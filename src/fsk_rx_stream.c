/*
 * fsk_rx_stream — live G4FSK receiver: continuous s16le mono audio on
 * stdin (e.g. `arecord -D default -f S16_LE -r 48000 -c 1`) -> payload
 * bytes on stdout, released per frame once its CRC passes.
 *
 * Burst framing: a CAUSAL DDC + single-pass LPF reference filter (same
 * taps as the offline chain) feeds an energy gate.  The gate closes when
 * the profile stays below 15% of the in-burst peak for 20 ms — the TX
 * inter-frame gap is 40 ms, so bursts never run together.  Each closed
 * burst's RAW samples (plus look-back and tail margins) are then re-run
 * through the exact offline chain (zero-phase double pass) and the shared
 * detector rx_run(), on a self-contained slice.
 *
 * Sample-rate ppm is handled per burst instead of globally: the same
 * coarse grid as offline ({0, +-300 .. +-1200}) plus a fine +-100/+-50
 * cascade refinement (escalating to the remaining offline-grid points
 * whenever a burst decoded nothing: garbage qn ranking, see
 * decode_range), tried newest-winner-first with EARLY EXIT at the first
 * trial that decodes cleanly (ok>0, bad==0).  Clock drift between bursts
 * is absorbed by the per-burst MF/prefix re-alignment, so no whole-file
 * rescan (and no EOF) is ever needed.
 *
 * Buffering: a flat sliding window (memmove-trimmed), never larger than
 * about two max-size frames; everything older than the decoded burst is
 * discarded after each decode.  stdout carries payload bytes only; all
 * diagnostics go to stderr.
 *
 * Exit code: 0 when at least one frame was decoded, 1 otherwise — the
 * same rule for a clean EOF and for Ctrl-C (SIGINT drains, summarizes).
 */
#include "common.h"
#include "rx_core.h"

#include <complex.h>
#include <errno.h>
#include <getopt.h>
#include <liquid/liquid.h>
#include <math.h>
#include <poll.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define GATE_Q       32u     /* pending-burst queue; >= 4x max closes per read chunk */
#define STREAM_REJECT_LIMIT 8u  /* see rx_run(): bound garbage sweeps per slice */
#define GATE_FLOOR   1e-4    /* absolute gate floor on mean |y|^2 (~ -40 dBFS)      */
#define READ_CHUNK   8192u   /* samples per read()                                   */

static volatile sig_atomic_t g_stop = 0;
static int g_gate_dbg = 0;   /* -vv gate open/close trace */
static void on_sigint(int sig) { (void)sig; g_stop = 1; }

static double nowd(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [options] < in.s16 > out.bin\n"
        "  e.g.  arecord -D default -f S16_LE -r 48000 -c 1 | %s\n"
        "\n"
        "Reads a continuous raw signed-16-bit little-endian mono stream on\n"
        "stdin, decodes G4FSK frames as their bursts complete and writes\n"
        "the payload bytes to stdout immediately (stderr for diagnostics).\n"
        "Bytes are released per frame, whole, only after the CRC passes;\n"
        "latency per frame is payload time + 40 ms gap + decode.\n"
        "Exit code: 0 if at least one frame was decoded, else 1 (also on\n"
        "EOF/Ctrl-C).\n"
        "\n"
        "PHY options MUST match the transmitter (defaults shown):\n"
        "  -r, --rate N       sample rate Hz              (48000)\n"
        "  -c, --fc F         carrier frequency Hz        (4000)\n"
        "  -k, --k N          samples per symbol, even    (20)\n"
        "  -H, --h-index F    modulation index h          (0.25)\n"
        "  -b, --bt F         Gaussian BT product          (0.8)\n"
        "  -d, --delay N      cpfsk filter delay symbols  (1)\n"
        "  -P, --preamble N   preamble symbols             (72, max 65536)\n"
        "      --pulse P      gmsk | rcos | square        (gmsk)\n"
        "      --max-frames N exit after N good frames\n"
        "      --idle-exit S  exit after S seconds with no carrier\n"
        "      --trial-budget S  wall-clock ppm-sweep budget per burst\n"
        "                         (0=default: burst audio time, min 2 s;\n"
        "                         -1=unlimited)\n"
        "  -v, --verbose      per-frame detail on stderr (-vv tracing)\n"
        "  -h, --help         this text\n",
        prog, prog);
}

/* ---- energy gate on the causal reference filter --------------------- */
typedef struct {
    nco_crcf     nco;
    firfilt_crcf lpfc;
    float       *hist;                 /* |f|^2 ring of win samples */
    unsigned     win, hpos, arm;
    double       acc, pmax, burst_pk, decay;
    int          armed, activity;
    unsigned     below;
    size_t       open_abs;
    size_t       q_open[GATE_Q], q_close[GATE_Q];
    unsigned     qn;
} gate;

static void gate_advance(gate *g, const int16_t *s, size_t n, size_t abs_base)
{
    for (size_t j = 0; j < n; j++) {
        float complex z, f;
        nco_crcf_mix_down(g->nco, (float complex)((float)s[j] / 32768.0f), &z);
        nco_crcf_step(g->nco);
        firfilt_crcf_execute_one(g->lpfc, z, &f);
        double e = (double)(crealf(f) * crealf(f) + cimagf(f) * cimagf(f));
        g->acc += e - g->hist[g->hpos];
        g->hist[g->hpos] = (float)e;
        if (++g->hpos == g->win) g->hpos = 0;
        double prof = g->acc / (double)g->win;
        if (prof < 0.0) prof = 0.0;

        /* long-window tracker (tau ~ 0.5 s) for OPENING: adapts to changing
         * levels within a second; FLOOR keeps a silent noisy line closed */
        g->pmax = prof > g->pmax ? prof : g->pmax * g->decay;
        double eth = 0.15 * g->pmax;
        if (eth < GATE_FLOOR) eth = GATE_FLOOR;
        if (prof >= eth) g->activity = 1;

        if (g->armed) {
            if (prof >= eth) {
                g->armed = 0; g->below = 0;
                g->open_abs = abs_base + j;
                g->burst_pk = prof;
            }
        } else {
            if (prof > g->burst_pk) g->burst_pk = prof;
            double ethc = 0.15 * g->burst_pk;    /* CLOSING hysteresis: a
                                                  * weak burst following a
                                                  * strong one still closes,
                                                  * and idle noise after a
                                                  * burst cannot reopen */
            if (ethc < GATE_FLOOR) ethc = GATE_FLOOR;
            if (prof < ethc) {
                if (++g->below >= g->arm) {
                    if (g_gate_dbg)
                        fprintf(stderr, "rx: gate open=%zu close=%zu (a=%zu)\n",
                                g->open_abs, abs_base + j + 1u, abs_base);
                    if (g->qn < GATE_Q) {
                        g->q_open[g->qn]  = g->open_abs;
                        g->q_close[g->qn] = abs_base + j + 1u;
                        g->qn++;
                    } else {
                        fprintf(stderr, "rx: burst queue full, gap ignored\n");
                    }
                    g->armed = 1; g->below = 0;
                }
            } else g->below = 0;
        }
    }
}

/* ---- decode context -------------------------------------------------- */
typedef struct {
    nsfsk_params p;
    int verbose;
    unsigned nf, back, tail, retain;
    size_t cap;
    float complex *tmpl; float et;
    float complex *y;                    /* scratch, cap entries       */
    int16_t *ring;
    size_t abs_lo;                       /* absolute index of ring[0]  */
    double last_ppm;
    double trial_budget;           /* s; 0 = auto (burst audio time, min 2),
                                    * <0 = unlimited */
    unsigned max_frames;           /* 0 = unlimited; enforced in deliver() */
    int pending_byte;              /* odd trailing byte of a split s16 sample */
    unsigned ok_tot, bad_tot, rej_tot;
    size_t bytes_tot;
} rx_ctx;

/* Trials on one burst slice; returns the winning result (caller frees). */
static rx_result decode_range(rx_ctx *X, const int16_t *raw, size_t n)
{
    rx_result best = {0};
    if (n < (size_t)X->p.k * (X->p.preamble + 8u)) return best;
    if (nsfsk_rx_filter_slice(&X->p, raw, n, X->y)) return best;

    int vb = X->verbose > 1 ? X->verbose - 1 : X->verbose;
    double cand[14]; unsigned nc = 0;
    cand[nc++] = X->last_ppm;                       /* newest winner first */
    const int coarse[8] = {0, 300, -300, 600, -600, 900, -900, 1200};
    for (unsigned i = 0; i < 8u; i++)
        if ((double)coarse[i] != X->last_ppm) cand[nc++] = (double)coarse[i];
    /* -1200 only if last_ppm wasn't it either */
    if (-1200.0 != X->last_ppm) cand[nc++] = -1200.0;

    double best_qn = -1.0, ppm_best = X->last_ppm;
    int accepted = 0;
    /* wall-clock budget for the trial sweep: a garbage burst must never
     * stall the reader far enough behind a live capture that audio is
     * overrun; trials already spent are kept (best-so-far still wins) */
    double budget = X->trial_budget > 0.0 ? X->trial_budget
                  : (X->trial_budget < 0.0 ? 0.0
                     : (double)n / (double)X->p.rate);
    if (budget > 0.0 && budget < 2.0) budget = 2.0;
    int    over_budget = 0;
    unsigned ntrial = 0;

    /* one ppm trial + bookkeeping: accept on first fully clean pass */
#define TRIAL(PPM)                                                        \
    do {                                                                  \
        double ppm_ = (PPM);                                              \
        double t_trial = nowd();                                          \
        if (g_stop) { rx_result_free(&best);                        \
                      goto trials_done; }                            \
        if (budget > 0.0 && nowd() - t_burst > budget) {                  \
            over_budget = 1; goto trials_done; }                          \
        ntrial++;                                                         \
        rx_result cur = {0};                                              \
        if (ppm_ == 0.0) {                                                \
            rx_run(X->y, n, X->nf, X->tmpl, X->et, X->p, vb,              \
                   STREAM_REJECT_LIMIT, 0u, 0u, &cur);                    \
        } else {                                                          \
            size_t N2;                                                    \
            float complex *y2 =                                           \
                resample_lin(X->y, n, 1.0 + ppm_ * 1e-6, &N2);            \
            if (y2 != NULL) {                                             \
                rx_run(y2, N2, X->nf, X->tmpl, X->et, X->p, vb,           \
                       STREAM_REJECT_LIMIT, 0u, 0u, &cur);                \
                free(y2);                                                 \
            }                                                             \
        }                                                                 \
        if (X->verbose > 1)                                               \
            fprintf(stderr, "rx: trial ppm=%+.0f -> ok=%u bad=%u (%.0f ms)\n", \
                    ppm_, cur.ok, cur.bad, 1000.0 * (nowd() - t_trial));  \
        if (cur.ok > 0 && cur.bad == 0) {          /* early exit */       \
            rx_result_free(&best);                                        \
            best = cur; ppm_best = ppm_; accepted = 1;                    \
            goto trials_done;                                             \
        }                                                                 \
        double qn = cur.qual / (double)(cur.ndec + 1u);                   \
        if (cur.ok > best.ok || (cur.ok == best.ok && qn > best_qn)) {    \
            rx_result_free(&best);                                        \
            best = cur; best_qn = qn; ppm_best = ppm_;                    \
        } else {                                                          \
            rx_result_free(&cur);                                         \
        }                                                                 \
    } while (0)                                                           

    double t_burst = nowd();
    for (unsigned t = 0; t < nc; t++) TRIAL(cand[t]);
    /* the offline RX cascades its fine grid and bootstraps it with frames
     * decoded elsewhere in the file; a lone long burst needs guaranteed
     * coverage: the true correction sits <= 150 ppm from one of the two
     * coarse neighbours, so probe both sides from the coarse winner */
    if (!accepted && !over_budget) {
        for (unsigned f = 0; f < 6u && !accepted; f++) {
            static const double fine[6] = {-150.0, -100.0, -50.0,
                                           50.0, 100.0, 150.0};
            TRIAL(ppm_best + fine[f]);
        }
        /* Nothing decoded anywhere: the qn ranking that recentred that
         * fine sweep was garbage, so it may have probed around the wrong
         * coarse winner (observed: a -100 ppm channel lost between the
         * 0 and -300 coarse points while a garbage +800 winner ranked
         * first).  The offline RX dropped adaptivity for exactly this
         * reason (fsk_rx.c dense-grid comment); do the same here as an
         * escalation -- the remaining offline-grid points give fixed
         * coverage of the +-300 coarse gaps.  Cost is bounded and only
         * paid by bursts that decoded nothing anyway. */
        static const double esc[14] = { -250, -200, -150, -100, -50,
                                        50, 100, 150, 200, 250,
                                        450, -450, 1000, -1000 };
        for (unsigned g = 0; g < 14u && !accepted; g++)
            TRIAL(esc[g]);
    }
#undef TRIAL
trials_done:
    if (over_budget && X->verbose)
        fprintf(stderr, "rx: trial budget %.1fs hit after %u trial(s)\n",
                budget, ntrial);
    X->last_ppm = ppm_best;
    return best;
}

static void deliver(rx_ctx *X, rx_result *R, size_t stamp_abs)
{
    /* --max-frames is exact: a burst (or the final drain) may contain more
     * frames than the remaining budget, so deliver only the first `want`
     * of its recorded frames (fr[] stores each payload's offset/length). */
    size_t nbytes = R->nbytes;
    unsigned ok = R->ok;
    if (X->max_frames && X->ok_tot + ok > X->max_frames) {
        unsigned want = X->max_frames - X->ok_tot;
        if (want > R->nfr) want = R->nfr;   /* OOM bookkeeping degradation */
        ok = want;
        if (want == 0)
            nbytes = 0;
        else if (want < R->nfr)
            nbytes = R->fr[want].off;       /* start of the next frame */
        else
            nbytes = R->fr[R->nfr - 1u].off + R->fr[R->nfr - 1u].len;
    }
    if (nbytes) { fwrite(R->buf, 1, nbytes, stdout); fflush(stdout); }
    X->ok_tot   += ok;
    X->bad_tot  += R->bad;
    X->bytes_tot += nbytes;
    if (R->ok == 0 && R->bad == 0) X->rej_tot++;
    if (X->verbose)
        fprintf(stderr, "rx: burst @%.1f ms ok=%u bad=%u byte(s)=%zu ppm=%+.0f\n",
                1000.0 * (double)stamp_abs / (double)X->p.rate,
                ok, R->bad, nbytes, X->last_ppm);
    rx_result_free(R);
}

int main(int argc, char **argv)
{
    rx_ctx X;
    memset(&X, 0, sizeof X);
    nsfsk_params_default(&X.p);

    static struct option longopts[] = {
        {"rate",       required_argument, 0, 'r'},
        {"fc",         required_argument, 0, 'c'},
        {"k",          required_argument, 0, 'k'},
        {"h-index",    required_argument, 0, 'H'},
        {"bt",         required_argument, 0, 'b'},
        {"delay",      required_argument, 0, 'd'},
        {"preamble",   required_argument, 0, 'P'},
        {"pulse",      required_argument, 0,  1 },
        {"max-frames", required_argument, 0,  2 },
        {"idle-exit",  required_argument, 0,  3 },
        {"trial-budget", required_argument, 0, 4 },
        {"verbose",    no_argument,       0, 'v'},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    unsigned max_frames = 0;
    double idle_exit = 0.0;
    int opt;
    while ((opt = getopt_long_only(argc, argv, "r:c:k:H:b:d:P:vh",
                                   longopts, NULL)) != -1) {
        switch (opt) {
        case 'r': X.p.rate = (unsigned)strtoul(optarg, NULL, 0); break;
        case 'c': X.p.fc = strtof(optarg, NULL); break;
        case 'k': X.p.k = (unsigned)strtoul(optarg, NULL, 0); break;
        case 'H': X.p.h = strtof(optarg, NULL); break;
        case 'b': X.p.bt = strtof(optarg, NULL); break;
        case 'd': X.p.delay = (unsigned)strtoul(optarg, NULL, 0); break;
        case 'P': X.p.preamble = (unsigned)strtoul(optarg, NULL, 0); break;
        case  1 : if (nsfsk_parse_pulse(optarg, &X.p.pulse)) {
                      fprintf(stderr, "error: bad --pulse '%s'\n", optarg);
                      return 1;
                  }
                  break;
        case  2 : max_frames = (unsigned)strtoul(optarg, NULL, 0); break;
        case  3 : idle_exit = strtod(optarg, NULL); break;
        case  4 : X.trial_budget = strtod(optarg, NULL); break;
        case 'v': X.verbose++; break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    const char *err = nsfsk_params_validate(&X.p);
    if (err) { fprintf(stderr, "error: %s\n", err); return 1; }

    g_gate_dbg = X.verbose > 1;
    X.max_frames = max_frames;
    X.pending_byte = -1;

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_sigint;                 /* no SA_RESTART: read()
                                                * returns EINTR on Ctrl-C */
    sigaction(SIGINT, &sa, NULL);

    /* ---- window geometry (see header comment) ---- */
    const nsfsk_params *p = &X.p;
    X.nf = nsfsk_lpf_taps(p);
    unsigned win = (unsigned)fminf((float)(p->rate / 50),
                                   0.5f * (float)(p->preamble * p->k));
    if (win < p->k) win = p->k;                /* must match rx_run() */
    unsigned arm  = p->rate / 50u;             /* 20 ms; TX gap is 40 ms */
    X.back   = win + 4u * p->k;                /* rx_run's MF look-back */
    X.tail   = (X.nf - 1u) / 2u + (2u * p->delay + 3u + 2u) * p->k;
    X.retain = X.back + arm + (p->preamble + 4u) * p->k + 2u * X.nf
             + (X.nf - 1u) / 2u;
    size_t maxframe = (size_t)(p->preamble + NSFSK_SYNC_LEN
                               + 4u * (NSFSK_HDR_LEN + NSFSK_MAX_PAYLOAD
                                       + NSFSK_CRC_LEN)) * p->k;
    /* worst-case occupancy: the ring holds up to retain + 2 read chunks
     * when a gate opens (the idle slide runs once per main-loop pass), and
     * open -> close+tail spans one max frame PLUS the TX trailing pad PLUS
     * the gate-close overshoot (causal LPF fill ~ win + 20 ms arm); then
     * one more chunk of headroom for the loop's read-space check.  A cap
     * that only covered maxframe overflowed mid-burst after any idle lead
     * (real audio path: silence before TX) and force-decoded frame 1 away;
     * pinned by stream.sh idle-lead*. */
    X.cap = X.retain + 2u * READ_CHUNK
          + maxframe + (size_t)nsfsk_tx_pad(p) * p->k + 2u * win + 2u * arm
          + X.tail + READ_CHUNK;

    X.ring = malloc(X.cap * sizeof(int16_t));
    X.y    = malloc(X.cap * sizeof(float complex));
    if (X.ring == NULL || X.y == NULL) {
        fprintf(stderr, "rx: out of memory (window %zu samples)\n", X.cap);
        return 1;
    }

    /* matched-filter template (identical to offline rx) */
    unsigned tn = NSFSK_TMPL_SYMS * p->k;
    X.tmpl = malloc(tn * sizeof(float complex));
    if (X.tmpl == NULL) { fprintf(stderr, "rx: out of memory\n"); return 1; }
    if (nsfsk_build_template(p, NSFSK_TMPL_SYM0(X.p), NSFSK_TMPL_SYMS,
                             X.tmpl) != tn) {
        fprintf(stderr, "error: template build failed\n");
        return 1;
    }
    X.et = 0.0f;
    for (unsigned i = 0; i < tn; i++) X.et += crealf(X.tmpl[i] * conjf(X.tmpl[i]));
    if (nsfsk_rx_cutoff(p) >= 2.0f * p->fc - nsfsk_rx_cutoff(p)) {
        fprintf(stderr, "error: fc=%.0f too low to reject the mixer image "
                        "for these PHY params\n", p->fc);
        return 1;
    }

    gate g;
    memset(&g, 0, sizeof g);
    g.win = win; g.arm = arm;
    g.decay = exp(-2.0 / (double)p->rate);     /* pmax tau ~ 0.5 s */
    g.hist = calloc(win, sizeof(float));
    if (g.hist == NULL) { fprintf(stderr, "rx: out of memory\n"); return 1; }
    g.armed = 1;
    g.nco = nco_crcf_create(LIQUID_NCO);
    g.lpfc = firfilt_crcf_create_kaiser(X.nf, nsfsk_rx_cutoff(p)
                                        / (float)p->rate, 80.0f, 0.0f);
    if (g.nco == NULL || g.lpfc == NULL) {
        fprintf(stderr, "rx: gate filter allocation failed\n");
        return 1;
    }
    nco_crcf_set_frequency(g.nco, 2.0f * (float)M_PI * p->fc / (float)p->rate);

    size_t consumed = 0;                       /* samples read since start */
    size_t safe_after = 0;                     /* absolute: every frame that
                                                * STARTS before this has been
                                                * delivered or discarded */
    int eof = 0;
    double t_act = nowd();

    while (!g_stop) {
        size_t len = consumed - X.abs_lo;
        if (X.cap - len < READ_CHUNK) {
            /* pathologic: carrier without a 20 ms gap (or decoder wedged on
             * noise) filled the window; decode everything held, reset */
            if (X.verbose)
                fprintf(stderr, "rx: window full with no gap, forcing decode\n");
            rx_result R = decode_range(&X, X.ring, len);
            deliver(&X, &R, X.abs_lo);
            safe_after = consumed;
            size_t nl = consumed > X.retain ? consumed - X.retain : 0u;
            /* an open burst is the one thing we must not trim away: keep
             * from its head so it still decodes when it closes, provided
             * keeping it leaves real room to grow */
            if (!g.armed) {
                size_t keep = g.open_abs > X.back ? g.open_abs - X.back : 0u;
                if (keep >= X.abs_lo && keep < nl &&
                    consumed - keep + 2u * READ_CHUNK <= X.cap)
                    nl = keep;
            }
            /* drop only the queued bursts whose slice data got trimmed */
            unsigned q0 = 0;
            while (q0 < g.qn &&
                   (g.q_open[q0] > X.back ? g.q_open[q0] - X.back : 0u) < nl)
                q0++;
            for (unsigned i = q0; i < g.qn; i++) {
                g.q_open[i - q0]  = g.q_open[i];
                g.q_close[i - q0] = g.q_close[i];
            }
            g.qn -= q0;
            memmove(X.ring, X.ring + (nl - X.abs_lo),
                    (consumed - nl) * sizeof(int16_t));
            X.abs_lo = nl;
            len = consumed - X.abs_lo;
        }
        /* If the producer stalls, a plain blocking read() would never let the
         * idle timer fire; poll for the remaining idle budget first.  When
         * idle_exit is 0 this is skipped and the read stays blocking. */
        if (idle_exit > 0.0) {
            double rem = idle_exit - (nowd() - t_act);
            if (rem <= 0.0) {
                if (X.verbose) fprintf(stderr, "rx: idle, exiting\n");
                break;
            }
            if (rem > 3600.0) rem = 3600.0;
            struct pollfd pfd;
            pfd.fd = 0; pfd.events = POLLIN; pfd.revents = 0;
            int pr = poll(&pfd, 1, (int)(rem * 1000.0) + 1);
            if (pr == 0) {
                if (X.verbose) fprintf(stderr, "rx: idle, exiting\n");
                break;
            }
            if (pr < 0) {
                if (errno == EINTR) { if (g_stop) break; }
                else { fprintf(stderr, "rx: poll error: %s\n", strerror(errno)); break; }
            }
        }
        size_t want = X.cap - len < READ_CHUNK ? X.cap - len : READ_CHUNK;
        if (want == 0) {
            fprintf(stderr, "rx: internal window full\n");
            break;
        }
        /* read bytes, not samples: a pipe may split an s16 sample across
         * reads, and dropping the odd byte would desynchronize the whole
         * stream.  Carry it into the next read (pending_byte). */
        unsigned char bytes[READ_CHUNK * 2];
        size_t boff = 0;
        if (X.pending_byte >= 0) { bytes[0] = (unsigned char)X.pending_byte; boff = 1; }
        ssize_t r = read(0, bytes + boff, want * sizeof(int16_t) - boff);
        size_t got;
        if (r > 0) {
            X.pending_byte = -1;
            size_t total = (size_t)r + boff;
            if (total & 1u) { X.pending_byte = bytes[total - 1u]; total--; }
            memcpy(X.ring + len, bytes, total);
            got = total / 2u;
        } else if (r == 0) {
            eof = 1; got = 0;
            X.pending_byte = -1;          /* trailing half sample unusable */
        } else {
            if (errno == EINTR) { if (g_stop) break; continue; }
            fprintf(stderr, "rx: read error: %s\n", strerror(errno));
            eof = 1; got = 0;
            X.pending_byte = -1;
        }
        if (got) {
            gate_advance(&g, X.ring + len, got, consumed);
            consumed += got;
            if (g.activity) { t_act = nowd(); g.activity = 0; }
        }

        /* decode every burst whose close+tail samples have arrived; the
         * ring is trimmed to `retain` samples after each decode */
        while (g.qn && consumed >= g.q_close[0] + X.tail) {
            size_t lo = g.q_open[0] > X.back ? g.q_open[0] - X.back : 0u;
            if (lo < X.abs_lo) lo = X.abs_lo;
            size_t hi = consumed < g.q_close[0] + X.tail
                      ? consumed : g.q_close[0] + X.tail;
            rx_result R = decode_range(&X, X.ring + (lo - X.abs_lo), hi - lo);
            deliver(&X, &R, g.q_open[0]);
            /* the last frame symbol is at least arm samples before the
             * gate closed, so nothing decoded can reappear ahead of this */
            if (g.q_close[0] > arm + 1u && g.q_close[0] - arm + 1u > safe_after)
                safe_after = g.q_close[0] - arm + 1u;
            /* pop */
            for (unsigned i = 1; i < g.qn; i++) {
                g.q_open[i - 1u] = g.q_open[i];
                g.q_close[i - 1u] = g.q_close[i];
            }
            g.qn--;
            /* trim: keep the newest `retain` samples (covers the next
             * burst's look-back, see header math) */
            size_t nl = hi > X.retain ? hi - X.retain : 0u;
            if (safe_after < nl) safe_after = nl;      /* never un-deliver */
            memmove(X.ring, X.ring + (nl - X.abs_lo),
                    (consumed - nl) * sizeof(int16_t));
            X.abs_lo = nl;
        }

        /* idle between bursts (armed, queue drained): nothing can start
         * before `consumed`, so slide the window down to `retain` samples;
         * without this, a long silence fills the ring and the force-decode
         * path would fire mid-burst and destroy the next frame's head */
        if (g.armed && g.qn == 0 && consumed - X.abs_lo > X.retain) {
            size_t nl = consumed - X.retain;
            if (nl > X.abs_lo) {
                memmove(X.ring, X.ring + (nl - X.abs_lo),
                        (consumed - nl) * sizeof(int16_t));
                X.abs_lo = nl;
                if (safe_after < nl) safe_after = nl;
            }
        }

        if (eof || g_stop) break;
        if (X.max_frames && X.ok_tot >= X.max_frames) break;
        if (idle_exit > 0.0 && nowd() - t_act > idle_exit) {
            if (X.verbose) fprintf(stderr, "rx: idle, exiting\n");
            break;
        }
    }

    /* drain: a complete final burst may sit past the last gate close */
    if (consumed > X.abs_lo) {
        size_t lo = safe_after > X.abs_lo ? safe_after : X.abs_lo;
        if (consumed - lo >= (size_t)p->k * (p->preamble + 8u)) {
            rx_result R = decode_range(&X, X.ring + (lo - X.abs_lo),
                                       consumed - lo);
            deliver(&X, &R, lo);
        }
    }

    fprintf(stderr, "rx: %u frame(s) ok, %u bad, %zu byte(s) decoded, "
                    "%u burst(s) rejected; last ppm=%+.0f; rate=%u fc=%.0f "
                    "k=%u h=%.3f bt=%.2f delay=%u pulse=%s\n",
            X.ok_tot, X.bad_tot, X.bytes_tot, X.rej_tot, X.last_ppm,
            p->rate, p->fc, p->k, p->h, p->bt, p->delay,
            p->pulse == LIQUID_CPFSK_GMSK ? "gmsk" :
            p->pulse == LIQUID_CPFSK_RCOS_FULL ? "rcos" : "square");

    free(X.ring); free(X.y); free(X.tmpl); free(g.hist);
    nco_crcf_destroy(g.nco);
    firfilt_crcf_destroy(g.lpfc);
    return X.ok_tot > 0 ? 0 : 1;
}
