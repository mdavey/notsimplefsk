/*
 * fsk_rx — receive G4FSK frames: s16le mono audio on stdin -> payload bytes
 * on stdout.  PHY options must match the transmitter (see --help).
 *
 * Pipeline (offline, whole input buffered): see src/rx_core.h for the
 * DDC + LPF, energy gate, matched filter, prefix/CFO alignment sweep and
 * differential CRC-gated decode — all shared with fsk_rx_stream.
 *
 * Sample-rate ppm error is handled offline by resample trials over the
 * whole capture (pick the run with most frames, then best per-decision
 * quality); total detection lag is absorbed by the lag search.  No AGC
 * (scale tracking covers gain error for detection purposes).
 */
#include "common.h"
#include "rx_core.h"

#include <complex.h>
#include <getopt.h>
#include <liquid/liquid.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Offline acquisition budget: a continuous undecodable burst (tone, fading,
 * reverb) otherwise triggers a full MF/prefix sweep after every energy dip.
 * Calibrated against roundtrip/stream and the p2 pathological cases:
 *  - valid single frames need up to 3 acquisitions (periodic preamble lobes);
 *    some noisy/multipath valid frames need 4 within a burst,
 *  - a garbage burst can use 4 before the 40 ms gap, then the real frame 1,
 *  - no-decode p2 tone/fading trials used 38-64 acquisitions without a cap.
 * Keep the old reject_limit (64): lowering it breaks valid alignments before
 * the acquisition cap can do its job. */
#define OFFLINE_REJECT_LIMIT   64u
#define OFFLINE_MAX_ACQ_TOTAL  5u
#define OFFLINE_MAX_ACQ_BURST  4u
/* The budget is only needed for expensive (large-k) sweeps.  At k<=24 a
 * full no-decode trial costs a fraction of a second even with many
 * acquisitions, and capping it regresses stochastic/ppm-ramp baselines
 * (e.g. p0 rayleigh-fd1 needs >8 acquisitions before its first good frame). */
#define OFFLINE_CAP_MIN_K      32u

static unsigned env_cap(const char *name, unsigned def)
{
    const char *s = getenv(name);
    char *end = NULL;
    unsigned long v;
    if (s == NULL || *s == '\0') return def;
    v = strtoul(s, &end, 0);
    if (end == s || *end != '\0' || v > 1000000UL) {
        fprintf(stderr, "rx: bad %s '%s' (using %u)\n", name, s, def);
        return def;
    }
    return (unsigned)v;
}

/* Cross-trial merge window.  Copies of one frame can be recovered at
 * different lag offsets per ppm trial; the observed starts differ by ~8k
 * (p2 echo) and ~5.8k (p4 noise).  Cover the lag sweep with a 2x margin,
 * but never exceed 25% of the minimum legal frame spacing: two genuinely
 * distinct frames (even with identical payloads) are at least that far
 * apart, so adjacent frames cannot be merged. */
static double merge_tolerance(const nsfsk_params *p)
{
    double min_frame_syms = (double)p->preamble + NSFSK_SYNC_LEN
                          + (double)NSFSK_WIRE_SYMS;
    double min_spacing = min_frame_syms * (double)p->k;
    double lag_max = 2.0 * (double)p->delay + 3.0;
    double lag_bound = 2.0 * lag_max * (double)p->k
                     + 4.0 * (double)p->k + 128.0;
    double floor = 4.0 * (double)p->k + 128.0;
    double cap = 0.25 * min_spacing;
    return fmin(fmax(floor, lag_bound), cap);
}

/* merged view of CRC-good frames collected across all ppm trials.
 * t0 is the frame start mapped back to the un-resampled capture
 * (resample_lin places input position P at output P/ratio, so the original
 * position is s0*ratio), so copies of one frame from different trials
 * collapse to the same t0 while distinct frames stay one frame period
 * apart. */
typedef struct { double t0; unsigned len; const unsigned char *p; } gfrm;
static int cmp_gfrm(const void *a, const void *b)
{
    const gfrm *x = (const gfrm *)a, *y = (const gfrm *)b;
    return (x->t0 > y->t0) - (x->t0 < y->t0);
}

static void usage(const char *prog){
    fprintf(stderr,
        "usage: %s [options] < in.s16 > out.bin\n"
        "\n"
        "Reads raw signed-16-bit little-endian mono audio on stdin, finds\n"
        "G4FSK frames and writes the recovered payload bytes to stdout.\n"
        "Diagnostics go to stderr.\n"
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
        "      --fec N        RS parity bytes             (16|32|64|84, 32)\n"
        "  -v, --verbose      per-frame detail on stderr (-vv tracing)\n"
        "  -h, --help         this text\n",
        prog);
}

int main(int argc, char **argv)
{
    nsfsk_params p;
    nsfsk_params_default(&p);
    int verbose = 0;

    static struct option longopts[] = {
        {"rate",     required_argument, 0, 'r'},
        {"fc",       required_argument, 0, 'c'},
        {"k",        required_argument, 0, 'k'},
        {"h-index",  required_argument, 0, 'H'},
        {"bt",       required_argument, 0, 'b'},
        {"delay",    required_argument, 0, 'd'},
        {"preamble", required_argument, 0, 'P'},
        {"pulse",    required_argument, 0,  1 },
        {"fec",      required_argument, 0,  2 },
        {"verbose",  no_argument,       0, 'v'},
        {"help",     no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long_only(argc, argv, "r:c:k:H:b:d:P:vh",
                                   longopts, NULL)) != -1) {
        switch (opt) {
        case 'r': p.rate = (unsigned)strtoul(optarg, NULL, 0); break;
        case 'c': p.fc = strtof(optarg, NULL); break;
        case 'k': p.k = (unsigned)strtoul(optarg, NULL, 0); break;
        case 'H': p.h = strtof(optarg, NULL); break;
        case 'b': p.bt = strtof(optarg, NULL); break;
        case 'd': p.delay = (unsigned)strtoul(optarg, NULL, 0); break;
        case 'P': p.preamble = (unsigned)strtoul(optarg, NULL, 0); break;
        case 1:   if (nsfsk_parse_pulse(optarg, &p.pulse)) {
                      fprintf(stderr, "error: bad --pulse '%s'\n", optarg);
                      return 1;
                  }
                  break;
        case 2: {
            char *end = NULL;
            unsigned long v = strtoul(optarg, &end, 0);
            if (end == optarg || *end != '\0') {
                fprintf(stderr, "error: bad --fec '%s'\n", optarg);
                return 1;
            }
            p.fec = (unsigned)v;
            break;
        }
        case 'v': verbose++; break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    const char *err = nsfsk_params_validate(&p);
    if (err) { fprintf(stderr, "error: %s\n", err); return 1; }

    unsigned offline_reject = env_cap("NSFSK_RX_REJECT_LIMIT",
                                      OFFLINE_REJECT_LIMIT);
    unsigned acq_total_default = p.k >= OFFLINE_CAP_MIN_K
                               ? OFFLINE_MAX_ACQ_TOTAL : 0u;
    unsigned acq_burst_default = p.k >= OFFLINE_CAP_MIN_K
                               ? OFFLINE_MAX_ACQ_BURST : 0u;
    unsigned offline_acq_total = env_cap("NSFSK_RX_ACQ_TOTAL",
                                         acq_total_default);
    unsigned offline_acq_burst = env_cap("NSFSK_RX_ACQ_BURST",
                                         acq_burst_default);

    /* ---- read all of stdin ---- */
    size_t cap = 1u << 16, n16 = 0;
    int16_t *in = malloc(cap * sizeof(int16_t));
    if (in == NULL) { fprintf(stderr, "rx: out of memory\n"); return 1; }
    for (;;) {
        if (n16 == cap) {
            if (cap > 1u << 26) { fprintf(stderr, "error: input too large\n"); return 1; }
            int16_t *ni = realloc(in, cap * 2u * sizeof(int16_t));
            if (ni == NULL) {
                fprintf(stderr, "rx: out of memory\n");
                free(in);
                return 1;
            }
            in = ni;
            cap *= 2;
        }
        size_t got = fread(in + n16, sizeof(int16_t), cap - n16, stdin);
        n16 += got;
        if (got == 0) break;
    }
    size_t N = n16;                         /* fread reads whole samples */
    if (N < (size_t)p.k * (p.preamble + 8u)) {
        fprintf(stderr, "rx: input too short for one frame\n");
        return 1;
    }

    unsigned nf = nsfsk_lpf_taps(&p);
    float complex *y = malloc(N * sizeof(float complex));
    if (y == NULL) { fprintf(stderr, "rx: out of memory\n"); return 1; }
    if (nsfsk_rx_filter_slice(&p, in, N, y)) return 1;

    /* ---- matched-filter template from steady preamble ---- */
    unsigned tn = NSFSK_TMPL_SYMS * p.k;
    float complex *tmpl = malloc(tn * sizeof(float complex));
    if (tmpl == NULL) { fprintf(stderr, "rx: out of memory\n"); return 1; }
    if (nsfsk_build_template(&p, NSFSK_TMPL_SYM0(p), NSFSK_TMPL_SYMS, tmpl) != tn) {
        fprintf(stderr, "error: template build failed\n");
        return 1;
    }
    float et = 0.0f;
    for (unsigned i = 0; i < tn; i++) et += crealf(tmpl[i] * conjf(tmpl[i]));

    /* ---- sample-rate drift trials (offline timing recovery) ----
     * Constant ppm clock error walks the decision instants across the
     * symbol (~k/2 samples = death); resample the whole capture on a ppm
     * grid and keep the run with most frames, then the best mean
     * |phase|-per-decision quality.
     *
     * Every trial's result is KEPT and the CRC-good frames are MERGED
     * positionally at the end: per-frame sensitivity to residual rate
     * error differs by several samples between frames (within-frame slip
     * plus per-burst MF refit luck), so no single global ppm is
     * guaranteed to decode all frames of a multi-frame capture.  A
     * CRC-verified frame from any trial is safe to deliver; copies of one
     * frame from different trials are rescaled back to the original
     * timebase and deduped by position + payload (see the merge below). */
    static rx_result tr[23];
    double ppm_t[23];
    double best_qn = -1.0, ppm_best = 0.0;
    unsigned ok_best = 0;
    /* Fixed dense grid (no adaptive fine stage: with every trial at
     * ok=0 the qn ranking is unreliable and recenters the fine grid on
     * garbage, e.g. missing a +-100ppm channel entirely between the
     * 300-step coarse points).  Per-frame union across trials makes
     * coverage of the grid the thing that matters. */
    const int grid[23] = { 0, 50, -50, 100, -100, 150, -150, 200, -200,
                           250, -250, 300, -300, 450, -450, 600, -600,
                           900, -900, 1000, -1000, 1200, -1200 };
    const char *pfix = getenv("NSFSK_PPM");      /* dev knob: single trial */
    unsigned ntrial = pfix ? 1u : 23u;
    if (pfix) {
        char *end = NULL;
        double v = strtod(pfix, &end);
        /* keep 1 + ppm*1e-6 strictly positive/finite for resample_lin() */
        if (end == pfix || *end != '\0' || !isfinite(v) ||
            v <= -900000.0 || v >= 900000.0) {
            fprintf(stderr, "rx: bad NSFSK_PPM '%s'\n", pfix);
            return 1;
        }
        ppm_best = v;
    }
    for (unsigned t = 0; t < ntrial; t++) {
        double ppm;
        if (pfix)            ppm = atof(pfix);
        else                 ppm = (double)grid[t];
        ppm_t[t] = ppm;
        rx_result *cur = &tr[t];
        /* reject_limit: MF-strong-but-prefix-garbage candidates cost a full
         * acquisition sweep each; unbounded (0) lets a bt<0.5 capture with a
         * CRC-failing frame thrash the sweep every 4k samples for minutes.
         * After the cap, idle to the next silence and re-arm. */
        if (ppm == 0.0) {
            rx_run(y, N, nf, tmpl, et, p, verbose > 1 ? verbose - 1 : verbose,
                   offline_reject, offline_acq_total, offline_acq_burst, cur);
        } else {
            size_t N2;
            float complex *y2 = resample_lin(y, N, 1.0 + ppm * 1e-6, &N2);
            if (y2 != NULL) {
                rx_run(y2, N2, nf, tmpl, et, p,
                       verbose > 1 ? verbose - 1 : verbose,
                       offline_reject, offline_acq_total, offline_acq_burst,
                       cur);
                free(y2);
            }
        }
        double qn = cur->qual / (double)(cur->ndec + 1u);
        if (verbose)
            fprintf(stderr, "rx: ppm trial %+6.0f -> ok=%u bad=%u acq=%u qn=%.4f\n",
                    ppm, cur->ok, cur->bad, cur->acq, qn);
        if (cur->ok > ok_best || (cur->ok == ok_best && qn > best_qn)) {
            ok_best = cur->ok; best_qn = qn; ppm_best = ppm;
        }
    }

    /* ---- merge the trials ----
     * Every CRC-good frame from every trial is mapped back to the original
     * (un-resampled) timebase, then duplicates are identified by BOTH the
     * rescaled position and the payload.  Position alone is not enough:
     * the cross-trial spread grows with the capture length and eventually
     * exceeds the spacing of short frames, which silently discarded real
     * frames (2026-10-29: 400 one-byte frames -> 386 delivered, 520 -> 405). */
    size_t na = 0;
    for (unsigned t = 0; t < ntrial; t++) na += tr[t].nfr;
    gfrm *all = na ? malloc(na * sizeof *all) : NULL;
    if (na && all == NULL) {
        fprintf(stderr, "rx: out of memory merging trials\n");
        for (unsigned t = 0; t < ntrial; t++) rx_result_free(&tr[t]);
        free(in); free(y); free(tmpl);
        return 1;
    }
    size_t na_fill = 0;
    for (unsigned t = 0; t < ntrial; t++) {
        const double ratio = 1.0 + ppm_t[t] * 1e-6;
        for (unsigned i = 0; i < tr[t].nfr; i++) {
            const rx_frame *f = &tr[t].fr[i];
            /* resample_lin places input position P at output P/ratio, so the
             * un-resampled position is s0 * ratio */
            all[na_fill].t0  = (double)f->s0 * ratio;
            all[na_fill].len = f->len;
            all[na_fill].p   = tr[t].buf + f->off;
            na_fill++;
        }
    }
    if (na_fill) qsort(all, na_fill, sizeof all[0], cmp_gfrm);
    /* The lag sweep can recover the same frame at different offsets in
     * different ppm trials; see merge_tolerance().  The 25% spacing cap
     * keeps adjacent identical-payload frames distinct. */
    const double tol = merge_tolerance(&p);
    unsigned nout = 0, ndup = 0; size_t nbytes = 0;
    double acc_t0 = 0.0; unsigned acc_len = 0; const unsigned char *acc_p = NULL;
    for (size_t i = 0; i < na_fill; i++) {
        /* copies of one frame from different trials have the same payload
         * and land within refit jitter of each other after rescaling; two
         * genuinely different frames are at least a frame period apart and
         * (even for identical payloads) cannot satisfy both tests. */
        if (nout && fabs(all[i].t0 - acc_t0) < tol &&
            all[i].len == acc_len &&
            (acc_len == 0 || memcmp(all[i].p, acc_p, acc_len) == 0)) {
            ndup++;
            continue;
        }
        if (all[i].len) fwrite(all[i].p, 1, all[i].len, stdout);
        nbytes += all[i].len;
        acc_t0 = all[i].t0; acc_len = all[i].len; acc_p = all[i].p;
        nout++;
    }
    if (nout) fflush(stdout);
    if (verbose)
        fprintf(stderr, "rx: merge tol=%.1f sample(s), suppressed %u "
                        "duplicate copy/copies\n", tol, ndup);
    fprintf(stderr, "rx: %u frame(s) merged from %u ok / %zu byte(s) decoded "
                    "(ppm=%+.0f); rate=%u fc=%.0f k=%u h=%.3f bt=%.2f "
                    "delay=%u pulse=%s\n",
            nout, ok_best, nbytes, ppm_best,
            p.rate, p.fc, p.k, p.h, p.bt, p.delay,
            p.pulse == LIQUID_CPFSK_GMSK ? "gmsk" :
            p.pulse == LIQUID_CPFSK_RCOS_FULL ? "rcos" : "square");

    free(all);
    for (unsigned t = 0; t < ntrial; t++) rx_result_free(&tr[t]);
    free(in); free(y); free(tmpl);
    return nout > 0 ? 0 : 1;
}
