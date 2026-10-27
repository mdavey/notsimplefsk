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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* merged view of CRC-good frames collected across all ppm trials */
typedef struct { size_t s0; unsigned len; const unsigned char *p; } gfrm;
static int cmp_gfrm(const void *a, const void *b)
{
    const gfrm *x = (const gfrm *)a, *y = (const gfrm *)b;
    return (x->s0 > y->s0) - (x->s0 < y->s0);
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
        "  -P, --preamble N   preamble symbols             (72)\n"
        "      --pulse P      gmsk | rcos | square        (gmsk)\n"
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
        case 'v': verbose++; break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    const char *err = nsfsk_params_validate(&p);
    if (err) { fprintf(stderr, "error: %s\n", err); return 1; }

    /* ---- read all of stdin ---- */
    size_t cap = 1u << 16, n16 = 0;
    int16_t *in = malloc(cap * sizeof(int16_t));
    for (;;) {
        if (n16 == cap) {
            if (cap > 1u << 26) { fprintf(stderr, "error: input too large\n"); return 1; }
            cap *= 2;
            in = realloc(in, cap * sizeof(int16_t));
        }
        size_t got = fread(in + n16, sizeof(int16_t), cap - n16, stdin);
        n16 += got;
        if (got == 0) break;
    }
    size_t N = n16 - (n16 & 1u);            /* whole samples only */
    if (N < (size_t)p.k * (p.preamble + 8u)) {
        fprintf(stderr, "rx: input too short for one frame\n");
        return 1;
    }

    unsigned nf = nsfsk_lpf_taps(&p);
    float complex *y = malloc(N * sizeof(float complex));
    if (nsfsk_rx_filter_slice(&p, in, N, y)) return 1;

    /* ---- matched-filter template from steady preamble ---- */
    unsigned tn = NSFSK_TMPL_SYMS * p.k;
    float complex *tmpl = malloc(tn * sizeof(float complex));
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
     * frame from different trials land within pos*2400e-6 samples of each
     * other and are deduped. */
    static rx_result tr[23];
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
    if (pfix) ppm_best = atof(pfix);
    for (unsigned t = 0; t < ntrial; t++) {
        double ppm;
        if (pfix)            ppm = atof(pfix);
        else                 ppm = (double)grid[t];
        rx_result *cur = &tr[t];
        if (ppm == 0.0) {
            rx_run(y, N, nf, tmpl, et, p, verbose > 1 ? verbose - 1 : verbose,
                   0, cur);
        } else {
            size_t N2;
            float complex *y2 = resample_lin(y, N, 1.0 + ppm * 1e-6, &N2);
            rx_run(y2, N2, nf, tmpl, et, p,
                   verbose > 1 ? verbose - 1 : verbose, 0, cur);
            free(y2);
        }
        double qn = cur->qual / (double)(cur->ndec + 1u);
        if (verbose)
            fprintf(stderr, "rx: ppm trial %+6.0f -> ok=%u bad=%u qn=%.4f\n",
                    ppm, cur->ok, cur->bad, qn);
        if (cur->ok > ok_best || (cur->ok == ok_best && qn > best_qn)) {
            ok_best = cur->ok; best_qn = qn; ppm_best = ppm;
        }
    }

    static gfrm all[23 * RX_MAXFRAMES];
    unsigned na = 0;
    for (unsigned t = 0; t < ntrial; t++)
        for (unsigned i = 0; i < tr[t].nfr && na < 23u * RX_MAXFRAMES; i++) {
            all[na].s0  = tr[t].fr[i].s0;
            all[na].len = tr[t].fr[i].len;
            all[na].p   = tr[t].buf + tr[t].fr[i].off;
            na++;
        }
    qsort(all, na, sizeof all[0], cmp_gfrm);
    unsigned nout = 0; size_t nbytes = 0, acc_s0 = 0;
    for (unsigned i = 0; i < na; i++) {
        /* same frame seen by several trials: positions agree to within
         * the trial rescale (0.24%) plus refit jitter; a different real
         * frame never starts closer than one short burst (~2.5k syms) */
        if (nout && all[i].s0 < acc_s0 + 1200u + (size_t)(acc_s0 * 0.0025))
            continue;
        fwrite(all[i].p, 1, all[i].len, stdout);
        nbytes += all[i].len; acc_s0 = all[i].s0; nout++;
    }
    if (nout) fflush(stdout);
    fprintf(stderr, "rx: %u frame(s) merged from %u ok / %zu byte(s) decoded "
                    "(ppm=%+.0f); rate=%u fc=%.0f k=%u h=%.3f bt=%.2f "
                    "delay=%u pulse=%s\n",
            nout, ok_best, nbytes, ppm_best,
            p.rate, p.fc, p.k, p.h, p.bt, p.delay,
            p.pulse == LIQUID_CPFSK_GMSK ? "gmsk" :
            p.pulse == LIQUID_CPFSK_RCOS_FULL ? "rcos" : "square");

    for (unsigned t = 0; t < ntrial; t++) free(tr[t].buf);
    free(in); free(y); free(tmpl);
    return nout > 0 ? 0 : 1;
}
