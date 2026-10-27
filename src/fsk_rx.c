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

static void usage(const char *prog)
{
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
     * |phase|-per-decision quality. */
    rx_result best = {0};
    double best_qn = -1.0, ppm_best = 0.0;
    const int coarse[9] = {0, 300, -300, 600, -600, 900, -900, 1200, -1200};
    const int fine[4]   = {-100, -50, 50, 100};
    for (unsigned t = 0; t < 13u; t++) {
        double ppm;
        if (t < 9u) ppm = (double)coarse[t];
        else        ppm = ppm_best + (double)fine[t - 9u];
        rx_result cur = {0};
        if (ppm == 0.0) {
            rx_run(y, N, nf, tmpl, et, p, verbose > 1 ? verbose - 1 : verbose,
                   0, &cur);
        } else {
            size_t N2;
            float complex *y2 = resample_lin(y, N, 1.0 + ppm * 1e-6, &N2);
            rx_run(y2, N2, nf, tmpl, et, p,
                   verbose > 1 ? verbose - 1 : verbose, 0, &cur);
            free(y2);
        }
        double qn = cur.qual / (double)(cur.ndec + 1u);
        if (cur.ok > best.ok || (cur.ok == best.ok && qn > best_qn)) {
            rx_result old = best;
            best = cur; best_qn = qn; ppm_best = ppm;
            free(old.buf);
        } else {
            free(cur.buf);
        }
    }

    if (best.nbytes) { fwrite(best.buf, 1, best.nbytes, stdout); fflush(stdout); }
    fprintf(stderr, "rx: %u frame(s) ok, %u bad, %zu byte(s) decoded "
                    "(ppm=%+.0f); rate=%u fc=%.0f k=%u h=%.3f bt=%.2f "
                    "delay=%u pulse=%s\n",
            best.ok, best.bad, best.nbytes, ppm_best,
            p.rate, p.fc, p.k, p.h, p.bt, p.delay,
            p.pulse == LIQUID_CPFSK_GMSK ? "gmsk" :
            p.pulse == LIQUID_CPFSK_RCOS_FULL ? "rcos" : "square");

    free(best.buf); free(in); free(y); free(tmpl);
    return best.ok > 0 ? 0 : 1;
}
