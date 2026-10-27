/*
 * channel — test impairment channel: s16le mono in -> impaired s16le out.
 *
 * Impairments (applied in this order):
 *   -d  integer sample delay (zeros prepended)
 *   -D  fractional sample delay (linear interpolation, resampling walk)
 *   -r  sample-rate offset in ppm (output stretched by 1 + r*1e-6)
 *   -g  gain in dB
 *   -f  carrier frequency offset in Hz (analytic-shift via Hilbert FIR)
 *   -n  additive white Gaussian noise sigma, full-scale units (e.g. 0.01)
 *   -s  RNG seed (default 1)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double gauss_pair_spare;
static int    gauss_pair_have;
static double gauss_next(unsigned long *st)
{
    /* xorshift64* + Box-Muller */
    double u1, u2;
    if (gauss_pair_have) { gauss_pair_have = 0; return gauss_pair_spare; }
    do {
        *st ^= *st << 13; *st ^= *st >> 7; *st ^= *st << 17;
        u1 = (double)(*st % 1000000000u) / 1000000000.0 + 1e-9;
        *st ^= *st << 13; *st ^= *st >> 7; *st ^= *st << 17;
        u2 = (double)(*st % 1000000000u) / 1000000000.0;
    } while (u1 <= 1e-9);
    gauss_pair_spare = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
    gauss_pair_have  = 1;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [-d delay] [-D fdelay] [-r ppm] [-g dB] [-f Hz] [-n sigma] "
        "[-s seed] [-R rate] < in.s16 > out.s16\n", prog);
}

int main(int argc, char **argv)
{
    long   delay = 0;
    long   rate = 48000;
    double fdelay = 0.0, ppm = 0.0, gain_db = 0.0, cfo = 0.0, sigma = 0.0;
    unsigned long seed = 1;
    int opt;
    while ((opt = getopt(argc, argv, "d:D:r:g:f:n:s:R:h")) != -1) {
        switch (opt) {
        case 'R': rate = strtol(optarg, NULL, 0); break;
        case 'd': delay = strtol(optarg, NULL, 0); break;
        case 'D': fdelay = strtod(optarg, NULL); break;
        case 'r': ppm = strtod(optarg, NULL); break;
        case 'g': gain_db = strtod(optarg, NULL); break;
        case 'f': cfo = strtod(optarg, NULL); break;
        case 'n': sigma = strtod(optarg, NULL); break;
        case 's': seed = strtoul(optarg, NULL, 0); break;
        default: usage(argv[0]); return 1;
        }
    }
    if (seed == 0) seed = 1;

    size_t cap = 1u << 16, n = 0;
    int16_t *in = malloc(cap * sizeof(int16_t));
    for (;;) {
        if (n == cap) {
            if (cap > 1u << 26) { fprintf(stderr, "channel: input too large\n"); return 1; }
            cap *= 2; in = realloc(in, cap * sizeof(int16_t));
        }
        size_t got = fread(in + n, sizeof(int16_t), cap - n, stdin);
        n += got;
        if (!got) break;
    }
    n -= n & 1u;

    /* ---- resampling walk: integer+fractional delay and ppm rate error ---- */
    double ratio = 1.0 + ppm * 1e-6;               /* out n -> in pos n*ratio - delay */
    size_t m = (size_t)(((double)n + (double)delay) / ratio);
    float *x = malloc((m + 4) * sizeof(float));
    for (size_t i = 0; i < m; i++) {
        double pos = (double)i * ratio - (double)delay - fdelay;
        if (pos < 0.0 || pos > (double)(n - 1)) { x[i] = 0.0f; continue; }
        size_t i0 = (size_t)pos;
        double f = pos - (double)i0;
        x[i] = (float)((1.0 - f) * (double)in[i0] + f * (double)in[i0 + 1]);
    }

    /* ---- gain ---- */
    double amp = pow(10.0, gain_db / 20.0);
    for (size_t i = 0; i < m; i++) x[i] *= (float)amp;

    /* ---- CFO via windowed-FIR Hilbert transform ---- */
    if (cfo != 0.0) {
        const int HL = 127, c = (HL - 1) / 2;      /* odd length, odd-sym h */
        double h[127];
        for (int j = 0; j < HL; j++) {
            int k = j - c;
            double v = 0.0;
            if (k != 0 && (k & 1)) v = 2.0 / (M_PI * (double)k);   /* ideal hilbert */
            /* Hann window */
            v *= 0.5 - 0.5 * cos(2.0 * M_PI * (double)j / (double)(HL - 1));
            h[j] = v;
        }
        float *xi = malloc((m + HL) * sizeof(float));
        memcpy(xi + c, x, m * sizeof(float));      /* zero-padded */
        double th = 0.0, dtheta = 2.0 * M_PI * cfo / (double)rate;
        for (size_t i = 0; i < m; i++) {
            double hb = 0.0;
            for (int j = 0; j < HL; j++) hb += h[j] * (double)xi[i + j];
            /* real{ (x + j*hb) * e^{j th} } */
            x[i] = (float)((double)x[i] * cos(th) - hb * sin(th));
            th += dtheta;
            if (th > 2.0 * M_PI) th -= 2.0 * M_PI;
        }
        free(xi);
    }

    /* ---- AWGN + clip ---- */
    if (sigma > 0.0)
        for (size_t i = 0; i < m; i++)
            x[i] += (float)(sigma * 32768.0 * gauss_next(&seed));

    int16_t *out = malloc((m + 8) * sizeof(int16_t));
    for (size_t i = 0; i < m; i++) {
        double v = x[i];
        if (v > 32767.0) v = 32767.0;
        if (v < -32768.0) v = -32768.0;
        out[i] = (int16_t)lround(v);
    }
    fwrite(out, sizeof(int16_t), m, stdout);
    fprintf(stderr, "channel: %zu -> %zu samples (delay=%ld fdelay=%+.2f "
                    "ppm=%+.1f gain=%+.1fdB cfo=%+.1fHz sigma=%g)\n",
            n, m, delay, fdelay, ppm, gain_db, cfo, sigma);
    return 0;
}
