/* probe12: dump per-symbol differential phase of one frame in a capture,
 * for s0 offsets around the RX-selected alignment.  Usage:
 *   probe12 in.s16 rate fc k delay s0 lag n0 n1
 * Prints w  a_hat  phi/sc for w in [n0,n1).  DDC + zero-phase LPF copied
 * from fsk_rx. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <liquid/liquid.h>
#include "common.h"

int main(int argc, char **argv)
{
    if (argc < 9) { fprintf(stderr, "usage: probe12 in rate fc k delay s0 lag n0 n1\n"); return 1; }
    nsfsk_params p; nsfsk_params_default(&p);
    p.rate  = (unsigned)strtoul(argv[1], 0, 0);
    p.fc    = strtof(argv[2], 0);
    p.k     = (unsigned)strtoul(argv[3], 0, 0);
    p.delay = (unsigned)strtoul(argv[4], 0, 0);
    long s0   = atol(argv[5]);
    long lag  = atol(argv[6]);
    long n0   = atol(argv[7]), n1 = atol(argv[8]);

    size_t cap = 1u << 16, n16 = 0;
    int16_t *in = malloc(cap * sizeof(int16_t));
    for (;;) {
        if (n16 == cap) { cap *= 2; in = realloc(in, cap * sizeof(int16_t)); }
        size_t got = fread(in + n16, sizeof(int16_t), cap - n16, stdin);
        n16 += got; if (!got) break;
    }
    size_t N = n16 - (n16 & 1u);

    unsigned nf = nsfsk_lpf_taps(&p);
    float rs = (float)p.rate / (float)p.k;
    float dev = p.h * rs * 3.0f / 2.0f, cutoff = dev + 0.85f * rs;
    nco_crcf nco = nco_crcf_create(LIQUID_NCO);
    firfilt_crcf lpf = firfilt_crcf_create_kaiser(nf, cutoff / (float)p.rate, 80.0f, 0.0f);
    nco_crcf_set_frequency(nco, 2.0f * (float)M_PI * p.fc / (float)p.rate);
    float complex *y = malloc(N * sizeof(float complex));
    for (size_t i = 0; i < N; i++) {
        float complex z, f;
        nco_crcf_mix_down(nco, (float complex)((float)in[i] / 32768.0f), &z);
        nco_crcf_step(nco);
        firfilt_crcf_execute_one(lpf, z, &f);
        y[i] = f;
    }
    firfilt_crcf_reset(lpf);
    for (size_t i = N; i-- > 0; ) { float complex f; firfilt_crcf_execute_one(lpf, y[i], &f); y[i] = f; }

    /* calibrate sc over prefix symbols TMPL_SYM0..D0 like the RX does */
    unsigned D0 = p.preamble + NSFSK_SYNC_LEN;
    float num = 0, den = 0;
    for (unsigned w = NSFSK_TMPL_SYM0(p); w < D0; w++) {
        size_t a = (size_t)s0 + (size_t)(w + lag) * p.k - 1u;
        size_t b = a + p.k;
        float phi = cargf(y[b] * conjf(y[a]));
        float av = 2.0f * (float)nsfsk_prefix_sym(&p, w) + 1.0f - (float)NSFSK_M;
        num += phi * av; den += av * av;
    }
    float sc = num / den;
    fprintf(stderr, "probe12: nf=%u cutoff=%.0f sc=%.6f (ideal %.6f) s0=%ld lag=%ld\n",
            nf, cutoff, sc, (float)M_PI * p.h, s0, lag);
    for (long w = n0; w < n1; w++) {
        size_t a = (size_t)s0 + (size_t)(w + lag) * p.k - 1u;
        size_t b = a + p.k;
        if (b >= N) break;
        float phi = cargf(y[b] * conjf(y[a]));
        float v = phi / sc;
        unsigned s = (v < -2.0f) ? 0u : (v < 0.0f) ? 1u : (v < 2.0f) ? 2u : 3u;
        float amag = cabsf(y[b]);
        printf("%ld exp=%u got=%u v=%+.3f |y|=%.5f\n", w,
               w < (long)D0 ? nsfsk_prefix_sym(&p, (unsigned)w) : 0xffu, s, v, amag);
    }
    return 0;
}
