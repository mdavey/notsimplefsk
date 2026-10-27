/*
 * probe_lpf — sweep kaiser LPF (cutoff, length) through the real RX
 * front-end chain: cpfskmod -> upconvert(real) -> DDC -> kaiser LPF ->
 * full-symbol differential decisions.  Reports clean BER per config.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <liquid/liquid.h>

#include "common.h"

static unsigned detect_phi(float phi, float scale)
{
    float v = phi / scale;
    return (v < -2.0f) ? 0u : (v < 0.0f) ? 1u : (v < 2.0f) ? 2u : 3u;
}

int main(void)
{
    nsfsk_params p;
    nsfsk_params_default(&p);
    const unsigned Nsym = 3000u;
    const size_t   N    = (size_t)Nsym * p.k;
    const float fs = (float)p.rate;

    cpfskmod mod = cpfskmod_create(NSFSK_BPS, p.h, p.k, p.delay, p.bt, p.pulse);
    float complex *bb = malloc(N * sizeof(float complex));
    float complex *tmp = malloc(p.k * sizeof(float complex));
    float complex *y  = malloc(N * sizeof(float complex));
    unsigned char *txs = malloc(Nsym);
    srand48(99);
    for (unsigned i = 0; i < Nsym; i++) txs[i] = (unsigned char)(lrand48() & 3u);
    for (unsigned i = 0; i < Nsym; i++) {
        cpfskmod_modulate(mod, txs[i], tmp);
        memcpy(bb + (size_t)i * p.k, tmp, p.k * sizeof(float complex));
    }

    /* DDC result of the real passband signal: half the baseband plus the
     * negative-frequency image at -2*fc */
    for (size_t n = 0; n < N; n++) {
        float complex img = 0.5f * conjf(bb[n]) *
            cexpf(-_Complex_I * 4.0f * (float)M_PI * p.fc * (float)n / fs);
        y[n] = 0.5f * bb[n] + img;
    }

    float scale = (float)M_PI * p.h * (float)(p.k - 1u) / (float)p.k;
    printf("%8s %5s %9s %10s\n", "cutoff", "n", "BER", "delay");
    float cutoffs[] = {1200.f, 1600.f, 2000.f, 2400.f, 2800.f, 3200.f};
    /* only lengths whose group delay is an integer number of symbols,
     * so decisions stay at symbol boundaries */
    unsigned lens[] = {41u, 81u, 121u, 161u, 201u, 281u};
    for (unsigned ci = 0; ci < sizeof(cutoffs)/sizeof(cutoffs[0]); ci++) {
        for (unsigned li = 0; li < sizeof(lens)/sizeof(lens[0]); li++) {
            if (((lens[li] - 1u) / 2u) % p.k) continue;
            firfilt_crcf lpf = firfilt_crcf_create_kaiser(
                lens[li], cutoffs[ci] / fs, 80.0f, 0.0f);
            const float *hbuf = firfilt_crcf_get_coefficients(lpf);
            unsigned err = 0, tot = 0;
            /* steady region after filter delay */
            size_t d = lens[li] / 2u;
            unsigned a = p.k ? (unsigned)(d / p.k) : 0u;   /* symbol delay */
            for (unsigned i = (unsigned)(lens[li] / p.k) + 2u * p.delay + a + 4u;
                 (i + 1u) * p.k <= N; i++) {
                float complex fa = 0, fb = 0;
                for (unsigned t = 0; t < lens[li]; t++) {
                    size_t i0 = (size_t)i * p.k - 1u, i1 = (size_t)(i+1u)*p.k - 1u;
                    if (i0 >= t) fa += hbuf[t] * y[i0 - t];
                    if (i1 >= t) fb += hbuf[t] * y[i1 - t];
                }
                float phi = cargf(fb * conjf(fa));
                if (detect_phi(phi, scale) != txs[i - 1u - a]) err++;
                tot++;
            }
            printf("%8.0f %5u %9.2e %10zu\n", cutoffs[ci], lens[li],
                   (double)err / (double)tot, d);
            firfilt_crcf_destroy(lpf);
        }
    }
    return 0;
}
