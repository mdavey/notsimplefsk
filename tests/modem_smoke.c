/*
 * modem_smoke — validate the notsimplefsk PHY loopback:
 *   cpfskmod (TX) -> [AWGN] -> own full-symbol differential detector (RX)
 *
 * usage: modem_smoke [k] [h] [bt] [delay] [pulse]   (defaults 20 0.25 0.8 1 gmsk)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <liquid/liquid.h>

#include "common.h"

static double gauss_noise(double *spare, int *have)
{
    double u1, u2;
    if (*have) { *have = 0; return *spare; }
    do { u1 = drand48(); } while (u1 < 1e-12);
    u2 = drand48();
    *spare = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
    *have  = 1;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* decision on one full-symbol phase increment */
static unsigned detect_phi(float phi, float scale)
{
    float v = phi / scale;              /* ~ a in {-3,-1,+1,+3} */
    return (v < -2.0f) ? 0u : (v < 0.0f) ? 1u : (v < 2.0f) ? 2u : 3u;
}

int main(int argc, char **argv)
{
    unsigned k = (argc > 1) ? (unsigned)atoi(argv[1]) : 20u;
    float    h = (argc > 2) ? (float)atof(argv[2]) : 0.25f;
    float   bt = (argc > 3) ? (float)atof(argv[3]) : 0.8f;
    unsigned m = (argc > 4) ? (unsigned)atoi(argv[4]) : 1u;
    int   pulse = LIQUID_CPFSK_GMSK;
    if (argc > 5 && nsfsk_parse_pulse(argv[5], &pulse)) {
        fprintf(stderr, "bad pulse '%s'\n", argv[5]);
        return 1;
    }
    const unsigned N = 40000u;
    const char *pname = pulse == LIQUID_CPFSK_GMSK ? "gmsk" :
                        pulse == LIQUID_CPFSK_RCOS_FULL ? "rcos" : "square";

    printf("loopback: bps=%u k=%u h=%.3f bt=%.2f delay=%u pulse=%s\n",
           NSFSK_BPS, k, h, bt, m, pname);

    cpfskmod mod = cpfskmod_create(NSFSK_BPS, h, k, m, bt, pulse);
    if (mod == NULL) { fprintf(stderr, "FAILED cpfskmod_create\n"); return 1; }

    unsigned char *tx = (unsigned char *)malloc(N);
    float complex *stream = malloc((size_t)(N + 2) * k * sizeof(float complex));
    float complex *y = (float complex *)malloc(k * sizeof(float complex));
    srand48(12345);
    for (unsigned i = 0; i < N; i++) tx[i] = (unsigned char)(lrand48() & 3u);

    cpfskmod_reset(mod);
    for (unsigned i = 0; i < N; i++) {
        cpfskmod_modulate(mod, tx[i], y);
        memcpy(stream + (size_t)i * k, y, k * sizeof(float complex));
    }

    /* phase increment spans k-1 samples between window edges */
    float scale = (float)M_PI * h * (float)(k - 1u) / (float)k;

    /* ---- clean BER, sweep lag (total loop delay) ---- */
    double best = 1.0; int bl = -1;
    for (unsigned lag = 0; lag <= 2u * m + 3u; lag++) {
        unsigned err = 0, tot = 0;
        for (unsigned i = 2u * m + 2u + lag; (i + 1u) * k <= N * k; i++) {
            float phi = cargf(stream[(size_t)(i + 1) * k - 1] *
                              conjf(stream[(size_t)i * k - 1]));
            if (detect_phi(phi, scale) != tx[i - lag]) err++;
            tot++;
        }
        double ber = (double)err / (double)tot;
        if (ber < best) { best = ber; bl = (int)lag; }
    }
    printf("clean loopback: BER=%.3e at lag=%d  %s\n", best, bl,
           best == 0.0 ? "OK" : "FAIL");

    /* ---- AWGN sweep at best lag ---- */
    double spare = 0.0;
    int have = 0;
    const double ebn0_db[] = { 6.0, 9.0, 12.0, 15.0, 18.0, 21.0 };
    printf("\nAWGN (Es=1, sigma^2/dim = 1/(4*Eb/N0)):\n");
    for (unsigned s = 0; s < sizeof(ebn0_db)/sizeof(ebn0_db[0]); s++) {
        double snr = pow(10.0, ebn0_db[s] / 10.0);
        float  sig = (float)sqrt(1.0 / (4.0 * snr));
        size_t ntot = (size_t)N * k;
        float complex *noz = malloc(ntot * sizeof(float complex));
        have = 0;
        for (size_t j = 0; j < ntot; j++)
            noz[j] = stream[j] +
                sig * ((float)gauss_noise(&spare, &have) +
                       _Complex_I * (float)gauss_noise(&spare, &have));
        unsigned err = 0, tot = 0;
        unsigned ulag = (unsigned)(bl < 0 ? 0 : bl);
        for (unsigned i = 2u * m + 2u + ulag; (i + 1u) * k <= N * k; i++) {
            float phi = cargf(noz[(size_t)(i + 1) * k - 1] *
                              conjf(noz[(size_t)i * k - 1]));
            if (detect_phi(phi, scale) != tx[i - ulag]) err++;
            tot++;
        }
        printf("  Eb/N0=%4.1f dB   BER=%10.6f\n",
               ebn0_db[s], (double)err / (double)tot);
        free(noz);
    }

    free(tx); free(stream); free(y);
    cpfskmod_destroy(mod);
    return best == 0.0 ? 0 : 1;
}
