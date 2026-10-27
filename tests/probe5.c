/* probe5: full (dn, lag) prefix-calibration table like rx_core's alignment
 * sweep, annotated with true data-symbol errors for each candidate.
 * Usage: probe5 rate fc k h bt delay payfile frameidx < in.s16
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <liquid/liquid.h>
#include "common.h"
#include "rx_core.h"

static float wrap_pi(float a)
{
    while (a > (float)M_PI)  a -= 2.0f * M_PI;
    while (a < -(float)M_PI) a += 2.0f * M_PI;
    return a;
}
static float unwrap_ph(float ph)
{
    if      (ph >  (float)M_PI - 0.03f) ph -= 2.0f * M_PI;
    else if (ph < -(float)M_PI + 0.03f) ph += 2.0f * M_PI;
    return ph;
}
static unsigned detect_phi(float phi, float scale)
{
    float v = phi / scale;
    return (v < -2.0f) ? 0u : (v < 0.0f) ? 1u : (v < 2.0f) ? 2u : 3u;
}

int main(int argc, char **argv)
{
    if (argc < 8) { fprintf(stderr, "usage: probe5 rate fc k h bt delay payfile frameidx < in.s16\n"); return 1; }
    nsfsk_params p; nsfsk_params_default(&p);
    p.rate  = (unsigned)strtoul(argv[1], 0, 0);
    p.fc    = strtof(argv[2], 0);
    p.k     = (unsigned)strtoul(argv[3], 0, 0);
    p.h     = strtof(argv[4], 0);
    p.bt    = strtof(argv[5], 0);
    p.delay = (unsigned)strtoul(argv[6], 0, 0);
    long fidx = atol(argv[8]);

    size_t pcap = 1u << 20; unsigned char *pay = malloc(pcap);
    FILE *f = fopen(argv[7], "rb");
    size_t npay = fread(pay, 1, pcap, f); fclose(f);
    unsigned gap = p.rate / 25u;
    unsigned cap_syms = NSFSK_MAX_FRAME_SYMS(p);
    unsigned char *syms = malloc(cap_syms);
    long soff = 0; size_t off = 0; long fi = 0; unsigned nsym = 0;
    while (off < npay || fi == 0) {
        unsigned plen = (unsigned)(npay - off);
        if (plen > NSFSK_MAX_PAYLOAD) plen = NSFSK_MAX_PAYLOAD;
        nsym = nsfsk_build_frame_syms(&p, 0, pay + off, plen, syms, cap_syms);
        if (fi == fidx) break;
        soff += (long)(nsym + nsfsk_tx_pad(&p)) * (long)p.k + (long)gap;
        off += plen; fi++;
    }

    size_t cap = 1u << 22, n16 = 0;
    int16_t *in = malloc(cap * sizeof(int16_t));
    for (;;) {
        if (n16 == cap) { cap *= 2; in = realloc(in, cap * sizeof(int16_t)); }
        size_t got = fread(in + n16, sizeof(int16_t), cap - n16, stdin);
        n16 += got; if (!got) break;
    }
    size_t N = n16 - (n16 & 1u);
    float complex *y = malloc(N * sizeof(float complex));
    nsfsk_rx_filter_slice(&p, in, N, y);

    const unsigned D0 = p.preamble + NSFSK_SYNC_LEN;
    const float scale = (float)M_PI * p.h * (float)(p.k - 1u) / (float)p.k;
    fprintf(stderr, "frame %ld true s0=%ld nsym=%u\n", fidx, soff, nsym);

    /* RX-selected best_n was 86680 for frame 2 => s0_ref = best_n - (delay+4)*k.
     * Sweep dn from -220..80 around the same s0_ref. */
    long s0_ref = soff + 140 - (long)(p.delay + 4u) * (long)p.k;  /* mimics MF peak +140 */
    for (long dn = -220; dn <= 80; dn += 1) {
        for (unsigned lag = 0; lag <= 2u * p.delay + 3u; lag++) {
            long s0 = s0_ref + dn;
            if (s0 < 1) continue;
            float num = 0, den = 0; unsigned mis = 0, tot = 0;
            static float phi_v[512]; static float a_v[512];
            for (unsigned w = p.delay + 4u; w < D0; w++) {
                size_t hi = (size_t)s0 + (size_t)(w + lag + 1u) * p.k - 1u;
                size_t lo = (size_t)s0 + (size_t)(w + lag) * p.k - 1u;
                if (hi >= N) break;
                float ph = unwrap_ph(cargf(y[hi] * conjf(y[lo])));
                float a = 2.0f * (float)nsfsk_prefix_sym(&p, w) + 1.0f - (float)NSFSK_M;
                phi_v[tot] = ph; a_v[tot] = a; tot++;
                num += ph * a; den += a * a;
            }
            float sc = den > 0 ? num / den : 0;
            for (unsigned q = 0; q < tot; q++)
                if (detect_phi(phi_v[q], sc) !=
                    nsfsk_prefix_sym(&p, p.delay + 4u + q)) mis++;
            if (mis > 2u || sc <= 0.1f) continue;
            /* data errors at this candidate: fixed calibrated scale + 1 adaptive pass */
            float dphi = 0, fsc = sc;
            unsigned dmis = 0;
            for (unsigned got = 0; D0 + got < nsym; got++) {
                size_t hi = (size_t)s0 + (size_t)(D0 + got + lag + 1u) * p.k - 1u;
                size_t lo = (size_t)s0 + (size_t)(D0 + got + lag) * p.k - 1u;
                if (hi >= N) break;
                float ph = unwrap_ph(cargf(y[hi] * conjf(y[lo])) - dphi);
                unsigned s = detect_phi(ph, fsc);
                float a_hat = 2.0f * (float)s + 1.0f - (float)NSFSK_M;
                float err = wrap_pi(ph - a_hat * fsc);
                dphi += 0.01f * err;
                fsc  += 0.005f * err * a_hat / 3.0f;
                if (s != syms[D0 + got]) dmis++;
            }
            fprintf(stderr, "dn=%+4ld s0=%ld lag=%u mis=%u sc=%.4f  dmis=%u\n",
                    dn, s0, lag, mis, sc, dmis);
        }
    }
    return 0;
}
