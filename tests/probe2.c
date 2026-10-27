/* probe2: effective-window-start sweep, single window (no trio), dmis curve */
#include "rx_core.c"

/* single-window tone bank: correlate [start, start+k) interpolated */
static unsigned one_win(const float complex *y, size_t N, double start,
                        unsigned k, const float complex *be, float *q)
{
    long b = (long)floor(start);
    float fr = (float)(start - (double)b);
    if (b < 0 || b + (long)k + 1 > (long)N) return BANK_OOB;
    float complex z[NSFSK_M], u[NSFSK_M];
    float e = 0.0f;
    for (unsigned s = 0u; s < NSFSK_M; s++) { z[s] = 0.0f; u[s] = 1.0f; }
    for (long n = 0; n < (long)k; n++) {
        const float complex v = (1.0f - fr) * y[b + n] + fr * y[b + n + 1L];
        e += crealf(v * conjf(v));
        for (unsigned s = 0u; s < NSFSK_M; s++) z[s] += v * u[s];
        for (unsigned s = 0u; s < NSFSK_M; s++) u[s] *= be[s];
    }
    float eb = -1.0f, gbest = -1.0f;
    unsigned best = 0u;
    for (unsigned s = 0u; s < NSFSK_M; s++) {
        float v = crealf(z[s] * conjf(z[s]));
        if (v > eb) { eb = v; best = s; }
        float g = v / ((float)k * e);
        if (g > gbest) gbest = g;
    }
    if (q) *q = gbest;
    return best;
}

int main(int argc, char **argv)
{
    nsfsk_params p;
    nsfsk_params_default(&p);
    p.h = 0.5f; p.fc = 5000.0f;
    if (argc > 1) p.h = (float)atof(argv[1]);
    if (argc > 2) p.fc = (float)atof(argv[2]);

    size_t cap = 1u << 16, n16 = 0;
    int16_t *in = malloc(cap * sizeof(int16_t));
    for (;;) {
        if (n16 == cap) { cap *= 2; in = realloc(in, cap * sizeof(int16_t)); }
        size_t got = fread(in + n16, sizeof(int16_t), cap - n16, stdin);
        n16 += got; if (!got) break;
    }
    size_t N = n16 - (n16 & 1u);
    unsigned nf = nsfsk_lpf_taps(&p);
    float complex *y = malloc(N * sizeof(float complex));
    if (nsfsk_rx_filter_slice(&p, in, N, y)) return 1;

    uint8_t payload[65536];
    FILE *f = fopen("/tmp/big/p1024.bin", "rb");
    size_t len = fread(payload, 1, sizeof payload, f); fclose(f);
    unsigned msglen = NSFSK_HDR_LEN + len + NSFSK_CRC_LEN;
    uint8_t *msg = calloc(msglen, 1);
    msg[0] = 0; msg[1] = (len >> 8) & 0xff; msg[2] = len & 0xff;
    memcpy(msg + NSFSK_HDR_LEN, payload, len);
    uint8_t crc[4];
    nsfsk_crc32(msg, NSFSK_HDR_LEN + len, crc);
    memcpy(msg + NSFSK_HDR_LEN + len, crc, 4);
    unsigned nsym = 4u * msglen;
    unsigned char *dsym = malloc(nsym);
    nsfsk_data_to_symbols(msg, msglen, dsym);

    const unsigned D0 = p.preamble + NSFSK_SYNC_LEN;
    const float aa[NSFSK_M] = { -3.0f, -1.0f, 1.0f, 3.0f };
    float complex bank_e[NSFSK_M];
    for (unsigned s = 0u; s < NSFSK_M; s++)
        bank_e[s] = cexpf(-I * (float)M_PI * p.h * aa[s] / (float)p.k);

    /* first data symbol's window starts somewhere around sample E0;
     * total frame start ~140 samples + D0 symbols... we don't know s0
     * absolutely, so sweep E0 within [120..220+40] coarsely then fine:
     * E0 = 140 + D0*20 + e, e in -30..+30 step 0.5 */
    double base = 17.0 + (double)D0 * (double)p.k;
    for (double e = -30.0; e <= 30.0; e += 1.0) {
        unsigned dmis = 0; float qsum = 0;
        for (unsigned w = 0; w < 400u; w++) {
            float qf;
            unsigned s = one_win(y, N, base + e + (double)w * (double)p.k,
                                 p.k, bank_e, &qf);
            if (s == BANK_OOB) break;
            if (s != dsym[w]) dmis++;
            qsum += qf;
        }
        printf("e=%+6.1f dmis=%3u/400 q=%.3f\n", e, dmis, qsum / 400.0);
    }
    return 0;
}
