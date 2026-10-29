/* probe: find which alignment/offset lets bank_decide reproduce wire syms */
#include "rx_core.c"

struct R { double mis; double off; double fr; };

static int cmp_d (const void *a, const void *b) {
    return ((const struct R*)a)->mis - ((const struct R*)b)->mis;
}

int main(int argc, char **argv)
{
    nsfsk_params p;
    nsfsk_params_default(&p);
    p.h = 0.5f; p.fc = 5000.0f;
    if (argc > 1) p.h = (float)atof(argv[1]);
    if (argc > 2) p.fc = (float)atof(argv[2]);

    /* payload */
    size_t cap = 1u << 16, n16 = 0;
    int16_t *in = malloc(cap * sizeof(int16_t));
    if (in == NULL) { fprintf(stderr, "probe_bank: out of memory\n"); return 1; }
    for (;;) {
        if (n16 == cap) {
            int16_t *ni = realloc(in, cap * 2u * sizeof(int16_t));
            if (ni == NULL) { fprintf(stderr, "probe_bank: out of memory\n"); free(in); return 1; }
            in = ni; cap *= 2;
        }
        size_t got = fread(in + n16, sizeof(int16_t), cap - n16, stdin);
        n16 += got; if (!got) break;
    }
    size_t N = n16;
    float complex *y = malloc(N * sizeof(float complex));
    if (y == NULL) { fprintf(stderr, "probe_bank: out of memory\n"); return 1; }
    if (nsfsk_rx_filter_slice(&p, in, N, y)) return 1;

    /* true wire: msg = [type,len16] + payload + crc */
    uint8_t payload[65536];
    FILE *f = fopen("/tmp/big/p1024.bin", "rb");
    if (f == NULL) { fprintf(stderr, "probe_bank: cannot open /tmp/big/p1024.bin\n"); return 1; }
    size_t len = fread(payload, 1, sizeof payload, f); fclose(f);
    unsigned msglen = NSFSK_HDR_LEN + (unsigned)len + NSFSK_CRC_LEN;
    uint8_t *msg = calloc(msglen, 1);
    if (msg == NULL) { fprintf(stderr, "probe_bank: out of memory\n"); return 1; }
    msg[0] = 0; msg[1] = (len >> 8) & 0xff; msg[2] = len & 0xff;
    memcpy(msg + NSFSK_HDR_LEN, payload, len);
    uint8_t crc[4];
    nsfsk_crc32(msg, NSFSK_HDR_LEN + len, crc);
    memcpy(msg + NSFSK_HDR_LEN + len, crc, 4);
    unsigned nsym = 4u * msglen;
    unsigned char *dsym = malloc(nsym);
    if (dsym == NULL) { fprintf(stderr, "probe_bank: out of memory\n"); return 1; }
    nsfsk_data_to_symbols(msg, msglen, dsym);

    const unsigned D0 = p.preamble + NSFSK_SYNC_LEN;
    const float aa[NSFSK_M] = { -3.0f, -1.0f, 1.0f, 3.0f };
    float complex bank_e[NSFSK_M];
    for (unsigned s = 0u; s < NSFSK_M; s++)
        bank_e[s] = cexpf(-I * (float)M_PI * p.h * aa[s] / (float)p.k);

    /* true frame waveform starts at s0_true (preamble first sym sample).
     * sweep integer s0c 0..40 and fractional window off -10..+10 */
    struct R res[41 * 21]; int n = 0;
    for (long s0c = 0; s0c <= 40; s0c++) {
        for (long off = -10; off <= 10; off++) {
            double pos0 = (double)s0c + (double)off;
            unsigned mis = 0, tot = 0;
            for (unsigned w = NSFSK_TMPL_SYM0(p); w < D0; w++) {
                unsigned s = bank_decide(y, N, pos0 + (double)w * (double)p.k,
                                         p.k, bank_e, NULL);
                if (s == BANK_OOB) break;
                tot++;
                if (s != nsfsk_prefix_sym(&p, w)) mis++;
            }
            if (tot < 60u) continue;
            /* data compare (first 200 syms for speed) */
            unsigned dmis = 0;
            for (unsigned w = 0; w < 200u; w++) {
                unsigned s = bank_decide(y, N,
                    pos0 + (double)(D0 + w) * (double)p.k, p.k, bank_e, NULL);
                if (s == BANK_OOB) break;
                if (s != dsym[w]) {
                    dmis++;
                    if (s0c == 13 && off == 10)
                        fprintf(stderr,
                            "ERR w=%u got=%u want=%u prev=%u pprev=%u\n",
                            w, s, dsym[w],
                            w ? dsym[w-1] : 255, w > 1 ? dsym[w-2] : 255);
                }
            }
            res[n++] = (struct R){ (double)mis + (double)dmis, (double)off, (double)s0c };
            if (mis == 0 || dmis == 0)
                fprintf(stderr, "s0c=%ld off=%+ld: pmis=%u dmis=%u/200\n",
                        s0c, off, mis, dmis);
            if (mis <= 3 && (s0c + off >= 21 && s0c + off <= 25))
                fprintf(stderr, "cell s0c=%ld off=%+ld: pmis=%u dmis=%u\n",
                        s0c, off, mis, dmis);
        }
    }
    qsort(res, n, sizeof res[0], cmp_d);
    for (int i = 0; i < n && i < 8; i++)
        printf("best: pmis+dmis=%2.0f s0c=%2.0f off=%+.0f\n",
               res[i].mis, res[i].fr, res[i].off);
    return 0;
}
