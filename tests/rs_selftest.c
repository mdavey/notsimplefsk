/*
 * rs_selftest — validate src/rs.c on all four parity tiers.
 *
 * Pins the properties the frame codec relies on:
 *   - clean blocks decode as clean (0 corrections);
 *   - any e <= t byte errors (every single position, and random bursts of
 *     0..t) are always FULLY corrected (data and parity restored);
 *   - > t errors NEVER silently reproduce the original codeword: the
 *     decoder returns -1, or corrects to a DIFFERENT (valid) codeword,
 *     which the frame CRC then rejects.  Return count is always <= t.
 *   - the shortened code over the padded region behaves like the full one
 *     (errors sweep across the data/parity boundary).
 */
#include "rs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned long rng_state = 0x2545F491uL;

static unsigned rnd(unsigned n)
{
    /* xorshift32-ish LCG mix, deterministic across platforms */
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return (unsigned)(rng_state >> 16) % n;
}

static unsigned nfail;

static void fill(uint8_t *b, unsigned n)
{
    for (unsigned i = 0; i < n; i++) b[i] = (uint8_t)rnd(256u);
}

/* inject ne errors at DISTINCT positions, returns ne */
static unsigned inject(uint8_t *blk, unsigned ne)
{
    static unsigned used[NSFSK_RS_N];
    static unsigned stamp;
    unsigned done = 0;

    stamp++;
    while (done < ne) {
        unsigned pos = rnd(NSFSK_RS_N);
        if (used[pos] == stamp) continue;
        used[pos] = stamp;
        blk[pos] ^= (uint8_t)(rnd(255u) + 1u);
        done++;
    }
    return ne;
}

static void trial_ok(unsigned P, unsigned ne)
{
    unsigned k = NSFSK_RS_N - P;
    uint8_t data[255], blk[255], out[255], par[NSFSK_RS_PMAX];

    fill(data, k);
    nsfsk_rs_encode(P, data, par);
    memcpy(blk, data, k);
    memcpy(blk + k, par, P);
    if (ne > 0) inject(blk, ne);

    int ret = nsfsk_rs_decode(P, blk, out);
    if (ret != (int)ne) {
        printf("FAIL P=%u e=%u: ret=%d\n", P, ne, ret);
        nfail++;
        return;
    }
    if (memcmp(out, data, k) || memcmp(blk + k, par, P)) {
        printf("FAIL P=%u e=%u: not fully corrected\n", P, ne);
        nfail++;
    }
}

static void trial_over(unsigned P, unsigned ne)
{
    unsigned k = NSFSK_RS_N - P, t = P / 2u;
    uint8_t data[255], blk[255], out[255], par[NSFSK_RS_PMAX];

    fill(data, k);
    nsfsk_rs_encode(P, data, par);
    memcpy(blk, data, k);
    memcpy(blk + k, par, P);
    inject(blk, ne);

    int ret = nsfsk_rs_decode(P, blk, out);
    if (ret >= 0) {
        if ((unsigned)ret > t) {
            printf("FAIL P=%u e=%u: ret=%d exceeds t\n", P, ne, ret);
            nfail++;
        }
        /* > t errors: landing back on the ORIGINAL codeword would mean the
         * decoder silently passed a wrong frame on to the CRC */
        if (memcmp(out, data, k) == 0 && memcmp(blk + k, par, P) == 0) {
            printf("FAIL P=%u e=%u: miscorrected back to original\n", P, ne);
            nfail++;
        }
    }
}

int main(void)
{
    static const unsigned tiers[4] = { 16u, 32u, 64u, 84u };

    for (unsigned ti = 0; ti < 4u; ti++) {
        unsigned P = tiers[ti], t = P / 2u, k = NSFSK_RS_N - P;
        uint8_t data[255], blk[255], out[255], par[NSFSK_RS_PMAX];

        /* single-error sweep across the whole 255-byte block: pins every
         * position mapping (data region, shortened/padded region, parity) */
        for (unsigned pos = 0; pos < NSFSK_RS_N; pos++) {
            fill(data, k);
            nsfsk_rs_encode(P, data, par);
            memcpy(blk, data, k);
            memcpy(blk + k, par, P);
            blk[pos] ^= (uint8_t)(rnd(255u) + 1u);
            int ret = nsfsk_rs_decode(P, blk, out);
            if (ret != 1 || memcmp(out, data, k) ||
                memcmp(blk + k, par, P)) {
                printf("FAIL P=%u pos=%u: single error not corrected\n",
                       P, pos);
                nfail++;
            }
        }

        /* clean blocks */
        for (unsigned n = 0; n < 20u; n++) trial_ok(P, 0);
        /* every correctable error count, incl. exactly t */
        for (unsigned e = 1; e <= t; e++)
            for (unsigned n = 0; n < 40u; n++) trial_ok(P, e);
        /* beyond t: never silently wrong */
        for (unsigned e = t + 1u; e <= t + 8u; e++)
            for (unsigned n = 0; n < 50u; n++) trial_over(P, e);

        printf("rs tier P=%-2u t=%-2u: ok (failures so far %u)\n",
               P, t, nfail);
    }

    if (nfail) {
        printf("rs_selftest: %u FAILURE(S)\n", nfail);
        return 1;
    }
    printf("rs_selftest: all tiers passed\n");
    return 0;
}
