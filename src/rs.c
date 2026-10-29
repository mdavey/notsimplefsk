#include "rs.h"

#include <stdlib.h>
#include <string.h>

/* GF(256) with primitive polynomial x^8 + x^4 + x^3 + x^2 + 1 (0x11d). */
static uint8_t gf_exp[512];   /* alpha^i; doubled so a+b never needs mod 255 */
static uint8_t gf_log[256];

/* Generator polynomials for the four tiers, computed once:
 * g(x) = prod_{i=1..P} (x + alpha^i);  rs_gen[t][j] = coeff of x^(P-1-j)
 * (descending order, leading 1 omitted) -> direct LFSR encoder taps. */
static const unsigned rs_tiers[4] = { 16u, 32u, 64u, 84u };
static uint8_t rs_gen[4][NSFSK_RS_PMAX];
static int rs_ready = 0;

static uint8_t gmul(uint8_t a, uint8_t b)
{
    return (a && b) ? gf_exp[gf_log[a] + gf_log[b]] : 0;
}

static uint8_t ginv(uint8_t a)
{
    return gf_exp[255u - gf_log[a]];
}

static void rs_init(void)
{
    unsigned i, t;
    int x = 1;

    for (i = 0; i < 255u; i++) {
        gf_exp[i]  = (uint8_t)x;
        gf_log[x]  = (uint8_t)i;
        x <<= 1;
        if (x & 0x100) x ^= 0x11d;
    }
    for (i = 255u; i < 512u; i++) gf_exp[i] = gf_exp[i - 255u];

    for (t = 0; t < 4u; t++) {
        unsigned P = rs_tiers[t], d;
        uint8_t ga[NSFSK_RS_PMAX + 1];      /* ascending: ga[j] = x^j coeff */

        ga[0] = 1;
        for (d = 1; d <= P; d++) ga[d] = 0;
        for (i = 1; i <= P; i++) {          /* multiply by (x + alpha^i) */
            uint8_t a = gf_exp[i];
            for (d = P; d >= 1; d--)
                ga[d] = (uint8_t)(ga[d - 1] ^ gmul(ga[d], a));
            ga[0] = gmul(ga[0], a);
        }
        for (i = 0; i < P; i++) rs_gen[t][i] = ga[P - 1 - i];
    }
    rs_ready = 1;
}

static const uint8_t *rs_genpoly(unsigned P)
{
    unsigned t;
    if (!rs_ready) rs_init();
    for (t = 0; t < 4u; t++)
        if (rs_tiers[t] == P) return rs_gen[t];
    return NULL;
}

void nsfsk_rs_encode(unsigned P, const uint8_t *data, uint8_t *parity)
{
    const uint8_t *gd = rs_genpoly(P);
    unsigned k = NSFSK_RS_N - P, i, j;

    if (gd == NULL) abort();                /* tier not validated */

    memset(parity, 0, P);
    for (i = 0; i < k; i++) {               /* LFSR: parity = data*x^P mod g */
        uint8_t fb = data[i] ^ parity[0];
        memmove(parity, parity + 1, P - 1u);
        parity[P - 1u] = 0;
        for (j = 0; j < P; j++) parity[j] ^= gmul(fb, gd[j]);
    }
}

/* S[1..P] over block[0..254] (block[j] = coeff of x^(254-j)) */
static void syndromes(const uint8_t *block, unsigned P, uint8_t *S)
{
    unsigned i, j;
    for (i = 1; i <= P; i++) {
        uint8_t a = gf_exp[i], v = 0;
        for (j = 0; j < NSFSK_RS_N; j++)
            v = (uint8_t)(gmul(v, a) ^ block[j]);
        S[i] = v;
    }
}

/* Horner evaluation of an ascending coefficient array (c[0..deg]) at x */
static uint8_t geval(const uint8_t *c, unsigned deg, uint8_t x)
{
    uint8_t v = 0;
    unsigned i;
    for (i = deg + 1u; i-- > 0u; )
        v = (uint8_t)(gmul(v, x) ^ c[i]);
    return v;
}

int nsfsk_rs_decode(unsigned P, uint8_t *block, uint8_t *out_data)
{
    unsigned t = P / 2u, k = NSFSK_RS_N - P, i, j;
    uint8_t S[NSFSK_RS_PMAX + 1];
    uint8_t C[NSFSK_RS_PMAX + 1], B[NSFSK_RS_PMAX + 1], T[NSFSK_RS_PMAX + 1];
    uint8_t Om[NSFSK_RS_PMAX + 1], Dp[NSFSK_RS_PMAX + 1];
    unsigned L, m, ncor;

    if (rs_genpoly(P) == NULL) return -1;

    syndromes(block, P, S);
    for (i = 1; i <= P; i++)
        if (S[i]) break;
    if (i > P) {                            /* clean block */
        if (out_data != block) memcpy(out_data, block, k);
        return 0;
    }

    /* Berlekamp-Massey: C(x) = error locator (ascending), L = its degree */
    memset(C, 0, sizeof C);
    memset(B, 0, sizeof B);
    C[0] = 1;
    B[0] = 1;
    L = 0;
    m = 1;
    {
        uint8_t b = 1;
        long N;
        for (N = 0; N < (long)P; N++) {
            uint8_t d = S[N + 1], coef;
            for (i = 1; i <= L; i++) d ^= gmul(C[i], S[N + 1 - i]);
            if (d == 0) { m++; continue; }
            memcpy(T, C, sizeof C);
            coef = gmul(d, ginv(b));
            for (i = 0; i + m <= P; i++) C[i + m] ^= gmul(coef, B[i]);
            if (2u * L <= (unsigned)N) {
                memcpy(B, T, sizeof B);
                b = d;
                L = (unsigned)N + 1u - L;
                m = 1;
            } else {
                m++;
            }
        }
    }
    if (L > t) return -1;                   /* needs more than t corrections */

    /* Omega(x) = S(x) C(x) mod x^P, S(x) = sum_{i=1..P} S[i] x^(i-1) */
    memset(Om, 0, sizeof Om);
    for (i = 0; i < P; i++)
        for (j = 0; i + j < P; j++)
            Om[i + j] ^= gmul(S[i + 1], C[j]);

    /* formal derivative C'(x) (odd-degree terms survive in characteristic 2) */
    memset(Dp, 0, sizeof Dp);
    for (i = 1; i <= L; i += 2) Dp[i - 1] = C[i];

    /* Chien search + Forney magnitudes: root at alpha^-j (exponent j) means
     * an error in block[254 - j]; e = Omega(alpha^-j)/C'(alpha^-j) */
    ncor = 0;
    for (j = 0; j < NSFSK_RS_N; j++) {
        uint8_t xi = gf_exp[255u - j];      /* alpha^-j */
        uint8_t e, dv;
        if (geval(C, L, xi) != 0) continue;
        dv = geval(Dp, L - 1u > 0u ? L - 1u : 0u, xi);
        if (dv == 0) return -1;             /* degenerate locator */
        e = gmul(geval(Om, P - 1u, xi), ginv(dv));
        block[NSFSK_RS_N - 1u - j] ^= e;
        ncor++;
    }
    if (ncor != L) return -1;               /* repeated/missing roots */

    /* never emit an unverified correction */
    syndromes(block, P, S);
    for (i = 1; i <= P; i++)
        if (S[i]) return -1;

    if (out_data != block) memcpy(out_data, block, k);
    return (int)ncor;
}
