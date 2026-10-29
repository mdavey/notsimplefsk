#include "common.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <liquid/liquid.h>

const unsigned char nsfsk_gray[NSFSK_M]     = { 0, 1, 3, 2 };
const unsigned char nsfsk_gray_inv[NSFSK_M] = { 0, 1, 3, 2 };

/* non-periodic by design: >=15/24 symbol mismatches vs any shift 1..8,
 * which resolves the alignment ambiguity created by the periodic 0,3
 * preamble (and by the period-4 sync of earlier drafts) */
const unsigned char nsfsk_sync_sym[NSFSK_SYNC_LEN] = {
    3,1,1,0, 2,1,0,0, 1,3,2,3, 3,3,3,1, 1,0,0,0, 1,2,0,3
};

void nsfsk_params_default(nsfsk_params *p)
{
    p->rate     = 48000u;
    p->fc       = 4000.0f;
    p->k        = 20u;
    p->h        = 0.25f;
    p->bt       = 0.8f;
    p->delay    = 1u;
    p->pulse    = LIQUID_CPFSK_GMSK;
    p->preamble = 72u;
    p->fec      = NSFSK_FEC_DEFAULT;
}

unsigned nsfsk_max_payload(const nsfsk_params *p)
{
    return NSFSK_WIRE_BYTES - NSFSK_HDR_LEN - NSFSK_CRC_LEN - p->fec;
}

int nsfsk_parse_pulse(const char *s, int *out)
{
    if      (!strcmp(s, "gmsk"))   *out = LIQUID_CPFSK_GMSK;
    else if (!strcmp(s, "rcos"))   *out = LIQUID_CPFSK_RCOS_FULL;
    else if (!strcmp(s, "square")) *out = LIQUID_CPFSK_SQUARE;
    else return -1;
    return 0;
}

const char *nsfsk_params_validate(const nsfsk_params *p)
{
    if (p->rate < 8000u)                   return "rate must be >= 8000";
    if (p->k < 4u || (p->k % 2u) || p->k > 160u)
        return "k must be even and within 4..160";
    if (!(p->h > 0.0f))                    return "h must be > 0";
    if (!(p->bt > 0.0f) || p->bt > 1.0f)   return "bt must be in (0,1]";
    if (p->delay == 0u || p->delay > 8u)   return "delay must be within 1..8";
    if (p->preamble > NSFSK_MAX_PREAMBLE)
        return "preamble too long (max 65536 symbols)";
    if (p->fec != 16u && p->fec != 32u && p->fec != 64u && p->fec != 84u)
        return "fec must be one of 16, 32, 64, 84";
    if (p->preamble < NSFSK_TMPL_SYM0(*p) + NSFSK_TMPL_SYMS)
        return "preamble too short for sync template";
    if (!(p->fc > 0.0f) || p->fc >= 0.49f * (float)p->rate)
        return "fc must be in (0, 0.5*rate)";

    /* h*(M-1) < 1 decodes with the per-symbol phase-advance detector;
     * >= 1 switches to the non-coherent tone bank (h >= 1 fully
     * orthogonal).  So h is bounded only by what the RX DDC+LPF can fit:
     * cutoff = dev + 0.85*Rs must stay clear of DC (mixer image at
     * -2*fc) on the low side and of Nyquist (image at rate-2*fc) on the
     * high side — mirrors nsfsk_rx_filter_slice()'s image check. */
    float rs  = (float)p->rate / (float)p->k;
    float fd  = p->h * rs * (float)(NSFSK_M - 1u) / 2.0f;  /* peak deviation */
    float occ = 2.0f * (fd + rs);           /* Carson occupied bandwidth */
    float cut = fd + 0.85f * rs;            /* RX LPF cutoff */
    if (p->fc <= cut)
        return "fc too low for h/k: RX LPF cannot pass the FM band and "
               "reject the mixer image (raise fc or lower h)";
    if (0.5f * (float)p->rate - p->fc <= cut)
        return "fc too high for h/k: RX LPF cannot pass the FM band and "
               "reject the upper mixer image (lower fc or h)";
    if (p->fc - 0.5f * occ < 100.0f || p->fc + 0.5f * occ > 0.5f * (float)p->rate - 100.0f)
        fprintf(stderr,
                "warning: occupied band %.0f..%.0f Hz is close to the band edges\n",
                p->fc - 0.5f * occ, p->fc + 0.5f * occ);
    if (p->pulse == LIQUID_CPFSK_GMSK && p->bt < 0.5f)
        fprintf(stderr,
                "warning: bt=%.3f leaves an ISI error floor for the simple "
                "differential detector (0.8 recommended; <0.5 needs MLSE)\n",
                p->bt);
    return NULL;
}

unsigned char nsfsk_prefix_sym(const nsfsk_params *p, unsigned i)
{
    if (i < p->preamble)
        return (i & 1u) ? 3u : 0u;
    i -= p->preamble;
    if (i < NSFSK_SYNC_LEN)
        return nsfsk_sync_sym[i];
    return 0u;
}

/* Payload whitening: additive m-sequence scrambler on the data bits before
 * NRZI.  Without it the payload statistics reach the wire verbatim: NRZI
 * only inserts transitions on 1 bits, so constant data (0x00, 0xFF, a
 * repeated byte, JPEG run tables) holds the dibit fixed and the TX emits
 * one tone for hundreds of symbols -- measured: 640 zero bytes = 2560
 * identical symbols, MtWellington.jpg max run 51.  Additive (local PN, not
 * feedback of the transmitted bits) means zero error multiplication and
 * bit-position determinism, which rx_core's header-first prefix decode
 * requires: every nsfsk_symbols_to_data() call re-arms the PN exactly as
 * the TX did at the start of the data section, both sides implicitly
 * resyncing on the sync table.  PN = Fibonacci LFSR with primitive
 * x^23 + x^18 + 1 (ADSL G.992 lineage, period verified 2^23-1), all-ones
 * seed; the period (2^23-1 bits) far exceeds the fixed 1020-symbol data
 * section, so one run of the sequence covers any frame and single-symbol
 * runs are bounded by the m-sequence's own max run (~23 bits -> ~12
 * symbols). */
#define PN_MASK 0x7fffffu
#define PN_SEED 0x7fffffu

static unsigned pn_bit(unsigned *reg)
{
    unsigned out = (*reg >> 22) & 1u;
    *reg = ((*reg << 1) | (out ^ ((*reg >> 17) & 1u))) & PN_MASK;
    return out;
}

void nsfsk_data_to_symbols(const uint8_t *data, unsigned nbytes,
                           unsigned char *sym)
{
    unsigned nsym = 0, i;
    unsigned e_prev = 0, acc = 0, nb = 0, pn = PN_SEED;

    for (i = 0; i < nbytes; i++) {
        int b;
        for (b = 7; b >= 0; b--) {
            unsigned d = (data[i] >> b) & 1u;
            unsigned e = e_prev ^ d ^ pn_bit(&pn);
            e_prev = e;
            acc = (acc << 1) | e;
            if (++nb == 2u) {
                sym[nsym++] = nsfsk_gray[acc];
                acc = 0;
                nb  = 0;
            }
        }
    }
}

void nsfsk_symbols_to_data(const unsigned char *sym, unsigned nsym,
                           uint8_t *data)
{
    unsigned nbyte = 0, i;
    unsigned e_prev = 0, acc = 0, nbits = 0, pn = PN_SEED;

    for (i = 0; i < nsym; i++) {
        unsigned v = nsfsk_gray_inv[sym[i] & (NSFSK_M - 1u)];
        int b;
        for (b = 1; b >= 0; b--) {
            unsigned e = (v >> b) & 1u;
            unsigned d = e ^ e_prev ^ pn_bit(&pn);
            e_prev = e;
            acc = (acc << 1) | d;
            if (++nbits == 8u) {
                data[nbyte++] = (uint8_t)acc;
                acc   = 0;
                nbits = 0;
            }
        }
    }
}

void nsfsk_crc32(const uint8_t *buf, unsigned n, uint8_t out[4])
{
    unsigned int crc = crc_generate_key(LIQUID_CRC_32,
                                        (unsigned char *)(uintptr_t)buf, n);
    out[0] = (uint8_t)(crc >> 24);
    out[1] = (uint8_t)(crc >> 16);
    out[2] = (uint8_t)(crc >>  8);
    out[3] = (uint8_t)(crc      );
}

void nsfsk_frame_encode(const nsfsk_params *p, uint8_t type,
                        const uint8_t *payload, unsigned plen, uint8_t *wire)
{
    unsigned k = NSFSK_WIRE_BYTES - p->fec;

    if (plen > nsfsk_max_payload(p)) plen = nsfsk_max_payload(p);
    wire[0] = type;
    wire[1] = (uint8_t)plen;
    if (plen) memcpy(wire + NSFSK_HDR_LEN, payload, plen);
    nsfsk_crc32(wire, NSFSK_HDR_LEN + plen, wire + NSFSK_HDR_LEN + plen);
    memset(wire + NSFSK_HDR_LEN + plen + NSFSK_CRC_LEN, 0,
           k - (NSFSK_HDR_LEN + plen + NSFSK_CRC_LEN));
    nsfsk_rs_encode(p->fec, wire, wire + k);
}

int nsfsk_frame_decode(const nsfsk_params *p, const uint8_t *wire,
                       uint8_t *type, uint8_t *payload, unsigned *plen,
                       unsigned *fixed)
{
    uint8_t blk[NSFSK_WIRE_BYTES], crc[4];
    unsigned len, maxp = nsfsk_max_payload(p);
    int cor;

    memcpy(blk, wire, NSFSK_WIRE_BYTES);
    cor = nsfsk_rs_decode(p->fec, blk, blk);
    if (cor < 0) {                      /* uncorrectable: reject outright */
        if (fixed != NULL) *fixed = 0;
        return 0;
    }
    if (fixed != NULL) *fixed = (unsigned)cor;
    len = blk[1];
    if (len > maxp)
        return 0;
    nsfsk_crc32(blk, NSFSK_HDR_LEN + len, crc);
    if (memcmp(crc, blk + NSFSK_HDR_LEN + len, NSFSK_CRC_LEN) != 0)
        return 0;
    *type = blk[0];
    *plen = len;
    if (len) memcpy(payload, blk + NSFSK_HDR_LEN, len);
    return 1;
}

unsigned nsfsk_build_frame_syms(const nsfsk_params *p, uint8_t type,
                                const uint8_t *payload, unsigned plen,
                                unsigned char *out, unsigned out_cap)
{
    uint8_t wire[NSFSK_WIRE_BYTES];
    size_t need, i;

    if (plen > nsfsk_max_payload(p))
        return 0;

    need = (size_t)p->preamble + NSFSK_SYNC_LEN + NSFSK_WIRE_SYMS;
    if (need > (size_t)out_cap || need > UINT_MAX)
        return 0;

    nsfsk_frame_encode(p, type, payload, plen, wire);

    for (i = 0; i < (size_t)p->preamble + NSFSK_SYNC_LEN; i++)
        out[i] = nsfsk_prefix_sym(p, (unsigned)i);

    nsfsk_data_to_symbols(wire, NSFSK_WIRE_BYTES,
                          out + (size_t)p->preamble + NSFSK_SYNC_LEN);
    return (unsigned)need;
}

unsigned nsfsk_lpf_taps(const nsfsk_params *p)
{
    float rs     = (float)p->rate / (float)p->k;
    float dev    = p->h * rs * (float)(NSFSK_M - 1u) / 2.0f;
    float cutoff = dev + 0.85f * rs;
    float img    = 2.0f * p->fc - cutoff;
    float trans  = img - cutoff;                 /* stop-band room */
    float trans2 = 2.0f * (cutoff - 1.3f * dev); /* flat-passband room */
    if (trans <= 0.0f) trans = trans2;           /* degenerate: RX rejects fc anyway */
    float tmin = fminf(trans, trans2);
    if (!(tmin > 1.0f)) tmin = 1.0f;

    /* kaiser As~80: transition ~ 4.5*fs/n; round so group delay = int
     * symbols.  Clamp the ideal length BEFORE converting to unsigned:
     * a degenerate fc (fc only infinitesimally above the cutoff) gives a
     * tiny transition and would otherwise make the float->unsigned cast
     * undefined/wrapped. */
    double nf0 = 4.5 * (double)p->rate / (double)tmin;
    if (!(nf0 > 1.0)) nf0 = 1.0;
    if (nf0 > 2401.0) nf0 = 2401.0;
    unsigned nf = (unsigned)((nf0 - 1.0) / (2.0 * (double)p->k) + 1.0)
                * 2u * p->k + 1u;
    if (nf < 2u * p->k + 1u) nf = 2u * p->k + 1u;
    if (nf > 2401u) nf = 2401u;
    nf |= 1u;
    return nf;
}

unsigned nsfsk_tx_pad(const nsfsk_params *p)
{
    /* the zero-phase LPF double pass reads (nf-1)/2 samples past every
     * sample; wire symbols closer to the buffer end than that are filtered
     * against (nonexistent) zeros and their amplitudes decay into the
     * trailer.  Pad so the last real symbol plus the detector's worst-case
     * lag still sits a half-filter-length clear of the end. */
    unsigned nf  = nsfsk_lpf_taps(p);
    unsigned reach = (nf - 1u) / 2u;                    /* samples */
    unsigned lag_max = 2u * p->delay + 3u;              /* rx lag sweep */
    return (reach + p->k - 1u) / p->k + lag_max + 3u;
}

unsigned nsfsk_build_template(const nsfsk_params *p, unsigned sym0,
                              unsigned nsyms, float complex *out)
{
    cpfskmod mod = cpfskmod_create(NSFSK_BPS, p->h, p->k, p->delay,
                                   p->bt, p->pulse);
    float complex *y;
    unsigned i;

    if (mod == NULL) {
        fprintf(stderr, "error: cpfskmod_create() failed (see liquid errors)\n");
        return 0;
    }

    y = (float complex *)malloc(p->k * sizeof(float complex));
    if (y == NULL) {
        fprintf(stderr, "error: out of memory\n");
        cpfskmod_destroy(mod);
        return 0;
    }
    for (i = 0; i < sym0 + nsyms; i++) {
        cpfskmod_modulate(mod, nsfsk_prefix_sym(p, i), y);
        if (i >= sym0)
            memcpy(out + (size_t)(i - sym0) * p->k, y,
                   p->k * sizeof(float complex));
    }

    free(y);
    cpfskmod_destroy(mod);
    return nsyms * p->k;
}
