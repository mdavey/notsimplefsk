/*
 * notsimplefsk — shared definitions for the G4FSK modem pair (fsk_tx / fsk_rx)
 *
 * PHY: 4-level continuous-phase FSK with Gaussian pulse shaping (G4FSK),
 *      TX via liquid-dsp cpfskmod (bps=2, LIQUID_CPFSK_GMSK by default).
 *      RX implements its own full-symbol differential detector
 *      (arg(y[t+k]·conj(y[t])) -> nearest pi*h*{+/-1,+/-3}); liquid's
 *      cpfskdem is NOT used: its M>2 non-coherent path in this build is
 *      unreliable (partial-response ISI, unfixed filter betas).
 *      NOTE: liquid's GMSK pulse spans ~3 symbols at any m; simple
 *      differential detection needs bt >= ~0.8 (default).  bt=rcos/square
 *      pulses work down to any bt (full response).
 *
 * Symbol -> frequency (relative to carrier fc):
 *      f(s) = (2s + 1 - M) * h * Rs / 2  [Hz]  M = 4,  Rs = rate/k
 *      h=0.25 -> tones at {-0.375,-0.125,+0.125,+0.375} * Rs
 *
 *      Non-coherent (differential) demod requires principal-value
 *      resolution of the per-symbol phase step: pi*h*(M-1) < pi,
 *      i.e. h < 1/(M-1) = 1/3 for 4FSK.  Default h = 0.25.
 *      Adjacent tone spacing = h * Rs.
 *
 * Frame (on the wire):
 *      [ preamble: 0,3,0,3,... x P ][ sync: 24-sym shift-distinguishing table ]
 *      [ Gray( NRZI( PN( type | len_hi | len_lo | payload | crc32_be ) ) ) ]
 *      [ pad: 0,3,... x nsfsk_tx_pad()  -- detector lag + RX LPF edge reach ]
 *
 *  - type  : 1 byte, arbitrary
 *  - len   : payload length, 16-bit big-endian, 0..NSFSK_MAX_PAYLOAD
 *  - crc32 : LIQUID_CRC_32 over type|len|payload, big-endian trailer
 *  - PN    : additive whitening with a 2^23-1 m-sequence (x^23+x^18+1,
 *            all-ones seed, re-armed per frame at the data section) so
 *            payload runs (zeros, repeated bytes) never become single-tone
 *            symbol runs on the wire; descramble = same local PN, no error
 *            multiplication, prefix-deterministic for the header-first decode
 *  - NRZI  : NRZ-M over the (whitened) bit stream, MSB-first, state = 0 at
 *            first bit:
 *              enc[n] = enc[n-1] XOR data[n]   (transition on 1)
 *            cancels 180-degree phase flips of the CPM waveform
 *  - Gray  : dibit -> symbol index via gray(v) = v ^ (v>>1) = {0,1,3,2}
 */
#ifndef NSFSK_COMMON_H
#define NSFSK_COMMON_H

#include <stdint.h>
#include <complex.h>

#define NSFSK_BPS           2u
#define NSFSK_M             (1u << NSFSK_BPS)
#define NSFSK_MAX_PAYLOAD   1024u
#define NSFSK_HDR_LEN       3u                 /* type + len16 */
#define NSFSK_CRC_LEN       4u

#define NSFSK_SYNC_REPS     6u
#define NSFSK_SYNC_LEN      (4u * NSFSK_SYNC_REPS)   /* 24 symbols */

/* Upper bound on the preamble, enforced by nsfsk_params_validate().  Keeps
 * every frame-size computation well inside unsigned/size_t range and stops
 * "-P -1"/huge values from turning into multi-GB allocations or OOB
 * pointer arithmetic. */
#define NSFSK_MAX_PREAMBLE  65536u

/* matched-filter template: 32 symbols taken from the steady preamble,
 * starting at prefix symbol (delay + 4) so the Gaussian filter is settled */
#define NSFSK_TMPL_SYMS     32u
#define NSFSK_TMPL_SYM0(p)  ((p).delay + 4u)

/* worst-case wire symbols per frame (size_t so the arithmetic cannot wrap) */
#define NSFSK_MAX_FRAME_SYMS(p) \
    ((size_t)(p).preamble + NSFSK_SYNC_LEN \
     + 4u * (size_t)(NSFSK_HDR_LEN + NSFSK_MAX_PAYLOAD + NSFSK_CRC_LEN))

typedef struct {
    unsigned rate;        /* sample rate, Hz                              */
    float    fc;          /* carrier frequency, Hz                        */
    unsigned k;           /* samples per symbol (even, 4..160)            */
    float    h;           /* CPM modulation index (> 0)                   */
    float    bt;          /* Gaussian bandwidth-time product, (0,1]       */
    unsigned delay;       /* filter delay, symbols 1..8                   */
    int      pulse;       /* LIQUID_CPFSK_GMSK (default) | RCOS_FULL | SQUARE */
    unsigned preamble;    /* preamble length, symbols                     */
} nsfsk_params;

void          nsfsk_params_default(nsfsk_params *p);
/* returns NULL when valid, else a short reason string; prints warnings */
const char   *nsfsk_params_validate(const nsfsk_params *p);
/* parse --pulse argument: "gmsk" | "rcos" | "square"; returns 0 on success */
int           nsfsk_parse_pulse(const char *s, int *out);

extern const unsigned char nsfsk_gray[NSFSK_M];         /* dibit -> sym   */
extern const unsigned char nsfsk_gray_inv[NSFSK_M];     /* sym -> dibit   */
extern const unsigned char nsfsk_sync_sym[NSFSK_SYNC_LEN];

/* wire symbol at index i of the raw prefix (preamble + sync) */
unsigned char nsfsk_prefix_sym(const nsfsk_params *p, unsigned i);

/* bytes -> wire symbols over the data section: nsym = 4*nbytes */
void nsfsk_data_to_symbols(const uint8_t *data, unsigned nbytes,
                           unsigned char *sym);
/* wire symbols -> bytes over the data section: nbytes = nsym/4 */
void nsfsk_symbols_to_data(const unsigned char *sym, unsigned nsym,
                           uint8_t *data);

/* CRC-32 (liquid LIQUID_CRC_32) over buf[0..n-1], big-endian out */
void nsfsk_crc32(const uint8_t *buf, unsigned n, uint8_t out[4]);

/* Build a complete frame's wire symbols (prefix + data section).
 * Returns number of symbols written, or 0 on error. */
unsigned nsfsk_build_frame_syms(const nsfsk_params *p, uint8_t type,
                                const uint8_t *payload, unsigned plen,
                                unsigned char *out, unsigned out_cap);

/* Generate the complex baseband waveform of prefix symbols [sym0, sym0+nsyms)
 * exactly as the transmitter emits them at frame start (cpfskmod fresh from
 * reset).  Returns number of samples written (nsyms*k). */
unsigned nsfsk_build_template(const nsfsk_params *p, unsigned sym0,
                              unsigned nsyms, float complex *out);

/* RX front-end LPF length (kaiser taps, group delay (nf-1)/2): shared rule
 * so the TX trailing pad can cover the RX filter edge reach. */
unsigned nsfsk_lpf_taps(const nsfsk_params *p);

/* TX trailing pad in symbols: covers detector lag (lag_max+1), the RX
 * zero-phase LPF edge reach (~(nf-1)/2 samples), and a margin. */
unsigned nsfsk_tx_pad(const nsfsk_params *p);

#endif /* NSFSK_COMMON_H */
