/*
 * rx_core — shared G4FSK receive detector for fsk_rx and fsk_rx_stream.
 *
 * The pieces are deliberately window-shaped: everything operates on a raw
 * s16le window or an already-filtered complex window, so the offline RX
 * (one window = the whole capture) and the streaming RX (one window = one
 * burst plus margins) share the identical signal chain.  A window is
 * self-contained because nsfsk_lpf_taps() rounds the tap count to an
 * integer-symbol group delay and the double pass removes it, and because
 * every consumer of y is invariant to a constant complex rotation (MF
 * magnitude sums, differential phases), so a slice's DDC needs no phase
 * continuity with what came before it.
 *
 * Design rationale: docs/reference.md.
 */
#ifndef NSFSK_RX_CORE_H
#define NSFSK_RX_CORE_H

#include "common.h"

#include <complex.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/* Decoded frames appended to buf; ok/bad count RS+CRC gate passes/failures;
 * qual/ndec accumulate per-decision quality for ppm trial selection. */
/* Per-frame record so the offline RX can MERGE frames across ppm trials:
 * each frame's sensitivity to residual sample-rate error is several
 * samples shifted from its neighbours (per-burst MF refit luck + slip),
 * so no single global ppm is guaranteed to decode every frame; a
 * CRC-verified frame from any trial can be delivered.  The record list
 * grows on demand (an arbitrary capture can contain any number of
 * frames); duplicates across trials are identified by rescaled position
 * + payload in fsk_rx.c, never by a position window alone.
 * fix = wire bytes RS had to correct (0 = the block arrived clean), so the
 * merged frame list can also report how much of the delivery is FEC work. */
typedef struct { size_t s0, off; unsigned len, fix; } rx_frame;

typedef struct {
    unsigned char *buf; size_t cap, nbytes;
    unsigned ok, bad; double qual; size_t ndec;
    /* FEC accounting: ok_rs counts accepted frames that only cleared the
     * RS+CRC gate because RS corrected something, rs_bits is the wire bits
     * RS fixed in accepted frames (8 per corrected byte).  lost[] holds the
     * start sample of every COMPLETE codeword the RS/CRC gate rejected
     * (candidates pruned on the 8-symbol header peek, or cut off by the end
     * of the window, are 'bad' but were never a whole block, so they are not
     * lost frames); the offline RX dedupes those positions across ppm trials
     * the same way it dedupes delivered frames. */
    unsigned ok_rs;
    unsigned long long rs_bits;
    unsigned acq;              /* expensive candidate acquisitions */
    rx_frame *fr; unsigned nfr, fr_cap;
    size_t *lost; unsigned nlost, lost_cap;
    /* -vv diagnostics: per-symbol bank quality of the last decoded frame */
    float dbg_q[4200]; unsigned dbg_n;
} rx_result;

/* Append n bytes to R->buf (n may be 0).  Returns 0 on success, -1 on
 * allocation failure (R is left unchanged). */
int rx_result_append(rx_result *R, const unsigned char *b, size_t n);

/* Release everything owned by R and zero it (safe on a zeroed/empty
 * result and safe to call twice). */
void rx_result_free(rx_result *R);

/* Shared FEC report line pair (stderr in both receivers): frames attempted /
 * recovered / lost and how many recovered frames needed RS, then the same
 * in wire bits (recovered + unrecovered = total, percentages, plus the bits
 * RS actually fixed).  ok = accepted frames, ok_rs = accepted frames that
 * needed RS, lost = rejected complete codewords, rs_bits = bits RS fixed. */
void rx_fec_report(FILE *f, unsigned ok, unsigned ok_rs, unsigned lost,
                   unsigned long long rs_bits);

/* complex linear-interpolation resampler: out[i] = in[i * ratio] */
float complex *resample_lin(const float complex *x, size_t n,
                            double ratio, size_t *m);

/* baseband cutoff of the receive LPF (dev + 0.85*Rs) */
float nsfsk_rx_cutoff(const nsfsk_params *p);

/* DDC at fc + zero-phase double-pass Kaiser LPF: raw s16le window
 * in[0..n-1] -> y[0..n-1] (caller allocates y with n entries).
 * Returns 0, or -1 when fc cannot reject the mixer image (message printed
 * on stderr). */
int nsfsk_rx_filter_slice(const nsfsk_params *p,
                          const int16_t *in, size_t n, float complex *y);

/* One detector pass over (possibly resampled) y: energy gate, sub-window
 * matched filter, (timing offset, detector lag, CFO) prefix sweep, and
 * differential decode with CRC gating; payloads appended to R.
 * reject_limit caps consecutive-ish prefix-misrejects (a mid-join garbage
 * burst otherwise costs a full MF+prefix sweep every 4k samples); at the
 * cap the pass disarms and idles to the next gap, so a later real frame
 * in the same window still decodes.
 * max_acq_total caps expensive acquisitions per pass (0 = unlimited); the
 * allowance grows by 8 per CRC-good frame so a productive pass can finish.
 * max_acq_burst caps acquisitions in one energy burst; at the cap the pass
 * idles to the next silence and re-arms, preserving a real frame after a
 * garbage burst while stopping repeated dips in a continuous tone/reverb. */
void rx_run(const float complex *y, size_t N, unsigned nf,
            const float complex *tmpl, float et,
            const nsfsk_params p, int verbose, unsigned reject_limit,
            unsigned max_acq_total, unsigned max_acq_burst,
            rx_result *R);

#endif /* NSFSK_RX_CORE_H */
