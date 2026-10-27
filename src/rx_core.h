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

/* Decoded frames appended to buf; ok/bad count CRC passes/failures;
 * qual/ndec accumulate per-decision quality for ppm trial selection. */
typedef struct {
    unsigned char *buf; size_t cap, nbytes;
    unsigned ok, bad; double qual; size_t ndec;
} rx_result;

void rx_result_append(rx_result *R, const unsigned char *b, size_t n);

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
 * in the same window still decodes.  0 = unlimited (offline RX). */
void rx_run(const float complex *y, size_t N, unsigned nf,
            const float complex *tmpl, float et,
            const nsfsk_params p, int verbose, unsigned reject_limit,
            rx_result *R);

#endif /* NSFSK_RX_CORE_H */
