/*
 * rx_core — shared G4FSK receive detector (see rx_core.h).
 * Moved verbatim from fsk_rx.c; behavior must stay bit-identical there.
 */
#include "rx_core.h"

#include <complex.h>
#include <liquid/liquid.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* same detector as modem_smoke: phi over k-1 samples between window edges */
static unsigned detect_phi(float phi, float scale)
{
    float v = phi / scale;              /* ~ a in {-3,-1,+1,+3} */
    return (v < -2.0f) ? 0u : (v < 0.0f) ? 1u : (v < 2.0f) ? 2u : 3u;
}

static float wrap_pi(float a)
{
    while (a > (float)M_PI)  a -= 2.0f * (float)M_PI;
    while (a < -(float)M_PI) a += 2.0f * (float)M_PI;
    return a;
}

/* With a CFO whose sign matches the tone, the outer a=+-3 symbols reach
 * |phi| = 3*unit + |dphi_cfo| ~ 2.9 rad and ISI pushes some past +-pi,
 * where cargf wraps them onto the opposite cluster (+-3 flips).  After the
 * CFO hypothesis is subtracted, legitimate ph never exceeds ~+-2.7 rad,
 * while a wrapped one lands beyond pi + |dphi_cfo| - excess, so a +-2pi
 * correction above pi-0.03 restores it. */
static float unwrap_ph(float ph)
{
    if      (ph >  (float)M_PI - 0.03f) ph -= 2.0f * (float)M_PI;
    else if (ph < -(float)M_PI + 0.03f) ph += 2.0f * (float)M_PI;
    return ph;
}

/* decision at boundary i of the frame (samples s0+i*k-1, s0+(i+1)*k-1),
 * raw phase increment or NAN when out of bounds */
static float boundary_phi(const float complex *y, size_t N, size_t s0,
                          unsigned k, unsigned i)
{
    size_t hi = s0 + (size_t)(i + 1u) * k - 1u;
    size_t lo = s0 + (size_t)i * k - 1u;
    if (hi >= N) return NAN;
    return cargf(y[hi] * conjf(y[lo]));
}

/* Non-coherent tone-bank detector, used when h*(M-1) >= 0.95 where the
 * phase-advance detector structurally cannot work (the outer +-3 tones
 * rotate past +-pi per symbol and wrap onto the inner cluster).
 * Decide the symbol whose head sits at fractional sample position pos:
 * correlate the k-sample window starting at pos (linear-interpolated
 * reads for fractional pos) against the M tone frequencies a*PI*h/k per
 * sample, a in {-3,-1,+1,+3}; pick the peak bin and report its best
 * normalized coherent gain |z|^2 / (k * window energy) ~ 1 when the
 * window is aligned with a tone (all k samples add coherently), falling
 * as boundaries mix symbols.  NOT the energy share: that stays high for
 * any window containing ANY strong tone, so it cannot rank ppm trials.
 * *q carries that gain (ppm-trial quality metric in place of
 * ph/fsc*a_hat).  BANK_OOB when the window is out of bounds.
 * Timing is NOT recovered here.  ±k/4 micro-window schemes (peak-bin,
 * sum-energy NDA, with or without feedback on the picked offset) were
 * all tried and removed: inside a tone run every candidate window scores
 * identically, so any per-symbol window selection flips at noise and the
 * flipped window straddles a symbol boundary = coin-toss decisions on
 * later symbols at earlier indices (sporadic outer<->inner tone errors,
 * high qmean, CRC death).  The static lattice (bpos += k, resampled per
 * ppm trial) plus per-frame union across trials (fsk_rx.c merge) is the
 * timing architecture; measured dmis=0/400 on clean frames at h=0.5,
 * 0.667 and 1.0 with this single-window detector.
 * Textbook non-coherent M-FSK decide, fully orthogonal at h=1; CFO needs
 * no correction (magnitude decision; Dirichlet null sits h*Rs away). */
#define BANK_OOB 0xFFFFFFFFu
static unsigned bank_decide(const float complex *y, size_t N, double pos,
                            unsigned k, const float complex *be, float *q)
{
    long base = (long)floor(pos);
    float  fr = (float)(pos - (double)base);
    if (base < 0L || base + (long)k + 1L > (long)N) return BANK_OOB;
    float complex z[NSFSK_M], u[NSFSK_M];
    float e = 0.0f;
    for (unsigned s = 0u; s < NSFSK_M; s++) { z[s] = 0.0f; u[s] = 1.0f; }
    for (long n = 0; n < (long)k; n++) {
        const float complex v = (1.0f - fr) * y[base + n]
                              + fr * y[base + n + 1L];
        e += crealf(v * conjf(v));
        for (unsigned s = 0u; s < NSFSK_M; s++) z[s] += v * u[s];
        for (unsigned s = 0u; s < NSFSK_M; s++) u[s] *= be[s];
    }
    if (!(e > 0.0f)) return BANK_OOB;
    float eb = -1.0f, gbest = -1.0f;
    unsigned best = 0u;
    for (unsigned s = 0u; s < NSFSK_M; s++) {
        float v = crealf(z[s] * conjf(z[s]));
        if (v > eb) { eb = v; best = s; }
        float g = v / ((float)k * e);
        if (g > gbest) gbest = g;
    }
    if (q) *q = gbest > 0.0f ? gbest : 0.0f;
    return best;
}

/* complex linear-interpolation resampler: out[i] = in[i * ratio].
 * Returns NULL with *m = 0 for degenerate input (n == 0, non-finite or
 * non-positive ratio) instead of producing an invalid size_t cast. */
float complex *resample_lin(const float complex *x, size_t n,
                            double ratio, size_t *m)
{
    if (n == 0 || !(ratio > 0.0) || !isfinite(ratio)) { *m = 0; return NULL; }
    size_t mo = (size_t)(((double)n - 1.0) / ratio) + 1u;
    float complex *o = malloc((mo + 4) * sizeof(float complex));
    if (o == NULL) { *m = 0; return NULL; }
    for (size_t i = 0; i < mo; i++) {
        double pos = (double)i * ratio;
        if (pos > (double)(n - 1)) { o[i] = 0.0f; continue; }
        size_t i0 = (size_t)pos;
        double f = pos - (double)i0;
        size_t i1 = (i0 + 1u < n) ? i0 + 1u : i0;
        o[i] = (float)(1.0 - f) * x[i0] + (float)f * x[i1];
    }
    *m = mo;
    return o;
}

int rx_result_append(rx_result *R, const unsigned char *b, size_t n)
{
    if (n == 0) return 0;
    if (R->nbytes + n > R->cap) {
        size_t cap = R->cap ? R->cap : 65536u;
        while (cap < R->nbytes + n) {
            if (cap > (SIZE_MAX / 2u)) return -1;
            cap *= 2u;
        }
        unsigned char *nb = realloc(R->buf, cap);
        if (nb == NULL) return -1;
        R->buf = nb;
        R->cap = cap;
    }
    memcpy(R->buf + R->nbytes, b, n);
    R->nbytes += n;
    return 0;
}

void rx_result_free(rx_result *R)
{
    free(R->buf);
    free(R->fr);
    free(R->lost);
    memset(R, 0, sizeof *R);
}

/* Record a rejected complete codeword's position (stats only: an allocation
 * failure just undercounts, the RX keeps running). */
static void rx_result_note_lost(rx_result *R, size_t s0)
{
    if (R->nlost == R->lost_cap) {
        unsigned ncap = R->lost_cap ? R->lost_cap * 2u : 64u;
        size_t *nl = realloc(R->lost, (size_t)ncap * sizeof *nl);
        if (nl == NULL) return;
        R->lost = nl;
        R->lost_cap = ncap;
    }
    R->lost[R->nlost++] = s0;
}

void rx_fec_report(FILE *f, unsigned ok, unsigned ok_rs, unsigned lost,
                   unsigned long long rs_bits)
{
    /* One data frame is exactly one shortened RS(255, 255-P) codeword, so
     * the bit totals are frame counts times the wire width; a recovered
     * frame is CRC-verified and therefore counted error-free even when RS
     * had to fix bits inside it (those are the RS-corrected subtotal). */
    const unsigned long long bf  = (unsigned long long)NSFSK_WIRE_BYTES * 8ull;
    const unsigned long long rec = (unsigned long long)ok * bf;
    const unsigned long long tot = (unsigned long long)(ok + lost) * bf;
    const unsigned long long unre = tot - rec;

#define PCT(x, d) ((d) ? 100.0 * (double)(x) / (double)(d) : 0.0)
    fprintf(f, "rx: FEC frames: %u attempted, %u recovered, %u lost, "
               "%u (%.2f%% of recovered) recovered via RS\n",
            ok + lost, ok, lost, ok_rs, PCT(ok_rs, ok));
    fprintf(f, "rx: FEC bits: %llu total, %llu recovered (%.2f%%), "
               "%llu unrecovered (%.2f%%), %llu bit(s)/frame, "
               "%llu RS-corrected (%.2f%% of total)\n",
            tot, rec, PCT(rec, tot), unre, PCT(unre, tot), bf,
            rs_bits, PCT(rs_bits, tot));
#undef PCT
}

float nsfsk_rx_cutoff(const nsfsk_params *p)
{
    float rs  = (float)p->rate / (float)p->k;
    float dev = p->h * rs * (float)(NSFSK_M - 1u) / 2.0f;      /* peak Hz dev */
    return dev + 0.85f * rs;
}

int nsfsk_rx_filter_slice(const nsfsk_params *p,
                          const int16_t *in, size_t n, float complex *y)
{
    if (n == 0) return 0;
    /* ---- DDC + LPF ----
     * After mixing down by fc the negative-frequency image (conjugate
     * baseband) lands at -2*fc; the LPF must pass the whole FM band
     * (dev + shaping, ~0.65*Rs skirt) yet stop that image, and its group
     * delay MUST be an integer number of symbols or the fixed decision
     * instants sample mid-symbol and read a blend of two tones
     * (verified in tests/probe_lpf.c: <2400 Hz cutoff or non-integer
     * symbol delay wrecks BER entirely). */
    float cutoff = nsfsk_rx_cutoff(p);
    float img    = 2.0f * p->fc - cutoff;                      /* image near edge */
    if (img <= cutoff) {
        fprintf(stderr, "error: fc=%.0f too low to reject the mixer image "
                        "for these PHY params\n", p->fc);
        return -1;
    }
    /* taps from the shared rule (also sizes the TX trailing pad):  1) stop
     * the image,  2) passband flat to ~1.3*dev over the FM tones */
    unsigned nf = nsfsk_lpf_taps(p);

    nco_crcf     nco = nco_crcf_create(LIQUID_NCO);
    firfilt_crcf lpf = firfilt_crcf_create_kaiser(nf, cutoff / (float)p->rate,
                                                  80.0f, 0.0f);
    if (nco == NULL || lpf == NULL) {
        fprintf(stderr, "error: RX filter allocation failed\n");
        if (nco) nco_crcf_destroy(nco);
        if (lpf) firfilt_crcf_destroy(lpf);
        return -1;
    }
    nco_crcf_set_frequency(nco, 2.0f * (float)M_PI * p->fc / (float)p->rate);

    for (size_t i = 0; i < n; i++) {
        float complex z, f;
        nco_crcf_mix_down(nco, (float complex)((float)in[i] / 32768.0f), &z);
        nco_crcf_step(nco);
        firfilt_crcf_execute_one(lpf, z, &f);
        y[i] = f;
    }
    /* zero-phase second pass (filtfilt-style): |H|^2 magnitude with NO group
     * delay.  A causal pass would shift the frame by (nf-1)/2 samples,
     * requiring the TX pad to exceed that delay (it does not for wide LPFs
     * at small k) and forcing the delay to an integer symbol count; both
     * constraints vanish with the double pass. */
    firfilt_crcf_reset(lpf);
    for (size_t i = n; i-- > 0; ) {
        float complex f;
        firfilt_crcf_execute_one(lpf, y[i], &f);
        y[i] = f;
    }

    nco_crcf_destroy(nco);
    firfilt_crcf_destroy(lpf);
    return 0;
}

/* ==================== one detector pass over (possibly resampled) y === */
/* Incoherent sub-window matched-filter alignment score of the template at
 * y[n..n+tn), same metric as the burst search (et division omitted: only
 * used for relative thresholds). */
static float mf_score(const float complex *y, const float complex *tmpl,
                      size_t n, unsigned sub, unsigned nwsub, float et,
                      float complex *rr)
{
    float accm = 0.0f, e = 0.0f;
    for (unsigned w = 0; w < nwsub; w++) {
        float complex r = 0.0f;
        for (unsigned t = w * sub; t < (w + 1u) * sub; t++) {
            float yr = crealf(y[n + t]), yi = cimagf(y[n + t]);
            r += tmpl[t] * (yr - I * yi);
            e += yr * yr + yi * yi;
        }
        if (rr) rr[w] = r;
        accm += cabsf(r);
    }
    return e > 0.0f ? (accm * accm) / (e * et) : 0.0f;
}

void rx_run(const float complex *y, size_t N, unsigned nf,
            const float complex *tmpl, float et,
            const nsfsk_params p, int verbose, unsigned reject_limit,
            unsigned max_acq_total, unsigned max_acq_burst,
            rx_result *R)
{
    /* ---- energy profile (moving average over win samples) ---- */
    if (N == 0) return;
    unsigned win = (unsigned)fminf((float)(p.rate / 50),
                                   0.5f * (float)(p.preamble * p.k));
    if (win < p.k) win = p.k;
    float *prof = malloc(N * sizeof(float));
    if (prof == NULL) { fprintf(stderr, "rx: out of memory\n"); return; }
    double acc = 0.0;
    for (size_t i = 0; i < N; i++) {
        acc += cabsf(y[i]) * cabsf(y[i]);
        if (i >= win) acc -= cabsf(y[i - win]) * cabsf(y[i - win]);
        prof[i] = (float)((i + 1u < win) ? acc / (double)(i + 1u)
                                         : acc / (double)win);
    }
    float pmax = 0.0f;
    for (size_t i = 0; i < N; i++) if (prof[i] > pmax) pmax = prof[i];
    if (!(pmax > 0.0f)) { free(prof); return; }

    const unsigned tn = NSFSK_TMPL_SYMS * p.k;
    const unsigned D0      = p.preamble + NSFSK_SYNC_LEN;  /* first data sym */
    const size_t   max_syn = NSFSK_WIRE_SYMS; /* fixed 1020 data symbols    */
    const unsigned max_pl  = nsfsk_max_payload(&p);
    unsigned char *wire    = malloc((size_t)D0 + max_syn); /* full frame syms */
    uint8_t  *blk  = malloc(NSFSK_WIRE_BYTES);  /* descrambled codeword     */
    uint8_t  *pay  = malloc(NSFSK_WIRE_BYTES);  /* delivered payload        */
    if (wire == NULL || blk == NULL || pay == NULL) {
        fprintf(stderr, "rx: out of memory\n");
        free(prof); free(wire); free(blk); free(pay);
        return;
    }
    /* The decision window is k samples wide (boundary_phi compares
     * samples i*k-1 and (i+1)*k-1), so the nominal phase advance per
     * unit symbol value is pi*h, NOT pi*h*(k-1)/k.  This is only a
     * prior: the prefix sweep calibrates the real scale per candidate. */
    const float scale = (float)M_PI * p.h;
    const unsigned lag_max = 2u * p.delay + 3u;
    /* h*(M-1) >= 0.95: outer tone phase advances within 0.05*pi of the
     * cargf wrap (unwrap_ph shaves legitimate readings; h=0.33 was already
     * losing the +-3 tones there) -> switch to the tone bank */
    const unsigned bank = (p.h * (float)(NSFSK_M - 1u) >= 0.95f);
    float complex bank_e[NSFSK_M];
    {
        const float aa[NSFSK_M] = { -3.0f, -1.0f, 1.0f, 3.0f };
        for (unsigned s = 0u; s < NSFSK_M; s++)
            bank_e[s] = cexpf(-I * (float)M_PI * p.h * aa[s] / (float)p.k);
    }
    const unsigned silence_arm = p.rate / 50u;  /* 20 ms; TX gap is 40 ms */
    const float ethresh = 0.15f * pmax;

    float *cand_phi = malloc((size_t)D0 * sizeof(float));
    float *cand_a   = malloc((size_t)D0 * sizeof(float));
    if (cand_phi == NULL || cand_a == NULL) {
        fprintf(stderr, "rx: out of memory\n");
        free(prof); free(wire); free(blk); free(pay);
        free(cand_phi); free(cand_a);
        return;
    }

    /* dn-sweep span (frame-start hypotheses around the MF peak) and the
     * prefilter score array over it; see SWEEP_DN below */
    const long  dn_lo  = -(long)(p.preamble * p.k + win + 2u * nf);
    const long  dn_hi  = 4L * (long)p.k;
    const size_t dn_span = (size_t)(dn_hi - dn_lo) + 1u;
    float *dn_sc = malloc(dn_span * sizeof(float));
    if (dn_sc == NULL) {
        fprintf(stderr, "rx: out of memory\n");
        free(prof); free(wire); free(blk); free(pay);
        free(cand_phi); free(cand_a);
        return;
    }
    long   cnd[1024];

    int armed = 1;
    unsigned below = silence_arm;                /* start as if silent */
    unsigned rejects = 0;
    unsigned burst_acq = 0;
    size_t i = win < N / 2u ? win : N / 2u;

    while (i + (size_t)(p.preamble * p.k) < N) {
        if (!armed) {
            if (prof[i] < ethresh) {
                if (++below >= silence_arm) {
                    armed = 1; below = 0; burst_acq = 0;
                }
            } else below = 0;
            i++;
            continue;
        }
        if (prof[i] < ethresh) { below = 0; i++; continue; }
        /* work budgets: the total cap stops this ppm trial; the burst cap
         * idles to the next silence so a real later frame can re-arm. */
        if (max_acq_total && R->acq >= max_acq_total + 8u * R->ok)
            break;
        if (max_acq_burst && burst_acq >= max_acq_burst) {
            armed = 0; below = 0; i++;
            continue;
        }
        if (verbose > 1)
            fprintf(stderr, "rx: burst candidate @%zu (prof=%.3g pmax=%.3g)\n",
                    i, prof[i], pmax);

        /* --- burst at energy rise i: incoherent sub-window matched filter
         * (4 x TMPL_SYMS/4 symbols; each sub-window nulls at ~300 Hz so the
         *  summed magnitude tracks alignment for any in-tolerance CFO) --- */
        /* the energy gate lags the burst onset by up to win samples (plus
         * the LPF delay), so the MF window must reach back that far */
        size_t back = (size_t)win + 4u * p.k;
        size_t n_lo = i > back ? i - back : 0u;
        size_t n_hi = i + (size_t)(p.preamble + 8u) * p.k;
        if (n_hi > N - tn) n_hi = N - tn;
        const unsigned nwsub = 8u;                 /* 4-symbol sub-windows: */
        const unsigned sub = tn / nwsub;           /* Dirichlet null at fs/sub */
        /* MF scan stride: must stay under the SCORE-LOBE half width, which
         * the periodic 0,3 preamble shrinks to ~k/5 samples at large h
         * (h=1.0 k=20 measured: peak +-4, ripple plateau out to +-20 - a
         * stride-8 grid aliased the real peak away and gate-skipped frame 1,
         * 2026-10-28).  Widest safe power of two with MS/2 inside the lobe:
         * MS = floor_pow2(k/4), capped at 8 for the large-k CFO case. */
        unsigned MS = 1u;
        while (MS < 8u && (MS << 1u) * 4u <= p.k) MS <<= 1u;
        size_t best_n = 0;
        float  best_sc = 0.0f;
        float complex best_rr[8], rr[8];
        {
            /* strided scan + full-resolution refine bracketing the stride
             * max: MS/2 stays inside the score lobe (rule above), so the
             * stride max lands in the winning lobe and refine resolves
             * the exact peak */
            size_t nb = n_lo; float nbv = 0.0f;
            for (size_t n = n_lo; n + tn <= n_hi; n += MS) {
                float sc = mf_score(y, tmpl, n, sub, nwsub, et, NULL);
                if (sc > nbv) { nbv = sc; nb = n; }
            }
            size_t r0 = (nb > n_lo + MS / 2u) ? nb - MS / 2u : n_lo;
            size_t r1 = nb + MS / 2u;
            for (size_t n = r0; n + tn <= n_hi && n <= r1; n++) {
                float sc = mf_score(y, tmpl, n, sub, nwsub, et, rr);
                if (sc > best_sc) { best_sc = sc; best_n = n;
                                    memcpy(best_rr, rr, sizeof rr); }
            }
        }
        if (verbose > 1)
            fprintf(stderr, "rx: mf best_n=%zu score=%.3f cfo_stage0\n",
                    best_n, best_sc);
        if (best_sc < 0.25f) {                      /* no real frame here */
            i += (size_t)(p.preamble * p.k);
            armed = 0; below = 0;
            continue;
        }
        R->acq++;
        burst_acq++;

        /* --- (frame-start hypothesis, detector lag) chosen by mismatch
         * count against the known prefix.  The periodic 0,3 preamble gives
         * an equivalent MF peak at every +-2k offset, so dn sweeps the whole
         * preamble range; only the non-periodic sync table resolves it. --- */
        long best_mis = (long)1 << 30, best_dnb = 0, best_lag = 0;
        float best_dsc = scale, best_cfo = 0.0f;
        /* CFO estimate: phase slope across the coherent sub-window MF sums
         * at best_n.  Template-aware, hence valid even for the periodic 0,3
         * preamble, whose self-autocorrelation at any pattern-period lag
         * folds the +-3*dev tone pair and mirrors/locks the phase (observed:
         * true +260 Hz estimated -262).  Ambiguity +-rate/(2*sub) = +-300 Hz
         * is resolved by the g hypotheses below; the prefix decisions need
         * the CFO removed first (a constant dphi eats the +-1-unit margin). */
        float cfo_c = 0.0f;
        {
            float complex cacc = 0.0f;
            for (unsigned w = 0u; w + 1u < nwsub; w++)
                cacc += conjf(best_rr[w]) * best_rr[w + 1u];
            if (cabsf(cacc) > 0.0f)
                cfo_c = -cargf(cacc) * (float)p.rate /
                        (2.0f * (float)M_PI * (float)sub);
        }
        const float cgmax = (float)p.rate / (8.0f * (float)p.k);
        long s0_ref = (long)best_n - (long)(NSFSK_TMPL_SYM0(p) * p.k);
        /* --- MF-prefilter the sample-granular dn sweep ---
         * every dn that can be the frame start has its template window
         * (s0 + SYM0*k) at best_n + dn, and a REAL start scores near the
         * span best even under CFO/ISI (the 4-symbol sub-windows stay
         * coherent for any in-tolerance offset).  Evaluating the prefix on
         * all span (~15k) offsets cost O(span*k) per ppm trial - seconds at
         * k=160, starving the live decoder (acoustic session 2026-10-28);
         * strided NMS peaks >= 40% of the span best, expanded by
         * +-MS/2+4 samples each, reproduce the full sweep's winner at
         * ~50-100x less work. */
        float scb = 0.0f;
        size_t ns = 0;
        for (long d = 0; dn_lo + d <= dn_hi; d += (long)MS, ns++) {
            long npos = (long)best_n + dn_lo + d;
            dn_sc[ns] = (npos >= 0 && (size_t)npos + tn <= N)
                      ? mf_score(y, tmpl, (size_t)npos, sub, nwsub, et, NULL)
                      : 0.0f;
            if (dn_sc[ns] > scb) scb = dn_sc[ns];
        }
        unsigned ncnd = 0;
        size_t npk_dbg = 0;
        if (scb > 0.0f) {
            const float thr = 0.4f * scb;
            long wq = (long)p.k / (2L * (long)MS); if (wq < 1) wq = 1;
            long pk[64]; size_t npk = 0;   /* npk mirrored to npk_dbg below */
            for (size_t js = 0; js < ns; js++) {
                if (dn_sc[js] < thr) continue;
                long q0 = (long)js > wq ? (long)js - wq : 0;
                long q1 = (long)js + wq + 1 < (long)ns
                        ? (long)js + wq + 1 : (long)ns;
                int ismax = 1;
                for (long q = q0; q < q1; q++)
                    if (dn_sc[q] > dn_sc[js]) { ismax = 0; break; }
                if (!ismax) continue;
                if (npk < 64u) pk[npk++] = (long)js;
                else {                       /* evict current worst */
                    unsigned mi = 0;
                    for (unsigned q = 1; q < 64u; q++)
                        if (dn_sc[pk[q]] < dn_sc[pk[mi]]) mi = q;
                    if (dn_sc[js] > dn_sc[pk[mi]]) pk[mi] = (long)js;
                }
            }
            for (size_t a = 1; a < npk; a++) {      /* best score first:
                tightens the mis-prune inside the sweep */
                long kj = pk[a]; float kv = dn_sc[kj]; size_t t = a;
                while (t > 0 && dn_sc[pk[t - 1u]] < kv) { pk[t] = pk[t-1u]; t--; }
                pk[t] = kj;
            }
            const long W = (long)MS / 2 + 4;       /* stride max is within
                MS/2 of the true peak; a few extra for prefix-vs-MF skew */
            for (size_t a = 0; a < npk; a++)
                for (long d = pk[a] * (long)MS - W; d <= pk[a] * (long)MS + W; d++) {
                    if (d < 0 || (unsigned long)d >= dn_span || ncnd >= 1024u)
                        continue;
                    int dup = 0;
                    for (unsigned q = 0; q < ncnd; q++)
                        if (cnd[q] == d) { dup = 1; break; }
                    if (!dup) cnd[ncnd++] = d;
                }
            /* scan in ASCENDING dn: the original full sweep resolved
             * equal-mis ties toward the SMALLEST dn (bank updates only on
             * strict <, non-bank on a higher calibrated scale).  Score-first
             * ordering flips those ties onto periodic-ridge impostors:
             * mis=0 against the periodic preamble, data-section garbage
             * (h1.0 battery pin, 2026-10-28). */
            for (unsigned a = 1; a < ncnd; a++) {
                long kj = cnd[a]; unsigned t = a;
                while (t > 0 && cnd[t - 1u] > cnd[t]) { cnd[t] = cnd[t-1u]; t--; }
                cnd[t] = kj;
            }
            npk_dbg = npk;
        }
        unsigned nswp = ncnd ? ncnd : dn_span;
        if (verbose > 1)
            fprintf(stderr, "rx: dn prefilter span=%zu stride=%u scb=%.3f "
                    "npk=%zu ncnd=%u best_n=%zu\n",
                    dn_span, MS, scb, npk_dbg, ncnd, best_n);
        for (unsigned si = 0; si < nswp; si++) {
            long dn = ncnd ? dn_lo + cnd[si] : dn_lo + (long)si;
            long s0c = s0_ref + dn;      /* sample-granular sweep: the MF
                                          * peak need not sit on a symbol
                                          * boundary under noise */
            if (s0c < 1) continue;       /* boundary_phi needs s0-1 >= 0;
                                          * frame start can be < k when the
                                          * burst begins at sample 0 */
            for (unsigned lg = 0; lg <= lag_max; lg++) {
                if (bank) {
                    /* magnitude decisions need no scale/CFO calibration:
                     * just count prefix mismatches.  Prune once mis can
                     * no longer beat best (wrong alignments give away
                     * themselves within a few symbols).  Scan from the
                     * sync-table END backwards: the preamble half of the
                     * prefix is 0,3-periodic and matches wrong dn/lag
                     * shifts for dozens of symbols, so a forward scan
                     * pays the full prefix per candidate (k=160 acoustic
                     * profile, 2026-10-28); the shift-distinguishing
                     * sync rejects losers within ~8 decisions. */
                    unsigned mis = 0, tot = 0;
                    for (unsigned wi = D0 - NSFSK_TMPL_SYM0(p); wi-- > 0u; ) {
                        unsigned w = NSFSK_TMPL_SYM0(p) + wi;
                        unsigned s = bank_decide(y, N,
                            (double)s0c + (double)(w + lg) * (double)p.k,
                            p.k, bank_e, NULL);
                        if (s == BANK_OOB) break;
                        tot++;
                        if (s != nsfsk_prefix_sym(&p, w)) mis++;
                        if ((long)mis > best_mis && tot >= 8u) break;
                    }
                    if (tot <= (D0 - NSFSK_TMPL_SYM0(p)) / 2u) continue;
                    if ((long)mis < best_mis) {
                        best_mis = mis; best_dnb = dn; best_lag = lg;
                        best_dsc = scale; best_cfo = 0.0f;
                    }
                    continue;
                }
                /* collect the known-prefix phase increments, then calibrate
                 * a per-candidate scale lambda against the expected symbol
                 * values a in {-3,-1,+1,+3} (LPF/ISI/bt all shift the
                 * effective phi amplitude) */
                float *phi_buf = cand_phi, *a_buf = cand_a;
                unsigned tot = 0;
                for (unsigned w = NSFSK_TMPL_SYM0(p); w < D0; w++) {
                    float phi = boundary_phi(y, N, (size_t)s0c, p.k, w + lg);
                    if (isnan(phi)) break;
                    phi_buf[tot] = phi;
                    a_buf[tot] = 2.0f * (float)nsfsk_prefix_sym(&p, w)
                               + 1.0f - (float)NSFSK_M;
                    tot++;
                }
                if (tot <= (D0 - NSFSK_TMPL_SYM0(p)) / 2u) continue;
                float den = 0.0f;
                for (unsigned q = 0; q < tot; q++)
                    den += a_buf[q] * a_buf[q];
                if (!(den > 0.0f)) continue;
                for (int g = 0; g < 3; g++) {
                    float cg = cfo_c + (g == 0 ? 0.0f :
                                (g == 1 ? 2.0f * cgmax : -2.0f * cgmax));
                    float dphi_g = 2.0f * (float)M_PI * cg
                                 * (float)p.k / (float)p.rate;
                    /* calibrate the decision scale on THIS hypothesis
                     * (num/den estimates the scale itself), then count
                     * prefix mismatches */
                    float num = 0.0f;
                    for (unsigned q = 0; q < tot; q++)
                        num += unwrap_ph(phi_buf[q] - dphi_g) * a_buf[q];
                    float sc = num / den;
                    if (!(sc > 0.1f) || sc > 3.0f) continue;
                    unsigned mis = 0;
                    for (unsigned q = 0; q < tot; q++)
                        if (detect_phi(unwrap_ph(phi_buf[q] - dphi_g), sc) !=
                            nsfsk_prefix_sym(&p, NSFSK_TMPL_SYM0(p) + q))
                            mis++;
                    if ((long)mis < best_mis ||
                        ((long)mis == best_mis && sc > best_dsc + 0.002f)) {
                        /* Tie-break on calibrated scale: at bt<0.8 the ISI-
                         * narrowed eye makes prefix mis count blind to +-1
                         * sample alignment (every class can tie at mis=1),
                         * but straddling a symbol boundary averages two
                         * tones and pulls the LS phase amplitude down
                         * (0.480 vs 0.501 at bt=0.5, tests/probe5.c), while
                         * the data loop inherits ~0.3% uncorrectable errors
                         * from the 1-sample-early lattice.  Max sc = eye
                         * centre.  A wrong g (CFO +-600 Hz) never reaches
                         * the tie pool: it blows up mis long before it can
                         * inflate sc (dphi_g error is +-pi/2 here). */
                        best_mis = mis; best_dnb = dn; best_lag = lg;
                        best_dsc = sc; best_cfo = cg;
                    }
                }
            }
        }
        if (verbose > 1)
            fprintf(stderr, "rx: prefix mism=%ld/%u doff=%ld s0=%ld lag=%ld "
                    "cfo=%+.0fHz\n",
                    best_mis, D0 - NSFSK_TMPL_SYM0(p), best_dnb,
                    s0_ref + best_dnb, best_lag, best_cfo);
        if (best_mis * 20 > (long)(D0 - NSFSK_TMPL_SYM0(p))) { /* >5%: rejects weak junk so a bad candidate cannot consume the frame */
            if (reject_limit && ++rejects >= reject_limit) {
                armed = 0; below = 0;            /* unaligned garbage (a
                                                  * mid-join burst costs a
                                                  * full sweep every 4k
                                                  * samples): idle to the
                                                  * gap and re-arm */
                i++;
            } else {
                i += 4u * (size_t)p.k;
            }
            continue;                            /* stay armed, rescan */
        }
        size_t s0 = (size_t)(s0_ref + best_dnb);
        float cfo = best_cfo;

        /* --- decode the fixed 1020 data symbols, RS + CRC gate ---
         * The data section is one 255-byte shortened RS codeword, so the
         * symbol count is known before decoding; only an 8-symbol sanity
         * peek (type+len: len <= 249-fec) runs before the full pass, and
         * wrong-alignment candidates die at the RS+CRC gate instead.
         * best_dsc is calibrated on the alternating 0,3 preamble, whose
         * heavy same-pulse ISI pulls it ~25% below the pattern-average
         * per-unit phase of random payload data.  The decision-directed
         * scale loop (fsc) converges to the data's own average unit so the
         * +-2*fsc thresholds stay in the true midpoints between the +-1
         * and +-3 phase clusters.  The first frame's early symbols sit in
         * the loop's convergence window, so a failed CRC re-decodes the
         * frame with fsc seeded from the previous pass' converged scale. */
        R->dbg_n = 0;
        float dphi = 2.0f * (float)M_PI * cfo * (float)p.k / (float)p.rate;
        float fsc  = best_dsc;
        float retry_fsc = best_dsc;
        float fsc_lo = 0.45f * best_dsc, fsc_hi = 1.6f * (float)M_PI * p.h
                     * (float)(p.k - 1u) / (float)p.k;
        unsigned need_syn = 4u * NSFSK_HDR_LEN;      /* type + len: 8 syms */
        unsigned got = 0, iter, frameok = 0, hdr_bad;
        unsigned len = 0, fix = 0;
        uint8_t ftype = 0;
        for (iter = 0; iter < (bank ? 1u : 3u); iter++) {
            if (iter) { dphi = 2.0f * (float)M_PI * cfo * (float)p.k
                                   / (float)p.rate;
                        fsc = retry_fsc; }
            double q = 0.0; size_t nq = 0;
            /* static symbol lattice (timing comes from the ppm trials,
             * see bank_decide) */
            double bpos = (double)s0
                        + (double)(D0 + best_lag) * (double)p.k;
            got = 0; hdr_bad = 0;
            for (; got < max_syn; got++) {
                unsigned s;
                if (bank) {
                    float qf = 0.25f;
                    s = bank_decide(y, N, bpos, p.k, bank_e, &qf);
                    if (s == BANK_OOB) break;
                    bpos += (double)p.k;
                    q += (double)qf;      /* ppm-trial quality metric */                    if (verbose && R->dbg_n < 4200u)
                        R->dbg_q[R->dbg_n++] = qf;
                } else {
                float phi = boundary_phi(y, N, s0, p.k, D0 + got + best_lag);
                if (isnan(phi)) break;
                float ph = unwrap_ph(phi - dphi);
                s = detect_phi(ph, fsc);
                /* decision-directed phase + scale tracking */
                float a_hat = 2.0f * (float)s + 1.0f - (float)NSFSK_M;
                float err   = wrap_pi(ph - a_hat * fsc);
                dphi += 0.01f * err;
                fsc  += 0.005f * err * a_hat / 3.0f;
                if (fsc < fsc_lo) fsc = fsc_lo;
                if (fsc > fsc_hi) fsc = fsc_hi;
                /* per-decision quality metric for the ppm trials */
                q += (double)(ph / fsc) * (double)a_hat;
                }
                nq++;
                wire[D0 + got] = (unsigned char)s;

                /* cheap early prune: the 2-byte header must carry a legal
                 * length (weaker than the old 24-sym header: junk
                 * alignments now fall to the RS+CRC gate below) */
                if (got + 1u == need_syn) {
                    nsfsk_symbols_to_data(wire + D0, need_syn, blk);
                    if (verbose > 1)
                        fprintf(stderr, "rx: hdr type=%u len=%u\n",
                                blk[0], blk[1]);
                    if (blk[1] > max_pl) { hdr_bad = 1; break; }
                }
            }
            R->qual += q; R->ndec += nq;
            if (hdr_bad || got < max_syn) break;

            nsfsk_symbols_to_data(wire + D0, max_syn, blk);
            if (nsfsk_frame_decode(&p, blk, &ftype, pay, &len, &fix)) {
                frameok = 1;
                break;
            }
            if (got >= 40u) retry_fsc = fsc;    /* seed the next pass */
        }
        i = s0 + (size_t)(D0 + got) * p.k + p.k;   /* resume after frame */
        armed = 0; below = 0;

        if (!frameok) {
            if (verbose) {
                fprintf(stderr, "rx: frame @%.1f ms RS/CRC REJECT",
                        1000.0 * (double)s0 / (double)p.rate);
                if (verbose > 1 && R->dbg_n > 0) {
                    float mn = 9; size_t mi = 0; double s = 0;
                    for (unsigned t = 0; t < R->dbg_n; t++) {
                        s += R->dbg_q[t];
                        if (R->dbg_q[t] < mn) { mn = R->dbg_q[t]; mi = t; }
                    }
                    fprintf(stderr, " qmean=%.3f qmin=%.3f@%zu/%u (%.0f%%) "
                                    "firstq=[%.2f %.2f %.2f %.2f]",
                            s / R->dbg_n, mn, mi, R->dbg_n,
                            100.0 * mi / R->dbg_n,
                            R->dbg_q[0], R->dbg_q[1], R->dbg_q[2], R->dbg_q[3]);
                }
                R->dbg_n = 0;
            }
            /* Only a full codeword that RS/CRC turned down is a LOST frame:
             * a candidate pruned on the header peek or cut short by the end
             * of the window never presented 2040 wire bits. */
            if (!hdr_bad && got >= max_syn) rx_result_note_lost(R, s0);
            R->bad++; continue;
        }

        if (rx_result_append(R, pay, len)) {
            fprintf(stderr, "rx: out of memory appending payload\n");
            R->bad++;
            continue;
        }
        if (R->nfr == R->fr_cap) {
            unsigned ncap = R->fr_cap ? R->fr_cap * 2u : 64u;
            rx_frame *nfr = realloc(R->fr, (size_t)ncap * sizeof *nfr);
            if (nfr != NULL) { R->fr = nfr; R->fr_cap = ncap; }
        }
        if (R->nfr < R->fr_cap) {
            R->fr[R->nfr].s0  = s0;
            R->fr[R->nfr].off = R->nbytes - len;  /* payload start in buf */
            R->fr[R->nfr].len = len;
            R->fr[R->nfr].fix = fix;
            R->nfr++;
        } else {
            fprintf(stderr, "rx: out of memory recording frame\n");
        }
        R->ok++;
        R->rs_bits += 8ull * fix;
        if (fix) R->ok_rs++;
        if (verbose)
            fprintf(stderr, "rx: frame @%.1f ms type=%u len=%u "
                            "cfo=%+.1f Hz lag=%u mism=%u rs=%u\n",
                    1000.0 * (double)s0 / (double)p.rate, ftype, len,
                    cfo, (int)best_lag, (int)best_mis, fix);
    }

    free(prof); free(wire); free(blk); free(pay);
    free(cand_phi); free(cand_a); free(dn_sc);
}
