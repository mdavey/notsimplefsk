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

/* complex linear-interpolation resampler: out[i] = in[i * ratio] */
float complex *resample_lin(const float complex *x, size_t n,
                            double ratio, size_t *m)
{
    size_t mo = (size_t)(((double)n - 1.0) / ratio) + 1u;
    float complex *o = malloc((mo + 4) * sizeof(float complex));
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

void rx_result_append(rx_result *R, const unsigned char *b, size_t n)
{
    if (R->nbytes + n > R->cap) {
        size_t cap = R->cap ? R->cap : 65536u;
        while (cap < R->nbytes + n) cap *= 2u;
        R->buf = realloc(R->buf, cap);
        R->cap = cap;
    }
    memcpy(R->buf + R->nbytes, b, n);
    R->nbytes += n;
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
void rx_run(const float complex *y, size_t N, unsigned nf,
            const float complex *tmpl, float et,
            const nsfsk_params p, int verbose, unsigned reject_limit,
            rx_result *R)
{
    /* ---- energy profile (moving average over win samples) ---- */
    unsigned win = (unsigned)fminf((float)(p.rate / 50),
                                   0.5f * (float)(p.preamble * p.k));
    if (win < p.k) win = p.k;
    float *prof = malloc(N * sizeof(float));
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
    const unsigned max_msg = NSFSK_HDR_LEN + NSFSK_MAX_PAYLOAD + NSFSK_CRC_LEN;
    const size_t   max_syn = (size_t)4u * max_msg;
    unsigned char *wire    = malloc(D0 + max_syn);         /* full frame syms */
    uint8_t  *msg  = malloc(max_msg);
    const float scale = (float)M_PI * p.h * (float)(p.k - 1u) / (float)p.k;
    const unsigned lag_max = 2u * p.delay + 3u;
    const unsigned silence_arm = p.rate / 50u;  /* 20 ms; TX gap is 40 ms */
    const float ethresh = 0.15f * pmax;

    float *cand_phi = malloc((size_t)D0 * sizeof(float));
    float *cand_a   = malloc((size_t)D0 * sizeof(float));

    int armed = 1;
    unsigned below = silence_arm;                /* start as if silent */
    unsigned rejects = 0;
    size_t i = win < N / 2u ? win : N / 2u;

    while (i + (size_t)(p.preamble * p.k) < N) {
        if (!armed) {
            if (prof[i] < ethresh) {
                if (++below >= silence_arm) { armed = 1; below = 0; }
            } else below = 0;
            i++;
            continue;
        }
        if (prof[i] < ethresh) { below = 0; i++; continue; }
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
        size_t best_n = 0;
        float  best_sc = 0.0f;
        float complex best_rr[8], rr[8];
        for (size_t n = n_lo; n + tn <= n_hi; n++) {
            float accm = 0.0f, e = 0.0f;
            for (unsigned w = 0; w < nwsub; w++) {
                float complex r = 0.0f;
                for (unsigned t = w * sub; t < (w + 1u) * sub; t++) {
                    r += tmpl[t] * conjf(y[n + t]);
                    e += cabsf(y[n + t]) * cabsf(y[n + t]);
                }
                rr[w] = r;
                accm += cabsf(r);
            }
            float sc = (e > 0.0f) ? (accm * accm) / (e * et) : 0.0f;
            if (sc > best_sc) { best_sc = sc; best_n = n;
                                memcpy(best_rr, rr, sizeof rr); }
        }
        if (verbose > 1)
            fprintf(stderr, "rx: mf best_n=%zu score=%.3f cfo_stage0\n",
                    best_n, best_sc);
        if (best_sc < 0.25f) {                      /* no real frame here */
            i += (size_t)(p.preamble * p.k);
            armed = 0; below = 0;
            continue;
        }

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
        for (long dn = -(long)(p.preamble * p.k + win + 2u * nf);
             dn <= (long)(4 * p.k); dn++) {
            long s0c = s0_ref + dn;      /* sample-granular sweep: the MF
                                          * peak need not sit on a symbol
                                          * boundary under noise */
            if (s0c < 1) continue;       /* boundary_phi needs s0-1 >= 0;
                                          * frame start can be < k when the
                                          * burst begins at sample 0 */
            for (unsigned lg = 0; lg <= lag_max; lg++) {
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
                    if ((long)mis < best_mis) {
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

        /* --- decode data-section symbols, header first, CRC gate ---
         * best_dsc is calibrated on the alternating 0,3 preamble, whose
         * heavy same-pulse ISI pulls it ~25% below the pattern-average
         * per-unit phase of random payload data.  The decision-directed
         * scale loop (fsc) converges to the data's own average unit so the
         * +-2*fsc thresholds stay in the true midpoints between the +-1
         * and +-3 phase clusters.  The first frame's early symbols sit in
         * the loop's convergence window, so a failed CRC re-decodes the
         * frame with fsc seeded from the previous pass' converged scale. */
        float dphi = 2.0f * (float)M_PI * cfo * (float)p.k / (float)p.rate;
        float fsc  = best_dsc;
        float retry_fsc = best_dsc;
        float fsc_lo = 0.45f * best_dsc, fsc_hi = 1.6f * (float)M_PI * p.h
                     * (float)(p.k - 1u) / (float)p.k;
        unsigned need_syn = 4u * NSFSK_HDR_LEN;      /* type + len16 */
        unsigned got = 0, need_total = 0, iter, crcok = 0;
        unsigned len = 0;
        uint8_t crc[4];
        for (iter = 0; iter < 3u; iter++) {
            if (iter) { dphi = 2.0f * (float)M_PI * cfo * (float)p.k
                                   / (float)p.rate;
                        fsc = retry_fsc; }
            double q = 0.0; size_t nq = 0;
            got = 0; need_total = 0;
            for (; got < max_syn; got++) {
                float phi = boundary_phi(y, N, s0, p.k, D0 + got + best_lag);
                if (isnan(phi)) break;
                float ph = unwrap_ph(phi - dphi);
                unsigned s = detect_phi(ph, fsc);
                /* decision-directed phase + scale tracking */
                float a_hat = 2.0f * (float)s + 1.0f - (float)NSFSK_M;
                float err   = wrap_pi(ph - a_hat * fsc);
                dphi += 0.01f * err;
                fsc  += 0.005f * err * a_hat / 3.0f;
                if (fsc < fsc_lo) fsc = fsc_lo;
                if (fsc > fsc_hi) fsc = fsc_hi;
                /* per-decision quality metric for the ppm trials */
                q += (double)(ph / fsc) * (double)a_hat;
                nq++;
                wire[D0 + got] = (unsigned char)s;

                if (need_total == 0 && got + 1u >= need_syn) {
                    nsfsk_symbols_to_data(wire + D0, need_syn, msg);
                    unsigned hl = ((unsigned)msg[1] << 8) | msg[2];
                    if (verbose > 1)
                        fprintf(stderr, "rx: hdr type=%u len=%u\n",
                                msg[0], hl);
                    if (hl > NSFSK_MAX_PAYLOAD) break;
                    need_total = 4u * (NSFSK_HDR_LEN + hl + NSFSK_CRC_LEN);
                }
                if (need_total && got + 1u >= need_total) break;
            }
            R->qual += q; R->ndec += nq;
            if (need_total == 0 || got + 1u < need_total) break;

            nsfsk_symbols_to_data(wire + D0, need_total, msg);
            len = ((unsigned)msg[1] << 8) | msg[2];
            nsfsk_crc32(msg, NSFSK_HDR_LEN + len, crc);
            if (memcmp(crc, msg + NSFSK_HDR_LEN + len, 4) == 0) {
                crcok = 1;
                break;
            }
            if (got >= 40u) retry_fsc = fsc;    /* seed the next pass */
        }
        i = s0 + (size_t)(D0 + got) * p.k + p.k;   /* resume after frame */
        armed = 0; below = 0;

        if (!crcok) { R->bad++; continue; }

        rx_result_append(R, msg + NSFSK_HDR_LEN, len);
        R->ok++;
        if (verbose)
            fprintf(stderr, "rx: frame @%.1f ms type=%u len=%u "
                            "cfo=%+.1f Hz lag=%u mism=%u\n",
                    1000.0 * (double)s0 / (double)p.rate, msg[0], len,
                    cfo, (int)best_lag, (int)best_mis);
    }

    free(prof); free(wire); free(msg);
    free(cand_phi); free(cand_a);
}
