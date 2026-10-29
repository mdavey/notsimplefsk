/*
 * channel — test impairment channel: s16le mono in -> impaired s16le out.
 *
 * Backward-compatible short options (same semantics as the original tool):
 *   -d  integer sample delay (zeros prepended)
 *   -D  fractional sample delay (linear interpolation, resampling walk)
 *   -r  sample-rate offset in ppm (output stretched by 1 + r*1e-6)
 *   -g  gain in dB
 *   -f  carrier frequency offset in Hz (analytic shift)
 *   -n  additive white Gaussian noise sigma, full-scale units (e.g. 0.01)
 *   -s  RNG seed (default 1)
 *
 * Additional long options model a more realistic audio path:
 *   --multipath MS:DB[,MS:DB...]  static taps (delay ms, gain dB, optional phase deg)
 *   --reverb RT60:WET_DB          direct path + synthetic diffuse tail
 *   --fade-doppler HZ             time-varying flat/frequency-selective fading
 *   --fade-rician K_DB            Rician K on the direct tap (default Rayleigh)
 *   --fade-spread MS              exponential power-delay profile spread
 *   --fade-taps N                 fading taps (default 4 when spread > 0)
 *   --snr-db DB                   AWGN from active-signal RMS (exclusive with -n)
 *   --noise-color white|pink|brown
 *   --bandpass LO:HI              RBJ biquad high-pass then low-pass
 *   --hum HZ:DBFS                 hum fundamental + 2nd/3rd harmonics
 *   --tone HZ:DBFS                in-band interferer (repeatable)
 *   --clicks RATE:AMP_DBFS:DECAY_MS
 *   --bursts RATE:LEN_MS:AMP_DBFS
 *   --dropouts RATE:LEN_MS:MODE   MODE=zero|drop
 *   --drop-at MS:LEN_MS[:MODE]    deterministic dropout (repeatable)
 *   --gain-steps MS:DB[,MS:DB...] piecewise gain changes
 *   --ppm-ramp START:END          linear sample-clock ramp
 *   --agc TARGET_DBFS:ATTACK_MS:RELEASE_MS[:MAX_GAIN_DB]
 *   --soft-clip DRIVE
 *   --quantize BITS               8..16
 *   --dc OFFSET                   full-scale DC offset
 *   --normalize                   output RMS -> input active RMS
 *   --report                      extra measured statistics on stderr
 *
 * Processing order (a model, not exact physics):
 *   clock/delay -> [analytic: multipath/fading/CFO] -> bandpass -> gain/steps
 *   -> hum/tones -> noise -> clicks/bursts -> AGC -> soft clip/DC/quantise
 *   -> dropouts -> s16 clip/write.
 *
 * All stochastic components are deterministic from -s.  stdout is always
 * pure s16le; diagnostics go to stderr.
 */
#include <complex.h>
#include <getopt.h>
#include <liquid/liquid.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_MULTIPATH  32
#define MAX_TONES       8
#define MAX_DROP_AT    32
#define MAX_GAIN_STEPS 32
#define MAX_REVERB_TAPS 64
#define MAX_STATIC_TAPS (MAX_MULTIPATH + MAX_REVERB_TAPS + 2)
#define MAX_FADE_TAPS   8
#define FADE_OSC       16
#define MAX_FFT_SIZE   (1u << 24)      /* 16 Mi samples: ~128 MB per FFT buffer */
#define MAX_SAMPLES    (1u << 26)      /* existing input cap */

/* ------------------------------------------------------------------ */
/* small deterministic RNG (xorshift64* + Box-Muller)                  */
/* ------------------------------------------------------------------ */
typedef struct {
    uint64_t s;
    double   spare;
    int      have_spare;
} rng_t;

static uint64_t splitmix64(uint64_t z)
{
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void rng_init(rng_t *r, uint64_t seed)
{
    uint64_t z = splitmix64(seed ^ 0x243f6a8885a308d3ULL);
    r->s = z ? z : 1u;
    r->spare = 0.0;
    r->have_spare = 0;
}

static uint64_t rng_u64(rng_t *r)
{
    uint64_t x = r->s;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    r->s = x;
    return x * 2685821657736338717ULL;
}

static double rng_double(rng_t *r)
{
    return (double)(rng_u64(r) >> 11) * (1.0 / 9007199254740992.0);
}

static double rng_gauss(rng_t *r)
{
    double u1, u2;
    if (r->have_spare) {
        r->have_spare = 0;
        return r->spare;
    }
    do { u1 = rng_double(r); } while (u1 <= 1e-15);
    u2 = rng_double(r);
    r->spare = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
    r->have_spare = 1;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ------------------------------------------------------------------ */
/* allocation helpers                                                  */
/* ------------------------------------------------------------------ */
static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1u);
    if (p == NULL) { fprintf(stderr, "channel: out of memory\n"); exit(1); }
    return p;
}

static void *xcalloc(size_t n, size_t sz)
{
    void *p = calloc(n ? n : 1u, sz ? sz : 1u);
    if (p == NULL) { fprintf(stderr, "channel: out of memory\n"); exit(1); }
    return p;
}

static void *xrealloc(void *p, size_t n)
{
    void *q = realloc(p, n ? n : 1u);
    if (q == NULL) { fprintf(stderr, "channel: out of memory\n"); exit(1); }
    return q;
}

/* ------------------------------------------------------------------ */
/* parsed impairment spec                                             */
/* ------------------------------------------------------------------ */
typedef struct {
    double delay_ms;
    double gain_db;
    double phase_rad;
} static_tap;

typedef struct {
    double hz;
    double db;
} tone_spec;

typedef struct {
    double at_ms;
    double len_ms;
    int    mode;                    /* 0 = zero, 1 = drop samples */
} drop_spec;

typedef struct {
    double at_ms;
    double db;
} gain_step;

typedef struct {
    double rate_hz;
    double attack_ms;
    double release_ms;
    double target_db;
    double max_gain_db;
} agc_spec;

/* fade state */
typedef struct {
    int      enabled;
    int      ntaps;
    int      nosc;
    double  *delay_ms;
    double  *weight;
    double  *power;
    double   doppler;
    double   rician_db;
    int      have_rician;
    double   los_amp;
    double   ray_amp;
    double complex *osc;            /* ntaps*nosc */
    double complex *rot;            /* ntaps*nosc */
    double complex *los;            /* ntaps */
    double   min_db, max_db;
    int      stats;
} fade_state;

/* ------------------------------------------------------------------ */
/* CLI parsing helpers                                                 */
/* ------------------------------------------------------------------ */
static int parse_double(const char *s, double *out)
{
    char *end = NULL;
    double v;
    if (s == NULL || *s == '\0') return -1;
    v = strtod(s, &end);
    if (end == s || *end != '\0' || !isfinite(v)) return -1;
    *out = v;
    return 0;
}

static int parse_long_strict(const char *s, long *out)
{
    char *end = NULL;
    long v;
    if (s == NULL || *s == '\0') return -1;
    v = strtol(s, &end, 0);
    if (end == s || *end != '\0') return -1;
    *out = v;
    return 0;
}

static int parse_mode(const char *s, int *mode)
{
    if (strcmp(s, "zero") == 0) { *mode = 0; return 0; }
    if (strcmp(s, "drop") == 0) { *mode = 1; return 0; }
    return -1;
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [options] < in.s16 > out.s16\n"
        "\n"
        "Backward-compatible short options:\n"
        "  -d delay        integer sample delay\n"
        "  -D fdelay       fractional sample delay\n"
        "  -r ppm          sample-rate offset in ppm\n"
        "  -g dB           gain\n"
        "  -f Hz           carrier frequency offset\n"
        "  -n sigma        AWGN sigma, full-scale units\n"
        "  -s seed         RNG seed (default 1)\n"
        "  -R rate         sample rate (default 48000)\n"
        "\n"
        "Long options:\n"
        "  --multipath MS:DB[,MS:DB...]      static taps (optional :PHASE_DEG)\n"
        "  --reverb RT60:WET_DB              diffuse room tail\n"
        "  --fade-doppler HZ                 enable fading (0 = static complex)\n"
        "  --fade-rician K_DB                Rician K on the direct tap\n"
        "  --fade-spread MS                  frequency-selective power spread\n"
        "  --fade-taps N                     fading taps (1..8)\n"
        "  --snr-db DB                       AWGN from active RMS (not with -n)\n"
        "  --noise-color white|pink|brown\n"
        "  --bandpass LO:HI                  RBJ high-pass then low-pass\n"
        "  --hum HZ:DBFS                     hum + 2nd/3rd harmonics\n"
        "  --tone HZ:DBFS                    in-band interferer (repeatable)\n"
        "  --clicks RATE:AMP_DBFS:DECAY_MS   Poisson clicks\n"
        "  --bursts RATE:LEN_MS:AMP_DBFS     Poisson noise bursts\n"
        "  --dropouts RATE:LEN_MS:MODE       Poisson dropouts (zero/drop)\n"
        "  --drop-at MS:LEN_MS[:MODE]        deterministic dropout (repeatable)\n"
        "  --gain-steps MS:DB[,MS:DB...]     piecewise gain changes\n"
        "  --ppm-ramp START:END              linear clock ramp\n"
        "  --agc TARGET_DBFS:ATTACK_MS:RELEASE_MS[:MAX_GAIN_DB]\n"
        "  --soft-clip DRIVE                 tanh saturation\n"
        "  --quantize BITS                   8..16 bit output\n"
        "  --dc OFFSET                       full-scale DC offset\n"
        "  --normalize                       scale output to input active RMS\n"
        "  --report                          measured stats on stderr\n"
        "  -h, --help\n",
        prog);
}

static void add_static_tap(static_tap *taps, int *n, double delay_ms,
                           double gain_db, double phase_deg)
{
    if (*n >= MAX_STATIC_TAPS) {
        fprintf(stderr, "channel: too many static taps (max %d)\n", MAX_STATIC_TAPS);
        exit(1);
    }
    taps[*n].delay_ms = delay_ms;
    taps[*n].gain_db = gain_db;
    taps[*n].phase_rad = phase_deg * M_PI / 180.0;
    (*n)++;
}

static int cmp_drop(const void *a, const void *b)
{
    double x = ((const drop_spec *)a)->at_ms;
    double y = ((const drop_spec *)b)->at_ms;
    return (x > y) - (x < y);
}

static int cmp_gain_step(const void *a, const void *b)
{
    double x = ((const gain_step *)a)->at_ms;
    double y = ((const gain_step *)b)->at_ms;
    return (x > y) - (x < y);
}

/* ------------------------------------------------------------------ */
/* resampling walk                                                     */
/* ------------------------------------------------------------------ */
static double interp_real(const int16_t *in, size_t n, double pos)
{
    size_t i0;
    double f;
    if (n == 0 || pos < 0.0 || pos > (double)(n - 1)) return 0.0;
    i0 = (size_t)pos;
    f = pos - (double)i0;
    if (i0 + 1u >= n) return (double)in[i0] / 32768.0;
    return ((1.0 - f) * (double)in[i0] + f * (double)in[i0 + 1u]) / 32768.0;
}

/* Constant ppm path preserves the original output-length rule. */
static double *resample_constant(const int16_t *in, size_t n, double ratio,
                                 long delay, double fdelay, size_t *out_n)
{
    double md = ((double)n + (double)delay) / ratio;
    size_t m = (n > 0 && md > 0.0) ? (size_t)md : 0u;
    double *x;

    if (m > MAX_SAMPLES) {
        fprintf(stderr, "channel: input too large after resampling\n");
        return NULL;
    }
    x = xmalloc(m * sizeof(double));
    for (size_t i = 0; i < m; i++) {
        double pos = (double)i * ratio - (double)delay - fdelay;
        x[i] = interp_real(in, n, pos);
    }
    *out_n = m;
    return x;
}

/* Linear ppm ramp; pos advances by the instantaneous ratio. */
static double *resample_ramp(const int16_t *in, size_t n, double ppm0,
                             double ppm1, long delay, double fdelay,
                             size_t *out_n)
{
    size_t cap = n + 1024u, m = 0;
    double *x = xmalloc(cap * sizeof(double));
    double pos = -(double)delay - fdelay;

    while (n > 0 && pos <= (double)(n - 1)) {
        double frac, ppm, ratio;
        if (m == cap) {
            if (cap > MAX_SAMPLES) {
                fprintf(stderr, "channel: input too large after resampling\n");
                free(x);
                return NULL;
            }
            cap *= 2u;
            x = xrealloc(x, cap * sizeof(double));
        }
        frac = pos <= 0.0 ? 0.0 : (pos >= (double)(n - 1) ? 1.0 : pos / (double)(n - 1));
        ppm = ppm0 + (ppm1 - ppm0) * frac;
        ratio = 1.0 + ppm * 1e-6;
        if (!(ratio > 0.0)) {
            fprintf(stderr, "channel: invalid ppm ramp ratio\n");
            free(x);
            return NULL;
        }
        x[m++] = interp_real(in, n, pos);
        pos += ratio;
    }
    *out_n = m;
    return x;
}

/* ------------------------------------------------------------------ */
/* FFT analytic signal                                                 */
/* ------------------------------------------------------------------ */
static size_t next_pow2(size_t x)
{
    size_t p = 1u;
    while (p < x) {
        if (p > (MAX_FFT_SIZE >> 1)) return 0;
        p <<= 1;
    }
    return p;
}

/* Returns a float-complex analytic signal of length nfft.  The real part is
 * forced back to the original samples so only the Hilbert imaginary part is
 * taken from the FFT.  liquid's backward FFT is unscaled (verified), so the
 * inverse is divided by nfft. */
static float complex *make_analytic(const double *x, size_t n, size_t nfft)
{
    float complex *a = (float complex *)fft_malloc(nfft * sizeof(float complex));
    float complex *b = (float complex *)fft_malloc(nfft * sizeof(float complex));
    fftplan fwd = NULL, bwd = NULL;
    int rc = 0;

    if (a == NULL || b == NULL) {
        fprintf(stderr, "channel: out of memory (FFT)\n");
        if (a) fft_free(a);
        if (b) fft_free(b);
        return NULL;
    }
    for (size_t i = 0; i < nfft; i++) {
        a[i] = (i < n) ? (float)x[i] : 0.0f;
        b[i] = 0.0f;
    }
    fwd = fft_create_plan((unsigned)nfft, a, b, LIQUID_FFT_FORWARD, 0);
    if (fwd == NULL) { fprintf(stderr, "channel: FFT plan failed\n"); rc = -1; goto done; }
    if (fft_execute(fwd) != 0) { fprintf(stderr, "channel: FFT failed\n"); rc = -1; goto done; }

    /* one-sided spectrum: double positive frequencies, zero negatives */
    for (size_t k = 1; k < nfft / 2u; k++)
        b[k] *= 2.0f;
    for (size_t k = 1; k < nfft / 2u; k++)
        b[nfft - k] = 0.0f;

    bwd = fft_create_plan((unsigned)nfft, b, a, LIQUID_FFT_BACKWARD, 0);
    if (bwd == NULL) { fprintf(stderr, "channel: inverse FFT plan failed\n"); rc = -1; goto done; }
    if (fft_execute(bwd) != 0) { fprintf(stderr, "channel: inverse FFT failed\n"); rc = -1; goto done; }

    for (size_t i = 0; i < nfft; i++) {
        float re = crealf(a[i]) / (float)nfft;
        float im = cimagf(a[i]) / (float)nfft;
        a[i] = (i < n ? (float)x[i] : 0.0f) + _Complex_I * im;
        (void)re;
    }

done:
    if (fwd) fft_destroy_plan(fwd);
    if (bwd) fft_destroy_plan(bwd);
    fft_free(b);
    if (rc != 0) { fft_free(a); return NULL; }
    return a;
}

static double complex analytic_at(const float complex *a, size_t n, double pos)
{
    size_t i0;
    double f;
    if (n == 0 || pos < 0.0 || pos > (double)(n - 1)) return 0.0;
    i0 = (size_t)pos;
    f = pos - (double)i0;
    if (i0 + 1u >= n) return (double complex)a[i0];
    return (1.0 - f) * (double complex)a[i0] + f * (double complex)a[i0 + 1u];
}

/* ------------------------------------------------------------------ */
/* biquads / AGC                                                       */
/* ------------------------------------------------------------------ */
typedef struct {
    double b0, b1, b2, a1, a2;
    double z1, z2;
} biquad;

static void biquad_reset(biquad *f)
{
    f->z1 = f->z2 = 0.0;
}

static void biquad_lowpass(biquad *f, double fs, double fc, double q)
{
    double w0 = 2.0 * M_PI * fc / fs;
    double cw = cos(w0), sw = sin(w0);
    double alpha = sw / (2.0 * q);
    double a0 = 1.0 + alpha;
    f->b0 = (1.0 - cw) * 0.5 / a0;
    f->b1 = (1.0 - cw) / a0;
    f->b2 = f->b0;
    f->a1 = -2.0 * cw / a0;
    f->a2 = (1.0 - alpha) / a0;
    biquad_reset(f);
}

static void biquad_highpass(biquad *f, double fs, double fc, double q)
{
    double w0 = 2.0 * M_PI * fc / fs;
    double cw = cos(w0), sw = sin(w0);
    double alpha = sw / (2.0 * q);
    double a0 = 1.0 + alpha;
    f->b0 = (1.0 + cw) * 0.5 / a0;
    f->b1 = -(1.0 + cw) / a0;
    f->b2 = f->b0;
    f->a1 = -2.0 * cw / a0;
    f->a2 = (1.0 - alpha) / a0;
    biquad_reset(f);
}

static double biquad_run(biquad *f, double x)
{
    double y = f->b0 * x + f->z1;
    f->z1 = f->b1 * x - f->a1 * y + f->z2;
    f->z2 = f->b2 * x - f->a2 * y;
    return y;
}

static double active_rms(const double *x, size_t n)
{
    double peak = 0.0, th, sum = 0.0;
    size_t cnt = 0;
    if (n == 0) return 0.0;
    for (size_t i = 0; i < n; i++) {
        double a = fabs(x[i]);
        if (a > peak) peak = a;
    }
    th = peak * 0.01;
    if (th < 1e-12) th = 1e-12;
    for (size_t i = 0; i < n; i++) {
        double v = x[i];
        if (fabs(v) >= th) { sum += v * v; cnt++; }
    }
    if (cnt == 0) {
        for (size_t i = 0; i < n; i++) sum += x[i] * x[i];
        cnt = n;
    }
    return cnt ? sqrt(sum / (double)cnt) : 0.0;
}

/* ------------------------------------------------------------------ */
/* fade state helpers                                                  */
/* ------------------------------------------------------------------ */
static void fade_init(fade_state *f, rng_t *rng, double rate,
                      double doppler, double rician_db, int have_rician,
                      double spread_ms, int ntaps_req)
{
    int ntaps = 1;
    double k_lin, power_sum = 0.0;

    memset(f, 0, sizeof(*f));
    f->doppler = doppler;
    f->rician_db = rician_db;
    f->have_rician = have_rician;
    f->nosc = FADE_OSC;
    f->min_db = 1e9;
    f->max_db = -1e9;

    if (spread_ms > 0.0) {
        ntaps = ntaps_req > 0 ? ntaps_req : 4;
        if (ntaps < 1) ntaps = 1;
        if (ntaps > MAX_FADE_TAPS) ntaps = MAX_FADE_TAPS;
    }
    f->ntaps = ntaps;
    f->enabled = 1;
    f->delay_ms = xcalloc((size_t)ntaps, sizeof(double));
    f->weight = xcalloc((size_t)ntaps, sizeof(double));
    f->power = xcalloc((size_t)ntaps, sizeof(double));
    f->osc = xcalloc((size_t)ntaps * (size_t)f->nosc, sizeof(double complex));
    f->rot = xcalloc((size_t)ntaps * (size_t)f->nosc, sizeof(double complex));
    f->los = xcalloc((size_t)ntaps, sizeof(double complex));

    for (int j = 0; j < ntaps; j++) {
        f->delay_ms[j] = (double)j * spread_ms;
        f->power[j] = spread_ms > 0.0 ? exp(-(double)j) : 1.0;
        power_sum += f->power[j];
    }
    if (power_sum <= 0.0) power_sum = 1.0;
    for (int j = 0; j < ntaps; j++) {
        f->power[j] /= power_sum;
        f->weight[j] = sqrt(f->power[j]);
    }

    if (have_rician) {
        k_lin = pow(10.0, rician_db / 10.0);
        f->los_amp = sqrt(k_lin / (k_lin + 1.0));
        f->ray_amp = sqrt(1.0 / (k_lin + 1.0));
    } else {
        f->los_amp = 0.0;
        f->ray_amp = 1.0;
    }

    for (int j = 0; j < ntaps; j++) {
        if (j == 0 && f->los_amp > 0.0) {
            double ph = 2.0 * M_PI * rng_double(rng);
            f->los[j] = f->los_amp * cexp(_Complex_I * ph);
        }
        for (int n = 0; n < f->nosc; n++) {
            double alpha = 2.0 * M_PI * rng_double(rng);
            double ph = 2.0 * M_PI * rng_double(rng);
            double omega = 2.0 * M_PI * doppler * cos(alpha) / rate;
            f->osc[j * f->nosc + n] = cexp(_Complex_I * ph);
            f->rot[j * f->nosc + n] = cexp(_Complex_I * omega);
        }
    }
}

static void fade_free(fade_state *f)
{
    free(f->delay_ms);
    free(f->weight);
    free(f->power);
    free(f->osc);
    free(f->rot);
    free(f->los);
}

static void fade_step(fade_state *f)
{
    for (int j = 0; j < f->ntaps; j++)
        for (int n = 0; n < f->nosc; n++)
            f->osc[j * f->nosc + n] *= f->rot[j * f->nosc + n];
}

/* complex fading gain for one tap at the current time */
static double complex fade_gain(const fade_state *f, int tap)
{
    double complex sum = 0.0;
    for (int n = 0; n < f->nosc; n++)
        sum += f->osc[tap * f->nosc + n];
    return f->los[tap] + f->ray_amp * sum / sqrt((double)f->nosc);
}

static void fade_renorm(fade_state *f)
{
    for (int j = 0; j < f->ntaps; j++) {
        for (int n = 0; n < f->nosc; n++) {
            double complex *z = &f->osc[j * f->nosc + n];
            double a = cabs(*z);
            if (a > 1e-12) *z /= a;
            else *z = 1.0;
        }
    }
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */
int main(int argc, char **argv)
{
    long   delay = 0;
    long   rate_l = 48000;
    double fdelay = 0.0, ppm = 0.0, gain_db = 0.0, cfo = 0.0, sigma = 0.0;
    unsigned long seed = 1;
    int    sigma_set = 0, ppm_set = 0;
    int    snr_set = 0;
    double snr_db = 0.0;
    const char *noise_color = "white";
    int    report = 0, normalize = 0;
    double band_lo = 0.0, band_hi = 0.0;
    int    bandpass = 0, band_hi_nyq = 0;
    double hum_hz = 0.0, hum_db = 0.0;
    int    hum_set = 0;
    tone_spec tones[MAX_TONES];
    int    ntones = 0;
    double click_rate = 0.0, click_amp_db = 0.0, click_decay_ms = 0.0;
    int    clicks_set = 0;
    double burst_rate = 0.0, burst_len_ms = 0.0, burst_amp_db = 0.0;
    int    bursts_set = 0;
    double dropout_rate = 0.0, dropout_len_ms = 0.0;
    int    dropout_mode = 0, dropouts_set = 0;
    drop_spec drop_at[MAX_DROP_AT];
    int    ndrop_at = 0;
    gain_step gain_steps[MAX_GAIN_STEPS];
    int    ngain_steps = 0;
    double ppm_ramp0 = 0.0, ppm_ramp1 = 0.0;
    int    ppm_ramp_set = 0;
    agc_spec agc;
    int    agc_set = 0;
    double soft_clip_drive = 0.0;
    int    quantize_bits = 16;
    double dc_offset = 0.0;
    static_tap staps[MAX_STATIC_TAPS];
    int    nstaps = 0;
    double reverb_rt60 = 0.0, reverb_wet_db = 0.0;
    int    reverb_set = 0;
    int    fade_doppler_set = 0, fade_taps_set = 0;
    double fade_doppler = 0.0, fade_rician_db = 0.0, fade_spread_ms = 0.0;
    int    fade_rician_set = 0, fade_taps_req = 0;
    int    opt;
    int    longidx = 0;

    enum {
        OPT_MULTIPATH = 1000, OPT_REVERB, OPT_FADE_DOPPLER, OPT_FADE_RICIAN,
        OPT_FADE_SPREAD, OPT_FADE_TAPS, OPT_SNR_DB, OPT_NOISE_COLOR,
        OPT_BANDPASS, OPT_HUM, OPT_TONE, OPT_CLICKS, OPT_BURSTS,
        OPT_DROPOUTS, OPT_DROP_AT, OPT_GAIN_STEPS, OPT_PPM_RAMP,
        OPT_AGC, OPT_SOFT_CLIP, OPT_QUANTIZE, OPT_DC, OPT_NORMALIZE,
        OPT_REPORT
    };
    static struct option longopts[] = {
        {"multipath",   required_argument, 0, OPT_MULTIPATH},
        {"reverb",      required_argument, 0, OPT_REVERB},
        {"fade-doppler",required_argument, 0, OPT_FADE_DOPPLER},
        {"fade-rician", required_argument, 0, OPT_FADE_RICIAN},
        {"fade-spread", required_argument, 0, OPT_FADE_SPREAD},
        {"fade-taps",   required_argument, 0, OPT_FADE_TAPS},
        {"snr-db",      required_argument, 0, OPT_SNR_DB},
        {"noise-color", required_argument, 0, OPT_NOISE_COLOR},
        {"bandpass",    required_argument, 0, OPT_BANDPASS},
        {"hum",         required_argument, 0, OPT_HUM},
        {"tone",        required_argument, 0, OPT_TONE},
        {"clicks",      required_argument, 0, OPT_CLICKS},
        {"bursts",      required_argument, 0, OPT_BURSTS},
        {"dropouts",    required_argument, 0, OPT_DROPOUTS},
        {"drop-at",     required_argument, 0, OPT_DROP_AT},
        {"gain-steps",  required_argument, 0, OPT_GAIN_STEPS},
        {"ppm-ramp",    required_argument, 0, OPT_PPM_RAMP},
        {"agc",         required_argument, 0, OPT_AGC},
        {"soft-clip",   required_argument, 0, OPT_SOFT_CLIP},
        {"quantize",    required_argument, 0, OPT_QUANTIZE},
        {"dc",          required_argument, 0, OPT_DC},
        {"normalize",   no_argument,       0, OPT_NORMALIZE},
        {"report",      no_argument,       0, OPT_REPORT},
        {"help",        no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    memset(&agc, 0, sizeof(agc));
    agc.max_gain_db = 20.0;

    while ((opt = getopt_long(argc, argv, "d:D:r:g:f:n:s:R:h",
                              longopts, &longidx)) != -1) {
        switch (opt) {
        case 'R': {
            long v;
            if (parse_long_strict(optarg, &v) || v <= 0 || v > 1000000) {
                fprintf(stderr, "channel: bad rate '%s'\n", optarg); return 1;
            }
            rate_l = v;
            break;
        }
        case 'd':
            if (parse_long_strict(optarg, &delay) || delay < 0) {
                fprintf(stderr, "channel: bad delay '%s'\n", optarg); return 1;
            }
            break;
        case 'D':
            if (parse_double(optarg, &fdelay)) {
                fprintf(stderr, "channel: bad fdelay '%s'\n", optarg); return 1;
            }
            break;
        case 'r':
            if (parse_double(optarg, &ppm)) {
                fprintf(stderr, "channel: bad ppm '%s'\n", optarg); return 1;
            }
            ppm_set = 1;
            break;
        case 'g':
            if (parse_double(optarg, &gain_db)) {
                fprintf(stderr, "channel: bad gain '%s'\n", optarg); return 1;
            }
            break;
        case 'f':
            if (parse_double(optarg, &cfo)) {
                fprintf(stderr, "channel: bad cfo '%s'\n", optarg); return 1;
            }
            break;
        case 'n':
            if (parse_double(optarg, &sigma) || sigma < 0.0) {
                fprintf(stderr, "channel: bad noise sigma '%s'\n", optarg); return 1;
            }
            sigma_set = 1;
            break;
        case 's':
            seed = (unsigned long)strtoul(optarg, NULL, 0);
            if (seed == 0) seed = 1;
            break;
        case 'h': usage(argv[0]); return 0;
        case OPT_MULTIPATH: {
            char *copy = strdup(optarg);
            char *save1 = NULL;
            if (copy == NULL) { fprintf(stderr, "channel: out of memory\n"); return 1; }
            for (char *tok = strtok_r(copy, ",", &save1); tok;
                 tok = strtok_r(NULL, ",", &save1)) {
                char *save2 = NULL;
                char *f1 = strtok_r(tok, ":", &save2);
                char *f2 = strtok_r(NULL, ":", &save2);
                char *f3 = strtok_r(NULL, ":", &save2);
                double dm, db, ph = 0.0;
                if (f1 == NULL || f2 == NULL || parse_double(f1, &dm) ||
                    parse_double(f2, &db) || dm < 0.0 ||
                    (f3 != NULL && parse_double(f3, &ph))) {
                    fprintf(stderr, "channel: bad multipath '%s'\n", optarg);
                    free(copy); return 1;
                }
                if (nstaps >= MAX_MULTIPATH) {
                    fprintf(stderr, "channel: too many multipath taps (max %d)\n", MAX_MULTIPATH);
                    free(copy); return 1;
                }
                add_static_tap(staps, &nstaps, dm, db, ph);
            }
            free(copy);
            break;
        }
        case OPT_REVERB: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            if (f1 == NULL || f2 == NULL ||
                parse_double(f1, &reverb_rt60) || parse_double(f2, &reverb_wet_db) ||
                reverb_rt60 <= 0.0) {
                fprintf(stderr, "channel: bad reverb '%s'\n", optarg);
                free(copy); return 1;
            }
            reverb_set = 1;
            free(copy);
            break;
        }
        case OPT_FADE_DOPPLER:
            if (parse_double(optarg, &fade_doppler) || fade_doppler < 0.0) {
                fprintf(stderr, "channel: bad fade doppler '%s'\n", optarg); return 1;
            }
            fade_doppler_set = 1;
            break;
        case OPT_FADE_RICIAN:
            if (parse_double(optarg, &fade_rician_db)) {
                fprintf(stderr, "channel: bad fade rician '%s'\n", optarg); return 1;
            }
            fade_rician_set = 1;
            break;
        case OPT_FADE_SPREAD:
            if (parse_double(optarg, &fade_spread_ms) || fade_spread_ms < 0.0) {
                fprintf(stderr, "channel: bad fade spread '%s'\n", optarg); return 1;
            }
            break;
        case OPT_FADE_TAPS: {
            long v;
            if (parse_long_strict(optarg, &v) || v < 1 || v > MAX_FADE_TAPS) {
                fprintf(stderr, "channel: bad fade taps '%s'\n", optarg); return 1;
            }
            fade_taps_req = (int)v;
            fade_taps_set = 1;
            break;
        }
        case OPT_SNR_DB:
            if (parse_double(optarg, &snr_db)) {
                fprintf(stderr, "channel: bad snr '%s'\n", optarg); return 1;
            }
            snr_set = 1;
            break;
        case OPT_NOISE_COLOR:
            if (strcmp(optarg, "white") && strcmp(optarg, "pink") &&
                strcmp(optarg, "brown")) {
                fprintf(stderr, "channel: bad noise color '%s'\n", optarg); return 1;
            }
            noise_color = optarg;
            break;
        case OPT_BANDPASS: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            if (f1 == NULL || f2 == NULL ||
                parse_double(f1, &band_lo)) {
                fprintf(stderr, "channel: bad bandpass '%s'\n", optarg);
                free(copy); return 1;
            }
            if (strcmp(f2, "nyq") == 0) band_hi_nyq = 1;
            else if (parse_double(f2, &band_hi)) {
                fprintf(stderr, "channel: bad bandpass '%s'\n", optarg);
                free(copy); return 1;
            }
            bandpass = 1;
            free(copy);
            break;
        }
        case OPT_HUM: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            if (f1 == NULL || f2 == NULL ||
                parse_double(f1, &hum_hz) || parse_double(f2, &hum_db) ||
                hum_hz <= 0.0) {
                fprintf(stderr, "channel: bad hum '%s'\n", optarg);
                free(copy); return 1;
            }
            hum_set = 1;
            free(copy);
            break;
        }
        case OPT_TONE: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            if (f1 == NULL || f2 == NULL || ntones >= MAX_TONES ||
                parse_double(f1, &tones[ntones].hz) ||
                parse_double(f2, &tones[ntones].db) ||
                tones[ntones].hz < 0.0) {
                fprintf(stderr, "channel: bad tone '%s'\n", optarg);
                free(copy); return 1;
            }
            ntones++;
            free(copy);
            break;
        }
        case OPT_CLICKS: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            char *f3 = f2 ? strtok_r(NULL, ":", &save) : NULL;
            if (f1 == NULL || f2 == NULL || f3 == NULL ||
                parse_double(f1, &click_rate) ||
                parse_double(f2, &click_amp_db) ||
                parse_double(f3, &click_decay_ms) ||
                click_rate < 0.0 || click_decay_ms <= 0.0) {
                fprintf(stderr, "channel: bad clicks '%s'\n", optarg);
                free(copy); return 1;
            }
            clicks_set = 1;
            free(copy);
            break;
        }
        case OPT_BURSTS: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            char *f3 = f2 ? strtok_r(NULL, ":", &save) : NULL;
            if (f1 == NULL || f2 == NULL || f3 == NULL ||
                parse_double(f1, &burst_rate) ||
                parse_double(f2, &burst_len_ms) ||
                parse_double(f3, &burst_amp_db) ||
                burst_rate < 0.0 || burst_len_ms <= 0.0) {
                fprintf(stderr, "channel: bad bursts '%s'\n", optarg);
                free(copy); return 1;
            }
            bursts_set = 1;
            free(copy);
            break;
        }
        case OPT_DROPOUTS: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            char *f3 = f2 ? strtok_r(NULL, ":", &save) : NULL;
            if (f1 == NULL || f2 == NULL || f3 == NULL ||
                parse_double(f1, &dropout_rate) ||
                parse_double(f2, &dropout_len_ms) ||
                parse_mode(f3, &dropout_mode) ||
                dropout_rate < 0.0 || dropout_len_ms <= 0.0) {
                fprintf(stderr, "channel: bad dropouts '%s'\n", optarg);
                free(copy); return 1;
            }
            dropouts_set = 1;
            free(copy);
            break;
        }
        case OPT_DROP_AT: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            char *f3 = f2 ? strtok_r(NULL, ":", &save) : NULL;
            drop_spec d;
            d.mode = 0;
            if (f1 == NULL || f2 == NULL || ndrop_at >= MAX_DROP_AT ||
                parse_double(f1, &d.at_ms) || parse_double(f2, &d.len_ms) ||
                d.at_ms < 0.0 || d.len_ms <= 0.0 ||
                (f3 != NULL && parse_mode(f3, &d.mode))) {
                fprintf(stderr, "channel: bad drop-at '%s'\n", optarg);
                free(copy); return 1;
            }
            drop_at[ndrop_at++] = d;
            free(copy);
            break;
        }
        case OPT_GAIN_STEPS: {
            char *copy = strdup(optarg);
            char *save1 = NULL;
            if (copy == NULL) { fprintf(stderr, "channel: out of memory\n"); return 1; }
            for (char *tok = strtok_r(copy, ",", &save1); tok;
                 tok = strtok_r(NULL, ",", &save1)) {
                char *save2 = NULL;
                char *f1 = strtok_r(tok, ":", &save2);
                char *f2 = strtok_r(NULL, ":", &save2);
                gain_step g;
                if (f1 == NULL || f2 == NULL || ngain_steps >= MAX_GAIN_STEPS ||
                    parse_double(f1, &g.at_ms) || parse_double(f2, &g.db) ||
                    g.at_ms < 0.0) {
                    fprintf(stderr, "channel: bad gain-steps '%s'\n", optarg);
                    free(copy); return 1;
                }
                gain_steps[ngain_steps++] = g;
            }
            free(copy);
            break;
        }
        case OPT_PPM_RAMP: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            if (f1 == NULL || f2 == NULL ||
                parse_double(f1, &ppm_ramp0) || parse_double(f2, &ppm_ramp1) ||
                1.0 + ppm_ramp0 * 1e-6 <= 0.0 ||
                1.0 + ppm_ramp1 * 1e-6 <= 0.0) {
                fprintf(stderr, "channel: bad ppm-ramp '%s'\n", optarg);
                free(copy); return 1;
            }
            ppm_ramp_set = 1;
            free(copy);
            break;
        }
        case OPT_AGC: {
            char *copy = strdup(optarg);
            char *save = NULL;
            char *f1 = copy ? strtok_r(copy, ":", &save) : NULL;
            char *f2 = f1 ? strtok_r(NULL, ":", &save) : NULL;
            char *f3 = f2 ? strtok_r(NULL, ":", &save) : NULL;
            char *f4 = f3 ? strtok_r(NULL, ":", &save) : NULL;
            if (f1 == NULL || f2 == NULL || f3 == NULL ||
                parse_double(f1, &agc.target_db) ||
                parse_double(f2, &agc.attack_ms) ||
                parse_double(f3, &agc.release_ms) ||
                agc.attack_ms <= 0.0 || agc.release_ms <= 0.0) {
                fprintf(stderr, "channel: bad agc '%s'\n", optarg);
                free(copy); return 1;
            }
            if (f4 != NULL &&
                (parse_double(f4, &agc.max_gain_db) || agc.max_gain_db < 0.0)) {
                fprintf(stderr, "channel: bad agc max gain '%s'\n", optarg);
                free(copy); return 1;
            }
            agc_set = 1;
            free(copy);
            break;
        }
        case OPT_SOFT_CLIP:
            if (parse_double(optarg, &soft_clip_drive) || soft_clip_drive <= 0.0) {
                fprintf(stderr, "channel: bad soft-clip '%s'\n", optarg); return 1;
            }
            break;
        case OPT_QUANTIZE: {
            long v;
            if (parse_long_strict(optarg, &v) || v < 8 || v > 16) {
                fprintf(stderr, "channel: bad quantize '%s'\n", optarg); return 1;
            }
            quantize_bits = (int)v;
            break;
        }
        case OPT_DC:
            if (parse_double(optarg, &dc_offset)) {
                fprintf(stderr, "channel: bad dc '%s'\n", optarg); return 1;
            }
            break;
        case OPT_NORMALIZE: normalize = 1; break;
        case OPT_REPORT: report = 1; break;
        default: usage(argv[0]); return 1;
        }
    }

    /* ---- validate combinations ---- */
    if (sigma_set && snr_set) {
        fprintf(stderr, "channel: -n and --snr-db are mutually exclusive\n");
        return 1;
    }
    if (ppm_set && ppm_ramp_set) {
        fprintf(stderr, "channel: -r and --ppm-ramp are mutually exclusive\n");
        return 1;
    }
    if ((fade_rician_set || fade_spread_ms != 0.0 || fade_taps_set) &&
        !fade_doppler_set) {
        fprintf(stderr, "channel: fading options require --fade-doppler\n");
        return 1;
    }
    if (fade_taps_set && fade_spread_ms <= 0.0) {
        fprintf(stderr, "channel: --fade-taps requires --fade-spread > 0\n");
        return 1;
    }
    if (fade_rician_set && (fade_rician_db < -100.0 || fade_rician_db > 100.0)) {
        fprintf(stderr, "channel: --fade-rician out of range (-100..100 dB)\n");
        return 1;
    }
    if (bandpass) {
        if (band_hi_nyq) band_hi = (double)rate_l * 0.5;
        if (band_lo < 0.0 || band_hi <= band_lo ||
            band_hi > (double)rate_l * 0.5 + 1e-9) {
            fprintf(stderr, "channel: invalid bandpass %.1f:%.1f at rate %ld\n",
                    band_lo, band_hi, rate_l);
            return 1;
        }
    }
    if (hum_set && (hum_hz <= 0.0 || hum_hz >= (double)rate_l * 0.5)) {
        fprintf(stderr, "channel: invalid hum frequency\n");
        return 1;
    }
    for (int i = 0; i < ntones; i++) {
        if (tones[i].hz >= (double)rate_l * 0.5) {
            fprintf(stderr, "channel: tone %.1f Hz above Nyquist\n", tones[i].hz);
            return 1;
        }
    }
    if (agc_set && (agc.target_db > 0.0 || agc.target_db < -160.0)) {
        fprintf(stderr, "channel: invalid AGC target\n");
        return 1;
    }
    if (quantize_bits < 8 || quantize_bits > 16) {
        fprintf(stderr, "channel: quantize must be 8..16\n");
        return 1;
    }
    if (dc_offset < -10.0 || dc_offset > 10.0) {
        fprintf(stderr, "channel: invalid dc offset\n");
        return 1;
    }
    if (reverb_set && reverb_rt60 < 0.001) {
        fprintf(stderr, "channel: reverb RT60 too small\n");
        return 1;
    }
    if (seed == 0) seed = 1;

    double rate = (double)rate_l;
    double ratio = 1.0 + ppm * 1e-6;
    if (!(ratio > 0.0)) {
        fprintf(stderr, "channel: invalid ppm\n");
        return 1;
    }

    /* ---- read all input ---- */
    size_t cap = 1u << 16, n = 0;
    int16_t *in = xmalloc(cap * sizeof(int16_t));
    for (;;) {
        if (n == cap) {
            if (cap > MAX_SAMPLES) {
                fprintf(stderr, "channel: input too large\n");
                free(in);
                return 1;
            }
            cap *= 2u;
            in = xrealloc(in, cap * sizeof(int16_t));
        }
        size_t got = fread(in + n, sizeof(int16_t), cap - n, stdin);
        n += got;
        if (!got) break;
    }

    /* ---- resampling walk ---- */
    size_t m = 0;
    double *x = NULL;
    if (n > 0) {
        x = ppm_ramp_set
            ? resample_ramp(in, n, ppm_ramp0, ppm_ramp1, delay, fdelay, &m)
            : resample_constant(in, n, ratio, delay, fdelay, &m);
        if (x == NULL) { free(in); return 1; }
    } else {
        x = xmalloc(sizeof(double));
    }

    /* ---- static taps + reverb ---- */
    if (reverb_set) {
        rng_t rev;
        double spacing, tail_end;
        int nt;
        rng_init(&rev, (uint64_t)seed ^ 0x6a09e667f3bcc909ULL);
        add_static_tap(staps, &nstaps, 0.0, 0.0, 0.0);
        tail_end = reverb_rt60 * 0.5;
        if (tail_end > 0.5) tail_end = 0.5;
        nt = (int)(tail_end / 0.005);
        if (nt < 1) nt = 1;
        if (nt > MAX_REVERB_TAPS) nt = MAX_REVERB_TAPS;
        spacing = tail_end / (double)nt;
        for (int j = 1; j <= nt; j++) {
            double t = (double)j * spacing;
            double db = reverb_wet_db - 60.0 * t / reverb_rt60;
            double ph = 2.0 * M_PI * rng_double(&rev);
            add_static_tap(staps, &nstaps, t * 1000.0, db, ph * 180.0 / M_PI);
        }
    }

    /* ---- fading init ---- */
    fade_state fade;
    memset(&fade, 0, sizeof(fade));
    if (fade_doppler_set) {
        rng_t frng;
        rng_init(&frng, (uint64_t)seed ^ 0x510e527fade682d1ULL);
        fade_init(&fade, &frng, rate, fade_doppler,
                  fade_rician_db, fade_rician_set, fade_spread_ms,
                  fade_taps_set ? fade_taps_req : 0);
    }

    /* ---- analytic channel ---- */
    int use_analytic = (cfo != 0.0) || (nstaps > 0) || fade.enabled;
    double max_delay_ms = 0.0;
    for (int i = 0; i < nstaps; i++)
        if (staps[i].delay_ms > max_delay_ms) max_delay_ms = staps[i].delay_ms;
    for (int i = 0; i < fade.ntaps; i++)
        if (fade.delay_ms[i] > max_delay_ms) max_delay_ms = fade.delay_ms[i];
    size_t max_delay = (size_t)ceil(max_delay_ms * rate / 1000.0);
    if (max_delay > MAX_SAMPLES) {
        fprintf(stderr, "channel: channel delay too large\n");
        free(x); free(in); fade_free(&fade); return 1;
    }
    size_t out_n = (m == 0) ? 0 : m + max_delay;
    double *y = xmalloc((out_n ? out_n : 1u) * sizeof(double));

    if (use_analytic && m > 0) {
        size_t nfft = next_pow2(out_n + 2u);
        float complex *ana;
        if (nfft == 0 || nfft > MAX_FFT_SIZE) {
            fprintf(stderr, "channel: input too long for analytic channel\n");
            free(y); free(x); free(in); fade_free(&fade); return 1;
        }
        ana = make_analytic(x, m, nfft);
        if (ana == NULL) {
            free(y); free(x); free(in); fade_free(&fade); return 1;
        }

        double complex cfo_z = 1.0;
        double complex cfo_rot = cexp(_Complex_I * 2.0 * M_PI * cfo / rate);
        int renorm_ctr = 0;

        for (size_t i = 0; i < out_n; i++) {
            double complex acc = 0.0;
            if (nstaps == 0 && !fade.enabled)
                acc = analytic_at(ana, m, (double)i);
            for (int t = 0; t < nstaps; t++) {
                double pos = (double)i - staps[t].delay_ms * rate / 1000.0;
                double amp = pow(10.0, staps[t].gain_db / 20.0);
                double complex g = amp * cexp(_Complex_I * staps[t].phase_rad);
                acc += g * analytic_at(ana, m, pos);
            }
            if (fade.enabled) {
                for (int t = 0; t < fade.ntaps; t++) {
                    double pos = (double)i - fade.delay_ms[t] * rate / 1000.0;
                    double complex g = fade_gain(&fade, t);
                    double a = cabs(g);
                    double db = 20.0 * log10(a > 1e-12 ? a : 1e-12);
                    if (db < fade.min_db) fade.min_db = db;
                    if (db > fade.max_db) fade.max_db = db;
                    acc += fade.weight[t] * g * analytic_at(ana, m, pos);
                }
                fade_step(&fade);
                if (++renorm_ctr >= 4096) { fade_renorm(&fade); renorm_ctr = 0; }
            }
            acc *= cfo_z;
            y[i] = creal(acc);
            cfo_z *= cfo_rot;
            if ((i & 1023u) == 0) {
                double a = cabs(cfo_z);
                if (a > 1e-12) cfo_z /= a;
            }
        }
        fft_free(ana);
    } else {
        for (size_t i = 0; i < m; i++) y[i] = x[i];
    }

    /* ---- bandpass ---- */
    if (bandpass) {
        /* fourth-order Butterworth as two cascaded biquads */
        static const double bq[2] = { 0.5411961001461971, 1.3065629648763765 };
        biquad hp[2], lp[2];
        if (band_lo > 0.0) {
            biquad_highpass(&hp[0], rate, band_lo, bq[0]);
            biquad_highpass(&hp[1], rate, band_lo, bq[1]);
            for (size_t i = 0; i < out_n; i++)
                y[i] = biquad_run(&hp[1], biquad_run(&hp[0], y[i]));
        }
        if (band_hi < rate * 0.5) {
            biquad_lowpass(&lp[0], rate, band_hi, bq[0]);
            biquad_lowpass(&lp[1], rate, band_hi, bq[1]);
            for (size_t i = 0; i < out_n; i++)
                y[i] = biquad_run(&lp[1], biquad_run(&lp[0], y[i]));
        }
    }

    /* ---- gain + gain steps ---- */
    {
        double g = pow(10.0, gain_db / 20.0);
        size_t si = 0;
        if (ngain_steps > 1) qsort(gain_steps, (size_t)ngain_steps, sizeof(gain_steps[0]), cmp_gain_step);
        for (size_t i = 0; i < out_n; i++) {
            while (si < (size_t)ngain_steps &&
                   gain_steps[si].at_ms * rate / 1000.0 <= (double)i) {
                g = pow(10.0, gain_steps[si].db / 20.0);
                si++;
            }
            y[i] *= g;
        }
    }

    /* ---- hum + tones ---- */
    if (hum_set) {
        double amp = pow(10.0, hum_db / 20.0);
        for (int h = 1; h <= 3; h++) {
            double a = h == 1 ? amp : (h == 2 ? amp * 0.5 : amp * 0.25);
            double th = 0.0, dth = 2.0 * M_PI * hum_hz * (double)h / rate;
            for (size_t i = 0; i < out_n; i++) {
                y[i] += a * sin(th);
                th += dth;
            }
        }
    }
    for (int t = 0; t < ntones; t++) {
        double amp = pow(10.0, tones[t].db / 20.0);
        double th = 0.0, dth = 2.0 * M_PI * tones[t].hz / rate;
        for (size_t i = 0; i < out_n; i++) {
            y[i] += amp * sin(th);
            th += dth;
        }
    }

    /* ---- noise ---- */
    double sig_rms = active_rms(y, out_n);
    double noise_sigma = sigma;
    if (snr_set) {
        if (sig_rms > 1e-12)
            noise_sigma = sig_rms / pow(10.0, snr_db / 20.0);
        else
            noise_sigma = 0.0;
    }
    if (noise_sigma > 0.0 && out_n > 0) {
        rng_t nrng;
        double *noise = xmalloc(out_n * sizeof(double));
        double mean = 0.0, rms = 0.0;
        rng_init(&nrng, (uint64_t)seed ^ 0x1f83d9abfb41bd6bULL);
        if (strcmp(noise_color, "white") == 0) {
            for (size_t i = 0; i < out_n; i++) noise[i] = rng_gauss(&nrng);
        } else if (strcmp(noise_color, "pink") == 0) {
            double b0 = 0, b1 = 0, b2 = 0, b3 = 0, b4 = 0, b5 = 0, b6 = 0;
            for (size_t i = 0; i < out_n; i++) {
                double w = rng_gauss(&nrng);
                b0 = 0.99886 * b0 + w * 0.0555179;
                b1 = 0.99332 * b1 + w * 0.0750759;
                b2 = 0.96900 * b2 + w * 0.1538520;
                b3 = 0.86650 * b3 + w * 0.3104856;
                b4 = 0.55000 * b4 + w * 0.5329522;
                b5 = -0.7616 * b5 - w * 0.0168980;
                noise[i] = b0 + b1 + b2 + b3 + b4 + b5 + b6 + w * 0.5362;
                b6 = w * 0.115926;
            }
        } else {
            double b = 0.0;
            for (size_t i = 0; i < out_n; i++) {
                b = 0.99 * b + 0.1 * rng_gauss(&nrng);
                noise[i] = b;
            }
        }
        for (size_t i = 0; i < out_n; i++) mean += noise[i];
        mean /= (double)out_n;
        for (size_t i = 0; i < out_n; i++) {
            noise[i] -= mean;
            rms += noise[i] * noise[i];
        }
        rms = sqrt(rms / (double)out_n);
        if (rms > 1e-12) {
            double scale = noise_sigma / rms;
            for (size_t i = 0; i < out_n; i++) y[i] += noise[i] * scale;
        }
        free(noise);
    }

    /* ---- clicks ---- */
    if (clicks_set && click_rate > 0.0 && out_n > 0) {
        rng_t crng;
        double peak = pow(10.0, click_amp_db / 20.0);
        double decay = click_decay_ms * rate / 1000.0;
        size_t tail = (size_t)(5.0 * decay) + 1u;
        double t = 0.0;
        rng_init(&crng, (uint64_t)seed ^ 0x5be0cd19137e2179ULL);
        for (;;) {
            double u = 1.0 - rng_double(&crng);
            t += -log(u) / click_rate;
            size_t idx = (size_t)(t * rate);
            if (idx >= out_n) break;
            double sign = rng_double(&crng) < 0.5 ? -1.0 : 1.0;
            size_t lim = idx + tail < out_n ? idx + tail : out_n;
            for (size_t k = idx; k < lim; k++)
                y[k] += sign * peak * exp(-(double)(k - idx) / decay);
        }
    }

    /* ---- noise bursts ---- */
    if (bursts_set && burst_rate > 0.0 && out_n > 0) {
        rng_t brng;
        double amp = pow(10.0, burst_amp_db / 20.0);
        size_t len = (size_t)(burst_len_ms * rate / 1000.0);
        if (len < 1) len = 1;
        double t = 0.0;
        rng_init(&brng, (uint64_t)seed ^ 0x3c6ef372fe94f82bULL);
        for (;;) {
            double u = 1.0 - rng_double(&brng);
            t += -log(u) / burst_rate;
            size_t idx = (size_t)(t * rate);
            if (idx >= out_n) break;
            size_t lim = idx + len < out_n ? idx + len : out_n;
            for (size_t k = idx; k < lim; k++)
                y[k] += amp * rng_gauss(&brng);
        }
    }

    /* ---- AGC ---- */
    if (agc_set && out_n > 0) {
        double target = pow(10.0, agc.target_db / 20.0);
        double att = exp(-1.0 / (agc.attack_ms * 0.001 * rate));
        double rel = exp(-1.0 / (agc.release_ms * 0.001 * rate));
        double maxg = pow(10.0, agc.max_gain_db / 20.0);
        double env = fabs(y[0]);
        for (size_t i = 0; i < out_n; i++) {
            double a = fabs(y[i]);
            double coeff = a > env ? att : rel;
            env = coeff * env + (1.0 - coeff) * a;
            double g = target / (env > 1e-12 ? env : 1e-12);
            if (g > maxg) g = maxg;
            y[i] *= g;
        }
    }

    /* ---- optional normalization ---- */
    if (normalize && out_n > 0) {
        double in_rms = active_rms(x, m);
        double out_rms = active_rms(y, out_n);
        if (in_rms > 1e-12 && out_rms > 1e-12) {
            double g = in_rms / out_rms;
            for (size_t i = 0; i < out_n; i++) y[i] *= g;
        }
    }

    /* ---- soft clip / DC / quantize ---- */
    if (soft_clip_drive > 0.0) {
        double denom = tanh(soft_clip_drive);
        for (size_t i = 0; i < out_n; i++)
            y[i] = tanh(soft_clip_drive * y[i]) / denom;
    }
    if (dc_offset != 0.0)
        for (size_t i = 0; i < out_n; i++) y[i] += dc_offset;
    if (quantize_bits < 16) {
        double levels = (double)(1u << (quantize_bits - 1));
        for (size_t i = 0; i < out_n; i++)
            y[i] = round(y[i] * levels) / levels;
    }

    /* ---- dropouts ---- */
    drop_spec *events = NULL;
    size_t nevents = 0, evcap = 0;
    for (int i = 0; i < ndrop_at; i++) {
        if (nevents == evcap) {
            evcap = evcap ? evcap * 2u : 16u;
            events = xrealloc(events, evcap * sizeof(*events));
        }
        events[nevents++] = drop_at[i];
    }
    if (dropouts_set && dropout_rate > 0.0 && out_n > 0) {
        rng_t drng;
        double t = 0.0;
        rng_init(&drng, (uint64_t)seed ^ 0x13198a2e03707344ULL);
        for (;;) {
            double u = 1.0 - rng_double(&drng);
            t += -log(u) / dropout_rate;
            size_t idx = (size_t)(t * rate);
            if (idx >= out_n) break;
            if (nevents == evcap) {
                evcap = evcap ? evcap * 2u : 16u;
                events = xrealloc(events, evcap * sizeof(*events));
            }
            events[nevents].at_ms = (double)idx * 1000.0 / rate;
            events[nevents].len_ms = dropout_len_ms;
            events[nevents].mode = dropout_mode;
            nevents++;
        }
    }
    unsigned ndrop_applied = 0;
    if (nevents > 0) {
        unsigned char *keep = NULL;
        int has_drop = 0;
        qsort(events, nevents, sizeof(*events), cmp_drop);
        for (size_t e = 0; e < nevents; e++) {
            size_t start = (size_t)(events[e].at_ms * rate / 1000.0);
            size_t len = (size_t)(events[e].len_ms * rate / 1000.0);
            if (start >= out_n) continue;
            if (len > out_n - start) len = out_n - start;
            ndrop_applied++;
            if (events[e].mode == 0) {
                memset(y + start, 0, len * sizeof(double));
            } else {
                if (keep == NULL) {
                    keep = xmalloc(out_n);
                    memset(keep, 1, out_n);
                }
                memset(keep + start, 0, len);
                has_drop = 1;
            }
        }
        if (has_drop) {
            size_t w = 0;
            size_t old_n = out_n;
            for (size_t i = 0; i < old_n; i++)
                if (keep[i]) y[w++] = y[i];
            out_n = w;
        }
        free(keep);
    }
    free(events);

    /* ---- count hard clips before the report/write pass ---- */
    size_t clips = 0;
    for (size_t i = 0; i < out_n; i++)
        if (y[i] > 1.0 || y[i] < -1.0) clips++;

    /* ---- report ---- */
    if (report) {
        double peak = 0.0;
        for (size_t i = 0; i < out_n; i++) {
            double a = fabs(y[i]);
            if (a > peak) peak = a;
        }
        double out_rms = active_rms(y, out_n);
        double snr_est = (noise_sigma > 0.0 && sig_rms > 0.0)
                       ? 20.0 * log10(sig_rms / noise_sigma) : INFINITY;
        fprintf(stderr,
            "channel-report: in_rms=%.6f out_rms=%.6f peak_dbfs=%.2f "
            "snr_db=%.2f sigma=%.6f fade_min_db=%.2f fade_max_db=%.2f "
            "dropouts=%u clips=%zu\n",
            active_rms(x, m), out_rms,
            peak > 1e-12 ? 20.0 * log10(peak) : -240.0, snr_est,
            noise_sigma, fade.enabled ? fade.min_db : 0.0,
            fade.enabled ? fade.max_db : 0.0, ndrop_applied, clips);
    }

    /* ---- write s16 ---- */
    int16_t *out = xmalloc((out_n ? out_n : 1u) * sizeof(int16_t));
    for (size_t i = 0; i < out_n; i++) {
        double v = y[i];
        if (v > 1.0) v = 1.0;
        else if (v < -1.0) v = -1.0;
        long q = lround(v * 32768.0);
        if (q > 32767) q = 32767;
        if (q < -32768) q = -32768;
        out[i] = (int16_t)q;
    }
    fwrite(out, sizeof(int16_t), out_n, stdout);

    fprintf(stderr,
        "channel: %zu -> %zu samples (delay=%ld fdelay=%+.2f ppm=%+.1f "
        "gain=%+.1fdB cfo=%+.1fHz sigma=%g clips=%zu drops=%u)\n",
        n, out_n, delay, fdelay, ppm, gain_db, cfo, noise_sigma, clips,
        ndrop_applied);

    free(out);
    free(y);
    if (x) free(x);
    free(in);
    fade_free(&fade);
    return 0;
}
