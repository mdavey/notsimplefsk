/*
 * fsk_tx — G4FSK frame transmitter: bytes on stdin, s16le mono audio on stdout
 *
 * Frames: stdin is chunked into payloads of at most 1024 bytes; each becomes
 * one frame [preamble][sync][Gray(NRZI(type|len16|payload|crc32))].
 * See common.h for the wire format and PHY notes.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <liquid/liquid.h>

#include "common.h"

static void usage(const char *prog)
{
    fprintf(stderr,
"G4FSK transmitter: raw bytes on stdin -> s16le mono audio on stdout\n"
"\n"
"usage: %s [options] < data > audio.s16\n"
"\n"
"  -t, --type N       frame type byte, 0..255 (default 0)\n"
"  -g, --gain dB      peak level in dBFS (default -12)\n"
"  -r, --rate HZ      sample rate (default 48000)\n"
"  -c, --fc HZ        carrier frequency (default 4000)\n"
"  -k, --k N          samples/symbol, even (default 20)\n"
"  -H, --h H          CPM modulation index (default 0.25)\n"
"  -b, --bt BT        Gaussian bandwidth-time product (default 0.8)\n"
"  -d, --delay N      CPM filter delay, symbols (default 1)\n"
"  -P, --preamble N   preamble length, symbols (default 72, max 65536)\n"
"      --pulse STR    gmsk | rcos | square (default gmsk)\n"
"  -h, --help         this text\n"
"\n"
"phy: tones at fc + {-3,-1,1,3}*h*Rs/2, Rs=rate/k; 2 bits/symbol\n"
"default: 2400 baud, 4800 bps, occupied ~6.6 kHz (Carson)\n"
"NOTE: TX/RX must use identical -r/-c/-k/-H/-b/-d/--pulse settings\n",
            prog);
}

int main(int argc, char **argv)
{
    nsfsk_params p;
    nsfsk_params_default(&p);
    long type = 0;
    float gain_db = -12.0f;

    static struct option longopts[] = {
        {"type",     required_argument, 0, 't'},
        {"gain",     required_argument, 0, 'g'},
        {"rate",     required_argument, 0, 'r'},
        {"fc",       required_argument, 0, 'c'},
        {"k",        required_argument, 0, 'k'},
        {"h-index",  required_argument, 0, 'H'},
        {"bt",       required_argument, 0, 'b'},
        {"delay",    required_argument, 0, 'd'},
        {"preamble", required_argument, 0, 'P'},
        {"pulse",    required_argument, 0,  1 },
        {"help",     no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long_only(argc, argv, "t:g:r:c:k:H:b:d:P:h",
                                   longopts, NULL)) != -1) {
        switch (opt) {
        case 't': {
            char *end = NULL;
            type = strtol(optarg, &end, 0);
            if (end == optarg || *end != '\0') {
                fprintf(stderr, "error: bad --type '%s'\n", optarg);
                return 1;
            }
            break;
        }
        case 'g': {
            char *end = NULL;
            gain_db = strtof(optarg, &end);
            if (end == optarg || *end != '\0' || !isfinite(gain_db)) {
                fprintf(stderr, "error: bad --gain '%s'\n", optarg);
                return 1;
            }
            break;
        }
        case 'r': p.rate = (unsigned)strtoul(optarg, NULL, 0); break;
        case 'c': p.fc = strtof(optarg, NULL); break;
        case 'k': p.k = (unsigned)strtoul(optarg, NULL, 0); break;
        case 'H': p.h = strtof(optarg, NULL); break;
        case 'b': p.bt = strtof(optarg, NULL); break;
        case 'd': p.delay = (unsigned)strtoul(optarg, NULL, 0); break;
        case 'P': p.preamble = (unsigned)strtoul(optarg, NULL, 0); break;
        case 1:   if (nsfsk_parse_pulse(optarg, &p.pulse)) {
                      fprintf(stderr, "error: bad --pulse '%s'\n", optarg);
                      return 1;
                  }
                  break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    if (type < 0 || type > 255) {
        fprintf(stderr, "error: type must be 0..255\n");
        return 1;
    }
    const char *err = nsfsk_params_validate(&p);
    if (err) { fprintf(stderr, "error: %s\n", err); return 1; }

    /* read all of stdin */
    size_t cap = 65536, n = 0;
    unsigned char *buf = malloc(cap);
    if (buf == NULL) { fprintf(stderr, "error: out of memory\n"); return 1; }
    for (;;) {
        if (n == cap) {
            if (cap > 64u * 1024u * 1024u) {
                fprintf(stderr, "error: input too large\n");
                return 1;
            }
            unsigned char *nb = realloc(buf, cap * 2u);
            if (nb == NULL) {
                fprintf(stderr, "error: out of memory\n");
                free(buf);
                return 1;
            }
            buf = nb;
            cap *= 2;
        }
        size_t got = fread(buf + n, 1, cap - n, stdin);
        n += got;
        if (got == 0) break;
    }

    cpfskmod mod = cpfskmod_create(NSFSK_BPS, p.h, p.k, p.delay, p.bt, p.pulse);
    if (mod == NULL) { fprintf(stderr, "error: cpfskmod_create failed\n"); return 1; }
    nco_crcf up = nco_crcf_create(LIQUID_NCO);
    nco_crcf_set_frequency(up, 2.0f * (float)M_PI * p.fc / (float)p.rate);

    float amp = powf(10.0f, gain_db / 20.0f);
    /* preamble is bounded by nsfsk_params_validate(), so this fits unsigned */
    unsigned cap_syms = (unsigned)NSFSK_MAX_FRAME_SYMS(p);
    unsigned char *syms = malloc(cap_syms);
    float complex *yk = malloc(p.k * sizeof(float complex));
    int16_t *out = malloc(8192 * sizeof(int16_t));
    if (syms == NULL || yk == NULL || out == NULL) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }
    size_t outn = 0;
    unsigned gap = p.rate / 25u;            /* 40 ms silence between frames */

    size_t off = 0;
    unsigned nframes = 0;
    double total_samples = 0.0;
    while (off < n || nframes == 0) {
        unsigned plen = (unsigned)(n - off);
        if (plen > NSFSK_MAX_PAYLOAD) plen = NSFSK_MAX_PAYLOAD;

        unsigned nsym = nsfsk_build_frame_syms(&p, (uint8_t)type,
                                               buf + off, plen, syms, cap_syms);
        if (nsym == 0) { fprintf(stderr, "error: build_frame_syms failed\n"); return 1; }

        if (nframes > 0) {                   /* inter-frame gap */
            for (unsigned i = 0; i < gap; i++) {
                out[outn++] = 0;
                if (outn == 8192) { fwrite(out, sizeof(int16_t), outn, stdout); outn = 0; }
            }
            total_samples += gap;
        }

        cpfskmod_reset(mod);
        /* trailing pad: covers the detector lag AND the RX zero-phase LPF
         * edge reach so the CRC trailer is never filtered against the
         * buffer end (values irrelevant; RX stops at the CRC) */
        unsigned pad = nsfsk_tx_pad(&p);
        for (unsigned s = 0; s < nsym + pad; s++) {
            cpfskmod_modulate(mod, s < nsym ? syms[s] : (s & 1u) ? 3u : 0u, yk);
            for (unsigned j = 0; j < p.k; j++) {
                float complex z;
                nco_crcf_mix_up(up, yk[j], &z);
                nco_crcf_step(up);
                float v = crealf(z) * amp;
                int32_t q = (int32_t)lrintf(v * 32767.0f);
                if (q > 32767) q = 32767;
                if (q < -32768) q = -32768;
                out[outn++] = (int16_t)q;
                if (outn == 8192) { fwrite(out, sizeof(int16_t), outn, stdout); outn = 0; }
            }
            total_samples += p.k;
        }

        off += plen;
        nframes++;
    }
    if (outn) fwrite(out, sizeof(int16_t), outn, stdout);
    fflush(stdout);

    fprintf(stderr, "tx: %u frame(s), %zu byte(s), %.0f samples (%.3f s) "
                    "rate=%u fc=%.0f k=%u h=%.3f bt=%.2f delay=%u pulse=%s\n",
            nframes, n, total_samples, total_samples / p.rate,
            p.rate, p.fc, p.k, p.h, p.bt, p.delay,
            p.pulse == LIQUID_CPFSK_GMSK ? "gmsk" :
            p.pulse == LIQUID_CPFSK_RCOS_FULL ? "rcos" : "square");

    free(buf); free(syms); free(yk); free(out);
    cpfskmod_destroy(mod);
    nco_crcf_destroy(up);
    return 0;
}
