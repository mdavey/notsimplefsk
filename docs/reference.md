# notsimplefsk — G4FSK modem reference

Two programs move arbitrary bytes over audio using 4-level continuous-phase
FSK with Gaussian pulse shaping (G4FSK), built on liquid-dsp:

    bin/fsk_tx   bytes (stdin)  -> s16le mono audio (stdout)
    bin/fsk_rx   s16le audio    -> recovered bytes (stdout)

All diagnostics go to stderr; stdout is always pure payload/audio so the
programs compose in pipelines. RX PHY options must match the TX exactly.

## Build

    make            # bin/{fsk_tx,fsk_rx,fsk_rx_stream,modem_smoke,channel,rs_selftest}
    make test       # fast gate: rs_selftest + roundtrip.sh + stream.sh
                    #            + channel_selftest.py
                    #            + rx_speed.sh (pathological + live-noise RX)
    make stress     # realistic impairment/PHY matrix (quick, opt-in)
    make stress-full # full matrix, 4 seeds (long; use --update-baseline)
    make speed-full # opt-in 2 KB pathological RX timing case

Requires libliquid >= 1.6 (`-lliquid`) and python3 for the stress harness.

## Quick start

    ./bin/fsk_tx < secret.bin > audio.s16
    ./bin/fsk_rx < audio.s16 > got.bin
    cmp secret.bin got.bin

Audition/playback: `sox audio.s16 -r 48000 -c 1 audio.wav`.
The channel simulator adds impairments for testing:

    ./bin/fsk_tx < in.bin | ./bin/channel -f 200 -r 300 -n 0.01 | ./bin/fsk_rx > out.bin

## Wire format (per frame)

    [ preamble 0,3,0,3,... x P ][ sync 24-sym table ][ Gray(NRZI(PN( codeword ))) ][ pad ]

The data section is always exactly ONE shortened Reed–Solomon codeword:
RS(255, 255−P) over GF(256), prim poly 0x11d, generator roots α..α^P.
The 255-byte block is independent of tier and payload, so its wire image is
always 1020 symbols (4 syms/byte at 2 bits/sym):

    codeword (255 B) = [ type 1 | len 1 | payload 0..249−P | crc32_be 4 | zeros to k ][ parity P ]
    data region k = 255 − P; RS parity P bytes; t = P/2 corrected bytes/block

    type    1 byte, arbitrary
    len     payload length, single byte (max payload = 249 − P)
    crc32   liquid LIQUID_CRC_32 over type|len|payload — INSIDE the codeword.
            RS correction + CRC is the frame gate: an uncorrectable block or a
            CRC mismatch rejects the whole frame (all-or-nothing, R->bad++).
            Padding lives inside the RS data region and is excluded from the
            CRC, so RX never discovers lengths or padding.
    parity  systematic RS parity over the full k-byte data region
    PN      additive whitening, m-sequence LFSR x^23+x^18+1, all-ones seed,
            re-armed per frame at the data section (period 2^23−1 bits ≫
            2040); keeps payload statistics off the air -- un-whitened, runs
            of zeros modulate as one pure tone; descramble is the same local
            PN (no error multiplication)
    NRZI    enc[n] = enc[n-1] XOR d[n] (d = whitened bit), state 0 at first
            bit (kills the 180-deg ambiguity of non-coherent detection)
    Gray    dibit -> symbol index v ^ (v>>1) = {0,1,3,2}: one wrong symbol
            corrupts at most ONE codeword byte, matching the RS byte-error
            model exactly
    pad     nsfsk_tx_pad() dummy symbols: covers detector lag and the RX
            zero-phase LPF edge reach
    sync    non-periodic 24-symbol table; the periodic 0,3 preamble has an
            MF sidelobe every +-2 symbols, only the sync table resolves it

### FEC tiers (`--fec P`, TX and RX must match)

    P    t corrected B/frame   max payload B    RS rate
    16           8                 233          239/255
    32 (dfl)    16                 217          223/255  (~7/8)
    64          32                 185          191/255
    84          42                 165          171/255

The decoder (src/rs.c, ~200 lines: syndromes → Berlekamp–Massey → Chien →
Forney) re-verifies the syndromes after correcting and returns −1 unless the
result is syndrome-clean — it never emits a bogus correction; the CRC then
gates the frame anyway.  Liquid-dsp's RS API is NOT used: the `.so` shipped
on the dev box contains only "libfec not installed" stubs (verified via
nm/strings).

Tier nesting (pinned fecmis-lo/-hi): every tier's generator roots are the
consecutive set a..a^P, so a higher-tier codeword is a VALID codeword of
every lower tier: an RX with --fec BELOW the TX tier decodes byte-exact and
loses only correction budget.  An RX tier ABOVE the TX tier sees a
non-codeword and rejects (RS+CRC gate).  RX --fec must therefore be ≤ TX.

Why one fixed block per frame: GF(256) caps any codeword at 255 B, so a
single-block frame gives a data section whose length is known BEFORE
decoding (no length discovery, no per-length decode paths); frame airtime
becomes tier-independent; and the decode loop is one straight path.  The
cost: tiny messages now occupy a full frame too (~3.7 s at 300 Bd where the
old 124-sym minimum cost ~0.4 s), and bulk transfers lose ~25 % airtime to
per-frame preamble/sync versus the old 1024-byte frames.  Decode-side benefit: ≤ t corrupt bytes per frame (short
acoustic dropouts, the bt ≤ 0.5 ISI pattern floor) turn from a lost frame
into a clean decode — pinned by the fecwin cases in roundtrip.sh, which the
pre-FEC binary fails on identical channel output.

Airtime per frame (1116 syms incl. preamble+sync; tier-independent):

    300 Bd  (k=160 @ 48k)   ~3.7 s   any message <= max payload
    1200 Bd (k=40  @ 48k)   ~0.93 s
    2400 Bd (k=20  @ 48k)   ~0.47 s

Inputs larger than the tier's max payload are chunked exactly as before:
TX emits one frame per ≤ (249−P)-byte chunk of stdin separated by 40 ms of
silence (the RX energy gate re-arms on it); larger payloads become multiple
frames and the consumer concatenates the delivered payloads in order.

## PHY parameters (identical flags on both programs)

    -r rate   sample rate Hz      (default 48000, min 8000)
    -c fc     carrier Hz          (default 4000, < 0.49*rate)
    -k N      samples per symbol  (default 20, even, 4..160)   Rs = rate/k
    -H h      modulation index    (default 0.25; 0..1 feasible range, see
                                  detector notes below)
    -b bt     Gaussian BT         (default 0.8; 0.65+ robust.  At 0.5 the
                                  ISI pattern floor makes ~1-2% of frames
                                  unrecoverable even noiseless - some
                                  waveforms admit no error-free decision
                                  alignment at all (probe5), MLSE would be
                                  needed; < 0.5 always needs MLSE)
    -d N      CPM filter delay, symbols (default 1, 1..8)
    -P N      preamble symbols    (default 72, max 65536)
       --pulse gmsk|rcos|square   (default gmsk)
       --fec N RS parity bytes 16|32|64|84 (default 32, see Wire format)

Tone frequencies: fc + (2s+1-M)*h*Rs/2, i.e. fc + {-3,-1,+1,+3} * h*Rs/2.
Defaults: 2400 baud, 4800 bps gross, Carson occupied bandwidth 6.6 kHz.
At fec=32 a max-size frame carries 217 payload bytes in 1116 symbols
(2232 bits): ~78 % end-to-end payload efficiency (RS alone costs 32/255).

The receiver picks its symbol detector from h: below h*(M-1) = 0.95 it
decodes the phase advance across each symbol (works for any h*(M-1) < 2,
but the +-3 tones must not wrap past +-pi: h = 1/3 for M = 4); at or above
that threshold it switches to a non-coherent tone bank (one k-sample
complex correlation per tone, peak bin wins - fully orthogonal at h = 1,
immune to CFO since the decision is a magnitude). h = 1 (e.g. fc=7000) and
h = 1 with k=40 both round-trip byte-exact, clean and through
CFO/noise/gain/ppm impairment sweeps.

Bandwidth: Carson occupied width = (3h + 2) * Rs (peak deviation
1.5*h*Rs). The RX validates both band edges against the mixer filter
(cut = 1.5*h*Rs + 0.85*Rs must fit below fc and below rate/2 - fc). At
Rs=2400: h=0.25 (6.6 kHz) fits fc=4000 (default); h=0.5 (8.4 kHz) wants
fc >= ~4300; h=1.0 (12 kHz) wants fc >= ~6100 or a lower symbol rate
(k=40, Rs=1200 fits fc=4000-5000).

## Achieved tolerances (tests/roundtrip.sh, all default PHY unless noted)

    payload sizes        0, 1..217, multi-frame 218..5000 (fec=32 chunks)
    FEC tiers            small + max payload round-trip at P=16/64/84;
                         6..40 ms deterministic mid-frame dropouts repaired
                         byte-exact (tier 16/32/84), > t damage rejects
    CFO                  +-295 Hz full tolerance (estimator range rate/(8k))
    sample-rate error    +-1000 ppm (offline: fixed 23-point trial grid
                         0..+-1200 ppm, CRC-good frames merged across
                         trials - no single ppm need decode all frames)
    delay                +-integer samples and fractional (linear interp)
    gain                 -30 dB .. +6 dB (clipping) — scale loop absorbs
    AWGN                 sigma = 0.03 FS passes (single + combos)
    combos               CFO+delay+ppm+noise+gain together, CFO+ppm+noise,
                         CFO+frac-delay+noise across 5-frame streams
    pattern payloads     zeros / repeated byte / JPEG chunks round-trip,
                         incl. multi-frame at bt=0.65 (whitening must
                         absorb every single-tone run pathology)
    PHY sweep            k 16..64, bt 0.65/1.0, delay 1..8, rcos/square,
                         rate 8k/16k/22050/44100, h 0.3..1.0 (incl. tone
                         bank + ppm/impairment cases), fc 8k, preamble 120

Physically rejected by design (hard error or guaranteed CRC gate, never a
false decode): fc <= dev + 0.85*Rs at the given k on EITHER band side
(mixer image cannot be filtered: at fc=4000 and the default h=0.25 this
rejects k <= 14; k=16 is the first usable value), non-integer-friendly
band edges. The former h < 1/3 cap is gone: it was a limit of the
phase-advance detector only, now bypassed by the tone bank.

Envelope note (probe-verified PRE-FEC with variable 1024-byte frames; RS
correction plausibly moved these margins — re-probe before quoting):
every PAIR of
{ppm +-300, CFO 100, noise 0.01, gain} passes; stacking ppm+CFO+noise on
top of a signal pulled to -24 dBFS fails (offline RX fails it too — it is
a decoder envelope limit, not a streaming regression).  The stress matrix
below covers multipath, fading, reverb, AGC, dropouts, colored/burst noise,
band-limiting, quantisation and longer clock drift; `tests/stress_baseline.json`
is the authoritative record of what the current RX actually achieves.

## Stress / impairment testing

`tests/stress.py` is a reproducibility-first matrix over PHY profiles and
channel impairments.  It is deliberately opt-in: `make test` stays the fast
commit gate, while `make stress` runs the curated quick matrix and
`make stress-full` runs every profile x scenario x 4 fixed seeds.

    make stress                 # quick, ~1-2 min
    make stress-full            # full, ~6 min on a 20-core machine
    python3 tests/stress.py --list
    python3 tests/stress.py --quick --offline --filter 'p0-default/rayleigh'
    python3 tests/stress.py --full --strict          # XFAIL/NEW are failures
    python3 tests/stress.py --full --paced           # extra 5x-real-time pipe feeds
    python3 tests/stress.py --full --update-baseline # after reviewing results

Every case is `TX -> channel -> RX` and both receivers are exercised by
default.  Payloads are deterministic, frames are matched in order, and the
result rejects corrupt, phantom, duplicate, reordered, truncated, timed-out
and crashed output.  Mild deterministic scenarios are gated: every frame must
arrive.  Stochastic cases are compared against `tests/stress_baseline.json`
(per-seed frame counts); known-bad cases are reported as XFAIL, and a case
that gets worse than its baseline is FAIL.  Missing baseline entries are NEW;
`--update-baseline` requires `--full` and refuses to write if any gated case
fails.  Hard and gated failures save the payload/audio/decoded streams and
exact commands in a printed temp directory; use `--keep-failures DIR` to keep
XFAIL/regression artifacts too and to choose the directory.

Most XFAILs in the current baseline are stochastic cases (Rayleigh,
frequency-selective fading, diffuse reverb, strong interference) where the
current envelope simply loses frames; the baseline records the expected count.
The old p2/p4 cross-trial duplicate cases are now gated regressions: the merge
window is derived from the lag sweep and capped at 25% of the minimum legal
frame spacing, so copies of one frame collapse while adjacent identical frames
remain distinct (`-v` reports the merge tolerance and suppressed copies).

Offline no-decode work is bounded by a per-burst and per-trial acquisition
budget (`src/fsk_rx.c`; `-v` prints `acq=` per trial).  A continuous tone,
fade or reverb tail cannot trigger dozens of MF/prefix sweeps; a real frame
after the normal 40 ms gap re-arms the gate.  The budget is applied to
expensive `k >= 32` sweeps; lower-k sweeps are cheap and are left uncapped
so sparse-frame stochastic/ppm-ramp cases keep their baseline.  The
pathological p2 cases run in single-digit or low-tens of seconds; the full
2 KB (10-frame) p2 case is opt-in via `make speed-full`.

`rx_speed.sh` also has a live-readiness case: 5 s of full-scale random noise
followed by a 40 ms gap and one max-payload p1-highrate frame, fed strictly at
1x.  The stream RX must decode that packet and finish within the feed
duration plus a small margin, proving it does not fall behind on pure noise.

`tests/channel.c` now models a fuller audio path (short options remain
backward-compatible):

    --multipath MS:DB[,MS:DB...]      static echoes (optional :PHASE_DEG)
    --reverb RT60:WET_DB              direct + diffuse room tail
    --fade-doppler HZ                 flat/frequency-selective fading
    --fade-rician K_DB                Rician K on the direct tap
    --fade-spread MS --fade-taps N    exponential power-delay profile
    --snr-db DB                       wideband SNR from active-signal RMS
    --noise-color white|pink|brown
    --bandpass LO:HI                  fourth-order Butterworth
    --hum HZ:DBFS  --tone HZ:DBFS     electrical hum / in-band interferer
    --clicks RATE:AMP_DBFS:DECAY_MS
    --bursts RATE:LEN_MS:AMP_DBFS
    --dropouts RATE:LEN_MS:MODE       Poisson zero/drop
    --drop-at MS:LEN_MS[:MODE]        deterministic dropout
    --gain-steps MS:DB[,MS:DB...]     level steps
    --ppm-ramp START:END              linear clock drift
    --agc TARGET_DBFS:ATTACK_MS:RELEASE_MS[:MAX_GAIN_DB]
    --soft-clip DRIVE --quantize BITS --dc OFFSET
    --report                          measured RMS/SNR/fade/clip stats

`--snr-db` is wideband (full-Nyquist) SNR referenced to the active signal RMS
at the noise-injection point; the in-band detector SNR is higher because the
RX low-pass removes out-of-band noise.  All stochastic parts are seeded by
`-s`; the same command gives byte-identical output.

Real captures can be added under `tests/fixtures/<name>/`: `case.json`
(`phy`, `rate`, `rx`, `expect=all|subsequence|none`, `min_frames`,
`timeout_s`), `audio.raw` (s16le mono; `.raw` is used because `.gitignore`
excludes `*.s16`), and `payload.bin`.  No fixtures are checked in.

## Streaming reception (fsk_rx_stream)

    arecord -D pipewire -f S16_LE -r 48000 -c 1 -B 2000000 | ./bin/fsk_rx_stream > out.bin

Decodes as audio arrives — no EOF, no whole-file rescan. An energy gate
(causal DDC + single-pass LPF reference) closes 20 ms after the carrier
drops (the TX inter-frame gap is 40 ms); the burst's RAW samples are then
re-run through the exact offline chain (zero-phase double pass) and the
shared rx_run() detector on a self-contained slice. Payload bytes go to
stdout per frame, whole, only after the CRC passes: latency per frame =
fixed frame airtime (1116 syms, payload-independent: 0.47 s at the default
PHY) + 40 ms gap + decode (typ. 20-100 ms; full-grid worst case a
few hundred ms, still under the arrival time of the next frame).

Differences from the offline receiver:

  * ppm trials run per burst: the coarse grid {0, +-300, +-600, +-900,
    +-1200} ordered with the previous burst's winner first, then (only if
    nothing decoded) a +-150/-100/-50 sweep centred on the coarse winner.
    Early exit on the first fully clean trial, so steady state is one
    trial per burst. Inter-burst clock drift is absorbed by the per-burst
    MF/prefix re-alignment, so the per-burst grid matches the offline
    +-1200 ppm tolerance.
  * sliding flat window (~2 max-size frames) trimmed to `retain` samples
    after every decode and slid forward during silence (no frame can
    start before the read cursor).
  * joining mid-frame skips the partial burst (bounded reject sweeps,
    then the gate idles to the next gap) and locks onto the first
    complete frame; frames cut by EOF never deliver.
  * exit code: 0 if >= 1 frame decoded, 1 otherwise — same rule for
    clean EOF and Ctrl-C. Options: --max-frames N, --idle-exit SEC.
  * single-threaded: decode runs between reads; the pipe plus arecord's
    -B buffer absorb decode stalls (a decode never exceeds one frame's
    arrival time on normal hardware).
  * regression: tests/stream.sh (fast + paced feeds, impairments, joins,
    idle gaps, --max-frames) is part of make test;
    tests/stream_loopback.sh is the manual real-device e2e (snd-aloop or
    two hosts) and is deliberately not run by make test.

## Detector design notes

Full rationale lives in the source comments; the load-bearing facts:

1. RX is a full-symbol differential detector (arg(y[t+(k-1)]·conj(y[t])) ->
   nearest pi*h*{+-1,+-3}), not liquid's cpfskdem (unreliable for M>2).
2. The LPF after DDC runs twice (filtfilt-style) so its group delay is zero
   and only |H| shapes the signal; its half-length sets the TX trailing pad
   (`nsfsk_lpf_taps` / `nsfsk_tx_pad` share the rule).
3. CFO comes from the phase slope of consecutive coherent sub-window matched
   filter sums. Preamble autocorrelation does NOT work: the alternating 0,3
   pattern folds the outer tones and mirrors the CFO phase.
4. Decision thresholds use a decision-directed scale loop: the prefix-
   calibrated scale measures the alternating-preamble ISI (~25 % low vs.
   random data), which would clip the +-1 clusters. Frames failing RS+CRC
   re-decode (<= 3 passes) seeded from the converged scale — kept for the
   scale-convergence case, not for error repair (that is RS's job now).
5. After the CFO hypothesis is removed, |phi| > pi - 0.03 triggers a +-2*pi
   wrap-undo; the outer tones plus CFO sit close to pi by the non-coherent
   limit 3*pi*h*(k-1)/k + 2*pi*|cfo|*(k-1)/rate < pi.

## Files

    src/common.{h,c}   wire format + frame codec, NRZI/Gray/CRC, params,
                       LPF+pad rule
    src/rs.{h,c}       shortened systematic RS(255,255-P), GF(256)/0x11d
    src/rx_core.{h,c}  shared detector: filter slice + rx_run() burst scan
    src/fsk_tx.c       transmitter
    src/fsk_rx.c       offline receiver (pipeline comment has the design)
    src/fsk_rx_stream.c live receiver (arecord | ...) + energy gate
    tests/rs_selftest.c RS codec self-test, all tiers (make test)
    tests/roundtrip.sh regression battery (make test)
    tests/stream.sh    streaming regression battery (make test)
    tests/channel_selftest.py fast channel property tests (make test)
    tests/rx_speed.sh  offline no-decode + live-noise RX checks (make test)
    tests/stress.py    impairment/PHY matrix + baseline (make stress)
    tests/stress_baseline.json checked-in per-seed frame counts (XFAIL tracking)
    tests/fixtures/    optional real-capture fixtures (none checked in)
    tests/stream_loopback.sh manual audio-device e2e (NOT in make test)
    tests/channel.c    impairment simulator: old short flags plus multipath,
                       reverb, fading, colored/burst noise, bandpass, AGC,
                       dropouts, clock ramp, quantisation, report
    tests/modem_smoke.c baseband loopback BER sanity check
    tests/probe_lpf.c  LPF config sweep diagnostic
    tests/probe12.c    per-symbol differential phase dump (RX debugging)
