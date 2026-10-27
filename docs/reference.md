# notsimplefsk — G4FSK modem reference

Two programs move arbitrary bytes over audio using 4-level continuous-phase
FSK with Gaussian pulse shaping (G4FSK), built on liquid-dsp:

    bin/fsk_tx   bytes (stdin)  -> s16le mono audio (stdout)
    bin/fsk_rx   s16le audio    -> recovered bytes (stdout)

All diagnostics go to stderr; stdout is always pure payload/audio so the
programs compose in pipelines. RX PHY options must match the TX exactly.

## Build

    make            # bin/{fsk_tx,fsk_rx,modem_smoke,channel}
    make test       # full regression battery (tests/roundtrip.sh, ~60 cases)

Requires libliquid >= 1.6 (`-lliquid`).

## Quick start

    ./bin/fsk_tx < secret.bin > audio.s16
    ./bin/fsk_rx < audio.s16 > got.bin
    cmp secret.bin got.bin

Audition/playback: `sox audio.s16 -r 48000 -c 1 audio.wav`.
The channel simulator adds impairments for testing:

    ./bin/fsk_tx < in.bin | ./bin/channel -f 200 -r 300 -n 0.01 | ./bin/fsk_rx > out.bin

## Wire format (per frame)

    [ preamble 0,3,0,3,... x P ][ sync 24-sym table ][ Gray(NRZI( data )) ][ pad ]

    data   = type(1) | len16_be(2) | payload(0..1024) | crc32_be(4)
    crc32  = liquid LIQUID_CRC_32 over type|len|payload
    NRZI   = enc[n] = enc[n-1] XOR d[n], state 0 at first bit (kills 180-deg
             ambiguity of non-coherent detection)
    Gray   = dibit -> symbol index  v ^ (v>>1)  = {0,1,3,2}
    pad    = nsfsk_tx_pad() dummy symbols: covers the detector lag and the
             RX zero-phase LPF edge reach so the CRC trailer is never
             filtered against the end of the buffer
    sync   = non-periodic 24-symbol table; the periodic 0,3 preamble has an
             MF sidelobe every +-2 symbols, only the sync table resolves it

TX emits one frame per <= 1024-byte chunk of stdin; larger inputs become
multiple frames separated by 40 ms silence (the RX energy gate re-arms on it).

## PHY parameters (identical flags on both programs)

    -r rate   sample rate Hz      (default 48000, min 8000)
    -c fc     carrier Hz          (default 4000, < 0.49*rate)
    -k N      samples per symbol  (default 20, even, 4..160)   Rs = rate/k
    -H h      modulation index    (default 0.25, must have h*(M-1) < 1)
    -b bt     Gaussian BT         (default 0.8; < 0.5 has an ISI floor)
    -d N      CPM filter delay, symbols (default 1, 1..8)
    -P N      preamble symbols    (default 72)
       --pulse gmsk|rcos|square   (default gmsk)

Tone frequencies: fc + (2s+1-M)*h*Rs/2, i.e. fc + {-3,-1,+1,+3} * h*Rs/2.
Defaults: 2400 baud, 4800 bps, Carson occupied bandwidth 6.6 kHz, ~97 %
payload efficiency at 1024-byte frames.

## Achieved tolerances (tests/roundtrip.sh, all default PHY unless noted)

    payload sizes        0, 1..1023, 1024, multi-frame 1025..5000
    CFO                  +-295 Hz full tolerance (estimator range rate/(8k))
    sample-rate error    +-1000 ppm (offline resample trials +-1200 + fine)
    delay                +-integer samples and fractional (linear interp)
    gain                 -30 dB .. +6 dB (clipping) — scale loop absorbs
    AWGN                 sigma = 0.03 FS passes (single + combos)
    combos               CFO+delay+ppm+noise+gain together, CFO+ppm+noise,
                         CFO+frac-delay+noise across 5-frame streams
    PHY sweep            k 16..64, bt 0.5/1.0, delay 1..8, rcos/square,
                         rate 8k/16k/22050/44100, h 0.3, fc 8k, preamble 120

Physically rejected by design (hard error or guaranteed CRC gate, never a
false decode): fc <= dev + 0.85*Rs at the given k (mixer image cannot be
filtered: k <= 10 @ fc=4000), h >= 1/3, non-integer-friendly band edges.

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
   random data), which would clip the +-1 clusters. Failed CRCs re-decode
   (<= 3 passes) seeded from the converged scale.
5. After the CFO hypothesis is removed, |phi| > pi - 0.03 triggers a +-2*pi
   wrap-undo; the outer tones plus CFO sit close to pi by the non-coherent
   limit 3*pi*h*(k-1)/k + 2*pi*|cfo|*(k-1)/rate < pi.

## Files

    src/common.{h,c}   wire format, NRZI/Gray/CRC, params, LPF+pad rule
    src/fsk_tx.c       transmitter
    src/fsk_rx.c       receiver (pipeline comment at top has the full design)
    tests/roundtrip.sh regression battery (make test)
    tests/channel.c    impairment simulator: -d int delay, -D fractional,
                       -r ppm, -g dB gain, -f Hz CFO (Hilbert SSB), -n sigma
                       AWGN, -s seed, -R rate
    tests/modem_smoke.c baseband loopback BER sanity check
    tests/probe_lpf.c  LPF config sweep diagnostic
    tests/probe12.c    per-symbol differential phase dump (RX debugging)
