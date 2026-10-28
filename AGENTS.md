# Notes for agents

## Audio loopback tests: the "pipewire" device is not one fixed path

`make loopback`, `make acoustic` and `tests/stream_loopback.sh` /
`tests/stream_acoustic.sh` pass the literal device name `pipewire` to both
`aplay` and `arecord`. That resolves to whatever the PipeWire *default*
sink/source happen to be at run time — and those change between sessions
and can even flip mid-session:

- sometimes a **virtual** sink/source pair (digital monitor routing,
  transparent, no room acoustics, levels always good);
- sometimes the **physical** path (speaker -> room -> microphone, e.g.
  MiniFuse 2 Pro input with miniDSP 2x4HD as sink), where gain staging,
  room noise and anything else playing on the audio system matter.

That is why the first loopback battery run after the 2026-10-28 rx_core
work failed: the active route was carrying a steady 4 kHz sine at ~-6
dBFS (the captured "FSK" had constant instantaneous frequency - the
modulation never made it through), and it passed again as soon as the
route was clean. Same run, same binary, different audio graph.

Practical consequences when a loopback/acoustic test fails:

1. A/B against a known-good binary on a saved capture before touching RX
   code - if the *old* binary fails identically, it is the audio path,
   not the code.
2. Check what the captured signal actually is: real CPFSK is
   constant-envelope with the zero-crossing rate swinging between tone
   bins; a constant-amplitude, constant-frequency section means external
   audio or a processing device, not a decoder bug.
3. Levels: a good acoustic capture is around -17 dBFS signal over a
   ~-45 dBFS floor. Near-zero SNR means sync will never lock.
4. `tests/stream_acoustic.sh CAP=/tmp/cap.raw ...` keeps the raw capture
   for offline analysis; the tail frame of a multi-frame acoustic run is
   the marginal one - rerun once before believing a regression.
