#!/bin/sh
# Manual ACOUSTIC coupling test for fsk_rx_stream: loudspeaker -> room ->
# microphone.  NOT part of make test (physical path; user runs by hand).
#
#   PHY="-r 48000 -c 1000 -k 160 -H 1.0 -b 1.0 -d 4 -P 72" \
#       tests/stream_acoustic.sh [payload_bytes]        (default 1100)
#
# Requires: aplay's device is audible from the mic that arecord reads.
# Defaults use the PipeWire default sink/source (-D pipewire).
#
# Gain staging rules learned 2026-10-28 (300 baud, h=1.0):
#  - room noise floor must sit well under the signal; at ~4 dB SNR nothing
#    locks (sync never found), at ~40 dB it decodes fine
#  - a -6 dBFS capture limiter is harmless (constant-envelope CPFSK clips
#    clean; pinned by sim "clip@16k" reasoning) - do not fear gain
#  - sink and source are separate async USB clocks: expect a real ppm
#    offset (observed +300); the grid covers +-1200
#  - keep the mic CLOSE to the speaker: at 3.3 ms symbols, room reflections
#    from 1-2 m distances are 1 symbol period of destructive ISI
set -u
cd "$(dirname "$0")/.." || exit 1

command -v arecord >/dev/null 2>&1 || { echo "SKIP: no arecord"; exit 0; }
command -v aplay   >/dev/null 2>&1 || { echo "SKIP: no aplay";   exit 0; }

PHY=${PHY:--c 1000 -k 160 -H 1.0 -b 1.0 -d 4}
RATE=${RATE:-48000}      # must match the PHY -r (default 48000)
N=${1:-1100}
DEV_REC=${DEV_REC:-pipewire}
DEV_PLAY=${DEV_PLAY:-pipewire}
RXV=${RXV:--v}
PAY=/tmp/nsfsk-aco.bin
AUD=/tmp/nsfsk-aco-audio.raw
OUT=/tmp/nsfsk-aco-out.bin
LOG=/tmp/nsfsk-aco.log
CAP=${CAP:-}              # set to a path to also keep the raw capture s16le

if [ ! -f "$PAY" ] || [ "$(stat -c%s "$PAY")" -ne "$N" ]; then
    LC_ALL=C gawk -v n="$N" 'BEGIN{srand(11);for(i=0;i<n;i++)printf "%c",int(rand()*256)}' > "$PAY"
fi
./bin/fsk_tx $PHY < "$PAY" > "$AUD" 2>/dev/null || { echo "FAIL: tx"; exit 1; }
ABYTES=$(stat -c%s "$AUD")
SECS=$(( ABYTES / (2 * RATE) + 6 ))

echo "acoustic: $DEV_PLAY -> $DEV_REC, $N bytes, $((ABYTES/(2*RATE)))s audio, ${SECS}s window [$PHY]"

# capture first, TX second (the proven lead; mid-stream joins on PipeWire
# graphs drop frames - see stream_loopback.sh header)
if [ -n "$CAP" ]; then
    timeout $SECS arecord -q -D "$DEV_REC" -f S16_LE -r "$RATE" -c 1 \
        2>/dev/null | tee "$CAP" | ./bin/fsk_rx_stream $RXV $PHY > "$OUT" 2>"$LOG" &
else
    timeout $SECS arecord -q -D "$DEV_REC" -f S16_LE -r "$RATE" -c 1 \
        2>/dev/null | ./bin/fsk_rx_stream $RXV $PHY > "$OUT" 2>"$LOG" &
fi
RXJOB=$!
sleep 1
aplay -q -D "$DEV_PLAY" -f S16_LE -r "$RATE" -c 1 - < "$AUD" || \
    { kill $RXJOB 2>/dev/null; echo "FAIL: aplay"; exit 1; }
wait $RXJOB 2>/dev/null
rc=$?

if [ $rc -eq 0 ] && cmp -s "$PAY" "$OUT"; then
    echo "ACOUSTIC PASS: $(stat -c%s "$OUT") of $N bytes byte-exact"
    tail -1 "$LOG" | sed 's/^/  /'
else
    echo "ACOUSTIC FAIL (rc=$rc): got $(stat -c%s "$OUT") of $N bytes"
    tail -3 "$LOG" | sed 's/^/  /'
    exit 1
fi
