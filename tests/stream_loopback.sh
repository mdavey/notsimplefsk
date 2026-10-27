#!/bin/sh
# Manual end-to-end loopback test for fsk_rx_stream — NOT part of make
# test (it needs a real audio path).
#
# Two hosts (true acoustic/cable path):
#   A:  cat payload.bin | ./bin/fsk_tx | aplay -q -f S16_LE -r 48000 -c 1 -
#   B:  arecord -q -D default -f S16_LE -r 48000 -c 1 -B 2000000 \
#           | ./bin/fsk_rx_stream > out.bin
#   then: cmp payload.bin out.bin
#
# One host, two devices (e.g. line-out of one sound card into line-in of
# another, or a loopback mux — modprobe snd-aloop if that's the path):
#   ./tests/stream_loopback.sh [capture-device] [playback-device] [payload-file]
#   capture default: sysdefault:CARD=M2      (arecord)
#   playback default: sysdefault:CARD=m2x4HD (aplay)
# List names with: arecord -l / aplay -l
#
# PHY flags must match on both sides if you override the defaults.
set -u
cd "$(dirname "$0")/.." || exit 1

command -v arecord >/dev/null 2>&1 || { echo "SKIP: no arecord"; exit 0; }
command -v aplay   >/dev/null 2>&1 || { echo "SKIP: no aplay";   exit 0; }

DEV_REC=${1:-sysdefault:CARD=M2}
DEV_PLAY=${2:-sysdefault:CARD=m2x4HD}
PAY=${3:-/tmp/nsfsk-loop.bin}
OUT=/tmp/nsfsk-loop-out.bin
LOG=/tmp/nsfsk-loop.log
SECS=15

[ -f "$PAY" ] || gawk -v n=5000 \
    'BEGIN{srand(7);for(i=0;i<n;i++)printf "%c",int(rand()*256)}' > "$PAY"

# capture SECS of the loopback device, decoding live
timeout $SECS arecord -q -D "$DEV_REC" -f S16_LE -r 48000 -c 1 -B 2000000 2>/dev/null \
    | ./bin/fsk_rx_stream > "$OUT" 2>"$LOG" &
RXJOB=$!

sleep 1                                   # let arecord spin up
echo "loopback: $DEV_PLAY -> $DEV_REC, $(stat -c%s "$PAY") byte(s)"
if ! ./bin/fsk_tx < "$PAY" | aplay -q -D "$DEV_PLAY" -f S16_LE -r 48000 -c 1 -; then
    kill $RXJOB 2>/dev/null
    echo "FAIL: aplay into $DEV_PLAY failed (aplay -l to list devices)"
    exit 1
fi
wait $RXJOB 2>/dev/null
echo "rx: $(tail -1 "$LOG")"
if cmp -s "$PAY" "$OUT"; then
    echo "LOOPBACK PASS ($(stat -c%s "$OUT") bytes)"
else
    echo "LOOPBACK FAIL: got $(stat -c%s "$OUT") of $(stat -c%s "$PAY") bytes"
    exit 1
fi
