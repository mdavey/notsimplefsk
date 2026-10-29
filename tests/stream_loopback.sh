#!/bin/sh
# Manual end-to-end loopback test for fsk_rx_stream — NOT part of make
# test (it needs a real audio path; aplay/arecord/PipeWire scheduling is
# exactly the thing under test, so no simulated substitute covers it).
#
# One host, one audio path (virtual sink + its source under PipeWire,
# snd-aloop, or line-out into line-in):
#   ./tests/stream_loopback.sh [capture-device] [playback-device] [payload-file]
#   capture default: pipewire   (arecord -D pipewire)
#   playback default: pipewire  (aplay -D pipewire)
# List names with: arecord -l / aplay -l
#
# Passes (REPS full + 1 join, REPS env, default 2):
#   full: payload -> fsk_tx -> aplay -> [loopback] -> arecord
#         -> fsk_rx_stream; stdout must equal the payload byte-exact.
#   join: aplay starts FIRST; the RX starts ~1 s into the stream, misses
#         frame 1 mid-air and must still deliver frames.  Connecting a
#         capture stream to a running PipeWire graph settles over a few
#         seconds and can swallow WHOLE frames anywhere in the stream
#         (observed: frames 1-3 eaten by settle, frame 4 dropped by a
#         mid-capture graph resync).  Contiguity is therefore NOT asserted;
#         the hard contract is: >= 2 frames delivered, the LAST payload
#         frame arrives, and the delivered bytes are an exact in-order
#         subsequence of payload frames (nothing corrupted, nothing
#         phantom).  Contiguous end-to-end delivery is pass_full's job.
#
# Two hosts (true acoustic/cable path): run the two halves from the
# header of the git history, or compare against a capture file directly.
# PHY flags must match on both sides if you override the defaults.
# RATE (sample rate used by arecord/aplay, default 48000) must match the
# PHY -r; e.g. RATE=44100 PHY="-r 44100 -c 3800 ...".
set -u
cd "$(dirname "$0")/.." || exit 1

command -v arecord >/dev/null 2>&1 || { echo "SKIP: no arecord"; exit 0; }
command -v aplay   >/dev/null 2>&1 || { echo "SKIP: no aplay";   exit 0; }

DEV_REC=${1:-pipewire}
DEV_PLAY=${2:-pipewire}
PAY=${3:-/tmp/nsfsk-loop.bin}
AUD=/tmp/nsfsk-loop-audio.raw
OUT=/tmp/nsfsk-loop-out.bin
LOG=/tmp/nsfsk-loop.log
REPS=${REPS:-2}
RXV=${RXV:--v}           # -v per-frame lines; -vv adds gate + ppm trace
RATE=${RATE:-48000}      # must match the PHY -r (default 48000)
PHY=${PHY:-}             # extra PHY flags for both fsk_tx and fsk_rx_stream

# LC_ALL=C: gawk "%c" emits 2 bytes for codes 128..255 in a UTF-8 locale
[ -f "$PAY" ] || LC_ALL=C gawk -v n=6144 \
    'BEGIN{srand(7);for(i=0;i<n;i++)printf "%c",int(rand()*256)}' > "$PAY"
: > "$LOG"
./bin/fsk_tx $PHY < "$PAY" > "$AUD" 2>/dev/null || { echo "FAIL: tx"; exit 1; }
PAYLEN=$(stat -c%s "$PAY")
ABYTES=$(stat -c%s "$AUD")
# capture window: audio seconds (2 B/sample) + RX spin-up + buffer drain
SECS=$(( ABYTES / (2 * RATE) + 4 ))

echo "loopback: $DEV_PLAY -> $DEV_REC, $PAYLEN payload bytes, $((ABYTES/(2*RATE)))s audio, ${SECS}s window"

pass_full() {
    timeout $SECS arecord -q -D "$DEV_REC" -f S16_LE -r "$RATE" -c 1 -B 2000000 \
        2>/dev/null | ./bin/fsk_rx_stream $RXV $PHY > "$OUT" 2>>"$LOG" &
    RXJOB=$!
    sleep 1                                   # let arecord spin up
    if ! aplay -q -D "$DEV_PLAY" -f S16_LE -r "$RATE" -c 1 - < "$AUD"; then
        kill $RXJOB 2>/dev/null
        echo "    aplay into $DEV_PLAY failed (aplay -l to list devices)"
        return 2
    fi
    wait $RXJOB 2>/dev/null
    rc=$?
    got=$(stat -c%s "$OUT")
    [ $rc -eq 0 ] || return 1
    cmp -s "$PAY" "$OUT"
}

pass_join() {
    aplay -q -D "$DEV_PLAY" -f S16_LE -r "$RATE" -c 1 - < "$AUD" 2>/dev/null &
    APJOB=$!
    sleep 1                                   # join mid frame 1
    timeout $SECS arecord -q -D "$DEV_REC" -f S16_LE -r "$RATE" -c 1 -B 2000000 \
        2>/dev/null | ./bin/fsk_rx_stream $RXV $PHY > "$OUT" 2>>"$LOG"
    rc=$?
    wait $APJOB 2>/dev/null
    got=$(stat -c%s "$OUT")
    [ $rc -eq 0 ] && [ "$got" -ge 2048 ] || return 1
    # last payload frame must have arrived intact (final frame may be short)
    nfr=$(( (PAYLEN + 1023) / 1024 ))
    lsz=$(( PAYLEN - 1024*(nfr-1) ))
    tail -c $lsz "$PAY" > "$OUT.exp"
    tail -c $lsz "$OUT" | cmp -s - "$OUT.exp" || return 1
    # greedy: OUT must be an in-order concatenation of whole payload frames
    qi=0; pi=0
    while [ "$qi" -lt "$got" ]; do
        k=$pi; found=0
        while [ "$k" -lt "$nfr" ]; do
            fsz=1024; [ "$k" -eq $((nfr-1)) ] && fsz=$(( PAYLEN - 1024*(nfr-1) ))
            dd if="$PAY" bs=1024 skip=$k count=1 2>/dev/null | head -c $fsz > "$OUT.pk"
            dd if="$OUT" bs=1 skip=$qi count=$fsz 2>/dev/null > "$OUT.qi"
            if cmp -s "$OUT.pk" "$OUT.qi"; then found=1; break; fi
            k=$((k+1))
        done
        [ $found -eq 1 ] || return 1
        qi=$(( qi + fsz )); pi=$(( k + 1 ))
    done
    return 0
}

FAILS=0
i=1
while [ "$i" -le "$REPS" ]; do
    if pass_full; then
        echo "PASS full#$i ($PAYLEN bytes byte-exact)"
    else
        echo "FAIL full#$i (got $(stat -c%s "$OUT") of $PAYLEN bytes)"
        tail -2 "$LOG" | sed 's/^/    /'
        FAILS=$((FAILS+1))
    fi
    i=$((i+1))
    sleep 0.5                                 # let the sink go quiet again
done
if [ "$PAYLEN" -lt 2048 ]; then
    echo "SKIP join (payload < 2 frames)"
elif pass_join; then
    echo "PASS join (mid-stream RX start, $(stat -c%s "$OUT") bytes, in-order subsequence)"
else
    echo "FAIL join (got $(stat -c%s "$OUT") of $PAYLEN bytes)"
    tail -2 "$LOG" | sed 's/^/    /'
    FAILS=$((FAILS+1))
fi

echo "rx last: $(grep '^rx: [0-9]* frame' "$LOG" | tail -1)"
[ "$FAILS" -eq 0 ] || { echo "LOOPBACK FAIL ($FAILS pass(es), log: $LOG)"; exit 1; }
echo "LOOPBACK PASS ($((REPS+1)) passes)"
