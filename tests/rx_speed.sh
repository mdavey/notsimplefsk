#!/bin/sh
# Offline RX runtime regression checks for the no-decode acquisition budget.
#
# These are the p2 cases that used to spend minutes in repeated MF/prefix
# sweeps (strong tone, frequency-selective fade, diffuse reverb) plus the
# 2 KB Rayleigh case (10 max-payload frames).  A garbage burst followed by
# a real frame after the
# normal 40 ms gap must still decode, proving the per-burst cap re-arms.
set -u
cd "$(dirname "$0")/.." || exit 1
TX=./bin/fsk_tx
RX=./bin/fsk_rx
CH=./bin/channel
TMP=${TMPDIR:-/tmp}/nsfsk-speed.$$
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

# Defaults must be exercised regardless of a developer shell.
unset NSFSK_PPM NSFSK_RX_ACQ_TOTAL NSFSK_RX_ACQ_BURST NSFSK_RX_REJECT_LIMIT

PHY="-r 48000 -c 5000 -k 40 -H 0.5 -b 1.0 -d 4 -P 72"
FAILS=0

pay() {
    LC_ALL=C gawk -v n="$1" -v s="$(printf %s "$2" | cksum | cut -d' ' -f1)" '
        BEGIN{srand(s);for(i=0;i<n;i++)printf "%c",int(rand()*256)}'
}

# feed_rt1 <file> : strict 1x producer (0.5 s audio per 0.5 s), like a live RF
# capture.  If the decoder falls behind, the pipe backpressures and the total
# wall time stretches past the audio duration.
feed_rt1() {
    total=$(stat -c%s "$1"); n=0
    while [ $((n * 48000)) -lt "$total" ]; do
        dd if="$1" bs=48000 skip=$n count=1 2>/dev/null
        n=$((n+1))
        sleep 0.5
    done
}

# run_speed <name> <nbytes> <max_acq> <timeout_s> [channel args...]
# Cases 1-4 do not require a particular decode result (the stress baseline
# owns that); they require a bounded run and bounded acquisition work.
run_speed() {
    name=$1 n=$2 maxacq=$3 tmo=$4
    shift 4
    pay "$n" "$name" > "$TMP/p"
    if ! $TX $PHY -g -20 < "$TMP/p" > "$TMP/a" 2>"$TMP/txe"; then
        echo "FAIL $name (tx)"; sed 's/^/    /' "$TMP/txe" | tail -2
        FAILS=$((FAILS+1)); return
    fi
    if ! $CH "$@" -s 1 --report < "$TMP/a" > "$TMP/c" 2>"$TMP/che"; then
        echo "FAIL $name (channel)"; sed 's/^/    /' "$TMP/che" | tail -2
        FAILS=$((FAILS+1)); return
    fi
    t0=$(date +%s)
    timeout "$tmo" $RX $PHY -v < "$TMP/c" > "$TMP/o" 2>"$TMP/rxe"
    rc=$?
    t1=$(date +%s)
    if [ "$rc" -eq 124 ]; then
        echo "FAIL $name (timeout ${tmo}s)"
        FAILS=$((FAILS+1)); return
    fi
    if [ "$rc" -ne 0 ] && [ "$rc" -ne 1 ]; then
        echo "FAIL $name (rx rc=$rc)"
        tail -2 "$TMP/rxe" | sed 's/^/    /'
        FAILS=$((FAILS+1)); return
    fi
    acq=$(grep -o 'acq=[0-9]*' "$TMP/rxe" | cut -d= -f2 |
          awk '{s+=$1} END{print s+0}')
    if [ "$acq" -le 0 ]; then
        echo "FAIL $name (no acq diagnostics)"
        FAILS=$((FAILS+1)); return
    fi
    if [ "$acq" -gt "$maxacq" ]; then
        echo "FAIL $name (acq=$acq > $maxacq)"
        FAILS=$((FAILS+1)); return
    fi
    echo "PASS $name (acq=$acq, $((t1-t0))s)"
}

echo "== offline no-decode acquisition budget (p2) =="
run_speed "p2-strong-tone-1k" 1024 128 30 --tone 4300:-20 --snr-db 30
run_speed "p2-sel-fade-1k"    1024 128 30 --fade-doppler 2 --fade-spread 2 \
          --fade-taps 4 --snr-db 25
run_speed "p2-reverb-1k"      1024 128 30 --reverb 0.3:-6 --snr-db 30
if [ "${RX_SPEED_FULL:-0}" = "1" ]; then
    # 10 max-payload frames ~ the old 8 KB single-codec-era audio length
    run_speed "p2-rayleigh-2k" 2170 160 90 --fade-doppler 5 --snr-db 20
else
    echo "SKIP p2-rayleigh-2k (set RX_SPEED_FULL=1 to run)"
fi

echo "== garbage burst, silence, then a valid frame =="
pay 217 "garbage-then-valid" > "$TMP/p"    # exactly one frame, as originally
$TX $PHY -g -20 < "$TMP/p" > "$TMP/frame" 2>/dev/null
head -c 48000 /dev/zero | $CH --tone 4300:-20 --snr-db 30 -s 1 \
    > "$TMP/garbage" 2>/dev/null       # 0.5 s of in-band garbage
head -c 3840 /dev/zero > "$TMP/gap"    # 40 ms silence
cat "$TMP/garbage" "$TMP/gap" "$TMP/frame" > "$TMP/in"
t0=$(date +%s)
timeout 30 $RX $PHY -v < "$TMP/in" > "$TMP/o" 2>"$TMP/rxe"
rc=$?
t1=$(date +%s)
acq=$(grep -o 'acq=[0-9]*' "$TMP/rxe" | cut -d= -f2 |
      awk '{s+=$1} END{print s+0}')
if [ "$rc" -eq 0 ] && cmp -s "$TMP/p" "$TMP/o" && [ "$acq" -le 128 ]; then
    echo "PASS garbage-then-valid (acq=$acq, $((t1-t0))s)"
else
    echo "FAIL garbage-then-valid (rc=$rc acq=$acq, got $(stat -c%s "$TMP/o")B)"
    tail -3 "$TMP/rxe" | sed 's/^/    /'
    FAILS=$((FAILS+1))
fi

echo "== live random noise then one packet (p1, high baud) =="
PHY1="-r 48000 -c 8000 -k 16 -H 0.25 -b 0.8 -d 1 -P 72"
pay 217 "live-noise-packet" > "$TMP/p"    # one max-payload frame
$TX $PHY1 -g -12 < "$TMP/p" > "$TMP/frame" 2>/dev/null
head -c $((48000 * 2 * 5)) /dev/urandom > "$TMP/noise"
head -c 3840 /dev/zero > "$TMP/gap"
cat "$TMP/noise" "$TMP/gap" "$TMP/frame" > "$TMP/in"
feedsec=$(( $(stat -c%s "$TMP/in") / 96000 ))
t0=$(date +%s)
( feed_rt1 "$TMP/in" | timeout 60 $RX $PHY1 > "$TMP/o" 2>"$TMP/rxe" )
rc=$?
t1=$(date +%s)
el=$((t1-t0))
if [ "$rc" -eq 0 ] && cmp -s "$TMP/p" "$TMP/o" && [ "$el" -le $((feedsec + 4)) ]; then
    echo "PASS live-noise-packet (${el}s/${feedsec}s audio, $(stat -c%s "$TMP/o")B)"
else
    echo "FAIL live-noise-packet (rc=$rc ${el}s/${feedsec}s, got $(stat -c%s "$TMP/o")B)"
    tail -3 "$TMP/rxe" | sed 's/^/    /'
    FAILS=$((FAILS+1))
fi

echo
[ "$FAILS" -eq 0 ] || { echo "RX SPEED FAIL ($FAILS)"; exit 1; }
echo "RX SPEED PASS"
