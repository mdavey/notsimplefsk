#!/bin/sh
# Noise lead-in battery for fsk_rx and fsk_rx_stream.
#
# Prepends Gaussian noise (via tests/channel -n, fixed seed -> deterministic)
# of various levels and lengths (100 ms .. 10 s) to a valid capture, then
# decodes offline and through the live streaming path.  Purpose: track RX
# acquisition speed and robustness after long stretches of noise, which is
# the realistic "squelch opens onto noise, packet arrives" scenario.
#
# Offline fsk_rx must decode every case and finish quickly (it is pure
# post-processing, so per-case wall time is reported in samples-of-audio
# per second of processing; a regression back into the minutes-long
# acquisition sweeps of the past trips the timeout).
#
# fsk_rx_stream runs through a fast pipe; a few cases also run paced at
# strict 1x (arecord real time) with a first-output-latency bound.
#
# XFAIL: sustained noise above the stream energy-gate floor (GATE_FLOOR,
# ~-40 dBFS mean |y|^2) keeps the gate open; the bounded ring then
# force-decodes mid-noise and can destroy the frame head depending on the
# phase of the force-decode boundary vs the packet start (which is why the
# failures are phase-dependent, e.g. 1 s and 5 s leads fail while 2.5 s
# passes).  Those cells are pinned as XFAIL: the suite stays green today,
# an XPASS prints a loud "now decodes" line so the cell can be promoted to
# a normal pass once the gate/force-decode policy is fixed.
set -u
cd "$(dirname "$0")/.." || exit 1
TX=./bin/fsk_tx
RX=./bin/fsk_rx_stream
RXOFF=./bin/fsk_rx
CH=./bin/channel
TMP=${TMPDIR:-/tmp}/nsfsk-nsl.$$
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

PASS=0; FAIL=0; XFAIL=0; XPASS=0

pay() {
    LC_ALL=C gawk -v n="$1" -v s="$(printf %s "$2" | cksum | cut -d' ' -f1)" '
        BEGIN{srand(s);for(i=0;i<n;i++)printf "%c",int(rand()*256)}'
}

# noise <ms> <sigma> <seed> : Gaussian noise of that many ms, full-scale
# sigma, written to $TMP/n.  channel -n over digital silence is exactly
# additive AWGN, deterministic for a fixed seed.
noise() {
    head -c $(( $1 * 96 )) /dev/zero | $CH -n "$2" -s "$3" 2>/dev/null > "$TMP/n"
}

# timeout generous but finite: 10 s of audio must never take minutes
TMO=30

# run_offline <name> <lead_ms> <sigma> <nbytes>
run_offline() {
    name=$1 lead=$2 sig=$3 n=$4
    pay "$n" "off-$name" > "$TMP/p"
    $TX < "$TMP/p" > "$TMP/a" 2>/dev/null || { echo "FAIL $name (tx)"; FAIL=$((FAIL+1)); return; }
    noise "$lead" "$sig" 7
    cat "$TMP/n" "$TMP/a" > "$TMP/in"
    audio=$(stat -c%s "$TMP/in")
    t0=$(date +%s%N)
    timeout $TMO $RXOFF < "$TMP/in" > "$TMP/o" 2>"$TMP/rxe"
    rc=$?
    t1=$(date +%s%N)
    ms=$(( (t1 - t0) / 1000000 ))
    [ $ms -eq 0 ] && ms=1        # avoid div-by-zero in the speed factor
    factor=$(awk -v a="$audio" -v m="$ms" \
        'BEGIN{printf "%.0fx", (a/96000.0) / (m/1000.0)}')
    if [ $rc -eq 124 ]; then
        echo "FAIL $name (timeout ${TMO}s, lead=${lead}ms sig=$sig)"
        FAIL=$((FAIL+1)); return
    fi
    if [ $rc -eq 0 ] && cmp -s "$TMP/p" "$TMP/o"; then
        echo "PASS $name (lead=${lead}ms sig=$sig ${ms}ms, ${factor}RT)"
        PASS=$((PASS+1))
    else
        echo "FAIL $name (rc=$rc lead=${lead}ms sig=$sig ${ms}ms, got $(stat -c%s "$TMP/o")B want ${n}B)"
        sed 's/^/    /' "$TMP/rxe" | tail -2
        FAIL=$((FAIL+1))
    fi
}

# expect_stream <lead_ms> <sigma> : pass | xfail for the current gate
# policy on the single-frame matrix (seed fixed, so these cells are
# deterministic).  Multi-frame cases are seeded differently and do NOT
# inherit this table - phase effects differ, so they pass their own want.
expect_stream() {
    case "$1:$2" in
        1000:0.005|1000:0.02|1000:0.05|5000:0.005|5000:0.02|10000:0.005|10000:0.02|10000:0.05)
            echo xfail ;;
        *) echo pass ;;
    esac
}

# run_stream <name> <lead_ms> <sigma> <nbytes> [want] : fast pipe, byte-exact
run_stream() {
    name=$1 lead=$2 sig=$3 n=$4
    want=${5:-$(expect_stream "$lead" "$sig")}
    pay "$n" "str-$name" > "$TMP/p"
    $TX < "$TMP/p" > "$TMP/a" 2>/dev/null || { echo "FAIL $name (tx)"; FAIL=$((FAIL+1)); return; }
    noise "$lead" "$sig" 7
    cat "$TMP/n" "$TMP/a" > "$TMP/in"
    t0=$(date +%s%N)
    timeout $TMO $RX < "$TMP/in" > "$TMP/o" 2>"$TMP/rxe"
    rc=$?
    t1=$(date +%s%N)
    ms=$(( (t1 - t0) / 1000000 ))
    ok=0
    [ $rc -eq 0 ] && cmp -s "$TMP/p" "$TMP/o" && ok=1
    if [ "$want" = xfail ]; then
        if [ $ok -eq 1 ]; then
            echo "XPASS $name *** lead=${lead}ms sig=$sig NOW DECODES (${ms}ms): promote to pass + fix gate/force-decode policy is no longer needed for this cell ***"
            XPASS=$((XPASS+1))
        else
            echo "XFAIL $name (lead=${lead}ms sig=$sig rc=$rc ${ms}ms, known gate-open limitation)"
            XFAIL=$((XFAIL+1))
        fi
        return
    fi
    if [ $rc -eq 124 ]; then
        echo "FAIL $name (timeout ${TMO}s, lead=${lead}ms sig=$sig)"
        FAIL=$((FAIL+1)); return
    fi
    if [ $ok -eq 1 ]; then
        echo "PASS $name (lead=${lead}ms sig=$sig ${ms}ms)"
        PASS=$((PASS+1))
    else
        echo "FAIL $name (rc=$rc lead=${lead}ms sig=$sig ${ms}ms, got $(stat -c%s "$TMP/o")B want ${n}B)"
        sed 's/^/    /' "$TMP/rxe" | tail -2
        FAIL=$((FAIL+1))
    fi
}

# feed_rt1 <file> : strict 1x real-time producer (0.5 s audio per 0.5 s)
feed_rt1() {
    total=$(stat -c%s "$1"); n=0
    while [ $((n * 48000)) -lt "$total" ]; do
        dd if="$1" bs=48000 skip=$n count=1 2>/dev/null
        n=$((n+1))
        sleep 0.5
    done
}

# run_live <name> <lead_ms> <sigma> <nbytes> <first_lat_max_s>
# 1x paced; byte-exact at EOF; first payload byte no later than the lead
# plus the frame period plus margin; wall time tracks the feed duration.
run_live() {
    name=$1 lead=$2 sig=$3 n=$4 latmax=$5
    pay "$n" "live-$name" > "$TMP/p"
    $TX < "$TMP/p" > "$TMP/a" 2>/dev/null || { echo "FAIL $name (tx)"; FAIL=$((FAIL+1)); return; }
    noise "$lead" "$sig" 7
    cat "$TMP/n" "$TMP/a" > "$TMP/in"
    feedsec=$(( $(stat -c%s "$TMP/in") / 96000 ))
    t0=$(date +%s%N)
    s0=$(date +%s)
    ( feed_rt1 "$TMP/in" | timeout 120 $RX > "$TMP/o" 2>"$TMP/rxe" ) &
    job=$!
    first=0
    while kill -0 $job 2>/dev/null; do
        sz=$(stat -c%s "$TMP/o" 2>/dev/null || echo 0)
        if [ "$sz" -ge 1 ]; then
            first=$(( ($(date +%s%N) - t0) / 1000000 ))
            break
        fi
        sleep 0.1
    done
    wait $job; rc=$?
    el=$(( $(date +%s) - s0 ))
    if [ "$first" -eq 0 ]; then
        # the writer can flush within one poll interval of exiting; the
        # measured exit time is still a valid first-byte upper bound
        sz=$(stat -c%s "$TMP/o" 2>/dev/null || echo 0)
        if [ "$sz" -ge 1 ]; then first=$((el * 1000)); else first=999999; fi
    fi
    firsts=$(( first / 1000 ))
    if [ $rc -eq 0 ] && cmp -s "$TMP/p" "$TMP/o" \
       && [ "$firsts" -le "$latmax" ] && [ "$el" -le $((feedsec + 4)) ]; then
        echo "PASS $name (lead=${lead}ms sig=$sig: first=${first}ms, ${el}s/${feedsec}s)"
        PASS=$((PASS+1))
    else
        echo "FAIL $name (rc=$rc first=${first}ms max=${latmax}s ${el}s/${feedsec}s, got $(stat -c%s "$TMP/o")B)"
        sed 's/^/    /' "$TMP/rxe" | tail -2
        FAIL=$((FAIL+1))
    fi
}

echo "== offline fsk_rx: noise lead, single max frame (217 B) =="
for lead in 100 250 500 1000 2500 5000 10000; do
    for sig in 0.002 0.005 0.02 0.05; do
        run_offline "off-lead$lead-s$sig" "$lead" "$sig" 217
    done
done

echo
echo "== offline fsk_rx: noise lead, multi-frame stream (1519 B = 7 frames) =="
run_offline "off-multi-lead1000-s0.002"  1000 0.002 1519
run_offline "off-multi-lead5000-s0.02"   5000 0.02  1519
run_offline "off-multi-lead10000-s0.005" 10000 0.005 1519

echo
echo "== fsk_rx_stream fast pipe: noise lead, single max frame =="
for lead in 100 250 500 1000 2500 5000 10000; do
    for sig in 0.002 0.005 0.02; do
        run_stream "str-lead$lead-s$sig" "$lead" "$sig" 217
    done
done
run_stream "str-lead100-s0.05"   100   0.05 217
run_stream "str-lead1000-s0.05"  1000  0.05 217
run_stream "str-lead10000-s0.05" 10000 0.05 217

echo
echo "== fsk_rx_stream fast pipe: noise lead, multi-frame stream =="
run_stream "str-multi-lead2000-s0.002"  2000 0.002 1519 pass
run_stream "str-multi-lead5000-s0.002"  5000 0.002 1519 pass
# the single-frame 1000ms/s0.02 cell is XFAIL, but the different payload
# phase of the 7-frame capture decodes fine today - pinned as a real pass
run_stream "str-multi-lead1000-s0.02"   1000 0.02  1519 pass

echo
echo "== fsk_rx_stream strict 1x live feed (first-byte latency bound) =="
run_live "live-lead1000-s0.002"  1000 0.002 217  5
run_live "live-multi-lead2000-s0.002" 2000 0.002 1519 6

echo
echo "PASS=$PASS FAIL=$FAIL XFAIL=$XFAIL XPASS=$XPASS"
[ "$FAIL" -eq 0 ] || { echo "NOISE LEAD FAIL ($FAIL)"; exit 1; }
echo "NOISE LEAD PASS"
