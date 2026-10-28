#!/bin/sh
# notsimplefsk streaming RX regression battery (fsk_rx_stream)
# run from the repo root:  sh tests/stream.sh   (after `make`)
#
# Same byte-exact contract as roundtrip.sh, but through the live path:
# bursts close on the energy gate, ppm trials run per burst with early
# exit.  "rt" cases pace the feed with dd+sleep so decode and production
# interleave through a real pipe.
set -u
cd "$(dirname "$0")/.." || exit 1
TX=./bin/fsk_tx
RX=./bin/fsk_rx_stream
CH=./bin/channel
TMP=${TMPDIR:-/tmp}/nsfsk-st.$$
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

PASS=0; FAIL=0
PHY=""                    # extra PHY flags, must match tx+rx

# pay <nbytes> <seedname> : deterministic pseudo-random payload
pay() {
    LC_ALL=C gawk -v n="$1" -v s="$(printf %s "$2" | cksum | cut -d' ' -f1)" '
        BEGIN{srand(s);for(i=0;i<n;i++)printf "%c",int(rand()*256)}'
}

# run <name> <nbytes> [channel args...] : fast pipe, bytes in == bytes out
run() {
    name=$1 n=$2; shift 2
    pay "$n" "$name" > "$TMP/p"
    if ! $TX $PHY < "$TMP/p" > "$TMP/a" 2>"$TMP/txe"; then
        echo "FAIL $name  (tx errored)"; sed 's/^/    /' "$TMP/txe"
        FAIL=$((FAIL+1)); return
    fi
    if [ $# -gt 0 ]; then
        $CH "$@" < "$TMP/a" > "$TMP/b" 2>/dev/null
        mv "$TMP/b" "$TMP/a"
    fi
    timeout 120 $RX $PHY < "$TMP/a" > "$TMP/o" 2>"$TMP/rxe"
    rc=$?
    if [ $rc -eq 0 ] && cmp -s "$TMP/p" "$TMP/o"; then
        echo "PASS $name"; PASS=$((PASS+1))
    else
        echo "FAIL $name  (rc=$rc)"; sed 's/^/    /' "$TMP/rxe" | tail -3
        FAIL=$((FAIL+1))
    fi
}

# feed_rt <file> : ~real-time paced producer (0.2 s sleep per 0.5 s chunk)
feed_rt() {
    total=$(stat -c%s "$1"); n=0
    while [ $((n * 96000)) -lt "$total" ]; do
        dd if="$1" bs=96000 skip=$n count=1 2>/dev/null
        n=$((n+1))
        sleep 0.2
    done
}

# runrt <name> <nbytes> : paced feed through a pipe
runrt() {
    name=$1 n=$2
    pay "$n" "$name" > "$TMP/p"
    $TX $PHY < "$TMP/p" > "$TMP/a" 2>/dev/null || { echo "FAIL $name (tx)"; FAIL=$((FAIL+1)); return; }
    ( feed_rt "$TMP/a" | timeout 180 $RX $PHY > "$TMP/o" 2>"$TMP/rxe" )
    rc=$?
    if [ $rc -eq 0 ] && cmp -s "$TMP/p" "$TMP/o"; then
        echo "PASS $name (rt)"; PASS=$((PASS+1))
    else
        echo "FAIL $name (rt)  (rc=$rc)"; sed 's/^/    /' "$TMP/rxe" | tail -3
        FAIL=$((FAIL+1))
    fi
}

# expect <name> <rx> args: compare $TMP/o to $TMP/x, check exit code
expect() {
    name=$1 rc=$2 want=$3
    if [ "$rc" -eq "$want" ] && cmp -s "$TMP/x" "$TMP/o"; then
        echo "PASS $name"; PASS=$((PASS+1))
    else
        echo "FAIL $name  (rc=$rc want=$want got=$(stat -c%s "$TMP/o")B)"
        FAIL=$((FAIL+1))
    fi
}

echo "== clean sizes (fast pipe) =="
for n in 0 1 64 1024 1025 3000 10000; do run "size$n" $n; done

echo "== paced real-time feeds =="
runrt rt3000 3000
runrt rt10000 10000

echo "== impairments =="
run "ppm+1000" 1500 -r 1000 -s 3
run "ppm-1000" 1500 -r -1000 -s 3
# +-100 ppm sits between the 300-step coarse points: only full coverage
# (winner-adjacent fine or the nothing-decoded escalation) finds it
run "ppm+100"  1500 -r 100 -s 3
run "ppm-100"  1500 -r -100 -s 3
run "cfo+295"  1500 -f 295 -s 5
run "cfo-295"  1500 -f -295 -s 5
run "noise.01" 1500 -n 0.01 -s 7
run "noise.02" 1500 -n 0.02 -s 7
run "gain-30"  1500 -g -30 -s 9
run "gain+6"   1500 -g 6 -s 9
run "delay7"   1500 -d 7 -s 2
run "fdelay.5" 1500 -D 0.5 -s 4
run "combo"    2500 -f 260 -d 777 -r 100 -n 0.01 -g -10 -s 11

echo "== PHY variants =="
PHY="-k 16";        run "k16"    1500
PHY="-k 40";        run "k40"    1500
PHY="-b 0.5";       run "bt0.5"  1500
PHY="--pulse rcos"; run "rcos"   1500
PHY="-d 4";         run "delay4" 1500
PHY=""

echo "== modulation index sweep (tone bank) =="
PHY="-H 0.5 -c 5000";  run "h0.5"    1500
PHY="-H 0.5 -c 5000";  run "h0.5ppm" 1500 -r 100 -s 3
PHY="-H 1.0 -c 7000";  run "h1.0"    1500
PHY="-H 1.0 -c 7000";  run "h1.0imp" 1500 -f 250 -n 0.005 -g -20 -s 5
PHY="-H 1.0 -k 40";    run "h1.0k40" 1500
PHY=""

echo "== mid-stream join =="
pay 1500 join1 > "$TMP/p"
$TX < "$TMP/p" > "$TMP/a" 2>/dev/null
# start inside frame 1's preamble: frame 1 lost, frame 2 (476 B) must arrive
tail -c +1001 "$TMP/a" | timeout 120 $RX > "$TMP/o" 2>/dev/null
rc=$?; tail -c 476 "$TMP/p" > "$TMP/x"; expect join-preamble $rc 0
# start inside frame 1's payload: frame 1 lost, nothing complete -> exit 1
head -c 110000 "$TMP/a" | timeout 120 $RX > "$TMP/o" 2>/dev/null
rc=$?; : > "$TMP/x"; expect join-midframe $rc 1

echo "== idle gap between streams =="
pay 700 g1 > "$TMP/p1"; pay 900 g2 > "$TMP/p2"
$TX < "$TMP/p1" > "$TMP/a1" 2>/dev/null
$TX < "$TMP/p2" > "$TMP/a2" 2>/dev/null
head -c 192000 /dev/zero > "$TMP/z"           # 2 s of digital silence
cat "$TMP/p1" "$TMP/p2" > "$TMP/x"
cat "$TMP/a1" "$TMP/z" "$TMP/a2" | timeout 120 $RX > "$TMP/o" 2>/dev/null
expect idle-gap $? 0

echo "== EOF mid-burst =="
head -c 5000 "$TMP/a" > "$TMP/t"              # a = 2-frame stream above
timeout 60 $RX < "$TMP/t" > "$TMP/o" 2>/dev/null
rc=$?; : > "$TMP/x"; expect eof-midframe $rc 1

echo "== --max-frames =="
pay 3000 mf > "$TMP/p"
$TX < "$TMP/p" > "$TMP/a" 2>/dev/null
timeout 60 $RX --max-frames 1 < "$TMP/a" > "$TMP/o" 2>/dev/null
rc=$?; head -c 1024 "$TMP/p" > "$TMP/x"; expect max-frames $rc 0

echo
echo "PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ] || exit 1
