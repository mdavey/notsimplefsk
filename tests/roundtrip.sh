#!/bin/sh
# notsimplefsk end-to-end regression battery
# run from the repo root:  sh tests/roundtrip.sh   (after `make`)
#
# pass case  : bytes in == bytes out (exit 0 per case, counted)
# reject case: PHY pair is physically impossible; we require only that the
#              RX never mis-delivers the payload (CRC gate must hold).
set -u
cd "$(dirname "$0")/.." || exit 1
TX=./bin/fsk_tx
RX=./bin/fsk_rx
CH=./bin/channel
TMP=${TMPDIR:-/tmp}/nsfsk-rt.$$
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

PASS=0; FAIL=0
PHY=""                    # extra PHY flags, must match tx+rx

# pay <nbytes> <seedname> : deterministic pseudo-random payload
pay() {
    LC_ALL=C gawk -v n="$1" -v s="$(printf %s "$2" | cksum | cut -d' ' -f1)" '
        BEGIN{srand(s);for(i=0;i<n;i++)printf "%c",int(rand()*256)}'
}

# run <name> <nbytes> [channel args...]
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
    $RX $PHY < "$TMP/a" > "$TMP/o" 2>"$TMP/rxe"
    if cmp -s "$TMP/p" "$TMP/o"; then
        echo "PASS $name"; PASS=$((PASS+1))
    else
        echo "FAIL $name"; sed 's/^/    /' "$TMP/rxe" | tail -3
        FAIL=$((FAIL+1))
    fi
}

# rej <name> [tx/rx flags...] : must NOT deliver the payload
rej() {
    name=$1; shift
    PHY="$*"
    pay 64 "$name" > "$TMP/p"
    $TX $PHY < "$TMP/p" > "$TMP/a" 2>/dev/null
    if $RX $PHY < "$TMP/a" > "$TMP/o" 2>/dev/null && cmp -s "$TMP/p" "$TMP/o"; then
        echo "FAIL $name  (should not decode)"; FAIL=$((FAIL+1))
    else
        echo "REJECT $name (expected)"; PASS=$((PASS+1))
    fi
    PHY=""
}

echo "== clean sizes (single frames) =="
for n in 0 1 3 17 255 1023 1024; do run "size$n" $n; done

echo "== multi-frame streams =="
for n in 1025 2048 3000 5000; do run "multi$n" $n; done

echo "== PHY sweep =="
PHY="-k 16";              run "k16"    300
PHY="-k 24";              run "k24"    300
PHY="-k 40";              run "k40"    300
PHY="-k 64";              run "k64"    300
PHY="-b 0.5";             run "bt0.5"  300
PHY="-b 1.0";             run "bt1.0"  300
PHY="-d 2";               run "delay2" 300
PHY="-d 4";               run "delay4" 300
PHY="-d 8";               run "delay8" 300
PHY="--pulse rcos";       run "rcos"   300
PHY="--pulse square";     run "square" 300
PHY="-c 8000";            run "fc8k"   300
PHY="-H 0.3";             run "h0.3"   300
PHY="-P 120";             run "preamble120" 300
PHY="-r 44100 -c 3800";   run "rate44100" 300
PHY="-r 22050 -c 3500";   run "rate22050" 300
PHY="-r 16000 -c 3200";   run "rate16000" 300
PHY="-r 8000  -c 2600";   run "rate8000"  300
PHY="-k 16 -d 3 -b 1.0";  run "k16mix" 500
PHY=""

echo "== impairments (single) =="
run "cfo+100"  300 -f  100
run "cfo-100"  300 -f -100
run "cfo+295"  300 -f  295
run "cfo-295"  300 -f -295
run "ppm+300"  300 -r  300
run "ppm-300"  300 -r -300
run "ppm+1000" 300 -r  1000
run "ppm-1000" 300 -r -1000
run "delay100" 300 -d 100
run "delay1000" 300 -d 1000
run "fdelay.5" 300 -D 0.5
run "gain-30"  300 -g -30
run "gain+6"   300 -g 6
run "noise.01" 300 -n 0.01
run "noise.03" 300 -n 0.03

echo "== impairment combinations =="
PHY=""; run "multi+cfo300"  2048 -f 300
run "multi+ppm1000" 2048 -r 1000
run "cfo200+noise"  300 -f 200 -n 0.01
run "fdelay+cfo"    300 -D 0.5 -f 150
run "ppm600+delay"  300 -r 600 -d 250
run "bigcombo"      300 -f 260 -d 777 -r 100 -n 0.01 -g -10
PHY="-k 24"; run "k24cfo200" 300 -f 200
PHY="-r 44100 -c 3800"; run "r44+cfo300" 300 -f 300
PHY=""; run "multi5k+cfo" 5000 -f 120
run "evilcombo" 5000 -f 295 -D 0.5 -r 600 -n 0.01

echo "== expected rejects (impossible fc/image combos) =="
rej "k8@fc4000"     -k 8
rej "k6@fc4000"     -k 6
rej "k4@fc4000"     -k 4
rej "h0.4"          -H 0.4
rej "fc1000@22050"  -r 22050 -c 1000

echo
echo "PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
