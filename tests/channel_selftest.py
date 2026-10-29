#!/usr/bin/env python3
"""Fast property/self-tests for bin/channel.

Run from the repo root (make test does this): python3 tests/channel_selftest.py
Only the Python standard library is used.  The tests are deliberately about
channel properties, not decoder performance; the stress matrix covers that.
"""
import math
import os
import re
import struct
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHANNEL = os.environ.get("CHANNEL_BIN", os.path.join(ROOT, "bin", "channel"))
RATE = 48000
FAILS = []


def run(args, samples, timeout=30):
    data = struct.pack("<%dh" % len(samples), *samples)
    p = subprocess.run([CHANNEL] + args, input=data, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, timeout=timeout)
    if p.returncode != 0:
        raise RuntimeError("channel %s failed: %s" %
                           (" ".join(args), p.stderr.decode(errors="replace")))
    out = list(struct.unpack("<%dh" % (len(p.stdout) // 2), p.stdout))
    return out, p.stderr.decode(errors="replace")


def sine(n, freq, amp, phase=0.0):
    return [int(round(amp * 32767.0 *
                      math.sin(2.0 * math.pi * freq * i / RATE + phase)))
            for i in range(n)]


def rms(x):
    return math.sqrt(sum(float(v) * v for v in x) / len(x)) if x else 0.0


def db(x):
    return 20.0 * math.log10(max(x, 1e-15))


def expect(name, cond, detail=""):
    if cond:
        print("PASS %s" % name)
    else:
        print("FAIL %s%s" % (name, ("  " + detail) if detail else ""))
        FAILS.append(name)


def report_fields(err):
    m = re.search(r"channel-report: (.*)", err)
    if not m:
        return {}
    out = {}
    for part in m.group(1).split():
        if "=" in part:
            k, v = part.split("=", 1)
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def main():
    if not os.path.exists(CHANNEL):
        print("SKIP: %s not built" % CHANNEL)
        return 0

    # 1. determinism of all stochastic components
    x = sine(9600, 1000, 0.2)
    args = ["--fade-doppler", "3", "--fade-rician", "6", "--snr-db", "20",
            "--clicks", "3:-20:5", "--dropouts", "1:10:zero", "-s", "7"]
    a, _ = run(args, x)
    b, _ = run(args, x)
    expect("deterministic stochastic channel", a == b)

    # 2. integer + fractional delay
    imp = [0] * 4800
    imp[100] = 16384
    y, _ = run(["-d", "10", "-D", "0.5"], imp)
    expect("delay length", len(y) == 4810, "len=%d" % len(y))
    pk = sorted(range(len(y)), key=lambda i: -abs(y[i]))[:3]
    expect("fractional delay position",
           set(pk[:2]) == {110, 111} and abs(y[110] - 8192) <= 2 and
           abs(y[111] - 8192) <= 2,
           "peaks=%s" % [(i, y[i]) for i in pk])

    # 3. static multipath impulse response
    imp0 = [0] * 4800
    imp0[0] = 16384
    y, _ = run(["--multipath", "0:0,1:-6"], imp0)
    expect("multipath length", len(y) == 4848, "len=%d" % len(y))
    expect("multipath taps",
           abs(y[0] - 16384) <= 2 and
           abs(y[48] - 16384 * 10 ** (-6.0 / 20.0)) <= 32,
           "y0=%d y48=%d" % (y[0], y[48]))

    # 4. SNR calibration (signal + noise only, no analytic path)
    n = 48000
    x = sine(n, 1000, 0.25)
    y, err = run(["--snr-db", "10", "-s", "7", "--report"], x)
    sig = rms(x) / 32768.0
    noise = rms([a - b for a, b in zip(y, x)]) / 32768.0
    got = db(sig / max(noise, 1e-15))
    expect("snr-db calibration", abs(got - 10.0) <= 1.5,
           "got %.2f dB" % got)
    rf = report_fields(err)
    expect("report has sigma", "sigma" in rf and rf["sigma"] > 0.0)

    # 5. bandpass pass/stop
    x = sine(9600, 1000, 0.25)
    y, _ = run(["--bandpass", "200:4000"], x)
    gain_db = db(rms(y) / max(rms(x), 1.0))
    expect("bandpass passband", abs(gain_db) <= 1.0, "gain %.2f dB" % gain_db)
    x = sine(9600, 8000, 0.25)
    y, _ = run(["--bandpass", "200:4000"], x)
    att_db = db(rms(y) / max(rms(x), 1.0))
    expect("bandpass stopband", att_db <= -20.0, "gain %.2f dB" % att_db)

    # 6. fading: finite, roughly unit power, deterministic, varies
    x = sine(48000, 1000, 0.25)
    y, err = run(["--fade-doppler", "2", "--report"], x)
    ratio = rms(y) / max(rms(x), 1.0)
    rf = report_fields(err)
    expect("fading finite", all(abs(v) <= 32768 for v in y))
    expect("fading mean power", 0.2 <= ratio <= 3.0, "ratio %.2f" % ratio)
    expect("fading varies", rf.get("fade_min_db", 0) < rf.get("fade_max_db", 0),
           "fade=%s..%s" % (rf.get("fade_min_db"), rf.get("fade_max_db")))

    # 7. AGC pushes a level step toward the target
    n = 48000
    x = [int(round(0.02 * 32767.0 * math.sin(2 * math.pi * 1000 * i / RATE)))
         if i < n // 2 else
         int(round(0.6 * 32767.0 * math.sin(2 * math.pi * 1000 * i / RATE)))
         for i in range(n)]
    y, _ = run(["--agc", "-18:5:500"], x)
    quiet_in = rms(x[:n // 2])
    quiet_out = rms(y[:n // 2])
    loud_in = rms(x[n // 2:])
    loud_out = rms(y[n // 2:])
    expect("agc raises quiet", quiet_out > quiet_in * 1.5,
           "in=%.0f out=%.0f" % (quiet_in, quiet_out))
    expect("agc lowers loud", loud_out < loud_in * 0.9,
           "in=%.0f out=%.0f" % (loud_in, loud_out))
    expect("agc no overflow", max(abs(v) for v in y) <= 32768)

    # 8. quantisation to 8 bits -> multiples of 256 in s16
    x = sine(9600, 1000, 0.25)
    y, _ = run(["--quantize", "8"], x)
    bad = [v for v in y if abs(v) % 256 != 0]
    expect("quantize 8", not bad, "bad=%d" % len(bad))

    # 8b. optional RMS normalisation restores the input level
    x = sine(9600, 1000, 0.10)
    y, _ = run(["-g", "-6", "--normalize"], x)
    ratio = rms(y) / max(rms(x), 1.0)
    expect("normalize", 0.85 <= ratio <= 1.15, "ratio %.2f" % ratio)

    # 9. deterministic dropout (zero) in place
    x = sine(9600, 1000, 0.25)
    y, _ = run(["--drop-at", "10:5:zero"], x)
    expect("dropout length zero", len(y) == len(x), "len=%d" % len(y))
    expect("dropout zeros", all(v == 0 for v in y[480:720]) and
           any(v != 0 for v in y[720:1200]))
    y2, _ = run(["--drop-at", "10:5:drop"], x)
    expect("dropout remove", len(y2) == len(x) - 240, "len=%d" % len(y2))

    # 10. old short flags still produce finite, sane output
    x = sine(9600, 1000, 0.25)
    for args in (["-d", "10"], ["-D", "0.5"], ["-r", "300"],
                 ["-g", "-6"], ["-f", "200"], ["-n", "0.01"], ["-s", "3"]):
        y, _ = run(args, x)
        ok = bool(y) and all(abs(v) <= 32768 for v in y)
        if args[0] == "-f":
            ok = ok and rms(y) > 0.05 * 32768
        expect("old flag %s" % " ".join(args), ok, "len=%d" % len(y))

    print()
    if FAILS:
        print("CHANNEL SELFTEST FAIL (%d)" % len(FAILS))
        return 1
    print("CHANNEL SELFTEST PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
