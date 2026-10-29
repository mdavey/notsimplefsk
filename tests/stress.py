#!/usr/bin/env python3
"""Realistic impairment/PHY stress matrix for notsimplefsk.

Runs TX -> channel impairments -> RX (offline and streaming) over a curated
matrix, checks frame-sequence integrity, and compares stochastic cases against
a checked-in per-seed baseline.

Usage:
  python3 tests/stress.py --quick
  python3 tests/stress.py --full --update-baseline
  python3 tests/stress.py --filter 'p0-default/rayleigh' --keep-failures /tmp/f
  python3 tests/stress.py --list

No third-party dependencies.  See docs/reference.md for the channel options
and the baseline workflow.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
TX = os.path.join(BIN, "fsk_tx")
RX_OFF = os.path.join(BIN, "fsk_rx")
RX_STREAM = os.path.join(BIN, "fsk_rx_stream")
CHANNEL = os.path.join(BIN, "channel")
BASELINE = os.path.join(ROOT, "tests", "stress_baseline.json")
FIXTURES = os.path.join(ROOT, "tests", "fixtures")

PHYS = {
    "p0-default": ["-r", "48000", "-c", "4000", "-k", "20", "-H", "0.25",
                   "-b", "0.8", "-d", "1", "-P", "72"],
    "p1-highrate": ["-r", "48000", "-c", "8000", "-k", "16", "-H", "0.25",
                    "-b", "0.8", "-d", "1", "-P", "72"],
    "p2-robust": ["-r", "48000", "-c", "5000", "-k", "40", "-H", "0.5",
                  "-b", "1.0", "-d", "4", "-P", "72"],
    "p3-tonebank": ["-r", "48000", "-c", "7000", "-k", "20", "-H", "1.0",
                    "-b", "0.8", "-d", "1", "-P", "72"],
    "p4-acoustic": ["-r", "48000", "-c", "1000", "-k", "160", "-H", "1.0",
                    "-b", "1.0", "-d", "4", "-P", "72"],
    "p5-rcos": ["-r", "48000", "-c", "4000", "-k", "20", "-H", "0.25",
                "-b", "0.8", "-d", "1", "-P", "72", "--pulse", "rcos"],
    "p6-square": ["-r", "48000", "-c", "4000", "-k", "20", "-H", "0.25",
                  "-b", "0.8", "-d", "1", "-P", "72", "--pulse", "square"],
    "p7-lowrate": ["-r", "8000", "-c", "2600", "-k", "20", "-H", "0.25",
                   "-b", "0.8", "-d", "1", "-P", "72"],
    "p8-pre120": ["-r", "48000", "-c", "4000", "-k", "20", "-H", "0.25",
                  "-b", "0.8", "-d", "1", "-P", "120"],
    "p9-bt065": ["-r", "48000", "-c", "4000", "-k", "20", "-H", "0.25",
                 "-b", "0.65", "-d", "1", "-P", "72"],
}

# Each scenario is a dict:
#   channel   channel args (seed is added by the runner)
#   payload   bytes
#   tx_gain   TX -g (TX-only; never passed to RX)
#   gated     True => every frame must be delivered
#   quick     included in --quick
#   phys      optional whitelist
#   tolerance allowed baseline shortfall (default 0)
#   allow_clip  do not fail when the channel clips
#   split     optional first-part size for an inserted idle gap
#   gap_s     inserted silence in seconds (split scenarios)
SCENARIOS = {
    # ---- deterministic / gated ----
    "clean": dict(channel=[], payload=3000, tx_gain=-12, gated=True,
                  quick=True),
    "awgn25": dict(channel=["--snr-db", "25"], payload=3000, tx_gain=-12,
                   gated=True, quick=True),
    "awgn20": dict(channel=["--snr-db", "20"], payload=3000, tx_gain=-12,
                   gated=True, quick=True),
    "cfo+295": dict(channel=["-f", "295"], payload=3000, tx_gain=-12,
                    gated=True, quick=True,
                    phys=["p0-default", "p1-highrate", "p3-tonebank",
                          "p5-rcos", "p6-square", "p8-pre120", "p9-bt065"]),
    "cfo-295": dict(channel=["-f", "-295"], payload=3000, tx_gain=-12,
                    gated=True, quick=True,
                    phys=["p0-default", "p1-highrate", "p3-tonebank",
                          "p5-rcos", "p6-square", "p8-pre120", "p9-bt065"]),
    "cfo+30": dict(channel=["-f", "30"], payload=3000, tx_gain=-12,
                   gated=True, quick=True,
                   phys=["p2-robust", "p4-acoustic", "p7-lowrate"]),
    "cfo-30": dict(channel=["-f", "-30"], payload=3000, tx_gain=-12,
                   gated=True, quick=True,
                   phys=["p2-robust", "p4-acoustic", "p7-lowrate"]),
    "ppm+300": dict(channel=["-r", "300"], payload=3000, tx_gain=-12,
                    gated=True, quick=True),
    "ppm-300": dict(channel=["-r", "-300"], payload=3000, tx_gain=-12,
                    gated=True, quick=True),
    "ppm+1000": dict(channel=["-r", "1000"], payload=3000, tx_gain=-12,
                     gated=False, quick=True, tolerance=0),
    "ppm-1000": dict(channel=["-r", "-1000"], payload=3000, tx_gain=-12,
                     gated=False, quick=True, tolerance=0),
    "ppm-ramp": dict(channel=["--ppm-ramp", "-300:300"], payload=3000,
                     tx_gain=-12, gated=False, quick=True, tolerance=0),
    "echo-0.1ms": dict(channel=["--multipath", "0:0,0.1:-6"], payload=3000,
                       tx_gain=-20, gated=True, quick=True,
                       known_hard_phys=["p2-robust"]),
    "echo-1ms-20": dict(channel=["--multipath", "0:0,1:-20"], payload=3000,
                        tx_gain=-20, gated=True, quick=True),
    "bp-dc-hum": dict(channel=["--bandpass", "200:8000", "--dc", "0.02",
                               "--hum", "50:-30"],
                      payload=3000, tx_gain=-12, gated=True, quick=True),
    "quant8-soft": dict(channel=["--quantize", "8", "--soft-clip", "2"],
                        payload=3000, tx_gain=-12, gated=True, quick=True),
    "agc-level": dict(channel=["--agc", "-18:5:500", "--gain-steps",
                               "0:0,3000:-12,6000:0"],
                      payload=3000, tx_gain=-12, gated=False, quick=True,
                      tolerance=0),
    "long-drift": dict(channel=["--ppm-ramp", "-300:300"], payload=32768,
                       tx_gain=-12, gated=True, quick=False,
                       phys=["p0-default", "p3-tonebank", "p5-rcos"],
                       timeout=300),
    "idle-gap": dict(channel=[], payload=8192, tx_gain=-12, gated=True,
                     quick=False, split=4096, gap_s=2.0,
                     phys=["p0-default", "p3-tonebank", "p5-rcos"],
                     timeout=300),

    # ---- deterministic but outside the current RX envelope (XFAIL) ----
    "echo-1ms-6": dict(channel=["--multipath", "0:0,1:-6"], payload=3000,
                       tx_gain=-20, gated=False, quick=True, tolerance=0),
    "strong-tone": dict(channel=["--tone", "4300:-20", "--snr-db", "30"],
                        payload=3000, tx_gain=-20, gated=False, quick=True,
                        tolerance=0, skip_phys=["p2-robust"]),
    "clip+12": dict(channel=["-g", "12"], payload=3000, tx_gain=-12,
                    gated=False, quick=False, allow_clip=True, tolerance=0),
    "agc-noise": dict(channel=["--agc", "-18:5:500", "--snr-db", "20"],
                      payload=3000, tx_gain=-20, gated=False, quick=True,
                      tolerance=0),
    "reverb-heavy": dict(channel=["--reverb", "0.3:-6"], payload=8192,
                         tx_gain=-20, gated=False, quick=False, tolerance=1,
                         skip_phys=["p2-robust"]),

    # ---- stochastic / baseline-tracked ----
    "rician10-fd1": dict(channel=["--fade-doppler", "1", "--fade-rician", "10",
                                  "--snr-db", "30"],
                         payload=8192, tx_gain=-20, gated=False, quick=True,
                         tolerance=1),
    "rayleigh-fd1": dict(channel=["--fade-doppler", "1", "--snr-db", "30"],
                         payload=8192, tx_gain=-20, gated=False, quick=True,
                         tolerance=1),
    "rayleigh-fd5": dict(channel=["--fade-doppler", "5", "--snr-db", "20"],
                         payload=8192, tx_gain=-20, gated=False, quick=False,
                         tolerance=1),
    "sel-fade": dict(channel=["--fade-doppler", "2", "--fade-spread", "2",
                              "--fade-taps", "4", "--snr-db", "25"],
                     payload=8192, tx_gain=-20, gated=False, quick=True,
                     tolerance=1, skip_phys=["p2-robust"]),
    "reverb-mild": dict(channel=["--reverb", "0.3:-20", "--snr-db", "30"],
                        payload=8192, tx_gain=-20, gated=False, quick=True,
                        tolerance=1),
    "clicks": dict(channel=["--clicks", "2:-20:5", "--snr-db", "25"],
                   payload=8192, tx_gain=-20, gated=False, quick=True,
                   tolerance=1),
    "bursts": dict(channel=["--bursts", "1:20:-20", "--snr-db", "25"],
                   payload=8192, tx_gain=-20, gated=False, quick=False,
                   tolerance=1),
    "drop-mid": dict(channel=["--dropouts", "0.2:20:drop", "--snr-db", "30"],
                     payload=8192, tx_gain=-20, gated=False, quick=True,
                     tolerance=1),
    "low-snr12": dict(channel=["--snr-db", "12"], payload=8192, tx_gain=-12,
                      gated=False, quick=True, tolerance=1),
    # Known offline RX issue (not fixed in this step): at k=160 with noise,
    # fsk_rx can merge the same frame from two ppm trials and emit it twice.
    "p4-dup-awgn25": dict(channel=["--snr-db", "25"], payload=1536,
                          tx_gain=-12, gated=False, quick=True,
                          phys=["p4-acoustic"], tolerance=0,
                          known_hard=True, no_cap=True),
    "realistic-acoustic": dict(
        channel=["--bandpass", "200:4000", "--multipath", "1:-6,5:-12",
                 "--fade-doppler", "2", "--fade-rician", "12",
                 "--hum", "50:-30", "--snr-db", "25",
                 "--agc", "-18:5:500", "--quantize", "8",
                 "-r", "300", "-f", "50"],
        payload=2048, tx_gain=-24, gated=False, quick=True,
        phys=["p4-acoustic"], tolerance=1),
}

QUICK_PHYS = ["p0-default", "p3-tonebank", "p4-acoustic", "p5-rcos",
              "p6-square"]
SLOW_PHYS = {"p4-acoustic", "p7-lowrate"}


def phy_seed(seed: int) -> int:
    z = (seed + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return z ^ (z >> 31)


class Prng:
    def __init__(self, seed):
        self.s = phy_seed(seed) or 1

    def u64(self):
        x = self.s
        x ^= (x >> 12) & 0xFFFFFFFFFFFFFFFF
        x ^= (x << 25) & 0xFFFFFFFFFFFFFFFF
        x ^= x >> 27
        self.s = x
        return (x * 2685821657736338717) & 0xFFFFFFFFFFFFFFFF


def make_payload(n, key):
    digest = hashlib.sha256(key.encode()).digest()
    r = Prng(int.from_bytes(digest[:8], "little"))
    out = bytearray()
    while len(out) < n:
        out += struct.pack("<Q", r.u64())
    return bytes(out[:n])


def split_frames(payload):
    return [payload[i:i + 1024] for i in range(0, len(payload), 1024)]


def payload_size_for(phy, scenario):
    n = int(scenario["payload"])
    if scenario.get("no_cap"):
        return n
    if phy == "p4-acoustic":
        n = min(n, 1024)
    elif phy in SLOW_PHYS:
        n = min(n, 2048)
    if phy == "p2-robust" and not scenario["gated"]:
        n = min(n, 1024)
    return n


def phy_rate(phy):
    args = PHYS[phy]
    for i, a in enumerate(args):
        if a == "-r" and i + 1 < len(args):
            return int(args[i + 1])
    return 48000


def run_cmd(args, data, timeout):
    try:
        p = subprocess.run(args, input=data, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, timeout=timeout)
        return p.returncode, p.stdout, p.stderr, False
    except subprocess.TimeoutExpired as e:
        return None, e.stdout or b"", e.stderr or b"", True
    except OSError as e:
        return None, b"", str(e).encode(), False


def parse_report(err):
    m = re.search(rb"channel-report: (.*)", err)
    if not m:
        return {}
    out = {}
    for part in m.group(1).split():
        if b"=" in part:
            k, v = part.split(b"=", 1)
            try:
                out[k.decode()] = float(v)
            except ValueError:
                out[k.decode()] = v.decode(errors="replace")
    return out


def match_frames(out, frames):
    """Greedy in-order frame matching; skips are allowed, corruption is not."""
    oi = 0
    ei = 0
    delivered = []
    total = len(frames)
    while oi < len(out):
        if ei < total and out.startswith(frames[ei], oi):
            delivered.append(ei)
            oi += len(frames[ei])
            ei += 1
            continue
        found = None
        for j in range(ei + 1, total):
            if out.startswith(frames[j], oi):
                found = j
                break
        if found is None:
            return dict(hard=True,
                        reason="unmatched output at byte %d (of %d)" %
                               (oi, len(out)),
                        delivered=delivered, lost=list(range(total)),
                        total=total, prefix=0)
        delivered.append(found)
        oi += len(frames[found])
        ei = found + 1
    lost = [i for i in range(total) if i not in set(delivered)]
    prefix = 0
    seen = set(delivered)
    while prefix in seen:
        prefix += 1
    return dict(hard=False, delivered=delivered, lost=lost, total=total,
                prefix=prefix)


def is_gated(scenario, phy):
    return bool(scenario["gated"] and
                phy not in scenario.get("gate_exclude", []) and
                phy not in scenario.get("known_hard_phys", []))


def status_for(scenario, phy, key, base_cases, delivered, total, hard):
    if hard:
        return "FAIL"
    if is_gated(scenario, phy):
        return "PASS" if delivered == total else "FAIL"
    if key not in base_cases:
        return "NEW"
    exp = int(base_cases[key].get("frames", 0))
    tol = int(base_cases[key].get("tolerance", scenario.get("tolerance", 0)))
    if delivered < exp - tol:
        return "FAIL"
    if delivered < total:
        return "XFAIL"
    return "PASS"


def make_case_id(phy, scenario, mode, seed, payload):
    return "%s/%s/%s/s%d/b%d" % (phy, scenario, mode, seed, payload)


def task_id(phy, scenario, seed):
    return "%s/%s/s%d" % (phy, scenario, seed)


def run_task(task, ctx):
    phy_id, scenario_id, seed = task
    phy = PHYS[phy_id]
    sc = SCENARIOS[scenario_id]
    nbytes = payload_size_for(phy_id, sc)
    payload = make_payload(nbytes,
                           "%s/%s/%d/%d" % (phy_id, scenario_id, seed, nbytes))
    timeout = sc.get("timeout", 180)
    artifacts = {}
    commands = {}

    parts = [payload]
    if sc.get("split"):
        cut = int(sc["split"])
        parts = [payload[:cut], payload[cut:]]
    expected_frames = []
    audio_parts = []
    for pi, part in enumerate(parts):
        tx_args = [TX] + phy + ["-g", str(sc["tx_gain"])]
        rc, out, err, tmo = run_cmd(tx_args, part, 120)
        commands["tx%d" % pi] = " ".join(tx_args)
        if tmo or rc != 0:
            artifacts["payload"] = payload
            return dict(task=task, error="tx failed/rc=%s%s" %
                        (rc, " timeout" if tmo else ""),
                        artifacts=artifacts, commands=commands)
        audio_parts.append(out)
        expected_frames.extend(split_frames(part))

    if len(audio_parts) == 1:
        tx_audio = audio_parts[0]
    else:
        gap = int(sc.get("gap_s", 0.0) * phy_rate(phy_id))
        tx_audio = audio_parts[0] + b"\x00" * (2 * gap) + audio_parts[1]

    ch_args = [CHANNEL] + sc["channel"] + ["-s", str(seed), "--report"]
    commands["channel"] = " ".join(ch_args)
    rc, ch_out, ch_err, tmo = run_cmd(ch_args, tx_audio, timeout)
    if tmo or rc != 0:
        artifacts.update(payload=payload, tx=tx_audio, channel_stderr=ch_err)
        return dict(task=task, error="channel failed/rc=%s%s" %
                    (rc, " timeout" if tmo else ""),
                    artifacts=artifacts, commands=commands)
    rep = parse_report(ch_err)
    clips = int(rep.get("clips", 0))
    if clips > 0 and not sc.get("allow_clip", False):
        artifacts.update(payload=payload, tx=tx_audio, channel=ch_out,
                         channel_stderr=ch_err)
        return dict(task=task, error="unintended clipping (%d samples)" % clips,
                    artifacts=artifacts, commands=commands, report=rep)

    results = {}
    any_short = False
    rx_art = {}
    for mode in ctx["modes"]:
        rx = RX_OFF if mode == "off" else RX_STREAM
        rx_args = [rx] + phy
        commands["rx_" + mode] = " ".join(rx_args)
        rc, out, err, tmo = run_cmd(rx_args, ch_out, timeout)
        m = match_frames(out, expected_frames)
        hard = False
        reason = ""
        if tmo:
            hard, reason = True, "rx timeout"
        elif m["hard"]:
            hard, reason = True, m["reason"]
        else:
            want_rc = 0 if len(m["delivered"]) > 0 else 1
            if rc != want_rc:
                hard = True
                reason = "rx rc=%s want %d" % (rc, want_rc)
        if len(m["delivered"]) < len(expected_frames):
            any_short = True
        results[mode] = dict(
            rc=rc, hard=hard, reason=reason,
            delivered=len(m["delivered"]), total=len(expected_frames),
            lost=m["lost"], prefix=m["prefix"],
            stdout_len=len(out), stdout=out if (hard or ctx["keep_failures"]) else None,
            stderr=err.decode(errors="replace")[-4000:],
            mode=mode,
        )
        if hard or ctx["keep_failures"]:
            rx_art["rx_%s_stdout" % mode] = out
            rx_art["rx_%s_stderr" % mode] = err
    hard_any = any(r["hard"] for r in results.values())
    keep = (any_short and ctx["keep_failures"]) or hard_any or (
        any_short and is_gated(sc, phy_id))
    if keep:
        artifacts.update(payload=payload, tx=tx_audio, channel=ch_out,
                         channel_stderr=ch_err)
        artifacts.update(rx_art)
        if ctx.get("failure_root") and ctx["keep_failures"]:
            save_artifacts(artifacts, commands,
                           {"key": task_id(phy_id, scenario_id, seed)},
                           ctx["failure_root"])
            artifacts = {}
    return dict(task=task, results=results, artifacts=artifacts,
                commands=commands, report=rep)


def case_key(phy, scenario, mode, seed, nbytes):
    return make_case_id(phy, scenario, mode, seed, nbytes)


def build_tasks(args, modes):
    phys = list(PHYS)
    if args.phy:
        rx = re.compile(args.phy)
        phys = [p for p in phys if rx.search(p)]
    if args.quick:
        phys = [p for p in phys if p in QUICK_PHYS]
    scenarios = []
    for name, sc in SCENARIOS.items():
        if args.quick and not sc.get("quick", False):
            continue
        if args.scenario and not re.search(args.scenario, name):
            continue
        if sc.get("phys") and not any(p in sc["phys"] for p in phys):
            continue
        scenarios.append(name)
    tasks = []
    seeds = args.seeds if args.seeds else ([1, 2] if args.quick else [1, 2, 3, 4])
    frx = re.compile(args.filter) if args.filter else None
    for phy in phys:
        for sid in scenarios:
            sc = SCENARIOS[sid]
            if sc.get("phys") and phy not in sc["phys"]:
                continue
            if phy in sc.get("skip_phys", []):
                continue
            for seed in seeds:
                if frx is not None:
                    nb = payload_size_for(phy, sc)
                    ids = [case_key(phy, sid, m, seed, nb) for m in modes]
                    if not any(frx.search(x) for x in ids):
                        continue
                tasks.append((phy, sid, seed))
    return tasks


def print_list(tasks):
    for t in tasks:
        print(task_id(*t))
    print("tasks=%d" % len(tasks))


def matrix_hash():
    blob = json.dumps({"phys": PHYS, "scenarios": SCENARIOS},
                      sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def load_baseline(path):
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_baseline(path, cases):
    doc = {
        "schema": 1,
        "suite": "stress-v1",
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "git_commit": git_commit(),
        "matrix_hash": matrix_hash(),
        "cases": {k: cases[k] for k in sorted(cases)},
    }
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=1, sort_keys=True)
        f.write("\n")
    os.replace(tmp, path)


def git_commit():
    try:
        p = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"],
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        return p.stdout.decode().strip()
    except OSError:
        return ""


def save_artifacts(artifacts, commands, result, root):
    if not artifacts:
        return None
    key = result.get("key", "case").replace("/", "_")
    d = os.path.join(root, key)
    os.makedirs(d, exist_ok=True)
    for name, data in artifacts.items():
        if isinstance(data, bytes):
            with open(os.path.join(d, name + ".bin"), "wb") as f:
                f.write(data)
            continue
        if isinstance(data, str):
            with open(os.path.join(d, name + ".txt"), "w", encoding="utf-8") as f:
                f.write(data)
    safe = {k: v for k, v in result.items()
            if k not in ("artifacts", "stdout")}
    with open(os.path.join(d, "case.json"), "w", encoding="utf-8") as f:
        json.dump({"commands": commands, "result": safe}, f, indent=1,
                  default=str)
    return d


def classify_results(raw, ctx):
    results = []
    base_cases = ctx["base_cases"]
    for item in raw:
        if "error" in item:
            phy, scenario, seed = item["task"]
            sc = SCENARIOS[scenario]
            nbytes = payload_size_for(phy, sc)
            for mode in ctx["modes"]:
                key = case_key(phy, scenario, mode, seed, nbytes)
                results.append(dict(
                    key=key, phy=phy, scenario=scenario, mode=mode, seed=seed,
                    status="FAIL", reason=item["error"], hard=True,
                    delivered=0, total=len(split_frames(make_payload(nbytes,
                        "%s/%s/%d/%d" % (phy, scenario, seed, nbytes)))),
                    commands=item.get("commands", {}),
                    artifacts=item.get("artifacts", {})))
            continue
        phy, scenario, seed = item["task"]
        sc = SCENARIOS[scenario]
        nbytes = payload_size_for(phy, sc)
        for mode, r in item["results"].items():
            key = case_key(phy, scenario, mode, seed, nbytes)
            base = base_cases.get(key)
            tol = sc.get("tolerance", 0)
            if base is not None:
                tol = int(base.get("tolerance", tol))
            known_hard = (sc.get("known_hard", False) or
                          phy in sc.get("known_hard_phys", []))
            if r["hard"] and known_hard:
                worse = False
                if base is not None:
                    bexp = int(base.get("frames", 0))
                    btol = int(base.get("tolerance", tol))
                    worse = r["delivered"] < bexp - btol
                status = "FAIL" if worse else "XFAIL"
                reason = (("known hard failure regressed: " if worse else
                           "known hard failure: ") + r["reason"])
            else:
                status = status_for(sc, phy, key, base_cases, r["delivered"],
                                    r["total"], r["hard"])
                reason = r["reason"]
            if status == "FAIL" and not r["hard"]:
                if is_gated(sc, phy):
                    reason = "gated: %d/%d frames delivered" % (
                        r["delivered"], r["total"])
                else:
                    exp = int(base.get("frames", 0)) if base else 0
                    reason = "baseline regression: %d/%d (expected >= %d)" % (
                        r["delivered"], r["total"], exp - tol)
            if status == "NEW":
                reason = "no baseline entry"
            results.append(dict(
                key=key, phy=phy, scenario=scenario, mode=mode, seed=seed,
                status=status, reason=reason, hard=r["hard"],
                delivered=r["delivered"], total=r["total"],
                lost=r["lost"], prefix=r["prefix"], rc=r["rc"],
                stdout_len=r["stdout_len"], stderr=r["stderr"],
                commands=item.get("commands", {}),
                artifacts=item.get("artifacts", {}),
                report=item.get("report", {})))
    return results



PACED_CASES = [("p0-default", "clean"), ("p0-default", "awgn20"),
               ("p0-default", "ppm+300"), ("p0-default", "drop-mid")]


def run_paced_case(phy_id, scenario_id, seed):
    """Feed a stream receiver at ~5x real time through a pipe."""
    phy = PHYS[phy_id]
    sc = SCENARIOS[scenario_id]
    nbytes = payload_size_for(phy_id, sc)
    payload = make_payload(nbytes, "paced/%s/%s/%d/%d" %
                          (phy_id, scenario_id, seed, nbytes))
    tx_args = [TX] + phy + ["-g", str(sc["tx_gain"])]
    rc, tx_audio, tx_err, tmo = run_cmd(tx_args, payload, 120)
    if tmo or rc != 0:
        return dict(key="paced/%s/%s/stream/s%d" % (phy_id, scenario_id, seed),
                    phy=phy_id, scenario="paced:" + scenario_id,
                    mode="stream", seed=seed, status="FAIL",
                    reason="tx failed", hard=True, delivered=0,
                    total=len(split_frames(payload)), commands={},
                    artifacts={})
    ch_args = [CHANNEL] + sc["channel"] + ["-s", str(seed), "--report"]
    rc, ch_out, ch_err, tmo = run_cmd(ch_args, tx_audio, sc.get("timeout", 180))
    if tmo or rc != 0:
        return dict(key="paced/%s/%s/stream/s%d" % (phy_id, scenario_id, seed),
                    phy=phy_id, scenario="paced:" + scenario_id,
                    mode="stream", seed=seed, status="FAIL",
                    reason="channel failed", hard=True, delivered=0,
                    total=len(split_frames(payload)),
                    commands={"channel": " ".join(ch_args)},
                    artifacts={"payload": payload, "tx": tx_audio,
                               "channel_stderr": ch_err})
    frames = split_frames(payload)
    tmpd = tempfile.mkdtemp(prefix="nsfsk-paced-")
    out_f = open(os.path.join(tmpd, "out.bin"), "wb")
    err_f = open(os.path.join(tmpd, "err.txt"), "wb")
    rx_args = [RX_STREAM] + phy
    p = subprocess.Popen(rx_args, stdin=subprocess.PIPE, stdout=out_f,
                         stderr=err_f)
    rate = phy_rate(phy_id)
    chunk = max(2, int(rate * 2 * 0.1))       # 0.1 s of s16le
    try:
        for off in range(0, len(ch_out), chunk):
            p.stdin.write(ch_out[off:off + chunk])
            p.stdin.flush()
            time.sleep(0.02)                   # ~5x real time
        p.stdin.close()
        p.wait(timeout=len(ch_out) / (2 * rate) / 5 + 30)
        rc = p.returncode
    except (BrokenPipeError, subprocess.TimeoutExpired):
        p.kill()
        rc = None
        p.wait()
    out_f.close()
    err_f.close()
    with open(os.path.join(tmpd, "out.bin"), "rb") as f:
        out = f.read()
    with open(os.path.join(tmpd, "err.txt"), "rb") as f:
        err = f.read()
    shutil.rmtree(tmpd, ignore_errors=True)
    m = match_frames(out, frames)
    hard = m["hard"] or rc != (0 if m["delivered"] else 1)
    delivered = len(m["delivered"])
    if hard:
        ok = False
    elif sc["gated"]:
        ok = delivered == len(frames)
    else:
        ok = delivered >= 1
    return dict(key="paced/%s/%s/stream/s%d" % (phy_id, scenario_id, seed),
                phy=phy_id, scenario="paced:" + scenario_id,
                mode="stream", seed=seed,
                status="PASS" if ok else "FAIL",
                reason=("" if ok else
                        (m.get("reason", "") or
                         "paced: %d/%d frames, rc=%s" %
                         (delivered, len(frames), rc))),
                hard=bool(hard), delivered=delivered, total=len(frames),
                commands={"rx": " ".join(rx_args)},
                artifacts=(dict(payload=payload, tx=tx_audio, channel=ch_out,
                                rx_stdout=out, rx_stderr=err)
                           if not ok else {}),
                stderr=err.decode(errors="replace")[-4000:], rc=rc,
                stdout_len=len(out))

def run_fixtures(args, ctx):
    runs = []
    fdir = args.fixtures
    if args.no_fixtures or not os.path.isdir(fdir):
        return runs
    for name in sorted(os.listdir(fdir)):
        cdir = os.path.join(fdir, name)
        cfgp = os.path.join(cdir, "case.json")
        if not os.path.isfile(cfgp):
            continue
        try:
            with open(cfgp, "r", encoding="utf-8") as f:
                cfg = json.load(f)
        except (OSError, ValueError) as e:
            runs.append(dict(key="fixture/%s" % name, phy="-", scenario=name,
                             mode="-", seed=0, status="FAIL",
                             reason="bad case.json: %s" % e, hard=True,
                             delivered=0, total=0, commands={},
                             artifacts={}))
            continue
        phy = cfg.get("phy", "")
        if isinstance(phy, str):
            phy = phy.split()
        audio_path = os.path.join(cdir, cfg.get("audio", "audio.raw"))
        payload_path = os.path.join(cdir, cfg.get("payload", "payload.bin"))
        try:
            with open(audio_path, "rb") as f:
                audio = f.read()
        except OSError as e:
            runs.append(dict(key="fixture/%s" % name, phy="-", scenario=name,
                             mode="-", seed=0, status="FAIL",
                             reason="missing audio: %s" % e, hard=True,
                             delivered=0, total=0, commands={},
                             artifacts={}))
            continue
        expect = cfg.get("expect", "all")
        payload = b""
        if expect != "none":
            try:
                with open(payload_path, "rb") as f:
                    payload = f.read()
            except OSError as e:
                runs.append(dict(key="fixture/%s" % name, phy="-",
                                 scenario=name, mode="-", seed=0,
                                 status="FAIL", reason="missing payload: %s" % e,
                                 hard=True, delivered=0, total=0,
                                 commands={}, artifacts={}))
                continue
        frames = split_frames(payload) if expect != "none" else []
        for mode in cfg.get("rx", ["offline", "stream"]):
            m = "off" if mode in ("off", "offline") else "stream"
            rx = RX_OFF if m == "off" else RX_STREAM
            rx_args = [rx] + phy
            rc, out, err, tmo = run_cmd(rx_args, audio,
                                        cfg.get("timeout_s", 120))
            match = match_frames(out, frames) if frames else dict(
                hard=False, delivered=[], lost=[], total=0, prefix=0)
            hard = tmo or match["hard"]
            reason = "timeout" if tmo else (match.get("reason", "") if match["hard"] else "")
            delivered = len(match["delivered"])
            if expect == "none":
                ok = (not tmo) and rc == 1 and len(out) == 0
            elif expect == "subsequence":
                ok = ((not hard) and delivered >= int(cfg.get("min_frames", 1))
                      and rc == 0)
            else:
                ok = ((not hard) and delivered == len(frames) and rc == 0)
            status = "PASS" if ok else "FAIL"
            if not ok and not reason:
                reason = "%s: %d/%d frames, rc=%s" % (expect, delivered,
                                                      len(frames), rc)
            runs.append(dict(key="fixture/%s/%s" % (name, m), phy="-",
                             scenario="fixture:" + name, mode=m, seed=0,
                             status=status, reason=reason, hard=bool(hard),
                             delivered=delivered, total=len(frames),
                             commands={"rx": " ".join(rx_args)},
                             artifacts={"audio": audio} if not ok else {},
                             stderr=err.decode(errors="replace")[-4000:],
                             rc=rc, stdout_len=len(out)))
    return runs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--quick", action="store_true", help="curated matrix")
    mode.add_argument("--full", action="store_true", help="full matrix")
    rx = ap.add_mutually_exclusive_group()
    rx.add_argument("--offline", action="store_true", help="only fsk_rx")
    rx.add_argument("--stream", action="store_true", help="only fsk_rx_stream")
    ap.add_argument("--both", action="store_true", help="both receivers (default)")
    ap.add_argument("--jobs", type=int, default=min(16, os.cpu_count() or 4))
    ap.add_argument("--seed", dest="seed_list", default=None,
                    help="comma-separated seeds")
    ap.add_argument("--scenario", default=None, help="scenario regex")
    ap.add_argument("--phy", default=None, help="PHY id regex")
    ap.add_argument("--filter", default=None, help="case id regex")
    ap.add_argument("--payload", type=int, default=None,
                    help="override payload bytes (for experiments)")
    ap.add_argument("--timeout", type=int, default=None,
                    help="override per-case timeout seconds")
    ap.add_argument("--baseline", default=BASELINE)
    ap.add_argument("--update-baseline", action="store_true")
    ap.add_argument("--strict", action="store_true",
                    help="treat XFAIL/NEW as failure")
    ap.add_argument("--keep-failures", default=None, metavar="DIR")
    ap.add_argument("--json", dest="json_out", default=None, metavar="FILE")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--fixtures", default=FIXTURES)
    ap.add_argument("--no-fixtures", action="store_true")
    ap.add_argument("--paced", action="store_true",
                    help="also run a few paced stream cases (full only)")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if args.full:
        args.quick = False
    else:
        args.quick = True
    if args.jobs < 1:
        ap.error("--jobs must be positive")
    if args.timeout is not None and args.timeout <= 0:
        ap.error("--timeout must be positive")
    if args.seed_list:
        try:
            args.seeds = [int(x) for x in args.seed_list.split(",") if x]
        except ValueError:
            ap.error("--seed must be comma-separated integers")
        if not args.seeds:
            ap.error("--seed must not be empty")
    else:
        args.seeds = None
    if args.payload is not None:
        if args.payload <= 0:
            ap.error("--payload must be positive")
        for sc in SCENARIOS.values():
            sc["payload"] = args.payload
    if args.timeout is not None:
        for sc in SCENARIOS.values():
            sc["timeout"] = args.timeout

    modes = ["off", "stream"]
    if args.offline:
        modes = ["off"]
    elif args.stream:
        modes = ["stream"]

    tasks = build_tasks(args, modes)
    if args.list:
        print_list(tasks)
        return 0

    for p in (TX, RX_OFF, RX_STREAM, CHANNEL):
        if not os.path.exists(p):
            print("missing %s (run make)" % p, file=sys.stderr)
            return 2

    baseline = load_baseline(args.baseline)
    base_cases = baseline.get("cases", {}) if baseline else {}
    if baseline and baseline.get("matrix_hash") != matrix_hash():
        print("warning: baseline matrix hash differs; run --full --update-baseline",
              file=sys.stderr)
    if baseline is None:
        print("warning: no baseline at %s; stochastic cases will be NEW" %
              args.baseline, file=sys.stderr)
    if args.update_baseline and not args.full:
        print("--update-baseline requires --full", file=sys.stderr)
        return 2
    if args.update_baseline and (args.filter or args.scenario or args.phy):
        print("--update-baseline requires the unfiltered full matrix",
              file=sys.stderr)
        return 2
    if args.update_baseline and (args.offline or args.stream):
        print("--update-baseline requires both receivers", file=sys.stderr)
        return 2
    if args.update_baseline and (args.payload is not None or args.paced):
        print("--update-baseline does not accept --payload/--paced",
              file=sys.stderr)
        return 2

    ctx = dict(modes=modes, keep_failures=bool(args.keep_failures),
               failure_root=args.keep_failures, base_cases=base_cases)
    raw = []
    t0 = time.time()
    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(run_task, t, ctx): t for t in tasks}
        for fut in concurrent.futures.as_completed(futs):
            try:
                raw.append(fut.result())
            except Exception as e:  # defensive: a worker bug must not hang the run
                t = futs[fut]
                raw.append(dict(task=t, error="worker exception: %s" % e,
                                artifacts={}, commands={}))
            done += 1
            if args.verbose:
                print("  [%d/%d] %s" % (done, len(tasks), task_id(*futs[fut])),
                      file=sys.stderr)
            elif done % 25 == 0:
                print("  %d/%d cases" % (done, len(tasks)), file=sys.stderr)
    results = classify_results(raw, ctx)
    if args.paced:
        if not args.full:
            print("warning: --paced requires --full; skipped", file=sys.stderr)
        else:
            for phy_id, sid in PACED_CASES:
                for seed in (args.seeds or [1, 2]):
                    results.append(run_paced_case(phy_id, sid, seed))
    results.extend(run_fixtures(args, ctx))
    elapsed = time.time() - t0

    # ---- summary ----
    counts = {}
    for r in results:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    groups = {}
    for r in results:
        g = (r["scenario"], r["phy"], r["mode"])
        d = groups.setdefault(g, dict(cases=0, frames=0, total=0,
                                      status=set()))
        d["cases"] += 1
        d["frames"] += r.get("delivered", 0)
        d["total"] += r.get("total", 0)
        d["status"].add(r["status"])
    print()
    print("== stress summary (%.1fs, %d workers) ==" % (elapsed, args.jobs))
    print("%-24s %-14s %-4s %6s %10s  %s" %
          ("scenario", "phy", "rx", "cases", "frames", "result"))
    for g in sorted(groups):
        d = groups[g]
        st = d["status"]
        label = ("FAIL" if "FAIL" in st else
                 "XFAIL" if "XFAIL" in st else
                 "NEW" if "NEW" in st else "PASS")
        print("%-24s %-14s %-4s %6d %4d/%-5d  %s" %
              (g[0], g[1], g[2], d["cases"], d["frames"], d["total"], label))
    bad = [r for r in results if r["status"] == "FAIL"]
    xf = [r for r in results if r["status"] == "XFAIL"]
    new = [r for r in results if r["status"] == "NEW"]
    if bad:
        print("\nFAILED cases (%d):" % len(bad))
        for r in bad[:80]:
            print("  %s: %s" % (r["key"], r["reason"]))
        if len(bad) > 80:
            print("  ... %d more" % (len(bad) - 80))
    if xf:
        print("\nXFAIL cases (%d, known/stochastic):" % len(xf))
        for r in xf[:30]:
            print("  %s: %d/%d" % (r["key"], r["delivered"], r["total"]))
        if len(xf) > 30:
            print("  ... %d more" % (len(xf) - 30))
    if new:
        print("\nNEW cases without baseline (%d)" % len(new))
    print("\nTotals: " + ", ".join("%s=%d" % (k, counts[k])
                                   for k in sorted(counts)))

    # ---- failure artifacts ----
    fail_root = args.keep_failures
    if fail_root is None and (bad or args.keep_failures):
        fail_root = tempfile.mkdtemp(prefix="nsfsk-stress-")
    if fail_root:
        os.makedirs(fail_root, exist_ok=True)
        saved = 0
        for r in results:
            if (r["status"] == "FAIL" or
                    (args.keep_failures and r["status"] == "XFAIL")):
                if r.get("artifacts"):
                    d = save_artifacts(r["artifacts"], r.get("commands", {}),
                                       r, fail_root)
                    if d:
                        saved += 1
        if saved:
            print("artifacts: %s" % fail_root)

    # ---- optional JSON ----
    if args.json_out:
        safe_results = [{k: v for k, v in r.items()
                         if k not in ("artifacts", "stdout")}
                        for r in results]
        with open(args.json_out, "w", encoding="utf-8") as f:
            json.dump({"results": safe_results, "elapsed_s": elapsed,
                       "counts": counts}, f, indent=1, default=str)
        print("json: %s" % args.json_out)

    # ---- baseline update ----
    if args.update_baseline:
        if bad:
            print("refusing to update baseline: %d FAIL case(s)" % len(bad),
                  file=sys.stderr)
            return 1
        cases = {}
        for r in results:
            sc = SCENARIOS.get(r["scenario"])
            if sc is None:
                continue
            cases[r["key"]] = {
                "frames": r["delivered"],
                "total": r["total"],
                "tolerance": sc.get("tolerance", 0),
                "mode": "gate" if is_gated(sc, r["phy"]) else "xfail",
                "known_hard": bool(sc.get("known_hard", False) or
                                   r["phy"] in sc.get("known_hard_phys", [])),
            }
        write_baseline(args.baseline, cases)
        print("baseline updated: %s (%d cases)" % (args.baseline, len(cases)))
        return 0

    if args.strict and (xf or new):
        return 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
