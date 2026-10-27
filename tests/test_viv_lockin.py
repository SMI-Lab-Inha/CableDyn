#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Cross-flow VIV lock-in of a taut line against MoorDyn-F.

A 100 m line (50 segments, D 0.1 m, 10 kg/m, EA 5e7 N, about 50 kN) between two points at 50 m
depth in a uniform 0.15 m/s cross current, lift coefficient 0.8 with MoorDyn's band (dF 0.08,
cF 0.18): the shedding frequency cF U/D = 0.27 Hz sits on the first vertical mode, and the line
locks in. MoorDyn-F (OpenFAST v5 moordyn_driver, dtM 5e-4 s) gives a mid-span A/D of 0.662 at
0.293 Hz over 150-300 s. Gates, over 60-120 s:

* at dtM 0.01 s the mid-span A/D is within 3 % of MoorDyn-F's and the response frequency within
  one spectral bin (1/60 Hz) of it, with the quarter-span amplitude of the first mode;
* at dtM 0.002 s, a step that resolves the mesh zig-zag mode, the line stays locked in with the
  same amplitude (within 2 %): with the nodes' own kinematics in the synchronisation model the
  zig-zag mode was pumped there and the line lost lock-in after about 70 s;
* the same line as a finite-EI cable (EI 1 N m^2, a string to within 1e-6 of its geometric
  stiffness: the cubic-Hermite route with the Gauss-integrated lift) locks in with the EI = 0
  line's mid-span and quarter-span amplitudes (within 1 %) and frequency, at dtM 0.01 s and at
  0.002 s (no smoothing needed on the Hermite route);
* the same finite-EI cable on the multibody march (End A Fixed, and a weighted Connect point on
  two chains 300 m away, which puts the deck on that march) locks in as on the Hermite route
  (within 1 % and one bin);
* a bending-stiff, neutrally buoyant cable (8.1 kg/m, EA 5e9 N, EI 5e7 N m^2, about 51 kN:
  beta = pi^2 EI/(T L^2) ~ 1) against the EI = 0 string tensioned to T (1 + beta), which has
  the same first-mode frequency and shape: at U 0.21 m/s both lock into mode 1 with the same
  mid-span A/D (within 3 %) and frequency; at U 0.30 m/s the string jumps to mode 2 (mid-span
  A/D below 0.15) while the cable, whose bending raises f2/f1 to about 3.15, stays in mode 1
  (mid-span A/D above 0.6).

usage: test_viv_lockin.py --exe <cabledyn> --output <dir>
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from cabledyn import read_output  # noqa: E402

FAILURES: list[str] = []
MD_AD, MD_F = 0.662, 0.293  # MoorDyn-F reference, mid span
D = 0.1


def check(ok: bool, message: str) -> None:
    print(("PASS: " if ok else "FAIL: ") + message)
    if not ok:
        FAILURES.append(message)


CURRENT = """MoorDyn v2 water kinematics file
uniform cross-flow current for the VIV lock-in gate
--------------------------- WAVES -------------------------------------
0                    WaveKinMod  - type of wave input
""                   WaveKinFile - file containing wave elevation time series
0.000000E+00  dtWave      - time step to use in setting up wave kinematics grid (s)
0                    WaveDir     - wave heading (deg)
2                                - X wave input type
-10, 110, 2                      - X wave grid point data
2                                - Y wave input type
-10, 10, 2                       - Y wave grid point data
2                                - Z wave input type
-100, 0, 2                       - Z wave grid point data
--------------------------- CURRENT -------------------------------------
1                    CurrentMod  - type of current input
z-depth     x-current      y-current
(m)           (m/s)         (m/s)
0.0    0.0  {u}
-100.0 0.0  {u}
--------------------- need this line ------------------
"""


def deck(dtm: float, ei: float = 0.0, ea: float = 5.0e7, mass: float = 10.0, tension: float = 5.0e4,
         waterkin: str = "current.dat") -> str:
    length = 100.0 / (1.0 + tension / ea)
    return f"""--------------------- CableDyn Input File ---------------------------------
VIV lock-in gate: taut line in a uniform cross current
----------------------- LINE TYPES ------------------------------------------
TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF
(name)     (m)    (kg/m)   (N)       (N-s/-)   (-)  (-)   (-)   (-)   (-)   (-)  (-)   (-)
riser      0.10   {mass}     {ea}     -0.8      {ei}  1.2   1.0   0.0   0.0   0.8  0.08  0.18
----------------------- POINTS ----------------------------------------------
ID   Attachment   X       Y     Z      M    V     CdA   CA
(#)  (-)          (m)     (m)   (m)    (kg) (m^3) (m^2) (-)
1    Coupled      0.0     0.0   -50.0  0    0     0     0
2    Fixed        100.0   0.0   -50.0  0    0     0     0
----------------------- LINES -----------------------------------------------
ID   LineType   AttachA  AttachB  UnstrLen  NumSegs  Outputs
(#)  (name)     (#)      (#)      (m)       (-)      (-)
1    riser      1        2        {length:.6f}  50       -
---------------------- OPTIONS ----------------------------------------------
{dtm}    dtM
100.0    WtrDpth
1025.0   rhoW
9.80665  g
{waterkin}  WaterKin
0.25     VIVrmsWindow
120.0    TMax
---------------------- OUTPUTS -----------------------------------------------
L1N25PZ
L1N12PZ
END
------------------------- need this line --------------------------------------
"""


def multibody(text: str) -> str:
    """The deck with End A Fixed and, 300 m away, a 5 t Connect point hung on two chains: a deck
    with a dynamic point runs on the multibody march."""
    for old, new in (
        ("1    Coupled      0.0     0.0   -50.0", "1    Fixed        0.0     0.0   -50.0"),
        ("riser      0.10", "chain      0.20   250.0    1.0e9     -1.0      0.0   1.2   1.0   0.0   0.0"
                           "   0.0  0.08  0.18\nriser      0.10"),
        ("2    Fixed        100.0   0.0   -50.0  0    0     0     0",
         "2    Fixed        100.0   0.0   -50.0  0    0     0     0\n"
         "3    Fixed        0.0     300.0 -20.0  0    0     0     0\n"
         "4    Connect      50.0    300.0 -60.0  5000 0     0     0\n"
         "5    Fixed        100.0   300.0 -20.0  0    0     0     0"),
        ("  50       -\n", "  50       -\n2    chain      3        4        70.0  20       -\n"
                                "3    chain      5        4        70.0  20       -\n"),
    ):
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    return text


def response(path: Path) -> tuple[float, float, float]:
    """Mid A/D, mid peak frequency, quarter A/D over 60-120 s."""
    hist = read_output(path)
    t = hist.column(hist.channels[0])
    mask = t >= 60.0
    mid = hist.column("L1N25PZ")[mask]
    qtr = hist.column("L1N12PZ")[mask]
    mid = mid - mid.mean()
    dt = float(np.median(np.diff(t[mask])))
    freq = np.fft.rfftfreq(mid.size, dt)
    power = np.abs(np.fft.rfft(mid * np.hanning(mid.size)))
    return (
        math.sqrt(2.0) * float(mid.std()) / D,
        float(freq[1:][np.argmax(power[1:])]),
        math.sqrt(2.0) * float((qtr - qtr.mean()).std()) / D,
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    out = Path(args.output)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    for u in (0.15, 0.21, 0.30):
        (out / f"current_{u}.dat").write_text(CURRENT.format(u=u), encoding="utf-8", newline="\n")
    # the stiff cable and its equal-frequency string (neutrally buoyant: no sag stiffening)
    beta = math.pi**2 * 5.0e7 / (5.0e4 * 100.0**2)
    stiff = {"ei": 5.0e7, "ea": 5.0e9, "mass": 8.1}
    string = {"ea": 5.0e7, "mass": 8.1, "tension": 5.0e4 * (1.0 + beta)}
    runs = [
        (0.01, "viv_dt0.01", deck(0.01, waterkin="current_0.15.dat")),
        (0.002, "viv_dt0.002", deck(0.002, waterkin="current_0.15.dat")),
        ("fei", "viv_fei", deck(0.01, ei=1.0, waterkin="current_0.15.dat")),
        ("fei_fine", "viv_fei_dt0.002", deck(0.002, ei=1.0, waterkin="current_0.15.dat")),
        ("fei_mb", "viv_fei_multibody", multibody(deck(0.01, ei=1.0, waterkin="current_0.15.dat"))),
    ]
    for u in (0.21, 0.30):
        runs.append((("stiff", u), f"viv_stiff_u{u}", deck(0.01, waterkin=f"current_{u}.dat", **stiff)))
        runs.append((("string", u), f"viv_string_u{u}", deck(0.01, waterkin=f"current_{u}.dat", **string)))
    for _, name, text in runs:
        (out / f"{name}.dat").write_text(text, encoding="utf-8", newline="\n")
    env = dict(os.environ, OMP_NUM_THREADS="1")

    def run(name: str) -> subprocess.CompletedProcess[str]:
        cmd = [args.exe, f"{name}.dat", name]
        return subprocess.run(cmd, cwd=out, capture_output=True, text=True, env=env)

    # the decks are independent: run them side by side (a checked Debug build is ~10x slower)
    with ThreadPoolExecutor(max_workers=max(1, min(4, os.cpu_count() or 1))) as pool:
        procs = list(pool.map(run, [name for _, name, _ in runs]))
    results = {}
    for (key, name, _), proc in zip(runs, procs):
        if proc.returncode != 0:
            print(proc.stdout[-2000:], proc.stderr[-2000:])
            check(False, f"{name}: the run completes (exit {proc.returncode})")
            continue
        results[key] = response(out / f"{name}.out")
    if 0.01 in results:
        ad, f, qtr = results[0.01]
        check(abs(ad - MD_AD) <= 0.03 * MD_AD,
              f"dtM 0.01 s: mid-span A/D {ad:.3f}, MoorDyn-F {MD_AD} (within 3 %)")
        check(abs(f - MD_F) <= 1.0 / 60.0,
              f"dtM 0.01 s: response frequency {f:.3f} Hz, MoorDyn-F {MD_F} Hz (within 1/60 Hz)")
        ratio = math.sin(math.pi * 22.0 / 100.0)
        check(abs(qtr / ad - ratio) <= 0.05,
              f"dtM 0.01 s: quarter/mid amplitude {qtr / ad:.3f}, first mode {ratio:.3f}")
    if 0.01 in results and 0.002 in results:
        a1, a2 = results[0.01][0], results[0.002][0]
        check(abs(a2 - a1) <= 0.02 * a1,
              f"dtM 0.002 s stays locked in: mid-span A/D {a2:.3f} against {a1:.3f} at 0.01 s (within 2 %)")
    if 0.01 in results and "fei" in results:
        (a0, f0, q0), (a1, f1, q1) = results[0.01], results["fei"]
        check(abs(a1 - a0) <= 0.01 * a0 and abs(q1 - q0) <= 0.01 * q0,
              f"finite-EI cable: mid/quarter A/D {a1:.3f}/{q1:.3f}, EI = 0 line {a0:.3f}/{q0:.3f} (within 1 %)")
        check(abs(f1 - f0) <= 1.0 / 60.0,
              f"finite-EI cable: response frequency {f1:.3f} Hz, EI = 0 line {f0:.3f} Hz (within 1/60 Hz)")
    if 0.01 in results and "fei_fine" in results:
        (a0, f0, q0), (a1, f1, q1) = results[0.01], results["fei_fine"]
        check(abs(a1 - a0) <= 0.01 * a0 and abs(q1 - q0) <= 0.01 * q0 and abs(f1 - f0) <= 1.0 / 60.0,
              f"finite-EI cable at dtM 0.002 s: mid/quarter A/D {a1:.3f}/{q1:.3f} at {f1:.3f} Hz, EI = 0 "
              f"line at 0.01 s {a0:.3f}/{q0:.3f} at {f0:.3f} Hz (within 1 % and one bin)")
    if "fei" in results and "fei_mb" in results:
        (a0, f0, q0), (a1, f1, q1) = results["fei"], results["fei_mb"]
        check(abs(a1 - a0) <= 0.01 * a0 and abs(q1 - q0) <= 0.01 * q0 and abs(f1 - f0) <= 1.0 / 60.0,
              f"finite-EI cable on the multibody march: mid/quarter A/D {a1:.3f}/{q1:.3f} at {f1:.3f} Hz, "
              f"Hermite route {a0:.3f}/{q0:.3f} at {f0:.3f} Hz (within 1 % and one bin)")
    if ("stiff", 0.21) in results and ("string", 0.21) in results:
        (a1, f1, _), (a0, f0, _) = results[("stiff", 0.21)], results[("string", 0.21)]
        check(abs(a1 - a0) <= 0.03 * a0 and abs(f1 - f0) <= 1.0 / 60.0,
              f"U 0.21 m/s: stiff cable A/D {a1:.3f} at {f1:.3f} Hz, equal-frequency string {a0:.3f} "
              f"at {f0:.3f} Hz (mode 1, within 3 % and one bin)")
    if ("stiff", 0.30) in results and ("string", 0.30) in results:
        a1, a0 = results[("stiff", 0.30)][0], results[("string", 0.30)][0]
        check(a1 > 0.6 and a0 < 0.15,
              f"U 0.30 m/s: stiff cable stays in mode 1 (mid-span A/D {a1:.3f}), the string jumps to "
              f"mode 2 (mid-span A/D {a0:.3f})")
    print(f"{len(FAILURES)} failure(s)")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
