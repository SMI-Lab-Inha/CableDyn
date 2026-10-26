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
  zig-zag mode was pumped there and the line lost lock-in after about 70 s.

usage: test_viv_lockin.py --exe <cabledyn> --output <dir>
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
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
0.0    0.0  0.15
-100.0 0.0  0.15
--------------------- need this line ------------------
"""


def deck(dtm: float) -> str:
    length = 100.0 / (1.0 + 5.0e4 / 5.0e7)
    return f"""--------------------- CableDyn Input File ---------------------------------
VIV lock-in gate: taut line in a uniform cross current
----------------------- LINE TYPES ------------------------------------------
TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF
(name)     (m)    (kg/m)   (N)       (N-s/-)   (-)  (-)   (-)   (-)   (-)   (-)  (-)   (-)
riser      0.10   10.0     5.0e7     -0.8      0.0  1.2   1.0   0.0   0.0   0.8  0.08  0.18
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
current.dat  WaterKin
0.25     VIVrmsWindow
120.0    TMax
---------------------- OUTPUTS -----------------------------------------------
L1N25PZ
L1N12PZ
END
------------------------- need this line --------------------------------------
"""


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
    (out / "current.dat").write_text(CURRENT, encoding="utf-8", newline="\n")
    results = {}
    for dtm in (0.01, 0.002):
        name = f"viv_dt{dtm}"
        (out / f"{name}.dat").write_text(deck(dtm), encoding="utf-8", newline="\n")
        proc = subprocess.run([args.exe, f"{name}.dat", name], cwd=out, capture_output=True, text=True)
        if proc.returncode != 0:
            print(proc.stdout[-2000:], proc.stderr[-2000:])
            check(False, f"dtM {dtm}: the run completes (exit {proc.returncode})")
            continue
        results[dtm] = response(out / f"{name}.out")
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
    print(f"{len(FAILURES)} failure(s)")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
