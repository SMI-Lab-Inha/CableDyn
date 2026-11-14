#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Cross-flow VIV of the SJTU truncated steel catenary riser in a 0.2 m/s uniform current.

The model test of Wang et al. (2014), as used by Thorsen et al. (2017, Sec. 5.1) to verify the
time-domain synchronisation model that CableDyn's cross-flow VIV implements: writes the deck, runs
it, and compares the r.m.s. of the dynamic cross-flow strain along the riser with the 25 measured
values (Thorsen et al. 2017, Fig. 17) and the dominant strain frequency 8.25 m from the lower end
with the measured 1.00 Hz.

    python wang_scr_uniform_current.py path/to/cabledyn.exe WORKDIR [--segments 200]
        [--dt 0.005] [--tmax 150] [--t0 30] [--summary DIR]

Sources:
  J. Wang, S. Fu, R. Baarholm, J. Wu, C.M. Larsen (2014), Fatigue damage of a steel catenary
  riser from vortex-induced vibration caused by vessel motions, Marine Structures 39, 131-156,
  https://doi.org/10.1016/j.marstruc.2014.07.002 (Table 2: the model riser)
  M.J. Thorsen, S. Saevik, C.M. Larsen (2017), Non-linear time domain analysis of cross-flow
  vortex-induced vibrations, Marine Structures 51, 134-151,
  https://doi.org/10.1016/j.marstruc.2016.10.007 (Table 1, Sec. 5.1, Fig. 17: the measured
  r.m.s. strains, read from the vector figure)
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# Thorsen et al. (2017) Table 1 / Wang et al. (2014) Table 2: the truncated model riser. The
# outer coating, PPR pipe and copper core are solid: no internal water is added to the 0.69 kg/m.
LENGTH = 23.71  # m
X_TOP = 21.0425  # m, horizontal distance of the top joint from the lower end
H_TOP = 9.0  # m, height of the top joint above the fake seabed
DIAMETER = 0.024  # m
MASS = 0.69  # kg/m
EI = 10.5  # N m^2
EA = 6.66e5  # N
# Thorsen et al. (2017) Sec. 5.1: seabed springs of 1 (kN/m)/m in compression, friction 0.2;
# their coefficient set for this case: Cv 0.7 (Re about 3000), f0 0.18, df 0.08, CD 1.2, CM 2
SEABED_K_PER_LENGTH = 1000.0  # N/m per m of riser
FRICTION = 0.2
CV, F0, DF, CD, CA = 0.7, 0.18, 0.08, 1.2, 1.0
CURRENT = 0.2  # m/s, uniform, in the riser plane from the lower end towards the top
WATER_DEPTH = (
    9.5  # m: the top joint 0.5 m below the surface (not given; the current is uniform)
)
MEASURED_FREQUENCY = (
    1.00  # Hz, 8.25 m from the lower end (Thorsen et al. 2017, Fig. 18)
)
FREQUENCY_STATION = 8.25  # m from the lower end


def write_deck(case: Path, segments: int, dt: float, tmax: float) -> Path:
    case.mkdir(parents=True, exist_ok=True)
    (case / "current_profile.txt").write_text(
        f"z ux uy uz\n(m) (m/s) (m/s) (m/s)\n{-WATER_DEPTH:.2f} {CURRENT} 0 0\n0.0 {CURRENT} 0 0\n",
        encoding="utf-8",
    )
    outs = "".join(f"CurvX1N{j}\nCurvY1N{j}\n" for j in range(1, segments + 2))
    text = f"""--------------------- CableDyn Input File ---------------------------------
Wang et al. (2014) truncated SCR, 0.2 m/s uniform current (Thorsen et al. 2017, case 1)
----------------------- LINE TYPES ------------------------------------------
TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t Cl dF cF
(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-) (-) (-) (-)
scr {DIAMETER} {MASS} {EA:.6e} 0.0 {EI} {CD} 0.0 {CA} 0.0 {CV} {DF} {F0}
----------------------- POINTS ----------------------------------------------
ID Attachment X Y Z M V CdA CA
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Coupled {X_TOP} 0.0 {-(WATER_DEPTH - H_TOP):.4f} 0 0 0 0
2 Fixed 0.0 0.0 {-WATER_DEPTH:.4f} 0 0 0 0
----------------------- LINES -----------------------------------------------
ID LineType AttachA AttachB UnstrLen NumSegs Outputs
(-) (-) (-) (-) (m) (-) (-)
1 scr 1 2 {LENGTH} {segments} -
---------------------- OPTIONS ----------------------------------------------
1025.0 rhoW
9.80665 g
{WATER_DEPTH} WtrDpth
{SEABED_K_PER_LENGTH / DIAMETER:.6e} kBot
{FRICTION} frictionMu
{dt:.6f} dtM
{round(tmax / dt) * dt:.6f} TMax
1 Currents
---------------------- OUTPUTS -----------------------------------------------
{outs}END
------------------------- need this line --------------------------------------
"""
    deck = case / "wang_scr.dat"
    deck.write_text(text, encoding="utf-8")
    return deck


def cross_flow_strain(out: Path, t0: float):
    """Arc length from the lower end, r.m.s. dynamic cross-flow strain, the strain histories.

    The cross-flow (out-of-plane) curvature is the section component with no static mean (the
    in-plane one carries the catenary curvature); the strain is kappa D/2 at the outer surface,
    as Fig. 17 reports it.
    """
    lines = out.read_text(encoding="utf-8").splitlines()
    head = lines[1].split()
    a = np.loadtxt(lines[3:])
    t = a[:, 0]
    keep = t > t0
    nodes = sorted({int(c[7:]) for c in head if c.startswith("CurvX1N")})
    cx = np.array([a[keep, head.index(f"CurvX1N{j}")] for j in nodes])
    cy = np.array([a[keep, head.index(f"CurvY1N{j}")] for j in nodes])
    cf = cx if np.abs(cx.mean(axis=1)).mean() < np.abs(cy.mean(axis=1)).mean() else cy
    strain = (cf - cf.mean(axis=1, keepdims=True)) * DIAMETER / 2.0
    s = LENGTH - (np.array(nodes) - 1) * LENGTH / (len(nodes) - 1)
    order = np.argsort(s)
    return (
        s[order],
        np.sqrt((strain[order] ** 2).mean(axis=1)),
        strain[order],
        t[keep][1] - t[keep][0],
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("exe")
    ap.add_argument("workdir")
    ap.add_argument("--segments", type=int, default=200)
    ap.add_argument("--dt", type=float, default=0.005)
    ap.add_argument("--tmax", type=float, default=150.0)
    ap.add_argument("--t0", type=float, default=30.0)
    ap.add_argument(
        "--summary",
        default=None,
        help="directory for the comparison CSV and summary JSON",
    )
    a = ap.parse_args()
    # the run starts in the case folder: resolve the executable against the invocation directory
    # (or PATH, for a bare name) first
    found = a.exe if Path(a.exe).is_file() else shutil.which(a.exe)
    if found is None:
        raise SystemExit(f"CableDyn executable not found: {a.exe}")
    exe = Path(found).resolve()
    case = Path(a.workdir) / f"wang_scr_n{a.segments}_dt{a.dt:g}"
    deck = write_deck(case, a.segments, a.dt, a.tmax)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(
        [str(exe), deck.name, "out"],
        cwd=case,
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    if proc.returncode != 0:
        raise SystemExit(f"CableDyn failed ({proc.returncode}):\n{proc.stderr[-2000:]}")
    s, rms, strain, dt = cross_flow_strain(case / "out.out", a.t0)
    with open(HERE / "measured_rms_strain.csv", encoding="utf-8") as fh:
        measured = [
            (float(r["s_m"]), float(r["rms_strain"])) for r in csv.DictReader(fh)
        ]
    sm = np.array([m[0] for m in measured])
    em = np.array([m[1] for m in measured])
    at = np.interp(sm, s, rms)
    j = int(np.argmin(np.abs(s - FREQUENCY_STATION)))
    spec = np.abs(np.fft.rfft(strain[j] * np.hanning(strain.shape[1])))
    freq = np.fft.rfftfreq(strain.shape[1], dt)
    f_dom = float(freq[1 + int(np.argmax(spec[1:]))])
    summary = {
        "segments": a.segments,
        "dt_s": a.dt,
        "analysis_window_s": [a.t0, a.tmax],
        "frequency_resolution_hz": float(freq[1]),
        "max_rms_strain": float(rms.max()),
        "measured_max_rms_strain": float(em.max()),
        "max_ratio": float(rms.max() / em.max()),
        "mean_ratio_at_stations": float(np.mean(at / em)),
        "rms_log_error_at_stations": float(np.sqrt(np.mean(np.log(at / em) ** 2))),
        "dominant_frequency_hz": f_dom,
        "measured_dominant_frequency_hz": MEASURED_FREQUENCY,
        "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
    }
    print(json.dumps(summary, indent=2))
    if a.summary:
        out = Path(a.summary)
        out.mkdir(parents=True, exist_ok=True)
        with open(
            out / "rms_strain_comparison.csv", "w", newline="", encoding="utf-8"
        ) as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["s_m", "measured_rms_strain", "cabledyn_rms_strain", "ratio"])
            for x, m, c in zip(sm, em, at, strict=True):
                w.writerow([f"{x:.4f}", f"{m:.4e}", f"{c:.4e}", f"{c / m:.4f}"])
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n"
        )


if __name__ == "__main__":
    main()
