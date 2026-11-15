#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Cross-flow and in-line VIV of the Delta Flume vertical tension riser, blind-trial cases 1-9.

The nine stepped-current cases of the blind-prediction trial of Chaplin et al. (2005b): the lower
5.94 m of the 13.12 m riser in a uniform current of 0.16-0.95 m/s, the rest in still water.
Writes a deck per case, sets its unstretched length so that the static top tension in the current
is the case's measured mean, runs it, and compares with the measured response: the dominant
cross-flow mode, the largest cross-flow amplitude, the largest mean in-line deflection and the
overall standard deviations of the in-line and cross-flow curvature.

    python delta_flume_viv.py path/to/cabledyn.exe WORKDIR [--cases 1-9] [--segments 64]
        [--dt-scale 1] [--tmax 40] [--t0 15] [--spring 11.4e3] [--jobs N] [--summary DIR]

Sources:
  J.R. Chaplin, P.W. Bearman, F.J. Huera Huarte, R.J. Pattenden (2005a), Laboratory measurements
  of vortex-induced vibrations of a vertical tension riser in a stepped current, Journal of Fluids
  and Structures 21(1), 3-24, https://doi.org/10.1016/j.jfluidstructs.2005.04.010 (the riser;
  Eq. 4, the drag amplification)
  J.R. Chaplin, P.W. Bearman, Y. Cheng, E. Fontaine, J.M.R. Graham, K. Herfjord, F.J. Huera
  Huarte, M. Isherwood, K. Lambrakos, C.M. Larsen, J.R. Meneghini, G. Moe, R.J. Pattenden,
  M.S. Triantafyllou, R.H.J. Willden (2005b), Blind predictions of laboratory measurements of
  vortex-induced vibrations of a tension riser, Journal of Fluids and Structures 21(1), 25-40,
  https://doi.org/10.1016/j.jfluidstructs.2005.05.016 (the blind-trial cases; the measured
  results are those of the authors' presentation of the trial, Trondheim, October 2004)
  The project's pages at the University of Southampton (archived): the mechanical properties
  (mass, structural damping, axial and top-spring stiffness) and the current zone
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# The riser (Chaplin et al. 2005a, Sec. 2; the Southampton mechanical-properties page): 13.12 m
# between pinned joints, D 28 mm, EI 29.9 N m^2, mass ratio 3.0 with the internal water
# (1.847 kg/m; the page gives 1.85 kg/m and the 12.1 N/m submerged weight), EA 5.88 MN,
# structural damping 0.33 % of critical, top spring 11.4 kN/m (the page gives 11.4 or 38.1).
LENGTH = 13.12  # m
DIAMETER = 0.028  # m
EI = 29.9  # N m^2
EA_RISER = 5.88e6  # N
MASS_RATIO = 3.0
RHO_W = 1000.0  # kg/m^3, fresh water
GRAVITY = 9.80665  # m/s^2
AREA = math.pi * DIAMETER**2 / 4.0
MASS = MASS_RATIO * RHO_W * AREA  # kg/m
SUBMERGED_WEIGHT = (MASS - RHO_W * AREA) * GRAVITY  # N/m
ZETA = 0.0033  # structural damping ratio, as bending damping of a pure bending mode at F_REF
F_REF = 2.0  # Hz (the frequency of the measured damping is not given)
CURRENT_HEIGHT = 5.94  # m: the current zone above the lower joint, still water above
Z_BOTTOM = -14.0  # m, lower joint; the top joint at Z_BOTTOM + LENGTH
WATER_DEPTH = 20.0  # m

# Blind-trial cases (presentation slide 3): towing speed [m/s] and mean top tension [N]
CASES: dict[int, tuple[float, float]] = {
    1: (0.16, 405.0),
    2: (0.21, 407.0),
    3: (0.31, 457.0),
    4: (0.40, 506.0),
    5: (0.54, 598.0),
    6: (0.60, 670.0),
    7: (0.70, 743.0),
    8: (0.85, 923.0),
    9: (0.95, 1002.0),
}
# Measured dominant cross-flow mode numbers (slide 23)
MEASURED_MODES = {1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 6, 7: 7, 8: 7, 9: 7}

# VIV set, not fitted to these measurements: synchronisation lift Cl 0.8 in the band
# cF 0.18 +- dF 0.08; in-line force ClX 0.85, alpha 0.15 (Ulveseter et al. 2018); Cd 1.2, Ca 1.
CD, CA, CL, DF, CF = 1.2, 1.0, 0.8, 0.08, 0.18
CLX, ALPHA = 0.85, 0.15
CV = 0.8  # the VIV AMPLITUDE table holds Cv = Cl constant
VIV_RMS_WINDOW = 0.25  # s, the same averaging time at every dtM


def drag_amplification(ad: float) -> float:
    """Chaplin et al. (2005a, Eq. 4): Cd/Cd0 = 1 + 2.3 s + 1.7 s^2, s = (A/D)/sqrt(2)."""
    s = ad / math.sqrt(2.0)
    return 1.0 + 2.3 * s + 1.7 * s * s


def effective_ea(spring: float) -> float:
    """The top spring in series with the riser, as one axial stiffness of the line."""
    return LENGTH / (1.0 / spring + LENGTH / EA_RISER)


def time_step(u: float, scale: float) -> float:
    """The largest round step at or below 1/200 of the shedding period D/(cF U), times scale."""
    target = DIAMETER / (CF * u) / 200.0
    base = max(s for s in (0.005, 0.0025, 0.002, 0.001, 0.0005) if s <= target)
    return base * scale


def deck_text(
    case: int, l0: float, segments: int, dt: float, tmax: float, ea: float
) -> str:
    u, t_top = CASES[case]
    amp = "".join(
        f"riser {ad:.2f} {CV} {drag_amplification(ad):.5f}\n"
        for ad in (0.1 * k for k in range(16))
    )
    return f"""--------------------- CableDyn Input File ---------------------------------
Delta Flume riser, blind-trial case {case}: U {u} m/s, top tension {t_top:g} N
----------------------- LINE TYPES ------------------------------------------
TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t Cl dF cF
(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-) (-) (-) (-)
riser {DIAMETER} {MASS:.6f} {ea:.6e} 0.0 {EI}|-{ZETA}|{F_REF} {CD} 0.0 {CA} 0.0 {CL} {DF} {CF}
----------------------- POINTS ----------------------------------------------
ID Attachment X Y Z M V CdA CA
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Coupled 0.0 0.0 {Z_BOTTOM + LENGTH:.6f} 0 0 0 0
2 Fixed 0.0 0.0 {Z_BOTTOM:.6f} 0 0 0 0
----------------------- LINES -----------------------------------------------
ID LineType AttachA AttachB UnstrLen NumSegs Outputs
(-) (-) (-) (-) (m) (-) (-)
1 riser 1 2 {l0:.9f} {segments} p
--------------------- VIV ----------------------------------------------
TypeName ClX alpha
(-) (-) (-)
riser {CLX} {ALPHA}
--------------------- VIV AMPLITUDE ------------------------------------
TypeName AD Cv CdAmp
(-) (-) (-) (-)
{amp}--------------------- VIV DAMPING ------------------------------------
TypeName a b c
(-) (-) (-) (-)
riser 0 0 1
---------------------- OPTIONS ----------------------------------------------
{RHO_W} rhoW
{GRAVITY} g
{WATER_DEPTH} WtrDpth
{dt:.6f} dtM
{round(tmax / dt) * dt:.6f} TMax
1 Currents
{VIV_RMS_WINDOW} VIVrmsWindow
---------------------- OUTPUTS -----------------------------------------------
FairTen1
END
------------------------- need this line --------------------------------------
"""


def write_current(folder: Path, u: float) -> None:
    """Uniform U over the lower CURRENT_HEIGHT, still water above (a 2 cm linear step)."""
    zc = Z_BOTTOM + CURRENT_HEIGHT
    (folder / "current_profile.txt").write_text(
        "z ux uy uz\n(m) (m/s) (m/s) (m/s)\n"
        f"{-WATER_DEPTH:.3f} {u} 0 0\n{zc - 0.01:.3f} {u} 0 0\n"
        f"{zc + 0.01:.3f} 0 0 0\n0.0 0 0 0\n",
        encoding="utf-8",
    )


def run(exe: Path, folder: Path, stem: str) -> None:
    env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run(
        [str(exe), f"{stem}.dat", stem],
        cwd=folder,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
        check=False,
    )
    (folder / f"{stem}.console.txt").write_text(
        proc.stdout + proc.stderr, encoding="utf-8"
    )
    if proc.returncode != 0:
        raise SystemExit(f"{folder / stem}: CableDyn exited with {proc.returncode}")


def static_top_tension(folder: Path, stem: str) -> float:
    rows = (folder / f"{stem}.static.out").read_text(encoding="utf-8").splitlines()
    head = rows[1].split()
    first = rows[3].split()  # node 1, End A, the top joint
    return float(first[head.index("Tension")])


def tuned_length(exe: Path, folder: Path, case: int, segments: int, ea: float):
    """Unstretched length whose static top tension in the current is the case's mean tension.

    The start is the still-water length, the exact integral of ds / (1 + T(s)/EA) with T falling
    by the submerged weight down the riser; the current's bow then raises the tension, and each
    iteration corrects the length by the tension error over EA.
    """
    t_target = CASES[case][1]
    t_bottom = t_target - SUBMERGED_WEIGHT * LENGTH
    l0 = ea / SUBMERGED_WEIGHT * math.log((ea + t_target) / (ea + t_bottom))
    for _ in range(12):
        stem = "static"
        (folder / f"{stem}.dat").write_text(
            deck_text(case, l0, segments, 0.005, 0.0, ea), encoding="utf-8"
        )
        run(exe, folder, stem)
        t = static_top_tension(folder, stem)
        if abs(t - t_target) < 0.05:
            return l0, t
        l0 *= 1.0 + (t - t_target) / ea
    raise SystemExit(f"case {case}: the static top tension did not converge ({t} N)")


def response(folder: Path, stem: str, t0: float) -> dict:
    """Response metrics of a run, statistics over t > t0.

    x is in-line (along the current), y cross-flow, both from the line between the joints.
    """
    path = folder / f"{stem}.Line1.p.out"
    a = np.loadtxt(path, skiprows=2)
    keep = a[:, 0] > t0
    t = a[keep, 0]
    xyz = a[keep, 1:].reshape(len(t), -1, 3)
    x, y, z = xyz[:, :, 0], xyz[:, :, 1], xyz[:, :, 2]
    # height over the chord, 0 at the lower joint
    s = (z.mean(axis=0) - Z_BOTTOM) / LENGTH
    half_pp = (y.max(axis=0) - y.min(axis=0)) / 2.0 / DIAMETER
    x_mean = x.mean(axis=0) / DIAMETER
    # dominant cross-flow mode: the largest r.m.s. projection on the pinned-pinned sines
    yd = y - y.mean(axis=0)
    o = np.argsort(s)
    so = s[o]
    proj = [
        float(
            np.sqrt(
                np.mean(
                    np.trapezoid(yd[:, o] * np.sin(n * np.pi * so), so, axis=1) ** 2
                )
            )
        )
        for n in range(1, 16)
    ]
    mode = int(np.argmax(proj)) + 1
    # dominant cross-flow frequency a quarter of the length above the lower joint
    j = int(np.argmin(np.abs(s - 0.25)))
    sig = yd[:, j] * np.hanning(len(t))
    freq = np.fft.rfftfreq(len(t), t[1] - t[0])
    f_cf = float(freq[1 + int(np.argmax(np.abs(np.fft.rfft(sig))[1:]))])
    # curvature from second differences of the node positions (equal elements along the arc);
    # overall standard deviation = r.m.s. over the interior nodes of the standard deviation
    ds = LENGTH / (xyz.shape[1] - 1)
    kx = (x[:, 2:] - 2.0 * x[:, 1:-1] + x[:, :-2]) / ds**2
    ky = (y[:, 2:] - 2.0 * y[:, 1:-1] + y[:, :-2]) / ds**2
    tension = np.loadtxt(folder / f"{stem}.out", skiprows=2)
    tension = tension[tension[:, 0] > t0, 1]
    return {
        "mode": mode,
        "f_cf_hz": f_cf,
        "frequency_resolution_hz": float(freq[1]),
        "max_cf_half_pp_D": float(half_pp.max()),
        "max_il_mean_D": float(x_mean.max()),
        "std_curv_il_kD": float(np.sqrt(np.mean(np.var(kx, axis=0))) * DIAMETER),
        "std_curv_cf_kD": float(np.sqrt(np.mean(np.var(ky, axis=0))) * DIAMETER),
        "mean_top_tension_N": float(tension.mean()),
        "std_top_tension_N": float(tension.std()),
    }


def measured() -> tuple[dict, dict]:
    with open(HERE / "measured_summary.csv", encoding="utf-8") as fh:
        summ = {
            (r["quantity"], int(r["case"])): float(r["value"])
            for r in csv.DictReader(fh)
        }
    env = {}
    for case in CASES:
        with open(
            HERE / "envelopes" / f"case{case}_envelopes.csv", encoding="utf-8"
        ) as fh:
            rows = list(csv.DictReader(fh))
        env[case] = max((float(r["y_max_D"]) - float(r["y_min_D"])) / 2.0 for r in rows)
    return summ, env


def parse_cases(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    if not out or any(c not in CASES for c in out):
        raise SystemExit(f"--cases: cases are 1-9, got {text!r}")
    return sorted(set(out))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("exe")
    ap.add_argument("workdir")
    ap.add_argument("--cases", default="1-9", help="for example 1-9, 5 or 5,9")
    ap.add_argument("--segments", type=int, default=64)
    ap.add_argument(
        "--dt-scale",
        type=float,
        default=1.0,
        help="factor of each case's time step (1: 1/200 of the shedding period, rounded down)",
    )
    ap.add_argument("--tmax", type=float, default=40.0)
    ap.add_argument(
        "--t0", type=float, default=15.0, help="start of the statistics [s]"
    )
    ap.add_argument(
        "--spring", type=float, default=11.4e3, help="top spring stiffness [N/m]"
    )
    ap.add_argument("--jobs", type=int, default=min(9, os.cpu_count() or 1))
    ap.add_argument(
        "--summary",
        default=None,
        help="directory for the comparison CSV and summary JSON",
    )
    a = ap.parse_args()
    found = a.exe if Path(a.exe).is_file() else shutil.which(a.exe)
    if found is None:
        raise SystemExit(f"CableDyn executable not found: {a.exe}")
    exe = Path(found).resolve()
    cases = parse_cases(a.cases)
    ea = effective_ea(a.spring)
    tag = f"n{a.segments}_dt{a.dt_scale:g}_k{a.spring:g}"

    def one(case: int) -> dict:
        u, _ = CASES[case]
        folder = Path(a.workdir).resolve() / f"case{case}_{tag}"
        folder.mkdir(parents=True, exist_ok=True)
        write_current(folder, u)
        l0, t_static = tuned_length(exe, folder, case, a.segments, ea)
        dt = time_step(u, a.dt_scale)
        stem = f"case{case}"
        (folder / f"{stem}.dat").write_text(
            deck_text(case, l0, a.segments, dt, a.tmax, ea), encoding="utf-8"
        )
        run(exe, folder, stem)
        out = {
            "case": case,
            "dt_s": dt,
            "unstretched_length_m": l0,
            "static_top_tension_N": t_static,
        }
        out.update(response(folder, stem, a.t0))
        return out

    with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as pool:
        results = list(pool.map(one, cases))

    summ, env = measured()

    def pct(model: float, meas: float) -> float:
        return round(100.0 * (model / meas - 1.0), 1)

    rows = []
    for r in results:
        c = r["case"]
        u, t_top = CASES[c]
        rows.append(
            {
                "case": c,
                "U_m_per_s": u,
                "top_tension_N": t_top,
                "dt_s": r["dt_s"],
                "mode_measured": MEASURED_MODES[c],
                "mode_cabledyn": r["mode"],
                "f_cf_cabledyn_hz": round(r["f_cf_hz"], 3),
                "max_cf_D_measured": round(summ[("max_cf_disp_D", c)], 4),
                "max_cf_half_pp_D_measured_envelope": round(env[c], 4),
                "max_cf_half_pp_D_cabledyn": round(r["max_cf_half_pp_D"], 4),
                "max_cf_error_percent": pct(
                    r["max_cf_half_pp_D"], summ[("max_cf_disp_D", c)]
                ),
                "max_cf_envelope_error_percent": pct(r["max_cf_half_pp_D"], env[c]),
                "max_il_mean_D_measured": round(summ[("max_il_mean_D", c)], 3),
                "max_il_mean_D_cabledyn": round(r["max_il_mean_D"], 3),
                "max_il_mean_error_percent": pct(
                    r["max_il_mean_D"], summ[("max_il_mean_D", c)]
                ),
                "std_curv_cf_kD_measured": f"{summ[('std_curv_cf_kD', c)]:.4e}",
                "std_curv_cf_kD_cabledyn": f"{r['std_curv_cf_kD']:.4e}",
                "std_curv_cf_error_percent": pct(
                    r["std_curv_cf_kD"], summ[("std_curv_cf_kD", c)]
                ),
                "std_curv_il_kD_measured": f"{summ[('std_curv_il_kD', c)]:.4e}",
                "std_curv_il_kD_cabledyn": f"{r['std_curv_il_kD']:.4e}",
                "std_curv_il_error_percent": pct(
                    r["std_curv_il_kD"], summ[("std_curv_il_kD", c)]
                ),
                "mean_top_tension_N_cabledyn": round(r["mean_top_tension_N"], 1),
            }
        )
        print(
            f"case {c}: mode {r['mode']}/{MEASURED_MODES[c]}, f_CF {r['f_cf_hz']:.2f} Hz, "
            f"A_CF/D {r['max_cf_half_pp_D']:.3f}/{summ[('max_cf_disp_D', c)]:.3f}, "
            f"IL mean/D {r['max_il_mean_D']:.2f}/{summ[('max_il_mean_D', c)]:.2f}, "
            f"sd kD CF {r['std_curv_cf_kD']:.2e}/{summ[('std_curv_cf_kD', c)]:.2e}, "
            f"IL {r['std_curv_il_kD']:.2e}/{summ[('std_curv_il_kD', c)]:.2e}"
        )

    def span(key: str) -> list[float]:
        values = [r[key] for r in rows]
        return [min(values), max(values)]

    il_err = [r["max_il_mean_error_percent"] for r in rows]
    summary = {
        "cases": cases,
        "segments": a.segments,
        "dt_scale": a.dt_scale,
        "tmax_s": a.tmax,
        "statistics_after_s": a.t0,
        "top_spring_N_per_m": a.spring,
        "effective_axial_stiffness_N": round(ea, 1),
        "modes_right": sum(r["mode_cabledyn"] == r["mode_measured"] for r in rows),
        "modes_within_one": sum(
            abs(r["mode_cabledyn"] - r["mode_measured"]) <= 1 for r in rows
        ),
        "il_mean_within_17_percent": sum(abs(e) <= 17.0 for e in il_err),
        "il_mean_error_percent_range": span("max_il_mean_error_percent"),
        "max_cf_error_percent_range": span("max_cf_error_percent"),
        "max_cf_envelope_error_percent_range": span("max_cf_envelope_error_percent"),
        "std_curv_cf_error_percent_range": span("std_curv_cf_error_percent"),
        "std_curv_il_error_percent_range": span("std_curv_il_error_percent"),
        "per_case": rows,
        "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
    }
    print(json.dumps({k: v for k, v in summary.items() if k != "per_case"}, indent=2))
    if a.summary:
        out = Path(a.summary)
        out.mkdir(parents=True, exist_ok=True)
        with open(
            out / "response_comparison.csv", "w", newline="", encoding="utf-8"
        ) as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n"
        )


if __name__ == "__main__":
    main()
