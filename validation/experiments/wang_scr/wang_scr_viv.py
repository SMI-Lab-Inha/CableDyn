#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Cross-flow VIV of the SJTU truncated steel catenary riser: uniform current and heave.

The model test of Wang et al. (2014) as Thorsen et al. (2017) use it to verify the time-domain
synchronisation model that CableDyn's cross-flow VIV implements, in its two parts:

- current: a 0.2 m/s uniform in-plane current (Thorsen et al. 2017, Sec. 5.1): the r.m.s. of
  the dynamic cross-flow strain along the riser against the 25 measured values (their Fig. 17)
  and the dominant strain frequency 8.25 m from the lower end against the measured 1.00 Hz;
- heave: still water, the upper end moved along the measured top trajectory (Sec. 5.2, period
  5.96 s, their Fig. 19): the r.m.s. cross-flow strain band-passed from 4 times the motion
  frequency to 50 Hz, as Wang et al. (2014) process the measurement, against the 25 measured
  values of their Fig. 22, with and without a VIV FRAME filter.

Writes a deck per run, runs them in parallel (a run whose folder holds the same input files and
executable is reused) and compares. The default runs are the committed record: the two primary
runs, their convergence runs (half the time step, twice the elements) and the frame-filter
sensitivity of the heave case.

    python wang_scr_viv.py path/to/cabledyn.exe WORKDIR [--runs NAME,...] [--jobs N]
        [--keep] [--summary DIR]

Sources:
  J. Wang, S. Fu, R. Baarholm, J. Wu, C.M. Larsen (2014), Fatigue damage of a steel catenary
  riser from vortex-induced vibration caused by vessel motions, Marine Structures 39, 131-156,
  https://doi.org/10.1016/j.marstruc.2014.07.002 (Table 2: the model riser; Sec. 3.1: the
  band-pass of the measured strain)
  M.J. Thorsen, S. Saevik, C.M. Larsen (2017), Non-linear time domain analysis of cross-flow
  vortex-induced vibrations, Marine Structures 51, 134-151,
  https://doi.org/10.1016/j.marstruc.2016.10.007 (Table 1; Secs. 5.1 and 5.2; Figs. 17, 19
  and 22, read by extract_thorsen2017.py)
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

# Thorsen et al. (2017) Table 1 / Wang et al. (2014) Table 2: the truncated model riser. The
# outer coating, PPR pipe and copper core are solid: no internal water is added to the 0.69 kg/m.
# Wang et al. give EA 6.67e5 N; Thorsen et al.'s 6.66e5 N is used.
LENGTH = 23.71  # m
X_TOP = 21.0425  # m, horizontal distance of the top joint from the lower end
H_TOP = 9.0  # m, height of the top joint above the fake seabed
DIAMETER = 0.024  # m
MASS = 0.69  # kg/m
EI = 10.5  # N m^2
EA = 6.66e5  # N
# Thorsen et al. (2017) Sec. 5.1: seabed springs of 1 (kN/m)/m in compression, friction 0.2
SEABED_K_PER_LENGTH = 1000.0  # N/m per m of riser
FRICTION = 0.2
# Their coefficient sets: C_v 0.7 in the current (Re about 3000), 0.8 in heave (Re about 4000);
# f0 0.18, df 0.08, cross-flow only; C_D 1.2, C_M 2 (C_a 1)
CV = {"current": 0.7, "heave": 0.8}
F0, DF, CD, CA = 0.18, 0.08, 1.2, 1.0
CURRENT = 0.2  # m/s, uniform, in the riser plane from the lower end towards the top
WATER_DEPTH = 9.5  # m: the top joint 0.5 m below the surface (not given)
FREQUENCY_STATION = 8.25  # m from the lower end
MEASURED_FREQUENCY = {
    "current": 1.00,
    "heave": 1.0,
}  # Hz: Fig. 18; "close to 1 Hz" (Fig. 23)
RAMP_PERIODS = 2.0  # the top motion is ramped in over two periods
BAND_HIGH = 50.0  # Hz, the upper end of Wang et al.'s band-pass
# The two response branches of the heave case: spectral shares of the band-passed strain 8.25 m
# from the lower end
BRANCHES = {"0.80-0.95": (0.80, 0.95), "0.95-1.10": (0.95, 1.10)}
PART = {
    "current": {"tmax": 150.0, "t0": 30.0, "fig": "thorsen2017_fig17.csv"},
    "heave": {"tmax": 200.0, "t0": 40.0, "fig": "thorsen2017_fig22.csv"},
}
ANALYSIS = "1"  # part of the reuse key: bump when analyse() changes

# name: (part, role, segments, dt, Tf); the first of each part is its primary run
RUNS: dict[str, tuple[str, str, int, float, float]] = {
    "current": ("current", "primary", 200, 0.005, 0.0),
    "current_dt2": ("current", "convergence: half the time step", 200, 0.0025, 0.0),
    "current_n400": ("current", "convergence: twice the elements", 400, 0.005, 0.0),
    "heave": ("heave", "primary (no VIV FRAME)", 200, 0.005, 0.0),
    "heave_dt2": ("heave", "convergence: half the time step", 200, 0.0025, 0.0),
    "heave_n400": ("heave", "convergence: twice the elements", 400, 0.005, 0.0),
    "heave_tf0.35": ("heave", "VIV FRAME sensitivity", 200, 0.005, 0.35),
    "heave_tf0.5": ("heave", "VIV FRAME sensitivity", 200, 0.005, 0.5),
    "heave_tf0.75": ("heave", "VIV FRAME sensitivity", 200, 0.005, 0.75),
    "heave_tf1": ("heave", "VIV FRAME sensitivity", 200, 0.005, 1.0),
    "heave_tf0.5_dt2": (
        "heave",
        "VIV FRAME sensitivity, half the time step",
        200,
        0.0025,
        0.5,
    ),
}


def top_motion() -> dict:
    """The sinusoids fitted to the measured top motion (extract_thorsen2017.py)."""
    return json.loads((HERE / "thorsen2017_fig19_fit.json").read_text(encoding="utf-8"))


def deck_files(part: str, segments: int, dt: float, tf: float) -> dict[str, str]:
    tmax = round(PART[part]["tmax"] / dt) * dt
    z_top = -(WATER_DEPTH - H_TOP)
    frame = (
        "--------------------- VIV FRAME ----------------------------------------\n"
        f"TypeName Tf\n(-) (s)\nscr {tf:g}\n"
        if tf > 0
        else ""
    )
    drive = "1 Currents" if part == "current" else "motion.txt motionFile"
    nodes = segments + 1
    # The cross-flow curvature at every node; the in-plane one at every tenth, to check which
    # component is which (analyse())
    outs = "FairTen1\n" + "".join(f"CurvY1N{j}\n" for j in range(1, nodes + 1))
    outs += "".join(f"CurvX1N{j}\n" for j in range(1, nodes + 1, 10))
    title = (
        f"{CURRENT} m/s uniform current (Thorsen et al. 2017, case 1)"
        if part == "current"
        else "heave-induced VIV in still water (Thorsen et al. 2017, case 2)"
    )
    files = {
        "wang.dat": f"""--------------------- CableDyn Input File ---------------------------------
Wang et al. (2014) truncated SCR, {title}
----------------------- LINE TYPES ------------------------------------------
TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t Cl dF cF
(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-) (-) (-) (-)
scr {DIAMETER} {MASS} {EA:.6e} 0.0 {EI} {CD} 0.0 {CA} 0.0 {CV[part]} {DF} {F0}
----------------------- POINTS ----------------------------------------------
ID Attachment X Y Z M V CdA CA
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Coupled {X_TOP} 0.0 {z_top:.4f} 0 0 0 0
2 Fixed 0.0 0.0 {-WATER_DEPTH:.4f} 0 0 0 0
----------------------- LINES -----------------------------------------------
ID LineType AttachA AttachB UnstrLen NumSegs Outputs
(-) (-) (-) (-) (m) (-) (-)
1 scr 1 2 {LENGTH} {segments} -
{frame}---------------------- OPTIONS ----------------------------------------------
1025.0 rhoW
9.80665 g
{WATER_DEPTH} WtrDpth
{SEABED_K_PER_LENGTH / DIAMETER:.6e} kBot
{FRICTION} frictionMu
{dt:.6f} dtM
{tmax:.6f} TMax
{drive}
---------------------- OUTPUTS -----------------------------------------------
{outs}END
------------------------- need this line --------------------------------------
"""
    }
    if part == "current":
        files["current_profile.txt"] = (
            f"z ux uy uz\n(m) (m/s) (m/s) (m/s)\n{-WATER_DEPTH:.2f} {CURRENT} 0 0\n"
            f"0.0 {CURRENT} 0 0\n"
        )
    else:
        files["motion.txt"] = motion_text(dt, tmax, z_top)
    return files


def motion_text(dt: float, tmax: float, z_top: float) -> str:
    """The top joint's x and z: the fitted sinusoids about the static position, ramped in over
    RAMP_PERIODS periods with a half-cosine, with their exact velocities and accelerations."""
    fit = top_motion()
    w = 2.0 * math.pi / fit["period_s"]
    ramp = RAMP_PERIODS * fit["period_s"]
    rows = ["# time point x y z vx vy vz ax ay az\n"]
    for k in range(round(tmax / dt) + 1):
        t = k * dt
        if t < ramp:
            q = math.pi / ramp
            r, dr, ddr = (
                0.5 * (1.0 - math.cos(q * t)),
                0.5 * q * math.sin(q * t),
                0.5 * q * q * math.cos(q * t),
            )
        else:
            r, dr, ddr = 1.0, 0.0, 0.0
        vals = []
        for comp, base in (("x", X_TOP), ("z", z_top)):
            amp, ph = fit[comp]["amp_m"], fit[comp]["phase_rad"]
            s, c = math.sin(w * t + ph), math.cos(w * t + ph)
            vals.append(
                (
                    base + amp * r * s,
                    amp * (dr * s + r * w * c),
                    amp * (ddr * s + 2.0 * dr * w * c - r * w * w * s),
                )
            )
        (x, vx, ax), (z, vz, az) = vals
        rows.append(
            f"{t:.6f} 1 {x:.9e} 0 {z:.9e} {vx:.9e} 0 {vz:.9e} {ax:.9e} 0 {az:.9e}\n"
        )
    return "".join(rows)


def share(spec: np.ndarray, f: np.ndarray, lo: float, hi: float) -> float | None:
    total = spec.sum()
    return float(spec[(f >= lo) & (f < hi)].sum() / total) if total > 0 else None


def analyse(out: Path, part: str, t0: float) -> dict:
    """r.m.s. cross-flow strain along the riser (arc length s from the lower end) and the
    spectrum of the strain 8.25 m from the lower end, over t > t0.

    The section axes of a line not restrained in torsion start at End A with x the projection
    of global +z normal to the line, so CurvX is the in-plane curvature (the catenary's) and
    CurvY the cross-flow one; this is checked on the static means. The strain is kappa D/2 at
    the outer surface, as Figs. 17 and 22 report it. Current: the dynamic part (mean removed).
    Heave: band-passed from 4 times the motion frequency to 50 Hz, as Wang et al. (2014, Sec.
    3.1) process the measurement.
    """
    with open(out, encoding="utf-8") as fh:
        fh.readline()
        head = fh.readline().split()
        fh.readline()
        a = np.loadtxt(fh)
    t = a[:, 0]
    keep = t > t0
    dt = float(t[1] - t[0])
    nodes = sorted(int(c[7:]) for c in head if c.startswith("CurvY1N"))
    cy = np.array([a[keep, head.index(f"CurvY1N{j}")] for j in nodes])
    sampled = sorted(int(c[7:]) for c in head if c.startswith("CurvX1N"))
    cx_mean = np.abs(
        [a[keep, head.index(f"CurvX1N{j}")].mean() for j in sampled]
    ).mean()
    cy_mean = np.abs(cy[[nodes.index(j) for j in sampled]].mean(axis=1)).mean()
    if not cx_mean > 10.0 * cy_mean:
        raise SystemExit(
            f"{out}: CurvX is not the in-plane curvature ({cx_mean}, {cy_mean})"
        )
    strain = (cy - cy.mean(axis=1, keepdims=True)) * DIAMETER / 2.0
    f = np.fft.rfftfreq(strain.shape[1], dt)
    band = None
    if part == "heave":
        band = (4.0 / top_motion()["period_s"], BAND_HIGH)
        sp = np.fft.rfft(strain, axis=1)
        sp[:, (f < band[0]) | (f > band[1])] = 0.0
        strain = np.fft.irfft(sp, strain.shape[1], axis=1)
    s = LENGTH - (np.array(nodes) - 1) * LENGTH / (len(nodes) - 1)
    order = np.argsort(s)
    s, strain = s[order], strain[order]
    rms = np.sqrt((strain**2).mean(axis=1))
    j = int(np.argmin(np.abs(s - FREQUENCY_STATION)))
    spec = np.abs(np.fft.rfft(strain[j] * np.hanning(strain.shape[1]))) ** 2
    inner = (
        (rms[1:-1] > rms[:-2]) & (rms[1:-1] > rms[2:]) & (rms[1:-1] > 0.5 * rms.max())
    )
    res = {
        "s_m": s.tolist(),
        "rms_strain": rms.tolist(),
        "frequency_station_s_m": float(s[j]),
        "dominant_frequency_hz": float(f[1 + int(np.argmax(spec[1:]))]),
        "frequency_resolution_hz": float(f[1]),
        "rms_peaks": int(inner.sum()),
        "mean_top_tension_N": float(a[keep, head.index("FairTen1")].mean()),
    }
    if band is not None:
        res["band_hz"] = list(band)
        for name, (lo, hi) in BRANCHES.items():
            res[f"share_{name}_hz"] = share(spec, f, lo, hi)
    return res


def run_case(
    exe: Path, folder: Path, files: dict, analyse_out, key: str, keep: bool, sha: str
):
    """Write the input files, run the deck and analyse it. A folder that holds the same input
    files run by the same executable is not run again: its result is reused when the analysis
    key is also the same, and its output is analysed again if it was kept."""
    folder.mkdir(parents=True, exist_ok=True)
    result, out, ran = (
        folder / "result.json",
        folder / "wang.out",
        folder / "executable_sha256.txt",
    )
    same = all(
        (folder / n).is_file() and (folder / n).read_text(encoding="utf-8") == v
        for n, v in files.items()
    )
    if same and result.is_file():
        old = json.loads(result.read_text(encoding="utf-8"))
        if old.pop("reuse_key", None) == f"{sha}:{key}":
            return old
    if not (
        same and out.is_file() and ran.is_file() and ran.read_text().strip() == sha
    ):
        # A run that fails or is interrupted must not leave the previous model's result next
        # to the new inputs, where the next attempt would reuse it.
        result.unlink(missing_ok=True)
        ran.unlink(missing_ok=True)
        out.unlink(missing_ok=True)
        for n, v in files.items():
            (folder / n).write_text(v, encoding="utf-8", newline="\n")
        proc = subprocess.run(
            [str(exe), "wang.dat", "wang"],
            cwd=folder,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            env=dict(os.environ, OMP_NUM_THREADS="1"),
            check=False,
        )
        (folder / "console.txt").write_text(proc.stdout + proc.stderr, encoding="utf-8")
        if proc.returncode != 0:
            raise SystemExit(f"{folder}: CableDyn exited with {proc.returncode}")
        ran.write_text(sha + "\n", encoding="utf-8")
    res = analyse_out(out)
    result.write_text(
        json.dumps({**res, "reuse_key": f"{sha}:{key}"}, indent=1) + "\n",
        encoding="utf-8",
    )
    if not keep:
        out.unlink()
    return res


def figure(name: str) -> dict[str, np.ndarray]:
    """The measured points and Thorsen et al.'s simulated line of Fig. 17 or 22."""
    data: dict[str, list] = {"measured": [], "thorsen2017": []}
    with open(HERE / name, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            data[r["series"]].append((float(r["s_m"]), float(r["rms_strain"])))
    return {k: np.array(sorted(v)) for k, v in data.items()}


def mean_or_none(values) -> float | None:
    values = list(values)
    return float(np.mean(values)) if values else None


def scores(s: np.ndarray, rms: np.ndarray, measured: np.ndarray) -> dict:
    """r.m.s. strain against the measured stations: the largest over the span, the station
    ratios, and their means over the lower and upper halves of the span."""
    at = np.interp(measured[:, 0], s, rms)
    ratio = at / measured[:, 1]
    lower = measured[:, 0] < LENGTH / 2.0
    return {
        "max_rms_strain": float(rms.max()),
        "s_of_max_m": float(s[int(np.argmax(rms))]),
        "max_ratio": float(rms.max() / measured[:, 1].max()),
        "mean_ratio_at_stations": float(ratio.mean()),
        "rms_log_error_at_stations": float(np.sqrt(np.mean(np.log(ratio) ** 2))),
        "mean_ratio_lower_half": mean_or_none(ratio[lower]),
        "mean_ratio_upper_half": mean_or_none(ratio[~lower]),
        "station_ratios": ratio.tolist(),
        "station_values": at.tolist(),
    }


def rounded(x):
    """Six significant digits, so the record does not carry the last bits of the arithmetic."""
    if isinstance(x, float):
        return float(f"{x:.6g}")
    if isinstance(x, dict):
        return {k: rounded(v) for k, v in x.items()}
    if isinstance(x, list):
        return [rounded(v) for v in x]
    return x


def fmt(value: float | None, spec: str = ".3f") -> str:
    """A score for the console; n/a when the selection holds nothing for it."""
    return "n/a" if value is None else format(value, spec)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("exe")
    ap.add_argument("workdir")
    ap.add_argument(
        "--runs", default="all", help="comma list of " + ", ".join(RUNS) + ", or all"
    )
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--keep", action="store_true", help="keep the output tables")
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
    sha = hashlib.sha256(exe.read_bytes()).hexdigest()
    names = (
        list(RUNS)
        if a.runs == "all"
        else [n.strip() for n in a.runs.split(",") if n.strip()]
    )
    if not names or any(n not in RUNS for n in names):
        raise SystemExit(f"--runs: choose from {', '.join(RUNS)}, got {a.runs!r}")
    work = Path(a.workdir).resolve()

    def one(name: str) -> dict:
        part, _, segments, dt, tf = RUNS[name]
        t0 = PART[part]["t0"]
        folder = work / f"{part}_n{segments}_dt{dt:g}_tf{tf:g}"
        files = deck_files(part, segments, dt, tf)
        return run_case(
            exe,
            folder,
            files,
            lambda out: analyse(out, part, t0),
            f"{part}:{t0:g}:{ANALYSIS}",
            a.keep,
            sha,
        )

    with ThreadPoolExecutor(max_workers=max(1, min(a.jobs, len(names)))) as pool:
        results = dict(zip(names, pool.map(one, names), strict=True))

    rows, summary = [], {}
    for part in PART:
        fig = figure(PART[part]["fig"])
        measured, th = fig["measured"], fig["thorsen2017"]
        th_scores = scores(th[:, 0], th[:, 1], measured)
        block = {
            "coefficients": {
                "Cv": CV[part],
                "f0": F0,
                "df": DF,
                "CD": CD,
                "CM": 1.0 + CA,
            },
            "tmax_s": PART[part]["tmax"],
            "statistics_after_s": PART[part]["t0"],
            "measured_max_rms_strain": float(measured[:, 1].max()),
            "measured_dominant_frequency_hz": MEASURED_FREQUENCY[part],
            "thorsen2017": {
                k: v
                for k, v in th_scores.items()
                if k not in ("station_ratios", "station_values")
            },
            "runs": {},
        }
        if part == "heave":
            fit = top_motion()
            block["top_motion"] = {
                "period_s": fit["period_s"],
                "x_amplitude_m": fit["x"]["amp_m"],
                "z_amplitude_m": fit["z"]["amp_m"],
                "x_lag_behind_z_deg": math.degrees(
                    fit["z"]["phase_rad"] - fit["x"]["phase_rad"]
                ),
                "ramp_periods": RAMP_PERIODS,
            }
        for name in (n for n in RUNS if n in results and RUNS[n][0] == part):
            _, role, segments, dt, tf = RUNS[name]
            r = results[name]
            sc = scores(np.array(r["s_m"]), np.array(r["rms_strain"]), measured)
            run = {
                "role": role,
                "segments": segments,
                "dt_s": dt,
                "viv_frame_Tf_s": tf,
                **{
                    k: v
                    for k, v in sc.items()
                    if k not in ("station_ratios", "station_values")
                },
                **{k: v for k, v in r.items() if k not in ("s_m", "rms_strain")},
            }
            block["runs"][name] = run
            th_at = np.interp(measured[:, 0], th[:, 0], th[:, 1])
            for k, (sm, em) in enumerate(measured):
                rows.append(
                    {
                        "run": name,
                        "part": part,
                        "segments": segments,
                        "dt_s": dt,
                        "viv_frame_Tf_s": tf,
                        "station": k + 1,
                        "s_m": f"{sm:.4f}",
                        "measured_rms_strain": f"{em:.4e}",
                        "cabledyn_rms_strain": f"{sc['station_values'][k]:.4e}",
                        "cabledyn_ratio": f"{sc['station_ratios'][k]:.4f}",
                        "thorsen2017_rms_strain": f"{th_at[k]:.4e}",
                        "thorsen2017_ratio": f"{th_at[k] / em:.4f}",
                    }
                )
            print(
                f"{name:16s} max r.m.s. / measured {fmt(run['max_ratio'])} at s "
                f"{run['s_of_max_m']:.2f} m; station mean {fmt(run['mean_ratio_at_stations'])} "
                f"(lower half {fmt(run['mean_ratio_lower_half'])}, upper "
                f"{fmt(run['mean_ratio_upper_half'])}), log error "
                f"{fmt(run['rms_log_error_at_stations'])}; f({FREQUENCY_STATION} m) "
                f"{run['dominant_frequency_hz']:.3f} Hz"
                + "".join(
                    f"; share {b} Hz {fmt(run.get(f'share_{b}_hz'), '.2f')}"
                    for b in BRANCHES
                    if part == "heave"
                )
            )
        if block["runs"]:
            primary = next(n for n in RUNS if RUNS[n][0] == part)
            base = block["runs"].get(primary)
            for name, run in block["runs"].items():
                if name != primary and base is not None:
                    run["max_ratio_change_from_primary"] = (
                        run["max_ratio"] / base["max_ratio"] - 1.0
                    )
            ths = block["thorsen2017"]
            print(
                f"{part}: Thorsen et al. (2017) max {fmt(ths['max_ratio'])}, station mean "
                f"{fmt(ths['mean_ratio_at_stations'])} (lower {fmt(ths['mean_ratio_lower_half'])}, "
                f"upper {fmt(ths['mean_ratio_upper_half'])}), log error "
                f"{fmt(ths['rms_log_error_at_stations'])}"
            )
            summary[part] = block
    summary["executable_sha256"] = sha
    summary = rounded(summary)
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
