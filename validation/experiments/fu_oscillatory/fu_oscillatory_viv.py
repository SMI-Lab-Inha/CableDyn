#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Cross-flow VIV of the SJTU flexible cylinder moved harmonically through still water.

The tests of Fu et al. (2014): a 4 m horizontal flexible cylinder (D 24 mm, pretension 500 N)
carried by a towing carriage on a sinusoidal path across its axis, A_m 0.68 m (KC 178) at three
periods, U_R,max = 2 pi A_m/(T f1 D) of 4, 6.5 and 7.9. Writes a deck per run, runs them in
parallel (a run whose folder holds the same input files and executable is reused) and compares
the cross-flow A/D at mid-span (station CF4) against the instantaneous reduced velocity U_R on
the acceleration and deceleration halves of each half motion cycle with the measured paths of
their Figs. 15-17 (read by extract_fu2014.py), for three published VIV coefficient sets; and the
share of each half-cycle the envelope spends above sqrt(2)/2 of its largest value with the
"lock-in" shares the article gives (about 17, 34 and 50 %).

The default runs are the committed record: every set and case with a VIV FRAME filter of 0.2 s,
their convergence runs at U_R,max 6.5 (half the time step, twice the elements) and the
sensitivity to the filter time (none, 0.1 and 0.3 s).

    python fu_oscillatory_viv.py path/to/cabledyn.exe WORKDIR [--runs NAME,...] [--jobs N]
        [--keep] [--summary DIR]

Sources:
  S. Fu, J. Wang, R. Baarholm, J. Wu, C.M. Larsen (2014), Features of vortex-induced vibration
  in oscillatory flow, Journal of Offshore Mechanics and Arctic Engineering 136(1), 011801,
  https://doi.org/10.1115/1.4025759 (Table 1: the cylinder; Tables 2 and 3: the stations and
  the test matrix; Secs. 4.3 and 4.6; Eqs. 19 and 20; Figs. 15-17 and their captions)
  J.V. Ulveseter, M.J. Thorsen, S. Saevik, C.M. Larsen (2019), Simulating fundamental and higher
  harmonic VIV of slender structures, Applied Ocean Research 90, 101856,
  https://doi.org/10.1016/j.apor.2019.101856 (Table 1)
  J.V. Ulveseter, M.J. Thorsen, S. Saevik, C.M. Larsen (2018), Time domain simulation of riser
  VIV in current and irregular waves, Marine Structures 60, 241-260,
  https://doi.org/10.1016/j.marstruc.2018.04.001 (Table 3)
  J.R. Chaplin, P.W. Bearman, F.J. Huera Huarte, R.J. Pattenden (2005), Laboratory measurements
  of vortex-induced vibrations of a vertical tension riser in a stepped current, Journal of
  Fluids and Structures 21(1), 3-24, https://doi.org/10.1016/j.jfluidstructs.2005.04.010 (Eq. 4)
  M.L. Facchinetti, E. de Langre, F. Biolley (2004), Coupling of structure and wake oscillators
  in vortex-induced vibrations, Journal of Fluids and Structures 19(2), 123-140,
  https://doi.org/10.1016/j.jfluidstructs.2003.12.004
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
GRAVITY = 9.80665  # m/s^2
RHO = 1000.0  # kg/m^3, fresh water (a towing tank; not stated)
# Fu et al. (2014) Table 1
LENGTH = 4.0  # m
DIAMETER = 0.024  # m
MASS = 0.69  # kg/m
EI = 10.5  # N m^2
EA = 6.67e5  # N
PRETENSION = 500.0  # N
F1 = 2.68  # Hz, their calculated first natural frequency in water, the frequency of U_R
CD, CA = (
    1.2,
    1.0,
)  # C_a 1 is the added mass of the article's natural frequencies (its Eq. 19)
AMPLITUDE = 0.68  # m, the largest A_m of Table 3, KC = 2 pi A_m/D = 178
DEPTH = 2.0  # m, the cylinder's submergence (not stated)
WATER_DEPTH = 4.0  # m
SEGMENTS, DT = 40, 0.001  # the axial period (8 ms) in 8 steps
# The three KC 178 cases of Figs. 15-17: U_R,max -> the period of the captions [s]
CASES = {"4": 16.5, "6.5": 10.2, "7.9": 8.45}
RAMP_PERIODS = 1.0  # the motion is ramped in over one period
CYCLES = 8  # motion cycles after the ramp
SKIP_PERIODS = 2.0  # statistics from two periods on (the ramp and one cycle)
BAND = (0.5, 25.0)  # Hz, the band-pass of the cross-flow displacement
UR_STEP = 0.25  # the U_R grid of the comparison
LOCK_IN = (
    math.sqrt(2.0) / 2.0
)  # Sec. 4.6: the envelope above sqrt(2)/2 of its largest value
MEASURED_LOCK_IN = {"4": 0.17, "6.5": 0.34, "7.9": 0.50}  # Sec. 4.6, "about"
ANALYSIS = "1"  # part of the reuse key: bump when an analysis changes


def chaplin_amplification(ad: float) -> float:
    """Chaplin et al. (2005, Eq. 4): Cd/Cd0 = 1 + 2.3 s + 1.7 s^2, s = (A/D)/sqrt(2)."""
    s = ad / math.sqrt(2.0)
    return 1.0 + 2.3 * s + 1.7 * s * s


# The coefficient sets of the FLOTANT and NDP cases: LINE TYPES Cl dF cF; VIV row (ClX, alpha,
# Chh); VIV AMPLITUDE with Chaplin's drag amplification; VIV DAMPING (a, b, c); WAKE OSCILLATOR
# (CL0, St, eps, A).
SETS: dict[str, dict] = {
    "uls19": {
        "label": "Ulveseter et al. (2019), Table 1, Chh 1.5",
        "cl": (0.85, 0.064, 0.144),
        "viv": (0.75, 0.15, 1.5),
        "amp": False,
        "damping": None,
        "wake": None,
    },
    "b1": {
        "label": "generic set (Delta Flume case)",
        "cl": (0.8, 0.08, 0.18),
        "viv": (0.85, 0.15, None),
        "amp": True,
        "damping": (0, 0, 1),
        "wake": None,
    },
    "wake": {
        "label": "Facchinetti et al. (2004) wake oscillator",
        "cl": (0, 0.08, 0.18),
        "viv": None,
        "amp": False,
        "damping": None,
        "wake": (0.3, 0.2, 0.3, 12),
    },
}

TF_PRIMARY = 0.2  # s
# name -> (set, U_R,max, role, segments, dt, Tf)
RUNS: dict[str, tuple] = {}
for _s in SETS:
    for _c in CASES:
        RUNS[f"{_s}_ur{_c}"] = (_s, _c, "primary", SEGMENTS, DT, TF_PRIMARY)
for _s in SETS:
    RUNS[f"{_s}_ur6.5_dt"] = (
        _s,
        "6.5",
        "half the time step",
        SEGMENTS,
        DT / 2.0,
        TF_PRIMARY,
    )
    RUNS[f"{_s}_ur6.5_n80"] = (
        _s,
        "6.5",
        "twice the elements",
        2 * SEGMENTS,
        DT,
        TF_PRIMARY,
    )
for _tf in (0.0, 0.1, 0.3):
    for _s in SETS:
        for _c in CASES:
            RUNS[f"{_s}_ur{_c}_tf{_tf:g}"] = (
                _s,
                _c,
                "VIV FRAME sensitivity",
                SEGMENTS,
                DT,
                _tf,
            )


def viv_sections(hset: str, tf: float) -> str:
    s = SETS[hset]
    text = ""
    if s["viv"]:
        clx, alpha, chh = s["viv"]
        cols, vals = "ClX alpha", f"{clx} {alpha}"
        if chh is not None:
            cols, vals = cols + " Chh", vals + f" {chh}"
        text += (
            "--------------------- VIV ----------------------------------------------\n"
            f"TypeName {cols}\n(-){' (-)' * len(cols.split())}\npipe {vals}\n"
        )
    if s["amp"]:
        text += (
            "--------------------- VIV AMPLITUDE ------------------------------------\n"
            "TypeName AD Cv CdAmp\n(-) (-) (-) (-)\n"
            + "".join(
                f"pipe {0.1 * k:.2f} {s['cl'][0]} {chaplin_amplification(0.1 * k):.5f}\n"
                for k in range(16)
            )
        )
    if s["damping"]:
        a, b, c = s["damping"]
        text += (
            "--------------------- VIV DAMPING ------------------------------------\n"
            f"TypeName a b c\n(-) (-) (-) (-)\npipe {a} {b} {c}\n"
        )
    if s["wake"]:
        text += (
            "--------------------- WAKE OSCILLATOR ----------------------------------\n"
            "TypeName CL0 St eps A\n(-) (-) (-) (-) (-)\n"
            f"pipe {' '.join(f'{v:g}' for v in s['wake'])}\n"
        )
    if tf > 0.0:
        text += (
            "--------------------- VIV FRAME ----------------------------------------\n"
            f"TypeName Tf\n(-) (s)\npipe {tf:g}\n"
        )
    return text


def ramp(t: float, period: float) -> tuple[float, float, float]:
    """The cosine ramp over RAMP_PERIODS periods and its first two derivatives."""
    tr = RAMP_PERIODS * period
    if t >= tr:
        return 1.0, 0.0, 0.0
    w = math.pi / tr
    return (
        0.5 * (1.0 - math.cos(w * t)),
        0.5 * w * math.sin(w * t),
        0.5 * w * w * math.cos(w * t),
    )


def deck_files(
    hset: str, case: str, segments: int, dt: float, tf: float
) -> dict[str, str]:
    period = CASES[case]
    tmax = round(period * (RAMP_PERIODS + CYCLES) / dt) * dt
    l0 = LENGTH / (1.0 + PRETENSION / EA)
    mid = segments // 2 + 1
    cl = " ".join(f"{v:g}" for v in SETS[hset]["cl"])
    deck = f"""--------------------- CableDyn Input File ---------------------------------
Fu et al. (2014) flexible cylinder in oscillatory flow: KC 178, U_R,max {case} (T {period:g} s), set {hset}
----------------------- LINE TYPES ------------------------------------------
TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t Cl dF cF
(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-) (-) (-) (-)
pipe {DIAMETER} {MASS} {EA:.6e} 0.0 {EI} {CD} 0.0 {CA} 0.0 {cl}
----------------------- POINTS ----------------------------------------------
ID Attachment X Y Z M V CdA CA
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Coupled 0.0 0.0 {-DEPTH} 0 0 0 0
2 Coupled {LENGTH} 0.0 {-DEPTH} 0 0 0 0
----------------------- LINES -----------------------------------------------
ID LineType AttachA AttachB UnstrLen NumSegs Outputs
(-) (-) (-) (-) (m) (-) (-)
1 pipe 1 2 {l0:.9f} {segments} -
{viv_sections(hset, tf)}---------------------- OPTIONS ----------------------------------------------
{RHO} rhoW
{GRAVITY} g
{WATER_DEPTH} WtrDpth
{dt:.6f} dtM
{tmax:.6f} TMax
motion.txt motionFile
---------------------- OUTPUTS -----------------------------------------------
FairTen1
L1N{mid}pz
END
------------------------- need this line --------------------------------------
"""
    w = 2.0 * math.pi / period
    lines = ["# time point x y z vx vy vz ax ay az"]
    for k in range(int(round(tmax / dt)) + 1):
        t = k * dt
        r, dr, ddr = ramp(t, period)
        s, c = math.sin(w * t), math.cos(w * t)
        y = AMPLITUDE * r * s
        vy = AMPLITUDE * (dr * s + r * w * c)
        ay = AMPLITUDE * (ddr * s + 2.0 * dr * w * c - r * w * w * s)
        for pid, x0 in ((1, 0.0), (2, LENGTH)):
            lines.append(
                f"{t:.6f} {pid} {x0:.9e} {y:.9e} {-DEPTH:.9e} 0 {vy:.9e} 0 0 {ay:.9e} 0"
            )
    return {"fu.dat": deck, "motion.txt": "\n".join(lines) + "\n"}


def read_output(out: Path) -> tuple[list[str], np.ndarray]:
    lines = out.read_text(encoding="utf-8").splitlines()
    return lines[1].split(), np.loadtxt(lines[3:], ndmin=2)


def envelope(x: np.ndarray, dt: float) -> np.ndarray:
    """Hilbert envelope of x band-passed over BAND (by FFT)."""
    n = len(x)
    f = np.fft.fftfreq(n, dt)
    sp = np.fft.fft(x - x.mean())
    sp[(np.abs(f) < BAND[0]) | (np.abs(f) > BAND[1])] = 0.0
    h = np.zeros(n)
    h[0] = 1.0
    h[1 : (n + 1) // 2] = 2.0
    if n % 2 == 0:
        h[n // 2] = 1.0
    return np.abs(np.fft.ifft(sp * h))


def analyse(out: Path, case: str, segments: int) -> dict:
    """A/D at mid-span on the U_R grid, acceleration and deceleration paths, averaged over the
    half-cycles after SKIP_PERIODS (the last half period dropped, where the envelope ends); the
    lock-in share of each half-cycle."""
    period = CASES[case]
    head, a = read_output(out)
    t = a[:, 0]
    dt = float(t[1] - t[0])
    env = envelope(a[:, head.index(f"L1N{segments // 2 + 1}pz")], dt) / DIAMETER
    w = 2.0 * math.pi / period
    ur = np.abs(AMPLITUDE * w * np.cos(w * t)) / (F1 * DIAMETER)
    accel = np.sin(2.0 * w * t) < 0.0  # |cos| rising
    keep = (t >= SKIP_PERIODS * period) & (t <= t[-1] - 0.5 * period)
    grid = np.arange(UR_STEP, float(ur.max()) + 1e-9, UR_STEP)
    paths = {}
    for name, sel in (("acceleration", accel), ("deceleration", ~accel)):
        vals = []
        for u0 in grid:
            b = keep & sel & (np.abs(ur - u0) <= 0.5 * UR_STEP)
            vals.append(float(env[b].mean()) if b.any() else None)
        paths[name] = vals
    shares, peaks = [], []
    k = int(math.ceil(2.0 * SKIP_PERIODS))
    while (k + 1) * period / 2.0 <= t[-1] - 0.5 * period + 1e-9:
        e = env[(t >= k * period / 2.0) & (t < (k + 1) * period / 2.0)]
        peaks.append(float(e.max()))
        shares.append(float(np.mean(e > LOCK_IN * e.max())))
        k += 1
    return {
        "U_R": [float(u) for u in grid],
        **paths,
        "lock_in_share": float(np.mean(shares)),
        "half_cycle_peak_AD": float(np.mean(peaks)),
        "half_cycles": len(shares),
        "steps": int(len(t)),
    }


def run_case(
    exe: Path, folder: Path, files: dict, analyse_fn, key: str, keep: bool, sha: str
) -> dict:
    """Write the deck, run it and analyse it. A folder that holds the same input files run by
    the same executable is not run again: its result is reused when the analysis key is also
    the same, and its output is analysed again if it was kept."""
    folder.mkdir(parents=True, exist_ok=True)
    result, out, ran = (
        folder / "result.json",
        folder / "fu.out",
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
            [str(exe), "fu.dat", "fu"],
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
    res = analyse_fn(out)
    result.write_text(
        json.dumps({**res, "reuse_key": f"{sha}:{key}"}, indent=1) + "\n",
        encoding="utf-8",
    )
    if not keep:
        out.unlink()
    return res


def measured_paths() -> dict[str, dict[str, dict[float, float]]]:
    """{U_R,max: {path: {U_R: A/D}}} of Figs. 15-17 (extract_fu2014.py)."""
    data: dict = {}
    with open(HERE / "fu2014_hysteresis.csv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            case = f"{float(r['U_R_max']):g}"
            data.setdefault(case, {}).setdefault(r["path"], {})[float(r["U_R"])] = (
                float(r["A_D"])
            )
    return data


def compare(res: dict, meas: dict) -> tuple[dict, list]:
    """Scores of a run against the measured paths of its case: the r.m.s. A/D error over both
    paths, the mean hysteresis gap (deceleration less acceleration) over the U_R where both
    paths exist, and the largest A/D."""
    model = {
        p: dict(zip(res["U_R"], res[p], strict=True))
        for p in ("acceleration", "deceleration")
    }
    errors, gaps_m, gaps_c, rows = [], [], [], []
    for path in ("acceleration", "deceleration"):
        for u, am in sorted(meas[path].items()):
            ac = model[path].get(u)
            rows.append((path, u, am, ac))
            if ac is not None:
                errors.append(ac - am)
    for u, am_acc in sorted(meas["acceleration"].items()):
        am_dec = meas["deceleration"].get(u)
        ac_acc, ac_dec = model["acceleration"].get(u), model["deceleration"].get(u)
        if None not in (am_dec, ac_acc, ac_dec):
            gaps_m.append(am_dec - am_acc)
            gaps_c.append(ac_dec - ac_acc)
    largest = [v for p in model.values() for v in p.values() if v is not None]
    return {
        "rms_AD_error": float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
        "mean_AD_bias": float(np.mean(errors)) if errors else None,
        "mean_hysteresis_gap": float(np.mean(gaps_c)) if gaps_c else None,
        "measured_mean_hysteresis_gap": float(np.mean(gaps_m)) if gaps_m else None,
        "largest_AD": max(largest) if largest else None,
        "measured_largest_AD": max(v for p in meas.values() for v in p.values()),
        "points": len(errors),
    }, rows


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
        "--runs", default="all", help="comma list of run names (see RUNS), or all"
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
        hset, case, _, segments, dt, tf = RUNS[name]
        folder = work / f"{hset}_ur{case}_n{segments}_dt{dt:g}_tf{tf:g}"
        return run_case(
            exe,
            folder,
            deck_files(hset, case, segments, dt, tf),
            lambda out: analyse(out, case, segments),
            f"{case}:{segments}:{ANALYSIS}",
            a.keep,
            sha,
        )

    with ThreadPoolExecutor(max_workers=max(1, min(a.jobs, len(names)))) as pool:
        results = dict(zip(names, pool.map(one, names), strict=True))

    meas = measured_paths()
    rows, runs = [], {}
    for name in (n for n in RUNS if n in results):
        hset, case, role, segments, dt, tf = RUNS[name]
        res = results[name]
        sc, pts = compare(res, meas[case])
        runs[name] = {
            "set": hset,
            "U_R_max": float(case),
            "period_s": CASES[case],
            "role": role,
            "segments": segments,
            "dt_s": dt,
            "viv_frame_Tf_s": tf,
            **sc,
            "lock_in_share": res["lock_in_share"],
            "measured_lock_in_share": MEASURED_LOCK_IN[case],
            "half_cycle_peak_AD": res["half_cycle_peak_AD"],
            "half_cycles": res["half_cycles"],
        }
        for path, u, am, ac in pts:
            rows.append(
                {
                    "run": name,
                    "set": hset,
                    "U_R_max": case,
                    "viv_frame_Tf_s": f"{tf:g}",
                    "segments": segments,
                    "dt_s": f"{dt:g}",
                    "path": path,
                    "U_R": f"{u:g}",
                    "measured_AD": f"{am:.4f}",
                    "cabledyn_AD": "" if ac is None else f"{ac:.4f}",
                }
            )
        r = runs[name]
        print(
            f"{name:22s} A/D r.m.s. error {fmt(r['rms_AD_error'])} (bias {fmt(r['mean_AD_bias'], '+.3f')}); "
            f"hysteresis gap {fmt(r['mean_hysteresis_gap'], '+.3f')} [measured "
            f"{fmt(r['measured_mean_hysteresis_gap'], '+.3f')}]; largest {fmt(r['largest_AD'])} "
            f"[{fmt(r['measured_largest_AD'])}]; lock-in share {r['lock_in_share']:.2f} "
            f"[{MEASURED_LOCK_IN[case]:.2f}]"
        )
    for name, r in runs.items():
        base = runs.get(f"{r['set']}_ur{r['U_R_max']:g}")
        if (
            r["role"] != "primary"
            and base is not None
            and base["largest_AD"]
            and r["largest_AD"]
        ):
            r["largest_AD_change_from_primary"] = (
                r["largest_AD"] / base["largest_AD"] - 1.0
            )
    summary = rounded(
        {
            "sets": {k: v["label"] for k, v in SETS.items()},
            "cases": {
                c: {
                    "period_s": p,
                    "U_R_max_from_period": 2
                    * math.pi
                    * AMPLITUDE
                    / (p * F1 * DIAMETER),
                }
                for c, p in CASES.items()
            },
            "amplitude_m": AMPLITUDE,
            "KC": 2.0 * math.pi * AMPLITUDE / DIAMETER,
            "statistics": {
                "skip_periods": SKIP_PERIODS,
                "cycles_after_ramp": CYCLES,
                "band_hz": list(BAND),
                "U_R_bin": UR_STEP,
                "lock_in_level": LOCK_IN,
            },
            "runs": runs,
            "executable_sha256": sha,
        }
    )
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
