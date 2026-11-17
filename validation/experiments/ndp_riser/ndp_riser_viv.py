#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Cross-flow VIV of the NDP 38 m high-mode riser in uniform flow, against the measurements.

The bare riser of the Norwegian Deepwater Programme tests (Trim et al. 2005), towed horizontally
at 0.3-2.4 m/s. Writes a deck per coefficient set and tow speed, runs them in parallel and
compares the maximum r.m.s. cross-flow strain along the riser (fundamental and total) with the
measured values of Thorsen et al. (2016, Fig. 8a), and the dominant cross-flow frequency and the
spatial mean of the cross-flow sigma/D with Trim et al. (2005, Figs. 18 and 16).

    python ndp_riser_viv.py path/to/cabledyn.exe WORKDIR [--sets b1,and24b,uls19,th16,th16c]
        [--speeds 0.5,1.0,1.5,2.0,2.4 | all] [--segments 200] [--dt 0.001] [--tmax 25]
        [--t0 10] [--tension 4000] [--kax 3500] [--jobs N] [--keep] [--summary DIR]

Sources:
  A.D. Trim, H. Braaten, H. Lie, M.A. Tognarelli (2005), Experimental investigation of
  vortex-induced vibration of long marine risers, Journal of Fluids and Structures 21(3),
  335-361, https://doi.org/10.1016/j.jfluidstructs.2005.07.014 (the rig, Sec. 2.1; Figs. 16
  and 18)
  E. Fontaine, H. Marcollo, K. Vandiver, M. Triantafyllou, C. Larsen, M. Tognarelli,
  Y. Constantinides, O. Oakley (2011), Reliability based factors of safety for VIV fatigue using
  NDP riser high mode VIV tests, OMAE2011, Vol. 7, 511-522, https://doi.org/10.1115/OMAE2011-49820
  (Table 1: the riser properties)
  M.J. Thorsen, S. Saevik, C.M. Larsen (2016), Time domain simulation of vortex-induced
  vibrations in stationary and oscillating flows, Journal of Fluids and Structures 61, 1-19,
  https://doi.org/10.1016/j.jfluidstructs.2015.11.006 (Table 1, Eq. 1, Figs. 7 and 8)
  J.V. Ulveseter, M.J. Thorsen, S. Saevik, C.M. Larsen (2019), Simulating fundamental and higher
  harmonic VIV of slender structures, Applied Ocean Research 90, 101856,
  https://doi.org/10.1016/j.apor.2019.101856 (Table 1, Sec. 4.1)
  M.L. Andersen, S. Saevik, J. Wu, B.J. Leira, H. Langseth (2024b), Simulating vortex-induced
  vibrations in sheared current by using an empirical time-domain model with adaptive
  parameters, Applied Ocean Research 153, 104284, https://doi.org/10.1016/j.apor.2024.104284
  (Table 3)
  J.R. Chaplin, P.W. Bearman, F.J. Huera Huarte, R.J. Pattenden (2005), Laboratory measurements
  of vortex-induced vibrations of a vertical tension riser in a stepped current, Journal of
  Fluids and Structures 21(1), 3-24, https://doi.org/10.1016/j.jfluidstructs.2005.04.010 (Eq. 4)
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

# The riser (Fontaine et al. 2011, Table 1; Thorsen et al. 2016, Table 1): 38.00 m between
# universal joints, OD 27 mm, EI 598.8 N m^2, EA 8.19e6 N, 0.933 kg/m flooded (as tested).
LENGTH = 38.0  # m
DIAMETER = 0.027  # m
EI = 598.8  # N m^2
EA_RISER = 8.19e6  # N
MASS = 0.933  # kg/m
# Fresh water: the mass ratio 1.62 of Fontaine et al. (2011, Table 1) is 0.933 kg/m over the
# mass of the displaced water at 1006 kg/m^3 (1.63 at 1000, 1.59 at 1025).
RHO_W = 1000.0  # kg/m^3
GRAVITY = 9.80665  # m/s^2
# The rig (Trim et al. 2005, Sec. 2.1): clump weights on inclined pendulum arms set the tension,
# 4-6 kN over the test speeds (not given per test). Modelled as an axial spring in series with
# the riser, T0 in still water; the default k_ax is derived from that range (see the README).
TENSION0 = 4000.0  # N
K_AX = 3500.0  # N/m
DEPTH = 5.0  # m, the riser's depth (not given; the current is uniform)
WATER_DEPTH = 10.0  # m
VIV_RMS_WINDOW = 0.5  # s: 500 steps of the default dtM, held fixed when dt changes
SPEEDS = [round(0.1 * k, 1) for k in range(3, 25)]  # the 22 uniform-flow tests [m/s]
DEFAULT_SPEEDS = "0.5,1.0,1.5,2.0,2.4"


def chaplin_amplification(ad: float) -> float:
    """Chaplin et al. (2005, Eq. 4): Cd/Cd0 = 1 + 2.3 s + 1.7 s^2, s = (A/D)/sqrt(2)."""
    s = ad / math.sqrt(2.0)
    return 1.0 + 2.3 * s + 1.7 * s * s


def thorsen_cv() -> list[tuple[float, float]]:
    """Thorsen et al. (2016, Fig. 7) C_v(y0/D), from extract_thorsen2016.py."""
    with open(HERE / "thorsen2016_fig7_cv.csv", encoding="utf-8") as fh:
        return [(float(r["y0_over_D"]), float(r["Cv"])) for r in csv.DictReader(fh)]


# Coefficient sets: LINE TYPES Cd, Cl dF cF; VIV row (ClX, alpha, Chh); VIV AMPLITUDE table;
# VIV DAMPING (a, b, c). Each as its source gives it; "fitted" says whether the source fitted
# it to these NDP tests.
SETS: dict[str, dict] = {
    "b1": {
        "label": "generic set (Delta Flume case)",
        "fitted": "no",
        "cd": 1.2,
        "cl": (0.8, 0.08, 0.18),
        "viv": (0.85, 0.15, None),
        "amp": "chaplin",
        "damping": (0, 0, 1),
    },
    "and24b": {
        "label": "Andersen et al. (2024b), Table 3, cross-flow average",
        "fitted": "sheared-flow tests of this riser",
        "cd": 1.16,
        "cl": (1.28, 0.066, 0.131),
        "viv": None,
        "amp": None,
        "damping": None,
    },
    "uls19": {
        "label": "Ulveseter et al. (2019), Table 1, Chh 1.5",
        "fitted": "these uniform-flow tests",
        "cd": 1.2,
        "cl": (0.85, 0.064, 0.144),
        "viv": (0.75, 0.15, 1.5),
        "amp": None,
        "damping": None,
    },
    "th16": {
        "label": "Thorsen et al. (2016), Fig. 7 and Eq. 1",
        "fitted": "these uniform and sheared-flow tests",
        "cd": 1.2,
        "cl": (1.0, 0.08, 0.18),  # the C_v table replaces Cl
        "viv": None,
        "amp": "thorsen",
        "damping": (0.31, 0.89, 0),
    },
    "th16c": {
        "label": "th16 with Morison's linear cross-flow damping kept (variant)",
        "fitted": "variant of th16",
        "cd": 1.2,
        "cl": (1.0, 0.08, 0.18),
        "viv": None,
        "amp": "thorsen",
        "damping": (0.31, 0.89, 1),
    },
}


def effective_ea(kax: float) -> float:
    """The rig spring in series with the riser, as one axial stiffness of the line."""
    return LENGTH / (1.0 / kax + LENGTH / EA_RISER)


def deck_text(
    hset: str,
    u: float,
    segments: int,
    dt: float,
    tmax: float,
    tension: float,
    kax: float,
) -> str:
    s = SETS[hset]
    ea = effective_ea(kax)
    l0 = LENGTH / (1.0 + tension / ea)  # unstretched length: T0 in still water
    cl, df, cf = s["cl"]
    viv = ""
    if s["viv"]:
        clx, alpha, chh = s["viv"]
        cols, vals = "ClX alpha", f"{clx} {alpha}"
        if chh is not None:
            cols, vals = cols + " Chh", vals + f" {chh}"
        viv += (
            "--------------------- VIV ----------------------------------------------\n"
            f"TypeName {cols}\n(-){' (-)' * len(cols.split())}\nndp {vals}\n"
        )
    if s["amp"] == "chaplin":  # C_v = Cl constant, Chaplin's drag amplification
        rows = "".join(
            f"ndp {0.1 * k:.2f} {cl} {chaplin_amplification(0.1 * k):.5f}\n"
            for k in range(16)
        )
    elif s["amp"] == "thorsen":
        rows = "".join(f"ndp {ad:.2f} {c:.4f} 1.0\n" for ad, c in thorsen_cv())
    if s["amp"]:
        viv += (
            "--------------------- VIV AMPLITUDE ------------------------------------\n"
            f"TypeName AD Cv CdAmp\n(-) (-) (-) (-)\n{rows}"
        )
    if s["damping"]:
        a, b, c = s["damping"]
        viv += (
            "--------------------- VIV DAMPING ------------------------------------\n"
            f"TypeName a b c\n(-) (-) (-) (-)\nndp {a} {b} {c}\n"
        )
    nodes = segments + 1
    outs = "FairTen1\n" + "".join(
        f"CurvX1N{j}\nCurvY1N{j}\n" for j in range(1, nodes + 1)
    )
    outs += "".join(f"L1N{j}py\nL1N{j}pz\n" for j in range(1, nodes + 1, 4))
    return f"""--------------------- CableDyn Input File ---------------------------------
NDP 38 m riser, uniform {u} m/s, T0 {tension:g} N, k_ax {kax:g} N/m, set {hset}
----------------------- LINE TYPES ------------------------------------------
TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t Cl dF cF
(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-) (-) (-) (-)
ndp {DIAMETER} {MASS} {ea:.6e} 0.0 {EI} {s["cd"]} 0.0 1.0 0.0 {cl} {df} {cf}
----------------------- POINTS ----------------------------------------------
ID Attachment X Y Z M V CdA CA
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Coupled 0.0 0.0 {-DEPTH} 0 0 0 0
2 Fixed {LENGTH} 0.0 {-DEPTH} 0 0 0 0
----------------------- LINES -----------------------------------------------
ID LineType AttachA AttachB UnstrLen NumSegs Outputs
(-) (-) (-) (-) (m) (-) (-)
1 ndp 1 2 {l0:.9f} {segments} -
{viv}---------------------- OPTIONS ----------------------------------------------
{RHO_W} rhoW
{GRAVITY} g
{WATER_DEPTH} WtrDpth
{dt:.6f} dtM
{round(tmax / dt) * dt:.6f} TMax
1 Currents
{VIV_RMS_WINDOW} VIVrmsWindow
---------------------- OUTPUTS -----------------------------------------------
{outs}END
------------------------- need this line --------------------------------------
"""


# Search bands of the dominant frequencies, in f D/U: the cross-flow one covers the
# synchronisation band of every set (0.065-0.26) and excludes a slow swing of the riser below it
# and the third harmonic above it; the in-line one is twice that.
CF_BAND = (0.05, 0.30)
IL_BAND = (0.10, 0.60)
ANALYSIS = "1"  # part of the reuse key: bump when analyse() changes


def spectrum_peak(x: np.ndarray, dt: float, lo: float, hi: float) -> float:
    """Peak of the spectra summed over the rows, between lo and hi [Hz]."""
    xm = x - x.mean(axis=-1, keepdims=True)
    sp = (np.abs(np.fft.rfft(xm, axis=-1)) ** 2).sum(axis=0)
    f = np.fft.rfftfreq(xm.shape[-1], dt)
    sp[(f < lo) | (f > hi)] = 0.0
    return float(f[int(np.argmax(sp))])


def bandpass(x: np.ndarray, dt: float, lo: float, hi: float) -> np.ndarray:
    xm = x - x.mean(axis=-1, keepdims=True)
    f = np.fft.rfftfreq(xm.shape[-1], dt)
    sp = np.fft.rfft(xm, axis=-1)
    sp[..., (f < lo) | (f > hi)] = 0.0
    return np.fft.irfft(sp, xm.shape[-1], axis=-1)


def analyse(out: Path, u: float, t0: float) -> dict:
    """Reduce a run to the compared quantities (statistics over t > t0).

    Cross-flow is vertical (z) and in-line horizontal (y). The section axis x of a line not
    restrained in torsion starts at End A as the projection of global +z normal to the line and
    is carried along it, so CurvX is the cross-flow curvature (the static bow tilts it by
    under 2 degrees) and CurvY the in-line one. The dominant frequencies are the peaks of the
    curvature spectra summed over the nodes, inside CF_BAND and IL_BAND. The strain is
    kappa D/2 at the outer surface; its fundamental part is band-passed to 0.5-1.5 times the
    dominant cross-flow frequency, as Thorsen et al. (2016) filter the measured strain. sigma/D
    is band-passed the same way about the dominant frequency of its direction, as the measured
    modal reconstruction holds the vibration; the raw values keep any slow swing.
    """
    with open(out, encoding="utf-8") as fh:
        fh.readline()
        head = fh.readline().split()
        fh.readline()
        a = np.loadtxt(fh)
    t = a[:, 0]
    k = t > t0
    dt = float(t[1] - t[0])
    s = u / DIAMETER  # f D/U = 1 at s Hz
    nodes = sorted({int(c[7:]) for c in head if c.startswith("CurvX1N")})
    cf = np.array([a[k, head.index(f"CurvX1N{j}")] for j in nodes])
    il = np.array([a[k, head.index(f"CurvY1N{j}")] for j in nodes])
    f_cf = spectrum_peak(cf, dt, CF_BAND[0] * s, CF_BAND[1] * s)
    f_il = spectrum_peak(il, dt, IL_BAND[0] * s, IL_BAND[1] * s)
    eps_cf = (cf - cf.mean(axis=1, keepdims=True)) * DIAMETER / 2.0
    eps_il = (il - il.mean(axis=1, keepdims=True)) * DIAMETER / 2.0
    rms_cf = np.sqrt((eps_cf**2).mean(axis=1))
    rms_cf_f = np.sqrt((bandpass(eps_cf, dt, 0.5 * f_cf, 1.5 * f_cf) ** 2).mean(axis=1))
    rms_il = np.sqrt((eps_il**2).mean(axis=1))
    pn = sorted(
        {int(c[3:-2]) for c in head if c.startswith("L1N") and c.endswith("pz")}
    )
    inner = pn[1:-1]  # the end nodes are pinned
    pz = np.array([a[k, head.index(f"L1N{j}pz")] for j in inner])
    py = np.array([a[k, head.index(f"L1N{j}py")] for j in inner])
    sd_cf = bandpass(pz, dt, 0.5 * f_cf, 1.5 * f_cf).std(axis=1) / DIAMETER
    sd_il = bandpass(py, dt, 0.5 * f_il, 1.5 * f_il).std(axis=1) / DIAMETER
    ten = a[k, head.index("FairTen1")]
    return {
        "mean_tension_N": float(ten.mean()),
        "f_cf_hz": f_cf,
        "f_il_hz": f_il,
        "fD_over_U_cf": f_cf / s,
        "f_cf_displacement_hz": spectrum_peak(pz, dt, CF_BAND[0] * s, CF_BAND[1] * s),
        "f_cf_displacement_above_0p3hz": spectrum_peak(pz, dt, 0.3, np.inf),
        "f_il_displacement_above_0p3hz": spectrum_peak(py, dt, 0.3, np.inf),
        "frequency_resolution_hz": 1.0 / (k.sum() * dt),
        "max_rms_strain_cf_fundamental": float(rms_cf_f.max()),
        "max_rms_strain_cf_total": float(rms_cf.max()),
        "max_rms_strain_il_total": float(rms_il.max()),
        "sigma_over_D_cf_mean": float(sd_cf.mean()),
        "sigma_over_D_il_mean": float(sd_il.mean()),
        "sigma_over_D_cf_raw_mean": float((pz.std(axis=1) / DIAMETER).mean()),
        "sigma_over_D_il_raw_mean": float((py.std(axis=1) / DIAMETER).mean()),
    }


def run_case(
    exe: Path, folder: Path, text: str, u: float, t0: float, keep: bool, sha: str
) -> dict:
    """Write the deck, run it and analyse it. A folder that holds the same deck run by the
    same executable is not run again: its result is reused when the statistics window and
    the analysis are also the same, and its output is analysed again if it was kept."""
    folder.mkdir(parents=True, exist_ok=True)
    deck, result, out = folder / "ndp.dat", folder / "result.json", folder / "ndp.out"
    key = f"{sha}:{t0:g}:{ANALYSIS}"
    same = deck.is_file() and deck.read_text(encoding="utf-8") == text
    if same and result.is_file():
        old = json.loads(result.read_text(encoding="utf-8"))
        if old.pop("reuse_key", None) == key:
            return old
    ran = folder / "executable_sha256.txt"
    if not (
        same and out.is_file() and ran.is_file() and ran.read_text().strip() == sha
    ):
        # A run that fails or is interrupted must not leave the previous model's result next
        # to the new inputs, where the next attempt would reuse it.
        result.unlink(missing_ok=True)
        ran.unlink(missing_ok=True)
        out.unlink(missing_ok=True)
        deck.write_text(text, encoding="utf-8")
        (folder / "current_profile.txt").write_text(
            f"z ux uy uz\n(m) (m/s) (m/s) (m/s)\n{-WATER_DEPTH:.2f} 0 {u} 0\n0.0 0 {u} 0\n",
            encoding="utf-8",
        )
        env = dict(os.environ, OMP_NUM_THREADS="1")
        proc = subprocess.run(
            [str(exe), deck.name, "ndp"],
            cwd=folder,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            env=env,
            check=False,
        )
        (folder / "console.txt").write_text(proc.stdout + proc.stderr, encoding="utf-8")
        if proc.returncode != 0:
            raise SystemExit(f"{folder}: CableDyn exited with {proc.returncode}")
        ran.write_text(sha + "\n", encoding="utf-8")
    res = analyse(out, u, t0)
    result.write_text(
        json.dumps({**res, "reuse_key": key}, indent=1) + "\n", encoding="utf-8"
    )
    if not keep:
        out.unlink()
    return res


def read_csv(name: str) -> list[dict]:
    with open(HERE / name, encoding="utf-8") as fh:
        return list(csv.DictReader(r for r in fh if not r.startswith("#")))


def measured() -> dict:
    """The measurements, keyed by speed rounded to 0.1 m/s."""
    m: dict = {"fig8": {}, "fig16": {}, "fig18": {}}
    for r in read_csv("thorsen2016_fig8_strain.csv"):
        if r["panel"] == "uniform":
            key = (r["series"], round(float(r["U"]), 1))
            m["fig8"][key] = float(r["strain_1e-4"]) * 1e-4
    for r in read_csv("measured_trim2005_fig16.csv"):
        m["fig16"][round(float(r["U"]), 1)] = (float(r["CF"]), float(r["IL"]))
    for r in read_csv("measured_trim2005_fig18.csv"):
        m["fig18"][round(float(r["U"]), 1)] = (float(r["CF"]), float(r["IL"]))
    return m


def rms_log(ratios: list[float]) -> float | None:
    return float(np.sqrt(np.mean(np.log(ratios) ** 2))) if ratios else None


def parse_speeds(text: str) -> list[float]:
    if text == "all":
        return SPEEDS
    out = sorted({round(float(v), 1) for v in text.split(",")})
    if not out or any(v not in SPEEDS for v in out):
        raise SystemExit(
            f"--speeds: the tests are 0.3, 0.4, ..., 2.4 m/s, got {text!r}"
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("exe")
    ap.add_argument("workdir")
    ap.add_argument(
        "--sets", default=",".join(SETS), help="comma list of " + ", ".join(SETS)
    )
    ap.add_argument(
        "--speeds",
        default=DEFAULT_SPEEDS,
        help="comma list of tow speeds [m/s], or all",
    )
    ap.add_argument("--segments", type=int, default=200)
    ap.add_argument("--dt", type=float, default=0.001)
    ap.add_argument("--tmax", type=float, default=25.0)
    ap.add_argument(
        "--t0", type=float, default=10.0, help="start of the statistics [s]"
    )
    ap.add_argument(
        "--tension", type=float, default=TENSION0, help="T0 in still water [N]"
    )
    ap.add_argument("--kax", type=float, default=K_AX, help="rig spring [N/m]")
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
    sets = [s.strip() for s in a.sets.split(",") if s.strip()]
    if not sets or any(s not in SETS for s in sets):
        raise SystemExit(f"--sets: choose from {', '.join(SETS)}, got {a.sets!r}")
    speeds = parse_speeds(a.speeds)
    tag = f"n{a.segments}_dt{a.dt:g}_T{a.tension:g}_k{a.kax:g}_t{a.tmax:g}"
    jobs = [(s, u) for s in sets for u in speeds]
    sha = hashlib.sha256(exe.read_bytes()).hexdigest()

    def one(job: tuple[str, float]) -> dict:
        s, u = job
        folder = Path(a.workdir).resolve() / f"{s}_u{u:.1f}_{tag}"
        text = deck_text(s, u, a.segments, a.dt, a.tmax, a.tension, a.kax)
        res = run_case(exe, folder, text, u, a.t0, a.keep, sha)
        return {"set": s, "U_m_per_s": u, **res}

    with ThreadPoolExecutor(max_workers=max(1, min(a.jobs, len(jobs)))) as pool:
        results = list(pool.map(one, jobs))

    m = measured()
    rows, per_set = [], {}
    for s in sets:
        rs = [r for r in results if r["set"] == s]
        ratio, ratio_u5, ratio_th, ratio_th_u5, ratio_total = [], [], [], [], []
        fdu, fdu_meas, sd, sd_meas, sdi, sdi_meas = [], [], [], [], [], []
        for r in rs:
            u = r["U_m_per_s"]
            ef = m["fig8"][("exp_fundamental", u)]
            et = m["fig8"][("exp_total", u)]
            th = m["fig8"][("sim", u)]
            q = r["max_rms_strain_cf_fundamental"] / ef
            ratio.append(q)
            ratio_th.append(th / ef)
            ratio_total.append(r["max_rms_strain_cf_total"] / et)
            if u >= 0.5:
                ratio_u5.append(q)
                ratio_th_u5.append(th / ef)
            row = {
                "set": s,
                "U_m_per_s": u,
                "mean_tension_N": round(r["mean_tension_N"], 0),
                "strain_fund_measured": f"{ef:.4e}",
                "strain_fund_cabledyn": f"{r['max_rms_strain_cf_fundamental']:.4e}",
                "strain_fund_ratio": round(q, 3),
                "strain_fund_thorsen2016_sim": f"{th:.4e}",
                "strain_total_measured": f"{et:.4e}",
                "strain_total_cabledyn": f"{r['max_rms_strain_cf_total']:.4e}",
                "f_cf_cabledyn_hz": round(r["f_cf_hz"], 3),
                "f_il_cabledyn_hz": round(r["f_il_hz"], 3),
                "fD_U_cf_cabledyn": round(r["fD_over_U_cf"], 4),
                "f_cf_measured_hz": "",
                "f_il_measured_hz": "",
                "sigma_D_cf_cabledyn": round(r["sigma_over_D_cf_mean"], 3),
                "sigma_D_il_cabledyn": round(r["sigma_over_D_il_mean"], 3),
                "sigma_D_cf_measured": "",
                "sigma_D_il_measured": "",
            }
            if u in m["fig18"]:
                fc, fi = m["fig18"][u]
                row["f_cf_measured_hz"], row["f_il_measured_hz"] = fc, fi
                fdu.append(r["fD_over_U_cf"])
                fdu_meas.append(fc * DIAMETER / u)
            if u in m["fig16"]:
                c, i = m["fig16"][u]
                row["sigma_D_cf_measured"], row["sigma_D_il_measured"] = c, i
                sd.append(r["sigma_over_D_cf_mean"])
                sd_meas.append(c)
                sdi.append(r["sigma_over_D_il_mean"])
                sdi_meas.append(i)
            rows.append(row)
        per_set[s] = {
            "label": SETS[s]["label"],
            "fitted_to": SETS[s]["fitted"],
            "speeds": [r["U_m_per_s"] for r in rs],
            "strain_fund_ratio_median": float(np.median(ratio)),
            "strain_fund_ratio_range": [float(min(ratio)), float(max(ratio))],
            "strain_fund_rms_log_error": rms_log(ratio),
            "strain_fund_rms_log_error_U_ge_0p5": rms_log(ratio_u5),
            "thorsen2016_sim_rms_log_error_same_speeds": rms_log(ratio_th),
            "thorsen2016_sim_rms_log_error_U_ge_0p5": rms_log(ratio_th_u5),
            "strain_total_ratio_median": float(np.median(ratio_total)),
            "strain_total_rms_log_error": rms_log(ratio_total),
            "fD_U_cf_mean_at_fig18_speeds": float(np.mean(fdu)) if fdu else None,
            "fD_U_cf_measured_mean_at_fig18_speeds": (
                float(np.mean(fdu_meas)) if fdu else None
            ),
            "fD_U_cf_mean_all_speeds": float(np.mean([r["fD_over_U_cf"] for r in rs])),
            "sigma_D_cf_mean": float(np.mean(sd)) if sd else None,
            "sigma_D_cf_measured_mean": float(np.mean(sd_meas)) if sd else None,
            "sigma_D_il_mean": float(np.mean(sdi)) if sdi else None,
            "sigma_D_il_measured_mean": float(np.mean(sdi_meas)) if sdi else None,
            "mean_tension_N_range": [
                round(min(r["mean_tension_N"] for r in rs)),
                round(max(r["mean_tension_N"] for r in rs)),
            ],
        }
        ps = per_set[s]
        print(
            f"{s}: strain fund/measured median {ps['strain_fund_ratio_median']:.2f} "
            f"({ps['strain_fund_ratio_range'][0]:.2f}-{ps['strain_fund_ratio_range'][1]:.2f}), "
            f"rms log {ps['strain_fund_rms_log_error']:.3f}"
            + (
                f" (U >= 0.5: {ps['strain_fund_rms_log_error_U_ge_0p5']:.3f}; Thorsen's sim "
                f"{ps['thorsen2016_sim_rms_log_error_U_ge_0p5']:.3f})"
                if ratio_u5
                else ""
            )
            + f"; fD/U {ps['fD_U_cf_mean_all_speeds']:.3f}"
            + (
                f" (Fig. 18 speeds {ps['fD_U_cf_mean_at_fig18_speeds']:.3f} vs "
                f"{ps['fD_U_cf_measured_mean_at_fig18_speeds']:.3f})"
                if fdu
                else ""
            )
            + (
                f"; sigma/D CF {ps['sigma_D_cf_mean']:.3f} vs {ps['sigma_D_cf_measured_mean']:.3f}"
                if sd
                else ""
            )
            + f"; T {ps['mean_tension_N_range'][0]}-{ps['mean_tension_N_range'][1]} N"
        )
    summary = {
        "segments": a.segments,
        "dt_s": a.dt,
        "tmax_s": a.tmax,
        "statistics_after_s": a.t0,
        "tension0_N": a.tension,
        "rig_spring_N_per_m": a.kax,
        "effective_axial_stiffness_N": round(effective_ea(a.kax), 1),
        "per_set": per_set,
        "executable_sha256": sha,
    }
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
