#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""VIV of the Southampton untensioned catenary model riser with high curvature in steady current.

The tests of Chaplin & King (2018): a 5.36 m, 56 mm riser of very low bending stiffness, pinned
on the floor of a 2 m deep flume and 0.2 m above the surface 4.30 m away, tensioned only by its
weight and the drag, at two mass ratios (submerged weight 2.9 and 20.9 N/m) in currents up to
1.4 m/s in both directions: positive with the upstream end on the floor, negative with it at
the surface. Writes a deck per run, runs them in parallel (a run whose folder holds the same
input files and executable is reused) and compares with the measurements of the article (read
by extract_chaplin2018.py):
  the mean top tension and its standard deviation against U (Figs. 10 and 13);
  the dominant cross-flow frequency against U (Fig. 9, positive; Fig. 15(a), negative) and its
  spread along the suspended span (uniform in the tests);
  the inclination of the straight part in negative current (Fig. 12(a)) and the drag
  coefficient it implies by the article's Eq. 9, as a comparison only;
for three published VIV coefficient sets and without VIV (C_d 1.2).

The default runs are the committed record: every set at U = 0.2-1.4 m/s in steps of 0.1 m/s for
both mass ratios and directions, the still-water runs, and the convergence runs and the runs
without structural damping and without floor friction (README).

    python chaplin_catenary_viv.py path/to/cabledyn.exe WORKDIR [--runs NAME,...] [--jobs N]
        [--keep] [--summary DIR]

Sources:
  J.R. Chaplin, R. King (2018), Laboratory measurements of the vortex-induced vibrations of an
  untensioned catenary riser with high curvature, Journal of Fluids and Structures 79, 26-38,
  https://doi.org/10.1016/j.jfluidstructs.2018.01.008 (Table 1; Secs. 2-4; Eqs. 9-10;
  Figs. 9, 10, 12, 13 and 15)
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
RHO = 1000.0  # kg/m^3: fresh water (not stated; Table 1's W and m* agree at 1000 for m* 1.12)
# Chaplin & King (2018), Table 1
LENGTH = 5.36  # m, between the pin ends
DIAMETER = 0.056  # m
EI = 1.3  # N m^2, estimated
EA = 7.9e5  # N, estimated
ZETA = 0.03  # structural damping in air, "about 3 % of critical"
WEIGHT = {"1.12": 2.9, "1.84": 20.9}  # N/m submerged, by mass ratio m*
DEPTH = 2.0  # m, the lower pin on the floor
TOP = 0.2  # m, the upper pin above the surface
SPAN = 4.30  # m, horizontal distance between the pins
CD, CA = 1.2, 1.0
FRICTION = (
    0.5  # floor friction coefficient (not given: a silicone sheath on the flume floor)
)
SEGMENTS, DT = 54, 0.0025  # about 0.1 m elements
# The current is ramped up linearly from still water over RAMP s: a static solve in the full
# current finds no equilibrium for some light-riser cases in negative current, where the riser
# lifts off the floor, so every run starts from the same still-water state.
RAMP = 10.0  # s
TMAX, T0 = 110.0, 30.0  # s: run length, statistics from T0
TMAX_NONE, T0_NONE = 40.0, 30.0  # s: the runs without VIV settle to a steady state
WELCH = 20.0  # s, the segment of the averaged spectra (0.05 Hz resolution)
# Hz, where the dominant cross-flow frequency is sought: in the displacement spectrum from 0.3 Hz
# (below the lowest measured frequency, 0.35 Hz; a light riser lifted in negative current also
# wanders slowly sideways), in the acceleration spectrum from 0.15 Hz
BAND = {0: (0.3, 6.0), 4: (0.15, 6.0)}
SPEEDS = [round(0.1 * k, 1) for k in range(2, 15)]  # m/s
DIRECTIONS = {"pos": 1.0, "neg": -1.0}
# Rayleigh damping 3 % at the ends of the band of the measured cross-flow frequencies
RAYLEIGH = (0.3, 3.0)  # Hz
# The points compared with the measured frequencies, by their mean elevation [m]: the dots of
# the insets of Fig. 9 (two near the top, (a); one low down, (b)), at their share of the height
# between the pins (0.65 and 0.37; 0.19), and mid-depth for Fig. 15 ("around mid-depth").
PROBES = {
    "upper_a": -DEPTH + 0.65 * (DEPTH + TOP),
    "upper_b": -DEPTH + 0.37 * (DEPTH + TOP),
    "lower": -DEPTH + 0.19 * (DEPTH + TOP),
    "mid_depth": -0.5 * DEPTH,
}
ANALYSIS = "2"  # part of the reuse key: bump when an analysis changes


def chaplin_amplification(ad: float) -> float:
    """Chaplin et al. (2005, Eq. 4): Cd/Cd0 = 1 + 2.3 s + 1.7 s^2, s = (A/D)/sqrt(2)."""
    s = ad / math.sqrt(2.0)
    return 1.0 + 2.3 * s + 1.7 * s * s


# The coefficient sets of the FLOTANT, NDP and V701 cases: LINE TYPES Cl dF cF; VIV row (ClX,
# alpha, Chh); VIV AMPLITUDE with Chaplin's drag amplification; VIV DAMPING (a, b, c); WAKE
# OSCILLATOR (CL0, St, eps, A). none: no VIV, C_d 1.2.
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
    "none": {
        "label": "no VIV, C_d 1.2",
        "cl": (0, 0.08, 0.18),
        "viv": None,
        "amp": False,
        "damping": None,
        "wake": None,
    },
}


def viv_sections(hset: str) -> str:
    s = SETS[hset]
    text = ""
    if s["viv"]:
        clx, alpha, chh = s["viv"]
        cols, vals = "ClX alpha", f"{clx} {alpha}"
        if chh is not None:
            cols, vals = cols + " Chh", vals + f" {chh}"
        text += (
            "--------------------- VIV ----------------------------------------------\n"
            f"TypeName {cols}\n(-){' (-)' * len(cols.split())}\nriser {vals}\n"
        )
    if s["amp"]:
        text += (
            "--------------------- VIV AMPLITUDE ------------------------------------\n"
            "TypeName AD Cv CdAmp\n(-) (-) (-) (-)\n"
            + "".join(
                f"riser {0.1 * k:.2f} {s['cl'][0]} {chaplin_amplification(0.1 * k):.5f}\n"
                for k in range(16)
            )
        )
    if s["damping"]:
        a, b, c = s["damping"]
        text += (
            "--------------------- VIV DAMPING ------------------------------------\n"
            f"TypeName a b c\n(-) (-) (-) (-)\nriser {a} {b} {c}\n"
        )
    if s["wake"]:
        text += (
            "--------------------- WAKE OSCILLATOR ----------------------------------\n"
            "TypeName CL0 St eps A\n(-) (-) (-) (-) (-)\n"
            f"riser {' '.join(f'{v:g}' for v in s['wake'])}\n"
        )
    return text


def mass(mstar: str) -> float:
    """Mass per metre with the internal water, from the submerged weight of Table 1."""
    return RHO * math.pi * DIAMETER**2 / 4.0 + WEIGHT[mstar] / GRAVITY


def deck_files(
    hset: str,
    mstar: str,
    direction: str,
    u: float,
    segments: int,
    dt: float,
    damping: bool,
    friction: float = FRICTION,
) -> dict[str, str]:
    """The deck: End A the upper pin, End B the lower pin on the floor at x = 0; the current
    along +x (positive: the upstream end on the floor) or -x (negative)."""
    ux = DIRECTIONS[direction] * u
    cl = " ".join(f"{v:g}" for v in SETS[hset]["cl"])
    ei = f"{EI}|-{ZETA}|{RAYLEIGH[0]}|{RAYLEIGH[1]}" if damping else f"{EI}"
    outs = ["FairTen1"] + [f"L1N{j}p{c}" for j in range(1, segments + 2) for c in "xyz"]
    tmax = TMAX_NONE if hset == "none" else TMAX
    deck = f"""--------------------- CableDyn Input File ---------------------------------
Chaplin & King (2018) catenary riser, m* {mstar}, U {ux:+.2f} m/s, set {hset}
----------------------- LINE TYPES ------------------------------------------
TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t Cl dF cF
(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-) (-) (-) (-)
riser {DIAMETER} {mass(mstar):.6f} {EA:.6e} 0.0 {ei} {CD} 0.0 {CA} 0.0 {cl}
----------------------- POINTS ----------------------------------------------
ID Attachment X Y Z M V CdA CA
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Coupled {SPAN} 0.0 {TOP} 0 0 0 0
2 Fixed 0.0 0.0 {-DEPTH} 0 0 0 0
----------------------- LINES -----------------------------------------------
ID LineType AttachA AttachB UnstrLen NumSegs Outputs
(-) (-) (-) (-) (m) (-) (-)
1 riser 1 2 {LENGTH} {segments} -
{viv_sections(hset)}---------------------- OPTIONS ----------------------------------------------
{RHO} rhoW
{GRAVITY} g
{DEPTH} WtrDpth
{friction:g} frictionMu
{dt:.6f} dtM
{tmax:.6f} TMax
2 Currents
---------------------- OUTPUTS -----------------------------------------------
{chr(10).join(outs)}
END
------------------------- need this line --------------------------------------
"""
    # a depth profile in time (Currents 2): uniform over the depth, zero at t = 0, U from RAMP
    # on; the last row past the end of the run (the table repeats past its last sample)
    profile = (
        "Chaplin & King (2018) flume current\n"
        "uniform over the depth, ramped from still water\n\n"
        "depths (m)\n"
        f"{-DEPTH - 0.5:.2f} 0.00\n"
        "t (s), ux at each depth (m/s)\n"
        + "".join(
            f"{t:g} {v:.4f} {v:.4f}\n"
            for t, v in (
                (0.0, 0.0),
                (RAMP, ux),
                (RAMP * math.ceil(tmax / RAMP + 1.0), ux),
            )
        )
    )
    return {"ck.dat": deck, "current_profile_dynamic.txt": profile}


def read_output(out: Path) -> tuple[list[str], np.ndarray]:
    with open(out, encoding="utf-8") as fh:
        fh.readline()
        head = fh.readline().split()
        fh.readline()
        a = np.loadtxt(fh)
    return head, a


def stations() -> list[float]:
    """Arc lengths from the upper pin of the twelve accelerometers: 0.4 m apart (Sec. 2), taken
    as centred on the line (0.48 m from each pin; the article gives no positions)."""
    first = 0.5 * (LENGTH - 11 * 0.4)
    return [first + 0.4 * k for k in range(12)]


def peak_frequency(y: np.ndarray, dt: float, power: int) -> float:
    """The peak of the averaged spectrum (Hann segments of WELCH s, half overlapping) of y,
    weighted by omega^power (4: the acceleration spectrum of a displacement), within its BAND."""
    n = int(round(WELCH / dt))
    win = np.hanning(n)
    f = np.fft.rfftfreq(n, dt)
    sp = np.zeros(len(f))
    for k in range(0, len(y) - n + 1, n // 2):
        seg = y[k : k + n] - y[k : k + n].mean()
        sp += np.abs(np.fft.rfft(seg * win)) ** 2
    sp *= (2.0 * np.pi * f) ** power
    lo, hi = BAND[power]
    m = (f >= lo) & (f <= hi)
    return float(f[m][np.argmax(sp[m])])


def highpassed_rms(y: np.ndarray, dt: float) -> float:
    """r.m.s. of y above the lower edge of the displacement BAND (the vibration without the
    slow sideways wandering)."""
    f = np.fft.rfftfreq(len(y), dt)
    sp = np.fft.rfft(y - y.mean())
    sp[f < BAND[0][0]] = 0.0
    return float(np.fft.irfft(sp, len(y)).std())


def analyse(out: Path, hset: str, direction: str, u: float, segments: int) -> dict:
    """Mean and standard deviation of the top tension; the mean profile; at each accelerometer
    station and probe the cross-flow A/D_rms (all of it and above 0.3 Hz) and the peak
    frequencies of its displacement and acceleration spectra; the inclination of the mean
    profile at mid-depth."""
    head, a = read_output(out)
    t = a[:, 0]
    dt = float(t[1] - t[0])
    t0 = T0 if hset != "none" else T0_NONE
    k = t > t0 + 1e-9
    nodes = range(1, segments + 2)
    col = {h: i for i, h in enumerate(head)}
    ten = a[k, col["FairTen1"]]
    xm = np.array([a[k, col[f"L1N{j}px"]].mean() for j in nodes])
    zm = np.array([a[k, col[f"L1N{j}pz"]].mean() for j in nodes])
    s0 = np.linspace(0.0, LENGTH, segments + 1)  # unstretched arc length of each node
    grounded = zm < -DEPTH + 0.01
    res: dict = {
        "mean_top_tension_N": float(ten.mean()),
        "sd_top_tension_N": float(ten.std()),
        "suspended_length_m": float(s0[~grounded].max()) if (~grounded).any() else 0.0,
        "touchdown_x_m": float(xm[~grounded][-1]) if grounded.any() else None,
    }
    # inclination of the mean profile to the horizontal at mid-depth (z = -DEPTH/2)
    j = int(np.argmin(np.abs(zm + 0.5 * DEPTH)))
    dz = zm[j - 1] - zm[j + 1]
    dx = abs(xm[j - 1] - xm[j + 1])
    res["psi_mid_depth_deg"] = float(math.degrees(math.atan2(abs(dz), dx)))
    res["x_mid_depth_m"] = float(xm[j])
    st = []
    for s in stations():
        n = int(round(s / LENGTH * segments)) + 1
        y = a[k, col[f"L1N{n}py"]]
        row = {
            "s_m": round(s, 3),
            "node": n,
            "x_m": float(xm[n - 1]),
            "z_m": float(zm[n - 1]),
            "suspended": bool(not grounded[n - 1]),
            "AD_rms": float(y.std() / DIAMETER),
            "AD_rms_above_0p3hz": highpassed_rms(y, dt) / DIAMETER,
        }
        if hset != "none":
            row["f_acc_hz"] = peak_frequency(y, dt, 4)
            row["f_disp_hz"] = peak_frequency(y, dt, 0)
        st.append(row)
    res["stations"] = st
    for name, zp in PROBES.items():
        free = np.where(~grounded)[0]
        n = int(free[np.argmin(np.abs(zm[free] - zp))]) + 1
        y = a[k, col[f"L1N{n}py"]]
        p = {"node": n, "x_m": float(xm[n - 1]), "z_m": float(zm[n - 1])}
        p["AD_rms"] = float(y.std() / DIAMETER)
        p["AD_rms_above_0p3hz"] = highpassed_rms(y, dt) / DIAMETER
        if hset != "none":
            p["f_acc_hz"] = peak_frequency(y, dt, 4)
            p["f_disp_hz"] = peak_frequency(y, dt, 0)
        res[name] = p
    res["frequency_resolution_hz"] = 1.0 / WELCH
    res["steps"] = int(len(t))
    return res


def run_case(
    exe: Path, folder: Path, files: dict, analyse_fn, key: str, keep: bool, sha: str
) -> dict:
    """Write the deck, run it and analyse it. A folder that holds the same input files run by
    the same executable is not run again: its result is reused when the analysis key is also
    the same, and its output is analysed again if it was kept."""
    folder.mkdir(parents=True, exist_ok=True)
    result, out, ran = (
        folder / "result.json",
        folder / "ck.out",
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
            [str(exe), "ck.dat", "ck"],
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


# name -> (set, m*, direction, U, role, segments, dt, structural damping, floor friction)
RUNS: dict[str, tuple] = {}
for _m in WEIGHT:
    RUNS[f"none_m{_m}_still"] = (
        "none",
        _m,
        "pos",
        0.0,
        "still water",
        SEGMENTS,
        DT,
        True,
        FRICTION,
    )
for _s in SETS:
    for _m in WEIGHT:
        for _d in DIRECTIONS:
            for _u in SPEEDS:
                RUNS[f"{_s}_m{_m}_{_d}{_u:g}"] = (
                    _s,
                    _m,
                    _d,
                    _u,
                    "primary",
                    SEGMENTS,
                    DT,
                    True,
                    FRICTION,
                )
CONVERGENCE = (("1.84", "pos", 0.8), ("1.12", "neg", 1.2))
for _s in ("uls19", "b1", "wake"):
    for _m, _d, _u in CONVERGENCE:
        _base = f"{_s}_m{_m}_{_d}{_u:g}"
        RUNS[_base + "_dt"] = (
            _s,
            _m,
            _d,
            _u,
            "half the time step",
            SEGMENTS,
            DT / 2,
            True,
            FRICTION,
        )
        RUNS[_base + "_n108"] = (
            _s,
            _m,
            _d,
            _u,
            "twice the elements",
            2 * SEGMENTS,
            DT,
            True,
            FRICTION,
        )
for _s in ("uls19", "b1", "wake"):
    for _m in WEIGHT:
        for _d in DIRECTIONS:
            for _u in (0.4, 0.8, 1.2):
                _base = f"{_s}_m{_m}_{_d}{_u:g}"
                RUNS[_base + "_undamped"] = (
                    _s,
                    _m,
                    _d,
                    _u,
                    "no structural damping",
                    SEGMENTS,
                    DT,
                    False,
                    FRICTION,
                )
                RUNS[_base + "_frictionless"] = (
                    _s,
                    _m,
                    _d,
                    _u,
                    "frictionless floor",
                    SEGMENTS,
                    DT,
                    True,
                    0.0,
                )


def measured() -> dict:
    """{(figure, m*): {U: [(value, sd), ...]}} of the extracted points (extract_chaplin2018.py)."""
    data: dict = {}
    with open(HERE / "measured_chaplin2018.csv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            key = (r["figure"], r["m_star"])
            sd = float(r["sd"]) if r["sd"] else None
            data.setdefault(key, {}).setdefault(round(float(r["U_mps"]), 2), []).append(
                (float(r["value"]), sd)
            )
    return data


def implied_cd(psi_deg: float, u: float, w: float) -> float | None:
    """The drag coefficient of a straight line at the inclination psi by the article's Eq. 9,
    0.5 rho Cd d (U sin psi)^2 = W cos psi."""
    p = math.radians(psi_deg)
    if u <= 0.0 or p <= 0.0:
        return None
    return 2.0 * w * math.cos(p) / (RHO * DIAMETER * (u * math.sin(p)) ** 2)


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


def nearest(values: list[float], v: float) -> float | None:
    return min(values, key=lambda q: abs(q - v)) if values else None


def compare_run(name: str, res: dict, meas: dict, still: dict) -> dict:
    """One run against the measured points at its speed."""
    hset, mstar, direction, u, role, segments, dt, damp, mu = RUNS[name]
    fig_t = "10" if direction == "pos" else "13"
    row: dict = {
        "run": name,
        "set": hset,
        "role": role,
        "m_star": mstar,
        "direction": direction,
        "U_mps": u,
        "segments": segments,
        "dt_s": dt,
        "structural_damping": int(damp),
        "floor_friction": mu,
        "model_mean_T_N": res["mean_top_tension_N"],
        "model_sd_T_N": res["sd_top_tension_N"],
        "model_still_T_N": still[mstar],
    }
    t = meas.get((fig_t, mstar), {}).get(round(u, 2), [])
    row["measured_mean_T_N"] = t[0][0] if t else None
    row["measured_sd_T_N"] = t[0][1] if t else None
    t0 = meas.get((fig_t, mstar), {}).get(0.0, [])
    row["measured_still_T_N"] = t0[0][0] if t0 else None
    row["model_psi_mid_depth_deg"] = res["psi_mid_depth_deg"]
    row["model_cd_eq9"] = implied_cd(res["psi_mid_depth_deg"], u, WEIGHT[mstar])
    p = meas.get(("12a", mstar), {}).get(round(u, 2), []) if direction == "neg" else []
    row["measured_psi_deg"] = p[0][0] if p else None
    row["measured_cd_eq9"] = implied_cd(p[0][0], u, WEIGHT[mstar]) if p else None
    if hset != "none":
        if direction == "pos":
            probes = {"upper": ("upper_a", "upper_b"), "lower": ("lower",)}
            figs = {"upper": "9a", "lower": "9b"}
        else:
            probes, figs = {"mid": ("mid_depth",)}, {"mid": "15a"}
        for where, names in probes.items():
            # the dominant frequency: the peak of the displacement spectrum (the article's
            # acceleration and video-displacement peaks agree, its Fig. 4); the acceleration
            # spectrum's peak as a second value
            row[f"model_f_{where}_hz"] = float(
                np.mean([res[q]["f_disp_hz"] for q in names])
            )
            row[f"model_f_acc_{where}_hz"] = float(
                np.mean([res[q]["f_acc_hz"] for q in names])
            )
            pts = [
                v for v, _ in meas.get((figs[where], mstar), {}).get(round(u, 2), [])
            ]
            row[f"measured_f_{where}_hz"] = ";".join(f"{v:.3f}" for v in sorted(pts))
            row[f"nearest_measured_f_{where}_hz"] = nearest(
                pts, row[f"model_f_{where}_hz"]
            )
        fs = [
            s["f_disp_hz"] for s in res["stations"] if s["suspended"] and s["z_m"] < 0.0
        ]
        row["model_f_span_min_hz"] = min(fs) if fs else None
        row["model_f_span_median_hz"] = float(np.median(fs)) if fs else None
        row["model_f_span_max_hz"] = max(fs) if fs else None
        row["model_AD_rms_max"] = max(s["AD_rms"] for s in res["stations"])
        row["model_AD_rms_above_0p3hz_max"] = max(
            s["AD_rms_above_0p3hz"] for s in res["stations"]
        )
    return row


def scores(rows: list[dict]) -> dict:
    """Scores of a set at one mass ratio and direction over the primary speeds."""
    out: dict = {}
    t = [
        (r["model_mean_T_N"], r["measured_mean_T_N"])
        for r in rows
        if r["measured_mean_T_N"]
    ]
    if t:
        e = np.array([m / x - 1.0 for m, x in t])
        out["T_rms_rel_error"] = float(np.sqrt(np.mean(e**2)))
        out["T_mean_rel_error"] = float(np.mean(e))
        out["T_points"] = len(t)
    inc = [
        (
            r["model_mean_T_N"] - r["model_still_T_N"],
            r["measured_mean_T_N"] - r["measured_still_T_N"],
        )
        for r in rows
        if r["measured_mean_T_N"] and r["measured_still_T_N"]
    ]
    if inc:
        out["dT_rms_error_N"] = float(np.sqrt(np.mean([(m - x) ** 2 for m, x in inc])))
        out["dT_mean_error_N"] = float(np.mean([m - x for m, x in inc]))
    sd = [
        r["model_sd_T_N"] / r["measured_sd_T_N"] for r in rows if r["measured_sd_T_N"]
    ]
    if sd:
        out["sd_T_ratio_median"] = float(np.median(sd))
        out["sd_T_points"] = len(sd)
    psi = [
        r["model_psi_mid_depth_deg"] - r["measured_psi_deg"]
        for r in rows
        if r["measured_psi_deg"] is not None
    ]
    if psi:
        out["psi_rms_error_deg"] = float(np.sqrt(np.mean(np.square(psi))))
        out["psi_mean_error_deg"] = float(np.mean(psi))
        out["cd_eq9_model_median"] = float(
            np.median([r["model_cd_eq9"] for r in rows if r["model_cd_eq9"]])
        )
        out["cd_eq9_measured_median"] = float(
            np.median([r["measured_cd_eq9"] for r in rows if r["measured_cd_eq9"]])
        )
    for where in ("upper", "lower", "mid"):
        k = f"model_f_{where}_hz"
        if not rows or k not in rows[0]:
            continue
        # runs whose acceleration spectrum peaks at a higher harmonic of the displacement's
        out[f"f_{where}_acc_peak_at_harmonic"] = int(
            sum(r[f"model_f_acc_{where}_hz"] >= 1.5 * r[k] for r in rows)
        )
        out[f"f_{where}_model_fd_over_U_median"] = float(
            np.median([r[k] * DIAMETER / r["U_mps"] for r in rows])
        )
        meas_fd = [
            float(v) * DIAMETER / r["U_mps"]
            for r in rows
            for v in r[f"measured_f_{where}_hz"].split(";")
            if v
        ]
        if meas_fd:
            out[f"f_{where}_measured_fd_over_U_median"] = float(np.median(meas_fd))
            out[f"f_{where}_measured_points"] = len(meas_fd)
        if where == "mid":
            # on the normal velocity U sin(psi) of the straight part (Fig. 15(b)): the model's
            # inclination at mid-depth, the measured one of Fig. 12(a)
            out["f_mid_model_fd_over_Un_median"] = float(
                np.median(
                    [
                        r[k]
                        * DIAMETER
                        / (
                            r["U_mps"]
                            * math.sin(math.radians(r["model_psi_mid_depth_deg"]))
                        )
                        for r in rows
                    ]
                )
            )
            mn = [
                float(v)
                * DIAMETER
                / (r["U_mps"] * math.sin(math.radians(r["measured_psi_deg"])))
                for r in rows
                if r["measured_psi_deg"]
                for v in r["measured_f_mid_hz"].split(";")
                if v
            ]
            if mn:
                out["f_mid_measured_fd_over_Un_median"] = float(np.median(mn))
        sel = [r for r in rows if r[f"nearest_measured_f_{where}_hz"]]
        if sel:
            ratio = [r[k] / r[f"nearest_measured_f_{where}_hz"] for r in sel]
            out[f"f_{where}_ratio_to_nearest_measured_median"] = float(np.median(ratio))
            out[f"f_{where}_within_10pct"] = int(
                sum(abs(q - 1.0) <= 0.1 for q in ratio)
            )
            out[f"f_{where}_points"] = len(sel)
    sp = [
        (r["model_f_span_max_hz"] - r["model_f_span_min_hz"])
        / r["model_f_span_median_hz"]
        for r in rows
        if r.get("model_f_span_median_hz")
    ]
    if sp:
        out["f_span_spread_median"] = float(np.median(sp))
    return out


def change(r: dict, base: dict) -> dict:
    """Relative change of a sensitivity or convergence run from its primary run."""

    def rel(k):
        if r.get(k) is None or not base.get(k):
            return None
        return r[k] / base[k] - 1.0

    ch = {
        "role": r["role"],
        "mean_T": rel("model_mean_T_N"),
        "sd_T": rel("model_sd_T_N"),
    }
    ch["AD_rms_max"] = rel("model_AD_rms_max")
    for w in ("upper", "lower", "mid"):
        if f"model_f_{w}_hz" in r:
            ch[f"f_{w}"] = rel(f"model_f_{w}_hz")
    return ch


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
    # the still-water runs are the reference of the tension increments: always run
    names = [n for n in RUNS if n.endswith("_still") and n not in names] + names
    work = Path(a.workdir).resolve()

    def one(name: str) -> dict:
        hset, mstar, direction, u, _, segments, dt, damp, mu = RUNS[name]
        folder = work / (
            f"{hset}_m{mstar}_{direction}{u:g}_n{segments}_dt{dt:g}_mu{mu:g}"
            + ("" if damp else "_undamped")
        )
        return run_case(
            exe,
            folder,
            deck_files(hset, mstar, direction, u, segments, dt, damp, mu),
            lambda out: analyse(out, hset, direction, u, segments),
            f"{segments}:{ANALYSIS}",
            a.keep,
            sha,
        )

    with ThreadPoolExecutor(max_workers=max(1, min(a.jobs, len(names)))) as pool:
        results = dict(zip(names, pool.map(one, names), strict=True))

    meas = measured()
    still = {
        RUNS[n][1]: results[n]["mean_top_tension_N"]
        for n in results
        if n.endswith("_still")
    }
    rows = [compare_run(n, results[n], meas, still) for n in RUNS if n in results]
    by = {r["run"]: r for r in rows}
    groups: dict = {}
    for r in rows:
        if r["role"] == "primary":
            groups.setdefault(f"{r['set']}_m{r['m_star']}_{r['direction']}", []).append(
                r
            )
    summary_scores = {k: scores(v) for k, v in groups.items()}
    for k, sc in summary_scores.items():
        freq = "; ".join(
            f"f_{w} fd/U {fmt(sc.get(f'f_{w}_model_fd_over_U_median'))} [nearest measured "
            f"ratio {fmt(sc.get(f'f_{w}_ratio_to_nearest_measured_median'), '.2f')}]"
            for w in ("upper", "lower", "mid")
            if f"f_{w}_model_fd_over_U_median" in sc
        )
        print(
            f"{k:18s} T r.m.s. rel. error {fmt(sc.get('T_rms_rel_error'))} "
            f"(mean {fmt(sc.get('T_mean_rel_error'), '+.3f')}); increment r.m.s. "
            f"{fmt(sc.get('dT_rms_error_N'), '.2f')} N; sd ratio "
            f"{fmt(sc.get('sd_T_ratio_median'), '.2f')}; {freq}"
        )
    changes = {
        r["run"]: change(r, by[r["run"].rsplit("_", 1)[0]])
        for r in rows
        if r["role"] not in ("primary", "still water")
        and r["run"].rsplit("_", 1)[0] in by
    }
    summary = rounded(
        {
            "sets": {k: v["label"] for k, v in SETS.items()},
            "model": {
                "length_m": LENGTH,
                "diameter_m": DIAMETER,
                "EI_Nm2": EI,
                "EA_N": EA,
                "mass_kg_per_m": {m: mass(m) for m in WEIGHT},
                "submerged_weight_N_per_m": WEIGHT,
                "rayleigh_damping": {
                    "zeta": ZETA,
                    "f1_hz": RAYLEIGH[0],
                    "f2_hz": RAYLEIGH[1],
                },
                "floor_friction": FRICTION,
                "segments": SEGMENTS,
                "dt_s": DT,
                "tmax_s": TMAX,
                "statistics_from_s": T0,
                "current_ramp_s": RAMP,
                "welch_segment_s": WELCH,
                "probe_elevations_m": PROBES,
            },
            "still_water_top_tension_N": {
                m: {
                    "model": still[m],
                    "measured_fig10": meas.get(("10", m), {}).get(0.0, [(None,)])[0][0],
                    "measured_fig13": meas.get(("13", m), {}).get(0.0, [(None,)])[0][0],
                }
                for m in still
            },
            "scores": summary_scores,
            "sensitivity_and_convergence": changes,
            "executable_sha256": sha,
        }
    )
    if a.summary:
        out = Path(a.summary)
        out.mkdir(parents=True, exist_ok=True)
        fields = list(dict.fromkeys(k for r in rows for k in r))

        def cell(v):
            if v is None:
                return ""
            if isinstance(v, float):
                return f"{v:.6g}"
            return v

        with open(
            out / "response_comparison.csv", "w", newline="", encoding="utf-8"
        ) as fh:
            w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
            w.writeheader()
            for r in rows:
                w.writerow({k: cell(r.get(k)) for k in fields})
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n"
        )


if __name__ == "__main__":
    main()
