#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""VIV of the FLOTANT OD 48 mm research cable towed at 200 kg, and its forced-oscillation RAOs.

The tests of de Wilde, Heijnsdijk & de Ridder (2021, FLOTANT Deliverable D5.5, MARIN): a 6.4 m
OD 48 mm cable hanging from the towing carriage with a 200 kg clump weight, towed at 0.025-1.0
m/s, and the same cable forced by a hexapod in air and in water. Writes a deck per case, runs
them in parallel and compares:
  tow  cross-flow and in-line A/D_rms (sigma/OD) at 1/2L, 3/4L and L against the measured
       values of Appendix 3, for three published VIV coefficient sets;
  rao  the RAOs of the stations in air (Appendix 1) and in water (Appendix 2), the structural
       check of the model: in air the cable section, in water three models of the lower end.

    python flotant_viv.py path/to/cabledyn.exe WORKDIR [--parts tow,rao] [--sets uls19,b1,wake]
        [--speeds all | 0.3,0.475,...] [--end identified|table|free] [--dt 0.002] [--refine 1]
        [--tmax 60] [--t0 20] [--jobs N] [--keep] [--summary DIR]

Sources:
  J. de Wilde, R. Heijnsdijk, E.-J. de Ridder (2021), D5.5 - Report on VIV (hydrodynamic)
  behaviour, FLOTANT deliverable 210430-FLT-WP5_D-5-5_V0 (MARIN report 30995-53-DWB),
  https://flotantproject.eu/wp-content/uploads/2021/11/210430_FLT-D5.5_VIV_test_v0.pdf
  (the set-up, Sec. 2; processing, Sec. 3.4; FEM input, Secs. 4.5 and 5.4; Tables 6-1 and 6-2;
  Appendices 1-3)
  J.V. Ulveseter, M.J. Thorsen, S. Saevik, C.M. Larsen (2019), Simulating fundamental and higher
  harmonic VIV of slender structures, Applied Ocean Research 90, 101856,
  https://doi.org/10.1016/j.apor.2019.101856 (Table 1, Sec. 4.1)
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
# Fresh water: the report does not state the density; Table 2-2's 2.8 kg/m and 9.7 N/m in
# water at OD 48 mm imply 1000 kg/m^3.
RHO = {"air": 1.225, "water": 1000.0}
OD = 0.048  # m (Table 2-1)
F1_UR = 1.77  # Hz, the frequency of Ur in Appendix 3 (Sec. 6.3)
W_HIGHPASS = (
    3.5  # rad/s, the tow-test filter of the double-integrated accelerations (Table 3-4)
)

# The structure. The line follows the report's FEM input tables (fem_air.csv, fem_water.csv)
# node for node; the cable section and the lower end in water are identified from the RAOs
# (see the README).
EA_TABLE = 1.0e8  # N, the tables' EA
EI_STEEL = 1.0e6  # N m^2, steel links (rigid next to the cable)
EI_HINGE = 1.0  # N m^2, the Cardan joints and the steel wire
CABLE_EI = 1650.0  # N m^2, identified from the in-air RAOs (the tables give 3000)
CABLE_BB = 130.0  # N m^2 s, bending damping, identified from the in-air RAOs
K_AX = (
    2.0e4  # N/m, the counterweight (air) or clump (water) as an axial spring in series
)
# The identified lower end in water: the link chain and load cell resist bending (CHAIN_EI) and
# the clump, a point mass on a short rigid segment, is held laterally by a spring END_K and
# damper END_C, modelled as a clamped cantilever post of length POST_LEN (K = 3 EI_p / l^3).
CHAIN_EI = 100.0  # N m^2
END_K = 32500.0  # N/m
END_C = 20.0  # N s/m
POST_LEN = 0.30  # m
CLUMP_SEG = 0.10  # m
POST_TENSION = (
    1000.0  # N, the static tension the post and segment lengths are stretched for
)
CLUMP_MASS = 200.0  # kg
CLUMP_WEIGHT = 1710.5  # N in water (Sec. 5.4 table, the clump node)
CLUMP_CDA = 0.05  # m^2, assumed (not given)
CD, CA = 1.2, 1.0  # the cable and the hardware
ENDS = ("identified", "table", "free")
SPEEDS = [round(0.025 * k, 3) for k in range(1, 41)]  # Table 6-2 [m/s]
STATIONS = {  # table node of each accelerometer -> station
    "water": {46: "1_2L", 67: "3_4L", 88: "L"},
    "air": {85: "1_8L", 75: "1_4L", 54: "1_2L", 33: "3_4L"},
}
# The forced-oscillation cases compared: the measured frequencies of Appendices 1 and 2 at
# 5 and 10 mm in a band about the first modes. The line is axisymmetric, so one direction (X)
# of forcing stands for both; it is compared with the surge and the sway tests.
RAO_CASES = {
    ("air", 5): [round(1.0 + 0.1 * k, 1) for k in range(16)],
    ("air", 10): [round(1.0 + 0.1 * k, 1) for k in range(11)],
    ("water", 5): [round(1.0 + 0.1 * k, 1) for k in range(21)],
}
RAO_SCORE_BAND = {"air": (1.4, 2.0), "water": (1.0, 3.0)}  # Hz
RAO_WINDOW = 10.0  # s, the harmonic fit at the end of each forced run
UR_MIN = 3.0  # the compared tow tests: from the onset the report states (Sec. 6.3)
RELEASE_FRACTION = 0.7  # release: the first test above the peak below 70 % of the peak
ANALYSIS = "1"  # part of the reuse key: bump when an analysis changes


def chaplin_amplification(ad: float) -> float:
    """Chaplin et al. (2005, Eq. 4): Cd/Cd0 = 1 + 2.3 s + 1.7 s^2, s = (A/D)/sqrt(2)."""
    s = ad / math.sqrt(2.0)
    return 1.0 + 2.3 * s + 1.7 * s * s


# Coefficient sets on the cable: LINE TYPES Cl dF cF; VIV row (ClX, alpha, Chh); VIV AMPLITUDE
# with Chaplin's drag amplification; VIV DAMPING (a, b, c); WAKE OSCILLATOR (CL0, St, eps, A).
SETS: dict[str, dict] = {
    "uls19": {
        "label": "Ulveseter et al. (2019), Table 1, Chh 1.5",
        "fitted": "the NDP riser in uniform flow; blind here",
        "cl": (0.85, 0.064, 0.144),
        "viv": (0.75, 0.15, 1.5),
        "amp": False,
        "damping": None,
        "wake": None,
    },
    "b1": {
        "label": "generic set (Delta Flume case)",
        "fitted": "no; blind",
        "cl": (0.8, 0.08, 0.18),
        "viv": (0.85, 0.15, None),
        "amp": True,
        "damping": (0, 0, 1),
        "wake": None,
    },
    "wake": {
        "label": "Facchinetti et al. (2004) wake oscillator",
        "fitted": "rigid-cylinder data; blind",
        "cl": (0, 0.08, 0.18),
        "viv": None,
        "amp": False,
        "damping": None,
        "wake": (0.3, 0.2, 0.3, 12),
    },
}


def load_table(medium: str) -> list[dict]:
    with open(HERE / f"fem_{medium}.csv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ("pos", "x", "y", "z", "mass", "weight", "tension", "EA", "EI"):
            r[k] = float(r[k]) if r[k] not in ("", "None") else None
    return rows


def clump_volume() -> float:
    return (CLUMP_MASS * GRAVITY - CLUMP_WEIGHT) / (RHO["water"] * GRAVITY)


def line_runs(medium: str, end: str) -> tuple[list[dict], list, float]:
    """The line as runs of equal table elements, End A (top) first.

    Each table element (0.0765 m) is one element of the run labelled by its lower table node
    (Wire, Steel, Joint, Cable, Chain, Sensor), with the node's lumped mass over the element
    and, in water, the diameter of the displaced volume of the table's weight (the cable: OD).
    The unstretched lengths give the table's static tensions with the effective axial stiffness
    EA_eff = L / (1/k_ax + L/EA). Returns (rows, runs, EA_eff); a run is
    [(label, mass per m, diameter, EI), length, elements].
    """
    rows = load_table(medium)
    ltot = abs(rows[-1]["z"] - rows[0]["z"])
    ea_eff = ltot / (1.0 / K_AX + ltot / EA_TABLE)
    elems = []
    for k in range(1, len(rows)):
        src, lab = rows[k], rows[k]["element"]
        if k == len(rows) - 1:
            # the element into the end point (the hexapod joint in air, the clump in water):
            # the end node's lumped mass is not line mass
            prev = rows[k - 1]
            if lab in ("Clump", "Hex", "Top"):
                lab = prev["element"]
            src = dict(prev)
            if lab == "Joint":
                src["mass"] = 0.1
                src["weight"] = 0.1 * GRAVITY * (1.0 if medium == "air" else 0.87)
        dz = abs(rows[k]["z"] - rows[k - 1]["z"])
        l0 = dz / (1.0 + src["tension"] / ea_eff)
        m = src["mass"] / dz
        if lab == "Cable":
            diam = OD
        elif medium == "water":
            vol = max(
                (m * GRAVITY - src["weight"] / dz) / (RHO["water"] * GRAVITY), 0.0
            )
            diam = math.sqrt(4.0 * vol / math.pi) if vol else 0.05
        else:
            diam = 0.05
        ei = (
            CABLE_EI
            if lab == "Cable"
            else (EI_HINGE if lab in ("Joint", "Wire") else EI_STEEL)
        )
        if end == "identified" and medium == "water" and lab in ("Chain", "Sensor"):
            ei = CHAIN_EI
        elems.append((lab, l0, m, diam, ei))
    if end == "identified" and medium == "water":
        # the clump as a point load on a light rigid segment, then the post to a still anchor
        stretch = 1.0 + POST_TENSION / ea_eff
        elems += [("Clumpseg", CLUMP_SEG / 2.0 / stretch, 1.0, 0.05, EI_STEEL)] * 2
        elems += [
            ("Post", POST_LEN / 4.0 / stretch, 1.0, 0.02, END_K * POST_LEN**3 / 3.0)
        ] * 4
    # consecutive elements of one label form a run; mass and diameter are length-weighted means
    # (the table's 0.076/0.077 m spacing would otherwise split every run)
    groups: list[list] = []
    for e in elems:
        if groups and groups[-1][0][0] == e[0]:
            groups[-1].append(e)
        else:
            groups.append([e])
    runs = []
    for g in groups:
        lg = sum(e[1] for e in g)
        m = sum(e[2] * e[1] for e in g) / lg
        dsq = sum(e[3] ** 2 * e[1] for e in g) / lg  # keeps the displaced volume
        runs.append([(g[0][0], m, math.sqrt(dsq), g[0][4]), lg, len(g)])
    return rows, runs, ea_eff


def station_nodes(medium: str, refine: int) -> dict[str, int]:
    """CableDyn node of each station: table node k is node k * refine + 1."""
    return {st: k * refine + 1 for k, st in STATIONS[medium].items()}


def deck(
    medium: str,
    end: str,
    hset: str | None,
    current: float,
    motion: tuple[float, float] | None,
    dt: float,
    tmax: float,
    refine: int,
) -> dict[str, str]:
    """The deck and its input files: {file name: text}."""
    rows, runs, ea_eff = line_runs(medium, end)
    s = SETS[hset] if hset else None
    types, secs, cable = [], [], []
    for i, ((lab, m, diam, ei), l0, n) in enumerate(runs):
        tname = f"{lab.lower()}{i}"
        cl = "0 0.08 0.18"
        cd, ca = CD, CA
        if lab == "Post":
            bb = END_C * POST_LEN**3 / 3.0
            eicol = f"{ei:.9g}|{bb:.9g}"
            cd, ca = 0.0, 0.0
        elif lab == "Cable":
            eicol = f"{CABLE_EI:g}|{CABLE_BB:g}"
            if s:
                cl = " ".join(f"{v:g}" for v in s["cl"])
            cable.append(tname)
        else:
            eicol = f"{ei:g}"
        types.append(
            f"{tname} {diam:.6f} {m:.6f} {ea_eff:.6e} 0.0 {eicol} {cd} 0.0 {ca} 0.0 {cl}"
        )
        secs.append(f"1 {tname} {l0:.9f} {n * refine}")
    z_top, z_bot = rows[0]["z"], rows[-1]["z"]
    extra = ""
    if end == "identified" and medium == "water":
        arc = 0.0
        for (lab, *_), l0, _n in runs:
            if lab == "Clumpseg":
                arc += 0.5 * l0
                break
            arc += l0
        extra += (
            "--------------------- ATTACHMENTS ----------------------------------\n"
            "LineID  ArcLength  Mass  Volume  CdA  Ca\n(-)  (m)  (kg)  (m^3)  (m^2)  (-)\n"
            f"1  {arc:.9f}  {CLUMP_MASS}  {clump_volume():.6f}  {CLUMP_CDA}  1.0\n"
            "--------------------- END CONNECTIONS -------------------------------\n"
            "LineID  End  Stiffness  EzX  EzY  EzZ\n(-)     (-)  (N-m/rad)  (-)  (-)  (-)\n"
            "1       B    Rigid      0    0    -1\n"
        )
        z_bot -= CLUMP_SEG + POST_LEN
    if s and s["viv"]:
        clx, alpha, chh = s["viv"]
        cols, vals = "ClX alpha", f"{clx} {alpha}"
        if chh is not None:
            cols, vals = cols + " Chh", vals + f" {chh}"
        extra += (
            "--------------------- VIV ----------------------------------------------\n"
            f"TypeName {cols}\n(-){' (-)' * len(cols.split())}\n"
            + "".join(f"{c} {vals}\n" for c in cable)
        )
    if s and s["amp"]:
        extra += (
            "--------------------- VIV AMPLITUDE ------------------------------------\n"
            "TypeName AD Cv CdAmp\n(-) (-) (-) (-)\n"
            + "".join(
                f"{c} {0.1 * k:.2f} {s['cl'][0]} {chaplin_amplification(0.1 * k):.5f}\n"
                for c in cable
                for k in range(16)
            )
        )
    if s and s["damping"]:
        a, b, c3 = s["damping"]
        extra += (
            "--------------------- VIV DAMPING ------------------------------------\n"
            "TypeName a b c\n(-) (-) (-) (-)\n"
            + "".join(f"{c} {a} {b} {c3}\n" for c in cable)
        )
    if s and s["wake"]:
        extra += (
            "--------------------- WAKE OSCILLATOR ----------------------------------\n"
            "TypeName CL0 St eps A\n(-) (-) (-) (-) (-)\n"
            + "".join(f"{c} {' '.join(f'{v:g}' for v in s['wake'])}\n" for c in cable)
        )
    bottom = f"2 Coupled 0.0 0.0 {z_bot:.6f} 0 0 0 0"
    if end == "free" and medium == "water":
        bottom = f"2 Free 0.0 0.0 {z_bot:.6f} {CLUMP_MASS} {clump_volume():.6f} {CLUMP_CDA} 1.0"
    nodes = station_nodes(medium, refine)
    outs = ["FairTen1"] + [
        f"L1N{j + d}p{c}" for j in nodes.values() for d in (-1, 0, 1) for c in "xyz"
    ]
    files: dict[str, str] = {}
    if current:
        files["current_profile.txt"] = (
            f"z ux uy uz\n(m) (m/s) (m/s) (m/s)\n{z_bot - 2.0:.3f} {current} 0 0\n0.0 {current} 0 0\n"
        )
    title = f"FLOTANT D5.5 OD 48 mm cable, {medium}, end {end}" + (
        f", tow {current} m/s, set {hset}" if current else f", forced {motion}"
    )
    files[
        "flotant.dat"
    ] = f"""--------------------- CableDyn Input File ---------------------------------
{title}
----------------------- LINE TYPES ------------------------------------------
TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t Cl dF cF
(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-) (-) (-) (-)
{chr(10).join(types)}
----------------------- POINTS ----------------------------------------------
ID Attachment X Y Z M V CdA CA
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Coupled 0.0 0.0 {z_top:.6f} 0 0 0 0
{bottom}
----------------------- LINES -----------------------------------------------
ID NodeA NodeB Outputs
(-) (-) (-) (-)
1 1 2 -
----------------------- SECTIONS --------------------------------------------
LineID LineType Length NumSegs
(-) (-) (m) (-)
{chr(10).join(secs)}
{extra}---------------------- OPTIONS ----------------------------------------------
{RHO[medium]} rhoW
{GRAVITY} g
{abs(z_bot) + 2.0:.3f} WtrDpth
{dt:.6f} dtM
{tmax:.6f} TMax
{"1 Currents" + chr(10) if current else ""}{"motion.txt" if motion else "0"} motionFile
---------------------- OUTPUTS -----------------------------------------------
{chr(10).join(outs)}
END
------------------------- need this line --------------------------------------
"""
    if motion:
        files["motion.txt"] = motion_text(medium, end, z_top, z_bot, motion, dt, tmax)
    return files


def motion_text(medium, end, z_top, z_bot, motion, dt, tmax) -> str:
    """The hexapod (the bottom point in air, the top point in water) moves in x as A sin(w t),
    reached over a 5-period cosine ramp; the other Coupled point is still."""
    amp, f = motion
    hexapod = 2 if medium == "air" else 1
    points = (
        ((1, z_top),)
        if (end == "free" and medium == "water")
        else ((1, z_top), (2, z_bot))
    )
    w, tr = 2.0 * math.pi * f, 5.0 / f
    lines = ["# time point x y z vx vy vz ax ay az"]
    for k in range(round(tmax / dt) + 1):
        t = k * dt
        if t < tr:
            r = 0.5 * (1.0 - math.cos(math.pi * t / tr))
            dr = 0.5 * math.pi / tr * math.sin(math.pi * t / tr)
            ddr = 0.5 * (math.pi / tr) ** 2 * math.cos(math.pi * t / tr)
        else:
            r, dr, ddr = 1.0, 0.0, 0.0
        sn, cs = math.sin(w * t), math.cos(w * t)
        for pid, z0 in points:
            x = vx = ax = 0.0
            if pid == hexapod:
                x = amp * r * sn
                vx = amp * (dr * sn + r * w * cs)
                ax = amp * (ddr * sn + 2.0 * dr * w * cs - r * w * w * sn)
            lines.append(
                f"{t:.6f} {pid} {x:.9e} {0.0:.9e} {z0:.9e} {vx:.9e} {0.0:.9e} "
                f"{0.0:.9e} {ax:.9e} {0.0:.9e} {0.0:.9e}"
            )
    return "\n".join(lines) + "\n"


def read_output(out: Path) -> tuple[list[str], np.ndarray]:
    with open(out, encoding="utf-8") as fh:
        fh.readline()
        head = fh.readline().split()
        fh.readline()
        a = np.loadtxt(fh)
    return head, a


def highpass(x: np.ndarray, theta: np.ndarray | None, dt: float) -> np.ndarray:
    """The displacement as the report derives it from an accelerometer that tilts with the
    cable: x + g theta / omega^2 (the gravity component on the tilted sensor axis, with theta
    the slope dx/dz; theta None leaves it out), high-passed at W_HIGHPASS rad/s."""
    f = np.fft.rfftfreq(len(x), dt)
    w = 2.0 * np.pi * f
    sp = np.fft.rfft(x - x.mean())
    keep = w >= W_HIGHPASS
    if theta is not None:
        sp[keep] += GRAVITY * np.fft.rfft(theta - theta.mean())[keep] / w[keep] ** 2
    sp[~keep] = 0.0
    return np.fft.irfft(sp, len(x))


def slope(a: np.ndarray, head: list[str], j: int, c: str) -> np.ndarray:
    """dx/dz at node j from its neighbours (central difference)."""
    xm, xp = a[:, head.index(f"L1N{j - 1}p{c}")], a[:, head.index(f"L1N{j + 1}p{c}")]
    zm, zp = a[:, head.index(f"L1N{j - 1}pz")], a[:, head.index(f"L1N{j + 1}pz")]
    return (xm - xp) / (zm - zp)


def analyse_tow(out: Path, u: float, t0: float, refine: int) -> dict:
    """A/D_rms = sigma / OD of each station and direction over t > t0, as measured (with the
    tilt term) and of the displacement itself; the dominant cross-flow frequency at 1/2L."""
    head, a = read_output(out)
    k = a[:, 0] > t0
    a = a[k]
    dt = float(a[1, 0] - a[0, 0])
    res: dict = {"mean_top_tension_N": float(a[:, head.index("FairTen1")].mean())}
    for st, j in station_nodes("water", refine).items():
        for comp, c in (("IL", "x"), ("CF", "y")):
            x = a[:, head.index(f"L1N{j}p{c}")]
            res[f"AD_{comp}_{st}"] = float(
                highpass(x, slope(a, head, j, c), dt).std() / OD
            )
            res[f"AD_{comp}_{st}_displacement"] = float(
                highpass(x, None, dt).std() / OD
            )
    j = station_nodes("water", refine)["1_2L"]
    y = highpass(a[:, head.index(f"L1N{j}py")], None, dt)
    f = np.fft.rfftfreq(len(y), dt)
    fp = float(f[int(np.argmax(np.abs(np.fft.rfft(y))))])
    res["f_cf_1_2L_hz"] = fp
    res["fD_over_U_cf_1_2L"] = fp * OD / u
    res["frequency_resolution_hz"] = float(1.0 / (len(y) * dt))
    return res


def harmonic(t: np.ndarray, x: np.ndarray, w: float) -> complex:
    m = np.column_stack([np.sin(w * t), np.cos(w * t), np.ones_like(t)])
    c, *_ = np.linalg.lstsq(m, x, rcond=None)
    return complex(c[0], c[1])


def analyse_rao(out: Path, medium: str, amp: float, f: float, refine: int) -> dict:
    """RAO of each station along the forcing over the last RAO_WINDOW s: the harmonic of
    x + g theta / omega^2 (as the accelerometers measure it) over the hexapod amplitude."""
    head, a = read_output(out)
    t = a[:, 0]
    k = t > t[-1] - RAO_WINDOW
    w = 2.0 * math.pi * f
    res = {}
    for st, j in station_nodes(medium, refine).items():
        x = a[k, head.index(f"L1N{j}px")]
        th = slope(a[k], head, j, "x")
        res[f"RAO_{st}"] = abs(harmonic(t[k], x + GRAVITY * th / w**2, w)) / amp
        res[f"RAO_{st}_displacement"] = abs(harmonic(t[k], x, w)) / amp
    return res


def run_case(
    exe: Path, folder: Path, files: dict, analyse, key: str, keep: bool, sha: str
) -> dict:
    """Write the deck, run it and analyse it. A folder that holds the same input files run by
    the same executable is not run again: its result is reused when the analysis key is also
    the same, and its output is analysed again if it was kept."""
    folder.mkdir(parents=True, exist_ok=True)
    result, out, ran = (
        folder / "result.json",
        folder / "flotant.out",
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
        for n, v in files.items():
            (folder / n).write_text(v, encoding="utf-8", newline="\n")
        proc = subprocess.run(
            [str(exe), "flotant.dat", "flotant"],
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
    res = analyse(out)
    result.write_text(
        json.dumps({**res, "reuse_key": f"{sha}:{key}"}, indent=1) + "\n",
        encoding="utf-8",
    )
    if not keep:
        out.unlink()
    return res


def read_csv(name: str) -> list[dict]:
    with open(HERE / name, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def rms_log(ratios: list[float]) -> float | None:
    return float(np.sqrt(np.mean(np.log(ratios) ** 2))) if ratios else None


def peak_and_release(points: list[tuple[float, float]]) -> dict:
    """The largest A/D, its Ur, and the release: the first Ur above the peak at which A/D
    falls below RELEASE_FRACTION of the peak (None when it never does)."""
    points = sorted(points)
    i = max(range(len(points)), key=lambda k: points[k][1])
    ur_p, ad_p = points[i]
    rel = next((ur for ur, ad in points[i + 1 :] if ad < RELEASE_FRACTION * ad_p), None)
    return {"peak_AD": ad_p, "peak_Ur": ur_p, "release_Ur": rel}


def parse_speeds(text: str) -> list[float]:
    if text == "all":
        return SPEEDS
    out = sorted({round(float(v), 3) for v in text.split(",")})
    if not out or any(v not in SPEEDS for v in out):
        raise SystemExit(
            f"--speeds: the tests are 0.025, 0.050, ..., 1.000 m/s, got {text!r}"
        )
    return out


def compare_tow(results: list[dict], sets: list[str]) -> tuple[list[dict], dict]:
    meas: dict = {}
    for r in read_csv("measured_tow_ad.csv"):
        comp = "CF" if r["component"] == "Y" else "IL"
        meas[(r["station"], comp, float(r["V_test"]))] = (
            float(r["AD_rms"]),
            float(r["Ur"]),
            r["direction"],
        )
    rows, per_set = [], {}
    stations = list(STATIONS["water"].values())
    for s in sets:
        rs = sorted((r for r in results if r["set"] == s), key=lambda r: r["U_m_per_s"])
        ratios: dict = {}
        for r in rs:
            u = r["U_m_per_s"]
            ur = u / (OD * F1_UR)
            row = {"set": s, "U_m_per_s": u, "Ur": round(ur, 3)}
            for comp in ("CF", "IL"):
                for st in stations:
                    mv = meas.get((st, comp, u))
                    model = r[f"AD_{comp}_{st}"]
                    row[f"AD_{comp}_{st}_measured"] = mv[0] if mv else ""
                    row[f"AD_{comp}_{st}_cabledyn"] = round(model, 4)
                    row[f"AD_{comp}_{st}_cabledyn_displacement"] = round(
                        r[f"AD_{comp}_{st}_displacement"], 4
                    )
                    if mv:
                        row["direction"] = mv[2]
                    if mv and ur >= UR_MIN:
                        ratios.setdefault((comp, st), []).append(model / mv[0])
                        ratios.setdefault((comp, st, mv[2]), []).append(model / mv[0])
            row["f_cf_1_2L_hz"] = round(r["f_cf_1_2L_hz"], 3)
            row["fD_U_cf_1_2L"] = round(r["fD_over_U_cf_1_2L"], 4)
            row["mean_top_tension_N"] = round(r["mean_top_tension_N"], 1)
            rows.append(row)
        model_pts = [(r["U_m_per_s"] / (OD * F1_UR), r["AD_CF_1_2L"]) for r in rs]
        meas_pts = [
            (mv[1], mv[0])
            for (st, comp, u), mv in meas.items()
            if st == "1_2L" and comp == "CF" and u in [r["U_m_per_s"] for r in rs]
        ]
        pk = peak_and_release(model_pts)
        f_at_peak = next(
            r["f_cf_1_2L_hz"]
            for r in rs
            if abs(r["U_m_per_s"] / (OD * F1_UR) - pk["peak_Ur"]) < 1e-9
        )
        cf_all = [q for st in stations for q in ratios.get(("CF", st), [])]
        il_all = [q for st in stations for q in ratios.get(("IL", st), [])]
        per_set[s] = {
            "label": SETS[s]["label"],
            "fitted_to": SETS[s]["fitted"],
            "speeds": [r["U_m_per_s"] for r in rs],
            "compared_tests": f"Ur >= {UR_MIN:g}",
            "cf_rms_log_error_all_stations": rms_log(cf_all),
            "cf_ratio_median_all_stations": float(np.median(cf_all))
            if cf_all
            else None,
            "il_rms_log_error_all_stations": rms_log(il_all),
            "il_ratio_median_all_stations": float(np.median(il_all))
            if il_all
            else None,
            "cf_rms_log_error": {
                st: rms_log(ratios.get(("CF", st), [])) for st in stations
            },
            "il_rms_log_error": {
                st: rms_log(ratios.get(("IL", st), [])) for st in stations
            },
            "cf_rms_log_error_1_2L_by_direction": {
                d: rms_log(ratios.get(("CF", "1_2L", d), []))
                for d in ("forwards", "backwards")
            },
            "cf_1_2L_peak": {**pk, "f_hz_at_peak": f_at_peak},
            "cf_1_2L_peak_measured": peak_and_release(meas_pts) if meas_pts else None,
            "mean_top_tension_N_range": [
                round(min(r["mean_top_tension_N"] for r in rs)),
                round(max(r["mean_top_tension_N"] for r in rs)),
            ],
        }
        p, q = per_set[s]["cf_1_2L_peak"], per_set[s]["cf_1_2L_peak_measured"]
        print(
            f"{s}: CF r.m.s. log error {per_set[s]['cf_rms_log_error_all_stations']:.3f} "
            f"(1/2L {per_set[s]['cf_rms_log_error']['1_2L']:.3f}), IL "
            f"{per_set[s]['il_rms_log_error_all_stations']:.3f}; 1/2L peak {p['peak_AD']:.3f} at Ur "
            f"{p['peak_Ur']:.2f} ({p['f_hz_at_peak']:.2f} Hz), release Ur {p['release_Ur']}"
            + (
                f" | measured {q['peak_AD']:.3f} at {q['peak_Ur']:.2f}, release {q['release_Ur']}"
                if q
                else ""
            )
        )
    return rows, per_set


def compare_rao(results: list[dict]) -> tuple[list[dict], dict]:
    meas: dict = {}
    for r in read_csv("measured_rao.csv"):
        meas[
            (r["medium"], r["motion"], r["station"], int(r["amp_mm"]), float(r["f_hz"]))
        ] = float(r["rao"])
    rows, summary = [], {}
    for r in sorted(
        results, key=lambda r: (r["medium"], r["end"], r["amp_mm"], r["f_hz"])
    ):
        for st in STATIONS[r["medium"]].values():
            row = {k: r[k] for k in ("medium", "end", "amp_mm", "f_hz")}
            row["station"] = st
            row["rao_cabledyn"] = round(r[f"RAO_{st}"], 4)
            row["rao_cabledyn_displacement"] = round(r[f"RAO_{st}_displacement"], 4)
            for motion, name in (("X", "surge"), ("Y", "sway")):
                row[f"rao_measured_{name}"] = meas.get(
                    (r["medium"], motion, st, r["amp_mm"], r["f_hz"]), ""
                )
            rows.append(row)
    groups = sorted({(r["medium"], r["end"], r["amp_mm"]) for r in rows})
    for medium, end, amp in groups:
        lo, hi = RAO_SCORE_BAND[medium]
        g = [
            r
            for r in rows
            if (r["medium"], r["end"], r["amp_mm"]) == (medium, end, amp)
        ]
        entry = {}
        for name in ("surge", "sway"):
            err = [
                math.log(r["rao_cabledyn"] / r[f"rao_measured_{name}"])
                for r in g
                if r[f"rao_measured_{name}"] != ""
                and lo - 1e-9 <= r["f_hz"] <= hi + 1e-9
            ]
            half = sorted(
                (r for r in g if r["station"] == "1_2L"), key=lambda r: r["f_hz"]
            )
            mh = [
                (r["f_hz"], r[f"rao_measured_{name}"])
                for r in half
                if r[f"rao_measured_{name}"] != ""
            ]
            entry[name] = {
                "rms_log_error": float(np.sqrt(np.mean(np.square(err))))
                if err
                else None,
                "band_hz": [lo, hi],
                "peaks_1_2L_cabledyn": local_peaks(
                    [(r["f_hz"], r["rao_cabledyn"]) for r in half]
                ),
                "peaks_1_2L_measured": local_peaks(mh),
            }
        summary[f"{medium}_{end}_{amp}mm"] = entry
        print(
            f"rao {medium:5s} {end:10s} {amp:2d} mm: r.m.s. log error surge {entry['surge']['rms_log_error']:.3f}, "
            f"sway {entry['sway']['rms_log_error']:.3f}; 1/2L peaks {entry['surge']['peaks_1_2L_cabledyn']} "
            f"(surge test {entry['surge']['peaks_1_2L_measured']}, sway {entry['sway']['peaks_1_2L_measured']})"
        )
    return rows, summary


def local_peaks(pts: list[tuple[float, float]], n: int = 2) -> list[list[float]]:
    """The n highest local maxima [f, RAO] of a sampled curve, in frequency order."""
    pts = sorted(pts)
    i = [
        k
        for k in range(1, len(pts) - 1)
        if pts[k][1] >= pts[k - 1][1] and pts[k][1] >= pts[k + 1][1]
    ]
    top = sorted(i, key=lambda k: -pts[k][1])[:n]
    return [[pts[k][0], round(pts[k][1], 3)] for k in sorted(top)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("exe")
    ap.add_argument("workdir")
    ap.add_argument("--parts", default="tow,rao", help="comma list of tow, rao")
    ap.add_argument(
        "--sets", default=",".join(SETS), help="comma list of " + ", ".join(SETS)
    )
    ap.add_argument(
        "--speeds", default="all", help="comma list of tow speeds [m/s], or all"
    )
    ap.add_argument(
        "--end",
        default="identified",
        choices=ENDS,
        help="the lower end of the tow model",
    )
    ap.add_argument("--dt", type=float, default=0.002)
    ap.add_argument("--refine", type=int, default=1, help="elements per table element")
    ap.add_argument("--tmax", type=float, default=60.0, help="tow run length [s]")
    ap.add_argument(
        "--t0", type=float, default=20.0, help="start of the tow statistics [s]"
    )
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--keep", action="store_true", help="keep the output tables")
    ap.add_argument(
        "--summary",
        default=None,
        help="directory for the comparison CSVs and summary JSON",
    )
    a = ap.parse_args()
    found = a.exe if Path(a.exe).is_file() else shutil.which(a.exe)
    if found is None:
        raise SystemExit(f"CableDyn executable not found: {a.exe}")
    exe = Path(found).resolve()
    sha = hashlib.sha256(exe.read_bytes()).hexdigest()
    parts = [p.strip() for p in a.parts.split(",") if p.strip()]
    if not parts or any(p not in ("tow", "rao") for p in parts):
        raise SystemExit(f"--parts: choose from tow, rao, got {a.parts!r}")
    sets = [s.strip() for s in a.sets.split(",") if s.strip()]
    if not sets or any(s not in SETS for s in sets):
        raise SystemExit(f"--sets: choose from {', '.join(SETS)}, got {a.sets!r}")
    speeds = parse_speeds(a.speeds)
    tag = f"dt{a.dt:g}_r{a.refine}"
    work = Path(a.workdir).resolve()
    jobs = []
    if "tow" in parts:
        tmax = round(a.tmax / a.dt) * a.dt
        for s in sets:
            for u in speeds:
                files = deck("water", a.end, s, u, None, a.dt, tmax, a.refine)
                folder = work / f"tow_{s}_{a.end}_v{u:.3f}_{tag}_t{a.tmax:g}"
                key = f"tow:{a.t0:g}:{ANALYSIS}"

                def an(out, u=u):
                    return analyse_tow(out, u, a.t0, a.refine)

                jobs.append(
                    ({"part": "tow", "set": s, "U_m_per_s": u}, folder, files, an, key)
                )
    if "rao" in parts:
        for (medium, amp), freqs in RAO_CASES.items():
            for end in ENDS if medium == "water" else ("table",):
                for f in freqs:
                    tmax = round((5.0 / f + 25.0) / a.dt) * a.dt
                    files = deck(
                        medium, end, None, 0.0, (amp / 1000.0, f), a.dt, tmax, a.refine
                    )
                    folder = work / f"rao_{medium}_{end}_a{amp}mm_f{f:.2f}_{tag}"

                    def an(out, medium=medium, amp=amp, f=f):
                        return analyse_rao(out, medium, amp / 1000.0, f, a.refine)

                    info = {
                        "part": "rao",
                        "medium": medium,
                        "end": end,
                        "amp_mm": amp,
                        "f_hz": f,
                    }
                    jobs.append(
                        (info, folder, files, an, f"rao:{RAO_WINDOW:g}:{ANALYSIS}")
                    )

    def one(job):
        info, folder, files, an, key = job
        return {**info, **run_case(exe, folder, files, an, key, a.keep, sha)}

    with ThreadPoolExecutor(max_workers=max(1, min(a.jobs, len(jobs)))) as pool:
        results = list(pool.map(one, jobs))
    summary: dict = {
        "dt_s": a.dt,
        "elements_per_table_element": a.refine,
        "tow": None,
        "rao": None,
        "executable_sha256": sha,
    }
    tow_rows = rao_rows = None
    if "tow" in parts:
        tow_rows, per_set = compare_tow(
            [r for r in results if r["part"] == "tow"], sets
        )
        summary["tow"] = {
            "end": a.end,
            "tmax_s": a.tmax,
            "statistics_after_s": a.t0,
            "highpass_rad_s": W_HIGHPASS,
            "Ur_frequency_hz": F1_UR,
            "release_fraction_of_peak": RELEASE_FRACTION,
            "per_set": per_set,
        }
    if "rao" in parts:
        rao_rows, summary["rao"] = compare_rao(
            [r for r in results if r["part"] == "rao"]
        )
    if a.summary:
        out = Path(a.summary)
        out.mkdir(parents=True, exist_ok=True)
        for name, rows in (
            ("response_comparison.csv", tow_rows),
            ("rao_comparison.csv", rao_rows),
        ):
            if rows:
                with open(out / name, "w", newline="", encoding="utf-8") as fh:
                    fields = list(dict.fromkeys(k for r in rows for k in r))
                    w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
                    w.writeheader()
                    w.writerows(rows)
        (out / "summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n"
        )


if __name__ == "__main__":
    main()
