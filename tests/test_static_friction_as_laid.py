#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Static seabed friction toward the as-laid position (OPTIONS staticsFriction asLaid, set in
every deck below; the default is none).

A grounded run on a seabed sloping down toward its anchor, with friction mu = 1 well above the
slope, must hold in statics where it was laid; with the default frictionless statics (none) the
run slides down the slope and folds back past its anchor (EI = 0), or is rejected as too long
for its span (finite EI). Every deck runs through the standalone
driver:

* ``slope_ei0`` (EI = 0 chain, 5 degree slope, 700 m line over a 450 m span): converges with
  the grounded run on its side of the anchor, and an independent nodal force balance (element
  tensions, lumped weight, seabed normal along the local normal) leaves at every grounded node
  a horizontal friction force within mu N, mobilised against the slope; its twin with
  ``none staticsFriction`` still folds past the anchor (the frictionless static of before);
* ``slope_ei0_rest``: the same deck marched for 2 s from the static state stays at rest (every
  output row reproduces the static end forces to 1e-5 of the fairlead tension): the friction
  anchors of the march continue the static springs, with no jump at t = 0;
* ``slope_fei`` (finite-EI cable on a 10 degree slope, too long for its span) converges, and
  marched for 1 s stays at rest likewise;
* ``flat_*`` (a grounded chain on a flat seabed): laid at the tension its run carries
  (``asLaidTension``), the friction is not mobilised and the end tensions equal the
  frictionless ones (FairTen to 1e-4, AnchTen to 0.5 %, see the check); laid slack (the
  default), friction holds the run back and the anchor tension falls by at most mu w times the
  grounded length;
* input checks: ``layAzimuth<ID>`` of a line that does not exist and an unknown policy are
  rejected.

usage: test_static_friction_as_laid.py --exe <cabledyn> --output <dir>
"""

from __future__ import annotations

import argparse
import math
import os
import subprocess
import sys
from pathlib import Path

G, RHO = 9.80665, 1025.0
CHAIN = (0.216, 286.56, 1.22976e9)  # diameter [m], mass [kg/m], EA [N]
REST_TOL = 1.0e-5
FAILURES: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILURES.append(msg)


def slope_grid(path: Path, depth0: float, slope_deg: float) -> None:
    """A seabed z = -(depth0 + tan(slope) x): planar in x, falling toward +x."""
    t = math.tan(math.radians(slope_deg))
    rows = []
    for x in range(-300, 1000, 100):
        for y in (-300.0, 300.0):
            rows.append("%.1f %.1f %.12g" % (x, y, depth0 + t * x))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")


def deck(title: str, ei: float, length: float, nseg: int, fair: tuple, anch: tuple, seabed: str,
         extra: str = "", tmax: float = 0.0, mu: float = 1.0) -> str:
    d, m, ea = CHAIN
    return (
        "--------------------- CableDyn Input File ---\n"
        f"{title}\n"
        "--------------------- LINE TYPES ---\n"
        "TypeName Diam MassDenInAir EA BA EI Cd_n Cd_t Ca_n Ca_t\n"
        "(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)\n"
        f"chain {d} {m} {ea} -1 {ei} 1.37 0.64 1 0\n"
        "--------------------- POINTS ---\n"
        "ID Type X Y Z Mass Vol CdA Ca\n(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)\n"
        f"1 Coupled {fair[0]} {fair[1]} {fair[2]} 0 0 0 0\n"
        f"2 Fixed {anch[0]} {anch[1]} {anch[2]} 0 0 0 0\n"
        "--------------------- LINES ---\n"
        "ID NodeA NodeB Outputs\n(-) (-) (-) (-)\n1 1 2 pt\n"
        "--------------------- SECTIONS ---\n"
        "LineID LineType Length NumSegs\n(-) (-) (m) (-)\n"
        f"1 chain {length} {nseg}\n"
        "--------------------- OPTIONS ---\n"
        f"{G} g\n{RHO} rhoW\n{seabed}\n100000.0 kBot\n1.0e4 cBot\n0.05 dtM\n{tmax} TMax\n"
        f"{mu} frictionMu\nasLaid staticsFriction\n{extra}"
        "--------------------- OUTPUTS ---\n"
        '"FairTen1"\n"AnchTen1"\n"Point1Fx"\n"Point1Fy"\n"Point1Fz"\n"Point2Fx"\n"Point2Fy"\n"Point2Fz"\n'
        "--------------------- need this line ---\n"
    )


def run(exe: str, cdir: Path, name: str, text: str) -> tuple:
    cdir.mkdir(parents=True, exist_ok=True)
    (cdir / (name + ".dat")).write_text(text, encoding="utf-8", newline="\n")
    env = dict(os.environ)
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    p = subprocess.run([exe, name + ".dat", name], cwd=cdir, env=env, capture_output=True, text=True,
                       timeout=600)
    return p.returncode, p.stdout + p.stderr


def table(path: Path) -> list:
    rows = [r for r in path.read_text(encoding="utf-8").splitlines() if r.strip() and not r.startswith("#")]
    names = rows[0].split()
    return [dict(zip(names, map(float, r.split()))) for r in rows[1:]]


def first_row(path: Path) -> list:
    rows = [r for r in path.read_text(encoding="utf-8").splitlines() if r.strip() and not r.startswith("#")]
    return [float(v) for v in rows[1].split()[1:]]


def grounded_balance(cdir: Path, name: str, length: float, nseg: int, mu: float, depth0: float,
                     slope_deg: float) -> tuple:
    """Largest friction utilisation |f|/(mu N) over grounded nodes and the total friction."""
    d, m, _ = CHAIN
    w = (m - RHO * math.pi / 4.0 * d * d) * G
    l0 = length / nseg
    p = first_row(cdir / (name + ".Line1.p.out"))
    pos = [p[3 * i:3 * i + 3] for i in range(len(p) // 3)]
    ten = first_row(cdir / (name + ".Line1.t.out"))
    gx = -math.tan(math.radians(slope_deg))
    nrm = [-gx / math.hypot(1.0, gx), 0.0, 1.0 / math.hypot(1.0, gx)]
    worst, total = 0.0, 0.0
    for i in range(1, len(pos) - 1):
        f = [0.0, 0.0, -w * l0]
        for e, j in ((i - 1, i - 1), (i, i + 1)):
            u = [pos[j][k] - pos[i][k] for k in range(3)]
            n = math.sqrt(sum(x * x for x in u))
            for k in range(3):
                f[k] += ten[e] * u[k] / n
        zf = -(depth0 + math.tan(math.radians(slope_deg)) * pos[i][0])
        normal = -f[2] / nrm[2]
        if pos[i][2] > zf + 0.5 or not normal > 1.0e-6 * w * l0:
            continue
        fr = math.hypot(f[0] + normal * nrm[0], f[1])
        total += fr
        worst = max(worst, fr / (mu * normal))
    past = max(p_[0] for p_ in pos) - pos[-1][0]  # beyond the anchor (+x) along the span
    return worst, total, past


def at_rest(cdir: Path, static_name: str, rest_name: str) -> float:
    s0 = table(cdir / (static_name + ".out"))[0]
    scale = max(abs(s0["FairTen1"]), 1.0)
    drift = 0.0
    for row in table(cdir / (rest_name + ".out")):
        drift = max(drift, max(abs(row[n] - s0[n]) for n in s0 if n != "Time(s)"))
    return drift / scale


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    root = Path(a.output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    d, m, ea = CHAIN
    w = (m - RHO * math.pi / 4.0 * d * d) * G

    # --- EI = 0 run on a 5 degree slope (the R3 reproducer) ---
    cdir = root / "slope_ei0"
    cdir.mkdir(parents=True, exist_ok=True)
    depth0, slope = 260.6301014133342, 5.0
    slope_grid(cdir / "bathy.txt", depth0, slope)
    sb = "bathy.txt bathymetryFile"
    base = dict(ei=0.0, length=700.0, nseg=70, fair=(0, 0, -10), anch=(450, 0, -300), seabed=sb)
    rc, log = run(a.exe, cdir, "laid", deck("slope_ei0", **base))
    check(rc == 0, "slope_ei0 converges (rc %d)" % rc)
    if rc == 0:
        worst, total, past = grounded_balance(cdir, "laid", 700.0, 70, 1.0, depth0, slope)
        check(past < 0.5, "slope_ei0: the grounded run stays on its side of the anchor (%.3f m past it)" % past)
        check(worst <= 1.0 + 1.0e-3, "slope_ei0: friction within mu N at every grounded node (max %.4f)" % worst)
        check(worst > 0.5 * math.tan(math.radians(slope)),
              "slope_ei0: friction holds the slope load (max utilisation %.4f)" % worst)
    rc, log = run(a.exe, cdir, "none", deck("slope_ei0 none", **base, extra="none staticsFriction\n"))
    check(rc == 0, "slope_ei0 with staticsFriction none converges (rc %d)" % rc)
    if rc == 0:
        _, total, past = grounded_balance(cdir, "none", 700.0, 70, 1.0, depth0, slope)
        check(past > 10.0, "slope_ei0 frictionless static still folds past the anchor (%.2f m)" % past)
    rc, log = run(a.exe, cdir, "rest", deck("slope_ei0 rest", **base, tmax=2.0))
    check(rc == 0, "slope_ei0 at-rest march runs (rc %d)" % rc)
    if rc == 0:
        drift = at_rest(cdir, "laid", "rest")
        check(drift <= REST_TOL, "slope_ei0 stays at rest from the static state (drift %.2e of FairTen)" % drift)

    # --- finite-EI cable on a 10 degree slope, too long for its span ---
    cdir = root / "slope_fei"
    cdir.mkdir(parents=True, exist_ok=True)
    depth0, slope = 220.65285868119076, 10.0
    slope_grid(cdir / "bathy.txt", depth0, slope)
    fei = dict(ei=1.0e4, length=650.0, nseg=140, fair=(0, 0, -10), anch=(450, 0, -300), seabed=sb)
    rc, log = run(a.exe, cdir, "laid", deck("slope_fei", **fei))
    check(rc == 0, "slope_fei converges (rc %d)%s" % (rc, "" if rc == 0 else ": " + log.strip()[-300:]))
    if rc == 0:
        rc, log = run(a.exe, cdir, "rest", deck("slope_fei rest", **fei, tmax=1.0))
        check(rc == 0, "slope_fei at-rest march runs (rc %d)" % rc)
        if rc == 0:
            drift = at_rest(cdir, "laid", "rest")
            check(drift <= REST_TOL, "slope_fei stays at rest from the static state (drift %.2e)" % drift)

    # --- flat seabed: friction not mobilised when laid at the run's tension ---
    cdir = root / "flat"
    flat = dict(ei=0.0, length=850.0, nseg=85, fair=(-58, 0, -14), anch=(-837.6, 0, -200),
                seabed="200.0 WtrDpth", mu=0.5)
    rc0, _ = run(a.exe, cdir, "none", deck("flat none", **flat, extra="none staticsFriction\n"))
    check(rc0 == 0, "flat frictionless static converges")
    if rc0 == 0:
        r0 = table(cdir / "none.out")[0]
        h = abs(r0["Point2Fx"])  # the grounded run carries the anchor's horizontal tension
        rc, _ = run(a.exe, cdir, "taut", deck("flat taut", **flat, extra="%.9g asLaidTension\n" % h))
        check(rc == 0, "flat laid at the run tension converges")
        if rc == 0:
            r1 = table(cdir / "taut.out")[0]
            # FairTen to 1e-4. AnchTen to 0.5 %: the anchor rests on the seabed surface and the
            # laid run one static penetration below it, so the run is drawn a fraction of a
            # millimetre toward the anchor and its friction springs add a small tension there
            # (the anchor-level effect Orcina's "OrcaFlex Friction Effects" note describes).
            for k, tol in (("FairTen1", 1.0e-4), ("AnchTen1", 5.0e-3)):
                rel = abs(r1[k] - r0[k]) / r0[k]
                check(rel < tol, "flat laid at the run tension: %s equals the frictionless %.6g N (rel %.1e)"
                      % (k, r0[k], rel))
        rc, _ = run(a.exe, cdir, "slack", deck("flat slack", **flat))
        check(rc == 0, "flat laid slack converges")
        if rc == 0:
            r2 = table(cdir / "slack.out")[0]
            drop = r0["AnchTen1"] - r2["AnchTen1"]
            check(0.0 < drop <= 0.5 * w * 850.0,
                  "flat laid slack: friction lowers AnchTen by %.4g N (within mu w L)" % drop)

    # --- input checks ---
    cdir = root / "inputs"
    rc, log = run(a.exe, cdir, "badid", deck("bad id", **flat, extra="30.0 layAzimuth9\n"))
    check(rc != 0 and "names no LINE" in log, "layAzimuth9 without a line 9 is rejected")
    rc, log = run(a.exe, cdir, "badpolicy", deck("bad policy", **flat, extra="restart staticsFriction\n"))
    check(rc != 0 and "asLaid or none" in log, "an unknown staticsFriction policy is rejected")
    rc, log = run(a.exe, cdir, "goodid", deck("good id", **flat, extra="180.0 layAzimuth1\n0 asLaidTension1\n"))
    check(rc == 0, "layAzimuth1 / asLaidTension1 rows are accepted")

    print("static friction as laid: %d failure(s)" % len(FAILURES))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
