#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Still-water natural frequencies of the Delta Flume riser of Chaplin et al. (2005).

Writes a CableDyn deck for each of the five test series (Table 2 top tensions), runs the modal
analysis about the static state, and compares the first seven frequencies with Table 2.

    python chaplin_natural_frequencies.py path/to/cabledyn.exe WORKDIR [--segments 64]
        [--ca 1.0] [--ea 5.5e6] [--summary DIR]

Source: J.R. Chaplin, P.W. Bearman, F.J. Huera Huarte, R.J. Pattenden (2005), Laboratory
measurements of vortex-induced vibrations of a vertical tension riser in a stepped current,
Journal of Fluids and Structures 21(1), 3-24, https://doi.org/10.1016/j.jfluidstructs.2005.04.010
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path

# Section 2 of the paper: length, diameter, bending stiffness, mass ratio including the internal
# water, submerged weight (12.1 N/m, reproduced by the mass ratio below), universal joints at
# both ends. The tank was filled with water, so the whole riser is submerged.
LENGTH = 13.12  # m
DIAMETER = 0.028  # m
EI = 29.9  # N m^2
MASS_RATIO = 3.0  # mass / displaced mass, internal water included
RHO_W = 1000.0  # kg/m^3, fresh water (assumed)
GRAVITY = 9.80665  # m/s^2
AREA = math.pi * DIAMETER**2 / 4.0
MASS = MASS_RATIO * RHO_W * AREA  # kg/m
SUBMERGED_WEIGHT = (MASS - RHO_W * AREA) * GRAVITY  # N/m
Z_BOTTOM = -14.0  # m; any depth that keeps the riser submerged
WATER_DEPTH = 20.0  # m

# Table 2: still-water top tension [N] and the first seven natural frequencies [Hz] of each series
TABLE2: dict[str, tuple[float, tuple[float, ...]]] = {
    "A": (389.0, (0.405, 0.873, 1.326, 1.804, 2.278, 2.797, 3.389)),
    "B": (813.0, (0.633, 1.329, 1.954, 2.625, 3.312, 4.028, 4.753)),
    "C": (1185.0, (0.782, 1.423, 2.394, 3.203, 4.024, 4.875, 5.704)),
    "D": (1546.0, (0.899, 1.815, 2.741, 3.614, 4.561, 5.525, 6.439)),
    "E": (1925.0, (1.037, 2.046, 3.084, 4.099, 5.127, 6.222, 7.298)),
}
N_MODES = 7


def unstretched_length(top_tension: float, ea: float) -> float:
    """Unstretched length that gives ``top_tension`` on the vertical riser of length LENGTH.

    The effective tension falls linearly with the submerged weight from the top joint down,
    and an element ds of the stretched riser has the unstretched length ds / (1 + T/EA).
    """
    t_bottom = top_tension - SUBMERGED_WEIGHT * LENGTH
    # exact integral of ds / (1 + T(s)/EA) with T linear in s
    return ea / SUBMERGED_WEIGHT * math.log((ea + top_tension) / (ea + t_bottom))


def deck_text(
    series: str, top_tension: float, segments: int, ca: float, ea: float
) -> str:
    l0 = unstretched_length(top_tension, ea)
    return f"""--------------------- CableDyn Input File ---------------------------------
Delta Flume riser of Chaplin et al. (2005), series {series}: still-water top tension {top_tension:g} N
----------------------- LINE TYPES ------------------------------------------
TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t
(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-)
riser {DIAMETER} {MASS:.6f} {ea:.6e} 0.0 {EI} 1.2 0.0 {ca} 0.0
----------------------- POINTS ----------------------------------------------
ID Attachment X Y Z M V CdA CA
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Coupled 0.0 0.0 {Z_BOTTOM + LENGTH:.6f} 0 0 0 0
2 Fixed 0.0 0.0 {Z_BOTTOM:.6f} 0 0 0 0
----------------------- LINES -----------------------------------------------
ID LineType AttachA AttachB UnstrLen NumSegs Outputs
(-) (-) (-) (-) (m) (-) (-)
1 riser 1 2 {l0:.9f} {segments} -
---------------------- OPTIONS ----------------------------------------------
{RHO_W} rhoW
{GRAVITY} g
{WATER_DEPTH} WtrDpth
{2 * (N_MODES + 1)} nModes
0.01 dtM
0.0 TMax
---------------------- OUTPUTS -----------------------------------------------
FairTen1
AnchTen1
END
------------------------- need this line --------------------------------------
"""


def frequencies(path: Path) -> list[float]:
    """The distinct frequencies of a straight line's .modes.out (each in both planes)."""
    lines = path.read_text(encoding="utf-8").splitlines()
    start = lines.index("# Natural frequencies") + 3
    values = []
    for line in lines[start:]:
        if not line.strip() or line.startswith("#"):
            break
        values.append(float(line.split()[2]))
    pairs = values[0::2]
    for a, b in zip(values[0::2], values[1::2]):
        if not math.isclose(a, b, rel_tol=1e-6):
            raise SystemExit(f"{path.name}: the two planes differ ({a} and {b} Hz)")
    return pairs


def top_tension(path: Path) -> float:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines()[2:]]
    return float(rows[0][1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("executable", type=Path)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--segments", type=int, default=64)
    parser.add_argument(
        "--ca", type=float, default=1.0, help="normal added-mass coefficient"
    )
    parser.add_argument(
        "--ea", type=float, default=5.5e6, help="axial stiffness [N] (not given)"
    )
    parser.add_argument(
        "--summary", type=Path, help="write the committed summary files here"
    )
    args = parser.parse_args()

    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    rows = []
    for series, (tension, measured) in TABLE2.items():
        stem = f"chaplin_{series}"
        (work / f"{stem}.dat").write_text(
            deck_text(series, tension, args.segments, args.ca, args.ea),
            encoding="utf-8",
        )
        run = subprocess.run(
            [str(args.executable.resolve()), f"{stem}.dat", stem],
            cwd=work,
            env=env,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
        (work / f"{stem}.console.txt").write_text(
            run.stdout + run.stderr, encoding="utf-8"
        )
        if run.returncode != 0:
            raise SystemExit(f"{stem}: CableDyn exited with {run.returncode}")
        computed = frequencies(work / f"{stem}.modes.out")[:N_MODES]
        t_top = top_tension(work / f"{stem}.out")
        for mode, (f_c, f_m) in enumerate(zip(computed, measured), start=1):
            rows.append(
                {
                    "series": series,
                    "top_tension_N": tension,
                    "computed_top_tension_N": round(t_top, 3),
                    "mode": mode,
                    "measured_Hz": f_m,
                    "cabledyn_Hz": round(f_c, 5),
                    "error_percent": round(100.0 * (f_c / f_m - 1.0), 2),
                }
            )

    errors = [abs(float(r["error_percent"])) for r in rows]
    # Series C mode 2 (Table 2: 1.423 Hz) breaks the near-harmonic pattern of every other series
    # (f2/f1 of 1.82 against 1.97-2.16): reported, and also left out of a second statistic
    regular = [
        abs(float(r["error_percent"]))
        for r in rows
        if (r["series"], r["mode"]) != ("C", 2)
    ]
    signed = [
        float(r["error_percent"]) for r in rows if (r["series"], r["mode"]) != ("C", 2)
    ]
    summary = {
        "source": "Chaplin et al. (2005), J. Fluids Struct. 21, 3-24, Table 2",
        "doi": "10.1016/j.jfluidstructs.2005.04.010",
        "segments": args.segments,
        "added_mass_coefficient": args.ca,
        "axial_stiffness_N": args.ea,
        "mass_kg_per_m": round(MASS, 5),
        "submerged_weight_N_per_m": round(SUBMERGED_WEIGHT, 3),
        "frequencies_compared": len(rows),
        "mean_abs_error_percent": round(sum(errors) / len(errors), 2),
        "max_abs_error_percent": round(max(errors), 2),
        "without_C2": {
            "frequencies_compared": len(regular),
            "mean_error_percent": round(sum(signed) / len(signed), 2),
            "mean_abs_error_percent": round(sum(regular) / len(regular), 2),
            "max_abs_error_percent": round(max(regular), 2),
        },
        "executable_sha256": hashlib.sha256(args.executable.read_bytes()).hexdigest(),
    }
    for row in rows:
        print(
            f"{row['series']} mode {row['mode']}: measured {row['measured_Hz']:.3f} Hz, "
            f"CableDyn {row['cabledyn_Hz']:.3f} Hz, {row['error_percent']:+.1f} %"
        )
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "executable_sha256"}, indent=2
        )
    )
    if args.summary is not None:
        args.summary.mkdir(parents=True, exist_ok=True)
        with (args.summary / "natural_frequencies.csv").open(
            "w", newline="", encoding="utf-8"
        ) as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        (args.summary / "natural_frequencies_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n"
        )


if __name__ == "__main__":
    main()
