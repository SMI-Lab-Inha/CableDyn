#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""A lazy wave restrained in torsion starts in an in-plane current.

The lazy wave of ``tests/data/lazy_wave_torsion_attachments.dat`` (25 buoyancy modules, Rigid in
torsion at both ends) in a 0.5 m/s current along its plane: its untwisted static state, solved in
the vertical plane, is unstable out of that plane on its mesh (one negative stiffness mode), which
the same deck without torsion accepts with a note. The twist stage added no instability, yet the
torsion stage refused it: at pretwist 0 "torsion continuation stalled at Phi = 0 (target 0): the
twisted equilibrium is unstable", at 360 deg an unstable equilibrium after a failed Newton solve.

Gate, on the driver's own output: at pretwist 0 and 360 deg the deck starts and runs (it stopped
before), the torsion note reports the twist solved, and the instability is named as the untwisted
state's own.

usage: test_torsion_current_start.py --exe <cabledyn> --data <tests/data> --output <dir>
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def deck_text(base: str, pretwist: str) -> str:
    lines = base.split("\n")
    out = []
    for line in lines:
        f = line.split()
        if f and f[0] == "bare" and len(f) == 14:
            f[9] = "0.12"
            line = " ".join(f)
        elif f[:3] == ["4", "B", "Rigid"]:
            line = " ".join(f[:-1] + [pretwist])
        elif f[1:] == ["dtM"]:
            line = "0.05 dtM"
        elif f[1:] == ["TMax"]:
            line = "0.5 TMax\nuniform 0.5 0 0 current"
        out.append(line)
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    base = (Path(args.data) / "lazy_wave_torsion_attachments.dat").read_text(
        encoding="utf-8"
    )
    failures = []
    for pretwist in ("0", "360"):
        deck = out / f"lw_current_twist{pretwist}.dat"
        deck.write_text(deck_text(base, pretwist), encoding="utf-8", newline="\n")
        if "uniform 0.5 0 0 current" not in deck.read_text(encoding="utf-8"):
            failures.append("the deck edit did not place the current")
            continue
        p = subprocess.run(
            [args.exe, str(deck), str(out / f"lw_current_twist{pretwist}")],
            capture_output=True,
            text=True,
            check=False,
        )
        log = p.stdout + p.stderr
        solved = "torsion solved" in log
        named = "keeps the instability of its untwisted state" in log
        print(
            f"pretwist {pretwist}: exit {p.returncode}, solved {solved}, untwisted instability named {named}"
        )
        if p.returncode != 0:
            failures.append(f"pretwist {pretwist}: the deck stopped: {log[-400:]}")
        if not solved:
            failures.append(f"pretwist {pretwist}: no 'torsion solved' note")
        if not named:
            failures.append(
                f"pretwist {pretwist}: the untwisted state's instability is not named"
            )
        if "stalled" in log:
            failures.append(
                f"pretwist {pretwist}: a stalled torsion continuation is reported"
            )
    for f in failures:
        print("FAIL: " + f)
    if not failures:
        print("PASS: a lazy wave restrained in torsion starts in an in-plane current")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
