#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Coupled linearization A/B: CableDyn (CompMooring=5) vs stock MoorDyn (=3).

Twin still-water IEA-15MW UMaine runs, identical except for the mooring module,
each linearized once at the same time with module-level output (LinOutMod).

WHAT IS GATED -- CableDyn's quasi-static dYdu stiffness block (force rows vs
translation-displacement columns, matched by ROW/COLUMN LABELS) must satisfy
the invariants of a true mooring tangent stiffness at equilibrium:
  * symmetry (a conservative static tangent);
  * negative displacement diagonals (restoring);
  * spread-symmetry: a collective surge produces (near-)zero net sway force
    on this y-symmetric three-line spread.
Measured on this model the block is symmetric to machine noise and its
collective surge stiffness is stable across operating points (~75 kN/m).

WHAT IS REPORTED, NOT GATED -- the same block from stock MoorDyn's
zero-frequency reduction K_md = D - C A^-1 B. Measured on this catenary +
seabed model, MoorDyn's own FD linearization is contact-noisy: its block is
asymmetric (~14 %) and violates the spread's y-symmetry by orders of
magnitude, so it is NOT a usable reference here (the numbers are printed as
evidence). The rigorous cross-code check is the finite-displacement settled
A/B (displaced-start static solves in both codes, `openfast_settled_ab.py`).

Fail-closed: missing .lin files, unparseable sections, empty matched blocks
(e.g. a module whose coupled-mesh inputs did not enter the linearization set),
and non-finite entries all abort with a named error.

Usage:
    python validation/scripts/openfast_lin_ab.py --exe <openfast.exe> --case-dir <dir>
        [--ref-cd Main_settle_cd5.fst] [--ref-md Main_settle_md3.fst]
        [--t-lin 60.0] [--gate-pct 5.0] [--skip-run]
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

LIN_BLOCK = re.compile(r"^-+\s*LINEARIZATION\s*-+", re.IGNORECASE)


def die(msg: str) -> None:
    print(f"FAIL: {msg}")
    sys.exit(1)


def swap_keys(src: Path, dst: Path, swaps: dict[str, str], required: bool, what: str) -> set[str]:
    """Copy a value-first OpenFAST input file, replacing the value field of the
    named keys. The key is the LAST word before ' - ' (LinTimes carries a
    comma-separated value list)."""
    lines = src.read_text().splitlines(keepends=True)
    out: list[str] = []
    seen: set[str] = set()
    for ln in lines:
        m = re.match(r"^\s*[^-\s].*?(\w+)\s+-\s", ln)
        key = m.group(1) if m else None
        if key in swaps:
            rest = ln[ln.index(key):]
            out.append(f"{swaps[key]:<22} {rest}")
            seen.add(key)
        else:
            out.append(ln)
    missing = set(swaps) - seen
    if required and missing:
        die(f"{src.name}: {what} keys not found: {sorted(missing)}")
    dst.write_text("".join(out))
    return seen


def prepare_lin_deck(src: Path, dst: Path, t_lin: float) -> None:
    """Copy an .fst deck, enabling one linearization snapshot at t_lin, and derive
    a linearization-compatible HydroDyn twin (no DFT excitation, no convolution
    radiation, no Newman 2nd-order -- HydroDyn refuses to linearize those)."""
    swap_keys(src, dst, {
        "TMax": f"{t_lin:.1f}",
        "Linearize": "True",
        "CalcSteady": "False",
        "NLinTimes": "1",
        "LinTimes": f"{t_lin:.6f}",
        # 2 = ALL module inputs/outputs: under the "standard" set (1) neither
        # mooring module's coupled MESH enters the linearization set (measured:
        # n_u = 0 for stock MoorDyn too), and the impedance needs the mesh columns.
        "LinInputs": "2",
        "LinOutputs": "2",
        "LinOutJac": "False",
        "LinOutMod": "True",
    }, required=True, what="linearization")

    text = dst.read_text()
    m = re.search(r'"([^"]+)"\s+HydroFile', text)
    if not m:
        die(f"{dst.name}: HydroFile entry not found")
    hydro_src = (dst.parent / m.group(1)).resolve()
    hydro_dst = dst.parent / (hydro_src.stem + "_lin" + hydro_src.suffix)
    seen = swap_keys(hydro_src, hydro_dst, {
        "ExctnMod": "0",
        "RdtnMod": "0",
        "NewmanApp": "0",
    }, required=False, what="HydroDyn linearization")
    if not {"ExctnMod", "RdtnMod"} & seen:
        die(f"{hydro_src.name}: neither ExctnMod nor RdtnMod found -- cannot make "
            f"the HydroDyn deck linearization-compatible")
    dst.write_text(text.replace(m.group(1), hydro_dst.name))


def run_case(exe: Path, fst: Path) -> None:
    proc = subprocess.run(
        [str(exe), str(fst.name)],
        cwd=fst.parent,
        capture_output=True,
        text=True,
        timeout=3600,
    )
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-25:])
        die(f"openfast failed on {fst.name} (rc={proc.returncode}):\n{tail}")


def parse_module_lin(path: Path) -> dict:
    """Parse an OpenFAST module .lin file: sizes, labels, and A/B/C/D matrices."""
    if not path.is_file():
        die(f"missing linearization file: {path}")
    text = path.read_text().splitlines()

    def read_int(pattern: str) -> int:
        for ln in text:
            m = re.search(pattern, ln)
            if m:
                return int(m.group(1))
        return 0

    n_x = read_int(r"Number of continuous states:\s*(\d+)")
    n_u = read_int(r"Number of inputs:\s*(\d+)")
    n_y = read_int(r"Number of outputs:\s*(\d+)")

    def read_labels(header: str, count: int) -> list[str]:
        # Layout: the header line, a column-title row, a dashes underline row, then
        # one entry per row starting with its integer index. Take ONLY rows whose
        # first token is an integer -- a positional skip miscounts (measured: the
        # dashes row became label 1 and every row/column selection shifted by one,
        # silently extracting a shifted block from a CORRECT matrix).
        labels: list[str] = []
        for i, ln in enumerate(text):
            if header in ln:
                j = i + 1
                while j < len(text) and len(labels) < count:
                    row = text[j].strip()
                    j += 1
                    if not row:
                        continue
                    tok = row.split()[0]
                    if not tok.lstrip("+-").isdigit():
                        continue
                    labels.append(row)
                break
        if len(labels) != count:
            die(f"{path.name}: expected {count} labels under '{header}', got {len(labels)}")
        return labels

    u_labels = read_labels("Order of inputs", n_u) if n_u else []
    y_labels = read_labels("Order of outputs", n_y) if n_y else []

    def read_matrix(name: str, rows: int, cols: int) -> np.ndarray:
        if rows == 0 or cols == 0:
            return np.zeros((rows, cols))
        pat = re.compile(rf"^{name}:\s*\d+\s*x\s*\d+")
        for i, ln in enumerate(text):
            if pat.match(ln.strip()):
                # Wide rows WRAP across physical lines in .lin files: accumulate
                # exactly rows*cols floats (never count physical lines), stopping
                # on the first non-numeric token (the next section header).
                flat: list[float] = []
                j = i + 1
                while j < len(text) and len(flat) < rows * cols:
                    toks = text[j].split()
                    j += 1
                    if not toks:
                        continue
                    try:
                        flat.extend(float(v) for v in toks)
                    except ValueError:
                        break
                if len(flat) != rows * cols:
                    die(f"{path.name}: {name} yielded {len(flat)} values, expected {rows * cols}")
                mat = np.array(flat).reshape(rows, cols)
                if not np.all(np.isfinite(mat)):
                    die(f"{path.name}: {name} contains non-finite entries")
                return mat
        die(f"{path.name}: matrix {name} ({rows}x{cols}) not found")
        raise AssertionError

    return {
        "n_x": n_x,
        "u_labels": u_labels,
        "y_labels": y_labels,
        "A": read_matrix("A", n_x, n_x),
        "B": read_matrix("B", n_x, n_u),
        "C": read_matrix("C", n_y, n_x),
        "D": read_matrix("D", n_y, n_u),
    }


def match_indices(labels: list[str], must: list[str]) -> list[int]:
    idx = [i for i, s in enumerate(labels) if all(m.lower() in s.lower() for m in must)]
    return idx


def static_impedance(mod: dict) -> np.ndarray:
    """K = D - C A^+ B (reduces to D when the module carries no states).

    A^+ is the pseudo-inverse: a lumped-mass mooring on a FRICTIONLESS seabed has
    neutral horizontal modes on its grounded nodes (zero restoring force), so A is
    singular and a naive solve fills the block with junk (measured: asymmetric,
    sign-flipped diagonals). The minimum-norm solution projects the neutral modes
    out -- they stay put under a static boundary displacement, which is exactly the
    quasi-static physics."""
    if mod["n_x"] == 0:
        return mod["D"]
    A = mod["A"]
    sv = np.linalg.svd(A, compute_uv=False)
    cond = sv[0] / max(sv[-1], 1e-300)
    if cond > 1e10:
        print(f"note: state matrix cond ~ {cond:.2e} (neutral modes); using the "
              f"pseudo-inverse zero-frequency reduction")
        X = np.linalg.lstsq(A, mod["B"], rcond=1e-8)[0]
    else:
        X = np.linalg.solve(A, mod["B"])
    return mod["D"] - mod["C"] @ X


def block(mod: dict, K: np.ndarray, row_must: list[str], col_must: list[str], what: str) -> np.ndarray:
    rows = match_indices(mod["y_labels"], row_must)
    cols = match_indices(mod["u_labels"], col_must)
    if not rows or not cols:
        die(f"{what}: no labels matched rows={row_must} ({len(rows)}) cols={col_must} ({len(cols)}) "
            f"-- the module's coupled mesh did not enter the linearization set")
    return K[np.ix_(rows, cols)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--ref-cd", default="Main_settle_cd5.fst")
    ap.add_argument("--ref-md", default="Main_settle_md3.fst")
    ap.add_argument("--t-lin", type=float, default=60.0)
    # accepted for compatibility with the other A/B scripts; the pass criteria here are the
    # invariants above, not a percentage
    ap.add_argument("--gate-pct", type=float, default=5.0)
    ap.add_argument("--skip-run", action="store_true")
    args = ap.parse_args()

    case = args.case_dir
    cd_fst = case / "Main_lin_cd5.fst"
    md_fst = case / "Main_lin_md3.fst"
    prepare_lin_deck(case / args.ref_cd, cd_fst, args.t_lin)
    prepare_lin_deck(case / args.ref_md, md_fst, args.t_lin)

    if not args.skip_run:
        for fst in (cd_fst, md_fst):
            print(f"running {fst.name} ...")
            run_case(args.exe, fst)

    cd = parse_module_lin(case / "Main_lin_cd5.1.CD.lin")
    md = parse_module_lin(case / "Main_lin_md3.1.MD.lin")
    print(f"CD: n_x={cd['n_x']} n_u={len(cd['u_labels'])} n_y={len(cd['y_labels'])}")
    print(f"MD: n_x={md['n_x']} n_u={len(md['u_labels'])} n_y={len(md['y_labels'])}")
    if cd["n_x"] != 0:
        die("CableDyn module reported continuous states under linearization "
            "(the dYdu-only reduction must register Nx = 0)")

    K_cd = static_impedance(cd)
    K_md = static_impedance(md)

    # The gated block: CD's quasi-static stiffness invariants.
    kb_cd = block(cd, K_cd, ["force"], ["translation displacement"], "CD stiffness")
    n = kb_cd.shape[0]
    sym = np.linalg.norm(kb_cd - kb_cd.T) / max(np.linalg.norm(kb_cd), 1e-30)
    print(f"CD stiffness block {kb_cd.shape}: symmetry |K-K^T|/|K| = {sym:.2e} (gate 1e-3)")
    diag = np.diag(kb_cd)
    print(f"CD displacement diagonals: min {diag.min():.4g}, max {diag.max():.4g} N/m (gate: all < 0)")
    e = np.zeros(n)
    e[0::3] = 1.0  # collective +x on every fairlead
    fx = float(kb_cd[0::3, :] @ e @ np.ones(n // 3))
    fy = float(kb_cd[1::3, :] @ e @ np.ones(n // 3))
    print(f"CD collective surge stiffness: dFx = {fx:.5g} N/m, net sway leak |dFy/dFx| = "
          f"{abs(fy) / max(abs(fx), 1e-30):.2e} (gate 1e-2 on this y-symmetric spread)")
    ok = sym <= 1e-3 and bool(np.all(diag < 0.0)) and abs(fy) <= 1e-2 * abs(fx)

    # Reported, not gated: stock MoorDyn's zero-frequency reduction on the same
    # block. Measured contact-noisy on catenary + seabed decks (asymmetric,
    # y-symmetry-violating) -- printed as evidence, not used as a reference.
    kb_md = block(md, K_md, ["force"], ["translation displacement"], "MD stiffness")
    sym_md = np.linalg.norm(kb_md - kb_md.T) / max(np.linalg.norm(kb_md), 1e-30)
    fx_md = float(kb_md[0::3, :] @ e @ np.ones(n // 3))
    fy_md = float(kb_md[1::3, :] @ e @ np.ones(n // 3))
    print(f"MD (informational): collective surge dFx = {fx_md:.5g} N/m, asymmetry {sym_md:.2f}, "
          f"sway leak |dFy/dFx| = {abs(fy_md) / max(abs(fx_md), 1e-30):.2f} -- a symmetric spread "
          f"must have ~zero leak; MoorDyn's -lin FD is contact-noisy here and is not the reference")

    if not ok:
        die("CD stiffness block violates the quasi-static invariants (see above)")
    print("PASS")


if __name__ == "__main__":
    main()
