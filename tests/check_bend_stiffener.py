# SPDX-License-Identifier: Apache-2.0
"""Bend stiffener against the planar elastica (driver_bend_stiffener in CMakeLists.txt).

A power cable leaves a hang-off clamped horizontal at End A through a 5 m polyurethane bend
stiffener (a 0.5 m cylindrical root and a 4.5 m cone) and hangs to a fixed End B 52 m away and
30 m lower. The driver solves the static equilibrium with the stiffener given as a STIFFENERS
table. The reference is an independent solution of the planar extensible elastica,

    x' = (1 + e) cos(theta),  z' = (1 + e) sin(theta),  theta' = M / EI(s),
    M' = -(1 + e) (cos(theta) n_z - sin(theta) n_x),   n' = (0, w(s)),
    e = (n_x cos(theta) + n_z sin(theta)) / EA,

with EI(s) the cable's EI plus E pi/64 (OD(s)^4 - ID^4) of the sleeve and w(s) the submerged
weight of cable and sleeve per unstretched length, solved by collocation (scipy solve_bvp)
with the clamped direction and position at End A, and the position and zero moment of the
pinned End B. Over the stiffener and the 3 m beyond it the cable's own bend moment (the
static profile's BendMoment, EI_cable kappa) must follow the reference within 2 % of its root
value, and the node positions within 1 cm; at the root and 2 m along, the cable's and the
stiffener's moments (BendMom1N<J> + StfMom1N<J>) must add up to the reference section moment
within 2 % of the root moment; and FairDefl1, the angle between the end force and the clamp
axis, must follow the reference within 0.05 deg.

The curvature over the same 8 m must also follow OrcaFlex 11.6d, with the stiffener given as a
line attachment (tests/data/orcaflex_bend_stiffener_curvature.csv), within 1.5 % of the root
curvature; this part needs no scipy.

The initialisation notes are checked too: none on this deck; the resolution note when the 5 m
stiffener spans 3 elements (2 m elements); the pinned-end note when End A has no END
CONNECTIONS row.

Usage: check_bend_stiffener.py <cabledyn executable>
Exit 77 (skipped) when scipy is not installed and the notes pass.
"""

import math
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

try:
    from scipy.integrate import solve_bvp
except ImportError:  # pragma: no cover - environment dependent
    solve_bvp = None

MOMENT_TOLERANCE = 0.02
POSITION_TOLERANCE = 0.01
ANGLE_TOLERANCE = 0.05
ORCAFLEX_TOLERANCE = 0.015
ORCAFLEX_CSV = pathlib.Path(__file__).resolve().parent / "data" / "orcaflex_bend_stiffener_curvature.csv"
WINDOW = 8.0

G = 9.80665
RHO_W = 1025.0
D_CABLE = 0.20
M_CABLE = 60.0
EA = 5.0e8
EI_CABLE = 2.0e4
LENGTH = 60.0
END_B = (51.90, -30.0)
# root side arc, tip side arc, root OD, tip OD (outward from End A)
SEGMENTS = [(0.0, 0.5, 0.60, 0.60), (0.5, 5.0, 0.60, 0.22)]
ID_STF = 0.20
E_STF = 5.0e7
RHO_STF = 1150.0

DECK = f"""--------------------- CableDyn Input File ------------------------------------
Bend stiffener check: hang-off clamped horizontal, static
--------------------- LINE TYPES ---------------------------------------------
TypeName Diam MassDenInAir EA BA/-zeta EI Cd_n Cd_t Ca_n Ca_t
(-) (m) (kg/m) (N) (N-s/-) (N-m^2) (-) (-) (-) (-)
cable {D_CABLE} {M_CABLE} {EA} -1.0 {EI_CABLE} 1.2 0.0 1.0 0.0
--------------------- POINTS -------------------------------------------------
ID Type X Y Z Mass Vol CdA Ca
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Coupled 0.0 0.0 -20.0 0 0 0 0
2 Fixed {END_B[0]} 0.0 {-20.0 + END_B[1]} 0 0 0 0
--------------------- LINES --------------------------------------------------
ID NodeA NodeB Outputs
(-) (-) (-) (-)
1 1 2 -
--------------------- SECTIONS -----------------------------------------------
LineID LineType Length NumSegs
(-) (-) (m) (-)
1 cable 8.0 160
1 cable 52.0 52
--------------------- STIFFENERS ---------------------------------------------
LineID End Length RootOD TipOD InnerD E Density
(-) (-) (m) (m) (m) (m) (Pa) (kg/m^3)
""" + "".join(
    f"1 A {s1 - s0} {od0} {od1} {ID_STF} {E_STF} {RHO_STF}\n" for s0, s1, od0, od1 in SEGMENTS
) + """--------------------- END CONNECTIONS ----------------------------------------
LineID End Stiffness EzX EzY EzZ
(-) (-) (N-m/rad) (-) (-) (-)
1 A Rigid 1.0 0.0 0.0
--------------------- OPTIONS ------------------------------------------------
9.80665 g
1025.0 rhoW
0.05 dtM
0.05 TMax
--------------------- OUTPUTS ------------------------------------------------
"FairTen1"
"FairDefl1"
"BendMom1N1"
"StfMom1N1"
"BendMom1N41"
"StfMom1N41"
END
"""


def outer_diameter(s: np.ndarray) -> np.ndarray:
    od = np.full_like(s, np.nan)
    for s0, s1, od0, od1 in SEGMENTS:
        inside = (s >= s0) & (s <= s1)
        od[inside] = od0 + (s[inside] - s0) / (s1 - s0) * (od1 - od0)
    return od


def section(s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Bending stiffness and submerged weight per unstretched length at arc s from End A."""
    od = outer_diameter(s)
    sleeve = ~np.isnan(od)
    ei = np.full_like(s, EI_CABLE)
    w = np.full_like(s, (M_CABLE - RHO_W * math.pi / 4.0 * D_CABLE**2) * G)
    area = math.pi / 4.0 * (od[sleeve] ** 2 - ID_STF**2)
    ei[sleeve] += E_STF * math.pi / 64.0 * (od[sleeve] ** 4 - ID_STF**4)
    w[sleeve] += (RHO_STF - RHO_W) * area * G
    return ei, w


FORCE_SCALE = 1.0e5
MOMENT_SCALE = 1.0e4


def reference(s_nodes: np.ndarray, x_seed: np.ndarray, z_seed: np.ndarray):
    """Planar extensible elastica along the unstretched arc, y = (x, z, theta, M, n_x, n_z), with
    the moment and the force scaled by MOMENT_SCALE and FORCE_SCALE. The seed is the shape only:
    the driver's node positions; moment zero and a tension of FORCE_SCALE along the line."""

    def rhs(s, y):
        x, z, th, m, nx, nz = y
        ei, w = section(s)
        c, sn = np.cos(th), np.sin(th)
        stretch = 1.0 + FORCE_SCALE * (nx * c + nz * sn) / EA
        return np.vstack([stretch * c, stretch * sn, MOMENT_SCALE * m / ei,
                          -stretch * FORCE_SCALE / MOMENT_SCALE * (c * nz - sn * nx),
                          np.zeros_like(s), w / FORCE_SCALE])

    def bc(ya, yb):
        return np.array([ya[0], ya[1], ya[2], yb[0] - END_B[0], yb[1] - END_B[1], yb[3]])

    # a dense start mesh: the collocation residual is met at 1e-4 without refinement, which
    # would split the mesh at the steps of the sleeve's weight and stiffness at its tip
    mesh = np.unique(np.round(np.concatenate([s_nodes, np.linspace(0.0, WINDOW, 1601),
                                              np.linspace(0.0, LENGTH, 1201)]), 9))
    xs = np.interp(mesh, s_nodes, x_seed)
    zs = np.interp(mesh, s_nodes, z_seed)
    th = np.unwrap(np.arctan2(np.gradient(zs, mesh), np.gradient(xs, mesh)))
    y0 = np.vstack([xs, zs, th, np.zeros_like(mesh), np.cos(th[-1]) * np.ones_like(mesh), np.sin(th)])
    sol = solve_bvp(rhs, bc, mesh, y0, tol=1.0e-4, max_nodes=50000)
    if not sol.success:
        raise RuntimeError(f"reference elastica did not converge: {sol.message}")
    return sol


RESOLUTION_NOTE = "stiffener at End A spans"
PINNED_NOTE = "the stiffener at End A sleeves a pinned end"


def solve(exe: str, tmp: str, name: str, deck_text: str):
    """Run the driver on a deck; return (exit code, stdout + stderr, static profile rows, the
    output channels at t = 0 by name)."""
    deck = pathlib.Path(tmp) / f"{name}.dat"
    deck.write_text(deck_text)
    root = str(pathlib.Path(tmp) / name)
    run = subprocess.run([exe, str(deck), root], capture_output=True, text=True)
    static = pathlib.Path(root + ".static.out")
    rows = static.read_text().splitlines()[3:] if static.exists() else []
    channels = {}
    series = pathlib.Path(root + ".out")
    if series.exists():
        lines = series.read_text().splitlines()
        channels = dict(zip(lines[1].split(), (float(v) for v in lines[2].split())))
    return run.returncode, run.stdout + run.stderr, rows, channels


def check_notes(exe: str, tmp: str, log: str) -> bool:
    ok = RESOLUTION_NOTE not in log and PINNED_NOTE not in log
    print(f"{'ok  ' if ok else 'FAIL'} no stiffener note on the resolved, clamped deck")
    coarse = DECK.replace("1 cable 8.0 160", "1 cable 8.0 4")
    rc, out, _, _ = solve(exe, tmp, "coarse", coarse)
    good = rc == 0 and RESOLUTION_NOTE in out and "spans 3 elements" in out
    print(f"{'ok  ' if good else 'FAIL'} resolution note with 2 m elements")
    ok = ok and good
    start = DECK.index("--------------------- END CONNECTIONS")
    pinned = DECK[:start] + DECK[DECK.index("--------------------- OPTIONS"):]
    # (the deflection angle needs a held end)
    rc, out, _, _ = solve(exe, tmp, "pinned", pinned.replace('"FairDefl1"\n', ""))
    good = rc == 0 and PINNED_NOTE in out
    print(f"{'ok  ' if good else 'FAIL'} pinned-end note without END CONNECTIONS")
    return ok and good


def check_orcaflex(arc: np.ndarray, curvature: np.ndarray) -> bool:
    """The curvature over the first WINDOW m against OrcaFlex 11.6d with its stiffener given as a
    line attachment (tests/data/orcaflex_bend_stiffener_curvature.csv, from
    validation/scripts/orcaflex_bend_stiffener.py), within ORCAFLEX_TOLERANCE of the root value."""
    table = np.loadtxt(ORCAFLEX_CSV, delimiter=",", comments="#", skiprows=2)
    window = arc <= WINDOW
    reference = np.interp(arc[window], table[:, 0], table[:, 1])
    error = np.max(np.abs(curvature[window] - reference)) / table[0, 1]
    ok = error <= ORCAFLEX_TOLERANCE
    print(f"{'ok  ' if ok else 'FAIL'} curvature against OrcaFlex 11.6d: largest difference "
          f"{100 * error:.2f} % of the root curvature {table[0, 1]:.4f} 1/m")
    return ok


def main() -> int:
    exe = sys.argv[1]
    with tempfile.TemporaryDirectory() as tmp:
        rc, log, rows, out = solve(exe, tmp, "stiffener", DECK)
        if rc != 0:
            print(log[-3000:])
            print("FAIL: the driver did not solve the stiffener deck")
            return 1
        notes_ok = check_notes(exe, tmp, log)
    data = np.array([[float(v) for v in r.split()[2:9]] for r in rows if r.strip()])
    orcaflex_ok = check_orcaflex(data[:, 0], data[:, 5])
    if solve_bvp is None:
        print("SKIP: scipy is not installed (the elastica comparison needs solve_bvp)")
        return 77 if notes_ok and orcaflex_ok else 1
    arc, x, z, moment = data[:, 0], data[:, 1], data[:, 3] + 20.0, data[:, 6]
    s_nodes = np.cumsum(np.r_[0.0, np.full(160, 8.0 / 160), np.full(52, 1.0)])
    sol = reference(s_nodes, x, z)
    ref = sol.sol(s_nodes)
    ref[3] *= MOMENT_SCALE
    section_moment = np.abs(ref[3])
    cable_moment = EI_CABLE * section_moment / section(s_nodes)[0]
    window = s_nodes <= WINDOW
    cable_error = np.max(np.abs(np.abs(moment[window]) - cable_moment[window])) / cable_moment[0]
    position_error = np.max(np.hypot(x - ref[0], z - ref[1]))
    print(f"root cable moment: driver {moment[0]:.1f} N m, elastica {cable_moment[0]:.1f} N m")
    print(f"largest cable-moment difference over the first {WINDOW} m: {100 * cable_error:.2f} % of the root value")
    print(f"largest node position difference: {1000 * position_error:.2f} mm")
    total_error = 0.0
    for node in (1, 41):
        total = out[f"BendMom1N{node}"] + out[f"StfMom1N{node}"]
        total_error = max(total_error, abs(total - section_moment[node - 1]) / section_moment[0])
        print(f"node {node}: cable + stiffener {total:.1f} N m, elastica section moment {section_moment[node - 1]:.1f} N m")
    end_force = sol.sol(np.array([0.0]))[4:6, 0]
    deflection = math.degrees(math.atan2(abs(end_force[1]), end_force[0]))
    angle_error = abs(out["FairDefl1"] - deflection)
    print(f"FairDefl1: driver {out['FairDefl1']:.4f} deg, elastica {deflection:.4f} deg")
    ok = (cable_error <= MOMENT_TOLERANCE and position_error <= POSITION_TOLERANCE and notes_ok
          and total_error <= MOMENT_TOLERANCE and angle_error <= ANGLE_TOLERANCE and orcaflex_ok)
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
