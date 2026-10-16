# SPDX-License-Identifier: Apache-2.0
"""The layouts behind the EI = 0 fold-note tests (driver_ei0_fold_note_* in CMakeLists.txt).

Those tests check that a note is printed or not; this one checks that the static layouts are
the ones the notes are about, so that a silent no-note test cannot pass on a layout that never
reached past its anchor. The anchor of each chain is at x = 450 m, its fairlead at x = 0, and
the note threshold is 1 % of the 450 m span plus 0.5 m, so a grounded node past the anchor by
more than 5 m (x > 455 m) is a run carried past it:

* ei0_hairpin_current.dat (a 2 m/s current, frictionless seabed): the run reaches past the
  anchor, and the driver notes it;
* ei0_hairpin_current_friction.dat (the same with frictionMu 0.5): friction holds the run
  short of the threshold, so there is nothing to note;
* ei0_aslaid_behind_anchor.dat (staticsFriction asLaid, layAzimuth 0, away from the
  fairlead): the run does lie past the anchor, laid there by the deck, which the note skips.

Usage: check_fold_layout.py <cabledyn executable> <tests/data directory>
"""

import pathlib
import subprocess
import sys
import tempfile

ANCHOR_X = 450.0
THRESHOLD = 0.01 * 450.0 + 0.5
FLOOR_Z = -100.0
GROUNDED_TOL = 0.05
NODES = 53


def node_positions(exe, deck):
    """The static x and z of every node of line 1 (deck node 1 = End A ... 53 = End B)."""
    text = deck.read_text(encoding="utf-8")
    head = text.split("--------------------- OUTPUTS")[0]
    channels = [f"L1N{k}px" for k in range(1, NODES + 1)] + [f"L1N{k}pz" for k in range(1, NODES + 1)]
    body = head + "--------------------- OUTPUTS ---\n" + "\n".join(channels) + "\n--------------------- END ---\n"
    with tempfile.TemporaryDirectory() as tmp:
        path = pathlib.Path(tmp) / deck.name
        path.write_text(body, encoding="utf-8", newline="\n")
        run = subprocess.run([exe, str(path), "layout"], cwd=tmp, capture_output=True, text=True)
        if run.returncode != 0:
            raise SystemExit(f"{deck.name}: driver failed ({run.returncode}):\n{run.stdout}\n{run.stderr}")
        rows = (pathlib.Path(tmp) / "layout.out").read_text(encoding="utf-8").splitlines()
    values = [float(v) for v in rows[2].split()[1:]]
    return values[:NODES], values[NODES:]


def farthest_grounded_x(exe, deck):
    xs, zs = node_positions(exe, deck)
    grounded = [x for x, z in zip(xs, zs) if z <= FLOOR_Z + GROUNDED_TOL]
    if not grounded:
        raise SystemExit(f"{deck.name}: no grounded node")
    return max(grounded)


def main():
    exe, data = str(pathlib.Path(sys.argv[1]).resolve()), pathlib.Path(sys.argv[2]).resolve()
    failures = []
    x = farthest_grounded_x(exe, data / "ei0_hairpin_current.dat")
    print(f"frictionless current: grounded run reaches x = {x:.3f} m")
    if not x > ANCHOR_X + THRESHOLD:
        failures.append("the frictionless current deck does not carry the run past its anchor")
    x = farthest_grounded_x(exe, data / "ei0_hairpin_current_friction.dat")
    print(f"current with frictionMu: grounded run reaches x = {x:.3f} m")
    if not x <= ANCHOR_X + THRESHOLD:
        failures.append("frictionMu does not hold the run short of its anchor")
    x = farthest_grounded_x(exe, data / "ei0_aslaid_behind_anchor.dat")
    print(f"as laid behind the anchor: grounded run reaches x = {x:.3f} m")
    if not x > ANCHOR_X + THRESHOLD:
        failures.append("the as-laid run does not lie past its anchor, so the skipped note is not exercised")
    for f in failures:
        print("FAIL: " + f)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
