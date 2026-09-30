# SPDX-License-Identifier: Apache-2.0
"""In-process object queries: lines, points, bodies, rods and output channels.

The parity gates run each deck through the standalone driver and compare its
written channels, per-line node positions and segment tensions with the
in-process queries at the same time. They need ``CABLEDYN_TEST_DRIVER`` (or
``CABLEDYN_DRIVER``) to name a built driver and are skipped otherwise. The
other gates need only the shared library.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("cabledyn._lib")  # loads the shared library or skips

import cabledyn
from cabledyn import _lib
from cabledyn import objects as objects_module
from cabledyn.errors import CableDynError
from cabledyn.model import CableDyn

REPO = Path(__file__).resolve().parents[2]
EXAMPLES = REPO / "examples"
_DRIVER = (
    os.environ.get("CABLEDYN_TEST_DRIVER", "").strip()
    or os.environ.get("CABLEDYN_DRIVER", "").strip()
)
needs_driver = pytest.mark.skipif(
    not _DRIVER or not Path(_DRIVER).is_file(),
    reason="set CABLEDYN_TEST_DRIVER to a built native driver to run the parity gates",
)


def _outputs(text: str, channels: list[str]) -> str:
    """Append OUTPUTS channels to a deck's existing OUTPUTS section."""
    block = "".join(f'"{name}"\n' for name in channels)
    marker = re.search(r"^-+\s*need this line.*$", text, flags=re.MULTILINE)
    assert marker is not None
    return text[: marker.start()] + block + text[marker.start() :]


def _held_chain(tmp_path: Path, tmax: float = 1.0) -> Path:
    text = (EXAMPLES / "dynamic_chain_held.dat").read_text(encoding="utf-8")
    text = text.replace("1     2       1       -", "1     2       1       pt")
    text = text.replace("10.0         TMax", f"{tmax}         TMax")
    text = _outputs(
        text,
        [
            "L1N5px",
            "L1N5pz",
            "L1N5vz",
            "L1N5az",
            "Ten1N5",
            "Curv1N5",
            "BendMom1N5",
            "L1N5Dec",
            "L1N5Azi",
            "Point1Fx",
            "Point1Fz",
            "Point2pz",
        ],
    )
    path = tmp_path / "held_chain.dat"
    path.write_text(text, encoding="utf-8")
    return path


def _lazy_wave(tmp_path: Path) -> Path:
    text = (EXAMPLES / "lazy_wave_modes.dat").read_text(encoding="utf-8")
    text = text.replace("1   1      2      -", "1   1      2      pt")
    text = text.replace("10       nModes", "0        nModes")
    text = _outputs(text, ["L1N96px", "L1N96pz", "Ten1N96", "L1N96Dec"])
    path = tmp_path / "lazy_wave.dat"
    path.write_text(text, encoding="utf-8")
    return path


_BODY_DECK = """Rigid6 buoy for the object-query gates (still water, no body hydrodynamics)
--- LINE TYPES ---
TypeName Diam MassDenInAir EA BA EI Cd_n Cd_t Ca_n Ca_t
(-) (m) (kg/m) (N) (N-s/-) (N-m^2) (-) (-) (-) (-)
poly 0.12 15.0 5.0e7 -1.0 0.0 1.2 0.2 1.0 0.0
--- BODIES ---
ID Type X Y Z Roll Pitch Yaw Mass Vol C33 C44 C55 CdA Ca Ixx Iyy Izz
(-) (-) (m) (m) (m) (deg) (deg) (deg) (kg) (m3) (N/m) (Nm) (Nm) (m2) (-) (kgm2) (kgm2) (kgm2)
1 Rigid6 0.0 0.0 -20.0 0.0 0.0 0.0 2.0e4 40.0 0.0 0.0 0.0 0.0 0.0 3.5e4 3.5e4 3.5e4
--- POINTS ---
ID Type X Y Z Mass Vol CdA Ca
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
1 Body1 1.5 0.0 -2.0 0 0 0 0
2 Body1 -0.75 1.299 -2.0 0 0 0 0
3 Body1 -0.75 -1.299 -2.0 0 0 0 0
4 Fixed 40.0 0.0 -100.0 0 0 0 0
5 Fixed -20.0 34.641 -100.0 0 0 0 0
6 Fixed -20.0 -34.641 -100.0 0 0 0 0
--- LINES ---
ID NodeA NodeB Outputs
(-) (-) (-) (-)
1 1 4 pt
2 2 5 -
3 3 6 -
--- SECTIONS ---
LineID LineType Length NumSegs
(-) (-) (m) (-)
1 poly 86.85 20
2 poly 86.85 20
3 poly 86.85 20
--- OPTIONS ---
9.80665 g
1025.0 rhoW
100.0 WtrDpth
1.0e5 kBot
1.0e4 cBot
0.01 dtM
0.1 TMax
--- OUTPUTS ---
"FairTen1"
"AnchTen2"
"Body1Px"
"Body1Pz"
"Body1Rx"
"Body1Vz"
"Body1Fz"
"Body1Mx"
"Point1pz"
"Point4Fz"
--- need this line ---
"""

_ROD_DECK = """Free rod for the object-query gates
--- LINE TYPES ---
Name Diam Mass EA BA EI Cdn Cdt Can Cat
(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)
line 0.10 80.0 1.0e5 0.0 0.0 1.0 0.0 1.0 0.0
--- ROD TYPES ---
Name Diam Mass Cd Ca CdEnd CaEnd
(-) (m) (kg/m) (-) (-) (-) (-)
rodmat 0.20 50.0 1.0 1.0 0.0 0.0
--- RODS ---
ID RodType Type XA YA ZA XB YB ZB NumSegs Outputs
(-) (-) (-) (m) (m) (m) (m) (m) (m) (-) (-)
1 rodmat Free 0.0 0.0 -4.0 0.0 0.0 -2.0 2 -
--- POINTS ---
ID Type X Y Z Mass Vol CdA Ca
(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)
1 Fixed 2.0 0.0 -5.0 0.0 0.0 0.0 0.0
2 Rod1A 0.0 0.0 -4.0 0.0 0.0 0.0 0.0
3 Rod1B 0.0 0.0 -2.0 0.0 0.0 0.0 0.0
4 Coupled 0.0 0.0 -1.0 0.0 0.0 0.0 0.0
--- LINES ---
ID NodeA NodeB Outputs
(-) (-) (-) (-)
1 2 1 -
2 3 4 -
--- SECTIONS ---
LineID LineType Length NumSegs
(-) (-) (m) (-)
1 line 3.5 3
2 line 1.2 2
--- OPTIONS ---
9.80665 g
1025.0 rhoW
10.0 WtrDpth
1.0e5 kBot
0.005 dtM
0.05 TMax
--- OUTPUTS ---
"Rod1Pz"
"Rod1Rx"
"Rod1N2Pz"
"Rod1TenA"
"Rod1Fz"
--- need this line ---
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _assert_close(
    actual: np.ndarray, expected: np.ndarray, what: str, *, peak: bool = False
) -> None:
    """The driver writes 8 significant digits (ES15.7); compare at that precision.

    ``peak`` scales by the largest magnitude instead: derived quantities of a stiff
    finite-EI line (curvature, bend moment) amplify the round-off difference between
    the driver's and the embedded static solution.
    """
    top = max(1.0, float(np.max(np.abs(expected))))
    scale = np.full_like(expected, top) if peak else np.maximum(np.abs(expected), 1.0e-6 * top)
    err = np.max(np.abs(actual - expected) / scale)
    assert err < 1.0e-6, f"{what}: max relative difference {err:.3e}"


def _driver_run(deck: Path, tmp_path: Path) -> cabledyn.DriverResult:
    return cabledyn.CableDynDriver(_DRIVER).run(deck, tmp_path / (deck.stem + "_run"))


def _check_row(model: CableDyn, result: cabledyn.DriverResult, time: float) -> None:
    main = result.read_main()
    row = int(np.argmin(np.abs(main.values[:, 0] - time)))
    assert abs(main.values[row, 0] - time) < 1.0e-9
    channels = list(main.channels[1:])
    _assert_close(model.channels(channels), main.values[row, 1:], f"channels at t={time}")
    for line_id in result.line_ids:
        line = model.line(line_id)
        with_p = result.read_line_positions(line_id)
        with_t = result.read_line_tensions(line_id)
        prow = int(np.argmin(np.abs(with_p.values[:, 0] - time)))
        _assert_close(line.node_positions().ravel(), with_p.values[prow, 1:], "node positions")
        trow = int(np.argmin(np.abs(with_t.values[:, 0] - time)))
        # Relative to the line's peak tension: a stiff element's tension amplifies the
        # round-off difference between the driver's and the embedded static solution.
        expected = with_t.values[trow, 1:]
        err = np.max(np.abs(line.segment_tensions() - expected)) / np.max(np.abs(expected))
        assert err < 1.0e-6, f"segment tensions: max difference {err:.3e} of the peak"


# --- parity with the standalone driver's output files -------------------------


@needs_driver
def test_point_route_matches_driver_outputs_at_every_sampled_step(tmp_path):
    deck = _held_chain(tmp_path)
    result = _driver_run(deck, tmp_path)
    with CableDyn(deck) as model:
        _check_row(model, result, 0.0)
        for step in range(20):
            model.step_held(0.05)
            if (step + 1) % 5 == 0:
                _check_row(model, result, model.time)


@needs_driver
def test_finite_ei_line_matches_driver_outputs_at_the_initial_state(tmp_path):
    deck = _lazy_wave(tmp_path)
    result = _driver_run(deck, tmp_path)
    with CableDyn(deck) as model:
        (line,) = model.lines
        assert line.finite_ei
        _check_row(model, result, 0.0)
        # the static profile carries curvature and bend moment for every node
        static = result.read_static()
        profile = np.asarray(static.values)
        cols = list(static.channels)
        curvature = profile[:, cols.index("Curvature")]
        moment = profile[:, cols.index("BendMoment")]
        _assert_close(line.curvature(), curvature, "curvature", peak=True)
        _assert_close(line.bend_moment(), moment, "bend moment", peak=True)


@needs_driver
def test_body_and_rod_decks_match_driver_outputs_at_the_initial_state(tmp_path):
    for name, text in (("body.dat", _BODY_DECK), ("rod.dat", _ROD_DECK)):
        deck = _write(tmp_path, name, text)
        result = _driver_run(deck, tmp_path)
        with CableDyn(deck) as model:
            _check_row(model, result, 0.0)


# --- library-only gates --------------------------------------------------------


def test_views_equal_the_channel_evaluator_bit_for_bit(tmp_path):
    deck = _held_chain(tmp_path)
    with CableDyn(deck) as model:
        model.step_held(0.05)
        (line,) = model.lines
        assert (line.index, line.id, line.name, line.n_nodes, line.n_segments) == (
            0,
            1,
            "Line1",
            42,
            41,
        )
        assert not line.finite_ei
        pos = line.node_positions()
        vel = line.node_velocities()
        acc = line.node_accelerations()
        for k in (1, 7, 42):
            for c, axis in enumerate("xyz"):
                assert pos[k - 1, c] == model.channel(f"L1N{k}p{axis}")
                assert vel[k - 1, c] == model.channel(f"L1N{k}v{axis}")
                assert acc[k - 1, c] == model.channel(f"L1N{k}a{axis}")
            assert line.node_tensions()[k - 1] == model.channel(f"Ten1N{k}")
            assert line.curvature()[k - 1] == model.channel(f"Curv1N{k}")
            assert line.declination()[k - 1] == model.channel(f"L1N{k}Dec")
            assert line.azimuth()[k - 1] == model.channel(f"L1N{k}Azi")
        assert np.all(line.bend_moment() == 0.0)
        assert line.fairlead_tension() == model.channel("FairTen1") == line.node_tensions()[0]
        assert line.anchor_tension() == model.channel("AnchTen1") == line.node_tensions()[-1]
        seg = line.segment_tensions()
        assert seg.shape == (41,) and np.all(seg > 0.0)
        arc = line.arc_length()
        assert arc[0] == 0.0 and arc[-1] == pytest.approx(410.0, rel=1.0e-2)

        fixed, coupled = model.points
        assert (fixed.id, fixed.kind, coupled.id, coupled.kind) == (1, "fixed", 2, "coupled")
        np.testing.assert_array_equal(fixed.position(), [400.0, 0.0, -50.0])
        np.testing.assert_array_equal(coupled.position(), pos[0])
        np.testing.assert_array_equal(fixed.position(), pos[-1])
        np.testing.assert_array_equal(coupled.velocity(), [0.0, 0.0, 0.0])
        force = fixed.force()
        assert force[0] == model.channel("Point1Fx")
        assert np.linalg.norm(force) == pytest.approx(line.anchor_tension(), rel=1.0e-12)
        assert model.bodies == () and model.rods == ()


def test_out_buffers_are_filled_in_place_and_validated(tmp_path):
    deck = _held_chain(tmp_path)
    with CableDyn(deck) as model:
        line = model.line(1)
        buf = np.empty((42, 3))
        assert line.node_positions(out=buf) is buf
        np.testing.assert_array_equal(buf, line.node_positions())
        seg = np.empty(41)
        assert line.segment_tensions(out=seg) is seg
        p3 = np.empty(3)
        assert model.point(1).force(out=p3) is p3
        with pytest.raises(ValueError, match="shape"):
            line.node_positions(out=np.empty((3, 42)))
        with pytest.raises(ValueError, match="C-contiguous"):
            line.node_positions(out=np.empty((42, 6))[:, ::2])
        with pytest.raises(ValueError, match="unknown line quantity"):
            line.values("strain")  # type: ignore[arg-type]


def test_finite_ei_views_report_bending(tmp_path):
    deck = _lazy_wave(tmp_path)
    with CableDyn(deck) as model:
        (line,) = model.lines
        # 256 deck segments, refined by the tensile-safety pass
        assert line.finite_ei and line.n_nodes >= 257
        curv = line.curvature()
        moment = line.bend_moment()
        assert curv.max() > 0.01
        # EI is uniform (1.99e4 N.m^2) along this cable
        np.testing.assert_allclose(moment, 1.99e4 * curv, rtol=1.0e-9, atol=1.0e-9)
        assert line.bend_moment()[95] == model.channel("BendMom1N96")
        np.testing.assert_allclose(line.node_positions()[0], [5.0, 0.0, -14.0], atol=1e-9)
        np.testing.assert_allclose(line.node_positions()[-1], [125.0, 0.0, -80.0], atol=1e-9)
        model.step_held(0.05)
        assert model.time == pytest.approx(0.05)
        assert np.all(np.isfinite(line.segment_tensions()))


def test_body_views(tmp_path):
    deck = _write(tmp_path, "body.dat", _BODY_DECK)
    with CableDyn(deck) as model:
        (body,) = model.bodies
        assert (body.id, body.name) == (1, "Body1")
        pose = body.pose()
        assert pose[2] == model.channel("Body1Pz")
        assert pose[3] == model.channel("Body1Rx")
        wrench = body.wrench()
        assert wrench[2] == model.channel("Body1Fz")
        assert wrench[3] == model.channel("Body1Mx")
        assert body.velocity()[2] == model.channel("Body1Vz")
        assert body.acceleration().shape == (6,)
        rot = body.rotation_matrix()
        np.testing.assert_allclose(rot @ rot.T, np.eye(3), atol=1e-12)
        assert len(model.lines) == 3 and len(model.points) >= 6
        assert model.body(1) is body
        model.step_held(0.01)
        assert np.all(np.isfinite(body.pose()))


def test_rod_views(tmp_path):
    deck = _write(tmp_path, "rod.dat", _ROD_DECK)
    with CableDyn(deck) as model:
        (rod,) = model.rods
        assert (rod.id, rod.n_nodes, rod.n_segments) == (1, 3, 2)
        nodes = rod.node_positions()
        np.testing.assert_array_equal(rod.end_a(), nodes[0])
        np.testing.assert_array_equal(rod.end_b(), nodes[-1])
        assert nodes[2, 2] == model.channel("Rod1N2Pz")
        assert rod.pose()[2] == model.channel("Rod1Pz")
        assert rod.pose()[5] == 0.0
        assert rod.wrench()[2] == model.channel("Rod1Fz")
        assert rod.velocity().shape == (6,)
        assert np.linalg.norm(rod.axis()) == pytest.approx(1.0)
        assert model.rod(1) is rod
        with pytest.raises(KeyError):
            model.rod(7)


def test_channel_errors_and_stale_views(tmp_path):
    deck = _held_chain(tmp_path)
    model = CableDyn(deck)
    for bad in ("", "x" * 65, "L1N99px", "Line1", "FairTen9", "L1N5q", "Body1Pz", "TDP1s"):
        with pytest.raises((CableDynError, ValueError)):
            model.channel(bad)
    with pytest.raises(TypeError):
        model.channel(3)  # type: ignore[arg-type]
    with pytest.raises(KeyError):
        model.line(5)
    line = model.lines[0]
    model.init_deck(deck)
    assert model.time == 0.0
    with pytest.raises(CableDynError, match="stale"):
        line.node_positions()
    fresh = model.lines[0]
    assert fresh is not line and fresh.n_nodes == 42
    model.close()
    with pytest.raises(CableDynError, match="stale"):
        fresh.node_tensions()
    assert "Line1" in repr(fresh)


def test_unbound_handle_and_raw_lines(monkeypatch, tmp_path):
    with CableDyn() as model:
        with pytest.raises(CableDynError, match="NOT_INITIALIZED"):
            model.lines  # noqa: B018
    # an old library without the object queries fails with a clear message
    monkeypatch.setattr(_lib, "abi_minor", lambda: 0)
    with CableDyn(_held_chain(tmp_path)) as model:
        with pytest.raises(RuntimeError, match="predates the object queries"):
            model.points  # noqa: B018


def test_abi_minor_is_reported():
    assert _lib.abi_minor() == _lib.OBJECT_QUERY_ABI_MINOR == cabledyn.abi_minor() == 1
    assert objects_module.OBJ_ROD == 4


def test_finite_ei_deck_needs_the_coupling_step(tmp_path):
    text = (EXAMPLES / "lazy_wave_modes.dat").read_text(encoding="utf-8")
    text = text.replace("10       nModes", "0        nModes")
    kept = [row for row in text.splitlines() if " dtM " not in row and " TMax " not in row]
    deck = _write(tmp_path, "no_dtm.dat", "\n".join(kept) + "\n")
    with pytest.raises(CableDynError, match="finite-EI lines, BODIES or RODS needs dtM"):
        CableDyn(deck)


def test_channel_tokens_must_be_ascii(tmp_path):
    with CableDyn(_held_chain(tmp_path)) as model:
        with pytest.raises(ValueError, match=r"1..64 ASCII characters"):
            model.channel("FairTen\u00b91")


def test_zero_length_rod_axis_is_a_value_error(tmp_path, monkeypatch):
    deck = _write(tmp_path, "rod.dat", _ROD_DECK)
    with CableDyn(deck) as model:
        (rod,) = model.rods
        monkeypatch.setattr(type(rod), "node_positions", lambda self, out=None: np.zeros((3, 3)))
        with pytest.raises(ValueError, match="zero length"):
            rod.axis()
