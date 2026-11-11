# SPDX-License-Identifier: Apache-2.0
"""VIV response diagnostics (cabledyn.viv) on synthetic line motions with known answers."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import pytest

from cabledyn import LineNodeHistory, VIVResponse, current_steps, viv_response
from cabledyn.results import TimeHistory

L = 100.0
D = 0.1
N_NODES = 41
F = 0.5  # cross-flow frequency [Hz]
DT = 0.01
T = 40.0  # 20 cross-flow periods
FLOW = (0.0, 0.5, 0.0)


def _history(displacement, *, time=None, depth=-50.0, lateral=0.0) -> LineNodeHistory:
    """A straight line along +x at y = lateral, z = depth with displacement(t, s) -> (dy, dz)."""
    t = np.arange(0.0, T, DT) if time is None else time
    s = np.linspace(0.0, L, N_NODES)
    dy, dz = displacement(t[:, None], s[None, :])
    x = np.broadcast_to(s, dy.shape)
    xyz = np.stack((x, lateral + dy, depth + dz), axis=2).reshape(t.size, -1)
    channels = ["Time(s)"] + [f"Node{n}{a}(m)" for n in range(1, N_NODES + 1) for a in "XYZ"]
    return LineNodeHistory(
        Path("synthetic.Line1.p.out"), "", tuple(channels), None, np.column_stack((t, xyz))
    )


def _standing(t, s, a_cf=0.08, a_il=0.02, phase_cf=0.3, phase_il=1.0):
    """Cross-flow (z) mode 1 at F, in-line (y) mode 2 at 2 F."""
    dz = a_cf * np.sin(np.pi * s / L) * np.cos(2.0 * np.pi * F * t + phase_cf)
    dy = a_il * np.sin(2.0 * np.pi * s / L) * np.cos(4.0 * np.pi * F * t + phase_il)
    return dy, dz


def test_local_frame_follows_the_tangent_and_the_flow():
    response = viv_response(_history(_standing), flow=FLOW, diameter=D)
    assert isinstance(response, VIVResponse)
    assert np.allclose(response.tangent, [1.0, 0.0, 0.0])
    assert np.allclose(response.in_line_axis, [0.0, 1.0, 0.0])
    # (e_x, e_y, t) right-handed: cross-flow = t x in-line = +z
    assert np.allclose(response.cross_flow_axis, [0.0, 0.0, 1.0])
    assert np.allclose(response.normal_speed, 0.5)
    assert response.arc_length[-1] == pytest.approx(L)


def test_amplitudes_frequencies_and_the_nondimensional_frequency():
    response = viv_response(_history(_standing), flow=FLOW, diameter=D)
    mid, quarter = N_NODES // 2, N_NODES // 4
    assert response.rms_over_diameter("cross_flow")[mid] == pytest.approx(0.08 / math.sqrt(2) / D)
    assert response.amplitude_over_diameter("cross_flow")[mid] == pytest.approx(0.8, rel=1e-3)
    assert response.amplitude_over_diameter("in_line")[quarter] == pytest.approx(0.2, rel=1e-3)
    assert response.frequency("cross_flow") == pytest.approx(F, abs=1e-3)
    assert response.frequency("in_line") == pytest.approx(2.0 * F, abs=1e-3)
    assert response.frequencies("cross_flow")[mid] == pytest.approx(F, abs=1e-3)
    # the pinned end nodes do not move
    assert math.isnan(response.frequencies("cross_flow")[0])
    assert np.allclose(response.nondimensional_frequency(), F * D / 0.5, rtol=1e-3)


def test_modal_amplitudes_find_the_modes():
    response = viv_response(_history(_standing), flow=FLOW, diameter=D)
    cross = response.modal_amplitudes("cross_flow", n_modes=6)
    assert cross.dominant_mode == 1
    assert cross.rms_over_diameter[0] == pytest.approx(0.08 / math.sqrt(2) / D, rel=1e-9)
    assert np.all(cross.rms_over_diameter[1:] < 1e-9)
    assert cross.residual < 1e-9
    assert cross.coefficients.shape == (response.time.size, 6)
    inline = response.modal_amplitudes("in_line")
    assert inline.dominant_mode == 2
    s = response.arc_length
    custom = response.modal_amplitudes("cross_flow", shapes=[np.sin(2 * np.pi * s / L)])
    assert custom.mode_numbers == (1,) and custom.residual == pytest.approx(1.0, rel=1e-6)


def _main_output(displacement, nodes) -> TimeHistory:
    """The L1N<J>p[xyz] channels of some nodes of the _history line, as a main output."""
    whole = _history(displacement)
    channels = ["Time(s)"]
    columns = [whole.time]
    for node in nodes:
        for axis in "XYZ":
            channels.append(f"L1N{node}p{axis.lower()}")
            columns.append(whole.column(f"Node{node}{axis}(m)"))
    return TimeHistory(Path("run.out"), "", tuple(channels), None, np.column_stack(columns))


def test_default_shapes_need_the_whole_line():
    interior = viv_response(_main_output(_standing, range(5, 38)), flow=FLOW, diameter=D, line_id=1)
    assert not interior.whole_line
    # the sines would vanish at nodes 5 and 37 instead of the line's ends: refused
    with pytest.raises(ValueError, match="the default shapes need the whole line"):
        interior.modal_amplitudes("cross_flow")
    # shapes at the nodes held are fitted as given
    s = np.linspace(0.0, L, N_NODES)[4:37]
    fit = interior.modal_amplitudes("cross_flow", shapes=[np.sin(np.pi * s / L)])
    assert fit.rms_over_diameter[0] == pytest.approx(0.08 / math.sqrt(2) / D, rel=1e-9)
    assert fit.residual < 1e-9
    assert viv_response(_history(_standing), flow=FLOW, diameter=D).whole_line


def test_dependent_shapes_are_refused():
    # ten uniform nodes: sin(9 pi s/S) vanishes at every node, sin(10 pi s/S) = -sin(8 pi s/S)
    n = 10
    t = np.arange(0.0, T, DT)
    s = np.linspace(0.0, L, n)
    dz = 0.05 * np.sin(8.0 * np.pi * s / L)[None, :] * np.cos(2.0 * np.pi * F * t)[:, None]
    xyz = np.stack((np.broadcast_to(s, dz.shape), 0.0 * dz, -50.0 + dz), axis=2)
    channels = ["Time(s)"] + [f"Node{k}{a}(m)" for k in range(1, n + 1) for a in "XYZ"]
    history = LineNodeHistory(
        Path("coarse.Line1.p.out"),
        "",
        tuple(channels),
        None,
        np.column_stack((t, xyz.reshape(t.size, -1))),
    )
    response = viv_response(history, flow=FLOW, diameter=D)
    with pytest.raises(ValueError, match=r"10 shapes are not independent .*\(rank 8\)"):
        response.modal_amplitudes("cross_flow", n_modes=10)
    eight = response.modal_amplitudes("cross_flow", n_modes=8)
    assert eight.dominant_mode == 8
    assert eight.rms_over_diameter[7] == pytest.approx(0.05 / math.sqrt(2) / D, rel=1e-9)
    assert eight.residual < 1e-9


def test_travelling_wave_index_is_nan_without_motion():
    still = viv_response(_history(lambda t, s: (0.0 * t * s, 0.0 * t * s)), flow=FLOW, diameter=D)
    assert math.isnan(still.frequency("cross_flow"))
    assert math.isnan(still.travelling_wave_index("cross_flow"))
    assert math.isnan(still.travelling_wave_index("in_line"))
    along = viv_response(_history(_standing), flow=(1.0, 0.0, 0.0), diameter=D)
    assert math.isnan(along.travelling_wave_index("cross_flow"))
    # still nodes at coordinates with no exact binary form keep a round-off offset from their
    # mean: they have no frequency, on a still line and at the pinned ends of a moving one
    for depth in (-50.1, -0.3, -1234.567):
        calm = viv_response(
            _history(lambda t, s: (0.0 * t * s, 0.0 * t * s), depth=depth),
            flow=FLOW,
            diameter=D,
        )
        assert math.isnan(calm.frequency("cross_flow"))
        assert np.all(np.isnan(calm.frequencies("cross_flow")))
        assert math.isnan(calm.travelling_wave_index("cross_flow"))
        moving = viv_response(_history(_standing, depth=depth), flow=FLOW, diameter=D)
        ends = moving.frequencies("cross_flow")[[0, -1]]
        assert np.all(np.isnan(ends))
        assert moving.frequency("cross_flow") == pytest.approx(F, abs=1e-3)
        # cross-flow motion only, at a lateral coordinate with no exact binary form: no in-line
        # motion, so no orbit phase
        flat = viv_response(
            _history(
                lambda t, s: (
                    0.0 * t * s,
                    0.08 * np.sin(np.pi * s / L) * np.cos(2 * np.pi * F * t),
                ),
                depth=depth,
                lateral=-0.3,
            ),
            flow=FLOW,
            diameter=D,
        )
        assert math.isnan(flat.phase_difference(N_NODES // 2 + 1))
        assert np.all(flat.complex_shape("in_line", 2.0 * F) == 0.0)


def test_main_output_nodes_keep_their_numbers(tmp_path):
    whole = viv_response(_history(_standing), flow=FLOW, diameter=D)
    held = (5, 10, 21, 37)
    some = viv_response(_main_output(_standing, held), flow=FLOW, diameter=D, line_id=1)
    assert some.node_numbers == held
    assert whole.node_numbers == tuple(range(1, N_NODES + 1))
    # a node is taken by its number: node 21 is the mid-span node of the whole line
    for a, b in zip(some.orbit(21), whole.orbit(21), strict=True):
        assert np.array_equal(a, b)
    assert some.phase_difference(10) == pytest.approx(whole.phase_difference(10), abs=1e-9)
    with pytest.raises(ValueError, match=r"one of the nodes held, \[5, 10, 21, 37\]"):
        some.orbit(1)
    with pytest.raises(ValueError, match="from 1 to 41"):
        whole.orbit(42)
    with (some.export(tmp_path / "some.csv")).open(encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    assert [row[0] for row in rows[1:]] == [str(node) for node in held]


def test_travelling_wave_index_separates_standing_and_travelling_waves():
    standing = viv_response(_history(_standing), flow=FLOW, diameter=D)
    assert standing.travelling_wave_index("cross_flow") < 1e-6

    def travelling(t, s):
        dz = 0.05 * np.cos(2.0 * np.pi * F * t - 2.0 * np.pi * 4.0 * s / L)
        return np.zeros_like(dz), dz

    wave = viv_response(_history(travelling), flow=FLOW, diameter=D)
    assert wave.travelling_wave_index("cross_flow") > 0.95


def test_complex_shape_and_phase_difference():
    response = viv_response(_history(_standing), flow=FLOW, diameter=D)
    mid, quarter = N_NODES // 2, N_NODES // 4
    shape = response.complex_shape("cross_flow")
    assert abs(shape[mid]) == pytest.approx(0.08, rel=1e-3)
    assert np.angle(shape[mid]) == pytest.approx(0.3, abs=1e-3)
    # theta = phi_x - 2 phi_y = 1.0 - 0.6 rad
    assert response.phase_difference(quarter + 1) == pytest.approx(math.degrees(0.4), abs=0.1)
    x, y = response.orbit(quarter + 1)
    assert x.shape == y.shape == response.time.shape
    assert math.isnan(response.phase_difference(1))


def test_flow_along_the_line_leaves_no_frame():
    response = viv_response(_history(_standing), flow=(0.4, 0.0, 0.0), diameter=D)
    assert np.all(np.isnan(response.in_line_axis)) and np.all(response.normal_speed == 0.0)
    assert np.all(np.isnan(response.frequencies("cross_flow")))
    assert math.isnan(response.frequency("cross_flow"))
    assert np.all(np.isnan(response.nondimensional_frequency()))
    with pytest.raises(ValueError, match="fewer nodes with a flow frame"):
        response.modal_amplitudes()


def test_window_and_export(tmp_path):
    response = viv_response(_history(_standing), flow=FLOW, diameter=D, start=10.0, stop=30.0)
    assert response.time[0] == pytest.approx(10.0) and response.time[-1] <= 30.0
    target = response.export(tmp_path / "viv.csv")
    with target.open(encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    assert rows[0][:3] == ["Node", "ArcLength_[m]", "NormalSpeed_[m/s]"]
    assert len(rows) == N_NODES + 1
    assert rows[1][7] == ""  # the pinned End A node has no frequency
    with pytest.raises(FileExistsError):
        response.export(target)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"flow": (0.0, 0.0)}, "finite, non-zero 3-vector"),
        ({"flow": (0.0, 0.0, 0.0)}, "finite, non-zero 3-vector"),
        ({"diameter": 0.0}, "diameter must be finite and positive"),
        ({"start": 39.97}, "at least eight samples"),
        ({"line_id": 1}, "line_id does not apply"),
    ],
)
def test_arguments_are_validated(kwargs, message):
    arguments = {"flow": FLOW, "diameter": D, **kwargs}
    with pytest.raises(ValueError, match=message):
        viv_response(_history(_standing), **arguments)


def test_non_uniform_sampling_is_refused():
    time = np.cumsum(np.full(100, DT))
    time[50:] += 0.003
    with pytest.raises(ValueError, match="uniform sampling"):
        viv_response(_history(_standing, time=time), flow=FLOW, diameter=D)


def test_methods_validate_direction_and_node():
    response = viv_response(_history(_standing), flow=FLOW, diameter=D)
    with pytest.raises(ValueError, match="direction must be one of"):
        response.rms_over_diameter("axial")
    with pytest.raises(ValueError, match="node must be an integer from 1 to 41"):
        response.orbit(0)
    with pytest.raises(ValueError, match="n_modes must be a positive integer"):
        response.modal_amplitudes(n_modes=0)
    with pytest.raises(ValueError, match="shapes must have shape"):
        response.modal_amplitudes(shapes=np.ones((2, 3)))


def test_current_steps_schedule_and_file(tmp_path):
    steps = current_steps([0.1, 0.2, 0.1], hold=20.0, ramp=5.0, direction=(0.0, 2.0, 0.0))
    np.testing.assert_allclose(steps.starts, [0.0, 25.0, 50.0])
    assert steps.duration == pytest.approx(70.0)
    np.testing.assert_allclose(steps.direction, [0.0, 1.0, 0.0])
    assert steps.windows(15.0) == ((15.0, 20.0), (40.0, 45.0), (65.0, 70.0))
    path = steps.write(tmp_path, depths=(-100.0, 0.0))
    assert path.name == "current_profile_dynamic.txt"
    lines = path.read_text(encoding="ascii").splitlines()
    assert lines[4].split() == ["-100", "0"]
    rows = np.array([[float(v) for v in line.split()] for line in lines[6:]])
    # t ux(z1 z2) uy(z1 z2) uz(z1 z2): held, ramped over 5 s, held
    np.testing.assert_allclose(rows[:, 0], [0.0, 20.0, 25.0, 45.0, 50.0, 70.0])
    np.testing.assert_allclose(rows[:, 3], [0.1, 0.1, 0.2, 0.2, 0.1, 0.1])
    assert np.all(rows[:, [1, 2, 5, 6]] == 0.0) and np.all(rows[:, 3] == rows[:, 4])
    with pytest.raises(FileExistsError):
        steps.write(tmp_path)
    with pytest.raises(ValueError, match="settle"):
        steps.windows(20.0)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"speeds": []}, "speeds must be"),
        ({"speeds": [0.1, 0.0]}, "speeds must be"),
        ({"ramp": 0.0}, "hold and ramp"),
        ({"hold": 22.0}, "whole number of ramps"),
        ({"direction": (0.0, 0.0, 0.0)}, "direction must be"),
    ],
)
def test_current_steps_validation(kwargs, message):
    arguments = {"speeds": [0.1, 0.2], "hold": 20.0, "ramp": 5.0, **kwargs}
    with pytest.raises(ValueError, match=message):
        current_steps(**arguments)


def test_current_steps_responses_take_each_hold():
    steps = current_steps([0.5, 0.25], hold=15.0, ramp=5.0)
    responses = steps.responses(_history(_standing), diameter=D, settle=5.0)
    assert len(responses) == 2
    assert responses[1].time[0] == pytest.approx(25.0)
    assert np.allclose(responses[0].normal_speed, 0.5) and np.allclose(
        responses[1].normal_speed, 0.25
    )
