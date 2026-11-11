# SPDX-License-Identifier: Apache-2.0
"""Vortex-induced vibration (VIV) response of a line from its node positions.

:func:`viv_response` resolves the motion of every node of a line about its mean
position into the local flow frame of the node and reports the quantities a VIV
study describes a response by: the cross-flow and in-line amplitude along the
span, the response frequencies and the non-dimensional frequency they give
against the flow, the modal content, how far the response is a travelling wave,
and the in-line/cross-flow orbits.

The local frame of a node follows the mean shape over the analysed window: the
tangent ``t`` from End A towards End B (central differences of the mean node
positions), the in-line direction ``e_x`` along the part of the flow velocity
normal to the tangent, and the cross-flow direction ``e_y = t x e_x``, so that
``(e_x, e_y, t)`` is right-handed. A node whose tangent is parallel to the flow
has no normal flow and no frame; its values are ``nan``.

The positions come from a CableDyn per-line position file
(``<root>.Line<L>.p.out``, a :class:`~cabledyn.LineNodeHistory`), a MoorDyn
per-line file with node positions (a :class:`~cabledyn.MoorDynLineHistory`), or
the ``L<L>N<J>px``/``py``/``pz`` channels of a main output for the nodes the
deck lists. The frequency analysis needs uniform sampling.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from cabledyn._csv import number, write_csv
from cabledyn._optional import pyplot
from cabledyn.formats import MoorDynLineHistory
from cabledyn.ranges import node_position_channels
from cabledyn.results import LineNodeHistory, TimeHistory
from cabledyn.spectra import _sampling
from cabledyn.touchdown import _positions

__all__ = ["CurrentSteps", "ModalAmplitudes", "VIVResponse", "current_steps", "viv_response"]

_DIRECTIONS = ("cross_flow", "in_line")


def _readonly(values: npt.ArrayLike, dtype: npt.DTypeLike = np.float64) -> npt.NDArray[Any]:
    array = np.array(values, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _direction(direction: str) -> str:
    if direction not in _DIRECTIONS:
        raise ValueError(f"direction must be one of {_DIRECTIONS}, not {direction!r}")
    return direction


def _peak_frequency(frequency: npt.NDArray[Any], power: npt.NDArray[Any]) -> float:
    """Frequency of the largest power above zero, refined by a parabola through its neighbours."""
    if power.size < 3 or not np.any(power[1:] > 0.0):
        return math.nan
    k = int(np.argmax(power[1:])) + 1
    if 1 <= k < power.size - 1:
        left, middle, right = power[k - 1], power[k], power[k + 1]
        curvature = left - 2.0 * middle + right
        if curvature < 0.0:
            shift = 0.5 * (left - right) / curvature
            return float(frequency[k] + shift * (frequency[1] - frequency[0]))
    return float(frequency[k])


@dataclass(frozen=True)
class ModalAmplitudes:
    """The modal content of one direction of a VIV response.

    Attributes
    ----------
    direction : str
        ``"cross_flow"`` or ``"in_line"``.
    mode_numbers : tuple[int, ...]
        One-based numbers of the shapes fitted.
    coefficients : numpy.ndarray
        Read-only ``(n_samples, n_modes)`` least-squares weights of the shapes
        in the displacement at every sample, in metres per unit shape.
    rms_over_diameter : numpy.ndarray
        Read-only ``(n_modes,)`` RMS of each weight over the window, divided by
        the diameter.
    residual : float
        RMS of the displacement the shapes leave unexplained, over the RMS of
        the displacement (0 when the shapes reproduce it).
    """

    direction: str
    mode_numbers: tuple[int, ...]
    coefficients: npt.NDArray[Any]
    rms_over_diameter: npt.NDArray[Any]
    residual: float

    @property
    def dominant_mode(self) -> int:
        """Number of the shape with the largest RMS weight."""
        return self.mode_numbers[int(np.argmax(self.rms_over_diameter))]


@dataclass(frozen=True)
class VIVResponse:
    """The VIV response of one line over an analysis window.

    Created by :func:`viv_response`. Nodes are in End-A-to-End-B order of the
    source; node arrays have one entry per node, time arrays one per sample.

    Attributes
    ----------
    time : numpy.ndarray
        Read-only ``(n_samples,)`` sample times of the window, in seconds.
    arc_length : numpy.ndarray
        Read-only ``(n_nodes,)`` chord length of the mean shape from End A, m.
    mean_position : numpy.ndarray
        Read-only ``(n_nodes, 3)`` mean node positions over the window, m.
    tangent, in_line_axis, cross_flow_axis : numpy.ndarray
        Read-only ``(n_nodes, 3)`` unit vectors of each node's local frame
        (``nan`` rows where the flow is along the tangent).
    normal_speed : numpy.ndarray
        Read-only ``(n_nodes,)`` speed of the flow normal to the mean tangent,
        m/s.
    in_line, cross_flow : numpy.ndarray
        Read-only ``(n_samples, n_nodes)`` displacement about the mean position
        along the node's in-line and cross-flow axes, m.
    diameter : float
        Diameter the amplitudes are divided by, m.
    flow : numpy.ndarray
        Read-only ``(3,)`` flow velocity the frame is taken against, m/s.
    whole_line : bool
        Whether the nodes are the whole line, End A to End B (a per-line
        file), so that ``arc_length`` spans it; a main output may hold only
        some of its nodes.
    node_numbers : tuple of int
        The CableDyn number of each node (1 = End A): ``1 .. n_nodes`` for a
        per-line file (MoorDyn's node ``J`` is ``J + 1``), the ``J`` of the
        ``L<L>N<J>p`` channels for a main output. The methods that take a
        node take this number. Empty means ``1 .. n_nodes``.
    """

    time: npt.NDArray[Any]
    arc_length: npt.NDArray[Any]
    mean_position: npt.NDArray[Any]
    tangent: npt.NDArray[Any]
    in_line_axis: npt.NDArray[Any]
    cross_flow_axis: npt.NDArray[Any]
    normal_speed: npt.NDArray[Any]
    in_line: npt.NDArray[Any]
    cross_flow: npt.NDArray[Any]
    diameter: float
    flow: npt.NDArray[Any]
    whole_line: bool = True
    node_numbers: tuple[int, ...] = ()

    @property
    def node_count(self) -> int:
        """Number of nodes."""
        return int(self.arc_length.size)

    @property
    def sampling_interval(self) -> float:
        """Uniform sampling interval of the window, s."""
        return float(self.time[1] - self.time[0])

    def displacement(self, direction: str) -> npt.NDArray[Any]:
        """Return the ``(n_samples, n_nodes)`` displacement of ``"cross_flow"`` or ``"in_line"``."""
        return self.cross_flow if _direction(direction) == "cross_flow" else self.in_line

    def rms_over_diameter(self, direction: str = "cross_flow") -> npt.NDArray[Any]:
        """Return the ``(n_nodes,)`` RMS displacement over the diameter.

        For a harmonic response the amplitude is ``sqrt(2)`` times the RMS.
        """
        values = self.displacement(direction)
        return _readonly(np.sqrt(np.mean(values * values, axis=0)) / self.diameter)

    def amplitude_over_diameter(self, direction: str = "cross_flow") -> npt.NDArray[Any]:
        """Return the ``(n_nodes,)`` half peak-to-peak displacement over the diameter."""
        values = self.displacement(direction)
        return _readonly(0.5 * (values.max(axis=0) - values.min(axis=0)) / self.diameter)

    def _moving_displacement(self, direction: str) -> npt.NDArray[Any]:
        """The displacement of a direction with still nodes set to exactly zero.

        A node whose displacement stays within round-off of its coordinates (a range below
        1e-12 of its distance from the origin plus the diameter) is still: its round-off
        offset from the mean, windowed or transformed, would otherwise read as motion at
        every frequency, with an arbitrary phase.
        """
        values = np.array(self.displacement(direction), dtype=np.float64)
        scale = np.linalg.norm(self.mean_position, axis=1) + self.diameter
        with np.errstate(invalid="ignore"):
            still = (values.max(axis=0) - values.min(axis=0)) <= 1.0e-12 * scale
        values[:, still] = 0.0
        return values

    def _power(self, direction: str) -> tuple[npt.NDArray[Any], npt.NDArray[Any]]:
        """Frequencies and the ``(n_freq, n_nodes)`` Hann-windowed periodogram of a direction.

        Still nodes (:meth:`_moving_displacement`) have no power.
        """
        values = self._moving_displacement(direction)
        window = np.hanning(self.time.size)
        spectrum = np.fft.rfft(values * window[:, None], axis=0)
        frequency = np.fft.rfftfreq(self.time.size, self.sampling_interval)
        return frequency, np.abs(spectrum) ** 2

    def frequencies(self, direction: str = "cross_flow") -> npt.NDArray[Any]:
        """Return the ``(n_nodes,)`` dominant response frequency of each node, Hz.

        The peak of the node's Hann-windowed periodogram above zero frequency,
        refined by a parabola through the bins beside it; ``nan`` for a node
        without a frame or without motion.
        """
        frequency, power = self._power(direction)
        result = np.full(self.node_count, math.nan)
        for node in range(self.node_count):
            if np.all(np.isfinite(power[:, node])):
                result[node] = _peak_frequency(frequency, power[:, node])
        return _readonly(result)

    def frequency(self, direction: str = "cross_flow") -> float:
        """Return the dominant response frequency of the line, Hz.

        The peak of the sum of the nodes' periodograms (``nan`` without motion).
        """
        frequency, power = self._power(direction)
        finite = np.all(np.isfinite(power), axis=0)
        if not np.any(finite):
            return math.nan
        return _peak_frequency(frequency, power[:, finite].sum(axis=1))

    def nondimensional_frequency(self, direction: str = "cross_flow") -> npt.NDArray[Any]:
        """Return ``f D / U_n`` at every node, with ``f`` the line's dominant frequency.

        For the cross-flow response this is the quantity a synchronisation
        model's band ``cF +/- dF`` bounds, and its reciprocal the reduced
        velocity ``U_n/(f D)``.
        """
        with np.errstate(divide="ignore", invalid="ignore"):
            result = self.frequency(direction) * self.diameter / self.normal_speed
        return _readonly(np.where(self.normal_speed > 0.0, result, math.nan))

    def modal_amplitudes(
        self,
        direction: str = "cross_flow",
        *,
        n_modes: int = 10,
        shapes: npt.ArrayLike | None = None,
    ) -> ModalAmplitudes:
        """Fit mode shapes to the displacement of a direction at every sample.

        Parameters
        ----------
        direction : str
            ``"cross_flow"`` or ``"in_line"``.
        n_modes : int
            Number of the default shapes ``sin(n pi s / S)``, ``n = 1 ..
            n_modes``, of arc length ``s`` from End A over the mean chord
            length ``S`` (the modes of a taut string with fixed ends). Ignored
            when ``shapes`` is given.
        shapes : array_like, optional
            ``(n_modes, n_nodes)`` scalar shapes of the direction at the nodes
            to fit instead, for example a component of the run's mode shapes.

        Returns
        -------
        ModalAmplitudes
            The weights of the shapes at every sample (least squares over the
            nodes with a frame), their RMS over the diameter and the share of
            the displacement they leave unexplained.

        Raises
        ------
        ValueError
            If the shapes are not finite, have the wrong number of nodes, or
            fewer nodes with a frame than shapes remain; if the shapes are not
            independent at those nodes (a sine that vanishes at every node, or
            two that alias each other on a coarse output); or if the default
            shapes are asked of a response that does not hold the whole line.
        """
        values = self.displacement(direction)
        if shapes is None:
            if isinstance(n_modes, bool) or int(n_modes) != n_modes or n_modes < 1:
                raise ValueError("n_modes must be a positive integer")
            if not self.whole_line:
                # the sines vanish at the line's ends, which the nodes of a main output may
                # not include: fitted over the span of the nodes it holds, they would be wrong
                raise ValueError(
                    "the default shapes need the whole line, End A to End B: analyse the "
                    "per-line .p.out file, or give shapes at the nodes held"
                )
            span = float(self.arc_length[-1])
            if not span > 0.0:
                raise ValueError("the mean shape has no length")
            numbers = np.arange(1, int(n_modes) + 1)
            basis = np.sin(np.pi * numbers[:, None] * self.arc_length[None, :] / span)
        else:
            basis = np.array(shapes, dtype=np.float64)
            if basis.ndim != 2 or basis.shape[1] != self.node_count or basis.shape[0] < 1:
                raise ValueError(f"shapes must have shape (n_modes, {self.node_count})")
            if not np.all(np.isfinite(basis)):
                raise ValueError("shapes must be finite")
            numbers = np.arange(1, basis.shape[0] + 1)
        framed = np.all(np.isfinite(values), axis=0)
        if int(np.count_nonzero(framed)) < basis.shape[0]:
            raise ValueError("fewer nodes with a flow frame than shapes to fit")
        design = basis[:, framed].T
        target = values[:, framed].T
        coefficients, _, rank, _ = np.linalg.lstsq(design, target, rcond=None)
        if int(rank) < basis.shape[0]:
            # dependent shapes at the nodes held (a mode that vanishes at every node, or two
            # that alias each other): their weights are not determined by the motion
            raise ValueError(
                f"the {basis.shape[0]} shapes are not independent at the nodes with a flow frame "
                f"(rank {int(rank)}): fit fewer modes, or use a finer output"
            )
        left = target - design @ coefficients
        total = float(np.sqrt(np.mean(target * target)))
        residual = float(np.sqrt(np.mean(left * left))) / total if total > 0.0 else 0.0
        rms = np.sqrt(np.mean(coefficients * coefficients, axis=1)) / self.diameter
        return ModalAmplitudes(
            direction=_direction(direction),
            mode_numbers=tuple(int(n) for n in numbers),
            coefficients=_readonly(coefficients.T),
            rms_over_diameter=_readonly(rms),
            residual=residual,
        )

    def complex_shape(
        self, direction: str = "cross_flow", frequency: float | None = None
    ) -> npt.NDArray[Any]:
        """Return the ``(n_nodes,)`` complex Fourier shape of a direction at one frequency.

        ``c_j = (2/W) sum_k w_k u_j(t_k) exp(-2 pi i f t_k)`` with the Hann
        window ``w`` of sum ``W``, so that a node moving as ``a cos(2 pi f t +
        phase)`` has ``c_j = a exp(i phase)``; ``frequency`` defaults to the
        line's dominant frequency of the direction. A still node has ``c_j = 0``.
        """
        f = self.frequency(direction) if frequency is None else float(frequency)
        if not (math.isfinite(f) and f >= 0.0):
            raise ValueError("the frequency must be finite and non-negative")
        window = np.hanning(self.time.size)
        phasor = window * np.exp(-2j * np.pi * f * self.time)
        shape: npt.NDArray[Any] = np.asarray(
            2.0 * (phasor @ self._moving_displacement(direction)) / window.sum(),
            dtype=np.complex128,
        )
        shape.setflags(write=False)
        return shape

    def travelling_wave_index(
        self, direction: str = "cross_flow", frequency: float | None = None
    ) -> float:
        """Return how far the response of a direction is a travelling wave, from 0 to 1.

        The ratio of the smaller to the larger singular value of the real and
        imaginary parts of :meth:`complex_shape` (the travelling index of the
        complex orthogonal decomposition, Feeny 2008): 0 for a standing wave,
        whose nodes move in or out of phase, and 1 for a wave of constant
        amplitude travelling over whole wavelengths. ``nan`` without motion.
        """
        if frequency is None and not math.isfinite(self.frequency(direction)):
            return math.nan  # no motion, or no node with a flow frame
        shape = self.complex_shape(direction, frequency)
        framed = np.isfinite(shape)
        parts = np.vstack((shape[framed].real, shape[framed].imag))
        if parts.shape[1] < 2:
            return math.nan
        singular = np.linalg.svd(parts, compute_uv=False)
        return float(singular[1] / singular[0]) if singular[0] > 0.0 else math.nan

    def orbit(self, node: int) -> tuple[npt.NDArray[Any], npt.NDArray[Any]]:
        """Return the in-line and cross-flow displacement of a node over the window.

        ``node`` is the node's number (:attr:`node_numbers`).
        """
        index = self._node_index(node)
        return self.in_line[:, index], self.cross_flow[:, index]

    def phase_difference(self, node: int, frequency: float | None = None) -> float:
        """Return the in-line phase less twice the cross-flow phase at a node, in degrees.

        ``theta = arg(c_x(2 f)) - 2 arg(c_y(f))`` in ``[0, 360)``, with ``c``
        the node's complex Fourier coefficients (:meth:`complex_shape`) and
        ``f`` the cross-flow frequency (default: the line's dominant one): the
        phase that sets the shape and direction of a figure-of-eight orbit.
        ``nan`` for a node without a frame or without motion at those
        frequencies.
        """
        index = self._node_index(node)
        f = self.frequency("cross_flow") if frequency is None else float(frequency)
        if not math.isfinite(f):
            return math.nan
        cy = self.complex_shape("cross_flow", f)[index]
        cx = self.complex_shape("in_line", 2.0 * f)[index]
        if not (np.isfinite(cx) and np.isfinite(cy)) or abs(cx) == 0.0 or abs(cy) == 0.0:
            return math.nan
        return float(np.degrees(np.angle(cx) - 2.0 * np.angle(cy)) % 360.0)

    def _numbers(self) -> tuple[int, ...]:
        """The node numbers, ``1 .. n_nodes`` when none were recorded."""
        return self.node_numbers or tuple(range(1, self.node_count + 1))

    def _node_index(self, node: int) -> int:
        numbers = self._numbers()
        if isinstance(node, bool) or int(node) != node or int(node) not in numbers:
            if numbers == tuple(range(1, self.node_count + 1)):
                raise ValueError(f"node must be an integer from 1 to {self.node_count}")
            raise ValueError(f"node must be one of the nodes held, {list(numbers)}")
        return numbers.index(int(node))

    def export(self, path: str | os.PathLike[str], *, overwrite: bool = False) -> Path:
        """Atomically write the node values as a unit-labelled CSV, one row per node.

        Columns: node number (:attr:`node_numbers`), arc length, normal flow speed,
        cross-flow and in-line RMS
        and half peak-to-peak displacement over the diameter, the nodes'
        dominant frequencies and the non-dimensional cross-flow frequency.
        Fields of a node without a frame are empty.

        Raises
        ------
        FileExistsError
            If ``path`` exists and ``overwrite`` is false.
        """
        columns = (
            self.arc_length,
            self.normal_speed,
            self.rms_over_diameter("cross_flow"),
            self.rms_over_diameter("in_line"),
            self.amplitude_over_diameter("cross_flow"),
            self.amplitude_over_diameter("in_line"),
            self.frequencies("cross_flow"),
            self.frequencies("in_line"),
            self.nondimensional_frequency("cross_flow"),
        )

        numbers = self._numbers()

        def row(index: int) -> tuple[str, ...]:
            fields = tuple(
                number(column[index]) if math.isfinite(float(column[index])) else ""
                for column in columns
            )
            return (str(numbers[index]), *fields)

        return write_csv(
            path,
            (
                "Node",
                "ArcLength_[m]",
                "NormalSpeed_[m/s]",
                "CrossFlowRMSOverD_[-]",
                "InLineRMSOverD_[-]",
                "CrossFlowAmplitudeOverD_[-]",
                "InLineAmplitudeOverD_[-]",
                "CrossFlowFrequency_[Hz]",
                "InLineFrequency_[Hz]",
                "CrossFlowNondimensionalFrequency_[-]",
            ),
            (row(index) for index in range(self.node_count)),
            overwrite=overwrite,
        )

    def plot_envelope(self, *, ax: Any = None) -> Any:
        """Plot the cross-flow and in-line RMS displacement over the diameter along the line."""
        plt = pyplot()
        if ax is None:
            _, ax = plt.subplots()
        ax.plot(self.arc_length, self.rms_over_diameter("cross_flow"), label="cross-flow")
        ax.plot(self.arc_length, self.rms_over_diameter("in_line"), label="in-line")
        ax.set_xlabel("Arc length from End A (m)")
        ax.set_ylabel("RMS displacement / D (-)")
        ax.legend()
        return ax

    def plot_orbit(self, node: int, *, ax: Any = None) -> Any:
        """Plot a node's cross-flow displacement against its in-line displacement, over D."""
        plt = pyplot()
        if ax is None:
            _, ax = plt.subplots()
        x, y = self.orbit(node)
        ax.plot(x / self.diameter, y / self.diameter)
        ax.set_xlabel("In-line displacement / D (-)")
        ax.set_ylabel("Cross-flow displacement / D (-)")
        ax.set_aspect("equal", adjustable="datalim")
        return ax


def viv_response(
    history: TimeHistory,
    *,
    flow: npt.ArrayLike,
    diameter: float,
    start: float | None = None,
    stop: float | None = None,
    line_id: int | None = None,
    uniform_rtol: float = 1.0e-6,
) -> VIVResponse:
    """Resolve a line's node motion into its VIV response.

    Parameters
    ----------
    history : TimeHistory
        A :class:`~cabledyn.LineNodeHistory`, a :class:`~cabledyn.MoorDynLineHistory`
        with node positions, or a main output with the line's
        ``L<L>N<J>p{x,y,z}`` channels.
    flow : array_like
        ``(3,)`` flow velocity the local frames are taken against, m/s: the
        current at the line (a uniform one, or a representative velocity of a
        sheared one).
    diameter : float
        Diameter the amplitudes are divided by, m.
    start, stop : float, optional
        The window to analyse, s (default: the whole record); a lock-in study
        leaves out the build-up.
    line_id : int, optional
        The line of a main output (required there, refused for a per-line file).
    uniform_rtol : float
        Relative tolerance of the uniform-sampling check.

    Returns
    -------
    VIVResponse
        The displacement in the local flow frames and the response measures.

    Raises
    ------
    ValueError
        If the flow is not a finite non-zero 3-vector, the diameter is not
        positive, the window holds fewer than eight samples, the sampling is
        not uniform, or the mean shape has coincident nodes.
    """
    velocity = np.array(flow, dtype=np.float64)
    if velocity.shape != (3,) or not np.all(np.isfinite(velocity)) or not np.any(velocity):
        raise ValueError("flow must be a finite, non-zero 3-vector")
    if not (math.isfinite(float(diameter)) and float(diameter) > 0.0):
        raise ValueError("diameter must be finite and positive")
    time, positions = _positions(history, line_id)
    whole_line = isinstance(history, (LineNodeHistory, MoorDynLineHistory))
    if whole_line or line_id is None:
        numbers = tuple(range(1, positions.shape[1] + 1))
    else:
        # the L<L>N<J>p nodes _positions took, in the same (ascending) order
        numbers = tuple(node_position_channels(history, line_id))
    keep = np.ones(time.size, dtype=bool)
    if start is not None:
        keep &= time >= float(start)
    if stop is not None:
        keep &= time <= float(stop)
    time = time[keep]
    positions = positions[keep]
    if time.size < 8:
        raise ValueError("the window must hold at least eight samples")
    _sampling(time, rtol=float(uniform_rtol), atol=0.0)
    if positions.shape[1] < 2:
        raise ValueError("the line must have at least two nodes")
    mean = positions.mean(axis=0)
    chords = np.diff(mean, axis=0)
    lengths = np.linalg.norm(chords, axis=1)
    if np.any(lengths <= 0.0):
        raise ValueError("the mean shape has coincident nodes")
    tangent = np.empty_like(mean)
    tangent[0] = chords[0]
    tangent[-1] = chords[-1]
    tangent[1:-1] = mean[2:] - mean[:-2]
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    normal = velocity[None, :] - (tangent @ velocity)[:, None] * tangent
    speed = np.linalg.norm(normal, axis=1)
    framed = speed > 1.0e-9 * float(np.linalg.norm(velocity))
    in_axis = np.full_like(mean, math.nan)
    in_axis[framed] = normal[framed] / speed[framed, None]
    cross_axis = np.cross(tangent, in_axis)
    motion = positions - mean[None, :, :]
    in_line = np.einsum("tnk,nk->tn", motion, in_axis)
    cross_flow = np.einsum("tnk,nk->tn", motion, cross_axis)
    return VIVResponse(
        time=_readonly(time),
        arc_length=_readonly(np.concatenate(([0.0], np.cumsum(lengths)))),
        mean_position=_readonly(mean),
        tangent=_readonly(tangent),
        in_line_axis=_readonly(in_axis),
        cross_flow_axis=_readonly(cross_axis),
        normal_speed=_readonly(np.where(framed, speed, 0.0)),
        in_line=_readonly(in_line),
        cross_flow=_readonly(cross_flow),
        diameter=float(diameter),
        flow=_readonly(velocity),
        whole_line=whole_line,
        node_numbers=numbers,
    )


@dataclass(frozen=True)
class CurrentSteps:
    """A uniform current held at a sequence of speeds, for a stepped-current VIV study.

    Created by :func:`current_steps`. The current is held at ``speeds[k]`` for
    ``hold`` seconds, then ramps linearly over ``ramp`` seconds to the next
    speed, so one run continues each step from the state the previous one
    left (the response branch an up-and-down sweep follows). :meth:`write`
    gives the deck its MoorDyn-C ``Currents 2`` file.

    Attributes
    ----------
    speeds : numpy.ndarray
        Read-only ``(n_steps,)`` current speeds, m/s.
    hold, ramp : float
        Hold time of each speed and ramp time between speeds, s.
    direction : numpy.ndarray
        Read-only ``(3,)`` unit direction of the current.
    """

    speeds: npt.NDArray[Any]
    hold: float
    ramp: float
    direction: npt.NDArray[Any]

    @property
    def starts(self) -> npt.NDArray[Any]:
        """Read-only ``(n_steps,)`` start time of each hold, s."""
        return _readonly(np.arange(self.speeds.size) * (self.hold + self.ramp))

    @property
    def duration(self) -> float:
        """End of the last hold, s: the deck's ``TMax`` (a longer run repeats the record)."""
        return float(self.starts[-1] + self.hold)

    def windows(self, settle: float) -> tuple[tuple[float, float], ...]:
        """Return the analysis window of every hold, from ``settle`` seconds into it to its end."""
        if not (math.isfinite(float(settle)) and 0.0 <= float(settle) < self.hold):
            raise ValueError("settle must be finite, non-negative and shorter than the hold")
        return tuple((float(t) + float(settle), float(t) + self.hold) for t in self.starts)

    def write(
        self,
        folder: str | os.PathLike[str],
        *,
        depths: tuple[float, float] = (-1.0e4, 0.0),
        overwrite: bool = False,
    ) -> Path:
        """Write ``current_profile_dynamic.txt``, the deck's ``Currents 2`` file, into ``folder``.

        The current is uniform between the two ``depths`` (z, m, increasing) and
        held beyond them. The file's smallest time interval is the ramp, the time
        step MoorDyn-C and CableDyn resample it to, so the holds stay exact.

        Raises
        ------
        FileExistsError
            If the file exists and ``overwrite`` is false.
        ValueError
            If the depths are not two finite, increasing elevations.
        """
        low, high = (float(z) for z in depths)
        if not (math.isfinite(low) and math.isfinite(high) and low < high):
            raise ValueError("depths must be two finite, increasing elevations")
        target = Path(folder).expanduser().resolve() / "current_profile_dynamic.txt"
        if target.exists() and not overwrite:
            raise FileExistsError(f"output already exists: {target}")
        rows: list[str] = []
        for k, (start, speed) in enumerate(zip(self.starts, self.speeds, strict=True)):
            velocity = float(speed) * self.direction
            values = " ".join(f"{number(v)} {number(v)}" for v in velocity)
            times = (float(start), float(start) + self.hold) if k else (0.0, self.hold)
            rows.extend(f"{number(t)} {values}" for t in times)
        lines = [
            "--------------------- MoorDyn dynamic currents File -------------------------",
            "a uniform current held at a sequence of speeds (cabledyn.current_steps)",
            "line 5 lists the depths (m); each later row is t ux(z1 z2) uy(z1 z2) uz(z1 z2)",
            "------------------------------------------------------------------------------",
            f"{number(low)} {number(high)}",
            "time ux... uy... uz...",
            *rows,
        ]
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("\n".join(lines) + "\n", encoding="ascii")
        return target

    def responses(
        self,
        history: TimeHistory,
        *,
        diameter: float,
        settle: float,
        line_id: int | None = None,
    ) -> tuple[VIVResponse, ...]:
        """Return the :func:`viv_response` of every hold, from ``settle`` seconds into it."""
        return tuple(
            viv_response(
                history,
                flow=float(speed) * self.direction,
                diameter=diameter,
                start=start,
                stop=stop,
                line_id=line_id,
            )
            for speed, (start, stop) in zip(self.speeds, self.windows(settle), strict=True)
        )


def current_steps(
    speeds: npt.ArrayLike,
    *,
    hold: float,
    ramp: float,
    direction: npt.ArrayLike = (0.0, 1.0, 0.0),
) -> CurrentSteps:
    """Return a stepped current schedule.

    Parameters
    ----------
    speeds : array_like
        Positive current speeds, m/s, in the order they are held; an
        up-and-down sweep lists them up and back.
    hold : float
        Hold time of each speed, s; a whole number of ramps.
    ramp : float
        Positive ramp time between speeds, s.
    direction : array_like
        ``(3,)`` direction of the current (normalised), default ``+y``.

    Raises
    ------
    ValueError
        If a speed is not finite and positive, the times are invalid, the hold
        is not a whole number of ramps, or the direction is not a finite,
        non-zero 3-vector.
    """
    values = np.array(speeds, dtype=np.float64).reshape(-1)
    if values.size < 1 or not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("speeds must be one or more finite, positive values")
    hold_s, ramp_s = float(hold), float(ramp)
    if not (math.isfinite(hold_s) and hold_s > 0.0 and math.isfinite(ramp_s) and ramp_s > 0.0):
        raise ValueError("hold and ramp must be finite and positive")
    ratio = hold_s / ramp_s
    if abs(ratio - round(ratio)) > 1.0e-9 * ratio:
        raise ValueError("hold must be a whole number of ramps (the file's resampling step)")
    unit = np.array(direction, dtype=np.float64)
    if unit.shape != (3,) or not np.all(np.isfinite(unit)) or not np.any(unit):
        raise ValueError("direction must be a finite, non-zero 3-vector")
    return CurrentSteps(
        speeds=_readonly(values),
        hold=hold_s,
        ramp=ramp_s,
        direction=_readonly(unit / np.linalg.norm(unit)),
    )
