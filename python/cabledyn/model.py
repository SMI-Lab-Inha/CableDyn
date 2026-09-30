# SPDX-License-Identifier: Apache-2.0
"""The CableDyn model class: handle lifecycle + NumPy coupling boundary."""

from __future__ import annotations

import contextlib
import ctypes
import math
import numbers
import os
import threading
from collections.abc import Callable
from ctypes import POINTER, byref, c_bool, c_double, c_int, c_void_p
from pathlib import Path
from types import TracebackType
from typing import Any

import numpy as np
import numpy.typing as npt

from cabledyn import _lib as _lib_module
from cabledyn._lib import _lib
from cabledyn.errors import CableDynError, ConvergenceError
from cabledyn.objects import Body, Line, Point, Rod, build_views, eval_channel

__all__ = ["Body", "CableDyn", "CableDynError", "ConvergenceError", "Line", "Point", "Rod"]

FloatArray = npt.NDArray[np.float64]
ArrayInput = npt.ArrayLike

# The C library retains up to 1024 characters per handle; +1 for the NUL.
_MESSAGE_BUFFER = 1025


def _as_real_array(x: ArrayInput, name: str) -> FloatArray:
    """Return ``x`` as a contiguous ``float64`` array, rejecting non-real input.

    Complex values would otherwise lose their imaginary part silently, and strings
    or objects would be parsed or fail deep inside NumPy.
    """
    src = np.asarray(x)
    if src.dtype.kind == "c":
        raise TypeError(f"{name} must be real-valued; got a complex array")
    if src.dtype.kind not in "biuf":
        raise TypeError(f"{name} must be a real numeric array; got dtype {src.dtype}")
    return np.ascontiguousarray(src, dtype=np.float64)


def _as_real_scalar(value: object, name: str, *, positive: bool) -> float:
    """Return ``value`` as a finite ``float``; ``positive`` also requires ``> 0``."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise TypeError(f"{name} must be a real number; got {type(value).__name__}")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite; got {result}")
    if positive and result <= 0.0:
        raise ValueError(f"{name} must be positive; got {result}")
    if not positive and result < 0.0:
        raise ValueError(f"{name} must be non-negative; got {result}")
    return result


def _as_dof_array(x: ArrayInput, n: int, name: str) -> FloatArray:
    """Return ``x`` as a contiguous flat ``float64`` array of ``n`` values.

    Accepts a flat ``(n,)`` array or a column-per-point ``(3, n // 3)`` array.
    """
    arr = _as_real_array(x, name)
    if arr.ndim == 2 and arr.shape[0] == 3:
        arr = np.ascontiguousarray(arr.reshape(-1, order="F"))
    if arr.ndim != 1 or arr.size != n:
        raise ValueError(
            f"{name} must be flat ({n},) or column-per-point (3, {n // 3}); got shape {np.shape(x)}"
        )
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite (no NaN or Inf)")
    return arr


def _dbl_ptr(arr: FloatArray) -> Any:
    return arr.ctypes.data_as(POINTER(c_double))


def _encode_path(deck: str | os.PathLike[str]) -> bytes:
    """Encode a deck path for the Fortran runtime's narrow-character ``OPEN``."""
    text = os.fspath(deck)
    if "\0" in text:
        raise ValueError("deck path must not contain a NUL character")
    if os.name == "nt":
        # gfortran and ifx open files through the active ANSI code page.
        try:
            return text.encode("mbcs", errors="strict")
        except UnicodeEncodeError as exc:
            raise ValueError(
                f"deck path cannot be represented in the active Windows code page: {text!r}; "
                "move the deck to a path using representable characters"
            ) from exc
    return os.fsencode(text)


class CableDyn:
    """One CableDyn solver instance behind the C ABI (context-manager friendly).

    Parameters
    ----------
    deck : str | os.PathLike | None
        A MoorDyn-style sectioned ``.dat`` deck (the deck format in the
        documentation, https://cabledyn.readthedocs.io).
        When given, the static initial condition is solved immediately.
        ``None`` creates an unbound handle for a later :meth:`init_deck`.

    The coupled exchange follows the coupling boundary described in the
    documentation (https://cabledyn.readthedocs.io): kinematics in
    (position, velocity, and acceleration at the coupled points, 3 DOFs per
    point, metres and seconds), loads out (the force the lines exert at those
    DOFs, newtons). Arrays may be flat ``(n_dof,)`` with interleaved ``x, y, z``
    values, or column-per-point ``(3, n_points)``; results are always flat.
    A ``(3, 3)`` array is read column-per-point, so pass three points flat when
    in doubt.

    Calls on one instance are serialized by an internal lock. Independent
    instances may be used from different threads.

    The objects of an initialized model are available as :attr:`lines`,
    :attr:`points`, :attr:`bodies` and :attr:`rods` (see
    :mod:`cabledyn.objects`); any output channel can be evaluated with
    :meth:`channel`. Both read the committed state at the time of the call.
    """

    def __init__(self, deck: str | os.PathLike[str] | None = None) -> None:
        self._lock = threading.RLock()
        self._handle = c_void_p()
        self._lib = _lib
        # Bumped by every (re)initialization and by close(): object views compare it
        # to detect a stale handle.
        self._generation = 0
        self._views: tuple[tuple[Line, ...], tuple[Point, ...], tuple[Body, ...], tuple[Rod, ...]]
        self._views = ((), (), (), ())
        self._views_valid = False
        self._time = 0.0
        self._deck: Path | None = None
        status = c_int()
        _lib.CableDyn_Create(byref(self._handle), byref(status))
        if status.value != 0:
            raise CableDynError("CableDyn_Create", status.value, "")
        if deck is not None:
            try:
                self.init_deck(deck)
            except BaseException:
                self.close()
                raise

    # -- lifecycle ------------------------------------------------------------

    def close(self) -> None:
        """Release the handle (idempotent)."""
        with self._lock:
            if self._handle:
                status = c_int()
                _lib.CableDyn_Close(byref(self._handle), byref(status))
                self._handle = c_void_p()
            self._generation += 1
            self._views_valid = False

    def __enter__(self) -> CableDyn:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def __del__(self) -> None:  # best-effort; close() is the contract
        with contextlib.suppress(Exception):
            self.close()

    def _check(self, call: str, status: int) -> None:
        if status != 0:
            raise CableDynError(call, status, self._last_error())

    def _last_error(self) -> str:
        buf = ctypes.create_string_buffer(_MESSAGE_BUFFER)
        _lib.CableDyn_GetLastError(self._handle, buf, len(buf))
        return buf.value.decode("utf-8", errors="replace").strip()

    def last_error(self) -> str:
        """The library's last diagnostic message for this handle."""
        with self._lock:
            return self._last_error()

    # -- initialization -------------------------------------------------------

    def init_deck(self, deck: str | os.PathLike[str]) -> None:
        """Initialize from a deck file: parse, build, and solve the static IC.

        A failed initialization keeps any previously initialized model.
        """
        _encode_path(deck)  # rejects a NUL before any file-system call
        # the deck path as the rest of the package reads it ('~' expanded)
        resolved = Path(deck).expanduser().resolve()
        path = _encode_path(resolved)
        with self._lock:
            status = c_int()
            _lib.CableDyn_InitDeck(self._handle, path, len(path), byref(status))
            self._check("CableDyn_InitDeck", status.value)
            self._generation += 1
            self._views_valid = False
            self._time = 0.0
            self._deck = resolved

    @property
    def initialized(self) -> bool:
        """Whether the handle holds an initialized model."""
        with self._lock:
            return bool(_lib.CableDyn_IsInitialized(self._handle))

    # -- structure queries ----------------------------------------------------

    def _query(self, fn: Callable[..., int], call: str) -> int:
        with self._lock:
            status = c_int()
            n = fn(self._handle, byref(status))
            self._check(call, status.value)
            return int(n)

    @property
    def n_coupled_dof(self) -> int:
        """Coupled DOF count (3 per host-driven point)."""
        return self._query(_lib.CableDyn_NCoupledDOF, "CableDyn_NCoupledDOF")

    @property
    def n_points(self) -> int:
        """Stored system point count (the fluid-field column count)."""
        return self._query(_lib.CableDyn_NPoints, "CableDyn_NPoints")

    @property
    def n_lines(self) -> int:
        """Owned line count."""
        return self._query(_lib.CableDyn_NLines, "CableDyn_NLines")

    @property
    def deck(self) -> Path | None:
        """Absolute path of the deck the model was initialized from, if any."""
        with self._lock:
            return self._deck

    @property
    def time(self) -> float:
        """Simulation time [s]: 0 after initialization, advanced by every completed step."""
        with self._lock:
            return self._time

    # -- object queries -------------------------------------------------------

    def _object_views(
        self,
    ) -> tuple[tuple[Line, ...], tuple[Point, ...], tuple[Body, ...], tuple[Rod, ...]]:
        with self._lock:
            if not self._views_valid:
                if _lib_module.abi_minor() < _lib_module.OBJECT_QUERY_ABI_MINOR:
                    # not a C status: the loaded library lacks the entry points
                    raise RuntimeError(
                        "the loaded CableDyn library predates the object queries "
                        "(ABI 1 minor extension 1); rebuild it from this checkout"
                    )
                self._views = build_views(self)
                self._views_valid = True
            return self._views

    @property
    def lines(self) -> tuple[Line, ...]:
        """The model's lines in ascending deck id (see :class:`~cabledyn.objects.Line`)."""
        return self._object_views()[0]

    @property
    def points(self) -> tuple[Point, ...]:
        """The model's points in solver order (see :class:`~cabledyn.objects.Point`)."""
        return self._object_views()[1]

    @property
    def bodies(self) -> tuple[Body, ...]:
        """The model's Rigid6 bodies (see :class:`~cabledyn.objects.Body`)."""
        return self._object_views()[2]

    @property
    def rods(self) -> tuple[Rod, ...]:
        """The model's rods (see :class:`~cabledyn.objects.Rod`)."""
        return self._object_views()[3]

    @staticmethod
    def _by_id(items: tuple[Any, ...], object_id: int, kind: str) -> Any:
        for item in items:
            if item.id == object_id:
                return item
        raise KeyError(f"no {kind} with deck id {object_id}")

    def line(self, line_id: int) -> Line:
        """The line with deck id ``line_id``."""
        result: Line = self._by_id(self.lines, line_id, "line")
        return result

    def point(self, point_id: int) -> Point:
        """The point with deck id ``point_id``."""
        result: Point = self._by_id(self.points, point_id, "point")
        return result

    def body(self, body_id: int) -> Body:
        """The body with deck id ``body_id``."""
        result: Body = self._by_id(self.bodies, body_id, "body")
        return result

    def rod(self, rod_id: int) -> Rod:
        """The rod with deck id ``rod_id``."""
        result: Rod = self._by_id(self.rods, rod_id, "rod")
        return result

    def channel(self, token: str) -> float:
        """Evaluate one output channel at the committed state.

        ``token`` uses the deck OUTPUTS vocabulary (``FairTen1``, ``L2N5px``,
        ``Curv1N3``, ``Point4Fz``, ``Body1Pz``, ``Rod2TenA``, ...) and returns the
        value the standalone driver writes for that channel at the same state.
        """
        with self._lock:
            self._object_views()
            return eval_channel(self, token)

    def channels(self, tokens: list[str] | tuple[str, ...]) -> FloatArray:
        """Evaluate several output channels; returns ``(len(tokens),)``."""
        with self._lock:
            return np.array([self.channel(token) for token in tokens], dtype=np.float64)

    def step_held(self, dt: float) -> int:
        """Advance one step with every coupled point held at its current position.

        The coupled kinematics are the current positions with zero velocity and
        acceleration, the standalone driver's treatment of a deck without vessel
        motion. Returns the Newton iteration count, as :meth:`step`.
        """
        with self._lock:
            q, _, _ = self.get_coupled_motion()
            zero = np.zeros_like(q)
            return self.step(dt, q, zero, zero)

    # -- the coupling boundary ------------------------------------------------

    def get_coupled_motion(self) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Current coupled kinematics ``(q, v, a)``, each shaped ``(n_dof,)``."""
        with self._lock:
            n = self.n_coupled_dof
            q = np.zeros(n)
            v = np.zeros(n)
            a = np.zeros(n)
            status = c_int()
            _lib.CableDyn_GetCoupledMotion(
                self._handle, _dbl_ptr(q), _dbl_ptr(v), _dbl_ptr(a), n, byref(status)
            )
            self._check("CableDyn_GetCoupledMotion", status.value)
            return q, v, a

    def update_states(self, q: ArrayInput, v: ArrayInput, a: ArrayInput) -> None:
        """Transfer coupled kinematics without advancing time (direct feedthrough)."""
        with self._lock:
            n = self.n_coupled_dof
            qa = _as_dof_array(q, n, "q")
            va = _as_dof_array(v, n, "v")
            aa = _as_dof_array(a, n, "a")
            status = c_int()
            _lib.CableDyn_UpdateStates(
                self._handle, _dbl_ptr(qa), _dbl_ptr(va), _dbl_ptr(aa), n, byref(status)
            )
            self._check("CableDyn_UpdateStates", status.value)

    def step(self, dt: float, q: ArrayInput, v: ArrayInput, a: ArrayInput) -> int:
        """Advance one implicit step of ``dt`` with the coupled kinematics at t+dt.

        Returns the Newton iteration count. A failed step raises
        :class:`CableDynError` and leaves the state at ``t``. A step that ends
        without meeting the Newton tolerance raises :class:`ConvergenceError`;
        the model has then already advanced to ``t + dt`` with its best iterate.
        ``dt`` must be a finite positive real number (a string, complex, or
        boolean value raises :class:`TypeError`).
        """
        step_dt = _as_real_scalar(dt, "dt", positive=True)
        with self._lock:
            n = self.n_coupled_dof
            qa = _as_dof_array(q, n, "q")
            va = _as_dof_array(v, n, "v")
            aa = _as_dof_array(a, n, "a")
            conv = c_bool()
            stall = c_bool()
            nit = c_int()
            status = c_int()
            _lib.CableDyn_Step(
                self._handle,
                step_dt,
                _dbl_ptr(qa),
                _dbl_ptr(va),
                _dbl_ptr(aa),
                n,
                byref(conv),
                byref(stall),
                byref(nit),
                byref(status),
            )
            self._check("CableDyn_Step", status.value)
            self._time += step_dt
            if not conv.value:
                detail = self._last_error() or (
                    f"Newton tolerance not met after {nit.value} iterations"
                    + (" (line search stalled)" if stall.value else "")
                )
                raise ConvergenceError(
                    "CableDyn_Step", detail, n_iter=int(nit.value), stalled=bool(stall.value)
                )
            return int(nit.value)

    def calc_output(self) -> FloatArray:
        """The coupled reaction loads at the current committed state, ``(n_dof,)`` N."""
        with self._lock:
            n = self.n_coupled_dof
            loads = np.zeros(n)
            status = c_int()
            _lib.CableDyn_CalcOutput(self._handle, _dbl_ptr(loads), n, byref(status))
            self._check("CableDyn_CalcOutput", status.value)
            return loads

    def update_point_fluid_fields(
        self,
        fluid_velocity: ArrayInput,
        fluid_acceleration: ArrayInput,
        waterline_z: ArrayInput,
        fluid_density: float,
    ) -> None:
        """Prescribe external per-point fluid kinematics (the CFD-coupling surface).

        ``fluid_velocity`` and ``fluid_acceleration`` are ``(3, n_points)`` or
        flat ``(3 * n_points,)``; ``waterline_z`` is ``(n_points,)``;
        ``fluid_density`` is a finite non-negative real number.
        """
        density = _as_real_scalar(fluid_density, "fluid_density", positive=False)
        with self._lock:
            npt_ = self.n_points
            fv = _as_dof_array(fluid_velocity, 3 * npt_, "fluid_velocity")
            fa = _as_dof_array(fluid_acceleration, 3 * npt_, "fluid_acceleration")
            wl = _as_real_array(waterline_z, "waterline_z")
            if wl.ndim != 1 or wl.size != npt_:
                raise ValueError(f"waterline_z must have {npt_} entries")
            if not np.all(np.isfinite(wl)):
                raise ValueError("waterline_z must be finite (no NaN or Inf)")
            status = c_int()
            _lib.CableDyn_UpdatePointFluidFields(
                self._handle,
                _dbl_ptr(fv),
                _dbl_ptr(fa),
                _dbl_ptr(wl),
                npt_,
                density,
                byref(status),
            )
            self._check("CableDyn_UpdatePointFluidFields", status.value)
