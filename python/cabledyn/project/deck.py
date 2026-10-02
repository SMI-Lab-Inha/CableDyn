# SPDX-License-Identifier: Apache-2.0
"""Deck adapters: read a CableDyn deck into a project, and write a project as a deck.

Both directions go through :class:`cabledyn.DeckModel` and
:class:`cabledyn.DeckFile`, so the deck grammar and the native deck rules
live in one place. :class:`DeckReader` builds the object graph from a
:class:`~cabledyn.DeckModel`; :class:`DeckWriter` builds a fresh
:class:`~cabledyn.DeckModel` from the graph and renders it with
:meth:`~cabledyn.DeckModel.to_text`. The writer returns an :class:`IdMap`
between objects and deck ids.

Deck ids, keyword spellings, option order and option descriptions are kept
as *hints* (``deck_hints``) on the objects, so reading a deck and writing it
back gives :meth:`DeckModel.load(deck).to_text() <cabledyn.DeckModel.to_text>`
exactly; deck comments and column alignment are not kept. Option rows the
object model does not type are kept verbatim as
:class:`~cabledyn.project.ExtraOption` rows.

The deck vocabulary is kept as it is: ``Vessel`` points and bodies, and the
``vesselMotion``, ``vesselRAO`` and ``vesselRef`` options, map to the
floater classes and properties of the object model.
"""

from __future__ import annotations

import os
import re
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from cabledyn import builder as _b
from cabledyn.builder import DeckModel, DeckReferenceError
from cabledyn.deck_file import (
    _ROD_ATTACHMENT_ALIASES,
    DeckFile,
    _native_int,
    _option_key,
)
from cabledyn.errors import DeckFormatError
from cabledyn.project.base import ModelObject
from cabledyn.project.bodies import Body, Floater, Point3Body, Rigid6Body, Turbine
from cabledyn.project.environment import (
    BathymetrySeabed,
    CurrentModel,
    FloaterMotionRecord,
    FloaterRAO,
    JonswapWave,
    MotionFile,
    MotionSource,
    MultiTrainSea,
    NoCurrent,
    NoMotion,
    NoWaves,
    OchiHubbleWave,
    ProfileCurrent,
    RegularWave,
    SpectrumWave,
    UniformCurrent,
    WaveModel,
    WaveTrain,
)
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.lines import (
    Attachment,
    BuoyancyModule,
    ClumpWeight,
    EndConnection,
    Line,
    Section,
)
from cabledyn.project.points import (
    BodyPoint,
    ConnectPoint,
    FixedPoint,
    FloaterPoint,
    FreePoint,
    LineEndTarget,
    Point,
    RodEndpoint,
    RodPoint,
    TurbinePoint,
)
from cabledyn.project.project import Project
from cabledyn.project.rods import Rod
from cabledyn.project.rows import (
    Control,
    EquivalentBuoyancy,
    ExternalLoad,
    Failure,
    NamedChannel,
    ObjectChannel,
    SyropeIC,
)
from cabledyn.project.settings import ExtraOption
from cabledyn.project.types import (
    BendingModel,
    GenericLineType,
    LinearAxial,
    LineType,
    RodType,
    SyropeAxial,
    ViscoelasticAxial,
)

__all__ = ["DeckExportError", "DeckReader", "DeckWriter", "IdMap", "typed_option"]

_TRUE = frozenset({"true", "t", "yes", "y", "on", "1"})
_FALSE = frozenset({"false", "f", "no", "n", "off", "0"})
_DISABLED = frozenset({"0", "none"})
_PINNED = frozenset({"pinned", "free", "zero"})
_RIGID = frozenset({"rigid", "infinity", "inf"})
_SPECTRA = {
    "pm": "pm",
    "piersonmoskowitz": "pm",
    "pierson-moskowitz": "pm",
    "issc": "issc",
    "bretschneider": "bretschneider",
    "torsethaugen": "torsethaugen",
}
_OCHI = frozenset({"ochihubble", "ochi-hubble", "ochi_hubble"})
_MOTION_KEYS = {
    "motionfile": MotionFile,
    "vesselmotion": FloaterMotionRecord,
    "vesselrao": FloaterRAO,
}
_MOTION_KEYWORDS = {MotionFile: "motionFile", FloaterMotionRecord: "vesselMotion"}
_SEA_WATER_DENSITY = 1025.0


class DeckExportError(ValueError):
    """The project cannot be written as a deck.

    Attributes
    ----------
    obj : ModelObject | None
        The object that could not be written, if known.
    """

    def __init__(self, message: str, obj: ModelObject | None = None) -> None:
        super().__init__(message)
        self.obj = obj


# --------------------------------------------------------------------------- option codecs


def _fnum(token: str) -> float:
    return _b._f(token)


def _num(value: float) -> str:
    return repr(float(value))


@dataclass(frozen=True)
class _Scalar:
    identity: str
    owner: str
    attr: str
    codec: str
    keyword: str


_SCALARS: tuple[_Scalar, ...] = (
    _Scalar("gravity", "environment", "gravity", "float", "g"),
    _Scalar("water_density", "environment", "water_density", "float", "rhoW"),
    _Scalar("water_depth", "environment", "water_depth", "float", "WtrDpth"),
    _Scalar("kbot", "seabed", "stiffness", "float", "kBot"),
    _Scalar("cbot", "seabed", "damping", "float", "cBot"),
    _Scalar("friction", "seabed", "friction", "friction", "frictionMu"),
    _Scalar("friction_axial", "seabed", "friction_axial", "float", "frictionMuAxial"),
    _Scalar("friction_lateral", "seabed", "friction_lateral", "float", "frictionMuLateral"),
    _Scalar("bodyic", "settings", "body_initial_condition", "text", "bodyIC"),
    _Scalar("icmode", "settings", "initial_condition_mode", "text", "ICmode"),
    _Scalar("bodywetting", "settings", "body_wetting", "text", "bodyWetting"),
    _Scalar("bodyhydro", "settings", "body_hydrodynamics", "text", "bodyHydro"),
    _Scalar("rodhydro", "settings", "rod_hydrodynamics", "text", "rodHydro"),
    _Scalar("bodyscheme", "settings", "body_scheme", "text", "bodyScheme"),
    _Scalar("bodysubstep", "settings", "body_substep", "text", "bodySubstep"),
    _Scalar("rhoinf", "settings", "spectral_radius", "float", "rhoInf"),
    _Scalar("max_strain", "settings", "max_strain", "float", "maxStrain"),
    _Scalar("modified_newton", "settings", "modified_newton", "bool", "modified_newton"),
    _Scalar(
        "cable_load_feedback", "settings", "cable_load_feedback", "bool", "cable_load_feedback"
    ),
    _Scalar("adaptive_mesh", "settings", "adaptive_mesh", "bool", "adaptive_mesh"),
    _Scalar("cable_statics", "settings", "cable_statics", "text", "cable_statics"),
    _Scalar("alpha_force_blend", "settings", "alpha_force_blend", "bool", "alpha_force_blend"),
    _Scalar("tensile_safety", "settings", "tensile_safety", "text", "tensile_safety"),
    _Scalar(
        "tensile_strain_tolerance",
        "settings",
        "tensile_strain_tolerance",
        "float",
        "tensile_strain_tolerance",
    ),
    _Scalar(
        "recovery_max_substeps", "settings", "recovery_max_substeps", "int", "recovery_max_substeps"
    ),
    _Scalar(
        "axial_quadrature_order",
        "settings",
        "axial_quadrature_order",
        "int",
        "axial_quadrature_order",
    ),
    _Scalar(
        "bending_quadrature_order",
        "settings",
        "bending_quadrature_order",
        "int",
        "bending_quadrature_order",
    ),
    _Scalar("dtm", "settings", "time_step", "float", "dtM"),
    _Scalar("tmax", "settings", "duration", "float", "TMax"),
    _Scalar("range_start", "settings", "range_start", "float", "RangeStart"),
    _Scalar("n_modes", "settings", "mode_count", "int", "nModes"),
    _Scalar("tscheme", "settings", "time_scheme", "text", "tScheme"),
    _Scalar("dtic", "settings", "ic_time_step", "float", "dtIC"),
    _Scalar("tmaxic", "settings", "ic_duration", "float", "TmaxIC"),
    _Scalar("cdscaleic", "settings", "ic_drag_scale", "float", "CdScaleIC"),
    _Scalar("threshic", "settings", "ic_threshold", "float", "threshIC"),
    _Scalar("writelog", "settings", "write_log", "int", "WriteLog"),
    _Scalar("dtout", "settings", "output_interval", "float", "dtOut"),
    _Scalar("ramp_time", "environment", "ramp_time", "float", "rampTime"),
    _Scalar("wave_seed", "environment", "wave_seed", "int", "WaveSeed"),
    _Scalar("wave_spreading", "environment", "wave_spreading", "float", "WaveSpreading"),
    _Scalar("wave_directions", "environment", "wave_directions", "int", "WaveDirections"),
    _Scalar("wave_components", "environment", "wave_components", "int", "WaveComponents"),
    _Scalar("stream_order", "environment", "stream_order", "int", "StreamOrder"),
    _Scalar("dtwave", "environment", "wave_time_step", "float", "dtWave"),
    _Scalar("waterkin", "environment", "water_kinematics", "text", "WaterKin"),
    _Scalar("currents", "environment", "current_mode", "int", "Currents"),
)
_SCALAR_BY_ID = {item.identity: item for item in _SCALARS}


def _decode_scalar(codec: str, values: tuple[str, ...]) -> Any:
    if len(values) != 1:
        raise ValueError("expected one value")
    token = values[0]
    if codec == "float":
        return _fnum(token)
    if codec == "friction":
        return 0.0 if token.lower() == "none" else _fnum(token)
    if codec == "int":
        return _native_int(token)
    if codec == "bool":
        lowered = token.lower()
        if lowered in _TRUE:
            return True
        if lowered in _FALSE:
            return False
        raise ValueError(f"{token!r} is not a logical value")
    return token


def _encode_scalar(codec: str, value: Any) -> tuple[str, ...]:
    if codec in {"float", "friction"}:
        return (_num(value),)
    if codec == "int":
        return (str(int(value)),)
    if codec == "bool":
        return ("True" if value else "False",)
    return (str(value),)


def _floats(values: Iterable[str]) -> list[float]:
    return [_fnum(token) for token in values]


def _decode_wave(values: tuple[str, ...], train: bool) -> WaveModel:
    mode = values[0].lower()
    numbers = _floats(values[1:])
    if mode == "none" and not numbers and not train:
        return NoWaves()
    if mode in {"airy", "stream", "dean"} and len(numbers) == 3:
        if train and mode != "airy":
            raise ValueError("only airy regular trains")
        return RegularWave(theory=mode, height=numbers[0], period=numbers[1], direction=numbers[2])
    if mode == "jonswap" and len(numbers) == 4:
        return JonswapWave(
            significant_height=numbers[0],
            peak_period=numbers[1],
            peak_enhancement=numbers[2],
            direction=numbers[3],
        )
    if mode in _SPECTRA and len(numbers) == 3:
        return SpectrumWave(
            spectrum=_SPECTRA[mode],
            significant_height=numbers[0],
            peak_period=numbers[1],
            direction=numbers[2],
        )
    if mode in _OCHI and len(numbers) == 7:
        return OchiHubbleWave(
            height_1=numbers[0],
            period_1=numbers[1],
            shape_1=numbers[2],
            height_2=numbers[3],
            period_2=numbers[4],
            shape_2=numbers[5],
            direction=numbers[6],
        )
    raise ValueError(f"unrecognised wave row {values!r}")


def _encode_wave(wave: WaveModel) -> tuple[str, ...]:
    if isinstance(wave, NoWaves):
        return ("none",)
    if isinstance(wave, RegularWave):
        return (wave.theory, _num(wave.height), _num(wave.period), _num(wave.direction))
    if isinstance(wave, JonswapWave):
        return (
            "jonswap",
            _num(wave.significant_height),
            _num(wave.peak_period),
            _num(wave.peak_enhancement),
            _num(wave.direction),
        )
    if isinstance(wave, SpectrumWave):
        return (
            wave.spectrum,
            _num(wave.significant_height),
            _num(wave.peak_period),
            _num(wave.direction),
        )
    if isinstance(wave, OchiHubbleWave):
        values = (
            wave.height_1,
            wave.period_1,
            wave.shape_1,
            wave.height_2,
            wave.period_2,
            wave.shape_2,
            wave.direction,
        )
        return ("ochihubble", *(_num(value) for value in values))
    raise DeckExportError(f"{wave.label()} cannot be written as a waves row", wave)


def _decode_train(values: tuple[str, ...]) -> WaveTrain:
    mode = values[0].lower()
    if mode == "airy":
        return WaveTrain(wave=_decode_wave(values, True))
    *head, spread = values
    train = WaveTrain(wave=_decode_wave(tuple(head), True))
    train.spreading = _fnum(spread)
    return train


def _encode_train(train: WaveTrain) -> tuple[str, ...]:
    tokens = _encode_wave(train.wave)
    if isinstance(train.wave, RegularWave):
        return tokens
    return (*tokens, _num(train.spreading or 0.0))


def _decode_current(values: tuple[str, ...]) -> CurrentModel:
    mode = values[0].lower()
    numbers = _floats(values[1:])
    if mode == "none" and not numbers:
        return NoCurrent()
    if mode == "uniform" and len(numbers) == 3:
        return UniformCurrent(velocity=numbers)
    if mode == "profile" and len(numbers) == 8:
        return ProfileCurrent(
            z_1=numbers[0], velocity_1=numbers[1:4], z_2=numbers[4], velocity_2=numbers[5:8]
        )
    raise ValueError(f"unrecognised current row {values!r}")


def _encode_current(current: CurrentModel) -> tuple[str, ...]:
    if isinstance(current, UniformCurrent):
        return ("uniform", *(_num(v) for v in current.velocity))
    if isinstance(current, ProfileCurrent):
        return (
            "profile",
            _num(current.z_1),
            *(_num(v) for v in current.velocity_1),
            _num(current.z_2),
            *(_num(v) for v in current.velocity_2),
        )
    if isinstance(current, NoCurrent):
        return ("none",)
    raise DeckExportError(f"{current.label()} cannot be written as a current row", current)


def _decode_vector(values: tuple[str, ...]) -> tuple[float, float, float]:
    parts = "|".join(values).split("|")
    if len(parts) != 3:
        raise ValueError("expected x|y|z")
    x, y, z = _floats(parts)
    return (x, y, z)


def _encode_vector(value: tuple[float, float, float]) -> tuple[str, ...]:
    return ("|".join(_num(v) for v in value),)


_SOLVER = (
    ("solver_relative_tolerance", "float"),
    ("solver_absolute_tolerance", "float"),
    ("solver_max_iterations", "int"),
    ("solver_backtracks", "int"),
)


# --------------------------------------------------------------------------- id map


@dataclass
class IdMap:
    """Correspondence between objects and deck ids, names and text rows.

    Attributes
    ----------
    keys : dict[str, tuple[str, int | str]]
        Object uid to ``(kind, deck id or name)``; kinds are ``line_type``,
        ``rod_type``, ``body``, ``rod``, ``turbine``, ``point``, ``line``
        and ``external_load``.
    objects : dict[tuple[str, int | str], ModelObject]
        The reverse mapping (names in lower case).
    rows : dict[int, ModelObject]
        Zero-based line number of the written deck text to the object that
        produced the row (filled by :meth:`DeckWriter.to_text`).
    """

    keys: dict[str, tuple[str, int | str]] = field(default_factory=dict)
    objects: dict[tuple[str, int | str], ModelObject] = field(default_factory=dict)
    rows: dict[int, ModelObject] = field(default_factory=dict)

    def add(self, obj: ModelObject, kind: str, key: int | str) -> None:
        """Record that ``obj`` is written as ``kind`` ``key``."""
        self.keys[obj.uid] = (kind, key)
        self.objects[(kind, key.lower() if isinstance(key, str) else key)] = obj

    def deck_id(self, obj: ModelObject) -> int | str | None:
        """Return the deck id (or name) of ``obj``, or ``None``."""
        found = self.keys.get(obj.uid)
        return None if found is None else found[1]

    def object_for(self, kind: str, key: int | str) -> ModelObject | None:
        """Return the object written as ``kind`` ``key`` (names case-insensitive)."""
        return self.objects.get((kind, key.lower() if isinstance(key, str) else key))

    def object_at_row(self, line: int) -> ModelObject | None:
        """Return the object behind zero-based deck text line ``line``, if any."""
        return self.rows.get(line)


# --------------------------------------------------------------------------- reader


class DeckReader:
    """Build a :class:`~cabledyn.project.Project` from a CableDyn deck.

    Examples
    --------
    >>> from cabledyn.project import DeckReader
    >>> project = DeckReader.read("examples/spread_3line_chain.dat")  # doctest: +SKIP
    """

    @classmethod
    def read(cls, path: str | os.PathLike[str], *, caller_driven: bool = False) -> Project:
        """Read and validate a deck file.

        Parameters
        ----------
        path : str | os.PathLike
            Deck file.
        caller_driven : bool
            Validate for the OpenFAST coupling instead of the standalone driver.

        Returns
        -------
        Project
            The project.

        Raises
        ------
        DeckFormatError
            If the deck violates the native deck contract.
        """
        return cls.from_model(DeckModel.load(path, caller_driven=caller_driven))

    @classmethod
    def from_text(
        cls,
        text: str,
        *,
        path: str | os.PathLike[str] = "deck.dat",
        caller_driven: bool = False,
    ) -> Project:
        """Validate deck text and return its project (see :meth:`read`)."""
        return cls.from_model(DeckModel.from_text(text, path=path, caller_driven=caller_driven))

    @classmethod
    def from_model(cls, model: DeckModel) -> Project:
        """Build a project from a :class:`~cabledyn.DeckModel`.

        Parameters
        ----------
        model : DeckModel
            The deck model (it is not modified).

        Returns
        -------
        Project
            The project, with deck ids and spellings kept as hints.
        """
        return _Reader(model).run()


class _Reader:
    def __init__(self, model: DeckModel) -> None:
        self.model = model
        self.project = Project()
        self.line_types: dict[int, LineType] = {}
        self.rod_types: dict[int, RodType] = {}
        self.bodies: dict[int, Body] = {}
        self.rods: dict[int, Rod] = {}
        self.points: dict[int, Point] = {}
        self.lines: dict[int, Line] = {}

    def run(self) -> Project:
        model, project = self.model, self.project
        project.name = model.path.stem
        project.title = model.title
        project.caller_driven = model.caller_driven
        project.deck_path = str(model.path)
        for line_type in model.line_types:
            self.line_type(line_type)
        for rod_type in model.rod_types:
            self.rod_type(rod_type)
        for body in model.bodies:
            self.body(body)
        for rod in model.rods:
            self.rod(rod)
        for turbine in model.turbines:
            project.turbines.append(
                Turbine(
                    f"Turbine {turbine.id}",
                    number=turbine.id,
                    position=(turbine.x, turbine.y, turbine.z),
                    platform_displacement=turbine.ptfm or (),
                )
            )
        for point in model.points:
            self.point(point)
        for line in model.lines:
            self.line(line)
        for index, connection in enumerate(model.end_connections):
            self.end_connection(connection, index)
        for eq in model.equivalent_buoyancy:
            project.equivalent_buoyancy.append(
                EquivalentBuoyancy(
                    f"Equivalent {eq.line_type.name}",
                    line_type=self.line_types[id(eq.line_type)],
                    diameter=eq.diam,
                    submerged_weight=eq.submerged_weight,
                )
            )
        for index, attachment in enumerate(model.attachments):
            self.attachment(attachment, index)
        for history in model.syrope_ic:
            project.syrope_history.append(
                SyropeIC(
                    lines=self._lines(history.lines),
                    max_tension=history.tmax0,
                    mean_tension=history.tmean0,
                )
            )
        for number, failure in enumerate(model.failures, start=1):
            project.failures.append(
                Failure(
                    f"Failure {number}",
                    point=self.points[id(failure.point)],
                    lines=self._lines(failure.lines),
                    time=failure.fail_time,
                    tension=failure.fail_tension,
                )
            )
        for control in model.controls:
            project.controls.append(
                Control(
                    f"Control {control.channel}",
                    channel=control.channel,
                    lines=self._lines(control.lines),
                )
            )
        for load in model.external_loads:
            self.external_load(load)
        self.options()
        self.outputs()
        return project

    # ------------------------------------------------------------------ libraries

    def _lines(self, rows: Iterable[_b.Line]) -> list[Line]:
        return [self.lines[id(line)] for line in rows]

    def line_type(self, row: _b.LineType) -> None:
        axial: LinearAxial | ViscoelasticAxial | SyropeAxial
        if isinstance(row.ea, _b.SyropeEA):
            axial = SyropeAxial(settings_file=row.ea.settings, alpha=row.ea.alpha, beta=row.ea.beta)
        elif isinstance(row.ea, tuple):
            if len(row.ea) == 2:
                axial = ViscoelasticAxial(static_stiffness=row.ea[0], dynamic_stiffness=row.ea[1])
            elif len(row.ea) == 3:
                axial = ViscoelasticAxial(
                    static_stiffness=row.ea[0], alpha_mbl=row.ea[1], beta=row.ea[2]
                )
            else:
                raise DeckFormatError(f"line type {row.name!r}: unsupported EA {row.ea!r}")
        else:
            axial = LinearAxial(stiffness=row.ea)
        if isinstance(row.ba, tuple) and len(row.ba) > 2:
            raise DeckFormatError(f"line type {row.name!r}: unsupported BA {row.ba!r}")
        axial.set_deck_damping(row.ba)
        bending = BendingModel(
            bending_stiffness=row.ei,
            shear_stiffness=row.gas,
            torsional_stiffness=row.gj,
            rotary_inertia_axial=row.irt,
            rotary_inertia_normal=row.irn,
        )
        item = GenericLineType(
            row.name,
            diameter=row.diam,
            mass_per_length=row.mass,
            axial=axial,
            bending=bending,
            drag_normal=row.cdn,
            drag_axial=row.cdt,
            added_mass_normal=row.can,
            added_mass_axial=row.cat,
        )
        self.project.line_types.append(item)
        self.line_types[id(row)] = item

    def rod_type(self, row: _b.RodType) -> None:
        item = RodType(
            row.name,
            diameter=row.diam,
            mass_per_length=row.mass,
            drag=row.cd,
            added_mass=row.ca,
            end_drag=row.cd_end,
            end_added_mass=row.ca_end,
            axial_drag=row.cd_ax,
            axial_added_mass=row.ca_ax,
        )
        self.project.rod_types.append(item)
        self.rod_types[id(row)] = item

    # ------------------------------------------------------------------ objects

    def body(self, row: _b.Body | _b.MoorDynBody) -> None:
        kind = row.type.lower()
        item: Body
        if kind in {"rigid6", "free"}:
            item = Rigid6Body(f"Body {row.id}")
        elif kind == "point3":
            item = Point3Body(f"Body {row.id}")
        elif kind in {"coupled", "cpld", "vessel", "ves"}:
            item = Floater(f"Body {row.id}")
            item.kind = "coupled" if kind in {"coupled", "cpld"} else "floater"
        else:
            raise DeckFormatError(f"body {row.id}: unsupported type {row.type!r}")
        item.position = (row.x, row.y, row.z)
        item.orientation = (row.roll, row.pitch, row.yaw)
        item.mass = row.mass
        item.volume = row.volume
        if isinstance(row, _b.MoorDynBody):
            item.row_format = "moordyn"
            item.centre_of_gravity = _tuple(row.cg)
            item.inertia = _tuple(row.inertia)
            item.drag_area = _tuple(row.cda)
            item.added_mass = _tuple(row.ca)
        else:
            item.row_format = "cabledyn"
            item.heave_stiffness = row.c33
            item.roll_stiffness = row.c44
            item.pitch_stiffness = row.c55
            item.drag_area = (row.cda,)
            item.added_mass = (row.ca,)
            item.inertia = row.inertia or ()
        item.deck_hints.update(deck_id=row.id, deck_type=row.type)
        self.project.bodies.append(item)
        self.bodies[id(row)] = item

    def rod(self, row: _b.Rod) -> None:
        kind = _ROD_ATTACHMENT_ALIASES.get(row.type.lower(), row.type.lower())
        attachment = {
            "vessel": "floater",
            "bodypinned": "body_pinned",
            "bodypin": "body_pinned",
        }.get(kind, kind)
        item = Rod(
            f"Rod {row.id}",
            rod_type=self.rod_types[id(row.rod_type)],
            attachment=attachment,
            body=None if row.body is None else self.bodies[id(row.body)],
            position_a=row.end_a,
            position_b=row.end_b,
            segments=row.num_segs,
            outputs=row.outputs,
        )
        item.deck_hints.update(deck_id=row.id, deck_type=row.type)
        self.project.rods.append(item)
        self.rods[id(row)] = item

    def _endpoint(self, end: _b.RodEnd) -> RodEndpoint:
        return self.rods[id(end.rod)].endpoint(end.end)

    def point(self, row: _b.Point) -> None:
        kind = row.type.lower()
        item: Point
        if kind == "fixed":
            item = FixedPoint()
        elif kind in {"coupled", "vessel"}:
            item = FloaterPoint(kind="coupled" if kind == "coupled" else "floater")
        elif kind == "free":
            item = FreePoint()
        elif kind == "connect":
            item = ConnectPoint()
        elif kind == "body" and row.body is not None:
            item = BodyPoint(body=self.bodies[id(row.body)])
        elif kind == "rod" and row.rod_end is not None:
            item = RodPoint(rod_end=self._endpoint(row.rod_end))
        elif kind in {"turbine", "t"} and row.turbine is not None:
            farm = {turbine.number: turbine for turbine in self.project.turbines}
            if row.turbine in farm:
                item = TurbinePoint(turbine=farm[row.turbine])
            else:
                item = TurbinePoint(turbine_number=row.turbine)
        else:
            raise DeckFormatError(f"point {row.id}: unsupported type {row.type!r}")
        item.name = f"Point {row.id}"
        item.position = (row.x, row.y, row.z)
        item.mass = row.mass
        item.volume = row.volume
        item.drag_area = row.cda
        item.added_mass = row.ca
        item.deck_hints.update(deck_id=row.id, deck_type=row.type)
        self.project.points.append(item)
        self.points[id(row)] = item

    def _target(self, end: _b.Point | _b.RodEnd) -> LineEndTarget:
        if isinstance(end, _b.RodEnd):
            return self._endpoint(end)
        return self.points[id(end)]

    def line(self, row: _b.Line) -> None:
        item = Line(
            f"Line {row.id}",
            end_a=self._target(row.end_a),
            end_b=self._target(row.end_b),
            outputs=row.outputs,
        )
        for index, section in enumerate(row.sections, start=1):
            item.sections.append(
                Section(
                    f"Section {index}",
                    line_type=self.line_types[id(section.line_type)],
                    length=section.length,
                    segments=section.num_segs,
                )
            )
        item.deck_hints["deck_id"] = row.id
        if row.stock_row:
            item.deck_hints["stock_row"] = True
        self.project.lines.append(item)
        self.lines[id(row)] = item

    def end_connection(self, row: _b.EndConnection, index: int) -> None:
        item = EndConnection(f"End {row.end}", end=row.end, direction=row.direction)
        if isinstance(row.stiffness, str):
            item.rotation = "pinned" if row.stiffness.lower() in _PINNED else "rigid"
            item.deck_hints["rotation_token"] = row.stiffness
        else:
            item.rotation = "stiffness"
            item.rotational_stiffness = row.stiffness
        if row.torsion_stiffness is not None:
            if isinstance(row.torsion_stiffness, str):
                lowered = row.torsion_stiffness.lower()
                item.torsion = "free" if lowered in _PINNED else "rigid"
                item.deck_hints["torsion_token"] = row.torsion_stiffness
            else:
                item.torsion = "stiffness"
                item.torsional_stiffness = row.torsion_stiffness
            item.normal = row.normal
            item.pretwist = row.pretwist
        item.deck_hints["row"] = index
        self.lines[id(row.line)].end_connections.append(item)

    def attachment(self, row: _b.Attachment, index: int) -> None:
        density = _option_density(self.model, _SEA_WATER_DENSITY)
        cls = BuoyancyModule if row.volume * density > row.mass else ClumpWeight
        item: Attachment = cls(
            mass=row.mass,
            volume=row.volume,
            drag_area=row.cda,
            added_mass=row.ca,
            axial_drag_area=row.cdax,
        )
        if isinstance(row.arc_length, tuple):
            item.arc_length, item.pitch, item.last_arc_length = row.arc_length
        else:
            item.arc_length = row.arc_length
        line = self.lines[id(row.line)]
        item.name = f"{cls.type_label} {len(line.attachments) + 1}"
        item.deck_hints["row"] = index
        line.attachments.append(item)

    def external_load(self, row: _b.ExternalLoad) -> None:
        item = ExternalLoad(
            f"Load {row.id}",
            body=self.bodies[id(row.body)],
            axes="body" if row.csys.upper() == "L" else "global",
            force=_tuple(row.force),
            linear_damping=_tuple(row.blin),
            quadratic_damping=_tuple(row.bquad),
        )
        item.deck_hints["csys"] = row.csys
        self.project.external_loads.append(item)

    # ------------------------------------------------------------------ options

    def options(self) -> None:
        project = self.project
        rows = list(self.model.options)
        hints: dict[str, Any] = {}
        typed: dict[str, int] = {}
        motion_rows: list[int] = []
        trains: list[int] = []
        for position, row in enumerate(rows):
            identity = _option_key(row.keyword)
            if identity in _MOTION_KEYS:
                if row.values and row.values[0].lower() not in _DISABLED:
                    motion_rows.append(position)
                continue
            if identity == "wavetrain":
                trains.append(position)
            else:
                typed[identity] = position
        if trains:
            typed.pop("waves", None)
        accepted: set[int] = set()
        for identity, position in typed.items():
            row = rows[position]
            values = tuple(row.values)
            try:
                encoded = self._apply(identity, values)
            except (TypeError, ValueError):
                continue
            if encoded is None:
                continue
            accepted.add(position)
            hints[identity] = _hint(position, row, encoded)
        if motion_rows:
            position = motion_rows[-1]
            row = rows[position]
            source = _MOTION_KEYS[_option_key(row.keyword)](file=row.values[0])
            reference = project.motion.reference
            project.motion = source
            source.reference = reference
            accepted.add(position)
            hints["motion"] = _hint(position, row, tuple(row.values))
        if trains:
            sea = MultiTrainSea()
            for position in trains:
                row = rows[position]
                try:
                    train = _decode_train(tuple(row.values))
                except (TypeError, ValueError):
                    continue
                train.name = f"Train {len(sea.trains) + 1}"
                train.deck_hints["option"] = _hint(position, row, _encode_train(train))
                sea.trains.append(train)
                accepted.add(position)
            if sea.trains:
                project.environment.waves = sea
        for position, row in enumerate(rows):
            if position in accepted:
                continue
            extra = ExtraOption(row.keyword, keyword=row.keyword, values=row.values)
            extra.note = row.description
            extra.deck_hints["position"] = position
            project.settings.extra_options.append(extra)
        project.deck_hints["options"] = hints

    def _apply(self, identity: str, values: tuple[str, ...]) -> tuple[str, ...] | None:
        project = self.project
        scalar = _SCALAR_BY_ID.get(identity)
        if scalar is not None:
            value = _decode_scalar(scalar.codec, values)
            owner = _owner(project, scalar.owner)
            owner.set_value(scalar.attr, value)
            return _encode_scalar(scalar.codec, value)
        if identity == "bathymetry":
            old = project.seabed
            seabed = BathymetrySeabed(file=values[0])
            for name in ("stiffness", "damping", "friction", "friction_axial", "friction_lateral"):
                seabed.set_value(name, old.get_value(name))
            project.seabed = seabed
            return values
        if identity == "waves":
            wave = _decode_wave(values, False)
            project.environment.waves = wave
            return _encode_wave(wave)
        if identity == "current":
            current = _decode_current(values)
            project.environment.current = current
            return _encode_current(current)
        if identity == "vesselref":
            reference = _decode_vector(values)
            project.motion.reference = reference
            return _encode_vector(reference)
        if identity == "dynamic_solver":
            if len(values) not in (4, 5):
                raise ValueError("four or five values")
            settings = project.settings
            decoded = [_decode_scalar(codec, (values[i],)) for i, (_, codec) in enumerate(_SOLVER)]
            for (name, _), value in zip(_SOLVER, decoded, strict=True):
                settings.set_value(name, value)
            if len(values) == 5:
                settings.solver_spectral_radius = _fnum(values[4])
            return _encode_solver(settings)
        return None

    def outputs(self) -> None:
        targets: dict[tuple[str, int], ModelObject] = {}
        for line in self.lines.values():
            targets[("line", int(line.deck_hints["deck_id"]))] = line
        for point in self.points.values():
            targets[("point", int(point.deck_hints["deck_id"]))] = point
        for body in self.bodies.values():
            targets[("body", int(body.deck_hints["deck_id"]))] = body
        for rod in self.rods.values():
            targets[("rod", int(rod.deck_hints["deck_id"]))] = rod
        channels = self.project.outputs.channels
        for name in self.model.outputs:
            reference = _b._channel_reference(name)
            target = None if reference is None else targets.get(reference[:2])
            if reference is not None and target is not None:
                match = reference[2]
                channels.append(
                    ObjectChannel(
                        name, target=target, quantity=match.group(1), qualifier=match.group(3)
                    )
                )
            else:
                channels.append(NamedChannel(name, channel=name))


def _option_density(model: DeckModel, fallback: float) -> float:
    value = model.options.value("water_density")
    if value is None:
        return fallback
    try:
        return _fnum(value)
    except ValueError:
        return fallback


def _tuple(value: float | tuple[float, ...]) -> tuple[float, ...]:
    return tuple(value) if isinstance(value, tuple) else (float(value),)


def _multi(value: tuple[float, ...]) -> float | tuple[float, ...]:
    return value[0] if len(value) == 1 else tuple(value)


def _hint(position: int, row: _b.Option, encoded: tuple[str, ...]) -> dict[str, Any]:
    return {
        "position": position,
        "keyword": row.keyword,
        "values": list(row.values),
        "description": row.description,
        "encoded": list(encoded),
    }


def _owner(project: Project, name: str) -> ModelObject:
    owner: ModelObject = project.get_value(name)
    return owner


def _encode_solver(settings: Any) -> tuple[str, ...]:
    tokens = [_encode_scalar(codec, settings.get_value(name))[0] for name, codec in _SOLVER]
    if settings.solver_spectral_radius is not None:
        tokens.append(_num(settings.solver_spectral_radius))
    return tuple(tokens)


# --------------------------------------------------------------------------- writer

_BODY_TOKENS = {
    (Rigid6Body, "cabledyn"): ("Rigid6", frozenset({"rigid6"})),
    (Rigid6Body, "moordyn"): ("Free", frozenset({"free"})),
    (Point3Body, "cabledyn"): ("Point3", frozenset({"point3"})),
    (Point3Body, "moordyn"): ("Point3", frozenset({"point3"})),
}
_ROD_TOKENS = {
    "free": "Free",
    "fixed": "Fixed",
    "pinned": "Pinned",
    "coupled": "Coupled",
    "floater": "Vessel",
    "body": "Body",
    "body_pinned": "BodyPinned",
}
_TABLE_SECTIONS = (
    "LINE TYPES",
    "BODIES",
    "ROD TYPES",
    "RODS",
    "TURBINES",
    "POINTS",
    "LINES",
    "SYROPE IC",
    "SECTIONS",
    "END CONNECTIONS",
    "EQUIVALENT BUOYANCY",
    "ATTACHMENTS",
    "FAILURE",
    "CONTROL",
    "EXTERNAL LOADS",
)


def _hinted(obj: ModelObject, canonical: str, accepted: frozenset[str]) -> str:
    """Return the deck token hint of ``obj`` if it still means the same, else ``canonical``."""
    hint = obj.deck_hints.get("deck_type")
    if isinstance(hint, str) and hint.lower() in accepted:
        return hint
    return canonical


@dataclass
class _Row:
    position: float
    keyword: str
    values: tuple[str, ...]
    description: str | None
    owner: ModelObject
    group: str | None


def _valid_hint(hint: Any) -> bool:
    """Whether an option hint (possibly from a hand-edited file) has the expected shape."""
    return (
        isinstance(hint, dict)
        and isinstance(hint.get("position"), (int, float))
        and not isinstance(hint.get("position"), bool)
        and isinstance(hint.get("keyword"), str)
        and isinstance(hint.get("values"), list)
        and all(isinstance(token, str) for token in hint["values"])
        and isinstance(hint.get("encoded"), list)
        and isinstance(hint.get("description"), (str, type(None)))
    )


# ``wavetrain`` rows add up (every row is a train), so a kept one never overrides.
_SPECIAL_OPTIONS = frozenset({"waves", "current", "vesselref", "bathymetry", "dynamic_solver"})


def typed_option(keyword: str) -> str | None:
    """Return the typed option an ``OPTIONS`` keyword sets, or ``None``.

    ``motionFile``, ``vesselMotion`` and ``vesselRAO`` all set ``"motion"``.
    """
    identity = _option_key(keyword)
    if identity in _MOTION_KEYS:
        return "motion"
    if identity in _SCALAR_BY_ID or identity in _SPECIAL_OPTIONS:  # one row wins
        return identity
    return None


class DeckWriter:
    """Write a :class:`~cabledyn.project.Project` as a CableDyn deck.

    Relative side-file paths resolve from the folder of the project's
    ``deck_path``; with no ``deck_path`` they resolve from the working folder
    (:class:`~cabledyn.project.ProjectStore` anchors them at the project
    file). Points, lines, bodies and rods keep their deck ids where they can;
    external loads are numbered 1, 2, ... in collection order, as the deck
    requires.

    Parameters
    ----------
    project : Project
        The project. Only ``physics`` properties are written.

    Examples
    --------
    >>> from cabledyn.project import DeckReader, DeckWriter
    >>> project = DeckReader.read("examples/spread_3line_chain.dat")  # doctest: +SKIP
    >>> text = DeckWriter(project).to_text()  # doctest: +SKIP
    """

    def __init__(self, project: Project) -> None:
        self.project = project
        self.id_map = IdMap()
        self._sections: dict[str, list[ModelObject]] = {}
        self._deck_bodies: dict[int, Any] = {}
        self._deck_rods: dict[int, _b.Rod] = {}
        self._deck_points: dict[int, _b.Point] = {}
        self._deck_lines: dict[int, _b.Line] = {}

    # ------------------------------------------------------------------ public

    def to_model(self) -> tuple[DeckModel, IdMap]:
        """Build a fresh :class:`~cabledyn.DeckModel` from the project.

        Returns
        -------
        tuple[DeckModel, IdMap]
            The deck model and the id map.

        Raises
        ------
        DeckExportError
            If an object cannot be expressed as a deck row (a missing
            reference, a bad value); ``obj`` names it.
        """
        self.id_map = IdMap()
        self._sections = {name: [] for name in (*_TABLE_SECTIONS, "OPTIONS", "OUTPUTS")}
        project = self.project
        model = DeckModel.new(
            title=project.title,
            path=project.deck_path or "deck.dat",
            caller_driven=project.caller_driven,
        )
        self._ids()
        self._libraries(model)
        self._bodies_and_rods(model)
        self._points(model)
        self._lines(model)
        self._rows(model)
        self._options(model)
        self._outputs(model)
        return model, self.id_map

    def to_text(self, *, validate: bool = True) -> str:
        """Return the canonical deck text of the project.

        Parameters
        ----------
        validate : bool
            Check the text with the native deck rules (:class:`cabledyn.DeckFile`).

        Raises
        ------
        DeckExportError
            If an object cannot be written.
        DeckFormatError
            If ``validate`` is true and the deck violates the native rules.
        """
        model, _ = self.to_model()
        text = self._render(model)
        if validate:
            DeckFile.from_text(text, path=model.path, caller_driven=model.caller_driven)
        return text

    def write(self, path: str | os.PathLike[str], *, overwrite: bool = False) -> Path:
        """Validate the project and write it as a deck file.

        Relative ancillary paths are rebased to the target folder (see
        :meth:`cabledyn.DeckModel.save`).

        Returns
        -------
        pathlib.Path
            The written file.
        """
        model, _ = self.to_model()
        self._render(model)
        return model.save(path, overwrite=overwrite)

    def validate(self) -> list[Issue]:
        """Write the deck and apply the native deck rules (validation layer 3).

        Returns
        -------
        list[Issue]
            Layer-3 issues, mapped to objects where the failing row is known.
        """
        try:
            model, _ = self.to_model()
            text = self._render(model)
        except DeckExportError as exc:
            return [Issue(Severity.ERROR, str(exc), exc.obj or self.project, None, 3)]
        label = "<project deck>"
        try:
            DeckFile.from_text(
                text, path=model.path, caller_driven=model.caller_driven, label=label
            )
        except DeckFormatError as exc:
            message = str(exc)
            obj: ModelObject = self.project
            match = re.match(rf"^{re.escape(label)}:(\d+): ?(.*)$", message, re.DOTALL)
            if match is not None:
                obj = self.id_map.object_at_row(int(match.group(1)) - 1) or self.project
                message = match.group(2)
            elif message.startswith(f"{label}: "):
                message = message[len(label) + 2 :]
            return [Issue(Severity.ERROR, message, obj, None, 3)]
        return []

    # ------------------------------------------------------------------ helpers

    def _render(self, model: DeckModel) -> str:
        try:
            text = model.to_text(validate=False)
        except (DeckReferenceError, TypeError, ValueError) as exc:
            raise DeckExportError(str(exc)) from exc
        self._map_rows(text)
        return text

    def _map_rows(self, text: str) -> None:
        rows = self.id_map.rows
        rows.clear()
        current: list[ModelObject] | None = None
        skip = 0
        count = 0
        for index, line in enumerate(text.splitlines()):
            if line.startswith("---------------------"):
                name = line.strip("- ").strip()
                current = self._sections.get(name)
                skip = 2 if name in _TABLE_SECTIONS else 0
                count = 0
                continue
            if current is None:
                continue
            if skip:
                skip -= 1
                continue
            if count < len(current):
                rows[index] = current[count]
            count += 1

    def _fail(self, obj: ModelObject, message: str) -> DeckExportError:
        return DeckExportError(f"{obj.label()}: {message}", obj)

    def _call(self, obj: ModelObject, action: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
        try:
            return action(*args, **kwargs)
        except DeckExportError:
            raise
        except (DeckReferenceError, KeyError, TypeError, ValueError) as exc:
            raise self._fail(obj, str(exc)) from exc

    def _ids(self) -> None:
        project = self.project
        for kind, items in (
            ("body", list(project.bodies)),
            ("rod", list(project.rods)),
            ("point", list(project.points)),
            ("line", list(project.lines)),
        ):
            used: set[int] = set()
            assigned: dict[int, int] = {}
            for position, item in enumerate(items):
                hint = item.deck_hints.get("deck_id")
                if isinstance(hint, int) and not isinstance(hint, bool) and hint >= 1:
                    if hint in used:
                        continue
                    used.add(hint)
                    assigned[position] = hint
            following = max(used, default=0)
            for position, item in enumerate(items):
                if position not in assigned:
                    following += 1
                    assigned[position] = following
                self.id_map.add(item, kind, assigned[position])

    def _id(self, obj: ModelObject) -> int:
        found = self.id_map.deck_id(obj)
        assert isinstance(found, int)  # every numbered object gets an id in _ids
        return found

    def _libraries(self, model: DeckModel) -> None:
        for item in self.project.line_types:
            self._call(item, self._line_type, model, item)
            self.id_map.add(item, "line_type", item.name)
            self._sections["LINE TYPES"].append(item)
        for rod_type in self.project.rod_types:
            self._call(
                rod_type,
                model.add_rod_type,
                rod_type.name,
                diam=rod_type.diameter,
                mass=rod_type.mass_per_length,
                cd=rod_type.drag,
                ca=rod_type.added_mass,
                cd_end=rod_type.end_drag,
                ca_end=rod_type.end_added_mass,
                cd_ax=rod_type.axial_drag,
                ca_ax=rod_type.axial_added_mass,
            )
            self.id_map.add(rod_type, "rod_type", rod_type.name)
            self._sections["ROD TYPES"].append(rod_type)

    def _line_type(self, model: DeckModel, item: LineType) -> None:
        axial = item.axial
        ea: float | tuple[float, ...] | _b.SyropeEA
        if isinstance(axial, LinearAxial):
            ea = axial.stiffness
        elif isinstance(axial, ViscoelasticAxial):
            if axial.dynamic_stiffness is not None:
                ea = (axial.static_stiffness, axial.dynamic_stiffness)
            elif axial.alpha_mbl is not None and axial.beta is not None:
                ea = (axial.static_stiffness, axial.alpha_mbl, axial.beta)
            else:
                raise self._fail(item, "the viscoelastic model needs its dynamic stiffness")
        elif isinstance(axial, SyropeAxial):
            if not axial.settings_file:
                raise self._fail(item, "the Syrope model needs a settings file")
            ea = _b.SyropeEA(axial.settings_file, axial.alpha, axial.beta)
        else:
            raise self._fail(item, f"{axial.label()} cannot be written")
        bending = item.bending
        model.add_line_type(
            item.name,
            diam=item.diameter,
            mass=item.mass_per_length,
            ea=ea,
            ba=axial.deck_damping(),
            ei=bending.bending_stiffness,
            cdn=item.drag_normal,
            cdt=item.drag_axial,
            can=item.added_mass_normal,
            cat=item.added_mass_axial,
            gas=bending.shear_stiffness,
            gj=bending.torsional_stiffness,
            irt=bending.rotary_inertia_axial,
            irn=bending.rotary_inertia_normal,
        )

    def _bodies_and_rods(self, model: DeckModel) -> None:
        self._deck_bodies = {}
        for body in self.project.bodies:
            self._deck_bodies[id(body)] = self._call(body, self._body, model, body)
            self._sections["BODIES"].append(body)
        self._deck_rods = {}
        for rod in self.project.rods:
            self._deck_rods[id(rod)] = self._call(rod, self._rod, model, rod)
            self._sections["RODS"].append(rod)
        for turbine in self.project.turbines:
            ptfm = turbine.platform_displacement or None
            self._call(turbine, model.add_turbine, turbine.number, *turbine.position, ptfm=ptfm)
            self.id_map.add(turbine, "turbine", turbine.number)
            self._sections["TURBINES"].append(turbine)

    def _body(self, model: DeckModel, body: Body) -> Any:
        if isinstance(body, Floater):
            if body.kind == "coupled":
                token = _hinted(body, "Coupled", frozenset({"coupled", "cpld"}))
            else:
                token = _hinted(body, "Vessel", frozenset({"vessel", "ves"}))
        else:
            canonical, accepted = _BODY_TOKENS[(type(body), body.row_format)]
            token = _hinted(body, canonical, accepted)
        x, y, z = body.position
        roll, pitch, yaw = body.orientation
        if body.row_format == "moordyn":
            return model.add_moordyn_body(
                self._id(body),
                token,
                x,
                y,
                z,
                roll=roll,
                pitch=pitch,
                yaw=yaw,
                mass=body.mass,
                cg=_multi(body.centre_of_gravity),
                inertia=_multi(body.inertia) if body.inertia else 0.0,
                volume=body.volume,
                cda=_multi(body.drag_area),
                ca=_multi(body.added_mass),
            )
        if len(body.drag_area) != 1 or len(body.added_mass) != 1:
            raise self._fail(body, "the CableDyn layout takes one drag area and one added mass")
        inertia: tuple[float, float, float] | None = None
        if body.inertia:
            if len(body.inertia) != 3:
                raise self._fail(body, "the CableDyn layout takes three inertias")
            inertia = (body.inertia[0], body.inertia[1], body.inertia[2])
        return model.add_body(
            self._id(body),
            token,
            x,
            y,
            z,
            roll=roll,
            pitch=pitch,
            yaw=yaw,
            mass=body.mass,
            volume=body.volume,
            c33=body.heave_stiffness,
            c44=body.roll_stiffness,
            c55=body.pitch_stiffness,
            cda=body.drag_area[0],
            ca=body.added_mass[0],
            inertia=inertia,
        )

    def _rod(self, model: DeckModel, rod: Rod) -> _b.Rod:
        if rod.rod_type is None:
            raise self._fail(rod, "has no rod type")
        if rod.rod_type.parent is not self.project:
            raise self._fail(rod, "uses a rod type outside the project")
        canonical = _ROD_TOKENS[rod.attachment]
        accepted = frozenset(
            token
            for token in {canonical.lower(), *_ROD_ATTACHMENT_ALIASES, "bodypin"}
            if _rod_kind(token) == rod.attachment
        )
        token = _hinted(rod, canonical, accepted)
        body = None
        if rod.body is not None:
            body = self._deck_bodies.get(id(rod.body))
            if body is None:
                raise self._fail(rod, "its body is not part of the project")
        rod_type = model.rod_types[rod.rod_type.name]
        return model.add_rod(
            self._id(rod),
            rod_type,
            token,
            rod.position_a,
            rod.position_b,
            rod.segments,
            outputs=rod.outputs,
            body=body,
        )

    def _rod_end(self, model: DeckModel, endpoint: RodEndpoint) -> _b.RodEnd:
        rod = endpoint.parent
        if not isinstance(rod, Rod) or id(rod) not in self._deck_rods:
            raise DeckExportError(f"{endpoint.label()} is not part of the project", endpoint)
        return model.rod_end(self._deck_rods[id(rod)], endpoint.end)

    def _points(self, model: DeckModel) -> None:
        self._deck_points = {}
        for point in self.project.points:
            self._deck_points[id(point)] = self._call(point, self._point, model, point)
            self._sections["POINTS"].append(point)

    def _point(self, model: DeckModel, point: Point) -> _b.Point:
        body = None
        rod_end = None
        turbine = None
        if isinstance(point, FixedPoint):
            token = _hinted(point, "Fixed", frozenset({"fixed"}))
        elif isinstance(point, FloaterPoint):
            if point.kind == "coupled":
                token = _hinted(point, "Coupled", frozenset({"coupled"}))
            else:
                token = _hinted(point, "Vessel", frozenset({"vessel"}))
        elif isinstance(point, FreePoint):
            token = _hinted(point, "Free", frozenset({"free"}))
        elif isinstance(point, ConnectPoint):
            token = _hinted(point, "Connect", frozenset({"connect"}))
        elif isinstance(point, BodyPoint):
            token = _hinted(point, "Body", frozenset({"body"}))
            if point.body is None or id(point.body) not in self._deck_bodies:
                raise self._fail(point, "is not attached to a body of the project")
            body = self._deck_bodies[id(point.body)]
        elif isinstance(point, RodPoint):
            token = _hinted(point, "Rod", frozenset({"rod"}))
            if point.rod_end is None:
                raise self._fail(point, "is not attached to a rod end")
            rod_end = self._rod_end(model, point.rod_end)
        elif isinstance(point, TurbinePoint):
            token = _hinted(point, "Turbine", frozenset({"turbine", "t"}))
            turbine = point.number()
            if turbine is None:
                raise self._fail(point, "has no turbine")
            if point.turbine is not None and point.turbine.parent is not self.project:
                raise self._fail(point, "refers to a turbine outside the project")
        else:
            raise self._fail(point, "has no deck point type")
        x, y, z = point.position
        return model.add_point(
            self._id(point),
            token,
            x,
            y,
            z,
            mass=point.mass,
            volume=point.volume,
            cda=point.drag_area,
            ca=point.added_mass,
            body=body,
            rod_end=rod_end,
            turbine=turbine,
        )

    def _end(self, model: DeckModel, line: Line, target: LineEndTarget | None) -> Any:
        if target is None:
            raise self._fail(line, "has an unattached end")
        if isinstance(target, RodEndpoint):
            return self._rod_end(model, target)
        if id(target) not in self._deck_points:
            raise self._fail(line, f"is attached to {target.label()}, outside the project")
        return self._deck_points[id(target)]

    def _lines(self, model: DeckModel) -> None:
        self._deck_lines = {}
        for line in self.project.lines:
            deck_line = self._call(line, self._line, model, line)
            self._deck_lines[id(line)] = deck_line
            self._sections["LINES"].append(line)
            if not (deck_line.stock_row and len(deck_line.sections) == 1):
                self._sections["SECTIONS"].extend(line.sections)

    def _line(self, model: DeckModel, line: Line) -> _b.Line:
        deck_line = model.add_line(
            self._id(line),
            self._end(model, line, line.end_a),
            self._end(model, line, line.end_b),
            outputs=line.outputs,
            stock_row=bool(line.deck_hints.get("stock_row", False)),
        )
        for section in line.sections:
            if section.line_type is None:
                raise self._fail(section, "has no line type")
            if section.line_type.parent is not self.project:
                raise self._fail(section, "uses a line type outside the project")
            model.add_section(
                deck_line,
                model.line_types[section.line_type.name],
                section.length,
                section.segments,
            )
        return deck_line

    def _line_of(self, obj: ModelObject, line: Line | None) -> _b.Line:
        if line is None or id(line) not in self._deck_lines:
            raise self._fail(obj, "refers to a line outside the project")
        return self._deck_lines[id(line)]

    def _lines_of(self, obj: ModelObject, lines: Iterable[Line]) -> list[_b.Line]:
        return [self._line_of(obj, line) for line in lines]

    def _rows(self, model: DeckModel) -> None:
        project = self.project
        for row in project.syrope_history:
            lines = self._lines_of(row, row.lines)
            self._call(row, model.add_syrope_ic, lines, row.max_tension, row.mean_tension)
            self._sections["SYROPE IC"].append(row)
        connections = _ordered(
            (row, line) for line in project.lines for row in line.end_connections
        )
        for connection, line in connections:
            self._call(connection, self._end_connection, model, connection, line)
            self._sections["END CONNECTIONS"].append(connection)
        for eq in project.equivalent_buoyancy:
            line_type = eq.line_type
            if line_type is None or line_type.parent is not project:
                raise self._fail(eq, "needs a line type of the project")
            deck_type = model.line_types[line_type.name]
            self._call(
                eq, model.add_equivalent_buoyancy, deck_type, eq.diameter, eq.submerged_weight
            )
            self._sections["EQUIVALENT BUOYANCY"].append(eq)
        attachments = _ordered((row, line) for line in project.lines for row in line.attachments)
        for attachment, line in attachments:
            self._call(attachment, self._attachment, model, attachment, line)
            self._sections["ATTACHMENTS"].append(attachment)
        for failure in project.failures:
            point = failure.point
            if point is None or id(point) not in self._deck_points:
                raise self._fail(failure, "needs a point of the project")
            self._call(
                failure,
                model.add_failure,
                self._deck_points[id(point)],
                self._lines_of(failure, failure.lines),
                fail_time=failure.time,
                fail_tension=failure.tension,
            )
            self._sections["FAILURE"].append(failure)
        for control in project.controls:
            lines = self._lines_of(control, control.lines)
            self._call(control, model.add_control, control.channel, lines)
            self._sections["CONTROL"].append(control)
        for number, load in enumerate(project.external_loads, start=1):
            if load.body is None or id(load.body) not in self._deck_bodies:
                raise self._fail(load, "needs a body of the project")
            hint = load.deck_hints.get("csys")
            letter = "L" if load.axes == "body" else "G"
            csys = hint if isinstance(hint, str) and hint.upper() == letter else letter
            self._call(
                load,
                model.add_external_load,
                number,
                self._deck_bodies[id(load.body)],
                csys=csys,
                force=_multi(load.force),
                blin=_multi(load.linear_damping),
                bquad=_multi(load.quadratic_damping),
            )
            self.id_map.add(load, "external_load", number)
            self._sections["EXTERNAL LOADS"].append(load)

    def _end_connection(self, model: DeckModel, row: EndConnection, line: Line) -> None:
        stiffness: float | str
        if row.rotation == "stiffness":
            stiffness = row.rotational_stiffness
        else:
            accepted = _PINNED if row.rotation == "pinned" else _RIGID
            hint = row.deck_hints.get("rotation_token")
            canonical = "Pinned" if row.rotation == "pinned" else "Rigid"
            stiffness = hint if isinstance(hint, str) and hint.lower() in accepted else canonical
        torsion: float | str | None = None
        if row.torsion == "stiffness":
            torsion = row.torsional_stiffness
        elif row.torsion in {"free", "rigid"}:
            accepted = _PINNED if row.torsion == "free" else _RIGID
            hint = row.deck_hints.get("torsion_token")
            canonical = "Free" if row.torsion == "free" else "Rigid"
            torsion = hint if isinstance(hint, str) and hint.lower() in accepted else canonical
        model.add_end_connection(
            self._line_of(row, line),
            row.end,
            stiffness,
            row.direction,
            torsion_stiffness=torsion,
            normal=row.normal if torsion is not None else None,
            pretwist=row.pretwist if torsion is not None else None,
        )

    def _attachment(self, model: DeckModel, row: Attachment, line: Line) -> None:
        arc: float | tuple[float, float, float] = row.arc_length
        if row.pitch is not None or row.last_arc_length is not None:
            if row.pitch is None or row.last_arc_length is None:
                raise self._fail(row, "a series needs both pitch and last_arc_length")
            arc = (row.arc_length, row.pitch, row.last_arc_length)
        model.add_attachment(
            self._line_of(row, line),
            arc,
            mass=row.mass,
            volume=row.volume,
            cda=row.drag_area,
            ca=row.added_mass,
            cdax=row.axial_drag_area,
        )

    # ------------------------------------------------------------------ options

    def _options(self, model: DeckModel) -> None:
        project = self.project
        stored = project.deck_hints.get("options")
        hints: dict[str, Any] = stored if isinstance(stored, dict) else {}
        rows: list[_Row] = []
        fresh = 1.0e9

        def hinted(
            hint: Any, keyword: str, encoded: tuple[str, ...], owner: ModelObject, group: str
        ) -> None:
            nonlocal fresh
            if _valid_hint(hint):
                same = list(encoded) == hint["encoded"]
                rows.append(
                    _Row(
                        float(hint["position"]),
                        hint["keyword"] or keyword,
                        tuple(hint["values"]) if same else encoded,
                        hint["description"],
                        owner,
                        group,
                    )
                )
            else:
                fresh += 1.0
                rows.append(_Row(fresh, keyword, encoded, None, owner, group))

        def add(identity: str, keyword: str, encoded: tuple[str, ...], owner: ModelObject) -> None:
            hinted(hints.get(identity), keyword, encoded, owner, identity.split(":")[0])

        for scalar in _SCALARS:
            owner = _owner(project, scalar.owner)
            value = owner.get_value(scalar.attr)
            if value is not None:
                add(scalar.identity, scalar.keyword, _encode_scalar(scalar.codec, value), owner)
        seabed = project.seabed
        if isinstance(seabed, BathymetrySeabed) and seabed.file:
            add("bathymetry", "bathymetryFile", (seabed.file,), seabed)
        settings = project.settings
        if settings.solver_row_given():
            add("dynamic_solver", "dynamic_solver", _encode_solver(settings), settings)
        environment = project.environment
        waves = environment.waves
        if isinstance(waves, MultiTrainSea):
            for train in waves.trains:
                encoded = _encode_train(train)
                hinted(train.deck_hints.get("option"), "wavetrain", encoded, train, "wavetrain")
        elif not isinstance(waves, NoWaves) or "waves" in hints:
            add("waves", "waves", _encode_wave(waves), environment)
        current = environment.current
        if not isinstance(current, NoCurrent) or "current" in hints:
            add("current", "current", _encode_current(current), environment)
        motion = project.motion
        self._motion(motion, hints, add)
        if motion.reference is not None:
            add("vesselref", "vesselRef", _encode_vector(motion.reference), motion)
        # Kept rows go where they were; a kept row of a typed option is moved
        # before every typed row of that option, so the typed property wins.
        last_extra: dict[str, float] = {}
        for extra in settings.extra_options:
            position = extra.deck_hints.get("position")
            if isinstance(position, int) and not isinstance(position, bool):
                order = float(position)
            else:
                fresh += 1.0
                order = fresh
            group = typed_option(extra.keyword)
            rows.append(_Row(order, extra.keyword, tuple(extra.values), extra.note, extra, None))
            if group is not None:
                last_extra[group] = max(order, last_extra.get(group, order))
        for row in rows:
            if row.group is not None and row.group in last_extra:
                row.position = max(row.position, last_extra[row.group] + 0.5)
        rows.sort(key=lambda row: row.position)
        for row in rows:
            self._call(
                row.owner, model.options.add, row.keyword, *row.values, description=row.description
            )
            self._sections["OPTIONS"].append(row.owner)

    def _motion(
        self,
        motion: MotionSource,
        hints: dict[str, Any],
        add: Callable[[str, str, tuple[str, ...], ModelObject], None],
    ) -> None:
        if isinstance(motion, NoMotion):
            return
        file = motion.get_value("file")
        if not file:
            raise self._fail(motion, "needs a file")
        keyword = "vesselRAO" if isinstance(motion, FloaterRAO) else _MOTION_KEYWORDS[type(motion)]
        hint = hints.get("motion")
        same_keyword = (
            isinstance(hint, dict)
            and _valid_hint(hint)
            and _option_key(str(hint["keyword"])) == keyword.lower()
        )
        add("motion" if same_keyword else "motion:changed", keyword, (file,), motion)

    def _outputs(self, model: DeckModel) -> None:
        for channel in self.project.outputs.channels:
            if isinstance(channel, ObjectChannel):
                target = channel.target
                deck_id = None if target is None else self.id_map.deck_id(target)
                if not isinstance(deck_id, int):
                    raise self._fail(channel, "refers to an object without a deck id")
                name = channel.channel_name(deck_id)
            else:
                name = channel.channel_name(None)
            self._call(channel, model.outputs.add, name)
            self._sections["OUTPUTS"].append(channel)


def _rod_kind(token: str) -> str:
    kind = _ROD_ATTACHMENT_ALIASES.get(token, token)
    return {
        "vessel": "floater",
        "bodypinned": "body_pinned",
        "bodypin": "body_pinned",
    }.get(kind, kind)


def _ordered(pairs: Iterable[tuple[ModelObject, Line]]) -> list[tuple[Any, Line]]:
    items = list(pairs)
    keyed = []
    for sequence, (obj, line) in enumerate(items):
        row = obj.deck_hints.get("row")
        keyed.append(
            ((float(row) if isinstance(row, int) else 1.0e9 + sequence), sequence, obj, line)
        )
    keyed.sort(key=lambda item: (item[0], item[1]))
    return [(obj, line) for _, _, obj, line in keyed]
