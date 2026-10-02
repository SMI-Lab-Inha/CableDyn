# SPDX-License-Identifier: Apache-2.0
"""Physical dimensions, units, and display-unit preferences.

Every quantity in :mod:`cabledyn.project` is stored in the unit the deck uses:
SI units, except angles, which are stored in degrees (the deck convention). A
:class:`Dimension` names what a quantity measures and its storage unit, a
:class:`Unit` converts between a display unit and the storage unit, and a
:class:`UnitSystem` holds the display unit chosen for each dimension.
Conversions are linear (``stored = value * scale``); no supported unit has an
offset.

Examples
--------
>>> from cabledyn.project.units import FORCE, UnitSystem, convert
>>> convert(2.5, "kN", "N")
2500.0
>>> UnitSystem.engineering().format(1.5e6, FORCE)
'1500 kN'
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

__all__ = [
    "ACCELERATION",
    "ANGLE",
    "AREA",
    "AXIAL_DAMPING",
    "BENDING_STIFFNESS",
    "DENSITY",
    "DIMENSIONLESS",
    "DIMENSIONS",
    "FORCE",
    "FORCE_PER_LENGTH",
    "FREQUENCY",
    "INERTIA",
    "LENGTH",
    "LINEAR_DAMPING",
    "MASS",
    "MASS_PER_LENGTH",
    "PRESSURE_PER_LENGTH",
    "PRESSURE_TIME_PER_LENGTH",
    "QUADRATIC_DAMPING",
    "ROTARY_INERTIA_PER_LENGTH",
    "ROTATIONAL_STIFFNESS",
    "STIFFNESS",
    "TIME",
    "UNITS",
    "VELOCITY",
    "VOLUME",
    "Dimension",
    "Unit",
    "UnitError",
    "UnitSystem",
    "convert",
    "parse_quantity",
    "unit",
    "units_for",
]


class UnitError(ValueError):
    """A unit is unknown or does not measure the expected dimension."""


@dataclass(frozen=True)
class Dimension:
    """A physical dimension and its SI unit symbol.

    Attributes
    ----------
    key : str
        Stable identifier, used in project files and preferences.
    label : str
        Human-readable name.
    si : str
        Symbol of the unit in which values are stored (SI; degrees for angles).
    """

    key: str
    label: str
    si: str


DIMENSIONLESS = Dimension("dimensionless", "Dimensionless", "-")
LENGTH = Dimension("length", "Length", "m")
AREA = Dimension("area", "Area", "m^2")
VOLUME = Dimension("volume", "Volume", "m^3")
MASS = Dimension("mass", "Mass", "kg")
TIME = Dimension("time", "Time", "s")
FREQUENCY = Dimension("frequency", "Frequency", "Hz")
ANGLE = Dimension("angle", "Angle", "deg")
VELOCITY = Dimension("velocity", "Velocity", "m/s")
ACCELERATION = Dimension("acceleration", "Acceleration", "m/s^2")
DENSITY = Dimension("density", "Density", "kg/m^3")
MASS_PER_LENGTH = Dimension("mass_per_length", "Mass per length", "kg/m")
FORCE = Dimension("force", "Force", "N")
FORCE_PER_LENGTH = Dimension("force_per_length", "Force per length", "N/m")
STIFFNESS = Dimension("stiffness", "Stiffness", "N/m")
ROTATIONAL_STIFFNESS = Dimension("rotational_stiffness", "Rotational stiffness", "N m/rad")
BENDING_STIFFNESS = Dimension("bending_stiffness", "Bending stiffness", "N m^2")
AXIAL_DAMPING = Dimension("axial_damping", "Axial damping", "N s")
LINEAR_DAMPING = Dimension("linear_damping", "Linear damping", "N s/m")
QUADRATIC_DAMPING = Dimension("quadratic_damping", "Quadratic damping", "N s^2/m^2")
INERTIA = Dimension("inertia", "Moment of inertia", "kg m^2")
ROTARY_INERTIA_PER_LENGTH = Dimension(
    "rotary_inertia_per_length", "Rotary inertia per length", "kg m"
)
PRESSURE_PER_LENGTH = Dimension("pressure_per_length", "Pressure per length", "Pa/m")
PRESSURE_TIME_PER_LENGTH = Dimension(
    "pressure_time_per_length", "Pressure-time per length", "Pa s/m"
)

DIMENSIONS: dict[str, Dimension] = {
    dim.key: dim
    for dim in (
        DIMENSIONLESS,
        LENGTH,
        AREA,
        VOLUME,
        MASS,
        TIME,
        FREQUENCY,
        ANGLE,
        VELOCITY,
        ACCELERATION,
        DENSITY,
        MASS_PER_LENGTH,
        FORCE,
        FORCE_PER_LENGTH,
        STIFFNESS,
        ROTATIONAL_STIFFNESS,
        BENDING_STIFFNESS,
        AXIAL_DAMPING,
        LINEAR_DAMPING,
        QUADRATIC_DAMPING,
        INERTIA,
        ROTARY_INERTIA_PER_LENGTH,
        PRESSURE_PER_LENGTH,
        PRESSURE_TIME_PER_LENGTH,
    )
}
"""Every known dimension, by key."""


@dataclass(frozen=True)
class Unit:
    """A display unit: ``si_value = value * scale``.

    Attributes
    ----------
    symbol : str
        Unit symbol. A symbol may measure several dimensions (``N/m`` is a
        force per length or a stiffness); then
        :func:`~cabledyn.project.units.unit` needs the dimension.
    dimension : Dimension
        What the unit measures.
    scale : float
        SI value of one unit.
    """

    symbol: str
    dimension: Dimension
    scale: float

    def to_si(self, value: float) -> float:
        """Convert a value in this unit to SI."""
        return value * self.scale

    def from_si(self, value: float) -> float:
        """Convert an SI value to this unit."""
        return value / self.scale


def _units() -> tuple[Unit, ...]:
    table: list[tuple[str, Dimension, float]] = [
        ("-", DIMENSIONLESS, 1.0),
        ("%", DIMENSIONLESS, 0.01),
        ("m", LENGTH, 1.0),
        ("mm", LENGTH, 1.0e-3),
        ("cm", LENGTH, 1.0e-2),
        ("km", LENGTH, 1.0e3),
        ("in", LENGTH, 0.0254),
        ("ft", LENGTH, 0.3048),
        ("m^2", AREA, 1.0),
        ("mm^2", AREA, 1.0e-6),
        ("cm^2", AREA, 1.0e-4),
        ("m^3", VOLUME, 1.0),
        ("l", VOLUME, 1.0e-3),
        ("kg", MASS, 1.0),
        ("t", MASS, 1.0e3),
        ("s", TIME, 1.0),
        ("min", TIME, 60.0),
        ("h", TIME, 3600.0),
        ("Hz", FREQUENCY, 1.0),
        ("rad/s", FREQUENCY, 1.0 / (2.0 * math.pi)),
        ("deg", ANGLE, 1.0),
        ("rad", ANGLE, 180.0 / math.pi),
        ("m/s", VELOCITY, 1.0),
        ("kn", VELOCITY, 1852.0 / 3600.0),
        ("m/s^2", ACCELERATION, 1.0),
        ("g0", ACCELERATION, 9.80665),
        ("kg/m^3", DENSITY, 1.0),
        ("t/m^3", DENSITY, 1.0e3),
        ("kg/m", MASS_PER_LENGTH, 1.0),
        ("t/m", MASS_PER_LENGTH, 1.0e3),
        ("N", FORCE, 1.0),
        ("kN", FORCE, 1.0e3),
        ("MN", FORCE, 1.0e6),
        ("GN", FORCE, 1.0e9),
        ("tf", FORCE, 9.80665e3),
        ("N/m", FORCE_PER_LENGTH, 1.0),
        ("kN/m", FORCE_PER_LENGTH, 1.0e3),
        ("N/m", STIFFNESS, 1.0),
        ("kN/m", STIFFNESS, 1.0e3),
        ("MN/m", STIFFNESS, 1.0e6),
        ("N m/rad", ROTATIONAL_STIFFNESS, 1.0),
        ("kN m/rad", ROTATIONAL_STIFFNESS, 1.0e3),
        ("MN m/rad", ROTATIONAL_STIFFNESS, 1.0e6),
        ("N m/deg", ROTATIONAL_STIFFNESS, 180.0 / math.pi),
        ("N m^2", BENDING_STIFFNESS, 1.0),
        ("kN m^2", BENDING_STIFFNESS, 1.0e3),
        ("MN m^2", BENDING_STIFFNESS, 1.0e6),
        ("N s", AXIAL_DAMPING, 1.0),
        ("kN s", AXIAL_DAMPING, 1.0e3),
        ("N s/m", LINEAR_DAMPING, 1.0),
        ("kN s/m", LINEAR_DAMPING, 1.0e3),
        ("N s^2/m^2", QUADRATIC_DAMPING, 1.0),
        ("kN s^2/m^2", QUADRATIC_DAMPING, 1.0e3),
        ("kg m^2", INERTIA, 1.0),
        ("t m^2", INERTIA, 1.0e3),
        ("kg m", ROTARY_INERTIA_PER_LENGTH, 1.0),
        ("Pa/m", PRESSURE_PER_LENGTH, 1.0),
        ("kPa/m", PRESSURE_PER_LENGTH, 1.0e3),
        ("Pa s/m", PRESSURE_TIME_PER_LENGTH, 1.0),
        ("kPa s/m", PRESSURE_TIME_PER_LENGTH, 1.0e3),
    ]
    return tuple(Unit(symbol, dim, scale) for symbol, dim, scale in table)


UNITS: tuple[Unit, ...] = _units()
"""Every known unit. A symbol may measure more than one dimension (``N/m``)."""


def unit(symbol: str, dimension: Dimension | None = None) -> Unit:
    """Return the unit with ``symbol``.

    Parameters
    ----------
    symbol : str
        Unit symbol, for example ``"kN"``.
    dimension : Dimension | None
        The dimension it must measure; needed only when the symbol is
        shared by several dimensions (``N/m``: force per length or stiffness).

    Returns
    -------
    Unit
        The unit.

    Raises
    ------
    UnitError
        If the symbol is unknown, does not measure ``dimension``, or is
        ambiguous without one.
    """
    found = [item for item in UNITS if item.symbol == symbol]
    if dimension is not None:
        found = [item for item in found if item.dimension == dimension]
    if not found:
        where = "" if dimension is None else f" for {dimension.label.lower()}"
        raise UnitError(f"unknown unit {symbol!r}{where}")
    if len(found) > 1:
        raise UnitError(f"unit {symbol!r} is ambiguous; give its dimension")
    return found[0]


def units_for(dimension: Dimension) -> list[Unit]:
    """Return the units that measure ``dimension``, SI unit first."""
    found = [item for item in UNITS if item.dimension == dimension]
    found.sort(key=lambda item: item.scale != 1.0)
    return found


def convert(value: float, source: str, target: str, dimension: Dimension | None = None) -> float:
    """Convert ``value`` from unit ``source`` to unit ``target``.

    Parameters
    ----------
    value : float
        Value in ``source`` units.
    source, target : str
        Unit symbols of the same dimension.
    dimension : Dimension | None
        The dimension, needed only for shared symbols such as ``N/m``.

    Returns
    -------
    float
        Value in ``target`` units.

    Raises
    ------
    UnitError
        If a unit is unknown or the dimensions differ.
    """
    a, b = unit(source, dimension), unit(target, dimension)
    if a.dimension != b.dimension:
        raise UnitError(f"cannot convert {source} ({a.dimension.label}) to {target}")
    if a is b:
        return float(value)
    return b.from_si(a.to_si(float(value)))


_QUANTITY = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(.*?)\s*$")


def parse_quantity(text: str, dimension: Dimension, default_unit: str | None = None) -> float:
    """Parse text such as ``"850 m"`` or ``"1.2 MN"`` into an SI value.

    Parameters
    ----------
    text : str
        A number optionally followed by a unit symbol.
    dimension : Dimension
        The expected dimension.
    default_unit : str | None
        Unit of a bare number; SI when ``None``.

    Returns
    -------
    float
        The SI value.

    Raises
    ------
    UnitError
        If the text is not a number with a unit of ``dimension``.
    """
    match = _QUANTITY.match(text)
    if match is None:
        raise UnitError(f"{text!r} is not a number with an optional unit")
    value = float(match.group(1))
    symbol = match.group(2) or default_unit or dimension.si
    return unit(symbol, dimension).to_si(value)


class UnitSystem:
    """Display-unit preferences: one unit per dimension (SI where unset).

    Parameters
    ----------
    choices : Mapping[str, str] | None
        Dimension key to unit symbol.

    Raises
    ------
    UnitError
        If a dimension or unit is unknown, or a unit measures another
        dimension.
    """

    def __init__(self, choices: Mapping[str, str] | None = None) -> None:
        self._choices: dict[str, str] = {}
        for key, symbol in (choices or {}).items():
            self.set(key, symbol)

    @classmethod
    def si(cls) -> UnitSystem:
        """Return SI display units for every dimension."""
        return cls()

    @classmethod
    def engineering(cls) -> UnitSystem:
        """Return offshore-engineering display units (kN, t, kN m^2, ...)."""
        return cls(
            {
                "force": "kN",
                "mass": "t",
                "force_per_length": "kN/m",
                "stiffness": "kN/m",
                "bending_stiffness": "kN m^2",
                "rotational_stiffness": "kN m/rad",
                "axial_damping": "kN s",
                "inertia": "t m^2",
            }
        )

    def set(self, dimension: str, symbol: str) -> None:
        """Choose the display unit of a dimension.

        Parameters
        ----------
        dimension : str
            Dimension key, for example ``"force"``.
        symbol : str
            Unit symbol of that dimension.
        """
        if dimension not in DIMENSIONS:
            raise UnitError(f"unknown dimension {dimension!r}")
        unit(symbol, DIMENSIONS[dimension])
        self._choices[dimension] = symbol

    def unit_for(self, dimension: Dimension) -> Unit:
        """Return the display unit chosen for ``dimension``."""
        return unit(self._choices.get(dimension.key, dimension.si), dimension)

    def to_display(self, value: float, dimension: Dimension) -> float:
        """Convert an SI value to the display unit of ``dimension``."""
        return self.unit_for(dimension).from_si(value)

    def from_display(self, value: float, dimension: Dimension) -> float:
        """Convert a value in the display unit of ``dimension`` to SI."""
        return self.unit_for(dimension).to_si(value)

    def format(self, value: float, dimension: Dimension, digits: int = 6) -> str:
        """Format an SI value in its display unit, for example ``'850 m'``.

        Parameters
        ----------
        value : float
            SI value.
        dimension : Dimension
            Its dimension.
        digits : int
            Significant digits.

        Returns
        -------
        str
            Number and unit symbol (no symbol for dimensionless values).
        """
        shown = self.unit_for(dimension)
        number = f"{shown.from_si(value):.{digits}g}"
        return number if shown.symbol == "-" else f"{number} {shown.symbol}"

    def to_dict(self) -> dict[str, str]:
        """Return the non-SI choices, for a project file."""
        return dict(self._choices)

    @classmethod
    def from_dict(cls, data: Mapping[str, str]) -> UnitSystem:
        """Rebuild preferences written by :meth:`to_dict`."""
        return cls(data)

    def dimensions(self) -> Iterable[str]:
        """Return the dimension keys with a non-SI choice."""
        return tuple(self._choices)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, UnitSystem) and self._choices == other._choices

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return f"UnitSystem({self._choices!r})"
