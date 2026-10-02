# SPDX-License-Identifier: Apache-2.0
"""Units layer of the project object model."""

from __future__ import annotations

import math

import pytest

from cabledyn.project.units import (
    ANGLE,
    DIMENSIONLESS,
    FORCE,
    FORCE_PER_LENGTH,
    FREQUENCY,
    LENGTH,
    MASS,
    STIFFNESS,
    UNITS,
    UnitError,
    UnitSystem,
    convert,
    parse_quantity,
    unit,
    units_for,
)


def test_convert_between_units_of_one_dimension() -> None:
    assert convert(2.5, "kN", "N") == 2500.0
    assert convert(1.0, "ft", "m") == pytest.approx(0.3048)
    assert convert(180.0, "deg", "rad") == pytest.approx(math.pi)
    assert convert(1.0, "Hz", "rad/s") == pytest.approx(2.0 * math.pi)
    assert convert(3.0, "m", "m") == 3.0


def test_convert_refuses_mixed_dimensions_and_unknown_units() -> None:
    with pytest.raises(UnitError, match="cannot convert"):
        convert(1.0, "m", "kg")
    with pytest.raises(UnitError, match="unknown unit"):
        convert(1.0, "furlong", "m")


def test_shared_symbols_need_a_dimension() -> None:
    with pytest.raises(UnitError, match="ambiguous"):
        unit("kN/m")
    assert unit("kN/m", STIFFNESS).dimension == STIFFNESS
    assert unit("kN/m", FORCE_PER_LENGTH).scale == 1.0e3
    assert convert(2.0, "kN/m", "N/m", STIFFNESS) == 2000.0
    with pytest.raises(UnitError, match="for mass"):
        unit("m", MASS)


def test_unit_round_trip_and_listing() -> None:
    tonne = unit("t")
    assert tonne.to_si(2.0) == 2000.0
    assert tonne.from_si(2000.0) == 2.0
    lengths = units_for(LENGTH)
    assert lengths[0].symbol == "m"
    assert {item.symbol for item in lengths} >= {"mm", "km", "ft"}
    assert all(item.dimension in {u.dimension for u in UNITS} for item in UNITS)


def test_parse_quantity() -> None:
    assert parse_quantity("850 m", LENGTH) == 850.0
    assert parse_quantity("1.2 MN", FORCE) == pytest.approx(1.2e6)
    assert parse_quantity("12", LENGTH, "mm") == pytest.approx(0.012)
    assert parse_quantity("-3e2", LENGTH) == -300.0
    assert parse_quantity("0.5", DIMENSIONLESS) == 0.5
    with pytest.raises(UnitError, match="not a number"):
        parse_quantity("about ten", LENGTH)
    with pytest.raises(UnitError):
        parse_quantity("10 kg", LENGTH)


def test_unit_system_display_and_preferences() -> None:
    si = UnitSystem.si()
    assert si.format(850.0, LENGTH) == "850 m"
    assert si.format(0.3, DIMENSIONLESS) == "0.3"
    eng = UnitSystem.engineering()
    assert eng.format(1.5e6, FORCE) == "1500 kN"
    assert eng.to_display(2000.0, MASS) == 2.0
    assert eng.from_display(2.0, MASS) == 2000.0
    assert eng.unit_for(ANGLE).symbol == "deg"
    assert "force" in eng.dimensions()
    copy = UnitSystem.from_dict(eng.to_dict())
    assert copy == eng
    assert copy != si
    assert repr(copy).startswith("UnitSystem(")
    with pytest.raises(UnitError, match="unknown dimension"):
        si.set("colour", "m")
    with pytest.raises(UnitError):
        si.set("length", "kg")
    si.set("frequency", "rad/s")
    assert si.to_display(1.0, FREQUENCY) == pytest.approx(2.0 * math.pi)
