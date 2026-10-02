# SPDX-License-Identifier: Apache-2.0
"""Numeric deck limits shared by the property checks and the deck reader.

The limits are not restated here: each entry refers to the admissible-range
constant that :class:`cabledyn.DeckFile` applies, so a property check
(validation layer 1) and the native deck rules (layer 3) cannot drift apart.
Layer 1 is an early, per-field hint; the deck rules remain the authority.

The table is the hook for a machine-readable deck schema: when the deck
reader takes its limits from a schema file, :data:`DECK_LIMITS` is the one
place to read it from.
"""

from __future__ import annotations

from dataclasses import dataclass

from cabledyn import deck_file as _deck

__all__ = ["DECK_LIMITS", "DeckLimit", "deck_limit"]


@dataclass(frozen=True)
class DeckLimit:
    """An admissible range of one deck field.

    Attributes
    ----------
    minimum, maximum : float | None
        Inclusive bounds, or ``None`` for no bound.
    absolute : bool
        Apply the bounds to the magnitude of the value.
    exclusive_minimum : bool
        The minimum itself is not admissible.
    """

    minimum: float | None = None
    maximum: float | None = None
    absolute: bool = False
    exclusive_minimum: bool = False

    def violation(self, value: float) -> str | None:
        """Return a message if ``value`` is outside the range, else ``None``."""
        checked = abs(value) if self.absolute else value
        what = "magnitude" if self.absolute else "value"
        if self.minimum is not None:
            low = checked <= self.minimum if self.exclusive_minimum else checked < self.minimum
            if low:
                relation = "greater than" if self.exclusive_minimum else "at least"
                return f"{what} must be {relation} {self.minimum:g}"
        if self.maximum is not None and checked > self.maximum:
            return f"{what} must be at most {self.maximum:g}"
        return None


_PERIOD_LOW, _PERIOD_HIGH = _deck._NATIVE_WAVE_PERIOD_RANGE

DECK_LIMITS: dict[str, DeckLimit] = {
    "coordinate": DeckLimit(maximum=_deck._NATIVE_MAX_COORDINATE, absolute=True),
    "coefficient": DeckLimit(maximum=_deck._NATIVE_MAX_COEFFICIENT, absolute=True),
    "small_coefficient": DeckLimit(maximum=_deck._NATIVE_MAX_SMALL_COEFFICIENT, absolute=True),
    "gravity": DeckLimit(maximum=_deck._NATIVE_MAX_GRAVITY),
    "water_density": DeckLimit(maximum=_deck._NATIVE_MAX_WATER_DENSITY),
    "water_depth": DeckLimit(maximum=_deck._NATIVE_MAX_COORDINATE),
    "line_diameter": DeckLimit(
        minimum=_deck._NATIVE_MIN_DIAMETER, maximum=_deck._NATIVE_MAX_DIAMETER
    ),
    "line_mass": DeckLimit(maximum=_deck._NATIVE_MAX_ABS_MASS, absolute=True),
    "axial_stiffness": DeckLimit(maximum=_deck._NATIVE_MAX_EA),
    "section_length": DeckLimit(
        minimum=0.0, maximum=_deck._NATIVE_MAX_SECTION_LENGTH, exclusive_minimum=True
    ),
    "segments": DeckLimit(minimum=1, maximum=_deck._NATIVE_MAX_SEGMENTS),
    "wave_height": DeckLimit(
        minimum=0.0, maximum=_deck._NATIVE_MAX_WAVE_HEIGHT, exclusive_minimum=True
    ),
    "wave_period": DeckLimit(minimum=_PERIOD_LOW, maximum=_PERIOD_HIGH),
    "fluid_speed": DeckLimit(maximum=_deck._NATIVE_MAX_FLUID_SPEED, absolute=True),
}
"""Shared limits by key, each taken from the native deck reader's constant."""


def deck_limit(key: str) -> DeckLimit:
    """Return the shared limit ``key``.

    Parameters
    ----------
    key : str
        A key of :data:`DECK_LIMITS`.

    Returns
    -------
    DeckLimit
        The limit.

    Raises
    ------
    KeyError
        If the key is unknown.
    """
    return DECK_LIMITS[key]
