# SPDX-License-Identifier: Apache-2.0
"""Validation issues reported by the three validation layers."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cabledyn.project.base import ModelObject

__all__ = ["Issue", "Severity"]


class Severity(str, Enum):
    """How serious an :class:`Issue` is."""

    ERROR = "error"
    """The model cannot be written as a valid deck, or is inconsistent."""
    WARNING = "warning"
    """The model is valid but probably not what was meant."""
    INFO = "info"
    """A note, for example an option kept verbatim."""


@dataclass(frozen=True)
class Issue:
    """One validation finding, attached to an object and optionally a property.

    Attributes
    ----------
    severity : Severity
        Error, warning, or note.
    message : str
        Human-readable message.
    obj : ModelObject | None
        The object concerned, or ``None`` for a model-wide finding.
    prop : str | None
        The property concerned, if any.
    layer : int
        ``1`` property checks, ``2`` object invariants, ``3`` native deck rules.
    """

    severity: Severity
    message: str
    obj: ModelObject | None = None
    prop: str | None = None
    layer: int = 1

    @property
    def is_error(self) -> bool:
        """Whether the issue is an error."""
        return self.severity is Severity.ERROR

    def __str__(self) -> str:
        where = "" if self.obj is None else f"{self.obj.label()}: "
        field = "" if self.prop is None else f"{self.prop}: "
        return f"{self.severity.value}: {where}{field}{self.message}"
