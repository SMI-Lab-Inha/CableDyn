# SPDX-License-Identifier: Apache-2.0
"""Rods: rigid cylinders with two ends that lines and points can attach to."""

from __future__ import annotations

from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.bodies import Body
from cabledyn.project.descriptors import (
    Child,
    Choice,
    Integer,
    OptionalColour,
    Ref,
    Text,
    Vec3,
)
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.points import RodEndpoint
from cabledyn.project.types import RodType
from cabledyn.project.units import LENGTH

__all__ = ["Rod"]

_BODY_KINDS = frozenset({"body", "body_pinned"})


@model_type("rod")
class Rod(ModelObject):
    """A rod (one deck ``RODS`` row).

    ``attachment`` ``coupled`` and ``floater`` correspond to the deck
    types ``Coupled`` and ``Vessel``; ``body`` and ``body_pinned`` need
    ``body``, and the end positions are then in body axes.
    """

    type_label = "Rod"
    abstract = False
    rod_type = Ref(RodType, group="Type")
    attachment = Choice(
        ("free", "fixed", "pinned", "coupled", "floater", "body", "body_pinned"),
        group="Attachment",
    )
    body = Ref(Body, required=False, group="Attachment")
    position_a = Vec3(LENGTH, limit="coordinate", group="Geometry", doc="End A position.")
    position_b = Vec3(
        LENGTH, (0.0, 0.0, -10.0), limit="coordinate", group="Geometry", doc="End B position."
    )
    segments = Integer(1, minimum=0, group="Discretisation", doc="0 for a zero-length rod.")
    outputs = Text("-", group="Outputs", doc="Output flags: '-' or 'p'.")
    end_a = Child(RodEndpoint, lambda: RodEndpoint.make("A"), read_only=True)
    end_b = Child(RodEndpoint, lambda: RodEndpoint.make("B"), read_only=True)
    colour = OptionalColour(group="Appearance")

    def required_references(self) -> frozenset[str]:
        """A body or body_pinned rod cannot do without its body."""
        found = super().required_references()
        return found | {"body"} if self.attachment in _BODY_KINDS else found

    def endpoint(self, end: str) -> RodEndpoint:
        """Return the endpoint object of ``end`` (``"A"`` or ``"B"``)."""
        if end.upper() not in {"A", "B"}:
            raise ValueError(f"end must be A or B, got {end!r}")
        return self.end_a if end.upper() == "A" else self.end_b

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if (self.attachment in _BODY_KINDS) != (self.body is not None):
            found.append(
                Issue(
                    Severity.ERROR,
                    "a body or body_pinned rod needs its body, and only those have one",
                    self,
                    "body",
                    2,
                )
            )
        return found
