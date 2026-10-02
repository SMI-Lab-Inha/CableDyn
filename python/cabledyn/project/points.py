# SPDX-License-Identifier: Apache-2.0
"""Points and the things a line end can attach to.

A line end refers to a :class:`LineEndTarget`: a :class:`Point` of any kind,
or one end of a rod (:class:`RodEndpoint`). The deck's ``POINTS`` types map
to the point classes; ``Coupled`` and ``Vessel`` points are
:class:`FloaterPoint` objects.
"""

from __future__ import annotations

from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.bodies import Body, Turbine
from cabledyn.project.descriptors import (
    Choice,
    OptionalColour,
    OptionalInteger,
    OptionalQuantity,
    Quantity,
    Ref,
    Role,
    Vec3,
)
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.units import AREA, DIMENSIONLESS, LENGTH, MASS, VOLUME

__all__ = [
    "BodyPoint",
    "ConnectPoint",
    "FixedPoint",
    "FloaterPoint",
    "FreePoint",
    "LineEndTarget",
    "Point",
    "RodEndpoint",
    "RodPoint",
    "TurbinePoint",
]

_A = Role.APPEARANCE


class LineEndTarget(ModelObject):
    """Base class of everything a line end can attach to."""

    type_label = "Attachment target"
    abstract = True


@model_type("rod.endpoint")
class RodEndpoint(LineEndTarget):
    """End A or End B of a rod; owned by the rod, never removed on its own."""

    type_label = "Rod end"
    abstract = False
    end = Choice(("A", "B"), read_only=True, role=Role.META, group="Identity")

    @classmethod
    def make(cls, end: str) -> RodEndpoint:
        """Return a new endpoint for ``end`` (``"A"`` or ``"B"``)."""
        item = cls(f"End {end}")
        item._values["end"] = type(item).end.coerce(end)
        return item

    def label(self) -> str:
        rod = self.parent
        owner = "" if rod is None else f"{rod.label()} "
        return f"{owner}End {self.end}"


class Point(LineEndTarget):
    """Base class of the points: position and lumped hydrodynamic properties."""

    type_label = "Point"
    abstract = True
    position = Vec3(LENGTH, limit="coordinate", group="Position")
    mass = Quantity(MASS, 0.0, limit="coefficient", group="Properties")
    volume = Quantity(VOLUME, 0.0, limit="coefficient", group="Properties", doc="Displaced volume.")
    drag_area = Quantity(AREA, 0.0, limit="coefficient", group="Properties", doc="CdA.")
    added_mass = Quantity(DIMENSIONLESS, 0.0, limit="small_coefficient", group="Properties")
    connector = Choice(
        (
            "none",
            "shackle",
            "h_link",
            "triplate",
            "subsea_connector",
            "bell_mouth",
            "hang_off",
        ),
        role=_A,
        group="Appearance",
        doc="Hardware drawn at the point.",
    )
    colour = OptionalColour(group="Appearance")


@model_type("point.fixed")
class FixedPoint(Point):
    """A point fixed in space, such as an anchor (deck ``Fixed``)."""

    type_label = "Fixed point"
    anchor = Choice(
        ("none", "drag", "suction_pile", "driven_pile", "gravity"),
        role=_A,
        group="Appearance",
        doc="Anchor drawn at the point.",
    )
    anchor_size = OptionalQuantity(LENGTH, minimum=0.0, role=_A, group="Appearance")


@model_type("point.floater")
class FloaterPoint(Point):
    """A point moved with the floater (deck ``Coupled`` or ``Vessel``)."""

    type_label = "Floater point"
    kind = Choice(("coupled", "floater"), group="Motion", doc="Deck point type Coupled or Vessel.")


@model_type("point.free")
class FreePoint(Point):
    """A free point solved for (deck ``Free``)."""

    type_label = "Free point"


@model_type("point.connect")
class ConnectPoint(Point):
    """A free connection point joining several lines (deck ``Connect``)."""

    type_label = "Connection point"


@model_type("point.turbine")
class TurbinePoint(Point):
    """A point carried by a FAST.Farm turbine (deck ``Turbine<J>``/``T<J>``).

    The point refers to its :class:`~cabledyn.project.Turbine` object, so it
    follows a renumbered turbine and is removed with it. In a coupled farm
    deck without a ``TURBINES`` table the turbine is known only by its number,
    given in ``turbine_number`` instead. The position is relative to the
    turbine's reference.
    """

    type_label = "Turbine point"
    turbine = Ref(Turbine, required=False, group="Farm", doc="The turbine.")
    turbine_number = OptionalInteger(
        minimum=1,
        maximum=10000,
        group="Farm",
        doc="Turbine number J when the deck has no TURBINES row for it.",
    )

    def number(self) -> int | None:
        """Return the turbine number J written to the deck, or ``None``."""
        return self.turbine.number if self.turbine is not None else self.turbine_number

    def required_references(self) -> frozenset[str]:
        found = super().required_references()
        return found | {"turbine"} if self.turbine_number is None else found

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if (self.turbine is None) == (self.turbine_number is None):
            found.append(
                Issue(
                    Severity.ERROR,
                    "give either the turbine or a turbine number",
                    self,
                    "turbine",
                    2,
                )
            )
        return found


@model_type("point.body")
class BodyPoint(Point):
    """A point fixed to a body; the position is in body axes (deck ``Body<N>``)."""

    type_label = "Body point"
    body = Ref(Body, group="Attachment")


@model_type("point.rod")
class RodPoint(Point):
    """A point at a rod end (deck ``Rod<N>A``/``Rod<N>B``)."""

    type_label = "Rod-end point"
    rod_end = Ref(RodEndpoint, group="Attachment")
