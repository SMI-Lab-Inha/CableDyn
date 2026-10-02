# SPDX-License-Identifier: Apache-2.0
"""Lines and what they own: sections, end connections, attachments and ancillaries."""

from __future__ import annotations

import itertools

from cabledyn.project.appearance import ModuleGeometry
from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.descriptors import (
    Bool,
    Child,
    Children,
    Choice,
    Integer,
    OptionalColour,
    OptionalQuantity,
    OptionalVec3,
    Quantity,
    Ref,
    Role,
    Text,
    Vec3,
)
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.points import LineEndTarget
from cabledyn.project.types import LineType
from cabledyn.project.units import (
    ANGLE,
    AREA,
    DIMENSIONLESS,
    LENGTH,
    MASS,
    ROTATIONAL_STIFFNESS,
    VOLUME,
)

__all__ = [
    "Attachment",
    "BendRestrictor",
    "BendStiffener",
    "BuoyancyModule",
    "ClumpWeight",
    "EndConnection",
    "Line",
    "LineAncillary",
    "Section",
    "TouchdownProtection",
]

_A = Role.APPEARANCE


@model_type("section")
class Section(ModelObject):
    """A length of one line type (one deck ``SECTIONS`` row), from End A."""

    type_label = "Section"
    abstract = False
    line_type = Ref(LineType, group="Section")
    length = Quantity(LENGTH, 100.0, limit="section_length", group="Section", doc="Unstretched.")
    segments = Integer(10, limit="segments", group="Section", doc="Number of elements.")
    colour = OptionalColour(group="Appearance", doc="Overrides the line type's colour.")


@model_type("end_connection")
class EndConnection(ModelObject):
    """The bending (and optional torsion) boundary at one end of a finite-EI line.

    One deck ``END CONNECTIONS`` row. ``rotation`` ``stiffness`` uses
    ``rotational_stiffness``; ``pinned`` and ``rigid`` are the deck
    keywords. ``torsion`` ``none`` writes the six-column row.
    """

    type_label = "End connection"
    abstract = False
    end = Choice(("A", "B"), group="Connection")
    rotation = Choice(("stiffness", "pinned", "rigid"), group="Bending")
    rotational_stiffness = Quantity(ROTATIONAL_STIFFNESS, 0.0, minimum=0.0, group="Bending")
    direction = Vec3(
        DIMENSIONLESS, (0.0, 0.0, -1.0), group="Bending", doc="Reference direction, End A to B."
    )
    torsion = Choice(("none", "stiffness", "free", "rigid"), group="Torsion")
    torsional_stiffness = Quantity(ROTATIONAL_STIFFNESS, 0.0, minimum=0.0, group="Torsion")
    normal = OptionalVec3(DIMENSIONLESS, group="Torsion", doc="Zero-twist reference normal.")
    pretwist = OptionalQuantity(ANGLE, group="Torsion", doc="Roll of the end frame.")

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if self.direction == (0.0, 0.0, 0.0):
            found.append(Issue(Severity.ERROR, "must not be zero", self, "direction", 2))
        if self.torsion == "none":
            if self.normal is not None or self.pretwist is not None:
                found.append(
                    Issue(Severity.ERROR, "normal and pretwist need a torsion", self, "torsion", 2)
                )
        elif self.normal is None:
            found.append(Issue(Severity.ERROR, "a torsion needs a normal", self, "normal", 2))
        return found


class Attachment(ModelObject):
    """Base class of discrete items on a finite-EI line (deck ``ATTACHMENTS``).

    The location is an arc length from End A, or a series ``arc_length``,
    ``arc_length + pitch``, ... up to ``last_arc_length``.
    """

    type_label = "Attachment"
    abstract = True
    arc_length = Quantity(LENGTH, 0.0, minimum=0.0, group="Location", doc="From End A.")
    pitch = OptionalQuantity(LENGTH, group="Location", doc="Spacing of a series.")
    last_arc_length = OptionalQuantity(LENGTH, group="Location", doc="Last item of a series.")
    mass = Quantity(MASS, 0.0, group="Properties", doc="Dry mass of one item.")
    volume = Quantity(VOLUME, 0.0, group="Properties", doc="Displaced volume of one item.")
    drag_area = Quantity(AREA, 0.0, group="Properties", doc="Normal drag area.")
    added_mass = Quantity(DIMENSIONLESS, 0.0, group="Properties")
    axial_drag_area = OptionalQuantity(AREA, group="Properties")

    def is_series(self) -> bool:
        """Whether the attachment is a series of identical items."""
        return self.pitch is not None

    def locations(self) -> list[float]:
        """Return the arc length of every item (one value for a single item)."""
        if self.pitch is None or self.last_arc_length is None or self.pitch <= 0.0:
            return [self.arc_length]
        found: list[float] = []
        count = int((self.last_arc_length - self.arc_length) / self.pitch + 1.0e-9) + 1
        for index in range(max(count, 0)):
            found.append(self.arc_length + index * self.pitch)
        return found

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if (self.pitch is None) != (self.last_arc_length is None):
            found.append(
                Issue(Severity.ERROR, "give pitch and last_arc_length together", self, "pitch", 2)
            )
        line = self.parent
        if isinstance(line, Line):
            total = line.length
            for location in self.locations()[:1] + self.locations()[-1:]:
                if location > total:
                    found.append(
                        Issue(
                            Severity.ERROR,
                            f"arc length {location:g} m is beyond the line ({total:g} m)",
                            self,
                            "arc_length",
                            2,
                        )
                    )
                    break
        return found


@model_type("attachment.buoyancy_module")
class BuoyancyModule(Attachment):
    """A buoyancy module (or a series of them)."""

    type_label = "Buoyancy module"
    geometry = Child(ModuleGeometry, role=_A)

    def geometric_volume(self) -> float | None:
        """Return the volume of the drawn shells, or ``None`` if not fully given."""
        g = self.geometry
        if g.outer_diameter is None or g.length is None:
            return None
        bore = g.bore or 0.0
        return 3.141592653589793 / 4.0 * (g.outer_diameter**2 - bore**2) * g.length

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        drawn = self.geometric_volume()
        if drawn is not None and self.volume > 0.0 and abs(drawn - self.volume) > 0.1 * self.volume:
            found.append(
                Issue(
                    Severity.WARNING,
                    f"the drawn module volume {drawn:.3g} m^3 differs from the physics "
                    f"volume {self.volume:.3g} m^3 by more than 10 %",
                    self,
                    "geometry",
                    2,
                )
            )
        return found


@model_type("attachment.clump_weight")
class ClumpWeight(Attachment):
    """A clump weight (or a series of them)."""

    type_label = "Clump weight"
    style = Choice(("cylinder", "chain_clump"), role=_A, group="Appearance")
    clump_diameter = OptionalQuantity(LENGTH, minimum=0.0, role=_A, group="Appearance")
    clump_length = OptionalQuantity(LENGTH, minimum=0.0, role=_A, group="Appearance")


class LineAncillary(ModelObject):
    """Base class of line ancillaries drawn over an arc-length range.

    Ancillaries are cosmetic: the stiffness and mass they add are modelled
    with sections and attachments. They have no deck rows.
    """

    type_label = "Ancillary"
    abstract = True
    start = Quantity(LENGTH, 0.0, minimum=0.0, role=_A, group="Location", doc="From End A.")
    colour = OptionalColour(group="Appearance")

    def extent(self) -> tuple[float, float]:
        """Return the covered arc-length range ``(start, end)``."""
        raise NotImplementedError

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        low, high = self.extent()
        if high < low:
            found.append(Issue(Severity.ERROR, "ends before it starts", self, "start", 2))
        line = self.parent
        if isinstance(line, Line) and high > line.length:
            found.append(Issue(Severity.ERROR, "extends beyond the line", self, "start", 2))
        return found


@model_type("ancillary.bend_stiffener")
class BendStiffener(LineAncillary):
    """A conical bend stiffener at a hang-off."""

    type_label = "Bend stiffener"
    length = Quantity(LENGTH, 5.0, minimum=0.0, role=_A, group="Shape")
    root_diameter = Quantity(LENGTH, 0.8, minimum=0.0, role=_A, group="Shape")
    tip_diameter = Quantity(LENGTH, 0.25, minimum=0.0, role=_A, group="Shape")
    profile = Choice(("linear", "curved"), role=_A, group="Shape")
    flange = Bool(True, role=_A, group="Shape", doc="Steel root flange.")

    def extent(self) -> tuple[float, float]:
        return (self.start, self.start + self.length)


@model_type("ancillary.bend_restrictor")
class BendRestrictor(LineAncillary):
    """Interlocking vertebrae limiting the bend radius."""

    type_label = "Bend restrictor"
    end = Quantity(LENGTH, 5.0, minimum=0.0, role=_A, group="Location")
    element_length = Quantity(LENGTH, 0.3, minimum=0.0, role=_A, group="Shape")
    outer_diameter = Quantity(LENGTH, 0.5, minimum=0.0, role=_A, group="Shape")
    locking_radius = OptionalQuantity(LENGTH, minimum=0.0, role=_A, group="Shape")

    def extent(self) -> tuple[float, float]:
        return (self.start, self.end)


@model_type("ancillary.touchdown_protection")
class TouchdownProtection(LineAncillary):
    """Half-shell protection sleeves over the touchdown zone."""

    type_label = "Touchdown protection"
    end = Quantity(LENGTH, 20.0, minimum=0.0, role=_A, group="Location")
    segment_length = Quantity(LENGTH, 1.0, minimum=0.0, role=_A, group="Shape")
    outer_diameter = Quantity(LENGTH, 0.4, minimum=0.0, role=_A, group="Shape")

    def extent(self) -> tuple[float, float]:
        return (self.start, self.end)


@model_type("line")
class Line(ModelObject):
    """A line from End A (fairlead side) to End B (anchor side).

    One deck ``LINES`` row with its ordered sections.
    """

    type_label = "Line"
    abstract = False
    end_a = Ref(LineEndTarget, group="Ends", doc="End A attachment.")
    end_b = Ref(LineEndTarget, group="Ends", doc="End B attachment.")
    sections = Children(Section, min_items=1, group="Sections")
    end_connections = Children(EndConnection, group="Ends")
    attachments = Children(Attachment, group="Attachments")
    ancillaries = Children(LineAncillary, role=_A, group="Ancillaries")
    outputs = Text("-", group="Outputs", doc="Output flags: '-', or p, t, r combined.")
    colour = OptionalColour(group="Appearance")

    @property
    def length(self) -> float:
        """Total unstretched length, the sum of the section lengths."""
        return float(sum(section.length for section in self.sections))

    @property
    def segment_count(self) -> int:
        """Total number of elements."""
        return sum(section.segments for section in self.sections)

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if self.end_a is not None and self.end_a is self.end_b:
            found.append(Issue(Severity.ERROR, "both ends attach to one object", self, "end_b", 2))
        seen: set[str] = set()
        for row in self.end_connections:
            if row.end in seen:
                found.append(
                    Issue(Severity.ERROR, f"End {row.end} has two end connections", row, "end", 2)
                )
            seen.add(row.end)
        ranges = sorted((item.extent(), item.uid) for item in self.ancillaries)
        for (first, _), (second, uid) in itertools.pairwise(ranges):
            if second[0] < first[1]:
                other = next(item for item in self.ancillaries if item.uid == uid)
                found.append(Issue(Severity.ERROR, "overlaps another ancillary", other, "start", 2))
        return found
