# SPDX-License-Identifier: Apache-2.0
"""Bodies and FAST.Farm turbines.

A body is one deck ``BODIES`` row. The deck's row types map to classes:
``Rigid6`` (and the MoorDyn ``Free`` row) to :class:`Rigid6Body`, ``Point3``
to :class:`Point3Body`, and ``Coupled``/``Vessel`` (with the MoorDyn aliases)
to :class:`Floater`. The term *floater* is used throughout the object model;
the deck keeps its own ``Vessel`` keyword, which the deck adapter maps.
"""

from __future__ import annotations

from cabledyn.project.appearance import MeshAsset
from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.descriptors import (
    Child,
    Choice,
    Integer,
    OptionalColour,
    Quantity,
    Ref,
    Role,
    Vec3,
    Vector,
)
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.types import FloaterType
from cabledyn.project.units import (
    ANGLE,
    AREA,
    DIMENSIONLESS,
    INERTIA,
    LENGTH,
    MASS,
    ROTATIONAL_STIFFNESS,
    STIFFNESS,
    VOLUME,
)

__all__ = ["Body", "Floater", "Point3Body", "Rigid6Body", "Turbine"]

_ROW_FORMATS = ("cabledyn", "moordyn")


class Body(ModelObject):
    """Base class of the bodies: pose, mass properties and hydrodynamics.

    The deck has two ``BODIES`` row layouts. The CableDyn layout carries the
    hydrostatic restoring terms and three optional inertias; the MoorDyn
    layout carries a centre of gravity and multi-valued drag and added-mass
    entries. ``row_format`` records which one the body uses.
    """

    type_label = "Body"
    abstract = True
    position = Vec3(LENGTH, limit="coordinate", group="Pose", doc="Reference position.")
    orientation = Vec3(ANGLE, group="Pose", doc="Roll, pitch, yaw.")
    mass = Quantity(MASS, 0.0, group="Mass")
    volume = Quantity(VOLUME, 0.0, group="Hydrostatics", doc="Displaced volume.")
    centre_of_gravity = Vector(
        LENGTH, (0.0,), lengths=(1, 3), group="Mass", doc="CG z, or (x, y, z); MoorDyn layout."
    )
    inertia = Vector(
        INERTIA, (), lengths=(0, 1, 3), group="Mass", doc="Ixx, Iyy, Izz (one value: all three)."
    )
    heave_stiffness = Quantity(STIFFNESS, 0.0, group="Hydrostatics", doc="C33; CableDyn layout.")
    roll_stiffness = Quantity(ROTATIONAL_STIFFNESS, 0.0, group="Hydrostatics", doc="C44.")
    pitch_stiffness = Quantity(ROTATIONAL_STIFFNESS, 0.0, group="Hydrostatics", doc="C55.")
    drag_area = Vector(AREA, (0.0,), lengths=(1, 2, 3, 6), group="Hydrodynamics", doc="CdA.")
    added_mass = Vector(DIMENSIONLESS, (0.0,), lengths=(1, 3), group="Hydrodynamics", doc="Ca.")
    row_format = Choice(
        _ROW_FORMATS, role=Role.META, group="Deck", doc="BODIES row layout used by the deck."
    )
    colour = OptionalColour(group="Appearance")
    mesh = Child(MeshAsset, role=Role.APPEARANCE)

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if self.row_format == "cabledyn":
            for name in ("drag_area", "added_mass"):
                if len(self.get_value(name)) != 1:
                    found.append(
                        Issue(Severity.ERROR, "the CableDyn layout takes one value", self, name, 2)
                    )
            if len(self.inertia) not in (0, 3):
                found.append(
                    Issue(Severity.ERROR, "the CableDyn layout takes 0 or 3", self, "inertia", 2)
                )
        elif not self.inertia:
            found.append(Issue(Severity.ERROR, "the MoorDyn layout needs it", self, "inertia", 2))
        return found


@model_type("body.rigid6")
class Rigid6Body(Body):
    """A free six-degree-of-freedom rigid body (deck ``Rigid6``, MoorDyn ``Free``)."""

    type_label = "Rigid body"


@model_type("body.point3")
class Point3Body(Body):
    """A three-degree-of-freedom point body (deck ``Point3``)."""

    type_label = "Point body"


@model_type("body.floater")
class Floater(Body):
    """A floater moved by the coupling host or by the prescribed motion.

    The deck body types ``Coupled`` and ``Vessel`` map to ``kind``
    ``"coupled"`` and ``"floater"``; the solver moves both kinds alike.
    """

    type_label = "Floater"
    kind = Choice(("coupled", "floater"), group="Motion", doc="Deck body type Coupled or Vessel.")
    floater_type = Ref(FloaterType, required=False, role=Role.CONSTRUCTION, group="Type")


@model_type("turbine")
class Turbine(ModelObject):
    """A FAST.Farm turbine of a standalone farm deck (``TURBINES`` row)."""

    type_label = "Turbine"
    abstract = False
    number = Integer(1, minimum=1, maximum=10000, group="Farm", doc="Turbine number J.")
    position = Vec3(LENGTH, limit="coordinate", group="Farm", doc="Farm-global reference.")
    platform_displacement = Vector(
        DIMENSIONLESS,
        (),
        lengths=(0, 6),
        limit="coordinate",
        group="Farm",
        doc="Initial surge, sway, heave (m) and roll, pitch, yaw (deg); empty for none.",
    )
