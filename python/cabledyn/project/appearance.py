# SPDX-License-Identifier: Apache-2.0
"""Appearance sub-objects: cosmetic data for the 3D views.

Every property here has the ``appearance`` role: the deck writer never reads
it. The presentation view generates geometry from these properties together
with the physics properties of the owning object (for example a chain's
nominal diameter and its link proportions). ``None`` means "derive from the
physics data".
"""

from __future__ import annotations

from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.descriptors import (
    Bool,
    Choice,
    Colour,
    FilePath,
    Integer,
    OptionalColour,
    OptionalQuantity,
    OptionalText,
    Quantity,
    Role,
    Text,
    Vec3,
)
from cabledyn.project.units import ANGLE, DIMENSIONLESS, LENGTH

__all__ = [
    "CableAppearance",
    "ChainAppearance",
    "EnvironmentAppearance",
    "FibreRopeAppearance",
    "LineAppearance",
    "MeshAsset",
    "ModuleGeometry",
    "PlainLineAppearance",
    "SeabedAppearance",
    "WireRopeAppearance",
]

_A = Role.APPEARANCE


class LineAppearance(ModelObject):
    """Base class of the surface descriptions of a line type."""

    type_label = "Line appearance"
    abstract = True


@model_type("appearance.line.plain")
class PlainLineAppearance(LineAppearance):
    """A smooth tube in the line type's colour."""

    type_label = "Plain tube"
    roughness = Quantity(DIMENSIONLESS, 0.5, minimum=0.0, maximum=1.0, role=_A, group="Surface")


@model_type("appearance.line.chain")
class ChainAppearance(LineAppearance):
    """Chain links generated from the nominal diameter and link proportions."""

    type_label = "Chain links"
    finish = Choice(("black", "galvanised", "rusted", "fouled"), role=_A, group="Surface")
    link_length = OptionalQuantity(
        LENGTH, role=_A, group="Links", doc="Outer link length; derived from the diameter if unset."
    )
    link_width = OptionalQuantity(
        LENGTH, role=_A, group="Links", doc="Outer link width; derived from the diameter if unset."
    )
    section_end_links = Bool(
        False, role=_A, group="Links", doc="Draw enlarged links at section ends."
    )


@model_type("appearance.line.wire_rope")
class WireRopeAppearance(LineAppearance):
    """Helical strands, optionally under a polymer sheath."""

    type_label = "Wire rope"
    lay_length = OptionalQuantity(LENGTH, role=_A, group="Lay", doc="Strand lay length.")
    lay_direction = Choice(("right", "left"), role=_A, group="Lay")
    finish = Choice(("galvanised", "bright"), role=_A, group="Surface")
    sheath_colour = OptionalColour(role=_A, group="Sheath", doc="Unset: no sheath.")
    sheath_thickness = OptionalQuantity(LENGTH, role=_A, group="Sheath")


@model_type("appearance.line.fibre_rope")
class FibreRopeAppearance(LineAppearance):
    """A braided jacket with optional marker stripes."""

    type_label = "Fibre rope"
    jacket_colour = Colour("#f2f2f2", group="Jacket")
    marker_stripes = Integer(1, minimum=0, role=_A, group="Markers")
    marker_colour = Colour("#1f4e9c", group="Markers")
    lay_length = OptionalQuantity(LENGTH, role=_A, group="Markers", doc="Stripe lay length.")


@model_type("appearance.line.cable")
class CableAppearance(LineAppearance):
    """A served power-cable exterior with stripes, marker bands and an optional cutaway."""

    type_label = "Cable serving"
    serving_colour = Colour("#1a1a1a", group="Serving")
    lay_angle = Quantity(ANGLE, 15.0, minimum=0.0, maximum=90.0, role=_A, group="Serving")
    stripe_count = Integer(1, minimum=0, role=_A, group="Stripes")
    stripe_colour = Colour("#f2c200", group="Stripes")
    stripe_width = OptionalQuantity(LENGTH, role=_A, group="Stripes")
    marker_band_spacing = OptionalQuantity(
        LENGTH, role=_A, group="Markers", doc="Spacing of marker bands; unset for none."
    )
    printed_id = Text("", role=_A, group="Markers", doc="Identification text printed along it.")
    cutaway = Bool(False, role=_A, group="Close-up", doc="Show cores and armour in close-ups.")


@model_type("appearance.module_geometry")
class ModuleGeometry(ModelObject):
    """The shape of a buoyancy module (two half-shells around the line)."""

    type_label = "Module geometry"
    abstract = False
    outer_diameter = OptionalQuantity(
        LENGTH, role=_A, group="Shape", doc="Unset: derived from the module volume."
    )
    length = OptionalQuantity(LENGTH, role=_A, group="Shape")
    bore = OptionalQuantity(LENGTH, role=_A, group="Shape", doc="Unset: the line diameter.")
    end_radius = OptionalQuantity(LENGTH, role=_A, group="Shape")
    seam = Bool(True, role=_A, group="Detail", doc="Draw the half-shell seam.")
    straps = Integer(2, minimum=0, role=_A, group="Detail")
    strap_width = OptionalQuantity(LENGTH, role=_A, group="Detail")
    colour = Colour("#f2a900", group="Surface")


@model_type("appearance.mesh")
class MeshAsset(ModelObject):
    """A 3D mesh drawn for a body or a floater type (glTF first; OBJ, STL)."""

    type_label = "Mesh"
    abstract = False
    file = FilePath(role=_A, group="File", doc="Mesh file; unset for a generic shape.")
    generic_shape = Choice(
        ("box", "cylinder", "semi_submersible", "spar", "barge", "buoy"),
        role=_A,
        group="File",
        doc="Built-in shape drawn when no file is given.",
    )
    scale = Quantity(DIMENSIONLESS, 1.0, minimum=0.0, role=_A, group="Placement")
    origin = Vec3(LENGTH, role=_A, group="Placement", doc="Mesh origin in body axes.")
    orientation = Vec3(ANGLE, role=_A, group="Placement", doc="Mesh rotation in body axes.")
    file_units = Choice(("m", "mm", "cm", "ft", "in"), role=_A, group="Placement")
    material = OptionalText(role=_A, group="Surface", doc="Material preset; unset: the file's.")
    colour = OptionalColour(role=_A, group="Surface")
    waterline = Bool(False, role=_A, group="Detail", doc="Mark the still-water line.")


@model_type("appearance.environment")
class EnvironmentAppearance(ModelObject):
    """Sea, sky and underwater look of the presentation view."""

    type_label = "Environment appearance"
    abstract = False
    sea_colour = Colour("#1d4f6e", group="Sea")
    clarity = Quantity(DIMENSIONLESS, 0.5, minimum=0.0, maximum=1.0, role=_A, group="Sea")
    ripples = Quantity(DIMENSIONLESS, 0.5, minimum=0.0, maximum=1.0, role=_A, group="Sea")
    sun_azimuth = Quantity(ANGLE, 200.0, role=_A, group="Sky")
    sun_elevation = Quantity(ANGLE, 30.0, minimum=-90.0, maximum=90.0, role=_A, group="Sky")
    sky_map = FilePath(role=_A, group="Sky", doc="Environment map; unset for the default sky.")
    underwater_fog = Quantity(
        DIMENSIONLESS, 0.5, minimum=0.0, maximum=1.0, role=_A, group="Underwater"
    )


@model_type("appearance.seabed")
class SeabedAppearance(ModelObject):
    """Seabed texture and colouring of the presentation view."""

    type_label = "Seabed appearance"
    abstract = False
    texture = Choice(("sand", "silt", "gravel", "rock"), role=_A, group="Surface")
    colour = Colour("#9c8a6a", group="Surface")
    depth_colouring = Bool(False, role=_A, group="Surface", doc="Colour by depth.")
    trench_exaggeration = Quantity(
        DIMENSIONLESS, 1.0, minimum=0.0, role=_A, group="Detail", doc="Cosmetic, labelled."
    )
