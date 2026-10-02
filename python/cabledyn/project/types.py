# SPDX-License-Identifier: Apache-2.0
"""Type libraries: line types, rod types and floater types.

A line type carries three kinds of data on one object, kept apart by the
property roles: physics written to the deck (diameter, mass, the axial and
bending models, hydrodynamic coefficients), construction data (chain grade,
rope construction, cable build; not written to the deck today), and the
appearance sub-object used by the presentation view.
"""

from __future__ import annotations

from cabledyn.project.appearance import (
    CableAppearance,
    ChainAppearance,
    FibreRopeAppearance,
    LineAppearance,
    MeshAsset,
    PlainLineAppearance,
    WireRopeAppearance,
)
from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.descriptors import (
    Child,
    Choice,
    Colour,
    FilePath,
    Integer,
    OptionalQuantity,
    OptionalText,
    OptionalVec3,
    Quantity,
    Role,
    Strategy,
    Vector,
)
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.units import (
    AREA,
    AXIAL_DAMPING,
    BENDING_STIFFNESS,
    DIMENSIONLESS,
    FORCE,
    INERTIA,
    LENGTH,
    MASS,
    MASS_PER_LENGTH,
    ROTARY_INERTIA_PER_LENGTH,
    VOLUME,
)

__all__ = [
    "AxialModel",
    "BendingModel",
    "CableType",
    "ChainType",
    "FibreRopeType",
    "FloaterType",
    "GenericLineType",
    "HydroDatabase",
    "HydroDatabaseFile",
    "LineType",
    "LinearAxial",
    "NoHydroDatabase",
    "RodType",
    "SyropeAxial",
    "ViscoelasticAxial",
    "WireRopeType",
]

_C = Role.CONSTRUCTION
_HYDRO = "Hydrodynamics"


# --------------------------------------------------------------------------- axial models


class AxialModel(ModelObject):
    """Base class of the axial (tension-strain) models of a line type."""

    type_label = "Axial model"
    abstract = True
    damping = Vector(
        AXIAL_DAMPING,
        (0.0,),
        lengths=(1, 2),
        limit="coefficient",
        group="Axial",
        doc="Axial damping BA (negative: damping ratio); two values (static, dynamic).",
    )


@model_type("axial.linear")
class LinearAxial(AxialModel):
    """Linear-elastic axial stiffness EA."""

    type_label = "Linear EA"
    stiffness = Quantity(FORCE, 1.0e9, limit="axial_stiffness", group="Axial", doc="EA.")


@model_type("axial.viscoelastic")
class ViscoelasticAxial(AxialModel):
    """Two-spring viscoelastic stiffness: ``Es|Ed`` or ``Es|alphaMBL|vbeta``."""

    type_label = "Viscoelastic EA"
    static_stiffness = Quantity(
        FORCE, 1.0e8, limit="axial_stiffness", group="Axial", doc="Static stiffness Es."
    )
    dynamic_stiffness = OptionalQuantity(
        FORCE, limit="axial_stiffness", group="Axial", doc="Constant dynamic stiffness Ed."
    )
    alpha_mbl = OptionalQuantity(
        DIMENSIONLESS, group="Axial", doc="Load-dependent dynamic stiffness intercept (x MBL)."
    )
    beta = OptionalQuantity(
        DIMENSIONLESS, group="Axial", doc="Load-dependent dynamic stiffness slope."
    )

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        constant = self.dynamic_stiffness is not None
        loaded = (self.alpha_mbl is not None, self.beta is not None)
        if constant == any(loaded) or (any(loaded) and not all(loaded)):
            found.append(
                Issue(
                    Severity.ERROR,
                    "give either dynamic_stiffness or both alpha_mbl and beta",
                    self,
                    None,
                    2,
                )
            )
        return found


@model_type("axial.syrope")
class SyropeAxial(AxialModel):
    """Syrope working-curve stiffness with a tension-dependent fast spring."""

    type_label = "Syrope"
    settings_file = FilePath(group="Axial", doc="Syrope settings file.")
    alpha = Quantity(FORCE, 1.0e8, group="Axial", doc="Fast-spring stiffness intercept.")
    beta = Quantity(DIMENSIONLESS, 20.0, group="Axial", doc="Fast-spring slope on tension.")

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if not self.settings_file:
            found.append(Issue(Severity.ERROR, "needs a settings file", self, "settings_file", 2))
        return found


@model_type("bending")
class BendingModel(ModelObject):
    """Bending, shear and torsion data of a line type (the finite-EI columns)."""

    type_label = "Bending model"
    abstract = False
    bending_stiffness = Quantity(
        BENDING_STIFFNESS, 0.0, minimum=0.0, limit="coefficient", group="Bending", doc="EI."
    )
    shear_stiffness = OptionalQuantity(FORCE, group="Finite EI", doc="GAs.")
    torsional_stiffness = OptionalQuantity(BENDING_STIFFNESS, group="Finite EI", doc="GJ.")
    rotary_inertia_axial = OptionalQuantity(
        ROTARY_INERTIA_PER_LENGTH, group="Finite EI", doc="Axial rotary inertia Irt."
    )
    rotary_inertia_normal = OptionalQuantity(
        ROTARY_INERTIA_PER_LENGTH, group="Finite EI", doc="Normal rotary inertia Irn."
    )

    def extended(self) -> bool:
        """Whether the four optional finite-EI columns are given."""
        return self.shear_stiffness is not None

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        given = [
            value is not None
            for value in (
                self.shear_stiffness,
                self.torsional_stiffness,
                self.rotary_inertia_axial,
                self.rotary_inertia_normal,
            )
        ]
        if any(given) and not all(given):
            found.append(
                Issue(
                    Severity.ERROR,
                    "give all four of the shear, torsion and rotary-inertia values or none",
                    self,
                    None,
                    2,
                )
            )
        return found


# --------------------------------------------------------------------------- line types


class LineType(ModelObject):
    """Base class of the line types (one deck ``LINE TYPES`` row each)."""

    type_label = "Line type"
    abstract = True
    diameter = Quantity(
        LENGTH, 0.1, limit="line_diameter", group="Geometry", doc="Volume-equivalent diameter."
    )
    mass_per_length = Quantity(
        MASS_PER_LENGTH, 100.0, limit="line_mass", group="Geometry", doc="Dry mass per metre."
    )
    axial = Strategy(AxialModel, LinearAxial, group="Axial")
    bending = Child(BendingModel, group="Bending")
    drag_normal = Quantity(
        DIMENSIONLESS, 1.2, limit="small_coefficient", group=_HYDRO, doc="Normal drag Cd_n."
    )
    drag_axial = Quantity(
        DIMENSIONLESS, 0.0, limit="small_coefficient", group=_HYDRO, doc="Axial drag Cd_t."
    )
    added_mass_normal = Quantity(
        DIMENSIONLESS, 1.0, limit="small_coefficient", group=_HYDRO, doc="Normal Ca_n."
    )
    added_mass_axial = Quantity(
        DIMENSIONLESS, 0.0, limit="small_coefficient", group=_HYDRO, doc="Axial Ca_t."
    )
    breaking_load = OptionalQuantity(
        FORCE, minimum=0.0, role=_C, group="Strength", doc="Minimum breaking load (MBL)."
    )
    colour = Colour("#d9c21a", group="Appearance", doc="Colour in the engineering view.")
    appearance = Strategy(LineAppearance, PlainLineAppearance, role=Role.APPEARANCE)


@model_type("line_type.generic")
class GenericLineType(LineType):
    """A line type given only by its physics (the form a deck import produces)."""

    type_label = "Line type"


@model_type("line_type.chain")
class ChainType(LineType):
    """Mooring chain."""

    type_label = "Chain"
    grade = Choice(("R3", "R3S", "R4", "R4S", "R5", "R6"), role=_C, group="Construction")
    link = Choice(("studless", "studlink"), role=_C, group="Construction")
    nominal_diameter = OptionalQuantity(LENGTH, role=_C, group="Construction", doc="Bar d.")
    appearance = Strategy(LineAppearance, ChainAppearance, role=Role.APPEARANCE)


@model_type("line_type.wire_rope")
class WireRopeType(LineType):
    """Steel wire rope."""

    type_label = "Wire rope"
    construction = Choice(
        ("six_strand", "spiral_strand", "sheathed_spiral_strand"), role=_C, group="Construction"
    )
    nominal_diameter = OptionalQuantity(LENGTH, role=_C, group="Construction")
    appearance = Strategy(LineAppearance, WireRopeAppearance, role=Role.APPEARANCE)


@model_type("line_type.fibre_rope")
class FibreRopeType(LineType):
    """Synthetic fibre rope."""

    type_label = "Fibre rope"
    material = Choice(("polyester", "nylon", "hmpe", "other"), role=_C, group="Construction")
    construction = Choice(
        ("eight_strand", "twelve_strand", "parallel_strand"), role=_C, group="Construction"
    )
    nominal_diameter = OptionalQuantity(LENGTH, role=_C, group="Construction")
    appearance = Strategy(LineAppearance, FibreRopeAppearance, role=Role.APPEARANCE)


@model_type("line_type.cable")
class CableType(LineType):
    """Power or umbilical cable."""

    type_label = "Cable"
    cores = Integer(3, minimum=0, role=_C, group="Construction")
    armour_layers = Integer(1, minimum=0, role=_C, group="Construction")
    storage_bend_radius = OptionalQuantity(
        LENGTH, minimum=0.0, role=_C, group="Strength", doc="Minimum bend radius in storage."
    )
    dynamic_bend_radius = OptionalQuantity(
        LENGTH, minimum=0.0, role=_C, group="Strength", doc="Minimum bend radius in service."
    )
    appearance = Strategy(LineAppearance, CableAppearance, role=Role.APPEARANCE)


# --------------------------------------------------------------------------- rod types


@model_type("rod_type")
class RodType(ModelObject):
    """A rod type (one deck ``ROD TYPES`` row)."""

    type_label = "Rod type"
    abstract = False
    diameter = Quantity(LENGTH, 1.0, limit="line_diameter", group="Geometry")
    mass_per_length = Quantity(MASS_PER_LENGTH, 100.0, limit="line_mass", group="Geometry")
    drag = Quantity(DIMENSIONLESS, 0.6, limit="small_coefficient", group=_HYDRO, doc="Cd.")
    added_mass = Quantity(DIMENSIONLESS, 1.0, limit="small_coefficient", group=_HYDRO, doc="Ca.")
    end_drag = Quantity(DIMENSIONLESS, 0.0, limit="small_coefficient", group=_HYDRO)
    end_added_mass = Quantity(DIMENSIONLESS, 0.0, limit="small_coefficient", group=_HYDRO)
    axial_drag = OptionalQuantity(DIMENSIONLESS, limit="small_coefficient", group=_HYDRO)
    axial_added_mass = OptionalQuantity(DIMENSIONLESS, limit="small_coefficient", group=_HYDRO)
    colour = Colour("#b0b0b0", group="Appearance")

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if (self.axial_drag is None) != (self.axial_added_mass is None):
            found.append(
                Issue(Severity.ERROR, "give both axial coefficients or neither", self, None, 2)
            )
        return found


# --------------------------------------------------------------------------- floater types


class HydroDatabase(ModelObject):
    """Base class of a floater type's hydrodynamic database.

    This is the extension point for a time-domain floater solver: a database
    kind (added mass, radiation damping and retardation functions, first-order
    excitation and RAOs, second-order drift, hydrostatics, wind and current
    load coefficients, imported from a radiation-diffraction tool) registers
    as a subclass with :func:`cabledyn.project.model_type`. Today floaters
    move by prescribed motion or as rigid bodies, and no database is written
    to the deck.
    """

    type_label = "Hydrodynamic database"
    abstract = True


@model_type("hydro.none")
class NoHydroDatabase(HydroDatabase):
    """No hydrodynamic database."""

    type_label = "None"


@model_type("hydro.file")
class HydroDatabaseFile(HydroDatabase):
    """A reference to a hydrodynamic database file produced by another tool."""

    type_label = "Database file"
    file = FilePath(role=_C, group=_HYDRO)
    source_format = OptionalText(role=_C, group=_HYDRO, doc="Name of the producing tool or format.")
    water_depth = OptionalQuantity(LENGTH, role=_C, group=_HYDRO, doc="Depth of the analysis.")


@model_type("floater_type")
class FloaterType(ModelObject):
    """A shared description of a floater (hull), referenced by :class:`Floater` objects.

    The physical description and the hydrodynamic database are construction
    data: the deck holds each floater as a prescribed-motion or rigid body
    row, so none of this is written to the deck today.
    """

    type_label = "Floater type"
    abstract = False
    length = OptionalQuantity(LENGTH, minimum=0.0, role=_C, group="Dimensions")
    beam = OptionalQuantity(LENGTH, minimum=0.0, role=_C, group="Dimensions")
    draft = OptionalQuantity(LENGTH, minimum=0.0, role=_C, group="Dimensions")
    mass = OptionalQuantity(MASS, minimum=0.0, role=_C, group="Mass")
    centre_of_gravity = OptionalVec3(LENGTH, role=_C, group="Mass")
    inertia = OptionalVec3(INERTIA, role=_C, group="Mass", doc="Ixx, Iyy, Izz about the CG.")
    displaced_volume = OptionalQuantity(VOLUME, minimum=0.0, role=_C, group="Hydrostatics")
    waterplane_area = OptionalQuantity(AREA, minimum=0.0, role=_C, group="Hydrostatics")
    hydrodynamics = Strategy(HydroDatabase, NoHydroDatabase, role=_C, group=_HYDRO)
    mesh = Child(MeshAsset, role=Role.APPEARANCE)
