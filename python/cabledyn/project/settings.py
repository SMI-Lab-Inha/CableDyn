# SPDX-License-Identifier: Apache-2.0
"""Analysis settings: the solver ``OPTIONS`` of statics, dynamics and modal runs.

Every option property is optional: ``None`` writes no row and leaves the
solver's own default in force. Options the model does not type are kept
verbatim as :class:`ExtraOption` rows, so nothing in a deck is lost.
"""

from __future__ import annotations

from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.descriptors import (
    Children,
    OptionalBool,
    OptionalInteger,
    OptionalQuantity,
    OptionalText,
    Text,
    TextList,
)
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.units import DIMENSIONLESS, TIME

__all__ = ["AnalysisSettings", "ExtraOption"]

_RUN = "Run"
_STATICS = "Statics"
_DYNAMICS = "Dynamics"
_SOLVER = "Dynamic solver"
_BODIES = "Bodies and rods"
_COMPAT = "Compatibility"


@model_type("option.extra")
class ExtraOption(ModelObject):
    """An ``OPTIONS`` row kept verbatim (keyword, value tokens, description)."""

    type_label = "Option"
    abstract = False
    keyword = Text("", group="Option")
    values = TextList((), group="Option", doc="Value tokens as written.")
    note = OptionalText(group="Option", doc="Description written after ' - '.")

    def invariants(self) -> list[Issue]:
        """Warn when the row sets an option the model types.

        Such a row (a duplicate kept from a deck, or added by hand) never
        overrides the typed property: the deck writer places it before the
        typed row, so the typed value is the one the solver uses.
        """
        from cabledyn.project.deck import typed_option

        found = super().invariants()
        typed = typed_option(self.keyword)
        if typed is not None:
            found.append(
                Issue(
                    Severity.WARNING,
                    f"this row sets the typed option {typed!r}; the typed property wins, so "
                    "edit it there and remove this row",
                    self,
                    "keyword",
                    2,
                )
            )
        return found


@model_type("settings")
class AnalysisSettings(ModelObject):
    """Run, statics, dynamics, modal and compatibility options."""

    type_label = "Analysis settings"
    abstract = False
    time_step = OptionalQuantity(TIME, minimum=0.0, group=_RUN, doc="Solver time step (dtM).")
    duration = OptionalQuantity(TIME, minimum=0.0, group=_RUN, doc="Run duration (TMax).")
    range_start = OptionalQuantity(TIME, minimum=0.0, group=_RUN, doc="Range-graph start.")
    cable_statics = OptionalText(group=_STATICS, doc="Finite-EI static route order.")
    body_initial_condition = OptionalText(group=_STATICS, doc="Start of free bodies (bodyIC).")
    initial_condition_mode = OptionalText(group=_STATICS, doc="Initial-condition mode (ICmode).")
    modified_newton = OptionalBool(group=_STATICS, doc="Guarded tangent reuse.")
    adaptive_mesh = OptionalBool(group=_STATICS, doc="Automatic finite-EI refinement.")
    spectral_radius = OptionalQuantity(
        DIMENSIONLESS, minimum=0.0, maximum=1.0, group=_DYNAMICS, doc="rhoInf."
    )
    max_strain = OptionalQuantity(DIMENSIONLESS, minimum=0.0, group=_DYNAMICS)
    alpha_force_blend = OptionalBool(group=_DYNAMICS)
    tensile_safety = OptionalText(group=_DYNAMICS, doc="Element axial-force audit mode.")
    tensile_strain_tolerance = OptionalQuantity(DIMENSIONLESS, minimum=0.0, group=_DYNAMICS)
    recovery_max_substeps = OptionalInteger(group=_DYNAMICS)
    axial_quadrature_order = OptionalInteger(group=_DYNAMICS)
    bending_quadrature_order = OptionalInteger(group=_DYNAMICS)
    solver_relative_tolerance = OptionalQuantity(DIMENSIONLESS, group=_SOLVER)
    solver_absolute_tolerance = OptionalQuantity(DIMENSIONLESS, group=_SOLVER)
    solver_max_iterations = OptionalInteger(group=_SOLVER)
    solver_backtracks = OptionalInteger(group=_SOLVER)
    solver_spectral_radius = OptionalQuantity(
        DIMENSIONLESS, group=_SOLVER, doc="Optional sixth value of the solver row."
    )
    cable_load_feedback = OptionalBool(group=_DYNAMICS, doc="Return cable reactions to the host.")
    body_wetting = OptionalText(group=_BODIES)
    body_hydrodynamics = OptionalText(group=_BODIES)
    rod_hydrodynamics = OptionalText(group=_BODIES)
    body_scheme = OptionalText(group=_BODIES)
    body_substep = OptionalText(group=_BODIES)
    mode_count = OptionalInteger(group="Modal", doc="Lowest natural modes per line.")
    time_scheme = OptionalText(group=_COMPAT, doc="Accepted for compatibility.")
    ic_time_step = OptionalQuantity(TIME, group=_COMPAT, doc="Accepted for compatibility.")
    ic_duration = OptionalQuantity(TIME, group=_COMPAT, doc="Accepted for compatibility.")
    ic_drag_scale = OptionalQuantity(DIMENSIONLESS, group=_COMPAT)
    ic_threshold = OptionalQuantity(DIMENSIONLESS, group=_COMPAT)
    write_log = OptionalInteger(group=_COMPAT)
    output_interval = OptionalQuantity(TIME, group=_COMPAT)
    extra_options = Children(ExtraOption, group="Other options")

    _SOLVER_FIELDS = (
        "solver_relative_tolerance",
        "solver_absolute_tolerance",
        "solver_max_iterations",
        "solver_backtracks",
    )

    def solver_row_given(self) -> bool:
        """Whether the dynamic-solver values are given (the four required ones)."""
        return self.solver_relative_tolerance is not None

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        given = [self.get_value(name) is not None for name in self._SOLVER_FIELDS]
        if any(given) and not all(given):
            found.append(
                Issue(
                    Severity.ERROR,
                    "give the four dynamic-solver values together",
                    self,
                    "solver_relative_tolerance",
                    2,
                )
            )
        if self.solver_spectral_radius is not None and not all(given):
            found.append(
                Issue(
                    Severity.ERROR,
                    "the solver spectral radius needs the four solver values",
                    self,
                    "solver_spectral_radius",
                    2,
                )
            )
        return found
