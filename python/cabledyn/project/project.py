# SPDX-License-Identifier: Apache-2.0
"""The :class:`Project`: the root of the object graph."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Any, TypeVar

from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.bodies import Body, Turbine
from cabledyn.project.descriptors import (
    Bool,
    Child,
    Children,
    OptionalText,
    Role,
    Strategy,
    Text,
)
from cabledyn.project.environment import (
    Environment,
    FlatSeabed,
    MotionSource,
    NoMotion,
    Seabed,
)
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.lines import Line
from cabledyn.project.points import Point
from cabledyn.project.rods import Rod
from cabledyn.project.rows import (
    Control,
    EquivalentBuoyancy,
    ExternalLoad,
    Failure,
    OutputRequest,
    SyropeIC,
)
from cabledyn.project.settings import AnalysisSettings
from cabledyn.project.studio import StudioData
from cabledyn.project.types import FloaterType, LineType, RodType

__all__ = ["Project"]

M = TypeVar("M", bound=ModelObject)


@model_type("project")
class Project(ModelObject):
    """A CableDyn model: environment, settings, libraries, objects and GUI data.

    Examples
    --------
    >>> from cabledyn.project import (
    ...     FixedPoint, FloaterPoint, GenericLineType, Line, Project, Section)
    >>> project = Project("Single chain")
    >>> chain = project.line_types.append(GenericLineType("chain", diameter=0.252))
    >>> fairlead = project.points.append(FloaterPoint("fairlead", position=(0, 0, -14)))
    >>> anchor = project.points.append(FixedPoint("anchor", position=(400, 0, -50)))
    >>> line = project.lines.append(Line("mooring", end_a=fairlead, end_b=anchor))
    >>> _ = line.sections.append(Section(line_type=chain, length=410.0, segments=41))
    >>> line.length
    410.0
    """

    type_label = "Project"
    abstract = False
    title = Text("CableDyn model", role=Role.META, group="Project", doc="Deck title.")
    caller_driven = Bool(
        False,
        role=Role.META,
        group="Project",
        doc="Validate the deck for the OpenFAST coupling instead of the standalone driver.",
    )
    deck_path = OptionalText(
        role=Role.META,
        group="Project",
        doc="Nominal deck path; relative files (motion, bathymetry, ...) resolve from its folder.",
    )
    environment = Child(Environment)
    seabed = Strategy(Seabed, FlatSeabed)
    motion = Strategy(MotionSource, NoMotion, doc="Prescribed motion of a standalone run.")
    settings = Child(AnalysisSettings)
    line_types = Children(LineType, group="Libraries")
    rod_types = Children(RodType, group="Libraries")
    floater_types = Children(FloaterType, group="Libraries")
    bodies = Children(Body, group="Objects")
    rods = Children(Rod, group="Objects")
    turbines = Children(Turbine, group="Objects")
    points = Children(Point, group="Objects")
    lines = Children(Line, group="Objects")
    equivalent_buoyancy = Children(EquivalentBuoyancy, group="Line rows")
    syrope_history = Children(SyropeIC, group="Line rows")
    failures = Children(Failure, group="Events")
    controls = Children(Control, group="Events")
    external_loads = Children(ExternalLoad, group="Loads")
    outputs = Child(OutputRequest)
    studio = Child(StudioData, role=Role.META)

    def _setup(self) -> None:
        super()._setup()
        self.unknown_objects: list[dict[str, Any]] = []
        self.load_warnings: list[str] = []

    def find(self, uid: str) -> ModelObject | None:
        """Return the object with ``uid`` anywhere in the project, or ``None``."""
        for obj in self.walk():
            if obj.uid == uid:
                return obj
        return None

    def objects_of(self, kind: type[M]) -> Iterator[M]:
        """Yield every object of class ``kind`` in the project, depth first."""
        for obj in self.walk():
            if isinstance(obj, kind):
                yield obj

    def used_by(self, obj: ModelObject) -> list[tuple[ModelObject, str]]:
        """Return the objects (and property names) that refer to ``obj`` or its descendants.

        Referrers inside ``obj`` itself are not listed.
        """
        from cabledyn.project.references import referrers

        return referrers(self, obj)

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        found.extend(_unique_names(self.line_types, "line type"))
        found.extend(_unique_names(self.rod_types, "rod type"))
        numbers: dict[int, Turbine] = {}
        for turbine in self.turbines:
            if turbine.number in numbers:
                found.append(
                    Issue(Severity.ERROR, "turbine number is used twice", turbine, "number", 2)
                )
            numbers[turbine.number] = turbine
        if self.unknown_objects:
            found.append(
                Issue(
                    Severity.WARNING,
                    f"{len(self.unknown_objects)} object(s) of unknown type are kept but "
                    "not written to the deck",
                    self,
                    None,
                    2,
                )
            )
        return found

    def validate_all(self, *, deck: bool = True) -> list[Issue]:
        """Run every validation layer over the whole project.

        Parameters
        ----------
        deck : bool
            Also run layer 3 (write the deck and apply the native deck rules).

        Returns
        -------
        list[Issue]
            Every issue found.
        """
        from cabledyn.project.validation import validate_project

        return validate_project(self, deck=deck)


def _unique_names(items: Iterable[ModelObject], what: str) -> list[Issue]:
    found: list[Issue] = []
    seen: set[str] = set()
    for item in items:
        key = item.name.casefold()
        if not item.name:
            found.append(Issue(Severity.ERROR, f"a {what} needs a name", item, "name", 2))
        elif key in seen:
            found.append(
                Issue(Severity.ERROR, f"{what} name {item.name!r} is used twice", item, "name", 2)
            )
        seen.add(key)
    return found
