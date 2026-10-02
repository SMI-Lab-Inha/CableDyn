# SPDX-License-Identifier: Apache-2.0
"""GUI-only project data: groups, camera bookmarks, layouts and specifications.

None of this is written to a deck. The layouts, plot and report
specifications are opaque JSON placeholders owned by the GUI.
"""

from __future__ import annotations

from typing import Any

from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.descriptors import (
    Children,
    Choice,
    JsonValue,
    OptionalColour,
    Quantity,
    RefList,
    Role,
    Vec3,
)
from cabledyn.project.units import ANGLE, LENGTH, UnitSystem

__all__ = ["CameraBookmark", "Group", "StudioData"]

_M = Role.META


@model_type("studio.group")
class Group(ModelObject):
    """A user group of objects in the model browser."""

    type_label = "Group"
    abstract = False
    members = RefList(ModelObject, role=_M, group="Group")
    colour = OptionalColour(group="Group")


@model_type("studio.camera")
class CameraBookmark(ModelObject):
    """A saved camera of a 3D view."""

    type_label = "Camera"
    abstract = False
    view = Choice(("engineering", "presentation"), role=_M, group="Camera")
    target = Vec3(LENGTH, role=_M, group="Camera", doc="Point looked at.")
    distance = Quantity(LENGTH, 100.0, minimum=0.0, role=_M, group="Camera")
    azimuth = Quantity(ANGLE, 270.0, role=_M, group="Camera")
    elevation = Quantity(ANGLE, 20.0, minimum=-90.0, maximum=90.0, role=_M, group="Camera")
    roll = Quantity(ANGLE, 0.0, role=_M, group="Camera")
    projection = Choice(("perspective", "parallel"), role=_M, group="Camera")


@model_type("studio")
class StudioData(ModelObject):
    """Everything the GUI stores in a project besides the model itself."""

    type_label = "Studio data"
    abstract = False
    groups = Children(Group, role=_M)
    cameras = Children(CameraBookmark, role=_M)
    units = JsonValue({}, group="Preferences", doc="Display-unit choices by dimension.")
    layouts = JsonValue({}, group="Workspace", doc="Window layouts.")
    plot_specs = JsonValue([], group="Workspace", doc="Saved graph specifications.")
    report_specs = JsonValue([], group="Workspace", doc="Saved report specifications.")

    def unit_system(self) -> UnitSystem:
        """Return the project's display units."""
        data: dict[str, Any] = self.units or {}
        return UnitSystem.from_dict(data)

    def set_unit_system(self, system: UnitSystem) -> None:
        """Store display-unit choices (not an undoable command)."""
        self.units = system.to_dict()
