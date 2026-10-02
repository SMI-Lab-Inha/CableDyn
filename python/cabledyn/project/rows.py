# SPDX-License-Identifier: Apache-2.0
"""Model-wide rows: equivalent buoyancy, Syrope history, failures, control,
external loads, and output channels."""

from __future__ import annotations

from cabledyn.project.appearance import ModuleGeometry
from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.bodies import Body
from cabledyn.project.descriptors import (
    Child,
    Children,
    Choice,
    Integer,
    OptionalQuantity,
    Quantity,
    Ref,
    RefList,
    Role,
    Text,
    Vector,
)
from cabledyn.project.lines import Line
from cabledyn.project.points import Point
from cabledyn.project.types import LineType
from cabledyn.project.units import (
    FORCE,
    FORCE_PER_LENGTH,
    LENGTH,
    LINEAR_DAMPING,
    QUADRATIC_DAMPING,
    TIME,
)

__all__ = [
    "Channel",
    "Control",
    "EquivalentBuoyancy",
    "ExternalLoad",
    "Failure",
    "NamedChannel",
    "ObjectChannel",
    "OutputRequest",
    "SyropeIC",
]


@model_type("equivalent_buoyancy")
class EquivalentBuoyancy(ModelObject):
    """A line type given by its net submerged weight (smeared buoyancy).

    One deck ``EQUIVALENT BUOYANCY`` row; it rewrites the referenced type.
    """

    type_label = "Equivalent buoyancy"
    abstract = False
    line_type = Ref(LineType, group="Section")
    diameter = Quantity(LENGTH, 0.5, limit="line_diameter", group="Section")
    submerged_weight = Quantity(
        FORCE_PER_LENGTH, 0.0, group="Section", doc="Net submerged weight (negative: uplift)."
    )
    module_geometry = Child(ModuleGeometry, role=Role.APPEARANCE)
    module_pitch = OptionalQuantity(
        LENGTH, minimum=0.0, role=Role.APPEARANCE, group="Appearance", doc="Cosmetic spacing."
    )


@model_type("syrope_ic")
class SyropeIC(ModelObject):
    """The prior load history of Syrope lines (one deck ``SYROPE IC`` row)."""

    type_label = "Syrope history"
    abstract = False
    lines = RefList(Line, min_items=1, group="History")
    max_tension = Quantity(FORCE, 0.0, minimum=0.0, group="History", doc="Running maximum.")
    mean_tension = Quantity(FORCE, 0.0, group="History")


@model_type("failure")
class Failure(ModelObject):
    """Lines detaching from a point at a time or tension (one deck ``FAILURE`` row)."""

    type_label = "Failure"
    abstract = False
    point = Ref(Point, group="Failure")
    lines = RefList(Line, min_items=1, group="Failure")
    time = Quantity(TIME, 0.0, minimum=0.0, group="Trigger", doc="0 disables the time trigger.")
    tension = Quantity(
        FORCE, 0.0, minimum=0.0, group="Trigger", doc="0 disables the tension trigger."
    )


@model_type("control")
class Control(ModelObject):
    """Lines whose length follows a coupling-host control channel (deck ``CONTROL``)."""

    type_label = "Control"
    abstract = False
    channel = Integer(1, minimum=1, group="Control")
    lines = RefList(Line, min_items=1, group="Control")


@model_type("external_load")
class ExternalLoad(ModelObject):
    """A constant force and damping on a rigid body (deck ``EXTERNAL LOADS``)."""

    type_label = "External load"
    abstract = False
    body = Ref(Body, group="Load")
    axes = Choice(("global", "body"), group="Load", doc="Global axes or body axes.")
    force = Vector(FORCE, (0.0,), lengths=(1, 3), group="Load")
    linear_damping = Vector(LINEAR_DAMPING, (0.0,), lengths=(1, 3), group="Damping")
    quadratic_damping = Vector(QUADRATIC_DAMPING, (0.0,), lengths=(1, 3), group="Damping")


class Channel(ModelObject):
    """Base class of a requested output channel."""

    type_label = "Channel"
    abstract = True

    def channel_name(self, deck_id: int | None) -> str:
        """Return the deck channel name, given the target's deck id (if any)."""
        raise NotImplementedError


@model_type("channel.object")
class ObjectChannel(Channel):
    """A channel of one object, such as a line's fairlead tension.

    The deck name is ``quantity`` + the target's deck id + ``qualifier``
    (``FairTen`` + ``1``; ``Ten`` + ``2`` + ``N5``). The channel follows its
    target when the target is renumbered and is removed with it.
    """

    type_label = "Object channel"
    target = Ref(ModelObject, group="Channel")
    quantity = Text("", group="Channel", doc="Channel stem, for example FairTen.")
    qualifier = Text("", group="Channel", doc="Suffix after the id, for example N5 or pz.")

    def channel_name(self, deck_id: int | None) -> str:
        return f"{self.quantity}{'' if deck_id is None else deck_id}{self.qualifier}"


@model_type("channel.named")
class NamedChannel(Channel):
    """A channel given by its literal deck name."""

    type_label = "Named channel"
    channel = Text("", group="Channel")

    def channel_name(self, deck_id: int | None) -> str:
        return self.channel


@model_type("outputs")
class OutputRequest(ModelObject):
    """The requested output channels (the deck's ``OUTPUTS`` list)."""

    type_label = "Outputs"
    abstract = False
    channels = Children(Channel, group="Channels")
