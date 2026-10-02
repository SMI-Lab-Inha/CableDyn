# SPDX-License-Identifier: Apache-2.0
"""Change events published by model objects.

Every :class:`~cabledyn.project.ModelObject` owns an :class:`EventHub`. A
change is published on the changed object's hub and then on the hub of every
ancestor, so a subscriber on the :class:`~cabledyn.project.Project` sees every
change in the model. The events are plain, immutable records; a GUI bridge
re-emits them as toolkit signals.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from cabledyn.project.base import ModelObject

__all__ = [
    "ChildAdded",
    "ChildMoved",
    "ChildRemoved",
    "Event",
    "EventHub",
    "PropertyChanged",
    "ReferenceChanged",
    "Subscription",
]


@dataclass(frozen=True)
class Event:
    """Base class of every model event.

    Attributes
    ----------
    source : ModelObject
        The object whose state changed.
    """

    source: ModelObject


@dataclass(frozen=True)
class PropertyChanged(Event):
    """A value property (number, text, choice, child object) changed.

    Attributes
    ----------
    name : str
        Property name.
    old, new : object
        Previous and new value.
    """

    name: str
    old: Any
    new: Any


@dataclass(frozen=True)
class ReferenceChanged(Event):
    """A reference property (:class:`~cabledyn.project.Ref` or
    :class:`~cabledyn.project.RefList`) now points elsewhere.

    Attributes
    ----------
    name : str
        Property name.
    old, new : object
        Previous and new target (an object, a tuple of objects, or ``None``).
    """

    name: str
    old: Any
    new: Any


@dataclass(frozen=True)
class ChildAdded(Event):
    """An object was inserted into a child collection of ``source``.

    Attributes
    ----------
    collection : str
        Name of the :class:`~cabledyn.project.Children` property.
    child : ModelObject
        The inserted object.
    index : int
        Its position.
    """

    collection: str
    child: ModelObject
    index: int


@dataclass(frozen=True)
class ChildRemoved(Event):
    """An object was removed from a child collection of ``source``.

    Attributes
    ----------
    collection : str
        Name of the collection.
    child : ModelObject
        The removed object.
    index : int
        Its former position.
    """

    collection: str
    child: ModelObject
    index: int


@dataclass(frozen=True)
class ChildMoved(Event):
    """An object moved within a child collection of ``source``.

    Attributes
    ----------
    collection : str
        Name of the collection.
    child : ModelObject
        The moved object.
    old_index, new_index : int
        Former and new position.
    """

    collection: str
    child: ModelObject
    old_index: int
    new_index: int


Listener = Callable[[Event], None]


class Subscription:
    """Handle returned by :meth:`EventHub.subscribe`; call :meth:`cancel` to stop."""

    def __init__(self, hub: EventHub, listener: Listener) -> None:
        self._hub = hub
        self._listener = listener

    def cancel(self) -> None:
        """Stop delivering events to the listener (idempotent)."""
        self._hub._remove(self._listener)


class EventHub:
    """A list of listeners called synchronously, in subscription order."""

    def __init__(self) -> None:
        self._listeners: list[Listener] = []

    def subscribe(self, listener: Listener) -> Subscription:
        """Call ``listener(event)`` for every event published on this hub.

        Parameters
        ----------
        listener : Callable[[Event], None]
            The callback.

        Returns
        -------
        Subscription
            Handle to cancel the subscription.
        """
        self._listeners.append(listener)
        return Subscription(self, listener)

    def _remove(self, listener: Listener) -> None:
        for index, existing in enumerate(self._listeners):
            if existing is listener:
                del self._listeners[index]
                return

    def emit(self, event: Event) -> None:
        """Deliver ``event`` to every listener.

        Parameters
        ----------
        event : Event
            The event.
        """
        for listener in tuple(self._listeners):
            listener(event)

    def __len__(self) -> int:
        return len(self._listeners)
