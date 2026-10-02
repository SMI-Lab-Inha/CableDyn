# SPDX-License-Identifier: Apache-2.0
"""Ordered, observable collections of owned model objects."""

from __future__ import annotations

from collections.abc import Iterator
from typing import TYPE_CHECKING, Generic, TypeVar, overload

from cabledyn.project.events import ChildAdded, ChildMoved, ChildRemoved

if TYPE_CHECKING:
    from cabledyn.project.base import ModelObject

__all__ = ["ObjectCollection"]

ObjT = TypeVar("ObjT", bound="ModelObject")


class ObjectCollection(Generic[ObjT]):
    """An ordered list of objects owned by one model object.

    The collection is the value of a :class:`~cabledyn.project.Children`
    property. Inserting an object makes the owner its parent; removing it
    clears the parent. Every change publishes :class:`ChildAdded`,
    :class:`ChildRemoved`, or :class:`ChildMoved` on the owner. An object can
    belong to one collection at a time.

    Parameters
    ----------
    owner : ModelObject
        The owning object.
    name : str
        The property name of the collection on ``owner``.
    item_type : type
        The class every item must be an instance of.
    """

    def __init__(self, owner: ModelObject, name: str, item_type: type[ObjT]) -> None:
        self._owner = owner
        self._name = name
        self._type = item_type
        self._items: list[ObjT] = []

    @property
    def owner(self) -> ModelObject:
        """The owning object."""
        return self._owner

    @property
    def name(self) -> str:
        """The collection's property name on its owner."""
        return self._name

    @property
    def item_type(self) -> type[ObjT]:
        """The class of the items."""
        return self._type

    def __iter__(self) -> Iterator[ObjT]:
        return iter(tuple(self._items))

    def __len__(self) -> int:
        return len(self._items)

    def __bool__(self) -> bool:
        return bool(self._items)

    @overload
    def __getitem__(self, index: int) -> ObjT: ...

    @overload
    def __getitem__(self, index: slice) -> list[ObjT]: ...

    def __getitem__(self, index: int | slice) -> ObjT | list[ObjT]:
        return self._items[index]

    def __contains__(self, item: object) -> bool:
        return any(existing is item for existing in self._items)

    def index(self, item: ObjT) -> int:
        """Return the position of ``item`` (identity, not equality).

        Raises
        ------
        ValueError
            If ``item`` is not in the collection.
        """
        for position, existing in enumerate(self._items):
            if existing is item:
                return position
        raise ValueError(f"{item!r} is not in {self._name}")

    def find(self, name: str) -> ObjT | None:
        """Return the first item called ``name`` (case-insensitive), or ``None``."""
        wanted = name.casefold()
        for item in self._items:
            if item.name.casefold() == wanted:
                return item
        return None

    def by_uid(self, uid: str) -> ObjT | None:
        """Return the item with ``uid``, or ``None``."""
        for item in self._items:
            if item.uid == uid:
                return item
        return None

    def insert(self, index: int, item: ObjT) -> None:
        """Insert ``item`` before ``index`` (clamped to the collection).

        Raises
        ------
        TypeError
            If ``item`` has the wrong type.
        ValueError
            If ``item`` already has a parent, would own its own owner, or
            brings a uid the project already uses.
        """
        if not isinstance(item, self._type):
            raise TypeError(
                f"{self._name} holds {self._type.__name__} objects, got {type(item).__name__}"
            )
        if item.parent is not None:
            raise ValueError(f"{item.label()} already belongs to {item.parent.label()}")
        self._owner._check_adoption(item)
        self._owner.root()._adopt(item)
        position = max(0, min(index, len(self._items)))
        self._items.insert(position, item)
        item._attach(self._owner, self._name)
        self._owner.emit(ChildAdded(self._owner, self._name, item, position))

    def append(self, item: ObjT) -> ObjT:
        """Append ``item`` and return it."""
        self.insert(len(self._items), item)
        return item

    def extend(self, items: list[ObjT]) -> None:
        """Append several items in order."""
        for item in items:
            self.append(item)

    def remove(self, item: ObjT) -> int:
        """Remove ``item`` and return its former position.

        Raises
        ------
        ValueError
            If ``item`` is not in the collection.
        """
        position = self.index(item)
        self._owner.root()._release(item)
        del self._items[position]
        item._detach()
        self._owner.emit(ChildRemoved(self._owner, self._name, item, position))
        return position

    def move(self, item: ObjT, index: int) -> None:
        """Move ``item`` to position ``index`` (clamped)."""
        old = self.index(item)
        new = max(0, min(index, len(self._items) - 1))
        if new == old:
            return
        del self._items[old]
        self._items.insert(new, item)
        self._owner.emit(ChildMoved(self._owner, self._name, item, old, new))

    def clear(self) -> None:
        """Remove every item, last first."""
        for item in reversed(self._items[:]):
            self.remove(item)

    def __repr__(self) -> str:
        names = [item.name or item.uid[:8] for item in self._items]
        return f"<{self._name}: {names!r}>"
