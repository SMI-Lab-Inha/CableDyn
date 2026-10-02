# SPDX-License-Identifier: Apache-2.0
"""Reverse references ("used by") and the cascade rules of object removal.

Objects refer to each other through :class:`~cabledyn.project.Ref` and
:class:`~cabledyn.project.RefList` properties. Removing an object that is
still referred to either fails (:class:`ObjectInUseError`) or cascades:

* a referrer whose *required* reference loses its target is removed too;
* an optional reference is cleared;
* a reference list drops the target, and a referrer left with fewer targets
  than its ``min_items`` is removed;
* an owner whose child collection falls below its ``min_items`` (a line
  without sections) is removed;

and so on until nothing changes. Only objects held in a collection can be
removed; an object owned through a single :class:`~cabledyn.project.Child`
slot (a rod end, the environment) goes with its owner.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from cabledyn.project.base import ModelObject, collection_of
from cabledyn.project.descriptors import Children, Ref, RefList

__all__ = ["ObjectInUseError", "RemovalPlan", "plan_removal", "referrers"]


class ObjectInUseError(ValueError):
    """An object cannot be removed without cascading: other objects refer to it.

    Attributes
    ----------
    referrers : tuple[tuple[ModelObject, str], ...]
        The referring objects and property names.
    """

    def __init__(self, message: str, found: list[tuple[ModelObject, str]]) -> None:
        super().__init__(message)
        self.referrers = tuple(found)


def _ids(obj: ModelObject) -> set[int]:
    return {id(item) for item in obj.walk()}


def referrers(root: ModelObject, obj: ModelObject) -> list[tuple[ModelObject, str]]:
    """Return ``(referrer, property name)`` pairs that point into ``obj``.

    Parameters
    ----------
    root : ModelObject
        The object graph to search (normally the project).
    obj : ModelObject
        The referenced object; references to its descendants count too.

    Returns
    -------
    list[tuple[ModelObject, str]]
        Referrers outside ``obj``, in model order, each property once.
    """
    inside = _ids(obj)
    found: list[tuple[ModelObject, str]] = []
    for item in root.walk():
        if id(item) in inside:
            continue
        names: list[str] = []
        for name, target in item.references():
            if id(target) in inside and name not in names:
                names.append(name)
        found.extend((item, name) for name in names)
    return found


@dataclass
class RemovalPlan:
    """What removing objects does to the model.

    Attributes
    ----------
    removed : list[ModelObject]
        Objects to take out of their collections (none inside another).
    reference_updates : list[tuple[ModelObject, str, Any, Any]]
        ``(object, property, old value, new value)`` for surviving referrers.
    """

    removed: list[ModelObject] = field(default_factory=list)
    reference_updates: list[tuple[ModelObject, str, Any, Any]] = field(default_factory=list)


def _removable(obj: ModelObject) -> bool:
    return collection_of(obj) is not None


def plan_removal(
    root: ModelObject, targets: list[ModelObject], *, cascade: bool = False
) -> RemovalPlan:
    """Work out the cascade of removing ``targets`` from ``root``.

    Parameters
    ----------
    root : ModelObject
        The object graph (normally the project).
    targets : list[ModelObject]
        Objects to remove; each must sit in a collection.
    cascade : bool
        Allow dependants to be removed or updated; otherwise refuse when
        anything outside the targets refers to them.

    Returns
    -------
    RemovalPlan
        The objects to remove and the references to update.

    Raises
    ------
    ValueError
        If a target is not held in a collection.
    ObjectInUseError
        If ``cascade`` is false and a target is still referred to.
    """
    for target in targets:
        if not _removable(target):
            raise ValueError(f"{target.label()} is not in a collection and cannot be removed")
    if not cascade:
        found: list[tuple[ModelObject, str]] = []
        doomed_ids: set[int] = set()
        for target in targets:
            doomed_ids |= _ids(target)
        for target in targets:
            for item, name in referrers(root, target):
                if id(item) not in doomed_ids and (item, name) not in found:
                    found.append((item, name))
        if found:
            names = ", ".join(f"{item.label()}.{name}" for item, name in found)
            raise ObjectInUseError(
                f"{targets[0].label()} is still used by {names}; remove those first or cascade",
                found,
            )
    doomed: list[ModelObject] = []
    for target in targets:
        if all(item is not target for item in doomed):
            doomed.append(target)
    changed = True
    while changed:
        changed = False
        doomed_ids = set()
        for item in doomed:
            doomed_ids |= _ids(item)
        candidates: list[ModelObject] = []
        for item in doomed:
            owner = item.parent
            collection = collection_of(item)
            if owner is None or collection is None or id(owner) in doomed_ids:
                continue
            prop = type(owner).get_property(collection.name)
            if isinstance(prop, Children) and prop.min_items > 0:
                remaining = [x for x in collection if id(x) not in doomed_ids]
                if len(remaining) < prop.min_items and _removable(owner):
                    candidates.append(owner)
        for item in root.walk():
            if id(item) in doomed_ids or not _removable(item):
                continue
            for prop in type(item).properties():
                value = item._values[prop.name]
                if isinstance(prop, Ref) and prop.required:
                    if value is not None and id(value) in doomed_ids:
                        candidates.append(item)
                        break
                elif isinstance(prop, RefList) and prop.min_items > 0:
                    remaining = [x for x in value if id(x) not in doomed_ids]
                    if len(remaining) != len(value) and len(remaining) < prop.min_items:
                        candidates.append(item)
                        break
        for item in candidates:
            if id(item) not in doomed_ids and all(x is not item for x in doomed):
                doomed.append(item)
                changed = True
    doomed_ids = set()
    for item in doomed:
        doomed_ids |= _ids(item)
    plan = RemovalPlan()
    for item in root.walk():
        if id(item) in doomed_ids:
            continue
        for prop in type(item).properties():
            value = item._values[prop.name]
            if isinstance(prop, Ref) and value is not None and id(value) in doomed_ids:
                plan.reference_updates.append((item, prop.name, value, None))
            elif isinstance(prop, RefList):
                kept = tuple(x for x in value if id(x) not in doomed_ids)
                if len(kept) != len(value):
                    plan.reference_updates.append((item, prop.name, value, kept))
    for item in doomed:
        ancestor = item.parent
        nested = False
        while ancestor is not None:
            if any(ancestor is other for other in doomed):
                nested = True
                break
            ancestor = ancestor.parent
        if not nested:
            plan.removed.append(item)
    return plan
