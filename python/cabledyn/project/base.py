# SPDX-License-Identifier: Apache-2.0
"""The :class:`ModelObject` base class, the type registry, and serialisation.

Every object of a CableDyn project is a :class:`ModelObject`: it has a stable
``uid``, a ``name``, a ``parent``, typed properties declared with the
descriptors of :mod:`cabledyn.project.descriptors`, change events, validation,
and a JSON form. Concrete classes register under a stable type key with
:func:`model_type`, which the project file uses to rebuild them.
"""

from __future__ import annotations

import uuid
from collections.abc import Callable, Iterator
from typing import Any, ClassVar, TypeVar

from cabledyn.project.collection import ObjectCollection
from cabledyn.project.descriptors import (
    Child,
    Children,
    ModelObjectLike,
    Property,
    Ref,
    RefList,
    Role,
    Text,
    TextList,
)
from cabledyn.project.events import Event, EventHub, PropertyChanged, ReferenceChanged
from cabledyn.project.issues import Issue

__all__ = [
    "REGISTRY",
    "ModelObject",
    "SerialContext",
    "TypeRegistry",
    "UnknownTypeError",
    "model_type",
]

M = TypeVar("M", bound="ModelObject")


class UnknownTypeError(KeyError):
    """A type key is not registered."""


class TypeRegistry:
    """Maps stable type keys to model classes (and back).

    The project file stores each object's key; a plugin registers new object
    kinds here (see :func:`model_type`).
    """

    def __init__(self) -> None:
        self._by_key: dict[str, type[ModelObject]] = {}

    def register(self, key: str, cls: type[ModelObject]) -> None:
        """Register ``cls`` under ``key``.

        Raises
        ------
        ValueError
            If ``key`` is taken by another class.
        """
        existing = self._by_key.get(key)
        if existing is not None and existing is not cls:
            raise ValueError(f"type key {key!r} is already used by {existing.__name__}")
        self._by_key[key] = cls

    def get(self, key: str) -> type[ModelObject]:
        """Return the class registered under ``key``.

        Raises
        ------
        UnknownTypeError
            If no class is registered under ``key``.
        """
        try:
            return self._by_key[key]
        except KeyError as exc:
            raise UnknownTypeError(key) from exc

    def __contains__(self, key: object) -> bool:
        return key in self._by_key

    def keys(self) -> list[str]:
        """Return the registered keys, sorted."""
        return sorted(self._by_key)

    def subclasses_of(self, base: type[ModelObject]) -> list[type[ModelObject]]:
        """Return the registered concrete classes that are ``base`` or derive from it."""
        return [
            cls
            for cls in self._by_key.values()
            if issubclass(cls, base) and not cls.__dict__.get("abstract", False)
        ]


REGISTRY = TypeRegistry()
"""The global registry of model classes."""


def model_type(key: str) -> Callable[[type[M]], type[M]]:
    """Class decorator registering a model class under a stable type key.

    Parameters
    ----------
    key : str
        Stable key written to project files, for example ``"line"``.

    Returns
    -------
    Callable
        The decorator; it returns the class unchanged.

    Examples
    --------
    >>> from cabledyn.project import ModelObject, Quantity, model_type
    >>> @model_type("example.marker")
    ... class Marker(ModelObject):
    ...     height = Quantity(default=1.0)
    >>> Marker("m1").height
    1.0
    """

    def register(cls: type[M]) -> type[M]:
        cls.type_key = key
        REGISTRY.register(key, cls)
        return cls

    return register


class SerialContext:
    """State of one project-file read or write.

    It maps uids to rebuilt objects, defers references until every object
    exists, and collects objects of unknown type (kept verbatim) and warnings.

    Parameters
    ----------
    fresh_uids : bool
        Give rebuilt objects new uids (for copy and paste) instead of the
        stored ones.
    """

    def __init__(self, *, fresh_uids: bool = False) -> None:
        self.fresh_uids = fresh_uids
        self.objects: dict[str, ModelObject] = {}
        self.unknown: list[dict[str, Any]] = []
        self.warnings: list[str] = []
        self._pending: list[tuple[ModelObject, Ref[Any] | RefList[Any], Any]] = []

    def warn(self, message: str) -> None:
        """Record a warning."""
        self.warnings.append(message)

    def build(
        self,
        data: Any,
        *,
        owner: ModelObject | None = None,
        slot: str | None = None,
        index: int | None = None,
    ) -> ModelObject | None:
        """Rebuild one object from its dictionary, or keep it verbatim if unknown.

        Returns
        -------
        ModelObject | None
            The object, or ``None`` for an unknown type (recorded in
            ``unknown``).
        """
        if not isinstance(data, dict) or not isinstance(data.get("type"), str):
            self.warn(f"skipped a malformed object entry: {data!r:.80}")
            return None
        key = data["type"]
        if key not in REGISTRY:
            self.unknown.append(
                {
                    "owner": None if owner is None else owner.uid,
                    "slot": slot,
                    "index": index,
                    "data": data,
                }
            )
            self.warn(f"kept an object of unknown type {key!r} verbatim")
            return None
        cls = REGISTRY.get(key)
        obj = cls.__new__(cls)
        obj._setup()
        stored = data.get("uid")
        if isinstance(stored, str) and stored and not self.fresh_uids:
            obj._uid = stored
        if isinstance(stored, str) and stored:
            self.objects[stored] = obj
        hints = data.get("deck_hints")
        if isinstance(hints, dict):
            obj.deck_hints = dict(hints)
        values = data.get("properties", {})
        if not isinstance(values, dict):
            values = {}
        for prop in cls.properties():
            if prop.name not in values:
                continue
            try:
                value = prop.from_json(values[prop.name], self, obj)
            except (TypeError, ValueError) as exc:
                self.warn(f"{cls.__name__}.{prop.name}: kept the default ({exc})")
                continue
            if not isinstance(prop, Children):
                old = obj._values[prop.name]
                if isinstance(prop, Child) and old is not value:
                    old._detach()
                obj._values[prop.name] = value
        for name in values:
            if name not in cls._property_map():
                self.warn(f"{cls.__name__}: ignored unknown property {name!r}")
        return obj

    def defer(self, obj: ModelObject, prop: Ref[Any] | RefList[Any], data: Any) -> None:
        """Resolve a reference once every object exists."""
        self._pending.append((obj, prop, data))

    def lookup(self, uid: Any) -> ModelObject:
        """Return the rebuilt object with ``uid``.

        Raises
        ------
        KeyError
            If no rebuilt object has ``uid``.
        """
        if not isinstance(uid, str) or uid not in self.objects:
            raise KeyError(f"no object with uid {uid!r}")
        return self.objects[uid]

    def resolve(self) -> None:
        """Resolve every deferred reference; dangling ones become warnings."""
        for obj, prop, data in self._pending:
            try:
                obj._values[prop.name] = prop.resolve(data, self)
            except (KeyError, TypeError) as exc:
                self.warn(f"{obj.label()}: dropped reference {prop.name} ({exc})")
        self._pending.clear()


class ModelObject(ModelObjectLike):
    """Base class of every object in a CableDyn project.

    Parameters
    ----------
    name : str
        Object name, shown in the browser. Names need not be unique, except
        where the deck needs it (line and rod type names).
    **values : object
        Initial property values, by property name.

    Raises
    ------
    AttributeError
        If a keyword is not a settable property of the class.

    Notes
    -----
    Assigning a property (``line.outputs = "pt"``) validates its type and
    publishes a change event, but is not undoable. Interactive edits go
    through :mod:`cabledyn.project.commands` and a
    :class:`~cabledyn.project.commands.CommandStack`.
    """

    type_key: ClassVar[str] = "object"
    """Stable type key in project files (set by :func:`model_type`)."""
    type_label: ClassVar[str] = "Object"
    """Human-readable class name for the GUI."""
    abstract: ClassVar[bool] = True
    """Abstract classes are not offered by :meth:`Strategy.options`."""

    name = Text("", role=Role.META, group="Identity", doc="Name shown in the model browser.")
    description = Text("", multiline=True, role=Role.META, group="Identity", doc="Free-text notes.")
    tags = TextList((), role=Role.META, group="Identity", doc="User tags for grouping.")

    _properties_cache: ClassVar[dict[type, tuple[Property[Any], ...]]] = {}

    def __init__(self, name: str = "", /, **values: object) -> None:
        self._setup()
        if name:
            self._values["name"] = type(self).name.coerce(name)
        props = type(self)._property_map()
        for key, value in values.items():
            prop = props.get(key)
            if prop is None or isinstance(prop, Children):
                raise AttributeError(f"{type(self).__name__} has no settable property {key!r}")
            setattr(self, key, value)

    def _setup(self) -> None:
        self._uid = uuid.uuid4().hex
        self._parent: ModelObject | None = None
        self._slot: str | None = None
        self.events = EventHub()
        self.deck_hints: dict[str, Any] = {}
        self._values: dict[str, Any] = {}
        for prop in type(self).properties():
            self._values[prop.name] = prop.initial(self)

    # ------------------------------------------------------------------ introspection

    @classmethod
    def properties(cls) -> tuple[Property[Any], ...]:
        """Return the property descriptors, base-class properties first."""
        cached = ModelObject._properties_cache.get(cls)
        if cached is not None:
            return cached
        found: dict[str, Property[Any]] = {}
        for klass in reversed(cls.__mro__):
            for key, value in vars(klass).items():
                if isinstance(value, Property):
                    found[key] = value
        result = tuple(found.values())
        ModelObject._properties_cache[cls] = result
        return result

    @classmethod
    def _property_map(cls) -> dict[str, Property[Any]]:
        return {prop.name: prop for prop in cls.properties()}

    @classmethod
    def get_property(cls, name: str) -> Property[Any]:
        """Return the descriptor called ``name``.

        Raises
        ------
        KeyError
            If the class has no such property.
        """
        try:
            return cls._property_map()[name]
        except KeyError as exc:
            raise KeyError(f"{cls.__name__} has no property {name!r}") from exc

    @property
    def uid(self) -> str:
        """Stable unique identifier (32 hexadecimal characters)."""
        return self._uid

    @property
    def parent(self) -> ModelObject | None:
        """The owning object, or ``None`` for a root or a detached object."""
        return self._parent

    @property
    def slot(self) -> str | None:
        """Name of the parent's property that owns this object."""
        return self._slot

    def label(self) -> str:
        """Return a short description such as ``Line 'Mooring 1'``."""
        return f"{self.type_label} {self.name!r}" if self.name else self.type_label

    def root(self) -> ModelObject:
        """Return the top-most ancestor (the project, when attached)."""
        node: ModelObject = self
        while node._parent is not None:
            node = node._parent
        return node

    def get_value(self, name: str) -> Any:
        """Return the value of property ``name``."""
        type(self).get_property(name)
        return self._values[name]

    def set_value(self, name: str, value: object) -> None:
        """Set property ``name`` (as attribute assignment does)."""
        prop = type(self).get_property(name)
        prop.__set__(self, value)

    # ------------------------------------------------------------------ ownership

    def _attach(self, parent: ModelObject, slot: str) -> None:
        self._parent = parent
        self._slot = slot

    def _detach(self) -> None:
        self._parent = None
        self._slot = None

    def children(self) -> Iterator[ModelObject]:
        """Yield the directly owned sub-objects, in property order."""
        for prop in type(self).properties():
            if isinstance(prop, Children):
                yield from self._values[prop.name]
            elif isinstance(prop, Child):
                yield self._values[prop.name]

    def walk(self) -> Iterator[ModelObject]:
        """Yield this object and every owned descendant, depth first."""
        yield self
        for child in self.children():
            yield from child.walk()

    def references(self) -> list[tuple[str, ModelObject]]:
        """Return the outgoing references as ``(property name, target)`` pairs."""
        found: list[tuple[str, ModelObject]] = []
        for prop in type(self).properties():
            if isinstance(prop, Ref):
                target = self._values[prop.name]
                if target is not None:
                    found.append((prop.name, target))
            elif isinstance(prop, RefList):
                found.extend((prop.name, target) for target in self._values[prop.name])
        return found

    # ------------------------------------------------------------------ changes

    def _assign(self, prop: Property[Any], value: Any) -> None:
        if prop.read_only:
            raise AttributeError(f"{prop.name} is read-only")
        old = self._values[prop.name]
        if prop.same(old, value):
            return
        if isinstance(prop, Child):
            if value.parent is not None:
                raise ValueError(f"{value.label()} already belongs to {value.parent.label()}")
            old._detach()
            value._attach(self, prop.name)
        self._values[prop.name] = value
        event_type = ReferenceChanged if prop.is_reference else PropertyChanged
        self.emit(event_type(self, prop.name, old, value))

    def emit(self, event: Event) -> None:
        """Publish ``event`` on this object's hub and on every ancestor's hub."""
        node: ModelObject | None = self
        while node is not None:
            node.events.emit(event)
            node = node._parent

    # ------------------------------------------------------------------ validation

    def invariants(self) -> list[Issue]:
        """Return the layer-2 issues of this object (structural invariants).

        Subclasses extend this; the base checks references and collection
        sizes declared by the descriptors.
        """
        return []

    def validate(self, *, recursive: bool = False) -> list[Issue]:
        """Return the layer-1 and layer-2 issues of this object.

        Parameters
        ----------
        recursive : bool
            Also validate every owned descendant.

        Returns
        -------
        list[Issue]
            Property-check issues (layer 1), then reference and invariant
            issues (layer 2).
        """
        found: list[Issue] = []
        for prop in type(self).properties():
            found.extend(prop.check(self, self._values[prop.name]))
        found.extend(self.invariants())
        if recursive:
            for child in self.children():
                found.extend(child.validate(recursive=True))
        return found

    # ------------------------------------------------------------------ serialisation

    def to_dict(self, ctx: SerialContext | None = None) -> dict[str, Any]:
        """Return the JSON-compatible form of this object and its descendants.

        References are written as the target's uid.
        """
        context = ctx or SerialContext()
        data: dict[str, Any] = {
            "type": self.type_key,
            "uid": self._uid,
            "properties": {
                prop.name: prop.to_json(self._values[prop.name], context)
                for prop in type(self).properties()
            },
        }
        if self.deck_hints:
            data["deck_hints"] = _json_copy(self.deck_hints)
        return data

    @classmethod
    def from_dict(cls: type[M], data: dict[str, Any], *, fresh_uids: bool = False) -> M:
        """Rebuild an object (and its descendants) written by :meth:`to_dict`.

        References to objects outside ``data`` are dropped with a warning.

        Parameters
        ----------
        data : dict
            The dictionary.
        fresh_uids : bool
            Give the rebuilt objects new uids.

        Returns
        -------
        ModelObject
            The rebuilt object.

        Raises
        ------
        TypeError
            If ``data`` does not describe a ``cls`` object.
        """
        ctx = SerialContext(fresh_uids=fresh_uids)
        obj = ctx.build(data)
        if not isinstance(obj, cls):
            raise TypeError(f"data does not describe a {cls.__name__}")
        ctx.resolve()
        return obj

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.name!r} {self._uid[:8]}>"


def _json_copy(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_copy(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_copy(item) for item in value]
    return value


def collection_of(obj: ModelObject) -> ObjectCollection[Any] | None:
    """Return the collection that holds ``obj``, or ``None`` (a single child or a root)."""
    parent = obj.parent
    if parent is None or obj.slot is None:
        return None
    value = parent._values[obj.slot]
    return value if isinstance(value, ObjectCollection) else None
