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

    def unregister(self, key: str) -> None:
        """Remove the class registered under ``key`` (no-op if none)."""
        self._by_key.pop(key, None)

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
            cls for cls in self._by_key.values() if issubclass(cls, base) and not cls.is_abstract()
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
        REGISTRY.register(key, cls)
        cls.type_key = key
        return cls

    return register


class SerialContext:
    """State of one project-file read or write.

    It maps uids to rebuilt objects, defers references until every object
    exists, and collects objects of unknown type (kept verbatim), duplicate
    uids and warnings.

    Values this version cannot use (an unknown property, a value a property
    refuses, a reference to an object that is not rebuilt, a sub-object of
    unknown type) are kept raw on the object and written back unchanged by
    :meth:`ModelObject.to_dict` until the property is assigned.

    Parameters
    ----------
    fresh_uids : bool
        Give rebuilt objects new uids (for copy and paste) instead of the
        stored ones. References inside the rebuilt data still resolve.
    within : ModelObject | None
        An existing object graph (the project a copy is pasted into);
        references to its objects resolve to them.
    """

    def __init__(self, *, fresh_uids: bool = False, within: ModelObject | None = None) -> None:
        self.fresh_uids = fresh_uids
        self.within = within
        self.objects: dict[str, ModelObject] = {}
        self.unknown: list[dict[str, Any]] = []
        self.duplicates: list[str] = []
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
            self.warn(f"kept an object of unknown type {key!r} verbatim")
            if owner is not None and index is None:
                return None  # a sub-object slot: the owner keeps the raw data
            self.unknown.append(
                {
                    "owner": None if owner is None else owner.uid,
                    "slot": slot,
                    "index": index,
                    "data": data,
                }
            )
            return None
        cls = REGISTRY.get(key)
        obj = cls.__new__(cls)
        obj._setup()
        stored = data.get("uid")
        if isinstance(stored, str) and stored:
            if stored in self.objects:
                self.duplicates.append(stored)
                self.warn(f"uid {stored} is used by more than one object")
            else:
                self.objects[stored] = obj
                if not self.fresh_uids:
                    obj._uid = stored
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
                self.warn(
                    f"{cls.__name__}.{prop.name}: kept the default and the stored value "
                    f"for a later save ({exc})"
                )
                obj._preserved[prop.name] = values[prop.name]
                continue
            if not isinstance(prop, Children):
                old = obj._values[prop.name]
                if isinstance(prop, Child) and old is not value:
                    old._detach()
                obj._values[prop.name] = value
        for name in values:
            if name not in cls._property_map():
                self.warn(f"{cls.__name__}: kept unknown property {name!r} verbatim")
                obj._preserved[name] = values[name]
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
        if isinstance(uid, str) and uid in self.objects:
            return self.objects[uid]
        if isinstance(uid, str) and self.within is not None:
            for obj in self.within.walk():
                if obj.uid == uid:
                    return obj
        raise KeyError(f"no object with uid {uid!r}")

    def resolve(self) -> None:
        """Resolve every deferred reference; dangling ones become warnings."""
        for obj, prop, data in self._pending:
            try:
                obj._values[prop.name] = prop.resolve(data, self)
            except (KeyError, TypeError) as exc:
                self.warn(
                    f"{obj.label()}: reference {prop.name} is unresolved and kept for a "
                    f"later save ({exc})"
                )
                obj._preserved[prop.name] = data
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
    """Set ``abstract = True`` in a class body to mark that class (not its
    subclasses) abstract; read it with :meth:`is_abstract`."""

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

    @classmethod
    def is_abstract(cls) -> bool:
        """Whether this class itself is abstract (not offered as a strategy option)."""
        return bool(cls.__dict__.get("abstract", False))

    def _setup(self) -> None:
        self._uid = uuid.uuid4().hex
        self._preserved: dict[str, Any] = {}
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

    def _check_adoption(self, obj: ModelObject) -> None:
        """Refuse to make ``obj`` a descendant of this object if that forms a cycle."""
        node: ModelObject | None = self
        while node is not None:
            if node is obj:
                raise ValueError(f"{obj.label()} cannot own itself or one of its owners")
            node = node._parent

    def _adopt(self, obj: ModelObject) -> None:
        """Called on the root when ``obj`` joins its tree (the project tracks uids)."""

    def _release(self, obj: ModelObject) -> None:
        """Called on the root when ``obj`` leaves its tree."""

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

    def required_references(self) -> frozenset[str]:
        """Return the names of the references this object cannot do without now.

        A cascading removal that takes the target of one of them removes this
        object too. The base returns the :class:`~cabledyn.project.Ref`
        properties declared ``required``; subclasses add conditional ones.
        """
        return frozenset(
            prop.name for prop in type(self).properties() if isinstance(prop, Ref) and prop.required
        )

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
            self._check_adoption(value)
            root = self.root()
            root._release(old)
            try:
                root._adopt(value)
            except ValueError:
                root._adopt(old)
                raise
            old._detach()
            value._attach(self, prop.name)
        self._values[prop.name] = value
        self._preserved.pop(prop.name, None)
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

        References are written as the target's uid. Raw values kept from a
        project file (see :class:`SerialContext`) are written back unchanged.

        Raises
        ------
        TypeError
            If the class is not registered with :func:`model_type` (it would
            be read back as another class).
        """
        cls = type(self)
        if self.type_key not in REGISTRY or REGISTRY.get(self.type_key) is not cls:
            raise TypeError(
                f"{cls.__name__} is not registered with model_type and cannot be written"
            )
        context = ctx or SerialContext()
        data: dict[str, Any] = {
            "type": self.type_key,
            "uid": self._uid,
            "properties": {
                prop.name: prop.to_json(self._values[prop.name], context)
                for prop in type(self).properties()
            },
        }
        data["properties"].update(_json_copy(self._preserved))
        if self.deck_hints:
            data["deck_hints"] = _json_copy(self.deck_hints)
        return data

    @classmethod
    def from_dict(
        cls: type[M],
        data: dict[str, Any],
        *,
        fresh_uids: bool = True,
        within: ModelObject | None = None,
    ) -> M:
        """Rebuild an object (and its descendants) written by :meth:`to_dict`.

        The copy gets new uids by default, so it can be pasted next to the
        original. References to objects outside ``data`` are kept raw (see
        :class:`SerialContext`).

        Parameters
        ----------
        data : dict
            The dictionary.
        fresh_uids : bool
            Give the rebuilt objects new uids; ``False`` keeps the stored ones.
        within : ModelObject | None
            The graph the copy will join (normally the project): references to
            its objects resolve, so a pasted line keeps its end points.

        Returns
        -------
        ModelObject
            The rebuilt object.

        Raises
        ------
        TypeError
            If ``data`` does not describe a ``cls`` object.
        """
        ctx = SerialContext(fresh_uids=fresh_uids, within=within)
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
