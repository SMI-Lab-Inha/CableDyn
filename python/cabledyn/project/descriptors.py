# SPDX-License-Identifier: Apache-2.0
"""Typed, unit-aware, introspectable property descriptors.

A model class declares its data as class attributes::

    class Section(ModelObject):
        line_type = Ref(LineType)
        length = Quantity(LENGTH, 100.0, limit="section_length")
        segments = Integer(10, limit="segments")

Each descriptor knows its type, SI dimension, default, limits, role, group
and documentation. The descriptors drive the property panel, the project file,
validation layer 1, and the deck writer (which reads only ``physics``
properties). Values are stored in SI units.
"""

from __future__ import annotations

import copy
import math
import numbers
import re
from collections.abc import Callable, Sequence
from enum import Enum
from typing import TYPE_CHECKING, Any, Generic, TypeVar, overload

from cabledyn.project.collection import ObjectCollection
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.schema import deck_limit
from cabledyn.project.units import DIMENSIONLESS, Dimension

if TYPE_CHECKING:
    from cabledyn.project.base import ModelObject, SerialContext

__all__ = [
    "Bool",
    "Child",
    "Children",
    "Choice",
    "Colour",
    "FilePath",
    "Integer",
    "JsonValue",
    "OptionalBool",
    "OptionalColour",
    "OptionalInteger",
    "OptionalQuantity",
    "OptionalText",
    "OptionalVec3",
    "Property",
    "Quantity",
    "Ref",
    "RefList",
    "Role",
    "Strategy",
    "Text",
    "TextList",
    "Vec3",
    "Vector",
]

T = TypeVar("T")
ObjT = TypeVar("ObjT", bound="ModelObject")

Vec3Value = tuple[float, float, float]


class Role(str, Enum):
    """What a property describes, and whether the deck writer reads it."""

    PHYSICS = "physics"
    """Solver input, written to the deck."""
    CONSTRUCTION = "construction"
    """Construction or library data (chain grade, rope construction, floater
    hydrodynamics); not written to the deck today."""
    APPEARANCE = "appearance"
    """Cosmetic data for the 3D views; never written to the deck."""
    DERIVED = "derived"
    """Computed, read-only, shown for information."""
    META = "meta"
    """Identity and bookkeeping (name, description, tags, GUI data)."""


def _is_number(value: object) -> bool:
    return isinstance(value, numbers.Real) and not isinstance(value, bool)


def _same(a: object, b: object) -> bool:
    if isinstance(a, float) and isinstance(b, float):
        return a == b or (math.isnan(a) and math.isnan(b))
    if isinstance(a, tuple) and isinstance(b, tuple):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b, strict=True))
    if isinstance(a, ModelObjectLike) or isinstance(b, ModelObjectLike):
        return a is b
    return bool(a == b)


class ModelObjectLike:
    """Marker base of objects compared by identity (set on ``ModelObject``)."""


class Property(Generic[T]):
    """Base class of every property descriptor.

    Parameters
    ----------
    default : object
        Value of a new object.
    role : Role
        Physics, construction, appearance, derived, or meta.
    group : str
        Property-panel group, for example ``"Hydrodynamics"``.
    label : str
        Display label; derived from the attribute name when empty.
    doc : str
        One-sentence description for tooltips and generated documentation.
    read_only : bool
        Shown but not editable in the property panel.
    """

    kind = "property"
    """Short type name for editors and documentation."""

    def __init__(
        self,
        default: T,
        *,
        role: Role = Role.PHYSICS,
        group: str = "",
        label: str = "",
        doc: str = "",
        read_only: bool = False,
    ) -> None:
        self.default = default
        self.role = role
        self.group = group
        self.label = label
        self.doc = doc
        self.read_only = read_only
        self.name = ""

    def __set_name__(self, owner: type[Any], name: str) -> None:
        self.name = name
        if not self.label:
            self.label = name.replace("_", " ").capitalize()

    @overload
    def __get__(self, obj: None, owner: type[Any]) -> Property[T]: ...

    @overload
    def __get__(self, obj: ModelObject, owner: type[Any]) -> T: ...

    def __get__(self, obj: ModelObject | None, owner: type[Any]) -> Property[T] | T:
        if obj is None:
            return self
        value: T = obj._values[self.name]
        return value

    def __set__(self, obj: ModelObject, value: T) -> None:
        obj._assign(self, self.coerce(value))

    @property
    def is_reference(self) -> bool:
        """Whether the property refers to objects owned elsewhere."""
        return False

    @property
    def is_owner(self) -> bool:
        """Whether the property owns child objects."""
        return False

    @property
    def dimension(self) -> Dimension:
        """The SI dimension of the value (dimensionless for non-quantities)."""
        return DIMENSIONLESS

    def initial(self, obj: ModelObject) -> T:
        """Return the value of a new object (a copy of the default)."""
        return copy.deepcopy(self.default)

    def coerce(self, value: object) -> T:
        """Return ``value`` converted to the stored type.

        Raises
        ------
        TypeError
            If ``value`` has the wrong type.
        ValueError
            If ``value`` cannot be stored at all (for example a bad choice).
        """
        raise NotImplementedError

    def check(self, obj: ModelObject, value: T) -> list[Issue]:
        """Return the layer-1 issues of ``value`` (finiteness and limits)."""
        return []

    def same(self, a: T, b: T) -> bool:
        """Whether two values are equal for change detection."""
        return _same(a, b)

    def to_json(self, value: T, ctx: SerialContext) -> Any:
        """Return a JSON-compatible form of ``value``."""
        return value

    def from_json(self, data: Any, ctx: SerialContext, obj: ModelObject) -> T:
        """Rebuild a value written by :meth:`to_json`."""
        return self.coerce(data)

    def describe(self) -> dict[str, Any]:
        """Return a JSON-compatible description, for documentation and editors."""
        return {
            "name": self.name,
            "kind": self.kind,
            "label": self.label,
            "role": self.role.value,
            "group": self.group,
            "doc": self.doc,
            "dimension": self.dimension.key,
            "read_only": self.read_only,
        }

    def _issue(self, obj: ModelObject, message: str, *, layer: int = 1) -> Issue:
        return Issue(Severity.ERROR, message, obj, self.name, layer)

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.name!r}>"


class _NumberMixin:
    """Limits shared by :class:`Quantity` and :class:`Integer`."""

    minimum: float | None
    maximum: float | None
    limit: str | None

    def _limit_issues(self, prop: Property[Any], obj: ModelObject, value: float) -> list[Issue]:
        if not math.isfinite(value):
            return [prop._issue(obj, f"must be finite, got {value!r}")]
        found: list[Issue] = []
        if self.minimum is not None and value < self.minimum:
            found.append(prop._issue(obj, f"must be at least {self.minimum:g}"))
        if self.maximum is not None and value > self.maximum:
            found.append(prop._issue(obj, f"must be at most {self.maximum:g}"))
        if self.limit is not None:
            message = deck_limit(self.limit).violation(value)
            if message is not None:
                found.append(prop._issue(obj, message))
        return found


def _float(value: object, what: str) -> float:
    if not _is_number(value):
        raise TypeError(f"{what} must be a number, got {value!r}")
    return float(value)  # type: ignore[arg-type]


def _int(value: object, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        if isinstance(value, float) and value.is_integer():
            return int(value)
        raise TypeError(f"{what} must be an integer, got {value!r}")
    return int(value)


class Quantity(_NumberMixin, Property[float]):
    """A real number with a physical dimension, stored in SI units.

    Parameters
    ----------
    dimension : Dimension
        What the value measures.
    default : float
        Value of a new object.
    minimum, maximum : float | None
        Inclusive bounds checked by validation layer 1.
    limit : str | None
        Key of a shared deck limit (:data:`cabledyn.project.schema.DECK_LIMITS`).
    **kwargs
        Role, group, label, doc, read_only (see :class:`Property`).
    """

    kind = "quantity"

    def __init__(
        self,
        dimension: Dimension = DIMENSIONLESS,
        default: float = 0.0,
        *,
        minimum: float | None = None,
        maximum: float | None = None,
        limit: str | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(float(default), **kwargs)
        self._dimension = dimension
        self.minimum = minimum
        self.maximum = maximum
        self.limit = limit

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    def coerce(self, value: object) -> float:
        return _float(value, self.name)

    def check(self, obj: ModelObject, value: float) -> list[Issue]:
        return self._limit_issues(self, obj, value)

    def describe(self) -> dict[str, Any]:
        info = super().describe()
        info.update(minimum=self.minimum, maximum=self.maximum, limit=self.limit)
        return info


class OptionalQuantity(_NumberMixin, Property[float | None]):
    """A :class:`Quantity` that may be unset (``None``).

    Parameters are those of :class:`Quantity`; the default is ``None``.
    """

    kind = "quantity"

    def __init__(
        self,
        dimension: Dimension = DIMENSIONLESS,
        default: float | None = None,
        *,
        minimum: float | None = None,
        maximum: float | None = None,
        limit: str | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(None if default is None else float(default), **kwargs)
        self._dimension = dimension
        self.minimum = minimum
        self.maximum = maximum
        self.limit = limit

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    def coerce(self, value: object) -> float | None:
        return None if value is None else _float(value, self.name)

    def check(self, obj: ModelObject, value: float | None) -> list[Issue]:
        return [] if value is None else self._limit_issues(self, obj, value)


class Integer(_NumberMixin, Property[int]):
    """An integer.

    Parameters
    ----------
    default : int
        Value of a new object.
    minimum, maximum : int | None
        Inclusive bounds.
    limit : str | None
        Key of a shared deck limit.
    **kwargs
        See :class:`Property`.
    """

    kind = "integer"

    def __init__(
        self,
        default: int = 0,
        *,
        minimum: int | None = None,
        maximum: int | None = None,
        limit: str | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(default, **kwargs)
        self.minimum = minimum
        self.maximum = maximum
        self.limit = limit

    def coerce(self, value: object) -> int:
        return _int(value, self.name)

    def check(self, obj: ModelObject, value: int) -> list[Issue]:
        return self._limit_issues(self, obj, float(value))


class OptionalInteger(_NumberMixin, Property[int | None]):
    """An :class:`Integer` that may be unset (``None``)."""

    kind = "integer"

    def __init__(
        self,
        default: int | None = None,
        *,
        minimum: int | None = None,
        maximum: int | None = None,
        limit: str | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(default, **kwargs)
        self.minimum = minimum
        self.maximum = maximum
        self.limit = limit

    def coerce(self, value: object) -> int | None:
        return None if value is None else _int(value, self.name)

    def check(self, obj: ModelObject, value: int | None) -> list[Issue]:
        return [] if value is None else self._limit_issues(self, obj, float(value))


class Bool(Property[bool]):
    """A Boolean flag."""

    kind = "bool"

    def __init__(self, default: bool = False, **kwargs: Any) -> None:
        super().__init__(default, **kwargs)

    def coerce(self, value: object) -> bool:
        if not isinstance(value, bool):
            raise TypeError(f"{self.name} must be True or False, got {value!r}")
        return value


class OptionalBool(Property[bool | None]):
    """A Boolean flag that may be unset (``None``: the native default applies)."""

    kind = "bool"

    def __init__(self, default: bool | None = None, **kwargs: Any) -> None:
        super().__init__(default, **kwargs)

    def coerce(self, value: object) -> bool | None:
        if value is not None and not isinstance(value, bool):
            raise TypeError(f"{self.name} must be True, False, or None, got {value!r}")
        return value


class Choice(Property[str]):
    """One of a fixed set of strings.

    Parameters
    ----------
    choices : Sequence[str]
        The admissible values.
    default : str
        Value of a new object; the first choice when empty.
    **kwargs
        See :class:`Property`.
    """

    kind = "choice"

    def __init__(self, choices: Sequence[str], default: str = "", **kwargs: Any) -> None:
        if not choices:
            raise ValueError("a Choice needs at least one choice")
        super().__init__(default or choices[0], **kwargs)
        self.choices = tuple(choices)
        if self.default not in self.choices:
            raise ValueError(f"default {self.default!r} is not one of {self.choices}")

    def coerce(self, value: object) -> str:
        if not isinstance(value, str):
            raise TypeError(f"{self.name} must be a string, got {value!r}")
        if value not in self.choices:
            raise ValueError(f"{self.name} must be one of {', '.join(self.choices)}; got {value!r}")
        return value

    def describe(self) -> dict[str, Any]:
        info = super().describe()
        info["choices"] = list(self.choices)
        return info


class Text(Property[str]):
    """A free-text string (single line unless ``multiline``)."""

    kind = "text"

    def __init__(self, default: str = "", *, multiline: bool = False, **kwargs: Any) -> None:
        super().__init__(default, **kwargs)
        self.multiline = multiline

    def coerce(self, value: object) -> str:
        if not isinstance(value, str):
            raise TypeError(f"{self.name} must be a string, got {value!r}")
        return value


class OptionalText(Property[str | None]):
    """A string that may be unset (``None``)."""

    kind = "text"

    def __init__(self, default: str | None = None, **kwargs: Any) -> None:
        super().__init__(default, **kwargs)

    def coerce(self, value: object) -> str | None:
        if value is not None and not isinstance(value, str):
            raise TypeError(f"{self.name} must be a string or None, got {value!r}")
        return value


class FilePath(OptionalText):
    """A file path, relative to the project's deck folder or absolute; may be unset."""

    kind = "file"


class TextList(Property[tuple[str, ...]]):
    """An ordered tuple of strings (for example tags)."""

    kind = "text_list"

    def __init__(self, default: Sequence[str] = (), **kwargs: Any) -> None:
        super().__init__(tuple(default), **kwargs)

    def coerce(self, value: object) -> tuple[str, ...]:
        if isinstance(value, str) or not isinstance(value, (list, tuple)):
            raise TypeError(f"{self.name} must be a sequence of strings, got {value!r}")
        if not all(isinstance(item, str) for item in value):
            raise TypeError(f"{self.name} must hold only strings")
        return tuple(value)

    def to_json(self, value: tuple[str, ...], ctx: SerialContext) -> Any:
        return list(value)


_COLOUR = re.compile(r"^#[0-9a-fA-F]{6}$")


def _colour(value: object, what: str) -> str:
    if not isinstance(value, str) or _COLOUR.match(value) is None:
        raise ValueError(f"{what} must be a colour '#rrggbb', got {value!r}")
    return value.lower()


class Colour(Property[str]):
    """An sRGB colour written ``#rrggbb``; an appearance property by default."""

    kind = "colour"

    def __init__(self, default: str = "#808080", **kwargs: Any) -> None:
        kwargs.setdefault("role", Role.APPEARANCE)
        super().__init__(_colour(default, "default"), **kwargs)

    def coerce(self, value: object) -> str:
        return _colour(value, self.name)


class OptionalColour(Property[str | None]):
    """A :class:`Colour` that may be unset (``None``: inherit)."""

    kind = "colour"

    def __init__(self, default: str | None = None, **kwargs: Any) -> None:
        kwargs.setdefault("role", Role.APPEARANCE)
        super().__init__(default, **kwargs)

    def coerce(self, value: object) -> str | None:
        return None if value is None else _colour(value, self.name)


def _floats(value: object, what: str) -> tuple[float, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{what} must be a sequence of numbers, got {value!r}")
    return tuple(_float(item, what) for item in value)


class Vector(_NumberMixin, Property[tuple[float, ...]]):
    """A tuple of real numbers of one dimension, with admissible lengths.

    Used for deck fields that take one or several values (``a|b|c``), such
    as a body's drag areas.

    Parameters
    ----------
    dimension : Dimension
        Dimension of every component.
    default : Sequence[float]
        Value of a new object.
    lengths : Sequence[int]
        Admissible numbers of components.
    limit : str | None
        Shared deck limit applied to every component.
    **kwargs
        See :class:`Property`.
    """

    kind = "vector"

    def __init__(
        self,
        dimension: Dimension = DIMENSIONLESS,
        default: Sequence[float] = (0.0,),
        *,
        lengths: Sequence[int] = (1, 3),
        limit: str | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(tuple(float(v) for v in default), **kwargs)
        self._dimension = dimension
        self.lengths = tuple(lengths)
        self.minimum = None
        self.maximum = None
        self.limit = limit

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    def coerce(self, value: object) -> tuple[float, ...]:
        values = _floats(value, self.name)
        if len(values) not in self.lengths:
            sizes = " or ".join(str(n) for n in self.lengths)
            raise ValueError(f"{self.name} must have {sizes} values, got {len(values)}")
        return values

    def check(self, obj: ModelObject, value: tuple[float, ...]) -> list[Issue]:
        found: list[Issue] = []
        for item in value:
            found.extend(self._limit_issues(self, obj, item))
        return found

    def to_json(self, value: tuple[float, ...], ctx: SerialContext) -> Any:
        return list(value)

    def describe(self) -> dict[str, Any]:
        info = super().describe()
        info["lengths"] = list(self.lengths)
        return info


class Vec3(_NumberMixin, Property[Vec3Value]):
    """A three-component vector ``(x, y, z)`` of one dimension."""

    kind = "vec3"

    def __init__(
        self,
        dimension: Dimension = DIMENSIONLESS,
        default: Sequence[float] = (0.0, 0.0, 0.0),
        *,
        limit: str | None = None,
        **kwargs: Any,
    ) -> None:
        x, y, z = (float(v) for v in default)
        super().__init__((x, y, z), **kwargs)
        self._dimension = dimension
        self.minimum = None
        self.maximum = None
        self.limit = limit

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    def coerce(self, value: object) -> Vec3Value:
        values = _floats(value, self.name)
        if len(values) != 3:
            raise ValueError(f"{self.name} must have three values, got {len(values)}")
        return (values[0], values[1], values[2])

    def check(self, obj: ModelObject, value: Vec3Value) -> list[Issue]:
        found: list[Issue] = []
        for item in value:
            found.extend(self._limit_issues(self, obj, item))
        return found

    def to_json(self, value: Vec3Value, ctx: SerialContext) -> Any:
        return list(value)


class OptionalVec3(_NumberMixin, Property[Vec3Value | None]):
    """A :class:`Vec3` that may be unset (``None``)."""

    kind = "vec3"

    def __init__(
        self,
        dimension: Dimension = DIMENSIONLESS,
        default: Sequence[float] | None = None,
        *,
        limit: str | None = None,
        **kwargs: Any,
    ) -> None:
        start: Vec3Value | None = None
        if default is not None:
            x, y, z = (float(v) for v in default)
            start = (x, y, z)
        super().__init__(start, **kwargs)
        self._dimension = dimension
        self.minimum = None
        self.maximum = None
        self.limit = limit

    @property
    def dimension(self) -> Dimension:
        return self._dimension

    def coerce(self, value: object) -> Vec3Value | None:
        if value is None:
            return None
        values = _floats(value, self.name)
        if len(values) != 3:
            raise ValueError(f"{self.name} must have three values, got {len(values)}")
        return (values[0], values[1], values[2])

    def check(self, obj: ModelObject, value: Vec3Value | None) -> list[Issue]:
        found: list[Issue] = []
        for item in value or ():
            found.extend(self._limit_issues(self, obj, item))
        return found

    def to_json(self, value: Vec3Value | None, ctx: SerialContext) -> Any:
        return None if value is None else list(value)


class JsonValue(Property[Any]):
    """Opaque JSON-compatible data (GUI layouts, plot specifications, ...).

    The value is deep-copied on assignment; it is never written to a deck.
    """

    kind = "json"

    def __init__(self, default: Any = None, **kwargs: Any) -> None:
        kwargs.setdefault("role", Role.META)
        super().__init__(default, **kwargs)

    def coerce(self, value: object) -> Any:
        return copy.deepcopy(value)


class Ref(Property[ObjT | None], Generic[ObjT]):
    """A reference to another object of the model (not to a name or an id).

    Renaming the target never breaks the reference. Removing a referenced
    object either fails or cascades (see
    :class:`cabledyn.project.commands.RemoveObject`): a referrer whose
    ``required`` reference loses its target is removed with it, and an
    optional reference is cleared.

    Parameters
    ----------
    target : type
        The class (or common base class) of admissible targets.
    required : bool
        A missing target is an error (validation layer 2).
    **kwargs
        See :class:`Property`.
    """

    kind = "ref"

    def __init__(self, target: type[ObjT], *, required: bool = True, **kwargs: Any) -> None:
        super().__init__(None, **kwargs)
        self.target = target
        self.required = required

    @property
    def is_reference(self) -> bool:
        return True

    def initial(self, obj: ModelObject) -> ObjT | None:
        return None

    def coerce(self, value: object) -> ObjT | None:
        if value is not None and not isinstance(value, self.target):
            raise TypeError(
                f"{self.name} must refer to a {self.target.__name__}, got {type(value).__name__}"
            )
        return value

    def same(self, a: ObjT | None, b: ObjT | None) -> bool:
        return a is b

    def check(self, obj: ModelObject, value: ObjT | None) -> list[Issue]:
        if value is None:
            if self.required:
                return [self._issue(obj, "refers to nothing", layer=2)]
            return []
        if value.root() is not obj.root():
            return [self._issue(obj, f"refers to {value.label()}, outside this model", layer=2)]
        return []

    def to_json(self, value: ObjT | None, ctx: SerialContext) -> Any:
        return None if value is None else value.uid

    def from_json(self, data: Any, ctx: SerialContext, obj: ModelObject) -> ObjT | None:
        if data is not None:
            ctx.defer(obj, self, data)
        return None

    def resolve(self, data: Any, ctx: SerialContext) -> ObjT | None:
        """Return the target of a uid written by :meth:`to_json`."""
        return self.coerce(ctx.lookup(data))

    def describe(self) -> dict[str, Any]:
        info = super().describe()
        info.update(target=self.target.__name__, required=self.required)
        return info


class RefList(Property[tuple[ObjT, ...]], Generic[ObjT]):
    """An ordered tuple of references to other objects.

    Parameters
    ----------
    target : type
        Class of admissible targets.
    min_items : int
        Fewest targets; a removal that leaves fewer removes the referrer
        when cascading.
    **kwargs
        See :class:`Property`.
    """

    kind = "ref_list"

    def __init__(self, target: type[ObjT], *, min_items: int = 0, **kwargs: Any) -> None:
        super().__init__((), **kwargs)
        self.target = target
        self.min_items = min_items

    @property
    def is_reference(self) -> bool:
        return True

    def initial(self, obj: ModelObject) -> tuple[ObjT, ...]:
        return ()

    def coerce(self, value: object) -> tuple[ObjT, ...]:
        if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
            raise TypeError(f"{self.name} must be a sequence of objects, got {value!r}")
        items = tuple(value)
        for item in items:
            if not isinstance(item, self.target):
                raise TypeError(
                    f"{self.name} must refer to {self.target.__name__} objects, "
                    f"got {type(item).__name__}"
                )
        return items

    def same(self, a: tuple[ObjT, ...], b: tuple[ObjT, ...]) -> bool:
        return len(a) == len(b) and all(x is y for x, y in zip(a, b, strict=True))

    def check(self, obj: ModelObject, value: tuple[ObjT, ...]) -> list[Issue]:
        found: list[Issue] = []
        if len(value) < self.min_items:
            found.append(self._issue(obj, f"needs at least {self.min_items} item(s)", layer=2))
        for item in value:
            if item.root() is not obj.root():
                found.append(
                    self._issue(obj, f"refers to {item.label()}, outside this model", layer=2)
                )
        return found

    def to_json(self, value: tuple[ObjT, ...], ctx: SerialContext) -> Any:
        return [item.uid for item in value]

    def from_json(self, data: Any, ctx: SerialContext, obj: ModelObject) -> tuple[ObjT, ...]:
        if data:
            ctx.defer(obj, self, data)
        return ()

    def resolve(self, data: Any, ctx: SerialContext) -> tuple[ObjT, ...]:
        """Return the targets of uids written by :meth:`to_json`."""
        return self.coerce([ctx.lookup(uid) for uid in data])

    def describe(self) -> dict[str, Any]:
        info = super().describe()
        info.update(target=self.target.__name__, min_items=self.min_items)
        return info


class Child(Property[ObjT], Generic[ObjT]):
    """One owned sub-object, created with the owner and replaceable as a whole.

    Parameters
    ----------
    base : type
        Class of the sub-object (subclasses are admissible).
    factory : Callable[[], ModelObject] | None
        Creates the default sub-object; ``base`` itself when ``None``.
    **kwargs
        See :class:`Property`.
    """

    kind = "child"

    def __init__(
        self, base: type[ObjT], factory: Callable[[], ObjT] | None = None, **kwargs: Any
    ) -> None:
        super().__init__(None, **kwargs)  # type: ignore[arg-type]
        self.base = base
        self.factory: Callable[[], ObjT] = factory or base

    @property
    def is_owner(self) -> bool:
        return True

    def initial(self, obj: ModelObject) -> ObjT:
        child = self.factory()
        child._attach(obj, self.name)
        return child

    def coerce(self, value: object) -> ObjT:
        if not isinstance(value, self.base):
            raise TypeError(
                f"{self.name} must be a {self.base.__name__}, got {type(value).__name__}"
            )
        return value

    def same(self, a: ObjT, b: ObjT) -> bool:
        return a is b

    def to_json(self, value: ObjT, ctx: SerialContext) -> Any:
        return value.to_dict(ctx)

    def from_json(self, data: Any, ctx: SerialContext, obj: ModelObject) -> ObjT:
        built = ctx.build(data, owner=obj, slot=self.name)
        if built is None or not isinstance(built, self.base):
            return self.initial(obj)
        built._attach(obj, self.name)
        return built

    def describe(self) -> dict[str, Any]:
        info = super().describe()
        info["base"] = self.base.__name__
        return info


class Strategy(Child[ObjT], Generic[ObjT]):
    """A polymorphic sub-object such as an axial model or a wave model.

    Any registered subclass of ``base`` can be swapped in (see
    :class:`cabledyn.project.commands.ReplaceStrategy`); :meth:`options` lists
    them for an editor.
    """

    kind = "strategy"

    def options(self) -> list[type[ModelObject]]:
        """Return the registered concrete classes admissible here."""
        from cabledyn.project.base import REGISTRY

        return REGISTRY.subclasses_of(self.base)


class Children(Property[ObjectCollection[ObjT]], Generic[ObjT]):
    """An owned, ordered, observable collection of sub-objects.

    The collection itself cannot be replaced; insert, remove, and move items
    through it (or through commands, for undo).

    Parameters
    ----------
    item_type : type
        Class of the items.
    min_items : int
        Fewest items; a cascading removal that leaves fewer removes the owner.
    **kwargs
        See :class:`Property`.
    """

    kind = "children"

    def __init__(self, item_type: type[ObjT], *, min_items: int = 0, **kwargs: Any) -> None:
        super().__init__(None, **kwargs)  # type: ignore[arg-type]
        self.item_type = item_type
        self.min_items = min_items

    @property
    def is_owner(self) -> bool:
        return True

    def __set__(self, obj: ModelObject, value: ObjectCollection[ObjT]) -> None:
        raise AttributeError(f"{self.name} is a collection; insert or remove its items")

    def initial(self, obj: ModelObject) -> ObjectCollection[ObjT]:
        return ObjectCollection(obj, self.name, self.item_type)

    def coerce(self, value: object) -> ObjectCollection[ObjT]:
        raise AttributeError(f"{self.name} is a collection; insert or remove its items")

    def same(self, a: ObjectCollection[ObjT], b: ObjectCollection[ObjT]) -> bool:
        return a is b

    def check(self, obj: ModelObject, value: ObjectCollection[ObjT]) -> list[Issue]:
        if len(value) < self.min_items:
            return [self._issue(obj, f"needs at least {self.min_items} item(s)", layer=2)]
        return []

    def to_json(self, value: ObjectCollection[ObjT], ctx: SerialContext) -> Any:
        return [item.to_dict(ctx) for item in value]

    def from_json(self, data: Any, ctx: SerialContext, obj: ModelObject) -> ObjectCollection[ObjT]:
        collection: ObjectCollection[ObjT] = obj._values[self.name]
        for index, item in enumerate(data):
            built = ctx.build(item, owner=obj, slot=self.name, index=index)
            if built is None:
                continue
            if not isinstance(built, self.item_type):
                ctx.warn(f"{self.name}: dropped a {type(built).__name__}, not a valid item")
                continue
            collection._items.append(built)
            built._attach(obj, self.name)
        return collection

    def describe(self) -> dict[str, Any]:
        info = super().describe()
        info.update(item_type=self.item_type.__name__, min_items=self.min_items)
        return info
