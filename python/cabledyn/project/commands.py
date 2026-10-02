# SPDX-License-Identifier: Apache-2.0
"""Undoable commands and the command stack.

Every interactive change to a project is a :class:`Command` pushed onto a
:class:`CommandStack`, which executes it and records it for undo and redo.
Undoing a command restores the project exactly (its
:meth:`~cabledyn.project.ModelObject.to_dict` is unchanged), and object
identity is preserved, so references to removed-and-restored objects stay
valid.

Examples
--------
>>> from cabledyn.project import CommandStack, FixedPoint, Project
>>> from cabledyn.project.commands import AddObject, SetProperty
>>> project = Project()
>>> stack = CommandStack()
>>> anchor = FixedPoint("anchor")
>>> stack.push(AddObject(project.points, anchor))
>>> stack.push(SetProperty(anchor, "position", (400.0, 0.0, -50.0)))
>>> stack.undo()
>>> anchor.position
(0.0, 0.0, 0.0)
"""

from __future__ import annotations

import contextlib
import copy
from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from cabledyn.project.base import ModelObject, collection_of
from cabledyn.project.collection import ObjectCollection
from cabledyn.project.descriptors import Child, Children
from cabledyn.project.lines import Line, Section
from cabledyn.project.points import LineEndTarget
from cabledyn.project.references import RemovalPlan, plan_removal

if TYPE_CHECKING:
    from cabledyn.project.project import Project

__all__ = [
    "AddObject",
    "ApplyLibraryItem",
    "Command",
    "CommandStack",
    "ImportDeck",
    "MacroCommand",
    "MoveInCollection",
    "Reconnect",
    "RemoveObject",
    "ReplaceStrategy",
    "SetProperty",
    "SplitSection",
    "StackChanged",
]


class Command:
    """Base class of an undoable change.

    Subclasses implement :meth:`do` and :meth:`undo`; :meth:`do` is called
    again for a redo and must give the same result.

    Attributes
    ----------
    text : str
        Description for the Edit menu ("Set length of Section").
    """

    text = "Edit"

    def do(self) -> None:
        """Apply the change."""
        raise NotImplementedError

    def undo(self) -> None:
        """Revert the change."""
        raise NotImplementedError

    def merge_with(self, other: Command) -> bool:
        """Absorb ``other`` (already applied) into this command, if compatible.

        Returns
        -------
        bool
            Whether ``other`` was absorbed; the stack then discards it.
        """
        return False

    def affected(self) -> tuple[ModelObject, ...]:
        """Return the objects the command changes."""
        return ()

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.text!r}>"


class SetProperty(Command):
    """Set one property of one object.

    Parameters
    ----------
    obj : ModelObject
        The object.
    name : str
        Property name.
    value : object
        New value (SI units); checked against the property type at once.
    merge : bool
        Merge with an immediately following mergeable ``SetProperty`` of the
        same object and property (one undo step for a drag).
    text : str | None
        Menu text; generated when ``None``.

    Raises
    ------
    KeyError
        If the object has no such property.
    TypeError, ValueError
        If the value does not fit the property.
    AttributeError
        If the property is read-only or a collection.
    """

    def __init__(
        self,
        obj: ModelObject,
        name: str,
        value: object,
        *,
        merge: bool = False,
        text: str | None = None,
    ) -> None:
        self.obj = obj
        self.prop = type(obj).get_property(name)
        if self.prop.read_only or isinstance(self.prop, Children):
            raise AttributeError(f"{name} cannot be set")
        self.new = self.prop.coerce(value)
        self.old: Any = obj._values[name]
        self._applied = False
        self.merge = merge
        self.text = text or f"Set {self.prop.label.lower()} of {obj.label()}"

    def do(self) -> None:
        if not self._applied:
            self.old = self.obj._values[self.prop.name]
            self._applied = True
        self.obj._assign(self.prop, self.new)

    def undo(self) -> None:
        self.obj._assign(self.prop, self.old)

    def merge_with(self, other: Command) -> bool:
        if (
            self.merge
            and isinstance(other, SetProperty)
            and other.merge
            and other.obj is self.obj
            and other.prop is self.prop
        ):
            self.new = other.new
            return True
        return False

    def affected(self) -> tuple[ModelObject, ...]:
        return (self.obj,)


class ReplaceStrategy(SetProperty):
    """Replace a strategy sub-object (axial model, wave model, seabed, ...).

    Values of properties the old and new classes share are copied to the new
    object unless ``copy_common`` is false.

    Parameters
    ----------
    obj : ModelObject
        The owner.
    name : str
        The :class:`~cabledyn.project.Strategy` (or :class:`~cabledyn.project.Child`)
        property.
    new : ModelObject | type
        The new sub-object, or a class to instantiate.
    copy_common : bool
        Copy shared property values from the current sub-object.
    """

    def __init__(
        self,
        obj: ModelObject,
        name: str,
        new: ModelObject | type[ModelObject],
        *,
        copy_common: bool = True,
    ) -> None:
        prop = type(obj).get_property(name)
        if not isinstance(prop, Child):
            raise TypeError(f"{name} is not a sub-object property")
        replacement = new() if isinstance(new, type) else new
        if copy_common:
            current = obj._values[name]
            shared = {p.name for p in type(current).properties()}
            for target in type(replacement).properties():
                if target.name in shared and not target.is_owner and not target.read_only:
                    with contextlib.suppress(TypeError, ValueError):
                        value = target.coerce(current._values[target.name])
                        replacement._values[target.name] = value
        super().__init__(
            obj, name, replacement, text=f"Change {prop.label.lower()} of {obj.label()}"
        )


class Reconnect(SetProperty):
    """Attach one end of a line to another point or rod end.

    Parameters
    ----------
    line : Line
        The line.
    end : str
        ``"A"`` or ``"B"``.
    target : LineEndTarget
        The new attachment.
    """

    def __init__(self, line: Line, end: str, target: LineEndTarget) -> None:
        letter = end.upper()
        if letter not in {"A", "B"}:
            raise ValueError(f"end must be A or B, got {end!r}")
        super().__init__(
            line,
            "end_a" if letter == "A" else "end_b",
            target,
            text=f"Connect End {letter} of {line.label()} to {target.label()}",
        )


class AddObject(Command):
    """Insert a new object into a collection.

    Parameters
    ----------
    collection : ObjectCollection
        The target collection, for example ``project.points``.
    obj : ModelObject
        The new object (without a parent).
    index : int | None
        Position; appended when ``None``.
    """

    def __init__(
        self, collection: ObjectCollection[Any], obj: ModelObject, index: int | None = None
    ) -> None:
        if not isinstance(obj, collection.item_type):
            raise TypeError(
                f"{collection.name} holds {collection.item_type.__name__}, not {type(obj).__name__}"
            )
        if obj.parent is not None:
            raise ValueError(f"{obj.label()} already belongs to {obj.parent.label()}")
        self.collection = collection
        self.obj = obj
        self.index = index
        self.text = f"Add {obj.label()}"

    def do(self) -> None:
        if self.index is None:
            self.index = len(self.collection)
        self.collection.insert(self.index, self.obj)

    def undo(self) -> None:
        self.collection.remove(self.obj)

    def affected(self) -> tuple[ModelObject, ...]:
        return (self.collection.owner, self.obj)


class RemoveObject(Command):
    """Remove objects, refusing or cascading when they are still referred to.

    The cascade (see :mod:`cabledyn.project.references`) is worked out when
    the command is first applied and captured, so undo restores every removed
    dependant and every cleared reference.

    Parameters
    ----------
    root : ModelObject
        The object graph, normally the project.
    objects : ModelObject | Sequence[ModelObject]
        The object(s) to remove.
    cascade : bool
        Remove or update dependants instead of refusing.

    Raises
    ------
    ObjectInUseError
        When applied without ``cascade`` to a referenced object.
    """

    def __init__(
        self,
        root: ModelObject,
        objects: ModelObject | Sequence[ModelObject],
        *,
        cascade: bool = False,
    ) -> None:
        self.root = root
        self.targets = [objects] if isinstance(objects, ModelObject) else list(objects)
        self.cascade = cascade
        self.plan: RemovalPlan | None = None
        self._positions: list[tuple[ObjectCollection[Any], ModelObject, int]] = []
        first = self.targets[0].label() if self.targets else "nothing"
        self.text = f"Remove {first}" if len(self.targets) < 2 else f"Remove {len(self.targets)}"

    def do(self) -> None:
        if self.plan is None:
            self.plan = plan_removal(self.root, self.targets, cascade=self.cascade)
        for obj, name, _old, new in self.plan.reference_updates:
            obj._assign(type(obj).get_property(name), new)
        self._positions = []
        for obj in self.plan.removed:
            collection = collection_of(obj)
            assert collection is not None
            self._positions.append((collection, obj, collection.remove(obj)))

    def undo(self) -> None:
        for collection, obj, index in reversed(self._positions):
            collection.insert(index, obj)
        assert self.plan is not None
        for obj, name, old, _new in reversed(self.plan.reference_updates):
            obj._assign(type(obj).get_property(name), old)

    def affected(self) -> tuple[ModelObject, ...]:
        if self.plan is None:
            return tuple(self.targets)
        return tuple(self.plan.removed) + tuple(obj for obj, *_ in self.plan.reference_updates)


class MoveInCollection(Command):
    """Move an object to another position within its collection.

    Parameters
    ----------
    obj : ModelObject
        The object.
    index : int
        New position.
    """

    def __init__(self, obj: ModelObject, index: int) -> None:
        collection = collection_of(obj)
        if collection is None:
            raise ValueError(f"{obj.label()} is not in a collection")
        self.collection = collection
        self.obj = obj
        self.index = index
        self.old = collection.index(obj)
        self.text = f"Move {obj.label()}"

    def do(self) -> None:
        self.old = self.collection.index(self.obj)
        self.collection.move(self.obj, self.index)

    def undo(self) -> None:
        self.collection.move(self.obj, self.old)

    def affected(self) -> tuple[ModelObject, ...]:
        return (self.obj,)


class MacroCommand(Command):
    """Several commands applied and undone as one.

    If a command fails while the macro is applied, the commands already
    applied are undone and the error is raised.

    Parameters
    ----------
    text : str
        Menu text.
    commands : Sequence[Command]
        The commands, in order.
    """

    def __init__(self, text: str, commands: Sequence[Command] = ()) -> None:
        self.text = text
        self.commands = list(commands)

    def do(self) -> None:
        done: list[Command] = []
        try:
            for command in self.commands:
                command.do()
                done.append(command)
        except Exception:
            for command in reversed(done):
                command.undo()
            raise

    def undo(self) -> None:
        for command in reversed(self.commands):
            command.undo()

    def affected(self) -> tuple[ModelObject, ...]:
        found: list[ModelObject] = []
        for command in self.commands:
            for obj in command.affected():
                if all(obj is not other for other in found):
                    found.append(obj)
        return tuple(found)


class SplitSection(MacroCommand):
    """Split a section in two at an arc length within it.

    The elements are shared in proportion to the lengths (at least one each).

    Parameters
    ----------
    section : Section
        The section; it must belong to a line.
    at : float
        Distance from the section's start, in m, strictly inside it.

    Raises
    ------
    ValueError
        If ``at`` is not inside the section or the section has one element.
    """

    def __init__(self, section: Section, at: float) -> None:
        collection = collection_of(section)
        if collection is None:
            raise ValueError("the section does not belong to a line")
        if not 0.0 < at < section.length:
            raise ValueError(f"split point {at!r} is not inside the section")
        if section.segments < 2:
            raise ValueError("a section with one element cannot be split")
        first = max(1, min(section.segments - 1, round(section.segments * at / section.length)))
        second = Section(
            section.name,
            line_type=section.line_type,
            length=section.length - at,
            segments=section.segments - first,
            colour=section.colour,
        )
        super().__init__(
            f"Split {section.label()}",
            [
                SetProperty(section, "length", at),
                SetProperty(section, "segments", first),
                AddObject(collection, second, collection.index(section) + 1),
            ],
        )
        self.new_section = second


class ApplyLibraryItem(MacroCommand):
    """Set several properties of an object from a library entry.

    Parameters
    ----------
    obj : ModelObject
        The object, typically a line type.
    values : Mapping[str, object]
        Property values by name (SI units).
    text : str | None
        Menu text.
    """

    def __init__(
        self, obj: ModelObject, values: Mapping[str, object], text: str | None = None
    ) -> None:
        super().__init__(
            text or f"Apply library data to {obj.label()}",
            [SetProperty(obj, name, value) for name, value in values.items()],
        )


class ImportDeck(Command):
    """Add the objects of a deck (or of another project) to a project.

    Libraries, objects and rows are appended; the environment, seabed,
    motion and settings of the imported model (and its option spellings and
    order) replace the project's when ``replace_environment`` is true. Name
    clashes of line or rod types are left for validation to report. The
    objects move out of ``source`` when the command is first applied.

    Parameters
    ----------
    project : Project
        The receiving project.
    source : Project
        A project built by :class:`~cabledyn.project.DeckReader`.
    replace_environment : bool
        Also take the environment, seabed, motion, settings and option hints.
    """

    _SINGLE = ("environment", "seabed", "motion", "settings")

    def __init__(
        self, project: Project, source: Project, *, replace_environment: bool = True
    ) -> None:
        self.project = project
        self.source = source
        self.replace_environment = replace_environment
        self.text = f"Import {source.title}"
        self._inner: MacroCommand | None = None
        self._hints: tuple[Any, Any] | None = None

    def _build(self) -> MacroCommand:
        project, source = self.project, self.source
        commands: list[Command] = []
        for prop in type(source).properties():
            if not isinstance(prop, Children):
                continue
            target = project._values[prop.name]
            for item in list(source._values[prop.name]):
                source._values[prop.name].remove(item)
                commands.append(AddObject(target, item))
        for channel in list(source.outputs.channels):
            source.outputs.channels.remove(channel)
            commands.append(AddObject(project.outputs.channels, channel))
        if self.replace_environment:
            for name in self._SINGLE:
                replacement = source._values[name]
                replacement._detach()
                source._values[name] = type(source).get_property(name).initial(source)
                commands.append(SetProperty(project, name, replacement))
            self._hints = (
                project.deck_hints.get("options"),
                copy.deepcopy(source.deck_hints.get("options")),
            )
        return MacroCommand(self.text, commands)

    def do(self) -> None:
        if self._inner is None:
            self._inner = self._build()
        self._inner.do()
        if self._hints is not None:
            self._set_hints(self._hints[1])

    def undo(self) -> None:
        assert self._inner is not None
        self._inner.undo()
        if self._hints is not None:
            self._set_hints(self._hints[0])

    def _set_hints(self, hints: Any) -> None:
        if hints is None:
            self.project.deck_hints.pop("options", None)
        else:
            self.project.deck_hints["options"] = copy.deepcopy(hints)

    def affected(self) -> tuple[ModelObject, ...]:
        return (self.project,) if self._inner is None else (self.project, *self._inner.affected())


@dataclass(frozen=True)
class StackChanged:
    """Published by :class:`CommandStack` after every change of its state.

    Attributes
    ----------
    action : str
        ``"push"``, ``"merge"``, ``"undo"``, ``"redo"``, ``"clean"``, or ``"clear"``.
    command : Command | None
        The command concerned.
    """

    action: str
    command: Command | None


class CommandStack:
    """Undo and redo history of commands.

    Parameters
    ----------
    limit : int
        Most commands kept (oldest dropped first); ``0`` keeps all.
    """

    def __init__(self, limit: int = 0) -> None:
        self._commands: list[Command] = []
        self._index = 0
        self._clean: int | None = 0
        self._limit = limit
        self._macro: list[MacroCommand] = []
        self._busy = False
        self._listeners: list[Callable[[StackChanged], None]] = []

    # ------------------------------------------------------------------ state

    @property
    def index(self) -> int:
        """Number of applied commands in the history."""
        return self._index

    def __len__(self) -> int:
        return len(self._commands)

    @property
    def can_undo(self) -> bool:
        """Whether there is a command to undo."""
        return self._index > 0 and not self._macro

    @property
    def can_redo(self) -> bool:
        """Whether there is a command to redo."""
        return self._index < len(self._commands) and not self._macro

    @property
    def undo_text(self) -> str:
        """Text of the command :meth:`undo` would revert ("" if none)."""
        return self._commands[self._index - 1].text if self.can_undo else ""

    @property
    def redo_text(self) -> str:
        """Text of the command :meth:`redo` would apply ("" if none)."""
        return self._commands[self._index].text if self.can_redo else ""

    @property
    def is_clean(self) -> bool:
        """Whether the state equals the last :meth:`set_clean` (for example a save)."""
        return self._clean == self._index

    def commands(self) -> tuple[Command, ...]:
        """Return the recorded commands, oldest first."""
        return tuple(self._commands)

    def subscribe(self, listener: Callable[[StackChanged], None]) -> Callable[[], None]:
        """Call ``listener`` after every stack change; returns an unsubscribe function."""
        self._listeners.append(listener)

        def cancel() -> None:
            if listener in self._listeners:
                self._listeners.remove(listener)

        return cancel

    def _notify(self, action: str, command: Command | None) -> None:
        for listener in tuple(self._listeners):
            listener(StackChanged(action, command))

    # ------------------------------------------------------------------ actions

    def _run(self, action: Callable[[], None]) -> None:
        if self._busy:
            raise RuntimeError(
                "the command stack is applying a command; a change listener must not push, "
                "undo or redo (schedule it after the current change instead)"
            )
        self._busy = True
        try:
            action()
        finally:
            self._busy = False

    def push(self, command: Command) -> None:
        """Apply ``command`` and record it (inside a macro: add it to the macro).

        Raises
        ------
        RuntimeError
            If called while a command is being applied (from a change listener).
        Exception
            Whatever the command raises; nothing is recorded then.
        """
        self._run(command.do)
        self._record(command)

    def _record(self, command: Command) -> None:
        if self._macro:
            self._macro[-1].commands.append(command)
            return
        del self._commands[self._index :]
        if self._clean is not None and self._clean > self._index:
            self._clean = None
        top = self._commands[-1] if self._commands else None
        if top is not None and self._clean != self._index and top.merge_with(command):
            self._notify("merge", top)
            return
        self._commands.append(command)
        self._index += 1
        if self._limit and len(self._commands) > self._limit:
            del self._commands[0]
            self._index -= 1
            self._clean = None if not self._clean else self._clean - 1
        self._notify("push", command)

    def undo(self) -> None:
        """Revert the last applied command (no-op if none)."""
        if not self.can_undo:
            return
        command = self._commands[self._index - 1]
        self._run(command.undo)
        self._index -= 1
        self._notify("undo", command)

    def redo(self) -> None:
        """Re-apply the next undone command (no-op if none)."""
        if not self.can_redo:
            return
        command = self._commands[self._index]
        self._run(command.do)
        self._index += 1
        self._notify("redo", command)

    def set_clean(self) -> None:
        """Mark the current state as saved."""
        self._clean = self._index
        self._notify("clean", None)

    def clear(self) -> None:
        """Forget the history (the model is unchanged)."""
        self._commands.clear()
        self._index = 0
        self._clean = 0
        self._notify("clear", None)

    @contextlib.contextmanager
    def macro(self, text: str) -> Iterator[MacroCommand]:
        """Group the commands pushed in a ``with`` block into one undo step.

        If the block raises, the commands pushed in it are undone and the
        error propagates.

        Parameters
        ----------
        text : str
            Menu text of the group.
        """
        group = MacroCommand(text)
        self._macro.append(group)
        try:
            yield group
        except BaseException:
            self._macro.pop()
            group.undo()
            raise
        self._macro.pop()
        if group.commands:
            self._record(group)
