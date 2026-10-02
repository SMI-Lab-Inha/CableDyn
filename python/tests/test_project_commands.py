# SPDX-License-Identifier: Apache-2.0
"""Commands, the undo stack, and reference semantics on rename and removal."""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

import pytest

from cabledyn.project import (
    AddObject,
    ApplyLibraryItem,
    BathymetrySeabed,
    Body,
    BodyPoint,
    ChainType,
    Command,
    CommandStack,
    Control,
    DeckReader,
    EndConnection,
    Failure,
    FixedPoint,
    FlatSeabed,
    FloaterPoint,
    GenericLineType,
    Group,
    ImportDeck,
    Line,
    MacroCommand,
    MoveInCollection,
    ObjectChannel,
    ObjectInUseError,
    Project,
    Reconnect,
    RemoveObject,
    ReplaceStrategy,
    Rigid6Body,
    Rod,
    RodPoint,
    RodType,
    Section,
    SetProperty,
    SplitSection,
    StackChanged,
    SyropeIC,
    ViscoelasticAxial,
    plan_removal,
    referrers,
)
from cabledyn.project.references import RemovalPlan

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def build() -> dict[str, Any]:
    project = Project("demo")
    chain = project.line_types.append(GenericLineType("chain", diameter=0.2))
    wire = project.line_types.append(GenericLineType("wire", diameter=0.1))
    fairlead = project.points.append(FloaterPoint("fairlead", position=(0, 0, -14)))
    anchor = project.points.append(FixedPoint("anchor", position=(400, 0, -50)))
    line = project.lines.append(Line("mooring", end_a=fairlead, end_b=anchor))
    s1 = line.sections.append(Section("s1", line_type=chain, length=200.0, segments=20))
    s2 = line.sections.append(Section("s2", line_type=wire, length=210.0, segments=21))
    second = project.lines.append(Line("second", end_a=fairlead, end_b=anchor))
    second.sections.append(Section("only", line_type=wire, length=410.0, segments=41))
    project.failures.append(Failure(point=fairlead, lines=[line, second]))
    project.controls.append(Control(lines=[line]))
    project.syrope_history.append(SyropeIC(lines=[second]))
    project.outputs.channels.append(ObjectChannel(target=line, quantity="FairTen"))
    group = project.studio.groups.append(Group("all", members=[line, second, anchor]))
    return {
        "project": project,
        "chain": chain,
        "wire": wire,
        "fairlead": fairlead,
        "anchor": anchor,
        "line": line,
        "second": second,
        "s1": s1,
        "s2": s2,
        "group": group,
    }


def check_undo_redo(project: Project, command: Command, stack: CommandStack | None = None) -> None:
    stack = stack or CommandStack()
    before = project.to_dict()
    stack.push(command)
    after = project.to_dict()
    assert after != before
    stack.undo()
    assert project.to_dict() == before
    stack.redo()
    assert project.to_dict() == after
    stack.undo()
    assert project.to_dict() == before


# --------------------------------------------------------------------------- basic commands


def test_set_property_and_merging() -> None:
    m = build()
    project, s1 = m["project"], m["s1"]
    check_undo_redo(project, SetProperty(s1, "length", 250.0))
    stack = CommandStack()
    stack.push(SetProperty(s1, "length", 201.0, merge=True))
    stack.push(SetProperty(s1, "length", 202.0, merge=True))
    stack.push(SetProperty(s1, "length", 203.0, merge=True))
    assert len(stack) == 1 and s1.length == 203.0
    stack.push(SetProperty(s1, "segments", 30, merge=True))
    assert len(stack) == 2
    stack.undo()
    stack.undo()
    assert s1.length == 200.0
    with pytest.raises(AttributeError, match="cannot be set"):
        SetProperty(project, "lines", [])
    with pytest.raises(TypeError):
        SetProperty(s1, "length", "long")
    with pytest.raises(KeyError):
        SetProperty(s1, "nothing", 1)
    assert SetProperty(s1, "length", 1.0).affected() == (s1,)
    assert "Set length of Section 's1'" in SetProperty(s1, "length", 1.0).text


def test_rename_keeps_every_reference() -> None:
    m = build()
    project, chain, fairlead = m["project"], m["chain"], m["fairlead"]
    stack = CommandStack()
    stack.push(SetProperty(chain, "name", "studless R4"))
    stack.push(SetProperty(fairlead, "name", "hang-off"))
    assert m["s1"].line_type is chain and m["line"].end_a is fairlead
    assert project.failures[0].point is fairlead
    stack.undo()
    assert fairlead.name == "fairlead"


def test_add_move_and_reconnect() -> None:
    m = build()
    project = m["project"]
    extra = FixedPoint("extra")
    check_undo_redo(project, AddObject(project.points, extra, 0))
    stack = CommandStack()
    stack.push(AddObject(project.points, extra))
    assert project.points[-1] is extra
    assert AddObject(project.points, FixedPoint()).affected()[0] is project
    with pytest.raises(TypeError, match="holds Point"):
        AddObject(project.points, Line())
    with pytest.raises(ValueError, match="already belongs"):
        AddObject(project.points, extra)
    check_undo_redo(project, MoveInCollection(extra, 0))
    assert MoveInCollection(extra, 0).affected() == (extra,)
    with pytest.raises(ValueError, match="not in a collection"):
        MoveInCollection(project.environment, 0)
    check_undo_redo(project, Reconnect(m["line"], "b", extra))
    check_undo_redo(project, Reconnect(m["line"], "A", extra))
    with pytest.raises(ValueError, match="A or B"):
        Reconnect(m["line"], "C", extra)


def test_replace_strategy_copies_shared_values() -> None:
    m = build()
    project = m["project"]
    project.seabed.stiffness = 3.0e5
    project.seabed.friction = 0.6
    command = ReplaceStrategy(project, "seabed", BathymetrySeabed)
    check_undo_redo(project, command)
    CommandStack().push(command)
    assert isinstance(project.seabed, BathymetrySeabed)
    assert project.seabed.stiffness == 3.0e5 and project.seabed.friction == 0.6
    chain = m["chain"]
    chain.axial.damping = 3.0
    check_undo_redo(project, ReplaceStrategy(chain, "axial", ViscoelasticAxial()))
    plain = ReplaceStrategy(chain, "axial", ViscoelasticAxial(), copy_common=False)
    assert plain.new.damping == 0.0
    with pytest.raises(TypeError, match="not a sub-object"):
        ReplaceStrategy(chain, "diameter", FlatSeabed)


def test_split_section_and_library_item() -> None:
    m = build()
    project, s1, line = m["project"], m["s1"], m["line"]
    command = SplitSection(s1, 50.0)
    check_undo_redo(project, command)
    CommandStack().push(command)
    assert [s.length for s in line.sections] == [50.0, 150.0, 210.0]
    assert [s.segments for s in line.sections] == [5, 15, 21]
    assert command.new_section.line_type is m["chain"]
    with pytest.raises(ValueError, match="not inside"):
        SplitSection(s1, 60.0)
    with pytest.raises(ValueError, match="does not belong"):
        SplitSection(Section(length=10.0, segments=4), 5.0)
    one = line.sections.append(Section(line_type=m["chain"], length=10.0, segments=1))
    with pytest.raises(ValueError, match="one element"):
        SplitSection(one, 5.0)
    chain = project.line_types.append(ChainType("R4"))
    item = ApplyLibraryItem(chain, {"diameter": 0.23, "mass_per_length": 290.0, "grade": "R4"})
    check_undo_redo(project, item)
    assert chain in item.affected()


def test_macro_command_rolls_back_on_failure() -> None:
    m = build()
    project, s1 = m["project"], m["s1"]

    class Boom(Command):
        def do(self) -> None:
            raise RuntimeError("boom")

    before = project.to_dict()
    macro = MacroCommand("bad", [SetProperty(s1, "length", 1.0), Boom()])
    with pytest.raises(RuntimeError):
        CommandStack().push(macro)
    assert project.to_dict() == before
    with pytest.raises(NotImplementedError):
        Command().do()
    with pytest.raises(NotImplementedError):
        Command().undo()
    assert Command().affected() == () and not Command().merge_with(Command())
    assert repr(Command()) == "<Command 'Edit'>"


# --------------------------------------------------------------------------- removal


def test_remove_refuses_while_referenced() -> None:
    m = build()
    project, fairlead = m["project"], m["fairlead"]
    with pytest.raises(ObjectInUseError) as raised:
        CommandStack().push(RemoveObject(project, fairlead))
    referrers_found = {(obj.name, name) for obj, name in raised.value.referrers}
    assert ("mooring", "end_a") in referrers_found and ("second", "end_a") in referrers_found
    assert len(project.points) == 2
    with pytest.raises(ValueError, match="not in a collection"):
        plan_removal(project, [project.environment])
    lonely = project.points.append(FixedPoint("lonely"))
    check_undo_redo(project, RemoveObject(project, lonely))


def test_cascading_removal_of_a_point() -> None:
    m = build()
    project, fairlead = m["project"], m["fairlead"]
    command = RemoveObject(project, fairlead, cascade=True)
    check_undo_redo(project, command)
    CommandStack().push(command)
    assert len(project.lines) == 0
    assert len(project.failures) == 0 and len(project.controls) == 0
    assert len(project.syrope_history) == 0 and len(project.outputs.channels) == 0
    assert m["group"].members == (m["anchor"],)
    assert fairlead in command.affected()


def test_cascading_removal_of_a_line_type_removes_emptied_lines() -> None:
    m = build()
    project, wire = m["project"], m["wire"]
    stack = CommandStack()
    before = project.to_dict()
    stack.push(RemoveObject(project, wire, cascade=True))
    assert [line.name for line in project.lines] == ["mooring"]
    assert [s.name for s in m["line"].sections] == ["s1"]
    assert project.failures[0].lines == (m["line"],)
    assert len(project.syrope_history) == 0
    stack.undo()
    assert project.to_dict() == before
    assert m["second"].parent is project and m["s2"].line_type is wire


def test_removal_of_a_section_and_a_rod() -> None:
    m = build()
    project = m["project"]
    check_undo_redo(project, RemoveObject(project, [m["s1"], m["s2"]], cascade=True))
    rod_type = project.rod_types.append(RodType("rt"))
    rod = project.rods.append(Rod("r", rod_type=rod_type))
    tip = project.points.append(RodPoint("tip", rod_end=rod.end_b))
    tether = project.lines.append(Line("tether", end_a=rod.end_a, end_b=m["anchor"]))
    tether.sections.append(Section(line_type=m["chain"], length=10.0, segments=2))
    assert {obj.name for obj, _ in referrers(project, rod)} == {"tip", "tether"}
    with pytest.raises(ObjectInUseError):
        plan_removal(project, [rod])
    plan = plan_removal(project, [rod], cascade=True)
    assert isinstance(plan, RemovalPlan)
    assert {obj.name for obj in plan.removed} == {"r", "tip", "tether"}
    check_undo_redo(project, RemoveObject(project, rod, cascade=True))
    del tip


def test_removal_of_a_body_with_optional_references() -> None:
    project = Project()
    body = project.bodies.append(Rigid6Body("b"))
    rod_type = project.rod_types.append(RodType("rt"))
    rod = project.rods.append(Rod("r", rod_type=rod_type, attachment="free", body=body))
    point = project.points.append(BodyPoint("p", body=body))
    plan = plan_removal(project, [body], cascade=True)
    assert plan.removed == [body, point]
    assert plan.reference_updates == [(rod, "body", body, None)]
    check_undo_redo(project, RemoveObject(project, body, cascade=True))
    assert isinstance(body, Body)


def test_removing_nested_targets_once() -> None:
    m = build()
    project = m["project"]
    plan = plan_removal(project, [m["line"], m["s1"], m["line"]], cascade=True)
    assert plan.removed[0] is m["line"]
    assert m["s1"] not in plan.removed
    command = RemoveObject(project, [])
    assert command.text == "Remove nothing"
    assert RemoveObject(project, [m["s1"], m["s2"]]).text == "Remove 2"
    end = m["line"].end_connections.append(EndConnection())
    check_undo_redo(project, RemoveObject(project, end))


# --------------------------------------------------------------------------- stack


def test_stack_state_and_notifications() -> None:
    m = build()
    s1 = m["s1"]
    stack = CommandStack()
    seen: list[StackChanged] = []
    cancel = stack.subscribe(seen.append)
    assert not stack.can_undo and not stack.can_redo
    assert stack.undo_text == "" and stack.redo_text == ""
    stack.undo()
    stack.redo()
    assert stack.is_clean
    stack.push(SetProperty(s1, "length", 1.0))
    assert not stack.is_clean and stack.can_undo and "length" in stack.undo_text
    stack.push(SetProperty(s1, "length", 2.0, merge=True))
    stack.set_clean()
    stack.push(SetProperty(s1, "length", 3.0, merge=True))
    assert len(stack) == 3  # no merge across the clean point
    stack.undo()
    assert stack.is_clean and stack.redo_text
    stack.undo()
    stack.undo()
    stack.push(SetProperty(s1, "segments", 3))
    assert not stack.is_clean and len(stack) == 1 and stack.index == 1
    assert [type(c) for c in stack.commands()] == [SetProperty]
    stack.clear()
    assert len(stack) == 0 and stack.is_clean
    cancel()
    cancel()
    stack.push(SetProperty(s1, "segments", 4))
    actions = [event.action for event in seen]
    assert actions[:4] == ["push", "push", "clean", "push"] and actions[-1] == "clear"
    assert "undo" in actions


def test_stack_limit() -> None:
    m = build()
    s1 = m["s1"]
    stack = CommandStack(limit=2)
    stack.set_clean()
    for value in (1.0, 2.0, 3.0):
        stack.push(SetProperty(s1, "length", value))
    assert len(stack) == 2 and not stack.is_clean
    stack.undo()
    stack.undo()
    assert s1.length == 1.0 and not stack.can_undo
    stack = CommandStack(limit=1)
    stack.push(SetProperty(s1, "length", 4.0))
    stack.set_clean()
    stack.push(SetProperty(s1, "length", 5.0))
    assert not stack.is_clean


def test_macro_context_groups_and_rolls_back() -> None:
    m = build()
    project, s1, s2 = m["project"], m["s1"], m["s2"]
    stack = CommandStack()
    before = project.to_dict()
    with stack.macro("Lengthen") as group:
        stack.push(SetProperty(s1, "length", 300.0))
        stack.push(SetProperty(s2, "length", 300.0))
        assert not stack.can_undo
    assert len(stack) == 1 and stack.undo_text == "Lengthen" and len(group.commands) == 2
    assert set(group.affected()) == {s1, s2}
    stack.undo()
    assert project.to_dict() == before
    with pytest.raises(RuntimeError), stack.macro("fails"):
        stack.push(SetProperty(s1, "length", 1.0))
        raise RuntimeError("stop")
    assert s1.length == 200.0 and len(stack) == 1
    with stack.macro("empty"):
        pass
    assert len(stack) == 1


def test_import_deck_command() -> None:
    project = Project()
    source = DeckReader.read(EXAMPLES / "spread_3line_chain.dat")
    counts = (len(source.lines), len(source.points), len(source.outputs.channels))
    environment = source.environment
    command = ImportDeck(project, source)
    before = project.to_dict()
    stack = CommandStack()
    stack.push(command)
    assert (len(project.lines), len(project.points), len(project.outputs.channels)) == counts
    assert project.environment is environment
    stack.undo()
    assert project.to_dict() == before
    stack.redo()
    assert project.lines[0].end_a in project.points
    keep = Project()
    keep_env = keep.environment
    stack.push(
        ImportDeck(
            keep, DeckReader.read(EXAMPLES / "spread_3line_chain.dat"), replace_environment=False
        )
    )
    assert keep.environment is keep_env


def test_random_command_sequences_undo_to_the_start() -> None:
    rng = random.Random(7)
    m = build()
    project = m["project"]
    stack = CommandStack()
    start = project.to_dict()
    for _ in range(200):
        lines = list(project.lines)
        sections = [s for line in lines for s in line.sections]
        choice = rng.randrange(5)
        if choice == 0 and sections:
            stack.push(SetProperty(rng.choice(sections), "length", rng.uniform(5.0, 300.0)))
        elif choice == 1:
            stack.push(AddObject(project.points, FixedPoint(position=(rng.random(), 0, -10))))
        elif choice == 2 and len(project.points) > 2:
            stack.push(RemoveObject(project, project.points[-1], cascade=True))
        elif choice == 3 and lines:
            stack.push(MoveInCollection(rng.choice(lines), rng.randrange(len(lines))))
        elif choice == 4 and lines:
            target = rng.choice(list(project.points))
            line = rng.choice(lines)
            if target is not line.end_b:
                stack.push(Reconnect(line, "A", target))
        if rng.random() < 0.2:
            stack.undo()
    while stack.can_undo:
        stack.undo()
    assert project.to_dict() == start
