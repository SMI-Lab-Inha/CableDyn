# SPDX-License-Identifier: Apache-2.0
"""Object model of cabledyn.project: descriptors, objects, events, invariants, serialisation."""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import pytest

from cabledyn.project import (
    REGISTRY,
    AnalysisSettings,
    BathymetrySeabed,
    BendingModel,
    BendRestrictor,
    BendStiffener,
    BodyPoint,
    BuoyancyModule,
    CableType,
    ChainType,
    ChildAdded,
    ChildMoved,
    ChildRemoved,
    ClumpWeight,
    EndConnection,
    Environment,
    Event,
    FibreRopeType,
    FixedPoint,
    Floater,
    FloaterPoint,
    FloaterType,
    GenericLineType,
    Group,
    HydroDatabaseFile,
    JonswapWave,
    Line,
    LinearAxial,
    ModelObject,
    MultiTrainSea,
    NoHydroDatabase,
    NoWaves,
    ObjectChannel,
    Project,
    PropertyChanged,
    Quantity,
    ReferenceChanged,
    RegularWave,
    Rigid6Body,
    Rod,
    RodType,
    Role,
    Section,
    Severity,
    StudioData,
    SyropeAxial,
    TouchdownProtection,
    Turbine,
    UnitSystem,
    ViscoelasticAxial,
    WaveTrain,
    WireRopeType,
    model_type,
)
from cabledyn.project.base import SerialContext, TypeRegistry, UnknownTypeError, collection_of
from cabledyn.project.descriptors import (
    Bool,
    Choice,
    Colour,
    Integer,
    JsonValue,
    OptionalBool,
    OptionalColour,
    OptionalInteger,
    OptionalQuantity,
    OptionalText,
    OptionalVec3,
    Text,
    TextList,
    Vec3,
    Vector,
)
from cabledyn.project.issues import Issue
from cabledyn.project.schema import DECK_LIMITS, DeckLimit, deck_limit
from cabledyn.project.units import LENGTH


class Widget(ModelObject):
    """A test class using every scalar descriptor."""

    type_label = "Widget"
    length = Quantity(LENGTH, 2.0, minimum=0.0, maximum=10.0, limit="coordinate")
    width = OptionalQuantity(LENGTH, minimum=1.0, maximum=3.0)
    count = Integer(3, minimum=1, maximum=5)
    spare = OptionalInteger(minimum=0)
    flag = Bool(False)
    maybe = OptionalBool()
    kind = Choice(("a", "b"))
    label_text = Text("x")
    note = OptionalText()
    colour = Colour("#112233")
    tint = OptionalColour()
    offset = Vec3(LENGTH)
    target = OptionalVec3(LENGTH, limit="coordinate")
    values = Vector(LENGTH, (1.0,), lengths=(1, 3))
    keys = TextList(("k",))
    blob = JsonValue({"a": 1})


@pytest.fixture(autouse=True, scope="module")
def _registered_widget() -> Iterator[None]:
    """Register the test class for this module only (no leak into other tests)."""
    model_type("test.widget")(Widget)
    yield
    REGISTRY.unregister("test.widget")


def _events(obj: ModelObject) -> list[Event]:
    seen: list[Event] = []
    obj.events.subscribe(seen.append)
    return seen


def chain_project() -> tuple[Project, GenericLineType, FloaterPoint, FixedPoint, Line]:
    project = Project("demo")
    chain = project.line_types.append(GenericLineType("chain", diameter=0.2))
    fairlead = project.points.append(FloaterPoint("fairlead", position=(0, 0, -14)))
    anchor = project.points.append(FixedPoint("anchor", position=(400, 0, -50)))
    line = project.lines.append(Line("mooring", end_a=fairlead, end_b=anchor))
    line.sections.append(Section("s1", line_type=chain, length=410.0, segments=41))
    return project, chain, fairlead, anchor, line


# --------------------------------------------------------------------------- descriptors


def test_descriptor_defaults_and_coercion() -> None:
    w = Widget("w")
    assert (w.length, w.width, w.count, w.flag, w.kind) == (2.0, None, 3, False, "a")
    assert w.offset == (0.0, 0.0, 0.0) and w.values == (1.0,) and w.keys == ("k",)
    w.length = 3
    assert isinstance(w.length, float)
    w.count = 4.0  # an integral float is accepted
    assert w.count == 4
    w.offset = [1, 2, 3]  # type: ignore[assignment]
    assert w.offset == (1.0, 2.0, 3.0)
    w.colour = "#AABBCC"
    assert w.colour == "#aabbcc"
    w.values = (1.0, 2.0, 3.0)
    w.target = (1, 1, 1)  # type: ignore[assignment]
    assert w.target == (1.0, 1.0, 1.0)
    w.target = None
    w.blob = {"b": [1, 2]}
    assert w.blob == {"b": [1, 2]}


@pytest.mark.parametrize(
    ("name", "value", "error"),
    [
        ("length", "3", TypeError),
        ("length", True, TypeError),
        ("count", 2.5, TypeError),
        ("count", True, TypeError),
        ("flag", 1, TypeError),
        ("maybe", "yes", TypeError),
        ("kind", "c", ValueError),
        ("kind", 1, TypeError),
        ("label_text", 1, TypeError),
        ("note", 1, TypeError),
        ("colour", "red", ValueError),
        ("tint", "#12", ValueError),
        ("offset", (1, 2), ValueError),
        ("offset", "abc", TypeError),
        ("target", (1, 2), ValueError),
        ("values", (1.0, 2.0), ValueError),
        ("keys", "abc", TypeError),
        ("keys", (1,), TypeError),
        ("width", "x", TypeError),
        ("spare", "1", TypeError),
    ],
)
def test_descriptor_rejects_bad_values(name: str, value: object, error: type) -> None:
    w = Widget()
    with pytest.raises(error):
        setattr(w, name, value)


def test_layer_one_checks() -> None:
    w = Widget()
    assert w.validate() == []
    w.length = 20.0
    w.width = 0.5
    w.count = 9
    w.spare = -1
    w.target = (2.0e6, 0.0, 0.0)
    messages = {(issue.prop, issue.message) for issue in w.validate()}
    assert ("length", "must be at most 10") in messages
    assert ("width", "must be at least 1") in messages
    assert ("count", "must be at most 5") in messages
    assert ("spare", "must be at least 0") in messages
    assert ("target", "magnitude must be at most 1e+06") in messages
    for bad in (float("nan"), float("inf")):
        with pytest.raises(ValueError, match="finite"):
            w.length = bad
        with pytest.raises(ValueError, match="finite"):
            w.values = (bad,)
    w.length = -1.0
    assert any(i.message == "must be at least 0" for i in w.validate())
    issue = w.validate()[0]
    assert issue.layer == 1 and issue.is_error
    assert str(issue).startswith("error: Widget: length:")


def test_shared_deck_limits() -> None:
    assert deck_limit("segments").violation(0) == "value must be at least 1"
    assert deck_limit("section_length").violation(0.0) == "value must be greater than 0"
    assert deck_limit("line_mass").violation(-2.0e6) == "magnitude must be at most 1e+06"
    assert deck_limit("wave_period").violation(5.0) is None
    assert DeckLimit().violation(1.0e30) is None
    assert set(DECK_LIMITS) >= {"coordinate", "gravity", "line_diameter"}


def test_property_introspection() -> None:
    props = {prop.name: prop for prop in Widget.properties()}
    assert list(props)[:3] == ["name", "description", "tags"]
    length = props["length"]
    info = length.describe()
    assert info["dimension"] == "length" and info["maximum"] == 10.0 and info["kind"] == "quantity"
    assert props["kind"].describe()["choices"] == ["a", "b"]
    assert props["values"].describe()["lengths"] == [1, 3]
    assert Widget.get_property("count").label == "Count"
    with pytest.raises(KeyError, match="no property"):
        Widget.get_property("nothing")
    assert repr(length) == "<Quantity 'length'>"
    assert Widget.__dict__["length"].role is Role.PHYSICS
    assert props["colour"].role is Role.APPEARANCE
    line_props = {prop.name: prop.describe() for prop in Line.properties()}
    assert line_props["end_a"]["target"] == "LineEndTarget"
    assert line_props["sections"]["min_items"] == 1
    assert {p.name: p.describe() for p in Project.properties()}["seabed"]["base"] == "Seabed"


def test_strategy_options_list_concrete_classes() -> None:
    options = ChainType.get_property("axial").options()  # type: ignore[attr-defined]
    assert {cls.__name__ for cls in options} >= {"LinearAxial", "ViscoelasticAxial", "SyropeAxial"}
    hydro = FloaterType.get_property("hydrodynamics").options()  # type: ignore[attr-defined]
    assert NoHydroDatabase in hydro and HydroDatabaseFile in hydro


def test_children_property_cannot_be_replaced() -> None:
    project = Project()
    with pytest.raises(AttributeError):
        project.lines = None  # type: ignore[assignment]
    with pytest.raises(AttributeError):
        Project.get_property("lines").coerce([])
    with pytest.raises(AttributeError, match="no settable property"):
        Project(lines=[])
    with pytest.raises(AttributeError, match="no settable property"):
        Widget(unknown=1)


def test_read_only_property() -> None:
    rod = Rod()
    with pytest.raises(AttributeError, match="read-only"):
        rod.end_a.end = "B"
    assert rod.endpoint("a") is rod.end_a and rod.endpoint("B") is rod.end_b
    with pytest.raises(ValueError, match="A or B"):
        rod.endpoint("C")
    assert rod.end_a.label() == "Rod End A"
    assert rod.end_a.parent is rod
    detached = rod.end_a.__class__.make("B")
    assert detached.label() == "End B"


# --------------------------------------------------------------------------- objects


def test_object_identity_and_tree() -> None:
    project, chain, fairlead, anchor, line = chain_project()
    assert len(project.uid) == 32 and project.uid != chain.uid
    assert chain.parent is project and chain.slot == "line_types"
    assert line.sections[0].root() is project
    assert project.find(line.sections[0].uid) is line.sections[0]
    assert project.find("missing") is None
    assert list(project.objects_of(Section)) == [line.sections[0]]
    assert chain in list(project.children())
    assert line.references() == [("end_a", fairlead), ("end_b", anchor)]
    assert project.used_by(chain) == [(line.sections[0], "line_type")]
    assert project.used_by(fairlead) == [(line, "end_a")]
    assert line.length == 410.0 and line.segment_count == 41
    assert project.get_value("title") == "CableDyn model"
    project.set_value("title", "renamed")
    assert project.title == "renamed"
    assert repr(chain).startswith("<GenericLineType 'chain'")
    assert chain.label() == "Line type 'chain'"
    assert Section().label() == "Section"
    assert collection_of(chain) is project.line_types
    assert collection_of(project.environment) is None
    assert collection_of(project) is None


def test_events_bubble_to_the_project() -> None:
    project, chain, fairlead, _anchor, line = chain_project()
    seen = _events(project)
    chain.diameter = 0.3
    chain.diameter = 0.3  # unchanged: no event
    line.end_b = fairlead
    extra = project.points.append(FixedPoint("extra"))
    project.points.move(extra, 0)
    project.points.move(extra, 0)  # no move
    project.points.remove(extra)
    kinds = [type(event) for event in seen]
    assert kinds == [PropertyChanged, ReferenceChanged, ChildAdded, ChildMoved, ChildRemoved]
    changed = seen[0]
    assert isinstance(changed, PropertyChanged)
    assert (changed.source, changed.name, changed.old, changed.new) == (chain, "diameter", 0.2, 0.3)
    moved = seen[3]
    assert isinstance(moved, ChildMoved) and (moved.old_index, moved.new_index) == (2, 0)
    assert extra.parent is None


def test_subscription_cancel() -> None:
    w = Widget()
    seen: list[Event] = []
    sub = w.events.subscribe(seen.append)
    assert len(w.events) == 1
    w.count = 2
    sub.cancel()
    sub.cancel()
    w.count = 4
    assert len(seen) == 1 and len(w.events) == 0


def test_collection_api() -> None:
    project = Project()
    a = project.points.append(FixedPoint("A"))
    b = FixedPoint("b")
    project.points.insert(-5, b)
    assert project.points[0] is b and project.points[0:2] == [b, a]
    assert project.points.find("a") is a and project.points.find("zz") is None
    assert project.points.by_uid(a.uid) is a and project.points.by_uid("x") is None
    assert a in project.points and 1 not in project.points
    assert bool(project.points) and not project.lines
    assert project.points.owner is project and project.points.name == "points"
    assert project.points.item_type.__name__ == "Point"
    assert "A" in repr(project.points)
    with pytest.raises(TypeError, match="holds Point"):
        project.points.append(Line())  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="already belongs"):
        project.points.append(a)
    with pytest.raises(ValueError, match="is not in"):
        project.lines.index(Line())
    project.points.extend([FixedPoint("c")])
    project.points.clear()
    assert len(project.points) == 0 and a.parent is None


def test_child_slot_replacement() -> None:
    project = Project()
    old = project.environment
    new = Environment()
    project.environment = new
    assert new.parent is project and old.parent is None
    with pytest.raises(ValueError, match="already belongs"):
        Project().environment = new
    with pytest.raises(TypeError):
        project.environment = Project()  # type: ignore[assignment]
    project.seabed = BathymetrySeabed(file="b.txt")
    assert project.seabed.parent is project


def test_line_types_keep_physics_construction_and_appearance_apart() -> None:
    chain = ChainType("R4 chain", grade="R4", nominal_diameter=0.12)
    roles = {prop.name: prop.role for prop in ChainType.properties()}
    assert roles["diameter"] is Role.PHYSICS
    assert roles["grade"] is Role.CONSTRUCTION
    assert roles["appearance"] is Role.APPEARANCE
    assert type(chain.appearance).__name__ == "ChainAppearance"
    assert type(WireRopeType().appearance).__name__ == "WireRopeAppearance"
    assert type(FibreRopeType().appearance).__name__ == "FibreRopeAppearance"
    assert type(CableType().appearance).__name__ == "CableAppearance"
    assert type(GenericLineType().appearance).__name__ == "PlainLineAppearance"
    names = [prop.name for prop in ChainType.properties()]
    assert names.index("appearance") < names.index("grade")


# --------------------------------------------------------------------------- invariants


def _messages(obj: ModelObject, *, recursive: bool = False) -> list[str]:
    return [issue.message for issue in obj.validate(recursive=recursive)]


def test_reference_invariants() -> None:
    line = Line()
    messages = _messages(line)
    assert messages.count("refers to nothing") == 2
    assert "needs at least 1 item(s)" in messages
    project, _chain, fairlead, _anchor, real = chain_project()
    outsider = FixedPoint("elsewhere")
    real.end_b = outsider
    assert any("outside this model" in m for m in _messages(real))
    real.end_b = fairlead
    assert "both ends attach to one object" in _messages(real)
    group = Group(members=[outsider])
    project.studio.groups.append(group)
    assert any("outside this model" in m for m in _messages(group))
    from cabledyn.project import SyropeIC

    assert "needs at least 1 item(s)" in _messages(SyropeIC())


def test_line_invariants() -> None:
    _project, _chain, _fairlead, _anchor, line = chain_project()
    line.end_connections.append(EndConnection(end="A"))
    line.end_connections.append(EndConnection(end="A"))
    module = BuoyancyModule(arc_length=500.0)
    line.attachments.append(module)
    series = ClumpWeight(arc_length=10.0, pitch=10.0, last_arc_length=40.0)
    line.attachments.append(series)
    assert series.is_series() and series.locations() == [10.0, 20.0, 30.0, 40.0]
    assert not module.is_series() and module.locations() == [500.0]
    line.ancillaries.append(BendStiffener(start=0.0, length=5.0))
    line.ancillaries.append(TouchdownProtection(start=3.0, end=8.0))
    line.ancillaries.append(BendRestrictor(start=20.0, end=10.0))
    line.ancillaries.append(TouchdownProtection(start=400.0, end=500.0))
    messages = _messages(line, recursive=True)
    assert "End A has two end connections" in messages
    assert any("beyond the line" in m for m in messages)
    assert "overlaps another ancillary" in messages
    assert "ends before it starts" in messages
    assert "extends beyond the line" in messages
    half = ClumpWeight(pitch=5.0)
    assert "give pitch and last_arc_length together" in _messages(half)
    assert half.locations() == [0.0]


def test_module_volume_warning() -> None:
    module = BuoyancyModule(volume=1.0)
    assert module.geometric_volume() is None
    module.geometry.outer_diameter = 1.0
    module.geometry.length = 1.0
    issues = module.validate()
    assert issues and issues[0].severity is Severity.WARNING
    module.geometry.length = 1.27324
    assert module.validate() == []


def test_end_connection_invariants() -> None:
    row = EndConnection(direction=(0, 0, 0))
    assert "must not be zero" in _messages(row)
    row = EndConnection(normal=(1, 0, 0))
    assert "normal and pretwist need a torsion" in _messages(row)
    row = EndConnection(torsion="free")
    assert "a torsion needs a normal" in _messages(row)
    row.normal = (1.0, 0.0, 0.0)
    assert _messages(row) == []


def test_type_invariants() -> None:
    assert _messages(ViscoelasticAxial()) == [
        "give either dynamic_stiffness or both alpha_mbl and beta"
    ]
    assert _messages(ViscoelasticAxial(dynamic_stiffness=2e8)) == []
    assert _messages(ViscoelasticAxial(alpha_mbl=10.0, beta=0.2)) == []
    assert _messages(ViscoelasticAxial(alpha_mbl=10.0))
    assert _messages(SyropeAxial()) == ["needs a settings file"]
    assert _messages(BendingModel(shear_stiffness=1.0))
    full = BendingModel(
        shear_stiffness=1.0,
        torsional_stiffness=1.0,
        rotary_inertia_axial=1.0,
        rotary_inertia_normal=1.0,
    )
    assert full.extended() and _messages(full) == []
    assert _messages(RodType(axial_drag=1.0)) == ["give both axial coefficients or neither"]
    assert _messages(BathymetrySeabed()) == ["needs a bathymetry file"]
    assert _messages(LinearAxial(stiffness=1e16)) == ["value must be at most 1e+15"]


def test_body_and_rod_invariants() -> None:
    body = Rigid6Body(drag_area=(1.0, 2.0, 3.0))
    messages = _messages(body)
    assert "the CableDyn layout takes one value" in messages
    body.inertia = (1.0,)
    assert "the CableDyn layout takes 0 or 3" in _messages(body)
    moordyn = Rigid6Body(row_format="moordyn")
    assert "the MoorDyn layout needs it" in _messages(moordyn)
    rod = Rod(attachment="body")
    assert any("needs its body" in m for m in _messages(rod))


def test_settings_and_wave_invariants() -> None:
    settings = AnalysisSettings(solver_relative_tolerance=1e-8)
    assert "give the four dynamic-solver values together" in _messages(settings)
    settings = AnalysisSettings(solver_spectral_radius=0.4)
    assert "the solver spectral radius needs the four solver values" in _messages(settings)
    assert _messages(WaveTrain(wave=NoWaves())) == ["a train needs a wave or spectrum"]
    assert _messages(WaveTrain(wave=RegularWave(theory="stream"))) == [
        "a regular train must be airy"
    ]
    assert _messages(WaveTrain(wave=JonswapWave())) == []
    assert "needs at least 1 item(s)" in _messages(MultiTrainSea())


def test_project_invariants() -> None:
    project, _chain, *_ = chain_project()
    project.line_types.append(GenericLineType("CHAIN"))
    project.line_types.append(GenericLineType(""))
    project.rod_types.append(RodType("r"))
    project.turbines.append(Turbine(number=1))
    project.turbines.append(Turbine(number=1))
    project.unknown_objects.append({"data": {}})
    messages = _messages(project)
    assert "line type name 'CHAIN' is used twice" in messages
    assert "a line type needs a name" in messages
    assert "turbine number is used twice" in messages
    assert any("unknown type" in m for m in messages)


# --------------------------------------------------------------------------- serialisation


def test_to_dict_and_from_dict_round_trip() -> None:
    project, chain, _fairlead, _anchor, line = chain_project()
    floater_type = project.floater_types.append(FloaterType("semi", draft=20.0))
    floater_type.hydrodynamics = HydroDatabaseFile(file="semi.hdb", source_format="BEM")
    floater = project.bodies.append(Floater("FOWT", floater_type=floater_type))
    project.points.append(BodyPoint("hang-off", body=floater))
    project.outputs.channels.append(ObjectChannel(target=line, quantity="FairTen"))
    project.studio.groups.append(Group("g", members=[line, chain]))
    project.environment.waves = JonswapWave()
    chain.deck_hints["note"] = ("kept", 1)
    data = project.to_dict()
    copy = Project.from_dict(data, fresh_uids=False)
    assert copy.to_dict() == data
    assert copy.uid == project.uid
    assert copy.lines[0].end_a is copy.points[0]
    assert copy.studio.groups[0].members == (copy.lines[0], copy.line_types[0])
    assert copy.bodies[0].floater_type is copy.floater_types[0]  # type: ignore[attr-defined]
    assert copy.line_types[0].deck_hints["note"] == ["kept", 1]
    fresh = Project.from_dict(data, fresh_uids=True)
    assert fresh.uid != project.uid
    with pytest.raises(TypeError, match="does not describe"):
        Line.from_dict(data)


def test_from_dict_tolerates_unknown_and_bad_entries() -> None:
    project, _chain, *_ = chain_project()
    data = project.to_dict()
    data["properties"]["points"].append({"type": "plugin.buoy", "uid": "u1", "properties": {}})
    data["properties"]["points"].append("garbage")
    data["properties"]["lines"][0]["properties"]["outputs"] = 5
    data["properties"]["lines"][0]["properties"]["bogus"] = 1
    data["properties"]["lines"][0]["properties"]["end_b"] = "no-such-uid"
    data["properties"]["line_types"].append(
        {"type": "section", "uid": "s", "properties": {"length": 1.0}}
    )
    data["properties"]["environment"] = {"type": "plugin.env", "properties": {}}
    data["properties"]["outputs"]["properties"] = "bad"
    ctx = SerialContext()
    rebuilt = ctx.build(data)
    ctx.resolve()
    assert isinstance(rebuilt, Project)
    assert len(rebuilt.points) == 2
    assert rebuilt.lines[0].outputs == "-"
    assert rebuilt.lines[0].end_b is None
    assert isinstance(rebuilt.environment, Environment)
    assert [entry["data"]["type"] for entry in ctx.unknown] == ["plugin.buoy"]
    written = rebuilt.to_dict()["properties"]
    assert written["environment"] == {"type": "plugin.env", "properties": {}}
    assert written["lines"][0]["properties"]["bogus"] == 1
    assert written["lines"][0]["properties"]["outputs"] == 5
    assert written["lines"][0]["properties"]["end_b"] == "no-such-uid"
    rebuilt.lines[0].outputs = "pt"
    assert rebuilt.to_dict()["properties"]["lines"][0]["properties"]["outputs"] == "pt"
    text = " | ".join(ctx.warnings)
    for fragment in ("unknown type", "malformed", "kept the default", "unknown property"):
        assert fragment in text
    assert "reference end_b is unresolved" in text
    assert "not a valid item" in text
    with pytest.raises(KeyError):
        ctx.lookup(5)


def test_registry() -> None:
    registry = TypeRegistry()
    registry.register("x", Widget)
    registry.register("x", Widget)
    with pytest.raises(ValueError, match="already used"):
        registry.register("x", Project)
    assert registry.get("x") is Widget and "x" in registry and registry.keys() == ["x"]
    with pytest.raises(UnknownTypeError):
        registry.get("y")
    assert REGISTRY.get("line") is Line
    assert "test.widget" in REGISTRY


def test_studio_data_and_units() -> None:
    studio = StudioData()
    assert studio.unit_system() == UnitSystem.si()
    studio.set_unit_system(UnitSystem.engineering())
    assert studio.unit_system() == UnitSystem.engineering()


def test_project_validate_all_runs_every_layer() -> None:
    project, *_ = chain_project()
    project.environment.water_depth = 50.0
    assert project.validate_all() == []
    project.lines[0].sections[0].segments = 0
    issues = project.validate_all()
    assert issues and all(issue.layer == 1 for issue in issues)
    assert isinstance(issues[0], Issue)


def test_floater_type_extension_slot(request: pytest.FixtureRequest) -> None:
    class TabulatedDatabase(HydroDatabaseFile):
        """A plugin database kind."""

    model_type("test.hydro.tabulated")(TabulatedDatabase)
    request.addfinalizer(lambda: REGISTRY.unregister("test.hydro.tabulated"))

    floater_type = FloaterType()
    assert isinstance(floater_type.hydrodynamics, NoHydroDatabase)
    floater_type.hydrodynamics = TabulatedDatabase(file="x")
    data: dict[str, Any] = floater_type.to_dict()
    assert data["properties"]["hydrodynamics"]["type"] == "test.hydro.tabulated"
    assert isinstance(FloaterType.from_dict(data).hydrodynamics, TabulatedDatabase)
    assert TabulatedDatabase in FloaterType.get_property("hydrodynamics").options()  # type: ignore[attr-defined]


def test_descriptor_edge_cases() -> None:
    from cabledyn.project.descriptors import Children, Ref, RefList, _same

    with pytest.raises(ValueError, match="at least one choice"):
        Choice(())
    with pytest.raises(ValueError, match="is not one of"):
        Choice(("a",), "b")
    assert OptionalVec3(LENGTH, (1, 2, 3)).default == (1.0, 2.0, 3.0)
    assert OptionalQuantity(LENGTH).dimension is LENGTH
    assert OptionalVec3(LENGTH).dimension is LENGTH
    assert Vector(LENGTH).dimension is LENGTH
    ref = Line.get_property("end_a")
    assert isinstance(ref, Ref)
    with pytest.raises(TypeError, match="must refer to a LineEndTarget"):
        ref.coerce(Section())
    refs = Group.get_property("members")
    assert isinstance(refs, RefList)
    with pytest.raises(TypeError, match="sequence of objects"):
        refs.coerce("abc")
    from cabledyn.project import SyropeIC

    with pytest.raises(TypeError, match="Line objects"):
        SyropeIC.get_property("lines").coerce([Section()])
    assert refs.describe()["min_items"] == 0
    assert refs.from_json([], SerialContext(), Group()) == ()
    children = Project.get_property("lines")
    assert isinstance(children, Children) and children.is_owner
    project = Project()
    assert children.same(project.lines, project.lines)
    a, b = Section(), Section()
    assert _same(a, a) and not _same(a, b)
    assert _same(float("nan"), float("nan"))
    assert not _same((1.0,), (1.0, 2.0))


def test_misc_edges() -> None:
    from cabledyn.project import RemoveObject

    project, _chain, _fairlead, _anchor, line = chain_project()
    command = RemoveObject(project, line)
    assert command.affected() == (line,)
    hub_owner = Widget()
    first = hub_owner.events.subscribe(lambda event: None)
    second = hub_owner.events.subscribe(lambda event: None)
    second.cancel()
    first.cancel()
    assert len(hub_owner.events) == 0
    assert line.references()[0][0] == "end_a"
    line.end_a = None
    assert [name for name, _ in line.references()] == ["end_b"]
