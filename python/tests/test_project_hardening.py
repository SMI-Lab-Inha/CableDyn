# SPDX-License-Identifier: Apache-2.0
"""Regressions for the hardening review of cabledyn.project (one test per finding class)."""

from __future__ import annotations

import json
import math
import os
import zipfile
from pathlib import Path
from typing import Any

import pytest

from cabledyn.builder import DeckModel
from cabledyn.errors import DeckFormatError
from cabledyn.project import (
    REGISTRY,
    AddObject,
    Attachment,
    BendStiffener,
    BuoyancyModule,
    ChainType,
    Command,
    CommandStack,
    DeckExportError,
    DeckReader,
    DeckWriter,
    ExternalLoad,
    ExtraOption,
    FixedPoint,
    Group,
    ImportDeck,
    Line,
    LinearAxial,
    ModelObject,
    MultiTrainSea,
    Project,
    ProjectFormatError,
    ProjectStore,
    RemoveObject,
    Rod,
    RodEndpoint,
    RodType,
    Section,
    SetProperty,
    Severity,
    TurbinePoint,
    UnitSystem,
    WaveTrain,
    model_type,
    validate_project,
)
from cabledyn.project import store as store_module
from cabledyn.project.settings import AnalysisSettings
from cabledyn.project.units import AXIAL_DAMPING, DIMENSIONLESS
from cabledyn.project.validation import errors

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples"


def _errors(project: Project) -> list[str]:
    return [str(issue) for issue in errors(validate_project(project))]


# --------------------------------------------------------------------------- P1 ids


def test_external_loads_are_numbered_in_order_after_edits() -> None:
    first = DeckReader.read(EXAMPLES / "als_volturnus_line_break_tension.dat")
    first.external_loads.append(ExternalLoad("second", body=first.bodies[0], force=(1.0, 0.0, 0.0)))
    deck = EXAMPLES / "als_volturnus_line_break_tension.dat"
    project = DeckReader.from_text(DeckWriter(first).to_text(), path=deck)
    body = project.bodies[0]
    assert len(project.external_loads) == 2 and _errors(project) == []
    stack = CommandStack()
    stack.push(RemoveObject(project, project.external_loads[0]))
    assert _errors(project) == []
    stack.push(AddObject(project.external_loads, ExternalLoad("third", body=body)))
    model, id_map = DeckWriter(project).to_model()
    assert [row.id for row in model.external_loads] == [1, 2]
    assert [id_map.deck_id(load) for load in project.external_loads] == [1, 2]
    assert _errors(project) == []


# --------------------------------------------------------------------------- P2


def test_a_kept_option_never_overrides_the_typed_property() -> None:
    project = DeckReader.read(EXAMPLES / "chain_catenary_r3_100m.dat")
    project.settings.time_step = 0.1
    project.settings.duration = 1.0
    project.settings.extra_options.append(ExtraOption(keyword="dtM", values=("0.5",)))
    issues = validate_project(project)
    warning = [i for i in issues if i.severity is Severity.WARNING and i.prop == "keyword"]
    assert warning and "typed property wins" in warning[0].message
    assert errors(issues) == []
    again = DeckReader.from_text(DeckWriter(project).to_text())
    assert again.settings.time_step == 0.1
    project.settings.extra_options.append(ExtraOption(keyword="MyOption", values=("1",)))
    assert ExtraOption(keyword="wavetrain").validate() == []
    assert ExtraOption(keyword="vesselMotion").validate()


def test_removing_a_body_removes_the_rods_that_need_it() -> None:
    project = DeckReader.read(EXAMPLES / "mixed_body_rods_points.dat")
    body_rods = [rod for rod in project.rods if rod.attachment in {"body", "body_pinned"}]
    assert body_rods
    body = body_rods[0].body
    assert body is not None
    before = project.to_dict()
    stack = CommandStack()
    stack.push(RemoveObject(project, body, cascade=True))
    assert all(rod.body is not None for rod in project.rods if rod.attachment.startswith("body"))
    assert all(issue.layer != 2 for issue in validate_project(project, deck=False))
    stack.undo()
    assert project.to_dict() == before
    free = Rod(attachment="free")
    assert free.required_references() == frozenset({"rod_type"})
    assert Rod(attachment="body").required_references() == frozenset({"rod_type", "body"})


def test_group_membership_never_blocks_a_removal() -> None:
    project = Project()
    anchor = project.points.append(FixedPoint("anchor"))
    group = project.studio.groups.append(Group("moorings", members=[anchor]))
    assert project.used_by(anchor) == []
    stack = CommandStack()
    stack.push(RemoveObject(project, anchor))
    assert group.members == ()
    stack.undo()
    assert group.members == (anchor,)
    describe = Group.get_property("members").describe()
    assert describe["weak"] is True


def test_uids_stay_unique() -> None:
    project = DeckReader.read(EXAMPLES / "spread_3line_chain.dat")
    line = project.lines[0]
    clone = Line.from_dict(line.to_dict(), within=project)
    assert clone.uid != line.uid and clone.end_a is line.end_a
    assert clone.sections[0].line_type is line.sections[0].line_type
    project.lines.append(clone)
    DeckWriter(project).to_model()
    same = Line.from_dict(line.to_dict(), fresh_uids=False)
    with pytest.raises(ValueError, match="already used"):
        project.lines.append(same)
    assert same.parent is None and project.find(line.uid) is line
    environment = project.environment
    clash = type(environment)()
    clash._uid = line.uid
    with pytest.raises(ValueError, match="already used"):
        project.environment = clash
    assert project.environment is environment and project.find(environment.uid) is environment
    same_uid = type(environment).from_dict(environment.to_dict(), fresh_uids=False)
    project.environment = same_uid  # replacing an object by a copy of itself is fine
    assert project.find(same_uid.uid) is same_uid
    project.points[0]._uid = project.points[1].uid
    assert project.duplicate_uids() == [project.points[1].uid]
    document = ProjectStore.to_json(project)
    points = document["project"]["properties"]["points"]
    points[1]["uid"] = points[0]["uid"]
    with pytest.raises(ProjectFormatError, match="more than one object"):
        ProjectStore.from_json(document)


def test_ownership_cycles_are_refused() -> None:
    sea = MultiTrainSea()
    train = sea.trains.append(WaveTrain())
    with pytest.raises(ValueError, match="cannot own itself"):
        train.wave = sea
    inner = WaveTrain()
    nested = MultiTrainSea()
    inner.wave = nested
    nested.trains.append(WaveTrain())
    with pytest.raises(ValueError, match="cannot own itself"):
        nested.trains.append(inner)
    assert sea.root() is sea and inner.root() is inner


def test_axial_damping_keeps_its_units_apart() -> None:
    axial = LinearAxial()
    axial.set_deck_damping(-0.8)
    assert (axial.damping_mode, axial.damping_ratio) == ("ratio", 0.8)
    assert axial.deck_damping() == -0.8
    assert LinearAxial.get_property("damping_ratio").dimension is DIMENSIONLESS
    assert LinearAxial.get_property("damping").dimension is AXIAL_DAMPING
    engineering = UnitSystem.engineering()
    assert engineering.to_display(axial.damping_ratio, DIMENSIONLESS) == 0.8
    axial.set_deck_damping(-0.0)
    assert math.copysign(1.0, axial.deck_damping()) < 0.0  # type: ignore[arg-type]
    axial.set_deck_damping((3.0e4, 2.0e5))
    assert axial.deck_damping() == (3.0e4, 2.0e5) and axial.damping_mode == "damping"
    model = DeckModel.new()
    model.add_line_type("odd", diam=0.1, mass=1.0, ea=1.0, ba=(1.0, 2.0, 3.0))
    with pytest.raises(DeckFormatError, match="unsupported BA"):
        DeckReader.from_model(model)


def test_the_project_deck_is_runnable_from_the_project_folder(tmp_path: Path) -> None:
    deck = EXAMPLES / "torsion_lazy_wave_hangoff_twist.dat"
    project = DeckReader.read(deck)
    folder = ProjectStore.save(project, tmp_path / "sub" / "twist.cdproj", folder=True)
    model = DeckModel.load(folder / "model.dat")
    motion = model.options.value("motionFile")
    assert motion is not None and (folder / motion).is_file()
    document = json.loads((folder / "project.json").read_text(encoding="utf-8"))
    assert document["deck"]["missing_files"] == []
    stored = document["project"]["properties"]["deck_path"]
    same_drive = os.path.splitdrive(str(tmp_path))[0] == os.path.splitdrive(str(deck))[0]
    assert Path(stored).is_absolute() != same_drive
    loaded = ProjectStore.load(folder)
    assert Path(loaded.deck_path or "") == deck.resolve()
    archive = ProjectStore.save(project, tmp_path / "zip" / "twist.cdproj")
    with zipfile.ZipFile(archive) as opened:
        text = opened.read("model.dat").decode("utf-8")
    zipped = DeckModel.from_text(text, path=archive.parent / "model.dat")
    zipped_motion = zipped.options.value("motionFile")
    assert zipped_motion is not None and (archive.parent / zipped_motion).is_file()
    project.motion.set_value("file", "no_such_motion.txt")
    target = ProjectStore.save(project, tmp_path / "missing.cdproj", folder=True)
    document = json.loads((target / "project.json").read_text(encoding="utf-8"))
    assert document["deck"]["missing_files"]
    fresh = Project("new")
    fresh.points.append(FixedPoint("p"))
    ProjectStore.save(fresh, tmp_path / "fresh.cdproj", folder=True)


def test_project_file_side_paths(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from types import SimpleNamespace

    from cabledyn.project.deck import DeckWriter as Writer

    work = tmp_path / "work"
    work.mkdir()
    (work / "deck.dat").write_text(
        (EXAMPLES / "spread_3line_chain.dat").read_text(encoding="utf-8"), encoding="utf-8"
    )
    project = DeckReader.read(work / "deck.dat")
    target = ProjectStore.save(project, work / "p.cdproj")
    document = json.loads(zipfile.ZipFile(target).read("project.json"))
    assert document["project"]["properties"]["deck_path"] == "deck.dat"
    moved = tmp_path / "moved"
    moved.mkdir()
    (moved / "p.cdproj").write_bytes(target.read_bytes())
    assert Path(ProjectStore.load(moved / "p.cdproj").deck_path or "") == moved / "deck.dat"
    syrope = DeckReader.read(EXAMPLES / "syrope_polyester_mooring.dat")
    folder = ProjectStore.save(syrope, tmp_path / "syrope.cdproj", folder=True)
    text = (folder / "model.dat").read_text(encoding="utf-8")
    assert "SYROPE:" in text
    record = SimpleNamespace(keyword="WaterKin", values=("SeaState",))
    numeric = SimpleNamespace(keyword="WaterKin", values=("3",))
    path = SimpleNamespace(keyword="WaterKin", values=("kin.dat",))
    off = SimpleNamespace(keyword="motionFile", values=("0",))
    fake = SimpleNamespace(line_types=[], options=[record, numeric, path, off])
    assert store_module._side_files(fake) == ["kin.dat"]  # type: ignore[arg-type]

    def fail(self: Any, target: Any) -> None:
        raise OSError("no access")

    monkeypatch.setattr(store_module._deck.DeckFile, "_rebase_relative_inputs", fail)
    folder = ProjectStore.save(project, tmp_path / "fails.cdproj", folder=True)
    document = json.loads((folder / "project.json").read_text(encoding="utf-8"))
    assert document["deck"] == {"written": False, "reason": "no access"}
    assert Writer is DeckWriter


def test_project_files_are_strict_json(tmp_path: Path) -> None:
    project = Project()
    with pytest.raises(ValueError, match="JSON data"):
        project.studio.layouts = {"k": {1, 2}}
    with pytest.raises(ValueError, match="JSON data"):
        project.studio.layouts = {"k": float("nan")}
    with pytest.raises(ValueError, match="finite"):
        project.settings.time_step = float("nan")
    project.studio.layouts = {"k": [1, 2]}
    target = ProjectStore.save(project, tmp_path / "p.cdproj", folder=True)
    text = (target / "project.json").read_text(encoding="utf-8")
    (target / "project.json").write_text(text.replace('"layouts": {', '"bad": NaN, "layouts": {'))
    with pytest.raises(ProjectFormatError, match="non-standard number NaN"):
        ProjectStore.load(target)


def test_model_dat_needs_every_validation_layer(tmp_path: Path) -> None:
    project = DeckReader.read(EXAMPLES / "spread_3line_chain.dat")
    line = project.lines[0]
    line.ancillaries.append(BendStiffener(start=0.0, length=5.0))
    line.ancillaries.append(BendStiffener(start=2.0, length=5.0))
    target = ProjectStore.save(project, tmp_path / "p.cdproj", folder=True)
    document = json.loads((target / "project.json").read_text(encoding="utf-8"))
    assert document["deck"]["written"] is False and "overlaps" in document["deck"]["reason"]
    assert not (target / "model.dat").exists()


# --------------------------------------------------------------------------- P3


def test_unknown_objects_keep_their_place_and_references(tmp_path: Path) -> None:
    project = DeckReader.read(EXAMPLES / "spread_3line_chain.dat")
    document = ProjectStore.to_json(project)
    points = document["project"]["properties"]["points"]
    retyped = points[1]
    retyped["type"] = "plugin.special_point"
    loaded = ProjectStore.from_json(document)
    assert loaded.lines[0].end_b is None or loaded.lines[0].end_b.uid != retyped["uid"]
    again = ProjectStore.to_json(loaded)
    assert again["project"]["properties"]["points"][1] == retyped
    assert again["unknown_objects"] == []
    line_dicts = again["project"]["properties"]["lines"]
    ends = [item["properties"][end] for item in line_dicts for end in ("end_a", "end_b")]
    assert retyped["uid"] in ends
    document["unknown_objects"] = [{"owner": "gone", "slot": "points", "index": 0, "data": {}}]
    kept = ProjectStore.from_json(document)
    assert ProjectStore.to_json(kept)["unknown_objects"][0]["owner"] == "gone"


def test_a_listener_cannot_push_during_a_change() -> None:
    project = Project()
    point = project.points.append(FixedPoint("p"))
    stack = CommandStack()
    failures: list[Exception] = []

    def listener(event: Any) -> None:
        try:
            stack.push(SetProperty(point, "volume", 2.0))
        except RuntimeError as exc:
            failures.append(exc)

    project.events.subscribe(listener)
    stack.push(SetProperty(point, "mass", 1.0))
    assert failures and point.volume == 0.0
    assert [type(c).__name__ for c in stack.commands()] == ["SetProperty"]


def test_failed_undo_and_redo_keep_the_stack_consistent() -> None:
    class Fragile(Command):
        def __init__(self) -> None:
            self.fail_undo = False
            self.fail_do = False

        def do(self) -> None:
            if self.fail_do:
                raise RuntimeError("do")

        def undo(self) -> None:
            if self.fail_undo:
                raise RuntimeError("undo")

    stack = CommandStack()
    command = Fragile()
    stack.push(command)
    command.fail_undo = True
    with pytest.raises(RuntimeError):
        stack.undo()
    assert stack.index == 1 and stack.can_undo
    command.fail_undo = False
    stack.undo()
    command.fail_do = True
    with pytest.raises(RuntimeError):
        stack.redo()
    assert stack.index == 0 and stack.can_redo


def test_malformed_hints_are_ignored(tmp_path: Path) -> None:
    project = DeckReader.read(EXAMPLES / "chain_catenary_r3_100m.dat")
    options = project.deck_hints["options"]
    first = next(iter(options))
    options[first] = {"keyword": "g", "values": ["9.8"]}
    options["motion"] = "garbage"
    project.settings.extra_options.append(ExtraOption(keyword="WriteLog", values=("0",)))
    project.settings.extra_options[-1].deck_hints["position"] = True
    assert _errors(project) == []
    ProjectStore.save(project, tmp_path / "p.cdproj")
    train_sea = MultiTrainSea()
    train = train_sea.trains.append(WaveTrain())
    train.deck_hints["option"] = {"position": "x"}
    project.environment.waves = train_sea
    assert "wavetrain" in DeckWriter(project).to_text(validate=False)


def test_loader_errors_are_project_format_errors(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    folder = tmp_path / "bad.cdproj"
    folder.mkdir()
    (folder / "project.json").write_bytes(b"\xff\xfe{")
    with pytest.raises(ProjectFormatError, match="UTF-8"):
        ProjectStore.load(folder)
    (folder / "project.json").write_text("[" * 100000 + "]" * 100000, encoding="utf-8")
    with pytest.raises(ProjectFormatError, match=r"nested too deeply|not a CableDyn"):
        ProjectStore.load(folder)
    good = ProjectStore.save(Project(), tmp_path / "good.cdproj")
    monkeypatch.setattr(store_module, "_MAX_DOCUMENT_BYTES", 10)
    with pytest.raises(ProjectFormatError, match="too large"):
        ProjectStore.load(good)
    good_folder = ProjectStore.save(Project(), tmp_path / "good_folder.cdproj", folder=True)
    with pytest.raises(ProjectFormatError, match="too large"):
        ProjectStore.load(good_folder)


def test_an_interrupted_save_keeps_the_previous_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = ProjectStore.save(Project("first"), tmp_path / "p.cdproj")
    before = target.read_bytes()

    def broken(source: Any, destination: Any) -> None:
        raise OSError("disk full")

    monkeypatch.setattr(os, "replace", broken)
    with pytest.raises(OSError, match="disk full"):
        ProjectStore.save(Project("second"), target, overwrite=True)
    assert target.read_bytes() == before
    assert [p.name for p in tmp_path.iterdir()] == ["p.cdproj"]


def test_unregistered_subclasses_are_not_written_as_their_parent() -> None:
    class MyChain(ChainType):
        """Not registered."""

    with pytest.raises(TypeError, match="not registered"):
        MyChain().to_dict()

    class Other(ModelObject):
        """Collides with a registered key."""

    with pytest.raises(ValueError, match="already used"):
        model_type("line")(Other)
    assert Other.type_key == "object"
    assert "line" in REGISTRY


def test_abstract_flag_is_per_class() -> None:
    assert Attachment.is_abstract() and not BuoyancyModule.is_abstract()
    assert not FixedPoint.is_abstract()


def test_dynamic_solver_with_spectral_radius() -> None:
    project = DeckReader.read(EXAMPLES / "chain_catenary_r3_100m.dat")
    model, _ = DeckWriter(project).to_model()
    model.options.set("dynamic_solver", "1e-8", "1e-14", "30", "12", "0.4")
    read = DeckReader.from_model(model)
    assert read.settings.solver_spectral_radius == 0.4
    assert "dynamic_solver 1e-8 1e-14 30 12 0.4" in DeckWriter(read).to_text(validate=False)


def test_describe_reports_numeric_limits() -> None:
    from cabledyn.project import Point

    assert Point.get_property("position").dimension.key == "length"
    assert AnalysisSettings.get_property("spectral_radius").describe()["maximum"] == 1.0
    assert Rod.get_property("segments").describe()["minimum"] == 0
    assert Section.get_property("segments").describe()["limit"] == "segments"


def test_rod_endpoints_cannot_be_replaced() -> None:
    rod = Rod()
    with pytest.raises(AttributeError):
        SetProperty(rod, "end_b", RodEndpoint.make("B"))
    with pytest.raises(AttributeError):
        rod.end_b = RodEndpoint.make("B")


def test_import_deck_moves_objects_only_when_applied() -> None:
    project = Project()
    source = DeckReader.read(EXAMPLES / "spread_3line_chain.dat")
    lines = len(source.lines)
    command = ImportDeck(project, source)
    assert len(source.lines) == lines and command.affected() == (project,)
    stack = CommandStack()
    stack.push(command)
    assert project.deck_hints["options"] and len(project.lines) == lines
    text = DeckWriter(project).to_text()
    assert text == DeckModel.load(EXAMPLES / "spread_3line_chain.dat").to_text().replace(
        DeckModel.load(EXAMPLES / "spread_3line_chain.dat").title, project.title
    )
    stack.undo()
    assert "options" not in project.deck_hints and len(project.lines) == 0
    stack.redo()
    assert len(project.lines) == lines


def test_turbine_points_follow_their_turbine() -> None:
    deck = ROOT / "validation" / "bodies" / "cases" / "V-S1m" / "cabledyn.dat"
    project = DeckReader.read(deck)
    points = [p for p in project.points if isinstance(p, TurbinePoint)]
    assert points and all(p.turbine is not None for p in points)
    turbine = points[0].turbine
    assert turbine is not None
    turbine.number = 7
    assert "Turbine7" in DeckWriter(project).to_text(validate=False) or "T7" in DeckWriter(
        project
    ).to_text(validate=False)
    stack = CommandStack()
    stack.push(RemoveObject(project, turbine, cascade=True))
    assert all(p.turbine is not turbine for p in project.objects_of(TurbinePoint))
    farm = DeckModel.new()
    farm.add_point(1, "Turbine3", 0.0, 0.0, -10.0)
    read = DeckReader.from_model(farm)
    assert read.points[0].turbine_number == 3  # type: ignore[attr-defined]
    loose = TurbinePoint(turbine_number=3)
    assert loose.number() == 3 and loose.validate() == []
    assert TurbinePoint().validate()
    project.points.append(TurbinePoint("orphan", turbine=type(turbine)(number=9)))
    with pytest.raises(DeckExportError, match="outside the project"):
        DeckWriter(project).to_text(validate=False)
    project.points[-1].turbine = None  # type: ignore[attr-defined]
    with pytest.raises(DeckExportError, match="has no turbine"):
        DeckWriter(project).to_text(validate=False)


def test_rod_types_must_belong_to_the_project() -> None:
    project = DeckReader.read(EXAMPLES / "mixed_body_rods_points.dat")
    project.rods[0].rod_type = RodType(project.rod_types[0].name)
    with pytest.raises(DeckExportError, match="outside the project"):
        DeckWriter(project).to_text(validate=False)


# --------------------------------------------------------------------------- oracles


def _strip_hints(project: Project) -> None:
    for obj in project.walk():
        obj.deck_hints.clear()


def _semantic(project: Project) -> Any:
    """The project's graph without uids (references become positions), names and hints."""
    graph = Project.from_dict(project.to_dict())
    _strip_hints(graph)
    data = graph.to_dict()
    order: dict[str, int] = {obj.uid: index for index, obj in enumerate(graph.walk())}

    def clean(node: Any) -> Any:
        if isinstance(node, dict):
            return {
                key: clean(value)
                for key, value in node.items()
                if key not in {"uid", "name", "deck_path"}
            }
        if isinstance(node, list):
            return [clean(value) for value in node]
        if isinstance(node, str) and node in order:
            return order[node]
        return node

    return clean(data["properties"])


@pytest.mark.parametrize(
    "name",
    [
        "spread_3line_chain.dat",
        "lazy_wave_buoyancy_modules.dat",
        "mixed_body_rods_points.dat",
        "chain_two_train_sea.dat",
        "lazy_wave_vessel_rao.dat",
        "torsion_lazy_wave_hangoff_twist.dat",
        "als_volturnus_line_break_tension.dat",
        "syrope_polyester_mooring.dat",
    ],
)
def test_hint_free_write_keeps_the_semantics(name: str) -> None:
    project = DeckReader.read(EXAMPLES / name)
    reference = _semantic(project)
    _strip_hints(project)
    text = DeckWriter(project).to_text()
    stripped = DeckReader.from_text(text, path=EXAMPLES / name)
    assert _semantic(stripped) == reference


HAND_DECK = """--------------------- CableDyn Input File ------------------------------------
hand-written check
--------------------- LINE TYPES ---------------------------------------------
TypeName Diam MassDenInAir EA BA/-zeta EI Cd_n Cd_t Ca_n Ca_t
(-) (m) (kg/m) (N) (N-s/-) (N-m^2) (-) (-) (-) (-)
chain 0.126 105.6 7.518e8 -0.8 0.0 1.37 0.64 1.0 0.0
--------------------- POINTS -------------------------------------------------
ID Type X Y Z Mass Vol CdA Ca
(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)
7 Vessel 1.0 2.0 -5.0 0 0 0 0
3 Fixed 100.0 0.0 -30.0 0 0 0 0
--------------------- LINES --------------------------------------------------
ID NodeA NodeB Outputs
(-) (-) (-) (-)
4 7 3 pt
--------------------- SECTIONS -----------------------------------------------
LineID LineType Length NumSegs
(-) (-) (m) (-)
4 chain 60.0 6
4 chain 50.0 5
--------------------- OPTIONS ------------------------------------------------
30.0 WtrDpth
1025 rhoW
0.01 dtM
1.0 TMax
uniform 0.5 0.0 0.0 current
--------------------- OUTPUTS ------------------------------------------------
"FairTen4"
--------------------- need this line -----------------------------------------
"""


def test_a_hand_written_deck_reads_to_the_expected_objects() -> None:
    project = DeckReader.from_text(HAND_DECK)
    chain = project.line_types[0]
    assert (chain.name, chain.diameter, chain.mass_per_length) == ("chain", 0.126, 105.6)
    assert isinstance(chain.axial, LinearAxial) and chain.axial.stiffness == 7.518e8
    assert (chain.axial.damping_mode, chain.axial.damping_ratio) == ("ratio", 0.8)
    assert (chain.drag_normal, chain.drag_axial, chain.added_mass_normal) == (1.37, 0.64, 1.0)
    fairlead, anchor = project.points
    assert type(fairlead).__name__ == "FloaterPoint" and fairlead.kind == "floater"  # type: ignore[attr-defined]
    assert fairlead.position == (1.0, 2.0, -5.0) and anchor.position == (100.0, 0.0, -30.0)
    line = project.lines[0]
    assert (line.end_a, line.end_b, line.outputs) == (fairlead, anchor, "pt")
    assert [(s.length, s.segments) for s in line.sections] == [(60.0, 6), (50.0, 5)]
    assert line.length == 110.0
    assert project.environment.water_depth == 30.0 and project.environment.water_density == 1025
    assert project.environment.current.velocity == (0.5, 0.0, 0.0)  # type: ignore[attr-defined]
    channel = project.outputs.channels[0]
    assert channel.target is line  # type: ignore[attr-defined]
    _strip_hints(project)
    text = DeckWriter(project).to_text()
    assert "-0.8" in text and "FairTen1" in text and "Vessel" in text


@pytest.mark.parametrize("seed", range(20))
def test_generated_graphs_survive_write_and_read(seed: int) -> None:
    """Seeded random spreads of chain lines read back to the same object graph."""
    import random

    from cabledyn.project import FloaterPoint, GenericLineType

    rng = random.Random(seed)
    project = Project(f"generated {seed}")
    project.environment.water_depth = 100.0
    types = [
        project.line_types.append(
            GenericLineType(
                f"type {index}",
                diameter=rng.uniform(0.05, 0.2),
                mass_per_length=rng.uniform(20.0, 200.0),
                axial=LinearAxial(stiffness=rng.uniform(1e8, 1e9)),
            )
        )
        for index in range(rng.randint(1, 3))
    ]
    for index in range(rng.randint(1, 4)):
        angle = rng.uniform(0.0, 2.0 * math.pi)
        fairlead = project.points.append(
            FloaterPoint(position=(10 * math.cos(angle), 10 * math.sin(angle), -10.0))
        )
        reach = rng.uniform(200.0, 400.0)
        anchor = project.points.append(
            FixedPoint(position=(reach * math.cos(angle), reach * math.sin(angle), -100.0))
        )
        line = project.lines.append(Line(f"line {index}", end_a=fairlead, end_b=anchor))
        for _ in range(rng.randint(1, 3)):
            line.sections.append(
                Section(
                    line_type=rng.choice(types),
                    length=rng.uniform(80.0, 200.0),
                    segments=rng.randint(4, 20),
                )
            )
    text = DeckWriter(project).to_text()
    again = DeckReader.from_text(text)
    assert _semantic(again) == _semantic(project)
    assert DeckWriter(again).to_text() == text
