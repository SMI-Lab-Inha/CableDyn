# SPDX-License-Identifier: Apache-2.0
"""Deck adapters: round trips over every repository deck, and the writer's mapping rules."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from cabledyn import builder as _b
from cabledyn.builder import DeckModel
from cabledyn.errors import DeckFormatError
from cabledyn.project import (
    AxialModel,
    BathymetrySeabed,
    BendingModel,
    BodyPoint,
    BuoyancyModule,
    ChainType,
    ClumpWeight,
    CommandStack,
    ConnectPoint,
    Control,
    CurrentModel,
    DeckExportError,
    DeckReader,
    DeckWriter,
    EndConnection,
    EquivalentBuoyancy,
    ExternalLoad,
    ExtraOption,
    Failure,
    FixedPoint,
    Floater,
    FloaterMotionRecord,
    FloaterPoint,
    FloaterRAO,
    FreePoint,
    GenericLineType,
    JonswapWave,
    Line,
    LinearAxial,
    MotionFile,
    MultiTrainSea,
    NamedChannel,
    NoCurrent,
    NoMotion,
    NoWaves,
    ObjectChannel,
    OchiHubbleWave,
    Point,
    Point3Body,
    ProfileCurrent,
    Project,
    RegularWave,
    RemoveObject,
    Rigid6Body,
    Rod,
    RodPoint,
    RodType,
    Section,
    SpectrumWave,
    SyropeAxial,
    SyropeIC,
    Turbine,
    TurbinePoint,
    UniformCurrent,
    ViscoelasticAxial,
    WaveModel,
    WaveTrain,
    model_type,
    validate_project,
)
from cabledyn.project.deck import (
    _decode_current,
    _decode_scalar,
    _decode_train,
    _decode_vector,
    _decode_wave,
    _encode_current,
    _encode_wave,
)
from cabledyn.project.validation import errors

ROOT = Path(__file__).resolve().parents[2]


@model_type("test.wave.odd")
class OddWave(WaveModel):
    """A wave model without a deck form."""


@model_type("test.current.odd")
class OddCurrent(CurrentModel):
    """A current model without a deck form."""


DECK_FOLDERS = ("examples", "tests/data", "validation")


def _candidates() -> list[Path]:
    found: list[Path] = []
    for folder in DECK_FOLDERS:
        found.extend(sorted((ROOT / folder).rglob("*.dat")))
    return found


def _load(path: Path) -> DeckModel | None:
    for caller_driven in (False, True):
        try:
            return DeckModel.load(path, caller_driven=caller_driven)
        except DeckFormatError:
            continue
    return None


CANDIDATES = _candidates()


def normalise(project: Project) -> Any:
    """Return the project's dictionary with uids replaced by their traversal order."""
    data = project.to_dict()
    uids: dict[str, str] = {}

    def collect(node: Any) -> None:
        if isinstance(node, dict):
            if "uid" in node and "type" in node:
                uids[node["uid"]] = f"#{len(uids)}"
            for value in node.values():
                collect(value)
        elif isinstance(node, list):
            for value in node:
                collect(value)

    def substitute(node: Any) -> Any:
        if isinstance(node, dict):
            return {key: substitute(value) for key, value in node.items()}
        if isinstance(node, list):
            return [substitute(value) for value in node]
        return uids.get(node, node) if isinstance(node, str) else node

    collect(data)
    return substitute(data)


def test_the_round_trip_covers_the_repository_decks() -> None:
    accepted = [path for path in CANDIDATES if _load(path) is not None]
    rejected = sorted(path.name for path in CANDIDATES if path not in accepted)
    # Not decks (Syrope tables, a harmonics table, a plot table) and two MoorDyn
    # driver decks without the dtM/TMax that a standalone CableDyn run requires.
    assert rejected == [
        "cable_200m.dat",
        "lozon_coupled_harmonics.dat",
        "marker_curvature_plot.dat",
        "mooring_80m.dat",
        "syrope_owc.dat",
        "syrope_settings.dat",
    ]
    assert len(accepted) >= 99


@pytest.mark.parametrize("path", CANDIDATES, ids=lambda p: p.relative_to(ROOT).as_posix())
def test_deck_round_trip(path: Path, tmp_path: Path) -> None:
    model = _load(path)
    if model is None:
        pytest.skip("not a deck the deck model accepts")
    canonical = model.to_text()
    project = DeckReader.from_model(model)
    text = DeckWriter(project).to_text()
    assert text == canonical
    again = DeckReader.from_text(text, path=model.path, caller_driven=model.caller_driven)
    assert normalise(again) == normalise(project)
    assert DeckWriter(again).to_text() == canonical
    assert errors(validate_project(project)) == []


# --------------------------------------------------------------------------- reader details


def test_reader_maps_deck_vocabulary_to_classes() -> None:
    project = DeckReader.read(ROOT / "examples" / "mixed_body_rods_points.dat")
    kinds = {type(point).__name__ for point in project.points}
    assert {"BodyPoint", "FixedPoint"} <= kinds
    assert project.bodies and project.rods
    assert project.title
    assert Path(project.deck_path or "").name == "mixed_body_rods_points.dat"
    lazy = DeckReader.read(ROOT / "examples" / "lazy_wave_buoyancy_modules.dat")
    attachments = [a for line in lazy.lines for a in line.attachments]
    assert attachments and all(isinstance(a, BuoyancyModule) for a in attachments)
    rao = DeckReader.read(ROOT / "examples" / "lazy_wave_vessel_rao.dat")
    assert isinstance(rao.motion, FloaterRAO) and rao.motion.reference is not None
    record = DeckReader.read(ROOT / "examples" / "lazy_wave_vessel_motion.dat")
    assert isinstance(record.motion, (FloaterMotionRecord, MotionFile))
    trains = DeckReader.read(ROOT / "examples" / "chain_two_train_sea.dat")
    assert isinstance(trains.environment.waves, MultiTrainSea)
    assert len(trains.environment.waves.trains) == 2
    clump = DeckReader.read(ROOT / "examples" / "clump_weight_free_point.dat")
    assert any(isinstance(p, FreePoint) for p in clump.points)


def test_reader_keeps_unknown_and_shadowed_options_verbatim() -> None:
    original = DeckModel.load(ROOT / "examples" / "chain_catenary_r3_100m.dat")
    original.options.add("dtM", 0.5)
    original.options.add("dtM", 0.25, description="the effective step")
    original.options.add("motionFile", "0")
    original.options.add("modified_newton", "maybe")
    original.options.add("waves", "airy", "1", "8")
    project = DeckReader.from_model(original)
    extras = [(row.keyword, tuple(row.values)) for row in project.settings.extra_options]
    assert ("dtM", ("0.5",)) in extras
    assert ("motionFile", ("0",)) in extras
    assert ("modified_newton", ("maybe",)) in extras
    assert ("waves", ("airy", "1", "8")) in extras
    assert project.settings.time_step == 0.25
    assert isinstance(project.motion, NoMotion)
    assert DeckWriter(project).to_text(validate=False) == original.to_text(validate=False)
    source = DeckModel.load(ROOT / "examples" / "chain_catenary_r3_100m.dat")
    source.options.add("motionFile", "a.txt")
    source.options.add("vesselMotion", "b.txt")
    source.options.add("wavetrain", "airy", "1", "8", "0")
    source.options.add("wavetrain", "airy", "x", "8", "0")
    project = DeckReader.from_model(source)
    assert isinstance(project.motion, FloaterMotionRecord)
    assert isinstance(project.environment.waves, MultiTrainSea)
    assert len(project.environment.waves.trains) == 1
    assert DeckWriter(project).to_text(validate=False) == source.to_text(validate=False)


def test_codecs() -> None:
    assert _decode_scalar("friction", ("none",)) == 0.0
    assert _decode_scalar("bool", ("yes",)) is True
    assert _decode_scalar("bool", ("off",)) is False
    with pytest.raises(ValueError, match="logical"):
        _decode_scalar("bool", ("maybe",))
    with pytest.raises(ValueError, match="one value"):
        _decode_scalar("float", ("1", "2"))
    assert _decode_scalar("text", ("SeaState",)) == "SeaState"
    assert isinstance(_decode_wave(("none",), False), NoWaves)
    for values, cls in (
        (("airy", "2", "8", "0"), RegularWave),
        (("jonswap", "2", "8", "3.3", "0"), JonswapWave),
        (("PiersonMoskowitz", "2", "8", "0"), SpectrumWave),
        (("ochi-hubble", "2", "12", "3", "1", "7", "1", "0"), OchiHubbleWave),
    ):
        wave = _decode_wave(values, False)
        assert isinstance(wave, cls)
        assert _decode_wave(_encode_wave(wave), False).to_dict()["properties"] == {
            **wave.to_dict()["properties"]
        }
    with pytest.raises(ValueError, match="unrecognised"):
        _decode_wave(("airy", "1"), False)
    with pytest.raises(ValueError, match="only airy"):
        _decode_wave(("stream", "1", "8", "0"), True)
    train = _decode_train(("jonswap", "3", "7", "3.3", "0", "6"))
    assert train.spreading == 6.0
    assert _decode_train(("airy", "1", "8", "0")).spreading is None
    current = _decode_current(("profile", "-30", "0.2", "0", "0", "0", "0.5", "0", "0"))
    assert isinstance(current, ProfileCurrent)
    assert _encode_current(current)[0] == "profile"
    assert isinstance(_decode_current(("uniform", "1", "0", "0")), UniformCurrent)
    assert isinstance(_decode_current(("none",)), NoCurrent)
    with pytest.raises(ValueError, match="unrecognised"):
        _decode_current(("uniform", "1"))
    assert _decode_vector(("1|2|3",)) == (1.0, 2.0, 3.0)
    with pytest.raises(ValueError, match=r"x\|y\|z"):
        _decode_vector(("1|2",))
    with pytest.raises(DeckExportError):
        _encode_wave(MultiTrainSea())
    with pytest.raises(DeckExportError):
        _encode_current(Project())  # type: ignore[arg-type]


# --------------------------------------------------------------------------- writer


def full_project() -> Project:
    """A project using most object kinds, written by hand (no hints)."""
    p = Project("hand built", title="hand built")
    chain = p.line_types.append(
        GenericLineType("chain", diameter=0.1, mass_per_length=50.0, drag_normal=1.2)
    )
    chain.axial = LinearAxial(stiffness=5.0e8, damping=(-1.0,))
    cable = p.line_types.append(ChainType("cable", diameter=0.2, mass_per_length=80.0))
    cable.bending = BendingModel(bending_stiffness=1.0e4)
    rod_type = p.rod_types.append(RodType("pipe", diameter=1.0, mass_per_length=100.0))
    body = p.bodies.append(
        Rigid6Body(
            "buoy", position=(50.0, 0.0, -20.0), mass=1.0e4, volume=20.0, inertia=(1e5, 1e5, 1e5)
        )
    )
    rod = p.rods.append(
        Rod(
            "spar",
            rod_type=rod_type,
            attachment="pinned",
            position_a=(0, 50, -5),
            position_b=(0, 50, -25),
            segments=4,
        )
    )
    fairlead = p.points.append(FloaterPoint("fairlead", position=(0.0, 0.0, -10.0)))
    anchor = p.points.append(FixedPoint("anchor", position=(300.0, 0.0, -100.0)))
    on_body = p.points.append(BodyPoint("pad", body=body, position=(0.0, 0.0, -1.0)))
    p.points.append(RodPoint("tip", rod_end=rod.end_b))
    line = p.lines.append(Line("main", end_a=fairlead, end_b=anchor, outputs="pt"))
    line.sections.append(Section(line_type=chain, length=150.0, segments=15))
    line.sections.append(Section(line_type=chain, length=180.0, segments=18))
    riser = p.lines.append(Line("riser", end_a=on_body, end_b=anchor))
    riser.sections.append(Section(line_type=chain, length=280.0, segments=28))
    tether = p.lines.append(Line("tether", end_a=rod.end_a, end_b=anchor))
    tether.sections.append(Section(line_type=chain, length=320.0, segments=32))
    p.environment.water_depth = 100.0
    p.environment.gravity = 9.81
    p.settings.time_step = 0.01
    p.settings.duration = 1.0
    p.outputs.channels.append(ObjectChannel(target=line, quantity="FairTen"))
    p.outputs.channels.append(ObjectChannel(target=on_body, quantity="Point", qualifier="pz"))
    p.outputs.channels.append(NamedChannel(channel="AnchTen1"))
    del cable
    return p


def test_writer_builds_a_valid_deck_from_scratch() -> None:
    project = full_project()
    writer = DeckWriter(project)
    text = writer.to_text()
    assert "FairTen1" in text and "Point3pz" in text
    assert "R1A" in text and "Rod1B" in text and "Body1" in text

    assert "9.81" in text and "WtrDpth" in text
    id_map = writer.id_map
    model, _ = DeckWriter(project).to_model()
    assert id_map.deck_id(project.lines[1]) == 2
    assert id_map.object_for("line_type", "CHAIN") is project.line_types[0]
    assert id_map.object_for("point", 4) is project.points[3]
    assert id_map.deck_id(Line()) is None
    row = text.splitlines().index(next(t for t in text.splitlines() if t.startswith("3 ")))
    assert id_map.object_at_row(row) is not None
    again = DeckReader.from_text(text)
    assert DeckWriter(again).to_text() == text
    assert len(model.lines) == 3


def test_writer_numbering_after_edits() -> None:
    project = DeckReader.read(ROOT / "examples" / "spread_3line_chain.dat")
    original = DeckWriter(project).to_text()
    CommandStack().push(RemoveObject(project, project.lines[0], cascade=True))
    extra = project.points.append(FixedPoint("new anchor", position=(10.0, 10.0, -20.0)))
    _, id_map = DeckWriter(project).to_model()
    ids = [id_map.deck_id(line) for line in project.lines]
    assert ids == [2, 3]
    assert (
        id_map.deck_id(extra)
        == max(id_map.deck_id(point) or 0 for point in project.points if point is not extra) + 1
    )
    project.points[0].deck_hints["deck_id"] = project.points[1].deck_hints["deck_id"]
    _, id_map = DeckWriter(project).to_model()
    assert len({id_map.deck_id(point) for point in project.points}) == len(project.points)
    project.line_types[0].name = "renamed chain"
    text = DeckWriter(project).to_text()
    assert '"renamed chain"' in text and text != original


def test_writer_follows_edited_option_values() -> None:
    project = DeckReader.read(ROOT / "examples" / "lazy_wave_vessel_rao.dat")
    text = DeckWriter(project).to_text()
    assert "vesselRAO" in text
    rao = project.motion
    assert isinstance(rao, FloaterRAO)
    project.motion = MotionFile(file=rao.file, reference=rao.reference)
    text = DeckWriter(project).to_text(validate=False)
    assert "motionFile" in text and "vesselRAO" not in text
    project.environment.water_density = 1030.0
    text = DeckWriter(project).to_text(validate=False)
    assert "1030.0" in text
    project.motion = NoMotion()
    text = DeckWriter(project).to_text(validate=False)
    assert "motionFile" not in text and "vesselRef" not in text


def test_writer_writes_every_option_kind() -> None:
    project = full_project()
    env, settings = project.environment, project.settings
    env.waves = JonswapWave(significant_height=3.0, peak_period=9.0)
    env.current = UniformCurrent(velocity=(0.3, 0.0, 0.0))
    env.ramp_time = 10.0
    env.wave_seed = 3
    settings.time_step = 0.02
    settings.duration = 10.0
    settings.modified_newton = True
    settings.solver_relative_tolerance = 1.0e-8
    settings.solver_absolute_tolerance = 1.0e-14
    settings.solver_max_iterations = 30
    settings.solver_backtracks = 12
    settings.solver_spectral_radius = 0.4
    project.seabed.stiffness = 2.0e5
    settings.extra_options.append(ExtraOption(keyword="WriteLog", values=("0",), note="log"))
    text = DeckWriter(project).to_text()
    for fragment in (
        "jonswap 3.0 9.0 3.3 0.0 waves",
        "uniform 0.3 0.0 0.0 current",
        "dynamic_solver 1e-08 1e-14 30 12 0.4",
        "True modified_newton",
        "0 WriteLog - log",
    ):
        assert fragment in text
    env.waves = MultiTrainSea()
    env.waves.trains.append(WaveTrain(wave=JonswapWave(), spreading=4.0))
    env.waves.trains.append(WaveTrain(wave=RegularWave(height=1.0)))
    text = DeckWriter(project).to_text(validate=False)
    assert "wavetrain" in text and "waves" not in text.split("OPTIONS")[1].split("\n")[1]
    env.waves = OchiHubbleWave()
    env.current = ProfileCurrent()
    project.seabed = BathymetrySeabed(file="bathy.txt")
    env.water_depth = None
    project.motion = FloaterMotionRecord(file="motion.txt", reference=(0.0, 0.0, 0.0))
    text = DeckWriter(project).to_text(validate=False)
    project.seabed.friction = 0.5
    text = DeckWriter(project).to_text(validate=False)
    for fragment in (
        "ochihubble",
        "profile",
        "bathymetryFile",
        "vesselMotion",
        "vesselRef",
        "0.5 frictionMu",
    ):
        assert fragment in text


def test_writer_body_rod_and_row_variants() -> None:
    project = full_project()
    chain = project.line_types[0]
    line = project.lines[0]
    body = project.bodies[0]
    project.bodies.append(Point3Body("p3", position=(0.0, -50.0, -10.0)))
    project.bodies.append(Floater("hull", kind="coupled"))
    project.bodies.append(Floater("hull2", kind="floater", row_format="moordyn", inertia=(1.0,)))
    project.bodies.append(
        Rigid6Body("md", row_format="moordyn", inertia=(1.0, 2.0, 3.0), drag_area=(1.0, 2.0))
    )
    project.turbines.append(Turbine(number=2, platform_displacement=(1, 0, 0, 0, 0, 0)))
    project.points.append(TurbinePoint("tp", turbine=2))
    project.points.append(ConnectPoint("c"))
    project.points.append(FloaterPoint("v", kind="floater"))
    line.end_connections.append(EndConnection(end="A", rotation="pinned"))
    line.end_connections.append(
        EndConnection(end="B", rotation="rigid", torsion="rigid", normal=(1, 0, 0), pretwist=5.0)
    )
    project.lines[1].end_connections.append(
        EndConnection(
            end="A",
            rotation="stiffness",
            rotational_stiffness=1e5,
            torsion="free",
            normal=(1, 0, 0),
        )
    )
    project.lines[2].end_connections.append(
        EndConnection(end="A", torsion="stiffness", torsional_stiffness=3.0, normal=(0, 1, 0))
    )
    line.attachments.append(BuoyancyModule(arc_length=10.0, volume=1.0, mass=10.0))
    line.attachments.append(
        ClumpWeight(
            arc_length=20.0, pitch=10.0, last_arc_length=40.0, mass=5.0, axial_drag_area=1.0
        )
    )
    project.equivalent_buoyancy.append(EquivalentBuoyancy(line_type=chain, diameter=0.3))
    project.syrope_history.append(SyropeIC(lines=[line]))
    project.failures.append(Failure(point=project.points[0], lines=[line], time=5.0))
    project.controls.append(Control(channel=2, lines=[line]))
    project.external_loads.append(ExternalLoad(body=body, axes="body", force=(1.0, 2.0, 3.0)))
    rod = project.rods[0]
    for attachment in ("free", "fixed", "coupled", "floater"):
        rod.attachment = attachment
        DeckWriter(project).to_model()
    rod.attachment = "body_pinned"
    rod.body = body
    model, _ = DeckWriter(project).to_model()
    text = model.to_text(validate=False)
    for fragment in (
        "Point3",
        "Coupled",
        "Vessel",
        "Free",
        "Turbine2",
        "Connect",
        "Pinned",
        "Rigid",
        "20.0:10.0:40.0",
        "SYROPE IC",
        "FAILURE",
        "CONTROL",
        "EXTERNAL LOADS",
        "Body1Pinned",
        "EQUIVALENT BUOYANCY",
        "Free  ",
    ):
        assert fragment in text, fragment


def test_writer_axial_models() -> None:
    project = full_project()
    chain = project.line_types[0]
    chain.axial = ViscoelasticAxial(static_stiffness=1e8, dynamic_stiffness=2e8, damping=(1.0, 2.0))
    assert "100000000.0|200000000.0" in DeckWriter(project).to_text(validate=False)
    chain.axial = ViscoelasticAxial(static_stiffness=1e8, alpha_mbl=10.0, beta=0.2)
    assert "100000000.0|10.0|0.2" in DeckWriter(project).to_text(validate=False)
    chain.axial = SyropeAxial(settings_file="settings.dat")
    assert "SYROPE:settings.dat|" in DeckWriter(project).to_text(validate=False)
    chain.bending = BendingModel(
        bending_stiffness=1.0,
        shear_stiffness=2.0,
        torsional_stiffness=3.0,
        rotary_inertia_axial=4.0,
        rotary_inertia_normal=5.0,
    )
    assert "GAs" in DeckWriter(project).to_text(validate=False)


@pytest.mark.parametrize(
    "breaker",
    [
        lambda p: setattr(p.lines[0], "end_a", None),
        lambda p: setattr(p.lines[0].sections[0], "line_type", None),
        lambda p: setattr(p.lines[0].sections[0], "line_type", GenericLineType("x")),
        lambda p: setattr(p.points[2], "body", None),
        lambda p: setattr(p.points[2], "body", Rigid6Body()),
        lambda p: setattr(p.points[3], "rod_end", None),
        lambda p: setattr(p.points[3], "rod_end", Rod().end_a),
        lambda p: setattr(p.rods[0], "rod_type", None),
        lambda p: setattr(p.lines[0], "end_b", FixedPoint()),
        lambda p: setattr(p.bodies[0], "drag_area", (1.0, 2.0, 3.0)),
        lambda p: setattr(p.bodies[0], "inertia", (1.0,)),
        lambda p: setattr(p.line_types[0], "axial", ViscoelasticAxial()),
        lambda p: setattr(p.line_types[0], "axial", SyropeAxial()),
        lambda p: setattr(p.line_types[0], "name", "bad#name"),
        lambda p: setattr(p, "motion", MotionFile()),
        lambda p: p.outputs.channels.append(ObjectChannel(target=p.line_types[0])),
        lambda p: p.outputs.channels.append(NamedChannel(channel="FairTen1")),
        lambda p: p.external_loads.append(ExternalLoad()),
        lambda p: p.failures.append(Failure(lines=[p.lines[0]])),
        lambda p: p.controls.append(Control(lines=[Line()])),
        lambda p: p.equivalent_buoyancy.append(EquivalentBuoyancy()),
        lambda p: p.lines[0].attachments.append(ClumpWeight(pitch=1.0)),
        lambda p: setattr(p.rods[0], "attachment", "body"),
        lambda p: setattr(p.rods[0], "body", Rigid6Body()),
        lambda p: setattr(p.environment, "waves", OddWave()),
        lambda p: setattr(p.environment, "current", OddCurrent()),
    ],
)
def test_writer_reports_the_failing_object(breaker: Any) -> None:
    project = full_project()
    breaker(project)
    with pytest.raises(DeckExportError) as raised:
        DeckWriter(project).to_text(validate=False)
    issues = DeckWriter(project).validate()
    assert len(issues) == 1 and issues[0].layer == 3
    assert raised.value.obj is None or raised.value.obj is issues[0].obj


def test_layer_three_maps_native_errors_to_objects() -> None:
    project = full_project()
    free = project.points.append(FreePoint("negative", mass=-5.0))
    issues = DeckWriter(project).validate()
    assert len(issues) == 1 and issues[0].obj is free and issues[0].layer == 3
    project.points.remove(free)
    project.seabed = BathymetrySeabed(file="bathy.txt")
    issues = DeckWriter(project).validate()
    assert len(issues) == 1 and issues[0].obj is project
    assert "WtrDpth" in issues[0].message or "bathymetry" in issues[0].message.lower()
    found = validate_project(project)
    assert [issue.layer for issue in found] == [3]
    assert validate_project(project, deck=False) == []


def test_write_to_file(tmp_path: Path) -> None:
    project = full_project()
    target = DeckWriter(project).write(tmp_path / "deck.dat")
    assert DeckReader.read(target).title == "hand built"
    with pytest.raises(FileExistsError):
        DeckWriter(project).write(target)
    DeckWriter(project).write(target, overwrite=True)


@model_type("test.axial.odd")
class OddAxial(AxialModel):
    """An axial model subclass the writer does not know."""


class OddPoint(Point):
    """A point kind without a deck type (not registered)."""


def test_reader_paths_of_an_in_memory_model() -> None:
    project = full_project()
    line = project.lines[0]
    body = project.bodies[0]
    project.bodies.append(Point3Body("p3", mass=1.0))
    project.bodies.append(Floater("hull", kind="coupled", mass=1.0))
    project.bodies.append(Floater("hull2", kind="floater", mass=1.0))
    line.end_connections.append(EndConnection(end="A", rotational_stiffness=5.0))
    line.end_connections.append(
        EndConnection(end="B", torsion="stiffness", torsional_stiffness=2.0, normal=(1, 0, 0))
    )
    line.attachments.append(ClumpWeight(arc_length=12.0, mass=100.0))
    project.controls.append(Control(channel=3, lines=[line]))
    project.external_loads.append(ExternalLoad(body=body))
    project.outputs.channels.append(NamedChannel(channel="FairTen9"))
    project.settings.extra_options.append(ExtraOption(keyword="rhoW", values=("abc",)))
    project.settings.solver_relative_tolerance = 1e-8
    project.settings.solver_absolute_tolerance = 1e-12
    project.settings.solver_max_iterations = 20
    project.settings.solver_backtracks = 5
    model, _ = DeckWriter(project).to_model()
    model.options.add("dynamic_solver", 1, 2, 3)
    model.options._rows.append(_b.Option("futurekeyword", ("1",)))
    again = DeckReader.from_model(model)
    assert [type(b).__name__ for b in again.bodies] == [
        "Rigid6Body",
        "Point3Body",
        "Floater",
        "Floater",
    ]
    assert again.lines[0].end_connections[0].rotation == "stiffness"
    assert again.lines[0].end_connections[1].torsion == "stiffness"
    assert isinstance(again.lines[0].attachments[0], ClumpWeight)
    assert again.controls[0].channel == 3
    assert again.outputs.channels[-1].name == "FairTen9"
    assert isinstance(again.outputs.channels[-1], NamedChannel)
    keywords = [row.keyword for row in again.settings.extra_options]
    assert keywords == ["dynamic_solver", "rhoW", "dynamic_solver", "futurekeyword"]
    with pytest.raises(DeckExportError, match="unknown option keyword"):
        DeckWriter(again).to_text(validate=False)
    again.settings.extra_options.remove(again.settings.extra_options[-1])
    assert "dynamic_solver 1 2 3" in DeckWriter(again).to_text(validate=False)


def test_reader_rejects_rows_it_cannot_type() -> None:
    model = DeckModel.new()
    model.add_line_type("odd", diam=0.1, mass=1.0, ea=(1.0, 2.0, 3.0, 4.0))
    with pytest.raises(DeckFormatError, match="unsupported EA"):
        DeckReader.from_model(model)
    model = DeckModel.new()
    model.add_body(1, "Weird", 0.0, 0.0, 0.0)
    with pytest.raises(DeckFormatError, match="unsupported type"):
        DeckReader.from_model(model)
    model = DeckModel.new()
    model.add_point(1, "Weird", 0.0, 0.0, 0.0)
    with pytest.raises(DeckFormatError, match="unsupported type"):
        DeckReader.from_model(model)


def test_writer_refuses_unknown_subclasses() -> None:
    project = full_project()
    project.line_types[0].axial = OddAxial()
    with pytest.raises(DeckExportError, match="cannot be written"):
        DeckWriter(project).to_text(validate=False)
    project = full_project()
    project.points.append(OddPoint("odd"))
    with pytest.raises(DeckExportError, match="no deck point type"):
        DeckWriter(project).to_text(validate=False)


def test_reader_bathymetry_and_density_fallback() -> None:
    project = full_project()
    project.environment.water_depth = None
    project.seabed = BathymetrySeabed(file="bathy.txt", stiffness=3.0e5, friction=0.2)
    project.lines[0].attachments.append(ClumpWeight(arc_length=5.0, mass=1.0, volume=1.0))
    model, _ = DeckWriter(project).to_model()
    again = DeckReader.from_model(model)
    seabed = again.seabed
    assert isinstance(seabed, BathymetrySeabed)
    assert (seabed.file, seabed.stiffness, seabed.friction) == ("bathy.txt", 3.0e5, 0.2)
    assert isinstance(again.lines[0].attachments[0], BuoyancyModule)
