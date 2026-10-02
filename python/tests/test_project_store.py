# SPDX-License-Identifier: Apache-2.0
"""The .cdproj project file."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path
from typing import Any

import pytest

from cabledyn.project import (
    SCHEMA_VERSION,
    CameraBookmark,
    ChainType,
    DeckReader,
    DeckWriter,
    FixedPoint,
    Group,
    Line,
    Project,
    ProjectFormatError,
    ProjectStore,
    UnitSystem,
)
from cabledyn.project import store as store_module

EXAMPLE = Path(__file__).resolve().parents[2] / "examples" / "lazy_wave_buoyancy_modules.dat"


def decorated() -> Project:
    project = DeckReader.read(EXAMPLE)
    chain = project.line_types.append(ChainType("R4 chain", grade="R4", nominal_diameter=0.1))
    chain.appearance.set_value("finish", "rusted")
    project.lines[0].colour = "#ff0000"
    project.studio.groups.append(Group("cables", members=list(project.lines)))
    project.studio.cameras.append(CameraBookmark("overview", azimuth=45.0))
    project.studio.layouts = {"main": {"docks": ["browser", "properties"]}}
    project.studio.set_unit_system(UnitSystem.engineering())
    return project


@pytest.mark.parametrize("folder", [False, True])
def test_save_and_load_keep_the_whole_graph(tmp_path: Path, folder: bool) -> None:
    project = decorated()
    target = ProjectStore.save(project, tmp_path / "model.cdproj", folder=folder)
    loaded = ProjectStore.load(target)
    assert loaded.to_dict() == project.to_dict()
    assert loaded.load_warnings == []
    assert DeckWriter(loaded).to_text() == DeckWriter(project).to_text()
    assert loaded.studio.unit_system() == UnitSystem.engineering()
    assert loaded.studio.groups[0].members[0] is loaded.lines[0]
    if folder:
        assert (target / "model.dat").read_text(encoding="utf-8") == DeckWriter(project).to_text()
        assert (target / "assets").is_dir()
        assert json.loads((target / "results.json").read_text(encoding="utf-8")) == {"runs": []}
    else:
        with zipfile.ZipFile(target) as archive:
            names = set(archive.namelist())
        assert {"project.json", "model.dat", "results.json", "assets/"} <= names


def test_save_refuses_to_overwrite(tmp_path: Path) -> None:
    project = decorated()
    target = ProjectStore.save(project, tmp_path / "a.cdproj")
    with pytest.raises(FileExistsError):
        ProjectStore.save(project, target)
    ProjectStore.save(project, target, overwrite=True)
    folder = ProjectStore.save(project, tmp_path / "b.cdproj", folder=True)
    project.lines[0].end_a = None
    ProjectStore.save(project, folder, folder=True, overwrite=True)
    assert not (folder / "model.dat").exists()
    document = json.loads((folder / "project.json").read_text(encoding="utf-8"))
    assert document["deck"]["written"] is False
    assert ProjectStore.load(folder).lines[0].end_a is None


def test_invalid_project_still_saves(tmp_path: Path) -> None:
    project = Project("draft")
    project.points.append(FixedPoint("lonely"))
    project.lines.append(Line("unfinished"))
    target = ProjectStore.save(project, tmp_path / "draft.cdproj")
    with zipfile.ZipFile(target) as archive:
        assert "model.dat" not in archive.namelist()
        document = json.loads(archive.read("project.json"))
    assert document["schema_version"] == SCHEMA_VERSION
    assert document["format"] == "cabledyn-project"
    assert ProjectStore.load(target).to_dict() == project.to_dict()


def test_unknown_object_types_are_kept(tmp_path: Path) -> None:
    project = decorated()
    document = ProjectStore.to_json(project)
    document["project"]["properties"]["points"].append(
        {"type": "plugin.special_point", "uid": "abc", "properties": {"size": 3}}
    )
    loaded = ProjectStore.from_json(document)
    assert len(loaded.points) == len(project.points)
    assert loaded.unknown_objects[0]["data"]["type"] == "plugin.special_point"
    assert any("unknown type" in warning for warning in loaded.load_warnings)
    again = ProjectStore.from_json(ProjectStore.to_json(loaded))
    assert again.unknown_objects == loaded.unknown_objects
    target = ProjectStore.save(loaded, tmp_path / "kept.cdproj")
    assert ProjectStore.load(target).unknown_objects == loaded.unknown_objects


def test_bad_documents(tmp_path: Path) -> None:
    with pytest.raises(ProjectFormatError, match="not a CableDyn project document"):
        ProjectStore.from_json({"format": "other"})
    with pytest.raises(ProjectFormatError, match="invalid schema version"):
        ProjectStore.from_json({"format": "cabledyn-project", "schema_version": "1"})
    with pytest.raises(ProjectFormatError, match="reads up to"):
        ProjectStore.from_json({"format": "cabledyn-project", "schema_version": 99})
    with pytest.raises(ProjectFormatError, match="holds no project"):
        ProjectStore.from_json(
            {"format": "cabledyn-project", "schema_version": 1, "project": {"type": "line"}}
        )
    with pytest.raises(FileNotFoundError):
        ProjectStore.load(tmp_path / "missing.cdproj")
    junk = tmp_path / "junk.cdproj"
    junk.write_text("not a zip", encoding="utf-8")
    with pytest.raises(ProjectFormatError, match="not a CableDyn project file"):
        ProjectStore.load(junk)
    empty = tmp_path / "empty.cdproj"
    empty.mkdir()
    with pytest.raises(ProjectFormatError, match=r"has no project\.json"):
        ProjectStore.load(empty)
    (empty / "project.json").write_text("{broken", encoding="utf-8")
    with pytest.raises(ProjectFormatError, match="not valid JSON"):
        ProjectStore.load(empty)


def test_migrations_run_in_order(monkeypatch: pytest.MonkeyPatch) -> None:
    project = decorated()
    document = ProjectStore.to_json(project)
    seen: list[int] = []

    def upgrade(data: dict[str, Any]) -> dict[str, Any]:
        seen.append(data["schema_version"])
        return {**data, "schema_version": data["schema_version"] + 1}

    monkeypatch.setattr(store_module, "SCHEMA_VERSION", 2)
    monkeypatch.setitem(store_module._MIGRATIONS, 1, upgrade)
    loaded = ProjectStore.from_json(document)
    assert seen == [1]
    assert loaded.to_dict() == project.to_dict()
