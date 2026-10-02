# SPDX-License-Identifier: Apache-2.0
"""The ``.cdproj`` project file.

A project file is a zip archive (or, for version control, an unzipped
folder) holding:

``project.json``
    The whole object graph (physics, construction, appearance and GUI data)
    with uids, SI values (angles in degrees) and a schema version. It is
    strict JSON: non-finite numbers are never written.
``model.dat``
    The deck written from the project, written only when the project passes
    every validation layer. Relative side-file paths (motion, bathymetry,
    WaterKin, Syrope files) are rewritten to resolve from the project folder
    (a folder project) or from the folder that holds the archive (a zip
    project), and the fixed-name kinematics files the solver reads from the
    deck folder are copied into a folder project. Side files that cannot be
    found are listed in ``project.json`` under ``deck.missing_files``.
``assets/``
    Meshes, textures and other files the project uses (reserved).
``results.json``
    References to solver runs (reserved; empty for now).

The project's ``deck_path`` is stored relative to the folder that holds the
project file when it can be. Objects and values this version does not know
are kept and written back unchanged (see
:class:`~cabledyn.project.SerialContext`).
"""

from __future__ import annotations

import io
import json
import os
import tempfile
import zipfile
from collections.abc import Callable
from pathlib import Path, PureWindowsPath
from typing import Any

from cabledyn import deck_file as _deck
from cabledyn.errors import DeckFormatError
from cabledyn.project.base import SerialContext
from cabledyn.project.deck import DeckExportError, DeckWriter
from cabledyn.project.project import Project

__all__ = ["FORMAT", "SCHEMA_VERSION", "ProjectFormatError", "ProjectStore"]

FORMAT = "cabledyn-project"
"""The ``format`` value of ``project.json``."""
SCHEMA_VERSION = 1
"""The schema version written by this release."""

_PROJECT = "project.json"
_DECK = "model.dat"
_RESULTS = "results.json"
_ASSETS = "assets/"
_MAX_DOCUMENT_BYTES = 256 * 1024 * 1024
_RELATIVE_MARK = "deck_path_relative_to"

Migration = Callable[[dict[str, Any]], dict[str, Any]]
_MIGRATIONS: dict[int, Migration] = {}
"""Upgrades from schema version ``n`` to ``n + 1``, applied in turn on load."""


class ProjectFormatError(ValueError):
    """A project file is not a readable CableDyn project, or is from a newer release."""


def _package_version() -> str:
    from cabledyn import __version__

    return str(__version__)


def _reject_constant(token: str) -> Any:
    raise ProjectFormatError(f"{_PROJECT} holds the non-standard number {token}")


def _is_absolute(value: str) -> bool:
    return Path(value).is_absolute() or PureWindowsPath(value).is_absolute()


def _find_object(node: Any, uid: str) -> dict[str, Any] | None:
    stack = [node]
    while stack:
        item = stack.pop()
        if isinstance(item, dict):
            if item.get("uid") == uid and isinstance(item.get("properties"), dict):
                return item
            stack.extend(item.values())
        elif isinstance(item, list):
            stack.extend(item)
    return None


def _side_files(deck: _deck.DeckFile) -> list[str]:
    """Return the side-file paths a deck reads (Syrope settings and path options)."""
    found: list[str] = []
    for row in deck.line_types:
        ea = row.tokens[3]
        if ea.lower().startswith("syrope:"):
            path = ea[len("SYROPE:") :].partition("|")[0]
            if path:
                found.append(path)
    for record in deck.options:
        keyword = record.keyword.lower()
        if keyword not in _deck._PATH_OPTIONS or len(record.values) != 1:
            continue
        value = record.values[0]
        group = _deck._PATH_OPTION_GROUP[keyword]
        if value.lower() in _deck._DISABLED_MOTION_FILE_VALUES and group != "bathymetry":
            continue
        if group == "waterkin":
            if value.lower() in {"none", "seastate"}:
                continue
            try:
                _deck._native_float(value)
            except ValueError:
                pass
            else:
                continue
        found.append(value)
    return found


def _atomic_write(target: Path, data: bytes) -> None:
    handle, temporary = tempfile.mkstemp(prefix=".cdproj-", dir=target.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(data)
        os.replace(temporary, target)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


class ProjectStore:
    """Save and load ``.cdproj`` project files."""

    @staticmethod
    def to_json(project: Project, base: str | os.PathLike[str] | None = None) -> dict[str, Any]:
        """Return the ``project.json`` document of ``project``.

        Objects of unknown type are put back at their place in their
        collection.

        Parameters
        ----------
        project : Project
            The project.
        base : str | os.PathLike | None
            Folder that holds the project file; ``deck_path`` is stored
            relative to it when possible.

        Returns
        -------
        dict
            JSON-compatible data.
        """
        ctx = SerialContext()
        graph = project.to_dict(ctx)
        orphans: list[dict[str, Any]] = []
        placed = sorted(
            (entry for entry in project.unknown_objects if isinstance(entry, dict)),
            key=lambda entry: int(entry["index"]) if isinstance(entry.get("index"), int) else -1,
        )
        for entry in placed:
            owner_uid, slot, index = entry.get("owner"), entry.get("slot"), entry.get("index")
            owner = _find_object(graph, owner_uid) if isinstance(owner_uid, str) else None
            items = None if owner is None else owner["properties"].get(slot)
            if isinstance(items, list) and isinstance(index, int):
                items.insert(min(max(index, 0), len(items)), entry.get("data"))
            else:
                orphans.append(entry)
        document: dict[str, Any] = {
            "format": FORMAT,
            "schema_version": SCHEMA_VERSION,
            "generator": f"cabledyn {_package_version()}",
            "project": graph,
            "unknown_objects": orphans,
        }
        deck_path = project.deck_path
        if base is not None and deck_path and _is_absolute(deck_path):
            try:
                relative = os.path.relpath(deck_path, Path(base).resolve())
            except ValueError:
                relative = None
            if relative is not None:
                graph["properties"]["deck_path"] = Path(relative).as_posix()
                document[_RELATIVE_MARK] = "project_folder"
        return document

    @staticmethod
    def from_json(data: Any, base: str | os.PathLike[str] | None = None) -> Project:
        """Rebuild a project from a ``project.json`` document.

        Parameters
        ----------
        data : object
            The parsed document.
        base : str | os.PathLike | None
            Folder that holds the project file (resolves a relative
            ``deck_path``).

        Raises
        ------
        ProjectFormatError
            If the document is not a project, is from a newer schema, or
            uses one uid for several objects.
        """
        if not isinstance(data, dict) or data.get("format") != FORMAT:
            raise ProjectFormatError("not a CableDyn project document")
        version = data.get("schema_version")
        if not isinstance(version, int) or isinstance(version, bool) or version < 1:
            raise ProjectFormatError(f"invalid schema version {version!r}")
        if version > SCHEMA_VERSION:
            raise ProjectFormatError(
                f"the project uses schema version {version}; this release reads up to "
                f"{SCHEMA_VERSION}"
            )
        while version < SCHEMA_VERSION:
            data = _MIGRATIONS[version](data)
            version += 1
        ctx = SerialContext()
        project = ctx.build(data.get("project"))
        if not isinstance(project, Project):
            raise ProjectFormatError("the document holds no project")
        if ctx.duplicates:
            raise ProjectFormatError(
                f"uid {ctx.duplicates[0]} is used by more than one object; the file is damaged"
            )
        ctx.resolve()
        kept = data.get("unknown_objects", [])
        project.unknown_objects = [*ctx.unknown, *(kept if isinstance(kept, list) else [])]
        project.load_warnings = list(ctx.warnings)
        deck_path = project.deck_path
        if data.get(_RELATIVE_MARK) and base is not None and deck_path:
            project._values["deck_path"] = str((Path(base) / deck_path).resolve())
        return project

    @staticmethod
    def _deck(project: Project, destination: Path, folder: bool) -> tuple[str | None, Any]:
        """Return the rebased deck text (or ``None``) and its ``deck`` record."""
        from cabledyn.project.validation import errors, validate_project

        problems = errors(validate_project(project))
        if problems:
            reasons = "; ".join(str(issue) for issue in problems[:3])
            return None, {"written": False, "reason": reasons, "errors": len(problems)}
        writer = DeckWriter(project)
        try:
            model, _ = writer.to_model()
            if not project.deck_path:
                model.path = destination / _DECK
            text = writer._render(model)
            deck = _deck.DeckFile.from_text(
                text, path=model.path, caller_driven=model.caller_driven
            )
            if destination != model.path.parent:
                deck._rebase_relative_inputs(destination)
                if folder:
                    plan = _deck._plan_companion_copies([deck], destination, overwrite=True)
                    _deck._copy_companions(plan)
        except (DeckExportError, DeckFormatError, FileExistsError, OSError) as exc:
            return None, {"written": False, "reason": str(exc)}
        missing = [
            value
            for value in _side_files(deck)
            if not (Path(value) if _is_absolute(value) else destination / value).exists()
        ]
        return deck.text(), {"written": True, "file": _DECK, "missing_files": missing}

    @classmethod
    def save(
        cls,
        project: Project,
        path: str | os.PathLike[str],
        *,
        folder: bool = False,
        overwrite: bool = False,
    ) -> Path:
        """Write ``project`` to a project file.

        Each file is written to a temporary name and then moved into place,
        so an interrupted save leaves the previous file intact.

        Parameters
        ----------
        project : Project
            The project.
        path : str | os.PathLike
            Target, normally ending in ``.cdproj``.
        folder : bool
            Write an unzipped folder instead of a zip archive.
        overwrite : bool
            Replace an existing file or folder contents.

        Returns
        -------
        pathlib.Path
            The written path.

        Raises
        ------
        FileExistsError
            If ``path`` exists and ``overwrite`` is false.
        """
        target = Path(path).expanduser().resolve()
        if target.exists() and not overwrite:
            raise FileExistsError(f"{target} already exists")
        target.parent.mkdir(parents=True, exist_ok=True)
        destination = target if folder else target.parent
        if folder:
            target.mkdir(exist_ok=True)
        document = cls.to_json(project, target.parent)
        deck_text, document["deck"] = cls._deck(project, destination, folder)
        payload = json.dumps(document, indent=1, allow_nan=False).encode("utf-8")
        results = json.dumps({"runs": []}, indent=1).encode("utf-8")
        if folder:
            _atomic_write(target / _PROJECT, payload)
            _atomic_write(target / _RESULTS, results)
            (target / _ASSETS).mkdir(exist_ok=True)
            deck_file = target / _DECK
            if deck_text is not None:
                _atomic_write(deck_file, deck_text.encode("utf-8"))
            elif deck_file.exists():
                deck_file.unlink()
            return target
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr(_PROJECT, payload)
            archive.writestr(_RESULTS, results)
            archive.writestr(zipfile.ZipInfo(_ASSETS), "")
            if deck_text is not None:
                archive.writestr(_DECK, deck_text)
        _atomic_write(target, buffer.getvalue())
        return target

    @classmethod
    def load(cls, path: str | os.PathLike[str]) -> Project:
        """Read a project file (zip archive or folder).

        Returns
        -------
        Project
            The project; ``Project.load_warnings`` lists anything kept but
            not understood.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.
        ProjectFormatError
            If the file is not a readable project.
        """
        source = Path(path).expanduser().resolve()
        try:
            if source.is_dir():
                document = source / _PROJECT
                if not document.is_file():
                    raise ProjectFormatError(f"{source} has no {_PROJECT}")
                if document.stat().st_size > _MAX_DOCUMENT_BYTES:
                    raise ProjectFormatError(f"{source}: {_PROJECT} is too large")
                text = document.read_bytes().decode("utf-8")
            elif source.is_file():
                try:
                    with zipfile.ZipFile(source) as archive:
                        if archive.getinfo(_PROJECT).file_size > _MAX_DOCUMENT_BYTES:
                            raise ProjectFormatError(f"{source}: {_PROJECT} is too large")
                        text = archive.read(_PROJECT).decode("utf-8")
                except (zipfile.BadZipFile, KeyError) as exc:
                    raise ProjectFormatError(f"{source} is not a CableDyn project file") from exc
            else:
                raise FileNotFoundError(f"{source} does not exist")
            data = json.loads(text, parse_constant=_reject_constant)
            return cls.from_json(data, source.parent)
        except UnicodeDecodeError as exc:
            raise ProjectFormatError(f"{source}: {_PROJECT} is not UTF-8 text") from exc
        except json.JSONDecodeError as exc:
            raise ProjectFormatError(f"{source}: {_PROJECT} is not valid JSON") from exc
        except RecursionError as exc:
            raise ProjectFormatError(f"{source}: {_PROJECT} is nested too deeply") from exc
