# SPDX-License-Identifier: Apache-2.0
"""The ``.cdproj`` project file.

A project file is a zip archive (or, for version control, an unzipped
folder) holding:

``project.json``
    The whole object graph (physics, construction, appearance and GUI data)
    with uids, SI values and a schema version.
``model.dat``
    The deck written from the project, so the driver can run it directly;
    left out (and noted in ``project.json``) when the project is not yet a
    valid deck.
``assets/``
    Meshes, textures and other files the project uses (reserved).
``results.json``
    References to solver runs (reserved; empty for now).

Objects of a type this version does not know are kept verbatim and written
back unchanged.
"""

from __future__ import annotations

import io
import json
import os
import zipfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

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

Migration = Callable[[dict[str, Any]], dict[str, Any]]
_MIGRATIONS: dict[int, Migration] = {}
"""Upgrades from schema version ``n`` to ``n + 1``, applied in turn on load."""


class ProjectFormatError(ValueError):
    """A project file is not a CableDyn project or is from a newer release."""


def _package_version() -> str:
    from cabledyn import __version__

    return str(__version__)


class ProjectStore:
    """Save and load ``.cdproj`` project files."""

    @staticmethod
    def to_json(project: Project) -> dict[str, Any]:
        """Return the ``project.json`` document of ``project``.

        Returns
        -------
        dict
            JSON-compatible data.
        """
        ctx = SerialContext()
        return {
            "format": FORMAT,
            "schema_version": SCHEMA_VERSION,
            "generator": f"cabledyn {_package_version()}",
            "project": project.to_dict(ctx),
            "unknown_objects": list(project.unknown_objects),
        }

    @staticmethod
    def from_json(data: Any) -> Project:
        """Rebuild a project from a ``project.json`` document.

        Raises
        ------
        ProjectFormatError
            If the document is not a project or its schema is newer than
            this release.
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
        ctx.resolve()
        kept = data.get("unknown_objects", [])
        project.unknown_objects = [*(kept if isinstance(kept, list) else []), *ctx.unknown]
        project.load_warnings = list(ctx.warnings)
        return project

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
        document = cls.to_json(project)
        deck_text: str | None = None
        try:
            deck_text = DeckWriter(project).to_text()
            document["deck"] = {"written": True, "file": _DECK}
        except (DeckExportError, DeckFormatError) as exc:
            document["deck"] = {"written": False, "reason": str(exc)}
        payload = json.dumps(document, indent=1, sort_keys=False, allow_nan=True)
        results = json.dumps({"runs": []}, indent=1)
        if folder:
            target.mkdir(parents=True, exist_ok=True)
            (target / _PROJECT).write_text(payload, encoding="utf-8")
            (target / _RESULTS).write_text(results, encoding="utf-8")
            (target / _ASSETS).mkdir(exist_ok=True)
            deck_file = target / _DECK
            if deck_text is not None:
                deck_file.write_text(deck_text, encoding="utf-8")
            elif deck_file.exists():
                deck_file.unlink()
            return target
        target.parent.mkdir(parents=True, exist_ok=True)
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr(_PROJECT, payload)
            archive.writestr(_RESULTS, results)
            archive.writestr(zipfile.ZipInfo(_ASSETS), "")
            if deck_text is not None:
                archive.writestr(_DECK, deck_text)
        target.write_bytes(buffer.getvalue())
        return target

    @classmethod
    def load(cls, path: str | os.PathLike[str]) -> Project:
        """Read a project file (zip archive or folder).

        Returns
        -------
        Project
            The project; ``Project.load_warnings`` lists anything
            skipped, and objects of unknown type are kept for the next save.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.
        ProjectFormatError
            If the file is not a readable project.
        """
        source = Path(path).expanduser().resolve()
        if source.is_dir():
            document = source / _PROJECT
            if not document.is_file():
                raise ProjectFormatError(f"{source} has no {_PROJECT}")
            text = document.read_text(encoding="utf-8")
        elif source.is_file():
            try:
                with zipfile.ZipFile(source) as archive:
                    text = archive.read(_PROJECT).decode("utf-8")
            except (zipfile.BadZipFile, KeyError) as exc:
                raise ProjectFormatError(f"{source} is not a CableDyn project file") from exc
        else:
            raise FileNotFoundError(f"{source} does not exist")
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ProjectFormatError(f"{source}: {_PROJECT} is not valid JSON") from exc
        return cls.from_json(data)
