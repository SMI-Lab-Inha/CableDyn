# SPDX-License-Identifier: Apache-2.0
"""The three validation layers of a project.

1. **Property checks**: type, finiteness and limits of each value, the
   numeric limits shared with the deck reader (:mod:`cabledyn.project.schema`).
2. **Object invariants**: references resolve, required collections are not
   empty, arc locations lie on their line, names the deck needs are unique.
3. **Native deck rules**: the project is written as a deck and checked by
   :class:`cabledyn.DeckFile`; a failure is mapped back to the object whose
   row failed. No deck rule is restated in the object model.
"""

from __future__ import annotations

from collections.abc import Iterable

from cabledyn.project.deck import DeckWriter
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.project import Project

__all__ = ["errors", "validate_project"]


def validate_project(project: Project, *, deck: bool = True) -> list[Issue]:
    """Validate a project with every layer.

    Parameters
    ----------
    project : Project
        The project.
    deck : bool
        Run layer 3. It runs only when layers 1 and 2 report no error,
        because a deck cannot be written from a structurally broken model.

    Returns
    -------
    list[Issue]
        Layer-1 and layer-2 issues of every object, then layer-3 issues.
    """
    found = project.validate(recursive=True)
    if deck and not any(issue.is_error for issue in found):
        found.extend(DeckWriter(project).validate())
    return found


def errors(issues: Iterable[Issue]) -> list[Issue]:
    """Return the issues that are errors."""
    return [issue for issue in issues if issue.severity is Severity.ERROR]
