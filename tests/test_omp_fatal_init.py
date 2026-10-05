#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Every OpenMP parallel region of the library and the driver prepares its threads for the
abnormal-end report: its first statement (the first statement of the loop body, for a
PARALLEL DO) is ``CALL CD_Fatal_Thread_Init()``, so a stack overflow on a worker thread is
reported like one on the main thread (src/cabledyn_fatal.c)."""

from __future__ import annotations

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIRECTIVE = re.compile(r"^\s*!\$omp\s+parallel\b(?P<rest>.*)$", re.IGNORECASE)
DO_CLAUSE = re.compile(r"^\s*do\b", re.IGNORECASE)
DO_STATEMENT = re.compile(r"^\s*(?:\w+\s*:\s*)?do\b", re.IGNORECASE)
INIT = re.compile(r"^\s*call\s+cd_fatal_thread_init\s*(?:\(\s*\))?\s*(?:!.*)?$", re.IGNORECASE)


def sources() -> list[Path]:
    """The Fortran sources of the core library and the driver (not the OpenFAST tree)."""
    files = sorted((ROOT / "src").glob("*.f90")) + sorted((ROOT / "app").glob("*.f90"))
    return [path for path in files if "openfast" not in path.parent.name.lower()]


def statements(lines: list[str], start: int):
    """Yield (index, text) of the statements after ``start``, skipping comments and blanks."""
    for index in range(start, len(lines)):
        text = lines[index].strip()
        if not text or text.startswith("!"):
            continue
        yield index, lines[index]


def missing_inits(lines: list[str]) -> tuple[int, list[int]]:
    """Count the parallel regions and list the 1-based lines of those without the call."""
    found, missing = 0, []
    index = 0
    while index < len(lines):
        match = DIRECTIVE.match(lines[index])
        if match is None:
            index += 1
            continue
        found += 1
        first = index
        loop = DO_CLAUSE.match(match.group("rest")) is not None
        while lines[index].rstrip().endswith("&"):
            index += 1
        following = statements(lines, index + 1)
        nxt = next(following, None)
        if loop:
            if nxt is None or not DO_STATEMENT.match(nxt[1]):
                missing.append(first + 1)
                index += 1
                continue
            nxt = next(following, None)
        if nxt is None or not INIT.match(nxt[1]):
            missing.append(first + 1)
        index += 1
    return found, missing


class OmpFatalInitTest(unittest.TestCase):
    def test_every_parallel_region_prepares_its_threads(self) -> None:
        total = 0
        problems = []
        for path in sources():
            found, missing = missing_inits(path.read_text(encoding="utf-8").splitlines())
            total += found
            problems += [f"{path.relative_to(ROOT)}:{line}" for line in missing]
        self.assertGreater(total, 0, "no OpenMP parallel region found: the scan is broken")
        self.assertEqual(
            problems,
            [],
            "OpenMP parallel regions whose first statement is not CALL CD_Fatal_Thread_Init()",
        )

    def test_the_scan_catches_a_region_without_the_call(self) -> None:
        good = [
            "    !$OMP PARALLEL DO DEFAULT(SHARED) &",
            "    !$OMP PRIVATE(i)",
            "    DO i = 1, n",
            "      ! comment",
            "      CALL CD_Fatal_Thread_Init() ! note",
            "    END DO",
            "    !$omp parallel",
            "    call cd_fatal_thread_init()",
            "    !$omp end parallel",
        ]
        self.assertEqual(missing_inits(good), (2, []))
        bad = ["  !$OMP PARALLEL DO", "  DO i = 1, n", "    x(i) = 0", "  END DO",
               "  !$OMP PARALLEL", "  x = 1"]
        self.assertEqual(missing_inits(bad), (2, [1, 5]))


if __name__ == "__main__":
    unittest.main()
