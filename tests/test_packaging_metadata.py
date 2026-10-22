# SPDX-License-Identifier: Apache-2.0
"""The published distribution (pyproject.toml) and the development package (python/pyproject.toml)
describe the same project: their [project] tables agree on everything but the readme and licence
files (paths relative to each file) and the classifiers (the platform wheels name their operating
systems and Fortran; the pure package is OS independent)."""

from __future__ import annotations

import unittest
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        tomllib = None

ROOT = Path(__file__).resolve().parents[1]
PER_FILE = {"readme", "license-files", "classifiers"}


def _project(path: Path) -> dict:
    with path.open("rb") as fh:
        return tomllib.load(fh)["project"]


@unittest.skipIf(tomllib is None, "needs tomllib (Python 3.11+) or tomli")
class PackagingMetadataTest(unittest.TestCase):
    def setUp(self):
        self.published = _project(ROOT / "pyproject.toml")
        self.development = _project(ROOT / "python" / "pyproject.toml")

    def test_the_two_projects_agree(self):
        keys = (set(self.published) | set(self.development)) - PER_FILE
        for key in sorted(keys):
            with self.subTest(key=key):
                self.assertEqual(self.published.get(key), self.development.get(key))

    def test_per_file_fields_point_at_existing_files(self):
        self.assertTrue((ROOT / self.published["readme"]).is_file())
        self.assertTrue((ROOT / "python" / self.development["readme"]).is_file())
        for name in self.published["license-files"]:
            self.assertTrue((ROOT / name).is_file(), name)
        for name in self.development["license-files"]:
            self.assertTrue((ROOT / "python" / name).is_file(), name)

    def test_the_wheel_classifiers_differ_only_in_the_platforms(self):
        published = set(self.published["classifiers"])
        development = set(self.development["classifiers"])
        self.assertEqual(development - published, {"Operating System :: OS Independent"})
        self.assertEqual(
            published - development,
            {
                "Operating System :: Microsoft :: Windows",
                "Operating System :: POSIX :: Linux",
                "Operating System :: MacOS",
                "Programming Language :: Fortran",
            },
        )


if __name__ == "__main__":
    unittest.main()
