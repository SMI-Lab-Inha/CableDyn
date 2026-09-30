#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Check that CableDyn release identifiers agree across public interfaces."""

from __future__ import annotations

import re
import runpy
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(relative_path: str) -> str:
    return (ROOT / relative_path).read_text(encoding="utf-8")


def capture(relative_path: str, pattern: str, label: str) -> str:
    match = re.search(pattern, read(relative_path), flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"{relative_path}: cannot find {label}")
    return match.group(1)


def main() -> int:
    try:
        version = capture(
            "CMakeLists.txt",
            r"project\(CableDynCore VERSION (\d+\.\d+\.\d+)",
            "CMake project version",
        )
        citation_version = capture("CITATION.cff", r"^version:\s*([^\s]+)$", "citation version")
        release_date = capture(
            "CITATION.cff", r"^date-released:\s*(\d{4}-\d{2}-\d{2})$", "release date"
        )
        abi_version = capture(
            "src/CableDyn_CAPI.f90",
            r"CD_C_ABI_VERSION\s*=\s*(\d+)_C_INT",
            "C ABI version",
        )
        abi_minor = capture(
            "src/CableDyn_CAPI.f90",
            r"CD_C_ABI_MINOR\s*=\s*(\d+)_C_INT",
            "C ABI minor level",
        )
    except ValueError as error:
        print(f"release metadata: FAIL\n- {error}", file=sys.stderr)
        return 1

    major, minor, patch = version.split(".")
    # The manual derives its version from CMakeLists.txt; evaluate its configuration.
    sphinx_conf = runpy.run_path(str(ROOT / "doc" / "conf.py"))

    checks = {
        "CITATION.cff version": citation_version == version,
        "Python project version": f'version = "{version}"' in read("python/pyproject.toml"),
        "Python runtime version": f'__version__ = "{version}"' in read("python/cabledyn/__init__.py"),
        "Fortran banner version": f"CableDyn  v{version}" in read("src/CableDyn_Banner.f90"),
        "C ABI version string": f"CableDyn {version} C-ABI {abi_version}" in read(
            "src/CableDyn_CAPI.f90"
        ),
        "C ABI numeric version": all(
            token in read("src/CableDyn_CAPI.f90")
            for token in (
                f"CD_C_VERSION_MAJOR = {major}_C_INT",
                f"CD_C_VERSION_MINOR = {minor}_C_INT",
                f"CD_C_VERSION_PATCH = {patch}_C_INT",
            )
        ),
        "C header numeric version": all(
            token in read("src/CableDyn_CAPI.h")
            for token in (
                f"CABLEDYN_CAPI_VERSION_MAJOR {major}",
                f"CABLEDYN_CAPI_VERSION_MINOR {minor}",
                f"CABLEDYN_CAPI_VERSION_PATCH {patch}",
            )
        ),
        "C header ABI minor level": f"CABLEDYN_CAPI_ABI_MINOR {abi_minor}" in read(
            "src/CableDyn_CAPI.h"
        ),
        "Python ABI minor level": f"OBJECT_QUERY_ABI_MINOR = {abi_minor}" in read(
            "python/cabledyn/_lib.py"
        ),
        "C header ABI version": f"CABLEDYN_CAPI_ABI_VERSION {abi_version}" in read(
            "src/CableDyn_CAPI.h"
        ),
        "Python ABI version": f"SUPPORTED_ABI = {abi_version}" in read(
            "python/cabledyn/_lib.py"
        ),
        "Fortran C API test version": (
            f"major == {major}_C_INT" in read("tests/test_c_api.f90")
            and f"minor == {minor}_C_INT" in read("tests/test_c_api.f90")
            and f"patch == {patch}_C_INT" in read("tests/test_c_api.f90")
            and f"abi_version == {abi_version}_C_INT" in read("tests/test_c_api.f90")
        ),
        "OpenFAST module version": f"'v{version}', '{release_date}'" in read(
            "src/openfast/CableDyn_OF.f90"
        ),
        "Sphinx version": sphinx_conf["version"] == version,
        "Sphinx release": sphinx_conf["release"] == version,
        "README stable version": f"stable release is `v{version}`" in read("README.md"),
        "manual stable version": f"``v{version}`` is the current stable release" in read(
            "doc/index.rst"
        ),
        "changelog release heading": f"## [{version}] - {release_date}" in read("CHANGELOG.md"),
        "OpenFAST integration version": f"CableDyn v{version})" in read(
            "integration/openfast/README.md"
        ),
        "OpenFAST example banner": f"Running CableDyn (v{version}," in read(
            "examples/openfast/README.md"
        ),
        "Python licence metadata": 'license-files = ["LICENSE"]' in read(
            "python/pyproject.toml"
        ),
        "Python licence file": (ROOT / "python" / "LICENSE").is_file(),
    }
    failures = [label for label, passed in checks.items() if not passed]
    if failures:
        print("release metadata: FAIL", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print(f"release metadata: PASS (CableDyn v{version}, {release_date})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
