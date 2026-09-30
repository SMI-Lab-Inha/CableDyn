<!-- SPDX-License-Identifier: Apache-2.0 -->

# CableDyn integration for OpenFAST v5

This directory is for developers who build OpenFAST, maintained by NLR (National Laboratory of
the Rockies, formerly NREL), with CableDyn as `CompMooring = 5`. Users
who only run coupled cases should use the released `openfast.exe` and the
[OpenFAST guide](../../doc/openfast.rst).

The patches in [`patches/`](patches/) apply to the public OpenFAST upstream commit
`2895884d2be01862173c88d70f86b358d2f1a50a` (the OpenFAST v5.0.0 release, used by
CableDyn v0.1.0). They add only the OpenFAST host integration; the CableDyn sources come from
this repository.

| File | Purpose |
| --- | --- |
| `patches/*.patch` | Patch series, one file per patch; `0006` and `0016`–`0020` add new CableDyn-owned files, the others modify upstream files |
| `prepare_openfast_source.ps1` | Checks the upstream revision and a clean worktree, then applies the series |
| `integrate_openfast_vs_solution.ps1` | Adds CableDyn to the Intel Fortran Visual Studio solution of the Windows release |
| `vs-build/CableDyn.vfproj.in` | Visual Studio project template used by that script |
| `openfast_build_gate.ps1` | Builds a configured OpenFAST tree and runs a coupled `CompMooring = 5` smoke case |

## Prepare a source tree

```powershell
git clone https://github.com/OpenFAST/openfast.git openfast-cabledyn
git -C openfast-cabledyn checkout 2895884d2be01862173c88d70f86b358d2f1a50a
./integration/openfast/prepare_openfast_source.ps1 -OpenFASTRoot ./openfast-cabledyn
```

The script stops if the revision differs, the worktree carries other changes, or a patch does
not apply. On a checkout where the whole series is already applied it reports that and changes
nothing, so it can be rerun. Configure the patched tree with CMake and pass `-DCABLEDYN_ROOT=<CableDyn
checkout>` (optionally `-DCABLEDYN_BUILD_DIR=<CableDyn build>` to link a prebuilt core). Then
check the build:

```powershell
./integration/openfast/openfast_build_gate.ps1 -BuildDir <openfast-build> -CaseFst <case.fst>
```

The script fails unless `openfast.exe` builds, the smoke case initialises the CableDyn module, and
the output is finite to the last row. An `openfast.exe` built through CMake does not embed the
UTF-8 active-code-page manifest, so on Windows it may not open input files in folders whose names
the system code page cannot spell; the static release build embeds it.

For the static Windows release, pass a clean, unpatched checkout to
`release/build_static_windows.ps1`; it applies the series itself (see
[Installation](../../doc/installation.rst)).

## License of modified files

The patches modify OpenFAST source files licensed under the Apache License, Version 2.0. The
changed files are modified files in the sense of Section 4(b) of that license: upstream copyright
and license notices are retained, and the changed files and the nature of the changes are
identified here and in the top-level [`NOTICE`](../../NOTICE).
