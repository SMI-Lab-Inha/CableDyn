<!-- SPDX-License-Identifier: Apache-2.0 -->

# CableDyn v0.4.0 release record

This record lists what was built and tested for the v0.4.0 release: the source commit,
toolchains, dependency pins, test results and asset hashes. Case definitions and pass criteria
are in [`VALIDATION.md`](../VALIDATION.md). The OpenFAST and r-test pins and the reference
LAPACK are those of [v0.3.0](RELEASE_0_3_0.md); the OpenFAST patch set grows from 20 to 22.

## Source and pins

| Item | Value |
| --- | --- |
| CableDyn | `0.4.0` (2026-11-04), C ABI 1, minor extension 3 |
| Source | The release-preparation tree (`release/0.4.0` at `476e2c5`): dev with #50 to #85 and the version, documentation, test-build, driver output-guard and Python changes of the release. The assets and the Windows gates were produced from this tree; afterwards only this record changed, which no build reads. The Debug CTest, the kernel sweep and the bodies suite below used `edd0d62`, the Release CTest `79760e0` (the driver output guard); the later commits change only the Python package, whose suite was re-run at `476e2c5` |
| OpenFAST base | `2895884d2be01862173c88d70f86b358d2f1a50a`; the 22 patches of [`integration/openfast/patches`](../integration/openfast/patches) applied to a clean checkout by `prepare_openfast_source.ps1` |
| OpenFAST r-test | `dd5feaaaa500ba7283140107806300d551cff0a7` |
| Reference LAPACK/BLAS | 3.12.1, release-tag archive with SHA-256 `2ca6407a001a474d4d4d35f3a61550156050c48016d949f0da0529c0aa052422`, built with IFX for x86-64-v3 (`/fp:precise /arch:CORE-AVX2`) |

## Toolchains

| Item | Version |
| --- | --- |
| GNU Fortran (test builds) | conda-forge GCC 16.2.0, with OpenBLAS 0.3.34 and LAPACK from the conda environment |
| Intel Fortran (Windows release) | IFX 2025.3.2, `/fp:precise /arch:CORE-AVX2 /Qfma-` (value-safe floating point, x86-64-v3, no contraction into fused multiply-adds), `/Qauto` (thread-safe local storage) |
| Microsoft C/C++ | Visual Studio 2022 (release workflow) |
| Python wheels | scikit-build-core; cibuildwheel v4.3.3 (Linux `manylinux_2_28`, macOS arm64); IFX static build on Windows (`release/build_windows_wheel.ps1`) |
| Python and packaging | Python 3.13.17, build 1.6.1, Twine 7.0.0 (release workflow) |

## Test results

| Gate | Result |
| --- | --- |
| GNU Release CTest, including the slow tier (Windows, local; a fresh build directory at `edd0d62`, rebuilt at `79760e0`) | 359 of 359 passed at both; 0 compiler warnings |
| GNU Debug CTest (`-fcheck=all`), including the slow tier (Windows, local, fresh build directory) | 359 of 359 passed; 0 compiler warnings |
| Round-off fragility (GNU Release build of the release tree on the OpenBLAS kernels Prescott, Nehalem, Sandybridge, Haswell, SkylakeX and Zen, through `OPENBLAS_CORETYPE`) | on each kernel the full CTest passes 359 of 359, the two kernel-pinned tests (`openfast_touchdown`, `hermite_node_reaction`) pass unpinned, and all 53 shipped example decks run with no fallback, stall or relaxed-tolerance note and the same exit status |
| Python package tests (local, native library and driver of the release tree) | 2167 passed, 47 skipped (platform- and tool-dependent cases) |
| CI matrix of #85 (the last change before the release tree) | Linux and Windows, Release and Debug; Python on three platforms with coverage and strict type checks |
| Sphinx documentation (`-W -n`) | built with no warnings |
| pre-commit, all files | passed |
| Release metadata | `check_release_metadata.py` passed (CableDyn v0.4.0, 2026-11-04) |
| Release-tree audit | worktree and `git archive` of the release tree: 921 files, the same SHA-256 |
| Experimental comparisons ([`experiments`](experiments/README.md)) | Holcombe (12 cases) and Bergdahl (30 primary and 8 refinement cases) re-run with a GNU Release build of the release-preparation tree (`cabledyn.exe` SHA-256 `d6603f989f41b01e7921dc52c1c9e580182208c99b45e6ebdd93377b4cbc289f`). Holcombe reproduces every summary and metric byte for byte. Bergdahl changes because the upper line end lies above the still-water level and 0.4.0 carries the dry weight of that part: RMSE 1.583 to 1.577 N, MAPE 3.13 % to 3.11 %. The Bergdahl waveform comparison was not re-run (its measured series was not available); its files and metrics are those of v0.3.0, as `provenance.json` states |
| Bodies suite (102 metrics) | 102 of 102 pass with a GNU Release build of the release tree |
| OpenFAST linearisation (`validation/scripts/openfast_lin_ab.py`, IEA-15 MW UMaine twin decks and a mixed cable deck) | passed with the release `openfast.exe` (single precision): CableDyn mooring block symmetric to 7.2e-6; surge stiffness 7.721e4 N/m against MoorDyn-F's 7.716e4 N/m; fairlead blocks to 1.1e-3 and the platform 6×6 stiffness to 1.0e-3 of MoorDyn-F's quasi-static reduction; the lowest nine line-mode frequencies to 0.3 %; the mixed cable deck's point stiffness to its static solves to 1.1e-4 |
| Wall time of the shipped examples, v0.3.0 vs dev during the release checks (GNU, same host, one thread, OpenBLAS Haswell kernel) | over the 51 decks both run, the geometric-mean ratio is 0.938 and the total falls from 213.0 to 196.5 s |

## Windows executable gates

Every gate of `build_static_windows.ps1` passed in the release workflow (`release-windows-static.yml`,
run 41709562272). The smoke runs used a system-only `PATH`.

| Gate | Result |
| --- | --- |
| Compiler | IFX 2025.3.2 (minimum 2025.3) |
| OpenFAST solution build | one pass: 34 projects succeeded, 0 failed |
| `CableDyn_driver.exe` imports | `KERNEL32.dll`, `SHELL32.dll` and `imagehlp.dll` only; stack reserve 268,435,456 bytes; `--version` reports `Build x86-64-v3 (AVX2)` |
| `openfast.exe` imports | `KERNEL32.dll` and `imagehlp.dll` only; `-v` reports `OpenFAST-v5.0.0-CableDyn-0.4.0`; reports `CableDyn build: x86-64-v3 (AVX2).` |
| Processor check | both executables stop at start with exit 3 and the message on a simulated missing x86-64-v3 feature |
| Standalone smoke | the shallow-chain static deck (`chain_catenary_shallow_30m.dat`) converges |
| Memory growth over a 10x longer run | 2.6 -> 2.6 MiB (`dynamic_chain_current.dat`), 20.7 -> 20.7 MiB (`lozon_gomex80_power_cable.dat`); limit +16 MiB |
| Exit contract | a run states the exit contract after the banner; a refused input ends stderr with the closing line |
| `CompMooring = 5` smoke | CableDyn banner, `dtM` 0.025 s (one glue step), equilibrium fairlead tensions 2.436 MN, normal termination |

## Clean-machine test

A self-contained check ran the release assets from a folder whose path holds spaces, an accented
letter and Hangul, with a system-only `PATH`. It checks the checksums of the four files and the
manifest, the start of both executables (both name 0.4.0), every shipped example deck (53) with
the release `CableDyn_driver.exe`, three `CompMooring = 5` runs (a 70 s JONSWAP sea with mixed
chains and a cable, a 400 s line failure, and a checkpoint at 30 s with a restart whose last
10 s of output, 400 rows, must be bit-identical), and the Python wheel in a new virtual
environment, where it imports as 0.4.0 and runs the release driver. All 12 checks passed.

## Python distribution

The PyPI distribution is built by `.github/workflows/wheels.yml`; each wheel is installed into a
clean environment and smoke-tested on its own platform (the library and the driver it carries
load and run a deck, and the driver reports the expected build). The preflight run (41709568758)
built the files below; the release run rebuilds them from the tag before the upload. The 0.4.0
dry run through TestPyPI installed the Windows wheel into a new environment and passed the same
smoke test.

| File | SHA-256 (preflight) |
| --- | --- |
| `cabledyn-0.4.0.tar.gz` | `ad0e5e10fd7d6c8ed2788bc58f6d2328500bad16eb1d39984521b461a4a777a8` |
| `cabledyn-0.4.0-py3-none-win_amd64.whl` | `bd08167a311ed1098dbf63052547d922275b46b055c31a220a80d1fcc819e3dc` |
| `cabledyn-0.4.0-py3-none-macosx_15_0_arm64.whl` | `390759cb3b15fc2f86a4d6d21f6d9ae24fc765d54f91be9dac608e4fe866bf7a` |
| `cabledyn-0.4.0-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl` | `d35868bdda8a07bdc052b1372b60e88bf788087e3c5eca01400d357fed328839` |

## Assets

| File | SHA-256 |
| --- | --- |
| `CableDyn_driver.exe` | `470c07fbcbd293b0b38012eb3a55ee919604736b6e6eef13c7900a77fe251502` |
| `openfast.exe` | `02ca8ffd36ecc9a5d47b58df6c028878683dd8e41837bb84c9bc22116980fa77` |
| `cabledyn-0.4.0-py3-none-any.whl` | `c43a37581377699e0337f2178953d33905477ce4beff810dbbb6cfcfb07ea4a9` |
| `cabledyn-0.4.0.tar.gz` | `20398e8cee8fe35547b5f7064470ff9172adeb2dcc9af8ffff1a9013dd6c1bfa` |
