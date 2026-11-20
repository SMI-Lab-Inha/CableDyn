<!-- SPDX-License-Identifier: Apache-2.0 -->

# CableDyn v0.4.3 release record

This record lists what was built and tested for the v0.4.3 release: the source commit,
toolchains, dependency pins, test results and asset hashes. Case definitions and pass criteria
are in [`VALIDATION.md`](../VALIDATION.md). The OpenFAST and r-test pins, the 23 OpenFAST
patches and the reference LAPACK are those of [v0.4.2](RELEASE_0_4_2.md).

## Source and pins

| Item | Value |
| --- | --- |
| CableDyn | `0.4.3` (2026-11-20), C ABI 1, minor extension 4 |
| Source | The release tree (`release/0.4.3` at `706ab87`, tree `ea6b1f0`): dev with #144 to #188, the version and the experiment records of the release. The final assets were built from this tree; the gates below ran on it, or on the version commit `1b1c469` before it, which has the same sources (the later commit changes only records and documentation). This record is added afterwards |
| OpenFAST base | `2895884d2be01862173c88d70f86b358d2f1a50a`; the 23 patches of [`integration/openfast/patches`](../integration/openfast/patches) applied to a clean checkout by `prepare_openfast_source.ps1` |
| OpenFAST r-test | `dd5feaaaa500ba7283140107806300d551cff0a7` |
| Reference LAPACK/BLAS | 3.12.1, release-tag archive with SHA-256 `2ca6407a001a474d4d4d35f3a61550156050c48016d949f0da0529c0aa052422`, built with IFX for x86-64-v3 (`/fp:precise /arch:CORE-AVX2`) |

## Toolchains

| Item | Version |
| --- | --- |
| GNU Fortran (test builds) | conda-forge GCC 16.2.0, with OpenBLAS 0.3.33 and LAPACK from the conda environment |
| Intel Fortran (Windows release) | IFX 2025.3.2, `/fp:precise /arch:CORE-AVX2 /Qfma-` (value-safe floating point, x86-64-v3, no contraction into fused multiply-adds), `/Qauto` (thread-safe local storage) |
| Microsoft C/C++ | Visual Studio 2022 (release workflow) |
| Python wheels | scikit-build-core; cibuildwheel v4.3.3 (Linux `manylinux_2_28`, macOS arm64); IFX static build on Windows (`release/build_windows_wheel.ps1`) |
| Python and packaging | Python 3.13, build 1.6.1, Twine 7.0.0 (release workflow); Python 3.15 with the same build and Twine locally |

## Test results

| Gate | Result |
| --- | --- |
| GNU Release CTest, including the slow tier, on five OpenBLAS kernels (Prescott, Nehalem, Sandybridge, Haswell and Zen through `OPENBLAS_CORETYPE`; Windows, local, fresh build directory) | 433 of 434 on each kernel before the experiment rerun. The one failure is `validation_records`, which needs provenance records that name the release; it passes on the release tree. 0 compiler warnings |
| GNU Debug CTest (`-fcheck=all`), including the slow tier, with the tests' own timeouts (Windows, local, fresh build directory) | 433 of 434, with the same `validation_records` row; 0 compiler warnings |
| Example fallback scan (`tests/check_example_fallbacks.py --coretype`, each kernel confirmed by name) | all 66 shipped decks, including those of the `moordynC_*` folders, run on each of the five kernels with no fallback note |
| Python package tests (local, native library and driver of the release tree) | 3232 passed, 47 skipped (platform- and tool-dependent cases) |
| Python distributions (`python -m build python`, `twine check`, locally and in the release workflow) | both pass; each carries the Apache-2.0 licence; the wheel holds only `cabledyn` and its metadata |
| CI matrix of the release tree (#189 at `706ab87`) | Linux and Windows, Release and Debug (the pull-request Debug tier, without the slow and long Debug tests, which ran locally above); the Prescott round-off job with the example fallback scan; Python 3.10 and 3.13 on three platforms with coverage and strict type checks; the wheel builds; Sphinx and Read the Docs; all passed. Codex review of the release tree: no findings, no review threads |
| Sphinx documentation (`-W --keep-going`, CI; Read the Docs build of #189) | built with no warnings |
| pre-commit, all files | passed |
| Release metadata | `check_release_metadata.py` passed (CableDyn v0.4.3, 2026-11-20) |
| Release-tree audit | Git archive of `706ab87`: 1086 files, SHA-256 `2ce14b2d64a524b676472f307bf30e9bc0a5c35fae751888c13c5350c1b13575` |
| Link check (README, top-level Markdown, manual sources) | every relative link resolves; every external URL answers, except publisher pages that refuse automated requests, whose DOIs all resolve at doi.org, and the `v0.4.3` release page, which exists once the release is published |
| Experimental comparisons ([`experiments`](experiments/README.md)) | Holcombe (12 cases) and Bergdahl (30 primary and 8 refinement cases) rerun with a GNU Release build of the release sources (`cabledyn.exe` SHA-256 `84a1f1ff4fda48cca59afcf7313b023425e45d7229a1e7db5199310e52fb14de`). Both reproduce their summaries, metrics, refinement figures and the hashed Holcombe static profiles byte for byte. The Bergdahl waveform comparison was not rerun (its measured series was not available); its files and metrics are those of v0.3.0, as `provenance.json` states |
| Bodies suite (102 metrics) | 102 of 102 pass with the same GNU Release build |
| OpenFAST build gate (`integration/openfast/openfast_build_gate.ps1`, a GNU double-precision OpenFAST of the release sources; the IEA-15 MW mixed `CompMooring = 5` case in still air, as no turbulent wind file is on the host) | passed: the module compiles, names itself `OpenFAST-v5.0.0-CableDyn-0.4.3` and runs 60 s with finite output |
| OpenFAST linearisation (`validation/scripts/openfast_lin_ab.py`, IEA-15 MW UMaine twin decks with the r-test `MD_Shared` turbine-1 inputs prepared as `examples/openfast/README.md` describes, and the mixed cable deck) | passed with the release `openfast.exe` (single precision), with the figures of v0.4.2. The CableDyn mooring block is symmetric to 7.2e-6. Surge stiffness is 7.721e4 N/m against MoorDyn-F's 7.716e4 N/m. The dominant fairlead-block entries agree to 1.1e-3, and every entry above 1 % of its block to 2.9e-3. The platform 6×6 stiffness agrees to 1.0e-3 with MoorDyn-F's quasi-static reduction, the lowest nine line-mode frequencies to 0.3 %, and the mixed cable deck's point stiffness to its static solves to 1.1e-4. It also passed with the double-precision build (symmetry 7.2e-6, 6×6 to 1.0e-3, cable deck to 1.0e-4) |
| Shipped examples, v0.4.2 against the release sources: outputs (GNU Release builds and the two release drivers, every output file compared) | 65 decks run with both. Every output file of 64 decks is byte-identical. The cross-flow VIV example `viv_taut_line_crossflow.dat` changes (instantaneous node lift force by up to 29 % of its column's largest value, fairlead tension by up to 0.34 %), as the CHANGELOG states: its phase now converges at second order in the time step (#171) and moves toward the fine-step response. The 66th deck, `viv_taut_cable_inline.dat`, uses the in-line `VIV` section, which v0.4.2 refuses |
| Shipped examples, v0.4.2 against the release sources: wall time (same host, one thread, OpenBLAS Haswell kernel, best of two runs in alternating order, shared host) | GNU Release builds: geometric-mean ratio 0.981 over the 65 decks, total 132.2 s to 132.0 s; the two torsion decks take 22 to 24 % less (#187), and the four decks that are slower by more than 5 % run in 0.03 s or less. Release drivers: ratio 1.013, total 90.2 s to 90.8 s; the five largest ratios (1.09 to 1.19) are on decks that run in 0.03 s or less, and the two VolturnUS-S line-break decks take 1.0 and 1.4 % more (24.2 s) |
| L5 wall-clock records ([`performance`](performance/README.md)) | not re-measured. The records of `VALIDATION.md` and `validation/performance` are those of the v0.4.2 drivers: the raw results of the external rows were not available on the release host. A control run of the CableDyn rows (`performance_records.py run`, one thread, median of three runs after a warm-up, on a shared host about 20 % busy) timed the v0.4.2 release driver before and after the v0.4.3 preflight driver. The 0.4.3 driver takes 1.449 / 2.574 / 5.461 s on the lazy-wave spans against 1.433–1.447 / 2.538–2.566 / 5.345–5.437 s, 4.714 s on the one-hour VolturnUS-S chain against 4.619–4.620 s, and 0.207 s on the WD0050 chain against 0.205 s. That is 0.1 to 2.0 % slower than the nearer 0.4.2 run, and every timed model reproduces its reference response as recorded |

The coupled IEA-15 MW throughput figures of `VALIDATION.md` (September 2026, another host) were
not re-measured: the turbulent wind file of those runs is not on the release host. They say so
in place. The coupled speed-ups of this release (#184, #185, #187) are measured in their
CHANGELOG entries.

## Windows executable gates

Every gate of `build_static_windows.ps1` passed in the release workflow (`release-windows-static.yml`,
run 44749642176, on `706ab87`; the preflight run 44735202382 on `1b1c469` passed the same gates). The
smoke runs used a system-only `PATH`.

| Gate | Result |
| --- | --- |
| Compiler | IFX 2025.3.2 (minimum 2025.3) |
| OpenFAST solution build | one pass: 34 projects succeeded, 0 failed |
| `CableDyn_driver.exe` imports | `KERNEL32.dll`, `SHELL32.dll` and `imagehlp.dll` only; stack reserve 268,435,456 bytes; `--version` reports `Build x86-64-v3 (AVX2)` |
| `openfast.exe` imports | `KERNEL32.dll` and `imagehlp.dll` only; `-v` reports `OpenFAST-v5.0.0-CableDyn-0.4.3` |
| Processor check | both executables stop at start with exit 3 and the message on a simulated missing x86-64-v3 feature, in the workflow and again on the downloaded assets |
| Standalone smoke | the shallow-chain static deck (`chain_catenary_shallow_30m.dat`) converges |
| Memory growth over a 10x longer run | 2.7 -> 2.7 MiB (`dynamic_chain_current.dat`), 21.6 -> 21.7 MiB (`lozon_gomex80_power_cable.dat`); limit +16 MiB |
| Exit contract | a run states the exit contract after the banner; a refused input ends stderr with the closing line |
| `CompMooring = 5` smoke | CableDyn banner, `dtM` 0.025 s (one glue step), equilibrium fairlead tensions 2.436 / 2.436 / 2.436 MN, the same as v0.4.2 to every printed digit, normal termination |

## Clean-machine test

A self-contained check ran the release assets from a folder whose path holds spaces, an
accented letter and Hangul, with a system-only `PATH`. It checks:

- the checksums of the four files and the manifest;
- that both executables start and name 0.4.3, and that both stop at start on a simulated missing
  processor feature;
- every shipped example deck (66, including the `moordynC_*` folders, each run from its own
  folder) with the release `CableDyn_driver.exe`;
- four `CompMooring = 5` runs;
- the Python wheel in a new conda environment, where it imports as 0.4.3 and runs the release
  driver.

The four runs are:

- a 70 s JONSWAP sea with mixed chains and a cable;
- a 400 s line failure;
- a checkpoint at 30 s with a restart whose last 10 s of output (400 rows) must be
  bit-identical;
- a checkpoint between two mooring steps (`dtM` = 2 `DT`, at 30.025 s) whose restart must
  match the rest of the run (200 rows) bit for bit.

The model inputs were rebuilt from the r-test revision above. All 14 checks passed on the
preflight assets and on the final assets.

## Python distribution

The PyPI distribution is built by `.github/workflows/wheels.yml`; each wheel is installed into a
clean environment and smoke-tested on its own platform (the library and the driver it carries
load and run a deck, and the driver reports the expected build). The preflight run (44751227094, on `706ab87`)
built the files below and uploaded them to TestPyPI; the release run rebuilds them from the tag
before the upload. The Windows wheel installed from TestPyPI into a new conda environment has the preflight checksum and passed `release/smoke_test_wheel.py` (library, driver, build x86-64-v3).

| File | SHA-256 (preflight) |
| --- | --- |
| `cabledyn-0.4.3.tar.gz` | `6c85c3aab297487ca9944bed42ebfe1269f063ac5d9431debe105d9bbf8c7a25` |
| `cabledyn-0.4.3-py3-none-win_amd64.whl` | `0b0682262f90ff50431e18e81343f9021a53d52b5ee4be69f9b069b6a29db258` |
| `cabledyn-0.4.3-py3-none-macosx_15_0_arm64.whl` | `e01e156d47764752fd07b447e7e5ad17c2352a0ce3b117ce2d52b06b2315da5a` |
| `cabledyn-0.4.3-py3-none-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl` | `7a363cfc8a1bbd407f1cf73e67650abc7e56a2c615fb7554a1cbcd6e7e27b2e9` |

## Assets

| File | SHA-256 |
| --- | --- |
| `CableDyn_driver.exe` | `5f301f064a69f030b3b9ac7618fce4cbb351883b2be8a2ab2ad88db12cb7aeaa` |
| `openfast.exe` | `8493360a5c30025047a3078a66ee7a7ebe17a8bc0689ad14b0da3a32bd7497fb` |
| `cabledyn-0.4.3-py3-none-any.whl` | `94232ecd0def2b793cb4d5d5a198163f546e9a511925b136e6c67e40724d4498` |
| `cabledyn-0.4.3.tar.gz` | `e9fbcead9867babdc7b392a0a9b454dd5dbef541439045d609b56b583ecc0302` |
