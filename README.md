<!-- SPDX-License-Identifier: Apache-2.0 -->

# CableDyn

[![Release](https://img.shields.io/github/v/release/SMI-Lab-Inha/CableDyn?sort=semver)](https://github.com/SMI-Lab-Inha/CableDyn/releases/latest)
[![Fortran CI](https://github.com/SMI-Lab-Inha/CableDyn/actions/workflows/fortran.yml/badge.svg)](https://github.com/SMI-Lab-Inha/CableDyn/actions/workflows/fortran.yml)
[![Python CI](https://github.com/SMI-Lab-Inha/CableDyn/actions/workflows/python.yml/badge.svg)](https://github.com/SMI-Lab-Inha/CableDyn/actions/workflows/python.yml)
[![Documentation](https://readthedocs.org/projects/cabledyn/badge/?version=latest)](https://cabledyn.readthedocs.io/en/latest/)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.oceaneng.2026.128332-blue.svg)](https://doi.org/10.1016/j.oceaneng.2026.128332)

CableDyn is an open-source finite-element solver for the static and dynamic
analysis of offshore mooring lines and dynamic power cables. It is written in
Fortran 2018 and runs as a standalone program, as a mooring module of OpenFAST,
maintained by NLR (National Laboratory of the Rockies, formerly NREL)
(`CompMooring = 5`), through a C interface, or from Python. It reads input decks
in the MoorDyn v2 format and solves them with two line formulations:

- **Tension-only lines** (`EI = 0`): chains, wire ropes, and synthetic-fibre
  moorings.
- **Finite-bending lines** (finite-EI): cubic-Hermite elements for dynamic power
  cables, including lazy-wave configurations, with continuous centreline curvature.

The current stable release is `v0.2.0`. Interfaces are
versioned and may change before `1.0.0`.

## Capabilities

| Area | Features |
|---|---|
| Statics | Catenary seeding, Newton solution, load continuation, mesh sequencing, adaptive meshing |
| Dynamics | Implicit generalised-alpha integration, prescribed support motion, consistent mass, automatic step subdivision |
| Environment | Buoyancy, Morison drag, added mass, Froude–Krylov loading, regular and irregular waves, current profiles, seabed contact with friction |
| Materials | Linear axial stiffness, Kelvin–Voigt damping, viscoelastic (series-Kelvin; MoorDyn `ElasticMod` 2/3), Syrope polyester model |
| Topology | Fixed, coupled, free, and connected points; clump weights; rigid bodies and rods; end rotational stiffness; active line-length control; line failure |
| Outputs | Effective tension, position, velocity, acceleration, curvature, bending moment, end angles, contact response, solver diagnostics |
| Interfaces | Standalone driver, OpenFAST and FAST.Farm, C ABI, Python package |

Section and keyword names follow MoorDyn v2; unsupported combinations stop with
a named error rather than running approximately.

## Quick start

1. Download from the
   [latest release](https://github.com/SMI-Lab-Inha/CableDyn/releases/latest):

   | Asset | Contents |
   |---|---|
   | `CableDyn_driver.exe` | Standalone solver for Windows x64 |
   | `openfast.exe` | OpenFAST with CableDyn available as `CompMooring = 5` |
   | `cabledyn-0.2.0-py3-none-any.whl` | Python package (pre- and post-processing, batch runs) |
   | `cabledyn-0.2.0.tar.gz` | Python source distribution |
   | `SHA256SUMS.txt` | Checksums for the assets above |

   Both executables are statically linked, need no installer, compiler, or runtime
   DLLs, and run on Windows 10 (version 1903 or later) and Windows 11, from folders
   named in any script. A turbine controller used by an OpenFAST model still needs
   its own `DISCON.dll`. The example decks are in the `examples` folder of the
   release's **Source code (zip)** archive: unzip it and copy `examples` next to the
   executables.

2. Verify the download (compare with `SHA256SUMS.txt`):

   ```powershell
   Get-FileHash CableDyn_driver.exe -Algorithm SHA256
   ```

3. Solve an example deck. The driver finds the static equilibrium, then runs the
   dynamics when the deck sets `dtM` and `TMax`:

   ```powershell
   .\CableDyn_driver.exe examples\wd0050_chain.dat wd0050
   ```

   Results are written to `wd0050.out` and related tables in the current folder.
   The driver does not create folders: an output root such as `results\wd0050`
   needs `New-Item -ItemType Directory -Force results | Out-Null` first.

4. Optionally, script runs from Python. `cabledyn-run` needs to find the driver:
   pass `--executable`, set the `CABLEDYN_DRIVER` environment variable, or put the
   driver on `PATH`. A relative output root is resolved against the deck's folder, so
   this run writes `examples\wd0050.out`:

   ```powershell
   python -m pip install "cabledyn-0.2.0-py3-none-any.whl[post]"
   cabledyn-run examples\wd0050_chain.dat wd0050 --executable .\CableDyn_driver.exe
   ```

   A source build names the driver `cabledyn` (see
   [Build from source](#build-from-source)); the release executable is
   `CableDyn_driver.exe`.

For OpenFAST, set `CompMooring = 5` in the `.fst` file and point `MooringFile` at
a CableDyn deck; ready-to-run cases are in
[`examples/openfast`](https://cabledyn.readthedocs.io/en/latest/tutorial_openfast.html). Every shipped deck is listed
in [`examples/`](https://cabledyn.readthedocs.io/en/latest/examples.html).

## Documentation

The manual at <https://cabledyn.readthedocs.io> covers installation, tutorials,
the deck format, options, outputs, theory, OpenFAST coupling, the C API, Python,
validation, and troubleshooting.

## Build from source

```powershell
conda env create -f environment.yml
conda activate cabledyn
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build -L fortran -LE slow --output-on-failure
```

The driver is `build/cabledyn` (`build\bin\cabledyn.exe` with the conda toolchain
on Windows); the release
ships the same program as `CableDyn_driver.exe`.

See the [development guide](https://cabledyn.readthedocs.io/en/latest/development.html)
for the full test suite, the Python package, and the OpenFAST build; the source
layout is described in the architecture notes in the repository.

## Validation

CableDyn is verified against analytical solutions and compared with MoorDyn,
OrcaFlex, two physical experiments (Holcombe lazy-wave, Bergdahl dynamic chain),
the IEA 15 MW VolturnUS-S reference systems, and coupled OpenFAST runs. Results,
acceptance criteria, and reproduction commands are in
the [validation record](https://cabledyn.readthedocs.io/en/latest/validation.html); data and scripts are in
[`validation/`](https://cabledyn.readthedocs.io/en/latest/validation.html).

## Repository layout

| Directory | Contents |
|---|---|
| `app/` | Standalone driver (`CableDyn_driver.exe`) |
| `src/` | Solver core, C interface, OpenFAST module |
| `python/` | Python package |
| `tests/` | CTest unit, integration, and validation tests |
| `doc/` | Sphinx manual |
| `examples/` | Standalone and OpenFAST input decks |
| `integration/openfast/` | OpenFAST patch series and build scripts |
| `release/` | Release checks and the static Windows build |
| `validation/` | Validation records, reference data, and scripts |

## Citation

If you use CableDyn in academic work, please cite:

> J. H. Seo, J. Lim, K. Shim, J. Song. CableDyn: Curvature-resolving implicit
> finite-element analysis of mooring lines and dynamic power cables for floating
> offshore wind. *Ocean Engineering* 368 (Part 2), 128332, 2026.
> <https://doi.org/10.1016/j.oceaneng.2026.128332>

To identify the software version, also cite the release (CableDyn 0.2.0).
The [citation file](https://cabledyn.readthedocs.io/en/latest/citing.html) contains both entries (GitHub **Cite this
repository**); BibTeX is on the manual's
[How to cite](https://cabledyn.readthedocs.io/en/latest/citing.html) page.

## Contributing and support

Questions, bug reports, and feature requests are welcome in
[GitHub Issues](https://github.com/SMI-Lab-Inha/CableDyn/issues); the
support notes in the repository say where each belongs. Read the
[contributing guide](https://cabledyn.readthedocs.io/en/latest/contributing.html)
before opening a pull request. Everyone taking part is expected to follow the
[Code of Conduct](https://github.com/SMI-Lab-Inha/CableDyn?tab=coc-ov-file). Report
security issues privately as described in the
[security policy](https://github.com/SMI-Lab-Inha/CableDyn/security/policy).

If you extend or fix CableDyn in your own copy, please consider contributing the
change here, so that every user benefits and the validation record stays in one
place.

CableDyn is research software. Users are responsible for model selection,
convergence checks, and compliance with the design codes that govern their work.

Maintainer: Prof. Jae Hoon Seo, Department of Naval Architecture and Ocean
Engineering, Inha University, Republic of Korea, <jaehoon.seo@inha.ac.kr>

## Acknowledgements

CableDyn builds on the work of others, and we are grateful to the developers and
maintainers of MoorDyn (Hall et al.), whose open input format CableDyn reads and
which serves as a comparison reference; of OpenFAST (NLR), which hosts CableDyn as
a mooring module; and of OrcaFlex, the industry reference used for many of the comparisons
in the [validation record](https://cabledyn.readthedocs.io/en/latest/validation.html). OrcaFlex is a product of Orcina Ltd.

## License and name

CableDyn is distributed under the Apache License 2.0. See [LICENSE](LICENSE) and
[NOTICE](NOTICE).

"CableDyn" identifies the official project in this repository. Modified versions
are welcome under the licence, but, as Section 6 of the Apache License 2.0
provides, they may not use the name in a way that implies they are the official
CableDyn or endorsed by its maintainers; describing a work as "based on CableDyn"
is fine.

This software is developed independently at the Department of Naval Architecture
and Ocean Engineering, Inha University. It is not affiliated with, derived from,
or endorsed by the EPSRC-funded research project (2022–2025) of the same name
led by the University of Edinburgh, "CableDyn: Subsea Power Cable Dynamics Under
Complex Ocean Environment" ([cabledyn.eng.ed.ac.uk](https://cabledyn.eng.ed.ac.uk)).
