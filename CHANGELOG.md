<!-- SPDX-License-Identifier: Apache-2.0 -->

# Changelog

All notable changes to CableDyn are recorded in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.1] - 2026-10-09

### Added

- Torsion of finite-EI lines in the standalone driver (condensed, quasi-static, uniform torque):
  `END CONNECTIONS` columns `TorsStiffness NxX NxY NxZ [Pretwist]`, a `motionFile` roll column,
  the torque returned to Rigid6 bodies, channels `Torq<L>N<J>`, `Twist<L>N<J>` and `Twist<L>`,
  range-graph columns, `DeckModel` and channel support in Python, and the example
  `torsion_lazy_wave_hangoff_twist.dat`. Validated in VALIDATION.md (Torsion). Routes outside
  its scope, coupled OpenFAST (NLR) runs included, stop with a named error; results without
  torsion are unchanged.

### Changed

- Python: NumPy arrays in the public functions carry explicit element types, and the package
  type-checks under `mypy --strict` on Python 3.10 and on current NumPy. Tests that need the
  native shared library are skipped with a reason when it has not been built.
- The test suite runs in parallel (`ctest -j`) and holds on hosted CI runners: the catenary
  reference solve converges from a robust start on every platform, and the driver-path tests
  check the case-insensitive spelling of an output root only on Windows.

### Fixed

- The standalone driver no longer ends without a word when its process is interrupted or hits a
  fatal fault. Ctrl+C, Ctrl+Break, closing the console window, `SIGINT`, `SIGTERM` and `SIGHUP`
  are reported as `CableDyn_driver: stopped by <cause>`, and an access violation, a stack
  overflow or another fatal fault as `CableDyn_driver: fatal error: <cause>`. Both reports give
  the simulated time of the last committed step, and the exit status is unchanged. On Windows a
  fault report also names the module and offset of the fault and, for an access violation, the
  address read or written. Previously a
  stack overflow in a Windows GNU build ended the process with no message. Every thread that
  enters an OpenMP parallel region keeps room to report its own stack overflow, and on Linux
  and macOS a hardware fault reaches the handler that was there before (the Fortran runtime's
  backtrace, or the default action and its core dump) with its original address and context,
  and a previous handler runs under its own flags and signal mask; a fault signal sent by
  another process (`kill`) is passed on as sent.
  A signal the driver was started with ignored (for example under `nohup`) stays ignored.
- Every non-zero exit the driver makes itself now ends stderr with the closing line
  `CableDyn_driver: ended with exit code <n>`. A process ended from outside (`taskkill /F`,
  *End task*, `kill -9`) runs none of its own code and cannot report anything; the missing
  closing line now identifies it, since `taskkill /F` leaves exit code 1, the same code as a
  refused input. The driver states this at start, after the banner, in an `Exit status:` line.
  `CableDynDriver.run`, `cabledyn-run` and `cabledyn-study` report such a run as ended early,
  with the time its output reached, instead of relaying the start-up log as the error; with a
  driver that does not write the `Exit status:` line (0.1.0 and earlier), a failure is
  reported with the driver's own text as before. The documentation explains how to recognise each kind of early end and warns that
  `taskkill /IM CableDyn_driver.exe /F` ends every CableDyn run on the computer.
- The static solve of a taut, neutrally buoyant finite-EI line (mass per length equal to
  the displaced mass) no longer depends on the sign of the round-off weight: a line whose
  total weight is below 1e-12 of its axial stiffness is seeded as a straight, uniformly
  stretched line. Previously the catenary seed could fail to close on some platforms.
- Windows builds with gfortran: the `modal` test could stop with a stack overflow on some
  CPUs, because it kept 1.7 MB of dense matrices on the stack. Large test arrays and the
  initialisation workspaces of the coupled aggregate and the C API are now allocated on the
  heap, every executable of such a build reserves a 64 MiB stack, and the new `stack` tests
  run the library on a 1 MiB thread stack. Results are unchanged.
- The static Windows `openfast.exe` compiles the CableDyn adapter with `/heap-arrays:1024`,
  like the CableDyn core, so its run-time-sized arrays no longer use the
  `openfast.exe` stack.
- Builds with gfortran 16 are free of `-Wuninitialized` warnings; the reported values (the
  bounds of unallocated components) were never read.
- The EPUB edition of the manual no longer contains the `.nojekyll` marker file.
- Python: `DeckFile` messages for an `END CONNECTIONS` row name the line, the end and the
  offending value, and the unit of a main-output channel written with leading zeros in its ids
  (`Ten01N04`, `Torq03N4`) is recognised, as the driver accepts such names.
- Torsion: an imposed static twist `Pretwist(B) − Pretwist(A)` above 1000 turns (for example
  a value given in the wrong units) now stops with a message naming the line and the value.
  Before, the static twist ramp could run for millions of stages.
- Torsion: installing a torsion description on a dynamic cable now invalidates its step
  snapshot, as a new end connection does, so a later restore cannot rewind to a state taken
  under the previous loads; a rejected description leaves the cable and its snapshot as they
  were, and ending a cable module also clears its torsion frame.

## [0.1.0] - 2026-10-01

First public release. Deck keywords, Python interfaces, and the C ABI are
versioned but may change before `1.0.0`. The validation evidence for every item
below is in
[VALIDATION.md](https://github.com/SMI-Lab-Inha/CableDyn/blob/main/VALIDATION.md).

### Added

#### Lines and materials

- Geometrically nonlinear tension-only (`EI = 0`) elements for chains, wire ropes,
  and fibre ropes, and cubic-Hermite finite-bending elements for dynamic power
  cables, with continuous centreline curvature and no rotational degrees of
  freedom.
- Linear axial stiffness, Kelvin–Voigt damping, series-Kelvin viscoelasticity, and
  the Syrope polyester working-curve model with load-history initialisation. A
  `SYROPE IC` line starts in the static equilibrium solved on its `Tmax0` working
  curve; `Tmean0` is accepted for MoorDyn compatibility but does not enter the
  state.
- Discrete attachments on finite-EI cables (`ATTACHMENTS`): buoyancy modules and
  clumps (mass, displaced volume, normal and axial drag, added mass) at single arc
  lengths or `first:pitch:last` series, loaded by currents and waves in statics and
  dynamics, standalone and coupled. `EQUIVALENT BUOYANCY` remains the smeared
  alternative.
- Topology: fixed, coupled, free, and connected points; clump weights; rigid bodies
  and rods; finite-stiffness, clamped (`Rigid`), or moment-free cable-end
  connections; line failure; and active unstretched-length control for
  tension-only lines.

#### Statics

- Newton solution with catenary seeding, load continuation, mesh sequencing,
  adaptive meshing, section-interface preservation, and flat or gridded seabed
  contact. The initial condition of every run is the static equilibrium, solved
  directly.
- Tension-only lines are seeded from the exact elastic catenary with
  frictionless-seabed complementarity, by one initialiser shared by the static and
  dynamic routes. Slack, taut, vertical, and finely meshed lines, soft or sloped
  seabeds, and anchors resting on the seabed are covered; compressed branches are
  rejected.
- Finite-EI cables are solved from the exact `EI = 0` catenary by continuation in
  `EI`, with energy-minimising stages and a physical-branch audit (no self-contact,
  stable tangent); `cable_statics` (default `continuation`) selects whether the
  continuation or the mesh-sequenced route runs first. A mesh that cannot resolve
  the equilibrium curvature is refined, graded toward rotationally restrained ends,
  and reported, and a mesh finer than `sqrt(EI/EA)` is named before the solve; a solve
  that stalls on such a mesh stops within seconds. A line too long for its span to hang
  as a catenary (over a flat or uniformly sloped seabed) is solved in the equilibrium its
  bending stiffness holds, reached by walking the anchor back into place from where the
  catenary exists. Over a flat seabed in still water, when that walk fails or reaches a
  state that loops through itself, the bow behind the fairlead is traced instead: the
  anchor starts beyond the fairlead, the grounded run is folded back to its side once it
  reaches it, and the anchor is walked into place. When neither is reached the solve stops
  naming the geometry. With `tensile_safety True` the static refinement of a
  compression stops, naming it, once the compression is the reaction of a held end or a
  refined polish cannot follow it, and a walked bow compressed beyond a tenth of its
  peak tension is rejected before any refinement.
- The frictionless contact of a finite-EI line with a sloped structured seabed acts
  along the surface normal, as for tension-only lines, in statics and dynamics: the
  tension of a grounded run changes along the slope by the weight component.
- A steady current (uniform or a depth profile) enters the static equilibrium,
  including the out-of-plane offset. A finite-EI line whose grounded run the current
  drives into compression is solved out of its plane; a current that leaves a line
  no tension-only equilibrium stops the solve with the compression named; in the
  plane of the line that diagnosis is not repeated on a refined mesh or at a looser
  tolerance, and the drag of the current is evaluated in parallel. A cable
  laid on the seabed with friction reaches its equilibrium in any current direction.
  A very strong current on a frictionless seabed can fold the line back past its
  anchor (a hairpin); that stable equilibrium is kept and a note names it. In a current
  along the plane of the line the equilibrium is selected on the deck mesh coarsened to
  the bending length, so a refined `NumSegs` reaches the same equilibrium, and a planar
  equilibrium that is unstable out of its plane gives way to the stable 3D one.
- Free `Rigid6` bodies, free rods, and `Free`/`Connect` points start at their static
  equilibrium (`bodyIC`, default `static`; `deck bodyIC` keeps the deck pose). The
  solve converges from poor deck poses and never returns an unstable equilibrium.
- A failed initialisation names the line or object, the seeds and solver stages
  tried, and the reason.

#### Dynamics

- Implicit generalised-α integration with consistent structural, hydrodynamic,
  contact, and prescribed-motion terms, and internal step subdivision that never
  commits a partial state.
- Finite-EI cables blend the internal forces of the two step ends
  (`alpha_force_blend`, default `True`); `False alpha_force_blend` selects the
  configuration blend of the journal article.
- Rigid6 bodies, `Point3` buoys, free and pinned rods, and `Free`/`Connect` points
  step monolithically with their lines (`bodyScheme`, default `monolithic`): one
  implicit step with a Newton iteration on the object and point accelerations, each
  line returning its end reaction and condensed end stiffness. Bodies and rods are
  integrated at second order without numerical dissipation. By default
  (`accuracy bodySubstep`) the step is sub-divided to the phase accuracy of the
  stiffest body-mooring mode; `none bodySubstep` takes `dtM` as given.
  `staggered bodyScheme` keeps the staggered predictor-corrector for one release;
  decks with a `motionFile` or a `FAILURE` section, and the OpenFAST and C-API
  coupled routes, use it.
- A finite-EI cable end on a Rigid6 body point or a rod end can be clamped or carry a
  rotational spring; the connection turns with the object and returns its moment to
  it, in statics and dynamics.
- Line modal analysis (`nModes`): natural frequencies and mode shapes of every line
  about its static equilibrium, for `EI = 0` and cubic-Hermite lines, written to
  `<root>.modes.out`, by a banded solver that computes only the requested modes.
- The finite-EI tensile audit (`tensile_safety`) judges the element-mean axial force
  averaged over three elements, the force that enters equilibrium, in the dynamic step
  and in the static contract alike, instead of the pointwise resultant, whose
  oscillation grows with the element length. A coarse mesh whose tensions and
  curvatures are converged is no longer read as compressed or refined for it (the Gulf
  of Mexico 80 m lazy wave under heave passes at 48 elements), and a genuine compression is still
  rejected on every mesh.
- A plausibility guard (`maxStrain`, default `0.5`) stops a run whose `EI = 0` lines
  reach a non-finite state or an element strain beyond the bound.
- After its first steps, a coupled step makes no heap allocation.

#### Environment

- Buoyancy, Morison drag, added mass, Froude–Krylov loading, axial hydrodynamic
  loading, wetting, and current profiles. Finite-EI buoyancy at the free surface
  follows the partial immersion of the circular section.
- Waves: Airy; random-phase `jonswap` seas (`WaveSeed`); ISSC/Pierson–Moskowitz,
  Torsethaugen, and Ochi–Hubble spectra; cos-2s directional spreading
  (`WaveSpreading`, `WaveDirections`); up to 16 superposed wave trains
  (`wavetrain`); Dean stream-function waves (`stream`, `StreamOrder`); and a
  half-cosine wave ramp (`rampTime`).
- MoorDyn-C water kinematics: `WaveKin 3` (`wave_elevation.txt`), `WaveKin 7`
  (`wave_frequencies.txt`), and `Currents 1` (`current_profile.txt`). The other
  MoorDyn-C modes fail closed by name.
- Seabed contact on flat and sloped bathymetry (along the surface normal) with the
  same law in statics and dynamics, and stick-slip seabed friction on `EI = 0` and
  finite-EI lines, isotropic (`frictionMu`) or anisotropic (`frictionMuAxial`,
  `frictionMuLateral`); rods and bodies use a velocity-regularised law. The friction
  state is committed with the step, the snapshots, and the OpenFAST checkpoint.

#### Bodies, rods, and points

- Rigid6 bodies with weight at the centre of gravity, buoyancy at the centre of
  buoyancy, the full rigid-body mass matrix, and hydrostatic restoring (`C44`/`C55`)
  about the body's own axes. MoorDyn v2 14-column `BODIES` rows are read. Buoyancy,
  drag, fluid inertia, and added mass follow the submerged fraction of an equivalent
  sphere (`bodyWetting sphere`, default); `bodyWetting moordyn` and
  `bodyHydro moordyn` reproduce MoorDyn's always-wet body.
- Rigid rods (MoorDyn `ROD TYPES`/`RODS`) with loads integrated over `NumSegs`
  segments, exact wet cross-sections at any tilt, and MoorDyn's end added mass, end
  drag, and end fluid inertia, plus the axial Froude–Krylov force on the end caps.
  `ROD TYPES` columns 6–7 are MoorDyn's `CdEnd CaEnd`; the axial side coefficients
  `CdAx CaAx` are optional columns 8–9. `rodHydro moordyn` reproduces MoorDyn's
  surface-piercing rod moment.
- Rod topologies: rods fixed to a body (`Body<N>`), pinned rods (`Pinned`), rods
  pinned to a body (`Body<N>Pinned`), zero-length rods, the MoorDyn `RODS` attachment
  aliases and `R<N>A`/`R<N>B` line attachments, and decks of bodies and rods without
  lines. One standalone march covers any mix of bodies, rods, points, `EI = 0` lines,
  and finite-EI cables.
- External loads (MoorDyn-F `EXTERNAL LOADS`): constant force and linear or
  quadratic translational damping on Rigid6 bodies, in the global or body frame.
- A `FAILURE` row may detach a line from a free Rigid6 body (standalone), for
  accidental-limit-state runs; the detached end moves as a free point.
- Standalone farm decks: a `TURBINES` section runs FAST.Farm `Turbine<J>` decks in
  the standalone driver, with per-turbine motion records in the `motionFile`.

#### Prescribed motion

- `motionFile` records for coupled points, bodies, and rods.
- 6-DOF vessel motion (`vesselMotion`, about `vesselRef`): position and roll, pitch,
  and yaw in the OrcaFlex `Rz·Ry·Rx` convention or a unit quaternion, with velocities
  and accelerations. Every `Coupled`/`Vessel` point moves rigidly with the vessel,
  and a finite or `Rigid` end connection turns with it.
- RAO-driven vessel motion (`vesselRAO`): a displacement RAO table in the OrcaFlex
  convention, driven by the deck's linear sea (regular, spectral, spread, wave-train,
  or WaterKin) with the same component phases and ramp as the line kinematics.

#### Outputs

- `FairTen`/`AnchTen` and the end-node `Ten<L>N<J>` channels report the line-end
  force, as MoorDyn and OrcaFlex do: the end element's force with axial damping and
  bending shear, plus the end node's share of weight, seabed contact, and drag,
  without the end node's inertia. At rest it is the static end reaction on every
  route. Finite-EI node tensions are segment tensions, as on `EI = 0` lines;
  `.elements.out` keeps the pointwise extrema.
- Body, rod, and point channels with MoorDyn-F names (`Body<N>…`, `Rod<N>…`,
  `Point<P>Fx/Fy/Fz/FH`), standalone and in OpenFAST.
- Range graphs (`LINES` `Outputs` flag `r`): the minimum, maximum, and mean of
  effective tension, curvature, bend moment, declination, and seabed clearance at
  every node, in `<root>.Line<L>.range.out`; `RangeStart` leaves out a start-up
  transient.
- Touchdown channels `TDP<L>s`, `TDP<L>x/y/z`, `TDP<L>Lay`, and `TDP<L>Exc` for a
  line grounded at one end, on every route including OpenFAST.

#### Standalone driver

- `CableDyn_driver` with documented input, output, exit-code, and convergence
  behaviour. Inputs are range-checked on every row, and errors name the deck line.
  Every number in a deck row, an option, or an auxiliary file must be finite and
  not subnormal, input magnitudes are bounded where larger ones would overflow the
  solver (listed in the deck format reference), and line attachments are unquoted
  point ids or rod ends read whole. Non-fatal notes flag an axial start-up shock from
  a prescribed motion, a slack line below its anchor on a deck without a seabed, an
  over-soft `kBot`, time steps or meshes too coarse for a finite-EI line, and a
  static layout folded back past its anchor.
- The initialisation report, standalone and coupled to OpenFAST, prints the fairlead
  end-force vector with its own inclination and the line tangent angles on a separate
  `line tangent` row; the two angles differ by
  the end node's share of the distributed load and, on a finite-EI line, the end
  shear.
- Arguments and file names are read as Unicode. The Windows release executables
  run with UTF-8 as their code page, so decks, outputs, and working folders may be
  named in any script whatever the system locale. The driver opens exactly the named
  file or refuses it with a reason (over-long paths and reserved device names, and,
  in a GNU source build, names the ANSI code page cannot spell). It never overwrites
  an input file and locks its output root.

#### OpenFAST module

- `CompMooring = 5`: mixed mooring and power-cable systems, SeaState kinematics, an
  independent CableDyn time step (the largest multiple of the glue `DT` not exceeding
  it), correction iterations, initial platform offset (`PtfmInit`),
  checkpoint/restart, quasi-static linearisation, active line control, and FAST.Farm
  shared moorings.
- `Coupled`/`Vessel` bodies and rods follow the platform as 6-DOF nodes of the
  coupling mesh and return the force and moment of their lines, hydrostatics,
  SeaState Morison loads, contact, external loads, and inertia. Coupled rods sample
  SeaState along their segments. Free Rigid6 bodies and rods share the deck with
  them and start at their static equilibrium with the platform at `PtfmInit`.
- Discrete cable attachments run on the coupled route with the host's SeaState
  kinematics.
- Linearisation fails closed on a deck with Rigid6 bodies or rods.
- The coupled initialisation report and `<root>.CD.static.out` give the static
  equilibrium, solved before the SeaState kinematics act. The requested `OUTPUTS`
  printed at initialisation, the first `<root>.CD.out` row and the first OpenFAST
  output row give the t = 0 values, which include the SeaState load at t = 0 on the
  line end nodes; the console labels each state.

#### C interface

- Opaque model handles, checked lifecycle calls, fluid-kinematics exchange, motion
  input, load output, and C ABI version 1. Initialisation is thread-safe, and steps
  on distinct handles run concurrently. `CableDyn_Close` releases the calling
  thread's OpenMP worker pool. OpenBLAS uses one thread per calling thread by
  default (`CABLEDYN_BLAS_THREADS`); a Windows GNU build that cannot load OpenBLAS
  fails with a message naming the library and each location tried.
- ABI 1 minor extension 1 (`CableDyn_GetAbiMinor`) adds object queries:
  `CableDyn_NObjects`, `CableDyn_GetObjectInfo`, `CableDyn_GetLineValues`,
  `CableDyn_GetPointState`, `CableDyn_GetBodyState`, `CableDyn_GetRodState`, and
  `CableDyn_EvalChannel`, which read through the evaluator behind the output
  channels.
- `CableDyn_InitDeck` runs decks with finite-EI lines, `Point3` buoys, free `Rigid6`
  bodies, and free, fixed, or pinned rods on the coupled aggregate (still water).

#### Python package (`cabledyn`)

- Native coupling through the C interface and standalone execution
  (`cabledyn-run`), with atomic result replacement on `overwrite=True`.
- Deck validation and loss-aware editing (`cabledyn-deck`), parameter grids, and
  batch studies (`cabledyn-study`) with provenance manifests.
- `DeckModel`: an object model to build and edit decks, with reference-safe rename
  and remove, validated through `DeckFile`.
- Readers for CableDyn, OpenFAST (text and binary), and MoorDyn outputs, including
  the output files of coupled runs. A channel that an OpenFAST OutList repeats is
  kept, the repeats renamed `<name>_2`, `<name>_3`, ... with a warning.
- Post-processing (`cabledyn-post`): run comparison, resampling and zero-phase
  filtering, spectra and coherence, rainflow counting and damage-equivalent ranges,
  block and up-crossing extremes with Gumbel and Weibull fits, and line geometry.
- Fatigue and design checks for moorings and dynamic power cables: S-N and T-N
  curves with the published DNV-RP-C203, DNV-OS-E301, and API RP 2SK constants,
  Palmgren-Miner damage with an optional Goodman correction, damage along the arc,
  lifetime damage over the sea states of a batch study, MBR bend checks, and tension
  checks against MBL.
- Range graphs, touchdown histories, results at a time as arc-length profiles, seabed
  and line-to-line clearance, and per-channel and per-line summary tables.
- In-process object views (`CableDyn.lines`, `points`, `bodies`, `rods`) and
  `channel()` for any output channel.
- `Snapshots`: time-indexed geometry of every object with the seabed and free
  surface, recorded in-process or collected from driver output files, stored as
  `.npz`, and played by `animate()` in Matplotlib.

#### Examples, validation, and distribution

- Examples: catenary, taut, and semi-taut moorings; synthetic ropes; connected
  objects; bodies and rods; regular and irregular forcing; vessel motion; the IEA
  15 MW VolturnUS-S platform; and the 80, 200, and 800 m reference lazy-wave power
  cables. The 80 m cable under prescribed heave runs on 64 elements, within 0.5 % of a
  1024-element solution in fairlead tension and peak curvature, in under a second.
- Validation record: analytical tests, comparisons with MoorDyn and OrcaFlex, the
  Holcombe lazy-wave and Bergdahl dynamic-chain experiments, coupled OpenFAST tests,
  refinement studies, performance measurements, and the reproduction of the journal
  article
  ([validation/PAPER_REPRODUCTION.md](https://github.com/SMI-Lab-Inha/CableDyn/blob/main/validation/PAPER_REPRODUCTION.md)).
  OrcaFlex comparisons are kept as the summary values the tests check; the scripts
  that regenerate the full OrcaFlex output for licence holders are included.
- Acknowledgements of MoorDyn, OpenFAST, and OrcaFlex in the README and the manual.
- Statically linked Windows executables (`CableDyn_driver.exe`, `openfast.exe`)
  built with Intel IFX against the reference LAPACK/BLAS; Python wheel and source
  distribution; SHA-256 checksums; a Sphinx user manual; and cross-platform
  continuous integration.

### Fixed

Corrections made during the pre-release review.

- `EXTERNAL LOADS` reads the MoorDyn-F column order `ID Object Fext Blin Bquad CSys`.
  Rows with `CSys` in the third column are still read, and several rows on one object
  add up. IDs run 1, 2, 3, … in row order. A row on a `Point3` body stops with an
  error (external loads apply to Rigid6 bodies only).
- On the coupled OpenFAST route, a deck `waves` option stops with an error, because the
  host's SeaState supplies the waves. A deck `current` is kept as a steady current only on a
  single-turbine, pure `EI = 0` deck without Rigid6 bodies or rods whose SeaState carries no
  waves or current; otherwise it stops with an error. `DeckFile(caller_driven=True)` in the
  Python package applies the parts of this rule that the deck alone decides.
- The `EI = 0` element tension behind the `Ten` channels, range graphs and line tension
  files includes the axial Kelvin–Voigt (`BA`) share, as the end forces already did.
- On sloped bathymetry, `EI = 0` seabed damping and the friction capacity act along the
  floor normal, as on the finite-EI path; a level floor is unchanged.
- The static drag continuation follows its path in one consistent metric near a limit
  point.
- A coupled step that fails after the body and rod update returns bodies and rods,
  as well as lines, to the start of the step. `CableDyn_Step` reports
  `converged = stalled = false` and `n_iter = 0` for any step that does not commit, on
  every route. `CableDyn_EvalChannel` accepts a buffer length with an embedded `NUL`.
- Checkpoint restart and linearisation rebuild from a temporary deck copy written
  beside the `MooringFile` and removed after use, so the files the deck names
  (bathymetry, WaterKin, Syrope tables) resolve as in the original run. The OpenFAST
  module says at initialisation that it provides no VTK visualisation meshes.
- A very large `WaveComponents` × `WaveDirections` request is rejected instead of
  overflowing the component count.
- `OPTION nModes` on a deck that mixes `EI = 0` and finite-EI lines stops with a named
  error; that route writes no modal files.
- MoorDyn friction keywords (`mu_kT`, `mu_kA`, `mc`, `cv`, `FricDamp`,
  `StatDynFricScale`), `CoupledPinned`/`VesselPinned` rods, `N0` line-node channels and
  numbers with a sign inside the mantissa fail with a named error. A relative
  `WaveKinFile` is resolved against the deck folder first, as in MoorDyn-F, then the
  WaterKin file's folder.
- Python: `tension_check` takes `strength_factor` (default 0.95, exported as
  `DNV_OS_E301_STRENGTH_FACTOR`); an `analysis` or `consequence_class` argument that
  does not belong to the chosen standard is rejected.
- Python: `RangeGraph.line_id` may be `None`, and `SeabedClearance` has a `line_id`;
  `PowerSpectrum.window` is stored in canonical form (`"hann"`, `"boxcar"`); MoorDyn-F
  `Node<J>ax/ay/az` channels are read.
- Python: `generate_deck_cases` and `DeckModel.save(rebase=True)` copy the MoorDyn-C
  kinematics files (`wave_elevation.txt`, `wave_frequencies.txt`, `current_profile.txt`)
  next to the new decks, and `cases.json` records them under `companion_files`.
- Python: `DeckFile` rejects the decks the solver rejects: the `WaveKin 3/7` and
  `Currents 1` rules, the `nModes` rules, output-channel and `FAILURE` rules, a
  non-positive `dtWave`, a JONSWAP `gamma` out of range, and an attachment beyond the
  line length.
- Documentation: every command in the manual, the examples and the validation records
  creates its output folder first; the source-build driver name (`cabledyn`) and the
  release name (`CableDyn_driver.exe`) are both given; the validation record quotes the
  values of the shipped example decks and states which tension (end element or
  `FairTen` end force) each comparison uses.

[Unreleased]: https://github.com/SMI-Lab-Inha/CableDyn/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/SMI-Lab-Inha/CableDyn/releases/tag/v0.1.1
[0.1.0]: https://github.com/SMI-Lab-Inha/CableDyn/releases/tag/v0.1.0
