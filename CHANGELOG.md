<!-- SPDX-License-Identifier: Apache-2.0 -->

# Changelog

All notable changes to CableDyn are recorded in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Seabed friction in the static solve, toward the as-laid position (OrcaFlex's *as laid*
  statics seabed friction policy, opt-in): `asLaid staticsFriction`. Every grounded node is
  held by a friction spring toward the position at which the line was laid from its lower end
  (the anchor, or the fairlead when that is lower) along the lay azimuth on the seabed, with the regularised Coulomb law of the dynamics
  (isotropic or axial/lateral, capped at `mu N`), in still water and in a current, on `EI = 0`
  and finite-EI lines, standalone, coupled OpenFAST and through the C API; the march continues
  the static springs, so a held run stays at rest. New OPTIONS `staticsFriction none|asLaid`,
  `layAzimuth` and `layAzimuth<ID>` (default: the horizontal direction of each line's chord,
  from its lower end),
  `asLaidTension` and `asLaidTension<ID>` (default 0). The default stays `none`, as in MoorDyn:
  frictionless statics in still water (friction from the still-water equilibrium in a current),
  conservative for anchor loads, and every deck without the new rows runs as before. With
  `asLaid`, a 700 m chain on a 5° slope with `frictionMu 1`, which the frictionless static
  folds 24.0 m back past its anchor, holds where it was laid (friction up to 0.087 `mu N`, the
  slope; FairTen 714.2 → 622.9 kN; held for 2 s to 1.5e-10 of FairTen), and a 650 m finite-EI
  cable on a 10° slope, rejected as too long for its span, solves. As an illustration of the
  effect on anchor loads, `asLaid` in `iea15mw_umaine_mixed_cabledyn.dat` lowers the chains'
  AnchTen from 1351.6 to 368.2 kN (FairTen +0.16 %) and the laid cable's from 1230 N to 0.36 N:
  runs laid slack and pulled taut are held back by friction. A run laid at the tension it
  carries (`asLaidTension`) mobilises none: on a flat seabed FairTen equals the frictionless
  one to 6e-6 and AnchTen to 0.4 %. A finite-EI line laid past the point under its fairlead has
  no static solution yet (its touchdown hook is not reached); a lay azimuth across the span
  solves it. A lay azimuth across the vertical plane of the chord, or a bathymetry sloping across
  it under the laid run, solves a finite-EI line in 3D from the start (so a march from it stays
  at rest), and a finite-EI line is laid from its lower end like an `EI = 0` one; a line whose
  fairlead lies below its other end is laid, and seeded, from the fairlead. With discrete
  `ATTACHMENTS` the reference sinks each node under its discrete load (its share of the line's
  weight plus the attachment's net weight), so a module node rests on the laid path like the rest
  of the run instead of taking its sink for a friction deflection. On a mesh the finite-EI
  statics refines (adaptive refinement, the bending-length and stability refinements) or
  coarsens, the reference is laid on that mesh itself, under its own lumped weights and nodal
  seabed stiffnesses, and the march takes it on the final mesh: a run refined 16x on a seabed
  curved along it and sloping across it, with a section change on the grounded run, rests where
  its twin declared on the refined mesh does (to solver precision), where a reference
  interpolated from the deck nodes moved it by up to 1 mm and reversed its lateral anchor force. The policy is rejected
  where the statics carry no seabed friction: on decks with bodies, rods or dynamic
  (`Connect`/`Free`) points (`frictionMu` itself is rejected on the latter and accepted on the
  default monolithic body and rod march), and on a finite-EI deck whose line has a moving End B
  (a line shared between two turbines, which takes the compatibility route). The Python `DeckFile`, `DeckWriter` and `DeckModel` take the new
  rows. See the theory page (Friction in statics).
- Run-time profile of the standalone driver. With the environment variable
  `CABLEDYN_PROFILE=1`, the driver prints at the end of the run how its wall time divides
  between input parsing, initialisation and statics, and the dynamic march. The march is split
  into the fluid sampled by the driver, the line, cable and object steps, the checks, the
  outputs, and the progress report. The profile also gives the committed steps and the time per
  step, the tension-only line steps and their Newton iterations, and (outside the mixed
  `EI = 0` + finite-EI route) the finite-EI cable steps with their residual, tangent and solve
  counts and times. Every output file is the same with and without it. See the manual's
  command-line page (Run-time profile).
- A note when a finite-EI line starts its dynamics from a static state that is unstable out of
  its vertical plane. A line solved in that plane (free in torsion, no current or end direction
  across it) has its out-of-plane DOFs held in the statics and released in the dynamics. On
  elements longer than half the bending length the touchdown can make that state unstable out of
  the plane, which the planar solve does not see. One stability test with the dynamic
  constraints now runs at the hand-off (0.35 ms on an 87-element lazy wave, 6 ms on 2048
  elements; on a line with `ATTACHMENTS`, on the discrete equilibrium after the attachment
  polish) and the note names the line, its mesh and the factor that meets the bending length,
  and recommends a finer mesh or `adaptive_mesh True`. Results do not change.
- Torsion of finite-EI lines in coupled OpenFAST (`CompMooring = 5`) and FAST.Farm runs. The
  End A frame of a line restrained in torsion turns with the platform orientation, the torque
  moment is returned to the platform, `PtfmInit` twists the line in the static initialisation,
  the `Torq`/`Twist` channels are written, and the unwrapped twist is restored on a checkpoint
  restart. The C API and the standalone deck that mixes `EI = 0` and finite-EI lines without a
  body, which move the line ends by translation only, still refuse torsion columns by name, and
  so does OpenFAST linearisation of a line restrained in torsion. See the manual's OpenFAST page
  (Line torsion).

### Changed

- Faster wave kinematics of directional and multi-train seas (`wavetrain` rows, spectral
  `waves` rows, `WaveSpreading`). Each component's phase cosine is now computed once for the
  surface elevation and the kinematics, and its depth profile with one exponential instead of
  two. `examples/chain_two_train_sea.dat` runs 17 % faster with the release build. The formulas
  are unchanged. With gfortran the results are bit-identical. With the Intel compiler of the
  release, which had computed the second cosine through a fused sine/cosine evaluation, they
  can differ in the last printed digit.
- Faster monolithic body steps with a finite-EI cable clamped or elastically connected to a
  body or rod. The cable's finite-difference columns of the junction Jacobian, each of which
  re-steps the whole cable, are now kept and reused across steps like the mooring lines'
  condensed tangent. They are rebuilt when 20 steps old, after a step that needed more than two
  outer iterations, when the step size changes, and from the third outer iteration of a step.
  `examples/buoy_clamped_cable.dat` runs 66 % faster with the release build, with the same
  outer iterations per step. The junction converges to the same tolerance, so results differ
  only within it, and by no more than the tolerance already allowed against a tightly converged
  run: tensions by about 2e-6 relative and the body position by under 1 µm.
- The OpenFAST linearisation check (`validation/scripts/openfast_lin_ab.py`) now compares
  CableDyn's quasi-static mooring stiffness entry by entry with stock MoorDyn's zero-frequency
  reduction `D - C A^-1 B`, and fails above 0.5 %. It writes the `.lin` files with the
  full-precision `OutFmt` `ES24.15E3`, lets MoorDyn settle its initial condition, and
  linearises at t = 0. Its default case is the `examples/openfast` IEA-15 UMaine twin. On
  that case the per-fairlead stiffness agrees to 0.07–0.10 % (0.29 % for the transverse
  term), and the platform 6 × 6 diagonals agree to 0.10 % or better. The earlier record had
  called MoorDyn's reduction unusable; with the 4-digit default `OutFmt` the same case
  reproduces that asymmetric reduction, so the check now refuses a reference written below
  full precision. VALIDATION.md and the OpenFAST page of the manual are updated.
- The FAST.Farm shared-mooring comparison with stock MoorDyn now covers 600 s of the r-test
  `MD_Shared` farm (`validation/scripts/fastfarm_shared_ab.py`): tension means agree to 0.02 %,
  extremes to 0.25 %, and platform motions to millimetres. The migration guide now says that on
  a line with no `Fixed` end `FairTen` and `AnchTen` name the opposite ends from MoorDyn's.
- The power cable of `examples/iea15mw_umaine_mixed_cabledyn.dat` is meshed to the
  half-bending-length rule: 248 elements (2 m, 2.5 m and 2.5 m in its three sections) instead
  of 87 (4 m, 4 m and 10 m). On the 10 m grounded elements the touchdown was under-resolved and
  the static state, solved in the vertical plane of the cable, was unstable out of that plane.
  The static fairlead tension of the cable changes from 11 898 N to 11 922 N (11 921 N on a
  1 m mesh), its anchor tension from 1218 N to 1230 N, and the touchdown curvature from
  0.058 1/m to 0.079 1/m (0.080 1/m on the 1 m mesh). The curvature output moves with the
  mesh: `Curv4N20`/`BendMom4N20` (76 m from the hang-off) become `Curv4N38`/`BendMom4N38`
  (74 m, the peak of the sag bend). In the one-hour full-turbine run of the validation record
  the hang-off tension statistics are unchanged (11.5 ± 1.1 kN), the touchdown curvature peaks
  at 0.109 1/m instead of 0.075 1/m, and the anchor tension stays at its static 1.23 kN, where
  the coarse mesh showed stick-slip peaks up to 4.1 kN. The validation record and the tutorials
  are updated.

### Fixed

- The development environment (`environment.yml`) excludes OpenBLAS 0.3.34, which caused the
  access violation in the modal analysis tests listed under Known issues for 0.1.1. On AMD
  Zen 4 and Zen 5 CPUs that expose AVX-512, OpenBLAS 0.3.34 enlarges the GEMM blocking,
  and its Zen and Haswell DGEMM kernels then overrun their stack buffer for an inner
  dimension above 128. This happens in builds without the AVX-512 kernels, such as the
  conda-forge Windows package, and in any build when `OPENBLAS_CORETYPE` forces one of these
  kernels (OpenMathLib/OpenBLAS#6013 and #6021). The dense reference eigensolve
  (`DSYGV`) of the modal tests reached it. On hosted CI runners the fault appeared only on
  those runners whose AMD CPU exposed AVX-512. The solver itself calls only banded LAPACK
  routines and a 4 × 4 eigensolve, none of which passes DGEMM an inner dimension above 128,
  and the released Windows executables use reference LAPACK, so results are unchanged. An existing environment takes the exclusion with
  `conda env update -f environment.yml`, followed by a fresh build directory.
- The anchor tension of a finite-EI cable laid on the seabed no longer jumps when the line moves
  out of its plane. The seabed contact acts at the nodes, which carry no moment, so at an
  anchor resting on the floor the end moment of the end element's weight turned the free end
  tangent into the floor: the element dipped below the seabed between its nodes (27 mm for the
  power cable of `examples/iea15mw_umaine_mixed_cabledyn.dat` on its former 10 m grounded
  elements) and carried axial forces from −47 kN to +127 kN about a line tension of 1.2 kN.
  Once the hang-off swayed, the element buckled sideways; with seabed friction the anchor tension then locked at about 2.5 times its
  static value for the rest of the run (1.22 kN to 3.1 kN in a coupled OpenFAST run while the
  platform drifted toward the anchor), and without friction the end nodes wandered by up to
  0.18 m. The seabed now carries the end element's weight with its end moment at a held line end
  resting on the floor (along the floor normal on a sloped bathymetry), in the statics and the
  dynamics alike, in the share in which the element lies on the floor: the contact force of
  its inner node against half the element's weight, through a smooth step that is zero when
  that node is clear of the floor and one from half the weight up. Lines whose held end does
  not rest on the seabed, or whose end element rises from it (a taut line anchored on the
  floor), are unchanged. In the
  shipped finite-EI examples with an anchor on the seabed the end tensions change by at most
  0.34 N (0.004 %). Two decks of the statics robustness matrix change outcome. A lazy wave
  in a 1 m/s current on a frictional bed (`LW_1015__cur_fr`) had been reported as a hairpin
  past its anchor only at its exact seabed stiffness. It now fails closed: the friction-held
  touchdown goes into compression at about 0.33 m/s, before any tension-only equilibrium is
  reached from the still-water layout, and a slow current ramp in the dynamics with
  `tensile_safety` stops at the same point. A short line in a 1.9 m/s current on a
  frictionless bed (`CUR_1013`) had failed closed. It now finds the equilibrium with the
  grounded run turned sideways on the seabed: when the refinement of the catenary
  continuation stops on a compression in a current, its failure now names the current limit,
  so the driver retries out of the plane.
- A finite-EI line solved in 3D (restrained in torsion, or with an end connection direction or
  a current out of the vertical plane of its chord) on a deck mesh too coarse at touchdown or at a section
  change no longer stops with "the converged state is an unstable equilibrium". The stable
  equilibrium had read as a saddle: elements longer than the bending length cannot follow
  the bend, the axial force at their Gauss points oscillates into compression, and that
  compression lowered the out-of-plane stiffness until the stability test failed. Such a
  state is now refined, from the mesh the route returned, to the mesh whose elements are at
  most half the bending length: the catenary-continuation route is repeated there, as for a
  failed polish, and the state of the mesh-sequenced route is prolongated there and polished.
  This also covers a mesh that an adaptive or safety refinement left too coarse; a mesh that
  already meets the rule is not refined again. A note names both meshes. An equilibrium that
  is still unstable on the refined mesh stops as before. Decks that solved before give the
  same results. `tests/data/lazy_wave_touchdown_torsion_coarse.dat` is the regression case,
  solved with each route first.
- A finite-EI line over a bathymetry that is not a plane along its chord no longer fails its
  statics on a few centimetres of roughness. The catenary seed is laid on the plane fitted to
  the bathymetry along the chord, which need not pass through the anchor, and it was rejected
  as soon as a node dipped about 1e-6 of the line length below that plane. The
  catenary-continuation route then stopped at once, and the mesh-sequenced fallback often ended
  on a kinked branch ("unresolved/kinked equilibrium"). In a sweep of a chain on a rough
  seabed, 39 of 96 runs on a 25 m survey grid failed this way, at every roughness from 2 cm to
  0.5 m. The seed is now measured against the actual bed: a node deeper below it than its
  static penalty penetration (submerged weight over `kBot` × `Diam` × length) is lifted onto
  the bed, a shallower one is left to the contact solve, and the continuation, which queries
  the full surface, settles the contact. The equilibrium is audited, tested for stability and
  checked for kinks as before. A line too long for the fitted plane is also walked into place
  from a retreated anchor over such a bathymetry (the moved anchor kept at its height above
  the bed under it); as the bed can hold the excess length in its dips, a failed walk there
  is not reported as proof that the line is too long, and the mesh-sequenced route still runs.
  All 192 finite-EI runs of the sweep now solve, within 1 % of the fairlead force on a mesh
  four times finer and 0.8 % of the `EI = 0` line. Decks over a flat seabed or a planar
  bathymetry are unchanged, and so are all shipped examples and test decks. A deck over a
  non-planar bathymetry that solved before can now take the catenary-continuation route
  instead of the mesh-sequenced one, or start from a seed lifted out of a bump: in the sweep
  29 of 153 such runs changed, the fairlead force by at most 0.03 % and the anchor tension by
  at most 2.2 % (on 0.5 m roughness, where the anchor tension of the same deck differs by up
  to 13 % on a mesh four times finer). The decks `ROUGH_1001` to `ROUGH_1003` of the cable
  statics matrix (with mesh twins) are the regression cases.
- A line resting in a trench, a valley or a pit of a gridded bathymetry now has a static
  equilibrium. Along every interior grid line the bilinear surface has a crease, and a node on
  it took the normal of one adjacent cell, so a node that had to rest in a concave crease was
  pushed back across it from either side and had no equilibrium. An `EI = 0` chain whose
  frictionless tension falls to zero inside a trench (its excess length resting slack at the
  bottom) failed with "no admissible static equilibrium", and so did a grounded run on a
  two-dimensional rough seabed, which slides sideways into a valley along a grid line; a
  finite-EI line on such a seabed could fail its first dynamic step. Each crease is now
  rounded within a band of half-width 5 % of the smaller grid spacing next to its line, by
  the cubic fillet `J (d - |u|)^3 / (6 d^2)` of the slope jump `J` across the line, applied as
  a tensor product in both directions, so that the surface is twice continuously
  differentiable. It moves the surface by at most `|J| d / 6` along one line (1.7 cm for a
  slope jump of 0.2 on a 10 m grid; at a crossing of two lines the two roundings add), and a
  grid that is exactly planar, and every patch without a crease (one that would move the
  surface by less than 1 µm), keep the bilinear surface bit for bit; the rounding errors of a
  plane written to two or three decimals are themselves creases and are rounded by a fraction
  of a millimetre. In a band the elastic contact force is the gradient of the penalty
  potential of the normal penetration (up to the twist of a bilinear cell, which the tangent
  leaves out as before), with its Hessian as the tangent, in the statics and the dynamics
  alike. Every bathymetry query uses
  this surface, the Python `Bathymetry` and its seabed clearance included. A `Fixed` point at
  the tabulated depth of a creased grid node rests on the seabed and is placed on the rounded
  surface. In the statics campaign the `EI = 0` decks on bathymetry solve 40 of 43 (32
  before; the three left have anchors below the seabed), a 30 m trench now solves at every
  depth from 2 m to 30 m (5 of 11 before), with the anchor tension `w` times the trench depth
  as for a chain piled at the bottom (the fairlead tension within 0.15 % of the same trench
  sampled with a level bottom cell, on a mesh four times finer), and the `EI = 0` rough-seabed
  sweep solves 192 of 192 (191). Decks over a flat seabed or a planar bathymetry are unchanged, and so are all shipped
  examples and test decks. Decks over creased bathymetries that solved before change: the
  fairlead tension by at most 0.04 % (0.4 % for an `EI = 0` run that settles in another valley
  of a two-dimensional rough seabed), the anchor tension by up to 8.5 %. The rounding removes
  horizontal reactions of the one-sided normals: on a 2 m trench the two ends' horizontal
  forces now agree to 0.02 % (1.2 % before), and the anchor tension of a finite-EI line on a
  rough two-dimensional seabed is half as far from the one on a mesh four times finer. The
  decks `trench_touchdown_pile`, `trench_grounded_run_pile` and `rough2d_seabed` of the `EI = 0`
  statics robustness matrix are the regression cases, with an independent nodal force balance
  on the rounded surface. When the anchor walk of a finite-EI line over a non-planar bathymetry
  reaches a state that is rejected, the failure no longer states that the line is too long for
  its span to hang as a catenary, which is proven over a planar seabed only.

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
  address read or written. Previously a stack overflow in a Windows GNU build ended the process
  with no message. Every thread that enters an OpenMP parallel region keeps room to report its
  own stack overflow, and on Linux and macOS a hardware fault reaches the handler that was there
  before (the Fortran runtime's backtrace, or the default action and its core dump) with its
  original address and context, and a previous handler runs under its own flags and signal mask;
  a fault signal sent by another process (`kill`) is passed on as sent. A signal the driver was
  started with ignored (for example under `nohup`) stays ignored.
- Every non-zero exit the driver makes itself now ends stderr with the closing line
  `CableDyn_driver: ended with exit code <n>`. A process ended from outside (`taskkill /F`, *End
  task*, `kill -9`) runs none of its own code and cannot report anything; the missing closing
  line now identifies it, since `taskkill /F` leaves exit code 1, the same code as a refused
  input. The driver states this at start, after the banner, in an `Exit status:` line.
  `CableDynDriver.run`, `cabledyn-run` and `cabledyn-study` report such a run as ended early,
  with the time its output reached, instead of relaying the start-up log as the error; with a
  driver that does not write the `Exit status:` line (0.1.0 and earlier), a failure is reported
  with the driver's own text as before. The documentation explains how to recognise each kind of
  early end and warns that `taskkill /IM CableDyn_driver.exe /F` ends every CableDyn run on the
  computer. A run the Fortran runtime ends with its own error (for example when memory runs out)
  is reported as ended with a Fortran runtime error, quoting the runtime's line.
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
  were, and ending a cable module also clears its torsion frame. End connections set after
  the torsion, and a torsion director that is not the rigid connection direction, are refused
  by name.
- `END CONNECTIONS` row errors now name the line, the end and the offending value.

### Known issues

- An intermittent access violation in the modal analysis tests has been seen on some hosted
  CI runners with GNU builds that use OpenBLAS 0.3.34. The cause is under investigation. The
  released Windows driver uses reference LAPACK and is not affected.

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
