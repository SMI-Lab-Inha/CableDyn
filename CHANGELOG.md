<!-- SPDX-License-Identifier: Apache-2.0 -->

# Changelog

All notable changes to CableDyn are recorded in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Finite-EI lines with two moving ends on the cubic-Hermite route: a line whose End B is a
  `Coupled`/`Vessel` point (a cable shared by two floaters, a FAST.Farm shared line) is driven at
  both ends standalone, in OpenFAST, FAST.Farm and the C API (End B on a moving column of its
  own), with `END CONNECTIONS` and torsion at both ends (both end frames turn with their
  parents) and every Hermite feature. A `Coupled` End B held at rest writes the `Fixed` End B's
  files byte for byte; a symmetric cable whose ends move apart stays mirror-symmetric to 4e-11 m;
  a rigidly translating line keeps its
  end loads to 4e-5 N of 1.6e4 N; driving End B matches the reversed line driven at End A to
  1e-3 N, with a turning clamped End B to 2e-5 N and 4e-7 N m; equal rolls of both end frames
  keep the torque to 3e-6 N m and opposite rolls change it by `GJ (2 phi)/L` to 5e-6 N m; a
  FAST.Farm aggregate driven by prescribed turbine motions reports the standalone deck's end
  tensions to 7e-3 N of 1.9e4 N; restarts are bit-identical. The compatibility (Cosserat) route
  these decks ran on is retired: on four shared-line decks it failed to converge or reached a
  compressive static state (-6.1e8 N), and the Hermite route runs all four.
- A finite-EI line between two objects on the standalone multibody march: End B on a Rigid6
  `Body<N>` point or a rod end as well as End A, both ends driven by their objects and their end
  forces, connection moments and torques returned to them, with `END CONNECTIONS` and torsion at
  both ends (`bodyScheme monolithic`; the staggered scheme takes pinned ends). The monolithic step
  re-steps such a cable at every outer iteration and couples the junction coordinates of its two
  objects through its finite-difference end-to-end stiffness in one Newton matrix (on a stiff taut
  cable between light bodies, dropping that coupling took 9.2 mean outer iterations and 300 step
  halvings over 100 steps, against 3.9 and none); a pinned cable with one moving end is coupled at
  the predictor whichever end that is, and an end on a fixed rod's end is the `Fixed` end there
  (byte for byte). `bodyIC static` balances both objects on the cubic-Hermite cable, its passes
  accelerated by Anderson mixing: a stiff (EI 1e6 N m^2) clamped cable between two unequal moored
  bodies, out of their plane, converges in 8 passes, written in either direction, to the same
  equilibrium. Two identical bodies stay mirror images to the output precision, a body joined to
  an immobile body (1e15 kg) moves as on a `Fixed` End B to 8e-9, two free bodies keep the total
  linear momentum to 3e-5 and the energy to 0.2 %, a cable twisted 90 deg returns equal and
  opposite torques (`GJ Phi/L` at the start, their sum below 2e-10 N m), and a moored pair starts
  and stays at rest. An End B on a `Free`/`Connect` point, a `Point3` buoy or a `Coupled`/`Vessel`
  point stays refused there by name (this march drives no `Coupled` point), and so are, by name,
  torsion at an End B on a rod end and a line restrained in torsion at a rod End A with a moving
  End B, which no route of this build solves.
- C API host rotations (ABI 1.3): `CableDyn_InitDeckRotating` initialises a handle whose host
  passes each coupled point's rotation (`CableDyn_UpdateRotations`: the global-to-body DCM and
  the angular velocity and acceleration) with every update and step, so a line restrained in
  torsion and a `Rigid` end connection turn with the host, as on the coupled OpenFAST route; a
  torsion deck, refused by `CableDyn_InitDeck`/`InitDeckAt`, runs there.
  `CableDyn_GetCoupledMoments` returns the moments the lines apply to the host. On a straight
  line rolled at 60 deg/s the coupled torque is the quasi-static `GJ phi / L` to 3e-6 N m, and
  `SaveState`/`RestoreState` carry the rotations. Python: `CableDyn.init_deck_rotating`,
  `update_rotations`, `coupled_moments`.
- `ATTACHMENTS` on a finite-EI line restrained in torsion: the torsion stage of the static start
  runs on the discrete equilibrium, after the attachment polish, with the modules' net weights
  and current drag (to a fixed point) as nodal loads. A lazy wave with 25 buoyancy modules and
  one turn of pretwist, in still water and in a cross-flow current, keeps its torque, twist and
  fairlead tension over a still 10 s run to the output precision (a torsion stage before the
  polish would ring by 3e-3 of the torque); neutral modules reproduce the bare line's GJ Phi / L.
- `CoupledPinned`/`VesselPinned` bodies in OpenFAST, as MoorDyn-F's coupled pinned body: the
  body reference point is a node of the coupling mesh that takes the platform position only, the
  body turns about it with its own three rotations, and the node returns the body's net load
  less its inertia, with no moment. `PtfmInit` moves and turns the body, `bodyIC static` balances
  it about its node, and the checkpoint state and correction rewind carry it. A dry compound
  pendulum keeps the elliptic-integral period to 1.6e-4 at `dt = T/100`, a body with a
  `CoupledPinned` rod's mass and inertia follows the rod to 2e-12, and against MoorDyn-F
  (`inertialF 1`) a pendulum swinging through 120 deg agrees to 0.09 deg (0.15 deg when its node
  surges). Outside a coupled host, without inertias, and in linearisation the body fails closed
  by name.
- Active line control (`CONTROL`, OpenFAST `CableDeltaL`/`CableDeltaLdot`) of a Syrope
  polyester segment, which stopped at initialisation before. As in MoorDyn-F and MoorDyn-C the
  segment's slow stretch carries over the length change (its slow strain is read at the new
  length) and the running maximum is kept. On a taut 100 m line whose fairlead segment is paid
  out by 0.5 m, hauled in to -1 m and paid back, the tension follows the rest curve of its load
  history (working curve, OWC, then the working curve of the new maximum) to 2.5e-4 along the
  ramps and 6e-7 at rest; after a length step it is MoorDyn's at the carried stretch to 1e-15
  (carrying the strain instead would be off by 0.4 to 2.4 %); a restart in the ramp continues
  bit for bit; a polyester-chain line with a Syrope fairlead segment and a
  chain-polyester-chain line run.
- Active line control of a viscoelastic segment (`Es|Ed` and `Es|alphaMBL|vbeta`), which
  stopped at initialisation before. As in MoorDyn-F and MoorDyn-C the slow-branch stretch
  `dl_1` is left as it is and read at the new length. On the same taut line and cycle the
  tension follows the static law `Es` times the strain to 1.7e-4 along the ramps (4.6e-4
  load-dependent) and 4e-9 at rest; after a length step it is MoorDyn's at the carried stretch
  to 4e-16 (carrying the strain instead would be off by 0.4 to 4.6 %); a restart in the ramp
  continues bit for bit.
- OpenFAST linearisation of finite-EI lines restrained in torsion: the twist angles of the nodes
  the end frames do not hold, and their rates, are continuous states (`Cable 1 node 2 twist,
  rad`); the coupled point's orientation turns the torsion frame of the line's coupled
  (host-driven) end, its angular velocity and acceleration reach the torsional damping and the
  held end section's inertia, and the fairlead
  moment carries the torque. A straight line's torsional modes from the eigenvalues of `A` are
  the lumped chain's to 5e-14 (the continuum's `nπ/L √(GJ/Irn)` less the discretisation error,
  1.6e-4 for the first mode at 50 elements) and its roll stiffness `D − C A⁻¹ B` is `GJ/L` to
  1.5e-14; a 90° twisted lazy wave matches the deck's static solutions for a rolled and a moved
  fairlead to 4.9e-4. In a `Linearize = True` run every input Jacobian, those of the time march
  included, steps the coupled-point rotation by 1e-4 rad instead of 0.1 rad, which changes the
  rotation columns of a deck with a hang-off spring by up to 0.17 %.
- Prescribed motion and torsion on the standalone mixed deck (`EI = 0` and finite-EI lines between
  `Fixed`/`Coupled`/`Vessel` points, run on the coupled aggregate): the driver marches a
  `motionFile` or `vesselMotion` as on the cubic-Hermite route, row 1 installed before `t = 0`
  (also with `TMax = 0`). A vessel turns the end connection and torsion frames of the cables on
  its points, and the roll column imposes the twist of a line restrained in torsion, which these
  decks no longer refuse.
  Under vessel roll, pitch and heave, or a surge with the roll column, a clamped line gives the
  single-line route's torque, twist and node positions to the printed eight digits and its end
  tensions and bend moment within 8.2e-6 of their range (the coupled cables' Newton tolerance; with
  the single-line route's, every channel to the printed digits). These figures are for motions
  that start from rest; a start with a jump in velocity differs by several percent of the range
  (3 to 11 % in the cases tried) under the coupled Newton tolerance. A straight line rolled
  slowly gives `Torq = GJ Phi / L` to 5e-6 N m. A static mixed run (`TMax = 0`) with torsion
  takes `False alpha_force_blend`, as the single-line route and the Python reader do.
  `vesselRAO` stays refused with the deck waves. A line between two moving points there takes the
  roll columns of both ends, `Phi = Pretwist(B) − Pretwist(A) − roll(A) − roll(B)`: rolled
  oppositely (a rigid roll about its axis) a line twisted 90 deg keeps its twist exactly, and rolled
  equally it twists by twice the roll to 5e-7 deg, with `Torq = GJ Phi / L` at both ends to
  5e-6 N m.
- OpenFAST linearisation with the line dynamics, as MoorDyn-F's: the continuous states are the
  interior-node positions of every `EI = 0` line (MoorDyn's names and order) and the free
  degrees of freedom of every finite-EI cable (node positions and tangents), then their
  velocities, with `A`, `B`, `C` and `D` by central differences of the accelerations solved at the
  perturbed state. On the UMaine chains the stiffness `D - C A^-1 B` matches MoorDyn-F's to 0.1 %
  and the lowest nine line-mode frequencies to 0.3 %; with a lazy-wave power cable every coupled
  point matches the deck's own static solutions to 3e-4, also with a hang-off spring. The
  reduction is the static stiffness at an operating point whose coupled points do not
  accelerate: the coupled load includes the line's inertia, so at the accelerating `t = 0`
  platform the lazy wave's force per unit rotation reads 1.2 % off, and 2.7e-4 with the platform
  DOFs disabled (`doc/openfast.rst`). A cable end connection with a numeric stiffness now
  linearises. The viscoelastic slow-branch stretch
  and the Syrope slow strain are first-order states, so the reduction gives these ropes' relaxed
  static stiffness (each matches its deck's static solutions to 7e-5; held, they would read the
  chain 0.6 % and the polyester 2.9 times too stiff). Free/Connect points, a `Rigid` cable end
  connection, VIV lines, bodies, rods, buoys and `CONTROL` fail closed with named messages. This
  replaces the quasi-static linearisation without states.
- Composite Syrope lines: a Syrope polyester section between plain `EI = 0` sections, such as a
  chain-polyester-chain mooring, in one LINE. The static start puts every Syrope element on the
  rest curve of its load history with its own secant EA (the tension varies along a sagging
  section), so the dynamics start in the static equilibrium; the chains keep their `BA` damping,
  and `SYROPE IC` applies to every Syrope section of the line. On a 400 m polyester section
  between chains every element starts on its OWC to 2e-11, and a still run keeps the fairlead
  tension. Example `composite_chain_syrope_chain.dat`, on a moving vessel.
- Torsional inertia and damping of finite-EI lines: in the dynamics a line restrained in torsion
  carries a twist angle at every node, against directors that follow the node tangents by the
  smallest rotation over each step, with the polar rotary inertia `Irn` of its LINE TYPES row and
  the torsional damping `Bt` of a new third part of the `GJ` column (`GJ|K|Bt`, negative for a
  damping ratio). The torque now varies along the line (`Torq<L>N<J>` is the node's) and
  travels as a torsional wave; the statics stay condensed, and in the limit of a negligible
  inertia the dynamics reproduce the condensed torsion. The step solves the angles with the
  cable in one band of half-width 13; the twist state joins the snapshot, the recovering step
  and the OpenFAST checkpoint mirror. Decks with torsion whose `Irn` was a placeholder (the
  column was unused before) need the section's real value. The Python deck reader and builder
  take `GJ|K|Bt` (builder `gj=(GJ, K, Bt)`).
- Active line control (`CONTROL`) of a finite-EI cable restrained in torsion, which was refused:
  the torsional compliance, the nodal twist's `GJ`/l and `Bt`/l, the torque–tension coupling and
  the polar inertia take the new lengths, a `-ζ` torsional damping keeps its `ζ`, and the damping
  acts on the rate of the twist per length. A taut line paid out and hauled in carries
  `GJ Φ/L(t) + K (H − L)/L` to 1.5e-12 and its tension to 4e-6; a constant pay-out carries the
  Kelvin–Voigt torque to 6e-6 N m of its 0.64 N m damping share; a pay-out under a ringing twist
  balances the energy to 8e-5 of the work done; restarts mid-ramp (snapshot, checkpoint mirror,
  a coupled deck) are bit-identical. The checkpoint mirror of a line with the nodal twist gains
  one value, its committed torque.
- Modal analysis (`nModes`) of finite-EI lines restrained in torsion: the modes carry the twist
  angle of every node, as the dynamics do, with the lumped polar inertia and the second variation
  of the twist energy at the static torque. A clamped straight line's torsional modes match
  `nπ/L √(GJ/Irn)` within the lumped chain's own error (2.6e-4 on the first), elastic ends and
  the torque–tension coupling their continuum solutions to 1e-3, and under a torque its lateral
  modes match the clamped Greenhill cable to 3e-7; eliminating the angles reproduces the
  condensed static tangent to 1e-16. `<root>.modes.out` adds `PolarShare` and `dTwist` when a
  line carries twist.
- Torsion of a finite-EI line restrained at one end only, the other free to twist (`TorsStiffness`
  `Free` there, or no row): solved where a line restrained at both ends is and its LINE TYPES
  rows give `GJ`, instead of being ignored. The free end carries no torque and turns with the
  twist of the line: statically `Twist<L> = −Σ (K/GJ)(l − L)` under `GJ|K` (End B relative to
  End A; to 1e-8 on a stretched straight cable, the moving end free or the anchored one), and in
  the dynamics with torsional waves that reflect at it. A straight line held rigidly at the other end rings at the
  clamped–free frequencies `(2n − 1) π/(2L) √(GJ/Irn)` (3e-4 and 3e-3 for n = 1, 2), and a roll
  of the restrained end gives the continuum's free-end turn `A/cos kL` and reaction
  `GJ A k tan kL` to 0.5 % below and above the first mode, in the standalone driver and on the
  coupled aggregate. The `Torq`/`Twist` channels and range columns work on such a line, and
  its `CurvX`/`CurvY` axes follow its twist. Results change on such decks: the dynamics carry
  the twist (on a swaying lazy-wave cable with `GJ|K` the tensions move by up to 1.5e-4, 7e-8
  in a planar heave), and a line restrained at End A on a Rigid6 body joins the static body
  balance of the clamped ends: bending-Pinned there, its body's static moment now balances
  (Body1My −103.9 N m to 1e-6 on a buoy, FairTen1 126885 to 126857 N, 2.2e-4, Body1Ry
  0.3004 to 0.2984 deg; the march then differs by up to 2.6 % at dtM 0.05 s, 0.5 % at
  0.0125 s); clamped there, by about 1e-7. OpenFAST checkpoints of such decks written before
  this change do not reload: the cable's state mirror now carries the nodal twist, and the
  rebuild stops at a mirror-extent check (`... mirror extends past the state vector`). A line
  restrained at one end that misses a rule (no `GJ`, two moving ends, the configuration blend,
  a restrained End A on a free, pinned or Coupled rod end or a non-Rigid6 body) is solved
  without torsion as before, with a note naming the rule; a solved one's restrained `Pretwist`
  is held within 1000 turns.
- Torsion at a rod end: a line restrained in torsion may hang from an end of a `Fixed` rod (an
  I-tube or J-tube), whose frame holds the End A torsion frame, with `Nx` in the rod's deck
  axes. Two turns on an inclined J-tube give `GJ Phi / L` to 1.2e-9, the same torque as from a
  held Coupled point, and `Rod<N>M` the torque along the line axis; the line written anchor
  first (MoorDyn's order) runs bit for bit as written rod first. On a rod fixed to a rolling Rigid6
  body the torque follows `GJ (Phi - theta)/L` to 4e-11, and `Rod<N>M` of a rod fixed to a body
  now includes the connection moment of a cable clamped at its end (it read 0). Free, pinned
  and Coupled rods refuse the torsion by name, and so does `staggered bodyScheme` on a body or
  rod (it drove the cable end by translation, so the torsion frame did not turn); a line between
  two rods is refused as ending on a rod end, not as having two moving ends.
- The torsion routes combined: a line restrained in torsion at one end only takes `ATTACHMENTS`,
  `CONTROL`, `nModes` and a `Fixed` rod end, which left it unsolved in torsion. With modules of
  net weights 14.6, −9.1 and 5.8 kN on a vertical cable its free end turns by `−(K/GJ)(H − L)`
  to 8.5e-9 and its mid-line node with the module tensions to 3.9e-4 (a first-order mesh
  error); a paid-out fairlead turns it by `−(K/GJ)(H − L(t))` to 4e-11 with no torque; its modes
  and the eigenvalues of the OpenFAST `A` are the clamped–free chain's to 8e-10 and 2e-14 and
  the quarter-wave `(2n − 1)π/(2L) √(GJ/Irn)` within the lumped-inertia error, with no roll
  stiffness. The modal twist Hessian no longer takes a free end for an elastic one of infinite
  stiffness. `ATTACHMENTS` on a line restrained at both ends run on the mixed `EI = 0` +
  finite-EI route and on a `Fixed` rod end as on the single-line route (to the printed digits).
- Torque–tension coupling of finite-EI cables: the 14-column LINE TYPES row's `GJ` column takes
  `GJ|K`, with `K` the coupling of an armoured cable or wire rope (`T = EA ε + K τ`,
  `M = GJ τ + K ε`). The cable carries the axial stiffness `EA − K²/GJ` of a section free to twist.
  A line restrained in torsion adds the stretch's share of the twist, `−Σ (K/GJ)(l − L)`, to its
  condensed torque, in statics and dynamics and on every route that carries torsion. The node
  and element tensions, the element profile, the twist channels and the static and dynamic
  `tensile_safety` audits include the coupling, so a stretched cable that opposite-sense torsion
  puts into compression is caught. A stretched straight cable reproduces
  `T = EA ε`, `M = K ε` between rigid ends and `T = (EA − K²/GJ) ε` with a free end. `GJ|0` is the
  uncoupled cable bit for bit. The Python deck reader and builder take `GJ|K` (builder
  `gj=(GJ, K)`).
- Section-frame curvature channels `CurvX<L>N<J>` and `CurvY<L>N<J>` on finite-EI lines: the
  curvature vector on the cable's cross-section axes, for the fatigue of an armour wire or
  conductor at a fixed point of the section. The axes start at End A from the end connection
  (the torsion normal, or the static start's upper normal), turn with a clamped hang-off, are
  carried along the line by parallel transport and turn with the twist of a line restrained in
  torsion; `CurvX^2 + CurvY^2 = Curv^2`. On a line twisted 90 degrees the carried axes reach
  End B within 3e-6 rad of its pretwisted normal. `EI = 0` lines and finite-EI lines with two
  moving ends refuse the channels.
- Cross-flow vortex-induced vibration of `EI = 0` lines in the standalone driver: MoorDyn v2's
  LINE TYPES columns `Cl` or `Cl dF cF` switch on the synchronisation model of Thorsen et al.,
  with the same phase start and running-RMS window (`VIVrmsWindow`, default 500 steps of
  `dtM`). On a taut 100 m line in a uniform current it follows MoorDyn-F's lock-in amplitude
  and frequency to two to three digits through mode-1 and mode-3 lock-in, at any `dtM` from
  5e-4 to 1e-2 s. The phases start from the deck's End A as MoorDyn's do, also on a line
  written anchor first, which the reader turns around, and at the deck nodes of a refined
  finite-EI mesh. MoorDyn-C's own VIV case, the catenary riser of Thorsen et al., reproduces
  its mean tension and tension range and the measured lock-in frequency (period within 1 %).
  The Python deck tools read and keep the columns.
- Cross-flow VIV of finite-EI cables in the standalone driver: the same LINE TYPES columns on an
  `EI > 0` type give a distributed lift, Gauss-integrated over each cubic-Hermite element with
  exact tangents, and a shedding phase at each node driven along the node's own tangent. A taut
  line with `EI` near zero reproduces the `EI = 0` lock-in to three digits, at any `dtM` from
  5e-4 to 1e-2 s without smoothing; a bending-stiff line matches the string of equal first-mode
  frequency in mode-1 lock-in and stays in mode 1 where that string jumps to mode 2. The
  multibody march carries it too (a cable on a Rigid6 body or rod, or beside free points), with
  the same lock-in. A finite-EI line with two moving ends refuses `Cl > 0`.
- LINES Outputs flag `V` (MoorDyn's): the cross-flow VIV lift force at the nodes of a line,
  written to `<root>.Line<L>.V.out` on dynamic runs of the standalone routes that carry VIV. On
  a finite-EI line each element's lift is lumped to its nodes; at the interior nodes of a taut
  line in lock-in it matches the `EI = 0` line's within 0.1 %. Upper case only, as in MoorDyn,
  whose `v` (node velocities) CableDyn does not write.
- Cross-flow VIV on the coupled routes: OpenFAST (`CompMooring = 5`), FAST.Farm and the
  standalone runner of a mixed `EI = 0` / finite-EI deck carry the lift on both line kinds, with
  the host's fluid field. The phases, running mean squares and averaging time are part of the
  checkpoint state and of the correction rewind, so a restart continues bit for bit. A riser and
  a cable in a 0.15 m/s host current lock in at A/D 0.665 and 0.292 Hz, MoorDyn-F's 0.662 and
  0.293 Hz.
- Bend stiffeners (`STIFFENERS` table): a tapered polyurethane sleeve over an end of a finite-EI
  line, given as conical segments from the root outward (outer diameters, bore, modulus,
  density). Each element it covers takes the compliance average of the cable-plus-sleeve
  bending stiffness and the sleeve's mass, buoyancy and drag; standalone and coupled, with the
  Python reader and builder. A hang-off clamped through a 5 m stiffener follows an independent
  solution of the planar elastica within 0.2 % of the root moment, and its curvature follows
  OrcaFlex 11.6d within 0.8 %. Initialisation notes name a stiffener spanning fewer than 10
  elements and one at a pinned end. Example `lazy_wave_bend_stiffener.dat`.
- Output channels `StfMom<L>N<J>`, the share of the section moment a bend stiffener carries,
  and `FairDefl<L>`/`AnchDefl<L>`, the deflection of the end force from the direction a held
  END CONNECTIONS row declares (with the end tension, a stiffener's load case).
- MoorDyn-C time-varying and gridded currents on the standalone routes: `Currents 2` (the depth
  profile in time `current_profile_dynamic.txt`) and `Currents 5` (the grid in x, y, z and time
  `current_profile_4d.txt`), read and interpolated as MoorDyn-C v2.7.1 does (resampled at the
  smallest time interval, repeating past the last sample, held beyond the grid, no fluid
  acceleration) at every node's position and time, the static start included. At 739 points and
  times the current is identical to MoorDyn-C v2.7.1's; a constant table equals `Currents 1` and a
  grid uniform in x and y equals `Currents 2`; on a 70 m rope the tension change of a time-varying
  and of a sheared current matches MoorDyn-C's within 0.9 %, and the static start in a sheared grid
  holds the tensions to 1.4e-7 of the change the current makes. The C API, without a clock,
  refuses a table that varies in time.

### Changed

- The macOS wheel needs macOS 15 or newer on Apple silicon (0.3.0: macOS 14). It is built on
  the macOS 15 runner, as GitHub retires the macOS 14 runner image on 2 November 2026, and its
  Homebrew gfortran and OpenBLAS support no older macOS. On macOS 14, build from source.
- With a `motionFile` or `vesselMotion`, the single-line finite-EI route writes `.static.out`
  before it installs the first motion row, so the file holds the static solve; before, it showed
  the fairlead already moved to row 1. A deck whose row 1 is off the deck fairlead changes
  `.static.out` only (a 2 cm offset changed the fairlead x by 0.02 m and its tension from 3.5 MN
  back to the static 29 kN); the example decks are unchanged.
- Standalone `TURBINES` decks read the `END CONNECTIONS` direction `Ez` and torsion normal `Nx`
  of an end on a `Turbine<J>` point in that turbine's axes, turned by its `PtfmInit` rotation at
  the static start, as the FAST.Farm aggregate does (they were global). A deck whose turbines
  start rotated and whose finite-EI lines carry end connections or torsion changes; on a clamped
  line twisted 90 deg between turbines yawed +15 and -8 deg the standalone run now matches the
  FAST.Farm aggregate to 8e-9 in every channel. Decks without a `PtfmInit` rotation, or without
  end connections and torsion, are unchanged.
- Standalone `TURBINES` motion rows turn the `END CONNECTIONS` and torsion frames of a finite-EI
  line end on a `Turbine<J>` point with the row's rotation, angular velocity and acceleration (End A
  and End B), as `vesselMotion` and the OpenFAST and FAST.Farm hosts do; they were held at their
  start orientation while the end moved. A `TURBINES` deck whose rows rotate the turbines and whose
  finite-EI lines carry end connections or torsion changes: a clamped line twisted 90 deg between
  turbines that roll 10 deg in opposite senses over 1 s ends at 613.6 N m of torque (788.1 N m
  before; the tensions change by under 5e-7) and matches the FAST.Farm aggregate to 3.4e-8 in every channel.
  No example deck uses `TURBINES` rows; the examples are unchanged.
- `TMax = 0` checks the static start of any dynamic deck as it stands: a `motionFile` or
  `vesselMotion` record may run past `TMax` (its later rows are checked on the `dtM` grid and
  not used, so a long record also serves a shorter run), and `RangeStart` is accepted at
  `TMax = 0`, where the range graphs hold the `t = 0` row. A standalone deck without `TMax` that
  asks for the march now says so: its error ends "to check a dynamic deck's static start only,
  give dtM with TMax = 0".
- Deck waves from a WaterKin elevation history (`WaveKinMod 1`) run two to four times faster.
  Below the surface, the velocity and acceleration sums stop once the components still to come,
  too short to reach that depth, cannot change them beyond round-off. On a 1024-element lazy-wave
  cable over 300 s, a record at `dtWave` 0.25 s (600 components) runs in 103 s instead of 226 s,
  and one at 0.1 s in 189 s instead of about 660 s, with output identical to the printed digits.
  A `jonswap` sea changes at round-off.
- An elevation record (WaterKin `WaveKinMod 1`) keeps only its Fourier components up to the new
  `WvHiCOff` option, in rad/s as in SeaState, by default 0.5 Hz (MoorDyn-C's fixed cut for
  `WaveKin 3`, which now follows the option too); `none WvHiCOff` keeps them up to the Nyquist
  frequency of `dtWave`, as MoorDyn-F and earlier CableDyn did. A one-hour record at `dtWave`
  0.1 s drops from 18,000 components to 1,800. On the 300 s lazy-wave case both records then run
  in 72 s, and every output channel's standard deviation moves by at most 0.002 %. A record with
  no component at or below the cut-off is refused by name.
- Coupled OpenFAST runs start from the equilibrium in SeaState's steady current (`CurrMod 1`),
  as standalone decks start in their `current` option, instead of from the still-water shape.
  The current is evaluated at SeaState's vertical grid levels for the static solves only; the
  dynamics take the SeaState field alone. On the IEA-15 MW mixed example in a 0.5 m/s cross
  current, the lazy-wave cable starts 29 m across at its deepest monitored node and moves less
  than 0.6 m over the first 30 s, where the still-water start has reached a third of that offset.
  `CurrMod 2` and InflowWind currents keep the still-water start, with a note.
- Syrope polyester lines take the seabed, Morison drag and added mass, deck current and waves,
  a coupled host's fluid field, a `bathymetryFile` and VIV, as every `EI = 0` line does and as
  MoorDyn's do; before, they ran without hydrodynamics and refused current, waves, host fluid
  and bathymetry. The static initial condition includes the current's drag. A taut polyester
  line whose fairlead sways 2 m in water matches MoorDyn-C's fairlead tension within 11 N on a
  45 kN range (1.03 MN mean); without drag and added mass the range is 5 % larger.
- Under a bend stiffener, `BendMom<L>N<J>` (and the `BendMoment` columns of the static profile,
  range graphs and element file) report the cable's own moment, the line type's EI times the
  curvature; lines without stiffeners are unchanged.
- A run whose outputs name line nodes (`L<L>N<J>…`, `Ten<L>N<J>` and the other node channels, or a
  LINES `Outputs` flag) notes that CableDyn numbers line nodes from 1 at End A while MoorDyn
  numbers them from 0, so a MoorDyn deck's `N<J>` is CableDyn's `N<J+1>`, or `N<NumSegs+1-J>` on
  a line written anchor-first. The numbering itself is unchanged; the outputs and migration pages
  state it.
- The run banner, README, manual, citation and package descriptions call CableDyn what it is: a
  geometrically exact finite-element solver for the statics and dynamics of mooring lines and
  dynamic power cables (its strain and curvature are the exact large-deformation measures), in
  place of "geometrically nonlinear cable & mooring dynamics".

### Fixed

- A rod numbered as a body no longer joins the body's attachment points. A `Rod<N>A/B` end point
  keeps its rod's id where a `Body<N>` point keeps its body's, and the body's point list took
  every point with its id: a rod end POINT row listed before the body's points became the body's
  first point, and the list overran its allocation (also in 0.3.0). A free body on three legs with
  a pinned arm rod numbered 1, its tether on an explicit `Rod1B` POINT listed first, now writes
  the outputs of the same deck with the arm as rod 3 bit for bit; before, the arm turned 172 deg
  in pitch, the body 1.06 deg, and a leg tension moved 3.0 kN in 1 s.
- The static start of a neutral finite-EI line whose length equals its span (straight at zero
  strain, so with neither weight nor tension) failed with neutral `ATTACHMENTS` or with
  attachments that only carry a drag area: the attachment polish stalled at the round-off of the
  axial forces, 7.6e-8 N, against a tolerance of 1e-8 times the solver's 1 N force floor, at
  every load step. For a line whose chord is at least its length the force scale of the
  convergence test now includes that round-off floor over the tolerance
  (`16 eps max(EA/L0) |x| / tol`, as on the `EI = 0` path), so these lines start straight and
  unloaded, also between two moving ends, restrained in torsion at one or both ends, and with a
  chord 1e-8 longer than the line (its 0.1 N tension to 2.2e-7 N); slack lines and lines that
  carry weight or tension above the floor are unchanged.
- The Newton tangent of the finite-EI static solve in a depth-profile current or a current grid
  carries the current's gradient: the drag samples the current at the node positions, so it moves
  with them. The drag part of the tangent now matches central differences of the residual to
  8e-10 (before off by 5e-2 to 0.2 in sheared profiles). Newton takes up to a third fewer
  iterations on a light cable in a strong shear (16 instead of 23); the equilibria are unchanged
  within the solver tolerance.
- The static start of an `EI = 0` line held above the still-water level kept the buoyancy of its
  dry part, while the march (with `WtrDpth` or a bathymetry) gives that part its dry weight. The
  start was not the march's equilibrium: a chain hung from 5 m above the surface onto the
  seabed moved by 1e-3 of its fairlead load in the first step, and the coupled start differed
  from the standalone one by 2.4e-4. The static solve now carries the dry weight (a fixed point
  on the wetted lengths): the start holds to 1e-10, the two routes agree to 1.4e-8, and a
  vertical chain's end loads sum to its weight to 1.5e-8 (before 1.4 % low). The force balance
  of `Connect` and `Free` points and free bodies and rods takes the same dry weight: a 50 t
  `Connect` point hung by two chains from 10 m above the surface holds its position to the
  output's 8 digits over 2 s (before it sank 1.1e-2 m). There a finite-EI line, carried without
  its bending stiffness, takes its dry weight with or without a seabed, as its cubic-Hermite
  cable does: a stiff body with a light cable pinned 5 m above the surface on a deck without a
  seabed holds its heave to 7e-6 m over 3 s (before it moved 2.8e-3 m), and an `EI = 0` chain
  between that body 5 m above the surface and a `Connect` clump holds both to the output's 8
  digits (before 5.5e-3 m). Lines below the surface are unchanged.
- `ATTACHMENTS` above the still-water level kept their buoyancy in the static start of a
  finite-EI line, while the march loads them with their weight alone. A 500 kg, 0.2 m^3 module
  3 m above the surface on a cable hung from 5 m above it moved the fairlead load by up to
  6.1 % in the first second, and as much on the same cable restrained in torsion. The
  attachment polish and the torsion stage after it now give each attachment the load of its side
  of the surface (a fixed point that repeats the solve until none changes side); both starts
  hold to the output's 8 digits. An attachment that its weight in air takes below the surface
  and its weight in water above it (a net-buoyant module that would float) has no static
  position and fails closed, naming the attachment that crossed the surface.
- A bending-pinned finite-EI line with `ATTACHMENTS` from a Rigid6 body or a rod end left its
  attachments out of the body static balance, which carries such a line as a cable without
  bending stiffness: on a moored body with a 2 t clump 50 m along a pinned cable the body started
  at the pose of the line without the clump and drifted by up to 2.9e-2 m and 0.20 deg in 5 s.
  Such a line now takes the static passes of a connected end, which balance the body on the
  cubic-Hermite cable's end loads; the body holds its pose to the output's 8 digits. Their
  passes are Anderson-mixed, as those of a line between two objects, and the core cable of the
  passes carries the line's `ATTACHMENTS` smeared, so the passes correct only the lumped,
  side-of-the-surface and bending differences: on a free spar rod with a 1500 kg clump on a cable
  pinned at its upper end the plain passes on a core cable without the clump diverged (the deck
  stopped), and Anderson-mixed ones converged in 26 passes or stopped by round-off (a 1-ulp
  change of an anchor coordinate; the Linux build stopped). They now converge in 6, the same
  for those 1-ulp changes, and the rod holds its pose to 8 digits. The mixing drops history
  columns that add less than 1e-6 of a new direction (1e-10 before), which had let round-off
  columns drive the mixed loads to 1e6 N. Across a sweep of the rod deck (clump 0 to 4 t, cable
  50 to 70 m, rod mass, clump position, start offsets, 1-ulp changes) 42 of 45 decks converge
  (34 before), the three left with no stable equilibrium found on either core cable.
- A bending-pinned finite-EI line on a Rigid6 body or a rod end without `ATTACHMENTS` was left to
  the body static solve's cable without bending stiffness, whose end force lacks the bending's
  share: a free spar rod with a 55 m pinned cable (EI 2e4 N m^2) started 5.7e-2 m off its
  equilibrium, drifted 7e-2 m in 3 s and moved a leg tension by 2.3 %. Such a line now takes the
  static passes too (4 here); the rod holds its pose to the output's 8 digits.
- In OpenFAST, Free Rigid6 bodies and Free, Pinned and CoupledPinned rods started their first
  step from zero acceleration: the module restores its state mirror after the output call at
  `t = 0`, which marked the objects' accelerations as committed before the start accelerations
  were solved. They are now solved when the module initialises (and, before the first step, at a
  committed step-start transfer and in new host fluid fields, so a SeaState current or wave starts
  them in that field), and the line ends on a CoupledPinned object move with them.
  An object in its static equilibrium changes at round-off (1e-12 m and 1e-6 deg over 10 s); one
  started out of balance (`bodyIC deck`) changes by the first step's velocity error, 0.02 m/s and
  0.05 deg for a hanging Pinned rod and up to 0.6 m/s for bodies and rods released far out of
  balance. The committed host kinematics at the first step are those of the initial state, so a
  host acceleration at `t = 0` enters neither version.
- MoorDyn-F's `inertialF` option was an unknown keyword; `inertialF 1` (CableDyn's behaviour)
  is now accepted, and `0` and `2` are rejected by name.
- A Syrope line held at its running maximum (the fairlead pulled out and held) under
  `False alpha_force_blend` stopped with "Syrope backward-Euler state solve did not converge":
  next to the top sample of the working curve the Newton step's difference quotient straddled
  the kink and crept along one side of the bracket. The solve now finishes by bisection there,
  and where a light rope's steep residual (BA_s|BA_d = 1e4|1e3) changes sign between two
  adjacent reals it returns that root to machine precision. Over 440,000 fuzzed states every
  solve that converged before is bit-identical, and the 26 that failed now converge.
- The static start of an `EI = 0` line now carries a deck current that is zero along the
  straight chord between its ends but not where the line sags (a depth profile or a current
  grid): it was solved in still water and the march drifted toward the current's equilibrium.
  On a 70 m rope with its ends at z = -20 m in 0.8 m/s below z = -32 m the start now holds to
  machine precision (it drifted 192 N in 4 s) and agrees with the finite-EI route within 0.14 %.
  The static tangent takes the current's gradient along the line, so the Newton iteration also
  converges in a current sheared across the line. A later `current` row now discards the levels
  of an earlier profile row.
- The modal analysis of finite-EI lines ignored the end connections and the discrete
  attachments: a clamped end was solved pinned (the first twelve frequencies of the lazy wave
  of `torsion_lazy_wave_hangoff_twist.dat`, untwisted, rise by 0.1 % to 3.8 % with its clamps)
  and buoyancy modules had no mass (with `nModes`, the first six frequencies of
  `lazy_wave_buoyancy_modules.dat` drop by 13 % to 29 %). The positions and translations of a
  line with a heading are written in the deck's axes, not the line's.
- A finite-EI line restrained in torsion whose bending-pinned End A is on a body or a rod is
  refused on `bodyScheme staggered`, as a clamped End A is: that scheme passes cable ends as point
  forces and held the torsion frame still while the body turned (a body could spin through +-180
  deg); `bodyScheme monolithic` turns the frame with the body.
- A deck with a WaterKin elevation history (`WaveKinMod 1`) at `TMax = 0`, a check of the static
  start, reduces the record over its own length instead of failing for want of a positive
  `TMax`.
- Reading a WaterKin `WaveKinFile` no longer slows with the square of its length: the rows grow
  by doubling instead of copying the whole record per row. A 10 s run on a 168,000-row record
  took 107 s and now takes 3.5 s.
- A weightless or near-neutral taut finite-EI line (a neutrally buoyant cable, a tendon) failed
  its static solve, even straight along its chord: the convergence test scaled the residual by
  the line's weight alone, far below the round-off of its axial forces. A taut line's scale now
  includes the tension its stretched chord implies; slack lines are unchanged.
- The anchor force of a finite-EI line resting on a sloped seabed (`bathymetryFile`). The floor
  took the whole downward part of the end force, so a cable laid down the slope from its anchor
  lost the down-slope pull: on a 10° slope `AnchTen` read 1.5 % low and `AnchFz`/`Point<P>Fz`
  read 0. The floor now takes only the part pressing into it along its normal, and on a
  frictional seabed the along-floor part of the end node's weight up to the friction it can
  hold; flat seabeds are unchanged.
- OpenFAST's native `<root>.CD.out` wrote a channel value below 1e-99 or from 1e100 on as a
  field of asterisks (`ES15.6E2`): the round-off noise of a zero component, such as the y
  position or tangent azimuth of a line in the x-z plane, broke every reader. Channel values now
  have a three-digit exponent (`ES15.6E3`, same width). The mixed `EI = 0` + finite-EI `.out`,
  the modal results, the static driver's mesh file and the OpenFAST module driver, which wrote
  such values without the exponent letter, do too.
- `openfast -restart` with a finite-EI cable continues the run bit for bit. The rebuild moved the
  cable's fairlead to the restored input mesh, whose last bits differ from the boundary the
  uninterrupted run committed, and the cable's `.CD.out` channels drifted by up to 1e-6
  (relative) within seconds. The checkpoint now carries the cable's committed boundary.
- A finite-EI line longer than its span with no net weight, or too little for its bending
  stiffness (0.001 N/m on a 120 m line with `EI` 1e4 N m²), failed its static solve: both routes
  reached the compressed straight strut or stalled at the round-off of the axial forces. Such a
  line is now solved heavier and carried back to its own weight in the vertical plane of its
  chord. A weightless line bows as the Euler elastica (end thrust and bow to 3e-5 of the exact
  solution); nothing in the deck chooses the side of its bow, so it is solved bowed downward and
  a note says so. With its ends one above the other it is still refused, by name.
- A line along +X reported its tangent azimuth (`L<L>N<J>Azi`, the static profile's `Azimuth`)
  at 0 or about 360 from one row to the next as it swayed a fraction of a degree to either
  side, or by the sign of a round-off component: in OpenFAST about half the rows of a cable in
  the x-z plane read 359.9x. A line whose chord heads within 45° of +X now reports its azimuths
  in (−180, 180], every other line in [0, 360) as before.

## [0.3.0] - 2026-10-23

### Added

- `pip install cabledyn`: the package is on PyPI, with wheels for Windows x64 and Linux x86-64
  (processors with AVX2) and macOS on Apple silicon that carry the solver library and the
  standalone driver, so it runs decks without the release executables; other platforms build
  from the source distribution. The release wheels are built and smoke-tested in CI and published
  through PyPI Trusted Publishing.
- OpenFAST with `dtM` larger than `DT`: between CableDyn solves the coupled loads are extrapolated
  from the last two solves instead of held (deck `OPTIONS predictLoads`, default `True`). On the
  IEA-15 MW VolturnUS-S at `dtM = 0.1 s` and `DT = 0.025 s` the platform motions differ from a
  `dtM = DT` run by at most 0.12 % of their standard deviation, against 0.42 % held. Runs at
  `dtM = DT` (the default) do not change.
- Initialisation notes, none of which changes results: a `dtM` that resolves a line's axial
  period in fewer than six steps, with the `dtM` that resolves it; and a finite-EI line whose
  grounded elements lock more seabed friction than 15 % of the touchdown tension, with the
  grounded element length that resolves the tension of the laid part.
- `validation/scripts/openfast_supercycle_ab.py`, which scores coupled runs at a coarser `dtM`
  against `dtM = DT` (`VALIDATION.md`, supercycling verification).
- The run-time profile (`CABLEDYN_PROFILE=1`) splits the line and cable steps into their kernels
  and assembly sections, and the OpenFAST `CalcOutput` and SeaState sampling into their parts.

### Changed

- Speed: the Windows release executables are compiled for processors with AVX2 (x86-64-v3:
  Intel Haswell and AMD Excavator or newer), which makes the driver 13-33 % faster on the
  benchmark decks; on an older processor they stop at start with a message naming the missing
  instructions. In OpenFAST, `CalcOutput` no longer evaluates cable reactions it discarded and
  assembles a fairlead reaction from the elements at its node (a 600 s IEA-15 MW run with a
  248-element power cable: 227 s to 194 s). Results change at round-off only.
- A static solve in a current that seabed friction cannot hold fails closed sooner (22 s to
  17.7 s on a 1000 m lazy wave), with the same message.
- The manual explains the mesh dependence of the negative `BA` (damping ratio) form and of
  node-lumped seabed friction, and the measured trade-off of `modified_newton`, which stays off.

### Fixed

- Finite-EI node tension, curvature and bend moment at touchdown understated the last suspended
  node's tension by half an element's floor reaction (5.3 % on a 10 m lazy-wave mesh); the coarse
  mesh now follows a four times finer one within 3 %. End forces and the solution do not change.
- OpenFAST restarts with `dtM` larger than `DT` now continue the run bit for bit (the checkpoint
  carries the held SeaState fields). A checkpoint written by an earlier CableDyn cannot be
  restarted with this one.
- A restore (`CableDyn_RestoreState`, OpenFAST correction iterations) of a finite-EI cable with
  seabed friction now continues bit for bit like the step it repeats.
- The FAQ and the OpenFAST tutorial gave 0.1 s as the default `dtM` in OpenFAST; it is the glue
  `DT`.

## [0.2.0] - 2026-10-20

### Added

- Seabed friction in the static solve toward the as-laid position (`staticsFriction asLaid`,
  with `layAzimuth` and `asLaidTension`), as in OrcaFlex. The default stays `none`, so existing
  decks run as before. A 700 m chain on a 5° slope that the frictionless statics folds 24 m past
  its anchor now holds where it was laid.
- Active tensioning (`CONTROL`) of finite-EI cables in OpenFAST, alone or in a deck with
  `EI = 0` lines.
- Sub-stepping of tension-only (`EI = 0`) lines through fast segment rotation and slack-to-taut
  snaps, on by default (`ei0_substeps auto|off`, `ei0_substep_angle`, `ei0_substep_max`).
  It brings snap-load peaks at coarse time steps closer to converged runs; `off` restores the
  previous results.
- Torsion of finite-EI lines in coupled OpenFAST (`CompMooring = 5`) and FAST.Farm runs.
- `Pinned` and `CoupledPinned`/`VesselPinned` rods in coupled OpenFAST runs.
- Seabed contact for `Free`/`Connect` points and `Point3` buoys: a heavy clump no longer sinks
  through the seabed.
- Anchor uplift channels `AnchFz<L>`, `AnchFH<L>` and `TDP<L>Grd` for drag-anchor checks.
- C ABI 1 minor extension 2 for hosts that drive `Coupled` points (for example CFD solvers):
  `CableDyn_InitDeckAt`, `CableDyn_SaveState`/`CableDyn_RestoreState` and
  `CableDyn_GetCoupledLayout`, with matching Python methods. Existing functions are unchanged.
- `OutputRecovery equilibrium|element` option and the `CD_HermiteCable_Dyn_Curvature_Exact` and
  `CD_HFMF_Curvature_Exact` functions (see Changed).
- Run-time profile of the standalone driver (`CABLEDYN_PROFILE=1`).
- New diagnostics, none of which changes results:
  - run audit of `EI = 0` lines in the standalone driver: tension spikes, fast segment turns
    and unresolved snaps, and a note for `rhoInf` above 0.95;
  - a note when the static state of a finite-EI line is unstable out of its vertical plane;
  - a note when the static layout of an `EI = 0` line folds back on the seabed, naming the
    cause and the remedy;
  - notes for a line stretched far beyond its elastic range and for a section too coarse for
    its line;
  - a warning when line nodes of a coupled OpenFAST run lie outside the SeaState
    wave-kinematics grid.

### Changed

- Finite-EI node tension, curvature and bend moment (and `.t.out` segment tensions) are now
  recovered from nodal equilibrium. On the lazy-wave example the largest tension error falls from
  1552 N to 60 N on its former 87-element mesh, and oscillations near touchdown on a coarse mesh
  are removed. End tensions and positions are unchanged. `OutputRecovery element` restores the
  previous definitions.
- `EI = 0` lines step with a generalised-α force blend, as finite-EI lines do. Tension spikes of
  lines that go slack are gone: a light cable under 3 m heave reported 4.2 MN where converged
  runs give 0.093 MN. Tension statistics of the shipped examples change by at most 0.23 % (mean)
  and 0.9 % (maximum). `False alpha_force_blend` selects the former scheme.
- Active tensioning of `EI = 0` lines follows MoorDyn-F through length ramps: fairlead tensions
  within 0.13 % instead of up to 6.3 %. A command outside MoorDyn-F's bounds stops the run
  instead of printing a warning, and the command acts on the segment at the host-driven
  (fairlead-side) end, as in MoorDyn-F.
- In OpenFAST, a deck without `dtM` now solves every glue step (`dtM = DT`) instead of about
  every 0.1 s, so the mooring loads reach the platform without delay. A deck `dtM` that is a
  multiple of `DT` still trades this for speed (see Known issues).
- The touchdown point (`TDP<L>`) now moves continuously instead of in steps of one element.
- Finite-EI statics refine the mesh locally instead of as a whole, which cuts solve time (for
  example 10.0 s to 0.33 s on a 462-element lazy wave). Decks that do not refine are unchanged.
- The power cable of `examples/iea15mw_umaine_mixed_cabledyn.dat` is re-meshed from 87 to 248
  elements. Its static fairlead tension changes by 0.2 %.
- Faster runs: directional and multi-train wave kinematics (17 %), finite-EI cables clamped to
  bodies or rods (66 %), and decks of `EI = 0` lines (13 %). Results can differ in the last digit
  with the Intel compiler.
- The OpenFAST linearisation check and the FAST.Farm shared-mooring comparison against MoorDyn
  are extended: per-fairlead stiffness agrees to 0.1 % (0.29 % for the transverse term) and
  shared-mooring tension means to 0.02 %.

### Fixed

- Node output channels of a finite-EI line whose mesh the static solve refined now read the
  deck node named, not the same index on the refined mesh.
- Static solves that failed or gave the wrong answer:
  - creased gridded bathymetry (trenches, valleys) now has a static equilibrium;
  - a finite-EI line over a rough, non-planar bathymetry no longer fails on a few centimetres
    of roughness;
  - a 3D finite-EI solve on a coarse mesh no longer stops with "unstable equilibrium";
  - a finite-EI line with `ATTACHMENTS` in a cross current on a frictional seabed now solves.
- The anchor tension of a finite-EI cable laid on the seabed no longer jumps when the line
  moves out of its plane.
- Explicit sub-stepping of free and pinned rods and Rigid6 bodies now accounts for seabed
  contact, free-surface crossing and decks without lines, and seabed friction of rods and
  bodies is applied implicitly (no more chatter on a frictional bed).
- `modified_newton` on `EI = 0` lines rebuilds the tangent as documented.
- Junction end tangent of a tension-only line step that ends right after a fresh tangent
  evaluation.
- Error messages are no longer truncated, and now name the cause and workaround for neutral
  buoyancy, unsupported couplings and anchors below the seabed.
- The development environment excludes OpenBLAS 0.3.34, which crashed the modal tests on some
  AMD CPUs. Results are unchanged.

### Known issues

- In OpenFAST, a deck `dtM` larger than `DT` holds the loads between CableDyn's steps, so they
  reach the platform about `dtM/2` late. The lag matters for stiff (taut or tension-leg)
  moorings and for tension statistics; check such results at a smaller `dtM`. See the note on
  the manual's OpenFAST page.

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

[Unreleased]: https://github.com/SMI-Lab-Inha/CableDyn/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/SMI-Lab-Inha/CableDyn/releases/tag/v0.3.0
[0.2.0]: https://github.com/SMI-Lab-Inha/CableDyn/releases/tag/v0.2.0
[0.1.1]: https://github.com/SMI-Lab-Inha/CableDyn/releases/tag/v0.1.1
[0.1.0]: https://github.com/SMI-Lab-Inha/CableDyn/releases/tag/v0.1.0
