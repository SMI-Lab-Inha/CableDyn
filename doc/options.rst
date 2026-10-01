.. SPDX-License-Identifier: Apache-2.0

OPTIONS reference and defaults
==============================

This page lists every keyword the ``OPTIONS`` section accepts, with every alias. Each keyword is
described with the same fields: value type and unit, default, valid values, the route that
consumes it, and an example record. The runnable :file:`examples/cabledyn_options_reference.dat`
shows every syntax form, including commented mutually exclusive alternatives. Every maintained
example repeats the common defaults that affect its route and links back to this page. A value
labelled *absent* is not represented by a magic number: omission itself is meaningful.

Record syntax
-------------

* **Scalar records** use the OpenFAST order ``value keyword``. Anything after the keyword is
  commentary, so ``0.05 dtM - Time step (s)`` and ``0.05 dtM time step`` are both valid.
* **Positional records** (``current``, ``waves``) put the model name first and the keyword last:
  ``airy 2.0 8.0 0.0 waves``. The ``dynamic_solver`` record puts its keyword first.
* **Descriptions.** In every OPTIONS record, the first ``-`` that has whitespace on both sides
  starts a description; it and everything after it are ignored. Negative numbers such as ``-1.0``
  are unaffected.
* **Keywords** are case-insensitive. An unknown keyword is fatal
  (``unknown OPTION keyword "<name>"``); CableDyn never silently ignores a misspelling.
* **Repeated keywords.** Records are applied in order, so the last record for a keyword wins
  (for example, a subsequent ``0 motionFile`` row disables an earlier path).
* **Numbers.** Every numeric value must be one plain number with a ``.`` decimal point
  (``1e-3``, ``0.05``, ``-50``). A value containing ``/``, ``,`` or ``;``, or written as a repeat
  count such as ``2*0``, is rejected with the deck line number instead of being read partially;
  ``1/20 dtM`` is an error, never ``dtM = 1``. Integer-valued options (for example
  ``recovery_max_substeps``, the quadrature orders, ``nModes``, ``StreamOrder`` and ``WaveSeed``)
  reject non-integral values.
* **Logical values** (for example ``modified_newton``, ``adaptive_mesh``, ``alpha_force_blend``
  and ``cable_load_feedback``) accept,
  case-insensitively, ``True``/``T``/``yes``/``y``/``on``/``1`` and
  ``False``/``F``/``no``/``n``/``off``/``0``. Any other spelling is an error.
* **File paths** are single whitespace-free tokens, optionally quoted, resolved relative to the
  deck's directory. They may contain ``/`` or ``\`` separators and may be longer than the
  64-character identifier limit. A ``#`` or ``!`` starts a comment even inside a path.

In the tables, *Standalone* means ``CableDyn_driver.exe``; *OpenFAST* means a ``CompMooring = 5``
(or ``Mod_SharedMooring = 5``) coupled run in OpenFAST, maintained by NLR (National Laboratory of
the Rockies, formerly NREL).

Environment and seabed
----------------------

.. list-table::
   :header-rows: 1
   :widths: 16 12 11 31 18 12

   * - Keyword and aliases
     - Type, unit
     - Default
     - Valid values and rules
     - Route
     - Example
   * - ``g`` / ``gravity``
     - real, m/s²
     - ``9.80665``
     - finite, > 0
     - both; OpenFAST host ``Gravity`` overrides the deck
     - ``9.80665 g``
   * - ``rhoW`` / ``WtrDnsty`` / ``water_density``
     - real, kg/m³
     - ``1025``
     - finite, > 0; also used by ``EQUIVALENT BUOYANCY`` conversion
     - both; OpenFAST host ``WtrDens`` overrides the deck
     - ``1025.0 rhoW``
   * - ``WtrDpth`` / ``water_depth``
     - real, m
     - absent: no flat seabed
     - finite, > 0; flat seabed at ``z = -WtrDpth``. Mutually exclusive with a
       ``bathymetryFile`` row. Required by ``waves``
     - both; OpenFAST uses host ``WtrDpth`` unless the deck names a ``bathymetryFile``
     - ``200.0 WtrDpth``
   * - ``bathymetryFile`` / ``bathymetry_file`` / ``seafloorFile`` / ``seafloor_file``
     - path
     - absent
     - non-empty. The file holds ``x y depth`` rows (plain numbers, any order, comments allowed)
       forming one complete rectangular grid of at least 2 × 2 points with no duplicate
       ``x y`` pair; depths finite and > 0
     - standalone static ``EI = 0``, independent-line ``EI = 0`` dynamic, ``Connect``/``Free``
       point-system, Hermite finite-EI dynamic, free/fixed rod, and Rigid6 decks (not Syrope
       lines); in OpenFAST the deck file overrides the host flat bed on supported routes
       (coupled Rigid6 decks require a flat ``WtrDpth``)
     - ``site.bty bathymetryFile``
   * - ``kBot`` / ``kb``
     - real, Pa/m (N/m³)
     - ``1.0e5``
     - finite, > 0 when a seabed is present. Mapped to each line node and rod contact station by
       local diameter × tributary length; a Rigid6 reference point uses a fixed 1 m² area
     - both
     - ``1.0e5 kBot``
   * - ``cBot`` / ``cb``
     - real, Pa·s/m (N·s/m³)
     - ``1.0e4``
     - finite, ≥ 0 when a seabed is present; scaled like ``kBot`` and acts only on contacted
       nodes moving downward
     - both
     - ``1.0e4 cBot``
   * - ``frictionMu`` / ``mu`` / ``frictionCoefficient``
     - real, –
     - ``0`` (off)
     - finite, ≥ 0, or ``none``. A positive value needs ``WtrDpth`` or ``bathymetryFile``, and
       ``dtM``/``TMax`` in the standalone driver. Line nodes carry stick-slip friction springs; in
       a current the static solve holds the line with friction springs from its still-water laid
       shape, and the march continues them. The isotropic shorthand of
       ``frictionMuAxial``/``frictionMuLateral``. Not supported (the deck is rejected) with a
       ``FAILURE`` section, on a deck with ``Connect``/``Free`` points or ``Point3`` buoys, or
       on a deck whose lines attach to ``Rigid6`` bodies or rods unless it runs on the default
       standalone monolithic multibody march: ``staggered bodyScheme``, a ``motionFile``,
       ``Coupled``/``Vessel`` rods and every coupled OpenFAST deck with bodies or rods reject it
     - both, within the limits in this row
     - ``0.5 frictionMu``
   * - ``frictionMuAxial`` / ``frictionMuLateral`` (aliases in the rules)
     - real, –
     - ``frictionMu``
     - finite, ≥ 0. Anisotropic seabed friction of line nodes (the axial/normal pair of
       OrcaFlex, by Orcina): the coefficient along the line and across it. A row omitted takes
       ``frictionMu``; both must then be positive, or both zero. An equal pair is the isotropic
       law bit-for-bit. Rods and bodies, which have no line axis, use the lateral coefficient;
       the finite-EI two-moving-end route rejects an unequal pair. The route limits of
       ``frictionMu`` apply. Aliases: ``frictionMu_axial``, ``mu_axial``, ``muaxial`` for the
       axial coefficient; ``frictionMu_lateral``, ``mu_lateral``, ``mulateral``,
       ``frictionMuNormal``, ``frictionMu_normal`` for the lateral one. See :doc:`theory`
     - both, within the limits of ``frictionMu``
     - ``0.2 frictionMuAxial``

Clock and initial condition
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 16 12 11 31 18 12

   * - Keyword and aliases
     - Type, unit
     - Default
     - Valid values and rules
     - Route
     - Example
   * - ``bodyIC``
     - text
     - ``static``
     - ``static``: free ``Rigid6`` bodies and free rods, jointly with ``Free``/``Connect``
       points, start at their static equilibrium (weight, buoyancy, hydrostatic restoring,
       steady-current drag, seabed contact, and the attached lines re-solved). A body with no
       reachable or restrained equilibrium stops with an error naming it. ``deck``: they start
       at their deck pose, for example to release a body in a free-decay test. Coupled runs
       solve it with the host's coupled points at their initial (``PtfmInit``) pose
     - both
     - ``deck bodyIC``
   * - ``bodyWetting``
     - text
     - ``sphere``
     - ``sphere``: a ``Rigid6`` body's buoyancy, drag, fluid inertia, and added mass scale with
       the submerged fraction of an equivalent sphere at its centre of buoyancy (a dry body
       carries none). ``moordyn``: fully wet at any elevation, as in MoorDyn. Bodies with
       ``C33``/``C44``/``C55`` are always fully wet
     - both
     - ``moordyn bodyWetting``
   * - ``rodHydro``
     - text
     - ``exact``
     - ``exact``: a surface-piercing rod carries the hydrostatic moment of its displaced volume
       (waterplane second moment :math:`\pi d^4/64`). ``moordyn``: adds MoorDyn's
       :math:`\rho g (\pi d^4/32)\sin\varphi\cos\varphi` while End A is below the surface and End B
       above, reproducing MoorDyn-C's pitch stiffness
     - both
     - ``moordyn rodHydro``
   * - ``bodyHydro``
     - text
     - ``morison``
     - ``morison``: ``Rigid6`` bodies carry the fluid inertia :math:`\rho V (1 + C_a)\dot{u}`.
       ``moordyn``: they do not, as in MoorDyn
     - both
     - ``moordyn bodyHydro``
   * - ``bodyScheme``
     - text
     - ``monolithic``
     - ``monolithic``: free ``Rigid6`` bodies, ``Point3`` buoys, free and pinned rods and
       ``Free``/``Connect`` points step in one implicit generalised-alpha step with their
       ``EI = 0`` lines (a Newton iteration on the object accelerations around the line steps);
       no sub-stepping is needed for stability. A deck without lines keeps its explicit
       central-difference march unless it names ``monolithic bodyScheme``. ``staggered``: the
       staggered predictor-corrector, sub-stepped for the stiffest body-mooring mode. Decks with
       a ``motionFile`` or a ``FAILURE`` section, and coupled runs, use ``staggered`` whatever
       the value. See :doc:`theory`
     - standalone
     - ``staggered bodyScheme``
   * - ``bodySubstep``
     - text
     - ``accuracy``
     - ``accuracy``: the monolithic step divides each coupling step so that the estimated
       stiffest body-mooring frequency satisfies :math:`\omega\,\Delta t \le 0.28`, the
       accuracy of the staggered scheme's own sub-steps. ``none``: it takes ``dtM`` as given
       (stable at any ``dtM``; the body-mooring modes are then resolved only as far as ``dtM``
       allows). Ignored by ``staggered bodyScheme``, which always sub-steps
     - standalone
     - ``none bodySubstep``
   * - ``dtM`` / ``dt``
     - real, s
     - absent (static-only)
     - finite, > 0. In the standalone driver ``dtM`` and ``TMax`` are given together or not at
       all
     - standalone step; OpenFAST targets ``0.1`` s when absent, rounds it to at least one whole
       glue ``DT``, and reports the value used
     - ``0.05 dtM``
   * - ``TMax``
     - real, s
     - absent (static-only)
     - finite, ≥ 0 and an integer multiple of ``dtM``. ``TMax = 0`` runs the dynamic
       initialisation and writes the ``t = 0`` row
     - standalone only; OpenFAST uses the host ``TMax``. Do not copy the host duration into a
       coupled deck
     - ``600.0 TMax``
   * - ``RangeStart`` / ``range_start``
     - real, s
     - ``0``
     - finite, ≥ 0 and ≤ ``TMax``. First output time accumulated into the range graphs
       (``<out_root>.Line<L>.range.out``, ``LINES`` ``Outputs`` flag ``r``), so a start-up
       transient is left out. Needs ``TMax`` and at least one line with the ``r`` flag
     - standalone only
     - ``100.0 RangeStart``

Nonlinear solver
----------------

Static-solver tolerances are built in and are not deck options; keys such as ``staticRelTol`` are
rejected as unknown. The dynamic Newton controls are set with ``dynamic_solver``.

.. list-table::
   :header-rows: 1
   :widths: 16 12 11 31 18 12

   * - Keyword and aliases
     - Type, unit
     - Default
     - Valid values and rules
     - Route
     - Example
   * - ``rhoInf`` / ``rho_inf``
     - real, –
     - ``0.4``
     - finite, in ``[0, 1]``; generalised-alpha high-frequency spectral radius
     - both
     - ``0.4 rhoInf``
   * - ``maxStrain`` / ``max_strain``
     - real, –
     - ``0.5``
     - finite, ≥ 0; ``0`` disables the strain bound. Plausibility guard checked after every
       committed dynamic step of the ``EI = 0`` lines: a non-finite nodal state, or any element
       stretched beyond this engineering strain (``0.5`` = 50 %, several times the breaking
       strain of steel, polyester, and nylon lines), stops the run with exit code 2 and a message
       naming the line, element, strain, estimated elastic tension, and time. It catches a
       numerically unstable march whose Newton residual still converges; a physical run never
       reaches the default. Finite-EI cables are guarded by ``tensile_safety`` instead
     - both (``EI = 0`` lines)
     - ``0.5 maxStrain``
   * - ``dynamic_solver rel abs max_iter backtracks [rhoInf]``
     - keyword-first record: two reals, two integers, optional real
     - ``1e-8 1e-14 30 12``
     - ``rel`` and ``abs`` finite and > 0; ``max_iter`` an integer ≥ 1; ``backtracks`` an integer
       ≥ 0; the optional sixth token sets ``rhoInf`` (same range). ``EI = 0`` dynamic solves
       use all four controls. A standalone finite-EI (Hermite) deck uses only ``rel`` and
       ``max_iter`` for each step (``abs`` and ``backtracks`` are ignored on this route) and
       converges its static equilibrium to at least ``min(1e-6, rel)``. Finite-EI cables in
       OpenFAST and in mixed ``EI = 0`` + finite-EI standalone decks do not read ``rel``,
       ``abs``, ``max_iter`` or ``backtracks``: they use a fixed step tolerance (relative
       ``5e-3``, 200 Newton iterations) and a ``1e-6`` static tolerance; only the optional
       ``rhoInf`` token applies to them. Use controls qualified for
       the mesh and load case; the maintained 952-element Gulf of Maine cable uses
       ``1e-4 1e-14 100 12``
     - both (see rules)
     - ``dynamic_solver 1e-8 1e-14 30 12``
   * - ``modified_newton`` / ``modifiednewton``
     - logical
     - ``False``
     - ``True`` enables guarded within-step tangent reuse on every dynamic solver the deck builds.
       It does not relax convergence tolerances and should be qualified against full Newton
     - both
     - ``False modified_newton``
   * - ``recovery_max_substeps`` / ``recoverymaxsubsteps``
     - integer, –
     - ``1024``
     - integer from 4 to 65536. Maximum internal subdivisions after a stalled tension-only step
       or a failed or under-resolved finite-bending interval. The nominal host time grid and
       output times are unchanged; prescribed position, velocity, and acceleration follow one
       C2 quintic trajectory inside the interval. The standalone completion summary reports the
       recovered intervals and the largest subdivision count used
     - both
     - ``1024 recovery_max_substeps``

Finite-EI cable controls
------------------------

.. list-table::
   :header-rows: 1
   :widths: 16 12 11 31 18 12

   * - Keyword and aliases
     - Type, unit
     - Default
     - Valid values and rules
     - Route
     - Example
   * - ``adaptive_mesh`` / ``adaptivemesh``
     - logical
     - ``False``
     - Permits finite-EI refinement when the static curvature diagnosis or mesh-scale contact
       chatter is fragile. Multiple adequately sampled grounded runs are allowed and are not
       treated as chatter merely because the contact topology has several islands. ``False``
       keeps the declared ``NumSegs`` authoritative but does not disable the folded-element
       safety check, nor the automatic refinement of a deck mesh that cannot resolve the
       curvature of the equilibrium (reported). When ``True``, reaching the refinement cap while
       the mesh remains fragile is reported
     - both (finite-EI lines)
     - ``True adaptive_mesh``
   * - ``alpha_force_blend`` / ``alphaforceblend``
     - logical
     - ``True``
     - Generalised-α blend of the finite-EI dynamics. ``True`` blends the forces of the
       two step ends; ``False`` evaluates the forces at the blended configuration
       :math:`q_{\alpha_f}` (Chung and Hulbert), whose blended tangent vectors are shortened
       while the line rotates and raise the mean tension at large ``dtM``: on the Lozon Gulf of
       Mexico 80 m cable in 3 m surge the mean hang-off tension is 12511 / 10097 / 9507 / 9360 N at
       ``dtM`` 0.1 / 0.05 / 0.025 / 0.0125 s with ``False`` and 9424 / 9322 / 9314 / 9311 N
       with ``True``. A run prints a note when a nodal tangent turns by more than 0.25 deg
       (``True``) or 0.1 deg (``False``) in one step. A dynamic run with a torsional line
       (``END CONNECTIONS`` ``TorsStiffness`` at both ends) needs ``True``
     - standalone and coupled (finite-EI lines)
     - ``False alpha_force_blend``
   * - ``cable_statics`` / ``cablestatics``
     - mode
     - ``continuation``
     - ``continuation`` solves a finite-EI cable from its exact ``EI = 0`` catenary by
       continuation in ``EI`` and falls back to the mesh-sequenced route; ``sequenced`` runs
       the mesh-sequenced route (coarse section-preserving hierarchy, ``EI`` and buoyancy
       continuation, prolongation) first and the continuation route as its fallback. Both
       results are audited; on the Lozon cables the two routes agree to the solver
       tolerance (see :doc:`theory`)
     - both (finite-EI lines)
     - ``sequenced cable_statics``
   * - ``axial_quadrature_order`` / ``axialquadratureorder`` / ``axial_quadrature``
     - integer, –
     - ``4``
     - integer from 1 to 6; Gauss order for the finite-EI axial energy. A value different from
       the bending order gives selective integration
     - both (finite-EI lines)
     - ``4 axial_quadrature_order``
   * - ``bending_quadrature_order`` / ``bendingquadratureorder`` / ``bending_quadrature``
     - integer, –
     - ``4``
     - integer from 1 to 6; Gauss order for the finite-EI curvature energy
     - both (finite-EI lines)
     - ``4 bending_quadrature_order``
   * - ``tensile_safety`` / ``tensilesafety``
     - mode
     - ``False``
     - ``True``/``T``/``yes``/``y``/``on``/``1``/``error`` reject axial compression outside the
       tolerance, judged on the element-mean axial force averaged over three elements (see
       :doc:`theory`); ``warn``/``warning``/``monitor`` commit the state, count accepted
       integration-step events, and report the worst force, threshold, element and time after a
       standalone run (the printed ``xi`` is always ``0.5``, the centre of the three-element
       window, because the audit judges element means); ``False``/``F``/``no``/``n``/``off``/``0``
       disable the audit (case-insensitive). A response diagnostic, not a tension-only
       constitutive law
     - both (finite-EI lines)
     - ``warn tensile_safety``
   * - ``tensile_strain_tolerance`` / ``tensilestraintolerance``
     - real, –
     - ``2e-6``
     - finite, ≥ 0; strain band used by the tensile audit. Changes require an engineering basis
       and should be recorded with the run settings
     - both (finite-EI lines)
     - ``2e-6 tensile_strain_tolerance``
   * - ``cable_load_feedback`` / ``cableloadfeedback``
     - logical
     - ``True``
     - ``False`` keeps the cable march, prescribed platform kinematics, host SeaState fields, and
       response channels active, but returns zero cable force and moment to the host: a
       controlled one-way comparison, not a physical operating configuration
     - OpenFAST (coupled finite-EI cables)
     - ``True cable_load_feedback``

Ambient fluid and prescribed motion
-----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 16 12 11 31 18 12

   * - Keyword and aliases
     - Type, unit
     - Default
     - Valid values and rules
     - Route
     - Example
   * - ``current``
     - positional record, m and m/s
     - ``none``
     - ``none current``; ``uniform vx vy vz current``; or the two-level
       ``profile z1 vx1 vy1 vz1 z2 vx2 vy2 vz2 current``. Values finite; the two profile levels
       may be given in either order, are sorted by ``z``, and must differ. The static initial
       condition includes the steady drag of the current on the line at rest. Requires
       ``dtM``/``TMax`` in the standalone driver; a coupled deck takes its clock from the host.
       Use a WaterKin ``CurrentMod 1`` table for more levels; a deck with both is rejected
     - standalone; in OpenFAST, a single-turbine pure ``EI = 0`` deck without Rigid6 bodies or
       rods keeps it as a steady current when SeaState carries no waves or current. It is
       rejected when SeaState carries waves or current (double counting), in FAST.Farm, and on a
       deck with finite-EI cables, Rigid6 bodies or rods. On a pure ``EI = 0`` deck, a WaterKin
       ``CurrentMod 1`` table is the deck-side current that combines with a SeaState field
     - ``uniform 0.5 0 0 current``
   * - ``waves`` / ``wave``
     - positional record, m, s, deg
     - ``none``
     - ``none waves``; ``airy H T direction waves``; ``stream H T direction waves`` (alias
       ``dean``: the regular nonlinear stream-function wave, see ``StreamOrder``);
       ``jonswap Hs Tp gamma direction waves``; ``pm Hs Tp direction waves`` (aliases ``issc``,
       ``bretschneider``); ``torsethaugen Hs Tp direction waves``; or
       ``ochihubble Hs1 Tp1 lambda1 Hs2 Tp2 lambda2 direction waves`` (see :doc:`theory`).
       Height and period finite and > 0; ``gamma`` finite and ≥ 1; direction finite, in
       degrees. ``WaveSpreading`` spreads a spectral row. Requires ``dtM``/``TMax`` and ``WtrDpth``
       (a ``bathymetryFile`` does not satisfy it). Not combinable with WaterKin
       ``WaveKinMod 1``. A ``jonswap`` sea is synthesised from 200 long-crested components
       over ``[0.2, 5]`` times the peak frequency: one component per equal-width frequency
       bin, placed at a random frequency inside its bin with a random phase (see
       ``WaveSeed``), and scaled so the discrete spectrum gives exactly ``Hs``. The
       frequencies are not equally spaced, so the record does not repeat
     - standalone only; always rejected in a coupled OpenFAST or FAST.Farm deck, because the
       host SeaState supplies the waves
     - ``airy 2.0 8.0 0.0 waves``
   * - ``StreamOrder`` / ``stream_order``
     - integer, –
     - ``0`` (20 terms)
     - ``0`` or an integer from 2 to 60: the number of Fourier terms of a ``stream`` wave. A
       wave higher than 0.8 d or steeper than H/L = 0.142 is rejected, and a solution that
       does not converge fails with a named error
     - standalone ``stream`` waves
     - ``30 StreamOrder``
   * - ``nModes`` / ``n_modes`` / ``modes``
     - integer, –
     - ``0`` (off)
     - integer from 0 to 1000. Modal analysis of every line about its static equilibrium: the
       ``nModes`` lowest natural frequencies and mode shapes of ``K φ = ω² M φ`` (static
       tangent stiffness, structural plus added mass, ends held), written to
       ``<root>.modes.out`` before any dynamic march. Decks of lines between Fixed and
       Coupled/Vessel points on a flat ``WtrDpth`` seabed: all ``EI = 0`` lines (static or
       dynamic deck), or all finite-EI lines on the cubic-Hermite route (dynamic deck, every
       End B Fixed). Other decks are rejected; a deck that mixes ``EI = 0`` and finite-EI
       lines stops with ``OPTION nModes is not supported for mixed EI=0/finite-EI decks;
       remove it or set it to 0``; a deck with a torsional line stops with ``modal analysis
       (OPTION nModes) of a line with torsion is not yet supported``. The banded solver has no
       line-length limit: a
       1024-element finite-EI cable takes a few seconds. See :doc:`theory`
     - standalone
     - ``10 nModes``
   * - ``WaveSeed`` / ``wave_seed``
     - integer, –
     - ``1``
     - integer from 1 to 2147483646. Seed of the ``jonswap`` component frequencies and phases:
       the same seed always gives the same sea on every platform, and a different seed gives an
       independent realisation. The stream is the Park–Miller MINSTD generator (multiplier
       48271, modulus 2\ :sup:`31` − 1) started from the seed scrambled by three rounds of a
       31-bit xorshift and one MINSTD step; each component draws its frequency offset, then its
       phase. Run several seeds for extreme-value statistics
     - standalone spectral waves and wave trains (train *i* uses seed + 7919 (*i* − 1))
     - ``7 WaveSeed``
   * - ``wavetrain`` / ``wave_train``
     - positional record, m, s, deg, –
     - none
     - one wave train of a multi-train sea; the rows add up (at most 16):
       ``airy H T direction wavetrain``, ``jonswap Hs Tp gamma direction s wavetrain``,
       ``pm Hs Tp direction s wavetrain``, ``torsethaugen Hs Tp direction s wavetrain``,
       ``ochihubble Hs1 Tp1 lambda1 Hs2 Tp2 lambda2 direction s wavetrain``. ``s`` is the
       cos-2s spreading exponent in [0, 1000]; ``0`` gives a long-crested train. Each train
       has its own heading. Not combinable with a ``waves`` row or ``WaveSpreading``. Same
       requirements as ``waves``
     - standalone; in OpenFAST rejected like ``waves``
     - ``jonswap 4 9 3.3 0 4 wavetrain``
   * - ``WaveSpreading`` / ``wave_spreading``
     - real, –
     - ``0``
     - cos-2s exponent in [0, 1000] of the spectral ``waves`` row: the energy is spread over
       ``direction ± 90°`` with ``D(θ) = K(s) cos^2s(θ − direction)``. ``0`` keeps the sea
       long-crested. Rejected with a regular (``airy``) wave, without waves, or with
       ``wavetrain`` rows
     - standalone spectral ``waves``
     - ``4 WaveSpreading``
   * - ``WaveDirections`` / ``wave_directions``
     - integer, –
     - ``9``
     - integer ≥ 1; equal-angle direction bins of a spread train (each with its own
       frequency set). Long-crested trains use one direction
     - standalone spread seas
     - ``15 WaveDirections``
   * - ``WaveComponents`` / ``wave_components``
     - integer, –
     - ``200``
     - integer ≥ 2; frequency components per train and direction. A ``jonswap`` waves row
       with the default 200 and no spreading keeps the original JONSWAP synthesis; any other
       value moves it to the spectral-sea synthesis. At most 100000 components in all
     - standalone spectral waves
     - ``400 WaveComponents``
   * - ``rampTime`` / ``ramp_time`` / ``tRamp``
     - real, s
     - ``0``
     - finite, ≥ 0; ``0`` disables the ramp. Scales every wave amplitude by the half-cosine
       ``r(t) = (1 − cos(πt/rampTime))/2`` from still water at ``t = 0`` to full height at
       ``t = rampTime``; the fluid acceleration includes the ramp rate, so the ramped field is
       kinematically consistent. The current is not ramped: the static initial condition already
       carries it, so ramping it would start the line out of equilibrium. Recommended for wave
       runs, which otherwise switch the full sea on at ``t = 0`` against a still-water static
       shape; one to two peak periods is typical
     - standalone ``waves``, ``wavetrain``, and WaterKin ``WaveKinMod 1`` waves
     - ``20.0 rampTime``
   * - ``motionFile``
     - path, or ``0``/``none``
     - absent
     - an active path requires ``dtM`` and ``TMax``; ``0`` or case-insensitive ``none`` disables
       it. On a deck with a line restrained in torsion at both ends, an optional twelfth column
       rolls that line's End A frame (degrees, 0 at ``t = 0``). See the file grammar in
       :doc:`file_formats`
     - standalone ``Coupled``/``Vessel`` points, prescribed rods, and Rigid6 bodies; rejected
       on ``Connect``/``Free`` point-system, FAILURE, and mixed decks; forbidden in OpenFAST,
       where the host owns coupled motion
     - ``motion.dat motionFile``
   * - ``vesselMotion``
     - path, or ``0``/``none``
     - absent
     - a 6-DOF vessel record (19 values per row with roll/pitch/yaw, or 20 with a unit
       quaternion) on the ``dtM`` grid; every ``Coupled``/``Vessel`` point moves rigidly with the
       vessel from its deck position, and a finite or ``Rigid`` END CONNECTION at such a point
       turns with the vessel. An alternative to ``motionFile`` and ``vesselRAO`` (only one may
       be active); requires ``dtM`` and ``TMax``. Grammar in :doc:`driver_format`
     - standalone line decks accepted with a ``motionFile``; rejected on rod, Rigid6 and
       ``TURBINES`` decks and in OpenFAST
     - ``vessel.dat vesselMotion``
   * - ``vesselRAO``
     - path, or ``0``/``none``
     - absent
     - a displacement RAO table (amplitude and phase lag per DOF, per period and relative
       heading). The vessel motion is the RAO response to the deck linear waves (``airy``,
       spectral rows with ``WaveSpreading``, ``wavetrain`` rows, or a WaterKin ``WaveKinMod 1``
       file), each component at its own heading and phase, with ``rampTime``; a ``stream`` wave
       is rejected. Otherwise as ``vesselMotion``
     - as ``vesselMotion``, with deck waves
     - ``rao.dat vesselRAO``
   * - ``vesselRef``
     - ``x|y|z``, m
     - ``0|0|0``
     - three finite values: the vessel reference point at its reference pose, which is the
       rotation centre of ``vesselMotion``/``vesselRAO`` and the RAO origin (the point at which
       the wave phase is referenced)
     - ``vesselMotion``, ``vesselRAO``
     - ``90.0|0.0|-14.0 vesselRef``
   * - ``WaterKin`` / ``WaveKin``
     - ``0``/``none``, ``3``, ``7``, path, or ``SEASTATE``
     - ``0`` (no WaterKin policy)
     - ``0`` or ``none`` means still water from this record; ``3`` and ``7`` are the MoorDyn-C
       modes that read ``wave_elevation.txt`` and ``wave_frequencies.txt`` in the deck folder; a
       value containing a letter (other than an exponent ``e``) is a MoorDyn-F WaterKin
       filename; ``SEASTATE`` selects the host field. The other MoorDyn-C modes fail closed by
       name. Mode rules are in :ref:`moordyn-c-wavekin` and :ref:`waterkin-file-modes`
     - ``3``/``7`` and file ``WaveKinMod 1``: standalone decks only; ``SEASTATE``,
       ``WaveKinMod 2``, and ``CurrentMod 2``: OpenFAST only
     - ``waterkin.dat WaterKin``
   * - ``dtWave``
     - positive number [s]
     - ``0.25``
     - the step at which ``WaveKin 3`` resamples ``wave_elevation.txt`` (MoorDyn-C)
     - standalone
     - ``0.1 dtWave``
   * - ``Currents``
     - ``0`` or ``1``
     - ``0``
     - ``1`` reads the MoorDyn-C steady profile ``current_profile.txt`` in the deck folder; the
       other MoorDyn-C modes fail closed by name. ``current`` and a WaterKin table are the other
       sources; two sources in one deck are rejected
     - standalone
     - ``1 Currents``

MoorDyn compatibility-only keywords
-----------------------------------

These names exist so a migrated MoorDyn deck fails less often at the syntax boundary. They do
**not** tune CableDyn physics and should not be added to a new CableDyn model.

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Keyword
     - Value
     - CableDyn behaviour
   * - ``tScheme``
     - any token
     - Accepted and ignored; CableDyn always uses implicit generalised-alpha.
   * - ``dtIC``, ``TmaxIC``, ``CdScaleIC``, ``threshIC``
     - plain number
     - Parsed strictly and ignored; CableDyn uses a Newton static initial condition rather than
       MoorDyn dynamic relaxation.
   * - ``WriteLog``
     - plain number
     - Parsed strictly and ignored; reporting uses the normal CableDyn/OpenFAST diagnostics.
   * - ``dtOut``
     - plain number
     - Parsed strictly and ignored. Standalone output follows ``dtM``; coupled ``.CD.out``
       follows committed CableDyn solves; OpenFAST's main table follows host ``DT_Out``.
   * - ``mu_kT``, ``mu_kA``
     - any
     - Rejected by name: the MoorDyn-F seabed friction coefficients map to ``frictionMu`` with
       ``frictionMuLateral`` (for ``mu_kT``) and ``frictionMuAxial`` (for ``mu_kA``).
   * - ``mc``, ``cv``, ``FricDamp``, ``StatDynFricScale``
     - any
     - Rejected by name: CableDyn applies one regularised Coulomb friction set by
       ``frictionMu``, ``frictionMuLateral`` and ``frictionMuAxial``, with no static-to-kinetic
       ratio or friction damping. Remove the row.

Other MoorDyn constructs that fail by name rather than as unknown input: ``CoupledPinned`` and
``VesselPinned`` rods, and ``CoupledPinned`` bodies. Line-node output channels start at node 1
(``N0`` is rejected), and a number with a sign inside its mantissa (such as ``1+2``) is rejected
instead of being read up to the sign.

Water-kinematics sources
------------------------

A WaterKin file keeps its independent ``WaveKinMod`` and ``CurrentMod`` selectors. Host modes
require the OpenFAST aggregate and a compatible SeaState field; the standalone driver rejects them
with ``WaterKin WaveKinMod 2/SEASTATE is coupled-only``. The detailed supported matrix is in
:ref:`waterkin-file-modes` in :doc:`driver_format`.

With no WaterKin policy, OpenFAST uses its complete SeaState field when one exists. Explicit file
selectors are authoritative: for example, file ``WaveKinMod=0`` and ``CurrentMod=1`` means no host
waves and only the file current profile. A coupled deck cannot carry inline ``waves`` or
``wavetrain`` rows: the host SeaState supplies the waves, so such a deck is rejected at
initialisation. An inline ``current`` row is kept as a steady current only on a single-turbine,
pure ``EI = 0`` deck without Rigid6 bodies or rods whose SeaState carries no waves or current; it
is rejected when SeaState carries waves or current, in FAST.Farm, and on a deck with finite-EI
cables, Rigid6 bodies or rods.

Standalone versus OpenFAST checklist
------------------------------------

For a standalone static calculation, omit both ``dtM`` and ``TMax`` or set a deliberate pair with
``TMax=0`` when exercising the dynamic initialisation route. For standalone dynamics, provide
both and ensure ``TMax/dtM`` is integral. ``current``, ``waves``, ``motionFile``, a positive
``frictionMu``, rods, bodies, ``Connect``/``Free`` points, and finite-EI sections (outside a mixed
deck) all need the ``dtM``/``TMax`` pair. File paths are resolved relative to the deck and may use
ordinary nested paths longer than the 64-character deck-identifier limit.

For ``CompMooring = 5``, OpenFAST owns gravity, density, flat depth, platform motion, run duration,
and its output clock. The CableDyn deck owns line properties, mesh, contact coefficients,
``rhoInf``, Newton policy, optional ``dtM``, and CableDyn output channels. The console prints the
actual ``dtM`` and its integer ratio to glue ``DT``; archive that line with production results.
