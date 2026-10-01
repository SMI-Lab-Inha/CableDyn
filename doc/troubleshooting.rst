.. SPDX-License-Identifier: Apache-2.0

Troubleshooting
===============

CableDyn is designed to **fail closed with a clear message** rather than return a quiet wrong
answer. This page maps the messages you will actually see — a rejected deck, a solve that will
not converge, a combination outside the coupled boundary of OpenFAST (maintained by NLR, the
National Laboratory of the Rockies, formerly NREL) — to the reason and the supported path
forward. Questions about results and settings are answered in :doc:`faq`.

Exit codes at a glance
----------------------

.. list-table::
   :header-rows: 1
   :widths: 12 40 48

   * - Code
     - Meaning
     - First thing to check
   * - ``0``
     - completed; the completion line
       ``CableDyn_driver: converged run written to <out_root>.out`` is on stdout
     - —
   * - ``1``
     - the run was refused or its input could not be used: unknown option, wrong argument
       count, empty or over-long argument, output root equal to the deck, unwritable output
       directory, unreadable or unparseable deck or auxiliary file, out-of-range value,
       unsupported feature, or an output file that cannot be opened during the run
     - the stderr message names the file, line, keyword, or feature — see below
   * - ``2``
     - an ``EI = 0`` static line or a dynamic step did not converge, a finite-EI static
       initialisation failed or was rejected, or a result became non-finite. A message
       starting ``plausibility guard:`` means a converged dynamic step produced an element
       strain beyond ``maxStrain`` (default 50 %) — a numerically unstable march, not a
       physical load. A Windows GNU source build also stops with code 2 when it cannot load
       its LAPACK library (``openblas.dll``)
     - geometry, mesh, time step, and continuation — see :ref:`Non-convergence <non-convergence>`;
       for the plausibility guard, reduce ``dtM`` first; for a missing library, see the last
       section of this page

The complete stream and exit-code contract is in :doc:`standalone_driver`. Error messages are on
**stderr**. If you redirected it away (``2>/dev/null``), rerun without the redirect.

Reading an error message
------------------------

A row-level parse error names the file and line and quotes the offending row::

   CableDyn_DeckDriver: deck line 42: column 3 value "500/2" contains "/", ",", ";" or a repeat count "n*"; ... [row: 1 chain 500/2 ...]

Auxiliary files use their own label in place of ``deck`` (for example ``motionFile line 7`` or
``bathymetry file line 3``). Validation errors that are not tied to one row name the object
instead: a duplicate or undefined id is quoted in the message (``duplicate POINT id 4``,
``LINE 2 references undefined POINT id 9``).

.. _deck-rejected:

A deck is rejected (exit 1)
---------------------------

Fail-closed at parse is a *feature*: it means the deck named something CableDyn will not
silently reinterpret. The tables below are grouped by message class.

Command line and files
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Message
     - Cause and fix
   * - ``unknown option "--x"``
     - only ``-h``/``--help`` and ``-v``/``--version`` (as the first argument) are options; any
       other argument starting with ``-`` is refused before the deck is read.
   * - usage printed on stderr, no other message
     - the driver needs exactly two arguments: ``<deck.dat> <out_root>``.
   * - ``argument N is longer than 4096 characters`` / ``command-line argument N is empty``
     - shorten the path (run from the case directory and use relative paths) or supply the
       missing argument.
   * - ``would overwrite the input deck``
     - the output root is the deck's own stem, so ``<out_root>.out`` or ``<out_root>.static.out``
       would be the deck. Choose a different root; do not append ``.out`` to it.
   * - ``cannot write output files at "<root>"``
     - the output directory does not exist or is not writable. Create it first; the driver does
       not create directories.
   * - ``the name has characters the Windows ANSI code page cannot represent``
     - only a driver built from source with the GNU toolchain gives this: the name is outside
       the system code page and the volume has no 8.3 short name for it. Use the release
       ``CableDyn_driver.exe``, which accepts names in any script, or rename the folder.
   * - a file in a folder with accented or non-Latin characters is not found, with an older
       or self-built executable
     - the release executables run with UTF-8 as their Windows code page and open such folders
       (Windows 10 version 1903 or later). An executable built without
       ``app/utf8_code_page.manifest`` uses the system ANSI code page; run it from a folder whose
       name that code page can spell.
   * - ``cannot open deck`` / ``cannot open motionFile`` / ``cannot open bathymetry file`` /
       ``cannot open the WaterKin file`` / ``cannot open Syrope settings file`` /
       ``cannot open Syrope OWC table``
     - a referenced file is missing or unreadable. Paths inside a deck are resolved relative to
       the deck's directory (the OWC table relative to the Syrope settings file), not to the
       terminal's current directory. A ``#`` or ``!`` in a path starts a comment and truncates
       it.
   * - ``cannot open output file`` (during the run)
     - an output file became unwritable after the start-up check (for example, it is open and
       locked in another program). Close it and rerun.

Records, tokens, and sections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Message
     - Cause and fix
   * - ``deck line N: column C value "500/2" contains "/", ",", ";" or a repeat count``
     - a table value that is not one plain number. CableDyn never reads ``500/2`` as ``500``,
       ``400,0`` as ``400``, or ``2*200`` as two values. Write the number with a ``.`` decimal
       point; quote a text value (for example a type name) that must contain such a character.
       The same rule applies to numeric ``OPTIONS`` values (``1/20 dtM`` is rejected as
       ``malformed numeric OPTION value``) and to ``motionFile``, bathymetry, WaterKin, and
       Syrope table rows.
   * - ``deck line N is longer than 512 characters``
     - a record whose text outside a ``#``/``!``/``--`` comment does not fit the record buffer.
       Split a long ``OUTPUTS`` list over several rows or move commentary behind ``#``. Long
       comment lines are accepted. Auxiliary files report the same error with their own label.
   * - ``unknown deck section "<NAME>"``
     - a record containing ``---`` is always read as a section header. Use one of the header
       names listed in :doc:`driver_format`, and keep ``---`` out of data rows and descriptions.
   * - ``LINE TYPES row needs 10 cols ...`` / ``POINTS row needs 5 ... or 9 ...`` /
       ``LINES row needs 3 ..., 4 ..., or the stock 7 cols ...`` /
       ``SECTIONS row needs 4 cols``
     - a row with the wrong number of columns. A stray comment without ``#``/``!`` or an
       unquoted name with a space adds a column.
   * - ``malformed LINES row: attachment X must be an unquoted point id or rod end (R<N>A,
       R<N>B)``
     - a LINES attachment is neither a point id nor a rod end. Rod ends can be named directly
       (``R1A``, ``Rod1B``); CableDyn creates the end point. A body attachment needs a POINT of
       type ``Body<N>``; attach the line to that point id.
   * - ``OUTPUTS channel name "..." is longer than 64 characters``
     - shorten the channel; valid channel names are far shorter than the limit.
   * - ``unknown OPTION keyword "..."``
     - a typo in an ``OPTIONS`` keyword, or a keyword that does not exist (for example
       ``staticRelTol``: static tolerances are built in; dynamic Newton controls use
       ``dynamic_solver``). Check the spelling against :doc:`options`.
   * - ``unknown OPTION keyword "2.0"`` on a ``current``/``waves`` row
     - text after the terminal ``waves``/``current`` keyword that is not introduced by a spaced
       ``-`` adds tokens and hides the positional form. Write ``airy 2.0 8.0 0.0 waves -
       description`` or move the text behind ``#``.
   * - ``malformed logical OPTION value`` / ``malformed tensile_safety value``
     - use one of the spellings listed in :doc:`options` (``True``/``False``; ``warn`` for
       ``tensile_safety``).
   * - ``ICmode is not an option``
     - CableDyn always solves the static equilibrium directly; delete the ``ICmode`` row.
   * - a channel like ``Point2px_raw`` is rejected
     - an ``OUTPUTS`` channel that does not match a supported form, names an unknown id, or has
       trailing text. See the exact channel grammar in :doc:`outputs`.

Clock, motion, and environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Message
     - Cause and fix
   * - ``dynamic run requires both dtM and TMax options``
     - the standalone driver needs both or neither. Add the missing one.
   * - ``OPTIONS TMax must be an integer multiple of dtM``
     - choose ``TMax = n·dtM`` exactly (for example ``dtM = 0.05``, ``TMax = 600``).
   * - ``standalone OPTION current requires dtM and TMax`` (likewise ``waves``, ``motionFile``,
       ``frictionMu``, ``RODS``, ``Rigid6 BODY``, ``Point3 BODY``, ``Connect``/``Free`` points,
       ``finite-EI section (EI>0)``)
     - these features act only in a dynamic march. Add ``dtM``/``TMax``, or remove the feature for
       a static calculation. A mixed ``EI = 0`` + finite-EI deck may run static without them. A
       coupled (caller-driven) deck takes its clock from the host, so ``current`` and ``frictionMu``
       need no ``dtM`` there.
   * - ``OPTION waves requires WtrDpth``
     - waves need a flat depth; a ``bathymetryFile`` alone does not satisfy the wave model.
   * - ``seabed friction is not supported on decks with Connect/Free points, bodies or rods on
       this route (coupled Rigid6/ROD decks, point decks and the staggered body scheme)``
     - these routes have no seabed friction. Remove ``frictionMu``, or run a standalone body or
       rod deck on the default ``monolithic`` body scheme.
   * - ``OPTION nModes is not supported for mixed EI=0/finite-EI decks; remove it or set it to 0``
     - a deck that mixes ``EI = 0`` and finite-EI lines writes no modal files. Run the modal
       analysis on a deck whose lines are all ``EI = 0`` or all finite-EI.
   * - motionFile ``time is outside the dtM/TMax grid``
     - every row time must be ``k·dtM`` for ``0 ≤ k ≤ TMax/dtM``. Resample the history onto the
       ``dtM`` grid.
   * - ``motionFile must provide every Coupled/Vessel point at every dtM time`` (or ``every
       eligible prescribed point``)
     - the file is missing at least one (point, time) pair. Every prescribed point needs a row at
       every grid time, including ``t = 0`` and ``t = TMax``.
   * - motionFile ``rows must be: time point_id x y z vx vy vz ax ay az``
     - a row has fewer than eleven plain numbers, or a header line is not commented. Comment
       header lines with ``#``.
   * - motionFile ``... which is not Coupled/Vessel`` / ``duplicate row for point``
     - the ``point_id`` is not a prescribed point, or a (point, time) pair appears twice.
   * - ``Connect/Free dynamic-point deck does not support motionFile coupling`` /
       ``mixed-deck motionFile is not supported``
     - these routes hold their coupled endpoints. Remove ``motionFile`` or separate the model.
   * - ``bathymetry file must define one complete rectangular grid`` / ``grid is incomplete`` /
       ``duplicate x/y entries`` / ``needs at least a 2x2 grid`` / ``depths must be positive``
     - the ``x y depth`` rows must cover every combination of the distinct ``x`` and ``y``
       values exactly once (any order), with at least two of each and positive depths.
   * - ``OPTIONS WtrDpth and bathymetryFile are mutually exclusive``
     - keep one seabed definition.
   * - ``Wheeler stretching is undefined: the wave trough reaches the seabed``
     - the wave trough is at or below the seabed (``depth + eta <= 0``). Reduce the wave height or
       increase ``WtrDpth``; the model is not valid in that sea state.
   * - ``WaterKin WaveKinMod 2/SEASTATE is coupled-only``
     - ``SEASTATE``, ``WaveKinMod 2``, and ``CurrentMod 2`` take kinematics from the OpenFAST
       SeaState field. In the standalone driver use ``WaveKinMod 0``/``1`` and ``CurrentMod 0``/
       ``1``, or native ``waves``/``current`` rows.
   * - ``WaterKin WaveKinMod 1 is supported only by the standalone driver``
     - an elevation-history file cannot drive a coupled run or a mixed standalone deck. Put the
       sea state in SeaState (coupled) or use a single-family standalone deck.
   * - ``WaterKin WaveKinMod 1 requires ...`` / ``WaveKinFile ...``
     - ``WaveKinMod 1`` needs a ``WaveKinFile`` of at least four finite ``time elevation`` rows
       starting at ``t = 0`` with increasing times, a positive ``dtWave``, and a positive
       ``TMax``; it cannot be combined with a ``waves`` OPTION or ``CurrentMod 2``.
   * - ``... BOTH a current OPTION and a WaterKin current (double-counting)``
     - keep either the ``current`` row or the WaterKin ``CurrentMod 1`` table.
   * - a WaterKin wave/current combination is rejected
     - the requested host and file sources cannot be separated by the nodewise fluid sampler
       without double counting. Use independent sources that :doc:`capabilities` lists, or move
       the complete environment into OpenFAST SeaState.

Model data
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Message
     - Cause and fix
   * - ``has a non-finite property`` / ``POINT has non-finite position or load data`` /
       ``BODY has non-finite data`` / ``motionFile contains a non-finite value``
     - a ``NaN`` or ``Inf`` reached the deck, usually from a generating script. Fix the source
       value; CableDyn does not substitute defaults.
   * - ``LINE End A fairlead (NodeA) must not be below End B anchor (NodeB)``
     - End A is the fairlead/upper end. A stock anchor-first row (``Fixed`` NodeA, fairlead NodeB)
       is swapped automatically; this message means the fairlead point itself lies below the
       anchor. Check the point coordinates.
   * - ``LINE End A (NodeA) must be Coupled/Vessel/Body<N> ...`` / ``LINE End B (NodeB) must be a
       Fixed (anchor) point``
     - in a deck without dynamic points every line runs from a fairlead to an anchor. Use
       ``Connect``/``Free`` points for line-to-line junctions.
   * - ``finite-EI section (EI>0) requires a dynamic deck with dtM and TMax``
     - a single-family finite-EI deck is solved by the dynamic finite-EI workflow. Add
       ``dtM``/``TMax`` (``TMax = 0`` gives the static configuration only).
   * - ``is net-buoyant (dry mass <= displaced water mass)``
     - an ``EI = 0`` section whose dry mass ≤ displaced mass (a buoyant arch) is outside the EI = 0
       route — a lazy-wave arch needs the finite-EI path. Give the buoyant section an ``EI > 0``
       line type, or specify it by target submerged weight (the ``EQUIVALENT BUOYANCY`` block in
       :doc:`driver_format`).
   * - ``the viscoelastic dynamic stiffness Ed (EA = Es|Ed) must be finite and exceed the static
       Es``
     - ``Ed`` must be strictly greater than ``Es``; equal values make the series spring infinite.
       For ``Es|alphaMBL|vbeta`` both extra values must be positive. Viscoelastic types need
       ``EI = 0``.
   * - ``a LINE TYPES BA entry cannot have more bar-separated values than its EA entry``
     - write ``BA`` as one value for a plain EA and at most ``Bs|Bd`` for a viscoelastic EA.
   * - ``a Syrope line does not support ...`` / ``a Syrope line must be a single section`` /
       ``must be taut at init`` / ``has no static output``
     - a Syrope line is one ``EI = 0`` section, taut at initialisation, marched with
       ``dtM``/``TMax``, without ``bathymetryFile``, deck current/waves, or host fluid loads (a
       flat ``WtrDpth`` is accepted). Adjust the model accordingly.
   * - ``Syrope alpha and beta ... must be positive`` / ``the Syrope BA column needs two parts`` /
       ``the Syrope settings file needs OWC, WCType, k1, and k2 rows``
     - see the Syrope column rules in :doc:`driver_format`: ``EA = SYROPE:<file>|alpha|beta``
       with positive ``alpha``/``beta``, ``BA = BA_s|BA_d`` non-negative with a positive sum.
   * - ``SYROPE IC references undefined LINE id``
     - the ``SYROPE IC`` section must come after ``LINES``; also check the line id.
   * - ``Point3 BODY requires Mass > 0 ...`` / ``Point3 BODY currently supports exactly one
       Body<N> attachment point``
     - a ``Point3`` body needs positive mass, non-negative ``Vol``/``CdA``/``Ca``, and exactly
       one ``Body<N>`` point.
   * - ``Rigid6 BODY requires Mass/Ixx/Iyy/Izz > 0 ...`` / ``BODIES row needs 15 cols plus
       optional Ixx Iyy Izz``
     - a ``Rigid6`` row needs the three trailing inertias (18 columns), all positive.
   * - ``Rigid6 dynamic deck cannot mix ...`` / ``Rod dynamic deck cannot mix ...``
     - each standalone route owns one object family. Split the model or use the combinations
       listed in :doc:`capabilities`.
   * - ``a ROD may have at most one Rod<N>A point and one Rod<N>B point`` /
       ``Coupled/Vessel ROD requires a motionFile (or a coupled host)`` / ``ROD endpoints must
       define a positive length``
     - remove the duplicate rod-end POINT row (the end points are created automatically when a
       line names the rod end), supply ``motionFile`` rows for both ends of a prescribed rod, and
       give the rod distinct end coordinates.
   * - ``FAST.Farm Turbine<J> decks are consumed by the CompMooring=5 aggregate path only``
     - ``Turbine<J>`` points use turbine-local coordinates. Run the deck through FAST.Farm, or
       replace them with ``Coupled`` points in global coordinates for a standalone run.
   * - ``END CONNECTIONS requires a finite-EI line`` / ``non-pinned END CONNECTIONS require a
       finite-EI line with Fixed End B``
     - bending end connections apply only to ``EI > 0`` lines whose End B is an anchor.
   * - ``a clamped or elastic END CONNECTION on a body or rod needs bodyScheme monolithic`` /
       ``the END CONNECTION direction on a rod end must be parallel to the rod axis``
     - remove ``staggered bodyScheme``; on a rod end give the direction along the rod axis
       (either sense).

Torsion
~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Message
     - Cause and fix
   * - ``END CONNECTIONS row needs 6 columns ..., or 10 or 11 with the torsion columns``
     - a torsion row is ``LineID End Stiffness EzX EzY EzZ TorsStiffness NxX NxY NxZ
       [Pretwist]``; see :doc:`driver_format`.
   * - ``END CONNECTIONS torsional stiffness must be finite and non-negative, Free, or Rigid`` /
       ``torsion reference normal (NxX NxY NxZ) must be non-zero`` / ``must not be parallel to
       the direction Ez``
     - give ``Free``, ``Rigid`` or a stiffness in N·m/rad, and a reference normal that is not
       along ``Ez``.
   * - ``line <L> is torsionally restrained at both ends: its LINE TYPES row ... must give an
       explicit GJ > 0``
     - torsion has no ``EI/1.3`` default: use the 14-column ``LINE TYPES`` row with ``GJ`` (and
       positive ``GAs``, ``Irt``, ``Irn``) for every section of the line.
   * - ``torsional END CONNECTIONS ... require a finite-EI line`` / ``require a finite-EI line
       with Fixed End B``
     - torsion is solved on cubic-Hermite lines whose End B is an anchor.
   * - ``a torsional END CONNECTION on a rod end is not supported`` / ``on a body needs a Rigid6
       body`` / ``torsion is not combined with ATTACHMENTS in this build``
     - outside the torsion scope (:doc:`capabilities`). Move End A to a ``Fixed``,
       ``Coupled``/``Vessel`` point or a Rigid6 body, or represent the modules as a smeared
       buoyancy section.
   * - ``torsion is not yet supported in coupled OpenFAST runs`` / ``torsion is not yet supported
       in a deck that mixes EI = 0 and finite-EI lines without a BODY, nor in coupled OpenFAST or
       FAST.Farm runs``
     - run the torsional line standalone (a mixed deck with a body runs on the multibody route,
       which supports torsion), or set ``TorsStiffness Free`` in every ``END CONNECTIONS`` row.
   * - ``a body holding a line restrained in torsion turns by more than 90 deg in one step even
       after 6 step halvings``
     - the body turns faster than one step can follow the twist; reduce ``dtM``.
   * - ``torsion in a dynamic run needs the force-blended generalised-alpha (OPTION
       alpha_force_blend True)`` / ``modal analysis (OPTION nModes) of a line with torsion is not
       yet supported``
     - remove ``False alpha_force_blend`` or ``nModes`` from a deck with torsion.
   * - ``OUTPUT "Torq...": line <L> is not torsionally restrained at both ends ..., so it carries
       no torque``
     - a ``Torq``/``Twist`` channel needs a line with a torsional restraint at both ends.
   * - ``motionFile: point <id> has a non-zero roll column, but no line restrained in torsion at
       both ends ... has its moving (non-Fixed) end there`` / ``the roll column of point <id> must
       be 0 at t = 0`` / ``gives the roll column (12th) on some rows only`` / ``the roll column
       (12th, degrees) must be a finite number``
     - the roll column drives the frame of a torsional line's moving end, starts from 0 (put a
       constant twist in ``Pretwist``), and appears on every row of the point or on none
       (:doc:`file_formats`).

FAILURE and CONTROL sections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - Message
     - Cause and fix
   * - ``FAILURE row needs 5 cols`` / ``FAILURE row has trailing content``
     - write ``FailID Point Lines FailTime FailTen`` with the line list comma-separated and no
       spaces (``1,2``).
   * - ``FAILURE FailIDs must be sequential from 1``
     - number the rows 1, 2, 3, … in deck order.
   * - ``FAILURE ... lists line L, which is not attached to point P``
     - every listed line must end at the failing point.
   * - ``FAILURE needs FailTime > 0 or FailTen > 0``
     - a row with both triggers at zero could never fire.
   * - ``FAILURE at a rod end (R<n>A/R<n>B) is not supported``
     - failures act on POINTS; name the point id (``P<n>`` or an integer).
   * - ``FAILURE requires a dynamic deck`` / ``FAILURE is not supported on finite-EI decks`` /
       ``... on ROD decks`` / ``FAILURE deck does not support motionFile coupling``
       / ``... seabed friction``
     - FAILURE runs on the ``EI = 0`` point-system and Rigid6 dynamic routes without rods,
       ``motionFile``, or ``frictionMu``.
   * - ``FAILURE on a Rigid6 BODY deck does not support Connect/Free or Point3 body points``
     - a FAILURE deck with Rigid6 bodies runs on the staggered body scheme, which takes the
       bodies and their lines only; model a clump or buoy on the line as part of a body.
   * - ``this entry point does not support the CONTROL section``
     - the standalone driver does not consume CONTROL. Line control is an OpenFAST
       ``CompMooring = 5`` feature for ``EI = 0`` decks; remove the section for a standalone run.
   * - ``CONTROL row has trailing content`` / ``line L is assigned to more than one CONTROL
       channel``
     - write ``ChannelID Lines`` with a comma list without spaces, and assign each line once.

.. _non-convergence:

A solve does not converge
-------------------------

An ``EI = 0`` static line that does not converge, or a dynamic step that does not converge after
the recovery substeps, ends the run with exit code ``2``; the ``.out`` holds the converged prefix
for inspection only. A finite-EI static initialisation failure (``finite-EI cable static solve
failed`` or ``unresolved/kinked finite-EI equilibrium``) also ends the run with exit code ``2``.

The ``EI = 0`` catenary path and the finite-EI static solve are both robust from a trivial seed,
but a few geometries are genuinely hard:

- **A large grounded touchdown on the** ``EI = 0`` **path.** A long, near-horizontal grounded run can
  make the linear-Lagrange tangent singular regardless of seed. If the case is a lazy-wave or
  bending-dominated cable, use the finite-EI path (``EI > 0``), which is built for it.
- **A finite-EI lazy-wave cable at too coarse a mesh.** An under-resolved arch or touchdown
  carries a large element ``h·κ`` and can stall the dynamic Newton. The resolution metric flags
  this; the builder can auto-refine a fragile mesh when the deck opts in (``adaptive_mesh True``).
  This remains true when a large net-buoyant mesh first uses coarse-to-fine sequencing: the
  converged sequenced state is still diagnosed and, when necessary, warm-start refined. Otherwise
  **increase** ``NumSegs`` on the sag-bend and hang-off sections.
- **Several separated finite-EI contact regions.** These can be physical: a cable can bridge a
  bathymetric rise, or a buoyant section can separate two broad grounded heavy sections. CableDyn
  records the island count but does not reject that topology. Adaptive refinement is requested
  only when contact classification toggles over one interior node (an isolated contact node or
  isolated suspended gap), which is a mesh-scale signal; refine the local sections if that
  signal remains active.
- **A residual-converged folded branch.** CableDyn checks the element rotation measure after every
  finite-EI sequence attempt, every final deck-built equilibrium, and the public automatic-mesh
  route. ``max(h*kappa) > 0.7`` is rejected by name even if Newton's residual is small. This is a
  mesh/branch safety failure, not permission to relax the tolerance. Increase resolution in the
  reported high-curvature section and inspect the static profile. The check samples element
  interiors adaptively, so a cubic-element fold cannot hide between output nodes.
- **Odd, prime, or section-wise nonuniform meshes.** These are supported; do not change
  ``NumSegs`` merely to make the total divisible by two. The nested hierarchy retains physical
  section interfaces and interpolates at the true cumulative rest-length coordinate. Its
  suspended-span bending-scale preflight checks the actual longest proposed coarse element, so a
  sparse long section is not over-coarsened merely because the deck also contains many short
  elements. Flat-bed installed lines with an already supported anchor retain their grounded-tail
  branch homotopy and are still accepted only after exact-mesh residual and kink checks.
- **A very slack cable on a frictionless flat bed.** If its rest length exceeds the vertical-plus-
  horizontal L-shaped path between the endpoints, a smooth single-touchdown equilibrium may not
  exist: the surplus must fold or form another contact region. CableDyn rejects a resulting
  element-localized wad through the curvature safety check. This is not fixed by relaxing Newton
  tolerances; revise the physical support/contact model or geometry.
- **A twisted finite-EI line.** ``torsion continuation stalled at Phi = ... (target ...)``
  reports the imposed twist the static ramp reached: beyond it the line has no nearby static
  equilibrium on the path, typically because a loop is forming (hockling), which needs a
  dynamic analysis and is not resolved as self-contact. ``finite-EI torsion static solve ended on
  an unstable equilibrium`` means the descent from a buckled (unstable) twisted state did not
  find a stable one. ``the twist moved more than pi/2 from its committed value in one step``
  stops a dynamic step that would change the twist by more than a quarter turn (the unwrapping
  of the twist cannot be trusted beyond it): reduce ``dtM`` or slow the imposed roll.
  ``tangent turns by more than 120 degrees`` names an element or node pair where the twist of
  the centreline is no longer reliable: refine the mesh there.
- **A cold start at full load.** The static solve continues the load in stages; if you have
  hand-tuned a case into a bad basin, remove the tuning and let the default continuation run
  from the catenary/arch seed.

If a *dynamic* step fails after a good static IC, the usual cause is a time step too large for a
transient event (a snap) — but note the implicit integrator's whole point is large steps, so
first confirm the static IC itself is smooth (plot the ``<out_root>.static.out`` curvature).
The failure message reports the simulated interval; ``recovery_max_substeps`` raises the
internal subdivision ceiling, and ``dynamic_solver`` sets the dynamic Newton tolerance and
iteration budget (see :doc:`options`). Qualify any change against a smaller ``dtM``.

Warnings and notes that do not stop a run
-----------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 42 58

   * - ``Note: line <L> is torsionally restrained at one end only; ... no torsion is solved``
     - with the other end free to twist the line carries no torque, so the torsion columns of the
       restrained end have no effect. Restrain both ends to solve torsion.
   * - ``Note: line <L>: torsion solved, imposed twist ... deg, torque ... N m (<n> twist
       stage(s), <m> buckling descent(s))``
     - the static twist stage of a torsional line. A non-zero descent count means the straight or
       untwisted branch was unstable at the imposed twist and the solve moved to a buckled shape;
       check the shape and the ``Twist<L>`` channel.

   * - Message
     - Cause and fix
   * - ``WARNING: line type ... sinks ... m into the seabed ...; kBot is soft for this line``
     - a grounded run of this type would sink deeper than its own diameter into the penalty
       seabed. Raise ``kBot``; the default ``1.0e5`` embeds a heavy chain by about half its
       diameter, which is normal.
   * - ``WARNING: LINE ... hangs ... m below its Fixed anchor and the deck has no seabed``
     - the deck has no ``WtrDpth`` or ``bathymetryFile``, so a slack mooring hangs through the
       missing seabed. This is typical of an OpenFAST ``MooringFile``, which takes the depth from
       the host; add ``WtrDpth`` for a standalone run.
   * - ``Note: motionFile point ... starts moving at ... m/s``
     - the lines start at rest, so an end that starts at speed sends an axial shock down the line
       and can put its upper part into compression for several steps. A harmonic started at
       ``t = 0`` begins at its peak velocity. Start the motion from rest or ramp its amplitude in
       over one or two periods, as ``examples/data/lozon/generate_gomex80_heave.py`` does.
   * - ``Note: line ...: the static layout carries the grounded run ... m past its anchor, where
       it folds back (a hairpin) ...``
     - a very strong current on a frictionless seabed has pushed the line past its anchor, and
       it folds back to it (see :doc:`theory`). The layout is a valid, stable equilibrium of the
       deck as written, but not the usual one. If the real seabed holds the line, declare
       ``frictionMu``; otherwise check the line length and the current, and inspect
       ``<out_root>.static.out`` before using the layout.
   * - ``Note: line ...: the equilibrium in the current is held in the vertical plane of the line; it
       is unstable out of that plane ...``
     - in a current along the plane of the line the static state is a planar equilibrium that a
       lateral disturbance would leave (the line is held in the plane only by the symmetry), and
       no out-of-plane equilibrium was found. The state is valid in the plane; a cross-flow
       component or seabed friction (``frictionMu``) gives the line a stable layout.
   * - ``Note: line ...: a nodal tangent turned up to ... deg in one step``
     - the time-step audit of a finite-EI line: above 0.25° per step (0.1° with
       ``False alpha_force_blend``) the mean tension can err by about 1 %. Reduce ``dtM``.
   * - ``Note: line ...: elements as short as ... are below sqrt(EI/EA)``
     - the finite-EI deck mesh is finer than the axial-bending length :math:`\sqrt{EI/EA}`:
       the axial stiffness of each element swamps its bending stiffness and the static
       Newton system loses the precision to converge. An 80 m lazy wave (6.5 mm limit)
       solves at 8 mm elements and not at 4 to 6 mm, where the solve stops within seconds
       and names the mesh. Coarsen ``NumSegs`` to at least the length the note gives.
   * - ``Note: line ... mesh too coarse for its bends``
     - the element-mean axial force of a finite-EI line, averaged with its neighbours, is
       compressive, and the pointwise resultant oscillates well beyond it. Refine ``NumSegs`` in
       the reported region, or set ``tensile_safety True`` to reject the state. A pointwise dip
       at a seabed touchdown with tensile neighbouring elements is not reported.
   * - curvature near a hang-off changes with ``dtM`` or the mesh
     - end curvature is a boundary-layer quantity. Run a ``dtM`` study (halve it until the value
       changes by less than a few percent), grade the mesh toward the end, and report the value
       1–2 m from a pinned end, not at the pinned node (see :doc:`modeling`, step 5).

An OpenFAST ``CompMooring = 5`` run fails at init
-------------------------------------------------

The module covers the scoped MoorDyn-F coupled surface and **fails closed** outside it. If init
aborts with a CableDyn fatal error, check whether the deck uses one of these combinations:

.. list-table::
   :header-rows: 1
   :widths: 46 54

   * - Not supported on the coupled route
     - Use instead
   * - active ServoDyn control (``CONTROL``) on a deck with a **finite-EI** cable, or in FAST.Farm
     - ``EI = 0`` line control **is** supported in a single-turbine model; for a controlled cable,
       use stock MoorDyn (``CompMooring = 3``) in the same binary
   * - deck ``waves``/``wavetrain`` OPTIONS, and a deck ``current`` on a deck with finite-EI
       cables, Rigid6 bodies or rods, in FAST.Farm, or with a SeaState that carries waves or
       current
     - put the sea state and current in SeaState, or a current table in a WaterKin file
       (``CurrentMod 1``); a steady deck ``current`` is kept on a single-turbine, pure
       ``EI = 0`` deck in a SeaState without waves or current
   * - a finite-EI cable in the same deck as Rigid6 bodies, rods, or ``Connect``/``Free`` points
     - model the cable and the body/point system in separate decks or use the supported
       combinations in :doc:`capabilities`
   * - OpenFAST linearisation with a platform-relative cable ``END CONNECTIONS`` row
       (``formal OpenFAST linearization with a platform-relative cable end connection is not
       supported``) or with a ``CONTROL`` section
     - linearise a pinned-hang-off, uncontrolled variant of the deck; time-domain runs support both
   * - a ``Rigid`` end-connection direction reversed by exactly 180 degrees in one host update
     - drive the platform through intermediate orientations
   * - VIV (``compVIV``)
     - use a solver with a separately validated VIV model
   * - a host/file fluid combination named as inseparable
     - select host wave/current and file wave/current independently as documented, or place the
       whole field in SeaState
   * - a ``Point3`` restoring model beyond buoyancy and supported fluid loads
     - represent the restoring physics with supported line/body terms or use stock MoorDyn

Finite-EI rotational end coupling **is** supported in the time domain on the coupled route: an
``END CONNECTIONS`` direction at a coupled hang-off rotates with the OpenFAST platform
orientation, the mesh angular velocity and acceleration enter the implicit residual, and the
connection moment is returned to OpenFAST on the load mesh (see :doc:`openfast`). In the
standalone driver an end-connection direction stays fixed in global axes on a ``Fixed`` or
``Coupled`` point driven by a ``motionFile``, and turns with the vessel under ``vesselMotion``
or ``vesselRAO``.

Everything CableDyn does support end-to-end — coupled time-domain analysis, ``PtfmInit``,
``NumCrctn > 0``, SeaState kinematics, checkpoint-restart, linearisation, rigid bodies/rods, and
``Mod_SharedMooring = 5`` — is listed in :doc:`capabilities`.

The standalone binary will not launch (Windows)
-----------------------------------------------

Use the release ``CableDyn_driver.exe`` or CableDyn-enabled ``openfast.exe``. Those two
non-OpenMP x64 files are statically linked and require no adjacent compiler, BLAS, OpenMP,
or C/C++ runtime DLL. A locally built GNU/OpenMP executable may still require its selected
toolchain libraries. A GNU build loads ``openblas.dll`` from its own directory when it starts;
if the library is missing, the driver stops with exit status ``2`` and a message naming the
library, each location tried, and the Windows load error. Keep ``openblas.dll`` next to the
executable (``cmake --install`` stages it) or run from the activated environment. A missing
``DISCON.dll`` is different: the OpenFAST model explicitly
requested a ServoDyn controller, so copy that model-specific controller with the model or
disable it. See :doc:`installation` (standalone binaries).

If your issue is not here, see :doc:`faq`; the :doc:`driver_format` row for the feature is the
authoritative statement of what is supported, and :doc:`validation` records exactly which cases
are validated.
