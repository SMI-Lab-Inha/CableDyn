.. SPDX-License-Identifier: Apache-2.0

OpenFAST with ``CompMooring = 5``
=================================

This page is the reference for running CableDyn inside OpenFAST, maintained by NLR (National
Laboratory of the Rockies, formerly NREL). For a step-by-step first run of
the IEA-15MW VolturnUS-S model with the release executable, follow :doc:`tutorial_openfast`.
Every CableDyn keyword, its default, and whether OpenFAST or the deck owns it are listed in
:doc:`options`; units, frames, signs, and the tension and angle definitions are in
:doc:`conventions`.

CableDyn is statically linked into the CableDyn-enabled ``openfast.exe`` (OpenFAST v5.0.0).
There is no DLL to copy and no second process: OpenFAST calls the module in memory and exchanges
platform motion, line loads, SeaState fields, checkpoints, and linearisation data through its
normal glue code. A stock OpenFAST executable does **not** recognise ``CompMooring = 5``; use the
executable from the CableDyn release or build the integration from ``integration/openfast/``
(:doc:`coupling`). Stock MoorDyn (``CompMooring = 3``) remains available in the same binary.

Required files
--------------

A coupled case contains:

* the CableDyn-enabled ``openfast.exe``;
* the top-level ``.fst`` and every turbine sub-model it references (ElastoDyn, SeaState,
  HydroDyn and its WAMIT data, AeroDyn, InflowWind, ServoDyn, …);
* one CableDyn ``MooringFile`` deck, plus any file that deck references (WaterKin, bathymetry,
  Syrope settings); and
* a controller DLL such as ``DISCON.dll`` only if the selected ServoDyn input names one. It
  belongs to the turbine model, not to CableDyn.

Relative paths inside the ``.fst`` are resolved by OpenFAST from the ``.fst`` folder; paths inside
the CableDyn deck are resolved from the deck's folder.

Configure the ``.fst``
----------------------

.. code-block:: text

         5   CompMooring     - {0=None; 1=MAP++; 2=FEAMooring; 3=MoorDyn; 4=OrcaFlex; 5=CableDyn}
   "CableDyn_UMaine.dat"  MooringFile

``CompMooring = 3`` with a MoorDyn file selects stock MoorDyn instead; nothing else in the model
has to change for an A/B comparison.

What OpenFAST owns and what the deck owns
-----------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Owner
     - Quantities
   * - OpenFAST
     - ``TMax``, glue step ``DT``, ``TStart``, ``DT_Out``, output format; gravity, water density,
       and flat water depth (they override the deck values before statics); platform position,
       velocity, and acceleration; the SeaState wave and current field
   * - CableDyn deck
     - line types, points, lines, sections, end connections, seabed contact (``kBot``, ``cBot``,
       ``frictionMu``, optional ``bathymetryFile``), solver settings, ``rhoInf``, ``dtM``,
       WaterKin source selection, and the ``OUTPUTS`` list

``motionFile`` and deck ``waves`` and ``wavetrain`` OPTIONS are rejected in a coupled deck,
because the host supplies motion and waves. A deck ``current`` is kept as a steady current only on
a single-turbine, pure ``EI = 0`` deck without Rigid6 bodies or rods when SeaState carries no waves
or current; it is rejected when SeaState carries waves or current, in FAST.Farm, and on a deck with
finite-EI cables, Rigid6 bodies or rods. On a pure ``EI = 0`` deck without Rigid6 bodies or rods,
a depth-dependent current from a WaterKin ``CurrentMod 1`` table combines with the SeaState field
(see *Fluid source*). ``Coupled``
(or ``Vessel``) points are fairleads given in platform
axes; ``Fixed`` points are anchors in global axes.

Time step
~~~~~~~~~

``DT`` is the OpenFAST glue step; ``dtM`` is CableDyn's committed step. With no deck ``dtM``,
CableDyn solves every glue step (``dtM = DT``). A deck ``dtM`` that is a whole multiple of ``DT``
solves once per that many glue steps, which is faster. Between commits the coupled loads are
extrapolated from the last two commits (deck ``OPTIONS predictLoads``, the default) or, with
``predictLoads False``, held at the last commit, which delays them (see the note below). The
console reports the result, for
example ``CableDyn time step dtM = 2.50000E-02 s (1 x glue DT; deck dtM)``. A deck ``dtM``
larger than ``DT`` but not a whole multiple of it is never coarsened: CableDyn uses the largest
multiple of ``DT`` that does not exceed it and reports a warning with both values (``dtM = 0.025``
with ``DT = 0.01`` runs at 0.02 s), matching MoorDyn, which only ever shortens its step to fit the
coupling interval. A deck ``dtM`` below ``DT`` (typical of MoorDyn decks written for explicit
integration) runs at ``dtM = DT``, reported on the console. ``TMax`` must end
on a commit boundary. For fatigue, snap, contact, or local finite-EI response compare at least
two ``dtM`` values (for example 0.1, 0.05, and 0.025 s) and demonstrate convergence of ranges,
peaks, and damage. If CableDyn needs a smaller step than the glue step, reduce ``DT`` as well.
``DT_Out`` changes only the sampling of the main OpenFAST table.

.. note::

   With ``dtM`` larger than ``DT``, the loads CableDyn returns at a glue step a fraction
   :math:`\tau = (t - t_c)/\mathrm{dtM}` past its last commit :math:`t_c` are extrapolated from
   that commit and the one before, :math:`F = F_c + \tau\,(F_c - F_p)`. Their error is of
   second order in ``dtM`` (:math:`\tau(1+\tau)\,\mathrm{dtM}^2 \ddot F/2`). Held instead
   (``predictLoads False``), the platform receives the mooring load on average half a ``dtM``
   late, and a spring force delayed this way acts on the platform as a negative damping of about
   the mooring stiffness times ``dtM/2``.

   On the IEA-15 MW VolturnUS-S in a turbulent 600 s DLC-1.6 environment at ``dtM = 0.1 s``
   (four glue steps of 0.025 s), the platform surge, sway, roll and pitch time series differ from
   the ``dtM = DT`` run by at most 0.12 % of their standard deviation with the extrapolation, on
   the chain mooring as on the chain and power-cable mooring, where the held loads give up to
   0.42 %. The chain-only run took a fifth less time (116 s against 141 s for ``dtM = DT``
   on the same machine). In the three-seed DLC-1.1 comparison
   of :doc:`validation`, made with held loads, the platform motions with ``dtM`` above ``DT``
   agree with MoorDyn within 0.6 % in mean and standard deviation, but the fairlead tension's
   standard deviation differs by about 3 % from CableDyn's own solution at every glue step.

   The extrapolation does not change how the lines themselves are integrated: each line still
   steps on ``dtM``, so its own resolution (see the ``dtM`` advisory in :doc:`modeling`) decides
   its tension accuracy, and the output channels and the ``.CD.out`` record the committed
   samples. For stiff taut, semi-taut or tension-leg moorings, and for any fatigue or
   tension-range result, check a smaller ``dtM`` (down to ``dtM = DT``) before relying on it.

Fluid source
~~~~~~~~~~~~

With no WaterKin policy (``0 WaterKin``, the default), CableDyn samples the complete enabled
SeaState wave and current field at its line nodes and refreshes it every CableDyn step. A
WaterKin file can select file waves, a file depth-dependent current, host waves, or host
current independently where the field is separable. Unsupported or double-counted combinations
fail during initialisation. The selector matrix is in :doc:`options` and :doc:`driver_format`.

SeaState grid
~~~~~~~~~~~~~

SeaState computes the wave and current kinematics on a grid set in its input file and
interpolates them at the line nodes: ``X_HalfWidth`` and ``Y_HalfWidth`` about the origin,
``Z_Depth`` below the still-water level, and ``NX``, ``NY``, ``NZ`` grid points (``NZ`` is
spaced more finely near the surface). A node outside the grid gets the kinematics of the
nearest grid boundary point, so its fluid load is wrong: beyond a half-width the waves are
out of phase, and below ``Z_Depth`` the wave and current kinematics are those of the grid
bottom, which are too strong. A SeaState file written for the platform alone usually covers
only the platform. The r-test ``MD_Shared`` file used in :doc:`tutorial_openfast` has
``X_HalfWidth = Y_HalfWidth = 100 m``, ``Z_Depth = 40 m`` and ``NZ = 2``, while the anchors lie
838 m from the platform on a 200 m seabed. Size the grid for the lines:

* ``X_HalfWidth`` and ``Y_HalfWidth`` beyond the farthest anchor or cable end, with margin for
  the platform offset;
* ``Z_Depth`` equal to the water depth (``"default"``); nodes on the seabed are then inside the
  grid; and
* ``NX``, ``NY`` and ``NZ`` raised with the extents to keep the grid spacing: a horizontal
  spacing of about a tenth of the peak wavelength or less. A long-crested sea (``WaveDirMod =
  0``) that travels along x does not vary in y, so ``NY = 2`` is enough for it.

SeaState holds the whole time history of the grid in memory, so memory and initialisation time
grow with (2 ``NX`` − 1)(2 ``NY`` − 1) ``NZ`` ``WaveTMax`` / ``WaveDT``. For the
IEA-15MW case of :doc:`tutorial_openfast` (Hs 6 m, Tp 12 s, ``WaveTMax = 4800 s``,
``WaveDT = 0.2 s``), ``X_HalfWidth = Y_HalfWidth = 900``, the default ``Z_Depth``, ``NX = 37``
(25 m spacing), ``NY = 2`` and ``NZ = 11`` needs about 4.6 GB. Halving the spacing (``NX = 73``)
or using ``NZ = 16`` changes the 10–70 s platform and fairlead-tension statistics by less than
0.1 %. The ``MD_Shared`` grid gives a line-1 fairlead tension standard deviation 18 % too high
(0.136 MN instead of 0.115 MN).

CableDyn compares its nodes with the grid at the static equilibrium and at every CableDyn step.
At the first node found outside it prints one warning and continues:

.. code-block:: text

   WARNING: CableDyn mooring module: 147 of the 153 CableDyn sampling nodes in the water lie
   outside the SeaState wave-kinematics grid at t = 0 s (warned once). SeaState clamps them to
   the grid boundary, so their fluid kinematics and loads are not physical. Farthest: x = -837.6
   m (grid |x| <= X_HalfWidth = 100 m); y = 725.38 m (grid |y| <= Y_HalfWidth = 100 m); z =
   -200.01 m (grid z >= -Z_Depth = -40 m). In the SeaState input file enlarge X_HalfWidth,
   Y_HalfWidth and Z_Depth (Z_Depth at most WtrDpth + MSL2SWL, its default) to cover the lines
   with margin for their motion, and raise NX, NY and NZ to keep the grid spacing.

Only nodes in the water are checked: a node above the free surface carries no fluid load,
whatever kinematics SeaState would give it, so a suspended cable or jumper above the water may
leave the grid. A node counts as in the water at or below the free-surface elevation sampled
for it, which also decides its wetting in CableDyn (the still-water level when the host waves
are off), and, without wave stretching (``WaveStMod = 0``), at or below the still-water level,
above which SeaState gives zero kinematics without reading the grid. Only the directions in
which nodes leave the grid are listed. Waves are checked in all three
directions and a SeaState current (uniform in x and y) only in depth; under Wheeler stretching
(``WaveStMod = 3``) the depth checked is the one SeaState interpolates at, the node depth mapped
linearly from the range between the seabed and the local wave elevation onto the range between
the seabed and the still-water level (with the elevation SeaState sampled, also when only the
host current is used). A dynamic current field (an InflowWind-driven current) is read on its
own grid; when that field is gridded, a separate one-time warning reports nodes outside its
vertical extent or lateral half-width. SeaState measures its field from the still-water level, and the
mooring samples it at positions relative to the mean sea level, so a run with host waves or
current from SeaState refuses a nonzero ``MSL2SWL`` at initialisation (a deck that takes
both from WaterKin does not query SeaState and is not affected). Excursions smaller than
1 % of the grid spacing at the boundary are ignored, and so are nodes below a grid that reaches
the seabed (``Z_Depth = WtrDpth + MSL2SWL``): they are in the seabed, where the kinematics of the
grid bottom are the right ones. FAST.Farm runs, which use still water, and SeaState fields
without waves or current are not checked.

Run and confirm
---------------

.. code-block:: powershell

   C:\CableDyn\openfast.exe IEA-15-UMaine_CompMooring5_CableDyn.fst
   if ($LASTEXITCODE -ne 0) { throw "OpenFAST failed with $LASTEXITCODE" }

Before trusting a run, confirm in the console:

#. ``Running CableDyn (v0.2.0, ...)`` and ``Parsing CableDyn input file: <deck>``;
#. the ``dtM`` line and the model inventory (``Created CableDyn model: ... line object(s)``);
#. for every line, the converged fairlead effective tension, force vector, inclination,
   declination, and azimuth at the static equilibrium — printed unconditionally, whether or not
   ``OUTPUTS`` is present — followed by the line's static-initialisation notes, if any
   (``Note: line ...``, for example a folded ``EI = 0`` layout; see :doc:`troubleshooting`);
#. ``Requested CableDyn OUTPUTS at t = 0 s``, the deck's channels in the first output row. With
   SeaState waves or current the header adds ``SeaState kinematics at t = 0``: the equilibrium is
   solved before the kinematics act, so a line-end tension at t = 0 can differ slightly from the
   equilibrium value above by the kinematic load on the end node. ``<root>.CD.static.out`` holds
   the equilibrium; and
#. ``OpenFAST terminated normally.``

The force is the load the line exerts on its End A (the fairlead), directed toward End B; the
angle definitions are in :ref:`conventions:Angles at line ends`. A tutorial run with the full
console text is in :doc:`tutorial_openfast`.

Output channels and files
-------------------------

CableDyn channels are selected in the deck's ``OUTPUTS`` section, one quoted channel per row,
closed by ``END`` — not in the ``.fst`` ``OutList``:

.. code-block:: text

   --------------------------- OUTPUTS -------------------------------------------
   "FairTen1"
   "AnchTen1"
   "FairIncl1"
   "AnchIncl1"
   END

The channel vocabulary (tensions, end angles, node tension, curvature, bend moment, position,
velocity, acceleration, declination, azimuth, the touchdown point, and the anchor uplift
channels ``AnchFz<L>``, ``AnchFH<L>`` and ``TDP<L>Grd``) and units are in :doc:`outputs`.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - File
     - Contents
   * - ``<RootName>.out`` / ``.outb``
     - the OpenFAST table at ``DT_Out``, including the CableDyn channels (held between CableDyn
       commits). With ``TStart = 0`` its first row equals the printed t = 0 ``OUTPUTS`` values.
   * - ``<RootName>.CD.out``
     - CableDyn-owned: the ``OUTPUTS`` channels at ``t = 0`` and once per committed ``dtM``
       solve, never repeated glue-rate values
   * - ``<RootName>.CD.static.out``
     - CableDyn-owned, always written: the static-equilibrium profile of every line, one row per
       node (arc length, coordinates, effective tension, curvature, bend moment, declination,
       inclination, azimuth) — the coupled range graph, at the state of the printed line summary
   * - ``<RootName>.CD.rst.dat``, ``<RootName>.CD.lin.dat``
     - temporary CableDyn-owned copies of the deck, written in the ``MooringFile`` folder only
       when a run restarts from a checkpoint or linearises (see *Lifecycle features*), and
       removed after use

With ``WrVTK > 0`` OpenFAST writes its VTK visualisation files as usual, but CableDyn provides no
mooring or cable meshes to them, and says so at initialisation: the lines do not appear in the
VTK output. Use the node-position
channels, the ``.CD.static.out`` profile or the Python snapshot tools (:doc:`python`) to view
line geometry.

Example files
-------------

.. list-table::
   :header-rows: 1
   :widths: 48 52

   * - File (in ``examples/``)
     - Purpose
   * - ``openfast/IEA-15-UMaine_CompMooring5_CableDyn.fst``
     - IEA-15MW UMaine VolturnUS-S template selecting CableDyn
   * - ``openfast/CableDyn_UMaine.dat``
     - three-line R4 chain mooring deck, ``dtM = DT``
   * - ``openfast/IEA-15-UMaine_CompMooring3_MoorDyn.fst`` + ``openfast/MoorDyn_UMaine.dat``
     - the like-for-like stock-MoorDyn twin; ``SeaState WaterKin`` gives its lines the same
       SeaState wave kinematics CableDyn samples
   * - ``iea15mw_umaine_mixed_cabledyn.dat``
     - the three chains plus one finite-EI lazy-wave power cable with a grounded tail
   * - ``iea15mw_umaine_openfast_cabledyn.dat``
     - an alternative OpenFAST-only mooring deck for the same platform

The turbine sub-models are not part of this repository. :doc:`tutorial_openfast` gives the exact
steps to fetch them from OpenFAST's ``r-test`` (``glue-codes/fast-farm/MD_Shared`` at revision
``dd5feaaaa500ba7283140107806300d551cff0a7``), to reset that farm model's turbine-1 initial
pose to a single turbine at the origin, and to enlarge its SeaState grid to cover the moorings
(see *SeaState grid*).

Mixed moorings and power cables
-------------------------------

One deck may combine ``EI = 0`` chains and finite-EI (``EI > 0``) power cables; each line uses
its own element. In ``iea15mw_umaine_mixed_cabledyn.dat`` line 4 runs from a platform hang-off
to a fixed termination on the 200 m seabed; its bare/buoyant/bare sections form the lazy wave and
a grounded tail (149 of 249 nodes on the seabed at equilibrium). ``WtrDpth``, ``kBot``, ``cBot``,
and ``frictionMu`` define contact; a structured ``bathymetryFile`` may replace the flat floor.
CableDyn finds touchdown during static initialisation. Inspect ``<RootName>.CD.static.out``
(curvature, grounded length, section transitions, minimum bend radius) before the time history.

Line-end bending connections
----------------------------

A finite-EI line may declare an ``END CONNECTIONS`` row at End A or End B. A numeric stiffness
defines an isotropic rotational spring, ``Pinned`` leaves the tangent direction free, and
``Rigid`` enforces the declared direction while keeping the axial stretch free. At a coupled
hang-off the declared direction is given in platform axes and follows the OpenFAST point-mesh
orientation; its angular velocity and acceleration enter the implicit residual and tangent, and
the connection moment is returned to OpenFAST. Fixed-end directions are global. A rigid
direction reversal of exactly 180° in one host update is rejected, because the shortest rotation
is not unique. Formal linearisation (``Linearize = True``) with a platform-relative end
connection fails during initialisation. Row grammar: :doc:`driver_format`.

Line torsion
------------

A finite-EI line restrained in torsion at both ends (``END CONNECTIONS`` torsion columns, see
:ref:`theory-torsion`) runs coupled as it does standalone. Its End A frame is given in platform
axes and turns with the OpenFAST point-mesh orientation, whatever the bending connection, so a
platform roll, pitch or yaw about the line's end direction twists the line. The torque moment is
returned to OpenFAST with the connection moment, the ``Torq<L>N<J>`` and ``Twist<L>``
channels are available, ``PtfmInit`` twists the line in the static initialisation, and the
unwrapped twist is part of the checkpoint state. End B must be a ``Fixed`` point. Within one
step the line may turn by more than 90° about its end direction, because the step is cut into
substeps; an output probe moves the frame without a step, is limited to 90° from the
committed twist, and stops by name beyond it (a restart first restores the checkpoint's twist,
so it is not limited this way). Formal
linearisation of a line restrained in torsion fails during initialisation: the torque's
rotational stiffness has no column in the translation-only quasi-static reduction. The C API
and the standalone deck that mixes ``EI = 0`` and finite-EI lines without a body move the line
ends by translation only and refuse torsion columns.

Lifecycle features
------------------

The same module takes part in OpenFAST's checkpoint, linearisation, farm and control operations
without an auxiliary process:

* **Checkpoint/restart** serialises the committed line, constitutive, connection, body, and rod
  state needed to continue the trajectory, with the SeaState fields the committed step holds and,
  with ``dtM`` larger than ``DT``, the loads of the solve before it (the load predictor), so a
  restart continues the run bit for bit at any ``dtM``. A checkpoint is accepted only on a committed CableDyn
  step, so choose ``ChkptTime`` as a whole multiple of ``dtM``; a checkpoint written between two
  commits stops the restart with a message that says so. On restart the module writes a
  temporary deck copy, ``<RootName>.CD.rst.dat``, beside the ``MooringFile``, rebuilds the model
  from it and removes it; files the deck names (``bathymetryFile``, a WaterKin file, Syrope
  tables) resolve against the ``MooringFile`` folder as in the original run. Keep restart files with the exact executable and complete input set.
* **Linearisation** publishes the quasi-static mooring load derivative at a committed boundary.
  It rebuilds the model from a temporary deck copy, ``<RootName>.CD.lin.dat``, written beside the
  ``MooringFile`` and removed after use, so side files resolve as on restart. Set ``dtM = DT`` and linearise at the initial
  operating point (``CalcSteady = False``,
  ``LinTimes = 0``, platform at its still-water equilibrium, rotor speed and pitch at their
  steady values). Contact switching makes derivatives near touchdown non-smooth. OpenFAST writes
  the ``.lin`` matrices in the ``.fst`` ``OutFmt``; to compare them with MoorDyn's
  ``D - C A^-1 B``, use a full-precision format such as ``"ES24.15E3"``, because the usual
  4-digit format cannot resolve that reduction. A deck with
  Rigid6 bodies, rods or ``Point3`` buoys is rejected at initialisation: the quasi-static
  reduction would hold them at their deck pose instead of re-solving their equilibrium (a
  ``Free`` point with the buoy's ``Mass``, ``Vol``, ``CdA`` and ``Ca`` is re-solved).
* **FAST.Farm** selects CableDyn with ``Mod_SharedMooring = 5`` for shared lines and
  connections. Gravity, water density, depth, and ``TMax`` are host-owned and must agree across
  turbines, otherwise initialisation fails.
* **Line failures** detach a line end on a time or tension trigger and replay correctly after a
  restart. Verify the trigger and the post-failure topology in a dedicated case first.
* **Active tensioning** maps ServoDyn ``CableDeltaL``/``CableDeltaLdot`` control onto the
  fairlead-side segment of ``EI = 0`` lines and finite-EI cables, the segment at the line's end
  on a host-driven point, or at the ``NodeB`` written in the deck when both or neither ends are
  (:doc:`driver_format`, CONTROL). On a finite-EI cable the command scales every cable element
  of that segment, and needs the force-blended generalised-α (the default). On both kinds a
  command outside MoorDyn-F's bounds, from minus half to plus one ``UnstrLen``/N of the line
  (the controlled segment's initial length on a line of one section), stops the run with a
  message that names the line, as it stops MoorDyn-F; so does a command that would make the
  controlled segment of a line of several sections non-positive. FAST.Farm decks reject
  CONTROL.
* **Rigid6 bodies, Point3 bodies, and rigid rods** are translated at the deck boundary; their
  supported loads and restrictions in coupled runs are listed in :doc:`capabilities`. Free
  Rigid6 bodies, free rods and pinned rods start at their static equilibrium with the platform
  at ``PtfmInit`` (``bodyIC static``, the default), as in standalone runs; one that carries no
  line is not part of that solve and starts at rest in its deck pose. ``Coupled``/``Vessel``
  rods are platform-borne: each is a node of the coupling mesh at its End A that follows the
  platform position and orientation and returns a force and a moment (MoorDyn-F's coupled rod,
  including the rod's own and added-mass inertia). ``CoupledPinned``/``VesselPinned`` rods are
  pinned to the platform: the node at End A follows the platform position only, the rod turns
  about it with its own rotation (the node orientation and angular rates are ignored), and the
  node returns a force with zero moment (MoorDyn-F's coupled pinned rod). That force is the
  rod's net load less its inertia, with the angular acceleration solved for the node's current
  acceleration; a CalcOutput at an off-step acceleration therefore returns the rod's effective
  inertia about the pin rather than its full mass (the rotation, the angular velocity and the
  held fluid stay those of the committed step). ``PtfmInit`` moves and turns such a rod with
  the platform before the static solve balances it about its pin. ``Pinned`` rods turn about a
  fixed End A. ``Coupled``/``Vessel`` bodies are platform-borne the same way as ``Coupled`` rods,
  with the node at the body reference point, and may share a deck with ``Free`` bodies.
  ``CoupledPinned`` bodies are not supported on this route.

Each is a separate qualification: a successful ordinary run does not qualify restart,
linearisation, control, failure, or farm behaviour for a new model.

Two different Jacobians
~~~~~~~~~~~~~~~~~~~~~~~

During time marching, the OpenFAST input–output solver asks for the instantaneous output
derivative with the cable state frozen; CableDyn evaluates it by symmetric perturbation of the
coupled input while preserving the committed state. In a ``Linearize = True`` run CableDyn
registers no dynamic states and reports the re-equilibrated quasi-static stiffness instead.
These are different derivatives and must not be interchanged. ``NumCrctn`` controls glue
correction passes; it does not change which derivative is used.

Performance
-----------

Measure coupled performance with the complete OpenFAST model, not with standalone line runs:
OpenFAST coupling iterations, SeaState interpolation, controller, and output work are part of
the cost. Report the ``Time Ratio (Sim/CPU)`` from a run that **terminated normally**, with
identical executable, inputs, thread counts, wind/wave realisation, and output channels between
the cases you compare, and with the default solver settings (full Newton, default dynamic
tolerances, ``rhoInf`` default of :ref:`conventions:Time integration parameter`, no adaptive
remeshing). ``dtM`` is the main cost lever; choose it by convergence, not speed. When
OpenFAST/CableDyn OpenMP threading is active, keep the BLAS library single-threaded, and compare
against one thread: small line sets can lose time to parallel overhead.

Common mistakes
---------------

* **CompMooring rejected during input validation:** a stock ``openfast.exe`` is being used.
* **Mooring file not found:** the path is resolved from the ``.fst`` folder, not from the
  executable's folder.
* **Missing DISCON.dll:** the ServoDyn input requests a controller; obtain the matching DLL or
  set ``CompServo = 0`` for a controller-free run.
* **Missing turbine files:** the repository provides mooring decks and templates, not the turbine
  model; see :doc:`tutorial_openfast`.
* **Implausible initial tensions** (for example hundreds of MN): the platform's initial pose in
  ElastoDyn does not match the fairlead and anchor layout of the deck.
* **Initialisation fails closed:** read the complete message and consult :doc:`troubleshooting`;
  do not bypass a rejected feature combination.
