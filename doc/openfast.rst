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
CableDyn targets 0.1 s rounded to a whole number of glue steps and holds its loads (zero-order
hold) between commits. ``dtM = DT`` solves every glue step. The console reports the result, for
example ``CableDyn time step dtM = 2.50000E-02 s (1 x glue DT; deck dtM )``. A deck ``dtM``
larger than ``DT`` but not a whole multiple of it is never coarsened: CableDyn uses the largest
multiple of ``DT`` that does not exceed it and reports a warning with both values (``dtM = 0.025``
with ``DT = 0.01`` runs at 0.02 s), matching MoorDyn, which only ever shortens its step to fit the
coupling interval. A deck ``dtM`` below ``DT`` (typical of MoorDyn decks written for explicit
integration) runs at ``dtM = DT``, reported on the console. ``TMax`` must end
on a commit boundary. For fatigue, snap, contact, or local finite-EI response compare at least
two ``dtM`` values (for example 0.1, 0.05, and 0.025 s) and demonstrate convergence of ranges,
peaks, and damage. If CableDyn needs a smaller step than the glue step, reduce ``DT`` as well.
``DT_Out`` changes only the sampling of the main OpenFAST table.

Fluid source
~~~~~~~~~~~~

With no WaterKin policy (``0 WaterKin``, the default), CableDyn samples the complete enabled
SeaState wave and current field at its line nodes and refreshes it every CableDyn step. A
WaterKin file can select file waves, a file depth-dependent current, host waves, or host
current independently where the field is separable. Unsupported or double-counted combinations
fail during initialisation. The selector matrix is in :doc:`options` and :doc:`driver_format`.

Run and confirm
---------------

.. code-block:: powershell

   C:\CableDyn\openfast.exe IEA-15-UMaine_CompMooring5_CableDyn.fst
   if ($LASTEXITCODE -ne 0) { throw "OpenFAST failed with $LASTEXITCODE" }

Before trusting a run, confirm in the console:

#. ``Running CableDyn (v0.1.1, ...)`` and ``Parsing CableDyn input file: <deck>``;
#. the ``dtM`` line and the model inventory (``Created CableDyn model: ... line object(s)``);
#. for every line, the converged fairlead effective tension, force vector, inclination,
   declination, and azimuth at the static equilibrium — printed unconditionally, whether or not
   ``OUTPUTS`` is present;
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
velocity, acceleration, declination, azimuth) and units are in :doc:`outputs`.

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
``dd5feaaaa500ba7283140107806300d551cff0a7``) and to reset that farm model's turbine-1 initial
pose to a single turbine at the origin.

Mixed moorings and power cables
-------------------------------

One deck may combine ``EI = 0`` chains and finite-EI (``EI > 0``) power cables; each line uses
its own element. In ``iea15mw_umaine_mixed_cabledyn.dat`` line 4 runs from a platform hang-off
to a fixed termination on the 200 m seabed; its bare/buoyant/bare sections form the lazy wave and
a grounded tail (38 of 88 nodes on the seabed at equilibrium). ``WtrDpth``, ``kBot``, ``cBot``,
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

Lifecycle features
------------------

The same module takes part in OpenFAST's checkpoint, linearisation, farm and control operations
without an auxiliary process:

* **Checkpoint/restart** serialises the committed line, constitutive, connection, body, and rod
  state needed to continue the trajectory. A checkpoint is accepted only on a committed CableDyn
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
  steady values). Contact switching makes derivatives near touchdown non-smooth. A deck with
  Rigid6 bodies or rods is rejected at initialisation: the quasi-static reduction would hold
  them at their deck pose instead of re-solving their equilibrium.
* **FAST.Farm** selects CableDyn with ``Mod_SharedMooring = 5`` for shared lines and
  connections. Gravity, water density, depth, and ``TMax`` are host-owned and must agree across
  turbines, otherwise initialisation fails.
* **Line failures** detach a line end on a time or tension trigger and replay correctly after a
  restart. Verify the trigger and the post-failure topology in a dedicated case first.
* **Active tensioning** maps ServoDyn ``CableDeltaL`` control onto ``EI = 0`` lines; finite-EI
  lines reject it.
* **Rigid6 bodies, Point3 bodies, and rigid rods** are translated at the deck boundary; their
  supported loads and restrictions in coupled runs are listed in :doc:`capabilities`. Free
  Rigid6 bodies and rods start at their static equilibrium with the platform at ``PtfmInit``
  (``bodyIC static``, the default), as in standalone runs. ``Coupled``/``Vessel`` rods are
  platform-borne: each is a node of the coupling mesh at its End A that follows the platform
  position and orientation and returns a force and a moment (MoorDyn-F's coupled rod, including
  the rod's own and added-mass inertia). ``Coupled``/``Vessel`` bodies are platform-borne the
  same way, with the node at the body reference point, and may share a deck with ``Free``
  bodies. ``CoupledPinned`` rods and bodies, and ``Pinned`` rods, are not supported on this route.

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
