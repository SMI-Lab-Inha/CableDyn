.. SPDX-License-Identifier: Apache-2.0

Coupling: OpenFAST and CFD
==========================

CableDyn is one caller-agnostic solver core behind a single coupling boundary, with a thin
shell per target. The core knows nothing about its caller; that separation is the central
architectural decision. The formal contract — kinematics in, loads out, with frames, units,
and cadence — is specified in :doc:`coupling_boundary`.

OpenFAST module (``CompMooring = 5``)
-------------------------------------

CableDyn has a mooring module for OpenFAST v5 (OpenFAST is maintained by NLR, the National
Laboratory of the Rockies, formerly NREL), selected with ``CompMooring = 5``, that coexists with
stock MoorDyn (``CompMooring = 3``) in an OpenFAST build carrying the CableDyn module. The
IEA-15MW VolturnUS-S UMaine model runs coupled end to end with ``EI = 0`` catenary moorings and
a finite-EI dynamic power cable in one deck (:doc:`tutorial_openfast`; the comparisons are in
:doc:`validation`).

At initialisation the module announces itself the standard OpenFAST way, identifies the parsed deck
inventory, and prints an unconditional converged-static fairlead report for every line. The report
includes effective tension, the global force on the fairlead, and the inclination, declination, and
azimuth of the line end as defined in :doc:`conventions`. It is independent of the deck ``OUTPUTS``
selection, while requested channels are printed in a separate block. This makes the active model and
angle convention visible directly in an OpenFAST run log. The module covers MoorDyn-F's coupled
integration surface: coupled time-domain analysis at a 0.1 s CableDyn step by default (with deck
``dtM`` control), a nonzero initial platform displacement (``PtfmInit``), glue correction iterations
(``NumCrctn > 0``), SeaState wave/current kinematics on the moorings and finite-EI cables,
checkpoint-restart, quasi-static linearisation (``dYdu``), FAST.Farm farm-level shared moorings
(``Mod_SharedMooring = 5``), coupled ``Rigid6`` 6-DOF rigid bodies (buoyancy,
``C33``/``C44``/``C55`` hydrostatic restoring, restart state, and SeaState drag / Froude-Krylov /
added-mass), coupled ``Point3`` buoys (buoyancy and SeaState drag / Froude-Krylov / added-mass), and
coupled rigid rods with host-driven endpoints, correction rewind, checkpoint state, and two-point
Gauss SeaState sampling.

The input/output Jacobian has two deliberately different contracts. During an ordinary nonlinear
OpenFAST run, CableDyn returns the partial direct-feedthrough derivative at fixed committed internal
state, which is the derivative required by OpenFAST's algebraic input/output solve. During formal
linearisation, where CableDyn exposes a zero-state quasi-static representation, it instead returns
the re-equilibrated zero-frequency mooring stiffness. Substituting the latter into the runtime glue
solve is mathematically incorrect; changing ``NumCrctn`` does not change that derivative contract.

.. admonition:: Remaining integration surface
   :class: warning

   Outside that surface the module fails closed — a clear fatal error, never a silent wrong
   answer — on:

   - active ServoDyn control (``CableDeltaL``) of a finite-EI cable, or any ``CONTROL`` section
     in FAST.Farm (``EI = 0`` line control is supported in a single-turbine model),
   - coupled ``Rigid6`` bodies and rods in a FAST.Farm deck,
   - deck ``waves``/``wavetrain`` OPTIONS (the host SeaState supplies the waves), and a deck
     ``current`` except as a steady current on a single-turbine, pure ``EI = 0`` deck without
     Rigid6 bodies or rods in a SeaState without waves or current (a WaterKin ``CurrentMod 1``
     table remains available),
   - VIV (``compVIV``), and
   - fluid-source combinations that the nodewise SeaState field cannot separate without double
     counting.

   ``Point3`` buoyancy and fluid loads are supported, but an imported restoring model that is not
   represented by those terms remains rejected. See :doc:`capabilities` for the route matrix. For
   an unsupported workflow, use stock MoorDyn (``CompMooring = 3``) in the same binary.

Standalone binaries
~~~~~~~~~~~~~~~~~~~~

The Windows release carries ``CableDyn_driver.exe`` and a CableDyn-enabled
``openfast.exe`` as non-OpenMP, statically linked x64 files. Neither needs an adjacent
compiler, BLAS, OpenMP, or C/C++ runtime DLL. A model can still request an external
ServoDyn controller such as ``DISCON.dll``; that is a model-specific plugin rather than a
dependency of CableDyn or ``openfast.exe``. The release build is described in
:doc:`installation`; end-user setup and a complete coupled walkthrough are in :doc:`openfast`.
The release compiler is pinned to IFX 2025.3 or newer, and both CableDyn and stock-MoorDyn
OpenFAST paths must pass separate smoke cases. A successful CableDyn case alone does not qualify
the executable as a drop-in OpenFAST replacement. The packaged OpenFAST executable is built from
the upstream Visual Studio ``Release|x64`` solution target (the same build family used for the
official Windows release), with CableDyn installed as a normal solution module; it is not a
CMake-built approximation of that executable.

C binding (``CableDyn_CAPI``)
-----------------------------

The second shell is a standalone C-style binding (``CableDyn_CAPI``), couplable to CFD solvers
such as STAR-CCM+ and OpenFOAM the way MoorDyn-C is. It consumes the same coupling boundary as
the OpenFAST shell, so the two shells are one core with a thin adapter per target, not two
codebases.

The coupling boundary
---------------------

Both shells pass coupled-point kinematics and nodal or point fluid fields in and structural
loads out across one defined boundary, so the line's hydrodynamic response is computed
identically regardless of caller. Loads follow the sign convention of :doc:`conventions` (force
exerted by the cable on the coupled object, global frame). The full contract is in
:doc:`coupling_boundary`.

For user-facing setup rather than architecture, use :doc:`openfast` (coupled execution),
:doc:`capi` (embedding), or :doc:`python` (automation and in-process access).
