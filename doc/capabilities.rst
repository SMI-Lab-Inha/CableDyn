.. SPDX-License-Identifier: Apache-2.0

Capabilities and route selection
================================

CableDyn deliberately exposes several execution routes over one solver core. A feature can be
implemented in the core yet unavailable on a particular route because that caller does not provide
the required motion, fluid field, or restart state. Use this page to choose the route before
writing a deck. The evidence behind every checked item is recorded in :doc:`validation`; rejected
combinations fail closed with a named error.

Execution routes
----------------

.. list-table::
   :header-rows: 1
   :widths: 19 25 30 26

   * - Route
     - Boundary owner
     - Best use
     - Start here
   * - Released ``CableDyn_driver.exe``
     - deck coordinates, held ends, ``motionFile``, ``vesselMotion``, or ``vesselRAO``
     - static design, property checkout, prescribed-motion dynamics, waves/current studies
     - :doc:`standalone_driver`
   * - Source-build ``cabledyn``
       (``build/cabledyn``; ``build\bin\cabledyn.exe`` with the conda toolchain on Windows)
     - same deck route as the release driver
     - source builds, tests, and research automation
     - :doc:`installation`
   * - CableDyn-enabled ``openfast.exe`` (OpenFAST, maintained by NLR, the National Laboratory of
       the Rockies, formerly NREL)
     - OpenFAST platform/body motion and SeaState
     - full turbine simulations, DLCs, mixed mooring/power-cable systems, FAST.Farm
     - :doc:`openfast`
   * - Python executable wrapper
     - a Python process launching the standalone driver
     - parameter sweeps, managed output, post-processing
     - :doc:`python`
   * - C ABI / in-process Python API
     - the embedding program
     - CFD and custom co-simulation with explicit kinematics-in/loads-out ownership
     - :doc:`capi`

Line and constitutive capability
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Model
     - Status
     - Scope
   * - Linear ``EI = 0`` chain, wire, polyester, nylon
     - Supported
     - static and implicit dynamics; taut, semi-taut, and grounded catenary configurations
   * - Composite line
     - Supported
     - one line assembled from ordered sections with independent type, length, and mesh density
   * - Finite-EI bending cable
     - Supported
     - production cubic-Hermite static/dynamic path for suspended and flat/structured
       seabed-contact cables in standalone and OpenFAST-coupled use; the uncommon
       two-moving-end standalone topology uses a separate compatibility solver
   * - Mixed ``EI = 0`` + finite-EI deck
     - Supported with route limits
     - a mixed deck of lines between held or coupled points only (no bodies, rods, or
       ``Free``/``Connect`` points) runs in OpenFAST and standalone as a static or held-end
       dynamic run, and rejects a standalone ``motionFile`` or deck waves; a mixed deck with
       bodies, rods, or ``Free``/``Connect`` points runs standalone on the multibody march, with
       moving objects and deck waves/current, but not with a ``motionFile``
   * - Viscoelastic rope (MoorDyn ``ElasticMod`` 2/3)
     - Supported
     - per-element series-Kelvin state (standard linear solid only when ``BA_D = 0``),
       including load-dependent dynamic stiffness; ``EI = 0`` lines
   * - Syrope polyester
     - Supported with limits
     - single-section taut dynamic line with OWC and two history states; composite, hydro,
       finite-EI, current, and wave combinations fail closed by name
   * - Torsion, hockling, nonlinear cross-section laws
     - Not provided
     - the production elements carry no torsion; the secondary Cosserat path
       (:doc:`solver_paths`) is not a production route
   * - Discrete attachments (``ATTACHMENTS``: buoyancy modules, clumps)
     - Supported
     - lumped at nodes of a finite-EI cable, standalone and coupled; rejected on the two-moving-end
       compatibility solver
   * - Modal analysis (``nModes``)
     - Standalone only
     - natural frequencies and mode shapes of each line about its static equilibrium (all
       ``EI = 0`` lines, or all finite-EI lines on the Hermite route; flat seabed)
   * - VIV (``compVIV``)
     - Not implemented
     - use another validated model; CableDyn rejects the request

Loads and environment
---------------------

.. list-table::
   :header-rows: 1
   :widths: 28 20 52

   * - Feature
     - Status
     - Important ownership rule
   * - Gravity and displaced-volume buoyancy
     - Supported
     - SI input; water density may come from the deck or OpenFAST environment
   * - Morison drag, added mass, Froude–Krylov
     - Supported
     - line and supported rigid-object fields use the same nodewise fluid contract
   * - Flat seabed contact and friction
     - Supported
     - declare ``WtrDpth``, ``kBot``, ``cBot`` and optional friction; supported for coupled
       finite-EI touchdown cables as well as EI = 0 lines
   * - Anisotropic seabed friction
     - Supported with limits
     - ``frictionMuAxial``/``frictionMuLateral`` (OrcaFlex axial and normal coefficients), within
       the friction limits listed in :doc:`options`
   * - Structured bathymetry
     - Supported with limits
     - mutually exclusive with flat ``WtrDpth``; coupled finite-EI cables query the global
       surface through their chord frame; see :doc:`driver_format`
   * - Uniform/profile current
     - Supported
     - standalone deck source; coupled runs use explicitly selected host/file sources
   * - Airy and JONSWAP waves
     - Supported
     - standalone deck source or OpenFAST SeaState, never an accidental double source
   * - Regular nonlinear (Dean stream-function) waves
     - Standalone only
     - deck ``stream H T dir waves``, ``StreamOrder``
   * - ISSC/Pierson-Moskowitz, Torsethaugen and Ochi-Hubble spectra, cos-2s spreading,
       multiple wave trains
     - Standalone only
     - deck ``waves``/``wavetrain`` rows with ``WaveSpreading``; see :doc:`options`
   * - Prescribed vessel motion (``vesselMotion``, ``vesselRAO``)
     - Standalone only
     - a 6-DOF vessel record, or the RAO response to the deck waves, moves every
       ``Coupled``/``Vessel`` point of a line deck; see :doc:`options`
   * - Range graphs (``r`` line output flag, ``RangeStart``)
     - Standalone only
     - minimum, maximum, and mean along each line over the run; a coupled OpenFAST deck rejects
       the flag
   * - WaterKin CurrentMod-1 file table
     - Supported
     - bit-identical to the equivalent inline depth profile
   * - WaterKin WaveKinMod-1 history
     - Standalone only
     - resampled over ``TMax``; caller-driven use fails closed because the aggregate has no
       self-driven per-step wave refresh
   * - WaterKin WaveKinMod-2 / in-file ``SEASTATE``
     - Coupled only
     - consumes the host SeaState through the OpenFAST shell and requires that field to exist

Topology and coupled objects
----------------------------

CableDyn's native objects are finite-element lines and boundary attachments. ``POINTS``, ``BODIES``,
and ``RODS`` are MoorDyn-compatible deck records that the driver translates; they are not an
alternative discretisation of the line. Supported translations include fixed/coupled/vessel
attachments, free or connecting masses and clumps, ``Point3`` buoys, ``Rigid6`` bodies, and coupled
rigid rods. In OpenFAST, ``Coupled``/``Vessel`` bodies and rods are platform-borne 6-DOF nodes that
return a force and a moment; the C API and Python run free bodies and rods in still water. The
coupled body/rod routes carry checkpoint state and SeaState drag, Froude–Krylov, and added mass.
``Rigid6`` also carries the documented hydrostatic restoring terms. Consult :doc:`driver_format`
before using an imported object because unsupported fields fail closed.

OpenFAST coverage
-----------------

The ``CompMooring = 5`` shell covers normal coupled time-domain stepping with its own ``dtM``,
nonzero ``PtfmInit``, correction iterations, SeaState fields, checkpoint/restart, scoped
quasi-static ``dYdu`` linearisation, line failures, active-tension control for ``EI = 0`` lines,
and FAST.Farm shared moorings. Restart and linearisation do not arise under FAST.Farm, which
does not call those module operations. Active ServoDyn control of a finite-EI cable
remains rejected. See :doc:`openfast` for setup and :doc:`coupling_boundary` for the interface
contract.

How to interpret this page
--------------------------

``Supported`` means the route is implemented and covered by regression tests, not that every
conceivable combination is valid. Physical source ownership still matters. When two requested
features cannot be separated without double counting—for example a host field whose private
current cannot be split from its wave field—CableDyn stops rather than guessing. The fatal message
and :doc:`troubleshooting` identify the conflicting sources.
