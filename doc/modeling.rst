.. SPDX-License-Identifier: Apache-2.0

Modelling workflow
==================

This is the practical engineering workflow for building a CableDyn model from geometry and line
data. It complements the exact grammar in :doc:`driver_format`: this page explains *what to choose
and check*, while the format reference defines every accepted field.

1. Choose the owner of motion and water kinematics
---------------------------------------------------

Start with :doc:`capabilities`. Use the standalone driver when endpoints are fixed, held, or
prescribed (``motionFile``, ``vesselMotion``, or ``vesselRAO``) and the deck owns waves/current.
Use OpenFAST, maintained by NLR (National Laboratory of the Rockies, formerly NREL), when the
turbine or platform owns endpoint motion and SeaState. Use the C ABI when another solver owns
the coupling loop. Do not carry a standalone motion source into an OpenFAST deck; do not declare a duplicate
fluid source unless the selector combination is explicitly supported.

2. Assemble the minimum engineering data
----------------------------------------

For every line collect:

* global End A and End B coordinates in metres;
* unstretched length and section transition locations;
* outer hydrodynamic diameter, dry mass per unit length, displaced volume convention, ``EA``,
  ``BA`` and ``EI``;
* normal/tangential drag and added-mass coefficients with their provenance;
* water depth or bathymetry, water density, gravity, and contact/friction assumptions;
* the motion, current, and wave time bases for a dynamic run; and
* the output quantities and physical locations needed for acceptance.

Record whether stiffness values are static, dynamic, mean-load-dependent, or working-curve data.
Do not substitute a dynamic rope modulus for ``EA`` without documenting the loading regime.

3. Build topology in public End-A to End-B order
-------------------------------------------------

End A is the upper/fairlead end and End B is the lower/anchor end. A ``LINE`` references its two
attachments; its ``SECTIONS`` rows then run from End A toward End B. Keep a physical cable or
mooring leg as one line even when its construction changes. A chain–polyester–chain leg is one
line with three sections, not three line objects joined merely to represent material boundaries.

Use a connection object only when there is a real force-balanced attachment such as a clump,
float, junction, or body. This distinction makes output identities, failure behaviour, and mesh
convergence unambiguous.

4. Define line types and constitutive response
-----------------------------------------------

``EI = 0`` selects the positions-only mooring element. ``EI > 0`` selects bending-cable behaviour
on the production cubic-Hermite route. The ``EA`` and ``BA`` tokens also select synthetic-rope
models:

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - Response
     - Input form
     - Use
   * - Linear
     - ``EA`` and ``BA``/``-zeta``
     - chain, wire, or a documented secant-modulus rope approximation
   * - Viscoelastic (series-Kelvin)
     - ``Es|Ed`` and ``Bs|Bd``
     - polyester/nylon with distinct slow and wave-frequency stiffness
   * - Viscoelastic, load-dependent
     - ``Es|alphaMBL|vbeta`` and ``Bs|Bd``
     - MoorDyn ``ElasticMod = 3`` behaviour
   * - Syrope
     - ``SYROPE:settings|alpha|beta`` and ``BA_s|BA_d``
     - supported single-section taut dynamic polyester use case

Use :doc:`examples` for runnable instances and :doc:`theory` for the state equations (the
viscoelastic model is a four-parameter series-Kelvin solid that reduces to a standard linear
solid only when ``BA_D = 0``). Preserve
the units shown in the ``LINE TYPES`` table; the parser does not infer or convert vendor units.

5. Choose and converge the mesh
-------------------------------

``NumSegs`` is section-local. Place section boundaries at physical property discontinuities, then
refine regions with high curvature or rapidly varying load: hang-off and bend-stiffener zones,
buoyancy-module transitions, sag bends, touchdown, and short weighted attachments. Avoid using a
material boundary solely as a substitute for mesh refinement.

Perform a mesh study on engineering observables, not node identity. At minimum compare endpoint
tension and the peak/position of curvature or touchdown response. Refine until the change is below
the project tolerance. Cross-code comparisons should interpolate by normalised arc length as
described in :doc:`validation`.

**Curvature at a line end.** The curvature near a hang-off or other held end is a boundary-layer
quantity: it varies over a length of order :math:`\sqrt{EI/T}` and is far more sensitive to the
time step and the local mesh than the interior sag- or hog-bend peak. On the reference lazy-wave
cables in a 3 m, 12 s heave with a pinned hang-off, the sag-bend peak is converged at
``dtM = 0.05 s``, but the curvature 1 m from the hang-off converges only at ``dtM ≤ 0.025 s``
(800 m cable) and ``dtM ≤ 0.0125 s`` (200 m cable), with the mesh graded toward the end. At
``dtM = 0.05 s`` both generalised-α blends overstate it, by a factor 1.6–2.3 at 1 m from the end.
Once converged, the values 0.5, 1, and 2 m from the end change by less than 2 % under a further
two- or four-fold local refinement. When end curvature matters:

- halve ``dtM`` until the end curvature changes by less than a few percent, and grade the mesh
  toward the end (a short, finely meshed section next to it);
- report the curvature 1–2 m from a pinned end, never the value at the pinned node: a pinned end
  carries no moment, so the exact curvature there is zero, and the computed nodal value is a
  discretisation residue (at ``dtM = 0.05 s`` up to 37 times the converged value 1 m away);
- for bend-stiffener or hang-off design, model the stiffener as a section and the end as a
  clamped or rotational-spring ``END CONNECTIONS`` row (see :doc:`driver_format`); a pinned end
  does not represent the moment a stiffener carries.

**Snap loads on slack lines.** When a slack line snaps taut it excites axial waves up to the
segment scale, near :math:`2c/l` with :math:`c = \sqrt{EA/m}` and :math:`l` the segment length.
The implicit integrator stays stable at any ``dtM``, but the snap tension converges in ``dtM`` only
once ``dtM`` resolves those modes (roughly :math:`\omega_\text{max}\,\mathrm{dtM} \lesssim 1`). Above
that, the tension history changes by a few tenths of a percent between time-step levels without
settling; means, envelopes and body motion are unaffected. On a buoy leg meshed with 320 segments
the limit is about 30 µs, and at 80 segments about 0.125 ms. An explicit code resolves these modes
only because its stability limit forces a step of that size; the implicit solver of OrcaFlex
(Orcina) shows a
similar plateau (about 0.5 %). When snap-tension histories matter, use the coarsest mesh that
resolves the line shape and refine ``dtM`` toward :math:`l/(2c)`, rather than refining both.

6. Establish the static initial condition
-----------------------------------------

Begin with the smallest model that should equilibrate:

#. omit ``dtM``/``TMax`` for an ``EI = 0`` static-only checkout;
#. verify endpoint order, length, submerged weight, and seabed elevation;
#. request endpoint tensions and the static along-arc profile;
#. inspect contact, touchdown, curvature, and symmetry; then
#. add constitutive states, connections, environmental forcing, and dynamics one at a time.

CableDyn solves the static equilibrium directly with Newton's method (with continuation as a
fallback) rather than by artificial drag-scaled relaxation. A converged
solver is necessary but not sufficient: reject a result with the wrong touchdown side, impossible
tension gradient, penetrated seabed, broken symmetry, or unresolved curvature peak.

7. Configure dynamics and forcing
----------------------------------

Set ``dtM`` and ``TMax`` for standalone time marching. ``rhoInf`` controls high-frequency
generalised-alpha dissipation; use the documented default unless a time-step study justifies a
change. With the default ``modified_newton = False`` every Newton iteration rebuilds the
tangent, except that a standalone finite-EI step on smooth motion starts from the previous
step's factorised tangent (:doc:`theory`). ``modified_newton = True`` also reuses the tangent
within a step, with guarded refresh, and converges to the same residual tolerance. For
performance comparisons keep the default.

Choose ``dtM`` by convergence of peaks, phase, rainflow ranges, and accumulated damage—not only by
solver success. Numerical stability at a large implicit step is not evidence of fatigue accuracy.
In OpenFAST the default target is ``dtM = 0.1 s``; loads and assembled OpenFAST channels are
zero-order-held between CableDyn solves. The separate ``.CD.out`` contains only genuine committed
``dtM`` samples, while ``DT_Out`` changes the assembled OpenFAST recording cadence only. Neither
file convention substitutes for a time-step convergence study, which should cover the forcing
bandwidth, the shortest element transit time, contact switching, and constitutive relaxation
times. For prescribed motion, row 1 defines the actual :math:`t=0` position, velocity, and
acceleration and must be consistent with the intended initial state.

8. Request auditable outputs
-----------------------------

Use one quoted channel per row:

.. code-block:: text

   --------------------- OUTPUTS ------------------------------------------
   "FairTen1"
   "AnchTen1"
   "Curv1N20"
   "BendMom1N20"

For static design, retain ``<out_root>.static.out`` and plot tension/curvature against arc length.
For dynamics, request endpoint loads plus channels at every fatigue-critical region. Keep units in
the analysis script explicit: CableDyn tension is N. See :doc:`outputs` for the complete grammar.

9. Validate in layers
----------------------

Use an evidence ladder appropriate to the decision:

* analytical catenary or simple suspended-line checks for geometry and submerged weight;
* mesh/time-step studies for discretisation error;
* like-for-like MoorDyn or OrcaFlex comparisons on common physical observables;
* standalone prescribed-motion tests before full OpenFAST coupling; and
* coupled A/B runs with identical turbine, wind, wave, current, and output windows.

Archive the executable version/hash, complete transitive input set, resolved toolchain, command,
solver log, and acceptance thresholds. Do not claim validation from a visually plausible trace.

Use-case map
------------

.. list-table::
   :header-rows: 1
   :widths: 28 34 38

   * - Use case
     - Recommended starting deck
     - Key checks
   * - Grounded chain catenary
     - ``wd0050_chain.dat``
     - touchdown, fairlead/anchor tension, seabed penetration
   * - Taut or semi-taut synthetic mooring
     - ``polyester_catenary_mooring.dat`` or ``semitaut_chain_polyester.dat``
     - stiffness convention, pretension, axial strain
   * - Dynamic-stiffness rope
     - ``ve_polyester_dynamic_waves.dat``
     - slow/dynamic state, relaxation time, cyclic tension
   * - Syrope polyester
     - ``syrope_polyester_mooring.dat``
     - OWC/settings files, initial history, supported-combination limits
   * - Current/wave-loaded mooring
     - ``dynamic_chain_current.dat`` / ``dynamic_chain_waves.dat``
     - source ownership, water depth, time-step and phase convergence
   * - Clump or line junction
     - ``clump_weight_free_point.dat`` / ``connect_weighted_point.dat``
     - attachment force balance and motion
   * - Composite chain/wire or chain/rope leg
     - ``composite_chain_wire.dat`` / ``composite_chain_poly_chain.dat``
     - section ordering and transition mesh
   * - IEA-15MW spread mooring
     - ``spread_3line_chain.dat`` (:doc:`tutorial_spread`)
     - symmetry, platform-restoring response
   * - Single VolturnUS-S mooring line
     - ``iea15mw_volturnus_mooring.dat``
     - single-line MoorDyn/OrcaFlex comparison
   * - Mixed mooring + finite-EI power cable
     - ``iea15mw_umaine_mixed_cabledyn.dat``
     - OpenFAST route, 200 m touchdown/contact, SeaState ownership, curvature and platform loads

Study layout and reproducibility
--------------------------------

Keep immutable inputs separate from generated results:

.. code-block:: text

   study/
   ├── inputs/       # decks and every transitively referenced file
   ├── scripts/      # commands, manifests, and post-processing
   ├── reference/    # analytical/cross-code evidence with provenance
   └── results/      # regenerated outputs and logs

Version-control inputs and scripts. Record ``CableDyn_driver.exe --version`` (and preferably the
release SHA-256) with every result. A conda environment file is useful but does not freeze resolved
package builds; archive ``conda list`` or the compiler/library versions when bitwise reproduction
matters.
