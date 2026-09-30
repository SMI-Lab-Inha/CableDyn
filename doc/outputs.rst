.. SPDX-License-Identifier: Apache-2.0

Output files and channels
=========================

This page is the reference for every result file CableDyn writes, the exact layout of each, and
the complete vocabulary of the deck ``OUTPUTS`` section. It covers the standalone driver
(``CableDyn_driver <deck.dat> <out_root>``, see :doc:`cli`) and the coupled module of OpenFAST,
maintained by NLR (National Laboratory of the Rockies, formerly NREL) (``CompMooring = 5``, see
:doc:`openfast`). Deck syntax is in :doc:`driver_format`; the input
files a deck can reference are in :doc:`file_formats`.

Every file is **replaced** on each run (``STATUS='REPLACE'``); CableDyn never appends to a
previous result. All physical quantities are SI: metres, seconds, newtons, radians only where
stated (angles in the output files are always degrees).

.. contents:: On this page
   :local:
   :depth: 2

Standalone driver files
-----------------------

The driver names every file from the ``<out_root>`` argument. Which files appear depends on the
*route* the deck selects, which is decided by the deck content:

.. list-table::
   :header-rows: 1
   :widths: 30 14 14 14 14 14

   * - Deck
     - ``.out``
     - ``.static.out``
     - ``.elements.out``
     - ``.Line<L>.{p,t}.out``
     - ``.Rod<R>.p.out``
   * - Static only (no ``dtM``/``TMax``), all lines ``EI = 0``
     - one row, ``t = 0``
     - yes
     - --
     - static tables
     - --
   * - Dynamic, all lines ``EI = 0``: single lines, and Connect/Free points, Rigid6 bodies
       and rods when they are not on the multibody march (``staggered bodyScheme``, a
       ``motionFile``, a ``FAILURE`` section, ``Coupled``/``Vessel`` rods)
     - time series
     - **no**
     - --
     - time series
     - rod decks only
   * - Dynamic multibody march (bodies, rods or ``Connect``/``Free`` points with their lines,
       ``EI = 0`` or finite-EI; the *Multibody march* of :doc:`driver_format`)
     - time series
     - no
     - --
     - time series
     - yes
   * - Dynamic bodies and rods without lines
     - time series
     - no
     - --
     - --
     - --
   * - Dynamic, all lines finite-EI, every End B on a ``Fixed`` point (cubic-Hermite route)
     - time series
     - yes (the converged initial state)
     - yes
     - time series
     - --
   * - Dynamic, all lines finite-EI, some End B not ``Fixed``
     - time series
     - no
     - --
     - time series
     - --
   * - Mixed ``EI = 0`` and finite-EI lines (static or dynamic)
     - ``t = 0`` row, plus time series when ``TMax > 0``
     - yes
     - --
     - rejected (fail-closed)
     - --

The ``EI = 0`` dynamic routes never write ``<out_root>.static.out``: use a separate static-only
run of the same deck (remove ``dtM`` and ``TMax``) when the static range table is needed. A
line or rod file is written only when the corresponding ``Outputs`` flag is set. The range graph
``<out_root>.Line<L>.range.out`` of the ``Outputs`` flag ``r`` is written on every route of the
table that has lines, the mixed route included.

Common layout
~~~~~~~~~~~~~

All standalone files are plain ASCII tables:

* fields are separated by a single **TAB** character; numeric fields are also right-justified
  inside a fixed width, so the files align in a terminal and parse with any whitespace splitter;
* the **time** column is Fortran ``ES25.16E3`` (25 characters, 17 significant digits, a
  three-digit exponent, e.g. ``  1.2500000000000000E+001``) so that time stamps survive long
  runs with small ``dtM`` without rounding;
* every other real value is Fortran ``ES15.7E3`` (15 characters, 8 significant digits, e.g.
  `` 5.0988032E+005``, three-digit exponent). The one exception is the channel columns of the
  mixed ``EI = 0`` + finite-EI route's ``.out``, written as ``ES15.7`` (the same 8 significant
  digits with a two-digit exponent, e.g. `` 5.0988032E+05``);
* integer key columns (``LineID``, ``Node``, ``Element``, ``Segment``) are written without
  padding (``I0``);
* a non-finite channel value is never written: a dynamic row that would contain ``NaN`` or
  ``Inf`` stops the run with a named error and exit code 2.

``<out_root>.out`` — the main table
-----------------------------------

Written by every route. Layout:

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - Line
     - Content
   * - 1
     - title line beginning with ``#`` (see below)
   * - 2
     - header row: ``Time(s)`` followed by each ``OUTPUTS`` channel token **exactly as written
       in the deck** (original spelling and case)
   * - 3 …
     - one data row per output time

There is **no units row**: units are fixed by the channel vocabulary below (the Python reader
:func:`cabledyn.read_output` attaches them automatically). A deck without an ``OUTPUTS``
section produces a valid file containing only the ``Time(s)`` column.

The title line identifies the route that produced the file:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Route
     - Title line
   * - static
     - ``# CableDyn driver output (static IC; converged=T)`` (``converged=F`` when a line did
       not converge — the file is then written for inspection and the driver exits with 2)
   * - ``EI = 0`` single-line dynamics
     - ``# CableDyn driver output (dynamic; convergence quality is reported by status message)``
   * - Connect/Free-point dynamics
     - ``# CableDyn driver output (dynamic points; convergence quality is reported by status
       message)``
   * - Rigid6 bodies
     - ``# CableDyn driver output (Rigid6 dynamic; convergence quality is reported by status
       message)``
   * - rods
     - ``# CableDyn driver output (rod dynamic; convergence quality is reported by status
       message)``
   * - finite-EI (cubic-Hermite)
     - ``# CableDyn driver output (finite-EI dynamic; production cubic-Hermite route)``
   * - finite-EI, two moving ends
     - ``# CableDyn driver output (finite-EI dynamic; convergence quality in the status
       message)``
   * - mixed ``EI = 0`` + finite-EI
     - ``# CableDyn driver output (mixed aggregate dynamic; convergence quality is reported by
       status message)``
   * - multibody march
     - ``# CableDyn driver output (multibody dynamic; convergence quality is reported by status
       message)``
   * - bodies and rods without lines
     - ``# CableDyn driver output (bodies and rods without lines, dynamic)``

Rows: a static run writes one row at ``t = 0``. A dynamic run writes the initial state at
``t = 0`` and then one row after every committed step, at ``t = k·dtM`` for
``k = 1 … NINT(TMax/dtM)``. The deck keyword ``dtOut`` is accepted for MoorDyn compatibility
but has no effect: the output cadence is always ``dtM``.

If a dynamic step does not converge, the march stops, the file keeps only the rows of committed
steps, and the driver exits with code 2. Such a file is a diagnostic record, not an accepted
response history.

Example (static route)::

   # CableDyn driver output (static IC; converged=T)
   Time(s)	FairTen1	AnchTen1	FairIncl1	AnchIncl1
     0.0000000000000000E+000	 5.0988032E+005	 1.9241379E+005	 6.7242193E+001	-6.8947041E-001

``<out_root>.static.out`` — static range table
----------------------------------------------

The along-arc static configuration, one row per node of every line (a *range graph* of the static
state). It is independent of the ``OUTPUTS`` selection. On a dynamic cubic-Hermite or mixed
run it records the converged initial state that the march starts from.

Layout: title line, header row, **units row**, then data. The title reads
``CableDyn static configuration profile (deformed arc; public node order EndA -> EndB)``
(static route), ``CableDyn standalone static configuration (deformed arc; public node order
End A -> End B)`` (cubic-Hermite route) or ``CableDyn coupled static configuration (deformed
arc; public node order End A -> End B)`` (mixed route and OpenFAST). The columns are identical
on every route:

.. list-table::
   :header-rows: 1
   :widths: 18 10 72

   * - Column
     - Unit
     - Meaning
   * - ``LineID``
     - ``(-)``
     - deck ``LINES`` id (keys the rows when several lines are present)
   * - ``Node``
     - ``(-)``
     - node number, 1 = End A … N = End B
   * - ``ArcLength``
     - ``(m)``
     - cumulative **deformed** chord length from End A
   * - ``X`` ``Y`` ``Z``
     - ``(m)``
     - node position in the global frame (Z up, ``Z = 0`` at the still-water line)
   * - ``Tension``
     - ``(N)``
     - nodal effective tension, the segment tension (same definition as ``Ten<L>N<J>``)
   * - ``Curvature``
     - ``(1/m)``
     - nodal curvature (same definition as ``Curv<L>N<J>``)
   * - ``BendMoment``
     - ``(N.m)``
     - bend moment (same definition as ``BendMom<L>N<J>``); written as 0 by the ``EI = 0``
       static route
   * - ``Declination``
     - ``(deg)``
     - axial-tangent angle from +Z (0 = up, 90 = horizontal, 180 = down)
   * - ``Inclination``
     - ``(deg)``
     - signed inclination below horizontal, ``Declination − 90``
   * - ``Azimuth``
     - ``(deg)``
     - axial-tangent azimuth from +X toward +Y, in [0, 360)

Plot ``Curvature`` or ``Tension`` against ``ArcLength`` for a fatigue or strength check.

.. note::

   Node tensions are segment tensions; ``.elements.out`` holds pointwise field values. On a
   finite-EI line the ``Tension`` column and ``Ten<L>N<J>`` report, at an interior node, the
   length-weighted average of the element-mean axial forces of the two neighbouring elements
   (the tension that enters equilibrium, as the segment tension of an ``EI = 0`` line). The
   ``MinimumAxialResultant`` and ``MaximumAxialResultant`` columns of ``.elements.out`` are
   extrema of the pointwise field :math:`EA\,(|\mathbf r'| - 1)`. With a stiff ``EA`` that
   field oscillates about the element mean wherever an element cannot follow the line -- at a
   nodal seabed-contact kink near touchdown and at a section junction -- and can dip below
   zero there while the line is tensile. Curvature and bend moment are pointwise in both files.

``<out_root>.elements.out`` — Hermite element extrema
-----------------------------------------------------

Written by the cubic-Hermite route together with ``.static.out``. For every element of every
finite-EI line it records the extrema of the **continuous** element fields at the converged
initial state — located by searching the whole cubic field, not only nodes or quadrature points
— so peak curvature between nodes is not missed. Element 1 is at End A; the local coordinate
``xi`` runs from 0 at the End-A side of the element to 1 at its End-B side.

Layout: title line ``CableDyn continuous Hermite-element extrema (public order End A -> End
B)``, header row, units row, data (12 columns):

.. list-table::
   :header-rows: 1
   :widths: 26 10 64

   * - Column
     - Unit
     - Meaning
   * - ``LineID``
     - ``(-)``
     - deck ``LINES`` id
   * - ``Element``
     - ``(-)``
     - element number, 1 = End A
   * - ``ReferenceArcStart``
     - ``(m)``
     - unstretched arc length from End A to the element's End-A side
   * - ``ReferenceArcEnd``
     - ``(m)``
     - unstretched arc length from End A to the element's End-B side
   * - ``PeakXi``
     - ``(-)``
     - local coordinate of the curvature maximum
   * - ``PeakReferenceArc``
     - ``(m)``
     - unstretched arc length of the curvature maximum
   * - ``PeakCurvature``
     - ``(1/m)``
     - maximum curvature over the element
   * - ``BendMomentAtPeak``
     - ``(N.m)``
     - ``EI × PeakCurvature``
   * - ``MinimumAxialResultant``
     - ``(N)``
     - minimum axial force resultant over the element (negative = compression)
   * - ``MinimumXi``
     - ``(-)``
     - local coordinate of that minimum
   * - ``MaximumAxialResultant``
     - ``(N)``
     - maximum axial force resultant over the element
   * - ``MaximumXi``
     - ``(-)``
     - local coordinate of that maximum

Per-line files (``LINES`` ``Outputs`` flag)
-------------------------------------------

Set the ``Outputs`` column of a ``LINES`` row to any combination of the letters below
(case-insensitive, e.g. ``ptr``); ``-`` or an empty field requests nothing. Any other letter is a
parse error. Mixed ``EI = 0`` + finite-EI decks reject the ``p`` and ``t`` flags (use main-file
channels instead) and accept ``r``.

.. list-table::
   :header-rows: 1
   :widths: 10 40 50

   * - Flag
     - File
     - Contents
   * - ``p``
     - ``<out_root>.Line<L>.p.out``
     - node positions, End A → End B
   * - ``t``
     - ``<out_root>.Line<L>.t.out``
     - segment (element) tensions, End A → End B
   * - ``r``
     - ``<out_root>.Line<L>.range.out``
     - range graph: minimum, maximum and mean over the run of the node tension, curvature, bend
       moment, declination and seabed clearance, one row per node (see below)

``<L>`` is the deck line id. Each file starts with a ``#`` comment line (e.g.
``# CableDyn static line node positions; public node order EndA -> EndB``) and a header row;
there is no units row — units are in the column names.

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - File
     - Static run (one row per node / segment)
     - Dynamic run (one row per output time)
   * - ``.p.out``
     - ``Node  X(m)  Y(m)  Z(m)``
     - ``Time(s)  Node1X(m)  Node1Y(m)  Node1Z(m)  Node2X(m) …  Node<N>Z(m)``
   * - ``.t.out``
     - ``Segment  Tension(N)``
     - ``Time(s)  Segment1Tension(N)  Segment2Tension(N) …  Segment<N−1>Tension(N)``

Dynamic per-line rows are written at the same times as ``.out``.

``<out_root>.Line<L>.range.out`` — range graphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The OrcaFlex range graph of line ``L``, accumulated by the solver during the run: at every node,
the minimum, maximum and mean over the output times of the range window. The window holds every
``.out`` row with ``t ≥ RangeStart`` (OPTION ``RangeStart``, default 0, see :doc:`options`), so a
start-up transient can be excluded; a static run has one sample, ``t = 0``. The values are the
node channels of the same run, ``Ten<L>N<J>``, ``Curv<L>N<J>``, ``BendMom<L>N<J>`` and
``L<L>N<J>Dec``, evaluated at every node: the minimum and maximum equal those of the channel time
histories exactly, and the mean is their arithmetic mean. The file is written when the run
completes; a run that stops early leaves none.

Layout: title line, header row, **units row**, then one row per node, End A first. The title reads
``CableDyn range graph (line <L>; <n> samples from t = <t0> s to t = <t1> s; public node order
End A -> End B)``.

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Column
     - Unit
     - Meaning
   * - ``Node``
     - ``(-)``
     - node number, 1 = End A … N = End B
   * - ``ArcLength``
     - ``(m)``
     - cumulative deformed chord length from End A at the first sample (``t = 0``, the static
       initial state; the ``ArcLength`` of ``.static.out``)
   * - ``TensionMin`` ``TensionMax`` ``TensionMean``
     - ``(N)``
     - effective tension, as ``Ten<L>N<J>`` (end nodes: the ``FairTen``/``AnchTen`` end force)
   * - ``CurvatureMin`` ``CurvatureMax`` ``CurvatureMean``
     - ``(1/m)``
     - curvature, as ``Curv<L>N<J>``
   * - ``BendMomentMin`` ``BendMomentMax`` ``BendMomentMean``
     - ``(N.m)``
     - bend moment, as ``BendMom<L>N<J>`` (0 on ``EI = 0`` lines)
   * - ``DeclinationMin`` ``DeclinationMax`` ``DeclinationMean``
     - ``(deg)``
     - declination of the axial tangent, as ``L<L>N<J>Dec``
   * - ``ClearanceMin`` ``ClearanceMax`` ``ClearanceMean``
     - ``(m)``
     - seabed clearance: node ``z`` minus the seabed elevation below the node (flat ``WtrDpth``
       or the ``bathymetryFile`` surface); negative where the node penetrates the penalty seabed.
       Present only when the deck defines a seabed

The accumulation keeps three numbers per node and quantity and makes no heap allocation per step.
Sampling a line evaluates its node channels once per output row; on the 1024-element cubic-Hermite
lazy-wave example this adds about 3 % to the step time. A coupled OpenFAST run does not write range
files and rejects the ``r`` flag; use the ``TDP<L>`` and node channels there.

``<out_root>.modes.out`` — natural frequencies and mode shapes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Written when the deck sets ``nModes`` (see :doc:`options`), before any dynamic march. The file
holds two tab-separated tables, each introduced by a ``#`` comment line, a header row and a units
row:

* ``# Natural frequencies``: ``LineID  Mode  Frequency  Period  Omega`` (Hz, s, rad/s), the
  ``nModes`` lowest modes of every line in ascending order;
* ``# Mode shapes``: ``LineID  Mode  Node  X  Y  Z  dX  dY  dZ``, one row per mode and node.
  ``X Y Z`` is the static node position and ``dX dY dZ`` the nodal translation of the mode,
  scaled so that the largest nodal displacement of the mode is 1. Nodes are numbered End A first,
  as in the other line files; the held end nodes have zero displacement.

Per-rod files (``RODS`` ``Outputs`` flag)
-----------------------------------------

For a dynamic rod deck, set a rod's ``Outputs`` field to ``p`` to write
``<out_root>.Rod<R>.p.out`` (``<R>`` = deck rod id). ``p`` is the only rod flag; any other
letter is a parse error. Layout: ``#`` comment line, header row
``Time(s)  EndAX(m)  EndAY(m)  EndAZ(m)  EndBX(m)  EndBY(m)  EndBZ(m)``, then one row per
output time.

``OUTPUTS`` channel vocabulary
------------------------------

Channel names are listed in the deck ``OUTPUTS`` section. Each row may hold one or several
names separated by whitespace, commas or tabs, optionally quoted with ``"`` or ``'`` (the
OpenFAST ``OutList`` style). Maintained decks use one double-quoted name per row::

   ---------------------- OUTPUTS -----------------------------------------
   "FairTen1"
   "AnchTen1"
   "Ten1N10"
   "L1N10pz"
   -------------------------------------------------------------------------

Matching is **case-insensitive** (``fairten1`` = ``FairTen1``), but the header row repeats the
token as written. A name may be at most 64 characters. Each channel may be requested **once**:
the deck is rejected (exit code ``1``, naming the channel, the earlier spelling and the deck
line) when two names select the same quantity -- the same name in any case or across rows, a
numeric id written with leading zeros (``Ten1N02`` = ``Ten1N2``), or an alias spelling
(``Con<P>p{x,y,z}`` = ``Point<P>p{x,y,z}``, ``FairAngle<L>`` = ``FairDecl<L>``,
``AnchAngle<L>`` = ``AnchDecl<L>``). ``<L>`` and ``<P>`` are deck
``LINES`` and ``POINTS`` ids (not array positions); ``<J>`` is a node number on line ``<L>``,
1 = End A to ``N`` = End B, where ``N`` = 1 + the total ``NumSegs`` of the line's sections.

Line-end channels
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 56 20

   * - Channel
     - Meaning
     - Unit
   * - ``FairTen<L>``
     - line-end tension at End A (fairlead) of line ``L``: on an EI = 0 line, the magnitude of
       the force the line actually exerts on its End-A point -- the end element's tension
       including axial damping, plus the end node's lumped loads (its share of the submerged
       weight, seabed contact and drag at the actual relative velocity), as MoorDyn and
       OrcaFlex report it. It equals the coupled load on that point less the end node's
       inertia, and at rest it is the static end reaction, so the static, ``TMax = 0`` and
       dynamic routes agree at the initial condition. A finite-EI line reports the same
       quantity: the end element's internal end force (axial and bending shear, axial
       damping) plus the end node's share of the submerged weight, buoyancy lost above the
       surface, seabed contact with its damping, and ambient-fluid drag at the actual
       velocity; an end resting on the seabed leaves the end node's weight to the floor. With
       bending it is the end-force magnitude, not the axial effective tension:
       where the end shear is significant it exceeds OrcaFlex's end ``Effective tension`` and
       matches the magnitude of its end ``GX``/``GZ`` force. Seabed friction at a grounded end
       node is part of the ``EI = 0`` end force but not of the finite-EI (cubic-Hermite) end
       force, which leaves the end node's friction, like its inertia, in the coupled load; the
       two agree whenever the end node is off the seabed or friction is off. On the
       two-moving-end finite-EI route (End B not ``Fixed``) the value is instead the signed
       chord-stretch axial tension ``EA (chord/L0 − 1)`` of the end element, without end shear,
       axial damping or the end node's lumped loads
     - N
   * - ``AnchTen<L>``
     - line-end tension at End B (anchor) of line ``L``, defined as for ``FairTen``. On a
       grounded run it is close to the horizontal tension. A ``Fixed`` anchor at the seabed
       depth sits slightly above the grounded chain, which sinks into the penalty seabed, so the
       anchor also carries part of its end node's weight (for example 192.4 kN against a
       horizontal tension of 191.5 kN in :doc:`tutorial_catenary`)
     - N
   * - ``FairIncl<L>`` / ``AnchIncl<L>``
     - signed inclination of the End-A / End-B tangent below horizontal
       (``Decl − 90``: 0 = horizontal, +90 = pointing down, −90 = pointing up)
     - deg
   * - ``FairDecl<L>`` / ``AnchDecl<L>``
     - declination of the End-A / End-B tangent from +Z (0 = up, 90 = horizontal, 180 = down)
     - deg
   * - ``FairAngle<L>`` / ``AnchAngle<L>``
     - aliases of ``FairDecl<L>`` / ``AnchDecl<L>``
     - deg

A grounded anchor segment reads ``AnchIncl ≈ 0``; the small non-zero value reports the actual
orientation of the last element, which CableDyn does not force to horizontal.

Line-node channels
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 56 20

   * - Channel
     - Meaning
     - Unit
   * - ``Ten<L>N<J>``
     - effective tension at node ``J``: an interior node takes the mean of its two neighbouring
       element (segment) tensions -- on a finite-EI line their element-mean axial forces,
       weighted by element length; an end node reports the line-end force of ``FairTen`` /
       ``AnchTen``. On the two-moving-end finite-EI route every value is the chord-stretch
       tension of ``FairTen`` there: the adjacent element at an end node and the unweighted
       mean of the two neighbouring elements at an interior node
     - N
   * - ``Curv<L>N<J>``
     - curvature at node ``J``. ``EI = 0`` lines: discrete curvature of the circle through
       the node and its two neighbours (end nodes take the adjacent interior value).
       Cubic-Hermite lines: exact curvature of the continuous centreline, taking the larger of
       the two one-sided element values at an interior node
     - 1/m
   * - ``BendMom<L>N<J>``
     - bend moment ``EI × κ`` measured relative to the stress-free reference shape. The
       cubic-Hermite reference is straight, so there it equals ``EI × Curv<L>N<J>``
       (larger one-sided value at an interior node); 0 on ``EI = 0`` lines
     - N·m
   * - ``L<L>N<J>px`` / ``py`` / ``pz``
     - node position component
     - m
   * - ``L<L>N<J>vx`` / ``vy`` / ``vz``
     - node velocity component (0 on a static run)
     - m/s
   * - ``L<L>N<J>ax`` / ``ay`` / ``az``
     - node acceleration component (0 on a static run)
     - m/s²
   * - ``L<L>N<J>Dec``
     - declination of the node's axial tangent from +Z (0 = up, 90 = horizontal, 180 = down)
     - deg
   * - ``L<L>N<J>Azi``
     - azimuth of the node's axial tangent from +X toward +Y, in [0, 360)
     - deg

Touchdown channels
~~~~~~~~~~~~~~~~~~

For a line that rests on the seabed at one end, the touchdown point (TDP) at every output time.
A node is grounded when its centreline is at most 1e-6 m (the height at which the seabed contact
law engages) above the seabed; the grounded end is the end grounded in the initial state. Walking
from that end, the TDP lies between the last grounded node and the next one, where the centreline
crosses that height (linear interpolation of the clearance). The definitions are those of
:func:`cabledyn.touchdown_history` with ``tolerance = 1e-6``.

.. list-table::
   :header-rows: 1
   :widths: 24 56 20

   * - Channel
     - Meaning
     - Unit
   * - ``TDP<L>s``
     - arc length of the TDP from End A (deformed chord length along the nodes)
     - m
   * - ``TDP<L>x`` / ``y`` / ``z``
     - TDP position
     - m
   * - ``TDP<L>Lay``
     - layback: horizontal distance from the TDP to the suspended end
     - m
   * - ``TDP<L>Exc``
     - TDP excursion: horizontal displacement of the TDP from its initial (``t = 0``) position,
       projected on the initial horizontal direction from the TDP to the suspended end (positive
       toward the suspended end)
     - m

The channels are available on every route, OpenFAST included, when the deck has a seabed
(``WtrDpth`` or ``bathymetryFile``). A line whose initial state is grounded at both ends or at
neither stops the run with exit code 1 naming the line. When the grounded end later lifts off,
the channels report that end node; when the whole line rests on the seabed, the suspended end
node.

Point channels
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 56 20

   * - Channel
     - Meaning
     - Unit
   * - ``Point<P>px`` / ``py`` / ``pz``
     - position component of point ``P`` (for a prescribed point on a ``motionFile`` run, the
       prescribed position at that time)
     - m
   * - ``Con<P>px`` / ``py`` / ``pz``
     - MoorDyn v1 spelling, identical to ``Point<P>p…``
     - m
   * - ``Point<P>Fx`` / ``Fy`` / ``Fz`` / ``FH``
     - resultant of the forces of the lines and finite-EI cables attached to point ``P`` (each
       the ``FairTen``/``AnchTen`` end force): the load on a shared anchor, or the line load a
       free point balances; ``FH`` is the horizontal magnitude. Available on every route
       (static, ``EI = 0``, finite-EI, multibody and OpenFAST) except the two-moving-end
       finite-EI route, which rejects it by name
     - N

Body and rod channels
~~~~~~~~~~~~~~~~~~~~~

The MoorDyn-F names, for Rigid6 bodies and for rods (free, fixed, pinned, prescribed, fixed or
pinned to a body), on every route that carries them, the OpenFAST ``CompMooring = 5`` route
included. Loads are evaluated at the committed state of the output time.

.. list-table::
   :header-rows: 1
   :widths: 30 52 18

   * - Channel
     - Meaning
     - Unit
   * - ``Body<N>Px`` / ``Py`` / ``Pz``
     - reference-point position
     - m
   * - ``Body<N>Rx`` / ``Ry`` / ``Rz``
     - attitude, x-y'-z'' Euler angles of the deck convention
     - deg
   * - ``Body<N>Vx`` … ``Vz``, ``RVx`` … ``RVz``
     - reference-point velocity; angular velocity
     - m/s, deg/s
   * - ``Body<N>Ax`` … ``Az``, ``RAx`` … ``RAz``
     - reference-point acceleration; angular acceleration
     - m/s², deg/s²
   * - ``Body<N>Fx`` / ``Fy`` / ``Fz``, ``Mx`` / ``My`` / ``Mz``
     - net external force and moment on the body about its reference point (global axes):
       weight, buoyancy and hydrostatic restoring, Morison and external loads, seabed contact,
       the loads of its fixed rods, the attached line and cable end forces, and the pin forces
       of the rods pinned to it. Zero for a free body at rest in equilibrium
     - N, N·m
   * - ``Rod<N>Px`` / ``Py`` / ``Pz``, ``Vx`` … ``Vz``, ``Ax`` … ``Az``
     - End A position, velocity and acceleration
     - m, m/s, m/s²
   * - ``Rod<N>Rx`` / ``Ry``
     - MoorDyn's roll and pitch of the rod axis: its tilt φ from the vertical times
       ``-sin β`` and ``cos β``, β the heading of the axis
     - deg
   * - ``Rod<N>RVx`` … ``RVz``, ``RAx`` … ``RAz``
     - angular velocity and acceleration
     - deg/s, deg/s²
   * - ``Rod<N>Fx`` / ``Fy`` / ``Fz``, ``Mx`` / ``My`` / ``Mz``
     - net external force on the rod and its moment about End A: weight, buoyancy, Morison and
       seabed loads and the line and cable end forces (not the pin reaction of a pinned rod,
       whose ``M`` therefore vanishes at a static equilibrium)
     - N, N·m
   * - ``Rod<N>TenA`` / ``TenB``
     - magnitude of the summed line and cable end force at End A / End B
     - N
   * - ``Rod<N>Sub``
     - submerged fraction of the rod length, below the local waterline (still water: z = 0)
     - –
   * - ``Rod<N>N<k>Px`` / ``Py`` / ``Pz``
     - position of rod node ``k`` (0 = End A, ``NumSegs`` = End B)
     - m

A zero-length rod (``NumSegs`` 0) is modelled as a point; its ``Rod<N>`` channels are rejected
by name, and its motion is reported by the ``Point<P>`` channels of that point.

The axial tangent of every orientation channel points from End A toward End B (OrcaFlex's node
``Ez`` axis). Curvature, declination and azimuth are evaluated from the solved geometry on
every route (static, ``EI = 0`` dynamic, finite-EI dynamic, rod and Rigid6 decks).

Validation of channel names
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every name is checked while the deck is parsed, before any solve, and a bad name stops the run
with exit code 1. A name is rejected when:

* it matches none of the forms above, or carries trailing text (``Point2px_raw``,
  ``Point2pzz``, ``FairTen1x``);
* it references an unknown line or point id, or a node number larger than the line's node
  count, or is a ``TDP<L>`` name with a suffix other than ``s``, ``x``, ``y``, ``z``, ``Lay`` or
  ``Exc``;
* on a mixed ``EI = 0`` + finite-EI deck (and in OpenFAST), it is a ``Point<P>`` channel of a
  ``Coupled``/``Vessel`` point attached **only** to finite-EI lines — that point is not part of
  the ``EI = 0`` point system that serves point channels. Use the cable's
  ``L<L>N<J>p{x,y,z}`` channel instead. ``Fixed`` points remain valid because they never move;
* in an OpenFAST ``CompMooring = 5`` run, it is longer than OpenFAST's 20-character channel
  header (``ChanLen``). The coupled initialisation stops with an error naming the channel
  rather than writing a truncated, possibly duplicate header.

Console output
--------------

The standalone driver writes a short initialisation report after the static solve. On the
``EI = 0``, finite-EI and static routes it goes to **stderr** (after the identity banner)::

     Parsing CableDyn input file: examples/chain_catenary_r3_100m.dat
      Created CableDyn model: 1 line object(s), 2 point(s), 1 section(s) [EI=0: 1, finite-EI: 0].
      Initial conditions: Newton static equilibrium with load continuation completed.
      Fairlead convention: force is on End A toward End B; inclinations are signed below horizontal.
      Line 1 fairlead effective tension:  5.09880E+005 N
         force [Fx, Fy, Fz]: [ 1.91549E+005,  0.00000E+000, -4.72532E+005] N, inclination=   67.934 deg
         line tangent: inclination=   67.242 deg, declination=  157.242 deg, azimuth=    0.000 deg
     CableDyn initialization completed.

``force`` is the force on the End-A node directed toward End B, with its inclination below the
horizontal. ``line tangent`` gives the direction of the line itself at End A. The two angles
differ slightly: the end force also carries the end node's share of the distributed load (weight,
drag, seabed reaction) and, on a finite-EI line, the end shear. On a finite-EI deck whose End B
is not ``Fixed`` the line is labelled ``axial force component`` instead: that route reports the
axial part of the end resultant only.
On a mixed ``EI = 0`` + finite-EI deck the report is shorter and goes to **stdout**::

     CableDyn mixed standalone aggregate: <n> line(s) [<n0> EI=0, <n1> finite-EI].
     Parsed <np> point(s) and <ns> section row(s).
     Static equilibrium fairlead results:
       Line <L>: FairTen=<T> N, tangent inclination=<deg> deg, force=(<Fx> <Fy> <Fz> ) N

Dynamic runs then print progress to **stdout** at 5 % intervals of the march, with the ETA
estimated from the average wall time per committed step::

     Dynamic simulation: 2000 step(s), simulated duration      200.000 s, dtM =  1.00000E-01 s.
     Progress:    5.0% | t =       10.000 s | elapsed 000:00:02 | ETA 000:00:38

Two further stdout records may follow a completed march: ``Recovery audit: …`` when the
integrator subdivided nominal steps (with the number of subdivided intervals and the maximum
sub-step count used), and, on cubic-Hermite decks with the tensile monitor in warning mode,
``Tensile monitor: line <L> …`` summarising compression events (``tensile_safety warn`` in
``OPTIONS``, see :doc:`options`). The final stdout line of a
successful run is described in :doc:`cli`.

In OpenFAST
-----------

Under ``CompMooring = 5`` the selected channels flow through OpenFAST's normal output system
(``WriteOutput``), and CableDyn additionally writes its own files, named from the OpenFAST
output root:

.. list-table::
   :header-rows: 1
   :widths: 36 64

   * - File
     - Contents
   * - ``<RootName>.CD.static.out``
     - the static range table of the converged initialisation (same columns as
       ``.static.out`` above; title ``CableDyn coupled static configuration …``), at the
       equilibrium the console summary reports, before any SeaState kinematics act. Always
       written, whatever the ``OUTPUTS`` selection.
   * - ``<RootName>.CD.out``
     - time history of the deck ``OUTPUTS`` channels, written only when at least one channel is
       selected
   * - ``<FarmRoot>.FarmCD.static.out`` / ``<FarmRoot>.FarmCD.out``
     - the same two files for a FAST.Farm run, named from the FAST.Farm output root
   * - ``<RootName>.CD.rst.dat`` / ``<RootName>.CD.lin.dat``
     - temporary copies of the CableDyn deck, written in the ``MooringFile`` folder when a run
       restarts from a checkpoint or linearises; the module rebuilds its model from them and
       removes them after use (see :doc:`openfast`)

``<RootName>.CD.out`` layout — note that it differs from the standalone ``.out``:

* **no title line**;
* header row ``Time`` followed by the channel tokens;
* a **units row**: ``(s)`` followed by each channel's unit — ``(N)`` tensions and forces,
  ``(deg)`` angles, ``(m)`` positions, ``(m/s)`` velocities, ``(m/s2)`` accelerations,
  ``(deg/s)`` / ``(deg/s2)`` body and rod angular rates, ``(1/m)`` curvature, ``(N.m)`` bend
  moment and body/rod moments, ``(-)`` ``Rod<N>Sub``;
* data rows: time as ``ES25.16E3`` (full double precision), values as ``ES15.6E2``, TAB-separated.

The first data row is the static initialisation at ``t = 0``; every further row is written after
each committed CableDyn step, at every ``dtM`` boundary. It is independent of OpenFAST's ``DT_Out``
and never repeats held values between
CableDyn steps. If OpenFAST corrects a step (predictor–corrector iterations), the provisional row
is replaced rather than duplicated; the file is flushed after every row, so it is readable during
a run and intact after a crash.

Select CableDyn channels in the ``OUTPUTS`` section of the CableDyn ``MooringFile``; they are not
members of the top-level ``.fst`` ``OutList``. The ``TDP<L>`` channels are evaluated from the
coupled state against the touchdown reference of the converged initialisation; the range files of
the ``LINES`` flag ``r`` are a standalone output, and a coupled deck with that flag is rejected. For
a single-turbine run the module also prints
every selected channel's static value to the screen after initialisation. A deck with no
``OUTPUTS`` section is valid but publishes no channels.

Reading and post-processing
---------------------------

:func:`cabledyn.read_output` reads every file on this page strictly and attaches the units
defined above, returning time-history or static-profile objects that plot ``FairTen1`` against
time, ``Curvature`` or ``Tension`` against ``ArcLength``, or the static centreline.
:func:`cabledyn.read_range_graphs` returns the envelopes of a ``.range.out`` file as
:class:`cabledyn.RangeGraph` objects. Dynamic
``Line<L>.p.out`` and ``Line<L>.t.out`` files additionally expose interpolated snapshots and
per-segment envelopes. See :doc:`python` and :doc:`api_python`.

For time-history channels, :meth:`cabledyn.TimeHistory.fatigue` provides weighted rainflow
cycles and an uncorrected damage-equivalent **range** (the Wöhler exponent and reference cycle
count or frequency are mandatory inputs); :meth:`cabledyn.TimeHistory.spectrum` a one-sided
Welch power spectral density with spectral moments; and :meth:`cabledyn.TimeHistory.coherence`
the magnitude-squared coherence of two channels. These are derived Python results; they never
modify the native files.

Some pyDatView versions do not recognise the multi-line ``.static.out`` table directly. Export
one line to a normalised table with units embedded in the header:

.. code-block:: powershell

   cabledyn-post export case.static.out line4_static.csv --line 4

The conversion preserves values and row order, does not add a time column, and refuses to
overwrite an existing file unless ``--overwrite`` is given.
