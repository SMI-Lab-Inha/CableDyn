.. SPDX-License-Identifier: Apache-2.0

Python API reference
====================

.. module:: cabledyn

This page documents every public name of the :mod:`cabledyn` package. Everything is imported
from the top-level package, for example ``from cabledyn import CableDynDriver, read_output``.
The workflows that use these objects are described in :doc:`python`.

Two groups of names have different installation requirements:

* The standalone-driver API (running decks, reading and analysing results, editing decks, batch
  studies, and all exception classes) is pure Python. It needs NumPy and, to run analyses, the
  ``CableDyn_driver`` executable, but never the CableDyn shared library.
* The in-process coupling API (:class:`CableDyn`, :func:`abi_version`, :func:`abi_minor`,
  :func:`library_path`, and :func:`version_string`) loads the CableDyn shared library the first
  time one of these names is accessed. These names are deliberately left out of
  ``cabledyn.__all__``, so
  ``from cabledyn import *`` works in an installation without the shared library.

Plotting methods import matplotlib and ``to_dataframe`` imports pandas only when they are
called; both are optional dependencies. Every plotting method returns the matplotlib axes it
drew on and draws on a new figure unless an ``ax`` is supplied. The exception is
:func:`animate`, which returns a matplotlib ``FuncAnimation``.

Unless stated otherwise, quantities are in SI units: metres, seconds, kilograms, and newtons.
Result objects are immutable, and the NumPy arrays they hold are read-only.

Running the standalone driver
-----------------------------

.. autoclass:: cabledyn.CableDynDriver
   :members:

.. autoclass:: cabledyn.DriverResult
   :members:

Batch studies
-------------

.. autofunction:: cabledyn.run_study

.. autoclass:: cabledyn.StudyResult
   :members:

.. autoclass:: cabledyn.StudyCaseResult
   :members:

Parameter sweeps
----------------

:func:`parameter_grid` builds the case specifications that :func:`generate_deck_cases`
consumes, from lists of values for each deck selector.

.. autofunction:: cabledyn.parameter_grid

Reading results
---------------

:func:`read_output` returns the most specific table class for the file it reads:
:class:`TimeHistory` for the main output, :class:`LineNodeHistory` and
:class:`LineSegmentHistory` for per-line dynamic outputs, :class:`StaticProfile` for the static
profile, and :class:`OutputTable` for any other numeric table. :class:`TimeHistory` and
:class:`StaticProfile` are subclasses of :class:`OutputTable`; :class:`LineNodeHistory` and
:class:`LineSegmentHistory` are subclasses of :class:`TimeHistory`, so every table method is
available on them.

.. autofunction:: cabledyn.read_output

.. autoclass:: cabledyn.OutputTable
   :members:

.. autoclass:: cabledyn.TimeHistory
   :members:

.. autoclass:: cabledyn.LineNodeHistory
   :members:

.. autoclass:: cabledyn.LineSegmentHistory
   :members:

.. autoclass:: cabledyn.StaticProfile
   :members:

.. autoclass:: cabledyn.ChannelStatistics
   :members:

.. autoclass:: cabledyn.StaticLineSummary
   :members:

.. autoclass:: cabledyn.SpatialStatistics
   :members:

Reading other codes' outputs
----------------------------

These readers return the same table classes as :func:`read_output`, so statistics, fatigue,
spectra, plots, comparisons, and exports work unchanged on results of OpenFAST (maintained by
NLR, the National Laboratory of the Rockies, formerly NREL) and MoorDyn and on the files of an
OpenFAST run that uses CableDyn as its mooring module. :func:`read_table`
chooses the reader from the file name. Malformed files raise :exc:`OutputFormatError`.

.. autofunction:: cabledyn.read_table

.. autofunction:: cabledyn.read_openfast_output

.. autofunction:: cabledyn.read_moordyn_output

.. autofunction:: cabledyn.read_moordyn_line

.. autoclass:: cabledyn.MoorDynLineHistory
   :members:

.. autofunction:: cabledyn.read_coupled_run

.. autoclass:: cabledyn.CoupledRun
   :members:

Comparing runs
--------------

:func:`compare_histories` aligns two :class:`TimeHistory` tables on one time grid, without
extrapolating beyond the interval they share, and reports error metrics and statistic changes
for each channel. Differences are candidate minus reference, in the unit of the channel.

.. autofunction:: cabledyn.compare_histories

.. autoclass:: cabledyn.HistoryComparison
   :members:

.. autoclass:: cabledyn.ChannelComparison
   :members:

Signal processing
-----------------

Each function returns a new table of the same class as its input, with the same channels,
units, source path, and title; the input is never modified. The filters require uniform
sampling, so resample a variable-step record first.

.. autofunction:: cabledyn.resample

.. autofunction:: cabledyn.moving_average

.. autofunction:: cabledyn.fft_filter

Fatigue analysis
----------------

:meth:`TimeHistory.fatigue` is the usual entry point. The functions below apply the same
rainflow counting and damage-equivalent-range calculation to any scalar sequence.

.. autofunction:: cabledyn.rainflow_cycles

.. autofunction:: cabledyn.cycle_histogram

.. autofunction:: cabledyn.damage_equivalent_range

.. autoclass:: cabledyn.FatigueResult
   :members:

.. autoclass:: cabledyn.RainflowCycle
   :members:

.. autoclass:: cabledyn.RainflowHistogram
   :members:

Fatigue damage
--------------

A :class:`FatigueCurve` gives cycles to failure for a stress range (S--N, MPa) or a normalised
tension range (T--N). The factory functions return curves with published constants; the source is
in each docstring and in :attr:`FatigueCurve.source`. :func:`channel_damage` and
:func:`damage_along_arc` apply Palmgren--Miner damage to rainflow cycles, and
:func:`lifetime_fatigue` weights the damage of several sea states over a design life.

.. autoclass:: cabledyn.FatigueCurve
   :members:

.. autofunction:: cabledyn.dnv_rp_c203_curve

.. autofunction:: cabledyn.dnv_os_e301_curve

.. autofunction:: cabledyn.api_rp_2sk_curve

.. autofunction:: cabledyn.chain_nominal_area

.. autofunction:: cabledyn.miner_damage

.. autofunction:: cabledyn.channel_damage

.. autoclass:: cabledyn.ChannelDamage
   :members:

.. autofunction:: cabledyn.cable_stress

.. autofunction:: cabledyn.damage_along_arc

.. autoclass:: cabledyn.DamageProfile
   :members:

.. autofunction:: cabledyn.study_sea_state_damage

.. autoclass:: cabledyn.SeaStateDamage
   :members:

.. autofunction:: cabledyn.lifetime_fatigue

.. autoclass:: cabledyn.LifetimeFatigue
   :members:

Range graphs
------------

.. autoclass:: cabledyn.RangeGraph
   :members:

.. autofunction:: cabledyn.read_range_graphs

.. autofunction:: cabledyn.read_range_graph

.. autofunction:: cabledyn.node_range_graph

.. autofunction:: cabledyn.static_range_graph

.. autofunction:: cabledyn.element_range_graph

Design checks
-------------

.. autofunction:: cabledyn.bend_check

.. autoclass:: cabledyn.BendLimits
   :members:

.. autoclass:: cabledyn.BendCheck
   :members:

.. autofunction:: cabledyn.mbr_utilisation

.. autofunction:: cabledyn.tension_check

.. py:data:: DNV_OS_E301_STRENGTH_FACTOR
   :type: float
   :value: 0.95

   DNV-OS-E301 characteristic strength of a new chain or steel-wire-rope line body,
   ``S_C = 0.95 S_mbs``: the default ``strength_factor`` of :func:`cabledyn.tension_check`.
   Other line bodies take the factor of their design basis.

.. autoclass:: cabledyn.TensionCheck
   :members:

.. py:data:: API_RP_2SK_SAFETY_FACTORS
   :type: dict[tuple[str, str], float]

   API RP 2SK (3rd edition, 2005) factors of safety on the maximum line tension, keyed by
   (condition, analysis): intact 2.00 (quasi-static) and 1.67 (dynamic), damaged (one line
   broken) 1.43 and 1.25, transient 1.18 and 1.05.

.. py:data:: DNV_OS_E301_PARTIAL_FACTORS
   :type: dict[tuple[str, int], tuple[float, float]]

   DNV-OS-E301 (Ch.2 Sec.2) partial safety factors (gamma_mean, gamma_dyn) on the mean and
   dynamic tension, keyed by (limit state, consequence class): ULS class 1 (1.10, 1.50), ULS
   class 2 (1.40, 2.10), ALS class 1 (1.00, 1.10), ALS class 2 (1.00, 1.25).

Spectral analysis
-----------------

:meth:`TimeHistory.spectrum` and :meth:`TimeHistory.coherence` are the usual entry points. The
functions below apply the same Welch estimators to arrays held in memory.

.. autofunction:: cabledyn.power_spectrum

.. autofunction:: cabledyn.magnitude_squared_coherence

.. autoclass:: cabledyn.PowerSpectrum
   :members:

.. autoclass:: cabledyn.SpectralPeak
   :members:

.. autoclass:: cabledyn.CoherenceResult
   :members:

Extreme values
--------------

:func:`block_maxima` and :func:`upcrossing_maxima` sample extremes from a record;
:func:`fit_gumbel` and :func:`fit_weibull` fit a distribution to them. The fits describe the
sample they are given: they do not check independence or stationarity, do not extrapolate a
block length, and give no confidence interval.

.. autofunction:: cabledyn.block_maxima

.. autofunction:: cabledyn.upcrossing_maxima

.. autofunction:: cabledyn.fit_gumbel

.. autoclass:: cabledyn.GumbelFit
   :members:

.. autofunction:: cabledyn.fit_weibull

.. autoclass:: cabledyn.WeibullFit
   :members:

Line geometry
-------------

:func:`line_geometry` derives arc length, inclination, and a discrete curvature estimate from
the node coordinates of one line. These are post-processing estimates from the output nodes;
the solver's own curvature, tension, and angle channels, where written, are authoritative.

.. autofunction:: cabledyn.line_geometry

.. autoclass:: cabledyn.LineGeometry
   :members:

.. autoclass:: cabledyn.Touchdown
   :members:

.. autofunction:: cabledyn.touchdown_history

.. autoclass:: cabledyn.TouchdownHistory
   :members:

Deck editing and case generation
--------------------------------

:class:`DeckFile` edits an existing deck while preserving everything it does not change.
:class:`DeckWriter` builds a new deck from Python values. :func:`generate_deck_cases` writes a
family of edited decks for a batch study. It and :meth:`DeckModel.save` (with ``rebase=True``)
rewrite relative ancillary paths for the new folder and copy the fixed-name MoorDyn-C
kinematics files (``wave_elevation.txt`` for ``3 WaveKin``, ``wave_frequencies.txt`` for
``7 WaveKin``, ``current_profile.txt`` for ``1 Currents``), which the solver reads from the deck's
own folder, next to the new decks (:doc:`python`).

.. autoclass:: cabledyn.DeckFile
   :members:

.. autoclass:: cabledyn.deck_file.DeckRecord
   :members:

.. autoclass:: cabledyn.deck_file.OptionRecord
   :members:

.. autofunction:: cabledyn.generate_deck_cases

.. autoclass:: cabledyn.GeneratedCase
   :members:

.. autoclass:: cabledyn.DeckWriter(title="cabledyn deck")
   :members:

Deck object model
-----------------

.. autoclass:: cabledyn.DeckModel
   :members:

.. autoexception:: cabledyn.DeckReferenceError

.. autoclass:: cabledyn.SyropeEA
   :members:

The row and object classes (``LineType``, ``Line``, ``Section``, ``Point``, ``Body``, ``Rod``,
``Failure``, ``EndConnection``, ...) are dataclasses in ``cabledyn.builder``.

Snapshots and animation
-----------------------

.. autoclass:: cabledyn.Snapshots
   :members:

.. autoclass:: cabledyn.Recorder
   :members:

.. autofunction:: cabledyn.animation.record

.. autofunction:: cabledyn.animate

.. autoclass:: cabledyn.Seabed
   :members:

.. autoclass:: cabledyn.WaterSurface
   :members:

Results at a time
-----------------

.. autofunction:: cabledyn.line_positions

.. autoclass:: cabledyn.LinePositions
   :members:

.. autofunction:: cabledyn.line_field

.. autoclass:: cabledyn.LineField
   :members:

.. autofunction:: cabledyn.available_quantities

.. autofunction:: cabledyn.profile_at

.. autofunction:: cabledyn.profiles_at

.. autoclass:: cabledyn.ArcProfile
   :members:

.. autofunction:: cabledyn.line_range_graph

Clearance
---------

.. autofunction:: cabledyn.seabed_clearance

.. autoclass:: cabledyn.SeabedClearance
   :members:

.. autofunction:: cabledyn.line_clearance

.. autoclass:: cabledyn.LineClearance
   :members:

.. autofunction:: cabledyn.clearance_matrix

.. autoclass:: cabledyn.ClearanceMatrix
   :members:

.. autofunction:: cabledyn.segment_distance

.. autofunction:: cabledyn.read_bathymetry

.. autoclass:: cabledyn.Bathymetry
   :members:

Summary tables
--------------

.. autofunction:: cabledyn.channel_summary

.. autofunction:: cabledyn.line_summary

.. autoclass:: cabledyn.ChannelSummary
   :members:

.. autoclass:: cabledyn.LineSummary
   :members:

.. autoclass:: cabledyn.SummaryTable
   :members:

In-process coupling
-------------------

These names load the CableDyn shared library on first access; see :doc:`python` for how the
library is located. Array ordering, frames, and units at the coupling boundary follow
:doc:`coupling_boundary`, and the underlying C functions are described in :doc:`capi`.

.. py:class:: CableDyn(deck=None)

   One CableDyn solver instance behind the C ABI. The class is a context manager: leaving a
   ``with`` block calls :meth:`close`.

   :param deck: Path of a sectioned ``.dat`` deck; see :doc:`driver_format`. When given, the
      deck is parsed and its static initial condition is solved immediately, as by
      :meth:`init_deck`. ``None`` (the default) creates an uninitialised instance for a later
      :meth:`init_deck` call.
   :type deck: str | os.PathLike | None
   :raises CableDynError: if the handle cannot be created or the deck fails to initialise.
      A failed initialisation releases the handle before the exception propagates.
   :raises ValueError: if the deck path contains a NUL character or, on Windows, cannot be
      represented in the active ANSI code page.

   The coupled exchange is kinematics in, loads out. Kinematics are the position (m), velocity
   (m/s), and acceleration (m/s\ :sup:`2`) of every coupled point, three degrees of freedom per
   point. Loads are the forces (N) the lines exert on those degrees of freedom. Input arrays may
   be flat, shaped ``(n_coupled_dof,)`` with interleaved ``x, y, z`` values per point, or
   column-per-point, shaped ``(3, n_points)``. Returned arrays are always flat ``float64``
   arrays. A ``(3, 3)`` input is read column-per-point, so pass three points flat when in doubt.

   Calls on one instance are serialised by an internal lock. Independent instances may be used
   from different threads; their deck initialisations are serialised inside the library, their
   steps run concurrently. ``dt`` and ``fluid_density`` must be finite real numbers (``dt``
   positive, ``fluid_density`` non-negative) and arrays real-valued; other values raise
   ``TypeError`` or ``ValueError`` before the library is called.

   .. py:method:: close()

      Release the native handle. Calling it again is harmless.

   .. py:method:: last_error()

      Return the library's most recent diagnostic message for this instance.

      :returns: The message, or an empty string when there is none.
      :rtype: str

   .. py:method:: init_deck(deck)

      Initialise from a deck file: parse it, build the model, and solve the static initial
      condition.

      :param deck: Path of the ``.dat`` deck. Relative ancillary files named in the deck are
         resolved as described in :doc:`driver_format`.
      :type deck: str | os.PathLike
      :raises CableDynError: if parsing, model construction, or the static solve fails. A
         failed initialisation keeps any previously initialised model.
      :raises ValueError: if the path contains a NUL character or, on Windows, cannot be
         represented in the active ANSI code page.

   .. py:property:: initialized
      :type: bool

      Whether the instance holds an initialised model.

   .. py:property:: n_coupled_dof
      :type: int

      Number of coupled degrees of freedom, three per host-driven point.

      :raises CableDynError: if the query fails, for example on an uninitialised model.

   .. py:property:: n_points
      :type: int

      Number of stored system points. This is the number of columns expected by
      :meth:`update_point_fluid_fields`.

      :raises CableDynError: if the query fails, for example on an uninitialised model.

   .. py:property:: n_lines
      :type: int

      Number of lines in the model.

      :raises CableDynError: if the query fails, for example on an uninitialised model.

   .. py:method:: get_coupled_motion()

      Return the current kinematics of the coupled points.

      :returns: ``(q, v, a)``: positions (m), velocities (m/s), and accelerations
         (m/s\ :sup:`2`), each a flat array of shape ``(n_coupled_dof,)``.
      :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
      :raises CableDynError: if the call fails, for example on an uninitialised model.

   .. py:method:: update_states(q, v, a)

      Transfer coupled-point kinematics to the model without advancing time. Use it when a
      host needs loads for trial kinematics before committing a step.

      :param q: Coupled-point positions (m).
      :type q: array_like
      :param v: Coupled-point velocities (m/s).
      :type v: array_like
      :param a: Coupled-point accelerations (m/s\ :sup:`2`).
      :type a: array_like
      :raises ValueError: if an array is neither flat ``(n_coupled_dof,)`` nor
         ``(3, n_coupled_dof // 3)``.
      :raises CableDynError: if the library rejects the call.

   .. py:method:: step(dt, q, v, a)

      Advance one implicit time step of length ``dt`` with the coupled-point kinematics
      prescribed at the end of the step, ``t + dt``.

      :param dt: Step length (s).
      :type dt: float
      :param q: Coupled-point positions at ``t + dt`` (m).
      :type q: array_like
      :param v: Coupled-point velocities at ``t + dt`` (m/s).
      :type v: array_like
      :param a: Coupled-point accelerations at ``t + dt`` (m/s\ :sup:`2`).
      :type a: array_like
      :returns: The number of Newton iterations used.
      :rtype: int
      :raises ValueError: if an array has the wrong shape.
      :raises CableDynError: if the step fails outright. The model is left at time ``t``.
      :raises ConvergenceError: if the step ends without meeting the Newton tolerance. The model
         has then already advanced to ``t + dt`` with its best iterate, so retrying the same
         step would advance time twice.

   .. py:method:: calc_output()

      Return the coupled reaction loads at the current committed state.

      :returns: Forces (N) the lines exert on the coupled degrees of freedom, a flat array of
         shape ``(n_coupled_dof,)``.
      :rtype: numpy.ndarray
      :raises CableDynError: if the call fails.

   .. py:method:: update_point_fluid_fields(fluid_velocity, fluid_acceleration, waterline_z, fluid_density)

      Prescribe external fluid kinematics at every system point.
      This is the interface for coupling to an external flow solver. Call it before
      :meth:`step`.

      :param fluid_velocity: Fluid velocity at each point (m/s), shaped ``(3, n_points)`` or
         flat ``(3 * n_points,)``.
      :type fluid_velocity: array_like
      :param fluid_acceleration: Fluid acceleration at each point (m/s\ :sup:`2`), with the
         same shape as ``fluid_velocity``.
      :type fluid_acceleration: array_like
      :param waterline_z: Global z coordinate of the free surface (m) at each point, shape
         ``(n_points,)``.
      :type waterline_z: array_like
      :param fluid_density: Fluid density (kg/m\ :sup:`3`).
      :type fluid_density: float
      :raises ValueError: if an array has the wrong shape.
      :raises CableDynError: if the library rejects the fields.

   .. py:property:: lines
      :type: tuple[cabledyn.objects.Line, ...]

      The model's lines in ascending deck id.

   .. py:property:: points
      :type: tuple[cabledyn.objects.Point, ...]

      The model's points in solver order.

   .. py:property:: bodies
      :type: tuple[cabledyn.objects.Body, ...]

      The model's Rigid6 bodies.

   .. py:property:: rods
      :type: tuple[cabledyn.objects.Rod, ...]

      The model's rods.

   .. py:method:: line(line_id)
                  point(point_id)
                  body(body_id)
                  rod(rod_id)

      Return the object with the given deck id.

      :raises KeyError: if there is no such object.

   .. py:method:: channel(token)

      Evaluate one output channel at the committed state and return its value. ``token`` uses
      the deck ``OUTPUTS`` vocabulary (:doc:`outputs`), for example ``FairTen1``, ``L2N5px``,
      ``Curv1N3``, ``Point4Fz``, ``Body1Pz``, or ``Rod2TenA``; the value equals what the
      standalone driver writes for that channel at the same state.

      :rtype: float
      :raises ValueError: if the token is empty or longer than 64 characters.
      :raises CableDynError: if the token is malformed or names an object the model lacks.

   .. py:method:: channels(tokens)

      Evaluate several channels; returns a ``(len(tokens),)`` array.

   .. py:method:: step_held(dt)

      Advance one step with every coupled point held at its current position (zero velocity
      and acceleration); returns the Newton iteration count, as :meth:`step`.

   .. py:property:: time
      :type: float

      Simulation time (s): 0 after initialisation, advanced by every completed step.

   .. py:property:: deck
      :type: pathlib.Path | None

      Absolute path of the deck the model was initialised from.

Object views
~~~~~~~~~~~~

The classes returned by :attr:`CableDyn.lines`, :attr:`CableDyn.points`,
:attr:`CableDyn.bodies` and :attr:`CableDyn.rods`. They live in :mod:`cabledyn.objects`, which
imports without the shared library.

.. automodule:: cabledyn.objects
   :no-members:

.. autoclass:: cabledyn.objects.Line
   :members:
   :inherited-members:

.. autoclass:: cabledyn.objects.Point
   :members:
   :inherited-members:

.. autoclass:: cabledyn.objects.Body
   :members:
   :inherited-members:

.. autoclass:: cabledyn.objects.Rod
   :members:
   :inherited-members:

.. currentmodule:: cabledyn

.. py:function:: abi_version()

   Return the C ABI version of the loaded shared library. Importing the in-process API fails
   with :exc:`OSError` when the library's ABI version differs from the one this package
   supports.

   :rtype: int

.. py:function:: abi_minor()

   Return the extension level of ABI 1 of the loaded library: 1 when it has the object queries
   behind :attr:`CableDyn.lines` and :meth:`CableDyn.channel`, 0 for an older library.

   :rtype: int

.. py:function:: library_path()

   Return the absolute path of the loaded CableDyn shared library.

   :rtype: str

.. py:function:: version_string()

   Return the loaded library's human-readable version string.

   :rtype: str

Exceptions
----------

All exception classes are importable without the shared library, so in-process errors can be
caught in any installation. The hierarchy is:

* :exc:`DriverError` (a :exc:`RuntimeError`): standalone-driver and result-file errors, with
  subclasses :exc:`DriverNotFoundError`, :exc:`DriverExecutionError`, and
  :exc:`OutputFormatError`.
* :exc:`DeckFormatError` and :exc:`StudyFormatError` (both :exc:`ValueError`): malformed input
  decks and study manifests.
* :exc:`DeckReferenceError` (a :exc:`ValueError`): a :class:`DeckModel` edit that would leave
  an object referenced by other deck objects, or a row that names an object of another model.
* :exc:`StudyOutputError` (an :exc:`OSError`): a study ran its cases but could not write its
  summary files.
* :exc:`CableDynError` (a :exc:`RuntimeError`): a failed in-process call, with subclass
  :exc:`ConvergenceError`.

.. autoexception:: cabledyn.DriverError

.. autoexception:: cabledyn.DriverNotFoundError

.. autoexception:: cabledyn.DriverExecutionError

.. autoexception:: cabledyn.OutputFormatError

.. autoexception:: cabledyn.DeckFormatError

.. autoexception:: cabledyn.StudyFormatError

.. autoexception:: cabledyn.StudyOutputError

.. py:exception:: CableDynError(call, status, detail)

   A CableDyn C API call failed.

   .. py:attribute:: call
      :type: str

      Name of the C function that failed, for example ``"CableDyn_Step"``.

   .. py:attribute:: status
      :type: int

      The C ABI status code: 1 ``BAD_HANDLE``, 2 ``ALLOC_FAIL``, 3 ``BAD_INPUT``,
      4 ``SOLVE_FAIL``, or 5 ``NOT_INITIALIZED``.

   .. py:attribute:: detail
      :type: str

      The library's diagnostic message for the handle, possibly empty.

.. py:exception:: ConvergenceError(call, detail, *, n_iter, stalled)

   An implicit step finished without meeting the Newton tolerance. A subclass of
   :exc:`CableDynError` with :attr:`~CableDynError.status` 4 (``SOLVE_FAIL``).

   The model has already advanced to ``t + dt`` and holds the best iterate found, so retrying
   the same step would advance time twice.

   .. py:attribute:: n_iter
      :type: int

      Number of Newton iterations performed.

   .. py:attribute:: stalled
      :type: bool

      Whether the line search stalled.
