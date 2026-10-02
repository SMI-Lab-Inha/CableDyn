.. SPDX-License-Identifier: Apache-2.0

Project object model
====================

.. module:: cabledyn.project

:mod:`cabledyn.project` describes a complete CableDyn model as a graph of plain Python objects.
It is the model layer of the planned graphical front end and is also usable from scripts. It
has no GUI-toolkit dependency.

The package is under development: its classes and the project-file schema may change before
they are used by a released front end.

Design
------

Objects and properties
   Every object is a :class:`ModelObject` with a stable ``uid``, a ``name``, a ``parent`` and
   typed properties declared with descriptors (:class:`Quantity`, :class:`Integer`,
   :class:`Choice`, :class:`Vec3`, :class:`Ref`, :class:`Children`, :class:`Strategy`, ...).
   Quantities are stored in SI units (angles in degrees, as in the deck); non-finite numbers
   are refused. Each property has a role: ``physics`` (written to the deck), ``construction``
   (library and construction data such as a chain grade or a floater's hydrodynamic database;
   not written today), ``appearance`` (for the 3D views; never written), ``derived`` or
   ``meta``. :meth:`ModelObject.properties` lists the descriptors, which carry label, group,
   unit dimension, limits and documentation for editors. A value whose meaning changes with its
   sign in the deck is split: the axial damping is either a damping in N s or a dimensionless
   damping ratio (:class:`AxialModel`).

Identity and ownership
   Every object has one owner and a uid that is unique in its project: inserting an object whose
   uid is already used, or making an object own one of its owners, is refused.
   :meth:`ModelObject.from_dict` gives a copy fresh uids by default; pass ``within=project``
   so the copy keeps its references to the project's objects (paste).

References, renaming and removal
   Objects refer to each other by reference, never by name or deck id, so renaming cannot break
   a link. Removing a referenced object either fails with :class:`ObjectInUseError` or cascades
   (:func:`plan_removal`): a referrer whose required reference is lost is removed with it
   (:meth:`ModelObject.required_references`, which can depend on the referrer's state, such as
   a body rod's body), an optional reference is cleared, and a line left without sections is
   removed. Weak references, such as group memberships, never block a removal and are cleared.
   :meth:`Project.used_by` answers "used by" questions.

Events and commands
   A change publishes :class:`PropertyChanged`, :class:`ReferenceChanged`,
   :class:`ChildAdded`, :class:`ChildRemoved` or :class:`ChildMoved` on the object and on every
   ancestor, so one subscription on the :class:`Project` sees the whole model. Interactive edits
   are :class:`Command` objects pushed onto a :class:`CommandStack`, which provides undo, redo,
   merging of drags, grouping (:meth:`CommandStack.macro`) and a clean (saved) marker. Undo
   restores the project's :meth:`~ModelObject.to_dict` exactly. A change listener must not push,
   undo or redo while a command is being applied; the stack refuses it.

Validation
   Three layers (:func:`validate_project`): property checks with the numeric limits shared with
   the deck reader (:mod:`cabledyn.project.schema`), object invariants (references resolve,
   lines have sections, arc locations lie on the line, names the deck needs are unique), and the
   native deck rules, applied by writing the deck and reading it with :class:`cabledyn.DeckFile`.
   A deck error is mapped back to the object whose row failed. No deck rule is restated in the
   object model.

Units
   :class:`UnitSystem` holds the display unit chosen for each :class:`Dimension`; values are
   converted only for display. :func:`cabledyn.project.units.parse_quantity` reads text such as
   ``"1.2 MN"``.

Floaters
   The object model uses the term *floater*: :class:`Floater` bodies, :class:`FloaterPoint`,
   the shared :class:`FloaterType` library entry, and the motion sources
   :class:`FloaterMotionRecord` and :class:`FloaterRAO`. The deck keeps its own ``Vessel``
   keywords (``Vessel`` points and bodies, ``vesselMotion``, ``vesselRAO``, ``vesselRef``), which
   the deck adapters map. The ``hydrodynamics`` property of :class:`FloaterType` is the slot for
   a hydrodynamic database: new database kinds register as :class:`HydroDatabase` subclasses.

Deck adapters and the project file
----------------------------------

:class:`DeckReader` reads a deck through :class:`cabledyn.DeckModel` and builds the object
graph; :class:`DeckWriter` builds a fresh :class:`~cabledyn.DeckModel` and renders it. Deck ids,
keyword spellings, the order and descriptions of ``OPTIONS`` rows are kept as hints, so a
deck read and written back gives exactly ``DeckModel.load(deck).to_text()``. Comments and
column alignment of the original file are not kept. ``OPTIONS`` rows the model does not type are
kept verbatim as :class:`ExtraOption` objects. An extra row that sets a typed option (a duplicate
kept from a deck, for example) never overrides the typed property: the writer places it before
the typed row, and validation warns about it. The writer returns an :class:`IdMap` between
objects and deck ids.

.. code-block:: python

   from cabledyn.project import CommandStack, DeckReader, DeckWriter, ProjectStore
   from cabledyn.project.commands import SetProperty

   project = DeckReader.read("examples/lazy_wave_buoyancy_modules.dat")
   stack = CommandStack()
   stack.push(SetProperty(project.lines[0].sections[0], "length", 120.0))
   stack.undo()
   ProjectStore.save(project, "lazy_wave.cdproj")
   text = DeckWriter(project).to_text()

:class:`ProjectStore` writes a ``.cdproj`` zip archive (or an unzipped folder) with
``project.json`` (the whole graph, physics and appearance, as strict JSON with a schema
version), ``model.dat``, ``assets/`` and ``results.json``. ``model.dat`` is written only when
the project passes every validation layer. Its relative side-file paths (motion, bathymetry,
WaterKin and Syrope files) are rewritten to resolve from the project folder (a folder project)
or from the folder that holds the archive (a zip project), and side files that cannot be found
are listed in ``project.json``. Each file is replaced atomically. Objects, properties and
references the reader does not understand are kept and written back unchanged.

Extending the model
-------------------

A new object kind is a :class:`ModelObject` subclass registered with :func:`model_type`:

.. code-block:: python

   from cabledyn.project import HydroDatabase, model_type
   from cabledyn.project.descriptors import FilePath

   @model_type("myplugin.hydro.tabulated")
   class TabulatedDatabase(HydroDatabase):
       """Hydrodynamic coefficients read from a table."""

       table = FilePath()

API reference
-------------

Base classes, descriptors and events
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: cabledyn.project.ModelObject
   :members:

.. autoclass:: cabledyn.project.ObjectCollection
   :members:

.. autofunction:: cabledyn.project.model_type

.. autoclass:: cabledyn.project.TypeRegistry
   :members:

.. autoexception:: cabledyn.project.UnknownTypeError

.. autoclass:: cabledyn.project.SerialContext
   :members:

.. autoclass:: cabledyn.project.Role
   :members:

.. autoclass:: cabledyn.project.Property
   :members: describe, coerce, check

.. autoclass:: cabledyn.project.Quantity
.. autoclass:: cabledyn.project.OptionalQuantity
.. autoclass:: cabledyn.project.Integer
.. autoclass:: cabledyn.project.OptionalInteger
.. autoclass:: cabledyn.project.Bool
.. autoclass:: cabledyn.project.OptionalBool
.. autoclass:: cabledyn.project.Choice
.. autoclass:: cabledyn.project.Text
.. autoclass:: cabledyn.project.OptionalText
.. autoclass:: cabledyn.project.FilePath
.. autoclass:: cabledyn.project.TextList
.. autoclass:: cabledyn.project.Colour
.. autoclass:: cabledyn.project.OptionalColour
.. autoclass:: cabledyn.project.Vector
.. autoclass:: cabledyn.project.Vec3
.. autoclass:: cabledyn.project.OptionalVec3
.. autoclass:: cabledyn.project.JsonValue
.. autoclass:: cabledyn.project.Ref
.. autoclass:: cabledyn.project.RefList
.. autoclass:: cabledyn.project.Child
.. autoclass:: cabledyn.project.Strategy
   :members: options
.. autoclass:: cabledyn.project.Children

.. autoclass:: cabledyn.project.EventHub
   :members:
.. autoclass:: cabledyn.project.Subscription
   :members:
.. autoclass:: cabledyn.project.Event
.. autoclass:: cabledyn.project.PropertyChanged
.. autoclass:: cabledyn.project.ReferenceChanged
.. autoclass:: cabledyn.project.ChildAdded
.. autoclass:: cabledyn.project.ChildRemoved
.. autoclass:: cabledyn.project.ChildMoved

Project and model objects
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: cabledyn.project.Project
   :members: find, objects_of, used_by, validate_all

.. autoclass:: cabledyn.project.Environment
.. autoclass:: cabledyn.project.WaveModel
.. autoclass:: cabledyn.project.NoWaves
.. autoclass:: cabledyn.project.RegularWave
.. autoclass:: cabledyn.project.JonswapWave
.. autoclass:: cabledyn.project.SpectrumWave
.. autoclass:: cabledyn.project.OchiHubbleWave
.. autoclass:: cabledyn.project.MultiTrainSea
.. autoclass:: cabledyn.project.WaveTrain
.. autoclass:: cabledyn.project.CurrentModel
.. autoclass:: cabledyn.project.NoCurrent
.. autoclass:: cabledyn.project.UniformCurrent
.. autoclass:: cabledyn.project.ProfileCurrent
.. autoclass:: cabledyn.project.Seabed
.. autoclass:: cabledyn.project.FlatSeabed
.. autoclass:: cabledyn.project.BathymetrySeabed
.. autoclass:: cabledyn.project.MotionSource
.. autoclass:: cabledyn.project.NoMotion
.. autoclass:: cabledyn.project.MotionFile
.. autoclass:: cabledyn.project.FloaterMotionRecord
.. autoclass:: cabledyn.project.FloaterRAO
.. autoclass:: cabledyn.project.AnalysisSettings
.. autoclass:: cabledyn.project.ExtraOption

.. autoclass:: cabledyn.project.LineType
.. autoclass:: cabledyn.project.GenericLineType
.. autoclass:: cabledyn.project.ChainType
.. autoclass:: cabledyn.project.WireRopeType
.. autoclass:: cabledyn.project.FibreRopeType
.. autoclass:: cabledyn.project.CableType
.. autoclass:: cabledyn.project.AxialModel
.. autoclass:: cabledyn.project.LinearAxial
.. autoclass:: cabledyn.project.ViscoelasticAxial
.. autoclass:: cabledyn.project.SyropeAxial
.. autoclass:: cabledyn.project.BendingModel
.. autoclass:: cabledyn.project.RodType
.. autoclass:: cabledyn.project.FloaterType
.. autoclass:: cabledyn.project.HydroDatabase
.. autoclass:: cabledyn.project.NoHydroDatabase
.. autoclass:: cabledyn.project.HydroDatabaseFile

.. autoclass:: cabledyn.project.Body
.. autoclass:: cabledyn.project.Rigid6Body
.. autoclass:: cabledyn.project.Point3Body
.. autoclass:: cabledyn.project.Floater
.. autoclass:: cabledyn.project.Turbine
.. autoclass:: cabledyn.project.Rod
   :members: endpoint
.. autoclass:: cabledyn.project.LineEndTarget
.. autoclass:: cabledyn.project.RodEndpoint
.. autoclass:: cabledyn.project.Point
.. autoclass:: cabledyn.project.FixedPoint
.. autoclass:: cabledyn.project.FloaterPoint
.. autoclass:: cabledyn.project.FreePoint
.. autoclass:: cabledyn.project.ConnectPoint
.. autoclass:: cabledyn.project.TurbinePoint
.. autoclass:: cabledyn.project.BodyPoint
.. autoclass:: cabledyn.project.RodPoint

.. autoclass:: cabledyn.project.Line
   :members: length, segment_count
.. autoclass:: cabledyn.project.Section
.. autoclass:: cabledyn.project.EndConnection
.. autoclass:: cabledyn.project.Attachment
   :members: is_series, locations
.. autoclass:: cabledyn.project.BuoyancyModule
.. autoclass:: cabledyn.project.ClumpWeight
.. autoclass:: cabledyn.project.LineAncillary
.. autoclass:: cabledyn.project.BendStiffener
.. autoclass:: cabledyn.project.BendRestrictor
.. autoclass:: cabledyn.project.TouchdownProtection
.. autoclass:: cabledyn.project.EquivalentBuoyancy
.. autoclass:: cabledyn.project.SyropeIC
.. autoclass:: cabledyn.project.Failure
.. autoclass:: cabledyn.project.Control
.. autoclass:: cabledyn.project.ExternalLoad
.. autoclass:: cabledyn.project.OutputRequest
.. autoclass:: cabledyn.project.Channel
.. autoclass:: cabledyn.project.ObjectChannel
.. autoclass:: cabledyn.project.NamedChannel

Appearance and GUI data
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: cabledyn.project.LineAppearance
.. autoclass:: cabledyn.project.PlainLineAppearance
.. autoclass:: cabledyn.project.ChainAppearance
.. autoclass:: cabledyn.project.WireRopeAppearance
.. autoclass:: cabledyn.project.FibreRopeAppearance
.. autoclass:: cabledyn.project.CableAppearance
.. autoclass:: cabledyn.project.ModuleGeometry
.. autoclass:: cabledyn.project.MeshAsset
.. autoclass:: cabledyn.project.EnvironmentAppearance
.. autoclass:: cabledyn.project.SeabedAppearance
.. autoclass:: cabledyn.project.StudioData
   :members: unit_system, set_unit_system
.. autoclass:: cabledyn.project.Group
.. autoclass:: cabledyn.project.CameraBookmark

Commands and references
~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: cabledyn.project.CommandStack
   :members:
.. autoclass:: cabledyn.project.StackChanged
.. autoclass:: cabledyn.project.Command
   :members:
.. autoclass:: cabledyn.project.SetProperty
.. autoclass:: cabledyn.project.ReplaceStrategy
.. autoclass:: cabledyn.project.Reconnect
.. autoclass:: cabledyn.project.AddObject
.. autoclass:: cabledyn.project.RemoveObject
.. autoclass:: cabledyn.project.MoveInCollection
.. autoclass:: cabledyn.project.MacroCommand
.. autoclass:: cabledyn.project.SplitSection
.. autoclass:: cabledyn.project.ApplyLibraryItem
.. autoclass:: cabledyn.project.ImportDeck
.. autofunction:: cabledyn.project.plan_removal
.. autofunction:: cabledyn.project.referrers
.. autoclass:: cabledyn.project.RemovalPlan
.. autoexception:: cabledyn.project.ObjectInUseError

Validation and units
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: cabledyn.project.validate_project
.. autoclass:: cabledyn.project.Issue
   :members: is_error
.. autoclass:: cabledyn.project.Severity
   :members:
.. autoclass:: cabledyn.project.UnitSystem
   :members:
.. autoclass:: cabledyn.project.Dimension
.. autoclass:: cabledyn.project.Unit
   :members:
.. autoexception:: cabledyn.project.UnitError
.. autofunction:: cabledyn.project.units.convert
.. autofunction:: cabledyn.project.units.parse_quantity
.. autofunction:: cabledyn.project.units.unit
.. autofunction:: cabledyn.project.units.units_for
.. autoclass:: cabledyn.project.schema.DeckLimit
   :members:
.. autofunction:: cabledyn.project.schema.deck_limit

.. autodata:: cabledyn.project.schema.DECK_LIMITS
   :no-value:

Adapters
~~~~~~~~

.. autoclass:: cabledyn.project.DeckReader
   :members:
.. autoclass:: cabledyn.project.DeckWriter
   :members:
.. autoclass:: cabledyn.project.IdMap
   :members:
.. autoexception:: cabledyn.project.DeckExportError
.. autoclass:: cabledyn.project.ProjectStore
   :members:
.. autoexception:: cabledyn.project.ProjectFormatError

Modules
~~~~~~~

Every public name above is importable from :mod:`cabledyn.project`; the classes are defined in
these modules.

.. py:module:: cabledyn.project.descriptors
   :synopsis: Property descriptors.
.. py:module:: cabledyn.project.events
   :synopsis: Change events.
.. py:module:: cabledyn.project.commands
   :synopsis: Commands and the command stack.
.. py:module:: cabledyn.project.references
   :synopsis: Reverse references and removal cascades.
.. py:module:: cabledyn.project.validation
   :synopsis: The validation layers.
.. py:module:: cabledyn.project.units
   :synopsis: Dimensions, units and display preferences.
.. py:module:: cabledyn.project.schema
   :synopsis: Deck limits shared with the deck reader.
