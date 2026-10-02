# SPDX-License-Identifier: Apache-2.0
"""Environment, sea state, current, seabed, and prescribed motion.

These objects hold the deck's environment ``OPTIONS``. An option property
left unset (``None``) is not written, so the solver's own default applies;
the defaults are deliberately not restated here.
"""

from __future__ import annotations

from cabledyn.project.appearance import EnvironmentAppearance, SeabedAppearance
from cabledyn.project.base import ModelObject, model_type
from cabledyn.project.descriptors import (
    Child,
    Children,
    Choice,
    FilePath,
    OptionalInteger,
    OptionalQuantity,
    OptionalText,
    OptionalVec3,
    Quantity,
    Role,
    Strategy,
    Vec3,
)
from cabledyn.project.issues import Issue, Severity
from cabledyn.project.units import (
    ACCELERATION,
    ANGLE,
    DENSITY,
    DIMENSIONLESS,
    LENGTH,
    PRESSURE_PER_LENGTH,
    PRESSURE_TIME_PER_LENGTH,
    TIME,
    VELOCITY,
)

__all__ = [
    "BathymetrySeabed",
    "CurrentModel",
    "Environment",
    "FlatSeabed",
    "FloaterMotionRecord",
    "FloaterRAO",
    "JonswapWave",
    "MotionFile",
    "MotionSource",
    "MultiTrainSea",
    "NoCurrent",
    "NoMotion",
    "NoWaves",
    "OchiHubbleWave",
    "ProfileCurrent",
    "RegularWave",
    "Seabed",
    "SpectrumWave",
    "UniformCurrent",
    "WaveModel",
    "WaveTrain",
]

_SEA = "Sea state"


# --------------------------------------------------------------------------- waves


class WaveModel(ModelObject):
    """Base class of the wave models (the deck's ``waves`` and ``wavetrain`` rows)."""

    type_label = "Waves"
    abstract = True


@model_type("wave.none")
class NoWaves(WaveModel):
    """Still water."""

    type_label = "No waves"


@model_type("wave.regular")
class RegularWave(WaveModel):
    """A regular wave: linear (Airy) or nonlinear stream-function theory."""

    type_label = "Regular wave"
    theory = Choice(
        ("airy", "stream", "dean"), group=_SEA, doc="Airy, or Dean stream-function theory."
    )
    height = Quantity(LENGTH, 2.0, limit="wave_height", group=_SEA, doc="Wave height H.")
    period = Quantity(TIME, 10.0, limit="wave_period", group=_SEA, doc="Wave period T.")
    direction = Quantity(ANGLE, 0.0, limit="coordinate", group=_SEA, doc="Propagation direction.")


@model_type("wave.jonswap")
class JonswapWave(WaveModel):
    """A JONSWAP sea state."""

    type_label = "JONSWAP spectrum"
    significant_height = Quantity(LENGTH, 2.0, limit="wave_height", group=_SEA, doc="Hs.")
    peak_period = Quantity(TIME, 10.0, limit="wave_period", group=_SEA, doc="Tp.")
    peak_enhancement = Quantity(DIMENSIONLESS, 3.3, minimum=1.0, group=_SEA, doc="Gamma.")
    direction = Quantity(ANGLE, 0.0, limit="coordinate", group=_SEA)


@model_type("wave.spectrum")
class SpectrumWave(WaveModel):
    """A two-parameter spectral sea state (Pierson-Moskowitz family, Torsethaugen)."""

    type_label = "Spectral sea"
    spectrum = Choice(("pm", "issc", "bretschneider", "torsethaugen"), group=_SEA)
    significant_height = Quantity(LENGTH, 2.0, limit="wave_height", group=_SEA, doc="Hs.")
    peak_period = Quantity(TIME, 10.0, limit="wave_period", group=_SEA, doc="Tp.")
    direction = Quantity(ANGLE, 0.0, limit="coordinate", group=_SEA)


@model_type("wave.ochi_hubble")
class OchiHubbleWave(WaveModel):
    """A two-peak Ochi-Hubble sea state."""

    type_label = "Ochi-Hubble spectrum"
    height_1 = Quantity(LENGTH, 2.0, limit="wave_height", group=_SEA, doc="Hs of the first peak.")
    period_1 = Quantity(TIME, 12.0, limit="wave_period", group=_SEA, doc="Tp of the first peak.")
    shape_1 = Quantity(DIMENSIONLESS, 3.0, group=_SEA, doc="Shape parameter of the first peak.")
    height_2 = Quantity(LENGTH, 1.0, minimum=0.0, group=_SEA, doc="Hs of the second peak.")
    period_2 = Quantity(TIME, 7.0, limit="wave_period", group=_SEA, doc="Tp of the second peak.")
    shape_2 = Quantity(DIMENSIONLESS, 1.0, group=_SEA, doc="Shape parameter of the second peak.")
    direction = Quantity(ANGLE, 0.0, limit="coordinate", group=_SEA)


@model_type("wave.train")
class WaveTrain(ModelObject):
    """One wave train of a multi-train sea (one deck ``wavetrain`` row)."""

    type_label = "Wave train"
    abstract = False
    wave = Strategy(WaveModel, RegularWave, group=_SEA, doc="The train's wave or spectrum.")
    spreading = OptionalQuantity(
        DIMENSIONLESS, group=_SEA, doc="cos-2s spreading exponent (spectral trains)."
    )

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if isinstance(self.wave, (NoWaves, MultiTrainSea)):
            found.append(Issue(Severity.ERROR, "a train needs a wave or spectrum", self, "wave", 2))
        elif isinstance(self.wave, RegularWave) and self.wave.theory != "airy":
            found.append(Issue(Severity.ERROR, "a regular train must be airy", self, "wave", 2))
        return found


@model_type("wave.multi_train")
class MultiTrainSea(WaveModel):
    """A sea made of several wave trains (swell and wind sea, for example)."""

    type_label = "Multi-train sea"
    trains = Children(WaveTrain, min_items=1, group=_SEA)


# --------------------------------------------------------------------------- current


class CurrentModel(ModelObject):
    """Base class of the current models (the deck's ``current`` row)."""

    type_label = "Current"
    abstract = True


@model_type("current.none")
class NoCurrent(CurrentModel):
    """No current."""

    type_label = "No current"


@model_type("current.uniform")
class UniformCurrent(CurrentModel):
    """A depth-uniform current."""

    type_label = "Uniform current"
    velocity = Vec3(VELOCITY, (0.5, 0.0, 0.0), limit="fluid_speed", group="Current")


@model_type("current.profile")
class ProfileCurrent(CurrentModel):
    """A current varying linearly between two depths."""

    type_label = "Two-level current"
    z_1 = Quantity(LENGTH, -50.0, limit="coordinate", group="Current", doc="First level z.")
    velocity_1 = Vec3(VELOCITY, limit="fluid_speed", group="Current")
    z_2 = Quantity(LENGTH, 0.0, limit="coordinate", group="Current", doc="Second level z.")
    velocity_2 = Vec3(VELOCITY, (0.5, 0.0, 0.0), limit="fluid_speed", group="Current")


# --------------------------------------------------------------------------- environment


@model_type("environment")
class Environment(ModelObject):
    """Water, gravity, waves and current."""

    type_label = "Environment"
    abstract = False
    gravity = OptionalQuantity(ACCELERATION, limit="gravity", group="Constants")
    water_density = OptionalQuantity(DENSITY, limit="water_density", group="Constants")
    water_depth = OptionalQuantity(
        LENGTH,
        limit="water_depth",
        group="Water",
        doc="Flat-seabed depth; unset with a bathymetry file or a host-supplied depth.",
    )
    waves = Strategy(WaveModel, NoWaves, group=_SEA)
    current = Strategy(CurrentModel, NoCurrent, group="Current")
    ramp_time = OptionalQuantity(TIME, minimum=0.0, group=_SEA, doc="Start-up ramp of the waves.")
    wave_seed = OptionalInteger(group=_SEA, doc="Random-phase seed of spectral seas.")
    wave_spreading = OptionalQuantity(DIMENSIONLESS, group=_SEA, doc="cos-2s exponent.")
    wave_directions = OptionalInteger(group=_SEA, doc="Direction bins of a spread sea.")
    wave_components = OptionalInteger(group=_SEA, doc="Frequency components per direction.")
    stream_order = OptionalInteger(group=_SEA, doc="Fourier terms of a stream-function wave.")
    wave_time_step = OptionalQuantity(TIME, group=_SEA, doc="Resampling step of wave files.")
    water_kinematics = OptionalText(
        group="Kinematics", doc="Kinematics source: a file, a host selector, or a mode number."
    )
    current_mode = OptionalInteger(group="Kinematics", doc="Current-profile file mode.")
    appearance = Child(EnvironmentAppearance, role=Role.APPEARANCE)


# --------------------------------------------------------------------------- seabed


class Seabed(ModelObject):
    """Base class of the seabed descriptions; holds the contact properties."""

    type_label = "Seabed"
    abstract = True
    stiffness = OptionalQuantity(
        PRESSURE_PER_LENGTH, limit="coefficient", group="Contact", doc="Penalty stiffness base."
    )
    damping = OptionalQuantity(
        PRESSURE_TIME_PER_LENGTH, limit="coefficient", group="Contact", doc="Normal damping base."
    )
    friction = OptionalQuantity(DIMENSIONLESS, minimum=0.0, group="Friction")
    friction_axial = OptionalQuantity(DIMENSIONLESS, minimum=0.0, group="Friction")
    friction_lateral = OptionalQuantity(DIMENSIONLESS, minimum=0.0, group="Friction")
    appearance = Child(SeabedAppearance, role=Role.APPEARANCE)


@model_type("seabed.flat")
class FlatSeabed(Seabed):
    """A horizontal seabed at the environment's water depth."""

    type_label = "Flat seabed"


@model_type("seabed.bathymetry")
class BathymetrySeabed(Seabed):
    """A seabed from a structured x/y/depth bathymetry file."""

    type_label = "Bathymetry"
    file = FilePath(group="Bathymetry", doc="Bathymetry file.")

    def invariants(self) -> list[Issue]:
        found = super().invariants()
        if not self.file:
            found.append(Issue(Severity.ERROR, "needs a bathymetry file", self, "file", 2))
        return found


# --------------------------------------------------------------------------- motion


class MotionSource(ModelObject):
    """Base class of the prescribed-motion sources of a standalone run.

    The deck keywords ``motionFile``, ``vesselMotion``, ``vesselRAO`` and
    ``vesselRef`` map to :class:`MotionFile`, :class:`FloaterMotionRecord`,
    :class:`FloaterRAO` and ``reference``.
    """

    type_label = "Prescribed motion"
    abstract = True
    reference = OptionalVec3(
        LENGTH, group="Motion", doc="Floater reference point, rotation centre and RAO origin."
    )


@model_type("motion.none")
class NoMotion(MotionSource):
    """No prescribed motion."""

    type_label = "No prescribed motion"


@model_type("motion.file")
class MotionFile(MotionSource):
    """Point, rod-end and body motion from a time-series file."""

    type_label = "Motion file"
    file = FilePath(group="Motion")


@model_type("motion.floater_record")
class FloaterMotionRecord(MotionSource):
    """A 6-DOF floater record moving every floater-borne point rigidly."""

    type_label = "Floater motion record"
    file = FilePath(group="Motion")


@model_type("motion.floater_rao")
class FloaterRAO(MotionSource):
    """Floater motion as the RAO response to the waves."""

    type_label = "Floater RAO"
    file = FilePath(group="Motion")
