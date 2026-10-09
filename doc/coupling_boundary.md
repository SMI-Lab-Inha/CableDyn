<!-- SPDX-License-Identifier: Apache-2.0 -->

# Coupling boundary

This page is for developers coupling CableDyn to a host code: OpenFAST (maintained by NLR, the
National Laboratory of the Rockies, formerly NREL), a CFD solver, or a custom co-simulation. It
defines what crosses the boundary, in which frame, and at what cadence. For
usage, see the [C API guide](capi.rst), the [Python package](python.rst), and the
[OpenFAST guide](openfast.rst).

CableDyn is one solver core with two coupling shells: the OpenFAST module (`CompMooring = 5`, the
MoorDyn-F counterpart) and an ISO C binding (the MoorDyn-C counterpart). Both exchange the same
quantities with the core: **kinematics in, loads out**. The core has no caller-specific logic; it
cannot tell whether its motion comes from an OpenFAST mesh or a C caller, or whether its fluid
field comes from its own wave model or from a host. Units, frames, and the load sign convention
are those of [Conventions](conventions.rst).

## Exchanged quantities

**Coupled degrees of freedom.** The exchange is translational. A handle exposes `n` coupled DOFs:
the `x, y, z` components of every coupled point (fairleads, vessel points, rigid-body attachment
points), interleaved as xyz triples in the point order fixed at initialisation. A point shared by
several lines appears once. `CableDyn_NCoupledDOF` returns `n`.

**Kinematics in.** Position, velocity, and acceleration of the coupled DOFs in the global frame.
For a step from `t` to `t + dt` they are the prescribed values at `t + dt`; the prescribed
acceleration enters the dynamic residual through the consistent mass.

**Fluid fields in.**

- Line nodes take their fluid kinematics from the deck environment (waves and current) in
  standalone and C API use, and from the host field in OpenFAST, where SeaState is sampled at the
  structural nodes and held over the step.
- Dynamic `Free`/`Connect` points with hydrodynamic properties accept an external point-level
  field through `CableDyn_UpdatePointFluidFields`: velocity and acceleration (3 × `n_point`),
  local surface elevation (`n_point`), and one fluid density. A point above its surface
  elevation, or with zero density, carries no fluid load.

**Loads out.** The force the cable system exerts on each coupled point: `n` values in the global
frame. It is the reaction `−(M·a + f_int − f_ext)` at the prescribed DOFs, so it includes the
line's inertia and hydrodynamic reaction at the coupled node (see [Theory](theory.rst)). The
OpenFAST shell also maps point loads to a rigid-body wrench `[Fx, Fy, Fz, Mx, My, Mz]` about the
platform reference point, and returns the end-connection moments of finite-EI lines on its load
mesh.

**Load derivatives out.** At a supplied coupled state, `CableDyn_CalcOutputDerivatives` returns
four n-by-n matrices:

| Matrix | Definition |
| --- | --- |
| `dload_dq` | ∂load/∂q, with the free DOFs of every line statically condensed at the current effective mass |
| `dload_dv` | ∂load/∂v, condensed the same way |
| `dload_da` | ∂load/∂a = −(M_cc − M_cf M_ff⁻¹ M_fc), where M includes added mass |
| `added_mass` | −`dload_da` |

All four are analytic, assembled over the whole system, and written column-major:
`J[i + j*n] = ∂load_i/∂x_j`. They include the structural, damping, drag, added-mass, and seabed
contributions. The `eps_fd` argument must be non-negative and does not affect the result. The
call leaves the handle at the supplied operating point; a caller that needs its previous state
must save and restore it.

**Two derivative contracts in OpenFAST.** During a nonlinear OpenFAST run the shell returns the
direct-feedthrough derivative at fixed committed internal state, which the host's input–output
solve requires. During linearisation it returns the re-equilibrated zero-frequency mooring
stiffness. The two are not interchangeable.

## Cadence

The caller owns the coupling step `dt` and advances the core once per step over `[t, t + dt]`.
Internal step subdivision (see [Theory](theory.rst)) is invisible to the caller: the step lands
on `t + dt` and the host time grid is unchanged. In OpenFAST the module may advance at its own
`dtM`, an integer multiple of the glue step, holding its committed loads between advances.

## C ABI (`CableDyn_CAPI.h`)

The public header declares the complete interface of C ABI version 1, including the object
queries of its minor extension 1 (ABI 1.1; check `CableDyn_GetAbiMinor() >= 1` before calling
them). Signatures and examples are in the [C API guide](capi.rst).

| Group | Functions |
| --- | --- |
| Version | `CableDyn_GetVersion`, `CableDyn_GetVersionString` |
| Lifecycle | `CableDyn_Create`, `CableDyn_Close`, `CableDyn_IsInitialized` |
| Initialisation | `CableDyn_InitDeck` (a supported `EI = 0` point-system deck, through the deck parser and static solve, or a deck with finite-EI lines, free bodies or rods on the coupled aggregate); `CableDyn_InitLine` (one meshed `EI = 0` line); `CableDyn_InitLines` (several lines with a compact coupled-DOF map) |
| Stepping | `CableDyn_UpdateStates` (set coupled kinematics and recompute internal accelerations), `CableDyn_UpdatePointFluidFields`, `CableDyn_Step` |
| Output | `CableDyn_CalcOutput`, `CableDyn_CalcOutputDerivatives`, `CableDyn_GetCoupledMotion` |
| Queries | `CableDyn_NCoupledDOF`, `CableDyn_NPoints`, `CableDyn_NLines`, `CableDyn_GetLastError` |
| Object queries (ABI 1.1) | `CableDyn_GetAbiMinor`, `CableDyn_NObjects`, `CableDyn_GetObjectInfo`, `CableDyn_GetLineValues`, `CableDyn_GetPointState`, `CableDyn_GetBodyState`, `CableDyn_GetRodState`, `CableDyn_EvalChannel`: the model's lines, points, bodies and rods, their committed state, and any `OUTPUTS` channel |

Rules shared by every entry point:

- Arrays are caller-owned and contiguous. Positions, velocities, and fluid fields are interleaved
  xyz triples; connectivity and fixed-DOF indices are 1-based.
- Calls with an `err_stat` argument return `CD_C_OK`, `CD_C_BAD_HANDLE`, `CD_C_ALLOC_FAIL`,
  `CD_C_BAD_INPUT`, `CD_C_SOLVE_FAIL`, or `CD_C_NOT_INITIALIZED`. Diagnostics are kept per handle
  (up to 1024 characters) and copied by `CableDyn_GetLastError`.
- `CableDyn_NPoints` returns the column count expected by `CableDyn_UpdatePointFluidFields`; it
  can exceed `NCoupledDOF/3` when a deck has output-only fixed or coupled points.
- Initialisation is transactional: a failed call leaves the handle as it was (uninitialised,
  or holding its previous model).

The interface ships as the versioned shared library `libcabledyn` (`cabledyn.dll` in an Intel
Fortran and MSVC build; CMake target `cabledyn_shared`) with the installed header; it reports
CableDyn 0.1.1 and C ABI version 1, minor extension 1. The Python package uses the same library.

### Concurrency

A process-local C11 atomic mutex protects the live-handle registry, so `Create` and `Close` are
thread-safe independently of OpenMP. When copies of one handle are closed concurrently, exactly
one close succeeds; the others are rejected without touching freed storage. A handle must not be
used by two threads at once, or while another thread closes it.

## OpenFAST module

The OpenFAST shell follows the MoorDyn-F module form (`Init`, `UpdateStates`, `CalcOutput`,
`End`): coupled-point kinematics arrive on an input point mesh and loads return on an output
point mesh. It supports:

- nonzero `PtfmInit` and correction iterations;
- checkpoint/restart, with continuous states (including constitutive and rigid-body states)
  stored in the OpenFAST state record;
- quasi-static linearisation;
- FAST.Farm shared moorings;
- line failures and active tension control of `EI = 0` lines;
- WaterKin files with independent selection of host wave and host current. SeaState's standard
  steady current is reconstructed on the same vertical nodes and weights SeaState uses, so it can
  be removed or kept exactly; a host field that cannot be separated is rejected.

The route matrix is in [Capabilities](capabilities.rst) and setup in [OpenFAST](openfast.rst).

## Module map

| Module | Role |
| --- | --- |
| `CableDyn_System` | multi-line owner: coupled-DOF map, points, bodies, aggregate loads and derivatives |
| `CableDyn_OpenFAST` | lifecycle shell over the system: update, step, output, derivatives, body wrench |
| `CableDyn_OpenFAST_Types`, `CableDyn_OpenFAST_Mesh` | registry-style state records and the point-mesh adapter |
| `CableDyn_OpenFAST_FMF`, `CableDyn_OpenFAST_HermiteFMF` | module-form lifecycle for `EI = 0` systems and finite-EI cables |
| `CableDyn_OpenFAST_Aggregate` | mixed-line, body, rod, and FAST.Farm aggregate used by `CompMooring = 5` |
| `CableDyn_RigidKinematics` | platform pose/rate/acceleration to point kinematics, and point loads to a wrench |
| `CableDyn_CAPI` | ISO C binding with opaque handles over the same shell |

The full source map is in the [source module map](api_reference.md).
