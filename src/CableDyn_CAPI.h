/* File: src/CableDyn_CAPI.h
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
 *
 * Public C declarations for the CableDyn ISO C binding.
 * The reference is doc/capi.rst (the C API reference). This header mirrors the
 * BIND(C) procedures implemented in src/CableDyn_CAPI.f90.
 *
 * Conventions shared by every entry point:
 *
 * - Scalar out-parameters (err_stat, converged, stalled, n_iter, the version
 *   components, and the void** of Create/Close) must point to valid storage.
 *   Only the handle, array, and string arguments are checked for NULL.
 * - Arrays are caller-owned and contiguous. Nodal positions, velocities, and
 *   fluid fields are interleaved xyz triples. Connectivity and fixed-DOF
 *   indices are 1-based. The Jacobians written by
 *   CableDyn_CalcOutputDerivatives are n-by-n column-major:
 *   J[i + j*n] = d(load_i)/d(q_j).
 * - Diagnostics are retained per handle (up to 1024 characters) and copied by
 *   CableDyn_GetLastError, truncated to message_len - 1 and NUL-terminated.
 * - Create and Close are thread-safe. InitDeck, InitLine, and InitLines are
 *   serialised process-wide (they parse files through the Fortran runtime), so
 *   concurrent initialisations of distinct handles succeed one after another.
 *   The remaining per-handle calls run concurrently on distinct handles. A
 *   handle must not be used by two threads at once, and must not be used while
 *   another thread closes it.
 * - The first Create applies the BLAS thread policy (one OpenBLAS thread unless
 *   the environment variable CABLEDYN_BLAS_THREADS says otherwise) and, in the
 *   Windows GNU build, loads openblas.dll. When that runtime cannot be loaded,
 *   InitDeck, InitLine and InitLines return CD_C_ALLOC_FAIL and GetLastError names
 *   the library, the locations tried and the Windows load error.
 */
#ifndef CABLEDYN_CAPI_H
#define CABLEDYN_CAPI_H

#include <stdbool.h>

#define CD_C_OK 0
#define CD_C_BAD_HANDLE 1
#define CD_C_ALLOC_FAIL 2
#define CD_C_BAD_INPUT 3
#define CD_C_SOLVE_FAIL 4
#define CD_C_NOT_INITIALIZED 5

#define CABLEDYN_CAPI_VERSION_MAJOR 0
#define CABLEDYN_CAPI_VERSION_MINOR 1
#define CABLEDYN_CAPI_VERSION_PATCH 0
#define CABLEDYN_CAPI_ABI_VERSION 1

#ifdef __cplusplus
extern "C" {
#endif

void CableDyn_GetVersion(int *major, int *minor, int *patch, int *abi_version);
void CableDyn_GetVersionString(char *version, int version_len);

void CableDyn_Create(void **handle, int *err_stat);
void CableDyn_Close(void **handle, int *err_stat);

void CableDyn_InitDeck(void *handle, const char *deck_path, int deck_path_len,
                       int *err_stat);

void CableDyn_InitLine(void *handle, int n_nodes, int n_elem, const double *q0,
                       const double *v0, const int *elem_conn, const double *l0,
                       const double *ea, const double *rho_a,
                       const int *fixed_dofs, int n_fixed, int tension_only,
                       double rho_inf, int *err_stat);

void CableDyn_InitLines(void *handle, int n_lines, const int *n_nodes,
                        const int *n_elem, const int *n_fixed, const double *q0,
                        const double *v0, const int *elem_conn,
                        const double *l0, const double *ea,
                        const double *rho_a, const int *fixed_dofs,
                        const int *coupled_map, int n_coupled_map,
                        int tension_only, double rho_inf, int *err_stat);

void CableDyn_UpdateStates(void *handle, const double *q_coupled,
                           const double *v_coupled, const double *a_coupled,
                           int n_coupled_dof, int *err_stat);

void CableDyn_UpdatePointFluidFields(void *handle, const double *fluid_velocity,
                                     const double *fluid_acceleration,
                                     const double *waterline_z, int n_point,
                                     double fluid_density, int *err_stat);

void CableDyn_Step(void *handle, double dt, const double *q_coupled,
                   const double *v_coupled, const double *a_coupled,
                   int n_coupled_dof, bool *converged, bool *stalled,
                   int *n_iter, int *err_stat);

void CableDyn_CalcOutput(void *handle, double *coupled_loads,
                         int n_coupled_dof, int *err_stat);

void CableDyn_CalcOutputDerivatives(void *handle, const double *q_coupled,
                                    const double *v_coupled,
                                    const double *a_coupled,
                                    int n_coupled_dof, double eps_fd,
                                    double *coupled_loads, double *dload_dq,
                                    double *dload_dv, double *dload_da,
                                    double *added_mass, int *err_stat);

void CableDyn_GetCoupledMotion(void *handle, double *q_coupled,
                               double *v_coupled, double *a_coupled,
                               int n_coupled_dof, int *err_stat);

void CableDyn_GetLastError(void *handle, char *message, int message_len);

bool CableDyn_IsInitialized(void *handle);
int CableDyn_NCoupledDOF(void *handle, int *err_stat);
/* Stored system point count: the column count CableDyn_UpdatePointFluidFields
 * expects. May exceed CableDyn_NCoupledDOF/3 when the deck has unbound
 * (output-only) Fixed/Coupled points. */
int CableDyn_NPoints(void *handle, int *err_stat);
int CableDyn_NLines(void *handle, int *err_stat);

/* ---- ABI 1, minor extension 1: in-process object queries ------------------
 *
 * Read the committed state of an initialized handle at any step. Every value
 * comes from the evaluator behind the deck OUTPUTS channels and the per-line
 * .Line<L>.p.out/.t.out files, so a query returns what the matching output
 * channel reports for the same state. Object indices are 0-based positions in
 * the handle's inventory: lines in ascending deck id, points in system order,
 * bodies and rods in solver order. Ids are deck ids (raw-line handles number
 * their lines 1..n). Buffers are caller-owned; where noted, NULL skips an
 * output. A deck with finite-EI lines, BODIES or RODS runs on the coupled
 * aggregate (it needs dtM); other decks keep the point route. Check
 * CableDyn_GetAbiMinor() >= 1 before using these on a library of unknown
 * vintage. */
#define CABLEDYN_CAPI_ABI_MINOR 1

#define CD_C_OBJ_LINE 1
#define CD_C_OBJ_POINT 2
#define CD_C_OBJ_BODY 3
#define CD_C_OBJ_ROD 4

/* Line quantities (the codes of the line-node output channels). */
#define CD_C_LINE_POSITION 1        /* 3 per node, xyz interleaved [m] */
#define CD_C_LINE_VELOCITY 2        /* 3 per node [m/s] */
#define CD_C_LINE_TENSION 3         /* 1 per node, Ten<L>N<k> [N] */
#define CD_C_LINE_ACCELERATION 4    /* 3 per node [m/s^2] */
#define CD_C_LINE_CURVATURE 5       /* 1 per node [1/m] */
#define CD_C_LINE_BEND_MOMENT 6     /* 1 per node [N.m], 0 when EI=0 */
#define CD_C_LINE_DECLINATION 7     /* 1 per node [deg] */
#define CD_C_LINE_AZIMUTH 8         /* 1 per node [deg] */
#define CD_C_LINE_SEGMENT_TENSION 9 /* 1 per segment, .t.out order [N] */

int CableDyn_GetAbiMinor(void);
int CableDyn_NObjects(void *handle, int kind, int *err_stat);
/* Lines: subtype 1 = finite-EI, 0 = EI=0. Points: n_nodes 1, subtype 1 Fixed,
 * 2 Coupled, 3 Free, 4 Connect. Bodies: n_nodes 1. Rods: n_nodes = NumSegs+1. */
void CableDyn_GetObjectInfo(void *handle, int kind, int index, int *id,
                            int *n_nodes, int *subtype, int *err_stat);
/* Nodes run End A -> End B; n_out must equal the value count exactly. */
void CableDyn_GetLineValues(void *handle, int index, int quantity, double *out,
                            int n_out, int *err_stat);
/* pos, vel, force: 3 doubles each (NULL skips). force is the resultant of
 * the forces the attached lines exert on the point (Point<P>F). */
void CableDyn_GetPointState(void *handle, int index, double *pos, double *vel,
                            double *force, int *err_stat);
/* 6 doubles each (NULL skips): pose [x y z rx ry rz] (m; deg, x-y'-z''),
 * vel [vx vy vz wx wy wz] (m/s; deg/s), acc (m/s^2; deg/s^2), wrench
 * [Fx Fy Fz Mx My Mz] about the reference point (N; N.m). */
void CableDyn_GetBodyState(void *handle, int index, double *pose, double *vel,
                           double *acc, double *wrench, int *err_stat);
/* nodes: 3 * n_nodes doubles, node 0 = End A; pose [x y z rx ry 0] of End A
 * with the axis roll and pitch from vertical; vel and wrench about End A as
 * for a body. NULL skips an output. */
void CableDyn_GetRodState(void *handle, int index, double *nodes, int n_nodes,
                          double *pose, double *vel, double *wrench,
                          int *err_stat);
/* Evaluate one deck OUTPUTS token ("FairTen1", "L2N5px", "Curv1N3",
 * "Point4Fz", "Body1Pz", "Rod2TenA", "TDP1s", ...) at the committed state. */
void CableDyn_EvalChannel(void *handle, const char *token, int token_len,
                          double *value, int *err_stat);

#ifdef __cplusplus
}
#endif

#endif
