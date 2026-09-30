! File: tests/test_model.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_model
  !! Lifecycle gate for CableDyn_Model, the persistent Fortran core API described
  !! in ARCHITECTURE.md and doc/coupling_boundary.md. The test keeps to the
  !! validated EI=0 dynamics path and checks product-level behaviour: init computes a
  !! consistent state, query reports coupled DOFs, a no-load fixed-end line remains
  !! at rest, prescribed fairlead motion advances through the model boundary, and
  !! End releases the model cleanly.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_Bathymetry, ONLY: CD_BathymetryType, CD_Init_Bathymetry, CD_End_Bathymetry, CD_BATHY_OK
  USE CableDyn_Dynamic, ONLY: GenAlphaConfig
  USE CableDyn_SeabedContact, ONLY: CD_SEABED_CONTACT_BLEND
  USE CableDyn_Line, ONLY: CD_LineType, CD_LineSection, CD_Build_Line_Mesh, CD_Nodal_Seabed_Stiffness, CD_LINE_OK
  USE CableDyn_Static, ONLY: CableSolverConfig
  USE CableDyn_Model, ONLY: CD_ModelType, CD_Init_Model, CD_Init_Line_Model, CD_Step_Model, &
                            CD_Set_Model_Friction_Anchors, CD_Get_Model_Friction_Anchors, &
                            CD_Step_Model_Recovering, CD_End_Model, &
                            CD_Model_NCoupledDOF, CD_Get_Model_State, CD_Model_Is_Initialized, &
                            CD_Model_NDOF, CD_Model_NElem, CD_Get_Model_CoupledDofs, &
                            CD_Get_Model_CoupledMotion, CD_Get_Model_Tension, CD_Get_Model_EndForces, &
                            CD_Calc_Model_CoupledLoads, CD_Calc_Model_CoupledKinematicDerivatives, &
                            CD_Calc_Model_CoupledAccelDerivative, &
                            CD_Update_Model_Hydro_Fields, &
                            CD_Update_Model_External_Loads, CD_Update_Model_CoupledMotion, &
                            CD_Recompute_Model_Acceleration, &
                            CD_Copy_Model, &
                            CD_MODEL_OK, CD_MODEL_BADINPUT, CD_MODEL_SOLVEFAIL, CD_MODEL_NOT_INITIALIZED, &
                            CD_Update_Model_SegmentLength
  USE CableDyn_Loads, ONLY: CD_Assemble_Distributed_Load
  USE CableDyn_Assemble, ONLY: CD_Assemble_Cable_Mass
  USE CableDyn_Hydro, ONLY: CD_Cable_Added_Mass_Matrix
  USE CableDyn_Linalg, ONLY: CD_Solve_Dense_As_Banded_Multiple
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN, IEEE_IS_FINITE

  USE, INTRINSIC :: IEEE_EXCEPTIONS, ONLY: IEEE_USUAL, IEEE_GET_HALTING_MODE, IEEE_SET_HALTING_MODE, IEEE_SET_FLAG
  IMPLICIT NONE
  LOGICAL :: fp_halt(3)   ! saved IEEE halting modes around deliberately non-finite inputs

  INTEGER :: nfail
  nfail = 0

  CALL case_lifecycle_static_step()
  CALL case_segment_length_update()
  CALL case_copy_model_preserves_state()
  CALL case_prescribed_step_updates_state()
  CALL case_recovery_common_path()
  CALL case_recovery_subdivision_path()
  CALL case_coupled_load_sign()
  CALL case_coupled_motion_update()
  CALL case_external_load_update()
  CALL case_recompute_acceleration_after_load_update()
  CALL case_seabed_coupled_load()
  CALL case_seabed_normal_damping_load()
  CALL case_seabed_normal_damping_tangent_fd()
  CALL case_seabed_friction_load()
  CALL case_seabed_friction_tangent_fd()
  CALL case_damping_model_load()
  CALL case_morison_drag_model_load()
  CALL case_morison_drag_element_coefficients()
  CALL case_froude_krylov_model_load()
  CALL case_buoyancy_recovery_model_load()
  CALL case_hydro_field_update()
  CALL case_hydro_damping_kinematic_derivatives_fd()
  CALL case_prescribed_hydro_added_mass_step()
  CALL case_banded_model_load_matches_dense_structured_bathymetry()
  CALL case_added_mass_model_accel()
  CALL case_accel_derivative_band_parity()
  CALL case_line_model_init_from_sections()
  CALL case_grounded_line_model_with_seabed()
  CALL case_line_model_fail_closed()
  CALL case_end_force_is_coupled_load()
  CALL case_fail_closed()

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: CableDyn_Model persistent lifecycle'

CONTAINS

  INCLUDE 'nan_max_abs.inc'

  SUBROUTINE straight_cable(ne, conn, l0, ea, rho_a, q0, fixed)
    !! Unit-spaced, pre-tensioned cable with both endpoints prescribed.
    INTEGER, INTENT(IN) :: ne
    INTEGER, INTENT(OUT) :: conn(2, ne), fixed(6)
    REAL(wp), INTENT(OUT) :: l0(ne), ea(ne), rho_a(ne), q0(3*(ne + 1))
    INTEGER :: e, i, ndof
    ndof = 3*(ne + 1)
    DO e = 1, ne
      conn(:, e) = [e, e + 1]
      l0(e) = 0.95_wp
      ea(e) = 1.0e6_wp
      rho_a(e) = 10.0_wp
    END DO
    DO i = 1, ne + 1
      q0(3*i - 2:3*i) = [REAL(i - 1, wp), 0.0_wp, 0.0_wp]
    END DO
    fixed = [1, 2, 3, ndof - 2, ndof - 1, ndof]
  END SUBROUTINE straight_cable

  SUBROUTINE case_lifecycle_static_step()
    !! A straight pretensioned line at rest with no external load remains exactly at
    !! rest after one model step, and the lifecycle/coupled-DOF queries are valid.
    INTEGER, PARAMETER :: NE = 4, NN = NE + 1, NDOF = 3*NN
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(6), coupled(6), es, n_iter, nc, nd, nelem
    REAL(wp) :: l0(NE), ea(NE), rho_a(NE), q0(NDOF), v0(NDOF), f_ext(NDOF)
    REAL(wp) :: q(NDOF), v(NDOF), a(NDOF)
    REAL(wp) :: loads(6)
    REAL(wp) :: dload_dq(6, 6), dload_dv(6, 6), dload_da(6, 6)
    LOGICAL :: conv, stalled
    CHARACTER(200) :: em

    CALL straight_cable(NE, conn, l0, ea, rho_a, q0, fixed)
    v0 = 0.0_wp
    f_ext = 0.0_wp
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK .AND. CD_Model_Is_Initialized(model), 'lifecycle:init')
    nc = CD_Model_NCoupledDOF(model, es, em)
    CALL require(es == CD_MODEL_OK .AND. nc == SIZE(fixed), 'lifecycle:ncoupled')
    nd = CD_Model_NDOF(model, es, em)
    CALL require(es == CD_MODEL_OK .AND. nd == NDOF, 'lifecycle:ndof')
    nelem = CD_Model_NElem(model, es, em)
    CALL require(es == CD_MODEL_OK .AND. nelem == NE, 'lifecycle:nelem')
    CALL CD_Get_Model_CoupledDofs(model, coupled, es, em)
    CALL require(es == CD_MODEL_OK .AND. ALL(coupled == fixed), 'lifecycle:coupled-map')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'lifecycle:loads')
    CALL require(ALLOCATED(model%dynamic_workspace%R) .AND. ALLOCATED(model%dynamic_workspace%load_eval), &
                 'lifecycle:load-query-workspace-owned')
    CALL CD_Calc_Model_CoupledKinematicDerivatives(model, dload_dq, dload_dv, es, em)
    CALL require(es == CD_MODEL_OK, 'lifecycle:kinematic-derivatives')
    CALL CD_Calc_Model_CoupledAccelDerivative(model, dload_da, es, em)
    CALL require(es == CD_MODEL_OK, 'lifecycle:accel-derivative')
    CALL require(ALLOCATED(model%dynamic_workspace%Kt_eval) .AND. &
                 ALLOCATED(model%dynamic_workspace%eff_free) .AND. &
                 ALLOCATED(model%dynamic_workspace%eff), 'lifecycle:derivative-workspace-owned')

    CALL CD_Step_Model(model, 0.02_wp, conv, stalled, n_iter, es, em)
    CALL require(es == CD_MODEL_OK .AND. conv .AND. .NOT. stalled, 'lifecycle:step')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK, 'lifecycle:get')
    CALL require(nan_max_abs(q - q0) < 1.0e-12_wp, 'lifecycle:q-at-rest')
    CALL require(nan_max_abs(v) < 1.0e-12_wp, 'lifecycle:v-at-rest')
    CALL require(nan_max_abs(a) < 1.0e-12_wp, 'lifecycle:a-at-rest')

    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK .AND. .NOT. CD_Model_Is_Initialized(model), 'lifecycle:end')
    CALL require(.NOT. ALLOCATED(model%dynamic_workspace%R) .AND. &
                 .NOT. ALLOCATED(model%dynamic_workspace%load_eval), 'lifecycle:load-query-workspace-released')
    CALL require(.NOT. ALLOCATED(model%dynamic_workspace%Kt_eval) .AND. &
                 .NOT. ALLOCATED(model%dynamic_workspace%eff_free) .AND. &
                 .NOT. ALLOCATED(model%dynamic_workspace%eff), 'lifecycle:derivative-workspace-released')
  END SUBROUTINE case_lifecycle_static_step

  SUBROUTINE case_segment_length_update()
    !! Active line control (the CtrlChan mechanism): CD_Update_Model_SegmentLength
    !! must leave the model INDISTINGUISHABLE from one BUILT with the new length --
    !! the stale-cache killer (any l0-derived quantity cached at init and not
    !! refreshed breaks the bitwise step equivalence). Also: paying out a taut
    !! segment drops its tension (sign), and the validation battery fails closed.
    INTEGER, PARAMETER :: NE = 4, NN = NE + 1, NDOF = 3*NN
    TYPE(CD_ModelType) :: ma, mb
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(6), es, n_iter
    REAL(wp) :: l0(NE), l0b(NE), ea(NE), rho_a(NE), q0(NDOF), v0(NDOF), f_ext(NDOF)
    REAL(wp) :: qa(NDOF), va(NDOF), aa(NDOF), qb(NDOF), vb(NDOF), ab(NDOF)
    REAL(wp) :: ten0(NE), ten1(NE), tenb(NE), loads_a(6), loads_b(6)
    LOGICAL :: conv, stalled
    CHARACTER(200) :: em

    CALL straight_cable(NE, conn, l0, ea, rho_a, q0, fixed)
    v0 = 0.0_wp
    f_ext = 0.0_wp
    ! a small transverse push on an interior node so the twin step is a real
    ! dynamic trajectory, not a fixed point
    f_ext(3*2 - 1) = 2.0_wp
    l0b = l0
    l0b(NE) = l0(NE) + 0.02_wp

    CALL CD_Init_Model(ma, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK, 'seglen:init-a')
    CALL CD_Init_Model(mb, q0, v0, conn, l0b, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK, 'seglen:init-b')

    ! sign: paying out the taut end segment drops its tension, others unchanged
    CALL CD_Get_Model_Tension(ma, ten0, es, em)
    CALL require(es == CD_MODEL_OK, 'seglen:tension-0')
    CALL CD_Update_Model_SegmentLength(ma, NE, l0b(NE), 0.0_wp, es, em)
    CALL require(es == CD_MODEL_OK, 'seglen:update: '//TRIM(em))
    CALL CD_Get_Model_Tension(ma, ten1, es, em)
    CALL require(es == CD_MODEL_OK, 'seglen:tension-1')
    CALL require(ten1(NE) < ten0(NE), 'seglen:payout-drops-end-tension')
    CALL require(nan_max_abs(ten1(1:NE - 1) - ten0(1:NE - 1)) <= 0.0_wp, 'seglen:other-segments-untouched')

    ! equivalence: the updated model matches the fresh-built twin bitwise --
    ! tensions, coupled loads (the re-assembled mass cache), and a full dynamic
    ! step under the same trajectory
    CALL CD_Get_Model_Tension(mb, tenb, es, em)
    CALL require(nan_max_abs(ten1 - tenb) <= 0.0_wp, 'seglen:tension-equals-fresh-build')
    CALL CD_Calc_Model_CoupledLoads(ma, loads_a, es, em)
    CALL require(es == CD_MODEL_OK, 'seglen:loads-a')
    CALL CD_Calc_Model_CoupledLoads(mb, loads_b, es, em)
    CALL require(es == CD_MODEL_OK, 'seglen:loads-b')
    CALL require(nan_max_abs(loads_a - loads_b) <= 0.0_wp, 'seglen:coupled-loads-equal-fresh-build')
    CALL CD_Step_Model(ma, 0.02_wp, conv, stalled, n_iter, es, em)
    CALL require(es == CD_MODEL_OK .AND. conv, 'seglen:step-a')
    CALL CD_Step_Model(mb, 0.02_wp, conv, stalled, n_iter, es, em)
    CALL require(es == CD_MODEL_OK .AND. conv, 'seglen:step-b')
    CALL CD_Get_Model_State(ma, qa, va, aa, es, em)
    CALL CD_Get_Model_State(mb, qb, vb, ab, es, em)
    CALL require(nan_max_abs(qa - qb) <= 0.0_wp .AND. nan_max_abs(va - vb) <= 0.0_wp .AND. &
                 nan_max_abs(aa - ab) <= 0.0_wp, 'seglen:step-equals-fresh-build-bitwise')

    ! DISTRIBUTED load tracks the length: with a stored per-length load, the
    ! updated model's f_ext must equal a fresh build's assembly at the new length
    ! (each end node carries 0.5*l0*w), proven by bitwise coupled-loads and step
    ! equality on WEIGHTED twins -- a regression against the paid-out segment
    ! running with the old net weight.
    BLOCK
      TYPE(CD_ModelType) :: mw0, mw1
      REAL(wp) :: wpl(3, NE), fw0(NDOF), fw1(NDOF), lw_a(6), lw_b(6)
      REAL(wp) :: qw0(NDOF), vw0(NDOF), aw0(NDOF), qw1(NDOF), vw1(NDOF), aw1(NDOF)
      wpl = 0.0_wp
      wpl(3, :) = -7.5_wp
      ! Assemble through the PRODUCTION assembler (never a hand-rolled twin): the
      ! bitwise fresh-build claim needs the exact arithmetic the model uses.
      CALL CD_Assemble_Distributed_Load(conn, l0, wpl, fw0, es, em)
      CALL require(es == 0, 'seglen:weighted-assemble-base')
      CALL CD_Assemble_Distributed_Load(conn, l0b, wpl, fw1, es, em)
      CALL require(es == 0, 'seglen:weighted-assemble-fresh')
      CALL CD_Init_Model(mw0, q0, v0, conn, l0, ea, rho_a, .FALSE., fw0, fixed, cfg, es, em, &
                         dist_load_per_length=wpl)
      CALL require(es == CD_MODEL_OK, 'seglen:init-weighted-base')
      CALL CD_Init_Model(mw1, q0, v0, conn, l0b, ea, rho_a, .FALSE., fw1, fixed, cfg, es, em, &
                         dist_load_per_length=wpl)
      CALL require(es == CD_MODEL_OK, 'seglen:init-weighted-fresh')
      CALL CD_Update_Model_SegmentLength(mw0, NE, l0b(NE), 0.0_wp, es, em)
      CALL require(es == CD_MODEL_OK, 'seglen:weighted-update: '//TRIM(em))
      CALL require(nan_max_abs(mw0%f_ext - mw1%f_ext) <= 0.0_wp, 'seglen:weighted-fext-bitwise')
      CALL require(nan_max_abs(mw0%l0 - mw1%l0) <= 0.0_wp, 'seglen:weighted-l0-bitwise')
      CALL require(nan_max_abs(mw0%a - mw1%a) <= 0.0_wp, 'seglen:weighted-accel-bitwise')
      CALL CD_Calc_Model_CoupledLoads(mw0, lw_a, es, em)
      CALL CD_Calc_Model_CoupledLoads(mw1, lw_b, es, em)
      CALL require(nan_max_abs(lw_a - lw_b) <= 0.0_wp, 'seglen:weighted-coupled-loads-equal-fresh-build')
      CALL CD_Step_Model(mw0, 0.02_wp, conv, stalled, n_iter, es, em)
      CALL require(es == CD_MODEL_OK, 'seglen:weighted-step-0')
      CALL CD_Step_Model(mw1, 0.02_wp, conv, stalled, n_iter, es, em)
      CALL require(es == CD_MODEL_OK, 'seglen:weighted-step-1')
      CALL CD_Get_Model_State(mw0, qw0, vw0, aw0, es, em)
      CALL CD_Get_Model_State(mw1, qw1, vw1, aw1, es, em)
      CALL require(nan_max_abs(qw0 - qw1) <= 0.0_wp .AND. nan_max_abs(vw0 - vw1) <= 0.0_wp .AND. &
                   nan_max_abs(aw0 - aw1) <= 0.0_wp, 'seglen:weighted-step-equals-fresh-build-bitwise')
      CALL CD_End_Model(mw0, es, em)
      CALL CD_End_Model(mw1, es, em)
    END BLOCK

    ! Runtime EXTERNAL-LOAD updates must survive a later control command: the
    ! segment-length apply reconstructs f_ext from the distributed-load remainder,
    ! so the external-load path must maintain that remainder -- pre-fix a control
    ! command silently WIPED any post-init external load. Bitwise against a fresh
    ! build at the new length carrying the same extra point load.
    BLOCK
      TYPE(CD_ModelType) :: mx0, mx1
      REAL(wp) :: wplx(3, NE), fx0(NDOF), fx1(NDOF), lx_a(6), lx_b(6)
      wplx = 0.0_wp
      wplx(3, :) = -7.5_wp
      CALL CD_Assemble_Distributed_Load(conn, l0, wplx, fx0, es, em)
      CALL require(es == 0, 'seglen:extload-assemble-base')
      CALL CD_Init_Model(mx0, q0, v0, conn, l0, ea, rho_a, .FALSE., fx0, fixed, cfg, es, em, &
                         dist_load_per_length=wplx)
      CALL require(es == CD_MODEL_OK, 'seglen:extload-init-base')
      ! runtime external-load update: a point load bump at an interior node
      fx0(3*2 - 1) = fx0(3*2 - 1) + 3.0_wp
      CALL CD_Update_Model_External_Loads(mx0, fx0, es, em)
      CALL require(es == CD_MODEL_OK, 'seglen:extload-update: '//TRIM(em))
      CALL CD_Update_Model_SegmentLength(mx0, NE, l0b(NE), 0.0_wp, es, em)
      CALL require(es == CD_MODEL_OK, 'seglen:extload-control: '//TRIM(em))
      ! fresh twin: assembled at the NEW length plus the same bump
      CALL CD_Assemble_Distributed_Load(conn, l0b, wplx, fx1, es, em)
      CALL require(es == 0, 'seglen:extload-assemble-fresh')
      fx1(3*2 - 1) = fx1(3*2 - 1) + 3.0_wp
      CALL CD_Init_Model(mx1, q0, v0, conn, l0b, ea, rho_a, .FALSE., fx1, fixed, cfg, es, em, &
                         dist_load_per_length=wplx)
      CALL require(es == CD_MODEL_OK, 'seglen:extload-init-fresh')
      CALL require(nan_max_abs(mx0%f_ext - mx1%f_ext) <= 0.0_wp, 'seglen:extload-survives-control-bitwise')
      CALL CD_Calc_Model_CoupledLoads(mx0, lx_a, es, em)
      CALL CD_Calc_Model_CoupledLoads(mx1, lx_b, es, em)
      CALL require(nan_max_abs(lx_a - lx_b) <= 0.0_wp, 'seglen:extload-coupled-loads-equal-fresh')
      CALL CD_End_Model(mx0, es, em)
      CALL CD_End_Model(mx1, es, em)
    END BLOCK

    ! The payout RATE must reach the PRODUCTION (banded) damping path: two
    ! ba-damped twins, one with l0_dot = 0.3 at the SAME length, must diverge
    ! within one step -- a regression against the banded element kernel ignoring
    ! l0_dot (which left the twins bit-equal).
    BLOCK
      TYPE(CD_ModelType) :: md0, md1
      REAL(wp) :: ba(NE), q0d(NDOF), v0d(NDOF)
      REAL(wp) :: qd0(NDOF), vd0(NDOF), ad0(NDOF), qd1(NDOF), vd1(NDOF), ad1(NDOF)
      ba = 5.0_wp
      q0d = q0
      v0d = 0.0_wp
      v0d(3*2 - 2) = 0.1_wp   ! axial motion so the damping term is live
      CALL CD_Init_Model(md0, q0d, v0d, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, ba=ba)
      CALL require(es == CD_MODEL_OK, 'seglen:init-damped-0')
      CALL CD_Init_Model(md1, q0d, v0d, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, ba=ba)
      CALL require(es == CD_MODEL_OK, 'seglen:init-damped-1')
      CALL CD_Update_Model_SegmentLength(md1, NE, l0(NE), 0.3_wp, es, em)
      CALL require(es == CD_MODEL_OK, 'seglen:rate-only-update')
      CALL CD_Step_Model(md0, 0.02_wp, conv, stalled, n_iter, es, em)
      CALL require(es == CD_MODEL_OK, 'seglen:damped-step-0')
      CALL CD_Step_Model(md1, 0.02_wp, conv, stalled, n_iter, es, em)
      CALL require(es == CD_MODEL_OK, 'seglen:damped-step-1')
      CALL CD_Get_Model_State(md0, qd0, vd0, ad0, es, em)
      CALL CD_Get_Model_State(md1, qd1, vd1, ad1, es, em)
      CALL require(nan_max_abs(qd1 - qd0) + nan_max_abs(vd1 - vd0) > 0.0_wp, &
                   'seglen:payout-rate-reaches-the-production-step')
      CALL CD_End_Model(md0, es, em)
      CALL CD_End_Model(md1, es, em)
    END BLOCK

    ! Seabed-DAMPING models without the cn/kn ratio must NOT enable the length
    ! recompute (a zero-default ratio would silently wipe cn on the first update);
    ! the update fails closed by name instead.
    BLOCK
      TYPE(CD_ModelType) :: ms
      REAL(wp) :: kn_s(NN), cn_s(NN), diam_s(NE)
      kn_s = 1.0e4_wp
      cn_s = 1.0e3_wp
      diam_s = 0.1_wp
      CALL CD_Init_Model(ms, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                         seabed_z_floor=-5.0_wp, seabed_kn=kn_s, seabed_cn=cn_s, &
                         seabed_kbot=1.0e5_wp, seabed_contact_diameter=diam_s)
      CALL require(es == CD_MODEL_OK, 'seglen:init-seabed-no-ratio')
      CALL CD_Update_Model_SegmentLength(ms, NE, 1.0_wp, 0.0_wp, es, em)
      CALL require(es /= CD_MODEL_OK, 'seglen:seabed-damping-without-ratio-fails-closed')
      CALL CD_End_Model(ms, es, em)
      ! REUSE: a model that previously carried recompute state, re-initialized as a
      ! seabed model WITHOUT the inputs, must fail the update closed -- not inherit
      ! the stale flag and read a deallocated diameter array.
      CALL CD_Init_Model(ms, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                         seabed_z_floor=-5.0_wp, seabed_kn=kn_s, seabed_cn=cn_s, &
                         seabed_kbot=1.0e5_wp, seabed_contact_diameter=diam_s, seabed_cn_over_kn=0.1_wp)
      CALL require(es == CD_MODEL_OK, 'seglen:init-recompute-enabled')
      CALL CD_Update_Model_SegmentLength(ms, NE, 1.0_wp, 0.0_wp, es, em)
      CALL require(es == CD_MODEL_OK, 'seglen:recompute-enabled-update-ok: '//TRIM(em))
      CALL CD_End_Model(ms, es, em)
      CALL CD_Init_Model(ms, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                         seabed_z_floor=-5.0_wp, seabed_kn=kn_s)
      CALL require(es == CD_MODEL_OK, 'seglen:reuse-init-without-recompute')
      CALL CD_Update_Model_SegmentLength(ms, NE, 1.0_wp, 0.0_wp, es, em)
      CALL require(es /= CD_MODEL_OK, 'seglen:reused-model-without-inputs-fails-closed')
      CALL CD_End_Model(ms, es, em)
    END BLOCK

    ! fail-closed battery
    CALL CD_Update_Model_SegmentLength(ma, 0, 1.0_wp, 0.0_wp, es, em)
    CALL require(es /= CD_MODEL_OK, 'seglen:bad-elem-fails')
    CALL CD_Update_Model_SegmentLength(ma, NE, 0.0_wp, 0.0_wp, es, em)
    CALL require(es /= CD_MODEL_OK, 'seglen:zero-l0-fails')
    ! deliberately non-finite or overflowing input: must not halt a trapping build
    CALL IEEE_GET_HALTING_MODE(IEEE_USUAL, fp_halt)
    CALL IEEE_SET_HALTING_MODE(IEEE_USUAL, .FALSE.)
    CALL CD_Update_Model_SegmentLength(ma, NE, IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN), 0.0_wp, es, em)
    CALL IEEE_SET_FLAG(IEEE_USUAL, .FALSE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_USUAL, fp_halt)
    CALL require(es /= CD_MODEL_OK, 'seglen:nan-l0-fails')

    CALL CD_End_Model(ma, es, em)
    CALL CD_End_Model(mb, es, em)
    CALL CD_Update_Model_SegmentLength(ma, 1, 1.0_wp, 0.0_wp, es, em)
    CALL require(es /= CD_MODEL_OK, 'seglen:uninitialized-fails')
  END SUBROUTINE case_segment_length_update

  SUBROUTINE case_copy_model_preserves_state()
    !! Checked model copies preserve persistent state and leave transient solver
    !! workspaces independent of the source object.
    INTEGER, PARAMETER :: NE = 3, NN = NE + 1, NDOF = 3*NN
    TYPE(CD_ModelType) :: src, dst
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(6), es
    REAL(wp) :: l0(NE), ea(NE), rho_a(NE), q0(NDOF), v0(NDOF), f_ext(NDOF)
    REAL(wp) :: q_src(NDOF), v_src(NDOF), a_src(NDOF), q_dst(NDOF), v_dst(NDOF), a_dst(NDOF)
    REAL(wp) :: load_src(6), load_dst(6)
    CHARACTER(200) :: em

    CALL straight_cable(NE, conn, l0, ea, rho_a, q0, fixed)
    v0 = 0.0_wp
    f_ext = 0.0_wp
    f_ext(5) = -25.0_wp

    CALL CD_Init_Model(src, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK .AND. CD_Model_Is_Initialized(src), 'copy:init-src: '//TRIM(em))
    CALL CD_Calc_Model_CoupledLoads(src, load_src, es, em)
    CALL require(es == CD_MODEL_OK .AND. ALLOCATED(src%dynamic_workspace%R), 'copy:source-loads: '//TRIM(em))

    CALL CD_Copy_Model(src, dst, es, em)
    CALL require(es == CD_MODEL_OK .AND. CD_Model_Is_Initialized(dst), 'copy:model: '//TRIM(em))
    CALL require(.NOT. ALLOCATED(dst%dynamic_workspace%R), 'copy:workspace-not-shared')
    CALL CD_Get_Model_State(src, q_src, v_src, a_src, es, em)
    CALL require(es == CD_MODEL_OK, 'copy:get-src: '//TRIM(em))
    CALL CD_Get_Model_State(dst, q_dst, v_dst, a_dst, es, em)
    CALL require(es == CD_MODEL_OK, 'copy:get-dst: '//TRIM(em))
    CALL require(nan_max_abs(q_dst - q_src) < 1.0e-12_wp .AND. &
                 nan_max_abs(v_dst - v_src) < 1.0e-12_wp .AND. &
                 nan_max_abs(a_dst - a_src) < 1.0e-12_wp, 'copy:state')
    CALL CD_Calc_Model_CoupledLoads(dst, load_dst, es, em)
    CALL require(es == CD_MODEL_OK, 'copy:dst-loads: '//TRIM(em))
    CALL require(nan_max_abs(load_dst - load_src) < 1.0e-12_wp, 'copy:loads')

    ! Copy again into the now-initialized dst: the commit path must release the previous
    ! model and reproduce the source state, not leak or corrupt the already-initialized target.
    CALL CD_Copy_Model(src, dst, es, em)
    CALL require(es == CD_MODEL_OK .AND. CD_Model_Is_Initialized(dst), 'copy:recopy-into-initialized: '//TRIM(em))
    CALL CD_Get_Model_State(dst, q_dst, v_dst, a_dst, es, em)
    CALL require(es == CD_MODEL_OK .AND. nan_max_abs(q_dst - q_src) < 1.0e-12_wp .AND. &
                 nan_max_abs(a_dst - a_src) < 1.0e-12_wp, 'copy:recopy-state')
    CALL CD_Calc_Model_CoupledLoads(dst, load_dst, es, em)
    CALL require(es == CD_MODEL_OK .AND. nan_max_abs(load_dst - load_src) < 1.0e-12_wp, 'copy:recopy-loads')

    CALL CD_End_Model(dst, es, em)
    CALL require(es == CD_MODEL_OK, 'copy:end-dst')
    CALL CD_End_Model(src, es, em)
    CALL require(es == CD_MODEL_OK, 'copy:end-src')
  END SUBROUTINE case_copy_model_preserves_state

  SUBROUTINE case_prescribed_step_updates_state()
    !! Prescribed endpoint motion passes through the model boundary and persists in
    !! the stored state after the step.
    INTEGER, PARAMETER :: NE = 4, NN = NE + 1, NDOF = 3*NN
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(6), es, n_iter
    REAL(wp) :: l0(NE), ea(NE), rho_a(NE), q0(NDOF), v0(NDOF), f_ext(NDOF)
    REAL(wp) :: pq(NDOF), pv(NDOF), pa(NDOF), q(NDOF), v(NDOF), a(NDOF)
    LOGICAL :: conv, stalled
    CHARACTER(200) :: em

    CALL straight_cable(NE, conn, l0, ea, rho_a, q0, fixed)
    v0 = 0.0_wp
    f_ext = 0.0_wp
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK, 'prescribed:init')

    pq = q0
    pv = 0.0_wp
    pa = 0.0_wp
    pq(NDOF) = q0(NDOF) + 0.02_wp
    pv(NDOF) = 1.0_wp
    CALL CD_Step_Model(model, 0.02_wp, conv, stalled, n_iter, es, em, &
                       prescribed_q=pq, prescribed_v=pv, prescribed_a=pa)
    CALL require(es == CD_MODEL_OK .AND. conv, 'prescribed:step')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK, 'prescribed:get')
    CALL require(ABS(q(NDOF) - pq(NDOF)) < 1.0e-12_wp, 'prescribed:q-stored')
    CALL require(ABS(v(NDOF) - pv(NDOF)) < 1.0e-12_wp, 'prescribed:v-stored')
    CALL require(ABS(a(NDOF) - pa(NDOF)) < 1.0e-12_wp, 'prescribed:a-stored')
    CALL require(nan_max_abs(v(4:NDOF - 3)) > 1.0e-3_wp, 'prescribed:interior-excited')

    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'prescribed:end')
  END SUBROUTINE case_prescribed_step_updates_state

  SUBROUTINE case_recovery_common_path()
    !! A nominal interval that converges directly remains bit-identical through
    !! the recovery entry point.  This protects the ordinary production path
    !! while the integration cases exercise the rare subdivision branch.
    INTEGER, PARAMETER :: NE = 4, NN = NE + 1, NDOF = 3*NN
    TYPE(CD_ModelType) :: plain, recovering
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(6), es, n_iter, substeps_used
    REAL(wp) :: l0(NE), ea(NE), rho_a(NE), q0(NDOF), v0(NDOF), f_ext(NDOF)
    REAL(wp) :: pq(NDOF), pv(NDOF), pa(NDOF)
    LOGICAL :: conv, stalled
    CHARACTER(240) :: em

    CALL straight_cable(NE, conn, l0, ea, rho_a, q0, fixed)
    v0 = 0.0_wp
    f_ext = 0.0_wp
    CALL CD_Init_Model(plain, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK, 'recovery-common:init: '//TRIM(em))
    CALL CD_Copy_Model(plain, recovering, es, em)
    CALL require(es == CD_MODEL_OK, 'recovery-common:copy: '//TRIM(em))
    ! Recovery owns a dedicated interval snapshot and must not overwrite the
    ! rollback buffers reserved for enclosing multi-line transactions.
    recovering%q_rollback = 101.0_wp
    recovering%v_rollback = 102.0_wp
    recovering%a_rollback = 103.0_wp

    pq = q0
    pv = 0.0_wp
    pa = 0.0_wp
    pq(NDOF) = q0(NDOF) + 0.02_wp
    pv(NDOF) = 1.0_wp
    CALL CD_Step_Model(plain, 0.02_wp, conv, stalled, n_iter, es, em, &
                       prescribed_q=pq, prescribed_v=pv, prescribed_a=pa)
    CALL require(es == CD_MODEL_OK .AND. conv .AND. .NOT. stalled, 'recovery-common:plain: '//TRIM(em))
    CALL CD_Step_Model_Recovering(recovering, 0.02_wp, conv, stalled, n_iter, es, em, &
                                  prescribed_q=pq, prescribed_v=pv, prescribed_a=pa, &
                                  substeps_used=substeps_used)
    CALL require(es == CD_MODEL_OK .AND. conv .AND. .NOT. stalled, &
                 'recovery-common:recovering: '//TRIM(em))
    CALL require(nan_max_abs(plain%q - recovering%q) <= 0.0_wp .AND. &
                 nan_max_abs(plain%v - recovering%v) <= 0.0_wp .AND. &
                 nan_max_abs(plain%a - recovering%a) <= 0.0_wp, &
                 'recovery-common:bit-identical')
    CALL require(substeps_used == 1, 'recovery-common:reports-direct-step')
    CALL require(nan_max_abs(recovering%q_rollback - 101.0_wp) <= 0.0_wp .AND. &
                 nan_max_abs(recovering%v_rollback - 102.0_wp) <= 0.0_wp .AND. &
                 nan_max_abs(recovering%a_rollback - 103.0_wp) <= 0.0_wp, &
                 'recovery-common:enclosing-rollback-workspace-preserved')

    CALL CD_End_Model(recovering, es, em)
    CALL CD_End_Model(plain, es, em)
  END SUBROUTINE case_recovery_common_path

  SUBROUTINE case_recovery_subdivision_path()
    !! Find a bounded, deterministic prescribed-motion increment for which the
    !! deliberately small nonlinear-iteration allowance rejects the full
    !! interval but accepts the same C2 trajectory after subdivision.  This
    !! exercises the rare recovery branch without relying on a long simulation or
    !! compiler-sensitive near-stall time history.
    INTEGER, PARAMETER :: NE = 4, NN = NE + 1, NDOF = 3*NN
    REAL(wp), PARAMETER :: dts(4) = [0.02_wp, 0.05_wp, 0.10_wp, 0.20_wp]
    REAL(wp), PARAMETER :: offsets(6) = [0.02_wp, 0.05_wp, 0.10_wp, 0.20_wp, 0.35_wp, 0.50_wp]
    TYPE(CD_ModelType) :: plain, recovering
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(6), es, n_iter, substeps_used, idt, idz
    REAL(wp) :: l0(NE), ea(NE), rho_a(NE), q0(NDOF), v0(NDOF), f_ext(NDOF)
    REAL(wp) :: pq(NDOF), pv(NDOF), pa(NDOF), q(NDOF), v(NDOF), a(NDOF)
    REAL(wp) :: found_dt, found_offset
    INTEGER :: found_substeps
    LOGICAL :: conv, stalled, found, direct_soft_failure
    CHARACTER(240) :: em

    CALL straight_cable(NE, conn, l0, ea, rho_a, q0, fixed)
    v0 = 0.0_wp
    f_ext = 0.0_wp
    cfg%max_iter = 3
    found = .FALSE.
    found_dt = 0.0_wp; found_offset = 0.0_wp; found_substeps = 0
    DO idt = 1, SIZE(dts)
      DO idz = 1, SIZE(offsets)
        CALL CD_Init_Model(plain, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
        CALL require(es == CD_MODEL_OK, 'recovery-subdivision:init: '//TRIM(em))
        IF (es /= CD_MODEL_OK) RETURN
        CALL CD_Copy_Model(plain, recovering, es, em)
        CALL require(es == CD_MODEL_OK, 'recovery-subdivision:copy: '//TRIM(em))
        IF (es /= CD_MODEL_OK) RETURN

        pq = q0
        pv = 0.0_wp
        pa = 0.0_wp
        pq(NDOF) = q0(NDOF) + offsets(idz)
        CALL CD_Step_Model(plain, dts(idt), conv, stalled, n_iter, es, em, &
                           prescribed_q=pq, prescribed_v=pv, prescribed_a=pa)
        direct_soft_failure = es == CD_MODEL_OK .AND. .NOT. conv
        IF (direct_soft_failure) THEN
          CALL CD_Step_Model_Recovering(recovering, dts(idt), conv, stalled, n_iter, es, em, &
                                        prescribed_q=pq, prescribed_v=pv, prescribed_a=pa, &
                                        max_substeps=256, substeps_used=substeps_used)
          IF (es == CD_MODEL_OK .AND. conv .AND. .NOT. stalled .AND. substeps_used > 1) THEN
            found = .TRUE.
            found_dt = dts(idt); found_offset = offsets(idz); found_substeps = substeps_used
            CALL CD_Get_Model_State(recovering, q, v, a, es, em)
            CALL require(es == CD_MODEL_OK .AND. ALL(IEEE_IS_FINITE(q)) .AND. &
                         ALL(IEEE_IS_FINITE(v)) .AND. ALL(IEEE_IS_FINITE(a)), &
                         'recovery-subdivision:finite committed state')
            CALL require(ABS(q(NDOF) - pq(NDOF)) <= 0.0_wp .AND. &
                         ABS(v(NDOF) - pv(NDOF)) <= 0.0_wp .AND. &
                         ABS(a(NDOF) - pa(NDOF)) <= 0.0_wp, &
                         'recovery-subdivision:lands on prescribed endpoint')
          END IF
        END IF
        CALL CD_End_Model(recovering, es, em)
        CALL CD_End_Model(plain, es, em)
        IF (found) EXIT
      END DO
      IF (found) EXIT
    END DO
    CALL require(found, 'recovery-subdivision:soft full-step rejection recovered by C2 subdivision')
    IF (found) WRITE (*, '(A,F7.4,A,F7.4,A,I0)') '  recovery branch gate: dt=', found_dt, &
      ' s, endpoint offset=', found_offset, ' m, substeps=', found_substeps
  END SUBROUTINE case_recovery_subdivision_path

  SUBROUTINE case_coupled_load_sign()
    !! Coupled loads are the force exerted by the cable on the coupled object. For a
    !! one-element line stretched from x=0 to x=1.2 with L0=1 and EA=100, T=20:
    !! the cable pulls the left support in +x and the right support in -x.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), loads(6), loads_bad(5)
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, 0.0_wp, 1.2_wp, 0.0_wp, 0.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [5.0_wp]
    f_ext = 0.0_wp

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .TRUE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK, 'loads:init')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'loads:calc')
    CALL require(ABS(loads(1) - 20.0_wp) < 1.0e-12_wp, 'loads:left-x')
    CALL require(ABS(loads(4) + 20.0_wp) < 1.0e-12_wp, 'loads:right-x')
    CALL require(nan_max_abs(loads([2, 3, 5, 6])) < 1.0e-12_wp, 'loads:transverse-zero')
    CALL CD_Calc_Model_CoupledLoads(model, loads_bad, es, em)
    CALL require(es == CD_MODEL_BADINPUT, 'loads:bad-shape')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'loads:end')
  END SUBROUTINE case_coupled_load_sign

  SUBROUTINE case_coupled_motion_update()
    !! Compact coupled-DOF kinematics updates use the same ordering as
    !! CD_Get_Model_CoupledDofs and are reflected by coupled-load extraction.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), loads(6)
    REAL(wp) :: q_c(6), v_c(6), a_c(6), q_out(6), v_out(6), a_out(6), bad(5)
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [5.0_wp]
    f_ext = 0.0_wp
    q_c = q0
    v_c = 0.0_wp
    a_c = 0.0_wp

    CALL CD_Update_Model_CoupledMotion(model, q_c, v_c, a_c, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED, 'motion-update:uninitialized')
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .TRUE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK, 'motion-update:init')
    q_c(4) = 1.2_wp
    CALL CD_Update_Model_CoupledMotion(model, q_c, v_c, a_c, es, em)
    CALL require(es == CD_MODEL_OK, 'motion-update:set')
    CALL CD_Get_Model_CoupledMotion(model, q_out, v_out, a_out, es, em)
    CALL require(es == CD_MODEL_OK, 'motion-update:get')
    CALL require(nan_max_abs(q_out - q_c) < 1.0e-12_wp, 'motion-update:q-roundtrip')
    CALL require(nan_max_abs(v_out - v_c) < 1.0e-12_wp, 'motion-update:v-roundtrip')
    CALL require(nan_max_abs(a_out - a_c) < 1.0e-12_wp, 'motion-update:a-roundtrip')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'motion-update:loads')
    CALL require(ABS(loads(1) - 20.0_wp) < 1.0e-12_wp, 'motion-update:left-x')
    CALL require(ABS(loads(4) + 20.0_wp) < 1.0e-12_wp, 'motion-update:right-x')
    bad = 0.0_wp
    CALL CD_Get_Model_CoupledMotion(model, bad, v_c, a_c, es, em)
    CALL require(es == CD_MODEL_BADINPUT, 'motion-update:get-bad-shape')
    CALL CD_Update_Model_CoupledMotion(model, bad, v_c, a_c, es, em)
    CALL require(es == CD_MODEL_BADINPUT, 'motion-update:bad-shape')
    q_c(2) = IEEE_VALUE(q_c(2), IEEE_QUIET_NAN)
    CALL CD_Update_Model_CoupledMotion(model, q_c, v_c, a_c, es, em)
    CALL require(es == CD_MODEL_BADINPUT, 'motion-update:nonfinite')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'motion-update:end')
  END SUBROUTINE case_coupled_motion_update

  SUBROUTINE case_external_load_update()
    !! Runtime force updates are reflected immediately by the coupled-load
    !! boundary. With an unstretched fixed one-element line, coupled loads equal
    !! the model-owned external nodal loads.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), loads(6), bad(5)
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [5.0_wp]
    f_ext = 0.0_wp

    CALL CD_Update_Model_External_Loads(model, f_ext, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED, 'force-update:uninitialized')
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK, 'force-update:init')
    f_ext = [12.0_wp, 0.0_wp, 0.0_wp, -3.0_wp, 4.0_wp, 0.0_wp]
    CALL CD_Update_Model_External_Loads(model, f_ext, es, em)
    CALL require(es == CD_MODEL_OK, 'force-update:set')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'force-update:loads')
    CALL require(nan_max_abs(loads - f_ext) < 1.0e-12_wp, 'force-update:loads-match')
    bad = 0.0_wp
    CALL CD_Update_Model_External_Loads(model, bad, es, em)
    CALL require(es == CD_MODEL_BADINPUT, 'force-update:bad-shape')
    f_ext(2) = IEEE_VALUE(f_ext(2), IEEE_QUIET_NAN)
    CALL CD_Update_Model_External_Loads(model, f_ext, es, em)
    CALL require(es == CD_MODEL_BADINPUT, 'force-update:nonfinite')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'force-update:end')
  END SUBROUTINE case_external_load_update

  SUBROUTINE case_recompute_acceleration_after_load_update()
    !! Runtime load changes immediately refresh acceleration. The explicit
    !! recompute API remains idempotent for callers that request it directly.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(3), es
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), q(6), v(6), a(6)
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3]
    q0 = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [3.0_wp]   ! free-node consistent mass = rho_a*L/3 = 1
    f_ext = 0.0_wp

    CALL CD_Recompute_Model_Acceleration(model, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED, 'accel-recompute:uninitialized')
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK, 'accel-recompute:init')
    f_ext(4) = 3.0_wp
    CALL CD_Update_Model_External_Loads(model, f_ext, es, em)
    CALL require(es == CD_MODEL_OK, 'accel-recompute:update-load')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK .AND. ABS(a(4) - 3.0_wp) < 1.0e-12_wp, 'accel-recompute:after-update')
    CALL CD_Recompute_Model_Acceleration(model, es, em)
    CALL require(es == CD_MODEL_OK, 'accel-recompute:call')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK .AND. ABS(a(4) - 3.0_wp) < 1.0e-12_wp, 'accel-recompute:after')
    f_ext = 0.0_wp
    CALL CD_Update_Model_External_Loads(model, f_ext, es, em)
    CALL require(es == CD_MODEL_OK, 'accel-recompute:clear-load')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    a(1:3) = [2.0_wp, 0.0_wp, 0.0_wp]
    CALL CD_Update_Model_CoupledMotion(model, q(1:3), v(1:3), a(1:3), es, em)
    CALL require(es == CD_MODEL_OK, 'accel-recompute:set-support-accel')
    CALL CD_Recompute_Model_Acceleration(model, es, em)
    CALL require(es == CD_MODEL_OK, 'accel-recompute:support-call')
    CALL require(ALLOCATED(model%dynamic_workspace%M) .AND. ALLOCATED(model%dynamic_workspace%R_free), &
                 'accel-recompute:workspace-owned')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK .AND. ABS(a(1) - 2.0_wp) < 1.0e-12_wp .AND. &
                 ABS(a(4) + 1.0_wp) < 1.0e-12_wp, 'accel-recompute:support-inertia')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'accel-recompute:end')
    CALL require(.NOT. ALLOCATED(model%dynamic_workspace%M) .AND. &
                 .NOT. ALLOCATED(model%dynamic_workspace%R_free), 'accel-recompute:workspace-released')
  END SUBROUTINE case_recompute_acceleration_after_load_update

  SUBROUTINE case_seabed_coupled_load()
    !! A model-owned seabed contributor participates in coupled-load extraction.
    !! The element is at rest length, node 2 is 0.1 m below the floor with k_n=1000,
    !! so the upward contact load passed through the model boundary is +100 N.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), kn(2), loads(6)
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, -0.1_wp]
    v0 = 0.0_wp
    l0 = [SQRT(1.01_wp)]
    ea = [100.0_wp]
    rho_a = [5.0_wp]
    f_ext = 0.0_wp
    kn = [1000.0_wp, 1000.0_wp]

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       seabed_z_floor=0.0_wp, seabed_kn=kn)
    CALL require(es == CD_MODEL_OK, 'seabed-load:init')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-load:calc')
    CALL require(ABS(loads(6) - 100.0_wp) < 1.0e-10_wp, 'seabed-load:node2-z')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-load:end')
  END SUBROUTINE case_seabed_coupled_load

  SUBROUTINE case_seabed_normal_damping_load()
    !! Contact normal damping contributes only while a contacted node moves downward:
    !! Fz = kn*(z_floor - z) - cn*vz for vz < 0.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp) :: q0(6), v0(6), f_ext(6), l0(1), ea(1), rho_a(1), kn(2), cn(2), loads(6)
    CHARACTER(160) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, -0.1_wp, 1.0_wp, 0.0_wp, -0.1_wp]
    v0 = 0.0_wp
    v0(3) = -2.0_wp
    v0(6) = 0.1_wp
    f_ext = 0.0_wp
    l0 = 1.0_wp
    ea = 1.0e6_wp
    rho_a = 10.0_wp
    kn = 100.0_wp
    cn = 100.0_wp

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       seabed_z_floor=0.0_wp, seabed_kn=kn, seabed_cn=cn)
    CALL require(es == CD_MODEL_OK, 'seabed-damp:init')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-damp:loads')
    CALL require(ABS(loads(3) - 210.0_wp) < 1.0e-10_wp, 'seabed-damp:downward-node')
    CALL require(ABS(loads(6) - 10.0_wp) < 1.0e-10_wp, 'seabed-damp:resting-node')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-damp:end')
  END SUBROUTINE case_seabed_normal_damping_load

  SUBROUTINE case_seabed_normal_damping_tangent_fd()
    !! A downward-moving node in the C1 touchdown blend has a continuous normal
    !! damping force and a position tangent. The reduced coupled-load derivative
    !! must match finite differences so the dynamic Newton path does not see a
    !! hard damping active-set switch at touchdown.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es, j
    REAL(wp) :: q0(6), v0(6), f_ext(6), l0(1), ea(1), rho_a(1), kn(2), cn(2)
    REAL(wp) :: q_c(6), v_c(6), a_c(6), q_p(6), q_m(6), v_p(6), v_m(6)
    REAL(wp) :: loads_p(6), loads_m(6), jq(6, 6), jv(6, 6), fdq(6, 6), fdv(6, 6)
    REAL(wp) :: hq, hv, err_q, err_v
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, -0.5_wp*CD_SEABED_CONTACT_BLEND, &
          1.0_wp, 0.0_wp, -0.1_wp]
    v0 = 0.0_wp
    v0(3) = -2.0_wp
    v0(6) = 0.1_wp
    f_ext = 0.0_wp
    l0 = [SQRT((q0(4) - q0(1))**2 + (q0(6) - q0(3))**2)]
    ea = [1.0e3_wp]
    rho_a = [10.0_wp]
    kn = [100.0_wp, 100.0_wp]
    cn = [100.0_wp, 100.0_wp]
    hq = 1.0e-8_wp
    hv = 1.0e-6_wp

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       seabed_z_floor=0.0_wp, seabed_kn=kn, seabed_cn=cn)
    CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:init: '//TRIM(em))
    CALL CD_Get_Model_CoupledMotion(model, q_c, v_c, a_c, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:get-motion')
    CALL CD_Calc_Model_CoupledKinematicDerivatives(model, jq, jv, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:analytic: '//TRIM(em))

    DO j = 1, 6
      q_p = q_c
      q_m = q_c
      q_p(j) = q_p(j) + hq
      q_m(j) = q_m(j) - hq
      CALL CD_Update_Model_CoupledMotion(model, q_p, v_c, a_c, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:q-plus')
      CALL CD_Calc_Model_CoupledLoads(model, loads_p, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:q-plus-loads')
      CALL CD_Update_Model_CoupledMotion(model, q_m, v_c, a_c, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:q-minus')
      CALL CD_Calc_Model_CoupledLoads(model, loads_m, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:q-minus-loads')
      fdq(:, j) = (loads_p - loads_m)/(2.0_wp*hq)

      v_p = v_c
      v_m = v_c
      v_p(j) = v_p(j) + hv
      v_m(j) = v_m(j) - hv
      CALL CD_Update_Model_CoupledMotion(model, q_c, v_p, a_c, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:v-plus')
      CALL CD_Calc_Model_CoupledLoads(model, loads_p, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:v-plus-loads')
      CALL CD_Update_Model_CoupledMotion(model, q_c, v_m, a_c, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:v-minus')
      CALL CD_Calc_Model_CoupledLoads(model, loads_m, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:v-minus-loads')
      fdv(:, j) = (loads_p - loads_m)/(2.0_wp*hv)
    END DO

    err_q = nan_max_abs(fdq - jq)
    err_v = nan_max_abs(fdv - jv)
    IF (err_v >= 2.0e-6_wp) WRITE (*, '(A,ES12.4)') 'seabed-damp-fd max dload/dv error = ', err_v
    CALL require(err_q < 2.0e-2_wp, 'seabed-damp-fd:dload-dq')
    CALL require(err_v < 2.0e-6_wp, 'seabed-damp-fd:dload-dv')
    CALL CD_Update_Model_CoupledMotion(model, q_c, v_c, a_c, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:restore')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-damp-fd:end')
  END SUBROUTINE case_seabed_normal_damping_tangent_fd

  SUBROUTINE case_seabed_friction_load()
    !! Stick-slip friction springs: anchored at the initial positions they carry no force;
    !! stretched below the capacity mu*N (N = spring penetration + one-sided normal damper)
    !! the load is -k d, beyond it -mu*N d/|d|.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp) :: q0(6), v0(6), f_ext(6), l0(1), ea(1), rho_a(1), kn(2), cn(2), loads(6), anchors(2, 2)
    CHARACTER(160) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, -0.1_wp, 1.0_wp, 0.0_wp, -0.1_wp]
    v0 = 0.0_wp
    v0(1) = 3.0_wp
    v0(2) = 4.0_wp
    v0(3) = -2.0_wp
    f_ext = 0.0_wp
    l0 = 1.0_wp
    ea = 1.0e6_wp
    rho_a = 10.0_wp
    kn = 100.0_wp
    cn = 100.0_wp

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       seabed_z_floor=0.0_wp, seabed_kn=kn, seabed_cn=cn, seabed_mu=0.5_wp)
    CALL require(es == CD_MODEL_OK, 'seabed-fric:init')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-fric:loads')
    CALL require(ABS(loads(1)) + ABS(loads(2)) <= 0.0_wp, 'seabed-fric:unstretched')
    CALL require(ABS(loads(3) - 210.0_wp) < 1.0e-10_wp, 'seabed-fric:normal')
    ! stick: |k d| = 100*0.5 = 50 below mu*N = 105
    CALL CD_Get_Model_Friction_Anchors(model, anchors, es, em)
    CALL require(es == CD_MODEL_OK .AND. nan_max_abs(anchors(:, 1) - q0(1:2)) <= 0.0_wp, 'seabed-fric:anchor0')
    anchors(:, 1) = [-0.3_wp, -0.4_wp]
    CALL CD_Set_Model_Friction_Anchors(model, anchors, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-fric:set-stick')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(ABS(loads(1) + 30.0_wp) < 1.0e-10_wp, 'seabed-fric:stick-x')
    CALL require(ABS(loads(2) + 40.0_wp) < 1.0e-10_wp, 'seabed-fric:stick-y')
    ! slip: |k d| = 500 above mu*N = 105
    anchors(:, 1) = [-3.0_wp, -4.0_wp]
    CALL CD_Set_Model_Friction_Anchors(model, anchors, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-fric:set-slip')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(ABS(loads(1) + 0.5_wp*210.0_wp*0.6_wp) < 1.0e-10_wp, 'seabed-fric:slip-x')
    CALL require(ABS(loads(2) + 0.5_wp*210.0_wp*0.8_wp) < 1.0e-10_wp, 'seabed-fric:slip-y')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-fric:end')
  END SUBROUTINE case_seabed_friction_load

  SUBROUTINE case_seabed_friction_tangent_fd()
    !! Sloped-bathymetry seabed friction derivatives match finite differences.
    TYPE(CD_ModelType) :: model
    TYPE(CD_BathymetryType) :: bathy
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es, j, ii, jj
    REAL(wp) :: q0(6), v0(6), a0(6), q_p(6), q_m(6), v_p(6), v_m(6), loads_p(6), loads_m(6)
    REAL(wp) :: f_ext(6), l0(1), ea(1), rho_a(1), kn(2), cn(2), jq(6, 6), jv(6, 6)
    REAL(wp) :: fdq(6, 6), fdv(6, 6), bx(2), by(2), depth(2, 2), h
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, -10.20_wp, 1.0_wp, 0.5_wp, -10.10_wp]
    v0 = [0.7_wp, -0.3_wp, -0.2_wp, -0.4_wp, 0.6_wp, -0.1_wp]
    a0 = 0.0_wp
    f_ext = 0.0_wp
    l0 = [1.20_wp]
    ea = [30.0_wp]
    rho_a = [1.0_wp]
    kn = [40.0_wp, 60.0_wp]
    cn = [5.0_wp, 7.0_wp]
    bx = [-1.0_wp, 2.0_wp]
    by = [-1.0_wp, 2.0_wp]
    DO jj = 1, 2
      DO ii = 1, 2
        depth(ii, jj) = 10.0_wp - 0.20_wp*bx(ii) + 0.10_wp*by(jj)
      END DO
    END DO
    CALL CD_Init_Bathymetry(bathy, bx, by, depth, es, em)
    CALL require(es == CD_BATHY_OK, 'seabed-fric-fd:bathy')
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       bathymetry=bathy, seabed_kn=kn, seabed_cn=cn, seabed_mu=0.35_wp)
    CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:init: '//TRIM(em))
    ! node 1 sliding (anchor far behind), node 2 sticking (anchor 1 mm away)
    BLOCK
      REAL(wp) :: anchors(2, 2)
      anchors(:, 1) = q0(1:2) - [0.9_wp, -0.6_wp]
      anchors(:, 2) = q0(4:5) - [1.0e-3_wp, 0.5e-3_wp]
      CALL CD_Set_Model_Friction_Anchors(model, anchors, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:anchors')
    END BLOCK
    CALL CD_Calc_Model_CoupledKinematicDerivatives(model, jq, jv, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:analytic: '//TRIM(em))

    h = 1.0e-6_wp
    DO j = 1, 6
      q_p = q0
      q_m = q0
      q_p(j) = q_p(j) + h
      q_m(j) = q_m(j) - h
      CALL CD_Update_Model_CoupledMotion(model, q_p, v0, a0, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:q-plus')
      CALL CD_Calc_Model_CoupledLoads(model, loads_p, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:q-load-plus')
      CALL CD_Update_Model_CoupledMotion(model, q_m, v0, a0, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:q-minus')
      CALL CD_Calc_Model_CoupledLoads(model, loads_m, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:q-load-minus')
      fdq(:, j) = (loads_p - loads_m)/(2.0_wp*h)

      v_p = v0
      v_m = v0
      v_p(j) = v_p(j) + h
      v_m(j) = v_m(j) - h
      CALL CD_Update_Model_CoupledMotion(model, q0, v_p, a0, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:v-plus')
      CALL CD_Calc_Model_CoupledLoads(model, loads_p, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:v-load-plus')
      CALL CD_Update_Model_CoupledMotion(model, q0, v_m, a0, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:v-minus')
      CALL CD_Calc_Model_CoupledLoads(model, loads_m, es, em)
      CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:v-load-minus')
      fdv(:, j) = (loads_p - loads_m)/(2.0_wp*h)
    END DO
    CALL require(nan_max_abs(fdq - jq) < 5.0e-5_wp, 'seabed-fric-fd:dload-dq')
    CALL require(nan_max_abs(fdv - jv) < 5.0e-5_wp, 'seabed-fric-fd:dload-dv')
    CALL CD_Update_Model_CoupledMotion(model, q0, v0, a0, es, em)
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'seabed-fric-fd:end')
    CALL CD_End_Bathymetry(bathy)
  END SUBROUTINE case_seabed_friction_tangent_fd

  SUBROUTINE case_damping_model_load()
    !! Model-owned BA damping contributes to both the initial acceleration and the
    !! coupled load boundary. A rest-length bar with node 2 free and moving +x has
    !! damping force -BA*v/L0 on node 2 and +BA*v/L0 on the fixed support. Because
    !! the element uses consistent mass, the fixed support also reports the inertial
    !! reaction from the free-node acceleration through the off-diagonal mass block.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(3), es
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), ba(1)
    REAL(wp) :: q(6), v(6), a(6), loads(3)
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3]
    q0 = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
    v0 = 0.0_wp
    v0(4) = 2.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [3.0_wp]   ! free-node consistent mass = rho_a*L0/3 = 1
    f_ext = 0.0_wp
    ba = [10.0_wp]

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, ba=ba)
    CALL require(es == CD_MODEL_OK, 'damping:init')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK, 'damping:get')
    CALL require(ABS(a(4) + 20.0_wp) < 1.0e-10_wp, 'damping:free-node-accel')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'damping:loads')
    CALL require(ABS(loads(1) - 30.0_wp) < 1.0e-10_wp, 'damping:fixed-node-load')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'damping:end')
  END SUBROUTINE case_damping_model_load

  SUBROUTINE case_morison_drag_model_load()
    !! Model-owned Morison current drag contributes to the coupled-load boundary.
    !! This mirrors the hydro primitive's one fully wet horizontal element case:
    !! U_y=2 m/s, rho=1000 kg/m^3, D=0.2 m, Cdn=1.2 gives 240 N per node.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), loads(6)
    REAL(wp) :: fluid(3, 2), waterline(2)
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, -1.0_wp, 1.0_wp, 0.0_wp, -1.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [5.0_wp]
    f_ext = 0.0_wp
    fluid = 0.0_wp
    fluid(2, :) = 2.0_wp
    waterline = 0.0_wp

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       fluid_velocity=fluid, drag_waterline_z=waterline, drag_rho=1000.0_wp, &
                       drag_diameter=0.2_wp, drag_cdn=1.2_wp, drag_cdt=0.5_wp)
    CALL require(es == CD_MODEL_OK, 'drag:init')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'drag:loads')
    CALL require(ABS(loads(2) - 240.0_wp) < 1.0e-10_wp, 'drag:node1-y')
    CALL require(ABS(loads(5) - 240.0_wp) < 1.0e-10_wp, 'drag:node2-y')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'drag:end')
  END SUBROUTINE case_morison_drag_model_load

  SUBROUTINE case_morison_drag_element_coefficients()
    !! Section-wise Morison coefficients are owned at element resolution. Two fully
    !! wet unit elements with D=[0.2, 0.4] and U_y=2 m/s give nodal y-loads
    !! [240, 720, 480] N after consistent half-element assembly.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 2), fixed(9), es
    REAL(wp) :: q0(9), v0(9), l0(2), ea(2), rho_a(2), f_ext(9), loads(9)
    REAL(wp) :: fluid(3, 3), waterline(3), diam(2), cdn(2), cdt(2)
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    conn(:, 2) = [2, 3]
    fixed = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    q0 = [0.0_wp, 0.0_wp, -1.0_wp, 1.0_wp, 0.0_wp, -1.0_wp, 2.0_wp, 0.0_wp, -1.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp, 1.0_wp]
    ea = [100.0_wp, 100.0_wp]
    rho_a = [5.0_wp, 5.0_wp]
    f_ext = 0.0_wp
    fluid = 0.0_wp
    fluid(2, :) = 2.0_wp
    waterline = 0.0_wp
    diam = [0.2_wp, 0.4_wp]
    cdn = [1.2_wp, 1.2_wp]
    cdt = [0.5_wp, 0.5_wp]

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       fluid_velocity=fluid, drag_waterline_z=waterline, drag_rho=1000.0_wp, &
                       drag_diameter_elem=diam, drag_cdn_elem=cdn, drag_cdt_elem=cdt)
    CALL require(es == CD_MODEL_OK, 'drag-elem:init: '//TRIM(em))
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'drag-elem:loads')
    CALL require(ABS(loads(2) - 240.0_wp) < 1.0e-10_wp, 'drag-elem:node1-y')
    CALL require(ABS(loads(5) - 720.0_wp) < 1.0e-10_wp, 'drag-elem:node2-y')
    CALL require(ABS(loads(8) - 480.0_wp) < 1.0e-10_wp, 'drag-elem:node3-y')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'drag-elem:end')
  END SUBROUTINE case_morison_drag_element_coefficients

  SUBROUTINE case_froude_krylov_model_load()
    !! Model-owned Froude-Krylov/fluid-inertia loading contributes to coupled loads.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp), PARAMETER :: pi = 3.141592653589793238462643383279502884197_wp
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), loads(6)
    REAL(wp) :: accel(3, 2), waterline(2), base
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, -1.0_wp, 1.0_wp, 0.0_wp, -1.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [5.0_wp]
    f_ext = 0.0_wp
    accel(:, 1) = [1.0_wp, 0.0_wp, 2.0_wp]
    accel(:, 2) = [1.0_wp, 0.0_wp, 2.0_wp]
    waterline = 0.0_wp
    base = 1025.0_wp*0.25_wp*pi*0.252_wp*0.252_wp

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       fluid_acceleration=accel, fk_waterline_z=waterline, fk_rho=1025.0_wp, &
                       fk_diameter=0.252_wp, fk_can=1.0_wp, fk_cat=0.0_wp)
    CALL require(es == CD_MODEL_OK, 'fk-model:init')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'fk-model:loads')
    CALL require(ABS(loads(1) - 0.5_wp*base) < 1.0e-10_wp, 'fk-model:node1-x')
    CALL require(ABS(loads(3) - 2.0_wp*base) < 1.0e-10_wp, 'fk-model:node1-z')
    CALL require(ABS(loads(4) - 0.5_wp*base) < 1.0e-10_wp, 'fk-model:node2-x')
    CALL require(ABS(loads(6) - 2.0_wp*base) < 1.0e-10_wp, 'fk-model:node2-z')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'fk-model:end')
  END SUBROUTINE case_froude_krylov_model_load

  SUBROUTINE case_buoyancy_recovery_model_load()
    !! Model-owned buoyancy recovery removes over-applied submerged weight on dry
    !! nodes. For a fully dry one-element line, each node gets -0.5 rho A g.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp), PARAMETER :: pi = 3.141592653589793238462643383279502884197_wp
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), loads(6), waterline(2), half_buoy
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 1.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [5.0_wp]
    f_ext = 0.0_wp
    waterline = 0.0_wp
    half_buoy = 0.5_wp*1025.0_wp*0.25_wp*pi*0.252_wp*0.252_wp*9.80665_wp

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       buoyancy_waterline_z=waterline, buoyancy_rho=1025.0_wp, &
                       buoyancy_diameter=0.252_wp, buoyancy_gravity=9.80665_wp)
    CALL require(es == CD_MODEL_OK, 'buoyancy-model:init')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'buoyancy-model:loads')
    CALL require(ABS(loads(3) + half_buoy) < 1.0e-10_wp, 'buoyancy-model:node1-z')
    CALL require(ABS(loads(6) + half_buoy) < 1.0e-10_wp, 'buoyancy-model:node2-z')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'buoyancy-model:end')
  END SUBROUTINE case_buoyancy_recovery_model_load

  SUBROUTINE case_hydro_field_update()
    !! The model hydro update API refreshes nodal environmental fields without
    !! rebuilding the model. Updating current speed from 1 to 2 m/s quadruples drag.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(6), es
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), loads(6)
    REAL(wp) :: fluid(3, 2), waterline(2)
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3, 4, 5, 6]
    q0 = [0.0_wp, 0.0_wp, -1.0_wp, 1.0_wp, 0.0_wp, -1.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [5.0_wp]
    f_ext = 0.0_wp
    fluid = 0.0_wp
    fluid(2, :) = 1.0_wp
    waterline = 0.0_wp

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       fluid_velocity=fluid, drag_waterline_z=waterline, drag_rho=1000.0_wp, &
                       drag_diameter=0.2_wp, drag_cdn=1.2_wp, drag_cdt=0.5_wp)
    CALL require(es == CD_MODEL_OK, 'hydro-update:init')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK .AND. ABS(loads(2) - 60.0_wp) < 1.0e-10_wp, 'hydro-update:before')
    fluid(2, :) = 2.0_wp
    CALL CD_Update_Model_Hydro_Fields(model, es, em, fluid_velocity=fluid)
    CALL require(es == CD_MODEL_OK, 'hydro-update:update')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK .AND. ABS(loads(2) - 240.0_wp) < 1.0e-10_wp, 'hydro-update:after')
    fluid(2, :) = 3.0_wp
    waterline(2) = IEEE_VALUE(waterline(2), IEEE_QUIET_NAN)
    CALL CD_Update_Model_Hydro_Fields(model, es, em, fluid_velocity=fluid, drag_waterline_z=waterline)
    CALL require(es == CD_MODEL_BADINPUT, 'hydro-update:rejects-bad-compound-update')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK .AND. ABS(loads(2) - 240.0_wp) < 1.0e-10_wp, &
                 'hydro-update:reject-preserves-previous-fields')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'hydro-update:end')
  END SUBROUTINE case_hydro_field_update

  SUBROUTINE case_hydro_damping_kinematic_derivatives_fd()
    !! Coupled-load derivatives for a pretensioned, hydro-loaded line match
    !! finite differences after each support perturbation recomputes free-node
    !! acceleration. This guards the tangent used by near-taut dynamic coupling.
    INTEGER, PARAMETER :: NE = 2, NN = NE + 1, NDOF = 3*NN, NC = 6
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(NC), es, j
    REAL(wp) :: q0(NDOF), v0(NDOF), f_ext(NDOF), l0(NE), ea(NE), rho_a(NE), ba(NE)
    REAL(wp) :: fluid(3, NN), accel(3, NN), waterline(NN)
    REAL(wp) :: q_c(NC), v_c(NC), a_c(NC), q_p(NC), q_m(NC), v_p(NC), v_m(NC)
    REAL(wp) :: loads_p(NC), loads_m(NC), jq(NC, NC), jv(NC, NC), fdq(NC, NC), fdv(NC, NC)
    REAL(wp) :: h, err_q, err_v, scale_q, scale_v
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    conn(:, 2) = [2, 3]
    fixed = [1, 2, 3, 7, 8, 9]
    q0 = [0.0_wp, 0.0_wp, -40.0_wp, 48.0_wp, 0.6_wp, -42.0_wp, &
          96.0_wp, -0.3_wp, -39.0_wp]
    v0 = [0.20_wp, -0.10_wp, 0.05_wp, -0.15_wp, 0.25_wp, -0.05_wp, &
          0.10_wp, 0.20_wp, -0.02_wp]
    f_ext = 0.0_wp
    l0 = [47.5_wp, 48.0_wp]
    ea = [1.0e5_wp, 1.2e5_wp]
    rho_a = [85.0_wp, 90.0_wp]
    ba = [55.0_wp, 70.0_wp]
    fluid = 0.0_wp
    fluid(1, :) = [0.20_wp, 0.30_wp, 0.25_wp]
    fluid(2, :) = [0.60_wp, 0.55_wp, 0.70_wp]
    accel = 0.0_wp
    accel(1, :) = [0.03_wp, 0.04_wp, 0.02_wp]
    accel(3, :) = [0.10_wp, 0.08_wp, 0.12_wp]
    waterline = 0.0_wp

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, ba=ba, &
                       fluid_velocity=fluid, drag_waterline_z=waterline, drag_rho=1025.0_wp, &
                       drag_diameter=0.09_wp, drag_cdn=1.15_wp, drag_cdt=0.25_wp, &
                       fluid_acceleration=accel, fk_waterline_z=waterline, fk_rho=1025.0_wp, &
                       fk_diameter=0.09_wp, fk_can=1.0_wp, fk_cat=0.0_wp, &
                       added_mass_waterline_z=waterline, added_mass_rho=1025.0_wp, &
                       added_mass_diameter=0.09_wp, added_mass_can=1.0_wp, added_mass_cat=0.0_wp, &
                       buoyancy_waterline_z=waterline, buoyancy_rho=1025.0_wp, &
                       buoyancy_diameter=0.09_wp, buoyancy_gravity=9.80665_wp)
    CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:init: '//TRIM(em))
    CALL CD_Recompute_Model_Acceleration(model, es, em)
    CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:recompute')
    CALL CD_Get_Model_CoupledMotion(model, q_c, v_c, a_c, es, em)
    CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:get-motion')
    CALL CD_Calc_Model_CoupledKinematicDerivatives(model, jq, jv, es, em)
    CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:analytic: '//TRIM(em))
    CALL require(ALLOCATED(model%load_force_work) .AND. ALLOCATED(model%load_jq_work) .AND. &
                 ALLOCATED(model%load_jv_work), 'hydro-deriv-fd:reuses-load-workspace')

    h = 1.0e-6_wp
    DO j = 1, NC
      q_p = q_c
      q_m = q_c
      q_p(j) = q_p(j) + h
      q_m(j) = q_m(j) - h
      CALL CD_Update_Model_CoupledMotion(model, q_p, v_c, a_c, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:q-plus-update')
      CALL CD_Recompute_Model_Acceleration(model, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:q-plus-recompute')
      CALL CD_Calc_Model_CoupledLoads(model, loads_p, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:q-plus-loads')
      CALL CD_Update_Model_CoupledMotion(model, q_m, v_c, a_c, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:q-minus-update')
      CALL CD_Recompute_Model_Acceleration(model, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:q-minus-recompute')
      CALL CD_Calc_Model_CoupledLoads(model, loads_m, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:q-minus-loads')
      fdq(:, j) = (loads_p - loads_m)/(2.0_wp*h)

      v_p = v_c
      v_m = v_c
      v_p(j) = v_p(j) + h
      v_m(j) = v_m(j) - h
      CALL CD_Update_Model_CoupledMotion(model, q_c, v_p, a_c, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:v-plus-update')
      CALL CD_Recompute_Model_Acceleration(model, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:v-plus-recompute')
      CALL CD_Calc_Model_CoupledLoads(model, loads_p, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:v-plus-loads')
      CALL CD_Update_Model_CoupledMotion(model, q_c, v_m, a_c, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:v-minus-update')
      CALL CD_Recompute_Model_Acceleration(model, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:v-minus-recompute')
      CALL CD_Calc_Model_CoupledLoads(model, loads_m, es, em)
      CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:v-minus-loads')
      fdv(:, j) = (loads_p - loads_m)/(2.0_wp*h)
    END DO

    err_q = nan_max_abs(fdq - jq)
    err_v = nan_max_abs(fdv - jv)
    scale_q = MAX(1.0_wp, MAX(nan_max_abs(fdq), nan_max_abs(jq)))
    scale_v = MAX(1.0_wp, MAX(nan_max_abs(fdv), nan_max_abs(jv)))
    CALL require(err_q/scale_q < 2.0e-4_wp, 'hydro-deriv-fd:dload-dq')
    CALL require(err_v/scale_v < 2.0e-4_wp, 'hydro-deriv-fd:dload-dv')

    CALL CD_Update_Model_CoupledMotion(model, q_c, v_c, a_c, es, em)
    CALL CD_Recompute_Model_Acceleration(model, es, em)
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'hydro-deriv-fd:end')
    CALL require(.NOT. ALLOCATED(model%load_force_work) .AND. .NOT. ALLOCATED(model%load_jq_work) .AND. &
                 .NOT. ALLOCATED(model%load_jv_work), 'hydro-deriv-fd:releases-load-workspace')
  END SUBROUTINE case_hydro_damping_kinematic_derivatives_fd

  SUBROUTINE case_prescribed_hydro_added_mass_step()
    !! Regression for the dense dynamic branch used by moving-fairlead,
    !! still-water hydrodynamic models: prescribed endpoint motion, Morison drag,
    !! and configuration-dependent added mass must converge, keep finite state,
    !! honor the prescribed support exactly, and differ from the dry-mass response.
    INTEGER, PARAMETER :: NE = 4, NN = NE + 1, NDOF = 3*NN, NC = 6
    TYPE(CD_ModelType) :: wet_model, dry_model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(NC), es, n_iter
    REAL(wp) :: q0(NDOF), v0(NDOF), f_ext(NDOF), l0(NE), ea(NE), rho_a(NE)
    REAL(wp) :: fluid(3, NN), waterline(NN), q_target(NDOF), v_target(NDOF), a_target(NDOF)
    REAL(wp) :: q_wet(NDOF), v_wet(NDOF), a_wet(NDOF), q_dry(NDOF), v_dry(NDOF), a_dry(NDOF)
    REAL(wp) :: dt
    LOGICAL :: converged, stalled
    CHARACTER(200) :: em

    CALL straight_cable(NE, conn, l0, ea, rho_a, q0, fixed)
    f_ext = 0.0_wp
    v0 = 0.0_wp
    fluid = 0.0_wp
    waterline = 5.0_wp
    dt = 0.02_wp
    cfg%max_iter = 80
    cfg%armijo_max_backtracks = 14

    CALL CD_Init_Model(wet_model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       fluid_velocity=fluid, drag_waterline_z=waterline, drag_rho=1025.0_wp, &
                       drag_diameter=0.16_wp, drag_cdn=1.2_wp, drag_cdt=0.25_wp, &
                       added_mass_waterline_z=waterline, added_mass_rho=1025.0_wp, &
                       added_mass_diameter=0.16_wp, added_mass_can=1.0_wp, added_mass_cat=0.0_wp)
    CALL require(es == CD_MODEL_OK, 'prescribed-hydro:init-wet: '//TRIM(em))
    CALL CD_Init_Model(dry_model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       fluid_velocity=fluid, drag_waterline_z=waterline, drag_rho=1025.0_wp, &
                       drag_diameter=0.16_wp, drag_cdn=1.2_wp, drag_cdt=0.25_wp)
    CALL require(es == CD_MODEL_OK, 'prescribed-hydro:init-dry: '//TRIM(em))

    q_target = q0
    v_target = 0.0_wp
    a_target = 0.0_wp
    q_target(NDOF - 1) = q0(NDOF - 1) + 0.015_wp
    q_target(NDOF) = q0(NDOF) + 0.025_wp
    v_target(NDOF - 1) = 0.75_wp
    v_target(NDOF) = 1.25_wp

    CALL CD_Step_Model(wet_model, dt, converged, stalled, n_iter, es, em, &
                       prescribed_q=q_target, prescribed_v=v_target, prescribed_a=a_target)
    CALL require(es == CD_MODEL_OK .AND. converged .AND. .NOT. stalled, 'prescribed-hydro:wet-step: '//TRIM(em))
    CALL require(ALLOCATED(wet_model%dynamic_workspace%M_add_band) .AND. &
                 ALLOCATED(wet_model%dynamic_workspace%dMa_a_dq_band), 'prescribed-hydro:banded-added-mass-path')
    CALL CD_Get_Model_State(wet_model, q_wet, v_wet, a_wet, es, em)
    CALL require(es == CD_MODEL_OK, 'prescribed-hydro:get-wet')
    CALL require(ALL(IEEE_IS_FINITE(q_wet)) .AND. ALL(IEEE_IS_FINITE(v_wet)) .AND. ALL(IEEE_IS_FINITE(a_wet)), &
                 'prescribed-hydro:wet-finite')
    CALL require(nan_max_abs(q_wet(fixed) - q_target(fixed)) < 1.0e-12_wp, 'prescribed-hydro:q-support')
    CALL require(nan_max_abs(v_wet(fixed) - v_target(fixed)) < 1.0e-12_wp, 'prescribed-hydro:v-support')
    CALL require(nan_max_abs(a_wet(fixed) - a_target(fixed)) < 1.0e-12_wp, 'prescribed-hydro:a-support')

    CALL CD_Step_Model(dry_model, dt, converged, stalled, n_iter, es, em, &
                       prescribed_q=q_target, prescribed_v=v_target, prescribed_a=a_target)
    CALL require(es == CD_MODEL_OK .AND. converged .AND. .NOT. stalled, 'prescribed-hydro:dry-step: '//TRIM(em))
    CALL CD_Get_Model_State(dry_model, q_dry, v_dry, a_dry, es, em)
    CALL require(es == CD_MODEL_OK, 'prescribed-hydro:get-dry')
    CALL require(nan_max_abs(q_wet(4:NDOF - 3) - q_dry(4:NDOF - 3)) > 1.0e-8_wp, &
                 'prescribed-hydro:added-mass-affects-interior')

    CALL CD_End_Model(wet_model, es, em)
    CALL require(es == CD_MODEL_OK, 'prescribed-hydro:end-wet')
    CALL CD_End_Model(dry_model, es, em)
    CALL require(es == CD_MODEL_OK, 'prescribed-hydro:end-dry')
  END SUBROUTINE case_prescribed_hydro_added_mass_step

  SUBROUTINE case_banded_model_load_matches_dense_structured_bathymetry()
    !! Dense and banded model-owned load paths must agree for the local physics
    !! contributors the banded path duplicates: structured-bathymetry seabed
    !! contact/damping/friction, BA damping, Morison drag, Froude-Krylov, and
    !! buoyancy recovery. The dense side carries zero added-mass coefficients only
    !! to force the dense callback route while preserving the same physics.
    TYPE(CD_ModelType) :: band_model, dense_model
    TYPE(CD_BathymetryType) :: bathy
    TYPE(GenAlphaConfig) :: cfg
    INTEGER, PARAMETER :: NE = 2, NN = 3, NDOF = 9
    INTEGER :: conn(2, NE), fixed(6), es, n_iter_band, n_iter_dense
    REAL(wp) :: l0(NE), ea(NE), rho_a(NE), q0(NDOF), v0(NDOF), f_ext(NDOF)
    REAL(wp) :: kn(NN), cn(NN), ba(NE), fluid(3, NN), accel(3, NN), waterline(NN)
    REAL(wp) :: xgrid(3), ygrid(3), depth(3, 3)
    REAL(wp) :: qb(NDOF), vb(NDOF), ab(NDOF), qd(NDOF), vd(NDOF), ad(NDOF)
    LOGICAL :: conv_band, stalled_band, conv_dense, stalled_dense
    CHARACTER(240) :: em

    conn = RESHAPE([1, 2, 2, 3], [2, NE])
    fixed = [1, 2, 3, 7, 8, 9]
    l0 = [0.95_wp, 0.95_wp]
    ea = [1.0e5_wp, 1.0e5_wp]
    rho_a = [12.0_wp, 12.0_wp]
    q0 = [0.0_wp, 0.0_wp, -100.05_wp, 1.0_wp, 0.0_wp, -100.10_wp, &
          2.0_wp, 0.0_wp, -100.05_wp]
    v0 = 0.0_wp
    v0(4:6) = [0.12_wp, -0.08_wp, -0.04_wp]
    f_ext = 0.0_wp
    kn = [5000.0_wp, 5000.0_wp, 5000.0_wp]
    cn = [120.0_wp, 120.0_wp, 120.0_wp]
    ba = [15.0_wp, 15.0_wp]
    fluid = 0.0_wp
    fluid(2, :) = [0.8_wp, 0.9_wp, 0.8_wp]
    accel = 0.0_wp
    accel(3, :) = [0.12_wp, 0.14_wp, 0.12_wp]
    waterline = 0.0_wp
    cfg%max_iter = 40

    xgrid = [-10.0_wp, 0.0_wp, 10.0_wp]
    ygrid = [-10.0_wp, 0.0_wp, 10.0_wp]
    depth = RESHAPE([100.20_wp, 100.00_wp, 100.35_wp, &
                     100.05_wp, 100.12_wp, 100.48_wp, &
                     100.32_wp, 100.18_wp, 100.71_wp], SHAPE(depth))
    CALL CD_Init_Bathymetry(bathy, xgrid, ygrid, depth, es, em)
    CALL require(es == CD_BATHY_OK, 'band-dense:bathy-init')

    CALL CD_Init_Model(band_model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       bathymetry=bathy, seabed_kn=kn, seabed_cn=cn, seabed_mu=0.35_wp, ba=ba, &
                       fluid_velocity=fluid, drag_waterline_z=waterline, drag_rho=1025.0_wp, &
                       drag_diameter=0.18_wp, drag_cdn=1.1_wp, drag_cdt=0.2_wp, &
                       fluid_acceleration=accel, fk_waterline_z=waterline, fk_rho=1025.0_wp, &
                       fk_diameter=0.18_wp, fk_can=1.0_wp, fk_cat=0.1_wp, &
                       buoyancy_waterline_z=waterline, buoyancy_rho=1025.0_wp, &
                       buoyancy_diameter=0.18_wp, buoyancy_gravity=9.80665_wp)
    CALL require(es == CD_MODEL_OK, 'band-dense:band-init: '//TRIM(em))
    CALL CD_Init_Model(dense_model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       bathymetry=bathy, seabed_kn=kn, seabed_cn=cn, seabed_mu=0.35_wp, ba=ba, &
                       fluid_velocity=fluid, drag_waterline_z=waterline, drag_rho=1025.0_wp, &
                       drag_diameter=0.18_wp, drag_cdn=1.1_wp, drag_cdt=0.2_wp, &
                       fluid_acceleration=accel, fk_waterline_z=waterline, fk_rho=1025.0_wp, &
                       fk_diameter=0.18_wp, fk_can=1.0_wp, fk_cat=0.1_wp, &
                       buoyancy_waterline_z=waterline, buoyancy_rho=1025.0_wp, &
                       buoyancy_diameter=0.18_wp, buoyancy_gravity=9.80665_wp, &
                       added_mass_waterline_z=waterline, added_mass_rho=1025.0_wp, &
                       added_mass_diameter=0.18_wp, added_mass_can=0.0_wp, added_mass_cat=0.0_wp)
    CALL require(es == CD_MODEL_OK, 'band-dense:dense-init: '//TRIM(em))

    CALL CD_Step_Model(band_model, 0.001_wp, conv_band, stalled_band, n_iter_band, es, em)
    CALL require(es == CD_MODEL_OK .AND. conv_band .AND. .NOT. stalled_band, 'band-dense:band-step: '//TRIM(em))
    CALL CD_Step_Model(dense_model, 0.001_wp, conv_dense, stalled_dense, n_iter_dense, es, em)
    CALL require(es == CD_MODEL_OK .AND. conv_dense .AND. .NOT. stalled_dense, 'band-dense:dense-step: '//TRIM(em))
    CALL CD_Get_Model_State(band_model, qb, vb, ab, es, em)
    CALL require(es == CD_MODEL_OK, 'band-dense:band-get')
    CALL CD_Get_Model_State(dense_model, qd, vd, ad, es, em)
    CALL require(es == CD_MODEL_OK, 'band-dense:dense-get')
    CALL require(nan_max_abs(qb - qd) < 1.0e-10_wp, 'band-dense:q-match')
    CALL require(nan_max_abs(vb - vd) < 1.0e-10_wp, 'band-dense:v-match')
    CALL require(nan_max_abs(ab - ad) < 1.0e-8_wp, 'band-dense:a-match')

    CALL CD_End_Model(band_model, es, em)
    CALL require(es == CD_MODEL_OK, 'band-dense:band-end')
    CALL CD_End_Model(dense_model, es, em)
    CALL require(es == CD_MODEL_OK, 'band-dense:dense-end')
    CALL CD_End_Bathymetry(bathy)
  END SUBROUTINE case_banded_model_load_matches_dense_structured_bathymetry

  SUBROUTINE case_added_mass_model_accel()
    !! Model-owned added mass augments the free-DOF inertia in initial
    !! acceleration. With node 1 fixed, structural free z-mass is rho_a*L/3=1;
    !! the fully wet normal added mass contributes rho*pi*D^2/4/3 at node 2.
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, 1), fixed(3), es
    REAL(wp), PARAMETER :: pi = 3.141592653589793238462643383279502884197_wp
    REAL(wp) :: q0(6), v0(6), l0(1), ea(1), rho_a(1), f_ext(6), waterline(2)
    REAL(wp) :: q(6), v(6), a(6), expected
    CHARACTER(200) :: em

    conn(:, 1) = [1, 2]
    fixed = [1, 2, 3]
    q0 = [0.0_wp, 0.0_wp, -1.0_wp, 1.0_wp, 0.0_wp, -1.0_wp]
    v0 = 0.0_wp
    l0 = [1.0_wp]
    ea = [100.0_wp]
    rho_a = [3.0_wp]
    f_ext = 0.0_wp
    f_ext(6) = 1.0_wp
    waterline = 0.0_wp
    expected = 1.0_wp/(1.0_wp + 1025.0_wp*pi*0.252_wp*0.252_wp/4.0_wp/3.0_wp)

    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       added_mass_waterline_z=waterline, added_mass_rho=1025.0_wp, &
                       added_mass_diameter=0.252_wp, added_mass_can=1.0_wp, added_mass_cat=0.0_wp)
    CALL require(es == CD_MODEL_OK, 'added-mass:init')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK, 'added-mass:get')
    CALL require(ABS(a(6) - expected) < 1.0e-12_wp, 'added-mass:free-z-accel')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'added-mass:end')
  END SUBROUTINE case_added_mass_model_accel

  SUBROUTINE case_accel_derivative_band_parity()
    !! The banded coupled acceleration derivative equals the dense Schur complement
    !! -(M_cc - M_cf M_ff^-1 M_fc) with M the structural plus added mass, on a line
    !! that crosses the waterline, has uneven segments and an interior coupled node.
    !! Repeated queries return the cached result; a new position (added mass) or a
    !! new segment length invalidates it.
    INTEGER, PARAMETER :: NE = 12, NN = NE + 1, NDOF = 3*NN, NC = 9
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(NC), es, e, i
    REAL(wp) :: l0(NE), ea(NE), rho_a(NE), q0(NDOF), v0(NDOF), f_ext(NDOF), waterline(NN)
    REAL(wp) :: da(NC, NC), da2(NC, NC), ref(NC, NC), qc(NC), zc(NC)
    CHARACTER(200) :: em

    DO e = 1, NE
      conn(:, e) = [e, e + 1]
      l0(e) = 0.9_wp + 0.03_wp*REAL(MOD(e, 4), wp)
      ea(e) = 1.0e6_wp
      rho_a(e) = 8.0_wp + REAL(MOD(e, 3), wp)
    END DO
    DO i = 1, NN
      q0(3*i - 2:3*i) = [REAL(i - 1, wp), 0.1_wp*SIN(REAL(i, wp)), -3.0_wp + 0.5_wp*REAL(i - 1, wp)]
    END DO
    ! both ends and interior node 7 are coupled
    fixed = [1, 2, 3, 19, 20, 21, NDOF - 2, NDOF - 1, NDOF]
    v0 = 0.0_wp
    f_ext = 0.0_wp
    waterline = 0.0_wp
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       added_mass_waterline_z=waterline, added_mass_rho=1025.0_wp, &
                       added_mass_diameter=0.15_wp, added_mass_can=1.0_wp, added_mass_cat=0.5_wp)
    CALL require(es == CD_MODEL_OK, 'accel-band:init')

    CALL CD_Calc_Model_CoupledAccelDerivative(model, da, es, em)
    CALL require(es == CD_MODEL_OK, 'accel-band:query')
    CALL accel_dense_reference(conn, l0, rho_a, waterline, fixed, q0, ref)
    CALL require(nan_max_abs(da - ref) <= 1.0e-12_wp*nan_max_abs(ref), 'accel-band:dense-parity')
    CALL CD_Calc_Model_CoupledAccelDerivative(model, da2, es, em)
    CALL require(es == CD_MODEL_OK .AND. ALL(da2 <= da .AND. da2 >= da), 'accel-band:cached-repeat')

    ! a coupled-motion update moves the positions, so the added mass changes
    qc = q0(fixed)
    zc = 0.0_wp
    qc(3) = qc(3) - 0.4_wp
    qc(6) = qc(6) + 0.7_wp
    CALL CD_Update_Model_CoupledMotion(model, qc, zc, zc, es, em)
    CALL require(es == CD_MODEL_OK, 'accel-band:update-motion')
    q0(fixed) = qc
    CALL CD_Calc_Model_CoupledAccelDerivative(model, da, es, em)
    CALL accel_dense_reference(conn, l0, rho_a, waterline, fixed, q0, ref)
    CALL require(es == CD_MODEL_OK .AND. nan_max_abs(da - ref) <= 1.0e-12_wp*nan_max_abs(ref), &
                 'accel-band:position-invalidates')

    ! a segment-length change (line control) changes the structural mass
    CALL CD_Update_Model_SegmentLength(model, 3, 1.2_wp, 0.0_wp, es, em)
    CALL require(es == CD_MODEL_OK, 'accel-band:update-l0')
    l0(3) = 1.2_wp
    CALL CD_Calc_Model_CoupledAccelDerivative(model, da, es, em)
    CALL accel_dense_reference(conn, l0, rho_a, waterline, fixed, q0, ref)
    CALL require(es == CD_MODEL_OK .AND. nan_max_abs(da - ref) <= 1.0e-12_wp*nan_max_abs(ref), &
                 'accel-band:length-invalidates')
    CALL CD_End_Model(model, es, em)
  END SUBROUTINE case_accel_derivative_band_parity

  SUBROUTINE accel_dense_reference(conn, l0, rho_a, waterline, fixed, q, out)
    !! Dense -(M_cc - M_cf M_ff^-1 M_fc) with M the structural plus added mass
    !! (rho 1025, D 0.15, Can 1, Cat 0.5), summed over the free DOFs in order.
    INTEGER, INTENT(IN) :: conn(:, :), fixed(:)
    REAL(wp), INTENT(IN) :: l0(:), rho_a(:), waterline(:), q(:)
    REAL(wp), INTENT(OUT) :: out(:, :)
    REAL(wp), ALLOCATABLE :: M(:, :), M_add(:, :), mff(:, :), x(:, :)
    INTEGER, ALLOCATABLE :: free(:)
    INTEGER :: n, nc, nf, g, r, c, k, es
    CHARACTER(200) :: em

    n = SIZE(q)
    nc = SIZE(fixed)
    ALLOCATE (M(n, n), M_add(n, n), mff(n - nc, n - nc), x(n - nc, nc), free(n - nc))
    CALL CD_Assemble_Cable_Mass(conn, l0, rho_a, M, es, em)
    CALL CD_Cable_Added_Mass_Matrix(q, conn, l0, waterline, 1025.0_wp, 0.15_wp, 1.0_wp, 0.5_wp, M_add, es, em)
    M = M + M_add
    nf = 0
    DO g = 1, n
      IF (ANY(fixed == g)) CYCLE
      nf = nf + 1
      free(nf) = g
    END DO
    mff = M(free, free)
    x = M(free, fixed)
    CALL CD_Solve_Dense_As_Banded_Multiple(mff, x, es, em)
    DO c = 1, nc
      DO r = 1, nc
        out(r, c) = -M(fixed(r), fixed(c))
        DO k = 1, nf
          out(r, c) = out(r, c) + M(fixed(r), free(k))*x(k, c)
        END DO
      END DO
    END DO
  END SUBROUTINE accel_dense_reference

  SUBROUTINE case_line_model_init_from_sections()
    !! Product-level initialisation from line types + sections: a suspended taut
    !! line-object builds a static IC, stores a dynamic model, reports the endpoint
    !! coupled DOFs, and remains at static equilibrium for one no-prescribed step.
    TYPE(CD_ModelType) :: model
    TYPE(CD_LineType) :: lts(1)
    TYPE(CD_LineSection) :: secs(1)
    TYPE(CableSolverConfig) :: static_cfg
    TYPE(GenAlphaConfig) :: dynamic_cfg
    REAL(wp), PARAMETER :: factors(4) = [0.25_wp, 0.5_wp, 0.75_wp, 1.0_wp]
    REAL(wp), PARAMETER :: anchor(3) = [0.0_wp, 0.0_wp, 5.0_wp]
    REAL(wp), PARAMETER :: fairlead(3) = [60.0_wp, 0.0_wp, 15.0_wp]
    REAL(wp), ALLOCATABLE :: q(:), v(:), a(:)
    INTEGER :: es, nc, nd, nelem, ndof, n_iter
    LOGICAL :: conv, stalled
    CHARACTER(240) :: em

    lts(1) = CD_LineType(ea=5.0e7_wp, mass_per_length=30.0_wp, diameter=0.08_wp)
    secs(1) = CD_LineSection(line_type=1, length=58.0_wp, n_segments=20)
    CALL CD_Init_Line_Model(model, anchor, fairlead, lts, secs, 9.81_wp, 1025.0_wp, .FALSE., &
                            static_cfg, dynamic_cfg, factors, es, em)
    CALL require(es == CD_MODEL_OK .AND. CD_Model_Is_Initialized(model), 'line-init:ok')
    nc = CD_Model_NCoupledDOF(model, es, em)
    CALL require(es == CD_MODEL_OK .AND. nc == 6, 'line-init:ncoupled-endpoints')
    nd = CD_Model_NDOF(model, es, em)
    nelem = CD_Model_NElem(model, es, em)
    CALL require(es == CD_MODEL_OK .AND. nd == 63 .AND. nelem == 20, 'line-init:dimensions')

    ndof = 3*(20 + 1)
    ALLOCATE (q(ndof), v(ndof), a(ndof))
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK, 'line-init:get')
    CALL require(nan_max_abs(q(1:3) - anchor) < 1.0e-9_wp, 'line-init:anchor')
    CALL require(nan_max_abs(q(ndof - 2:ndof) - fairlead) < 1.0e-9_wp, 'line-init:fairlead')
    CALL require(nan_max_abs(v) < 1.0e-12_wp, 'line-init:vzero')
    CALL require(nan_max_abs(a) < 1.0e-6_wp, 'line-init:consistent-accel')

    CALL CD_Step_Model(model, 0.01_wp, conv, stalled, n_iter, es, em)
    CALL require(es == CD_MODEL_OK .AND. conv .AND. .NOT. stalled, 'line-init:step')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK, 'line-init:get-after-step')
    CALL require(nan_max_abs(q(1:3) - anchor) < 1.0e-9_wp, 'line-init:anchor-after-step')
    CALL require(nan_max_abs(q(ndof - 2:ndof) - fairlead) < 1.0e-9_wp, 'line-init:fairlead-after-step')

    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'line-init:end')
  END SUBROUTINE case_line_model_init_from_sections

  SUBROUTINE case_grounded_line_model_with_seabed()
    !! A grounded composite line can now initialise and step through CD_Model with
    !! the same seabed contributor used by the static IC. This is the first dynamic
    !! model gate for a grounded OrcaFlex-style line object.
    TYPE(CD_ModelType) :: model
    TYPE(CD_LineType) :: lts(2)
    TYPE(CD_LineSection) :: secs(2)
    TYPE(CableSolverConfig) :: static_cfg
    TYPE(GenAlphaConfig) :: dynamic_cfg
    INTEGER, ALLOCATABLE :: conn(:, :)
    REAL(wp), ALLOCATABLE :: l0(:), ea(:), mpl(:), dia(:), kn(:), q(:), v(:), a(:)
    REAL(wp), PARAMETER :: factors(4) = [0.25_wp, 0.5_wp, 0.75_wp, 1.0_wp]
    REAL(wp), PARAMETER :: anchor(3) = [0.0_wp, 0.0_wp, 0.0_wp]
    REAL(wp), PARAMETER :: fairlead(3) = [60.0_wp, 0.0_wp, 35.0_wp]
    INTEGER :: es, ndof, n_iter
    LOGICAL :: conv, stalled
    CHARACTER(240) :: em

    lts(1) = CD_LineType(ea=5.0e6_wp, mass_per_length=50.0_wp, diameter=0.10_wp)
    lts(2) = CD_LineType(ea=8.0e6_wp, mass_per_length=20.0_wp, diameter=0.08_wp)
    secs(1) = CD_LineSection(line_type=1, length=40.0_wp, n_segments=8)
    secs(2) = CD_LineSection(line_type=2, length=40.0_wp, n_segments=12)

    CALL CD_Build_Line_Mesh(secs, lts, .FALSE., conn, l0, ea, mpl, dia, es, em)
    CALL require(es == CD_LINE_OK, 'grounded:mesh')
    CALL CD_Nodal_Seabed_Stiffness(1.0e5_wp, dia, l0, kn, es, em)
    CALL require(es == CD_LINE_OK, 'grounded:kn')

    CALL CD_Init_Line_Model(model, anchor, fairlead, lts, secs, 9.81_wp, 1025.0_wp, .FALSE., &
                            static_cfg, dynamic_cfg, factors, es, em, &
                            seabed_z_floor=0.0_wp, seabed_kn=kn)
    CALL require(es == CD_MODEL_OK .AND. CD_Model_Is_Initialized(model), 'grounded:init')
    ndof = CD_Model_NDOF(model, es, em)
    CALL require(es == CD_MODEL_OK .AND. ndof == 63, 'grounded:ndof')
    ALLOCATE (q(ndof), v(ndof), a(ndof))
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK, 'grounded:get')
    CALL require(MINVAL(q(6:ndof - 3:3)) > -0.1_wp, 'grounded:no-through-seabed')
    CALL require(MINVAL(q(6:ndof - 3:3)) < 0.01_wp, 'grounded:touches-down')

    CALL CD_Step_Model(model, 0.005_wp, conv, stalled, n_iter, es, em)
    CALL require(es == CD_MODEL_OK .AND. conv .AND. .NOT. stalled, 'grounded:step')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK, 'grounded:get-after-step')
    CALL require(nan_max_abs(q(1:3) - anchor) < 1.0e-9_wp, 'grounded:anchor-after-step')
    CALL require(nan_max_abs(q(ndof - 2:ndof) - fairlead) < 1.0e-9_wp, 'grounded:fairlead-after-step')

    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK, 'grounded:end')
  END SUBROUTINE case_grounded_line_model_with_seabed

  SUBROUTINE case_end_force_is_coupled_load()
    !! FairTen/AnchTen source: CD_Get_Model_EndForces is the actual force the line exerts on
    !! its held ends, axial damping and drag at the actual velocity included. With every
    !! acceleration zero (a line at its static state whose fairlead is given a velocity but
    !! no acceleration) it equals the coupled load exactly, and it differs from the at-rest
    !! end force by the damping and drag of the moving end element.
    TYPE(CD_ModelType) :: model
    TYPE(CD_LineType) :: lts(1)
    TYPE(CD_LineSection) :: secs(1)
    TYPE(CableSolverConfig) :: static_cfg
    TYPE(GenAlphaConfig) :: dynamic_cfg
    REAL(wp), PARAMETER :: factors(4) = [0.25_wp, 0.5_wp, 0.75_wp, 1.0_wp]
    REAL(wp) :: ba(30), fluid(3, 31), wl(31), qc(6), vc(6), ac(6), loads(6), f_a(3), f_b(3), r_a(3), r_b(3)
    INTEGER :: es
    CHARACTER(240) :: em

    lts(1) = CD_LineType(ea=5.0e8_wp, mass_per_length=80.0_wp, diameter=0.1_wp)
    secs(1) = CD_LineSection(line_type=1, length=300.0_wp, n_segments=30)
    ba = 2.0e6_wp
    fluid = 0.0_wp
    wl = 0.0_wp
    CALL CD_Init_Line_Model(model, [250.0_wp, 0.0_wp, -150.0_wp], [0.0_wp, 0.0_wp, -10.0_wp], lts, secs, &
                            9.80665_wp, 1025.0_wp, .FALSE., static_cfg, dynamic_cfg, factors, es, em, ba=ba, &
                            fluid_velocity=fluid, drag_waterline_z=wl, drag_rho=1025.0_wp, drag_diameter=0.1_wp, &
                            drag_cdn=1.2_wp, drag_cdt=0.5_wp, dynamic_tension_only=.TRUE.)
    CALL require(es == CD_MODEL_OK, 'endforce:init: '//TRIM(em))
    IF (es /= CD_MODEL_OK) RETURN
    CALL CD_Get_Model_EndForces(model, r_a, r_b, es, em)
    CALL require(es == CD_MODEL_OK, 'endforce:rest')
    CALL CD_Get_Model_CoupledMotion(model, qc, vc, ac, es, em)
    CALL require(es == CD_MODEL_OK, 'endforce:motion')
    vc = 0.0_wp
    vc(4:6) = [0.8_wp, 0.3_wp, -0.5_wp]
    ac = 0.0_wp
    CALL CD_Update_Model_CoupledMotion(model, qc, vc, ac, es, em)
    CALL require(es == CD_MODEL_OK, 'endforce:update')
    CALL CD_Get_Model_EndForces(model, f_a, f_b, es, em)
    CALL require(es == CD_MODEL_OK, 'endforce:moving')
    CALL CD_Calc_Model_CoupledLoads(model, loads, es, em)
    CALL require(es == CD_MODEL_OK, 'endforce:coupled')
    CALL require(nan_max_abs(loads(1:3) - f_a) <= 1.0e-8_wp*NORM2(f_a) .AND. &
                 nan_max_abs(loads(4:6) - f_b) <= 1.0e-8_wp*NORM2(f_b), 'endforce:equals-coupled-load-at-zero-accel')
    CALL require(NORM2(f_b - r_b) > 1.0e-3_wp*NORM2(r_b), 'endforce:includes-damping-and-drag')
    CALL require(nan_max_abs(f_a - r_a) <= 1.0e-8_wp*NORM2(r_a), 'endforce:still-anchor-unchanged')
    CALL CD_End_Model(model, es, em)
  END SUBROUTINE case_end_force_is_coupled_load

  SUBROUTINE case_line_model_fail_closed()
    !! The line-object model API rejects unsupported/invalid physics at the boundary:
    !! finite-EI sections and impossible static initialisation do not leave a half
    !! initialised model behind.
    TYPE(CD_ModelType) :: model
    TYPE(CD_LineType) :: lts(1)
    TYPE(CD_LineSection) :: secs(1)
    TYPE(CableSolverConfig) :: static_cfg
    TYPE(GenAlphaConfig) :: dynamic_cfg
    REAL(wp), PARAMETER :: factors(4) = [0.25_wp, 0.5_wp, 0.75_wp, 1.0_wp]
    REAL(wp), PARAMETER :: bad_factors(2) = [0.25_wp, 0.5_wp]
    REAL(wp), PARAMETER :: bad_kn(2) = [1000.0_wp, 1000.0_wp]
    INTEGER :: es
    CHARACTER(240) :: em

    lts(1) = CD_LineType(ea=5.0e7_wp, mass_per_length=30.0_wp, diameter=0.08_wp, ei=1.0e3_wp)
    secs(1) = CD_LineSection(line_type=1, length=58.0_wp, n_segments=20)
    CALL CD_Init_Line_Model(model, [0.0_wp, 0.0_wp, 5.0_wp], [60.0_wp, 0.0_wp, 15.0_wp], &
                            lts, secs, 9.81_wp, 1025.0_wp, .FALSE., static_cfg, dynamic_cfg, &
                            factors, es, em)
    CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'line-fail:finite-ei')

    lts(1) = CD_LineType(ea=5.0e7_wp, mass_per_length=30.0_wp, diameter=0.08_wp)
    secs(1) = CD_LineSection(line_type=1, length=58.0_wp, n_segments=20)
    CALL CD_Init_Line_Model(model, [0.0_wp, 0.0_wp, 5.0_wp], [60.0_wp, 0.0_wp, 15.0_wp], &
                            lts, secs, 9.81_wp, 1025.0_wp, .FALSE., static_cfg, dynamic_cfg, &
                            bad_factors, es, em)
    CALL require(es == CD_MODEL_SOLVEFAIL .AND. .NOT. CD_Model_Is_Initialized(model), 'line-fail:bad-factors')

    CALL CD_Init_Line_Model(model, [0.0_wp, 0.0_wp, 5.0_wp], [60.0_wp, 0.0_wp, 15.0_wp], &
                            lts, secs, 9.81_wp, 1025.0_wp, .FALSE., static_cfg, dynamic_cfg, &
                            factors, es, em, seabed_z_floor=0.0_wp, seabed_kn=bad_kn)
    CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'line-fail:bad-seabed-kn')
  END SUBROUTINE case_line_model_fail_closed

  SUBROUTINE case_fail_closed()
    !! Uninitialised queries and malformed input are rejected with explicit model
    !! ErrStat values, not silently accepted.
    INTEGER, PARAMETER :: NE = 1, NN = NE + 1, NDOF = 3*NN
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    INTEGER :: conn(2, NE), fixed(6), coupled(6), es, n_iter, nc, nd, nelem
    REAL(wp) :: l0(NE), ea(NE), rho_a(NE), q0(NDOF), v0(NDOF), f_ext(NDOF)
    REAL(wp) :: q(NDOF), v(NDOF), a(NDOF), q_save(NDOF), tension(NE)
    LOGICAL :: conv, stalled
    CHARACTER(200) :: em

    nc = CD_Model_NCoupledDOF(model, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED .AND. nc == 0, 'fail:ncoupled-uninit')
    nd = CD_Model_NDOF(model, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED .AND. nd == 0, 'fail:ndof-uninit')
    nelem = CD_Model_NElem(model, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED .AND. nelem == 0, 'fail:nelem-uninit')
    coupled = 99
    CALL CD_Get_Model_CoupledDofs(model, coupled, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED .AND. MAXVAL(ABS(coupled)) == 0, 'fail:coupled-map-uninit')
    q = 123.0_wp
    v = 456.0_wp
    a = 789.0_wp
    CALL CD_Get_Model_CoupledMotion(model, q(1:6), v(1:6), a(1:6), es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED .AND. nan_max_abs(q(1:6)) < 1.0e-12_wp .AND. &
                 nan_max_abs(v(1:6)) < 1.0e-12_wp .AND. nan_max_abs(a(1:6)) < 1.0e-12_wp, &
                 'fail:coupled-motion-uninit')
    tension = 123.0_wp
    CALL CD_Get_Model_Tension(model, tension, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED .AND. nan_max_abs(tension) < 1.0e-12_wp, 'fail:tension-uninit')
    CALL CD_Calc_Model_CoupledLoads(model, q(1:6), es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED, 'fail:loads-uninit')
    CALL CD_Step_Model(model, 0.01_wp, conv, stalled, n_iter, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED, 'fail:step-uninit')

    CALL straight_cable(NE, conn, l0, ea, rho_a, q0, fixed)
    v0 = 0.0_wp
    f_ext = 0.0_wp
    l0(1) = -1.0_wp
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'fail:bad-l0')

    q = 123.0_wp
    v = 456.0_wp
    a = 789.0_wp
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_NOT_INITIALIZED .AND. nan_max_abs(q) < 1.0e-12_wp .AND. &
                 nan_max_abs(v) < 1.0e-12_wp .AND. nan_max_abs(a) < 1.0e-12_wp, 'fail:get-uninit')

    l0(1) = 1.0_wp
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                       seabed_z_floor=0.0_wp)
    CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'fail:partial-seabed')

    BLOCK
      REAL(wp) :: ba_bad(1)
      ba_bad = [-1.0_wp]
      CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, ba=ba_bad)
      CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'fail:negative-ba')
    END BLOCK
    BLOCK
      REAL(wp) :: fluid_bad(3, 2)
      fluid_bad = 0.0_wp
      CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                         fluid_velocity=fluid_bad)
      CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'fail:partial-drag')
    END BLOCK
    BLOCK
      REAL(wp) :: fluid_ok(3, 2), waterline(2), diam(1), cdn(1)
      fluid_ok = 0.0_wp
      waterline = 0.0_wp
      diam = [0.2_wp]
      cdn = [1.2_wp]
      CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                         fluid_velocity=fluid_ok, drag_waterline_z=waterline, drag_rho=1000.0_wp, &
                         drag_diameter_elem=diam, drag_cdn_elem=cdn)
      CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'fail:partial-drag-elem')
    END BLOCK
    BLOCK
      REAL(wp) :: waterline(2)
      waterline = 0.0_wp
      CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                         added_mass_waterline_z=waterline, added_mass_rho=1025.0_wp, &
                         added_mass_diameter=0.2_wp, added_mass_can=-1.0_wp, added_mass_cat=0.0_wp)
      CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'fail:negative-added-mass')
    END BLOCK
    BLOCK
      REAL(wp) :: accel_bad(3, 2)
      accel_bad = 0.0_wp
      CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                         fluid_acceleration=accel_bad)
      CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'fail:partial-fk')
    END BLOCK
    BLOCK
      REAL(wp) :: waterline(2)
      waterline = 0.0_wp
      CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em, &
                         buoyancy_waterline_z=waterline, buoyancy_rho=1025.0_wp, &
                         buoyancy_diameter=0.2_wp, buoyancy_gravity=-9.81_wp)
      CALL require(es == CD_MODEL_BADINPUT .AND. .NOT. CD_Model_Is_Initialized(model), 'fail:negative-buoyancy-g')
    END BLOCK

    l0(1) = 1.0_wp
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_OK .AND. CD_Model_Is_Initialized(model), 'fail:valid-before-bad-reinit')
    CALL CD_Get_Model_State(model, q_save, v, a, es, em)
    CALL require(es == CD_MODEL_OK, 'fail:valid-before-bad-reinit-state')
    l0(1) = -1.0_wp
    CALL CD_Init_Model(model, q0, v0, conn, l0, ea, rho_a, .FALSE., f_ext, fixed, cfg, es, em)
    CALL require(es == CD_MODEL_BADINPUT .AND. CD_Model_Is_Initialized(model), 'fail:bad-reinit-preserves-model')
    CALL CD_Get_Model_State(model, q, v, a, es, em)
    CALL require(es == CD_MODEL_OK .AND. nan_max_abs(q - q_save) < 1.0e-12_wp, 'fail:bad-reinit-preserves-state')
    CALL CD_End_Model(model, es, em)
    CALL require(es == CD_MODEL_OK .AND. .NOT. CD_Model_Is_Initialized(model), 'fail:bad-reinit-cleanup')
  END SUBROUTINE case_fail_closed

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. cond) THEN
      WRITE (*, '(A,A,A)') 'MISMATCH [', label, ']: condition false'
      nfail = nfail + 1
    END IF
  END SUBROUTINE require

END PROGRAM test_model
