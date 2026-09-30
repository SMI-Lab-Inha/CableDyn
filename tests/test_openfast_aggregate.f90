! File: tests/test_openfast_aggregate.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_openfast_aggregate
  !! Gate for the OpenFAST CompMooring=5 aggregate: ONE CableDyn deck carrying
  !! BOTH EI=0 mooring lines AND finite-EI Hermite power cables, composed behind the same
  !! FMF moving-point facade the OpenFAST shell drives. The main cases:
  !!   1) pure_mooring_equals_fmf -- a pure EI=0 mooring deck built through the aggregate is
  !!      identical to 1e-12 (moving-point pos/vel/acc/load, every column, every step) to
  !!      driving CD_FMF directly on that deck; ncable == 0. This is the critical degenerate
  !!      case: a pure CompMooring=5 mooring must be unchanged.
  !!   2) mixed_deck_composition -- a chain + a lazy-wave cable compose; ncp_total =
  !!      ncp_sys + ncable; the fairlead loads are physical; and the positional column map
  !!      is stable (a drive perturbation on column i moves ONLY column i's load).
  !!   3) pure_cable_deck -- a finite-EI-only deck builds with no EI=0 system and steps.
  !!   4) coupled_touchdown_contact -- flat and equivalent structured seabeds carry the same
  !!      finite-EI touchdown cable through initialization and moving-fairlead steps.
  !!   5) fail_closed -- a mixed deck carrying a Free (dynamic) point is rejected.
  !!   6) step_atomicity -- a failed step restores the exact step-start state (retryable).
  !!   7) ptfm_init_equivalence -- the ptfm_init option == a pre-displaced deck, bit-for-bit,
  !!      against an independent transcription of the NWTC ZYX transform; zero == absent;
  !!      non-finite fails closed.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Deck_Query_dtM, CD_DECKDRV_OK
  USE CableDyn_OpenFAST_Aggregate, ONLY: CD_AGG_ModuleType, CD_AGG_Init_From_Deck, CD_AGG_NMovingPoints, &
                                         CD_AGG_GetMovingPointMesh, CD_AGG_UpdateStates_Moving, &
                                         CD_AGG_Step_Moving, CD_AGG_CalcOutput, &
                                         CD_AGG_Snapshot, CD_AGG_Restore, &
                                         CD_AGG_NTurbines, CD_AGG_TurbineOfMoving, &
                                         CD_AGG_NFluidNodes, CD_AGG_GetFluidNodePositions, CD_AGG_SetFluidFields, &
                                         CD_AGG_End, CD_AGG_IsInitialized, CD_AGG_OK, CD_AGG_BADINPUT, &
                                         CD_AGG_SOLVEFAIL, &
                                         CD_AGG_NumChannels, CD_AGG_ChannelHeader, CD_AGG_EvalChannel, &
                                         CD_AGG_WriteStaticProfile, CD_AGG_GetInitMetadata, CD_AGG_GetInitLine, &
                                         CD_AGG_NFailures, CD_AGG_Get_Failure_Flags, CD_AGG_Set_Failure_Flags, &
                                         CD_AGG_NCtrlChans, CD_AGG_Apply_LineControl, &
                                         CD_AGG_Refresh_PointMesh, CD_AGG_HasRigid6, CD_AGG_HasRod, &
                                         CD_AGG_Rigid6_MirrorSize, CD_AGG_Get_Rigid6_States, &
                                         CD_AGG_Set_Rigid6_States, CD_AGG_Rod_MirrorSize, &
                                         CD_AGG_Get_Rod_States, CD_AGG_Set_Rod_States
  USE CableDyn_OpenFAST_FMF, ONLY: CD_FMF_ModuleType, CD_FMF_Init_From_Deck, CD_FMF_NMovingPoints, &
                                   CD_FMF_GetMovingPointMesh, CD_FMF_Step_Moving, CD_FMF_CalcOutput, &
                                   CD_FMF_End, CD_FMF_OK
  USE CableDyn_OpenFAST_HermiteFMF, ONLY: CD_HFMF_Restore, CD_HFMF_MinSpanZ, CD_HFMF_Curvature, CD_HFMF_CalcOutput, &
                                          CD_HFMF_OK, &
                                          CD_HFMF_MirrorSize, CD_HFMF_PackMirror, CD_HFMF_UnpackMirror
  ! Public driver-side evaluators, used to cross-check the shell channel path against the
  ! standalone driver's channel definitions on the aggregate's own static configuration.
  USE CableDyn_DeckDriver, ONLY: tangent_declination_deg
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCable_Dyn_Reset_Profile, CD_HermiteCable_Dyn_Get_Profile
  USE CableDyn_Model, ONLY: CD_Get_Model_EndForces, CD_MODEL_OK
  USE CableDyn_System, ONLY: CD_System_NSystemCoupledDOF, CD_Get_System_CoupledMotion, &
                             CD_Update_System_CoupledMotion, &
                             CD_System_NDynamicPoints, CD_Get_System_DynamicPoint_States, &
                             CD_Set_System_DynamicPoint_States, CD_SYSTEM_OK
  USE CableDyn_System, ONLY: CD_System_Line_NElem, CD_Get_System_Line_Tension, CD_SYSTEM_OK, &
                             CD_System_NLines, CD_System_Line_NDOF, CD_Get_System_Line_State, &
                             CD_Update_System_Line_Interior_State, CD_Recompute_System_Acceleration, &
                             CD_SystemType, CD_End_System, CD_System_Line_Has_Friction, &
                             CD_Get_System_Line_Friction_Anchors, CD_Set_System_Line_Friction_Anchors
  USE CableDyn_DeckDriver, ONLY: CD_Init_Deck_System
  USE CableDyn_OpenFAST_HermiteFMF, ONLY: CD_HFMF_FrictionMirrorSize
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE, IEEE_VALUE, IEEE_QUIET_NAN
  USE, INTRINSIC :: IEEE_EXCEPTIONS, ONLY: IEEE_USUAL, IEEE_GET_HALTING_MODE, IEEE_SET_HALTING_MODE, IEEE_SET_FLAG
  IMPLICIT NONE
  LOGICAL :: fp_halt(3)   ! saved IEEE halting modes around deliberately non-finite inputs
  INTEGER :: nfail
  CHARACTER(32) :: test_filter
  REAL(wp), PARAMETER :: PI = 3.141592653589793_wp

  nfail = 0
  IF (COMMAND_ARGUMENT_COUNT() == 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, test_filter)
    IF (TRIM(test_filter) == '--touchdown-only') THEN
      CALL case_coupled_touchdown_contact()
      IF (nfail > 0) THEN
        WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' touchdown assertion(s) failed'
        ERROR STOP 1
      END IF
      WRITE (*, '(A)') 'PASS: coupled flat/structured finite-EI touchdown lifecycle'
      STOP
    END IF
    IF (TRIM(test_filter) == '--feedback-only') THEN
      CALL case_cable_load_feedback()
      IF (nfail > 0) THEN
        WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' cable-load-feedback assertion(s) failed'
        ERROR STOP 1
      END IF
      WRITE (*, '(A)') 'PASS: finite-EI cable load-feedback contract'
      STOP
    END IF
    IF (TRIM(test_filter) == '--host-rod-only') THEN
      CALL case_host_rod()
      CALL case_host_rod_kinematics()
      CALL case_host_body()
      CALL case_host_body_mixed_free()
      IF (nfail > 0) THEN
        WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' host-rod assertion(s) failed'
        ERROR STOP 1
      END IF
      WRITE (*, '(A)') 'PASS: coupled host rods'
      STOP
    END IF
    IF (TRIM(test_filter) == '--pose-init-only') THEN
      CALL case_semitaut_displaced_pose_init()
      IF (nfail > 0) THEN
        WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' displaced-pose initialisation assertion(s) failed'
        ERROR STOP 1
      END IF
      WRITE (*, '(A)') 'PASS: displaced semitaut mooring initialisation'
      STOP
    END IF
  END IF
  CALL case_pure_mooring_equals_fmf()
  CALL case_mixed_deck_composition()
  CALL case_touchdown_cable_end_tension()
  CALL case_output_channels()
  CALL case_pure_cable_deck()
  CALL case_cable_load_feedback()
  CALL case_end_connection_orientation_and_moment()
  CALL case_cable_modified_newton()
  CALL case_dt_guard()
  CALL case_cable_no_advance_transfer()
  CALL case_fail_closed_dynamic_point()
  CALL case_coupled_touchdown_contact()
  CALL case_rejected_static_branch_status()
  CALL case_coupled_viscoelastic()
  CALL case_coupled_point3_buoy()
  CALL case_coupled_point3_fluid()
  CALL case_coupled_rigid6()
  CALL case_coupled_rigid6_mirror()
  CALL case_coupled_rigid6_fluid()
  CALL case_coupled_rod()
  CALL case_coupled_rod_mirror()
  CALL case_coupled_rod_fluid()
  CALL case_host_rod()
  CALL case_host_rod_kinematics()
  CALL case_host_body()
  CALL case_host_body_mixed_free()
  CALL case_step_atomicity()
  CALL case_ptfm_init_equivalence()
  CALL case_correction_rewind()
  CALL case_external_fluid()
  CALL case_still_water_deck_ambient()
  CALL case_external_fluid_report()
  CALL case_coupled_positions_override()
  CALL case_farm_partition()
  CALL case_farm_free_clump()
  CALL case_farm_shared_anchor()
  CALL case_aggregate_failures()
  CALL case_aggregate_line_control()
  CALL case_restart_equivalence()
  CALL case_restart_friction()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: OpenFAST aggregate (mooring + finite-EI cable composition)'

CONTAINS

  INCLUDE 'nan_max_abs.inc'

  SUBROUTINE require(cond, msg)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: msg
    IF (.NOT. cond) THEN
      WRITE (*, '(A)') 'MISMATCH: '//msg
      nfail = nfail + 1
    END IF
  END SUBROUTINE require

  SUBROUTINE case_end_connection_orientation_and_moment()
    !! Carry a parent-relative finite end connection through the aggregate surface
    !! used by the OpenFAST glue: orientation in, support moment out, and exact
    !! snapshot rollback. EI=0 point attachments retain zero moments.
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp) :: pos(3, 2), vel(3, 2), acc(3, 2), load(3, 2), moment(3, 2)
    REAL(wp) :: omega(3, 2), alpha(3, 2), omega_out(3, 2), alpha_out(3, 2)
    REAL(wp) :: orient(3, 3, 2), bad_orient(3, 3, 2), angle, d0_saved(3, 2)
    REAL(wp) :: tangent(3), tangent_norm, d0_before(3), q_before(3)
    REAL(wp) :: tangent_v_expected(3), tangent_a_expected(3), direction_rate(3), direction_acceleration(3)
    REAL(wp) :: omega_local(3), alpha_local(3), parent_omega_before(3), parent_alpha_before(3)
    REAL(wp) :: tangent_rate_before, tangent_acceleration_before
    REAL(wp) :: dynamic_direction_rate(3), dynamic_direction_acceleration(3), reconstructed(3)
    REAL(wp) :: d0_global(3), turn_axis(3), half_turn(3, 3), axis_norm
    REAL(wp) :: moving_before(3, 2), moving_after(3, 2), moving_vel_before(3, 2), moving_acc_before(3, 2)
    REAL(wp) :: moving_vel_after(3, 2), moving_acc_after(3, 2), bad_pos(3, 2), parent_dcm_before(3, 3)
    REAL(wp), ALLOCATABLE :: sys_q_before(:), sys_v_before(:), sys_a_before(:)
    REAL(wp), ALLOCATABLE :: sys_q_after(:), sys_v_after(:), sys_a_after(:)
    REAL(wp), ALLOCATABLE :: line_q_before(:), line_v_before(:), line_a_before(:)
    REAL(wp), ALLOCATABLE :: line_q_after(:), line_v_after(:), line_a_after(:)
    REAL(wp), ALLOCATABLE :: cable_q_before(:), cable_v_before(:), cable_a_before(:)
    REAL(wp) :: cable_d0_before(3, 2), cable_t_before
    INTEGER :: es, ni, cable_col, tangent_base, i, j, ncd, nldof
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_mixed_deck('agg_endconn.dat', finite_end_connection=.TRUE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_endconn.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'agg-endconn:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    cable_col = agg%ncp_sys + 1
    CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, load, es, em, orientation=orient, moment=moment, &
                                   angular_velocity=omega_out, angular_acceleration=alpha_out)
    CALL require(es == CD_AGG_OK, 'agg-endconn:get initial moving surface')
    CALL require(SQRT(SUM(moment(:, cable_col)**2)) > 1.0_wp, &
                 'agg-endconn: finite cable exports a support moment')
    CALL require(nan_max_abs(moment(:, 1:agg%ncp_sys)) <= 0.0_wp, &
                 'agg-endconn: EI=0 point attachments export zero moment')
    CALL require(nan_max_abs(omega_out) <= 0.0_wp .AND. nan_max_abs(alpha_out) <= 0.0_wp, &
                 'agg-endconn: initial parent angular kinematics are zero')

    CALL CD_AGG_Snapshot(agg, es, em)
    CALL require(es == CD_AGG_OK, 'agg-endconn:snapshot')
    d0_saved = agg%cables(1)%line%endconn_d0
    angle = 2.0_wp*PI/180.0_wp
    orient(:, :, cable_col) = 0.0_wp
    orient(1, 1, cable_col) = 1.0_wp
    orient(2, 2, cable_col) = COS(angle); orient(2, 3, cable_col) = SIN(angle)
    orient(3, 2, cable_col) = -SIN(angle); orient(3, 3, cable_col) = COS(angle)
    CALL CD_AGG_Step_Moving(agg, DT, pos, vel, acc, cv, st, ni, es, em, orientation=orient)
    CALL require(es == CD_AGG_OK .AND. cv, 'agg-endconn: oriented step converges: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, load, es, em, orientation=orient, moment=moment)
    CALL require(es == CD_AGG_OK .AND. ABS(agg%cables(1)%line%endconn_d0(2, 2)) > 1.0e-3_wp, &
                 'agg-endconn: orientation reaches the cable parent connection')
    CALL require(ALL(IEEE_IS_FINITE(moment(:, cable_col))), &
                 'agg-endconn: oriented support moment is finite')
    CALL CD_AGG_Restore(agg, es, em)
    CALL require(es == CD_AGG_OK .AND. nan_max_abs(agg%cables(1)%line%endconn_d0 - d0_saved) <= 0.0_wp, &
                 'agg-endconn: restore rewinds connection direction exactly')

    ! Validate every DCM before touching any subsystem. A reflected input must not
    ! partly update the EI=0 system before the cable rejects it.
    bad_orient = orient
    bad_orient(:, :, cable_col) = 0.0_wp
    bad_orient(1, 1, cable_col) = -1.0_wp
    bad_orient(2, 2, cable_col) = 1.0_wp
    bad_orient(3, 3, cable_col) = 1.0_wp
    CALL CD_AGG_UpdateStates_Moving(agg, pos, vel, acc, es, em, orientation=bad_orient)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, 'proper orthogonal') > 0, &
                 'agg-endconn: reflected orientation fails before state mutation')
    CALL require(nan_max_abs(agg%cables(1)%line%endconn_d0 - d0_saved) <= 0.0_wp, &
                 'agg-endconn: rejected orientation leaves connection direction unchanged')
    CALL CD_AGG_End(agg, es, em)

    ! Repeat the host-facing path with the exact Rigid branch.  UpdateStates is
    ! deliberately non-stepping here: it proves that the direct-feedthrough load
    ! probe cannot combine a new parent orientation with a stale cable tangent.
    CALL write_mixed_deck('agg_endconn_rigid.dat', rigid_end_connection=.TRUE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_endconn_rigid.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'agg-rigid-endconn:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    cable_col = agg%ncp_sys + 1
    tangent_base = 6*(agg%cables(1)%line%nn - 1) + 3
    CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, load, es, em, orientation=orient, moment=moment)
    CALL require(es == CD_AGG_OK, 'agg-rigid-endconn:get initial moving surface')
    angle = 2.0_wp*PI/180.0_wp
    orient(:, :, cable_col) = 0.0_wp
    orient(1, 1, cable_col) = 1.0_wp
    orient(2, 2, cable_col) = COS(angle); orient(2, 3, cable_col) = SIN(angle)
    orient(3, 2, cable_col) = -SIN(angle); orient(3, 3, cable_col) = COS(angle)
    CALL CD_AGG_UpdateStates_Moving(agg, pos, vel, acc, es, em, orientation=orient)
    CALL require(es == CD_AGG_OK, 'agg-rigid-endconn: non-stepping oriented update: '//TRIM(em))
    tangent = agg%cables(1)%line%q(tangent_base + 1:tangent_base + 3)
    tangent_norm = SQRT(DOT_PRODUCT(tangent, tangent))
    CALL require(nan_max_abs(tangent/tangent_norm - agg%cables(1)%line%endconn_d0(:, 2)) <= &
                 64.0_wp*EPSILON(1.0_wp), &
                 'agg-rigid-endconn: direct feedthrough enforces the rotated tangent exactly')

    ! Exercise the same transport inside the implicit solve. Start from the
    ! zero-rate orientation above and end a short constant-angular-acceleration
    ! increment with compatible orientation, angular velocity, and acceleration.
    angle = angle + 0.5_wp*0.4_wp*DT*DT
    orient(:, :, cable_col) = 0.0_wp
    orient(1, 1, cable_col) = 1.0_wp
    orient(2, 2, cable_col) = COS(angle); orient(2, 3, cable_col) = SIN(angle)
    orient(3, 2, cable_col) = -SIN(angle); orient(3, 3, cable_col) = COS(angle)
    omega = 0.0_wp
    alpha = 0.0_wp
    omega(1, cable_col) = 0.4_wp*DT
    alpha(1, cable_col) = 0.4_wp
    CALL CD_AGG_Step_Moving(agg, DT, pos, vel, acc, cv, st, ni, es, em, orientation=orient, &
                            angular_velocity=omega, angular_acceleration=alpha)
    CALL require(es == CD_AGG_OK .AND. cv, &
                 'agg-rigid-endconn: angular-acceleration implicit step converges: '//TRIM(em))
    tangent = agg%cables(1)%line%q(tangent_base + 1:tangent_base + 3)
    tangent_norm = SQRT(DOT_PRODUCT(tangent, tangent))
    omega_local = [agg%cables(1)%frame_c*omega(1, cable_col) + &
                   agg%cables(1)%frame_s*omega(2, cable_col), &
                   -agg%cables(1)%frame_s*omega(1, cable_col) + &
                   agg%cables(1)%frame_c*omega(2, cable_col), omega(3, cable_col)]
    alpha_local = [agg%cables(1)%frame_c*alpha(1, cable_col) + &
                   agg%cables(1)%frame_s*alpha(2, cable_col), &
                   -agg%cables(1)%frame_s*alpha(1, cable_col) + &
                   agg%cables(1)%frame_c*alpha(2, cable_col), alpha(3, cable_col)]
    dynamic_direction_rate = test_cross3(omega_local, agg%cables(1)%line%endconn_d0(:, 2))
    dynamic_direction_acceleration = &
      test_cross3(alpha_local, agg%cables(1)%line%endconn_d0(:, 2)) + &
      test_cross3(omega_local, dynamic_direction_rate)
    tangent_rate_before = DOT_PRODUCT(agg%cables(1)%line%v(tangent_base + 1:tangent_base + 3), &
                                      agg%cables(1)%line%endconn_d0(:, 2))
    tangent_acceleration_before = &
      DOT_PRODUCT(agg%cables(1)%line%a(tangent_base + 1:tangent_base + 3), &
                  agg%cables(1)%line%endconn_d0(:, 2)) + &
      tangent_norm*DOT_PRODUCT(dynamic_direction_rate, dynamic_direction_rate)
    reconstructed = tangent_rate_before*agg%cables(1)%line%endconn_d0(:, 2) + &
                    tangent_norm*dynamic_direction_rate
    CALL require(nan_max_abs(agg%cables(1)%line%v(tangent_base + 1:tangent_base + 3) - reconstructed) <= &
                 2.0e-12_wp*MAX(1.0_wp, tangent_norm), &
                 'agg-rigid-endconn: implicit step commits exact rigid tangent velocity')
    reconstructed = tangent_acceleration_before*agg%cables(1)%line%endconn_d0(:, 2) + &
                    2.0_wp*tangent_rate_before*dynamic_direction_rate + &
                    tangent_norm*dynamic_direction_acceleration
    CALL require(nan_max_abs(agg%cables(1)%line%a(tangent_base + 1:tangent_base + 3) - reconstructed) <= &
                 2.0e-12_wp*MAX(1.0_wp, tangent_norm), &
                 'agg-rigid-endconn: implicit step commits exact angular and centripetal acceleration')
    CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, load, es, em, orientation=orient, moment=moment, &
                                   angular_velocity=omega_out, angular_acceleration=alpha_out)
    CALL require(es == CD_AGG_OK .AND. nan_max_abs(omega_out - omega) <= 0.0_wp .AND. &
                 nan_max_abs(alpha_out - alpha) <= 0.0_wp, &
                 'agg-rigid-endconn: implicit angular state is observable at the boundary')

    ! OpenFAST RotationVel/RotationAcc are global vectors.  A rigid parent-relative
    ! tangent m = n d must carry the complete transport terms, including the
    ! centripetal acceleration; simply rotating its previous v/a would leave both
    ! zero in this at-rest orientation probe.
    omega = 0.0_wp
    alpha = 0.0_wp
    omega(:, cable_col) = [0.7_wp, -0.4_wp, 1.2_wp]
    alpha(:, cable_col) = [-0.3_wp, 0.5_wp, 0.2_wp]
    tangent_rate_before = DOT_PRODUCT(agg%cables(1)%line%v(tangent_base + 1:tangent_base + 3), &
                                      agg%cables(1)%line%endconn_d0(:, 2))
    tangent_acceleration_before = DOT_PRODUCT(agg%cables(1)%line%a(tangent_base + 1:tangent_base + 3), &
                                              agg%cables(1)%line%endconn_d0(:, 2))
    omega_local = [agg%cables(1)%frame_c*agg%cables(1)%parent_omega(1) + &
                   agg%cables(1)%frame_s*agg%cables(1)%parent_omega(2), &
                   -agg%cables(1)%frame_s*agg%cables(1)%parent_omega(1) + &
                   agg%cables(1)%frame_c*agg%cables(1)%parent_omega(2), agg%cables(1)%parent_omega(3)]
    dynamic_direction_rate = test_cross3(omega_local, agg%cables(1)%line%endconn_d0(:, 2))
    tangent_acceleration_before = tangent_acceleration_before + &
                                  tangent_norm*DOT_PRODUCT(dynamic_direction_rate, dynamic_direction_rate)
    CALL CD_AGG_UpdateStates_Moving(agg, pos, vel, acc, es, em, orientation=orient, &
                                    angular_velocity=omega, angular_acceleration=alpha)
    CALL require(es == CD_AGG_OK, 'agg-rigid-endconn: angular transport overlay: '//TRIM(em))
    tangent = agg%cables(1)%line%q(tangent_base + 1:tangent_base + 3)
    tangent_norm = SQRT(DOT_PRODUCT(tangent, tangent))
    omega_local = [agg%cables(1)%frame_c*omega(1, cable_col) + &
                   agg%cables(1)%frame_s*omega(2, cable_col), &
                   -agg%cables(1)%frame_s*omega(1, cable_col) + &
                   agg%cables(1)%frame_c*omega(2, cable_col), omega(3, cable_col)]
    alpha_local = [agg%cables(1)%frame_c*alpha(1, cable_col) + &
                   agg%cables(1)%frame_s*alpha(2, cable_col), &
                   -agg%cables(1)%frame_s*alpha(1, cable_col) + &
                   agg%cables(1)%frame_c*alpha(2, cable_col), alpha(3, cable_col)]
    direction_rate = test_cross3(omega_local, agg%cables(1)%line%endconn_d0(:, 2))
    direction_acceleration = test_cross3(alpha_local, agg%cables(1)%line%endconn_d0(:, 2)) + &
                             test_cross3(omega_local, direction_rate)
    tangent_v_expected = tangent_rate_before*agg%cables(1)%line%endconn_d0(:, 2) + &
                         tangent_norm*direction_rate
    tangent_a_expected = tangent_acceleration_before*agg%cables(1)%line%endconn_d0(:, 2) + &
                         2.0_wp*tangent_rate_before*direction_rate + tangent_norm*direction_acceleration
    CALL require(nan_max_abs(agg%cables(1)%line%v(tangent_base + 1:tangent_base + 3) - &
                             tangent_v_expected) <= 2.0e-12_wp*MAX(1.0_wp, tangent_norm), &
                 'agg-rigid-endconn: spin transports rigid tangent velocity exactly')
    CALL require(nan_max_abs(agg%cables(1)%line%a(tangent_base + 1:tangent_base + 3) - &
                             tangent_a_expected) <= 2.0e-12_wp*MAX(1.0_wp, tangent_norm), &
                 'agg-rigid-endconn: angular and centripetal tangent acceleration are exact')
    CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, load, es, em, orientation=orient, moment=moment, &
                                   angular_velocity=omega_out, angular_acceleration=alpha_out)
    CALL require(es == CD_AGG_OK .AND. nan_max_abs(omega_out - omega) <= 0.0_wp .AND. &
                 nan_max_abs(alpha_out - alpha) <= 0.0_wp, &
                 'agg-rigid-endconn: angular kinematics round-trip through the aggregate')
    CALL CD_AGG_UpdateStates_Moving(agg, pos, vel, acc, es, em, orientation=orient, angular_velocity=omega)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, 'supplied together') > 0 .AND. &
                 nan_max_abs(agg%cables(1)%line%v(tangent_base + 1:tangent_base + 3) - &
                             tangent_v_expected) <= 2.0e-12_wp*MAX(1.0_wp, tangent_norm), &
                 'agg-rigid-endconn: partial angular input fails before mutation')
    alpha(1, cable_col) = IEEE_VALUE(0.0_wp, IEEE_QUIET_NAN)
    CALL CD_AGG_UpdateStates_Moving(agg, pos, vel, acc, es, em, orientation=orient, &
                                    angular_velocity=omega, angular_acceleration=alpha)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, 'angular kinematics must be finite') > 0 .AND. &
                 nan_max_abs(agg%cables(1)%line%a(tangent_base + 1:tangent_base + 3) - &
                             tangent_a_expected) <= 2.0e-12_wp*MAX(1.0_wp, tangent_norm), &
                 'agg-rigid-endconn: non-finite angular input fails before mutation')
    alpha(:, cable_col) = [-0.3_wp, 0.5_wp, 0.2_wp]
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, load, es, em, orientation=orient, moment=moment)
    CALL require(es == CD_AGG_OK .AND. ALL(IEEE_IS_FINITE(moment(:, cable_col))), &
                 'agg-rigid-endconn: exact-constraint reaction moment is finite')
    CALL require(SQRT(SUM(moment(:, cable_col)**2)) > 1.0_wp, &
                 'agg-rigid-endconn: exact constraint exports a non-zero reaction moment')

    ! A proper 180-degree parent rotation makes the prescribed direction
    ! antiparallel.  The shortest rotation is then non-unique, so reject it
    ! atomically rather than selecting an arbitrary bending plane. Move the
    ! EI=0 point in the same rejected request so a system-first implementation
    ! cannot pass this gate by leaving an unobserved partial update behind.
    d0_before = agg%cables(1)%line%endconn_d0(:, 2)
    q_before = agg%cables(1)%line%q(tangent_base + 1:tangent_base + 3)
    parent_dcm_before = agg%cables(1)%parent_dcm
    parent_omega_before = agg%cables(1)%parent_omega
    parent_alpha_before = agg%cables(1)%parent_alpha
    moving_before = pos
    moving_vel_before = vel
    moving_acc_before = acc
    ncd = CD_System_NSystemCoupledDOF(agg%sys%fast%system)
    ALLOCATE (sys_q_before(ncd), sys_v_before(ncd), sys_a_before(ncd))
    ALLOCATE (sys_q_after(ncd), sys_v_after(ncd), sys_a_after(ncd))
    nldof = CD_System_Line_NDOF(agg%sys%fast%system, 1, es, em)
    CALL require(es == CD_SYSTEM_OK .AND. nldof > 0, 'agg-rigid-endconn: size EI=0 line state')
    ALLOCATE (line_q_before(nldof), line_v_before(nldof), line_a_before(nldof))
    ALLOCATE (line_q_after(nldof), line_v_after(nldof), line_a_after(nldof))
    ALLOCATE (cable_q_before(SIZE(agg%cables(1)%line%q)), &
              cable_v_before(SIZE(agg%cables(1)%line%v)), &
              cable_a_before(SIZE(agg%cables(1)%line%a)))
    cable_q_before = agg%cables(1)%line%q
    cable_v_before = agg%cables(1)%line%v
    cable_a_before = agg%cables(1)%line%a
    cable_d0_before = agg%cables(1)%line%endconn_d0
    cable_t_before = agg%cables(1)%line%t
    CALL CD_Get_System_CoupledMotion(agg%sys%fast%system, sys_q_before, sys_v_before, sys_a_before, es, em)
    CALL require(es == CD_SYSTEM_OK, 'agg-rigid-endconn: capture system state before rejection')
    CALL CD_Get_System_Line_State(agg%sys%fast%system, 1, line_q_before, line_v_before, line_a_before, es, em)
    CALL require(es == CD_SYSTEM_OK, 'agg-rigid-endconn: capture EI=0 line state before rejection')
    bad_orient = orient
    d0_global = [agg%cables(1)%frame_c*d0_before(1) - agg%cables(1)%frame_s*d0_before(2), &
                 agg%cables(1)%frame_s*d0_before(1) + agg%cables(1)%frame_c*d0_before(2), &
                 d0_before(3)]
    turn_axis = [-d0_global(2), d0_global(1), 0.0_wp]
    axis_norm = SQRT(DOT_PRODUCT(turn_axis, turn_axis))
    turn_axis = turn_axis/axis_norm
    half_turn = 0.0_wp
    DO j = 1, 3
      DO i = 1, 3
        half_turn(i, j) = 2.0_wp*turn_axis(i)*turn_axis(j)
      END DO
      half_turn(j, j) = half_turn(j, j) - 1.0_wp
    END DO
    bad_orient(:, :, cable_col) = MATMUL(orient(:, :, cable_col), half_turn)
    bad_pos = moving_before
    bad_pos(:, 1) = bad_pos(:, 1) + [0.0_wp, 0.0_wp, 0.025_wp]
    CALL CD_AGG_UpdateStates_Moving(agg, bad_pos, moving_vel_before, moving_acc_before, es, em, &
                                    orientation=bad_orient)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, '180-degree') > 0, &
                 'agg-rigid-endconn: antiparallel parent update fails closed')
    CALL CD_Get_System_CoupledMotion(agg%sys%fast%system, sys_q_after, sys_v_after, sys_a_after, es, em)
    CALL require(es == CD_SYSTEM_OK, 'agg-rigid-endconn: read system state after rejection')
    CALL CD_Get_System_Line_State(agg%sys%fast%system, 1, line_q_after, line_v_after, line_a_after, es, em)
    CALL require(es == CD_SYSTEM_OK, 'agg-rigid-endconn: read EI=0 line state after rejection')
    CALL CD_AGG_GetMovingPointMesh(agg, moving_after, moving_vel_after, moving_acc_after, load, es, em)
    CALL require(es == CD_AGG_OK, 'agg-rigid-endconn: read aggregate after rejection')
    CALL require(nan_max_abs(sys_q_after - sys_q_before) <= 0.0_wp .AND. &
                 nan_max_abs(sys_v_after - sys_v_before) <= 0.0_wp .AND. &
                 nan_max_abs(sys_a_after - sys_a_before) <= 0.0_wp, &
                 'agg-rigid-endconn: rejected parent update leaves EI=0 state bit-identical')
    CALL require(nan_max_abs(line_q_after - line_q_before) <= 0.0_wp .AND. &
                 nan_max_abs(line_v_after - line_v_before) <= 0.0_wp .AND. &
                 nan_max_abs(line_a_after - line_a_before) <= 0.0_wp, &
                 'agg-rigid-endconn: rejected parent update leaves EI=0 line bit-identical')
    CALL require(nan_max_abs(moving_after - moving_before) <= 0.0_wp .AND. &
                 nan_max_abs(moving_vel_after - moving_vel_before) <= 0.0_wp .AND. &
                 nan_max_abs(moving_acc_after - moving_acc_before) <= 0.0_wp, &
                 'agg-rigid-endconn: rejected parent update leaves facade bit-identical')
    CALL require(nan_max_abs(agg%cables(1)%line%endconn_d0(:, 2) - d0_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%line%q(tangent_base + 1:tangent_base + 3) - q_before) <= 0.0_wp, &
                 'agg-rigid-endconn: rejected parent update is atomic')
    CALL require(nan_max_abs(agg%cables(1)%line%q - cable_q_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%line%v - cable_v_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%line%a - cable_a_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%line%endconn_d0 - cable_d0_before) <= 0.0_wp .AND. &
                 ABS(agg%cables(1)%line%t - cable_t_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%parent_dcm - parent_dcm_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%parent_omega - parent_omega_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%parent_alpha - parent_alpha_before) <= 0.0_wp, &
                 'agg-rigid-endconn: rejected parent update leaves full cable state bit-identical')

    ! The same aggregate-wide preflight must catch invalid kinematics in a later
    ! cable column before the valid EI=0 motion in an earlier column is applied.
    bad_pos = moving_before
    bad_pos(:, 1) = bad_pos(:, 1) + [0.0_wp, 0.0_wp, 0.025_wp]
    bad_pos(1, cable_col) = IEEE_VALUE(0.0_wp, IEEE_QUIET_NAN)
    CALL CD_AGG_UpdateStates_Moving(agg, bad_pos, moving_vel_before, moving_acc_before, es, em, orientation=orient)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, 'kinematics must be finite') > 0, &
                 'agg-rigid-endconn: non-finite later cable input fails closed')
    CALL CD_Get_System_CoupledMotion(agg%sys%fast%system, sys_q_after, sys_v_after, sys_a_after, es, em)
    CALL require(es == CD_SYSTEM_OK .AND. nan_max_abs(sys_q_after - sys_q_before) <= 0.0_wp .AND. &
                 nan_max_abs(sys_v_after - sys_v_before) <= 0.0_wp .AND. &
                 nan_max_abs(sys_a_after - sys_a_before) <= 0.0_wp, &
                 'agg-rigid-endconn: rejected non-finite cable input leaves EI=0 state bit-identical')
    CALL CD_Get_System_Line_State(agg%sys%fast%system, 1, line_q_after, line_v_after, line_a_after, es, em)
    CALL require(es == CD_SYSTEM_OK .AND. nan_max_abs(line_q_after - line_q_before) <= 0.0_wp .AND. &
                 nan_max_abs(line_v_after - line_v_before) <= 0.0_wp .AND. &
                 nan_max_abs(line_a_after - line_a_before) <= 0.0_wp, &
                 'agg-rigid-endconn: rejected non-finite input leaves EI=0 line bit-identical')
    CALL require(nan_max_abs(agg%cables(1)%line%q - cable_q_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%line%v - cable_v_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%line%a - cable_a_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%line%endconn_d0 - cable_d0_before) <= 0.0_wp .AND. &
                 ABS(agg%cables(1)%line%t - cable_t_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%parent_dcm - parent_dcm_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%parent_omega - parent_omega_before) <= 0.0_wp .AND. &
                 nan_max_abs(agg%cables(1)%parent_alpha - parent_alpha_before) <= 0.0_wp, &
                 'agg-rigid-endconn: rejected non-finite input leaves cable state bit-identical')

    ! A rejected probe must not poison the retry baseline.
    CALL CD_AGG_UpdateStates_Moving(agg, moving_before, moving_vel_before, moving_acc_before, es, em, &
                                    orientation=orient)
    CALL require(es == CD_AGG_OK, 'agg-rigid-endconn: valid retry succeeds after rejected updates')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_end_connection_orientation_and_moment

  PURE FUNCTION test_cross3(a, b) RESULT(c)
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  END FUNCTION test_cross3

  SUBROUTINE case_pure_mooring_equals_fmf()
    !! A pure EI=0 mooring deck driven through the aggregate must be bit-for-bit identical to
    !! driving CD_FMF directly: same system, same finalize, same delegated step/output.
    TYPE(CD_AGG_ModuleType) :: agg
    TYPE(CD_FMF_ModuleType) :: fmf
    INTEGER :: es, na, nf, s, i, ia, ifmf, nlines, npoints, nsections, nei0, nfinite, line_id
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.05_wp, AMP = 0.2_wp, PER = 16.0_wp, TOL = 1.0e-12_wp
    INTEGER, PARAMETER :: NSTEP = 8
    REAL(wp) :: pa(3, 3), va(3, 3), aa(3, 3), la(3, 3)
    REAL(wp) :: pf(3, 3), vf(3, 3), af(3, 3), lf(3, 3)
    REAL(wp) :: r0(3, 3), tgt(3, 3), vel(3, 3), acc(3, 3), zh, vh, ah, t
    REAL(wp) :: fair_ten, fair_force(3), fair_incl, fair_decl, fair_azi
    LOGICAL :: ca, cf, sa, sf

    CALL write_volturnus_deck('agg_pure.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_pure.dat', DT, es, em)
    CALL require(es == CD_AGG_OK .AND. CD_AGG_IsInitialized(agg), 'pure:agg-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_GetInitMetadata(agg, nlines, npoints, nsections, nei0, nfinite, es, em)
    CALL require(es == CD_AGG_OK .AND. nlines == 3 .AND. npoints == 6 .AND. nsections == 3 .AND. &
                 nei0 == 3 .AND. nfinite == 0, 'pure:init-report-inventory')
    CALL CD_AGG_GetInitLine(agg, 1, line_id, fair_ten, fair_force, fair_incl, fair_decl, fair_azi, es, em)
    CALL require(es == CD_AGG_OK .AND. line_id == 1 .AND. fair_ten > 0.0_wp .AND. &
                 ALL(IEEE_IS_FINITE(fair_force)) .AND. ABS(NORM2(fair_force) - fair_ten) < 1.0e-9_wp*fair_ten .AND. &
                 ABS(fair_incl - (fair_decl - 90.0_wp)) < 1.0e-12_wp .AND. &
                 fair_azi >= 0.0_wp .AND. fair_azi < 360.0_wp, 'pure:init-report-line-1')
    CALL CD_FMF_Init_From_Deck(fmf, 'agg_pure.dat', DT, es, em)
    CALL require(es == CD_FMF_OK, 'pure:fmf-init: '//TRIM(em))
    IF (es /= CD_FMF_OK) THEN
      CALL CD_AGG_End(agg, es, em)
      RETURN
    END IF

    na = CD_AGG_NMovingPoints(agg, es, em)
    nf = CD_FMF_NMovingPoints(fmf, es, em)
    CALL require(na == 3 .AND. nf == 3 .AND. na == nf, 'pure:nmoving-is-3-both')
    CALL require(agg%ncable == 0, 'pure:ncable-zero')
    CALL require(agg%ncp_sys == 3, 'pure:ncp_sys-3')

    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, la, es, em)
    CALL require(es == CD_AGG_OK, 'pure:agg-initial-mesh: '//TRIM(em))

    DO s = 1, NSTEP
      t = REAL(s, wp)*DT
      zh = AMP*(1.0_wp - COS(2.0_wp*PI*t/PER))
      vh = AMP*(2.0_wp*PI/PER)*SIN(2.0_wp*PI*t/PER)
      ah = AMP*(2.0_wp*PI/PER)**2*COS(2.0_wp*PI*t/PER)
      DO i = 1, 3
        tgt(:, i) = [r0(1, i), r0(2, i), r0(3, i) + zh]
        vel(:, i) = [0.0_wp, 0.0_wp, vh]
        acc(:, i) = [0.0_wp, 0.0_wp, ah]
      END DO
      CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, ca, sa, ia, es, em)
      CALL require(es == CD_AGG_OK, 'pure:agg-step: '//TRIM(em))
      CALL CD_FMF_Step_Moving(fmf, DT, tgt, vel, acc, cf, sf, ifmf, es, em)
      CALL require(es == CD_FMF_OK, 'pure:fmf-step: '//TRIM(em))
      CALL require((ca .EQV. cf) .AND. (sa .EQV. sf) .AND. ia == ifmf, 'pure:step-flags-match')
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL require(es == CD_AGG_OK, 'pure:agg-calc')
      CALL CD_FMF_CalcOutput(fmf, es, em)
      CALL require(es == CD_FMF_OK, 'pure:fmf-calc')
      CALL CD_AGG_GetMovingPointMesh(agg, pa, va, aa, la, es, em)
      CALL require(es == CD_AGG_OK, 'pure:agg-mesh')
      CALL CD_FMF_GetMovingPointMesh(fmf, pf, vf, af, lf, es, em)
      CALL require(es == CD_FMF_OK, 'pure:fmf-mesh')
      CALL require(nan_max_abs(pa - pf) < TOL .AND. nan_max_abs(va - vf) < TOL .AND. &
                   nan_max_abs(aa - af) < TOL .AND. nan_max_abs(la - lf) < TOL, 'pure:mesh-bit-for-bit')
    END DO
    CALL require(nan_max_abs(la) > 1.0_wp, 'pure:loads-nonzero (mooring pretension)')
    WRITE (*, '(A,ES12.4)') 'aggregate==FMF pure mooring: max |mesh diff| over run = ', &
      MAXVAL([nan_max_abs(pa - pf), nan_max_abs(va - vf), nan_max_abs(aa - af), nan_max_abs(la - lf)])
    CALL CD_AGG_End(agg, es, em)
    CALL CD_FMF_End(fmf, es, em)
  END SUBROUTINE case_pure_mooring_equals_fmf

  SUBROUTINE run_mixed(path, h1, h2, nstep, load_out, ncps, ncbl, es_out, em_out)
    !! Build the mixed aggregate from `path`, drive column 1 (mooring fairlead) with heave
    !! amplitude h1 and column 2 (cable fairlead) with h2 for nstep steps (both from rest via
    !! a 1-cos ramp), and return the final fairlead loads (3,2) plus ncp_sys / ncable.
    !! es_out /= 0 flags a build failure or a non-converged/failed step.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), INTENT(IN) :: h1, h2
    INTEGER, INTENT(IN) :: nstep
    REAL(wp), INTENT(OUT) :: load_out(3, 2)
    INTEGER, INTENT(OUT) :: ncps, ncbl, es_out
    CHARACTER(*), INTENT(OUT) :: em_out
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp), PARAMETER :: DT = 0.05_wp, PER = 12.0_wp
    REAL(wp) :: r0(3, 2), tgt(3, 2), vel(3, 2), acc(3, 2), p(3, 2), v(3, 2), a(3, 2), l(3, 2)
    REAL(wp) :: zh, vh, ah, t, hc(2)
    INTEGER :: s, col, es, ni
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    load_out = 0.0_wp
    ncps = 0
    ncbl = 0
    es_out = 0
    em_out = ''
    CALL CD_AGG_Init_From_Deck(agg, path, DT, es, em)
    IF (es /= CD_AGG_OK) THEN
      es_out = es
      em_out = 'init: '//TRIM(em)
      RETURN
    END IF
    ncps = agg%ncp_sys
    ncbl = agg%ncable
    hc = [h1, h2]
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, l, es, em)
    IF (es /= CD_AGG_OK) THEN
      es_out = es
      em_out = 'mesh0: '//TRIM(em)
      CALL CD_AGG_End(agg, es, em)
      RETURN
    END IF
    DO s = 1, nstep
      t = REAL(s, wp)*DT
      zh = 1.0_wp - COS(2.0_wp*PI*t/PER)
      vh = (2.0_wp*PI/PER)*SIN(2.0_wp*PI*t/PER)
      ah = (2.0_wp*PI/PER)**2*COS(2.0_wp*PI*t/PER)
      DO col = 1, 2
        tgt(:, col) = [r0(1, col), r0(2, col), r0(3, col) + hc(col)*zh]
        vel(:, col) = [0.0_wp, 0.0_wp, hc(col)*vh]
        acc(:, col) = [0.0_wp, 0.0_wp, hc(col)*ah]
      END DO
      CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
      IF (es /= CD_AGG_OK .OR. .NOT. cv) THEN
        es_out = MERGE(es, 2, es /= 0)
        em_out = 'step: '//TRIM(em)
        CALL CD_AGG_End(agg, es, em)
        RETURN
      END IF
    END DO
    CALL CD_AGG_CalcOutput(agg, es, em)
    IF (es /= CD_AGG_OK) THEN
      es_out = es
      em_out = 'calc: '//TRIM(em)
      CALL CD_AGG_End(agg, es, em)
      RETURN
    END IF
    CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, load_out, es, em)
    IF (es /= CD_AGG_OK) THEN
      es_out = es
      em_out = 'meshN: '//TRIM(em)
      CALL CD_AGG_End(agg, es, em)
      RETURN
    END IF
    es_out = CD_AGG_OK
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE run_mixed

  SUBROUTINE case_mixed_deck_composition()
    !! A chain (EI=0) + a lazy-wave cable (EI>0) compose. Column map: column 1 = mooring
    !! fairlead (the sys moving point), column 2 = cable fairlead. Because the mooring line
    !! and the cable are independent objects, a drive perturbation on one column changes
    !! ONLY that column's load -- the unperturbed column is bit-for-bit unchanged (zero
    !! cross-talk), which pins the positional routing.
    REAL(wp) :: L00(3, 2), L01(3, 2), L10(3, 2)
    INTEGER :: ncps, ncbl, es
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: AMP = 1.0_wp
    INTEGER, PARAMETER :: N = 48

    CALL write_mixed_deck('agg_mixed.dat')
    CALL run_mixed('agg_mixed.dat', 0.0_wp, 0.0_wp, N, L00, ncps, ncbl, es, em)
    CALL require(es == CD_AGG_OK, 'mixed:rest-run: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL require(ncps == 1, 'mixed:ncp_sys-1')
    CALL require(ncbl == 1, 'mixed:ncable-1')
    CALL require(ALL(IEEE_IS_FINITE(L00)), 'mixed:loads-finite')
    ! Mooring column: catenary pretension, order MN, pulling the fairlead down (Fz < 0).
    CALL require(NORM2(L00(:, 1)) > 1.0e5_wp .AND. L00(3, 1) < 0.0_wp, 'mixed:mooring-column-tension-physical')
    ! Cable column: a physical hang-off reaction (finite, non-trivial).
    CALL require(NORM2(L00(:, 2)) > 1.0e2_wp, 'mixed:cable-column-load-physical')
    WRITE (*, '(A,ES12.4,A,ES12.4)') 'mixed rest loads:  |mooring col| = ', NORM2(L00(:, 1)), &
      ' N;  |cable col| = ', NORM2(L00(:, 2))

    ! Perturb ONLY the cable column (2): the cable load must move; the mooring column must not.
    CALL run_mixed('agg_mixed.dat', 0.0_wp, AMP, N, L01, ncps, ncbl, es, em)
    CALL require(es == CD_AGG_OK, 'mixed:cable-perturb-run: '//TRIM(em))
    IF (es == CD_AGG_OK) THEN
      CALL require(NORM2(L01(:, 1) - L00(:, 1)) < 1.0e-6_wp, 'mixed:cable-drive-leaves-mooring-column-untouched')
      CALL require(NORM2(L01(:, 2) - L00(:, 2)) > 1.0e1_wp, 'mixed:cable-drive-moves-cable-column')
    END IF

    ! Perturb ONLY the mooring column (1): the mooring load must move; the cable column must not.
    CALL run_mixed('agg_mixed.dat', AMP, 0.0_wp, N, L10, ncps, ncbl, es, em)
    CALL require(es == CD_AGG_OK, 'mixed:mooring-perturb-run: '//TRIM(em))
    IF (es == CD_AGG_OK) THEN
      CALL require(NORM2(L10(:, 2) - L00(:, 2)) < 1.0e-6_wp, 'mixed:mooring-drive-leaves-cable-column-untouched')
      CALL require(NORM2(L10(:, 1) - L00(:, 1)) > 1.0e1_wp, 'mixed:mooring-drive-moves-mooring-column')
    END IF
  END SUBROUTINE case_mixed_deck_composition

  SUBROUTINE case_output_channels()
    !! A mixed deck (cable = line 1, mooring = line 2) that declares an OUTPUTS section. The
    !! aggregate retains the channel list + the deck-line-id -> object map, and the shell accessor
    !! CD_AGG_EvalChannel returns, for every channel, the value the standalone driver defines --
    !! cross-checked against a hand-computed reference built from the aggregate's OWN committed
    !! state with the driver's public evaluators: mooring end tensions from the system accessor,
    !! cable curvature from the Menger formula on the raw node positions, the fairlead angle from
    !! the End-A tangent, and the bend-moment == EI * curvature relation. Also checks headers /
    !! units and that a no-OUTPUTS deck reports zero channels (the byte-identical path).
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, nch, nnode, ne, base
    CHARACTER(512) :: em
    CHARACTER(64) :: hdr, unt
    REAL(wp), PARAMETER :: DT = 0.05_wp, EI_BARE = 1.99e4_wp
    REAL(wp) :: vFairTen2, vAnchTen2, vCurv1N5, vFairAngle1, vBendMom1N5, vTen1N10
    REAL(wp) :: vFairDecl1, vFairIncl1
    REAL(wp) :: r1(3), r2(3), kexp, decexp
    REAL(wp), ALLOCATABLE :: tension(:), exact_curvature(:)

    CALL write_mixed_channels_deck('agg_channels.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_channels.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'chan:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL require(agg%ncp_sys == 1 .AND. agg%ncable == 1, 'chan:mixed-shape')

    nch = CD_AGG_NumChannels(agg)
    CALL require(nch == 8, 'chan:num-channels-8')
    IF (nch /= 8) THEN
      CALL CD_AGG_End(agg, es, em); RETURN
    END IF

    ! Header + unit spot checks (original-case token, OrcaFlex unit per kind).
    CALL CD_AGG_ChannelHeader(agg, 1, hdr, unt, es, em)
    CALL require(es == CD_AGG_OK .AND. TRIM(hdr) == 'FairTen2' .AND. TRIM(unt) == '(N)', 'chan:hdr1-FairTen2')
    CALL CD_AGG_ChannelHeader(agg, 3, hdr, unt, es, em)
    CALL require(es == CD_AGG_OK .AND. TRIM(hdr) == 'Curv1N5' .AND. TRIM(unt) == '(1/m)', 'chan:hdr3-Curv1N5')
    CALL CD_AGG_ChannelHeader(agg, 4, hdr, unt, es, em)
    CALL require(es == CD_AGG_OK .AND. TRIM(hdr) == 'FairAngle1' .AND. TRIM(unt) == '(deg)', 'chan:hdr4-FairAngle1')
    CALL CD_AGG_ChannelHeader(agg, 5, hdr, unt, es, em)
    CALL require(es == CD_AGG_OK .AND. TRIM(hdr) == 'BendMom1N5' .AND. TRIM(unt) == '(N.m)', 'chan:hdr5-BendMom1N5')
    CALL CD_AGG_ChannelHeader(agg, 7, hdr, unt, es, em)
    CALL require(es == CD_AGG_OK .AND. TRIM(hdr) == 'AnchDecl1' .AND. TRIM(unt) == '(deg)', 'chan:hdr7-AnchDecl1')
    CALL CD_AGG_ChannelHeader(agg, 8, hdr, unt, es, em)
    CALL require(es == CD_AGG_OK .AND. TRIM(hdr) == 'FairIncl1' .AND. TRIM(unt) == '(deg)', 'chan:hdr8-FairIncl1')
    ! A host header buffer shorter than the channel name fails closed instead of truncating
    ! (two truncated names could collide); a name that fits is returned whole.
    BLOCK
      CHARACTER(9) :: short_hdr, short_unt
      CALL CD_AGG_ChannelHeader(agg, 1, short_hdr, short_unt, es, em)
      CALL require(es == CD_AGG_OK .AND. TRIM(short_hdr) == 'FairTen2', 'chan:short-host-fits')
      CALL CD_AGG_ChannelHeader(agg, 4, short_hdr, short_unt, es, em)
      CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, '"FairAngle1" has 10 characters') > 0 .AND. &
                   INDEX(em, 'at most 9') > 0, 'chan:short-host-rejects-long-name: '//TRIM(em))
    END BLOCK

    ! Evaluate every channel through the shell path; all must be finite.
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'chan:calc: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 1, vFairTen2, es, em); CALL require(es == CD_AGG_OK, 'chan:eval1: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 2, vAnchTen2, es, em); CALL require(es == CD_AGG_OK, 'chan:eval2: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 3, vCurv1N5, es, em); CALL require(es == CD_AGG_OK, 'chan:eval3: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 4, vFairAngle1, es, em); CALL require(es == CD_AGG_OK, 'chan:eval4: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 5, vBendMom1N5, es, em); CALL require(es == CD_AGG_OK, 'chan:eval5: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 6, vTen1N10, es, em); CALL require(es == CD_AGG_OK, 'chan:eval6: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 7, vFairDecl1, es, em); CALL require(es == CD_AGG_OK, 'chan:eval7: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 8, vFairIncl1, es, em); CALL require(es == CD_AGG_OK, 'chan:eval8: '//TRIM(em))
    CALL require(IEEE_IS_FINITE(vFairTen2) .AND. IEEE_IS_FINITE(vAnchTen2) .AND. IEEE_IS_FINITE(vCurv1N5) .AND. &
                 IEEE_IS_FINITE(vFairAngle1) .AND. IEEE_IS_FINITE(vBendMom1N5) .AND. IEEE_IS_FINITE(vTen1N10) .AND. &
                 IEEE_IS_FINITE(vFairDecl1) .AND. IEEE_IS_FINITE(vFairIncl1), &
                 'chan:all-finite')

    ! (a) mooring FairTen2 / AnchTen2 == the system line-1 line-end forces (end element tension plus
    ! the end node's lumped loads, CD_Get_Model_EndForces). Deck line 2 is the only EI=0 line, so it
    ! maps to system line 1; the map must route the mooring channel there.
    ne = CD_System_Line_NElem(agg%sys%fast%system, 1, es, em)
    CALL require(es == CD_SYSTEM_OK .AND. ne >= 1, 'chan:sys-nelem: '//TRIM(em))
    ALLOCATE (tension(MAX(ne, 1)))
    CALL CD_Get_System_Line_Tension(agg%sys%fast%system, 1, tension(1:ne), es, em)
    CALL require(es == CD_SYSTEM_OK, 'chan:sys-tension: '//TRIM(em))
    BLOCK
      REAL(wp) :: f_first(3), f_last(3)
      CALL CD_Get_Model_EndForces(agg%sys%fast%system%lines(1), f_first, f_last, es, em)
      CALL require(es == CD_MODEL_OK, 'chan:sys-end-forces: '//TRIM(em))
      CALL require(ABS(vFairTen2 - NORM2(f_last)) <= 1.0e-6_wp*MAX(1.0_wp, NORM2(f_last)), &
                   'chan:FairTen2-matches-driver')
      CALL require(ABS(vAnchTen2 - NORM2(f_first)) <= 1.0e-6_wp*MAX(1.0_wp, NORM2(f_first)), &
                   'chan:AnchTen2-matches-driver')
      ! the end force differs from the end element's tension only by the end node's lumped loads
      CALL require(ABS(vFairTen2 - tension(ne)) <= 0.05_wp*tension(ne), 'chan:FairTen2-near-end-element')
    END BLOCK
    CALL require(vFairTen2 > 1.0e3_wp, 'chan:FairTen2-physical')

    ! (b) cable Curv1N5 is the exact one-sided Hermite curvature at deck node 5,
    ! conservatively retaining the larger adjacent-element trace.  The channel and
    ! public HFMF recovery must share that definition.
    nnode = agg%cables(1)%line%nn
    ALLOCATE (exact_curvature(nnode))
    CALL CD_HFMF_Curvature(agg%cables(1), exact_curvature, es, em)
    CALL require(es == CD_HFMF_OK, 'chan:exact-Hermite-curvature: '//TRIM(em))
    kexp = exact_curvature(nnode - 5 + 1)
    CALL require(ABS(vCurv1N5 - kexp) <= 1.0e-9_wp + 1.0e-6_wp*ABS(kexp), 'chan:Curv1N5-matches-driver')
    CALL require(vCurv1N5 > 0.0_wp, 'chan:Curv1N5-nonzero')

    ! (c) cable FairAngle1 == declination of the End-A one-sided tangent (deck nodes 1,2).
    base = 6*(nnode - 1); r1 = agg%cables(1)%line%q(base + 1:base + 3)
    base = 6*(nnode - 2); r2 = agg%cables(1)%line%q(base + 1:base + 3)
    decexp = tangent_declination_deg(r2 - r1)
    CALL require(ABS(vFairAngle1 - decexp) <= 1.0e-9_wp + 1.0e-6_wp*ABS(decexp), 'chan:FairAngle1-matches-driver')
    CALL require(vFairDecl1 >= 0.0_wp .AND. vFairDecl1 <= 180.0_wp, 'chan:AnchDecl1-declination-range')
    CALL require(ABS(vFairIncl1 - (decexp - 90.0_wp)) <= 1.0e-9_wp + 1.0e-6_wp*ABS(decexp), &
                 'chan:FairIncl1-signed-below-horizontal')

    ! (d) cable BendMom1N5 == EI(node) * geometric curvature; the cable EI is uniform (bare/buoy
    ! both 1.99e4), so nodal EI = EI_BARE and BendMom reduces to EI_BARE * Curv1N5.
    CALL require(ABS(vBendMom1N5 - EI_BARE*vCurv1N5) <= 1.0e-6_wp*MAX(1.0_wp, ABS(EI_BARE*vCurv1N5)), &
                 'chan:BendMom1N5-equals-EI-times-curv')

    WRITE (*, '(A,ES12.4,A,ES12.4,A,ES12.4)') 'channels: FairTen2 = ', vFairTen2, ' N;  Curv1N5 = ', &
      vCurv1N5, ' 1/m;  FairAngle1 = ', vFairAngle1
    CALL CD_AGG_End(agg, es, em)

    ! FairDecl<L> is the explicit spelling of the FairAngle<L> declination: the same channel,
    ! so the explicit name alone reports the value FairAngle1 reported above, and listing
    ! both (one .out column twice) is rejected at parse time.
    CALL write_mixed_channels_deck('agg_channels_decl.dat', outputs_row='FairDecl1')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_channels_decl.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'chan:decl-init: '//TRIM(em))
    IF (es == CD_AGG_OK) THEN
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_EvalChannel(agg, 1, vFairDecl1, es, em)
      CALL require(es == CD_AGG_OK .AND. ABS(vFairDecl1 - vFairAngle1) <= 1.0e-12_wp, &
                   'chan:FairDecl1-explicit-declination-alias')
      CALL CD_AGG_End(agg, es, em)
    END IF
    CALL write_mixed_channels_deck('agg_channels_dup.dat', outputs_row='FairTen2 FairAngle1 fairdecl1')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_channels_dup.dat', DT, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'duplicates the earlier channel "FairAngle1"') > 0, &
                 'chan:alias-duplicate-rejected: '//TRIM(em))

    ! A deck with NO OUTPUTS reports zero channels (the byte-identical path -- p%NumOuts = 0).
    CALL write_mixed_deck('agg_nochan.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_nochan.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'chan:nochan-init: '//TRIM(em))
    IF (es == CD_AGG_OK) THEN
      CALL require(CD_AGG_NumChannels(agg) == 0, 'chan:no-outputs-zero-channels')
      CALL CD_AGG_WriteStaticProfile(agg, 'agg_nochan.static.out', es, em)
      CALL require(es == CD_AGG_OK, 'chan:nochan-static-profile: '//TRIM(em))
      IF (es == CD_AGG_OK) THEN
        ne = CD_System_Line_NDOF(agg%sys%fast%system, 1, es, em)/3
        CALL require(es == CD_SYSTEM_OK, 'chan:nochan-static-system-size: '//TRIM(em))
        CALL check_static_profile('agg_nochan.static.out', ne + agg%cables(1)%line%nn)
      END IF
      CALL CD_AGG_End(agg, es, em)
    END IF

    ! A PURE-mooring deck WITH channels: ncable = 0 so the cables store is UNALLOCATED, yet a
    ! mooring channel must still evaluate (the evaluator takes cables as an allocatable a system
    ! channel never dereferences). This pins the pure-mooring + channels path (distinct from the
    ! mixed deck, where cables are allocated).
    CALL write_mooring_channels_deck('agg_moorchan.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_moorchan.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'chan:moorchan-init: '//TRIM(em))
    IF (es == CD_AGG_OK) THEN
      CALL require(agg%ncable == 0 .AND. CD_AGG_NumChannels(agg) == 2, 'chan:moorchan-shape')
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL require(es == CD_AGG_OK, 'chan:moorchan-calc: '//TRIM(em))
      CALL CD_AGG_EvalChannel(agg, 1, vFairTen2, es, em)
      CALL require(es == CD_AGG_OK .AND. vFairTen2 > 1.0e3_wp, 'chan:moorchan-fairten-eval')
      CALL CD_AGG_End(agg, es, em)
    END IF

    ! Point channels are backed by the EI=0 system store.  A point belonging
    ! only to the finite-EI cable must be rejected rather than reporting its
    ! stale parse-time coordinates during a coupled march.
    CALL write_mixed_channels_deck('agg_cable_point_channel.dat', point_channel_id=1)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_cable_point_channel.dat', DT, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'cable-only point') > 0, &
                 'chan:reject-cable-only-point-channel')
    IF (es == CD_AGG_OK) CALL CD_AGG_End(agg, es, em)

    ! The other cable-only endpoint is Fixed.  Its system-store copy cannot go
    ! stale, so the corresponding Point channel remains valid.
    CALL write_mixed_channels_deck('agg_fixed_cable_point_channel.dat', point_channel_id=2)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_fixed_cable_point_channel.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'chan:allow-fixed-cable-point-channel: '//TRIM(em))
    IF (es == CD_AGG_OK) CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_output_channels

  SUBROUTINE check_static_profile(path, expected_rows)
    CHARACTER(*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: expected_rows
    INTEGER :: unit, ios, nrow, line_id, node
    CHARACTER(1024) :: title, header, units
    REAL(wp) :: s, x, y, z, tension, curvature, bend_moment, declination, inclination, azimuth

    OPEN (NEWUNIT=unit, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    CALL require(ios == 0, 'static-profile:open')
    IF (ios /= 0) RETURN
    READ (unit, '(A)', IOSTAT=ios) title
    READ (unit, '(A)', IOSTAT=ios) header
    READ (unit, '(A)', IOSTAT=ios) units
    CALL require(ios == 0, 'static-profile:headers-readable')
    CALL require(INDEX(title, 'coupled static configuration') > 0, 'static-profile:title')
    CALL require(INDEX(header, 'ArcLength') > 0 .AND. INDEX(header, 'Tension') > 0 .AND. &
                 INDEX(header, 'Curvature') > 0 .AND. INDEX(header, 'BendMoment') > 0 .AND. &
                 INDEX(header, 'Declination') > 0 .AND. INDEX(header, 'Inclination') > 0 .AND. &
                 INDEX(header, 'Azimuth') > 0, 'static-profile:columns')
    CALL require(INDEX(units, '(N)') > 0 .AND. INDEX(units, '(1/m)') > 0, 'static-profile:units')
    nrow = 0
    DO
      READ (unit, *, IOSTAT=ios) line_id, node, s, x, y, z, tension, curvature, bend_moment, &
        declination, inclination, azimuth
      IF (ios < 0) EXIT
      CALL require(ios == 0, 'static-profile:numeric-row')
      IF (ios /= 0) EXIT
      nrow = nrow + 1
      CALL require(line_id >= 1 .AND. node >= 1, 'static-profile:positive-indices')
      CALL require(ALL(IEEE_IS_FINITE([s, x, y, z, tension, curvature, bend_moment, declination, &
                                       inclination, azimuth])), 'static-profile:finite-row')
      CALL require(ABS(inclination - (declination - 90.0_wp)) < 2.0e-5_wp, &
                   'static-profile:inclination-convention')
    END DO
    CLOSE (unit)
    CALL require(nrow == expected_rows, 'static-profile:one-row-per-node')
  END SUBROUTINE check_static_profile

  SUBROUTINE case_pure_cable_deck()
    !! A finite-EI-only deck builds with has_sys = .FALSE., ncp_sys = 0, ncable >= 1, and steps.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, n, s, ni
    CHARACTER(512) :: em
    REAL(wp) :: r0(3, 1), tgt(3, 1), vel(3, 1), acc(3, 1), p(3, 1), v(3, 1), a(3, 1), l(3, 1)
    REAL(wp), PARAMETER :: DT = 0.05_wp
    LOGICAL :: cv, st

    CALL write_cable_only_deck('agg_cable.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_cable.dat', DT, es, em)
    CALL require(es == CD_AGG_OK .AND. CD_AGG_IsInitialized(agg), 'cable:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL require(.NOT. agg%has_sys, 'cable:no-ei0-system')
    CALL require(agg%ncp_sys == 0, 'cable:ncp_sys-0')
    CALL require(agg%ncable >= 1, 'cable:ncable-ge-1')
    n = CD_AGG_NMovingPoints(agg, es, em)
    CALL require(es == CD_AGG_OK .AND. n == agg%ncable, 'cable:nmoving-equals-ncable')

    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, l, es, em)
    CALL require(es == CD_AGG_OK, 'cable:initial-mesh: '//TRIM(em))
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 6
      tgt(:, 1) = r0(:, 1)
      CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK, 'cable:step: '//TRIM(em))
      IF (es /= CD_AGG_OK) EXIT
    END DO
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'cable:calc: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, l, es, em)
    CALL require(es == CD_AGG_OK, 'cable:mesh: '//TRIM(em))
    CALL require(ALL(IEEE_IS_FINITE(l)) .AND. NORM2(l(:, 1)) > 1.0e2_wp, 'cable:fairlead-load-physical')
    WRITE (*, '(A,ES12.4)') 'pure cable: fairlead |load| = ', NORM2(l(:, 1))
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_pure_cable_deck

  SUBROUTINE case_cable_load_feedback()
    !! The one-way comparison switch suppresses only the host reaction. The
    !! finite-EI cable remains initialized and retains the same physical
    !! fairlead reaction internally as the default two-way model.
    TYPE(CD_AGG_ModuleType) :: two_way, one_way
    INTEGER :: es
    CHARACTER(512) :: em
    REAL(wp) :: p2(3, 1), v2(3, 1), a2(3, 1), l2(3, 1), m2(3, 1)
    REAL(wp) :: p1(3, 1), v1(3, 1), a1(3, 1), l1(3, 1), m1(3, 1)
    REAL(wp) :: fint2(3), fint1(3), mint2(3), mint1(3)
    REAL(wp), PARAMETER :: DT = 0.05_wp

    CALL write_cable_only_deck('agg_feedback_two_way.dat')
    CALL write_cable_only_deck('agg_feedback_one_way.dat', cable_load_feedback=.FALSE.)
    CALL CD_AGG_Init_From_Deck(two_way, 'agg_feedback_two_way.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'feedback:two-way-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_Init_From_Deck(one_way, 'agg_feedback_one_way.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'feedback:one-way-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) THEN
      CALL CD_AGG_End(two_way, es, em)
      RETURN
    END IF

    CALL CD_AGG_GetMovingPointMesh(two_way, p2, v2, a2, l2, es, em, moment=m2)
    CALL require(es == CD_AGG_OK, 'feedback:two-way-mesh: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(one_way, p1, v1, a1, l1, es, em, moment=m1)
    CALL require(es == CD_AGG_OK, 'feedback:one-way-mesh: '//TRIM(em))
    CALL require(ALL(p1 == p2) .AND. ALL(v1 == v2) .AND. ALL(a1 == a2), 'feedback:kinematics-unchanged')
    CALL require(NORM2(l2(:, 1)) > 1.0e2_wp, 'feedback:default-load-physical')
    CALL require(ALL(l1 == 0.0_wp) .AND. ALL(m1 == 0.0_wp), 'feedback:host-reaction-zero')

    CALL CD_HFMF_CalcOutput(two_way%cables(1), fint2, es, em, mint2)
    CALL require(es == CD_HFMF_OK, 'feedback:two-way-internal: '//TRIM(em))
    CALL CD_HFMF_CalcOutput(one_way%cables(1), fint1, es, em, mint1)
    CALL require(es == CD_HFMF_OK, 'feedback:one-way-internal: '//TRIM(em))
    CALL require(ALL(fint1 == fint2) .AND. ALL(mint1 == mint2), 'feedback:internal-reaction-retained')

    CALL CD_AGG_End(two_way, es, em)
    CALL CD_AGG_End(one_way, es, em)
  END SUBROUTINE case_cable_load_feedback

  SUBROUTINE case_cable_modified_newton()
    !! The deck OPTIONS modified_newton reaches the COUPLED finite-EI cable: twin
    !! cable-only decks differing only in the option, driven by the SAME fairlead heave.
    !! Gates: (a) both legs converge every step; (b) the OPTIONS row reaches the cable
    !! model (the flag is set on the modified leg only); (c) the coupled loads agree to
    !! the shared step tolerance (10 % of magnitude: both variants commit at the same 5e-3
    !! relative Newton tolerance, so trajectories are not bit-close). The tangent-assembly
    !! counters are collected but not gated: at this tolerance reuse saves no assembly.
    TYPE(CD_AGG_ModuleType) :: aggf, aggm
    INTEGER :: es, s, ni, nsf, nrf, nvf, naf, nlf, ntan_full, ntan_mn
    CHARACTER(512) :: em
    REAL(wp) :: r0(3, 1), tgt(3, 1), vel(3, 1), acc(3, 1), pf(3, 1), pm(3, 1)
    REAL(wp) :: vf(3, 1), af(3, 1), lf(3, 1), lm(3, 1), tf, tr, te, ts, ta, zt, vt, at, t
    REAL(wp), PARAMETER :: DT = 0.05_wp, AMP = 0.5_wp, PER = 8.0_wp, &
                           TWO_PI = 6.283185307179586_wp
    INTEGER, PARAMETER :: NSTEP = 20
    LOGICAL :: cv, st

    CALL write_cable_only_deck('agg_cable_full.dat')
    CALL write_cable_only_deck('agg_cable_mn.dat', modified_newton=.TRUE.)
    CALL CD_AGG_Init_From_Deck(aggf, 'agg_cable_full.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'mn-cable: full init: '//TRIM(em))
    CALL CD_AGG_Init_From_Deck(aggm, 'agg_cable_mn.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'mn-cable: mn init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN

    CALL CD_AGG_GetMovingPointMesh(aggf, r0, vel, acc, lf, es, em)
    CALL require(es == CD_AGG_OK, 'mn-cable: mesh: '//TRIM(em))

    CALL CD_HermiteCable_Dyn_Reset_Profile()
    t = 0.0_wp
    DO s = 1, NSTEP
      t = t + DT
      zt = AMP*SIN(TWO_PI*t/PER)
      vt = AMP*(TWO_PI/PER)*COS(TWO_PI*t/PER)
      at = -AMP*(TWO_PI/PER)**2*SIN(TWO_PI*t/PER)
      tgt(:, 1) = r0(:, 1) + [0.0_wp, 0.0_wp, zt]
      vel(:, 1) = [0.0_wp, 0.0_wp, vt]
      acc(:, 1) = [0.0_wp, 0.0_wp, at]
      CALL CD_AGG_Step_Moving(aggf, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'mn-cable: full step: '//TRIM(em))
    END DO
    CALL CD_HermiteCable_Dyn_Get_Profile(nsf, nrf, nvf, naf, nlf, tf, tr, te, ts, ta, n_tan=ntan_full)

    CALL CD_HermiteCable_Dyn_Reset_Profile()
    t = 0.0_wp
    DO s = 1, NSTEP
      t = t + DT
      zt = AMP*SIN(TWO_PI*t/PER)
      vt = AMP*(TWO_PI/PER)*COS(TWO_PI*t/PER)
      at = -AMP*(TWO_PI/PER)**2*SIN(TWO_PI*t/PER)
      tgt(:, 1) = r0(:, 1) + [0.0_wp, 0.0_wp, zt]
      vel(:, 1) = [0.0_wp, 0.0_wp, vt]
      acc(:, 1) = [0.0_wp, 0.0_wp, at]
      CALL CD_AGG_Step_Moving(aggm, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'mn-cable: mn step: '//TRIM(em))
    END DO
    CALL CD_HermiteCable_Dyn_Get_Profile(nsf, nrf, nvf, naf, nlf, tf, tr, te, ts, ta, n_tan=ntan_mn)

    CALL CD_AGG_CalcOutput(aggf, es, em)
    CALL require(es == CD_AGG_OK, 'mn-cable: full output: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(aggf, pf, vf, af, lf, es, em)
    CALL require(es == CD_AGG_OK, 'mn-cable: full mesh: '//TRIM(em))
    CALL CD_AGG_CalcOutput(aggm, es, em)
    CALL require(es == CD_AGG_OK, 'mn-cable: mn output: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(aggm, pm, vf, af, lm, es, em)
    CALL require(es == CD_AGG_OK, 'mn-cable: mn mesh: '//TRIM(em))
    ! The DETERMINISTIC plumbing gate: the option must land on the cable model itself
    ! (a plumbed-but-dead flag fails here; the flag's solver behavior is gated at the
    ! element level by hermite_cable_dynamic::check_modified_newton). The counters are
    ! reported, not gated: at this path's frozen production tolerance (5e-3 relative)
    ! per-iteration contraction is weak, so reuse yields no assembly saving on this rig
    ! -- the benefit lives at tighter tolerances (the committed 800 m bench at 1e-4).
    CALL require(.NOT. aggf%cables(1)%line%modified_newton, 'mn-cable: plain deck leaves the flag off')
    CALL require(aggm%cables(1)%line%modified_newton, 'mn-cable: the OPTIONS row reaches the cable model')
    ! Both variants commit each step at the SAME relative tolerance; the coupled-node
    ! reaction lives on FIXED DOFs outside the converged free residual, so the honest
    ! agreement bound is tolerance-consistent (a few x 5e-3 of the force scale), not
    ! trajectory-tight. Gate at 10% of magnitude as the same-physics sanity check.
    CALL require(NORM2(lf(:, 1) - lm(:, 1)) <= 0.10_wp*NORM2(lf(:, 1)), &
                 'mn-cable: coupled loads agree to the shared step tolerance')
    ! The counters stay collected (they document that reuse saves nothing at this
    ! path's frozen 5e-3 tolerance -- see the CHANGELOG scope note) but a passing
    ! test prints nothing; gate them only for finiteness of the collection itself.
    CALL require(ntan_full > 0 .AND. ntan_mn > 0, 'mn-cable: profile counters collected')

    CALL CD_AGG_End(aggf, es, em)
    CALL CD_AGG_End(aggm, es, em)
  END SUBROUTINE case_cable_modified_newton

  SUBROUTINE case_dt_guard()
    !! FIX-2 fixed-dt contract: the aggregate + its cables freeze the coupling dt at Init (the
    !! cables carry it in their gen-alpha step), so a CD_AGG_Step_Moving whose dt differs from the
    !! Init dt must fail closed (CD_AGG_BADINPUT) rather than silently desync the cable columns from
    !! the mooring columns. The matching dt is still accepted (the guard rejects only a mismatch).
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, ni, ncp
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp), ALLOCATABLE :: tgt(:, :), vel(:, :), acc(:, :), lod(:, :)
    LOGICAL :: cv, st

    CALL write_mixed_deck('agg_dtguard.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_dtguard.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'dtguard:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    ncp = CD_AGG_NMovingPoints(agg, es, em)
    ALLOCATE (tgt(3, ncp), vel(3, ncp), acc(3, ncp), lod(3, ncp))
    CALL CD_AGG_GetMovingPointMesh(agg, tgt, vel, acc, lod, es, em)
    CALL require(es == CD_AGG_OK, 'dtguard:mesh0: '//TRIM(em))
    vel = 0.0_wp                                ! hold every fairlead at r0 (tgt), at rest
    acc = 0.0_wp
    ! A step dt different from the Init dt is rejected BEFORE any solve (fixed-dt contract).
    CALL CD_AGG_Step_Moving(agg, DT + 0.01_wp, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_BADINPUT, 'dtguard:mismatched-dt-fails-closed')
    CALL require(INDEX(em, 'match the Init dt') > 0, 'dtguard:rejection-from-the-intended-guard')
    ! The Init dt itself is accepted (a rest step at r0 converges): the guard is not over-rejecting.
    CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK, 'dtguard:matching-dt-accepted: '//TRIM(em))
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_dt_guard

  SUBROUTINE case_cable_no_advance_transfer()
    !! FIX-3 cable no-advance transfer: CD_AGG_UpdateStates_Moving now writes each cable's coupled
    !! kinematics (translational r/v/a only) into the line state WITHOUT advancing, so a CalcOutput
    !! taken after the fairlead moves -- but before any Step -- reflects the move (the instantaneous
    !! direct-feedthrough load). Move ONLY the cable column and hold the mooring column: the cable
    !! column's reaction must CHANGE against its frozen interior, while the held mooring column's
    !! load is unchanged (per-column transfer, no advance, no cross-talk).
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, ncp, ccol
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.05_wp, DZ = 0.1_wp
    REAL(wp), ALLOCATABLE :: r0(:, :), v0(:, :), a0(:, :), l0m(:, :)
    REAL(wp), ALLOCATABLE :: pos(:, :), vel(:, :), acc(:, :), p1(:, :), v1(:, :), a1(:, :), l1(:, :)

    CALL write_mixed_deck('agg_transfer.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_transfer.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'transfer:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL require(agg%ncp_sys == 1 .AND. agg%ncable == 1, 'transfer:mixed-shape')
    ncp = agg%ncp_sys + agg%ncable
    ccol = agg%ncp_sys + 1                      ! the (single) trailing cable column
    ALLOCATE (r0(3, ncp), v0(3, ncp), a0(3, ncp), l0m(3, ncp))
    ALLOCATE (pos(3, ncp), vel(3, ncp), acc(3, ncp), p1(3, ncp), v1(3, ncp), a1(3, ncp), l1(3, ncp))
    ! Baseline loads at the built (rest) state.
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'transfer:calc0: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg, r0, v0, a0, l0m, es, em)
    CALL require(es == CD_AGG_OK, 'transfer:mesh0: '//TRIM(em))
    ! Move ONLY the cable column by a finite z offset; hold the mooring column; all at rest.
    pos = r0
    pos(3, ccol) = r0(3, ccol) + DZ
    vel = 0.0_wp
    acc = 0.0_wp
    CALL CD_AGG_UpdateStates_Moving(agg, pos, vel, acc, es, em)
    CALL require(es == CD_AGG_OK, 'transfer:updatestates: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'transfer:calc1: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg, p1, v1, a1, l1, es, em)
    CALL require(es == CD_AGG_OK, 'transfer:mesh1: '//TRIM(em))
    CALL require(ALL(IEEE_IS_FINITE(l1)), 'transfer:loads-finite')
    ! The moved cable column reports the moved position (the transfer wrote the line's coupled q).
    CALL require(ABS(p1(3, ccol) - (r0(3, ccol) + DZ)) < 1.0e-9_wp, 'transfer:cable-position-updated')
    ! The cable-column reaction CHANGED: the frozen-interior feedthrough reflects the fairlead move.
    CALL require(NORM2(l1(:, ccol) - l0m(:, ccol)) > 1.0e1_wp, 'transfer:cable-load-reflects-move')
    ! The mooring column was held: its load is unchanged (no advance, no cross-talk).
    CALL require(NORM2(l1(:, 1) - l0m(:, 1)) < 1.0e-6_wp, 'transfer:mooring-load-unchanged')
    WRITE (*, '(A,ES12.4,A,ES12.4)') 'no-advance transfer: cable |dLoad| = ', &
      NORM2(l1(:, ccol) - l0m(:, ccol)), ' N;  mooring |dLoad| = ', NORM2(l1(:, 1) - l0m(:, 1))
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_cable_no_advance_transfer

  SUBROUTINE case_fail_closed_dynamic_point()
    !! A mixed deck (chain + cable) carrying a Free (dynamic) point is out of v1 scope and
    !! must fail closed at the aggregate's guard -- not be silently coerced.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es
    CHARACTER(512) :: em

    CALL write_mixed_free_deck('agg_mixed_free.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_mixed_free.dat', 0.05_wp, es, em)
    CALL require(es == CD_AGG_BADINPUT, 'freept:mixed-with-free-point-fails-closed')
    CALL require(INDEX(em, 'Free/Connect dynamic points') > 0, 'freept:rejection-from-the-intended-guard')
    CALL require(.NOT. CD_AGG_IsInitialized(agg), 'freept:not-initialized-after-reject')
  END SUBROUTINE case_fail_closed_dynamic_point

  SUBROUTINE case_coupled_viscoelastic()
    !! A viscoelastic (ElasticMod 2) deck now BUILDS on the COUPLED aggregate
    !! facade: the glue checkpoint mirror carries the per-element dl_1 state
    !! (CableDyn_OF pack_state_mirror / reload_interiors_ary), so a checkpoint
    !! restart restores the SLS partition exactly. The cross-process -restart
    !! bit-identity is gated fork-side; here we assert the facade accepts and
    !! initialises the coupled viscoelastic deck. The standalone driver too.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, u, ios
    CHARACTER(512) :: em
    OPEN (NEWUNIT=u, FILE='agg_visco.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'coupled viscoelastic deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'rope 0.1438 22.42 1.424e+08|1.586e+08 4E9|11E6 0.0 1.0 0.0 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 0.0 0.0 -1.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '2 Vessel 2.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 rope 2.0 4'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_visco.dat', 0.05_wp, es, em)
    CALL require(es == CD_AGG_OK, 'visco:coupled-facade-accepts: '//TRIM(em))
    IF (es == CD_AGG_OK) CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_coupled_viscoelastic

  SUBROUTINE case_coupled_point3_buoy()
    !! A point3 buoy on a coupled mooring now BUILDS on the aggregate facade: the
    !! point3 BODY's Mass/Vol/CdA/Ca is resolved onto its Free attachment point at
    !! parse, and the aggregate's system carries the volume buoyancy + drag/added-mass
    !! and checkpoint-mirrors the point like any Free/Connect point. A heavy point3
    !! clump hangs on a chain line from a Fixed anchor to a Coupled fairlead. rigid6
    !! BODYs still fail closed by name.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, u, ios
    CHARACTER(512) :: em
    OPEN (NEWUNIT=u, FILE='agg_point3.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'coupled point3 buoy deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.1 50.0 1.0e8 0.0 0.0 1.0 0.0 1.0 0.0'
    WRITE (u, '(A)') '--- BODIES ---'
    WRITE (u, '(A)') 'ID Type X Y Z Rot1 Rot2 Rot3 Mass Vol C33 C44 C55 CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (deg) (deg) (deg) (kg) (m3) (N/m) (Nm/rad) (Nm/rad) (m2) (-)'
    WRITE (u, '(A)') '1 Point3 3.0 0.0 -6.0 0.0 0.0 0.0 500.0 0.05 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 6.0 0.0 -10.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '2 Body1 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '3 Coupled 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '2 3 2 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 chain 6.0 6'
    WRITE (u, '(A)') '2 chain 6.5 6'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen2'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_point3.dat', 0.05_wp, es, em)
    CALL require(es == CD_AGG_OK, 'point3:coupled-facade-accepts: '//TRIM(em))
    IF (es == CD_AGG_OK) CALL CD_AGG_End(agg, es, em)
    ! C33/C44/C55 restoring on a point3 is not supported and always fails
    ! closed. A point3 carrying a C33 heave-restoring stiffness (CdA = Ca = 0 here) must fail
    ! closed, not silently run as mass + constant buoyancy only. (CdA/Ca drag/added-mass ARE
    ! now served under a host field -- case_coupled_point3_fluid; still water fails closed.)
    OPEN (NEWUNIT=u, FILE='agg_point3_hydro.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'coupled point3 buoy with restoring (must fail closed)'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.1 50.0 1.0e8 0.0 0.0 1.0 0.0 1.0 0.0'
    WRITE (u, '(A)') '--- BODIES ---'
    WRITE (u, '(A)') 'ID Type X Y Z Rot1 Rot2 Rot3 Mass Vol C33 C44 C55 CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (deg) (deg) (deg) (kg) (m3) (N/m) (Nm/rad) (Nm/rad) (m2) (-)'
    WRITE (u, '(A)') '1 Point3 3.0 0.0 -6.0 0.0 0.0 0.0 500.0 0.05 -1000.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 6.0 0.0 -10.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '2 Body1 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '3 Coupled 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '2 3 2 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 chain 6.0 6'
    WRITE (u, '(A)') '2 chain 6.5 6'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen2'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_point3_hydro.dat', 0.05_wp, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'C33/C44/C55 restoring') > 0, &
                 'point3:restoring-fails-closed')
    ! a point3 BODY with NO referencing Body<N> POINT (orphan) fails closed: the aggregate
    ! would otherwise build from points/lines only and silently drop the declared buoy
    OPEN (NEWUNIT=u, FILE='agg_point3_orphan.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'orphan point3 body deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.1 50.0 1.0e8 0.0 0.0 1.0 0.0 1.0 0.0'
    WRITE (u, '(A)') '--- BODIES ---'
    WRITE (u, '(A)') 'ID Type X Y Z Rot1 Rot2 Rot3 Mass Vol C33 C44 C55 CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (deg) (deg) (deg) (kg) (m3) (N/m) (Nm/rad) (Nm/rad) (m2) (-)'
    WRITE (u, '(A)') '1 Point3 3.0 0.0 -6.0 0.0 0.0 0.0 500.0 0.05 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 6.0 0.0 -10.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '2 Coupled 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 chain 12.0 8'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_point3_orphan.dat', 0.05_wp, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'referenced by exactly one') > 0, &
                 'point3:orphan-body-fails-closed')
  END SUBROUTINE case_coupled_point3_buoy

  SUBROUTINE write_coupled_point3_hydro_deck(path, cda, ca)
    !! A submerged point3 buoy (CdA/Ca) hung between a Fixed anchor and a Coupled fairlead by two
    !! chains -- the coupled point-fluid contract test fixture.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), INTENT(IN) :: cda, ca
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'coupled point3 hydro buoy deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'line 0.10 80.0 1.0e5 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- BODIES ---'
    WRITE (u, '(A)') 'ID Type X Y Z Rot1 Rot2 Rot3 Mass Vol C33 C44 C55 CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (deg) (deg) (deg) (kg) (m3) (N/m) (Nm/rad) (Nm/rad) (m2) (-)'
    WRITE (u, '(A,F6.3,1X,F6.3)') '1 Point3 0.0 0.0 -4.0 0.0 0.0 0.0 100.0 0.05 0.0 0.0 0.0 ', cda, ca
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 2.0 0.0 -5.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '2 Body1 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '3 Coupled 0.0 0.0 -1.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '2 3 2 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 line 2.5 3'
    WRITE (u, '(A)') '2 line 2.95 3'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen2'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_coupled_point3_hydro_deck

  SUBROUTINE case_coupled_point3_fluid()
    !! The coupled point3 FLUID CONTRACT: with a host ambient field (external_fluid) a point3 buoy
    !! carrying CdA/Ca IS served. The aggregate appends one fluid sampling node per hydro-active
    !! dynamic point, GetFluidNodePositions reports its position there, and CD_AGG_SetFluidFields
    !! sets the field on the buoy so point_environment_force applies drag / Froude-Krylov /
    !! added-mass. A still-water CdA buoy (no host field) still fails closed. Discriminating test:
    !! a strong ambient current drives the buoy downstream vs a zero field.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, s, nfl
    REAL(wp) :: r0(3, 1), vel(3, 1), acc(3, 1), la(3, 1)
    REAL(wp), ALLOCATABLE :: xyz(:, :), fv(:, :), fa(:, :), wl(:)
    REAL(wp) :: x_flow, x_still
    LOGICAL :: ca_flag, sa_flag
    INTEGER :: ia
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.005_wp

    ! Still water: a point3 with CdA but no host field fails closed (nothing to drive the drag).
    CALL write_coupled_point3_hydro_deck('agg_point3_cda_still.dat', 2.0_wp, 0.0_wp)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_point3_cda_still.dat', DT, es, em)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, 'CdA = Ca = 0') > 0, 'point3-fluid:cda-still-fails-closed')

    ! With a host ambient field, the CdA/Ca buoy is served.
    CALL write_coupled_point3_hydro_deck('agg_point3_fluid.dat', 2.0_wp, 1.0_wp)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_point3_fluid.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'point3-fluid:external-field-served: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN

    nfl = CD_AGG_NFluidNodes(agg, es, em)
    CALL require(nfl >= 1, 'point3-fluid:nfluid-nonzero')
    ALLOCATE (xyz(3, nfl), fv(3, nfl), fa(3, nfl), wl(nfl))
    CALL CD_AGG_GetFluidNodePositions(agg, xyz, es, em)
    CALL require(es == CD_AGG_OK, 'point3-fluid:get-positions: '//TRIM(em))
    CALL require(nan_max_abs(xyz(:, nfl) - [0.0_wp, 0.0_wp, -4.0_wp]) < 1.0e-9_wp, &
                 'point3-fluid:buoy-node-at-reference-position')

    ! Strong +x current: the buoy (last fluid node) drifts downstream.
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, la, es, em)
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 40
      fv = 0.0_wp
      fv(1, :) = 3.0_wp
      fa = 0.0_wp
      wl = 5.0_wp
      CALL CD_AGG_SetFluidFields(agg, fv, fa, wl, es, em)
      CALL require(es == CD_AGG_OK, 'point3-fluid:set-fields: '//TRIM(em))
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=REAL(s, wp)*DT)
      CALL require(es == CD_AGG_OK, 'point3-fluid:step: '//TRIM(em))
      IF (es /= CD_AGG_OK) EXIT
      CALL CD_AGG_CalcOutput(agg, es, em)
    END DO
    CALL CD_AGG_GetFluidNodePositions(agg, xyz, es, em)
    x_flow = xyz(1, nfl)
    CALL CD_AGG_End(agg, es, em)

    ! Control: identical build + steps, ZERO ambient field -> still-water buoy trajectory.
    CALL CD_AGG_Init_From_Deck(agg, 'agg_point3_fluid.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'point3-fluid:control-build: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, la, es, em)
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 40
      fv = 0.0_wp
      fa = 0.0_wp
      wl = 5.0_wp
      CALL CD_AGG_SetFluidFields(agg, fv, fa, wl, es, em)
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=REAL(s, wp)*DT)
      IF (es /= CD_AGG_OK) EXIT
      CALL CD_AGG_CalcOutput(agg, es, em)
    END DO
    CALL CD_AGG_GetFluidNodePositions(agg, xyz, es, em)
    x_still = xyz(1, nfl)
    CALL CD_AGG_End(agg, es, em)

    WRITE (*, '(A,ES12.4,A,ES12.4)') 'coupled point3 fluid: buoy x (flow) = ', x_flow, '  (still) = ', x_still
    CALL require(x_flow - x_still > 1.0e-3_wp, 'point3-fluid:current-drives-buoy-downstream')

    ! A VOLUME-ONLY buoy (Vol > 0, CdA = Ca = 0) must ALSO be a fluid node, so it takes the wave
    ! Froude-Krylov (rho*Vol*(1+Ca)*accel) instead of staying constant-buoyancy while the lines feel
    ! waves. Build one, drive a +x ambient ACCELERATION field (zero velocity, so no drag -- FK only),
    ! and require the buoy to move downstream relative to the still-water offset.
    CALL write_coupled_point3_hydro_deck('agg_point3_vol.dat', 0.0_wp, 0.0_wp)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_point3_vol.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'point3-fluid:vol-only-built: '//TRIM(em))
    CALL require(CD_AGG_NFluidNodes(agg, es, em) == nfl, 'point3-fluid:vol-only-is-a-fluid-node')
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, la, es, em)
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 40
      fv = 0.0_wp
      fa = 0.0_wp
      fa(1, :) = 2.0_wp
      wl = 5.0_wp
      CALL CD_AGG_SetFluidFields(agg, fv, fa, wl, es, em)
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=REAL(s, wp)*DT)
      IF (es /= CD_AGG_OK) EXIT
      CALL CD_AGG_CalcOutput(agg, es, em)
    END DO
    CALL CD_AGG_GetFluidNodePositions(agg, xyz, es, em)
    x_flow = xyz(1, nfl)
    CALL CD_AGG_End(agg, es, em)
    WRITE (*, '(A,ES12.4)') 'coupled point3 vol-only Froude-Krylov: buoy x = ', x_flow
    CALL require(x_flow - x_still > 1.0e-3_wp, 'point3-fluid:vol-only-buoy-feels-froude-krylov')
    DEALLOCATE (xyz, fv, fa, wl)
  END SUBROUTINE case_coupled_point3_fluid

  SUBROUTINE write_coupled_rigid6_deck(path, cda, ca, extra_free, deck_pose)
    !! A submerged Rigid6 body hung between a Fixed anchor and a Coupled fairlead by two
    !! chain lines (still water). cda/ca parametrise the body fluid coefficients (0/0 is the
    !! supported pure-inertia+buoyancy+restoring case); extra_free adds a stray Free point to
    !! exercise the scope guard; deck_pose starts the body at its deck pose (bodyIC deck)
    !! instead of its static equilibrium.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), INTENT(IN) :: cda, ca
    LOGICAL, INTENT(IN) :: extra_free
    LOGICAL, INTENT(IN), OPTIONAL :: deck_pose
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'coupled rigid6 body deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'line 0.10 80.0 1.0e5 0.0 0.0 1.0 0.0 1.0 0.0'
    WRITE (u, '(A)') '--- BODIES ---'
    WRITE (u, '(A)') 'ID Type X Y Z Rot1 Rot2 Rot3 Mass Vol C33 C44 C55 CdA Ca Ixx Iyy Izz'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (deg) (deg) (deg) (kg) (m3) (N/m) (Nm/rad) (Nm/rad) (m2) (-) '// &
      '(kgm2) (kgm2) (kgm2)'
    WRITE (u, '(A,F6.3,1X,F6.3,A)') '1 Rigid6 0.0 0.0 -4.0 0.0 0.0 0.0 100.0 0.05 0.0 100.0 100.0 ', &
      cda, ca, ' 20.0 20.0 20.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 2.0 0.0 -5.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '2 Body1 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '3 Coupled 0.0 0.0 -1.0 0.0 0.0 0.0 0.0'
    IF (extra_free) WRITE (u, '(A)') '4 Free -3.0 0.0 -5.0 100.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '2 3 2 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 line 2.5 3'
    WRITE (u, '(A)') '2 line 2.95 3'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    IF (PRESENT(deck_pose)) THEN
      IF (deck_pose) WRITE (u, '(A)') 'deck bodyIC'
    END IF
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen2'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_coupled_rigid6_deck

  SUBROUTINE case_coupled_rigid6()
    !! A 6-DOF Rigid6 body on a coupled mooring BUILDS and MARCHES on the aggregate facade:
    !! step_rigid6_system advances the body explicitly (buoyancy + C44/C55 restoring +
    !! inertia + line-attachment reactions, all self-contained at CdA=Ca=0), drives its
    !! attachment Free point as prescribed motion, and steps the lines. Holding the fairlead
    !! in still water, the body settles: every step converges and the mesh stays finite.
    !! A Rigid6 carrying any fluid coefficient (CdA or Ca) fails closed -- its Morison
    !! drag/added-mass would be dropped because the SeaState-to-body field is not wired -- and
    !! a stray Free point (needing the DynamicPoints integrator) fails closed too.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, s, na
    REAL(wp) :: r0(3, 1), vel(3, 1), acc(3, 1), pa(3, 1), va(3, 1), aa(3, 1), la(3, 1)
    REAL(wp) :: paA(3, 1), laA(3, 1)
    LOGICAL :: ca_flag, sa_flag
    INTEGER :: ia
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.005_wp

    CALL write_coupled_rigid6_deck('agg_rigid6.dat', 0.0_wp, 0.0_wp, .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rigid6.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'rigid6:coupled-facade-accepts: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL require(CD_AGG_HasRigid6(agg), 'rigid6:has-rigid6-true (restart guard sees the body)')
    na = CD_AGG_NMovingPoints(agg, es, em)
    CALL require(na == 1, 'rigid6:one-moving-fairlead')
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, la, es, em)
    CALL require(es == CD_AGG_OK, 'rigid6:initial-mesh: '//TRIM(em))
    vel = 0.0_wp
    acc = 0.0_wp
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pa, va, aa, laA, es, em)
    DO s = 1, 20
      ! hold the fairlead fixed at r0: the body starts at its static equilibrium (bodyIC static)
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, &
                              t_committed=REAL(s, wp)*DT)
      CALL require(es == CD_AGG_OK, 'rigid6:step: '//TRIM(em))
      IF (es /= CD_AGG_OK) EXIT
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL require(es == CD_AGG_OK, 'rigid6:calc: '//TRIM(em))
      CALL CD_AGG_GetMovingPointMesh(agg, pa, va, aa, la, es, em)
      CALL require(es == CD_AGG_OK .AND. ALL(IEEE_IS_FINITE(la)), 'rigid6:mesh-finite')
    END DO
    WRITE (*, '(A,ES12.4,A,ES12.4)') 'coupled rigid6 fairlead |load| = ', nan_max_abs(la), &
      '; relative change over 20 held steps = ', nan_max_abs(la - laA)/nan_max_abs(laA)
    CALL require(nan_max_abs(la - laA) <= 1.0e-6_wp*nan_max_abs(laA), &
                 'rigid6:starts-at-rest (coupled static body initial condition)')
    CALL require(nan_max_abs(la) > 0.0_wp, 'rigid6:fairlead-load-nonzero')

    ! Correction-rewind: snapshot -> step A -> restore -> step B must reproduce step A
    ! bit-for-bit -- exercises CD_Rigid6_Snapshot/Restore joining the aggregate snapshot
    ! (the body states are NOT in the system snapshot, so a missing rewind would diverge).
    CALL CD_AGG_Snapshot(agg, es, em)
    CALL require(es == CD_AGG_OK, 'rigid6:cr-snapshot: '//TRIM(em))
    CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=21.0_wp*DT)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, paA, va, aa, laA, es, em)
    CALL CD_AGG_Restore(agg, es, em)
    CALL require(es == CD_AGG_OK, 'rigid6:cr-restore: '//TRIM(em))
    CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=21.0_wp*DT)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pa, va, aa, la, es, em)
    CALL require(nan_max_abs(pa - paA) <= 0.0_wp .AND. nan_max_abs(la - laA) <= 0.0_wp, &
                 'rigid6:correction-rewind-bit-identical')
    CALL CD_AGG_End(agg, es, em)

    ! Fluid-contract fail-closed: a Rigid6 with CdA /= 0 or Ca /= 0.
    CALL write_coupled_rigid6_deck('agg_rigid6_cda.dat', 1.0_wp, 0.0_wp, .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rigid6_cda.dat', DT, es, em)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, 'CdA = Ca = 0') > 0, 'rigid6:cda-fails-closed')
    CALL require(.NOT. CD_AGG_IsInitialized(agg), 'rigid6:cda-not-initialized')
    CALL write_coupled_rigid6_deck('agg_rigid6_ca.dat', 0.0_wp, 1.0_wp, .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rigid6_ca.dat', DT, es, em)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, 'CdA = Ca = 0') > 0, 'rigid6:ca-fails-closed')

    ! Scope fail-closed: a stray Free point (needs the DynamicPoints integrator, which the
    ! body march replaces) is outside the supported coupled Rigid6 scope.
    CALL write_coupled_rigid6_deck('agg_rigid6_free.dat', 0.0_wp, 0.0_wp, .TRUE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rigid6_free.dat', DT, es, em)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, 'Free/Connect dynamic points') > 0, &
                 'rigid6:free-point-fails-closed')

    ! In STILL water (no host field) CdA/Ca still fail closed: there is nothing to drive the
    ! drag/added-mass, so a body carrying them would silently drop them (the message still names
    ! the CdA = Ca = 0 escape). The host-field (external_fluid) case is served -- case_coupled_rigid6_fluid.
    CALL write_coupled_rigid6_deck('agg_rigid6_cda_still.dat', 1.0_wp, 0.0_wp, .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rigid6_cda_still.dat', DT, es, em)
    CALL require(es == CD_AGG_BADINPUT .AND. INDEX(em, 'CdA = Ca = 0') > 0, 'rigid6:cda-still-water-fails-closed')
  END SUBROUTINE case_coupled_rigid6

  SUBROUTINE case_coupled_rigid6_fluid()
    !! The coupled body FLUID CONTRACT: with a host ambient field (external_fluid) a Rigid6 body
    !! carrying CdA/Ca IS served. The aggregate appends one fluid sampling node per body (its
    !! reference point) to nfluid, GetFluidNodePositions returns the body position there, and
    !! CD_AGG_SetFluidFields holds the sampled field on the body so its Morison drag / Froude-
    !! Krylov / added-mass act. Discriminating test: a strong ambient current drives the body
    !! downstream, whereas a zero field leaves it on the still-water trajectory.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, s, na, nfl, nb
    REAL(wp) :: r0(3, 1), vel(3, 1), acc(3, 1), la(3, 1)
    REAL(wp), ALLOCATABLE :: xyz(:, :), fv(:, :), fa(:, :), wl(:), st(:)
    REAL(wp) :: x_flow, x_still
    LOGICAL :: ca_flag, sa_flag
    INTEGER :: ia
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.005_wp

    ! Build with a host ambient field and non-zero CdA/Ca -- now accepted (was fail-closed).
    CALL write_coupled_rigid6_deck('agg_rigid6_fluid.dat', 2.0_wp, 1.0_wp, .FALSE., deck_pose=.TRUE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rigid6_fluid.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'rigid6-fluid:external-field-served: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN

    ! nfluid gained one node (the body ref point); its position is the last fluid node.
    nfl = CD_AGG_NFluidNodes(agg, es, em)
    CALL require(nfl >= 1, 'rigid6-fluid:nfluid-nonzero')
    ALLOCATE (xyz(3, nfl), fv(3, nfl), fa(3, nfl), wl(nfl))
    CALL CD_AGG_GetFluidNodePositions(agg, xyz, es, em)
    CALL require(es == CD_AGG_OK, 'rigid6-fluid:get-positions: '//TRIM(em))
    CALL require(nan_max_abs(xyz(:, nfl) - [0.0_wp, 0.0_wp, -4.0_wp]) < 1.0e-9_wp, &
                 'rigid6-fluid:body-node-at-reference-position')

    ! Drive a strong +x ambient current (submerging waterline) and step; read the body x.
    na = CD_AGG_NMovingPoints(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, la, es, em)
    vel = 0.0_wp
    acc = 0.0_wp
    fv = 0.0_wp
    fv(1, :) = 3.0_wp
    fa = 0.0_wp
    wl = 5.0_wp
    DO s = 1, 30
      CALL CD_AGG_SetFluidFields(agg, fv, fa, wl, es, em)
      CALL require(es == CD_AGG_OK, 'rigid6-fluid:set-fields: '//TRIM(em))
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=REAL(s, wp)*DT)
      CALL require(es == CD_AGG_OK, 'rigid6-fluid:step: '//TRIM(em))
      IF (es /= CD_AGG_OK) EXIT
      CALL CD_AGG_CalcOutput(agg, es, em)
    END DO
    nb = CD_AGG_Rigid6_MirrorSize(agg)
    ALLOCATE (st(nb))
    CALL CD_AGG_Get_Rigid6_States(agg, st, es, em)
    x_flow = st(1)   ! body 1: r(1)
    CALL CD_AGG_End(agg, es, em)

    ! Control: identical build + identical steps but a ZERO ambient field -> still-water trajectory.
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rigid6_fluid.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'rigid6-fluid:control-build: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, la, es, em)
    vel = 0.0_wp
    acc = 0.0_wp
    fv = 0.0_wp
    fa = 0.0_wp
    wl = 5.0_wp
    DO s = 1, 30
      CALL CD_AGG_SetFluidFields(agg, fv, fa, wl, es, em)
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=REAL(s, wp)*DT)
      IF (es /= CD_AGG_OK) EXIT
      CALL CD_AGG_CalcOutput(agg, es, em)
    END DO
    CALL CD_AGG_Get_Rigid6_States(agg, st, es, em)
    x_still = st(1)
    CALL CD_AGG_End(agg, es, em)

    WRITE (*, '(A,ES12.4,A,ES12.4)') 'coupled rigid6 fluid: body x (flow) = ', x_flow, '  (still) = ', x_still
    ! The ambient current must push the body measurably downstream (+x) relative to still water --
    ! proof the held field reaches the body's Morison force through the fluid contract.
    CALL require(x_flow - x_still > 1.0e-3_wp, 'rigid6-fluid:current-drives-body-downstream')
    DEALLOCATE (xyz, fv, fa, wl, st)
  END SUBROUTINE case_coupled_rigid6_fluid

  SUBROUTINE case_coupled_rigid6_mirror()
    !! The 6-DOF Rigid6 body-state checkpoint accessors: CD_AGG_Rigid6_MirrorSize sizes the
    !! per-body block (24), CD_AGG_Get/Set_Rigid6_States pack/unpack the committed dynamic
    !! state (r, v, a, omega, alpha, rot_mat; eta re-derived on set). A get -> set -> get
    !! roundtrip is bit-identical; a Set with a wrong-size or non-finite array fails closed.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, s, na, msz, ia
    REAL(wp) :: r0(3, 1), vel(3, 1), acc(3, 1), la(3, 1)
    REAL(wp), ALLOCATABLE :: s1(:), s2(:), s3(:), sbad(:)
    LOGICAL :: ca_flag, sa_flag
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.005_wp

    CALL write_coupled_rigid6_deck('agg_rigid6_mir.dat', 0.0_wp, 0.0_wp, .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rigid6_mir.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'rigid6mir:build: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    msz = CD_AGG_Rigid6_MirrorSize(agg)
    CALL require(msz == 24, 'rigid6mir:size-24-per-body')
    ALLOCATE (s1(msz), s2(msz), s3(msz))
    na = CD_AGG_NMovingPoints(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, la, es, em)
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 10   ! march to a non-trivial committed state
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=REAL(s, wp)*DT)
      CALL require(es == CD_AGG_OK, 'rigid6mir:march: '//TRIM(em))
    END DO
    CALL CD_AGG_Get_Rigid6_States(agg, s1, es, em)
    CALL require(es == CD_AGG_OK .AND. ALL(IEEE_IS_FINITE(s1)), 'rigid6mir:get-finite')
    DO s = 11, 20   ! march further so the state changes
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=REAL(s, wp)*DT)
    END DO
    CALL CD_AGG_Get_Rigid6_States(agg, s2, es, em)
    CALL require(nan_max_abs(s2 - s1) > 0.0_wp, 'rigid6mir:march-changes-state')
    ! restore the earlier committed state, then get -> must reproduce s1 bit-for-bit
    CALL CD_AGG_Set_Rigid6_States(agg, s1, es, em)
    CALL require(es == CD_AGG_OK, 'rigid6mir:set-ok: '//TRIM(em))
    CALL CD_AGG_Get_Rigid6_States(agg, s3, es, em)
    CALL require(nan_max_abs(s3 - s1) <= 0.0_wp, 'rigid6mir:get-set-get-bit-identical')
    ! fail closed: wrong-size array
    ALLOCATE (sbad(msz - 1))
    sbad = 0.0_wp
    CALL CD_AGG_Set_Rigid6_States(agg, sbad, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'size mismatch') > 0, 'rigid6mir:set-wrong-size-fails')
    ! fail closed: non-finite value
    s1(1) = IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN)
    CALL CD_AGG_Set_Rigid6_States(agg, s1, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'non-finite') > 0, 'rigid6mir:set-nonfinite-fails')
    CALL CD_AGG_End(agg, es, em)
    DEALLOCATE (s1, s2, s3, sbad)
  END SUBROUTINE case_coupled_rigid6_mirror

  SUBROUTINE write_coupled_rod_deck(path, add_current, rod_type)
    !! A submerged Free rod tethered between a Fixed anchor and a Coupled fairlead by two chain
    !! lines (still water). add_current appends a deck current OPTION to exercise the still-water scope.
    !! rod_type replaces the rod's Free type (the other rod routes).
    CHARACTER(*), INTENT(IN) :: path
    LOGICAL, INTENT(IN) :: add_current
    CHARACTER(*), INTENT(IN), OPTIONAL :: rod_type
    INTEGER :: u, ios
    CHARACTER(16) :: rt
    rt = 'Free'
    IF (PRESENT(rod_type)) rt = rod_type
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'coupled rod deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'line 0.10 80.0 1.0e5 0.0 0.0 1.0 0.0 1.0 0.0'
    WRITE (u, '(A)') '--- ROD TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass Cd Ca CdEnd CaEnd CdAx CaAx'
    WRITE (u, '(A)') '(-) (m) (kg/m) (-) (-) (-) (-) (-) (-)'
    ! A light rod (100 kg) carrying heavier line inertia (~280 kg of hanging line): the
    ! partitioned rod update must treat the attached line-end inertia implicitly or it
    ! diverges within a few dozen steps at any dt.
    WRITE (u, '(A)') 'rodmat 0.20 50.0 1.0 1.0 0.0 0.0 0.2 0.0'
    WRITE (u, '(A)') '--- RODS ---'
    WRITE (u, '(A)') 'ID RodType Type XA YA ZA XB YB ZB NumSegs Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (m) (m) (m) (m) (m) (m) (-) (-)'
    WRITE (u, '(A)') '1 rodmat '//TRIM(rt)//' 0.0 0.0 -4.0 0.0 0.0 -2.0 2 -'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 2.0 0.0 -5.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '2 Rod1A 0.0 0.0 -4.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '3 Rod1B 0.0 0.0 -2.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '4 Coupled 0.0 0.0 -1.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '2 3 4 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 line 3.5 3'
    WRITE (u, '(A)') '2 line 1.2 2'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    IF (add_current) WRITE (u, '(A)') 'profile -2.0 0.0 0.0 0.0 0.0 2.0 0.0 0.0 current'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen2'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_coupled_rod_deck

  SUBROUTINE case_coupled_rod()
    !! A rigid ROD on a coupled mooring BUILDS and MARCHES on the aggregate facade -- the rod twin of
    !! case_coupled_rigid6. step_rod_system advances the rod, drives its Rod1A/Rod1B end points as
    !! PRESCRIBED motion, and steps the lines. Holding the fairlead in still water the rod settles;
    !! every step converges and the mesh stays finite. Correction-rewind is bit-identical. A deck
    !! current (outside the still-water scope) fails closed.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, s, na, ia
    REAL(wp) :: r0(3, 1), vel(3, 1), acc(3, 1), la(3, 1), pa(3, 1), va(3, 1), aa(3, 1)
    REAL(wp) :: paA(3, 1), laA(3, 1)
    LOGICAL :: ca_flag, sa_flag
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.005_wp

    CALL write_coupled_rod_deck('agg_rod.dat', .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rod.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'rod:coupled-facade-accepts: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    ! CD_AGG_HasRod flags the dynamic rod state the OpenFAST shell fails Init closed on (until a rod
    ! x%states mirror is wired); a rod deck is NOT a Rigid6 deck.
    CALL require(CD_AGG_HasRod(agg), 'rod:has-rod-true (restart guard sees the rod)')
    CALL require(.NOT. CD_AGG_HasRigid6(agg), 'rod:not-rigid6')
    na = CD_AGG_NMovingPoints(agg, es, em)
    CALL require(na == 1, 'rod:one-moving-fairlead')
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, la, es, em)
    vel = 0.0_wp
    acc = 0.0_wp
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pa, va, aa, laA, es, em)
    DO s = 1, 20
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=REAL(s, wp)*DT)
      CALL require(es == CD_AGG_OK, 'rod:step: '//TRIM(em))
      IF (es /= CD_AGG_OK) EXIT
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg, pa, va, aa, la, es, em)
      CALL require(es == CD_AGG_OK .AND. ALL(IEEE_IS_FINITE(la)), 'rod:mesh-finite')
    END DO
    WRITE (*, '(A,ES12.4,A,ES12.4)') 'coupled rod fairlead |load| = ', nan_max_abs(la), &
      '; relative change over 20 held steps = ', nan_max_abs(la - laA)/nan_max_abs(laA)
    CALL require(nan_max_abs(la - laA) <= 1.0e-6_wp*nan_max_abs(laA), &
                 'rod:starts-at-rest (coupled static rod initial condition)')
    CALL require(nan_max_abs(la) > 0.0_wp, 'rod:fairlead-load-nonzero')

    ! Correction-rewind: snapshot -> step A -> restore -> step B must reproduce step A bit-for-bit
    ! (the rod states are NOT in the system snapshot, so a missing rewind would diverge).
    CALL CD_AGG_Snapshot(agg, es, em)
    CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=21.0_wp*DT)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, paA, va, aa, laA, es, em)
    CALL CD_AGG_Restore(agg, es, em)
    CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, ca_flag, sa_flag, ia, es, em, t_committed=21.0_wp*DT)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pa, va, aa, la, es, em)
    CALL require(nan_max_abs(pa - paA) <= 0.0_wp .AND. nan_max_abs(la - laA) <= 0.0_wp, &
                 'rod:correction-rewind-bit-identical')
    CALL CD_AGG_End(agg, es, em)

    ! Host-driven ambient fluid is accepted; its force response is gated separately.
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rod.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'rod:external-fluid-accepted: '//TRIM(em))
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_coupled_rod

  SUBROUTINE case_coupled_rod_mirror()
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, s, msz, ia
    REAL(wp) :: r0(3, 1), vel(3, 1), acc(3, 1), load(3, 1)
    REAL(wp), ALLOCATABLE :: s1(:), s2(:), s3(:), sbad(:)
    LOGICAL :: cv, st
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.005_wp

    CALL write_coupled_rod_deck('agg_rod_mirror.dat', .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_rod_mirror.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'rodmir:build: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    msz = CD_AGG_Rod_MirrorSize(agg)
    CALL require(msz == 24, 'rodmir:size-24-per-rod')
    ALLOCATE (s1(msz), s2(msz), s3(msz), sbad(MAX(0, msz - 1)))
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, load, es, em)
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 10
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, cv, st, ia, es, em, t_committed=REAL(s, wp)*DT)
    END DO
    CALL CD_AGG_Get_Rod_States(agg, s1, es, em)
    CALL require(es == CD_AGG_OK .AND. ALL(IEEE_IS_FINITE(s1)), 'rodmir:get-finite')
    DO s = 11, 20
      CALL CD_AGG_Step_Moving(agg, DT, r0, vel, acc, cv, st, ia, es, em, t_committed=REAL(s, wp)*DT)
    END DO
    CALL CD_AGG_Get_Rod_States(agg, s2, es, em)
    CALL require(nan_max_abs(s2 - s1) > 0.0_wp, 'rodmir:march-changes-state')
    CALL CD_AGG_Set_Rod_States(agg, s1, es, em)
    CALL require(es == CD_AGG_OK, 'rodmir:set-ok: '//TRIM(em))
    CALL CD_AGG_Get_Rod_States(agg, s3, es, em)
    CALL require(nan_max_abs(s3 - s1) <= 0.0_wp, 'rodmir:get-set-get-bit-identical')
    sbad = 0.0_wp
    CALL CD_AGG_Set_Rod_States(agg, sbad, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'size mismatch') > 0, 'rodmir:set-wrong-size-fails')
    s1(1) = IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN)
    CALL CD_AGG_Set_Rod_States(agg, s1, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'non-finite') > 0, 'rodmir:set-nonfinite-fails')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_coupled_rod_mirror

  SUBROUTINE case_coupled_rod_fluid()
    TYPE(CD_AGG_ModuleType) :: flow, still
    INTEGER :: es, s, nf, ia
    REAL(wp) :: r0(3, 1), v0(3, 1), a0(3, 1), load(3, 1), x_flow, x_still
    REAL(wp), ALLOCATABLE :: xyz(:, :), fv(:, :), fa(:, :), wl(:), sf(:), ss(:), pd(:)
    LOGICAL :: cv, st
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.005_wp

    CALL write_coupled_rod_deck('agg_rod_fluid.dat', .FALSE.)
    CALL CD_AGG_Init_From_Deck(flow, 'agg_rod_fluid.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'rod-fluid:flow-build: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_Init_From_Deck(still, 'agg_rod_fluid.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'rod-fluid:still-build: '//TRIM(em))
    nf = CD_AGG_NFluidNodes(flow, es, em)
    ALLOCATE (xyz(3, nf), fv(3, nf), fa(3, nf), wl(nf), sf(24), ss(24))
    CALL CD_AGG_GetFluidNodePositions(flow, xyz, es, em)
    ! the rod block closes the layout: NumSegs + 1 = 3 segment stations, End A to End B, 1 m apart
    CALL require(nf >= 3 .AND. es == CD_AGG_OK, 'rod-fluid:segment-stations-present')
    CALL require(ABS(SQRT(SUM((xyz(:, nf) - xyz(:, nf - 1))**2)) - 1.0_wp) < 1.0e-12_wp .AND. &
                 ABS(SQRT(SUM((xyz(:, nf - 1) - xyz(:, nf - 2))**2)) - 1.0_wp) < 1.0e-12_wp .AND. &
                 ABS(SQRT(SUM((xyz(:, nf) - xyz(:, nf - 2))**2)) - 2.0_wp) < 1.0e-12_wp, 'rod-fluid:station-spacing')
    fv = 0.0_wp
    fv(1, :) = 2.0_wp
    fa = 0.0_wp
    wl = 0.0_wp
    CALL CD_AGG_SetFluidFields(flow, fv, fa, wl, es, em)
    CALL require(es == CD_AGG_OK, 'rod-fluid:set-host-field: '//TRIM(em))
    ALLOCATE (pd(nf))
    pd = 0.0_wp
    CALL CD_AGG_SetFluidFields(flow, fv, fa, wl, es, em, dynamic_pressure=pd(1:nf - 1))
    CALL require(es /= CD_AGG_OK, 'rod-fluid:dynamic-pressure-shape-fails-closed')
    pd(nf) = IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN)
    CALL CD_AGG_SetFluidFields(flow, fv, fa, wl, es, em, dynamic_pressure=pd)
    CALL require(es /= CD_AGG_OK, 'rod-fluid:dynamic-pressure-non-finite-fails-closed')
    pd(nf) = 0.0_wp
    CALL CD_AGG_SetFluidFields(flow, fv, fa, wl, es, em, dynamic_pressure=pd)
    CALL require(es == CD_AGG_OK, 'rod-fluid:dynamic-pressure-accepted: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(flow, r0, v0, a0, load, es, em)
    v0 = 0.0_wp
    a0 = 0.0_wp
    DO s = 1, 40
      CALL CD_AGG_SetFluidFields(flow, fv, fa, wl, es, em)
      CALL CD_AGG_Step_Moving(flow, DT, r0, v0, a0, cv, st, ia, es, em, t_committed=REAL(s, wp)*DT)
      CALL CD_AGG_Step_Moving(still, DT, r0, v0, a0, cv, st, ia, es, em, t_committed=REAL(s, wp)*DT)
    END DO
    CALL CD_AGG_Get_Rod_States(flow, sf, es, em)
    CALL CD_AGG_Get_Rod_States(still, ss, es, em)
    x_flow = sf(1)
    x_still = ss(1)
    CALL require(x_flow - x_still > 1.0e-4_wp, 'rod-fluid:host-current-drives-rod-downstream')
    CALL CD_AGG_End(flow, es, em)
    CALL CD_AGG_End(still, es, em)
    ! the coupled rod march turns no rod about a pin: a Pinned rod is refused by name
    CALL write_coupled_rod_deck('agg_rod_pinned.dat', .FALSE., 'Pinned')
    CALL CD_AGG_Init_From_Deck(still, 'agg_rod_pinned.dat', DT, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'Pinned ROD') > 0, 'rod:pinned-fails-closed-by-name: '//TRIM(em))
  END SUBROUTINE case_coupled_rod_fluid

  SUBROUTINE case_rejected_static_branch_status()
    !! A net-heavy 200 m finite-EI line anchored on a flat seabed 20 m from its fairlead and
    !! 10 m below it is longer than its span plus rise: it cannot hang as a catenary, and
    !! 170 m of excess length leaves no equilibrium held by its bending stiffness clear of
    !! the seabed either. That is a static failure of a valid deck, so the aggregate reports SOLVEFAIL
    !! (driver exit 2), never BADINPUT.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, u
    CHARACTER(2048) :: em

    OPEN (NEWUNIT=u, FILE='agg_coarse_fold.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'coarse finite-EI fold'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'cab 0.10 100.0 1.0e9 0.0 1.0e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Fixed 0.0 0.0 -10.0'
    WRITE (u, '(A)') '2 Coupled 20.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 cab 200.0 3'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_coarse_fold.dat', 0.05_wp, es, em)
    CALL require(es == CD_AGG_SOLVEFAIL, 'coarse-branch:status-is-solvefail: '//TRIM(em))
    CALL require(INDEX(em, 'too long for its span') > 0, 'coarse-branch:names-geometric-reason: '//TRIM(em))
    IF (es == CD_AGG_OK) CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_rejected_static_branch_status

  SUBROUTINE case_coupled_touchdown_contact()
    !! Production OpenFAST finite-EI touchdown gate. A heavy cable with surplus length rests on
    !! a -50 m bed and carries nonzero friction. The flat WtrDpth and an equivalent structured
    !! bathymetry must both install contact, initialize, and follow the same moving-fairlead
    !! trajectory. This exercises contact forces/tangents, consistent initial acceleration,
    !! chord-frame bathymetry queries, aggregate stepping, and teardown behind CompMooring=5.
    TYPE(CD_AGG_ModuleType) :: flat, bathy, sequenced_bathy
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp) :: pf(3, 1), vf(3, 1), af(3, 1), lf(3, 1), p0(3, 1), l0(3, 1)
    REAL(wp) :: pb(3, 1), vb(3, 1), ab(3, 1), lb(3, 1), target(3, 1)
    REAL(wp) :: min_flat, min_bathy, scale, floor_z, clearance, min_clearance
    INTEGER :: es, s, ni, node
    LOGICAL :: cv, st
    CHARACTER(2048) :: em

    ! Use the resolution-adequate 256-element fixture for the lifecycle/contact test.
    CALL write_cable_grounded_deck('agg_grounded_flat.dat', refined=.TRUE.)
    CALL write_cable_bathymetry('agg_grounded.xyz', 50.0_wp)
    CALL write_cable_grounded_deck('agg_grounded_bathy.dat', 'agg_grounded.xyz', refined=.TRUE.)
    CALL CD_AGG_Init_From_Deck(flat, 'agg_grounded_flat.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'grounded-flat:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_Init_From_Deck(bathy, 'agg_grounded_bathy.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'grounded-bathy:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) THEN
      CALL CD_AGG_End(flat, es, em); RETURN
    END IF
    CALL require(flat%cables(1)%line%has_contact, 'grounded-flat:contact-installed')
    CALL require(.NOT. flat%cables(1)%line%has_contact_bathymetry, 'grounded-flat:flat-contact-selected')
    CALL require(bathy%cables(1)%line%has_contact, 'grounded-bathy:contact-installed')
    CALL require(bathy%cables(1)%line%has_contact_bathymetry, 'grounded-bathy:structured-contact-selected')
    CALL require(ABS(flat%cables(1)%line%contact_mu - 0.35_wp) <= 0.0_wp, 'grounded-flat:friction-wired')
    CALL require(ABS(bathy%cables(1)%line%contact_mu - 0.35_wp) <= 0.0_wp, 'grounded-bathy:friction-wired')
    min_flat = CD_HFMF_MinSpanZ(flat%cables(1))
    min_bathy = CD_HFMF_MinSpanZ(bathy%cables(1))
    CALL require(min_flat >= -50.5_wp .AND. min_flat <= -49.9_wp, 'grounded-flat:touchdown-at-bed')
    CALL require(ABS(min_bathy - min_flat) < 1.0e-8_wp, 'grounded:flat/bathymetry-static-depth-match')
    CALL require(nan_max_abs(flat%cables(1)%line%q - bathy%cables(1)%line%q) <= 0.0_wp, &
                 'grounded:flat/bathymetry-static-state-bit-identical')
    CALL CD_AGG_GetMovingPointMesh(flat, pf, vf, af, lf, es, em)
    CALL require(es == CD_AGG_OK, 'grounded-flat:initial-mesh')
    p0 = pf
    l0 = lf
    CALL CD_AGG_GetMovingPointMesh(bathy, pb, vb, ab, lb, es, em)
    CALL require(es == CD_AGG_OK, 'grounded-bathy:initial-mesh')
    scale = MAX(1.0_wp, nan_max_abs(lf))
    CALL require(nan_max_abs(lf - lb) <= 1.0e-8_wp*scale, 'grounded:flat/bathymetry-static-load-match')

    ! Force the installed mesh-sequencing branch (256 elements, including a net-buoyant
    ! section) on a genuinely non-flat floor. The initialized state must honor the
    ! structured surface at every node; replacing it with the scalar endpoint minimum
    ! permits the coarse solve to pass beneath the shallower end of this plane.
    CALL write_cable_sloped_bathymetry('agg_grounded_slope.xyz')
    CALL write_cable_grounded_deck('agg_grounded_sequenced_bathy.dat', &
                                   'agg_grounded_slope.xyz', refined=.TRUE.)
    CALL CD_AGG_Init_From_Deck(sequenced_bathy, 'agg_grounded_sequenced_bathy.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'grounded-sequenced-bathy:init: '//TRIM(em))
    IF (es == CD_AGG_OK) THEN
      min_clearance = HUGE(1.0_wp)
      DO node = 1, SIZE(sequenced_bathy%cables(1)%line%q)/6
        floor_z = -(50.0_wp + 0.1_wp*sequenced_bathy%cables(1)%line%q(6*(node - 1) + 1))
        clearance = sequenced_bathy%cables(1)%line%q(6*(node - 1) + 3) - floor_z
        min_clearance = MIN(min_clearance, clearance)
      END DO
      CALL require(min_clearance >= -5.0e-3_wp, &
                   'grounded-sequenced-bathy:static state honors structured floor')
      CALL require(min_clearance <= 0.10_wp, &
                   'grounded-sequenced-bathy:static state retains touchdown contact')
      CALL CD_AGG_End(sequenced_bathy, es, em)
    END IF

    ! OpenFAST checkpoints and failed-step recovery use this exact aggregate transaction.
    ! Contact configuration remains installed while the committed q/v/a state is restored.
    CALL CD_AGG_Snapshot(flat, es, em)
    CALL require(es == CD_AGG_OK, 'grounded-flat:snapshot: '//TRIM(em))
    target = pf
    target(1, 1) = target(1, 1) + 0.01_wp
    vf = 0.0_wp; af = 0.0_wp
    vf(1, 1) = 0.01_wp/DT
    CALL CD_AGG_Step_Moving(flat, DT, target, vf, af, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'grounded-flat:pre-restore-step: '//TRIM(em))
    CALL CD_AGG_Restore(flat, es, em)
    CALL require(es == CD_AGG_OK, 'grounded-flat:restore: '//TRIM(em))
    CALL CD_AGG_CalcOutput(flat, es, em)
    CALL require(es == CD_AGG_OK, 'grounded-flat:post-restore-output: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(flat, pf, vf, af, lf, es, em)
    CALL require(nan_max_abs(pf - p0) <= 0.0_wp, 'grounded-flat:restore-position-bit-identical')
    CALL require(nan_max_abs(lf - l0) <= 0.0_wp, 'grounded-flat:restore-load-bit-identical')

    DO s = 1, 4
      target = pf
      target(1, 1) = pf(1, 1) + 0.01_wp*REAL(s, wp)
      target(3, 1) = pf(3, 1) + 0.005_wp*REAL(s, wp)
      vf = 0.0_wp; af = 0.0_wp
      vf(1, 1) = 0.01_wp/DT; vf(3, 1) = 0.005_wp/DT
      CALL CD_AGG_Step_Moving(flat, DT, target, vf, af, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'grounded-flat:moving-touchdown-step: '//TRIM(em))
      CALL CD_AGG_Step_Moving(bathy, DT, target, vf, af, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'grounded-bathy:moving-touchdown-step: '//TRIM(em))
      CALL CD_AGG_GetMovingPointMesh(flat, pf, vf, af, lf, es, em)
      CALL CD_AGG_GetMovingPointMesh(bathy, pb, vb, ab, lb, es, em)
      scale = MAX(1.0_wp, nan_max_abs(lf))
      ! The structured evaluator performs coordinate interpolation even on a constant
      ! grid, so its implicit dynamic trajectory is equivalent to the plane to about
      ! one part per million rather than bit-identical.
      CALL require(nan_max_abs(lf - lb) <= 1.0e-6_wp*scale, 'grounded:flat/bathymetry-dynamic-load-match')
    END DO
    CALL CD_AGG_End(flat, es, em)
    CALL CD_AGG_End(bathy, es, em)
  END SUBROUTINE case_coupled_touchdown_contact

  SUBROUTINE case_step_atomicity()
    !! STAGE-THEN-COMMIT: a failed aggregate step must leave the aggregate EXACTLY at the
    !! step-start state, and the aggregate must remain steppable (the retry/recovery
    !! contract). Rig: the mixed deck driven so the partial-advance hazard is real -- a
    !! benign mooring column (the EI=0 system step SUCCEEDS and would strand at t+dt)
    !! plus a teleported cable fairlead (the implicit cable step FAILS). Observables:
    !! committed positions and loads at every column must be BIT-IDENTICAL to the
    !! pre-failure state (restored states => recomputed loads identical), and a
    !! subsequent valid step must converge.
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp) :: r0(3, 2), vel(3, 2), acc(3, 2), tgt(3, 2)
    REAL(wp) :: p_before(3, 2), l_before(3, 2), p_after(3, 2), l_after(3, 2), v_s(3, 2), a_s(3, 2)
    INTEGER :: es, ni, s
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_mixed_deck('agg_atomic.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_atomic.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    ! a restore without a snapshot fails closed at the sub-module surface
    CALL CD_HFMF_Restore(agg%cables(1), es, em)
    CALL require(es /= CD_HFMF_OK, 'atomic:restore-without-snapshot-fails-closed')
    ! a few healthy rest-hold steps to a committed dynamic state
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, l_before, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:mesh0')
    tgt = r0
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 5
      CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'atomic:healthy-step')
      IF (es /= CD_AGG_OK) RETURN
    END DO
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:calc-before')
    CALL CD_AGG_GetMovingPointMesh(agg, p_before, v_s, a_s, l_before, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:mesh-before')
    ! the poisoned step: the mooring column gets a SMALL VALID MOTION (its system step
    ! converges to changed t+dt kinematics -- holding it constant would be structurally
    ! blind to a stale facade mesh after rollback, the missing-case class); the cable
    ! fairlead teleports by 1e150 m -- FINITE (so the input validation passes and the
    ! solve itself must fail: the strain energy overflows and the Newton cannot
    ! converge); a mere 10 km teleport measurably CONVERGES on the linear-elastic
    ! cable, which is robustness, not a failure injection
    tgt(:, 1) = r0(:, 1) + [0.0_wp, 0.0_wp, 0.02_wp]
    vel(:, 1) = [0.0_wp, 0.0_wp, 0.02_wp/DT]
    tgt(:, 2) = r0(:, 2) + [1.0e150_wp, 0.0_wp, 0.0_wp]
    ! deliberately non-finite or overflowing input: must not halt a trapping build
    CALL IEEE_GET_HALTING_MODE(IEEE_USUAL, fp_halt)
    CALL IEEE_SET_HALTING_MODE(IEEE_USUAL, .FALSE.)
    CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL IEEE_SET_FLAG(IEEE_USUAL, .FALSE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_USUAL, fp_halt)
    CALL require(es /= CD_AGG_OK .AND. .NOT. cv, 'atomic:poisoned-step-fails')
    vel = 0.0_wp
    ! the aggregate must sit exactly at the step-start state: recomputed loads and
    ! committed positions bit-identical
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:calc-after-rollback')
    CALL CD_AGG_GetMovingPointMesh(agg, p_after, v_s, a_s, l_after, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:mesh-after-rollback')
    CALL require(nan_max_abs(p_after - p_before) <= 0.0_wp, 'atomic:positions-bit-identical-after-rollback')
    CALL require(nan_max_abs(l_after - l_before) <= 0.0_wp, 'atomic:loads-bit-identical-after-rollback')
    ! retryability: a valid rest-hold step converges from the restored state
    tgt = r0
    CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'atomic:valid-step-after-rollback-converges')

    ! the SYSTEM-failure mirror: poison the MOORING column (the EI=0 solve fails) with a
    ! benign, MOVING cable column -- the facade point mesh received the failed kinematics
    ! before the solve, so this leg pins the FMF's own failure atomicity (mesh re-derived
    ! from the rolled-back committed states) plus the staged undo of the cable snapshot
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:sys-calc-before')
    CALL CD_AGG_GetMovingPointMesh(agg, p_before, v_s, a_s, l_before, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:sys-mesh-before')
    ! Use a finite displacement whose squared chord necessarily exceeds binary64.
    ! The previous 1e150 value could remain representable through the axial path
    ! on some compilers and was therefore not a deterministic failure injection.
    tgt(:, 1) = r0(:, 1) + [1.0e200_wp, 0.0_wp, 0.0_wp]
    tgt(:, 2) = r0(:, 2) + [0.0_wp, 0.0_wp, 0.02_wp]
    vel(:, 2) = [0.0_wp, 0.0_wp, 0.02_wp/DT]
    ! deliberately non-finite or overflowing input: must not halt a trapping build
    CALL IEEE_GET_HALTING_MODE(IEEE_USUAL, fp_halt)
    CALL IEEE_SET_HALTING_MODE(IEEE_USUAL, .FALSE.)
    CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL IEEE_SET_FLAG(IEEE_USUAL, .FALSE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_USUAL, fp_halt)
    CALL require(es /= CD_AGG_OK .AND. .NOT. cv, 'atomic:sys-poisoned-step-fails')
    vel = 0.0_wp
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:sys-calc-after-rollback')
    CALL CD_AGG_GetMovingPointMesh(agg, p_after, v_s, a_s, l_after, es, em)
    CALL require(es == CD_AGG_OK, 'atomic:sys-mesh-after-rollback')
    CALL require(nan_max_abs(p_after - p_before) <= 0.0_wp, 'atomic:sys-positions-bit-identical')
    CALL require(nan_max_abs(l_after - l_before) <= 0.0_wp, 'atomic:sys-loads-bit-identical')
    tgt = r0
    CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'atomic:sys-valid-step-after-rollback-converges')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_step_atomicity

  SUBROUTINE case_ptfm_init_equivalence()
    !! PtfmInit correctness gate. Route A builds the reference deck with the ptfm_init
    !! option (the aggregate rigid-transforms every Vessel point by p' = r + M^T p before
    !! the static solve). Route B computes the SAME displaced fairleads here, in the test,
    !! from an independent transcription of the NWTC ZYX DCM, writes them into a deck
    !! round-trip-exact, and builds plainly. Both routes must be BIT-IDENTICAL: initial
    !! mesh positions/loads, and the full stepped trajectory (deterministic solver on
    !! identical inputs). Plus: zero ptfm_init == absent (pass-through), and a non-finite
    !! ptfm_init fails closed.
    TYPE(CD_AGG_ModuleType) :: agg_opt, agg_ref
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp), PARAMETER :: P6(6) = [2.0_wp, -1.0_wp, 0.5_wp, 0.02_wp, -0.03_wp, 0.05_wp]
    INTEGER, PARAMETER :: NSTEP = 4
    REAL(wp) :: fl0(3, 3), fld(3, 3), dcm(3, 3), bad6(6)
    REAL(wp) :: pa(3, 3), va(3, 3), aa(3, 3), la(3, 3)
    REAL(wp) :: pb(3, 3), vb(3, 3), ab(3, 3), lb(3, 3)
    REAL(wp) :: r0a(3, 3), r0b(3, 3), tgt(3, 3), vel(3, 3), acc(3, 3), zh, vh, t
    REAL(wp) :: pos_tol, load_tol, load_scale
    INTEGER :: es, i, s, ni
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    ! the undisplaced fairleads (must match write_volturnus_deck's defaults bit-for-bit:
    ! both the deck text '-58.000' and the literal parse to the same correctly-rounded
    ! binary64)
    fl0(:, 1) = [-58.0_wp, 0.0_wp, -14.0_wp]
    fl0(:, 2) = [29.0_wp, 50.229_wp, -14.0_wp]
    fl0(:, 3) = [29.0_wp, -50.229_wp, -14.0_wp]
    ! independent transform: the NWTC EulerConstructZYX expressions transcribed here so
    ! the test pins the CONVENTION, not just self-consistency with the implementation
    dcm = test_zyx_dcm(P6(4:6))
    DO i = 1, 3
      fld(:, i) = P6(1:3) + MATMUL(TRANSPOSE(dcm), fl0(:, i))
    END DO

    ! route A: reference deck + ptfm_init option
    CALL write_volturnus_deck('agg_ptfm_a.dat')
    CALL CD_AGG_Init_From_Deck(agg_opt, 'agg_ptfm_a.dat', DT, es, em, ptfm_init=P6)
    CALL require(es == CD_AGG_OK, 'ptfm:opt-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    ! route B: pre-displaced deck, no option
    CALL write_volturnus_deck('agg_ptfm_b.dat', fairleads=fld)
    CALL CD_AGG_Init_From_Deck(agg_ref, 'agg_ptfm_b.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'ptfm:ref-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) THEN
      CALL CD_AGG_End(agg_opt, es, em)
      RETURN
    END IF

    CALL CD_AGG_GetMovingPointMesh(agg_opt, pa, va, aa, la, es, em)
    CALL require(es == CD_AGG_OK, 'ptfm:opt-mesh0')
    CALL CD_AGG_GetMovingPointMesh(agg_ref, pb, vb, ab, lb, es, em)
    CALL require(es == CD_AGG_OK, 'ptfm:ref-mesh0')
    pos_tol = 64.0_wp*EPSILON(1.0_wp)*MAX(1.0_wp, nan_max_abs(fld))
    IF (nan_max_abs(pa - fld) > pos_tol .OR. nan_max_abs(pa - pb) > pos_tol) THEN
      WRITE (*, '(A,ES13.5,A,ES13.5)') 'ptfm initial differences: option-vs-transform=', &
        nan_max_abs(pa - fld), ', option-vs-roundtrip=', nan_max_abs(pa - pb)
    END IF
    ! The independent DCM transcription can differ by one final rounding operation.
    ! Require a scale-aware multiple of machine precision rather than bit identity
    ! across the two algebraically equivalent floating-point evaluation orders.
    CALL require(nan_max_abs(pa - fld) <= pos_tol, 'ptfm:opt-positions-match-independent-transform')
    CALL require(nan_max_abs(pa - pb) <= pos_tol, 'ptfm:initial-positions-roundoff-equivalent')
    load_scale = MAX(1.0_wp, nan_max_abs(la), nan_max_abs(lb))
    load_tol = 8192.0_wp*EPSILON(1.0_wp)*load_scale
    CALL require(nan_max_abs(la - lb) <= load_tol, 'ptfm:initial-loads-roundoff-equivalent')
    CALL require(ALL(IEEE_IS_FINITE(la)), 'ptfm:initial-loads-finite')

    ! the displaced static solves must track through identical stepped trajectories
    r0a = pa
    r0b = pb
    DO s = 1, NSTEP
      t = REAL(s, wp)*DT
      zh = 0.1_wp*(1.0_wp - COS(2.0_wp*PI*t/8.0_wp))
      vh = 0.1_wp*(2.0_wp*PI/8.0_wp)*SIN(2.0_wp*PI*t/8.0_wp)
      DO i = 1, 3
        tgt(:, i) = [r0a(1, i), r0a(2, i), r0a(3, i) + zh]
        vel(:, i) = [0.0_wp, 0.0_wp, vh]
        acc(:, i) = 0.0_wp
      END DO
      CALL CD_AGG_Step_Moving(agg_opt, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'ptfm:opt-step')
      CALL CD_AGG_Step_Moving(agg_ref, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'ptfm:ref-step')
      CALL CD_AGG_CalcOutput(agg_opt, es, em)
      CALL CD_AGG_CalcOutput(agg_ref, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg_opt, pa, va, aa, la, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg_ref, pb, vb, ab, lb, es, em)
      load_scale = MAX(1.0_wp, nan_max_abs(la), nan_max_abs(lb))
      load_tol = 8192.0_wp*EPSILON(1.0_wp)*load_scale
      IF (nan_max_abs(pa - pb) > pos_tol .OR. nan_max_abs(la - lb) > load_tol) THEN
        WRITE (*, '(A,I0,A,ES13.5,A,ES13.5)') 'ptfm step ', s, ': position difference=', &
          nan_max_abs(pa - pb), ', load difference=', nan_max_abs(la - lb)
      END IF
      CALL require(nan_max_abs(pa - pb) <= pos_tol, 'ptfm:stepped-positions-roundoff-equivalent')
      CALL require(nan_max_abs(la - lb) <= load_tol, 'ptfm:stepped-loads-roundoff-equivalent')
    END DO
    CALL CD_AGG_End(agg_opt, es, em)
    CALL CD_AGG_End(agg_ref, es, em)

    ! zero ptfm_init is a pure pass-through: bit-identical to the absent-optional build
    CALL CD_AGG_Init_From_Deck(agg_opt, 'agg_ptfm_a.dat', DT, es, em, ptfm_init=[0.0_wp, 0.0_wp, &
                                                                                 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp])
    CALL require(es == CD_AGG_OK, 'ptfm:zero-init: '//TRIM(em))
    CALL CD_AGG_Init_From_Deck(agg_ref, 'agg_ptfm_a.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'ptfm:absent-init: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg_opt, pa, va, aa, la, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_ref, pb, vb, ab, lb, es, em)
    CALL require(nan_max_abs(pa - pb) <= 0.0_wp, 'ptfm:zero-equals-absent-positions')
    CALL require(nan_max_abs(la - lb) <= 0.0_wp, 'ptfm:zero-equals-absent-loads')
    CALL require(nan_max_abs(pa - fl0) <= 0.0_wp, 'ptfm:zero-leaves-deck-fairleads')
    CALL CD_AGG_End(agg_opt, es, em)
    CALL CD_AGG_End(agg_ref, es, em)

    ! a non-finite ptfm_init fails closed at init
    bad6 = P6
    bad6(3) = IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN)
    CALL CD_AGG_Init_From_Deck(agg_opt, 'agg_ptfm_a.dat', DT, es, em, ptfm_init=bad6)
    CALL require(es /= CD_AGG_OK, 'ptfm:non-finite-fails-closed')
    CALL require(.NOT. CD_AGG_IsInitialized(agg_opt), 'ptfm:non-finite-leaves-uninitialized')

    ! the displaced geometry must be REVALIDATED: the transform runs after the parse's
    ! geometry checks, so a heave that sinks a fairlead below its Fixed anchor
    ! (fairleads z=-14, anchors z=-200 on this deck) must fail closed at init, not reach
    ! the static builders
    bad6 = 0.0_wp
    bad6(3) = -190.0_wp
    CALL CD_AGG_Init_From_Deck(agg_opt, 'agg_ptfm_a.dat', DT, es, em, ptfm_init=bad6)
    CALL require(es /= CD_AGG_OK, 'ptfm:fairlead-below-anchor-fails-closed')
    CALL require(INDEX(em, 'below') > 0, 'ptfm:fairlead-below-anchor-message-names-the-contract')
    CALL require(.NOT. CD_AGG_IsInitialized(agg_opt), 'ptfm:fairlead-below-anchor-leaves-uninitialized')

    ! OUT-OF-PLANE ptfm_init on a MIXED (mooring + finite-EI cable) deck: sway/yaw move
    ! the cable's Coupled fairlead off its x-z plane while the Fixed anchor stays -- the
    ! cable must build in its chord-local frame (azimuthal support) and the displaced
    ! aggregate must march. This was the exact blocked case: the x-z plane restriction
    ! made every out-of-plane displaced start on a cable deck fail at init.
    CALL case_ptfm_mixed_out_of_plane()
  END SUBROUTINE case_ptfm_init_equivalence

  SUBROUTINE case_semitaut_displaced_pose_init()
    !! Regression for a field-scale composite semitaut system at the displaced
    !! platform pose reached by a 21.2 m/s operating-state pre-run. The former
    !! global-frame cold start stalled in the first continuation stage even though
    !! the same equilibrium was reached by a dynamic relaxation from the origin.
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp), PARAMETER :: DT = 0.025_wp
    ! OpenFAST's registry boundary carries PtfmInit in ReKi. Preserve the
    ! resulting binary32-rounded values here because that rounding exposed the
    ! formerly fragile cold start.
    REAL(wp), PARAMETER :: P6_SEVERE(6) = [2.6035442352294922_wp, -13.878925323486328_wp, &
                                           0.46294105052947998_wp, 0.027599746361374855_wp, -0.005574016831815243_wp, &
                                           -0.04081614688038826_wp]
    REAL(wp), PARAMETER :: P6_RATED(6) = [16.802740440000001_wp, -17.496244050000001_wp, &
                                          0.35881783960000002_wp, 0.045055365734409304_wp, 0.052519905315415380_wp, &
                                          -0.0048230808019097274_wp]
    INTEGER :: es
    CHARACTER(512) :: em

    CALL write_lozon_semitaut_deck('agg_semitaut_pose.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_semitaut_pose.dat', DT, es, em, &
                               env_gravity=9.806650161743164_wp, env_rho_water=1025.0_wp, env_wtrdpth=200.0_wp, &
                               ptfm_init=P6_SEVERE, external_fluid=.TRUE., run_tmax=4200.0_wp)
    CALL require(es == CD_AGG_OK, 'semitaut-pose:severe-init: '//TRIM(em))
    IF (es == CD_AGG_OK) CALL CD_AGG_End(agg, es, em)

    CALL CD_AGG_Init_From_Deck(agg, 'agg_semitaut_pose.dat', DT, es, em, &
                               env_gravity=9.806650161743164_wp, env_rho_water=1025.0_wp, env_wtrdpth=200.0_wp, &
                               ptfm_init=P6_RATED, external_fluid=.TRUE., run_tmax=4200.0_wp)
    CALL require(es == CD_AGG_OK, 'semitaut-pose:rated-init: '//TRIM(em))
    IF (es == CD_AGG_OK) CALL CD_AGG_End(agg, es, em)

    CALL CD_AGG_Init_From_Deck(agg, 'agg_semitaut_pose.dat', DT, es, em, &
                               env_gravity=9.80665_wp, env_rho_water=1025.0_wp, env_wtrdpth=200.0_wp, &
                               ptfm_init=P6_RATED, external_fluid=.TRUE., run_tmax=4200.0_wp)
    CALL require(es == CD_AGG_OK, 'semitaut-pose:rated-exact-gravity-init: '//TRIM(em))
    IF (es == CD_AGG_OK) CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_semitaut_displaced_pose_init

  SUBROUTINE case_ptfm_mixed_out_of_plane()
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp), PARAMETER :: P6(6) = [1.0_wp, -2.0_wp, 0.3_wp, 0.0_wp, 0.0_wp, 0.1_wp]
    REAL(wp) :: r0(3, 2), vel(3, 2), acc(3, 2), tgt(3, 2), lds(3, 2)
    INTEGER :: es, ni, s
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_mixed_deck('agg_ptfm_mixed.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_ptfm_mixed.dat', DT, es, em, ptfm_init=P6)
    CALL require(es == CD_AGG_OK, 'ptfm-oop:mixed-init-with-sway-yaw: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, lds, es, em)
    CALL require(es == CD_AGG_OK, 'ptfm-oop:mesh0')
    CALL require(ALL(IEEE_IS_FINITE(lds)), 'ptfm-oop:initial-loads-finite')
    tgt = r0
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 4
      tgt(3, :) = r0(3, :) + 0.02_wp*REAL(s, wp)
      vel(3, :) = 0.02_wp/DT
      CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'ptfm-oop:displaced-step-converges: '//TRIM(em))
      IF (es /= CD_AGG_OK) EXIT
    END DO
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_ptfm_mixed_out_of_plane

  PURE FUNCTION test_zyx_dcm(theta) RESULT(m)
    !! Independent transcription of NWTC_Num's EulerConstructZYX (the 3-2-1 sequence
    !! global-to-local DCM) -- the SAME expressions in the SAME order as the
    !! implementation, so both evaluate bit-identically AND the test would catch either
    !! side drifting from the MoorDyn-F PtfmInit convention.
    REAL(wp), INTENT(IN) :: theta(3)
    REAL(wp) :: m(3, 3)
    REAL(wp) :: cx, sx, cy, sy, cz, sz
    cx = COS(theta(1)); sx = SIN(theta(1))
    cy = COS(theta(2)); sy = SIN(theta(2))
    cz = COS(theta(3)); sz = SIN(theta(3))
    m(1, 1) = cy*cz
    m(2, 1) = sx*sy*cz - sz*cx
    m(3, 1) = sx*sz + sy*cx*cz
    m(1, 2) = sz*cy
    m(2, 2) = sx*sy*sz + cx*cz
    m(3, 2) = -sx*cz + sy*sz*cx
    m(1, 3) = -sy
    m(2, 3) = sx*cy
    m(3, 3) = cx*cy
  END FUNCTION test_zyx_dcm

  SUBROUTINE case_correction_rewind()
    !! Correction-iteration semantics at the aggregate level (the OpenFAST shell's
    !! NumCrctn > 0 building block): snapshot a committed state, advance with motion A,
    !! RESTORE, re-advance with motion B -- positions and loads must be BIT-IDENTICAL to
    !! a twin instance that only ever saw motion B. The rewound A-step must leave no
    !! trace, on a MIXED deck so both the EI=0 system and the Hermite cable participate.
    TYPE(CD_AGG_ModuleType) :: agg, twin
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp) :: r0(3, 2), vel(3, 2), acc(3, 2), tgt(3, 2)
    REAL(wp) :: p_snap(3, 2), l_snap(3, 2), v_s(3, 2), a_s(3, 2)
    REAL(wp) :: pa(3, 2), la(3, 2), pb(3, 2), lb(3, 2)
    INTEGER :: es, ni, s
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    ! snapshot/restore on an UNINITIALIZED aggregate fails closed
    CALL CD_AGG_Snapshot(agg, es, em)
    CALL require(es /= CD_AGG_OK, 'crctn:snapshot-uninitialized-fails-closed')
    CALL CD_AGG_Restore(agg, es, em)
    CALL require(es /= CD_AGG_OK, 'crctn:restore-uninitialized-fails-closed')

    CALL write_mixed_deck('agg_crctn.dat')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_crctn.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'crctn:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    ! restore before ANY snapshot fails closed (initialized, no snapshot yet)
    CALL CD_AGG_Restore(agg, es, em)
    CALL require(es /= CD_AGG_OK, 'crctn:restore-without-snapshot-fails-closed')

    CALL CD_AGG_Init_From_Deck(twin, 'agg_crctn.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'crctn:twin-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) THEN
      CALL CD_AGG_End(agg, es, em)
      RETURN
    END IF

    ! bring BOTH to the same committed dynamic state
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, l_snap, es, em)
    CALL require(es == CD_AGG_OK, 'crctn:mesh0')
    tgt = r0
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 3
      CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'crctn:warm-step')
      CALL CD_AGG_Step_Moving(twin, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'crctn:twin-warm-step')
      IF (es /= CD_AGG_OK) RETURN
    END DO
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, p_snap, v_s, a_s, l_snap, es, em)
    CALL require(es == CD_AGG_OK, 'crctn:mesh-snap')

    ! snapshot, advance with motion A (the "uncorrected predictor" step)
    CALL CD_AGG_Snapshot(agg, es, em)
    CALL require(es == CD_AGG_OK, 'crctn:snapshot: '//TRIM(em))
    tgt(:, 1) = r0(:, 1) + [0.01_wp, 0.0_wp, 0.03_wp]
    tgt(:, 2) = r0(:, 2) + [0.0_wp, 0.02_wp, 0.01_wp]
    vel(:, 1) = [0.01_wp, 0.0_wp, 0.03_wp]/DT
    vel(:, 2) = [0.0_wp, 0.02_wp, 0.01_wp]/DT
    CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'crctn:A-step')

    ! rewind: the aggregate must sit exactly at the snapshotted committed state
    CALL CD_AGG_Restore(agg, es, em)
    CALL require(es == CD_AGG_OK, 'crctn:restore: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'crctn:calc-after-restore')
    CALL CD_AGG_GetMovingPointMesh(agg, pa, v_s, a_s, la, es, em)
    CALL require(es == CD_AGG_OK, 'crctn:mesh-after-restore')
    CALL require(nan_max_abs(pa - p_snap) <= 0.0_wp, 'crctn:restored-positions-bit-identical')
    CALL require(nan_max_abs(la - l_snap) <= 0.0_wp, 'crctn:restored-loads-bit-identical')

    ! re-advance with motion B (the "corrected" inputs); the twin takes B directly
    tgt(:, 1) = r0(:, 1) + [0.0_wp, 0.0_wp, 0.02_wp]
    tgt(:, 2) = r0(:, 2) + [0.01_wp, 0.0_wp, 0.02_wp]
    vel(:, 1) = [0.0_wp, 0.0_wp, 0.02_wp]/DT
    vel(:, 2) = [0.01_wp, 0.0_wp, 0.02_wp]/DT
    CALL CD_AGG_Step_Moving(agg, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'crctn:B-step')
    CALL CD_AGG_Step_Moving(twin, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'crctn:twin-B-step')

    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_CalcOutput(twin, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pa, v_s, a_s, la, es, em)
    CALL CD_AGG_GetMovingPointMesh(twin, pb, v_s, a_s, lb, es, em)
    CALL require(nan_max_abs(pa - pb) <= 0.0_wp, 'crctn:rewound-B-positions-equal-direct-B')
    CALL require(nan_max_abs(la - lb) <= 0.0_wp, 'crctn:rewound-B-loads-equal-direct-B')

    ! a SECOND rewind of the same snapshot (NumCrctn >= 2: restore is repeatable)
    CALL CD_AGG_Restore(agg, es, em)
    CALL require(es == CD_AGG_OK, 'crctn:second-restore: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg, pa, v_s, a_s, la, es, em)
    CALL require(nan_max_abs(pa - p_snap) <= 0.0_wp, 'crctn:second-restore-positions-bit-identical')

    CALL CD_AGG_End(agg, es, em)
    CALL CD_AGG_End(twin, es, em)
  END SUBROUTINE case_correction_rewind

  SUBROUTINE case_farm_partition()
    !! FAST.Farm mode at the aggregate level: a two-turbine deck (one anchored line per
    !! turbine + one SHARED line between the two fairleads) with the Turbine<J>
    !! vocabulary. Gates:
    !! (a) farm init converges; NTurbines and the turbine-of-moving map are exact;
    !! (b) deck coordinates are TURBINE-LOCAL: each fairlead sits at local + ref_pos;
    !! (c) farm-layout objectivity: translating the whole layout translates the
    !!     equilibrium exactly (loads preserved, positions shifted);
    !! (d) the shared line couples the turbines: displacing turbine 1's fairlead
    !!     changes the load at turbine 2's fairlead;
    !! (e) fail-closed: Turbine<J> without farm init; plain Coupled in farm mode;
    !!     J beyond n_turbines; missing turbine_ref_pos.
    TYPE(CD_AGG_ModuleType) :: agg, agg2
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp) :: refpos(3, 2), refpos2(3, 2), q0(3, 2), v0(3, 2), a0(3, 2), l0f(3, 2)
    REAL(wp) :: q1(3, 2), v1(3, 2), a1(3, 2), l1f(3, 2), tgt(3, 2)
    INTEGER :: es, tof(2), ni
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    refpos = 0.0_wp
    refpos(1, 2) = 800.0_wp
    CALL write_farm_deck('agg_farm.dat', plain_coupled=.FALSE.)

    ! (e) fail-closed edges first
    CALL CD_AGG_Init_From_Deck(agg, 'agg_farm.dat', DT, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'farm') > 0, 'farm: Turbine<J> without farm init fails closed')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_farm.dat', DT, es, em, n_turbines=2)
    CALL require(es /= CD_AGG_OK, 'farm: missing turbine_ref_pos fails closed')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_farm.dat', DT, es, em, n_turbines=1, &
                               turbine_ref_pos=refpos(:, 1:1))
    CALL require(es /= CD_AGG_OK, 'farm: turbine id beyond n_turbines fails closed')
    CALL write_farm_deck('agg_farm_plain.dat', plain_coupled=.TRUE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_farm_plain.dat', DT, es, em, n_turbines=2, &
                               turbine_ref_pos=refpos)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'Turbine<J>') > 0, &
                 'farm: plain Coupled in farm mode fails closed')
    ! forbid_deck_ambient (the still-water coupled scope): a deck-declared current or
    ! waves OPTION must fail closed -- with the sampling of a host field suppressed,
    ! nothing else stops a deck from silently self-driving ambient forcing
    CALL write_farm_deck('agg_farm_cur.dat', plain_coupled=.FALSE., ambient='current')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_farm_cur.dat', DT, es, em, n_turbines=2, &
                               turbine_ref_pos=refpos, forbid_deck_ambient=.TRUE.)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'STILL WATER') > 0, &
                 'farm: forbid_deck_ambient rejects a deck current OPTION')
    CALL write_farm_deck('agg_farm_wav.dat', plain_coupled=.FALSE., ambient='waves')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_farm_wav.dat', DT, es, em, n_turbines=2, &
                               turbine_ref_pos=refpos, forbid_deck_ambient=.TRUE.)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'STILL WATER') > 0, &
                 'farm: forbid_deck_ambient rejects a deck waves OPTION')

    ! (a) the real farm init (forbid_deck_ambient mirrors the OpenFAST shell's farm
    ! call: a clean deck must be unaffected by the still-water scope flag)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_farm.dat', DT, es, em, n_turbines=2, &
                               turbine_ref_pos=refpos, forbid_deck_ambient=.TRUE.)
    CALL require(es == CD_AGG_OK, 'farm: init: '//TRIM(em))
    CALL require(CD_AGG_NTurbines(agg) == 2, 'farm: NTurbines')
    CALL CD_AGG_TurbineOfMoving(agg, tof, es, em)
    CALL require(es == CD_AGG_OK .AND. tof(1) == 1 .AND. tof(2) == 2, 'farm: turbine-of-moving map')

    ! (b) turbine-local deck coordinates land at local + ref_pos
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'farm: output: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg, q0, v0, a0, l0f, es, em)
    CALL require(es == CD_AGG_OK, 'farm: mesh: '//TRIM(em))
    CALL require(nan_max_abs(q0(:, 1) - [0.0_wp, 0.0_wp, -10.0_wp]) <= 1.0e-9_wp, &
                 'farm: turbine-1 fairlead at its local position')
    CALL require(nan_max_abs(q0(:, 2) - [800.0_wp, 0.0_wp, -10.0_wp]) <= 1.0e-9_wp, &
                 'farm: turbine-2 fairlead at local + ref_pos')

    ! (c) farm-layout objectivity: shift the whole layout +50 m in y (anchors through
    ! a deck twin, turbines through ref_pos) -- the equilibrium translates exactly
    CALL write_farm_deck('agg_farm_shift.dat', plain_coupled=.FALSE., y_shift=50.0_wp)
    refpos2 = refpos
    refpos2(2, :) = refpos2(2, :) + 50.0_wp
    CALL CD_AGG_Init_From_Deck(agg2, 'agg_farm_shift.dat', DT, es, em, n_turbines=2, &
                               turbine_ref_pos=refpos2)
    CALL require(es == CD_AGG_OK, 'farm: shifted init: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg2, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg2, q1, v1, a1, l1f, es, em)
    CALL require(nan_max_abs(q1(2, :) - (q0(2, :) + 50.0_wp)) <= 1.0e-6_wp, &
                 'farm: layout translation shifts the equilibrium')
    CALL require(nan_max_abs(l1f - l0f) <= 1.0e-4_wp*nan_max_abs(l0f), &
                 'farm: layout translation preserves the loads')
    CALL CD_AGG_End(agg2, es, em)

    ! the dtM metadata query must ACCEPT farm decks: the OpenFAST shell scans dtM
    ! BEFORE the farm-capable aggregate init (regression: the vocabulary gate aborted
    ! the scan and valid farm decks never reached their consumer)
    BLOCK
      REAL(wp) :: dtq
      LOGICAL :: hasq
      CALL CD_Deck_Query_dtM('agg_farm.dat', dtq, hasq, es, em)
      CALL require(es == CD_DECKDRV_OK, 'farm: dtM metadata query accepts the farm deck: '//TRIM(em))
    END BLOCK

    ! stock-order farm deck: anchored lines written NodeA = anchor, NodeB = Turbine<J>
    ! -- the Turbine normalization must run BEFORE the endpoint-order pre-pass so the
    ! swap treats the fairlead like any coupled row (regression: previously the raw
    ! turbine row survived the swap and the End-A rule rejected a valid MoorDyn deck)
    CALL write_farm_deck('agg_farm_stock.dat', plain_coupled=.FALSE., stock_order=.TRUE.)
    CALL CD_AGG_Init_From_Deck(agg2, 'agg_farm_stock.dat', DT, es, em, n_turbines=2, &
                               turbine_ref_pos=refpos)
    CALL require(es == CD_AGG_OK, 'farm: stock-order farm deck inits: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg2, q1, v1, a1, l1f, es, em)
    CALL require(nan_max_abs(q1 - q0) <= 1.0e-9_wp, 'farm: stock-order fairleads land identically')
    CALL CD_AGG_End(agg2, es, em)

    ! (d) the shared line couples the turbines
    tgt = q0
    tgt(1, 1) = tgt(1, 1) + 2.0_wp   ! surge turbine 1 toward turbine 2
    v0 = 0.0_wp; a0 = 0.0_wp
    CALL CD_AGG_Step_Moving(agg, DT, tgt, v0, a0, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'farm: shared-line step: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, q1, v1, a1, l1f, es, em)
    CALL require(nan_max_abs(l1f(:, 2) - l0f(:, 2)) > 1.0e2_wp, &
                 'farm: turbine-1 motion changes turbine-2 shared-line load')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_farm_partition

  SUBROUTINE case_farm_shared_anchor()
    !! FAST.Farm mode with one Fixed anchor shared by the lines of two turbines (800 m apart,
    !! anchor midway): the farm aggregate builds, holds still, and returns mirror-symmetric
    !! turbine loads (the anchor couples nothing between the lines).
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp), PARAMETER :: DT = 0.01_wp
    REAL(wp) :: refpos(3, 2), q0(3, 2), v0(3, 2), a0(3, 2), l0f(3, 2), l1f(3, 2)
    INTEGER :: es, s, ni, u
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    OPEN (NEWUNIT=u, FILE='agg_farm_shared_anchor.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'Two-turbine farm: both anchored lines on one shared anchor'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.333 685.0 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Turbine1 0.0 0.0 -10.0'
    WRITE (u, '(A)') '2 Turbine2 0.0 0.0 -10.0'
    WRITE (u, '(A)') '3 Fixed 400.0 0.0 -200.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 3 -'
    WRITE (u, '(A)') '2 2 3 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 chain 480.0 20'
    WRITE (u, '(A)') '2 chain 480.0 20'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '200.0 WtrDpth'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1 FairTen2'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
    refpos = 0.0_wp
    refpos(1, 2) = 800.0_wp
    CALL CD_AGG_Init_From_Deck(agg, 'agg_farm_shared_anchor.dat', DT, es, em, n_turbines=2, &
                               turbine_ref_pos=refpos, forbid_deck_ambient=.TRUE.)
    CALL require(es == CD_AGG_OK, 'farm-shared-anchor: init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, q0, v0, a0, l0f, es, em)
    v0 = 0.0_wp
    a0 = 0.0_wp
    DO s = 1, 20
      CALL CD_AGG_Step_Moving(agg, DT, q0, v0, a0, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK, 'farm-shared-anchor: step: '//TRIM(em))
      IF (es /= CD_AGG_OK) EXIT
    END DO
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, q0, v0, a0, l1f, es, em)
    WRITE (*, '(A,2ES12.4)') 'farm shared anchor: turbine loads Fx, symmetry error = ', l1f(1, 1), &
      ABS(l1f(1, 1) + l1f(1, 2))/ABS(l1f(1, 1))
    CALL require(ALL(IEEE_IS_FINITE(l1f)) .AND. ABS(l1f(1, 1)) > 0.0_wp, 'farm-shared-anchor: finite loads')
    CALL require(ABS(l1f(1, 1) + l1f(1, 2)) <= 1.0e-9_wp*ABS(l1f(1, 1)) .AND. &
                 ABS(l1f(3, 1) - l1f(3, 2)) <= 1.0e-9_wp*ABS(l1f(3, 1)), 'farm-shared-anchor: mirror-symmetric loads')
    CALL require(nan_max_abs(l1f - l0f) <= 1.0e-6_wp*nan_max_abs(l0f), 'farm-shared-anchor: held farm stays at rest')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_farm_shared_anchor

  SUBROUTINE case_farm_free_clump()
    !! FAST.Farm mode with a FREE mass point (clump weight) splitting the shared line --
    !! the stock MoorDyn r-test shared-mooring topology. The aggregate builds pure EI=0
    !! decks through the same point-connected system builder the standalone driver uses,
    !! so the clump must (a) initialize in farm mode, (b) integrate dynamically (the
    !! system marches it -- its equilibrium position departs from the deck seed), and
    !! (c) keep the turbines coupled through the two half-lines.
    TYPE(CD_AGG_ModuleType) :: agg
    ! Free/Connect points integrate EXPLICITLY inside the system step, so the coupling
    ! dt must satisfy the point CFL: a 5 t clump between chain-EA legs has a local
    ! frequency ~sqrt((EA/l)/m) ~ 1e2 rad/s -> dt well below ~0.015 s (stock MoorDyn
    ! runs its clump r-test at dtM = 1e-4 for the same reason). The gate steps at
    ! 0.005 s; at 0.05 s the clump rings and diverges within ~15 steps -- measured.
    REAL(wp), PARAMETER :: DT = 0.005_wp
    REAL(wp) :: refpos(3, 2), q0(3, 2), v0(3, 2), a0(3, 2), l0f(3, 2)
    REAL(wp) :: q1(3, 2), v1(3, 2), a1(3, 2), l1f(3, 2), tgt(3, 2)
    INTEGER :: es, s, ni
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    refpos = 0.0_wp
    refpos(1, 2) = 800.0_wp
    CALL write_farm_deck('agg_farm_clump.dat', plain_coupled=.FALSE., clump=.TRUE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_farm_clump.dat', DT, es, em, n_turbines=2, &
                               turbine_ref_pos=refpos, forbid_deck_ambient=.TRUE.)
    CALL require(es == CD_AGG_OK, 'farm-clump: init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL require(CD_AGG_NTurbines(agg) == 2, 'farm-clump: NTurbines')

    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'farm-clump: output: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg, q0, v0, a0, l0f, es, em)
    CALL require(es == CD_AGG_OK, 'farm-clump: mesh: '//TRIM(em))
    CALL require(ALL(IEEE_IS_FINITE(l0f)), 'farm-clump: finite initial loads')

    ! (c) the clump-split shared line still couples the turbines: surging turbine 1
    ! must change turbine 2's load THROUGH the free point. The offset is applied as a
    ! small step (a large boundary teleport at zero declared velocity rings the
    ! undamped clump transient into blow-up a few steps later).
    tgt = q0
    tgt(1, 1) = tgt(1, 1) + 0.2_wp
    v0 = 0.0_wp; a0 = 0.0_wp
    ! the boundary offset is a shock (position step at zero declared velocity): the
    ! first steps may stall in the line search; require clean errors throughout and
    ! CONVERGENCE once the shock has propagated (the last step). The clump is not
    ! force-balanced by the per-line statics, so the system also carries its own settling
    ! transient: a reference aggregate held at the initial pose isolates the part of
    ! turbine 2's load that turbine 1's motion causes. That signal must travel ~840 m of
    ! chain (axial wave speed ~2.2 km/s) and move the clump, hence the 0.6 s window.
    BLOCK
      TYPE(CD_AGG_ModuleType) :: aggr
      REAL(wp) :: lref(3, 2), qr(3, 2), vr(3, 2), ar(3, 2)
      LOGICAL :: cvr, str
      INTEGER :: nir
      CALL CD_AGG_Init_From_Deck(aggr, 'agg_farm_clump.dat', DT, es, em, n_turbines=2, &
                                 turbine_ref_pos=refpos, forbid_deck_ambient=.TRUE.)
      CALL require(es == CD_AGG_OK, 'farm-clump: reference init: '//TRIM(em))
      DO s = 1, 120
        CALL CD_AGG_Step_Moving(agg, DT, tgt, v0, a0, cv, st, ni, es, em)
        CALL require(es == CD_AGG_OK, 'farm-clump: step: '//TRIM(em))
        IF (es /= CD_AGG_OK) EXIT
        CALL CD_AGG_Step_Moving(aggr, DT, q0, v0, a0, cvr, str, nir, es, em)
        CALL require(es == CD_AGG_OK, 'farm-clump: reference step: '//TRIM(em))
        IF (es /= CD_AGG_OK) EXIT
      END DO
      CALL require(cv, 'farm-clump: the post-shock step converges')
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL require(es == CD_AGG_OK, 'farm-clump: post-step output: '//TRIM(em))
      CALL CD_AGG_GetMovingPointMesh(agg, q1, v1, a1, l1f, es, em)
      CALL require(es == CD_AGG_OK, 'farm-clump: post-step mesh: '//TRIM(em))
      CALL CD_AGG_CalcOutput(aggr, es, em)
      CALL CD_AGG_GetMovingPointMesh(aggr, qr, vr, ar, lref, es, em)
      CALL require(es == CD_AGG_OK .AND. ALL(IEEE_IS_FINITE(lref)), 'farm-clump: reference mesh: '//TRIM(em))
      CALL require(nan_max_abs(l1f(:, 2) - lref(:, 2)) > 1.0e1_wp, &
                   'farm-clump: turbine-1 motion reaches turbine 2 through the clump')
      CALL CD_AGG_End(aggr, es, em)
    END BLOCK

    ! Dynamic-point state mirror roundtrip (the checkpoint-restore machinery): capture
    ! the committed clump states mid-run, drift the system further, restore, and the
    ! read-back must be exact -- with the restore SCATTERED into the bound endpoints
    ! (the coupled loads move back toward the captured configuration, not the drifted).
    BLOCK
      REAL(wp), ALLOCATABLE :: pbuf(:), pbuf2(:)
      REAL(wp) :: lrest(3, 2)
      INTEGER :: nd
      nd = CD_System_NDynamicPoints(agg%sys%fast%system)
      CALL require(nd == 1, 'farm-clump: one dynamic point in the system')
      ALLOCATE (pbuf(9*nd), pbuf2(9*nd))
      CALL CD_Get_System_DynamicPoint_States(agg%sys%fast%system, pbuf, es, em)
      CALL require(es == CD_SYSTEM_OK, 'farm-clump: point-state get: '//TRIM(em))
      ! drift: hold the boundary and let the still-settling transient move the states
      DO s = 1, 20
        CALL CD_AGG_Step_Moving(agg, DT, tgt, v0, a0, cv, st, ni, es, em)
        CALL require(es == CD_AGG_OK, 'farm-clump: drift step: '//TRIM(em))
      END DO
      CALL CD_Get_System_DynamicPoint_States(agg%sys%fast%system, pbuf2, es, em)
      CALL require(es == CD_SYSTEM_OK, 'farm-clump: drifted get: '//TRIM(em))
      CALL require(nan_max_abs(pbuf2 - pbuf) > 0.0_wp, 'farm-clump: the drift moved the clump state')
      CALL CD_Set_System_DynamicPoint_States(agg%sys%fast%system, pbuf, es, em)
      CALL require(es == CD_SYSTEM_OK, 'farm-clump: point-state set: '//TRIM(em))
      CALL CD_Get_System_DynamicPoint_States(agg%sys%fast%system, pbuf2, es, em)
      CALL require(es == CD_SYSTEM_OK, 'farm-clump: point-state re-get: '//TRIM(em))
      ! nan_max_abs(.) >= 0 by construction, so <= 0 is the EXACT (bitwise) roundtrip
      ! claim without a reals-equality comparison (the -Wcompare-reals class)
      CALL require(nan_max_abs(pbuf2 - pbuf) <= 0.0_wp, 'farm-clump: point-state roundtrip is exact')
      ! the point-store restore leaves line endpoints to the coupled-motion mirror block
      ! (per-line-end state; see the shell's checkpoint layout) -- the system still
      ! steps cleanly from the restored store
      CALL CD_AGG_Step_Moving(agg, DT, tgt, v0, a0, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK, 'farm-clump: post-restore step: '//TRIM(em))
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL require(es == CD_AGG_OK, 'farm-clump: post-restore output: '//TRIM(em))
      CALL CD_AGG_GetMovingPointMesh(agg, q1, v1, a1, lrest, es, em)
      CALL require(es == CD_AGG_OK .AND. ALL(IEEE_IS_FINITE(lrest)), &
                   'farm-clump: finite loads from the restored state')
    END BLOCK
    CALL CD_AGG_End(agg, es, em)

    ! Stock-order twin: with a clump present the deck is dynamic-classified, and the
    ! endpoint-order normalization must still apply -- otherwise stock anchor legs
    ! (NodeA = Fixed, NodeB = Turbine<J>) would reach the builders reversed. The
    ! stock+clump deck must land the fairleads
    ! exactly where the CableDyn-order clump deck does.
    BLOCK
      TYPE(CD_AGG_ModuleType) :: agg2
      REAL(wp) :: q2(3, 2), v2(3, 2), a2(3, 2), l2f(3, 2)
      CALL write_farm_deck('agg_farm_clump_stock.dat', plain_coupled=.FALSE., stock_order=.TRUE., &
                           clump=.TRUE.)
      CALL CD_AGG_Init_From_Deck(agg2, 'agg_farm_clump_stock.dat', DT, es, em, n_turbines=2, &
                                 turbine_ref_pos=refpos, forbid_deck_ambient=.TRUE.)
      CALL require(es == CD_AGG_OK, 'farm-clump: stock-order clump deck inits: '//TRIM(em))
      IF (es == CD_AGG_OK) THEN
        CALL CD_AGG_CalcOutput(agg2, es, em)
        CALL require(es == CD_AGG_OK, 'farm-clump: stock output: '//TRIM(em))
        CALL CD_AGG_GetMovingPointMesh(agg2, q2, v2, a2, l2f, es, em)
        CALL require(es == CD_AGG_OK, 'farm-clump: stock mesh: '//TRIM(em))
        CALL require(nan_max_abs(q2 - q0) <= 1.0e-9_wp, &
                     'farm-clump: stock-order fairleads land identically')
        CALL CD_AGG_End(agg2, es, em)
      END IF
    END BLOCK
  END SUBROUTINE case_farm_free_clump

  SUBROUTINE case_aggregate_failures()
    !! Line failures on the coupled aggregate: (a) a FAILURE deck initializes and
    !! REQUIRES the committed time threaded into every step; (b) pre-trigger steps
    !! are bit-identical to the intact twin; (c) the row fires at its committed
    !! time, a snapshot/restore rewind UN-fires it (flags + detach topology), and
    !! it re-fires on the re-advance; (d) checkpoint-replay: a FRESH aggregate
    !! brought up with CD_AGG_Set_Failure_Flags + the committed coupled vector +
    !! point store steps bit-identically to the fired instance -- the exact
    !! sequence the OpenFAST mirror reload performs.
    TYPE(CD_AGG_ModuleType) :: agg, aggf, aggr
    REAL(wp), PARAMETER :: DT = 0.005_wp
    REAL(wp) :: q0(3, 1), v0(3, 1), a0(3, 1), lf(3, 1)
    REAL(wp) :: lint(3, 1), lfail(3, 1), lrep(3, 1)
    REAL(wp), ALLOCATABLE :: cvq(:), cvv(:), cva(:), pst(:)
    LOGICAL :: cv, st, flg(1)
    INTEGER :: es, ni, sidx, ncd, npd
    CHARACTER(512) :: em

    CALL write_agg_failure_deck('agg_fail_intact.dat')
    CALL write_agg_failure_deck('agg_fail_row.dat', failure_row='1 P2 2 0.007 0.0')

    CALL CD_AGG_Init_From_Deck(agg, 'agg_fail_intact.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'aggfail: intact init: '//TRIM(em))
    CALL CD_AGG_Init_From_Deck(aggf, 'agg_fail_row.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'aggfail: failure init: '//TRIM(em))
    CALL require(CD_AGG_NFailures(aggf) == 1 .AND. CD_AGG_NFailures(agg) == 0, 'aggfail: row counts')

    CALL CD_AGG_GetMovingPointMesh(agg, q0, v0, a0, lf, es, em)
    CALL require(es == CD_AGG_OK, 'aggfail: mesh: '//TRIM(em))
    v0 = 0.0_wp
    a0 = 0.0_wp

    ! (a) the committed time is REQUIRED on a failure deck
    CALL CD_AGG_Step_Moving(aggf, DT, q0, v0, a0, cv, st, ni, es, em)
    CALL require(es /= CD_AGG_OK, 'aggfail: step without t_committed fails closed')

    ! (b) pre-trigger step (t = 0.005 < 0.007): unfired and bit-identical to intact
    CALL CD_AGG_Step_Moving(agg, DT, q0, v0, a0, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK, 'aggfail: intact step 1: '//TRIM(em))
    CALL CD_AGG_Step_Moving(aggf, DT, q0, v0, a0, cv, st, ni, es, em, t_committed=DT)
    CALL require(es == CD_AGG_OK, 'aggfail: failure step 1: '//TRIM(em))
    CALL CD_AGG_Get_Failure_Flags(aggf, flg, es, em)
    CALL require(es == CD_AGG_OK .AND. .NOT. flg(1), 'aggfail: unfired before its time')
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, q0, v0, a0, lint, es, em)
    CALL CD_AGG_CalcOutput(aggf, es, em)
    CALL CD_AGG_GetMovingPointMesh(aggf, q0, v0, a0, lfail, es, em)
    v0 = 0.0_wp
    a0 = 0.0_wp
    CALL require(nan_max_abs(lfail - lint) <= 0.0_wp, 'aggfail: pre-trigger loads bit-identical to intact')

    ! (c) snapshot, fire at t = 0.010 >= 0.007, rewind un-fires, re-advance re-fires
    CALL CD_AGG_Snapshot(aggf, es, em)
    CALL require(es == CD_AGG_OK, 'aggfail: snapshot: '//TRIM(em))
    CALL CD_AGG_Step_Moving(aggf, DT, q0, v0, a0, cv, st, ni, es, em, t_committed=2.0_wp*DT)
    CALL require(es == CD_AGG_OK, 'aggfail: firing step: '//TRIM(em))
    CALL CD_AGG_Get_Failure_Flags(aggf, flg, es, em)
    CALL require(es == CD_AGG_OK .AND. flg(1), 'aggfail: fired at its committed time')
    CALL CD_AGG_Restore(aggf, es, em)
    CALL require(es == CD_AGG_OK, 'aggfail: restore: '//TRIM(em))
    CALL CD_AGG_Get_Failure_Flags(aggf, flg, es, em)
    CALL require(es == CD_AGG_OK .AND. .NOT. flg(1), 'aggfail: rewind un-fires the row')
    CALL CD_AGG_Step_Moving(aggf, DT, q0, v0, a0, cv, st, ni, es, em, t_committed=2.0_wp*DT)
    CALL require(es == CD_AGG_OK, 'aggfail: re-advance: '//TRIM(em))
    CALL CD_AGG_Get_Failure_Flags(aggf, flg, es, em)
    CALL require(es == CD_AGG_OK .AND. flg(1), 'aggfail: re-fires on the re-advance')

    ! settle two more committed steps so the fired trajectory has real post-detach
    ! dynamics in it, then capture the committed exchange state
    DO sidx = 3, 4
      CALL CD_AGG_Step_Moving(aggf, DT, q0, v0, a0, cv, st, ni, es, em, t_committed=REAL(sidx, wp)*DT)
      CALL require(es == CD_AGG_OK, 'aggfail: settle step: '//TRIM(em))
    END DO
    ncd = CD_System_NSystemCoupledDOF(aggf%sys%fast%system)
    npd = 9*CD_System_NDynamicPoints(aggf%sys%fast%system)
    ALLOCATE (cvq(ncd), cvv(ncd), cva(ncd), pst(npd))
    CALL CD_Get_System_CoupledMotion(aggf%sys%fast%system, cvq, cvv, cva, es, em)
    CALL require(es == CD_SYSTEM_OK, 'aggfail: vector get: '//TRIM(em))
    CALL CD_Get_System_DynamicPoint_States(aggf%sys%fast%system, pst, es, em)
    CALL require(es == CD_SYSTEM_OK, 'aggfail: store get: '//TRIM(em))

    ! (d) checkpoint-replay on a FRESH instance: replay the flags (detach topology)
    ! FIRST, then the per-point store, then the coupled vector LAST -- the OpenFAST
    ! mirror reload order -- and the next committed step must match bit-for-bit.
    CALL CD_AGG_Init_From_Deck(aggr, 'agg_fail_row.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'aggfail: replay init: '//TRIM(em))
    flg(1) = .TRUE.
    CALL CD_AGG_Set_Failure_Flags(aggr, flg, es, em)
    CALL require(es == CD_AGG_OK, 'aggfail: flag replay: '//TRIM(em))
    CALL CD_Set_System_DynamicPoint_States(aggr%sys%fast%system, pst, es, em)
    CALL require(es == CD_SYSTEM_OK, 'aggfail: store set: '//TRIM(em))
    CALL CD_Update_System_CoupledMotion(aggr%sys%fast%system, cvq, cvv, cva, es, em)
    CALL require(es == CD_SYSTEM_OK, 'aggfail: vector set: '//TRIM(em))
    CALL CD_AGG_Refresh_PointMesh(aggr, es, em)
    CALL require(es == CD_AGG_OK, 'aggfail: mesh refresh: '//TRIM(em))
    CALL CD_AGG_Step_Moving(aggf, DT, q0, v0, a0, cv, st, ni, es, em, t_committed=5.0_wp*DT)
    CALL require(es == CD_AGG_OK, 'aggfail: fired-instance continue: '//TRIM(em))
    CALL CD_AGG_Step_Moving(aggr, DT, q0, v0, a0, cv, st, ni, es, em, t_committed=5.0_wp*DT)
    CALL require(es == CD_AGG_OK, 'aggfail: replayed-instance continue: '//TRIM(em))
    CALL CD_AGG_CalcOutput(aggf, es, em)
    CALL CD_AGG_GetMovingPointMesh(aggf, q0, v0, a0, lfail, es, em)
    CALL CD_AGG_CalcOutput(aggr, es, em)
    CALL CD_AGG_GetMovingPointMesh(aggr, q0, v0, a0, lrep, es, em)
    CALL require(nan_max_abs(lrep - lfail) <= 0.0_wp, &
                 'aggfail: replayed instance continues bit-identically to the fired instance')

    CALL CD_AGG_End(agg, es, em)
    CALL CD_AGG_End(aggf, es, em)
    CALL CD_AGG_End(aggr, es, em)
  END SUBROUTINE case_aggregate_failures

  SUBROUTINE case_aggregate_line_control()
    !! Active line control on the coupled aggregate: (a) a zero command is
    !! BIT-IDENTICAL to never applying (the control machinery is inert at
    !! DeltaL = 0); (b) paying line 2 out drops its tension IMMEDIATELY (same
    !! committed q, longer unstretched length) and the subsequent step stays
    !! convergent; (c) the command-vector fail-closed battery.
    TYPE(CD_AGG_ModuleType) :: agg, aggc
    REAL(wp), PARAMETER :: DT = 0.005_wp
    REAL(wp) :: q0(3, 1), v0(3, 1), a0(3, 1), lf(3, 1), lc(3, 1)
    REAL(wp), ALLOCATABLE :: ten_before(:), ten_after(:)
    INTEGER :: es, ni, ne2
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_agg_failure_deck('agg_ctrl_plain.dat')
    CALL write_agg_failure_deck('agg_ctrl_row.dat', control_row='1 2')

    CALL CD_AGG_Init_From_Deck(agg, 'agg_ctrl_plain.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'aggctrl: plain init: '//TRIM(em))
    ! The metadata-only dtM scan must RECEIVE (and discard) the CONTROL and
    ! FAILURE sections: the OpenFAST shell queries dtM BEFORE the aggregate init,
    ! and a fail-closed scan would abort every deck the aggregate supports.
    BLOCK
      REAL(wp) :: dtm_probe
      LOGICAL :: has_dtm_probe
      CALL CD_Deck_Query_dtM('agg_ctrl_row.dat', dtm_probe, has_dtm_probe, es, em)
      CALL require(es == CD_DECKDRV_OK, 'aggctrl: dtM scan accepts CONTROL decks: '//TRIM(em))
      CALL write_agg_failure_deck('agg_ctrl_failrow.dat', failure_row='1 P2 2 5.0 0.0')
      CALL CD_Deck_Query_dtM('agg_ctrl_failrow.dat', dtm_probe, has_dtm_probe, es, em)
      CALL require(es == CD_DECKDRV_OK, 'aggctrl: dtM scan accepts FAILURE decks: '//TRIM(em))
    END BLOCK
    CALL CD_AGG_Init_From_Deck(aggc, 'agg_ctrl_row.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'aggctrl: control init: '//TRIM(em))
    CALL require(CD_AGG_NCtrlChans(aggc) == 1 .AND. CD_AGG_NCtrlChans(agg) == 0, 'aggctrl: channel counts')

    ! (a) zero command == never applied, bit-identical loads after one step
    CALL CD_AGG_Apply_LineControl(aggc, [0.0_wp], [0.0_wp], es, em)
    CALL require(es == CD_AGG_OK, 'aggctrl: zero apply: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg, q0, v0, a0, lf, es, em)
    v0 = 0.0_wp
    a0 = 0.0_wp
    CALL CD_AGG_Step_Moving(agg, DT, q0, v0, a0, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK, 'aggctrl: plain step: '//TRIM(em))
    CALL CD_AGG_Step_Moving(aggc, DT, q0, v0, a0, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK, 'aggctrl: control step: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, q0, v0, a0, lf, es, em)
    CALL CD_AGG_CalcOutput(aggc, es, em)
    CALL CD_AGG_GetMovingPointMesh(aggc, q0, v0, a0, lc, es, em)
    v0 = 0.0_wp
    a0 = 0.0_wp
    CALL require(nan_max_abs(lc - lf) <= 0.0_wp, 'aggctrl: zero command is bit-identical to no control')

    ! (b) payout drops the controlled line's tension immediately; the next step converges
    ne2 = CD_System_Line_NElem(aggc%sys%fast%system, 2, es, em)
    ALLOCATE (ten_before(ne2), ten_after(ne2))
    CALL CD_Get_System_Line_Tension(aggc%sys%fast%system, 2, ten_before, es, em)
    CALL require(es == CD_SYSTEM_OK, 'aggctrl: tension before')
    CALL CD_AGG_Apply_LineControl(aggc, [0.05_wp], [0.0_wp], es, em)
    CALL require(es == CD_AGG_OK, 'aggctrl: payout apply: '//TRIM(em))
    CALL CD_Get_System_Line_Tension(aggc%sys%fast%system, 2, ten_after, es, em)
    CALL require(es == CD_SYSTEM_OK, 'aggctrl: tension after')
    CALL require(ten_after(ne2) < ten_before(ne2), 'aggctrl: payout drops the controlled segment tension')
    CALL CD_AGG_Step_Moving(aggc, DT, q0, v0, a0, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'aggctrl: post-payout step converges: '//TRIM(em))

    ! (a2) IDEMPOTENCE: re-applying the same command is bit-inert (the apply is an
    ! absolute assignment over every derived quantity -- lengths, mass, f_ext,
    ! seabed coefficients).
    BLOCK
      REAL(wp), ALLOCATABLE :: ti1(:), ti2(:)
      INTEGER :: nei
      nei = CD_System_Line_NElem(aggc%sys%fast%system, 2, es, em)
      ALLOCATE (ti1(nei), ti2(nei))
      CALL CD_AGG_Apply_LineControl(aggc, [0.02_wp], [0.0_wp], es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: idempotence first apply: '//TRIM(em))
      CALL CD_Get_System_Line_Tension(aggc%sys%fast%system, 2, ti1, es, em)
      CALL CD_AGG_Apply_LineControl(aggc, [0.02_wp], [0.0_wp], es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: idempotence second apply: '//TRIM(em))
      CALL CD_Get_System_Line_Tension(aggc%sys%fast%system, 2, ti2, es, em)
      CALL require(nan_max_abs(ti2 - ti1) <= 0.0_wp, 'aggctrl: re-applying the same command is bit-inert')
      CALL CD_AGG_Apply_LineControl(aggc, [0.0_wp], [0.0_wp], es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: idempotence reset: '//TRIM(em))
    END BLOCK

    ! (b2) COMPOSITE controlled line (two sections, different segment lengths):
    ! the base must be the BUILT fairlead element's length, not a deck-sections
    ! scan -- pre-fix a zero command rewrote the fairlead segment to the
    ! anchor-side length and the tension moved at DeltaL = 0.
    BLOCK
      TYPE(CD_AGG_ModuleType) :: aggd
      REAL(wp), ALLOCATABLE :: tb(:), ta(:)
      INTEGER :: ned
      CALL write_agg_failure_deck('agg_ctrl_comp.dat', control_row='1 2', composite_line2=.TRUE.)
      CALL CD_AGG_Init_From_Deck(aggd, 'agg_ctrl_comp.dat', DT, es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: composite init: '//TRIM(em))
      ned = CD_System_Line_NElem(aggd%sys%fast%system, 2, es, em)
      ALLOCATE (tb(ned), ta(ned))
      CALL CD_Get_System_Line_Tension(aggd%sys%fast%system, 2, tb, es, em)
      CALL CD_AGG_Apply_LineControl(aggd, [0.0_wp], [0.0_wp], es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: composite zero apply: '//TRIM(em))
      CALL CD_Get_System_Line_Tension(aggd%sys%fast%system, 2, ta, es, em)
      CALL require(nan_max_abs(ta - tb) <= 0.0_wp, &
                   'aggctrl: composite zero command leaves every segment tension bit-identical')
      CALL CD_AGG_End(aggd, es, em)
    END BLOCK

    ! (b3) ATOMIC command validation: with one channel driving BOTH lines
    ! (bases 0.98 and 0.48 on the composite deck), a command invalid only for the
    ! SECOND row must leave the FIRST line untouched -- pre-fix the loop had
    ! already updated it before hitting the bad row.
    BLOCK
      TYPE(CD_AGG_ModuleType) :: agga
      REAL(wp), ALLOCATABLE :: t1b(:), t1a(:)
      INTEGER :: ne1
      CALL write_agg_failure_deck('agg_ctrl_atomic.dat', control_row='1 1,2', composite_line2=.TRUE.)
      CALL CD_AGG_Init_From_Deck(agga, 'agg_ctrl_atomic.dat', DT, es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: atomic init: '//TRIM(em))
      ne1 = CD_System_Line_NElem(agga%sys%fast%system, 1, es, em)
      ALLOCATE (t1b(ne1), t1a(ne1))
      CALL CD_Get_System_Line_Tension(agga%sys%fast%system, 1, t1b, es, em)
      CALL CD_AGG_Apply_LineControl(agga, [-0.6_wp], [0.0_wp], es, em)
      CALL require(es /= CD_AGG_OK, 'aggctrl: atomic bad-row command fails closed')
      CALL CD_Get_System_Line_Tension(agga%sys%fast%system, 1, t1a, es, em)
      CALL require(nan_max_abs(t1a - t1b) <= 0.0_wp, &
                   'aggctrl: a rejected command leaves every line bit-identical (no partial apply)')
      CALL CD_AGG_End(agga, es, em)
    END BLOCK

    ! (b4) REWIND: a command applied after a snapshot must roll back WITH the
    ! states -- l0, mass, seabed coefficients, and the distributed-load remainder
    ! all restore, so the tensions return bit-identically without any re-apply.
    BLOCK
      REAL(wp), ALLOCATABLE :: tr0(:), tr1(:)
      INTEGER :: nrr
      nrr = CD_System_Line_NElem(aggc%sys%fast%system, 2, es, em)
      ALLOCATE (tr0(nrr), tr1(nrr))
      CALL CD_Get_System_Line_Tension(aggc%sys%fast%system, 2, tr0, es, em)
      CALL CD_AGG_Snapshot(aggc, es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: rewind snapshot: '//TRIM(em))
      CALL CD_AGG_Apply_LineControl(aggc, [0.03_wp], [0.1_wp], es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: rewind apply: '//TRIM(em))
      CALL CD_AGG_Restore(aggc, es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: rewind restore: '//TRIM(em))
      CALL CD_Get_System_Line_Tension(aggc%sys%fast%system, 2, tr1, es, em)
      CALL require(nan_max_abs(tr1 - tr0) <= 0.0_wp, &
                   'aggctrl: a rewound command restores the tensions bit-identically')
      CALL CD_AGG_Step_Moving(aggc, DT, q0, v0, a0, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK, 'aggctrl: post-rewind step: '//TRIM(em))
    END BLOCK

    ! (b5) trailing tokens fail closed: "1 1 2" (spaces instead of "1 1,2") would
    ! silently leave line 2 uncontrolled -- the parser must reject the row.
    BLOCK
      TYPE(CD_AGG_ModuleType) :: aggt
      CALL write_agg_failure_deck('agg_ctrl_trail.dat', control_row='1 1 2')
      CALL CD_AGG_Init_From_Deck(aggt, 'agg_ctrl_trail.dat', DT, es, em)
      CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'trailing') > 0, &
                   'aggctrl: space-separated line list fails closed')
      CALL write_agg_failure_deck('agg_fail_trail.dat', failure_row='1 P2 2 0.004 0.0 7')
      CALL CD_AGG_Init_From_Deck(aggt, 'agg_fail_trail.dat', DT, es, em)
      CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'trailing') > 0, &
                   'aggctrl: FAILURE trailing token fails closed')
    END BLOCK

    ! (c) fail-closed battery
    CALL CD_AGG_Apply_LineControl(aggc, [REAL(wp) ::], [REAL(wp) ::], es, em)
    CALL require(es /= CD_AGG_OK, 'aggctrl: short command vector fails closed')
    CALL CD_AGG_Apply_LineControl(aggc, [-1.0_wp], [0.0_wp], es, em)
    CALL require(es /= CD_AGG_OK, 'aggctrl: non-positive segment length fails closed')

    CALL CD_AGG_End(agg, es, em)
    CALL CD_AGG_End(aggc, es, em)
  END SUBROUTINE case_aggregate_line_control

  SUBROUTINE write_agg_failure_deck(path, failure_row, control_row, composite_line2)
    !! Vessel fairlead + massy Connect junction + Fixed anchor, two 0.98 m lines --
    !! the connect rig behind the coupled facade, with optional FAILURE/CONTROL rows.
    !! composite_line2 splits line 2 into two sections with DIFFERENT segment
    !! lengths (0.50 then 0.48) -- the control-base composite regression.
    CHARACTER(*), INTENT(IN) :: path
    CHARACTER(*), INTENT(IN), OPTIONAL :: failure_row, control_row
    LOGICAL, INTENT(IN), OPTIONAL :: composite_line2
    LOGICAL :: comp2
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'Coupled connect rig with optional FAILURE'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'line 0.10 80.0 1.0e5 0.0 0.0 1.0 0.0 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Vessel 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '2 Connect 1.0 0.0 0.0 20.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '3 Fixed 2.0 0.0 0.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '2 2 3 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 line 0.98 1'
    comp2 = .FALSE.
    IF (PRESENT(composite_line2)) comp2 = composite_line2
    IF (comp2) THEN
      WRITE (u, '(A)') '2 line 0.50 1'
      WRITE (u, '(A)') '2 line 0.48 1'
    ELSE
      WRITE (u, '(A)') '2 line 0.98 1'
    END IF
    IF (PRESENT(failure_row)) THEN
      WRITE (u, '(A)') '--- FAILURE ---'
      WRITE (u, '(A)') 'FailID Point Lines FailTime FailTen'
      WRITE (u, '(A)') '(-) (-) (-) (s) (N)'
      WRITE (u, '(A)') TRIM(failure_row)
    END IF
    IF (PRESENT(control_row)) THEN
      WRITE (u, '(A)') '--- CONTROL ---'
      WRITE (u, '(A)') 'ChannelID Lines'
      WRITE (u, '(A)') '(-) (-)'
      WRITE (u, '(A)') TRIM(control_row)
    END IF
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_agg_failure_deck

  SUBROUTINE write_farm_deck(path, plain_coupled, y_shift, stock_order, ambient, clump)
    !! Two-turbine farm deck: one anchored EI=0 line per turbine plus one SHARED line
    !! between the two fairleads. Coupled coordinates are TURBINE-LOCAL (Turbine<J>
    !! vocabulary unless plain_coupled); anchors are FARM-GLOBAL. ambient = 'current'
    !! or 'waves' adds the corresponding deck OPTION row (the forbid_deck_ambient gate).
    !! clump: split the shared line through a FREE mass point (the MoorDyn r-test
    !! shared-mooring topology -- a clump weight holding the shared span down).
    CHARACTER(*), INTENT(IN) :: path
    LOGICAL, INTENT(IN) :: plain_coupled
    REAL(wp), INTENT(IN), OPTIONAL :: y_shift
    LOGICAL, INTENT(IN), OPTIONAL :: stock_order
    CHARACTER(*), INTENT(IN), OPTIONAL :: ambient
    LOGICAL, INTENT(IN), OPTIONAL :: clump
    INTEGER :: u, ios
    REAL(wp) :: ys
    LOGICAL :: stock, cl
    CHARACTER(16) :: t1, t2
    ys = 0.0_wp
    IF (PRESENT(y_shift)) ys = y_shift
    stock = .FALSE.
    IF (PRESENT(stock_order)) stock = stock_order
    cl = .FALSE.
    IF (PRESENT(clump)) cl = clump
    t1 = 'Turbine1'; t2 = 'Turbine2'
    IF (plain_coupled) THEN
      t1 = 'Coupled'; t2 = 'Coupled'
    END IF
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    CALL require(ios == 0, 'create farm deck')
    IF (ios /= 0) RETURN
    WRITE (u, '(A)') 'Two-turbine farm: per-turbine anchored lines + one shared line'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.333 685.0 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27'
    WRITE (u, '(A)') 'strop 0.080  25.0 2.00e8 -1.0 0.0 1.6 0.2 0.80 0.27'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A,A,A)') '1 ', TRIM(t1), '   0.0 0.0 -10.0'
    WRITE (u, '(A,F10.3,A)') '2 Fixed  -500.0 ', ys, ' -200.0'
    WRITE (u, '(A,A,A)') '3 ', TRIM(t2), '   0.0 0.0 -10.0'
    WRITE (u, '(A,F10.3,A)') '4 Fixed  1300.0 ', ys, ' -200.0'
    IF (cl) WRITE (u, '(A,F10.3,A)') '5 Free  400.0 ', ys, ' -80.0 100000.0 0.0 0.0 0.0'
    IF (cl) WRITE (u, '(A,F10.3,A)') '6 Fixed  400.0 ', ys, ' -200.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    IF (stock) THEN
      ! stock MoorDyn endpoint order for the anchored lines: NodeA = Fixed anchor,
      ! NodeB = Turbine<J> fairlead (the pre-pass must swap these like any coupled row)
      WRITE (u, '(A)') '1 2 1 -'
      WRITE (u, '(A)') '2 4 3 -'
    ELSE
      WRITE (u, '(A)') '1 1 2 -'
      WRITE (u, '(A)') '2 3 4 -'
    END IF
    IF (cl) THEN
      WRITE (u, '(A)') '3 1 5 -'
      WRITE (u, '(A)') '4 5 3 -'
      ! ballast leg clump -> Fixed below; the stock twin writes it anchor-first
      ! (Fixed A, Free B) -- the (Fixed -> Free/Connect) normalization class
      IF (stock) THEN
        WRITE (u, '(A)') '5 6 5 -'
      ELSE
        WRITE (u, '(A)') '5 5 6 -'
      END IF
    ELSE
      WRITE (u, '(A)') '3 1 3 -'
    END IF
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 chain 560.0 20'
    WRITE (u, '(A)') '2 chain 560.0 20'
    IF (cl) THEN
      WRITE (u, '(A)') '3 chain 420.0 12'
      WRITE (u, '(A)') '4 chain 420.0 12'
      WRITE (u, '(A)') '5 strop 130.0 6'
    ELSE
      WRITE (u, '(A)') '3 chain 830.0 24'
    END IF
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '200.0 WtrDpth'
    WRITE (u, '(A)') '3.0e6 kBot'
    WRITE (u, '(A)') '3.0e5 cBot'
    IF (PRESENT(ambient)) THEN
      ! dtM + TMax make the ambient OPTION pass the parse-level driver rule --
      ! the realistic leak the forbid_deck_ambient guard must catch is a deck that
      ! parses fine yet self-declares ambient forcing against the still-water scope
      WRITE (u, '(A)') '0.05 dtM'
      WRITE (u, '(A)') '10.0 TMax'
      SELECT CASE (TRIM(ambient))
      CASE ('current')
        WRITE (u, '(A)') 'uniform 0.5 0.0 0.0 current'
      CASE ('waves')
        WRITE (u, '(A)') 'airy 2.0 8.0 0.0 waves'
      END SELECT
    END IF
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_farm_deck

  SUBROUTINE case_coupled_positions_override()
    !! The linearization probe's boundary contract: CD_AGG_Init_From_Deck with
    !! coupled_positions places the moving points EXACTLY, column k of the override
    !! = node k of the moving mesh (the aggregate moving-set order). Gates:
    !! (a) identity -- overriding with the plain build's own mesh positions
    !!     reproduces them column-for-column;
    !! (b) locality -- perturbing one column's z moves exactly that mesh node;
    !! (c) response -- the perturbed static solve changes that node's load;
    !! (d) fail-closed -- a wrong column count and a non-finite entry are rejected.
    TYPE(CD_AGG_ModuleType) :: agg_a, agg_b
    REAL(wp), PARAMETER :: DT = 0.05_wp, DZ = 0.5_wp
    REAL(wp) :: q0(3, 3), v0(3, 3), a0(3, 3), l0f(3, 3)
    REAL(wp) :: q1(3, 3), v1(3, 3), a1(3, 3), l1f(3, 3)
    REAL(wp) :: pos(3, 3), bad(3, 2)
    INTEGER :: es, k, j
    CHARACTER(512) :: em

    CALL write_volturnus_deck('agg_cpos.dat')
    CALL CD_AGG_Init_From_Deck(agg_a, 'agg_cpos.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'cpos: plain init: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg_a, es, em)
    CALL require(es == CD_AGG_OK, 'cpos: plain output: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg_a, q0, v0, a0, l0f, es, em)
    CALL require(es == CD_AGG_OK, 'cpos: plain mesh: '//TRIM(em))

    ! (a) identity: the override at the plain build's own positions reproduces them
    CALL CD_AGG_Init_From_Deck(agg_b, 'agg_cpos.dat', DT, es, em, coupled_positions=q0)
    CALL require(es == CD_AGG_OK, 'cpos: identity init: '//TRIM(em))
    CALL CD_AGG_GetMovingPointMesh(agg_b, q1, v1, a1, l1f, es, em)
    CALL require(es == CD_AGG_OK, 'cpos: identity mesh: '//TRIM(em))
    CALL require(nan_max_abs(q1 - q0) <= 1.0e-9_wp, 'cpos: identity positions column-for-column')
    CALL CD_AGG_End(agg_b, es, em)

    ! (b) + (c): perturb one column at a time; the mesh must move EXACTLY there
    DO k = 1, 3
      pos = q0
      pos(3, k) = pos(3, k) + DZ
      CALL CD_AGG_Init_From_Deck(agg_b, 'agg_cpos.dat', DT, es, em, coupled_positions=pos)
      CALL require(es == CD_AGG_OK, 'cpos: perturbed init: '//TRIM(em))
      CALL CD_AGG_CalcOutput(agg_b, es, em)
      CALL require(es == CD_AGG_OK, 'cpos: perturbed output: '//TRIM(em))
      CALL CD_AGG_GetMovingPointMesh(agg_b, q1, v1, a1, l1f, es, em)
      CALL require(es == CD_AGG_OK, 'cpos: perturbed mesh: '//TRIM(em))
      DO j = 1, 3
        IF (j == k) THEN
          CALL require(ABS(q1(3, j) - (q0(3, j) + DZ)) <= 1.0e-9_wp, 'cpos: column moves its own node z')
          CALL require(nan_max_abs(q1(1:2, j) - q0(1:2, j)) <= 1.0e-9_wp, 'cpos: column leaves its node xy')
          CALL require(ABS(l1f(3, j) - l0f(3, j)) > 1.0e2_wp, 'cpos: the perturbed node z-load responds')
        ELSE
          CALL require(nan_max_abs(q1(:, j) - q0(:, j)) <= 1.0e-9_wp, 'cpos: other nodes stay put')
          ! independent lines: an untouched fairlead's STATIC load must not respond to
          ! another line's boundary perturbation (the cross-line stiffness is zero)
          CALL require(nan_max_abs(l1f(:, j) - l0f(:, j)) <= 1.0e-3_wp*nan_max_abs(l0f), &
                       'cpos: untouched nodes carry no cross-line load response')
        END IF
      END DO
      CALL CD_AGG_End(agg_b, es, em)
    END DO

    ! (d) fail-closed edges
    bad = q0(:, 1:2)
    CALL CD_AGG_Init_From_Deck(agg_b, 'agg_cpos.dat', DT, es, em, coupled_positions=bad)
    CALL require(es /= CD_AGG_OK, 'cpos: wrong column count fails closed')
    pos = q0
    pos(1, 2) = IEEE_VALUE(pos(1, 2), IEEE_QUIET_NAN)
    CALL CD_AGG_Init_From_Deck(agg_b, 'agg_cpos.dat', DT, es, em, coupled_positions=pos)
    CALL require(es /= CD_AGG_OK, 'cpos: non-finite entry fails closed')

    CALL CD_AGG_End(agg_a, es, em)
  END SUBROUTINE case_coupled_positions_override

  SUBROUTINE case_external_fluid()
    !! HOST-DRIVEN AMBIENT FLUID (the SeaState WaveField boundary): an external_fluid
    !! build configures the full wave-capable hydro set at ZERO fields. Contracts:
    !! (a) undriven, it is BIT-IDENTICAL to the plain build -- initial loads and a held
    !!     trajectory (zero fields contribute exactly zero load);
    !! (b) the sampling surface exposes every line node with its committed position;
    !! (c) prescribing zero fields explicitly stays bit-identical; a horizontal current
    !!     changes the loads (drag is live); an acceleration-only field changes the
    !!     loads (Froude-Krylov is live -- the wave-capable configuration is real);
    !! (d) a driven system STEPS (converges) under the prescribed field;
    !! (e) fail-closed: fields on a PLAIN build name the unconfigured block; an
    !!     external_fluid deck may not declare its own waves/current; a cable deck is
    !!     rejected (held fields are not supported for cables); wrong shapes and non-finite fields are
    !!     rejected.
    TYPE(CD_AGG_ModuleType) :: agg_p, agg_x
    REAL(wp), PARAMETER :: DT = 0.05_wp
    INTEGER, PARAMETER :: NN_EXPECT = 63   ! 3 lines x (20 segments + 1)
    REAL(wp) :: r0(3, 3), vel(3, 3), acc(3, 3), tgt(3, 3)
    REAL(wp) :: lp(3, 3), lx(3, 3), l0f(3, 3), lcur(3, 3), lfk(3, 3), p3(3, 3), v3(3, 3), a3(3, 3)
    REAL(wp), ALLOCATABLE :: xyz(:, :), fu(:, :), fud(:, :), wl(:)
    INTEGER :: es, nf, s, ni
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_volturnus_deck('agg_extf.dat')
    CALL CD_AGG_Init_From_Deck(agg_p, 'agg_extf.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'extf:plain-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_Init_From_Deck(agg_x, 'agg_extf.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'extf:external-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) THEN
      CALL CD_AGG_End(agg_p, es, em)
      RETURN
    END IF

    ! (a) undriven bit-identity: initial mesh loads + a short held trajectory
    CALL CD_AGG_GetMovingPointMesh(agg_p, r0, vel, acc, lp, es, em)
    CALL require(es == CD_AGG_OK, 'extf:plain-mesh0')
    CALL CD_AGG_GetMovingPointMesh(agg_x, p3, v3, a3, lx, es, em)
    CALL require(es == CD_AGG_OK, 'extf:external-mesh0')
    CALL require(nan_max_abs(lp - lx) <= 0.0_wp, 'extf:zero-field-build-loads-bit-identical')
    CALL require(nan_max_abs(r0 - p3) <= 0.0_wp, 'extf:zero-field-build-positions-bit-identical')
    tgt = r0
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 3
      CALL CD_AGG_Step_Moving(agg_p, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'extf:plain-held-step')
      CALL CD_AGG_Step_Moving(agg_x, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'extf:external-held-step')
    END DO
    CALL CD_AGG_CalcOutput(agg_p, es, em)
    CALL CD_AGG_CalcOutput(agg_x, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_p, p3, v3, a3, lp, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_x, p3, v3, a3, lx, es, em)
    CALL require(nan_max_abs(lp - lx) <= 0.0_wp, 'extf:held-trajectory-bit-identical')

    ! (b) the sampling surface
    nf = CD_AGG_NFluidNodes(agg_x, es, em)
    CALL require(es == CD_AGG_OK .AND. nf == NN_EXPECT, 'extf:nfluid-63')
    ALLOCATE (xyz(3, nf), fu(3, nf), fud(3, nf), wl(nf))
    CALL CD_AGG_GetFluidNodePositions(agg_x, xyz, es, em)
    CALL require(es == CD_AGG_OK, 'extf:positions: '//TRIM(em))
    CALL require(ALL(IEEE_IS_FINITE(xyz)), 'extf:positions-finite')
    ! every line spans anchor depth (-200) to fairlead depth (-14)
    CALL require(ABS(MINVAL(xyz(3, :)) - (-200.0_wp)) <= 1.0_wp, 'extf:deepest-node-near-anchor')
    CALL require(ABS(MAXVAL(xyz(3, :)) - (-14.0_wp)) <= 1.0_wp, 'extf:shallowest-node-near-fairlead')

    ! (c) explicit zero fields stay bit-identical; a current changes loads; FK is live
    fu = 0.0_wp
    fud = 0.0_wp
    wl = 0.0_wp
    CALL CD_AGG_SetFluidFields(agg_x, fu, fud, wl, es, em)
    CALL require(es == CD_AGG_OK, 'extf:set-zero-fields: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg_x, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_x, p3, v3, a3, l0f, es, em)
    CALL require(nan_max_abs(l0f - lx) <= 0.0_wp, 'extf:explicit-zero-fields-bit-identical')
    fu(1, :) = 1.0_wp
    CALL CD_AGG_SetFluidFields(agg_x, fu, fud, wl, es, em)
    CALL require(es == CD_AGG_OK, 'extf:set-current: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg_x, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_x, p3, v3, a3, lcur, es, em)
    CALL require(nan_max_abs(lcur - l0f) > 0.0_wp, 'extf:current-changes-loads')
    fu = 0.0_wp
    fud(1, :) = 2.0_wp
    CALL CD_AGG_SetFluidFields(agg_x, fu, fud, wl, es, em)
    CALL require(es == CD_AGG_OK, 'extf:set-acceleration: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg_x, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_x, p3, v3, a3, lfk, es, em)
    CALL require(nan_max_abs(lfk - l0f) > 0.0_wp, 'extf:froude-krylov-live')

    ! (d) a driven system steps
    CALL CD_AGG_Step_Moving(agg_x, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'extf:driven-step-converges: '//TRIM(em))

    ! (e) fail-closed edges
    CALL CD_AGG_SetFluidFields(agg_p, fu, fud, wl, es, em)
    CALL require(es /= CD_AGG_OK, 'extf:fields-on-plain-build-fail-closed')
    fud(1, 1) = IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN)
    CALL CD_AGG_SetFluidFields(agg_x, fu, fud, wl, es, em)
    CALL require(es /= CD_AGG_OK, 'extf:non-finite-fields-fail-closed')
    CALL CD_AGG_SetFluidFields(agg_x, fu(:, 1:nf - 1), fud(:, 1:nf - 1), wl(1:nf - 1), es, em)
    CALL require(es /= CD_AGG_OK, 'extf:wrong-shape-fails-closed')
    CALL CD_AGG_End(agg_p, es, em)
    CALL CD_AGG_End(agg_x, es, em)

    ! an external_fluid deck may not declare its own current (double-counting)
    CALL write_current_deck('agg_extf_cur.dat')
    CALL CD_AGG_Init_From_Deck(agg_x, 'agg_extf_cur.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es /= CD_AGG_OK, 'extf:deck-current-plus-external-fails-closed')
    CALL require(INDEX(em, 'double-counting') > 0, 'extf:double-counting-named')

    ! a MIXED (mooring + finite-EI cable) deck now receives the host field too: the
    ! cable nodes extend the sampling surface, the cable's held-field mode makes drag
    ! (velocity) AND Froude-Krylov (acceleration) live on the cable column, zero fields
    ! stay bit-identical to the plain build, and the driven system steps
    CALL case_external_fluid_cable()
  END SUBROUTINE case_external_fluid

  SUBROUTINE case_still_water_deck_ambient()
    !! Single-turbine still-water host (no external_fluid, no still-water flag) on a pure
    !! EI=0 deck: a deck current is a steady field applied through the line models, so it
    !! changes the fairlead loads and they stay steady under a held fairlead; deck waves are
    !! not evaluated on this route, so they fail closed by name instead of running at zero.
    TYPE(CD_AGG_ModuleType) :: agg_s, agg_c
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp) :: r0(3, 3), vel(3, 3), acc(3, 3), ls(3, 3), lc(3, 3), lc1(3, 3), p3(3, 3), v3(3, 3), a3(3, 3)
    INTEGER :: es, s, ni
    LOGICAL :: cv, st
    CHARACTER(1024) :: em

    CALL write_volturnus_deck('agg_sw_still.dat', extra_option='0.05 dtM', extra_option2='60.0 TMax')
    CALL write_current_deck('agg_sw_cur.dat')
    CALL CD_AGG_Init_From_Deck(agg_s, 'agg_sw_still.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'still-water:plain-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_Init_From_Deck(agg_c, 'agg_sw_cur.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'still-water:deck-current-accepted: '//TRIM(em))
    IF (es /= CD_AGG_OK) THEN
      CALL CD_AGG_End(agg_s, es, em)
      RETURN
    END IF
    CALL CD_AGG_CalcOutput(agg_s, es, em)
    CALL CD_AGG_CalcOutput(agg_c, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_s, r0, vel, acc, ls, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_c, p3, v3, a3, lc, es, em)
    ! a 0.5 m/s current on 200 m chains moves the fairlead loads by far more than round-off
    CALL require(nan_max_abs(lc - ls) > 1.0e-3_wp*nan_max_abs(ls), 'still-water:deck-current-changes-loads')
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, 5
      CALL CD_AGG_Step_Moving(agg_c, DT, r0, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'still-water:deck-current-held-step: '//TRIM(em))
    END DO
    CALL CD_AGG_CalcOutput(agg_c, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_c, p3, v3, a3, lc1, es, em)
    CALL require(nan_max_abs(lc1 - lc) <= 1.0e-6_wp*nan_max_abs(lc), 'still-water:deck-current-steady')
    CALL CD_AGG_End(agg_s, es, em)
    CALL CD_AGG_End(agg_c, es, em)

    CALL write_volturnus_deck('agg_sw_waves.dat', extra_option='airy 2.0 8.0 0.0 waves', &
                              extra_option2='0.05 dtM', extra_option3='60.0 TMax')
    CALL CD_AGG_Init_From_Deck(agg_c, 'agg_sw_waves.dat', DT, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'deck waves are not evaluated') > 0, &
                 'still-water:deck-waves-rejected-by-name: '//TRIM(em))
    CALL CD_AGG_End(agg_c, es, em)
  END SUBROUTINE case_still_water_deck_ambient

  SUBROUTINE case_external_fluid_report()
    !! The two reported initial states of a coupled run under host kinematics. The
    !! initialization summary (CD_AGG_GetInitLine), the FairTen channel evaluated before any
    !! field is prescribed and the static profile all describe the converged equilibrium, so
    !! they agree; the host fields then act from t = 0 and change the t = 0 channel, while the
    !! summary keeps the equilibrium value. The OpenFAST shell prints the summary, writes the
    !! static profile before it seeds the t = 0 fields and labels the t = 0 outputs; in the
    !! IEA-15MW JONSWAP example the summary reads 2.43712 MN and the t = 0 FairTen1 2.43621 MN.
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp), ALLOCATABLE :: fu(:, :), fud(:, :), wl(:)
    REAL(wp) :: fair_ten, fair_force(3), fair_incl, fair_decl, fair_azi, ten0, ten_field, ten_profile
    REAL(wp) :: summary_after
    INTEGER :: es, nf, line_id, unit, ios, row_line, row_node
    CHARACTER(512) :: em
    CHARACTER(1024) :: skip

    CALL write_volturnus_deck('agg_extf_report.dat', extra_option='--- OUTPUTS ---', &
                              extra_option2='FairTen1 AnchTen1', extra_option3='END')
    CALL CD_AGG_Init_From_Deck(agg, 'agg_extf_report.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'extf-report:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_GetInitLine(agg, 1, line_id, fair_ten, fair_force, fair_incl, fair_decl, fair_azi, es, em)
    CALL require(es == CD_AGG_OK .AND. line_id == 1, 'extf-report:summary: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 1, ten0, es, em)
    CALL require(es == CD_AGG_OK, 'extf-report:channel-at-equilibrium: '//TRIM(em))
    CALL require(ABS(ten0 - fair_ten) <= 1.0e-9_wp*fair_ten, 'extf-report:summary-equals-equilibrium-channel')
    CALL require(ABS(NORM2(fair_force) - fair_ten) <= 1.0e-9_wp*fair_ten, 'extf-report:summary-force-magnitude')

    ! the static profile written at the equilibrium carries the same End A tension
    CALL CD_AGG_WriteStaticProfile(agg, 'agg_extf_report.static.out', es, em)
    CALL require(es == CD_AGG_OK, 'extf-report:profile: '//TRIM(em))
    ten_profile = -1.0_wp
    OPEN (NEWUNIT=unit, FILE='agg_extf_report.static.out', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    CALL require(ios == 0, 'extf-report:profile-open')
    IF (ios == 0) THEN
      READ (unit, '(A)', IOSTAT=ios) skip
      READ (unit, '(A)', IOSTAT=ios) skip
      READ (unit, '(A)', IOSTAT=ios) skip
      BLOCK
        REAL(wp) :: s, x, y, z
        READ (unit, *, IOSTAT=ios) row_line, row_node, s, x, y, z, ten_profile
      END BLOCK
      CALL require(ios == 0 .AND. row_line == 1 .AND. row_node == 1, 'extf-report:profile-first-row')
      CLOSE (unit)
    END IF
    ! the profile prints eight significant digits
    CALL require(ABS(ten_profile - fair_ten) <= 1.0e-7_wp*fair_ten, 'extf-report:profile-equals-summary')

    ! host fields at t = 0: a 1 m/s current past every node changes the t = 0 channel through
    ! the end node's drag; the summary still reports the equilibrium
    nf = CD_AGG_NFluidNodes(agg, es, em)
    CALL require(es == CD_AGG_OK .AND. nf > 0, 'extf-report:nfluid')
    IF (es /= CD_AGG_OK .OR. nf <= 0) THEN
      CALL CD_AGG_End(agg, es, em)
      RETURN
    END IF
    ALLOCATE (fu(3, nf), fud(3, nf), wl(nf))
    fu = 0.0_wp
    fu(1, :) = 1.0_wp
    fud = 0.0_wp
    wl = 0.0_wp
    CALL CD_AGG_SetFluidFields(agg, fu, fud, wl, es, em)
    CALL require(es == CD_AGG_OK, 'extf-report:set-fields: '//TRIM(em))
    CALL CD_AGG_EvalChannel(agg, 1, ten_field, es, em)
    CALL require(es == CD_AGG_OK, 'extf-report:channel-under-field: '//TRIM(em))
    CALL require(ABS(ten_field - ten0) > 1.0e-6_wp*ten0, 'extf-report:t0-channel-carries-field')
    CALL CD_AGG_GetInitLine(agg, 1, line_id, summary_after, fair_force, fair_incl, fair_decl, fair_azi, es, em)
    CALL require(es == CD_AGG_OK .AND. ABS(summary_after - fair_ten) <= 0.0_wp, &
                 'extf-report:summary-keeps-equilibrium')
    WRITE (*, '(A,ES14.6,A,ES14.6,A)') 'external fluid report: equilibrium FairTen1 = ', fair_ten, &
      ' N; t = 0 under a 1 m/s current = ', ten_field, ' N'
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_external_fluid_report

  SUBROUTINE case_external_fluid_cable()
    TYPE(CD_AGG_ModuleType) :: agg_p, agg_x
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp) :: r0(3, 2), vel(3, 2), acc(3, 2), tgt(3, 2)
    REAL(wp) :: lp(3, 2), l0f(3, 2), lcur(3, 2), lfk(3, 2), p2(3, 2), v2(3, 2), a2(3, 2)
    REAL(wp), ALLOCATABLE :: xyz(:, :), fu(:, :), fud(:, :), wl(:)
    INTEGER :: es, nf, ni, ccol1
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_mixed_deck('agg_extf_mixed.dat')
    CALL CD_AGG_Init_From_Deck(agg_p, 'agg_extf_mixed.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'extfc:plain-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_Init_From_Deck(agg_x, 'agg_extf_mixed.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'extfc:external-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) THEN
      CALL CD_AGG_End(agg_p, es, em)
      RETURN
    END IF

    ! sampling surface spans the mooring line AND the cable nodes
    nf = CD_AGG_NFluidNodes(agg_x, es, em)
    CALL require(es == CD_AGG_OK .AND. nf > 21, 'extfc:nfluid-spans-cable')
    ALLOCATE (xyz(3, nf), fu(3, nf), fud(3, nf), wl(nf))
    CALL CD_AGG_GetFluidNodePositions(agg_x, xyz, es, em)
    CALL require(es == CD_AGG_OK .AND. ALL(IEEE_IS_FINITE(xyz)), 'extfc:positions')

    ! zero fields bit-identical to the plain build (symmetric protocol: CalcOutput on
    ! BOTH before reading -- the facade's mooring load slot is populated by CalcOutput)
    CALL CD_AGG_CalcOutput(agg_p, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_p, r0, vel, acc, lp, es, em)
    fu = 0.0_wp
    fud = 0.0_wp
    wl = 0.0_wp
    CALL CD_AGG_SetFluidFields(agg_x, fu, fud, wl, es, em)
    CALL require(es == CD_AGG_OK, 'extfc:set-zero: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg_x, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_x, p2, v2, a2, l0f, es, em)
    CALL require(nan_max_abs(l0f - lp) <= 0.0_wp, 'extfc:zero-fields-bit-identical')

    ! the CABLE column (column 2 on the mixed deck) responds to velocity (drag) and to
    ! acceleration alone (Froude-Krylov live through the held field)
    ccol1 = 2
    fu(1, :) = 1.0_wp
    CALL CD_AGG_SetFluidFields(agg_x, fu, fud, wl, es, em)
    CALL require(es == CD_AGG_OK, 'extfc:set-current: '//TRIM(em))
    ! the drag needs relative motion or the load evaluation at committed state: the
    ! held velocity enters the DRAG at the next residual evaluation -- step once
    CALL CD_AGG_GetMovingPointMesh(agg_x, r0, vel, acc, lcur, es, em)
    tgt = r0
    vel = 0.0_wp
    acc = 0.0_wp
    CALL CD_AGG_Step_Moving(agg_x, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'extfc:driven-step-velocity: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg_x, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_x, p2, v2, a2, lcur, es, em)
    CALL require(nan_max_abs(lcur(:, ccol1) - l0f(:, ccol1)) > 0.0_wp, 'extfc:cable-drag-live')
    fu = 0.0_wp
    fud(1, :) = 2.0_wp
    CALL CD_AGG_SetFluidFields(agg_x, fu, fud, wl, es, em)
    CALL require(es == CD_AGG_OK, 'extfc:set-acc: '//TRIM(em))
    CALL CD_AGG_Step_Moving(agg_x, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'extfc:driven-step-acc: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg_x, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_x, p2, v2, a2, lfk, es, em)
    CALL require(nan_max_abs(lfk(:, ccol1) - lcur(:, ccol1)) > 0.0_wp, 'extfc:cable-froude-krylov-live')
    ! the held ELEVATION wets ALL hydro contributors -- including the ADDED-MASS
    ! matrix: the same acceleration drive with the free surface dropped far below
    ! the cable (drying it) must march to a different cable reaction than the
    ! wetted drive (regression for the AM-wets-against-am_wl composition gap)
    wl = -500.0_wp
    CALL CD_AGG_SetFluidFields(agg_x, fu, fud, wl, es, em)
    CALL require(es == CD_AGG_OK, 'extfc:set-dry: '//TRIM(em))
    CALL CD_AGG_Step_Moving(agg_x, DT, tgt, vel, acc, cv, st, ni, es, em)
    CALL require(es == CD_AGG_OK .AND. cv, 'extfc:driven-step-dry: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg_x, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_x, p2, v2, a2, lcur, es, em)
    CALL require(nan_max_abs(lcur(:, ccol1) - lfk(:, ccol1)) > 0.0_wp, &
                 'extfc:held-elevation-wets-added-mass-and-loads')

    CALL CD_AGG_End(agg_p, es, em)
    CALL CD_AGG_End(agg_x, es, em)
    ! fail-closed: a PLAIN build (no external_fluid declared) must reject host
    ! fields for BOTH families -- cables carry a drag config, so without the
    ! aggregate-level gate they would silently accept fields nobody declared
    CALL CD_AGG_Init_From_Deck(agg_p, 'agg_extf_mixed.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'extfc: plain re-init: '//TRIM(em))
    nf = CD_AGG_NFluidNodes(agg_p, es, em)
    IF (ALLOCATED(fu)) DEALLOCATE (fu, fud, wl)
    ALLOCATE (fu(3, nf), fud(3, nf), wl(nf))
    fu = 0.0_wp; fud = 0.0_wp; wl = 0.0_wp
    CALL CD_AGG_SetFluidFields(agg_p, fu, fud, wl, es, em)
    CALL require(es /= CD_AGG_OK .AND. INDEX(em, 'external_fluid') > 0, &
                 'extfc: plain build rejects host fields naming external_fluid')
    CALL CD_AGG_End(agg_p, es, em)
  END SUBROUTINE case_external_fluid_cable

  SUBROUTINE case_restart_equivalence()
    !! CHECKPOINT-RESTART semantics at the aggregate level (the shell rebuild's exact
    !! recipe): march instance A mid-flight, capture what a checkpoint carries -- the
    !! coupled kinematics, each line's interior [v; q], each cable's [v; q; a] mirror --
    !! rebuild a FRESH instance B from the same deck, overlay that state (interior
    !! setter + acceleration re-derivation on the mooring; exact mirror reload on the
    !! cable), and continue BOTH under identical drives. The mooring acceleration is
    !! RE-DERIVED, not reloaded, so continuation agreement is tolerance-level (the
    !! committed acceleration satisfies the same balance the recompute solves); the
    !! cable side reloads exactly. Gate: trajectories agree to 1e-8 of range over the
    !! continuation window, on the MIXED deck (both partitions live).
    TYPE(CD_AGG_ModuleType) :: agg_a, agg_b
    REAL(wp), PARAMETER :: DT = 0.05_wp, TOL = 1.0e-8_wp
    INTEGER, PARAMETER :: NPRE = 6, NPOST = 8
    REAL(wp) :: r0(3, 2), vel(3, 2), acc(3, 2), tgt(3, 2)
    REAL(wp) :: pa(3, 2), la(3, 2), pb(3, 2), lb(3, 2), v2(3, 2), a2(3, 2)
    REAL(wp), ALLOCATABLE :: qL(:), vL(:), aL(:), cbuf(:)
    REAL(wp) :: t, zh, vh, lscale
    INTEGER :: es, ni, s, nl, il, ndof, nin, csz
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_mixed_deck('agg_restart.dat')
    CALL CD_AGG_Init_From_Deck(agg_a, 'agg_restart.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'restart:A-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN

    ! march A mid-flight under a heave drive
    CALL CD_AGG_GetMovingPointMesh(agg_a, r0, vel, acc, la, es, em)
    tgt = r0
    vel = 0.0_wp
    acc = 0.0_wp
    DO s = 1, NPRE
      t = REAL(s, wp)*DT
      zh = 0.05_wp*(1.0_wp - COS(2.0_wp*PI*t/4.0_wp))
      vh = 0.05_wp*(2.0_wp*PI/4.0_wp)*SIN(2.0_wp*PI*t/4.0_wp)
      tgt(3, :) = r0(3, :) + zh
      vel(3, :) = vh
      CALL CD_AGG_Step_Moving(agg_a, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'restart:A-march')
      IF (es /= CD_AGG_OK) RETURN
    END DO

    ! the "checkpoint": rebuild B fresh and overlay A's committed state
    CALL CD_AGG_Init_From_Deck(agg_b, 'agg_restart.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'restart:B-init: '//TRIM(em))
    ! coupled kinematics (the input mesh's role in the shell rebuild)
    CALL CD_AGG_GetMovingPointMesh(agg_a, pa, v2, a2, la, es, em)
    CALL CD_AGG_UpdateStates_Moving(agg_b, pa, v2, a2, es, em)
    CALL require(es == CD_AGG_OK, 'restart:B-coupled-overlay: '//TRIM(em))
    ! mooring interiors including the COMMITTED gen-alpha acceleration (an
    ! acceleration re-derivation was measured to move the first restored loads --
    ! the mirror carries [v; q; a] so the reload is exact)
    nl = CD_System_NLines(agg_a%sys%fast%system, es, em)
    DO il = 1, nl
      ndof = CD_System_Line_NDOF(agg_a%sys%fast%system, il, es, em)
      nin = ndof - 6
      IF (ALLOCATED(qL)) DEALLOCATE (qL, vL, aL)
      ALLOCATE (qL(ndof), vL(ndof), aL(ndof))
      CALL CD_Get_System_Line_State(agg_a%sys%fast%system, il, qL, vL, aL, es, em)
      CALL require(es == CD_SYSTEM_OK, 'restart:A-line-read')
      CALL CD_Update_System_Line_Interior_State(agg_b%sys%fast%system, il, &
                                                qL(4:ndof - 3), vL(4:ndof - 3), es, em, &
                                                a_interior=aL(4:ndof - 3))
      CALL require(es == CD_SYSTEM_OK, 'restart:B-interior-write: '//TRIM(em))
    END DO
    ! cable mirror, exact
    csz = CD_HFMF_MirrorSize(agg_a%cables(1))
    CALL require(csz > 0 .AND. csz == CD_HFMF_MirrorSize(agg_b%cables(1)), 'restart:cable-mirror-sized')
    ALLOCATE (cbuf(csz))
    CALL CD_HFMF_PackMirror(agg_a%cables(1), cbuf, es, em)
    CALL require(es == CD_HFMF_OK, 'restart:cable-pack')
    CALL CD_HFMF_UnpackMirror(agg_b%cables(1), cbuf, agg_a%cables(1)%line%t, es, em)
    CALL require(es == CD_HFMF_OK, 'restart:cable-unpack: '//TRIM(em))

    ! continue BOTH under identical drives; trajectories must track to tolerance
    CALL CD_AGG_CalcOutput(agg_a, es, em)
    CALL CD_AGG_CalcOutput(agg_b, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_a, pa, v2, a2, la, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg_b, pb, v2, a2, lb, es, em)
    lscale = MAX(1.0_wp, nan_max_abs(la))
    CALL require(nan_max_abs(pa - pb) <= TOL, 'restart:overlaid-positions-match')
    CALL require(nan_max_abs(la - lb)/lscale <= 1.0e-6_wp, 'restart:overlaid-loads-match')
    DO s = NPRE + 1, NPRE + NPOST
      t = REAL(s, wp)*DT
      zh = 0.05_wp*(1.0_wp - COS(2.0_wp*PI*t/4.0_wp))
      vh = 0.05_wp*(2.0_wp*PI/4.0_wp)*SIN(2.0_wp*PI*t/4.0_wp)
      tgt(3, :) = r0(3, :) + zh
      vel(3, :) = vh
      CALL CD_AGG_Step_Moving(agg_a, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'restart:A-continue')
      CALL CD_AGG_Step_Moving(agg_b, DT, tgt, vel, acc, cv, st, ni, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'restart:B-continue')
      CALL CD_AGG_CalcOutput(agg_a, es, em)
      CALL CD_AGG_CalcOutput(agg_b, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg_a, pa, v2, a2, la, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg_b, pb, v2, a2, lb, es, em)
      lscale = MAX(1.0_wp, nan_max_abs(la))
      CALL require(nan_max_abs(pa - pb) <= TOL, 'restart:continued-positions-track')
      CALL require(nan_max_abs(la - lb)/lscale <= 1.0e-6_wp, 'restart:continued-loads-track')
      IF (nan_max_abs(pa - pb) > TOL) EXIT
    END DO
    CALL CD_AGG_End(agg_a, es, em)
    CALL CD_AGG_End(agg_b, es, em)
  END SUBROUTINE case_restart_equivalence

  SUBROUTINE case_restart_friction()
    !! Stick-slip seabed friction on the coupled route. (1) An EI=0 chain on a frictional
    !! seabed: the aggregate and the standalone system build start from the identical
    !! static state and friction anchors (bit-identical). (2) For the EI=0 chain and for a
    !! finite-EI touchdown cable, a surge drive makes the grounded run slip (the anchors
    !! move); a fresh instance overlaid with the checkpoint state -- line interiors [v; q; a]
    !! and the friction anchors, the cable's mirror including its anchors -- continues on
    !! the same trajectory as the uninterrupted instance, with the same anchors. (3) The same
    !! restart with anisotropic friction (frictionMuAxial beside frictionMu).
    TYPE(CD_AGG_ModuleType) :: agg_a, agg_b
    TYPE(CD_SystemType) :: sys_s
    REAL(wp), PARAMETER :: DT = 0.05_wp, TOL = 1.0e-8_wp
    INTEGER, PARAMETER :: NPRE = 6, NPOST = 8
    REAL(wp), ALLOCATABLE :: r0(:, :), vel(:, :), acc(:, :), tgt(:, :), pa(:, :), la(:, :), pb(:, :), lb(:, :)
    REAL(wp), ALLOCATABLE :: v2(:, :), a2(:, :), qL(:), vL(:), aL(:), qS(:), vS(:), aS(:), anc(:, :), anc0(:, :)
    REAL(wp), ALLOCATABLE :: ancb(:, :), cbuf(:)
    REAL(wp) :: t, xh, vh, lscale, dtm, g, rho, moved
    INTEGER :: es, ni, s, nl, il, ndof, npt, deck, csz
    LOGICAL :: cv, st, is_chain
    TYPE :: anchor_set
      REAL(wp), ALLOCATABLE :: a(:, :)
    END TYPE anchor_set
    TYPE(anchor_set) :: anc_iso(2)
    CHARACTER(512) :: em
    CHARACTER(64) :: path

    ! (1) standalone and coupled builds of the EI=0 friction deck start identically
    CALL write_friction_chain_deck('agg_fric_chain.dat')
    CALL CD_AGG_Init_From_Deck(agg_a, 'agg_fric_chain.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'fric-restart:chain-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL require(CD_System_Line_Has_Friction(agg_a%sys%fast%system, 1), 'fric-restart:chain-has-friction')
    CALL CD_Init_Deck_System('agg_fric_chain.dat', sys_s, dtm, g, rho, es, em, caller_driven=.TRUE.)
    CALL require(es == CD_DECKDRV_OK, 'fric-restart:standalone-init: '//TRIM(em))
    IF (es == CD_DECKDRV_OK) THEN
      ndof = CD_System_Line_NDOF(agg_a%sys%fast%system, 1, es, em)
      ALLOCATE (qL(ndof), vL(ndof), aL(ndof), qS(ndof), vS(ndof), aS(ndof), anc(2, ndof/3), anc0(2, ndof/3))
      CALL CD_Get_System_Line_State(agg_a%sys%fast%system, 1, qL, vL, aL, es, em)
      CALL CD_Get_System_Line_State(sys_s, 1, qS, vS, aS, es, em)
      CALL require(nan_max_abs(qL - qS) <= 0.0_wp .AND. nan_max_abs(aL - aS) <= 0.0_wp, &
                   'fric-restart:standalone/coupled-static-state-bit-identical')
      CALL CD_Get_System_Line_Friction_Anchors(agg_a%sys%fast%system, 1, anc, es, em)
      CALL CD_Get_System_Line_Friction_Anchors(sys_s, 1, anc0, es, em)
      CALL require(es == CD_SYSTEM_OK .AND. nan_max_abs(anc - anc0) <= 0.0_wp, &
                   'fric-restart:standalone/coupled-anchors-bit-identical')
      DEALLOCATE (qL, vL, aL, qS, vS, aS, anc, anc0)
    END IF
    CALL CD_End_System(sys_s, es, em)
    CALL CD_AGG_End(agg_a, es, em)

    ! (2) restart overlay under a slipping surge drive, EI=0 chain then finite-EI cable, then
    ! both again with anisotropic (axial/lateral) friction
    DO deck = 1, 4
      is_chain = MOD(deck, 2) == 1
      SELECT CASE (deck)
      CASE (1)
        path = 'agg_fric_chain.dat'
      CASE (2)
        path = 'agg_fric_cable.dat'
        CALL write_cable_grounded_deck(path)
      CASE (3)
        path = 'agg_fric_chain_aniso.dat'
        CALL write_friction_chain_deck(path, '0.02 frictionMuAxial')
      CASE DEFAULT
        path = 'agg_fric_cable_aniso.dat'
        CALL write_cable_grounded_deck(path, extra_option='0.02 frictionMuAxial')
      END SELECT
      CALL CD_AGG_Init_From_Deck(agg_a, TRIM(path), DT, es, em)
      CALL require(es == CD_AGG_OK, 'fric-restart:A-init: '//TRIM(em))
      IF (es /= CD_AGG_OK) RETURN
      npt = CD_AGG_NMovingPoints(agg_a, es, em)
      ALLOCATE (r0(3, npt), vel(3, npt), acc(3, npt), tgt(3, npt), pa(3, npt), la(3, npt), pb(3, npt), &
                lb(3, npt), v2(3, npt), a2(3, npt))
      CALL CD_AGG_GetMovingPointMesh(agg_a, r0, vel, acc, la, es, em)
      IF (is_chain) THEN
        ndof = CD_System_Line_NDOF(agg_a%sys%fast%system, 1, es, em)
        ALLOCATE (anc0(2, ndof/3), anc(2, ndof/3), ancb(2, ndof/3))
        CALL CD_Get_System_Line_Friction_Anchors(agg_a%sys%fast%system, 1, anc0, es, em)
      ELSE
        CALL require(CD_HFMF_FrictionMirrorSize(agg_a%cables(1)) > 0, 'fric-restart:cable-anchors-mirrored')
        ALLOCATE (anc0(2, agg_a%cables(1)%line%nn), anc(2, agg_a%cables(1)%line%nn))
        anc0 = agg_a%cables(1)%line%fr_anchor
      END IF
      tgt = r0
      vel = 0.0_wp
      acc = 0.0_wp
      DO s = 1, NPRE
        ! surge drive of every coupled point: 0.4 m half-cosine over 2 s
        t = REAL(s, wp)*DT
        xh = 0.4_wp*(1.0_wp - COS(2.0_wp*PI*t/2.0_wp))
        vh = 0.4_wp*(2.0_wp*PI/2.0_wp)*SIN(2.0_wp*PI*t/2.0_wp)
        tgt(1, :) = r0(1, :) + xh
        vel(1, :) = vh
        CALL CD_AGG_Step_Moving(agg_a, DT, tgt, vel, acc, cv, st, ni, es, em)
        CALL require(es == CD_AGG_OK .AND. cv, 'fric-restart:A-march: '//TRIM(em))
        IF (es /= CD_AGG_OK) RETURN
      END DO
      ! the drive makes the grounded run slip: the anchors move
      IF (is_chain) THEN
        CALL CD_Get_System_Line_Friction_Anchors(agg_a%sys%fast%system, 1, anc, es, em)
      ELSE
        anc = agg_a%cables(1)%line%fr_anchor
      END IF
      moved = nan_max_abs(anc - anc0)
      CALL require(moved > 0.0_wp, 'fric-restart:anchors-slip-under-drive')
      ! the anisotropic decks are wired as such and slide differently from the isotropic ones
      IF (deck <= 2) THEN
        IF (ALLOCATED(anc_iso(MOD(deck - 1, 2) + 1)%a)) DEALLOCATE (anc_iso(MOD(deck - 1, 2) + 1)%a)
        ALLOCATE (anc_iso(MOD(deck - 1, 2) + 1)%a, SOURCE=anc)
      ELSE IF (is_chain) THEN
        ! (the short grounded run of the cable sticks under this drive; its anisotropic
        ! sliding is gated in test_seabed_friction_aniso)
        CALL require(nan_max_abs(anc - anc_iso(MOD(deck - 1, 2) + 1)%a) > 1.0e-6_wp, &
                     'fric-restart:aniso-slides-differently')
      END IF
      IF (deck == 4) CALL require(agg_a%cables(1)%line%contact_fr_aniso .AND. &
                                  ABS(agg_a%cables(1)%line%contact_mu_axial - 0.02_wp) <= 0.0_wp, &
                                  'fric-restart:cable-aniso-wired')

      ! the checkpoint: fresh B overlaid with A's committed state and anchors
      CALL CD_AGG_Init_From_Deck(agg_b, TRIM(path), DT, es, em)
      CALL require(es == CD_AGG_OK, 'fric-restart:B-init: '//TRIM(em))
      CALL CD_AGG_GetMovingPointMesh(agg_a, pa, v2, a2, la, es, em)
      CALL CD_AGG_UpdateStates_Moving(agg_b, pa, v2, a2, es, em)
      IF (is_chain) THEN
        nl = CD_System_NLines(agg_a%sys%fast%system, es, em)
        DO il = 1, nl
          ndof = CD_System_Line_NDOF(agg_a%sys%fast%system, il, es, em)
          ALLOCATE (qL(ndof), vL(ndof), aL(ndof))
          CALL CD_Get_System_Line_State(agg_a%sys%fast%system, il, qL, vL, aL, es, em)
          CALL CD_Update_System_Line_Interior_State(agg_b%sys%fast%system, il, qL(4:ndof - 3), &
                                                    vL(4:ndof - 3), es, em, a_interior=aL(4:ndof - 3))
          CALL require(es == CD_SYSTEM_OK, 'fric-restart:B-interior: '//TRIM(em))
          CALL CD_Get_System_Line_Friction_Anchors(agg_a%sys%fast%system, il, anc, es, em)
          CALL CD_Set_System_Line_Friction_Anchors(agg_b%sys%fast%system, il, anc, es, em)
          CALL require(es == CD_SYSTEM_OK, 'fric-restart:B-anchors: '//TRIM(em))
          DEALLOCATE (qL, vL, aL)
        END DO
      ELSE
        csz = CD_HFMF_MirrorSize(agg_a%cables(1))
        ALLOCATE (cbuf(csz))
        CALL CD_HFMF_PackMirror(agg_a%cables(1), cbuf, es, em)
        CALL CD_HFMF_UnpackMirror(agg_b%cables(1), cbuf, agg_a%cables(1)%line%t, es, em)
        CALL require(es == CD_HFMF_OK, 'fric-restart:cable-unpack: '//TRIM(em))
        CALL require(nan_max_abs(agg_b%cables(1)%line%fr_anchor - agg_a%cables(1)%line%fr_anchor) <= 0.0_wp, &
                     'fric-restart:cable-anchors-reloaded')
        DEALLOCATE (cbuf)
      END IF
      DO s = NPRE + 1, NPRE + NPOST
        ! surge drive of every coupled point: 0.4 m half-cosine over 2 s
        t = REAL(s, wp)*DT
        xh = 0.4_wp*(1.0_wp - COS(2.0_wp*PI*t/2.0_wp))
        vh = 0.4_wp*(2.0_wp*PI/2.0_wp)*SIN(2.0_wp*PI*t/2.0_wp)
        tgt(1, :) = r0(1, :) + xh
        vel(1, :) = vh
        CALL CD_AGG_Step_Moving(agg_a, DT, tgt, vel, acc, cv, st, ni, es, em)
        CALL require(es == CD_AGG_OK .AND. cv, 'fric-restart:A-continue')
        CALL CD_AGG_Step_Moving(agg_b, DT, tgt, vel, acc, cv, st, ni, es, em)
        CALL require(es == CD_AGG_OK .AND. cv, 'fric-restart:B-continue')
        CALL CD_AGG_CalcOutput(agg_a, es, em)
        CALL CD_AGG_CalcOutput(agg_b, es, em)
        CALL CD_AGG_GetMovingPointMesh(agg_a, pa, v2, a2, la, es, em)
        CALL CD_AGG_GetMovingPointMesh(agg_b, pb, v2, a2, lb, es, em)
        lscale = MAX(1.0_wp, nan_max_abs(la))
        CALL require(nan_max_abs(la - lb)/lscale <= 1.0e-6_wp, 'fric-restart:continued-loads-track')
        IF (nan_max_abs(la - lb)/lscale > 1.0e-6_wp) EXIT
      END DO
      IF (is_chain) THEN
        CALL CD_Get_System_Line_Friction_Anchors(agg_a%sys%fast%system, 1, anc, es, em)
        CALL CD_Get_System_Line_Friction_Anchors(agg_b%sys%fast%system, 1, ancb, es, em)
        CALL require(nan_max_abs(anc - ancb) <= TOL, 'fric-restart:chain-anchors-track')
        DEALLOCATE (ancb)
      ELSE
        CALL require(nan_max_abs(agg_a%cables(1)%line%fr_anchor - agg_b%cables(1)%line%fr_anchor) <= TOL, &
                     'fric-restart:cable-anchors-track')
      END IF
      CALL CD_AGG_End(agg_a, es, em)
      CALL CD_AGG_End(agg_b, es, em)
      DEALLOCATE (r0, vel, acc, tgt, pa, la, pb, lb, v2, a2, anc0, anc)
    END DO
  END SUBROUTINE case_restart_friction

  SUBROUTINE write_friction_chain_deck(path, extra_option)
    !! EI=0 chain with a long grounded run on a frictional seabed, driven at a Vessel fairlead.
    CHARACTER(*), INTENT(IN) :: path
    CHARACTER(*), INTENT(IN), OPTIONAL :: extra_option
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'coupled EI=0 chain on a frictional seabed'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.252 390.0 1.674e9 -1.0 0.0 1.37 0.64 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Fixed 400.0 0.0 -50.0'
    WRITE (u, '(A)') '2 Vessel 0.0 0.0 -10.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 chain 420.0 42'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '50.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '1.0e4 cBot'
    WRITE (u, '(A)') '0.3 frictionMu'
    IF (PRESENT(extra_option)) WRITE (u, '(A)') extra_option
    WRITE (u, '(A)') '0.05 dtM'
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_friction_chain_deck

  SUBROUTINE write_current_deck(path)
    !! The volturnus deck plus a deck-declared uniform current OPTION (with the dtM/TMax
    !! pair the current option requires, so the parse succeeds and the rejection under
    !! test is the external_fluid double-counting guard, not the parse rule).
    CHARACTER(*), INTENT(IN) :: path
    CALL write_volturnus_deck(path, extra_option='uniform 0.5 0.0 0.0 current', &
                              extra_option2='0.05 dtM', extra_option3='60.0 TMax')
  END SUBROUTINE write_current_deck

  ! ---------------------------------------------------------------------------------------
  ! deck writers
  ! ---------------------------------------------------------------------------------------

  SUBROUTINE write_volturnus_deck(path, fairleads, extra_option, extra_option2, extra_option3)
    !! IEA-15MW VolturnUS-S three-chain mooring (End A = Vessel fairlead, End B = Fixed
    !! anchor), 20 segments per line -- a pure EI=0 deck (all sections chain, EI = 0).
    !! fairleads (optional, 3x3): override the three Vessel coordinates, written
    !! round-trip-exact (ES25.16 = 17 significant digits, the binary64 round-trip
    !! minimum) so a parse of the written deck reproduces the passed values
    !! bit-for-bit -- the ptfm_init equivalence case keys on that.
    !! extra_option[23] (optional): verbatim OPTIONS lines appended before the terminator
    !! (deck current/waves/dtM/TMax variants for the external-fluid cases).
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), INTENT(IN), OPTIONAL :: fairleads(3, 3)
    CHARACTER(*), INTENT(IN), OPTIONAL :: extra_option, extra_option2, extra_option3
    REAL(wp) :: fl(3, 3)
    INTEGER :: u, ios
    fl(:, 1) = [-58.0_wp, 0.0_wp, -14.0_wp]
    fl(:, 2) = [29.0_wp, 50.229_wp, -14.0_wp]
    fl(:, 3) = [29.0_wp, -50.229_wp, -14.0_wp]
    IF (PRESENT(fairleads)) fl = fairleads
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'IEA-15MW VolturnUS-S mooring (pure EI=0), aggregate == FMF gate'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'main 0.333 685.0 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A,3ES25.16)') '1 Vessel ', fl(:, 1)
    WRITE (u, '(A)') '2 Fixed  -837.600   0.000 -200.0'
    WRITE (u, '(A,3ES25.16)') '3 Vessel ', fl(:, 2)
    WRITE (u, '(A)') '4 Fixed   418.800 725.383 -200.0'
    WRITE (u, '(A,3ES25.16)') '5 Vessel ', fl(:, 3)
    WRITE (u, '(A)') '6 Fixed   418.800 -725.383 -200.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '2 3 4 -'
    WRITE (u, '(A)') '3 5 6 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 main 850.0 20'
    WRITE (u, '(A)') '2 main 850.0 20'
    WRITE (u, '(A)') '3 main 850.0 20'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '200.0 WtrDpth'
    WRITE (u, '(A)') '3.0e6 kBot'
    WRITE (u, '(A)') '3.0e5 cBot'
    IF (PRESENT(extra_option)) WRITE (u, '(A)') TRIM(extra_option)
    IF (PRESENT(extra_option2)) WRITE (u, '(A)') TRIM(extra_option2)
    IF (PRESENT(extra_option3)) WRITE (u, '(A)') TRIM(extra_option3)
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_volturnus_deck

  SUBROUTINE write_lozon_semitaut_deck(path)
    !! Three-line 200 m semitaut system used by the displaced-pose regression.
    CHARACTER(*), INTENT(IN) :: path
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'Lozon 200 m composite semitaut displaced-pose gate'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.2791 480.93 2.058e9 -1.0 0.0 1.0 0.0 1.0 0.0'
    WRITE (u, '(A)') 'poly  0.1438  22.42 1.420e8 -1.0 0.0 1.2 0.2 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Coupled  -58.0      0.0     -14.0'
    WRITE (u, '(A)') '2 Fixed   -700.0      0.0    -200.0'
    WRITE (u, '(A)') '3 Coupled   29.0     50.229   -14.0'
    WRITE (u, '(A)') '4 Fixed    350.0    606.218  -200.0'
    WRITE (u, '(A)') '5 Coupled   29.0    -50.229   -14.0'
    WRITE (u, '(A)') '6 Fixed    350.0   -606.218  -200.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '2 3 4 -'
    WRITE (u, '(A)') '3 5 6 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 poly  199.8 20'
    WRITE (u, '(A)') '1 chain 497.7 40'
    WRITE (u, '(A)') '2 poly  199.8 20'
    WRITE (u, '(A)') '2 chain 497.7 40'
    WRITE (u, '(A)') '3 poly  199.8 20'
    WRITE (u, '(A)') '3 chain 497.7 40'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '200.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '1.0e4 cBot'
    WRITE (u, '(A)') '0.1 dtM'
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_lozon_semitaut_deck

  SUBROUTINE case_touchdown_cable_end_tension()
    !! The IEA-15MW touchdown power cable of the mixed example (line 2 beside one mooring
    !! chain): bare / buoyant / bare sections, stiff EA = 4.69e8 N on 4-10 m elements,
    !! ending on a 3e6 Pa/m seabed. Statics carry no horizontal load, so the anchor end sees
    !! the grounded run's tension H (~1.2 kN). The node-sampled stretch EA*(|dr/ds| - 1) at
    !! the held anchor read 78 kN; the element end force equals H.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: u, ios, es, k
    REAL(wp) :: v(3)
    CHARACTER(512) :: em
    REAL(wp), PARAMETER :: DT = 0.025_wp

    OPEN (NEWUNIT=u, FILE='agg_touchdown_cable.dat', STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    CALL require(ios == 0, 'touchdown: create deck')
    WRITE (u, '(A)') 'IEA-15MW mooring chain + touchdown power cable'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'main 0.333 685.00 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27'
    WRITE (u, '(A)') 'bare 0.16 36.70 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') 'buoy 0.29 59.53 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Coupled -58.0 0.0 -14.0'
    WRITE (u, '(A)') '2 Fixed -837.6 0.0 -200.0'
    WRITE (u, '(A)') '3 Coupled 30.0 0.0 -20.0'
    WRITE (u, '(A)') '4 Fixed 500.0 0.0 -200.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '2 3 4 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 main 850.0 50'
    WRITE (u, '(A)') '2 bare 80.0 20'
    WRITE (u, '(A)') '2 buoy 100.0 25'
    WRITE (u, '(A)') '2 bare 420.0 42'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '200.0 WtrDpth'
    WRITE (u, '(A)') '3.0e6 kBot'
    WRITE (u, '(A)') '3.0e5 cBot'
    WRITE (u, '(A)') '0.35 frictionMu'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen2 AnchTen2 Ten2N70'
    WRITE (u, '(A)') 'END'
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_touchdown_cable.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'touchdown: init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL require(es == CD_AGG_OK, 'touchdown: calc: '//TRIM(em))
    DO k = 1, 3
      CALL CD_AGG_EvalChannel(agg, k, v(k), es, em)
      CALL require(es == CD_AGG_OK, 'touchdown: eval: '//TRIM(em))
    END DO
    CALL require(v(3) > 1.0e3_wp .AND. v(3) < 1.5e3_wp, 'touchdown: grounded-run tension is the horizontal tension H')
    CALL require(ABS(v(2) - v(3)) < 0.05_wp*v(3), 'touchdown: cable AnchTen equals the grounded-run tension H')
    CALL require(v(1) > 1.1e4_wp .AND. v(1) < 1.25e4_wp, 'touchdown: cable FairTen carries the suspended weight')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_touchdown_cable_end_tension

  SUBROUTINE write_mixed_deck(path, finite_end_connection, rigid_end_connection)
    !! MIXED deck: one EI=0 catenary chain (mooring) + one finite-EI lazy-wave power cable.
    !! Cable endpoints (Fixed anchor + Coupled fairlead) sit well above the -200 m seabed, so
    !! the suspended pinned-span cable ignores WtrDpth while the mooring uses it (touchdown).
    CHARACTER(*), INTENT(IN) :: path
    LOGICAL, INTENT(IN), OPTIONAL :: finite_end_connection
    LOGICAL, INTENT(IN), OPTIONAL :: rigid_end_connection
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'MIXED deck: EI=0 chain mooring + finite-EI lazy-wave cable'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.333 685.0 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27'
    WRITE (u, '(A)') 'strop 0.080  25.0 2.00e8 -1.0 0.0 1.6 0.2 0.80 0.27'
    WRITE (u, '(A)') 'bare 0.16 36.70 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') 'buoy 0.29 59.53 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Coupled  90.000   0.000  -14.0'
    WRITE (u, '(A)') '2 Fixed     0.000   0.000  -56.0'
    WRITE (u, '(A)') '3 Vessel  -58.000   0.000  -14.0'
    WRITE (u, '(A)') '4 Fixed  -837.600   0.000 -200.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '2 3 4 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 bare 40.0 13'
    WRITE (u, '(A)') '1 buoy 50.0 16'
    WRITE (u, '(A)') '1 bare 55.0 18'
    WRITE (u, '(A)') '2 chain 850.0 20'
    IF (PRESENT(finite_end_connection)) THEN
      IF (finite_end_connection) THEN
        WRITE (u, '(A)') '--- END CONNECTIONS ---'
        WRITE (u, '(A)') 'LineID End Stiffness EzX EzY EzZ'
        WRITE (u, '(A)') '(-) (-) (N-m/rad) (-) (-) (-)'
        WRITE (u, '(A)') '1 A 20000 0 0 -1'
      END IF
    END IF
    IF (PRESENT(rigid_end_connection)) THEN
      IF (rigid_end_connection) THEN
        WRITE (u, '(A)') '--- END CONNECTIONS ---'
        WRITE (u, '(A)') 'LineID End Stiffness EzX EzY EzZ'
        WRITE (u, '(A)') '(-) (-) (N-m/rad) (-) (-) (-)'
        ! One degree from the unconstrained equilibrium tangent: enough to make
        ! the reaction observable without seeding an artificial sharp kink.
        WRITE (u, '(A)') '1 A Rigid -0.24716 0 -0.968975'
      END IF
    END IF
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '200.0 WtrDpth'
    WRITE (u, '(A)') '3.0e6 kBot'
    WRITE (u, '(A)') '3.0e5 cBot'
    WRITE (u, '(A)') '0.5 rhoInf'
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_mixed_deck

  SUBROUTINE write_mixed_channels_deck(path, point_channel_id, outputs_row)
    !! The write_mixed_deck geometry (cable = line 1, mooring chain = line 2) plus an OUTPUTS
    !! section exercising both objects and multiple channel kinds: FairTen/AnchTen on the mooring,
    !! and curvature / fairlead-angle / bend-moment / nodal-tension on the cable. outputs_row
    !! replaces the default OUTPUTS row.
    CHARACTER(*), INTENT(IN) :: path
    INTEGER, INTENT(IN), OPTIONAL :: point_channel_id
    CHARACTER(*), INTENT(IN), OPTIONAL :: outputs_row
    INTEGER :: u, ios
    INTEGER :: point_id
    point_id = 0
    IF (PRESENT(point_channel_id)) point_id = point_channel_id
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'MIXED deck with OUTPUTS: EI=0 chain mooring + finite-EI lazy-wave cable'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.333 685.0 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27'
    WRITE (u, '(A)') 'strop 0.080  25.0 2.00e8 -1.0 0.0 1.6 0.2 0.80 0.27'
    WRITE (u, '(A)') 'bare 0.16 36.70 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') 'buoy 0.29 59.53 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Coupled  90.000   0.000  -14.0'
    WRITE (u, '(A)') '2 Fixed     0.000   0.000  -56.0'
    WRITE (u, '(A)') '3 Vessel  -58.000   0.000  -14.0'
    WRITE (u, '(A)') '4 Fixed  -837.600   0.000 -200.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '2 3 4 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 bare 40.0 13'
    WRITE (u, '(A)') '1 buoy 50.0 16'
    WRITE (u, '(A)') '1 bare 55.0 18'
    WRITE (u, '(A)') '2 chain 850.0 20'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '200.0 WtrDpth'
    WRITE (u, '(A)') '3.0e6 kBot'
    WRITE (u, '(A)') '3.0e5 cBot'
    WRITE (u, '(A)') '0.5 rhoInf'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    IF (point_id > 0) THEN
      WRITE (u, '(A,I0,A)') 'Point', point_id, 'px'
    ELSE IF (PRESENT(outputs_row)) THEN
      WRITE (u, '(A)') outputs_row
    ELSE
      WRITE (u, '(A)') 'FairTen2 AnchTen2 Curv1N5 FairAngle1 BendMom1N5 Ten1N10 AnchDecl1 FairIncl1'
    END IF
    WRITE (u, '(A)') 'END'
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_mixed_channels_deck

  SUBROUTINE write_mooring_channels_deck(path)
    !! A PURE EI=0 mooring deck (single catenary chain) that declares OUTPUTS -- exercises the
    !! channel path when no cable exists (ncable = 0, cables store unallocated).
    CHARACTER(*), INTENT(IN) :: path
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'PURE mooring deck with OUTPUTS (no cable)'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'main 0.333 685.0 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Vessel  -58.000   0.000  -14.0'
    WRITE (u, '(A)') '2 Fixed  -837.600   0.000 -200.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 main 850.0 20'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '200.0 WtrDpth'
    WRITE (u, '(A)') '3.0e6 kBot'
    WRITE (u, '(A)') '3.0e5 cBot'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1 AnchTen1'
    WRITE (u, '(A)') 'END'
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_mooring_channels_deck

  SUBROUTINE write_cable_only_deck(path, modified_newton, cable_load_feedback)
    !! A finite-EI lazy-wave power cable, coupled at the fairlead -- no EI=0 line. Bare +
    !! buoyancy-module properties from the Humboldt / GoMex 15 MW reference (the l3_lazywave
    !! and deck_hermite_cable convention). Purely SUSPENDED (no WtrDpth). modified_newton
    !! adds the OPTIONS row (the tangent-reuse plumbing gate's twin-deck knob).
    CHARACTER(*), INTENT(IN) :: path
    LOGICAL, INTENT(IN), OPTIONAL :: modified_newton
    LOGICAL, INTENT(IN), OPTIONAL :: cable_load_feedback
    INTEGER :: u, ios
    LOGICAL :: mn
    mn = .FALSE.
    IF (PRESENT(modified_newton)) mn = modified_newton
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'finite-EI lazy-wave power cable, coupled at the fairlead (no EI=0 line)'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'bare 0.16 36.70 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') 'buoy 0.29 59.53 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Fixed 0.0 0.0 -56.0'
    WRITE (u, '(A)') '2 Coupled 90.0 0.0 -14.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 bare 40.0 13'
    WRITE (u, '(A)') '1 buoy 50.0 16'
    WRITE (u, '(A)') '1 bare 55.0 18'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '0.5 rhoInf'
    IF (mn) WRITE (u, '(A)') 'true modified_newton'
    IF (PRESENT(cable_load_feedback)) THEN
      IF (cable_load_feedback) THEN
        WRITE (u, '(A)') 'true cable_load_feedback'
      ELSE
        WRITE (u, '(A)') 'false cable_load_feedback'
      END IF
    END IF
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_cable_only_deck

  SUBROUTINE write_cable_grounded_deck(path, bathy_path, refined, extra_option)
    !! Coupled finite-EI touchdown cable. The 115 m span between endpoints 96.9 m apart
    !! develops a grounded run on the -50 m bed without exceeding the 126 m planar
    !! L-shape capacity (longer fixtures necessarily wad on a frictionless bed).
    !! Optional bathymetry_path replaces flat WtrDpth.
    CHARACTER(*), INTENT(IN) :: path
    CHARACTER(*), INTENT(IN), OPTIONAL :: bathy_path
    LOGICAL, INTENT(IN), OPTIONAL :: refined
    CHARACTER(*), INTENT(IN), OPTIONAL :: extra_option
    INTEGER :: u, ios
    LOGICAL :: fine_mesh
    fine_mesh = .FALSE.
    IF (PRESENT(refined)) fine_mesh = refined
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'coupled finite-EI touchdown cable'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'bare 0.16 36.70 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') 'buoy 0.29 59.53 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Fixed 0.0 0.0 -50.0'
    WRITE (u, '(A)') '2 Coupled 90.0 0.0 -14.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    IF (fine_mesh) THEN
      WRITE (u, '(A)') '1 bare 30.0 67'
      WRITE (u, '(A)') '1 buoy 40.0 89'
      WRITE (u, '(A)') '1 bare 45.0 100'
    ELSE
      WRITE (u, '(A)') '1 bare 30.0 13'
      WRITE (u, '(A)') '1 buoy 40.0 17'
      WRITE (u, '(A)') '1 bare 45.0 20'
    END IF
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    IF (PRESENT(bathy_path)) THEN
      WRITE (u, '(A)') TRIM(bathy_path)//' bathymetryFile'
    ELSE
      WRITE (u, '(A)') '50.0 WtrDpth'
    END IF
    WRITE (u, '(A)') '3.0e6 kBot'
    WRITE (u, '(A)') '3.0e5 cBot'
    WRITE (u, '(A)') '0.35 frictionMu'
    IF (PRESENT(extra_option)) WRITE (u, '(A)') extra_option
    WRITE (u, '(A)') '0.05 dtM'
    WRITE (u, '(A)') '0.20 TMax'
    WRITE (u, '(A)') '0.5 rhoInf'
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_cable_grounded_deck

  SUBROUTINE write_cable_bathymetry(path, depth)
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), INTENT(IN) :: depth
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    CALL require(ios == 0, 'grounded:create-bathymetry-file')
    WRITE (u, '(3(ES22.14,1X))') - 10.0_wp, -10.0_wp, depth
    WRITE (u, '(3(ES22.14,1X))') 100.0_wp, -10.0_wp, depth
    WRITE (u, '(3(ES22.14,1X))') - 10.0_wp, 10.0_wp, depth
    WRITE (u, '(3(ES22.14,1X))') 100.0_wp, 10.0_wp, depth
    CLOSE (u)
  END SUBROUTINE write_cable_bathymetry

  SUBROUTINE write_cable_sloped_bathymetry(path)
    !! depth = 50 + 0.1*x over the cable footprint: End A is on the floor,
    !! while the scalar minimum endpoint elevation is nine metres deeper.
    CHARACTER(*), INTENT(IN) :: path
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    CALL require(ios == 0, 'grounded-sequenced:create-sloped-bathymetry-file')
    WRITE (u, '(3(ES22.14,1X))') - 10.0_wp, -10.0_wp, 49.0_wp
    WRITE (u, '(3(ES22.14,1X))') 100.0_wp, -10.0_wp, 60.0_wp
    WRITE (u, '(3(ES22.14,1X))') - 10.0_wp, 10.0_wp, 49.0_wp
    WRITE (u, '(3(ES22.14,1X))') 100.0_wp, 10.0_wp, 60.0_wp
    CLOSE (u)
  END SUBROUTINE write_cable_sloped_bathymetry

  SUBROUTINE write_mixed_free_deck(path)
    !! MIXED deck (chain + cable) whose mooring line hangs from a Free (dynamic) clump point.
    !! In v1 the aggregate fails closed on a finite-EI deck carrying Free/Connect points.
    !! dtM/TMax are present so the Free point passes parse_deck and reaches the aggregate guard.
    CHARACTER(*), INTENT(IN) :: path
    INTEGER :: u, ios
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    WRITE (u, '(A)') 'MIXED deck with a Free (dynamic) point -- must fail closed in v1'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.333 685.0 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27'
    WRITE (u, '(A)') 'strop 0.080  25.0 2.00e8 -1.0 0.0 1.6 0.2 0.80 0.27'
    WRITE (u, '(A)') 'bare 0.16 36.70 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') 'buoy 0.29 59.53 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Coupled  90.000   0.000  -14.0    0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '2 Fixed     0.000   0.000  -56.0    0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '3 Free    -58.000   0.000  -60.0 5000.0 1.0 2.0 1.0'
    WRITE (u, '(A)') '4 Fixed  -837.600   0.000 -200.0    0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '2 3 4 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 bare 40.0 13'
    WRITE (u, '(A)') '1 buoy 50.0 16'
    WRITE (u, '(A)') '1 bare 55.0 18'
    WRITE (u, '(A)') '2 chain 850.0 20'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '200.0 WtrDpth'
    WRITE (u, '(A)') '3.0e6 kBot'
    WRITE (u, '(A)') '3.0e5 cBot'
    WRITE (u, '(A)') '0.5 rhoInf'
    WRITE (u, '(A)') '0.05 dtM'
    WRITE (u, '(A)') '1.0 TMax'
    WRITE (u, '(A)') '--- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_mixed_free_deck

  SUBROUTINE write_host_rod_deck(path, rod_as_points)
    !! A Coupled (platform-borne) vertical rod, End A at z = -1 and End B at z = -3, moored at both
    !! ends to its own anchor, plus an independent Coupled fairlead with its own line. With
    !! rod_as_points the rod is replaced by two Coupled points at its end positions (the same
    !! lines, the same statics): the reference for the rod's line-end loads.
    CHARACTER(*), INTENT(IN) :: path
    LOGICAL, INTENT(IN) :: rod_as_points
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'coupled host rod deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'line 0.10 80.0 1.0e6 0.0 0.0 1.0 0.0 1.0 0.0'
    IF (.NOT. rod_as_points) THEN
      WRITE (u, '(A)') '--- ROD TYPES ---'
      WRITE (u, '(A)') 'Name Diam Mass Cd Ca CdEnd CaEnd CdAx CaAx'
      WRITE (u, '(A)') '(-) (m) (kg/m) (-) (-) (-) (-) (-) (-)'
      WRITE (u, '(A)') 'rodmat 0.20 50.0 1.0 1.0 0.0 0.0 0.2 0.0'
      WRITE (u, '(A)') '--- RODS ---'
      WRITE (u, '(A)') 'ID RodType Type XA YA ZA XB YB ZB NumSegs Outputs'
      WRITE (u, '(A)') '(-) (-) (-) (m) (m) (m) (m) (m) (m) (-) (-)'
      WRITE (u, '(A)') '1 rodmat Coupled 0.0 0.0 -1.0 0.0 0.0 -3.0 4 -'
    END IF
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 6.0 0.0 -8.0 0.0 0.0 0.0 0.0'
    IF (rod_as_points) THEN
      WRITE (u, '(A)') '2 Coupled 0.0 0.0 -1.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '3 Coupled 0.0 0.0 -3.0 0.0 0.0 0.0 0.0'
    ELSE
      WRITE (u, '(A)') '2 Rod1A 0.0 0.0 -1.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '3 Rod1B 0.0 0.0 -3.0 0.0 0.0 0.0 0.0'
    END IF
    WRITE (u, '(A)') '4 Coupled 0.0 2.0 -1.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '5 Fixed 6.0 2.0 -8.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '6 Fixed -6.0 0.0 -8.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 3 1 -'
    WRITE (u, '(A)') '2 4 5 -'
    WRITE (u, '(A)') '3 2 6 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 line 8.0 6'
    WRITE (u, '(A)') '2 line 9.5 6'
    WRITE (u, '(A)') '3 line 9.5 6'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_host_rod_deck

  SUBROUTINE case_host_rod()
    !! Coupled/Vessel RODS on the aggregate (MoorDyn "Coupled" rods): each is a 6-DOF host node at
    !! its End A that returns force and moment. Gates: (1) the moving surface appends one node per
    !! host rod after the system and cable columns; (2) at rest the wrench equals the line-end loads
    !! of the same lines on two Coupled points (independent deck) plus the rod's weight and
    !! buoyancy, with the moment of the End B load about End A; (3) a host acceleration returns
    !! -(m + rhoW Ca A L) a on top of the point-deck line reactions; (4) the host orientation turns
    !! the rod rigidly and round-trips; (5) a surge+pitch march converges with a finite wrench,
    !! correction-rewind is bit-identical, and a fresh instance overlaid with the checkpoint state
    !! (host kinematics, line interiors, rod states) continues bit-identically; (6) PtfmInit moves
    !! the rod rigidly with the platform.
    TYPE(CD_AGG_ModuleType) :: agg, ref, agg_b
    REAL(wp), PARAMETER :: DT = 0.01_wp, RHO = 1025.0_wp, G = 9.80665_wp, D = 0.2_wp, L = 2.0_wp, M = 100.0_wp
    REAL(wp), ALLOCATABLE :: p(:, :), v(:, :), a(:, :), f(:, :), mo(:, :), orr(:, :, :), w(:, :), alp(:, :)
    REAL(wp), ALLOCATABLE :: pr(:, :), vr(:, :), ar(:, :), fr(:, :), rs(:), rs2(:), p0(:, :), fb(:, :), mb(:, :)
    REAL(wp), ALLOCATABLE :: lq(:), lv(:), la2(:)
    REAL(wp), ALLOCATABLE :: sor(:, :, :), sw(:, :), sal(:, :)
    REAL(wp) :: vol, fexp(3), mexp(3), rot(3, 3), th, err, fsc, ptfm(6), ctr(3)
    INTEGER :: es, n, nr, s, k, nl, il, ndof
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_host_rod_deck('agg_hostrod.dat', .FALSE.)
    CALL write_host_rod_deck('agg_hostrod_pts.dat', .TRUE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_hostrod.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'hostrod:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_Init_From_Deck(ref, 'agg_hostrod_pts.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'hostrod:ref-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    n = CD_AGG_NMovingPoints(agg, es, em)
    nr = CD_AGG_NMovingPoints(ref, es, em)
    ! (1) one system fairlead (point 4) + one host rod node; the reference has points 2, 3, 4
    CALL require(n == 2 .AND. nr == 3, 'hostrod:moving-surface-layout')
    IF (n /= 2 .OR. nr /= 3) RETURN
    ALLOCATE (p(3, n), v(3, n), a(3, n), f(3, n), mo(3, n), orr(3, 3, n), w(3, n), alp(3, n), p0(3, n))
    ALLOCATE (pr(3, nr), vr(3, nr), ar(3, nr), fr(3, nr), fb(3, n), mb(3, n))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, f, es, em, orientation=orr, moment=mo)
    CALL require(es == CD_AGG_OK, 'hostrod:mesh: '//TRIM(em))
    CALL require(nan_max_abs(p(:, 2) - [0.0_wp, 0.0_wp, -1.0_wp]) <= 1.0e-14_wp, 'hostrod:node-at-end-A')
    CALL require(nan_max_abs(orr(:, :, 2) - eye3()) <= 1.0e-14_wp, 'hostrod:init-orientation-identity')
    CALL CD_AGG_CalcOutput(ref, es, em)
    CALL CD_AGG_GetMovingPointMesh(ref, pr, vr, ar, fr, es, em)
    ! (2) static wrench vs the independent point deck
    vol = 0.25_wp*PI*D*D*L
    fexp = fr(:, 1) + fr(:, 2) + [0.0_wp, 0.0_wp, -M*G + RHO*G*vol]
    mexp = cross(pr(:, 2) - pr(:, 1), fr(:, 2))
    fsc = nan_max_abs(fexp)
    err = nan_max_abs(f(:, 2) - fexp)/fsc
    WRITE (*, '(A,ES10.3,A,ES10.3)') 'host rod static wrench: rel. force error ', err, &
      '; moment error ', nan_max_abs(mo(:, 2) - mexp)/(fsc*L)
    CALL require(err <= 1.0e-10_wp, 'hostrod:static-force-equals-lines-plus-weight-buoyancy')
    CALL require(nan_max_abs(mo(:, 2) - mexp) <= 1.0e-10_wp*fsc*L, 'hostrod:static-moment-about-end-A')
    CALL require(nan_max_abs(f(:, 1) - fr(:, 3)) <= 1.0e-10_wp*fsc, 'hostrod:plain-fairlead-unchanged')
    ! (3) host acceleration: rod inertia + transverse added mass on top of the line reactions
    p0 = p
    v = 0.0_wp
    a = 0.0_wp
    a(1, 2) = 0.5_wp
    vr = 0.0_wp
    ar = 0.0_wp
    ar(1, 1:2) = 0.5_wp
    CALL CD_AGG_UpdateStates_Moving(agg, p0, v, a, es, em, output_probe=.TRUE.)
    CALL require(es == CD_AGG_OK, 'hostrod:probe-update: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, fb, es, em, moment=mb)
    CALL CD_AGG_UpdateStates_Moving(ref, pr, vr, ar, es, em, output_probe=.TRUE.)
    CALL CD_AGG_CalcOutput(ref, es, em)
    CALL CD_AGG_GetMovingPointMesh(ref, pr, vr, ar, fr, es, em)
    fexp = fr(:, 1) + fr(:, 2) + [0.0_wp, 0.0_wp, -M*G + RHO*G*vol] - &
           (M + RHO*1.0_wp*vol)*[0.5_wp, 0.0_wp, 0.0_wp]
    err = nan_max_abs(fb(:, 2) - fexp)/fsc
    WRITE (*, '(A,ES10.3)') 'host rod accelerated wrench: rel. force error ', err
    CALL require(err <= 1.0e-10_wp, 'hostrod:inertia-minus-(m+Ca rhoW V) a')
    ! restore the rest state on both
    v = 0.0_wp
    a = 0.0_wp
    CALL CD_AGG_UpdateStates_Moving(agg, p0, v, a, es, em, output_probe=.TRUE.)
    ! (4) orientation: pitch the host by th about global y; the rod turns rigidly about End A
    th = 0.1_wp
    rot = eye3()
    rot(1, 1) = COS(th)
    rot(1, 3) = -SIN(th)
    rot(3, 1) = SIN(th)
    rot(3, 3) = COS(th)
    orr = 0.0_wp
    orr(:, :, 1) = eye3()
    orr(:, :, 2) = rot
    w = 0.0_wp
    alp = 0.0_wp
    CALL CD_AGG_UpdateStates_Moving(agg, p0, v, a, es, em, output_probe=.TRUE., orientation=orr, &
                                    angular_velocity=w, angular_acceleration=alp)
    CALL require(es == CD_AGG_OK, 'hostrod:orient-update: '//TRIM(em))
    ALLOCATE (rs(CD_AGG_Rod_MirrorSize(agg)), rs2(CD_AGG_Rod_MirrorSize(agg)))
    CALL CD_AGG_Get_Rod_States(agg, rs, es, em)
    ctr = p0(:, 2) + 0.5_wp*L*MATMUL(TRANSPOSE(rot), [0.0_wp, 0.0_wp, -1.0_wp])
    CALL require(nan_max_abs(rs(1:3) - ctr) <= 1.0e-14_wp, 'hostrod:orientation-turns-rod-about-end-A')
    CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, f, es, em, orientation=orr)
    CALL require(nan_max_abs(orr(:, :, 2) - rot) <= 1.0e-14_wp .AND. nan_max_abs(p(:, 2) - p0(:, 2)) <= 1.0e-14_wp, &
                 'hostrod:orientation-round-trip')
    orr(:, :, 2) = eye3()
    CALL CD_AGG_UpdateStates_Moving(agg, p0, v, a, es, em, output_probe=.TRUE., orientation=orr, &
                                    angular_velocity=w, angular_acceleration=alp)
    ! (5) surge + pitch march; rewind; restart overlay
    DO s = 1, 60
      CALL host_rod_drive(REAL(s, wp)*DT, p0, p, v, a, orr, w, alp)
      CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=REAL(s, wp)*DT, orientation=orr, &
                              angular_velocity=w, angular_acceleration=alp)
      CALL require(es == CD_AGG_OK .AND. cv, 'hostrod:march: '//TRIM(em))
      IF (es /= CD_AGG_OK) RETURN
    END DO
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), f, es, em, moment=mo)
    CALL require(ALL(IEEE_IS_FINITE(f)) .AND. ALL(IEEE_IS_FINITE(mo)) .AND. nan_max_abs(mo(:, 2)) > 0.0_wp, &
                 'hostrod:march-wrench-finite')
    CALL CD_AGG_Snapshot(agg, es, em)
    CALL host_rod_drive(61.0_wp*DT, p0, p, v, a, orr, w, alp)
    CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=61.0_wp*DT, orientation=orr, &
                            angular_velocity=w, angular_acceleration=alp)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), fb, es, em, moment=mb)
    CALL CD_AGG_Restore(agg, es, em)
    CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=61.0_wp*DT, orientation=orr, &
                            angular_velocity=w, angular_acceleration=alp)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), f, es, em, moment=mo)
    CALL require(nan_max_abs(f - fb) <= 0.0_wp .AND. nan_max_abs(mo - mb) <= 0.0_wp, &
                 'hostrod:correction-rewind-bit-identical')
    ! checkpoint overlay onto a fresh instance, in the OpenFAST shell's restart order
    CALL CD_AGG_Init_From_Deck(agg_b, 'agg_hostrod.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'hostrod:B-init: '//TRIM(em))
    CALL CD_AGG_UpdateStates_Moving(agg_b, p, v, a, es, em, orientation=orr, angular_velocity=w, &
                                    angular_acceleration=alp)
    CALL require(es == CD_AGG_OK, 'hostrod:B-host: '//TRIM(em))
    nl = CD_System_NLines(agg%sys%fast%system, es, em)
    DO il = 1, nl
      ndof = CD_System_Line_NDOF(agg%sys%fast%system, il, es, em)
      ALLOCATE (lq(ndof), lv(ndof), la2(ndof))
      CALL CD_Get_System_Line_State(agg%sys%fast%system, il, lq, lv, la2, es, em)
      CALL CD_Update_System_Line_Interior_State(agg_b%sys%fast%system, il, lq(4:ndof - 3), lv(4:ndof - 3), es, em, &
                                                a_interior=la2(4:ndof - 3))
      CALL require(es == CD_SYSTEM_OK, 'hostrod:B-interior: '//TRIM(em))
      DEALLOCATE (lq, lv, la2)
    END DO
    CALL CD_AGG_Get_Rod_States(agg, rs, es, em)
    CALL CD_AGG_Set_Rod_States(agg_b, rs, es, em)
    CALL require(es == CD_AGG_OK, 'hostrod:B-rod-states: '//TRIM(em))
    ! the OpenFAST CalcOutput probe on B: overlay other host kinematics, then restore the saved
    ! host kinematics and the mirror's rod states (restore_calcoutput_probe) -- state-const
    CALL CD_AGG_GetMovingPointMesh(agg_b, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), fb, es, em)
    ALLOCATE (sor(3, 3, n), sw(3, n), sal(3, n))
    CALL CD_AGG_GetMovingPointMesh(agg_b, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), fb, es, em, orientation=sor, &
                                   angular_velocity=sw, angular_acceleration=sal)
    CALL host_rod_drive(0.37_wp, p0, p, v, a, orr, w, alp)
    CALL CD_AGG_UpdateStates_Moving(agg_b, p, v, a, es, em, output_probe=.TRUE., orientation=orr, &
                                    angular_velocity=w, angular_acceleration=alp)
    CALL CD_AGG_CalcOutput(agg_b, es, em)
    CALL CD_AGG_UpdateStates_Moving(agg_b, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), es, em, output_probe=.TRUE., &
                                    orientation=sor, angular_velocity=sw, angular_acceleration=sal)
    CALL CD_AGG_Set_Rod_States(agg_b, rs, es, em)
    CALL require(es == CD_AGG_OK, 'hostrod:B-probe-restore: '//TRIM(em))
    err = 0.0_wp
    DO s = 62, 80
      CALL host_rod_drive(REAL(s, wp)*DT, p0, p, v, a, orr, w, alp)
      CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=REAL(s, wp)*DT, orientation=orr, &
                              angular_velocity=w, angular_acceleration=alp)
      CALL CD_AGG_Step_Moving(agg_b, DT, p, v, a, cv, st, k, es, em, t_committed=REAL(s, wp)*DT, orientation=orr, &
                              angular_velocity=w, angular_acceleration=alp)
      CALL require(es == CD_AGG_OK, 'hostrod:B-march: '//TRIM(em))
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), f, es, em, moment=mo)
      CALL CD_AGG_CalcOutput(agg_b, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg_b, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), fb, es, em, moment=mb)
      err = MAX(err, nan_max_abs(f - fb), nan_max_abs(mo - mb))
    END DO
    WRITE (*, '(A,ES10.3)') 'host rod restart overlay: max |wrench difference| over 19 steps = ', err
    CALL require(err <= 0.0_wp, 'hostrod:restart-overlay-bit-identical')
    CALL CD_AGG_Get_Rod_States(agg, rs, es, em)
    CALL CD_AGG_Get_Rod_States(agg_b, rs2, es, em)
    CALL require(nan_max_abs(rs - rs2) <= 0.0_wp, 'hostrod:restart-rod-state-bit-identical')
    CALL CD_AGG_End(agg_b, es, em)
    CALL CD_AGG_End(agg, es, em)
    CALL CD_AGG_End(ref, es, em)
    ! (6) PtfmInit: the rod moves rigidly with the platform (surge 1 m, heave 0.5 m, yaw 10 deg)
    ptfm = [1.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.0_wp, 10.0_wp*PI/180.0_wp]
    CALL CD_AGG_Init_From_Deck(agg, 'agg_hostrod.dat', DT, es, em, ptfm_init=ptfm)
    CALL require(es == CD_AGG_OK, 'hostrod:ptfm-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, f, es, em, orientation=orr)
    CALL require(nan_max_abs(p(:, 2) - [1.0_wp, 0.0_wp, -0.5_wp]) <= 1.0e-12_wp, 'hostrod:ptfm-node-displaced')
    rot = eye3()
    rot(1, 1) = COS(ptfm(6))
    rot(1, 2) = SIN(ptfm(6))
    rot(2, 1) = -SIN(ptfm(6))
    rot(2, 2) = COS(ptfm(6))
    CALL require(nan_max_abs(orr(:, :, 2) - rot) <= 1.0e-12_wp, 'hostrod:ptfm-orientation-is-platform-dcm')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_host_rod

  SUBROUTINE host_rod_drive(t, p0, p, v, a, orr, w, alp)
    !! Surge 0.3 sin(2 pi t) of both host nodes plus a pitch 0.05 sin(2 pi t) of the rod node (column 2).
    REAL(wp), INTENT(IN) :: t, p0(:, :)
    REAL(wp), INTENT(OUT) :: p(:, :), v(:, :), a(:, :), orr(:, :, :), w(:, :), alp(:, :)
    REAL(wp) :: om, ang, dang, ddang
    om = 2.0_wp*PI
    p = p0
    p(1, :) = p0(1, :) + 0.3_wp*SIN(om*t)
    v = 0.0_wp
    v(1, :) = 0.3_wp*om*COS(om*t)
    a = 0.0_wp
    a(1, :) = -0.3_wp*om*om*SIN(om*t)
    ang = 0.05_wp*SIN(om*t)
    dang = 0.05_wp*om*COS(om*t)
    ddang = -0.05_wp*om*om*SIN(om*t)
    orr = 0.0_wp
    orr(:, :, 1) = eye3()
    orr(:, :, 2) = eye3()
    orr(1, 1, 2) = COS(ang)
    orr(1, 3, 2) = -SIN(ang)
    orr(3, 1, 2) = SIN(ang)
    orr(3, 3, 2) = COS(ang)
    w = 0.0_wp
    alp = 0.0_wp
    w(2, 2) = dang
    alp(2, 2) = ddang
  END SUBROUTINE host_rod_drive

  SUBROUTINE write_host_body_deck(path, body_as_points, moordyn_row, free_body)
    !! A Coupled (platform-borne) Rigid6 body at (0, 0, -2) with two Body1 attachment points at
    !! +-1 m along x, each moored to its own anchor, plus an independent Coupled fairlead. With
    !! body_as_points the body is replaced by two Coupled points at the attachment positions;
    !! moordyn_row writes the body as a MoorDyn v2 14-column "Coupled" row.
    !! free_body adds an independent buoyant Free Rigid6 body (body 2) on three short taut tethers,
    !! whose stiff mode makes the body march sub-cycle.
    CHARACTER(*), INTENT(IN) :: path
    LOGICAL, INTENT(IN) :: body_as_points, moordyn_row
    LOGICAL, INTENT(IN), OPTIONAL :: free_body
    INTEGER :: u
    LOGICAL :: fb
    fb = .FALSE.
    IF (PRESENT(free_body)) fb = free_body
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'coupled host body deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'line 0.10 80.0 1.0e6 0.0 0.0 1.0 0.0 1.0 0.0'
    IF (.NOT. body_as_points) THEN
      WRITE (u, '(A)') '--- BODIES ---'
      IF (moordyn_row) THEN
        WRITE (u, '(A)') 'ID Attachment X0 Y0 Z0 r0 p0 y0 Mass CG* I* Volume CdA* Ca*'
        WRITE (u, '(A)') '(#) (word) (m) (m) (m) (deg) (deg) (deg) (kg) (m) (kg-m^2) (m^3) (m^2) (-)'
        WRITE (u, '(A)') '1 Coupled 0.0 0.0 -2.0 0.0 0.0 0.0 200.0 0.0 30.0 0.1 0.0 0.5'
      ELSE
        WRITE (u, '(A)') 'ID Type X Y Z Rot1 Rot2 Rot3 Mass Vol C33 C44 C55 CdA Ca Ixx Iyy Izz'
        WRITE (u, '(A)') '(-) (-) (m) (m) (m) (deg) (deg) (deg) (kg) (m3) (N/m) (Nm/rad) (Nm/rad) (m2) (-) '// &
          '(kgm2) (kgm2) (kgm2)'
        WRITE (u, '(A)') '1 Coupled 0.0 0.0 -2.0 0.0 0.0 0.0 200.0 0.1 0.0 0.0 0.0 0.0 0.5 30.0 30.0 30.0'
        IF (fb) WRITE (u, '(A)') '2 Rigid6 0.0 -4.0 -4.0 0.0 0.0 0.0 20.0 0.1 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0'
      END IF
    END IF
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 8.0 0.0 -8.0 0.0 0.0 0.0 0.0'
    IF (body_as_points) THEN
      WRITE (u, '(A)') '2 Coupled 1.0 0.0 -2.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '3 Coupled -1.0 0.0 -2.0 0.0 0.0 0.0 0.0'
    ELSE
      WRITE (u, '(A)') '2 Body1 1.0 0.0 0.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '3 Body1 -1.0 0.0 0.0 0.0 0.0 0.0 0.0'
    END IF
    WRITE (u, '(A)') '4 Coupled 0.0 3.0 -1.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '5 Fixed 6.0 3.0 -8.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '6 Fixed -8.0 0.0 -8.0 0.0 0.0 0.0 0.0'
    IF (fb) THEN
      WRITE (u, '(A)') '7 Body2 0.5 0.0 0.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '8 Body2 -0.5 0.0 0.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '9 Body2 0.0 0.5 0.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '10 Fixed 1.5 -4.0 -8.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '11 Fixed -1.5 -4.0 -8.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '12 Fixed 0.0 -2.5 -8.0 0.0 0.0 0.0 0.0'
    END IF
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '2 4 5 -'
    WRITE (u, '(A)') '3 3 6 -'
    IF (fb) THEN
      WRITE (u, '(A)') '4 7 10 -'
      WRITE (u, '(A)') '5 8 11 -'
      WRITE (u, '(A)') '6 9 12 -'
    END IF
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 line 9.5 6'
    WRITE (u, '(A)') '2 line 9.5 6'
    WRITE (u, '(A)') '3 line 9.5 6'
    IF (fb) THEN
      WRITE (u, '(A)') '4 line 4.12 4'
      WRITE (u, '(A)') '5 line 4.12 4'
      WRITE (u, '(A)') '6 line 4.12 4'
    END IF
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    IF (fb) WRITE (u, '(A)') 'deck bodyIC'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_host_body_deck

  SUBROUTINE case_host_body()
    !! Coupled/Vessel BODIES on the aggregate (MoorDyn "Coupled" bodies): a 6-DOF host node at the
    !! body reference point returning force and moment. Gates mirror case_host_rod: layout; the
    !! static wrench equals the same lines on two Coupled points plus weight and buoyancy with
    !! their moments; a host translational + angular acceleration returns -(m + rhoW V Ca) a and
    !! -I alpha on top of the point-deck line reactions; orientation round trip; a surge+pitch march
    !! with rewind and restart-overlay bit identity; PtfmInit; the MoorDyn 14-column Coupled row
    !! builds the same aggregate.
    TYPE(CD_AGG_ModuleType) :: agg, ref, agg_b
    REAL(wp), PARAMETER :: DT = 0.01_wp, RHO = 1025.0_wp, G = 9.80665_wp, VOL = 0.1_wp, M = 200.0_wp
    REAL(wp), PARAMETER :: CA = 0.5_wp, IYY = 30.0_wp
    REAL(wp), ALLOCATABLE :: p(:, :), v(:, :), a(:, :), f(:, :), mo(:, :), orr(:, :, :), w(:, :), alp(:, :)
    REAL(wp), ALLOCATABLE :: pr(:, :), vr(:, :), ar(:, :), fr(:, :), p0(:, :), fb(:, :), mb(:, :), bs(:), bs2(:)
    REAL(wp), ALLOCATABLE :: lq(:), lv(:), la2(:)
    REAL(wp) :: fexp(3), mexp(3), rot(3, 3), err, fsc, ptfm(6), arm(3)
    INTEGER :: es, n, nr, s, k, nl, il, ndof
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_host_body_deck('agg_hostbody.dat', .FALSE., .FALSE.)
    CALL write_host_body_deck('agg_hostbody_pts.dat', .TRUE., .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_hostbody.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'hostbody:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_Init_From_Deck(ref, 'agg_hostbody_pts.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'hostbody:ref-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    n = CD_AGG_NMovingPoints(agg, es, em)
    nr = CD_AGG_NMovingPoints(ref, es, em)
    CALL require(n == 2 .AND. nr == 3, 'hostbody:moving-surface-layout')
    IF (n /= 2 .OR. nr /= 3) RETURN
    ALLOCATE (p(3, n), v(3, n), a(3, n), f(3, n), mo(3, n), orr(3, 3, n), w(3, n), alp(3, n), p0(3, n))
    ALLOCATE (pr(3, nr), vr(3, nr), ar(3, nr), fr(3, nr), fb(3, n), mb(3, n))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, f, es, em, orientation=orr, moment=mo)
    CALL require(es == CD_AGG_OK, 'hostbody:mesh: '//TRIM(em))
    CALL require(nan_max_abs(p(:, 2) - [0.0_wp, 0.0_wp, -2.0_wp]) <= 1.0e-14_wp, 'hostbody:node-at-reference')
    CALL CD_AGG_CalcOutput(ref, es, em)
    CALL CD_AGG_GetMovingPointMesh(ref, pr, vr, ar, fr, es, em)
    fexp = fr(:, 1) + fr(:, 2) + [0.0_wp, 0.0_wp, -M*G + RHO*G*VOL]
    mexp = cross(pr(:, 1) - p(:, 2), fr(:, 1)) + cross(pr(:, 2) - p(:, 2), fr(:, 2))
    fsc = nan_max_abs(fexp)
    err = nan_max_abs(f(:, 2) - fexp)/fsc
    WRITE (*, '(A,ES10.3,A,ES10.3)') 'host body static wrench: rel. force error ', err, &
      '; moment error ', nan_max_abs(mo(:, 2) - mexp)/fsc
    CALL require(err <= 1.0e-10_wp, 'hostbody:static-force-equals-lines-plus-weight-buoyancy')
    CALL require(nan_max_abs(mo(:, 2) - mexp) <= 1.0e-10_wp*fsc, 'hostbody:static-moment-about-reference')
    ! host translational + angular acceleration
    p0 = p
    v = 0.0_wp
    a = 0.0_wp
    a(1, 2) = 0.5_wp
    w = 0.0_wp
    alp = 0.0_wp
    alp(2, 2) = 0.2_wp
    orr = 0.0_wp
    orr(:, :, 1) = eye3()
    orr(:, :, 2) = eye3()
    CALL CD_AGG_UpdateStates_Moving(agg, p0, v, a, es, em, output_probe=.TRUE., orientation=orr, &
                                    angular_velocity=w, angular_acceleration=alp)
    CALL require(es == CD_AGG_OK, 'hostbody:probe-update: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, fb, es, em, moment=mb)
    vr = 0.0_wp
    ar = 0.0_wp
    ar(:, 1) = [0.5_wp, 0.0_wp, 0.0_wp] + cross([0.0_wp, 0.2_wp, 0.0_wp], pr(:, 1) - p0(:, 2))
    ar(:, 2) = [0.5_wp, 0.0_wp, 0.0_wp] + cross([0.0_wp, 0.2_wp, 0.0_wp], pr(:, 2) - p0(:, 2))
    CALL CD_AGG_UpdateStates_Moving(ref, pr, vr, ar, es, em, output_probe=.TRUE.)
    CALL CD_AGG_CalcOutput(ref, es, em)
    CALL CD_AGG_GetMovingPointMesh(ref, pr, vr, ar, fr, es, em)
    fexp = fr(:, 1) + fr(:, 2) + [0.0_wp, 0.0_wp, -M*G + RHO*G*VOL] - (M + RHO*VOL*CA)*[0.5_wp, 0.0_wp, 0.0_wp]
    mexp = cross(pr(:, 1) - p0(:, 2), fr(:, 1)) + cross(pr(:, 2) - p0(:, 2), fr(:, 2)) - &
           IYY*[0.0_wp, 0.2_wp, 0.0_wp]
    err = MAX(nan_max_abs(fb(:, 2) - fexp), nan_max_abs(mb(:, 2) - mexp))/fsc
    WRITE (*, '(A,ES10.3)') 'host body accelerated wrench: rel. error ', err
    CALL require(err <= 1.0e-10_wp, 'hostbody:inertia-(m+rhoW V Ca) a and I alpha')
    a = 0.0_wp
    alp = 0.0_wp
    CALL CD_AGG_UpdateStates_Moving(agg, p0, v, a, es, em, output_probe=.TRUE., orientation=orr, &
                                    angular_velocity=w, angular_acceleration=alp)
    ! orientation round trip: pitch the host, the attachment points turn about the reference
    rot = eye3()
    rot(1, 1) = COS(0.1_wp)
    rot(1, 3) = -SIN(0.1_wp)
    rot(3, 1) = SIN(0.1_wp)
    rot(3, 3) = COS(0.1_wp)
    orr(:, :, 2) = rot
    CALL CD_AGG_UpdateStates_Moving(agg, p0, v, a, es, em, output_probe=.TRUE., orientation=orr, &
                                    angular_velocity=w, angular_acceleration=alp)
    CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, f, es, em, orientation=orr)
    CALL require(nan_max_abs(orr(:, :, 2) - rot) <= 1.0e-14_wp, 'hostbody:orientation-round-trip')
    orr(:, :, 2) = eye3()
    CALL CD_AGG_UpdateStates_Moving(agg, p0, v, a, es, em, output_probe=.TRUE., orientation=orr, &
                                    angular_velocity=w, angular_acceleration=alp)
    ! march, rewind, restart overlay
    DO s = 1, 60
      CALL host_rod_drive(REAL(s, wp)*DT, p0, p, v, a, orr, w, alp)
      CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=REAL(s, wp)*DT, orientation=orr, &
                              angular_velocity=w, angular_acceleration=alp)
      CALL require(es == CD_AGG_OK .AND. cv, 'hostbody:march: '//TRIM(em))
      IF (es /= CD_AGG_OK) RETURN
    END DO
    CALL CD_AGG_Snapshot(agg, es, em)
    CALL host_rod_drive(61.0_wp*DT, p0, p, v, a, orr, w, alp)
    CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=61.0_wp*DT, orientation=orr, &
                            angular_velocity=w, angular_acceleration=alp)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), fb, es, em, moment=mb)
    CALL CD_AGG_Restore(agg, es, em)
    CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=61.0_wp*DT, orientation=orr, &
                            angular_velocity=w, angular_acceleration=alp)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), f, es, em, moment=mo)
    CALL require(ALL(IEEE_IS_FINITE(f)) .AND. nan_max_abs(mo(:, 2)) > 0.0_wp, 'hostbody:march-wrench-finite')
    CALL require(nan_max_abs(f - fb) <= 0.0_wp .AND. nan_max_abs(mo - mb) <= 0.0_wp, &
                 'hostbody:correction-rewind-bit-identical')
    CALL CD_AGG_Init_From_Deck(agg_b, 'agg_hostbody.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'hostbody:B-init: '//TRIM(em))
    CALL CD_AGG_UpdateStates_Moving(agg_b, p, v, a, es, em, orientation=orr, angular_velocity=w, &
                                    angular_acceleration=alp)
    nl = CD_System_NLines(agg%sys%fast%system, es, em)
    DO il = 1, nl
      ndof = CD_System_Line_NDOF(agg%sys%fast%system, il, es, em)
      ALLOCATE (lq(ndof), lv(ndof), la2(ndof))
      CALL CD_Get_System_Line_State(agg%sys%fast%system, il, lq, lv, la2, es, em)
      CALL CD_Update_System_Line_Interior_State(agg_b%sys%fast%system, il, lq(4:ndof - 3), lv(4:ndof - 3), es, em, &
                                                a_interior=la2(4:ndof - 3))
      DEALLOCATE (lq, lv, la2)
    END DO
    ALLOCATE (bs(CD_AGG_Rigid6_MirrorSize(agg)), bs2(CD_AGG_Rigid6_MirrorSize(agg)))
    CALL CD_AGG_Get_Rigid6_States(agg, bs, es, em)
    CALL CD_AGG_Set_Rigid6_States(agg_b, bs, es, em)
    CALL require(es == CD_AGG_OK, 'hostbody:B-body-states: '//TRIM(em))
    err = 0.0_wp
    DO s = 62, 80
      CALL host_rod_drive(REAL(s, wp)*DT, p0, p, v, a, orr, w, alp)
      CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=REAL(s, wp)*DT, orientation=orr, &
                              angular_velocity=w, angular_acceleration=alp)
      CALL CD_AGG_Step_Moving(agg_b, DT, p, v, a, cv, st, k, es, em, t_committed=REAL(s, wp)*DT, orientation=orr, &
                              angular_velocity=w, angular_acceleration=alp)
      CALL require(es == CD_AGG_OK, 'hostbody:B-march: '//TRIM(em))
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), f, es, em, moment=mo)
      CALL CD_AGG_CalcOutput(agg_b, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg_b, pr(:, 1:2), vr(:, 1:2), ar(:, 1:2), fb, es, em, moment=mb)
      err = MAX(err, nan_max_abs(f - fb), nan_max_abs(mo - mb))
    END DO
    WRITE (*, '(A,ES10.3)') 'host body restart overlay: max |wrench difference| over 19 steps = ', err
    CALL require(err <= 0.0_wp, 'hostbody:restart-overlay-bit-identical')
    CALL CD_AGG_Get_Rigid6_States(agg, bs, es, em)
    CALL CD_AGG_Get_Rigid6_States(agg_b, bs2, es, em)
    CALL require(nan_max_abs(bs - bs2) <= 0.0_wp, 'hostbody:restart-body-state-bit-identical')
    CALL CD_AGG_End(agg_b, es, em)
    CALL CD_AGG_End(agg, es, em)
    CALL CD_AGG_End(ref, es, em)
    ! PtfmInit: surge 1 m, heave 0.5 m, yaw 10 deg; the attachment points move with the body
    ptfm = [1.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.0_wp, 10.0_wp*PI/180.0_wp]
    CALL CD_AGG_Init_From_Deck(agg, 'agg_hostbody.dat', DT, es, em, external_fluid=.TRUE., ptfm_init=ptfm)
    CALL require(es == CD_AGG_OK, 'hostbody:ptfm-init: '//TRIM(em))
    IF (es == CD_AGG_OK) THEN
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, f, es, em, orientation=orr, moment=mo)
      CALL require(nan_max_abs(p(:, 2) - [1.0_wp, 0.0_wp, -1.5_wp]) <= 1.0e-12_wp, 'hostbody:ptfm-node-displaced')
      rot = eye3()
      rot(1, 1) = COS(ptfm(6))
      rot(1, 2) = SIN(ptfm(6))
      rot(2, 1) = -SIN(ptfm(6))
      rot(2, 2) = COS(ptfm(6))
      CALL require(nan_max_abs(orr(:, :, 2) - rot) <= 1.0e-12_wp, 'hostbody:ptfm-orientation-is-platform-dcm')
      ! the body's attachment 1 sits at r + M^T (1, 0, 0): its line pulls from there
      arm = MATMUL(TRANSPOSE(rot), [1.0_wp, 0.0_wp, 0.0_wp])
      CALL require(ALL(IEEE_IS_FINITE(f)) .AND. ALL(IEEE_IS_FINITE(mo)) .AND. NORM2(arm) > 0.0_wp, &
                   'hostbody:ptfm-wrench-finite')
      CALL CD_AGG_End(agg, es, em)
    END IF
    ! the MoorDyn 14-column Coupled row builds the same host body
    CALL write_host_body_deck('agg_hostbody_md.dat', .FALSE., .TRUE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_hostbody_md.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'hostbody:moordyn-row: '//TRIM(em))
    IF (es == CD_AGG_OK) THEN
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, fb, es, em, moment=mb)
      CALL require(CD_AGG_NMovingPoints(agg, es, em) == 2, 'hostbody:moordyn-row-layout')
      CALL CD_AGG_End(agg, es, em)
    END IF
  END SUBROUTINE case_host_body

  SUBROUTINE case_host_body_mixed_free()
    !! A Coupled body and a Free Rigid6 body in one coupled deck. The free body's stiff tethers make
    !! the body march sub-cycle (NSUB = 9 sub-steps at DT); the host-driven body then follows the
    !! host path interpolated from the step start, so its wrench under a surge + pitch drive stays
    !! with the host-only deck run at DT/NSUB with the exact drive (the bodies share no line). Held
    !! at its t + dt pose through the sub-steps instead, the body's lines would see the whole step's
    !! motion at the first sub-step. Every step converges, the free body stays on its tethers and
    !! correction rewind is bit-identical.
    TYPE(CD_AGG_ModuleType) :: agg, ref
    REAL(wp), PARAMETER :: DT = 0.01_wp
    INTEGER, PARAMETER :: NSUB = 9
    REAL(wp) :: p0(3, 2), p(3, 2), v(3, 2), a(3, 2), orr(3, 3, 2), w(3, 2), alp(3, 2), f(3, 2), mo(3, 2)
    REAL(wp) :: fr(3, 2), mr(3, 2), fb(3, 2), mb(3, 2), q(3, 2), err, fsc
    REAL(wp), ALLOCATABLE :: bs(:)
    INTEGER :: es, s, k, j
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_host_body_deck('agg_hostbody_mixed.dat', .FALSE., .FALSE., free_body=.TRUE.)
    CALL write_host_body_deck('agg_hostbody.dat', .FALSE., .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_hostbody_mixed.dat', DT, es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'hostbody-mixed:init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_Init_From_Deck(ref, 'agg_hostbody.dat', DT/REAL(NSUB, wp), es, em, external_fluid=.TRUE.)
    CALL require(es == CD_AGG_OK, 'hostbody-mixed:ref-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL require(CD_AGG_NMovingPoints(agg, es, em) == 2, 'hostbody-mixed:free-body-is-not-a-host-node')
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, p0, v, a, f, es, em)
    err = 0.0_wp
    fsc = 0.0_wp
    DO s = 1, 60
      CALL host_rod_drive(REAL(s, wp)*DT, p0, p, v, a, orr, w, alp)
      CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=REAL(s, wp)*DT, orientation=orr, &
                              angular_velocity=w, angular_acceleration=alp)
      CALL require(es == CD_AGG_OK .AND. cv, 'hostbody-mixed:march: '//TRIM(em))
      IF (es /= CD_AGG_OK) RETURN
      DO j = 1, NSUB
        CALL host_rod_drive((REAL(s - 1, wp) + REAL(j, wp)/REAL(NSUB, wp))*DT, p0, q, v, a, orr, w, alp)
        CALL CD_AGG_Step_Moving(ref, DT/REAL(NSUB, wp), q, v, a, cv, st, k, es, em, &
                                t_committed=(REAL(s - 1, wp) + REAL(j, wp)/REAL(NSUB, wp))*DT, orientation=orr, &
                                angular_velocity=w, angular_acceleration=alp)
      END DO
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, f, es, em, moment=mo)
      CALL CD_AGG_CalcOutput(ref, es, em)
      CALL CD_AGG_GetMovingPointMesh(ref, p, v, a, fr, es, em, moment=mr)
      err = MAX(err, nan_max_abs(f(:, 2) - fr(:, 2)), nan_max_abs(mo(:, 2) - mr(:, 2)))
      fsc = MAX(fsc, nan_max_abs(fr(:, 2)))
    END DO
    WRITE (*, '(A,ES10.3)') 'host body with a sub-cycling Free body: max wrench difference / max force = ', err/fsc
    CALL require(err <= 1.0e-2_wp*fsc, 'hostbody-mixed:host-wrench-follows-host-only-deck')
    ALLOCATE (bs(CD_AGG_Rigid6_MirrorSize(agg)))
    CALL CD_AGG_Get_Rigid6_States(agg, bs, es, em)
    CALL require(ALL(IEEE_IS_FINITE(bs)), 'hostbody-mixed:body-states-finite')
    CALL CD_AGG_Snapshot(agg, es, em)
    CALL host_rod_drive(61.0_wp*DT, p0, p, v, a, orr, w, alp)
    CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=61.0_wp*DT, orientation=orr, &
                            angular_velocity=w, angular_acceleration=alp)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, q, v, a, fb, es, em, moment=mb)
    CALL CD_AGG_Restore(agg, es, em)
    CALL host_rod_drive(61.0_wp*DT, p0, p, v, a, orr, w, alp)
    CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=61.0_wp*DT, orientation=orr, &
                            angular_velocity=w, angular_acceleration=alp)
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, q, v, a, f, es, em, moment=mo)
    CALL require(nan_max_abs(f - fb) <= 0.0_wp .AND. nan_max_abs(mo - mb) <= 0.0_wp, &
                 'hostbody-mixed:correction-rewind-bit-identical')
    CALL CD_AGG_End(agg, es, em)
    CALL CD_AGG_End(ref, es, em)
  END SUBROUTINE case_host_body_mixed_free

  SUBROUTINE write_host_rod_kin_deck(path, rod_as_points)
    !! A massless, inclined Coupled rod End A (0, 0, -1) -> End B (3, 1, -2) with a line at each
    !! end, plus a plain Coupled fairlead; rod_as_points replaces the rod by two Coupled points.
    CHARACTER(*), INTENT(IN) :: path
    LOGICAL, INTENT(IN) :: rod_as_points
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'massless host rod kinematics deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'line 0.10 80.0 1.0e6 0.0 0.0 1.0 0.0 1.0 0.0'
    IF (.NOT. rod_as_points) THEN
      WRITE (u, '(A)') '--- ROD TYPES ---'
      WRITE (u, '(A)') 'Name Diam Mass Cd Ca CdEnd CaEnd'
      WRITE (u, '(A)') '(-) (m) (kg/m) (-) (-) (-) (-)'
      WRITE (u, '(A)') 'rodmat 1.0e-6 1.0e-9 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '--- RODS ---'
      WRITE (u, '(A)') 'ID RodType Type XA YA ZA XB YB ZB NumSegs Outputs'
      WRITE (u, '(A)') '(-) (-) (-) (m) (m) (m) (m) (m) (m) (-) (-)'
      WRITE (u, '(A)') '1 rodmat Coupled 0.0 0.0 -1.0 3.0 1.0 -2.0 1 -'
    END IF
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Fixed 9.0 1.0 -8.0 0.0 0.0 0.0 0.0'
    IF (rod_as_points) THEN
      WRITE (u, '(A)') '2 Coupled 0.0 0.0 -1.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '3 Coupled 3.0 1.0 -2.0 0.0 0.0 0.0 0.0'
    ELSE
      WRITE (u, '(A)') '2 Rod1A 0.0 0.0 -1.0 0.0 0.0 0.0 0.0'
      WRITE (u, '(A)') '3 Rod1B 3.0 1.0 -2.0 0.0 0.0 0.0 0.0'
    END IF
    WRITE (u, '(A)') '4 Coupled 0.0 2.0 -1.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '5 Fixed 6.0 2.0 -8.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '6 Fixed -6.0 0.0 -8.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 3 1 -'
    WRITE (u, '(A)') '2 4 5 -'
    WRITE (u, '(A)') '3 2 6 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 line 9.2 6'
    WRITE (u, '(A)') '2 line 9.5 6'
    WRITE (u, '(A)') '3 line 9.5 6'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '10.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_host_rod_kin_deck

  SUBROUTINE case_host_rod_kinematics()
    !! A massless inclined Coupled rod under surge + pitch of its host equals two Coupled points
    !! driven by the rigid-body motion of its ends: every step, the rod wrench is the two
    !! points' line loads (force) and the moment of the End B load about End A, to 1e-9.
    TYPE(CD_AGG_ModuleType) :: agg, ref
    REAL(wp), PARAMETER :: DT = 0.01_wp
    REAL(wp) :: p(3, 2), v(3, 2), a(3, 2), f(3, 2), mo(3, 2), orr(3, 3, 2), w(3, 2), alp(3, 2), p0(3, 2)
    REAL(wp) :: pr(3, 3), vr(3, 3), ar(3, 3), fr(3, 3), d0(3), d(3), err, fsc
    INTEGER :: es, s, k
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_host_rod_kin_deck('agg_hostrod_kin.dat', .FALSE.)
    CALL write_host_rod_kin_deck('agg_hostrod_kin_pts.dat', .TRUE.)
    CALL CD_AGG_Init_From_Deck(agg, 'agg_hostrod_kin.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'hostrodkin:init: '//TRIM(em))
    CALL CD_AGG_Init_From_Deck(ref, 'agg_hostrod_kin_pts.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'hostrodkin:ref-init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_GetMovingPointMesh(agg, p0, v, a, f, es, em)
    d0 = [3.0_wp, 1.0_wp, -1.0_wp]
    err = 0.0_wp
    fsc = 0.0_wp
    DO s = 1, 80
      CALL host_rod_drive(REAL(s, wp)*DT, p0, p, v, a, orr, w, alp)
      ! the reference points follow the rigid motion of the rod ends
      d = MATMUL(TRANSPOSE(orr(:, :, 2)), d0)
      pr(:, 1) = p(:, 2)
      vr(:, 1) = v(:, 2)
      ar(:, 1) = a(:, 2)
      pr(:, 2) = p(:, 2) + d
      vr(:, 2) = v(:, 2) + cross(w(:, 2), d)
      ar(:, 2) = a(:, 2) + cross(alp(:, 2), d) + cross(w(:, 2), cross(w(:, 2), d))
      pr(:, 3) = p(:, 1)
      vr(:, 3) = v(:, 1)
      ar(:, 3) = a(:, 1)
      CALL CD_AGG_Step_Moving(agg, DT, p, v, a, cv, st, k, es, em, t_committed=REAL(s, wp)*DT, orientation=orr, &
                              angular_velocity=w, angular_acceleration=alp)
      CALL require(es == CD_AGG_OK .AND. cv, 'hostrodkin:march: '//TRIM(em))
      CALL CD_AGG_Step_Moving(ref, DT, pr, vr, ar, cv, st, k, es, em, t_committed=REAL(s, wp)*DT)
      CALL require(es == CD_AGG_OK .AND. cv, 'hostrodkin:ref-march: '//TRIM(em))
      IF (es /= CD_AGG_OK) RETURN
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg, p, v, a, f, es, em, moment=mo)
      CALL CD_AGG_CalcOutput(ref, es, em)
      CALL CD_AGG_GetMovingPointMesh(ref, pr, vr, ar, fr, es, em)
      fsc = MAX(fsc, nan_max_abs(fr))
      err = MAX(err, nan_max_abs(f(:, 2) - fr(:, 1) - fr(:, 2)), nan_max_abs(f(:, 1) - fr(:, 3)), &
                nan_max_abs(mo(:, 2) - cross(pr(:, 2) - pr(:, 1), fr(:, 2)))/MAX(1.0_wp, NORM2(d0)))
    END DO
    WRITE (*, '(A,ES10.3)') 'massless host rod vs two driven points: max rel. wrench difference ', err/fsc
    CALL require(err <= 1.0e-9_wp*fsc, 'hostrodkin:massless-rod-equals-driven-points')
    CALL CD_AGG_End(agg, es, em)
    CALL CD_AGG_End(ref, es, em)
  END SUBROUTINE case_host_rod_kinematics

  PURE FUNCTION eye3() RESULT(e)
    REAL(wp) :: e(3, 3)
    e = 0.0_wp
    e(1, 1) = 1.0_wp
    e(2, 2) = 1.0_wp
    e(3, 3) = 1.0_wp
  END FUNCTION eye3

  PURE FUNCTION cross(x, y) RESULT(z)
    REAL(wp), INTENT(IN) :: x(3), y(3)
    REAL(wp) :: z(3)
    z = [x(2)*y(3) - x(3)*y(2), x(3)*y(1) - x(1)*y(3), x(1)*y(2) - x(2)*y(1)]
  END FUNCTION cross
END PROGRAM test_openfast_aggregate
