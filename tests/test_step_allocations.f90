! File: tests/test_step_allocations.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_step_allocations
  !! Allocation-count gate for the per-step advance of the OpenFAST aggregate route (the
  !! calls openfast.exe makes every step): after warm-up steps size every reusable
  !! workspace, a further march of fluid-field updates and moving-point steps (EI = 0 lines
  !! and a Hermite cable) must make no heap allocation at all. The count comes from
  !! malloc/calloc/realloc wrappers linked into this executable (tests/allocation_counter.c,
  !! GNU ld --wrap); it covers ALLOCATE, run-time-sized automatic arrays, and the array
  !! temporaries a compiler places on the heap. The snapshot is outside the counted window
  !! (its storage is sized once, at the first call). The output calls (loads, mesh, output
  !! channels) are counted and reported but not gated: the channel evaluation still
  !! allocates per call. With a fourth argument "ranges" the deck's range graphs (LINES flag r)
  !! are accumulated after every step inside the gated window, which must stay allocation-free
  !! too: the envelope accumulation, the touchdown evaluation, and the EI = 0 and Hermite line
  !! samples with their FairTen/AnchTen end-force recovery.
  !!
  !! With the single argument "torsion" the gate runs on a finite-EI cable with condensed torsion
  !! instead (CD_HFMF_UpdateStates): a line clamped at both ends and restrained in torsion, its
  !! coupled end swaying while the parent rolls about the line axis and the imposed twist
  !! changes, so every step evaluates the twist kernel, its Hessian and the bordered solve. The
  !! committed-state outputs (CD_HFMF_SetCoupledKinematics, CD_HFMF_CalcOutput with the connection
  !! moment, both end forces, the energy and the torque query) are gated as well.
  !!
  !! Usage: test_step_allocations <deck.dat> <water_depth_m> <n_steps> [ranges]
  !!        test_step_allocations torsion
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_LONG_LONG
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_EndConnection, ONLY: CD_ENDCONN_RIGID
  USE CableDyn_HermiteTorsion, ONLY: CD_HermiteTorsionType
  USE CableDyn_OpenFAST_HermiteFMF, ONLY: CD_HFMF_ModuleType, CD_HFMF_Init, CD_HFMF_Set_EndConnection, &
                                          CD_HFMF_Set_Torsion, CD_HFMF_UpdateStates, CD_HFMF_CalcOutput, CD_HFMF_End, &
                                          CD_HFMF_SetCoupledKinematics
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCable_Dyn_End_Force, CD_HermiteCable_Dyn_Energy, &
                                          CD_HermiteCable_Dyn_Torsion_State
  USE CableDyn_OpenFAST_Aggregate, ONLY: CD_AGG_ModuleType, CD_AGG_Init_From_Deck, CD_AGG_NMovingPoints, &
                                         CD_AGG_NFluidNodes, CD_AGG_NumChannels, CD_AGG_GetMovingPointMesh, &
                                         CD_AGG_Snapshot, CD_AGG_SetFluidFields, CD_AGG_Step_Moving, &
                                         CD_AGG_CalcOutput, CD_AGG_EvalChannel, CD_AGG_End, &
                                         CD_AGG_Range_Sample
  IMPLICIT NONE
  INTERFACE
    FUNCTION cd_test_allocation_count() BIND(C, NAME='cd_test_allocation_count') RESULT(n)
      IMPORT :: C_LONG_LONG
      INTEGER(C_LONG_LONG) :: n
    END FUNCTION cd_test_allocation_count
    SUBROUTINE cd_test_allocation_counting(on) BIND(C, NAME='cd_test_allocation_counting')
      IMPORT :: C_INT
      INTEGER(C_INT), VALUE :: on
    END SUBROUTINE cd_test_allocation_counting
  END INTERFACE
  INTEGER, PARAMETER :: WARMUP_STEPS = 40
  REAL(wp), PARAMETER :: DT = 0.025_wp, SURGE_AMPLITUDE = 0.5_wp, SURGE_PERIOD = 20.0_wp
  TYPE(CD_AGG_ModuleType) :: agg
  CHARACTER(1024) :: deck, arg
  CHARACTER(2048) :: em
  REAL(wp), ALLOCATABLE :: pos0(:, :), pos(:, :), vel(:, :), acc(:, :), loads(:, :), orient(:, :, :)
  REAL(wp), ALLOCATABLE :: fluid_v(:, :), fluid_a(:, :), waterline(:)
  REAL(wp) :: depth, value
  INTEGER :: es, n_steps, ncp, nfluid, nch, i, ios
  INTEGER(C_LONG_LONG) :: n_step_alloc, n_output_alloc
  LOGICAL :: ranges

  IF (COMMAND_ARGUMENT_COUNT() == 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, arg)
    IF (TRIM(arg) /= 'torsion') ERROR STOP 'a single argument must be "torsion"'
    CALL torsion_gate()
    STOP
  END IF
  IF (COMMAND_ARGUMENT_COUNT() /= 3 .AND. COMMAND_ARGUMENT_COUNT() /= 4) THEN
    PRINT '(A)', 'usage: test_step_allocations <deck.dat> <water_depth_m> <n_steps> [ranges]'
    ERROR STOP 2
  END IF
  ranges = COMMAND_ARGUMENT_COUNT() == 4
  IF (ranges) THEN
    CALL GET_COMMAND_ARGUMENT(4, arg)
    IF (TRIM(arg) /= 'ranges') ERROR STOP 'the fourth argument must be "ranges"'
  END IF
  CALL GET_COMMAND_ARGUMENT(1, deck)
  CALL GET_COMMAND_ARGUMENT(2, arg)
  READ (arg, *, IOSTAT=ios) depth
  IF (ios /= 0) ERROR STOP 'invalid water depth'
  CALL GET_COMMAND_ARGUMENT(3, arg)
  READ (arg, *, IOSTAT=ios) n_steps
  IF (ios /= 0 .OR. n_steps <= 0) ERROR STOP 'invalid step count'

  CALL CD_AGG_Init_From_Deck(agg, TRIM(deck), DT, es, em, env_gravity=9.80665_wp, env_rho_water=1025.0_wp, &
                             env_wtrdpth=depth, external_fluid=.TRUE., range_files=ranges)
  CALL check('initialization')
  IF (ranges .AND. .NOT. agg%ranges%active) ERROR STOP 'the deck requests no range graph'
  ncp = CD_AGG_NMovingPoints(agg, es, em)
  CALL check('moving-point count')
  nfluid = CD_AGG_NFluidNodes(agg, es, em)
  CALL check('fluid-node count')
  nch = CD_AGG_NumChannels(agg)
  IF (ncp <= 0) ERROR STOP 'deck has no moving points to drive'
  ALLOCATE (pos0(3, ncp), pos(3, ncp), vel(3, ncp), acc(3, ncp), loads(3, ncp), orient(3, 3, ncp))
  ALLOCATE (fluid_v(3, nfluid), fluid_a(3, nfluid), waterline(nfluid))
  CALL CD_AGG_GetMovingPointMesh(agg, pos0, vel, acc, loads, es, em)
  CALL check('initial mesh')
  orient = 0.0_wp
  DO i = 1, ncp
    orient(1, 1, i) = 1.0_wp
    orient(2, 2, i) = 1.0_wp
    orient(3, 3, i) = 1.0_wp
  END DO
  fluid_v = 0.0_wp
  fluid_a = 0.0_wp
  waterline = 0.0_wp

  n_step_alloc = 0
  n_output_alloc = 0
  CALL march(0, WARMUP_STEPS, .FALSE.)
  CALL march(WARMUP_STEPS, n_steps, .TRUE.)
  CALL CD_AGG_End(agg, es, em)
  CALL check('end')
  PRINT '(A,I0,A,I0,A,I0,A)', TRIM(deck)//': ', n_steps, ' counted steps, ', n_step_alloc, &
    ' heap allocations in the steps, ', n_output_alloc, ' in the outputs'
  IF (n_step_alloc /= 0) THEN
    PRINT '(A)', 'FAIL: the per-step aggregate advance allocated after warm-up'
    ERROR STOP 1
  END IF
  PRINT '(A)', 'PASS: no heap allocation in the per-step advance after warm-up'

CONTAINS

  SUBROUTINE march(first_step, count, counted)
    INTEGER, INTENT(IN) :: first_step, count
    LOGICAL, INTENT(IN) :: counted
    REAL(wp) :: t, w, ramp
    LOGICAL :: converged, stalled
    INTEGER :: k, n_iter, ich
    INTEGER(C_LONG_LONG) :: c0

    w = 2.0_wp*ACOS(-1.0_wp)/SURGE_PERIOD
    DO k = first_step + 1, first_step + count
      t = k*DT
      ramp = MIN(t/SURGE_PERIOD, 1.0_wp)
      pos = pos0
      vel = 0.0_wp
      acc = 0.0_wp
      pos(1, :) = pos0(1, :) + ramp*SURGE_AMPLITUDE*SIN(w*t)
      vel(1, :) = ramp*SURGE_AMPLITUDE*w*COS(w*t)
      acc(1, :) = -ramp*SURGE_AMPLITUDE*w*w*SIN(w*t)
      CALL CD_AGG_Snapshot(agg, es, em)
      CALL check('snapshot')
      c0 = cd_test_allocation_count()
      IF (counted) CALL cd_test_allocation_counting(1_C_INT)
      CALL CD_AGG_SetFluidFields(agg, fluid_v, fluid_a, waterline, es, em)
      CALL check('fluid fields')
      CALL CD_AGG_Step_Moving(agg, DT, pos, vel, acc, converged, stalled, n_iter, es, em, t_committed=t, &
                              orientation=orient)
      IF (ranges .AND. es == 0) CALL CD_AGG_Range_Sample(agg, t, es, em)
      CALL cd_test_allocation_counting(0_C_INT)
      CALL check('step')
      IF (.NOT. converged .OR. stalled) ERROR STOP 'aggregate step did not converge'
      n_step_alloc = n_step_alloc + (cd_test_allocation_count() - c0)
      c0 = cd_test_allocation_count()
      IF (counted) CALL cd_test_allocation_counting(1_C_INT)
      CALL CD_AGG_CalcOutput(agg, es, em)
      IF (es == 0) CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, loads, es, em)
      DO ich = 1, nch
        IF (es == 0) CALL CD_AGG_EvalChannel(agg, ich, value, es, em)
      END DO
      CALL cd_test_allocation_counting(0_C_INT)
      CALL check('output')
      n_output_alloc = n_output_alloc + (cd_test_allocation_count() - c0)
    END DO
  END SUBROUTINE march

  SUBROUTINE torsion_gate()
    !! Zero heap allocation in CD_HFMF_UpdateStates of a cable with condensed torsion.
    INTEGER, PARAMETER :: NE = 20, NN = NE + 1, NDOF = 6*NN, N_WARM = 40, N_COUNT = 200
    REAL(wp), PARAMETER :: LL = 10.0_wp, DTT = 0.01_wp
    TYPE(CD_HFMF_ModuleType) :: cab
    TYPE(CD_HermiteTorsionType) :: tors
    REAL(wp) :: q(NDOF), l0(NE), eye(3, 3), frame(3, 2), t, x(3), v(3), a(3), dcm(3, 3), f(3), m(3), psi
    REAL(wp) :: fe(3), ke, se, th, mt
    INTEGER :: k
    INTEGER(C_LONG_LONG) :: c0, n_alloc, n_out, n_energy
    n_out = 0
    n_energy = 0
    q = 0.0_wp
    DO k = 1, NN
      q(6*k - 5) = LL*REAL(k - 1, wp)/REAL(NE, wp)
      q(6*k - 2) = 1.0_wp
    END DO
    l0 = LL/REAL(NE, wp)
    eye = 0.0_wp
    eye(1, 1) = 1.0_wp
    eye(2, 2) = 1.0_wp
    eye(3, 3) = 1.0_wp
    CALL CD_HFMF_Init(cab, l0, [(1.0e5_wp, k=1, NE)], [(100.0_wp, k=1, NE)], [(1.0_wp, k=1, NE)], &
                      [(0.0_wp, k=1, NE)], q, [1, 2, 3, NDOF - 5, NDOF - 4, NDOF - 3], 0.0_wp, 0.0_wp, 0.8_wp, &
                      DTT, NN, es, em, max_iter=60, tol=1.0e-8_wp)
    CALL check('torsion cable init')
    CALL CD_HFMF_Set_EndConnection(cab, [0.0_wp, 0.0_wp], RESHAPE([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], &
                                                                  [3, 2]), es, em, &
                                   coupled_d0_parent=[1.0_wp, 0.0_wp, 0.0_wp], parent_orientation=eye, &
                                   connection_mode=[CD_ENDCONN_RIGID, CD_ENDCONN_RIGID])
    CALL check('torsion cable end connections')
    tors%active = .TRUE.
    tors%phi = 2.0_wp
    tors%ends = RESHAPE([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
                         0.0_wp, 0.0_wp, 1.0_wp], [3, 4])
    ALLOCATE (tors%gj(NE))
    tors%gj = 80.0_wp
    tors%has_theta = .TRUE.
    frame(:, 1) = [1.0_wp, 0.0_wp, 0.0_wp]
    frame(:, 2) = [0.0_wp, 0.0_wp, 1.0_wp]
    CALL CD_HFMF_Set_Torsion(cab, tors, es, em, coupled_frame_parent=frame, parent_orientation=eye)
    CALL check('torsion cable torsion')
    n_alloc = 0
    DO k = 1, N_WARM + N_COUNT
      t = DTT*REAL(k, wp)
      x = [LL, 0.05_wp*SIN(2.0_wp*t), 0.0_wp]
      v = [0.0_wp, 0.1_wp*COS(2.0_wp*t), 0.0_wp]
      a = [0.0_wp, -0.2_wp*SIN(2.0_wp*t), 0.0_wp]
      psi = 3.0_wp*t
      dcm = eye
      dcm(2, 2) = COS(psi)
      dcm(3, 3) = COS(psi)
      dcm(2, 3) = SIN(psi)
      dcm(3, 2) = -SIN(psi)
      c0 = cd_test_allocation_count()
      IF (k > N_WARM) CALL cd_test_allocation_counting(1_C_INT)
      CALL CD_HFMF_UpdateStates(cab, x, v, a, es, em, u_orientation=dcm, u_angular_velocity=[3.0_wp, 0.0_wp, 0.0_wp], &
                                u_angular_acceleration=[0.0_wp, 0.0_wp, 0.0_wp], u_twist=2.0_wp - psi + 0.3_wp*SIN(t))
      CALL cd_test_allocation_counting(0_C_INT)
      CALL check('torsion cable step')
      n_alloc = n_alloc + (cd_test_allocation_count() - c0)
      ! the outputs of the committed state and the non-stepping boundary write: the loads and
      ! connection moment, the end forces, the energy and the twist and torque queries
      c0 = cd_test_allocation_count()
      IF (k > N_WARM) CALL cd_test_allocation_counting(1_C_INT)
      CALL CD_HFMF_SetCoupledKinematics(cab, x, v, a, es, em, u_orientation=dcm, &
                                        u_angular_velocity=[3.0_wp, 0.0_wp, 0.0_wp], &
                                        u_angular_acceleration=[0.0_wp, 0.0_wp, 0.0_wp])
      IF (es == 0) CALL CD_HFMF_CalcOutput(cab, f, es, em, y_moment=m)
      IF (es == 0) CALL CD_HermiteCable_Dyn_End_Force(cab%line, 1, fe, es, em)
      IF (es == 0) CALL CD_HermiteCable_Dyn_End_Force(cab%line, NN, fe, es, em)
      IF (es == 0) CALL CD_HermiteCable_Dyn_Torsion_State(cab%line, th, mt, es, em)
      CALL cd_test_allocation_counting(0_C_INT)
      CALL check('torsion cable output')
      n_out = n_out + (cd_test_allocation_count() - c0)
      ! the energy diagnostic (not a per-step output; its mass product uses a run-time-sized
      ! temporary): only its torsion part is required to be allocation-free, so the count of
      ! the call may not exceed that of the same call without torsion (one temporary)
      c0 = cd_test_allocation_count()
      IF (k > N_WARM) CALL cd_test_allocation_counting(1_C_INT)
      CALL CD_HermiteCable_Dyn_Energy(cab%line, ke, se, es, em)
      CALL cd_test_allocation_counting(0_C_INT)
      CALL check('torsion cable energy')
      IF (k > N_WARM) n_energy = MAX(n_energy, cd_test_allocation_count() - c0)
    END DO
    CALL CD_HFMF_End(cab)
    PRINT '(A,I0,A,I0,A,I0,A)', 'torsion cable: ', N_COUNT, ' counted steps, ', n_alloc, &
      ' heap allocations in the steps, ', n_out, ' in the outputs'
    PRINT '(A,I0)', 'torsion cable energy diagnostic: heap allocations per call ', n_energy
    IF (n_alloc /= 0 .OR. n_out /= 0 .OR. n_energy > 1) THEN
      PRINT '(A)', 'FAIL: the per-step advance or the outputs of a cable with torsion allocated after warm-up'
      ERROR STOP 1
    END IF
    PRINT '(A)', 'PASS: no heap allocation in the per-step advance or outputs of a cable with torsion after warm-up'
  END SUBROUTINE torsion_gate

  SUBROUTINE check(stage)
    CHARACTER(*), INTENT(IN) :: stage

    IF (es /= 0) THEN
      PRINT '(A)', 'FAIL: '//stage//': '//TRIM(em)
      ERROR STOP 1
    END IF
  END SUBROUTINE check

END PROGRAM test_step_allocations
