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
  !! Usage: test_step_allocations <deck.dat> <water_depth_m> <n_steps> [ranges]
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_LONG_LONG
  USE CableDyn_Precision, ONLY: wp
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

  SUBROUTINE check(stage)
    CHARACTER(*), INTENT(IN) :: stage

    IF (es /= 0) THEN
      PRINT '(A)', 'FAIL: '//stage//': '//TRIM(em)
      ERROR STOP 1
    END IF
  END SUBROUTINE check

END PROGRAM test_step_allocations
