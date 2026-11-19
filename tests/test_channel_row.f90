! File: tests/test_channel_row.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_channel_row
  !! The deck channels of one output row evaluated together, as the OpenFAST shell's CalcOutput,
  !! its committed-step record, its linearization and the standalone mixed route evaluate them
  !! (CD_AGG_EvalChannels: the line-wide values a cable's channels read -- element twists, torsion
  !! state, end twists, node torques, section curvature -- once per row), against the same
  !! channels one at a time (CD_AGG_EvalChannel), bit for bit. The deck is coupled: an EI = 0
  !! line and two finite-EI cables restrained in torsion (nodal twist), the second with
  !! torque-tension coupling (GJ|K) and a torsional end spring, every fairlead surging and
  !! rolling with the platform; the channels mix twists, torques, bend moments, curvatures,
  !! tensions, kinematics and point forces of all three lines.
  !!   1. the static state, then every step of a rolling, surging march: a row's values must
  !!      follow the state (a line-wide value kept from the row before, the stale row, fails);
  !!      a second row at the same state repeats the first bit for bit;
  !!   2. a non-stepping transfer of other kinematics (the shell's output probe) changes the
  !!      state without a step: the row follows it.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_OpenFAST_Aggregate, ONLY: CD_AGG_ModuleType, CD_AGG_Init_From_Deck, CD_AGG_NMovingPoints, &
                                         CD_AGG_GetMovingPointMesh, CD_AGG_UpdateStates_Moving, &
                                         CD_AGG_Step_Moving, CD_AGG_CalcOutput, CD_AGG_NumChannels, &
                                         CD_AGG_EvalChannel, CD_AGG_EvalChannels, CD_AGG_End, CD_AGG_OK
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, RHOW = 1025.0_wp
  REAL(wp), PARAMETER :: DT = 0.05_wp, DEG = PI/180.0_wp
  INTEGER, PARAMETER :: NSTEP = 12
  INTEGER :: n_fail = 0

  CALL write_deck('chrow.dat')
  CALL check_rows()
  IF (n_fail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', n_fail, ' channel row check(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: channel rows (one row against one channel at a time, through a march and a '// &
    'probe transfer)'

CONTAINS

  SUBROUTINE require(cond, msg)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: msg
    IF (.NOT. cond) THEN
      n_fail = n_fail + 1
      WRITE (*, '(A,A)') 'FAILED: ', msg
    END IF
  END SUBROUTINE require

  LOGICAL FUNCTION same_bits(a, b) RESULT(same)
    !! Whether two arrays of the same size have the same bit patterns.
    REAL(wp), INTENT(IN) :: a(:), b(:)
    INTEGER :: k
    same = SIZE(a) == SIZE(b)
    IF (.NOT. same) RETURN
    DO k = 1, SIZE(a)
      IF (TRANSFER(a(k), 0_INT64) /= TRANSFER(b(k), 0_INT64)) THEN
        same = .FALSE.
        RETURN
      END IF
    END DO
  END FUNCTION same_bits

  SUBROUTINE write_deck(path)
    !! Three taut 100 m lines from fairleads at x = 0 to anchors at x = 100 m, 40 m deep, the cables
    !! neutrally buoyant: line 1 EI = 0 (twice the displaced mass); line 2 finite EI, clamped and
    !! restrained in torsion at both ends with a quarter turn of pretwist; line 3 finite EI under
    !! torque-tension coupling, clamped, rigid in torsion at End A and on a torsional spring at End B.
    CHARACTER(*), INTENT(IN) :: path
    INTEGER :: u
    REAL(wp) :: m
    CHARACTER(32) :: mtxt, rtxt
    m = RHOW*0.25_wp*PI*0.04_wp
    WRITE (mtxt, '(ES24.16)') m
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'channel row gate'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI GAs GJ Irt Irn Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (N) (Nm2) (kgm) (kgm) (-) (-) (-) (-)'
    WRITE (rtxt, '(ES24.16)') 2.0_wp*m
    WRITE (u, '(A)') 'rope 0.2 '//TRIM(rtxt)//' 1.0e7 0.0 0.0 0.0 0.0 0.0 0.0 1.2 0.0 1.0 0.0'
    WRITE (u, '(A)') 'cab 0.2 '//TRIM(mtxt)//' 1.0e7 0.0 1.0e6 1.0e8 5.0e4 1.0 2.0 1.2 0.0 1.0 0.0'
    WRITE (u, '(A)') 'arm 0.2 '//TRIM(mtxt)//' 1.0e7 0.0 1.0e6 1.0e8 5.0e4|4.0e3 1.0 2.0 1.2 0.0 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Coupled 0.0 -20.0 -40.0'
    WRITE (u, '(A)') '2 Fixed 100.0 -20.0 -40.0'
    WRITE (u, '(A)') '3 Coupled 0.0 0.0 -40.0'
    WRITE (u, '(A)') '4 Fixed 100.0 0.0 -40.0'
    WRITE (u, '(A)') '5 Coupled 0.0 20.0 -40.0'
    WRITE (u, '(A)') '6 Fixed 100.0 20.0 -40.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '2 3 4 -'
    WRITE (u, '(A)') '3 5 6 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 rope 99.5 10'
    WRITE (u, '(A)') '2 cab 99.5 12'
    WRITE (u, '(A)') '3 arm 99.5 12'
    WRITE (u, '(A)') '--- END CONNECTIONS ---'
    WRITE (u, '(A)') 'LineID End Stiffness EzX EzY EzZ TorsStiffness NxX NxY NxZ Pretwist'
    WRITE (u, '(A)') '(-) (-) (N-m/rad) (-) (-) (-) (N-m/rad) (-) (-) (-) (deg)'
    WRITE (u, '(A)') '2 A Rigid 1 0 0 Rigid 0 0 1 0'
    WRITE (u, '(A)') '2 B Rigid 1 0 0 Rigid 0 0 1 90'
    WRITE (u, '(A)') '3 A Rigid 1 0 0 Rigid 0 0 1 0'
    WRITE (u, '(A)') '3 B Rigid 1 0 0 2.0e4 0 0 1 45'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '0.05 dtM'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1'
    WRITE (u, '(A)') 'L1N4pz'
    WRITE (u, '(A)') 'L1N4ax'
    WRITE (u, '(A)') 'Point1Fx'
    WRITE (u, '(A)') 'Twist2'
    WRITE (u, '(A)') 'Twist2N1'
    WRITE (u, '(A)') 'Twist2N5'
    WRITE (u, '(A)') 'Twist2N12'
    WRITE (u, '(A)') 'Torq2N1'
    WRITE (u, '(A)') 'Torq2N7'
    WRITE (u, '(A)') 'BendMom2N3'
    WRITE (u, '(A)') 'BendMom2N8'
    WRITE (u, '(A)') 'Curv2N4'
    WRITE (u, '(A)') 'CurvX2N4'
    WRITE (u, '(A)') 'CurvY2N9'
    WRITE (u, '(A)') 'Ten2N6'
    WRITE (u, '(A)') 'L2N6vy'
    WRITE (u, '(A)') 'L2N6ay'
    WRITE (u, '(A)') 'FairTen2'
    WRITE (u, '(A)') 'Twist3'
    WRITE (u, '(A)') 'Twist3N2'
    WRITE (u, '(A)') 'Twist3N11'
    WRITE (u, '(A)') 'Torq3N1'
    WRITE (u, '(A)') 'Torq3N13'
    WRITE (u, '(A)') 'BendMom3N6'
    WRITE (u, '(A)') 'CurvX3N6'
    WRITE (u, '(A)') 'Ten3N6'
    WRITE (u, '(A)') 'AnchTen3'
    WRITE (u, '(A)') 'Point5Fz'
    WRITE (u, '(A)') '--- END ---'
    CLOSE (u)
  END SUBROUTINE write_deck

  FUNCTION roll_dcm(phi) RESULT(r)
    !! OpenFAST mesh orientation (global to local) of a roll phi about +x.
    REAL(wp), INTENT(IN) :: phi
    REAL(wp) :: r(3, 3)
    r = 0.0_wp
    r(1, 1) = 1.0_wp
    r(2, 2) = COS(phi)
    r(2, 3) = SIN(phi)
    r(3, 2) = -SIN(phi)
    r(3, 3) = COS(phi)
  END FUNCTION roll_dcm

  SUBROUTINE one_at_a_time(agg, values)
    TYPE(CD_AGG_ModuleType), INTENT(INOUT) :: agg
    REAL(wp), INTENT(OUT) :: values(:)
    INTEGER :: i, es
    CHARACTER(512) :: em
    DO i = 1, SIZE(values)
      CALL CD_AGG_EvalChannel(agg, i, values(i), es, em)
      CALL require(es == CD_AGG_OK, 'channel evaluates: '//TRIM(em))
    END DO
  END SUBROUTINE one_at_a_time

  SUBROUTINE row(agg, values)
    TYPE(CD_AGG_ModuleType), INTENT(INOUT) :: agg
    REAL(wp), INTENT(OUT) :: values(:)
    INTEGER :: es, n_done
    CHARACTER(512) :: em
    values = 0.0_wp
    CALL CD_AGG_EvalChannels(agg, values, n_done, es, em)
    CALL require(es == CD_AGG_OK .AND. n_done == SIZE(values), 'the row evaluates: '//TRIM(em))
  END SUBROUTINE row

  SUBROUTINE check_rows()
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER, PARAMETER :: NCP = 3
    ! the twist and torque channels away from End A (Twist2N1, the twist from End A to End A, is 0)
    INTEGER, PARAMETER :: TWIST_CH(10) = [5, 7, 8, 9, 10, 20, 21, 22, 23, 24]
    REAL(wp) :: pos(3, NCP), vel(3, NCP), acc(3, NCP), load(3, NCP), orient(3, 3, NCP)
    REAL(wp) :: omega(3, NCP), alpha(3, NCP), x0(3, NCP), pos2(3, NCP)
    REAL(wp), ALLOCATABLE :: v_row(:), v_one(:), v_again(:), v_prev(:), v_b(:)
    REAL(wp) :: phi, rate, amp, w, t
    INTEGER :: es, ni, k, nmp, nch, j
    LOGICAL :: cv, st, rows_ok, moved
    CHARACTER(512) :: em

    CALL CD_AGG_Init_From_Deck(agg, 'chrow.dat', DT, es, em, host_orientation=.TRUE.)
    CALL require(es == CD_AGG_OK, 'the channel row deck builds on the aggregate: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    nmp = CD_AGG_NMovingPoints(agg, es, em)
    CALL require(nmp == NCP .AND. agg%ncable == 2, 'the deck is one EI = 0 line and two coupled cables')
    IF (nmp /= NCP) RETURN
    nch = CD_AGG_NumChannels(agg)
    CALL require(nch == 29, 'the deck declares its 29 channels')
    ALLOCATE (v_row(nch), v_one(nch), v_again(nch), v_prev(nch), v_b(nch))

    ! 1. the static state, then a rolling, surging march
    CALL row(agg, v_row)
    CALL one_at_a_time(agg, v_one)
    CALL require(same_bits(v_row, v_one), '1: static state: the row equals the channels one at a time')
    CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, load, es, em, orientation=orient)
    x0 = pos
    rate = 20.0_wp*DEG
    amp = 0.5_wp
    w = 2.0_wp*PI/1.2_wp
    rows_ok = .TRUE.
    moved = .TRUE.
    v_prev = v_row
    DO k = 1, NSTEP
      t = DT*k
      phi = rate*t
      DO j = 1, NCP
        orient(:, :, j) = roll_dcm(phi)
        pos(:, j) = x0(:, j)
        pos(1, j) = x0(1, j) + amp*SIN(w*t)
        vel(:, j) = 0.0_wp
        vel(1, j) = amp*w*COS(w*t)
        acc(:, j) = 0.0_wp
        acc(1, j) = -amp*w*w*SIN(w*t)
        omega(:, j) = [rate, 0.0_wp, 0.0_wp]
        alpha(:, j) = 0.0_wp
      END DO
      CALL CD_AGG_Step_Moving(agg, DT, pos, vel, acc, cv, st, ni, es, em, orientation=orient, &
                              angular_velocity=omega, angular_acceleration=alpha)
      CALL require(es == CD_AGG_OK .AND. cv, '1: the march step converges: '//TRIM(em))
      IF (es /= CD_AGG_OK) EXIT
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL row(agg, v_row)
      CALL one_at_a_time(agg, v_one)
      CALL row(agg, v_again)
      rows_ok = rows_ok .AND. same_bits(v_row, v_one) .AND. same_bits(v_again, v_row)
      ! every twist and torque channel moves with the roll: a row value kept from the state before
      ! would show here
      moved = moved .AND. ALL(ABS(v_row(TWIST_CH) - v_prev(TWIST_CH)) > 0.0_wp)
      v_prev = v_row
    END DO
    CALL require(rows_ok, '1: every step: the row equals the channels one at a time, and a second row '// &
                 'at the same state repeats it, bit for bit')
    CALL require(moved, '1: the twist and torque channels change at every step (the rows follow the state)')
    WRITE (*, '(A,F10.4,A,F10.4,A)') '1 march: Twist2 ', v_row(5), ' deg, Twist3 ', v_row(20), ' deg'

    ! 2. a non-stepping transfer of other kinematics changes the state without a step
    pos2 = pos
    pos2(2, :) = pos(2, :) + 0.05_wp
    CALL CD_AGG_UpdateStates_Moving(agg, pos2, vel, acc, es, em, output_probe=.TRUE., orientation=orient, &
                                    angular_velocity=omega, angular_acceleration=alpha)
    CALL require(es == CD_AGG_OK, '2: the probe transfer applies: '//TRIM(em))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL row(agg, v_b)
    CALL one_at_a_time(agg, v_one)
    CALL require(same_bits(v_b, v_one), '2: after a probe transfer: the row equals the channels one at a time')
    CALL require(.NOT. same_bits(v_b, v_row), '2: the probe transfer changes the row (it follows the state)')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE check_rows
END PROGRAM test_channel_row
