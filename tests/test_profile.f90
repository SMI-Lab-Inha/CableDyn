! File: tests/test_profile.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_profile
  !! The run-time profile of the driver (CableDyn_Profile):
  !!   1. Unarmed, it records nothing.
  !!   2. Armed, the phases partition the elapsed time: each mark charges the time since the
  !!      previous one to the phase that was running, so their sum is the time from arming to
  !!      the last flush.
  !!   3. Time spent in a phase that never hands over (a run that stops with an error) stays
  !!      in that phase when the report flushes it, and a phase never entered (the
  !!      finalisation of such a run) stays exactly zero.
  !!   4. The EI = 0 line-step counters are exact when the lines of an OpenMP team update them
  !!      at once.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_Profile
  USE CableDyn_FatalReport, ONLY: CD_Fatal_Thread_Init
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int64
  IMPLICIT NONE
  INTEGER :: nfail

  nfail = 0
  CALL case_unarmed()
  CALL case_partition_and_failure()
  CALL case_concurrent_counts()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'test_profile: ', nfail, ' failure(s)'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'test_profile: all checks passed'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. cond) THEN
      WRITE (*, '(A)') 'MISMATCH ['//label//']'
      nfail = nfail + 1
    END IF
  END SUBROUTINE require

  SUBROUTINE busy(seconds)
    !! Spin for at least the given wall time.
    REAL(wp), INTENT(IN) :: seconds
    INTEGER(int64) :: c0, c, rate
    CALL SYSTEM_CLOCK(c0, rate)
    DO
      CALL SYSTEM_CLOCK(c)
      IF (REAL(c - c0, wp) >= seconds*REAL(rate, wp)) EXIT
    END DO
  END SUBROUTINE busy

  REAL(wp) FUNCTION total_seconds() RESULT(s)
    INTEGER :: k
    s = 0.0_wp
    DO k = 1, CD_PROF_NPHASE
      s = s + CD_Profile_Seconds(k)
    END DO
  END FUNCTION total_seconds

  SUBROUTINE case_unarmed()
    INTEGER(int64) :: march, steps, iterations
    INTEGER :: max_iterations
    CALL require(.NOT. CD_Profile_Armed(), 'unarmed:not-armed-at-start')
    CALL CD_Profile_Enter(CD_PROF_SOLVE)
    CALL busy(0.002_wp)
    CALL CD_Profile_Flush()
    CALL CD_Profile_Step()
    CALL CD_Profile_Count_Line_Step(4)
    CALL CD_Profile_Counts(march, steps, iterations, max_iterations)
    CALL require(total_seconds() <= 0.0_wp, 'unarmed:no-time')
    CALL require(march == 0 .AND. steps == 0 .AND. iterations == 0 .AND. max_iterations == 0, 'unarmed:no-counts')
  END SUBROUTINE case_unarmed

  SUBROUTINE case_partition_and_failure()
    INTEGER(int64) :: c0, c1, rate
    REAL(wp) :: elapsed
    CALL SYSTEM_CLOCK(c0, rate)
    CALL CD_Profile_Arm()
    CALL require(CD_Profile_Armed(), 'phases:armed')
    CALL busy(0.010_wp)                      ! input
    CALL CD_Profile_Enter(CD_PROF_INIT)
    CALL busy(0.020_wp)                      ! initialisation
    CALL CD_Profile_Enter(CD_PROF_SOLVE)
    CALL busy(0.030_wp)                      ! a step that stops with an error: no further mark
    CALL CD_Profile_Flush()                  ! the report of the failed run
    CALL SYSTEM_CLOCK(c1)
    elapsed = REAL(c1 - c0, wp)/REAL(rate, wp)
    CALL require(CD_Profile_Seconds(CD_PROF_INPUT) >= 0.010_wp, 'phases:input')
    CALL require(CD_Profile_Seconds(CD_PROF_INIT) >= 0.020_wp, 'phases:initialisation')
    CALL require(CD_Profile_Seconds(CD_PROF_SOLVE) >= 0.030_wp, 'phases:failed-step-stays-in-its-phase')
    CALL require(CD_Profile_Seconds(CD_PROF_FINAL) <= 0.0_wp, 'phases:finalisation-never-entered-is-zero')
    CALL require(total_seconds() >= 0.060_wp .AND. total_seconds() <= elapsed, 'phases:sum-is-the-elapsed-time')
    ! a second flush adds the new time to the same running phase
    CALL busy(0.005_wp)
    CALL CD_Profile_Flush()
    CALL require(CD_Profile_Seconds(CD_PROF_SOLVE) >= 0.035_wp .AND. CD_Profile_Seconds(CD_PROF_FINAL) <= 0.0_wp, &
                 'phases:flush-keeps-the-running-phase')
    ! re-arming starts a new profile
    CALL CD_Profile_Arm()
    CALL require(total_seconds() <= 0.0_wp, 'phases:re-armed-from-zero')
  END SUBROUTINE case_partition_and_failure

  SUBROUTINE case_concurrent_counts()
    INTEGER, PARAMETER :: N = 20000
    INTEGER(int64) :: march, steps, iterations
    INTEGER :: max_iterations, i
    CALL CD_Profile_Arm()
    !$OMP PARALLEL NUM_THREADS(4) DEFAULT(SHARED) PRIVATE(i)
    CALL CD_Fatal_Thread_Init()
    !$OMP DO SCHEDULE(STATIC)
    DO i = 1, N
      CALL CD_Profile_Count_Line_Step(MOD(i, 7))
    END DO
    !$OMP END DO
    !$OMP END PARALLEL
    CALL CD_Profile_Step()
    CALL CD_Profile_Counts(march, steps, iterations, max_iterations)
    CALL require(march == 1, 'counts:march-steps')
    CALL require(steps == N, 'counts:line-steps-exact')
    CALL require(iterations == SUM([(INT(MOD(i, 7), int64), i=1, N)]), 'counts:newton-iterations-exact')
    CALL require(max_iterations == 6, 'counts:most-iterations')
  END SUBROUTINE case_concurrent_counts

END PROGRAM test_profile
