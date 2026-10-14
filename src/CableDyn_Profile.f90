! File: src/CableDyn_Profile.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_Profile
  !! Opt-in wall-clock profile of a standalone driver run, requested with the environment
  !! variable CABLEDYN_PROFILE (as for the OpenFAST module).
  !!
  !! The driver arms it once at start. The run then attributes its wall time to phases by
  !! laps: CD_Profile_Lap(phase) adds the time since the previous lap to that phase, so the
  !! phases partition the run exactly. Each time march laps after its fluid update, its step,
  !! its checks and its output, and once per committed step (CD_Profile_Step). CD_Step_Model
  !! counts the EI = 0 line steps and their Newton iterations. Unarmed, every entry point
  !! returns after one test, so a run that does not ask for the profile, and a host that loads
  !! the library, are not affected. The state is process-global and serves one run: the laps
  !! are taken on the driver's thread, and the line-step counters are updated atomically
  !! because lines step concurrently.
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int64
  USE CableDyn_Precision, ONLY: wp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: CD_Profile_Requested, CD_Profile_Arm, CD_Profile_Armed, CD_Profile_Lap, CD_Profile_Step, &
            CD_Profile_Count_Line_Step, CD_Profile_Seconds, CD_Profile_Counts

  INTEGER, PARAMETER, PUBLIC :: CD_PROF_INPUT = 1     ! deck parsing and auxiliary input files
  INTEGER, PARAMETER, PUBLIC :: CD_PROF_INIT = 2      ! model build, static equilibrium, initial outputs
  INTEGER, PARAMETER, PUBLIC :: CD_PROF_FLUID = 3     ! march: deck waves and current sampled by the driver
  INTEGER, PARAMETER, PUBLIC :: CD_PROF_SOLVE = 4     ! march: line, cable and object steps
  INTEGER, PARAMETER, PUBLIC :: CD_PROF_CHECK = 5     ! march: plausibility checks and line failures
  INTEGER, PARAMETER, PUBLIC :: CD_PROF_OUTPUT = 6    ! march: output channels, line files, range graphs
  INTEGER, PARAMETER, PUBLIC :: CD_PROF_PROGRESS = 7  ! march: progress report and bookkeeping
  INTEGER, PARAMETER, PUBLIC :: CD_PROF_FINAL = 8     ! after the march: closing files, range files
  INTEGER, PARAMETER, PUBLIC :: CD_PROF_NPHASE = 8

  LOGICAL, SAVE :: armed = .FALSE.
  INTEGER(int64), SAVE :: rate = 0_int64, last = 0_int64
  INTEGER(int64), SAVE :: ticks(CD_PROF_NPHASE) = 0_int64
  INTEGER(int64), SAVE :: n_march_steps = 0_int64, n_line_steps = 0_int64, n_line_iter = 0_int64
  INTEGER, SAVE :: max_line_iter = 0

CONTAINS

  LOGICAL FUNCTION CD_Profile_Requested() RESULT(requested)
    !! TRUE when CABLEDYN_PROFILE is set to a non-empty value other than 0.
    CHARACTER(32) :: value
    INTEGER :: stat
    value = ''
    CALL GET_ENVIRONMENT_VARIABLE('CABLEDYN_PROFILE', value, STATUS=stat)
    requested = stat == 0 .AND. LEN_TRIM(value) > 0 .AND. TRIM(value) /= '0'
  END FUNCTION CD_Profile_Requested

  SUBROUTINE CD_Profile_Arm()
    !! Zero the phases and counters and start the clock: the first lap measures from here.
    ticks = 0_int64
    n_march_steps = 0_int64
    n_line_steps = 0_int64
    n_line_iter = 0_int64
    max_line_iter = 0
    CALL SYSTEM_CLOCK(last, rate)
    armed = rate > 0_int64
  END SUBROUTINE CD_Profile_Arm

  LOGICAL FUNCTION CD_Profile_Armed() RESULT(is_armed)
    is_armed = armed
  END FUNCTION CD_Profile_Armed

  SUBROUTINE CD_Profile_Lap(phase)
    !! Add the wall time since the previous lap to phase (a CD_PROF_* value).
    INTEGER, INTENT(IN) :: phase
    INTEGER(int64) :: now
    IF (.NOT. armed) RETURN
    IF (phase < 1 .OR. phase > CD_PROF_NPHASE) RETURN
    CALL SYSTEM_CLOCK(now)
    ticks(phase) = ticks(phase) + (now - last)
    last = now
  END SUBROUTINE CD_Profile_Lap

  SUBROUTINE CD_Profile_Step()
    !! Count a committed step of the march and lap its progress and bookkeeping.
    IF (.NOT. armed) RETURN
    n_march_steps = n_march_steps + 1_int64
    CALL CD_Profile_Lap(CD_PROF_PROGRESS)
  END SUBROUTINE CD_Profile_Step

  SUBROUTINE CD_Profile_Count_Line_Step(n_iter)
    !! Count one EI = 0 line step (a nominal step, a recovery substep or a retry) and its
    !! Newton iterations. Lines of a system step in an OpenMP team, so the updates are atomic.
    INTEGER, INTENT(IN) :: n_iter
    INTEGER(int64) :: iterations
    IF (.NOT. armed) RETURN
    iterations = INT(MAX(0, n_iter), int64)
    !$OMP ATOMIC UPDATE
    n_line_steps = n_line_steps + 1_int64
    !$OMP ATOMIC UPDATE
    n_line_iter = n_line_iter + iterations
    !$OMP ATOMIC UPDATE
    max_line_iter = MAX(max_line_iter, n_iter)
  END SUBROUTINE CD_Profile_Count_Line_Step

  REAL(wp) FUNCTION CD_Profile_Seconds(phase) RESULT(seconds)
    !! Wall time attributed to phase so far [s]; zero when the profile is not armed.
    INTEGER, INTENT(IN) :: phase
    seconds = 0.0_wp
    IF (.NOT. armed) RETURN
    IF (phase < 1 .OR. phase > CD_PROF_NPHASE) RETURN
    seconds = REAL(ticks(phase), wp)/REAL(rate, wp)
  END FUNCTION CD_Profile_Seconds

  SUBROUTINE CD_Profile_Counts(march_steps, line_steps, line_iterations, line_max_iterations)
    !! Committed march steps, EI = 0 line steps, their Newton iterations and the most in one step.
    INTEGER(int64), INTENT(OUT) :: march_steps, line_steps, line_iterations
    INTEGER, INTENT(OUT) :: line_max_iterations
    march_steps = n_march_steps
    line_steps = n_line_steps
    line_iterations = n_line_iter
    line_max_iterations = max_line_iter
  END SUBROUTINE CD_Profile_Counts

END MODULE CableDyn_Profile
