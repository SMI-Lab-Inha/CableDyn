! File: app/cabledyn.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM cabledyn
  !! The standalone CableDyn driver (CableDyn_driver): solve a MoorDyn-style `.dat` deck
  !! (doc/driver_format.md) and write a tabular `.out`.
  !!
  !!   CableDyn_driver <deck.dat> <out_root>
  !!
  !! The deck follows the OpenFAST-style ecosystem (sectioned text, classic
  !! vocabulary) with the OrcaFlex line model (a line is End A -> End B, built from
  !! ordered SECTIONS, each its own line type + mesh). The driver solves the static
  !! initial condition of every deck it accepts -- EI=0 and finite-EI lines, points,
  !! BODIES and RODS -- and, when dtM/TMax are present, marches the dynamics while
  !! writing the requested OUTPUTS at each row (doc/cli.rst).
  !!
  !! Exit codes: 0 converged, 1 bad/unparseable deck, 2 solve failure. Every non-zero exit
  !! the driver makes itself ends stderr with the closing line "CableDyn_driver: ended with
  !! exit code <n>" (stop_run); an interrupt or a fatal fault is reported by
  !! CableDyn_FatalReport with the simulated time reached.
  USE CableDyn_Precision, ONLY: wp, CD_ZERO
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_Deck_Query_dtM, &
                                 CD_Classify_Dynamic_Completion, CD_DECKDRV_OK, &
                                 CD_DECKDRV_BADINPUT, CD_DECKDRV_SOLVEFAIL, CD_DECK_PATHLEN
  USE CableDyn_PathIO, ONLY: CD_Get_Argument, CD_Native_Path, CD_Lock_Path, CD_PATH_OK, &
                             CD_LOCK_ACQUIRED, CD_LOCK_HELD
  USE CableDyn_OpenFAST_Aggregate, ONLY: CD_AGG_ModuleType, CD_AGG_Init_From_Deck, &
                                         CD_AGG_NMovingPoints, CD_AGG_GetMovingPointMesh, &
                                         CD_AGG_Step_Moving, CD_AGG_CalcOutput, CD_AGG_End, &
                                         CD_AGG_NumChannels, CD_AGG_ChannelHeader, &
                                         CD_AGG_EvalChannel, CD_AGG_WriteStaticProfile, &
                                         CD_AGG_GetInitMetadata, CD_AGG_GetInitLine, &
                                         CD_AGG_Range_Sample, CD_AGG_Range_Write, CD_AGG_OK, CD_AGG_BADINPUT
  USE CableDyn_Banner, ONLY: CD_Print_Banner
  USE CableDyn_Linalg, ONLY: CD_Blas_Runtime_Check, CD_LINALG_OK
  USE CableDyn_FatalReport, ONLY: CD_Fatal_Report_Install, CD_Fatal_Report_Time
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: error_unit, output_unit, int64
  IMPLICIT NONE
  TYPE :: StandaloneProgress
    INTEGER :: nstep = 0, stride = 1
    INTEGER(int64) :: start_count = 0_int64, count_rate = 0_int64
  END TYPE StandaloneProgress
  ! CD_AGG_Init_From_Deck owns dynamic workspaces and therefore requires a positive
  ! clock even when the executable requests only equilibrium. This internal value is
  ! never marched and does not enter the static residual or output equations.
  REAL(wp), PARAMETER :: MIXED_STATIC_INIT_DT = 1.0_wp
  CHARACTER(*), PARAMETER :: PROG = 'CableDyn_driver'
  ! Lock file beside the outputs that keeps a second run off the same output root.
  CHARACTER(*), PARAMETER :: LOCK_SUFFIX = '.cabledyn.lock'
  ! Longest suffix a run appends to the output root (".Line<id>.range.out" with a 10-digit id).
  INTEGER, PARAMETER :: MAX_OUTPUT_SUFFIX = 25
  ! Command-line paths as UTF-8; the native spellings are what the Fortran runtime opens.
  CHARACTER(4096) :: deck_path, out_root
  CHARACTER(:), ALLOCATABLE :: deck_native, root_native
  CHARACTER(CD_DECK_PATHLEN), ALLOCATABLE :: input_files(:)
  INTEGER :: nargs, ErrStat, n_ei0, n_finite_ei, iarg
  REAL(wp) :: dtM, tmax
  LOGICAL :: converged, has_dtm, has_tmax, has_motion_file, has_line_outputs, has_objects
  ! Sized for the aggregate and deck-driver messages it relays (512-character buffers
  ! plus deck-line context), so no diagnostic is cut short.
  CHARACTER(1024) :: ErrMsg

  nargs = COMMAND_ARGUMENT_COUNT()

  ! Version / help flags: print the identity banner (and, for help, the usage) and
  ! exit 0 -- the way MoorDyn's driver announces itself when run with no real work.
  IF (nargs >= 1) THEN
    CALL get_argument(1, deck_path)
    SELECT CASE (TRIM(deck_path))
    CASE ('-v', '-V', '-version', '-VERSION', '--version')
      CALL CD_Print_Banner()
      STOP 0, QUIET = .TRUE.
    CASE ('-h', '-H', '-help', '-HELP', '--help', '-?', '/?')
      CALL CD_Print_Banner()
      CALL print_usage(output_unit)
      STOP 0, QUIET = .TRUE.
    END SELECT
  END IF
  ! Any other dash-prefixed argument is an unknown option, never a deck or output path.
  DO iarg = 1, nargs
    CALL get_argument(iarg, deck_path)
    IF (deck_path(1:1) == '-') THEN
      WRITE (error_unit, '(A)') PROG//': unknown option "'//TRIM(deck_path)//'"'
      CALL print_usage(error_unit)
      CALL stop_run(1)
    END IF
  END DO

  IF (nargs /= 2) THEN
    CALL CD_Print_Banner(error_unit)
    CALL print_usage(error_unit)
    CALL stop_run(1)
  END IF
  CALL get_argument(1, deck_path)
  CALL get_argument(2, out_root)
  CALL prepare_paths(TRIM(deck_path), TRIM(out_root), deck_native, root_native)
  CALL check_output_root(TRIM(deck_path), deck_native, TRIM(out_root), root_native)

  ! Startup identity to stderr. Initialization/progress records and the completion line
  ! remain visible on stdout; automation must use the exit status rather than assuming
  ! stdout contains only one machine-readable record.
  CALL CD_Print_Banner(error_unit)
  ! The exit contract, stated once at start: a caller that finds this line and later no
  ! closing line knows the process was ended from outside (a driver older than this line
  ! writes neither, so its failures must not be read that way).
  WRITE (error_unit, '(A)') '  Exit status: every failure ends stderr with "'//PROG// &
    ': ended with exit code <n>".'
  ! The LAPACK runtime is loaded on first use in the Windows GNU build; a missing or
  ! incompatible library stops here with its diagnostic instead of as a solver failure.
  CALL CD_Blas_Runtime_Check(PROG, ErrStat, ErrMsg)
  IF (ErrStat /= CD_LINALG_OK) THEN
    WRITE (error_unit, '(A)') TRIM(ErrMsg)
    CALL stop_run(2)
  END IF
  ! From here on an interrupt or a fatal fault (an access violation, a stack overflow) is
  ! reported on stderr with the simulated time reached, instead of ending the process
  ! without a word. Installed after the LAPACK runtime has loaded, so its handlers are the
  ! process's own.
  CALL CD_Fatal_Report_Install()
  CALL CD_Deck_Query_dtM(TRIM(deck_path), dtM, has_dtm, ErrStat, ErrMsg, tmax=tmax, &
                         has_tmax=has_tmax, n_ei0=n_ei0, n_finite_ei=n_finite_ei, &
                         has_motion_file=has_motion_file, standalone_scan=.TRUE., &
                         has_line_outputs=has_line_outputs, input_files=input_files, has_objects=has_objects)
  IF (ErrStat == CD_DECKDRV_OK) THEN
    DO iarg = 1, SIZE(input_files)
      CALL check_input_not_output(TRIM(input_files(iarg)), TRIM(out_root), root_native, 'input file')
    END DO
    ! Taken before the first output file is written and held until the process ends.
    CALL lock_output_root(TRIM(out_root))
    ! EI=0 and finite-EI lines between held/moving points run on the coupled aggregate; with
    ! bodies, rods or Free/Connect points the deck driver's multibody march owns the run
    IF (n_ei0 > 0 .AND. n_finite_ei > 0 .AND. .NOT. has_objects) THEN
      CALL run_mixed_held_deck(TRIM(deck_path), root_native, dtM, has_dtm, tmax, has_tmax, &
                               has_motion_file, has_line_outputs, converged, ErrStat, ErrMsg)
    ELSE
      CALL CD_Run_Deck_Driver(TRIM(deck_path), root_native, converged, ErrStat, ErrMsg)
    END IF
  END IF
  IF (ErrStat /= CD_DECKDRV_OK) THEN
    WRITE (error_unit, '(A)') TRIM(ErrMsg)
    IF (ErrStat == CD_DECKDRV_BADINPUT) THEN
      CALL stop_run(1)
    ELSE
      CALL stop_run(2)
    END IF
  END IF
  ! Keep the process-level contract fail-closed even if a future driver route
  ! accidentally returns a non-converged flag without a matching error status.
  IF (.NOT. converged) THEN
    IF (LEN_TRIM(ErrMsg) > 0) THEN
      WRITE (error_unit, '(A)') TRIM(ErrMsg)
    ELSE
      WRITE (error_unit, '(A)') 'CableDyn_driver: run did not converge; output is for inspection only.'
    END IF
    CALL stop_run(2)
  END IF
  WRITE (output_unit, '(A,A,A)') PROG//': converged run written to ', TRIM(out_root), '.out'

CONTAINS

  SUBROUTINE stop_run(code)
    !! End the process with a non-zero exit code, after the diagnostic already written.
    !! The closing line is the last stderr record of every exit the driver makes itself, so
    !! a caller can tell such an exit from a process ended from outside (Task Manager,
    !! `taskkill /F`, `kill -9`), which runs no code of the driver and leaves no closing line.
    INTEGER, INTENT(IN) :: code
    WRITE (error_unit, '(A,I0)') PROG//': ended with exit code ', code
    FLUSH (error_unit)
    STOP code, QUIET = .TRUE.
  END SUBROUTINE stop_run

  SUBROUTINE get_argument(index, value)
    !! Command-line argument that fails closed when it cannot be retrieved whole
    !! (a path longer than the buffer would otherwise be silently truncated).
    !! The argument is UTF-8 (on Windows, taken from the wide command line).
    INTEGER, INTENT(IN) :: index
    CHARACTER(*), INTENT(OUT) :: value
    INTEGER :: arg_len, arg_stat
    CALL CD_Get_Argument(index, value, arg_len, arg_stat)
    IF (arg_stat == -1) THEN
      WRITE (error_unit, '(A,I0,A,I0,A)') PROG//': argument ', index, ' is longer than ', LEN(value), &
        ' characters'
      CALL stop_run(1)
    ELSE IF (arg_stat /= 0) THEN
      WRITE (error_unit, '(A,I0)') PROG//': cannot read command-line argument ', index
      CALL stop_run(1)
    END IF
    IF (arg_len < 1) THEN
      WRITE (error_unit, '(A,I0,A)') PROG//': command-line argument ', index, ' is empty'
      CALL stop_run(1)
    END IF
  END SUBROUTINE get_argument

  SUBROUTINE prepare_paths(deck, root, deck_spelling, root_spelling)
    !! The spellings the Fortran runtime opens for the UTF-8 deck path and output root
    !! (CableDyn_PathIO). A name that cannot be opened exactly -- a reserved Windows
    !! device, characters outside the ANSI code page with no 8.3 short name, or a
    !! path beyond the length limit -- stops the run here with the reason, so a
    !! look-alike file is never read or written instead.
    CHARACTER(*), INTENT(IN) :: deck, root
    CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: deck_spelling, root_spelling
    CHARACTER(512) :: why
    INTEGER :: stat

    CALL CD_Native_Path(deck, deck_spelling, stat, why)
    IF (stat /= CD_PATH_OK) THEN
      WRITE (error_unit, '(A)') PROG//': cannot read deck "'//deck//'": '//TRIM(why)
      CALL stop_run(1)
    END IF
    CALL CD_Native_Path(root, root_spelling, stat, why, for_output=.TRUE., reserve=MAX_OUTPUT_SUFFIX)
    IF (stat /= CD_PATH_OK) THEN
      WRITE (error_unit, '(A)') PROG//': cannot write output files at "'//root//'": '//TRIM(why)
      CALL stop_run(1)
    END IF
  END SUBROUTINE prepare_paths

  SUBROUTINE check_output_root(deck, deck_spelling, root, root_spelling)
    !! Validate the output root before any solve: no output file may be the deck itself,
    !! and the output directory must accept a new file. deck/root are the UTF-8 names,
    !! the *_spelling arguments the names the runtime opens.
    CHARACTER(*), INTENT(IN) :: deck, deck_spelling, root, root_spelling
    CHARACTER(:), ALLOCATABLE :: probe
    INTEGER :: unit, ios
    LOGICAL :: exists

    CALL check_input_not_output(deck, root, root_spelling, 'deck', deck_spelling)
    ! A probe file beside the outputs proves the directory exists and is writable.
    probe = root_spelling//'.write_check.tmp'
    INQUIRE (FILE=probe, EXIST=exists)
    OPEN (NEWUNIT=unit, FILE=probe, STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE (error_unit, '(A)') PROG//': cannot write output files at "'//root// &
        '" (check that the directory exists and is writable)'
      CALL stop_run(1)
    END IF
    IF (exists) THEN
      CLOSE (unit)
    ELSE
      CLOSE (unit, STATUS='DELETE')
    END IF
  END SUBROUTINE check_output_root

  SUBROUTINE check_input_not_output(input, root, root_spelling, what, input_spelling)
    !! Stop before any solve when the input file (the deck or a file the deck names)
    !! is one the run would overwrite. The run writes <root>.out, <root>.static.out,
    !! <root>.elements.out, <root>.Line<N>.p.out / .t.out / .range.out, <root>.Rod<N>.p.out, plus
    !! the transient <root>.write_check.tmp probe and the <root>.cabledyn.lock lock.
    !! Each fixed name is compared with the input by file identity; the numbered
    !! per-line and per-rod names are recognized from the input's own name and then
    !! confirmed by identity, so every spelling of one file (relative, absolute,
    !! differently cased on Windows) is caught.
    CHARACTER(*), INTENT(IN) :: input, root, root_spelling, what
    CHARACTER(*), INTENT(IN), OPTIONAL :: input_spelling
    ! The fixed output suffixes, blank-separated.
    CHARACTER(*), PARAMETER :: FIXED = '.out .static.out .elements.out .write_check.tmp '//LOCK_SUFFIX//' '
    CHARACTER(:), ALLOCATABLE :: spelling, base_in, base_root, suffix
    CHARACTER(512) :: why
    INTEGER :: stat, first, last
    LOGICAL :: clash

    IF (PRESENT(input_spelling)) THEN
      spelling = input_spelling
    ELSE
      CALL CD_Native_Path(input, spelling, stat, why)
      ! A name the runtime cannot spell is not one of the outputs, which it can.
      IF (stat /= CD_PATH_OK) RETURN
    END IF
    clash = .FALSE.
    first = 1
    DO WHILE (first < LEN(FIXED))
      last = first + INDEX(FIXED(first:), ' ') - 2
      clash = same_file(spelling, root_spelling//FIXED(first:last))
      IF (clash) EXIT
      first = last + 2
    END DO
    IF (.NOT. clash) THEN
      base_in = file_name(input)
      base_root = file_name(root)
      IF (LEN(base_in) > LEN(base_root)) THEN
        IF (lower(base_in(1:LEN(base_root))) == lower(base_root)) THEN
          suffix = base_in(LEN(base_root) + 1:)
          IF (is_numbered_output_suffix(lower(suffix))) clash = same_file(spelling, root_spelling//suffix)
        END IF
      END IF
    END IF
    IF (clash) THEN
      IF (what == 'deck') THEN
        WRITE (error_unit, '(A)') PROG//': output root "'//root//'" would overwrite the input deck'
      ELSE
        WRITE (error_unit, '(A)') PROG//': output root "'//root//'" would overwrite the '//what//' "'// &
          input//'" that the deck reads'
      END IF
      CALL stop_run(1)
    END IF
  END SUBROUTINE check_input_not_output

  PURE LOGICAL FUNCTION is_numbered_output_suffix(suffix) RESULT(yes)
    !! True for ".line<digits>.p.out", ".line<digits>.t.out", ".line<digits>.range.out" and
    !! ".rod<digits>.p.out" (lower case), the per-line and per-rod output names.
    CHARACTER(*), INTENT(IN) :: suffix
    INTEGER :: first, n
    yes = .FALSE.
    n = LEN(suffix)
    IF (n > 5) THEN
      IF (suffix(1:5) == '.line') THEN
        first = 6
      ELSE IF (suffix(1:4) == '.rod') THEN
        first = 5
      ELSE
        RETURN
      END IF
    ELSE
      RETURN
    END IF
    IF (suffix(1:5) == '.line' .AND. n - first + 1 >= 11) THEN
      IF (suffix(n - 9:n) == '.range.out' .AND. VERIFY(suffix(first:n - 10), '0123456789') == 0) THEN
        yes = .TRUE.
        RETURN
      END IF
    END IF
    IF (n - first + 1 < 7) RETURN
    IF (VERIFY(suffix(first:n - 6), '0123456789') /= 0) RETURN
    IF (suffix(n - 5:n) == '.p.out') THEN
      yes = .TRUE.
    ELSE IF (suffix(n - 5:n) == '.t.out') THEN
      yes = suffix(1:5) == '.line'
    END IF
  END FUNCTION is_numbered_output_suffix

  PURE FUNCTION file_name(path) RESULT(name)
    !! The last component of a path (after the last '/' or '\').
    CHARACTER(*), INTENT(IN) :: path
    CHARACTER(:), ALLOCATABLE :: name
    INTEGER :: p
    p = SCAN(path, '/\', BACK=.TRUE.)
    name = path(p + 1:)
  END FUNCTION file_name

  PURE FUNCTION lower(s) RESULT(out)
    !! ASCII lower case.
    CHARACTER(*), INTENT(IN) :: s
    CHARACTER(LEN(s)) :: out
    INTEGER :: i, c
    out = s
    DO i = 1, LEN(s)
      c = IACHAR(s(i:i))
      IF (c >= IACHAR('A') .AND. c <= IACHAR('Z')) out(i:i) = ACHAR(c + 32)
    END DO
  END FUNCTION lower

  SUBROUTINE lock_output_root(root)
    !! Hold <root>.cabledyn.lock for the rest of the run, so two runs writing one output
    !! root cannot interleave their files. The operating system releases the lock when
    !! the process ends, even when it is killed, so a stale lock never blocks a run.
    CHARACTER(*), INTENT(IN) :: root
    SELECT CASE (CD_Lock_Path(root//LOCK_SUFFIX))
    CASE (CD_LOCK_ACQUIRED)
      RETURN
    CASE (CD_LOCK_HELD)
      WRITE (error_unit, '(A)') PROG//': another CableDyn run is writing output root "'//root// &
        '" (it holds '//root//LOCK_SUFFIX//'); wait for it to finish or choose another output root'
    CASE DEFAULT
      WRITE (error_unit, '(A)') PROG//': cannot create the output lock file "'//root//LOCK_SUFFIX//'"'
    END SELECT
    CALL stop_run(1)
  END SUBROUTINE lock_output_root

  LOGICAL FUNCTION same_file(a, b) RESULT(same)
    !! True when both names reach one existing file. With a connected to a unit, an
    !! INQUIRE by the name b reports OPENED when b resolves to that same file, however
    !! it is spelled (relative, absolute, or differently cased on Windows); normalized
    !! spelling is compared as well.
    CHARACTER(*), INTENT(IN) :: a, b
    LOGICAL :: exists_a, exists_b, opened
    INTEGER :: unit, ios
    same = .FALSE.
    INQUIRE (FILE=a, EXIST=exists_a)
    INQUIRE (FILE=b, EXIST=exists_b)
    IF (.NOT. (exists_a .AND. exists_b)) RETURN
    same = normalized_path(a) == normalized_path(b)
    IF (same) RETURN
    OPEN (NEWUNIT=unit, FILE=a, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    opened = .FALSE.
    INQUIRE (FILE=b, OPENED=opened)
    CLOSE (unit)
    same = opened
  END FUNCTION same_file

  FUNCTION normalized_path(path) RESULT(norm)
    CHARACTER(*), INTENT(IN) :: path
    CHARACTER(LEN(path)) :: norm
    INTEGER :: i, c
    norm = ADJUSTL(path)
    DO i = 1, LEN_TRIM(norm)
      c = IACHAR(norm(i:i))
      IF (norm(i:i) == '\') THEN
        norm(i:i) = '/'
      ELSE IF (c >= IACHAR('A') .AND. c <= IACHAR('Z')) THEN
        norm(i:i) = ACHAR(c + 32)
      END IF
    END DO
    ! drop leading "./" segments so "deck.dat" and "./deck.dat" compare equal
    DO WHILE (LEN_TRIM(norm) > 2)
      IF (norm(1:2) /= './') EXIT
      norm = norm(3:)
    END DO
  END FUNCTION normalized_path

  SUBROUTINE run_mixed_held_deck(deck, root, dt, has_dt, run_time, has_run_time, has_motion, line_outputs, &
                                 run_converged, run_stat, run_msg)
    !! Standalone owner for a mixed {EI=0 + finite-EI} deck. It deliberately reuses
    !! the production OpenFAST aggregate and holds the initialized coupled endpoints.
    !! Prescribed reference-body motion and deck-owned ambient fields require separate
    !! transformations/updaters and therefore remain explicit fail-closed boundaries.
    CHARACTER(*), INTENT(IN) :: deck, root
    REAL(wp), INTENT(IN) :: dt, run_time
    LOGICAL, INTENT(IN) :: has_dt, has_run_time, has_motion, line_outputs
    LOGICAL, INTENT(OUT) :: run_converged
    INTEGER, INTENT(OUT) :: run_stat
    CHARACTER(*), INTENT(OUT) :: run_msg

    TYPE(CD_AGG_ModuleType), ALLOCATABLE :: aggregate   ! about 37 KB: on the heap, not the stack
    REAL(wp), ALLOCATABLE :: position(:, :), velocity(:, :), acceleration(:, :), load(:, :)
    REAL(wp) :: time, first_miss_time, last_miss_time, step_dt
    INTEGER :: es, n_moving, nstep, step, n_iter, miss_count, consecutive_miss, max_consecutive_miss
    INTEGER :: unit, ios
    LOGICAL :: step_converged, stalled, file_open
    CHARACTER(512) :: em
    CHARACTER(4096) :: static_path, dynamic_path
    TYPE(StandaloneProgress) :: progress

    run_converged = .FALSE.
    run_stat = CD_DECKDRV_BADINPUT
    run_msg = ''
    file_open = .FALSE.
    unit = 0
    IF (has_motion) THEN
      run_msg = 'CableDyn_driver: mixed-deck motionFile is not supported; platform-reference motion '// &
                'must be rigidly transformed to every coupled hang-off before marching'
      RETURN
    END IF
    IF (line_outputs) THEN
      run_msg = 'CableDyn_driver: mixed-deck LINE Outputs p/t are not supported; '// &
                'use main OUTPUTS channels or remove the per-line flags'
      RETURN
    END IF
    ! The mixed route writes no modal files: fail closed rather than ignore the request.
    IF (deck_requests_modes(deck)) THEN
      run_msg = 'CableDyn_driver: OPTION nModes is not supported for mixed EI=0/finite-EI decks; '// &
                'remove it or set it to 0'
      RETURN
    END IF
    IF (has_run_time .AND. run_time < CD_ZERO) THEN
      run_msg = 'CableDyn_driver: mixed-deck TMax must be non-negative'
      RETURN
    END IF
    nstep = 0
    step_dt = MIXED_STATIC_INIT_DT
    IF (has_dt .AND. dt > CD_ZERO) step_dt = dt
    IF (has_run_time .AND. run_time > CD_ZERO) THEN
      IF (.NOT. has_dt .OR. .NOT. (dt > CD_ZERO)) THEN
        run_msg = 'CableDyn_driver: a positive mixed-deck TMax requires a finite positive dtM'
        RETURN
      END IF
      IF (.NOT. (run_time/step_dt <= REAL(HUGE(1) - 1, wp))) THEN
        run_msg = 'CableDyn_driver: mixed-deck TMax/dtM needs more time steps than the driver can count; '// &
                  'increase dtM or reduce TMax'
        RETURN
      END IF
      nstep = NINT(run_time/step_dt)
    END IF

    mixed_run: BLOCK
      ALLOCATE (aggregate, STAT=ios)
      IF (ios /= 0) THEN
        run_stat = CD_DECKDRV_SOLVEFAIL
        run_msg = 'CableDyn_driver: cannot allocate the mixed aggregate'
        EXIT mixed_run
      END IF
      ! The aggregate is caller-driven and cannot refresh standalone wave/current
      ! OPTIONS. Reject such decks by name instead of silently holding zero fluid.
      CALL CD_AGG_Init_From_Deck(aggregate, deck, step_dt, es, em, forbid_deck_ambient=.TRUE., &
                                 run_tmax=MAX(CD_ZERO, run_time), range_files=.TRUE.)
      IF (es /= CD_AGG_OK) THEN
        CALL set_aggregate_error(es, em, run_stat, run_msg)
        EXIT mixed_run
      END IF
      CALL print_mixed_initialization(aggregate, run_stat, run_msg)
      IF (run_stat /= CD_DECKDRV_OK) EXIT mixed_run

      n_moving = CD_AGG_NMovingPoints(aggregate, es, em)
      IF (es /= CD_AGG_OK) THEN
        CALL set_aggregate_error(es, em, run_stat, run_msg)
        EXIT mixed_run
      END IF
      ALLOCATE (position(3, n_moving), velocity(3, n_moving), acceleration(3, n_moving), &
                load(3, n_moving), STAT=ios)
      IF (ios /= 0) THEN
        run_stat = CD_DECKDRV_SOLVEFAIL
        run_msg = 'CableDyn_driver: cannot allocate the mixed aggregate boundary workspace'
        EXIT mixed_run
      END IF
      CALL CD_AGG_GetMovingPointMesh(aggregate, position, velocity, acceleration, load, es, em)
      IF (es /= CD_AGG_OK) THEN
        CALL set_aggregate_error(es, em, run_stat, run_msg)
        EXIT mixed_run
      END IF
      ! A standalone deck with no motionFile means held coupling boundaries. Initial
      ! solver accelerations are not prescribed boundary accelerations and must not be
      ! fed back into the first step.
      velocity = CD_ZERO
      acceleration = CD_ZERO

      static_path = TRIM(root)//'.static.out'
      CALL CD_AGG_WriteStaticProfile(aggregate, TRIM(static_path), es, em)
      IF (es /= CD_AGG_OK) THEN
        CALL set_aggregate_error(es, em, run_stat, run_msg)
        EXIT mixed_run
      END IF
      dynamic_path = TRIM(root)//'.out'
      OPEN (NEWUNIT=unit, FILE=TRIM(dynamic_path), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) THEN
        run_stat = CD_DECKDRV_BADINPUT
        run_msg = 'CableDyn_driver: cannot open mixed dynamic output: '//TRIM(dynamic_path)
        EXIT mixed_run
      END IF
      file_open = .TRUE.
      WRITE (unit, '(A)', IOSTAT=ios) &
        '# CableDyn driver output (mixed aggregate dynamic; convergence quality is reported by status message)'
      IF (ios /= 0) THEN
        run_stat = CD_DECKDRV_BADINPUT
        run_msg = 'CableDyn_driver: cannot write mixed dynamic output title'
        EXIT mixed_run
      END IF
      CALL write_mixed_header(unit, aggregate, es, em)
      IF (es /= CD_AGG_OK) THEN
        CALL set_aggregate_error(es, em, run_stat, run_msg)
        EXIT mixed_run
      END IF
      CALL CD_AGG_CalcOutput(aggregate, es, em)
      IF (es /= CD_AGG_OK) THEN
        CALL set_aggregate_error(es, em, run_stat, run_msg)
        EXIT mixed_run
      END IF
      CALL write_mixed_row(unit, CD_ZERO, aggregate, es, em)
      IF (es == CD_AGG_OK) CALL CD_AGG_Range_Sample(aggregate, CD_ZERO, es, em)
      IF (es /= CD_AGG_OK) THEN
        CALL set_aggregate_error(es, em, run_stat, run_msg)
        EXIT mixed_run
      END IF

      miss_count = 0
      consecutive_miss = 0
      max_consecutive_miss = 0
      first_miss_time = CD_ZERO
      last_miss_time = CD_ZERO
      IF (nstep > 0) CALL start_standalone_progress(progress, nstep, step_dt)
      DO step = 1, nstep
        time = REAL(step, wp)*step_dt
        CALL CD_AGG_Step_Moving(aggregate, step_dt, position, velocity, acceleration, step_converged, stalled, &
                                n_iter, es, em, t_committed=time)
        IF (es /= CD_AGG_OK) THEN
          CALL set_aggregate_error(es, em, run_stat, run_msg)
          EXIT mixed_run
        END IF
        IF (.NOT. step_converged .OR. stalled) THEN
          miss_count = miss_count + 1
          consecutive_miss = consecutive_miss + 1
          max_consecutive_miss = MAX(max_consecutive_miss, consecutive_miss)
          IF (miss_count == 1) first_miss_time = time
          last_miss_time = time
          CALL CD_Classify_Dynamic_Completion('mixed aggregate dynamics', nstep, miss_count, &
                                              max_consecutive_miss, first_miss_time, last_miss_time, &
                                              run_converged, run_stat, run_msg)
          EXIT mixed_run
        ELSE
          consecutive_miss = 0
        END IF
        CALL CD_AGG_CalcOutput(aggregate, es, em)
        IF (es /= CD_AGG_OK) THEN
          CALL set_aggregate_error(es, em, run_stat, run_msg)
          EXIT mixed_run
        END IF
        CALL write_mixed_row(unit, time, aggregate, es, em)
        IF (es == CD_AGG_OK) CALL CD_AGG_Range_Sample(aggregate, time, es, em)
        IF (es /= CD_AGG_OK) THEN
          CALL set_aggregate_error(es, em, run_stat, run_msg)
          EXIT mixed_run
        END IF
        CALL update_standalone_progress(progress, step, time)
      END DO

      CALL CD_Classify_Dynamic_Completion('mixed aggregate dynamics', nstep, miss_count, &
                                          max_consecutive_miss, first_miss_time, last_miss_time, &
                                          run_converged, run_stat, run_msg)
      IF (nstep == 0) run_converged = .TRUE.
      IF (run_converged .AND. run_stat == CD_DECKDRV_OK) THEN
        CALL CD_AGG_Range_Write(aggregate, root, es, em)
        IF (es /= CD_AGG_OK) THEN
          run_converged = .FALSE.
          CALL set_aggregate_error(es, em, run_stat, run_msg)
        END IF
      END IF
    END BLOCK mixed_run

    IF (file_open) CLOSE (unit)
    IF (ALLOCATED(aggregate)) CALL CD_AGG_End(aggregate, es, em)
  END SUBROUTINE run_mixed_held_deck

  LOGICAL FUNCTION deck_requests_modes(deck) RESULT(requested)
    !! TRUE when the deck's OPTIONS section carries a nonzero nModes row (MoorDyn-style
    !! "<value> nModes"; the keyword spellings the deck parser accepts). The deck parser
    !! validates the value; this scan only detects the request.
    CHARACTER(*), INTENT(IN) :: deck
    CHARACTER(4096) :: line
    CHARACTER(64) :: t1, t2
    INTEGER :: u, ios
    LOGICAL :: in_options
    REAL(wp) :: x
    requested = .FALSE.
    in_options = .FALSE.
    OPEN (NEWUNIT=u, FILE=deck, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    DO
      READ (u, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      line = ADJUSTL(line)
      IF (line(1:3) == '---') THEN
        in_options = INDEX(ascii_lower(line), 'option') > 0
        CYCLE
      END IF
      IF (.NOT. in_options) CYCLE
      t1 = ''
      t2 = ''
      READ (line, *, IOSTAT=ios) t1, t2
      IF (ios /= 0) CYCLE
      SELECT CASE (ascii_lower(t2))
      CASE ('nmodes', 'n_modes', 'modes')
        READ (t1, *, IOSTAT=ios) x
        IF (ios == 0) THEN
          IF (ABS(x) > CD_ZERO) requested = .TRUE.
        END IF
      END SELECT
      IF (requested) EXIT
    END DO
    CLOSE (u)
  END FUNCTION deck_requests_modes

  FUNCTION ascii_lower(s) RESULT(r)
    !! ASCII lower case of s.
    CHARACTER(*), INTENT(IN) :: s
    CHARACTER(LEN(s)) :: r
    INTEGER :: k
    r = s
    DO k = 1, LEN(s)
      IF (s(k:k) >= 'A' .AND. s(k:k) <= 'Z') r(k:k) = ACHAR(IACHAR(s(k:k)) + 32)
    END DO
  END FUNCTION ascii_lower

  SUBROUTINE print_mixed_initialization(aggregate, stat, msg)
    TYPE(CD_AGG_ModuleType), INTENT(IN) :: aggregate
    INTEGER, INTENT(OUT) :: stat
    CHARACTER(*), INTENT(OUT) :: msg
    INTEGER :: n_lines, n_points, n_sections, n_ei0, n_finite, i, line_id, es
    REAL(wp) :: tension, force(3), inclination, declination, azimuth
    CHARACTER(512) :: em
    CHARACTER(1024) :: note

    CALL CD_AGG_GetInitMetadata(aggregate, n_lines, n_points, n_sections, n_ei0, n_finite, es, em)
    IF (es /= CD_AGG_OK) THEN
      CALL set_aggregate_error(es, em, stat, msg)
      RETURN
    END IF
    WRITE (output_unit, '(A,I0,A,I0,A,I0,A)') '  CableDyn mixed standalone aggregate: ', n_lines, &
      ' line(s) [', n_ei0, ' EI=0, ', n_finite, ' finite-EI].'
    WRITE (output_unit, '(A,I0,A,I0,A)') '  Parsed ', n_points, ' point(s) and ', n_sections, ' section row(s).'
    WRITE (output_unit, '(A)') '  Static equilibrium fairlead results:'
    DO i = 1, n_lines
      CALL CD_AGG_GetInitLine(aggregate, i, line_id, tension, force, inclination, declination, azimuth, es, em, &
                              note=note)
      IF (es /= CD_AGG_OK) THEN
        CALL set_aggregate_error(es, em, stat, msg)
        RETURN
      END IF
      WRITE (output_unit, '(A,I0,A,ES12.5,A,F8.3,A,3(ES12.5,1X),A)') '    Line ', line_id, &
        ': FairTen=', tension, ' N, tangent inclination=', inclination, ' deg, force=(', force, ') N'
      IF (LEN_TRIM(note) > 0) CALL print_notes(note)
    END DO
    FLUSH (output_unit)
    stat = CD_DECKDRV_OK
    msg = ''
  END SUBROUTINE print_mixed_initialization

  SUBROUTINE print_notes(text)
    !! Print the initialisation notes (one per text line) as a standalone run does.
    CHARACTER(*), INTENT(IN) :: text
    INTEGER :: a, b
    a = 1
    DO
      b = INDEX(text(a:), NEW_LINE('a'))
      IF (b == 0) THEN
        WRITE (output_unit, '(A)') '      Note: '//TRIM(text(a:))//'.'
        EXIT
      END IF
      WRITE (output_unit, '(A)') '      Note: '//text(a:a + b - 2)//'.'
      a = a + b
    END DO
  END SUBROUTINE print_notes

  SUBROUTINE write_mixed_header(unit, aggregate, stat, msg)
    INTEGER, INTENT(IN) :: unit
    TYPE(CD_AGG_ModuleType), INTENT(IN) :: aggregate
    INTEGER, INTENT(OUT) :: stat
    CHARACTER(*), INTENT(OUT) :: msg
    INTEGER :: i, ios
    CHARACTER(64) :: header, unit_label

    WRITE (unit, '(A)', ADVANCE='NO', IOSTAT=ios) 'Time(s)'
    IF (ios /= 0) THEN
      stat = CD_AGG_BADINPUT; msg = 'cannot write mixed output header'; RETURN
    END IF
    DO i = 1, CD_AGG_NumChannels(aggregate)
      CALL CD_AGG_ChannelHeader(aggregate, i, header, unit_label, stat, msg)
      IF (stat /= CD_AGG_OK) RETURN
      WRITE (unit, '(A,A)', ADVANCE='NO', IOSTAT=ios) CHAR(9), TRIM(header)
      IF (ios /= 0) THEN
        stat = CD_AGG_BADINPUT; msg = 'cannot write mixed output header'; RETURN
      END IF
    END DO
    WRITE (unit, '(A)', IOSTAT=ios) ''
    IF (ios /= 0) THEN
      stat = CD_AGG_BADINPUT; msg = 'cannot write mixed output header'; RETURN
    END IF
    stat = CD_AGG_OK
    msg = ''
  END SUBROUTINE write_mixed_header

  SUBROUTINE write_mixed_row(unit, time, aggregate, stat, msg)
    INTEGER, INTENT(IN) :: unit
    REAL(wp), INTENT(IN) :: time
    TYPE(CD_AGG_ModuleType), INTENT(IN) :: aggregate
    INTEGER, INTENT(OUT) :: stat
    CHARACTER(*), INTENT(OUT) :: msg
    INTEGER :: i, ios
    REAL(wp) :: value

    ! full-precision time column (the deck-driver writers use the same format)
    WRITE (unit, '(ES25.16E3)', ADVANCE='NO', IOSTAT=ios) time
    IF (ios /= 0) THEN
      stat = CD_AGG_BADINPUT; msg = 'cannot write mixed output row'; RETURN
    END IF
    DO i = 1, CD_AGG_NumChannels(aggregate)
      CALL CD_AGG_EvalChannel(aggregate, i, value, stat, msg)
      IF (stat /= CD_AGG_OK) RETURN
      WRITE (unit, '(A,ES15.7)', ADVANCE='NO', IOSTAT=ios) CHAR(9), value
      IF (ios /= 0) THEN
        stat = CD_AGG_BADINPUT; msg = 'cannot write mixed output row'; RETURN
      END IF
    END DO
    WRITE (unit, '(A)', IOSTAT=ios) ''
    IF (ios /= 0) THEN
      stat = CD_AGG_BADINPUT; msg = 'cannot write mixed output row'; RETURN
    END IF
    stat = CD_AGG_OK
    msg = ''
  END SUBROUTINE write_mixed_row

  SUBROUTINE set_aggregate_error(aggregate_stat, aggregate_msg, stat, msg)
    INTEGER, INTENT(IN) :: aggregate_stat
    CHARACTER(*), INTENT(IN) :: aggregate_msg
    INTEGER, INTENT(OUT) :: stat
    CHARACTER(*), INTENT(OUT) :: msg
    IF (aggregate_stat == CD_AGG_BADINPUT) THEN
      stat = CD_DECKDRV_BADINPUT
    ELSE
      stat = CD_DECKDRV_SOLVEFAIL
    END IF
    msg = aggregate_msg
  END SUBROUTINE set_aggregate_error

  SUBROUTINE start_standalone_progress(progress, nstep, dt)
    TYPE(StandaloneProgress), INTENT(OUT) :: progress
    INTEGER, INTENT(IN) :: nstep
    REAL(wp), INTENT(IN) :: dt
    progress%nstep = nstep
    progress%stride = MAX(1, (nstep + 19)/20)
    CALL SYSTEM_CLOCK(progress%start_count, progress%count_rate)
    WRITE (output_unit, '(A,I0,A,F12.3,A,ES12.5,A)') '  Dynamic simulation: ', nstep, &
      ' step(s), simulated duration ', REAL(nstep, wp)*dt, ' s, dtM = ', dt, ' s.'
    FLUSH (output_unit)
  END SUBROUTINE start_standalone_progress

  SUBROUTINE update_standalone_progress(progress, step, simulated_time)
    TYPE(StandaloneProgress), INTENT(IN) :: progress
    INTEGER, INTENT(IN) :: step
    REAL(wp), INTENT(IN) :: simulated_time
    INTEGER(int64) :: now
    REAL(wp) :: elapsed, remaining, percent
    CHARACTER(16) :: elapsed_text, remaining_text
    ! Every committed step: the time an abnormal end of the process reports.
    CALL CD_Fatal_Report_Time(simulated_time)
    IF (step < progress%nstep .AND. MOD(step, progress%stride) /= 0) RETURN
    CALL SYSTEM_CLOCK(now)
    elapsed = CD_ZERO
    IF (progress%count_rate > 0_int64 .AND. now >= progress%start_count) &
      elapsed = REAL(now - progress%start_count, wp)/REAL(progress%count_rate, wp)
    remaining = CD_ZERO
    IF (step < progress%nstep .AND. elapsed > CD_ZERO) &
      remaining = elapsed*REAL(progress%nstep - step, wp)/REAL(step, wp)
    percent = 100.0_wp*REAL(step, wp)/REAL(progress%nstep, wp)
    elapsed_text = format_wall_time(elapsed)
    remaining_text = format_wall_time(remaining)
    WRITE (output_unit, '(A,F6.1,A,F12.3,A,A,A,A)') '  Progress: ', percent, '% | t = ', &
      simulated_time, ' s | elapsed ', TRIM(elapsed_text), ' | ETA ', TRIM(remaining_text)
    FLUSH (output_unit)
  END SUBROUTINE update_standalone_progress

  PURE FUNCTION format_wall_time(seconds) RESULT(text)
    REAL(wp), INTENT(IN) :: seconds
    CHARACTER(16) :: text
    INTEGER :: total, hours, minutes, secs
    total = MAX(0, NINT(seconds))
    hours = total/3600
    minutes = MOD(total, 3600)/60
    secs = MOD(total, 60)
    WRITE (text, '(I3.3,A,I2.2,A,I2.2)') hours, ':', minutes, ':', secs
  END FUNCTION format_wall_time

  SUBROUTINE print_usage(unit)
    !! Command-line usage summary for the standalone CableDyn driver.
    INTEGER, INTENT(IN) :: unit
    WRITE (unit, '(A)') ''
    WRITE (unit, '(A)') 'Usage:'
    WRITE (unit, '(A)') '  CableDyn_driver <deck.dat> <out_root>   solve a deck; write <out_root>.out'
    WRITE (unit, '(A)') '  CableDyn_driver -v | --version          print version and exit'
    WRITE (unit, '(A)') '  CableDyn_driver -h | --help             print this help and exit'
    WRITE (unit, '(A)') ''
    WRITE (unit, '(A)') 'The deck follows MoorDyn v2 vocabulary with sectioned line objects'
    WRITE (unit, '(A)') '(a line runs End A -> End B through ordered SECTIONS). See examples/.'
  END SUBROUTINE print_usage

END PROGRAM cabledyn
