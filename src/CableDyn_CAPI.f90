! File: src/CableDyn_CAPI.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_CAPI
  !! ISO C binding for CableDyn lifecycle ownership and coupled exchange.
  !!
  !! Its reference is doc/capi.rst (the C API reference): C callers receive
  !! an opaque handle that owns the OpenFAST-facing CableDyn module object. The
  !! ABI exposes deck and raw-array line initialization plus update/step/output calls.
  !! Raw-line initialization remains explicit, while deck initialization preserves
  !! the point-system object graph used by OpenFAST-facing integrations.
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_BOOL, C_CHAR, C_DOUBLE, C_INT, C_NULL_CHAR, C_NULL_PTR, C_PTR, &
                                                                            C_ASSOCIATED, C_LOC, C_F_POINTER
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE CableDyn_DeckDriver, ONLY: CD_Init_Deck_System, CD_DECKDRV_OK, CD_DECKDRV_BADINPUT, CD_DECKDRV_SOLVEFAIL, &
                                 CD_Deck_Query_dtM, CD_DECK_NAMELEN, CD_Eval_Aggregate_Channel, &
                                 CD_Eval_Aggregate_Segment_Tensions, CD_Rigid6_Inventory, CD_Rod_Inventory, &
                                 CD_Channel_Token_Parses
  USE CableDyn_OpenFAST_Aggregate, ONLY: CD_AGG_ModuleType, CD_AGG_Init_From_Deck, CD_AGG_NMovingPoints, &
                                         CD_AGG_UpdateStates_Moving, CD_AGG_Step_Moving, CD_AGG_CalcOutput, &
                                         CD_AGG_GetMovingPointMesh, CD_AGG_End, CD_AGG_IsInitialized, CD_AGG_OK, &
                                         CD_AGG_BADINPUT, CD_AGG_NOT_INITIALIZED, CD_AGG_ALLOCFAIL
  USE CableDyn_Dynamic, ONLY: GenAlphaConfig
  USE CableDyn_Linalg, ONLY: CD_Blas_Runtime_Check, CD_LINALG_OK
  USE CableDyn_Model, ONLY: CD_ModelType, CD_Init_Model, CD_End_Model, CD_MODEL_OK, CD_MODEL_BADINPUT, &
                            CD_MODEL_SOLVEFAIL, CD_MODEL_ALLOCFAIL
  USE CableDyn_OpenFAST, ONLY: CD_FAST_ModuleType, CD_FAST_End, CD_FAST_Init_From_Models, &
                               CD_FAST_UpdateStates, CD_FAST_Step, CD_FAST_CalcOutput, &
                               CD_FAST_CalcOutputDerivatives, CD_FAST_GetCoupledMotion, CD_FAST_IsInitialized, &
                               CD_FAST_NCoupledDOF, CD_FAST_NPoints, CD_FAST_NLines, CD_FAST_OK, CD_FAST_BADINPUT, &
                               CD_FAST_NOT_INITIALIZED, CD_FAST_ALLOCFAIL, CD_FAST_UpdatePointFluidFields
  USE CableDyn_OpenFAST_HermiteFMF, ONLY: CD_HFMF_ModuleType
  USE CableDyn_Precision, ONLY: wp, CD_ZERO
  USE CableDyn_System, ONLY: CD_SystemType, CD_System_Line_NDOF
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CableDyn_Create
  PUBLIC :: CableDyn_Close
  PUBLIC :: CableDyn_InitDeck
  PUBLIC :: CableDyn_InitLine
  PUBLIC :: CableDyn_InitLines
  PUBLIC :: CableDyn_UpdateStates
  PUBLIC :: CableDyn_UpdatePointFluidFields
  PUBLIC :: CableDyn_Step
  PUBLIC :: CableDyn_CalcOutput
  PUBLIC :: CableDyn_CalcOutputDerivatives
  PUBLIC :: CableDyn_GetCoupledMotion
  PUBLIC :: CableDyn_GetLastError
  PUBLIC :: CableDyn_GetVersion
  PUBLIC :: CableDyn_GetVersionString
  PUBLIC :: CableDyn_IsInitialized
  PUBLIC :: CableDyn_NCoupledDOF
  PUBLIC :: CableDyn_NPoints
  PUBLIC :: CableDyn_NLines
  PUBLIC :: CableDyn_GetAbiMinor
  PUBLIC :: CableDyn_NObjects
  PUBLIC :: CableDyn_GetObjectInfo
  PUBLIC :: CableDyn_GetLineValues
  PUBLIC :: CableDyn_GetPointState
  PUBLIC :: CableDyn_GetBodyState
  PUBLIC :: CableDyn_GetRodState
  PUBLIC :: CableDyn_EvalChannel

  INTEGER(C_INT), PARAMETER :: CD_C_OK = 0_C_INT
  INTEGER(C_INT), PARAMETER :: CD_C_BAD_HANDLE = 1_C_INT
  INTEGER(C_INT), PARAMETER :: CD_C_ALLOC_FAIL = 2_C_INT
  INTEGER(C_INT), PARAMETER :: CD_C_BAD_INPUT = 3_C_INT
  INTEGER(C_INT), PARAMETER :: CD_C_SOLVE_FAIL = 4_C_INT
  INTEGER(C_INT), PARAMETER :: CD_C_NOT_INITIALIZED = 5_C_INT
  INTEGER(C_INT), PARAMETER :: CD_C_VERSION_MAJOR = 0_C_INT
  INTEGER(C_INT), PARAMETER :: CD_C_VERSION_MINOR = 1_C_INT
  INTEGER(C_INT), PARAMETER :: CD_C_VERSION_PATCH = 1_C_INT
  INTEGER(C_INT), PARAMETER :: CD_C_ABI_VERSION = 1_C_INT
  CHARACTER(*), PARAMETER :: CD_C_VERSION_STRING = 'CableDyn 0.1.1 C-ABI 1'
  !> Minor extension level of ABI 1: 1 adds the in-process object queries.
  INTEGER(C_INT), PARAMETER :: CD_C_ABI_MINOR = 1_C_INT
  !> Object kinds of CableDyn_NObjects / CableDyn_GetObjectInfo.
  INTEGER(C_INT), PARAMETER :: CD_C_OBJ_LINE = 1_C_INT, CD_C_OBJ_POINT = 2_C_INT, CD_C_OBJ_BODY = 3_C_INT, &
                               CD_C_OBJ_ROD = 4_C_INT
  !> Line quantities of CableDyn_GetLineValues (the codes of the line-node output channels).
  INTEGER(C_INT), PARAMETER :: CD_C_LINE_POSITION = 1_C_INT, CD_C_LINE_VELOCITY = 2_C_INT, &
                               CD_C_LINE_TENSION = 3_C_INT, CD_C_LINE_ACCELERATION = 4_C_INT, &
                               CD_C_LINE_CURVATURE = 5_C_INT, CD_C_LINE_BEND_MOMENT = 6_C_INT, &
                               CD_C_LINE_DECLINATION = 7_C_INT, CD_C_LINE_AZIMUTH = 8_C_INT, &
                               CD_C_LINE_SEGMENT_TENSION = 9_C_INT
  !> Capacity of the per-handle diagnostic retained for CableDyn_GetLastError.
  INTEGER, PARAMETER :: CD_C_MSG_LEN = 1024

  ! Live-handle integrity. A handle is valid only if its C_PTR is in the registry of
  ! handles this module created and has not closed (so a foreign, stale, or
  ! double-closed pointer is rejected WITHOUT being dereferenced) AND the target's
  ! magic field is intact (defence in depth once the registry says the deref is safe).
  ! Concurrency: registry operations run under a process-local C11 atomic mutex, and Close uses a
  ! single atomic check-and-remove (registry_claim) so a concurrent Close of the same
  ! or a copied handle has exactly one winner that frees -- no double-free. This
  ! synchronization is independent of the optional OpenMP performance runtime, so handle
  ! lifecycle stays thread-safe in serial/non-OpenMP builds too. (Calling
  ! a query on a handle while another thread is Closing it is a use-after-free the caller
  ! must avoid, as with any handle API; the registry defends against stale/foreign/
  ! double-closed handles, not against racing a live handle against its own teardown.)
  ! Residual limitation: because the ABI hands back a raw C_PTR, a stale copy whose
  ! freed address is later reused by a new Create would alias the live handle and pass
  ! -- inherent to a pointer ABI; an opaque integer-id ABI would remove it.
  INTEGER(INT64), PARAMETER :: CD_C_HANDLE_MAGIC = 3400658189_INT64    ! 0xCAB1ED0D, a distinctive sentinel
  TYPE(C_PTR), ALLOCATABLE, SAVE :: g_handle_registry(:)
  INTEGER, SAVE :: g_handle_count = 0

  INTERFACE
    SUBROUTINE cabledyn_registry_lock() BIND(C, NAME='cabledyn_registry_lock')
    END SUBROUTINE cabledyn_registry_lock
    SUBROUTINE cabledyn_registry_unlock() BIND(C, NAME='cabledyn_registry_unlock')
    END SUBROUTINE cabledyn_registry_unlock
    ! Process-wide blocking lock around deck parsing and handle initialisation
    ! (src/cabledyn_mutex.c). The Fortran runtime's file and formatted I/O is not
    ! safe to run from several threads at once, so InitDeck/InitLine/InitLines hold
    ! it; Step and the other per-handle calls do no I/O on their success paths and
    ! do not take it.
    SUBROUTINE cabledyn_input_lock() BIND(C, NAME='cabledyn_input_lock')
    END SUBROUTINE cabledyn_input_lock
    SUBROUTINE cabledyn_input_unlock() BIND(C, NAME='cabledyn_input_unlock')
    END SUBROUTINE cabledyn_input_unlock
    ! Limits OpenBLAS to one thread once per process and, where CableDyn loads
    ! OpenBLAS at run time, loads it with that pool size (src/cabledyn_blas.c).
    ! Returns 0 on success, 1 when that OpenBLAS runtime could not be loaded.
    FUNCTION cabledyn_blas_threads_once() RESULT(failed) BIND(C, NAME='cabledyn_blas_threads_once')
      IMPORT :: C_INT
      INTEGER(C_INT) :: failed
    END FUNCTION cabledyn_blas_threads_once
  END INTERFACE

  TYPE :: CD_C_Handle
    INTEGER(INT64) :: magic = 0_INT64
    TYPE(CD_FAST_ModuleType) :: fast
    REAL(wp), ALLOCATABLE :: q_work(:)
    REAL(wp), ALLOCATABLE :: v_work(:)
    REAL(wp), ALLOCATABLE :: a_work(:)
    REAL(wp), ALLOCATABLE :: loads_work(:)
    REAL(wp), ALLOCATABLE :: jq_work(:, :)
    REAL(wp), ALLOCATABLE :: jv_work(:, :)
    REAL(wp), ALLOCATABLE :: ja_work(:, :)
    REAL(wp), ALLOCATABLE :: madd_work(:, :)
    REAL(wp), ALLOCATABLE :: fluid_velocity_work(:, :)
    REAL(wp), ALLOCATABLE :: fluid_acceleration_work(:, :)
    REAL(wp), ALLOCATABLE :: waterline_work(:)
    ! Decks with BODIES or RODS run on the coupled aggregate (use_agg): its moving points are the
    ! coupled DOFs, its step is the deck dtM (agg_dt), and agg_t is its committed time.
    LOGICAL :: use_agg = .FALSE.
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp) :: agg_dt = CD_ZERO, agg_t = CD_ZERO
    REAL(wp), ALLOCATABLE :: agg_pos(:, :), agg_vel(:, :), agg_acc(:, :), agg_load(:, :)
    ! Point-route line inventory for the object queries: the deck-line-id -> system-line map
    ! (pt_obj, 0 = no such line; pt_is_cable all .FALSE.) the shared channel evaluator resolves
    ! through, and an empty cable store to pass it. Raw-line handles number their lines 1..n.
    LOGICAL, ALLOCATABLE :: pt_is_cable(:)
    INTEGER, ALLOCATABLE :: pt_obj(:)
    TYPE(CD_HFMF_ModuleType), ALLOCATABLE :: no_cables(:)
    CHARACTER(CD_C_MSG_LEN) :: last_error = ''
  END TYPE CD_C_Handle

CONTAINS

  SUBROUTINE CableDyn_GetVersion(major, minor, patch, abi_version) BIND(C, NAME='CableDyn_GetVersion')
    !! Return the CableDyn C ABI version tuple.
    INTEGER(C_INT), INTENT(OUT) :: major, minor, patch, abi_version

    major = CD_C_VERSION_MAJOR
    minor = CD_C_VERSION_MINOR
    patch = CD_C_VERSION_PATCH
    abi_version = CD_C_ABI_VERSION
  END SUBROUTINE CableDyn_GetVersion

  SUBROUTINE CableDyn_GetVersionString(version_ptr, version_len) BIND(C, NAME='CableDyn_GetVersionString')
    !! Copy a null-terminated human-readable CableDyn C ABI version string.
    TYPE(C_PTR), VALUE :: version_ptr
    INTEGER(C_INT), VALUE :: version_len

    CHARACTER(KIND=C_CHAR), POINTER :: version_c(:)

    IF (version_len < 1_C_INT .OR. .NOT. C_ASSOCIATED(version_ptr)) RETURN
    CALL C_F_POINTER(version_ptr, version_c, [INT(version_len)])
    CALL copy_message_to_c(CD_C_VERSION_STRING, version_c)
  END SUBROUTINE CableDyn_GetVersionString

  SUBROUTINE CableDyn_Create(handle, err_stat) BIND(C, NAME='CableDyn_Create')
    !! Allocate a CableDyn C handle. The returned handle is uninitialized until
    !! CableDyn_InitDeck, CableDyn_InitLine, or CableDyn_InitLines succeeds.
    TYPE(C_PTR), INTENT(OUT) :: handle
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    INTEGER :: stat

    handle = C_NULL_PTR
    err_stat = CD_C_OK
    ! Every CableDyn solve is a small banded system: keep a linked OpenBLAS from
    ! starting one worker (and one committed work buffer) per logical CPU. A run-time
    ! loaded BLAS library that cannot be loaded is reported by the Init calls, which
    ! have a handle to carry the diagnostic (see blas_runtime_ready).
    stat = INT(cabledyn_blas_threads_once())
    ALLOCATE (h, STAT=stat)
    IF (stat /= 0) THEN
      err_stat = CD_C_ALLOC_FAIL
      RETURN
    END IF
    h%magic = CD_C_HANDLE_MAGIC
    h%last_error = ''
    handle = C_LOC(h)
    CALL registry_add(handle, stat)
    IF (stat /= 0) THEN
      DEALLOCATE (h)
      handle = C_NULL_PTR
      err_stat = CD_C_ALLOC_FAIL
    END IF
  END SUBROUTINE CableDyn_Create

  SUBROUTINE CableDyn_InitDeck(handle, deck_path_ptr, deck_path_len, err_stat) BIND(C, NAME='CableDyn_InitDeck')
    !! Initialize from a supported MoorDyn-style .dat deck. This routes through the
    !! production deck parser and static-IC builder, and therefore inherits its
    !! fail-closed feature subset. Deck initialisations are serialised process-wide:
    !! concurrent calls (on distinct handles) run one after another.
    TYPE(C_PTR), VALUE :: handle, deck_path_ptr
    INTEGER(C_INT), VALUE :: deck_path_len
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    IF (.NOT. blas_runtime_ready(handle, 'CableDyn_InitDeck', err_stat)) RETURN
    CALL cabledyn_input_lock()
    CALL init_deck_serialised(handle, deck_path_ptr, deck_path_len, err_stat)
    CALL cabledyn_input_unlock()
  END SUBROUTINE CableDyn_InitDeck

  SUBROUTINE init_deck_serialised(handle, deck_path_ptr, deck_path_len, err_stat)
    !! CableDyn_InitDeck body; the caller holds the process-wide input lock.
    TYPE(C_PTR), VALUE :: handle, deck_path_ptr
    INTEGER(C_INT), VALUE :: deck_path_len
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    TYPE(CD_C_Handle), POINTER :: h
    TYPE(CD_FAST_ModuleType) :: new_fast
    CHARACTER(KIND=C_CHAR), POINTER :: path_c(:)
    CHARACTER(:), ALLOCATABLE :: deck_path
    INTEGER :: i, es, end_es, stat, n_objects, n_finite
    INTEGER, ALLOCATABLE :: line_ids(:)
    REAL(wp) :: dtm, gravity, rho_water
    LOGICAL :: has_dtm
    CHARACTER(CD_C_MSG_LEN) :: em, end_em

    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (deck_path_len < 1_C_INT .OR. .NOT. C_ASSOCIATED(deck_path_ptr)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitDeck: deck path pointer is null or length is non-positive', &
                      err_stat)
      RETURN
    END IF
    CALL C_F_POINTER(deck_path_ptr, path_c, [INT(deck_path_len)])
    ALLOCATE (CHARACTER(LEN=INT(deck_path_len)) :: deck_path, STAT=stat)
    IF (stat /= 0) THEN
      CALL set_status(h, CD_C_ALLOC_FAIL, 'CableDyn_InitDeck: deck path allocation failed', err_stat)
      RETURN
    END IF
    DO i = 1, INT(deck_path_len)
      IF (path_c(i) == C_NULL_CHAR) THEN
        IF (i == 1) THEN
          CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitDeck: deck path is empty', err_stat)
          RETURN
        END IF
        deck_path = deck_path(:i - 1)
        EXIT
      END IF
      deck_path(i:i) = path_c(i)
    END DO
    ! A deck with BODIES or RODS (standalone rules: a motionFile-free deck with dtM) runs on the
    ! coupled aggregate; every other deck keeps the point-system route and its diagnostics.
    ! A deck with finite-EI lines likewise runs on the aggregate (its Hermite cables).
    CALL CD_Deck_Query_dtM(TRIM(deck_path), dtm, has_dtm, es, em, standalone_scan=.TRUE., n_bodies_rods=n_objects, &
                           n_finite_ei=n_finite)
    IF (es == CD_DECKDRV_OK .AND. (n_objects > 0 .OR. n_finite > 0)) THEN
      CALL init_deck_aggregate(h, TRIM(deck_path), has_dtm, dtm, err_stat)
      RETURN
    END IF
    CALL CD_Init_Deck_System(TRIM(deck_path), new_fast%system, dtm, gravity, rho_water, es, em, line_ids=line_ids)
    IF (es /= CD_DECKDRV_OK) THEN
      CALL map_deck_status(es, err_stat)
      CALL store_status(h, err_stat, 'CableDyn_InitDeck: '//TRIM(em))
      CALL CD_FAST_End(new_fast, end_es, end_em)
      RETURN
    END IF
    CALL CD_FAST_End(h%fast, end_es, end_em)
    h%fast = new_fast
    CALL CD_FAST_End(new_fast, end_es, end_em)
    CALL drop_aggregate(h)
    CALL set_point_line_map(h, line_ids)
    CALL store_status(h, CD_C_OK, '')
    err_stat = CD_C_OK
  END SUBROUTINE init_deck_serialised

  SUBROUTINE init_deck_aggregate(h, deck_path, has_dtm, dtm, err_stat)
    !! Build a BODIES/RODS deck on the coupled aggregate (transactional: the handle keeps its
    !! previous model on failure). Coupled/Vessel rods and bodies are 6-DOF host nodes, which the
    !! translational C ABI cannot drive; the standalone deck rules already reject them without a
    !! motionFile, and the aggregate rejects a motionFile.
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    CHARACTER(*), INTENT(IN) :: deck_path
    LOGICAL, INTENT(IN) :: has_dtm
    REAL(wp), INTENT(IN) :: dtm
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    ! On the heap: InitDeck runs on the host's thread, whose stack size the host chooses.
    TYPE(CD_AGG_ModuleType), ALLOCATABLE :: new_agg
    INTEGER :: es, nm, stat
    REAL(wp), ALLOCATABLE :: pos(:, :), vel(:, :), acc(:, :), load(:, :)
    CHARACTER(CD_C_MSG_LEN) :: em

    IF (.NOT. (has_dtm .AND. dtm > CD_ZERO)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitDeck: a deck with finite-EI lines, BODIES or RODS '// &
                      'needs dtM, the fixed step of CableDyn_Step', err_stat)
      RETURN
    END IF
    ALLOCATE (new_agg, STAT=stat)
    IF (stat /= 0) THEN
      CALL set_status(h, CD_C_ALLOC_FAIL, 'CableDyn_InitDeck: aggregate allocation failed', err_stat)
      RETURN
    END IF
    CALL CD_AGG_Init_From_Deck(new_agg, deck_path, dtm, es, em)
    IF (es /= CD_AGG_OK) THEN
      CALL map_agg_status(es, err_stat)
      CALL store_status(h, err_stat, 'CableDyn_InitDeck: '//TRIM(em))
      CALL CD_AGG_End(new_agg, es, em)
      RETURN
    END IF
    nm = CD_AGG_NMovingPoints(new_agg, es, em)
    ! Allocate the new route's workspace before releasing the old route, so an
    ! allocation failure leaves the handle's previous model in place.
    ALLOCATE (pos(3, nm), vel(3, nm), acc(3, nm), load(3, nm), STAT=stat)
    IF (stat /= 0) THEN
      CALL CD_AGG_End(new_agg, es, em)
      CALL set_status(h, CD_C_ALLOC_FAIL, 'CableDyn_InitDeck: aggregate workspace allocation failed', err_stat)
      RETURN
    END IF
    CALL CD_FAST_End(h%fast, es, em)
    CALL drop_aggregate(h)
    CALL clear_c_workspace(h)
    CALL MOVE_ALLOC(pos, h%agg_pos)
    CALL MOVE_ALLOC(vel, h%agg_vel)
    CALL MOVE_ALLOC(acc, h%agg_acc)
    CALL MOVE_ALLOC(load, h%agg_load)
    CALL move_aggregate(new_agg, h%agg)
    CALL clear_point_line_map(h)
    h%use_agg = .TRUE.
    h%agg_dt = dtm
    h%agg_t = CD_ZERO
    CALL store_status(h, CD_C_OK, '')
    err_stat = CD_C_OK
  END SUBROUTINE init_deck_aggregate

  SUBROUTINE move_aggregate(src, dst)
    !! Take over an initialized aggregate (deep copy, then release the source).
    TYPE(CD_AGG_ModuleType), INTENT(INOUT) :: src, dst
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em
    CALL CD_AGG_End(dst, es, em)
    dst = src
    CALL CD_AGG_End(src, es, em)
  END SUBROUTINE move_aggregate

  SUBROUTINE drop_aggregate(h)
    !! Release the aggregate route of a handle (another initialization replaces it).
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em
    IF (h%use_agg .OR. CD_AGG_IsInitialized(h%agg)) CALL CD_AGG_End(h%agg, es, em)
    h%use_agg = .FALSE.
    h%agg_dt = CD_ZERO
    h%agg_t = CD_ZERO
    IF (ALLOCATED(h%agg_pos)) DEALLOCATE (h%agg_pos, h%agg_vel, h%agg_acc, h%agg_load)
  END SUBROUTINE drop_aggregate

  SUBROUTINE map_agg_status(agg_stat, err_stat)
    INTEGER, INTENT(IN) :: agg_stat
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    SELECT CASE (agg_stat)
    CASE (CD_AGG_OK)
      err_stat = CD_C_OK
    CASE (CD_AGG_BADINPUT)
      err_stat = CD_C_BAD_INPUT
    CASE (CD_AGG_NOT_INITIALIZED)
      err_stat = CD_C_NOT_INITIALIZED
    CASE (CD_AGG_ALLOCFAIL)
      err_stat = CD_C_ALLOC_FAIL
    CASE DEFAULT
      err_stat = CD_C_SOLVE_FAIL
    END SELECT
  END SUBROUTINE map_agg_status

  SUBROUTINE CableDyn_InitLine(handle, n_nodes, n_elem, q0_ptr, v0_ptr, elem_conn_ptr, l0_ptr, ea_ptr, rho_a_ptr, &
                               fixed_dofs_ptr, n_fixed, tension_only, rho_inf, err_stat) &
    BIND(C, NAME='CableDyn_InitLine')
    !! Initialize the handle from one already-meshed EI=0 line. Connectivity and
    !! fixed DOFs use CableDyn's 1-based Fortran numbering convention. Serialised
    !! process-wide with the other initialisation calls.
    TYPE(C_PTR), VALUE :: handle, q0_ptr, v0_ptr, elem_conn_ptr, l0_ptr, ea_ptr, rho_a_ptr, fixed_dofs_ptr
    INTEGER(C_INT), VALUE :: n_nodes, n_elem, n_fixed, tension_only
    REAL(C_DOUBLE), VALUE :: rho_inf
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    IF (.NOT. blas_runtime_ready(handle, 'CableDyn_InitLine', err_stat)) RETURN
    CALL cabledyn_input_lock()
    CALL init_line_serialised(handle, n_nodes, n_elem, q0_ptr, v0_ptr, elem_conn_ptr, l0_ptr, ea_ptr, rho_a_ptr, &
                              fixed_dofs_ptr, n_fixed, tension_only, rho_inf, err_stat)
    CALL cabledyn_input_unlock()
  END SUBROUTINE CableDyn_InitLine

  SUBROUTINE init_line_serialised(handle, n_nodes, n_elem, q0_ptr, v0_ptr, elem_conn_ptr, l0_ptr, ea_ptr, &
                                  rho_a_ptr, fixed_dofs_ptr, n_fixed, tension_only, rho_inf, err_stat)
    !! CableDyn_InitLine body; the caller holds the process-wide input lock.
    TYPE(C_PTR), VALUE :: handle, q0_ptr, v0_ptr, elem_conn_ptr, l0_ptr, ea_ptr, rho_a_ptr, fixed_dofs_ptr
    INTEGER(C_INT), VALUE :: n_nodes, n_elem, n_fixed, tension_only
    REAL(C_DOUBLE), VALUE :: rho_inf
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    TYPE(CD_C_Handle), POINTER :: h
    TYPE(CD_ModelType) :: model
    TYPE(GenAlphaConfig) :: cfg
    REAL(C_DOUBLE), POINTER :: q0_c(:), v0_c(:), l0_c(:), ea_c(:), rho_a_c(:)
    INTEGER(C_INT), POINTER :: elem_c(:), fixed_c(:)
    REAL(wp), ALLOCATABLE :: q0(:), v0(:), l0(:), ea(:), rho_a(:), f_ext(:)
    INTEGER, ALLOCATABLE :: elem_conn(:, :), fixed_dofs(:)
    INTEGER :: i, es, end_es, stat
    INTEGER(INT64) :: n_nodes_wide, n_elem_wide, n_fixed_wide, int_huge
    CHARACTER(CD_C_MSG_LEN) :: em, fast_em

    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (n_nodes < 2_C_INT .OR. n_elem < 1_C_INT .OR. n_fixed < 0_C_INT) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLine: invalid node, element, or fixed-DOF count', err_stat)
      RETURN
    END IF
    n_nodes_wide = INT(n_nodes, INT64)
    n_elem_wide = INT(n_elem, INT64)
    n_fixed_wide = INT(n_fixed, INT64)
    int_huge = INT(HUGE(i), INT64)
    IF (3_INT64*n_nodes_wide > int_huge .OR. 2_INT64*n_elem_wide > int_huge .OR. &
        n_fixed_wide > int_huge) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLine: input sizes exceed supported Fortran array bounds', &
                      err_stat)
      RETURN
    END IF
    IF (.NOT. (C_ASSOCIATED(q0_ptr) .AND. C_ASSOCIATED(v0_ptr) .AND. C_ASSOCIATED(elem_conn_ptr) .AND. &
               C_ASSOCIATED(l0_ptr) .AND. C_ASSOCIATED(ea_ptr) .AND. C_ASSOCIATED(rho_a_ptr))) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLine: required input pointer is null', err_stat)
      RETURN
    END IF
    IF (n_fixed > 0_C_INT .AND. .NOT. C_ASSOCIATED(fixed_dofs_ptr)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLine: fixed_dofs pointer is null with n_fixed > 0', err_stat)
      RETURN
    END IF

    CALL C_F_POINTER(q0_ptr, q0_c, [3*INT(n_nodes)])
    CALL C_F_POINTER(v0_ptr, v0_c, [3*INT(n_nodes)])
    CALL C_F_POINTER(elem_conn_ptr, elem_c, [2*INT(n_elem)])
    CALL C_F_POINTER(l0_ptr, l0_c, [INT(n_elem)])
    CALL C_F_POINTER(ea_ptr, ea_c, [INT(n_elem)])
    CALL C_F_POINTER(rho_a_ptr, rho_a_c, [INT(n_elem)])
    IF (n_fixed > 0_C_INT) CALL C_F_POINTER(fixed_dofs_ptr, fixed_c, [INT(n_fixed)])

    ALLOCATE (q0(3*INT(n_nodes)), v0(3*INT(n_nodes)), f_ext(3*INT(n_nodes)), &
              elem_conn(2, INT(n_elem)), l0(INT(n_elem)), ea(INT(n_elem)), rho_a(INT(n_elem)), &
              fixed_dofs(INT(n_fixed)), STAT=stat)
    IF (stat /= 0) THEN
      CALL set_status(h, CD_C_ALLOC_FAIL, 'CableDyn_InitLine: temporary input allocation failed', err_stat)
      RETURN
    END IF
    q0 = REAL(q0_c, wp)
    v0 = REAL(v0_c, wp)
    l0 = REAL(l0_c, wp)
    ea = REAL(ea_c, wp)
    rho_a = REAL(rho_a_c, wp)
    f_ext = CD_ZERO
    DO i = 1, INT(n_elem)
      elem_conn(:, i) = INT(elem_c(2*i - 1:2*i))
    END DO
    IF (n_fixed > 0_C_INT) fixed_dofs = INT(fixed_c)
    IF (.NOT. (rho_inf >= 0.0_C_DOUBLE .AND. rho_inf <= 1.0_C_DOUBLE)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLine: rho_inf must lie in [0, 1]', err_stat)
      RETURN
    END IF
    cfg%rho_inf = REAL(rho_inf, wp)

    CALL CD_Init_Model(model, q0, v0, elem_conn, l0, ea, rho_a, tension_only /= 0_C_INT, f_ext, fixed_dofs, &
                       cfg, es, em)
    IF (es /= CD_MODEL_OK) THEN
      CALL set_model_init_failure(h, 'CableDyn_InitLine', es, em, err_stat)
      CALL CD_End_Model(model, es, em)
      RETURN
    END IF
    CALL CD_FAST_Init_From_Models(h%fast, [model], es, fast_em)
    CALL CD_End_Model(model, end_es, em)
    CALL map_fast_status(es, err_stat)
    IF (err_stat == CD_C_OK) THEN
      CALL drop_aggregate(h)
      CALL set_point_line_map(h, [1])
      CALL store_status(h, err_stat, '')
    ELSE
      CALL store_status(h, err_stat, 'CableDyn_InitLine: '//TRIM(fast_em))
    END IF
  END SUBROUTINE init_line_serialised

  SUBROUTINE CableDyn_InitLines(handle, n_lines, n_nodes_ptr, n_elem_ptr, n_fixed_ptr, q0_ptr, v0_ptr, &
                                elem_conn_ptr, l0_ptr, ea_ptr, rho_a_ptr, fixed_dofs_ptr, &
                                coupled_map_ptr, n_coupled_map, tension_only, rho_inf, err_stat) &
    BIND(C, NAME='CableDyn_InitLines')
    !! Initialize from multiple already-meshed EI=0 lines. Per-line arrays are
    !! concatenated in line order. Connectivity and fixed DOFs are local to each
    !! line and use 1-based Fortran numbering. `coupled_map_ptr` may be null when
    !! `n_coupled_map == 0`, in which case the local coupled DOFs are not shared.
    !! Serialised process-wide with the other initialisation calls.
    TYPE(C_PTR), VALUE :: handle, n_nodes_ptr, n_elem_ptr, n_fixed_ptr, q0_ptr, v0_ptr
    TYPE(C_PTR), VALUE :: elem_conn_ptr, l0_ptr, ea_ptr, rho_a_ptr, fixed_dofs_ptr, coupled_map_ptr
    INTEGER(C_INT), VALUE :: n_lines, n_coupled_map, tension_only
    REAL(C_DOUBLE), VALUE :: rho_inf
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    IF (.NOT. blas_runtime_ready(handle, 'CableDyn_InitLines', err_stat)) RETURN
    CALL cabledyn_input_lock()
    CALL init_lines_serialised(handle, n_lines, n_nodes_ptr, n_elem_ptr, n_fixed_ptr, q0_ptr, v0_ptr, &
                               elem_conn_ptr, l0_ptr, ea_ptr, rho_a_ptr, fixed_dofs_ptr, &
                               coupled_map_ptr, n_coupled_map, tension_only, rho_inf, err_stat)
    CALL cabledyn_input_unlock()
  END SUBROUTINE CableDyn_InitLines

  SUBROUTINE init_lines_serialised(handle, n_lines, n_nodes_ptr, n_elem_ptr, n_fixed_ptr, q0_ptr, v0_ptr, &
                                   elem_conn_ptr, l0_ptr, ea_ptr, rho_a_ptr, fixed_dofs_ptr, &
                                   coupled_map_ptr, n_coupled_map, tension_only, rho_inf, err_stat)
    !! CableDyn_InitLines body; the caller holds the process-wide input lock.
    TYPE(C_PTR), VALUE :: handle, n_nodes_ptr, n_elem_ptr, n_fixed_ptr, q0_ptr, v0_ptr
    TYPE(C_PTR), VALUE :: elem_conn_ptr, l0_ptr, ea_ptr, rho_a_ptr, fixed_dofs_ptr, coupled_map_ptr
    INTEGER(C_INT), VALUE :: n_lines, n_coupled_map, tension_only
    REAL(C_DOUBLE), VALUE :: rho_inf
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    TYPE(CD_C_Handle), POINTER :: h
    INTEGER(C_INT), POINTER :: n_nodes_c(:), n_elem_c(:), n_fixed_c(:), elem_c(:), fixed_c(:), map_c(:)
    INTEGER(C_INT), ALLOCATABLE :: empty_fixed(:)
    REAL(C_DOUBLE), POINTER :: q0_c(:), v0_c(:), l0_c(:), ea_c(:), rho_a_c(:)
    TYPE(CD_ModelType), ALLOCATABLE :: models(:)
    TYPE(GenAlphaConfig) :: cfg
    INTEGER, ALLOCATABLE :: coupled_map(:)
    INTEGER :: i, es, stat, total_nodes, total_elem, total_fixed
    INTEGER :: dof_lo, dof_hi, elem_lo, elem_hi, elem_pair_lo, elem_pair_hi, fixed_lo, fixed_hi
    INTEGER(INT64) :: total_nodes_wide, total_elem_wide, total_fixed_wide
    INTEGER(INT64) :: max_nodes_wide, max_elem_wide, max_fixed_wide
    INTEGER(INT64) :: dof_off, elem_off, fixed_off, local_nodes, local_elem, local_fixed, local_dof
    INTEGER(INT64) :: dof_stop, elem_stop, elem_pair_stop, fixed_stop, int_huge
    CHARACTER(CD_C_MSG_LEN) :: em

    NULLIFY (fixed_c, map_c)
    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (n_lines < 1_C_INT .OR. .NOT. (C_ASSOCIATED(n_nodes_ptr) .AND. C_ASSOCIATED(n_elem_ptr) .AND. &
                                      C_ASSOCIATED(n_fixed_ptr))) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLines: invalid line count or null count pointer', err_stat)
      RETURN
    END IF
    CALL C_F_POINTER(n_nodes_ptr, n_nodes_c, [INT(n_lines)])
    CALL C_F_POINTER(n_elem_ptr, n_elem_c, [INT(n_lines)])
    CALL C_F_POINTER(n_fixed_ptr, n_fixed_c, [INT(n_lines)])
    IF (ANY(n_nodes_c < 2_C_INT) .OR. ANY(n_elem_c < 1_C_INT) .OR. ANY(n_fixed_c < 0_C_INT)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLines: invalid per-line node, element, or fixed-DOF count', &
                      err_stat)
      RETURN
    END IF
    total_nodes_wide = SUM(INT(n_nodes_c, INT64))
    total_elem_wide = SUM(INT(n_elem_c, INT64))
    total_fixed_wide = SUM(INT(n_fixed_c, INT64))
    int_huge = INT(HUGE(total_nodes), INT64)
    max_nodes_wide = int_huge/3_INT64
    max_elem_wide = int_huge/2_INT64
    max_fixed_wide = int_huge
    IF (total_nodes_wide > max_nodes_wide .OR. total_elem_wide > max_elem_wide .OR. &
        total_fixed_wide > max_fixed_wide) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLines: aggregate input sizes exceed supported bounds', &
                      err_stat)
      RETURN
    END IF
    total_nodes = INT(total_nodes_wide)
    total_elem = INT(total_elem_wide)
    total_fixed = INT(total_fixed_wide)
    IF (.NOT. (C_ASSOCIATED(q0_ptr) .AND. C_ASSOCIATED(v0_ptr) .AND. C_ASSOCIATED(elem_conn_ptr) .AND. &
               C_ASSOCIATED(l0_ptr) .AND. C_ASSOCIATED(ea_ptr) .AND. C_ASSOCIATED(rho_a_ptr))) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLines: required input pointer is null', err_stat)
      RETURN
    END IF
    IF (total_fixed > 0 .AND. .NOT. C_ASSOCIATED(fixed_dofs_ptr)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLines: fixed_dofs pointer is null with fixed DOFs present', &
                      err_stat)
      RETURN
    END IF
    CALL C_F_POINTER(q0_ptr, q0_c, [3*total_nodes])
    CALL C_F_POINTER(v0_ptr, v0_c, [3*total_nodes])
    CALL C_F_POINTER(elem_conn_ptr, elem_c, [2*total_elem])
    CALL C_F_POINTER(l0_ptr, l0_c, [total_elem])
    CALL C_F_POINTER(ea_ptr, ea_c, [total_elem])
    CALL C_F_POINTER(rho_a_ptr, rho_a_c, [total_elem])
    IF (total_fixed > 0) CALL C_F_POINTER(fixed_dofs_ptr, fixed_c, [total_fixed])
    IF (n_coupled_map < 0_C_INT) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLines: n_coupled_map must be non-negative', err_stat)
      RETURN
    END IF
    IF (n_coupled_map > 0_C_INT) THEN
      IF (.NOT. C_ASSOCIATED(coupled_map_ptr) .OR. &
          INT(n_coupled_map, INT64) /= 6_INT64*INT(n_lines, INT64)) THEN
        CALL set_status(h, CD_C_BAD_INPUT, &
                        'CableDyn_InitLines: coupled_map must have 6 entries per line when supplied', err_stat)
        RETURN
      END IF
      IF (ANY(n_fixed_c /= 6_C_INT)) THEN
        CALL set_status(h, CD_C_BAD_INPUT, &
                        'CableDyn_InitLines: shared coupled map requires six fixed DOFs per line', err_stat)
        RETURN
      END IF
      CALL C_F_POINTER(coupled_map_ptr, map_c, [INT(n_coupled_map)])
      ALLOCATE (coupled_map(INT(n_coupled_map)), STAT=stat)
      IF (stat /= 0) THEN
        CALL set_status(h, CD_C_ALLOC_FAIL, 'CableDyn_InitLines: coupled map allocation failed', err_stat)
        RETURN
      END IF
      coupled_map = INT(map_c)
    END IF

    ALLOCATE (models(INT(n_lines)), empty_fixed(0), STAT=stat)
    IF (stat /= 0) THEN
      CALL set_status(h, CD_C_ALLOC_FAIL, 'CableDyn_InitLines: model workspace allocation failed', err_stat)
      RETURN
    END IF
    IF (.NOT. (rho_inf >= 0.0_C_DOUBLE .AND. rho_inf <= 1.0_C_DOUBLE)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLines: rho_inf must lie in [0, 1]', err_stat)
      RETURN
    END IF
    cfg%rho_inf = REAL(rho_inf, wp)
    dof_off = 0_INT64
    elem_off = 0_INT64
    fixed_off = 0_INT64
    DO i = 1, INT(n_lines)
      local_nodes = INT(n_nodes_c(i), INT64)
      local_elem = INT(n_elem_c(i), INT64)
      local_fixed = INT(n_fixed_c(i), INT64)
      local_dof = 3_INT64*local_nodes
      dof_stop = dof_off + local_dof
      elem_stop = elem_off + local_elem
      elem_pair_stop = 2_INT64*elem_stop
      fixed_stop = fixed_off + local_fixed
      IF (dof_stop > int_huge .OR. elem_stop > int_huge .OR. elem_pair_stop > int_huge .OR. &
          fixed_stop > int_huge) THEN
        CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_InitLines: per-line slice exceeds supported bounds', err_stat)
        CALL end_temp_models(models)
        RETURN
      END IF
      dof_lo = INT(dof_off + 1_INT64)
      dof_hi = INT(dof_stop)
      elem_lo = INT(elem_off + 1_INT64)
      elem_hi = INT(elem_stop)
      elem_pair_lo = INT(2_INT64*elem_off + 1_INT64)
      elem_pair_hi = INT(elem_pair_stop)
      fixed_lo = INT(fixed_off + 1_INT64)
      fixed_hi = INT(fixed_stop)
      IF (n_fixed_c(i) > 0_C_INT) THEN
        CALL init_one_raw_model(models(i), q0_c(dof_lo:dof_hi), v0_c(dof_lo:dof_hi), &
                                elem_c(elem_pair_lo:elem_pair_hi), &
                                l0_c(elem_lo:elem_hi), ea_c(elem_lo:elem_hi), rho_a_c(elem_lo:elem_hi), &
                                fixed_c(fixed_lo:fixed_hi), &
                                tension_only /= 0_C_INT, cfg, es, em)
      ELSE
        CALL init_one_raw_model(models(i), q0_c(dof_lo:dof_hi), v0_c(dof_lo:dof_hi), &
                                elem_c(elem_pair_lo:elem_pair_hi), &
                                l0_c(elem_lo:elem_hi), ea_c(elem_lo:elem_hi), rho_a_c(elem_lo:elem_hi), &
                                empty_fixed, tension_only /= 0_C_INT, cfg, es, em)
      END IF
      IF (es /= CD_MODEL_OK) THEN
        CALL set_model_init_failure(h, 'CableDyn_InitLines', es, em, err_stat)
        CALL end_temp_models(models)
        RETURN
      END IF
      dof_off = dof_stop
      elem_off = elem_stop
      fixed_off = fixed_stop
    END DO
    IF (ALLOCATED(coupled_map)) THEN
      CALL CD_FAST_Init_From_Models(h%fast, models, es, em, coupled_dof_map=coupled_map)
    ELSE
      CALL CD_FAST_Init_From_Models(h%fast, models, es, em)
    END IF
    CALL end_temp_models(models)
    CALL map_fast_status(es, err_stat)
    IF (err_stat == CD_C_OK) THEN
      CALL drop_aggregate(h)
      CALL set_point_line_map(h, [(i, i=1, INT(n_lines))])
      CALL store_status(h, err_stat, '')
    ELSE
      CALL store_status(h, err_stat, 'CableDyn_InitLines: '//TRIM(em))
    END IF
  END SUBROUTINE init_lines_serialised

  SUBROUTINE CableDyn_UpdateStates(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, err_stat) &
    BIND(C, NAME='CableDyn_UpdateStates')
    !! Update coupled kinematics without advancing time.
    TYPE(C_PTR), VALUE :: handle, q_ptr, v_ptr, a_ptr
    INTEGER(C_INT), VALUE :: n_coupled_dof
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    REAL(C_DOUBLE), POINTER :: q_c(:), v_c(:), a_c(:)
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em

    CALL coupled_inputs(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, h, q_c, v_c, a_c, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      IF (n_coupled_dof > 0_C_INT) THEN
        h%agg_pos = RESHAPE(REAL(q_c, wp), SHAPE(h%agg_pos))
        h%agg_vel = RESHAPE(REAL(v_c, wp), SHAPE(h%agg_vel))
        h%agg_acc = RESHAPE(REAL(a_c, wp), SHAPE(h%agg_acc))
      END IF
      CALL CD_AGG_UpdateStates_Moving(h%agg, h%agg_pos, h%agg_vel, h%agg_acc, es, em)
      CALL map_agg_status(es, err_stat)
      CALL store_status(h, err_stat, em)
      RETURN
    END IF
    CALL ensure_c_vector_workspace(h, INT(n_coupled_dof), need_load=.FALSE., err_stat=err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (n_coupled_dof > 0_C_INT) THEN
      h%q_work(1:INT(n_coupled_dof)) = REAL(q_c, wp)
      h%v_work(1:INT(n_coupled_dof)) = REAL(v_c, wp)
      h%a_work(1:INT(n_coupled_dof)) = REAL(a_c, wp)
    END IF
    CALL CD_FAST_UpdateStates(h%fast, h%q_work(1:INT(n_coupled_dof)), h%v_work(1:INT(n_coupled_dof)), &
                              h%a_work(1:INT(n_coupled_dof)), es, em)
    CALL map_fast_status(es, err_stat)
    CALL store_status(h, err_stat, em)
  END SUBROUTINE CableDyn_UpdateStates

  SUBROUTINE CableDyn_UpdatePointFluidFields(handle, fluid_velocity_ptr, fluid_acceleration_ptr, waterline_z_ptr, &
                                             n_point, fluid_density, err_stat) &
    BIND(C, NAME='CableDyn_UpdatePointFluidFields')
    !! Update externally supplied lumped point-fluid kinematics on a point-initialized handle.
    TYPE(C_PTR), VALUE :: handle, fluid_velocity_ptr, fluid_acceleration_ptr, waterline_z_ptr
    INTEGER(C_INT), VALUE :: n_point
    REAL(C_DOUBLE), VALUE :: fluid_density
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    TYPE(CD_C_Handle), POINTER :: h
    REAL(C_DOUBLE), POINTER :: fluid_velocity_c(:), fluid_acceleration_c(:), waterline_z_c(:)
    INTEGER :: es, expected_point, p
    CHARACTER(CD_C_MSG_LEN) :: em

    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_UpdatePointFluidFields: point fluid fields are not '// &
                      'available on a deck with BODIES or RODS (still water)', err_stat)
      RETURN
    END IF
    IF (n_point < 0_C_INT) THEN
      CALL set_status(h, CD_C_BAD_INPUT, &
                      'CableDyn_UpdatePointFluidFields: point count must be non-negative', err_stat)
      RETURN
    END IF
    ! Size the fluid-field update by the stored point count, not the coupled-DOF count: a deck
    ! with unbound (output-only) Fixed/Coupled points has fewer coupled DOFs than stored points,
    ! and the system point fluid-field update requires one column per stored point.
    expected_point = CD_FAST_NPoints(h%fast, es, em)
    IF (es /= CD_FAST_OK) THEN
      CALL map_fast_status(es, err_stat)
      CALL store_status(h, err_stat, em)
      RETURN
    END IF
    IF (INT(n_point) /= expected_point) THEN
      CALL set_status(h, CD_C_BAD_INPUT, &
                      'CableDyn_UpdatePointFluidFields: point count does not match the system point count', err_stat)
      RETURN
    END IF
    IF (n_point > 0_C_INT .AND. &
        .NOT. (C_ASSOCIATED(fluid_velocity_ptr) .AND. C_ASSOCIATED(fluid_acceleration_ptr) .AND. &
               C_ASSOCIATED(waterline_z_ptr))) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_UpdatePointFluidFields: input pointer is null', err_stat)
      RETURN
    END IF

    CALL ensure_c_fluid_workspace(h, expected_point, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (expected_point > 0) THEN
      CALL C_F_POINTER(fluid_velocity_ptr, fluid_velocity_c, [3*expected_point])
      CALL C_F_POINTER(fluid_acceleration_ptr, fluid_acceleration_c, [3*expected_point])
      CALL C_F_POINTER(waterline_z_ptr, waterline_z_c, [expected_point])
      DO p = 1, expected_point
        h%fluid_velocity_work(:, p) = REAL(fluid_velocity_c(3*p - 2:3*p), wp)
        h%fluid_acceleration_work(:, p) = REAL(fluid_acceleration_c(3*p - 2:3*p), wp)
      END DO
      h%waterline_work(1:expected_point) = REAL(waterline_z_c, wp)
    END IF
    CALL CD_FAST_UpdatePointFluidFields(h%fast, h%fluid_velocity_work(:, 1:expected_point), &
                                        h%fluid_acceleration_work(:, 1:expected_point), &
                                        h%waterline_work(1:expected_point), &
                                        REAL(fluid_density, wp), es, em)
    CALL map_fast_status(es, err_stat)
    CALL store_status(h, err_stat, em)
  END SUBROUTINE CableDyn_UpdatePointFluidFields

  SUBROUTINE CableDyn_Step(handle, dt, q_ptr, v_ptr, a_ptr, n_coupled_dof, converged, stalled, n_iter, err_stat) &
    BIND(C, NAME='CableDyn_Step')
    !! Advance one coupled time step.
    TYPE(C_PTR), VALUE :: handle, q_ptr, v_ptr, a_ptr
    REAL(C_DOUBLE), VALUE :: dt
    INTEGER(C_INT), VALUE :: n_coupled_dof
    LOGICAL(C_BOOL), INTENT(OUT) :: converged, stalled
    INTEGER(C_INT), INTENT(OUT) :: n_iter, err_stat
    TYPE(CD_C_Handle), POINTER :: h
    REAL(C_DOUBLE), POINTER :: q_c(:), v_c(:), a_c(:)
    LOGICAL :: conv_f, stall_f
    INTEGER :: es, nit
    CHARACTER(CD_C_MSG_LEN) :: em

    converged = .FALSE._C_BOOL
    stalled = .FALSE._C_BOOL
    n_iter = 0_C_INT
    CALL coupled_inputs(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, h, q_c, v_c, a_c, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      IF (n_coupled_dof > 0_C_INT) THEN
        h%agg_pos = RESHAPE(REAL(q_c, wp), SHAPE(h%agg_pos))
        h%agg_vel = RESHAPE(REAL(v_c, wp), SHAPE(h%agg_vel))
        h%agg_acc = RESHAPE(REAL(a_c, wp), SHAPE(h%agg_acc))
      END IF
      CALL CD_AGG_Step_Moving(h%agg, REAL(dt, wp), h%agg_pos, h%agg_vel, h%agg_acc, conv_f, stall_f, nit, es, em, &
                              t_committed=h%agg_t + REAL(dt, wp))
      IF (es == CD_AGG_OK) h%agg_t = h%agg_t + REAL(dt, wp)
      converged = LOGICAL(conv_f, C_BOOL)
      stalled = LOGICAL(stall_f, C_BOOL)
      n_iter = INT(nit, C_INT)
      CALL map_agg_status(es, err_stat)
      CALL store_status(h, err_stat, em)
      CALL clear_failed_step_flags()
      RETURN
    END IF
    CALL ensure_c_vector_workspace(h, INT(n_coupled_dof), need_load=.FALSE., err_stat=err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (n_coupled_dof > 0_C_INT) THEN
      h%q_work(1:INT(n_coupled_dof)) = REAL(q_c, wp)
      h%v_work(1:INT(n_coupled_dof)) = REAL(v_c, wp)
      h%a_work(1:INT(n_coupled_dof)) = REAL(a_c, wp)
    END IF
    CALL CD_FAST_Step(h%fast, REAL(dt, wp), h%q_work(1:INT(n_coupled_dof)), h%v_work(1:INT(n_coupled_dof)), &
                      h%a_work(1:INT(n_coupled_dof)), conv_f, stall_f, nit, es, em)
    converged = LOGICAL(conv_f, C_BOOL)
    stalled = LOGICAL(stall_f, C_BOOL)
    n_iter = INT(nit, C_INT)
    CALL map_fast_status(es, err_stat)
    CALL store_status(h, err_stat, em)
    CALL clear_failed_step_flags()

  CONTAINS

    SUBROUTINE clear_failed_step_flags()
      !! A step that did not commit reports converged = stalled = false and n_iter = 0 on
      !! every route (the documented hard-failure contract).
      IF (err_stat /= CD_C_OK) THEN
        converged = .FALSE._C_BOOL
        stalled = .FALSE._C_BOOL
        n_iter = 0_C_INT
      END IF
    END SUBROUTINE clear_failed_step_flags
  END SUBROUTINE CableDyn_Step

  SUBROUTINE CableDyn_CalcOutput(handle, loads_ptr, n_coupled_dof, err_stat) BIND(C, NAME='CableDyn_CalcOutput')
    !! Copy coupled reaction loads into a C-owned array.
    TYPE(C_PTR), VALUE :: handle, loads_ptr
    INTEGER(C_INT), VALUE :: n_coupled_dof
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    REAL(C_DOUBLE), POINTER :: loads_c(:)
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em

    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (.NOT. valid_coupled_dof_count(h, n_coupled_dof, err_stat)) RETURN
    IF (n_coupled_dof > 0_C_INT .AND. .NOT. C_ASSOCIATED(loads_ptr)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_CalcOutput: coupled_loads pointer is null', err_stat)
      RETURN
    END IF
    CALL ensure_c_vector_workspace(h, INT(n_coupled_dof), need_load=.TRUE., err_stat=err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (n_coupled_dof > 0_C_INT) THEN
      CALL C_F_POINTER(loads_ptr, loads_c, [INT(n_coupled_dof)])
      loads_c = 0.0_C_DOUBLE
    END IF
    IF (h%use_agg) THEN
      CALL CD_AGG_CalcOutput(h%agg, es, em)
      IF (es == CD_AGG_OK) CALL CD_AGG_GetMovingPointMesh(h%agg, h%agg_pos, h%agg_vel, h%agg_acc, h%agg_load, es, em)
      IF (es == CD_AGG_OK .AND. n_coupled_dof > 0_C_INT) &
        loads_c = REAL(RESHAPE(h%agg_load, [INT(n_coupled_dof)]), C_DOUBLE)
      CALL map_agg_status(es, err_stat)
      CALL store_status(h, err_stat, em)
      RETURN
    END IF
    CALL CD_FAST_CalcOutput(h%fast, h%loads_work(1:INT(n_coupled_dof)), es, em)
    IF (es == CD_FAST_OK .AND. n_coupled_dof > 0_C_INT) &
      loads_c = REAL(h%loads_work(1:INT(n_coupled_dof)), C_DOUBLE)
    CALL map_fast_status(es, err_stat)
    CALL store_status(h, err_stat, em)
  END SUBROUTINE CableDyn_CalcOutput

  SUBROUTINE CableDyn_CalcOutputDerivatives(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, eps_fd, loads_ptr, &
                                            dload_dq_ptr, dload_dv_ptr, dload_da_ptr, added_mass_ptr, err_stat) &
    BIND(C, NAME='CableDyn_CalcOutputDerivatives')
    !! Return coupled loads plus reduced shell load derivatives for tight
    !! coupling. The q/v/a blocks are assembled at the Fortran model/system
    !! boundary; added mass is the negative acceleration derivative. Matrices are
    !! flat column-major arrays of size n*n.
    TYPE(C_PTR), VALUE :: handle, q_ptr, v_ptr, a_ptr, loads_ptr
    TYPE(C_PTR), VALUE :: dload_dq_ptr, dload_dv_ptr, dload_da_ptr, added_mass_ptr
    INTEGER(C_INT), VALUE :: n_coupled_dof
    REAL(C_DOUBLE), VALUE :: eps_fd
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    REAL(C_DOUBLE), POINTER :: q_c(:), v_c(:), a_c(:), loads_c(:)
    REAL(C_DOUBLE), POINTER :: jq_c(:), jv_c(:), ja_c(:), ma_c(:)
    INTEGER :: es, i, j, n
    CHARACTER(CD_C_MSG_LEN) :: em

    NULLIFY (q_c, v_c, a_c, loads_c, jq_c, jv_c, ja_c, ma_c)
    CALL coupled_inputs(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, h, q_c, v_c, a_c, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_CalcOutputDerivatives: load derivatives are not available on '// &
                      'a deck with BODIES or RODS', err_stat)
      RETURN
    END IF
    IF (n_coupled_dof > 0_C_INT .AND. &
        .NOT. (C_ASSOCIATED(loads_ptr) .AND. C_ASSOCIATED(dload_dq_ptr) .AND. C_ASSOCIATED(dload_dv_ptr) .AND. &
               C_ASSOCIATED(dload_da_ptr) .AND. C_ASSOCIATED(added_mass_ptr))) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_CalcOutputDerivatives: output pointer is null', err_stat)
      RETURN
    END IF
    n = INT(n_coupled_dof)
    IF (n > 0) THEN
      ! Guard the flat n*n shapes in INT64 before forming them: a corrupted/huge
      ! n_coupled_dof would otherwise wrap the default-integer n*n and map the C
      ! buffers at the wrong size. Reject before any C_F_POINTER / workspace alloc.
      IF (INT(n, INT64)*INT(n, INT64) > INT(HUGE(n), INT64)) THEN
        CALL set_status(h, CD_C_BAD_INPUT, &
                        'CableDyn_CalcOutputDerivatives: n_coupled_dof too large; n*n overflows the buffer index', &
                        err_stat)
        RETURN
      END IF
      CALL C_F_POINTER(loads_ptr, loads_c, [n])
      CALL C_F_POINTER(dload_dq_ptr, jq_c, [n*n])
      CALL C_F_POINTER(dload_dv_ptr, jv_c, [n*n])
      CALL C_F_POINTER(dload_da_ptr, ja_c, [n*n])
      CALL C_F_POINTER(added_mass_ptr, ma_c, [n*n])
      loads_c = 0.0_C_DOUBLE
      jq_c = 0.0_C_DOUBLE
      jv_c = 0.0_C_DOUBLE
      ja_c = 0.0_C_DOUBLE
      ma_c = 0.0_C_DOUBLE
    END IF
    CALL ensure_c_derivative_workspace(h, n, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (n > 0) THEN
      h%q_work(1:n) = REAL(q_c, wp)
      h%v_work(1:n) = REAL(v_c, wp)
      h%a_work(1:n) = REAL(a_c, wp)
    END IF
    CALL CD_FAST_CalcOutputDerivatives(h%fast, h%q_work(1:n), h%v_work(1:n), h%a_work(1:n), &
                                       REAL(eps_fd, wp), h%loads_work(1:n), h%jq_work(1:n, 1:n), &
                                       h%jv_work(1:n, 1:n), h%ja_work(1:n, 1:n), h%madd_work(1:n, 1:n), es, em)
    IF (es == CD_FAST_OK .AND. n > 0) THEN
      loads_c = REAL(h%loads_work(1:n), C_DOUBLE)
      DO j = 1, n
        DO i = 1, n
          jq_c(i + (j - 1)*n) = REAL(h%jq_work(i, j), C_DOUBLE)
          jv_c(i + (j - 1)*n) = REAL(h%jv_work(i, j), C_DOUBLE)
          ja_c(i + (j - 1)*n) = REAL(h%ja_work(i, j), C_DOUBLE)
          ma_c(i + (j - 1)*n) = REAL(h%madd_work(i, j), C_DOUBLE)
        END DO
      END DO
    END IF
    CALL map_fast_status(es, err_stat)
    CALL store_status(h, err_stat, em)
  END SUBROUTINE CableDyn_CalcOutputDerivatives

  SUBROUTINE CableDyn_GetCoupledMotion(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, err_stat) &
    BIND(C, NAME='CableDyn_GetCoupledMotion')
    !! Copy current coupled kinematics into C-owned arrays.
    TYPE(C_PTR), VALUE :: handle, q_ptr, v_ptr, a_ptr
    INTEGER(C_INT), VALUE :: n_coupled_dof
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    REAL(C_DOUBLE), POINTER :: q_c(:), v_c(:), a_c(:)
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em

    CALL coupled_outputs(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, h, q_c, v_c, a_c, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      IF (n_coupled_dof > 0_C_INT) THEN
        q_c = 0.0_C_DOUBLE
        v_c = 0.0_C_DOUBLE
        a_c = 0.0_C_DOUBLE
      END IF
      CALL CD_AGG_GetMovingPointMesh(h%agg, h%agg_pos, h%agg_vel, h%agg_acc, h%agg_load, es, em)
      IF (es == CD_AGG_OK .AND. n_coupled_dof > 0_C_INT) THEN
        q_c = REAL(RESHAPE(h%agg_pos, [INT(n_coupled_dof)]), C_DOUBLE)
        v_c = REAL(RESHAPE(h%agg_vel, [INT(n_coupled_dof)]), C_DOUBLE)
        a_c = REAL(RESHAPE(h%agg_acc, [INT(n_coupled_dof)]), C_DOUBLE)
      END IF
      CALL map_agg_status(es, err_stat)
      CALL store_status(h, err_stat, em)
      RETURN
    END IF
    CALL ensure_c_vector_workspace(h, INT(n_coupled_dof), need_load=.FALSE., err_stat=err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (n_coupled_dof > 0_C_INT) THEN
      q_c = 0.0_C_DOUBLE
      v_c = 0.0_C_DOUBLE
      a_c = 0.0_C_DOUBLE
    END IF
    CALL CD_FAST_GetCoupledMotion(h%fast, h%q_work(1:INT(n_coupled_dof)), h%v_work(1:INT(n_coupled_dof)), &
                                  h%a_work(1:INT(n_coupled_dof)), es, em)
    IF (es == CD_FAST_OK .AND. n_coupled_dof > 0_C_INT) THEN
      q_c = REAL(h%q_work(1:INT(n_coupled_dof)), C_DOUBLE)
      v_c = REAL(h%v_work(1:INT(n_coupled_dof)), C_DOUBLE)
      a_c = REAL(h%a_work(1:INT(n_coupled_dof)), C_DOUBLE)
    END IF
    CALL map_fast_status(es, err_stat)
    CALL store_status(h, err_stat, em)
  END SUBROUTINE CableDyn_GetCoupledMotion

  SUBROUTINE CableDyn_GetLastError(handle, message_ptr, message_len) BIND(C, NAME='CableDyn_GetLastError')
    !! Copy the last diagnostic message stored on the opaque handle into a
    !! C-owned character buffer. The buffer is always null-terminated when
    !! `message_len > 0` and `message_ptr` is associated.
    TYPE(C_PTR), VALUE :: handle, message_ptr
    INTEGER(C_INT), VALUE :: message_len

    TYPE(CD_C_Handle), POINTER :: h
    CHARACTER(KIND=C_CHAR), POINTER :: message_c(:)
    CHARACTER(CD_C_MSG_LEN) :: msg

    IF (message_len < 1_C_INT .OR. .NOT. C_ASSOCIATED(message_ptr)) RETURN
    CALL C_F_POINTER(message_ptr, message_c, [INT(message_len)])
    message_c = C_NULL_CHAR
    IF (validate_handle(handle, h) /= CD_C_OK) THEN
      msg = 'CableDyn C API: bad handle'
    ELSE
      msg = h%last_error
    END IF
    CALL copy_message_to_c(TRIM(msg), message_c)
  END SUBROUTINE CableDyn_GetLastError

  SUBROUTINE CableDyn_Close(handle, err_stat) BIND(C, NAME='CableDyn_Close')
    !! Release a CableDyn C handle. A null handle is an idempotent no-op
    !! (CD_C_OK). A non-null pointer that is not a live CableDyn handle (stale,
    !! double-closed, or foreign) is left untouched and reported CD_C_BAD_HANDLE --
    !! it is never dereferenced or deallocated. Only a validated live handle is
    !! de-registered, finalized, deallocated, and nulled.
    !!
    !! Closing a live handle also releases the OpenMP thread pool of the calling
    !! thread (outside a parallel region only; the next parallel region starts a new
    !! one). A pool left to the runtime is torn down when its host thread exits, and
    !! its threads then end detached: a process that exits while one of them is
    !! still ending can hang in the OpenMP runtime's exit handler (MinGW-w64 libgomp
    !! with winpthreads). Joining the pool here means a host thread that closes its
    !! handles before it ends leaves no such thread behind.
!$  USE OMP_LIB, ONLY: omp_pause_resource_all, omp_pause_soft
    TYPE(C_PTR), INTENT(INOUT) :: handle
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em
!$  INTEGER :: omp_stat

    err_stat = CD_C_OK
    IF (.NOT. C_ASSOCIATED(handle)) RETURN
    ! Atomically claim ownership: exactly one concurrent Close of this (or a copied)
    ! handle removes it from the registry and proceeds; the rest get .FALSE. and must
    ! not free. A foreign/stale/double-closed pointer also fails the claim.
    IF (.NOT. registry_claim(handle)) THEN
      err_stat = CD_C_BAD_HANDLE
      RETURN
    END IF
    ! We now exclusively own the target: the registry guaranteed it was a live handle
    ! we created, and the claim removed it so no other Close can reach it.
    CALL C_F_POINTER(handle, h)
    h%magic = 0_INT64
    CALL CD_FAST_End(h%fast, es, em)
    CALL drop_aggregate(h)
    CALL clear_c_workspace(h)
    DEALLOCATE (h)
    handle = C_NULL_PTR
    ! Nonzero inside an active parallel region, where the pool is left alone.
!$  omp_stat = omp_pause_resource_all(omp_pause_soft)
  END SUBROUTINE CableDyn_Close

  FUNCTION CableDyn_IsInitialized(handle) RESULT(is_initialized) BIND(C, NAME='CableDyn_IsInitialized')
    !! Return whether the handle owns an initialized CableDyn module.
    TYPE(C_PTR), VALUE :: handle
    LOGICAL(C_BOOL) :: is_initialized
    TYPE(CD_C_Handle), POINTER :: h

    is_initialized = .FALSE._C_BOOL
    IF (validate_handle(handle, h) /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      is_initialized = LOGICAL(CD_AGG_IsInitialized(h%agg), C_BOOL)
    ELSE
      is_initialized = LOGICAL(CD_FAST_IsInitialized(h%fast), C_BOOL)
    END IF
  END FUNCTION CableDyn_IsInitialized

  FUNCTION CableDyn_NCoupledDOF(handle, err_stat) RESULT(n) BIND(C, NAME='CableDyn_NCoupledDOF')
    !! Return aggregate coupled DOF count, or zero for a null/uninitialized handle.
    TYPE(C_PTR), VALUE :: handle
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    INTEGER(C_INT) :: n
    TYPE(CD_C_Handle), POINTER :: h
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em

    n = 0_C_INT
    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      n = INT(3*CD_AGG_NMovingPoints(h%agg, es, em), C_INT)
      CALL map_agg_status(es, err_stat)
      CALL store_status(h, err_stat, em)
      RETURN
    END IF
    n = INT(CD_FAST_NCoupledDOF(h%fast, es, em), C_INT)
    CALL map_fast_status(es, err_stat)
    CALL store_status(h, err_stat, em)
  END FUNCTION CableDyn_NCoupledDOF

  FUNCTION CableDyn_NPoints(handle, err_stat) RESULT(n) BIND(C, NAME='CableDyn_NPoints')
    !! Return the stored system point count, or zero for a null/uninitialized handle. This
    !! is the column count CableDyn_UpdatePointFluidFields expects, and can exceed
    !! CableDyn_NCoupledDOF/3 when the deck has unbound (output-only) Fixed/Coupled points.
    TYPE(C_PTR), VALUE :: handle
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    INTEGER(C_INT) :: n
    TYPE(CD_C_Handle), POINTER :: h
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em

    n = 0_C_INT
    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      n = INT(h%agg%n_points, C_INT)
      CALL store_status(h, CD_C_OK, '')
      RETURN
    END IF
    n = INT(CD_FAST_NPoints(h%fast, es, em), C_INT)
    CALL map_fast_status(es, err_stat)
    CALL store_status(h, err_stat, em)
  END FUNCTION CableDyn_NPoints

  FUNCTION CableDyn_NLines(handle, err_stat) RESULT(n) BIND(C, NAME='CableDyn_NLines')
    !! Return owned line count, or zero for a null/uninitialized handle.
    TYPE(C_PTR), VALUE :: handle
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    INTEGER(C_INT) :: n
    TYPE(CD_C_Handle), POINTER :: h
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em

    n = 0_C_INT
    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      n = INT(h%agg%n_lines, C_INT)
      CALL store_status(h, CD_C_OK, '')
      RETURN
    END IF
    n = INT(CD_FAST_NLines(h%fast, es, em), C_INT)
    CALL map_fast_status(es, err_stat)
    CALL store_status(h, err_stat, em)
  END FUNCTION CableDyn_NLines

  ! ---------------------------------------------------------------------------------------
  ! In-process object queries (C ABI 1, minor extension 1). Every value is read from the
  ! committed state through the evaluator the standalone driver's output channels use
  ! (CD_Eval_Aggregate_Channel) or through the state accessors its per-line files use, so a
  ! query reports the value the matching output channel reports at the same step. Object
  ! indices are 0-based positions in the handle's inventory; ids are deck ids.
  ! ---------------------------------------------------------------------------------------

  FUNCTION CableDyn_GetAbiMinor() RESULT(minor) BIND(C, NAME='CableDyn_GetAbiMinor')
    !! Minor extension level of C ABI 1 (1 = the object queries are present).
    INTEGER(C_INT) :: minor
    minor = CD_C_ABI_MINOR
  END FUNCTION CableDyn_GetAbiMinor

  FUNCTION CableDyn_NObjects(handle, kind, err_stat) RESULT(n) BIND(C, NAME='CableDyn_NObjects')
    !! Number of objects of one kind (CD_C_OBJ_LINE/POINT/BODY/ROD) in the handle's inventory.
    TYPE(C_PTR), VALUE :: handle
    INTEGER(C_INT), VALUE :: kind
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    INTEGER(C_INT) :: n
    TYPE(CD_C_Handle), POINTER :: h
    INTEGER, ALLOCATABLE :: ids(:), nnodes(:), sub(:)

    n = 0_C_INT
    IF (.NOT. query_ready(handle, h, 'CableDyn_NObjects', err_stat)) RETURN
    CALL object_inventory(h, INT(kind), ids, nnodes, sub, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    n = INT(SIZE(ids), C_INT)
    CALL store_status(h, CD_C_OK, '')
  END FUNCTION CableDyn_NObjects

  SUBROUTINE CableDyn_GetObjectInfo(handle, kind, index, id, n_nodes, subtype, err_stat) &
    BIND(C, NAME='CableDyn_GetObjectInfo')
    !! Deck id, node count and subtype of object `index` (0-based) of one kind. Lines: node
    !! count, subtype 1 for a finite-EI (bending) line, 0 for EI=0. Points: 1 node, subtype the
    !! point kind (1 Fixed, 2 Coupled, 3 Free, 4 Connect). Bodies: 1 node, subtype 0. Rods:
    !! NumSegs + 1 nodes, subtype 0.
    TYPE(C_PTR), VALUE :: handle
    INTEGER(C_INT), VALUE :: kind, index
    INTEGER(C_INT), INTENT(OUT) :: id, n_nodes, subtype, err_stat
    TYPE(CD_C_Handle), POINTER :: h
    INTEGER, ALLOCATABLE :: ids(:), nnodes(:), sub(:)

    id = 0_C_INT
    n_nodes = 0_C_INT
    subtype = 0_C_INT
    IF (.NOT. query_ready(handle, h, 'CableDyn_GetObjectInfo', err_stat)) RETURN
    CALL object_inventory(h, INT(kind), ids, nnodes, sub, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (index < 0_C_INT .OR. INT(index) >= SIZE(ids)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_GetObjectInfo: object index out of range', err_stat)
      RETURN
    END IF
    id = INT(ids(index + 1), C_INT)
    n_nodes = INT(nnodes(index + 1), C_INT)
    subtype = INT(sub(index + 1), C_INT)
    CALL store_status(h, CD_C_OK, '')
  END SUBROUTINE CableDyn_GetObjectInfo

  SUBROUTINE CableDyn_GetLineValues(handle, index, quantity, out_ptr, n_out, err_stat) &
    BIND(C, NAME='CableDyn_GetLineValues')
    !! One quantity along line `index` (0-based) in public node order End A -> End B.
    !! quantity: CD_C_LINE_POSITION / VELOCITY / ACCELERATION (3 per node, xyz interleaved),
    !! TENSION / CURVATURE / BEND_MOMENT / DECLINATION / AZIMUTH (1 per node), or
    !! SEGMENT_TENSION (1 per segment). n_out must equal the value count exactly.
    TYPE(C_PTR), VALUE :: handle, out_ptr
    INTEGER(C_INT), VALUE :: index, quantity, n_out
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    REAL(C_DOUBLE), POINTER :: out_c(:)
    INTEGER, ALLOCATABLE :: ids(:), nnodes(:), sub(:)
    INTEGER :: nnode, nvals, line_id, node, c, es
    REAL(wp) :: val
    CHARACTER(CD_DECK_NAMELEN) :: token, stem
    CHARACTER(CD_C_MSG_LEN) :: em
    CHARACTER(1) :: field
    CHARACTER(1), PARAMETER :: xyz(3) = ['x', 'y', 'z']

    IF (.NOT. query_ready(handle, h, 'CableDyn_GetLineValues', err_stat)) RETURN
    CALL object_inventory(h, INT(CD_C_OBJ_LINE), ids, nnodes, sub, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (index < 0_C_INT .OR. INT(index) >= SIZE(ids)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_GetLineValues: line index out of range', err_stat)
      RETURN
    END IF
    line_id = ids(index + 1)
    nnode = nnodes(index + 1)
    field = 'p'
    SELECT CASE (quantity)
    CASE (CD_C_LINE_POSITION, CD_C_LINE_VELOCITY, CD_C_LINE_ACCELERATION)
      nvals = 3*nnode
      IF (quantity == CD_C_LINE_VELOCITY) field = 'v'
      IF (quantity == CD_C_LINE_ACCELERATION) field = 'a'
    CASE (CD_C_LINE_TENSION, CD_C_LINE_CURVATURE, CD_C_LINE_BEND_MOMENT, CD_C_LINE_DECLINATION, &
          CD_C_LINE_AZIMUTH)
      nvals = nnode
    CASE (CD_C_LINE_SEGMENT_TENSION)
      nvals = nnode - 1
    CASE DEFAULT
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_GetLineValues: unknown quantity code', err_stat)
      RETURN
    END SELECT
    IF (INT(n_out) /= nvals .OR. .NOT. C_ASSOCIATED(out_ptr)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_GetLineValues: output buffer is null or its size does '// &
                      'not match the line', err_stat)
      RETURN
    END IF
    CALL C_F_POINTER(out_ptr, out_c, [nvals])

    IF (quantity == CD_C_LINE_SEGMENT_TENSION) THEN
      BLOCK
        REAL(wp) :: ten(nvals)
        IF (h%use_agg) THEN
          CALL CD_Eval_Aggregate_Segment_Tensions(line_id, h%agg%has_sys, h%agg%sys%fast%system, h%agg%cables, &
                                                  h%agg%line_is_cable, h%agg%line_obj_index, ten, es, em)
        ELSE
          CALL CD_Eval_Aggregate_Segment_Tensions(line_id, .TRUE., h%fast%system, h%no_cables, h%pt_is_cable, &
                                                  h%pt_obj, ten, es, em)
        END IF
        CALL map_deck_status(es, err_stat)
        IF (err_stat /= CD_C_OK) THEN
          CALL store_status(h, err_stat, 'CableDyn_GetLineValues: '//TRIM(em))
          RETURN
        END IF
        out_c = REAL(ten, C_DOUBLE)
      END BLOCK
      CALL store_status(h, CD_C_OK, '')
      RETURN
    END IF

    DO node = 1, nnode
      stem = TRIM(int_text(line_id))//'N'//TRIM(int_text(node))
      SELECT CASE (quantity)
      CASE (CD_C_LINE_POSITION, CD_C_LINE_VELOCITY, CD_C_LINE_ACCELERATION)
        DO c = 1, 3
          CALL eval_token(h, 'L'//TRIM(stem)//field//xyz(c), val, err_stat)
          IF (err_stat /= CD_C_OK) RETURN
          out_c(3*(node - 1) + c) = REAL(val, C_DOUBLE)
        END DO
        CYCLE
      CASE (CD_C_LINE_TENSION)
        token = 'Ten'//TRIM(stem)
      CASE (CD_C_LINE_CURVATURE)
        token = 'Curv'//TRIM(stem)
      CASE (CD_C_LINE_BEND_MOMENT)
        token = 'BendMom'//TRIM(stem)
      CASE (CD_C_LINE_DECLINATION)
        token = 'L'//TRIM(stem)//'Dec'
      CASE DEFAULT
        token = 'L'//TRIM(stem)//'Azi'
      END SELECT
      CALL eval_token(h, token, val, err_stat)
      IF (err_stat /= CD_C_OK) RETURN
      out_c(node) = REAL(val, C_DOUBLE)
    END DO
    CALL store_status(h, CD_C_OK, '')
  END SUBROUTINE CableDyn_GetLineValues

  SUBROUTINE CableDyn_GetPointState(handle, index, pos_ptr, vel_ptr, force_ptr, err_stat) &
    BIND(C, NAME='CableDyn_GetPointState')
    !! Position, velocity and line force (the resultant of the forces the attached lines exert
    !! on it, the Point<P>F channels) of point `index` (0-based). Each output is 3 doubles; a
    !! null pointer skips that output.
    TYPE(C_PTR), VALUE :: handle, pos_ptr, vel_ptr, force_ptr
    INTEGER(C_INT), VALUE :: index
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    REAL(C_DOUBLE), POINTER :: out_c(:)
    INTEGER, ALLOCATABLE :: ids(:), nnodes(:), sub(:)
    INTEGER :: c, pidx
    REAL(wp) :: val
    CHARACTER(1), PARAMETER :: xyz(3) = ['x', 'y', 'z']

    IF (.NOT. query_ready(handle, h, 'CableDyn_GetPointState', err_stat)) RETURN
    CALL object_inventory(h, INT(CD_C_OBJ_POINT), ids, nnodes, sub, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (index < 0_C_INT .OR. INT(index) >= SIZE(ids)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_GetPointState: point index out of range', err_stat)
      RETURN
    END IF
    pidx = INT(index) + 1
    IF (C_ASSOCIATED(pos_ptr)) THEN
      CALL C_F_POINTER(pos_ptr, out_c, [3])
      DO c = 1, 3
        CALL eval_token(h, 'Point'//TRIM(int_text(ids(pidx)))//'p'//xyz(c), val, err_stat)
        IF (err_stat /= CD_C_OK) RETURN
        out_c(c) = REAL(val, C_DOUBLE)
      END DO
    END IF
    IF (C_ASSOCIATED(vel_ptr)) THEN
      CALL C_F_POINTER(vel_ptr, out_c, [3])
      IF (h%use_agg) THEN
        out_c = REAL(h%agg%sys%fast%system%points(pidx)%v, C_DOUBLE)
      ELSE
        out_c = REAL(h%fast%system%points(pidx)%v, C_DOUBLE)
      END IF
    END IF
    IF (C_ASSOCIATED(force_ptr)) THEN
      CALL C_F_POINTER(force_ptr, out_c, [3])
      DO c = 1, 3
        CALL eval_token(h, 'Point'//TRIM(int_text(ids(pidx)))//'F'//xyz(c), val, err_stat)
        IF (err_stat /= CD_C_OK) RETURN
        out_c(c) = REAL(val, C_DOUBLE)
      END DO
    END IF
    CALL store_status(h, CD_C_OK, '')
  END SUBROUTINE CableDyn_GetPointState

  SUBROUTINE CableDyn_GetBodyState(handle, index, pose_ptr, vel_ptr, acc_ptr, wrench_ptr, err_stat) &
    BIND(C, NAME='CableDyn_GetBodyState')
    !! Rigid6 body `index` (0-based), each output 6 doubles (a null pointer skips it): pose
    !! [x, y, z, rx, ry, rz] (reference point [m], x-y'-z'' Euler angles [deg]); velocity
    !! [vx, vy, vz, wx, wy, wz] ([m/s], [deg/s]); acceleration ([m/s^2], [deg/s^2]); net external
    !! wrench [Fx, Fy, Fz, Mx, My, Mz] about the reference point ([N], [N.m]) -- the Body<N>
    !! P/R, V/RV, A/RA and F/M channels.
    TYPE(C_PTR), VALUE :: handle, pose_ptr, vel_ptr, acc_ptr, wrench_ptr
    INTEGER(C_INT), VALUE :: index
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    INTEGER, ALLOCATABLE :: ids(:), nnodes(:), sub(:)
    CHARACTER(CD_DECK_NAMELEN) :: prefix

    IF (.NOT. query_ready(handle, h, 'CableDyn_GetBodyState', err_stat)) RETURN
    CALL object_inventory(h, INT(CD_C_OBJ_BODY), ids, nnodes, sub, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (index < 0_C_INT .OR. INT(index) >= SIZE(ids)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_GetBodyState: body index out of range', err_stat)
      RETURN
    END IF
    prefix = 'Body'//TRIM(int_text(ids(index + 1)))
    CALL eval_six(h, prefix, 'P', 'R', .TRUE., pose_ptr, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    CALL eval_six(h, prefix, 'V', 'RV', .TRUE., vel_ptr, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    CALL eval_six(h, prefix, 'A', 'RA', .TRUE., acc_ptr, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    CALL eval_six(h, prefix, 'F', 'M', .TRUE., wrench_ptr, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    CALL store_status(h, CD_C_OK, '')
  END SUBROUTINE CableDyn_GetBodyState

  SUBROUTINE CableDyn_GetRodState(handle, index, nodes_ptr, n_nodes, pose_ptr, vel_ptr, wrench_ptr, err_stat) &
    BIND(C, NAME='CableDyn_GetRodState')
    !! Rod `index` (0-based); a null pointer skips that output. nodes: node positions End A
    !! (node 0) -> End B (node NumSegs), 3 doubles per node (n_nodes = NumSegs + 1). pose:
    !! [x, y, z, rx, ry, 0], End A [m] and the roll/pitch of the axis from the vertical [deg].
    !! vel: End A [vx, vy, vz, wx, wy, wz] ([m/s], [deg/s]). wrench: net [Fx, Fy, Fz, Mx, My,
    !! Mz] about End A ([N], [N.m]) -- the Rod<N>N<k>P, P/R, V/RV and F/M channels.
    TYPE(C_PTR), VALUE :: handle, nodes_ptr, pose_ptr, vel_ptr, wrench_ptr
    INTEGER(C_INT), VALUE :: index, n_nodes
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    REAL(C_DOUBLE), POINTER :: out_c(:)
    INTEGER, ALLOCATABLE :: ids(:), nnodes(:), sub(:)
    INTEGER :: node, c
    REAL(wp) :: val
    CHARACTER(CD_DECK_NAMELEN) :: prefix
    CHARACTER(1), PARAMETER :: xyz(3) = ['x', 'y', 'z']

    IF (.NOT. query_ready(handle, h, 'CableDyn_GetRodState', err_stat)) RETURN
    CALL object_inventory(h, INT(CD_C_OBJ_ROD), ids, nnodes, sub, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    IF (index < 0_C_INT .OR. INT(index) >= SIZE(ids)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_GetRodState: rod index out of range', err_stat)
      RETURN
    END IF
    prefix = 'Rod'//TRIM(int_text(ids(index + 1)))
    IF (C_ASSOCIATED(nodes_ptr)) THEN
      IF (INT(n_nodes) /= nnodes(index + 1)) THEN
        CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_GetRodState: n_nodes does not match the rod', err_stat)
        RETURN
      END IF
      CALL C_F_POINTER(nodes_ptr, out_c, [3*nnodes(index + 1)])
      DO node = 0, nnodes(index + 1) - 1
        DO c = 1, 3
          CALL eval_token(h, TRIM(prefix)//'N'//TRIM(int_text(node))//'P'//xyz(c), val, err_stat)
          IF (err_stat /= CD_C_OK) RETURN
          out_c(3*node + c) = REAL(val, C_DOUBLE)
        END DO
      END DO
    END IF
    CALL eval_six(h, prefix, 'P', 'R', .FALSE., pose_ptr, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    CALL eval_six(h, prefix, 'V', 'RV', .TRUE., vel_ptr, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    CALL eval_six(h, prefix, 'F', 'M', .TRUE., wrench_ptr, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    CALL store_status(h, CD_C_OK, '')
  END SUBROUTINE CableDyn_GetRodState

  SUBROUTINE CableDyn_EvalChannel(handle, token_ptr, token_len, value, err_stat) BIND(C, NAME='CableDyn_EvalChannel')
    !! Evaluate one OUTPUTS channel token (the deck OUTPUTS vocabulary, e.g. "FairTen1",
    !! "L2N5px", "Curv1N3", "Point4Fz", "Body1Pz", "Rod2TenA") at the committed state. A
    !! malformed token, an unknown object, or a channel family this handle cannot evaluate
    !! is CD_C_BAD_INPUT or CD_C_SOLVE_FAIL with the reason in CableDyn_GetLastError.
    !! token_len is the buffer length; an embedded NUL ends the token early. The token
    !! itself (up to the NUL) is at most 64 characters.
    TYPE(C_PTR), VALUE :: handle, token_ptr
    INTEGER(C_INT), VALUE :: token_len
    REAL(C_DOUBLE), INTENT(OUT) :: value
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    TYPE(CD_C_Handle), POINTER :: h
    CHARACTER(KIND=C_CHAR), POINTER :: tok_c(:)
    CHARACTER(CD_DECK_NAMELEN) :: token
    INTEGER :: i, n
    REAL(wp) :: val

    value = 0.0_C_DOUBLE
    IF (.NOT. query_ready(handle, h, 'CableDyn_EvalChannel', err_stat)) RETURN
    IF (token_len < 1_C_INT .OR. .NOT. C_ASSOCIATED(token_ptr)) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_EvalChannel: token pointer is null or its length is '// &
                      'not positive', err_stat)
      RETURN
    END IF
    CALL C_F_POINTER(token_ptr, tok_c, [INT(token_len)])
    token = ''
    n = INT(token_len)
    DO i = 1, INT(token_len)
      IF (tok_c(i) == C_NULL_CHAR) THEN
        n = i - 1
        EXIT
      END IF
      IF (i > CD_DECK_NAMELEN) THEN
        CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_EvalChannel: the channel token is longer than 64 '// &
                        'characters', err_stat)
        RETURN
      END IF
      token(i:i) = tok_c(i)
    END DO
    IF (n < 1) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_EvalChannel: empty channel token', err_stat)
      RETURN
    END IF
    IF (.NOT. CD_Channel_Token_Parses(token(1:n))) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn_EvalChannel: "'//token(1:n)// &
                      '" is not a supported output channel token', err_stat)
      RETURN
    END IF
    CALL eval_token(h, token(1:n), val, err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    value = REAL(val, C_DOUBLE)
    CALL store_status(h, CD_C_OK, '')
  END SUBROUTINE CableDyn_EvalChannel

  LOGICAL FUNCTION query_ready(handle, h, caller, err_stat) RESULT(ready)
    !! Validate the handle of an object query and require an initialized model.
    TYPE(C_PTR), VALUE :: handle
    TYPE(CD_C_Handle), POINTER, INTENT(OUT) :: h
    CHARACTER(*), INTENT(IN) :: caller
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    LOGICAL :: init

    ready = .FALSE.
    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (h%use_agg) THEN
      init = CD_AGG_IsInitialized(h%agg)
    ELSE
      init = CD_FAST_IsInitialized(h%fast)
    END IF
    IF (.NOT. init) THEN
      CALL set_status(h, CD_C_NOT_INITIALIZED, caller//': handle is not initialized', err_stat)
      RETURN
    END IF
    ready = .TRUE.
  END FUNCTION query_ready

  SUBROUTINE object_inventory(h, kind, ids, nnodes, sub, err_stat)
    !! The handle's objects of one kind: deck ids, node counts and subtypes. Lines are listed in
    !! ascending deck id, points in system order, bodies and rods in runtime order.
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER, INTENT(IN) :: kind
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ids(:), nnodes(:), sub(:)
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    INTEGER, ALLOCATABLE :: body_ids(:), rod_ids(:), rod_nsegs(:)

    err_stat = CD_C_OK
    SELECT CASE (kind)
    CASE (CD_C_OBJ_LINE)
      IF (h%use_agg) THEN
        CALL collect_lines(h%agg%line_is_cable, h%agg%line_obj_index)
      ELSE
        CALL collect_lines(h%pt_is_cable, h%pt_obj)
      END IF
    CASE (CD_C_OBJ_POINT)
      IF (.NOT. h%use_agg) THEN
        CALL collect_points(h%fast%system)
      ELSE IF (h%agg%has_sys) THEN
        CALL collect_points(h%agg%sys%fast%system)
      ELSE
        ALLOCATE (ids(0), nnodes(0), sub(0))
      END IF
    CASE (CD_C_OBJ_BODY)
      IF (h%use_agg .AND. h%agg%has_rigid6) THEN
        CALL CD_Rigid6_Inventory(h%agg%rigid6, body_ids, rod_ids, rod_nsegs)
        CALL MOVE_ALLOC(body_ids, ids)
      ELSE
        ALLOCATE (ids(0))
      END IF
      ALLOCATE (nnodes(SIZE(ids)), sub(SIZE(ids)))
      nnodes = 1
      sub = 0
    CASE (CD_C_OBJ_ROD)
      IF (h%use_agg .AND. h%agg%has_rigid6) THEN
        CALL CD_Rigid6_Inventory(h%agg%rigid6, body_ids, ids, rod_nsegs)
      ELSE IF (h%use_agg .AND. h%agg%has_rod) THEN
        CALL CD_Rod_Inventory(h%agg%rod, ids, rod_nsegs)
      ELSE
        ALLOCATE (ids(0), rod_nsegs(0))
      END IF
      ALLOCATE (nnodes(SIZE(ids)), sub(SIZE(ids)))
      nnodes = rod_nsegs + 1
      sub = 0
    CASE DEFAULT
      ALLOCATE (ids(0), nnodes(0), sub(0))
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn C API: unknown object kind', err_stat)
    END SELECT

  CONTAINS

    SUBROUTINE collect_lines(is_cable, obj_index)
      LOGICAL, ALLOCATABLE, INTENT(IN) :: is_cable(:)
      INTEGER, ALLOCATABLE, INTENT(IN) :: obj_index(:)
      INTEGER :: k, n, obj, es
      CHARACTER(CD_C_MSG_LEN) :: em
      IF (.NOT. (ALLOCATED(obj_index) .AND. ALLOCATED(is_cable))) THEN
        ALLOCATE (ids(0), nnodes(0), sub(0))
        RETURN
      END IF
      n = COUNT(obj_index > 0)
      ALLOCATE (ids(n), nnodes(n), sub(n))
      n = 0
      DO k = 1, SIZE(obj_index)
        obj = obj_index(k)
        IF (obj < 1) CYCLE
        n = n + 1
        ids(n) = k
        sub(n) = 0
        IF (is_cable(k)) THEN
          nnodes(n) = h%agg%cables(obj)%line%nn
          sub(n) = 1
        ELSE IF (h%use_agg) THEN
          nnodes(n) = CD_System_Line_NDOF(h%agg%sys%fast%system, obj, es, em)/3
        ELSE
          nnodes(n) = CD_System_Line_NDOF(h%fast%system, obj, es, em)/3
        END IF
      END DO
    END SUBROUTINE collect_lines

    SUBROUTINE collect_points(system)
      TYPE(CD_SystemType), INTENT(IN) :: system
      INTEGER :: n
      IF (.NOT. ALLOCATED(system%points)) THEN
        ALLOCATE (ids(0), nnodes(0), sub(0))
        RETURN
      END IF
      n = SIZE(system%points)
      ALLOCATE (ids(n), nnodes(n), sub(n))
      ids = system%points(:)%id
      nnodes = 1
      sub = system%points(:)%point_type
    END SUBROUTINE collect_points
  END SUBROUTINE object_inventory

  SUBROUTINE eval_token(h, token, val, err_stat)
    !! Evaluate one channel token on the handle's model (aggregate or point route).
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    CHARACTER(*), INTENT(IN) :: token
    REAL(wp), INTENT(OUT) :: val
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    INTEGER :: es
    CHARACTER(CD_C_MSG_LEN) :: em

    val = 0.0_wp
    IF (h%use_agg) THEN
      IF (h%agg%has_rigid6) THEN
        CALL CD_Eval_Aggregate_Channel(token, h%agg%has_sys, h%agg%sys%fast%system, h%agg%cables, &
                                       h%agg%line_is_cable, h%agg%line_obj_index, val, es, em, &
                                       cable_points=h%agg%cable_points, rigid_rt=h%agg%rigid6, ranges=h%agg%ranges)
      ELSE IF (h%agg%has_rod) THEN
        CALL CD_Eval_Aggregate_Channel(token, h%agg%has_sys, h%agg%sys%fast%system, h%agg%cables, &
                                       h%agg%line_is_cable, h%agg%line_obj_index, val, es, em, &
                                       cable_points=h%agg%cable_points, rod_rt=h%agg%rod, ranges=h%agg%ranges)
      ELSE
        CALL CD_Eval_Aggregate_Channel(token, h%agg%has_sys, h%agg%sys%fast%system, h%agg%cables, &
                                       h%agg%line_is_cable, h%agg%line_obj_index, val, es, em, &
                                       cable_points=h%agg%cable_points, ranges=h%agg%ranges)
      END IF
    ELSE
      IF (.NOT. (ALLOCATED(h%pt_obj) .AND. ALLOCATED(h%pt_is_cable))) THEN
        CALL set_status(h, CD_C_NOT_INITIALIZED, 'CableDyn C API: the handle has no line inventory', err_stat)
        RETURN
      END IF
      CALL CD_Eval_Aggregate_Channel(token, .TRUE., h%fast%system, h%no_cables, h%pt_is_cable, h%pt_obj, &
                                     val, es, em)
    END IF
    CALL map_deck_status(es, err_stat)
    IF (err_stat /= CD_C_OK) THEN
      CALL store_status(h, err_stat, TRIM(em))
    ELSE IF (.NOT. IEEE_IS_FINITE(val)) THEN
      CALL set_status(h, CD_C_SOLVE_FAIL, 'CableDyn C API: channel "'//TRIM(token)//'" is not finite', err_stat)
    END IF
  END SUBROUTINE eval_token

  SUBROUTINE eval_six(h, prefix, lin, ang, has_third_angle, out_ptr, err_stat)
    !! Fill a 6-vector [<lin>x, <lin>y, <lin>z, <ang>x, <ang>y, <ang>z] of object channels; a null
    !! pointer is a no-op. has_third_angle = .FALSE. (a rod's R, which has no yaw) sets the sixth
    !! entry to zero.
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    CHARACTER(*), INTENT(IN) :: prefix, lin, ang
    LOGICAL, INTENT(IN) :: has_third_angle
    TYPE(C_PTR), VALUE :: out_ptr
    INTEGER(C_INT), INTENT(OUT) :: err_stat
    REAL(C_DOUBLE), POINTER :: out_c(:)
    CHARACTER(1), PARAMETER :: xyz(3) = ['x', 'y', 'z']
    INTEGER :: c
    REAL(wp) :: val

    err_stat = CD_C_OK
    IF (.NOT. C_ASSOCIATED(out_ptr)) RETURN
    CALL C_F_POINTER(out_ptr, out_c, [6])
    out_c = 0.0_C_DOUBLE
    DO c = 1, 3
      CALL eval_token(h, TRIM(prefix)//lin//xyz(c), val, err_stat)
      IF (err_stat /= CD_C_OK) RETURN
      out_c(c) = REAL(val, C_DOUBLE)
    END DO
    DO c = 1, 3
      IF (c == 3 .AND. .NOT. has_third_angle) EXIT
      CALL eval_token(h, TRIM(prefix)//ang//xyz(c), val, err_stat)
      IF (err_stat /= CD_C_OK) RETURN
      out_c(3 + c) = REAL(val, C_DOUBLE)
    END DO
  END SUBROUTINE eval_six

  SUBROUTINE set_point_line_map(h, line_ids)
    !! Record the point route's line inventory: the deck-line-id -> system-line map (system
    !! line k carries deck id line_ids(k)) the shared channel evaluator resolves tokens through.
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER, INTENT(IN) :: line_ids(:)
    INTEGER :: k, maxid

    CALL clear_point_line_map(h)
    maxid = 1
    IF (SIZE(line_ids) > 0) maxid = MAX(1, MAXVAL(line_ids))
    ALLOCATE (h%pt_is_cable(maxid), h%pt_obj(maxid))
    h%pt_is_cable = .FALSE.
    h%pt_obj = 0
    DO k = 1, SIZE(line_ids)
      IF (line_ids(k) >= 1) h%pt_obj(line_ids(k)) = k
    END DO
  END SUBROUTINE set_point_line_map

  SUBROUTINE clear_point_line_map(h)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    IF (ALLOCATED(h%pt_is_cable)) DEALLOCATE (h%pt_is_cable)
    IF (ALLOCATED(h%pt_obj)) DEALLOCATE (h%pt_obj)
  END SUBROUTINE clear_point_line_map

  PURE FUNCTION int_text(i) RESULT(s)
    !! Decimal text of a non-negative id (no formatted I/O: the queries stay I/O-free).
    INTEGER, INTENT(IN) :: i
    CHARACTER(12) :: s
    INTEGER :: v, k
    CHARACTER(12) :: rev

    s = ''
    rev = ''
    v = MAX(0, i)
    k = 0
    DO
      k = k + 1
      rev(k:k) = ACHAR(IACHAR('0') + MOD(v, 10))
      v = v/10
      IF (v == 0 .OR. k == 12) EXIT
    END DO
    DO v = 1, k
      s(v:v) = rev(k - v + 1:k - v + 1)
    END DO
  END FUNCTION int_text

  INTEGER(C_INT) FUNCTION validate_handle(handle, h) RESULT(err_stat)
    !! The single handle validator used by every public entry point. A handle is
    !! accepted ONLY if (1) it is non-null, (2) it is in the live-handle registry --
    !! so a foreign, stale, or double-closed pointer is rejected before any
    !! dereference -- and (3) the target's magic word is intact. On any failure h is
    !! nullified and CD_C_BAD_HANDLE is returned.
    TYPE(C_PTR), VALUE :: handle
    TYPE(CD_C_Handle), POINTER, INTENT(OUT) :: h

    h => NULL()
    err_stat = CD_C_BAD_HANDLE
    IF (.NOT. C_ASSOCIATED(handle)) RETURN
    IF (.NOT. registry_contains(handle)) RETURN
    CALL C_F_POINTER(handle, h)
    IF (h%magic /= CD_C_HANDLE_MAGIC) THEN
      h => NULL()
      RETURN
    END IF
    err_stat = CD_C_OK
  END FUNCTION validate_handle

  SUBROUTINE registry_add(p, stat)
    !! Register a freshly created handle pointer as live. Grows the registry as
    !! needed. stat /= 0 only on an allocation failure.
    TYPE(C_PTR), INTENT(IN) :: p
    INTEGER, INTENT(OUT) :: stat
    TYPE(C_PTR), ALLOCATABLE :: tmp(:)
    stat = 0
    CALL cabledyn_registry_lock()
    IF (.NOT. ALLOCATED(g_handle_registry)) ALLOCATE (g_handle_registry(8), STAT=stat)
    IF (stat == 0) THEN
      IF (g_handle_count == SIZE(g_handle_registry)) THEN
        ALLOCATE (tmp(2*SIZE(g_handle_registry)), STAT=stat)
        IF (stat == 0) THEN
          tmp(1:g_handle_count) = g_handle_registry(1:g_handle_count)
          CALL MOVE_ALLOC(tmp, g_handle_registry)
        END IF
      END IF
    END IF
    IF (stat == 0) THEN
      g_handle_count = g_handle_count + 1
      g_handle_registry(g_handle_count) = p
    END IF
    CALL cabledyn_registry_unlock()
  END SUBROUTINE registry_add

  LOGICAL FUNCTION registry_contains(p) RESULT(found)
    !! True iff p is a currently live handle pointer. Pure lookup; no dereference.
    TYPE(C_PTR), INTENT(IN) :: p
    INTEGER :: i
    found = .FALSE.
    CALL cabledyn_registry_lock()
    IF (ALLOCATED(g_handle_registry)) THEN
      DO i = 1, g_handle_count
        IF (C_ASSOCIATED(g_handle_registry(i), p)) THEN
          found = .TRUE.
          EXIT
        END IF
      END DO
    END IF
    CALL cabledyn_registry_unlock()
  END FUNCTION registry_contains

  LOGICAL FUNCTION registry_claim(p) RESULT(claimed)
    !! Atomically check-and-remove p from the live-handle registry (swap-with-last),
    !! all inside one critical region. Returns .TRUE. to exactly ONE caller -- the one
    !! that removed it -- and .FALSE. to every other (foreign/stale/already-claimed).
    !! This is what makes Close's validate -> remove -> free safe against a concurrent
    !! Close of the same or a copied handle: only the winning claim frees the target,
    !! so there is no double-free window between validation and removal.
    TYPE(C_PTR), INTENT(IN) :: p
    INTEGER :: i
    claimed = .FALSE.
    CALL cabledyn_registry_lock()
    DO i = 1, g_handle_count
      IF (C_ASSOCIATED(g_handle_registry(i), p)) THEN
        g_handle_registry(i) = g_handle_registry(g_handle_count)
        g_handle_count = g_handle_count - 1
        claimed = .TRUE.
        EXIT
      END IF
    END DO
    CALL cabledyn_registry_unlock()
  END FUNCTION registry_claim

  SUBROUTINE coupled_inputs(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, h, q_c, v_c, a_c, err_stat)
    TYPE(C_PTR), VALUE :: handle, q_ptr, v_ptr, a_ptr
    INTEGER(C_INT), VALUE :: n_coupled_dof
    TYPE(CD_C_Handle), POINTER, INTENT(OUT) :: h
    REAL(C_DOUBLE), POINTER, INTENT(OUT) :: q_c(:), v_c(:), a_c(:)
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    CALL coupled_outputs(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, h, q_c, v_c, a_c, err_stat)
  END SUBROUTINE coupled_inputs

  SUBROUTINE coupled_outputs(handle, q_ptr, v_ptr, a_ptr, n_coupled_dof, h, q_c, v_c, a_c, err_stat)
    TYPE(C_PTR), VALUE :: handle, q_ptr, v_ptr, a_ptr
    INTEGER(C_INT), VALUE :: n_coupled_dof
    TYPE(CD_C_Handle), POINTER, INTENT(OUT) :: h
    REAL(C_DOUBLE), POINTER, INTENT(OUT) :: q_c(:), v_c(:), a_c(:)
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    err_stat = validate_handle(handle, h)
    IF (err_stat /= CD_C_OK) RETURN
    IF (.NOT. valid_coupled_dof_count(h, n_coupled_dof, err_stat)) RETURN
    IF (n_coupled_dof > 0_C_INT .AND. &
        .NOT. (C_ASSOCIATED(q_ptr) .AND. C_ASSOCIATED(v_ptr) .AND. C_ASSOCIATED(a_ptr))) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn C API: coupled kinematics pointer is null', err_stat)
      RETURN
    END IF
    IF (n_coupled_dof > 0_C_INT) THEN
      CALL C_F_POINTER(q_ptr, q_c, [INT(n_coupled_dof)])
      CALL C_F_POINTER(v_ptr, v_c, [INT(n_coupled_dof)])
      CALL C_F_POINTER(a_ptr, a_c, [INT(n_coupled_dof)])
    END IF
  END SUBROUTINE coupled_outputs

  SUBROUTINE map_fast_status(fast_stat, err_stat)
    INTEGER, INTENT(IN) :: fast_stat
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    SELECT CASE (fast_stat)
    CASE (CD_FAST_OK)
      err_stat = CD_C_OK
    CASE (CD_FAST_BADINPUT)
      err_stat = CD_C_BAD_INPUT
    CASE (CD_FAST_NOT_INITIALIZED)
      err_stat = CD_C_NOT_INITIALIZED
    CASE (CD_FAST_ALLOCFAIL)
      err_stat = CD_C_ALLOC_FAIL
    CASE DEFAULT
      err_stat = CD_C_SOLVE_FAIL
    END SELECT
  END SUBROUTINE map_fast_status

  SUBROUTINE map_deck_status(deck_stat, err_stat)
    INTEGER, INTENT(IN) :: deck_stat
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    SELECT CASE (deck_stat)
    CASE (CD_DECKDRV_OK)
      err_stat = CD_C_OK
    CASE (CD_DECKDRV_BADINPUT)
      err_stat = CD_C_BAD_INPUT
    CASE (CD_DECKDRV_SOLVEFAIL)
      err_stat = CD_C_SOLVE_FAIL
    CASE DEFAULT
      err_stat = CD_C_SOLVE_FAIL
    END SELECT
  END SUBROUTINE map_deck_status

  SUBROUTINE set_model_init_failure(h, context, model_stat, model_msg, err_stat)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    CHARACTER(*), INTENT(IN) :: context
    INTEGER, INTENT(IN) :: model_stat
    CHARACTER(*), INTENT(IN) :: model_msg
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    INTEGER(C_INT) :: c_stat

    SELECT CASE (model_stat)
    CASE (CD_MODEL_BADINPUT)
      c_stat = CD_C_BAD_INPUT
    CASE (CD_MODEL_ALLOCFAIL)
      c_stat = CD_C_ALLOC_FAIL
    CASE DEFAULT
      c_stat = CD_C_SOLVE_FAIL
    END SELECT
    CALL set_status(h, c_stat, TRIM(context)//': '//TRIM(model_msg), err_stat)
  END SUBROUTINE set_model_init_failure

  LOGICAL FUNCTION valid_coupled_dof_count(h, n_coupled_dof, err_stat) RESULT(ok)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER(C_INT), VALUE :: n_coupled_dof
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    INTEGER :: es, expected
    CHARACTER(CD_C_MSG_LEN) :: em

    ok = .FALSE.
    IF (n_coupled_dof < 0_C_INT) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn C API: coupled DOF count must be non-negative', err_stat)
      RETURN
    END IF
    IF (h%use_agg) THEN
      expected = 3*CD_AGG_NMovingPoints(h%agg, es, em)
      IF (es /= CD_AGG_OK) THEN
        CALL map_agg_status(es, err_stat)
        CALL store_status(h, err_stat, em)
        RETURN
      END IF
    ELSE
      expected = CD_FAST_NCoupledDOF(h%fast, es, em)
      IF (es /= CD_FAST_OK) THEN
        CALL map_fast_status(es, err_stat)
        CALL store_status(h, err_stat, em)
        RETURN
      END IF
    END IF
    IF (INT(n_coupled_dof) /= expected) THEN
      err_stat = CD_C_BAD_INPUT
      CALL store_status(h, err_stat, 'CableDyn C API: coupled DOF count does not match initialized handle')
      RETURN
    END IF
    err_stat = CD_C_OK
    CALL store_status(h, err_stat, '')
    ok = .TRUE.
  END FUNCTION valid_coupled_dof_count

  SUBROUTINE ensure_c_vector_workspace(h, n, need_load, err_stat)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(IN) :: need_load
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    INTEGER :: stat
    LOGICAL :: need_q_workspace, need_load_workspace

    err_stat = CD_C_OK
    IF (n < 0) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn C API: workspace size must be non-negative', err_stat)
      RETURN
    END IF
    need_q_workspace = .NOT. ALLOCATED(h%q_work)
    IF (.NOT. need_q_workspace) need_q_workspace = SIZE(h%q_work) < n
    IF (need_q_workspace) THEN
      CALL clear_c_vector_workspace(h)
      ALLOCATE (h%q_work(n), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
      ALLOCATE (h%v_work(n), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
      ALLOCATE (h%a_work(n), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
    END IF
    IF (n > 0) THEN
      h%q_work(1:n) = CD_ZERO
      h%v_work(1:n) = CD_ZERO
      h%a_work(1:n) = CD_ZERO
    END IF
    IF (need_load) THEN
      need_load_workspace = .NOT. ALLOCATED(h%loads_work)
      IF (.NOT. need_load_workspace) need_load_workspace = SIZE(h%loads_work) < n
      IF (need_load_workspace) THEN
        IF (ALLOCATED(h%loads_work)) DEALLOCATE (h%loads_work)
        ALLOCATE (h%loads_work(n), STAT=stat)
        IF (stat /= 0) THEN
          CALL fail_c_workspace_allocation(h, err_stat)
          RETURN
        END IF
      END IF
      IF (n > 0) h%loads_work(1:n) = CD_ZERO
    END IF
  END SUBROUTINE ensure_c_vector_workspace

  SUBROUTINE ensure_c_derivative_workspace(h, n, err_stat)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER, INTENT(IN) :: n
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    INTEGER :: stat
    LOGICAL :: need_derivative_workspace

    CALL ensure_c_vector_workspace(h, n, need_load=.TRUE., err_stat=err_stat)
    IF (err_stat /= CD_C_OK) RETURN
    need_derivative_workspace = .NOT. ALLOCATED(h%jq_work)
    IF (.NOT. need_derivative_workspace) THEN
      need_derivative_workspace = SIZE(h%jq_work, 1) < n .OR. SIZE(h%jq_work, 2) < n
    END IF
    IF (need_derivative_workspace) THEN
      CALL clear_c_derivative_workspace(h)
      ALLOCATE (h%jq_work(n, n), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
      ALLOCATE (h%jv_work(n, n), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
      ALLOCATE (h%ja_work(n, n), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
      ALLOCATE (h%madd_work(n, n), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
    END IF
    IF (n > 0) THEN
      h%jq_work(1:n, 1:n) = CD_ZERO
      h%jv_work(1:n, 1:n) = CD_ZERO
      h%ja_work(1:n, 1:n) = CD_ZERO
      h%madd_work(1:n, 1:n) = CD_ZERO
    END IF
  END SUBROUTINE ensure_c_derivative_workspace

  SUBROUTINE ensure_c_fluid_workspace(h, n_point, err_stat)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER, INTENT(IN) :: n_point
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    INTEGER :: stat
    LOGICAL :: need_fluid_workspace

    err_stat = CD_C_OK
    IF (n_point < 0) THEN
      CALL set_status(h, CD_C_BAD_INPUT, 'CableDyn C API: point workspace size must be non-negative', err_stat)
      RETURN
    END IF
    need_fluid_workspace = .NOT. ALLOCATED(h%fluid_velocity_work)
    IF (.NOT. need_fluid_workspace) need_fluid_workspace = SIZE(h%fluid_velocity_work, 2) < n_point
    IF (need_fluid_workspace) THEN
      CALL clear_c_fluid_workspace(h)
      ALLOCATE (h%fluid_velocity_work(3, n_point), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
      ALLOCATE (h%fluid_acceleration_work(3, n_point), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
      ALLOCATE (h%waterline_work(n_point), STAT=stat)
      IF (stat /= 0) THEN
        CALL fail_c_workspace_allocation(h, err_stat)
        RETURN
      END IF
    END IF
    IF (n_point > 0) THEN
      h%fluid_velocity_work(:, 1:n_point) = CD_ZERO
      h%fluid_acceleration_work(:, 1:n_point) = CD_ZERO
      h%waterline_work(1:n_point) = CD_ZERO
    END IF
  END SUBROUTINE ensure_c_fluid_workspace

  SUBROUTINE fail_c_workspace_allocation(h, err_stat)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    CALL clear_c_workspace(h)
    CALL set_status(h, CD_C_ALLOC_FAIL, 'CableDyn C API: workspace allocation failed', err_stat)
  END SUBROUTINE fail_c_workspace_allocation

  SUBROUTINE clear_c_workspace(h)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h

    CALL clear_c_vector_workspace(h)
    CALL clear_c_derivative_workspace(h)
    CALL clear_c_fluid_workspace(h)
  END SUBROUTINE clear_c_workspace

  SUBROUTINE clear_c_vector_workspace(h)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h

    IF (ALLOCATED(h%q_work)) DEALLOCATE (h%q_work)
    IF (ALLOCATED(h%v_work)) DEALLOCATE (h%v_work)
    IF (ALLOCATED(h%a_work)) DEALLOCATE (h%a_work)
    IF (ALLOCATED(h%loads_work)) DEALLOCATE (h%loads_work)
  END SUBROUTINE clear_c_vector_workspace

  SUBROUTINE clear_c_derivative_workspace(h)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h

    IF (ALLOCATED(h%jq_work)) DEALLOCATE (h%jq_work)
    IF (ALLOCATED(h%jv_work)) DEALLOCATE (h%jv_work)
    IF (ALLOCATED(h%ja_work)) DEALLOCATE (h%ja_work)
    IF (ALLOCATED(h%madd_work)) DEALLOCATE (h%madd_work)
  END SUBROUTINE clear_c_derivative_workspace

  SUBROUTINE clear_c_fluid_workspace(h)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h

    IF (ALLOCATED(h%fluid_velocity_work)) DEALLOCATE (h%fluid_velocity_work)
    IF (ALLOCATED(h%fluid_acceleration_work)) DEALLOCATE (h%fluid_acceleration_work)
    IF (ALLOCATED(h%waterline_work)) DEALLOCATE (h%waterline_work)
  END SUBROUTINE clear_c_fluid_workspace

  LOGICAL FUNCTION blas_runtime_ready(handle, caller, err_stat) RESULT(ready)
    !! Model initialisation gate for the LAPACK runtime: false, with CD_C_ALLOC_FAIL and
    !! the load diagnostic (library, locations tried, load error) stored on the handle,
    !! when the run-time loaded OpenBLAS is unavailable; false with the validation
    !! status for a bad handle. True otherwise, with err_stat untouched.
    TYPE(C_PTR), VALUE :: handle
    CHARACTER(*), INTENT(IN) :: caller
    INTEGER(C_INT), INTENT(INOUT) :: err_stat

    TYPE(CD_C_Handle), POINTER :: h
    INTEGER :: es
    INTEGER(C_INT) :: vs
    CHARACTER(CD_C_MSG_LEN) :: em

    ready = .FALSE.
    vs = validate_handle(handle, h)
    IF (vs /= CD_C_OK) THEN
      err_stat = vs
      RETURN
    END IF
    CALL CD_Blas_Runtime_Check(caller, es, em)
    IF (es /= CD_LINALG_OK) THEN
      CALL set_status(h, CD_C_ALLOC_FAIL, TRIM(em), err_stat)
      RETURN
    END IF
    ready = .TRUE.
  END FUNCTION blas_runtime_ready

  SUBROUTINE set_status(h, stat, msg, err_stat)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER(C_INT), INTENT(IN) :: stat
    CHARACTER(*), INTENT(IN) :: msg
    INTEGER(C_INT), INTENT(OUT) :: err_stat

    err_stat = stat
    CALL store_status(h, stat, msg)
  END SUBROUTINE set_status

  SUBROUTINE store_status(h, stat, msg)
    TYPE(CD_C_Handle), POINTER, INTENT(INOUT) :: h
    INTEGER(C_INT), INTENT(IN) :: stat
    CHARACTER(*), INTENT(IN) :: msg

    IF (stat == CD_C_OK) THEN
      h%last_error = ''
    ELSE IF (LEN_TRIM(msg) > 0) THEN
      h%last_error = TRIM(msg)
    ELSE
      h%last_error = default_status_message(stat)
    END IF
  END SUBROUTINE store_status

  FUNCTION default_status_message(stat) RESULT(msg)
    INTEGER(C_INT), INTENT(IN) :: stat
    CHARACTER(CD_C_MSG_LEN) :: msg

    SELECT CASE (stat)
    CASE (CD_C_BAD_HANDLE)
      msg = 'CableDyn C API: bad handle'
    CASE (CD_C_ALLOC_FAIL)
      msg = 'CableDyn C API: allocation failed'
    CASE (CD_C_BAD_INPUT)
      msg = 'CableDyn C API: bad input'
    CASE (CD_C_SOLVE_FAIL)
      msg = 'CableDyn C API: solve failed'
    CASE (CD_C_NOT_INITIALIZED)
      msg = 'CableDyn C API: handle is not initialized'
    CASE DEFAULT
      msg = ''
    END SELECT
  END FUNCTION default_status_message

  SUBROUTINE copy_message_to_c(msg, buffer)
    CHARACTER(*), INTENT(IN) :: msg
    CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: buffer(:)

    INTEGER :: i, ncopy

    IF (SIZE(buffer) < 1) RETURN
    buffer = C_NULL_CHAR
    ncopy = MIN(LEN_TRIM(msg), SIZE(buffer) - 1)
    DO i = 1, ncopy
      buffer(i) = msg(i:i)
    END DO
    buffer(ncopy + 1) = C_NULL_CHAR
  END SUBROUTINE copy_message_to_c

  SUBROUTINE init_one_raw_model(model, q0_c, v0_c, elem_c, l0_c, ea_c, rho_a_c, fixed_c, tension_only, cfg, &
                                ErrStat, ErrMsg)
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(C_DOUBLE), INTENT(IN) :: q0_c(:), v0_c(:), l0_c(:), ea_c(:), rho_a_c(:)
    INTEGER(C_INT), INTENT(IN) :: elem_c(:), fixed_c(:)
    LOGICAL, INTENT(IN) :: tension_only
    TYPE(GenAlphaConfig), INTENT(IN) :: cfg
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp), ALLOCATABLE :: q0(:), v0(:), l0(:), ea(:), rho_a(:), f_ext(:)
    INTEGER, ALLOCATABLE :: elem_conn(:, :), fixed_dofs(:)
    INTEGER :: e, n_elem, stat

    n_elem = SIZE(l0_c)
    ALLOCATE (q0(SIZE(q0_c)), v0(SIZE(v0_c)), f_ext(SIZE(q0_c)), &
              l0(n_elem), ea(n_elem), rho_a(n_elem), elem_conn(2, n_elem), fixed_dofs(SIZE(fixed_c)), STAT=stat)
    IF (stat /= 0) THEN
      ErrStat = CD_MODEL_ALLOCFAIL
      ErrMsg = 'CableDyn_CAPI: raw-model temporary allocation failed'
      RETURN
    END IF
    q0 = REAL(q0_c, wp)
    v0 = REAL(v0_c, wp)
    l0 = REAL(l0_c, wp)
    ea = REAL(ea_c, wp)
    rho_a = REAL(rho_a_c, wp)
    f_ext = CD_ZERO
    DO e = 1, n_elem
      elem_conn(:, e) = INT(elem_c(2*e - 1:2*e))
    END DO
    fixed_dofs = INT(fixed_c)
    CALL CD_Init_Model(model, q0, v0, elem_conn, l0, ea, rho_a, tension_only, f_ext, fixed_dofs, cfg, &
                       ErrStat, ErrMsg)
  END SUBROUTINE init_one_raw_model

  SUBROUTINE end_temp_models(models)
    TYPE(CD_ModelType), INTENT(INOUT) :: models(:)
    INTEGER :: i, es
    CHARACTER(120) :: em

    DO i = 1, SIZE(models)
      CALL CD_End_Model(models(i), es, em)
    END DO
  END SUBROUTINE end_temp_models

END MODULE CableDyn_CAPI
