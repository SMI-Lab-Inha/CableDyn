!**********************************************************************************************************************************
! SPDX-License-Identifier: Apache-2.0
! LICENSING
! Copyright (C) 2026 Jae Hoon Seo, SMI Lab, Inha University
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
MODULE CableDyn
   !! CableDyn: a first-class mooring and cable-dynamics module for OpenFAST
   !! (CompMooring = 5). CableDyn is an independent module, not a MoorDyn variant: the
   !! input file is CableDyn's, the line object follows the OrcaFlex model strictly
   !! (a LINE runs End A -> End B through ordered SECTIONS, each with its own line
   !! type, length, and mesh), the initial condition is a true Newton static
   !! equilibrium (no dynamic relaxation, no relaxation tuning), and time stepping
   !! is implicit generalised-alpha on the banded production core (one implicit
   !! step per coupling step; the module never changes the glue's DT). The deck's
   !! OPTIONS section carries CableDyn's own static- and dynamic-solver settings.
   !! Output channels follow the OrcaFlex vocabulary: the deck OUTPUTS list flows to the
   !! OpenFAST WriteOutput surface (per-line and per-node tensions, curvature, bend moment,
   !! the fairlead/anchor angle, node positions/velocities/accelerations) for both the EI=0
   !! mooring lines and the finite-EI cables.
   !!
   !! Exports the FAST-modularization-framework surface the glue couples against
   !! (CD_Init / CD_UpdateStates / CD_CalcOutput / CD_CalcContStateDeriv / CD_End + the
   !! four Jacobian routines) over the registry-generated CableDyn_Types.
   !!
   !! Coexists with stock MoorDyn in one binary: CompMooring = 3 selects MoorDyn,
   !! CompMooring = 5 selects this module (an A/B comparison is a one-digit switch).
   !! For migration, decks written in the MoorDyn v2 format also parse: each stock
   !! construct is routed into the corresponding CableDyn feature (bodies and rods,
   !! viscoelastic and Syrope ropes, line failures, active tension through the CONTROL
   !! section, SeaState water kinematics), or, where CableDyn has no counterpart (for
   !! example a CoupledPinned body or a Pinned rod on a line deck), fails closed with a
   !! named error.
   !!
   !! Supported coupled surface (each verified against stock MoorDyn or a named gate):
   !! Fixed/Coupled(Vessel) points, EI=0 mooring lines (BA damping + the Morison set),
   !! finite-EI power cables, penalty seabed, nonzero initial platform displacement
   !! (PtfmInit), the module's own mooring step (supercycling + zero-order-hold loads),
   !! glue correction iterations (NumCrctn > 0: snapshot-rewind re-advance), SeaState
   !! wave/current kinematics sampled at the line nodes (single-turbine only; see the
   !! farm scope at the use_extf guard), checkpoint restart (rebuild from the packed
   !! deck + state mirror, accepted on committed mooring-interval boundaries),
   !! linearization (quasi-static dYdu only -- no dynamic state linearization; the
   !! module registers zero continuous states and CD_CalcContStateDeriv stays fatal
   !! by design), and FAST.Farm farm-level moorings (FarmSize > 0: Turbine<J> decks,
   !! per-turbine coupled meshes, shared lines; still-water hydro; farm restart and
   !! farm linearization fail closed with named messages). Checkpoint restart and
   !! linearization re-parse a temporary copy of the deck, written next to the
   !! MooringFile (so its side files resolve as in the original run) and removed after
   !! use.
   !!
   !! State handling: the authoritative states live in the module-level instance
   !! (keyed by m%CDInst); x%states carries the mirror (per EI=0 mooring line the
   !! MoorDyn-layout interior-node velocities then positions, then per finite-EI cable
   !! its free-DOF velocities then positions) so checkpoint packing, correction-
   !! iteration rewind, and the ModVars x-variable metadata are framework-conformant.

   USE CableDyn_Types
   USE NWTC_Library
   USE CableDyn_Precision, ONLY: wp
   USE CableDyn_Hydro, ONLY: CD_Current_Profile_Velocity, CD_SeaState_Steady_Current, CD_HYDRO_OK
   USE CableDyn_DeckDriver, ONLY: CD_Deck_Query_dtM, CD_DECKDRV_OK
   USE CableDyn_OpenFAST_Aggregate, ONLY: CD_AGG_ModuleType, CD_AGG_Init_From_Deck, CD_AGG_UpdateStates_Moving, &
                                          CD_AGG_Step_Moving, CD_AGG_CalcOutput, CD_AGG_GetMovingPointMesh, &
                                          CD_AGG_Snapshot, CD_AGG_Restore, &
                                          CD_AGG_NMovingPoints, CD_AGG_NDynamicPoints, CD_AGG_Refresh_PointMesh, &
                                          CD_AGG_NFailures, &
                                          CD_AGG_HasRigid6, CD_AGG_HasRod, &
                                          CD_AGG_Rigid6_MirrorSize, CD_AGG_Get_Rigid6_States, &
                                          CD_AGG_Set_Rigid6_States, &
                                          CD_AGG_Rod_MirrorSize, CD_AGG_Get_Rod_States, CD_AGG_Set_Rod_States, &
                                          CD_AGG_NCtrlChans, CD_AGG_Apply_LineControl, &
                                          CD_AGG_Get_Failure_Flags, CD_AGG_Set_Failure_Flags, CD_AGG_End, &
                                          CD_AGG_IsInitialized, CD_AGG_OK, &
                                          CD_AGG_NTurbines, CD_AGG_TurbineOfMoving, &
                                          CD_AGG_NFluidNodes, CD_AGG_GetFluidNodePositions, CD_AGG_SetFluidFields, &
                                          CD_AGG_NumChannels, CD_AGG_ChannelHeader, CD_AGG_EvalChannel, &
                                          CD_AGG_GetInitMetadata, CD_AGG_GetInitLine, &
                                          CD_AGG_WriteStaticProfile
   USE CableDyn_OpenFAST_HermiteFMF, ONLY: CD_HFMF_MirrorSize, CD_HFMF_PackMirror, CD_HFMF_UnpackMirror, &
                                           CD_HFMF_FrictionMirrorSize, CD_HFMF_ForceMirrorSize, &
                                           CD_HFMF_OK
   USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCable_Dyn_Reset_Profile, &
                                            CD_HermiteCable_Dyn_Disable_Profile, &
                                            CD_HermiteCable_Dyn_Get_Profile, &
                                            CD_HermiteCable_Dyn_Recovery_Count, &
                                            CD_HermiteCable_Dyn_Recovery_Reset
   USE CableDyn_System, ONLY: CD_System_NLines, CD_System_Line_NDOF, CD_System_Line_NElem, &
                              CD_Get_System_DynamicPoint_States, CD_Set_System_DynamicPoint_States, &
                              CD_System_NSystemCoupledDOF, CD_Get_System_CoupledMotion, &
                              CD_Update_System_CoupledMotion, &
                              CD_Get_System_Line_State, CD_Update_System_Line_Interior_State, &
                              CD_System_Line_Has_Viscoelastic, CD_Get_System_Line_VE_Dl1, &
                              CD_Set_System_Line_VE_Dl1, &
                              CD_System_Line_Has_Syrope, CD_Get_System_Line_Syrope_State, &
                              CD_Set_System_Line_Syrope_State, &
                              CD_System_Line_Has_Friction, CD_Get_System_Line_Friction_Anchors, &
                              CD_Set_System_Line_Friction_Anchors, &
                              CD_Recompute_System_Acceleration, CD_SYSTEM_OK
   USE SeaSt_WaveField, ONLY: WaveField_GetNodeWaveKin
   USE SeaSt_WaveField_Types, ONLY: WaveMod_None

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER :: CableDyn_ProgDesc = ProgDesc( 'CableDyn', 'v0.1.0', '2026-10-01' )

   !> Module-level instance registry: the CableDyn solver objects cannot live inside the
   !> registry-generated types, so each initialized module holds an integer handle
   !> (m%CDInst) into this PRIVATE table. A handle that does not refer to an ACTIVE
   !> entry (e.g. after a checkpoint restore in a fresh process) fails closed.
   TYPE :: CD_OF_Instance
      TYPE(CD_AGG_ModuleType) :: agg
      INTEGER(IntKi)          :: ncp     = 0        ! number of coupled points (mesh nodes)
      INTEGER(IntKi)          :: nlines  = 0
      INTEGER(IntKi)          :: last_n  = -1       ! last completed glue step number (correction detection)
      INTEGER(IntKi)          :: last_solve_n = -1  ! glue-step index of the last actual mooring solve
      ! CableDyn's OWN time step (supercycling): the implicit mooring solve advances at
      ! dt_moor = nsuper * (glue DT), ONE step per nsuper glue steps, instead of being
      ! chained to the turbine's small DT. The deck's OPTIONS dtM sets it (clamped to at
      ! least the glue DT); absent a deck dtM the default target is 0.1 s. Between mooring
      ! solves CD_CalcOutput returns the last
      ! committed load state as a zero-order hold; feeding instantaneous endpoints into a
      ! lagged mooring interior creates artificial stiff-EA load chatter.
      INTEGER(IntKi)          :: nsuper  = 1        ! glue steps per mooring step
      REAL(wp)                :: dt_moor = 0.0_wp   ! the mooring step frozen into the aggregate
      LOGICAL                 :: active  = .FALSE.
      ! Persistent marching-path workspace (3, ncp): CD_UpdateStates and CD_CalcOutput run
      ! every glue step of a long simulation, so their kinematics/load scratch must not be
      ! heap-allocated per call. Grow-only: allocated on first use, reused thereafter.
      REAL(wp), ALLOCATABLE   :: ws_pos(:,:), ws_vel(:,:), ws_acc(:,:), ws_fld(:,:), ws_moment(:,:)
      REAL(wp), ALLOCATABLE   :: ws_omega(:,:), ws_alpha(:,:)
      REAL(wp), ALLOCATABLE   :: ws_orient(:,:,:)
      ! saved committed coupled kinematics for the state-const CalcOutput at nsuper = 1
      ! (save -> transfer u -> evaluate -> restore; see CD_CalcOutput)
      REAL(wp), ALLOCATABLE   :: sv_pos(:,:), sv_vel(:,:), sv_acc(:,:), sv_orient(:,:,:)
      REAL(wp), ALLOCATABLE   :: sv_omega(:,:), sv_alpha(:,:)
      ! Supercycled-CalcOutput cache: between mooring steps the committed state does not
      ! change, so held loads/channels are recomputed only after a new solve commit.
      LOGICAL                 :: loads_cached = .FALSE.
      INTEGER(IntKi)          :: loads_at_n = -1
      REAL(wp), ALLOCATABLE   :: ws_chan(:)
      ! Persistent interpolated-input container: CD_UpdateStates interpolates the coupled
      ! kinematics to the mooring interval end; the mesh copy is built ONCE and reused
      ! (ExtrapInterp overwrites the fields in place), never created/destroyed per step.
      TYPE(CD_InputType)      :: u_interp
      LOGICAL                 :: u_interp_ready = .FALSE.
      ! State-mirror packing scratch (pack_state_mirror runs once per actual mooring step;
      ! as locals its grow-only guards were dead logic -- automatic deallocation on return
      ! made every call re-allocate). Grow-only across calls here.
      REAL(wp), ALLOCATABLE   :: ws_qL(:), ws_vL(:), ws_aL(:), ws_cbuf(:)
      ! committed-state mirror captured before a state-const CalcOutput transfer;
      ! restored after (the boundary-only restore recomputes interior accelerations,
      ! which are gen-alpha COMMITTED values, not (q, v)-derivable)
      REAL(R8Ki), ALLOCATABLE :: ws_msave(:)
      ! FAST.Farm: turbine count (0 = plain single-turbine mode), the driving turbine
      ! per flat aggregate moving slot, that slot's node index within its turbine's
      ! mesh, the farm-layout reference offsets, and the per-turbine PtfmInit.
      INTEGER(IntKi)          :: nT = 0
      INTEGER(IntKi), ALLOCATABLE :: tmap(:), node_of_slot(:)
      REAL(R8Ki), ALLOCATABLE :: farm_ref(:,:), farm_ptfm(:,:)
      ! CableDyn-owned committed-step record. Unlike OpenFAST WriteOutput (sampled on the
      ! host DT_Out clock and therefore repeated under supercycling), this writes exactly
      ! one row per committed dtM solve. A correction re-evaluation overwrites the last row.
      ! FAST.Farm keeps its established <RootName>.out name; a single-turbine run writes
      ! <OpenFASTRoot>.CD.out. native_un < 0 means no deck time-history channels.
      INTEGER(IntKi)          :: native_un = -1
      INTEGER(IntKi)          :: native_wrote_n = -123456789
      REAL(DbKi)              :: native_wrote_t = 0.0_DbKi
      LOGICAL                 :: native_pending = .FALSE.
      CHARACTER(1024)         :: native_root = ''
      ! SeaState ambient-fluid sampling (host-driven wave/current kinematics): when the
      ! glue's WaveField declares waves or a current, the aggregate is built wave-capable
      ! (external_fluid) and every mooring advance samples the field at the line nodes'
      ! committed positions -- velocity, acceleration, and the local surface elevation --
      ! and prescribes them for the step (a per-mooring-step zero-order hold, the same
      ! weak-coupling cadence the loads-out side documents). Workspaces sized at Init.
      LOGICAL                 :: use_extfluid = .FALSE.
      INTEGER(IntKi)          :: nfluid = 0
      REAL(wp), ALLOCATABLE   :: fx_xyz(:,:), fx_vel(:,:), fx_acc(:,:), fx_wl(:), fx_pd(:)
      ! CHECKPOINT-RESTART identity: serial is the generation stamp minted at
      ! Init/rebuild and mirrored into other%restart_serial; t_commit is the time of the
      ! last committed mooring advance, mirrored into other%t_commit. A restored
      ! OtherState that does not match the live instance (stamp mismatch = fresh
      ! process; live clock AHEAD of the restored clock = same-process rewind) triggers
      ! rebuild_from_checkpoint.
      INTEGER(IntKi)          :: serial = 0
      REAL(DbKi)              :: t_commit = -1.0_DbKi
      ! Opt-in end-to-end profiler. CABLEDYN_PROFILE=1 arms these counters and
      ! the existing detailed Hermite profiler; ordinary production runs pay no
      ! SYSTEM_CLOCK overhead. Times are accumulated only for successful calls.
      LOGICAL                 :: prof_enabled = .FALSE.
      INTEGER(IntKi)          :: prof_n_update = 0, prof_n_solve = 0
      INTEGER(IntKi)          :: prof_n_calcout = 0, prof_n_fresh = 0, prof_n_wave = 0
      REAL(wp)                :: prof_t_init = 0.0_wp, prof_t_update = 0.0_wp
      REAL(wp)                :: prof_t_wave = 0.0_wp, prof_t_step = 0.0_wp
      REAL(wp)                :: prof_t_calcout = 0.0_wp, prof_t_loads = 0.0_wp
      REAL(wp)                :: prof_t_native = 0.0_wp
   END TYPE CD_OF_Instance

   TYPE(CD_OF_Instance), ALLOCATABLE, SAVE :: Inst(:)
   INTEGER(IntKi), SAVE :: g_serial = 0   ! monotonic generation counter (Init/rebuild mints)

   PUBLIC :: CD_Init
   PUBLIC :: CD_UpdateStates
   PUBLIC :: CD_CalcOutput
   PUBLIC :: CD_CalcContStateDeriv
   PUBLIC :: CD_End
   PUBLIC :: CD_JacobianPContState
   PUBLIC :: CD_JacobianPInput
   PUBLIC :: CD_JacobianPDiscState
   PUBLIC :: CD_JacobianPConstrState

CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE open_native_outputs( id, root, t_start, ErrStat, ErrMsg, write_profile )
      !! Create a self-contained CableDyn result set at the current committed state.
      !! A checkpoint rebuild must first release or adopt any process-local owner of
      !! the same path through restore_native_outputs. write_profile = .FALSE. keeps a
      !! static profile the caller has already written (CD_Init writes it at the
      !! equilibrium, before the host kinematics are seeded).
      INTEGER(IntKi),       INTENT(IN)    :: id
      CHARACTER(*),        INTENT(IN)    :: root
      REAL(DbKi),           INTENT(IN)    :: t_start
      INTEGER(IntKi),       INTENT(OUT)   :: ErrStat
      CHARACTER(*),         INTENT(OUT)   :: ErrMsg
      LOGICAL, OPTIONAL,    INTENT(IN)    :: write_profile

      character(*), parameter :: RoutineName = 'CableDyn:open_native_outputs'
      INTEGER(IntKi)          :: es, un, ic, nch, owner, ios
      CHARACTER(512)          :: em, iomsg
      CHARACTER(ChanLen)      :: hdr, unt
      REAL(wp)                :: value
      LOGICAL                 :: write_static

      ErrStat = ErrID_None
      ErrMsg = ''
      IF (id < 1_IntKi .OR. .NOT. ALLOCATED(Inst) .OR. id > SIZE(Inst)) THEN
         CALL SetErrStat( ErrID_Fatal, 'Invalid solver-instance handle for native output.', &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (Inst(id)%native_un > 0) CLOSE (Inst(id)%native_un)
      Inst(id)%native_un = -1_IntKi
      Inst(id)%native_wrote_n = -123456789_IntKi
      Inst(id)%native_wrote_t = MAX(t_start, 0.0_DbKi)
      Inst(id)%native_pending = .FALSE.
      Inst(id)%native_root = TRIM(root)

      write_static = .TRUE.
      IF (PRESENT(write_profile)) write_static = write_profile
      IF (write_static) THEN
         CALL CD_AGG_WriteStaticProfile( Inst(id)%agg, TRIM(root)//'.static.out', es, em )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      END IF

      nch = CD_AGG_NumChannels(Inst(id)%agg)
      IF (nch <= 0) RETURN
      ! One run-level path has one owner. A rebuilt state must adopt the existing
      ! unit through restore_native_outputs rather than attempting a second OPEN.
      DO owner = 1, SIZE(Inst)
         IF (owner == id) CYCLE
         IF (Inst(owner)%active .AND. Inst(owner)%native_un > 0 .AND. &
             TRIM(Inst(owner)%native_root) == TRIM(root)) THEN
            CALL SetErrStat( ErrID_Fatal, 'Native output path "'//TRIM(root)//'.out" is already owned '// &
                             'by another live CableDyn instance.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      END DO
      CALL GetNewUnit(un)
      ! Read/write access is deliberate. The glue can restore a registered state
      ! after a row has been committed, in which case restore_native_outputs must
      ! BACKSPACE and replace the abandoned row. OpenFOutFile opens write-only and
      ! therefore makes that standard record-positioning operation fail on gfortran.
      iomsg = ''
      OPEN (un, FILE=TRIM(root)//'.out', STATUS='REPLACE', FORM='FORMATTED', &
            ACTION='READWRITE', POSITION='REWIND', IOSTAT=ios, IOMSG=iomsg)
      IF (ios /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Cannot open native output file "'//TRIM(root)// &
                          '.out" for transactional read/write access. '//TRIM(iomsg), &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      WRITE (un, '(A)', ADVANCE='NO') 'Time'
      DO ic = 1, nch
         CALL CD_AGG_ChannelHeader(Inst(id)%agg, ic, hdr, unt, es, em)
         IF (es /= CD_AGG_OK) THEN
            CLOSE (un, STATUS='DELETE')
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         WRITE (un, '(A)', ADVANCE='NO') CHAR(9)//TRIM(hdr)
      END DO
      WRITE (un, '(A)') ''
      WRITE (un, '(A)', ADVANCE='NO') '(s)'
      DO ic = 1, nch
         CALL CD_AGG_ChannelHeader(Inst(id)%agg, ic, hdr, unt, es, em)
         IF (es /= CD_AGG_OK) THEN
            CLOSE (un, STATUS='DELETE')
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         WRITE (un, '(A)', ADVANCE='NO') CHAR(9)//TRIM(unt)
      END DO
      WRITE (un, '(A)') ''
      WRITE (un, '(ES25.16E3)', ADVANCE='NO') MAX(t_start, 0.0_DbKi)
      DO ic = 1, nch
         CALL CD_AGG_EvalChannel(Inst(id)%agg, ic, value, es, em)
         IF (es /= CD_AGG_OK) THEN
            CLOSE (un, STATUS='DELETE')
            CALL SetErrStat( ErrID_Fatal, 'Initial committed-step channel record: '//TRIM(em), &
                             ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         WRITE (un, '(A,ES15.6E2)', ADVANCE='NO') CHAR(9), value
      END DO
      WRITE (un, '(A)') ''
      FLUSH (un)
      Inst(id)%native_un = un
      ! At initialization no glue step has advanced. A rebuilt output begins at the
      ! checkpoint state and the next solve receives a non-negative glue-step index.
      Inst(id)%native_wrote_n = -1_IntKi
      Inst(id)%native_wrote_t = MAX(t_start, 0.0_DbKi)
   END SUBROUTINE open_native_outputs

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE restore_native_outputs( id, old_id, root, t_start, restart_serial, dt_moor, nsuper, &
                                      ErrStat, ErrMsg )
      !! Transfer the run-level native history from a superseded in-process solver
      !! instance to its checkpoint reconstruction. If the restored state rewinds
      !! time, remove the abandoned committed rows before handing over the still-open
      !! unit. A fresh-process restart has no live owner and creates a new history.
      INTEGER(IntKi), INTENT(IN)  :: id, old_id, restart_serial, nsuper
      CHARACTER(*),   INTENT(IN)  :: root
      REAL(DbKi),      INTENT(IN)  :: t_start, dt_moor
      INTEGER(IntKi), INTENT(OUT) :: ErrStat
      CHARACTER(*),   INTENT(OUT) :: ErrMsg

      character(*), parameter :: RoutineName = 'CableDyn:restore_native_outputs'
      INTEGER(IntKi)          :: ntrim, k, ios
      REAL(DbKi)              :: delta, tol
      LOGICAL                 :: live_owner, same_history

      ErrStat = ErrID_None
      ErrMsg = ''
      live_owner = .FALSE.
      same_history = .FALSE.
      IF (ALLOCATED(Inst)) THEN
         IF (old_id >= 1_IntKi .AND. old_id <= SIZE(Inst)) THEN
            live_owner = Inst(old_id)%active .AND. Inst(old_id)%native_un > 0 .AND. &
                         TRIM(Inst(old_id)%native_root) == TRIM(root)
            same_history = live_owner .AND. Inst(old_id)%serial == restart_serial
         END IF
      END IF
      IF (.NOT. live_owner) THEN
         CALL open_native_outputs( id, root, t_start, ErrStat, ErrMsg )
         RETURN
      END IF

      ! A fresh-process restart first passes through CD_Init, which creates a
      ! provisional t=0 stream. The restored generation stamp differs, so release
      ! that stream and recreate it at the checkpoint. It must not be treated as a
      ! same-process rewind because it contains no pre-checkpoint history.
      IF (.NOT. same_history) THEN
         CLOSE (Inst(old_id)%native_un, IOSTAT=ios)
         IF (ios /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not release the provisional native-output stream during restart.', &
                             ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         Inst(old_id)%native_un = -1_IntKi
         Inst(old_id)%native_root = ''
         Inst(old_id)%native_wrote_n = -123456789_IntKi
         Inst(old_id)%native_wrote_t = 0.0_DbKi
         Inst(old_id)%native_pending = .FALSE.
         CALL open_native_outputs( id, root, t_start, ErrStat, ErrMsg )
         RETURN
      END IF

      tol = 1.0E-7_DbKi*MAX(1.0_DbKi, ABS(t_start), ABS(Inst(old_id)%native_wrote_t))
      delta = Inst(old_id)%native_wrote_t - t_start
      IF (delta < -tol) THEN
         CALL SetErrStat( ErrID_Fatal, 'The live native history ends before the restored committed state.', &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      ntrim = MAX(0_IntKi, NINT(MAX(delta, 0.0_DbKi)/dt_moor, IntKi))
      IF (ABS(delta - REAL(ntrim, DbKi)*dt_moor) > tol) THEN
         CALL SetErrStat( ErrID_Fatal, 'The native-history rewind is not an integer number of CableDyn steps.', &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      DO k = 1, ntrim
         BACKSPACE (Inst(old_id)%native_un, IOSTAT=ios)
         IF (ios /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not remove abandoned native-output record '// &
                             TRIM(Num2LStr(k))//' of '//TRIM(Num2LStr(ntrim))//' while rewinding from t = '// &
                             TRIM(Num2LStr(Inst(old_id)%native_wrote_t))//' s to t = '// &
                             TRIM(Num2LStr(t_start))//' s.', &
                             ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      END DO
      IF (ntrim > 0_IntKi) THEN
         ENDFILE (Inst(old_id)%native_un, IOSTAT=ios)
         IF (ios == 0) BACKSPACE (Inst(old_id)%native_un, IOSTAT=ios)
         IF (ios /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not truncate the abandoned native-output branch.', &
                             ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      END IF

      Inst(id)%native_un = Inst(old_id)%native_un
      Inst(id)%native_root = Inst(old_id)%native_root
      Inst(id)%native_wrote_t = t_start
      IF (Inst(old_id)%native_wrote_n >= 0_IntKi) THEN
         Inst(id)%native_wrote_n = MAX(-1_IntKi, Inst(old_id)%native_wrote_n - ntrim*MAX(1_IntKi, nsuper))
      ELSE
         Inst(id)%native_wrote_n = Inst(old_id)%native_wrote_n
      END IF
      Inst(id)%native_pending = .FALSE.

      Inst(old_id)%native_un = -1_IntKi
      Inst(old_id)%native_root = ''
      Inst(old_id)%native_wrote_n = -123456789_IntKi
      Inst(old_id)%native_wrote_t = 0.0_DbKi
      Inst(old_id)%native_pending = .FALSE.
   END SUBROUTINE restore_native_outputs

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CD_Init(InitInp, u, p, x, xd, z, other, y, m, DTcoupling, InitOut, ErrStat, ErrMsg)
      TYPE(CD_InitInputType),       INTENT(IN   )  :: InitInp
      TYPE(CD_InputType),           INTENT(  OUT)  :: u
      TYPE(CD_ParameterType),       INTENT(  OUT)  :: p
      TYPE(CD_ContinuousStateType), INTENT(  OUT)  :: x
      TYPE(CD_DiscreteStateType),   INTENT(  OUT)  :: xd
      TYPE(CD_ConstraintStateType), INTENT(  OUT)  :: z
      TYPE(CD_OtherStateType),      INTENT(  OUT)  :: other
      TYPE(CD_OutputType),          INTENT(  OUT)  :: y
      TYPE(CD_MiscVarType),         INTENT(  OUT)  :: m
      REAL(DbKi),                   INTENT(INOUT)  :: DTcoupling
      TYPE(CD_InitOutputType),      INTENT(  OUT)  :: InitOut
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg

      character(*), parameter   :: RoutineName = 'CD_Init'
      INTEGER(IntKi)            :: ErrStat2, es
      CHARACTER(ErrMsgLen)      :: ErrMsg2
      CHARACTER(512)            :: em
      INTEGER(IntKi)            :: id, ncp, nlines, il, i, ndofl, nstates, ic, nch, nsuper, nglue
      INTEGER(IntKi)            :: nreportlines, npoints, nsections, nei0, nfinite, line_id
      INTEGER(IntKi)            :: clock_count, nT, prof_c0, prof_c1, prof_rate
      LOGICAL                   :: farm_has_ptfm, profile_requested
      REAL(wp), ALLOCATABLE     :: pos(:,:), vel(:,:), acc(:,:), fld(:,:)
      REAL(R8Ki)                :: OrientIdent(3,3)
      REAL(R8Ki)                :: OrMatInit(3,3), ptfm6(6), refpos(3)
      REAL(wp)                  :: cval, dtC, deck_dtm, dt_moor
      REAL(wp)                  :: fair_ten, fair_force(3), fair_incl, fair_decl, fair_azi, force_incl
      LOGICAL                   :: has_dtm, has_ptfm, use_extf
      CHARACTER(ChanLen)        :: chan_hdr, chan_unt
      CHARACTER(1024)           :: native_root
      REAL(wp), PARAMETER       :: DTM_DEFAULT = 0.1_wp

      ErrStat = ErrID_None
      ErrMsg  = ''
      profile_requested = profile_environment_requested()
      IF (profile_requested) THEN
         CALL SYSTEM_CLOCK(prof_c0, prof_rate)
      ELSE
         ! Disarm before any static/dynamic initialization work. This also covers an
         ! initialization that later fails and never reaches the live-instance setup.
         CALL CD_HermiteCable_Dyn_Disable_Profile()
      END IF

      ! ------------------------------------------------------------------------------
      ! The CableDyn identity, announced at OpenFAST initialization the way every
      ! OpenFAST module does (DispNVD one-line version banner, then WrScr identity
      ! lines; cf. MoorDyn.f90 Init). DispCopyrightLicense is deliberately NOT used --
      ! it hard-codes the NREL/Envision copyright, whereas CableDyn is (C) SMI Lab,
      ! Inha University under Apache-2.0. InitOut%Ver feeds the glue summary file.
      ! ------------------------------------------------------------------------------
      CALL DispNVD( CableDyn_ProgDesc )
      InitOut%Ver = CableDyn_ProgDesc
      CALL WrScr( '   CableDyn: geometrically nonlinear cable & mooring dynamics for floating wind.' )
      CALL WrScr( '   Author Prof. Jae Hoon Seo, Inha University.  License Apache-2.0.' )
      CALL WrScr( '   CompMooring = 5: mooring and dynamic power-cable loads are computed by CableDyn,' )
      CALL WrScr( '   an OpenFAST mooring module alongside MoorDyn (CompMooring = 3) in this binary.' )

      ! ------------------------------------------------------------------------------
      ! Scope guards (named fail-closed messages; see the module header)
      ! ------------------------------------------------------------------------------
      IF (InitInp%FarmSize > 0 .AND. InitInp%Linearize) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: linearization in FAST.Farm mode is '// &
                          'not supported (linearize single-turbine cases).', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      ! Linearize is accepted for the established force-only boundary: the module presents
      ! the MAP-shaped reduction -- ZERO
      ! continuous-state variables (CD_InitVars) and a QUASI-STATIC dYdu whose
      ! displacement columns come from scratch-aggregate static re-solves at the
      ! perturbed boundary (each init IS a static solve), so the linearized mooring
      ! is the re-equilibrated stiffness, not the frozen-interior response
      ! (CD_JacobianPInput). Velocity/acceleration/orientation columns are zero:
      ! frequency-dependent mooring damping is not representable without states. A
      ! platform-relative cable end connection is rejected after deck construction below,
      ! because it introduces orientation-dependent force and moment columns.
      IF (.NOT. InitInp%UsePrimaryInputFile) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: passed-in primary input data is not '// &
                          'supported; provide the MoorDyn-format deck as a file.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      ! Initial platform displacement (PtfmInit): the deck's Coupled/Vessel points are
      ! rigid-transformed BEFORE the build, so the true-Newton static equilibrium is
      ! solved at the DISPLACED pose -- correct by construction, never an implicit step
      ! from the undisplaced statics. Mirrors MoorDyn-F's convention exactly: the ZYX
      ! (3-2-1) DCM applied transposed, mesh reference positions kept at the UNDISPLACED
      ! deck coordinates, and the initial displacement carried in TranslationDisp (the
      ! glue's motion mapping anchors on reference positions).
      ptfm6 = 0.0_R8Ki
      nT = MAX(0_IntKi, InitInp%FarmSize)
      IF (nT > 0) THEN
         ! FAST.Farm: per-turbine PtfmInit + farm-layout reference offsets are consumed
         ! below through the aggregate's farm mode; the single-turbine ptfm path is off.
         IF (.NOT. ALLOCATED(InitInp%TurbineRefPos)) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: FAST.Farm mode needs '// &
                             'TurbineRefPos(3, FarmSize).', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         IF (SIZE(InitInp%TurbineRefPos, 1) /= 3 .OR. SIZE(InitInp%TurbineRefPos, 2) /= nT) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: TurbineRefPos must be shaped '// &
                             '(3, FarmSize).', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      ELSE IF (ALLOCATED(InitInp%PtfmInit)) THEN
         ! SIZE > 0 first: MAXVAL of a zero-size array is -HUGE (well-defined but obscure);
         ! an empty allocation reads explicitly as "no displacement declared".
         IF (SIZE(InitInp%PtfmInit, 1) >= 6 .AND. SIZE(InitInp%PtfmInit, 2) >= 1) THEN
            ptfm6 = REAL(InitInp%PtfmInit(1:6, 1), R8Ki)
         END IF
      END IF
      has_ptfm = MAXVAL(ABS(ptfm6)) > 0.0_R8Ki
      IF (has_ptfm) THEN
         WRITE (em, '(A,6(1X,ES24.16))') '   CableDyn initial platform pose [x y z roll pitch yaw]:', ptfm6
         CALL WrScr(TRIM(em))
      END IF

      ! Keep a packable copy of the deck in the parameters so a checkpoint carries it.
      CALL ProcessComFile( InitInp%FileName, p%DeckCopy, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! ------------------------------------------------------------------------------
      ! Build the CableDyn system from the deck: parse, feature-gate (the deck layer
      ! carries its own named rejections), TRUE Newton static equilibrium, models.
      ! ------------------------------------------------------------------------------
      id = new_instance()

      ! ------------------------------------------------------------------------------
      ! CableDyn's OWN time step. The implicit generalised-alpha solve can remain stable at a
      ! large mooring step, but stability does not establish time-domain or fatigue accuracy.
      ! The practical default targets 0.1 s, advancing at
      ! dt_moor = nsuper * DT (one implicit step per nsuper glue steps): the deck's OPTIONS
      ! dtM is the author's intent (rounded to the nearest integer multiple of the glue DT
      ! and clamped to at least one glue step -- a MoorDyn-migrated deck's explicit-RK4
      ! dtM like 0.001 s does not transfer to an implicit solver); absent a deck dtM, the
      ! target is DTM_DEFAULT. The CableDyn-owned .CD.out records only committed dtM rows;
      ! OpenFAST's assembled output still receives zero-order-held module channels at DT_Out.
      ! Between coarse mooring steps, CD_CalcOutput returns the last committed mooring load
      ! state as a zero-order hold; injecting instantaneous
      ! endpoint motion into a lagged line interior creates nonphysical stiff-EA chatter.
      ! ------------------------------------------------------------------------------
      dtC = REAL(DTcoupling, wp)
      CALL CD_Deck_Query_dtM( TRIM(InitInp%FileName), deck_dtm, has_dtm, es, em, &
                              env_gravity=REAL(InitInp%g, wp), env_rho_water=REAL(InitInp%rhoW, wp), &
                              env_wtrdpth=REAL(InitInp%WtrDepth, wp) )
      IF (es /= CD_DECKDRV_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (has_dtm) THEN
         ! Like MoorDyn-F (which shortens dtM so a whole number of its steps fills the
         ! coupling interval), an inexact deck dtM is only ever refined, never coarsened:
         ! use the largest multiple of the glue DT that does not exceed it, and say so as
         ! a WARNING stating both values. A dtM below DT is clamped to one glue step (the
         ! documented implicit-solver floor for MoorDyn-migrated explicit dtM values).
         nsuper = MAX(1_IntKi, FLOOR(deck_dtm/dtC*(1.0_wp + 1.0e-6_wp), IntKi))
         IF (deck_dtm < dtC*(1.0_wp - 1.0e-6_wp)) THEN
            CALL WrScr( '   CableDyn: deck dtM = '//TRIM(Num2LStr(deck_dtm))//' s is below the glue DT = '// &
                        TRIM(Num2LStr(dtC))//' s; the implicit mooring step uses dtM = DT.' )
         ELSE IF (ABS(REAL(nsuper, wp)*dtC - deck_dtm) > 1.0e-6_wp*MAX(dtC, deck_dtm)) THEN
            CALL SetErrStat( ErrID_Warn, 'CableDyn mooring module: deck dtM = '//TRIM(Num2LStr(deck_dtm))// &
                             ' s is not an integer multiple of the glue DT = '//TRIM(Num2LStr(dtC))// &
                             ' s; using dtM = '//TRIM(Num2LStr(REAL(nsuper, wp)*dtC))//' s ('// &
                             TRIM(Num2LStr(nsuper))//' x DT), the largest multiple not exceeding the '// &
                             'requested step. Set dtM to a multiple of DT to silence this warning.', &
                             ErrStat, ErrMsg, RoutineName )
            CALL WrScr( '   WARNING: CableDyn deck dtM = '//TRIM(Num2LStr(deck_dtm))//' s changed to '// &
                        TRIM(Num2LStr(REAL(nsuper, wp)*dtC))//' s ('//TRIM(Num2LStr(nsuper))// &
                        ' x glue DT = '//TRIM(Num2LStr(dtC))//' s).' )
         END IF
      ELSE
         nsuper = MAX(1_IntKi, NINT(DTM_DEFAULT/dtC, IntKi))
      END IF
      dt_moor = REAL(nsuper, wp)*dtC
      CALL WrScr( '   CableDyn time step dtM = '//TRIM(Num2LStr(dt_moor))//' s ('//TRIM(Num2LStr(nsuper))// &
                  ' x glue DT'//TRIM(MERGE('; deck dtM   ', '; the default', has_dtm))//')' )
      IF (nsuper > 1_IntKi) THEN
         CALL WrScr( '   CableDyn: coupled loads use a zero-order hold between dtM solves; the separate .CD.out '// &
                     'records committed CableDyn samples only. Verify dtM convergence for fatigue, snap, contact, '// &
                     'or local cable-response calculations; set deck dtM = DT for glue-rate solves.' )
      END IF
      ! The run must END on a mooring boundary: with supercycling, a trailing partial
      ! interval would take only skipped glue steps, leaving the committed mooring state --
      ! and the held (zero-order-hold) loads -- up to dtM stale at TMax, silently. Fail
      ! closed at Init with the fix named rather than finish on a stale state.
      IF (nsuper > 1_IntKi) THEN
         nglue = NINT(REAL(InitInp%Tmax, wp)/dtC, IntKi)
         IF (ABS(REAL(nglue, wp)*dtC - REAL(InitInp%Tmax, wp)) > 1.0e-6_wp*MAX(dtC, REAL(InitInp%Tmax, wp)) .OR. &
             MOD(nglue, nsuper) /= 0_IntKi) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: TMax = '//TRIM(Num2LStr(InitInp%Tmax))// &
                             ' s is not a whole number of CableDyn steps (dtM = '//TRIM(Num2LStr(dt_moor))// &
                             ' s = '//TRIM(Num2LStr(nsuper))//' x DT): the run would end inside a mooring '// &
                             'interval with the committed state up to dtM stale. Choose TMax as a multiple '// &
                             'of dtM, or set the deck OPTIONS dtM accordingly.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      END IF

      ! SeaState ambient fluid: when the glue hands a WaveField carrying waves or a
      ! current, the mooring receives it -- the aggregate is built WAVE-CAPABLE
      ! (external_fluid: the full drag/FK/added-mass/wetting hydro configuration at zero
      ! initial fields) and every mooring advance samples the field at the line nodes.
      ! A still-water field (WaveMod 0, no current) keeps the plain build, bit-for-bit.
      !
      ! FARM SCOPE: farm-level mooring runs in still water. The only WaveField reachable
      ! at the farm level is turbine 1's TURBINE-LOCAL field, whose spatial grid (SeaState
      ! X/Y_HalfWidth, typically ~100 m) does not cover a multi-kilometre farm mooring
      ! span -- sampling it beyond the grid returns clamped, nonphysical kinematics.
      ! Stock farm-level MoorDyn behaves the same way in practice (it holds the field
      ! pointer but samples nothing unless its deck declares WaterKin). Farm-wide wave
      ! kinematics on shared moorings needs a farm-domain field and is deferred with the
      ! scope stated here rather than silently wrong at the grid edge.
      use_extf = .FALSE.
      IF (ASSOCIATED(InitInp%WaveField) .AND. nT == 0_IntKi) THEN
         ! hasCurrField covers a dynamically-populated current field (e.g. an
         ! InflowWind-driven MHK current) that carries no CurrMod declaration --
         ! without it such a run would silently leave the mooring in still water.
         use_extf = wavefield_is_ambient( InitInp%WaveField )
         IF (use_extf) THEN
            ! the mooring water column must be the SeaState one (MoorDyn's own guard set)
            IF (ABS(InitInp%WaveField%WtrDpth - InitInp%WtrDepth) > 1.0E-3_ReKi) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: SeaState water depth ('// &
                                TRIM(Num2LStr(InitInp%WaveField%WtrDpth))//' m) does not match the mooring '// &
                                'water depth ('//TRIM(Num2LStr(InitInp%WtrDepth))//' m).', ErrStat, ErrMsg, &
                                RoutineName )
               RETURN
            END IF
         END IF
      END IF

      ! The host owns the environment: gravity, water density, and depth override
      ! whatever the deck declares (MoorDyn-faithful coupled behavior).
      CALL WrScr( '  Parsing CableDyn input file: '//TRIM(InitInp%FileName) )
      IF (nT > 0) THEN
         ! an allocated-but-empty/too-small PtfmInit means "no initial displacement"
         ! (the single-turbine path above treats that case as absent); the aggregate
         ! requires a PRESENT farm array to be exactly (6, n_turbines), so gate the
         ! pass-through on shape and omit the optional otherwise
         farm_has_ptfm = .FALSE.
         IF (ALLOCATED(InitInp%PtfmInit)) THEN
            farm_has_ptfm = SIZE(InitInp%PtfmInit, 1) >= 6 .AND. SIZE(InitInp%PtfmInit, 2) >= nT
         END IF
         ! forbid_deck_ambient: farm scope is STILL WATER (see the use_extf guard above),
         ! and that must hold against the deck too -- a deck-declared waves/current OPTION
         ! would silently self-drive ambient forcing the documented scope excludes.
         IF (farm_has_ptfm) THEN
            CALL CD_AGG_Init_From_Deck( Inst(id)%agg, TRIM(InitInp%FileName), dt_moor, es, em, &
                                        env_gravity=REAL(InitInp%g, wp), env_rho_water=REAL(InitInp%rhoW, wp), &
                                        env_wtrdpth=REAL(InitInp%WtrDepth, wp), external_fluid=use_extf, &
                                        n_turbines=nT, &
                                        turbine_ref_pos=REAL(InitInp%TurbineRefPos, wp), &
                                        farm_ptfm_init=REAL(InitInp%PtfmInit(1:6, 1:nT), wp), &
                                        forbid_deck_ambient=.TRUE., run_tmax=REAL(InitInp%Tmax, wp) )
         ELSE
            CALL CD_AGG_Init_From_Deck( Inst(id)%agg, TRIM(InitInp%FileName), dt_moor, es, em, &
                                        env_gravity=REAL(InitInp%g, wp), env_rho_water=REAL(InitInp%rhoW, wp), &
                                        env_wtrdpth=REAL(InitInp%WtrDepth, wp), external_fluid=use_extf, &
                                        n_turbines=nT, &
                                        turbine_ref_pos=REAL(InitInp%TurbineRefPos, wp), &
                                        forbid_deck_ambient=.TRUE., run_tmax=REAL(InitInp%Tmax, wp) )
         END IF
      ELSE
         CALL CD_AGG_Init_From_Deck( Inst(id)%agg, TRIM(InitInp%FileName), dt_moor, es, em, &
                                     env_gravity=REAL(InitInp%g, wp), env_rho_water=REAL(InitInp%rhoW, wp), &
                                     env_wtrdpth=REAL(InitInp%WtrDepth, wp), ptfm_init=REAL(ptfm6, wp), &
                                     external_fluid=use_extf, run_tmax=REAL(InitInp%Tmax, wp) )
      END IF
      IF (es /= CD_AGG_OK) THEN
         ! The slot stays inactive (it is activated only after full success, below), so this and
         ! every other early-return error path leaves it reusable by the next new_instance() --
         ! no explicit per-path release needed.
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      ! Rigid6 bodies and rods carry their own degrees of freedom, but the quasi-static
      ! linearization re-solves only the line and point equilibrium at each perturbed
      ! boundary and holds them at their deck pose, so its dYdu would silently omit
      ! their compliance. Fail closed until body-state linearization exists.
      IF (InitInp%Linearize .AND. (CD_AGG_HasRigid6( Inst(id)%agg ) .OR. CD_AGG_HasRod( Inst(id)%agg ))) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: formal OpenFAST linearization '// &
                          'with Rigid6 bodies or RODs is not supported; the quasi-static reduction '// &
                          'holds them at their deck pose instead of re-solving their equilibrium '// &
                          '(linearize a deck without BODIES/RODS).', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (InitInp%Linearize .AND. Inst(id)%agg%ncable > 0) THEN
         DO ic = 1, Inst(id)%agg%ncable
            IF (.NOT. Inst(id)%agg%cables(ic)%has_parent_endconn) CYCLE
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: formal OpenFAST linearization '// &
                             'with a platform-relative cable end connection is not supported; '// &
                             'the required orientation and connection-moment Jacobian columns '// &
                             'are not part of the current quasi-static reduction.', &
                             ErrStat, ErrMsg, RoutineName )
            RETURN
         END DO
      END IF
      CALL CD_AGG_GetInitMetadata( Inst(id)%agg, nreportlines, npoints, nsections, nei0, nfinite, es, em )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      Inst(id)%nsuper = nsuper
      Inst(id)%dt_moor = dt_moor
      Inst(id)%use_extfluid = use_extf
      ! CableDyn-owned result files. Single-turbine names use p%RootName =
      ! <OpenFASTRoot>.CD; FAST.Farm preserves its established FarmCD root.
      IF (nT > 0) THEN
         native_root = TRIM(InitInp%RootName)
      ELSE
         native_root = TRIM(InitInp%RootName)//'.CD'
      END IF
      ! The static profile records the converged static equilibrium, the state the
      ! initialization summary reports, so it is written before the SeaState kinematics
      ! are seeded below: they act from t = 0 and enter the t = 0 output row only.
      CALL CD_AGG_WriteStaticProfile( Inst(id)%agg, TRIM(native_root)//'.static.out', es, em )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (use_extf) THEN
         ! the WaveField pointer + sampling workspaces (sized once; no per-step allocation)
         p%WaveField => InitInp%WaveField
         Inst(id)%nfluid = CD_AGG_NFluidNodes( Inst(id)%agg, es, em )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         ALLOCATE (Inst(id)%fx_xyz(3, Inst(id)%nfluid), Inst(id)%fx_vel(3, Inst(id)%nfluid), &
                   Inst(id)%fx_acc(3, Inst(id)%nfluid), Inst(id)%fx_wl(Inst(id)%nfluid), &
                   Inst(id)%fx_pd(Inst(id)%nfluid), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate the ambient-fluid sampling workspace.', &
                             ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         CALL WrScr( '   CableDyn: SeaState wave/current kinematics drive the mooring hydro ('// &
                     TRIM(Num2LStr(Inst(id)%nfluid))//' sampling nodes, refreshed every mooring step).' )
         ! Seed the held fields at t = 0 (the static-IC positions): with supercycling
         ! (nsuper > 1) the first mooring advance happens only at the END of the first
         ! interval, and every CalcOutput before it would otherwise serve loads from the
         ! wave-capable-at-ZERO-fields build -- a still-water start under a declared sea
         ! state. Seeding after the static solve matches the stock-MoorDyn convention:
         ! the IC is the still-water equilibrium, wave kinematics act from t = 0.
         CALL sample_wavefield( Inst(id), p, m, 0.0_DbKi, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF

      ! The coupling mesh carries ONLY the host-prescribed (Coupled/Vessel) points.
      ! Fixed anchors stay inside the solver's boundary set, held at their initial
      ! positions -- they must never receive platform motion or return loads to it.
      ncp = CD_AGG_NMovingPoints( Inst(id)%agg, es, em )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
      END IF
      IF (Inst(id)%agg%has_sys) THEN
         nlines = CD_System_NLines( Inst(id)%agg%sys%fast%system, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
      ELSE
         nlines = 0                          ! pure finite-EI-cable deck: no EI=0 mooring system
      END IF
      IF (ncp < 1) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the deck declares no Coupled/Vessel '// &
                          'points; a mooring module needs at least one coupled fairlead.', &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      Inst(id)%ncp    = ncp
      Inst(id)%nlines = nlines

      ! ------------------------------------------------------------------------------
      ! Parameters (glue-visible bookkeeping only; the physics lives in the instance)
      ! ------------------------------------------------------------------------------
      p%RootName    = TRIM(InitInp%RootName)//'.CD'
      ! The MooringFile folder: checkpoint restart and linearization write their temporary
      ! deck copy there, so deck-relative side files resolve as in the original run.
      CALL GetPath( InitInp%FileName, p%PriPath )
      ! 0 for a single turbine, FarmSize under FAST.Farm (the checkpoint rebuild guard).
      p%nTurbines   = nT
      p%nLines      = nlines
      p%NumOuts     = 0
      p%g           = InitInp%g
      p%rhoW        = InitInp%rhoW
      p%WtrDpth     = InitInp%WtrDepth
      p%dtM0        = REAL(dt_moor, DbKi)
      p%dtCoupling  = REAL(DTcoupling, DbKi)
      p%PtfmInit    = REAL(ptfm6, ReKi)   ! the rebuild reapplies the displaced-start geometry
      p%Tmax        = InitInp%Tmax
      ! DTcoupling itself is deliberately left untouched (SubSteps = 1 in MV_AddModule).
      ! The module provides no visualization meshes: say so when the glue asks for them
      ! (WrVTK > 0) instead of leaving the mooring silently absent from the VTK output.
      IF (InitInp%VisMeshes) CALL WrScr( '   CableDyn: visualization meshes are not provided; the mooring '// &
                                         'lines and points are not written to the VTK output.' )

      ! ------------------------------------------------------------------------------
      ! Input/output meshes: one point mesh over the coupled points (all six motion
      ! fields), plus the Force/Moment sibling for the returned loads.
      ! ------------------------------------------------------------------------------
      ALLOCATE (pos(3, ncp), vel(3, ncp), acc(3, ncp), fld(3, ncp), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate coupled-point work arrays.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      CALL CD_AGG_GetMovingPointMesh( Inst(id)%agg, pos, vel, acc, fld, es, em )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
      END IF
      ! Pre-allocate the instance's persistent marching workspace here (ncp is final): the
      ! supercycled CalcOutput path reads state through these WITHOUT a mesh_to_arrays call,
      ! so they must exist before the first marching call, not grow on first use.
      ALLOCATE (Inst(id)%ws_pos(3, ncp), Inst(id)%ws_vel(3, ncp), Inst(id)%ws_acc(3, ncp), &
                Inst(id)%ws_fld(3, ncp), Inst(id)%ws_moment(3, ncp), Inst(id)%ws_orient(3, 3, ncp), &
                Inst(id)%ws_omega(3, ncp), Inst(id)%ws_alpha(3, ncp), &
                Inst(id)%sv_pos(3, ncp), Inst(id)%sv_vel(3, ncp), Inst(id)%sv_acc(3, ncp), &
                Inst(id)%sv_orient(3, 3, ncp), Inst(id)%sv_omega(3, ncp), Inst(id)%sv_alpha(3, ncp), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate the marching workspace.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (nT > 0) THEN
         ! FAST.Farm: one mesh pair per turbine over the flat aggregate moving surface;
         ! deck coupled coordinates are TURBINE-LOCAL (mesh reference positions), the
         ! per-turbine PtfmInit is carried in TranslationDisp, and the farm-layout
         ! reference offset is applied at the transfer boundary (MoorDyn's
         ! Position + TranslationDisp + TurbineRefPos convention).
         CALL build_farm_meshes( Inst(id), u, y, InitInp, nT, ncp, pos, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
      ELSE

      ALLOCATE (u%CoupledKinematics(1), y%CoupledLoads(1), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate mesh containers.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      CALL MeshCreate( BlankMesh=u%CoupledKinematics(1), IOS=COMPONENT_INPUT, Nnodes=ncp, &
                       TranslationDisp=.TRUE., Orientation=.TRUE., TranslationVel=.TRUE., &
                       RotationVel=.TRUE., TranslationAcc=.TRUE., RotationAcc=.TRUE., &
                       ErrStat=ErrStat2, ErrMess=ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! Mesh reference positions are the UNDISPLACED deck coordinates (the motion
      ! mapping anchors its lever arms on Position); the aggregate returns the coupled
      ! points at the DISPLACED, statics-solved pose, so with PtfmInit active the
      ! reference is recovered by the exact inverse rigid transform
      ! rRef = M (p - r) (p = r + M^T rRef), and the initial displacement p - rRef is
      ! carried in TranslationDisp -- MoorDyn-F's convention.
      IF (has_ptfm) THEN
         OrMatInit = EulerConstructZYX( ptfm6(4:6) )
      END IF
      CALL Eye( OrientIdent, ErrStat2, ErrMsg2 )
      DO i = 1, ncp
         IF (has_ptfm) THEN
            refpos = MATMUL( OrMatInit, REAL(pos(:, i), R8Ki) - ptfm6(1:3) )
         ELSE
            refpos = REAL(pos(:, i), R8Ki)
         END IF
         CALL MeshPositionNode( u%CoupledKinematics(1), i, REAL(refpos, ReKi), ErrStat2, ErrMsg2, &
                                Orient=OrientIdent )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL MeshConstructElement( u%CoupledKinematics(1), ELEMENT_POINT, ErrStat2, ErrMsg2, i )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL MeshCommit( u%CoupledKinematics(1), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      u%CoupledKinematics(1)%TranslationDisp = 0.0_ReKi
      u%CoupledKinematics(1)%TranslationVel  = 0.0_ReKi
      u%CoupledKinematics(1)%RotationVel     = 0.0_ReKi
      u%CoupledKinematics(1)%TranslationAcc  = 0.0_ReKi
      u%CoupledKinematics(1)%RotationAcc     = 0.0_ReKi
      IF (has_ptfm) THEN
         DO i = 1, ncp
            u%CoupledKinematics(1)%TranslationDisp(:, i) = REAL(pos(:, i), ReKi) &
               - u%CoupledKinematics(1)%Position(:, i)
            u%CoupledKinematics(1)%Orientation(:, :, i) = REAL(OrMatInit, ReKi)
         END DO
      END IF

      ! Active line control (CONTROL section): allocate the per-channel command
      ! inputs (zero = no control, bit-identical to an uncontrolled deck) and tell
      ! the glue which channels a line actually requested (CableCChanRqst drives
      ! the ServoDyn cable-control hookup, the MoorDyn convention).
      IF (CD_AGG_NCtrlChans( Inst(id)%agg ) > 0) THEN
         BLOCK
            INTEGER(IntKi) :: nchc, jc
            nchc = CD_AGG_NCtrlChans( Inst(id)%agg )
            ALLOCATE (u%DeltaL(nchc), u%DeltaLdot(nchc), InitOut%CableCChanRqst(nchc), STAT=ErrStat2)
            IF (ErrStat2 /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, 'Could not allocate the line-control channel inputs.', &
                                ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            IF (InitInp%Linearize) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: linearization with a CONTROL '// &
                                'section is not supported (the DeltaL input variables are not '// &
                                'registered); linearize an uncontrolled deck.', ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
            u%DeltaL = 0.0_ReKi
            u%DeltaLdot = 0.0_ReKi
            InitOut%CableCChanRqst = .FALSE.
            DO jc = 1, SIZE(Inst(id)%agg%ctrl_chan)
               InitOut%CableCChanRqst(Inst(id)%agg%ctrl_chan(jc)) = .TRUE.
            END DO
         END BLOCK
      END IF

      CALL MeshCopy( SrcMesh=u%CoupledKinematics(1), DestMesh=y%CoupledLoads(1), CtrlCode=MESH_SIBLING, &
                     IOS=COMPONENT_OUTPUT, Force=.TRUE., Moment=.TRUE., ErrStat=ErrStat2, ErrMess=ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      y%CoupledLoads(1)%Force  = 0.0_ReKi
      y%CoupledLoads(1)%Moment = 0.0_ReKi

      END IF   ! nT == 0 (plain single-turbine mesh path, bit-for-bit)

      ! ------------------------------------------------------------------------------
      ! Framework states: the MoorDyn-layout mirror (interior-node velocities then
      ! positions per line) + trivial xd/z/other.
      ! ------------------------------------------------------------------------------
      ALLOCATE (m%LineStateIs1(nlines), m%LineStateIsN(nlines), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate line-state index arrays.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      ! Coupled Rigid6 bodies are marched by the aggregate (CD_AGG_Step_Moving); their 6-DOF body
      ! state (pose/rates, CD_AGG_Rigid6_MirrorSize reals) is carried in the x%states checkpoint
      ! mirror below -- appended after the finite-EI cable block, following the same pattern. The body
      ! drives its attachment Free points, so restoring the 6-DOF state on restart makes the points
      ! consistent (they are re-derived from the restored pose on the next advance).
      nstates = 0
      DO il = 1, nlines
         ndofl = CD_System_Line_NDOF( Inst(id)%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         m%LineStateIs1(il) = nstates + 1
         nstates = nstates + 3*(ndofl - 6)          ! interior nodes only: [rd, r, rdd]
         ! a viscoelastic line appends its per-element dl_1 history (the series-Kelvin internal
         ! strain a restart cannot re-derive from q/v) after the interior block; this
         ! stays inside the line's LineStateIs1..IsN span (see pack_state_mirror)
         IF (CD_System_Line_Has_Viscoelastic( Inst(id)%agg%sys%fast%system, il )) THEN
            nstates = nstates + CD_System_Line_NElem( Inst(id)%agg%sys%fast%system, il, es, em )
            IF (es /= CD_SYSTEM_OK) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
         END IF
         ! a Syrope line appends TWO per-element blocks after the interior: the
         ! slow-spring static strain then the running-maximum tension (neither
         ! re-derivable from q/v). A line is viscoelastic OR Syrope, never both.
         IF (CD_System_Line_Has_Syrope( Inst(id)%agg%sys%fast%system, il )) THEN
            nstates = nstates + 2*CD_System_Line_NElem( Inst(id)%agg%sys%fast%system, il, es, em )
            IF (es /= CD_SYSTEM_OK) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
         END IF
         m%LineStateIsN(il) = nstates
      END DO
      ! finite-EI cables append their free-DOF [v_free; q_free; a_free] mirror after the mooring block
      DO ic = 1, Inst(id)%agg%ncable
         nstates = nstates + CD_HFMF_MirrorSize( Inst(id)%agg%cables(ic) )
      END DO
      ! coupled Rigid6 bodies append their 6-DOF body-state mirror after the cable block (0 if none)
      nstates = nstates + CD_AGG_Rigid6_MirrorSize( Inst(id)%agg )
      ! coupled rigid rods append the equivalent committed rigid-state mirror.
      nstates = nstates + CD_AGG_Rod_MirrorSize( Inst(id)%agg )
      ! EI=0 mooring lines on a frictional seabed append their stick-slip friction anchors,
      ! (x, y) per node, line by line (committed state a restart cannot re-derive)
      nstates = nstates + line_friction_mirror_size( Inst(id) )
      ! dynamic (Free/Connect) points append TWO trailing blocks: the compact
      ! coupled-motion vector [q; v; a] and the per-point [v; q; a] store.
      !
      ! INVARIANT (why one triple per shared point suffices): "two line ends bound to
      ! one point with different committed endpoint accelerations" is not a reachable
      ! committed state -- the dynamic-point step writes exactly one value per point
      ! into every bound endpoint each advance, and gather_mapped_motion FAILS CLOSED
      ! (1e-10) if any bound ends ever disagree. The same gather runs inside every
      ! normal step, so a state violating the invariant cannot even advance, let alone
      ! be packed. The vector's real payload is the COMMITTED boundary kinematics: the
      ! vessel-end accelerations the last advance committed differ at round-off from the
      ! restored input mesh, and on a point deck that round-off is amplified through the
      ! explicit-point window. Pure coupled decks keep the interiors-only mirror (their
      ! restart is exact from the host-restored boundary alone at output precision).
      IF (Inst(id)%agg%has_sys) THEN
         ! the committed boundary kinematics of the EI=0 lines: a restart takes the coupled
         ! endpoints' q/v from the restored input mesh, but their committed accelerations
         ! (the generalised-alpha a_n) only from here
         nstates = nstates + 3*CD_System_NSystemCoupledDOF( Inst(id)%agg%sys%fast%system )
      END IF
      IF (CD_AGG_NDynamicPoints( Inst(id)%agg ) > 0 .AND. Inst(id)%agg%has_sys) THEN
         nstates = nstates + 9*CD_AGG_NDynamicPoints( Inst(id)%agg )
         ! FAILURE rows: one fired flag per row at the mirror tail (0/1). A failure
         ! deck always carries dynamic points (its reserves are Free points), so this
         ! block lives inside the dynamic-point mirror gate.
         nstates = nstates + CD_AGG_NFailures( Inst(id)%agg )
      END IF
      m%Nx    = nstates
      m%Nxtra = nstates
      ALLOCATE (x%states(MAX(nstates, 1)), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate the state mirror.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      x%states = 0.0_DbKi
      CALL pack_state_mirror( Inst(id), m, x, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      xd%dummy    = 0.0_SiKi
      z%dummy     = 0.0_SiKi
      IF (g_serial == 0_IntKi) THEN
         ! first mint in this process: salt from the wall clock so a fresh process can
         ! never reproduce another process's serials (a -restart process would otherwise
         ! mint the same small integers the checkpointed run saved)
         CALL SYSTEM_CLOCK( COUNT=clock_count )
         g_serial = INT(MOD(ABS(clock_count), 8388608_IntKi), IntKi)*128_IntKi
      END IF
      g_serial = g_serial + 1_IntKi
      Inst(id)%serial = g_serial
      Inst(id)%t_commit = 0.0_DbKi          ! committed at the static IC (t = 0)
      other%restart_serial = g_serial
      other%t_commit = Inst(id)%t_commit

      ! ------------------------------------------------------------------------------
      ! Output channels: the deck OUTPUTS list in the OrcaFlex vocabulary, emitted through
      ! the OpenFAST WriteOutput surface (per-line tensions, curvature, fairlead angle, ...).
      ! The channels were validated when the aggregate was built (parse_deck routes each
      ! OUTPUTS token through check_channel against the deck lines/points and rejects an
      ! unknown token or an out-of-range line node with a named fatal). Here we size the
      ! framework arrays, label every channel with its header token + OrcaFlex unit, and
      ! evaluate each once at the static equilibrium -- which also confirms the whole eval
      ! path resolves (a channel the parser admits but the eval cannot serve, e.g. a Point
      ! channel with no mooring system, fails closed HERE rather than at the first output
      ! step). A deck with no OUTPUTS gives nch = 0, size-zero arrays, and behaviour
      ! byte-identical to the pre-channel module (p%NumOuts = 0, no CableDyn columns).
      ! ------------------------------------------------------------------------------
      nch = CD_AGG_NumChannels( Inst(id)%agg )
      p%NumOuts = nch
      ALLOCATE (y%WriteOutput(nch), InitOut%WriteOutputHdr(nch), InitOut%WriteOutputUnt(nch), &
                Inst(id)%ws_chan(nch), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate the output channel arrays.', &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      DO ic = 1, nch
         CALL CD_AGG_ChannelHeader( Inst(id)%agg, ic, chan_hdr, chan_unt, es, em )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         InitOut%WriteOutputHdr(ic) = chan_hdr
         InitOut%WriteOutputUnt(ic) = chan_unt
         CALL CD_AGG_EvalChannel( Inst(id)%agg, ic, cval, es, em )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: output channel "'//TRIM(chan_hdr)// &
                             '" could not be evaluated at the static equilibrium: '//TRIM(em), &
                             ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         y%WriteOutput(ic) = REAL(cval, ReKi)
      END DO

      CALL CD_InitVars( InitOut%Vars, u, p, x, y, m, Inst(id), InitInp%Linearize, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! The time history contains the selected deck channels and is therefore opened only
      ! when nch > 0. The static profile, unconditional and covering every line node
      ! regardless of OUTPUTS selection, was written at the equilibrium above.
      CALL open_native_outputs( id, TRIM(native_root), 0.0_DbKi, ErrStat2, ErrMsg2, write_profile=.FALSE. )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      IF (profile_requested) THEN
         CALL SYSTEM_CLOCK(prof_c1)
         CALL reset_instance_profile( id, .TRUE., &
                                      REAL(prof_c1 - prof_c0, wp)/REAL(prof_rate, wp), .TRUE. )
         CALL WrScr( '    CableDyn profiler armed by CABLEDYN_PROFILE.' )
      ELSE
         CALL reset_instance_profile( id, .FALSE., 0.0_wp, .TRUE. )
      END IF

      CALL WrScr( '   Created CableDyn model: '//TRIM(Num2LStr(nreportlines))//' line object(s), '// &
                  TRIM(Num2LStr(npoints))//' point(s), '//TRIM(Num2LStr(nsections))//' section(s) [EI=0: '// &
                  TRIM(Num2LStr(nei0))//', finite-EI: '//TRIM(Num2LStr(nfinite))//'].' )
      CALL WrScr( '   Initial conditions: Newton static equilibrium with load continuation completed.' )
      IF (use_extf) THEN
         ! The equilibrium is solved at zero SeaState fields and the SeaState kinematics act
         ! from t = 0, so the t = 0 outputs below differ from these equilibrium values by
         ! the hydrodynamic load of the t = 0 kinematics on the line end nodes.
         CALL WrScr( '   Line results below are at this equilibrium; SeaState kinematics act from t = 0.' )
      END IF
      CALL WrScr( '   Fairlead convention: force is on End A toward End B; inclinations are signed below horizontal.' )
      DO il = 1, nreportlines
         CALL CD_AGG_GetInitLine( Inst(id)%agg, il, line_id, fair_ten, fair_force, fair_incl, fair_decl, &
                                  fair_azi, es, em )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         CALL WrScr( '   Line '//TRIM(Num2LStr(line_id))//' fairlead effective tension: '// &
                     TRIM(Num2LStr(REAL(fair_ten, DbKi)))//' N' )
         ! The force carries the end node's share of the line load (and, finite-EI, the end
         ! shear), so its angle differs slightly from the line tangent's: both are printed.
         force_incl = ATAN2(-fair_force(3), NORM2(fair_force(1:2)))*REAL(R2D_D, wp)
         CALL WrScr( '      force [Fx, Fy, Fz]: ['//TRIM(Num2LStr(REAL(fair_force(1), DbKi)))//', '// &
                     TRIM(Num2LStr(REAL(fair_force(2), DbKi)))//', '// &
                     TRIM(Num2LStr(REAL(fair_force(3), DbKi)))//'] N, inclination='// &
                     TRIM(Num2LStr(REAL(force_incl, DbKi)))//' deg' )
         CALL WrScr( '      line tangent: inclination='//TRIM(Num2LStr(REAL(fair_incl, DbKi)))// &
                     ' deg, declination='//TRIM(Num2LStr(REAL(fair_decl, DbKi)))// &
                     ' deg, azimuth='//TRIM(Num2LStr(REAL(fair_azi, DbKi)))//' deg' )
      END DO
      CALL WrScr( '  CableDyn initialization completed.' )

      ! Initialization fully succeeded: bind the handle and activate the slot together. Until this
      ! point the slot was reserved but inactive, so any early-return error path above (including a
      ! failed report query) leaves it reusable rather than leaking a half-built active instance.
      ! Crucially m%CDInst is also left at its INTENT(OUT) default (0) on those failure paths, so a
      ! later init can reuse this slot without a failed module record aliasing the new instance id.
      m%CDInst        = id
      Inst(id)%active = .TRUE.

      ! In a single-turbine OpenFAST run, expose the converged t = 0 result immediately.
      ! These are the same values already placed in y%WriteOutput for the first host
      ! output row; printing them here makes a direct executable smoke test useful even
      ! before the user opens the OpenFAST .out/.outb file. FAST.Farm is excluded to
      ! avoid repeating a potentially long channel list for every turbine instance; its
      ! static row is written to <RootName>.FarmCD.out above.
      IF (nT == 0) THEN
         IF (nch > 0) THEN
            IF (use_extf) THEN
               CALL WrScr( '    Requested CableDyn OUTPUTS at t = 0 s (equilibrium pose, SeaState kinematics '// &
                           'at t = 0):' )
            ELSE
               CALL WrScr( '    Requested CableDyn OUTPUTS at t = 0 s (static equilibrium):' )
            END IF
            DO ic = 1, nch
               CALL WrScr( '      '//TRIM(InitOut%WriteOutputHdr(ic))//' = '// &
                           TRIM(Num2LStr(REAL(y%WriteOutput(ic), DbKi)))//' '// &
                           TRIM(InitOut%WriteOutputUnt(ic)) )
            END DO
            CALL WrScr( '    These values are also the t = 0 CableDyn columns in the OpenFAST output file.' )
         ELSE
            CALL WrScr( '    CableDyn: no deck OUTPUTS requested; no CableDyn result columns will be written.' )
         END IF
      END IF

   END SUBROUTINE CD_Init

!----------------------------------------------------------------------------------------------------------------------------------
   LOGICAL FUNCTION instance_matches( m, other ) RESULT(ok)
      !! TRUE iff m%CDInst refers to a live instance whose identity matches the
      !! (possibly checkpoint-restored) OtherState: same generation stamp AND a
      !! committed clock no further ahead than the restored one. A fresh process fails
      !! the stamp; a same-process rewind fails the clock; both route to
      !! rebuild_from_checkpoint.
      TYPE(CD_MiscVarType),    INTENT(IN) :: m
      TYPE(CD_OtherStateType), INTENT(IN) :: other
      ok = .FALSE.
      IF (.NOT. ALLOCATED(Inst)) RETURN
      IF (m%CDInst < 1 .OR. m%CDInst > SIZE(Inst)) RETURN
      IF (.NOT. Inst(m%CDInst)%active) RETURN
      IF (Inst(m%CDInst)%serial /= other%restart_serial) RETURN
      ! EQUALITY, not merely not-ahead: a fresh instance (t_commit = 0) must never be
      ! mistaken for a checkpointed one, even if the process-local serial collides.
      IF (ABS(Inst(m%CDInst)%t_commit - other%t_commit) > 1.0E-9_DbKi*MAX(1.0_DbKi, ABS(other%t_commit))) RETURN
      ok = .TRUE.
   END FUNCTION instance_matches

   SUBROUTINE rebuild_from_checkpoint( t, u1, p, x, other, m, ErrStat, ErrMsg )
      !! CHECKPOINT-RESTART: rebuild the module's internal solver instance from the
      !! registry-packed data alone -- the deck copy in p, the state mirror in x, the
      !! coupled mesh in u -- exactly the data an `openfast -restart` run unpacks. The
      !! deck is re-parsed and the static build re-run (structure, meshes, hydro
      !! configuration), then the committed dynamic state is overlaid: coupled endpoint
      !! kinematics from the input mesh, line interiors from the mirror with the
      !! acceleration re-derived (the EI=0 mirror is MoorDyn-layout [rd; r]), and each
      !! cable's [v; q; a] free-DOF mirror reloaded exactly with the restart clock.
      REAL(DbKi),                   INTENT(IN   ) :: t
      TYPE(CD_InputType),           INTENT(IN   ) :: u1
      TYPE(CD_ParameterType),       INTENT(IN   ) :: p
      TYPE(CD_ContinuousStateType), INTENT(IN   ) :: x
      TYPE(CD_OtherStateType),      INTENT(IN   ) :: other
      TYPE(CD_MiscVarType),         INTENT(INOUT) :: m
      INTEGER(IntKi),               INTENT(  OUT) :: ErrStat
      CHARACTER(*),                 INTENT(  OUT) :: ErrMsg

      character(*), parameter   :: RoutineName = 'CableDyn:rebuild_from_checkpoint'
      CHARACTER(2100)           :: rst_deck
      INTEGER(IntKi)            :: id, old_id, stale_id, i, es, ncp, nlines, il, ndofl, nin, k, ic, csz, koff
      LOGICAL                   :: use_extf, profile_requested, preserve_profile
      CHARACTER(512)            :: em
      INTEGER(IntKi)            :: ErrStat2
      CHARACTER(ErrMsgLen)      :: ErrMsg2

      ErrStat = ErrID_None
      ErrMsg  = ''
      IF (p%nTurbines > 0_IntKi) THEN
         ! FAST.Farm restart follows FAST.Farm's own checkpoint protocol, which CableDyn does not implement.
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: checkpoint restart in FAST.Farm '// &
                          'mode is not supported.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (p%DeckCopy%NumLines < 1) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: checkpoint rebuild requested but the '// &
                          'parameters carry no deck copy.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      ! A checkpoint written mid-mooring-interval carries interiors from the LAST
      ! committed solve (other%t_commit) while the glue would overlay the boundary at
      ! the newer t -- a state no continuous run ever had. Fail closed: restarts are
      ! valid only on committed interval boundaries (choose ChkptTime as a multiple
      ! of the mooring step dtM).
      IF (ABS(REAL(t, DbKi) - other%t_commit) > 1.0E-6_DbKi*MAX(1.0_DbKi, ABS(other%t_commit))) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the checkpoint was written at t = '// &
                          TRIM(Num2LStr(t))//' s, between mooring commits (last commit at '// &
                          TRIM(Num2LStr(other%t_commit))//' s). Choose ChkptTime as a multiple of the '// &
                          'mooring step dtM so the restored state is a committed interval boundary.', &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      ! the deck copy back to disk (the parser's input boundary is a file), in the
      ! MooringFile folder so deck-relative side files resolve as in the original run
      CALL write_deck_copy( p, '.rst.dat', rst_deck, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! Profiling is opt-in per host run, while the Hermite profiler is process-global.
      ! A same-process rewind should retain its accumulated diagnostics. A fresh
      ! restart (or a host that removed CABLEDYN_PROFILE between sequential cases)
      ! must explicitly arm/reset or disarm the global timer before rebuilding.
      old_id = m%CDInst
      profile_requested = profile_environment_requested()
      preserve_profile = .FALSE.
      IF (ALLOCATED(Inst)) THEN
         IF (old_id >= 1_IntKi .AND. old_id <= SIZE(Inst)) THEN
            IF (Inst(old_id)%active) THEN
               preserve_profile = Inst(old_id)%prof_enabled .AND. profile_requested .AND. &
                                  Inst(old_id)%serial == other%restart_serial
            END IF
         END IF
      END IF
      id = new_instance()
      IF (preserve_profile) THEN
         CALL copy_instance_profile( id, old_id )
      ELSE
         CALL reset_instance_profile( id, profile_requested, 0.0_wp, .TRUE. )
      END IF
      use_extf = .FALSE.
      IF (ASSOCIATED(p%WaveField)) use_extf = wavefield_is_ambient( p%WaveField )
      CALL CD_AGG_Init_From_Deck( Inst(id)%agg, TRIM(rst_deck), REAL(p%dtM0, wp), es, em, &
                                  env_gravity=REAL(p%g, wp), env_rho_water=REAL(p%rhoW, wp), &
                                  env_wtrdpth=REAL(p%WtrDpth, wp), external_fluid=use_extf, &
                                  ptfm_init=REAL(p%PtfmInit, wp), run_tmax=REAL(p%Tmax, wp) )
      CALL delete_file( rst_deck )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: checkpoint rebuild: '//TRIM(em), &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      ! structural consistency with the restored registry data (a different deck or a
      ! corrupted checkpoint must fail closed, not silently mis-map states)
      ncp = CD_AGG_NMovingPoints( Inst(id)%agg, es, em )
      IF (es /= CD_AGG_OK .OR. ncp /= u1%CoupledKinematics(1)%Nnodes) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: checkpoint rebuild: the rebuilt deck '// &
                          'has a different coupled-point count than the restored mesh.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (Inst(id)%agg%has_sys) THEN
         nlines = CD_System_NLines( Inst(id)%agg%sys%fast%system, es, em )
      ELSE
         nlines = 0
      END IF
      IF (nlines /= p%nLines) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: checkpoint rebuild: the rebuilt deck '// &
                          'has a different line count than the restored parameters.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      ! Dynamic (Free/Connect) point states ride two trailing mirror blocks: the compact
      ! coupled-motion vector (the committed endpoint/boundary kinematics -- the vessel
      ! ends' committed accelerations differ at round-off from the restored input mesh,
      ! which a point deck amplifies) and the per-point store. A deck whose layout
      ! changed since the checkpoint fails closed on the mirror-extent checks.

      Inst(id)%ncp    = ncp
      Inst(id)%nlines = nlines
      Inst(id)%nsuper = MAX(1_IntKi, NINT(p%dtM0/p%dtCoupling, IntKi))
      Inst(id)%dt_moor = REAL(p%dtM0, wp)
      ALLOCATE (Inst(id)%ws_pos(3, ncp), Inst(id)%ws_vel(3, ncp), Inst(id)%ws_acc(3, ncp), &
                Inst(id)%ws_fld(3, ncp), Inst(id)%ws_moment(3, ncp), Inst(id)%ws_orient(3, 3, ncp), &
                Inst(id)%ws_omega(3, ncp), Inst(id)%ws_alpha(3, ncp), &
                Inst(id)%sv_pos(3, ncp), Inst(id)%sv_vel(3, ncp), Inst(id)%sv_acc(3, ncp), &
                Inst(id)%sv_orient(3, 3, ncp), Inst(id)%sv_omega(3, ncp), Inst(id)%sv_alpha(3, ncp), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate the marching workspace (rebuild).', &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      Inst(id)%use_extfluid = use_extf
      IF (use_extf) THEN
         Inst(id)%nfluid = CD_AGG_NFluidNodes( Inst(id)%agg, es, em )
         ALLOCATE (Inst(id)%fx_xyz(3, Inst(id)%nfluid), Inst(id)%fx_vel(3, Inst(id)%nfluid), &
                   Inst(id)%fx_acc(3, Inst(id)%nfluid), Inst(id)%fx_wl(Inst(id)%nfluid), &
                   Inst(id)%fx_pd(Inst(id)%nfluid), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate the ambient-fluid workspace (rebuild).', &
                             ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      END IF
      ! the WriteOutput channel buffer (CalcOutput evaluates the deck channels into it
      ! every reporting step; Init sizes it from the aggregate's channel count)
      IF (CD_AGG_NumChannels(Inst(id)%agg) > 0) THEN
         ALLOCATE (Inst(id)%ws_chan(CD_AGG_NumChannels(Inst(id)%agg)), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate the channel workspace (rebuild).', &
                             ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      END IF

      ! overlay the committed dynamic state: coupled endpoints from the restored mesh...
      CALL mesh_to_arrays( u1%CoupledKinematics(1), ncp, Inst(id)%ws_pos, Inst(id)%ws_vel, &
                           Inst(id)%ws_acc, Inst(id)%ws_orient, Inst(id)%ws_omega, Inst(id)%ws_alpha, &
                           ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      CALL CD_AGG_UpdateStates_Moving( Inst(id)%agg, Inst(id)%ws_pos, Inst(id)%ws_vel, Inst(id)%ws_acc, es, em, &
                                       orientation=Inst(id)%ws_orient, angular_velocity=Inst(id)%ws_omega, &
                                       angular_acceleration=Inst(id)%ws_alpha )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: checkpoint rebuild: '//TRIM(em), &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      ! ...line interiors + cable free-DOF mirrors, exactly (committed accelerations
      ! reload from the mirror -- shared helper with the state-const CalcOutput)
      CALL reload_interiors_ary( Inst(id), m, x%states, REAL(t, wp), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! restore the ambient-fluid fields the reference's committed state carried: the
      ! last pre-checkpoint solve sampled the WaveField at ITS interval end (= the
      ! restart time), so sampling here reproduces those fields on the restored
      ! positions. Without this the first post-restart load evaluations would see still
      ! water and ring for a few steps.
      Inst(id)%active = .TRUE.       ! sample_wavefield reads through the live instance
      IF (use_extf) THEN
         CALL sample_wavefield( Inst(id), p, m, t, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL release_instance( id )
            RETURN
         END IF
      END IF

      ! Active line control: the rebuild re-initialized at the BASE segment lengths;
      ! apply the restored command HERE, on the rebuild transition itself, so every
      ! entry that can trigger a rebuild (UpdateStates -- including its supercycle
      ! early-return substeps -- and CalcOutput) resumes at the checkpointed lengths.
      IF (ALLOCATED(u1%DeltaL)) THEN
         CALL CD_AGG_Apply_LineControl( Inst(id)%agg, REAL(u1%DeltaL, wp), REAL(u1%DeltaLdot, wp), es, em )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: line-control apply after the '// &
                             'checkpoint rebuild failed: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
            CALL release_instance( id )
            RETURN
         END IF
      END IF
      ! A fresh-process restart first passes through CD_Init, whose provisional instance
      ! owns the output stream of the same root; the checkpointed m%CDInst need not be
      ! that slot. Find it so its stream is released and the slot freed with old_id.
      stale_id = 0_IntKi
      DO k = 1, SIZE(Inst)
         IF (k == id .OR. k == old_id) CYCLE
         IF (Inst(k)%active .AND. Inst(k)%native_un > 0 .AND. &
             TRIM(Inst(k)%native_root) == TRIM(p%RootName)) THEN
            stale_id = k
            EXIT
         END IF
      END DO
      IF (stale_id > 0_IntKi) THEN
         CLOSE (Inst(stale_id)%native_un)
         Inst(stale_id)%native_un = -1_IntKi
         Inst(stale_id)%native_root = ''
      END IF
      ! A same-process replacement adopts the run-level stream from its predecessor.
      ! Any rows beyond the restored commit are removed before the open unit changes
      ! owner. A fresh-process restart has no live unit and starts a new result file at
      ! the checkpoint state, matching the normal OpenFAST restart-output policy.
      CALL restore_native_outputs( id, old_id, TRIM(p%RootName), other%t_commit, other%restart_serial, &
                                   p%dtM0, Inst(id)%nsuper, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL release_instance( id )
         RETURN
      END IF
      ! The rebuilt instance supersedes the one the checkpoint replaced (and the provisional
      ! instance of a fresh process): release their aggregates and workspaces.
      IF (old_id >= 1_IntKi .AND. old_id <= SIZE(Inst) .AND. old_id /= id) CALL release_instance( old_id )
      IF (stale_id > 0_IntKi) CALL release_instance( stale_id )
      ! adopt the restored identity (the OtherState is authoritative across the restart)
      Inst(id)%serial = other%restart_serial
      Inst(id)%t_commit = other%t_commit
      Inst(id)%last_n = -1_IntKi
      Inst(id)%last_solve_n = -1_IntKi
      Inst(id)%loads_cached = .FALSE.
      m%CDInst = id
      CALL WrScr( '   CableDyn: rebuilt the solver instance from the checkpoint at t = '// &
                  TRIM(Num2LStr(t))//' s ('//TRIM(Num2LStr(ncp))//' coupled points, '// &
                  TRIM(Num2LStr(nlines))//' mooring lines, '// &
                  TRIM(Num2LStr(Inst(id)%agg%ncable))//' cables).' )
   END SUBROUTINE rebuild_from_checkpoint

   SUBROUTINE CD_UpdateStates( t, n, u, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg)
      REAL(DbKi)                      , INTENT(IN   ) :: t
      INTEGER(IntKi)                  , INTENT(IN   ) :: n
      TYPE(CD_InputType)              , INTENT(INOUT) :: u(:)
      REAL(DbKi)                      , INTENT(IN   ) :: t_array(:)
      TYPE(CD_ParameterType)          , INTENT(INOUT) :: p
      TYPE(CD_ContinuousStateType)    , INTENT(INOUT) :: x
      TYPE(CD_DiscreteStateType)      , INTENT(INOUT) :: xd
      TYPE(CD_ConstraintStateType)    , INTENT(INOUT) :: z
      TYPE(CD_OtherStateType)         , INTENT(INOUT) :: other
      TYPE(CD_MiscVarType)            , INTENT(INOUT) :: m
      INTEGER(IntKi)                  , INTENT(  OUT) :: ErrStat
      CHARACTER(*)                    , INTENT(  OUT) :: ErrMsg

      character(*), parameter   :: RoutineName = 'CD_UpdateStates'
      INTEGER(IntKi)            :: ErrStat2, es, id, niter, iu_t, i_u
      INTEGER(IntKi)            :: prof_c0, prof_c1, prof_rate, prof_phase0
      CHARACTER(ErrMsgLen)      :: ErrMsg2
      CHARACTER(512)            :: em
      LOGICAL                   :: conv, stalled

      ErrStat = ErrID_None
      ErrMsg  = ''
      ! CHECKPOINT-RESTART: a restored OtherState that does not match the live instance
      ! (fresh process, or a same-process rewind) rebuilds the solver from p + x + u.
      ! The overlay must use the input AT the restart time t -- u(1) here is the glue's
      ! PREDICTED input for the upcoming interval (t_array(1) = t + DT), one glue step
      ! AHEAD of the committed boundary; overlaying it would start the restored run from
      ! a boundary it never had. Pick the u(:) entry whose time is t.
      IF (.NOT. instance_matches(m, other)) THEN
         iu_t = 1_IntKi
         DO i_u = 1, SIZE(t_array)
            IF (ABS(t_array(i_u) - t) <= 1.0E-6_DbKi*MAX(1.0_DbKi, ABS(t))) THEN
               iu_t = i_u
               EXIT
            END IF
         END DO
         CALL rebuild_from_checkpoint( t, u(iu_t), p, x, other, m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF
      id = m%CDInst
      IF (Inst(id)%prof_enabled) THEN
         CALL SYSTEM_CLOCK(prof_c0, prof_rate)
         Inst(id)%prof_n_update = Inst(id)%prof_n_update + 1_IntKi
      END IF

      ! Correction iterations (NumCrctn > 0): the glue re-runs UpdateStates over the SAME
      ! interval with corrected inputs. A re-entry over a glue step where nothing was
      ! advanced (an intermediate supercycle call) is a no-op; a re-entry over a COMPLETED
      ! mooring step rewinds the aggregate to the snapshot taken at that step's start and
      ! falls through to re-advance with the corrected inputs -- exactly the
      ! copy-states-back semantics a registry-state module gets from the glue, expressed
      ! on the module's internal state.
      IF (n == Inst(id)%last_n) THEN
         IF (Inst(id)%last_solve_n /= n) THEN
            IF (Inst(id)%prof_enabled) THEN
               CALL SYSTEM_CLOCK(prof_c1)
               Inst(id)%prof_t_update = Inst(id)%prof_t_update + &
                  REAL(prof_c1 - prof_c0, wp)/REAL(prof_rate, wp)
            END IF
            RETURN
         END IF
         CALL CD_AGG_Restore( Inst(id)%agg, es, em )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: correction-iteration rewind failed: '// &
                             TRIM(em), ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      END IF

      ! Supercycling (CableDyn's own dtM): the mooring interval [.., (n+1)*DT] completes only
      ! on every nsuper-th glue step; an intermediate call records the glue step and returns
      ! without an implicit solve. CD_CalcOutput serves the last committed load state between
      ! mooring solves.
      IF (MOD(n + 1_IntKi, Inst(id)%nsuper) /= 0_IntKi) THEN
         Inst(id)%last_n = n
         IF (Inst(id)%prof_enabled) THEN
            CALL SYSTEM_CLOCK(prof_c1)
            Inst(id)%prof_t_update = Inst(id)%prof_t_update + &
               REAL(prof_c1 - prof_c0, wp)/REAL(prof_rate, wp)
         END IF
         RETURN
      END IF

      ! Interpolate the input to the END of the coupling interval (= the end of the mooring
      ! interval on this stepping call): the implicit step prescribes the coupled-point
      ! kinematics at t + dt. The interpolation container is built ONCE (a full mesh copy)
      ! and reused for the rest of the run -- ExtrapInterp overwrites its fields in place --
      ! so the stepping path never creates/destroys a mesh.
      IF (.NOT. Inst(id)%u_interp_ready) THEN
         CALL CD_CopyInput( u(1), Inst(id)%u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         Inst(id)%u_interp_ready = .TRUE.
      END IF
      CALL CD_Input_ExtrapInterp( u, t_array, Inst(id)%u_interp, t + p%dtCoupling, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat < AbortErrLev) THEN

         IF (Inst(id)%nT > 0) THEN
            CALL farm_gather_kinematics( Inst(id), Inst(id)%u_interp, ErrStat2, ErrMsg2 )
         ELSE
            CALL mesh_to_arrays( Inst(id)%u_interp%CoupledKinematics(1), Inst(id)%ncp, Inst(id)%ws_pos, &
                                 Inst(id)%ws_vel, Inst(id)%ws_acc, Inst(id)%ws_orient, Inst(id)%ws_omega, &
                                 Inst(id)%ws_alpha, ErrStat2, ErrMsg2 )
         END IF
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat < AbortErrLev) THEN
            ! Snapshot the step-start state before advancing: a later correction re-entry
            ! over this same n restores it and re-advances. Idempotent (buffers reused),
            ! and on the corrected re-advance it re-captures the SAME restored state.
            CALL CD_AGG_Snapshot( Inst(id)%agg, es, em )
            IF (es /= CD_AGG_OK) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: step-start snapshot failed: '// &
                                TRIM(em), ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
            ! SeaState ambient fluid: sample the WaveField at the line nodes' committed
            ! positions, at the END of this mooring interval (where the implicit step
            ! prescribes the boundary), and hold the fields over the step. The advance
            ! runs on the LAST glue substep of the interval, so the interval end is
            ! t + p%dtCoupling -- the SAME time the coupled-kinematics interpolation
            ! above uses (t + dt_moor from here would sit almost a full mooring step in
            ! the future under supercycling, phase-shifting the wave kinematics).
            ! Runs AFTER the snapshot: a correction re-entry restores the step-start
            ! state and re-samples at the same committed positions.
            IF (Inst(id)%use_extfluid) THEN
               IF (Inst(id)%prof_enabled) CALL SYSTEM_CLOCK(prof_phase0)
               CALL sample_wavefield( Inst(id), p, m, t + p%dtCoupling, ErrStat2, ErrMsg2 )
               IF (Inst(id)%prof_enabled) THEN
                  CALL SYSTEM_CLOCK(prof_c1)
                  Inst(id)%prof_n_wave = Inst(id)%prof_n_wave + 1_IntKi
                  Inst(id)%prof_t_wave = Inst(id)%prof_t_wave + &
                     REAL(prof_c1 - prof_phase0, wp)/REAL(prof_rate, wp)
               END IF
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN
            END IF
            ! Active line control: apply the interpolated DeltaL/DeltaLdot command to
            ! the controlled segments BEFORE the advance (the MoorDyn placement).
            ! Runs AFTER the snapshot so a correction re-entry restores the
            ! step-start lengths and re-applies the corrected command.
            IF (ALLOCATED(Inst(id)%u_interp%DeltaL)) THEN
               CALL CD_AGG_Apply_LineControl( Inst(id)%agg, REAL(Inst(id)%u_interp%DeltaL, wp), &
                                              REAL(Inst(id)%u_interp%DeltaLdot, wp), es, em )
               IF (es /= CD_AGG_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: line-control apply failed: '// &
                                   TRIM(em), ErrStat, ErrMsg, RoutineName )
                  RETURN
               END IF
            END IF
            IF (Inst(id)%prof_enabled) CALL SYSTEM_CLOCK(prof_phase0)
            CALL CD_AGG_Step_Moving( Inst(id)%agg, Inst(id)%dt_moor, Inst(id)%ws_pos, Inst(id)%ws_vel, &
                                     Inst(id)%ws_acc, conv, stalled, niter, es, em, &
                                     t_committed=REAL(t + p%dtCoupling, wp), orientation=Inst(id)%ws_orient, &
                                     angular_velocity=Inst(id)%ws_omega, angular_acceleration=Inst(id)%ws_alpha )
            IF (Inst(id)%prof_enabled) THEN
               CALL SYSTEM_CLOCK(prof_c1)
               Inst(id)%prof_n_solve = Inst(id)%prof_n_solve + 1_IntKi
               Inst(id)%prof_t_step = Inst(id)%prof_t_step + &
                  REAL(prof_c1 - prof_phase0, wp)/REAL(prof_rate, wp)
            END IF
            IF (es /= CD_AGG_OK .OR. .NOT. conv) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: implicit step did not converge at t = '// &
                                TRIM(Num2LStr(t))//' s: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
            ELSE
               Inst(id)%last_n = n
               Inst(id)%last_solve_n = n
               ! The committed state changed: a served ZOH cache keyed on last_solve_n
               ! alone would go stale on a correction re-advance at the SAME n, so the
               ! cache is invalidated on every successful advance (bit-identical for the
               ! normal path -- last_solve_n moved anyway -- and correct for corrections).
               Inst(id)%loads_cached = .FALSE.
               Inst(id)%native_pending = Inst(id)%native_un > 0_IntKi
               ! commit the restart clock: the state now committed corresponds to the END
               ! of this mooring interval. The advance runs on the LAST glue substep, so
               ! the interval end is t + dtCoupling (t + dt_moor from here overshoots by
               ! dt_moor - DT under supercycling -- the same arithmetic as the WaveField
               ! sampling time).
               Inst(id)%t_commit = t + p%dtCoupling
               other%t_commit = Inst(id)%t_commit
               CALL pack_state_mirror( Inst(id), m, x, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            END IF
         END IF

      END IF

      IF (Inst(id)%prof_enabled) THEN
         CALL SYSTEM_CLOCK(prof_c1)
         Inst(id)%prof_t_update = Inst(id)%prof_t_update + &
            REAL(prof_c1 - prof_c0, wp)/REAL(prof_rate, wp)
      END IF

   END SUBROUTINE CD_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE sample_wavefield( inst_e, p, m, t_sample, ErrStat, ErrMsg )
      !! Sample the SeaState WaveField at every mooring line node's committed position --
      !! fluid velocity (waves + current), fluid acceleration, and the local free-surface
      !! elevation -- and prescribe them to the aggregate for the coming advance. Nodes
      !! are forced in-water for kinematics (forceNodeInWater = .TRUE., the MoorDyn
      !! convention: the WETTING is decided by CableDyn's own waterline machinery from
      !! the sampled elevation, not by the sampler).
      TYPE(CD_OF_Instance),     INTENT(INOUT) :: inst_e
      TYPE(CD_ParameterType),   INTENT(IN   ) :: p
      TYPE(CD_MiscVarType),     INTENT(INOUT) :: m
      REAL(DbKi),               INTENT(IN   ) :: t_sample
      INTEGER(IntKi),           INTENT(  OUT) :: ErrStat
      CHARACTER(*),             INTENT(  OUT) :: ErrMsg

      character(*), parameter   :: RoutineName = 'sample_wavefield'
      INTEGER(IntKi)            :: i, l, es, nodeInWater, n_current_nodes
      REAL(SiKi)                :: zeta1, zeta2, zeta, pdyn, fv(3), fa(3), famcf(3)
      REAL(SiKi)                :: wave_zeta1, wave_zeta2, wave_zeta, wave_pdyn
      REAL(SiKi)                :: wave_fv(3), wave_fa(3), wave_famcf(3)
      REAL(wp)                  :: profile_fv(3), steady_fv_wp(3), grid_z
      REAL(SiKi)                :: steady_fv(3), grid_fv(3), grid_weight
      CHARACTER(512)            :: em
      INTEGER(IntKi)            :: ErrStat2
      CHARACTER(ErrMsgLen)      :: ErrMsg2

      ErrStat = ErrID_None
      ErrMsg  = ''
      IF (.NOT. ASSOCIATED(p%WaveField)) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the WaveField pointer is null but '// &
                          'SeaState kinematics were requested at Init.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      CALL CD_AGG_GetFluidNodePositions( inst_e%agg, inst_e%fx_xyz, es, em )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF ((inst_e%agg%host_wave_enabled .NEQV. inst_e%agg%host_current_enabled) .AND. &
          p%WaveField%Current_InitInput%CurrMod == 2_IntKi) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: WaterKin requests separate host wave/current '// &
                          'modes, but SeaState CurrMod 2 is a private user current and cannot be decomposed.', &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      DO i = 1, inst_e%nfluid
         IF (.NOT. inst_e%agg%host_wave_enabled .AND. .NOT. inst_e%agg%host_current_enabled) THEN
            ! An explicit WaterKin file may disable both host components even
            ! though the turbine owns a SeaState field. Do not sample it at all.
            fv = 0.0_SiKi
            fa = 0.0_SiKi
            zeta = 0.0_SiKi
            nodeInWater = 0_IntKi
         ELSE
            CALL WaveField_GetNodeWaveKin( p%WaveField, m%WaveField_m, t_sample, &
                                           REAL(inst_e%fx_xyz(:, i), ReKi), .TRUE., &
                                           inst_e%agg%host_current_enabled, nodeInWater, &
                                           zeta1, zeta2, zeta, pdyn, fv, fa, famcf, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) RETURN
         END IF
         steady_fv = 0.0_SiKi
         IF (nodeInWater /= 0_IntKi .AND. &
             (.NOT. inst_e%agg%host_wave_enabled .OR. .NOT. inst_e%agg%host_current_enabled)) THEN
            ! Reconstruct the current exactly where WaveField_GetNodeWaveKin
            ! samples it. For vertical/extrapolation stretching above SWL that
            ! is the single z=0 value. SeaState deliberately does NOT add
            ! PCurrV/Pz to PWaveVel0 (see Waves.f90), so the extrapolation term
            ! is wave-only and must not be subtracted as steady current here.
            ! Below SWL and for Wheeler stretching, retain the same four
            ! vertical grid nodes and interpolation weights as the sampler.
            n_current_nodes = 4
            IF (p%WaveField%WaveStMod < 3_IntKi .AND. inst_e%fx_xyz(3, i) > 0.0_wp) &
               n_current_nodes = 1
            DO l = 1, n_current_nodes
               IF (n_current_nodes == 1) THEN
                  grid_z = 0.0_wp
                  grid_weight = 1.0_SiKi
               ELSE
                  grid_z = REAL(p%WaveField%GridDepth, wp)* &
                           (SIN(REAL(p%WaveField%VolGridParams%pZero(4) + &
                                     m%WaveField_m%Indx(l, 4)*p%WaveField%VolGridParams%delta(4), wp)) - 1.0_wp)
                  grid_weight = SUM(m%WaveField_m%N4D(:, :, :, l))
               END IF
               CALL CD_SeaState_Steady_Current(grid_z, REAL(p%WaveField%EffWtrDpth, wp), &
                    INT(p%WaveField%Current_InitInput%CurrMod), &
                    REAL(p%WaveField%Current_InitInput%CurrSSV0, wp), &
                    REAL(p%WaveField%Current_InitInput%CurrSSDir, wp), &
                    REAL(p%WaveField%Current_InitInput%CurrNSRef, wp), &
                    REAL(p%WaveField%Current_InitInput%CurrNSV0, wp), &
                    REAL(p%WaveField%Current_InitInput%CurrNSDir, wp), &
                    REAL(p%WaveField%Current_InitInput%CurrDIV, wp), &
                    REAL(p%WaveField%Current_InitInput%CurrDIDir, wp), steady_fv_wp, es, em)
               IF (es /= CD_HYDRO_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: SeaState current decomposition: '// &
                                   TRIM(em), ErrStat, ErrMsg, RoutineName )
                  RETURN
               END IF
               grid_fv = REAL(steady_fv_wp, SiKi)
               steady_fv = steady_fv + grid_weight*grid_fv
            END DO
         END IF
         IF (.NOT. inst_e%agg%host_wave_enabled .AND. inst_e%agg%host_current_enabled) THEN
            ! The second sample contains waves plus the built-in steady current;
            ! subtract both, then restore the analytically reconstructed steady
            ! current. The remainder from the first sample is dynamic current.
            CALL WaveField_GetNodeWaveKin( p%WaveField, m%WaveField_m, t_sample, &
                                           REAL(inst_e%fx_xyz(:, i), ReKi), .TRUE., .FALSE., nodeInWater, &
                                           wave_zeta1, wave_zeta2, wave_zeta, wave_pdyn, wave_fv, wave_fa, &
                                           wave_famcf, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) RETURN
            fv = fv - wave_fv + steady_fv
            fa = fa - wave_fa
            zeta = 0.0_SiKi
         ELSE IF (.NOT. inst_e%agg%host_current_enabled) THEN
            ! fetchDynCurrent does not gate SeaState's built-in steady current.
            fv = fv - steady_fv
         END IF
         IF (ALLOCATED(inst_e%agg%host_current_profile_z)) THEN
            CALL CD_Current_Profile_Velocity(REAL(inst_e%fx_xyz(3, i), wp), &
                                             inst_e%agg%host_current_profile_z, &
                                             inst_e%agg%host_current_profile_velocity, &
                                             profile_fv, es, em)
            IF (es /= CD_HYDRO_OK) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: host-wave/file-current profile: '// &
                                TRIM(em), ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
            fv = fv + REAL(profile_fv, SiKi)
         END IF
         inst_e%fx_vel(:, i) = REAL(fv, wp)
         inst_e%fx_acc(:, i) = REAL(fa, wp)
         inst_e%fx_wl(i) = REAL(zeta, wp)
         ! the wave dynamic pressure (Pa) for the rod end caps: waves only, zero when the host
         ! wave component is disabled
         inst_e%fx_pd(i) = 0.0_wp
         IF (inst_e%agg%host_wave_enabled) inst_e%fx_pd(i) = REAL(pdyn, wp)
      END DO
      CALL CD_AGG_SetFluidFields( inst_e%agg, inst_e%fx_vel, inst_e%fx_acc, inst_e%fx_wl, es, em, &
                                  dynamic_pressure=inst_e%fx_pd )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
      END IF
   END SUBROUTINE sample_wavefield

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CD_CalcOutput( t, u, p, x, xd, z, other, y, m, ErrStat, ErrMsg )
      REAL(DbKi)                     , INTENT(IN   ) :: t
      TYPE( CD_InputType )           , INTENT(IN   ) :: u
      TYPE( CD_ParameterType )       , INTENT(IN   ) :: p
      TYPE( CD_ContinuousStateType ) , INTENT(IN   ) :: x
      TYPE( CD_DiscreteStateType )   , INTENT(IN   ) :: xd
      TYPE( CD_ConstraintStateType ) , INTENT(IN   ) :: z
      TYPE( CD_OtherStateType )      , INTENT(IN   ) :: other
      TYPE( CD_OutputType )          , INTENT(INOUT) :: y
      TYPE(CD_MiscVarType)           , INTENT(INOUT) :: m
      INTEGER(IntKi)                 , INTENT(INOUT) :: ErrStat
      CHARACTER(*)                   , INTENT(INOUT) :: ErrMsg

      character(*), parameter   :: RoutineName = 'CD_CalcOutput'
      INTEGER(IntKi)            :: ErrStat2, es, id, i
      INTEGER(IntKi)            :: prof_c0, prof_c1, prof_rate, prof_phase0
      CHARACTER(ErrMsgLen)      :: ErrMsg2
      CHARACTER(512)            :: em
      REAL(wp)                  :: cval
      LOGICAL                   :: fresh, state_overlaid

      ! Report this call's OWN status: clear at entry (ErrStat is INTENT(INOUT) and
      ! every SetErrStat below only ACCUMULATES, so on the all-success path an
      ! uninitialised incoming status would survive -- the FD-Jacobian probe passes
      ! its own uninitialised locals here). Mirrors every other CD_* entry point.
      ErrStat = ErrID_None
      ErrMsg  = ''

      ! CHECKPOINT-RESTART: the glue's first post-restore call can be CalcOutput, so this
      ! entry also rebuilds a non-matching instance (OtherState is INTENT(IN) here; the
      ! rebuild ADOPTS its identity without writing it).
      IF (.NOT. instance_matches(m, other)) THEN
         CALL rebuild_from_checkpoint( t, u, p, x, other, m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         ! (active line control is re-applied INSIDE rebuild_from_checkpoint, so
         ! this entry -- like every rebuild entry -- resumes at the checkpointed
         ! lengths; between advances CalcOutput stays state-const with the
         ! committed lengths, the supercycled zero-order-hold semantics)
      END IF
      id = m%CDInst
      IF (Inst(id)%prof_enabled) THEN
         CALL SYSTEM_CLOCK(prof_c0, prof_rate)
         Inst(id)%prof_n_calcout = Inst(id)%prof_n_calcout + 1_IntKi
      END IF
      state_overlaid = .FALSE.

      ! Evaluate the loads the mooring exerts on the coupling points. At nsuper = 1 (the
      ! mooring stepping at the glue rate) the CURRENT coupled kinematics are transferred
      ! first, so the loads carry the instantaneous direct feedthrough -- the original
      ! behavior, bit-for-bit. With SUPERCYCLING (nsuper > 1) the loads are instead the
      ! LAST COMMITTED mooring step's (a zero-order hold): transferring the instantaneous
      ! endpoint against the lagged interior would ring the stiff EA elastic feedthrough
      ! with a per-glue-step tension sawtooth, whereas the held loads lag by at most
      ! dt_moor. All scratch is the instance's persistent workspace.
      IF (Inst(id)%nsuper == 1_IntKi) THEN
         ! Gather the caller input before touching the aggregate.
         IF (Inst(id)%nT > 0) THEN
            CALL farm_gather_kinematics( Inst(id), u, ErrStat2, ErrMsg2 )
         ELSE
            CALL mesh_to_arrays( u%CoupledKinematics(1), Inst(id)%ncp, Inst(id)%ws_pos, Inst(id)%ws_vel, &
                                 Inst(id)%ws_acc, Inst(id)%ws_orient, Inst(id)%ws_omega, Inst(id)%ws_alpha, &
                                 ErrStat2, ErrMsg2 )
         END IF
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         fresh = .TRUE.
         ! At the glue rate, every call re-evaluates with the current kinematics
         ! transferred -- and then RESTORED. CalcOutput must be state-const: the glue
         ! calls it out-of-band with off-step inputs (checkpoint packing, extra output
         ! evaluations, probes), and committing those boundary kinematics would move the
         ! module's committed state. Save the committed coupled kinematics, transfer,
         ! evaluate, restore: the direct-feedthrough loads keep their semantics and the
         ! committed state is untouched on exit.
         CALL CD_AGG_GetMovingPointMesh( Inst(id)%agg, Inst(id)%sv_pos, Inst(id)%sv_vel, Inst(id)%sv_acc, &
                                         Inst(id)%ws_fld, es, em, orientation=Inst(id)%sv_orient, &
                                         angular_velocity=Inst(id)%sv_omega, angular_acceleration=Inst(id)%sv_alpha )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         ! ... and the FULL committed state. The output-probe overlay below changes only
         ! prescribed endpoints; the mirror remains the final authority for every
         ! committed generalized-alpha interior value.
         IF (.NOT. ALLOCATED(Inst(id)%ws_msave)) THEN
            ALLOCATE (Inst(id)%ws_msave(m%Nx), STAT=ErrStat2)
            IF (ErrStat2 /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, 'Could not allocate the CalcOutput state mirror.', &
                                ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
         END IF
         CALL pack_state_mirror_ary( Inst(id), m, Inst(id)%ws_msave, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         state_overlaid = .TRUE. ! cleanup is required even if an output probe fails part-way
         CALL CD_AGG_UpdateStates_Moving( Inst(id)%agg, Inst(id)%ws_pos, Inst(id)%ws_vel, Inst(id)%ws_acc, es, em, &
                                          output_probe=.TRUE., orientation=Inst(id)%ws_orient, &
                                          angular_velocity=Inst(id)%ws_omega, angular_acceleration=Inst(id)%ws_alpha )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
            CALL restore_calcoutput_probe( Inst(id), m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
      ELSE
         ! Supercycled: the committed state changes only when a mooring step completes, so the
         ! held loads/channels are re-evaluated only then (serving the cache is bit-identical
         ! to re-evaluating the unchanged state, at 1/nsuper of the evaluations).
         fresh = (.NOT. Inst(id)%loads_cached) .OR. (Inst(id)%loads_at_n /= Inst(id)%last_solve_n)
      END IF

      IF (fresh) THEN
         IF (Inst(id)%prof_enabled) THEN
            Inst(id)%prof_n_fresh = Inst(id)%prof_n_fresh + 1_IntKi
            CALL SYSTEM_CLOCK(prof_phase0)
         END IF
         CALL CD_AGG_CalcOutput( Inst(id)%agg, es, em )
         IF (es /= CD_AGG_OK) THEN
            IF (Inst(id)%prof_enabled) THEN
               CALL SYSTEM_CLOCK(prof_c1)
               Inst(id)%prof_t_loads = Inst(id)%prof_t_loads + &
                  REAL(prof_c1 - prof_phase0, wp)/REAL(prof_rate, wp)
            END IF
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
            IF (state_overlaid) THEN
               CALL restore_calcoutput_probe( Inst(id), m, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            END IF
            RETURN
         END IF
         CALL CD_AGG_GetMovingPointMesh( Inst(id)%agg, Inst(id)%ws_pos, Inst(id)%ws_vel, Inst(id)%ws_acc, &
                                         Inst(id)%ws_fld, es, em, orientation=Inst(id)%ws_orient, &
                                         moment=Inst(id)%ws_moment, angular_velocity=Inst(id)%ws_omega, &
                                         angular_acceleration=Inst(id)%ws_alpha )
         ! The aggregate CalcOutput above recovers EI=0 loads; finite-EI fairlead reactions are
         ! evaluated on demand by GetMovingPointMesh. Stop the phase timer only after both halves
         ! so mixed-run load recovery is not misattributed to the enclosing CalcOutput time.
         IF (Inst(id)%prof_enabled) THEN
            CALL SYSTEM_CLOCK(prof_c1)
            Inst(id)%prof_t_loads = Inst(id)%prof_t_loads + &
               REAL(prof_c1 - prof_phase0, wp)/REAL(prof_rate, wp)
         END IF
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
            IF (state_overlaid) THEN
               CALL restore_calcoutput_probe( Inst(id), m, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            END IF
            RETURN
         END IF
         ! Deck OUTPUTS channels are evaluated at the same aggregate state as the returned
         ! loads: with supercycling the last committed zero-order-hold mooring state, at
         ! nsuper = 1 the state including the current coupled kinematics.
         DO i = 1, p%NumOuts
            CALL CD_AGG_EvalChannel( Inst(id)%agg, i, cval, es, em )
            IF (es /= CD_AGG_OK) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: output channel evaluation failed: '// &
                                TRIM(em), ErrStat, ErrMsg, RoutineName )
               IF (state_overlaid) THEN
                  CALL restore_calcoutput_probe( Inst(id), m, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF
               RETURN
            END IF
            Inst(id)%ws_chan(i) = cval
         END DO
         Inst(id)%loads_cached = .TRUE.
         Inst(id)%loads_at_n = Inst(id)%last_solve_n
      END IF

      ! restore the committed coupled kinematics saved above (state-const CalcOutput),
      ! then reload the interiors from the captured mirror: the boundary transfer just
      ! re-derived interior accelerations, and only the mirror carries the committed
      ! gen-alpha values.
      IF (state_overlaid) THEN
         CALL restore_calcoutput_probe( Inst(id), m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF

      ! CableDyn-owned channel recording: one row per COMMITTED mooring advance, stamped with
      ! the commit time and evaluated HERE -- after the committed-state restore -- because
      ! at nsuper = 1 the fresh branch above intentionally evaluates ws_chan at the
      ! caller's CURRENT input (the documented instantaneous feedthrough an out-of-band
      ! probe/extra CalcOutput may drive with transient kinematics). The row must carry
      ! the committed state the stamp names, so its channels re-evaluate against the
      ! restored aggregate (at nsuper > 1 the committed state is live anyway and the
      ! re-evaluation is identical). Keyed on last_solve_n: exactly one row per advance,
      ! whichever call after it comes first; the t = 0 static row was written at Init.
      ! If a glue correction re-advances the same n, BACKSPACE replaces the provisional
      ! last record with the corrected committed result rather than duplicating the time.
      IF (Inst(id)%native_un > 0 .AND. Inst(id)%native_pending .AND. Inst(id)%last_solve_n >= 0_IntKi) THEN
         IF (Inst(id)%prof_enabled) CALL SYSTEM_CLOCK(prof_phase0)
         IF (Inst(id)%last_solve_n == Inst(id)%native_wrote_n) BACKSPACE (Inst(id)%native_un)
         WRITE (Inst(id)%native_un, '(ES25.16E3)', ADVANCE='NO') MAX(Inst(id)%t_commit, 0.0_DbKi)
         DO i = 1, p%NumOuts
            CALL CD_AGG_EvalChannel( Inst(id)%agg, i, cval, es, em )
            IF (es /= CD_AGG_OK) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: committed-step channel record: '//TRIM(em), &
                                ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            WRITE (Inst(id)%native_un, '(A,ES15.6E2)', ADVANCE='NO') CHAR(9), cval
         END DO
         WRITE (Inst(id)%native_un, '(A)') ''
         ! keep the record readable during the run (and intact after a crash)
         FLUSH (Inst(id)%native_un)
         Inst(id)%native_wrote_n = Inst(id)%last_solve_n
         Inst(id)%native_wrote_t = MAX(Inst(id)%t_commit, 0.0_DbKi)
         Inst(id)%native_pending = .FALSE.
         IF (Inst(id)%prof_enabled) THEN
            CALL SYSTEM_CLOCK(prof_c1)
            Inst(id)%prof_t_native = Inst(id)%prof_t_native + &
               REAL(prof_c1 - prof_phase0, wp)/REAL(prof_rate, wp)
         END IF
      END IF

      IF (Inst(id)%nT > 0) THEN
         DO i = 1, Inst(id)%nT
            y%CoupledLoads(i)%Force  = 0.0_ReKi
            y%CoupledLoads(i)%Moment = 0.0_ReKi
         END DO
         DO i = 1, Inst(id)%ncp
            y%CoupledLoads(Inst(id)%tmap(i))%Force(:, Inst(id)%node_of_slot(i)) = &
               REAL(Inst(id)%ws_fld(:, i), ReKi)
            y%CoupledLoads(Inst(id)%tmap(i))%Moment(:, Inst(id)%node_of_slot(i)) = &
               REAL(Inst(id)%ws_moment(:, i), ReKi)
         END DO
      ELSE
      DO i = 1, Inst(id)%ncp
         y%CoupledLoads(1)%Force(:, i) = REAL(Inst(id)%ws_fld(:, i), ReKi)
         y%CoupledLoads(1)%Moment(:, i) = REAL(Inst(id)%ws_moment(:, i), ReKi)
      END DO
      END IF
      DO i = 1, p%NumOuts
         y%WriteOutput(i) = REAL(Inst(id)%ws_chan(i), ReKi)
      END DO

      IF (Inst(id)%prof_enabled) THEN
         CALL SYSTEM_CLOCK(prof_c1)
         Inst(id)%prof_t_calcout = Inst(id)%prof_t_calcout + &
            REAL(prof_c1 - prof_c0, wp)/REAL(prof_rate, wp)
      END IF

   END SUBROUTINE CD_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE restore_calcoutput_probe( inst_e, m, ErrStat, ErrMsg )
      !! Restore every state surface after the output-probe
      !! overlay. This is used on both success and every post-overlay error exit.
      TYPE(CD_OF_Instance), INTENT(INOUT) :: inst_e
      TYPE(CD_MiscVarType), INTENT(INOUT) :: m
      INTEGER(IntKi), INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg
      INTEGER(IntKi) :: es, ErrStat2
      CHARACTER(512) :: em
      CHARACTER(ErrMsgLen) :: ErrMsg2

      ErrStat = ErrID_None
      ErrMsg = ''
      CALL CD_AGG_UpdateStates_Moving( inst_e%agg, inst_e%sv_pos, inst_e%sv_vel, inst_e%sv_acc, es, em, &
                                       output_probe=.TRUE., orientation=inst_e%sv_orient, &
                                       angular_velocity=inst_e%sv_omega, angular_acceleration=inst_e%sv_alpha )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn CalcOutput endpoint restore failed: '//TRIM(em), &
                          ErrStat, ErrMsg, 'restore_calcoutput_probe' )
      END IF
      CALL reload_interiors_ary( inst_e, m, inst_e%ws_msave, &
                                 REAL(MAX(inst_e%t_commit, 0.0_DbKi), wp), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'restore_calcoutput_probe' )
   END SUBROUTINE restore_calcoutput_probe

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CD_CalcContStateDeriv( t, u, p, x, xd, z, other, m, dxdt, ErrStat, ErrMsg )
      REAL(DbKi),                         INTENT(IN )    :: t
      TYPE(CD_InputType),                 INTENT(IN )    :: u
      TYPE(CD_ParameterType),             INTENT(IN )    :: p
      TYPE(CD_ContinuousStateType),       INTENT(IN )    :: x
      TYPE(CD_DiscreteStateType),         INTENT(IN )    :: xd
      TYPE(CD_ConstraintStateType),       INTENT(IN )    :: z
      TYPE(CD_OtherStateType),            INTENT(IN )    :: other
      TYPE(CD_MiscVarType),               INTENT(INOUT)  :: m
      TYPE(CD_ContinuousStateType),       INTENT(INOUT)  :: dxdt
      INTEGER(IntKi),                     INTENT( OUT)   :: ErrStat
      CHARACTER(*),                       INTENT( OUT)   :: ErrMsg

      character(*), parameter :: RoutineName = 'CD_CalcContStateDeriv'

      ErrStat = ErrID_None
      ErrMsg  = ''
      ! Only the linearization/GetOP paths reach this routine (the module marches
      ! implicitly through UpdateStates). CableDyn exposes no continuous states.
      CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: continuous-state derivatives are part of '// &
                       'the linearization surface; CableDyn exposes no continuous states to OpenFAST.', &
                       ErrStat, ErrMsg, RoutineName )
   END SUBROUTINE CD_CalcContStateDeriv

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CD_End(u, p, x, xd, z, other, y, m, ErrStat , ErrMsg)
      TYPE(CD_InputType) ,            INTENT(INOUT) :: u
      TYPE(CD_ParameterType) ,        INTENT(INOUT) :: p
      TYPE(CD_ContinuousStateType) ,  INTENT(INOUT) :: x
      TYPE(CD_DiscreteStateType) ,    INTENT(INOUT) :: xd
      TYPE(CD_ConstraintStateType) ,  INTENT(INOUT) :: z
      TYPE(CD_OtherStateType) ,       INTENT(INOUT) :: other
      TYPE(CD_OutputType) ,           INTENT(INOUT) :: y
      TYPE(CD_MiscVarType),           INTENT(INOUT) :: m
      INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat
      CHARACTER(*),                   INTENT(  OUT) :: ErrMsg

      character(*), parameter :: RoutineName = 'CD_End'
      INTEGER(IntKi)          :: ErrStat2
      INTEGER                 :: hn_step, hn_resid, hn_solve, hn_am, hn_alloc, hn_tan
      REAL(wp)                :: ht_step, ht_resid, ht_eff, ht_solve, ht_am, ht_tan
      CHARACTER(ErrMsgLen)    :: ErrMsg2

      ErrStat = ErrID_None
      ErrMsg  = ''

      IF (ALLOCATED(Inst)) THEN
         IF (m%CDInst >= 1 .AND. m%CDInst <= SIZE(Inst)) THEN
            IF (Inst(m%CDInst)%active) THEN
               IF (Inst(m%CDInst)%prof_enabled) THEN
                  CALL WrScr( '    CableDyn profile: init_s='// &
                              TRIM(Num2LStr(REAL(Inst(m%CDInst)%prof_t_init, DbKi)))// &
                              ' update_s='//TRIM(Num2LStr(REAL(Inst(m%CDInst)%prof_t_update, DbKi)))// &
                              ' wave_s='//TRIM(Num2LStr(REAL(Inst(m%CDInst)%prof_t_wave, DbKi)))// &
                              ' step_s='//TRIM(Num2LStr(REAL(Inst(m%CDInst)%prof_t_step, DbKi))) )
                  CALL WrScr( '    CableDyn profile: calcout_s='// &
                              TRIM(Num2LStr(REAL(Inst(m%CDInst)%prof_t_calcout, DbKi)))// &
                              ' loads_s='//TRIM(Num2LStr(REAL(Inst(m%CDInst)%prof_t_loads, DbKi)))// &
                              ' native_io_s='//TRIM(Num2LStr(REAL(Inst(m%CDInst)%prof_t_native, DbKi))) )
                  CALL WrScr( '    CableDyn profile counts: update='// &
                              TRIM(Num2LStr(Inst(m%CDInst)%prof_n_update))// &
                              ' solve='//TRIM(Num2LStr(Inst(m%CDInst)%prof_n_solve))// &
                              ' wave='//TRIM(Num2LStr(Inst(m%CDInst)%prof_n_wave))// &
                              ' calcout='//TRIM(Num2LStr(Inst(m%CDInst)%prof_n_calcout))// &
                              ' fresh='//TRIM(Num2LStr(Inst(m%CDInst)%prof_n_fresh)) )
                  IF (Inst(m%CDInst)%agg%ncable > 0) THEN
                     CALL CD_HermiteCable_Dyn_Get_Profile(hn_step, hn_resid, hn_solve, hn_am, hn_alloc, &
                                                          ht_step, ht_resid, ht_eff, ht_solve, ht_am, &
                                                          hn_tan, ht_tan)
                     CALL WrScr( '    CableDyn Hermite profile: step_s='//TRIM(Num2LStr(REAL(ht_step, DbKi)))// &
                                 ' residual_s='//TRIM(Num2LStr(REAL(ht_resid, DbKi)))// &
                                 ' tangent_s='//TRIM(Num2LStr(REAL(ht_tan, DbKi)))// &
                                 ' effective_s='//TRIM(Num2LStr(REAL(ht_eff, DbKi)))// &
                                 ' solve_s='//TRIM(Num2LStr(REAL(ht_solve, DbKi)))// &
                                 ' added_mass_s='//TRIM(Num2LStr(REAL(ht_am, DbKi))) )
                     CALL WrScr( '    CableDyn Hermite counts: step='//TRIM(Num2LStr(hn_step))// &
                                 ' residual='//TRIM(Num2LStr(hn_resid))// &
                                 ' tangent='//TRIM(Num2LStr(hn_tan))// &
                                 ' solve='//TRIM(Num2LStr(hn_solve))// &
                                 ' added_mass='//TRIM(Num2LStr(hn_am))// &
                                 ' alloc='//TRIM(Num2LStr(hn_alloc)) )
                     CALL WrScr( '    CableDyn Hermite recoveries: intervals='// &
                                 TRIM(Num2LStr(CD_HermiteCable_Dyn_Recovery_Count())) )
                  END IF
               END IF
               ! Propagate a teardown failure as SEVERE (not fatal): the glue's End path
               ! aggregates statuses and continues cleanup; swallowing it would hide it.
               CALL release_instance( m%CDInst, ErrStat2, ErrMsg2 )
               IF (ErrStat2 /= ErrID_None) CALL SetErrStat( ErrID_Severe, TRIM(ErrMsg2), ErrStat, ErrMsg, RoutineName )
            END IF
         END IF
      END IF

      ! Framework cleanup of everything the registry types allocated.
      CALL CD_DestroyInput( u, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CD_DestroyParam( p, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CD_DestroyContState( x, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CD_DestroyDiscState( xd, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CD_DestroyConstrState( z, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CD_DestroyOtherState( other, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CD_DestroyOutput( y, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL CD_DestroyMisc( m, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END SUBROUTINE CD_End

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CD_JacobianPInput(Vars, t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu)
      TYPE(ModVarsType),                  INTENT(IN   ) :: Vars
      REAL(DbKi),                         INTENT(IN   ) :: t
      TYPE(CD_InputType),                 INTENT(INOUT) :: u
      TYPE(CD_ParameterType),             INTENT(IN   ) :: p
      TYPE(CD_ContinuousStateType),       INTENT(IN   ) :: x
      TYPE(CD_DiscreteStateType),         INTENT(IN   ) :: xd
      TYPE(CD_ConstraintStateType),       INTENT(IN   ) :: z
      TYPE(CD_OtherStateType),            INTENT(IN   ) :: OtherState
      TYPE(CD_OutputType),                INTENT(INOUT) :: y
      TYPE(CD_MiscVarType),               INTENT(INOUT) :: m
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dYdu(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dXdu(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dXddu(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dZdu(:,:)

      character(*), parameter :: RoutineName = 'CD_JacobianPInput'
      INTEGER(IntKi)          :: ErrStat2, i, j, iCol
      CHARACTER(ErrMsgLen)    :: ErrMsg2

      ErrStat = ErrID_None
      ErrMsg  = ''
      ! The Failed() idiom merges these locals into ErrStat; seed them clean so a
      ! probe never reads a stale/uninitialised status (belt-and-suspenders with the
      ! CD_CalcOutput entry-clear).
      ErrStat2 = ErrID_None
      ErrMsg2  = ''

      ! Get OP values here (the FD reference point)
      CALL CD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 ); IF (Failed()) RETURN

      ! dY/du has TWO distinct contracts. During a normal nonlinear OpenFAST run the
      ! registered CableDyn states exist (Nx > 0), and the glue input/output solve needs
      ! the PARTIAL derivative g_u at fixed internal state. Re-equilibrating the line for
      ! those probes returns a total quasi-static derivative instead; inserting that into
      ! the runtime algebraic loop can make its Jacobian rank deficient. For formal
      ! linearization CableDyn deliberately registers no states (Nx = 0), so the useful
      ! reduction is the re-equilibrated, zero-frequency stiffness. Keep those paths
      ! explicit rather than using the linearization reduction in the runtime solve.
      IF (PRESENT(dYdu)) THEN
         IF (.NOT. ALLOCATED(dYdu)) THEN
            CALL AllocAry( dYdu, m%Jac%Ny, m%Jac%Nu, 'dYdu', ErrStat2, ErrMsg2 ); IF (Failed()) RETURN
         END IF
         dYdu = 0.0_R8Ki
         IF (m%Jac%Nx > 0) THEN
            ! Match stock MoorDyn's runtime Jacobian contract: perturb the current input,
            ! evaluate CalcOutput with the same x/xd/z/other, and central-difference the
            ! packed outputs. CD_CalcOutput's probe overlay restores the authoritative
            ! CableDyn state after every evaluation.
            CALL CD_CopyInput( u, m%u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); IF (Failed()) RETURN
            CALL CD_VarsPackInput( Vars, u, m%Jac%u )
            DO i = 1, SIZE(Vars%u)
               DO j = 1, Vars%u(i)%Num
                  iCol = Vars%u(i)%iLoc(1) + j - 1
                  CALL MV_Perturb( Vars%u(i), j, 1, m%Jac%u, m%Jac%u_perturb )
                  CALL CD_VarsUnpackInput( Vars, m%Jac%u_perturb, m%u_perturb )
                  CALL CD_CalcOutput( t, m%u_perturb, p, x, xd, z, OtherState, m%y_lin, m, ErrStat2, ErrMsg2 )
                  IF (Failed()) RETURN
                  CALL CD_VarsPackOutput( Vars, m%y_lin, m%Jac%y_pos )

                  CALL MV_Perturb( Vars%u(i), j, -1, m%Jac%u, m%Jac%u_perturb )
                  CALL CD_VarsUnpackInput( Vars, m%Jac%u_perturb, m%u_perturb )
                  CALL CD_CalcOutput( t, m%u_perturb, p, x, xd, z, OtherState, m%y_lin, m, ErrStat2, ErrMsg2 )
                  IF (Failed()) RETURN
                  CALL CD_VarsPackOutput( Vars, m%y_lin, m%Jac%y_neg )

                  CALL MV_ComputeCentralDiff( Vars%y, Vars%u(i)%Perturb, m%Jac%y_pos, m%Jac%y_neg, &
                                              dYdu(:, iCol) )
               END DO
            END DO
            ! The final negative probe leaves only output/cache workspaces at the
            ! perturbed input (the authoritative structural state was already restored).
            ! Re-evaluate the operating point so the callback is observationally neutral
            ! to the next glue evaluation as well as state neutral.
            CALL CD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
            IF (Failed()) RETURN
         ELSE
            ! Zero-state formal linearization: each displacement probe re-solves static
            ! equilibrium on a scratch aggregate, giving the reduced mooring stiffness.
            CALL quasistatic_dydu( Vars, u, p, m, dYdu, ErrStat2, ErrMsg2 ); IF (Failed()) RETURN
         END IF
      END IF

      ! dX/du: the implicit march has no x' = f(x, u) surface. Under linearization the
      ! module registers ZERO state variables (see CD_InitVars), so the honest dXdu is
      ! the empty (0 x Nu) block; a caller holding a state-registered instance (Nx > 0)
      ! is asking for a surface that does not exist -- fail loudly, never zero-fill.
      IF (PRESENT(dXdu)) THEN
         IF (m%Jac%Nx > 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: dXdu is undefined for the implicit '// &
                             'march with registered state variables; linearization must be requested at '// &
                             'Init so the module presents its zero-state (dYdu-only) reduction.', &
                             ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         IF (.NOT. ALLOCATED(dXdu)) THEN
            CALL AllocAry( dXdu, m%Jac%Nx, m%Jac%Nu, 'dXdu', ErrStat2, ErrMsg2 ); IF (Failed()) RETURN
         END IF
      END IF

      IF (PRESENT(dXddu)) THEN
         IF (ALLOCATED(dXddu)) DEALLOCATE (dXddu)
      END IF
      IF (PRESENT(dZdu)) THEN
         IF (ALLOCATED(dZdu)) DEALLOCATE (dZdu)
      END IF

   CONTAINS
      LOGICAL FUNCTION Failed()
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         Failed = ErrStat >= AbortErrLev
      END FUNCTION Failed
   END SUBROUTINE CD_JacobianPInput

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE quasistatic_dydu( Vars, u, p, m, dYdu, ErrStat, ErrMsg )
      !! The quasi-static dYdu displacement columns: for each coupled-mesh
      !! FieldTransDisp perturbation, re-solve the mooring static equilibrium at the
      !! perturbed boundary on a SCRATCH aggregate built from the deck copy (each
      !! CD_AGG_Init_From_Deck IS a static solve; coupled_positions places the moving
      !! points explicitly), read the loads and deck channels, and central-difference.
      !! The live instance is never touched -- the probe cannot perturb committed
      !! state by construction. All other input fields keep zero columns.
      TYPE(ModVarsType),        INTENT(IN   ) :: Vars
      TYPE(CD_InputType),       INTENT(IN   ) :: u
      TYPE(CD_ParameterType),   INTENT(IN   ) :: p
      TYPE(CD_MiscVarType),     INTENT(INOUT) :: m
      REAL(R8Ki),               INTENT(INOUT) :: dYdu(:,:)
      INTEGER(IntKi),           INTENT(  OUT) :: ErrStat
      CHARACTER(*),             INTENT(  OUT) :: ErrMsg

      character(*), parameter   :: RoutineName = 'quasistatic_dydu'
      INTEGER(IntKi)            :: ErrStat2, es, i, j, iCol, node, comp, ncp, k
      CHARACTER(ErrMsgLen)      :: ErrMsg2
      CHARACTER(2100)           :: lin_deck
      REAL(wp), ALLOCATABLE     :: op_pos(:,:), prb_pos(:,:)
      LOGICAL                   :: use_extf

      ErrStat = ErrID_None
      ErrMsg  = ''
      IF (p%DeckCopy%NumLines < 1) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the linearization probe needs the '// &
                          'deck copy in the parameters.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      ncp = u%CoupledKinematics(1)%Nnodes
      ALLOCATE (op_pos(3, ncp), prb_pos(3, ncp), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate the linearization probe positions.', &
                          ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      ! the operating-point boundary: the committed mesh positions at the snapshot
      DO i = 1, ncp
         op_pos(:, i) = REAL(u%CoupledKinematics(1)%Position(:, i), wp) + &
                        REAL(u%CoupledKinematics(1)%TranslationDisp(:, i), wp)
      END DO

      ! the deck copy back to disk once (beside the MooringFile, so its side files
      ! resolve); every probe init reuses it and it is removed afterwards
      CALL write_deck_copy( p, '.lin.dat', lin_deck, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      use_extf = .FALSE.
      IF (ASSOCIATED(p%WaveField)) use_extf = wavefield_is_ambient( p%WaveField )

      DO i = 1, SIZE(Vars%u)
         IF (Vars%u(i)%Field /= FieldTransDisp) CYCLE   ! vel/acc/orientation: zero columns
         DO j = 1, Vars%u(i)%Num
            iCol = Vars%u(i)%iLoc(1) + j - 1
            node = (j - 1)/3 + 1
            comp = MOD(j - 1, 3) + 1
            IF (node > ncp) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: linearization probe node index '// &
                                'exceeds the coupled mesh.', ErrStat, ErrMsg, RoutineName )
               CALL delete_file( lin_deck )
               RETURN
            END IF
            DO k = 1, 2   ! +perturb then -perturb
               prb_pos = op_pos
               IF (k == 1) THEN
                  prb_pos(comp, node) = prb_pos(comp, node) + REAL(Vars%u(i)%Perturb, wp)
               ELSE
                  prb_pos(comp, node) = prb_pos(comp, node) - REAL(Vars%u(i)%Perturb, wp)
               END IF
               CALL probe_loads( prb_pos, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL delete_file( lin_deck )
                  RETURN
               END IF
               IF (k == 1) THEN
                  CALL CD_VarsPackOutput( Vars, m%y_lin, m%Jac%y_pos )
               ELSE
                  CALL CD_VarsPackOutput( Vars, m%y_lin, m%Jac%y_neg )
               END IF
            END DO
            CALL MV_ComputeCentralDiff( Vars%y, Vars%u(i)%Perturb, m%Jac%y_pos, m%Jac%y_neg, dYdu(:, iCol) )
         END DO
      END DO
      CALL delete_file( lin_deck )

   CONTAINS

      SUBROUTINE probe_loads( pos, PErrStat, PErrMsg )
         !! One static probe: scratch-aggregate init at the given boundary (the init
         !! IS the static solve), loads + deck channels into m%y_lin, then release.
         REAL(wp),       INTENT(IN   ) :: pos(:,:)
         INTEGER(IntKi), INTENT(  OUT) :: PErrStat
         CHARACTER(*),   INTENT(  OUT) :: PErrMsg

         TYPE(CD_AGG_ModuleType) :: agg_p
         REAL(wp)                :: q_p(3, ncp), v_p(3, ncp), a_p(3, ncp), f_p(3, ncp)
         REAL(wp)                :: cval
         INTEGER(IntKi)          :: es_p, ich
         CHARACTER(512)          :: em_p

         PErrStat = ErrID_None
         PErrMsg  = ''
         CALL CD_AGG_Init_From_Deck( agg_p, TRIM(lin_deck), REAL(p%dtM0, wp), es_p, em_p, &
                                     env_gravity=REAL(p%g, wp), env_rho_water=REAL(p%rhoW, wp), &
                                     env_wtrdpth=REAL(p%WtrDpth, wp), external_fluid=use_extf, &
                                     coupled_positions=pos, run_tmax=REAL(p%Tmax, wp) )
         IF (es_p /= CD_AGG_OK) THEN
            PErrStat = ErrID_Fatal
            PErrMsg  = 'linearization static probe init failed: '//TRIM(em_p)
            RETURN
         END IF
         CALL CD_AGG_CalcOutput( agg_p, es_p, em_p )
         IF (es_p == CD_AGG_OK) CALL CD_AGG_GetMovingPointMesh( agg_p, q_p, v_p, a_p, f_p, es_p, em_p )
         IF (es_p /= CD_AGG_OK) THEN
            CALL CD_AGG_End( agg_p, es_p, em_p )
            PErrStat = ErrID_Fatal
            PErrMsg  = 'linearization static probe output failed: '//TRIM(em_p)
            RETURN
         END IF
         DO ich = 1, ncp
            m%y_lin%CoupledLoads(1)%Force(:, ich) = REAL(f_p(:, ich), ReKi)
            m%y_lin%CoupledLoads(1)%Moment(:, ich) = 0.0_ReKi
         END DO
         DO ich = 1, p%NumOuts
            CALL CD_AGG_EvalChannel( agg_p, ich, cval, es_p, em_p )
            IF (es_p /= CD_AGG_OK) THEN
               CALL CD_AGG_End( agg_p, es_p, em_p )
               PErrStat = ErrID_Fatal
               PErrMsg  = 'linearization static probe channel failed: '//TRIM(em_p)
               RETURN
            END IF
            m%y_lin%WriteOutput(ich) = REAL(cval, ReKi)
         END DO
         CALL CD_AGG_End( agg_p, es_p, em_p )
      END SUBROUTINE probe_loads

   END SUBROUTINE quasistatic_dydu

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CD_JacobianPContState(Vars, t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx)
      TYPE(ModVarsType),                  INTENT(IN   ) :: Vars
      REAL(DbKi),                         INTENT(IN   ) :: t
      TYPE(CD_InputType),                 INTENT(INOUT) :: u
      TYPE(CD_ParameterType),             INTENT(IN   ) :: p
      TYPE(CD_ContinuousStateType),       INTENT(INOUT) :: x
      TYPE(CD_DiscreteStateType),         INTENT(IN   ) :: xd
      TYPE(CD_ConstraintStateType),       INTENT(IN   ) :: z
      TYPE(CD_OtherStateType),            INTENT(IN   ) :: OtherState
      TYPE(CD_OutputType),                INTENT(INOUT) :: y
      TYPE(CD_MiscVarType),               INTENT(INOUT) :: m
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dYdx(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dXdx(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dXddx(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dZdx(:,:)

      character(*), parameter :: RoutineName = 'CD_JacobianPContState'

      INTEGER(IntKi)          :: ErrStat2
      CHARACTER(ErrMsgLen)    :: ErrMsg2

      ErrStat = ErrID_None
      ErrMsg  = ''
      ErrStat2 = ErrID_None
      ErrMsg2  = ''
      ! Under linearization the module registers ZERO state variables (the MAP-shaped
      ! dYdu-only reduction; see CD_InitVars), so every state Jacobian is an empty
      ! block. A state-registered instance (Nx > 0) has no x' = f(x, u) surface to
      ! differentiate -- fail loudly rather than fabricate one.
      IF (m%Jac%Nx > 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: state Jacobians are undefined for the '// &
                          'implicit march with registered state variables; linearization must be '// &
                          'requested at Init so the module presents its zero-state (dYdu-only) '// &
                          'reduction.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (PRESENT(dYdx)) THEN
         IF (.NOT. ALLOCATED(dYdx)) THEN
            CALL AllocAry( dYdx, m%Jac%Ny, m%Jac%Nx, 'dYdx', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) RETURN
         END IF
      END IF
      IF (PRESENT(dXdx)) THEN
         IF (.NOT. ALLOCATED(dXdx)) THEN
            CALL AllocAry( dXdx, m%Jac%Nx, m%Jac%Nx, 'dXdx', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) RETURN
         END IF
      END IF
      IF (PRESENT(dXddx)) THEN
         IF (ALLOCATED(dXddx)) DEALLOCATE (dXddx)
      END IF
      IF (PRESENT(dZdx)) THEN
         IF (ALLOCATED(dZdx)) DEALLOCATE (dZdx)
      END IF
   END SUBROUTINE CD_JacobianPContState

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
      REAL(DbKi),                         INTENT(IN   ) :: t
      TYPE(CD_InputType),                 INTENT(INOUT) :: u
      TYPE(CD_ParameterType),             INTENT(IN   ) :: p
      TYPE(CD_ContinuousStateType),       INTENT(IN   ) :: x
      TYPE(CD_DiscreteStateType),         INTENT(IN   ) :: xd
      TYPE(CD_ConstraintStateType),       INTENT(IN   ) :: z
      TYPE(CD_OtherStateType),            INTENT(IN   ) :: OtherState
      TYPE(CD_OutputType),                INTENT(INOUT) :: y
      TYPE(CD_MiscVarType),               INTENT(INOUT) :: m
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dYdxd(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dXdxd(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dXddxd(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dZdxd(:,:)

      ! No discrete states (matches stock: a no-op).
      ErrStat = ErrID_None
      ErrMsg  = ''
   END SUBROUTINE CD_JacobianPDiscState

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
      REAL(DbKi),                         INTENT(IN   ) :: t
      TYPE(CD_InputType),                 INTENT(INOUT) :: u
      TYPE(CD_ParameterType),             INTENT(IN   ) :: p
      TYPE(CD_ContinuousStateType),       INTENT(IN   ) :: x
      TYPE(CD_DiscreteStateType),         INTENT(IN   ) :: xd
      TYPE(CD_ConstraintStateType),       INTENT(IN   ) :: z
      TYPE(CD_OtherStateType),            INTENT(IN   ) :: OtherState
      TYPE(CD_OutputType),                INTENT(INOUT) :: y
      TYPE(CD_MiscVarType),               INTENT(INOUT) :: m
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dYdz(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dXdz(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dXddz(:,:)
      REAL(R8Ki), ALLOCATABLE, OPTIONAL,  INTENT(INOUT) :: dZdz(:,:)

      ! No constraint states (matches stock: a no-op).
      ErrStat = ErrID_None
      ErrMsg  = ''
   END SUBROUTINE CD_JacobianPConstrState

!----------------------------------------------------------------------------------------------------------------------------------
!  Private helpers
!----------------------------------------------------------------------------------------------------------------------------------

   !> Register the module variables (the ModVars surface MV_AddModule consumes):
   !> u = the coupled-kinematics mesh; y = the coupled-loads mesh; x = the line
   !> interior-node velocity/position mirror in the MoorDyn framework ordering
   !> (all velocities registered after all positions).
   SUBROUTINE CD_InitVars(Vars, u, p, x, y, m, inst_, Linearize, ErrStat, ErrMsg)
      TYPE(ModVarsType),               INTENT(OUT)    :: Vars
      TYPE(CD_InputType),              INTENT(INOUT)  :: u
      TYPE(CD_ParameterType),          INTENT(INOUT)  :: p
      TYPE(CD_ContinuousStateType),    INTENT(INOUT)  :: x
      TYPE(CD_OutputType),             INTENT(INOUT)  :: y
      TYPE(CD_MiscVarType),            INTENT(INOUT)  :: m
      TYPE(CD_OF_Instance),            INTENT(IN)     :: inst_
      LOGICAL,                         INTENT(IN)     :: Linearize
      INTEGER(IntKi),                  INTENT(OUT)    :: ErrStat
      CHARACTER(*),                    INTENT(OUT)    :: ErrMsg

      character(*), parameter   :: RoutineName = 'CD_InitVars'
      INTEGER(IntKi)            :: ErrStat2, es
      CHARACTER(ErrMsgLen)      :: ErrMsg2
      CHARACTER(512)            :: em
      INTEGER(IntKi)            :: il, i, j, N, ndofl, nin, ic, csz, nf, koff, nve
      CHARACTER(32)             :: LinStr
      CHARACTER(ChanLen)        :: whdr, wunt
      INTEGER(IntKi)            :: ivT
      CHARACTER(LinChanLen), ALLOCATABLE :: WrOutLinNames(:)
      CHARACTER(LinChanLen), ALLOCATABLE :: CblLinVel(:), CblLinPos(:), CblLinAcc(:)
      character(20), parameter  :: TransDispSuffix(*) = [' Px, m', ' Py, m', ' Pz, m']
      character(20), parameter  :: TransVelSuffix(*)  = [' Vx, m/s', ' Vy, m/s', ' Vz, m/s']

      ErrStat = ErrID_None
      ErrMsg  = ''

      IF (Linearize) THEN
         ! Linearization presents the module in the MAP shape: ZERO continuous-state
         ! variables. The module marches implicitly (generalised-alpha behind
         ! UpdateStates), so there is no x' = f(x, u) surface to expose, and the glue
         ! force-flags every registered x variable into the linearization set. With
         ! Nx = 0 the linearized mooring is the direct-feedthrough impedance dYdu --
         ! stiffness from the displacement columns, damping from the velocity columns,
         ! inertia from the acceleration columns: the quasi-static reduction MAP
         ! exposes, and what a FOWT linearization consumes from a mooring. The
         ! checkpoint and state-copy machinery is UNAFFECTED: pack_state_mirror and
         ! the registry state copies use x%states directly, never Vars%x.
         ALLOCATE (Vars%x(0), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate Vars%x.', ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
      ELSE

      ! position states (interior line nodes; the mirror layout stores [rd, r] per
      ! line, so positions start at LineStateIs1 + 3*(N-1)).
      ! NB the inline LinNames constructors below are SAFE: each element is TRIM(LinStr) (a
      ! fixed scalar) concatenated with a PARAMETER-array suffix element -- no function call
      ! inside the constructor. The heap-corruption class this file guards against (see the
      ! finite-EI cable block and the WriteOutput block, both built via explicit loops) is
      ! specifically an implied-do character array constructor whose element expression calls
      ! a function that performs internal I/O (Num2LStr does a WRITE); those must be pre-built.
      DO il = 1, inst_%nlines
         ndofl = CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         N = ndofl/3 - 1                       ! number of segments
         DO i = 0, N - 2
            LinStr = 'Line '//TRIM(Num2LStr(il))//' node '//TRIM(Num2LStr(i + 1))
            CALL MV_AddVar( Vars%x, LinStr, FieldTransDisp, DatLoc(CD_x_states), &
                            iAry=m%LineStateIs1(il) + 3*(N - 1) + 3*i, &
                            Num=3, Flags=VF_DerivOrder2, Perturb=0.05_R8Ki, &
                            LinNames=[(TRIM(LinStr)//TransDispSuffix(j), j=1, 3)] )
         END DO
      END DO

      ! velocity states
      DO il = 1, inst_%nlines
         ndofl = CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         N = ndofl/3 - 1
         DO i = 0, N - 2
            LinStr = 'Line '//TRIM(Num2LStr(il))//' node '//TRIM(Num2LStr(i + 1))
            CALL MV_AddVar( Vars%x, LinStr, FieldTransVel, DatLoc(CD_x_states), &
                            iAry=m%LineStateIs1(il) + 3*i, &
                            Num=3, Flags=VF_DerivOrder2, Perturb=0.1_R8Ki, &
                            LinNames=[(TRIM(LinStr)//TransVelSuffix(j), j=1, 3)] )
         END DO
      END DO

      ! acceleration states: the third mirror block per line reloads the committed
      ! generalised-alpha acceleration on checkpoint restart. Registered as ONE generic
      ! (FieldScalar) coverage block per line -- acceleration is not a linearization
      ! field, and the linearization registers no x variables; the goal is that Vars%x
      ! covers the WHOLE mirror so the framework pack/unpack sizing stays consistent
      ! (the Init-time size assertion below enforces exactly that).
      DO il = 1, inst_%nlines
         ndofl = CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         nin = ndofl - 6
         LinStr = 'Line '//TRIM(Num2LStr(il))
         ALLOCATE (CblLinVel(nin), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate line accel lin-name work array.', &
                             ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         DO j = 1, nin
            CblLinVel(j) = TRIM(LinStr)//' interior accel '//TRIM(Num2LStr(j))
         END DO
         CALL MV_AddVar( Vars%x, TRIM(LinStr)//' interior accel', FieldScalar, DatLoc(CD_x_states), &
                         iAry=m%LineStateIs1(il) + 2*nin, Num=nin, Perturb=0.1_R8Ki, &
                         LinNames=CblLinVel )
         DEALLOCATE (CblLinVel)
      END DO

      ! viscoelastic dl_1 states: the per-element series-Kelvin internal strain appended after the
      ! interior [rd; r; rdd] block for a viscoelastic line (mirror offset LineStateIs1 +
      ! 3*nin). Registered as a generic FieldScalar coverage block (not a linearization
      ! field) so Vars%x covers the WHOLE mirror -- the Init-time Vars%Nx == m%Nx assertion.
      DO il = 1, inst_%nlines
         IF (.NOT. CD_System_Line_Has_Viscoelastic( inst_%agg%sys%fast%system, il )) CYCLE
         ndofl = CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         nin = ndofl - 6
         nve = CD_System_Line_NElem( inst_%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         LinStr = 'Line '//TRIM(Num2LStr(il))
         ALLOCATE (CblLinVel(nve), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate line dl_1 lin-name work array.', &
                             ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         DO j = 1, nve
            CblLinVel(j) = TRIM(LinStr)//' viscoelastic dl_1 '//TRIM(Num2LStr(j))
         END DO
         CALL MV_AddVar( Vars%x, TRIM(LinStr)//' viscoelastic dl_1', FieldScalar, DatLoc(CD_x_states), &
                         iAry=m%LineStateIs1(il) + 3*nin, Num=nve, Perturb=0.1_R8Ki, &
                         LinNames=CblLinVel )
         DEALLOCATE (CblLinVel)
      END DO

      ! Syrope states: two per-element blocks (slow-spring static strain then running-
      ! maximum tension) appended after the interior block for a Syrope line (mirror
      ! offset LineStateIs1 + 3*nin; a line is viscoelastic OR Syrope, never both).
      ! Registered as a generic FieldScalar coverage block so Vars%x covers the WHOLE
      ! mirror -- the Init-time Vars%Nx == m%Nx assertion.
      DO il = 1, inst_%nlines
         IF (.NOT. CD_System_Line_Has_Syrope( inst_%agg%sys%fast%system, il )) CYCLE
         ndofl = CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         nin = ndofl - 6
         nve = CD_System_Line_NElem( inst_%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         LinStr = 'Line '//TRIM(Num2LStr(il))
         ALLOCATE (CblLinVel(2*nve), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate line Syrope lin-name work array.', &
                             ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         DO j = 1, nve
            CblLinVel(j)       = TRIM(LinStr)//' Syrope slow '//TRIM(Num2LStr(j))
            CblLinVel(nve + j) = TRIM(LinStr)//' Syrope tmax '//TRIM(Num2LStr(j))
         END DO
         CALL MV_AddVar( Vars%x, TRIM(LinStr)//' Syrope state', FieldScalar, DatLoc(CD_x_states), &
                         iAry=m%LineStateIs1(il) + 3*nin, Num=2*nve, Perturb=0.1_R8Ki, &
                         LinNames=CblLinVel )
         DEALLOCATE (CblLinVel)
      END DO

      ! finite-EI cable states: pack_state_mirror appends the free-DOF mirror
      ! [v_free; q_free] of each cable AFTER the EI=0 mooring block, at the x%states
      ! offset koff = LineStateIsN(nlines) (0 when there are no mooring lines). Register
      ! the SAME index range here so Vars%x covers the WHOLE mirror (the appended cable
      ! slots were previously invisible to the framework pack/unpack). The free DOFs of
      ! a cable mix a node translation r and tangent m -- they are NOT per-node
      ! 3-vectors -- so each cable contributes two GENERIC (FieldScalar) blocks, a
      ! velocity half then a position half, each sized to the free-DOF count. This is
      ! a size/index mirror (the linearization registers no x variables), so the
      ! goal is index coverage, not a physically-typed linearization surface.
      koff = 0
      IF (inst_%nlines > 0) koff = m%LineStateIsN(inst_%nlines)
      DO ic = 1, inst_%agg%ncable
         csz = CD_HFMF_MirrorSize( inst_%agg%cables(ic) )
         IF (csz <= 0) CYCLE                 ! uninitialised cable contributes no states
         ! free DOFs; mirror = [v_free; q_free; a_free] then the friction anchors (if any)
         nf = (csz - CD_HFMF_FrictionMirrorSize( inst_%agg%cables(ic) ) - &
               CD_HFMF_ForceMirrorSize( inst_%agg%cables(ic) ))/3
         LinStr = 'Cable '//TRIM(Num2LStr(ic))
         ! Build the per-DOF linearization names in an explicit pre-allocated loop, NOT an
         ! inline implied-do array constructor: an impure function (Num2LStr) evaluated inside
         ! a character array constructor mismanages compiler temporaries and corrupts the heap.
         ! This is the SAME failure class fixed for the WriteOutput var below -- fix the class,
         ! not the instance (the WriteOutput fix missed these two cable-var instances).
         ALLOCATE (CblLinVel(nf), CblLinPos(nf), CblLinAcc(nf), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate cable lin-name work arrays.', &
                             ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         DO j = 1, nf
            CblLinVel(j) = TRIM(LinStr)//' free-DOF vel '//TRIM(Num2LStr(j))
            CblLinPos(j) = TRIM(LinStr)//' free-DOF pos '//TRIM(Num2LStr(j))
            CblLinAcc(j) = TRIM(LinStr)//' free-DOF accel '//TRIM(Num2LStr(j))
         END DO
         ! velocity block -> x%states(koff+1 : koff+nf)
         CALL MV_AddVar( Vars%x, TRIM(LinStr)//' free-DOF vel', FieldScalar, DatLoc(CD_x_states), &
                         iAry=koff + 1, Num=nf, Flags=VF_DerivOrder2, DerivOrder=1, Perturb=0.1_R8Ki, &
                         LinNames=CblLinVel )
         ! position block -> x%states(koff+nf+1 : koff+2*nf)
         CALL MV_AddVar( Vars%x, TRIM(LinStr)//' free-DOF pos', FieldScalar, DatLoc(CD_x_states), &
                         iAry=koff + nf + 1, Num=nf, Flags=VF_DerivOrder2, DerivOrder=0, Perturb=0.05_R8Ki, &
                         LinNames=CblLinPos )
         ! acceleration block -> x%states(koff+2*nf+1 : koff+3*nf): the committed
         ! gen-alpha acceleration, carried for exact checkpoint restart (coverage-only,
         ! like the mooring interior-accel block above)
         CALL MV_AddVar( Vars%x, TRIM(LinStr)//' free-DOF accel', FieldScalar, DatLoc(CD_x_states), &
                         iAry=koff + 2*nf + 1, Num=nf, Perturb=0.1_R8Ki, &
                         LinNames=CblLinAcc )
         DEALLOCATE (CblLinVel, CblLinPos, CblLinAcc)
         ! stick-slip friction anchors (x, y per node) and the force-blend committed force
         ! (flag + force) after the acceleration block: coverage-only, like the acceleration block
         IF (csz > 3*nf) THEN
            ALLOCATE (CblLinVel(csz - 3*nf), STAT=ErrStat2)
            IF (ErrStat2 /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, 'Could not allocate cable state lin-name work array.', &
                                ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            DO j = 1, csz - 3*nf
               CblLinVel(j) = TRIM(LinStr)//' integrator state '//TRIM(Num2LStr(j))
            END DO
            CALL MV_AddVar( Vars%x, TRIM(LinStr)//' integrator state', FieldScalar, DatLoc(CD_x_states), &
                            iAry=koff + 3*nf + 1, Num=csz - 3*nf, Perturb=0.1_R8Ki, LinNames=CblLinVel )
            DEALLOCATE (CblLinVel)
         END IF
         koff = koff + csz
      END DO

      ! coupled Rigid6 6-DOF body state: coverage-only registration after the cable blocks (Vars%Nx
      ! must equal the mirror size m%Nx). Names built in a pre-allocated loop, NOT an implied-do with
      ! an impure Num2LStr inside an array constructor (that corrupts compiler temporaries -- the heap
      ! class fixed for the cable/WriteOutput vars).
      IF (CD_AGG_HasRigid6( inst_%agg )) THEN
         BLOCK
            INTEGER(IntKi) :: r6sz, j2
            r6sz = CD_AGG_Rigid6_MirrorSize( inst_%agg )
            IF (r6sz > 0) THEN
               ALLOCATE (CblLinVel(r6sz), STAT=ErrStat2)
               IF (ErrStat2 /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, 'Could not allocate Rigid6 lin-name work array.', &
                                   ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
               DO j2 = 1, r6sz
                  CblLinVel(j2) = 'Rigid6 body-state '//TRIM(Num2LStr(j2))
               END DO
               CALL MV_AddVar( Vars%x, 'Rigid6 body-state', FieldScalar, DatLoc(CD_x_states), &
                               iAry=koff + 1, Num=r6sz, Perturb=0.1_R8Ki, LinNames=CblLinVel )
               DEALLOCATE (CblLinVel)
               koff = koff + r6sz
            END IF
         END BLOCK
      END IF

      IF (CD_AGG_HasRod( inst_%agg )) THEN
         BLOCK
            INTEGER(IntKi) :: rdsz, j2
            rdsz = CD_AGG_Rod_MirrorSize( inst_%agg )
            IF (rdsz > 0) THEN
               ALLOCATE (CblLinVel(rdsz), STAT=ErrStat2)
               IF (ErrStat2 /= 0) THEN
                  CALL SetErrStat( ErrID_Fatal, 'Could not allocate rigid-rod lin-name work array.', &
                                   ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
               DO j2 = 1, rdsz
                  CblLinVel(j2) = 'Rigid rod state '//TRIM(Num2LStr(j2))
               END DO
               CALL MV_AddVar( Vars%x, 'Rigid rod state', FieldScalar, DatLoc(CD_x_states), &
                               iAry=koff + 1, Num=rdsz, Perturb=0.1_R8Ki, LinNames=CblLinVel )
               DEALLOCATE (CblLinVel)
               koff = koff + rdsz
            END IF
         END BLOCK
      END IF

      ! EI=0 line stick-slip friction anchors: coverage-only registration after the rod block.
      BLOCK
         INTEGER(IntKi) :: frsz, j2
         frsz = line_friction_mirror_size( inst_ )
         IF (frsz > 0) THEN
            ALLOCATE (CblLinVel(frsz), STAT=ErrStat2)
            IF (ErrStat2 /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, 'Could not allocate line friction lin-name work array.', &
                                ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            DO j2 = 1, frsz
               CblLinVel(j2) = 'Line friction anchor '//TRIM(Num2LStr(j2))
            END DO
            CALL MV_AddVar( Vars%x, 'Line friction anchor', FieldScalar, DatLoc(CD_x_states), &
                            iAry=koff + 1, Num=frsz, Perturb=0.1_R8Ki, LinNames=CblLinVel )
            DEALLOCATE (CblLinVel)
            koff = koff + frsz
         END IF
      END BLOCK

      ! the compact coupled-motion vector [q; v; a] (per-LINE-END slots) of every deck with
      ! EI=0 lines, then block 2 of dynamic (Free/Connect) points. Coverage-only
      ! registration like every trailing block -- Vars%Nx must equal the mirror size m%Nx.
      IF (inst_%agg%has_sys) THEN
         BLOCK
            INTEGER(IntKi) :: ncd2
            ncd2 = CD_System_NSystemCoupledDOF( inst_%agg%sys%fast%system )
            ALLOCATE (CblLinVel(ncd2), CblLinPos(ncd2), CblLinAcc(ncd2), STAT=ErrStat2)
            IF (ErrStat2 /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, 'Could not allocate coupled-vector lin-name work arrays.', &
                                ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            DO j = 1, ncd2
               CblLinPos(j) = 'CoupledVec pos '//TRIM(Num2LStr(j))
               CblLinVel(j) = 'CoupledVec vel '//TRIM(Num2LStr(j))
               CblLinAcc(j) = 'CoupledVec accel '//TRIM(Num2LStr(j))
            END DO
            CALL MV_AddVar( Vars%x, 'CoupledVec pos', FieldScalar, DatLoc(CD_x_states), &
                            iAry=koff + 1, Num=ncd2, Perturb=0.05_R8Ki, LinNames=CblLinPos )
            CALL MV_AddVar( Vars%x, 'CoupledVec vel', FieldScalar, DatLoc(CD_x_states), &
                            iAry=koff + ncd2 + 1, Num=ncd2, Perturb=0.1_R8Ki, LinNames=CblLinVel )
            CALL MV_AddVar( Vars%x, 'CoupledVec accel', FieldScalar, DatLoc(CD_x_states), &
                            iAry=koff + 2*ncd2 + 1, Num=ncd2, Perturb=0.1_R8Ki, LinNames=CblLinAcc )
            DEALLOCATE (CblLinVel, CblLinPos, CblLinAcc)
            koff = koff + 3*ncd2
         END BLOCK
      END IF
      ! dynamic (Free/Connect) points, block 2: the trailing per-point [v; q; a] mirror block.
      ! Coverage-only index registration (like the accel blocks above) so Vars%Nx equals
      ! the mirror size m%Nx -- without it the drift check below kills EVERY dynamic-
      ! point deck at Init, not just restarts. LinNames use the same explicit
      ! pre-allocated loop as the cable blocks (the impure-function-in-array-constructor
      ! heap-corruption class).
      DO j = 1, CD_AGG_NDynamicPoints( inst_%agg )
         LinStr = 'DynPoint '//TRIM(Num2LStr(j))
         ALLOCATE (CblLinVel(3), CblLinPos(3), CblLinAcc(3), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate dynamic-point lin-name work arrays.', &
                             ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         DO ic = 1, 3
            CblLinVel(ic) = TRIM(LinStr)//' vel '//TRIM(Num2LStr(ic))
            CblLinPos(ic) = TRIM(LinStr)//' pos '//TRIM(Num2LStr(ic))
            CblLinAcc(ic) = TRIM(LinStr)//' accel '//TRIM(Num2LStr(ic))
         END DO
         CALL MV_AddVar( Vars%x, TRIM(LinStr)//' vel', FieldScalar, DatLoc(CD_x_states), &
                         iAry=koff + 1, Num=3, Flags=VF_DerivOrder2, DerivOrder=1, Perturb=0.1_R8Ki, &
                         LinNames=CblLinVel )
         CALL MV_AddVar( Vars%x, TRIM(LinStr)//' pos', FieldScalar, DatLoc(CD_x_states), &
                         iAry=koff + 4, Num=3, Flags=VF_DerivOrder2, DerivOrder=0, Perturb=0.05_R8Ki, &
                         LinNames=CblLinPos )
         CALL MV_AddVar( Vars%x, TRIM(LinStr)//' accel', FieldScalar, DatLoc(CD_x_states), &
                         iAry=koff + 7, Num=3, Perturb=0.1_R8Ki, &
                         LinNames=CblLinAcc )
         DEALLOCATE (CblLinVel, CblLinPos, CblLinAcc)
         koff = koff + 9
      END DO
      ! FAILURE fired flags: coverage-only registration of the mirror tail (0/1 per
      ! row) so Vars%Nx equals m%Nx -- the drift-guard class. Same explicit
      ! pre-allocated LinNames loop as every block above.
      IF (CD_AGG_NFailures( inst_%agg ) > 0) THEN
         BLOCK
            INTEGER(IntKi) :: nfl2
            nfl2 = CD_AGG_NFailures( inst_%agg )
            ALLOCATE (CblLinPos(nfl2), STAT=ErrStat2)
            IF (ErrStat2 /= 0) THEN
               CALL SetErrStat( ErrID_Fatal, 'Could not allocate failure-flag lin-name work array.', &
                                ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            DO j = 1, nfl2
               CblLinPos(j) = 'FailFlag '//TRIM(Num2LStr(j))
            END DO
            CALL MV_AddVar( Vars%x, 'FailFlag', FieldScalar, DatLoc(CD_x_states), &
                            iAry=koff + 1, Num=nfl2, Perturb=0.1_R8Ki, &
                            LinNames=CblLinPos )
            DEALLOCATE (CblLinPos)
            koff = koff + nfl2
         END BLOCK
      END IF

      END IF   ! .NOT. Linearize (the x-variable registration)

      ! input variables: the coupled-kinematics mesh
      ALLOCATE (Vars%u(0), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate Vars%u.', ErrStat, ErrMsg, RoutineName ); RETURN
      END IF
      ! one mesh variable per turbine (FAST.Farm); plain mode is the 1-mesh case
      DO ivT = 1, MAX(1, inst_%nT)
         CALL MV_AddMeshVar( Vars%u, "CoupledKinematics", MotionFields, &
                             DatLoc(CD_u_CoupledKinematics, ivT), &
                             Mesh=u%CoupledKinematics(ivT), &
                             Perturbs=[0.05_R8Ki, &   ! FieldTransDisp
                                       0.1_R8Ki, &    ! FieldOrientation
                                       0.1_R8Ki, &    ! FieldTransVel
                                       0.1_R8Ki, &    ! FieldAngularVel
                                       0.1_R8Ki, &    ! FieldTransAcc
                                       0.1_R8Ki] )    ! FieldAngularAcc
      END DO

      ! output variables: the coupled-loads mesh, then the WriteOutput channels
      ALLOCATE (Vars%y(0), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate Vars%y.', ErrStat, ErrMsg, RoutineName ); RETURN
      END IF
      DO ivT = 1, MAX(1, inst_%nT)
         CALL MV_AddMeshVar( Vars%y, "CoupledLoads", LoadFields, &
                             DatLoc(CD_y_CoupledLoads, ivT), &
                             Mesh=y%CoupledLoads(ivT) )
      END DO

      ! WriteOutput channels (the OrcaFlex-vocabulary deck OUTPUTS), registered as a scalar
      ! array so linearization / visualization see them -- the same y-variable MoorDyn adds for
      ! its WriteOutput. Registered ONLY when the deck declares channels, so a channel-free deck
      ! keeps Vars%y carrying just the loads mesh (byte-identical to the pre-channel module) and
      ! MV_InitVarsJac counts no extra y-slots. The linearization names ("Header, Unit") are built
      ! into an explicit array first (the aggregate owns the header/unit; InitOut is out of scope).
      IF (p%NumOuts > 0) THEN
         ALLOCATE (WrOutLinNames(p%NumOuts), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate WriteOutput lin-name work array.', &
                             ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         DO j = 1, p%NumOuts
            CALL CD_AGG_ChannelHeader( inst_%agg, j, whdr, wunt, es, em )
            WrOutLinNames(j) = TRIM(whdr)//', '//TRIM(wunt)
         END DO
         CALL MV_AddVar( Vars%y, "WriteOutput", FieldScalar, DatLoc(CD_y_WriteOutput), &
                         Flags=VF_WriteOut, Num=p%NumOuts, LinNames=WrOutLinNames )
         DEALLOCATE (WrOutLinNames)
      END IF

      ! finalize the variable metadata + the Jacobian workspace in m%Jac
      CALL MV_InitVarsJac( Vars, m%Jac, Linearize, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! Self-consistency guard: the registered continuous-state variables MUST cover
      ! EXACTLY the state mirror x%states(1:m%Nx) that pack_state_mirror writes.
      ! MV_InitVarsJac has set Vars%Nx = sum(Vars%x%Num); if a future change grows the
      ! mirror (a new state block) without a matching MV_AddVar here (or shrinks it),
      ! the packed x-variable size and m%Nx drift apart -- fail LOUDLY at init instead
      ! of silently leaving mirror slots invisible to pack/unpack/checkpoint.
      IF (.NOT. Linearize .AND. Vars%Nx /= m%Nx) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the registered continuous-state '// &
                          'variable size ('//TRIM(Num2LStr(Vars%Nx))//') does not equal the state '// &
                          'mirror size m%Nx ('//TRIM(Num2LStr(m%Nx))//'); the ModVars x-registration '// &
                          'and the pack_state_mirror layout have drifted apart.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      ! linearization workspaces used by the FD Jacobian
      CALL CD_CopyInput( u, m%u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL CD_CopyOutput( y, m%y_lin, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL CD_CopyContState( x, m%x_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL CD_CopyContState( x, m%dxdt_lin, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END SUBROUTINE CD_InitVars

   SUBROUTINE reload_interiors_ary(inst_, m, states, t_clock, ErrStat, ErrMsg)
      !! Reload the committed solver state from a mirror array (the checkpoint
      !! x%states layout): interior [rd; r; rdd] per mooring line -- the acceleration
      !! block restores the gen-alpha COMMITTED acceleration exactly (a re-derivation
      !! would move the first restored load evaluation) -- and each cable's
      !! [v; q; a] free-DOF mirror. Shared by the checkpoint rebuild and the
      !! state-const CalcOutput restore.
      TYPE(CD_OF_Instance),         INTENT(INOUT) :: inst_
      TYPE(CD_MiscVarType),         INTENT(IN)    :: m
      REAL(R8Ki),                   INTENT(IN)    :: states(:)
      REAL(wp),                     INTENT(IN)    :: t_clock
      INTEGER(IntKi),               INTENT(OUT)   :: ErrStat
      CHARACTER(*),                 INTENT(OUT)   :: ErrMsg

      character(*), parameter :: RoutineName = 'CableDyn:reload_interiors'
      INTEGER(IntKi)          :: es, il, ndofl, nin, k, ic, csz, koff, npd, ncd, nfl, nve
      CHARACTER(512)          :: em

      ErrStat = ErrID_None
      ErrMsg  = ''
      DO il = 1, inst_%nlines
         ndofl = CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         nin = ndofl - 6
         k = m%LineStateIs1(il)
         ! a one-segment line (NumSegs = 1, n_dof = 6) has NO interior DOFs: nothing
         ! to reload for the interior, and the model-level interior update rightly
         ! rejects zero-length arrays -- guard the interior block, but STILL reload a
         ! viscoelastic dl_1 tail below (a 1-segment line has one element with dl_1)
         IF (nin > 0) THEN
            IF (k + 3*nin - 1 > SIZE(states)) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: line '//TRIM(Num2LStr(il))// &
                                ' state-mirror extent does not match the mesh.', ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
            IF (ALLOCATED(inst_%ws_qL)) THEN
               IF (SIZE(inst_%ws_qL) < nin) DEALLOCATE (inst_%ws_qL, inst_%ws_vL, inst_%ws_aL)
            END IF
            IF (.NOT. ALLOCATED(inst_%ws_qL)) ALLOCATE (inst_%ws_qL(ndofl), inst_%ws_vL(ndofl), inst_%ws_aL(ndofl))
            inst_%ws_vL(1:nin) = REAL(states(k:k + nin - 1), wp)
            inst_%ws_qL(1:nin) = REAL(states(k + nin:k + 2*nin - 1), wp)
            inst_%ws_aL(1:nin) = REAL(states(k + 2*nin:k + 3*nin - 1), wp)
            CALL CD_Update_System_Line_Interior_State( inst_%agg%sys%fast%system, il, &
                                                       inst_%ws_qL(1:nin), inst_%ws_vL(1:nin), es, em, &
                                                       a_interior=inst_%ws_aL(1:nin) )
            IF (es /= CD_SYSTEM_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
         END IF
         ! viscoelastic dl_1 reload: the per-element tail packed after the interior block
         ! (mirrors pack_state_mirror; independent of nin so 1-segment lines restore too)
         IF (CD_System_Line_Has_Viscoelastic( inst_%agg%sys%fast%system, il )) THEN
            nve = CD_System_Line_NElem( inst_%agg%sys%fast%system, il, es, em )
            IF (es /= CD_SYSTEM_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            IF (k + 3*nin + nve - 1 > SIZE(states)) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: line '//TRIM(Num2LStr(il))// &
                                ' viscoelastic dl_1 mirror extends past the state vector.', &
                                ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
            BLOCK
               REAL(wp) :: dl1(nve)
               dl1 = REAL(states(k + 3*nin:k + 3*nin + nve - 1), wp)
               CALL CD_Set_System_Line_VE_Dl1( inst_%agg%sys%fast%system, il, dl1, es, em )
               IF (es /= CD_SYSTEM_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
            END BLOCK
         END IF
         ! Syrope state reload: two per-element tail blocks (slow then tmax) after the
         ! interior block (mirrors pack_state_mirror; independent of nin)
         IF (CD_System_Line_Has_Syrope( inst_%agg%sys%fast%system, il )) THEN
            nve = CD_System_Line_NElem( inst_%agg%sys%fast%system, il, es, em )
            IF (es /= CD_SYSTEM_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            IF (k + 3*nin + 2*nve - 1 > SIZE(states)) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: line '//TRIM(Num2LStr(il))// &
                                ' Syrope state mirror extends past the state vector.', ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
            BLOCK
               REAL(wp) :: slw(nve), tmx(nve)
               slw = REAL(states(k + 3*nin:k + 3*nin + nve - 1), wp)
               tmx = REAL(states(k + 3*nin + nve:k + 3*nin + 2*nve - 1), wp)
               CALL CD_Set_System_Line_Syrope_State( inst_%agg%sys%fast%system, il, slw, tmx, es, em )
               IF (es /= CD_SYSTEM_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
            END BLOCK
         END IF
      END DO
      koff = 0
      IF (inst_%nlines > 0) koff = m%LineStateIsN(inst_%nlines)
      DO ic = 1, inst_%agg%ncable
         csz = CD_HFMF_MirrorSize( inst_%agg%cables(ic) )
         IF (csz <= 0) CYCLE
         IF (koff + csz > SIZE(states)) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: cable '//TRIM(Num2LStr(ic))// &
                             ' mirror extends past the state vector.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         IF (ALLOCATED(inst_%ws_cbuf)) THEN
            IF (SIZE(inst_%ws_cbuf) /= csz) DEALLOCATE (inst_%ws_cbuf)
         END IF
         IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(csz))
         inst_%ws_cbuf = REAL(states(koff + 1:koff + csz), wp)
         CALL CD_HFMF_UnpackMirror( inst_%agg%cables(ic), inst_%ws_cbuf, t_clock, es, em )
         IF (es /= CD_HFMF_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         koff = koff + csz
      END DO
      ! coupled Rigid6 6-DOF body state: restore after the cable block, before the dynamic points
      ! (the body drives its attachment points, so the restored pose re-derives them next advance)
      IF (CD_AGG_HasRigid6( inst_%agg )) THEN
         csz = CD_AGG_Rigid6_MirrorSize( inst_%agg )
         IF (csz > 0) THEN
            IF (koff + csz > SIZE(states)) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the Rigid6 body-state mirror '// &
                                'extends past the state vector.', ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            IF (ALLOCATED(inst_%ws_cbuf)) THEN
               IF (SIZE(inst_%ws_cbuf) < csz) DEALLOCATE (inst_%ws_cbuf)
            END IF
            IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(csz))
            inst_%ws_cbuf(1:csz) = REAL(states(koff + 1:koff + csz), wp)
            CALL CD_AGG_Set_Rigid6_States( inst_%agg, inst_%ws_cbuf(1:csz), es, em )
            IF (es /= CD_AGG_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            koff = koff + csz
         END IF
      END IF
      ! Coupled rigid-rod committed state follows the Rigid6 block. Set scatters
      ! Rod<N>A/B kinematics immediately, keeping restart-time output consistent.
      IF (CD_AGG_HasRod( inst_%agg )) THEN
         csz = CD_AGG_Rod_MirrorSize( inst_%agg )
         IF (csz > 0) THEN
            IF (koff + csz > SIZE(states)) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the rigid-rod state mirror '// &
                                'extends past the state vector.', ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            IF (ALLOCATED(inst_%ws_cbuf)) THEN
               IF (SIZE(inst_%ws_cbuf) < csz) DEALLOCATE (inst_%ws_cbuf)
            END IF
            IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(csz))
            inst_%ws_cbuf(1:csz) = REAL(states(koff + 1:koff + csz), wp)
            CALL CD_AGG_Set_Rod_States( inst_%agg, inst_%ws_cbuf(1:csz), es, em )
            IF (es /= CD_AGG_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            koff = koff + csz
         END IF
      END IF
      ! EI=0 line stick-slip friction anchors, line by line after the rod block
      IF (line_friction_mirror_size( inst_ ) > 0) THEN
         DO il = 1, inst_%nlines
            IF (.NOT. CD_System_Line_Has_Friction( inst_%agg%sys%fast%system, il )) CYCLE
            ndofl = CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )
            IF (es /= CD_SYSTEM_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            csz = 2*(ndofl/3)
            IF (koff + csz > SIZE(states)) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: line '//TRIM(Num2LStr(il))// &
                                ' friction-anchor mirror extends past the state vector.', ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
            BLOCK
               REAL(wp) :: anch(2, ndofl/3)
               anch = RESHAPE(REAL(states(koff + 1:koff + csz), wp), [2, ndofl/3])
               CALL CD_Set_System_Line_Friction_Anchors( inst_%agg%sys%fast%system, il, anch, es, em )
               IF (es /= CD_SYSTEM_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
            END BLOCK
            koff = koff + csz
         END DO
      END IF
      ! dynamic (Free/Connect) points: the per-point store FIRST (its self-sufficient
      ! scatter positions the point slots), then the compact coupled-motion vector LAST
      ! -- the CHECKPOINTED ENDPOINT TRUTH must win on every coupled slot, so it is
      ! applied after any store scatter could touch them (at pack time store and vector
      ! point slots are bit-equal by the step's commit contract, but the reload order
      ! guarantees the vector's values regardless).
      npd = 9*CD_AGG_NDynamicPoints( inst_%agg )
      IF (npd == 0 .AND. inst_%agg%has_sys) THEN
         ! no dynamic points: restore the committed boundary kinematics (the restored input
         ! mesh carries the coupled positions and velocities, not the committed accelerations)
         ncd = CD_System_NSystemCoupledDOF( inst_%agg%sys%fast%system )
         IF (koff + 3*ncd > SIZE(states)) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the coupled-motion mirror extends '// &
                             'past the state vector.', ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         IF (ALLOCATED(inst_%ws_cbuf)) THEN
            IF (SIZE(inst_%ws_cbuf) < 3*ncd) DEALLOCATE (inst_%ws_cbuf)
         END IF
         IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(3*ncd))
         inst_%ws_cbuf(1:3*ncd) = REAL(states(koff + 1:koff + 3*ncd), wp)
         CALL CD_Update_System_CoupledMotion( inst_%agg%sys%fast%system, inst_%ws_cbuf(1:ncd), &
                                              inst_%ws_cbuf(ncd + 1:2*ncd), inst_%ws_cbuf(2*ncd + 1:3*ncd), es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         koff = koff + 3*ncd
      END IF
      IF (npd > 0 .AND. inst_%agg%has_sys) THEN
         ncd = CD_System_NSystemCoupledDOF( inst_%agg%sys%fast%system )
         nfl = CD_AGG_NFailures( inst_%agg )
         IF (koff + 3*ncd + npd + nfl > SIZE(states)) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the dynamic-point mirror extends '// &
                             'past the state vector.', ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         ! FAILURE flags replay FIRST (they live at the mirror tail): the detach
         ! topology must match the checkpointed exchange layout BEFORE any state
         ! block is written -- the restored vector/store were packed on the
         ! post-detach system, and writing them through a pre-detach map would put
         ! committed values on the wrong slots.
         IF (nfl > 0) THEN
            BLOCK
               LOGICAL :: flg(nfl)
               flg = states(koff + 3*ncd + npd + 1:koff + 3*ncd + npd + nfl) > 0.5_DbKi
               CALL CD_AGG_Set_Failure_Flags( inst_%agg, flg, es, em )
               IF (es /= CD_AGG_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
            END BLOCK
         END IF
         IF (ALLOCATED(inst_%ws_cbuf)) THEN
            IF (SIZE(inst_%ws_cbuf) < MAX(3*ncd, npd)) DEALLOCATE (inst_%ws_cbuf)
         END IF
         IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(MAX(3*ncd, npd)))
         inst_%ws_cbuf(1:npd) = REAL(states(koff + 3*ncd + 1:koff + 3*ncd + npd), wp)
         CALL CD_Set_System_DynamicPoint_States( inst_%agg%sys%fast%system, inst_%ws_cbuf(1:npd), es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         inst_%ws_cbuf(1:3*ncd) = REAL(states(koff + 1:koff + 3*ncd), wp)
         CALL CD_Update_System_CoupledMotion( inst_%agg%sys%fast%system, inst_%ws_cbuf(1:ncd), &
                                              inst_%ws_cbuf(ncd + 1:2*ncd), inst_%ws_cbuf(2*ncd + 1:3*ncd), es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         ! the FMF facade's point mesh caches the coupled surface; the moving-step path
         ! seeds from it, so re-derive it from the just-restored system state (else the
         ! next step would feed pre-restore Free/Connect values back into the system)
         CALL CD_AGG_Refresh_PointMesh( inst_%agg, es, em )
         IF (es /= CD_AGG_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
      END IF
   END SUBROUTINE reload_interiors_ary

   SUBROUTINE pack_state_mirror(inst_, m, x, ErrStat, ErrMsg)
      TYPE(CD_OF_Instance),         INTENT(INOUT) :: inst_
      TYPE(CD_MiscVarType),         INTENT(IN)    :: m
      TYPE(CD_ContinuousStateType), INTENT(INOUT) :: x
      INTEGER(IntKi),               INTENT(OUT)   :: ErrStat
      CHARACTER(*),                 INTENT(OUT)   :: ErrMsg
      CALL pack_state_mirror_ary(inst_, m, x%states, ErrStat, ErrMsg)
   END SUBROUTINE pack_state_mirror

   INTEGER(IntKi) FUNCTION line_friction_mirror_size(inst_) RESULT(n)
      !! Mirror slots of the EI=0 lines' stick-slip friction anchors: (x, y) per node of
      !! every line on a frictional seabed (0 without such lines).
      TYPE(CD_OF_Instance), INTENT(IN) :: inst_
      INTEGER(IntKi) :: il, es
      CHARACTER(256) :: em
      n = 0
      IF (.NOT. inst_%agg%has_sys) RETURN
      DO il = 1, inst_%nlines
         IF (.NOT. CD_System_Line_Has_Friction( inst_%agg%sys%fast%system, il )) CYCLE
         n = n + 2*(CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )/3)
      END DO
   END FUNCTION line_friction_mirror_size

   !> Copy the instance states into the framework mirror x%states: per EI=0 mooring
   !> line, the interior-node velocities then the interior-node positions (the MoorDyn
   !> layout; end nodes are excluded -- they are prescribed or fixed); then, after the
   !> mooring block, each finite-EI cable contributes its free-DOF velocities then
   !> positions [v_free; q_free]. The module stays self-authoritative: it reloads its
   !> state from the mirror only on a checkpoint rebuild or the state-const CalcOutput
   !> restore (reload_interiors_ary).
   SUBROUTINE pack_state_mirror_ary(inst_, m, states, ErrStat, ErrMsg)
      !! Pack the committed solver state (interior [rd; r; rdd] per mooring line,
      !! [v; q; a] free DOFs per cable) into the given mirror array -- the checkpoint
      !! x%states layout. Also used by the state-const CalcOutput to capture the
      !! committed state before an out-of-band transfer.
      TYPE(CD_OF_Instance),         INTENT(INOUT) :: inst_
      TYPE(CD_MiscVarType),         INTENT(IN)    :: m
      REAL(R8Ki),                   INTENT(INOUT) :: states(:)
      INTEGER(IntKi),               INTENT(OUT)   :: ErrStat
      CHARACTER(*),                 INTENT(OUT)   :: ErrMsg

      character(*), parameter :: RoutineName = 'CableDyn:pack_state_mirror'
      INTEGER(IntKi)          :: es, il, ndofl, nin, k, ic, csz, koff, npd, ncd, nfl, nve
      CHARACTER(512)          :: em

      ErrStat = ErrID_None
      ErrMsg  = ''
      DO il = 1, inst_%nlines
         ndofl = CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         ! grow-only INSTANCE scratch: this routine runs once per mooring step, so the
         ! buffers must persist across calls (as locals they re-allocated every call).
         IF (ALLOCATED(inst_%ws_qL)) THEN
            IF (SIZE(inst_%ws_qL) /= ndofl) DEALLOCATE (inst_%ws_qL, inst_%ws_vL, inst_%ws_aL)
         END IF
         IF (.NOT. ALLOCATED(inst_%ws_qL)) ALLOCATE (inst_%ws_qL(ndofl), inst_%ws_vL(ndofl), inst_%ws_aL(ndofl))
         CALL CD_Get_System_Line_State( inst_%agg%sys%fast%system, il, inst_%ws_qL, inst_%ws_vL, inst_%ws_aL, es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         nin = ndofl - 6                       ! interior DOFs (drop both end nodes)
         k = m%LineStateIs1(il)
         states(k:k + nin - 1)           = REAL(inst_%ws_vL(4:ndofl - 3), DbKi)
         states(k + nin:k + 2*nin - 1)   = REAL(inst_%ws_qL(4:ndofl - 3), DbKi)
         ! the committed gen-alpha acceleration: carried so a checkpoint restart reloads
         ! the integrator state EXACTLY (re-deriving it moves the first restored loads)
         states(k + 2*nin:k + 3*nin - 1) = REAL(inst_%ws_aL(4:ndofl - 3), DbKi)
         ! viscoelastic dl_1 history: the per-element tail of this line's block, sized
         ! to match the Init nstates layout (LineStateIs1..IsN spans it)
         IF (CD_System_Line_Has_Viscoelastic( inst_%agg%sys%fast%system, il )) THEN
            nve = CD_System_Line_NElem( inst_%agg%sys%fast%system, il, es, em )
            IF (es /= CD_SYSTEM_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            BLOCK
               REAL(wp) :: dl1(nve)
               CALL CD_Get_System_Line_VE_Dl1( inst_%agg%sys%fast%system, il, dl1, es, em )
               IF (es /= CD_SYSTEM_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
               states(k + 3*nin:k + 3*nin + nve - 1) = REAL(dl1, DbKi)
            END BLOCK
         END IF
         ! Syrope state: two per-element tail blocks (slow then tmax), sized to match
         ! the Init nstates layout (LineStateIs1..IsN spans it)
         IF (CD_System_Line_Has_Syrope( inst_%agg%sys%fast%system, il )) THEN
            nve = CD_System_Line_NElem( inst_%agg%sys%fast%system, il, es, em )
            IF (es /= CD_SYSTEM_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            BLOCK
               REAL(wp) :: slw(nve), tmx(nve)
               CALL CD_Get_System_Line_Syrope_State( inst_%agg%sys%fast%system, il, slw, tmx, es, em )
               IF (es /= CD_SYSTEM_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
               states(k + 3*nin:k + 3*nin + nve - 1)             = REAL(slw, DbKi)
               states(k + 3*nin + nve:k + 3*nin + 2*nve - 1)     = REAL(tmx, DbKi)
            END BLOCK
         END IF
      END DO
      ! finite-EI cables: append per-cable free-DOF [v_free; q_free; a_free] after the mooring block
      koff = 0
      IF (inst_%nlines > 0) koff = m%LineStateIsN(inst_%nlines)
      DO ic = 1, inst_%agg%ncable
         csz = CD_HFMF_MirrorSize( inst_%agg%cables(ic) )
         IF (csz <= 0) CYCLE
         IF (ALLOCATED(inst_%ws_cbuf)) THEN
            IF (SIZE(inst_%ws_cbuf) /= csz) DEALLOCATE (inst_%ws_cbuf)
         END IF
         IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(csz))
         CALL CD_HFMF_PackMirror( inst_%agg%cables(ic), inst_%ws_cbuf, es, em )
         IF (es /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         states(koff + 1:koff + csz) = REAL(inst_%ws_cbuf, DbKi)
         koff = koff + csz
      END DO
      ! coupled Rigid6 6-DOF body state: pack after the cable block (koff carries the running offset)
      IF (CD_AGG_HasRigid6( inst_%agg )) THEN
         csz = CD_AGG_Rigid6_MirrorSize( inst_%agg )
         IF (csz > 0) THEN
            IF (ALLOCATED(inst_%ws_cbuf)) THEN
               IF (SIZE(inst_%ws_cbuf) < csz) DEALLOCATE (inst_%ws_cbuf)
            END IF
            IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(csz))
            CALL CD_AGG_Get_Rigid6_States( inst_%agg, inst_%ws_cbuf(1:csz), es, em )
            IF (es /= CD_AGG_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            states(koff + 1:koff + csz) = REAL(inst_%ws_cbuf(1:csz), DbKi)
            koff = koff + csz
         END IF
      END IF
      IF (CD_AGG_HasRod( inst_%agg )) THEN
         csz = CD_AGG_Rod_MirrorSize( inst_%agg )
         IF (csz > 0) THEN
            IF (ALLOCATED(inst_%ws_cbuf)) THEN
               IF (SIZE(inst_%ws_cbuf) < csz) DEALLOCATE (inst_%ws_cbuf)
            END IF
            IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(csz))
            CALL CD_AGG_Get_Rod_States( inst_%agg, inst_%ws_cbuf(1:csz), es, em )
            IF (es /= CD_AGG_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            states(koff + 1:koff + csz) = REAL(inst_%ws_cbuf(1:csz), DbKi)
            koff = koff + csz
         END IF
      END IF
      ! EI=0 line stick-slip friction anchors, line by line after the rod block
      IF (line_friction_mirror_size( inst_ ) > 0) THEN
         DO il = 1, inst_%nlines
            IF (.NOT. CD_System_Line_Has_Friction( inst_%agg%sys%fast%system, il )) CYCLE
            ndofl = CD_System_Line_NDOF( inst_%agg%sys%fast%system, il, es, em )
            IF (es /= CD_SYSTEM_OK) THEN
               CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            csz = 2*(ndofl/3)
            BLOCK
               REAL(wp) :: anch(2, ndofl/3)
               CALL CD_Get_System_Line_Friction_Anchors( inst_%agg%sys%fast%system, il, anch, es, em )
               IF (es /= CD_SYSTEM_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
               states(koff + 1:koff + csz) = REAL(RESHAPE(anch, [csz]), DbKi)
            END BLOCK
            koff = koff + csz
         END DO
      END IF
      ! dynamic (Free/Connect) points: the compact coupled-motion vector [q; v; a]
      ! (the committed boundary/endpoint kinematics -- see the sizing comment at Init),
      ! then the per-point store
      npd = 9*CD_AGG_NDynamicPoints( inst_%agg )
      IF (npd == 0 .AND. inst_%agg%has_sys) THEN
         ! no dynamic points: the committed boundary kinematics alone
         ncd = CD_System_NSystemCoupledDOF( inst_%agg%sys%fast%system )
         IF (koff + 3*ncd > SIZE(states)) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the coupled-motion mirror extends '// &
                             'past the state vector.', ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         IF (ALLOCATED(inst_%ws_cbuf)) THEN
            IF (SIZE(inst_%ws_cbuf) < 3*ncd) DEALLOCATE (inst_%ws_cbuf)
         END IF
         IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(3*ncd))
         CALL CD_Get_System_CoupledMotion( inst_%agg%sys%fast%system, inst_%ws_cbuf(1:ncd), &
                                           inst_%ws_cbuf(ncd + 1:2*ncd), inst_%ws_cbuf(2*ncd + 1:3*ncd), es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         states(koff + 1:koff + 3*ncd) = REAL(inst_%ws_cbuf(1:3*ncd), DbKi)
         koff = koff + 3*ncd
      END IF
      IF (npd > 0 .AND. inst_%agg%has_sys) THEN
         ncd = CD_System_NSystemCoupledDOF( inst_%agg%sys%fast%system )
         IF (koff + 3*ncd + npd > SIZE(states)) THEN
            CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the dynamic-point mirror extends '// &
                             'past the state vector.', ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         IF (ALLOCATED(inst_%ws_cbuf)) THEN
            IF (SIZE(inst_%ws_cbuf) < MAX(3*ncd, npd)) DEALLOCATE (inst_%ws_cbuf)
         END IF
         IF (.NOT. ALLOCATED(inst_%ws_cbuf)) ALLOCATE (inst_%ws_cbuf(MAX(3*ncd, npd)))
         CALL CD_Get_System_CoupledMotion( inst_%agg%sys%fast%system, inst_%ws_cbuf(1:ncd), &
                                           inst_%ws_cbuf(ncd + 1:2*ncd), inst_%ws_cbuf(2*ncd + 1:3*ncd), es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         states(koff + 1:koff + 3*ncd) = REAL(inst_%ws_cbuf(1:3*ncd), DbKi)
         koff = koff + 3*ncd
         CALL CD_Get_System_DynamicPoint_States( inst_%agg%sys%fast%system, inst_%ws_cbuf(1:npd), es, em )
         IF (es /= CD_SYSTEM_OK) THEN
            CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
         END IF
         states(koff + 1:koff + npd) = REAL(inst_%ws_cbuf(1:npd), DbKi)
         koff = koff + npd
         ! FAILURE fired flags at the mirror tail (0/1 per row)
         nfl = CD_AGG_NFailures( inst_%agg )
         IF (nfl > 0) THEN
            IF (koff + nfl > SIZE(states)) THEN
               CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the failure-flag mirror extends '// &
                                'past the state vector.', ErrStat, ErrMsg, RoutineName ); RETURN
            END IF
            BLOCK
               LOGICAL :: flg(nfl)
               CALL CD_AGG_Get_Failure_Flags( inst_%agg, flg, es, em )
               IF (es /= CD_AGG_OK) THEN
                  CALL SetErrStat( ErrID_Fatal, TRIM(em), ErrStat, ErrMsg, RoutineName ); RETURN
               END IF
               states(koff + 1:koff + nfl) = MERGE(1.0_DbKi, 0.0_DbKi, flg)
            END BLOCK
         END IF
      END IF
   END SUBROUTINE pack_state_mirror_ary

   SUBROUTINE build_farm_meshes( inst_e, u, y, InitInp, nT, ncp, pos, ErrStat, ErrMsg )
      !! FAST.Farm mesh construction: one CoupledKinematics/CoupledLoads pair per
      !! turbine over the flat aggregate moving surface. Records the turbine map,
      !! each slot's node index within its turbine's mesh, the farm reference
      !! offsets, and the per-turbine PtfmInit. Mesh reference positions are the
      !! UNDISPLACED TURBINE-LOCAL deck coordinates (recovered by the exact inverse
      !! rigid transform from the displaced farm-global aggregate positions), with
      !! the initial displacement carried in TranslationDisp -- MoorDyn's convention.
      TYPE(CD_OF_Instance),     INTENT(INOUT) :: inst_e
      TYPE(CD_InputType),       INTENT(INOUT) :: u
      TYPE(CD_OutputType),      INTENT(INOUT) :: y
      TYPE(CD_InitInputType),   INTENT(IN   ) :: InitInp
      INTEGER(IntKi),           INTENT(IN   ) :: nT, ncp
      REAL(wp),                 INTENT(IN   ) :: pos(:,:)
      INTEGER(IntKi),           INTENT(  OUT) :: ErrStat
      CHARACTER(*),             INTENT(  OUT) :: ErrMsg

      character(*), parameter   :: RoutineName = 'CableDyn:build_farm_meshes'
      INTEGER(IntKi)            :: es, s, J, k, nk, ErrStat2
      CHARACTER(512)            :: em
      CHARACTER(ErrMsgLen)      :: ErrMsg2
      REAL(R8Ki)                :: OrMatJ(3,3), OrientIdent(3,3), refp(3), ptfmJ(6), ploc(3)
      LOGICAL                   :: hasp

      ErrStat = ErrID_None
      ErrMsg  = ''
      inst_e%nT = nT
      ALLOCATE (inst_e%tmap(ncp), inst_e%node_of_slot(ncp), inst_e%farm_ref(3, nT), &
                inst_e%farm_ptfm(6, nT), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate the farm maps.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      CALL CD_AGG_TurbineOfMoving( inst_e%agg, inst_e%tmap, es, em )
      IF (es /= CD_AGG_OK) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: '//TRIM(em), ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      IF (ANY(inst_e%tmap < 1) .OR. ANY(inst_e%tmap > nT)) THEN
         CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: the farm turbine map carries a slot '// &
                          'outside 1..FarmSize.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      inst_e%farm_ref = REAL(InitInp%TurbineRefPos, R8Ki)
      inst_e%farm_ptfm = 0.0_R8Ki
      hasp = ALLOCATED(InitInp%PtfmInit)
      IF (hasp) THEN
         IF (SIZE(InitInp%PtfmInit, 1) >= 6 .AND. SIZE(InitInp%PtfmInit, 2) >= nT) THEN
            inst_e%farm_ptfm = REAL(InitInp%PtfmInit(1:6, 1:nT), R8Ki)
         END IF
      END IF
      DO J = 1, nT
         k = 0
         DO s = 1, ncp
            IF (inst_e%tmap(s) == J) THEN
               k = k + 1
               inst_e%node_of_slot(s) = k
            END IF
         END DO
      END DO

      ALLOCATE (u%CoupledKinematics(nT), y%CoupledLoads(nT), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate the farm mesh containers.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF
      CALL Eye( OrientIdent, ErrStat2, ErrMsg2 )
      DO J = 1, nT
         nk = COUNT(inst_e%tmap == J)
         CALL MeshCreate( BlankMesh=u%CoupledKinematics(J), IOS=COMPONENT_INPUT, Nnodes=MAX(nk, 1), &
                          TranslationDisp=.TRUE., Orientation=.TRUE., TranslationVel=.TRUE., &
                          RotationVel=.TRUE., TranslationAcc=.TRUE., RotationAcc=.TRUE., &
                          ErrStat=ErrStat2, ErrMess=ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         ptfmJ = inst_e%farm_ptfm(:, J)
         OrMatJ = EulerConstructZYX( ptfmJ(4:6) )
         DO s = 1, ncp
            IF (inst_e%tmap(s) /= J) CYCLE
            ! farm-global aggregate position -> turbine-local displaced -> undisplaced
            ploc = REAL(pos(:, s), R8Ki) - inst_e%farm_ref(:, J)
            refp = MATMUL( OrMatJ, ploc - ptfmJ(1:3) )
            CALL MeshPositionNode( u%CoupledKinematics(J), inst_e%node_of_slot(s), REAL(refp, ReKi), &
                                   ErrStat2, ErrMsg2, Orient=OrientIdent )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL MeshConstructElement( u%CoupledKinematics(J), ELEMENT_POINT, ErrStat2, ErrMsg2, &
                                       inst_e%node_of_slot(s) )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
         IF (nk == 0) THEN
            ! a turbine with no moored objects still needs a committed 1-node mesh so the
            ! farm glue's per-turbine mappings stay well-formed; it carries zero loads
            CALL MeshPositionNode( u%CoupledKinematics(J), 1, REAL(inst_e%farm_ref(:, J), ReKi), &
                                   ErrStat2, ErrMsg2, Orient=OrientIdent )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL MeshConstructElement( u%CoupledKinematics(J), ELEMENT_POINT, ErrStat2, ErrMsg2, 1 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF
         CALL MeshCommit( u%CoupledKinematics(J), ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         u%CoupledKinematics(J)%TranslationDisp = 0.0_ReKi
         u%CoupledKinematics(J)%TranslationVel  = 0.0_ReKi
         u%CoupledKinematics(J)%RotationVel     = 0.0_ReKi
         u%CoupledKinematics(J)%TranslationAcc  = 0.0_ReKi
         u%CoupledKinematics(J)%RotationAcc     = 0.0_ReKi
         DO k = 1, MAX(nk, 1)
            u%CoupledKinematics(J)%Orientation(:, :, k) = REAL(OrMatJ, ReKi)
         END DO
         DO s = 1, ncp
            IF (inst_e%tmap(s) /= J) CYCLE
            u%CoupledKinematics(J)%TranslationDisp(:, inst_e%node_of_slot(s)) = &
               REAL(REAL(pos(:, s), R8Ki) - inst_e%farm_ref(:, J), ReKi) - &
               u%CoupledKinematics(J)%Position(:, inst_e%node_of_slot(s))
         END DO
         CALL MeshCopy( SrcMesh=u%CoupledKinematics(J), DestMesh=y%CoupledLoads(J), CtrlCode=MESH_SIBLING, &
                        IOS=COMPONENT_OUTPUT, Force=.TRUE., Moment=.TRUE., ErrStat=ErrStat2, ErrMess=ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         y%CoupledLoads(J)%Force  = 0.0_ReKi
         y%CoupledLoads(J)%Moment = 0.0_ReKi
      END DO
   END SUBROUTINE build_farm_meshes

   SUBROUTINE farm_gather_kinematics(inst_e, u, ErrStat, ErrMsg)
      !! FAST.Farm gather: every turbine's mesh into the flat aggregate arrays,
      !! farm-global (Position + TranslationDisp + TurbineRefPos -- MoorDyn's exact
      !! boundary convention). Fills the instance ws_pos/ws_vel/ws_acc.
      TYPE(CD_OF_Instance),  INTENT(INOUT) :: inst_e
      TYPE(CD_InputType),    INTENT(IN)    :: u
      INTEGER(IntKi),        INTENT(OUT)   :: ErrStat
      CHARACTER(*),          INTENT(OUT)   :: ErrMsg

      INTEGER(IntKi) :: s, J, k

      ErrStat = ErrID_None
      ErrMsg  = ''
      IF (SIZE(u%CoupledKinematics) < inst_e%nT) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'CableDyn:farm_gather: the input carries fewer meshes than FarmSize.'
         RETURN
      END IF
      DO s = 1, inst_e%ncp
         J = inst_e%tmap(s)
         k = inst_e%node_of_slot(s)
         inst_e%ws_pos(:, s) = REAL(u%CoupledKinematics(J)%Position(:, k), wp) + &
                               REAL(u%CoupledKinematics(J)%TranslationDisp(:, k), wp) + &
                               REAL(inst_e%farm_ref(:, J), wp)
         inst_e%ws_vel(:, s) = REAL(u%CoupledKinematics(J)%TranslationVel(:, k), wp)
         inst_e%ws_acc(:, s) = REAL(u%CoupledKinematics(J)%TranslationAcc(:, k), wp)
         inst_e%ws_orient(:, :, s) = REAL(u%CoupledKinematics(J)%Orientation(:, :, k), wp)
         inst_e%ws_omega(:, s) = REAL(u%CoupledKinematics(J)%RotationVel(:, k), wp)
         inst_e%ws_alpha(:, s) = REAL(u%CoupledKinematics(J)%RotationAcc(:, k), wp)
      END DO
   END SUBROUTINE farm_gather_kinematics

   !> Extract absolute coupled-point kinematics from the mesh: position = reference
   !> Position + TranslationDisp (single-turbine: TurbineRefPos = 0). GROW-ONLY on the
   !> passed workspace: the marching path (UpdateStates/CalcOutput, every glue step of a
   !> long run) hands in the instance's persistent ws_* arrays, so the allocation happens
   !> once and the per-step cost is the copy alone.
   SUBROUTINE mesh_to_arrays(mesh, ncp, pos, vel, acc, orient, omega, alpha, ErrStat, ErrMsg)
      TYPE(MeshType),        INTENT(IN)    :: mesh
      INTEGER(IntKi),        INTENT(IN)    :: ncp
      REAL(wp), ALLOCATABLE, INTENT(INOUT) :: pos(:,:), vel(:,:), acc(:,:)
      REAL(wp), ALLOCATABLE, INTENT(INOUT) :: orient(:,:,:)
      REAL(wp), ALLOCATABLE, INTENT(INOUT) :: omega(:,:), alpha(:,:)
      INTEGER(IntKi),        INTENT(OUT)   :: ErrStat
      CHARACTER(*),          INTENT(OUT)   :: ErrMsg

      INTEGER(IntKi) :: i, ErrStat2

      ErrStat = ErrID_None
      ErrMsg  = ''
      IF (ALLOCATED(pos)) THEN
         IF (SIZE(pos, 2) /= ncp) DEALLOCATE (pos, vel, acc, orient, omega, alpha)
      END IF
      IF (.NOT. ALLOCATED(pos)) THEN
         ALLOCATE (pos(3, ncp), vel(3, ncp), acc(3, ncp), orient(3, 3, ncp), &
                   omega(3, ncp), alpha(3, ncp), STAT=ErrStat2)
         IF (ErrStat2 /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = 'CableDyn:mesh_to_arrays: could not allocate kinematics arrays.'
            RETURN
         END IF
      END IF
      DO i = 1, ncp
         pos(:, i) = REAL(mesh%Position(:, i), wp) + REAL(mesh%TranslationDisp(:, i), wp)
         vel(:, i) = REAL(mesh%TranslationVel(:, i), wp)
         acc(:, i) = REAL(mesh%TranslationAcc(:, i), wp)
         orient(:, :, i) = REAL(mesh%Orientation(:, :, i), wp)
         omega(:, i) = REAL(mesh%RotationVel(:, i), wp)
         alpha(:, i) = REAL(mesh%RotationAcc(:, i), wp)
      END DO
   END SUBROUTINE mesh_to_arrays

   !> Release an instance slot: end its aggregate, destroy its interpolation input, free
   !> its workspaces, close its native output stream and mark it inactive (reusable).
   !> Status and message are optional so rebuild error paths can release silently.
   SUBROUTINE release_instance( k, ErrStat, ErrMsg )
      INTEGER(IntKi),           INTENT(IN)  :: k
      INTEGER(IntKi), OPTIONAL, INTENT(OUT) :: ErrStat
      CHARACTER(*),   OPTIONAL, INTENT(OUT) :: ErrMsg
      INTEGER(IntKi)       :: es, ErrStat2
      CHARACTER(512)       :: em
      CHARACTER(ErrMsgLen) :: ErrMsg2
      IF (PRESENT(ErrStat)) ErrStat = ErrID_None
      IF (PRESENT(ErrMsg)) ErrMsg = ''
      IF (.NOT. ALLOCATED(Inst)) RETURN
      IF (k < 1_IntKi .OR. k > SIZE(Inst)) RETURN
      IF (CD_AGG_IsInitialized( Inst(k)%agg )) THEN
         CALL CD_AGG_End( Inst(k)%agg, es, em )
         IF (es /= CD_AGG_OK) THEN
            IF (PRESENT(ErrStat)) ErrStat = ErrID_Severe
            IF (PRESENT(ErrMsg)) ErrMsg = 'CableDyn mooring module: '//TRIM(em)
         END IF
      END IF
      ! The persistent interpolation container owns a full mesh copy; destroy it
      ! through the registry (mesh components are not plain allocatables).
      IF (Inst(k)%u_interp_ready) THEN
         CALL CD_DestroyInput( Inst(k)%u_interp, ErrStat2, ErrMsg2 )
         Inst(k)%u_interp_ready = .FALSE.
      END IF
      IF (ALLOCATED(Inst(k)%fx_xyz)) DEALLOCATE (Inst(k)%fx_xyz, Inst(k)%fx_vel, Inst(k)%fx_acc, &
                                                 Inst(k)%fx_wl, Inst(k)%fx_pd)
      Inst(k)%use_extfluid = .FALSE.
      Inst(k)%nfluid = 0
      ! close the CableDyn-owned committed-step record (flushes the final rows)
      IF (Inst(k)%native_un > 0) THEN
         CLOSE (Inst(k)%native_un)
         Inst(k)%native_un = -1_IntKi
      END IF
      Inst(k)%native_root = ''
      Inst(k)%native_wrote_n = -123456789_IntKi
      Inst(k)%native_wrote_t = 0.0_DbKi
      Inst(k)%native_pending = .FALSE.
      Inst(k)%active = .FALSE.
   END SUBROUTINE release_instance

   !> TRUE when a SeaState wave field carries ambient kinematics for the mooring: waves,
   !> a declared current, or a dynamically populated current field (an MHK current
   !> driven by InflowWind carries no CurrMod declaration).
   LOGICAL FUNCTION wavefield_is_ambient( wf ) RESULT(ambient)
      TYPE(SeaSt_WaveFieldType), INTENT(IN) :: wf
      ambient = wf%WaveMod /= WaveMod_None .OR. wf%Current_InitInput%CurrMod /= 0_IntKi .OR. wf%hasCurrField
   END FUNCTION wavefield_is_ambient

   !> Write the packed deck copy to <MooringFile folder>/<output root name><suffix>, the
   !> file the checkpoint rebuild and the linearization probes re-parse. In the MooringFile
   !> folder, deck-relative side files (bathymetry, WaterKin, Syrope tables) resolve as in
   !> the original run. The caller removes it after use (delete_file).
   SUBROUTINE write_deck_copy( p, suffix, path, ErrStat, ErrMsg )
      TYPE(CD_ParameterType), INTENT(IN)  :: p
      CHARACTER(*),           INTENT(IN)  :: suffix
      CHARACTER(*),           INTENT(OUT) :: path
      INTEGER(IntKi),         INTENT(OUT) :: ErrStat
      CHARACTER(*),           INTENT(OUT) :: ErrMsg
      INTEGER(IntKi) :: un, i, cut
      cut = MAX(INDEX(p%RootName, '/', BACK=.TRUE.), INDEX(p%RootName, '\', BACK=.TRUE.))
      path = TRIM(p%PriPath)//TRIM(p%RootName(cut + 1:))//suffix
      CALL GetNewUnit( un )
      CALL OpenFOutFile( un, TRIM(path), ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN
      DO i = 1, p%DeckCopy%NumLines
         WRITE (un, '(A)') TRIM(p%DeckCopy%Lines(i))
      END DO
      CLOSE (un)
   END SUBROUTINE write_deck_copy

   !> Remove a temporary file (silently: a failed removal leaves only a stale copy).
   SUBROUTINE delete_file( path )
      CHARACTER(*), INTENT(IN) :: path
      INTEGER(IntKi) :: un, ios
      CALL GetNewUnit( un )
      OPEN (un, FILE=TRIM(path), STATUS='OLD', IOSTAT=ios)
      IF (ios == 0) CLOSE (un, STATUS='DELETE', IOSTAT=ios)
   END SUBROUTINE delete_file

   !> Reserve a slot in the module-level instance registry, left INACTIVE. CD_Init marks it
   !> active only once initialization fully succeeds; a failed init returns early on any of its
   !> many error paths, so leaving the slot inactive means the NEXT new_instance() reuses it
   !> (rather than accumulating half-built active slots). Reuse resets the slot, which also
   !> auto-deallocates any partial fmf left by a prior failed init on that slot.
   INTEGER(IntKi) FUNCTION new_instance() RESULT(id)
      TYPE(CD_OF_Instance), ALLOCATABLE :: tmp(:)
      TYPE(CD_OF_Instance) :: empty
      INTEGER(IntKi) :: i

      IF (.NOT. ALLOCATED(Inst)) ALLOCATE (Inst(0))
      DO i = 1, SIZE(Inst)
         IF (.NOT. Inst(i)%active) THEN
            ! Intrinsic assignment from a default-initialized object resets every
            ! component and deallocates partial failed-init storage.  Do not use a
            ! structure constructor here: several nested aggregate workspaces have
            ! PRIVATE components, which makes that constructor nonconforming and is
            ! correctly rejected by Intel Fortran even though gfortran accepts it.
            Inst(i) = empty
            id = i
            RETURN
         END IF
      END DO
      ALLOCATE (tmp(SIZE(Inst) + 1))
      tmp(1:SIZE(Inst)) = Inst
      CALL MOVE_ALLOC(tmp, Inst)
      id = SIZE(Inst)                     ! new slot is default-initialized (active = .FALSE.)
   END FUNCTION new_instance

   LOGICAL FUNCTION profile_environment_requested() RESULT(requested)
      CHARACTER(32) :: value
      INTEGER :: stat

      value = ''
      CALL GET_ENVIRONMENT_VARIABLE('CABLEDYN_PROFILE', value, STATUS=stat)
      requested = stat == 0 .AND. LEN_TRIM(value) > 0 .AND. TRIM(value) /= '0'
   END FUNCTION profile_environment_requested

   SUBROUTINE reset_instance_profile( id, enabled, init_seconds, reset_global )
      INTEGER(IntKi), INTENT(IN) :: id
      LOGICAL, INTENT(IN) :: enabled, reset_global
      REAL(wp), INTENT(IN) :: init_seconds

      Inst(id)%prof_enabled = enabled
      Inst(id)%prof_n_update = 0_IntKi
      Inst(id)%prof_n_solve = 0_IntKi
      Inst(id)%prof_n_calcout = 0_IntKi
      Inst(id)%prof_n_fresh = 0_IntKi
      Inst(id)%prof_n_wave = 0_IntKi
      Inst(id)%prof_t_init = init_seconds
      Inst(id)%prof_t_update = 0.0_wp
      Inst(id)%prof_t_wave = 0.0_wp
      Inst(id)%prof_t_step = 0.0_wp
      Inst(id)%prof_t_calcout = 0.0_wp
      Inst(id)%prof_t_loads = 0.0_wp
      Inst(id)%prof_t_native = 0.0_wp
      IF (reset_global) THEN
         CALL CD_HermiteCable_Dyn_Recovery_Reset()
         IF (enabled) THEN
            CALL CD_HermiteCable_Dyn_Reset_Profile()
         ELSE
            CALL CD_HermiteCable_Dyn_Disable_Profile()
         END IF
      END IF
   END SUBROUTINE reset_instance_profile

   SUBROUTINE copy_instance_profile( destination, source )
      INTEGER(IntKi), INTENT(IN) :: destination, source

      Inst(destination)%prof_enabled = Inst(source)%prof_enabled
      Inst(destination)%prof_n_update = Inst(source)%prof_n_update
      Inst(destination)%prof_n_solve = Inst(source)%prof_n_solve
      Inst(destination)%prof_n_calcout = Inst(source)%prof_n_calcout
      Inst(destination)%prof_n_fresh = Inst(source)%prof_n_fresh
      Inst(destination)%prof_n_wave = Inst(source)%prof_n_wave
      Inst(destination)%prof_t_init = Inst(source)%prof_t_init
      Inst(destination)%prof_t_update = Inst(source)%prof_t_update
      Inst(destination)%prof_t_wave = Inst(source)%prof_t_wave
      Inst(destination)%prof_t_step = Inst(source)%prof_t_step
      Inst(destination)%prof_t_calcout = Inst(source)%prof_t_calcout
      Inst(destination)%prof_t_loads = Inst(source)%prof_t_loads
      Inst(destination)%prof_t_native = Inst(source)%prof_t_native
   END SUBROUTINE copy_instance_profile

   !> Validate an instance handle (fails closed on 0 / stale handles, e.g. after a
   !> checkpoint restore in a fresh process, which does not rebuild the solver instance).
   LOGICAL FUNCTION valid_instance(id, ErrStat, ErrMsg, caller) RESULT(ok)
      INTEGER(IntKi), INTENT(IN)    :: id
      INTEGER(IntKi), INTENT(INOUT) :: ErrStat
      CHARACTER(*),   INTENT(INOUT) :: ErrMsg
      CHARACTER(*),   INTENT(IN)    :: caller

      ok = .FALSE.
      IF (ALLOCATED(Inst)) THEN
         IF (id >= 1 .AND. id <= SIZE(Inst)) THEN
            IF (Inst(id)%active) THEN
               ok = .TRUE.
               RETURN
            END IF
         END IF
      END IF
      CALL SetErrStat( ErrID_Fatal, 'CableDyn mooring module: no active solver instance for this module '// &
                       'data (restoring a checkpoint in a new process does not rebuild the solver; '// &
                       'restart the simulation from its input files).', ErrStat, ErrMsg, caller )
   END FUNCTION valid_instance

END MODULE CableDyn
