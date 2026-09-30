! File: src/CableDyn_OpenFAST_HermiteFMF.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_OpenFAST_HermiteFMF
  !! OpenFAST-style module-form lifecycle (Init / UpdateStates / CalcOutput / End) over the
  !! Hermite finite-EI DYNAMIC POWER CABLE (CableDyn_HermiteCableDynamic) -- the coupling
  !! boundary a floating-platform solver drives:
  !!
  !!   kinematics IN : the coupled endpoint's position / velocity / acceleration and the
  !!                   parent orientation / angular velocity / angular acceleration at
  !!                   t_{n+1} (the platform hang-off, prescribed through the
  !!                   generalised-alpha step's exact boundary machinery);
  !!   loads OUT     : the force the cable exerts ON the platform at that endpoint (the
  !!                   constraint reaction -[M a + f_int - f_ext] at the coupled node's
  !!                   translational DOFs, evaluated at the committed state).
  !!
  !! This is the finite-EI counterpart of the EI=0 mooring shell (CableDyn_OpenFAST_FMF):
  !! the same lifecycle contract, carried by the bending-cable element the lazy-wave
  !! geometry needs. Solver controls (Newton budget, tolerance) freeze at Init, mirroring
  !! the MoorDyn-F parameter pattern; the optional Morison drag / added-mass / wave
  !! configuration passes through the underlying model's setters before the first step.
  !! One line per module instance (the multi-line registry aggregation is the OpenFAST
  !! integration layer's job, out of scope here).
  USE CableDyn_Precision, ONLY: wp, CD_ZERO, CD_ONE, CD_All_Finite, CD_Is_Finite
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Set_Drag, CD_HermiteCable_Dyn_Set_AddedMass, &
                                          CD_HermiteCable_Dyn_Set_Axial_Damping, &
                                          CD_HermiteCable_Dyn_Set_Waves, CD_HermiteCable_Dyn_Set_Held_Fluid, &
                                          CD_HermiteCable_Dyn_Set_Contact, CD_HermiteCable_Dyn_Set_Attachments, &
                                          CD_HermiteCable_Dyn_Set_EndConnection, &
                                          CD_HermiteCable_Dyn_Recompute_Acceleration, &
                                          CD_HermiteCable_Dyn_Set_ModifiedNewton, &
                                          CD_HermiteCable_Dyn_Set_Tensile_Safety, &
                                          CD_HermiteCable_Dyn_Set_Tensile_Monitor, &
                                          CD_HermiteCable_Dyn_Get_Tensile_Diagnostics, &
                                          CD_HermiteCable_Dyn_Get_Recovery_Diagnostics, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_Step_Recovering, &
                                          CD_HermiteCable_Dyn_Reaction, &
                                          CD_HermiteCable_Dyn_EndConnection_Moment, &
                                          CD_HermiteCable_Dyn_Curvature, &
                                          CD_HermiteCable_Dyn_Snapshot, CD_HermiteCable_Dyn_Restore, &
                                          CD_HermiteCable_Dyn_End, CD_HCDYN_OK
  USE CableDyn_HermiteCable, ONLY: CD_HermiteCable_Shapes
  USE CableDyn_EndConnection, ONLY: CD_EndConn_Project, CD_ENDCONN_OK, CD_ENDCONN_BADINPUT, &
                                    CD_ENDCONN_PINNED, CD_ENDCONN_FINITE, CD_ENDCONN_RIGID
  USE CableDyn_Bathymetry, ONLY: CD_BathymetryType
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CD_HFMF_ModuleType
  PUBLIC :: CD_HFMF_Init
  PUBLIC :: CD_HFMF_Set_Drag
  PUBLIC :: CD_HFMF_Set_Axial_Damping
  PUBLIC :: CD_HFMF_Set_ModifiedNewton
  PUBLIC :: CD_HFMF_Set_Tensile_Safety
  PUBLIC :: CD_HFMF_Set_Tensile_Monitor
  PUBLIC :: CD_HFMF_Get_Tensile_Diagnostics
  PUBLIC :: CD_HFMF_Get_Recovery_Diagnostics
  PUBLIC :: CD_HFMF_Set_Recovery_Max_Substeps
  PUBLIC :: CD_HFMF_Set_AddedMass
  PUBLIC :: CD_HFMF_Set_Waves
  PUBLIC :: CD_HFMF_Set_Held_Fluid
  PUBLIC :: CD_HFMF_Set_Contact
  PUBLIC :: CD_HFMF_Set_Attachments
  PUBLIC :: CD_HFMF_Set_EndConnection
  PUBLIC :: CD_HFMF_Refresh_Acceleration
  PUBLIC :: CD_HFMF_NNodes
  PUBLIC :: CD_HFMF_GetNodePositions
  PUBLIC :: CD_HFMF_UpdateStates
  PUBLIC :: CD_HFMF_Snapshot
  PUBLIC :: CD_HFMF_Restore
  PUBLIC :: CD_HFMF_PreflightCoupledKinematics
  PUBLIC :: CD_HFMF_SetCoupledKinematics
  PUBLIC :: CD_HFMF_GetCoupledKinematics
  PUBLIC :: CD_HFMF_CalcOutput
  PUBLIC :: CD_HFMF_Curvature
  PUBLIC :: CD_HFMF_MinSpanZ
  PUBLIC :: CD_HFMF_MirrorSize
  PUBLIC :: CD_HFMF_FrictionMirrorSize
  PUBLIC :: CD_HFMF_ForceMirrorSize
  PUBLIC :: CD_HFMF_PackMirror
  PUBLIC :: CD_HFMF_UnpackMirror
  PUBLIC :: CD_HFMF_End
  INTEGER, PARAMETER, PUBLIC :: CD_HFMF_OK = 0, CD_HFMF_BADINPUT = 1, CD_HFMF_SOLVEFAIL = 2

  TYPE :: CD_HFMF_ModuleType
    !! Module-form owner of one Hermite dynamic power cable behind the FMF lifecycle.
    TYPE(CD_HermiteCableDynType) :: line
    INTEGER :: coupled_node = 0            ! node whose translations the platform drives
    INTEGER :: pdof(3) = 0                 ! that node's global r DOFs
    REAL(wp) :: dt = CD_ZERO               ! coupling step, frozen at Init
    INTEGER :: max_iter = 100              ! Newton budget per step, frozen at Init
    REAL(wp) :: tol = 1.0e-4_wp            ! Newton tolerance, frozen at Init
    INTEGER :: recovery_max_substeps = 1024 ! internal recovery cap; nominal coupling dt is unchanged
    ! Cable-frame heading (cos, sin of the rotation about +z from the CALLER's global
    ! frame to the cable's internal solve frame). The initial static builder uses a frame
    ! whose x-z plane contains the cable chord. Connection-enabled deck models leave the
    ! out-of-plane dynamic DOFs free so a rotating parent can produce a 3D response.
    ! An azimuthal cable stores its heading here and the FMF boundary rotates
    ! kinematics IN (global -> local) and loads OUT (local -> global), so every caller
    ! works purely in global coordinates. Identity (1, 0) leaves all paths bit-for-bit.
    REAL(wp) :: frame_c = 1.0_wp, frame_s = 0.0_wp
    ! A finite or rigid connection at the coupled endpoint is attached to the host, so its
    ! no-moment direction is expressed in parent axes and rotates with the input
    ! mesh orientation. The DCM convention is OpenFAST's global-to-parent matrix.
    LOGICAL :: has_parent_endconn = .FALSE.
    INTEGER :: parent_endconn_index = 0
    REAL(wp) :: parent_endconn_d0(3) = CD_ZERO
    REAL(wp) :: parent_dcm(3, 3) = RESHAPE([CD_ONE, CD_ZERO, CD_ZERO, &
                                            CD_ZERO, CD_ONE, CD_ZERO, &
                                            CD_ZERO, CD_ZERO, CD_ONE], [3, 3])
    ! OpenFAST supplies mesh angular velocity and acceleration in global axes.
    ! Retain them with the parent orientation so an exact rigid connection can
    ! transport its tangent q/v/a consistently through output probes, rewinds,
    ! and the next implicit step.
    REAL(wp) :: parent_omega(3) = CD_ZERO
    REAL(wp) :: parent_alpha(3) = CD_ZERO
    REAL(wp) :: snap_parent_dcm(3, 3) = RESHAPE([CD_ONE, CD_ZERO, CD_ZERO, &
                                                 CD_ZERO, CD_ONE, CD_ZERO, &
                                                 CD_ZERO, CD_ZERO, CD_ONE], [3, 3])
    REAL(wp) :: snap_parent_omega(3) = CD_ZERO
    REAL(wp) :: snap_parent_alpha(3) = CD_ZERO
    LOGICAL :: initialized = .FALSE.
  END TYPE CD_HFMF_ModuleType

CONTAINS

  SUBROUTINE CD_HFMF_Init(self, l0, EA, EI, rho_a, w, q_static, fixed_dofs, seabed_z, kn, &
                          rho_inf, dt, coupled_node, ErrStat, ErrMsg, max_iter, tol, frame_cs, &
                          axial_quadrature_order, bending_quadrature_order)
    !! Initialise the module: build the Hermite dynamic model on the static-equilibrium
    !! seed and bind the coupled endpoint. The coupled node's THREE translational DOFs
    !! must be in fixed_dofs (the platform drives them; the module rejects a free coupled
    !! endpoint rather than silently pinning it). max_iter / tol default to the production
    !! step controls (100, 1e-4) and freeze at Init. frame_cs = (cos, sin) of the cable's
    !! frame heading about +z (see CD_HFMF_ModuleType): q_static and fixed_dofs are given
    !! in the LOCAL solve frame; the boundary rotates. Defaults to identity. Initialisation
    !! builds into a local CANDIDATE and replaces self only on success (the EI=0 FMF
    !! wrapper's pattern), so a failed re-initialisation preserves an already-working
    !! module instead of clobbering it on entry.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: l0(:), EA(:), EI(:), rho_a(:), w(:), q_static(:)
    INTEGER, INTENT(IN) :: fixed_dofs(:)
    REAL(wp), INTENT(IN) :: seabed_z, kn, rho_inf, dt
    INTEGER, INTENT(IN) :: coupled_node
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: max_iter
    REAL(wp), INTENT(IN), OPTIONAL :: tol
    REAL(wp), INTENT(IN), OPTIONAL :: frame_cs(2)
    INTEGER, INTENT(IN), OPTIONAL :: axial_quadrature_order, bending_quadrature_order
    TYPE(CD_HFMF_ModuleType) :: candidate
    INTEGER :: es, s, gd
    CHARACTER(300) :: em

    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. CD_Is_Finite(dt) .OR. dt <= CD_ZERO) THEN
      CALL fail('dt must be finite and positive'); RETURN
    END IF
    IF (coupled_node < 1 .OR. coupled_node > SIZE(l0) + 1) THEN
      CALL fail('coupled_node out of range'); RETURN
    END IF
    DO s = 1, 3
      gd = 6*(coupled_node - 1) + s
      IF (.NOT. ANY(fixed_dofs == gd)) THEN
        CALL fail('the coupled node''s three translational DOFs must be in fixed_dofs'); RETURN
      END IF
    END DO
    IF (PRESENT(max_iter)) THEN
      IF (max_iter < 1) THEN
        CALL fail('max_iter must be >= 1'); RETURN
      END IF
      candidate%max_iter = max_iter
    END IF
    IF (PRESENT(tol)) THEN
      IF (.NOT. CD_Is_Finite(tol) .OR. tol <= CD_ZERO) THEN
        CALL fail('tol must be finite and positive'); RETURN
      END IF
      candidate%tol = tol
    END IF
    IF (PRESENT(frame_cs)) THEN
      IF (.NOT. CD_All_Finite(frame_cs)) THEN
        CALL fail('frame_cs must be finite'); RETURN
      END IF
      IF (ABS(frame_cs(1)**2 + frame_cs(2)**2 - 1.0_wp) > 1.0e-9_wp) THEN
        CALL fail('frame_cs must be a unit heading (cos^2 + sin^2 = 1)'); RETURN
      END IF
      candidate%frame_c = frame_cs(1)
      candidate%frame_s = frame_cs(2)
    END IF
    CALL CD_HermiteCable_Dyn_Init(candidate%line, l0, EA, EI, rho_a, w, q_static, fixed_dofs, &
                                  seabed_z, kn, rho_inf, es, em, axial_quadrature_order, &
                                  bending_quadrature_order)
    IF (es /= CD_HCDYN_OK) THEN
      CALL fail(TRIM(em)); RETURN
    END IF
    candidate%coupled_node = coupled_node
    DO s = 1, 3
      candidate%pdof(s) = 6*(coupled_node - 1) + s
    END DO
    candidate%dt = dt
    candidate%initialized = .TRUE.
    ! success: release any existing module, install the candidate, release its copy
    CALL CD_HFMF_End(self)
    self = candidate
    CALL CD_HFMF_End(candidate)

  CONTAINS
    SUBROUTINE fail(msg)
      CHARACTER(*), INTENT(IN) :: msg
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_Init: '//msg
    END SUBROUTINE fail
  END SUBROUTINE CD_HFMF_Init

  SUBROUTINE CD_HFMF_Set_Drag(self, rho_w, diam, cdn, cdt, waterline_z, current, ErrStat, ErrMsg, gravity)
    !! Pass-through to the model's Morison drag configuration (call between Init and the
    !! first UpdateStates). The configuration also sets the free surface above which the
    !! line regains its displaced-water buoyancy; gravity (optional, default 9.80665 m/s^2)
    !! must be the value the submerged weight was formed with.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: rho_w, diam(:), cdn(:), cdt(:), waterline_z, current(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: gravity
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Drag: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Set_Drag(self%line, rho_w, diam, cdn, cdt, waterline_z, current, es, em, &
                                      gravity=gravity)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Drag: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_Drag

  SUBROUTINE CD_HFMF_Set_Attachments(self, node, mass, volume, cda, cdax, ca, rho_w, gravity, ErrStat, ErrMsg)
    !! Install discrete point attachments lumped at cable nodes (internal node order,
    !! anchor = node 1); see CD_HermiteCable_Dyn_Set_Attachments.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: node(:)
    REAL(wp), INTENT(IN) :: mass(:), volume(:), cda(:), cdax(:), ca(:), rho_w, gravity
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Attachments: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Set_Attachments(self%line, node, mass, volume, cda, cdax, ca, rho_w, gravity, &
                                             es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_Attachments

  SUBROUTINE CD_HFMF_Set_Axial_Damping(self, BA, ErrStat, ErrMsg)
    !! Pass-through for resolved per-element axial Kelvin-Voigt coefficients [N s].
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: BA(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK
    ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_Set_Axial_Damping: module not initialised'
      RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Set_Axial_Damping(self%line, BA, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_Set_Axial_Damping: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_Axial_Damping

  SUBROUTINE CD_HFMF_Set_EndConnection(self, k_rot, d0_global, ErrStat, ErrMsg, &
                                       coupled_d0_parent, parent_orientation, connection_mode)
    !! Install pinned, finite-stiffness, or exact rigid end connections. End order follows the internal
    !! anchor-to-fairlead model: column 1 is node 1 and column 2 is the final node.
    !! Preferred directions arrive in the caller's global frame and are rotated
    !! into the cable solve frame here.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: k_rot(2), d0_global(3, 2)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: coupled_d0_parent(3), parent_orientation(3, 3)
    INTEGER, INTENT(IN), OPTIONAL :: connection_mode(2)
    INTEGER :: es, iend, coupled_iend, mode_eff(2)
    LOGICAL :: endpoint_connected
    REAL(wp) :: d0_local(3, 2)
    REAL(wp) :: parent_d0(3), parent_d0_norm, dcm(3, 3)
    CHARACTER(300) :: em

    ErrStat = CD_HFMF_OK
    ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_Set_EndConnection: module not initialised'
      RETURN
    END IF
    IF (PRESENT(coupled_d0_parent) .NEQV. PRESENT(parent_orientation)) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_Set_EndConnection: coupled_d0_parent and parent_orientation must be supplied together'
      RETURN
    END IF
    coupled_iend = 0
    IF (self%coupled_node == 1) coupled_iend = 1
    IF (self%coupled_node == self%line%nn) coupled_iend = 2
    dcm = identity_dcm()
    parent_d0 = CD_ZERO
    IF (PRESENT(connection_mode)) THEN
      mode_eff = connection_mode
    ELSE
      WHERE (k_rot > CD_ZERO)
        mode_eff = CD_ENDCONN_FINITE
      ELSEWHERE
        mode_eff = CD_ENDCONN_PINNED
      END WHERE
    END IF
    IF (PRESENT(coupled_d0_parent)) THEN
      ! Nested tests: Fortran does not short-circuit .OR., so mode_eff(0) must not be read.
      endpoint_connected = coupled_iend /= 0
      IF (endpoint_connected) endpoint_connected = mode_eff(coupled_iend) /= CD_ENDCONN_PINNED
      IF (.NOT. endpoint_connected) THEN
        ErrStat = CD_HFMF_BADINPUT
        ErrMsg = 'CD_HFMF_Set_EndConnection: parent-relative direction requires a connected coupled endpoint'
        RETURN
      END IF
      IF (.NOT. valid_parent_dcm(parent_orientation)) THEN
        ErrStat = CD_HFMF_BADINPUT
        ErrMsg = 'CD_HFMF_Set_EndConnection: parent_orientation must be a finite proper orthogonal DCM'
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(coupled_d0_parent)) THEN
        ErrStat = CD_HFMF_BADINPUT
        ErrMsg = 'CD_HFMF_Set_EndConnection: coupled parent direction must be finite'
        RETURN
      END IF
      parent_d0_norm = SQRT(DOT_PRODUCT(coupled_d0_parent, coupled_d0_parent))
      IF (.NOT. CD_Is_Finite(parent_d0_norm) .OR. parent_d0_norm <= SQRT(TINY(CD_ONE))) THEN
        ErrStat = CD_HFMF_BADINPUT
        ErrMsg = 'CD_HFMF_Set_EndConnection: coupled parent direction must be non-null'
        RETURN
      END IF
      dcm = parent_orientation
      parent_d0 = coupled_d0_parent/parent_d0_norm
    END IF
    DO iend = 1, 2
      d0_local(:, iend) = frame_to_local(self, d0_global(:, iend))
    END DO
    IF (PRESENT(coupled_d0_parent)) THEN
      d0_local(:, coupled_iend) = frame_to_local(self, MATMUL(TRANSPOSE(dcm), parent_d0))
    END IF
    CALL CD_HermiteCable_Dyn_Set_EndConnection(self%line, k_rot, d0_local, es, em, &
                                               connection_mode=mode_eff)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_Set_EndConnection: '//TRIM(em)
      RETURN
    END IF
    self%has_parent_endconn = PRESENT(coupled_d0_parent)
    self%parent_endconn_index = 0
    self%parent_endconn_d0 = CD_ZERO
    self%parent_dcm = identity_dcm()
    self%parent_omega = CD_ZERO
    self%parent_alpha = CD_ZERO
    self%snap_parent_dcm = self%parent_dcm
    self%snap_parent_omega = CD_ZERO
    self%snap_parent_alpha = CD_ZERO
    IF (PRESENT(coupled_d0_parent)) THEN
      self%parent_endconn_index = coupled_iend
      self%parent_endconn_d0 = parent_d0
      self%parent_dcm = dcm
      self%snap_parent_dcm = dcm
    END IF
  END SUBROUTINE CD_HFMF_Set_EndConnection

  SUBROUTINE CD_HFMF_Set_ModifiedNewton(self, enabled, ErrStat, ErrMsg)
    !! Pass-through to the model's modified-Newton (within-step tangent reuse) switch:
    !! one effective-tangent assembly + factorization per dynamic step, contraction-gated
    !! with a stale-direction refresh (see CD_HermiteCable_Dyn_Set_ModifiedNewton).
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    LOGICAL, INTENT(IN) :: enabled
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_ModifiedNewton: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Set_ModifiedNewton(self%line, enabled, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_ModifiedNewton: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_ModifiedNewton

  SUBROUTINE CD_HFMF_Set_Tensile_Safety(self, enabled, ErrStat, ErrMsg)
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    LOGICAL, INTENT(IN) :: enabled
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(200) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Tensile_Safety: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Set_Tensile_Safety(self%line, enabled, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Tensile_Safety: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_Tensile_Safety

  SUBROUTINE CD_HFMF_Set_Tensile_Monitor(self, mode, ErrStat, ErrMsg, strain_tolerance)
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: mode
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: strain_tolerance
    INTEGER :: es
    CHARACTER(200) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Tensile_Monitor: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Set_Tensile_Monitor(self%line, mode, es, em, strain_tolerance)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Tensile_Monitor: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_Tensile_Monitor

  SUBROUTINE CD_HFMF_Get_Tensile_Diagnostics(self, event_count, worst_force, worst_threshold, &
                                             worst_element, worst_xi, worst_time, ErrStat, ErrMsg)
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    INTEGER, INTENT(OUT) :: event_count, worst_element, ErrStat
    REAL(wp), INTENT(OUT) :: worst_force, worst_threshold, worst_xi, worst_time
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(200) :: em
    IF (.NOT. self%initialized) THEN
      event_count = 0; worst_element = 0
      worst_force = CD_ZERO; worst_threshold = CD_ZERO; worst_xi = CD_ZERO; worst_time = CD_ZERO
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Get_Tensile_Diagnostics: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Get_Tensile_Diagnostics(self%line, event_count, worst_force, &
                                                     worst_threshold, worst_element, worst_xi, worst_time, es, em)
    IF (es == CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_OK; ErrMsg = ''
    ELSE
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Get_Tensile_Diagnostics: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Get_Tensile_Diagnostics

  SUBROUTINE CD_HFMF_Get_Recovery_Diagnostics(self, event_count, max_substeps_used, ErrStat, ErrMsg)
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    INTEGER, INTENT(OUT) :: event_count, max_substeps_used, ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(200) :: em
    IF (.NOT. self%initialized) THEN
      event_count = 0; max_substeps_used = 0
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Get_Recovery_Diagnostics: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Get_Recovery_Diagnostics(self%line, event_count, max_substeps_used, es, em)
    IF (es == CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_OK; ErrMsg = ''
    ELSE
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Get_Recovery_Diagnostics: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Get_Recovery_Diagnostics

  SUBROUTINE CD_HFMF_Set_Recovery_Max_Substeps(self, max_substeps, ErrStat, ErrMsg)
    !! Set the maximum internal subdivisions attempted after a failed or under-resolved nominal
    !! coupling step. This does not change the host time grid or the committed interval endpoint.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    INTEGER, INTENT(IN) :: max_substeps
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_Set_Recovery_Max_Substeps: module not initialised'; RETURN
    END IF
    IF (max_substeps < 4 .OR. max_substeps > 65536) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_Set_Recovery_Max_Substeps: max_substeps must be in [4,65536]'; RETURN
    END IF
    self%recovery_max_substeps = max_substeps
  END SUBROUTINE CD_HFMF_Set_Recovery_Max_Substeps

  SUBROUTINE CD_HFMF_Set_AddedMass(self, rho_w, diam, can, cat, waterline_z, ErrStat, ErrMsg)
    !! Pass-through to the model's Morison added-mass configuration.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: rho_w, diam(:), can(:), cat(:), waterline_z
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_AddedMass: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Set_AddedMass(self%line, rho_w, diam, can, cat, waterline_z, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_AddedMass: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_AddedMass

  SUBROUTINE CD_HFMF_Set_Waves(self, height, period, direction_deg, depth, gravity, ErrStat, ErrMsg)
    !! Pass-through to the model's regular (Airy) wave-field configuration.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: height, period, direction_deg, depth, gravity
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Waves: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Set_Waves(self%line, height, period, direction_deg, depth, gravity, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Waves: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_Waves

  SUBROUTINE CD_HFMF_UpdateStates(self, u_pos, u_vel, u_acc, ErrStat, ErrMsg, u_orientation, &
                                  u_angular_velocity, u_angular_acceleration)
    !! Advance the cable one coupling step of dt with the coupled endpoint's kinematics
    !! at t_{n+1} prescribed to (u_pos, u_vel, u_acc) -- the FMF kinematics-in boundary.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: u_pos(3), u_vel(3), u_acc(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: u_orientation(3, 3)
    REAL(wp), INTENT(IN), OPTIONAL :: u_angular_velocity(3), u_angular_acceleration(3)
    INTEGER :: es
    REAL(wp) :: target_dcm(3, 3), target_d0(3, 2), target_omega(3), target_alpha(3)
    REAL(wp) :: target_d0_rate(3, 2), target_d0_acceleration(3, 2), omega_local(3), alpha_local(3)
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_UpdateStates: module not initialised'; RETURN
    END IF
    IF (.NOT. (CD_All_Finite(u_pos) .AND. CD_All_Finite(u_vel) .AND. &
               CD_All_Finite(u_acc))) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_UpdateStates: inputs must be finite'; RETURN
    END IF
    IF (PRESENT(u_angular_velocity) .NEQV. PRESENT(u_angular_acceleration)) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_UpdateStates: angular velocity and acceleration must be supplied together'
      RETURN
    END IF
    target_omega = CD_ZERO
    target_alpha = CD_ZERO
    IF (PRESENT(u_angular_velocity)) THEN
      IF (.NOT. (CD_All_Finite(u_angular_velocity) .AND. &
                 CD_All_Finite(u_angular_acceleration))) THEN
        ErrStat = CD_HFMF_BADINPUT
        ErrMsg = 'CD_HFMF_UpdateStates: angular kinematics must be finite'
        RETURN
      END IF
      target_omega = u_angular_velocity
      target_alpha = u_angular_acceleration
    END IF
    target_dcm = self%parent_dcm
    IF (PRESENT(u_orientation)) THEN
      IF (.NOT. valid_parent_dcm(u_orientation)) THEN
        ErrStat = CD_HFMF_BADINPUT
        ErrMsg = 'CD_HFMF_UpdateStates: orientation must be a finite proper orthogonal DCM'
        RETURN
      END IF
      target_dcm = u_orientation
    END IF
    target_d0 = self%line%endconn_d0
    target_d0_rate = CD_ZERO
    target_d0_acceleration = CD_ZERO
    IF (self%has_parent_endconn) THEN
      target_d0(:, self%parent_endconn_index) = frame_to_local(self, &
                                                               MATMUL(TRANSPOSE(target_dcm), self%parent_endconn_d0))
      IF (self%line%endconn_mode(self%parent_endconn_index) == CD_ENDCONN_RIGID) THEN
        omega_local = frame_to_local(self, target_omega)
        alpha_local = frame_to_local(self, target_alpha)
        target_d0_rate(:, self%parent_endconn_index) = &
          cross3(omega_local, target_d0(:, self%parent_endconn_index))
        target_d0_acceleration(:, self%parent_endconn_index) = &
          cross3(alpha_local, target_d0(:, self%parent_endconn_index)) + &
          cross3(omega_local, target_d0_rate(:, self%parent_endconn_index))
      END IF
    END IF
    ! Advance the finite-EI line over the prescribed interval. Internal substeps are used when the
    ! nominal interval is too large for the boundary-motion or environmental-load increment, while
    ! preserving the requested endpoint time.
    IF (self%has_parent_endconn) THEN
      CALL CD_HermiteCable_Dyn_Step_Recovering(self%line, self%dt, self%max_iter, self%tol, es, em, &
                                               pres_dofs=self%pdof, pres_q=frame_to_local(self, u_pos), &
                                               pres_v=frame_to_local(self, u_vel), pres_a=frame_to_local(self, u_acc), &
                                               max_substeps=self%recovery_max_substeps, &
                                               endconn_direction=target_d0, &
                                               endconn_direction_rate=target_d0_rate, &
                                               endconn_direction_acceleration=target_d0_acceleration)
    ELSE
      CALL CD_HermiteCable_Dyn_Step_Recovering(self%line, self%dt, self%max_iter, self%tol, es, em, &
                                               pres_dofs=self%pdof, pres_q=frame_to_local(self, u_pos), &
                                               pres_v=frame_to_local(self, u_vel), pres_a=frame_to_local(self, u_acc), &
                                               max_substeps=self%recovery_max_substeps)
    END IF
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_SOLVEFAIL; ErrMsg = 'CD_HFMF_UpdateStates: '//TRIM(em)
    ELSE IF (self%has_parent_endconn) THEN
      self%parent_dcm = target_dcm
      self%parent_omega = target_omega
      self%parent_alpha = target_alpha
    END IF
  END SUBROUTINE CD_HFMF_UpdateStates

  SUBROUTINE CD_HFMF_Snapshot(self, ErrStat, ErrMsg)
    !! Capture the cable's committed step state (q, v, a, t) for the aggregate
    !! stage-then-commit contract (delegates to CD_HermiteCable_Dyn_Snapshot).
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Snapshot: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Snapshot(self%line, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Snapshot: '//TRIM(em)
    ELSE
      self%snap_parent_dcm = self%parent_dcm
      self%snap_parent_omega = self%parent_omega
      self%snap_parent_alpha = self%parent_alpha
    END IF
  END SUBROUTINE CD_HFMF_Snapshot

  SUBROUTINE CD_HFMF_Restore(self, ErrStat, ErrMsg)
    !! Restore the cable's committed step state from the last Snapshot (fails closed
    !! without a valid snapshot).
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Restore: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Restore(self%line, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Restore: '//TRIM(em)
    ELSE
      self%parent_dcm = self%snap_parent_dcm
      self%parent_omega = self%snap_parent_omega
      self%parent_alpha = self%snap_parent_alpha
    END IF
  END SUBROUTINE CD_HFMF_Restore

  SUBROUTINE CD_HFMF_PreflightCoupledKinematics(self, u_pos, u_vel, u_acc, ErrStat, ErrMsg, u_orientation, &
                                                u_angular_velocity, u_angular_acceleration)
    !! Validate every deterministic failure condition in SetCoupledKinematics without
    !! changing the cable. Aggregate callers use this before updating any sibling
    !! subsystem, so a bad later cable column cannot leave an earlier column mutated.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    REAL(wp), INTENT(IN) :: u_pos(3), u_vel(3), u_acc(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: u_orientation(3, 3)
    REAL(wp), INTENT(IN), OPTIONAL :: u_angular_velocity(3), u_angular_acceleration(3)
    REAL(wp) :: target_dcm(3, 3), d0_target(3, 2), tangent_q(3), tangent_v(3), tangent_a(3)
    REAL(wp) :: target_omega(3), target_alpha(3)
    CHARACTER(300) :: em

    CALL prepare_coupled_kinematics(self, u_pos, u_vel, u_acc, target_dcm, target_omega, target_alpha, &
                                    d0_target, tangent_q, tangent_v, tangent_a, ErrStat, em, &
                                    u_orientation, u_angular_velocity, u_angular_acceleration)
    IF (ErrStat /= CD_HFMF_OK) THEN
      ErrMsg = 'CD_HFMF_PreflightCoupledKinematics: '//TRIM(em)
    ELSE
      ErrMsg = ''
    END IF
  END SUBROUTINE CD_HFMF_PreflightCoupledKinematics

  SUBROUTINE CD_HFMF_SetCoupledKinematics(self, u_pos, u_vel, u_acc, ErrStat, ErrMsg, u_orientation, &
                                          u_angular_velocity, u_angular_acceleration)
    !! Direct-feedthrough transfer: write ONLY the coupled node's three translational (r, v, a)
    !! DOFs into the line state, WITHOUT stepping the interior or advancing simulation time. Used
    !! between an UpdateStates_Moving (the host moves the fairlead) and the following CalcOutput so
    !! the reported reaction reflects the CURRENT coupled kinematics against the FROZEN interior
    !! (the last committed step's shape) -- the instantaneous feedthrough load, rather than the load
    !! at the last committed fairlead. It does not re-solve, does not touch the interior DOFs, and
    !! does not advance self%line%t; the next CD_HFMF_UpdateStates re-prescribes the boundary and
    !! performs the actual implicit advance.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: u_pos(3), u_vel(3), u_acc(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: u_orientation(3, 3)
    REAL(wp), INTENT(IN), OPTIONAL :: u_angular_velocity(3), u_angular_acceleration(3)
    REAL(wp) :: target_dcm(3, 3), d0_target(3, 2)
    REAL(wp) :: tangent_q(3), tangent_v(3), tangent_a(3)
    REAL(wp) :: target_omega(3), target_alpha(3)
    INTEGER :: base
    CHARACTER(300) :: em

    CALL prepare_coupled_kinematics(self, u_pos, u_vel, u_acc, target_dcm, target_omega, target_alpha, &
                                    d0_target, tangent_q, tangent_v, tangent_a, ErrStat, em, &
                                    u_orientation, u_angular_velocity, u_angular_acceleration)
    IF (ErrStat /= CD_HFMF_OK) THEN
      ErrMsg = 'CD_HFMF_SetCoupledKinematics: '//TRIM(em)
      RETURN
    END IF
    IF (self%has_parent_endconn) THEN
      IF (self%line%endconn_mode(self%parent_endconn_index) == CD_ENDCONN_RIGID) THEN
        IF (self%parent_endconn_index == 1) THEN
          base = 3
        ELSE
          base = 6*(self%line%nn - 1) + 3
        END IF
        self%line%q(base + 1:base + 3) = tangent_q
        self%line%v(base + 1:base + 3) = tangent_v
        self%line%a(base + 1:base + 3) = tangent_a
      END IF
      self%line%endconn_d0(:, self%parent_endconn_index) = d0_target(:, self%parent_endconn_index)
      self%parent_dcm = target_dcm
      self%parent_omega = target_omega
      self%parent_alpha = target_alpha
    END IF
    ! self%pdof holds the coupled node's three translational DOFs in the cable's local
    ! frame; the boundary rotates the caller's global kinematics in.
    self%line%q(self%pdof) = frame_to_local(self, u_pos)
    self%line%v(self%pdof) = frame_to_local(self, u_vel)
    self%line%a(self%pdof) = frame_to_local(self, u_acc)
    ErrStat = CD_HFMF_OK
    ErrMsg = ''
  END SUBROUTINE CD_HFMF_SetCoupledKinematics

  SUBROUTINE prepare_coupled_kinematics(self, u_pos, u_vel, u_acc, target_dcm, target_omega, target_alpha, &
                                        d0_target, tangent_q, tangent_v, tangent_a, ErrStat, ErrMsg, &
                                        u_orientation, u_angular_velocity, u_angular_acceleration)
    !! Prepare the complete non-stepping boundary update without mutating self. Keep
    !! every fallible operation here so Preflight and Set have identical acceptance.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    REAL(wp), INTENT(IN) :: u_pos(3), u_vel(3), u_acc(3)
    REAL(wp), INTENT(OUT) :: target_dcm(3, 3), target_omega(3), target_alpha(3), d0_target(3, 2)
    REAL(wp), INTENT(OUT) :: tangent_q(3), tangent_v(3), tangent_a(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: u_orientation(3, 3)
    REAL(wp), INTENT(IN), OPTIONAL :: u_angular_velocity(3), u_angular_acceleration(3)
    REAL(wp) :: dn, magnitude, magnitude_rate, magnitude_acceleration
    REAL(wp) :: old_direction(3), old_direction_rate(3), target_direction(3)
    REAL(wp) :: target_direction_rate(3), target_direction_acceleration(3)
    REAL(wp) :: omega_local(3), alpha_local(3), rotation_check(3, 3)
    INTEGER :: base, ecstat
    CHARACTER(200) :: ecmsg

    ErrStat = CD_HFMF_OK
    ErrMsg = ''
    target_dcm = identity_dcm()
    d0_target = CD_ZERO
    tangent_q = CD_ZERO
    tangent_v = CD_ZERO
    tangent_a = CD_ZERO
    target_omega = CD_ZERO
    target_alpha = CD_ZERO
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'module not initialised'
      RETURN
    END IF
    IF (.NOT. (CD_All_Finite(u_pos) .AND. CD_All_Finite(u_vel) .AND. &
               CD_All_Finite(u_acc))) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'inputs must be finite'
      RETURN
    END IF
    IF (PRESENT(u_angular_velocity) .NEQV. PRESENT(u_angular_acceleration)) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'angular velocity and acceleration must be supplied together'
      RETURN
    END IF
    IF (PRESENT(u_angular_velocity)) THEN
      IF (.NOT. (CD_All_Finite(u_angular_velocity) .AND. &
                 CD_All_Finite(u_angular_acceleration))) THEN
        ErrStat = CD_HFMF_BADINPUT
        ErrMsg = 'angular kinematics must be finite'
        RETURN
      END IF
      target_omega = u_angular_velocity
      target_alpha = u_angular_acceleration
    END IF
    target_dcm = self%parent_dcm
    IF (PRESENT(u_orientation)) THEN
      IF (.NOT. valid_parent_dcm(u_orientation)) THEN
        ErrStat = CD_HFMF_BADINPUT
        ErrMsg = 'orientation must be a finite proper orthogonal DCM'
        RETURN
      END IF
      target_dcm = u_orientation
    END IF
    d0_target = self%line%endconn_d0
    IF (.NOT. self%has_parent_endconn) RETURN

    d0_target(:, self%parent_endconn_index) = frame_to_local(self, &
                                                             MATMUL(TRANSPOSE(target_dcm), self%parent_endconn_d0))
    dn = SQRT(DOT_PRODUCT(d0_target(:, self%parent_endconn_index), &
                          d0_target(:, self%parent_endconn_index)))
    IF (.NOT. CD_Is_Finite(dn) .OR. dn <= SQRT(TINY(CD_ONE))) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'rotated end-connection direction is singular'
      RETURN
    END IF
    d0_target(:, self%parent_endconn_index) = d0_target(:, self%parent_endconn_index)/dn
    IF (self%line%endconn_mode(self%parent_endconn_index) /= CD_ENDCONN_RIGID) RETURN

    IF (self%parent_endconn_index == 1) THEN
      base = 3
    ELSE
      base = 6*(self%line%nn - 1) + 3
    END IF
    CALL CD_EndConn_Project(self%line%q(base + 1:base + 3), &
                            d0_target(:, self%parent_endconn_index), tangent_q, ecstat, ecmsg)
    IF (ecstat /= CD_ENDCONN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'rigid tangent projection failed: '//TRIM(ecmsg)
      RETURN
    END IF
    CALL direction_rotation(self%line%endconn_d0(:, self%parent_endconn_index), &
                            d0_target(:, self%parent_endconn_index), rotation_check, ecstat, ecmsg)
    IF (ecstat /= CD_ENDCONN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'rigid direction update failed: '//TRIM(ecmsg)
      RETURN
    END IF
    old_direction = self%line%endconn_d0(:, self%parent_endconn_index)
    target_direction = d0_target(:, self%parent_endconn_index)
    magnitude = SQRT(DOT_PRODUCT(tangent_q, tangent_q))
    omega_local = frame_to_local(self, self%parent_omega)
    old_direction_rate = cross3(omega_local, old_direction)
    magnitude_rate = DOT_PRODUCT(self%line%v(base + 1:base + 3), old_direction)
    magnitude_acceleration = DOT_PRODUCT(self%line%a(base + 1:base + 3), old_direction) + &
                             magnitude*DOT_PRODUCT(old_direction_rate, old_direction_rate)
    omega_local = frame_to_local(self, target_omega)
    alpha_local = frame_to_local(self, target_alpha)
    target_direction_rate = cross3(omega_local, target_direction)
    target_direction_acceleration = cross3(alpha_local, target_direction) + &
                                    cross3(omega_local, target_direction_rate)
    tangent_v = magnitude_rate*target_direction + magnitude*target_direction_rate
    tangent_a = magnitude_acceleration*target_direction + &
                2.0_wp*magnitude_rate*target_direction_rate + magnitude*target_direction_acceleration
    IF (.NOT. (CD_All_Finite(tangent_v) .AND. CD_All_Finite(tangent_a))) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'rigid tangent transport produced non-finite kinematics'
    END IF
  END SUBROUTINE prepare_coupled_kinematics

  SUBROUTINE CD_HFMF_GetCoupledKinematics(self, pos, vel, acc, ErrStat, ErrMsg, orientation, &
                                          angular_velocity, angular_acceleration)
    !! Read the coupled node's committed kinematics in the CALLER's global frame (the
    !! inverse of the SetCoupledKinematics rotation). Encapsulates the frame so no
    !! consumer needs to poke self%line%q directly.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    REAL(wp), INTENT(OUT) :: pos(3), vel(3), acc(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(OUT), OPTIONAL :: orientation(3, 3)
    REAL(wp), INTENT(OUT), OPTIONAL :: angular_velocity(3), angular_acceleration(3)
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    pos = CD_ZERO; vel = CD_ZERO; acc = CD_ZERO
    IF (PRESENT(orientation)) orientation = identity_dcm()
    IF (PRESENT(angular_velocity)) angular_velocity = CD_ZERO
    IF (PRESENT(angular_acceleration)) angular_acceleration = CD_ZERO
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_GetCoupledKinematics: module not initialised'; RETURN
    END IF
    pos = frame_to_global(self, self%line%q(self%pdof))
    vel = frame_to_global(self, self%line%v(self%pdof))
    acc = frame_to_global(self, self%line%a(self%pdof))
    IF (PRESENT(orientation)) orientation = self%parent_dcm
    IF (PRESENT(angular_velocity)) angular_velocity = self%parent_omega
    IF (PRESENT(angular_acceleration)) angular_acceleration = self%parent_alpha
  END SUBROUTINE CD_HFMF_GetCoupledKinematics

  SUBROUTINE CD_HFMF_Set_Held_Fluid(self, fluid_velocity, fluid_acceleration, waterline_z, ErrStat, ErrMsg)
    !! Prescribe the host-sampled ambient-fluid field on the cable's nodes: per-node
    !! fluid velocity, acceleration, and local free-surface elevation, held frozen until
    !! the next call (the per-step zero-order hold). Requires the drag config; mutually
    !! exclusive with the internal wave families (the model fails closed). The nodal
    !! vectors arrive in the caller's GLOBAL frame and rotate into the cable's local
    !! solve frame here (frame_to_local per column, exactly as the coupled kinematics
    !! do); the free-surface elevation is invariant under the about-+z heading.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: fluid_velocity(:, :), fluid_acceleration(:, :), waterline_z(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es, i
    LOGICAL :: has_heading
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Held_Fluid: module not initialised'; RETURN
    END IF
    IF (SIZE(fluid_velocity, 1) /= 3 .OR. SIZE(fluid_velocity, 2) /= self%line%nn .OR. &
        SIZE(fluid_acceleration, 1) /= 3 .OR. SIZE(fluid_acceleration, 2) /= self%line%nn) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_Set_Held_Fluid: fluid arrays must be shaped (3, nn)'; RETURN
    END IF
    has_heading = ABS(self%frame_s) > CD_ZERO .OR. ABS(self%frame_c - CD_ONE) > CD_ZERO
    IF (has_heading) THEN
      BLOCK
        REAL(wp) :: u_loc(3, self%line%nn), ud_loc(3, self%line%nn)
        DO i = 1, self%line%nn
          u_loc(:, i) = frame_to_local(self, fluid_velocity(:, i))
          ud_loc(:, i) = frame_to_local(self, fluid_acceleration(:, i))
        END DO
        CALL CD_HermiteCable_Dyn_Set_Held_Fluid(self%line, u_loc, ud_loc, waterline_z, es, em)
      END BLOCK
    ELSE
      CALL CD_HermiteCable_Dyn_Set_Held_Fluid(self%line, fluid_velocity, fluid_acceleration, waterline_z, es, em)
    END IF
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Held_Fluid: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_Held_Fluid

  SUBROUTINE CD_HFMF_Set_Contact(self, kn_node, cn_node, mu, ErrStat, ErrMsg, bathymetry)
    !! Configure standalone nodal seabed contact through the FMF owner. Bathymetry
    !! coordinates are global; the stored cable heading is forwarded for query rotation.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: kn_node(:), cn_node(:), mu
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    TYPE(CD_BathymetryType), INTENT(IN), OPTIONAL :: bathymetry
    INTEGER :: es
    CHARACTER(240) :: em

    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Contact: module not initialised'; RETURN
    END IF
    IF (PRESENT(bathymetry)) THEN
      CALL CD_HermiteCable_Dyn_Set_Contact(self%line, kn_node, cn_node, mu, es, em, &
                                           bathymetry=bathymetry, frame_cs=[self%frame_c, self%frame_s])
    ELSE
      CALL CD_HermiteCable_Dyn_Set_Contact(self%line, kn_node, cn_node, mu, es, em, &
                                           frame_cs=[self%frame_c, self%frame_s])
    END IF
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Set_Contact: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Set_Contact

  SUBROUTINE CD_HFMF_Refresh_Acceleration(self, ErrStat, ErrMsg)
    !! Recompute free-node acceleration after a direct-feedthrough boundary commit.
    !! The caller-prescribed coupled acceleration participates in the consistent-mass
    !! solve, including its off-diagonal coupling to the interior predictor state.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    REAL(wp) :: coupled_a(3)
    CHARACTER(240) :: em

    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Refresh_Acceleration: module not initialised'; RETURN
    END IF
    coupled_a = self%line%a(self%pdof)
    CALL CD_HermiteCable_Dyn_Recompute_Acceleration(self%line, es, em, &
                                                    pres_dofs=self%pdof, pres_a=coupled_a)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Refresh_Acceleration: '//TRIM(em); RETURN
    END IF
  END SUBROUTINE CD_HFMF_Refresh_Acceleration

  INTEGER FUNCTION CD_HFMF_NNodes(self) RESULT(n)
    !! Node count of the owned cable (the held-fluid sampling surface). Zero if
    !! uninitialised.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    n = 0
    IF (self%initialized) n = self%line%nn
  END FUNCTION CD_HFMF_NNodes

  SUBROUTINE CD_HFMF_GetNodePositions(self, xyz, ErrStat, ErrMsg)
    !! Current committed positions of every cable node -- where the host samples its
    !! wave field for CD_HFMF_Set_Held_Fluid. Positions are rotated from the cable's
    !! local solve frame to the caller's global frame (frame_to_global), so the host
    !! always samples at true global coordinates.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    REAL(wp), INTENT(OUT) :: xyz(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: i
    xyz = CD_ZERO
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_GetNodePositions: module not initialised'; RETURN
    END IF
    IF (SIZE(xyz, 1) /= 3 .OR. SIZE(xyz, 2) /= self%line%nn) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_GetNodePositions: array must be shaped (3, nn)'; RETURN
    END IF
    DO i = 1, self%line%nn
      xyz(:, i) = frame_to_global(self, self%line%q(6*(i - 1) + 1:6*(i - 1) + 3))
    END DO
  END SUBROUTINE CD_HFMF_GetNodePositions

  SUBROUTINE CD_HFMF_CalcOutput(self, y_fair, ErrStat, ErrMsg, y_moment)
    !! The FMF loads-out boundary: the force the cable exerts ON the platform at the
    !! coupled endpoint, at the current committed state (global coordinates, N).
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(OUT) :: y_fair(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(OUT), OPTIONAL :: y_moment(3)
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''; y_fair = CD_ZERO
    IF (PRESENT(y_moment)) y_moment = CD_ZERO
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_CalcOutput: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Reaction(self%line, self%coupled_node, y_fair, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_SOLVEFAIL; ErrMsg = 'CD_HFMF_CalcOutput: '//TRIM(em)
      RETURN
    END IF
    y_fair = frame_to_global(self, y_fair)
    IF (PRESENT(y_moment)) THEN
      CALL CD_HermiteCable_Dyn_EndConnection_Moment(self%line, self%coupled_node, y_moment, es, em)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = CD_HFMF_SOLVEFAIL
        ErrMsg = 'CD_HFMF_CalcOutput: '//TRIM(em)
        RETURN
      END IF
      y_moment = frame_to_global(self, y_moment)
    END IF
  END SUBROUTINE CD_HFMF_CalcOutput

  PURE FUNCTION identity_dcm() RESULT(dcm)
    REAL(wp) :: dcm(3, 3)
    dcm = CD_ZERO
    dcm(1, 1) = CD_ONE
    dcm(2, 2) = CD_ONE
    dcm(3, 3) = CD_ONE
  END FUNCTION identity_dcm

  PURE FUNCTION cross3(a, b) RESULT(c)
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  END FUNCTION cross3

  PURE SUBROUTINE direction_rotation(from_direction, to_direction, rotation, ErrStat, ErrMsg)
    !! Shortest proper rotation carrying one unit direction into another. Used
    !! only by the non-stepping rigid-boundary overlay, so tangent state remains
    !! kinematically aligned when OpenFAST probes loads at a new orientation.
    REAL(wp), INTENT(IN) :: from_direction(3), to_direction(3)
    REAL(wp), INTENT(OUT) :: rotation(3, 3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: a(3), b(3), v(3), skew(3, 3), n1, n2, c, s2

    rotation = identity_dcm()
    ErrStat = CD_ENDCONN_OK
    ErrMsg = ''
    n1 = SQRT(DOT_PRODUCT(from_direction, from_direction))
    n2 = SQRT(DOT_PRODUCT(to_direction, to_direction))
    IF (.NOT. CD_Is_Finite(n1) .OR. .NOT. CD_Is_Finite(n2) .OR. &
        n1 <= SQRT(TINY(CD_ONE)) .OR. n2 <= SQRT(TINY(CD_ONE))) THEN
      ErrStat = CD_ENDCONN_BADINPUT
      ErrMsg = 'directions must be finite and non-null'
      RETURN
    END IF
    a = from_direction/n1
    b = to_direction/n2
    c = MAX(-CD_ONE, MIN(CD_ONE, DOT_PRODUCT(a, b)))
    IF (c <= -CD_ONE + 64.0_wp*EPSILON(CD_ONE)) THEN
      ErrStat = CD_ENDCONN_BADINPUT
      ErrMsg = 'a 180-degree direction change has no unique shortest rotation'
      RETURN
    END IF
    v = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    s2 = DOT_PRODUCT(v, v)
    IF (s2 <= 64.0_wp*EPSILON(CD_ONE)) RETURN
    skew = RESHAPE([CD_ZERO, v(3), -v(2), -v(3), CD_ZERO, v(1), &
                    v(2), -v(1), CD_ZERO], [3, 3])
    rotation = rotation + skew + MATMUL(skew, skew)*(CD_ONE - c)/s2
  END SUBROUTINE direction_rotation

  PURE LOGICAL FUNCTION valid_parent_dcm(dcm) RESULT(valid)
    !! OpenFAST mesh orientations are proper global-to-local direction-cosine
    !! matrices. Reject reflections and matrices that would change vector norms.
    REAL(wp), INTENT(IN) :: dcm(3, 3)
    REAL(wp) :: gram(3, 3), det
    valid = .FALSE.
    IF (.NOT. CD_All_Finite(dcm)) RETURN
    gram = MATMUL(dcm, TRANSPOSE(dcm))
    det = dcm(1, 1)*(dcm(2, 2)*dcm(3, 3) - dcm(2, 3)*dcm(3, 2)) - &
          dcm(1, 2)*(dcm(2, 1)*dcm(3, 3) - dcm(2, 3)*dcm(3, 1)) + &
          dcm(1, 3)*(dcm(2, 1)*dcm(3, 2) - dcm(2, 2)*dcm(3, 1))
    valid = MAXVAL(ABS(gram - identity_dcm())) <= 1.0e-8_wp .AND. ABS(det - CD_ONE) <= 1.0e-8_wp
  END FUNCTION valid_parent_dcm

  PURE FUNCTION frame_is_identity(self) RESULT(is_id)
    !! Exact identity-heading test. Order comparisons preserve the bit-for-bit identity
    !! fast path without tripping -Wcompare-reals.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    LOGICAL :: is_id
    is_id = .NOT. (ABS(self%frame_s) > CD_ZERO .OR. ABS(self%frame_c - CD_ONE) > CD_ZERO)
  END FUNCTION frame_is_identity

  PURE FUNCTION frame_to_local(self, v) RESULT(vl)
    !! Rotate a caller-global vector into the cable's local solve frame (about +z by the
    !! stored heading). Identity heading returns v unchanged, preserving bit-for-bit
    !! behavior on every pre-existing (x-z plane) cable.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    REAL(wp), INTENT(IN) :: v(3)
    REAL(wp) :: vl(3)
    IF (frame_is_identity(self)) THEN
      vl = v
    ELSE
      vl = [self%frame_c*v(1) + self%frame_s*v(2), &
            -self%frame_s*v(1) + self%frame_c*v(2), v(3)]
    END IF
  END FUNCTION frame_to_local

  PURE FUNCTION frame_to_global(self, v) RESULT(vg)
    !! Rotate a cable-local vector back to the caller's global frame (the transpose of
    !! frame_to_local). Identity heading returns v unchanged.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    REAL(wp), INTENT(IN) :: v(3)
    REAL(wp) :: vg(3)
    IF (frame_is_identity(self)) THEN
      vg = v
    ELSE
      vg = [self%frame_c*v(1) - self%frame_s*v(2), &
            self%frame_s*v(1) + self%frame_c*v(2), v(3)]
    END IF
  END FUNCTION frame_to_global

  SUBROUTINE CD_HFMF_Curvature(self, curv_out, ErrStat, ErrMsg)
    !! Nodal centreline curvature at the current state (the fatigue report channel).
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(OUT) :: curv_out(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(300) :: em
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Curvature: module not initialised'; RETURN
    END IF
    CALL CD_HermiteCable_Dyn_Curvature(self%line, curv_out, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_Curvature: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HFMF_Curvature

  REAL(wp) FUNCTION CD_HFMF_MinSpanZ(self) RESULT(zmin)
    !! Minimum global z along the cubic-Hermite CENTRELINE (its static-seeded shape), computed
    !! EXACTLY per element from the cubic's stationary points -- NOT by fixed sampling. On each
    !! element z(xi) is a cubic on [0,1]; a cubic can dip below both endpoint nodes (and below a
    !! declared seabed) BETWEEN any finite set of samples, so its minimum is taken over the two
    !! endpoints {xi=0, xi=1} plus any real root of dz/dxi (a quadratic) that lies in (0,1) --
    !! exact for the cubic. The aggregate uses this to decide whether a suspended cable actually
    !! reaches a declared seabed (a grounded/touchdown cable the pinned kn=0 static cannot
    !! represent). Returns HUGE for an uninitialised module (no shape -> never flagged grounded).
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    INTEGER :: e, b1, b2, nc, ic
    REAL(wp) :: z1, mz1, z2, mz2, Lel
    REAL(wp) :: Aq, Bq, Cq, dz0, dz1, dzh, disc, sq, scal
    REAL(wp) :: cand(4), Hs(4), dHs(4), zc
    zmin = HUGE(0.0_wp)
    IF (.NOT. self%initialized) RETURN
    DO e = 1, self%line%ne
      b1 = 6*(e - 1)                          ! node e   base DOF (r,m): z at +3, tangent-z at +6
      b2 = 6*e                                ! node e+1 base DOF
      Lel = self%line%l0(e)
      z1 = self%line%q(b1 + 3); mz1 = self%line%q(b1 + 6)
      z2 = self%line%q(b2 + 3); mz2 = self%line%q(b2 + 6)
      ! z(xi) is a cubic; dz/dxi(xi) = Aq xi^2 + Bq xi + Cq is a quadratic. Recover A, B, C from
      ! dz/dxi at xi = 0, 1, 1/2:  Cq = dz(0);  Aq+Bq+Cq = dz(1);  Aq/4+Bq/2+Cq = dz(1/2).
      CALL CD_HermiteCable_Shapes(0.0_wp, Lel, Hs, dHs)
      dz0 = dHs(1)*z1 + dHs(2)*mz1 + dHs(3)*z2 + dHs(4)*mz2
      CALL CD_HermiteCable_Shapes(1.0_wp, Lel, Hs, dHs)
      dz1 = dHs(1)*z1 + dHs(2)*mz1 + dHs(3)*z2 + dHs(4)*mz2
      CALL CD_HermiteCable_Shapes(0.5_wp, Lel, Hs, dHs)
      dzh = dHs(1)*z1 + dHs(2)*mz1 + dHs(3)*z2 + dHs(4)*mz2
      Cq = dz0
      Aq = 2.0_wp*(dz1 + Cq - 2.0_wp*dzh)
      Bq = dz1 - Cq - Aq
      ! Candidate abscissae: the two endpoints always, plus real roots of dz/dxi strictly in (0,1).
      nc = 2
      cand(1) = 0.0_wp
      cand(2) = 1.0_wp
      scal = ABS(Aq) + ABS(Bq) + ABS(Cq)
      IF (scal > 0.0_wp) THEN
        IF (ABS(Aq) <= 1.0e-12_wp*scal) THEN
          ! (near-)linear dz/dxi: a single stationary point at -Cq/Bq when Bq is non-degenerate.
          IF (ABS(Bq) > 1.0e-12_wp*scal) CALL add_interior_root(-Cq/Bq, cand, nc)
        ELSE
          disc = Bq*Bq - 4.0_wp*Aq*Cq
          IF (disc >= 0.0_wp) THEN
            sq = SQRT(disc)
            CALL add_interior_root((-Bq + sq)/(2.0_wp*Aq), cand, nc)
            CALL add_interior_root((-Bq - sq)/(2.0_wp*Aq), cand, nc)
          END IF
        END IF
      END IF
      DO ic = 1, nc
        CALL CD_HermiteCable_Shapes(cand(ic), Lel, Hs, dHs)
        ! z(xi) = H1 z1 + H2 mz1 + H3 z2 + H4 mz2  (DOF order [r1, m1, r2, m2]; z is comp 3)
        zc = Hs(1)*z1 + Hs(2)*mz1 + Hs(3)*z2 + Hs(4)*mz2
        zmin = MIN(zmin, zc)
      END DO
    END DO

  CONTAINS
    SUBROUTINE add_interior_root(xi, xs, n)
      !! Append xi to the candidate list xs (count n) only when it is a strict interior stationary
      !! point (0,1); endpoints are already candidates, so a root at/beyond a boundary adds nothing.
      REAL(wp), INTENT(IN) :: xi
      REAL(wp), INTENT(INOUT) :: xs(:)
      INTEGER, INTENT(INOUT) :: n
      IF (xi > 0.0_wp .AND. xi < 1.0_wp) THEN
        n = n + 1
        xs(n) = xi
      END IF
    END SUBROUTINE add_interior_root
  END FUNCTION CD_HFMF_MinSpanZ

  INTEGER FUNCTION CD_HFMF_MirrorSize(self) RESULT(n)
    !! Number of continuous-state mirror slots this cable contributes to the OpenFAST state
    !! vector: the free (non-prescribed, non-fixed) DOF velocities, positions, and
    !! accelerations, [v_free; q_free; a_free], so n = 3 * (number of free DOFs). The
    !! acceleration block exists for CHECKPOINT-RESTART: the generalised-alpha step
    !! carries the committed acceleration, and unlike the EI=0 side (whose acceleration
    !! is re-derived by CD_Recompute_System_Acceleration) the cable reloads it exactly.
    !! A cable on a frictional seabed appends its stick-slip friction anchors, (x, y) per
    !! node (CD_HFMF_FrictionMirrorSize slots): they are committed state like q. Zero if
    !! uninitialised.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    n = 0
    IF (.NOT. self%initialized) RETURN
    IF (ALLOCATED(self%line%freemask)) n = 3*COUNT(self%line%freemask)
    n = n + CD_HFMF_FrictionMirrorSize(self) + CD_HFMF_ForceMirrorSize(self)
  END FUNCTION CD_HFMF_MirrorSize

  INTEGER FUNCTION CD_HFMF_ForceMirrorSize(self) RESULT(n)
    !! Mirror slots of the force-blended generalised-alpha integrator state beyond [v; q; a] of
    !! the free DOFs: a validity flag and the committed force F_n (ndof values), then the
    !! committed accelerations of the held DOFs; zero with the configuration blend. The
    !! force blend carries F_n and a_n of the prescribed end into the next step, so a restart
    !! reloads them instead of re-evaluating F_n and taking a_n from the host input mesh.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    n = 0
    IF (.NOT. self%initialized) RETURN
    IF (self%line%force_blend) n = 1 + self%line%ndof + (self%line%ndof - COUNT(self%line%freemask))
  END FUNCTION CD_HFMF_ForceMirrorSize

  INTEGER FUNCTION CD_HFMF_FrictionMirrorSize(self) RESULT(n)
    !! Mirror slots of the stick-slip friction anchors (2 per node), zero without friction.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    n = 0
    IF (.NOT. self%initialized) RETURN
    IF (self%line%fr_active) n = 2*self%line%nn
  END FUNCTION CD_HFMF_FrictionMirrorSize

  SUBROUTINE CD_HFMF_PackMirror(self, buf, ErrStat, ErrMsg)
    !! Write the cable's continuous state into buf as [v_free; q_free; a_free] -- the
    !! layout CD_HFMF_MirrorSize sizes. The mirror is written every mooring step and
    !! reloaded ONLY by a checkpoint-restart rebuild (CD_HFMF_UnpackMirror); the module
    !! stays self-authoritative in normal marching.
    TYPE(CD_HFMF_ModuleType), INTENT(IN) :: self
    REAL(wp), INTENT(OUT) :: buf(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: i, nf, k
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    buf = CD_ZERO
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_PackMirror: module not initialised'; RETURN
    END IF
    nf = COUNT(self%line%freemask)
    IF (SIZE(buf) /= CD_HFMF_MirrorSize(self)) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_PackMirror: buffer must be sized CD_HFMF_MirrorSize'
      RETURN
    END IF
    k = 0
    DO i = 1, self%line%ndof
      IF (self%line%freemask(i)) THEN
        k = k + 1
        buf(k) = self%line%v(i)           ! velocity block
        buf(nf + k) = self%line%q(i)      ! position block
        buf(2*nf + k) = self%line%a(i)    ! acceleration block (gen-alpha committed state)
      END IF
    END DO
    k = 3*nf
    ! stick-slip friction anchors, (x, y) per node, after the acceleration block
    IF (self%line%fr_active) THEN
      buf(k + 1:k + 2*self%line%nn) = RESHAPE(self%line%fr_anchor, [2*self%line%nn])
      k = k + 2*self%line%nn
    END IF
    ! force-blend committed force: flag (1 when the cache holds F at the committed state)
    IF (self%line%force_blend) THEN
      IF (self%line%fc_valid .AND. .NOT. (ABS(self%line%fc_t - self%line%t) > CD_ZERO)) THEN
        buf(k + 1) = 1.0_wp
        buf(k + 2:k + 1 + self%line%ndof) = self%line%fc_R
      END IF
      k = k + 1 + self%line%ndof
      DO i = 1, self%line%ndof
        IF (self%line%freemask(i)) CYCLE
        k = k + 1
        buf(k) = self%line%a(i)
      END DO
    END IF
  END SUBROUTINE CD_HFMF_PackMirror

  SUBROUTINE CD_HFMF_UnpackMirror(self, buf, t, ErrStat, ErrMsg)
    !! CHECKPOINT-RESTART reload: overwrite the cable's free-DOF state from a mirror
    !! written by CD_HFMF_PackMirror ([v_free; q_free; a_free]) and set the line clock
    !! to the restart time. The prescribed (coupled) and fixed DOFs are untouched -- the
    !! caller re-prescribes the coupled endpoint from the host mesh
    !! (CD_HFMF_SetCoupledKinematics) around this call. The module must already be
    !! initialised from the SAME deck (the rebuild re-runs the static build first).
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    REAL(wp), INTENT(IN) :: buf(:)
    REAL(wp), INTENT(IN) :: t
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: i, nf, k
    ErrStat = CD_HFMF_OK; ErrMsg = ''
    IF (.NOT. self%initialized) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_UnpackMirror: module not initialised'; RETURN
    END IF
    IF (.NOT. CD_Is_Finite(t)) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_UnpackMirror: restart time must be finite'; RETURN
    END IF
    nf = COUNT(self%line%freemask)
    IF (SIZE(buf) /= CD_HFMF_MirrorSize(self)) THEN
      ErrStat = CD_HFMF_BADINPUT
      ErrMsg = 'CD_HFMF_UnpackMirror: buffer must be sized CD_HFMF_MirrorSize'
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(buf)) THEN
      ErrStat = CD_HFMF_BADINPUT; ErrMsg = 'CD_HFMF_UnpackMirror: mirror must be finite'; RETURN
    END IF
    k = 0
    DO i = 1, self%line%ndof
      IF (self%line%freemask(i)) THEN
        k = k + 1
        self%line%v(i) = buf(k)
        self%line%q(i) = buf(nf + k)
        self%line%a(i) = buf(2*nf + k)
      END IF
    END DO
    k = 3*nf
    IF (self%line%fr_active) THEN
      self%line%fr_anchor = RESHAPE(buf(k + 1:k + 2*self%line%nn), [2, self%line%nn])
      self%line%ws_recovery_fr = self%line%fr_anchor
      k = k + 2*self%line%nn
    END IF
    self%line%t = t
    ! The reloaded acceleration is committed history: a first held fluid field after the
    ! reload (the host re-sampling its field on the restored state) must keep it.
    self%line%a_restored = .TRUE.
    ! the force-blend committed force reloads with the state it was evaluated at
    self%line%fc_valid = .FALSE.
    IF (self%line%force_blend) THEN
      IF (buf(k + 1) > 0.5_wp) THEN
        self%line%fc_R = buf(k + 2:k + 1 + self%line%ndof)
        self%line%fc_q = self%line%q
        self%line%fc_v = self%line%v
        self%line%fc_t = self%line%t
        self%line%fc_d0 = self%line%endconn_d0
        self%line%fc_valid = .TRUE.
      END IF
      k = k + 1 + self%line%ndof
      DO i = 1, self%line%ndof
        IF (self%line%freemask(i)) CYCLE
        k = k + 1
        self%line%a(i) = buf(k)
      END DO
    END IF
  END SUBROUTINE CD_HFMF_UnpackMirror

  SUBROUTINE CD_HFMF_End(self)
    !! Release the module and its line model.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: self
    CALL CD_HermiteCable_Dyn_End(self%line)
    self%coupled_node = 0
    self%pdof = 0
    self%dt = CD_ZERO
    self%max_iter = 100
    self%tol = 1.0e-4_wp
    self%recovery_max_substeps = 1024
    self%frame_c = CD_ONE
    self%frame_s = CD_ZERO
    self%has_parent_endconn = .FALSE.
    self%parent_endconn_index = 0
    self%parent_endconn_d0 = CD_ZERO
    self%parent_dcm = identity_dcm()
    self%parent_omega = CD_ZERO
    self%parent_alpha = CD_ZERO
    self%snap_parent_dcm = identity_dcm()
    self%snap_parent_omega = CD_ZERO
    self%snap_parent_alpha = CD_ZERO
    self%initialized = .FALSE.
  END SUBROUTINE CD_HFMF_End

END MODULE CableDyn_OpenFAST_HermiteFMF
