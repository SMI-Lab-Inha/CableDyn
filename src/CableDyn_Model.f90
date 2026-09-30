! File: src/CableDyn_Model.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_Model
  !! Persistent solver model for the Fortran production core. This module is the
  !! lifecycle layer over the EI=0 cable dynamics path described in ARCHITECTURE.md
  !! `src/` and doc/coupling_boundary.md: callers initialise a model once, advance it
  !! by steps, query its state, and release it through one Fortran-native API. The
  !! OpenFAST shell and the C API call this layer rather than owning element arrays
  !! directly.
  !!
  !! Scope of this module:
  !!  * positions-only EI=0 line dynamics,
  !!  * line-object initialisation from CableDyn_Line sections/types,
  !!  * external, hydrodynamic, seabed, and constitutive (viscoelastic / Syrope) loads,
  !!  * optional prescribed motion on the declared fixed/coupled DOFs,
  !!  * generalized-alpha integration via CableDyn_Dynamic.
  !!
  !! The element, assembly, and time-integration kernels live in CableDyn_Dynamic,
  !! CableDyn_Assemble, and CableDyn_Mesh; this module owns the object boundary the
  !! OpenFAST and C shells consume.
  USE CableDyn_Precision, ONLY: wp, CD_ZERO, CD_ONE, CD_All_Finite, CD_Is_Finite
  USE CableDyn_Line, ONLY: CD_LineType, CD_LineSection, CD_Build_Line_Mesh, CD_LINE_OK
  USE CableDyn_Loads, ONLY: CD_Submerged_Weight, CD_Assemble_Distributed_Load, CD_Seabed_Penalty_Load
  USE CableDyn_Line, ONLY: CD_Nodal_Seabed_Stiffness
  USE CableDyn_CableElem, ONLY: CD_Compute_Cable_Element
  USE CableDyn_SeabedContact, ONLY: CD_Seabed_Normal_Law, CD_Seabed_Normal_Contact, CD_Seabed_Contact_Gate, &
                                    CD_Seabed_Friction_Stick_Slip, CD_Seabed_Friction_Aniso, &
                                    CD_Seabed_Friction_Mu_Dir, CD_FRICTION_STICK_SLIP
  USE CableDyn_Bathymetry, ONLY: CD_BathymetryType, CD_Init_Bathymetry, CD_Bathymetry_Is_Initialized, &
                                 CD_Bathymetry_Floor, CD_Bathymetry_Floor_Gradient, &
                                 CD_Bathymetry_Seabed_Load, CD_End_Bathymetry, CD_BATHY_OK
  USE CableDyn_Static, ONLY: CableSolverConfig, CableCurrentLoad, CD_StaticLineReport, CD_Static_Line_Equilibrium, &
                             CD_Static_Friction_Anchors, &
                             CD_STATIC_OK
  USE CableDyn_Assemble, ONLY: CD_Assemble_Cable_Internal_Force, CD_Assemble_Cable_Mass, &
                               CD_Assemble_Cable_Tangent_Force, CD_Compute_Cable_Tension
  USE CableDyn_Damping, ONLY: CD_Cable_Axial_Damping_Load, CD_Cable_Axial_Damping_Force, &
                              CD_Axial_Damping_Element_Load, CD_Axial_Damping_Element_Force
  USE CableDyn_Viscoelastic, ONLY: CD_Viscoelastic_Params, CD_Viscoelastic_LoadDependent_EAD, &
                                   CD_Viscoelastic_Steady_State, CD_Viscoelastic_Element_Load, &
                                   CD_Viscoelastic_Element_Force, CD_Viscoelastic_Element_Tension, &
                                   CD_Viscoelastic_State_Advance, CD_VISCO_OK, CD_VISCO_BADINPUT, &
                                   CD_VISCO_TINY_LEN
  USE CableDyn_Syrope, ONLY: CD_SyropeType, CD_Syrope_Init, CD_Syrope_End, CD_Syrope_Is_Ready, &
                             CD_Syrope_Working_Curve, CD_Syrope_Rate_And_Tension, CD_Syrope_Branch_Divider, &
                             CD_Syrope_Element_Load, CD_Syrope_Check_Range, CD_SYROPE_OK, CD_SYROPE_NWC
  USE CableDyn_Hydro, ONLY: CD_Cable_Morison_Drag_Load, CD_Cable_Morison_Drag_Force, &
                            CD_Morison_Drag_Element_Load, CD_Morison_Drag_Element_Force, &
                            CD_Cable_Froude_Krylov_Load, CD_Cable_Froude_Krylov_Force, &
                            CD_Froude_Krylov_Element_Load, CD_Froude_Krylov_Element_Force, &
                            CD_Cable_Buoyancy_Recovery_Load, CD_Cable_Buoyancy_Recovery_Force, &
                            CD_Buoyancy_Recovery_Element_Load, CD_Buoyancy_Recovery_Element_Force, &
                            CD_Cable_Added_Mass, CD_Cable_Added_Mass_Matrix, CD_Cable_Added_Mass_Force, &
                            CD_Added_Mass_Element, CD_Added_Mass_Element_Matrix, CD_Added_Mass_Element_Force
  USE CableDyn_Dynamic, ONLY: GenAlphaConfig, CD_CableGenAlphaWorkspace, CD_Clear_GenAlpha_Workspace, &
                              CD_Cable_Initial_Acceleration, CD_Cable_Gen_Alpha_Step, CD_DYN_OK, CD_DYN_ALLOCFAIL
  USE CableDyn_Mesh, ONLY: CD_Partition_Free_Dofs, CD_Validate_Connectivity
  USE CableDyn_Linalg, ONLY: CD_Solve_Dense_As_Banded, CD_Solve_Dense_As_Banded_Multiple, CD_Factor_Banded, &
                             CD_Solve_Factored_Banded_Multiple
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CD_ModelType
  PUBLIC :: CD_Init_Model
  PUBLIC :: CD_Init_Line_Model
  PUBLIC :: CD_Step_Model
  PUBLIC :: CD_Step_Model_Recovering
  PUBLIC :: CD_End_Model
  PUBLIC :: CD_Model_NCoupledDOF
  PUBLIC :: CD_Model_NDOF
  PUBLIC :: CD_Model_NElem
  PUBLIC :: CD_Get_Model_CoupledDofs
  PUBLIC :: CD_Get_Model_CoupledMotion
  PUBLIC :: CD_Get_Model_State
  PUBLIC :: CD_Model_Has_Viscoelastic
  PUBLIC :: CD_Get_Model_VE_Dl1
  PUBLIC :: CD_Set_Model_VE_Dl1
  PUBLIC :: CD_Model_Has_Syrope
  PUBLIC :: CD_Model_Has_Friction
  PUBLIC :: CD_Get_Model_Friction_Anchors
  PUBLIC :: CD_Set_Model_Friction_Anchors
  PUBLIC :: CD_Set_Model_Friction_Axial
  PUBLIC :: CD_Get_Model_Syrope_State
  PUBLIC :: CD_Set_Model_Syrope_State
  PUBLIC :: CD_Get_Model_Tension
  PUBLIC :: CD_Get_Model_EndForces
  PUBLIC :: CD_Get_Model_EndNodeMassDiag
  PUBLIC :: CD_Get_Model_EndNodeAddedMass
  PUBLIC :: CD_Get_Model_EndNodeTangent
  PUBLIC :: CD_Set_Model_EndQuery
  PUBLIC :: CD_Get_Model_EndQuery
  PUBLIC :: CD_Update_Model_SegmentLength
  PUBLIC :: CD_Calc_Model_CoupledLoads
  PUBLIC :: CD_Calc_Model_CoupledKinematicDerivatives
  PUBLIC :: CD_Calc_Model_CoupledAccelDerivative, CD_Update_Model_CoupledAccelDerivative
  PUBLIC :: CD_Update_Model_Hydro_Fields
  PUBLIC :: CD_Update_Model_External_Loads
  PUBLIC :: CD_Update_Model_CoupledMotion
  PUBLIC :: CD_Update_Model_Interior_State
  PUBLIC :: CD_Recompute_Model_Acceleration
  PUBLIC :: CD_Model_Is_Initialized
  PUBLIC :: CD_Copy_Model

  INTEGER, PARAMETER, PUBLIC :: CD_MODEL_OK = 0
  INTEGER, PARAMETER, PUBLIC :: CD_MODEL_BADINPUT = 1
  INTEGER, PARAMETER, PUBLIC :: CD_MODEL_SOLVEFAIL = 2
  INTEGER, PARAMETER, PUBLIC :: CD_MODEL_NOT_INITIALIZED = 3
  INTEGER, PARAMETER, PUBLIC :: CD_MODEL_ALLOCFAIL = 4
  INTEGER, PARAMETER :: MODEL_DEFAULT_MAX_SUBSTEPS = 1024

  TYPE :: CD_ModelType
    !! Persistent positions-only EI=0 model state and immutable mesh/properties.
    LOGICAL :: initialized = .FALSE.
    LOGICAL :: tension_only = .FALSE.
    INTEGER :: n_dof = 0
    INTEGER :: n_elem = 0
    INTEGER, ALLOCATABLE :: elem_conn(:, :)
    INTEGER, ALLOCATABLE :: fixed_dofs(:)
    REAL(wp), ALLOCATABLE :: l0(:)
    REAL(wp), ALLOCATABLE :: ea(:)
    REAL(wp), ALLOCATABLE :: rho_a(:)
    ! Per-element unstretched-length rate (active line control; zero everywhere
    ! until CD_Update_Model_SegmentLength sets it). Feeds the axial-damping
    ! effective stretch rate.
    REAL(wp), ALLOCATABLE :: l0_dot(:)
    ! Per-element distributed load per unit LENGTH (3, n_elem) plus the
    ! NON-distributed remainder of the init f_ext, stored when the caller supplies
    ! the per-length load: a segment-length update REASSEMBLES the distributed part
    ! from the current lengths (f_ext = remainder + assemble(dist_load, l0)), which
    ! is bit-identical to a fresh build -- a nodal delta (0.5*dl*w added to an
    ! existing share) is NOT, in floating point. Unallocated = no distributed load
    ! to track (point loads and zero f_ext are length-independent).
    REAL(wp), ALLOCATABLE :: dist_load(:, :)
    REAL(wp), ALLOCATABLE :: f_ext_dist_base(:)
    ! Seabed nodal-coefficient recompute inputs (active line control): when the
    ! caller supplies them, a segment-length update rebuilds seabed_kn (and cn via
    ! the stored ratio) through CD_Nodal_Seabed_Stiffness at the current lengths --
    ! the SAME helper the deck used, so a fresh build matches bit-for-bit. Without
    ! them a length update on a seabed model fails closed (a controlled grounded
    ! segment would otherwise keep its old tributary contact area).
    LOGICAL :: has_seabed_recompute = .FALSE.
    REAL(wp) :: seabed_kbot = CD_ZERO
    REAL(wp) :: seabed_cn_over_kn = CD_ZERO
    REAL(wp), ALLOCATABLE :: seabed_diameter_elem(:)
    REAL(wp), ALLOCATABLE :: f_ext(:)
    REAL(wp), ALLOCATABLE :: q(:)
    REAL(wp), ALLOCATABLE :: v(:)
    REAL(wp), ALLOCATABLE :: a(:)
    REAL(wp), ALLOCATABLE :: q_step(:)
    REAL(wp), ALLOCATABLE :: v_step(:)
    REAL(wp), ALLOCATABLE :: a_step(:)
    REAL(wp), ALLOCATABLE :: q_prescribed_work(:)
    REAL(wp), ALLOCATABLE :: v_prescribed_work(:)
    REAL(wp), ALLOCATABLE :: a_prescribed_work(:)
    REAL(wp), ALLOCATABLE :: q_rollback(:)
    REAL(wp), ALLOCATABLE :: v_rollback(:)
    REAL(wp), ALLOCATABLE :: a_rollback(:)
    ! Dedicated nominal-interval snapshot for CD_Step_Model_Recovering.  These
    ! arrays are separate from q_rollback/v_rollback/a_rollback, which belong to
    ! the multi-line system and transactional load-update contracts.
    REAL(wp), ALLOCATABLE :: recovery_q0(:)
    REAL(wp), ALLOCATABLE :: recovery_v0(:)
    REAL(wp), ALLOCATABLE :: recovery_a0(:)
    REAL(wp), ALLOCATABLE :: f_ext_rollback(:)
    INTEGER, ALLOCATABLE :: free_map_work(:)
    TYPE(CD_CableGenAlphaWorkspace) :: dynamic_workspace
    REAL(wp), ALLOCATABLE :: load_force_work(:)
    REAL(wp), ALLOCATABLE :: load_jq_work(:, :)
    REAL(wp), ALLOCATABLE :: load_jv_work(:, :)
    LOGICAL :: has_seabed = .FALSE.
    LOGICAL :: has_bathymetry = .FALSE.
    REAL(wp) :: seabed_z_floor = CD_ZERO
    TYPE(CD_BathymetryType) :: bathymetry
    REAL(wp), ALLOCATABLE :: seabed_kn(:)
    LOGICAL :: has_seabed_damping = .FALSE.
    REAL(wp), ALLOCATABLE :: seabed_cn(:)
    LOGICAL :: has_seabed_friction = .FALSE.
    REAL(wp) :: seabed_mu = CD_ZERO
    ! Anisotropic friction (CD_Set_Model_Friction_Axial): seabed_mu is then the lateral
    ! coefficient and seabed_mu_axial the one along the line (CD_Seabed_Friction_Aniso).
    LOGICAL :: seabed_friction_aniso = .FALSE.
    REAL(wp) :: seabed_mu_axial = CD_ZERO
    ! Stick-slip seabed friction (has_seabed_friction): the committed anchor (x, y) of each
    ! node's friction spring (stiffness seabed_kn, capacity mu times the normal reaction),
    ! a STATE like q, with its system-rollback and step-recovery twins.
    REAL(wp), ALLOCATABLE :: fr_anchor(:, :), fr_anchor_rollback(:, :), recovery_fr_anchor(:, :)
    LOGICAL :: has_damping = .FALSE.
    REAL(wp), ALLOCATABLE :: ba(:)
    ! Viscoelastic axial model (MoorDyn ElasticMod 2, two Kelvin branches in series):
    ! per-element parameters (ve_ea_d == 0 marks a plain element) and the
    ! committed static-branch stretch state dl_1. model%ea keeps the relaxed
    ! composite everywhere it means STATIC stiffness (tension recovery of plain
    ! elements, coupled fscale); the DYNAMIC structural assembly and the
    ! coupled-load residual use ea_dyn (viscoelastic entries zeroed) with the
    ! series-Kelvin load contributor supplying the full elastic + damping tension for
    ! those elements, and ba_dyn likewise masks the plain axial damping there
    ! (the static-branch dashpot Bs lives inside the series law). ve_dt carries the
    ! active step size into the contributors; 0 = the instantaneous
    ! (committed-state) limit used by IC, coupled-load, and tension queries.
    ! ve_dl_1_rollback is the system-level rollback twin of q_rollback: a line
    ! whose step committed (advancing dl_1) can be rolled back when a sibling
    ! line or the coupled point solve fails.
    ! ve_mode: 0 = plain element, 2 = constant dynamic stiffness (Es|Ed),
    ! 3 = mean-load-dependent dynamic stiffness (Es|alphaMBL|vbeta). Mode 3
    ! evaluates EA_D from the COMMITTED dl_1 (lagged over the step, so the
    ! in-step elimination Jacobians stay exact for the lagged operator).
    LOGICAL :: has_viscoelastic = .FALSE.
    INTEGER, ALLOCATABLE :: ve_mode(:)
    REAL(wp), ALLOCATABLE :: ve_alpha_mbl(:)
    REAL(wp), ALLOCATABLE :: ve_vbeta(:)
    REAL(wp), ALLOCATABLE :: ve_ea_d(:)
    REAL(wp), ALLOCATABLE :: ve_ea_1(:)
    REAL(wp), ALLOCATABLE :: ve_ba(:)
    REAL(wp), ALLOCATABLE :: ve_ba_d(:)
    REAL(wp), ALLOCATABLE :: ve_dl_1(:)
    REAL(wp), ALLOCATABLE :: ve_dl_1_rollback(:)
    REAL(wp), ALLOCATABLE :: recovery_ve_dl_1(:)
    REAL(wp), ALLOCATABLE :: ea_dyn(:)
    REAL(wp), ALLOCATABLE :: ba_dyn(:)
    REAL(wp) :: ve_dt = CD_ZERO
    ! Syrope polyester working-curve constitutive (the SYROPE deck dialect): a
    ! load-history-dependent alternative to the series-Kelvin law. syrope_is(:) marks the Syrope
    ! elements; syrope_type(:) holds each element's constitutive data (OWC
    ! table + working-curve shape + spring/dashpot constants -- immutable after
    ! init). syrope_slow is the committed slow-spring static strain (the STATE)
    ! and syrope_tmax the running maximum mean tension (a second, monotone STATE
    ! whose growth regenerates the working curve on demand). Syrope elements
    ! share the ea_dyn/ba_dyn masking with the viscoelastic path (their elastic
    ! EA and plain damping are zeroed; the Syrope contributor supplies the
    ! tension, which has no direct velocity term). The _rollback twins mirror
    ! q_rollback for the system-level stage-then-commit contract.
    LOGICAL :: has_syrope = .FALSE.
    LOGICAL, ALLOCATABLE :: syrope_is(:)
    TYPE(CD_SyropeType), ALLOCATABLE :: syrope_type(:)
    REAL(wp), ALLOCATABLE :: syrope_slow(:)
    REAL(wp), ALLOCATABLE :: syrope_slow_rollback(:)
    REAL(wp), ALLOCATABLE :: syrope_tmax(:)
    REAL(wp), ALLOCATABLE :: syrope_tmax_rollback(:)
    REAL(wp), ALLOCATABLE :: recovery_syrope_slow(:)
    REAL(wp), ALLOCATABLE :: recovery_syrope_tmax(:)
    LOGICAL :: has_morison_drag = .FALSE.
    REAL(wp), ALLOCATABLE :: fluid_velocity(:, :)
    REAL(wp), ALLOCATABLE :: drag_waterline_z(:)
    REAL(wp) :: drag_rho = CD_ZERO
    REAL(wp) :: drag_diameter = CD_ZERO
    REAL(wp) :: drag_cdn = CD_ZERO
    REAL(wp) :: drag_cdt = CD_ZERO
    REAL(wp), ALLOCATABLE :: drag_diameter_elem(:)
    REAL(wp), ALLOCATABLE :: drag_cdn_elem(:)
    REAL(wp), ALLOCATABLE :: drag_cdt_elem(:)
    LOGICAL :: has_froude_krylov = .FALSE.
    REAL(wp), ALLOCATABLE :: fluid_acceleration(:, :)
    REAL(wp), ALLOCATABLE :: fk_waterline_z(:)
    REAL(wp) :: fk_rho = CD_ZERO
    REAL(wp) :: fk_diameter = CD_ZERO
    REAL(wp) :: fk_can = CD_ZERO
    REAL(wp) :: fk_cat = CD_ZERO
    REAL(wp), ALLOCATABLE :: fk_diameter_elem(:)
    REAL(wp), ALLOCATABLE :: fk_can_elem(:)
    REAL(wp), ALLOCATABLE :: fk_cat_elem(:)
    LOGICAL :: has_buoyancy_recovery = .FALSE.
    REAL(wp), ALLOCATABLE :: buoyancy_waterline_z(:)
    REAL(wp) :: buoyancy_rho = CD_ZERO
    REAL(wp) :: buoyancy_diameter = CD_ZERO
    REAL(wp) :: buoyancy_gravity = CD_ZERO
    REAL(wp), ALLOCATABLE :: buoyancy_diameter_elem(:)
    LOGICAL :: has_added_mass = .FALSE.
    REAL(wp), ALLOCATABLE :: added_mass_waterline_z(:)
    REAL(wp) :: added_mass_rho = CD_ZERO
    REAL(wp) :: added_mass_diameter = CD_ZERO
    REAL(wp) :: added_mass_can = CD_ZERO
    REAL(wp) :: added_mass_cat = CD_ZERO
    REAL(wp), ALLOCATABLE :: added_mass_diameter_elem(:)
    REAL(wp), ALLOCATABLE :: added_mass_can_elem(:)
    REAL(wp), ALLOCATABLE :: added_mass_cat_elem(:)
    TYPE(GenAlphaConfig) :: cfg
  END TYPE CD_ModelType

CONTAINS

  SUBROUTINE CD_Init_Model(model, q0, v0, elem_conn, l0, ea, rho_a, tension_only, &
                           f_ext, fixed_dofs, cfg, ErrStat, ErrMsg, seabed_z_floor, seabed_kn, seabed_cn, &
                           seabed_mu, ba, &
                           fluid_velocity, drag_waterline_z, drag_rho, drag_diameter, drag_cdn, drag_cdt, &
                           drag_diameter_elem, drag_cdn_elem, drag_cdt_elem, &
                           fluid_acceleration, fk_waterline_z, fk_rho, fk_diameter, fk_can, fk_cat, &
                           fk_diameter_elem, fk_can_elem, fk_cat_elem, &
                           buoyancy_waterline_z, buoyancy_rho, buoyancy_diameter, buoyancy_gravity, &
                           buoyancy_diameter_elem, &
                           added_mass_waterline_z, added_mass_rho, added_mass_diameter, added_mass_can, &
                           added_mass_cat, added_mass_diameter_elem, added_mass_can_elem, added_mass_cat_elem, &
                           bathymetry, &
                           dist_load_per_length, seabed_kbot, seabed_contact_diameter, seabed_cn_over_kn, &
                           ve_ea_d, ve_ba, ve_ba_d, ve_alpha_mbl, ve_vbeta, &
                           syrope_is, syrope_type, syrope_slow0, syrope_tmax0)
    !! Initialise a persistent EI=0 cable model from already-meshed arrays. The
    !! initial acceleration is computed consistently from the current structure,
    !! load, and fixed DOFs; callers do not provide a0. Existing model storage is
    !! released before the new model is populated.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: q0(:), v0(:)
    INTEGER, INTENT(IN) :: elem_conn(:, :)
    REAL(wp), INTENT(IN) :: l0(:), ea(:), rho_a(:)
    LOGICAL, INTENT(IN) :: tension_only
    REAL(wp), INTENT(IN) :: f_ext(:)
    INTEGER, INTENT(IN) :: fixed_dofs(:)
    TYPE(GenAlphaConfig), INTENT(IN) :: cfg
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_z_floor
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_kn(:)
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_cn(:)
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_mu
    REAL(wp), INTENT(IN), OPTIONAL :: ba(:)
    ! Per-element distributed load per unit length (3, n_elem): stored so runtime
    ! segment-length updates (active line control) can rescale f_ext exactly.
    REAL(wp), INTENT(IN), OPTIONAL :: dist_load_per_length(:, :)
    ! Seabed recompute inputs (see the type note): kBot, per-element contact
    ! diameter, and the cn/kn ratio.
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_kbot, seabed_cn_over_kn
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_contact_diameter(:)
    REAL(wp), INTENT(IN), OPTIONAL :: fluid_velocity(:, :)
    REAL(wp), INTENT(IN), OPTIONAL :: drag_waterline_z(:), drag_rho, drag_diameter, drag_cdn, drag_cdt
    REAL(wp), INTENT(IN), OPTIONAL :: drag_diameter_elem(:), drag_cdn_elem(:), drag_cdt_elem(:)
    REAL(wp), INTENT(IN), OPTIONAL :: fluid_acceleration(:, :)
    REAL(wp), INTENT(IN), OPTIONAL :: fk_waterline_z(:), fk_rho, fk_diameter, fk_can, fk_cat
    REAL(wp), INTENT(IN), OPTIONAL :: fk_diameter_elem(:), fk_can_elem(:), fk_cat_elem(:)
    REAL(wp), INTENT(IN), OPTIONAL :: buoyancy_waterline_z(:), buoyancy_rho, buoyancy_diameter, buoyancy_gravity
    REAL(wp), INTENT(IN), OPTIONAL :: buoyancy_diameter_elem(:)
    REAL(wp), INTENT(IN), OPTIONAL :: added_mass_waterline_z(:)
    REAL(wp), INTENT(IN), OPTIONAL :: added_mass_rho, added_mass_diameter, added_mass_can, added_mass_cat
    REAL(wp), INTENT(IN), OPTIONAL :: added_mass_diameter_elem(:), added_mass_can_elem(:), added_mass_cat_elem(:)
    TYPE(CD_BathymetryType), INTENT(IN), OPTIONAL :: bathymetry
    ! Viscoelastic per-element parameters (MoorDyn ElasticMod 2/3): ve_ea_d > 0
    ! marks a constant-Ed element (mode 2); ve_alpha_mbl > 0 with ve_vbeta > 0
    ! marks a mean-load-dependent element (mode 3); zero everywhere = plain.
    ! One element cannot carry both. ve_ba is the static-branch dashpot Bs and
    ! ve_ba_d the dynamic-branch Bd; the first three arrays must be supplied
    ! together (the mode-3 pair together too). For viscoelastic elements the
    ! plain axial damping ba is REPLACED by the series-Kelvin law (ve_ba carries the Bs the
    ! MoorDyn BA column would have supplied).
    REAL(wp), INTENT(IN), OPTIONAL :: ve_ea_d(:), ve_ba(:), ve_ba_d(:)
    REAL(wp), INTENT(IN), OPTIONAL :: ve_alpha_mbl(:), ve_vbeta(:)
    ! Syrope per-element data: syrope_is(e) marks a Syrope element,
    ! syrope_type(e) its (already-initialized) constitutive data, syrope_slow0
    ! the committed slow-strain IC (the caller sets it -- for a fixed point,
    ! the steady partition of the initial stretch), and syrope_tmax0 the initial
    ! running-maximum mean tension. All four are supplied together.
    LOGICAL, INTENT(IN), OPTIONAL :: syrope_is(:)
    TYPE(CD_SyropeType), INTENT(IN), OPTIONAL :: syrope_type(:)
    REAL(wp), INTENT(IN), OPTIONAL :: syrope_slow0(:), syrope_tmax0(:)

    INTEGER :: n_dof, n_elem, n_nodes, es, istat, e
    INTEGER, ALLOCATABLE :: free(:)
    TYPE(CD_ModelType) :: old_model
    CHARACTER(160) :: em
    LOGICAL :: use_seabed, use_seabed_damping, use_seabed_friction
    LOGICAL :: use_damping, use_drag, use_fk, use_buoyancy, use_added_mass, has_load
    LOGICAL :: use_viscoelastic, use_syrope
    LOGICAL :: scalar_drag, elem_drag, scalar_fk, elem_fk, scalar_buoyancy, elem_buoyancy
    LOGICAL :: scalar_added_mass, elem_added_mass
    REAL(wp) :: syrope_ws(CD_SYROPE_NWC), syrope_wt(CD_SYROPE_NWC), syrope_wsl(CD_SYROPE_NWC)

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    ! Drop transient solver scratch before the rollback snapshot. The dense load-contributor work
    ! matrices and the gen-alpha workspace are rebuilt on demand and are intentionally excluded from
    ! CD_Copy_Model; letting the intrinsic old_model = model snapshot deep-copy them would duplicate
    ! n_dof x n_dof arrays through an unchecked assignment on every reinitialization (and can abort
    ! under memory pressure before any checked allocation runs).
    CALL clear_model_load_contributor_workspace(model)
    CALL CD_Clear_GenAlpha_Workspace(model%dynamic_workspace)
    old_model = model

    n_dof = SIZE(q0)
    n_elem = SIZE(elem_conn, 2)
    n_nodes = n_dof/3
    IF (MOD(n_dof, 3) /= 0 .OR. n_dof < 6) THEN
      CALL fail(ErrStat, ErrMsg, 'q0 must be a positions-only state of shape (3 n_nodes), n_nodes >= 2')
      RETURN
    END IF
    IF (SIZE(v0) /= n_dof .OR. SIZE(f_ext) /= n_dof) THEN
      CALL fail(ErrStat, ErrMsg, 'q0, v0, and f_ext must share shape (3 n_nodes)')
      RETURN
    END IF
    IF (n_elem < 1 .OR. SIZE(elem_conn, 1) /= 2) THEN
      CALL fail(ErrStat, ErrMsg, 'elem_conn must have shape (2, n_elem), n_elem >= 1')
      RETURN
    END IF
    IF (SIZE(l0) /= n_elem .OR. SIZE(ea) /= n_elem .OR. SIZE(rho_a) /= n_elem) THEN
      CALL fail(ErrStat, ErrMsg, 'l0, ea, and rho_a must have shape (n_elem)')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(q0) .OR. .NOT. CD_All_Finite(v0) &
        .OR. .NOT. CD_All_Finite(f_ext)) THEN
      CALL fail(ErrStat, ErrMsg, 'q0, v0, and f_ext must be finite')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(l0) .OR. ANY(l0 <= CD_ZERO) &
        .OR. .NOT. CD_All_Finite(ea) .OR. ANY(ea <= CD_ZERO) &
        .OR. .NOT. CD_All_Finite(rho_a) .OR. ANY(rho_a <= CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'l0, ea, and rho_a must be finite and positive')
      RETURN
    END IF
    IF (ANY(elem_conn < 1) .OR. ANY(elem_conn > n_dof/3)) THEN
      CALL fail(ErrStat, ErrMsg, 'elem_conn contains a node index outside q0')
      RETURN
    END IF
    CALL CD_Partition_Free_Dofs(fixed_dofs, n_dof, free, es, em)
    IF (es /= 0) THEN
      CALL fail(ErrStat, ErrMsg, 'fixed_dofs rejected: '//TRIM(em))
      RETURN
    END IF
    IF (PRESENT(seabed_kn) .NEQV. (PRESENT(seabed_z_floor) .OR. PRESENT(bathymetry))) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed requires seabed_kn with either seabed_z_floor or bathymetry')
      RETURN
    END IF
    IF (PRESENT(seabed_z_floor) .AND. PRESENT(bathymetry)) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed_z_floor and bathymetry are mutually exclusive')
      RETURN
    END IF
    use_seabed = PRESENT(seabed_kn)
    use_seabed_damping = PRESENT(seabed_cn)
    use_seabed_friction = PRESENT(seabed_mu)
    use_damping = PRESENT(ba)
    use_drag = PRESENT(fluid_velocity) .OR. PRESENT(drag_waterline_z) .OR. PRESENT(drag_rho) .OR. &
               PRESENT(drag_diameter) .OR. PRESENT(drag_cdn) .OR. PRESENT(drag_cdt) .OR. &
               PRESENT(drag_diameter_elem) .OR. PRESENT(drag_cdn_elem) .OR. PRESENT(drag_cdt_elem)
    use_fk = PRESENT(fluid_acceleration) .OR. PRESENT(fk_waterline_z) .OR. PRESENT(fk_rho) .OR. &
             PRESENT(fk_diameter) .OR. PRESENT(fk_can) .OR. PRESENT(fk_cat) .OR. &
             PRESENT(fk_diameter_elem) .OR. PRESENT(fk_can_elem) .OR. PRESENT(fk_cat_elem)
    use_buoyancy = PRESENT(buoyancy_waterline_z) .OR. PRESENT(buoyancy_rho) .OR. &
                   PRESENT(buoyancy_diameter) .OR. PRESENT(buoyancy_gravity) .OR. &
                   PRESENT(buoyancy_diameter_elem)
    use_added_mass = PRESENT(added_mass_waterline_z) .OR. PRESENT(added_mass_rho) .OR. &
                     PRESENT(added_mass_diameter) .OR. PRESENT(added_mass_can) .OR. &
                     PRESENT(added_mass_cat) .OR. PRESENT(added_mass_diameter_elem) .OR. &
                     PRESENT(added_mass_can_elem) .OR. PRESENT(added_mass_cat_elem)
    scalar_drag = .FALSE.
    elem_drag = .FALSE.
    scalar_fk = .FALSE.
    elem_fk = .FALSE.
    scalar_buoyancy = .FALSE.
    elem_buoyancy = .FALSE.
    scalar_added_mass = .FALSE.
    elem_added_mass = .FALSE.
    IF (use_seabed) THEN
      IF (PRESENT(seabed_z_floor)) THEN
        IF (.NOT. CD_Is_Finite(seabed_z_floor)) THEN
          CALL fail(ErrStat, ErrMsg, 'seabed_z_floor must be finite')
          RETURN
        END IF
      ELSE
        IF (.NOT. CD_Bathymetry_Is_Initialized(bathymetry)) THEN
          CALL fail(ErrStat, ErrMsg, 'bathymetry seabed requires an initialized bathymetry object')
          RETURN
        END IF
      END IF
      IF (SIZE(seabed_kn) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'seabed_kn must have shape (n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(seabed_kn) .OR. ANY(seabed_kn <= CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'seabed_kn must be finite and positive')
        RETURN
      END IF
    END IF
    IF (use_seabed_damping) THEN
      IF (.NOT. use_seabed) THEN
        CALL fail(ErrStat, ErrMsg, 'seabed_cn requires seabed_z_floor and seabed_kn')
        RETURN
      END IF
      IF (SIZE(seabed_cn) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'seabed_cn must have shape (n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(seabed_cn) .OR. ANY(seabed_cn < CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'seabed_cn must be finite and non-negative')
        RETURN
      END IF
    END IF
    IF (use_seabed_friction) THEN
      IF (.NOT. use_seabed) THEN
        CALL fail(ErrStat, ErrMsg, 'seabed_mu requires seabed_z_floor and seabed_kn')
        RETURN
      END IF
      IF (.NOT. CD_Is_Finite(seabed_mu) .OR. seabed_mu < CD_ZERO) THEN
        CALL fail(ErrStat, ErrMsg, 'seabed_mu must be finite and non-negative')
        RETURN
      END IF
      use_seabed_friction = seabed_mu > CD_ZERO
    END IF
    IF (use_damping) THEN
      IF (SIZE(ba) /= n_elem) THEN
        CALL fail(ErrStat, ErrMsg, 'ba must have shape (n_elem)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(ba) .OR. ANY(ba < CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'ba must be finite and non-negative')
        RETURN
      END IF
    END IF
    IF (use_drag) THEN
      scalar_drag = PRESENT(drag_diameter) .AND. PRESENT(drag_cdn) .AND. PRESENT(drag_cdt)
      elem_drag = PRESENT(drag_diameter_elem) .AND. PRESENT(drag_cdn_elem) .AND. PRESENT(drag_cdt_elem)
      IF (.NOT. (PRESENT(fluid_velocity) .AND. PRESENT(drag_waterline_z) .AND. PRESENT(drag_rho) .AND. &
                 (scalar_drag .OR. elem_drag))) THEN
        CALL fail(ErrStat, ErrMsg, &
                  'Morison drag requires fluid_velocity, waterline_z, rho, and scalar or per-element D/Cdn/Cdt')
        RETURN
      END IF
      IF (scalar_drag .AND. elem_drag) THEN
        CALL fail(ErrStat, ErrMsg, 'Morison drag accepts either scalar or per-element coefficients, not both')
        RETURN
      END IF
      IF (SIZE(fluid_velocity, 1) /= 3 .OR. SIZE(fluid_velocity, 2) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_velocity must have shape (3, n_nodes)')
        RETURN
      END IF
      IF (SIZE(drag_waterline_z) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'drag_waterline_z must have shape (n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(fluid_velocity) .OR. .NOT. CD_All_Finite(drag_waterline_z)) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_velocity and drag_waterline_z must be finite')
        RETURN
      END IF
      IF (.NOT. (CD_Is_Finite(drag_rho) .AND. drag_rho > CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'Morison drag rho must be finite and positive')
        RETURN
      END IF
      IF (scalar_drag) THEN
        IF (.NOT. (CD_Is_Finite(drag_diameter) .AND. drag_diameter > CD_ZERO .AND. &
                   CD_Is_Finite(drag_cdn) .AND. drag_cdn >= CD_ZERO .AND. &
                   CD_Is_Finite(drag_cdt) .AND. drag_cdt >= CD_ZERO)) THEN
          CALL fail(ErrStat, ErrMsg, &
                    'Morison drag scalars must be finite with positive diameter and non-negative Cd')
          RETURN
        END IF
      ELSE
        CALL validate_positive_vector(drag_diameter_elem, n_elem, 'drag_diameter_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL validate_nonnegative_vector(drag_cdn_elem, n_elem, 'drag_cdn_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL validate_nonnegative_vector(drag_cdt_elem, n_elem, 'drag_cdt_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END IF
    END IF
    IF (use_fk) THEN
      scalar_fk = PRESENT(fk_diameter) .AND. PRESENT(fk_can) .AND. PRESENT(fk_cat)
      elem_fk = PRESENT(fk_diameter_elem) .AND. PRESENT(fk_can_elem) .AND. PRESENT(fk_cat_elem)
      IF (.NOT. (PRESENT(fluid_acceleration) .AND. PRESENT(fk_waterline_z) .AND. PRESENT(fk_rho) .AND. &
                 (scalar_fk .OR. elem_fk))) THEN
        CALL fail(ErrStat, ErrMsg, &
                  'Froude-Krylov requires fluid_acceleration, waterline_z, rho, and scalar or per-element D/Can/Cat')
        RETURN
      END IF
      IF (scalar_fk .AND. elem_fk) THEN
        CALL fail(ErrStat, ErrMsg, 'Froude-Krylov accepts either scalar or per-element coefficients, not both')
        RETURN
      END IF
      IF (SIZE(fluid_acceleration, 1) /= 3 .OR. SIZE(fluid_acceleration, 2) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_acceleration must have shape (3, n_nodes)')
        RETURN
      END IF
      IF (SIZE(fk_waterline_z) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'fk_waterline_z must have shape (n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(fluid_acceleration) .OR. .NOT. CD_All_Finite(fk_waterline_z)) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_acceleration and fk_waterline_z must be finite')
        RETURN
      END IF
      IF (.NOT. (CD_Is_Finite(fk_rho) .AND. fk_rho > CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'Froude-Krylov rho must be finite and positive')
        RETURN
      END IF
      IF (scalar_fk) THEN
        IF (.NOT. (CD_Is_Finite(fk_diameter) .AND. fk_diameter > CD_ZERO .AND. &
                   CD_Is_Finite(fk_can) .AND. fk_can >= CD_ZERO .AND. &
                   CD_Is_Finite(fk_cat) .AND. fk_cat >= CD_ZERO)) THEN
          CALL fail(ErrStat, ErrMsg, &
                    'Froude-Krylov scalars must be finite with positive diameter and non-negative Ca')
          RETURN
        END IF
      ELSE
        CALL validate_positive_vector(fk_diameter_elem, n_elem, 'fk_diameter_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL validate_nonnegative_vector(fk_can_elem, n_elem, 'fk_can_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL validate_nonnegative_vector(fk_cat_elem, n_elem, 'fk_cat_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END IF
    END IF
    IF (use_buoyancy) THEN
      scalar_buoyancy = PRESENT(buoyancy_diameter)
      elem_buoyancy = PRESENT(buoyancy_diameter_elem)
      IF (.NOT. (PRESENT(buoyancy_waterline_z) .AND. PRESENT(buoyancy_rho) .AND. &
                 PRESENT(buoyancy_gravity) .AND. (scalar_buoyancy .OR. elem_buoyancy))) THEN
        CALL fail(ErrStat, ErrMsg, &
                  'buoyancy recovery requires waterline_z, rho, gravity, and scalar or per-element diameter')
        RETURN
      END IF
      IF (scalar_buoyancy .AND. elem_buoyancy) THEN
        CALL fail(ErrStat, ErrMsg, 'buoyancy recovery accepts either scalar or per-element diameter, not both')
        RETURN
      END IF
      IF (SIZE(buoyancy_waterline_z) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'buoyancy_waterline_z must have shape (n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(buoyancy_waterline_z)) THEN
        CALL fail(ErrStat, ErrMsg, 'buoyancy_waterline_z must be finite')
        RETURN
      END IF
      IF (.NOT. (CD_Is_Finite(buoyancy_rho) .AND. buoyancy_rho > CD_ZERO .AND. &
                 CD_Is_Finite(buoyancy_gravity) .AND. buoyancy_gravity > CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'buoyancy recovery rho and gravity must be finite and positive')
        RETURN
      END IF
      IF (scalar_buoyancy) THEN
        IF (.NOT. (CD_Is_Finite(buoyancy_diameter) .AND. buoyancy_diameter > CD_ZERO)) THEN
          CALL fail(ErrStat, ErrMsg, 'buoyancy recovery diameter must be finite and positive')
          RETURN
        END IF
      ELSE
        CALL validate_positive_vector(buoyancy_diameter_elem, n_elem, 'buoyancy_diameter_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END IF
    END IF
    IF (use_added_mass) THEN
      scalar_added_mass = PRESENT(added_mass_diameter) .AND. PRESENT(added_mass_can) .AND. PRESENT(added_mass_cat)
      elem_added_mass = PRESENT(added_mass_diameter_elem) .AND. PRESENT(added_mass_can_elem) .AND. &
                        PRESENT(added_mass_cat_elem)
      IF (.NOT. (PRESENT(added_mass_waterline_z) .AND. PRESENT(added_mass_rho) .AND. &
                 (scalar_added_mass .OR. elem_added_mass))) THEN
        CALL fail(ErrStat, ErrMsg, &
                  'added mass requires waterline_z, rho, and scalar or per-element D/Can/Cat')
        RETURN
      END IF
      IF (scalar_added_mass .AND. elem_added_mass) THEN
        CALL fail(ErrStat, ErrMsg, 'added mass accepts either scalar or per-element coefficients, not both')
        RETURN
      END IF
      IF (SIZE(added_mass_waterline_z) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'added_mass_waterline_z must have shape (n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(added_mass_waterline_z)) THEN
        CALL fail(ErrStat, ErrMsg, 'added_mass_waterline_z must be finite')
        RETURN
      END IF
      IF (.NOT. (CD_Is_Finite(added_mass_rho) .AND. added_mass_rho > CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'added-mass rho must be finite and positive')
        RETURN
      END IF
      IF (scalar_added_mass) THEN
        IF (.NOT. (CD_Is_Finite(added_mass_diameter) .AND. added_mass_diameter > CD_ZERO .AND. &
                   CD_Is_Finite(added_mass_can) .AND. added_mass_can >= CD_ZERO .AND. &
                   CD_Is_Finite(added_mass_cat) .AND. added_mass_cat >= CD_ZERO)) THEN
          CALL fail(ErrStat, ErrMsg, &
                    'added-mass scalars must be finite with positive diameter and non-negative Ca')
          RETURN
        END IF
      ELSE
        CALL validate_positive_vector(added_mass_diameter_elem, n_elem, 'added_mass_diameter_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL validate_nonnegative_vector(added_mass_can_elem, n_elem, 'added_mass_can_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL validate_nonnegative_vector(added_mass_cat_elem, n_elem, 'added_mass_cat_elem', ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END IF
    END IF

    use_viscoelastic = PRESENT(ve_ea_d) .OR. PRESENT(ve_ba) .OR. PRESENT(ve_ba_d) .OR. &
                       PRESENT(ve_alpha_mbl) .OR. PRESENT(ve_vbeta)
    IF (use_viscoelastic) THEN
      IF (.NOT. (PRESENT(ve_ea_d) .AND. PRESENT(ve_ba) .AND. PRESENT(ve_ba_d))) THEN
        CALL fail(ErrStat, ErrMsg, 'viscoelastic parameters need all of ve_ea_d, ve_ba, ve_ba_d')
        RETURN
      END IF
      IF (PRESENT(ve_alpha_mbl) .NEQV. PRESENT(ve_vbeta)) THEN
        CALL fail(ErrStat, ErrMsg, 'the load-dependent pair ve_alpha_mbl and ve_vbeta must be '// &
                  'supplied together')
        RETURN
      END IF
      IF (SIZE(ve_ea_d) /= n_elem .OR. SIZE(ve_ba) /= n_elem .OR. SIZE(ve_ba_d) /= n_elem) THEN
        CALL fail(ErrStat, ErrMsg, 've_ea_d, ve_ba, and ve_ba_d must have shape (n_elem)')
        RETURN
      END IF
      IF (.NOT. (CD_All_Finite(ve_ea_d) .AND. CD_All_Finite(ve_ba) .AND. &
                 CD_All_Finite(ve_ba_d))) THEN
        CALL fail(ErrStat, ErrMsg, 'viscoelastic parameters must be finite')
        RETURN
      END IF
      IF (ANY(ve_ea_d < CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 've_ea_d entries must be zero (plain element) or positive')
        RETURN
      END IF
      IF (PRESENT(ve_alpha_mbl)) THEN
        IF (SIZE(ve_alpha_mbl) /= n_elem .OR. SIZE(ve_vbeta) /= n_elem) THEN
          CALL fail(ErrStat, ErrMsg, 've_alpha_mbl and ve_vbeta must have shape (n_elem)')
          RETURN
        END IF
        IF (.NOT. (CD_All_Finite(ve_alpha_mbl) .AND. CD_All_Finite(ve_vbeta))) THEN
          CALL fail(ErrStat, ErrMsg, 'viscoelastic load-dependent parameters must be finite')
          RETURN
        END IF
        IF (ANY(ve_alpha_mbl < CD_ZERO) .OR. ANY(ve_vbeta < CD_ZERO)) THEN
          CALL fail(ErrStat, ErrMsg, 've_alpha_mbl/ve_vbeta entries must be zero (not load-'// &
                    'dependent) or positive')
          RETURN
        END IF
        DO e = 1, n_elem
          IF ((ve_alpha_mbl(e) > CD_ZERO) .NEQV. (ve_vbeta(e) > CD_ZERO)) THEN
            CALL fail(ErrStat, ErrMsg, 'a load-dependent element needs BOTH ve_alpha_mbl > 0 '// &
                      'and ve_vbeta > 0')
            RETURN
          END IF
          IF (ve_alpha_mbl(e) > CD_ZERO .AND. ve_ea_d(e) > CD_ZERO) THEN
            CALL fail(ErrStat, ErrMsg, 'an element cannot carry both a constant dynamic stiffness '// &
                      '(ve_ea_d) and the load-dependent pair (ve_alpha_mbl/ve_vbeta)')
            RETURN
          END IF
        END DO
      END IF
      DO e = 1, n_elem
        IF (ve_ea_d(e) <= CD_ZERO .AND. .NOT. mode3_marked(ve_alpha_mbl, e) .AND. &
            ABS(ve_ba(e)) + ABS(ve_ba_d(e)) > CD_ZERO) THEN
          CALL fail(ErrStat, ErrMsg, 'a plain element (no ve_ea_d, no ve_alpha_mbl) cannot carry '// &
                    'series-Kelvin dashpots (ve_ba/ve_ba_d must be zero there)')
          RETURN
        END IF
      END DO
      ! all-zero arrays mean no viscoelastic elements at all
      use_viscoelastic = ANY(ve_ea_d > CD_ZERO)
      IF (PRESENT(ve_alpha_mbl)) use_viscoelastic = use_viscoelastic .OR. ANY(ve_alpha_mbl > CD_ZERO)
    END IF

    ! Syrope: the four arrays come together; syrope_type entries flagged by
    ! syrope_is must be initialized; a Syrope element cannot also be viscoelastic.
    use_syrope = PRESENT(syrope_is) .OR. PRESENT(syrope_type) .OR. PRESENT(syrope_slow0) .OR. &
                 PRESENT(syrope_tmax0)
    IF (use_syrope) THEN
      IF (.NOT. (PRESENT(syrope_is) .AND. PRESENT(syrope_type) .AND. PRESENT(syrope_slow0) .AND. &
                 PRESENT(syrope_tmax0))) THEN
        CALL fail(ErrStat, ErrMsg, 'Syrope needs all of syrope_is, syrope_type, syrope_slow0, syrope_tmax0')
        RETURN
      END IF
      IF (SIZE(syrope_is) /= n_elem .OR. SIZE(syrope_type) /= n_elem .OR. &
          SIZE(syrope_slow0) /= n_elem .OR. SIZE(syrope_tmax0) /= n_elem) THEN
        CALL fail(ErrStat, ErrMsg, 'the Syrope arrays must have shape (n_elem)')
        RETURN
      END IF
      IF (.NOT. (CD_All_Finite(syrope_slow0) .AND. CD_All_Finite(syrope_tmax0))) THEN
        CALL fail(ErrStat, ErrMsg, 'Syrope IC (slow0, tmax0) must be finite')
        RETURN
      END IF
      DO e = 1, n_elem
        IF (.NOT. syrope_is(e)) CYCLE
        IF (syrope_slow0(e) < CD_ZERO) THEN
          CALL fail(ErrStat, ErrMsg, 'Syrope slow0 must be non-negative on a Syrope element')
          RETURN
        END IF
        IF (.NOT. CD_Syrope_Is_Ready(syrope_type(e))) THEN
          CALL fail(ErrStat, ErrMsg, 'a Syrope element was flagged without initialized constitutive data')
          RETURN
        END IF
        IF (.NOT. (syrope_tmax0(e) > CD_ZERO)) THEN
          CALL fail(ErrStat, ErrMsg, 'Syrope tmax0 must be positive on a Syrope element')
          RETURN
        END IF
        CALL CD_Syrope_Working_Curve(syrope_type(e), syrope_tmax0(e), syrope_ws, syrope_wt, syrope_wsl, es, em)
        IF (es /= CD_SYROPE_OK) THEN
          CALL fail(ErrStat, ErrMsg, 'Syrope tmax0 is invalid: '//TRIM(em))
          RETURN
        END IF
        CALL CD_Syrope_Check_Range(syrope_type(e), syrope_slow0(e), es, em)
        IF (es /= CD_SYROPE_OK) THEN
          CALL fail(ErrStat, ErrMsg, 'Syrope slow0 is invalid: '//TRIM(em))
          RETURN
        END IF
        IF (use_viscoelastic) THEN
          IF (ve_ea_d(e) > CD_ZERO .OR. mode3_marked(ve_alpha_mbl, e)) THEN
            CALL fail(ErrStat, ErrMsg, 'an element cannot be both Syrope and viscoelastic')
            RETURN
          END IF
        END IF
      END DO
      use_syrope = ANY(syrope_is)
    END IF

    CALL CD_End_Model(model, es, em)
    model%n_dof = n_dof
    model%n_elem = n_elem
    model%tension_only = tension_only
    model%cfg = cfg
    ALLOCATE (model%elem_conn(2, n_elem), model%fixed_dofs(SIZE(fixed_dofs)), &
              model%l0(n_elem), model%ea(n_elem), model%ea_dyn(n_elem), model%rho_a(n_elem), model%l0_dot(n_elem), &
              model%f_ext(n_dof), &
              model%q(n_dof), model%v(n_dof), model%a(n_dof), &
              model%q_step(n_dof), model%v_step(n_dof), model%a_step(n_dof), &
              model%q_prescribed_work(n_dof), model%v_prescribed_work(n_dof), model%a_prescribed_work(n_dof), &
              model%q_rollback(n_dof), model%v_rollback(n_dof), model%a_rollback(n_dof), &
              model%recovery_q0(n_dof), model%recovery_v0(n_dof), model%recovery_a0(n_dof), &
              model%f_ext_rollback(n_dof), model%free_map_work(n_dof), STAT=istat)
    IF (istat /= 0) THEN
      CALL restore_after_init_alloc_fail(model, old_model, 'base model storage', ErrStat, ErrMsg)
      RETURN
    END IF
    model%elem_conn = elem_conn
    model%fixed_dofs = fixed_dofs
    model%l0 = l0
    model%l0_dot = CD_ZERO
    model%ea = ea
    model%ea_dyn = ea
    model%rho_a = rho_a
    model%f_ext = f_ext
    IF (PRESENT(dist_load_per_length)) THEN
      IF (SIZE(dist_load_per_length, 1) /= 3 .OR. SIZE(dist_load_per_length, 2) /= n_elem .OR. &
          .NOT. CD_All_Finite(dist_load_per_length)) THEN
        ! Restore the previous model, then report the bad input as such.
        CALL restore_after_init_alloc_fail(model, old_model, 'distributed-load shape', ErrStat, ErrMsg)
        ErrStat = CD_MODEL_BADINPUT
        ErrMsg = 'CableDyn_Model: dist_load_per_length must be finite with shape (3, n_elem)'
        RETURN
      END IF
      ALLOCATE (model%dist_load(3, n_elem), model%f_ext_dist_base(n_dof), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'distributed-load storage', ErrStat, ErrMsg)
        RETURN
      END IF
      model%dist_load = dist_load_per_length
      ! the non-distributed remainder: for the deck path (f_ext IS the assembled
      ! distributed load) this is exactly zero, so later reassembly reproduces a
      ! fresh build bit-for-bit
      CALL CD_Assemble_Distributed_Load(elem_conn, l0, dist_load_per_length, model%f_ext_dist_base, es, em)
      IF (es /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'distributed-load remainder', ErrStat, ErrMsg)
        ErrStat = CD_MODEL_BADINPUT
        ErrMsg = 'CableDyn_Model: distributed-load assembly failed: '//TRIM(em)
        RETURN
      END IF
      model%f_ext_dist_base = model%f_ext - model%f_ext_dist_base
    END IF
    IF (PRESENT(seabed_kbot) .AND. PRESENT(seabed_contact_diameter)) THEN
      BLOCK
        LOGICAL :: ratio_ok
        ! With seabed DAMPING active the recompute rebuilds cn from the ratio; a
        ! missing ratio would silently wipe the damping to zero on the first
        ! length update, so recompute stays DISABLED without it (the update then
        ! fails closed by name instead). Without damping the ratio is irrelevant.
        ratio_ok = .NOT. use_seabed_damping
        IF (PRESENT(seabed_cn_over_kn)) THEN
          IF (CD_Is_Finite(seabed_cn_over_kn) .AND. seabed_cn_over_kn >= CD_ZERO) ratio_ok = .TRUE.
        END IF
        IF (CD_Is_Finite(seabed_kbot) .AND. seabed_kbot > CD_ZERO .AND. &
            SIZE(seabed_contact_diameter) == n_elem .AND. &
            CD_All_Finite(seabed_contact_diameter) .AND. ALL(seabed_contact_diameter > CD_ZERO) .AND. &
            ratio_ok) THEN
          ALLOCATE (model%seabed_diameter_elem(n_elem), STAT=istat)
          IF (istat /= 0) THEN
            CALL restore_after_init_alloc_fail(model, old_model, 'seabed recompute storage', ErrStat, ErrMsg)
            RETURN
          END IF
          model%seabed_diameter_elem = seabed_contact_diameter
          model%has_seabed_recompute = .TRUE.
          model%seabed_kbot = seabed_kbot
          model%seabed_cn_over_kn = CD_ZERO
          IF (PRESENT(seabed_cn_over_kn)) THEN
            IF (CD_Is_Finite(seabed_cn_over_kn) .AND. seabed_cn_over_kn >= CD_ZERO) &
              model%seabed_cn_over_kn = seabed_cn_over_kn
          END IF
        END IF
      END BLOCK
    END IF
    model%q = q0
    model%v = v0
    model%a = CD_ZERO
    model%q_step = CD_ZERO
    model%v_step = CD_ZERO
    model%a_step = CD_ZERO
    model%q_prescribed_work = CD_ZERO
    model%v_prescribed_work = CD_ZERO
    model%a_prescribed_work = CD_ZERO
    model%q_rollback = CD_ZERO
    model%v_rollback = CD_ZERO
    model%a_rollback = CD_ZERO
    model%recovery_q0 = CD_ZERO
    model%recovery_v0 = CD_ZERO
    model%recovery_a0 = CD_ZERO
    model%f_ext_rollback = CD_ZERO
    model%free_map_work = 0
    ! The consistent mass is never stored densely: the step path applies it element by
    ! element (band assembly and matvec) and the coupled-derivative queries assemble it
    ! on demand. Validate the connectivity it is built from here, as the stored mass
    ! assembly did.
    ! validated once here: the connectivity never changes afterwards, so the model's own
    ! assembly calls pass topology_validated
    CALL CD_Validate_Connectivity(model%elem_conn, model%n_dof/3, model%n_elem, es, em)
    IF (es /= 0) THEN
      CALL CD_End_Model(model, ErrStat, ErrMsg)
      model = old_model
      ErrStat = CD_MODEL_SOLVEFAIL
      ErrMsg = 'CableDyn_Model: mesh connectivity rejected: '//TRIM(em)
      RETURN
    END IF
    model%has_seabed = use_seabed
    model%has_bathymetry = use_seabed .AND. PRESENT(bathymetry)
    model%has_seabed_damping = use_seabed_damping
    model%has_seabed_friction = use_seabed_friction
    model%has_damping = use_damping
    model%has_viscoelastic = use_viscoelastic
    model%has_syrope = use_syrope
    model%has_morison_drag = use_drag
    model%has_froude_krylov = use_fk
    model%has_buoyancy_recovery = use_buoyancy
    model%has_added_mass = use_added_mass
    IF (use_seabed) THEN
      IF (PRESENT(seabed_z_floor)) model%seabed_z_floor = seabed_z_floor
      IF (PRESENT(bathymetry)) model%bathymetry = bathymetry
      ALLOCATE (model%seabed_kn(n_nodes), source=seabed_kn, STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'seabed stiffness storage', ErrStat, ErrMsg)
        RETURN
      END IF
    END IF
    IF (use_seabed_damping) THEN
      ALLOCATE (model%seabed_cn(n_nodes), source=seabed_cn, STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'seabed damping storage', ErrStat, ErrMsg)
        RETURN
      END IF
    END IF
    IF (use_seabed_friction) THEN
      model%seabed_mu = seabed_mu
      ! Friction springs anchored at the initial positions (no force); a static solve with
      ! friction installs its spring forces through CD_Set_Model_Friction_Anchors.
      ALLOCATE (model%fr_anchor(2, n_nodes), model%fr_anchor_rollback(2, n_nodes), &
                model%recovery_fr_anchor(2, n_nodes), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'seabed friction storage', ErrStat, ErrMsg)
        RETURN
      END IF
      DO e = 1, n_nodes
        model%fr_anchor(:, e) = q0(3*e - 2:3*e - 1)
      END DO
      model%fr_anchor_rollback = model%fr_anchor
      model%recovery_fr_anchor = model%fr_anchor
    END IF
    IF (use_damping) THEN
      ALLOCATE (model%ba(n_elem), source=ba, STAT=istat)
      IF (istat == 0) ALLOCATE (model%ba_dyn(n_elem), source=ba, STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'axial damping storage', ErrStat, ErrMsg)
        RETURN
      END IF
    END IF
    IF (use_viscoelastic) THEN
      ALLOCATE (model%ve_mode(n_elem), model%ve_ea_d(n_elem), model%ve_ea_1(n_elem), &
                model%ve_alpha_mbl(n_elem), model%ve_vbeta(n_elem), model%ve_ba(n_elem), &
                model%ve_ba_d(n_elem), model%ve_dl_1(n_elem), model%ve_dl_1_rollback(n_elem), &
                model%recovery_ve_dl_1(n_elem), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'viscoelastic storage', ErrStat, ErrMsg)
        RETURN
      END IF
      model%ve_mode = 0
      model%ve_ea_d = ve_ea_d
      model%ve_ba = ve_ba
      model%ve_ba_d = ve_ba_d
      model%ve_ea_1 = CD_ZERO
      model%ve_alpha_mbl = CD_ZERO
      model%ve_vbeta = CD_ZERO
      IF (PRESENT(ve_alpha_mbl)) THEN
        model%ve_alpha_mbl = ve_alpha_mbl
        model%ve_vbeta = ve_vbeta
      END IF
      model%ve_dl_1 = CD_ZERO
      DO e = 1, n_elem
        IF (ve_ea_d(e) > CD_ZERO) THEN
          model%ve_mode(e) = 2
        ELSE IF (mode3_marked(ve_alpha_mbl, e)) THEN
          model%ve_mode(e) = 3
        ELSE
          CYCLE
        END IF
        ! IC: dl_1 starts at the STEADY-STATE partition of the initial stretch
        ! (the ld_1 = 0 point), at which the series-Kelvin tension equals the static
        ! composite EA*dl/l0 -- the static equilibrium the caller solved IS the
        ! dynamic fixed point. MoorDyn instead seeds dl_1 = lstr - l and
        ! relaxes it to this same partition over its TmaxIC window before
        ! t = 0, so the two codes enter production time in the same state.
        BLOCK
          INTEGER :: na, nb, it
          REAL(wp) :: lstr0, dl0, ea_d_eff, dl1, dl1_new
          na = elem_conn(1, e)
          nb = elem_conn(2, e)
          lstr0 = NORM2(q0(3*nb - 2:3*nb) - q0(3*na - 2:3*na))
          dl0 = lstr0 - l0(e)
          es = CD_VISCO_OK
          IF (model%ve_mode(e) == 2) THEN
            ea_d_eff = ve_ea_d(e)
            CALL CD_Viscoelastic_Steady_State(ea(e), ea_d_eff, dl0, model%ve_dl_1(e))
          ELSE
            ! mode 3: the partition is self-consistent (EA_D depends on dl_1);
            ! fixed-point iteration dl_1 <- dl0*(EA_D(dl_1) - EA)/EA_D(dl_1) --
            ! EA_D is monotone and bounded in dl_1, so the map contracts; fail
            ! closed if it does not settle
            dl1 = dl0
            ea_d_eff = CD_ZERO
            DO it = 1, 200
              CALL CD_Viscoelastic_LoadDependent_EAD(ea(e), l0(e), ve_alpha_mbl(e), ve_vbeta(e), &
                                                     dl1, ea_d_eff, es, em)
              IF (es /= CD_VISCO_OK) EXIT
              CALL CD_Viscoelastic_Steady_State(ea(e), ea_d_eff, dl0, dl1_new)
              IF (ABS(dl1_new - dl1) <= 1.0e-14_wp*MAX(1.0_wp, ABS(dl1_new))) THEN
                dl1 = dl1_new
                EXIT
              END IF
              dl1 = dl1_new
              IF (it == 200) THEN
                es = CD_VISCO_BADINPUT
                em = 'CableDyn_Viscoelastic: the load-dependent steady-state partition did not settle'
              END IF
            END DO
            model%ve_dl_1(e) = dl1
            ! EA_D(dl_1) is increasing in dl_1 with EA_D = alphaMBL at dl_1 <= 0, so
            ! the premise EA_D > EA holds over every admissible state (taut and
            ! slack) exactly when alphaMBL > EA; a check at the initial stretch
            ! alone would pass a taut IC and fail at the first slack segment.
            IF (.NOT. (ve_alpha_mbl(e) > ea(e))) THEN
              es = CD_VISCO_BADINPUT
              em = ''
              WRITE (em, '(A,I0,A,ES12.5,A,ES12.5,A)', IOSTAT=it) &
                'CableDyn_Viscoelastic: ElasticMod 3 element ', e, ' needs alphaMBL (', ve_alpha_mbl(e), &
                ' N) > static EA (', ea(e), ' N); a slack segment has EA_D = alphaMBL'
            END IF
          END IF
          ! the fatal battery (EA_D <= EA, negative dashpots, BA + BA_D = 0);
          ! mode 2 stores the derived EA_1, mode 3 derives both per evaluation
          ! from the committed dl_1 (the ElasticMod 3 premise is required above
          ! over the whole admissible dl_1 range, not only at the IC value)
          IF (es == CD_VISCO_OK) &
            CALL CD_Viscoelastic_Params(ea(e), ea_d_eff, ve_ba(e), ve_ba_d(e), model%ve_ea_1(e), es, em)
          IF (es /= CD_VISCO_OK) THEN
            CALL CD_End_Model(model, ErrStat, ErrMsg)
            model = old_model
            ErrStat = CD_MODEL_BADINPUT
            ErrMsg = 'CableDyn_Model: '//TRIM(em)
            RETURN
          END IF
          ! mode 3 has no stored EA_1 (state-dependent); keep only mode 2's
          IF (model%ve_mode(e) == 3) model%ve_ea_1(e) = CD_ZERO
        END BLOCK
        ! the series-Kelvin law replaces both the elastic EA and the plain axial damping of
        ! this element in the dynamic assembly
        model%ea_dyn(e) = CD_ZERO
        IF (ALLOCATED(model%ba_dyn)) model%ba_dyn(e) = CD_ZERO
      END DO
      model%ve_dl_1_rollback = model%ve_dl_1
      model%recovery_ve_dl_1 = model%ve_dl_1
    END IF
    IF (use_syrope) THEN
      ALLOCATE (model%syrope_is(n_elem), model%syrope_type(n_elem), model%syrope_slow(n_elem), &
                model%syrope_slow_rollback(n_elem), model%syrope_tmax(n_elem), &
                model%syrope_tmax_rollback(n_elem), model%recovery_syrope_slow(n_elem), &
                model%recovery_syrope_tmax(n_elem), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'Syrope storage', ErrStat, ErrMsg)
        RETURN
      END IF
      model%syrope_is = .FALSE.
      model%syrope_slow = CD_ZERO
      model%syrope_tmax = CD_ZERO
      DO e = 1, n_elem
        IF (.NOT. syrope_is(e)) CYCLE
        model%syrope_is(e) = .TRUE.
        ! deep-copy the (validated) constitutive data; the Syrope contributor
        ! supplies this element's tension, so its elastic EA and plain axial
        ! damping are masked out of the dynamic structural assembly
        model%syrope_type(e) = syrope_type(e)
        model%syrope_slow(e) = syrope_slow0(e)
        model%syrope_tmax(e) = syrope_tmax0(e)
        model%ea_dyn(e) = CD_ZERO
        IF (ALLOCATED(model%ba_dyn)) model%ba_dyn(e) = CD_ZERO
      END DO
      model%syrope_slow_rollback = model%syrope_slow
      model%syrope_tmax_rollback = model%syrope_tmax
      model%recovery_syrope_slow = model%syrope_slow
      model%recovery_syrope_tmax = model%syrope_tmax
    END IF
    IF (use_drag) THEN
      ALLOCATE (model%fluid_velocity(3, n_nodes), model%drag_waterline_z(n_nodes), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'drag field storage', ErrStat, ErrMsg)
        RETURN
      END IF
      model%fluid_velocity = fluid_velocity
      model%drag_waterline_z = drag_waterline_z
      model%drag_rho = drag_rho
      ALLOCATE (model%drag_diameter_elem(n_elem), model%drag_cdn_elem(n_elem), model%drag_cdt_elem(n_elem), &
                STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'drag coefficient storage', ErrStat, ErrMsg)
        RETURN
      END IF
      IF (scalar_drag) THEN
        model%drag_diameter = drag_diameter
        model%drag_cdn = drag_cdn
        model%drag_cdt = drag_cdt
        model%drag_diameter_elem = drag_diameter
        model%drag_cdn_elem = drag_cdn
        model%drag_cdt_elem = drag_cdt
      ELSE
        model%drag_diameter_elem = drag_diameter_elem
        model%drag_cdn_elem = drag_cdn_elem
        model%drag_cdt_elem = drag_cdt_elem
      END IF
    END IF
    IF (use_fk) THEN
      ALLOCATE (model%fluid_acceleration(3, n_nodes), model%fk_waterline_z(n_nodes), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'fluid-acceleration field storage', ErrStat, ErrMsg)
        RETURN
      END IF
      model%fluid_acceleration = fluid_acceleration
      model%fk_waterline_z = fk_waterline_z
      model%fk_rho = fk_rho
      ALLOCATE (model%fk_diameter_elem(n_elem), model%fk_can_elem(n_elem), model%fk_cat_elem(n_elem), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'fluid-inertia coefficient storage', ErrStat, ErrMsg)
        RETURN
      END IF
      IF (scalar_fk) THEN
        model%fk_diameter = fk_diameter
        model%fk_can = fk_can
        model%fk_cat = fk_cat
        model%fk_diameter_elem = fk_diameter
        model%fk_can_elem = fk_can
        model%fk_cat_elem = fk_cat
      ELSE
        model%fk_diameter_elem = fk_diameter_elem
        model%fk_can_elem = fk_can_elem
        model%fk_cat_elem = fk_cat_elem
      END IF
    END IF
    IF (use_buoyancy) THEN
      ALLOCATE (model%buoyancy_waterline_z(n_nodes), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'buoyancy waterline storage', ErrStat, ErrMsg)
        RETURN
      END IF
      model%buoyancy_waterline_z = buoyancy_waterline_z
      model%buoyancy_rho = buoyancy_rho
      model%buoyancy_gravity = buoyancy_gravity
      ALLOCATE (model%buoyancy_diameter_elem(n_elem), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'buoyancy diameter storage', ErrStat, ErrMsg)
        RETURN
      END IF
      IF (scalar_buoyancy) THEN
        model%buoyancy_diameter = buoyancy_diameter
        model%buoyancy_diameter_elem = buoyancy_diameter
      ELSE
        model%buoyancy_diameter_elem = buoyancy_diameter_elem
      END IF
    END IF
    IF (use_added_mass) THEN
      ALLOCATE (model%added_mass_waterline_z(n_nodes), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'added-mass waterline storage', ErrStat, ErrMsg)
        RETURN
      END IF
      model%added_mass_waterline_z = added_mass_waterline_z
      model%added_mass_rho = added_mass_rho
      ALLOCATE (model%added_mass_diameter_elem(n_elem), model%added_mass_can_elem(n_elem), &
                model%added_mass_cat_elem(n_elem), STAT=istat)
      IF (istat /= 0) THEN
        CALL restore_after_init_alloc_fail(model, old_model, 'added-mass coefficient storage', ErrStat, ErrMsg)
        RETURN
      END IF
      IF (scalar_added_mass) THEN
        model%added_mass_diameter = added_mass_diameter
        model%added_mass_can = added_mass_can
        model%added_mass_cat = added_mass_cat
        model%added_mass_diameter_elem = added_mass_diameter
        model%added_mass_can_elem = added_mass_can
        model%added_mass_cat_elem = added_mass_cat
      ELSE
        model%added_mass_diameter_elem = added_mass_diameter_elem
        model%added_mass_can_elem = added_mass_can_elem
        model%added_mass_cat_elem = added_mass_cat_elem
      END IF
    END IF

    has_load = model%has_seabed .OR. model%has_seabed_damping .OR. model%has_seabed_friction .OR. &
               model%has_damping .OR. model%has_viscoelastic .OR. model%has_syrope .OR. model%has_morison_drag .OR. &
               model%has_froude_krylov .OR. model%has_buoyancy_recovery
    IF (has_load .AND. model%has_added_mass) THEN
      CALL CD_Cable_Initial_Acceleration(model%q, model%v, model%elem_conn, model%l0, model%ea_dyn, &
                                         model%rho_a, model%tension_only, model%f_ext, &
                                         model%fixed_dofs, model%a, es, em, &
                                         load_force_proc=model_force, workspace=model%dynamic_workspace, &
                                         added_mass_band_proc=model_added_mass_band)
    ELSE IF (has_load) THEN
      CALL CD_Cable_Initial_Acceleration(model%q, model%v, model%elem_conn, model%l0, model%ea_dyn, &
                                         model%rho_a, model%tension_only, model%f_ext, &
                                         model%fixed_dofs, model%a, es, em, &
                                         load_force_proc=model_force, workspace=model%dynamic_workspace)
    ELSE IF (model%has_added_mass) THEN
      CALL CD_Cable_Initial_Acceleration(model%q, model%v, model%elem_conn, model%l0, model%ea_dyn, &
                                         model%rho_a, model%tension_only, model%f_ext, &
                                         model%fixed_dofs, model%a, es, em, &
                                         workspace=model%dynamic_workspace, added_mass_band_proc=model_added_mass_band)
    ELSE
      CALL CD_Cable_Initial_Acceleration(model%q, model%v, model%elem_conn, model%l0, model%ea_dyn, &
                                         model%rho_a, model%tension_only, model%f_ext, &
                                         model%fixed_dofs, model%a, es, em, workspace=model%dynamic_workspace)
    END IF
    IF (es /= CD_DYN_OK) THEN
      CALL CD_End_Model(model, ErrStat, ErrMsg)
      model = old_model
      ErrStat = model_status_from_dyn(es)
      ErrMsg = 'CableDyn_Model: initial acceleration failed: '//TRIM(em)
      RETURN
    END IF
    model%initialized = .TRUE.

  CONTAINS

    SUBROUTINE model_force(q, v, force, ErrStat, ErrMsg)
      REAL(wp), INTENT(IN) :: q(:), v(:)
      REAL(wp), INTENT(OUT) :: force(:)
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg
      CALL model_contributor_force(model, q, v, force, ErrStat, ErrMsg)
      ErrStat = as_dyn_callback_status(ErrStat)
    END SUBROUTINE model_force

    SUBROUTINE model_added_mass_band(q, accel, free, kl, ku, M_add_a, M_add_band, dMa_a_dq_band, ErrStat, ErrMsg, &
                                     need_tangent)
      REAL(wp), INTENT(IN) :: q(:), accel(:)
      INTEGER, INTENT(IN) :: free(:), kl, ku
      REAL(wp), INTENT(OUT) :: M_add_a(:), M_add_band(:, :), dMa_a_dq_band(:, :)
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg
      LOGICAL, INTENT(IN), OPTIONAL :: need_tangent
      CALL model_added_mass_banded(model, q, accel, free, kl, ku, M_add_a, M_add_band, dMa_a_dq_band, &
                                   ErrStat, ErrMsg, need_tangent=need_tangent)
      ErrStat = as_dyn_callback_status(ErrStat)
    END SUBROUTINE model_added_mass_band
  END SUBROUTINE CD_Init_Model

  SUBROUTINE CD_Init_Line_Model(model, anchor, fairlead, line_types, sections, gravity, rho_water, &
                                tension_only, static_cfg, dynamic_cfg, load_factors, ErrStat, ErrMsg, &
                                anchor_is_end_b, seabed_z_floor, seabed_kn, seabed_cn, seabed_mu, ba, &
                                fluid_velocity, drag_waterline_z, drag_rho, drag_diameter, drag_cdn, drag_cdt, &
                                drag_diameter_elem, drag_cdn_elem, drag_cdt_elem, &
                                fluid_acceleration, fk_waterline_z, fk_rho, fk_diameter, fk_can, fk_cat, &
                                fk_diameter_elem, fk_can_elem, fk_cat_elem, seabed_mu_axial, &
                                buoyancy_waterline_z, buoyancy_rho, buoyancy_diameter, buoyancy_gravity, &
                                buoyancy_diameter_elem, &
                                added_mass_waterline_z, added_mass_rho, added_mass_diameter, added_mass_can, &
                                added_mass_cat, added_mass_diameter_elem, added_mass_can_elem, added_mass_cat_elem, &
                                bathymetry, dynamic_tension_only, &
                                seabed_kbot, seabed_contact_diameter, seabed_cn_over_kn, &
                                ve_ea_d, ve_ba, ve_ba_d, ve_alpha_mbl, ve_vbeta, &
                                syrope_is, syrope_type, syrope_slow0, syrope_tmax0, &
                                current_profile_z, current_profile_velocity)
    !! Initialise a persistent EI=0 dynamic model from an OrcaFlex-style line object:
    !! line types + sections + end points. This routine builds the mesh, assembles
    !! submerged-weight loads, obtains a static initial condition by load continuation,
    !! and delegates the persistent-state population to CD_Init_Model.
    !!
    !! Optional seabed and axial-damping contributors are stored on the model and
    !! reused during initial acceleration, stepping, and coupled-load extraction so
    !! the static IC and the first dynamic step see the same load model.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: anchor(3), fairlead(3)
    TYPE(CD_LineType), INTENT(IN) :: line_types(:)
    TYPE(CD_LineSection), INTENT(IN) :: sections(:)
    REAL(wp), INTENT(IN) :: gravity, rho_water
    LOGICAL, INTENT(IN) :: tension_only
    ! dynamic_tension_only (default = tension_only): the COMMITTED dynamic model's
    ! compression convention, decoupled from the static IC's. The static continuation
    ! stays compression-capable where requested (a fully slack element has zero tangent
    ! stiffness, so tension-only STATICS on a slack span is singular), while the marched
    ! model can be tension-only -- MoorDyn's lumped-mass convention, the physical one
    ! for chains/ropes: a slack element exerts nothing and reports zero tension.
    TYPE(CableSolverConfig), INTENT(IN) :: static_cfg
    TYPE(GenAlphaConfig), INTENT(IN) :: dynamic_cfg
    REAL(wp), INTENT(IN) :: load_factors(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    LOGICAL, INTENT(IN), OPTIONAL :: anchor_is_end_b
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_z_floor
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_kn(:)
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_cn(:)
    ! Anisotropic seabed friction: the coefficient along the line (seabed_mu is then the
    ! lateral one); absent, or equal to seabed_mu, keeps the isotropic law.
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_mu_axial
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_mu
    REAL(wp), INTENT(IN), OPTIONAL :: ba(:)
    REAL(wp), INTENT(IN), OPTIONAL :: fluid_velocity(:, :)
    REAL(wp), INTENT(IN), OPTIONAL :: drag_waterline_z(:), drag_rho, drag_diameter, drag_cdn, drag_cdt
    REAL(wp), INTENT(IN), OPTIONAL :: drag_diameter_elem(:), drag_cdn_elem(:), drag_cdt_elem(:)
    REAL(wp), INTENT(IN), OPTIONAL :: fluid_acceleration(:, :)
    REAL(wp), INTENT(IN), OPTIONAL :: fk_waterline_z(:), fk_rho, fk_diameter, fk_can, fk_cat
    REAL(wp), INTENT(IN), OPTIONAL :: fk_diameter_elem(:), fk_can_elem(:), fk_cat_elem(:)
    REAL(wp), INTENT(IN), OPTIONAL :: buoyancy_waterline_z(:), buoyancy_rho, buoyancy_diameter, buoyancy_gravity
    REAL(wp), INTENT(IN), OPTIONAL :: buoyancy_diameter_elem(:)
    REAL(wp), INTENT(IN), OPTIONAL :: added_mass_waterline_z(:)
    REAL(wp), INTENT(IN), OPTIONAL :: added_mass_rho, added_mass_diameter, added_mass_can, added_mass_cat
    REAL(wp), INTENT(IN), OPTIONAL :: added_mass_diameter_elem(:), added_mass_can_elem(:), added_mass_cat_elem(:)
    TYPE(CD_BathymetryType), INTENT(IN), OPTIONAL :: bathymetry
    ! dynamic_tension_only sits LAST in the argument list: appending (never inserting)
    ! a new optional preserves every existing positional call's meaning (the API-
    ! stability rule for public module routines).
    LOGICAL, INTENT(IN), OPTIONAL :: dynamic_tension_only
    ! Seabed recompute inputs for active line control, passed through to
    ! CD_Init_Model (see the model type note).
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_kbot, seabed_cn_over_kn
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_contact_diameter(:)
    ! Viscoelastic per-element parameters, passed through to CD_Init_Model in
    ! the caller's FINAL element order (callers using anchor_is_end_b reversal
    ! must pre-reverse these alongside ba). The static IC below always solves
    ! with the composite ea (the series-Kelvin law relaxes to it), so only the dynamic model
    ! sees the split.
    REAL(wp), INTENT(IN), OPTIONAL :: ve_ea_d(:), ve_ba(:), ve_ba_d(:)
    REAL(wp), INTENT(IN), OPTIONAL :: ve_alpha_mbl(:), ve_vbeta(:)
    ! Syrope per-element data, passed through to CD_Init_Model in the caller's
    ! FINAL element order (pre-reversed alongside ba for anchor_is_end_b lines).
    ! The caller sets the composite ea to the OWC secant so the static IC solves
    ! at a strain consistent with the working-curve tension.
    LOGICAL, INTENT(IN), OPTIONAL :: syrope_is(:)
    TYPE(CD_SyropeType), INTENT(IN), OPTIONAL :: syrope_type(:)
    REAL(wp), INTENT(IN), OPTIONAL :: syrope_slow0(:), syrope_tmax0(:)
    ! Depth profile of the steady current (levels, velocities): the static initial condition
    ! then samples the current at each node's own elevation, as the dynamic march does, so the
    ! march starts from rest. Without it the static drag uses fluid_velocity as given.
    REAL(wp), INTENT(IN), OPTIONAL :: current_profile_z(:), current_profile_velocity(:, :)

    INTEGER, ALLOCATABLE :: elem_conn(:, :), fixed_dofs(:)
    REAL(wp), ALLOCATABLE :: l0(:), ea(:), mass_per_length(:), diameter(:)
    REAL(wp), ALLOCATABLE :: weight(:), elem_load(:, :), f_ext(:), q_static(:), v0(:)
    INTEGER :: n_elem, n_dof, e, es, istat
    LOGICAL :: dyn_to
    LOGICAL :: reverse_line, converged, use_current
    LOGICAL :: use_seabed, use_seabed_damping, use_seabed_friction, use_damping
    CHARACTER(200) :: em
    TYPE(CableCurrentLoad) :: cur
    TYPE(CD_StaticLineReport) :: report

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    reverse_line = .FALSE.
    IF (PRESENT(anchor_is_end_b)) reverse_line = anchor_is_end_b
    IF (PRESENT(seabed_kn) .NEQV. (PRESENT(seabed_z_floor) .OR. PRESENT(bathymetry))) THEN
      CALL fail(ErrStat, ErrMsg, 'line model seabed requires seabed_kn with either seabed_z_floor or bathymetry')
      RETURN
    END IF
    IF (PRESENT(seabed_z_floor) .AND. PRESENT(bathymetry)) THEN
      CALL fail(ErrStat, ErrMsg, 'line model seabed_z_floor and bathymetry are mutually exclusive')
      RETURN
    END IF
    use_seabed = PRESENT(seabed_kn)
    use_seabed_damping = PRESENT(seabed_cn)
    use_seabed_friction = PRESENT(seabed_mu)
    use_damping = PRESENT(ba)

    IF (.NOT. CD_All_Finite(anchor) .OR. .NOT. CD_All_Finite(fairlead)) THEN
      CALL fail(ErrStat, ErrMsg, 'anchor and fairlead must be finite')
      RETURN
    END IF
    IF (.NOT. CD_Is_Finite(gravity) .OR. gravity <= CD_ZERO) THEN
      CALL fail(ErrStat, ErrMsg, 'gravity must be finite and positive')
      RETURN
    END IF
    IF (.NOT. CD_Is_Finite(rho_water) .OR. rho_water <= CD_ZERO) THEN
      CALL fail(ErrStat, ErrMsg, 'rho_water must be finite and positive')
      RETURN
    END IF

    CALL CD_Build_Line_Mesh(sections, line_types, reverse_line, elem_conn, l0, ea, &
                            mass_per_length, diameter, es, em)
    IF (es /= CD_LINE_OK) THEN
      CALL fail(ErrStat, ErrMsg, 'line mesh build failed: '//TRIM(em))
      RETURN
    END IF
    n_elem = SIZE(l0)
    n_dof = 3*(n_elem + 1)
    ALLOCATE (weight(n_elem), elem_load(3, n_elem), f_ext(n_dof), q_static(n_dof), v0(n_dof), &
              STAT=istat)
    IF (istat /= 0) THEN
      CALL alloc_fail(ErrStat, ErrMsg, 'line model temporary allocation failed')
      RETURN
    END IF

    CALL CD_Submerged_Weight(mass_per_length, diameter, rho_water, gravity, weight, es, em)
    IF (es /= 0) THEN
      CALL fail(ErrStat, ErrMsg, 'submerged weight failed: '//TRIM(em))
      RETURN
    END IF
    DO e = 1, n_elem
      elem_load(:, e) = [CD_ZERO, CD_ZERO, -weight(e)]
    END DO
    CALL CD_Assemble_Distributed_Load(elem_conn, l0, elem_load, f_ext, es, em)
    IF (es /= 0) THEN
      CALL fail(ErrStat, ErrMsg, 'load assembly failed: '//TRIM(em))
      RETURN
    END IF
    IF (use_seabed) THEN
      IF (PRESENT(seabed_z_floor)) THEN
        IF (.NOT. CD_Is_Finite(seabed_z_floor)) THEN
          CALL fail(ErrStat, ErrMsg, 'seabed_z_floor must be finite')
          RETURN
        END IF
      ELSE
        IF (.NOT. CD_Bathymetry_Is_Initialized(bathymetry)) THEN
          CALL fail(ErrStat, ErrMsg, 'line model bathymetry must be initialized')
          RETURN
        END IF
      END IF
      IF (SIZE(seabed_kn) /= n_elem + 1) THEN
        CALL fail(ErrStat, ErrMsg, 'line model seabed_kn must have shape (n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(seabed_kn) .OR. ANY(seabed_kn <= CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'line model seabed_kn must be finite and positive')
        RETURN
      END IF
    END IF
    IF (use_seabed_damping) THEN
      IF (.NOT. use_seabed) THEN
        CALL fail(ErrStat, ErrMsg, 'line model seabed_cn requires seabed_z_floor and seabed_kn')
        RETURN
      END IF
      IF (SIZE(seabed_cn) /= n_elem + 1) THEN
        CALL fail(ErrStat, ErrMsg, 'line model seabed_cn must have shape (n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(seabed_cn) .OR. ANY(seabed_cn < CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'line model seabed_cn must be finite and non-negative')
        RETURN
      END IF
    END IF
    IF (use_seabed_friction) THEN
      IF (.NOT. use_seabed) THEN
        CALL fail(ErrStat, ErrMsg, 'line model seabed_mu requires seabed_z_floor and seabed_kn')
        RETURN
      END IF
      IF (.NOT. CD_Is_Finite(seabed_mu) .OR. seabed_mu < CD_ZERO) THEN
        CALL fail(ErrStat, ErrMsg, 'line model seabed_mu must be finite and non-negative')
        RETURN
      END IF
    END IF
    IF (use_damping) THEN
      IF (SIZE(ba) /= n_elem) THEN
        CALL fail(ErrStat, ErrMsg, 'line model ba must have shape (n_elem)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(ba) .OR. ANY(ba < CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'line model ba must be finite and non-negative')
        RETURN
      END IF
    END IF

    dyn_to = tension_only
    IF (PRESENT(dynamic_tension_only)) dyn_to = dynamic_tension_only

    ! Static initial condition: the shared robust line initializer (the same one the
    ! static-only deck route uses). The equilibrium obeys the committed model's constitutive
    ! law (dyn_to): a tension-only line gets its tension-only equilibrium, never a compressed
    ! strut. A steady current enters as Morison drag on the line at rest, with the node
    ! velocities and element hydrodynamics the dynamic model is built with, so the first
    ! dynamic step starts from equilibrium.
    use_current = .FALSE.
    IF (PRESENT(fluid_velocity) .AND. PRESENT(drag_rho)) THEN
      IF (SIZE(fluid_velocity, 1) == 3 .AND. SIZE(fluid_velocity, 2) == n_elem + 1) THEN
        IF (ANY(ABS(fluid_velocity) > CD_ZERO)) THEN
          IF (PRESENT(drag_diameter_elem) .AND. PRESENT(drag_cdn_elem) .AND. PRESENT(drag_cdt_elem)) THEN
            use_current = .TRUE.
            cur%elem_diameter = drag_diameter_elem
            cur%elem_cdn = drag_cdn_elem
            cur%elem_cdt = drag_cdt_elem
          ELSE IF (PRESENT(drag_diameter) .AND. PRESENT(drag_cdn) .AND. PRESENT(drag_cdt)) THEN
            use_current = .TRUE.
            cur%diameter = drag_diameter
            cur%cdn = drag_cdn
            cur%cdt = drag_cdt
          END IF
        END IF
      END IF
    END IF
    IF (use_current) THEN
      cur%rho = drag_rho
      cur%node_velocity = fluid_velocity
      IF (PRESENT(current_profile_z) .AND. PRESENT(current_profile_velocity)) THEN
        cur%profile_z = current_profile_z
        cur%profile_velocity = current_profile_velocity
      END IF
      IF (PRESENT(drag_waterline_z)) cur%node_waterline = drag_waterline_z
      ! A frictional seabed holds the line against the current from its still-water laid
      ! shape: solve the line without the current first; that equilibrium is the reference
      ! of the static friction springs, and their forces seed the stick-slip anchors below.
      IF (use_seabed .AND. use_seabed_friction) THEN
        IF (seabed_mu > CD_ZERO) THEN
          CALL CD_Static_Line_Equilibrium(anchor, fairlead, elem_conn, l0, ea, weight, f_ext, static_cfg, dyn_to, &
                                          q_static, converged, report, es, em, seabed_z_floor=seabed_z_floor, &
                                          seabed_kn=seabed_kn, bathymetry=bathymetry, load_factors=load_factors)
          IF (es /= CD_STATIC_OK .OR. .NOT. converged) THEN
            ErrStat = CD_MODEL_SOLVEFAIL
            ErrMsg = 'CableDyn_Model: still-water reference of the seabed friction failed: '//TRIM(em)// &
                     '; '//TRIM(report%message)
            RETURN
          END IF
          cur%friction_mu = seabed_mu
          IF (PRESENT(seabed_mu_axial)) cur%friction_mu_axial = seabed_mu_axial
          ALLOCATE (cur%friction_ref(2, n_elem + 1))
          DO istat = 1, n_elem + 1
            cur%friction_ref(:, istat) = q_static(3*istat - 2:3*istat - 1)
          END DO
        END IF
      END IF
      CALL CD_Static_Line_Equilibrium(anchor, fairlead, elem_conn, l0, ea, weight, f_ext, static_cfg, dyn_to, &
                                      q_static, converged, report, es, em, seabed_z_floor=seabed_z_floor, &
                                      seabed_kn=seabed_kn, bathymetry=bathymetry, current=cur, &
                                      load_factors=load_factors)
    ELSE
      CALL CD_Static_Line_Equilibrium(anchor, fairlead, elem_conn, l0, ea, weight, f_ext, static_cfg, dyn_to, &
                                      q_static, converged, report, es, em, seabed_z_floor=seabed_z_floor, &
                                      seabed_kn=seabed_kn, bathymetry=bathymetry, load_factors=load_factors)
    END IF
    IF (es /= CD_STATIC_OK) THEN
      ErrStat = CD_MODEL_SOLVEFAIL
      ErrMsg = 'CableDyn_Model: line static initialisation failed: '//TRIM(em)
      RETURN
    END IF
    IF (.NOT. converged) THEN
      ErrStat = CD_MODEL_SOLVEFAIL
      ErrMsg = 'CableDyn_Model: line static initialisation did not converge: '//TRIM(report%message)
      ! A line that has a still-water equilibrium but none in the current: name the limit.
      IF (use_current .AND. use_seabed) THEN
        BLOCK
          REAL(wp), ALLOCATABLE :: q_sw(:)
          TYPE(CD_StaticLineReport) :: report_sw
          LOGICAL :: conv_sw
          ALLOCATE (q_sw(SIZE(q_static)))
          CALL CD_Static_Line_Equilibrium(anchor, fairlead, elem_conn, l0, ea, weight, f_ext, static_cfg, dyn_to, &
                                          q_sw, conv_sw, report_sw, es, em, seabed_z_floor=seabed_z_floor, &
                                          seabed_kn=seabed_kn, bathymetry=bathymetry, load_factors=load_factors)
          IF (es == CD_STATIC_OK .AND. conv_sw) THEN
            IF (ALLOCATED(cur%friction_ref)) THEN
              ErrMsg = 'CableDyn_Model: the line has a still-water equilibrium but none in the current: '// &
                       'the seabed friction cannot hold the line against the current (it pushes the line '// &
                       'along the seabed or into the touchdown harder than the suspended span holds it); '// &
                       TRIM(report%message)
            ELSE
              ErrMsg = 'CableDyn_Model: the line has a still-water equilibrium but none in the current: '// &
                       'the current pushes the line along the frictionless seabed toward its anchor harder '// &
                       'than the suspended span holds it (declare OPTIONS frictionMu where the seabed holds '// &
                       'the run); '//TRIM(report%message)
            END IF
          END IF
        END BLOCK
      END IF
      RETURN
    END IF
    ALLOCATE (fixed_dofs(6), STAT=istat)
    IF (istat /= 0) THEN
      CALL alloc_fail(ErrStat, ErrMsg, 'line model fixed-dof allocation failed')
      RETURN
    END IF
    fixed_dofs = [1, 2, 3, n_dof - 2, n_dof - 1, n_dof]

    v0 = CD_ZERO
    IF (use_seabed .AND. use_damping) THEN
      CALL CD_Init_Model(model, q_static, v0, elem_conn, l0, ea, mass_per_length, dyn_to, &
                         f_ext, fixed_dofs, dynamic_cfg, ErrStat, ErrMsg, &
                         dist_load_per_length=elem_load, &
                         seabed_kbot=seabed_kbot, seabed_contact_diameter=seabed_contact_diameter, &
                         seabed_cn_over_kn=seabed_cn_over_kn, &
                         seabed_z_floor=seabed_z_floor, seabed_kn=seabed_kn, seabed_cn=seabed_cn, &
                         seabed_mu=seabed_mu, ba=ba, &
                         fluid_velocity=fluid_velocity, drag_waterline_z=drag_waterline_z, &
                         drag_rho=drag_rho, drag_diameter=drag_diameter, drag_cdn=drag_cdn, drag_cdt=drag_cdt, &
                         drag_diameter_elem=drag_diameter_elem, drag_cdn_elem=drag_cdn_elem, &
                         drag_cdt_elem=drag_cdt_elem, &
                         fluid_acceleration=fluid_acceleration, fk_waterline_z=fk_waterline_z, &
                         fk_rho=fk_rho, fk_diameter=fk_diameter, fk_can=fk_can, fk_cat=fk_cat, &
                         fk_diameter_elem=fk_diameter_elem, fk_can_elem=fk_can_elem, fk_cat_elem=fk_cat_elem, &
                         buoyancy_waterline_z=buoyancy_waterline_z, buoyancy_rho=buoyancy_rho, &
                         buoyancy_diameter=buoyancy_diameter, buoyancy_gravity=buoyancy_gravity, &
                         buoyancy_diameter_elem=buoyancy_diameter_elem, &
                         added_mass_waterline_z=added_mass_waterline_z, added_mass_rho=added_mass_rho, &
                         added_mass_diameter=added_mass_diameter, added_mass_can=added_mass_can, &
                         added_mass_cat=added_mass_cat, added_mass_diameter_elem=added_mass_diameter_elem, &
                         added_mass_can_elem=added_mass_can_elem, added_mass_cat_elem=added_mass_cat_elem, &
                         bathymetry=bathymetry, ve_ea_d=ve_ea_d, ve_ba=ve_ba, ve_ba_d=ve_ba_d, &
                         ve_alpha_mbl=ve_alpha_mbl, ve_vbeta=ve_vbeta, &
                         syrope_is=syrope_is, syrope_type=syrope_type, syrope_slow0=syrope_slow0, &
                         syrope_tmax0=syrope_tmax0)
    ELSE IF (use_seabed) THEN
      CALL CD_Init_Model(model, q_static, v0, elem_conn, l0, ea, mass_per_length, dyn_to, &
                         f_ext, fixed_dofs, dynamic_cfg, ErrStat, ErrMsg, &
                         dist_load_per_length=elem_load, &
                         seabed_kbot=seabed_kbot, seabed_contact_diameter=seabed_contact_diameter, &
                         seabed_cn_over_kn=seabed_cn_over_kn, &
                         seabed_z_floor=seabed_z_floor, seabed_kn=seabed_kn, seabed_cn=seabed_cn, &
                         seabed_mu=seabed_mu, &
                         fluid_velocity=fluid_velocity, drag_waterline_z=drag_waterline_z, &
                         drag_rho=drag_rho, drag_diameter=drag_diameter, drag_cdn=drag_cdn, drag_cdt=drag_cdt, &
                         drag_diameter_elem=drag_diameter_elem, drag_cdn_elem=drag_cdn_elem, &
                         drag_cdt_elem=drag_cdt_elem, &
                         fluid_acceleration=fluid_acceleration, fk_waterline_z=fk_waterline_z, &
                         fk_rho=fk_rho, fk_diameter=fk_diameter, fk_can=fk_can, fk_cat=fk_cat, &
                         fk_diameter_elem=fk_diameter_elem, fk_can_elem=fk_can_elem, fk_cat_elem=fk_cat_elem, &
                         buoyancy_waterline_z=buoyancy_waterline_z, buoyancy_rho=buoyancy_rho, &
                         buoyancy_diameter=buoyancy_diameter, buoyancy_gravity=buoyancy_gravity, &
                         buoyancy_diameter_elem=buoyancy_diameter_elem, &
                         added_mass_waterline_z=added_mass_waterline_z, added_mass_rho=added_mass_rho, &
                         added_mass_diameter=added_mass_diameter, added_mass_can=added_mass_can, &
                         added_mass_cat=added_mass_cat, added_mass_diameter_elem=added_mass_diameter_elem, &
                         added_mass_can_elem=added_mass_can_elem, added_mass_cat_elem=added_mass_cat_elem, &
                         bathymetry=bathymetry, ve_ea_d=ve_ea_d, ve_ba=ve_ba, ve_ba_d=ve_ba_d, &
                         ve_alpha_mbl=ve_alpha_mbl, ve_vbeta=ve_vbeta, &
                         syrope_is=syrope_is, syrope_type=syrope_type, syrope_slow0=syrope_slow0, &
                         syrope_tmax0=syrope_tmax0)
    ELSE IF (use_damping) THEN
      CALL CD_Init_Model(model, q_static, v0, elem_conn, l0, ea, mass_per_length, dyn_to, &
                         f_ext, fixed_dofs, dynamic_cfg, ErrStat, ErrMsg, ba=ba, &
                         ve_ea_d=ve_ea_d, ve_ba=ve_ba, ve_ba_d=ve_ba_d, ve_alpha_mbl=ve_alpha_mbl, ve_vbeta=ve_vbeta, &
                         syrope_is=syrope_is, syrope_type=syrope_type, syrope_slow0=syrope_slow0, &
                         syrope_tmax0=syrope_tmax0, &
                         dist_load_per_length=elem_load, &
                         seabed_kbot=seabed_kbot, seabed_contact_diameter=seabed_contact_diameter, &
                         seabed_cn_over_kn=seabed_cn_over_kn, &
                         fluid_velocity=fluid_velocity, drag_waterline_z=drag_waterline_z, &
                         drag_rho=drag_rho, drag_diameter=drag_diameter, drag_cdn=drag_cdn, drag_cdt=drag_cdt, &
                         drag_diameter_elem=drag_diameter_elem, drag_cdn_elem=drag_cdn_elem, &
                         drag_cdt_elem=drag_cdt_elem, &
                         fluid_acceleration=fluid_acceleration, fk_waterline_z=fk_waterline_z, &
                         fk_rho=fk_rho, fk_diameter=fk_diameter, fk_can=fk_can, fk_cat=fk_cat, &
                         fk_diameter_elem=fk_diameter_elem, fk_can_elem=fk_can_elem, fk_cat_elem=fk_cat_elem, &
                         buoyancy_waterline_z=buoyancy_waterline_z, buoyancy_rho=buoyancy_rho, &
                         buoyancy_diameter=buoyancy_diameter, buoyancy_gravity=buoyancy_gravity, &
                         buoyancy_diameter_elem=buoyancy_diameter_elem, &
                         added_mass_waterline_z=added_mass_waterline_z, added_mass_rho=added_mass_rho, &
                         added_mass_diameter=added_mass_diameter, added_mass_can=added_mass_can, &
                         added_mass_cat=added_mass_cat, added_mass_diameter_elem=added_mass_diameter_elem, &
                         added_mass_can_elem=added_mass_can_elem, added_mass_cat_elem=added_mass_cat_elem)
    ELSE
      CALL CD_Init_Model(model, q_static, v0, elem_conn, l0, ea, mass_per_length, dyn_to, &
                         f_ext, fixed_dofs, dynamic_cfg, ErrStat, ErrMsg, &
                         ve_ea_d=ve_ea_d, ve_ba=ve_ba, ve_ba_d=ve_ba_d, ve_alpha_mbl=ve_alpha_mbl, ve_vbeta=ve_vbeta, &
                         syrope_is=syrope_is, syrope_type=syrope_type, syrope_slow0=syrope_slow0, &
                         syrope_tmax0=syrope_tmax0, &
                         dist_load_per_length=elem_load, &
                         seabed_kbot=seabed_kbot, seabed_contact_diameter=seabed_contact_diameter, &
                         seabed_cn_over_kn=seabed_cn_over_kn, &
                         fluid_velocity=fluid_velocity, drag_waterline_z=drag_waterline_z, &
                         drag_rho=drag_rho, drag_diameter=drag_diameter, drag_cdn=drag_cdn, drag_cdt=drag_cdt, &
                         drag_diameter_elem=drag_diameter_elem, drag_cdn_elem=drag_cdn_elem, &
                         drag_cdt_elem=drag_cdt_elem, &
                         fluid_acceleration=fluid_acceleration, fk_waterline_z=fk_waterline_z, &
                         fk_rho=fk_rho, fk_diameter=fk_diameter, fk_can=fk_can, fk_cat=fk_cat, &
                         fk_diameter_elem=fk_diameter_elem, fk_can_elem=fk_can_elem, fk_cat_elem=fk_cat_elem, &
                         buoyancy_waterline_z=buoyancy_waterline_z, buoyancy_rho=buoyancy_rho, &
                         buoyancy_diameter=buoyancy_diameter, buoyancy_gravity=buoyancy_gravity, &
                         buoyancy_diameter_elem=buoyancy_diameter_elem, &
                         added_mass_waterline_z=added_mass_waterline_z, added_mass_rho=added_mass_rho, &
                         added_mass_diameter=added_mass_diameter, added_mass_can=added_mass_can, &
                         added_mass_cat=added_mass_cat, added_mass_diameter_elem=added_mass_diameter_elem, &
                         added_mass_can_elem=added_mass_can_elem, added_mass_cat_elem=added_mass_cat_elem)
    END IF
    IF (ErrStat /= CD_MODEL_OK) RETURN
    IF (PRESENT(seabed_mu_axial)) THEN
      IF (CD_Model_Has_Friction(model) .AND. seabed_mu_axial > CD_ZERO) THEN
        CALL CD_Set_Model_Friction_Axial(model, seabed_mu_axial, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END IF
    END IF
    ! Continue the static friction springs into the dynamics: the stick-slip anchors carry
    ! the static spring forces, so the march starts in the static equilibrium.
    IF (ALLOCATED(cur%friction_ref) .AND. CD_Model_Has_Friction(model)) THEN
      BLOCK
        REAL(wp), ALLOCATABLE :: anchors(:, :)
        ALLOCATE (anchors(2, n_elem + 1))
        IF (PRESENT(bathymetry)) THEN
          CALL CD_Static_Friction_Anchors(q_static, cur, seabed_kn, anchors, es, em, bathymetry=bathymetry)
        ELSE
          CALL CD_Static_Friction_Anchors(q_static, cur, seabed_kn, anchors, es, em, seabed_z_floor=seabed_z_floor)
        END IF
        IF (es /= CD_STATIC_OK) THEN
          es = CD_MODEL_SOLVEFAIL
        ELSE
          CALL CD_Set_Model_Friction_Anchors(model, anchors, es, em)
          IF (es == CD_MODEL_OK) CALL CD_Recompute_Model_Acceleration(model, es, em)
        END IF
        IF (es /= CD_MODEL_OK) THEN
          ErrStat = CD_MODEL_SOLVEFAIL
          ErrMsg = 'CableDyn_Model: seabed friction anchors failed: '//TRIM(em)
          RETURN
        END IF
      END BLOCK
    END IF
  END SUBROUTINE CD_Init_Line_Model

  SUBROUTINE CD_Step_Model(model, dt, converged, stalled, n_iter, ErrStat, ErrMsg, &
                           prescribed_q, prescribed_v, prescribed_a)
    !! Advance the persistent model by one generalized-alpha step. Optional
    !! prescribed_q/v/a arrays are full model-sized vectors; only entries in
    !! model%fixed_dofs are consumed by CableDyn_Dynamic. On successful or stalled
    !! dynamics returns, the model owns the returned state. On solver failure, the
    !! underlying step contract preserves the previous state.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: dt
    LOGICAL, INTENT(OUT) :: converged, stalled
    INTEGER, INTENT(OUT) :: n_iter, ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: prescribed_q(:), prescribed_v(:), prescribed_a(:)

    INTEGER :: es
    CHARACTER(160) :: em
    LOGICAL :: has_load

    converged = .FALSE.
    stalled = .FALSE.
    n_iter = 0
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (PRESENT(prescribed_q) .OR. PRESENT(prescribed_v) .OR. PRESENT(prescribed_a)) THEN
      IF (.NOT. (PRESENT(prescribed_q) .AND. PRESENT(prescribed_v) .AND. PRESENT(prescribed_a))) THEN
        CALL fail(ErrStat, ErrMsg, 'prescribed motion needs all of prescribed_q, prescribed_v, prescribed_a')
        RETURN
      END IF
      IF (SIZE(prescribed_q) /= model%n_dof .OR. SIZE(prescribed_v) /= model%n_dof &
          .OR. SIZE(prescribed_a) /= model%n_dof) THEN
        CALL fail(ErrStat, ErrMsg, 'prescribed_q/v/a must each have shape (3 n_nodes)')
        RETURN
      END IF
    END IF

    has_load = model%has_seabed .OR. model%has_seabed_damping .OR. model%has_seabed_friction .OR. &
               model%has_damping .OR. model%has_viscoelastic .OR. model%has_syrope .OR. model%has_morison_drag .OR. &
               model%has_froude_krylov .OR. model%has_buoyancy_recovery
    ! arm the viscoelastic contributors with the step size: in-step load
    ! evaluations use the backward-Euler-eliminated dl_1(q, v); it is reset to
    ! the instantaneous limit (0) on EVERY exit below so committed-state queries
    ! (coupled loads, tension) never see a stale step size
    model%ve_dt = dt
    IF (PRESENT(prescribed_q) .AND. has_load .AND. model%has_added_mass) THEN
      CALL CD_Cable_Gen_Alpha_Step(model%q, model%v, model%a, model%elem_conn, model%l0, model%ea_dyn, &
                                   model%rho_a, model%tension_only, model%f_ext, model%fixed_dofs, &
                                   dt, model%cfg, model%q_step, model%v_step, model%a_step, converged, stalled, &
                                   n_iter, &
                                   es, em, load_force_proc=model_force, load_band_proc=model_load_band, &
                                   added_mass_band_proc=model_added_mass_band, prescribed_q=prescribed_q, &
                                   prescribed_v=prescribed_v, prescribed_a=prescribed_a, &
                                   workspace=model%dynamic_workspace)
    ELSE IF (PRESENT(prescribed_q) .AND. has_load) THEN
      CALL CD_Cable_Gen_Alpha_Step(model%q, model%v, model%a, model%elem_conn, model%l0, model%ea_dyn, &
                                   model%rho_a, model%tension_only, model%f_ext, model%fixed_dofs, &
                                   dt, model%cfg, model%q_step, model%v_step, model%a_step, converged, stalled, &
                                   n_iter, &
                                   es, em, load_force_proc=model_force, load_band_proc=model_load_band, &
                                   prescribed_q=prescribed_q, prescribed_v=prescribed_v, &
                                   prescribed_a=prescribed_a, workspace=model%dynamic_workspace)
    ELSE IF (PRESENT(prescribed_q) .AND. model%has_added_mass) THEN
      CALL CD_Cable_Gen_Alpha_Step(model%q, model%v, model%a, model%elem_conn, model%l0, model%ea_dyn, &
                                   model%rho_a, model%tension_only, model%f_ext, model%fixed_dofs, &
                                   dt, model%cfg, model%q_step, model%v_step, model%a_step, converged, stalled, &
                                   n_iter, &
                                   es, em, added_mass_band_proc=model_added_mass_band, prescribed_q=prescribed_q, &
                                   prescribed_v=prescribed_v, prescribed_a=prescribed_a, &
                                   workspace=model%dynamic_workspace)
    ELSE IF (PRESENT(prescribed_q)) THEN
      CALL CD_Cable_Gen_Alpha_Step(model%q, model%v, model%a, model%elem_conn, model%l0, model%ea_dyn, &
                                   model%rho_a, model%tension_only, model%f_ext, model%fixed_dofs, &
                                   dt, model%cfg, model%q_step, model%v_step, model%a_step, converged, stalled, &
                                   n_iter, &
                                   es, em, prescribed_q=prescribed_q, prescribed_v=prescribed_v, &
                                   prescribed_a=prescribed_a, workspace=model%dynamic_workspace)
    ELSE IF (has_load .AND. model%has_added_mass) THEN
      CALL CD_Cable_Gen_Alpha_Step(model%q, model%v, model%a, model%elem_conn, model%l0, model%ea_dyn, &
                                   model%rho_a, model%tension_only, model%f_ext, model%fixed_dofs, &
                                   dt, model%cfg, model%q_step, model%v_step, model%a_step, converged, stalled, &
                                   n_iter, es, em, &
                                   load_force_proc=model_force, load_band_proc=model_load_band, &
                                   added_mass_band_proc=model_added_mass_band, workspace=model%dynamic_workspace)
    ELSE IF (has_load) THEN
      CALL CD_Cable_Gen_Alpha_Step(model%q, model%v, model%a, model%elem_conn, model%l0, model%ea_dyn, &
                                   model%rho_a, model%tension_only, model%f_ext, model%fixed_dofs, &
                                   dt, model%cfg, model%q_step, model%v_step, model%a_step, converged, stalled, &
                                   n_iter, es, em, &
                                   load_force_proc=model_force, load_band_proc=model_load_band, &
                                   workspace=model%dynamic_workspace)
    ELSE IF (model%has_added_mass) THEN
      CALL CD_Cable_Gen_Alpha_Step(model%q, model%v, model%a, model%elem_conn, model%l0, model%ea_dyn, &
                                   model%rho_a, model%tension_only, model%f_ext, model%fixed_dofs, &
                                   dt, model%cfg, model%q_step, model%v_step, model%a_step, converged, stalled, &
                                   n_iter, es, em, &
                                   added_mass_band_proc=model_added_mass_band, workspace=model%dynamic_workspace)
    ELSE
      CALL CD_Cable_Gen_Alpha_Step(model%q, model%v, model%a, model%elem_conn, model%l0, model%ea_dyn, &
                                   model%rho_a, model%tension_only, model%f_ext, model%fixed_dofs, &
                                   dt, model%cfg, model%q_step, model%v_step, model%a_step, converged, stalled, &
                                   n_iter, es, em, workspace=model%dynamic_workspace)
    END IF
    model%ve_dt = CD_ZERO
    IF (es /= CD_DYN_OK) THEN
      ErrStat = model_status_from_dyn(es)
      ErrMsg = 'CableDyn_Model: step failed: '//TRIM(em)
      converged = .FALSE.
      RETURN
    END IF
    ! ==== stateful constitutive commit: validate EVERY model (viscoelastic AND
    ! Syrope) BEFORE mutating ANY internal state, so a failure in one model cannot
    ! leave another model's state advanced past a step that never commits q/v/a.
    ! Each pass-1 is exhaustive enough that its pass-2 cannot fail mid-mutation, so
    ! the four passes together are one atomic commit. Each element reads/writes only
    ! its own state, and viscoelastic vs Syrope elements are disjoint segments, so
    ! the pass ordering does not change any success-path result (bit-identical). ====
    IF (model%has_viscoelastic) THEN
      ! Pass 1 (validate): every viscoelastic element length + lagged coefficients
      ! (deterministic from the COMMITTED dl_1), the only failure modes of the advance.
      BLOCK
        INTEGER :: e2, na, nb, es2
        REAL(wp) :: ea_d_eff, ea_1_eff
        CHARACTER(200) :: em2
        DO e2 = 1, model%n_elem
          IF (model%ve_mode(e2) == 0) CYCLE
          na = model%elem_conn(1, e2)
          nb = model%elem_conn(2, e2)
          IF (NORM2(model%q_step(3*nb - 2:3*nb) - model%q_step(3*na - 2:3*na)) <= &
              CD_VISCO_TINY_LEN) THEN
            CALL fail(ErrStat, ErrMsg, 'viscoelastic state advance found a collapsed element in the '// &
                      'converged step state')
            converged = .FALSE.
            RETURN
          END IF
          CALL model_viscoelastic_coeffs(model, e2, ea_d_eff, ea_1_eff, es2, em2)
          IF (es2 /= CD_MODEL_OK) THEN
            ErrStat = es2
            ErrMsg = em2
            converged = .FALSE.
            RETURN
          END IF
        END DO
      END BLOCK
    END IF
    IF (model%has_syrope) THEN
      ! Pass 1 (validate): every Syrope element (collapsed length, working-curve
      ! regeneration at the current AND any ratcheted maximum).
      BLOCK
        INTEGER :: e2, idx_scr(6)
        REAL(wp) :: q2b(6), v2b(6), slow_next, tmean
        DO e2 = 1, model%n_elem
          IF (.NOT. model%syrope_is(e2)) CYCLE
          CALL extract_element_state(model, e2, model%q_step, model%v_step, idx_scr, q2b, v2b)
          CALL model_syrope_state_next(model, e2, q2b, v2b, dt, slow_next, tmean, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) THEN
            converged = .FALSE.
            RETURN
          END IF
        END DO
      END BLOCK
    END IF
    ! Friction anchors of the converged state, mapped before the state is committed so a
    ! failure leaves the step-start state (and its anchors) intact; the passes below are
    ! validated and cannot fail.
    IF (model%has_seabed_friction) THEN
      CALL mapped_friction_anchors(model, model%q_step, model%v_step, .FALSE., ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) THEN
        converged = .FALSE.
        RETURN
      END IF
    END IF
    IF (model%has_viscoelastic) THEN
      ! Pass 2 (mutate): commit dl_1 from the SAME lagged coefficients the in-step
      ! loads used (mode 3 reads the still-committed dl_1). Validated above, so no
      ! failure can occur here mid-mutation.
      BLOCK
        INTEGER :: e2, na, nb, es2
        REAL(wp) :: dl1_next, ea_d_eff, ea_1_eff
        CHARACTER(200) :: em2
        DO e2 = 1, model%n_elem
          IF (model%ve_mode(e2) == 0) CYCLE
          CALL model_viscoelastic_coeffs(model, e2, ea_d_eff, ea_1_eff, es2, em2)
          IF (es2 /= CD_MODEL_OK) THEN
            ErrStat = es2
            ErrMsg = em2
            converged = .FALSE.
            RETURN
          END IF
          na = model%elem_conn(1, e2)
          nb = model%elem_conn(2, e2)
          CALL CD_Viscoelastic_State_Advance(model%q_step(3*na - 2:3*na), model%q_step(3*nb - 2:3*nb), &
                                             model%v_step(3*na - 2:3*na), model%v_step(3*nb - 2:3*nb), &
                                             model%l0(e2), ea_d_eff, ea_1_eff, &
                                             model%ve_ba(e2), model%ve_ba_d(e2), model%ve_dl_1(e2), dt, &
                                             dl1_next, es2, em2)
          IF (es2 /= CD_VISCO_OK) THEN
            CALL fail(ErrStat, ErrMsg, TRIM(em2))
            converged = .FALSE.
            RETURN
          END IF
          model%ve_dl_1(e2) = dl1_next
        END DO
      END BLOCK
    END IF
    IF (model%has_syrope) THEN
      ! Pass 2 (mutate): advance the committed slow strain and ratchet the running
      ! maximum. The element re-eval reads the still-committed slow/tmax (mutated
      ! only after), and each element writes only its own state. Validated above.
      BLOCK
        INTEGER :: e2, idx_scr(6)
        REAL(wp) :: q2b(6), v2b(6), tmean, slow_next
        DO e2 = 1, model%n_elem
          IF (.NOT. model%syrope_is(e2)) CYCLE
          CALL extract_element_state(model, e2, model%q_step, model%v_step, idx_scr, q2b, v2b)
          CALL model_syrope_state_next(model, e2, q2b, v2b, dt, slow_next, tmean, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) THEN
            converged = .FALSE.
            RETURN
          END IF
          model%syrope_slow(e2) = slow_next
          IF (tmean > model%syrope_tmax(e2)) model%syrope_tmax(e2) = tmean
        END DO
      END BLOCK
    END IF
    model%q = model%q_step
    model%v = model%v_step
    model%a = model%a_step
    ! The mapping repeats the queries validated above, so it cannot fail here.
    IF (model%has_seabed_friction) &
      CALL mapped_friction_anchors(model, model%q, model%v, .TRUE., ErrStat, ErrMsg)

  CONTAINS

    SUBROUTINE model_force(q, v, force, ErrStat, ErrMsg)
      REAL(wp), INTENT(IN) :: q(:), v(:)
      REAL(wp), INTENT(OUT) :: force(:)
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg
      CALL model_contributor_force(model, q, v, force, ErrStat, ErrMsg)
      ErrStat = as_dyn_callback_status(ErrStat)
    END SUBROUTINE model_force

    SUBROUTINE model_load_band(q, v, free, kl, ku, force, jac_q_band, jac_v_band, ErrStat, ErrMsg)
      REAL(wp), INTENT(IN) :: q(:), v(:)
      INTEGER, INTENT(IN) :: free(:), kl, ku
      REAL(wp), INTENT(OUT) :: force(:), jac_q_band(:, :), jac_v_band(:, :)
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg
      CALL model_contributor_load_banded(model, q, v, free, kl, ku, force, jac_q_band, jac_v_band, &
                                         ErrStat, ErrMsg)
    END SUBROUTINE model_load_band

    SUBROUTINE model_added_mass_band(q, accel, free, kl, ku, M_add_a, M_add_band, dMa_a_dq_band, ErrStat, ErrMsg, &
                                     need_tangent)
      REAL(wp), INTENT(IN) :: q(:), accel(:)
      INTEGER, INTENT(IN) :: free(:), kl, ku
      REAL(wp), INTENT(OUT) :: M_add_a(:), M_add_band(:, :), dMa_a_dq_band(:, :)
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg
      LOGICAL, INTENT(IN), OPTIONAL :: need_tangent
      CALL model_added_mass_banded(model, q, accel, free, kl, ku, M_add_a, M_add_band, dMa_a_dq_band, &
                                   ErrStat, ErrMsg, need_tangent=need_tangent)
      ErrStat = as_dyn_callback_status(ErrStat)
    END SUBROUTINE model_added_mass_band
  END SUBROUTINE CD_Step_Model

  SUBROUTINE CD_Step_Model_Recovering(model, dt, converged, stalled, n_iter, ErrStat, ErrMsg, &
                                      prescribed_q, prescribed_v, prescribed_a, max_substeps, substeps_used)
    !! Advance one EI=0 model over a nominal interval, recovering a clean
    !! nonlinear stall by internal temporal subdivision.  The nominal interval
    !! and its prescribed endpoint kinematics remain unchanged.  A C2 quintic
    !! joins the committed and target boundary position, velocity and
    !! acceleration, so the subdivision follows one kinematically consistent
    !! trajectory rather than a sequence of independently interpolated fields.
    !!
    !! Every recovery rung starts from the same complete committed state,
    !! including viscoelastic and Syrope history variables.  If no rung
    !! completes, that state is restored and the routine fails closed.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: dt
    LOGICAL, INTENT(OUT) :: converged, stalled
    INTEGER, INTENT(OUT) :: n_iter, ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: prescribed_q(:), prescribed_v(:), prescribed_a(:)
    INTEGER, INTENT(IN), OPTIONAL :: max_substeps
    INTEGER, INTENT(OUT), OPTIONAL :: substeps_used

    INTEGER :: recovery_cap, nsub, k, es, iter_sub, max_iter_seen
    LOGICAL :: has_pres, conv_sub, stall_sub, rung_ok
    REAL(wp) :: subdt, fraction
    REAL(wp), ALLOCATABLE :: pq(:), pv(:), pa(:)
    CHARACTER(300) :: em, first_em

    converged = .FALSE.
    stalled = .FALSE.
    n_iter = 0
    IF (PRESENT(substeps_used)) substeps_used = 0
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: recovery step requires an initialized model'
      RETURN
    END IF
    has_pres = PRESENT(prescribed_q) .AND. PRESENT(prescribed_v) .AND. PRESENT(prescribed_a)
    IF ((PRESENT(prescribed_q) .OR. PRESENT(prescribed_v) .OR. PRESENT(prescribed_a)) .AND. .NOT. has_pres) THEN
      CALL fail(ErrStat, ErrMsg, 'recovery step needs all of prescribed_q, prescribed_v and prescribed_a')
      RETURN
    END IF
    IF (has_pres) THEN
      IF (SIZE(prescribed_q) /= model%n_dof .OR. SIZE(prescribed_v) /= model%n_dof .OR. &
          SIZE(prescribed_a) /= model%n_dof) THEN
        CALL fail(ErrStat, ErrMsg, 'recovery prescribed_q/v/a must each have shape (3 n_nodes)')
        RETURN
      END IF
      IF (.NOT. (CD_All_Finite(prescribed_q) .AND. CD_All_Finite(prescribed_v) .AND. &
                 CD_All_Finite(prescribed_a))) THEN
        CALL fail(ErrStat, ErrMsg, 'recovery prescribed_q/v/a must be finite')
        RETURN
      END IF
    END IF
    recovery_cap = MODEL_DEFAULT_MAX_SUBSTEPS
    IF (PRESENT(max_substeps)) recovery_cap = max_substeps
    IF (recovery_cap < 4 .OR. recovery_cap > 65536) THEN
      CALL fail(ErrStat, ErrMsg, 'recovery max_substeps must be in [4,65536]')
      RETURN
    END IF
    IF (.NOT. (ALLOCATED(model%recovery_q0) .AND. ALLOCATED(model%recovery_v0) .AND. &
               ALLOCATED(model%recovery_a0))) THEN
      CALL fail(ErrStat, ErrMsg, 'recovery snapshot workspace is not initialized')
      RETURN
    END IF

    CALL save_recovery_state(model)
    IF (has_pres) THEN
      CALL CD_Step_Model(model, dt, converged, stalled, n_iter, ErrStat, ErrMsg, &
                         prescribed_q=prescribed_q, prescribed_v=prescribed_v, prescribed_a=prescribed_a)
    ELSE
      CALL CD_Step_Model(model, dt, converged, stalled, n_iter, ErrStat, ErrMsg)
    END IF
    max_iter_seen = n_iter
    IF (ErrStat /= CD_MODEL_OK) THEN
      CALL restore_recovery_state(model)
      RETURN
    END IF
    IF (converged) THEN
      IF (PRESENT(substeps_used)) substeps_used = 1
      RETURN
    END IF

    first_em = ErrMsg
    CALL restore_recovery_state(model)
    IF (has_pres) THEN
      ALLOCATE (pq(model%n_dof), pv(model%n_dof), pa(model%n_dof), STAT=es)
      IF (es /= 0) THEN
        ! The committed state was restored above; fail closed without a retry.
        converged = .FALSE.
        stalled = .TRUE.
        ErrStat = CD_MODEL_ALLOCFAIL
        ErrMsg = 'CableDyn_Model: recovery boundary workspace allocation failed'
        RETURN
      END IF
    END IF
    nsub = 4
    DO
      CALL restore_recovery_state(model)
      subdt = dt/REAL(nsub, wp)
      rung_ok = .TRUE.
      es = CD_MODEL_OK
      em = ''
      DO k = 1, nsub
        IF (has_pres) THEN
          fraction = REAL(k, wp)/REAL(nsub, wp)
          CALL quintic_model_boundary_state(model%recovery_q0, model%recovery_v0, model%recovery_a0, &
                                            prescribed_q, prescribed_v, prescribed_a, dt, fraction, pq, pv, pa)
          IF (k == nsub) THEN
            pq = prescribed_q
            pv = prescribed_v
            pa = prescribed_a
          END IF
          CALL CD_Step_Model(model, subdt, conv_sub, stall_sub, iter_sub, es, em, &
                             prescribed_q=pq, prescribed_v=pv, prescribed_a=pa)
        ELSE
          CALL CD_Step_Model(model, subdt, conv_sub, stall_sub, iter_sub, es, em)
        END IF
        max_iter_seen = MAX(max_iter_seen, iter_sub)
        IF (es /= CD_MODEL_OK .OR. .NOT. conv_sub) THEN
          rung_ok = .FALSE.
          EXIT
        END IF
      END DO
      IF (rung_ok) THEN
        converged = .TRUE.
        stalled = .FALSE.
        n_iter = max_iter_seen
        IF (PRESENT(substeps_used)) substeps_used = nsub
        ErrStat = CD_MODEL_OK
        ErrMsg = ''
        RETURN
      END IF
      IF (es /= CD_MODEL_OK) EXIT
      IF (nsub == recovery_cap) EXIT
      IF (nsub > recovery_cap/4) THEN
        nsub = recovery_cap
      ELSE
        nsub = 4*nsub
      END IF
    END DO

    CALL restore_recovery_state(model)
    converged = .FALSE.
    stalled = .TRUE.
    n_iter = max_iter_seen
    IF (es /= CD_MODEL_OK) THEN
      ErrStat = es
    ELSE
      ErrStat = CD_MODEL_SOLVEFAIL
    END IF
    IF (LEN_TRIM(em) == 0) em = first_em
    IF (LEN_TRIM(em) == 0) em = 'nonlinear solve did not converge'
    BLOCK
      CHARACTER(1024) :: wbuf
      INTEGER :: wios
      wbuf = ''
      WRITE (wbuf, '(A,ES12.4,A,I0,A)', IOSTAT=wios) 'CableDyn_Model: recovery failed over dt=', dt, &
        ' s at the subdivision cap ', recovery_cap, ': '//TRIM(em)
      ErrMsg = wbuf
    END BLOCK
  END SUBROUTINE CD_Step_Model_Recovering

  SUBROUTINE save_recovery_state(model)
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    model%recovery_q0 = model%q
    model%recovery_v0 = model%v
    model%recovery_a0 = model%a
    IF (model%has_viscoelastic) model%recovery_ve_dl_1 = model%ve_dl_1
    IF (model%has_syrope) THEN
      model%recovery_syrope_slow = model%syrope_slow
      model%recovery_syrope_tmax = model%syrope_tmax
    END IF
    IF (model%has_seabed_friction) model%recovery_fr_anchor = model%fr_anchor
  END SUBROUTINE save_recovery_state

  SUBROUTINE restore_recovery_state(model)
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    model%q = model%recovery_q0
    model%v = model%recovery_v0
    model%a = model%recovery_a0
    IF (model%has_viscoelastic) model%ve_dl_1 = model%recovery_ve_dl_1
    IF (model%has_syrope) THEN
      model%syrope_slow = model%recovery_syrope_slow
      model%syrope_tmax = model%recovery_syrope_tmax
    END IF
    IF (model%has_seabed_friction) model%fr_anchor = model%recovery_fr_anchor
    model%ve_dt = CD_ZERO
  END SUBROUTINE restore_recovery_state

  SUBROUTINE mapped_friction_anchors(model, q, v, commit, ErrStat, ErrMsg)
    !! Return mapping of the stick-slip friction springs at the state (q, v): a spring
    !! stretched beyond its capacity (mu times the normal reaction) moves its anchor to the
    !! capacity distance behind the node, so the committed force is the one the converged
    !! residual used; a node off the seabed (capacity zero) carries its anchor along. With
    !! commit false only the queries run (the step validates them before it commits any
    !! state); with commit true the anchors are updated.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:)
    LOGICAL, INTENT(IN) :: commit
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: i, ix
    REAL(wp) :: normal, cap, d(2), dlen, k, mu_d, q_d, gq_d(2)
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    DO i = 1, SIZE(model%fr_anchor, 2)
      ix = 3*i - 2
      normal = seabed_normal_reaction(model, q(ix), q(ix + 1), q(ix + 2), v(ix:ix + 2), i, ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
      IF (.NOT. commit) CYCLE
      IF (model%seabed_friction_aniso) THEN
        d = q(ix:ix + 1) - model%fr_anchor(:, i)
        CALL CD_Seabed_Friction_Mu_Dir(d, model%seabed_mu_axial, model%seabed_mu, friction_axis(q, i), &
                                       mu_d, q_d, gq_d)
        cap = mu_d*MAX(normal, CD_ZERO)
      ELSE
        cap = model%seabed_mu*MAX(normal, CD_ZERO)
      END IF
      k = model%seabed_kn(i)
      d = q(ix:ix + 1) - model%fr_anchor(:, i)
      dlen = SQRT(d(1)*d(1) + d(2)*d(2))
      IF (k*dlen > cap) model%fr_anchor(:, i) = q(ix:ix + 1) - (cap/(k*dlen))*d
    END DO
  END SUBROUTINE mapped_friction_anchors

  SUBROUTINE quintic_model_boundary_state(q0, v0, a0, q1, v1, a1, duration, u, q, v, a)
    REAL(wp), INTENT(IN) :: q0(:), v0(:), a0(:), q1(:), v1(:), a1(:), duration, u
    REAL(wp), INTENT(OUT) :: q(:), v(:), a(:)
    REAL(wp) :: u2, u3, u4, u5, t2

    u2 = u*u
    u3 = u2*u
    u4 = u3*u
    u5 = u4*u
    t2 = duration*duration
    q = (CD_ONE - 10.0_wp*u3 + 15.0_wp*u4 - 6.0_wp*u5)*q0 + &
        (u - 6.0_wp*u3 + 8.0_wp*u4 - 3.0_wp*u5)*duration*v0 + &
        0.5_wp*(u2 - 3.0_wp*u3 + 3.0_wp*u4 - u5)*t2*a0 + &
        (10.0_wp*u3 - 15.0_wp*u4 + 6.0_wp*u5)*q1 + &
        (-4.0_wp*u3 + 7.0_wp*u4 - 3.0_wp*u5)*duration*v1 + &
        0.5_wp*(u3 - 2.0_wp*u4 + u5)*t2*a1
    v = (-30.0_wp*u2 + 60.0_wp*u3 - 30.0_wp*u4)*q0/duration + &
        (CD_ONE - 18.0_wp*u2 + 32.0_wp*u3 - 15.0_wp*u4)*v0 + &
        (u - 4.5_wp*u2 + 6.0_wp*u3 - 2.5_wp*u4)*duration*a0 + &
        (30.0_wp*u2 - 60.0_wp*u3 + 30.0_wp*u4)*q1/duration + &
        (-12.0_wp*u2 + 28.0_wp*u3 - 15.0_wp*u4)*v1 + &
        (1.5_wp*u2 - 4.0_wp*u3 + 2.5_wp*u4)*duration*a1
    a = (-60.0_wp*u + 180.0_wp*u2 - 120.0_wp*u3)*q0/t2 + &
        (-36.0_wp*u + 96.0_wp*u2 - 60.0_wp*u3)*v0/duration + &
        (CD_ONE - 9.0_wp*u + 18.0_wp*u2 - 10.0_wp*u3)*a0 + &
        (60.0_wp*u - 180.0_wp*u2 + 120.0_wp*u3)*q1/t2 + &
        (-24.0_wp*u + 84.0_wp*u2 - 60.0_wp*u3)*v1/duration + &
        (3.0_wp*u - 12.0_wp*u2 + 10.0_wp*u3)*a1
  END SUBROUTINE quintic_model_boundary_state

  SUBROUTINE ensure_model_query_workspace(model, need_mass, need_added, need_tangent, need_reduced, need_rhs_cols, &
                                          ErrStat, ErrMsg)
    !! Grow reusable model-owned scratch for output/load derivative queries.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    LOGICAL, INTENT(IN) :: need_mass, need_added, need_tangent
    INTEGER, INTENT(IN) :: need_reduced, need_rhs_cols
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: istat, n

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (need_reduced < 0 .OR. need_rhs_cols < 0) THEN
      CALL fail(ErrStat, ErrMsg, 'query workspace reduced size must be non-negative')
      RETURN
    END IF
    n = model%n_dof

    IF (vector_needs_capacity(model%dynamic_workspace%R, n)) THEN
      IF (ALLOCATED(model%dynamic_workspace%R)) DEALLOCATE (model%dynamic_workspace%R)
      ALLOCATE (model%dynamic_workspace%R(n), STAT=istat)
      IF (istat /= 0) THEN
        CALL alloc_fail(ErrStat, ErrMsg, 'query residual workspace allocation failed')
        RETURN
      END IF
    END IF
    IF (vector_needs_capacity(model%dynamic_workspace%fint_eval, n)) THEN
      IF (ALLOCATED(model%dynamic_workspace%fint_eval)) DEALLOCATE (model%dynamic_workspace%fint_eval)
      ALLOCATE (model%dynamic_workspace%fint_eval(n), STAT=istat)
      IF (istat /= 0) THEN
        CALL alloc_fail(ErrStat, ErrMsg, 'query internal-force workspace allocation failed')
        RETURN
      END IF
    END IF
    IF (vector_needs_capacity(model%dynamic_workspace%load_eval, n)) THEN
      IF (ALLOCATED(model%dynamic_workspace%load_eval)) DEALLOCATE (model%dynamic_workspace%load_eval)
      ALLOCATE (model%dynamic_workspace%load_eval(n), STAT=istat)
      IF (istat /= 0) THEN
        CALL alloc_fail(ErrStat, ErrMsg, 'query load workspace allocation failed')
        RETURN
      END IF
    END IF

    IF (need_mass) THEN
      IF (matrix_needs_capacity(model%dynamic_workspace%M, n, n)) THEN
        IF (ALLOCATED(model%dynamic_workspace%M)) DEALLOCATE (model%dynamic_workspace%M)
        ALLOCATE (model%dynamic_workspace%M(n, n), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query mass workspace allocation failed')
          RETURN
        END IF
      END IF
    END IF

    IF (need_added) THEN
      IF (matrix_needs_capacity(model%dynamic_workspace%M_add_eval, n, n)) THEN
        IF (ALLOCATED(model%dynamic_workspace%M_add_eval)) DEALLOCATE (model%dynamic_workspace%M_add_eval)
        ALLOCATE (model%dynamic_workspace%M_add_eval(n, n), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query added-mass workspace allocation failed')
          RETURN
        END IF
      END IF
      IF (matrix_needs_capacity(model%dynamic_workspace%dMa_a_dq_eval, n, n)) THEN
        IF (ALLOCATED(model%dynamic_workspace%dMa_a_dq_eval)) DEALLOCATE (model%dynamic_workspace%dMa_a_dq_eval)
        ALLOCATE (model%dynamic_workspace%dMa_a_dq_eval(n, n), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query added-mass tangent workspace allocation failed')
          RETURN
        END IF
      END IF
    END IF

    IF (need_tangent) THEN
      IF (matrix_needs_capacity(model%dynamic_workspace%Kt_eval, n, n)) THEN
        IF (ALLOCATED(model%dynamic_workspace%Kt_eval)) DEALLOCATE (model%dynamic_workspace%Kt_eval)
        ALLOCATE (model%dynamic_workspace%Kt_eval(n, n), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query tangent workspace allocation failed')
          RETURN
        END IF
      END IF
      IF (matrix_needs_capacity(model%dynamic_workspace%jac_q_eval, n, n)) THEN
        IF (ALLOCATED(model%dynamic_workspace%jac_q_eval)) DEALLOCATE (model%dynamic_workspace%jac_q_eval)
        ALLOCATE (model%dynamic_workspace%jac_q_eval(n, n), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query dload/dq workspace allocation failed')
          RETURN
        END IF
      END IF
      IF (matrix_needs_capacity(model%dynamic_workspace%jac_v_eval, n, n)) THEN
        IF (ALLOCATED(model%dynamic_workspace%jac_v_eval)) DEALLOCATE (model%dynamic_workspace%jac_v_eval)
        ALLOCATE (model%dynamic_workspace%jac_v_eval(n, n), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query dload/dv workspace allocation failed')
          RETURN
        END IF
      END IF
      IF (vector_needs_capacity(model%dynamic_workspace%tension_eval, model%n_elem)) THEN
        IF (ALLOCATED(model%dynamic_workspace%tension_eval)) DEALLOCATE (model%dynamic_workspace%tension_eval)
        ALLOCATE (model%dynamic_workspace%tension_eval(model%n_elem), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query tension workspace allocation failed')
          RETURN
        END IF
      END IF
    END IF

    IF (need_reduced > 0) THEN
      IF (matrix_needs_capacity(model%dynamic_workspace%eff_free, need_reduced, need_reduced)) THEN
        IF (ALLOCATED(model%dynamic_workspace%eff_free)) DEALLOCATE (model%dynamic_workspace%eff_free)
        ALLOCATE (model%dynamic_workspace%eff_free(need_reduced, need_reduced), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query reduced-matrix workspace allocation failed')
          RETURN
        END IF
      END IF
    END IF

    IF (need_reduced > 0 .AND. need_rhs_cols > 0) THEN
      IF (matrix_needs_capacity(model%dynamic_workspace%eff, need_reduced, need_rhs_cols)) THEN
        IF (ALLOCATED(model%dynamic_workspace%eff)) DEALLOCATE (model%dynamic_workspace%eff)
        ALLOCATE (model%dynamic_workspace%eff(need_reduced, need_rhs_cols), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query rhs workspace allocation failed')
          RETURN
        END IF
      END IF
      IF (matrix_needs_capacity(model%dynamic_workspace%M_eff_eval, need_reduced, need_rhs_cols)) THEN
        IF (ALLOCATED(model%dynamic_workspace%M_eff_eval)) DEALLOCATE (model%dynamic_workspace%M_eff_eval)
        ALLOCATE (model%dynamic_workspace%M_eff_eval(need_reduced, need_rhs_cols), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'query secondary-rhs workspace allocation failed')
          RETURN
        END IF
      END IF
    END IF

  CONTAINS

    LOGICAL FUNCTION vector_needs_capacity(x, n_min) RESULT(needs)
      REAL(wp), ALLOCATABLE, INTENT(IN) :: x(:)
      INTEGER, INTENT(IN) :: n_min

      needs = .NOT. ALLOCATED(x)
      IF (.NOT. needs) needs = SIZE(x) < n_min
    END FUNCTION vector_needs_capacity

    LOGICAL FUNCTION matrix_needs_capacity(x, n1_min, n2_min) RESULT(needs)
      REAL(wp), ALLOCATABLE, INTENT(IN) :: x(:, :)
      INTEGER, INTENT(IN) :: n1_min, n2_min

      needs = .NOT. ALLOCATED(x)
      IF (.NOT. needs) needs = SIZE(x, 1) < n1_min .OR. SIZE(x, 2) < n2_min
    END FUNCTION matrix_needs_capacity
  END SUBROUTINE ensure_model_query_workspace

  SUBROUTINE fill_nodes_workspace(model, ErrStat)
    !! model%dynamic_workspace%nodes_eval = RESHAPE(model%q, [3, n_nodes]), in storage that
    !! is allocated once and reused.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    INTEGER :: n_nodes, k, istat

    ErrStat = CD_MODEL_OK
    n_nodes = model%n_dof/3
    IF (ALLOCATED(model%dynamic_workspace%nodes_eval)) THEN
      IF (SIZE(model%dynamic_workspace%nodes_eval, 2) /= n_nodes) DEALLOCATE (model%dynamic_workspace%nodes_eval)
    END IF
    IF (.NOT. ALLOCATED(model%dynamic_workspace%nodes_eval)) THEN
      ALLOCATE (model%dynamic_workspace%nodes_eval(3, n_nodes), STAT=istat)
      IF (istat /= 0) THEN
        ErrStat = CD_MODEL_ALLOCFAIL
        RETURN
      END IF
    END IF
    DO k = 1, n_nodes
      model%dynamic_workspace%nodes_eval(:, k) = model%q(3*k - 2:3*k)
    END DO
  END SUBROUTINE fill_nodes_workspace

  SUBROUTINE partition_model_free_workspace(model, n_free, ErrStat, ErrMsg)
    !! Fill model%free_map_work with free DOFs using model-owned marker scratch.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: n_free
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: i, dof, istat
    LOGICAL :: free_map_ready, need_marker

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    n_free = 0
    free_map_ready = ALLOCATED(model%free_map_work)
    IF (free_map_ready) free_map_ready = SIZE(model%free_map_work) >= model%n_dof
    IF (.NOT. free_map_ready) THEN
      CALL fail(ErrStat, ErrMsg, 'free-DOF partition requires initialized free-map workspace')
      RETURN
    END IF
    need_marker = .NOT. ALLOCATED(model%dynamic_workspace%dof_marker)
    IF (.NOT. need_marker) need_marker = SIZE(model%dynamic_workspace%dof_marker) < model%n_dof
    IF (need_marker) THEN
      IF (ALLOCATED(model%dynamic_workspace%dof_marker)) DEALLOCATE (model%dynamic_workspace%dof_marker)
      ALLOCATE (model%dynamic_workspace%dof_marker(model%n_dof), STAT=istat)
      IF (istat /= 0) THEN
        CALL alloc_fail(ErrStat, ErrMsg, 'free-DOF marker workspace allocation failed')
        RETURN
      END IF
    END IF
    model%dynamic_workspace%dof_marker(1:model%n_dof) = 0
    DO i = 1, SIZE(model%fixed_dofs)
      dof = model%fixed_dofs(i)
      IF (dof < 1 .OR. dof > model%n_dof) THEN
        CALL fail(ErrStat, ErrMsg, 'fixed DOF is out of range')
        RETURN
      END IF
      IF (model%dynamic_workspace%dof_marker(dof) /= 0) THEN
        CALL fail(ErrStat, ErrMsg, 'fixed DOFs contain duplicates')
        RETURN
      END IF
      model%dynamic_workspace%dof_marker(dof) = 1
    END DO
    DO dof = 1, model%n_dof
      IF (model%dynamic_workspace%dof_marker(dof) == 0) THEN
        n_free = n_free + 1
        model%free_map_work(n_free) = dof
      END IF
    END DO
  END SUBROUTINE partition_model_free_workspace

  SUBROUTINE CD_Get_Model_State(model, q, v, a, ErrStat, ErrMsg)
    !! Copy the current model state into caller-owned arrays.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: q(:), v(:), a(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    q = CD_ZERO
    v = CD_ZERO
    a = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (SIZE(q) /= model%n_dof .OR. SIZE(v) /= model%n_dof .OR. SIZE(a) /= model%n_dof) THEN
      CALL fail(ErrStat, ErrMsg, 'q, v, and a outputs must each have shape (3 n_nodes)')
      RETURN
    END IF
    q = model%q
    v = model%v
    a = model%a
  END SUBROUTINE CD_Get_Model_State

  PURE LOGICAL FUNCTION CD_Model_Has_Viscoelastic(model) RESULT(has)
    !! True if the model carries per-element series-Kelvin viscoelastic state.
    TYPE(CD_ModelType), INTENT(IN) :: model
    has = model%has_viscoelastic
  END FUNCTION CD_Model_Has_Viscoelastic

  SUBROUTINE CD_Get_Model_VE_Dl1(model, dl1, ErrStat, ErrMsg)
    !! Copy the per-element viscoelastic internal strain dl_1 (size n_elem) for the
    !! coupled checkpoint mirror. A non-viscoelastic model returns zeros (the mirror
    !! never queries one, but keep the accessor total).
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: dl1(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    dl1 = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (SIZE(dl1) /= model%n_elem) THEN
      CALL fail(ErrStat, ErrMsg, 'dl1 output must have shape (n_elem)')
      RETURN
    END IF
    IF (model%has_viscoelastic) dl1 = model%ve_dl_1
  END SUBROUTINE CD_Get_Model_VE_Dl1

  SUBROUTINE CD_Set_Model_VE_Dl1(model, dl1, ErrStat, ErrMsg)
    !! Overwrite the per-element viscoelastic internal strain dl_1 from a checkpoint
    !! mirror (the cross-process restart reload). Fails closed on a non-viscoelastic
    !! model, a shape mismatch, or a non-finite mirror.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: dl1(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (.NOT. model%has_viscoelastic) THEN
      CALL fail(ErrStat, ErrMsg, 'model carries no viscoelastic state to set')
      RETURN
    END IF
    IF (SIZE(dl1) /= model%n_elem) THEN
      CALL fail(ErrStat, ErrMsg, 'dl1 input must have shape (n_elem)')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(dl1)) THEN
      CALL fail(ErrStat, ErrMsg, 'dl1 mirror must be finite')
      RETURN
    END IF
    model%ve_dl_1 = dl1
    ! The rollback and recovery twins follow, so a later rollback cannot bring back the
    ! state this replaces (as CD_Set_Model_Friction_Anchors does).
    model%ve_dl_1_rollback = dl1
    model%recovery_ve_dl_1 = dl1
  END SUBROUTINE CD_Set_Model_VE_Dl1

  PURE LOGICAL FUNCTION CD_Model_Has_Syrope(model) RESULT(has)
    !! True if the model carries per-element Syrope (working-curve) internal state.
    TYPE(CD_ModelType), INTENT(IN) :: model
    has = model%has_syrope
  END FUNCTION CD_Model_Has_Syrope

  SUBROUTINE CD_Get_Model_Syrope_State(model, slow, tmax, ErrStat, ErrMsg)
    !! Copy the per-element Syrope committed states -- the slow-spring static strain
    !! slow and the running-maximum tension tmax (each size n_elem) -- for the coupled
    !! checkpoint mirror. A non-Syrope model returns zeros (total accessor).
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: slow(:), tmax(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    slow = CD_ZERO
    tmax = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (SIZE(slow) /= model%n_elem .OR. SIZE(tmax) /= model%n_elem) THEN
      CALL fail(ErrStat, ErrMsg, 'slow and tmax outputs must each have shape (n_elem)')
      RETURN
    END IF
    IF (model%has_syrope) THEN
      slow = model%syrope_slow
      tmax = model%syrope_tmax
    END IF
  END SUBROUTINE CD_Get_Model_Syrope_State

  SUBROUTINE CD_Set_Model_Syrope_State(model, slow, tmax, ErrStat, ErrMsg)
    !! Overwrite the per-element Syrope committed states (slow-spring static strain and
    !! running-maximum tension) from a checkpoint mirror (the cross-process restart
    !! reload). Fails closed on a non-Syrope model, a shape mismatch, or a non-finite
    !! mirror.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: slow(:), tmax(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: e, es
    REAL(wp) :: ws(CD_SYROPE_NWC), wt(CD_SYROPE_NWC), wsl(CD_SYROPE_NWC)
    CHARACTER(160) :: em
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (.NOT. model%has_syrope) THEN
      CALL fail(ErrStat, ErrMsg, 'model carries no Syrope state to set')
      RETURN
    END IF
    IF (SIZE(slow) /= model%n_elem .OR. SIZE(tmax) /= model%n_elem) THEN
      CALL fail(ErrStat, ErrMsg, 'slow and tmax inputs must each have shape (n_elem)')
      RETURN
    END IF
    IF (.NOT. (CD_All_Finite(slow) .AND. CD_All_Finite(tmax))) THEN
      CALL fail(ErrStat, ErrMsg, 'Syrope state mirror must be finite')
      RETURN
    END IF
    ! CD_Syrope_Working_Curve rejects a nonpositive running maximum, so fail closed at the
    ! reload boundary (a corrupted or edited checkpoint) rather than storing an invalid
    ! tmax that only trips the next load evaluation. Only Syrope elements carry a live tmax.
    IF (ANY(model%syrope_is .AND. .NOT. (tmax > CD_ZERO))) THEN
      CALL fail(ErrStat, ErrMsg, 'Syrope running-maximum tension (tmax) must be positive for every '// &
                'Syrope element')
      RETURN
    END IF
    IF (ANY(model%syrope_is .AND. slow < CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'Syrope slow-spring strain must be non-negative for every Syrope element')
      RETURN
    END IF
    DO e = 1, model%n_elem
      IF (.NOT. model%syrope_is(e)) CYCLE
      CALL CD_Syrope_Working_Curve(model%syrope_type(e), tmax(e), ws, wt, wsl, es, em)
      IF (es /= CD_SYROPE_OK) THEN
        CALL fail(ErrStat, ErrMsg, 'Syrope state mirror has invalid tmax: '//TRIM(em))
        RETURN
      END IF
      CALL CD_Syrope_Check_Range(model%syrope_type(e), slow(e), es, em)
      IF (es /= CD_SYROPE_OK) THEN
        CALL fail(ErrStat, ErrMsg, 'Syrope state mirror has invalid slow strain: '//TRIM(em))
        RETURN
      END IF
    END DO
    model%syrope_slow = slow
    model%syrope_tmax = tmax
    ! The rollback and recovery twins follow (see CD_Set_Model_VE_Dl1).
    model%syrope_slow_rollback = slow
    model%syrope_tmax_rollback = tmax
    model%recovery_syrope_slow = slow
    model%recovery_syrope_tmax = tmax
  END SUBROUTINE CD_Set_Model_Syrope_State

  SUBROUTINE CD_Set_Model_Friction_Axial(model, mu_axial, ErrStat, ErrMsg)
    !! Make the stick-slip seabed friction of a model anisotropic: the model's seabed_mu
    !! becomes the lateral coefficient and mu_axial the coefficient along the line (the
    !! horizontal projection of the node's chord tangent, re-evaluated at every residual
    !! call; the Jacobian omits its dependence on the positions). mu_axial equal to
    !! seabed_mu keeps the isotropic law bit-for-bit.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: mu_axial
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. (CD_Is_Finite(mu_axial) .AND. mu_axial > CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'the axial friction coefficient must be finite and positive')
      RETURN
    END IF
    IF (.NOT. model%has_seabed_friction) THEN
      CALL fail(ErrStat, ErrMsg, 'the model carries no seabed friction')
      RETURN
    END IF
    model%seabed_mu_axial = mu_axial
    model%seabed_friction_aniso = ABS(mu_axial - model%seabed_mu) > CD_ZERO
  END SUBROUTINE CD_Set_Model_Friction_Axial

  PURE FUNCTION friction_axis(q, node) RESULT(axis)
    !! Horizontal line axis at a node for anisotropic friction: the chord through its
    !! neighbours (one-sided at the ends), projected on the horizontal. The neighbours
    !! are nodes node-1 and node+1: the line models number their nodes along the line
    !! (CD_Build_Line_Mesh), which the anisotropic friction assumes.
    REAL(wp), INTENT(IN) :: q(:)
    INTEGER, INTENT(IN) :: node
    REAL(wp) :: axis(2)
    INTEGER :: nn, ia, ib
    nn = SIZE(q)/3
    ia = MAX(1, node - 1)
    ib = MIN(nn, node + 1)
    axis = q(3*ib - 2:3*ib - 1) - q(3*ia - 2:3*ia - 1)
  END FUNCTION friction_axis

  PURE LOGICAL FUNCTION CD_Model_Has_Friction(model) RESULT(has)
    !! True if the model carries stick-slip seabed friction (anchor state).
    TYPE(CD_ModelType), INTENT(IN) :: model
    has = model%initialized .AND. model%has_seabed_friction .AND. ALLOCATED(model%fr_anchor)
  END FUNCTION CD_Model_Has_Friction

  SUBROUTINE CD_Get_Model_Friction_Anchors(model, anchors, ErrStat, ErrMsg)
    !! Copy the committed stick-slip friction anchors, (2, n_nodes).
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: anchors(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    anchors = CD_ZERO
    IF (.NOT. CD_Model_Has_Friction(model)) THEN
      CALL fail(ErrStat, ErrMsg, 'the model carries no seabed friction state')
      RETURN
    END IF
    IF (SIZE(anchors, 1) /= 2 .OR. SIZE(anchors, 2) /= SIZE(model%fr_anchor, 2)) THEN
      CALL fail(ErrStat, ErrMsg, 'friction anchors must have shape (2, n_nodes)')
      RETURN
    END IF
    anchors = model%fr_anchor
  END SUBROUTINE CD_Get_Model_Friction_Anchors

  SUBROUTINE CD_Set_Model_Friction_Anchors(model, anchors, ErrStat, ErrMsg)
    !! Overwrite the committed stick-slip friction anchors, (2, n_nodes): from a static
    !! solve with friction (then call CD_Recompute_Model_Acceleration) or a checkpoint
    !! mirror. The rollback and recovery twins follow, so a later rollback cannot bring
    !! back the anchors this replaces.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: anchors(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. CD_Model_Has_Friction(model)) THEN
      CALL fail(ErrStat, ErrMsg, 'the model carries no seabed friction state')
      RETURN
    END IF
    IF (SIZE(anchors, 1) /= 2 .OR. SIZE(anchors, 2) /= SIZE(model%fr_anchor, 2)) THEN
      CALL fail(ErrStat, ErrMsg, 'friction anchors must have shape (2, n_nodes)')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(anchors)) THEN
      CALL fail(ErrStat, ErrMsg, 'friction anchors must be finite')
      RETURN
    END IF
    model%fr_anchor = anchors
    model%fr_anchor_rollback = anchors
    model%recovery_fr_anchor = anchors
  END SUBROUTINE CD_Set_Model_Friction_Anchors

  SUBROUTINE CD_Calc_Model_CoupledLoads(model, coupled_loads, ErrStat, ErrMsg)
    !! Return loads exerted BY the CableDyn model ON the externally coupled DOFs.
    !! The internal dynamic residual at prescribed DOFs is
    !!   R_c = M*a + f_int - f_ext,
    !! which is the support force needed on the cable. The equal-and-opposite force
    !! that the cable applies to the coupled object is therefore -R_c. This is the
    !! force-vector half of the OpenFAST/C kinematics-in/loads-out boundary.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(OUT) :: coupled_loads(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: es, istat, k
    CHARACTER(160) :: em

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    coupled_loads = CD_ZERO
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (SIZE(coupled_loads) /= SIZE(model%fixed_dofs)) THEN
      CALL fail(ErrStat, ErrMsg, 'coupled_loads output must have shape (n_coupled_dof)')
      RETURN
    END IF

    CALL ensure_model_query_workspace(model, need_mass=.FALSE., need_added=.FALSE., &
                                      need_tangent=.FALSE., need_reduced=0, need_rhs_cols=0, &
                                      ErrStat=ErrStat, ErrMsg=ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    IF (model%has_added_mass) THEN
      IF (.NOT. ALLOCATED(model%dynamic_workspace%M_add_a_eval)) THEN
        ALLOCATE (model%dynamic_workspace%M_add_a_eval(model%n_dof), STAT=istat)
      ELSE IF (SIZE(model%dynamic_workspace%M_add_a_eval) < model%n_dof) THEN
        DEALLOCATE (model%dynamic_workspace%M_add_a_eval)
        ALLOCATE (model%dynamic_workspace%M_add_a_eval(model%n_dof), STAT=istat)
      ELSE
        istat = 0
      END IF
      IF (istat /= 0) THEN
        CALL alloc_fail(ErrStat, ErrMsg, 'coupled-load added-mass force workspace allocation failed')
        RETURN
      END IF
    END IF
    ASSOCIATE (fint => model%dynamic_workspace%fint_eval(1:model%n_dof), &
               residual => model%dynamic_workspace%R(1:model%n_dof), &
               load => model%dynamic_workspace%load_eval(1:model%n_dof))
      CALL fill_nodes_workspace(model, es)
      IF (es /= CD_MODEL_OK) THEN
        CALL alloc_fail(ErrStat, ErrMsg, 'nodal workspace allocation failed')
        RETURN
      END IF
      CALL CD_Assemble_Cable_Internal_Force(model%dynamic_workspace%nodes_eval, model%elem_conn, &
                                            model%l0, model%ea_dyn, model%tension_only, fint, es, em, &
                                            topology_validated=.TRUE.)
      IF (es /= 0) THEN
        ErrStat = CD_MODEL_SOLVEFAIL
        ErrMsg = 'CableDyn_Model: internal-force assembly failed during coupled-load calculation: '//TRIM(em)
        RETURN
      END IF
      CALL model_structural_mass_product(model, model%a, residual)
      residual = residual + fint - model%f_ext
      IF (model%has_added_mass) THEN
        CALL model_added_mass_force(model, model%q, model%a, &
                                    model%dynamic_workspace%M_add_a_eval(1:model%n_dof), es, em)
        IF (es /= 0) THEN
          ErrStat = CD_MODEL_SOLVEFAIL
          ErrMsg = 'CableDyn_Model: added-mass assembly failed during coupled-load calculation: '//TRIM(em)
          RETURN
        END IF
        residual = residual + model%dynamic_workspace%M_add_a_eval(1:model%n_dof)
      END IF
      CALL model_contributor_force(model, model%q, model%v, load, es, em)
      IF (es /= 0) THEN
        ErrStat = CD_MODEL_SOLVEFAIL
        ErrMsg = 'CableDyn_Model: load assembly failed during coupled-load calculation: '//TRIM(em)
        RETURN
      END IF
      residual = residual - load
      DO k = 1, SIZE(model%fixed_dofs)
        coupled_loads(k) = -residual(model%fixed_dofs(k))
      END DO
    END ASSOCIATE
  END SUBROUTINE CD_Calc_Model_CoupledLoads

  SUBROUTINE CD_Calc_Model_CoupledKinematicDerivatives(model, dload_dq, dload_dv, ErrStat, ErrMsg)
    !! Analytic reduced coupled-load derivatives with respect to prescribed
    !! coupled positions and velocities. Free DOFs are statically condensed with
    !! the current effective mass, matching the acceleration recompute boundary.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(OUT) :: dload_dq(:, :), dload_dv(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: es, nc, n_free, i, j, row, col
    CHARACTER(160) :: em

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    dload_dq = CD_ZERO
    dload_dv = CD_ZERO
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    nc = SIZE(model%fixed_dofs)
    IF (SIZE(dload_dq, 1) /= nc .OR. SIZE(dload_dq, 2) /= nc .OR. &
        SIZE(dload_dv, 1) /= nc .OR. SIZE(dload_dv, 2) /= nc) THEN
      CALL fail(ErrStat, ErrMsg, 'dload_dq/dload_dv outputs must have shape (n_coupled_dof, n_coupled_dof)')
      RETURN
    END IF
    IF (nc == 0) RETURN

    CALL partition_model_free_workspace(model, n_free, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL ensure_model_query_workspace(model, need_mass=.TRUE., need_added=model%has_added_mass, &
                                      need_tangent=.TRUE., need_reduced=n_free, need_rhs_cols=nc, &
                                      ErrStat=ErrStat, ErrMsg=ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    ASSOCIATE (Kt => model%dynamic_workspace%Kt_eval(1:model%n_dof, 1:model%n_dof), &
               fint => model%dynamic_workspace%fint_eval(1:model%n_dof), &
               tension => model%dynamic_workspace%tension_eval(1:model%n_elem), &
               load => model%dynamic_workspace%load_eval(1:model%n_dof), &
               jac_q => model%dynamic_workspace%jac_q_eval(1:model%n_dof, 1:model%n_dof), &
               jac_v => model%dynamic_workspace%jac_v_eval(1:model%n_dof, 1:model%n_dof), &
               M => model%dynamic_workspace%M(1:model%n_dof, 1:model%n_dof))
      CALL fill_nodes_workspace(model, es)
      IF (es /= CD_MODEL_OK) THEN
        CALL alloc_fail(ErrStat, ErrMsg, 'nodal workspace allocation failed')
        RETURN
      END IF
      CALL CD_Assemble_Cable_Tangent_Force(model%dynamic_workspace%nodes_eval, model%elem_conn, model%l0, &
                                           model%ea_dyn, model%tension_only, Kt, fint, tension, es, em, &
                                           topology_validated=.TRUE.)
      IF (es /= 0) THEN
        ErrStat = CD_MODEL_SOLVEFAIL
        ErrMsg = 'CableDyn_Model: tangent assembly failed during kinematic derivative: '//TRIM(em)
        RETURN
      END IF
      CALL model_contributor_load(model, model%q, model%v, load, jac_q, jac_v, es, em)
      IF (es /= 0) THEN
        ! Preserve an allocation failure from the load-contributor workspace; only genuine
        ! solver failures collapse to CD_MODEL_SOLVEFAIL.
        IF (es == CD_MODEL_ALLOCFAIL) THEN
          ErrStat = CD_MODEL_ALLOCFAIL
        ELSE
          ErrStat = CD_MODEL_SOLVEFAIL
        END IF
        ErrMsg = 'CableDyn_Model: load Jacobian assembly failed during kinematic derivative: '//TRIM(em)
        RETURN
      END IF
      CALL CD_Assemble_Cable_Mass(model%elem_conn, model%l0, model%rho_a, M, es, em, topology_validated=.TRUE.)
      IF (es /= 0) THEN
        ErrStat = CD_MODEL_SOLVEFAIL
        ErrMsg = 'CableDyn_Model: mass assembly failed during kinematic derivative: '//TRIM(em)
        RETURN
      END IF
      IF (model%has_added_mass) THEN
        ASSOCIATE (M_add => model%dynamic_workspace%M_add_eval(1:model%n_dof, 1:model%n_dof), &
                   dMa_a_dq => model%dynamic_workspace%dMa_a_dq_eval(1:model%n_dof, 1:model%n_dof))
          CALL model_added_mass_full(model, model%q, model%a, M_add, dMa_a_dq, es, em)
          IF (es /= 0) THEN
            ErrStat = CD_MODEL_SOLVEFAIL
            ErrMsg = 'CableDyn_Model: added-mass Jacobian assembly failed during kinematic derivative: '//TRIM(em)
            RETURN
          END IF
          M = M + M_add
          Kt = Kt + dMa_a_dq
        END ASSOCIATE
      END IF
      Kt = Kt + jac_q

      DO j = 1, nc
        col = model%fixed_dofs(j)
        DO i = 1, nc
          row = model%fixed_dofs(i)
          dload_dq(i, j) = -Kt(row, col)
          dload_dv(i, j) = -jac_v(row, col)
        END DO
      END DO
      IF (n_free == 0) RETURN

      ASSOCIATE (Mff => model%dynamic_workspace%eff_free(1:n_free, 1:n_free), &
                 xq => model%dynamic_workspace%eff(1:n_free, 1:nc), &
                 xv => model%dynamic_workspace%M_eff_eval(1:n_free, 1:nc))
        DO j = 1, n_free
          col = model%free_map_work(j)
          DO i = 1, n_free
            row = model%free_map_work(i)
            Mff(i, j) = M(row, col)
          END DO
        END DO
        DO j = 1, nc
          col = model%fixed_dofs(j)
          DO i = 1, n_free
            row = model%free_map_work(i)
            xq(i, j) = Kt(row, col)
            xv(i, j) = jac_v(row, col)
          END DO
        END DO
        CALL CD_Solve_Dense_As_Banded_Multiple(Mff, xq, es, em)
        IF (es /= 0) THEN
          ErrStat = CD_MODEL_SOLVEFAIL
          ErrMsg = 'CableDyn_Model: position derivative condensation solve failed: '//TRIM(em)
          RETURN
        END IF
        CALL CD_Solve_Dense_As_Banded_Multiple(Mff, xv, es, em)
        IF (es /= 0) THEN
          ErrStat = CD_MODEL_SOLVEFAIL
          ErrMsg = 'CableDyn_Model: velocity derivative condensation solve failed: '//TRIM(em)
          RETURN
        END IF
        DO j = 1, nc
          DO i = 1, nc
            row = model%fixed_dofs(i)
            DO col = 1, n_free
              dload_dq(i, j) = dload_dq(i, j) + M(row, model%free_map_work(col))*xq(col, j)
              dload_dv(i, j) = dload_dv(i, j) + M(row, model%free_map_work(col))*xv(col, j)
            END DO
          END DO
        END DO
      END ASSOCIATE
    END ASSOCIATE
  END SUBROUTINE CD_Calc_Model_CoupledKinematicDerivatives

  SUBROUTINE CD_Update_Model_CoupledAccelDerivative(model, ErrStat, ErrMsg)
    !! Bring the coupled-load derivative with respect to prescribed coupled accelerations
    !! up to date in model%dynamic_workspace%ad_result (n_coupled_dof square; see
    !! CD_Calc_Model_CoupledAccelDerivative, which copies it out). For free DOFs eliminated
    !! by the consistent-mass solve, d(load_c)/d(a_c) = -(M_cc - M_cf M_ff^{-1} M_fc), with M
    !! including added mass.
    !! M is block-tridiagonal in node order, so it is assembled directly in band
    !! storage (structural and added mass summed per entry as in the dense form),
    !! M_ff is factored once and solved for all coupled columns together, and the
    !! result is cached until an input of the mass changes. Work and storage are
    !! O(n_dof); the workspace is reused, so a step allocates nothing.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: es, nc, n, n_elem, n_free, i, j, k, e, a, b, ia, ib, bwf, kl, ldab, r, c, g, fi, fj
    INTEGER :: idx(6), conn1(2, 1)
    REAL(wp) :: coeff, s, q2(6), l02(1), wl2(2), m2(6, 6)
    LOGICAL :: hit, ok
    CHARACTER(160) :: em

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    nc = SIZE(model%fixed_dofs)
    IF (nc == 0) RETURN
    n = model%n_dof
    n_elem = model%n_elem

    ASSOCIATE (w => model%dynamic_workspace)
      ! Cache hit: the mass depends only on the positions (through the added mass), the
      ! unstretched lengths, the mass per length and the fixed-DOF set; the hydrodynamic
      ! coefficients are set at initialisation, which clears this workspace, and a moved
      ! added-mass waterline clears it in CD_Update_Model_Hydro_Fields. Equality is
      ! tested as <= and >= so that a NaN never matches.
      hit = w%ad_valid
      IF (hit) hit = SIZE(w%ad_fixed) == nc .AND. SIZE(w%ad_l0) == n_elem .AND. SIZE(w%ad_q) == n
      IF (hit) hit = ALL(w%ad_fixed == model%fixed_dofs)
      IF (hit) hit = ALL(w%ad_l0 <= model%l0(1:n_elem) .AND. w%ad_l0 >= model%l0(1:n_elem))
      IF (hit) hit = ALL(w%ad_rho <= model%rho_a(1:n_elem) .AND. w%ad_rho >= model%rho_a(1:n_elem))
      IF (hit .AND. model%has_added_mass) hit = ALL(w%ad_q <= model%q(1:n) .AND. w%ad_q >= model%q(1:n))
      IF (hit) RETURN
      w%ad_valid = .FALSE.

      DO e = 1, n_elem
        IF (.NOT. (model%l0(e) > CD_ZERO .AND. model%rho_a(e) > CD_ZERO) .OR. &
            .NOT. CD_Is_Finite(model%l0(e)) .OR. .NOT. CD_Is_Finite(model%rho_a(e))) THEN
          ErrStat = CD_MODEL_SOLVEFAIL
          ErrMsg = 'CableDyn_Model: mass assembly failed during acceleration derivative: '// &
                   'l0 and rho_a must be positive and finite'
          RETURN
        END IF
      END DO
      ! half-bandwidth of the mass in the full node ordering
      bwf = 2
      DO e = 1, n_elem
        bwf = MAX(bwf, 3*ABS(model%elem_conn(1, e) - model%elem_conn(2, e)) + 2)
      END DO

      ok = ensure_i1(w%ad_finv, n)
      IF (ok) ok = ensure_i1(w%ad_fixed, nc)
      IF (ok) ok = ensure_r1(w%ad_q, n)
      IF (ok) ok = ensure_r1(w%ad_l0, n_elem)
      IF (ok) ok = ensure_r1(w%ad_rho, n_elem)
      IF (ok) ok = ensure_r2(w%ad_result, nc, nc)
      IF (ok) ok = ensure_r2(w%ad_full, 2*bwf + 1, n)
      IF (ok .AND. model%has_added_mass) ok = ensure_r2(w%ad_add, 2*bwf + 1, n)
      IF (.NOT. ok) THEN
        CALL alloc_fail(ErrStat, ErrMsg, 'acceleration-derivative workspace allocation failed')
        RETURN
      END IF

      w%ad_result = CD_ZERO
      ! free-DOF numbering in ascending global order (0 marks a fixed DOF)
      w%ad_finv = 1
      DO i = 1, nc
        g = model%fixed_dofs(i)
        IF (g < 1 .OR. g > n) THEN
          CALL fail(ErrStat, ErrMsg, 'fixed DOF is out of range')
          RETURN
        END IF
        IF (w%ad_finv(g) == 0) THEN
          CALL fail(ErrStat, ErrMsg, 'fixed DOFs contain duplicates')
          RETURN
        END IF
        w%ad_finv(g) = 0
      END DO
      n_free = 0
      DO g = 1, n
        IF (w%ad_finv(g) /= 0) THEN
          n_free = n_free + 1
          w%ad_finv(g) = n_free
        END IF
      END DO

      ! structural consistent mass, element by element: m/6 [[2I, I], [I, 2I]]
      w%ad_full = CD_ZERO
      DO e = 1, n_elem
        a = model%elem_conn(1, e)
        b = model%elem_conn(2, e)
        coeff = model%rho_a(e)*model%l0(e)/6.0_wp
        ia = 3*a - 2
        ib = 3*b - 2
        DO k = 0, 2
          w%ad_full(bwf + 1, ia + k) = w%ad_full(bwf + 1, ia + k) + 2.0_wp*coeff
          w%ad_full(bwf + 1 + ia - ib, ib + k) = w%ad_full(bwf + 1 + ia - ib, ib + k) + coeff
          w%ad_full(bwf + 1 + ib - ia, ia + k) = w%ad_full(bwf + 1 + ib - ia, ia + k) + coeff
          w%ad_full(bwf + 1, ib + k) = w%ad_full(bwf + 1, ib + k) + 2.0_wp*coeff
        END DO
      END DO
      ! added mass, assembled on its own and then summed entry by entry
      IF (model%has_added_mass) THEN
        w%ad_add = CD_ZERO
        conn1(:, 1) = [1, 2]
        DO e = 1, n_elem
          CALL extract_element_positions(model, e, model%q, idx, q2)
          CALL extract_nodal_scalar(model%added_mass_waterline_z, model, e, wl2)
          l02(1) = model%l0(e)
          CALL CD_Cable_Added_Mass_Matrix(q2, conn1, l02, wl2, model%added_mass_rho, &
                                          model%added_mass_diameter_elem(e), model%added_mass_can_elem(e), &
                                          model%added_mass_cat_elem(e), m2, es, em)
          IF (es /= 0) THEN
            ErrStat = CD_MODEL_SOLVEFAIL
            ErrMsg = 'CableDyn_Model: added-mass assembly failed during acceleration derivative: '//TRIM(em)
            RETURN
          END IF
          DO j = 1, 6
            DO i = 1, 6
              r = bwf + 1 + idx(i) - idx(j)
              w%ad_add(r, idx(j)) = w%ad_add(r, idx(j)) + m2(i, j)
            END DO
          END DO
        END DO
        w%ad_full = w%ad_full + w%ad_add
      END IF

      DO j = 1, nc
        c = model%fixed_dofs(j)
        DO i = 1, nc
          r = model%fixed_dofs(i)
          IF (ABS(r - c) <= bwf) w%ad_result(i, j) = -w%ad_full(bwf + 1 + r - c, c)
        END DO
      END DO

      IF (n_free > 0) THEN
        IF (.NOT. CD_All_Finite(w%ad_full)) THEN
          ErrStat = CD_MODEL_SOLVEFAIL
          ErrMsg = 'CableDyn_Model: acceleration derivative mass solve failed: the mass must be finite'
          RETURN
        END IF
        ! half-bandwidth of the free/free block in the free numbering
        kl = 0
        DO c = 1, n
          fj = w%ad_finv(c)
          IF (fj == 0) CYCLE
          DO r = MAX(1, c - bwf), MIN(n, c + bwf)
            fi = w%ad_finv(r)
            IF (fi /= 0) kl = MAX(kl, ABS(fi - fj))
          END DO
        END DO
        ldab = 3*kl + 1
        ok = ensure_r2(w%ad_ab, ldab, n_free)
        IF (ok) ok = ensure_r2(w%ad_x, n_free, nc)
        IF (ok) ok = ensure_i1(w%ad_ipiv, n_free)
        IF (.NOT. ok) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'acceleration-derivative workspace allocation failed')
          RETURN
        END IF
        w%ad_kl = kl
        ! M_ff in LAPACK general-band storage (kl = ku)
        w%ad_ab = CD_ZERO
        DO c = 1, n
          fj = w%ad_finv(c)
          IF (fj == 0) CYCLE
          DO r = MAX(1, c - bwf), MIN(n, c + bwf)
            fi = w%ad_finv(r)
            IF (fi /= 0) w%ad_ab(2*kl + 1 + fi - fj, fj) = w%ad_full(bwf + 1 + r - c, c)
          END DO
        END DO
        ! right-hand sides M_fc
        w%ad_x = CD_ZERO
        DO j = 1, nc
          c = model%fixed_dofs(j)
          DO r = MAX(1, c - bwf), MIN(n, c + bwf)
            fi = w%ad_finv(r)
            IF (fi /= 0) w%ad_x(fi, j) = w%ad_full(bwf + 1 + r - c, c)
          END DO
        END DO
        CALL CD_Factor_Banded(w%ad_ab, kl, kl, w%ad_ipiv, es, em, matrix_validated=.TRUE.)
        IF (es == 0) CALL CD_Solve_Factored_Banded_Multiple(w%ad_ab, kl, kl, w%ad_ipiv, w%ad_x, es, em)
        IF (es /= 0) THEN
          ErrStat = CD_MODEL_SOLVEFAIL
          ErrMsg = 'CableDyn_Model: acceleration derivative mass solve failed: '//TRIM(em)
          RETURN
        END IF
        ! + M_cf M_ff^{-1} M_fc, summed over the free columns in ascending order
        DO j = 1, nc
          DO i = 1, nc
            r = model%fixed_dofs(i)
            s = w%ad_result(i, j)
            DO g = MAX(1, r - bwf), MIN(n, r + bwf)
              fi = w%ad_finv(g)
              IF (fi /= 0) s = s + w%ad_full(bwf + 1 + r - g, g)*w%ad_x(fi, j)
            END DO
            w%ad_result(i, j) = s
          END DO
        END DO
      END IF

      w%ad_fixed = model%fixed_dofs
      w%ad_l0 = model%l0(1:n_elem)
      w%ad_rho = model%rho_a(1:n_elem)
      w%ad_q = model%q(1:n)
      w%ad_valid = .TRUE.
    END ASSOCIATE

  CONTAINS

    LOGICAL FUNCTION ensure_r1(x, n1) RESULT(ready)
      !! Give x exactly n1 entries, reallocating only when the length differs.
      REAL(wp), ALLOCATABLE, INTENT(INOUT) :: x(:)
      INTEGER, INTENT(IN) :: n1
      INTEGER :: istat

      ready = .TRUE.
      IF (ALLOCATED(x)) THEN
        IF (SIZE(x) == n1) RETURN
        DEALLOCATE (x)
      END IF
      ALLOCATE (x(n1), STAT=istat)
      ready = istat == 0
    END FUNCTION ensure_r1

    LOGICAL FUNCTION ensure_i1(x, n1) RESULT(ready)
      !! Give x exactly n1 entries, reallocating only when the length differs.
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: x(:)
      INTEGER, INTENT(IN) :: n1
      INTEGER :: istat

      ready = .TRUE.
      IF (ALLOCATED(x)) THEN
        IF (SIZE(x) == n1) RETURN
        DEALLOCATE (x)
      END IF
      ALLOCATE (x(n1), STAT=istat)
      ready = istat == 0
    END FUNCTION ensure_i1

    LOGICAL FUNCTION ensure_r2(x, n1, n2) RESULT(ready)
      !! Give x exactly the shape (n1, n2), reallocating only when the shape differs.
      REAL(wp), ALLOCATABLE, INTENT(INOUT) :: x(:, :)
      INTEGER, INTENT(IN) :: n1, n2
      INTEGER :: istat

      ready = .TRUE.
      IF (ALLOCATED(x)) THEN
        IF (SIZE(x, 1) == n1 .AND. SIZE(x, 2) == n2) RETURN
        DEALLOCATE (x)
      END IF
      ALLOCATE (x(n1, n2), STAT=istat)
      ready = istat == 0
    END FUNCTION ensure_r2
  END SUBROUTINE CD_Update_Model_CoupledAccelDerivative

  SUBROUTINE CD_Calc_Model_CoupledAccelDerivative(model, dload_da, ErrStat, ErrMsg)
    !! Analytic coupled-load derivative with respect to prescribed coupled
    !! accelerations. For free DOFs eliminated by the consistent-mass solve,
    !! d(load_c)/d(a_c) = -(M_cc - M_cf M_ff^{-1} M_fc), with M including added mass
    !! (computed by CD_Update_Model_CoupledAccelDerivative and copied out).
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(OUT) :: dload_da(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: nc

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    dload_da = CD_ZERO
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    nc = SIZE(model%fixed_dofs)
    IF (SIZE(dload_da, 1) /= nc .OR. SIZE(dload_da, 2) /= nc) THEN
      CALL fail(ErrStat, ErrMsg, 'dload_da output must have shape (n_coupled_dof, n_coupled_dof)')
      RETURN
    END IF
    IF (nc == 0) RETURN
    CALL CD_Update_Model_CoupledAccelDerivative(model, ErrStat, ErrMsg)
    IF (ErrStat == CD_MODEL_OK) dload_da = model%dynamic_workspace%ad_result
  END SUBROUTINE CD_Calc_Model_CoupledAccelDerivative

  SUBROUTINE CD_Update_Model_SegmentLength(model, elem, l0_new, l0_dot, ErrStat, ErrMsg)
    !! Active line control (the MoorDyn CtrlChan mechanism): set one element's
    !! unstretched length and its rate. The elastic tension and the per-step mass
    !! assembly read model%l0 directly, so the dynamics track the new length from
    !! the next step; the cached structural mass (the coupled-loads/inertia cache)
    !! is re-assembled here so every consumer sees one consistent length. BA stays
    !! frozen at its init resolution (MoorDyn's convention -- BA derives from the
    !! INITIAL UnstrLen). Callers owning point bindings (the system) must recompute
    !! any end-node mass-diagonal shares after this call.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    INTEGER, INTENT(IN) :: elem
    REAL(wp), INTENT(IN) :: l0_new, l0_dot
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(160) :: em
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (elem < 1 .OR. elem > SIZE(model%l0)) THEN
      CALL fail(ErrStat, ErrMsg, 'segment-length update: element index out of range')
      RETURN
    END IF
    IF (.NOT. (CD_Is_Finite(l0_new) .AND. CD_Is_Finite(l0_dot) .AND. l0_new > CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'segment-length update: l0 must be finite and positive, l0_dot finite')
      RETURN
    END IF
    IF (model%has_seabed .AND. .NOT. model%has_seabed_recompute) THEN
      CALL fail(ErrStat, ErrMsg, 'segment-length update on a seabed-contact model needs the seabed '// &
                'recompute scalars (kBot + contact diameter) supplied at init')
      RETURN
    END IF
    IF (model%has_viscoelastic) THEN
      IF (model%ve_mode(elem) /= 0) THEN
        ! the dl_1 state partitions the stretch against the INIT-time natural
        ! length; MoorDyn defines no rescaling law for a controlled
        ! viscoelastic segment, so this composition fails closed by name
        CALL fail(ErrStat, ErrMsg, 'segment-length update on a viscoelastic (ElasticMod > 1) element '// &
                  'is not supported: the dl_1 state has no length-rescaling law')
        RETURN
      END IF
    END IF
    IF (model%has_syrope) THEN
      IF (model%syrope_is(elem)) THEN
        ! the Syrope strain states partition against the INIT-time natural
        ! length; no rescaling law for a controlled Syrope segment, so this
        ! composition fails closed by name (same rule as the viscoelastic path)
        CALL fail(ErrStat, ErrMsg, 'segment-length update on a Syrope element is not supported: '// &
                  'the working-curve strain states have no length-rescaling law')
        RETURN
      END IF
    END IF
    ! The update is atomic: every length-dependent quantity is built from a trial
    ! length vector first, and the committed fields are mutated only once all builds
    ! succeeded. A failed acceleration recompute restores the entry state whole.
    BLOCK
      REAL(wp) :: l0_trial(SIZE(model%l0)), f_dist(SIZE(model%f_ext)), f_ext_entry(SIZE(model%f_ext))
      REAL(wp) :: l0_entry, l0_dot_entry
      REAL(wp), ALLOCATABLE :: kn_new(:), kn_entry(:), cn_entry(:)
      l0_trial = model%l0
      l0_trial(elem) = l0_new
      ! Distributed external load tracks the length: REASSEMBLE the distributed part
      ! from the current lengths onto the stored non-distributed remainder --
      ! bit-identical to a fresh build with the new length (a nodal delta is not,
      ! in floating point).
      IF (ALLOCATED(model%dist_load)) THEN
        CALL CD_Assemble_Distributed_Load(model%elem_conn, l0_trial, model%dist_load, f_dist, es, em)
        IF (es /= 0) THEN
          CALL fail(ErrStat, ErrMsg, 'segment-length update: distributed-load reassembly failed: '//TRIM(em))
          RETURN
        END IF
      END IF
      ! Seabed nodal coefficients track the tributary lengths: rebuild through the
      ! SAME helper the deck used, so a fresh build at the new length matches
      ! bit-for-bit; cn follows via the stored ratio.
      IF (model%has_seabed .AND. model%has_seabed_recompute) THEN
        CALL CD_Nodal_Seabed_Stiffness(model%seabed_kbot, model%seabed_diameter_elem, l0_trial, &
                                       kn_new, es, em)
        IF (es /= 0) THEN
          CALL fail(ErrStat, ErrMsg, 'segment-length update: seabed stiffness recompute failed: '//TRIM(em))
          RETURN
        END IF
      END IF
      ! Commit. The consistent mass follows l0 implicitly (it is applied element by
      ! element from the committed lengths, never stored).
      l0_entry = model%l0(elem)
      l0_dot_entry = model%l0_dot(elem)
      f_ext_entry = model%f_ext
      model%l0(elem) = l0_new
      model%l0_dot(elem) = l0_dot
      IF (ALLOCATED(model%dist_load)) model%f_ext = model%f_ext_dist_base + f_dist
      IF (model%has_seabed .AND. model%has_seabed_recompute) THEN
        CALL MOVE_ALLOC(model%seabed_kn, kn_entry)
        CALL MOVE_ALLOC(kn_new, model%seabed_kn)
        IF (ALLOCATED(model%seabed_cn)) THEN
          cn_entry = model%seabed_cn
          model%seabed_cn = model%seabed_kn*model%seabed_cn_over_kn
        END IF
      END IF
      ! The committed acceleration satisfied the OLD length's equation of motion;
      ! recompute it at the current q/v under the new length so the model is
      ! indistinguishable from one built with it (load exchanges and the next
      ! gen-alpha step start from a consistent state). A FULLY COUPLED model (no
      ! free DOFs -- every line in a point-connected system) carries only
      ! caller-prescribed accelerations, so there is nothing to recompute.
      IF (SIZE(model%fixed_dofs) < model%n_dof) THEN
        CALL CD_Recompute_Model_Acceleration(model, es, em)
        IF (es /= CD_MODEL_OK) THEN
          ! The recompute restored the acceleration; restore every other field.
          model%l0(elem) = l0_entry
          model%l0_dot(elem) = l0_dot_entry
          model%f_ext = f_ext_entry
          IF (ALLOCATED(kn_entry)) CALL MOVE_ALLOC(kn_entry, model%seabed_kn)
          IF (ALLOCATED(cn_entry)) model%seabed_cn = cn_entry
          CALL fail(ErrStat, ErrMsg, 'segment-length update: acceleration recompute failed: '//TRIM(em))
        END IF
      END IF
    END BLOCK
  END SUBROUTINE CD_Update_Model_SegmentLength

  SUBROUTINE CD_Get_Model_EndNodeMassDiag(model, coupled_slot, mass_diag, ErrStat, ErrMsg)
    !! The coupled end node's OWN diagonal share of the consistent structural mass:
    !! sum over elements attached to that node of rho_a*l0/3 (the 2m/6 block diagonal
    !! of CD_Assemble_Cable_Mass). This is the piece of the end-node inertia a point
    !! integrator can treat IMPLICITLY (added to the point's effective mass, with its
    !! lagged product removed from the coupled-loads reaction) so a light or massless
    !! point with attached lines integrates stably; the m/6 neighbour cross term stays
    !! in the lagged reaction. coupled_slot selects the endpoint triple in the model's
    !! compact coupled ordering (1 = fixed_dofs(1:3), 2 = fixed_dofs(4:6)), matching
    !! CD_LINE_END_A/B in the system point map.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: coupled_slot
    REAL(wp), INTENT(OUT) :: mass_diag
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: e, node, d0
    mass_diag = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (coupled_slot /= 1 .AND. coupled_slot /= 2) THEN
      CALL fail(ErrStat, ErrMsg, 'coupled_slot must be 1 (end A triple) or 2 (end B triple)')
      RETURN
    END IF
    IF (SIZE(model%fixed_dofs) < 3*coupled_slot) THEN
      CALL fail(ErrStat, ErrMsg, 'model does not expose the requested coupled endpoint triple')
      RETURN
    END IF
    d0 = model%fixed_dofs(3*coupled_slot - 2)
    IF (MOD(d0 - 1, 3) /= 0 .OR. &
        model%fixed_dofs(3*coupled_slot - 1) /= d0 + 1 .OR. &
        model%fixed_dofs(3*coupled_slot) /= d0 + 2) THEN
      CALL fail(ErrStat, ErrMsg, 'the requested coupled slot is not a contiguous per-node triple')
      RETURN
    END IF
    node = (d0 + 2)/3
    DO e = 1, SIZE(model%elem_conn, 2)
      IF (model%elem_conn(1, e) == node .OR. model%elem_conn(2, e) == node) THEN
        mass_diag = mass_diag + model%rho_a(e)*model%l0(e)/3.0_wp
      END IF
    END DO
  END SUBROUTINE CD_Get_Model_EndNodeMassDiag

  SUBROUTINE CD_Get_Model_EndNodeAddedMass(model, coupled_slot, ma_block, ErrStat, ErrMsg)
    !! The coupled end node's OWN 3x3 block of the configuration-dependent Morison
    !! added-mass matrix at the current model%q: the sum over elements attached to
    !! that node of the node's self block of CD_Added_Mass_Element_Matrix (the same
    !! element matrix, wetted fraction and orientation whose product M_add(q)*a enters
    !! CD_Calc_Model_CoupledLoads). A point integrator moves this block into its
    !! implicit mass together with the structural diagonal, cancelling the lagged
    !! product it carries in the coupled-loads reaction. Zero when the model has no
    !! added mass. coupled_slot follows CD_Get_Model_EndNodeMassDiag.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: coupled_slot
    REAL(wp), INTENT(OUT) :: ma_block(3, 3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: e, node, d0, idx(6), es
    REAL(wp) :: q2(6), wl2(2), m2(6, 6)
    CHARACTER(160) :: em
    ma_block = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (coupled_slot /= 1 .AND. coupled_slot /= 2) THEN
      CALL fail(ErrStat, ErrMsg, 'coupled_slot must be 1 (end A triple) or 2 (end B triple)')
      RETURN
    END IF
    IF (SIZE(model%fixed_dofs) < 3*coupled_slot) THEN
      CALL fail(ErrStat, ErrMsg, 'model does not expose the requested coupled endpoint triple')
      RETURN
    END IF
    d0 = model%fixed_dofs(3*coupled_slot - 2)
    IF (MOD(d0 - 1, 3) /= 0 .OR. &
        model%fixed_dofs(3*coupled_slot - 1) /= d0 + 1 .OR. &
        model%fixed_dofs(3*coupled_slot) /= d0 + 2) THEN
      CALL fail(ErrStat, ErrMsg, 'the requested coupled slot is not a contiguous per-node triple')
      RETURN
    END IF
    IF (.NOT. model%has_added_mass) RETURN
    node = (d0 + 2)/3
    DO e = 1, model%n_elem
      IF (model%elem_conn(1, e) /= node .AND. model%elem_conn(2, e) /= node) CYCLE
      CALL extract_element_positions(model, e, model%q, idx, q2)
      CALL extract_nodal_scalar(model%added_mass_waterline_z, model, e, wl2)
      CALL CD_Added_Mass_Element_Matrix(q2(1:3), q2(4:6), model%l0(e), wl2(1), wl2(2), &
                                        model%added_mass_rho, model%added_mass_diameter_elem(e), &
                                        model%added_mass_can_elem(e), model%added_mass_cat_elem(e), &
                                        m2, es, em)
      IF (es /= 0) THEN
        ErrStat = CD_MODEL_SOLVEFAIL
        ErrMsg = 'CableDyn_Model: end-node added-mass block failed: '//TRIM(em)
        ma_block = CD_ZERO
        RETURN
      END IF
      ! An element never joins a node to itself, so each attached element
      ! contributes exactly one of its two nodal self blocks.
      IF (model%elem_conn(1, e) == node) THEN
        ma_block = ma_block + m2(1:3, 1:3)
      ELSE
        ma_block = ma_block + m2(4:6, 4:6)
      END IF
    END DO
    IF (.NOT. CD_All_Finite(ma_block)) THEN
      CALL fail(ErrStat, ErrMsg, 'end-node added-mass block is non-finite')
      ma_block = CD_ZERO
    END IF
  END SUBROUTINE CD_Get_Model_EndNodeAddedMass

  SUBROUTINE CD_Set_Model_EndQuery(model, request, active, arm_guess, refresh, levels, level_dofs)
    !! Arm (request = .TRUE.) or disarm the junction query of the next CD_Step_Model: a
    !! converged step then leaves the alpha-level reaction rows of its fixed DOFs and the
    !! condensed dynamic end tangent (CD_Get_Model_EndQuery). active (one flag per fixed DOF,
    !! in model%fixed_dofs order) limits the tangent columns to the junction DOFs; arm_guess
    !! warm-starts the step from the previous converged solve of the same step.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    LOGICAL, INTENT(IN) :: request
    LOGICAL, INTENT(IN), OPTIONAL :: active(:)
    LOGICAL, INTENT(IN), OPTIONAL :: arm_guess
    !! refresh = .FALSE. keeps the condensed tangent of an earlier solve when one is valid
    LOGICAL, INTENT(IN), OPTIONAL :: refresh
    !! levels = [beta, gamma, alpha_m] of the objects that drive the fixed DOFs: their Newmark
    !! and inertia level replace the line's at those DOFs (the monolithic junction)
    REAL(wp), INTENT(IN), OPTIONAL :: levels(3)
    !! level_dofs (one flag per fixed DOF, in model%fixed_dofs order): the DOFs that take levels
    !! (all when absent); the others keep the line's own
    LOGICAL, INTENT(IN), OPTIONAL :: level_dofs(:)
    model%dynamic_workspace%end_request = request
    model%dynamic_workspace%end_levels = PRESENT(levels) .AND. request
    IF (PRESENT(levels)) THEN
      model%dynamic_workspace%end_beta = levels(1)
      model%dynamic_workspace%end_gamma = levels(2)
      model%dynamic_workspace%end_alpham = levels(3)
    END IF
    ! Reuse the array while its size is unchanged (one call per junction iteration).
    IF (PRESENT(level_dofs)) THEN
      IF (ALLOCATED(model%dynamic_workspace%end_level_dofs)) THEN
        IF (SIZE(model%dynamic_workspace%end_level_dofs) /= SIZE(level_dofs)) &
          DEALLOCATE (model%dynamic_workspace%end_level_dofs)
      END IF
      IF (.NOT. ALLOCATED(model%dynamic_workspace%end_level_dofs)) &
        ALLOCATE (model%dynamic_workspace%end_level_dofs(SIZE(level_dofs)))
      model%dynamic_workspace%end_level_dofs = level_dofs
    ELSE IF (ALLOCATED(model%dynamic_workspace%end_level_dofs)) THEN
      DEALLOCATE (model%dynamic_workspace%end_level_dofs)
    END IF
    model%dynamic_workspace%end_refresh = .TRUE.
    IF (PRESENT(refresh)) model%dynamic_workspace%end_refresh = refresh
    model%dynamic_workspace%guess_armed = .FALSE.
    IF (PRESENT(arm_guess)) model%dynamic_workspace%guess_armed = arm_guess .AND. request
    IF (PRESENT(active)) THEN
      IF (ALLOCATED(model%dynamic_workspace%end_active)) THEN
        IF (SIZE(model%dynamic_workspace%end_active) /= SIZE(active)) &
          DEALLOCATE (model%dynamic_workspace%end_active)
      END IF
      IF (.NOT. ALLOCATED(model%dynamic_workspace%end_active)) &
        ALLOCATE (model%dynamic_workspace%end_active(SIZE(active)))
      model%dynamic_workspace%end_active = active
    END IF
  END SUBROUTINE CD_Set_Model_EndQuery

  SUBROUTINE CD_Get_Model_EndQuery(model, reaction, tangent, valid)
    !! The junction query of the last CD_Step_Model (see CD_Set_Model_EndQuery): reaction(j)
    !! is the alpha-level residual row of fixed DOF j (the force the support applies to the
    !! line; the line loads the support with -reaction) and tangent(i, j) = d reaction_i /
    !! d q_{n+1, j} with the free DOFs condensed. valid is false when the last step did not
    !! produce them (not armed, not converged, or a dense-path step).
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: reaction(:), tangent(:, :)
    LOGICAL, INTENT(OUT) :: valid
    INTEGER :: n
    reaction = CD_ZERO
    tangent = CD_ZERO
    valid = model%dynamic_workspace%end_valid
    IF (.NOT. valid) RETURN
    n = SIZE(model%dynamic_workspace%end_r)
    IF (SIZE(reaction) /= n .OR. SIZE(tangent, 1) /= n .OR. SIZE(tangent, 2) /= n) THEN
      valid = .FALSE.
      RETURN
    END IF
    reaction = model%dynamic_workspace%end_r
    tangent = model%dynamic_workspace%end_k
  END SUBROUTINE CD_Get_Model_EndQuery

  SUBROUTINE CD_Get_Model_EndNodeTangent(model, coupled_slot, k_block, c_block, ErrStat, ErrMsg, kv_neighbour)
    !! The coupled end node's OWN 3x3 stiffness and damping blocks at the current
    !! model%q/v, summed over the elements attached to that node:
    !!   k_block = sum (EA_dyn/l0) t t^T + (max(T, 0)/l) (I - t t^T)   (elastic tangent),
    !!   c_block = sum (BA_dyn/l0) t t^T                             (axial damping).
    !! These are the self blocks of the element tangent and axial-damping Jacobian with
    !! the geometric term clamped to tension, so both blocks are symmetric positive
    !! semi-definite. A point integrator uses them to advance the end-segment spring
    !! and damper linearly implicitly; the reaction itself is still the exact coupled
    !! load, so the blocks only shape the implicit correction. The end node's own seabed
    !! contact is added: the normal penalty tangent to k_block and, while the node moves
    !! down, the normal damping to c_block, so a free line end resting on a stiff seabed
    !! is not limited by the explicit contact frequency. Viscoelastic, Syrope, drag and
    !! seabed friction contributions are not included (they stay in the reaction at the
    !! step-start state). coupled_slot follows CD_Get_Model_EndNodeMassDiag.
    !! kv_neighbour (optional) is sum K_e v_o over the same elements, v_o the velocity of each
    !! element's other node at model%v: the cross-block (-K_e) coupling that carries the
    !! neighbours along with the end node, so an implicit end-node update springs against the
    !! relative motion (v_end - v_o) rather than against neighbours held in place.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: coupled_slot
    REAL(wp), INTENT(OUT) :: k_block(3, 3), c_block(3, 3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(OUT), OPTIONAL :: kv_neighbour(3)
    INTEGER :: e, node, d0, idx(6), i
    REAL(wp) :: q2(6), chord(3), ell, tang(3), tt(3, 3), strain, material, tension, ke(3, 3), v_other(3)
    REAL(wp) :: z_floor, dfdx, dfdy, gap, fvec(3), jac3(3, 3), damp_force, damp_dgap, damp_dvz
    k_block = CD_ZERO
    c_block = CD_ZERO
    IF (PRESENT(kv_neighbour)) kv_neighbour = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (coupled_slot /= 1 .AND. coupled_slot /= 2) THEN
      CALL fail(ErrStat, ErrMsg, 'coupled_slot must be 1 (end A triple) or 2 (end B triple)')
      RETURN
    END IF
    IF (SIZE(model%fixed_dofs) < 3*coupled_slot) THEN
      CALL fail(ErrStat, ErrMsg, 'model does not expose the requested coupled endpoint triple')
      RETURN
    END IF
    d0 = model%fixed_dofs(3*coupled_slot - 2)
    IF (MOD(d0 - 1, 3) /= 0 .OR. &
        model%fixed_dofs(3*coupled_slot - 1) /= d0 + 1 .OR. &
        model%fixed_dofs(3*coupled_slot) /= d0 + 2) THEN
      CALL fail(ErrStat, ErrMsg, 'the requested coupled slot is not a contiguous per-node triple')
      RETURN
    END IF
    node = (d0 + 2)/3
    DO e = 1, model%n_elem
      IF (model%elem_conn(1, e) /= node .AND. model%elem_conn(2, e) /= node) CYCLE
      CALL extract_element_positions(model, e, model%q, idx, q2)
      chord = q2(4:6) - q2(1:3)
      ell = NORM2(chord)
      IF (.NOT. (CD_Is_Finite(ell) .AND. ell > CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'end-node tangent: attached element collapsed or non-finite')
        k_block = CD_ZERO
        c_block = CD_ZERO
        RETURN
      END IF
      tang = chord/ell
      DO i = 1, 3
        tt(:, i) = tang*tang(i)
      END DO
      strain = ell/model%l0(e) - CD_ONE
      material = model%ea_dyn(e)/model%l0(e)
      tension = model%ea_dyn(e)*strain
      IF (model%tension_only .AND. strain < CD_ZERO) material = CD_ZERO
      ke = material*tt + (MAX(tension, CD_ZERO)/ell)*(identity3() - tt)
      k_block = k_block + ke
      IF (PRESENT(kv_neighbour)) THEN
        IF (model%elem_conn(1, e) == node) THEN
          v_other = model%v(idx(4:6))
        ELSE
          v_other = model%v(idx(1:3))
        END IF
        kv_neighbour = kv_neighbour + MATMUL(ke, v_other)
      END IF
      IF (model%has_damping) c_block = c_block + (MAX(model%ba_dyn(e), CD_ZERO)/model%l0(e))*tt
    END DO
    IF (model%has_seabed) THEN
      ! the end node's seabed contact (its neighbour is the fixed seabed, so no kv term)
      IF (model%has_bathymetry) THEN
        CALL model_floor_gradient(model, model%q(d0), model%q(d0 + 1), z_floor, dfdx, dfdy, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) THEN
          k_block = CD_ZERO
          c_block = CD_ZERO
          IF (PRESENT(kv_neighbour)) kv_neighbour = CD_ZERO
          RETURN
        END IF
      ELSE
        z_floor = model%seabed_z_floor
        dfdx = CD_ZERO
        dfdy = CD_ZERO
      END IF
      gap = z_floor - model%q(d0 + 2)
      CALL CD_Seabed_Normal_Contact(gap, dfdx, dfdy, model%seabed_kn(node), fvec, jac3)
      k_block = k_block + jac3
      IF (model%has_seabed_damping) THEN
        BLOCK
          REAL(wp) :: gap_n, vn, nvec(3)
          LOGICAL :: sloped
          INTEGER :: cb
          CALL node_contact_normal(model, model%q(d0), model%q(d0 + 1), model%q(d0 + 2), model%v(d0:d0 + 2), &
                                   gap_n, vn, nvec, sloped, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) THEN
            k_block = CD_ZERO
            c_block = CD_ZERO
            IF (PRESENT(kv_neighbour)) kv_neighbour = CD_ZERO
            RETURN
          END IF
          CALL seabed_normal_damping_terms(gap_n, vn, model%seabed_cn(node), damp_force, damp_dgap, damp_dvz)
          IF (sloped) THEN
            DO cb = 1, 3
              c_block(:, cb) = c_block(:, cb) - damp_dvz*nvec*nvec(cb)
            END DO
          ELSE
            c_block(3, 3) = c_block(3, 3) - damp_dvz
          END IF
        END BLOCK
      END IF
    END IF
    IF (PRESENT(kv_neighbour)) THEN
      IF (.NOT. CD_All_Finite(kv_neighbour)) THEN
        CALL fail(ErrStat, ErrMsg, 'end-node neighbour coupling is non-finite')
        k_block = CD_ZERO
        c_block = CD_ZERO
        kv_neighbour = CD_ZERO
        RETURN
      END IF
    END IF
    IF (.NOT. (CD_All_Finite(k_block) .AND. CD_All_Finite(c_block))) THEN
      CALL fail(ErrStat, ErrMsg, 'end-node tangent blocks are non-finite')
      k_block = CD_ZERO
      c_block = CD_ZERO
      IF (PRESENT(kv_neighbour)) kv_neighbour = CD_ZERO
    END IF
  CONTAINS
    PURE FUNCTION identity3() RESULT(eye)
      REAL(wp) :: eye(3, 3)
      INTEGER :: j
      eye = CD_ZERO
      DO j = 1, 3
        eye(j, j) = CD_ONE
      END DO
    END FUNCTION identity3
  END SUBROUTINE CD_Get_Model_EndNodeTangent

  SUBROUTINE CD_Get_Model_EndForces(model, force_first, force_last, ErrStat, ErrMsg)
    !! Force the line exerts on the points holding its first and last node at the
    !! committed state (q, v): the end element's internal force (elastic, axial damping and
    !! any constitutive contributor) plus the end node's lumped external loads (submerged
    !! weight, buoyancy, seabed contact with its damping and friction, drag at the actual
    !! relative velocity, Froude-Krylov), i.e. f_ext + contributors - f_int at the node.
    !! This is the actual force on the attachment, as MoorDyn and OrcaFlex report the end
    !! tension; it equals the coupled load of the node (CD_Calc_Model_CoupledLoads) less the
    !! end node's structural and added-mass inertia. At rest it is the static initializer's
    !! end force. Its norm is the FairTen/AnchTen channel value.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: force_first(3), force_last(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp) :: load(9)
    REAL(wp) :: fint_first(3), fint_last(3), nodes6(6), Kt6(6, 6), fint6(6), Te
    INTEGER :: es, n, nn, e, a, b
    CHARACTER(160) :: em

    force_first = CD_ZERO
    force_last = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    n = model%n_dof
    nn = n/3
    ! Only the two end nodes are read, so only the elements that touch them are evaluated.
    ! Each end-node sum takes its terms in the element order of the whole-line assembly,
    ! so the result is bit-identical to assembling the full internal and load vectors.
    fint_first = CD_ZERO
    fint_last = CD_ZERO
    DO e = 1, model%n_elem
      a = model%elem_conn(1, e)
      b = model%elem_conn(2, e)
      IF (a /= 1 .AND. a /= nn .AND. b /= 1 .AND. b /= nn) CYCLE
      nodes6(1:3) = model%q(3*a - 2:3*a)
      nodes6(4:6) = model%q(3*b - 2:3*b)
      CALL CD_Compute_Cable_Element(nodes6, model%ea_dyn(e), model%l0(e), model%tension_only, Kt6, fint6, &
                                    Te, es, em)
      IF (es /= 0) THEN
        ErrStat = CD_MODEL_SOLVEFAIL
        ErrMsg = 'CableDyn_Model: internal-force assembly failed during end-force recovery: '// &
                 'CD_Assemble_Cable_Internal_Force: element '//TRIM(em)
        RETURN
      END IF
      IF (a == 1) fint_first = fint_first + fint6(1:3)
      IF (b == 1) fint_first = fint_first + fint6(4:6)
      IF (a == nn) fint_last = fint_last + fint6(1:3)
      IF (b == nn) fint_last = fint_last + fint6(4:6)
    END DO
    ! compact end-node load [first node, last node, other reached nodes]: no heap scratch
    CALL model_contributor_force(model, model%q, model%v, load, es, em, end_nodes_only=.TRUE.)
    IF (es /= 0) THEN
      ErrStat = CD_MODEL_SOLVEFAIL
      ErrMsg = 'CableDyn_Model: load assembly failed during end-force recovery: '//TRIM(em)
      RETURN
    END IF
    force_first = model%f_ext(1:3) + load(1:3) - fint_first
    force_last = model%f_ext(n - 2:n) + load(4:6) - fint_last
  END SUBROUTINE CD_Get_Model_EndForces

  SUBROUTINE CD_Get_Model_Tension(model, tension, ErrStat, ErrMsg)
    !! Recover current element axial tension from the model state.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: tension(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: es
    CHARACTER(160) :: em

    tension = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (SIZE(tension) /= model%n_elem) THEN
      CALL fail(ErrStat, ErrMsg, 'tension output must have shape (n_elem)')
      RETURN
    END IF
    CALL CD_Compute_Cable_Tension(RESHAPE(model%q, [3, model%n_dof/3]), model%elem_conn, model%l0, model%ea, &
                                  model%tension_only, tension, es, em, topology_validated=.TRUE.)
    IF (es /= 0) THEN
      ErrStat = CD_MODEL_SOLVEFAIL
      ErrMsg = 'CableDyn_Model: tension recovery failed: '//TRIM(em)
      RETURN
    END IF
    IF (model%has_damping) THEN
      ! The axial Kelvin-Voigt share BA*(strain rate) is part of the reported tension, as it
      ! is of the end forces (CD_Get_Model_EndForces): the same element law, projected on the
      ! chord. It also acts on a slack element. Viscoelastic and Syrope elements have
      ! ba_dyn = 0 here and report their own rate terms below.
      BLOCK
        INTEGER :: e, na, nb
        REAL(wp) :: f2(6), chord(3), clen
        DO e = 1, model%n_elem
          IF (.NOT. (model%ba_dyn(e) > CD_ZERO)) CYCLE
          na = model%elem_conn(1, e)
          nb = model%elem_conn(2, e)
          CALL CD_Axial_Damping_Element_Force(model%q(3*na - 2:3*na), model%q(3*nb - 2:3*nb), &
                                              model%v(3*na - 2:3*na), model%v(3*nb - 2:3*nb), &
                                              model%l0(e), model%ba_dyn(e), f2, es, em, &
                                              l0_dot=model%l0_dot(e))
          IF (es /= 0) THEN
            ErrStat = CD_MODEL_SOLVEFAIL
            ErrMsg = 'CableDyn_Model: axial damping tension recovery failed: '//TRIM(em)
            RETURN
          END IF
          chord = model%q(3*nb - 2:3*nb) - model%q(3*na - 2:3*na)
          clen = NORM2(chord)
          tension(e) = tension(e) + DOT_PRODUCT(f2(1:3), chord/clen)
        END DO
      END BLOCK
    END IF
    IF (model%has_viscoelastic) THEN
      ! viscoelastic elements report the series-Kelvin tension MagT + MagTd at the
      ! COMMITTED state and dl_1 (the elastic EA recovery above does not apply
      ! to them -- their dynamic elastic stiffness lives in the state)
      BLOCK
        INTEGER :: e, na, nb
        REAL(wp) :: ea_d_eff, ea_1_eff
        DO e = 1, model%n_elem
          IF (model%ve_mode(e) == 0) CYCLE
          CALL model_viscoelastic_coeffs(model, e, ea_d_eff, ea_1_eff, es, em)
          IF (es /= CD_MODEL_OK) THEN
            ErrStat = es
            ErrMsg = em
            RETURN
          END IF
          na = model%elem_conn(1, e)
          nb = model%elem_conn(2, e)
          CALL CD_Viscoelastic_Element_Tension(model%q(3*na - 2:3*na), model%q(3*nb - 2:3*nb), &
                                               model%v(3*na - 2:3*na), model%v(3*nb - 2:3*nb), &
                                               model%l0(e), ea_d_eff, ea_1_eff, &
                                               model%ve_ba(e), model%ve_ba_d(e), model%ve_dl_1(e), &
                                               tension(e), es, em)
          IF (es /= CD_VISCO_OK) THEN
            ErrStat = CD_MODEL_SOLVEFAIL
            ErrMsg = 'CableDyn_Model: viscoelastic tension recovery failed: '//TRIM(em)
            RETURN
          END IF
        END DO
      END BLOCK
    END IF
    IF (model%has_syrope) THEN
      ! Syrope elements report the working-curve tension at the COMMITTED state
      ! (the along-tangent component of the element force; the elastic EA
      ! recovery above does not apply -- their stiffness lives in the states)
      BLOCK
        INTEGER :: e, idx_scr(6)
        REAL(wp) :: q2b(6), v2b(6), f6(6), j6(6, 6), jv6(6, 6), tmean, dslow, chord(3)
        DO e = 1, model%n_elem
          IF (.NOT. model%syrope_is(e)) CYCLE
          CALL extract_element_state(model, e, model%q, model%v, idx_scr, q2b, v2b)
          CALL model_syrope_element(model, e, q2b, v2b, f6, j6, jv6, tmean, dslow, es, em)
          IF (es /= CD_MODEL_OK) THEN
            ErrStat = es
            ErrMsg = em
            RETURN
          END IF
          chord = q2b(4:6) - q2b(1:3)
          tension(e) = DOT_PRODUCT(f6(1:3), chord/NORM2(chord))
        END DO
      END BLOCK
    END IF
  END SUBROUTINE CD_Get_Model_Tension

  SUBROUTINE CD_Update_Model_Hydro_Fields(model, ErrStat, ErrMsg, fluid_velocity, drag_waterline_z, &
                                          fluid_acceleration, fk_waterline_z, added_mass_waterline_z, &
                                          buoyancy_waterline_z, fields_changed)
    !! Refresh nodal environmental fields owned by an initialized model. This is
    !! the small time-varying hydro boundary used by the deck drivers and
    !! OpenFAST shells: coefficients are fixed at initialization, while current,
    !! wave acceleration, free surface, and wetting fields can change every step.
    !! This routine does not recompute acceleration. Inside a time march the
    !! committed acceleration is the generalized-alpha ALGORITHMIC acceleration and
    !! must be preserved: replacing it with the equilibrium acceleration implied by
    !! the new fields destroys the scheme's high-frequency dissipation. Only an
    !! initial-condition station (before the first step) should follow a changed
    !! field with CD_Recompute_Model_Acceleration; fields_changed reports whether any
    !! supplied field differs from the stored one, so an unchanged (e.g. zero) field
    !! leaves the model bit-identical to one that was never refreshed.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: fluid_velocity(:, :)
    REAL(wp), INTENT(IN), OPTIONAL :: drag_waterline_z(:)
    REAL(wp), INTENT(IN), OPTIONAL :: fluid_acceleration(:, :)
    REAL(wp), INTENT(IN), OPTIONAL :: fk_waterline_z(:)
    REAL(wp), INTENT(IN), OPTIONAL :: added_mass_waterline_z(:)
    REAL(wp), INTENT(IN), OPTIONAL :: buoyancy_waterline_z(:)
    LOGICAL, INTENT(OUT), OPTIONAL :: fields_changed

    INTEGER :: n_nodes
    LOGICAL :: changed

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (PRESENT(fields_changed)) fields_changed = .FALSE.
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    n_nodes = model%n_dof/3
    IF (PRESENT(fluid_velocity)) THEN
      IF (.NOT. model%has_morison_drag) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_velocity update requires Morison drag configured at initialization')
        RETURN
      END IF
      IF (SIZE(fluid_velocity, 1) /= 3 .OR. SIZE(fluid_velocity, 2) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_velocity update must have shape (3, n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(fluid_velocity)) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_velocity update must be finite')
        RETURN
      END IF
    END IF
    IF (PRESENT(drag_waterline_z)) THEN
      IF (.NOT. model%has_morison_drag) THEN
        CALL fail(ErrStat, ErrMsg, 'drag_waterline_z update requires Morison drag configured at initialization')
        RETURN
      END IF
      CALL validate_vector_field(drag_waterline_z, n_nodes, 'drag_waterline_z', ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
    END IF
    IF (PRESENT(fluid_acceleration)) THEN
      IF (.NOT. model%has_froude_krylov) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_acceleration update requires Froude-Krylov configured at initialization')
        RETURN
      END IF
      IF (SIZE(fluid_acceleration, 1) /= 3 .OR. SIZE(fluid_acceleration, 2) /= n_nodes) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_acceleration update must have shape (3, n_nodes)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(fluid_acceleration)) THEN
        CALL fail(ErrStat, ErrMsg, 'fluid_acceleration update must be finite')
        RETURN
      END IF
    END IF
    IF (PRESENT(fk_waterline_z)) THEN
      IF (.NOT. model%has_froude_krylov) THEN
        CALL fail(ErrStat, ErrMsg, 'fk_waterline_z update requires Froude-Krylov configured at initialization')
        RETURN
      END IF
      CALL validate_vector_field(fk_waterline_z, n_nodes, 'fk_waterline_z', ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
    END IF
    IF (PRESENT(added_mass_waterline_z)) THEN
      IF (.NOT. model%has_added_mass) THEN
        CALL fail(ErrStat, ErrMsg, 'added_mass_waterline_z update requires added mass configured at initialization')
        RETURN
      END IF
      CALL validate_vector_field(added_mass_waterline_z, n_nodes, 'added_mass_waterline_z', ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
    END IF
    IF (PRESENT(buoyancy_waterline_z)) THEN
      IF (.NOT. model%has_buoyancy_recovery) THEN
        CALL fail(ErrStat, ErrMsg, &
                  'buoyancy_waterline_z update requires buoyancy recovery configured at initialization')
        RETURN
      END IF
      CALL validate_vector_field(buoyancy_waterline_z, n_nodes, 'buoyancy_waterline_z', ErrStat, ErrMsg)
    END IF
    IF (ErrStat /= CD_MODEL_OK) RETURN

    IF (PRESENT(fields_changed)) THEN
      changed = .FALSE.
      IF (PRESENT(fluid_velocity)) changed = changed .OR. field2_differs(model%fluid_velocity, fluid_velocity)
      IF (PRESENT(drag_waterline_z)) changed = changed .OR. field1_differs(model%drag_waterline_z, drag_waterline_z)
      IF (PRESENT(fluid_acceleration)) &
        changed = changed .OR. field2_differs(model%fluid_acceleration, fluid_acceleration)
      IF (PRESENT(fk_waterline_z)) changed = changed .OR. field1_differs(model%fk_waterline_z, fk_waterline_z)
      IF (PRESENT(added_mass_waterline_z)) &
        changed = changed .OR. field1_differs(model%added_mass_waterline_z, added_mass_waterline_z)
      IF (PRESENT(buoyancy_waterline_z)) &
        changed = changed .OR. field1_differs(model%buoyancy_waterline_z, buoyancy_waterline_z)
      fields_changed = changed
    END IF
    IF (PRESENT(fluid_velocity)) model%fluid_velocity = fluid_velocity
    IF (PRESENT(drag_waterline_z)) model%drag_waterline_z = drag_waterline_z
    IF (PRESENT(fluid_acceleration)) model%fluid_acceleration = fluid_acceleration
    IF (PRESENT(fk_waterline_z)) model%fk_waterline_z = fk_waterline_z
    IF (PRESENT(added_mass_waterline_z)) THEN
      ! The added mass of a partially wet element depends on the waterline: a moved one
      ! invalidates the cached acceleration derivative (keyed on q, l0, rho_a, fixed DOFs).
      IF (field1_differs(model%added_mass_waterline_z, added_mass_waterline_z)) &
        model%dynamic_workspace%ad_valid = .FALSE.
      model%added_mass_waterline_z = added_mass_waterline_z
    END IF
    IF (PRESENT(buoyancy_waterline_z)) model%buoyancy_waterline_z = buoyancy_waterline_z

  CONTAINS

    LOGICAL FUNCTION field1_differs(stored, new) RESULT(differs)
      !! An unallocated or differently shaped store counts as a change.
      REAL(wp), ALLOCATABLE, INTENT(IN) :: stored(:)
      REAL(wp), INTENT(IN) :: new(:)
      differs = .TRUE.
      IF (.NOT. ALLOCATED(stored)) RETURN
      IF (SIZE(stored) /= SIZE(new)) RETURN
      differs = ANY(ABS(stored - new) > CD_ZERO)
    END FUNCTION field1_differs

    LOGICAL FUNCTION field2_differs(stored, new) RESULT(differs)
      REAL(wp), ALLOCATABLE, INTENT(IN) :: stored(:, :)
      REAL(wp), INTENT(IN) :: new(:, :)
      differs = .TRUE.
      IF (.NOT. ALLOCATED(stored)) RETURN
      IF (SIZE(stored, 1) /= SIZE(new, 1) .OR. SIZE(stored, 2) /= SIZE(new, 2)) RETURN
      differs = ANY(ABS(stored - new) > CD_ZERO)
    END FUNCTION field2_differs
  END SUBROUTINE CD_Update_Model_Hydro_Fields

  SUBROUTINE CD_Update_Model_External_Loads(model, f_ext, ErrStat, ErrMsg)
    !! Replace the model-owned nodal external force vector. The update is consumed
    !! by subsequent coupled-load extraction and dynamic steps. The current
    !! acceleration is recomputed before return so direct load queries see a
    !! state consistent with the updated external force vector.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: f_ext(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (SIZE(f_ext) /= model%n_dof) THEN
      CALL fail(ErrStat, ErrMsg, 'f_ext update must have shape (3 n_nodes)')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(f_ext)) THEN
      CALL fail(ErrStat, ErrMsg, 'f_ext update must be finite')
      RETURN
    END IF
    IF (.NOT. (ALLOCATED(model%f_ext_rollback) .AND. ALLOCATED(model%a_rollback))) THEN
      CALL fail(ErrStat, ErrMsg, 'external-load update requires initialized rollback workspace')
      RETURN
    END IF
    model%f_ext_rollback = model%f_ext
    model%a_rollback = model%a
    model%f_ext = f_ext
    CALL CD_Recompute_Model_Acceleration(model, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) THEN
      model%f_ext = model%f_ext_rollback
      model%a = model%a_rollback
      RETURN
    END IF
    ! Maintain the distributed-load invariant f_ext = remainder + assemble(dist, l0):
    ! a later segment-length update RECONSTRUCTS f_ext from the remainder, so a
    ! runtime external-load update that skipped this would be silently WIPED by the
    ! next control command. Recomputed only on success (the rollback path above
    ! restores the committed pair untouched).
    IF (ALLOCATED(model%dist_load)) THEN
      BLOCK
        REAL(wp) :: f_dist(SIZE(model%f_ext))
        INTEGER :: es2
        CHARACTER(160) :: em2
        CALL CD_Assemble_Distributed_Load(model%elem_conn, model%l0, model%dist_load, f_dist, es2, em2)
        IF (es2 /= 0) THEN
          model%f_ext = model%f_ext_rollback
          model%a = model%a_rollback
          CALL fail(ErrStat, ErrMsg, 'external-load update: distributed-load remainder failed: '//TRIM(em2))
          RETURN
        END IF
        model%f_ext_dist_base = model%f_ext - f_dist
      END BLOCK
    END IF
  END SUBROUTINE CD_Update_Model_External_Loads

  SUBROUTINE CD_Update_Model_CoupledMotion(model, q_coupled, v_coupled, a_coupled, ErrStat, ErrMsg)
    !! Refresh the model-owned state on externally coupled DOFs without advancing
    !! the internal dynamics. The input order matches CD_Get_Model_CoupledDofs.
    !! This is the light kinematics-in half of the shell boundary used before
    !! CD_Calc_Model_CoupledLoads in loose-coupling or predictor workflows.
    !! Call CD_Recompute_Model_Acceleration after this update when support
    !! acceleration consistency is required before output/load derivatives.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: q_coupled(:), v_coupled(:), a_coupled(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: n

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    n = SIZE(model%fixed_dofs)
    IF (SIZE(q_coupled) /= n .OR. SIZE(v_coupled) /= n .OR. SIZE(a_coupled) /= n) THEN
      CALL fail(ErrStat, ErrMsg, 'coupled motion updates must have shape (n_coupled_dof)')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(q_coupled) .OR. .NOT. CD_All_Finite(v_coupled) .OR. &
        .NOT. CD_All_Finite(a_coupled)) THEN
      CALL fail(ErrStat, ErrMsg, 'coupled motion updates must be finite')
      RETURN
    END IF
    model%q(model%fixed_dofs) = q_coupled
    model%v(model%fixed_dofs) = v_coupled
    model%a(model%fixed_dofs) = a_coupled
  END SUBROUTINE CD_Update_Model_CoupledMotion

  SUBROUTINE CD_Update_Model_Interior_State(model, q_interior, v_interior, ErrStat, ErrMsg, a_interior)
    !! Overwrite the model's INTERIOR-node state (DOFs 4 .. n_dof-3; both endpoint nodes
    !! untouched) -- the checkpoint-restart reload boundary. The framework state mirror
    !! carries exactly these interior [v; q; a] blocks (endpoints are boundary: coupled
    !! ends are re-prescribed by the host mesh, anchors never move), so a restart
    !! rebuilds the model from the deck and reloads the interiors here. a_interior
    !! reloads the generalised-alpha COMMITTED acceleration exactly; when omitted, call
    !! CD_Recompute_Model_Acceleration afterwards (an equations-consistent re-derivation
    !! that differs from the committed algorithmic acceleration at the solve-tolerance
    !! scale -- measured to move the first restored load evaluation, which is why the
    !! mirror carries the acceleration).
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: q_interior(:), v_interior(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: a_interior(:)

    INTEGER :: nin

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    nin = model%n_dof - 6
    IF (nin <= 0) THEN
      CALL fail(ErrStat, ErrMsg, 'interior-state update needs a line with interior nodes')
      RETURN
    END IF
    IF (SIZE(q_interior) /= nin .OR. SIZE(v_interior) /= nin) THEN
      CALL fail(ErrStat, ErrMsg, 'interior-state updates must have shape (n_dof - 6)')
      RETURN
    END IF
    IF (.NOT. (CD_All_Finite(q_interior) .AND. CD_All_Finite(v_interior))) THEN
      CALL fail(ErrStat, ErrMsg, 'interior-state updates must be finite')
      RETURN
    END IF
    IF (PRESENT(a_interior)) THEN
      IF (SIZE(a_interior) /= nin) THEN
        CALL fail(ErrStat, ErrMsg, 'interior-state updates must have shape (n_dof - 6)')
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(a_interior)) THEN
        CALL fail(ErrStat, ErrMsg, 'interior-state updates must be finite')
        RETURN
      END IF
    END IF
    model%q(4:model%n_dof - 3) = q_interior
    model%v(4:model%n_dof - 3) = v_interior
    IF (PRESENT(a_interior)) model%a(4:model%n_dof - 3) = a_interior
  END SUBROUTINE CD_Update_Model_Interior_State

  SUBROUTINE CD_Recompute_Model_Acceleration(model, ErrStat, ErrMsg)
    !! Recompute the model-owned acceleration from the current q/v state,
    !! external loads, environmental loads, added mass, and fixed/coupled DOF map.
    !! Use this after runtime updates when the caller needs a consistent
    !! acceleration before a load exchange or diagnostic query.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: es
    CHARACTER(160) :: em
    LOGICAL :: has_load
    INTEGER :: nc
    LOGICAL :: prescribed_workspace_ready
    ! Entry acceleration: every failure below restores it whole, so a failed
    ! recompute leaves the committed state untouched.
    REAL(wp) :: a_entry(model%n_dof)

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    has_load = model%has_seabed .OR. model%has_seabed_damping .OR. model%has_seabed_friction .OR. &
               model%has_damping .OR. model%has_viscoelastic .OR. model%has_syrope .OR. model%has_morison_drag .OR. &
               model%has_froude_krylov .OR. model%has_buoyancy_recovery
    nc = SIZE(model%fixed_dofs)
    prescribed_workspace_ready = ALLOCATED(model%a_prescribed_work)
    IF (prescribed_workspace_ready) prescribed_workspace_ready = SIZE(model%a_prescribed_work) >= nc
    IF (.NOT. prescribed_workspace_ready) THEN
      CALL fail(ErrStat, ErrMsg, 'acceleration recompute requires initialized prescribed-acceleration workspace')
      RETURN
    END IF
    a_entry = model%a
    model%a_prescribed_work(1:nc) = model%a(model%fixed_dofs)
    IF (has_load .AND. model%has_added_mass) THEN
      CALL CD_Cable_Initial_Acceleration(model%q, model%v, model%elem_conn, model%l0, model%ea_dyn, &
                                         model%rho_a, model%tension_only, model%f_ext, &
                                         model%fixed_dofs, model%a, es, em, &
                                         load_force_proc=model_force, workspace=model%dynamic_workspace, &
                                         added_mass_band_proc=model_added_mass_band)
    ELSE IF (has_load) THEN
      CALL CD_Cable_Initial_Acceleration(model%q, model%v, model%elem_conn, model%l0, model%ea_dyn, &
                                         model%rho_a, model%tension_only, model%f_ext, &
                                         model%fixed_dofs, model%a, es, em, &
                                         load_force_proc=model_force, workspace=model%dynamic_workspace)
    ELSE IF (model%has_added_mass) THEN
      CALL CD_Cable_Initial_Acceleration(model%q, model%v, model%elem_conn, model%l0, model%ea_dyn, &
                                         model%rho_a, model%tension_only, model%f_ext, &
                                         model%fixed_dofs, model%a, es, em, &
                                         workspace=model%dynamic_workspace, added_mass_band_proc=model_added_mass_band)
    ELSE
      CALL CD_Cable_Initial_Acceleration(model%q, model%v, model%elem_conn, model%l0, model%ea_dyn, &
                                         model%rho_a, model%tension_only, model%f_ext, &
                                         model%fixed_dofs, model%a, es, em, workspace=model%dynamic_workspace)
    END IF
    IF (es /= CD_DYN_OK) THEN
      model%a = a_entry
      ErrStat = model_status_from_dyn(es)
      ErrMsg = 'CableDyn_Model: acceleration recompute failed: '//TRIM(em)
      RETURN
    END IF
    model%a(model%fixed_dofs) = model%a_prescribed_work(1:nc)
    CALL apply_prescribed_acceleration_coupling(model, model%a_prescribed_work(1:nc), ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) THEN
      model%a = a_entry
      RETURN
    END IF

  CONTAINS

    SUBROUTINE apply_prescribed_acceleration_coupling(model, prescribed_a, ErrStat, ErrMsg)
      TYPE(CD_ModelType), INTENT(INOUT) :: model
      REAL(wp), INTENT(IN) :: prescribed_a(:)
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg

      INTEGER :: es, n_free, i, j, row, col
      CHARACTER(160) :: em

      ErrStat = CD_MODEL_OK
      ErrMsg = ''
      IF (SIZE(prescribed_a) == 0) RETURN
      IF (MAXVAL(ABS(prescribed_a)) <= EPSILON(CD_ZERO)) RETURN
      CALL partition_model_free_work(model, n_free, ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
      IF (n_free == 0) RETURN
      CALL ensure_support_acceleration_workspace(model, n_free, ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
      ASSOCIATE (M => model%dynamic_workspace%M(1:model%n_dof, 1:model%n_dof), &
                 M_free => model%dynamic_workspace%eff_free(1:n_free, 1:n_free), &
                 rhs => model%dynamic_workspace%R_free(1:n_free))
        CALL CD_Assemble_Cable_Mass(model%elem_conn, model%l0, model%rho_a, M, es, em, topology_validated=.TRUE.)
        IF (es /= 0) THEN
          ErrStat = CD_MODEL_SOLVEFAIL
          ErrMsg = 'CableDyn_Model: support-inertia mass assembly failed: '//TRIM(em)
          RETURN
        END IF
        IF (model%has_added_mass) THEN
          ASSOCIATE (M_add => model%dynamic_workspace%M_add_eval(1:model%n_dof, 1:model%n_dof), &
                     dMa_a_dq => model%dynamic_workspace%dMa_a_dq_eval(1:model%n_dof, 1:model%n_dof))
            CALL model_added_mass_full(model, model%q, model%a, M_add, dMa_a_dq, es, em)
            IF (es /= 0) THEN
              ErrStat = CD_MODEL_SOLVEFAIL
              ErrMsg = 'CableDyn_Model: support-inertia added-mass assembly failed: '//TRIM(em)
              RETURN
            END IF
            M = M + M_add
          END ASSOCIATE
        END IF
        DO j = 1, n_free
          col = model%free_map_work(j)
          DO i = 1, n_free
            row = model%free_map_work(i)
            M_free(i, j) = M(row, col)
          END DO
        END DO
        rhs = CD_ZERO
        DO j = 1, SIZE(model%fixed_dofs)
          col = model%fixed_dofs(j)
          DO i = 1, n_free
            row = model%free_map_work(i)
            rhs(i) = rhs(i) + M(row, col)*prescribed_a(j)
          END DO
        END DO
        CALL CD_Solve_Dense_As_Banded(M_free, rhs, es, em)
        IF (es /= 0 .OR. .NOT. CD_All_Finite(rhs)) THEN
          ErrStat = CD_MODEL_SOLVEFAIL
          ErrMsg = 'CableDyn_Model: support-inertia solve failed: '//TRIM(em)
          RETURN
        END IF
        DO i = 1, n_free
          row = model%free_map_work(i)
          model%a(row) = model%a(row) - rhs(i)
        END DO
      END ASSOCIATE
    END SUBROUTINE apply_prescribed_acceleration_coupling

    SUBROUTINE partition_model_free_work(model, n_free, ErrStat, ErrMsg)
      TYPE(CD_ModelType), INTENT(INOUT) :: model
      INTEGER, INTENT(OUT) :: n_free
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg

      INTEGER :: i, dof
      LOGICAL :: free_map_ready, marker_ready

      ErrStat = CD_MODEL_OK
      ErrMsg = ''
      n_free = 0
      free_map_ready = ALLOCATED(model%free_map_work)
      IF (free_map_ready) free_map_ready = SIZE(model%free_map_work) >= model%n_dof
      IF (.NOT. free_map_ready) THEN
        CALL fail(ErrStat, ErrMsg, 'support-inertia coupling requires initialized free-map workspace')
        RETURN
      END IF
      marker_ready = ALLOCATED(model%dynamic_workspace%dof_marker)
      IF (marker_ready) marker_ready = SIZE(model%dynamic_workspace%dof_marker) >= model%n_dof
      IF (.NOT. marker_ready) THEN
        CALL fail(ErrStat, ErrMsg, 'support-inertia coupling requires initialized DOF marker workspace')
        RETURN
      END IF
      model%dynamic_workspace%dof_marker(1:model%n_dof) = 0
      DO i = 1, SIZE(model%fixed_dofs)
        dof = model%fixed_dofs(i)
        IF (dof < 1 .OR. dof > model%n_dof) THEN
          CALL fail(ErrStat, ErrMsg, 'support-inertia fixed DOF is out of range')
          RETURN
        END IF
        IF (model%dynamic_workspace%dof_marker(dof) /= 0) THEN
          CALL fail(ErrStat, ErrMsg, 'support-inertia fixed DOFs contain duplicates')
          RETURN
        END IF
        model%dynamic_workspace%dof_marker(dof) = 1
      END DO
      DO dof = 1, model%n_dof
        IF (model%dynamic_workspace%dof_marker(dof) == 0) THEN
          n_free = n_free + 1
          model%free_map_work(n_free) = dof
        END IF
      END DO
    END SUBROUTINE partition_model_free_work

    SUBROUTINE ensure_support_acceleration_workspace(model, n_free, ErrStat, ErrMsg)
      TYPE(CD_ModelType), INTENT(INOUT) :: model
      INTEGER, INTENT(IN) :: n_free
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg

      INTEGER :: istat, n
      LOGICAL :: need_workspace

      ErrStat = CD_MODEL_OK
      ErrMsg = ''
      n = model%n_dof
      need_workspace = .NOT. ALLOCATED(model%dynamic_workspace%M)
      IF (.NOT. need_workspace) need_workspace = SIZE(model%dynamic_workspace%M, 1) < n .OR. &
                                                 SIZE(model%dynamic_workspace%M, 2) < n
      IF (need_workspace) THEN
        IF (ALLOCATED(model%dynamic_workspace%M)) DEALLOCATE (model%dynamic_workspace%M)
        ALLOCATE (model%dynamic_workspace%M(n, n), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'support-inertia mass workspace allocation failed')
          RETURN
        END IF
      END IF
      need_workspace = .NOT. ALLOCATED(model%dynamic_workspace%eff_free)
      IF (.NOT. need_workspace) need_workspace = SIZE(model%dynamic_workspace%eff_free, 1) < n_free .OR. &
                                                 SIZE(model%dynamic_workspace%eff_free, 2) < n_free
      IF (need_workspace) THEN
        IF (ALLOCATED(model%dynamic_workspace%eff_free)) DEALLOCATE (model%dynamic_workspace%eff_free)
        ALLOCATE (model%dynamic_workspace%eff_free(n_free, n_free), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'support-inertia reduced-mass workspace allocation failed')
          RETURN
        END IF
      END IF
      need_workspace = .NOT. ALLOCATED(model%dynamic_workspace%R_free)
      IF (.NOT. need_workspace) need_workspace = SIZE(model%dynamic_workspace%R_free) < n_free
      IF (need_workspace) THEN
        IF (ALLOCATED(model%dynamic_workspace%R_free)) DEALLOCATE (model%dynamic_workspace%R_free)
        ALLOCATE (model%dynamic_workspace%R_free(n_free), STAT=istat)
        IF (istat /= 0) THEN
          CALL alloc_fail(ErrStat, ErrMsg, 'support-inertia rhs workspace allocation failed')
          RETURN
        END IF
      END IF
      IF (model%has_added_mass) THEN
        need_workspace = .NOT. ALLOCATED(model%dynamic_workspace%M_add_eval)
        IF (.NOT. need_workspace) need_workspace = SIZE(model%dynamic_workspace%M_add_eval, 1) < n .OR. &
                                                   SIZE(model%dynamic_workspace%M_add_eval, 2) < n
        IF (need_workspace) THEN
          IF (ALLOCATED(model%dynamic_workspace%M_add_eval)) DEALLOCATE (model%dynamic_workspace%M_add_eval)
          ALLOCATE (model%dynamic_workspace%M_add_eval(n, n), STAT=istat)
          IF (istat /= 0) THEN
            CALL alloc_fail(ErrStat, ErrMsg, 'support-inertia added-mass workspace allocation failed')
            RETURN
          END IF
        END IF
        need_workspace = .NOT. ALLOCATED(model%dynamic_workspace%dMa_a_dq_eval)
        IF (.NOT. need_workspace) need_workspace = SIZE(model%dynamic_workspace%dMa_a_dq_eval, 1) < n .OR. &
                                                   SIZE(model%dynamic_workspace%dMa_a_dq_eval, 2) < n
        IF (need_workspace) THEN
          IF (ALLOCATED(model%dynamic_workspace%dMa_a_dq_eval)) DEALLOCATE (model%dynamic_workspace%dMa_a_dq_eval)
          ALLOCATE (model%dynamic_workspace%dMa_a_dq_eval(n, n), STAT=istat)
          IF (istat /= 0) THEN
            CALL alloc_fail(ErrStat, ErrMsg, 'support-inertia added-mass tangent workspace allocation failed')
            RETURN
          END IF
        END IF
      END IF
    END SUBROUTINE ensure_support_acceleration_workspace

    SUBROUTINE model_force(q, v, force, ErrStat, ErrMsg)
      REAL(wp), INTENT(IN) :: q(:), v(:)
      REAL(wp), INTENT(OUT) :: force(:)
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg
      CALL model_contributor_force(model, q, v, force, ErrStat, ErrMsg)
      ErrStat = as_dyn_callback_status(ErrStat)
    END SUBROUTINE model_force

    SUBROUTINE model_added_mass_band(q, accel, free, kl, ku, M_add_a, M_add_band, dMa_a_dq_band, ErrStat, ErrMsg, &
                                     need_tangent)
      REAL(wp), INTENT(IN) :: q(:), accel(:)
      INTEGER, INTENT(IN) :: free(:), kl, ku
      REAL(wp), INTENT(OUT) :: M_add_a(:), M_add_band(:, :), dMa_a_dq_band(:, :)
      INTEGER, INTENT(OUT) :: ErrStat
      CHARACTER(*), INTENT(OUT) :: ErrMsg
      LOGICAL, INTENT(IN), OPTIONAL :: need_tangent
      CALL model_added_mass_banded(model, q, accel, free, kl, ku, M_add_a, M_add_band, dMa_a_dq_band, &
                                   ErrStat, ErrMsg, need_tangent=need_tangent)
      ErrStat = as_dyn_callback_status(ErrStat)
    END SUBROUTINE model_added_mass_band
  END SUBROUTINE CD_Recompute_Model_Acceleration

  SUBROUTINE CD_End_Model(model, ErrStat, ErrMsg)
    !! Release all allocatable model storage. Safe to call on an uninitialised model.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    IF (ALLOCATED(model%elem_conn)) DEALLOCATE (model%elem_conn)
    IF (ALLOCATED(model%fixed_dofs)) DEALLOCATE (model%fixed_dofs)
    IF (ALLOCATED(model%l0)) DEALLOCATE (model%l0)
    IF (ALLOCATED(model%ea)) DEALLOCATE (model%ea)
    IF (ALLOCATED(model%rho_a)) DEALLOCATE (model%rho_a)
    IF (ALLOCATED(model%f_ext)) DEALLOCATE (model%f_ext)
    IF (ALLOCATED(model%l0_dot)) DEALLOCATE (model%l0_dot)
    IF (ALLOCATED(model%dist_load)) DEALLOCATE (model%dist_load)
    IF (ALLOCATED(model%f_ext_dist_base)) DEALLOCATE (model%f_ext_dist_base)
    IF (ALLOCATED(model%seabed_diameter_elem)) DEALLOCATE (model%seabed_diameter_elem)
    ! reset the recompute state with its array: CD_Init_Model only ever SETS the
    ! flag, so a reused model would otherwise inherit .TRUE. without the inputs
    ! and the length update would read a deallocated diameter array
    model%has_seabed_recompute = .FALSE.
    model%seabed_kbot = CD_ZERO
    model%seabed_cn_over_kn = CD_ZERO
    IF (ALLOCATED(model%q)) DEALLOCATE (model%q)
    IF (ALLOCATED(model%v)) DEALLOCATE (model%v)
    IF (ALLOCATED(model%a)) DEALLOCATE (model%a)
    IF (ALLOCATED(model%q_step)) DEALLOCATE (model%q_step)
    IF (ALLOCATED(model%v_step)) DEALLOCATE (model%v_step)
    IF (ALLOCATED(model%a_step)) DEALLOCATE (model%a_step)
    IF (ALLOCATED(model%q_prescribed_work)) DEALLOCATE (model%q_prescribed_work)
    IF (ALLOCATED(model%v_prescribed_work)) DEALLOCATE (model%v_prescribed_work)
    IF (ALLOCATED(model%a_prescribed_work)) DEALLOCATE (model%a_prescribed_work)
    IF (ALLOCATED(model%q_rollback)) DEALLOCATE (model%q_rollback)
    IF (ALLOCATED(model%v_rollback)) DEALLOCATE (model%v_rollback)
    IF (ALLOCATED(model%a_rollback)) DEALLOCATE (model%a_rollback)
    IF (ALLOCATED(model%recovery_q0)) DEALLOCATE (model%recovery_q0)
    IF (ALLOCATED(model%recovery_v0)) DEALLOCATE (model%recovery_v0)
    IF (ALLOCATED(model%recovery_a0)) DEALLOCATE (model%recovery_a0)
    IF (ALLOCATED(model%f_ext_rollback)) DEALLOCATE (model%f_ext_rollback)
    IF (ALLOCATED(model%free_map_work)) DEALLOCATE (model%free_map_work)
    CALL CD_Clear_GenAlpha_Workspace(model%dynamic_workspace)
    CALL clear_model_load_contributor_workspace(model)
    IF (ALLOCATED(model%seabed_kn)) DEALLOCATE (model%seabed_kn)
    IF (ALLOCATED(model%seabed_cn)) DEALLOCATE (model%seabed_cn)
    CALL CD_End_Bathymetry(model%bathymetry)
    IF (ALLOCATED(model%ba)) DEALLOCATE (model%ba)
    IF (ALLOCATED(model%ba_dyn)) DEALLOCATE (model%ba_dyn)
    IF (ALLOCATED(model%ea_dyn)) DEALLOCATE (model%ea_dyn)
    IF (ALLOCATED(model%ve_mode)) DEALLOCATE (model%ve_mode)
    IF (ALLOCATED(model%ve_alpha_mbl)) DEALLOCATE (model%ve_alpha_mbl)
    IF (ALLOCATED(model%ve_vbeta)) DEALLOCATE (model%ve_vbeta)
    IF (ALLOCATED(model%ve_ea_d)) DEALLOCATE (model%ve_ea_d)
    IF (ALLOCATED(model%ve_ea_1)) DEALLOCATE (model%ve_ea_1)
    IF (ALLOCATED(model%ve_ba)) DEALLOCATE (model%ve_ba)
    IF (ALLOCATED(model%ve_ba_d)) DEALLOCATE (model%ve_ba_d)
    IF (ALLOCATED(model%ve_dl_1)) DEALLOCATE (model%ve_dl_1)
    IF (ALLOCATED(model%ve_dl_1_rollback)) DEALLOCATE (model%ve_dl_1_rollback)
    IF (ALLOCATED(model%recovery_ve_dl_1)) DEALLOCATE (model%recovery_ve_dl_1)
    IF (ALLOCATED(model%syrope_type)) THEN
      BLOCK
        INTEGER :: e
        DO e = 1, SIZE(model%syrope_type)
          CALL CD_Syrope_End(model%syrope_type(e))
        END DO
      END BLOCK
      DEALLOCATE (model%syrope_type)
    END IF
    IF (ALLOCATED(model%syrope_is)) DEALLOCATE (model%syrope_is)
    IF (ALLOCATED(model%syrope_slow)) DEALLOCATE (model%syrope_slow)
    IF (ALLOCATED(model%syrope_slow_rollback)) DEALLOCATE (model%syrope_slow_rollback)
    IF (ALLOCATED(model%syrope_tmax)) DEALLOCATE (model%syrope_tmax)
    IF (ALLOCATED(model%syrope_tmax_rollback)) DEALLOCATE (model%syrope_tmax_rollback)
    IF (ALLOCATED(model%recovery_syrope_slow)) DEALLOCATE (model%recovery_syrope_slow)
    IF (ALLOCATED(model%recovery_syrope_tmax)) DEALLOCATE (model%recovery_syrope_tmax)
    IF (ALLOCATED(model%fr_anchor)) DEALLOCATE (model%fr_anchor)
    IF (ALLOCATED(model%fr_anchor_rollback)) DEALLOCATE (model%fr_anchor_rollback)
    IF (ALLOCATED(model%recovery_fr_anchor)) DEALLOCATE (model%recovery_fr_anchor)
    IF (ALLOCATED(model%fluid_velocity)) DEALLOCATE (model%fluid_velocity)
    IF (ALLOCATED(model%drag_waterline_z)) DEALLOCATE (model%drag_waterline_z)
    IF (ALLOCATED(model%drag_diameter_elem)) DEALLOCATE (model%drag_diameter_elem)
    IF (ALLOCATED(model%drag_cdn_elem)) DEALLOCATE (model%drag_cdn_elem)
    IF (ALLOCATED(model%drag_cdt_elem)) DEALLOCATE (model%drag_cdt_elem)
    IF (ALLOCATED(model%fluid_acceleration)) DEALLOCATE (model%fluid_acceleration)
    IF (ALLOCATED(model%fk_waterline_z)) DEALLOCATE (model%fk_waterline_z)
    IF (ALLOCATED(model%fk_diameter_elem)) DEALLOCATE (model%fk_diameter_elem)
    IF (ALLOCATED(model%fk_can_elem)) DEALLOCATE (model%fk_can_elem)
    IF (ALLOCATED(model%fk_cat_elem)) DEALLOCATE (model%fk_cat_elem)
    IF (ALLOCATED(model%buoyancy_waterline_z)) DEALLOCATE (model%buoyancy_waterline_z)
    IF (ALLOCATED(model%buoyancy_diameter_elem)) DEALLOCATE (model%buoyancy_diameter_elem)
    IF (ALLOCATED(model%added_mass_waterline_z)) DEALLOCATE (model%added_mass_waterline_z)
    IF (ALLOCATED(model%added_mass_diameter_elem)) DEALLOCATE (model%added_mass_diameter_elem)
    IF (ALLOCATED(model%added_mass_can_elem)) DEALLOCATE (model%added_mass_can_elem)
    IF (ALLOCATED(model%added_mass_cat_elem)) DEALLOCATE (model%added_mass_cat_elem)
    model%initialized = .FALSE.
    model%tension_only = .FALSE.
    model%has_seabed = .FALSE.
    model%has_bathymetry = .FALSE.
    model%has_seabed_damping = .FALSE.
    model%has_seabed_friction = .FALSE.
    model%has_damping = .FALSE.
    model%has_viscoelastic = .FALSE.
    model%ve_dt = CD_ZERO
    model%has_syrope = .FALSE.
    model%has_morison_drag = .FALSE.
    model%has_froude_krylov = .FALSE.
    model%has_buoyancy_recovery = .FALSE.
    model%has_added_mass = .FALSE.
    model%seabed_z_floor = CD_ZERO
    model%seabed_mu = CD_ZERO
    model%seabed_friction_aniso = .FALSE.
    model%seabed_mu_axial = CD_ZERO
    model%drag_rho = CD_ZERO
    model%drag_diameter = CD_ZERO
    model%drag_cdn = CD_ZERO
    model%drag_cdt = CD_ZERO
    model%fk_rho = CD_ZERO
    model%fk_diameter = CD_ZERO
    model%fk_can = CD_ZERO
    model%fk_cat = CD_ZERO
    model%buoyancy_rho = CD_ZERO
    model%buoyancy_diameter = CD_ZERO
    model%buoyancy_gravity = CD_ZERO
    model%added_mass_rho = CD_ZERO
    model%added_mass_diameter = CD_ZERO
    model%added_mass_can = CD_ZERO
    model%added_mass_cat = CD_ZERO
    model%n_dof = 0
    model%n_elem = 0
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END SUBROUTINE CD_End_Model

  SUBROUTINE CD_Copy_Model(src, dst, ErrStat, ErrMsg)
    !! Deep-copy persistent model state with checked allocation. Transient solver
    !! workspaces are intentionally not copied; they are rebuilt on demand.
    TYPE(CD_ModelType), INTENT(IN) :: src
    TYPE(CD_ModelType), INTENT(INOUT) :: dst
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    TYPE(CD_ModelType) :: work
    INTEGER :: es
    CHARACTER(200) :: em

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    ! Build the copy into a temporary model and commit it into dst only after every checked
    ! allocation succeeds. dst is left untouched on any failure, so a caller's already-initialized
    ! model survives an out-of-memory copy intact. work is a local: its allocatable components are
    ! released automatically on every early return. The final assignment is the single commit point.

    work%tension_only = src%tension_only
    work%n_dof = src%n_dof
    work%n_elem = src%n_elem
    work%has_seabed = src%has_seabed
    work%has_bathymetry = src%has_bathymetry
    work%seabed_z_floor = src%seabed_z_floor
    work%has_seabed_damping = src%has_seabed_damping
    work%has_seabed_friction = src%has_seabed_friction
    work%seabed_mu = src%seabed_mu
    work%seabed_friction_aniso = src%seabed_friction_aniso
    work%seabed_mu_axial = src%seabed_mu_axial
    work%has_seabed_recompute = src%has_seabed_recompute
    work%seabed_kbot = src%seabed_kbot
    work%seabed_cn_over_kn = src%seabed_cn_over_kn
    work%has_damping = src%has_damping
    work%has_viscoelastic = src%has_viscoelastic
    work%ve_dt = src%ve_dt
    work%has_syrope = src%has_syrope
    work%has_morison_drag = src%has_morison_drag
    work%drag_rho = src%drag_rho
    work%drag_diameter = src%drag_diameter
    work%drag_cdn = src%drag_cdn
    work%drag_cdt = src%drag_cdt
    work%has_froude_krylov = src%has_froude_krylov
    work%fk_rho = src%fk_rho
    work%fk_diameter = src%fk_diameter
    work%fk_can = src%fk_can
    work%fk_cat = src%fk_cat
    work%has_buoyancy_recovery = src%has_buoyancy_recovery
    work%buoyancy_rho = src%buoyancy_rho
    work%buoyancy_diameter = src%buoyancy_diameter
    work%buoyancy_gravity = src%buoyancy_gravity
    work%has_added_mass = src%has_added_mass
    work%added_mass_rho = src%added_mass_rho
    work%added_mass_diameter = src%added_mass_diameter
    work%added_mass_can = src%added_mass_can
    work%added_mass_cat = src%added_mass_cat
    work%cfg = src%cfg

    CALL copy_int2(src%elem_conn, work%elem_conn, 'elem_conn', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_int1(src%fixed_dofs, work%fixed_dofs, 'fixed_dofs', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%l0, work%l0, 'l0', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ea, work%ea, 'ea', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%rho_a, work%rho_a, 'rho_a', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%l0_dot, work%l0_dot, 'l0_dot', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real2(src%dist_load, work%dist_load, 'dist_load', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%f_ext_dist_base, work%f_ext_dist_base, 'f_ext_dist_base', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%seabed_diameter_elem, work%seabed_diameter_elem, 'seabed_diameter_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%f_ext, work%f_ext, 'f_ext', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%q, work%q, 'q', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%v, work%v, 'v', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%a, work%a, 'a', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%q_step, work%q_step, 'q_step', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%v_step, work%v_step, 'v_step', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%a_step, work%a_step, 'a_step', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%q_prescribed_work, work%q_prescribed_work, 'q_prescribed_work', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%v_prescribed_work, work%v_prescribed_work, 'v_prescribed_work', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%a_prescribed_work, work%a_prescribed_work, 'a_prescribed_work', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%q_rollback, work%q_rollback, 'q_rollback', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%v_rollback, work%v_rollback, 'v_rollback', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%a_rollback, work%a_rollback, 'a_rollback', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%recovery_q0, work%recovery_q0, 'recovery_q0', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%recovery_v0, work%recovery_v0, 'recovery_v0', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%recovery_a0, work%recovery_a0, 'recovery_a0', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%f_ext_rollback, work%f_ext_rollback, 'f_ext_rollback', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_int1(src%free_map_work, work%free_map_work, 'free_map_work', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN

    CALL copy_real1(src%seabed_kn, work%seabed_kn, 'seabed_kn', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%seabed_cn, work%seabed_cn, 'seabed_cn', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    IF (src%has_bathymetry) THEN
      IF (.NOT. CD_Bathymetry_Is_Initialized(src%bathymetry)) THEN
        ErrStat = CD_MODEL_BADINPUT
        ErrMsg = 'CableDyn_Model: source model has uninitialized bathymetry'
        RETURN
      END IF
      CALL CD_Init_Bathymetry(work%bathymetry, src%bathymetry%x, src%bathymetry%y, src%bathymetry%depth, es, em)
      IF (es /= CD_BATHY_OK) THEN
        ErrStat = CD_MODEL_ALLOCFAIL
        ErrMsg = 'CableDyn_Model: bathymetry copy failed: '//TRIM(em)
        RETURN
      END IF
    END IF
    CALL copy_real1(src%ba, work%ba, 'ba', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ba_dyn, work%ba_dyn, 'ba_dyn', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ea_dyn, work%ea_dyn, 'ea_dyn', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_int1(src%ve_mode, work%ve_mode, 've_mode', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ve_alpha_mbl, work%ve_alpha_mbl, 've_alpha_mbl', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ve_vbeta, work%ve_vbeta, 've_vbeta', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ve_ea_d, work%ve_ea_d, 've_ea_d', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ve_ea_1, work%ve_ea_1, 've_ea_1', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ve_ba, work%ve_ba, 've_ba', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ve_ba_d, work%ve_ba_d, 've_ba_d', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ve_dl_1, work%ve_dl_1, 've_dl_1', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%ve_dl_1_rollback, work%ve_dl_1_rollback, 've_dl_1_rollback', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%recovery_ve_dl_1, work%recovery_ve_dl_1, 'recovery_ve_dl_1', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_syrope_types(src%syrope_type, work%syrope_type, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_logical1(src%syrope_is, work%syrope_is, 'syrope_is', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%syrope_slow, work%syrope_slow, 'syrope_slow', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%syrope_slow_rollback, work%syrope_slow_rollback, 'syrope_slow_rollback', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%syrope_tmax, work%syrope_tmax, 'syrope_tmax', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%syrope_tmax_rollback, work%syrope_tmax_rollback, 'syrope_tmax_rollback', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%recovery_syrope_slow, work%recovery_syrope_slow, 'recovery_syrope_slow', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%recovery_syrope_tmax, work%recovery_syrope_tmax, 'recovery_syrope_tmax', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real2(src%fr_anchor, work%fr_anchor, 'fr_anchor', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real2(src%fr_anchor_rollback, work%fr_anchor_rollback, 'fr_anchor_rollback', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real2(src%recovery_fr_anchor, work%recovery_fr_anchor, 'recovery_fr_anchor', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real2(src%fluid_velocity, work%fluid_velocity, 'fluid_velocity', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%drag_waterline_z, work%drag_waterline_z, 'drag_waterline_z', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%drag_diameter_elem, work%drag_diameter_elem, 'drag_diameter_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%drag_cdn_elem, work%drag_cdn_elem, 'drag_cdn_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%drag_cdt_elem, work%drag_cdt_elem, 'drag_cdt_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real2(src%fluid_acceleration, work%fluid_acceleration, 'fluid_acceleration', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%fk_waterline_z, work%fk_waterline_z, 'fk_waterline_z', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%fk_diameter_elem, work%fk_diameter_elem, 'fk_diameter_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%fk_can_elem, work%fk_can_elem, 'fk_can_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%fk_cat_elem, work%fk_cat_elem, 'fk_cat_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%buoyancy_waterline_z, work%buoyancy_waterline_z, 'buoyancy_waterline_z', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%buoyancy_diameter_elem, work%buoyancy_diameter_elem, 'buoyancy_diameter_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%added_mass_waterline_z, work%added_mass_waterline_z, 'added_mass_waterline_z', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%added_mass_diameter_elem, work%added_mass_diameter_elem, 'added_mass_diameter_elem', ErrStat, &
                    ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%added_mass_can_elem, work%added_mass_can_elem, 'added_mass_can_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL copy_real1(src%added_mass_cat_elem, work%added_mass_cat_elem, 'added_mass_cat_elem', ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN

    work%initialized = src%initialized
    ! Commit: release the caller's previous model, then move the fully-built copy in component by
    ! component. MOVE_ALLOC cannot fail, so the commit itself cannot abort under memory pressure --
    ! an intrinsic dst = work would re-allocate every allocatable component again without STAT.
    CALL CD_End_Model(dst, es, em)
    CALL commit_model_copy(work, dst)
  END SUBROUTINE CD_Copy_Model

  SUBROUTINE commit_model_copy(from, to)
    !! Move every persistent component of a freshly-built model copy into the just-ended destination
    !! without reallocating: scalars are assigned and allocatable arrays (plus the bathymetry grid)
    !! are transferred with MOVE_ALLOC. The transient solver scratch (dynamic_workspace, load_*_work)
    !! is intentionally not carried, matching CD_Copy_Model. The component list mirrors the type
    !! definition and CD_End_Model; every state- or load-bearing field must be moved here, and the
    !! model-copy regression test checks that none is dropped.
    TYPE(CD_ModelType), INTENT(INOUT) :: from
    TYPE(CD_ModelType), INTENT(INOUT) :: to

    to%initialized = from%initialized
    to%tension_only = from%tension_only
    to%n_dof = from%n_dof
    to%n_elem = from%n_elem
    to%has_seabed = from%has_seabed
    to%has_bathymetry = from%has_bathymetry
    to%seabed_z_floor = from%seabed_z_floor
    to%has_seabed_damping = from%has_seabed_damping
    to%has_seabed_friction = from%has_seabed_friction
    to%seabed_mu = from%seabed_mu
    to%seabed_friction_aniso = from%seabed_friction_aniso
    to%seabed_mu_axial = from%seabed_mu_axial
    to%has_seabed_recompute = from%has_seabed_recompute
    to%seabed_kbot = from%seabed_kbot
    to%seabed_cn_over_kn = from%seabed_cn_over_kn
    CALL MOVE_ALLOC(from%seabed_diameter_elem, to%seabed_diameter_elem)
    to%has_damping = from%has_damping
    to%has_viscoelastic = from%has_viscoelastic
    to%ve_dt = from%ve_dt
    to%has_syrope = from%has_syrope
    to%has_morison_drag = from%has_morison_drag
    to%drag_rho = from%drag_rho
    to%drag_diameter = from%drag_diameter
    to%drag_cdn = from%drag_cdn
    to%drag_cdt = from%drag_cdt
    to%has_froude_krylov = from%has_froude_krylov
    to%fk_rho = from%fk_rho
    to%fk_diameter = from%fk_diameter
    to%fk_can = from%fk_can
    to%fk_cat = from%fk_cat
    to%has_buoyancy_recovery = from%has_buoyancy_recovery
    to%buoyancy_rho = from%buoyancy_rho
    to%buoyancy_diameter = from%buoyancy_diameter
    to%buoyancy_gravity = from%buoyancy_gravity
    to%has_added_mass = from%has_added_mass
    to%added_mass_rho = from%added_mass_rho
    to%added_mass_diameter = from%added_mass_diameter
    to%added_mass_can = from%added_mass_can
    to%added_mass_cat = from%added_mass_cat
    to%cfg = from%cfg

    CALL MOVE_ALLOC(from%elem_conn, to%elem_conn)
    CALL MOVE_ALLOC(from%fixed_dofs, to%fixed_dofs)
    CALL MOVE_ALLOC(from%l0, to%l0)
    CALL MOVE_ALLOC(from%ea, to%ea)
    CALL MOVE_ALLOC(from%rho_a, to%rho_a)
    CALL MOVE_ALLOC(from%l0_dot, to%l0_dot)
    CALL MOVE_ALLOC(from%dist_load, to%dist_load)
    CALL MOVE_ALLOC(from%f_ext_dist_base, to%f_ext_dist_base)
    CALL MOVE_ALLOC(from%f_ext, to%f_ext)
    CALL MOVE_ALLOC(from%q, to%q)
    CALL MOVE_ALLOC(from%v, to%v)
    CALL MOVE_ALLOC(from%a, to%a)
    CALL MOVE_ALLOC(from%q_step, to%q_step)
    CALL MOVE_ALLOC(from%v_step, to%v_step)
    CALL MOVE_ALLOC(from%a_step, to%a_step)
    CALL MOVE_ALLOC(from%q_prescribed_work, to%q_prescribed_work)
    CALL MOVE_ALLOC(from%v_prescribed_work, to%v_prescribed_work)
    CALL MOVE_ALLOC(from%a_prescribed_work, to%a_prescribed_work)
    CALL MOVE_ALLOC(from%q_rollback, to%q_rollback)
    CALL MOVE_ALLOC(from%v_rollback, to%v_rollback)
    CALL MOVE_ALLOC(from%a_rollback, to%a_rollback)
    CALL MOVE_ALLOC(from%recovery_q0, to%recovery_q0)
    CALL MOVE_ALLOC(from%recovery_v0, to%recovery_v0)
    CALL MOVE_ALLOC(from%recovery_a0, to%recovery_a0)
    CALL MOVE_ALLOC(from%f_ext_rollback, to%f_ext_rollback)
    CALL MOVE_ALLOC(from%free_map_work, to%free_map_work)
    CALL MOVE_ALLOC(from%seabed_kn, to%seabed_kn)
    CALL MOVE_ALLOC(from%seabed_cn, to%seabed_cn)
    CALL MOVE_ALLOC(from%ba, to%ba)
    CALL MOVE_ALLOC(from%ba_dyn, to%ba_dyn)
    CALL MOVE_ALLOC(from%ea_dyn, to%ea_dyn)
    CALL MOVE_ALLOC(from%ve_mode, to%ve_mode)
    CALL MOVE_ALLOC(from%ve_alpha_mbl, to%ve_alpha_mbl)
    CALL MOVE_ALLOC(from%ve_vbeta, to%ve_vbeta)
    CALL MOVE_ALLOC(from%ve_ea_d, to%ve_ea_d)
    CALL MOVE_ALLOC(from%ve_ea_1, to%ve_ea_1)
    CALL MOVE_ALLOC(from%ve_ba, to%ve_ba)
    CALL MOVE_ALLOC(from%ve_ba_d, to%ve_ba_d)
    CALL MOVE_ALLOC(from%ve_dl_1, to%ve_dl_1)
    CALL MOVE_ALLOC(from%ve_dl_1_rollback, to%ve_dl_1_rollback)
    CALL MOVE_ALLOC(from%recovery_ve_dl_1, to%recovery_ve_dl_1)
    CALL MOVE_ALLOC(from%syrope_is, to%syrope_is)
    CALL MOVE_ALLOC(from%syrope_type, to%syrope_type)
    CALL MOVE_ALLOC(from%syrope_slow, to%syrope_slow)
    CALL MOVE_ALLOC(from%syrope_slow_rollback, to%syrope_slow_rollback)
    CALL MOVE_ALLOC(from%syrope_tmax, to%syrope_tmax)
    CALL MOVE_ALLOC(from%syrope_tmax_rollback, to%syrope_tmax_rollback)
    CALL MOVE_ALLOC(from%recovery_syrope_slow, to%recovery_syrope_slow)
    CALL MOVE_ALLOC(from%recovery_syrope_tmax, to%recovery_syrope_tmax)
    CALL MOVE_ALLOC(from%fr_anchor, to%fr_anchor)
    CALL MOVE_ALLOC(from%fr_anchor_rollback, to%fr_anchor_rollback)
    CALL MOVE_ALLOC(from%recovery_fr_anchor, to%recovery_fr_anchor)
    CALL MOVE_ALLOC(from%fluid_velocity, to%fluid_velocity)
    CALL MOVE_ALLOC(from%drag_waterline_z, to%drag_waterline_z)
    CALL MOVE_ALLOC(from%drag_diameter_elem, to%drag_diameter_elem)
    CALL MOVE_ALLOC(from%drag_cdn_elem, to%drag_cdn_elem)
    CALL MOVE_ALLOC(from%drag_cdt_elem, to%drag_cdt_elem)
    CALL MOVE_ALLOC(from%fluid_acceleration, to%fluid_acceleration)
    CALL MOVE_ALLOC(from%fk_waterline_z, to%fk_waterline_z)
    CALL MOVE_ALLOC(from%fk_diameter_elem, to%fk_diameter_elem)
    CALL MOVE_ALLOC(from%fk_can_elem, to%fk_can_elem)
    CALL MOVE_ALLOC(from%fk_cat_elem, to%fk_cat_elem)
    CALL MOVE_ALLOC(from%buoyancy_waterline_z, to%buoyancy_waterline_z)
    CALL MOVE_ALLOC(from%buoyancy_diameter_elem, to%buoyancy_diameter_elem)
    CALL MOVE_ALLOC(from%added_mass_waterline_z, to%added_mass_waterline_z)
    CALL MOVE_ALLOC(from%added_mass_diameter_elem, to%added_mass_diameter_elem)
    CALL MOVE_ALLOC(from%added_mass_can_elem, to%added_mass_can_elem)
    CALL MOVE_ALLOC(from%added_mass_cat_elem, to%added_mass_cat_elem)

    IF (from%has_bathymetry) THEN
      to%bathymetry%nx = from%bathymetry%nx
      to%bathymetry%ny = from%bathymetry%ny
      to%bathymetry%average_depth = from%bathymetry%average_depth
      to%bathymetry%minimum_depth = from%bathymetry%minimum_depth
      CALL MOVE_ALLOC(from%bathymetry%x, to%bathymetry%x)
      CALL MOVE_ALLOC(from%bathymetry%y, to%bathymetry%y)
      CALL MOVE_ALLOC(from%bathymetry%depth, to%bathymetry%depth)
    END IF
  END SUBROUTINE commit_model_copy

  SUBROUTINE copy_int1(src, dst, label, ErrStat, ErrMsg)
    INTEGER, ALLOCATABLE, INTENT(IN) :: src(:)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: dst(:)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: istat

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. ALLOCATED(src)) RETURN
    ALLOCATE (dst(SIZE(src)), SOURCE=src, STAT=istat)
    IF (istat /= 0) CALL copy_alloc_fail(label, ErrStat, ErrMsg)
  END SUBROUTINE copy_int1

  SUBROUTINE copy_logical1(src, dst, label, ErrStat, ErrMsg)
    LOGICAL, ALLOCATABLE, INTENT(IN) :: src(:)
    LOGICAL, ALLOCATABLE, INTENT(INOUT) :: dst(:)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: istat

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. ALLOCATED(src)) RETURN
    ALLOCATE (dst(SIZE(src)), SOURCE=src, STAT=istat)
    IF (istat /= 0) CALL copy_alloc_fail(label, ErrStat, ErrMsg)
  END SUBROUTINE copy_logical1

  SUBROUTINE copy_syrope_types(src, dst, ErrStat, ErrMsg)
    !! Deep-copy the per-element Syrope constitutive array. Each element's
    !! allocatable components (the OWC table + its slow-strain column) are
    !! carried by intrinsic derived-type assignment; the outer allocation is
    !! STAT-checked so an out-of-memory copy leaves the destination model whole.
    TYPE(CD_SyropeType), ALLOCATABLE, INTENT(IN) :: src(:)
    TYPE(CD_SyropeType), ALLOCATABLE, INTENT(INOUT) :: dst(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: istat, e, es2
    CHARACTER(160) :: em2

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. ALLOCATED(src)) RETURN
    ALLOCATE (dst(SIZE(src)), STAT=istat)
    IF (istat /= 0) THEN
      CALL copy_alloc_fail('syrope_type', ErrStat, ErrMsg)
      RETURN
    END IF
    ! Re-Init each populated element from the source's own raw inputs: this is
    ! a STAT-checked deep copy (CD_Syrope_Init allocates with STAT and fails
    ! closed), preserving the "dst survives an out-of-memory copy" contract that
    ! an unchecked intrinsic `dst(e) = src(e)` would break. Plain (uninitialized)
    ! entries stay default -- they carry no Syrope element.
    DO e = 1, SIZE(src)
      IF (.NOT. CD_Syrope_Is_Ready(src(e))) CYCLE
      CALL CD_Syrope_Init(dst(e), src(e)%owc_strain, src(e)%owc_tension, src(e)%wc_mod, &
                          src(e)%p1, src(e)%p2, src(e)%alpha, src(e)%beta, src(e)%c1, src(e)%c2, es2, em2)
      IF (es2 /= CD_SYROPE_OK) THEN
        CALL copy_alloc_fail('syrope_type', ErrStat, ErrMsg)
        RETURN
      END IF
    END DO
  END SUBROUTINE copy_syrope_types

  SUBROUTINE copy_int2(src, dst, label, ErrStat, ErrMsg)
    INTEGER, ALLOCATABLE, INTENT(IN) :: src(:, :)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: dst(:, :)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: istat

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. ALLOCATED(src)) RETURN
    ALLOCATE (dst(SIZE(src, 1), SIZE(src, 2)), SOURCE=src, STAT=istat)
    IF (istat /= 0) CALL copy_alloc_fail(label, ErrStat, ErrMsg)
  END SUBROUTINE copy_int2

  SUBROUTINE copy_real1(src, dst, label, ErrStat, ErrMsg)
    REAL(wp), ALLOCATABLE, INTENT(IN) :: src(:)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: dst(:)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: istat

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. ALLOCATED(src)) RETURN
    ALLOCATE (dst(SIZE(src)), SOURCE=src, STAT=istat)
    IF (istat /= 0) CALL copy_alloc_fail(label, ErrStat, ErrMsg)
  END SUBROUTINE copy_real1

  SUBROUTINE copy_real2(src, dst, label, ErrStat, ErrMsg)
    REAL(wp), ALLOCATABLE, INTENT(IN) :: src(:, :)
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: dst(:, :)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: istat

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. ALLOCATED(src)) RETURN
    ALLOCATE (dst(SIZE(src, 1), SIZE(src, 2)), SOURCE=src, STAT=istat)
    IF (istat /= 0) CALL copy_alloc_fail(label, ErrStat, ErrMsg)
  END SUBROUTINE copy_real2

  SUBROUTINE copy_alloc_fail(label, ErrStat, ErrMsg)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    ErrStat = CD_MODEL_ALLOCFAIL
    ErrMsg = 'CableDyn_Model: allocation failed while copying '//TRIM(label)
  END SUBROUTINE copy_alloc_fail

  SUBROUTINE restore_after_init_alloc_fail(model, old_model, label, ErrStat, ErrMsg)
    !! Restore a previously valid model if CD_Init_Model cannot allocate the new state.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    TYPE(CD_ModelType), INTENT(IN) :: old_model
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: es
    CHARACTER(160) :: em

    CALL CD_End_Model(model, es, em)
    model = old_model
    ErrStat = CD_MODEL_ALLOCFAIL
    ErrMsg = 'CableDyn_Model: allocation failed while initializing '//TRIM(label)
  END SUBROUTINE restore_after_init_alloc_fail

  INTEGER FUNCTION CD_Model_NCoupledDOF(model, ErrStat, ErrMsg) RESULT(n)
    !! Number of externally prescribed/coupled DOFs currently represented by the
    !! fixed_dofs list. This is the scalar queried by the C/OpenFAST shells before
    !! exchanging kinematics and loads.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    IF (.NOT. model%initialized) THEN
      n = 0
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    n = SIZE(model%fixed_dofs)
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END FUNCTION CD_Model_NCoupledDOF

  INTEGER FUNCTION CD_Model_NDOF(model, ErrStat, ErrMsg) RESULT(n)
    !! Total positions-only state size, 3*n_nodes. Shells use this before
    !! allocating q/v/a arrays for CD_Get_Model_State.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    IF (.NOT. model%initialized) THEN
      n = 0
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    n = model%n_dof
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END FUNCTION CD_Model_NDOF

  INTEGER FUNCTION CD_Model_NElem(model, ErrStat, ErrMsg) RESULT(n)
    !! Number of line elements owned by the model. Useful for shell-side diagnostics
    !! and output allocation.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    IF (.NOT. model%initialized) THEN
      n = 0
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    n = model%n_elem
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END FUNCTION CD_Model_NElem

  SUBROUTINE CD_Get_Model_CoupledDofs(model, coupled_dofs, ErrStat, ErrMsg)
    !! Copy the 1-based externally prescribed/coupled DOF map. The EI=0 model
    !! represents coupled endpoints through the fixed_dofs list; shells query this
    !! accessor rather than CD_ModelType internals.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(OUT) :: coupled_dofs(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    coupled_dofs = 0
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    IF (SIZE(coupled_dofs) /= SIZE(model%fixed_dofs)) THEN
      CALL fail(ErrStat, ErrMsg, 'coupled_dofs output must have shape (n_coupled_dof)')
      RETURN
    END IF
    coupled_dofs = model%fixed_dofs
  END SUBROUTINE CD_Get_Model_CoupledDofs

  SUBROUTINE CD_Get_Model_CoupledMotion(model, q_coupled, v_coupled, a_coupled, ErrStat, ErrMsg)
    !! Copy the current coupled-DOF state in the same compact ordering returned
    !! by CD_Get_Model_CoupledDofs and accepted by CD_Update_Model_CoupledMotion.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: q_coupled(:), v_coupled(:), a_coupled(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: n, i

    q_coupled = CD_ZERO
    v_coupled = CD_ZERO
    a_coupled = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_MODEL_NOT_INITIALIZED
      ErrMsg = 'CableDyn_Model: model is not initialized'
      RETURN
    END IF
    n = SIZE(model%fixed_dofs)
    IF (SIZE(q_coupled) /= n .OR. SIZE(v_coupled) /= n .OR. SIZE(a_coupled) /= n) THEN
      CALL fail(ErrStat, ErrMsg, 'coupled motion outputs must have shape (n_coupled_dof)')
      RETURN
    END IF
    DO i = 1, n
      q_coupled(i) = model%q(model%fixed_dofs(i))
      v_coupled(i) = model%v(model%fixed_dofs(i))
      a_coupled(i) = model%a(model%fixed_dofs(i))
    END DO
  END SUBROUTINE CD_Get_Model_CoupledMotion

  LOGICAL FUNCTION CD_Model_Is_Initialized(model) RESULT(is_initialized)
    !! Query-only helper for shells and tests.
    TYPE(CD_ModelType), INTENT(IN) :: model
    is_initialized = model%initialized
  END FUNCTION CD_Model_Is_Initialized

  SUBROUTINE ensure_model_load_contributor_workspace(model, n, ErrStat, ErrMsg)
    !! Grow reusable dense contributor scratch used by model_contributor_load.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: istat
    LOGICAL :: need_workspace

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (n < 0) THEN
      CALL fail(ErrStat, ErrMsg, 'load contributor workspace size must be non-negative')
      RETURN
    END IF
    need_workspace = .NOT. ALLOCATED(model%load_force_work)
    IF (.NOT. need_workspace) need_workspace = SIZE(model%load_force_work) < n
    IF (need_workspace) THEN
      CALL clear_model_load_contributor_workspace(model)
      ALLOCATE (model%load_force_work(n), STAT=istat)
      IF (istat /= 0) THEN
        CALL load_workspace_alloc_fail(model, ErrStat, ErrMsg)
        RETURN
      END IF
      ALLOCATE (model%load_jq_work(n, n), STAT=istat)
      IF (istat /= 0) THEN
        CALL load_workspace_alloc_fail(model, ErrStat, ErrMsg)
        RETURN
      END IF
      ALLOCATE (model%load_jv_work(n, n), STAT=istat)
      IF (istat /= 0) THEN
        CALL load_workspace_alloc_fail(model, ErrStat, ErrMsg)
        RETURN
      END IF
    END IF
  END SUBROUTINE ensure_model_load_contributor_workspace

  SUBROUTINE load_workspace_alloc_fail(model, ErrStat, ErrMsg)
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    CALL clear_model_load_contributor_workspace(model)
    ErrStat = CD_MODEL_ALLOCFAIL
    ErrMsg = 'CableDyn_Model: load contributor workspace allocation failed'
  END SUBROUTINE load_workspace_alloc_fail

  SUBROUTINE clear_model_load_contributor_workspace(model)
    TYPE(CD_ModelType), INTENT(INOUT) :: model

    IF (ALLOCATED(model%load_force_work)) DEALLOCATE (model%load_force_work)
    IF (ALLOCATED(model%load_jq_work)) DEALLOCATE (model%load_jq_work)
    IF (ALLOCATED(model%load_jv_work)) DEALLOCATE (model%load_jv_work)
  END SUBROUTINE clear_model_load_contributor_workspace

  PURE LOGICAL FUNCTION mode3_marked(alpha_mbl, e) RESULT(marked)
    !! Safe probe of the OPTIONAL load-dependent marker array (nested guard --
    !! never short-circuit on an absent OPTIONAL).
    REAL(wp), INTENT(IN), OPTIONAL :: alpha_mbl(:)
    INTEGER, INTENT(IN) :: e
    marked = .FALSE.
    IF (PRESENT(alpha_mbl)) marked = alpha_mbl(e) > CD_ZERO
  END FUNCTION mode3_marked

  SUBROUTINE model_viscoelastic_coeffs(model, e, ea_d_eff, ea_1_eff, ErrStat, ErrMsg)
    !! The effective series-Kelvin stiffness pair of one viscoelastic element: mode 2
    !! returns the stored constants; mode 3 evaluates MoorDyn's mean-load-
    !! dependent EA_D at the COMMITTED dl_1 (lagged over the step -- constant
    !! within it, so the elimination Jacobians remain exact for the lagged
    !! operator) and derives EA_1 from it. Fails closed if the mode-3 premise
    !! is violated at the current state.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: e
    REAL(wp), INTENT(OUT) :: ea_d_eff, ea_1_eff
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(200) :: em
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    ea_d_eff = CD_ZERO
    ea_1_eff = CD_ZERO
    SELECT CASE (model%ve_mode(e))
    CASE (2)
      ea_d_eff = model%ve_ea_d(e)
      ea_1_eff = model%ve_ea_1(e)
    CASE (3)
      CALL CD_Viscoelastic_LoadDependent_EAD(model%ea(e), model%l0(e), model%ve_alpha_mbl(e), &
                                             model%ve_vbeta(e), model%ve_dl_1(e), ea_d_eff, es, em)
      IF (es /= CD_VISCO_OK) THEN
        ErrStat = CD_MODEL_SOLVEFAIL
        ErrMsg = 'CableDyn_Model: '//TRIM(em)
        RETURN
      END IF
      ea_1_eff = ea_d_eff*model%ea(e)/(ea_d_eff - model%ea(e))
    CASE DEFAULT
      CALL fail(ErrStat, ErrMsg, 'viscoelastic coefficients requested for a plain element')
    END SELECT
  END SUBROUTINE model_viscoelastic_coeffs

  SUBROUTINE model_syrope_element(model, e, q2, v2, force, jac_q, jac_v, t_mean, dslow, ErrStat, ErrMsg)
    !! One Syrope element's production MoorDyn force and position/velocity
    !! Jacobians. Regenerates the working curve from the committed running maximum,
    !! selects the branch by the evaluated slow strain vs the branch divider at
    !! that maximum, and calls the element load. Inside a step (ve_dt > 0) the
    !! slow strain is eliminated by the same backward-Euler solve the commit makes;
    !! outside (ve_dt = 0) the committed state is used as is.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: e
    REAL(wp), INTENT(IN) :: q2(6), v2(6)
    REAL(wp), INTENT(OUT) :: force(6), jac_q(6, 6), jac_v(6, 6), t_mean, dslow
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: ws(CD_SYROPE_NWC), wt(CD_SYROPE_NWC), wsl(CD_SYROPE_NWC), slow_at
    LOGICAL :: on_wc
    INTEGER :: es
    CHARACTER(200) :: em
    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    t_mean = CD_ZERO
    dslow = CD_ZERO
    CALL CD_Syrope_Working_Curve(model%syrope_type(e), model%syrope_tmax(e), ws, wt, wsl, es, em)
    IF (es /= CD_SYROPE_OK) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(em)); RETURN
    END IF
    ! the branch follows the evaluated slow strain against the continuous divider
    slow_at = CD_Syrope_Branch_Divider(model%syrope_type(e), model%syrope_tmax(e))
    on_wc = model%syrope_slow(e) < slow_at
    CALL CD_Syrope_Element_Load(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), model%syrope_type(e), &
                                ws, wt, wsl, on_wc, model%syrope_slow(e), model%tension_only, &
                                force, jac_q, t_mean, dslow, es, em, jac_v=jac_v, dt=model%ve_dt, &
                                divider=slow_at)
    IF (es /= CD_SYROPE_OK) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(em)); RETURN
    END IF
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END SUBROUTINE model_syrope_element

  SUBROUTINE model_syrope_state_next(model, e, q2, v2, dt, slow_next, t_mean, ErrStat, ErrMsg)
    !! Evaluate the backward-Euler Syrope state candidate without mutation. The
    !! caller invokes this for every element before the atomic commit pass. The
    !! candidate comes from the element load's own implicit elimination, so the
    !! committed slow strain is bit-identical to the one the converged in-step
    !! force was evaluated at. A candidate beyond the OWC table fails the step.
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: e
    REAL(wp), INTENT(IN) :: q2(6), v2(6), dt
    REAL(wp), INTENT(OUT) :: slow_next, t_mean
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: ws(CD_SYROPE_NWC), wt(CD_SYROPE_NWC), wsl(CD_SYROPE_NWC)
    REAL(wp) :: chord(3), length, eps, slow_at, f6(6), jq6(6, 6), dslow
    LOGICAL :: on_wc
    INTEGER :: es
    CHARACTER(200) :: em

    slow_next = model%syrope_slow(e)
    t_mean = CD_ZERO
    chord = q2(4:6) - q2(1:3)
    length = NORM2(chord)
    IF (length <= CD_VISCO_TINY_LEN) THEN
      CALL fail(ErrStat, ErrMsg, 'Syrope state advance found a collapsed element')
      RETURN
    END IF
    IF (.NOT. (dt > CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'Syrope state advance needs dt > 0')
      RETURN
    END IF
    eps = (length - model%l0(e))/model%l0(e)
    CALL CD_Syrope_Working_Curve(model%syrope_type(e), model%syrope_tmax(e), ws, wt, wsl, es, em)
    IF (es /= CD_SYROPE_OK) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(em)); RETURN
    END IF
    slow_at = CD_Syrope_Branch_Divider(model%syrope_type(e), model%syrope_tmax(e))
    on_wc = model%syrope_slow(e) < slow_at
    CALL CD_Syrope_Element_Load(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), model%syrope_type(e), &
                                ws, wt, wsl, on_wc, model%syrope_slow(e), model%tension_only, &
                                f6, jq6, t_mean, dslow, es, em, dt=dt, slow_next=slow_next, divider=slow_at)
    IF (es /= CD_SYROPE_OK) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(em)); RETURN
    END IF
    ! the interpolation clamps at the last OWC row, so a state past it would run on
    ! with a silently wrong tension and unbounded slow-strain creep: fail closed
    ! (the total-strain bound applies on the OWC branch of the new state)
    IF (slow_next < slow_at) THEN
      CALL CD_Syrope_Check_Range(model%syrope_type(e), slow_next, es, em)
    ELSE
      CALL CD_Syrope_Check_Range(model%syrope_type(e), slow_next, es, em, eps=eps)
    END IF
    IF (es /= CD_SYROPE_OK) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(em)); RETURN
    END IF
    IF (t_mean > model%syrope_tmax(e)) THEN
      CALL CD_Syrope_Working_Curve(model%syrope_type(e), t_mean, ws, wt, wsl, es, em)
      IF (es /= CD_SYROPE_OK) THEN
        CALL fail(ErrStat, ErrMsg, TRIM(em)); RETURN
      END IF
    END IF
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END SUBROUTINE model_syrope_state_next

  SUBROUTINE model_contributor_load(model, q, v, force, jac_q, jac_v, ErrStat, ErrMsg)
    !! Assemble all configured non-structural model load contributors and tangents.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: e, idx(6)
    REAL(wp) :: q2(6), v2(6), f2(6), jq2(6, 6), jv2(6, 6), ea_d_eff, ea_1_eff, syr_tm, syr_ds

    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    IF (SIZE(v) /= SIZE(q)) THEN
      ErrStat = CD_MODEL_BADINPUT
      ErrMsg = 'CableDyn_Model: model load received inconsistent q/v shapes'
      RETURN
    END IF
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    CALL ensure_model_load_contributor_workspace(model, SIZE(q), ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    ASSOCIATE (f_tmp => model%load_force_work(1:SIZE(q)), &
               jq_tmp => model%load_jq_work(1:SIZE(q), 1:SIZE(q)), &
               jv_tmp => model%load_jv_work(1:SIZE(q), 1:SIZE(q)))
      IF (model%has_seabed) THEN
        f_tmp = CD_ZERO
        jq_tmp = CD_ZERO
        IF (model%has_bathymetry) THEN
          CALL CD_Bathymetry_Seabed_Load(model%bathymetry, reshape3(q, SIZE(q)), model%seabed_kn, &
                                         f_tmp, jq_tmp, ErrStat, ErrMsg)
        ELSE
          CALL CD_Seabed_Penalty_Load(reshape3(q, SIZE(q)), model%seabed_kn, model%seabed_z_floor, &
                                      f_tmp, jq_tmp, ErrStat, ErrMsg)
        END IF
        IF (ErrStat /= 0) THEN
          ErrStat = CD_MODEL_BADINPUT   ! a foreign module's status mapped to the model's
          RETURN
        END IF
        force = force + f_tmp
        jac_q = jac_q + jq_tmp
      END IF
      IF (model%has_seabed_damping) THEN
        f_tmp = CD_ZERO
        jv_tmp = CD_ZERO
        jq_tmp = CD_ZERO
        IF (model%has_bathymetry) THEN
          CALL bathymetry_normal_damping_load(model, q, v, f_tmp, jq_tmp, jv_tmp, ErrStat, ErrMsg)
        ELSE
          CALL seabed_normal_damping_load(q, v, model%seabed_z_floor, model%seabed_cn, f_tmp, jq_tmp, jv_tmp, &
                                          ErrStat, ErrMsg)
        END IF
        IF (ErrStat /= 0) RETURN
        force = force + f_tmp
        jac_q = jac_q + jq_tmp
        jac_v = jac_v + jv_tmp
      END IF
      IF (model%has_seabed_friction) THEN
        f_tmp = CD_ZERO
        jq_tmp = CD_ZERO
        jv_tmp = CD_ZERO
        CALL seabed_friction_load(model, q, v, f_tmp, jq_tmp, jv_tmp, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        force = force + f_tmp
        jac_q = jac_q + jq_tmp
        jac_v = jac_v + jv_tmp
      END IF
      IF (model%has_damping) THEN
        f_tmp = CD_ZERO
        jq_tmp = CD_ZERO
        jv_tmp = CD_ZERO
        CALL CD_Cable_Axial_Damping_Load(q, v, model%elem_conn, model%l0, model%ba_dyn, &
                                         f_tmp, jq_tmp, jv_tmp, ErrStat, ErrMsg, &
                                         l0_dot=model%l0_dot)
        IF (ErrStat /= 0) RETURN
        force = force + f_tmp
        jac_q = jac_q + jq_tmp
        jac_v = jac_v + jv_tmp
      END IF
      IF (model%has_viscoelastic) THEN
        DO e = 1, model%n_elem
          IF (model%ve_mode(e) == 0) CYCLE
          CALL model_viscoelastic_coeffs(model, e, ea_d_eff, ea_1_eff, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL extract_element_state(model, e, q, v, idx, q2, v2)
          CALL CD_Viscoelastic_Element_Load(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), &
                                            ea_d_eff, ea_1_eff, model%ve_ba(e), &
                                            model%ve_ba_d(e), model%ve_dl_1(e), model%ve_dt, &
                                            f2, jq2, jv2, ErrStat, ErrMsg)
          IF (ErrStat /= 0) RETURN
          CALL scatter_element_vector(force, idx, f2)
          CALL scatter_element_matrix(jac_q, idx, jq2)
          CALL scatter_element_matrix(jac_v, idx, jv2)
        END DO
      END IF
      IF (model%has_syrope) THEN
        DO e = 1, model%n_elem
          IF (.NOT. model%syrope_is(e)) CYCLE
          CALL extract_element_state(model, e, q, v, idx, q2, v2)
          CALL model_syrope_element(model, e, q2, v2, f2, jq2, jv2, syr_tm, syr_ds, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL scatter_element_vector(force, idx, f2)
          CALL scatter_element_matrix(jac_q, idx, jq2)
          CALL scatter_element_matrix(jac_v, idx, jv2)
        END DO
      END IF
      IF (model%has_morison_drag) THEN
        f_tmp = CD_ZERO
        jq_tmp = CD_ZERO
        jv_tmp = CD_ZERO
        CALL model_morison_drag_load(model, q, v, f_tmp, jq_tmp, jv_tmp, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        force = force + f_tmp
        jac_q = jac_q + jq_tmp
        jac_v = jac_v + jv_tmp
      END IF
      IF (model%has_froude_krylov) THEN
        f_tmp = CD_ZERO
        jq_tmp = CD_ZERO
        jv_tmp = CD_ZERO
        CALL model_froude_krylov_load(model, q, v, f_tmp, jq_tmp, jv_tmp, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        force = force + f_tmp
        jac_q = jac_q + jq_tmp
        jac_v = jac_v + jv_tmp
      END IF
      IF (model%has_buoyancy_recovery) THEN
        f_tmp = CD_ZERO
        jq_tmp = CD_ZERO
        jv_tmp = CD_ZERO
        CALL model_buoyancy_recovery_load(model, q, v, f_tmp, jq_tmp, jv_tmp, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        force = force + f_tmp
        jac_q = jac_q + jq_tmp
        jac_v = jac_v + jv_tmp
      END IF
    END ASSOCIATE
  END SUBROUTINE model_contributor_load

  SUBROUTINE model_contributor_load_banded(model, q, v, free, kl, ku, force, jac_q_band, jac_v_band, &
                                           ErrStat, ErrMsg)
    !! Assemble model-owned sparse-local load contributors into full force and
    !! reduced free/free banded residual-Jacobian blocks. This is the production
    !! Newton path for no-added-mass EI=0 dynamic lines, including flat and structured
    !! bathymetry contact/damping/friction.
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:)
    INTEGER, INTENT(IN) :: free(:), kl, ku
    REAL(wp), INTENT(OUT) :: force(:), jac_q_band(:, :), jac_v_band(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: n, ldab, e, i, ix, iy, iz, idx(6)
    REAL(wp) :: gap, normal, g(2), z_floor, dfdx, dfdy, tangent, dn_dv(3), dfc(2), nvec(3), vn
    LOGICAL :: sloped
    REAL(wp) :: d_fr(2), mu_c
    REAL(wp) :: damp_force, damp_dgap, damp_dvz
    REAL(wp) :: q2(6), v2(6), f2(6), jq2(6, 6), jv2(6, 6), ea_d_eff, ea_1_eff, syr_tm, syr_ds
    LOGICAL :: free_map_ready
    REAL(wp) :: fluid2(3, 2), accel2(3, 2), wl2(2), dfdv(2, 2)

    force = CD_ZERO
    jac_q_band = CD_ZERO
    jac_v_band = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    n = SIZE(q)
    ldab = 2*kl + ku + 1
    IF (SIZE(v) /= n .OR. SIZE(force) /= n) THEN
      CALL fail(ErrStat, ErrMsg, 'banded model load received inconsistent q/v/force shapes')
      RETURN
    END IF
    IF (kl < 0 .OR. ku < 0 .OR. SIZE(jac_q_band, 1) < ldab .OR. SIZE(jac_v_band, 1) < ldab .OR. &
        SIZE(jac_q_band, 2) /= SIZE(free) .OR. SIZE(jac_v_band, 2) /= SIZE(free)) THEN
      CALL fail(ErrStat, ErrMsg, 'banded model load received invalid band storage')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(q) .OR. .NOT. CD_All_Finite(v)) THEN
      CALL fail(ErrStat, ErrMsg, 'banded model load q/v must be finite')
      RETURN
    END IF

    free_map_ready = ALLOCATED(model%free_map_work)
    IF (free_map_ready) free_map_ready = SIZE(model%free_map_work) >= n
    IF (.NOT. free_map_ready) THEN
      CALL fail(ErrStat, ErrMsg, 'banded model load requires initialized free-map workspace')
      RETURN
    END IF
    model%free_map_work(1:n) = 0
    DO i = 1, SIZE(free)
      IF (free(i) < 1 .OR. free(i) > n) THEN
        CALL fail(ErrStat, ErrMsg, 'banded model load free DOF index out of range')
        RETURN
      END IF
      IF (model%free_map_work(free(i)) /= 0) THEN
        CALL fail(ErrStat, ErrMsg, 'banded model load duplicate free DOF')
        RETURN
      END IF
      model%free_map_work(free(i)) = i
    END DO

    IF (model%has_seabed) THEN
      DO i = 1, n/3
        ix = 3*i - 2
        iy = 3*i - 1
        iz = 3*i
        IF (model%has_bathymetry) THEN
          CALL model_floor_gradient(model, q(ix), q(iy), z_floor, dfdx, dfdy, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
        ELSE
          z_floor = model%seabed_z_floor
          dfdx = CD_ZERO
          dfdy = CD_ZERO
        END IF
        gap = z_floor - q(iz)
        IF (ABS(dfdx) > CD_ZERO .OR. ABS(dfdy) > CD_ZERO) THEN
          ! sloped seabed: frictionless contact along the surface normal (the static
          ! initializer uses the same law, so the initial condition is an equilibrium)
          BLOCK
            REAL(wp) :: fvec(3), jac3(3, 3)
            INTEGER :: ra, cb
            CALL CD_Seabed_Normal_Contact(gap, dfdx, dfdy, model%seabed_kn(i), fvec, jac3)
            force(ix:iz) = force(ix:iz) + fvec
            DO cb = 1, 3
              DO ra = 1, 3
                CALL scatter_scalar_band(jac_q_band, model%free_map_work, kl, ku, ix + ra - 1, ix + cb - 1, &
                                         jac3(ra, cb), ErrStat, ErrMsg)
                IF (ErrStat /= CD_MODEL_OK) RETURN
              END DO
            END DO
          END BLOCK
          CYCLE
        END IF
        CALL CD_Seabed_Normal_Law(gap, model%seabed_kn(i), normal, tangent)
        force(iz) = force(iz) + normal
        CALL scatter_scalar_band(jac_q_band, model%free_map_work, kl, ku, iz, iz, tangent, &
                                 ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END DO
    END IF

    IF (model%has_seabed_damping) THEN
      DO i = 1, n/3
        ix = 3*i - 2
        iy = 3*i - 1
        iz = 3*i
        ! damping along the floor normal, at the normal penetration and normal velocity
        CALL node_contact_normal(model, q(ix), q(iy), q(iz), v(ix:iz), gap, vn, nvec, sloped, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL seabed_normal_damping_terms(gap, vn, model%seabed_cn(i), damp_force, damp_dgap, damp_dvz)
        IF (.NOT. sloped) THEN
          force(iz) = force(iz) + damp_force
          CALL scatter_scalar_band(jac_q_band, model%free_map_work, kl, ku, iz, iz, damp_dgap, &
                                   ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL scatter_scalar_band(jac_v_band, model%free_map_work, kl, ku, iz, iz, -damp_dvz, &
                                   ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CYCLE
        END IF
        force(ix:iz) = force(ix:iz) + damp_force*nvec
        ! -dF/dq = damp_dgap n n^T (d gap_n/dq = -n); -dF/dv = -damp_dvz n n^T
        BLOCK
          INTEGER :: ra, cb
          DO cb = 1, 3
            DO ra = 1, 3
              CALL scatter_scalar_band(jac_q_band, model%free_map_work, kl, ku, ix + ra - 1, ix + cb - 1, &
                                       damp_dgap*nvec(ra)*nvec(cb), ErrStat, ErrMsg)
              IF (ErrStat /= CD_MODEL_OK) RETURN
              CALL scatter_scalar_band(jac_v_band, model%free_map_work, kl, ku, ix + ra - 1, ix + cb - 1, &
                                       -damp_dvz*nvec(ra)*nvec(cb), ErrStat, ErrMsg)
              IF (ErrStat /= CD_MODEL_OK) RETURN
            END DO
          END DO
        END BLOCK
      END DO
    END IF

    IF (model%has_seabed_friction) THEN
      DO i = 1, n/3
        ix = 3*i - 2
        iy = 3*i - 1
        iz = 3*i
        normal = seabed_normal_reaction(model, q(ix), q(iy), q(iz), v(ix:iz), i, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        IF (normal <= CD_ZERO) CYCLE
        ! stick-slip spring to the committed anchor; the load on the node is -f
        d_fr(1) = q(ix) - model%fr_anchor(1, i)
        d_fr(2) = q(iy) - model%fr_anchor(2, i)
        IF (model%seabed_friction_aniso) THEN
          ! dfc is then d(force)/d(normal), so the capacity chain factor is one
          CALL CD_Seabed_Friction_Aniso(CD_FRICTION_STICK_SLIP, d_fr, model%seabed_kn(i), normal, &
                                        model%seabed_mu_axial, model%seabed_mu, friction_axis(q, i), g, dfdv, dfc)
          mu_c = CD_ONE
        ELSE
          CALL CD_Seabed_Friction_Stick_Slip(d_fr, model%seabed_kn(i), &
                                             model%seabed_mu*normal, g, dfdv, dfc)
          mu_c = model%seabed_mu
        END IF
        force(ix) = force(ix) - g(1)
        force(iy) = force(iy) - g(2)
        CALL scatter_scalar_band(jac_q_band, model%free_map_work, kl, ku, ix, ix, dfdv(1, 1), ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL scatter_scalar_band(jac_q_band, model%free_map_work, kl, ku, ix, iy, dfdv(1, 2), ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL scatter_scalar_band(jac_q_band, model%free_map_work, kl, ku, iy, ix, dfdv(2, 1), ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL scatter_scalar_band(jac_q_band, model%free_map_work, kl, ku, iy, iy, dfdv(2, 2), ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        IF (ABS(dfc(1)) + ABS(dfc(2)) <= CD_ZERO) CYCLE
        ! sliding: the capacity follows the normal reaction
        CALL scatter_seabed_friction_position_band(model, q(ix), q(iy), q(iz), v(ix:iz), i, ix, dfc(1), jac_q_band, &
                                                   model%free_map_work, kl, ku, ErrStat, ErrMsg, mu_c)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL scatter_seabed_friction_position_band(model, q(ix), q(iy), q(iz), v(ix:iz), i, iy, dfc(2), jac_q_band, &
                                                   model%free_map_work, kl, ku, ErrStat, ErrMsg, mu_c)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL seabed_normal_velocity_gradient(model, q(ix), q(iy), q(iz), v(ix:iz), i, dn_dv, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        BLOCK
          INTEGER :: rr, cc
          DO rr = 1, 2
            DO cc = 1, 3
              ! a level floor has dn_dv = (0, 0, d): skip its structural zeros
              IF (cc < 3 .AND. .NOT. (ABS(dn_dv(cc)) > CD_ZERO)) CYCLE
              CALL scatter_scalar_band(jac_v_band, model%free_map_work, kl, ku, ix + rr - 1, ix + cc - 1, &
                                       mu_c*dn_dv(cc)*dfc(rr), ErrStat, ErrMsg)
              IF (ErrStat /= CD_MODEL_OK) RETURN
            END DO
          END DO
        END BLOCK
      END DO
    END IF

    DO e = 1, model%n_elem
      CALL extract_element_state(model, e, q, v, idx, q2, v2)
      IF (model%has_damping) THEN
        CALL CD_Axial_Damping_Element_Load(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), &
                                           model%ba_dyn(e), f2, jq2, jv2, ErrStat, ErrMsg, &
                                           l0_dot=model%l0_dot(e))
        IF (ErrStat /= 0) RETURN
        CALL scatter_element_vector(force, idx, f2)
        CALL scatter_element_matrix_band(jac_q_band, model%free_map_work, kl, ku, idx, jq2, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL scatter_element_matrix_band(jac_v_band, model%free_map_work, kl, ku, idx, jv2, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END IF
      IF (model%has_viscoelastic) THEN
        IF (model%ve_mode(e) /= 0) THEN
          CALL model_viscoelastic_coeffs(model, e, ea_d_eff, ea_1_eff, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL CD_Viscoelastic_Element_Load(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), &
                                            ea_d_eff, ea_1_eff, model%ve_ba(e), &
                                            model%ve_ba_d(e), model%ve_dl_1(e), model%ve_dt, &
                                            f2, jq2, jv2, ErrStat, ErrMsg)
          IF (ErrStat /= 0) RETURN
          CALL scatter_element_vector(force, idx, f2)
          CALL scatter_element_matrix_band(jac_q_band, model%free_map_work, kl, ku, idx, jq2, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL scatter_element_matrix_band(jac_v_band, model%free_map_work, kl, ku, idx, jv2, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
        END IF
      END IF
      IF (model%has_syrope) THEN
        IF (model%syrope_is(e)) THEN
          CALL model_syrope_element(model, e, q2, v2, f2, jq2, jv2, syr_tm, syr_ds, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL scatter_element_vector(force, idx, f2)
          CALL scatter_element_matrix_band(jac_q_band, model%free_map_work, kl, ku, idx, jq2, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL scatter_element_matrix_band(jac_v_band, model%free_map_work, kl, ku, idx, jv2, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
        END IF
      END IF
      IF (model%has_morison_drag) THEN
        CALL extract_nodal_vector(model%fluid_velocity, model, e, fluid2)
        CALL extract_nodal_scalar(model%drag_waterline_z, model, e, wl2)
        CALL CD_Morison_Drag_Element_Load(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), &
                                          fluid2(:, 1), fluid2(:, 2), wl2(1), wl2(2), model%drag_rho, &
                                          model%drag_diameter_elem(e), model%drag_cdn_elem(e), &
                                          model%drag_cdt_elem(e), f2, jq2, jv2, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        CALL scatter_element_vector(force, idx, f2)
        CALL scatter_element_matrix_band(jac_q_band, model%free_map_work, kl, ku, idx, jq2, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL scatter_element_matrix_band(jac_v_band, model%free_map_work, kl, ku, idx, jv2, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END IF
      IF (model%has_froude_krylov) THEN
        CALL extract_nodal_vector(model%fluid_acceleration, model, e, accel2)
        CALL extract_nodal_scalar(model%fk_waterline_z, model, e, wl2)
        CALL CD_Froude_Krylov_Element_Load(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), &
                                           accel2(:, 1), accel2(:, 2), wl2(1), wl2(2), model%fk_rho, &
                                           model%fk_diameter_elem(e), model%fk_can_elem(e), &
                                           model%fk_cat_elem(e), f2, jq2, jv2, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        CALL scatter_element_vector(force, idx, f2)
        CALL scatter_element_matrix_band(jac_q_band, model%free_map_work, kl, ku, idx, jq2, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL scatter_element_matrix_band(jac_v_band, model%free_map_work, kl, ku, idx, jv2, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END IF
      IF (model%has_buoyancy_recovery) THEN
        CALL extract_nodal_scalar(model%buoyancy_waterline_z, model, e, wl2)
        CALL CD_Buoyancy_Recovery_Element_Load(q2(1:3), q2(4:6), model%l0(e), wl2(1), wl2(2), &
                                               model%buoyancy_rho, model%buoyancy_diameter_elem(e), &
                                               model%buoyancy_gravity, f2, jq2, jv2, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        CALL scatter_element_vector(force, idx, f2)
        CALL scatter_element_matrix_band(jac_q_band, model%free_map_work, kl, ku, idx, jq2, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        CALL scatter_element_matrix_band(jac_v_band, model%free_map_work, kl, ku, idx, jv2, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
      END IF
    END DO

    IF (.NOT. CD_All_Finite(force) .OR. .NOT. CD_All_Finite(jac_q_band) .OR. &
        .NOT. CD_All_Finite(jac_v_band)) THEN
      CALL fail(ErrStat, ErrMsg, 'banded model load returned non-finite force/Jacobian')
    END IF
  END SUBROUTINE model_contributor_load_banded

  SUBROUTINE model_contributor_force(model, q, v, force, ErrStat, ErrMsg, end_nodes_only)
    !! Assemble all configured non-structural model load contributors without tangents.
    !! end_nodes_only: only the first and last node's entries are needed (end-force recovery);
    !! the other nodes and the elements not touching an end are skipped, and the two end entries
    !! are bit-identical to the full assembly (same terms in the same order). `force` is then the
    !! compact 9-vector [first node (3), last node (3), sum over the other nodes the end
    !! elements reach (3)], so the recovery needs no n_dof-long scratch; the third block only
    !! carries the finiteness check of the full assembly.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    LOGICAL, INTENT(IN), OPTIONAL :: end_nodes_only

    INTEGER :: n, e, i, ix, iy, iz, idx(6), nn
    LOGICAL :: ends
    REAL(wp) :: gap, z_floor, normal, tangent, damp_force, damp_dgap, damp_dvz, fr_f(2), fr_dd(2, 2), fr_dc(2)
    REAL(wp) :: d_fr(2)
    REAL(wp) :: q2(6), v2(6), f2(6), ea_d_eff, ea_1_eff, syr_tm, syr_ds, syr_jac(6, 6), syr_jac_v(6, 6)
    REAL(wp) :: fluid2(3, 2), accel2(3, 2), wl2(2)

    force = CD_ZERO
    n = SIZE(q)
    ends = .FALSE.
    IF (PRESENT(end_nodes_only)) ends = end_nodes_only
    IF (SIZE(v) /= n .OR. SIZE(force) /= MERGE(9, n, ends)) THEN
      ErrStat = CD_MODEL_BADINPUT
      ErrMsg = 'CableDyn_Model: model force received inconsistent q/v/force shapes'
      RETURN
    END IF
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    nn = n/3
    IF (.NOT. CD_All_Finite(q) .OR. .NOT. CD_All_Finite(v)) THEN
      CALL fail(ErrStat, ErrMsg, 'model force q/v must be finite')
      RETURN
    END IF

    IF (model%has_seabed) THEN
      DO i = 1, n/3
        IF (ends .AND. i /= 1 .AND. i /= nn) CYCLE
        ix = 3*i - 2
        iy = 3*i - 1
        iz = 3*i
        IF (model%has_bathymetry) THEN
          ! sloped seabed: frictionless contact along the surface normal
          BLOCK
            REAL(wp) :: fvec(3), jac3(3, 3), gx, gy
            CALL model_floor_gradient(model, q(ix), q(iy), z_floor, gx, gy, ErrStat, ErrMsg)
            IF (ErrStat /= CD_MODEL_OK) RETURN
            IF (ABS(gx) > CD_ZERO .OR. ABS(gy) > CD_ZERO) THEN
              CALL CD_Seabed_Normal_Contact(z_floor - q(iz), gx, gy, model%seabed_kn(i), fvec, jac3)
              force(slot(ix):slot(iz)) = force(slot(ix):slot(iz)) + fvec
              CYCLE
            END IF
          END BLOCK
        ELSE
          z_floor = model%seabed_z_floor
        END IF
        gap = z_floor - q(iz)
        CALL CD_Seabed_Normal_Law(gap, model%seabed_kn(i), normal, tangent)
        force(slot(iz)) = force(slot(iz)) + normal
      END DO
    END IF

    IF (model%has_seabed_damping) THEN
      DO i = 1, n/3
        IF (ends .AND. i /= 1 .AND. i /= nn) CYCLE
        ix = 3*i - 2
        iy = 3*i - 1
        iz = 3*i
        BLOCK
          REAL(wp) :: vn, nvec(3)
          LOGICAL :: sloped
          CALL node_contact_normal(model, q(ix), q(iy), q(iz), v(ix:iz), gap, vn, nvec, sloped, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL seabed_normal_damping_terms(gap, vn, model%seabed_cn(i), damp_force, damp_dgap, damp_dvz)
          IF (sloped) THEN
            force(slot(ix):slot(iz)) = force(slot(ix):slot(iz)) + damp_force*nvec
          ELSE
            force(slot(iz)) = force(slot(iz)) + damp_force
          END IF
        END BLOCK
      END DO
    END IF

    IF (model%has_seabed_friction) THEN
      DO i = 1, n/3
        IF (ends .AND. i /= 1 .AND. i /= nn) CYCLE
        ix = 3*i - 2
        iy = 3*i - 1
        iz = 3*i
        normal = seabed_normal_reaction(model, q(ix), q(iy), q(iz), v(ix:iz), i, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        IF (normal <= CD_ZERO) CYCLE
        d_fr(1) = q(ix) - model%fr_anchor(1, i)
        d_fr(2) = q(iy) - model%fr_anchor(2, i)
        IF (model%seabed_friction_aniso) THEN
          CALL CD_Seabed_Friction_Aniso(CD_FRICTION_STICK_SLIP, d_fr, model%seabed_kn(i), normal, &
                                        model%seabed_mu_axial, model%seabed_mu, friction_axis(q, i), &
                                        fr_f, fr_dd, fr_dc)
        ELSE
          CALL CD_Seabed_Friction_Stick_Slip(d_fr, model%seabed_kn(i), &
                                             model%seabed_mu*normal, fr_f, fr_dd, fr_dc)
        END IF
        force(slot(ix)) = force(slot(ix)) - fr_f(1)
        force(slot(iy)) = force(slot(iy)) - fr_f(2)
      END DO
    END IF

    DO e = 1, model%n_elem
      IF (ends) THEN
        IF (ALL(model%elem_conn(:, e) /= 1) .AND. ALL(model%elem_conn(:, e) /= nn)) CYCLE
      END IF
      CALL extract_element_state(model, e, q, v, idx, q2, v2)
      IF (ends) idx = [(slot(idx(i)), i=1, 6)]
      IF (model%has_damping) THEN
        CALL CD_Axial_Damping_Element_Force(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), &
                                            model%ba_dyn(e), f2, ErrStat, ErrMsg, &
                                            l0_dot=model%l0_dot(e))
        IF (ErrStat /= 0) RETURN
        CALL scatter_element_vector(force, idx, f2)
      END IF
      IF (model%has_viscoelastic) THEN
        IF (model%ve_mode(e) /= 0) THEN
          CALL model_viscoelastic_coeffs(model, e, ea_d_eff, ea_1_eff, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL CD_Viscoelastic_Element_Force(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), &
                                             ea_d_eff, ea_1_eff, model%ve_ba(e), &
                                             model%ve_ba_d(e), model%ve_dl_1(e), model%ve_dt, &
                                             f2, ErrStat, ErrMsg)
          IF (ErrStat /= 0) RETURN
          CALL scatter_element_vector(force, idx, f2)
        END IF
      END IF
      IF (model%has_syrope) THEN
        IF (model%syrope_is(e)) THEN
          CALL model_syrope_element(model, e, q2, v2, f2, syr_jac, syr_jac_v, syr_tm, syr_ds, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
          CALL scatter_element_vector(force, idx, f2)
        END IF
      END IF
      IF (model%has_morison_drag) THEN
        CALL extract_nodal_vector(model%fluid_velocity, model, e, fluid2)
        CALL extract_nodal_scalar(model%drag_waterline_z, model, e, wl2)
        CALL CD_Morison_Drag_Element_Force(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), &
                                           fluid2(:, 1), fluid2(:, 2), wl2(1), wl2(2), model%drag_rho, &
                                           model%drag_diameter_elem(e), model%drag_cdn_elem(e), &
                                           model%drag_cdt_elem(e), f2, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        CALL scatter_element_vector(force, idx, f2)
      END IF
      IF (model%has_froude_krylov) THEN
        CALL extract_nodal_vector(model%fluid_acceleration, model, e, accel2)
        CALL extract_nodal_scalar(model%fk_waterline_z, model, e, wl2)
        CALL CD_Froude_Krylov_Element_Force(q2(1:3), q2(4:6), v2(1:3), v2(4:6), model%l0(e), &
                                            accel2(:, 1), accel2(:, 2), wl2(1), wl2(2), model%fk_rho, &
                                            model%fk_diameter_elem(e), model%fk_can_elem(e), &
                                            model%fk_cat_elem(e), f2, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        CALL scatter_element_vector(force, idx, f2)
      END IF
      IF (model%has_buoyancy_recovery) THEN
        CALL extract_nodal_scalar(model%buoyancy_waterline_z, model, e, wl2)
        CALL CD_Buoyancy_Recovery_Element_Force(q2(1:3), q2(4:6), model%l0(e), wl2(1), wl2(2), &
                                                model%buoyancy_rho, model%buoyancy_diameter_elem(e), &
                                                model%buoyancy_gravity, f2, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
        CALL scatter_element_vector(force, idx, f2)
      END IF
    END DO
    IF (.NOT. CD_All_Finite(force)) CALL fail(ErrStat, ErrMsg, 'model force returned non-finite force')

  CONTAINS

    PURE INTEGER FUNCTION slot(k) RESULT(j)
      !! Storage index of global DOF k: itself in a full assembly, its compact slot in an
      !! end-node recovery.
      INTEGER, INTENT(IN) :: k
      INTEGER :: node
      j = k
      IF (.NOT. ends) RETURN
      node = (k + 2)/3
      j = k - 3*(node - 1)
      IF (node == nn) THEN
        j = j + 3
      ELSE IF (node /= 1) THEN
        j = j + 6
      END IF
    END FUNCTION slot
  END SUBROUTINE model_contributor_force

  SUBROUTINE seabed_normal_damping_load(q, v, z_floor, c_n, force, jac_q, jac_v, ErrStat, ErrMsg)
    !! Contact-only normal seabed damping with the shared C1 touchdown gate.
    !! For fully contacted downward motion, force_z = -c_n*v_z; through the
    !! touchdown blend the force and position tangent ramp continuously.
    REAL(wp), INTENT(IN) :: q(:), v(:), z_floor, c_n(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: i, iz, n_nodes
    REAL(wp) :: gap, damp_force, damp_dgap, damp_dvz

    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    CALL validate_seabed_damping(q, v, z_floor, c_n, force, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    IF (SIZE(jac_q, 1) /= SIZE(q) .OR. SIZE(jac_q, 2) /= SIZE(q) .OR. &
        SIZE(jac_v, 1) /= SIZE(q) .OR. SIZE(jac_v, 2) /= SIZE(q)) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed damping Jacobian shapes are inconsistent')
      RETURN
    END IF
    n_nodes = SIZE(q)/3
    DO i = 1, n_nodes
      iz = 3*i
      gap = z_floor - q(iz)
      CALL seabed_normal_damping_terms(gap, v(iz), c_n(i), damp_force, damp_dgap, damp_dvz)
      force(iz) = damp_force
      jac_q(iz, iz) = damp_dgap
      jac_v(iz, iz) = -damp_dvz
    END DO
  END SUBROUTINE seabed_normal_damping_load

  SUBROUTINE validate_seabed_damping(q, v, z_floor, c_n, force, ErrStat, ErrMsg)
    REAL(wp), INTENT(IN) :: q(:), v(:), z_floor, c_n(:)
    REAL(wp), INTENT(IN) :: force(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n_nodes

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (MOD(SIZE(q), 3) /= 0 .OR. SIZE(q) < 6) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed damping q must be a positions-only state')
      RETURN
    END IF
    n_nodes = SIZE(q)/3
    IF (SIZE(v) /= SIZE(q) .OR. SIZE(force) /= SIZE(q) .OR. SIZE(c_n) /= n_nodes) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed damping q/v/force/c_n shapes are inconsistent')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(q) .OR. .NOT. CD_All_Finite(v) .OR. &
        .NOT. CD_All_Finite(c_n) .OR. .NOT. CD_Is_Finite(z_floor)) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed damping inputs must be finite')
      RETURN
    END IF
    IF (ANY(c_n < CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed damping coefficients must be non-negative')
      RETURN
    END IF
  END SUBROUTINE validate_seabed_damping

  SUBROUTINE bathymetry_normal_damping_load(model, q, v, force, jac_q, jac_v, ErrStat, ErrMsg)
    !! One-sided vertical normal damping against a structured bathymetry floor,
    !! with the same C1 touchdown gate used for flat-floor contact.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: i, ix, iy, iz, n_nodes, cb
    REAL(wp) :: gap, vn, nvec(3), damp_force, damp_dgap, damp_dvz
    LOGICAL :: sloped

    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    CALL validate_bathymetry_damping(model, q, v, force, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    IF (SIZE(jac_q, 1) /= SIZE(q) .OR. SIZE(jac_q, 2) /= SIZE(q) .OR. &
        SIZE(jac_v, 1) /= SIZE(q) .OR. SIZE(jac_v, 2) /= SIZE(q)) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry damping Jacobian shapes are inconsistent')
      RETURN
    END IF
    n_nodes = SIZE(q)/3
    DO i = 1, n_nodes
      ix = 3*i - 2
      iy = ix + 1
      iz = ix + 2
      ! damping along the floor normal, at the normal penetration and normal velocity
      CALL node_contact_normal(model, q(ix), q(iy), q(iz), v(ix:iz), gap, vn, nvec, sloped, ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
      CALL seabed_normal_damping_terms(gap, vn, model%seabed_cn(i), damp_force, damp_dgap, damp_dvz)
      IF (.NOT. sloped) THEN
        force(iz) = damp_force
        jac_q(iz, iz) = damp_dgap
        jac_v(iz, iz) = -damp_dvz
        CYCLE
      END IF
      force(ix:iz) = damp_force*nvec
      DO cb = 1, 3
        jac_q(ix:iz, ix + cb - 1) = damp_dgap*nvec*nvec(cb)
        jac_v(ix:iz, ix + cb - 1) = -damp_dvz*nvec*nvec(cb)
      END DO
    END DO
  END SUBROUTINE bathymetry_normal_damping_load

  SUBROUTINE validate_bathymetry_damping(model, q, v, force, ErrStat, ErrMsg)
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:), force(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n_nodes

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%has_bathymetry .OR. .NOT. CD_Bathymetry_Is_Initialized(model%bathymetry)) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry damping requires initialized bathymetry contact')
      RETURN
    END IF
    IF (MOD(SIZE(q), 3) /= 0 .OR. SIZE(q) < 6) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry damping q must be a positions-only state')
      RETURN
    END IF
    n_nodes = SIZE(q)/3
    IF (SIZE(v) /= SIZE(q) .OR. SIZE(force) /= SIZE(q) .OR. SIZE(model%seabed_cn) /= n_nodes) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry damping q/v/force/c_n shapes are inconsistent')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(q) .OR. .NOT. CD_All_Finite(v) .OR. &
        .NOT. CD_All_Finite(model%seabed_cn) .OR. ANY(model%seabed_cn < CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry damping inputs must be finite with c_n >= 0')
      RETURN
    END IF
  END SUBROUTINE validate_bathymetry_damping

  SUBROUTINE seabed_friction_load(model, q, v, force, jac_q, jac_v, ErrStat, ErrMsg)
    !! Stick-slip seabed friction: each node in contact holds a horizontal spring of
    !! stiffness k_n to its committed anchor, capped at mu*N (CD_Seabed_Friction_Stick_Slip),
    !! where N is the active spring normal reaction plus one-sided normal damping.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: i, ix, iy, iz, n_nodes, r
    REAL(wp) :: normal, f(2), dfd(2, 2), dfc(2), dn_dq(3), dn_dv(3)
    REAL(wp) :: d_fr(2), mu_c

    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    dn_dq = CD_ZERO
    dn_dv = CD_ZERO
    CALL validate_seabed_friction(model, q, v, force, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    n_nodes = SIZE(q)/3
    DO i = 1, n_nodes
      ix = 3*i - 2
      iy = 3*i - 1
      iz = 3*i
      normal = seabed_normal_reaction(model, q(ix), q(iy), q(iz), v(ix:iz), i, ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
      IF (normal <= CD_ZERO) CYCLE
      d_fr(1) = q(ix) - model%fr_anchor(1, i)
      d_fr(2) = q(iy) - model%fr_anchor(2, i)
      IF (model%seabed_friction_aniso) THEN
        CALL CD_Seabed_Friction_Aniso(CD_FRICTION_STICK_SLIP, d_fr, model%seabed_kn(i), normal, &
                                      model%seabed_mu_axial, model%seabed_mu, friction_axis(q, i), f, dfd, dfc)
        mu_c = CD_ONE
      ELSE
        CALL CD_Seabed_Friction_Stick_Slip(d_fr, model%seabed_kn(i), &
                                           model%seabed_mu*normal, f, dfd, dfc)
        mu_c = model%seabed_mu
      END IF
      force(ix) = -f(1)
      force(iy) = -f(2)
      jac_q(ix:iy, ix:iy) = jac_q(ix:iy, ix:iy) + dfd
      IF (ABS(dfc(1)) + ABS(dfc(2)) <= CD_ZERO) CYCLE
      CALL seabed_normal_position_gradient(model, q(ix), q(iy), q(iz), v(ix:iz), i, dn_dq, ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
      CALL seabed_normal_velocity_gradient(model, q(ix), q(iy), q(iz), v(ix:iz), i, dn_dv, ErrStat, ErrMsg)
      IF (ErrStat /= CD_MODEL_OK) RETURN
      DO r = 1, 2
        jac_q(ix + r - 1, ix:iz) = jac_q(ix + r - 1, ix:iz) + mu_c*dfc(r)*dn_dq
        jac_v(ix + r - 1, ix:iz) = jac_v(ix + r - 1, ix:iz) + mu_c*dfc(r)*dn_dv
      END DO
    END DO
  END SUBROUTINE seabed_friction_load

  SUBROUTINE scatter_seabed_friction_position_band(model, x, y, z, vel, node, row, g_component, jac_q_band, &
                                                   free_map, kl, ku, ErrStat, ErrMsg, mu_c)
    !! mu_c: the factor that turns g_component into d(force)/d(normal) (seabed_mu for the
    !! isotropic law, whose g_component is d(force)/d(capacity); one for the anisotropic law).
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: x, y, z, vel(3), g_component, mu_c
    INTEGER, INTENT(IN) :: node, row, free_map(:), kl, ku
    REAL(wp), INTENT(INOUT) :: jac_q_band(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: ix, iy, iz
    REAL(wp) :: dn_dq(3)

    ix = 3*node - 2
    iy = ix + 1
    iz = ix + 2
    dn_dq = CD_ZERO
    CALL seabed_normal_position_gradient(model, x, y, z, vel, node, dn_dq, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL scatter_scalar_band(jac_q_band, free_map, kl, ku, row, ix, &
                             mu_c*dn_dq(1)*g_component, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL scatter_scalar_band(jac_q_band, free_map, kl, ku, row, iy, &
                             mu_c*dn_dq(2)*g_component, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL scatter_scalar_band(jac_q_band, free_map, kl, ku, row, iz, &
                             mu_c*dn_dq(3)*g_component, ErrStat, ErrMsg)
  END SUBROUTINE scatter_seabed_friction_position_band

  SUBROUTINE validate_seabed_friction(model, q, v, force, ErrStat, ErrMsg)
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:), force(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n_nodes

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%has_seabed) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed friction requires seabed contact')
      RETURN
    END IF
    IF (MOD(SIZE(q), 3) /= 0 .OR. SIZE(q) < 6) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed friction q must be a positions-only state')
      RETURN
    END IF
    n_nodes = SIZE(q)/3
    IF (SIZE(v) /= SIZE(q) .OR. SIZE(force) /= SIZE(q) .OR. SIZE(model%seabed_kn) /= n_nodes) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed friction state/contact shapes are inconsistent')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(q) .OR. .NOT. CD_All_Finite(v) .OR. &
        .NOT. CD_Is_Finite(model%seabed_mu) .OR. model%seabed_mu < CD_ZERO) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed friction inputs must be finite with mu >= 0')
      RETURN
    END IF
    IF (.NOT. ALLOCATED(model%fr_anchor)) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed friction anchors are not initialized')
      RETURN
    END IF
    IF (SIZE(model%fr_anchor, 2) /= n_nodes) THEN
      CALL fail(ErrStat, ErrMsg, 'seabed friction anchors must have shape (2, n_nodes)')
      RETURN
    END IF
  END SUBROUTINE validate_seabed_friction

  SUBROUTINE model_floor_gradient(model, x, y, z_floor, dfdx, dfdy, ErrStat, ErrMsg)
    !! Seabed elevation and gradient of the model's bathymetry, with the bathymetry status
    !! mapped to the model's (CD_MODEL_BADINPUT on a failed query).
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: x, y
    REAL(wp), INTENT(OUT) :: z_floor, dfdx, dfdy
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    CHARACTER(160) :: em
    CALL CD_Bathymetry_Floor_Gradient(model%bathymetry, x, y, z_floor, dfdx, dfdy, es, em)
    IF (es == CD_BATHY_OK) THEN
      ErrStat = CD_MODEL_OK
      ErrMsg = ''
    ELSE
      ErrStat = CD_MODEL_BADINPUT
      ErrMsg = 'CableDyn_Model: bathymetry gradient query failed: '//TRIM(em)
    END IF
  END SUBROUTINE model_floor_gradient

  SUBROUTINE node_contact_normal(model, x, y, z, vel, gap_n, vn, nvec, sloped, ErrStat, ErrMsg)
    !! Contact geometry of a node against the seabed: the penetration normal to the
    !! (locally planar) floor gap_n = (z_floor - z)/s, the velocity along the upward unit
    !! normal vn = v . n and n = (-dfdx, -dfdy, 1)/s, s = sqrt(1 + dfdx^2 + dfdy^2). On a
    !! level floor n = e_z, gap_n = z_floor - z and vn = v_z exactly (sloped = .FALSE.).
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: x, y, z, vel(3)
    REAL(wp), INTENT(OUT) :: gap_n, vn, nvec(3)
    LOGICAL, INTENT(OUT) :: sloped
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es
    REAL(wp) :: z_floor, dfdx, dfdy, sfac
    CHARACTER(160) :: em

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    gap_n = CD_ZERO
    vn = CD_ZERO
    nvec = [CD_ZERO, CD_ZERO, CD_ONE]
    sloped = .FALSE.
    dfdx = CD_ZERO
    dfdy = CD_ZERO
    IF (model%has_bathymetry) THEN
      CALL CD_Bathymetry_Floor_Gradient(model%bathymetry, x, y, z_floor, dfdx, dfdy, es, em)
      IF (es /= 0) THEN
        ErrStat = CD_MODEL_BADINPUT
        ErrMsg = 'CableDyn_Model: bathymetry gradient query failed: '//TRIM(em)
        RETURN
      END IF
    ELSE
      z_floor = model%seabed_z_floor
    END IF
    sloped = ABS(dfdx) > CD_ZERO .OR. ABS(dfdy) > CD_ZERO
    IF (sloped) THEN
      sfac = SQRT(CD_ONE + dfdx*dfdx + dfdy*dfdy)
      nvec = [-dfdx, -dfdy, CD_ONE]/sfac
      gap_n = (z_floor - z)/sfac
      vn = DOT_PRODUCT(vel, nvec)
    ELSE
      gap_n = z_floor - z
      vn = vel(3)
    END IF
  END SUBROUTINE node_contact_normal

  REAL(wp) FUNCTION seabed_normal_reaction(model, x, y, z, vel, node, ErrStat, ErrMsg) RESULT(normal)
    !! Magnitude of the seabed normal reaction (elastic law plus one-sided damping) at the
    !! normal penetration and normal velocity; the friction capacity is mu times it.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: x, y, z, vel(3)
    INTEGER, INTENT(IN) :: node
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: gap_n, vn, nvec(3), tangent, damping, dgap, dvz
    LOGICAL :: sloped

    normal = CD_ZERO
    CALL node_contact_normal(model, x, y, z, vel, gap_n, vn, nvec, sloped, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL CD_Seabed_Normal_Law(gap_n, model%seabed_kn(node), normal, tangent)
    IF (model%has_seabed_damping) THEN
      CALL seabed_normal_damping_terms(gap_n, vn, model%seabed_cn(node), damping, dgap, dvz)
      normal = normal + damping
    END IF
  END FUNCTION seabed_normal_reaction

  SUBROUTINE seabed_normal_position_gradient(model, x, y, z, vel, node, dn_dq, ErrStat, ErrMsg)
    !! d(normal reaction)/dq with the floor gradient (and so the normal) held over the node.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: x, y, z, vel(3)
    INTEGER, INTENT(IN) :: node
    REAL(wp), INTENT(OUT) :: dn_dq(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: gap_n, vn, nvec(3), normal, tangent, damping, damp_dgap, damp_dvz
    LOGICAL :: sloped

    dn_dq = CD_ZERO
    CALL node_contact_normal(model, x, y, z, vel, gap_n, vn, nvec, sloped, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL CD_Seabed_Normal_Law(gap_n, model%seabed_kn(node), normal, tangent)
    IF (model%has_seabed_damping) THEN
      CALL seabed_normal_damping_terms(gap_n, vn, model%seabed_cn(node), damping, damp_dgap, damp_dvz)
      tangent = tangent + damp_dgap
    END IF
    ! d(gap_n)/dq = -n
    dn_dq = -tangent*nvec
  END SUBROUTINE seabed_normal_position_gradient

  SUBROUTINE seabed_normal_velocity_gradient(model, x, y, z, vel, node, dn_dv, ErrStat, ErrMsg)
    !! d(normal reaction)/dv = (d damping/d vn) n.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: x, y, z, vel(3)
    INTEGER, INTENT(IN) :: node
    REAL(wp), INTENT(OUT) :: dn_dv(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: gap_n, vn, nvec(3), damping, damp_dgap, dn_dvn
    LOGICAL :: sloped

    dn_dv = CD_ZERO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (.NOT. model%has_seabed_damping) RETURN
    CALL node_contact_normal(model, x, y, z, vel, gap_n, vn, nvec, sloped, ErrStat, ErrMsg)
    IF (ErrStat /= CD_MODEL_OK) RETURN
    CALL seabed_normal_damping_terms(gap_n, vn, model%seabed_cn(node), damping, damp_dgap, dn_dvn)
    dn_dv = dn_dvn*nvec
  END SUBROUTINE seabed_normal_velocity_gradient

  PURE SUBROUTINE seabed_normal_damping_terms(gap, vz, c_n, force, dforce_dgap, dforce_dvz)
    !! One-sided normal seabed damping at normal penetration gap and normal velocity vz
    !! (vertical on a level floor): force = -c_n*gate(gap)*vz for vz < 0.
    REAL(wp), INTENT(IN) :: gap, vz, c_n
    REAL(wp), INTENT(OUT) :: force, dforce_dgap, dforce_dvz
    REAL(wp) :: gate, dgate_dgap

    force = CD_ZERO
    dforce_dgap = CD_ZERO
    dforce_dvz = CD_ZERO
    IF (c_n <= CD_ZERO .OR. vz >= CD_ZERO) RETURN
    CALL CD_Seabed_Contact_Gate(gap, gate, dgate_dgap)
    force = -c_n*gate*vz
    dforce_dgap = -c_n*dgate_dgap*vz
    dforce_dvz = -c_n*gate
  END SUBROUTINE seabed_normal_damping_terms

  SUBROUTINE model_morison_drag_load(model, q, v, force, jac_q, jac_v, ErrStat, ErrMsg)
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: e, conn1(2, 1), idx(6)
    REAL(wp) :: q2(6), v2(6), l02(1), fluid2(3, 2), wl2(2), f2(6), jq2(6, 6), jv2(6, 6)

    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    conn1(:, 1) = [1, 2]
    DO e = 1, model%n_elem
      CALL extract_element_state(model, e, q, v, idx, q2, v2)
      CALL extract_nodal_vector(model%fluid_velocity, model, e, fluid2)
      CALL extract_nodal_scalar(model%drag_waterline_z, model, e, wl2)
      l02(1) = model%l0(e)
      CALL CD_Cable_Morison_Drag_Load(q2, v2, conn1, l02, fluid2, wl2, model%drag_rho, &
                                      model%drag_diameter_elem(e), model%drag_cdn_elem(e), &
                                      model%drag_cdt_elem(e), f2, jq2, jv2, ErrStat, ErrMsg)
      IF (ErrStat /= 0) RETURN
      CALL scatter_element_vector(force, idx, f2)
      CALL scatter_element_matrix(jac_q, idx, jq2)
      CALL scatter_element_matrix(jac_v, idx, jv2)
    END DO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END SUBROUTINE model_morison_drag_load

  SUBROUTINE model_froude_krylov_load(model, q, v, force, jac_q, jac_v, ErrStat, ErrMsg)
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: e, conn1(2, 1), idx(6)
    REAL(wp) :: q2(6), v2(6), l02(1), accel2(3, 2), wl2(2), f2(6), jq2(6, 6), jv2(6, 6)

    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    conn1(:, 1) = [1, 2]
    DO e = 1, model%n_elem
      CALL extract_element_state(model, e, q, v, idx, q2, v2)
      CALL extract_nodal_vector(model%fluid_acceleration, model, e, accel2)
      CALL extract_nodal_scalar(model%fk_waterline_z, model, e, wl2)
      l02(1) = model%l0(e)
      CALL CD_Cable_Froude_Krylov_Load(q2, v2, conn1, l02, accel2, wl2, model%fk_rho, &
                                       model%fk_diameter_elem(e), model%fk_can_elem(e), &
                                       model%fk_cat_elem(e), f2, jq2, jv2, ErrStat, ErrMsg)
      IF (ErrStat /= 0) RETURN
      CALL scatter_element_vector(force, idx, f2)
      CALL scatter_element_matrix(jac_q, idx, jq2)
      CALL scatter_element_matrix(jac_v, idx, jv2)
    END DO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END SUBROUTINE model_froude_krylov_load

  SUBROUTINE model_buoyancy_recovery_load(model, q, v, force, jac_q, jac_v, ErrStat, ErrMsg)
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: e, conn1(2, 1), idx(6)
    REAL(wp) :: q2(6), v2(6), l02(1), wl2(2), f2(6), jq2(6, 6), jv2(6, 6)

    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    conn1(:, 1) = [1, 2]
    DO e = 1, model%n_elem
      CALL extract_element_state(model, e, q, v, idx, q2, v2)
      CALL extract_nodal_scalar(model%buoyancy_waterline_z, model, e, wl2)
      l02(1) = model%l0(e)
      CALL CD_Cable_Buoyancy_Recovery_Load(q2, v2, conn1, l02, wl2, model%buoyancy_rho, &
                                           model%buoyancy_diameter_elem(e), model%buoyancy_gravity, &
                                           f2, jq2, jv2, ErrStat, ErrMsg)
      IF (ErrStat /= 0) RETURN
      CALL scatter_element_vector(force, idx, f2)
      CALL scatter_element_matrix(jac_q, idx, jq2)
      CALL scatter_element_matrix(jac_v, idx, jv2)
    END DO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END SUBROUTINE model_buoyancy_recovery_load

  SUBROUTINE model_added_mass_full(model, q, accel, M_add, dMa_a_dq, ErrStat, ErrMsg)
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), accel(:)
    REAL(wp), INTENT(OUT) :: M_add(:, :), dMa_a_dq(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: e, conn1(2, 1), idx(6)
    REAL(wp) :: q2(6), a2(6), l02(1), wl2(2), m2(6, 6), dm2(6, 6)

    M_add = CD_ZERO
    dMa_a_dq = CD_ZERO
    conn1(:, 1) = [1, 2]
    DO e = 1, model%n_elem
      CALL extract_element_state(model, e, q, accel, idx, q2, a2)
      CALL extract_nodal_scalar(model%added_mass_waterline_z, model, e, wl2)
      l02(1) = model%l0(e)
      CALL CD_Cable_Added_Mass(q2, a2, conn1, l02, wl2, model%added_mass_rho, &
                               model%added_mass_diameter_elem(e), model%added_mass_can_elem(e), &
                               model%added_mass_cat_elem(e), m2, dm2, ErrStat, ErrMsg)
      IF (ErrStat /= 0) RETURN
      CALL scatter_element_matrix(M_add, idx, m2)
      CALL scatter_element_matrix(dMa_a_dq, idx, dm2)
    END DO
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END SUBROUTINE model_added_mass_full

  SUBROUTINE model_added_mass_banded(model, q, accel, free, kl, ku, M_add_a, M_add_band, dMa_a_dq_band, &
                                     ErrStat, ErrMsg, need_tangent)
    TYPE(CD_ModelType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: q(:), accel(:)
    INTEGER, INTENT(IN) :: free(:), kl, ku
    REAL(wp), INTENT(OUT) :: M_add_a(:), M_add_band(:, :), dMa_a_dq_band(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    LOGICAL, INTENT(IN), OPTIONAL :: need_tangent

    INTEGER :: e, i, n, ldab, idx(6)
    REAL(wp) :: q2(6), a2(6), wl2(2), m2(6, 6), dm2(6, 6), ma2(6)
    LOGICAL :: want_tangent, zero_accel, free_map_ready

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    want_tangent = .TRUE.
    IF (PRESENT(need_tangent)) want_tangent = need_tangent
    zero_accel = MAXVAL(ABS(accel)) <= EPSILON(CD_ZERO)
    M_add_a = CD_ZERO
    IF (want_tangent) THEN
      M_add_band = CD_ZERO
      dMa_a_dq_band = CD_ZERO
    END IF
    n = SIZE(q)
    ldab = 2*kl + ku + 1
    IF (SIZE(accel) /= n .OR. SIZE(M_add_a) /= n) THEN
      CALL fail(ErrStat, ErrMsg, 'banded added mass received inconsistent q/accel/output shapes')
      RETURN
    END IF
    IF (want_tangent .AND. &
        (kl < 0 .OR. ku < 0 .OR. SIZE(M_add_band, 1) < ldab .OR. SIZE(dMa_a_dq_band, 1) < ldab .OR. &
         SIZE(M_add_band, 2) /= SIZE(free) .OR. SIZE(dMa_a_dq_band, 2) /= SIZE(free))) THEN
      CALL fail(ErrStat, ErrMsg, 'banded added mass received invalid band storage')
      RETURN
    END IF
    free_map_ready = ALLOCATED(model%free_map_work)
    IF (free_map_ready) free_map_ready = SIZE(model%free_map_work) >= n
    IF (want_tangent .AND. .NOT. free_map_ready) THEN
      CALL fail(ErrStat, ErrMsg, 'banded added mass requires initialized free-map workspace')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(q) .OR. .NOT. CD_All_Finite(accel)) THEN
      CALL fail(ErrStat, ErrMsg, 'banded added mass q/accel must be finite')
      RETURN
    END IF

    IF (want_tangent) THEN
      model%free_map_work(1:n) = 0
      DO i = 1, SIZE(free)
        IF (free(i) < 1 .OR. free(i) > n) THEN
          CALL fail(ErrStat, ErrMsg, 'banded added mass free DOF index out of range')
          RETURN
        END IF
        IF (model%free_map_work(free(i)) /= 0) THEN
          CALL fail(ErrStat, ErrMsg, 'banded added mass duplicate free DOF')
          RETURN
        END IF
        model%free_map_work(free(i)) = i
      END DO
    END IF

    DO e = 1, model%n_elem
      CALL extract_element_state(model, e, q, accel, idx, q2, a2)
      CALL extract_nodal_scalar(model%added_mass_waterline_z, model, e, wl2)
      IF (want_tangent) THEN
        IF (zero_accel) THEN
          CALL CD_Added_Mass_Element_Matrix(q2(1:3), q2(4:6), model%l0(e), wl2(1), wl2(2), &
                                            model%added_mass_rho, model%added_mass_diameter_elem(e), &
                                            model%added_mass_can_elem(e), model%added_mass_cat_elem(e), &
                                            m2, ErrStat, ErrMsg)
          IF (ErrStat /= 0) RETURN
          ma2 = CD_ZERO
        ELSE
          CALL CD_Added_Mass_Element(q2(1:3), q2(4:6), a2(1:3), a2(4:6), model%l0(e), wl2(1), wl2(2), &
                                     model%added_mass_rho, model%added_mass_diameter_elem(e), &
                                     model%added_mass_can_elem(e), model%added_mass_cat_elem(e), &
                                     m2, dm2, ErrStat, ErrMsg)
          IF (ErrStat /= 0) RETURN
          ma2 = MATMUL(m2, a2)
        END IF
      ELSE
        CALL CD_Added_Mass_Element_Force(q2(1:3), q2(4:6), a2(1:3), a2(4:6), model%l0(e), wl2(1), wl2(2), &
                                         model%added_mass_rho, model%added_mass_diameter_elem(e), &
                                         model%added_mass_can_elem(e), model%added_mass_cat_elem(e), &
                                         ma2, ErrStat, ErrMsg)
        IF (ErrStat /= 0) RETURN
      END IF
      CALL scatter_element_vector(M_add_a, idx, ma2)
      IF (want_tangent) THEN
        CALL scatter_element_matrix_band(M_add_band, model%free_map_work, kl, ku, idx, m2, ErrStat, ErrMsg)
        IF (ErrStat /= CD_MODEL_OK) RETURN
        IF (.NOT. zero_accel) THEN
          CALL scatter_element_matrix_band(dMa_a_dq_band, model%free_map_work, kl, ku, idx, dm2, ErrStat, ErrMsg)
          IF (ErrStat /= CD_MODEL_OK) RETURN
        END IF
      END IF
    END DO
    IF (.NOT. CD_All_Finite(M_add_a) .OR. &
        (want_tangent .AND. (.NOT. CD_All_Finite(M_add_band) .OR. &
                             .NOT. CD_All_Finite(dMa_a_dq_band)))) THEN
      CALL fail(ErrStat, ErrMsg, 'banded added mass returned non-finite force/matrix/Jacobian')
    END IF
  END SUBROUTINE model_added_mass_banded

  PURE SUBROUTINE model_structural_mass_product(model, x, product)
    !! Apply the constant consistent translational mass in O(n_elem) without forming
    !! or traversing its dense storage. Each element is m/6 [[2I,I],[I,2I]].
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: x(:)
    REAL(wp), INTENT(OUT) :: product(:)

    INTEGER :: e, a, b, ia, ib
    REAL(wp) :: coeff

    product = CD_ZERO
    DO e = 1, model%n_elem
      a = model%elem_conn(1, e); b = model%elem_conn(2, e)
      ia = 3*a - 2; ib = 3*b - 2
      coeff = model%rho_a(e)*model%l0(e)/6.0_wp
      product(ia:ia + 2) = product(ia:ia + 2) + coeff*(2.0_wp*x(ia:ia + 2) + x(ib:ib + 2))
      product(ib:ib + 2) = product(ib:ib + 2) + coeff*(x(ia:ia + 2) + 2.0_wp*x(ib:ib + 2))
    END DO
  END SUBROUTINE model_structural_mass_product

  SUBROUTINE model_added_mass_force(model, q, accel, M_add_a, ErrStat, ErrMsg)
    !! Force-only added-mass application for load recovery. This avoids allocating,
    !! zeroing and multiplying an n_dof-square matrix on every CalcOutput call.
    TYPE(CD_ModelType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q(:), accel(:)
    REAL(wp), INTENT(OUT) :: M_add_a(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: e, idx(6)
    REAL(wp) :: q2(6), a2(6), wl2(2), ma2(6)

    M_add_a = CD_ZERO
    DO e = 1, model%n_elem
      CALL extract_element_positions(model, e, q, idx, q2)
      a2 = accel(idx)
      CALL extract_nodal_scalar(model%added_mass_waterline_z, model, e, wl2)
      CALL CD_Added_Mass_Element_Force(q2(1:3), q2(4:6), a2(1:3), a2(4:6), model%l0(e), wl2(1), wl2(2), &
                                       model%added_mass_rho, model%added_mass_diameter_elem(e), &
                                       model%added_mass_can_elem(e), model%added_mass_cat_elem(e), &
                                       ma2, ErrStat, ErrMsg)
      IF (ErrStat /= 0) RETURN
      CALL scatter_element_vector(M_add_a, idx, ma2)
    END DO
    IF (.NOT. CD_All_Finite(M_add_a)) THEN
      CALL fail(ErrStat, ErrMsg, 'added-mass force returned non-finite values')
      RETURN
    END IF
    ErrStat = CD_MODEL_OK
    ErrMsg = ''
  END SUBROUTINE model_added_mass_force

  PURE SUBROUTINE extract_element_state(model, e, q, v, idx, q2, v2)
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: e
    REAL(wp), INTENT(IN) :: q(:), v(:)
    INTEGER, INTENT(OUT) :: idx(6)
    REAL(wp), INTENT(OUT) :: q2(6), v2(6)

    INTEGER :: a, b

    a = model%elem_conn(1, e)
    b = model%elem_conn(2, e)
    idx = [3*a - 2, 3*a - 1, 3*a, 3*b - 2, 3*b - 1, 3*b]
    q2 = q(idx)
    v2 = v(idx)
  END SUBROUTINE extract_element_state

  PURE SUBROUTINE extract_element_positions(model, e, q, idx, q2)
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: e
    REAL(wp), INTENT(IN) :: q(:)
    INTEGER, INTENT(OUT) :: idx(6)
    REAL(wp), INTENT(OUT) :: q2(6)

    INTEGER :: a, b

    a = model%elem_conn(1, e)
    b = model%elem_conn(2, e)
    idx = [3*a - 2, 3*a - 1, 3*a, 3*b - 2, 3*b - 1, 3*b]
    q2 = q(idx)
  END SUBROUTINE extract_element_positions

  PURE SUBROUTINE extract_nodal_vector(field, model, e, field2)
    REAL(wp), INTENT(IN) :: field(:, :)
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: e
    REAL(wp), INTENT(OUT) :: field2(3, 2)

    field2(:, 1) = field(:, model%elem_conn(1, e))
    field2(:, 2) = field(:, model%elem_conn(2, e))
  END SUBROUTINE extract_nodal_vector

  PURE SUBROUTINE extract_nodal_scalar(field, model, e, field2)
    REAL(wp), INTENT(IN) :: field(:)
    TYPE(CD_ModelType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: e
    REAL(wp), INTENT(OUT) :: field2(2)

    field2(1) = field(model%elem_conn(1, e))
    field2(2) = field(model%elem_conn(2, e))
  END SUBROUTINE extract_nodal_scalar

  PURE SUBROUTINE scatter_element_vector(global, idx, local)
    REAL(wp), INTENT(INOUT) :: global(:)
    INTEGER, INTENT(IN) :: idx(6)
    REAL(wp), INTENT(IN) :: local(6)
    INTEGER :: i

    DO i = 1, 6
      global(idx(i)) = global(idx(i)) + local(i)
    END DO
  END SUBROUTINE scatter_element_vector

  PURE SUBROUTINE scatter_element_matrix(global, idx, local)
    REAL(wp), INTENT(INOUT) :: global(:, :)
    INTEGER, INTENT(IN) :: idx(6)
    REAL(wp), INTENT(IN) :: local(6, 6)
    INTEGER :: i, j

    DO j = 1, 6
      DO i = 1, 6
        global(idx(i), idx(j)) = global(idx(i), idx(j)) + local(i, j)
      END DO
    END DO
  END SUBROUTINE scatter_element_matrix

  SUBROUTINE scatter_element_matrix_band(ab, free_map, kl, ku, idx, local, ErrStat, ErrMsg)
    !! Scatter a 6x6 element residual-Jacobian block into reduced free/free
    !! LAPACK band storage. Fixed/prescribed rows or columns are skipped because
    !! Newton only solves the free block. This runs for every contributor of every
    !! element on each Jacobian evaluation, so ErrMsg is written only on failure.
    REAL(wp), INTENT(INOUT) :: ab(:, :)
    INTEGER, INTENT(IN) :: free_map(:), kl, ku, idx(6)
    REAL(wp), INTENT(IN) :: local(6, 6)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(INOUT) :: ErrMsg
    INTEGER :: i, j, row, fmap(6)

    ErrStat = CD_MODEL_OK
    IF (ANY(idx < 1) .OR. ANY(idx > SIZE(free_map))) THEN
      CALL fail(ErrStat, ErrMsg, 'band scatter global index out of range')
      RETURN
    END IF
    fmap = free_map(idx)
    DO j = 1, 6
      IF (fmap(j) == 0) CYCLE
      DO i = 1, 6
        IF (fmap(i) == 0) CYCLE
        row = kl + ku + 1 + fmap(i) - fmap(j)
        IF (row < 1 .OR. row > SIZE(ab, 1)) THEN
          CALL fail(ErrStat, ErrMsg, 'band scatter bandwidth too small')
          RETURN
        END IF
        ab(row, fmap(j)) = ab(row, fmap(j)) + local(i, j)
      END DO
    END DO
  END SUBROUTINE scatter_element_matrix_band

  SUBROUTINE scatter_scalar_band(ab, free_map, kl, ku, gi, gj, value, ErrStat, ErrMsg)
    !! Add one global matrix entry to reduced free/free LAPACK band storage.
    !! ErrMsg is written only on failure (this is a per-node hot path).
    REAL(wp), INTENT(INOUT) :: ab(:, :)
    INTEGER, INTENT(IN) :: free_map(:), kl, ku, gi, gj
    REAL(wp), INTENT(IN) :: value
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(INOUT) :: ErrMsg
    INTEGER :: fi, fj, row

    ErrStat = CD_MODEL_OK
    IF (gi < 1 .OR. gi > SIZE(free_map) .OR. gj < 1 .OR. gj > SIZE(free_map)) THEN
      CALL fail(ErrStat, ErrMsg, 'band scatter global index out of range')
      RETURN
    END IF
    fi = free_map(gi)
    fj = free_map(gj)
    IF (fi == 0 .OR. fj == 0) RETURN
    row = kl + ku + 1 + fi - fj
    IF (row < 1 .OR. row > SIZE(ab, 1)) THEN
      CALL fail(ErrStat, ErrMsg, 'band scatter bandwidth too small')
      RETURN
    END IF
    ab(row, fj) = ab(row, fj) + value
  END SUBROUTINE scatter_scalar_band

  PURE FUNCTION reshape3(q, n_dof) RESULT(nodes)
    !! Copy the flat positions-only state into a (3, n_nodes) node array.
    INTEGER, INTENT(IN) :: n_dof
    REAL(wp), INTENT(IN) :: q(n_dof)
    REAL(wp) :: nodes(3, n_dof/3)
    nodes = RESHAPE(q, [3, n_dof/3])
  END FUNCTION reshape3

  SUBROUTINE validate_vector_field(source, n, name, ErrStat, ErrMsg)
    !! Validate a model-owned nodal scalar field without mutating storage.
    REAL(wp), INTENT(IN) :: source(:)
    INTEGER, INTENT(IN) :: n
    CHARACTER(*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (SIZE(source) /= n) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(name)//' update must have shape (n_nodes)')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(source)) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(name)//' update must be finite')
    END IF
  END SUBROUTINE validate_vector_field

  SUBROUTINE validate_positive_vector(values, n, name, ErrStat, ErrMsg)
    !! Validate a finite positive element-wise coefficient vector.
    REAL(wp), INTENT(IN) :: values(:)
    INTEGER, INTENT(IN) :: n
    CHARACTER(*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (SIZE(values) /= n) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(name)//' must have shape (n_elem)')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(values) .OR. ANY(values <= CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(name)//' must be finite and positive')
    END IF
  END SUBROUTINE validate_positive_vector

  SUBROUTINE validate_nonnegative_vector(values, n, name, ErrStat, ErrMsg)
    !! Validate a finite non-negative element-wise coefficient vector.
    REAL(wp), INTENT(IN) :: values(:)
    INTEGER, INTENT(IN) :: n
    CHARACTER(*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    ErrStat = CD_MODEL_OK
    ErrMsg = ''
    IF (SIZE(values) /= n) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(name)//' must have shape (n_elem)')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(values) .OR. ANY(values < CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, TRIM(name)//' must be finite and non-negative')
    END IF
  END SUBROUTINE validate_nonnegative_vector

  SUBROUTINE fail(ErrStat, ErrMsg, msg)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    CHARACTER(*), INTENT(IN) :: msg
    ErrStat = CD_MODEL_BADINPUT
    ErrMsg = 'CableDyn_Model: '//msg
  END SUBROUTINE fail

  SUBROUTINE alloc_fail(ErrStat, ErrMsg, msg)
    !! Distinct allocation-failure reporter: resource exhaustion is not malformed input, so
    !! these paths return CD_MODEL_ALLOCFAIL and let deck/system/C-ABI callers map it to their
    !! own allocation status rather than collapsing it to bad input or solve failure.
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    CHARACTER(*), INTENT(IN) :: msg
    ErrStat = CD_MODEL_ALLOCFAIL
    ErrMsg = 'CableDyn_Model: '//msg
  END SUBROUTINE alloc_fail

  PURE INTEGER FUNCTION as_dyn_callback_status(model_es) RESULT(dyn_es)
    !! Translate a model load-callback status into the integrator's namespace so an allocation
    !! failure inside model_contributor_load/_force is forwarded as CD_DYN_ALLOCFAIL (and preserved
    !! by the integrator) rather than collapsed to a generic dynamic failure. Other statuses pass
    !! through unchanged (the integrator treats any nonzero callback status as a failure).
    INTEGER, INTENT(IN) :: model_es
    IF (model_es == CD_MODEL_ALLOCFAIL) THEN
      dyn_es = CD_DYN_ALLOCFAIL
    ELSE
      dyn_es = model_es
    END IF
  END FUNCTION as_dyn_callback_status

  PURE INTEGER FUNCTION model_status_from_dyn(dyn_es) RESULT(model_es)
    !! Map an integrator status back into the model namespace at the dynamic-step / acceleration
    !! wrap sites, preserving a propagated allocation failure as CD_MODEL_ALLOCFAIL.
    INTEGER, INTENT(IN) :: dyn_es
    IF (dyn_es == CD_DYN_ALLOCFAIL) THEN
      model_es = CD_MODEL_ALLOCFAIL
    ELSE
      model_es = CD_MODEL_SOLVEFAIL
    END IF
  END FUNCTION model_status_from_dyn

END MODULE CableDyn_Model
