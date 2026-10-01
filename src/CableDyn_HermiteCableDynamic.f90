! File: src/CableDyn_HermiteCableDynamic.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_HermiteCableDynamic
  !! Finite-EI DYNAMICS of a cubic-Hermite bending-cable line (CableDyn_HermiteCable).
  !!
  !! The line is a chain of 2-node elements sharing interior nodes; each node carries
  !! [r(3), m(3)] with m = dr/ds, so the global state is 6*n_nodes position-like DOFs. Because
  !! there are NO rotation DOFs the equations of motion are the STANDARD nonlinear structural
  !! dynamics system
  !!   M q_ddot + f_int(q) = f_ext(q, t),
  !! with a CONSTANT consistent mass M (the reference-length cubic-Hermite mass), the closed-form
  !! internal force/tangent from CD_HermiteCable_Element, distributed self-weight, and a penalty
  !! seabed reaction. None of the Cosserat finite-EI complications arise: no |theta| < pi chart,
  !! no additive-vs-spatial angular-velocity bookkeeping, no gyroscopic term, no SO(3) tangent
  !! Jacobians -- q_ddot and q live in the same flat vector space, so the generalised-alpha step
  !! is the textbook Chung & Hulbert (1993) form.
  !!
  !! This is the finite-EI power-cable / lazy-wave DYNAMIC path: it inherits the wall-free Hermite
  !! element that the static solver (CableDyn_HermiteCableStatic) uses, so it can carry the same
  !! buoyant-arch geometries the Cosserat dynamic model cannot. The static solve provides the
  !! at-rest equilibrium seed; this module advances it in time under gravity, seabed contact,
  !! inertia, and optional Morison quadratic drag (CD_HermiteCable_Dyn_Set_Drag) -- the damping that
  !! makes a driven dynamic cable physical, consistently distributed over the Hermite shape functions
  !! (CD_HermiteCable_Drag_Element) so the material-tangent DOFs carry drag too. Boundary DOFs listed
  !! in fixed_dofs are held at rest by default, or driven by a prescribed time history
  !! (a moving fairlead / hang-off heave) via the
  !! optional pres_* arguments to CD_HermiteCable_Dyn_Step. The Newton solve starts from the
  !! Newmark predictor (or, when the motion is smooth, a constant-acceleration predictor) and is
  !! globalised with a backtracking line search so it stays robust at the large implicit dt the
  !! stiff (EA >> EI) cable demands; the velocity-dependent drag Jacobian is folded into the effective
  !! tangent through the Newmark dv/dq coupling. By default (CD_HC_FORCE_BLEND_DEFAULT) the
  !! residual blends the forces, G = M a_am + (1-af) r(q_{n+1}) + af r(q_n)
  !! (CD_HermiteCable_Dyn_Set_ForceBlend), rather than evaluating them at the blended state as
  !! written below, and the factored tangent of a committed step can start the next one
  !! (cross-step reuse, CD_HermiteCable_Dyn_Set_TangentReuse); the constant-acceleration
  !! predictor and the reuse apply to the force blend only. A regular (Airy) wave field
  !! (CD_HermiteCable_Dyn_Set_Waves) feeds the drag's fluid velocity and drives the Froude-Krylov +
  !! fluid-inertia load, evaluated at the generalised-alpha intermediate time each step.
  !!
  !! Generalised-alpha (Chung & Hulbert 1993), spectral radius rho_inf in [0,1]:
  !!   alpha_m = (2 rho_inf - 1)/(rho_inf + 1),  alpha_f = rho_inf/(rho_inf + 1),
  !!   beta = 1/4 (1 - alpha_m + alpha_f)^2,     gamma = 1/2 - alpha_m + alpha_f.
  !! Newmark update in the unknown q_{n+1}:
  !!   a_{n+1} = (q_{n+1}-q_n)/(beta dt^2) - v_n/(beta dt) - (1/(2 beta) - 1) a_n,
  !!   v_{n+1} = v_n + dt((1-gamma) a_n + gamma a_{n+1}).
  !! Residual at the generalised-alpha intermediate point (f_ext may depend on velocity via drag):
  !!   G(q_{n+1}) = M((1-alpha_m) a_{n+1} + alpha_m a_n)
  !!             + f_int((1-alpha_f) q_{n+1} + alpha_f q_n) - f_ext(q_alpha, v_alpha),
  !! effective tangent  dG/dq_{n+1} = (1-alpha_m)/(beta dt^2) M + (1-alpha_f) K(q_alpha)
  !!                                + (1-alpha_f) gamma/(beta dt) Kv(q_alpha, v_alpha),
  !! where Kv = dR/dv is the drag velocity Jacobian and gamma/(beta dt) = d v_{n+1}/d q_{n+1}.
  USE CableDyn_Precision, ONLY: wp, CD_ZERO, CD_ONE, CD_All_Finite, CD_Is_Finite
  USE CableDyn_EndConnection, ONLY: CD_EndConn_Spring, CD_EndConn_Energy, &
                                    CD_EndConn_Reaction_Moment, CD_EndConn_Basis, &
                                    CD_EndConn_Project, CD_ENDCONN_OK, &
                                    CD_ENDCONN_PINNED, CD_ENDCONN_FINITE, CD_ENDCONN_RIGID
  USE CableDyn_HermiteCable, ONLY: CD_HermiteCable_Element, CD_HermiteCable_Mass, &
                                   CD_HermiteCable_Curvature, CD_HermiteCable_Shapes, &
                                   CD_HermiteCable_Gauss_Rule, &
                                   CD_HermiteCable_Axial_Resultant, &
                                   CD_HermiteCable_Axial_Strain_Lower_Bound, CD_HermiteCable_Dry_Buoyancy, &
                                   CD_HCABLE_OK
  USE CableDyn_Linalg, ONLY: CD_Factor_Banded, CD_Solve_Factored_Banded, CD_LINALG_OK
  USE CableDyn_Hydro, ONLY: CD_Morison_Drag_Per_Length_Jac, CD_Airy_Wave_Kinematics_Precomputed, &
                            CD_Component_Wave_Kinematics, CD_Solve_Dispersion_Wavenumber, CD_HYDRO_OK
  USE CableDyn_SeabedContact, ONLY: CD_Seabed_Normal_Law, CD_Seabed_Contact_Gate, CD_SEABED_CONTACT_BLEND, &
                                    CD_Seabed_Friction_Stick_Slip, CD_Seabed_Friction_Aniso, &
                                    CD_Seabed_Friction_Mu_Dir, CD_FRICTION_STICK_SLIP
  USE CableDyn_Bathymetry, ONLY: CD_BathymetryType, CD_Bathymetry_Is_Initialized, &
                                 CD_Bathymetry_Floor_Gradient, CD_End_Bathymetry, CD_BATHY_OK
  USE CableDyn_HermiteTorsion, ONLY: CD_HermiteTorsionType, CD_HermiteTorsion_Line, CD_HermiteTorsion_Unwrap, &
                                     CD_HermiteTorsion_Compliance, CD_HermiteTorsion_Validate, CD_HTORS_OK, &
                                     CD_HTORS_MAX_STEP
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64
!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_in_parallel
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CD_HermiteCableDynType
  PUBLIC :: CD_HermiteCable_Dyn_Init
  PUBLIC :: CD_HermiteCable_Dyn_Set_Drag
  PUBLIC :: CD_HermiteCable_Dyn_Set_Axial_Damping
  PUBLIC :: CD_HermiteCable_Dyn_Set_EndConnection
  PUBLIC :: CD_HermiteCable_Dyn_Set_Torsion
  PUBLIC :: CD_HermiteCable_Dyn_Torsion_State
  PUBLIC :: CD_HermiteCable_Dyn_EndConnection_Moment
  PUBLIC :: CD_HermiteCable_Dyn_Set_ModifiedNewton
  PUBLIC :: CD_HermiteCable_Dyn_Set_TangentReuse
  PUBLIC :: CD_HermiteCable_Dyn_Set_AdaptiveNewton
  PUBLIC :: CD_HermiteCable_Dyn_Set_Tensile_Safety
  PUBLIC :: CD_HermiteCable_Dyn_Set_Tensile_Monitor
  PUBLIC :: CD_HermiteCable_Dyn_Get_Tensile_Diagnostics
  PUBLIC :: CD_HermiteCable_Dyn_Get_Recovery_Diagnostics
  PUBLIC :: CD_HermiteCable_Dyn_Set_AddedMass
  PUBLIC :: CD_HermiteCable_Dyn_Set_Waves
  PUBLIC :: CD_HermiteCable_Dyn_Set_Held_Fluid
  PUBLIC :: CD_HermiteCable_Dyn_Set_Contact
  PUBLIC :: CD_HermiteCable_Dyn_Set_Attachments
  PUBLIC :: CD_HermiteCable_Attachment_Drag
  PUBLIC :: CD_HermiteCable_Dyn_Recompute_Acceleration
  PUBLIC :: CD_HermiteCable_Dyn_Set_Irregular_Waves
  PUBLIC :: CD_HermiteCable_Drag_Element
  PUBLIC :: CD_HermiteCable_Axial_Damping_Element
  PUBLIC :: CD_HermiteCable_Axial_Damping_Resultant
  PUBLIC :: CD_HermiteCable_AddedMass_Element
  PUBLIC :: CD_HermiteCable_FK_Element
  PUBLIC :: CD_HermiteCable_Dyn_Step
  PUBLIC :: CD_HermiteCable_Dyn_Step_Recovering
  PUBLIC :: CD_HermiteCable_Dyn_Recovery_Count, CD_HermiteCable_Dyn_Recovery_Reset
  PUBLIC :: CD_HermiteCable_Dyn_Snapshot
  PUBLIC :: CD_HermiteCable_Dyn_Set_Friction_Anchors
  PUBLIC :: CD_HermiteCable_Dyn_Set_Friction_Axial
  PUBLIC :: CD_HermiteCable_Dyn_Restore
  PUBLIC :: CD_HermiteCable_Dyn_Energy
  PUBLIC :: CD_HermiteCable_Dyn_Curvature
  PUBLIC :: CD_HermiteCable_Dyn_Reaction
  PUBLIC :: CD_HermiteCable_Dyn_End_Force
  PUBLIC :: CD_HermiteCable_Dyn_Set_ForceBlend
  PUBLIC :: CD_HermiteCable_Dyn_Max_Step_Rotation
  ! Default generalised-alpha blend of the Hermite dynamics: force blend.
  LOGICAL, PARAMETER, PUBLIC :: CD_HC_FORCE_BLEND_DEFAULT = .TRUE.
  PUBLIC :: CD_HermiteCable_Dyn_End
  PUBLIC :: CD_HermiteCable_Dyn_Reset_Profile, CD_HermiteCable_Dyn_Disable_Profile, &
            CD_HermiteCable_Dyn_Get_Profile
  PUBLIC :: CD_HermiteCable_Dyn_Profile_Enabled
  INTEGER, PARAMETER, PUBLIC :: CD_HCDYN_OK = 0, CD_HCDYN_BADINPUT = 1, CD_HCDYN_NOCONVERGE = 2
  INTEGER, PARAMETER, PUBLIC :: CD_HCDYN_TENSILE_OFF = 0
  INTEGER, PARAMETER, PUBLIC :: CD_HCDYN_TENSILE_WARN = 1
  INTEGER, PARAMETER, PUBLIC :: CD_HCDYN_TENSILE_ERROR = 2
  REAL(wp), PARAMETER :: CD_HCDYN_AXIAL_STRAIN_TOL = 2.0e-6_wp

  !! Diagnostic step-profile counters (performance instrumentation for the banded
  !! dynamic step): wall time and call counts attributed to the dynamic step's
  !! internal phases since the last reset, so a benchmark can split a run into
  !! element/hydro assembly (eval_res), dense effective-tangent formation + free-DOF
  !! reduction, the linear solve, and the frozen added-mass assembly, and can count
  !! the arrays allocated on the step path (zero unless a recovery rewind ran). Same contract as
  !! the Cosserat assembly counters: PROCESS-GLOBAL, single-run, single-threaded
  !! instrumentation -- not per-model and not synchronized; diagnostics only, never
  !! part of a numerical result.
  ! The line is a chain of 2-node elements with 6 DOF/node in nodal order, so every global
  ! operator (mass, stiffness, drag/FK Jacobians, effective tangent) couples DOFs at most 11
  ! apart: a constant half-bandwidth known a priori. All global matrices live in LAPACK
  ! general-band storage (DGBSV layout: entry (i, j) at ab(KL_D + KU_D + 1 + i - j, j), with
  ! KL_D extra top rows for the factorization fill) -- no dense matrices, no bandwidth scan.
  INTEGER, PARAMETER :: KL_D = 11, KU_D = 11
  INTEGER, PARAMETER :: LDAB_D = 2*KL_D + KU_D + 1

  ! Profiling is OFF by default: SYSTEM_CLOCK in the step/residual hot path is measurable
  ! overhead in a production run. CD_HermiteCable_Dyn_Reset_Profile ARMS it (a benchmark
  ! calls it before a timed run); production paths never call it, so they pay one branch.
  LOGICAL, SAVE :: prof_enabled = .FALSE.
  INTEGER, SAVE :: prof_n_step = 0        ! CD_HermiteCable_Dyn_Step calls
  INTEGER, SAVE :: prof_n_resid = 0       ! eval_res calls (predictor + accepted iterates + line-search trials)
  INTEGER, SAVE :: prof_n_tan = 0         ! tangent assemblies (once per continuing solve iteration)
  INTEGER, SAVE :: prof_n_solve = 0       ! effective-tangent linear solves (= Newton iterations)
  INTEGER, SAVE :: prof_n_am = 0          ! frozen added-mass assemblies (once per step when enabled)
  ! Arrays heap-allocated on the step path: the plain Step allocates none; the recovery
  ! rewind of CD_HermiteCable_Dyn_Step_Recovering allocates its interval buffers.
  INTEGER, SAVE :: prof_n_alloc = 0
  REAL(wp), SAVE :: prof_t_step = 0.0_wp  ! total Step wall (s)
  REAL(wp), SAVE :: prof_t_resid = 0.0_wp ! eval_res wall: element/drag/FK/contact assembly + inertia band matvec
  REAL(wp), SAVE :: prof_t_tan = 0.0_wp   ! tangent-assembly wall (element Hessian + drag/FK Jacobians + band scatters)
  REAL(wp), SAVE :: prof_t_eff = 0.0_wp   ! banded Keff triad + in-band Dirichlet wall
  REAL(wp), SAVE :: prof_t_solve = 0.0_wp ! banded linear-solve wall
  REAL(wp), SAVE :: prof_t_am = 0.0_wp    ! frozen added-mass assembly wall (band)

  ! Number of coupling intervals completed by internal generalised-alpha substeps after the nominal
  ! interval was too large for the prescribed boundary motion or environmental load. This counter is
  ! independent of prof_enabled so temporal-refinement studies can report the actual line advances.
  INTEGER, SAVE :: n_substep_recoveries = 0

  TYPE :: CD_HermiteCableDynType
    !! Persistent finite-EI Hermite dynamic model: geometry, gen-alpha constants, and state.
    INTEGER :: ne = 0, nn = 0, ndof = 0
    INTEGER :: axial_quadrature_order = 4, bending_quadrature_order = 4
    REAL(wp), ALLOCATABLE :: l0(:), EA(:), EI(:), rho_a(:), w(:)   ! per element
    ! Optional axial Kelvin-Voigt damping. BA is the axial constitutive dashpot
    ! coefficient [N s] after resolving the deck's direct-BA / negative-zeta convention.
    ! It is integrated over reference arc length with the same Hermite kinematics and
    ! axial quadrature as the elastic resultant. Off by default.
    LOGICAL :: has_axial_damping = .FALSE.
    REAL(wp), ALLOCATABLE :: BA(:)
    REAL(wp) :: seabed_z = CD_ZERO, kn = CD_ZERO
    REAL(wp), ALLOCATABLE :: trib(:)                              ! nodal tributary length
    ! Production standalone contact. contact_kn/contact_cn are nodal coefficients [N/m]
    ! and [N.s/m] of the deck's nodal seabed law. When absent, the scalar kn*tributary
    ! flat-bed term of CD_HermiteCable_Dyn_Init applies.
    LOGICAL :: has_contact = .FALSE., has_contact_bathymetry = .FALSE.
    REAL(wp) :: contact_mu = CD_ZERO, contact_frame_c = CD_ONE, contact_frame_s = CD_ZERO
    ! Anisotropic friction (CD_HermiteCable_Dyn_Set_Friction_Axial): contact_mu is then the
    ! lateral coefficient and contact_mu_axial the one along the nodal tangent.
    LOGICAL :: contact_fr_aniso = .FALSE.
    REAL(wp) :: contact_mu_axial = CD_ZERO
    REAL(wp), ALLOCATABLE :: contact_kn(:), contact_cn(:)
    TYPE(CD_BathymetryType) :: contact_bathymetry
    ! Stick-slip seabed friction (contact with mu > 0). Each node holds a horizontal spring
    ! of stiffness fr_k (the nodal normal stiffness) to its committed anchor
    ! fr_anchor(:, node); the spring force is capped at mu times the normal contact force
    ! and the anchor follows the node when the cap is reached (return mapping at the
    ! committed state). Set_Contact anchors the springs at the present positions;
    ! CD_HermiteCable_Dyn_Set_Friction_Anchors installs the force of a static solve with
    ! friction, which a spring then holds at rest. ws_recovery_fr is the step-start copy of
    ! the recovering step, snap_fr_anchor the aggregate snapshot.
    LOGICAL :: fr_active = .FALSE.
    REAL(wp), ALLOCATABLE :: fr_anchor(:, :), fr_k(:), ws_recovery_fr(:, :), snap_fr_anchor(:, :)
    REAL(wp), ALLOCATABLE :: Mgb(:, :)                            ! constant global mass (band)
    LOGICAL, ALLOCATABLE :: freemask(:), solve_mask(:)
    INTEGER, ALLOCATABLE :: fmap(:)                              ! free -> global DOF index
    REAL(wp), ALLOCATABLE :: q(:), v(:), a(:)                    ! dynamic state
    ! Step workspace, allocated ONCE at Init so a dynamic step performs ZERO heap
    ! allocations: the Newmark/gen-alpha state vectors, the residual/direction vectors,
    ! and the banded global operators (position tangent K, velocity tangent Kv, their
    ! line-search trial copies, the effective tangent, the total and added mass).
    REAL(wp), ALLOCATABLE :: ws_qn(:), ws_vn(:), ws_an(:), ws_qk(:), ws_ak(:)
    REAL(wp), ALLOCATABLE :: ws_recovery_q0(:), ws_recovery_v0(:), ws_recovery_a0(:)
    REAL(wp), ALLOCATABLE :: ws_qa(:), ws_aa(:), ws_va(:), ws_vk(:)
    REAL(wp), ALLOCATABLE :: ws_R(:), ws_G(:), ws_Gt(:), ws_inert(:), ws_qtrial(:), ws_dq(:)
    ! Force-blended generalised-alpha: the unblended residual R(q_{n+1}, v_{n+1}, t_{n+1}) of
    ! the current iterate (ws_Rk) and of a line-search trial (ws_Rt), and the force at the
    ! committed state F_n = R(q_n, v_n, t_n). F_n is cached from the previous step's converged
    ! evaluation and keyed on the exact state it was evaluated at (fc_q, fc_v, fc_t, fc_d0)
    ! plus fc_valid, which every load-configuration setter and Restore clears. A key mismatch
    ! (host-updated boundary, restored or externally written state) recomputes it, so the
    ! cache only ever replaces a bitwise-identical evaluation.
    REAL(wp), ALLOCATABLE :: ws_Rk(:), ws_Rt(:), fc_R(:), fc_q(:), fc_v(:)
    REAL(wp) :: fc_t = 0.0_wp, fc_d0(3, 2) = 0.0_wp
    LOGICAL :: fc_valid = .FALSE.
    ! Generalised-alpha blend: .TRUE. blends the forces (the default), .FALSE. evaluates the
    ! forces at the blended configuration q_alpha = (1-af) q_{n+1} + af q_n.
    LOGICAL :: force_blend = CD_HC_FORCE_BLEND_DEFAULT
    ! Largest turn of a nodal tangent over one committed step [rad] (run diagnostic).
    REAL(wp) :: max_step_rotation = 0.0_wp
    REAL(wp), ALLOCATABLE :: ws_Kb(:, :), ws_Kvb(:, :)
    REAL(wp), ALLOCATABLE :: ws_Keffb(:, :), ws_Msumb(:, :), ws_Mamb(:, :)
    INTEGER, ALLOCATABLE :: ws_ipiv(:)                    ! band-factorization pivots
    REAL(wp), ALLOCATABLE :: ws_elem_f(:, :), ws_elem_K(:, :, :), ws_elem_Kdrag(:, :, :)
    REAL(wp), ALLOCATABLE :: ws_elem_Kv(:, :, :), ws_elem_Kwave(:, :, :), ws_elem_Kheld(:, :, :)
    INTEGER, ALLOCATABLE :: ws_elem_es(:)
    CHARACTER(200), ALLOCATABLE :: ws_elem_em(:)
    REAL(wp) :: alpha_m = CD_ZERO, alpha_f = CD_ZERO, beta = CD_ZERO, gamma = CD_ZERO
    ! Optional Morison drag (transverse + tangential quadratic drag on the submerged cable, relative
    ! to a uniform fluid velocity). Off by default; wired in via CD_HermiteCable_Dyn_Set_Drag.
    LOGICAL :: has_drag = .FALSE.
    REAL(wp) :: rho_w = CD_ZERO, waterline_z = CD_ZERO, current(3) = CD_ZERO
    REAL(wp), ALLOCATABLE :: diam(:), cdn(:), cdt(:)             ! per element
    ! The drag configuration also defines the free surface and the displaced section
    ! (rho_w, diam): the part of the line above the surface regains its buoyancy
    ! rho_w*gravity*pi*diam^2/4 per reference length (see CD_HermiteCable_Dry_Buoyancy).
    REAL(wp) :: gravity = 9.80665_wp
    ! Optional Morison added mass (normal/tangential fluid-inertia coefficients on the submerged
    ! cable). Off by default; wired in via CD_HermiteCable_Dyn_Set_AddedMass. Kept independent of the
    ! drag configuration so either effect can be enabled alone.
    LOGICAL :: has_am = .FALSE.
    REAL(wp) :: am_rho = CD_ZERO, am_wl = CD_ZERO
    REAL(wp), ALLOCATABLE :: am_diam(:), am_can(:), am_cat(:)    ! per element
    ! Optional regular (Airy) wave field. Off by default; wired in via CD_HermiteCable_Dyn_Set_Waves.
    ! The wave velocity feeds the drag's fluid field (requires the drag config); the wave
    ! acceleration feeds the Froude-Krylov + fluid-inertia load with coefficients (1 + Can/Cat) from
    ! the added-mass config (requires it).
    LOGICAL :: has_wave = .FALSE.
    REAL(wp) :: wv_h = CD_ZERO, wv_om = CD_ZERO, wv_k = CD_ZERO, wv_depth = CD_ZERO, wv_dir = CD_ZERO
    ! Optional COMPONENT (irregular) wave field, wired in via
    ! CD_HermiteCable_Dyn_Set_Irregular_Waves: a caller-supplied long-crested component
    ! table (amplitude/frequency/wavenumber/phase per component) sharing wv_depth/wv_dir
    ! with the scalar field above. wv_ncomp = 0 -> the scalar regular field is active
    ! (when has_wave); wv_ncomp > 0 -> the component field is active and the scalar
    ! trio is cleared. The two are mutually exclusive by construction.
    INTEGER :: wv_ncomp = 0
    REAL(wp), ALLOCATABLE :: wv_ampc(:), wv_omc(:), wv_kc(:), wv_phc(:)
    ! exp(-2 wv_kc(i) wv_depth), fixed with the table (refresh_wave_table)
    REAL(wp), ALLOCATABLE :: wv_e2kd(:)
    ! Optional HELD ambient-fluid field, wired in via CD_HermiteCable_Dyn_Set_Held_Fluid:
    ! per-NODE fluid velocity/acceleration and local free-surface elevation prescribed by
    ! an external host (e.g. the OpenFAST shell sampling a SeaState WaveField) and held
    ! FROZEN over each step -- the weak-coupling cadence the coupled loads-out side
    ! documents. Consumed per ELEMENT as the two end nodes' average (the same
    ! granularity the host sampled at). Mutually exclusive with the internal wave
    ! families; requires the drag config (velocity) and adds Froude-Krylov forcing when
    ! the added-mass config is present (acceleration).
    LOGICAL :: has_held_field = .FALSE.
    ! Line-end bending connection at node 1 and node nn. Finite mode uses the
    ! isotropic spring; rigid mode is an exact tangent-direction constraint. The
    ! default pinned mode contributes nothing and remains bit-identical.
    LOGICAL :: has_endconn = .FALSE.
    INTEGER :: endconn_mode(2) = CD_ENDCONN_PINNED
    REAL(wp) :: endconn_k(2) = CD_ZERO
    REAL(wp) :: endconn_d0(3, 2) = CD_ZERO
    REAL(wp), ALLOCATABLE :: hf_u(:, :), hf_ud(:, :), hf_wl(:)   ! (3, nn), (3, nn), (nn)
    ! Discrete point attachments (clumps, buoyancy modules) lumped at nodes: node att_node(k)
    ! carries the net weight att_w(k) = (m - rho V) g [N], the buoyancy att_b(k) = rho g V
    ! [N] it loses above the free surface, the fluid-inertia coefficient att_fk(k) =
    ! rho V (1 + Ca) [kg], and the drag areas normal (att_cda) and axial (att_cdax) to the
    ! line tangent at the node [m^2] with the fluid density att_rho. Its mass
    ! plus added mass m + Ca rho V is part of the constant mass Mgb.
    LOGICAL :: has_attach = .FALSE.
    ! Set by a checkpoint-restart reload (CD_HFMF_UnpackMirror): the committed acceleration is
    ! a generalised-alpha history value, so the first held fluid field after the reload must
    ! preserve it instead of re-deriving an initial-condition acceleration.
    LOGICAL :: a_restored = .FALSE.
    INTEGER, ALLOCATABLE :: att_node(:)
    REAL(wp), ALLOCATABLE :: att_w(:), att_b(:), att_fk(:), att_cda(:), att_cdax(:)
    REAL(wp) :: att_rho = CD_ZERO
    ! Rollback copies of CD_HermiteCable_Dyn_Set_Held_Fluid, kept so a per-step field update
    ! allocates nothing: a (ndof), hf_u/hf_ud (3, nn), hf_wl (nn).
    REAL(wp), ALLOCATABLE :: ws_hf_a_old(:), ws_hf_u_old(:, :), ws_hf_ud_old(:, :), ws_hf_wl_old(:)
    ! Simulation time (s): 0 at Init, advanced by Step; wave loads evaluate at the generalised-alpha
    ! intermediate time within each step.
    REAL(wp) :: t = CD_ZERO
    ! Snapshot of the COMMITTED step state (q, v, a, t) for the aggregate stage-then-commit
    ! contract: a multi-object caller snapshots every sub-object before stepping and
    ! restores them all if any later sub-step fails, so a failed aggregate step can never
    ! leave partially advanced state. Buffers allocate on first Snapshot (amortised zero;
    ! never touched by Step itself).
    REAL(wp), ALLOCATABLE :: snap_q(:), snap_v(:), snap_a(:)
    REAL(wp) :: snap_t = CD_ZERO
    REAL(wp) :: snap_endconn_d0(3, 2) = CD_ZERO
    ! Condensed torsion (CD_HermiteCable_Dyn_Set_Torsion): end frames, GJ, imposed twist and the
    ! accepted Theta of the static solve. The residual, the end loads and the energy include it;
    ! the time step does not yet (a model with torsion refuses to step). snap_torsion_theta is
    ! the Theta of the last snapshot.
    TYPE(CD_HermiteTorsionType) :: torsion
    REAL(wp) :: snap_torsion_theta = CD_ZERO
    LOGICAL :: snap_valid = .FALSE.
    ! Adaptive-Newton warmup counters are step-history too: a staged step (aggregate
    ! stage-then-commit) that advances them and is then rolled back must restore them, or the
    ! warmup measurement/decision desyncs from the committed step count.
    INTEGER :: snap_na_steps_done = 0
    REAL(wp) :: snap_na_iter_sum = CD_ZERO
    LOGICAL :: snap_na_mn_active = .FALSE.
    INTEGER :: snap_tensile_event_count = 0, snap_tensile_worst_element = 0
    REAL(wp) :: snap_tensile_worst_force = CD_ZERO, snap_tensile_worst_threshold = CD_ZERO
    REAL(wp) :: snap_tensile_worst_xi = CD_ZERO, snap_tensile_worst_time = CD_ZERO
    REAL(wp) :: snap_max_step_rotation = CD_ZERO
    INTEGER :: snap_recovery_event_count = 0, snap_recovery_max_substeps_used = 0
    ! Modified Newton (tangent reuse): factor the effective tangent ONCE per step (at the
    ! predictor) and reuse the factorization for every subsequent iteration of that step;
    ! a line search that finds no improving step under the stale direction triggers ONE
    ! refresh (rebuild + refactor at the current iterate) before the safe damped fallback.
    ! Within-step reuse is off by default (cross-step reuse, below, is on). Wired in
    ! via CD_HermiteCable_Dyn_Set_ModifiedNewton (the deck OPTIONS modified_newton).
    LOGICAL :: modified_newton = .FALSE.
    ! Cross-step factor reuse (on by default; CD_HermiteCable_Dyn_Set_TangentReuse): the
    ! factored effective tangent of the last committed step (ws_Keffb/ws_ipiv) starts the next
    ! step's Newton when that step continues bitwise from the committed state (the force-cache
    ! key) at the same dt and its predictor test finds the motion smooth (the constant-
    ! acceleration guess runs). The factor is kept while each accepted iterate contracts the residual
    ! to <= 1/4 of the previous one (the modified-Newton gate),
    ! for at most three iterations, and is otherwise refreshed at the current iterate (full
    ! Newton for the rest of the step); a stale direction that finds no improving step refreshes
    ! once before the damped fallback, and a step that fails on the carried factor is re-solved
    ! from the Newmark guess on a fresh tangent. Only the Newton matrix is reused: every step
    ! converges on the exact residual (to 0.3 tol while on the carried factor, so the outputs
    ! stay as close to the converged solution as full Newton's), but the iterate it stops at depends
    ! on the factor's age, which no checkpoint carries, so a restart is exact only with reuse
    ! off. tr_valid marks a usable factor for step size tr_dt; any other writer of ws_Keffb,
    ! every setter and Restore clear it.
    LOGICAL :: tangent_reuse = .TRUE.
    ! Back-off: a step on the carried factor that still factorised at least as often as the
    ! last fresh step (xs_fresh_fact factorisations) saved nothing -- in a fast-changing
    ! geometry, Newton from the stale factor's iterate costs more than from the predictor. It
    ! makes the next xs_skip steps factor fresh; xs_backoff doubles on each consecutive such
    ! step (to XS_MAX_BACKOFF) and resets when a reuse step saves a factorisation.
    INTEGER :: xs_skip = 0, xs_backoff = 0, xs_fresh_fact = 0
    LOGICAL :: tr_valid = .FALSE.
    REAL(wp) :: tr_dt = CD_ZERO
    ! Adaptive modified Newton: run the first na_warmup steps in FULL Newton to measure the
    ! per-step iteration count, then enable tangent reuse for the rest of the run IFF the warmup
    ! averaged more than na_threshold iterations/step -- reuse pays on harder/larger systems
    ! and costs the refresh overhead on easy ones. Off by default (the
    ! explicit modified_newton flag still forces reuse); set via CD_HermiteCable_Dyn_Set_AdaptiveNewton.
    LOGICAL :: newton_adaptive = .FALSE.
    INTEGER :: na_warmup = 4
    REAL(wp) :: na_threshold = 5.0_wp
    INTEGER :: na_steps_done = 0
    REAL(wp) :: na_iter_sum = 0.0_wp
    LOGICAL :: na_mn_active = .FALSE.
    ! Continuous axial-resultant audit. The compatibility logical is true only for the
    ! hard-reject mode; tensile_mode also provides a non-fatal, counted warning mode for long runs.
    LOGICAL :: require_tensile = .FALSE.
    INTEGER :: tensile_mode = CD_HCDYN_TENSILE_OFF
    REAL(wp) :: tensile_strain_tolerance = CD_HCDYN_AXIAL_STRAIN_TOL
    INTEGER :: tensile_event_count = 0, tensile_worst_element = 0
    REAL(wp) :: tensile_worst_force = CD_ZERO, tensile_worst_threshold = CD_ZERO
    REAL(wp) :: tensile_worst_xi = CD_ZERO, tensile_worst_time = CD_ZERO
    INTEGER :: recovery_event_count = 0, recovery_max_substeps_used = 0
    LOGICAL :: initialized = .FALSE.
  END TYPE CD_HermiteCableDynType

CONTAINS

  LOGICAL FUNCTION CD_HermiteCable_Dyn_Profile_Enabled() RESULT(enabled)
    !! Process-global profiling counters are intentionally unsynchronised diagnostics.
    !! Aggregate callers query this flag to keep profiled multi-cable runs serial.
    enabled = prof_enabled
  END FUNCTION CD_HermiteCable_Dyn_Profile_Enabled

  SUBROUTINE CD_HermiteCable_Dyn_Reset_Profile()
    !! Zero the diagnostic step-profile counters and ARM the instrumentation (call before a
    !! timed run). Profiling stays off until this is called, so production runs never pay
    !! the SYSTEM_CLOCK overhead in the step/residual hot path.
    prof_enabled = .TRUE.
    prof_n_step = 0; prof_n_resid = 0; prof_n_solve = 0; prof_n_am = 0; prof_n_alloc = 0
    prof_n_tan = 0
    prof_t_step = 0.0_wp; prof_t_resid = 0.0_wp; prof_t_eff = 0.0_wp
    prof_t_tan = 0.0_wp
    prof_t_solve = 0.0_wp; prof_t_am = 0.0_wp
  END SUBROUTINE CD_HermiteCable_Dyn_Reset_Profile

  SUBROUTINE CD_HermiteCable_Dyn_Disable_Profile()
    !! Disarm process-global diagnostic timing. A host library can execute several
    !! sequential cases in one process; a non-profiled case must not inherit the
    !! timing overhead or serial multi-cable policy from an earlier profiled case.
    prof_enabled = .FALSE.
  END SUBROUTINE CD_HermiteCable_Dyn_Disable_Profile

  SUBROUTINE CD_HermiteCable_Dyn_Get_Profile(n_step, n_resid, n_solve, n_am, n_alloc, &
                                             t_step, t_resid, t_eff, t_solve, t_am, n_tan, t_tan)
    !! Read the diagnostic step-profile counters accumulated since the last reset.
    !! t_step is the total wall inside CD_HermiteCable_Dyn_Step; t_resid + t_eff +
    !! t_solve + t_am partition its instrumented interior (the remainder is
    !! predictor/commit bookkeeping and workspace allocation).
    INTEGER, INTENT(OUT) :: n_step, n_resid, n_solve, n_am, n_alloc
    REAL(wp), INTENT(OUT) :: t_step, t_resid, t_eff, t_solve, t_am
    ! n_tan/t_tan report the tangent assemblies separately from the residual evaluations,
    ! so the eval_res averages stay honest (OPTIONAL so callers may omit them).
    INTEGER, INTENT(OUT), OPTIONAL :: n_tan
    REAL(wp), INTENT(OUT), OPTIONAL :: t_tan
    n_step = prof_n_step; n_resid = prof_n_resid; n_solve = prof_n_solve
    n_am = prof_n_am; n_alloc = prof_n_alloc
    t_step = prof_t_step; t_resid = prof_t_resid; t_eff = prof_t_eff
    t_solve = prof_t_solve; t_am = prof_t_am
    IF (PRESENT(n_tan)) n_tan = prof_n_tan
    IF (PRESENT(t_tan)) t_tan = prof_t_tan
  END SUBROUTINE CD_HermiteCable_Dyn_Get_Profile

  SUBROUTINE CD_HermiteCable_Dyn_Init(model, l0, EA, EI, rho_a, w, seed, fixed_dofs, &
                                      seabed_z, kn, rho_inf, ErrStat, ErrMsg, &
                                      axial_quadrature_order, bending_quadrature_order)
    !! Initialise the Hermite dynamic model: validate, assemble the constant mass, store the
    !! at-rest state (v = 0), and solve the consistent initial acceleration M a0 = f_ext - f_int.
    !!
    !! l0/EA/EI/rho_a/w(ne) : per-element rest length, axial/bending stiffness, mass per length,
    !!                        signed submerged weight per length (>0 net-heavy pulls -z).
    !! seed(6*nn)           : initial [r,m] per node (nn = ne + 1); the static equilibrium seed.
    !! fixed_dofs(:)        : global DOF indices held fixed (v = a = 0 there).
    !! seabed_z, kn         : penalty seabed plane and stiffness (kn >= 0).
    !! rho_inf              : gen-alpha spectral radius in [0,1] (0.8 = light high-freq damping).
    TYPE(CD_HermiteCableDynType), INTENT(OUT) :: model
    REAL(wp), INTENT(IN) :: l0(:), EA(:), EI(:), rho_a(:), w(:), seed(:)
    INTEGER, INTENT(IN) :: fixed_dofs(:)
    REAL(wp), INTENT(IN) :: seabed_z, kn, rho_inf
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: axial_quadrature_order, bending_quadrature_order

    INTEGER :: ne, nn, ndof, i, e, es, nf, gi
    REAL(wp) :: Me(12, 12), rho
    CHARACTER(200) :: em2

    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    model%initialized = .FALSE.
    model%axial_quadrature_order = 4
    model%bending_quadrature_order = 4
    IF (PRESENT(axial_quadrature_order)) model%axial_quadrature_order = axial_quadrature_order
    IF (PRESENT(bending_quadrature_order)) model%bending_quadrature_order = bending_quadrature_order

    ne = SIZE(l0)
    nn = ne + 1
    ndof = 6*nn
    IF (ne < 1 .OR. SIZE(EA) /= ne .OR. SIZE(EI) /= ne .OR. SIZE(rho_a) /= ne .OR. SIZE(w) /= ne) THEN
      CALL fail('l0/EA/EI/rho_a/w must be same length ne >= 1'); RETURN
    END IF
    IF (SIZE(seed) /= ndof) THEN
      CALL fail('seed must be 6*(ne+1)'); RETURN
    END IF
    IF (.NOT. CD_Is_Finite(rho_inf) .OR. rho_inf < CD_ZERO .OR. rho_inf > CD_ONE) THEN
      CALL fail('rho_inf must be finite in [0,1]'); RETURN
    END IF
    IF (model%axial_quadrature_order < 1 .OR. model%axial_quadrature_order > 6 .OR. &
        model%bending_quadrature_order < 1 .OR. model%bending_quadrature_order > 6) THEN
      CALL fail('quadrature orders must lie in [1,6]'); RETURN
    END IF
    IF (kn < CD_ZERO .OR. .NOT. CD_Is_Finite(kn) .OR. .NOT. CD_Is_Finite(seabed_z)) THEN
      CALL fail('kn>=0 and finite seabed_z required'); RETURN
    END IF
    IF (.NOT. CD_All_Finite(seed) .OR. .NOT. CD_All_Finite(l0) .OR. &
        .NOT. CD_All_Finite(EA) .OR. .NOT. CD_All_Finite(EI) .OR. &
        .NOT. CD_All_Finite(rho_a) .OR. .NOT. CD_All_Finite(w)) THEN
      CALL fail('inputs must be finite'); RETURN
    END IF
    IF (ANY(l0 <= CD_ZERO) .OR. ANY(EA < CD_ZERO) .OR. ANY(EI < CD_ZERO) .OR. ANY(rho_a < CD_ZERO)) THEN
      CALL fail('l0>0, EA>=0, EI>=0, rho_a>=0 required'); RETURN
    END IF
    ! A dynamic model needs mass: an all-massless line gives a singular consistent mass, so the
    ! initial-acceleration solve M a0 = f_ext - f_int has no solution. Reject it with a clear message
    ! rather than letting the linear solve fail cryptically (partial masslessness that still leaves
    ! every free node with mass is allowed and is caught by the solve if it does turn singular).
    IF (ALL(rho_a <= CD_ZERO)) THEN
      CALL fail('a dynamic model needs positive mass (all rho_a = 0)'); RETURN
    END IF
    IF (ANY(fixed_dofs < 1) .OR. ANY(fixed_dofs > ndof)) THEN
      CALL fail('fixed_dofs out of range'); RETURN
    END IF

    model%ne = ne; model%nn = nn; model%ndof = ndof
    model%seabed_z = seabed_z; model%kn = kn
    ALLOCATE (model%l0(ne), model%EA(ne), model%EI(ne), model%rho_a(ne), model%w(ne))
    model%l0 = l0; model%EA = EA; model%EI = EI; model%rho_a = rho_a; model%w = w

    ! gen-alpha constants
    rho = rho_inf
    model%alpha_m = (2.0_wp*rho - CD_ONE)/(rho + CD_ONE)
    model%alpha_f = rho/(rho + CD_ONE)
    model%beta = 0.25_wp*(CD_ONE - model%alpha_m + model%alpha_f)**2
    model%gamma = 0.5_wp - model%alpha_m + model%alpha_f

    ! nodal tributary lengths (weight + seabed penalty)
    ALLOCATE (model%trib(nn))
    model%trib = CD_ZERO
    DO e = 1, ne
      model%trib(e) = model%trib(e) + 0.5_wp*l0(e)
      model%trib(e + 1) = model%trib(e + 1) + 0.5_wp*l0(e)
    END DO

    ! free/fixed partition
    ALLOCATE (model%freemask(ndof), model%solve_mask(ndof))
    model%freemask = .TRUE.
    DO i = 1, SIZE(fixed_dofs)
      model%freemask(fixed_dofs(i)) = .FALSE.
    END DO
    model%solve_mask = model%freemask
    nf = COUNT(model%freemask)
    ALLOCATE (model%fmap(nf))
    gi = 0
    DO i = 1, ndof
      IF (model%freemask(i)) THEN
        gi = gi + 1; model%fmap(gi) = i
      END IF
    END DO

    ! constant consistent mass, assembled once into band storage
    ALLOCATE (model%Mgb(LDAB_D, ndof))
    model%Mgb = CD_ZERO
    DO e = 1, ne
      CALL CD_HermiteCable_Mass(rho_a(e), l0(e), Me, es, em2)
      IF (es /= CD_HCABLE_OK) THEN
        CALL fail('mass assembly: '//TRIM(em2)); RETURN
      END IF
      CALL scatter12_band(model%Mgb, Me, e)
    END DO

    ! step workspace, allocated once (a dynamic step performs zero heap allocations)
    ALLOCATE (model%ws_qn(ndof), model%ws_vn(ndof), model%ws_an(ndof), model%ws_qk(ndof), &
              model%ws_ak(ndof), model%ws_qa(ndof), model%ws_aa(ndof), model%ws_va(ndof), &
              model%ws_vk(ndof), model%ws_R(ndof), model%ws_G(ndof), model%ws_Gt(ndof), &
              model%ws_inert(ndof), model%ws_qtrial(ndof), model%ws_dq(ndof), &
              model%ws_recovery_q0(ndof), model%ws_recovery_v0(ndof), model%ws_recovery_a0(ndof))
    ALLOCATE (model%ws_Rk(ndof), model%ws_Rt(ndof), model%fc_R(ndof), model%fc_q(ndof), model%fc_v(ndof))
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    model%force_blend = CD_HC_FORCE_BLEND_DEFAULT
    model%max_step_rotation = CD_ZERO
    ALLOCATE (model%ws_Kb(LDAB_D, ndof), model%ws_Kvb(LDAB_D, ndof), &
              model%ws_Keffb(LDAB_D, ndof), model%ws_Msumb(LDAB_D, ndof), &
              model%ws_Mamb(LDAB_D, ndof))
    ALLOCATE (model%ws_ipiv(ndof))
    ALLOCATE (model%ws_elem_f(12, ne), model%ws_elem_K(12, 12, ne), model%ws_elem_Kdrag(12, 12, ne), &
              model%ws_elem_Kv(12, 12, ne), model%ws_elem_Kwave(12, 12, ne), model%ws_elem_Kheld(12, 12, ne), &
              model%ws_elem_es(ne), model%ws_elem_em(ne))

    ! state: at-rest seed, v = 0, consistent a0
    ALLOCATE (model%q(ndof), model%v(ndof), model%a(ndof))
    model%q = seed
    model%v = CD_ZERO
    model%a = CD_ZERO

    ! Consistent initial acceleration M a0 = f_ext - f_int at the seed state (v = 0). Shared with
    ! Set_Drag so a0 is refreshed when drag is added (e.g. a current gives a non-zero rest drag).
    CALL solve_consistent_acceleration(model, es, em2)
    IF (es /= CD_HCDYN_OK) THEN
      CALL fail(TRIM(em2)); RETURN
    END IF

    model%initialized = .TRUE.

  CONTAINS
    SUBROUTINE fail(msg)
      CHARACTER(*), INTENT(IN) :: msg
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Init: '//msg
    END SUBROUTINE fail
  END SUBROUTINE CD_HermiteCable_Dyn_Init

  SUBROUTINE CD_HermiteCable_Dyn_Set_Contact(model, kn_node, cn_node, mu, ErrStat, ErrMsg, &
                                             bathymetry, frame_cs)
    !! Install the standalone deck's nodal seabed law on an initialized Hermite model.
    !! The force/Jacobian convention is that of the deck's nodal seabed law:
    !! C1 normal penalty, compression-only dashpot, and (mu > 0) stick-slip friction springs
    !! anchored at the present node positions (CD_Seabed_Friction_Stick_Slip).
    !! A structured bathymetry is queried in global coordinates; frame_cs maps the
    !! chord-local Hermite solve frame back to global x/y before each query.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: kn_node(:), cn_node(:), mu
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    TYPE(CD_BathymetryType), INTENT(IN), OPTIONAL :: bathymetry
    REAL(wp), INTENT(IN), OPTIONAL :: frame_cs(2)

    INTEGER :: es, istat, i_fr
    LOGICAL :: old_fr_active, old_fr_aniso
    REAL(wp), ALLOCATABLE :: old_fr_anchor(:, :), old_fr_k(:), old_fr_recovery(:, :), old_fr_snap(:, :)
    REAL(wp) :: cs(2), old_mu, old_fc, old_fs, old_mu_axial
    REAL(wp), ALLOCATABLE :: new_kn(:), new_cn(:), old_kn(:), old_cn(:), old_a(:)
    LOGICAL :: old_has_contact, old_has_bathymetry
    TYPE(CD_BathymetryType) :: new_bathymetry, old_bathymetry
    CHARACTER(200) :: em

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    ! A load-configuration change invalidates the cached committed force F_n.
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    IF (.NOT. model%initialized) THEN
      CALL fail_contact('model not initialized'); RETURN
    END IF
    IF (SIZE(kn_node) /= model%nn .OR. SIZE(cn_node) /= model%nn) THEN
      CALL fail_contact('kn_node/cn_node must have length n_nodes'); RETURN
    END IF
    IF (.NOT. CD_All_Finite(kn_node) .OR. ANY(kn_node <= CD_ZERO) .OR. &
        .NOT. CD_All_Finite(cn_node) .OR. ANY(cn_node < CD_ZERO) .OR. &
        .NOT. CD_Is_Finite(mu) .OR. mu < CD_ZERO) THEN
      CALL fail_contact('need finite kn_node > 0, cn_node >= 0, and mu >= 0'); RETURN
    END IF
    cs = [CD_ONE, CD_ZERO]
    IF (PRESENT(frame_cs)) THEN
      IF (.NOT. CD_All_Finite(frame_cs) .OR. &
          ABS(SUM(frame_cs*frame_cs) - CD_ONE) > 1.0e-9_wp) THEN
        CALL fail_contact('frame_cs must be a finite unit heading'); RETURN
      END IF
      cs = frame_cs
    END IF
    IF (PRESENT(bathymetry)) THEN
      IF (.NOT. CD_Bathymetry_Is_Initialized(bathymetry)) THEN
        CALL fail_contact('bathymetry is not initialized'); RETURN
      END IF
      new_bathymetry = bathymetry
    END IF
    ALLOCATE (new_kn(model%nn), new_cn(model%nn), old_a(model%ndof), STAT=istat)
    IF (istat /= 0) THEN
      CALL fail_contact('workspace allocation failed'); RETURN
    END IF
    new_kn = kn_node
    new_cn = cn_node

    ! Preserve the complete committed contact configuration until the acceleration
    ! refresh succeeds. Setters are allowed after initialization, so a failed
    ! replacement must not disable or partially overwrite a previously valid bed.
    old_has_contact = model%has_contact
    old_has_bathymetry = model%has_contact_bathymetry
    old_mu = model%contact_mu; old_fc = model%contact_frame_c; old_fs = model%contact_frame_s
    old_bathymetry = model%contact_bathymetry
    old_a = model%a
    CALL MOVE_ALLOC(model%contact_kn, old_kn)
    CALL MOVE_ALLOC(model%contact_cn, old_cn)
    CALL MOVE_ALLOC(new_kn, model%contact_kn)
    CALL MOVE_ALLOC(new_cn, model%contact_cn)
    model%contact_bathymetry = new_bathymetry
    model%has_contact_bathymetry = PRESENT(bathymetry)
    model%contact_mu = mu
    model%contact_frame_c = cs(1)
    model%contact_frame_s = cs(2)
    model%has_contact = .TRUE.
    ! Seabed friction (mu > 0): stick-slip springs anchored at the present node positions,
    ! i.e. carrying no force (CD_HermiteCable_Dyn_Set_Friction_Anchors installs the force
    ! of a static solve with friction). Workspace is allocated once per model.
    old_fr_active = model%fr_active
    old_fr_aniso = model%contact_fr_aniso
    old_mu_axial = model%contact_mu_axial
    IF (ALLOCATED(model%fr_anchor)) THEN
      old_fr_anchor = model%fr_anchor
      old_fr_k = model%fr_k
      old_fr_recovery = model%ws_recovery_fr
      old_fr_snap = model%snap_fr_anchor
    END IF
    ! A new coefficient resets the friction to isotropic; CD_HermiteCable_Dyn_Set_Friction_Axial
    ! (called after this setter) makes it anisotropic against the new lateral coefficient.
    model%contact_fr_aniso = .FALSE.
    model%contact_mu_axial = mu
    model%fr_active = mu > CD_ZERO
    IF (model%fr_active) THEN
      IF (.NOT. ALLOCATED(model%fr_anchor)) THEN
        ALLOCATE (model%fr_anchor(2, model%nn), model%fr_k(model%nn), model%ws_recovery_fr(2, model%nn), &
                  model%snap_fr_anchor(2, model%nn))
        model%snap_fr_anchor = CD_ZERO
      END IF
      model%fr_k = model%contact_kn
      DO i_fr = 1, model%nn
        model%fr_anchor(:, i_fr) = model%q(6*i_fr - 5:6*i_fr - 4)
      END DO
      model%ws_recovery_fr = model%fr_anchor
      ! The new anchors are the committed ones: a later Restore must not bring back
      ! the anchors of the replaced configuration.
      model%snap_fr_anchor = model%fr_anchor
    END IF
    ! The contact set is part of f_ext at t=0; refresh a0 so the first predictor
    ! starts from the same residual subsequently evaluated inside every Newton step.
    CALL solve_consistent_acceleration(model, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      IF (ALLOCATED(model%contact_kn)) DEALLOCATE (model%contact_kn)
      IF (ALLOCATED(model%contact_cn)) DEALLOCATE (model%contact_cn)
      CALL MOVE_ALLOC(old_kn, model%contact_kn)
      CALL MOVE_ALLOC(old_cn, model%contact_cn)
      model%contact_bathymetry = old_bathymetry
      model%has_contact = old_has_contact
      model%has_contact_bathymetry = old_has_bathymetry
      model%contact_mu = old_mu; model%contact_frame_c = old_fc; model%contact_frame_s = old_fs
      model%fr_active = old_fr_active
      model%contact_fr_aniso = old_fr_aniso
      model%contact_mu_axial = old_mu_axial
      IF (ALLOCATED(old_fr_anchor)) THEN
        model%fr_anchor = old_fr_anchor
        model%fr_k = old_fr_k
        model%ws_recovery_fr = old_fr_recovery
        model%snap_fr_anchor = old_fr_snap
      END IF
      model%a = old_a
      CALL fail_contact('initial acceleration refresh failed: '//TRIM(em)); RETURN
    END IF

  CONTAINS
    SUBROUTINE fail_contact(msg)
      CHARACTER(*), INTENT(IN) :: msg
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Contact: '//msg
    END SUBROUTINE fail_contact
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Contact

  SUBROUTINE CD_HermiteCable_Dyn_Set_Friction_Axial(model, mu_axial, ErrStat, ErrMsg)
    !! Make the stick-slip seabed friction anisotropic: contact_mu becomes the lateral
    !! coefficient and mu_axial the coefficient along the horizontal projection of the
    !! nodal tangent (CD_Seabed_Friction_Aniso). Call after CD_HermiteCable_Dyn_Set_Contact;
    !! mu_axial equal to contact_mu keeps the isotropic law bit-for-bit.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: mu_axial
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. (CD_Is_Finite(mu_axial) .AND. mu_axial > CD_ZERO)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Friction_Axial: the axial coefficient must be finite and positive'
      RETURN
    END IF
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Friction_Axial: model not initialised'
      RETURN
    END IF
    IF (.NOT. model%fr_active) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Friction_Axial: seabed contact with friction is required'
      RETURN
    END IF
    model%contact_mu_axial = mu_axial
    model%contact_fr_aniso = ABS(mu_axial - model%contact_mu) > CD_ZERO
    ! A friction-law change invalidates the cached committed force F_n and the carried factor.
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Friction_Axial

  SUBROUTINE CD_HermiteCable_Dyn_Set_Friction_Anchors(model, anchors, ErrStat, ErrMsg)
    !! Switch the seabed friction of a contact model with mu > 0 to elastic-plastic friction
    !! springs anchored at anchors(1:2, node) (x, y in the solve frame). The spring stiffness
    !! is the nodal normal contact stiffness. A static solve with seabed friction gives each
    !! node its spring force f_s; anchors = x - f_s/k then reproduce f_s at the committed
    !! state, so the march starts in the static equilibrium with every node sticking.
    !! The consistent initial acceleration is refreshed; failure leaves the model unchanged.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: anchors(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), ALLOCATABLE :: a_old(:), anchor_old(:, :)
    LOGICAL :: active_old
    INTEGER :: es
    CHARACTER(200) :: em

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      CALL fail_fr('model not initialised'); RETURN
    END IF
    IF (.NOT. model%has_contact .OR. .NOT. (model%contact_mu > CD_ZERO)) THEN
      CALL fail_fr('seabed contact with a positive friction coefficient is required'); RETURN
    END IF
    IF (SIZE(anchors, 1) /= 2 .OR. SIZE(anchors, 2) /= model%nn) THEN
      CALL fail_fr('anchors must be (2, n_nodes)'); RETURN
    END IF
    IF (.NOT. CD_All_Finite(anchors)) THEN
      CALL fail_fr('anchors must be finite'); RETURN
    END IF
    IF (.NOT. ALLOCATED(model%fr_anchor)) THEN
      ALLOCATE (model%fr_anchor(2, model%nn), model%fr_k(model%nn), model%ws_recovery_fr(2, model%nn), &
                model%snap_fr_anchor(2, model%nn))
      model%fr_anchor = CD_ZERO
      model%snap_fr_anchor = CD_ZERO
    END IF
    a_old = model%a
    anchor_old = model%fr_anchor
    active_old = model%fr_active
    model%fr_k = model%contact_kn
    model%fr_anchor = anchors
    model%fr_active = .TRUE.
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    CALL solve_consistent_acceleration(model, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      model%a = a_old
      model%fr_anchor = anchor_old
      model%fr_active = active_old
      CALL fail_fr('initial acceleration refresh failed: '//TRIM(em)); RETURN
    END IF
    ! The installed anchors are the committed ones: a later Restore must not bring back
    ! the anchors they replace.
    model%snap_fr_anchor = model%fr_anchor
    model%ws_recovery_fr = model%fr_anchor

  CONTAINS
    SUBROUTINE fail_fr(msg)
      CHARACTER(*), INTENT(IN) :: msg
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Friction_Anchors: '//msg
    END SUBROUTINE fail_fr
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Friction_Anchors

  SUBROUTINE CD_HermiteCable_Dyn_Set_Attachments(model, node, mass, volume, cda, cdax, ca, rho_w, gravity, &
                                                 ErrStat, ErrMsg)
    !! Install discrete point attachments (clumps, buoyancy modules) lumped at nodes. Each
    !! carries its weight and buoyancy (m - rho V) g downward, isotropic added mass
    !! Ca rho V (with the dry mass m, a constant nodal mass on the translational DOFs), the
    !! Froude-Krylov plus fluid-inertia force rho V (1 + Ca) du/dt and the quadratic drag
    !! rho [CdA |w_n| w_n + CdAx |w_t| w_t] / 2 of the relative velocity w = u - v split
    !! normal and axial to the line tangent at the node (CD_HermiteCable_Attachment_Drag),
    !! against the held fluid field at the node (still water without one). Above the free
    !! surface an attachment loses its buoyancy and fluid loads. Called once, after Init and
    !! never together with internal waves; the initial acceleration is refreshed. The
    !! kinetic energy (CD_HermiteCable_Dyn_Energy) then includes the attachment mass and its
    !! added mass Ca rho V.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER, INTENT(IN) :: node(:)
    REAL(wp), INTENT(IN) :: mass(:), volume(:), cda(:), cdax(:), ca(:), rho_w, gravity
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n, k, d, dof, es
    REAL(wp), ALLOCATABLE :: a_old(:), mgb_old(:, :)
    CHARACTER(200) :: em

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    n = SIZE(node)
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Set_Attachments: model not initialised'; RETURN
    END IF
    IF (model%has_attach) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Set_Attachments: attachments are already set'
      RETURN
    END IF
    IF (model%has_wave) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Attachments: attachments take the held fluid field, not internal waves'
      RETURN
    END IF
    IF (n < 1 .OR. SIZE(mass) /= n .OR. SIZE(volume) /= n .OR. SIZE(cda) /= n .OR. SIZE(ca) /= n .OR. &
        SIZE(cdax) /= n) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Attachments: node/mass/volume/cda/ca must share one length >= 1'
      RETURN
    END IF
    IF (ANY(node < 1) .OR. ANY(node > model%nn)) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Set_Attachments: node out of range'; RETURN
    END IF
    IF (.NOT. (CD_All_Finite(mass) .AND. CD_All_Finite(volume) .AND. CD_All_Finite(cda) .AND. &
               CD_All_Finite(cdax) .AND. CD_All_Finite(ca) .AND. CD_Is_Finite(rho_w) .AND. CD_Is_Finite(gravity))) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Set_Attachments: values must be finite'; RETURN
    END IF
    IF (ANY(mass < CD_ZERO) .OR. ANY(volume < CD_ZERO) .OR. ANY(cda < CD_ZERO) .OR. ANY(ca < CD_ZERO) .OR. &
        ANY(cdax < CD_ZERO) .OR. &
        rho_w <= CD_ZERO .OR. gravity <= CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Attachments: mass, volume, CdA and Ca must be >= 0; rho_w, g > 0'
      RETURN
    END IF
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    a_old = model%a
    mgb_old = model%Mgb
    ALLOCATE (model%att_node(n), model%att_w(n), model%att_b(n), model%att_fk(n), model%att_cda(n), &
              model%att_cdax(n))
    model%att_node = node
    model%att_w = (mass - rho_w*volume)*gravity
    model%att_b = rho_w*gravity*volume
    model%att_fk = rho_w*volume*(CD_ONE + ca)
    model%att_cda = cda
    model%att_cdax = cdax
    model%att_rho = rho_w
    DO k = 1, n
      DO d = 1, 3
        dof = 6*(node(k) - 1) + d
        model%Mgb(KL_D + KU_D + 1, dof) = model%Mgb(KL_D + KU_D + 1, dof) + mass(k) + ca(k)*rho_w*volume(k)
      END DO
    END DO
    model%has_attach = .TRUE.
    CALL solve_consistent_acceleration(model, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      DEALLOCATE (model%att_node, model%att_w, model%att_b, model%att_fk, model%att_cda, model%att_cdax)
      model%has_attach = .FALSE.
      model%Mgb = mgb_old
      model%a = a_old
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Attachments: initial acceleration refresh failed: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Attachments

  SUBROUTINE CD_HermiteCable_Attachment_Drag(rel_velocity, tangent, rho, cda, cdax, force, ErrStat, ErrMsg, &
                                             jac_rel, jac_tan)
    !! Quadratic drag of a line attachment: f = rho [CdA |w_n| w_n + CdAx |w_t| w_t] / 2 with
    !! w the relative fluid velocity split normal and axial to the line tangent (any non-zero
    !! vector; normalised internally). Optional Jacobians d f/d w and d f/d tangent, as a pair.
    !! Shared by the dynamic residual and the deck static solve, so both use one law.
    REAL(wp), INTENT(IN) :: rel_velocity(3), tangent(3), rho, cda, cdax
    REAL(wp), INTENT(OUT) :: force(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(OUT), OPTIONAL :: jac_rel(3, 3), jac_tan(3, 3)
    REAL(wp), PARAMETER :: PI_A = 3.14159265358979323846_wp
    ! The Morison per-length primitive with unit diameter is 1/2 rho Cdn |w_n| w_n +
    ! 1/2 rho pi Cdt |w_t| w_t, so Cdn = CdA and Cdt = CdAx/pi.
    IF (PRESENT(jac_rel) .AND. PRESENT(jac_tan)) THEN
      CALL CD_Morison_Drag_Per_Length_Jac(rel_velocity, tangent, rho, CD_ONE, cda, cdax/PI_A, force, &
                                          jac_rel, jac_tan, ErrStat, ErrMsg)
    ELSE
      CALL CD_Morison_Drag_Per_Length_Jac(rel_velocity, tangent, rho, CD_ONE, cda, cdax/PI_A, force, &
                                          ErrStat=ErrStat, ErrMsg=ErrMsg)
    END IF
    IF (ErrStat /= CD_HYDRO_OK) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'attachment drag: '//TRIM(ErrMsg)
    ELSE
      ErrStat = CD_HCDYN_OK
    END IF
  END SUBROUTINE CD_HermiteCable_Attachment_Drag

  SUBROUTINE assemble_attachments(model, q_cfg, v_cfg, want_residual, want_tangent, R, Kb, Kvb, ErrStat, ErrMsg)
    !! Attachment loads into R = f_int - f_ext, with the velocity tangent dR/dv and the
    !! drag's dependence on the nodal line tangent in dR/dq. The held fluid field is
    !! constant over the step.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: q_cfg(:), v_cfg(:)
    LOGICAL, INTENT(IN) :: want_residual, want_tangent
    REAL(wp), INTENT(INOUT) :: R(:), Kb(:, :), Kvb(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: k, base, a, b
    LOGICAL :: has_jac
    REAL(wp) :: rk(3), jrel(3, 3), jtan(3, 3)
    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    DO k = 1, SIZE(model%att_node)
      base = 6*(model%att_node(k) - 1)
      CALL attachment_residual(model, k, q_cfg, v_cfg, want_tangent, rk, has_jac, jrel, jtan, ErrStat, ErrMsg)
      IF (ErrStat /= CD_HCDYN_OK) RETURN
      IF (want_residual) R(base + 1:base + 3) = R(base + 1:base + 3) + rk
      IF (want_tangent .AND. has_jac) THEN
        ! w = u - v, so dR/dv = +df/dw; the tangent DOFs m enter through df/dm
        DO b = 1, 3
          DO a = 1, 3
            Kvb(KL_D + KU_D + 1 + a - b, base + b) = Kvb(KL_D + KU_D + 1 + a - b, base + b) + jrel(a, b)
            Kb(KL_D + KU_D + 1 + a - b - 3, base + b + 3) = Kb(KL_D + KU_D + 1 + a - b - 3, base + b + 3) - &
                                                            jtan(a, b)
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE assemble_attachments

  SUBROUTINE attachment_residual(model, k, q_cfg, v_cfg, want_jac, rk, has_jac, jrel, jtan, ErrStat, ErrMsg)
    !! Contribution rk of attachment k to R = f_int - f_ext at the translations of its node:
    !! weight, buoyancy, Froude-Krylov force and drag (inertia excluded). With want_jac, the
    !! drag Jacobians d f/d w and d f/d tangent are returned when has_jac is set.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: k
    REAL(wp), INTENT(IN) :: q_cfg(:), v_cfg(:)
    LOGICAL, INTENT(IN) :: want_jac
    REAL(wp), INTENT(OUT) :: rk(3), jrel(3, 3), jtan(3, 3)
    LOGICAL, INTENT(OUT) :: has_jac
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: j, base
    REAL(wp) :: u(3), ud(3), wrel(3), wl, f(3)
    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    rk = CD_ZERO; jrel = CD_ZERO; jtan = CD_ZERO; has_jac = .FALSE.
    j = model%att_node(k)
    base = 6*(j - 1)
    u = CD_ZERO; ud = CD_ZERO; wl = CD_ZERO
    IF (model%has_drag) THEN
      wl = model%waterline_z
      u = model%current
    ELSE IF (model%has_am) THEN
      wl = model%am_wl
    END IF
    IF (model%has_held_field) THEN
      u = u + model%hf_u(:, j); ud = model%hf_ud(:, j); wl = model%hf_wl(j)
    END IF
    IF (q_cfg(base + 3) > wl) THEN
      ! above the free surface: weight only, no buoyancy or fluid load
      rk(3) = model%att_w(k) + model%att_b(k)
      RETURN
    END IF
    rk(3) = model%att_w(k)
    rk = rk - model%att_fk(k)*ud
    IF (.NOT. (model%att_cda(k) > CD_ZERO .OR. model%att_cdax(k) > CD_ZERO)) RETURN
    wrel = u - v_cfg(base + 1:base + 3)
    IF (want_jac) THEN
      CALL CD_HermiteCable_Attachment_Drag(wrel, q_cfg(base + 4:base + 6), model%att_rho, model%att_cda(k), &
                                           model%att_cdax(k), f, ErrStat, ErrMsg, jrel, jtan)
      has_jac = .TRUE.
    ELSE
      CALL CD_HermiteCable_Attachment_Drag(wrel, q_cfg(base + 4:base + 6), model%att_rho, model%att_cda(k), &
                                           model%att_cdax(k), f, ErrStat, ErrMsg)
    END IF
    IF (ErrStat /= CD_HCDYN_OK) RETURN
    rk = rk - f
  END SUBROUTINE attachment_residual

  SUBROUTINE CD_HermiteCable_Dyn_Recompute_Acceleration(model, ErrStat, ErrMsg, pres_dofs, pres_a)
    !! Public lifecycle hook for a non-stepping boundary-state commit. Recompute the
    !! free-DOF acceleration from the complete installed load set at the model's current
    !! q/v/t. Optional prescribed fixed-DOF accelerations enter the consistent-mass RHS,
    !! including their off-diagonal coupling to free DOFs. Failure is atomic with respect
    !! to the previously committed acceleration.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: pres_dofs(:)
    REAL(wp), INTENT(IN), OPTIONAL :: pres_a(:)
    REAL(wp), ALLOCATABLE :: a_old(:)

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Recompute_Acceleration: model not initialized'; RETURN
    END IF
    a_old = model%a
    IF (PRESENT(pres_dofs) .AND. PRESENT(pres_a)) THEN
      CALL solve_consistent_acceleration(model, ErrStat, ErrMsg, pres_dofs=pres_dofs, pres_a=pres_a)
    ELSE IF (PRESENT(pres_dofs) .OR. PRESENT(pres_a)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Recompute_Acceleration: pres_dofs and pres_a must be supplied together'
      RETURN
    ELSE
      CALL solve_consistent_acceleration(model, ErrStat, ErrMsg)
    END IF
    IF (ErrStat /= CD_HCDYN_OK) THEN
      model%a = a_old
      ErrMsg = 'CD_HermiteCable_Dyn_Recompute_Acceleration: '//TRIM(ErrMsg)
    END IF
  END SUBROUTINE CD_HermiteCable_Dyn_Recompute_Acceleration

  SUBROUTINE CD_HermiteCable_Dyn_Set_ForceBlend(model, enabled, ErrStat, ErrMsg)
    !! Select the generalised-alpha blend. enabled = .TRUE. (the default): the residual
    !! blends the forces, G = M a_am + (1-af) r(q_{n+1}, v_{n+1}, t_{n+1}) + af r(q_n, v_n, t_n),
    !! with r_n reused from the previous converged step. enabled = .FALSE.: the forces are
    !! evaluated at the blended state (q_af, v_af, t_af), the classical form of Chung and
    !! Hulbert, whose blended Hermite tangents are shortened when the line rotates.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    LOGICAL, INTENT(IN) :: enabled
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_ForceBlend: model not initialised'
      RETURN
    END IF
    model%force_blend = enabled
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
  END SUBROUTINE CD_HermiteCable_Dyn_Set_ForceBlend

  SUBROUTINE CD_HermiteCable_Dyn_Set_ModifiedNewton(model, enabled, ErrStat, ErrMsg)
    !! Enable/disable modified Newton (within-step tangent reuse) on the (already
    !! initialised) dynamic model: the effective tangent is assembled and factored ONCE
    !! per step at the predictor and the factorization is reused for every subsequent
    !! Newton iteration of that step; a line search that cannot improve under the stale
    !! direction refreshes the tangent once before the safe damped fallback. Disabled
    !! (the default), a step factors a fresh tangent every iteration after any carried
    !! factor (cross-step reuse) is dropped. The committed state still
    !! satisfies the SAME convergence tolerance either way -- reuse changes the iteration
    !! path, not the equations or the convergence criterion.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    LOGICAL, INTENT(IN) :: enabled
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_ModifiedNewton: model not initialised'; RETURN
    END IF
    model%modified_newton = enabled
  END SUBROUTINE CD_HermiteCable_Dyn_Set_ModifiedNewton

  SUBROUTINE CD_HermiteCable_Dyn_Set_TangentReuse(model, enabled, ErrStat, ErrMsg)
    !! Enable (the default) or disable cross-step factor reuse (see the type note): a step
    !! that continues from the last committed state at the same dt starts Newton on that
    !! step's factored tangent and refreshes it when the contraction weakens. Every step
    !! still converges to the same tolerance. A host that restarts from a checkpoint and
    !! must reproduce the uninterrupted run disables it: the factor is not in the state.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    LOGICAL, INTENT(IN) :: enabled
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_TangentReuse: model not initialised'; RETURN
    END IF
    model%tangent_reuse = enabled
    model%tr_valid = .FALSE.
  END SUBROUTINE CD_HermiteCable_Dyn_Set_TangentReuse

  SUBROUTINE CD_HermiteCable_Dyn_Set_Tensile_Safety(model, enabled, ErrStat, ErrMsg)
    !! Require a physically tensile finite-bending response. The element-mean axial force,
    !! averaged over three elements, is checked after nonlinear convergence and before commit.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    LOGICAL, INTENT(IN) :: enabled
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    CALL CD_HermiteCable_Dyn_Set_Tensile_Monitor(model, &
                                                 MERGE(CD_HCDYN_TENSILE_ERROR, &
                                                       CD_HCDYN_TENSILE_OFF, enabled), &
                                                 ErrStat, ErrMsg, &
                                                 strain_tolerance=CD_HCDYN_AXIAL_STRAIN_TOL)
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Tensile_Safety

  SUBROUTINE CD_HermiteCable_Dyn_Set_Tensile_Monitor(model, mode, ErrStat, ErrMsg, strain_tolerance)
    !! Configure the axial-force audit of the converged step (the element-mean axial force,
    !! averaged over three elements; see CD_HermiteCable_Dyn_Step). ERROR is the
    !! tensile_safety=True hard gate. WARN commits the state and records one event per accepted
    !! integration step, allowing long runs to remain auditable without suppressing results.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER, INTENT(IN) :: mode
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: strain_tolerance
    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Tensile_Monitor: model not initialised'; RETURN
    END IF
    IF (mode < CD_HCDYN_TENSILE_OFF .OR. mode > CD_HCDYN_TENSILE_ERROR) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Tensile_Monitor: invalid mode'; RETURN
    END IF
    IF (PRESENT(strain_tolerance)) THEN
      IF (.NOT. CD_Is_Finite(strain_tolerance) .OR. strain_tolerance < CD_ZERO) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Tensile_Monitor: strain_tolerance must be finite and non-negative'
        RETURN
      END IF
      model%tensile_strain_tolerance = strain_tolerance
    END IF
    model%tensile_mode = mode
    model%require_tensile = mode == CD_HCDYN_TENSILE_ERROR
    model%tensile_event_count = 0
    model%tensile_worst_element = 0
    model%tensile_worst_force = CD_ZERO
    model%tensile_worst_threshold = CD_ZERO
    model%tensile_worst_xi = CD_ZERO
    model%tensile_worst_time = CD_ZERO
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Tensile_Monitor

  SUBROUTINE CD_HermiteCable_Dyn_Get_Tensile_Diagnostics(model, event_count, worst_force, &
                                                         worst_threshold, worst_element, worst_xi, &
                                                         worst_time, ErrStat, ErrMsg)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    INTEGER, INTENT(OUT) :: event_count, worst_element, ErrStat
    REAL(wp), INTENT(OUT) :: worst_force, worst_threshold, worst_xi, worst_time
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    event_count = 0; worst_element = 0
    worst_force = CD_ZERO; worst_threshold = CD_ZERO; worst_xi = CD_ZERO; worst_time = CD_ZERO
    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Get_Tensile_Diagnostics: model not initialised'; RETURN
    END IF
    event_count = model%tensile_event_count
    worst_force = model%tensile_worst_force
    worst_threshold = model%tensile_worst_threshold
    worst_element = model%tensile_worst_element
    worst_xi = model%tensile_worst_xi
    worst_time = model%tensile_worst_time
  END SUBROUTINE CD_HermiteCable_Dyn_Get_Tensile_Diagnostics

  SUBROUTINE CD_HermiteCable_Dyn_Get_Recovery_Diagnostics(model, event_count, max_substeps_used, &
                                                          ErrStat, ErrMsg)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    INTEGER, INTENT(OUT) :: event_count, max_substeps_used, ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    event_count = 0; max_substeps_used = 0; ErrStat = CD_HCDYN_OK; ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Get_Recovery_Diagnostics: model not initialised'; RETURN
    END IF
    event_count = model%recovery_event_count
    max_substeps_used = model%recovery_max_substeps_used
  END SUBROUTINE CD_HermiteCable_Dyn_Get_Recovery_Diagnostics

  SUBROUTINE CD_HermiteCable_Dyn_Set_AdaptiveNewton(model, enabled, ErrStat, ErrMsg, warmup, threshold)
    !! Enable/disable ADAPTIVE modified Newton: the solver runs the first `warmup` steps in
    !! full Newton, then latches tangent reuse ON for the rest of the run iff the warmup
    !! averaged more than `threshold` Newton iterations/step -- so reuse is used only where it
    !! pays (harder/larger systems) and full Newton keeps the easy cases cheap. Resets the
    !! measurement counters. Convergence tolerance is unchanged (reuse changes only the
    !! iteration path). warmup defaults to 4, threshold to 5.0. Disabling reverts to the
    !! explicit modified_newton flag.
    !!
    !! CAVEAT (off by default): the warmup-then-latch heuristic is NOT a reliable
    !! win and is not wired into any deck. Forced modified Newton (Set_ModifiedNewton) is the
    !! working reuse path -- it starts reuse at REST where the reference tangent is fresh.
    !! The adaptive latch instead switches ON during the harder transient, where the reference
    !! tangent is immediately stale, and can need more solves per step than both full and
    !! forced modified Newton; the short warmup also under-measures, so a genuinely stiff run
    !! may never latch. A correct auto-selector needs a rebuild-reference-on-activation or
    !! per-step disable-if-not-helping scheme, not a warmup threshold. Prefer
    !! Set_ModifiedNewton; use this only for experiments.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    LOGICAL, INTENT(IN) :: enabled
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: warmup
    REAL(wp), INTENT(IN), OPTIONAL :: threshold
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_AdaptiveNewton: model not initialised'; RETURN
    END IF
    ! Validate BOTH optional parameters before mutating the model: a rejected call must not
    ! leave a partial change (e.g. a valid warmup applied while an invalid threshold is rejected).
    IF (PRESENT(warmup)) THEN
      IF (warmup < 1) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_AdaptiveNewton: warmup must be >= 1'; RETURN
      END IF
    END IF
    IF (PRESENT(threshold)) THEN
      IF (.NOT. (threshold > CD_ZERO)) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_AdaptiveNewton: threshold must be > 0'; RETURN
      END IF
    END IF
    IF (PRESENT(warmup)) model%na_warmup = warmup
    IF (PRESENT(threshold)) model%na_threshold = threshold
    model%newton_adaptive = enabled
    model%na_steps_done = 0
    model%na_iter_sum = CD_ZERO
    model%na_mn_active = .FALSE.
  END SUBROUTINE CD_HermiteCable_Dyn_Set_AdaptiveNewton

  SUBROUTINE CD_HermiteCable_Dyn_Set_EndConnection(model, k_rot, d0, ErrStat, ErrMsg, connection_mode)
    !! Attach pinned, finite-stiffness, or exact rigid bending connections to
    !! an initialised dynamic model.
    !!
    !! Finite mode is an isotropic joint carrying rotational stiffness about a
    !! preferred no-moment direction. It enters as a configuration-dependent
    !! generalized load with a consistent symmetric tangent. Rigid mode instead
    !! eliminates the two transverse tangent coordinates in an orthonormal local
    !! basis; it is not a large-stiffness approximation.
    !!
    !! The spring restrains only the tangent DIRECTION: its generalized force satisfies
    !! f . d == 0 identically, so it never acts on |m|, which carries axial stretch.
    !! Clamping that component instead would impose an inextensibility over-constraint.
    !!
    !! k_rot(2)  : rotational stiffness [N m/rad] at end A (node 1) and end B (node nn).
    !!             A zero entry leaves that end pinned. Negative and non-finite values
    !!             are rejected. Rigid mode requires zero here.
    !! d0(3,2)   : preferred no-moment directions for those ends; each must be non-null and
    !!             is normalised here. Required for every non-pinned end.
    !! connection_mode(2), optional: Pinned, Finite, or Rigid. If omitted,
    !!             zero/positive stiffness selects Pinned/Finite.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: k_rot(2), d0(3, 2)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: connection_mode(2)
    INTEGER :: iend, base, mode_new(2), mode_old(2)
    REAL(wp) :: dn, d_new(3, 2), f_check(3), k_check(3, 3), basis(3, 3), projected(3)
    REAL(wp), ALLOCATABLE :: q_new(:), v_new(:), q_old(:), v_old(:), a_old(:)
    REAL(wp) :: k_old(2), d_old(3, 2)
    LOGICAL :: had_conn, snap_was_valid
    LOGICAL, ALLOCATABLE :: solve_mask_old(:)
    INTEGER :: es
    CHARACTER(200) :: em

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    ! A load-configuration change invalidates the cached committed force F_n.
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: model is not initialised'
      RETURN
    END IF

    IF (.NOT. CD_All_Finite(k_rot) .OR. ANY(k_rot < CD_ZERO)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: stiffnesses must be finite and non-negative'
      RETURN
    END IF

    ! Validate and normalise BEFORE touching the model, so a rejected call leaves no trace.
    d_new = CD_ZERO
    IF (PRESENT(connection_mode)) THEN
      mode_new = connection_mode
    ELSE
      WHERE (k_rot > CD_ZERO)
        mode_new = CD_ENDCONN_FINITE
      ELSEWHERE
        mode_new = CD_ENDCONN_PINNED
      END WHERE
    END IF
    ALLOCATE (q_new(model%ndof), v_new(model%ndof), q_old(model%ndof), &
              v_old(model%ndof), a_old(model%ndof), solve_mask_old(model%ndof))
    q_new = model%q
    v_new = model%v
    DO iend = 1, 2
      IF (mode_new(iend) < CD_ENDCONN_PINNED .OR. mode_new(iend) > CD_ENDCONN_RIGID) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: mode must be Pinned, Finite, or Rigid'
        RETURN
      END IF
      IF (mode_new(iend) == CD_ENDCONN_PINNED) THEN
        IF (k_rot(iend) > CD_ZERO) THEN
          ErrStat = CD_HCDYN_BADINPUT
          ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: Pinned mode requires zero stiffness'
          RETURN
        END IF
        CYCLE
      ELSE IF (mode_new(iend) == CD_ENDCONN_FINITE) THEN
        IF (k_rot(iend) <= CD_ZERO) THEN
          ErrStat = CD_HCDYN_BADINPUT
          ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: Finite mode requires positive stiffness'
          RETURN
        END IF
      ELSE IF (k_rot(iend) > CD_ZERO) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: Rigid mode is exact and requires zero stiffness'
        RETURN
      END IF
      IF (.NOT. CD_All_Finite(d0(:, iend))) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: preferred directions must be finite'
        RETURN
      END IF
      dn = SQRT(DOT_PRODUCT(d0(:, iend), d0(:, iend)))
      IF (.NOT. CD_Is_Finite(dn) .OR. dn <= SQRT(TINY(CD_ONE))) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: preferred direction must be non-null '// &
                 'for a non-pinned end'
        RETURN
      END IF
      d_new(:, iend) = d0(:, iend)/dn
      IF (iend == 1) THEN
        base = 3
      ELSE
        base = 6*(model%nn - 1) + 3
      END IF
      IF (mode_new(iend) == CD_ENDCONN_FINITE) THEN
        CALL CD_EndConn_Spring(model%q(base + 1:base + 3), d_new(:, iend), k_rot(iend), &
                               f_check, k_check, es, em)
        IF (es /= CD_ENDCONN_OK) THEN
          ErrStat = CD_HCDYN_BADINPUT
          ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: '//TRIM(em)
          RETURN
        END IF
      ELSE
        IF (ANY(.NOT. model%freemask(base + 1:base + 3))) THEN
          ErrStat = CD_HCDYN_BADINPUT
          ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: rigid tangent DOFs must not also be fixed'
          RETURN
        END IF
        CALL CD_EndConn_Basis(d_new(:, iend), basis, es, em)
        IF (es == CD_ENDCONN_OK) &
          CALL CD_EndConn_Project(q_new(base + 1:base + 3), d_new(:, iend), projected, es, em)
        IF (es /= CD_ENDCONN_OK) THEN
          ErrStat = CD_HCDYN_BADINPUT
          ErrMsg = 'CD_HermiteCable_Dyn_Set_EndConnection: '//TRIM(em)
          RETURN
        END IF
        q_new(base + 1:base + 3) = projected
        v_new(base + 1:base + 3) = DOT_PRODUCT(v_new(base + 1:base + 3), basis(:, 1))*basis(:, 1)
      END IF
    END DO

    ! Transactional, as Set_Drag is: the acceleration refresh below can fail after the
    ! parameters are installed, and the model must then be restored rather than left in a
    ! partially applied state.
    had_conn = model%has_endconn
    k_old = model%endconn_k
    d_old = model%endconn_d0
    mode_old = model%endconn_mode
    q_old = model%q
    v_old = model%v
    a_old = model%a
    solve_mask_old = model%solve_mask
    snap_was_valid = model%snap_valid

    model%endconn_d0 = d_new
    model%endconn_k = k_rot
    model%endconn_mode = mode_new
    model%has_endconn = ANY(mode_new /= CD_ENDCONN_PINNED)
    model%q = q_new
    model%v = v_new
    model%a = a_old
    model%solve_mask = model%freemask
    DO iend = 1, 2
      IF (mode_new(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (iend == 1) THEN
        base = 3
      ELSE
        base = 6*(model%nn - 1) + 3
      END IF
      model%a(base + 1:base + 3) = DOT_PRODUCT(model%a(base + 1:base + 3), &
                                               d_new(:, iend))*d_new(:, iend)
      model%solve_mask(base + 2:base + 3) = .FALSE.
    END DO
    model%snap_valid = .FALSE.

    ! Refresh the consistent acceleration now that the load set includes the connection.
    ! Init solved M a0 = f_ext - f_int without it; leaving a0 stale would launch the first
    ! step with an acceleration that omits the connection moment, so a state that is a true
    ! equilibrium WITH the connection would still move.
    CALL solve_consistent_acceleration(model, ErrStat, ErrMsg)
    IF (ErrStat /= CD_HCDYN_OK) THEN
      model%a = a_old
      model%q = q_old
      model%v = v_old
      model%endconn_k = k_old
      model%endconn_d0 = d_old
      model%endconn_mode = mode_old
      model%solve_mask = solve_mask_old
      model%has_endconn = had_conn
      model%snap_valid = snap_was_valid
    END IF
  END SUBROUTINE CD_HermiteCable_Dyn_Set_EndConnection

  SUBROUTINE prepare_endconn_directions(model, requested, prepared, ErrStat, ErrMsg)
    !! Validate and normalise a target set of preferred directions without mutating
    !! the committed model. Every non-pinned connection consumes a direction.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: requested(3, 2)
    REAL(wp), INTENT(OUT) :: prepared(3, 2)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: iend
    REAL(wp) :: dn

    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    prepared = model%endconn_d0
    IF (.NOT. model%has_endconn) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'moving end-connection directions require a configured connection'
      RETURN
    END IF
    DO iend = 1, 2
      IF (model%endconn_mode(iend) == CD_ENDCONN_PINNED) CYCLE
      IF (.NOT. CD_All_Finite(requested(:, iend))) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'moving end-connection directions must be finite'
        RETURN
      END IF
      dn = SQRT(DOT_PRODUCT(requested(:, iend), requested(:, iend)))
      IF (.NOT. CD_Is_Finite(dn) .OR. dn <= SQRT(TINY(CD_ONE))) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'moving end-connection direction must be non-null at each connected end'
        RETURN
      END IF
      prepared(:, iend) = requested(:, iend)/dn
    END DO
  END SUBROUTINE prepare_endconn_directions

  SUBROUTINE interpolate_endconn_directions(model, d0, d1, fraction, d, ErrStat, ErrMsg)
    !! Shortest-arc interpolation of the preferred direction over a coupling
    !! interval. This supplies the direction at the generalised-alpha load station
    !! and at each recovery substep without inventing a stiffness-dependent penalty.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: d0(3, 2), d1(3, 2), fraction
    REAL(wp), INTENT(OUT) :: d(3, 2)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: iend
    REAL(wp) :: dot01, theta, sintheta, dn

    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    d = d0
    IF (.NOT. CD_Is_Finite(fraction) .OR. fraction < CD_ZERO .OR. fraction > CD_ONE) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'end-connection interpolation fraction must lie in [0,1]'
      RETURN
    END IF
    DO iend = 1, 2
      IF (model%endconn_mode(iend) == CD_ENDCONN_PINNED) CYCLE
      dot01 = MAX(-CD_ONE, MIN(CD_ONE, DOT_PRODUCT(d0(:, iend), d1(:, iend))))
      IF (dot01 <= -CD_ONE + 64.0_wp*EPSILON(CD_ONE)) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'preferred direction rotates by 180 degrees in one step; subdivide the host motion'
        RETURN
      ELSE IF (dot01 >= CD_ONE - 64.0_wp*EPSILON(CD_ONE)) THEN
        d(:, iend) = (CD_ONE - fraction)*d0(:, iend) + fraction*d1(:, iend)
      ELSE
        theta = ACOS(dot01)
        sintheta = SIN(theta)
        d(:, iend) = SIN((CD_ONE - fraction)*theta)/sintheta*d0(:, iend) + &
                     SIN(fraction*theta)/sintheta*d1(:, iend)
      END IF
      dn = SQRT(DOT_PRODUCT(d(:, iend), d(:, iend)))
      IF (.NOT. CD_Is_Finite(dn) .OR. dn <= SQRT(TINY(CD_ONE))) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'preferred-direction interpolation is singular'
        RETURN
      END IF
      d(:, iend) = d(:, iend)/dn
    END DO
  END SUBROUTINE interpolate_endconn_directions

  SUBROUTINE prepare_endconn_derivatives(model, directions, rates, accelerations, ErrStat, ErrMsg, &
                                         requested_rates, requested_accelerations)
    !! Validate first and second derivatives of rigid unit preferred directions. For a
    !! unit vector d, d.d_dot = 0 and d.d_ddot = -d_dot.d_dot. Enforcing these
    !! identities at the API boundary prevents a host from silently changing the
    !! tangent magnitude through malformed angular kinematics.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: directions(3, 2)
    REAL(wp), INTENT(OUT) :: rates(3, 2), accelerations(3, 2)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: requested_rates(3, 2), requested_accelerations(3, 2)
    INTEGER :: iend
    REAL(wp) :: orthogonality, normal_identity, scale

    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    rates = CD_ZERO
    accelerations = CD_ZERO
    IF (PRESENT(requested_rates) .NEQV. PRESENT(requested_accelerations)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'direction rate and acceleration must be supplied together'
      RETURN
    END IF
    IF (.NOT. PRESENT(requested_rates)) RETURN
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (.NOT. (CD_All_Finite(requested_rates(:, iend)) .AND. &
                 CD_All_Finite(requested_accelerations(:, iend)))) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'rigid direction derivatives must be finite'
        RETURN
      END IF
      rates(:, iend) = requested_rates(:, iend)
      accelerations(:, iend) = requested_accelerations(:, iend)
      orthogonality = DOT_PRODUCT(directions(:, iend), rates(:, iend))
      normal_identity = DOT_PRODUCT(directions(:, iend), accelerations(:, iend)) + &
                        DOT_PRODUCT(rates(:, iend), rates(:, iend))
      scale = MAX(CD_ONE, SQRT(DOT_PRODUCT(rates(:, iend), rates(:, iend))), &
                  SQRT(DOT_PRODUCT(accelerations(:, iend), accelerations(:, iend))))
      IF (ABS(orthogonality) > 1.0e-9_wp*scale .OR. ABS(normal_identity) > 1.0e-9_wp*scale*scale) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'rigid direction derivatives violate unit-vector kinematics'
        RETURN
      END IF
    END DO
  END SUBROUTINE prepare_endconn_derivatives

  SUBROUTINE rigid_tangent_scalar_state(model, directions, q, v, a, magnitude, rate, acceleration)
    !! Recover n, n_dot, and n_ddot from the global tangent state m = n d.
    !! The transverse velocity determines d_dot; d.d_ddot = -|d_dot|^2 then
    !! separates the material magnitude acceleration from m_ddot.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: directions(3, 2), q(:), v(:), a(:)
    REAL(wp), INTENT(OUT) :: magnitude(2), rate(2), acceleration(2)
    INTEGER :: iend, base
    REAL(wp) :: direction_rate(3)
    magnitude = CD_ZERO
    rate = CD_ZERO
    acceleration = CD_ZERO
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (iend == 1) THEN
        base = 3
      ELSE
        base = 6*(model%nn - 1) + 3
      END IF
      magnitude(iend) = DOT_PRODUCT(q(base + 1:base + 3), directions(:, iend))
      rate(iend) = DOT_PRODUCT(v(base + 1:base + 3), directions(:, iend))
      direction_rate = (v(base + 1:base + 3) - rate(iend)*directions(:, iend))/magnitude(iend)
      acceleration(iend) = DOT_PRODUCT(a(base + 1:base + 3), directions(:, iend)) + &
                           magnitude(iend)*DOT_PRODUCT(direction_rate, direction_rate)
    END DO
  END SUBROUTINE rigid_tangent_scalar_state

  SUBROUTINE rigid_direction_derivatives_from_state(model, directions, q, v, a, rates, accelerations)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: directions(3, 2), q(:), v(:), a(:)
    REAL(wp), INTENT(OUT) :: rates(3, 2), accelerations(3, 2)
    REAL(wp) :: magnitude(2), magnitude_rate(2), magnitude_acceleration(2)
    INTEGER :: iend, base
    rates = CD_ZERO
    accelerations = CD_ZERO
    CALL rigid_tangent_scalar_state(model, directions, q, v, a, magnitude, magnitude_rate, magnitude_acceleration)
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (iend == 1) THEN
        base = 3
      ELSE
        base = 6*(model%nn - 1) + 3
      END IF
      rates(:, iend) = (v(base + 1:base + 3) - magnitude_rate(iend)*directions(:, iend))/magnitude(iend)
      accelerations(:, iend) = (a(base + 1:base + 3) - &
                                magnitude_acceleration(iend)*directions(:, iend) - &
                                2.0_wp*magnitude_rate(iend)*rates(:, iend))/magnitude(iend)
    END DO
  END SUBROUTINE rigid_direction_derivatives_from_state

  SUBROUTINE CD_HermiteCable_Dyn_Set_Drag(model, rho_w, diam, cdn, cdt, waterline_z, current, ErrStat, ErrMsg, &
                                          gravity)
    !! Enable per-element Morison quadratic drag on the (already initialised) dynamic model. The
    !! drag opposes the cable's velocity relative to a uniform fluid velocity `current`; a Gauss
    !! point above waterline_z is dry (no drag). This is the damping that makes a driven dynamic
    !! cable physical (an undamped cable builds up under sustained forcing).
    !!
    !! The same configuration defines the free surface for the self-weight: the part of the line
    !! above the surface displaces no water, so its buoyancy rho_w*gravity*pi*diam^2/4 per
    !! reference length is restored on top of the submerged weight w (exact dry-interval
    !! integration with its waterline-crossing tangent; CD_HermiteCable_Dry_Buoyancy).
    !!
    !! rho_w        : fluid mass density (kg/m^3), > 0.
    !! diam(ne)     : per-element hydrodynamic (drag and displacement) diameter, > 0.
    !! cdn, cdt(ne) : per-element normal / tangential drag coefficients, >= 0.
    !! waterline_z  : still-water level; a Gauss point with z above it is dry.
    !! current(3)   : uniform fluid velocity (0 = still water, drag purely damps the cable's motion).
    !! gravity      : optional gravitational acceleration for the dry-part buoyancy (m/s^2, > 0;
    !!                default 9.80665). Pass the value the submerged weight w was formed with.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: rho_w, diam(:), cdn(:), cdt(:), waterline_z, current(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: gravity

    REAL(wp), ALLOCATABLE :: diam_old(:), cdn_old(:), cdt_old(:), a_old(:)
    REAL(wp) :: rho_old, wl_old, cur_old(3), grav_old, grav_new
    LOGICAL :: had_drag

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    ! A load-configuration change invalidates the cached committed force F_n.
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Set_Drag: model not initialised'; RETURN
    END IF
    IF (SIZE(diam) /= model%ne .OR. SIZE(cdn) /= model%ne .OR. SIZE(cdt) /= model%ne) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Drag: diam/cdn/cdt must have length ne'; RETURN
    END IF
    IF (.NOT. CD_Is_Finite(rho_w) .OR. rho_w <= CD_ZERO .OR. .NOT. CD_Is_Finite(waterline_z) .OR. &
        .NOT. CD_All_Finite(current)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Drag: need rho_w>0, finite waterline/current'; RETURN
    END IF
    IF (.NOT. CD_All_Finite(diam) .OR. ANY(diam <= CD_ZERO) .OR. &
        .NOT. CD_All_Finite(cdn) .OR. ANY(cdn < CD_ZERO) .OR. &
        .NOT. CD_All_Finite(cdt) .OR. ANY(cdt < CD_ZERO)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Drag: need diam>0, cdn>=0, cdt>=0 (finite)'; RETURN
    END IF
    grav_new = 9.80665_wp
    IF (PRESENT(gravity)) grav_new = gravity
    IF (.NOT. CD_Is_Finite(grav_new) .OR. grav_new <= CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Drag: gravity must be finite and positive'; RETURN
    END IF
    ! One free surface: once waves are on, every hydro configuration must share the waterline
    ! (Set_Waves enforces this at enable time; this closes the reversed call order).
    IF (model%has_wave .AND. model%has_am) THEN
      IF (ABS(waterline_z - model%am_wl) > CD_ZERO) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Drag: waterline disagrees with the added-mass config while waves are on'
        RETURN
      END IF
    END IF

    ! Snapshot for rollback: Set_Drag is TRANSACTIONAL. The acceleration refresh below can fail
    ! after the parameters are installed (e.g. a finite but overflowing current makes the rest drag
    ! non-finite); the model must then be restored to its previous valid drag state and previous
    ! consistent acceleration -- not left with a partially applied configuration.
    had_drag = model%has_drag
    a_old = model%a
    rho_old = model%rho_w; wl_old = model%waterline_z; cur_old = model%current
    grav_old = model%gravity
    ! Allocate the snapshots unconditionally (zero when there was no prior drag) so every restore
    ! path reads a defined array.
    ALLOCATE (diam_old(model%ne), cdn_old(model%ne), cdt_old(model%ne))
    IF (had_drag) THEN
      diam_old = model%diam; cdn_old = model%cdn; cdt_old = model%cdt
    ELSE
      diam_old = CD_ZERO; cdn_old = CD_ZERO; cdt_old = CD_ZERO
    END IF

    IF (ALLOCATED(model%diam)) DEALLOCATE (model%diam)
    IF (ALLOCATED(model%cdn)) DEALLOCATE (model%cdn)
    IF (ALLOCATED(model%cdt)) DEALLOCATE (model%cdt)
    ALLOCATE (model%diam(model%ne), model%cdn(model%ne), model%cdt(model%ne))
    model%diam = diam; model%cdn = cdn; model%cdt = cdt
    model%rho_w = rho_w; model%waterline_z = waterline_z; model%current = current
    model%gravity = grav_new
    model%has_drag = .TRUE.

    ! Refresh the consistent acceleration now that the load set includes drag: at v = 0 a current
    ! exerts a non-zero rest drag, so the predictor for the first step must see it.
    CALL solve_consistent_acceleration(model, ErrStat, ErrMsg)
    IF (ErrStat /= CD_HCDYN_OK) THEN
      ! Roll back to the pre-call state: previous acceleration, previous drag configuration.
      model%a = a_old
      model%rho_w = rho_old; model%waterline_z = wl_old; model%current = cur_old
      model%gravity = grav_old
      IF (had_drag) THEN
        model%diam = diam_old; model%cdn = cdn_old; model%cdt = cdt_old
      ELSE
        DEALLOCATE (model%diam, model%cdn, model%cdt)
      END IF
      model%has_drag = had_drag
    END IF
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Drag

  SUBROUTINE CD_HermiteCable_Dyn_Set_Axial_Damping(model, BA, ErrStat, ErrMsg)
    !! Install the resolved per-element axial Kelvin-Voigt coefficients [N s].
    !! The constitutive viscous resultant is N_v = BA*strain_rate, where strain_rate
    !! follows the deformed centreline and material coordinate. An all-zero array
    !! disables the contribution exactly. The setter is transactional and refreshes
    !! the consistent acceleration because callers may replace the damping after the
    !! cable has acquired a non-zero velocity.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: BA(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp), ALLOCATABLE :: old_BA(:), new_BA(:), old_a(:)
    LOGICAL :: old_has
    INTEGER :: es, istat
    CHARACTER(300) :: em

    ErrStat = CD_HCDYN_OK
    ! A load-configuration change invalidates the cached committed force F_n.
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      CALL fail('model not initialised'); RETURN
    END IF
    IF (SIZE(BA) /= model%ne) THEN
      CALL fail('BA must have length ne'); RETURN
    END IF
    IF (.NOT. CD_All_Finite(BA) .OR. ANY(BA < CD_ZERO)) THEN
      CALL fail('BA must be finite and non-negative'); RETURN
    END IF

    ALLOCATE (new_BA(model%ne), old_a(model%ndof), STAT=istat)
    IF (istat /= 0) THEN
      CALL fail('workspace allocation failed'); RETURN
    END IF
    new_BA = BA
    old_a = model%a
    old_has = model%has_axial_damping
    CALL MOVE_ALLOC(model%BA, old_BA)
    IF (ANY(new_BA > CD_ZERO)) THEN
      CALL MOVE_ALLOC(new_BA, model%BA)
      model%has_axial_damping = .TRUE.
    ELSE
      model%has_axial_damping = .FALSE.
    END IF

    CALL solve_consistent_acceleration(model, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      IF (ALLOCATED(model%BA)) DEALLOCATE (model%BA)
      CALL MOVE_ALLOC(old_BA, model%BA)
      model%has_axial_damping = old_has
      model%a = old_a
      CALL fail('initial acceleration refresh failed: '//TRIM(em)); RETURN
    END IF

  CONTAINS
    SUBROUTINE fail(msg)
      CHARACTER(*), INTENT(IN) :: msg
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Axial_Damping: '//msg
    END SUBROUTINE fail
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Axial_Damping

  SUBROUTINE CD_HermiteCable_Dyn_Set_AddedMass(model, rho_w, diam, can, cat, waterline_z, ErrStat, ErrMsg)
    !! Enable per-element Morison added mass on the (already initialised) dynamic model: the
    !! submerged cable accelerates surrounding fluid, adding rho_w A (Can (I - t t^T) + Cat t t^T)
    !! per unit length to the inertia. The added mass shifts the wet natural frequencies down by
    !! sqrt(m/(m + m_a)) -- essential for the dynamic (fatigue) response of a submerged power cable.
    !! Within each step the added-mass matrix is FROZEN at the beginning-of-step configuration (an
    !! O(dt) approximation consistent with gen-alpha accuracy, exact at any at-rest fixed point), so
    !! the Newton solve keeps a constant mass. Transactional like Set_Drag: a failed refresh rolls
    !! back to the previous valid state. Independent of the drag configuration.
    !!
    !! rho_w       : fluid mass density (kg/m^3), > 0.
    !! diam(ne)    : per-element hydrodynamic diameter (m), > 0.
    !! can, cat(ne): per-element normal / tangential added-mass coefficients, >= 0.
    !! waterline_z : still-water level; a Gauss point above it carries no added mass.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: rho_w, diam(:), can(:), cat(:), waterline_z
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp), ALLOCATABLE :: d_old(:), can_old(:), cat_old(:), a_old(:)
    REAL(wp) :: rho_old, wl_old
    LOGICAL :: had_am

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    ! A load-configuration change invalidates the cached committed force F_n.
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Set_AddedMass: model not initialised'; RETURN
    END IF
    IF (SIZE(diam) /= model%ne .OR. SIZE(can) /= model%ne .OR. SIZE(cat) /= model%ne) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_AddedMass: diam/can/cat must have length ne'; RETURN
    END IF
    IF (.NOT. CD_Is_Finite(rho_w) .OR. rho_w <= CD_ZERO .OR. .NOT. CD_Is_Finite(waterline_z)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_AddedMass: need rho_w>0 and a finite waterline'; RETURN
    END IF
    IF (.NOT. CD_All_Finite(diam) .OR. ANY(diam <= CD_ZERO) .OR. &
        .NOT. CD_All_Finite(can) .OR. ANY(can < CD_ZERO) .OR. &
        .NOT. CD_All_Finite(cat) .OR. ANY(cat < CD_ZERO)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_AddedMass: need diam>0, can>=0, cat>=0 (finite)'; RETURN
    END IF
    ! One free surface: once waves are on, every hydro configuration must share the waterline
    ! (Set_Waves enforces this at enable time; this closes the reversed call order).
    IF (model%has_wave .AND. model%has_drag) THEN
      IF (ABS(waterline_z - model%waterline_z) > CD_ZERO) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_AddedMass: waterline disagrees with the drag config while waves are on'
        RETURN
      END IF
    END IF

    ! Snapshot for rollback (transactional; see Set_Drag).
    had_am = model%has_am
    a_old = model%a
    rho_old = model%am_rho; wl_old = model%am_wl
    ALLOCATE (d_old(model%ne), can_old(model%ne), cat_old(model%ne))
    IF (had_am) THEN
      d_old = model%am_diam; can_old = model%am_can; cat_old = model%am_cat
    ELSE
      d_old = CD_ZERO; can_old = CD_ZERO; cat_old = CD_ZERO
    END IF

    IF (ALLOCATED(model%am_diam)) DEALLOCATE (model%am_diam)
    IF (ALLOCATED(model%am_can)) DEALLOCATE (model%am_can)
    IF (ALLOCATED(model%am_cat)) DEALLOCATE (model%am_cat)
    ALLOCATE (model%am_diam(model%ne), model%am_can(model%ne), model%am_cat(model%ne))
    model%am_diam = diam; model%am_can = can; model%am_cat = cat
    model%am_rho = rho_w; model%am_wl = waterline_z
    model%has_am = .TRUE.

    ! Refresh the consistent acceleration with the total (structural + added) mass.
    CALL solve_consistent_acceleration(model, ErrStat, ErrMsg)
    IF (ErrStat /= CD_HCDYN_OK) THEN
      model%a = a_old
      model%am_rho = rho_old; model%am_wl = wl_old
      IF (had_am) THEN
        model%am_diam = d_old; model%am_can = can_old; model%am_cat = cat_old
      ELSE
        DEALLOCATE (model%am_diam, model%am_can, model%am_cat)
      END IF
      model%has_am = had_am
    END IF
  END SUBROUTINE CD_HermiteCable_Dyn_Set_AddedMass

  SUBROUTINE CD_HermiteCable_Dyn_Set_Waves(model, height, period, direction_deg, depth, gravity, &
                                           ErrStat, ErrMsg)
    !! Enable a regular (Airy) wave field on the (already initialised) dynamic model. The wave
    !! VELOCITY feeds the drag's relative-velocity fluid field (requires the drag configuration);
    !! the wave ACCELERATION drives the Froude-Krylov + fluid-inertia load per unit length
    !!   f = rho A [ (1 + Can) a_n + (1 + Cat) a_t ]
    !! with the added-mass configuration's rho/diameter/coefficients (requires it). At least one of
    !! the two hydro configurations must already be set, else the waves would silently do nothing.
    !! When both are set their waterline data must agree (one free surface). Kinematics use Wheeler
    !! stretching and are zero above the instantaneous surface. Transactional: a failed refresh
    !! rolls back. Waves make the wave loads time-dependent, so an at-rest state is NOT a fixed
    !! point once waves are on -- the model's simulation time starts at Init (t = 0) and advances
    !! with each step.
    !!
    !! height        : wave height H (m), > 0 (amplitude H/2).
    !! period        : wave period T (s), > 0.
    !! direction_deg : propagation direction in the x-y plane (deg from +x).
    !! depth         : water depth d (m), > 0 (seabed at waterline - d).
    !! gravity       : gravitational acceleration (m/s^2), > 0.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: height, period, direction_deg, depth, gravity
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp), PARAMETER :: PI_L = 3.14159265358979323846_wp
    REAL(wp), ALLOCATABLE :: a_old(:)
    REAL(wp) :: h_old, om_old, k_old, dep_old, dir_old, omega, k
    LOGICAL :: had_wave
    INTEGER :: es
    CHARACTER(200) :: em2

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    ! A load-configuration change invalidates the cached committed force F_n.
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Set_Waves: model not initialised'; RETURN
    END IF
    IF (model%has_held_field) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Waves: a held host fluid field is configured; an internal '// &
               'wave family on top of it would double-count the wave kinematics'; RETURN
    END IF
    IF (model%has_attach) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Waves: attachments take the held fluid field, not internal waves'
      RETURN
    END IF
    IF (.NOT. (model%has_drag .OR. model%has_am)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Waves: needs the drag and/or added-mass config first '// &
               '(waves would otherwise do nothing)'; RETURN
    END IF
    IF (model%has_drag .AND. model%has_am) THEN
      IF (ABS(model%waterline_z - model%am_wl) > CD_ZERO) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Waves: drag and added-mass waterlines disagree'; RETURN
      END IF
    END IF
    IF (.NOT. CD_Is_Finite(height) .OR. height <= CD_ZERO .OR. &
        .NOT. CD_Is_Finite(period) .OR. period <= CD_ZERO .OR. &
        .NOT. CD_Is_Finite(depth) .OR. depth <= CD_ZERO .OR. &
        .NOT. CD_Is_Finite(direction_deg) .OR. .NOT. CD_Is_Finite(gravity) .OR. gravity <= CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Waves: need height>0, period>0, depth>0, gravity>0 (finite)'; RETURN
    END IF
    omega = 2.0_wp*PI_L/period
    CALL CD_Solve_Dispersion_Wavenumber(omega, depth, gravity, k, es, em2)
    IF (es /= CD_HYDRO_OK) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Set_Waves: '//TRIM(em2); RETURN
    END IF

    ! Snapshot for rollback (transactional; see Set_Drag). A previously configured
    ! COMPONENT field is displaced by this regular field (mutual exclusivity), so its
    ! table is part of the snapshot too.
    had_wave = model%has_wave
    a_old = model%a
    h_old = model%wv_h; om_old = model%wv_om; k_old = model%wv_k
    dep_old = model%wv_depth; dir_old = model%wv_dir
    BLOCK
      INTEGER :: nc_old
      REAL(wp), ALLOCATABLE :: ampc_old(:), omc_old(:), kc_old(:), phc_old(:)
      nc_old = model%wv_ncomp
      IF (nc_old > 0) THEN
        ampc_old = model%wv_ampc; omc_old = model%wv_omc
        kc_old = model%wv_kc; phc_old = model%wv_phc
      ELSE
        ALLOCATE (ampc_old(0), omc_old(0), kc_old(0), phc_old(0))
      END IF

      model%wv_h = height; model%wv_om = omega; model%wv_k = k
      model%wv_depth = depth; model%wv_dir = direction_deg
      model%wv_ncomp = 0
      model%has_wave = .TRUE.

      ! Refresh the consistent acceleration under the wave loads at the current simulation time.
      CALL solve_consistent_acceleration(model, ErrStat, ErrMsg)
      IF (ErrStat /= CD_HCDYN_OK) THEN
        model%a = a_old
        model%wv_h = h_old; model%wv_om = om_old; model%wv_k = k_old
        model%wv_depth = dep_old; model%wv_dir = dir_old
        model%wv_ncomp = nc_old
        IF (nc_old > 0) THEN
          model%wv_ampc = ampc_old; model%wv_omc = omc_old
          model%wv_kc = kc_old; model%wv_phc = phc_old
        END IF
        model%has_wave = had_wave
        CALL refresh_wave_table(model)
      END IF
    END BLOCK
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Waves

  SUBROUTINE CD_HermiteCable_Dyn_Set_Irregular_Waves(model, amplitude, period, phase_deg, &
                                                     direction_deg, depth, gravity, ErrStat, ErrMsg)
    !! Enable a COMPONENT (irregular) wave field on the (already initialised) dynamic model: a
    !! caller-supplied long-crested linear component table -- amplitude a_i, period T_i, phase
    !! phi_i per component, one shared direction -- summed by CD_Component_Wave_Kinematics with
    !! total-eta Wheeler stretching (the OrcaFlex irregular-wave convention). The total wave
    !! VELOCITY feeds the drag's relative-velocity fluid field and the total wave ACCELERATION
    !! drives the Froude-Krylov + fluid-inertia load, exactly as the regular field above; the
    !! same hydro-configuration prerequisites and the one-free-surface rule apply. A component
    !! field displaces a previously configured regular field and vice versa (mutually
    !! exclusive). Transactional: a failed refresh rolls back.
    !!
    !! This is the deterministic component-matched irregular-sea path: an external tool
    !! synthesizes ONE spectrum discretization (e.g. JONSWAP) and feeds the SAME table to
    !! CableDyn and to the reference solver, so an irregular-sea parity case carries no
    !! stochastic seed ambiguity.
    !!
    !! amplitude(:)  : component amplitudes a_i (m), >= 0 finite, at least one > 0.
    !! period(:)     : component periods T_i (s), > 0 finite (same length).
    !! phase_deg(:)  : component phases phi_i (deg), finite (same length);
    !!                 eta_i = a_i cos(k_i x - omega_i t + phi_i).
    !! direction_deg : shared propagation direction in the x-y plane (deg from +x).
    !! depth         : water depth d (m), > 0 (seabed at waterline - d).
    !! gravity       : gravitational acceleration (m/s^2), > 0.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: amplitude(:), period(:), phase_deg(:)
    REAL(wp), INTENT(IN) :: direction_deg, depth, gravity
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp), PARAMETER :: PI_L = 3.14159265358979323846_wp
    REAL(wp), ALLOCATABLE :: a_old(:), ampc_old(:), omc_old(:), kc_old(:), phc_old(:)
    REAL(wp), ALLOCATABLE :: omc(:), kc(:), phc(:)
    REAL(wp) :: h_old, om_old, k_old, dep_old, dir_old
    LOGICAL :: had_wave
    INTEGER :: i, n, nc_old, es
    CHARACTER(200) :: em2

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    ! A load-configuration change invalidates the cached committed force F_n.
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Irregular_Waves: model not initialised'; RETURN
    END IF
    IF (model%has_held_field) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Irregular_Waves: a held host fluid field is configured; '// &
               'an internal wave family on top of it would double-count the wave kinematics'; RETURN
    END IF
    IF (model%has_attach) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Irregular_Waves: attachments take the held fluid field, '// &
               'not internal waves'
      RETURN
    END IF
    IF (.NOT. (model%has_drag .OR. model%has_am)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Irregular_Waves: needs the drag and/or added-mass '// &
               'config first (waves would otherwise do nothing)'; RETURN
    END IF
    IF (model%has_drag .AND. model%has_am) THEN
      IF (ABS(model%waterline_z - model%am_wl) > CD_ZERO) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Irregular_Waves: drag and added-mass waterlines disagree'
        RETURN
      END IF
    END IF
    n = SIZE(amplitude)
    IF (n < 1 .OR. SIZE(period) /= n .OR. SIZE(phase_deg) /= n) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Irregular_Waves: amplitude/period/phase_deg must share '// &
               'one length >= 1'; RETURN
    END IF
    IF (.NOT. (CD_All_Finite(amplitude) .AND. ALL(amplitude >= CD_ZERO) .AND. &
               ANY(amplitude > CD_ZERO) .AND. &
               CD_All_Finite(period) .AND. ALL(period > CD_ZERO) .AND. &
               CD_All_Finite(phase_deg))) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Irregular_Waves: need finite amplitudes >= 0 (one > 0), '// &
               'periods > 0, finite phases'; RETURN
    END IF
    IF (.NOT. CD_Is_Finite(depth) .OR. depth <= CD_ZERO .OR. &
        .NOT. CD_Is_Finite(direction_deg) .OR. .NOT. CD_Is_Finite(gravity) .OR. gravity <= CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Irregular_Waves: need depth>0, gravity>0, finite direction'
      RETURN
    END IF

    ALLOCATE (omc(n), kc(n), phc(n))
    DO i = 1, n
      omc(i) = 2.0_wp*PI_L/period(i)
      CALL CD_Solve_Dispersion_Wavenumber(omc(i), depth, gravity, kc(i), es, em2)
      IF (es /= CD_HYDRO_OK) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Irregular_Waves: '//TRIM(em2); RETURN
      END IF
      phc(i) = phase_deg(i)*PI_L/180.0_wp
    END DO

    ! Snapshot for rollback (transactional; see Set_Drag). Displaces any regular field.
    had_wave = model%has_wave
    a_old = model%a
    h_old = model%wv_h; om_old = model%wv_om; k_old = model%wv_k
    dep_old = model%wv_depth; dir_old = model%wv_dir
    nc_old = model%wv_ncomp
    IF (nc_old > 0) THEN
      ampc_old = model%wv_ampc; omc_old = model%wv_omc
      kc_old = model%wv_kc; phc_old = model%wv_phc
    ELSE
      ALLOCATE (ampc_old(0), omc_old(0), kc_old(0), phc_old(0))
    END IF

    model%wv_ampc = amplitude; model%wv_omc = omc; model%wv_kc = kc; model%wv_phc = phc
    model%wv_ncomp = n
    model%wv_h = CD_ZERO; model%wv_om = CD_ZERO; model%wv_k = CD_ZERO
    model%wv_depth = depth; model%wv_dir = direction_deg
    model%has_wave = .TRUE.
    CALL refresh_wave_table(model)

    ! Refresh the consistent acceleration under the wave loads at the current simulation time.
    CALL solve_consistent_acceleration(model, ErrStat, ErrMsg)
    IF (ErrStat /= CD_HCDYN_OK) THEN
      model%a = a_old
      model%wv_h = h_old; model%wv_om = om_old; model%wv_k = k_old
      model%wv_depth = dep_old; model%wv_dir = dir_old
      model%wv_ncomp = nc_old
      IF (nc_old > 0) THEN
        model%wv_ampc = ampc_old; model%wv_omc = omc_old
        model%wv_kc = kc_old; model%wv_phc = phc_old
      END IF
      model%has_wave = had_wave
      CALL refresh_wave_table(model)
    END IF
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Irregular_Waves

  SUBROUTINE refresh_wave_table(model)
    !! exp(-2 k_i depth) of the component wave table, once per table instead of per evaluation
    !! (the same expression CD_Component_Wave_Kinematics would evaluate).
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER :: i
    IF (model%wv_ncomp < 1) RETURN
    IF (ALLOCATED(model%wv_e2kd)) THEN
      IF (SIZE(model%wv_e2kd) /= model%wv_ncomp) DEALLOCATE (model%wv_e2kd)
    END IF
    IF (.NOT. ALLOCATED(model%wv_e2kd)) ALLOCATE (model%wv_e2kd(model%wv_ncomp))
    DO i = 1, model%wv_ncomp
      model%wv_e2kd(i) = EXP(-2.0_wp*(model%wv_kc(i)*model%wv_depth))
    END DO
  END SUBROUTINE refresh_wave_table

  SUBROUTINE CD_HermiteCable_Dyn_Set_Held_Fluid(model, fluid_velocity, fluid_acceleration, waterline_z, &
                                                ErrStat, ErrMsg)
    !! Prescribe the HELD ambient-fluid field: per-node fluid velocity, fluid
    !! acceleration, and local free-surface elevation, frozen over the coming step(s)
    !! until the next call (the external host's sampling cadence -- a per-step
    !! zero-order hold). The velocity adds to the drag's relative-velocity field (on
    !! top of any constant current from the drag config); the acceleration drives the
    !! Froude-Krylov + fluid-inertia load with the added-mass config's coefficients;
    !! the elevation drives the wetting of both. Consumed per element as the two end
    !! nodes' average -- the same spatial granularity the host sampled at. Requires the
    !! drag configuration; mutually exclusive with the internal wave families (a deck
    !! wave model AND a host field would double-count). The first call establishes a
    !! consistent initial acceleration. Later calls prescribe the load for the next
    !! interval without rewriting the committed generalised-alpha acceleration, which
    !! belongs to the preceding time station.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: fluid_velocity(:, :), fluid_acceleration(:, :), waterline_z(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: istat, es
    LOGICAL :: had_field
    CHARACTER(240) :: em
    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    ! A load-configuration change invalidates the cached committed force F_n. The carried step
    ! factor stays usable (a per-step field update only perturbs the drag and inertia tangent;
    ! the first field's consistent-acceleration solve drops it).
    model%fc_valid = .FALSE.
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: model not initialised'; RETURN
    END IF
    IF (.NOT. model%has_drag) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: requires the Morison drag configuration '// &
               '(CD_HermiteCable_Dyn_Set_Drag) -- the held velocity feeds its fluid field'; RETURN
    END IF
    IF (model%has_wave) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: an internal wave family is configured; a held '// &
               'host field on top of it would double-count the wave kinematics'; RETURN
    END IF
    IF (SIZE(fluid_velocity, 1) /= 3 .OR. SIZE(fluid_velocity, 2) /= model%nn .OR. &
        SIZE(fluid_acceleration, 1) /= 3 .OR. SIZE(fluid_acceleration, 2) /= model%nn .OR. &
        SIZE(waterline_z) /= model%nn) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: fields must be shaped (3, nn) / (nn)'; RETURN
    END IF
    IF (.NOT. (CD_All_Finite(fluid_velocity) .AND. CD_All_Finite(fluid_acceleration) .AND. &
               CD_All_Finite(waterline_z))) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: fields must be finite'; RETURN
    END IF
    ! Transactional update. The FIRST held field completes dynamic initialisation and
    ! therefore refreshes a0. A later field is the load over the NEXT interval; replacing
    ! a_n with the acceleration implied by that future load corrupts the Newmark history
    ! and can pump high-frequency energy into a cable driven by SeaState. Validate later
    ! fields at the committed q_n/v_n, but preserve the committed a_n exactly.
    had_field = model%has_held_field
    IF (had_field .AND. .NOT. (ALLOCATED(model%hf_u) .AND. ALLOCATED(model%hf_ud) .AND. &
                               ALLOCATED(model%hf_wl))) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: committed field storage is incomplete'; RETURN
    END IF
    IF (.NOT. had_field .AND. (ALLOCATED(model%hf_u) .OR. ALLOCATED(model%hf_ud) .OR. &
                               ALLOCATED(model%hf_wl))) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: uncommitted field storage is present'; RETURN
    END IF
    IF (had_field) THEN
      IF (SIZE(model%hf_u, 1) /= 3 .OR. SIZE(model%hf_u, 2) /= model%nn .OR. &
          SIZE(model%hf_ud, 1) /= 3 .OR. SIZE(model%hf_ud, 2) /= model%nn .OR. &
          SIZE(model%hf_wl) /= model%nn) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: committed field storage has invalid dimensions'; RETURN
      END IF
    END IF
    istat = 0
    IF (ALLOCATED(model%ws_hf_a_old)) THEN
      IF (SIZE(model%ws_hf_a_old) /= SIZE(model%a) .OR. SIZE(model%ws_hf_wl_old) /= model%nn) &
        DEALLOCATE (model%ws_hf_a_old, model%ws_hf_u_old, model%ws_hf_ud_old, model%ws_hf_wl_old)
    END IF
    IF (.NOT. ALLOCATED(model%ws_hf_a_old)) &
      ALLOCATE (model%ws_hf_a_old(SIZE(model%a)), model%ws_hf_u_old(3, model%nn), &
                model%ws_hf_ud_old(3, model%nn), model%ws_hf_wl_old(model%nn), STAT=istat)
    IF (istat /= 0) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: rollback storage allocation failed'; RETURN
    END IF
    model%ws_hf_a_old = model%a
    IF (had_field) THEN
      model%ws_hf_u_old = model%hf_u; model%ws_hf_ud_old = model%hf_ud; model%ws_hf_wl_old = model%hf_wl
    END IF
    IF (.NOT. ALLOCATED(model%hf_u)) THEN
      ALLOCATE (model%hf_u(3, model%nn), model%hf_ud(3, model%nn), model%hf_wl(model%nn), STAT=istat)
      IF (istat /= 0) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: field allocation failed'; RETURN
      END IF
    END IF
    model%hf_u = fluid_velocity
    model%hf_ud = fluid_acceleration
    model%hf_wl = waterline_z
    model%has_held_field = .TRUE.
    IF (.NOT. had_field .AND. .NOT. model%a_restored) THEN
      CALL solve_consistent_acceleration(model, ErrStat, ErrMsg)
    ELSE
      CALL assemble_residual(model, model%q, model%v, model%t, .TRUE., .FALSE., es, em)
      ErrStat = es
      ErrMsg = em
      IF (es == CD_HCDYN_OK .AND. model%force_blend) THEN
        ! This validation is the next step's committed force F_n = R(q_n, v_n, t_n) under the
        ! new field, the evaluation the step would otherwise repeat: cache it, keyed as usual.
        model%fc_R = model%ws_R
        model%fc_q = model%q
        model%fc_v = model%v
        model%fc_t = model%t
        model%fc_d0 = model%endconn_d0
        model%fc_valid = .TRUE.
      END IF
    END IF
    IF (ErrStat == CD_HCDYN_OK) model%a_restored = .FALSE.
    IF (ErrStat /= CD_HCDYN_OK) THEN
      model%a = model%ws_hf_a_old
      IF (had_field) THEN
        model%hf_u = model%ws_hf_u_old; model%hf_ud = model%ws_hf_ud_old; model%hf_wl = model%ws_hf_wl_old
        model%has_held_field = .TRUE.
      ELSE
        IF (ALLOCATED(model%hf_u)) DEALLOCATE (model%hf_u, model%hf_ud, model%hf_wl)
        model%has_held_field = .FALSE.
      END IF
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Held_Fluid: '//TRIM(ErrMsg)
    END IF
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Held_Fluid

  SUBROUTINE CD_HermiteCable_Dyn_Step(model, dt, max_iter, tol, ErrStat, ErrMsg, iters_out, res_out, &
                                      pres_dofs, pres_q, pres_v, pres_a, residual_history, history_count, &
                                      endconn_direction, endconn_direction_rate, endconn_direction_acceleration, &
                                      tensile_failed)
    !! Advance the state one generalised-alpha step of size dt. Newton-iterates the intermediate
    !! residual to a dimensionless relative tolerance, then commits q_{n+1}, v_{n+1}, a_{n+1}.
    !!
    !! Optional PRESCRIBED MOTION: pres_dofs(:) lists global DOFs (which MUST be in the model's
    !! fixed_dofs set) whose kinematics at t_{n+1} are prescribed to pres_q/pres_v/pres_a(:) --
    !! e.g. a moving fairlead/hang-off heave driving a lazy-wave cable. These DOFs are held at
    !! pres_q during the Newton solve and their prescribed q/v/a enter the free-DOF residual
    !! through the mass and stiffness coupling; all four arrays must be supplied together.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT), TARGET :: model
    REAL(wp), INTENT(IN) :: dt, tol
    INTEGER, INTENT(IN) :: max_iter
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(OUT), OPTIONAL :: iters_out
    REAL(wp), INTENT(OUT), OPTIONAL :: res_out
    INTEGER, INTENT(IN), OPTIONAL :: pres_dofs(:)
    REAL(wp), INTENT(IN), OPTIONAL :: pres_q(:), pres_v(:), pres_a(:)
    REAL(wp), INTENT(OUT), OPTIONAL :: residual_history(:)
    INTEGER, INTENT(OUT), OPTIONAL :: history_count
    REAL(wp), INTENT(IN), OPTIONAL :: endconn_direction(3, 2)
    REAL(wp), INTENT(IN), OPTIONAL :: endconn_direction_rate(3, 2), endconn_direction_acceleration(3, 2)
    !! .TRUE. when the step failed the tensile qualification of the converged field (ERROR
    !! mode), as opposed to a Newton failure; the state is not committed in either case.
    LOGICAL, INTENT(OUT), OPTIONAL :: tensile_failed

    INTEGER :: ndof, it, es, i, np, ip, ls, gdof, nsolve_step, axial_element, axial_es
    INTEGER :: axial_worst_element, iv, attempt, n_attempt
    REAL(wp) :: c_a, c_v, c_a2, am, af, bt, gm, rnorm, fscale, rnt, lambda, chosen, lam_ok, cw, t_alpha
    REAL(wp) :: e_bal, e_ma
    REAL(wp) :: axial_qe(12), axial_min
    REAL(wp) :: axial_threshold, axial_margin, axial_worst_margin, axial_worst_force, axial_worst_xi
    REAL(wp) :: axial_worst_threshold
    ! Axial audit window (elements e-1, e, e+1; see the audit after convergence).
    INTEGER :: win_elem(3)
    LOGICAL :: win_clear(3), win_done(3)
    REAL(wp) :: win_mean(3), win_len, win_sum, win_ea
    LOGICAL :: has_pres, scale_set, have_tangent, mn, mn_factored, mn_fresh, has_endconn_motion
    LOGICAL :: fc_hit, xs_active, factor_ok, use_ca, xs_start
    INTEGER :: n_fact
    INTEGER :: n_stale
    REAL(wp) :: step_endconn_target(3, 2), step_endconn_eval(3, 2)
    REAL(wp) :: step_endconn_rate(3, 2), step_endconn_acceleration(3, 2)
    REAL(wp) :: step_rigid_basis(3, 3, 2)
    REAL(wp) :: step_n0(2), step_ndot0(2), step_nddot0(2)
    ! Modified-Newton contraction threshold: an accepted iterate must cut the residual
    ! to <= MN_CONTRACT * rnorm for the reused factor to be kept. ONE constant feeds both
    ! the accept-time validate-and-build decision and the commit-time rebuild schedule --
    ! they must agree, or a weak accept could skip validation yet still get rebuilt at.
    REAL(wp), PARAMETER :: MN_CONTRACT = 0.25_wp
    ! Constant-acceleration predictor gate: relative force imbalance of the committed state.
    REAL(wp), PARAMETER :: PRED_BALANCE = 0.5_wp
    ! Cross-step reuse: at most this many accepted iterations on the carried factor per step;
    ! a solve still short of the tolerance then refreshes (quadratic convergence from there).
    INTEGER, PARAMETER :: XS_MAX_STALE = 3
    ! A carried-factor iterate contracting the residual by less than this abandons the attempt
    ! (see the abandon test in the Newton loop); between MN_CONTRACT and it, the factor is
    ! refreshed at the iterate as usual.
    REAL(wp), PARAMETER :: XS_ABANDON = 0.5_wp
    ! Longest run of fresh-factor steps after repeated refreshes (see xs_backoff).
    INTEGER, PARAMETER :: XS_MAX_BACKOFF = 32
    ! A carried factor converges linearly, so its iterate stops just inside the tolerance where a
    ! fresh tangent's quadratic last step lands far inside it: iterates on it must reach
    ! XS_TOL_FACTOR*tol, which keeps the outputs as close to the converged solution as full
    ! Newton's.
    REAL(wp), PARAMETER :: XS_TOL_FACTOR = 0.3_wp
    INTEGER(INT64) :: pc0, pc1, pcr, pt0, pt1
    CHARACTER(200) :: em2

    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    nsolve_step = 0
    IF (model%torsion%active) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step: torsion is supported in statics only in this build'
      RETURN
    END IF
    ! Clock scratch is defined even when profiling is off, purely so -Wmaybe-uninitialized
    ! can prove it (every read sits under the same prof_enabled guard as its write).
    pt0 = 0_INT64; pc0 = 0_INT64; pcr = 1_INT64
    IF (prof_enabled) THEN
      CALL SYSTEM_CLOCK(pt0, pcr)
      prof_n_step = prof_n_step + 1
    END IF
    IF (PRESENT(iters_out)) iters_out = 0
    IF (PRESENT(res_out)) res_out = CD_ZERO
    IF (PRESENT(history_count)) history_count = 0
    IF (PRESENT(residual_history)) residual_history = CD_ZERO
    IF (PRESENT(tensile_failed)) tensile_failed = .FALSE.
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Step: model not initialised'; RETURN
    END IF
    IF (.NOT. CD_Is_Finite(dt) .OR. dt <= CD_ZERO .OR. max_iter < 1 .OR. &
        .NOT. CD_Is_Finite(tol) .OR. tol <= CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step: need dt>0, max_iter>=1, finite tol>0'; RETURN
    END IF
    IF (PRESENT(residual_history)) THEN
      IF (SIZE(residual_history) < max_iter) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Step: residual_history must have at least max_iter entries'; RETURN
      END IF
    END IF

    ndof = model%ndof

    has_endconn_motion = PRESENT(endconn_direction)
    step_endconn_target = model%endconn_d0
    step_endconn_eval = model%endconn_d0
    step_endconn_rate = CD_ZERO
    step_endconn_acceleration = CD_ZERO
    IF (PRESENT(endconn_direction_rate) .NEQV. PRESENT(endconn_direction_acceleration)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step: direction rate and acceleration must be supplied together'
      RETURN
    END IF
    IF (PRESENT(endconn_direction_rate) .AND. .NOT. has_endconn_motion) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step: direction derivatives require endconn_direction'
      RETURN
    END IF
    IF (has_endconn_motion) THEN
      CALL prepare_endconn_directions(model, endconn_direction, step_endconn_target, es, em2)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es
        ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2)
        RETURN
      END IF
      CALL prepare_endconn_derivatives(model, step_endconn_target, step_endconn_rate, &
                                       step_endconn_acceleration, es, em2, endconn_direction_rate, &
                                       endconn_direction_acceleration)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es
        ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2)
        RETURN
      END IF
      CALL interpolate_endconn_directions(model, model%endconn_d0, step_endconn_target, &
                                          MERGE(CD_ONE, CD_ONE - model%alpha_f, model%force_blend), &
                                          step_endconn_eval, es, em2)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es
        ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2)
        RETURN
      END IF
    END IF
    CALL build_rigid_bases(model, step_endconn_target, step_rigid_basis, es, em2)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = es
      ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2)
      RETURN
    END IF

    ! Prescribed-motion validation: all four arrays together, sizes match, DOFs are fixed + finite.
    ! Reject ANY partial set -- including supplying pres_q/v/a while omitting pres_dofs, which would
    ! otherwise silently skip the drive and advance the boundary as held instead of failing the call.
    has_pres = PRESENT(pres_dofs) .AND. PRESENT(pres_q) .AND. PRESENT(pres_v) .AND. PRESENT(pres_a)
    np = 0
    IF ((PRESENT(pres_dofs) .OR. PRESENT(pres_q) .OR. PRESENT(pres_v) .OR. PRESENT(pres_a)) &
        .AND. .NOT. has_pres) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step: prescribed motion needs pres_dofs/q/v/a together'; RETURN
    END IF
    IF (has_pres) THEN
      np = SIZE(pres_dofs)
      IF (np < 1) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Step: pres_dofs must not be empty'; RETURN
      END IF
      IF (SIZE(pres_q) /= np .OR. SIZE(pres_v) /= np .OR. SIZE(pres_a) /= np) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Step: pres_q/v/a must match pres_dofs length'; RETURN
      END IF
      DO ip = 1, np
        gdof = pres_dofs(ip)
        IF (gdof < 1 .OR. gdof > ndof) THEN
          ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Step: pres_dof out of range'; RETURN
        END IF
        IF (model%freemask(gdof)) THEN
          ErrStat = CD_HCDYN_BADINPUT
          ErrMsg = 'CD_HermiteCable_Dyn_Step: prescribed DOF must be in fixed_dofs'; RETURN
        END IF
        IF (.NOT. CD_Is_Finite(pres_q(ip)) .OR. .NOT. CD_Is_Finite(pres_v(ip)) .OR. &
            .NOT. CD_Is_Finite(pres_a(ip))) THEN
          ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Step: prescribed q/v/a must be finite'; RETURN
        END IF
      END DO
    END IF
    am = model%alpha_m; af = model%alpha_f; bt = model%beta; gm = model%gamma
    ! Newmark acceleration-from-displacement coefficients.
    c_a = CD_ONE/(bt*dt*dt)              ! d a_{n+1} / d q_{n+1}
    c_v = CD_ONE/(bt*dt)
    c_a2 = CD_ONE/(2.0_wp*bt) - CD_ONE
    ! Velocity-tangent coupling: d(v_alpha)/d(q_{n+1}) = (1-af) * d(v_{n+1})/d(q_{n+1}) =
    ! (1-af) * gamma/(beta dt). Folds a velocity-dependent load (drag) into the effective tangent.
    cw = (CD_ONE - af)*gm/(bt*dt)
    ! FORCE-BLENDED generalised alpha (force_blend, the default): G = M a_alpha
    ! + (1-af) F(q_{n+1}, v_{n+1}, t_{n+1}) + af F(q_n, v_n, t_n). Blending the
    ! configurations instead evaluates the internal force at (1-af) q_{n+1} + af q_n, whose
    ! rotating Hermite tangents are shortened by the blend: a spurious stretch of order
    ! af (1-af) (omega dt)^2 EA that shifts the mean tension with dt. The force blend has
    ! the same second-order accuracy and tangent (1-af) K(q_{n+1}). Time-dependent (wave)
    ! loads evaluate at t_{n+1}, or at t_alpha for the configuration blend.
    IF (model%force_blend) THEN
      t_alpha = model%t + dt
    ELSE
      t_alpha = model%t + (CD_ONE - af)*dt
    END IF

    ! Workspace lives in the model (allocated once at Init): a dynamic step performs ZERO
    ! heap allocations. Local names below are the model's ws_* arrays.
    DO iv = 1, SIZE(model%ws_qn)
      model%ws_qn(iv) = model%q(iv)
      model%ws_vn(iv) = model%v(iv)
      model%ws_an(iv) = model%a(iv)
    END DO
    CALL rigid_tangent_scalar_state(model, model%endconn_d0, model%ws_qn, model%ws_vn, model%ws_an, &
                                    step_n0, step_ndot0, step_nddot0)
    ! Total mass for this step: structural + added mass FROZEN at the beginning-of-step
    ! configuration q_n (constant within the Newton -> no dM/dq term; O(dt) freeze error consistent
    ! with gen-alpha accuracy, exact at any at-rest fixed point). All band arrays share the
    ! DGBSV layout, so the sum is a plain array operation.
    CALL band_copy(model%Mgb, model%ws_Msumb)
    IF (model%has_am) THEN
      IF (prof_enabled) CALL SYSTEM_CLOCK(pc0)
      CALL assemble_added_mass(model, model%ws_qn, es, em2)
      IF (prof_enabled) THEN
        CALL SYSTEM_CLOCK(pc1)
        prof_n_am = prof_n_am + 1
        prof_t_am = prof_t_am + REAL(pc1 - pc0, wp)/REAL(pcr, wp)
      END IF
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2); RETURN
      END IF
      CALL band_add(model%ws_Mamb, model%ws_Msumb)
    END IF
    ! A fixed DOF that is NOT prescribed this step is a boundary HELD AT REST: zero its velocity and
    ! acceleration so the predictor keeps it at q_n and the commit leaves v = a = 0. Without this a
    ! DOF driven by pres_* in an earlier step and then released (omitted from pres_dofs) would keep
    ! drifting at its last prescribed rate, since Newton never touches fixed DOFs. Prescribed DOFs
    ! restore their own v_n/a_n history (needed for the gen-alpha blend and the pres_* override).
    DO i = 1, ndof
      IF (.NOT. model%freemask(i)) THEN
        model%ws_vn(i) = CD_ZERO; model%ws_an(i) = CD_ZERO
      END IF
    END DO
    IF (has_pres) THEN
      DO ip = 1, np
        model%ws_vn(pres_dofs(ip)) = model%v(pres_dofs(ip))
        model%ws_an(pres_dofs(ip)) = model%a(pres_dofs(ip))
      END DO
    END IF
    ! Committed force F_n = R(q_n, v_n, t_n) with the current load configuration: reused
    ! from the previous step's converged evaluation when its key matches the committed
    ! state exactly, evaluated here otherwise (one residual assembly).
    fc_hit = .FALSE.
    IF (model%force_blend) fc_hit = force_cache_matches()
    ! Cross-step factor reuse needs the step to continue bitwise from the state the carried
    ! factor's step committed (the force-cache key), at the same dt and without a moving rigid
    ! end frame (whose basis transform is baked into the factor). The carried factor is
    ! consumed here: only this step's commit re-validates ws_Keffb.
    xs_active = model%tangent_reuse .AND. model%tr_valid .AND. fc_hit .AND. .NOT. has_endconn_motion
    IF (xs_active) xs_active = .NOT. (ABS(dt - model%tr_dt) > CD_ZERO)
    factor_ok = xs_active
    model%tr_valid = .FALSE.
    IF (model%force_blend .AND. .NOT. fc_hit) THEN
      CALL assemble_residual(model, model%q, model%v, model%t, .TRUE., .FALSE., es, em2)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2); RETURN
      END IF
      DO iv = 1, SIZE(model%fc_R)
        model%fc_R(iv) = model%ws_R(iv)
      END DO
      ! Key the recomputed force on the state it describes.
      model%fc_q = model%q
      model%fc_v = model%v
      model%fc_t = model%t
      model%fc_d0 = model%endconn_d0
      model%fc_valid = .TRUE.
    END IF
    ! Newton predictor. The standard Newmark guess q_n + dt v_n + dt^2 (1/2 - beta) a_n (the
    ! acceleration map's base a_{n+1} = 0) is a far better start than q_n when the state is
    ! moving. When the committed acceleration is the physical one, the constant-acceleration
    ! guess a_{n+1} = a_n, q_n + dt v_n + dt^2/2 a_n, starts an order of magnitude closer: it
    ! saves an iteration per step and, at a loose tolerance, removes the a = 0 bias of the
    ! Newmark guess. It is NOT safe while an under-resolved (omega dt > 1) mode rings: its
    ! algorithmic acceleration is not a smooth signal, a loosely converged step keeps the
    ! extrapolated value and the mode is no longer damped, whereas the a = 0 guess damps it.
    ! The test uses only the committed state (restart-exact, no history): the extrapolation
    ! runs when the committed acceleration balances the committed force, |M a_n + F_n| <=
    ! PRED_BALANCE |M a_n| on the translational free DOFs (F_n = fc_R, force blend only). A
    ! failed constant-acceleration attempt re-solves the step from the Newmark guess, so
    ! robustness is never below the Newmark-predictor solve.
    n_attempt = 1
    use_ca = .FALSE.
    e_bal = CD_ZERO
    e_ma = CD_ZERO
    IF (model%force_blend) THEN
      CALL band_matvec(model%ws_Msumb, model%ws_an, model%ws_inert)
      DO iv = 1, ndof
        IF (.NOT. model%solve_mask(iv) .OR. MOD(iv - 1, 6) >= 3) CYCLE
        e_bal = e_bal + (model%ws_inert(iv) + model%fc_R(iv))**2
        e_ma = e_ma + model%ws_inert(iv)**2
      END DO
      use_ca = e_ma > CD_ZERO .AND. e_bal <= PRED_BALANCE**2*e_ma
    END IF
    ! The carried factor is used under the same smoothness test: through a transient it
    ! stops each step just inside the tolerance, and those states can lead the start-up of a
    ! stiff cable into a step no Newton converges.
    xs_active = xs_active .AND. use_ca
    IF (xs_active .AND. model%xs_skip > 0) THEN
      model%xs_skip = model%xs_skip - 1
      xs_active = .FALSE.
    END IF
    xs_start = xs_active
    n_fact = 0
    ! Each acceleration is dropped in turn on failure: a failed carried-factor solve retries
    ! the constant-acceleration guess on a fresh tangent (a stale factor can lead the line
    ! search off the smooth solution, while the guess itself is still the best start), and
    ! a failed constant-acceleration solve retries the Newmark guess on a fresh tangent,
    ! i.e. the plain solve. Dropping both at once restarts a step the constant-acceleration
    ! guess converges in one iteration from the a = 0 guess, which at a large dt can reach
    ! a different (spurious) root.
    IF (use_ca) n_attempt = 2
    IF (xs_active) n_attempt = 3
    predict: DO attempt = 1, n_attempt
    IF (attempt > 1) THEN
      ! The retry is a fresh solve: its iteration count is the step's.
      ErrStat = CD_HCDYN_OK
      ErrMsg = ''
      nsolve_step = 0
      IF (PRESENT(iters_out)) iters_out = 0
      IF (PRESENT(res_out)) res_out = CD_ZERO
      IF (PRESENT(history_count)) history_count = 0
      IF (PRESENT(residual_history)) residual_history = CD_ZERO
      IF (.NOT. (xs_start .AND. attempt == 2)) use_ca = .FALSE.
      xs_active = .FALSE.
    END IF
    n_stale = 0
    IF (use_ca) THEN
      DO iv = 1, SIZE(model%ws_qk)
        model%ws_qk(iv) = model%ws_qn(iv) + dt*model%ws_vn(iv) + 0.5_wp*dt*dt*model%ws_an(iv)
      END DO
    ELSE
      DO iv = 1, SIZE(model%ws_qk)
        model%ws_qk(iv) = model%ws_qn(iv) + dt*model%ws_vn(iv) + (0.5_wp - bt)*dt*dt*model%ws_an(iv)
      END DO
    END IF
    ! Hold the prescribed DOFs at their target position for the whole solve (Newton never touches
    ! them -- they are in fixed_dofs, hence not in fmap).
    IF (has_pres) THEN
      DO ip = 1, np
        model%ws_qk(pres_dofs(ip)) = pres_q(ip)
      END DO
    END IF
    CALL project_rigid_tangents(model, step_endconn_target, model%ws_qk, es, em2)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = es
      ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2)
      IF (attempt < n_attempt) CYCLE predict
      RETURN
    END IF

    ! Residual + tangent at the predictor, then damped (line-searched) Newton. The line search on
    ! the residual norm is what keeps the solve robust at the large implicit dt this element is
    ! meant to exploit (axial stiffness EA >> EI makes an undamped full step overshoot on a stiff
    ! transient) -- mirroring the damped Newton the static solver uses for the same geometry.
    ! The force scale fscale is FROZEN at the predictor so every rnorm in the step (predictor and
    ! all line-search trials) is divided by the same denominator -- otherwise a per-configuration
    ! fscale makes the line-search comparison rnt < rnorm inconsistent and rejects good steps.
    scale_set = .FALSE.
    have_tangent = .FALSE.
    mn = model%modified_newton
    ! Adaptive policy overrides the flag: full Newton during the warmup, then the decided mode.
    IF (model%newton_adaptive) THEN
      IF (model%na_steps_done < model%na_warmup) THEN
        mn = .FALSE.
      ELSE
        mn = model%na_mn_active
      END IF
    END IF
    ! A carried factor starts the iteration already factored (and stale).
    mn_factored = xs_active
    mn_fresh = .FALSE.
    CALL eval_res(model%ws_qk, .FALSE., rnorm, es, em2)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2)
      IF (attempt < n_attempt) CYCLE predict
      RETURN
    END IF

    DO it = 1, max_iter
      IF (PRESENT(residual_history)) residual_history(it) = rnorm
      IF (PRESENT(history_count)) history_count = it
      IF (PRESENT(iters_out)) iters_out = nsolve_step
      IF (PRESENT(res_out)) res_out = rnorm
      IF (rnorm < MERGE(XS_TOL_FACTOR*tol, tol, xs_active)) EXIT
      IF (it == max_iter) THEN
        ErrStat = CD_HCDYN_NOCONVERGE
        BLOCK
          ! Format into a local buffer: a record longer than the caller's ErrMsg
          ! must truncate, not abort the internal write.
          CHARACTER(1024) :: wbuf
          INTEGER :: wios
          wbuf = ''
          WRITE (wbuf, '(A,ES12.5,A,ES12.5)', IOSTAT=wios) &
            'CD_HermiteCable_Dyn_Step: Newton did not converge, ||G||=', rnorm, ' > tol=', tol
          ErrMsg = wbuf
        END BLOCK
        IF (attempt < n_attempt) CYCLE predict
        RETURN
      END IF
      ! Tangent assembly ONCE per solve iteration, at the alpha-blend states (ws_qa/ws_va)
      ! the accepted (or predictor) residual evaluation just stored -- the evaluation point
      ! the Newton iteration linearizes about, so rejected line-search trials and each
      ! step's CONVERGED accepted evaluation do not pay for tangents that are never
      ! factored. An accepted trial that will iterate again
      ! validates and builds its tangent inside the line search (have_tangent), so this
      ! fires only for the predictor iterate.
      !
      ! MODIFIED NEWTON (model%modified_newton): the tangent is assembled + factored once
      ! per STEP (at the predictor iterate) and the FACTORIZATION is reused for every
      ! subsequent iteration -- only the residual back-substitution repeats. The factored
      ! band + pivots persist untouched in ws_Keffb/ws_ipiv between iterations. mn_fresh
      ! records whether the current factor was built at the current iterate (the refresh
      ! path resets it), so a failed line search can tell a stale direction from a
      ! genuinely hard step. The stale-direction refresh retries WITHIN this named
      ! construct -- inside the SAME outer Newton iteration -- so it is never charged
      ! against max_iter (a refresh on the last allowed iteration must still be tried,
      ! not aborted by the budget check before its fresh direction ever runs).
      !
      ! CROSS-STEP REUSE (xs_active): the same machinery runs on the factor carried from the
      ! previous step, with the same MN_CONTRACT gate. Its first refresh ends the reuse for the
      ! rest of the step (full Newton, unless modified_newton keeps reusing within the step).
      mn_retry: DO
      IF (.NOT. ((mn .OR. xs_active) .AND. mn_factored)) THEN
        IF (.NOT. have_tangent) THEN
          IF (prof_enabled) CALL SYSTEM_CLOCK(pc0)
          CALL assemble_residual(model, model%ws_qa, model%ws_va, t_alpha, .FALSE., .TRUE., es, em2, &
                                 endconn_direction_eval=step_endconn_eval)
          IF (prof_enabled) THEN
            CALL SYSTEM_CLOCK(pc1)
            prof_n_tan = prof_n_tan + 1
            prof_t_tan = prof_t_tan + REAL(pc1 - pc0, wp)/REAL(pcr, wp)
          END IF
          IF (es /= CD_HCDYN_OK) THEN
            ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2)
            IF (attempt < n_attempt) CYCLE predict
            RETURN
          END IF
        END IF
        ! Effective tangent dG/dq_{n+1} = (1-am) c_a M + (1-af) K + cw Kv -- a plain array
        ! operation on the row-aligned band arrays. Dirichlet applies IN-BAND (row/column
        ! zeroed within the band, unit diagonal, zero residual slot), so the full-size banded
        ! solve returns exactly zero update at every fixed DOF -- no free-DOF reduction copy.
        factor_ok = .FALSE.
        IF (prof_enabled) CALL SYSTEM_CLOCK(pc0)
        CALL band_combine((CD_ONE - am)*c_a, model%ws_Msumb, CD_ONE - af, model%ws_Kb, cw, model%ws_Kvb, &
                          model%ws_Keffb)
        IF (has_endconn_motion) CALL add_rigid_transport_tangent()
        CALL transform_rigid_band(model, step_rigid_basis, model%ws_Keffb)
        CALL dirichlet_band(model%ws_Keffb, model%solve_mask)
        IF (prof_enabled) THEN
          CALL SYSTEM_CLOCK(pc1)
          prof_t_eff = prof_t_eff + REAL(pc1 - pc0, wp)/REAL(pcr, wp)
        END IF
        ! Factor with the model-owned pivot workspace (DGBTRF -- exactly what DGBSV runs
        ! internally, without its per-call pivot allocation).
        IF (prof_enabled) CALL SYSTEM_CLOCK(pc0)
        CALL CD_Factor_Banded(model%ws_Keffb, KL_D, KU_D, model%ws_ipiv, es, em2, matrix_validated=.TRUE.)
        IF (prof_enabled) THEN
          CALL SYSTEM_CLOCK(pc1)
          prof_t_solve = prof_t_solve + REAL(pc1 - pc0, wp)/REAL(pcr, wp)
        END IF
        IF (es /= CD_LINALG_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Dyn_Step: tangent factor: '//TRIM(em2)
          IF (attempt < n_attempt) CYCLE predict
          RETURN
        END IF
        mn_factored = .TRUE.
        mn_fresh = .TRUE.
        factor_ok = .TRUE.
        xs_active = .FALSE.
        n_fact = n_fact + 1
        ! The line-search build (if any) is now consumed by this factor. Clearing the
        ! flag matters under modified Newton: a weak-contraction accept builds its
        ! tangent, the scheduled rebuild factors it here, and a LATER refresh must not
        ! mistake that already-consumed assembly for one taken at the refresh iterate.
        have_tangent = .FALSE.
      END IF
      DO iv = 1, SIZE(model%ws_dq)
        model%ws_dq(iv) = -model%ws_G(iv)
      END DO
      DO i = 1, ndof
        IF (.NOT. model%solve_mask(i)) model%ws_dq(i) = CD_ZERO
      END DO
      IF (prof_enabled) CALL SYSTEM_CLOCK(pc0)
      CALL CD_Solve_Factored_Banded(model%ws_Keffb, KL_D, KU_D, model%ws_ipiv, model%ws_dq, es, em2, &
                                    factor_validated=.TRUE.)
      IF (prof_enabled) THEN
        CALL SYSTEM_CLOCK(pc1)
        prof_n_solve = prof_n_solve + 1
        prof_t_solve = prof_t_solve + REAL(pc1 - pc0, wp)/REAL(pcr, wp)
      END IF
      IF (es /= CD_LINALG_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Dyn_Step: tangent solve: '//TRIM(em2)
        IF (attempt < n_attempt) CYCLE predict
        RETURN
      END IF
      CALL transform_rigid_vector_to_global(model, step_rigid_basis, model%ws_dq)
      nsolve_step = nsolve_step + 1
      IF (PRESENT(iters_out)) iters_out = nsolve_step
      ! Backtracking line search: accept the first step length that reduces the residual (and gives
      ! a finite element geometry); if none does, take a small safe damped step so the iteration
      ! never commits an overshoot. max_iter then bounds a genuinely stuck solve honestly.
      ! dq(fixed) = 0 by construction, so the whole-vector trial leaves fixed DOFs untouched.
      chosen = -CD_ONE
      lam_ok = -CD_ONE                                    ! smallest step length with a finite eval
      lambda = CD_ONE
      DO ls = 1, 10
        DO iv = 1, SIZE(model%ws_qtrial)
          model%ws_qtrial(iv) = model%ws_qk(iv) + lambda*model%ws_dq(iv)
        END DO
        CALL eval_res(model%ws_qtrial, .TRUE., rnt, es, em2)
        IF (es == CD_HCDYN_OK) THEN
          lam_ok = lambda                                 ! loop shrinks lambda, so this ends smallest
          IF (rnt < rnorm) THEN
            ! Accepting a trial that will iterate again requires a FINITE tangent at it: a
            ! residual can stay finite where the Hessian/Jacobians overflow (near-degenerate
            ! geometry, extreme stiffness), and accepting such a point would only abort at
            ! the next factor. Validate-and-build here (the trial's alpha states are current)
            ! and treat an assembly failure exactly like a non-finite trial: keep halving.
            ! A CONVERGED accept skips the tangent -- it is never factored.
            ! Under MODIFIED NEWTON only a STRONG-contraction accept skips the
            ! validate-and-build (the existing factor is kept; no rebuild is scheduled
            ! at this point) -- that skip is the reuse saving. A WEAK-contraction accept
            ! (rnt > 0.25 rnorm, the same predicate the commit's contraction gate uses)
            ! schedules an immediate rebuild AT this point, so it validates-and-builds
            ! exactly like full Newton: committing a finite-residual/overflowing-tangent
            ! geometry and then rebuilding there would abort where full Newton halves
            ! to a usable trial.
            IF (rnt < MERGE(XS_TOL_FACTOR*tol, tol, xs_active)) THEN
              chosen = lambda; EXIT
            END IF
            IF ((mn .OR. xs_active) .AND. rnt <= MN_CONTRACT*rnorm) THEN
              chosen = lambda; EXIT
            END IF
            IF (prof_enabled) CALL SYSTEM_CLOCK(pc0)
            CALL assemble_residual(model, model%ws_qa, model%ws_va, t_alpha, .FALSE., .TRUE., es, em2, &
                                   endconn_direction_eval=step_endconn_eval)
            IF (prof_enabled) THEN
              CALL SYSTEM_CLOCK(pc1)
              prof_n_tan = prof_n_tan + 1
              prof_t_tan = prof_t_tan + REAL(pc1 - pc0, wp)/REAL(pcr, wp)
            END IF
            IF (es == CD_HCDYN_OK) THEN
              have_tangent = .TRUE.
              chosen = lambda; EXIT
            END IF
          END IF
        END IF
        lambda = 0.5_wp*lambda
      END DO
      IF (xs_active .AND. (chosen < CD_ZERO .OR. &
                           (rnt > XS_ABANDON*rnorm .AND. rnt >= XS_TOL_FACTOR*tol))) THEN
        ! The carried factor is kept only while it contracts like Newton. A weak or failed
        ! contraction abandons the attempt instead of refreshing at the current iterate: the
        ! iterates of a stale direction reduce the residual norm without staying near the
        ! solution, and a fresh Newton from there can crawl to a different root. The retry
        ! restarts from the constant-acceleration guess on a fresh tangent.
        ErrStat = CD_HCDYN_NOCONVERGE
        ErrMsg = 'CD_HermiteCable_Dyn_Step: carried factor stopped contracting'
        CYCLE predict
      END IF
      IF (chosen < CD_ZERO .AND. (mn .OR. xs_active) .AND. .NOT. mn_fresh) THEN
        ! The stale direction found no improving step. Before the safe damped fallback,
        ! refresh ONCE: rebuild + refactor at the CURRENT iterate and re-derive the
        ! direction, WITHIN this outer iteration (CYCLE mn_retry, not the Newton loop --
        ! the refresh is bookkeeping, not a step, and must not consume iteration budget).
        ! The rejected trials overwrote ws_qa/ws_va (eval_res stages the alpha states on
        ! every call), so re-evaluate at the current iterate first -- it restores
        ! ws_qa/ws_va/ws_G/rnorm to the iterate's own values, and the retry-top then
        ! assembles the fresh tangent at the right point. At most one refresh fires per
        ! iteration (the rebuild sets mn_fresh, which blocks this guard). If the fresh
        ! direction still cannot improve, the step is genuinely hard and the fallback
        ! below handles it exactly as full Newton would.
        CALL eval_res(model%ws_qk, .FALSE., rnorm, es, em2)
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2)
          IF (attempt < n_attempt) CYCLE predict
          RETURN
        END IF
        mn_factored = .FALSE.
        CYCLE mn_retry
      END IF
      EXIT mn_retry
      END DO mn_retry
      IF (chosen < CD_ZERO) THEN
        ! No step length reduced the residual (with a usable tangent): fall back to the
        ! SMALLEST step that at least gave a finite evaluation (a safe damped step). Using
        ! the actual smallest finite trial -- not a fixed 2^-8 -- means a step whose larger
        ! lengths are non-finite but whose smaller ones are valid is still taken, rather
        ! than spuriously reported as a failure. The fallback validates its tangent the
        ! same way, halving further if needed.
        IF (lam_ok < CD_ZERO) THEN
          ErrStat = CD_HCDYN_NOCONVERGE
          ErrMsg = 'CD_HermiteCable_Dyn_Step: line search: element non-finite at all step lengths'
          IF (attempt < n_attempt) CYCLE predict
          RETURN
        END IF
        lambda = lam_ok
        DO ls = 1, 10
          DO iv = 1, SIZE(model%ws_qtrial)
            model%ws_qtrial(iv) = model%ws_qk(iv) + lambda*model%ws_dq(iv)
          END DO
          CALL eval_res(model%ws_qtrial, .TRUE., rnt, es, em2)
          IF (es == CD_HCDYN_OK) THEN
            ! The fallback step rarely improves the residual, so under modified Newton the
            ! contraction gate will schedule a rebuild at it -- validate-and-build here
            ! exactly like full Newton (the mn saving never came from this rare path).
            IF (rnt < tol) EXIT                           ! converged accept: tangent never factored
            IF (prof_enabled) CALL SYSTEM_CLOCK(pc0)
            CALL assemble_residual(model, model%ws_qa, model%ws_va, t_alpha, .FALSE., .TRUE., es, em2, &
                                   endconn_direction_eval=step_endconn_eval)
            IF (prof_enabled) THEN
              CALL SYSTEM_CLOCK(pc1)
              prof_n_tan = prof_n_tan + 1
              prof_t_tan = prof_t_tan + REAL(pc1 - pc0, wp)/REAL(pcr, wp)
            END IF
            IF (es == CD_HCDYN_OK) THEN
              have_tangent = .TRUE.
              EXIT
            END IF
          END IF
          IF (ls == 10) THEN
            ErrStat = CD_HCDYN_NOCONVERGE
            ErrMsg = 'CD_HermiteCable_Dyn_Step: line search: no step length with a finite '// &
                     'residual and tangent'
            IF (attempt < n_attempt) CYCLE predict
            RETURN
          END IF
          lambda = 0.5_wp*lambda
        END DO
      END IF
      ! commit the accepted iterate + its residual/tangents; the step's reused factor is
      ! now at least one iterate old (stale), so a later stuck line search may refresh it
      IF (mn .OR. xs_active) THEN
        ! CONTRACTION-GATED reuse: keep the factor only while the accepted iterate still
        ! contracts the residual strongly (Newton-like). A stale tangent on a stiff cable
        ! otherwise degrades to slow linear convergence (measured: the pure-reuse variant
        ! plateaus above tol and burns the whole iteration budget), so a weak contraction
        ! schedules a rebuild at the next iteration -- the reuse then pays only in the
        ! well-conditioned phase where it is free.
        IF (rnt > MN_CONTRACT*rnorm) mn_factored = .FALSE.
        IF (xs_active) THEN
          n_stale = n_stale + 1
          IF (n_stale >= XS_MAX_STALE) mn_factored = .FALSE.
        END IF
        mn_fresh = .FALSE.
      END IF
      DO iv = 1, SIZE(model%ws_qk)
        model%ws_qk(iv) = model%ws_qtrial(iv)
      END DO
      DO iv = 1, SIZE(model%ws_G)
        model%ws_G(iv) = model%ws_Gt(iv)
      END DO
      DO iv = 1, SIZE(model%ws_Rk)
        model%ws_Rk(iv) = model%ws_Rt(iv)
      END DO
      rnorm = rnt
    END DO
    EXIT predict
    END DO predict

    ! Force blend: the committed force F_{n+1} enters every later step through af F_n, so the
    ! accepted state is taken one Newton correction beyond the tolerance, with the factor of
    ! the last iteration (one back-substitution and one residual evaluation, no new
    ! factorisation). The correction is kept only when it lowers the residual.
    IF (model%force_blend .AND. nsolve_step > 0 .AND. rnorm > CD_ZERO) THEN
      DO iv = 1, SIZE(model%ws_dq)
        model%ws_dq(iv) = -model%ws_G(iv)
      END DO
      DO i = 1, ndof
        IF (.NOT. model%solve_mask(i)) model%ws_dq(i) = CD_ZERO
      END DO
      CALL CD_Solve_Factored_Banded(model%ws_Keffb, KL_D, KU_D, model%ws_ipiv, model%ws_dq, es, em2, &
                                    factor_validated=.TRUE.)
      IF (es == CD_LINALG_OK) THEN
        CALL transform_rigid_vector_to_global(model, step_rigid_basis, model%ws_dq)
        DO iv = 1, SIZE(model%ws_qtrial)
          model%ws_qtrial(iv) = model%ws_qk(iv) + model%ws_dq(iv)
        END DO
        CALL eval_res(model%ws_qtrial, .TRUE., rnt, es, em2)
        IF (es == CD_HCDYN_OK .AND. rnt < rnorm) THEN
          DO iv = 1, SIZE(model%ws_qk)
            model%ws_qk(iv) = model%ws_qtrial(iv)
          END DO
          DO iv = 1, SIZE(model%ws_G)
            model%ws_G(iv) = model%ws_Gt(iv)
          END DO
          DO iv = 1, SIZE(model%ws_Rk)
            model%ws_Rk(iv) = model%ws_Rt(iv)
          END DO
          rnorm = rnt
          IF (PRESENT(res_out)) res_out = rnorm
        ELSE
          ! Restore the accepted iterate's staged states (the trial overwrote them).
          CALL eval_res(model%ws_qk, .FALSE., rnorm, es, em2)
          IF (es /= CD_HCDYN_OK) THEN
            ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_Step: '//TRIM(em2); RETURN
          END IF
        END IF
      END IF
    END IF

    ! Axial audit of the converged field. ERROR is the hard gate; WARN commits the converged
    ! state and records one event for this accepted step. The audited force is the one that
    ! enters equilibrium, as in the static branch audit: the element-mean axial force (the
    ! three-point Gauss mean of the signed resultant EA (|r'| - 1)), averaged over the element
    ! and its two neighbours weighted by length (clipped at the line ends). The pointwise
    ! trace of a stiff-EA cubic-Hermite element oscillates about that mean, and a single
    ! element mean dips next to a nodal seabed-contact reaction, by amounts that grow with
    ! the element length: judged pointwise, a coarse mesh whose tensions and curvatures are
    ! converged reads as compressed. A compressed span of three or more elements is never
    ! averaged away. The band is tensile_strain_tolerance times the length-weighted EA of the
    ! same three elements; a window whose elements all have a strain lower bound inside the
    ! band cannot fall below it and is skipped.
    IF (model%tensile_mode /= CD_HCDYN_TENSILE_OFF) THEN
      axial_worst_margin = HUGE(CD_ONE)
      axial_worst_force = CD_ZERO; axial_worst_threshold = CD_ZERO
      axial_worst_element = 0; axial_worst_xi = CD_ZERO
      ! Sliding window over elements e-1, e, e+1 (slot 1..3; element 0 = none).
      win_elem = 0
      win_clear = .TRUE.
      win_done = .FALSE.
      win_mean = CD_ZERO
      DO axial_element = 1, model%ne
        IF (axial_element == 1) THEN
          CALL window_load(2, 1)
          IF (model%ne >= 2) CALL window_load(3, 2)
        ELSE
          win_elem(1:2) = win_elem(2:3); win_clear(1:2) = win_clear(2:3)
          win_done(1:2) = win_done(2:3); win_mean(1:2) = win_mean(2:3)
          win_elem(3) = 0; win_clear(3) = .TRUE.; win_done(3) = .FALSE.
          IF (axial_element + 1 <= model%ne) CALL window_load(3, axial_element + 1)
        END IF
        IF (ALL(win_clear)) CYCLE
        win_len = CD_ZERO; win_sum = CD_ZERO; win_ea = CD_ZERO
        DO iv = 1, 3
          IF (win_elem(iv) == 0) CYCLE
          IF (.NOT. win_done(iv)) THEN
            CALL window_mean(iv, axial_es)
            IF (axial_es /= CD_HCABLE_OK) THEN
              ErrStat = CD_HCDYN_NOCONVERGE
              ErrMsg = 'CD_HermiteCable_Dyn_Step: tensile audit failed: '//TRIM(em2)
              IF (PRESENT(tensile_failed)) tensile_failed = .TRUE.
              RETURN
            END IF
          END IF
          win_len = win_len + model%l0(win_elem(iv))
          win_sum = win_sum + model%l0(win_elem(iv))*win_mean(iv)
          win_ea = win_ea + model%l0(win_elem(iv))*model%EA(win_elem(iv))
        END DO
        axial_min = win_sum/win_len
        axial_threshold = model%tensile_strain_tolerance*win_ea/win_len
        axial_margin = axial_min + axial_threshold
        IF (axial_margin < axial_worst_margin) THEN
          axial_worst_margin = axial_margin
          axial_worst_force = axial_min
          axial_worst_threshold = axial_threshold
          axial_worst_element = axial_element
          axial_worst_xi = 0.5_wp   ! the centre of the three-element window (its middle element)
        END IF
      END DO
      IF (axial_worst_margin < CD_ZERO) THEN
        IF (model%tensile_mode == CD_HCDYN_TENSILE_ERROR) THEN
          ErrStat = CD_HCDYN_NOCONVERGE
          BLOCK
            CHARACTER(1024) :: wbuf
            INTEGER :: wios
            wbuf = ''
            WRITE (wbuf, '(A,ES12.4,A,ES12.4,A,I0,A,F7.4,A)', IOSTAT=wios) &
              'CD_HermiteCable_Dyn_Step: axial compression ', axial_worst_force, ' N below -', &
              axial_worst_threshold, ' N (element-mean axial force over elements e-1..e+1) at element ', &
              axial_worst_element, ', xi=', axial_worst_xi, &
              '; increase tensile_strain_tolerance only with an engineering basis, or use warn mode'
            ErrMsg = wbuf
          END BLOCK
          IF (PRESENT(tensile_failed)) tensile_failed = .TRUE.
          RETURN
        END IF
        model%tensile_event_count = model%tensile_event_count + 1
        IF (model%tensile_worst_element == 0 .OR. &
            axial_worst_margin < model%tensile_worst_force + model%tensile_worst_threshold) THEN
          model%tensile_worst_force = axial_worst_force
          model%tensile_worst_threshold = axial_worst_threshold
          model%tensile_worst_element = axial_worst_element
          model%tensile_worst_xi = axial_worst_xi
          model%tensile_worst_time = model%t + dt
        END IF
      END IF
    END IF

    ! Adaptive modified Newton: this step CONVERGED (a non-converged step returned above and
    ! does not count). Accumulate the warmup's full-Newton solve counts and, at the end of
    ! the warmup, latch tangent reuse on for the rest of the run iff the warmup averaged more
    ! than na_threshold solves/step.
    IF (model%newton_adaptive .AND. model%na_steps_done < model%na_warmup) THEN
      model%na_iter_sum = model%na_iter_sum + REAL(nsolve_step, wp)
      model%na_steps_done = model%na_steps_done + 1
      IF (model%na_steps_done == model%na_warmup) THEN
        model%na_mn_active = (model%na_iter_sum/REAL(model%na_warmup, wp)) > model%na_threshold
      END IF
    ELSE IF (model%newton_adaptive) THEN
      model%na_steps_done = model%na_steps_done + 1
    END IF

    ! Commit q_{n+1}, a_{n+1}, v_{n+1}.
    DO iv = 1, SIZE(model%ws_ak)
      model%ws_ak(iv) = c_a*(model%ws_qk(iv) - model%ws_qn(iv)) - c_v*model%ws_vn(iv) - c_a2*model%ws_an(iv)
    END DO
    IF (has_pres) THEN
      DO ip = 1, np
        model%ws_ak(pres_dofs(ip)) = pres_a(ip)
      END DO
    END IF
    DO iv = 1, SIZE(model%q)
      model%q(iv) = model%ws_qk(iv)
    END DO
    DO iv = 1, SIZE(model%a)
      model%a(iv) = model%ws_ak(iv)
    END DO
    DO iv = 1, SIZE(model%v)
      model%v(iv) = model%ws_vn(iv) + dt*((CD_ONE - gm)*model%ws_an(iv) + gm*model%ws_ak(iv))
    END DO
    ! The prescribed DOFs carry their prescribed velocity exactly (not the Newmark blend).
    IF (has_pres) THEN
      DO ip = 1, np
        model%v(pres_dofs(ip)) = pres_v(ip)
      END DO
    END IF
    IF (has_endconn_motion) CALL impose_step_rigid_kinematics(model%q, model%v, model%a)
    IF (model%fr_active) CALL commit_friction_anchors(model)
    model%t = model%t + dt
    IF (has_endconn_motion) model%endconn_d0 = step_endconn_target
    ! Carry the factor of this step's last iteration to the next step.
    model%tr_valid = model%tangent_reuse .AND. factor_ok
    IF (xs_start) THEN
      IF (n_fact > 0 .AND. n_fact >= model%xs_fresh_fact) THEN
        model%xs_backoff = MIN(MAX(1, 2*model%xs_backoff), XS_MAX_BACKOFF)
        model%xs_skip = model%xs_backoff
      ELSE
        model%xs_backoff = 0
      END IF
    ELSE IF (use_ca) THEN
      ! a fresh step on smooth motion: the reference a reuse step must beat
      model%xs_fresh_fact = n_fact
    END IF
    model%tr_dt = dt
    ! The converged iterate's unblended force is F at the committed state: cache it, keyed
    ! on that state, as the next step's F_n.
    IF (model%force_blend) THEN
      DO iv = 1, SIZE(model%fc_R)
        model%fc_R(iv) = model%ws_Rk(iv)
      END DO
      model%fc_q = model%q
      model%fc_v = model%v
      model%fc_t = model%t
      model%fc_d0 = model%endconn_d0
      model%fc_valid = .TRUE.
    END IF
    ! Largest nodal tangent turn of the step (the time-step diagnostic).
    DO i = 1, model%nn
      model%max_step_rotation = MAX(model%max_step_rotation, &
                                    vector_angle(model%ws_qn(6*(i - 1) + 4:6*i), model%q(6*(i - 1) + 4:6*i)))
    END DO
    ! Step wall accumulates on the success path only; error paths abandon their partial time
    ! (the profile is a diagnostic of converged runs).
    IF (prof_enabled) THEN
      CALL SYSTEM_CLOCK(pt1)
      prof_t_step = prof_t_step + REAL(pt1 - pt0, wp)/REAL(pcr, wp)
    END IF

  CONTAINS

    SUBROUTINE window_load(slot, e)
      !! Axial audit window slot for element e: whether its strain lower bound lies inside the
      !! band (its element mean cannot then be below -tolerance*EA); the mean is computed on
      !! demand (window_mean).
      INTEGER, INTENT(IN) :: slot, e
      win_elem(slot) = e
      win_done(slot) = .FALSE.
      win_mean(slot) = CD_ZERO
      axial_qe(1:6) = model%ws_qk(6*(e - 1) + 1:6*e)
      axial_qe(7:12) = model%ws_qk(6*e + 1:6*(e + 1))
      win_clear(slot) = CD_HermiteCable_Axial_Strain_Lower_Bound(axial_qe, model%l0(e)) >= &
        -model%tensile_strain_tolerance
    END SUBROUTINE window_load

    SUBROUTINE window_mean(slot, esx)
      !! Element-mean axial force of the element in window slot: the three-point Gauss mean
      !! of the signed axial resultant (the static branch audit's element mean).
      INTEGER, INTENT(IN) :: slot
      INTEGER, INTENT(OUT) :: esx
      REAL(wp), PARAMETER :: GX(3) = [0.5_wp - 0.5_wp*SQRT(0.6_wp), 0.5_wp, 0.5_wp + 0.5_wp*SQRT(0.6_wp)]
      REAL(wp), PARAMETER :: GW(3) = [5.0_wp/18.0_wp, 8.0_wp/18.0_wp, 5.0_wp/18.0_wp]
      REAL(wp) :: ng
      INTEGER :: e, kg
      e = win_elem(slot)
      axial_qe(1:6) = model%ws_qk(6*(e - 1) + 1:6*e)
      axial_qe(7:12) = model%ws_qk(6*e + 1:6*(e + 1))
      win_mean(slot) = CD_ZERO
      DO kg = 1, 3
        CALL CD_HermiteCable_Axial_Resultant(axial_qe, model%l0(e), model%EA(e), GX(kg), ng, esx, em2)
        IF (esx /= CD_HCABLE_OK) RETURN
        win_mean(slot) = win_mean(slot) + GW(kg)*ng
      END DO
      win_done(slot) = .TRUE.
    END SUBROUTINE window_mean

    LOGICAL FUNCTION force_cache_matches() RESULT(ok)
      !! The cached F_n was evaluated at exactly the committed state and configuration.
      INTEGER :: ii
      ok = model%fc_valid
      IF (.NOT. ok) RETURN
      ok = .NOT. (ABS(model%fc_t - model%t) > CD_ZERO)
      IF (ok) ok = ALL(.NOT. (ABS(model%fc_d0 - model%endconn_d0) > CD_ZERO))
      DO ii = 1, ndof
        IF (.NOT. ok) RETURN
        ok = .NOT. (ABS(model%fc_q(ii) - model%q(ii)) > CD_ZERO .OR. &
                    ABS(model%fc_v(ii) - model%v(ii)) > CD_ZERO)
      END DO
    END FUNCTION force_cache_matches

    SUBROUTINE eval_res(qcur, use_trial, rn, esx, emx)
      !! gen-alpha residual G = M a_alpha + f_int - f_ext, its banded position/velocity
      !! tangents, and the dimensionless per-slot residual norm rn at the trial
      !! configuration qcur. Prescribed DOFs carry their prescribed acceleration/velocity;
      !! the rest follow the Newmark map from qcur. Results land in the model's OWN
      !! buffers -- ws_G (committed) or ws_Gt (line-search trial) selected by use_trial;
      !! evaluations are RESIDUAL-ONLY (the tangent is assembled once per solve iteration
      !! at the loop top, directly into the committed band buffers) -- through pointers into the
      !! host's model, so no model subobject is ever argument-associated alongside the
      !! model itself (which the Fortran aliasing rule would forbid). Uses the ws_* Newmark scratch and the host fscale.
      REAL(wp), INTENT(IN) :: qcur(:)
      LOGICAL, INTENT(IN) :: use_trial
      REAL(wp), INTENT(OUT) :: rn
      INTEGER, INTENT(OUT) :: esx
      CHARACTER(*), INTENT(OUT) :: emx
      REAL(wp), POINTER :: Gout(:)
      INTEGER :: ii, nd
      INTEGER(INT64) :: ec0, ec1
      REAL(wp) :: xref
      REAL(wp), PARAMETER :: RESID_ROUNDOFF = 16.0_wp
      REAL(wp), PARAMETER :: STATE_DIVERGED = 1.0e12_wp
      ec0 = 0_INT64      ! defined even with profiling off (see the step-entry note)
      IF (prof_enabled) THEN
        CALL SYSTEM_CLOCK(ec0)
        prof_n_resid = prof_n_resid + 1
      END IF
      IF (use_trial) THEN
        Gout => model%ws_Gt
      ELSE
        Gout => model%ws_G
      END IF
      ! A diverging Newton iterate (or line-search trial) is rejected before the element
      ! kernels see it: their strain and drag powers of a runaway configuration overflow
      ! (a trapped FPE, or Inf in Release) instead of reporting non-convergence.
      IF (.NOT. CD_All_Finite(qcur)) THEN
        esx = CD_HCDYN_NOCONVERGE
        emx = 'the Newton iterate is not finite'
        rn = CD_ZERO
        RETURN
      END IF
      IF (MAXVAL(ABS(qcur)) > STATE_DIVERGED) THEN
        esx = CD_HCDYN_NOCONVERGE
        emx = 'the Newton iterate diverged (a coordinate beyond 1e12 m)'
        rn = CD_ZERO
        RETURN
      END IF
      IF (.NOT. rigid_rays_valid(model, step_endconn_target, qcur)) THEN
        esx = CD_HCDYN_NOCONVERGE
        emx = 'rigid endpoint tangent left the positive prescribed ray'
        rn = CD_ZERO
        IF (prof_enabled) THEN
          CALL SYSTEM_CLOCK(ec1)
          prof_t_resid = prof_t_resid + REAL(ec1 - ec0, wp)/REAL(pcr, wp)
        END IF
        RETURN
      END IF
      CALL newmark_acceleration(c_a, qcur, model%ws_qn, c_v, model%ws_vn, c_a2, model%ws_an, model%ws_ak)
      IF (has_pres) THEN
        DO ii = 1, np
          model%ws_ak(pres_dofs(ii)) = pres_a(ii)
        END DO
      END IF
      ! Newmark velocity v_{n+1} from the current acceleration, with the prescribed boundary carrying
      ! its exact prescribed velocity (held boundaries have v_n = a_n = 0, so vk = 0 there).
      CALL newmark_velocity(dt, CD_ONE - gm, gm, model%ws_vn, model%ws_an, model%ws_ak, model%ws_vk)
      IF (has_pres) THEN
        DO ii = 1, np
          model%ws_vk(pres_dofs(ii)) = pres_v(ii)
        END DO
      END IF
      IF (has_endconn_motion) CALL impose_step_rigid_kinematics(qcur, model%ws_vk, model%ws_ak)
      ! Force blend: the configuration the force (and the loop-top tangent) is evaluated at is
      ! the iterate itself; only the inertia uses the alpha_m blend. Configuration blend: the
      ! forces are evaluated at the alpha_f-blended state.
      IF (model%force_blend) THEN
        CALL vector_copy(qcur, model%ws_qa)
        CALL vector_copy(model%ws_vk, model%ws_va)
      ELSE
        CALL vector_blend(CD_ONE - af, qcur, af, model%ws_qn, model%ws_qa)
        CALL vector_blend(CD_ONE - af, model%ws_vk, af, model%ws_vn, model%ws_va)
      END IF
      CALL vector_blend(CD_ONE - am, model%ws_ak, am, model%ws_an, model%ws_aa)
      CALL assemble_residual(model, model%ws_qa, model%ws_va, t_alpha, .TRUE., .FALSE., esx, emx, &
                             endconn_direction_eval=step_endconn_eval)
      IF (esx /= CD_HCDYN_OK) THEN
        rn = CD_ZERO
        IF (prof_enabled) THEN
          CALL SYSTEM_CLOCK(ec1)
          prof_t_resid = prof_t_resid + REAL(ec1 - ec0, wp)/REAL(pcr, wp)
        END IF
        RETURN
      END IF
      IF (use_trial) THEN
        CALL vector_copy(model%ws_R, model%ws_Rt)
      ELSE
        CALL vector_copy(model%ws_R, model%ws_Rk)
      END IF
      ! total (structural + frozen added) mass inertia via the band matvec
      CALL band_matvec(model%ws_Msumb, model%ws_aa, model%ws_inert)
      IF (model%force_blend) THEN
        CALL vector_sum3(model%ws_inert, CD_ONE - af, model%ws_R, af, model%fc_R, Gout)
      ELSE
        CALL vector_sum2(model%ws_inert, model%ws_R, Gout)
      END IF
      CALL transform_rigid_vector_to_local(model, step_rigid_basis, Gout)
      ! Freeze the force scale at the first (predictor) evaluation so it is constant across the
      ! step's line-search trials -- see the caller's note. The scale is taken over the FREE
      ! equations only: a fixed/prescribed boundary can carry a large support/contact/inertial
      ! reaction, and including those constrained entries would inflate the scale and let the free
      ! residual look converged prematurely. The convergence verdict is about the equations solved.
      IF (.NOT. scale_set) THEN
        fscale = CD_ONE
        DO ii = 1, ndof
          IF (.NOT. model%solve_mask(ii)) CYCLE
          fscale = MAX(fscale, ABS(model%ws_inert(ii)), ABS(model%ws_R(ii)))
        END DO
        ! Physical force-scale floor: the peak element gravity load |w|*l0. Without it, a step
        ! taken at or near rest (small inertia, a static seed whose f_int - f_ext residual is
        ! already small) collapses fscale to the CD_ONE floor, which turns the RELATIVE tolerance
        ! into an ABSOLUTE one far below the round-off floor of a stiff cable's O(EA) internal
        ! forces -- so the step can never converge however good the seed. The floor only binds when
        ! the motion is gentle; a vigorous step's inertia/residual already dominates it, so the
        ! validated moving-cable gates (whose fscale is inertia/residual-set) are unaffected.
        DO ii = 1, SIZE(model%w)
          fscale = MAX(fscale, ABS(model%w(ii))*model%l0(ii))
        END DO
        ! Round-off floor: an axial force EA*(|r'| - 1) is formed from coordinate differences,
        ! so its rounding error is about eps*|r|/l0*EA per element. On a finely meshed, stiff,
        ! long cable at rest (e.g. 1024 elements, EA ~ 5e8 N, coordinates ~ 1e3 m) that exceeds
        ! tol*|w|*l0 and a held step stagnates at the floor forever. Raise fscale so tol*fscale
        ! is never below RESID_ROUNDOFF times that error; the floor binds only there.
        xref = MAXVAL(model%l0)
        DO ii = 1, ndof
          IF (MOD(ii - 1, 6) < 3) xref = MAX(xref, ABS(qcur(ii)))
        END DO
        fscale = MAX(fscale, RESID_ROUNDOFF*EPSILON(CD_ONE)*MAXVAL(model%EA/model%l0)*xref/tol)
        scale_set = .TRUE.
      END IF
      rn = CD_ZERO
      DO ii = 1, ndof
        IF (.NOT. model%solve_mask(ii)) CYCLE
        IF (MOD(ii - 1, 6) >= 3) THEN                       ! material-tangent slot (m)
          nd = (ii - 1)/6 + 1
          rn = MAX(rn, ABS(Gout(ii))/(fscale*model%trib(nd)))
        ELSE                                                ! translational slot (r)
          rn = MAX(rn, ABS(Gout(ii))/fscale)
        END IF
      END DO
      IF (.NOT. CD_Is_Finite(rn)) THEN
        esx = CD_HCDYN_NOCONVERGE; emx = 'non-finite residual'
      END IF
      IF (prof_enabled) THEN
        CALL SYSTEM_CLOCK(ec1)
        prof_t_resid = prof_t_resid + REAL(ec1 - ec0, wp)/REAL(pcr, wp)
      END IF
    END SUBROUTINE eval_res

    SUBROUTINE impose_step_rigid_kinematics(q_state, v_state, a_state)
      !! Eliminate the two prescribed rotational tangent components while leaving
      !! the axial magnitude as the single solved coordinate. The scalar magnitude
      !! follows the same Newmark scheme as every free DOF; d, d_dot, and d_ddot
      !! are prescribed parent kinematics at the interval end.
      REAL(wp), INTENT(IN) :: q_state(:)
      REAL(wp), INTENT(INOUT) :: v_state(:), a_state(:)
      INTEGER :: iend, base
      REAL(wp) :: magnitude, magnitude_rate, magnitude_acceleration
      DO iend = 1, 2
        IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
        IF (iend == 1) THEN
          base = 3
        ELSE
          base = 6*(model%nn - 1) + 3
        END IF
        magnitude = DOT_PRODUCT(q_state(base + 1:base + 3), step_endconn_target(:, iend))
        magnitude_acceleration = c_a*(magnitude - step_n0(iend)) - &
                                 c_v*step_ndot0(iend) - c_a2*step_nddot0(iend)
        magnitude_rate = step_ndot0(iend) + &
                         dt*((CD_ONE - gm)*step_nddot0(iend) + gm*magnitude_acceleration)
        v_state(base + 1:base + 3) = magnitude_rate*step_endconn_target(:, iend) + &
                                     magnitude*step_endconn_rate(:, iend)
        a_state(base + 1:base + 3) = magnitude_acceleration*step_endconn_target(:, iend) + &
                                     2.0_wp*magnitude_rate*step_endconn_rate(:, iend) + &
                                     magnitude*step_endconn_acceleration(:, iend)
      END DO
    END SUBROUTINE impose_step_rigid_kinematics

    SUBROUTINE add_rigid_transport_tangent()
      !! Complete dG/dn for m = n d. The baseline effective matrix already
      !! contains the c_a*d and gamma/(beta*dt)*d contributions. Add the terms
      !! from d_dot and d_ddot as a rank-one update whose right vector is d;
      !! after the rigid basis transform this affects only the free axial column.
      INTEGER :: iend, base, row, j
      REAL(wp) :: rate_coefficient, extra_acceleration(3)
      rate_coefficient = gm/(bt*dt)
      DO iend = 1, 2
        IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
        IF (iend == 1) THEN
          base = 3
        ELSE
          base = 6*(model%nn - 1) + 3
        END IF
        extra_acceleration = 2.0_wp*rate_coefficient*step_endconn_rate(:, iend) + &
                             step_endconn_acceleration(:, iend)
        model%ws_dq = CD_ZERO
        model%ws_dq(base + 1:base + 3) = extra_acceleration
        CALL band_matvec(model%ws_Msumb, model%ws_dq, model%ws_inert)
        model%ws_inert = (CD_ONE - am)*model%ws_inert
        model%ws_dq = CD_ZERO
        model%ws_dq(base + 1:base + 3) = step_endconn_rate(:, iend)
        CALL band_matvec(model%ws_Kvb, model%ws_dq, model%ws_vk)
        model%ws_inert = model%ws_inert + (CD_ONE - af)*model%ws_vk
        DO j = 1, 3
          DO row = MAX(1, base + j - KU_D), MIN(model%ndof, base + j + KL_D)
            CALL set_dynamic_band_entry(model, model%ws_Keffb, row, base + j, &
                                        dynamic_band_entry(model, model%ws_Keffb, row, base + j) + &
                                        model%ws_inert(row)*step_endconn_target(j, iend))
          END DO
        END DO
      END DO
    END SUBROUTINE add_rigid_transport_tangent

  END SUBROUTINE CD_HermiteCable_Dyn_Step

  LOGICAL FUNCTION rigid_rays_valid(model, directions, state) RESULT(valid)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: directions(3, 2), state(:)
    INTEGER :: iend, base
    REAL(wp) :: magnitude
    valid = .TRUE.
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (iend == 1) THEN
        base = 3
      ELSE
        base = 6*(model%nn - 1) + 3
      END IF
      magnitude = DOT_PRODUCT(state(base + 1:base + 3), directions(:, iend))
      IF (.NOT. CD_Is_Finite(magnitude) .OR. magnitude <= SQRT(TINY(CD_ONE))) THEN
        valid = .FALSE.
        RETURN
      END IF
    END DO
  END FUNCTION rigid_rays_valid

  SUBROUTINE CD_HermiteCable_Dyn_Step_Recovering(model, dt, max_iter, tol, ErrStat, ErrMsg, &
                                                 pres_dofs, pres_q, pres_v, pres_a, max_substeps, &
                                                 endconn_direction, endconn_direction_rate, &
                                                 endconn_direction_acceleration)
    !! Advance one coupling interval of dt, recovering from a non-converged or temporally
    !! under-resolved full step by internal gen-alpha substepping. The common path is ONE call to
    !! CD_HermiteCable_Dyn_Step followed by an O(nn) geometric increment check. A NOCONVERGE return
    !! or an excessive accepted-state increment trips the safety net: restore the step-start
    !! state and re-advance the SAME interval as n substeps of dt/n. Prescribed position, velocity,
    !! and acceleration are joined by one quintic Hermite trajectory, so all three quantities remain
    !! kinematically consistent inside the coupling interval and land exactly on their t+dt targets.
    !!
    !! Subdivision reduces the boundary and load increment seen by each nonlinear solve and sharpens
    !! the Newmark predictor. The host-sampled fluid field and beginning-of-interval added mass remain
    !! fixed within each substep by definition. Landing at t+dt keeps the cable and mooring reactions
    !! at one physical time. Failure at the finest subdivision restores the interval-start state.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT), TARGET :: model
    REAL(wp), INTENT(IN) :: dt, tol
    INTEGER, INTENT(IN) :: max_iter
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: pres_dofs(:)
    REAL(wp), INTENT(IN), OPTIONAL :: pres_q(:), pres_v(:), pres_a(:)
    INTEGER, INTENT(IN), OPTIONAL :: max_substeps
    REAL(wp), INTENT(IN), OPTIONAL :: endconn_direction(3, 2)
    REAL(wp), INTENT(IN), OPTIONAL :: endconn_direction_rate(3, 2), endconn_direction_acceleration(3, 2)
    ! Default subdivision ladder: 4, 16, 64, 256, 1024 substeps. A caller may change the final cap;
    ! the geometric ladder itself is fixed, so the recovery cost of the default ladder is unaffected.
    INTEGER, PARAMETER :: DEFAULT_MAX_SUBSTEPS = 1024
    LOGICAL :: has_pres
    INTEGER :: es, nsub, k, ip, np, recovery_cap
    REAL(wp) :: frac, subdt
    REAL(wp), ALLOCATABLE :: q_save(:), v_save(:), a_save(:)
    REAL(wp), ALLOCATABLE :: q0(:), v0(:), a0(:), pq(:), pv(:), pa(:)
    REAL(wp) :: t_save
    REAL(wp) :: d0_save(3, 2), d0_target(3, 2), d0_sub(3, 2)
    REAL(wp) :: d0_rate_save(3, 2), d0_acceleration_save(3, 2)
    REAL(wp) :: d0_rate_target(3, 2), d0_acceleration_target(3, 2)
    REAL(wp) :: d0_rate_sub(3, 2), d0_acceleration_sub(3, 2)
    ! Adaptive-Newton warmup counters mutate inside every successful CD_HermiteCable_Dyn_Step, so they
    ! are part of the interval-start state that a subdivision rewind must restore: a rung that commits
    ! some substeps then diverges would otherwise leave the next rung (or the final failure) running
    ! from the restored physical state but with polluted na_* counters / a flipped modified-Newton mode.
    INTEGER :: na_steps_save
    REAL(wp) :: na_sum_save
    LOGICAL :: na_mn_save
    INTEGER :: tensile_events_save, tensile_element_save
    REAL(wp) :: tensile_force_save, tensile_threshold_save, tensile_xi_save, tensile_time_save
    ! The step-rotation diagnostic must not keep the turn of a step that is rewound.
    REAL(wp) :: rotation_save
    CHARACTER(300) :: em
    ! Tensile-qualification failure flag reported by each step (classified by flag, not by
    ! message text, so a short caller buffer cannot change the recovery decision).
    LOGICAL :: tensile_fail

    ! Match the public plain-step contract before touching model-owned snapshot workspace. A default
    ! constructed or already-ended model has no allocated state arrays and must fail closed, not trip
    ! an allocatable access runtime error in the recovery wrapper.
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: model not initialised'
      RETURN
    END IF
    IF (model%torsion%active) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: torsion is supported in statics only in this build'
      RETURN
    END IF

    has_pres = PRESENT(pres_dofs) .AND. PRESENT(pres_q) .AND. PRESENT(pres_v) .AND. PRESENT(pres_a)
    ! Reject a partial prescribed-motion set here too (not just in the underlying step): otherwise a
    ! caller that omits e.g. pres_a makes has_pres false and the common-path call below silently holds
    ! the boundary instead of surfacing the same bad input CD_HermiteCable_Dyn_Step rejects.
    IF ((PRESENT(pres_dofs) .OR. PRESENT(pres_q) .OR. PRESENT(pres_v) .OR. PRESENT(pres_a)) &
        .AND. .NOT. has_pres) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: prescribed motion needs pres_dofs/q/v/a together'
      RETURN
    END IF
    np = 0
    IF (has_pres) THEN
      np = SIZE(pres_dofs)
      IF (np < 1) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: pres_dofs must not be empty'
        RETURN
      END IF
    END IF
    recovery_cap = DEFAULT_MAX_SUBSTEPS
    IF (PRESENT(max_substeps)) recovery_cap = max_substeps
    IF (recovery_cap < 4 .OR. recovery_cap > 65536) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: max_substeps must be in [4,65536]'
      RETURN
    END IF

    d0_save = model%endconn_d0
    d0_target = d0_save
    d0_rate_target = CD_ZERO
    d0_acceleration_target = CD_ZERO
    IF (PRESENT(endconn_direction_rate) .NEQV. PRESENT(endconn_direction_acceleration)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: direction rate and acceleration must be supplied together'
      RETURN
    END IF
    IF (PRESENT(endconn_direction_rate) .AND. .NOT. PRESENT(endconn_direction)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: direction derivatives require endconn_direction'
      RETURN
    END IF
    IF (PRESENT(endconn_direction)) THEN
      CALL prepare_endconn_directions(model, endconn_direction, d0_target, ErrStat, ErrMsg)
      IF (ErrStat /= CD_HCDYN_OK) THEN
        ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: '//TRIM(ErrMsg)
        RETURN
      END IF
      CALL prepare_endconn_derivatives(model, d0_target, d0_rate_target, d0_acceleration_target, &
                                       ErrStat, ErrMsg, endconn_direction_rate, endconn_direction_acceleration)
      IF (ErrStat /= CD_HCDYN_OK) THEN
        ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: '//TRIM(ErrMsg)
        RETURN
      END IF
    END IF
    CALL rigid_direction_derivatives_from_state(model, d0_save, model%q, model%v, model%a, &
                                                d0_rate_save, d0_acceleration_save)

    ! Preserve the complete committed state before the full-step attempt in model-owned workspace.
    ! This remains allocation-free on the hot path and makes both genuine Newton failure and a later
    ! quality-gate/subdivision failure atomic. In particular, Step deliberately zeros the history of
    ! fixed DOFs omitted from pres_dofs, so ws_vn/ws_an are not a faithful pre-call snapshot.
    model%ws_recovery_q0 = model%q
    model%ws_recovery_v0 = model%v
    model%ws_recovery_a0 = model%a
    IF (model%fr_active) model%ws_recovery_fr = model%fr_anchor
    t_save = model%t
    na_steps_save = model%na_steps_done; na_sum_save = model%na_iter_sum; na_mn_save = model%na_mn_active
    tensile_events_save = model%tensile_event_count
    tensile_element_save = model%tensile_worst_element
    tensile_force_save = model%tensile_worst_force
    tensile_threshold_save = model%tensile_worst_threshold
    tensile_xi_save = model%tensile_worst_xi
    tensile_time_save = model%tensile_worst_time
    rotation_save = model%max_step_rotation
    ! Common path: one full step. It performs no allocation or snapshot.
    IF (has_pres) THEN
      IF (PRESENT(endconn_direction)) THEN
        CALL CD_HermiteCable_Dyn_Step(model, dt, max_iter, tol, ErrStat, ErrMsg, &
                                      pres_dofs=pres_dofs, pres_q=pres_q, pres_v=pres_v, pres_a=pres_a, &
                                      endconn_direction=d0_target, endconn_direction_rate=d0_rate_target, &
                                      endconn_direction_acceleration=d0_acceleration_target, &
                                      tensile_failed=tensile_fail)
      ELSE
        CALL CD_HermiteCable_Dyn_Step(model, dt, max_iter, tol, ErrStat, ErrMsg, &
                                      pres_dofs=pres_dofs, pres_q=pres_q, pres_v=pres_v, pres_a=pres_a, &
                                      tensile_failed=tensile_fail)
      END IF
    ELSE
      IF (PRESENT(endconn_direction)) THEN
        CALL CD_HermiteCable_Dyn_Step(model, dt, max_iter, tol, ErrStat, ErrMsg, &
                                      endconn_direction=d0_target, endconn_direction_rate=d0_rate_target, &
                                      endconn_direction_acceleration=d0_acceleration_target, &
                                      tensile_failed=tensile_fail)
      ELSE
        CALL CD_HermiteCable_Dyn_Step(model, dt, max_iter, tol, ErrStat, ErrMsg, tensile_failed=tensile_fail)
      END IF
    END IF
    IF (ErrStat == CD_HCDYN_OK) THEN
      IF (temporal_increment_ok(model)) RETURN
      ErrStat = CD_HCDYN_NOCONVERGE
      ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: converged full step exceeded the mesh-aware '// &
               'temporal increment limit'
    ELSE IF (tensile_fail) THEN
      ! Spatial subdivision, not temporal subdivision, is required when the
      ! converged interval-end field fails the tensile qualification.  Step has
      ! not committed its candidate state, so return immediately and preserve
      ! the authoritative interval-start state.
      em = ErrMsg
      ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: tensile qualification failed: '//TRIM(em)
      RETURN
    ELSE IF (ErrStat /= CD_HCDYN_NOCONVERGE) THEN
      RETURN
    END IF

    ! Recovery and final rollback always use the authoritative pre-call snapshot. The full step may
    ! either have committed a quality-rejected state or left model%q/v/a unchanged after Newton
    ! failure; using one anchor for both paths prevents their failure-atomic semantics from drifting.
    ALLOCATE (q_save(model%ndof), v_save(model%ndof), a_save(model%ndof))
    ALLOCATE (q0(np), v0(np), a0(np), pq(np), pv(np), pa(np))
    IF (prof_enabled) prof_n_alloc = prof_n_alloc + 9
    q_save = model%ws_recovery_q0
    v_save = model%ws_recovery_v0
    a_save = model%ws_recovery_a0
    IF (has_pres) THEN
      DO ip = 1, np
        q0(ip) = q_save(pres_dofs(ip))   ! committed boundary at interval start (local frame)
        v0(ip) = v_save(pres_dofs(ip))
        a0(ip) = a_save(pres_dofs(ip))
      END DO
    END IF

    em = ErrMsg
    nsub = 4
    DO
      subdt = dt/REAL(nsub, wp)
      ! Restore the interval-start state (an earlier rung may have committed some substeps before
      ! diverging on a later one) -- including the adaptive-Newton counters, so each rung starts from
      ! the same warmup state the interval began with.
      model%q = q_save; model%v = v_save; model%a = a_save; model%t = t_save
      model%endconn_d0 = d0_save
      IF (model%fr_active) model%fr_anchor = model%ws_recovery_fr
      model%na_steps_done = na_steps_save; model%na_iter_sum = na_sum_save; model%na_mn_active = na_mn_save
      model%tensile_event_count = tensile_events_save
      model%tensile_worst_element = tensile_element_save
      model%tensile_worst_force = tensile_force_save
      model%tensile_worst_threshold = tensile_threshold_save
      model%tensile_worst_xi = tensile_xi_save
      model%tensile_worst_time = tensile_time_save
      model%max_step_rotation = rotation_save
      es = CD_HCDYN_OK
      DO k = 1, nsub
        tensile_fail = .FALSE.
        IF (PRESENT(endconn_direction)) THEN
          frac = REAL(k, wp)/REAL(nsub, wp)
          CALL interpolate_endconn_kinematics(model, d0_save, d0_rate_save, d0_acceleration_save, &
                                              d0_target, d0_rate_target, d0_acceleration_target, dt, frac, &
                                              d0_sub, d0_rate_sub, d0_acceleration_sub, es, em)
          IF (es /= CD_HCDYN_OK) EXIT
        END IF
        IF (has_pres) THEN
          frac = REAL(k, wp)/REAL(nsub, wp)
          CALL quintic_boundary_state(q0, v0, a0, pres_q, pres_v, pres_a, dt, frac, pq, pv, pa)
          ! Avoid endpoint round-off in the polynomial basis. The interior samples remain one C2
          ! trajectory; the final sample is the caller's prescribed state bit for bit.
          IF (k == nsub) THEN
            pq = pres_q; pv = pres_v; pa = pres_a
          END IF
          IF (PRESENT(endconn_direction)) THEN
            CALL CD_HermiteCable_Dyn_Step(model, subdt, max_iter, tol, es, em, &
                                          pres_dofs=pres_dofs, pres_q=pq, pres_v=pv, pres_a=pa, &
                                          endconn_direction=d0_sub, endconn_direction_rate=d0_rate_sub, &
                                          endconn_direction_acceleration=d0_acceleration_sub, &
                                          tensile_failed=tensile_fail)
          ELSE
            CALL CD_HermiteCable_Dyn_Step(model, subdt, max_iter, tol, es, em, &
                                          pres_dofs=pres_dofs, pres_q=pq, pres_v=pv, pres_a=pa, &
                                          tensile_failed=tensile_fail)
          END IF
        ELSE
          IF (PRESENT(endconn_direction)) THEN
            CALL CD_HermiteCable_Dyn_Step(model, subdt, max_iter, tol, es, em, &
                                          endconn_direction=d0_sub, endconn_direction_rate=d0_rate_sub, &
                                          endconn_direction_acceleration=d0_acceleration_sub, &
                                          tensile_failed=tensile_fail)
          ELSE
            CALL CD_HermiteCable_Dyn_Step(model, subdt, max_iter, tol, es, em, tensile_failed=tensile_fail)
          END IF
        END IF
        IF (es == CD_HCDYN_OK) THEN
          IF (.NOT. temporal_increment_ok(model)) THEN
            es = CD_HCDYN_NOCONVERGE
            em = 'substep exceeded mesh-aware temporal increment limit'
          END IF
        END IF
        IF (es /= CD_HCDYN_OK) THEN
          IF (tensile_fail) THEN
            ! Some earlier substeps on this rung may already have committed.
            ! Restore the complete interval-start state before returning.
            model%q = q_save; model%v = v_save; model%a = a_save; model%t = t_save
            model%endconn_d0 = d0_save
            IF (model%fr_active) model%fr_anchor = model%ws_recovery_fr
            model%na_steps_done = na_steps_save
            model%na_iter_sum = na_sum_save
            model%na_mn_active = na_mn_save
            model%tensile_event_count = tensile_events_save
            model%tensile_worst_element = tensile_element_save
            model%tensile_worst_force = tensile_force_save
            model%tensile_worst_threshold = tensile_threshold_save
            model%tensile_worst_xi = tensile_xi_save
            model%tensile_worst_time = tensile_time_save
            model%max_step_rotation = rotation_save
            ErrStat = CD_HCDYN_NOCONVERGE
            ErrMsg = 'CD_HermiteCable_Dyn_Step_Recovering: tensile qualification failed during '// &
                     'internal subdivision: '//TRIM(em)
            RETURN
          END IF
          EXIT
        END IF
      END DO
      IF (es == CD_HCDYN_OK) THEN
        !$OMP ATOMIC UPDATE
        n_substep_recoveries = n_substep_recoveries + 1
        model%recovery_event_count = model%recovery_event_count + 1
        model%recovery_max_substeps_used = MAX(model%recovery_max_substeps_used, nsub)
        ErrStat = CD_HCDYN_OK; ErrMsg = ''
        RETURN
      END IF
      IF (nsub == recovery_cap) EXIT
      IF (nsub > recovery_cap/4) THEN
        nsub = recovery_cap
      ELSE
        nsub = 4*nsub
      END IF
    END DO

    ! Even the finest subdivision diverged: a genuine non-solution (the inexact-tangent residual
    ! plateau at an extreme wave crest, which no step size clears). Fail closed with the interval-start
    ! state restored -- physical state AND adaptive-Newton counters -- so the caller's rollback sees a
    ! consistent step-start.
    model%q = q_save; model%v = v_save; model%a = a_save; model%t = t_save
    model%endconn_d0 = d0_save
    IF (model%fr_active) model%fr_anchor = model%ws_recovery_fr
    model%na_steps_done = na_steps_save; model%na_iter_sum = na_sum_save; model%na_mn_active = na_mn_save
    model%tensile_event_count = tensile_events_save
    model%tensile_worst_element = tensile_element_save
    model%tensile_worst_force = tensile_force_save
    model%tensile_worst_threshold = tensile_threshold_save
    model%tensile_worst_xi = tensile_xi_save
    model%tensile_worst_time = tensile_time_save
    model%max_step_rotation = rotation_save
    ErrStat = CD_HCDYN_NOCONVERGE
    BLOCK
      ! Format into a local buffer: a record longer than the caller's ErrMsg
      ! must truncate, not abort the internal write.
      CHARACTER(1024) :: wbuf
      INTEGER :: wios
      wbuf = ''
      IF (has_pres) THEN
        WRITE (wbuf, '(A,ES12.4,A,ES12.4,A,I0,A,ES12.4,A)', IOSTAT=wios) &
          'Cable recovery failed: t=[', t_save, ', ', t_save + dt, '] s; cap=', recovery_cap, &
          '; max prescribed displacement increment=', &
          MAXVAL(ABS(pres_q - q0)), ': '//TRIM(em)
      ELSE
        WRITE (wbuf, '(A,ES12.4,A,ES12.4,A,I0,A)', IOSTAT=wios) &
          'Cable recovery failed: t=[', t_save, ', ', t_save + dt, '] s; cap=', recovery_cap, &
          ': '//TRIM(em)
      END IF
      ErrMsg = wbuf
    END BLOCK
  END SUBROUTINE CD_HermiteCable_Dyn_Step_Recovering

  SUBROUTINE quintic_boundary_state(q0, v0, a0, q1, v1, a1, duration, u, q, v, a)
    !! C2 trajectory through prescribed position, velocity, and acceleration at both ends of one
    !! coupling interval. u is normalised time in [0,1]; v and a are the exact first and second
    !! time derivatives of q, not independent interpolants.
    REAL(wp), INTENT(IN) :: q0(:), v0(:), a0(:), q1(:), v1(:), a1(:), duration, u
    REAL(wp), INTENT(OUT) :: q(:), v(:), a(:)
    REAL(wp) :: u2, u3, u4, u5, t2

    u2 = u*u; u3 = u2*u; u4 = u3*u; u5 = u4*u
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
  END SUBROUTINE quintic_boundary_state

  SUBROUTINE interpolate_endconn_kinematics(model, d0, v0, a0, d1, v1, a1, duration, fraction, &
                                            d, v, a, ErrStat, ErrMsg)
    !! Interpolate preferred directions during recovery. Finite springs retain
    !! their shortest-arc path; for a rigid end, a Cartesian quintic joins the
    !! complete endpoint states and analytic normalisation restores unit length
    !! plus its first two differential identities at every sample.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: d0(3, 2), v0(3, 2), a0(3, 2), d1(3, 2), v1(3, 2), a1(3, 2)
    REAL(wp), INTENT(IN) :: duration, fraction
    REAL(wp), INTENT(OUT) :: d(3, 2), v(3, 2), a(3, 2)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: iend
    REAL(wp) :: raw_d(3), raw_v(3), raw_a(3), scale, scale_rate, scale_acceleration

    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    CALL interpolate_endconn_directions(model, d0, d1, fraction, d, ErrStat, ErrMsg)
    IF (ErrStat /= CD_HCDYN_OK) RETURN
    v = CD_ZERO
    a = CD_ZERO
    DO iend = 1, 2
      ! A finite spring needs only its preferred direction in the residual and
      ! retains the established shortest-arc interpolation above. The complete
      ! C2 direction kinematics below are required only by a rigid tangent.
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      CALL quintic_boundary_state(d0(:, iend), v0(:, iend), a0(:, iend), &
                                  d1(:, iend), v1(:, iend), a1(:, iend), &
                                  duration, fraction, raw_d, raw_v, raw_a)
      scale = SQRT(DOT_PRODUCT(raw_d, raw_d))
      IF (.NOT. CD_Is_Finite(scale) .OR. scale <= SQRT(TINY(CD_ONE))) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'preferred-direction C2 interpolation is singular'
        RETURN
      END IF
      d(:, iend) = raw_d/scale
      scale_rate = DOT_PRODUCT(d(:, iend), raw_v)
      v(:, iend) = (raw_v - scale_rate*d(:, iend))/scale
      scale_acceleration = (DOT_PRODUCT(raw_v, raw_v) + DOT_PRODUCT(raw_d, raw_a) - &
                            scale_rate*scale_rate)/scale
      a(:, iend) = (raw_a - scale_acceleration*d(:, iend) - 2.0_wp*scale_rate*v(:, iend))/scale
    END DO
    IF (fraction >= CD_ONE) THEN
      d = d1
      v = v1
      a = a1
    END IF
  END SUBROUTINE interpolate_endconn_kinematics

  LOGICAL FUNCTION temporal_increment_ok(model) RESULT(ok)
    !! Cheap accepted-step quality gate. Newton convergence proves equilibrium at its discrete
    !! iterate, but not that one coupling interval resolved the nonlinear trajectory. Reject a step
    !! when the chord of any element (the relative position of its two nodes) changes by more
    !! than 20% of its reference length, or a nodal material-tangent direction turns more than
    !! 15 degrees, in one interval. Both measure deformation: a rigid translation of the line,
    !! however fast (a prescribed hang-off drive on short elements), changes neither and is
    !! integrated exactly by the full step, while a large relative nodal motion or a large
    !! local rotation still trips the gate. A separate
    !! geometry floor rejects a nodal |dr/ds| approaching zero: that state is a singular
    !! centreline parametrisation, not a large but admissible cable strain. The 0.01 floor
    !! still permits 99% compression and therefore acts only as a fail-fast corruption guard,
    !! not as a tensile-only constitutive qualification. The increment limits remain
    !! dimensionless and mesh-aware; they impose no cable-specific curvature or load ceiling.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), PARAMETER :: POS_FRAC_MAX = 0.20_wp
    REAL(wp), PARAMETER :: COS_TURN_MAX = 0.9659258262890683_wp ! cos(15 deg)
    REAL(wp), PARAMETER :: NORM_FLOOR = 1.0e-12_wp
    REAL(wp), PARAMETER :: MIN_MATERIAL_STRETCH = 1.0e-2_wp
    INTEGER :: i, gi
    REAL(wp) :: dr2, n0, n1, cosine

    ok = .TRUE.
    DO i = 1, model%ne
      gi = 6*(i - 1)
      ! Change of the element chord r_(i+1) - r_i over the interval.
      dr2 = SUM(((model%q(gi + 7:gi + 9) - model%q(gi + 1:gi + 3)) - &
                 (model%ws_qn(gi + 7:gi + 9) - model%ws_qn(gi + 1:gi + 3)))**2)
      IF (.NOT. (dr2 <= (POS_FRAC_MAX*model%l0(i))**2)) THEN
        ok = .FALSE.; RETURN
      END IF
    END DO
    DO i = 1, model%nn
      gi = 6*(i - 1)
      n0 = SQRT(SUM(model%ws_qn(gi + 4:gi + 6)**2))
      n1 = SQRT(SUM(model%q(gi + 4:gi + 6)**2))
      IF (.NOT. CD_Is_Finite(n0) .OR. .NOT. CD_Is_Finite(n1) .OR. &
          n0 <= NORM_FLOOR .OR. n1 < MIN_MATERIAL_STRETCH) THEN
        ok = .FALSE.; RETURN
      END IF
      cosine = DOT_PRODUCT(model%ws_qn(gi + 4:gi + 6), model%q(gi + 4:gi + 6))/(n0*n1)
      IF (cosine < COS_TURN_MAX) THEN
        ok = .FALSE.; RETURN
      END IF
    END DO
  END FUNCTION temporal_increment_ok

  PURE FUNCTION CD_HermiteCable_Dyn_Recovery_Count() RESULT(n)
    !! Coupling intervals that completed only via internal substepping since the last
    !! CD_HermiteCable_Dyn_Recovery_Reset (the finite-EI analogue of CD_System_Fallback_Count).
    INTEGER :: n
    n = n_substep_recoveries
  END FUNCTION CD_HermiteCable_Dyn_Recovery_Count

  SUBROUTINE CD_HermiteCable_Dyn_Recovery_Reset()
    !! Zero the substep-recovery tripwire (start of a measured window).
    n_substep_recoveries = 0
  END SUBROUTINE CD_HermiteCable_Dyn_Recovery_Reset

  SUBROUTINE CD_HermiteCable_Dyn_Energy(model, kinetic, strain, ErrStat, ErrMsg, connection)
    !! Kinetic energy 1/2 v^T M v with the consistent structural mass plus any discrete
    !! attachment's mass and added mass Ca rho V (the line's distributed added mass is fluid
    !! kinetic energy, excluded so the conservation diagnostics track the structure) and total
    !! elastic energy at the current state. The optional
    !! connection output reports the end-spring share already included in strain.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: kinetic, strain
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(OUT), OPTIONAL :: connection
    INTEGER :: e, es
    REAL(wp) :: qe(12), Ee, fe(12), Ke(12, 12), connection_energy, end_energy
    CHARACTER(200) :: em2
    ErrStat = CD_HCDYN_OK; ErrMsg = ''; kinetic = CD_ZERO; strain = CD_ZERO
    connection_energy = CD_ZERO
    IF (PRESENT(connection)) connection = CD_ZERO
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Energy: model not initialised'; RETURN
    END IF
    BLOCK
      REAL(wp) :: Mv(model%ndof)
      CALL band_matvec(model%Mgb, model%v, Mv)
      kinetic = 0.5_wp*DOT_PRODUCT(model%v, Mv)
    END BLOCK
    DO e = 1, model%ne
      CALL gather_elem(model, e, qe)
      CALL CD_HermiteCable_Element(qe, model%l0(e), model%EA(e), model%EI(e), Ee, fe, Ke, es, em2, &
                                   axial_quadrature_order=model%axial_quadrature_order, &
                                   bending_quadrature_order=model%bending_quadrature_order)
      IF (es /= CD_HCABLE_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Dyn_Energy: '//TRIM(em2); RETURN
      END IF
      strain = strain + Ee
    END DO
    IF (model%has_endconn) THEN
      IF (model%endconn_k(1) > CD_ZERO) THEN
        CALL CD_EndConn_Energy(model%q(4:6), model%endconn_d0(:, 1), model%endconn_k(1), &
                               end_energy, es, em2)
        IF (es /= CD_ENDCONN_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE
          ErrMsg = 'CD_HermiteCable_Dyn_Energy: node 1 end connection: '//TRIM(em2)
          RETURN
        END IF
        connection_energy = connection_energy + end_energy
      END IF
      IF (model%endconn_k(2) > CD_ZERO) THEN
        CALL CD_EndConn_Energy(model%q(6*(model%nn - 1) + 4:6*(model%nn - 1) + 6), &
                               model%endconn_d0(:, 2), model%endconn_k(2), end_energy, es, em2)
        IF (es /= CD_ENDCONN_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE
          ErrMsg = 'CD_HermiteCable_Dyn_Energy: final-node end connection: '//TRIM(em2)
          RETURN
        END IF
        connection_energy = connection_energy + end_energy
      END IF
      strain = strain + connection_energy
    END IF
    IF (model%torsion%active) THEN
      BLOCK
        REAL(wp) :: th, mt, gq(model%ndof), eg(6)
        CALL CD_HermiteCable_Dyn_Torsion_State(model, th, mt, ErrStat, em2, gq, eg)
        IF (ErrStat /= CD_HCDYN_OK) THEN
          ErrMsg = 'CD_HermiteCable_Dyn_Energy: '//TRIM(em2)
          RETURN
        END IF
        strain = strain + 0.5_wp*CD_HermiteTorsion_Compliance(model%torsion, model%l0)*mt*mt
      END BLOCK
    END IF
    IF (PRESENT(connection)) connection = connection_energy
  END SUBROUTINE CD_HermiteCable_Dyn_Energy

  SUBROUTINE CD_HermiteCable_Dyn_Curvature(model, curv_out, ErrStat, ErrMsg)
    !! Per-node centreline curvature at the current state (element-midpoint / endpoint maxima).
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: curv_out(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: a, es
    REAL(wp) :: qe(12), kL, kR
    CHARACTER(200) :: em2
    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Curvature: model not initialised'; RETURN
    END IF
    IF (SIZE(curv_out) /= model%nn) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Curvature: curv_out must be nn'; RETURN
    END IF
    DO a = 1, model%nn
      curv_out(a) = CD_ZERO
      ! A degenerate nodal tangent makes the curvature undefined; propagate that failure rather than
      ! silently reporting zero curvature -- a fatigue assessment reads this channel, so a hidden
      ! zero at a degenerate state would understate the peak sag-bend curvature.
      IF (a <= model%ne) THEN
        CALL gather_elem(model, a, qe)
        CALL CD_HermiteCable_Curvature(qe, model%l0(a), CD_ZERO, kR, es, em2)
        IF (es /= CD_HCABLE_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Dyn_Curvature: '//TRIM(em2); RETURN
        END IF
        curv_out(a) = kR
      END IF
      IF (a >= 2) THEN
        CALL gather_elem(model, a - 1, qe)
        CALL CD_HermiteCable_Curvature(qe, model%l0(a - 1), CD_ONE, kL, es, em2)
        IF (es /= CD_HCABLE_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Dyn_Curvature: '//TRIM(em2); RETURN
        END IF
        curv_out(a) = MAX(curv_out(a), kL)
      END IF
    END DO
  END SUBROUTINE CD_HermiteCable_Dyn_Curvature

  SUBROUTINE CD_HermiteCable_Dyn_Reaction(model, node, freac, ErrStat, ErrMsg)
    !! Constraint reaction at a fixed/prescribed node: the force the LINE exerts ON its
    !! support (fairlead/anchor) at the current committed state, in global coordinates.
    !!
    !! At a fixed translational DOF the equation of motion is not solved; the support
    !! supplies whatever generalized force closes it. With G = M a + f_int - f_ext
    !! evaluated at the committed (q, v, a, t), the support applies +G(dof) TO the line,
    !! so the line applies -G(dof) to the support -- the fairlead load an outer solver
    !! (the OpenFAST FMF coupling) needs. The total mass (structural + added) is
    !! REASSEMBLED here at the committed configuration from the model's current public
    !! fields rather than read from the ws_Msumb scratch: that scratch is only "last
    !! populated by Init / Step / setter refreshes", and a failed transactional hydro
    !! reconfiguration rolls the public fields back while leaving it partially
    !! refreshed -- a stale-mass hazard this output must not inherit.
    !! Meaningful at nodes whose translations are fixed or prescribed; at a free node the
    !! same expression returns the residual (~0 at convergence), not a support load.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER, INTENT(IN) :: node
    REAL(wp), INTENT(OUT) :: freac(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: es, gi
    CHARACTER(200) :: em2
    ErrStat = CD_HCDYN_OK; ErrMsg = ''; freac = CD_ZERO
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Reaction: model not initialised'; RETURN
    END IF
    IF (node < 1 .OR. node > model%nn) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_Reaction: node out of range'; RETURN
    END IF
    ! Residual-only assembly at the committed state: no tangent buffer is touched, so the
    ! committed K/Kv (which callers may still hold between steps) are preserved.
    CALL assemble_residual(model, model%q, model%v, model%t, .TRUE., .FALSE., es, em2)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_Reaction: '//TRIM(em2); RETURN
    END IF
    ! rebuild the total mass from the CURRENT fields at the committed configuration
    model%ws_Msumb = model%Mgb
    IF (model%has_am) THEN
      CALL assemble_added_mass(model, model%q, es, em2)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_Reaction: '//TRIM(em2); RETURN
      END IF
      model%ws_Msumb = model%ws_Msumb + model%ws_Mamb
    END IF
    CALL band_matvec(model%ws_Msumb, model%a, model%ws_inert)
    gi = 6*(node - 1)
    freac = -(model%ws_inert(gi + 1:gi + 3) + model%ws_R(gi + 1:gi + 3))
  END SUBROUTINE CD_HermiteCable_Dyn_Reaction

  SUBROUTINE CD_HermiteCable_Dyn_End_Force(model, node, force, ErrStat, ErrMsg)
    !! Force the line exerts on the point holding its end node (node = 1 or nn) at the
    !! committed state (q, v): the end element's internal force (elastic plus Kelvin-Voigt
    !! axial damping) plus the end node's share of the distributed external loads --
    !! consistent submerged weight, dry-part buoyancy above the free surface, seabed contact
    !! normal force with its damping, and the drag of the ambient fluid (current, waves or
    !! held host field) at the actual relative velocity, and the weight, buoyancy,
    !! Froude-Krylov force and drag of any discrete attachment on the end node -- i.e.
    !! f_ext - f_int at the node. This is the actual force on the attachment, the EI = 0 lines' CD_Get_Model_EndForces
    !! quantity; inertia (and seabed friction) stay in the coupled load
    !! (CD_HermiteCable_Dyn_Reaction). At rest it is the static end reaction. Read-only: no
    !! model buffer is written.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: node
    REAL(wp), INTENT(OUT) :: force(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), PARAMETER :: PI_L = 3.14159265358979323846_wp
    REAL(wp) :: qe(12), ve(12), fe(12), fdrag(12), fdamp(12), energy, fdry(4), kdry(4, 4), r3(3), gate, gate_gap
    REAL(wp) :: dry_wl, dry_b, dry_r, gap, normal, normal_gap, z_floor, xg, yg, gxg, gyg
    REAL(wp) :: gx, gy, sfac, nvec(3), vn, rk(3), jrel(3, 3), jtan(3, 3)
    INTEGER :: e, i0, iz, zslot, es, k
    LOGICAL :: dry, has_jac
    CHARACTER(200) :: em2

    ErrStat = CD_HCDYN_OK; ErrMsg = ''; force = CD_ZERO; gap = -HUGE(CD_ONE)
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_End_Force: model not initialised'; RETURN
    END IF
    IF (node /= 1 .AND. node /= model%nn) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Dyn_End_Force: node must be a line end'; RETURN
    END IF
    IF (node == 1) THEN
      e = 1; i0 = 0; zslot = 1
    ELSE
      e = model%ne; i0 = 6; zslot = 3
    END IF
    qe(1:6) = model%q(6*(e - 1) + 1:6*e)
    qe(7:12) = model%q(6*e + 1:6*(e + 1))
    ve(1:6) = model%v(6*(e - 1) + 1:6*e)
    ve(7:12) = model%v(6*e + 1:6*(e + 1))
    ! R = f_int - f_ext at the end node's translations, accumulated below.
    CALL CD_HermiteCable_Element(qe, model%l0(e), model%EA(e), model%EI(e), energy, fe, ErrStat=es, ErrMsg=em2, &
                                 axial_quadrature_order=model%axial_quadrature_order, &
                                 bending_quadrature_order=model%bending_quadrature_order)
    IF (es /= CD_HCABLE_OK) THEN
      ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Dyn_End_Force: '//TRIM(em2); RETURN
    END IF
    r3 = fe(i0 + 1:i0 + 3)
    IF (model%has_axial_damping) THEN
      CALL CD_HermiteCable_Axial_Damping_Element(qe, ve, model%l0(e), model%BA(e), fdamp, ErrStat=es, ErrMsg=em2, &
                                                 quadrature_order=model%axial_quadrature_order)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Dyn_End_Force: '//TRIM(em2); RETURN
      END IF
      r3 = r3 + fdamp(i0 + 1:i0 + 3)
    END IF
    r3(3) = r3(3) + model%w(e)*0.5_wp*model%l0(e)
    IF (model%has_drag .OR. model%has_am) THEN
      IF (model%has_drag) THEN
        dry_wl = model%waterline_z
        dry_b = model%rho_w*model%gravity*0.25_wp*PI_L*model%diam(e)**2
        dry_r = 0.5_wp*model%diam(e)
      ELSE
        dry_wl = model%am_wl
        dry_b = model%am_rho*model%gravity*0.25_wp*PI_L*model%am_diam(e)**2
        dry_r = 0.5_wp*model%am_diam(e)
      END IF
      IF (model%has_held_field) dry_wl = 0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1))
      CALL CD_HermiteCable_Dry_Buoyancy(qe([3, 6, 9, 12]), model%l0(e), dry_b, dry_wl, fdry, kdry, dry, &
                                        radius=dry_r)
      IF (dry) r3(3) = r3(3) + fdry(zslot)
    END IF
    iz = 6*(node - 1) + 3
    IF (model%has_contact) THEN
      z_floor = model%seabed_z
      gx = CD_ZERO
      gy = CD_ZERO
      IF (model%has_contact_bathymetry) THEN
        xg = model%contact_frame_c*model%q(iz - 2) - model%contact_frame_s*model%q(iz - 1)
        yg = model%contact_frame_s*model%q(iz - 2) + model%contact_frame_c*model%q(iz - 1)
        CALL CD_Bathymetry_Floor_Gradient(model%contact_bathymetry, xg, yg, z_floor, gxg, gyg, es, em2)
        IF (es /= CD_BATHY_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Dyn_End_Force: '//TRIM(em2); RETURN
        END IF
        gx = model%contact_frame_c*gxg + model%contact_frame_s*gyg
        gy = -model%contact_frame_s*gxg + model%contact_frame_c*gyg
      END IF
      ! frictionless: along the upward floor normal, at the normal penetration and velocity
      sfac = SQRT(CD_ONE + gx*gx + gy*gy)
      nvec = [-gx, -gy, CD_ONE]/sfac
      gap = (z_floor - model%q(iz))/sfac
      CALL CD_Seabed_Normal_Law(gap, model%contact_kn(node), normal, normal_gap)
      vn = DOT_PRODUCT(model%v(iz - 2:iz), nvec)
      IF (ALLOCATED(model%contact_cn) .AND. vn < CD_ZERO) THEN
        CALL CD_Seabed_Contact_Gate(gap, gate, gate_gap)
        normal = normal - model%contact_cn(node)*gate*vn
      END IF
      r3 = r3 - normal*nvec
    ELSE IF (model%kn > CD_ZERO) THEN
      gap = model%seabed_z - model%q(iz)
      IF (gap > CD_ZERO) r3(3) = r3(3) - model%kn*gap*model%trib(node)
    END IF
    IF (model%has_drag) THEN
      IF (model%has_wave .AND. model%wv_ncomp > 0) THEN
        CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                          model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                          fdrag, ErrStat=es, ErrMsg=em2, &
                                          wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=model%t, &
                                          wvc_amp=model%wv_ampc, wvc_om=model%wv_omc, &
                                          wvc_k=model%wv_kc, wvc_ph=model%wv_phc, &
                                          wvc_e2kd=model%wv_e2kd)
      ELSE IF (model%has_wave) THEN
        CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                          model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                          fdrag, ErrStat=es, ErrMsg=em2, &
                                          wv_h=model%wv_h, wv_om=model%wv_om, wv_k=model%wv_k, &
                                          wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=model%t)
      ELSE IF (model%has_held_field) THEN
        CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), &
                                          model%current + 0.5_wp*(model%hf_u(:, e) + model%hf_u(:, e + 1)), &
                                          0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1)), &
                                          model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                          fdrag, ErrStat=es, ErrMsg=em2)
      ELSE
        CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                          model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                          fdrag, ErrStat=es, ErrMsg=em2)
      END IF
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_End_Force: '//TRIM(em2); RETURN
      END IF
      r3 = r3 - fdrag(i0 + 1:i0 + 3)
    END IF
    IF (model%has_attach) THEN
      ! A clump or buoy on the end node loads the point holding it, exactly as it loads the
      ! end node in the residual (assemble_attachments).
      DO k = 1, SIZE(model%att_node)
        IF (model%att_node(k) /= node) CYCLE
        CALL attachment_residual(model, k, model%q, model%v, .FALSE., rk, has_jac, jrel, jtan, es, em2)
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_End_Force: '//TRIM(em2); RETURN
        END IF
        r3 = r3 + rk
      END DO
    END IF
    force = -r3
    ! Condensed torsion: R carries -M_t dTheta/dq, so the end force gains +M_t dTheta/dr_end.
    IF (model%torsion%active) THEN
      BLOCK
        REAL(wp) :: th, mt, gq(model%ndof), eg(6)
        CALL CD_HermiteCable_Dyn_Torsion_State(model, th, mt, es, em2, gq, eg)
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = es; ErrMsg = 'CD_HermiteCable_Dyn_End_Force: '//TRIM(em2); RETURN
        END IF
        force = force + mt*gq(6*(node - 1) + 1:6*(node - 1) + 3)
      END BLOCK
    END IF
    ! An end resting on the seabed (at or inside the contact blend) does not carry the end
    ! node's weight: the floor does, as it does for the grounded nodes beside it. The
    ! downward part of the end force is left to the seabed, so a grounded end reports the
    ! horizontal (grounded-run) tension.
    IF ((model%has_contact .OR. model%kn > CD_ZERO) .AND. .NOT. (gap <= -CD_SEABED_CONTACT_BLEND)) &
      force(3) = MAX(force(3), CD_ZERO)
  END SUBROUTINE CD_HermiteCable_Dyn_End_Force

  PURE REAL(wp) FUNCTION CD_HermiteCable_Dyn_Max_Step_Rotation(model) RESULT(angle)
    !! Largest angle [rad] through which a nodal tangent turned over one committed step
    !! since initialisation: the time-step resolution of the tangent rotation.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    angle = model%max_step_rotation
  END FUNCTION CD_HermiteCable_Dyn_Max_Step_Rotation

  PURE REAL(wp) FUNCTION vector_angle(a, b) RESULT(angle)
    !! Angle between two vectors (zero when either vanishes), robust near 0 and pi.
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    angle = ATAN2(NORM2(c), DOT_PRODUCT(a, b))
  END FUNCTION vector_angle

  SUBROUTINE CD_HermiteCable_Dyn_EndConnection_Moment(model, node, moment, ErrStat, ErrMsg)
    !! Moment exerted by the line-end connection on its supporting body. Finite
    !! mode evaluates the spring force; rigid mode recovers the exact constraint
    !! reaction from the committed equation of motion. Pinned mode returns zero.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER, INTENT(IN) :: node
    REAL(wp), INTENT(OUT) :: moment(3)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: iend, base, es
    REAL(wp) :: f_ec(3), k_ec(3, 3), direction(3), bend_moment(3)
    CHARACTER(200) :: em

    moment = CD_ZERO
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_EndConnection_Moment: model not initialised'
      RETURN
    END IF
    IF (node == 1) THEN
      iend = 1
    ELSE IF (node == model%nn) THEN
      iend = 2
    ELSE
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_EndConnection_Moment: node must be a cable endpoint'
      RETURN
    END IF
    IF (model%torsion%active) THEN
      ! Torque of the condensed torsion on the support: the generalised moment
      ! -dE_t/dw = M_t dTheta/dw of the end frame (its axial part is -+M_t d); the slaved
      ! tangent of a rigid end adds its share through the constraint reaction below.
      BLOCK
        REAL(wp) :: th, mt, gq(model%ndof), eg(6)
        CALL CD_HermiteCable_Dyn_Torsion_State(model, th, mt, es, em, gq, eg)
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = es
          ErrMsg = 'CD_HermiteCable_Dyn_EndConnection_Moment: '//TRIM(em)
          RETURN
        END IF
        moment = mt*eg(3*iend - 2:3*iend)
      END BLOCK
    END IF
    IF (.NOT. model%has_endconn) RETURN
    IF (model%endconn_mode(iend) == CD_ENDCONN_PINNED) RETURN

    base = 6*(node - 1) + 3
    IF (model%endconn_mode(iend) == CD_ENDCONN_FINITE) THEN
      CALL CD_EndConn_Spring(model%q(base + 1:base + 3), model%endconn_d0(:, iend), &
                             model%endconn_k(iend), f_ec, k_ec, es, em)
      IF (es /= CD_ENDCONN_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE
        ErrMsg = 'CD_HermiteCable_Dyn_EndConnection_Moment: '//TRIM(em)
        RETURN
      END IF
    ELSE
      ! G0 is the unconstrained dynamic residual. The constraint contribution
      ! required in the residual is -G0 in the transverse tangent plane; under
      ! virtual work, m cross that contribution is the line's moment on the
      ! supporting body, matching the finite-spring convention.
      CALL assemble_residual(model, model%q, model%v, model%t, .TRUE., .FALSE., es, em)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es
        ErrMsg = 'CD_HermiteCable_Dyn_EndConnection_Moment: '//TRIM(em)
        RETURN
      END IF
      model%ws_Msumb = model%Mgb
      IF (model%has_am) THEN
        CALL assemble_added_mass(model, model%q, es, em)
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = es
          ErrMsg = 'CD_HermiteCable_Dyn_EndConnection_Moment: '//TRIM(em)
          RETURN
        END IF
        model%ws_Msumb = model%ws_Msumb + model%ws_Mamb
      END IF
      CALL band_matvec(model%ws_Msumb, model%a, model%ws_inert)
      f_ec = -(model%ws_inert(base + 1:base + 3) + model%ws_R(base + 1:base + 3))
      direction = model%endconn_d0(:, iend)
      f_ec = f_ec - DOT_PRODUCT(f_ec, direction)*direction
    END IF
    CALL CD_EndConn_Reaction_Moment(model%q(base + 1:base + 3), f_ec, bend_moment, es, em)
    IF (es /= CD_ENDCONN_OK) THEN
      ErrStat = CD_HCDYN_NOCONVERGE
      ErrMsg = 'CD_HermiteCable_Dyn_EndConnection_Moment: '//TRIM(em)
      moment = CD_ZERO
      RETURN
    END IF
    moment = moment + bend_moment
  END SUBROUTINE CD_HermiteCable_Dyn_EndConnection_Moment

  SUBROUTINE CD_HermiteCable_Dyn_Set_Torsion(model, torsion, ErrStat, ErrMsg)
    !! Install the condensed torsion of a converged static solve (CD_HermiteCable_Static_Solve
    !! with torsion): the end frames and GJ in the model's solve frame, Phi and the accepted
    !! Theta. The residual, the support reaction, the end force, the connection moment (now
    !! including the torque) and the energy then include E_t = (Phi - Theta)**2/(2 C); the
    !! acceleration is recomputed from the complete load set. An inactive description removes
    !! the torsion. The time step refuses a model with torsion (statics only in this build).
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    TYPE(CD_HermiteTorsionType), INTENT(IN) :: torsion
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    TYPE(CD_HermiteTorsionType) :: saved
    REAL(wp) :: th, mt, eg(6)
    REAL(wp), ALLOCATABLE :: gq(:)
    INTEGER :: es
    CHARACTER(200) :: em
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Torsion: model not initialised'
      RETURN
    END IF
    saved = model%torsion
    IF (.NOT. torsion%active) THEN
      model%torsion = CD_HermiteTorsionType()
    ELSE
      CALL CD_HermiteTorsion_Validate(torsion, model%l0, es, em)
      IF (es /= CD_HTORS_OK) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Torsion: '//TRIM(em)
        RETURN
      END IF
      IF (.NOT. torsion%has_theta) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Torsion: the torsion state has no accepted Theta (solve the statics first)'
        RETURN
      END IF
      model%torsion = torsion
      ALLOCATE (gq(model%ndof))
      CALL CD_HermiteCable_Dyn_Torsion_State(model, th, mt, es, em, gq, eg)
      IF (es /= CD_HCDYN_OK) THEN
        model%torsion = saved
        ErrStat = es
        ErrMsg = 'CD_HermiteCable_Dyn_Set_Torsion: '//TRIM(em)
        RETURN
      END IF
      model%torsion%theta = th
      model%torsion%torque = mt
    END IF
    model%snap_torsion_theta = model%torsion%theta
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    CALL CD_HermiteCable_Dyn_Recompute_Acceleration(model, es, em)
    IF (es /= CD_HCDYN_OK) THEN
      model%torsion = saved
      ErrStat = es
      ErrMsg = 'CD_HermiteCable_Dyn_Set_Torsion: '//TRIM(em)
    END IF
  END SUBROUTINE CD_HermiteCable_Dyn_Set_Torsion

  SUBROUTINE CD_HermiteCable_Dyn_Torsion_State(model, theta, torque, ErrStat, ErrMsg, grad, end_grad, q_eval)
    !! Theta (unwrapped against the stored accepted value), the torque M_t = (Phi - Theta)/C and,
    !! optionally, dTheta/dq and dTheta/d[w_1, w_2] at the committed state (or at q_eval). Read
    !! only: the stored Theta is not advanced. A change of more than pi/2 from the stored value
    !! fails closed (a 2 pi branch could otherwise be lost).
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(OUT) :: theta, torque
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(OUT), OPTIONAL :: grad(:), end_grad(6)
    REAL(wp), INTENT(IN), OPTIONAL :: q_eval(:)
    REAL(wp), ALLOCATABLE :: g(:)
    REAL(wp) :: raw, eg(6), c
    INTEGER :: es
    CHARACTER(200) :: em
    theta = CD_ZERO
    torque = CD_ZERO
    IF (PRESENT(grad)) grad = CD_ZERO
    IF (PRESENT(end_grad)) end_grad = CD_ZERO
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. model%torsion%active) RETURN
    ALLOCATE (g(model%ndof))
    IF (PRESENT(q_eval)) THEN
      CALL CD_HermiteTorsion_Line(q_eval, model%l0, model%torsion%ends, raw, g, es, em, end_grad=eg, &
                                  quadrature_order=model%torsion%quadrature_order)
    ELSE
      CALL CD_HermiteTorsion_Line(model%q, model%l0, model%torsion%ends, raw, g, es, em, end_grad=eg, &
                                  quadrature_order=model%torsion%quadrature_order)
    END IF
    IF (es /= CD_HTORS_OK) THEN
      ErrStat = CD_HCDYN_NOCONVERGE
      ErrMsg = 'torsion: '//TRIM(em)
      RETURN
    END IF
    theta = CD_HermiteTorsion_Unwrap(raw, model%torsion%theta)
    IF (ABS(theta - model%torsion%theta) > CD_HTORS_MAX_STEP) THEN
      ErrStat = CD_HCDYN_NOCONVERGE
      ErrMsg = 'torsion: the twist moved more than pi/2 from its accepted value'
      theta = CD_ZERO
      RETURN
    END IF
    c = CD_HermiteTorsion_Compliance(model%torsion, model%l0)
    torque = (model%torsion%phi - theta)/c
    IF (PRESENT(grad)) grad = g
    IF (PRESENT(end_grad)) end_grad = eg
  END SUBROUTINE CD_HermiteCable_Dyn_Torsion_State

  SUBROUTINE CD_HermiteCable_Dyn_Snapshot(model, ErrStat, ErrMsg)
    !! Capture the committed step state (q, v, a, t) into the model-owned snapshot
    !! buffers (see the type note: the aggregate stage-then-commit contract). Buffers
    !! allocate on first use and are reused thereafter.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Snapshot: model not initialised'
      RETURN
    END IF
    IF (.NOT. ALLOCATED(model%snap_q)) THEN
      ALLOCATE (model%snap_q(model%ndof), model%snap_v(model%ndof), model%snap_a(model%ndof))
    END IF
    model%snap_q = model%q
    model%snap_v = model%v
    model%snap_a = model%a
    model%snap_t = model%t
    model%snap_endconn_d0 = model%endconn_d0
    model%snap_torsion_theta = model%torsion%theta
    IF (model%fr_active) model%snap_fr_anchor = model%fr_anchor
    model%snap_na_steps_done = model%na_steps_done
    model%snap_na_iter_sum = model%na_iter_sum
    model%snap_na_mn_active = model%na_mn_active
    model%snap_tensile_event_count = model%tensile_event_count
    model%snap_tensile_worst_element = model%tensile_worst_element
    model%snap_tensile_worst_force = model%tensile_worst_force
    model%snap_tensile_worst_threshold = model%tensile_worst_threshold
    model%snap_tensile_worst_xi = model%tensile_worst_xi
    model%snap_tensile_worst_time = model%tensile_worst_time
    model%snap_max_step_rotation = model%max_step_rotation
    model%snap_recovery_event_count = model%recovery_event_count
    model%snap_recovery_max_substeps_used = model%recovery_max_substeps_used
    model%snap_valid = .TRUE.
  END SUBROUTINE CD_HermiteCable_Dyn_Snapshot

  SUBROUTINE CD_HermiteCable_Dyn_Restore(model, ErrStat, ErrMsg)
    !! Restore the committed step state from the last Snapshot. Fails closed if no valid
    !! snapshot exists (a restore must never silently run on stale or absent state).
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    IF (.NOT. model%initialized) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Restore: model not initialised'
      RETURN
    END IF
    IF (.NOT. model%snap_valid) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Dyn_Restore: no valid snapshot'
      RETURN
    END IF
    model%q = model%snap_q
    model%v = model%snap_v
    model%a = model%snap_a
    model%t = model%snap_t
    model%endconn_d0 = model%snap_endconn_d0
    model%torsion%theta = model%snap_torsion_theta
    IF (model%fr_active) model%fr_anchor = model%snap_fr_anchor
    model%na_steps_done = model%snap_na_steps_done
    model%na_iter_sum = model%snap_na_iter_sum
    model%na_mn_active = model%snap_na_mn_active
    model%tensile_event_count = model%snap_tensile_event_count
    model%tensile_worst_element = model%snap_tensile_worst_element
    model%tensile_worst_force = model%snap_tensile_worst_force
    model%tensile_worst_threshold = model%snap_tensile_worst_threshold
    model%tensile_worst_xi = model%snap_tensile_worst_xi
    model%tensile_worst_time = model%snap_tensile_worst_time
    model%max_step_rotation = model%snap_max_step_rotation
    model%recovery_event_count = model%snap_recovery_event_count
    model%recovery_max_substeps_used = model%snap_recovery_max_substeps_used
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
  END SUBROUTINE CD_HermiteCable_Dyn_Restore

  SUBROUTINE CD_HermiteCable_Dyn_End(model)
    !! Release the model's persistent workspaces.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    IF (ALLOCATED(model%att_node)) DEALLOCATE (model%att_node, model%att_w, model%att_b, model%att_fk, &
                                               model%att_cda, model%att_cdax)
    model%has_attach = .FALSE.
    IF (ALLOCATED(model%l0)) DEALLOCATE (model%l0)
    IF (ALLOCATED(model%EA)) DEALLOCATE (model%EA)
    IF (ALLOCATED(model%EI)) DEALLOCATE (model%EI)
    IF (ALLOCATED(model%rho_a)) DEALLOCATE (model%rho_a)
    IF (ALLOCATED(model%w)) DEALLOCATE (model%w)
    IF (ALLOCATED(model%BA)) DEALLOCATE (model%BA)
    model%has_axial_damping = .FALSE.
    IF (ALLOCATED(model%trib)) DEALLOCATE (model%trib)
    IF (ALLOCATED(model%contact_kn)) DEALLOCATE (model%contact_kn)
    IF (ALLOCATED(model%contact_cn)) DEALLOCATE (model%contact_cn)
    IF (ALLOCATED(model%fr_anchor)) DEALLOCATE (model%fr_anchor, model%fr_k, model%ws_recovery_fr, &
                                                model%snap_fr_anchor)
    model%fr_active = .FALSE.
    CALL CD_End_Bathymetry(model%contact_bathymetry)
    model%has_contact = .FALSE.; model%has_contact_bathymetry = .FALSE.
    IF (ALLOCATED(model%Mgb)) DEALLOCATE (model%Mgb)
    IF (ALLOCATED(model%ws_qn)) THEN
      DEALLOCATE (model%ws_qn, model%ws_vn, model%ws_an, model%ws_qk, model%ws_ak, &
                  model%ws_qa, model%ws_aa, model%ws_va, model%ws_vk, model%ws_R, &
                  model%ws_G, model%ws_Gt, model%ws_inert, model%ws_qtrial, model%ws_dq, &
                  model%ws_recovery_q0, model%ws_recovery_v0, model%ws_recovery_a0)
      DEALLOCATE (model%ws_Rk, model%ws_Rt, model%fc_R, model%fc_q, model%fc_v)
      DEALLOCATE (model%ws_Kb, model%ws_Kvb, &
                  model%ws_Keffb, model%ws_Msumb, model%ws_Mamb)
      DEALLOCATE (model%ws_ipiv)
      DEALLOCATE (model%ws_elem_f, model%ws_elem_K, model%ws_elem_Kdrag, model%ws_elem_Kv, &
                  model%ws_elem_Kwave, model%ws_elem_Kheld, &
                  model%ws_elem_es, model%ws_elem_em)
    END IF
    IF (ALLOCATED(model%freemask)) DEALLOCATE (model%freemask)
    IF (ALLOCATED(model%solve_mask)) DEALLOCATE (model%solve_mask)
    IF (ALLOCATED(model%fmap)) DEALLOCATE (model%fmap)
    IF (ALLOCATED(model%q)) DEALLOCATE (model%q)
    IF (ALLOCATED(model%v)) DEALLOCATE (model%v)
    IF (ALLOCATED(model%a)) DEALLOCATE (model%a)
    IF (ALLOCATED(model%diam)) DEALLOCATE (model%diam)
    IF (ALLOCATED(model%cdn)) DEALLOCATE (model%cdn)
    IF (ALLOCATED(model%cdt)) DEALLOCATE (model%cdt)
    model%has_drag = .FALSE.
    IF (ALLOCATED(model%am_diam)) DEALLOCATE (model%am_diam)
    IF (ALLOCATED(model%am_can)) DEALLOCATE (model%am_can)
    IF (ALLOCATED(model%am_cat)) DEALLOCATE (model%am_cat)
    model%has_am = .FALSE.
    model%has_wave = .FALSE.
    model%wv_h = CD_ZERO; model%wv_om = CD_ZERO; model%wv_k = CD_ZERO
    model%wv_depth = CD_ZERO; model%wv_dir = CD_ZERO
    model%wv_ncomp = 0
    IF (ALLOCATED(model%wv_ampc)) DEALLOCATE (model%wv_ampc)
    IF (ALLOCATED(model%wv_omc)) DEALLOCATE (model%wv_omc)
    IF (ALLOCATED(model%wv_kc)) DEALLOCATE (model%wv_kc)
    IF (ALLOCATED(model%wv_e2kd)) DEALLOCATE (model%wv_e2kd)
    IF (ALLOCATED(model%wv_phc)) DEALLOCATE (model%wv_phc)
    IF (ALLOCATED(model%hf_u)) DEALLOCATE (model%hf_u, model%hf_ud, model%hf_wl)
    IF (ALLOCATED(model%ws_hf_a_old)) &
      DEALLOCATE (model%ws_hf_a_old, model%ws_hf_u_old, model%ws_hf_ud_old, model%ws_hf_wl_old)
    model%has_held_field = .FALSE.
    model%has_endconn = .FALSE.
    model%endconn_k = CD_ZERO
    model%endconn_d0 = CD_ZERO
    IF (ALLOCATED(model%snap_q)) DEALLOCATE (model%snap_q, model%snap_v, model%snap_a)
    model%fc_valid = .FALSE.
    model%tr_valid = .FALSE.
    model%snap_t = CD_ZERO
    model%snap_endconn_d0 = CD_ZERO
    model%torsion = CD_HermiteTorsionType()
    model%snap_torsion_theta = CD_ZERO
    model%snap_valid = .FALSE.
    model%require_tensile = .FALSE.
    model%tensile_mode = CD_HCDYN_TENSILE_OFF
    model%tensile_strain_tolerance = CD_HCDYN_AXIAL_STRAIN_TOL
    model%tensile_event_count = 0; model%tensile_worst_element = 0
    model%tensile_worst_force = CD_ZERO; model%tensile_worst_threshold = CD_ZERO
    model%tensile_worst_xi = CD_ZERO; model%tensile_worst_time = CD_ZERO
    model%recovery_event_count = 0; model%recovery_max_substeps_used = 0
    model%t = CD_ZERO
    model%ne = 0; model%nn = 0; model%ndof = 0
    model%initialized = .FALSE.
  END SUBROUTINE CD_HermiteCable_Dyn_End

  ! ------------------------------------------------------------------------------------------

  SUBROUTINE assemble_residual(model, q_cfg, v_cfg, t_eval, want_residual, want_tangent, ErrStat, ErrMsg, &
                               endconn_direction_eval)
    !! Residual R = f_int(q) - f_ext(q, v, t), banded position tangent Kb = dR/dq, and banded
    !! velocity tangent Kvb = dR/dv at configuration q_cfg with nodal velocity v_cfg at time
    !! t_eval (LAPACK band layout, see the module note), written into the MODEL'S OWN
    !! buffers -- ws_R when want_residual, ws_Kb/ws_Kvb when want_tangent. The two flags
    !! split the Newton's work honestly: residual-only for every evaluation (predictor and
    !! line-search trials, which never solve with a tangent) and tangent-only once per
    !! solve iteration at the accepted iterate -- the tangents land directly in the
    !! committed band buffers, so no trial tangent set and no per-iteration band copies
    !! exist. Buffer access goes through pointers into the single model dummy so no model
    !! subobject is argument-associated alongside the model itself (the Fortran
    !! aliasing/intent rule). f_ext = consistent
    !! self-weight (configuration-independent) + the buoyancy of any part above the drag
    !! configuration's free surface + penalty seabed reaction + optional Morison drag (the
    !! only velocity-dependent term, hence Kvb) + optional wave loads (drag against the wave-plus-
    !! current fluid velocity, and the Froude-Krylov + fluid-inertia load from the wave acceleration)
    !! evaluated at t_eval. The caller folds Kvb into the effective tangent through the Newmark
    !! velocity-displacement coupling dv/dq. q_cfg/v_cfg may be model workspace vectors;
    !! they are only read, never defined, so that association is legal.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT), TARGET :: model
    REAL(wp), INTENT(IN) :: q_cfg(:), v_cfg(:), t_eval
    LOGICAL, INTENT(IN) :: want_residual, want_tangent
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: endconn_direction_eval(3, 2)
    REAL(wp), POINTER :: R(:), Kb(:, :), Kvb(:, :)
    INTEGER :: ne, e, a, gi, es, omp_threads
    LOGICAL :: tangent_only
    REAL(wp) :: qe(12), ve(12), Ee, fe(12), Ke(12, 12), pen
    REAL(wp) :: fdrag(12), djq(12, 12), djv(12, 12), ffk(12), fjq(12, 12)
    REAL(wp) :: fdamp(12), kdq(12, 12), kdv(12, 12)
    REAL(wp) :: hf_mid(3)
    REAL(wp), PARAMETER :: PI_L = 3.14159265358979323846_wp
    REAL(wp) :: fdry(4), kdry(4, 4), dry_wl, dry_b, dry_r
    INTEGER :: zdof(4), b
    LOGICAL :: dry
    REAL(wp) :: endconn_d0_eval(3, 2)
    CHARACTER(200) :: em2

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    ne = model%ne
    omp_threads = 1
!$  omp_threads = MIN(4, omp_get_max_threads())
!$  IF (omp_in_parallel()) omp_threads = 1
    tangent_only = want_tangent .AND. .NOT. want_residual
    R => model%ws_R
    Kb => model%ws_Kb
    Kvb => model%ws_Kvb
    endconn_d0_eval = model%endconn_d0
    IF (PRESENT(endconn_direction_eval)) endconn_d0_eval = endconn_direction_eval
    IF (want_residual) R = CD_ZERO
    IF (want_tangent) THEN
      Kb = CD_ZERO; Kvb = CD_ZERO
    END IF

    ! Tangent-only element kernels are independent. Evaluate them in parallel into persistent
    ! buffers, then validate and scatter serially in the original element order. Residual and
    ! combined-output requests retain the lower-overhead serial path. A four-thread cap avoids
    ! tiny-kernel oversubscription while respecting a smaller host OpenMP thread limit.
    IF (tangent_only) THEN
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(qe) SCHEDULE(STATIC) NUM_THREADS(omp_threads) IF(ne >= 32)
      DO e = 1, ne
        CALL evaluate_tangent_element(model, e, q_cfg, v_cfg, t_eval, model%ws_elem_K(:, :, e), &
                                      model%ws_elem_Kdrag(:, :, e), model%ws_elem_Kv(:, :, e), &
                                      model%ws_elem_Kwave(:, :, e), model%ws_elem_Kheld(:, :, e), &
                                      model%ws_elem_es(e), model%ws_elem_em(e))
      END DO
      !$OMP END PARALLEL DO
      DO e = 1, ne
        IF (model%ws_elem_es(e) /= CD_HCABLE_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'element: '//TRIM(model%ws_elem_em(e)); RETURN
        END IF
        CALL scatter12_band(Kb, model%ws_elem_K(:, :, e), e)
      END DO
    ELSE
      DO e = 1, ne
        qe(1:3) = q_cfg(6*(e - 1) + 1:6*(e - 1) + 3)
        qe(4:6) = q_cfg(6*(e - 1) + 4:6*(e - 1) + 6)
        qe(7:9) = q_cfg(6*e + 1:6*e + 3)
        qe(10:12) = q_cfg(6*e + 4:6*e + 6)
        IF (want_tangent) THEN
          CALL CD_HermiteCable_Element(qe, model%l0(e), model%EA(e), model%EI(e), Ee, fe, Ke, es, em2, &
                                       symmetric_half=.TRUE., &
                                       axial_quadrature_order=model%axial_quadrature_order, &
                                       bending_quadrature_order=model%bending_quadrature_order)
        ELSE
          CALL CD_HermiteCable_Element(qe, model%l0(e), model%EA(e), model%EI(e), Ee, fe, ErrStat=es, ErrMsg=em2, &
                                       axial_quadrature_order=model%axial_quadrature_order, &
                                       bending_quadrature_order=model%bending_quadrature_order)
        END IF
        IF (es /= CD_HCABLE_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'element: '//TRIM(em2); RETURN
        END IF
        IF (want_residual) THEN
          DO a = 1, 12
            gi = elem_gdof(e, a)
            R(gi) = R(gi) + fe(a)
          END DO
        END IF
        IF (want_tangent) CALL scatter12_band(Kb, Ke, e)
      END DO
    END IF

    IF (want_residual) THEN
      ! consistent cubic-Hermite self-weight on z DOFs (R = fint - fext => -fext adds +w)
      DO e = 1, ne
        R(6*(e - 1) + 3) = R(6*(e - 1) + 3) + model%w(e)*0.5_wp*model%l0(e)
        R(6*(e - 1) + 6) = R(6*(e - 1) + 6) + model%w(e)*model%l0(e)*model%l0(e)/12.0_wp
        R(6*e + 3) = R(6*e + 3) + model%w(e)*0.5_wp*model%l0(e)
        R(6*e + 6) = R(6*e + 6) - model%w(e)*model%l0(e)*model%l0(e)/12.0_wp
      END DO
    END IF
    IF (model%has_attach) THEN
      CALL assemble_attachments(model, q_cfg, v_cfg, want_residual, want_tangent, R, Kb, Kvb, es, em2)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = TRIM(em2); RETURN
      END IF
    END IF

    ! Dry-part buoyancy recovery. w is the submerged weight; the part of an element above the
    ! free surface of the drag configuration regains its displaced-water buoyancy. The surface
    ! is the still waterline, or with a held host field the element's sampled elevation (the
    ! same surface the drag uses). Fully submerged elements return immediately. Without a
    ! drag configuration the added-mass configuration defines the surface and section, so
    ! the buoyancy is restored whenever either hydrodynamic configuration is present.
    IF (model%has_drag .OR. model%has_am) THEN
      DO e = 1, ne
        IF (model%has_drag) THEN
          dry_wl = model%waterline_z
          dry_b = model%rho_w*model%gravity*0.25_wp*PI_L*model%diam(e)**2
          dry_r = 0.5_wp*model%diam(e)
        ELSE
          dry_wl = model%am_wl
          dry_b = model%am_rho*model%gravity*0.25_wp*PI_L*model%am_diam(e)**2
          dry_r = 0.5_wp*model%am_diam(e)
        END IF
        IF (model%has_held_field) dry_wl = 0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1))
        zdof = [6*(e - 1) + 3, 6*(e - 1) + 6, 6*e + 3, 6*e + 6]
        CALL CD_HermiteCable_Dry_Buoyancy(q_cfg(zdof), model%l0(e), dry_b, dry_wl, fdry, kdry, dry, &
                                          radius=dry_r)
        IF (.NOT. dry) CYCLE
        IF (want_residual) R(zdof) = R(zdof) + fdry
        IF (want_tangent) THEN
          DO b = 1, 4
            DO a = 1, 4
              Kb(KL_D + KU_D + 1 + zdof(a) - zdof(b), zdof(b)) = &
                Kb(KL_D + KU_D + 1 + zdof(a) - zdof(b), zdof(b)) + kdry(a, b)
            END DO
          END DO
        END IF
      END DO
    END IF

    ! Standalone production contact uses the same nodal C1 normal/damper/friction
    ! law as the deck compatibility path, including structured-bed slope terms.
    ! Retain the scalar flat penalty for existing direct API callers.
    IF (model%has_contact) THEN
      CALL assemble_contact(model, q_cfg, v_cfg, want_residual, want_tangent, es, em2)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = 'contact: '//TRIM(em2); RETURN
      END IF
    ELSE IF (model%kn > CD_ZERO) THEN
      DO a = 1, model%nn
        gi = 6*(a - 1) + 3
        pen = model%seabed_z - q_cfg(gi)
        IF (pen > CD_ZERO) THEN
          IF (want_residual) R(gi) = R(gi) - model%kn*pen*model%trib(a)
          IF (want_tangent) Kb(KL_D + KU_D + 1, gi) = Kb(KL_D + KU_D + 1, gi) + model%kn*model%trib(a)
        END IF
      END DO
    END IF

    ! Scatter the element-private hydrodynamic tangents in the same phase and element
    ! order as the serial formulation, retaining its floating-point accumulation order.
    IF (tangent_only .AND. (model%has_drag .OR. model%has_axial_damping)) THEN
      DO e = 1, ne
        CALL scatter12_band(Kb, model%ws_elem_Kdrag(:, :, e), e)
        CALL scatter12_band(Kvb, model%ws_elem_Kv(:, :, e), e)
      END DO
    END IF
    IF (tangent_only .AND. model%has_wave .AND. model%has_am) THEN
      DO e = 1, ne
        CALL scatter12_band(Kb, model%ws_elem_Kwave(:, :, e), e)
      END DO
    END IF
    IF (tangent_only .AND. model%has_held_field .AND. model%has_am) THEN
      DO e = 1, ne
        CALL scatter12_band(Kb, model%ws_elem_Kheld(:, :, e), e)
      END DO
    END IF

    ! Axial Kelvin-Voigt damping is an internal constitutive force. It uses the same
    ! reference-arc Hermite interpolation and axial quadrature as the elastic resultant.
    IF (model%has_axial_damping .AND. .NOT. tangent_only) THEN
      DO e = 1, ne
        qe(1:3) = q_cfg(6*(e - 1) + 1:6*(e - 1) + 3); qe(4:6) = q_cfg(6*(e - 1) + 4:6*(e - 1) + 6)
        qe(7:9) = q_cfg(6*e + 1:6*e + 3); qe(10:12) = q_cfg(6*e + 4:6*e + 6)
        ve(1:3) = v_cfg(6*(e - 1) + 1:6*(e - 1) + 3); ve(4:6) = v_cfg(6*(e - 1) + 4:6*(e - 1) + 6)
        ve(7:9) = v_cfg(6*e + 1:6*e + 3); ve(10:12) = v_cfg(6*e + 4:6*e + 6)
        IF (want_tangent) THEN
          CALL CD_HermiteCable_Axial_Damping_Element(qe, ve, model%l0(e), model%BA(e), &
                                                     fdamp, kdq, kdv, es, em2, &
                                                     model%axial_quadrature_order)
        ELSE
          CALL CD_HermiteCable_Axial_Damping_Element(qe, ve, model%l0(e), model%BA(e), &
                                                     fdamp, ErrStat=es, ErrMsg=em2, &
                                                     quadrature_order=model%axial_quadrature_order)
        END IF
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = TRIM(em2); RETURN
        END IF
        IF (want_residual) THEN
          DO a = 1, 12
            gi = elem_gdof(e, a)
            R(gi) = R(gi) + fdamp(a)
          END DO
        END IF
        IF (want_tangent) THEN
          CALL scatter12_band(Kb, kdq, e)
          CALL scatter12_band(Kvb, kdv, e)
        END IF
      END DO
    END IF

    ! Consistent Morison quadratic drag on the submerged cable, distributed over the full cubic-
    ! Hermite shape functions (so the material-tangent DOFs, which move the element interior, also
    ! carry drag). The element returns the 12-DOF generalized drag LOAD fdrag (f_ext) and its force
    ! Jacobians djq = d(fdrag)/dq, djv = d(fdrag)/dv. R = f_int - f_ext, so R -= fdrag, K -= djq,
    ! Kv -= djv over all 12 element DOFs.
    IF (model%has_drag .AND. .NOT. tangent_only) THEN
      DO e = 1, ne
        qe(1:3) = q_cfg(6*(e - 1) + 1:6*(e - 1) + 3); qe(4:6) = q_cfg(6*(e - 1) + 4:6*(e - 1) + 6)
        qe(7:9) = q_cfg(6*e + 1:6*e + 3); qe(10:12) = q_cfg(6*e + 4:6*e + 6)
        ve(1:3) = v_cfg(6*(e - 1) + 1:6*(e - 1) + 3); ve(4:6) = v_cfg(6*(e - 1) + 4:6*(e - 1) + 6)
        ve(7:9) = v_cfg(6*e + 1:6*e + 3); ve(10:12) = v_cfg(6*e + 4:6*e + 6)
        IF (want_tangent) THEN
          IF (model%has_wave .AND. model%wv_ncomp > 0) THEN
            CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                              model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                              fdrag, djq, djv, es, em2, &
                                              wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=t_eval, &
                                              wvc_amp=model%wv_ampc, wvc_om=model%wv_omc, &
                                              wvc_k=model%wv_kc, wvc_ph=model%wv_phc, &
                                              wvc_e2kd=model%wv_e2kd)
          ELSE IF (model%has_wave) THEN
            CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                              model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                              fdrag, djq, djv, es, em2, &
                                              wv_h=model%wv_h, wv_om=model%wv_om, wv_k=model%wv_k, &
                                              wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=t_eval)
          ELSE IF (model%has_held_field) THEN
            ! HELD host field: the element's ambient velocity = the drag config's constant
            ! current + the two end nodes' held-velocity average; the element waterline =
            ! the two nodes' sampled-elevation average (the host's spatial granularity).
            CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), &
                                              model%current + 0.5_wp*(model%hf_u(:, e) + model%hf_u(:, e + 1)), &
                                              0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1)), &
                                              model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                              fdrag, djq, djv, es, em2)
          ELSE
            CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                              model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                              fdrag, djq, djv, es, em2)
          END IF
        ELSE
          IF (model%has_wave .AND. model%wv_ncomp > 0) THEN
            CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                              model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                              fdrag, ErrStat=es, ErrMsg=em2, &
                                              wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=t_eval, &
                                              wvc_amp=model%wv_ampc, wvc_om=model%wv_omc, &
                                              wvc_k=model%wv_kc, wvc_ph=model%wv_phc, &
                                              wvc_e2kd=model%wv_e2kd)
          ELSE IF (model%has_wave) THEN
            CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                              model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                              fdrag, ErrStat=es, ErrMsg=em2, &
                                              wv_h=model%wv_h, wv_om=model%wv_om, wv_k=model%wv_k, &
                                              wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=t_eval)
          ELSE IF (model%has_held_field) THEN
            CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), &
                                              model%current + 0.5_wp*(model%hf_u(:, e) + model%hf_u(:, e + 1)), &
                                              0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1)), &
                                              model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                              fdrag, ErrStat=es, ErrMsg=em2)
          ELSE
            CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                              model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                              fdrag, ErrStat=es, ErrMsg=em2)
          END IF
        END IF
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = es; ErrMsg = TRIM(em2); RETURN
        END IF
        IF (want_residual) THEN
          DO a = 1, 12
            gi = elem_gdof(e, a)
            R(gi) = R(gi) - fdrag(a)
          END DO
        END IF
        IF (want_tangent) THEN
          CALL scatter12_band(Kb, -djq, e)
          CALL scatter12_band(Kvb, -djv, e)
        END IF
      END DO
    END IF

    ! ---- rotational line-end connections ----
    ! A node-local generalized load on the end tangents: the 3x3 tangent block sits on the
    ! band diagonal (|gi - gj| <= 2 against KU_D = 11), so it cannot widen the band. The
    ! spring is a potential, so it enters the residual with the same sign as the internal
    ! force, and f . d == 0 keeps it off the stretch component by construction.
    IF (model%has_endconn) THEN
      CALL add_endconn_band(R, Kb, q_cfg, 1, model%endconn_k(1), endconn_d0_eval(:, 1), &
                            want_residual, want_tangent, es, em2)
      IF (es /= CD_ENDCONN_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE
        ErrMsg = 'end connection at node 1: '//TRIM(em2)
        RETURN
      END IF
      CALL add_endconn_band(R, Kb, q_cfg, ne + 1, model%endconn_k(2), endconn_d0_eval(:, 2), &
                            want_residual, want_tangent, es, em2)
      IF (es /= CD_ENDCONN_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE
        ErrMsg = 'end connection at final node: '//TRIM(em2)
        RETURN
      END IF
    END IF

    ! Froude-Krylov + fluid-inertia wave load (from the wave ACCELERATION), with the added-mass
    ! configuration's coefficients. An external load: R -= ffk, K -= fjq (geometric tangent chain;
    ! the wave-field spatial gradient is neglected in K -- see the element note).
    IF (model%has_wave .AND. model%has_am .AND. .NOT. tangent_only) THEN
      DO e = 1, ne
        qe(1:3) = q_cfg(6*(e - 1) + 1:6*(e - 1) + 3); qe(4:6) = q_cfg(6*(e - 1) + 4:6*(e - 1) + 6)
        qe(7:9) = q_cfg(6*e + 1:6*e + 3); qe(10:12) = q_cfg(6*e + 4:6*e + 6)
        IF (model%wv_ncomp > 0) THEN
          IF (want_tangent) THEN
            CALL CD_HermiteCable_FK_Element(qe, model%l0(e), model%am_wl, model%am_rho, model%am_diam(e), &
                                            model%am_can(e), model%am_cat(e), &
                                            wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=t_eval, &
                                            ffk=ffk, fjq=fjq, ErrStat=es, ErrMsg=em2, &
                                            wvc_amp=model%wv_ampc, wvc_om=model%wv_omc, &
                                            wvc_k=model%wv_kc, wvc_ph=model%wv_phc, &
                                            wvc_e2kd=model%wv_e2kd)
          ELSE
            CALL CD_HermiteCable_FK_Element(qe, model%l0(e), model%am_wl, model%am_rho, model%am_diam(e), &
                                            model%am_can(e), model%am_cat(e), &
                                            wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=t_eval, &
                                            ffk=ffk, ErrStat=es, ErrMsg=em2, &
                                            wvc_amp=model%wv_ampc, wvc_om=model%wv_omc, &
                                            wvc_k=model%wv_kc, wvc_ph=model%wv_phc, &
                                            wvc_e2kd=model%wv_e2kd)
          END IF
        ELSE
          IF (want_tangent) THEN
            CALL CD_HermiteCable_FK_Element(qe, model%l0(e), model%am_wl, model%am_rho, model%am_diam(e), &
                                            model%am_can(e), model%am_cat(e), &
                                            wv_h=model%wv_h, wv_om=model%wv_om, wv_k=model%wv_k, &
                                            wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=t_eval, &
                                            ffk=ffk, fjq=fjq, ErrStat=es, ErrMsg=em2)
          ELSE
            CALL CD_HermiteCable_FK_Element(qe, model%l0(e), model%am_wl, model%am_rho, model%am_diam(e), &
                                            model%am_can(e), model%am_cat(e), &
                                            wv_h=model%wv_h, wv_om=model%wv_om, wv_k=model%wv_k, &
                                            wv_depth=model%wv_depth, wv_dir=model%wv_dir, wv_t=t_eval, &
                                            ffk=ffk, ErrStat=es, ErrMsg=em2)
          END IF
        END IF
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = es; ErrMsg = TRIM(em2); RETURN
        END IF
        IF (want_residual) THEN
          DO a = 1, 12
            gi = elem_gdof(e, a)
            R(gi) = R(gi) - ffk(a)
          END DO
        END IF
        IF (want_tangent) CALL scatter12_band(Kb, -fjq, e)
      END DO
    END IF

    ! HELD-FIELD Froude-Krylov + fluid-inertia load: the host-sampled fluid acceleration
    ! (per-element two-node average, frozen over the step) with the added-mass
    ! configuration's coefficients -- the held-mode twin of the wave FK block above.
    IF (model%has_held_field .AND. model%has_am .AND. .NOT. tangent_only) THEN
      DO e = 1, ne
        qe(1:3) = q_cfg(6*(e - 1) + 1:6*(e - 1) + 3); qe(4:6) = q_cfg(6*(e - 1) + 4:6*(e - 1) + 6)
        qe(7:9) = q_cfg(6*e + 1:6*e + 3); qe(10:12) = q_cfg(6*e + 4:6*e + 6)
        IF (want_tangent) THEN
          hf_mid(1) = 0.5_wp*(model%hf_ud(1, e) + model%hf_ud(1, e + 1))
          hf_mid(2) = 0.5_wp*(model%hf_ud(2, e) + model%hf_ud(2, e + 1))
          hf_mid(3) = 0.5_wp*(model%hf_ud(3, e) + model%hf_ud(3, e + 1))
          CALL CD_HermiteCable_FK_Element(qe, model%l0(e), &
                                          0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1)), &
                                          model%am_rho, model%am_diam(e), &
                                          model%am_can(e), model%am_cat(e), &
                                          wv_depth=CD_ZERO, wv_dir=CD_ZERO, wv_t=CD_ZERO, &
                                          ffk=ffk, fjq=fjq, ErrStat=es, ErrMsg=em2, &
                                          hf_a=hf_mid)
        ELSE
          hf_mid(1) = 0.5_wp*(model%hf_ud(1, e) + model%hf_ud(1, e + 1))
          hf_mid(2) = 0.5_wp*(model%hf_ud(2, e) + model%hf_ud(2, e + 1))
          hf_mid(3) = 0.5_wp*(model%hf_ud(3, e) + model%hf_ud(3, e + 1))
          CALL CD_HermiteCable_FK_Element(qe, model%l0(e), &
                                          0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1)), &
                                          model%am_rho, model%am_diam(e), &
                                          model%am_can(e), model%am_cat(e), &
                                          wv_depth=CD_ZERO, wv_dir=CD_ZERO, wv_t=CD_ZERO, &
                                          ffk=ffk, ErrStat=es, ErrMsg=em2, &
                                          hf_a=hf_mid)
        END IF
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = es; ErrMsg = TRIM(em2); RETURN
        END IF
        IF (want_residual) THEN
          DO a = 1, 12
            gi = elem_gdof(e, a)
            R(gi) = R(gi) - ffk(a)
          END DO
        END IF
        IF (want_tangent) CALL scatter12_band(Kb, -fjq, e)
      END DO
    END IF

    ! Condensed torsion (statics only in this build): R <- R - M_t dTheta/dq at q_cfg.
    IF (model%torsion%active .AND. want_residual) THEN
      BLOCK
        REAL(wp) :: th, mt, gq(model%ndof)
        CALL CD_HermiteCable_Dyn_Torsion_State(model, th, mt, es, em2, gq, q_eval=q_cfg)
        IF (es /= CD_HCDYN_OK) THEN
          ErrStat = es; ErrMsg = TRIM(em2); RETURN
        END IF
        R = R - mt*gq
      END BLOCK
    END IF
  END SUBROUTINE assemble_residual

  SUBROUTINE add_endconn_band(R, Kb, q_cfg, node, k_rot, d0, want_residual, want_tangent, &
                              ErrStat, ErrMsg)
    !! Scatter one rotational end connection into the dynamic residual and band tangent.
    REAL(wp), INTENT(INOUT) :: R(:), Kb(:, :)
    REAL(wp), INTENT(IN) :: q_cfg(:), k_rot, d0(3)
    INTEGER, INTENT(IN) :: node
    LOGICAL, INTENT(IN) :: want_residual, want_tangent
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: f_ec(3), k_ec(3, 3)
    INTEGER :: base, i, j, gi, gj
    ErrStat = CD_ENDCONN_OK
    ErrMsg = ''
    IF (k_rot <= CD_ZERO) RETURN
    base = 6*(node - 1) + 3                      ! tangent DOFs are base+1 .. base+3
    CALL CD_EndConn_Spring(q_cfg(base + 1:base + 3), d0, k_rot, f_ec, k_ec, ErrStat, ErrMsg)
    IF (ErrStat /= CD_ENDCONN_OK) RETURN
    IF (want_residual) R(base + 1:base + 3) = R(base + 1:base + 3) + f_ec
    IF (want_tangent) THEN
      DO j = 1, 3
        gj = base + j
        DO i = 1, 3
          gi = base + i
          Kb(KL_D + KU_D + 1 + gi - gj, gj) = Kb(KL_D + KU_D + 1 + gi - gj, gj) + k_ec(i, j)
        END DO
      END DO
    END IF
  END SUBROUTINE add_endconn_band

  SUBROUTINE assemble_contact(model, q_cfg, v_cfg, want_residual, want_tangent, ErrStat, ErrMsg)
    !! Assemble nodal contact directly into the model-owned residual bands. The
    !! Jacobians are d(f_int-f_ext)/dq and d(f_int-f_ext)/dv.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: q_cfg(:), v_cfg(:)
    LOGICAL, INTENT(IN) :: want_residual, want_tangent
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: i, ix, iy, iz, es, jr, jc
    REAL(wp) :: xl, yl, xg, yg, z_floor, gxg, gyg, gx, gy
    REAL(wp) :: gap, normal, normal_gap, gate, gate_gap
    REAL(wp) :: dnormal_dgap, dnormal_dvz, dnormal_dq(3)
    ! Frictionless contact on a sloped floor acts along its upward normal nvec: normal
    ! penetration gap/sfac, normal velocity vn (sfac = 1 and nvec = e_z on a level floor).
    REAL(wp) :: sfac, nvec(3), vn, dnormal_dv(3)
    REAL(wp) :: d_fr(2), mu_c
    REAL(wp) :: ffr(2), dfr_dd(2, 2), dfr_dc(2)
    CHARACTER(200) :: em

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    IF (.NOT. model%has_contact .OR. .NOT. ALLOCATED(model%contact_kn) .OR. &
        .NOT. ALLOCATED(model%contact_cn)) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'contact configuration is incomplete'; RETURN
    END IF
    DO i = 1, model%nn
      ix = 6*i - 5; iy = ix + 1; iz = ix + 2
      gx = CD_ZERO; gy = CD_ZERO
      IF (model%has_contact_bathymetry) THEN
        xl = q_cfg(ix); yl = q_cfg(iy)
        xg = model%contact_frame_c*xl - model%contact_frame_s*yl
        yg = model%contact_frame_s*xl + model%contact_frame_c*yl
        CALL CD_Bathymetry_Floor_Gradient(model%contact_bathymetry, xg, yg, z_floor, gxg, gyg, es, em)
        IF (es /= CD_BATHY_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'bathymetry query failed: '//TRIM(em); RETURN
        END IF
        ! Chain the global floor gradient through the inverse local-to-global heading.
        gx = model%contact_frame_c*gxg + model%contact_frame_s*gyg
        gy = -model%contact_frame_s*gxg + model%contact_frame_c*gyg
      ELSE
        z_floor = model%seabed_z
      END IF

      sfac = SQRT(CD_ONE + gx*gx + gy*gy)
      nvec = [-gx, -gy, CD_ONE]/sfac
      gap = (z_floor - q_cfg(iz))/sfac
      CALL CD_Seabed_Normal_Law(gap, model%contact_kn(i), normal, normal_gap)
      CALL CD_Seabed_Contact_Gate(gap, gate, gate_gap)
      IF (normal <= CD_ZERO .AND. gate <= CD_ZERO) CYCLE
      vn = DOT_PRODUCT(v_cfg(ix:iz), nvec)
      dnormal_dgap = normal_gap
      dnormal_dvz = CD_ZERO
      IF (vn < CD_ZERO) THEN
        normal = normal - model%contact_cn(i)*gate*vn
        dnormal_dgap = dnormal_dgap - model%contact_cn(i)*gate_gap*vn
        dnormal_dvz = -model%contact_cn(i)*gate
      END IF
      ! d(gap/sfac)/dq = (gx, gy, -1)/sfac = -nvec (floor gradient held over the node)
      dnormal_dq = -dnormal_dgap*nvec
      dnormal_dv = dnormal_dvz*nvec

      IF (want_residual) model%ws_R(ix:iz) = model%ws_R(ix:iz) - normal*nvec
      IF (want_tangent) THEN
        DO jc = 1, 3
          DO jr = 1, 3
            CALL scatter_contact_band(model%ws_Kb, ix + jr - 1, ix + jc - 1, -nvec(jr)*dnormal_dq(jc))
            CALL scatter_contact_band(model%ws_Kvb, ix + jr - 1, ix + jc - 1, -nvec(jr)*dnormal_dv(jc))
          END DO
        END DO
      END IF

      IF (.NOT. model%fr_active) CYCLE
      ! Stick-slip friction spring to the committed anchor (CD_Seabed_Friction_Stick_Slip):
      ! capacity mu*normal; its dependence on the normal force enters the tangent.
      d_fr(1) = q_cfg(ix) - model%fr_anchor(1, i)
      d_fr(2) = q_cfg(iy) - model%fr_anchor(2, i)
      IF (model%contact_fr_aniso) THEN
        ! dfr_dc is then d(force)/d(normal): the capacity chain factor is one
        CALL CD_Seabed_Friction_Aniso(CD_FRICTION_STICK_SLIP, d_fr, model%fr_k(i), normal, &
                                      model%contact_mu_axial, model%contact_mu, q_cfg(ix + 3:ix + 4), &
                                      ffr, dfr_dd, dfr_dc)
        mu_c = CD_ONE
      ELSE
        CALL CD_Seabed_Friction_Stick_Slip(d_fr, model%fr_k(i), &
                                           model%contact_mu*normal, ffr, dfr_dd, dfr_dc)
        mu_c = model%contact_mu
      END IF
      IF (want_residual) model%ws_R(ix:iy) = model%ws_R(ix:iy) + ffr
      IF (want_tangent) THEN
        DO jr = 1, 2
          CALL scatter_contact_band(model%ws_Kb, ix + jr - 1, ix, dfr_dd(jr, 1) + &
                                    mu_c*dfr_dc(jr)*dnormal_dq(1))
          CALL scatter_contact_band(model%ws_Kb, ix + jr - 1, iy, dfr_dd(jr, 2) + &
                                    mu_c*dfr_dc(jr)*dnormal_dq(2))
          CALL scatter_contact_band(model%ws_Kb, ix + jr - 1, iz, mu_c*dfr_dc(jr)*dnormal_dq(3))
          DO jc = 1, 3
            CALL scatter_contact_band(model%ws_Kvb, ix + jr - 1, ix + jc - 1, mu_c*dfr_dc(jr)*dnormal_dv(jc))
          END DO
        END DO
      END IF
    END DO
  END SUBROUTINE assemble_contact

  SUBROUTINE commit_friction_anchors(model)
    !! Return mapping of the elastic-plastic friction springs at the committed state: a
    !! spring stretched beyond the capacity mu*normal moves its anchor to the capacity
    !! distance behind the node (a node off the seabed, capacity zero, carries its anchor
    !! with it), so the committed force is exactly the one of the converged residual.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER :: i, ix, iz, es
    REAL(wp) :: xl, yl, xg, yg, z_floor, gxg, gyg, gap, normal, normal_gap, gate, gate_gap
    REAL(wp) :: dxy(2), dlen, cap, mu_d, q_d, gq_d(2), gx, gy, sfac, nvec(3), vn
    CHARACTER(200) :: em
    DO i = 1, model%nn
      ix = 6*i - 5; iz = ix + 2
      gx = CD_ZERO
      gy = CD_ZERO
      IF (model%has_contact_bathymetry) THEN
        xl = model%q(ix); yl = model%q(ix + 1)
        xg = model%contact_frame_c*xl - model%contact_frame_s*yl
        yg = model%contact_frame_s*xl + model%contact_frame_c*yl
        CALL CD_Bathymetry_Floor_Gradient(model%contact_bathymetry, xg, yg, z_floor, gxg, gyg, es, em)
        IF (es /= CD_BATHY_OK) THEN
          ! Not reached in practice: the converged step's contact assembly queried this
          ! committed position and fails closed on the same error.
          z_floor = model%seabed_z
        ELSE
          gx = model%contact_frame_c*gxg + model%contact_frame_s*gyg
          gy = -model%contact_frame_s*gxg + model%contact_frame_c*gyg
        END IF
      ELSE
        z_floor = model%seabed_z
      END IF
      ! the normal force of assemble_contact (along the floor normal)
      sfac = SQRT(CD_ONE + gx*gx + gy*gy)
      nvec = [-gx, -gy, CD_ONE]/sfac
      gap = (z_floor - model%q(iz))/sfac
      CALL CD_Seabed_Normal_Law(gap, model%contact_kn(i), normal, normal_gap)
      CALL CD_Seabed_Contact_Gate(gap, gate, gate_gap)
      vn = DOT_PRODUCT(model%v(ix:iz), nvec)
      IF (vn < CD_ZERO) normal = normal - model%contact_cn(i)*gate*vn
      dxy = model%q(ix:ix + 1) - model%fr_anchor(:, i)
      IF (model%contact_fr_aniso) THEN
        CALL CD_Seabed_Friction_Mu_Dir(dxy, model%contact_mu_axial, model%contact_mu, model%q(ix + 3:ix + 4), &
                                       mu_d, q_d, gq_d)
        cap = mu_d*MAX(normal, CD_ZERO)
      ELSE
        cap = model%contact_mu*MAX(normal, CD_ZERO)
      END IF
      dlen = SQRT(dxy(1)*dxy(1) + dxy(2)*dxy(2))
      IF (model%fr_k(i)*dlen > cap) THEN
        model%fr_anchor(:, i) = model%q(ix:ix + 1) - (cap/(model%fr_k(i)*dlen))*dxy
      END IF
    END DO
  END SUBROUTINE commit_friction_anchors

  SUBROUTINE scatter_contact_band(ab, row_dof, col_dof, value)
    REAL(wp), INTENT(INOUT) :: ab(:, :)
    INTEGER, INTENT(IN) :: row_dof, col_dof
    REAL(wp), INTENT(IN) :: value
    INTEGER :: brow
    brow = KL_D + KU_D + 1 + row_dof - col_dof
    IF (brow >= 1 .AND. brow <= SIZE(ab, 1)) ab(brow, col_dof) = ab(brow, col_dof) + value
  END SUBROUTINE scatter_contact_band

  SUBROUTINE evaluate_tangent_element(model, e, q_cfg, v_cfg, t_eval, Kstruct, Kdrag, Kvel, Kwave, Kheld, &
                                      ErrStat, ErrMsg)
    !! Evaluate the complete element-local Newton package. Keeping structural, drag and
    !! fluid-inertia tangents in one worker avoids three OpenMP regions per iteration;
    !! the caller scatters these private blocks in deterministic element order.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: e
    REAL(wp), INTENT(IN) :: q_cfg(:), v_cfg(:), t_eval
    REAL(wp), INTENT(OUT) :: Kstruct(12, 12), Kdrag(12, 12), Kvel(12, 12)
    REAL(wp), INTENT(OUT) :: Kwave(12, 12), Kheld(12, 12)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp) :: qe(12), ve(12), energy, force(12), fload(12)
    REAL(wp) :: djq(12, 12), djv(12, 12), fjq(12, 12)
    REAL(wp) :: hf_mid(3)
    INTEGER :: es
    CHARACTER(200) :: em

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    qe(1:6) = q_cfg(6*(e - 1) + 1:6*e)
    qe(7:12) = q_cfg(6*e + 1:6*(e + 1))
    ve(1:6) = v_cfg(6*(e - 1) + 1:6*e)
    ve(7:12) = v_cfg(6*e + 1:6*(e + 1))
    Kdrag = CD_ZERO; Kvel = CD_ZERO; Kwave = CD_ZERO; Kheld = CD_ZERO

    CALL CD_HermiteCable_Element(qe, model%l0(e), model%EA(e), model%EI(e), energy, force, Kstruct, es, em, &
                                 symmetric_half=.TRUE., tangent_only=.TRUE., &
                                 axial_quadrature_order=model%axial_quadrature_order, &
                                 bending_quadrature_order=model%bending_quadrature_order)
    IF (es /= CD_HCABLE_OK) THEN
      ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'element: '//TRIM(em); RETURN
    END IF

    IF (model%has_axial_damping) THEN
      CALL CD_HermiteCable_Axial_Damping_Element(qe, ve, model%l0(e), model%BA(e), &
                                                 fload, djq, djv, es, em, &
                                                 model%axial_quadrature_order)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = TRIM(em); RETURN
      END IF
      Kdrag = Kdrag + djq
      Kvel = Kvel + djv
    END IF

    IF (model%has_drag) THEN
      IF (model%has_wave .AND. model%wv_ncomp > 0) THEN
        CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                          model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                          fload, djq, djv, es, em, wv_depth=model%wv_depth, &
                                          wv_dir=model%wv_dir, wv_t=t_eval, wvc_amp=model%wv_ampc, &
                                          wvc_om=model%wv_omc, wvc_k=model%wv_kc, wvc_ph=model%wv_phc, &
                                          wvc_e2kd=model%wv_e2kd)
      ELSE IF (model%has_wave) THEN
        CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                          model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                          fload, djq, djv, es, em, wv_h=model%wv_h, wv_om=model%wv_om, &
                                          wv_k=model%wv_k, wv_depth=model%wv_depth, wv_dir=model%wv_dir, &
                                          wv_t=t_eval)
      ELSE IF (model%has_held_field) THEN
        CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), &
                                          model%current + 0.5_wp*(model%hf_u(:, e) + model%hf_u(:, e + 1)), &
                                          0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1)), model%rho_w, &
                                          model%diam(e), model%cdn(e), model%cdt(e), fload, djq, djv, es, em)
      ELSE
        CALL CD_HermiteCable_Drag_Element(qe, ve, model%l0(e), model%current, model%waterline_z, &
                                          model%rho_w, model%diam(e), model%cdn(e), model%cdt(e), &
                                          fload, djq, djv, es, em)
      END IF
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = TRIM(em); RETURN
      END IF
      Kdrag = Kdrag - djq
      Kvel = Kvel - djv
    END IF

    IF (model%has_wave .AND. model%has_am) THEN
      IF (model%wv_ncomp > 0) THEN
        CALL CD_HermiteCable_FK_Element(qe, model%l0(e), model%am_wl, model%am_rho, model%am_diam(e), &
                                        model%am_can(e), model%am_cat(e), wv_depth=model%wv_depth, &
                                        wv_dir=model%wv_dir, wv_t=t_eval, ffk=fload, fjq=fjq, &
                                        ErrStat=es, ErrMsg=em, wvc_amp=model%wv_ampc, wvc_om=model%wv_omc, &
                                        wvc_k=model%wv_kc, wvc_ph=model%wv_phc, &
                                        wvc_e2kd=model%wv_e2kd)
      ELSE
        CALL CD_HermiteCable_FK_Element(qe, model%l0(e), model%am_wl, model%am_rho, model%am_diam(e), &
                                        model%am_can(e), model%am_cat(e), wv_h=model%wv_h, &
                                        wv_om=model%wv_om, wv_k=model%wv_k, wv_depth=model%wv_depth, &
                                        wv_dir=model%wv_dir, wv_t=t_eval, ffk=fload, fjq=fjq, &
                                        ErrStat=es, ErrMsg=em)
      END IF
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = TRIM(em); RETURN
      END IF
      Kwave = -fjq
    END IF

    IF (model%has_held_field .AND. model%has_am) THEN
      hf_mid(1) = 0.5_wp*(model%hf_ud(1, e) + model%hf_ud(1, e + 1))
      hf_mid(2) = 0.5_wp*(model%hf_ud(2, e) + model%hf_ud(2, e + 1))
      hf_mid(3) = 0.5_wp*(model%hf_ud(3, e) + model%hf_ud(3, e + 1))
      CALL CD_HermiteCable_FK_Element(qe, model%l0(e), 0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1)), &
                                      model%am_rho, model%am_diam(e), model%am_can(e), model%am_cat(e), &
                                      wv_depth=CD_ZERO, wv_dir=CD_ZERO, wv_t=CD_ZERO, ffk=fload, fjq=fjq, &
                                      ErrStat=es, ErrMsg=em, &
                                      hf_a=hf_mid)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = TRIM(em); RETURN
      END IF
      Kheld = -fjq
    END IF
  END SUBROUTINE evaluate_tangent_element

  SUBROUTINE CD_HermiteCable_Axial_Damping_Resultant(qe, ve, l0, BA, xi, resultant, ErrStat, ErrMsg)
    !! Signed viscous axial resultant N_v = BA*strain_rate at one material point.
    !! This is the damping share added to the elastic resultant for tension outputs.
    REAL(wp), INTENT(IN) :: qe(12), ve(12), l0, BA, xi
    REAL(wp), INTENT(OUT) :: resultant
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: H(4), dH(4), r_xi(3), v_xi(3), speed
    INTEGER :: s

    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    resultant = CD_ZERO
    IF (.NOT. CD_Is_Finite(l0) .OR. l0 <= CD_ZERO .OR. &
        .NOT. CD_Is_Finite(BA) .OR. BA < CD_ZERO .OR. &
        .NOT. CD_Is_Finite(xi) .OR. xi < CD_ZERO .OR. xi > CD_ONE) THEN
      CALL fail('need finite l0 > 0, BA >= 0, and xi in [0,1]'); RETURN
    END IF
    IF (.NOT. CD_All_Finite(qe) .OR. .NOT. CD_All_Finite(ve)) THEN
      CALL fail('qe and ve must be finite'); RETURN
    END IF
    IF (.NOT. (BA > CD_ZERO)) RETURN
    CALL CD_HermiteCable_Shapes(xi, l0, H, dH)
    r_xi = CD_ZERO
    v_xi = CD_ZERO
    DO s = 1, 4
      r_xi = r_xi + dH(s)*qe(3*(s - 1) + 1:3*s)
      v_xi = v_xi + dH(s)*ve(3*(s - 1) + 1:3*s)
    END DO
    speed = SQRT(SUM(r_xi*r_xi))
    IF (.NOT. CD_Is_Finite(speed) .OR. speed <= 1.0e-12_wp*MAX(CD_ONE, l0)) THEN
      CALL fail('collapsed or non-finite element tangent'); RETURN
    END IF
    resultant = BA*DOT_PRODUCT(r_xi/speed, v_xi)/l0
    IF (.NOT. CD_Is_Finite(resultant)) CALL fail('non-finite damping resultant')

  CONTAINS
    SUBROUTINE fail(msg)
      CHARACTER(*), INTENT(IN) :: msg
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Axial_Damping_Resultant: '//msg
      resultant = CD_ZERO
    END SUBROUTINE fail
  END SUBROUTINE CD_HermiteCable_Axial_Damping_Resultant

  SUBROUTINE CD_HermiteCable_Axial_Damping_Element(qe, ve, l0, BA, fdamp, Kq, Kv, ErrStat, ErrMsg, &
                                                   quadrature_order)
    !! Consistent axial Kelvin-Voigt contribution for one cubic-Hermite element.
    !! The material coordinate is reference arc length s0 = l0*xi. With
    !! r_xi = dr/dxi, t = r_xi/|r_xi| and v_xi = dv/dxi,
    !!
    !!   strain_rate = t dot v_xi / l0,
    !!   fdamp_a     = integral BA*strain_rate*dH_a*t dxi.
    !!
    !! fdamp is an INTERNAL generalized force, so v dot fdamp is non-negative.
    !! Kq = d(fdamp)/dq and Kv = d(fdamp)/dv are exact analytical Jacobians.
    REAL(wp), INTENT(IN) :: qe(12), ve(12), l0, BA
    REAL(wp), INTENT(OUT) :: fdamp(12)
    REAL(wp), INTENT(OUT), OPTIONAL :: Kq(12, 12), Kv(12, 12)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: quadrature_order

    REAL(wp) :: gp(6), gw(6), H(4), dH(4)
    REAL(wp) :: r_xi(3), v_xi(3), tangent(3), pperp(3), projector(3, 3)
    REAL(wp) :: speed, strain_rate, scale_q, scale_v
    INTEGER :: order, ig, s, sp, c, cp, ia, ib

    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    fdamp = CD_ZERO
    IF (PRESENT(Kq)) Kq = CD_ZERO
    IF (PRESENT(Kv)) Kv = CD_ZERO
    order = 4
    IF (PRESENT(quadrature_order)) order = quadrature_order
    IF (.NOT. CD_Is_Finite(l0) .OR. l0 <= CD_ZERO .OR. &
        .NOT. CD_Is_Finite(BA) .OR. BA < CD_ZERO) THEN
      CALL fail('need finite l0 > 0 and BA >= 0'); RETURN
    END IF
    IF (.NOT. CD_All_Finite(qe) .OR. .NOT. CD_All_Finite(ve)) THEN
      CALL fail('qe and ve must be finite'); RETURN
    END IF
    IF (order < 1 .OR. order > 6) THEN
      CALL fail('quadrature_order must lie in [1,6]'); RETURN
    END IF
    IF (.NOT. (BA > CD_ZERO)) RETURN

    CALL CD_HermiteCable_Gauss_Rule(order, gp, gw)
    DO ig = 1, order
      CALL CD_HermiteCable_Shapes(gp(ig), l0, H, dH)
      r_xi = CD_ZERO
      v_xi = CD_ZERO
      DO s = 1, 4
        r_xi = r_xi + dH(s)*qe(3*(s - 1) + 1:3*s)
        v_xi = v_xi + dH(s)*ve(3*(s - 1) + 1:3*s)
      END DO
      speed = SQRT(SUM(r_xi*r_xi))
      IF (.NOT. CD_Is_Finite(speed) .OR. speed <= 1.0e-12_wp*MAX(CD_ONE, l0)) THEN
        CALL fail('collapsed or non-finite element tangent at a quadrature point'); RETURN
      END IF
      tangent = r_xi/speed
      strain_rate = DOT_PRODUCT(tangent, v_xi)/l0
      projector = CD_ZERO
      DO c = 1, 3
        projector(c, c) = CD_ONE
      END DO
      DO c = 1, 3
        DO cp = 1, 3
          projector(c, cp) = projector(c, cp) - tangent(c)*tangent(cp)
        END DO
      END DO
      pperp = MATMUL(projector, v_xi)

      DO s = 1, 4
        DO c = 1, 3
          ia = 3*(s - 1) + c
          fdamp(ia) = fdamp(ia) + gw(ig)*BA*strain_rate*dH(s)*tangent(c)
          DO sp = 1, 4
            DO cp = 1, 3
              ib = 3*(sp - 1) + cp
              IF (PRESENT(Kv)) THEN
                scale_v = gw(ig)*BA*dH(s)*dH(sp)/l0
                Kv(ia, ib) = Kv(ia, ib) + scale_v*tangent(c)*tangent(cp)
              END IF
              IF (PRESENT(Kq)) THEN
                scale_q = gw(ig)*BA*dH(s)*dH(sp)/speed
                Kq(ia, ib) = Kq(ia, ib) + scale_q*(strain_rate*projector(c, cp) + &
                                                   tangent(c)*pperp(cp)/l0)
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO

    IF (.NOT. CD_All_Finite(fdamp)) THEN
      CALL fail('non-finite damping force or Jacobian')
      RETURN
    END IF
    IF (PRESENT(Kq)) THEN
      IF (.NOT. CD_All_Finite(Kq)) THEN
        CALL fail('non-finite damping force or Jacobian'); RETURN
      END IF
    END IF
    IF (PRESENT(Kv)) THEN
      IF (.NOT. CD_All_Finite(Kv)) THEN
        CALL fail('non-finite damping force or Jacobian'); RETURN
      END IF
    END IF

  CONTAINS
    SUBROUTINE fail(msg)
      CHARACTER(*), INTENT(IN) :: msg
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Axial_Damping_Element: '//msg
      fdamp = CD_ZERO
      IF (PRESENT(Kq)) Kq = CD_ZERO
      IF (PRESENT(Kv)) Kv = CD_ZERO
    END SUBROUTINE fail
  END SUBROUTINE CD_HermiteCable_Axial_Damping_Element

  SUBROUTINE CD_HermiteCable_Drag_Element(qe, ve, l0, fluid, waterline_z, rho, diam, cdn, cdt, &
                                          fdrag, djq, djv, ErrStat, ErrMsg, &
                                          wv_h, wv_om, wv_k, wv_depth, wv_dir, wv_t, &
                                          wvc_amp, wvc_om, wvc_k, wvc_ph, wvc_e2kd)
    !! Consistent Gauss-integrated Morison quadratic drag on one cubic-Hermite bending-cable element.
    !! qe(12), ve(12) = [r1, m1, r2, m2] nodal positions and velocities; the interior velocity and
    !! tangent use the SAME cubic-Hermite interpolation as the geometry, so the material-tangent (m)
    !! DOFs -- which move the element interior while the endpoints may be fixed -- also carry drag.
    !! Outputs the 12-DOF generalized drag LOAD fdrag = integral N^T f_perlen ds (f_ext on the cable)
    !! and its force Jacobians djq = d fdrag/d qe, djv = d fdrag/d ve. The per-length drag and its
    !! rel-velocity / tangent Jacobians come from the analytical Morison primitive, so djq/djv are
    !! exact (FD-verified) and vanish smoothly at zero relative velocity (fixed-point safe).
    !!
    !! Submergence is evaluated per Gauss point (a point with r_z above waterline_z contributes no
    !! drag). For a FULLY submerged element -- the production regime for a lazy-wave power cable or a
    !! mooring, which hangs entirely below the free surface -- every Gauss point is wet and the
    !! integral is exact. For an element that STRADDLES the waterline this is a bounded Gauss-sampling
    !! approximation (the wetted length is quantised to the abscissae and the crossing-point
    !! derivative is dropped); surface-piercing elements do not get an exact wet-portion
    !! integration, which fully-submerged cables do not need.
    !! OPTIONAL WAVES: supplying a wave family adds the wave velocity at each wet Gauss point
    !! (Wheeler-stretched, z measured from waterline_z) to the fluid velocity the drag sees.
    !! Exactly one family with the shared wv_depth/wv_dir/wv_t: the scalar REGULAR trio
    !! wv_h/wv_om/wv_k (Airy), or the COMPONENT quartet wvc_amp/wvc_om/wvc_k/wvc_ph (an
    !! irregular long-crested table evaluated by CD_Component_Wave_Kinematics with total-eta
    !! Wheeler stretching). The wave-field SPATIAL GRADIENT is neglected in djq --
    !! the standard explicit-field treatment (the residual is exact; the tangent omits a term of
    !! order k*H*omega, far below the structural stiffness), so Newton convergence is preserved.
    REAL(wp), INTENT(IN) :: qe(12), ve(12), l0, fluid(3), waterline_z, rho, diam, cdn, cdt
    REAL(wp), INTENT(OUT) :: fdrag(12)
    ! djq/djv are OPTIONAL as a pair: a residual-only evaluation (a line-search trial, which
    ! never solves with a tangent) omits both and skips the 12x12 Jacobian accumulation --
    ! the dominant share of the drag-element cost. The force is unchanged.
    REAL(wp), INTENT(OUT), OPTIONAL :: djq(12, 12), djv(12, 12)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: wv_h, wv_om, wv_k, wv_depth, wv_dir, wv_t
    REAL(wp), INTENT(IN), OPTIONAL :: wvc_amp(:), wvc_om(:), wvc_k(:), wvc_ph(:)
    ! exp(-2 k_i depth) of the component table (see CD_Component_Wave_Kinematics)
    REAL(wp), INTENT(IN), OPTIONAL :: wvc_e2kd(:)

    INTEGER, PARAMETER :: NGP = 3
    REAL(wp), PARAMETER :: G3 = 0.7745966692414834_wp   ! sqrt(3/5)
    REAL(wp) :: gp(NGP), gw(NGP), H(4), dH(4)
    REAL(wp) :: rpos(3), vvel(3), rp(3), speed, tangent(3), urel(3), fpl(3), Jrel(3, 3), Jtan(3, 3)
    REAL(wp) :: scale, eta_w, uw(3), aw(3)
    LOGICAL :: has_wv, has_shared, has_trio, has_cmp, any_wv_arg, want_jac
    INTEGER :: g, s, c, sp, cp, a, b, es
    CHARACTER(200) :: em2

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    want_jac = PRESENT(djq)
    IF (want_jac .NEQV. PRESENT(djv)) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Drag_Element: djq and djv must be supplied together (or neither)'; RETURN
    END IF
    fdrag = CD_ZERO
    IF (want_jac) THEN
      djq = CD_ZERO; djv = CD_ZERO
    END IF
    has_shared = PRESENT(wv_depth) .AND. PRESENT(wv_dir) .AND. PRESENT(wv_t)
    has_trio = PRESENT(wv_h) .AND. PRESENT(wv_om) .AND. PRESENT(wv_k)
    has_cmp = PRESENT(wvc_amp) .AND. PRESENT(wvc_om) .AND. PRESENT(wvc_k) .AND. PRESENT(wvc_ph)
    any_wv_arg = PRESENT(wv_h) .OR. PRESENT(wv_om) .OR. PRESENT(wv_k) .OR. PRESENT(wv_depth) .OR. &
                 PRESENT(wv_dir) .OR. PRESENT(wv_t) .OR. PRESENT(wvc_amp) .OR. PRESENT(wvc_om) .OR. &
                 PRESENT(wvc_k) .OR. PRESENT(wvc_ph)
    has_wv = has_shared .AND. (has_trio .NEQV. has_cmp)
    IF (any_wv_arg .AND. .NOT. has_wv) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Drag_Element: wave arguments must form exactly one family '// &
               '(shared depth/dir/t plus the regular trio OR the component quartet)'; RETURN
    END IF
    IF (has_wv .AND. has_cmp) THEN
      ! component-table validity BEFORE the wet/dry Gauss loop (config validity must not
      ! depend on the element's wet/dry state; the kernel re-checks per call)
      IF (SIZE(wvc_amp) < 1 .OR. SIZE(wvc_om) /= SIZE(wvc_amp) .OR. SIZE(wvc_k) /= SIZE(wvc_amp) .OR. &
          SIZE(wvc_ph) /= SIZE(wvc_amp)) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Drag_Element: component arrays must share one length >= 1'; RETURN
      END IF
      IF (.NOT. (CD_All_Finite(wvc_amp) .AND. ALL(wvc_amp >= CD_ZERO) .AND. &
                 CD_All_Finite(wvc_om) .AND. ALL(wvc_om > CD_ZERO) .AND. &
                 CD_All_Finite(wvc_k) .AND. ALL(wvc_k > CD_ZERO) .AND. &
                 CD_All_Finite(wvc_ph))) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_Drag_Element: component data must have finite positive omega/k, '// &
                 'non-negative amplitude, finite phase'; RETURN
      END IF
    END IF
    IF (.NOT. CD_All_Finite(qe) .OR. .NOT. CD_All_Finite(ve) .OR. .NOT. CD_Is_Finite(l0) .OR. &
        .NOT. CD_All_Finite(fluid) .OR. .NOT. CD_Is_Finite(waterline_z)) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Drag_Element: inputs must be finite'; RETURN
    END IF
    IF (l0 <= CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_Drag_Element: element length must be positive'; RETURN
    END IF
    ! Validate the hydrodynamic scalars BEFORE the wet/dry Gauss loop: a fully dry element cycles
    ! every Gauss point and would otherwise skip the only scalar check (inside the per-length
    ! primitive), silently returning OK for a malformed drag configuration. Config validity must not
    ! depend on the element's wet/dry state.
    IF (.NOT. CD_Is_Finite(rho) .OR. rho <= CD_ZERO .OR. .NOT. CD_Is_Finite(diam) .OR. diam <= CD_ZERO .OR. &
        .NOT. CD_Is_Finite(cdn) .OR. cdn < CD_ZERO .OR. .NOT. CD_Is_Finite(cdt) .OR. cdt < CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_Drag_Element: need rho>0, diam>0, cdn>=0, cdt>=0 (finite)'; RETURN
    END IF
    ! 3-point Gauss-Legendre on [0,1].
    gp = [0.5_wp*(CD_ONE - G3), 0.5_wp, 0.5_wp*(CD_ONE + G3)]
    gw = [5.0_wp/18.0_wp, 8.0_wp/18.0_wp, 5.0_wp/18.0_wp]

    DO g = 1, NGP
      CALL CD_HermiteCable_Shapes(gp(g), l0, H, dH)
      DO c = 1, 3
        rpos(c) = H(1)*qe(c) + H(2)*qe(3 + c) + H(3)*qe(6 + c) + H(4)*qe(9 + c)
        vvel(c) = H(1)*ve(c) + H(2)*ve(3 + c) + H(3)*ve(6 + c) + H(4)*ve(9 + c)
        rp(c) = dH(1)*qe(c) + dH(2)*qe(3 + c) + dH(3)*qe(6 + c) + dH(4)*qe(9 + c)
      END DO
      ! Wet/dry culling. Still water: a point above waterline_z is dry. With waves the free surface
      ! is the INSTANTANEOUS elevation eta(x, y, t): a point above still water but below a crest is
      ! wet (and carries drag against current + wave velocity); a point above eta is dry. The
      ! d(eta)/dq contribution is part of the documented neglected field gradient.
      IF (has_wv) THEN
        IF (has_cmp) THEN
          CALL CD_Component_Wave_Kinematics(rpos(1), rpos(2), rpos(3) - waterline_z, wv_t, &
                                            wv_depth, wv_dir, .TRUE., wvc_om, wvc_k, wvc_amp, wvc_ph, &
                                            eta_w, uw, aw, es, em2, inputs_validated=.TRUE., &
                                            exp_m2kd=wvc_e2kd)
        ELSE
          CALL CD_Airy_Wave_Kinematics_Precomputed(rpos(1), rpos(2), rpos(3) - waterline_z, wv_t, &
                                                   wv_h, wv_om, wv_k, wv_depth, wv_dir, .TRUE., &
                                                   eta_w, uw, aw, es, em2)
        END IF
        IF (es /= CD_HYDRO_OK) THEN
          ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Drag_Element: '//TRIM(em2); RETURN
        END IF
        IF (rpos(3) - waterline_z > eta_w) CYCLE       ! above the instantaneous surface: dry
      ELSE
        IF (rpos(3) > waterline_z) CYCLE               ! dry Gauss point: no drag, no Jacobian
      END IF
      speed = SQRT(rp(1)**2 + rp(2)**2 + rp(3)**2)     ! |dr/dxi| = deformed arc length per unit xi
      IF (speed <= 1.0e-12_wp) THEN
        ErrStat = CD_HCDYN_NOCONVERGE
        ErrMsg = 'CD_HermiteCable_Drag_Element: degenerate centreline (|dr/dxi| ~ 0)'; RETURN
      END IF
      tangent = rp/speed
      urel = fluid - vvel
      IF (has_wv) urel = urel + uw
      ! Pass the raw dr/dxi as the tangent vector: the primitive normalises it and returns Jtan =
      ! d f_perlen / d(dr/dxi), so d f_perlen/d qe folds in through d(dr/dxi)/d qe = dH directly.
      IF (want_jac) THEN
        CALL CD_Morison_Drag_Per_Length_Jac(urel, rp, rho, diam, cdn, cdt, fpl, Jrel, Jtan, es, em2)
      ELSE
        ! same primitive without the Jacobian pair: the force statements are identical, so
        ! the residual-only trial force is bit-for-bit the with-Jacobian force.
        CALL CD_Morison_Drag_Per_Length_Jac(urel, rp, rho, diam, cdn, cdt, fpl, ErrStat=es, ErrMsg=em2)
      END IF
      IF (es /= CD_HYDRO_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Drag_Element: '//TRIM(em2); RETURN
      END IF
      scale = gw(g)*speed                              ! ds = |dr/dxi| dxi (deformed cable length)
      DO s = 1, 4
        DO c = 1, 3
          a = 3*(s - 1) + c
          fdrag(a) = fdrag(a) + scale*H(s)*fpl(c)
        END DO
      END DO
      IF (.NOT. want_jac) CYCLE
      DO s = 1, 4
        DO c = 1, 3
          a = 3*(s - 1) + c
          DO sp = 1, 4
            DO cp = 1, 3
              b = 3*(sp - 1) + cp
              ! velocity: d urel/d ve(b) = -H(sp), scaled by the deformed length measure
              djv(a, b) = djv(a, b) - scale*H(s)*H(sp)*Jrel(c, cp)
              ! position: d[speed * f_perlen_c]/d qe(b) with d(dr/dxi)/d qe = dH(sp) --
              !   speed * Jtan(c,cp) * dH(sp)  (force direction/tangent)
              ! + f_perlen_c * d(speed)/d qe(b),  d(speed)/d qe(b) = tangent(cp) dH(sp)
              djq(a, b) = djq(a, b) + gw(g)*H(s)*dH(sp)*(speed*Jtan(c, cp) + fpl(c)*tangent(cp))
            END DO
          END DO
        END DO
      END DO
    END DO

    IF (.NOT. CD_All_Finite(fdrag)) THEN
      ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Drag_Element: non-finite drag load/Jacobian'
      RETURN
    END IF
    IF (want_jac) THEN
      IF (.NOT. CD_All_Finite(djq) .OR. .NOT. CD_All_Finite(djv)) THEN
        ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_Drag_Element: non-finite drag load/Jacobian'
      END IF
    END IF
  END SUBROUTINE CD_HermiteCable_Drag_Element

  SUBROUTINE CD_HermiteCable_AddedMass_Element(qe, l0, waterline_z, rho, diam, can, cat, Ma, ErrStat, ErrMsg)
    !! Consistent Gauss-integrated Morison ADDED MASS on one cubic-Hermite bending-cable element:
    !!   Ma = integral N^T m_a(t) N ds,   m_a(t) = rho A (Can (I - t t^T) + Cat t t^T),
    !! with A = pi d^2/4, t the unit tangent of the deformed centreline, and ds = |dr/dxi| dxi the
    !! DEFORMED length measure -- the same interpolation, wet/dry culling, and measure as the drag
    !! element, so the two hydro effects are mutually consistent. Symmetric positive semi-definite
    !! (zero for a dry element). The submergence note on CD_HermiteCable_Drag_Element applies.
    REAL(wp), INTENT(IN) :: qe(12), l0, waterline_z, rho, diam, can, cat
    REAL(wp), INTENT(OUT) :: Ma(12, 12)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    ! 4-point Gauss: the H_i H_j mass integrand is degree 6, which 4-pt integrates exactly (the same
    ! rule the consistent structural mass uses) -- 3-pt would under-integrate the added mass.
    INTEGER, PARAMETER :: NGP = 4
    REAL(wp), PARAMETER :: PI_L = 3.14159265358979323846_wp
    REAL(wp) :: gp(NGP), gw(NGP), H(4), dH(4)
    REAL(wp) :: rpos(3), rp(3), speed, tangent(3), ma3(3, 3), area, scale
    INTEGER :: g, s, c, sp, cp, a, b

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    Ma = CD_ZERO
    IF (.NOT. CD_All_Finite(qe) .OR. .NOT. CD_Is_Finite(l0) .OR. .NOT. CD_Is_Finite(waterline_z)) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_AddedMass_Element: inputs must be finite'; RETURN
    END IF
    IF (l0 <= CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_AddedMass_Element: element length must be positive'; RETURN
    END IF
    IF (.NOT. CD_Is_Finite(rho) .OR. rho <= CD_ZERO .OR. .NOT. CD_Is_Finite(diam) .OR. diam <= CD_ZERO .OR. &
        .NOT. CD_Is_Finite(can) .OR. can < CD_ZERO .OR. .NOT. CD_Is_Finite(cat) .OR. cat < CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_AddedMass_Element: need rho>0, diam>0, can>=0, cat>=0 (finite)'; RETURN
    END IF
    area = 0.25_wp*PI_L*diam*diam
    gp = [0.5_wp*(CD_ONE - 0.8611363115940526_wp), 0.5_wp*(CD_ONE - 0.3399810435848563_wp), &
          0.5_wp*(CD_ONE + 0.3399810435848563_wp), 0.5_wp*(CD_ONE + 0.8611363115940526_wp)]
    gw = 0.5_wp*[0.3478548451374538_wp, 0.6521451548625461_wp, &
                 0.6521451548625461_wp, 0.3478548451374538_wp]

    DO g = 1, NGP
      CALL CD_HermiteCable_Shapes(gp(g), l0, H, dH)
      DO c = 1, 3
        rpos(c) = H(1)*qe(c) + H(2)*qe(3 + c) + H(3)*qe(6 + c) + H(4)*qe(9 + c)
        rp(c) = dH(1)*qe(c) + dH(2)*qe(3 + c) + dH(3)*qe(6 + c) + dH(4)*qe(9 + c)
      END DO
      IF (rpos(3) > waterline_z) CYCLE                 ! dry Gauss point: no added mass
      speed = SQRT(rp(1)**2 + rp(2)**2 + rp(3)**2)
      IF (speed <= 1.0e-12_wp) THEN
        ErrStat = CD_HCDYN_NOCONVERGE
        ErrMsg = 'CD_HermiteCable_AddedMass_Element: degenerate centreline (|dr/dxi| ~ 0)'; RETURN
      END IF
      tangent = rp/speed
      ! m_a(t) = rho A (Can I + (Cat - Can) t t^T)
      DO cp = 1, 3
        DO c = 1, 3
          ma3(c, cp) = rho*area*(cat - can)*tangent(c)*tangent(cp)
        END DO
        ma3(cp, cp) = ma3(cp, cp) + rho*area*can
      END DO
      scale = gw(g)*speed                              ! ds = |dr/dxi| dxi (deformed cable length)
      DO s = 1, 4
        DO c = 1, 3
          a = 3*(s - 1) + c
          DO sp = 1, 4
            DO cp = 1, 3
              b = 3*(sp - 1) + cp
              Ma(a, b) = Ma(a, b) + scale*H(s)*H(sp)*ma3(c, cp)
            END DO
          END DO
        END DO
      END DO
    END DO

    ! Fail closed on an overflowed accumulation (finite but huge rho/diam/coefficients): a public
    ! element must never return CD_HCDYN_OK with a non-finite matrix.
    IF (.NOT. CD_All_Finite(Ma)) THEN
      ErrStat = CD_HCDYN_NOCONVERGE
      ErrMsg = 'CD_HermiteCable_AddedMass_Element: non-finite added-mass matrix'
    END IF
  END SUBROUTINE CD_HermiteCable_AddedMass_Element

  SUBROUTINE CD_HermiteCable_FK_Element(qe, l0, waterline_z, rho, diam, can, cat, &
                                        wv_h, wv_om, wv_k, wv_depth, wv_dir, wv_t, &
                                        ffk, fjq, ErrStat, ErrMsg, &
                                        wvc_amp, wvc_om, wvc_k, wvc_ph, hf_a, wvc_e2kd)
    !! Consistent Gauss-integrated FROUDE-KRYLOV + FLUID-INERTIA wave load on one cubic-Hermite
    !! bending-cable element: per wet Gauss point the Airy wave acceleration a_w (Wheeler-stretched,
    !! z measured from waterline_z) drives
    !!   f = rho A [ (1 + Can) a_w + ((1 + Cat) - (1 + Can)) (a_w . t) t ],
    !! distributed as ffk = integral N^T f ds over the deformed wet length (same interpolation,
    !! culling, and measure as the drag/added-mass elements). fjq = d ffk / d qe carries the
    !! GEOMETRIC tangent chain (through t(dr/dxi) and the length measure); the wave-field SPATIAL
    !! GRADIENT is neglected -- the standard explicit-field treatment (residual exact, tangent
    !! omits an O(k H omega^2) term far below the structural stiffness).
    !!
    !! FLUID-ACCELERATION SOURCES (exactly one): the scalar REGULAR trio wv_h/wv_om/wv_k,
    !! the COMPONENT quartet wvc_amp/wvc_om/wvc_k/wvc_ph (an irregular long-crested table
    !! evaluated by CD_Component_Wave_Kinematics with total-eta Wheeler stretching), or the
    !! HELD field hf_a -- a per-element CONSTANT fluid acceleration prescribed by an
    !! external host (sampled at the line nodes and node-averaged per element; frozen over
    !! the step, the documented weak-coupling cadence). With hf_a the caller's waterline_z
    !! already carries the local free-surface elevation, so the wet check uses it directly
    !! (eta_w = 0 relative) and wv_depth/wv_dir/wv_t are unread (pass zeros).
    !! wv_depth/wv_dir/wv_t are shared by both wave families.
    REAL(wp), INTENT(IN) :: qe(12), l0, waterline_z, rho, diam, can, cat
    REAL(wp), INTENT(IN), OPTIONAL :: wv_h, wv_om, wv_k
    REAL(wp), INTENT(IN) :: wv_depth, wv_dir, wv_t
    REAL(wp), INTENT(OUT) :: ffk(12)
    ! fjq is OPTIONAL: a residual-only evaluation (a line-search trial, which never solves
    ! with a tangent) omits it and skips the 12x12 Jacobian accumulation. ffk is unchanged.
    REAL(wp), INTENT(OUT), OPTIONAL :: fjq(12, 12)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: wvc_amp(:), wvc_om(:), wvc_k(:), wvc_ph(:)
    ! exp(-2 k_i depth) of the component table (see CD_Component_Wave_Kinematics)
    REAL(wp), INTENT(IN), OPTIONAL :: wvc_e2kd(:)
    REAL(wp), INTENT(IN), OPTIONAL :: hf_a(3)

    INTEGER, PARAMETER :: NGP = 3
    REAL(wp), PARAMETER :: G3 = 0.7745966692414834_wp   ! sqrt(3/5)
    REAL(wp), PARAMETER :: PI_L = 3.14159265358979323846_wp
    REAL(wp) :: gp(NGP), gw(NGP), H(4), dH(4)
    REAL(wp) :: rpos(3), rp(3), speed, tangent(3), eta_w, uw(3), aw(3), area, alpha
    REAL(wp) :: fpl(3), Jtan(3, 3), proj(3, 3), projk(3), dalpha
    LOGICAL :: has_trio, has_cmp, has_held, want_jac
    INTEGER :: g, s, c, sp, cp, a, b, k, es
    CHARACTER(200) :: em2

    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    want_jac = PRESENT(fjq)
    ffk = CD_ZERO
    IF (want_jac) fjq = CD_ZERO
    has_trio = PRESENT(wv_h) .AND. PRESENT(wv_om) .AND. PRESENT(wv_k)
    has_cmp = PRESENT(wvc_amp) .AND. PRESENT(wvc_om) .AND. PRESENT(wvc_k) .AND. PRESENT(wvc_ph)
    IF ((.NOT. has_trio) .AND. (PRESENT(wv_h) .OR. PRESENT(wv_om) .OR. PRESENT(wv_k))) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_FK_Element: the regular-wave trio must be supplied together'; RETURN
    END IF
    IF ((.NOT. has_cmp) .AND. (PRESENT(wvc_amp) .OR. PRESENT(wvc_om) .OR. PRESENT(wvc_k) .OR. &
                               PRESENT(wvc_ph))) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_FK_Element: the component-wave quartet must be supplied together'; RETURN
    END IF
    has_held = PRESENT(hf_a)
    IF (COUNT([has_trio, has_cmp, has_held]) /= 1) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_FK_Element: supply exactly one fluid-acceleration source '// &
               '(regular trio, component quartet, OR the held field hf_a)'; RETURN
    END IF
    IF (has_held) THEN
      IF (.NOT. CD_All_Finite(hf_a)) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_FK_Element: the held fluid acceleration must be finite'; RETURN
      END IF
    END IF
    IF (.NOT. CD_All_Finite(qe) .OR. .NOT. CD_Is_Finite(l0) .OR. .NOT. CD_Is_Finite(waterline_z) .OR. &
        .NOT. CD_Is_Finite(wv_t)) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_FK_Element: inputs must be finite'; RETURN
    END IF
    IF (l0 <= CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT; ErrMsg = 'CD_HermiteCable_FK_Element: element length must be positive'; RETURN
    END IF
    ! Validate all scalars BEFORE the wet/dry Gauss loop (config validity must not depend on the
    ! element's wet/dry state).
    IF (.NOT. CD_Is_Finite(rho) .OR. rho <= CD_ZERO .OR. .NOT. CD_Is_Finite(diam) .OR. diam <= CD_ZERO .OR. &
        .NOT. CD_Is_Finite(can) .OR. can < CD_ZERO .OR. .NOT. CD_Is_Finite(cat) .OR. cat < CD_ZERO) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'CD_HermiteCable_FK_Element: need rho>0, diam>0, can>=0, cat>=0 (finite)'; RETURN
    END IF
    IF (.NOT. has_held) THEN
      IF (.NOT. CD_Is_Finite(wv_depth) .OR. wv_depth <= CD_ZERO .OR. .NOT. CD_Is_Finite(wv_dir)) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_FK_Element: need wave depth > 0 and finite direction'; RETURN
      END IF
    END IF
    IF (has_held) THEN
      CONTINUE   ! the held field carries no wave configuration to validate
    ELSE IF (has_trio) THEN
      IF (.NOT. CD_Is_Finite(wv_h) .OR. wv_h <= CD_ZERO .OR. .NOT. CD_Is_Finite(wv_om) .OR. &
          wv_om <= CD_ZERO .OR. .NOT. CD_Is_Finite(wv_k) .OR. wv_k <= CD_ZERO) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_FK_Element: need wave height/omega/k > 0 (finite)'; RETURN
      END IF
    ELSE
      ! component-table validity BEFORE the wet/dry Gauss loop (same rule as the scalars;
      ! the kernel re-checks, but a fully dry element must not skip the config check)
      IF (SIZE(wvc_amp) < 1 .OR. SIZE(wvc_om) /= SIZE(wvc_amp) .OR. SIZE(wvc_k) /= SIZE(wvc_amp) .OR. &
          SIZE(wvc_ph) /= SIZE(wvc_amp)) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_FK_Element: component arrays must share one length >= 1'; RETURN
      END IF
      IF (.NOT. (CD_All_Finite(wvc_amp) .AND. ALL(wvc_amp >= CD_ZERO) .AND. &
                 CD_All_Finite(wvc_om) .AND. ALL(wvc_om > CD_ZERO) .AND. &
                 CD_All_Finite(wvc_k) .AND. ALL(wvc_k > CD_ZERO) .AND. &
                 CD_All_Finite(wvc_ph))) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'CD_HermiteCable_FK_Element: component data must have finite positive omega/k, '// &
                 'non-negative amplitude, finite phase'; RETURN
      END IF
    END IF
    area = 0.25_wp*PI_L*diam*diam
    gp = [0.5_wp*(CD_ONE - G3), 0.5_wp, 0.5_wp*(CD_ONE + G3)]
    gw = [5.0_wp/18.0_wp, 8.0_wp/18.0_wp, 5.0_wp/18.0_wp]

    DO g = 1, NGP
      CALL CD_HermiteCable_Shapes(gp(g), l0, H, dH)
      DO c = 1, 3
        rpos(c) = H(1)*qe(c) + H(2)*qe(3 + c) + H(3)*qe(6 + c) + H(4)*qe(9 + c)
        rp(c) = dH(1)*qe(c) + dH(2)*qe(3 + c) + dH(3)*qe(6 + c) + dH(4)*qe(9 + c)
      END DO
      ! Wet/dry culling on the INSTANTANEOUS free surface eta(x, y, t): a point above still water
      ! but below a crest is wet and carries the wave load; a point above eta is dry (the Wheeler
      ! kinematics are zero there anyway). d(eta)/dq is part of the neglected field gradient.
      ! With the HELD field the caller's waterline_z already carries the sampled local elevation,
      ! so the wet check is against it directly (eta_w = 0 relative).
      IF (has_held) THEN
        eta_w = CD_ZERO
        aw = hf_a
        es = CD_HYDRO_OK
      ELSE IF (has_cmp) THEN
        CALL CD_Component_Wave_Kinematics(rpos(1), rpos(2), rpos(3) - waterline_z, wv_t, &
                                          wv_depth, wv_dir, .TRUE., wvc_om, wvc_k, wvc_amp, wvc_ph, &
                                          eta_w, uw, aw, es, em2, inputs_validated=.TRUE., &
                                          exp_m2kd=wvc_e2kd)
      ELSE
        CALL CD_Airy_Wave_Kinematics_Precomputed(rpos(1), rpos(2), rpos(3) - waterline_z, wv_t, &
                                                 wv_h, wv_om, wv_k, wv_depth, wv_dir, .TRUE., &
                                                 eta_w, uw, aw, es, em2)
      END IF
      IF (es /= CD_HYDRO_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_FK_Element: '//TRIM(em2); RETURN
      END IF
      IF (rpos(3) - waterline_z > eta_w) CYCLE         ! above the instantaneous surface: dry
      speed = SQRT(rp(1)**2 + rp(2)**2 + rp(3)**2)
      IF (speed <= 1.0e-12_wp) THEN
        ErrStat = CD_HCDYN_NOCONVERGE
        ErrMsg = 'CD_HermiteCable_FK_Element: degenerate centreline (|dr/dxi| ~ 0)'; RETURN
      END IF
      tangent = rp/speed
      alpha = DOT_PRODUCT(aw, tangent)
      fpl = rho*area*((CD_ONE + can)*aw + (cat - can)*alpha*tangent)
      IF (want_jac) THEN
        ! Jtan = d fpl / d(dr/dxi): through d t/d rp = (I - t t^T)/speed. Column k uses the
        ! projected direction projk; d fpl = rho A (cat-can) (dalpha t + alpha projk),
        ! dalpha = aw . projk.
        DO k = 1, 3
          DO c = 1, 3
            proj(c, k) = -tangent(c)*tangent(k)/speed
          END DO
          proj(k, k) = proj(k, k) + CD_ONE/speed
        END DO
        DO k = 1, 3
          projk = proj(:, k)
          dalpha = DOT_PRODUCT(aw, projk)
          Jtan(:, k) = rho*area*(cat - can)*(dalpha*tangent + alpha*projk)
        END DO
      END IF
      DO s = 1, 4
        DO c = 1, 3
          a = 3*(s - 1) + c
          ffk(a) = ffk(a) + gw(g)*speed*H(s)*fpl(c)
        END DO
      END DO
      IF (.NOT. want_jac) CYCLE
      DO s = 1, 4
        DO c = 1, 3
          a = 3*(s - 1) + c
          DO sp = 1, 4
            DO cp = 1, 3
              b = 3*(sp - 1) + cp
              ! d[speed * fpl_c]/d qe(b): the tangent chain plus the length-measure derivative.
              fjq(a, b) = fjq(a, b) + gw(g)*H(s)*dH(sp)*(speed*Jtan(c, cp) + fpl(c)*tangent(cp))
            END DO
          END DO
        END DO
      END DO
    END DO

    IF (.NOT. CD_All_Finite(ffk)) THEN
      ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_FK_Element: non-finite wave load/Jacobian'
      RETURN
    END IF
    IF (want_jac) THEN
      IF (.NOT. CD_All_Finite(fjq)) THEN
        ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'CD_HermiteCable_FK_Element: non-finite wave load/Jacobian'
      END IF
    END IF
  END SUBROUTINE CD_HermiteCable_FK_Element

  SUBROUTINE assemble_added_mass(model, q_cfg, ErrStat, ErrMsg)
    !! Global consistent added-mass matrix at configuration q_cfg (zero when disabled),
    !! written into the model's own ws_Mamb band buffer through the single model dummy
    !! (no subobject argument-associated alongside the model -- the aliasing rule).
    !! q_cfg may be a model workspace vector; it is only read.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    REAL(wp), INTENT(IN) :: q_cfg(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: e, es
    REAL(wp) :: qe(12), Me(12, 12), wl_e
    CHARACTER(200) :: em2
    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    model%ws_Mamb = CD_ZERO
    IF (.NOT. model%has_am) RETURN
    DO e = 1, model%ne
      qe(1:3) = q_cfg(6*(e - 1) + 1:6*(e - 1) + 3); qe(4:6) = q_cfg(6*(e - 1) + 4:6*(e - 1) + 6)
      qe(7:9) = q_cfg(6*e + 1:6*e + 3); qe(10:12) = q_cfg(6*e + 4:6*e + 6)
      ! Held-field mode: the inertia matrix must wet against the SAME free surface the
      ! drag and Froude-Krylov use (the host-sampled elevation, node-averaged per
      ! element) -- against the init-time still-water am_wl a surface-piercing element
      ! would carry the wrong added mass exactly when the host drives waves.
      wl_e = model%am_wl
      IF (model%has_held_field) wl_e = 0.5_wp*(model%hf_wl(e) + model%hf_wl(e + 1))
      CALL CD_HermiteCable_AddedMass_Element(qe, model%l0(e), wl_e, model%am_rho, &
                                             model%am_diam(e), model%am_can(e), model%am_cat(e), Me, es, em2)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = TRIM(em2); RETURN
      END IF
      CALL scatter12_band(model%ws_Mamb, Me, e)
    END DO
    ! Every element matrix can be finite while two adjacent contributions
    ! overflow at their shared-node band entry.  Validate the assembled global
    ! operator before it is added to the structural mass or factorised.
    IF (.NOT. CD_All_Finite(model%ws_Mamb)) THEN
      ErrStat = CD_HCDYN_NOCONVERGE
      ErrMsg = 'assemble_added_mass: non-finite global added-mass matrix'
    END IF
  END SUBROUTINE assemble_added_mass

  SUBROUTINE build_rigid_bases(model, directions, bases, ErrStat, ErrMsg)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: directions(3, 2)
    REAL(wp), INTENT(OUT) :: bases(3, 3, 2)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: iend, es
    CHARACTER(200) :: em
    bases = CD_ZERO
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      CALL CD_EndConn_Basis(directions(:, iend), bases(:, :, iend), es, em)
      IF (es /= CD_ENDCONN_OK) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'invalid rigid end direction: '//TRIM(em)
        RETURN
      END IF
    END DO
  END SUBROUTINE build_rigid_bases

  SUBROUTINE project_rigid_tangents(model, directions, q, ErrStat, ErrMsg)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: directions(3, 2)
    REAL(wp), INTENT(INOUT) :: q(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: iend, base, es
    REAL(wp) :: projected(3)
    CHARACTER(200) :: em
    ErrStat = CD_HCDYN_OK
    ErrMsg = ''
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (iend == 1) THEN
        base = 3
      ELSE
        base = 6*(model%nn - 1) + 3
      END IF
      CALL CD_EndConn_Project(q(base + 1:base + 3), directions(:, iend), projected, es, em)
      IF (es /= CD_ENDCONN_OK) THEN
        ErrStat = CD_HCDYN_NOCONVERGE
        ErrMsg = 'rigid endpoint projection failed: '//TRIM(em)
        RETURN
      END IF
      q(base + 1:base + 3) = projected
    END DO
  END SUBROUTINE project_rigid_tangents

  SUBROUTINE transform_rigid_vector_to_local(model, bases, vector)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: bases(3, 3, 2)
    REAL(wp), INTENT(INOUT) :: vector(:)
    INTEGER :: iend, base
    REAL(wp) :: local_values(3)
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (iend == 1) THEN
        base = 3
      ELSE
        base = 6*(model%nn - 1) + 3
      END IF
      local_values = MATMUL(TRANSPOSE(bases(:, :, iend)), vector(base + 1:base + 3))
      vector(base + 1:base + 3) = local_values
    END DO
  END SUBROUTINE transform_rigid_vector_to_local

  SUBROUTINE transform_rigid_vector_to_global(model, bases, vector)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: bases(3, 3, 2)
    REAL(wp), INTENT(INOUT) :: vector(:)
    INTEGER :: iend, base
    REAL(wp) :: global_values(3)
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (iend == 1) THEN
        base = 3
      ELSE
        base = 6*(model%nn - 1) + 3
      END IF
      global_values = MATMUL(bases(:, :, iend), vector(base + 1:base + 3))
      vector(base + 1:base + 3) = global_values
    END DO
  END SUBROUTINE transform_rigid_vector_to_global

  PURE SUBROUTINE newmark_acceleration(c_a, q, q_n, c_v, v_n, c_a2, a_n, a)
    !! a = c_a*(q - q_n) - c_v*v_n - c_a2*a_n (a distinct from the inputs).
    REAL(wp), INTENT(IN) :: c_a, c_v, c_a2, q(:), q_n(:), v_n(:), a_n(:)
    REAL(wp), INTENT(INOUT) :: a(:)
    a = c_a*(q - q_n) - c_v*v_n - c_a2*a_n
  END SUBROUTINE newmark_acceleration

  PURE SUBROUTINE newmark_velocity(dt, w_n, w_k, v_n, a_n, a_k, v)
    !! v = v_n + dt*(w_n*a_n + w_k*a_k) (v distinct from the inputs).
    REAL(wp), INTENT(IN) :: dt, w_n, w_k, v_n(:), a_n(:), a_k(:)
    REAL(wp), INTENT(INOUT) :: v(:)
    v = v_n + dt*(w_n*a_n + w_k*a_k)
  END SUBROUTINE newmark_velocity

  PURE SUBROUTINE vector_copy(x, y)
    !! y = x (distinct arrays).
    REAL(wp), INTENT(IN) :: x(:)
    REAL(wp), INTENT(INOUT) :: y(:)
    y = x
  END SUBROUTINE vector_copy

  PURE SUBROUTINE vector_blend(a, x, b, y, out)
    !! out = a*x + b*y (out distinct from x and y).
    REAL(wp), INTENT(IN) :: a, b, x(:), y(:)
    REAL(wp), INTENT(INOUT) :: out(:)
    out = a*x + b*y
  END SUBROUTINE vector_blend

  PURE SUBROUTINE vector_sum2(x, y, out)
    !! out = x + y (out distinct from x and y).
    REAL(wp), INTENT(IN) :: x(:), y(:)
    REAL(wp), INTENT(INOUT) :: out(:)
    out = x + y
  END SUBROUTINE vector_sum2

  PURE SUBROUTINE vector_sum3(x, b, y, c, z, out)
    !! out = x + b*y + c*z (out distinct from x, y and z).
    REAL(wp), INTENT(IN) :: b, c, x(:), y(:), z(:)
    REAL(wp), INTENT(INOUT) :: out(:)
    out = x + b*y + c*z
  END SUBROUTINE vector_sum3

  PURE SUBROUTINE band_copy(x, y)
    !! y = x for band storage (distinct arrays).
    REAL(wp), INTENT(IN) :: x(:, :)
    REAL(wp), INTENT(INOUT) :: y(:, :)
    y = x
  END SUBROUTINE band_copy

  PURE SUBROUTINE band_add(x, y)
    !! y = y + x for band storage (distinct arrays).
    REAL(wp), INTENT(IN) :: x(:, :)
    REAL(wp), INTENT(INOUT) :: y(:, :)
    y = y + x
  END SUBROUTINE band_add

  PURE SUBROUTINE band_combine(a, x, b, y, c, z, out)
    !! out = a*x + b*y + c*z for band storage (out distinct from x, y, z).
    REAL(wp), INTENT(IN) :: a, b, c, x(:, :), y(:, :), z(:, :)
    REAL(wp), INTENT(INOUT) :: out(:, :)
    out = a*x + b*y + c*z
  END SUBROUTINE band_combine

  SUBROUTINE transform_rigid_band(model, bases, matrix)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: bases(3, 3, 2)
    REAL(wp), INTENT(INOUT) :: matrix(:, :)
    INTEGER :: iend, first, row, col, j
    REAL(wp) :: old_values(3), new_values(3)
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (iend == 1) THEN
        first = 4
      ELSE
        first = 6*(model%nn - 1) + 4
      END IF
      DO row = MAX(1, first - KU_D), MIN(model%ndof, first + 2 + KL_D)
        DO j = 1, 3
          old_values(j) = dynamic_band_entry(model, matrix, row, first + j - 1)
        END DO
        new_values = MATMUL(old_values, bases(:, :, iend))
        DO j = 1, 3
          CALL set_dynamic_band_entry(model, matrix, row, first + j - 1, new_values(j))
        END DO
      END DO
      DO col = MAX(1, first - KL_D), MIN(model%ndof, first + 2 + KU_D)
        DO j = 1, 3
          old_values(j) = dynamic_band_entry(model, matrix, first + j - 1, col)
        END DO
        new_values = MATMUL(TRANSPOSE(bases(:, :, iend)), old_values)
        DO j = 1, 3
          CALL set_dynamic_band_entry(model, matrix, first + j - 1, col, new_values(j))
        END DO
      END DO
    END DO
  END SUBROUTINE transform_rigid_band

  PURE REAL(wp) FUNCTION dynamic_band_entry(model, matrix, row, col) RESULT(value)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(IN) :: matrix(:, :)
    INTEGER, INTENT(IN) :: row, col
    value = CD_ZERO
    IF (row < 1 .OR. row > model%ndof .OR. col < 1 .OR. col > model%ndof) RETURN
    IF (row - col < -KU_D .OR. row - col > KL_D) RETURN
    value = matrix(KL_D + KU_D + 1 + row - col, col)
  END FUNCTION dynamic_band_entry

  SUBROUTINE set_dynamic_band_entry(model, matrix, row, col, value)
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    REAL(wp), INTENT(INOUT) :: matrix(:, :)
    INTEGER, INTENT(IN) :: row, col
    REAL(wp), INTENT(IN) :: value
    IF (row < 1 .OR. row > model%ndof .OR. col < 1 .OR. col > model%ndof) RETURN
    IF (row - col < -KU_D .OR. row - col > KL_D) RETURN
    matrix(KL_D + KU_D + 1 + row - col, col) = value
  END SUBROUTINE set_dynamic_band_entry

  SUBROUTINE solve_consistent_acceleration(model, ErrStat, ErrMsg, pres_dofs, pres_a)
    !! Set model%a to the consistent acceleration M a = f_ext(q,v) - f_int(q) on the free DOFs at the
    !! current state. Fixed DOFs keep a = 0 unless pres_dofs/pres_a supplies a prescribed acceleration;
    !! its M_fp*a_p contribution is then included in the free equations. Used at Init and whenever the
    !! load set changes (drag), so the gen-alpha predictor always starts from a physically consistent
    !! acceleration.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: model
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: pres_dofs(:)
    REAL(wp), INTENT(IN), OPTIONAL :: pres_a(:)
    INTEGER :: nf, np, i, ip, jp, gdof, es, iend, base
    LOGICAL :: has_pres
    REAL(wp) :: rigid_basis(3, 3, 2), old_acceleration(3), axial_acceleration
    CHARACTER(200) :: em2
    ErrStat = CD_HCDYN_OK; ErrMsg = ''
    ! This solve factors the mass into ws_Keffb: the carried step factor is gone.
    model%tr_valid = .FALSE.
    nf = COUNT(model%solve_mask)
    has_pres = PRESENT(pres_dofs) .AND. PRESENT(pres_a)
    IF ((PRESENT(pres_dofs) .OR. PRESENT(pres_a)) .AND. .NOT. has_pres) THEN
      ErrStat = CD_HCDYN_BADINPUT
      ErrMsg = 'solve_consistent_acceleration: pres_dofs and pres_a must be supplied together'; RETURN
    END IF
    np = 0
    IF (has_pres) THEN
      np = SIZE(pres_dofs)
      IF (SIZE(pres_a) /= np) THEN
        ErrStat = CD_HCDYN_BADINPUT
        ErrMsg = 'solve_consistent_acceleration: pres_a must match pres_dofs length'; RETURN
      END IF
      DO ip = 1, np
        gdof = pres_dofs(ip)
        IF (gdof < 1 .OR. gdof > model%ndof) THEN
          ErrStat = CD_HCDYN_BADINPUT
          ErrMsg = 'solve_consistent_acceleration: prescribed DOF out of range'; RETURN
        END IF
        IF (model%freemask(gdof)) THEN
          ErrStat = CD_HCDYN_BADINPUT
          ErrMsg = 'solve_consistent_acceleration: prescribed DOF must be in fixed_dofs'; RETURN
        END IF
        IF (.NOT. CD_Is_Finite(pres_a(ip))) THEN
          ErrStat = CD_HCDYN_BADINPUT
          ErrMsg = 'solve_consistent_acceleration: prescribed acceleration must be finite'; RETURN
        END IF
        DO jp = 1, ip - 1
          IF (pres_dofs(jp) == gdof) THEN
            ErrStat = CD_HCDYN_BADINPUT
            ErrMsg = 'solve_consistent_acceleration: prescribed DOFs must be unique'; RETURN
          END IF
        END DO
      END DO
    END IF
    ! Preserve the already committed acceleration of a moving rigid direction.
    ! Load setters may refresh the free acceleration between coupling steps; they
    ! must not erase the kinematic transverse acceleration established by the
    ! preceding constrained step.
    model%ws_vk = model%a
    model%a = CD_ZERO
    ! Assemble (and thereby VALIDATE) the residual even for a fully constrained model: skipping it
    ! would let Init accept a geometrically invalid seed, and Set_Drag report success on a drag
    ! configuration whose loads are non-finite -- deferring the failure to the next Step instead of
    ! failing closed (and rolling back) here. Only the mass solve is skipped when nf = 0.
    ! Uses the model's step workspace (allocated at Init before this is first called).
    CALL assemble_residual(model, model%q, model%v, model%t, .TRUE., .FALSE., es, em2)  ! ws_R = f_int - f_ext
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = es; ErrMsg = 'solve_consistent_acceleration: '//TRIM(em2); RETURN
    END IF
    ! Total mass = structural + added (the added-mass matrix at the current configuration).
    model%ws_Msumb = model%Mgb
    IF (model%has_am) THEN
      CALL assemble_added_mass(model, model%q, es, em2)
      IF (es /= CD_HCDYN_OK) THEN
        ErrStat = es; ErrMsg = 'solve_consistent_acceleration: '//TRIM(em2); RETURN
      END IF
      model%ws_Msumb = model%ws_Msumb + model%ws_Mamb
    END IF
    ! Build the fixed-DOF acceleration vector before the early return: a fully constrained
    ! model still has to retain a caller-prescribed boundary acceleration.
    model%ws_an = CD_ZERO
    IF (has_pres) model%ws_an(pres_dofs) = pres_a
    CALL build_rigid_bases(model, model%endconn_d0, rigid_basis, es, em2)
    IF (es /= CD_HCDYN_OK) THEN
      ErrStat = es; ErrMsg = 'solve_consistent_acceleration: '//TRIM(em2); RETURN
    END IF
    DO iend = 1, 2
      IF (model%endconn_mode(iend) /= CD_ENDCONN_RIGID) CYCLE
      IF (iend == 1) THEN
        base = 3
      ELSE
        base = 6*(model%nn - 1) + 3
      END IF
      old_acceleration = model%ws_vk(base + 1:base + 3)
      axial_acceleration = DOT_PRODUCT(old_acceleration, rigid_basis(:, 1, iend))
      model%ws_an(base + 1:base + 3) = old_acceleration - &
                                       axial_acceleration*rigid_basis(:, 1, iend)
    END DO
    IF (nf == 0) THEN
      model%a = model%ws_an
      RETURN
    END IF
    ! In-band Dirichlet on a mass copy + full-size banded solve. The free RHS is
    ! f_ext-f_int-M*a_prescribed; fixed RHS slots remain zero for the Dirichlet solve.
    model%ws_Keffb = model%ws_Msumb
    CALL band_matvec(model%ws_Msumb, model%ws_an, model%ws_inert)
    model%ws_dq = -model%ws_R - model%ws_inert
    CALL transform_rigid_vector_to_local(model, rigid_basis, model%ws_dq)
    CALL transform_rigid_band(model, rigid_basis, model%ws_Keffb)
    CALL dirichlet_band(model%ws_Keffb, model%solve_mask)
    DO i = 1, model%ndof
      IF (.NOT. model%solve_mask(i)) model%ws_dq(i) = CD_ZERO
    END DO
    CALL CD_Factor_Banded(model%ws_Keffb, KL_D, KU_D, model%ws_ipiv, es, em2)
    IF (es == CD_LINALG_OK) THEN
      CALL CD_Solve_Factored_Banded(model%ws_Keffb, KL_D, KU_D, model%ws_ipiv, model%ws_dq, es, em2)
    END IF
    IF (es /= CD_LINALG_OK) THEN
      ErrStat = CD_HCDYN_NOCONVERGE; ErrMsg = 'solve_consistent_acceleration: solve failed: '//TRIM(em2); RETURN
    END IF
    CALL transform_rigid_vector_to_global(model, rigid_basis, model%ws_dq)
    model%a = model%ws_dq + model%ws_an
    IF (has_pres) model%a(pres_dofs) = pres_a
  END SUBROUTINE solve_consistent_acceleration

  SUBROUTINE scatter12_band(abg, Ke, e)
    !! Scatter a 12x12 element matrix into LAPACK general-band storage for element e
    !! (entry (gi, gj) at abg(KL_D + KU_D + 1 + gi - gj, gj); the chain connectivity
    !! keeps |gi - gj| <= 11 = KU_D, so every entry is in band by construction).
    REAL(wp), INTENT(INOUT) :: abg(:, :)
    REAL(wp), INTENT(IN) :: Ke(12, 12)
    INTEGER, INTENT(IN) :: e
    INTEGER :: a, b, gi, gj
    DO a = 1, 12
      gi = elem_gdof(e, a)
      DO b = 1, 12
        gj = elem_gdof(e, b)
        abg(KL_D + KU_D + 1 + gi - gj, gj) = abg(KL_D + KU_D + 1 + gi - gj, gj) + Ke(a, b)
      END DO
    END DO
  END SUBROUTINE scatter12_band

  PURE SUBROUTINE band_matvec(abg, x, y)
    !! y = A x for A in the module's LAPACK general-band layout.
    REAL(wp), INTENT(IN) :: abg(:, :), x(:)
    REAL(wp), INTENT(OUT) :: y(:)
    INTEGER :: i, j, n
    n = SIZE(x)
    y = CD_ZERO
    DO j = 1, n
      DO i = MAX(1, j - KU_D), MIN(n, j + KL_D)
        y(i) = y(i) + abg(KL_D + KU_D + 1 + i - j, j)*x(j)
      END DO
    END DO
  END SUBROUTINE band_matvec

  PURE SUBROUTINE dirichlet_band(abg, freemask)
    !! Apply Dirichlet conditions in band form: zero each fixed DOF's row and column
    !! within the band and set a unit diagonal, so the full-size banded solve returns
    !! exactly zero update at every fixed DOF (with a zero RHS slot).
    REAL(wp), INTENT(INOUT) :: abg(:, :)
    LOGICAL, INTENT(IN) :: freemask(:)
    INTEGER :: j, gi, gj, n
    n = SIZE(freemask)
    DO j = 1, n
      IF (freemask(j)) CYCLE
      DO gi = MAX(1, j - KL_D), MIN(n, j + KU_D)          ! column j entries (rows gi)
        abg(KL_D + KU_D + 1 + gi - j, j) = CD_ZERO
      END DO
      DO gj = MAX(1, j - KU_D), MIN(n, j + KL_D)          ! row j entries (columns gj)
        abg(KL_D + KU_D + 1 + j - gj, gj) = CD_ZERO
      END DO
      abg(KL_D + KU_D + 1, j) = CD_ONE
    END DO
  END SUBROUTINE dirichlet_band

  SUBROUTINE gather_elem(model, e, qe)
    !! Gather element e's 12 DOFs from the model state.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    INTEGER, INTENT(IN) :: e
    REAL(wp), INTENT(OUT) :: qe(12)
    qe(1:3) = model%q(6*(e - 1) + 1:6*(e - 1) + 3)
    qe(4:6) = model%q(6*(e - 1) + 4:6*(e - 1) + 6)
    qe(7:9) = model%q(6*e + 1:6*e + 3)
    qe(10:12) = model%q(6*e + 4:6*e + 6)
  END SUBROUTINE gather_elem

  PURE INTEGER FUNCTION elem_gdof(e, a) RESULT(g)
    !! Map element-local DOF a (1..12) of element e to the global DOF index.
    INTEGER, INTENT(IN) :: e, a
    IF (a <= 6) THEN
      g = 6*(e - 1) + a
    ELSE
      g = 6*e + (a - 6)
    END IF
  END FUNCTION elem_gdof

END MODULE CableDyn_HermiteCableDynamic
