! File: tests/test_hermite_torsion_dynamic.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_hermite_torsion_dynamic
  !! Dynamic gates of condensed (Route A) torsion on the cubic-Hermite cable, on a 10 m line
  !! clamped (Rigid bending connections) at both ends with both ends restrained in torsion:
  !!   S  step response: an imposed-twist step (u_twist) on a straight line gives the quasi-static
  !!      torque GJ Phi / L at the very first step and every later one, with the centreline
  !!      unchanged. Route A has no torsional inertia, so the torque follows the imposed twist
  !!      without a torsional wave or overshoot: the dynamic and quasi-static answers coincide;
  !!   T  turning parent: the coupled end's torsion frame turns with the parent orientation about
  !!      the line axis through 1.5 turns while the imposed twist follows it, so Theta passes
  !!      several 2 pi branches (unwrapped, no slip) and the torque stays constant; turning the
  !!      parent alone changes the torque by (roll)/C;
  !!   E  energy: free vibration of a bowed, twisted line (rho_inf = 1, no drag): kinetic +
  !!      axial/bending + torsional energy (the torsional share exchanges with the bending)
  !!      stays within 0.5 % of the vibration energy over ten periods;
  !!   R  restart: a moving, rolling run restored from a snapshot, and one rebuilt from the
  !!      checkpoint mirror into a fresh module, continue bit-identically (Theta and torque
  !!      included) after the frame has turned through more than one turn;
  !!   B  the bordered (Sherman-Morrison) step converges as fast as the line without torsion
  !!      on the same swaying, rolling motion with bending-pinned torsional ends (in such steps
  !!      the rank-one term is small beside the inertia; it carries weight near buckling), and
  !!      a model with torsion refuses the configuration blend.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_EndConnection, ONLY: CD_ENDCONN_RIGID
  USE CableDyn_HermiteTorsion, ONLY: CD_HermiteTorsionType, CD_HermiteTorsion_Line, CD_HermiteTorsion_Unwrap, &
                                     CD_HTORS_OK
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_Energy, &
                                          CD_HermiteCable_Dyn_End, CD_HermiteCable_Dyn_Set_EndConnection, &
                                          CD_HermiteCable_Dyn_Set_Torsion, CD_HermiteCable_Dyn_Set_TangentReuse, &
                                          CD_HermiteCable_Dyn_Set_ForceBlend, CD_HermiteCable_Dyn_Torsion_State, &
                                          CD_HermiteCable_Dyn_Reset_Profile, CD_HermiteCable_Dyn_Get_Profile, &
                                          CD_HermiteCable_Dyn_Disable_Profile, CD_HCDYN_OK
  USE CableDyn_OpenFAST_HermiteFMF, ONLY: CD_HFMF_ModuleType, CD_HFMF_Init, CD_HFMF_Set_EndConnection, &
                                          CD_HFMF_Set_Torsion, CD_HFMF_UpdateStates, CD_HFMF_CalcOutput, &
                                          CD_HFMF_Snapshot, CD_HFMF_Restore, CD_HFMF_MirrorSize, &
                                          CD_HFMF_PackMirror, CD_HFMF_UnpackMirror, CD_HFMF_SetCoupledKinematics, &
                                          CD_HFMF_End, CD_HFMF_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
  REAL(wp), PARAMETER :: LL = 10.0_wp, EA = 1.0e5_wp, EI = 100.0_wp, GJ = 80.0_wp, RHOA = 1.0_wp
  INTEGER, PARAMETER :: NE = 20, NN = NE + 1, NDOF = 6*NN
  INTEGER :: nfail = 0

  CALL check_step_response()
  CALL check_turning_parent()
  CALL check_energy()
  CALL check_restart()
  CALL check_bordered_and_blend()

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: Hermite torsion dynamics (step, turning parent, energy, restart, bordered step)'

CONTAINS

  INCLUDE 'nan_max_abs.inc'

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. cond) THEN
      nfail = nfail + 1
      WRITE (*, '(A)') 'FAIL: '//label
    END IF
  END SUBROUTINE require

  FUNCTION roll_x(psi) RESULT(r)
    !! Global-to-parent DCM of a parent rolled by psi about +x (the transpose of the rotation).
    REAL(wp), INTENT(IN) :: psi
    REAL(wp) :: r(3, 3)
    r = 0.0_wp
    r(1, 1) = 1.0_wp
    r(2, 2) = COS(psi)
    r(3, 3) = COS(psi)
    r(2, 3) = SIN(psi)
    r(3, 2) = -SIN(psi)
  END FUNCTION roll_x

  SUBROUTINE line_seed(bow_y, bow_z, q)
    !! Straight line along +x (node 1 at the origin), optionally bowed by y = bow_y sin^2(pi s/L),
    !! z = bow_z sin^2(2 pi s/L) (end tangents stay along +x, as the clamps require).
    REAL(wp), INTENT(IN) :: bow_y, bow_z
    REAL(wp), INTENT(OUT) :: q(NDOF)
    INTEGER :: k
    REAL(wp) :: s
    q = 0.0_wp
    DO k = 1, NN
      s = LL*REAL(k - 1, wp)/REAL(NE, wp)
      q(6*k - 5) = s
      q(6*k - 4) = bow_y*SIN(PI*s/LL)**2
      q(6*k - 3) = bow_z*SIN(2.0_wp*PI*s/LL)**2
      q(6*k - 2) = 1.0_wp
      q(6*k - 1) = bow_y*(PI/LL)*SIN(2.0_wp*PI*s/LL)
      q(6*k) = bow_z*(2.0_wp*PI/LL)*SIN(4.0_wp*PI*s/LL)
    END DO
  END SUBROUTINE line_seed

  SUBROUTINE torsion_of(q, phi, tors, gj_line)
    !! Condensed torsion of the line: both ends clamped along +x with reference normal +z, GJ on
    !! every element, imposed twist phi, Theta of the seed q (branch nearest 0).
    REAL(wp), INTENT(IN) :: q(NDOF), phi
    TYPE(CD_HermiteTorsionType), INTENT(OUT) :: tors
    REAL(wp), INTENT(IN), OPTIONAL :: gj_line
    REAL(wp) :: g(NDOF), raw, l0(NE)
    INTEGER :: es
    CHARACTER(200) :: em
    l0 = LL/REAL(NE, wp)
    tors%active = .TRUE.
    tors%phi = phi
    tors%ends = RESHAPE([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
                         0.0_wp, 0.0_wp, 1.0_wp], [3, 4])
    ALLOCATE (tors%gj(NE))
    tors%gj = GJ
    IF (PRESENT(gj_line)) tors%gj = gj_line
    CALL CD_HermiteTorsion_Line(q, l0, tors%ends, raw, g, es, em)
    CALL require(es == CD_HTORS_OK, 'seed twist evaluates: '//TRIM(em))
    tors%theta = CD_HermiteTorsion_Unwrap(raw, 0.0_wp)
    tors%has_theta = .TRUE.
  END SUBROUTINE torsion_of

  SUBROUTINE build_cable(cab, q, phi, dt, tol, rho_inf, pinned, gj_line)
    !! HFMF cable on seed q: node NN coupled (driven), node 1 held; Rigid bending connections
    !! along +x at both ends (bending-pinned with pinned), the coupled one in parent axes
    !! (identity attitude); torsion with the coupled end's frame attached to the parent.
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: cab
    REAL(wp), INTENT(IN) :: q(NDOF), phi, dt, tol, rho_inf
    LOGICAL, INTENT(IN), OPTIONAL :: pinned
    REAL(wp), INTENT(IN), OPTIONAL :: gj_line
    TYPE(CD_HermiteTorsionType) :: tors
    REAL(wp) :: l0(NE), eye(3, 3), frame(3, 2)
    INTEGER :: es
    LOGICAL :: rigid
    CHARACTER(300) :: em
    l0 = LL/REAL(NE, wp)
    eye = roll_x(0.0_wp)
    rigid = .TRUE.
    IF (PRESENT(pinned)) rigid = .NOT. pinned
    CALL CD_HFMF_Init(cab, l0, [(EA, es=1, NE)], [(EI, es=1, NE)], [(RHOA, es=1, NE)], [(0.0_wp, es=1, NE)], q, &
                      [1, 2, 3, NDOF - 5, NDOF - 4, NDOF - 3], 0.0_wp, 0.0_wp, rho_inf, dt, NN, es, em, &
                      max_iter=60, tol=tol)
    CALL require(es == CD_HFMF_OK, 'cable init: '//TRIM(em))
    IF (rigid) THEN
      CALL CD_HFMF_Set_EndConnection(cab, [0.0_wp, 0.0_wp], RESHAPE([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], &
                                                                    [3, 2]), es, em, &
                                     coupled_d0_parent=[1.0_wp, 0.0_wp, 0.0_wp], parent_orientation=eye, &
                                     connection_mode=[CD_ENDCONN_RIGID, CD_ENDCONN_RIGID])
      CALL require(es == CD_HFMF_OK, 'cable end connections: '//TRIM(em))
    END IF
    CALL torsion_of(q, phi, tors, gj_line)
    frame(:, 1) = [1.0_wp, 0.0_wp, 0.0_wp]
    frame(:, 2) = [0.0_wp, 0.0_wp, 1.0_wp]
    CALL CD_HFMF_Set_Torsion(cab, tors, es, em, coupled_frame_parent=frame, parent_orientation=eye)
    CALL require(es == CD_HFMF_OK, 'cable torsion: '//TRIM(em))
  END SUBROUTINE build_cable

  ! ------------------------------------------------------------------------------------------
  SUBROUTINE check_step_response()
    TYPE(CD_HFMF_ModuleType) :: cab
    REAL(wp) :: q(NDOF), xa(3), phi, worst, bent, th, mt, f(3), m(3), mref
    INTEGER :: k, es
    CHARACTER(300) :: em
    CALL line_seed(0.0_wp, 0.0_wp, q)
    CALL build_cable(cab, q, 1.0_wp, 0.01_wp, 1.0e-8_wp, 0.8_wp)
    xa = q(NDOF - 5:NDOF - 3)
    worst = 0.0_wp
    bent = 0.0_wp
    DO k = 1, 40
      phi = 1.0_wp
      IF (k >= 3) phi = 3.0_wp                  ! step of 2 rad at t = 0.03 s
      CALL CD_HFMF_UpdateStates(cab, xa, [0.0_wp, 0.0_wp, 0.0_wp], [0.0_wp, 0.0_wp, 0.0_wp], es, em, u_twist=phi)
      CALL require(es == CD_HFMF_OK, 'step-response step: '//TRIM(em))
      IF (es /= CD_HFMF_OK) EXIT
      CALL CD_HermiteCable_Dyn_Torsion_State(cab%line, th, mt, es, em)
      mref = GJ*phi/LL
      worst = MAX(worst, ABS(mt - mref)/mref)
      bent = MAX(bent, MAXVAL(ABS(cab%line%q(2:NDOF:6))), MAXVAL(ABS(cab%line%q(3:NDOF:6))))
      IF (k == 3) THEN
        CALL CD_HFMF_CalcOutput(cab, f, es, em, y_moment=m)
        CALL require(es == CD_HFMF_OK .AND. ABS(m(1) + mref) <= 1.0e-9_wp*mref .AND. &
                     NORM2(m(2:3)) <= 1.0e-9_wp*mref, 'S: the coupled clamp carries -M_t along the axis')
      END IF
    END DO
    WRITE (*, '(A,ES10.3,A,ES10.3,A)') 'S step response: torque vs GJ Phi/L ', worst, &
      ' (every step, quasi-static from the first), transverse motion ', bent, ' m'
    CALL require(worst <= 1.0e-9_wp, 'S: torque equals the quasi-static GJ Phi / L at every step')
    CALL require(bent <= 1.0e-12_wp, 'S: a straight twisted line below the buckling torque stays straight')
    CALL CD_HFMF_End(cab)
  END SUBROUTINE check_step_response

  ! ------------------------------------------------------------------------------------------
  SUBROUTINE check_turning_parent()
    TYPE(CD_HFMF_ModuleType) :: cab
    REAL(wp) :: q(NDOF), xa(3), psi, th, mt, th0, err_th, err_mt
    INTEGER :: k, es
    INTEGER, PARAMETER :: NSTEP = 120
    CHARACTER(300) :: em
    CALL line_seed(0.0_wp, 0.0_wp, q)
    CALL build_cable(cab, q, 2.0_wp, 0.01_wp, 1.0e-8_wp, 0.8_wp)
    xa = q(NDOF - 5:NDOF - 3)
    th0 = cab%line%torsion%theta
    err_th = 0.0_wp
    err_mt = 0.0_wp
    ! 1.5 turns of the parent about the line axis; the imposed twist follows the roll
    DO k = 1, NSTEP
      psi = 3.0_wp*PI*REAL(k, wp)/REAL(NSTEP, wp)
      CALL CD_HFMF_UpdateStates(cab, xa, [0.0_wp, 0.0_wp, 0.0_wp], [0.0_wp, 0.0_wp, 0.0_wp], es, em, &
                                u_orientation=roll_x(psi), u_twist=2.0_wp - psi)
      CALL require(es == CD_HFMF_OK, 'T: turning-parent step: '//TRIM(em))
      IF (es /= CD_HFMF_OK) EXIT
      CALL CD_HermiteCable_Dyn_Torsion_State(cab%line, th, mt, es, em)
      err_th = MAX(err_th, ABS(th - (th0 - psi)))
      err_mt = MAX(err_mt, ABS(mt - GJ*2.0_wp/LL))
    END DO
    WRITE (*, '(A,F9.4,A,ES10.3,A,ES10.3)') 'T turning parent: Theta after 1.5 turns ', cab%line%torsion%theta, &
      ' rad, Theta error ', err_th, ', torque error ', err_mt
    CALL require(err_th <= 1.0e-10_wp, 'T: Theta follows the parent roll through 1.5 turns (unwrapped)')
    CALL require(err_mt <= 1.0e-9_wp, 'T: the torque stays GJ Phi_0 / L when the twist follows the roll')
    ! the parent alone turns on by 0.3 rad: the torque grows by 0.3/C
    CALL CD_HFMF_UpdateStates(cab, xa, [0.0_wp, 0.0_wp, 0.0_wp], [0.0_wp, 0.0_wp, 0.0_wp], es, em, &
                              u_orientation=roll_x(3.0_wp*PI + 0.3_wp))
    CALL CD_HermiteCable_Dyn_Torsion_State(cab%line, th, mt, es, em)
    WRITE (*, '(A,ES14.6,A,ES14.6)') 'T parent roll alone: torque ', mt, ' expected ', GJ*(2.0_wp + 0.3_wp)/LL
    CALL require(es == CD_HFMF_OK .AND. ABS(mt - GJ*(2.0_wp + 0.3_wp)/LL) <= 1.0e-9_wp, &
                 'T: a parent roll alone changes the torque by roll/C (the frame turn adds to the imposed twist)')
    CALL CD_HFMF_End(cab)
  END SUBROUTINE check_turning_parent

  ! ------------------------------------------------------------------------------------------
  SUBROUTINE energy_run(phi, band, e_vib, t_exchange)
    !! Free vibration of the bowed line clamped at both ends (no gravity, no fluid), rho_inf = 1,
    !! ten periods of the first bending mode; band = (max - min) of kinetic + elastic (with the
    !! torsional share) over the run, relative to the vibration energy above the straight
    !! twisted state; t_exchange = the swing of the torsional energy alone, same scale.
    REAL(wp), INTENT(IN) :: phi
    REAL(wp), INTENT(OUT) :: band, e_vib, t_exchange
    TYPE(CD_HermiteCableDynType) :: m
    TYPE(CD_HermiteTorsionType) :: tors
    REAL(wp) :: q(NDOF), l0(NE), ke, se, e0, emax, emin, t1, dt, th, mt, et, etmin, etmax, e_straight
    INTEGER :: k, es, nstep
    CHARACTER(300) :: em
    l0 = LL/REAL(NE, wp)
    CALL line_seed(0.05_wp, 0.03_wp, q)
    CALL CD_HermiteCable_Dyn_Init(m, l0, [(EA, k=1, NE)], [(EI, k=1, NE)], [(RHOA, k=1, NE)], [(0.0_wp, k=1, NE)], q, &
                                  [1, 2, 3, NDOF - 5, NDOF - 4, NDOF - 3], 0.0_wp, 0.0_wp, 1.0_wp, es, em)
    CALL require(es == CD_HCDYN_OK, 'E: init: '//TRIM(em))
    CALL CD_HermiteCable_Dyn_Set_EndConnection(m, [0.0_wp, 0.0_wp], &
                                               RESHAPE([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]), &
                                               es, em, connection_mode=[CD_ENDCONN_RIGID, CD_ENDCONN_RIGID])
    CALL require(es == CD_HCDYN_OK, 'E: end connections: '//TRIM(em))
    IF (phi > 0.0_wp) THEN
      CALL torsion_of(q, phi, tors)
      CALL CD_HermiteCable_Dyn_Set_Torsion(m, tors, es, em)
      CALL require(es == CD_HCDYN_OK, 'E: torsion: '//TRIM(em))
    END IF
    t1 = 2.0_wp*PI/(22.373_wp*SQRT(EI/(RHOA*LL**4)))
    dt = t1/400.0_wp
    nstep = 4000
    e_straight = 0.5_wp*(LL/GJ)*(GJ*phi/LL)**2
    CALL CD_HermiteCable_Dyn_Energy(m, ke, se, es, em)
    e0 = ke + se
    e_vib = e0 - e_straight
    emax = e0
    emin = e0
    etmax = -HUGE(1.0_wp)
    etmin = HUGE(1.0_wp)
    DO k = 1, nstep
      CALL CD_HermiteCable_Dyn_Step(m, dt, 40, 1.0e-10_wp, es, em)
      IF (es /= CD_HCDYN_OK) THEN
        CALL require(.FALSE., 'E: step: '//TRIM(em))
        EXIT
      END IF
      CALL CD_HermiteCable_Dyn_Energy(m, ke, se, es, em)
      emax = MAX(emax, ke + se)
      emin = MIN(emin, ke + se)
      IF (phi > 0.0_wp) THEN
        CALL CD_HermiteCable_Dyn_Torsion_State(m, th, mt, es, em)
        et = 0.5_wp*(LL/GJ)*mt*mt
        etmax = MAX(etmax, et)
        etmin = MIN(etmin, et)
      END IF
    END DO
    band = (emax - emin)/e_vib
    t_exchange = 0.0_wp
    IF (phi > 0.0_wp) t_exchange = (etmax - etmin)/e_vib
    CALL CD_HermiteCable_Dyn_End(m)
  END SUBROUTINE energy_run

  SUBROUTINE check_energy()
    REAL(wp) :: band0, band1, ev0, ev1, tx0, tx1
    CALL energy_run(0.0_wp, band0, ev0, tx0)
    CALL energy_run(4.0_wp, band1, ev1, tx1)
    WRITE (*, '(A,ES10.3,A,ES10.3,A,ES10.3,A)') 'E energy band over 10 periods (rho_inf = 1): twisted ', band1, &
      ', untwisted ', band0, ' of the vibration energy; torsional share swings by ', tx1, ' of it'
    CALL require(band1 <= 5.0e-3_wp, 'E: kinetic + elastic + torsional energy conserved within 0.5 % (rho_inf = 1)')
    CALL require(tx1 >= 1.0e-4_wp, 'E: the torsional energy exchanges with the bending (coupling active)')
  END SUBROUTINE check_energy

  ! ------------------------------------------------------------------------------------------
  SUBROUTINE drive(k, x, v, a, dcm, w, al, phi)
    !! Coupled-end drive of the restart runs: a lateral sway with a roll of the parent about the
    !! line axis (3 rad/s, more than a turn by step 250) and an imposed twist following it.
    INTEGER, INTENT(IN) :: k
    REAL(wp), INTENT(OUT) :: x(3), v(3), a(3), dcm(3, 3), w(3), al(3), phi
    REAL(wp), PARAMETER :: DT = 0.01_wp, AMP = 0.05_wp, OM = 2.0_wp, ROLL_RATE = 3.0_wp
    REAL(wp) :: t
    t = DT*REAL(k, wp)
    x = [LL, AMP*SIN(OM*t), 0.0_wp]
    v = [0.0_wp, AMP*OM*COS(OM*t), 0.0_wp]
    a = [0.0_wp, -AMP*OM*OM*SIN(OM*t), 0.0_wp]
    dcm = roll_x(ROLL_RATE*t)
    w = [ROLL_RATE, 0.0_wp, 0.0_wp]
    al = 0.0_wp
    phi = 1.0_wp - ROLL_RATE*t + 0.5_wp*SIN(OM*t)
  END SUBROUTINE drive

  SUBROUTINE march(cab, k0, k1, ok)
    TYPE(CD_HFMF_ModuleType), INTENT(INOUT) :: cab
    INTEGER, INTENT(IN) :: k0, k1
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp) :: x(3), v(3), a(3), dcm(3, 3), w(3), al(3), phi
    INTEGER :: k, es
    CHARACTER(300) :: em
    ok = .TRUE.
    DO k = k0 + 1, k1
      CALL drive(k, x, v, a, dcm, w, al, phi)
      CALL CD_HFMF_UpdateStates(cab, x, v, a, es, em, u_orientation=dcm, u_angular_velocity=w, &
                                u_angular_acceleration=al, u_twist=phi)
      IF (es /= CD_HFMF_OK) THEN
        WRITE (*, '(A,I0,A)') '  restart march step ', k, ': '//TRIM(em)
        ok = .FALSE.
        RETURN
      END IF
    END DO
  END SUBROUTINE march

  SUBROUTINE check_restart()
    INTEGER, PARAMETER :: K_SNAP = 250, K_END = 400
    TYPE(CD_HFMF_ModuleType) :: cab, fresh
    REAL(wp) :: q(NDOF), qa(NDOF), va(NDOF), tha, mta, th, mt, x(3), v(3), a(3), dcm(3, 3), w(3), al(3), phi
    REAL(wp) :: t_snap
    REAL(wp), ALLOCATABLE :: buf(:)
    LOGICAL :: ok
    INTEGER :: es
    CHARACTER(300) :: em
    CALL line_seed(0.0_wp, 0.0_wp, q)
    CALL build_cable(cab, q, 1.0_wp, 0.01_wp, 1.0e-6_wp, 0.8_wp)
    ! a restart is exact with cross-step factor reuse off (the factor's age is not state)
    CALL CD_HermiteCable_Dyn_Set_TangentReuse(cab%line, .FALSE., es, em)
    CALL march(cab, 0, K_SNAP, ok)
    CALL require(ok, 'R: march to the snapshot')
    CALL require(ABS(cab%line%torsion%theta) > 2.0_wp*PI, 'R: Theta has passed more than one turn at the snapshot')
    t_snap = cab%line%t
    CALL CD_HFMF_Snapshot(cab, es, em)
    CALL require(es == CD_HFMF_OK, 'R: snapshot: '//TRIM(em))
    ALLOCATE (buf(CD_HFMF_MirrorSize(cab)))
    CALL CD_HFMF_PackMirror(cab, buf, es, em)
    CALL require(es == CD_HFMF_OK, 'R: pack mirror: '//TRIM(em))
    CALL march(cab, K_SNAP, K_END, ok)
    CALL require(ok, 'R: uninterrupted march')
    qa = cab%line%q
    va = cab%line%v
    CALL CD_HermiteCable_Dyn_Torsion_State(cab%line, tha, mta, es, em)
    ! (1) restore the snapshot and re-run
    CALL CD_HFMF_Restore(cab, es, em)
    CALL require(es == CD_HFMF_OK, 'R: restore: '//TRIM(em))
    CALL march(cab, K_SNAP, K_END, ok)
    CALL CD_HermiteCable_Dyn_Torsion_State(cab%line, th, mt, es, em)
    CALL require(ok .AND. nan_max_abs(cab%line%q - qa) <= 0.0_wp .AND. nan_max_abs(cab%line%v - va) <= 0.0_wp .AND. &
                 .NOT. (ABS(th - tha) > 0.0_wp) .AND. .NOT. (ABS(mt - mta) > 0.0_wp), &
                 'R: snapshot/restore continues bit-identically (state, Theta, torque)')
    ! (2) a fresh module from the same static state, reloaded from the mirror
    CALL build_cable(fresh, q, 1.0_wp, 0.01_wp, 1.0e-6_wp, 0.8_wp)
    CALL CD_HermiteCable_Dyn_Set_TangentReuse(fresh%line, .FALSE., es, em)
    CALL drive(K_SNAP, x, v, a, dcm, w, al, phi)
    CALL CD_HFMF_UnpackMirror(fresh, buf, t_snap, es, em)
    CALL require(es == CD_HFMF_OK, 'R: unpack mirror: '//TRIM(em))
    CALL CD_HFMF_SetCoupledKinematics(fresh, x, v, a, es, em, u_orientation=dcm, u_angular_velocity=w, &
                                      u_angular_acceleration=al)
    CALL require(es == CD_HFMF_OK, 'R: re-prescribe the coupled end after the reload: '//TRIM(em))
    CALL march(fresh, K_SNAP, K_END, ok)
    CALL CD_HermiteCable_Dyn_Torsion_State(fresh%line, th, mt, es, em)
    WRITE (*, '(A,F10.4,A,ES10.3,A,ES10.3)') 'R restart: Theta ', th, ' rad, |dq| (mirror) ', &
      nan_max_abs(fresh%line%q - qa), ', |dTheta| ', ABS(th - tha)
    CALL require(ok .AND. nan_max_abs(fresh%line%q - qa) <= 0.0_wp .AND. nan_max_abs(fresh%line%v - va) <= 0.0_wp &
                 .AND. .NOT. (ABS(th - tha) > 0.0_wp), &
                 'R: a checkpoint-mirror restart continues bit-identically (state, Theta)')
    CALL CD_HFMF_End(cab)
    CALL CD_HFMF_End(fresh)
  END SUBROUTINE check_restart

  ! ------------------------------------------------------------------------------------------
  SUBROUTINE check_bordered_and_blend()
    !! The same swaying, rolling drive on a bowed line with and without torsion (a twist of about
    !! half the buckling twist): the bordered step keeps the Newton solve count of the untwisted
    !! line (a missing rank-one term would cost iterations), and torsion refuses the
    !! configuration blend.
    TYPE(CD_HFMF_ModuleType) :: cab
    TYPE(CD_HermiteCableDynType) :: m
    TYPE(CD_HermiteTorsionType) :: tors
    REAL(wp) :: q(NDOF), x(3), v(3), a(3), dcm(3, 3), w(3), al(3), phi, l0(NE), tt(6)
    INTEGER :: k, es, pass, nst, nres, nsol(2), nam, nal, ntan
    CHARACTER(300) :: em
    CALL line_seed(0.02_wp, 0.01_wp, q)
    l0 = LL/REAL(NE, wp)
    DO pass = 1, 2
      CALL build_cable(cab, q, 3.0_wp, 0.01_wp, 1.0e-8_wp, 0.8_wp, pinned=.TRUE.)
      IF (pass == 1) THEN
        tors%active = .FALSE.
        CALL CD_HermiteCable_Dyn_Set_Torsion(cab%line, tors, es, em)
      END IF
      CALL CD_HermiteCable_Dyn_Set_TangentReuse(cab%line, .FALSE., es, em)
      CALL CD_HermiteCable_Dyn_Reset_Profile()
      DO k = 1, 80
        CALL drive(k, x, v, a, dcm, w, al, phi)
        x(2) = 10.0_wp*x(2)
        v(2) = 10.0_wp*v(2)
        a(2) = 10.0_wp*a(2)
        IF (pass == 1) THEN
          CALL CD_HFMF_UpdateStates(cab, x, v, a, es, em, u_orientation=dcm, u_angular_velocity=w, &
                                    u_angular_acceleration=al)
        ELSE
          CALL CD_HFMF_UpdateStates(cab, x, v, a, es, em, u_orientation=dcm, u_angular_velocity=w, &
                                    u_angular_acceleration=al, u_twist=phi + 2.0_wp)
        END IF
        CALL require(es == CD_HFMF_OK, 'B: step: '//TRIM(em))
        IF (es /= CD_HFMF_OK) EXIT
      END DO
      CALL CD_HermiteCable_Dyn_Get_Profile(nst, nres, nsol(pass), nam, nal, tt(1), tt(2), tt(3), tt(4), tt(5), &
                                           ntan, tt(6))
      CALL CD_HermiteCable_Dyn_Disable_Profile()
      CALL CD_HFMF_End(cab)
    END DO
    WRITE (*, '(A,I0,A,I0,A)') 'B bordered step: ', nsol(2), ' Newton solves with torsion, ', nsol(1), &
      ' without (80 steps)'
    CALL require(nsol(2) <= nsol(1) + 80, 'B: the bordered step needs at most one more solve per step')
    ! the configuration blend is refused by name
    CALL CD_HermiteCable_Dyn_Init(m, l0, [(EA, k=1, NE)], [(EI, k=1, NE)], [(RHOA, k=1, NE)], [(0.0_wp, k=1, NE)], q, &
                                  [1, 2, 3, NDOF - 5, NDOF - 4, NDOF - 3], 0.0_wp, 0.0_wp, 0.8_wp, es, em)
    CALL CD_HermiteCable_Dyn_Set_EndConnection(m, [0.0_wp, 0.0_wp], &
                                               RESHAPE([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp], [3, 2]), &
                                               es, em, connection_mode=[CD_ENDCONN_RIGID, CD_ENDCONN_RIGID])
    CALL torsion_of(q, 1.0_wp, tors)
    CALL CD_HermiteCable_Dyn_Set_Torsion(m, tors, es, em)
    CALL CD_HermiteCable_Dyn_Set_ForceBlend(m, .FALSE., es, em)
    CALL CD_HermiteCable_Dyn_Step(m, 0.01_wp, 30, 1.0e-8_wp, es, em)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'force-blended') > 0, &
                 'B: torsion refuses the configuration blend by name (got: '//TRIM(em)//')')
    CALL CD_HermiteCable_Dyn_End(m)
  END SUBROUTINE check_bordered_and_blend

END PROGRAM test_hermite_torsion_dynamic
