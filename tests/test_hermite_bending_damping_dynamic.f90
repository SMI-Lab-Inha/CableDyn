! File: tests/test_hermite_bending_damping_dynamic.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_hermite_bending_damping_dynamic
  !! Free decay of the first mode of a weightless, pinned-pinned, tensioned finite-EI cable with
  !! bending Kelvin-Voigt damping (CD_HermiteCable_Dyn_Set_Bending_Damping), generalised-alpha with
  !! rho_inf = 1 (no numerical damping). For mode shape sin(k x), k = pi/L, the bending moment is
  !! EI k^2 a and the tension carries T k^2 a, so the mode's damping ratio is
  !!   zeta_1 = omega_1 Bb / (2 EI) * EI k^2 / (EI k^2 + T),   omega_1^2 = (EI k^4 + T k^2)/m:
  !! only the bending share of the modal stiffness is damped.
  !!   B1 a bending-dominated cable (share 0.96) decays at that zeta to 1 %;
  !!   B2 a tension-dominated cable (share 0.2) decays at that smaller zeta to 1 %;
  !!   B3 Bb = 0 leaves the undamped decay (|zeta| below 1e-5) and installing an all-zero array
  !!      keeps the configuration off;
  !!   B4 a negative coefficient or a wrong length is refused;
  !!   B5 a straight cable moving in its first bending mode with zero acceleration has no elastic
  !!      curvature: the node recovery reports none although the bending damping acts (its viscous
  !!      moment is not EI kappa), while the end force carries the damping.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_Set_Bending_Damping, &
                                          CD_HermiteCable_Dyn_Node_Recovery, CD_HermiteCable_Dyn_End_Force, &
                                          CD_HCDYN_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
  REAL(wp), PARAMETER :: LL = 20.0_wp, EA = 1.0e8_wp, EI = 1.0e5_wp, RHOA = 10.0_wp
  INTEGER, PARAMETER :: NE = 40, NN = NE + 1, NDOF = 6*NN, MID = NE/2 + 1
  INTEGER :: nfail = 0

  CALL decay_case('B1 bending-dominated', 1.0e2_wp, 0.02_wp)
  CALL decay_case('B2 tension-dominated', 1.0e4_wp, 0.02_wp)
  CALL decay_case('B3 undamped', 1.0e4_wp, 0.0_wp)
  CALL check_refusals()
  CALL check_recovery()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: Hermite bending damping dynamic decay'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE setup(model, tens, a0)
    TYPE(CD_HermiteCableDynType), INTENT(OUT) :: model
    REAL(wp), INTENT(IN) :: tens, a0
    REAL(wp) :: l0(NE), ea_(NE), ei_(NE), rho(NE), w(NE), q(NDOF), eps, lt, x
    INTEGER :: k, es
    CHARACTER(512) :: em
    eps = tens/EA; lt = (1.0_wp + eps)*LL
    l0 = LL/REAL(NE, wp); ea_ = EA; ei_ = EI; rho = RHOA; w = 0.0_wp
    DO k = 1, NN
      x = lt*REAL(k - 1, wp)/REAL(NE, wp)
      q(6*k - 5:6*k) = [x, 0.0_wp, a0*SIN(PI*x/lt), (1.0_wp + eps), 0.0_wp, (1.0_wp + eps)*a0*PI/lt*COS(PI*x/lt)]
    END DO
    CALL CD_HermiteCable_Dyn_Init(model, l0, ea_, ei_, rho, w, q, [1, 2, 3, 6*NE + 1, 6*NE + 2, 6*NE + 3], &
                                  -1.0e4_wp, 0.0_wp, 1.0_wp, es, em)
    CALL require(es == CD_HCDYN_OK, 'init: '//TRIM(em))
  END SUBROUTINE setup

  SUBROUTINE decay_case(label, tens, zeta_ref)
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: tens, zeta_ref
    TYPE(CD_HermiteCableDynType) :: model
    REAL(wp) :: kk, om1, share, bb, zeta_expect, dt, z, zprev, zprev2, peaks(200), tpk(200), delta, zeta_meas, res
    INTEGER :: n, k, np, es, it, nstep
    CHARACTER(512) :: em
    kk = PI/LL
    om1 = SQRT((EI*kk**4 + tens*kk**2)/RHOA)
    share = EI*kk**2/(EI*kk**2 + tens)
    bb = 2.0_wp*zeta_ref*EI/om1
    zeta_expect = om1*bb/(2.0_wp*EI)*share
    CALL setup(model, tens, 1.0e-3_wp)
    CALL CD_HermiteCable_Dyn_Set_Bending_Damping(model, [(bb, k=1, NE)], es, em)
    CALL require(es == CD_HCDYN_OK, label//': set damping: '//TRIM(em))
    dt = (2.0_wp*PI/om1)/200.0_wp
    nstep = 200*12
    np = 0; zprev = model%q(6*(MID - 1) + 3); zprev2 = zprev
    DO n = 1, nstep
      CALL CD_HermiteCable_Dyn_Step(model, dt, 50, 1.0e-10_wp, es, em, it, res)
      IF (es /= CD_HCDYN_OK) THEN
        CALL require(.FALSE., label//': step failed: '//TRIM(em)); RETURN
      END IF
      z = model%q(6*(MID - 1) + 3)
      ! a positive peak of the midspan deflection (parabolic refinement of the sampled peak)
      IF (n >= 2 .AND. zprev > zprev2 .AND. zprev >= z .AND. zprev > 0.0_wp .AND. np < SIZE(peaks)) THEN
        np = np + 1
        peaks(np) = zprev - 0.125_wp*(z - zprev2)**2/(z - 2.0_wp*zprev + zprev2)
        tpk(np) = dt*REAL(n - 1, wp)
      END IF
      zprev2 = zprev; zprev = z
    END DO
    CALL require(np >= 8, label//': enough peaks')
    IF (np < 8) RETURN
    delta = LOG(peaks(2)/peaks(np))/REAL(np - 2, wp)
    zeta_meas = delta/SQRT(4.0_wp*PI*PI + delta*delta)
    WRITE (*, '(A,A,F7.4,A,ES11.4,A,ES11.4,A,F8.3)') label, ': bending share ', share, ', zeta expected ', &
      zeta_expect, ', measured ', zeta_meas, ', period ratio ', (tpk(np) - tpk(2))/REAL(np - 2, wp)*om1/(2.0_wp*PI)
    IF (zeta_ref > 0.0_wp) THEN
      CALL require(ABS(zeta_meas - zeta_expect) <= 0.01_wp*zeta_expect, label//': decay matches zeta to 1 %')
    ELSE
      CALL require(ABS(zeta_meas) <= 1.0e-5_wp, label//': no damping without Bb')
    END IF
  END SUBROUTINE decay_case

  SUBROUTINE check_refusals()
    TYPE(CD_HermiteCableDynType) :: model
    INTEGER :: es, k
    CHARACTER(512) :: em
    CALL setup(model, 1.0e4_wp, 0.02_wp)
    CALL CD_HermiteCable_Dyn_Set_Bending_Damping(model, [(-1.0_wp, k=1, NE)], es, em)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'finite and non-negative') > 0, 'B4: negative Bb refused')
    CALL CD_HermiteCable_Dyn_Set_Bending_Damping(model, [(1.0_wp, k=1, NE - 1)], es, em)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'length ne') > 0, 'B4: wrong length refused')
    CALL CD_HermiteCable_Dyn_Set_Bending_Damping(model, [(0.0_wp, k=1, NE)], es, em)
    CALL require(es == CD_HCDYN_OK .AND. .NOT. model%has_bend_damping, 'B3: all-zero Bb keeps it off')
  END SUBROUTINE check_refusals

  SUBROUTINE check_recovery()
    TYPE(CD_HermiteCableDynType) :: model
    REAL(wp) :: curv(NN), bb, v0, lt, x, kdot, f_damped(3), f_plain(3)
    INTEGER :: es, k
    CHARACTER(512) :: em
    v0 = 0.3_wp
    bb = 2.0e3_wp
    CALL setup(model, 1.0e4_wp, 0.0_wp)
    lt = model%q(6*(NN - 1) + 1)
    CALL CD_HermiteCable_Dyn_Set_Bending_Damping(model, [(bb, k=1, NE)], es, em)
    CALL require(es == CD_HCDYN_OK, 'B5: set damping: '//TRIM(em))
    ! straight line, first-mode velocity z' = v0 sin(pi x/lt) and its slope rate, no acceleration
    model%v = 0.0_wp
    DO k = 1, NN
      x = model%q(6*(k - 1) + 1)
      model%v(6*(k - 1) + 3) = v0*SIN(PI*x/lt)
      model%v(6*(k - 1) + 6) = v0*PI/lt*COS(PI*x/lt)
    END DO
    model%v([1, 2, 3, 6*NE + 1, 6*NE + 2, 6*NE + 3]) = 0.0_wp
    model%a = 0.0_wp
    CALL CD_HermiteCable_Dyn_Node_Recovery(model, es, em, curvature=curv)
    CALL require(es == CD_HCDYN_OK, 'B5: node recovery: '//TRIM(em))
    kdot = v0*(PI/lt)**2
    WRITE (*, '(A,ES11.3,A,ES11.3)') 'B5: max recovered curvature ', MAXVAL(ABS(curv)), &
      ' 1/m; the viscous moment over EI would read ', bb*kdot/EI
    CALL require(MAXVAL(ABS(curv)) <= 1.0e-6_wp*bb*kdot/EI, 'B5: no viscous curvature in the node recovery')
    CALL CD_HermiteCable_Dyn_End_Force(model, 1, f_damped, es, em)
    model%has_bend_damping = .FALSE.
    CALL CD_HermiteCable_Dyn_End_Force(model, 1, f_plain, es, em)
    WRITE (*, '(A,ES11.3,A)') 'B5: the bending damping changes the end force by ', NORM2(f_damped - f_plain), ' N'
    CALL require(NORM2(f_damped - f_plain) > 0.0_wp, 'B5: the end force carries the damping')
  END SUBROUTINE check_recovery

END PROGRAM test_hermite_bending_damping_dynamic
