! File: tests/test_hermite_mass_damping_dynamic.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_hermite_mass_damping_dynamic
  !! Free decays of modes 1, 2 and 3 of a weightless, pinned-pinned finite-EI cable with Rayleigh
  !! damping: the mass-proportional transverse damping alpha (CD_HermiteCable_Dyn_Set_Mass_Damping)
  !! and the bending damping Bb (CD_HermiteCable_Dyn_Set_Bending_Damping), generalised-alpha with
  !! rho_inf = 1 (no numerical damping), 200 steps a period over 12 periods. For the mode
  !! sin(k x), k = n pi/L, of frequency omega^2 = (EI k^4 + T k^2)/m the damping ratio is
  !!   zeta = alpha/(2 omega) + omega Bb/(2 EI) * EI k^2/(EI k^2 + T),
  !! alpha acting on the whole modal mass (the cable's own) and Bb on the bending share.
  !!   R1 alpha alone on a tensioned cable: modes 1 and 3 decay at alpha/(2 omega) to 0.25 %;
  !!   R2 Bb alone on the same cable: modes 1 and 3 at the bending-share ratio to 0.25 %;
  !!   R3 both: the sum, to 0.25 %;
  !!   R4 the two-frequency pair (CD_HermiteCable_Rayleigh_Pair, the EI|-zeta|f1|f2 deck form) on
  !!      a cable without tension, f1 and f2 the frequencies of modes 1 and 3: both decay at the
  !!      target zeta to 0.25 %, and mode 2, between them, at the lower Rayleigh ratio to 0.25 %;
  !!   R5 a negative coefficient or a wrong length is refused, and an all-zero array keeps the
  !!      configuration off.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_HermiteCable, ONLY: CD_HermiteCable_Rayleigh_Pair
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_Set_Bending_Damping, &
                                          CD_HermiteCable_Dyn_Set_Mass_Damping, CD_HCDYN_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
  REAL(wp), PARAMETER :: LL = 20.0_wp, EA = 1.0e8_wp, EI = 1.0e5_wp, RHOA = 10.0_wp, TENS = 1.0e4_wp
  REAL(wp), PARAMETER :: ZT = 0.02_wp
  INTEGER, PARAMETER :: NE = 40, NN = NE + 1, NDOF = 6*NN, MID = NE/2 + 1, QUARTER = NE/4 + 1
  INTEGER :: nfail = 0
  REAL(wp) :: alpha, bb, f1, f3

  ! R1-R3: the tensioned cable; alpha and Bb each give zeta ZT to mode 1 in their own share
  alpha = 2.0_wp*ZT*omega(1, TENS)
  bb = 2.0_wp*ZT*EI/omega(1, TENS)
  CALL decay_case('R1 alpha alone', TENS, 1, alpha, 0.0_wp)
  CALL decay_case('R1 alpha alone', TENS, 3, alpha, 0.0_wp)
  CALL decay_case('R2 Bb alone', TENS, 1, 0.0_wp, bb)
  CALL decay_case('R2 Bb alone', TENS, 3, 0.0_wp, bb)
  CALL decay_case('R3 alpha and Bb', TENS, 1, alpha, bb)
  CALL decay_case('R3 alpha and Bb', TENS, 3, alpha, bb)
  ! R4: the pair holding ZT at the frequencies of modes 1 and 3 of the cable without tension
  f1 = omega(1, 0.0_wp)/(2.0_wp*PI)
  f3 = omega(3, 0.0_wp)/(2.0_wp*PI)
  CALL CD_HermiteCable_Rayleigh_Pair(ZT, f1, f3, EI, bb, alpha)
  CALL decay_case('R4 zeta at f1 and f3', 0.0_wp, 1, alpha, bb, ZT)
  CALL decay_case('R4 zeta at f1 and f3', 0.0_wp, 3, alpha, bb, ZT)
  CALL decay_case('R4 zeta at f1 and f3', 0.0_wp, 2, alpha, bb)
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: Hermite Rayleigh damping dynamic decay'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  PURE REAL(wp) FUNCTION omega(n, tens)
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(IN) :: tens
    REAL(wp) :: k
    k = REAL(n, wp)*PI/LL
    omega = SQRT((EI*k**4 + tens*k**2)/RHOA)
  END FUNCTION omega

  SUBROUTINE setup(model, tens, n, a0)
    !! The cable along x, pinned at both ends, at rest in its mode n of amplitude a0 (z).
    TYPE(CD_HermiteCableDynType), INTENT(OUT) :: model
    REAL(wp), INTENT(IN) :: tens, a0
    INTEGER, INTENT(IN) :: n
    REAL(wp) :: l0(NE), ea_(NE), ei_(NE), rho(NE), w(NE), q(NDOF), eps, lt, x, kn
    INTEGER :: k, es
    CHARACTER(512) :: em
    eps = tens/EA; lt = (1.0_wp + eps)*LL
    kn = REAL(n, wp)*PI/lt
    l0 = LL/REAL(NE, wp); ea_ = EA; ei_ = EI; rho = RHOA; w = 0.0_wp
    DO k = 1, NN
      x = lt*REAL(k - 1, wp)/REAL(NE, wp)
      q(6*k - 5:6*k) = [x, 0.0_wp, a0*SIN(kn*x), (1.0_wp + eps), 0.0_wp, (1.0_wp + eps)*a0*kn*COS(kn*x)]
    END DO
    CALL CD_HermiteCable_Dyn_Init(model, l0, ea_, ei_, rho, w, q, [1, 2, 3, 6*NE + 1, 6*NE + 2, 6*NE + 3], &
                                  -1.0e4_wp, 0.0_wp, 1.0_wp, es, em)
    CALL require(es == CD_HCDYN_OK, 'init: '//TRIM(em))
  END SUBROUTINE setup

  SUBROUTINE decay_case(label, tens, n, alpha_c, bb_c, zeta_target)
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: tens, alpha_c, bb_c
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(IN), OPTIONAL :: zeta_target
    TYPE(CD_HermiteCableDynType) :: model
    REAL(wp) :: kk, om, share, zeta_expect, dt, z, zprev, zprev2, peaks(200), tpk(200), delta, zeta_meas, res, sgn
    INTEGER :: m, k, np, es, it, node
    CHARACTER(512) :: em
    CHARACTER(16) :: tag
    kk = REAL(n, wp)*PI/LL
    om = omega(n, tens)
    share = EI*kk**2/(EI*kk**2 + tens)
    zeta_expect = alpha_c/(2.0_wp*om) + om*bb_c/(2.0_wp*EI)*share
    WRITE (tag, '(A,I0)') ', mode ', n
    CALL setup(model, tens, n, 1.0e-3_wp)
    IF (alpha_c > 0.0_wp) THEN
      CALL CD_HermiteCable_Dyn_Set_Mass_Damping(model, [(alpha_c, k=1, NE)], es, em)
      CALL require(es == CD_HCDYN_OK, label//TRIM(tag)//': set alpha: '//TRIM(em))
    END IF
    IF (bb_c > 0.0_wp) THEN
      CALL CD_HermiteCable_Dyn_Set_Bending_Damping(model, [(bb_c, k=1, NE)], es, em)
      CALL require(es == CD_HCDYN_OK, label//TRIM(tag)//': set Bb: '//TRIM(em))
    END IF
    ! an antinode of the mode: midspan for modes 1 and 3 (sign -1 for mode 3), the quarter for mode 2
    node = MID
    IF (n == 2) node = QUARTER
    sgn = SIGN(1.0_wp, SIN(kk*model%q(6*(node - 1) + 1)))
    dt = (2.0_wp*PI/om)/200.0_wp
    np = 0; zprev = sgn*model%q(6*(node - 1) + 3); zprev2 = zprev
    DO m = 1, 200*12
      CALL CD_HermiteCable_Dyn_Step(model, dt, 50, 1.0e-10_wp, es, em, it, res)
      IF (es /= CD_HCDYN_OK) THEN
        CALL require(.FALSE., label//TRIM(tag)//': step failed: '//TRIM(em)); RETURN
      END IF
      z = sgn*model%q(6*(node - 1) + 3)
      ! a positive peak of the antinode deflection (parabolic refinement of the sampled peak)
      IF (m >= 2 .AND. zprev > zprev2 .AND. zprev >= z .AND. zprev > 0.0_wp .AND. np < SIZE(peaks)) THEN
        np = np + 1
        peaks(np) = zprev - 0.125_wp*(z - zprev2)**2/(z - 2.0_wp*zprev + zprev2)
        tpk(np) = dt*REAL(m - 1, wp)
      END IF
      zprev2 = zprev; zprev = z
    END DO
    CALL require(np >= 8, label//TRIM(tag)//': enough peaks')
    IF (np < 8) RETURN
    delta = LOG(peaks(2)/peaks(np))/REAL(np - 2, wp)
    zeta_meas = delta/SQRT(4.0_wp*PI*PI + delta*delta)
    WRITE (*, '(A,A,F7.4,A,ES11.4,A,ES11.4,A,ES10.2,A,F8.5)') label//TRIM(tag), ': bending share ', share, &
      ', zeta expected ', zeta_expect, ', measured ', zeta_meas, ', rel. error ', &
      (zeta_meas - zeta_expect)/zeta_expect, ', period ratio ', (tpk(np) - tpk(2))/REAL(np - 2, wp)*om/(2.0_wp*PI)
    CALL require(ABS(zeta_meas - zeta_expect) <= 2.5e-3_wp*zeta_expect, &
                 label//TRIM(tag)//': decay matches zeta to 0.25 %')
    IF (PRESENT(zeta_target)) THEN
      WRITE (*, '(A,ES11.4,A,ES10.2)') label//TRIM(tag)//': target zeta ', zeta_target, ', measured rel. error ', &
        (zeta_meas - zeta_target)/zeta_target
      CALL require(ABS(zeta_meas - zeta_target) <= 2.5e-3_wp*zeta_target, &
                   label//TRIM(tag)//': the target zeta to 0.25 %')
    END IF
  END SUBROUTINE decay_case

  SUBROUTINE check_refusals()
    TYPE(CD_HermiteCableDynType) :: model
    INTEGER :: es, k
    CHARACTER(512) :: em
    CALL setup(model, TENS, 1, 0.02_wp)
    CALL CD_HermiteCable_Dyn_Set_Mass_Damping(model, [(-1.0_wp, k=1, NE)], es, em)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'finite and non-negative') > 0, 'R5: negative alpha refused')
    CALL CD_HermiteCable_Dyn_Set_Mass_Damping(model, [(1.0_wp, k=1, NE - 1)], es, em)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'length ne') > 0, 'R5: wrong length refused')
    CALL CD_HermiteCable_Dyn_Set_Mass_Damping(model, [(0.5_wp, k=1, NE)], es, em)
    CALL require(es == CD_HCDYN_OK .AND. model%has_mass_damping, 'R5: a positive alpha switches it on')
    CALL CD_HermiteCable_Dyn_Set_Mass_Damping(model, [(0.0_wp, k=1, NE)], es, em)
    CALL require(es == CD_HCDYN_OK .AND. .NOT. model%has_mass_damping, 'R5: all-zero alpha keeps it off')
  END SUBROUTINE check_refusals

END PROGRAM test_hermite_mass_damping_dynamic
