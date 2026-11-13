! File: tests/test_viv_wake_kernel.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_wake_kernel
  !! The wake-oscillator kernel CD_VIV_Wake_Advance (Facchinetti, de Langre and Biolley, J. Fluids
  !! Struct. 19, 2004, 123-140) against its published properties, with eps 0.3 and A 12:
  !!   K1 the free oscillator (no forcing) from q = 0.1 settles on the limit cycle q = 2 cos(om t):
  !!      amplitude 2 and period 2 pi/om within 1 %;
  !!   K2 forced at its own frequency by a structure moving as y0 D cos(om t), the lift
  !!      magnification K, half the first-harmonic amplitude of q, follows the closed form of their
  !!      Eq. 17, a one-harmonic balance, for y0 = 0.05-0.5;
  !!   K3 one call over an interval equals the same interval in two calls with the forcing split
  !!      at its midpoint, and refining the sub-steps converges at fourth order;
  !!   K4 a large eps (20, om 1) from (2, 0), whose damping rate eps om (q^2 - 1) is sixty times the
  !!      oscillation, stays finite in 0.1 s calls and matches 0.001 s calls within 1e-6 at t = 1.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_VIV, ONLY: CD_VIV_Wake_Advance, CD_VIV_WAKE_DEFAULT_EPS, CD_VIV_WAKE_DEFAULT_A
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
  REAL(wp), PARAMETER :: EPS = CD_VIV_WAKE_DEFAULT_EPS, A = CD_VIV_WAKE_DEFAULT_A
  REAL(wp), PARAMETER :: D = 0.1_wp, OM = 2.0_wp*PI ! 1 Hz wake frequency
  INTEGER :: nfail = 0

  CALL check_free()
  CALL check_forced()
  CALL check_steps()
  CALL check_stiff()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: wake-oscillator kernel'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE check_free()
    REAL(wp) :: s(2), dt, qmax, t, t_up(2), q_old
    INTEGER :: i, n_up
    s = [0.1_wp, 0.0_wp]
    dt = 0.01_wp
    qmax = 0.0_wp
    n_up = 0
    t_up = 0.0_wp
    DO i = 1, 8000 ! 80 periods
      q_old = s(1)
      CALL CD_VIV_Wake_Advance(s, OM, EPS, A/D, 0.0_wp, 0.0_wp, dt)
      t = REAL(i, wp)*dt
      IF (t < 70.0_wp) CYCLE
      qmax = MAX(qmax, ABS(s(1)))
      ! upward zero crossings, interpolated
      IF (q_old < 0.0_wp .AND. s(1) >= 0.0_wp .AND. n_up < 2) THEN
        n_up = n_up + 1
        t_up(n_up) = t - dt*s(1)/(s(1) - q_old)
      END IF
    END DO
    WRITE (*, '(A,F9.5,A,F9.5,A)') 'K1: limit-cycle amplitude ', qmax, ', period ', t_up(2) - t_up(1), ' s'
    CALL require(ABS(qmax/2.0_wp - 1.0_wp) <= 0.01_wp, 'K1: the free limit cycle has the amplitude 2')
    CALL require(n_up == 2 .AND. ABS(t_up(2) - t_up(1) - 1.0_wp) <= 0.01_wp, 'K1: and the wake period')
  END SUBROUTINE check_free

  REAL(wp) FUNCTION eq17(y0) RESULT(k)
    !! Facchinetti et al. (2004), Eq. 17: K = (X/36)^(1/3) + (4/(3X))^(1/3),
    !! X = 9 (A/eps) y0 + sqrt((9 (A/eps) y0)^2 - 48).
    REAL(wp), INTENT(IN) :: y0
    REAL(wp) :: x, b
    b = 9.0_wp*(A/EPS)*y0
    x = b + SQRT(b*b - 48.0_wp)
    k = (x/36.0_wp)**(1.0_wp/3.0_wp) + (4.0_wp/(3.0_wp*x))**(1.0_wp/3.0_wp)
  END FUNCTION eq17

  SUBROUTINE check_forced()
    !! Eq. 17 keeps the first harmonic only; the cubic damping term also makes a third harmonic,
    !! so the computed K exceeds it by a share that grows with the forcing (about 2 % at
    !! y0 = 0.5). The step is converged: halving it changes K by under 1e-4 of itself.
    REAL(wp), PARAMETER :: Y0S(5) = [0.05_wp, 0.1_wp, 0.2_wp, 0.3_wp, 0.5_wp]
    REAL(wp) :: k_num, k_ref, worst, k_half
    INTEGER :: j
    worst = 0.0_wp
    DO j = 1, SIZE(Y0S)
      k_num = forced_k(Y0S(j), 0.005_wp)
      k_ref = eq17(Y0S(j))
      worst = MAX(worst, ABS(k_num/k_ref - 1.0_wp))
      WRITE (*, '(A,F5.2,A,F8.4,A,F8.4)') 'K2: y0 ', Y0S(j), ': K ', k_num, ', Eq. 17 ', k_ref
    END DO
    k_half = forced_k(0.5_wp, 0.0025_wp)
    WRITE (*, '(A,ES10.3)') 'K2: change of K at y0 0.5 with half the step: ', ABS(k_half - k_num)
    CALL require(worst <= 0.02_wp, 'K2: the forced lift magnification follows Eq. 17 within 2 %')
    CALL require(ABS(k_half - k_num) <= 1.0e-4_wp*k_num, 'K2: with a converged step')
  END SUBROUTINE check_forced

  REAL(wp) FUNCTION forced_k(y0, dt) RESULT(k)
    !! Half the first-harmonic amplitude of q over periods 91-100 of the forcing y0 D cos(om t).
    REAL(wp), INTENT(IN) :: y0, dt
    REAL(wp) :: s(2), t, c, sn
    INTEGER :: i, n, n_keep
    n = NINT(100.0_wp/dt)
    n_keep = NINT(10.0_wp/dt)
    s = [2.0_wp, 0.0_wp]
    c = 0.0_wp
    sn = 0.0_wp
    DO i = 1, n
      t = REAL(i - 1, wp)*dt
      CALL CD_VIV_Wake_Advance(s, OM, EPS, A/D, -OM*OM*y0*D*COS(OM*t), -OM*OM*y0*D*COS(OM*(t + dt)), dt)
      IF (i > n - n_keep) THEN
        c = c + s(1)*COS(OM*(t + dt))
        sn = sn + s(1)*SIN(OM*(t + dt))
      END IF
    END DO
    k = 0.5_wp*(2.0_wp/REAL(n_keep, wp))*SQRT(c*c + sn*sn)
  END FUNCTION forced_k

  SUBROUTINE check_steps()
    REAL(wp) :: s1(2), s2(2), sc(2), sf(2), h, e1, e2
    h = 0.4_wp
    s1 = [1.3_wp, -4.0_wp]
    s2 = s1
    CALL CD_VIV_Wake_Advance(s1, OM, EPS, A/D, 0.2_wp, -0.6_wp, h)
    CALL CD_VIV_Wake_Advance(s2, OM, EPS, A/D, 0.2_wp, -0.2_wp, 0.5_wp*h)
    CALL CD_VIV_Wake_Advance(s2, OM, EPS, A/D, -0.2_wp, -0.6_wp, 0.5_wp*h)
    WRITE (*, '(A,ES10.3)') 'K3: one call against two halves: ', MAXVAL(ABS(s1 - s2))
    CALL require(MAXVAL(ABS(s1 - s2)) <= 1.0e-12_wp*MAXVAL(ABS(s1)), 'K3: an interval in two calls is one call')
    ! fourth order: calls short enough for one sub-step each (om h/N <= 0.1) set the step; 32 and
    ! 64 of them against 1024
    sf = march(1024, h)
    e1 = MAXVAL(ABS(march(32, h) - sf))
    e2 = MAXVAL(ABS(march(64, h) - sf))
    WRITE (*, '(A,2ES11.3,A,F7.2)') 'K3: errors with 32 and 64 steps ', e1, e2, ', ratio ', e1/e2
    CALL require(e1 < 1.0e-5_wp*MAXVAL(ABS(sf)), 'K3: the advance is accurate at om dt = 0.08')
    CALL require(e1/e2 > 12.0_wp, 'K3: and converges at fourth order (ratio near 16)')
    ! the kernel's own sub-steps in one call match the reference as closely
    sc = [1.3_wp, -4.0_wp]
    CALL CD_VIV_Wake_Advance(sc, OM, EPS, A/D, 0.2_wp, -0.6_wp, h)
    CALL require(MAXVAL(ABS(sc - sf)) < 1.0e-5_wp*MAXVAL(ABS(sf)), 'K3: one call sub-steps as finely')
  END SUBROUTINE check_steps

  FUNCTION march(n, h) RESULT(x)
    !! The K3 interval in n equal calls, the forcing linear from 0.2 to -0.6 over it.
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(IN) :: h
    REAL(wp) :: x(2)
    INTEGER :: k
    x = [1.3_wp, -4.0_wp]
    DO k = 1, n
      CALL CD_VIV_Wake_Advance(x, OM, EPS, A/D, 0.2_wp - 0.8_wp*REAL(k - 1, wp)/REAL(n, wp), &
                               0.2_wp - 0.8_wp*REAL(k, wp)/REAL(n, wp), h/REAL(n, wp))
    END DO
  END FUNCTION march

  SUBROUTINE check_stiff()
    !! The sub-steps follow the damping rate as well as the oscillation (eps 20, om 1).
    REAL(wp) :: sc(2), sf(2)
    INTEGER :: k
    sc = [2.0_wp, 0.0_wp]
    sf = sc
    DO k = 1, 10
      CALL CD_VIV_Wake_Advance(sc, 1.0_wp, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.1_wp)
    END DO
    DO k = 1, 1000
      CALL CD_VIV_Wake_Advance(sf, 1.0_wp, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.001_wp)
    END DO
    WRITE (*, '(A,2ES12.4,A,ES10.3)') 'K4: eps 20 at t = 1 s: ', sc, ', against 1 ms calls ', MAXVAL(ABS(sc - sf))
    CALL require(ALL(ABS(sc) < 10.0_wp) .AND. MAXVAL(ABS(sc - sf)) <= 1.0e-6_wp, 'K4: a large eps stays '// &
                 'resolved in long calls')
  END SUBROUTINE check_stiff

END PROGRAM test_viv_wake_kernel
