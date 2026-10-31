! File: src/CableDyn_WaveMath.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_WaveMath
  !! The sines, cosines and exponentials of the wave kinematics, evaluated in blocks.
  !!
  !! An irregular sea sums, at every point and time, one phase cosine and sine and one depth
  !! exponential per component: thousands of transcendental evaluations per node and step. The
  !! intrinsic functions become one scalar library call per element where the compiler may not
  !! use its vector library (Intel ifx under /fp:precise, as every ifx build compiles), and the
  !! MinGW-w64 library of gfortran reduces large arguments in x87 precision (errors of 10^4 to
  !! 10^5 ulp at |x| ~ 10^4). Here each block is plain IEEE arithmetic that compilers vectorise
  !! under any value-safe floating-point model (a value-unsafe one may fold the rounding to an
  !! integer below; the build compiles this file value-safe with every compiler):
  !!
  !! * sine and cosine together, from ONE argument reduction (they are always needed in pairs):
  !!   theta = n pi/2 + r with n the integer nearest 2 theta/pi and
  !!   r = ((theta - n P1) - n P2) - n P3, the Cody-Waite reduction with pi/2 split in three
  !!   parts whose products with |n| < 2^20 are exact (P1 and P2 carry 33 bits); sin r and cos r
  !!   by their Taylor series through r^17 and r^18 (|r| <= pi/4, truncation below 2e-19); the
  !!   quadrant n mod 4 selects and signs them;
  !! * exp(x) = 2^n exp(r) with n the integer nearest x/ln 2, r = (x - n L1) - n L2 (ln 2 in two
  !!   parts, |r| <= ln(2)/2), exp r by its Taylor series through r^13 (truncation below 1e-17)
  !!   and 2^n from its IEEE bits.
  !!
  !! Against quadruple precision (test_wave_math, and sweeps of 10^6 phases per decade) the
  !! results are within 2.4 ulp over |theta| <= 2^20 pi/2 (sine, cosine) and 1.1 ulp over
  !! -708 <= x <= 709 (exponential). Next to a zero of the sine or cosine, the double nearest
  !! n pi/2, the value is as small as theta's distance from the multiple (down to 6.2e-19) and the
  !! part of pi/2 below P1 + P2 + P3, 8.5e-32, adds |n| 8.5e-32 (at most 9e-26) to the error:
  !! up to 2.8e6 ulp of that tiny value, 4e-10 ulp of 1. A phase formed as k x - omega t + phi
  !! already carries a rounding error of |theta| ulp(1), far above both.
  !! Elements beyond the reduction (|theta| > 2^20 pi/2) or the normal range (x < -708,
  !! x > 709), and NaN, take the intrinsic, so a NaN argument gives NaN as the intrinsic does (the
  !! clamp to the range may turn it into a finite bound: Fortran leaves MIN and MAX of a NaN
  !! unspecified, and gfortran returns the other argument). The arrays are explicit-shape:
  !! inlined into a caller's loops, the blocks then still vectorise.
  USE CableDyn_Precision, ONLY: wp
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CD_Wave_SinCos, CD_Wave_Exp

  REAL(wp), PARAMETER :: TWO_OVER_PI = 0.63661977236758134308_wp
  ! pi/2 = P1 + P2 + P3 (+ 8.5e-32): P1 and P2 carry 33 significant bits (fdlibm's pio2_1, pio2_2)
  REAL(wp), PARAMETER :: P1 = 1.57079632673412561417e+00_wp
  REAL(wp), PARAMETER :: P2 = 6.07710050630396597660e-11_wp
  REAL(wp), PARAMETER :: P3 = 2.02226624871116645580e-21_wp
  REAL(wp), PARAMETER :: THETA_MAX = 1.6470993291652855e+06_wp   ! 2^20 pi/2
  REAL(wp), PARAMETER :: INV_LN2 = 1.44269504088896338700e+00_wp
  ! ln 2 = L1 + L2: L1 carries 32 significant bits (fdlibm's ln2_hi), so n L1 is exact for |n| < 2^21
  REAL(wp), PARAMETER :: L1 = 6.93147180369123816490e-01_wp
  REAL(wp), PARAMETER :: L2 = 1.90821492927058770002e-10_wp
  REAL(wp), PARAMETER :: X_LOW = -708.0_wp, X_HIGH = 709.0_wp
  ! (y + ROUNDER) - ROUNDER is y rounded to the nearest integer for |y| < 2^51, in plain IEEE
  ! arithmetic (ANINT can become a library call that blocks the vector loop)
  REAL(wp), PARAMETER :: ROUNDER = 6755399441055744.0_wp   ! 1.5 * 2^52
  ! Taylor coefficients 1/m!
  REAL(wp), PARAMETER :: F2 = 0.5_wp, F3 = 1.0_wp/6.0_wp, F4 = 1.0_wp/24.0_wp, F5 = 1.0_wp/120.0_wp, &
                         F6 = 1.0_wp/720.0_wp, F7 = 1.0_wp/5040.0_wp, F8 = 1.0_wp/40320.0_wp, &
                         F9 = 1.0_wp/362880.0_wp, F10 = 1.0_wp/3628800.0_wp, F11 = 1.0_wp/39916800.0_wp, &
                         F12 = 1.0_wp/479001600.0_wp, F13 = 1.0_wp/6227020800.0_wp, &
                         F14 = 1.0_wp/87178291200.0_wp, F15 = 1.0_wp/1307674368000.0_wp, &
                         F16 = 1.0_wp/20922789888000.0_wp, F17 = 1.0_wp/355687428096000.0_wp, &
                         F18 = 1.0_wp/6402373705728000.0_wp

CONTAINS

  PURE SUBROUTINE CD_Wave_SinCos(n, theta, s, c)
    !! s = sin(theta), c = cos(theta) for the n elements (explicit-shape arrays, so that the
    !! loop vectorises wherever the call is inlined).
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(IN) :: theta(n)
    REAL(wp), INTENT(OUT) :: s(n), c(n)
    INTEGER :: i, q
    REAL(wp) :: th, fn, r, r2, ps, pc, ss, cc
    LOGICAL :: far

    DO i = 1, n
      th = MAX(-THETA_MAX, MIN(THETA_MAX, theta(i)))
      fn = (th*TWO_OVER_PI + ROUNDER) - ROUNDER
      r = ((th - fn*P1) - fn*P2) - fn*P3
      r2 = r*r
      ! Horner's scheme, its innermost terms first
      ps = F13 + r2*(-F15 + r2*F17)
      ps = r + (r*r2)*(-F3 + r2*(F5 + r2*(-F7 + r2*(F9 + r2*(-F11 + r2*ps)))))
      pc = F12 + r2*(-F14 + r2*(F16 - r2*F18))
      pc = (1.0_wp - F2*r2) + (r2*r2)*(F4 + r2*(-F6 + r2*(F8 + r2*(-F10 + r2*pc))))
      ps = MERGE(ps, r, ABS(r) > 0.0_wp)   ! sin(-0) = -0 (ifx's SIGN would drop the sign)
      ! quadrant q = n mod 4: (sin, cos) = (ps, pc), (pc, -ps), (-ps, -pc), (-pc, ps)
      q = IAND(INT(fn), 3)
      ss = MERGE(pc, ps, BTEST(q, 0))
      cc = MERGE(ps, pc, BTEST(q, 0))
      s(i) = MERGE(-ss, ss, BTEST(q, 1))
      c(i) = MERGE(-cc, cc, BTEST(q + 1, 1))
    END DO
    ! beyond the reduction, or NaN (every comparison with a NaN is false)
    far = .FALSE.
    DO i = 1, n
      far = far .OR. .NOT. (ABS(theta(i)) <= THETA_MAX)
    END DO
    IF (.NOT. far) RETURN
    DO i = 1, n
      IF (.NOT. (ABS(theta(i)) <= THETA_MAX)) THEN
        s(i) = SIN(theta(i))
        c(i) = COS(theta(i))
      END IF
    END DO
  END SUBROUTINE CD_Wave_SinCos

  PURE SUBROUTINE CD_Wave_Exp(n, x, e)
    !! e = exp(x) for the n elements (explicit-shape arrays, as CD_Wave_SinCos).
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(IN) :: x(n)
    REAL(wp), INTENT(OUT) :: e(n)
    INTEGER :: i
    INTEGER(INT64) :: m
    REAL(wp) :: xc, fn, r, p
    LOGICAL :: outside

    DO i = 1, n
      xc = MAX(X_LOW, MIN(X_HIGH, x(i)))
      fn = (xc*INV_LN2 + ROUNDER) - ROUNDER
      r = (xc - fn*L1) - fn*L2
      p = F9 + r*(F10 + r*(F11 + r*(F12 + r*F13)))
      p = 1.0_wp + r*(1.0_wp + r*(F2 + r*(F3 + r*(F4 + r*(F5 + r*(F6 + r*(F7 + r*(F8 + r*p))))))))
      ! times 2^m from its IEEE bits, m in [-1022, 1023] for the clamped argument (SCALE can
      ! become a library call per element)
      m = INT(fn, INT64) + 1023_INT64
      e(i) = p*TRANSFER(ISHFT(m, 52), p)
    END DO
    ! outside the normal range, or NaN
    outside = .FALSE.
    DO i = 1, n
      outside = outside .OR. .NOT. (x(i) >= X_LOW .AND. x(i) <= X_HIGH)
    END DO
    IF (.NOT. outside) RETURN
    DO i = 1, n
      IF (.NOT. (x(i) >= X_LOW .AND. x(i) <= X_HIGH)) e(i) = EXP(x(i))
    END DO
  END SUBROUTINE CD_Wave_Exp

END MODULE CableDyn_WaveMath
