! File: tests/test_wave_math.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_wave_math
  !! Gate for the block kernels of CableDyn_WaveMath against an independent oracle, the
  !! quadruple-precision intrinsics (a double argument is exact in quadruple precision):
  !!   1. sine and cosine of 400,000 phases over |theta| <= 4e4 (the phases of a 3-hour sea)
  !!      within 3 ulp of the correctly rounded values, and the identities sin^2 + cos^2 = 1
  !!      to 4 ulp and sin(-theta) = -sin(theta), cos(-theta) = cos(theta) bit for bit;
  !!   2. the quadrant boundaries (theta = n pi/2 rounded, n up to 2^19) and tiny phases
  !!      (sin theta = theta, cos theta = 1 below 2^-27, signed zero kept);
  !!   3. exp over -745 <= x <= 709 within 2 ulp (in the normal range), the subnormal tail and
  !!      underflow to zero as the intrinsic gives them, exp(0) = 1;
  !!   4. phases beyond the reduction's range (|theta| > 2^20 pi/2) take the intrinsic (to 1e-6:
  !!      a clamped reduction would be off by O(1));
  !!   5. a NaN argument gives NaN (sine, cosine and exponential) and an infinite phase NaN, as
  !!      the intrinsics: the clamp to the reduction's range must not turn a NaN into a finite
  !!      value (gfortran's MIN and MAX return the non-NaN argument).
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: INT64
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN, IEEE_POSITIVE_INF, IEEE_IS_NAN
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_WaveMath, ONLY: CD_Wave_SinCos, CD_Wave_Exp
  IMPLICIT NONE
  INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(30)
  INTEGER, PARAMETER :: N = 400000
  INTEGER :: nfail, i
  REAL(wp), ALLOCATABLE :: th(:), s(:), c(:), x(:), e(:), s2(:), c2(:)
  REAL(wp) :: ulp_s, ulp_c, ulp_e, ulp_id, pi_wp
  REAL(qp) :: pi_q

  nfail = 0
  ALLOCATE (th(N), s(N), c(N), x(N), e(N), s2(N), c2(N))
  ! 1. phases over |theta| <= 4e4, with irrational spacing
  DO i = 1, N
    th(i) = -4.0e4_wp + 8.0e4_wp*REAL(i, wp)/REAL(N, wp) + 1.0e-7_wp*REAL(i, wp)*SQRT(2.0_wp)
  END DO
  CALL CD_Wave_SinCos(N, th, s, c)
  ulp_s = 0.0_wp
  ulp_c = 0.0_wp
  ulp_id = 0.0_wp
  DO i = 1, N
    ulp_s = MAX(ulp_s, ulp_error(s(i), SIN(REAL(th(i), qp))))
    ulp_c = MAX(ulp_c, ulp_error(c(i), COS(REAL(th(i), qp))))
    ulp_id = MAX(ulp_id, REAL(ABS(REAL(s(i), qp)**2 + REAL(c(i), qp)**2 - 1.0_qp), wp)/EPSILON(1.0_wp))
  END DO
  WRITE (*, '(A,F6.3,A,F6.3,A,F6.3,A)') '    sin, cos over |theta| <= 4e4: max error ', ulp_s, ' and ', ulp_c, &
    ' ulp; sin^2 + cos^2 - 1 within ', ulp_id, ' ulp'
  CALL require(ulp_s <= 3.0_wp .AND. ulp_c <= 3.0_wp, 'sine and cosine within 3 ulp')
  CALL require(ulp_id <= 4.0_wp, 'sin^2 + cos^2 = 1 within 4 ulp')
  CALL CD_Wave_SinCos(N, -th, s2, c2)
  CALL require(ALL(s2 == -s) .AND. ALL(c2 == c), 'sine odd and cosine even bit for bit')

  ! 2. quadrant boundaries and tiny phases
  pi_wp = 4.0_wp*ATAN(1.0_wp)
  pi_q = 4.0_qp*ATAN(1.0_qp)
  DO i = 1, 2000
    th(i) = REAL(REAL(MOD(i*262, 524288) - 262144, qp)*pi_q/2.0_qp, wp)
  END DO
  CALL CD_Wave_SinCos(2000, th, s, c)
  ulp_s = 0.0_wp
  DO i = 1, 2000
    ulp_s = MAX(ulp_s, abs_error(s(i), SIN(REAL(th(i), qp)), REAL(th(i), qp)), &
                abs_error(c(i), COS(REAL(th(i), qp)), REAL(th(i), qp)))
  END DO
  WRITE (*, '(A,ES9.2)') '    at the quadrant boundaries (|n| <= 2^18): max |error| / ulp(theta) ', ulp_s
  CALL require(ulp_s <= 1.0_wp, 'quadrant boundaries: error below ulp(theta)')
  th(1:4) = [0.0_wp, -0.0_wp, 1.0e-10_wp, -2.0_wp**(-30)]
  CALL CD_Wave_SinCos(4, th, s, c)
  ! the sign bits themselves: ifx's SIGN treats -0 as +0 by default (/assume:nominus0)
  CALL require(s(1) == 0.0_wp .AND. s(2) == 0.0_wp .AND. .NOT. BTEST(TRANSFER(s(1), 0_INT64), 63) .AND. &
               BTEST(TRANSFER(s(2), 0_INT64), 63), 'sin(+-0) = +-0')
  CALL require(ALL(c(1:4) == 1.0_wp) .AND. s(3) == th(3) .AND. s(4) == th(4), 'tiny phases: sin = theta, cos = 1')

  ! 3. exponential
  DO i = 1, N
    x(i) = -745.0_wp + 1454.0_wp*REAL(i, wp)/REAL(N, wp) + 1.0e-9_wp*REAL(i, wp)
  END DO
  x(N) = 0.0_wp
  CALL CD_Wave_Exp(N, x, e)
  ulp_e = 0.0_wp
  DO i = 1, N
    IF (x(i) >= -708.0_wp .AND. x(i) <= 709.0_wp) THEN
      ulp_e = MAX(ulp_e, ulp_error(e(i), EXP(REAL(x(i), qp))))
    ELSE
      CALL require(e(i) == EXP(x(i)), 'exp outside the normal range is the intrinsic')
    END IF
  END DO
  WRITE (*, '(A,F6.3,A)') '    exp over -708 <= x <= 709: max error ', ulp_e, ' ulp'
  CALL require(ulp_e <= 2.0_wp, 'exp within 2 ulp')
  CALL require(e(N) == 1.0_wp, 'exp(0) = 1')
  x(1:3) = [-745.2_wp, -800.0_wp, -720.0_wp]
  CALL CD_Wave_Exp(3, x, e)
  CALL require(e(1) == EXP(x(1)) .AND. e(2) == 0.0_wp .AND. e(3) == EXP(x(3)) .AND. e(3) > 0.0_wp, &
               'subnormal tail and underflow as the intrinsic')

  ! 4. beyond the reduction's range
  th(1:3) = [2.0e6_wp, -3.5e7_wp, 1.0e9_wp]
  CALL CD_Wave_SinCos(3, th, s, c)
  ulp_s = 0.0_wp
  DO i = 1, 3
    ulp_s = MAX(ulp_s, REAL(ABS(REAL(s(i), qp) - SIN(REAL(th(i), qp))), wp), &
                REAL(ABS(REAL(c(i), qp) - COS(REAL(th(i), qp))), wp))
  END DO
  WRITE (*, '(A,ES9.2)') '    beyond 2^20 pi/2 (the intrinsic): max |error| ', ulp_s
  CALL require(ulp_s <= 1.0e-6_wp, 'phases beyond 2^20 pi/2 are not reduced as near ones')

  ! 5. NaN and infinite arguments, alone and within a block of finite ones
  th(1:5) = [0.5_wp, IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN), 1.0e3_wp, -IEEE_VALUE(1.0_wp, IEEE_POSITIVE_INF), 2.0_wp]
  CALL CD_Wave_SinCos(5, th, s, c)
  WRITE (*, '(A,4ES11.3)') '    sin, cos of NaN and of -Inf: ', s(2), c(2), s(4), c(4)
  CALL require(IEEE_IS_NAN(s(2)) .AND. IEEE_IS_NAN(c(2)) .AND. IEEE_IS_NAN(s(4)) .AND. IEEE_IS_NAN(c(4)), &
               'sine and cosine of NaN and of an infinite phase are NaN')
  CALL require(ABS(REAL(s(3), qp) - SIN(REAL(th(3), qp))) <= REAL(4.0_wp*SPACING(1.0_wp), qp), &
               'finite phases of the same block unaffected')
  x(1:3) = [-1.0_wp, IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN), 2.0_wp]
  CALL CD_Wave_Exp(3, x, e)
  WRITE (*, '(A,ES11.3)') '    exp of NaN: ', e(2)
  CALL require(IEEE_IS_NAN(e(2)) .AND. .NOT. IEEE_IS_NAN(e(1)) .AND. .NOT. IEEE_IS_NAN(e(3)), 'exp of NaN is NaN')

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'test_wave_math: ', nfail, ' check(s) FAILED'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'test_wave_math: all checks passed'

CONTAINS

  REAL(wp) FUNCTION ulp_error(value, exact) RESULT(ulps)
    !! |value - exact| in units of the last place of the correctly rounded exact value.
    REAL(wp), INTENT(IN) :: value
    REAL(qp), INTENT(IN) :: exact
    REAL(wp) :: rounded
    rounded = REAL(exact, wp)
    IF (rounded == 0.0_wp) THEN
      ulps = REAL(ABS(REAL(value, qp) - exact), wp)/TINY(1.0_wp)
    ELSE
      ulps = REAL(ABS(REAL(value, qp) - exact)/REAL(SPACING(rounded), qp), wp)
    END IF
  END FUNCTION ulp_error

  REAL(wp) FUNCTION abs_error(value, exact, arg) RESULT(rel)
    !! |value - exact| over ulp(arg): near a zero of sine or cosine the phase's own rounding
    !! (half an ulp of theta) bounds what any evaluation can resolve.
    REAL(wp), INTENT(IN) :: value
    REAL(qp), INTENT(IN) :: exact, arg
    rel = REAL(ABS(REAL(value, qp) - exact)/REAL(SPACING(MAX(1.0_wp, REAL(ABS(arg), wp))), qp), wp)
  END FUNCTION abs_error

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. cond) THEN
      nfail = nfail + 1
      WRITE (*, '(A)') 'FAIL: '//label
    END IF
  END SUBROUTINE require

END PROGRAM test_wave_math
