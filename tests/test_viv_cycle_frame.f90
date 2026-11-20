! File: tests/test_viv_cycle_frame.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_cycle_frame
  !! The cycle-mean flow frame kernel (CD_VIV_Cycle_Init, CD_VIV_Cycle_Update, CD_VIV_Cycle_Valid):
  !! vbar(t) = (1/T) int_{t-T}^{t} v dt' + (v(t) - v(t-T))/2 over the window T, from a ring of
  !! CD_VIV_CYC_NB samples per window of the velocity's integral.
  !! * C1 a component of period T and its second harmonic (a cross-flow VIV at f = 1/T and an
  !!   in-line one at 2f) on a constant drop out: vbar is the constant within 1 % of the harmonics'
  !!   amplitude after the first window (measured 0.25 %; what is left is the sampled derivative of
  !!   v(t - T), the quadratic through three samples, of order (omega h)^2/6 with h = T/NB);
  !! * C2 a ramp v = a + b t passes without lag: vbar = v(t) to round-off (1e-12; the first-order
  !!   filter of time T lags it by b T); with the harmonics of C1 on it within 1 %;
  !! * C3 a quadratic v = c t^2 gives c (t^2 - T^2/6), the kernel's second-order error, within
  !!   c h^2 (the sampled integral of a cubic);
  !! * C4 C2 and C3 hold for steps longer than a sample interval (dt = 1.7 T/NB) and for a step that
  !!   does not divide T;
  !! * C5 a node at rest at v0 (Init) gives vbar = v0 at once; the state stays valid
  !!   (CD_VIV_Cycle_Valid) over a long run, and a state with a sample time outside [0, 1] or a
  !!   slot off the ring is invalid.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_VIV, ONLY: CD_VIV_Cycle_Init, CD_VIV_Cycle_Update, CD_VIV_Cycle_Valid, CD_VIV_NMS_CYC, CD_VIV_CYC_NB
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.141592653589793238_wp, TW = 0.4_wp
  INTEGER :: nfail
  REAL(wp) :: st(CD_VIV_NMS_CYC), vs(3), err
  nfail = 0

  ! C5 at rest
  CALL CD_VIV_Cycle_Init(st, [0.3_wp, -0.2_wp, 0.1_wp])
  CALL CD_VIV_Cycle_Update(st, [0.3_wp, -0.2_wp, 0.1_wp], [0.3_wp, -0.2_wp, 0.1_wp], 1.0e-3_wp, TW, vs)
  CALL require(MAXVAL(ABS(vs - [0.3_wp, -0.2_wp, 0.1_wp])) <= 1.0e-14_wp, 'C5: a node at rest keeps its velocity')

  ! C1 harmonics on a constant
  err = run_case(1, 1.0e-3_wp)
  WRITE (*, '(A,ES10.3)') '  C1 harmonics on a constant, max |vbar - v_slow|/amp: ', err
  CALL require(err <= 1.0e-2_wp, 'C1: a component of period T and its second harmonic drop out')
  ! C2 ramp
  err = run_case(2, 1.0e-3_wp)
  WRITE (*, '(A,ES10.3)') '  C2 ramp, max |vbar - v|: ', err
  CALL require(err <= 1.0e-12_wp, 'C2: a ramp passes without lag')
  err = run_case(4, 1.0e-3_wp)
  WRITE (*, '(A,ES10.3)') '  C2 ramp with harmonics, max |vbar - v_slow|/amp: ', err
  CALL require(err <= 1.0e-2_wp, 'C2: a ramp with harmonics on it passes without lag')
  ! C3 quadratic
  err = run_case(3, 1.0e-3_wp)
  WRITE (*, '(A,ES10.3)') '  C3 quadratic, max |vbar - c (t^2 - T^2/6)| / (c h^2): ', err
  CALL require(err <= 1.0_wp, 'C3: a quadratic gives the second-order error c T^2/6')
  ! C4 long steps that do not divide T
  err = run_case(2, 1.7_wp*TW/REAL(CD_VIV_CYC_NB, wp))
  CALL require(err <= 1.0e-12_wp, 'C4: a ramp passes without lag with steps longer than a sample interval')
  err = run_case(3, 1.7_wp*TW/REAL(CD_VIV_CYC_NB, wp))
  CALL require(err <= 1.0_wp, 'C4: the quadratic error holds with steps longer than a sample interval')

  ! C5 validity
  CALL CD_VIV_Cycle_Init(st, [0.0_wp, 0.0_wp, 0.0_wp])
  CALL require(CD_VIV_Cycle_Valid(st), 'C5: the initial state is valid')
  st(1) = 1.5_wp
  CALL require(.NOT. CD_VIV_Cycle_Valid(st), 'C5: a sample time beyond one interval is invalid')
  st(1) = 0.5_wp
  st(2) = REAL(CD_VIV_CYC_NB + 3, wp)
  CALL require(.NOT. CD_VIV_Cycle_Valid(st), 'C5: a slot off the ring is invalid')
  st(2) = 2.5_wp
  CALL require(.NOT. CD_VIV_Cycle_Valid(st), 'C5: a fractional slot is invalid')

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: the cycle-mean flow frame kernel'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. cond) THEN
      WRITE (*, '(A)') 'MISMATCH: '//label
      nfail = nfail + 1
    END IF
  END SUBROUTINE require

  PURE FUNCTION signal(kind, t) RESULT(v)
    !! The test velocity (three components) and, in v(4:6), its slow part (the expected vbar).
    INTEGER, INTENT(IN) :: kind
    REAL(wp), INTENT(IN) :: t
    REAL(wp) :: v(6), w, h(3)
    w = 2.0_wp*PI/TW
    h = [0.0_wp, 0.5_wp*SIN(w*t + 0.3_wp), 0.2_wp*COS(2.0_wp*w*t)]
    SELECT CASE (kind)
    CASE (1)
      v(4:6) = [0.25_wp, -0.1_wp, 0.05_wp]
      v(1:3) = v(4:6) + h
    CASE (2)
      v(4:6) = [0.25_wp + 0.03_wp*t, -0.1_wp - 0.02_wp*t, 0.05_wp]
      v(1:3) = v(4:6)
    CASE (3)
      v(1:3) = [0.02_wp*t*t, -0.01_wp*t*t, 0.0_wp]
      v(4:6) = v(1:3) - [0.02_wp, -0.01_wp, 0.0_wp]*TW*TW/6.0_wp
    CASE DEFAULT
      v(4:6) = [0.25_wp + 0.03_wp*t, -0.1_wp - 0.02_wp*t, 0.05_wp]
      v(1:3) = v(4:6) + h
    END SELECT
  END FUNCTION signal

  REAL(wp) FUNCTION run_case(kind, dt) RESULT(err)
    !! The largest error of vbar against the slow part over 20 windows, after the first two
    !! (a ramp or quadratic case starts from its value continued backwards in time exactly: the
    !! ring is primed with the integral of the signal over the window before t = 0).
    INTEGER, INTENT(IN) :: kind
    REAL(wp), INTENT(IN) :: dt
    REAL(wp) :: s(CD_VIV_NMS_CYC), t, v0(6), v1(6), out(3), scale, h
    INTEGER :: n, k
    CALL prime(kind, s)
    h = TW/REAL(CD_VIV_CYC_NB, wp)
    ! the harmonics' amplitude; the ramp absolute (round-off); the quadratic's c h^2
    SELECT CASE (kind)
    CASE (1, 4)
      scale = 0.5_wp
    CASE (3)
      scale = 0.02_wp*h*h
    CASE DEFAULT
      scale = 1.0_wp
    END SELECT
    err = 0.0_wp
    n = CEILING(20.0_wp*TW/dt)
    t = 0.0_wp
    DO k = 1, n
      v0 = signal(kind, t)
      v1 = signal(kind, t + dt)
      CALL CD_VIV_Cycle_Update(s, v0(1:3), v1(1:3), dt, TW, out)
      t = t + dt
      IF (t > 2.0_wp*TW) err = MAX(err, MAXVAL(ABS(out - v1(4:6)))/scale)
      CALL require(CD_VIV_Cycle_Valid(s), 'C5: the state stays valid')
    END DO
  END FUNCTION run_case

  SUBROUTINE prime(kind, s)
    !! The ring at t = 0 for the signal continued backwards: sample j (j - 1 intervals h back) holds
    !! (I(-(j - 1) h) - I(0))/h, the integral by Simpson's rule (exact for these polynomials and
    !! close for the harmonics, which the first windows flush anyway).
    INTEGER, INTENT(IN) :: kind
    REAL(wp), INTENT(OUT) :: s(CD_VIV_NMS_CYC)
    REAL(wp) :: h, a, b, va(6), vm(6), vb(6), acc(3)
    INTEGER :: j
    CALL CD_VIV_Cycle_Init(s, [0.0_wp, 0.0_wp, 0.0_wp])
    h = TW/REAL(CD_VIV_CYC_NB, wp)
    acc = 0.0_wp
    s(3:5) = 0.0_wp
    DO j = 2, CD_VIV_CYC_NB + 2
      b = -REAL(j - 2, wp)*h
      a = b - h
      va = signal(kind, a)
      vm = signal(kind, 0.5_wp*(a + b))
      vb = signal(kind, b)
      acc = acc - h*(va(1:3) + 4.0_wp*vm(1:3) + vb(1:3))/6.0_wp
      s(3 + 3*(j - 1):5 + 3*(j - 1)) = acc/h
    END DO
  END SUBROUTINE prime

END PROGRAM test_viv_cycle_frame
