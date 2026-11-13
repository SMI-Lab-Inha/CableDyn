! File: tests/test_viv_amplitude_kernel.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_amplitude_kernel
  !! The cross-flow amplitude estimate of Thorsen et al. (2016, Eq. 8 and 9; CD_VIV_Amplitude_Update)
  !! and the amplitude tables (CD_VIV_Table, CD_VIV_Check_Amp_Table):
  !!   K1 a sinusoidal cross-flow velocity y' = om y0 cos(om t + p) with a step that does not divide
  !!      the period: after the first half cycle the half-cycle amplitude is y0 within 3e-6 relative
  !!      (the trapezoidal (om dt)^2/12 is 1.9e-6),
  !!      halving the step cuts the error about four times (trapezoidal travel, interpolated
  !!      crossings), and the smoothed amplitude follows Eq. 9 (within 1e-9 of y0 after 400 steps);
  !!   K2 a velocity that touches zero without changing sign ends no half cycle, a velocity that
  !!      changes sign does, also through a sample at exactly zero (the travel is signed);
  !!   K3 a growing oscillation y0(t) = a t: each half-cycle value is the mean amplitude of its half
  !!      cycle (within 0.1 % from the sixth period, the extremes then within that of a t);
  !!   K4 the table interpolates linearly and holds its end values; a one-row table is constant;
  !!      a NaN argument stays inside the table (NaN back, or the one row);
  !!   K5 the table check refuses each bad table by name and accepts a good one.
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN, IEEE_IS_NAN
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_VIV, ONLY: CD_VIV_Amplitude_Update, CD_VIV_Table, CD_VIV_Check_Amp_Table, CD_VIV_AMP_MAXP, &
                          CD_VIV_AMP_KDT
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
  INTEGER :: nfail = 0

  CALL check_sinusoid()
  CALL check_touch()
  CALL check_growth()
  CALL check_table()
  CALL check_tables_validated()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: VIV amplitude estimate and tables'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE run_sinusoid(dt, nstep, y0, err, smooth_err)
    !! The largest relative error of the half-cycle amplitude over the run after the first two
    !! crossings, and the final error of the smoothed one.
    REAL(wp), INTENT(IN) :: dt, y0
    INTEGER, INTENT(IN) :: nstep
    REAL(wp), INTENT(OUT) :: err, smooth_err
    REAL(wp), PARAMETER :: om = 2.0_wp*PI/1.3_wp, ph = 0.4_wp
    REAL(wp) :: st(4), y_half
    INTEGER :: i, ncross
    st = 0.0_wp
    st(4) = om*y0*COS(ph)
    err = 0.0_wp
    ncross = 0
    y_half = 0.0_wp
    DO i = 1, nstep
      CALL CD_VIV_Amplitude_Update(st, om*y0*COS(om*REAL(i, wp)*dt + ph), dt)
      IF (st(2) /= y_half) THEN
        ncross = ncross + 1
        y_half = st(2)
        IF (ncross > 1) err = MAX(err, ABS(st(2)/y0 - 1.0_wp))
      END IF
    END DO
    smooth_err = ABS(st(3)/y0 - 1.0_wp)
    CALL require(ncross > 10, 'K1: the sinusoid crosses zero')
  END SUBROUTINE run_sinusoid

  SUBROUTINE check_sinusoid()
    REAL(wp) :: e1, e2, s1, s2
    CALL run_sinusoid(1.0e-3_wp, 13000, 0.07_wp, e1, s1)
    CALL run_sinusoid(0.5e-3_wp, 26000, 0.07_wp, e2, s2)
    WRITE (*, '(A,2ES11.3,A,F6.2)') 'K1: half-cycle error at dt, dt/2: ', e1, e2, '; ratio ', e1/e2
    CALL require(e1 <= 3.0e-6_wp, 'K1: the half-cycle amplitude is y0 within 3e-6')
    CALL require(e1/e2 > 3.0_wp .AND. e1/e2 < 5.0_wp, 'K1: the error is second order in the step')
    CALL require(s1 <= 1.0e-5_wp .AND. s2 <= 1.0e-5_wp, 'K1: the smoothed amplitude settles on y0')
    CALL check_smoothing()
  END SUBROUTINE check_sinusoid

  SUBROUTINE check_smoothing()
    !! Eq. 9 per step: y0s_n+1 = y0s_n + k_a dt (y0 - y0s_n) with the half-cycle amplitude held.
    REAL(wp) :: st(4), expect
    INTEGER :: i
    st = [0.0_wp, 0.5_wp, 0.0_wp, 1.0_wp]
    expect = 0.0_wp
    DO i = 1, 400
      CALL CD_VIV_Amplitude_Update(st, 1.0_wp, 1.0e-3_wp)
      expect = expect + CD_VIV_AMP_KDT*(0.5_wp - expect)
    END DO
    CALL require(st(3) == expect, 'K1: the smoothing is Eq. 9 step by step')
    CALL require(ABS(st(3) - 0.5_wp) <= 1.0e-9_wp, 'K1: and reaches the held amplitude')
    CALL require(st(2) == 0.5_wp, 'K1: no crossing, the half-cycle amplitude is held')
  END SUBROUTINE check_smoothing

  SUBROUTINE check_touch()
    REAL(wp) :: st(4)
    st = [0.0_wp, 0.25_wp, 0.25_wp, 1.0_wp]
    CALL CD_VIV_Amplitude_Update(st, 0.0_wp, 0.1_wp)
    CALL CD_VIV_Amplitude_Update(st, 1.0_wp, 0.1_wp)
    CALL require(st(2) == 0.25_wp, 'K2: touching zero without a sign change ends no half cycle')
    CALL require(ABS(st(1) - 0.1_wp) <= 1.0e-15_wp, 'K2: the travel keeps accumulating')
    CALL CD_VIV_Amplitude_Update(st, -1.0_wp, 0.1_wp)
    ! the crossing halfway through the step: travel 0.1 + 0.025 before it, half of it the amplitude
    CALL require(ABS(st(2) - 0.0625_wp) <= 1.0e-15_wp, 'K2: a sign change ends the half cycle at the crossing')
    CALL require(ABS(st(1) + 0.025_wp) <= 1.0e-15_wp, 'K2: the rest of the step starts the next one, signed')
    CALL check_zero_samples()
  END SUBROUTINE check_touch

  SUBROUTINE check_zero_samples()
    !! A triangular velocity sampled on its zeros, 0 1 2 1 0 -1 -2 -1 0 ... with dt 0.5: every
    !! reversal passes through a sample at zero, and each half cycle travels 2 (exact for the
    !! trapezoidal rule on a piecewise-linear velocity), so y0 = 1; a touch (1 0 1) ends none.
    REAL(wp), PARAMETER :: wave(8) = [1.0_wp, 2.0_wp, 1.0_wp, 0.0_wp, -1.0_wp, -2.0_wp, -1.0_wp, 0.0_wp]
    REAL(wp) :: st(4), s_before
    INTEGER :: i, n_half
    LOGICAL :: all_one
    st = 0.0_wp
    n_half = 0
    all_one = .TRUE.
    DO i = 1, 8*5
      s_before = st(1)
      CALL CD_VIV_Amplitude_Update(st, wave(MOD(i - 1, 8) + 1), 0.5_wp)
      ! a half cycle ends where the signed travel changes sign
      IF (s_before*st(1) < 0.0_wp) THEN
        n_half = n_half + 1
        all_one = all_one .AND. st(2) == 1.0_wp
      END IF
    END DO
    CALL require(n_half == 9 .AND. all_one, 'K2: each reversal through a sample at zero ends a half cycle, y0 = 1')
    CALL require(st(1) == -2.0_wp .AND. st(4) == 0.0_wp, 'K2: the run ends on a zero sample, its half cycle open')
    st = [0.75_wp, 0.5_wp, 0.5_wp, 0.0_wp]
    CALL CD_VIV_Amplitude_Update(st, 1.0_wp, 0.5_wp)
    CALL require(st(2) == 0.5_wp .AND. st(1) == 1.0_wp, 'K2: a touch through a zero sample (+, 0, +) ends none')
    st = [0.75_wp, 0.5_wp, 0.5_wp, 0.0_wp]
    CALL CD_VIV_Amplitude_Update(st, -1.0_wp, 0.5_wp)
    CALL require(st(2) == 0.375_wp .AND. st(1) == -0.25_wp, 'K2: a reversal through it (+, 0, -) ends one there')
  END SUBROUTINE check_zero_samples

  SUBROUTINE check_growth()
    !! y = a t sin(om t): y' = a sin(om t) + a t om cos(om t). Between crossings of y' the travel
    !! is |y(t2) - y(t1)|, half of which is the mean of the two extreme amplitudes.
    REAL(wp), PARAMETER :: om = 2.0_wp*PI, a = 0.01_wp, dt = 1.0e-4_wp
    REAL(wp) :: st(4), y_half, t, yd, t_prev_cross, y_prev_ext, err
    INTEGER :: i, ncross
    st = 0.0_wp
    st(4) = 0.0_wp
    y_half = 0.0_wp
    ncross = 0
    err = 0.0_wp
    t_prev_cross = 0.0_wp
    y_prev_ext = 0.0_wp
    DO i = 1, 100000
      t = REAL(i, wp)*dt
      yd = a*SIN(om*t) + a*t*om*COS(om*t)
      CALL CD_VIV_Amplitude_Update(st, yd, dt)
      IF (st(2) /= y_half) THEN
        ncross = ncross + 1
        y_half = st(2)
        ! the extreme just passed is at about the crossing time; its |y| = a t
        IF (ncross > 10) err = MAX(err, ABS(st(2)/(0.5_wp*(a*t + y_prev_ext)) - 1.0_wp))
        y_prev_ext = a*t
      END IF
    END DO
    WRITE (*, '(A,ES11.3)') 'K3: growing oscillation, largest error of the half-cycle mean amplitude: ', err
    CALL require(ncross > 15 .AND. err <= 1.0e-3_wp, &
                 'K3: each half-cycle value is the mean amplitude of the half cycle')
  END SUBROUTINE check_growth

  SUBROUTINE check_table()
    REAL(wp) :: xs(3), ys(3), nan, y
    xs = [0.0_wp, 0.4_wp, 0.8_wp]
    ys = [0.5_wp, 1.0_wp, 0.0_wp]
    CALL require(CD_VIV_Table(-1.0_wp, xs, ys) == 0.5_wp, 'K4: held below the first row')
    CALL require(CD_VIV_Table(0.0_wp, xs, ys) == 0.5_wp, 'K4: the first row')
    CALL require(ABS(CD_VIV_Table(0.2_wp, xs, ys) - 0.75_wp) <= 1.0e-15_wp, 'K4: linear between rows')
    CALL require(CD_VIV_Table(0.4_wp, xs, ys) == 1.0_wp, 'K4: an interior row')
    CALL require(ABS(CD_VIV_Table(0.7_wp, xs, ys) - 0.25_wp) <= 1.0e-15_wp, 'K4: linear on the falling side')
    CALL require(CD_VIV_Table(5.0_wp, xs, ys) == 0.0_wp, 'K4: held beyond the last row')
    CALL require(CD_VIV_Table(3.0_wp, xs(1:1), ys(1:1)) == 0.5_wp, 'K4: a one-row table is constant')
    ! A NaN amplitude (a diverged run) stays inside the table: NaN back, or the one row.
    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    y = CD_VIV_Table(nan, xs, ys)
    CALL require(IEEE_IS_NAN(y), 'K4: a NaN amplitude gives NaN')
    CALL require(CD_VIV_Table(nan, xs(1:1), ys(1:1)) == 0.5_wp, 'K4: a one-row table is constant at NaN')
  END SUBROUTINE check_table

  SUBROUTINE check_tables_validated()
    REAL(wp) :: amp(3, CD_VIV_AMP_MAXP, 2), cl(2)
    INTEGER :: n(2)
    CHARACTER(160) :: msg
    CALL good(amp, n, cl)
    CALL CD_VIV_Check_Amp_Table(amp, n, cl, msg)
    CALL require(LEN_TRIM(msg) == 0, 'K5: a good table passes: '//TRIM(msg))
    CALL good(amp, n, cl)
    amp(1, 2, 1) = 0.0_wp
    CALL bad(amp, n, cl, 'increase strictly', 'a repeated A/D')
    CALL good(amp, n, cl)
    amp(1, 1, 1) = -0.1_wp
    CALL bad(amp, n, cl, 'A/D >= 0', 'a negative A/D')
    CALL good(amp, n, cl)
    amp(2, 2, 1) = -0.1_wp
    CALL bad(amp, n, cl, 'C_v >= 0', 'a negative C_v')
    CALL good(amp, n, cl)
    amp(3, 1, 1) = 0.0_wp
    CALL bad(amp, n, cl, 'drag factor > 0', 'a zero drag factor')
    CALL good(amp, n, cl)
    cl(1) = 0.0_wp
    CALL bad(amp, n, cl, 'needs the cross-flow lift', 'a table on an element without lift')
    CALL good(amp, n, cl)
    n(1) = CD_VIV_AMP_MAXP + 1
    CALL bad(amp, n, cl, '0 to 32 rows', 'too many rows')
    CALL good(amp, n, cl)
    CALL CD_VIV_Check_Amp_Table(amp(:, 1:4, :), n, cl, msg)
    CALL require(INDEX(msg, 'shapes') > 0, 'K5: a table of the wrong shape is refused: '//TRIM(msg))
  END SUBROUTINE check_tables_validated

  SUBROUTINE good(amp, n, cl)
    REAL(wp), INTENT(OUT) :: amp(:, :, :), cl(:)
    INTEGER, INTENT(OUT) :: n(:)
    amp = 0.0_wp
    amp(:, 1, 1) = [0.0_wp, 0.8_wp, 1.0_wp]
    amp(:, 2, 1) = [0.5_wp, 0.0_wp, 2.0_wp]
    n = [2, 0]
    cl = [0.8_wp, 0.0_wp]
  END SUBROUTINE good

  SUBROUTINE bad(amp, n, cl, needle, what)
    REAL(wp), INTENT(IN) :: amp(:, :, :), cl(:)
    INTEGER, INTENT(IN) :: n(:)
    CHARACTER(*), INTENT(IN) :: needle, what
    CHARACTER(160) :: msg
    CALL CD_VIV_Check_Amp_Table(amp, n, cl, msg)
    CALL require(INDEX(msg, needle) > 0, 'K5: '//what//' is refused: '//TRIM(msg))
  END SUBROUTINE bad

END PROGRAM test_viv_amplitude_kernel
