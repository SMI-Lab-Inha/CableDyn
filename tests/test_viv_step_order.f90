! File: tests/test_viv_step_order.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_step_order
  !! Time-step convergence of a locked-in VIV response, the line and the synchronisation model
  !! coupled: a weightless taut line (20 m, 10 elements, D 0.1 m, 8.1 kg/m, EA 5e7 N, 3.8 kN; the
  !! finite-EI cable with EI 1e3 N m^2), pinned at both ends in a 0.3 m/s cross current, Cl 0.8,
  !! cF 0.18, dF 0.08, Cd 1.2, running RMS window 0.25 s. It locks in at about 0.65 Hz with a
  !! mid-span A/D of 0.57 (0.55 on the cable). The lock-in frequency is the mean shedding rate of
  !! the mid-span node over the second half of 200 nominal cycles D/(cF U) = 1.85 s, at n steps per
  !! cycle.
  !!
  !! The phase advances by the implicit mid-point rule on the mid-step motion, the velocity and
  !! acceleration the step means from the converged positions and velocities, and its lift is
  !! predicted from the last step's mid-point rate: the
  !! coupled response is second order in the step at any rho_inf. Sampling the committed
  !! generalised-alpha acceleration instead made it first order (that acceleration is the line's
  !! at t_n + (alpha_m - alpha_f) dt) and, at rho_inf = 1, where it carries an undamped alternating
  !! mode, it converged to a frequency 3.4 % low (over 400 cycles 3.9403 rad/s at 3200 steps per
  !! cycle against 4.0797). On that code every check below fails: the ratios were 2.44, 2.27, 2.18
  !! and 2.76, and the rho_inf = 0.4 frequency 3.1 % off the rho_inf = 1 one. Tolerances from a
  !! gfortran Release run, the measured value with each:
  !!
  !! S1 EI = 0 line, rho_inf = 1: the frequency changes from 100 to 200 and from 200 to 400 steps
  !!    per cycle in a ratio in [3, 5.5] (second order 4; measured 3.99), and at 200 steps per cycle
  !!    it is within 2e-4 of the 400-step one (measured 5.3e-5);
  !! S2 the same line at the deck default rho_inf = 0.4: second order the same way (measured 3.67),
  !!    and its frequency at 200 steps per cycle within 1e-4 of the rho_inf = 1 one (measured
  !!    4.3e-5): the numerical damping no longer moves the lock-in;
  !! S3 the wake oscillator on the EI = 0 line (St 0.2, eps 0.3, A 12, C_L0 0.3, rho_inf 0.8): the
  !!    mid-span frequency from upward zero crossings, second order the same way (measured 3.90);
  !! S4 the finite-EI cable, rho_inf = 0.8, 200 to 800 steps per cycle: second order the same way
  !!    (measured 4.39).
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_VIV, ONLY: CD_VIV_TWO_PI
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Set_Drag, CD_HermiteCable_Dyn_Set_VIV, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_End, CD_HCDYN_OK
  USE CableDyn_Dynamic, ONLY: GenAlphaConfig
  USE CableDyn_Model, ONLY: CD_ModelType, CD_Init_Model, CD_Set_Model_VIV, CD_Step_Model, CD_End_Model, &
                            CD_Recompute_Model_Acceleration, CD_MODEL_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: RHOW = 1025.0_wp, D = 0.1_wp, CF = 0.18_wp, DF = 0.08_wp, U = 0.3_wp
  INTEGER, PARAMETER :: NCYC = 200
  INTEGER :: nfail = 0
  REAL(wp) :: r(3), r04(3), ratio

  ! S1
  CALL sweep(.FALSE., 1.0_wp, .FALSE., [100, 200, 400], r)
  ratio = (r(2) - r(1))/(r(3) - r(2))
  WRITE (*, '(A,3F14.9,A,F7.3)') '  S1 EI = 0 line, rho_inf 1, 100/200/400 steps: ', r, ' rad/s, ratio ', ratio
  CALL require(ratio >= 3.0_wp .AND. ratio <= 5.5_wp, 'S1: the lock-in frequency is second order in the step')
  CALL require(ABS(r(2) - r(3)) <= 2.0e-4_wp*r(3), 'S1: at 200 steps per cycle it is within 2e-4 of 400')
  ! S2
  CALL sweep(.FALSE., 0.4_wp, .FALSE., [100, 200, 400], r04)
  ratio = (r04(2) - r04(1))/(r04(3) - r04(2))
  WRITE (*, '(A,3F14.9,A,F7.3,A,ES10.2)') '  S2 rho_inf 0.4: ', r04, ' rad/s, ratio ', ratio, &
    ', against rho_inf 1 at 200 steps ', ABS(r04(2) - r(2))/r(2)
  CALL require(ratio >= 3.0_wp .AND. ratio <= 5.5_wp, 'S2: second order at rho_inf 0.4')
  CALL require(ABS(r04(2) - r(2)) <= 1.0e-4_wp*r(2), 'S2: the lock-in does not depend on rho_inf')
  ! S3
  CALL sweep(.FALSE., 0.8_wp, .TRUE., [100, 200, 400], r)
  ratio = (r(2) - r(1))/(r(3) - r(2))
  WRITE (*, '(A,3F14.9,A,F7.3)') '  S3 wake oscillator, crossing frequency: ', r, ' rad/s, ratio ', ratio
  CALL require(ratio >= 3.0_wp .AND. ratio <= 5.5_wp, 'S3: the wake response is second order in the step')
  ! S4
  CALL sweep(.TRUE., 0.8_wp, .FALSE., [200, 400, 800], r)
  ratio = (r(2) - r(1))/(r(3) - r(2))
  WRITE (*, '(A,3F14.9,A,F7.3)') '  S4 finite-EI cable, rho_inf 0.8, 200/400/800 steps: ', r, ' rad/s, ratio ', ratio
  CALL require(ratio >= 3.0_wp .AND. ratio <= 5.5_wp, 'S4: the cable lock-in is second order in the step')

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: coupled VIV response second order in the time step'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE sweep(finite_ei, rho_inf, wake, nper, freq)
    !! The lock-in frequency at each of the steps per cycle nper.
    LOGICAL, INTENT(IN) :: finite_ei, wake
    REAL(wp), INTENT(IN) :: rho_inf
    INTEGER, INTENT(IN) :: nper(3)
    REAL(wp), INTENT(OUT) :: freq(3)
    INTEGER :: k
    DO k = 1, 3
      freq(k) = lockin(finite_ei, rho_inf, wake, nper(k))
    END DO
  END SUBROUTINE sweep

  REAL(wp) FUNCTION lockin(finite_ei, rho_inf, wake, nper) RESULT(freq)
    !! The mean shedding rate of the mid-span node over the second half of NCYC cycles (with the
    !! wake oscillator, the mean frequency of the upward zero crossings of its cross-flow
    !! displacement there), 0 on a failed set-up or step.
    LOGICAL, INTENT(IN) :: finite_ei, wake
    REAL(wp), INTENT(IN) :: rho_inf
    INTEGER, INTENT(IN) :: nper
    INTEGER, PARAMETER :: NE = 10, NN = NE + 1, MIDN = NE/2 + 1
    REAL(wp), PARAMETER :: SPAN = 20.0_wp, Z0 = -50.0_wp, EA = 5.0e7_wp, MLIN = 8.1_wp, TENS = 3.8e3_wp
    REAL(wp), PARAMETER :: TW = 0.25_wp
    TYPE(CD_HermiteCableDynType) :: hm
    TYPE(CD_ModelType) :: em0
    TYPE(GenAlphaConfig) :: cfg
    REAL(wp) :: seed(6*NN), l0(NE), eav(NE), eiv(NE), rhoa(NE), wv(NE), hd(NE), hn(NE), ht(NE), wk(3, NE)
    REAL(wp) :: cl(NE), cfv(NE), dfv(NE), dt, q0(3*NN), v0(3*NN), fext(3*NN), fluid(3, NN), wl(NN)
    REAL(wp) :: z, zold, ph, phold, phu, ph0, tc1, tcn, t
    INTEGER :: conn(2, NE), i, n, nstep, n0, es, n_iter, ncross
    CHARACTER(1024) :: em
    LOGICAL :: conv, stalled, ok
    freq = 0.0_wp
    dt = D/(CF*U)/REAL(nper, wp)
    nstep = nper*NCYC
    n0 = nstep/2
    cl = MERGE(0.3_wp, 0.8_wp, wake)
    cfv = CF
    dfv = DF
    wk(1, :) = 0.2_wp
    wk(2, :) = 0.3_wp
    wk(3, :) = 12.0_wp
    l0 = SPAN/REAL(NE, wp)/(1.0_wp + TENS/EA)
    eav = EA
    rhoa = MLIN
    IF (finite_ei) THEN
      seed = 0.0_wp
      DO i = 1, NN
        seed(6*(i - 1) + 1) = SPAN*REAL(i - 1, wp)/REAL(NE, wp)
        seed(6*(i - 1) + 3) = Z0
        seed(6*(i - 1) + 4) = 1.0_wp
      END DO
      eiv = 1.0e3_wp
      wv = 0.0_wp
      CALL CD_HermiteCable_Dyn_Init(hm, l0, eav, eiv, rhoa, wv, seed, [1, 2, 3, 6*NE + 1, 6*NE + 2, 6*NE + 3], &
                                    -1.0e4_wp, 0.0_wp, rho_inf, es, em)
      ok = es == CD_HCDYN_OK
      IF (ok) THEN
        hd = D
        hn = 1.2_wp
        ht = 0.0_wp
        CALL CD_HermiteCable_Dyn_Set_Drag(hm, RHOW, hd, hn, ht, 0.0_wp, [0.0_wp, U, 0.0_wp], es, em)
        ok = es == CD_HCDYN_OK
      END IF
      IF (ok) THEN
        IF (wake) THEN
          CALL CD_HermiteCable_Dyn_Set_VIV(hm, cl, cfv, dfv, TW, es, em, wake=wk)
        ELSE
          CALL CD_HermiteCable_Dyn_Set_VIV(hm, cl, cfv, dfv, TW, es, em)
        END IF
        ok = es == CD_HCDYN_OK
      END IF
      CALL require(ok, 'set-up of the finite-EI cable: '//TRIM(em))
      IF (.NOT. ok) RETURN
      phold = hm%viv_phi(MIDN)
    ELSE
      DO i = 1, NE
        conn(:, i) = [i, i + 1]
      END DO
      q0 = 0.0_wp
      DO i = 1, NN
        q0(3*i - 2) = SPAN*REAL(i - 1, wp)/REAL(NE, wp)
        q0(3*i) = Z0
      END DO
      v0 = 0.0_wp
      fext = 0.0_wp
      fluid = 0.0_wp
      fluid(2, :) = U
      wl = 0.0_wp
      cfg%rho_inf = rho_inf
      CALL CD_Init_Model(em0, q0, v0, conn, l0, eav, rhoa, .TRUE., fext, [1, 2, 3, 3*NN - 2, 3*NN - 1, 3*NN], cfg, &
                         es, em, fluid_velocity=fluid, drag_waterline_z=wl, drag_rho=RHOW, drag_diameter=D, &
                         drag_cdn=1.2_wp, drag_cdt=0.0_wp)
      ok = es == CD_MODEL_OK
      IF (ok) THEN
        hd = D
        IF (wake) THEN
          CALL CD_Set_Model_VIV(em0, RHOW, cl, hd, cfv, dfv, TW, es, em, wake=wk)
        ELSE
          CALL CD_Set_Model_VIV(em0, RHOW, cl, hd, cfv, dfv, TW, es, em)
        END IF
        ok = es == CD_MODEL_OK
      END IF
      IF (ok) THEN
        CALL CD_Recompute_Model_Acceleration(em0, es, em)
        ok = es == CD_MODEL_OK
      END IF
      CALL require(ok, 'set-up of the EI = 0 line: '//TRIM(em))
      IF (.NOT. ok) RETURN
      phold = em0%viv_phi(MIDN)
    END IF
    phu = phold
    ph0 = 0.0_wp
    zold = 0.0_wp
    ncross = 0
    tc1 = 0.0_wp
    tcn = 0.0_wp
    DO n = 1, nstep
      IF (finite_ei) THEN
        CALL CD_HermiteCable_Dyn_Step(hm, dt, 50, 1.0e-8_wp, es, em)
        ok = es == CD_HCDYN_OK
        ph = hm%viv_phi(MIDN)
        z = hm%q(6*(MIDN - 1) + 3) - Z0
      ELSE
        CALL CD_Step_Model(em0, dt, conv, stalled, n_iter, es, em)
        ok = es == CD_MODEL_OK .AND. conv
        ph = em0%viv_phi(MIDN)
        z = em0%q(3*MIDN) - Z0
      END IF
      CALL require(ok, 'step: '//TRIM(em))
      IF (.NOT. ok) EXIT
      phu = phu + (ph - phold) - CD_VIV_TWO_PI*ANINT((ph - phold)/CD_VIV_TWO_PI)
      phold = ph
      t = REAL(n, wp)*dt
      IF (n == n0) ph0 = phu
      IF (n > n0 .AND. zold <= 0.0_wp .AND. z > 0.0_wp) THEN
        IF (ncross == 0) tc1 = t - dt*z/(z - zold)
        tcn = t - dt*z/(z - zold)
        ncross = ncross + 1
      END IF
      zold = z
    END DO
    IF (ok) THEN
      IF (wake) THEN
        IF (ncross > 1) freq = CD_VIV_TWO_PI*REAL(ncross - 1, wp)/(tcn - tc1)
      ELSE
        freq = (phu - ph0)/(REAL(nstep - n0, wp)*dt)
      END IF
    END IF
    IF (finite_ei) THEN
      CALL CD_HermiteCable_Dyn_End(hm)
    ELSE
      CALL CD_End_Model(em0, es, em)
    END IF
  END FUNCTION lockin
END PROGRAM test_viv_step_order
