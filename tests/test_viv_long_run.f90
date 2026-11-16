! File: tests/test_viv_long_run.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_long_run
  !! Long runs of the cross-flow VIV synchronisation model (CableDyn_VIV): hundreds of cycles, as a
  !! VIV study runs them, at 100 steps per cycle. The shedding phase advances by Heun's rule
  !! (viv_step_begin / viv_step_end of both line kinds), so its error must grow linearly in time
  !! and fall with dt^2; the load powers must close the energy balance over the whole run; and a
  !! locked-in response must be stationary. Every tolerance below comes from a dt sweep (25, 50,
  !! 100, 200, 400 steps per cycle, gfortran Release) and is stated with the measured value.
  !!
  !! P  Prescribed motion, the finite-EI cable (10 elements, 20 m, every node's position held and
  !!    moved): all nodes move across a uniform flow U = 0.3 m/s (along y) as a rigid body,
  !!    z = z0 - A cos(om t), A = 1e-6 m (the phase does not depend on A; a small A keeps the
  !!    relative flow and so U and the cross-flow axis fixed to (A om/U)^2 ~ 1e-10). The running
  !!    mean squares start at those of the motion and average over 1000 s, so the motion phase is
  !!    phi_y = om t - pi/2 to the averaging ripple (about 1e-4 rad), and psi = phi_y - phi obeys
  !!    Adler's equation d psi/dt = a - b sin(psi), a = om - Om cF, b = Om dF, Om = 2 pi U/D,
  !!    whose solution is closed form: tan(psi/2) = (b + s tan(theta/2))/a, theta = theta0 + s t,
  !!    s = sqrt(a^2 - b^2) outside the band (a > b), and psi -> asin(a/b) inside it.
  !!    P1 outside the band (om 5.5 rad/s, slip period 2 pi/s = 4.3 s), 400 motion cycles: the
  !!       largest phase error over the 11 nodes against the closed form grows linearly: error(T)/
  !!       error(T/2) in [1.8, 2.2] (linear 2, quadratic 4; measured 2.014 and 1.988 at 50 and 100
  !!       steps per cycle);
  !!    P2 it is second order: error at 50 over error at 100 steps per cycle in [3.6, 4.4]
  !!       (measured 4.01), and at 100 steps per cycle at most 0.05 rad after 400 cycles (measured
  !!       0.0167, i.e. 4e-5 rad per cycle);
  !!    P3 inside the band (om = Om (cF + dF/2), psi* = pi/6), 50 steps per cycle: the phase lag
  !!       stays within 1e-4 rad of psi* at T/2 and T (measured 2.9e-5 at 25, 50 and 100 steps per
  !!       cycle: the averaging-ripple floor, not the step) and does not grow (error(T) <=
  !!       1.1 error(T/2), measured 1.0004), and the mean shedding rate over the second half is om
  !!       within 1e-8 (measured 5.8e-12).
  !!
  !! C  Lock-in, both line kinds: a weightless taut line (20 m, 10 elements, D 0.1 m, 8.1 kg/m,
  !!    EA 5e7 N, 3.8 kN, the finite-EI cable EI 1e3 N m^2, no added mass so that the structural
  !!    energy is the whole stored energy), pinned at both ends in the 0.3 m/s cross current, Cl
  !!    0.8, cF 0.18, dF 0.08, Cd 1.2, running RMS window 0.25 s (the examples' value), rho_inf = 1
  !!    (no numerical damping, so that the load powers are the whole energy change). It locks in at
  !!    about 0.63 Hz with a mid-span A/D of 0.56-0.57. 400 cycles of the nominal shedding period
  !!    D/(cF U) = 1.85 s (740 s) at 100 and at 200 steps per cycle.
  !!    C1 energy: E(t) - E(0) = int (P_viv + P_drag + P_damp) dt (trapezoidal), the mismatch at
  !!       most 3e-4 of the VIV work at 100 steps per cycle (measured 2.0e-4 finite EI, 1.8e-4
  !!       EI = 0: the lift takes the phase predicted from the previous step's mid-point rate, half
  !!       a step behind the step's own, a second-order lag), second order (the ratio from 100 to
  !!       200 steps per cycle in [3, 5], measured 3.98 and 3.67), and it accumulates with the work
  !!       and not faster: the mismatch fraction at T within 1.25 times that at T/2 (measured 1.00
  !!       and 1.00);
  !!    C2 statistics over the last three quarters of the run, its two halves: the mid-span RMS
  !!       A/D within 2e-3 and the mean VIV power within 1e-2 of each other (measured 7.3e-4 and
  !!       1.8e-3 finite EI, 4.1e-4 and 9.0e-4 EI = 0), the same bounds as the change of the
  !!       second half's values from 100 to 200 steps per cycle must meet (measured A/D 7.9e-5 and
  !!       6.0e-4, power 1.1e-3 and 3.3e-4), so a drift between the halves stays below the
  !!       discretisation error of the statistic;
  !!    C3 phase: the mean shedding rate (the lock-in frequency) of the mid-span node over the two
  !!       halves agrees within 2e-3 (measured 7.4e-5 finite EI, 1.9e-4 EI = 0), so the phase error
  !!       against a finer step grows linearly (a frequency offset) and not quadratically, and the
  !!       second half's rate changes by at most 5e-3 from 100 to 200 steps per cycle (measured
  !!       1.6e-4 finite EI, 1.9e-4 EI = 0). The coupled lock-in frequency converges at second
  !!       order (test viv_step_order).
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_VIV, ONLY: CD_VIV_TWO_PI
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Set_Drag, CD_HermiteCable_Dyn_Set_VIV, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_Recompute_Acceleration, &
                                          CD_HermiteCable_Dyn_Energy, CD_HermiteCable_Dyn_Load_Power, &
                                          CD_HermiteCable_Dyn_End, CD_HCDYN_OK
  USE CableDyn_Dynamic, ONLY: GenAlphaConfig
  USE CableDyn_Model, ONLY: CD_ModelType, CD_Init_Model, CD_Set_Model_VIV, CD_Step_Model, CD_End_Model, &
                            CD_Recompute_Model_Acceleration, CD_Get_Model_Load_Power, CD_MODEL_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
  REAL(wp), PARAMETER :: RHOW = 1025.0_wp, D = 0.1_wp, CF = 0.18_wp, DF = 0.08_wp, U = 0.3_wp
  REAL(wp), PARAMETER :: OMS = CD_VIV_TWO_PI*U/D
  INTEGER, PARAMETER :: NCYC = 400

  !> The results of one lock-in run (C): the energy mismatch over the VIV work at T/2 and T, the
  !> mid-span RMS A/D and mean VIV power over the two halves of the last three quarters, and the
  !> mean shedding rate of the mid-span node over the same halves.
  TYPE :: lockin_result
    LOGICAL :: ok = .FALSE.
    REAL(wp) :: closure(2) = 0.0_wp, ad(2) = 0.0_wp, pviv(2) = 0.0_wp, rate(2) = 0.0_wp
  END TYPE lockin_result

  INTEGER :: nfail = 0

  CALL check_prescribed()
  CALL check_lockin(.TRUE., 'finite-EI cable')
  CALL check_lockin(.FALSE., 'EI = 0 line')
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: long-run VIV phase, energy and statistics'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  PURE REAL(wp) FUNCTION wrap_pi(x) RESULT(y)
    !! x mapped to [-pi, pi].
    REAL(wp), INTENT(IN) :: x
    y = x - CD_VIV_TWO_PI*ANINT(x/CD_VIV_TWO_PI)
  END FUNCTION wrap_pi

  PURE REAL(wp) FUNCTION rel(x, y) RESULT(r)
    !! The relative difference of two statistics.
    REAL(wp), INTENT(IN) :: x, y
    r = ABS(x - y)/MAX(ABS(x), ABS(y), TINY(1.0_wp))
  END FUNCTION rel

  PURE REAL(wp) FUNCTION adler_psi(t, psi0, a, b) RESULT(psi)
    !! The continuous solution of d psi/dt = a - b sin(psi), a > b > 0, from psi(0) = psi0 in
    !! (-pi, pi): tan(psi/2) = (b + s tan(theta/2))/a, theta = theta0 + s t, s = sqrt(a^2 - b^2),
    !! continued over the turns of theta.
    REAL(wp), INTENT(IN) :: t, psi0, a, b
    REAL(wp) :: s, th, thr, k
    s = SQRT(a*a - b*b)
    th = 2.0_wp*ATAN((a*TAN(0.5_wp*psi0) - b)/s) + s*t
    k = ANINT(th/CD_VIV_TWO_PI)
    thr = th - CD_VIV_TWO_PI*k
    psi = CD_VIV_TWO_PI*k + 2.0_wp*ATAN((b + s*TAN(0.5_wp*thr))/a)
  END FUNCTION adler_psi

  SUBROUTINE check_prescribed()
    REAL(wp) :: e50(2), e100(2), el(2), rate_err
    LOGICAL :: ok
    CALL prescribed(5.5_wp, 50, e50, rate_err, ok)
    IF (.NOT. ok) RETURN
    CALL prescribed(5.5_wp, 100, e100, rate_err, ok)
    IF (.NOT. ok) RETURN
    WRITE (*, '(A,2ES11.3,A,2ES11.3)') '  P1/P2 outside the band, phase error at T/2, T: 50 steps ', e50, &
      ', 100 steps ', e100
    WRITE (*, '(A,3F8.4)') '  P1/P2 growth T/(T/2) at 50, 100 steps, order 50/100: ', e50(2)/e50(1), e100(2)/e100(1), &
      e50(2)/e100(2)
    CALL require(e50(2) >= 1.8_wp*e50(1) .AND. e50(2) <= 2.2_wp*e50(1) .AND. &
                 e100(2) >= 1.8_wp*e100(1) .AND. e100(2) <= 2.2_wp*e100(1), &
                 'P1: the phase error grows linearly in time')
    CALL require(e50(2) >= 3.6_wp*e100(2) .AND. e50(2) <= 4.4_wp*e100(2), 'P2: the phase error is second order')
    CALL require(e100(2) <= 0.05_wp, 'P2: at most 0.05 rad after 400 cycles at 100 steps per cycle')
    CALL prescribed(OMS*(CF + 0.5_wp*DF), 50, el, rate_err, ok)
    IF (.NOT. ok) RETURN
    WRITE (*, '(A,2ES11.3,A,ES11.3)') '  P3 inside the band, phase-lag error at T/2, T: ', el, &
      ', mean rate error ', rate_err
    CALL require(MAXVAL(el) <= 1.0e-4_wp .AND. el(2) <= 1.1_wp*el(1), 'P3: the locked phase lag holds')
    CALL require(rate_err <= 1.0e-8_wp, 'P3: the locked shedding rate is the motion frequency')
  END SUBROUTINE check_prescribed

  SUBROUTINE prescribed(om, nper, err, rate_err, ok)
    !! Every node of the straight cable moved across the flow as a rigid body at om, nper steps per
    !! motion cycle over NCYC cycles: err the largest phase error over the nodes at T/2 and T
    !! (outside the band against the closed form, inside it the phase lag against psi*), rate_err
    !! the relative error of the mean shedding rate over the second half inside the band.
    REAL(wp), INTENT(IN) :: om
    INTEGER, INTENT(IN) :: nper
    REAL(wp), INTENT(OUT) :: err(2), rate_err
    LOGICAL, INTENT(OUT) :: ok
    INTEGER, PARAMETER :: NE = 10, NN = NE + 1
    REAL(wp), PARAMETER :: SPAN = 20.0_wp, AMP = 1.0e-6_wp, Z0 = -50.0_wp, TW = 1.0e3_wp
    TYPE(CD_HermiteCableDynType) :: m
    REAL(wp) :: seed(6*NN), l0(NE), eav(NE), eiv(NE), rhoa(NE), wv(NE), hd(NE), hn(NE), ht(NE)
    REAL(wp) :: cl(NE), cfv(NE), dfv(NE), pq(3*NN), pv(3*NN), pa(3*NN), dt, t, a, b, phi0(NN), psi0(NN)
    REAL(wp) :: phu(NN), phold(NN), ex, e, ph_half(NN), psi_star
    INTEGER :: pdofs(3*NN), i, n, nstep, es
    CHARACTER(1024) :: em
    LOGICAL :: locked
    err = 0.0_wp; rate_err = 0.0_wp
    seed = 0.0_wp
    DO i = 1, NN
      seed(6*(i - 1) + 1) = SPAN*REAL(i - 1, wp)/REAL(NE, wp)
      seed(6*(i - 1) + 3) = Z0 - AMP
      seed(6*(i - 1) + 4) = 1.0_wp
      pdofs(3*i - 2:3*i) = [6*(i - 1) + 1, 6*(i - 1) + 2, 6*(i - 1) + 3]
    END DO
    ! a stiff, unstretched cable: the slopes (free) barely turn under the tiny lift and the drag
    l0 = SPAN/REAL(NE, wp); eav = 5.0e7_wp; eiv = 1.0e6_wp; rhoa = 8.1_wp; wv = 0.0_wp
    CALL CD_HermiteCable_Dyn_Init(m, l0, eav, eiv, rhoa, wv, seed, pdofs, -1.0e4_wp, 0.0_wp, 1.0_wp, es, em)
    ok = es == CD_HCDYN_OK
    CALL require(ok, 'P: init: '//TRIM(em))
    IF (.NOT. ok) RETURN
    hd = D; hn = 1.2_wp; ht = 0.0_wp
    CALL CD_HermiteCable_Dyn_Set_Drag(m, RHOW, hd, hn, ht, 0.0_wp, [0.0_wp, U, 0.0_wp], es, em)
    ok = es == CD_HCDYN_OK
    CALL require(ok, 'P: drag: '//TRIM(em))
    IF (.NOT. ok) RETURN
    cl = 1.0e-3_wp; cfv = CF; dfv = DF
    CALL CD_HermiteCable_Dyn_Set_VIV(m, cl, cfv, dfv, TW, es, em)
    ok = es == CD_HCDYN_OK
    CALL require(ok, 'P: VIV: '//TRIM(em))
    IF (.NOT. ok) RETURN
    CALL kin(0.0_wp, om, seed, AMP, Z0, NN, pq, pv, pa)
    CALL CD_HermiteCable_Dyn_Recompute_Acceleration(m, es, em, pres_dofs=pdofs, pres_a=pa)
    ok = es == CD_HCDYN_OK
    CALL require(ok, 'P: start acceleration: '//TRIM(em))
    IF (.NOT. ok) RETURN
    ! the running mean squares start at those of the motion, y' = A om sin(om t)
    m%viv_ms(1, :) = 0.5_wp*(AMP*om)**2
    m%viv_ms(2, :) = 0.5_wp*(AMP*om*om)**2
    a = om - OMS*CF
    b = OMS*DF
    locked = ABS(a) < b
    psi_star = ASIN(MAX(-1.0_wp, MIN(1.0_wp, a/b)))
    phi0 = m%viv_phi
    phu = phi0
    phold = phi0
    DO i = 1, NN
      psi0(i) = wrap_pi(-0.5_wp*PI - phi0(i))
    END DO
    dt = CD_VIV_TWO_PI/om/REAL(nper, wp)
    nstep = nper*NCYC
    DO n = 1, nstep
      t = REAL(n, wp)*dt
      CALL kin(t, om, seed, AMP, Z0, NN, pq, pv, pa)
      CALL CD_HermiteCable_Dyn_Step(m, dt, 50, 1.0e-10_wp, es, em, pres_dofs=pdofs, pres_q=pq, pres_v=pv, pres_a=pa)
      ok = es == CD_HCDYN_OK
      CALL require(ok, 'P: step: '//TRIM(em))
      IF (.NOT. ok) EXIT
      DO i = 1, NN
        phu(i) = phu(i) + wrap_pi(m%viv_phi(i) - phold(i))
        phold(i) = m%viv_phi(i)
      END DO
      IF (n /= nstep/2 .AND. n /= nstep) CYCLE
      e = 0.0_wp
      DO i = 1, NN
        IF (locked) THEN
          ex = wrap_pi(phu(i) - (om*t - 0.5_wp*PI - psi_star))
        ELSE
          ex = phu(i) - (phi0(i) + om*t - (adler_psi(t, psi0(i), a, b) - psi0(i)))
        END IF
        e = MAX(e, ABS(ex))
      END DO
      IF (n == nstep/2) THEN
        err(1) = e
        ph_half = phu
      ELSE
        err(2) = e
      END IF
    END DO
    IF (ok .AND. locked) rate_err = ABS(SUM(phu - ph_half)/REAL(NN, wp)/(0.5_wp*REAL(nstep, wp)*dt) - om)/om
    CALL CD_HermiteCable_Dyn_End(m)
  END SUBROUTINE prescribed

  SUBROUTINE kin(tt, om, seed, amp, z0, nn, pq, pv, pa)
    !! The prescribed rigid-body motion of every node at time tt.
    REAL(wp), INTENT(IN) :: tt, om, seed(:), amp, z0
    INTEGER, INTENT(IN) :: nn
    REAL(wp), INTENT(OUT) :: pq(:), pv(:), pa(:)
    INTEGER :: j
    DO j = 1, nn
      pq(3*j - 2:3*j) = [seed(6*(j - 1) + 1), 0.0_wp, z0 - amp*COS(om*tt)]
      pv(3*j - 2:3*j) = [0.0_wp, 0.0_wp, amp*om*SIN(om*tt)]
      pa(3*j - 2:3*j) = [0.0_wp, 0.0_wp, amp*om*om*COS(om*tt)]
    END DO
  END SUBROUTINE kin

  SUBROUTINE check_lockin(finite_ei, label)
    LOGICAL, INTENT(IN) :: finite_ei
    CHARACTER(*), INTENT(IN) :: label
    TYPE(lockin_result) :: r1, r2
    CALL lockin(finite_ei, 100, r1)
    IF (.NOT. r1%ok) RETURN
    CALL lockin(finite_ei, 200, r2)
    IF (.NOT. r2%ok) RETURN
    WRITE (*, '(3A,2ES11.3,A,2ES11.3)') '  C1 ', label, ': energy mismatch / VIV work at T/2, T, 100 steps ', &
      r1%closure, ', 200 steps ', r2%closure
    WRITE (*, '(3A,4F11.6)') '  C2 ', label, ': RMS A/D and VIV power [W] over the halves ', r1%ad/D, r1%pviv
    WRITE (*, '(3A,2F12.8,A,F12.8)') '  C3 ', label, ': mean shedding rate [rad/s] over the halves ', r1%rate, &
      ', at 200 steps ', r2%rate(2)
    WRITE (*, '(3A,5ES10.2)') '  ', label, ': halves A/D, power, rate; 100 to 200 steps A/D, power ', &
      rel(r1%ad(1), r1%ad(2)), rel(r1%pviv(1), r1%pviv(2)), rel(r1%rate(1), r1%rate(2)), &
      rel(r1%ad(2), r2%ad(2)), rel(r1%pviv(2), r2%pviv(2))
    CALL require(ABS(r1%closure(2)) <= 3.0e-4_wp, 'C1 '//label//': the load powers close the energy balance')
    CALL require(ABS(r1%closure(2)) >= 3.0_wp*ABS(r2%closure(2)) .AND. &
                 ABS(r1%closure(2)) <= 5.0_wp*ABS(r2%closure(2)), 'C1 '//label//': the mismatch is second order')
    CALL require(ABS(r1%closure(2)) <= 1.25_wp*ABS(r1%closure(1)), &
                 'C1 '//label//': the mismatch grows no faster than the work')
    CALL require(rel(r1%ad(1), r1%ad(2)) <= 2.0e-3_wp .AND. rel(r1%pviv(1), r1%pviv(2)) <= 1.0e-2_wp, &
                 'C2 '//label//': the response is stationary')
    CALL require(rel(r1%ad(2), r2%ad(2)) <= 2.0e-3_wp .AND. rel(r1%pviv(2), r2%pviv(2)) <= 1.0e-2_wp, &
                 'C2 '//label//': the statistics are converged in the step')
    CALL require(rel(r1%rate(1), r1%rate(2)) <= 2.0e-3_wp, 'C3 '//label//': the lock-in frequency does not drift')
    CALL require(rel(r1%rate(2), r2%rate(2)) <= 5.0e-3_wp, 'C3 '//label//': the lock-in frequency is converged')
  END SUBROUTINE check_lockin

  SUBROUTINE lockin(finite_ei, nper, res)
    !! The taut line locked in over NCYC nominal cycles at nper steps per cycle.
    LOGICAL, INTENT(IN) :: finite_ei
    INTEGER, INTENT(IN) :: nper
    TYPE(lockin_result), INTENT(OUT) :: res
    INTEGER, PARAMETER :: NE = 10, NN = NE + 1, MIDN = NE/2 + 1
    REAL(wp), PARAMETER :: SPAN = 20.0_wp, Z0 = -50.0_wp, EA = 5.0e7_wp, MLIN = 8.1_wp, TENS = 3.8e3_wp
    REAL(wp), PARAMETER :: TW = 0.25_wp
    TYPE(CD_HermiteCableDynType) :: hm
    TYPE(CD_ModelType) :: em0
    TYPE(GenAlphaConfig) :: cfg
    REAL(wp) :: seed(6*NN), l0(NE), eav(NE), eiv(NE), rhoa(NE), wv(NE), hd(NE), hn(NE), ht(NE)
    REAL(wp) :: cl(NE), cfv(NE), dfv(NE), dt, tv, e0, e_now, p_prev(3), p_now(3), work(3)
    REAL(wp) :: q0(3*NN), v0(3*NN), fext(3*NN), fluid(3, NN), wl(NN), z, ph, phold, phu, s2(2), sp(2), cnt(2)
    REAL(wp) :: ph_mark(3)
    INTEGER :: conn(2, NE), i, n, nstep, es, n_iter, n0, n1, h
    CHARACTER(1024) :: em
    LOGICAL :: conv, stalled, ok
    tv = D/(CF*U)
    dt = tv/REAL(nper, wp)
    nstep = nper*NCYC
    ! the statistics window: the last three quarters, in two halves (n0, n1] and (n1, nstep]
    n0 = nper*(NCYC/4)
    n1 = nper*(NCYC/4 + (NCYC - NCYC/4)/2)
    cl = 0.8_wp; cfv = CF; dfv = DF
    l0 = SPAN/REAL(NE, wp)/(1.0_wp + TENS/EA)
    eav = EA; rhoa = MLIN
    IF (finite_ei) THEN
      seed = 0.0_wp
      DO i = 1, NN
        seed(6*(i - 1) + 1) = SPAN*REAL(i - 1, wp)/REAL(NE, wp)
        seed(6*(i - 1) + 3) = Z0
        seed(6*(i - 1) + 4) = 1.0_wp
      END DO
      eiv = 1.0e3_wp; wv = 0.0_wp
      CALL CD_HermiteCable_Dyn_Init(hm, l0, eav, eiv, rhoa, wv, seed, [1, 2, 3, 6*NE + 1, 6*NE + 2, 6*NE + 3], &
                                    -1.0e4_wp, 0.0_wp, 1.0_wp, es, em)
      ok = es == CD_HCDYN_OK
      IF (ok) THEN
        hd = D; hn = 1.2_wp; ht = 0.0_wp
        CALL CD_HermiteCable_Dyn_Set_Drag(hm, RHOW, hd, hn, ht, 0.0_wp, [0.0_wp, U, 0.0_wp], es, em)
        ok = es == CD_HCDYN_OK
      END IF
      IF (ok) THEN
        CALL CD_HermiteCable_Dyn_Set_VIV(hm, cl, cfv, dfv, TW, es, em)
        ok = es == CD_HCDYN_OK
      END IF
      CALL require(ok, 'C: the finite-EI cable sets up: '//TRIM(em))
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
      v0 = 0.0_wp; fext = 0.0_wp
      fluid = 0.0_wp; fluid(2, :) = U; wl = 0.0_wp
      cfg%rho_inf = 1.0_wp
      CALL CD_Init_Model(em0, q0, v0, conn, l0, eav, rhoa, .TRUE., fext, [1, 2, 3, 3*NN - 2, 3*NN - 1, 3*NN], cfg, &
                         es, em, fluid_velocity=fluid, drag_waterline_z=wl, drag_rho=RHOW, drag_diameter=D, &
                         drag_cdn=1.2_wp, drag_cdt=0.0_wp)
      ok = es == CD_MODEL_OK
      IF (ok) THEN
        hd = D
        CALL CD_Set_Model_VIV(em0, RHOW, cl, hd, cfv, dfv, TW, es, em)
        ok = es == CD_MODEL_OK
      END IF
      IF (ok) THEN
        CALL CD_Recompute_Model_Acceleration(em0, es, em)
        ok = es == CD_MODEL_OK
      END IF
      CALL require(ok, 'C: the EI = 0 line sets up: '//TRIM(em))
      IF (.NOT. ok) RETURN
      phold = em0%viv_phi(MIDN)
    END IF
    phu = phold
    CALL energy(finite_ei, hm, em0, l0, EA, MLIN, e0)
    CALL powers(finite_ei, hm, em0, NN, p_prev)
    work = 0.0_wp
    s2 = 0.0_wp; sp = 0.0_wp; cnt = 0.0_wp
    ph_mark = 0.0_wp
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
      CALL require(ok, 'C: step: '//TRIM(em))
      IF (.NOT. ok) EXIT
      phu = phu + wrap_pi(ph - phold)
      phold = ph
      CALL powers(finite_ei, hm, em0, NN, p_now)
      work = work + 0.5_wp*dt*(p_prev + p_now)
      p_prev = p_now
      IF (n == n0) ph_mark(1) = phu
      IF (n == n1) ph_mark(2) = phu
      IF (n > n0) THEN
        h = MERGE(1, 2, n <= n1)
        s2(h) = s2(h) + z*z
        sp(h) = sp(h) + p_now(1)
        cnt(h) = cnt(h) + 1.0_wp
      END IF
      IF (n == nstep/2 .OR. n == nstep) THEN
        CALL energy(finite_ei, hm, em0, l0, EA, MLIN, e_now)
        res%closure(MERGE(1, 2, n == nstep/2)) = (e_now - e0 - SUM(work))/work(1)
      END IF
    END DO
    ph_mark(3) = phu
    IF (ok) THEN
      res%ok = .TRUE.
      res%ad = SQRT(2.0_wp*s2/cnt)
      res%pviv = sp/cnt
      res%rate(1) = (ph_mark(2) - ph_mark(1))/(REAL(n1 - n0, wp)*dt)
      res%rate(2) = (ph_mark(3) - ph_mark(2))/(REAL(nstep - n1, wp)*dt)
    END IF
    IF (finite_ei) THEN
      CALL CD_HermiteCable_Dyn_End(hm)
    ELSE
      CALL CD_End_Model(em0, es, em)
    END IF
  END SUBROUTINE lockin

  SUBROUTINE energy(finite_ei, hm, em0, l0, ea, mlin, e)
    !! The line's kinetic plus elastic energy: the finite-EI cable's own (consistent mass), or for
    !! the EI = 0 line the consistent translational mass (rho_a l0/6 [[2, 1], [1, 2]]) and the
    !! tension-only axial strain energy EA (l - l0)^2/(2 l0).
    LOGICAL, INTENT(IN) :: finite_ei
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: hm
    TYPE(CD_ModelType), INTENT(IN) :: em0
    REAL(wp), INTENT(IN) :: l0(:), ea, mlin
    REAL(wp), INTENT(OUT) :: e
    INTEGER :: k, ia, ib, ess
    REAL(wp) :: l, kin_, str_
    CHARACTER(256) :: emm
    IF (finite_ei) THEN
      CALL CD_HermiteCable_Dyn_Energy(hm, kin_, str_, ess, emm)
      CALL require(ess == CD_HCDYN_OK, 'C: energy: '//TRIM(emm))
      e = kin_ + str_
    ELSE
      e = 0.0_wp
      DO k = 1, SIZE(l0)
        ia = 3*k - 2; ib = 3*k + 1
        e = e + mlin*l0(k)/6.0_wp*(DOT_PRODUCT(em0%v(ia:ia + 2), em0%v(ia:ia + 2)) + &
                                   DOT_PRODUCT(em0%v(ia:ia + 2), em0%v(ib:ib + 2)) + &
                                   DOT_PRODUCT(em0%v(ib:ib + 2), em0%v(ib:ib + 2)))
        l = NORM2(em0%q(ib:ib + 2) - em0%q(ia:ia + 2))
        IF (l > l0(k)) e = e + 0.5_wp*ea*(l - l0(k))**2/l0(k)
      END DO
    END IF
  END SUBROUTINE energy

  SUBROUTINE powers(finite_ei, hm, em0, nn, p)
    !! The line totals of the VIV, drag and damping powers.
    LOGICAL, INTENT(IN) :: finite_ei
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: hm
    TYPE(CD_ModelType), INTENT(IN) :: em0
    INTEGER, INTENT(IN) :: nn
    REAL(wp), INTENT(OUT) :: p(3)
    REAL(wp) :: pv(nn), pd(nn), pm(nn)
    INTEGER :: ess
    CHARACTER(256) :: emm
    IF (finite_ei) THEN
      CALL CD_HermiteCable_Dyn_Load_Power(hm, pv, pd, pm, ess, emm)
    ELSE
      CALL CD_Get_Model_Load_Power(em0, pv, pd, pm, ess, emm)
    END IF
    CALL require(ess == 0, 'C: load power: '//TRIM(emm))
    p = [SUM(pv), SUM(pd), SUM(pm)]
  END SUBROUTINE powers

END PROGRAM test_viv_long_run
