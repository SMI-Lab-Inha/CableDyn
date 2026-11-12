! File: tests/test_hermite_load_power.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_hermite_load_power
  !! The load powers of a finite-EI cable (CD_HermiteCable_Dyn_Load_Power) close its energy
  !! balance: on a weightless, pinned-pinned, tensioned cable with no added mass the kinetic plus
  !! elastic energy changes by the work of the VIV loads, the drag and the structural damping,
  !!   E(t) - E(0) = int (P_viv + P_drag + P_damp) dt,
  !! integrated by the trapezoidal rule over generalised-alpha steps with rho_inf = 1 (no numerical
  !! damping) of 1/200 of the period of the motion the case resolves, so the two sides agree to
  !! the steps' second-order error.
  !!   L1 the free decay of the first bending mode with the bending dashpot Bb;
  !!   L2 the free decay of the first axial mode with the axial dashpot BA (its own time step);
  !!   L3 the bending decay in still water with Morison drag and no dashpots;
  !!   L4 a cable starting straight in a 0.5 m/s current with the cross-flow, in-line and
  !!      higher-harmonic (Chh 1) VIV loads,
  !!      the drag and the bending dashpot;
  !! each within 0.5 % of the largest work; the damping power is never positive at any node, and
  !! the VIV and drag powers are zero without those loads. L5 power arrays of the wrong size are
  !! refused.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_Set_Bending_Damping, &
                                          CD_HermiteCable_Dyn_Set_Axial_Damping, CD_HermiteCable_Dyn_Set_Drag, &
                                          CD_HermiteCable_Dyn_Set_VIV, CD_HermiteCable_Dyn_Energy, &
                                          CD_HermiteCable_Dyn_Load_Power, CD_HCDYN_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
  REAL(wp), PARAMETER :: LL = 20.0_wp, EA = 1.0e8_wp, EI = 1.0e5_wp, RHOA = 10.0_wp, TENS = 1.0e4_wp
  REAL(wp), PARAMETER :: DIAM = 0.1_wp, RHOW = 1025.0_wp
  INTEGER, PARAMETER :: NE = 40, NN = NE + 1, NDOF = 6*NN
  INTEGER, PARAMETER :: BENDING = 1, AXIAL = 2, DRAG = 3, VIV = 4
  INTEGER :: nfail = 0

  CALL balance_case('L1 bending dashpot decay', BENDING)
  CALL balance_case('L2 axial dashpot decay', AXIAL)
  CALL balance_case('L3 still-water drag decay', DRAG)
  CALL balance_case('L4 VIV in a current', VIV)
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: Hermite load power energy balance'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE setup(model, a0, u0, ok)
    !! The tensioned cable along x, pinned at both ends, at rest in its first bending mode of
    !! amplitude a0 (z) and its first axial mode of amplitude u0 (x), rho_inf = 1.
    TYPE(CD_HermiteCableDynType), INTENT(OUT) :: model
    REAL(wp), INTENT(IN) :: a0, u0
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp) :: l0(NE), ea_(NE), ei_(NE), rho(NE), w(NE), q(NDOF), eps, lt, x, s, c
    INTEGER :: k, es
    CHARACTER(512) :: em
    eps = TENS/EA; lt = (1.0_wp + eps)*LL
    l0 = LL/REAL(NE, wp); ea_ = EA; ei_ = EI; rho = RHOA; w = 0.0_wp
    DO k = 1, NN
      x = lt*REAL(k - 1, wp)/REAL(NE, wp)
      s = SIN(PI*x/lt)
      c = COS(PI*x/lt)
      q(6*k - 5:6*k) = [x + u0*s, 0.0_wp, a0*s, (1.0_wp + eps) + (1.0_wp + eps)*u0*PI/lt*c, 0.0_wp, &
                        (1.0_wp + eps)*a0*PI/lt*c]
    END DO
    CALL CD_HermiteCable_Dyn_Init(model, l0, ea_, ei_, rho, w, q, [1, 2, 3, 6*NE + 1, 6*NE + 2, 6*NE + 3], &
                                  -1.0e4_wp, 0.0_wp, 1.0_wp, es, em)
    ok = es == CD_HCDYN_OK
    CALL require(ok, 'init: '//TRIM(em))
  END SUBROUTINE setup

  SUBROUTINE balance_case(label, what)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: what
    TYPE(CD_HermiteCableDynType) :: model
    REAL(wp) :: kk, period, dt, ek, es0, e0, e1, p_prev(3), p_now(3), work(3), res, scale, pm_max
    REAL(wp) :: cl(NE), cf(NE), df(NE), clx(NE), alx(NE), chh(NE)
    INTEGER :: n, k, es, it
    CHARACTER(512) :: em
    LOGICAL :: zero_viv, zero_drag, ok
    kk = PI/LL
    IF (what == AXIAL) THEN
      period = 2.0_wp*LL/SQRT(EA/RHOA)
      CALL setup(model, 0.0_wp, 1.0e-4_wp, ok)
    ELSE
      period = 2.0_wp*PI/SQRT((EI*kk**4 + TENS*kk**2)/RHOA)
      CALL setup(model, MERGE(0.0_wp, 0.05_wp, what == VIV), 0.0_wp, ok)
    END IF
    IF (.NOT. ok) RETURN
    dt = period/200.0_wp
    IF (what == AXIAL) THEN
      CALL CD_HermiteCable_Dyn_Set_Axial_Damping(model, [(2.0e4_wp, k=1, NE)], es, em)
      CALL require(es == CD_HCDYN_OK, label//': axial damping: '//TRIM(em))
    END IF
    IF (what == BENDING .OR. what == VIV) THEN
      CALL CD_HermiteCable_Dyn_Set_Bending_Damping(model, [(2.0e2_wp, k=1, NE)], es, em)
      CALL require(es == CD_HCDYN_OK, label//': bending damping: '//TRIM(em))
    END IF
    IF (what == DRAG .OR. what == VIV) THEN
      CALL CD_HermiteCable_Dyn_Set_Drag(model, RHOW, [(DIAM, k=1, NE)], [(1.2_wp, k=1, NE)], [(0.0_wp, k=1, NE)], &
                                        1.0e3_wp, [0.0_wp, MERGE(0.5_wp, 0.0_wp, what == VIV), 0.0_wp], es, em)
      CALL require(es == CD_HCDYN_OK, label//': drag: '//TRIM(em))
    END IF
    IF (what == VIV) THEN
      cl = 0.8_wp; cf = 0.18_wp; df = 0.08_wp; clx = 0.85_wp; alx = 0.15_wp; chh = 1.0_wp
      CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=clx, il_alpha=alx, chh=chh)
      CALL require(es == CD_HCDYN_OK, label//': VIV: '//TRIM(em))
    END IF
    IF (es /= CD_HCDYN_OK) RETURN
    CALL CD_HermiteCable_Dyn_Energy(model, ek, es0, es, em)
    e0 = ek + es0
    work = 0.0_wp
    pm_max = -HUGE(1.0_wp)
    zero_viv = .TRUE.
    zero_drag = .TRUE.
    CALL powers(model, label, p_prev, pm_max, zero_viv, zero_drag)
    DO n = 1, 2400
      CALL CD_HermiteCable_Dyn_Step(model, dt, 50, 1.0e-10_wp, es, em, it, res)
      IF (es /= CD_HCDYN_OK) THEN
        CALL require(.FALSE., label//': step failed: '//TRIM(em)); RETURN
      END IF
      CALL powers(model, label, p_now, pm_max, zero_viv, zero_drag)
      work = work + 0.5_wp*dt*(p_prev + p_now)
      p_prev = p_now
    END DO
    CALL CD_HermiteCable_Dyn_Energy(model, ek, es0, es, em)
    e1 = ek + es0
    scale = MAXVAL(ABS(work))
    WRITE (*, '(2A,ES12.4,A,3ES12.4,A,ES10.2)') label, ': dE ', e1 - e0, ', work VIV/drag/damp ', work, &
      ', mismatch/largest ', ABS(e1 - e0 - SUM(work))/MAX(scale, TINY(1.0_wp))
    CALL require(scale > 0.0_wp, label//': the loads did work')
    IF (.NOT. (scale > 0.0_wp)) RETURN
    CALL require(ABS(e1 - e0 - SUM(work)) <= 5.0e-3_wp*scale, label//': the load work closes the energy balance')
    CALL require(pm_max <= 0.0_wp, label//': the damping power is never positive at a node')
    CALL require(zero_viv .EQV. what /= VIV, label//': VIV power only with VIV')
    CALL require(zero_drag .EQV. (what /= DRAG .AND. what /= VIV), label//': drag power only with drag')
    IF (what == BENDING .OR. what == AXIAL) CALL require(work(3) < 0.0_wp, label//': the dashpot dissipates')
    IF (what == DRAG) CALL require(work(2) < 0.0_wp, label//': still-water drag dissipates')
    IF (what == VIV) CALL require(work(1) > 0.0_wp, label//': the VIV loads put energy in')
  END SUBROUTINE balance_case

  SUBROUTINE powers(model, label, p, pm_max, zero_viv, zero_drag)
    !! The line totals of the three load powers, with the running checks: the largest node
    !! damping power, and whether the VIV and drag powers stayed zero at every node.
    TYPE(CD_HermiteCableDynType), INTENT(IN) :: model
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), INTENT(OUT) :: p(3)
    REAL(wp), INTENT(INOUT) :: pm_max
    LOGICAL, INTENT(INOUT) :: zero_viv, zero_drag
    REAL(wp) :: pv(NN), pd(NN), pm(NN)
    INTEGER :: esp
    CHARACTER(512) :: emp
    CALL CD_HermiteCable_Dyn_Load_Power(model, pv, pd, pm, esp, emp)
    CALL require(esp == CD_HCDYN_OK, label//': load power: '//TRIM(emp))
    p = [SUM(pv), SUM(pd), SUM(pm)]
    pm_max = MAX(pm_max, MAXVAL(pm))
    zero_viv = zero_viv .AND. ALL(pv == 0.0_wp)
    zero_drag = zero_drag .AND. ALL(pd == 0.0_wp)
  END SUBROUTINE powers

  SUBROUTINE check_refusals()
    TYPE(CD_HermiteCableDynType) :: model
    REAL(wp) :: pv(NN), pd(NN), pm(NN - 1)
    INTEGER :: es
    LOGICAL :: ok
    CHARACTER(512) :: em
    CALL setup(model, 0.0_wp, 0.0_wp, ok)
    IF (.NOT. ok) RETURN
    CALL CD_HermiteCable_Dyn_Load_Power(model, pv, pd, pm, es, em)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'one entry per node') > 0, 'L5: a wrong-size power array is refused')
  END SUBROUTINE check_refusals

END PROGRAM test_hermite_load_power
