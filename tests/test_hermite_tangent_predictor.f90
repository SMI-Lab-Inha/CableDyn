! File: tests/test_hermite_tangent_predictor.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_hermite_tangent_predictor
  !! The Newton predictor of a finite-EI step and the stretch mode of the nodal tangents
  !! (CableDyn_HermiteCableDynamic, dyn_step_core). A nodal tangent m = dr/ds carries the line's
  !! stiffest motion, its stretch (EA against the tangent's small consistent mass, here omega ~
  !! 8e3 rad/s, omega dt ~ 300 and 600); with little numerical damping it keeps ringing at the offset
  !! it starts from, and the committed velocity and acceleration of that mode are far larger than
  !! the motion it makes. The constant-acceleration and Newmark guesses extrapolated them and carried
  !! the tangents across zero, and Newton converged to the root with every tangent reversed (each
  !! element folded into a Z along its chord, 2.4e8 J of strain energy against 3.5 J), from which the
  !! march failed a few steps on. A guess that changes the length of a free tangent by more than a
  !! quarter of it now holds that tangent at its committed value, and a step that still ends with a
  !! tangent turned by 90 degrees or more fails.
  !!
  !! The line: a weightless taut line (20 m, 10 elements, D 0.1 m, 8.1 kg/m, EA 5e7 N, EI 1e3 N m^2,
  !! 3.8 kN, no added mass), pinned at both ends in a 0.3 m/s cross current, Cd 1.2, and the VIV lift
  !! (Cl 0.8, cF 0.18, dF 0.08, running RMS window 0.25 s), under which it locks in at about 0.63 Hz
  !! with a mid-span A/D of 0.56-0.57 (the lock-in case of the long-run VIV checks). Its seed is the
  !! straight line with unit tangents, while the elements are stretched by T/EA = 7.6e-5: the tangent
  !! stretch mode starts 7.6e-5 off its balance (|a| = 4.9e3 1/s^2 on the tangents). The reference
  !! is the same line seeded with tangents of length 1 + T/EA (balanced), which every predictor ran
  !! before and runs unchanged. 20 cycles of the nominal shedding period D/(cF U) = 1.85 s (37 s).
  !!
  !! Gates, for rho_inf 0.5, 0.8, 0.9 and 1 at 25 and 50 steps per cycle:
  !! * G1 every step of the unit-tangent seed converges (before: at rho_inf 0.8, 0.9 and 1 the run
  !!   failed at step 5, 13, 4 (25 steps per cycle) and 10, 7, 6 (50), each after steps on the
  !!   reversed root; at 25 steps per cycle and rho_inf 0.5 it ran through reversed roots to a
  !!   mid-span A/D of 5.7);
  !! * G2 its largest strain energy stays within 1.5 times the reference's (measured at most 1.20,
  !!   at 25 steps per cycle and rho_inf 1; the reversed root has 3e7 times it);
  !! * G3 its largest mid-span cross-flow displacement is the reference's within 2 % (measured at
  !!   most 0.75 %), so the ringing tangent stretch does not change the response;
  !! * G4 it takes at most 8 Newton iterations per step on average and no step more than 30 (the deck
  !!   default budget; measured at most 5.9 and 23, at 25 steps per cycle and rho_inf 1; before,
  !!   the first step alone took 21 at 50 steps per cycle and rho_inf 1, now 3);
  !! * G5 the same line without the lift at 50 steps per cycle and rho_inf 1 (before: failed at step
  !!   38 after steps on the reversed root): every step converges, the largest strain energy is
  !!   within 1.5 times the reference's (measured 1.09) and the mid-span in-line deflection at the end
  !!   within 1 % (measured to 4 digits);
  !! * G6 the same no-lift line kept in its plane (r_z and m_z fixed at every node, the planar
  !!   constraint pattern): the guard applies to the free tangent components of a partly
  !!   constrained node, so it runs as G5 (before: the guard skipped every node with a fixed
  !!   tangent component, and the reversed-root check stopped the run).
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Set_Drag, CD_HermiteCable_Dyn_Set_VIV, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_Energy, &
                                          CD_HermiteCable_Dyn_End, CD_HCDYN_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: RHOW = 1025.0_wp, D = 0.1_wp, CF = 0.18_wp, DF = 0.08_wp, U = 0.3_wp
  REAL(wp), PARAMETER :: EA = 5.0e7_wp, TENS = 3.8e3_wp
  INTEGER, PARAMETER :: NCYC = 20

  !> One run: whether every step converged, the largest strain energy, the largest mid-span
  !> cross-flow displacement, the mid-span in-line deflection at the end, the mean iterations.
  TYPE :: run_result
    LOGICAL :: ok = .FALSE.
    INTEGER :: nstep = 0
    INTEGER :: max_iters = 0
    REAL(wp) :: es_max = 0.0_wp, z_max = 0.0_wp, y_end = 0.0_wp, iters = 0.0_wp
  END TYPE run_result

  INTEGER :: nfail = 0, ir, ip
  REAL(wp), PARAMETER :: RHOS(4) = [0.5_wp, 0.8_wp, 0.9_wp, 1.0_wp]
  INTEGER, PARAMETER :: NPERS(2) = [25, 50]
  TYPE(run_result) :: unit, ref
  CHARACTER(64) :: lab

  DO ip = 1, SIZE(NPERS)
    DO ir = 1, SIZE(RHOS)
      WRITE (lab, '(A,F4.2,A,I0,A)') 'rho_inf ', RHOS(ir), ', ', NPERS(ip), ' steps per cycle'
      CALL run(NPERS(ip), RHOS(ir), .TRUE., .FALSE., 1.0_wp, unit)
      CALL run(NPERS(ip), RHOS(ir), .TRUE., .FALSE., 1.0_wp + TENS/EA, ref)
      CALL report(lab, unit, ref)
      CALL require(ref%ok, 'reference '//TRIM(lab)//': every step converges')
      CALL require(unit%ok, 'G1 '//TRIM(lab)//': every step of the unit-tangent seed converges')
      IF (.NOT. (unit%ok .AND. ref%ok)) CYCLE
      CALL require(unit%es_max <= 1.5_wp*ref%es_max, 'G2 '//TRIM(lab)//': the strain energy stays bounded')
      CALL require(ABS(unit%z_max - ref%z_max) <= 0.02_wp*ref%z_max, &
                   'G3 '//TRIM(lab)//': the cross-flow response is the balanced seed''s')
      CALL require(unit%iters <= 8.0_wp .AND. unit%max_iters <= 30, &
                   'G4 '//TRIM(lab)//': Newton converges in a few iterations')
    END DO
  END DO
  CALL run(50, 1.0_wp, .FALSE., .FALSE., 1.0_wp, unit)
  CALL run(50, 1.0_wp, .FALSE., .FALSE., 1.0_wp + TENS/EA, ref)
  CALL report('no lift, rho_inf 1.00, 50 steps per cycle', unit, ref)
  CALL no_lift_gates('G5 no lift', unit, ref)
  CALL run(50, 1.0_wp, .FALSE., .TRUE., 1.0_wp, unit)
  CALL run(50, 1.0_wp, .FALSE., .TRUE., 1.0_wp + TENS/EA, ref)
  CALL report('planar, no lift, rho_inf 1.00, 50 steps per cycle', unit, ref)
  CALL no_lift_gates('G6 planar no lift', unit, ref)
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: the predictor holds the ringing tangent stretch'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE no_lift_gates(label, a, b)
    CHARACTER(*), INTENT(IN) :: label
    TYPE(run_result), INTENT(IN) :: a, b
    CALL require(b%ok .AND. a%ok, label//': every step converges')
    IF (.NOT. (b%ok .AND. a%ok)) RETURN
    CALL require(a%es_max <= 1.5_wp*b%es_max, label//': the strain energy stays bounded')
    CALL require(ABS(a%y_end - b%y_end) <= 0.01_wp*ABS(b%y_end), &
                 label//': the in-line deflection is the balanced seed''s')
  END SUBROUTINE no_lift_gates

  SUBROUTINE report(label, a, b)
    CHARACTER(*), INTENT(IN) :: label
    TYPE(run_result), INTENT(IN) :: a, b
    WRITE (*, '(2X,A,A,I0,A,L1,A,3ES11.3,F6.2,I4,A,3ES11.3,F6.2,I4)') TRIM(label), ': unit seed steps ', a%nstep, &
      ' ok ', a%ok, ' (Es max, z max/D, y end, mean and max iterations) ', a%es_max, a%z_max/D, a%y_end, a%iters, &
      a%max_iters, &
      '; balanced ', b%es_max, b%z_max/D, b%y_end, b%iters, b%max_iters
  END SUBROUTINE report

  SUBROUTINE run(nper, rho_inf, with_lift, planar, m0, res)
    !! The taut line from its straight seed with tangents of length m0, over NCYC nominal cycles;
    !! planar fixes r_z and m_z at every node (the line moves in the plane of the current).
    INTEGER, INTENT(IN) :: nper
    REAL(wp), INTENT(IN) :: rho_inf, m0
    LOGICAL, INTENT(IN) :: with_lift, planar
    TYPE(run_result), INTENT(OUT) :: res
    INTEGER, PARAMETER :: NE = 10, NN = NE + 1, MIDN = NE/2 + 1
    REAL(wp), PARAMETER :: SPAN = 20.0_wp, Z0 = -50.0_wp, MLIN = 8.1_wp
    TYPE(CD_HermiteCableDynType) :: hm
    REAL(wp) :: seed(6*NN), l0(NE), eav(NE), eiv(NE), rhoa(NE), wv(NE), hd(NE), hn(NE), ht(NE)
    REAL(wp) :: cl(NE), cfv(NE), dfv(NE), dt, ek, es_
    INTEGER :: i, n, es, its, nit
    INTEGER, ALLOCATABLE :: fixed(:)
    CHARACTER(1024) :: em
    dt = D/(CF*U)/REAL(nper, wp)
    l0 = SPAN/REAL(NE, wp)/(1.0_wp + TENS/EA)
    eav = EA; rhoa = MLIN; eiv = 1.0e3_wp; wv = 0.0_wp
    seed = 0.0_wp
    DO i = 1, NN
      seed(6*(i - 1) + 1) = SPAN*REAL(i - 1, wp)/REAL(NE, wp)
      seed(6*(i - 1) + 3) = Z0
      seed(6*(i - 1) + 4) = m0
    END DO
    fixed = [1, 2, 3, 6*NE + 1, 6*NE + 2, 6*NE + 3]
    IF (planar) THEN
      fixed = [1, 2, 3, 6]
      DO i = 2, NN - 1
        fixed = [fixed, 6*(i - 1) + 3, 6*(i - 1) + 6]
      END DO
      fixed = [fixed, 6*NE + 1, 6*NE + 2, 6*NE + 3, 6*NE + 6]
    END IF
    CALL CD_HermiteCable_Dyn_Init(hm, l0, eav, eiv, rhoa, wv, seed, fixed, -1.0e4_wp, 0.0_wp, rho_inf, es, em)
    IF (es == CD_HCDYN_OK) THEN
      hd = D; hn = 1.2_wp; ht = 0.0_wp
      CALL CD_HermiteCable_Dyn_Set_Drag(hm, RHOW, hd, hn, ht, 0.0_wp, [0.0_wp, U, 0.0_wp], es, em)
    END IF
    IF (es == CD_HCDYN_OK .AND. with_lift) THEN
      cl = 0.8_wp; cfv = CF; dfv = DF
      CALL CD_HermiteCable_Dyn_Set_VIV(hm, cl, cfv, dfv, 0.25_wp, es, em)
    END IF
    CALL require(es == CD_HCDYN_OK, 'the line sets up: '//TRIM(em))
    IF (es /= CD_HCDYN_OK) RETURN
    nit = 0
    res%ok = .TRUE.
    DO n = 1, nper*NCYC
      CALL CD_HermiteCable_Dyn_Step(hm, dt, 50, 1.0e-8_wp, es, em, iters_out=its)
      IF (es /= CD_HCDYN_OK) THEN
        res%ok = .FALSE.
        WRITE (*, '(A,I0,A)') '    step ', n, ': '//TRIM(em)
        EXIT
      END IF
      res%nstep = n
      nit = nit + its
      res%max_iters = MAX(res%max_iters, its)
      CALL CD_HermiteCable_Dyn_Energy(hm, ek, es_, es, em)
      IF (es /= CD_HCDYN_OK) THEN
        res%ok = .FALSE.
        EXIT
      END IF
      res%es_max = MAX(res%es_max, es_)
      res%z_max = MAX(res%z_max, ABS(hm%q(6*(MIDN - 1) + 3) - Z0))
    END DO
    res%y_end = hm%q(6*(MIDN - 1) + 2)
    res%iters = REAL(nit, wp)/REAL(MAX(1, res%nstep), wp)
    CALL CD_HermiteCable_Dyn_End(hm)
  END SUBROUTINE run

END PROGRAM test_hermite_tangent_predictor
