! File: tests/test_viv_kernels.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_kernels
  !! Gate for the node kernels of the cross-flow VIV synchronisation model (CableDyn_VIV):
  !!   1) the lift equals MoorDyn's Lf = 1/4 rho d |Vp| Cl cos(phi) (q x Vp) (l_i + l_i+1) through
  !!      coef = sum (1/2 l0)(1/2 rho D Cl); it is normal to the tangent and to the normal flow,
  !!      |F| = coef |cos(phi)| U^2, and vanishes for flow along the line;
  !!   2) its residual Jacobians (-dF/dr', -dF/dr_prev, -dF/dr_next) match central differences;
  !!   3) the cross-flow velocity and acceleration are the components along e = q x v_n / U;
  !!   4) the phase rate is 2 pi (cF + dF sin(phi_y - phi)) U/D, with phi_y from the RMS-normalised
  !!      y', y'' (cF at phase lock, cF + dF a quarter period ahead);
  !!   5) the running mean square is MoorDyn's ((n_m - 1) old^2 + new^2)/n_m;
  !!   6) the phase wraps to [0, 2 pi), and a node whose neighbours coincide fails closed;
  !!   7) the per-length lift of a centreline point (the finite-EI Gauss points);
  !!   8) in-line VIV (Ulveseter et al. 2018): the node and per-length forces are
  !!      1/2 rho D C_vx cos(phi_x) U v_n along the normal flow, their Jacobians match central
  !!      differences, the in-line motion is the component along v_n/U, and the phase rate is
  !!      2 rate_y [1 + alpha sin(phi_xd - phi_x)] (twice the cross-flow rate at phase lock).
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_VIV, ONLY: CD_VIV_Node_Lift, CD_VIV_Crossflow_Motion, CD_VIV_Phase_Rate, CD_VIV_Mean_Square, &
                          CD_VIV_Wrap_Phase, CD_VIV_TWO_PI, CD_VIV_NM, CD_VIV_Lift_Per_Length, &
                          CD_VIV_Node_Inline, CD_VIV_Inline_Per_Length, CD_VIV_Inline_Phase_Rate, &
                          CD_VIV_Crossflow_Along
  IMPLICIT NONE
  REAL(wp), PARAMETER :: RHO = 1025.0_wp, DIAM = 0.1_wp, CL = 0.8_wp, L1 = 2.0_wp, L2 = 1.5_wp
  REAL(wp) :: rp(3), rn(3), rd(3), u(3), coef, phi, f(3), q(3), vr(3), vn(3), un, w(3), ref(3)
  REAL(wp) :: jv(3, 3), jp(3, 3), jn(3, 3), fd(3, 3), fp(3), fm(3), h, e(3), yd, ydd, un2, rate
  REAL(wp) :: ms, ms_md, rpv(3), kk, fl(3), jr(3, 3), jt(3, 3), xd, xdd, yd2, ydd2, un3, rdd(3)
  LOGICAL :: ok
  INTEGER :: nfail, j, k

  nfail = 0
  rp = [0.1_wp, -0.2_wp, -50.3_wp]
  rn = [3.6_wp, 0.4_wp, -49.8_wp]
  rd = [0.05_wp, -0.12_wp, 0.31_wp]
  u = [0.2_wp, 0.15_wp, -0.03_wp]
  coef = (0.5_wp*L1 + 0.5_wp*L2)*0.5_wp*RHO*DIAM*CL
  phi = 0.7_wp

  ! (1) MoorDyn's lift and its geometry
  CALL CD_VIV_Node_Lift(rp, rn, rd, u, coef, phi, f, ok)
  CALL require(ok, 'lift: frame defined')
  q = (rn - rp)/NORM2(rn - rp)
  vr = u - rd
  vn = vr - DOT_PRODUCT(q, vr)*q
  un = NORM2(vn)
  w = cross(q, vn)
  ref = 0.25_wp*RHO*DIAM*un*CL*COS(phi)*w*(L1 + L2)
  CALL require(NORM2(f - ref) <= 1.0e-12_wp*NORM2(ref), 'lift equals MoorDyn''s Lf')
  CALL require(ABS(DOT_PRODUCT(f, q)) <= 1.0e-12_wp*NORM2(f) .AND. ABS(DOT_PRODUCT(f, vn)) <= 1.0e-12_wp*NORM2(f)*un, &
               'lift normal to the tangent and to the normal flow')
  CALL require(ABS(NORM2(f) - coef*ABS(COS(phi))*un*un) <= 1.0e-12_wp*NORM2(f), '|lift| = coef |cos phi| U^2')
  CALL CD_VIV_Node_Lift(rp, rn, rd, rd + 0.3_wp*q, coef, phi, f, ok)
  CALL require(ok .AND. NORM2(f) <= 1.0e-14_wp, 'no lift for flow along the line')

  ! (2) Jacobians against central differences
  CALL CD_VIV_Node_Lift(rp, rn, rd, u, coef, phi, f, ok, jv, jp, jn)
  h = 1.0e-6_wp
  DO j = 1, 3
    CALL CD_VIV_Node_Lift(rp, rn, rd + h*unit(j), u, coef, phi, fp, ok)
    CALL CD_VIV_Node_Lift(rp, rn, rd - h*unit(j), u, coef, phi, fm, ok)
    fd(:, j) = -(fp - fm)/(2.0_wp*h)
  END DO
  CALL require(MAXVAL(ABS(jv - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), 'velocity Jacobian matches differences')
  DO j = 1, 3
    CALL CD_VIV_Node_Lift(rp + h*unit(j), rn, rd, u, coef, phi, fp, ok)
    CALL CD_VIV_Node_Lift(rp - h*unit(j), rn, rd, u, coef, phi, fm, ok)
    fd(:, j) = -(fp - fm)/(2.0_wp*h)
  END DO
  CALL require(MAXVAL(ABS(jp - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), 'previous-node Jacobian matches differences')
  DO j = 1, 3
    CALL CD_VIV_Node_Lift(rp, rn + h*unit(j), rd, u, coef, phi, fp, ok)
    CALL CD_VIV_Node_Lift(rp, rn - h*unit(j), rd, u, coef, phi, fm, ok)
    fd(:, j) = -(fp - fm)/(2.0_wp*h)
  END DO
  CALL require(MAXVAL(ABS(jn - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), 'next-node Jacobian matches differences')

  ! (3) cross-flow motion along e = q x v_n / U
  CALL CD_VIV_Crossflow_Motion(rp, rn, rd, [1.0_wp, 2.0_wp, -0.5_wp], u, yd, ydd, un2, ok)
  e = w/un
  CALL require(ok .AND. ABS(un2 - un) <= 1.0e-14_wp, 'cross-flow: normal flow speed')
  CALL require(ABS(yd - DOT_PRODUCT(rd, e)) <= 1.0e-14_wp .AND. &
               ABS(ydd - DOT_PRODUCT([1.0_wp, 2.0_wp, -0.5_wp], e)) <= 1.0e-14_wp, &
               'cross-flow velocity and acceleration along e')

  ! (4) phase rate: at phase lock (phi_y = phi) the band centre, a quarter period ahead cF + dF
  ! (y' = cos, y'' = -sin with unit RMS: phi_y is the motion phase)
  rate = CD_VIV_Phase_Rate(COS(1.1_wp), -SIN(1.1_wp), 1.0_wp, 1.0_wp, 1.1_wp, 0.15_wp, DIAM, 0.18_wp, &
                           0.08_wp)
  CALL require(ABS(rate - CD_VIV_TWO_PI*0.18_wp*0.15_wp/DIAM) <= 1.0e-12_wp*rate, 'rate at phase lock is cF U/D')
  rate = CD_VIV_Phase_Rate(COS(1.1_wp), -SIN(1.1_wp), 1.0_wp, 1.0_wp, 1.1_wp - 0.25_wp*CD_VIV_TWO_PI, &
                           0.15_wp, DIAM, 0.18_wp, 0.08_wp)
  CALL require(ABS(rate - CD_VIV_TWO_PI*0.26_wp*0.15_wp/DIAM) <= 1.0e-12_wp*rate, 'rate a quarter period ahead is (cF+dF) U/D')

  ! (5) MoorDyn's running mean square
  ms_md = 0.3_wp
  DO k = 1, 50
    ms = CD_VIV_Mean_Square(ms_md, SIN(0.1_wp*k), 1.0_wp/REAL(CD_VIV_NM, wp))
    CALL require(ABS(ms - (REAL(CD_VIV_NM - 1, wp)*ms_md + SIN(0.1_wp*k)**2)/REAL(CD_VIV_NM, wp)) <= 1.0e-15_wp, &
                 'mean square is MoorDyn''s running RMS update')
    ms_md = ms
  END DO

  ! (7) per-length lift at a centreline point (the finite-EI Gauss points): the node formula with
  ! the exact tangent, scaled to unit length, independent of |dr/dxi|; Jacobians vs differences
  rpv = 2.5_wp*(rn - rp)
  kk = 0.5_wp*RHO*DIAM*CL*COS(phi)
  CALL CD_VIV_Lift_Per_Length(u - rd, rpv, kk, fl, ok, jr, jt)
  CALL CD_VIV_Node_Lift(rp, rn, rd, u, 1.0_wp, phi, f, ok)
  CALL require(ok .AND. NORM2(fl - 0.5_wp*RHO*DIAM*CL*f) <= 1.0e-12_wp*NORM2(fl), &
               'per-length lift = 1/2 rho D Cl cos(phi) U (t x v_n)')
  DO j = 1, 3
    CALL CD_VIV_Lift_Per_Length(u - rd + h*unit(j), rpv, kk, fp, ok)
    CALL CD_VIV_Lift_Per_Length(u - rd - h*unit(j), rpv, kk, fm, ok)
    fd(:, j) = (fp - fm)/(2.0_wp*h)
  END DO
  CALL require(MAXVAL(ABS(jr - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), 'per-length lift: flow Jacobian matches differences')
  DO j = 1, 3
    CALL CD_VIV_Lift_Per_Length(u - rd, rpv + h*unit(j), kk, fp, ok)
    CALL CD_VIV_Lift_Per_Length(u - rd, rpv - h*unit(j), kk, fm, ok)
    fd(:, j) = (fp - fm)/(2.0_wp*h)
  END DO
  CALL require(MAXVAL(ABS(jt - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), &
               'per-length lift: tangent Jacobian matches differences')
  CALL CD_VIV_Lift_Per_Length(u - rd, 3.0_wp*rpv, kk, fp, ok)
  CALL require(NORM2(fp - fl) <= 1.0e-12_wp*NORM2(fl), 'per-length lift does not depend on |dr/dxi|')

  ! (8) in-line VIV: the node force along the normal flow and its Jacobians
  CALL CD_VIV_Node_Inline(rp, rn, rd, u, coef, phi, f, ok)
  ref = 0.25_wp*RHO*DIAM*un*CL*COS(phi)*vn*(L1 + L2)
  CALL require(ok .AND. NORM2(f - ref) <= 1.0e-12_wp*NORM2(ref), 'in-line: node force = coef cos(phi_x) U v_n')
  CALL require(ABS(DOT_PRODUCT(f, q)) <= 1.0e-12_wp*NORM2(f), 'in-line: node force normal to the tangent')
  CALL CD_VIV_Node_Inline(rp, rn, rd, rd + 0.3_wp*q, coef, phi, f, ok)
  CALL require(ok .AND. NORM2(f) <= 1.0e-14_wp, 'in-line: no force for flow along the line')
  CALL CD_VIV_Node_Inline(rp, rn, rd, u, coef, phi, f, ok, jv, jp, jn)
  DO j = 1, 3
    CALL CD_VIV_Node_Inline(rp, rn, rd + h*unit(j), u, coef, phi, fp, ok)
    CALL CD_VIV_Node_Inline(rp, rn, rd - h*unit(j), u, coef, phi, fm, ok)
    fd(:, j) = -(fp - fm)/(2.0_wp*h)
  END DO
  CALL require(MAXVAL(ABS(jv - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), 'in-line: velocity Jacobian matches differences')
  DO j = 1, 3
    CALL CD_VIV_Node_Inline(rp + h*unit(j), rn, rd, u, coef, phi, fp, ok)
    CALL CD_VIV_Node_Inline(rp - h*unit(j), rn, rd, u, coef, phi, fm, ok)
    fd(:, j) = -(fp - fm)/(2.0_wp*h)
  END DO
  CALL require(MAXVAL(ABS(jp - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), 'in-line: previous-node Jacobian matches differences')
  DO j = 1, 3
    CALL CD_VIV_Node_Inline(rp, rn + h*unit(j), rd, u, coef, phi, fp, ok)
    CALL CD_VIV_Node_Inline(rp, rn - h*unit(j), rd, u, coef, phi, fm, ok)
    fd(:, j) = -(fp - fm)/(2.0_wp*h)
  END DO
  CALL require(MAXVAL(ABS(jn - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), 'in-line: next-node Jacobian matches differences')
  CALL CD_VIV_Node_Inline(rp, rp, rd, u, coef, phi, f, ok)
  CALL require(.NOT. ok .AND. NORM2(f) <= 0.0_wp, 'in-line: coinciding neighbours fail closed with no force')

  ! in-line force per unit length: the node formula with the exact tangent, independent of |dr/dxi|
  CALL CD_VIV_Inline_Per_Length(u - rd, rpv, kk, fl, ok, jr, jt)
  CALL CD_VIV_Node_Inline(rp, rn, rd, u, 1.0_wp, phi, f, ok)
  CALL require(ok .AND. NORM2(fl - 0.5_wp*RHO*DIAM*CL*f) <= 1.0e-12_wp*NORM2(fl), &
               'in-line: per-length force = 1/2 rho D C_vx cos(phi_x) U v_n')
  DO j = 1, 3
    CALL CD_VIV_Inline_Per_Length(u - rd + h*unit(j), rpv, kk, fp, ok)
    CALL CD_VIV_Inline_Per_Length(u - rd - h*unit(j), rpv, kk, fm, ok)
    fd(:, j) = (fp - fm)/(2.0_wp*h)
  END DO
CALL require(MAXVAL(ABS(jr - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), 'in-line: per-length flow Jacobian matches differences')
  DO j = 1, 3
    CALL CD_VIV_Inline_Per_Length(u - rd, rpv + h*unit(j), kk, fp, ok)
    CALL CD_VIV_Inline_Per_Length(u - rd, rpv - h*unit(j), kk, fm, ok)
    fd(:, j) = (fp - fm)/(2.0_wp*h)
  END DO
  CALL require(MAXVAL(ABS(jt - fd)) <= 1.0e-6_wp*MAXVAL(ABS(fd)), &
               'in-line: per-length tangent Jacobian matches differences')
  CALL CD_VIV_Inline_Per_Length(u - rd, 3.0_wp*rpv, kk, fp, ok)
  CALL require(NORM2(fp - fl) <= 1.0e-12_wp*NORM2(fl), 'in-line: per-length force does not depend on |dr/dxi|')

  ! in-line motion along v_n/U, from both motion routines; the cross-flow outputs are unchanged
  rdd = [1.0_wp, 2.0_wp, -0.5_wp]
  CALL CD_VIV_Crossflow_Motion(rp, rn, rd, rdd, u, yd2, ydd2, un3, ok, xd, xdd)
  CALL require(ok .AND. yd2 == yd .AND. ydd2 == ydd, 'in-line: the cross-flow motion is unchanged')
  CALL require(ABS(xd - DOT_PRODUCT(rd, vn)/un) <= 1.0e-14_wp .AND. ABS(xdd - DOT_PRODUCT(rdd, vn)/un) <= 1.0e-14_wp, &
               'in-line: velocity and acceleration along v_n/U')
  CALL CD_VIV_Crossflow_Along(rn - rp, rd, rdd, u, yd2, ydd2, un3, ok, xd, xdd)
  CALL require(ok .AND. ABS(yd2 - yd) <= 1.0e-14_wp .AND. ABS(xd - DOT_PRODUCT(rd, vn)/un) <= 1.0e-14_wp .AND. &
               ABS(xdd - DOT_PRODUCT(rdd, vn)/un) <= 1.0e-14_wp, 'in-line: motion along a given tangent')

  ! in-line phase rate: twice the cross-flow rate at phase lock, (1 + alpha) times that a quarter
  ! period ahead; a node with no in-line motion takes phi_xd = 0
  rate = CD_VIV_Inline_Phase_Rate(COS(2.2_wp), -SIN(2.2_wp), 1.0_wp, 1.0_wp, 2.2_wp, 3.0_wp, 0.15_wp)
  CALL require(ABS(rate - 6.0_wp) <= 1.0e-12_wp*rate, 'in-line: rate at phase lock is twice the cross-flow rate')
  rate = CD_VIV_Inline_Phase_Rate(COS(2.2_wp), -SIN(2.2_wp), 1.0_wp, 1.0_wp, 2.2_wp - 0.25_wp*CD_VIV_TWO_PI, &
                                  3.0_wp, 0.15_wp)
  CALL require(ABS(rate - 6.0_wp*1.15_wp) <= 1.0e-12_wp*rate, 'in-line: rate a quarter period ahead is 2 (1+alpha) rate_y')
  rate = CD_VIV_Inline_Phase_Rate(0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp, 3.0_wp, 0.15_wp)
  CALL require(ABS(rate - 6.0_wp*(1.0_wp + 0.15_wp*SIN(-0.5_wp))) <= 1.0e-12_wp*rate, &
               'in-line: no in-line motion takes phi_xd = 0')

  ! (6) wrap and degenerate frame
  CALL require(CD_VIV_Wrap_Phase(-1.0e-300_wp) >= 0.0_wp .AND. CD_VIV_Wrap_Phase(-1.0e-300_wp) < CD_VIV_TWO_PI, &
               'a tiny negative phase wraps into [0, 2 pi)')
  CALL require(ABS(CD_VIV_Wrap_Phase(7.0_wp) - (7.0_wp - CD_VIV_TWO_PI)) <= 1.0e-15_wp, 'phase wraps by 2 pi')
  CALL CD_VIV_Node_Lift(rp, rp, rd, u, coef, phi, f, ok)
  CALL require(.NOT. ok .AND. NORM2(f) <= 0.0_wp, 'coinciding neighbours fail closed with no lift')

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: VIV synchronisation-model kernels'

CONTAINS

  PURE FUNCTION cross(a, b) RESULT(c)
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  END FUNCTION cross

  PURE FUNCTION unit(j) RESULT(e)
    INTEGER, INTENT(IN) :: j
    REAL(wp) :: e(3)
    e = 0.0_wp
    e(j) = 1.0_wp
  END FUNCTION unit

  SUBROUTINE require(cond, msg)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: msg
    IF (.NOT. cond) THEN
      WRITE (*, '(A)') 'MISMATCH: '//msg
      nfail = nfail + 1
    END IF
  END SUBROUTINE require

END PROGRAM test_viv_kernels
