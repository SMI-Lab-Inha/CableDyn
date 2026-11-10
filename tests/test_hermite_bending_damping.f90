! File: tests/test_hermite_bending_damping.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_hermite_bending_damping
  !! Kernel gates of the Kelvin-Voigt bending damping (CD_HermiteCable_Bending_Damping):
  !!   K1 Kv equals the central difference of fdamp in the rates, Kq that in the positions, on a
  !!      stretched, twisted 3D element and on a straight one;
  !!   K2 fdamp is linear in the rates, dissipative (v . fdamp >= 0) and zero for Bb = 0;
  !!   K3 zero for a rigid translation and for a rigid rotation whose axis is normal to the
  !!      tangent everywhere (a planar arc spun about the normal of its plane);
  !!   K4 the planar rate: a straight element bent in its plane at curvature rate kdot dissipates
  !!      Bb kdot^2 L, the dissipation function's closed form;
  !!   K5 bad inputs are refused by name.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_HermiteCable, ONLY: CD_HermiteCable_Bending_Damping, CD_HCABLE_OK
  IMPLICIT NONE
  INTEGER :: nfail = 0
  REAL(wp), PARAMETER :: L = 2.5_wp, BB = 37.0_wp

  CALL check_jacobians()
  CALL check_linear_dissipative()
  CALL check_rigid_motions()
  CALL check_planar_rate()
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: Hermite bending damping kernel'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE curved_state(q, v)
    !! A stretched, out-of-plane curved element with arbitrary nodal rates.
    REAL(wp), INTENT(OUT) :: q(12), v(12)
    q = [0.0_wp, 0.0_wp, 0.0_wp, 1.02_wp, 0.3_wp, -0.2_wp, &
         2.4_wp, 0.5_wp, -0.3_wp, 0.9_wp, -0.25_wp, 0.35_wp]
    v = [0.11_wp, -0.23_wp, 0.31_wp, 0.4_wp, -0.15_wp, 0.27_wp, &
         -0.19_wp, 0.08_wp, 0.42_wp, -0.33_wp, 0.21_wp, 0.12_wp]
  END SUBROUTINE curved_state

  SUBROUTINE fd_check(q, v, label)
    REAL(wp), INTENT(IN) :: q(12), v(12)
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp) :: f(12), Kq(12, 12), Kv(12, 12), fp(12), fm(12), Kq_fd(12, 12), Kv_fd(12, 12), dq(12), h
    REAL(wp) :: eq, ev
    INTEGER :: es, j
    CHARACTER(256) :: em
    CALL CD_HermiteCable_Bending_Damping(q, v, L, BB, f, Kq, Kv, es, em)
    CALL require(es == CD_HCABLE_OK, 'K1 '//label//': evaluation: '//TRIM(em))
    DO j = 1, 12
      h = 1.0e-6_wp*MAX(1.0_wp, ABS(q(j)))
      dq = 0.0_wp; dq(j) = h
      CALL CD_HermiteCable_Bending_Damping(q + dq, v, L, BB, fp, ErrStat=es, ErrMsg=em)
      CALL CD_HermiteCable_Bending_Damping(q - dq, v, L, BB, fm, ErrStat=es, ErrMsg=em)
      Kq_fd(:, j) = (fp - fm)/(2.0_wp*h)
      h = 1.0e-6_wp*MAX(1.0_wp, ABS(v(j)))
      dq = 0.0_wp; dq(j) = h
      CALL CD_HermiteCable_Bending_Damping(q, v + dq, L, BB, fp, ErrStat=es, ErrMsg=em)
      CALL CD_HermiteCable_Bending_Damping(q, v - dq, L, BB, fm, ErrStat=es, ErrMsg=em)
      Kv_fd(:, j) = (fp - fm)/(2.0_wp*h)
    END DO
    eq = MAXVAL(ABS(Kq - Kq_fd))/MAX(1.0e-12_wp, MAXVAL(ABS(Kq_fd)))
    ev = MAXVAL(ABS(Kv - Kv_fd))/MAX(1.0e-12_wp, MAXVAL(ABS(Kv_fd)))
    WRITE (*, '(A,2ES11.3)') 'K1 '//label//': relative |Kq - FD|, |Kv - FD| ', eq, ev
    CALL require(eq < 1.0e-6_wp, 'K1 '//label//': Kq matches the central difference')
    CALL require(ev < 1.0e-6_wp, 'K1 '//label//': Kv matches the central difference')
    CALL require(MAXVAL(ABS(Kv - TRANSPOSE(Kv))) <= 1.0e-12_wp*MAXVAL(ABS(Kv)), 'K1 '//label//': Kv symmetric')
  END SUBROUTINE fd_check

  SUBROUTINE check_jacobians()
    REAL(wp) :: q(12), v(12)
    CALL curved_state(q, v)
    CALL fd_check(q, v, 'curved 3D')
    ! straight, unstretched element along x (the zero-curvature state)
    q = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, L, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
    CALL fd_check(q, v, 'straight')
  END SUBROUTINE check_jacobians

  SUBROUTINE check_linear_dissipative()
    REAL(wp) :: q(12), v(12), f1(12), f2(12)
    INTEGER :: es
    CHARACTER(256) :: em
    CALL curved_state(q, v)
    CALL CD_HermiteCable_Bending_Damping(q, v, L, BB, f1, ErrStat=es, ErrMsg=em)
    CALL CD_HermiteCable_Bending_Damping(q, 2.0_wp*v, L, BB, f2, ErrStat=es, ErrMsg=em)
    CALL require(MAXVAL(ABS(f2 - 2.0_wp*f1)) <= 1.0e-12_wp*MAXVAL(ABS(f1)), 'K2: linear in the rates')
    CALL require(DOT_PRODUCT(v, f1) > 0.0_wp, 'K2: dissipative (v . f > 0)')
    CALL CD_HermiteCable_Bending_Damping(q, v, L, 0.0_wp, f1, ErrStat=es, ErrMsg=em)
    CALL require(es == CD_HCABLE_OK .AND. MAXVAL(ABS(f1)) == 0.0_wp, 'K2: zero for Bb = 0')
  END SUBROUTINE check_linear_dissipative

  SUBROUTINE check_rigid_motions()
    REAL(wp) :: q(12), v(12), f(12), scale, omega(3)
    INTEGER :: es, k
    CHARACTER(256) :: em
    CALL curved_state(q, v)
    ! rigid translation: every node moves with c, tangents fixed
    v = [0.3_wp, -0.7_wp, 0.2_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.3_wp, -0.7_wp, 0.2_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    CALL CD_HermiteCable_Bending_Damping(q, v, L, BB, f, ErrStat=es, ErrMsg=em)
    CALL require(MAXVAL(ABS(f)) <= 1.0e-12_wp*BB, 'K3: zero for a rigid translation')
    ! a planar arc in the x-z plane spun about y (normal to every tangent): v = omega x r, m' = omega x m
    q = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.4_wp, 2.3_wp, 0.0_wp, 0.6_wp, 0.8_wp, 0.0_wp, -0.3_wp]
    omega = [0.0_wp, 0.9_wp, 0.0_wp]
    DO k = 0, 3
      v(3*k + 1:3*k + 3) = cross(omega, q(3*k + 1:3*k + 3))
    END DO
    CALL CD_HermiteCable_Bending_Damping(q, v, L, BB, f, ErrStat=es, ErrMsg=em)
    ! reference: the same element bent at the same rate scale
    scale = BB*0.9_wp/L
    WRITE (*, '(A,ES11.3)') 'K3: rigid in-plane rotation, max |f| / (Bb omega / L) ', MAXVAL(ABS(f))/scale
    CALL require(MAXVAL(ABS(f)) <= 1.0e-12_wp*scale, 'K3: zero for a rigid rotation normal to the tangents')
  END SUBROUTINE check_rigid_motions

  SUBROUTINE check_planar_rate()
    !! A straight unstretched element along x bent in the x-z plane with rates that give a
    !! uniform curvature rate kdot (w(x) = 1/2 kdot x^2 t, end slopes kdot x): the dissipation
    !! v . f = 2 psi = Bb kdot^2 L at the straight state.
    REAL(wp) :: q(12), v(12), f(12), kdot, d
    INTEGER :: es
    CHARACTER(256) :: em
    kdot = 0.17_wp
    q = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, L, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
    v = 0.0_wp
    ! w = kdot x^2 / 2: node 2 rate kdot L^2/2, slope handle m = dr/ds -> rate kdot L
    v(9) = 0.5_wp*kdot*L*L
    v(12) = kdot*L
    CALL CD_HermiteCable_Bending_Damping(q, v, L, BB, f, ErrStat=es, ErrMsg=em)
    d = DOT_PRODUCT(v, f)
    WRITE (*, '(A,ES14.6,A,ES14.6)') 'K4: v . f ', d, ', Bb kdot^2 L ', BB*kdot*kdot*L
    CALL require(ABS(d - BB*kdot*kdot*L) <= 1.0e-10_wp*BB*kdot*kdot*L, 'K4: planar dissipation Bb kdot^2 L')
  END SUBROUTINE check_planar_rate

  SUBROUTINE check_refusals()
    REAL(wp) :: q(12), v(12), f(12)
    INTEGER :: es
    CHARACTER(256) :: em
    CALL curved_state(q, v)
    CALL CD_HermiteCable_Bending_Damping(q, v, L, -1.0_wp, f, ErrStat=es, ErrMsg=em)
    CALL require(es /= CD_HCABLE_OK .AND. INDEX(em, 'Bb must be non-negative') > 0, 'K5: negative Bb')
    CALL CD_HermiteCable_Bending_Damping(q, v, 0.0_wp, BB, f, ErrStat=es, ErrMsg=em)
    CALL require(es /= CD_HCABLE_OK .AND. INDEX(em, 'length must be positive') > 0, 'K5: zero length')
    q(4:6) = 0.0_wp; q(7:9) = 0.0_wp; q(10:12) = 0.0_wp
    CALL CD_HermiteCable_Bending_Damping(q, v, L, BB, f, ErrStat=es, ErrMsg=em)
    CALL require(es /= CD_HCABLE_OK .AND. INDEX(em, 'degenerate centreline') > 0, 'K5: collapsed element')
  END SUBROUTINE check_refusals

  PURE FUNCTION cross(a, b) RESULT(c)
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  END FUNCTION cross

END PROGRAM test_hermite_bending_damping
