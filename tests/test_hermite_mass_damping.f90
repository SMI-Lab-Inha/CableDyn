! File: tests/test_hermite_mass_damping.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_hermite_mass_damping
  !! Kernel gates of the mass-proportional (Rayleigh alpha) transverse damping
  !! (CD_HermiteCable_Mass_Damping) and of the two-frequency Rayleigh pair:
  !!   M1 Kv equals the central difference of fdamp in the rates, Kq that in the positions, on a
  !!      stretched 3D curved element and on a straight one; Kv is symmetric;
  !!   M2 fdamp is d psi / d v of the dissipation function 1/2 cm |v|^2 - 1/2 cm (t . v)^2 built
  !!      independently from the public shape functions, v . fdamp = 2 psi >= 0, and fdamp is
  !!      linear in the rates and zero for cm = 0;
  !!   M3 on a straight element along x the force is cm M v in y and z with M the consistent mass
  !!      matrix (CD_HermiteCable_Mass, per unit mass), zero in x, and an axial sliding is not
  !!      damped; a rigid transverse translation is (the absolute-velocity convention);
  !!   M4 bad inputs are refused by name;
  !!   M5 CD_HermiteCable_Rayleigh_Pair: alpha/(2 omega) + (Bb/EI) omega/2 equals zeta at f1 and
  !!      f2, and its minimum 2 zeta sqrt(f1 f2)/(f1 + f2) lies at sqrt(f1 f2).
  USE CableDyn_Precision, ONLY: wp
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
  USE CableDyn_HermiteCable, ONLY: CD_HermiteCable_Mass_Damping, CD_HermiteCable_Rayleigh_Pair, CD_HCABLE_OK, &
                                   CD_HermiteCable_Mass, CD_HermiteCable_Shapes, CD_HermiteCable_Gauss_Rule
  IMPLICIT NONE
  INTEGER :: nfail = 0
  REAL(wp), PARAMETER :: L = 2.5_wp, CM = 13.0_wp, PI = 3.14159265358979323846_wp

  CALL check_jacobians()
  CALL check_dissipation()
  CALL check_straight()
  CALL check_refusals()
  CALL check_pair()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: Hermite mass-proportional damping kernel'

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
    CALL CD_HermiteCable_Mass_Damping(q, v, L, CM, f, Kq, Kv, es, em)
    CALL require(es == CD_HCABLE_OK, 'M1 '//label//': evaluation: '//TRIM(em))
    DO j = 1, 12
      h = 1.0e-6_wp*MAX(1.0_wp, ABS(q(j)))
      dq = 0.0_wp; dq(j) = h
      CALL CD_HermiteCable_Mass_Damping(q + dq, v, L, CM, fp, ErrStat=es, ErrMsg=em)
      CALL CD_HermiteCable_Mass_Damping(q - dq, v, L, CM, fm, ErrStat=es, ErrMsg=em)
      Kq_fd(:, j) = (fp - fm)/(2.0_wp*h)
      h = 1.0e-6_wp*MAX(1.0_wp, ABS(v(j)))
      dq = 0.0_wp; dq(j) = h
      CALL CD_HermiteCable_Mass_Damping(q, v + dq, L, CM, fp, ErrStat=es, ErrMsg=em)
      CALL CD_HermiteCable_Mass_Damping(q, v - dq, L, CM, fm, ErrStat=es, ErrMsg=em)
      Kv_fd(:, j) = (fp - fm)/(2.0_wp*h)
    END DO
    eq = MAXVAL(ABS(Kq - Kq_fd))/MAX(1.0e-12_wp, MAXVAL(ABS(Kq_fd)))
    ev = MAXVAL(ABS(Kv - Kv_fd))/MAX(1.0e-12_wp, MAXVAL(ABS(Kv_fd)))
    WRITE (*, '(A,2ES11.3)') 'M1 '//label//': relative |Kq - FD|, |Kv - FD| ', eq, ev
    CALL require(eq < 1.0e-7_wp, 'M1 '//label//': Kq matches the central difference')
    CALL require(ev < 1.0e-7_wp, 'M1 '//label//': Kv matches the central difference')
    CALL require(MAXVAL(ABS(Kv - TRANSPOSE(Kv))) <= 1.0e-13_wp*MAXVAL(ABS(Kv)), 'M1 '//label//': Kv symmetric')
  END SUBROUTINE fd_check

  SUBROUTINE check_jacobians()
    REAL(wp) :: q(12), v(12)
    CALL curved_state(q, v)
    CALL fd_check(q, v, 'curved 3D')
    q = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, L, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
    CALL fd_check(q, v, 'straight')
  END SUBROUTINE check_jacobians

  REAL(wp) FUNCTION psi(q, v)
    !! 1/2 cm integral (|v|^2 - (t . v)^2) ds on the four-point rule, from the public shapes
    REAL(wp), INTENT(IN) :: q(12), v(12)
    REAL(wp) :: gp(6), gw(6), H(4), dH(4), rp(3), vv(3), t(3)
    INTEGER :: g, c
    CALL CD_HermiteCable_Gauss_Rule(4, gp, gw)
    psi = 0.0_wp
    DO g = 1, 4
      CALL CD_HermiteCable_Shapes(gp(g), L, H, dH)
      DO c = 1, 3
        rp(c) = dH(1)*q(c) + dH(2)*q(3 + c) + dH(3)*q(6 + c) + dH(4)*q(9 + c)
        vv(c) = H(1)*v(c) + H(2)*v(3 + c) + H(3)*v(6 + c) + H(4)*v(9 + c)
      END DO
      t = rp/NORM2(rp)
      psi = psi + 0.5_wp*CM*gw(g)*L*(DOT_PRODUCT(vv, vv) - DOT_PRODUCT(t, vv)**2)
    END DO
  END FUNCTION psi

  SUBROUTINE check_dissipation()
    REAL(wp) :: q(12), v(12), f(12), f2(12), fd(12), dv(12), h, p0
    INTEGER :: es, j
    CHARACTER(256) :: em
    CALL curved_state(q, v)
    CALL CD_HermiteCable_Mass_Damping(q, v, L, CM, f, ErrStat=es, ErrMsg=em)
    CALL require(es == CD_HCABLE_OK, 'M2: evaluation: '//TRIM(em))
    DO j = 1, 12
      h = 1.0e-5_wp
      dv = 0.0_wp; dv(j) = h
      fd(j) = (psi(q, v + dv) - psi(q, v - dv))/(2.0_wp*h)
    END DO
    p0 = psi(q, v)
    WRITE (*, '(A,ES11.3,A,ES11.3)') 'M2: relative |f - d psi/dv| ', MAXVAL(ABS(f - fd))/MAXVAL(ABS(fd)), &
      ', |v . f - 2 psi|/psi ', ABS(DOT_PRODUCT(v, f) - 2.0_wp*p0)/p0
    CALL require(MAXVAL(ABS(f - fd)) <= 1.0e-8_wp*MAXVAL(ABS(fd)), 'M2: the force is d psi / d v')
    CALL require(ABS(DOT_PRODUCT(v, f) - 2.0_wp*p0) <= 1.0e-12_wp*p0 .AND. p0 > 0.0_wp, &
                 'M2: v . f = 2 psi > 0 (dissipative)')
    CALL CD_HermiteCable_Mass_Damping(q, 2.0_wp*v, L, CM, f2, ErrStat=es, ErrMsg=em)
    CALL require(MAXVAL(ABS(f2 - 2.0_wp*f)) <= 1.0e-14_wp*MAXVAL(ABS(f)), 'M2: linear in the rates')
    CALL CD_HermiteCable_Mass_Damping(q, v, L, 0.0_wp, f2, ErrStat=es, ErrMsg=em)
    CALL require(es == CD_HCABLE_OK .AND. ALL(f2 == 0.0_wp), 'M2: zero for cm = 0')
  END SUBROUTINE check_dissipation

  SUBROUTINE check_straight()
    REAL(wp) :: q(12), v(12), f(12), Mm(12, 12), fm(12)
    INTEGER :: es, k
    CHARACTER(256) :: em
    q = [0.0_wp, 0.0_wp, 0.0_wp, 1.01_wp, 0.0_wp, 0.0_wp, 1.01_wp*L, 0.0_wp, 0.0_wp, 1.01_wp, 0.0_wp, 0.0_wp]
    v = [0.3_wp, -0.2_wp, 0.5_wp, 0.1_wp, 0.4_wp, -0.3_wp, -0.2_wp, 0.6_wp, 0.1_wp, 0.2_wp, -0.1_wp, 0.25_wp]
    CALL CD_HermiteCable_Mass_Damping(q, v, L, CM, f, ErrStat=es, ErrMsg=em)
    CALL CD_HermiteCable_Mass(1.0_wp, L, Mm, es, em)
    fm = CM*MATMUL(Mm, v)
    DO k = 0, 3
      fm(3*k + 1) = 0.0_wp
    END DO
    WRITE (*, '(A,ES11.3)') 'M3: straight element, relative |f - cm M v_transverse| ', &
      MAXVAL(ABS(f - fm))/MAXVAL(ABS(fm))
    CALL require(MAXVAL(ABS(f - fm)) <= 1.0e-13_wp*MAXVAL(ABS(fm)), 'M3: cm M v in y and z, zero in x')
    ! axial sliding: every node moves along x, the handles still
    v = [0.7_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    CALL CD_HermiteCable_Mass_Damping(q, v, L, CM, f, ErrStat=es, ErrMsg=em)
    CALL require(ALL(f == 0.0_wp), 'M3: an axial sliding is not damped')
    ! rigid transverse translation: damped as cm M v (the absolute-velocity convention)
    v = [0.0_wp, 0.4_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.4_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
    CALL CD_HermiteCable_Mass_Damping(q, v, L, CM, f, ErrStat=es, ErrMsg=em)
    WRITE (*, '(A,ES11.3,A,ES11.3)') 'M3: rigid transverse translation, total force ', f(2) + f(8), &
      ' N against cm L v ', CM*L*0.4_wp
    CALL require(ABS(f(2) + f(8) - CM*L*0.4_wp) <= 1.0e-13_wp*CM*L*0.4_wp, &
                 'M3: a rigid transverse translation is damped at cm L v')
  END SUBROUTINE check_straight

  SUBROUTINE check_refusals()
    REAL(wp) :: q(12), v(12), f(12)
    INTEGER :: es
    CHARACTER(256) :: em
    CALL curved_state(q, v)
    CALL CD_HermiteCable_Mass_Damping(q, v, L, -1.0_wp, f, ErrStat=es, ErrMsg=em)
    CALL require(es /= CD_HCABLE_OK .AND. INDEX(em, 'non-negative') > 0, 'M4: negative cm refused')
    CALL CD_HermiteCable_Mass_Damping(q, v, 0.0_wp, CM, f, ErrStat=es, ErrMsg=em)
    CALL require(es /= CD_HCABLE_OK .AND. INDEX(em, 'length must be positive') > 0, 'M4: zero length refused')
    v(3) = IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN)
    CALL CD_HermiteCable_Mass_Damping(q, v, L, CM, f, ErrStat=es, ErrMsg=em)
    CALL require(es /= CD_HCABLE_OK .AND. INDEX(em, 'finite') > 0 .AND. ALL(f == 0.0_wp), 'M4: NaN rate refused')
    q = 0.0_wp
    v = 0.0_wp
    CALL CD_HermiteCable_Mass_Damping(q, v, L, CM, f, ErrStat=es, ErrMsg=em)
    CALL require(es /= CD_HCABLE_OK .AND. INDEX(em, 'degenerate') > 0, 'M4: degenerate centreline refused')
  END SUBROUTINE check_refusals

  SUBROUTINE check_pair()
    REAL(wp), PARAMETER :: ZETA = 0.015_wp, F1 = 0.4_wp, F2 = 3.0_wp, EI = 2.0e4_wp
    REAL(wp) :: bb, alpha, z1, z2, zm, fmid
    CALL CD_HermiteCable_Rayleigh_Pair(ZETA, F1, F2, EI, bb, alpha)
    z1 = rayleigh_zeta(alpha, bb/EI, F1)
    z2 = rayleigh_zeta(alpha, bb/EI, F2)
    fmid = SQRT(F1*F2)
    zm = rayleigh_zeta(alpha, bb/EI, fmid)
    WRITE (*, '(A,2ES11.3,A,2ES12.4,A,ES12.4)') 'M5: Bb, alpha ', bb, alpha, '; zeta at f1, f2 ', z1, z2, &
      '; at sqrt(f1 f2) ', zm
    CALL require(ABS(z1 - ZETA) <= 1.0e-14_wp .AND. ABS(z2 - ZETA) <= 1.0e-14_wp, 'M5: zeta at f1 and f2')
    CALL require(ABS(zm - 2.0_wp*ZETA*fmid/(F1 + F2)) <= 1.0e-14_wp, 'M5: the minimum between them')
    CALL require(rayleigh_zeta(alpha, bb/EI, 0.99_wp*fmid) > zm .AND. rayleigh_zeta(alpha, bb/EI, 1.01_wp*fmid) > zm, &
                 'M5: sqrt(f1 f2) is the minimum')
  END SUBROUTINE check_pair

  PURE REAL(wp) FUNCTION rayleigh_zeta(alpha, beta, f)
    !! alpha/(2 omega) + beta omega/2 at omega = 2 pi f
    REAL(wp), INTENT(IN) :: alpha, beta, f
    rayleigh_zeta = alpha/(4.0_wp*PI*f) + beta*PI*f
  END FUNCTION rayleigh_zeta

END PROGRAM test_hermite_mass_damping
