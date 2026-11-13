! File: tests/test_viv_cf_damping_kernel.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_cf_damping_kernel
  !! The split drag of Thorsen et al. (2016) per unit length (CD_Morison_Drag_Split_Per_Length_Jac):
  !!   K1 a point at rest in a steady flow: the split force is Morison's, whatever cdy (w = U e_x);
  !!   K2 a cross-flow velocity y' across the flow: the cross-flow force is
  !!      1/2 rho D cdy |w_y| w_y (w_y = -y'), the in-line force Morison's with the full |w|, and
  !!      the blend beta = U^2/(U^2 + u_ref^2) mixes it with Morison's;
  !!   K3 no ambient normal flow (still water, or a flow along the tangent): Morison's drag, bit for
  !!      bit; and with cdy = 0 and cdy_lin = cdn, the split recombines Morison's drag in a flow;
  !!   K4 the Jacobians against central differences of the force: in the relative velocity and in
  !!      the (non-unit) tangent vector, with tangential relative flow, cdy /= cdn, a Morison-form
  !!      cross-flow term (cdy_lin) and U ~ u_ref;
  !!   K5 bad inputs refused by name.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_Hydro, ONLY: CD_Morison_Drag_Split_Per_Length_Jac, CD_Morison_Drag_Per_Length_Jac, CD_HYDRO_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: RHO = 1025.0_wp, DIAM = 0.1_wp, CDN = 1.2_wp, CDT = 0.3_wp, CDY = 0.55_wp
  INTEGER :: nfail = 0

  CALL check_at_rest()
  CALL check_crossflow()
  CALL check_no_flow()
  CALL check_jacobians()
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: split cross-flow drag kernel'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE check_at_rest()
    REAL(wp) :: amb(3), tg(3), fs(3), fm(3)
    INTEGER :: es
    CHARACTER(200) :: em
    amb = [0.3_wp, 0.8_wp, 0.1_wp]
    tg = [1.0_wp, 0.2_wp, -0.1_wp]
    CALL CD_Morison_Drag_Split_Per_Length_Jac(amb, amb, tg, RHO, DIAM, CDN, CDY, CDT, 0.01_wp, fs, &
                                              ErrStat=es, ErrMsg=em)
    CALL require(es == CD_HYDRO_OK, 'K1: split call: '//TRIM(em))
    CALL CD_Morison_Drag_Per_Length_Jac(amb, tg, RHO, DIAM, CDN, CDT, fm, ErrStat=es, ErrMsg=em)
    WRITE (*, '(A,ES10.2)') 'K1 at rest: |split - Morison|/|Morison| ', NORM2(fs - fm)/NORM2(fm)
    CALL require(NORM2(fs - fm) <= 1.0e-13_wp*NORM2(fm), 'K1: a point at rest feels Morison drag')
  END SUBROUTINE check_at_rest

  SUBROUTINE check_crossflow()
    REAL(wp) :: amb(3), tg(3), v(3), rel(3), f(3), ex(3), ey(3), w(3), wm, beta, fy, fx, uref, u
    INTEGER :: es
    CHARACTER(200) :: em
    u = 0.6_wp
    amb = [0.0_wp, u, 0.0_wp]
    tg = [2.0_wp, 0.0_wp, 0.0_wp]
    ex = [0.0_wp, 1.0_wp, 0.0_wp]
    ey = [0.0_wp, 0.0_wp, 1.0_wp]        ! t x e_x = x x y = z
    v = [0.0_wp, 0.1_wp, 0.35_wp]        ! in-line 0.1, cross-flow 0.35
    rel = amb - v
    uref = 0.2_wp
    CALL CD_Morison_Drag_Split_Per_Length_Jac(rel, amb, tg, RHO, DIAM, CDN, CDY, CDT, uref, f, ErrStat=es, ErrMsg=em)
    CALL require(es == CD_HYDRO_OK, 'K2: split call: '//TRIM(em))
    w = [0.0_wp, rel(2), rel(3)]
    wm = NORM2(w)
    beta = u*u/(u*u + uref*uref)
    fy = beta*0.5_wp*RHO*DIAM*CDY*ABS(rel(3))*rel(3) + (1.0_wp - beta)*0.5_wp*RHO*DIAM*CDN*wm*rel(3)
    fx = 0.5_wp*RHO*DIAM*CDN*wm*rel(2)
    WRITE (*, '(A,2ES12.4,A,2ES12.4)') 'K2 cross-flow: force (in line, across) ', DOT_PRODUCT(f, ex), &
      DOT_PRODUCT(f, ey), ' expected ', fx, fy
    CALL require(ABS(DOT_PRODUCT(f, ey) - fy) <= 1.0e-13_wp*ABS(fy), &
                 'K2: the cross-flow force is the blended cdy |w_y| w_y')
    CALL require(ABS(DOT_PRODUCT(f, ex) - fx) <= 1.0e-13_wp*ABS(fx), 'K2: the in-line force is Morison''s')
    CALL require(ABS(f(1)) <= 1.0e-12_wp, 'K2: no tangential force without tangential flow')
  END SUBROUTINE check_crossflow

  SUBROUTINE check_no_flow()
    REAL(wp) :: tg(3), rel(3), fs(3), fm(3), amb(3)
    INTEGER :: es
    CHARACTER(200) :: em
    tg = [0.3_wp, -1.0_wp, 0.4_wp]
    rel = [0.2_wp, 0.1_wp, -0.3_wp]
    amb = 0.0_wp
   CALL CD_Morison_Drag_Split_Per_Length_Jac(rel, amb, tg, RHO, DIAM, CDN, CDY, CDT, 0.01_wp, fs, ErrStat=es, ErrMsg=em)
    CALL CD_Morison_Drag_Per_Length_Jac(rel, tg, RHO, DIAM, CDN, CDT, fm, ErrStat=es, ErrMsg=em)
    CALL require(ALL(fs == fm), 'K3: no ambient flow, Morison drag bit for bit')
    amb = 0.7_wp*tg
   CALL CD_Morison_Drag_Split_Per_Length_Jac(rel, amb, tg, RHO, DIAM, CDN, CDY, CDT, 0.01_wp, fs, ErrStat=es, ErrMsg=em)
    CALL require(ALL(fs == fm), 'K3: a flow along the tangent, Morison drag bit for bit')
    amb = [0.1_wp, 0.6_wp, -0.2_wp]
    CALL CD_Morison_Drag_Split_Per_Length_Jac(rel, amb, tg, RHO, DIAM, CDN, 0.0_wp, CDT, 0.01_wp, fs, ErrStat=es, &
                                              ErrMsg=em, cdy_lin=CDN)
    WRITE (*, '(A,ES10.2)') 'K3 cdy 0, cdy_lin cdn: |split - Morison|/|Morison| ', NORM2(fs - fm)/NORM2(fm)
    CALL require(es == CD_HYDRO_OK .AND. NORM2(fs - fm) <= 1.0e-13_wp*NORM2(fm), &
                 'K3: cdy 0 and cdy_lin cdn recombine Morison''s drag')
  END SUBROUTINE check_no_flow

  SUBROUTINE check_jacobians()
    REAL(wp) :: rel(3), amb(3), tg(3), f(3), fp(3), fmn(3), jr(3, 3), jt(3, 3), h, err_r, err_t, scale, d(3)
    INTEGER :: es, k
    CHARACTER(200) :: em
    rel = [0.21_wp, 0.13_wp, -0.37_wp]
    amb = [0.05_wp, 0.18_wp, 0.02_wp]
    tg = [1.3_wp, -0.4_wp, 0.25_wp]
    CALL CD_Morison_Drag_Split_Per_Length_Jac(rel, amb, tg, RHO, DIAM, CDN, CDY, CDT, 0.12_wp, f, jr, jt, es, em, &
                                              cdy_lin=0.7_wp)
    CALL require(es == CD_HYDRO_OK, 'K4: split call with Jacobians: '//TRIM(em))
    h = 1.0e-6_wp
    err_r = 0.0_wp; err_t = 0.0_wp
    scale = MAX(MAXVAL(ABS(jr)), MAXVAL(ABS(jt)))
    DO k = 1, 3
      d = 0.0_wp; d(k) = h
      CALL CD_Morison_Drag_Split_Per_Length_Jac(rel + d, amb, tg, RHO, DIAM, CDN, CDY, CDT, 0.12_wp, fp, ErrStat=es, &
                                                ErrMsg=em, cdy_lin=0.7_wp)
      CALL CD_Morison_Drag_Split_Per_Length_Jac(rel - d, amb, tg, RHO, DIAM, CDN, CDY, CDT, 0.12_wp, fmn, ErrStat=es, &
                                                ErrMsg=em, cdy_lin=0.7_wp)
      err_r = MAX(err_r, MAXVAL(ABS((fp - fmn)/(2.0_wp*h) - jr(:, k))))
      CALL CD_Morison_Drag_Split_Per_Length_Jac(rel, amb, tg + d, RHO, DIAM, CDN, CDY, CDT, 0.12_wp, fp, ErrStat=es, &
                                                ErrMsg=em, cdy_lin=0.7_wp)
      CALL CD_Morison_Drag_Split_Per_Length_Jac(rel, amb, tg - d, RHO, DIAM, CDN, CDY, CDT, 0.12_wp, fmn, ErrStat=es, &
                                                ErrMsg=em, cdy_lin=0.7_wp)
      err_t = MAX(err_t, MAXVAL(ABS((fp - fmn)/(2.0_wp*h) - jt(:, k))))
    END DO
    WRITE (*, '(A,2ES10.2,A,ES10.2)') 'K4 Jacobians against central differences: rel, tangent ', err_r, err_t, &
      ' (scale ', scale, ')'
    CALL require(err_r <= 1.0e-7_wp*scale .AND. err_t <= 1.0e-7_wp*scale, 'K4: exact Jacobians')
  END SUBROUTINE check_jacobians

  SUBROUTINE check_refusals()
    REAL(wp) :: f(3)
    INTEGER :: es
    CHARACTER(200) :: em
    CALL CD_Morison_Drag_Split_Per_Length_Jac([0.1_wp, 0.0_wp, 0.0_wp], [0.0_wp, 0.5_wp, 0.0_wp], &
                                              [1.0_wp, 0.0_wp, 0.0_wp], &
                                              RHO, DIAM, CDN, -0.1_wp, CDT, 0.01_wp, f, ErrStat=es, ErrMsg=em)
    CALL require(es /= CD_HYDRO_OK .AND. INDEX(em, 'non-negative') > 0, 'K5: a negative cdy refused')
    CALL CD_Morison_Drag_Split_Per_Length_Jac([0.1_wp, 0.0_wp, 0.0_wp], [0.0_wp, 0.5_wp, 0.0_wp], &
                                              [1.0_wp, 0.0_wp, 0.0_wp], &
                                              RHO, DIAM, CDN, CDY, CDT, 0.0_wp, f, ErrStat=es, ErrMsg=em)
    CALL require(es /= CD_HYDRO_OK .AND. INDEX(em, 'u_ref positive') > 0, 'K5: a zero u_ref refused')
  END SUBROUTINE check_refusals

END PROGRAM test_viv_cf_damping_kernel
