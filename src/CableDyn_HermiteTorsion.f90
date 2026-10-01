! File: src/CableDyn_HermiteTorsion.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_HermiteTorsion
  !! Parallel-transport twist of a cubic-Hermite cable line, for condensed isotropic torsion.
  !!
  !! The Hermite cable (CableDyn_HermiteCable) carries [r(3), m(3)] per node, m = dr/ds, and has
  !! no twist coordinate. For an isotropic section without distributed torque the torque is
  !! uniform along a line, and the twist field condenses to one scalar per line:
  !!   E_t = (Phi - Theta(q))**2 / (2 C),   C = sum_e L_e / GJ_e,   M_t = (Phi - Theta) / C,
  !! where Phi is the imposed relative roll of the end frames and Theta(q) is the parallel-
  !! transport (smallest-rotation) holonomy of the centreline: the angle about d_B, right-hand
  !! rule, from the End-B reference normal n_B to the End-A normal n_A carried along the line.
  !! Because the Hermite curve is C1, Theta is the sum of
  !!   * one element holonomy per element (12 element DOFs),
  !!       h_e = + integral_e t1 . (a x b) / (|a|**2 (1 + t1 . a/|a|)) ds,
  !!     a = r', b = r'', t1 = m1/|m1| (Gauss-Legendre quadrature, four points by default);
  !!   * the smallest-rotation chain d_A -> t_0 -> ... -> t_N -> d_B through the unit node
  !!     tangents t_k = m_k/|m_k|.
  !! The value of the chain is evaluated by sequential composition. Its first and second
  !! derivatives are formed as a sum of local terms: every link j (u_j -> u_j+1) is given a gauge
  !! with a pole p_j = unit(u_j + u_j+1) frozen at the current configuration, so that the chain
  !! equals, up to a locally constant multiple of 2 pi,
  !!   sum_links Omega(p_j, u_j, u_j+1) - sum_vertices Omega(p_j-1, p_j, u_j) + end terms,
  !! with Omega(x, y, z) = 2 atan2(x . (y x z), 1 + x.y + y.z + z.x) the signed solid angle of a
  !! geodesic triangle. Every term depends on at most two neighbouring tangents, so the Hessian
  !! has the half-bandwidth 11 of the cable tangent. Gradient and Hessian are closed form.
  !!
  !! The end frames (d, n) enter through incremental rotations exp(w) about the current frame,
  !! w = 0: the routine also returns dTheta/dw and the w-w and w-m second derivatives.
  !!
  !! Theta is known modulo 2 pi. CD_HermiteTorsion_Accept selects the branch nearest to the
  !! previous accepted value and rejects changes larger than pi/2, so a 2 pi slip cannot pass
  !! unnoticed. All routines fail closed with named status codes: folds (1 + t1 . t < 0.5 at a
  !! Gauss point or across a chain link), non-finite values, and malformed input.
  USE CableDyn_Precision, ONLY: wp, CD_ZERO, CD_ONE, CD_All_Finite, CD_Is_Finite
  USE CableDyn_HermiteCable, ONLY: CD_HermiteCable_Gauss_Rule
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CD_HermiteTorsion_Element
  PUBLIC :: CD_HermiteTorsion_Line
  PUBLIC :: CD_HermiteTorsion_Fold_Check
  PUBLIC :: CD_HermiteTorsion_Unwrap
  PUBLIC :: CD_HermiteTorsion_Accept

  INTEGER, PARAMETER, PUBLIC :: CD_HTORS_OK = 0
  INTEGER, PARAMETER, PUBLIC :: CD_HTORS_BADINPUT = 1
  INTEGER, PARAMETER, PUBLIC :: CD_HTORS_FOLD = 2
  INTEGER, PARAMETER, PUBLIC :: CD_HTORS_NONFINITE = 3
  INTEGER, PARAMETER, PUBLIC :: CD_HTORS_STEP = 4
  !! Half-bandwidth of d2Theta/dq2 (equal to that of the cable tangent).
  INTEGER, PARAMETER, PUBLIC :: CD_HTORS_KBAND = 11
  !! Fold guard: 1 + t1 . t at every Gauss point and 1 + u_j . u_j+1 on every chain link.
  REAL(wp), PARAMETER, PUBLIC :: CD_HTORS_FOLD_MIN = 0.5_wp
  REAL(wp), PARAMETER, PUBLIC :: CD_HTORS_PI = 3.14159265358979323846264338327950288_wp
  !! Largest accepted change of Theta between two accepted states.
  REAL(wp), PARAMETER, PUBLIC :: CD_HTORS_MAX_STEP = 0.5_wp*CD_HTORS_PI

  INTEGER, PARAMETER :: NG_DEFAULT = 4, NG_MAX = 6
  REAL(wp), PARAMETER :: TWO_PI = 2.0_wp*CD_HTORS_PI
  REAL(wp), PARAMETER :: GEOM_TOL = 1.0e-12_wp
  REAL(wp), PARAMETER :: FRAME_TOL = 1.0e-8_wp

CONTAINS

  ! ------------------------------------------------------------------------------------------
  ! public routines
  ! ------------------------------------------------------------------------------------------

  PURE SUBROUTINE CD_HermiteTorsion_Element(qe, L, h, grad, ErrStat, ErrMsg, hess, fold_margin, quadrature_order)
    !! Element holonomy h_e of one 12-DOF element [r1, m1, r2, m2] of reference length L, with
    !! its closed-form gradient and (optionally) Hessian. fold_margin returns min(1 + t1 . t)
    !! over the Gauss points; below CD_HTORS_FOLD_MIN the routine fails with CD_HTORS_FOLD.
    REAL(wp), INTENT(IN) :: qe(12), L
    REAL(wp), INTENT(OUT) :: h, grad(12)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(OUT), OPTIONAL :: hess(12, 12)
    REAL(wp), INTENT(OUT), OPTIONAL :: fold_margin
    INTEGER, INTENT(IN), OPTIONAL :: quadrature_order

    REAL(wp) :: hk(12, 12), margin
    INTEGER :: order
    CHARACTER(200) :: em

    order = NG_DEFAULT
    IF (PRESENT(quadrature_order)) order = quadrature_order
    h = CD_ZERO
    grad = CD_ZERO
    IF (PRESENT(hess)) hess = CD_ZERO
    IF (PRESENT(fold_margin)) fold_margin = CD_ZERO
    IF (order < 1 .OR. order > NG_MAX) THEN
      ErrStat = CD_HTORS_BADINPUT
      ErrMsg = 'CD_HermiteTorsion_Element: quadrature order must lie in [1,6]'
      RETURN
    END IF
    IF (.NOT. (CD_All_Finite(qe) .AND. CD_Is_Finite(L))) THEN
      ErrStat = CD_HTORS_NONFINITE
      ErrMsg = 'CD_HermiteTorsion_Element: non-finite element DOFs or length'
      RETURN
    END IF
    IF (L <= GEOM_TOL) THEN
      ErrStat = CD_HTORS_BADINPUT
      ErrMsg = 'CD_HermiteTorsion_Element: element length must be positive'
      RETURN
    END IF
    CALL element_holonomy(qe, L, order, PRESENT(hess), h, grad, hk, margin, ErrStat, em)
    IF (PRESENT(fold_margin)) fold_margin = margin
    IF (ErrStat /= CD_HTORS_OK) THEN
      ErrMsg = 'CD_HermiteTorsion_Element: '//em
      h = CD_ZERO
      grad = CD_ZERO
      RETURN
    END IF
    IF (PRESENT(hess)) hess = hk
  END SUBROUTINE CD_HermiteTorsion_Element

  PURE SUBROUTINE CD_HermiteTorsion_Fold_Check(q, Le, ends, margin, ErrStat, ErrMsg, quadrature_order)
    !! Fold guard for one line: margin = min(1 + t1 . t) over every Gauss point of every element
    !! and min(1 + u_j . u_j+1) over every chain link d_A -> t_0 -> ... -> t_N -> d_B. Fails with
    !! CD_HTORS_FOLD (naming the element or link) when margin < CD_HTORS_FOLD_MIN, i.e. when the
    !! tangent turns by more than 120 degrees inside one element, between two nodes, or across an
    !! articulated end. The holonomy and its derivatives grow without bound as the margin -> 0.
    REAL(wp), INTENT(IN) :: q(:), Le(:), ends(3, 4)
    REAL(wp), INTENT(OUT) :: margin
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: quadrature_order
    INTEGER :: order

    order = NG_DEFAULT
    IF (PRESENT(quadrature_order)) order = quadrature_order
    CALL line_fold_check(q, Le, ends, order, 'CD_HermiteTorsion_Fold_Check', margin, ErrStat, ErrMsg)
  END SUBROUTINE CD_HermiteTorsion_Fold_Check

  PURE SUBROUTINE line_fold_check(q, Le, ends, order, who, margin, ErrStat, ErrMsg)
    REAL(wp), INTENT(IN) :: q(:), Le(:), ends(3, 4)
    INTEGER, INTENT(IN) :: order
    CHARACTER(*), INTENT(IN) :: who
    REAL(wp), INTENT(OUT) :: margin
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: nn, e, j, ig
    REAL(wp) :: gp(NG_MAX), gw(NG_MAX), t1(3), a(3), b(3), rho, sp, val
    REAL(wp), ALLOCATABLE :: u(:, :)
    CHARACTER(16) :: tag

    margin = CD_ZERO
    CALL validate_line(q, Le, ends, order, who, ErrStat, ErrMsg)
    IF (ErrStat /= CD_HTORS_OK) RETURN
    nn = SIZE(q)/6
    ALLOCATE (u(3, 0:nn + 1))
    CALL chain_vertices(q, ends, u)
    CALL CD_HermiteCable_Gauss_Rule(order, gp, gw)
    margin = HUGE(CD_ONE)
    DO e = 1, nn - 1
      rho = NORM2(q(6*e - 2:6*e))
      t1 = q(6*e - 2:6*e)/rho
      DO ig = 1, order
        CALL element_ab(q(6*e - 5:6*e + 6), Le(e), gp(ig), a, b)
        sp = NORM2(a)
        IF (sp <= GEOM_TOL) THEN
          ErrStat = CD_HTORS_BADINPUT
          WRITE (tag, '(I0)') e
          ErrMsg = who//': zero centreline tangent in element '//TRIM(tag)
          RETURN
        END IF
        val = CD_ONE + DOT_PRODUCT(t1, a)/sp
        IF (val < margin) margin = val
        IF (val < CD_HTORS_FOLD_MIN) THEN
          ErrStat = CD_HTORS_FOLD
          WRITE (tag, '(I0)') e
          ErrMsg = who//': tangent turns by more than 120 degrees inside element ' &
                   //TRIM(tag)//' (1 + t1.t < 0.5); refine the mesh'
          RETURN
        END IF
      END DO
    END DO
    DO j = 0, nn
      val = CD_ONE + DOT_PRODUCT(u(:, j), u(:, j + 1))
      IF (val < margin) margin = val
      IF (val < CD_HTORS_FOLD_MIN) THEN
        ErrStat = CD_HTORS_FOLD
        IF (j == 0) THEN
          ErrMsg = who//': End-A director and first node tangent differ by more than ' &
                   //'120 degrees (1 + d.t < 0.5)'
        ELSE IF (j == nn) THEN
          ErrMsg = who//': last node tangent and End-B director differ by more than ' &
                   //'120 degrees (1 + t.d < 0.5)'
        ELSE
          WRITE (tag, '(I0)') j
          ErrMsg = who//': node tangents turn by more than 120 degrees across element ' &
                   //TRIM(tag)//' (1 + t1.t2 < 0.5); refine the mesh'
        END IF
        RETURN
      END IF
    END DO
    ErrStat = CD_HTORS_OK
    ErrMsg = ''
  END SUBROUTINE line_fold_check

  PURE SUBROUTINE CD_HermiteTorsion_Line(q, Le, ends, theta_raw, grad, ErrStat, ErrMsg, hband, band_scale, &
                                         end_grad, end_hess, end_cross, fold_margin, quadrature_order)
    !! Theta, dTheta/dq and d2Theta/dq2 of one line.
    !!
    !!   q(6 n)          nodal DOFs [r_k, m_k], k = 1..n (n >= 2), as in the cable solvers
    !!   Le(n - 1)       element reference lengths
    !!   ends(3, 4)      columns d_A, n_A, d_B, n_B: unit end directors and unit reference
    !!                   normals (n perpendicular to d). A clamped end has d = its end tangent.
    !!   theta_raw       sum of element holonomies + chain angle in (-pi, pi]; Theta is defined
    !!                   modulo 2 pi -- pass it through CD_HermiteTorsion_Accept
    !!   grad(6 n)       dTheta/dq
    !!   hband           optional, INOUT: band_scale * d2Theta/dq2 is ADDED in general-band
    !!                   storage with half-bandwidth CD_HTORS_KBAND: entry (i, j) at
    !!                   hband(SIZE(hband, 1) - KBAND + i - j, j). This is the LAPACK DGBSV
    !!                   layout for SIZE(hband, 1) = 3 KBAND + 1 and the DGBMV layout for
    !!                   2 KBAND + 1. band_scale defaults to 1.
    !!   end_grad(6)     dTheta/d[w_A, w_B], w the rotation vector of an incremental rotation
    !!                   exp(w) of the end frame (d and n rotate together), at w = 0
    !!   end_hess(3,3,2) d2Theta/dw_A2 and d2Theta/dw_B2 (same parametrisation)
    !!   end_cross(3,3,2) d2Theta/dw_A dm_1 and d2Theta/dw_B dm_n (rows w, columns m). Theta has
    !!                   no other w coupling (no w_A-w_B, no w-r).
    !!   fold_margin     min(1 + t1 . t) over Gauss points and chain links
    !!
    !! On any failure every output is zero and hband is left unchanged.
    REAL(wp), INTENT(IN) :: q(:), Le(:), ends(3, 4)
    REAL(wp), INTENT(OUT) :: theta_raw, grad(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(INOUT), OPTIONAL :: hband(:, :)
    REAL(wp), INTENT(IN), OPTIONAL :: band_scale
    REAL(wp), INTENT(OUT), OPTIONAL :: end_grad(6), end_hess(3, 3, 2), end_cross(3, 3, 2)
    REAL(wp), INTENT(OUT), OPTIONAL :: fold_margin
    INTEGER, INTENT(IN), OPTIONAL :: quadrature_order

    INTEGER :: order, nn, n, e, j, k, i0, i, jj, row0
    LOGICAL :: want_hess
    REAL(wp) :: scale, margin, he, ge(12), hk(12, 12), v(3), chi, sum_h, eg(6), eh(3, 3, 2), ec(3, 3, 2)
    REAL(wp) :: om, g9(9), h9(9, 9), pj(3, 3), pk(3, 3), sd(3, 3), dmy
    REAL(wp), ALLOCATABLE :: u(:, :), pole(:, :), gu(:, :), huu(:, :, :), huv(:, :, :), band(:, :)
    CHARACTER(200) :: em

    theta_raw = CD_ZERO
    grad = CD_ZERO
    IF (PRESENT(end_grad)) end_grad = CD_ZERO
    IF (PRESENT(end_hess)) end_hess = CD_ZERO
    IF (PRESENT(end_cross)) end_cross = CD_ZERO
    IF (PRESENT(fold_margin)) fold_margin = CD_ZERO
    order = NG_DEFAULT
    IF (PRESENT(quadrature_order)) order = quadrature_order
    scale = CD_ONE
    IF (PRESENT(band_scale)) scale = band_scale
    want_hess = PRESENT(hband) .OR. PRESENT(end_hess) .OR. PRESENT(end_cross)

    CALL line_fold_check(q, Le, ends, order, 'CD_HermiteTorsion_Line', margin, ErrStat, ErrMsg)
    IF (ErrStat /= CD_HTORS_OK) RETURN
    n = SIZE(q)
    nn = n/6
    IF (SIZE(grad) /= n) THEN
      ErrStat = CD_HTORS_BADINPUT
      ErrMsg = 'CD_HermiteTorsion_Line: grad must have SIZE(q) entries'
      RETURN
    END IF
    IF (PRESENT(hband)) THEN
      IF (SIZE(hband, 1) < 2*CD_HTORS_KBAND + 1 .OR. SIZE(hband, 2) /= n) THEN
        ErrStat = CD_HTORS_BADINPUT
        ErrMsg = 'CD_HermiteTorsion_Line: hband must be at least (2*KBAND+1, SIZE(q))'
        RETURN
      END IF
    END IF
    IF (.NOT. CD_Is_Finite(scale)) THEN
      ErrStat = CD_HTORS_NONFINITE
      ErrMsg = 'CD_HermiteTorsion_Line: non-finite band_scale'
      RETURN
    END IF

    ALLOCATE (u(3, 0:nn + 1), pole(3, 0:nn), gu(3, 0:nn + 1), huu(3, 3, 0:nn + 1), huv(3, 3, 0:nn))
    ALLOCATE (band(2*CD_HTORS_KBAND + 1, n))
    row0 = CD_HTORS_KBAND + 1          ! diagonal row of the local band
    band = CD_ZERO
    CALL chain_vertices(q, ends, u)

    ! element holonomies
    sum_h = CD_ZERO
    DO e = 1, nn - 1
      i0 = 6*(e - 1)
      CALL element_holonomy(q(i0 + 1:i0 + 12), Le(e), order, want_hess, he, ge, hk, dmy, ErrStat, em)
      IF (ErrStat /= CD_HTORS_OK) THEN
        ErrMsg = 'CD_HermiteTorsion_Line: '//em
        grad = CD_ZERO
        RETURN
      END IF
      sum_h = sum_h + he
      grad(i0 + 1:i0 + 12) = grad(i0 + 1:i0 + 12) + ge
      IF (want_hess) THEN
        DO jj = 1, 12
          DO i = 1, 12
            band(row0 + i - jj, i0 + jj) = band(row0 + i - jj, i0 + jj) + hk(i, jj)
          END DO
        END DO
      END IF
    END DO

    ! smallest-rotation chain: value by sequential composition
    v = ends(:, 2)
    DO j = 0, nn
      v = smallest_rotation(u(:, j), u(:, j + 1), v)
    END DO
    chi = ATAN2(DOT_PRODUCT(ends(:, 3), cross(ends(:, 4), v)), DOT_PRODUCT(ends(:, 4), v))
    theta_raw = sum_h + chi

    ! chain derivatives: local solid-angle terms with frozen poles, w.r.t. the raw unit vertices
    gu = CD_ZERO
    huu = CD_ZERO
    huv = CD_ZERO
    DO j = 0, nn
      pole(:, j) = u(:, j) + u(:, j + 1)
      pole(:, j) = pole(:, j)/NORM2(pole(:, j))
    END DO
    DO j = 0, nn                                   ! links: + Omega(p_j, u_j, u_j+1)
      CALL solid_angle(pole(:, j), u(:, j), u(:, j + 1), want_hess, om, g9, h9)
      gu(:, j) = gu(:, j) + g9(4:6)
      gu(:, j + 1) = gu(:, j + 1) + g9(7:9)
      IF (want_hess) THEN
        huu(:, :, j) = huu(:, :, j) + h9(4:6, 4:6)
        huu(:, :, j + 1) = huu(:, :, j + 1) + h9(7:9, 7:9)
        huv(:, :, j) = huv(:, :, j) + h9(4:6, 7:9)
      END IF
    END DO
    DO j = 1, nn                                   ! interior vertices: - Omega(p_j-1, p_j, u_j)
      CALL solid_angle(pole(:, j - 1), pole(:, j), u(:, j), want_hess, om, g9, h9)
      gu(:, j) = gu(:, j) - g9(7:9)
      IF (want_hess) huu(:, :, j) = huu(:, :, j) - h9(7:9, 7:9)
    END DO
    ! end vertices: gauge change to a pole at the frozen end director, plus the end-frame twist
    CALL solid_angle(ends(:, 1), pole(:, 0), u(:, 0), want_hess, om, g9, h9)
    gu(:, 0) = gu(:, 0) - g9(7:9)
    IF (want_hess) huu(:, :, 0) = huu(:, :, 0) - h9(7:9, 7:9)
    CALL solid_angle(pole(:, nn), ends(:, 3), u(:, nn + 1), want_hess, om, g9, h9)
    gu(:, nn + 1) = gu(:, nn + 1) - g9(7:9)
    IF (want_hess) huu(:, :, nn + 1) = huu(:, :, nn + 1) - h9(7:9, 7:9)

    ! chain to the node tangent handles m_k (u = m/|m|)
    DO k = 1, nn
      i0 = 6*(k - 1) + 3
      CALL normalisation_jacobian(q(i0 + 1:i0 + 3), pj)
      grad(i0 + 1:i0 + 3) = grad(i0 + 1:i0 + 3) + MATMUL(pj, gu(:, k))
      IF (want_hess) THEN
        sd = normalised_hessian(q(i0 + 1:i0 + 3), gu(:, k), huu(:, :, k))
        CALL add_block(band, row0, i0, i0, sd)
        IF (k < nn) THEN
          CALL normalisation_jacobian(q(i0 + 7:i0 + 9), pk)
          sd = MATMUL(pj, MATMUL(huv(:, :, k), pk))
          CALL add_block(band, row0, i0, i0 + 6, sd)
          sd = TRANSPOSE(sd)
          CALL add_block(band, row0, i0 + 6, i0, sd)
        END IF
      END IF
    END DO

    ! chain to the end-frame rotations: d(w) = exp(w) d, dd/dw = -S(d), second order 1/2 w x (w x d)
    eg = CD_ZERO
    eh = CD_ZERO
    ec = CD_ZERO
    eg(1:3) = cross(ends(:, 1), gu(:, 0)) + ends(:, 1)
    eg(4:6) = cross(ends(:, 3), gu(:, nn + 1)) - ends(:, 3)
    IF (want_hess) THEN
      eh(:, :, 1) = rotation_hessian(ends(:, 1), gu(:, 0), huu(:, :, 0))
      eh(:, :, 2) = rotation_hessian(ends(:, 3), gu(:, nn + 1), huu(:, :, nn + 1))
      CALL normalisation_jacobian(q(4:6), pj)
      sd = skew(ends(:, 1))
      ec(:, :, 1) = MATMUL(sd, MATMUL(huv(:, :, 0), pj))
      CALL normalisation_jacobian(q(n - 2:n), pk)
      sd = skew(ends(:, 3))
      ec(:, :, 2) = MATMUL(sd, MATMUL(TRANSPOSE(huv(:, :, nn)), pk))
    END IF

    IF (.NOT. (CD_Is_Finite(theta_raw) .AND. CD_All_Finite(grad) .AND. CD_All_Finite(band) .AND. &
               CD_All_Finite(eg) .AND. CD_All_Finite(RESHAPE(eh, [18])) .AND. CD_All_Finite(RESHAPE(ec, [18])))) THEN
      ErrStat = CD_HTORS_NONFINITE
      ErrMsg = 'CD_HermiteTorsion_Line: non-finite twist or derivative'
      theta_raw = CD_ZERO
      grad = CD_ZERO
      RETURN
    END IF

    IF (PRESENT(hband)) THEN
      i0 = SIZE(hband, 1) - 2*CD_HTORS_KBAND - 1
      DO jj = 1, n
        DO i = MAX(1, jj - CD_HTORS_KBAND), MIN(n, jj + CD_HTORS_KBAND)
          hband(i0 + row0 + i - jj, jj) = hband(i0 + row0 + i - jj, jj) + scale*band(row0 + i - jj, jj)
        END DO
      END DO
    END IF
    IF (PRESENT(end_grad)) end_grad = eg
    IF (PRESENT(end_hess)) end_hess = eh
    IF (PRESENT(end_cross)) end_cross = ec
    IF (PRESENT(fold_margin)) fold_margin = margin
    ErrStat = CD_HTORS_OK
    ErrMsg = ''
  END SUBROUTINE CD_HermiteTorsion_Line

  ELEMENTAL REAL(wp) FUNCTION CD_HermiteTorsion_Unwrap(theta_raw, theta_prev) RESULT(theta)
    !! The 2 pi branch of theta_raw nearest to theta_prev: theta_raw + 2 pi round((prev - raw)/2 pi).
    !! Non-finite input returns theta_raw unchanged (the caller's finiteness check rejects it).
    REAL(wp), INTENT(IN) :: theta_raw, theta_prev
    theta = theta_raw
    IF (CD_Is_Finite(theta_raw) .AND. CD_Is_Finite(theta_prev)) &
      theta = theta_raw + TWO_PI*ANINT((theta_prev - theta_raw)/TWO_PI)
  END FUNCTION CD_HermiteTorsion_Unwrap

  PURE SUBROUTINE CD_HermiteTorsion_Accept(theta_prev, theta_raw, theta_new, ErrStat, ErrMsg, max_step)
    !! Unwrap a trial value against the previous accepted Theta and enforce the step limit
    !! |theta_new - theta_prev| <= max_step (default pi/2, a margin below the hard limit pi at
    !! which the nearest branch becomes ambiguous). On rejection theta_new = theta_prev and
    !! ErrStat = CD_HTORS_STEP: the caller must cut the step. Call it after every accepted
    !! Newton iterate and every accepted time step or load stage.
    REAL(wp), INTENT(IN) :: theta_prev, theta_raw
    REAL(wp), INTENT(OUT) :: theta_new
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: max_step
    REAL(wp) :: limit

    limit = CD_HTORS_MAX_STEP
    IF (PRESENT(max_step)) limit = max_step
    theta_new = theta_prev
    IF (.NOT. (CD_Is_Finite(theta_prev) .AND. CD_Is_Finite(theta_raw) .AND. CD_Is_Finite(limit))) THEN
      ErrStat = CD_HTORS_NONFINITE
      ErrMsg = 'CD_HermiteTorsion_Accept: non-finite twist value or step limit'
      RETURN
    END IF
    IF (limit <= CD_ZERO .OR. limit >= CD_HTORS_PI) THEN
      ErrStat = CD_HTORS_BADINPUT
      ErrMsg = 'CD_HermiteTorsion_Accept: step limit must lie in (0, pi)'
      RETURN
    END IF
    IF (ABS(CD_HermiteTorsion_Unwrap(theta_raw, theta_prev) - theta_prev) > limit) THEN
      ErrStat = CD_HTORS_STEP
      ErrMsg = 'CD_HermiteTorsion_Accept: twist change exceeds the step limit; cut the step'
      RETURN
    END IF
    theta_new = CD_HermiteTorsion_Unwrap(theta_raw, theta_prev)
    ErrStat = CD_HTORS_OK
    ErrMsg = ''
  END SUBROUTINE CD_HermiteTorsion_Accept

  ! ------------------------------------------------------------------------------------------
  ! element holonomy
  ! ------------------------------------------------------------------------------------------

  PURE SUBROUTINE element_ab(qe, L, xi, a, b)
    !! a = dr/ds and b = d2r/ds2 of the cubic-Hermite element at xi in [0, 1].
    REAL(wp), INTENT(IN) :: qe(12), L, xi
    REAL(wp), INTENT(OUT) :: a(3), b(3)
    REAL(wp) :: dh(4), ddh(4)
    INTEGER :: i
    CALL shape_derivatives(L, xi, dh, ddh)
    a = CD_ZERO
    b = CD_ZERO
    DO i = 1, 4
      a = a + dh(i)*qe(3*i - 2:3*i)
      b = b + ddh(i)*qe(3*i - 2:3*i)
    END DO
  END SUBROUTINE element_ab

  PURE SUBROUTINE shape_derivatives(L, x, dh, ddh)
    !! d/ds and d2/ds2 of the Hermite shape functions [H1, L H2, H3, L H4] at xi = x.
    REAL(wp), INTENT(IN) :: L, x
    REAL(wp), INTENT(OUT) :: dh(4), ddh(4)
    dh = [-6.0_wp*x + 6.0_wp*x*x, L*(CD_ONE - 4.0_wp*x + 3.0_wp*x*x), 6.0_wp*x - 6.0_wp*x*x, &
          L*(-2.0_wp*x + 3.0_wp*x*x)]/L
    ddh = [-6.0_wp + 12.0_wp*x, L*(-4.0_wp + 6.0_wp*x), 6.0_wp - 12.0_wp*x, L*(-2.0_wp + 6.0_wp*x)]/(L*L)
  END SUBROUTINE shape_derivatives

  PURE SUBROUTINE element_holonomy(qe, L, order, want_hess, h, grad, hess, margin, ErrStat, ErrMsg)
    !! h_e = L sum_g w_g f(a_g, b_g, t1), f = t1.(a x b) / (|a|**2 + |a| (t1.a)), t1 = m1/|m1|.
    !! f = N/D is differentiated in closed form in x = (a, b, t1) and chained to the 12 DOFs by
    !! the constant maps a, b = sum_i dh_i q_i, sum_i ddh_i q_i and by t1 = m1/|m1|.
    REAL(wp), INTENT(IN) :: qe(12), L
    INTEGER, INTENT(IN) :: order
    LOGICAL, INTENT(IN) :: want_hess
    REAL(wp), INTENT(OUT) :: h, grad(12), hess(12, 12), margin
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp) :: gp(NG_MAX), gw(NG_MAX), dh(4), ddh(4), m1(3), rho, t(3), a(3), b(3), s2, sp, ta
    REAL(wp) :: num, den, f, gn(9), gd(9), hn(9, 9), hd(9, 9), gf(9), hf(9, 9), w, eye(3, 3)
    REAL(wp) :: gt(3), htt(3, 3), xat(3, 3, 4), pm(3, 3), blk(3, 3)
    INTEGER :: ig, i, j

    h = CD_ZERO
    grad = CD_ZERO
    hess = CD_ZERO
    margin = CD_ZERO
    m1 = qe(4:6)
    rho = NORM2(m1)
    IF (rho <= GEOM_TOL) THEN
      ErrStat = CD_HTORS_BADINPUT
      ErrMsg = 'element holonomy: zero tangent handle at the first node'
      RETURN
    END IF
    t = m1/rho
    eye = identity3()
    gt = CD_ZERO
    htt = CD_ZERO
    xat = CD_ZERO
    margin = HUGE(CD_ONE)
    CALL CD_HermiteCable_Gauss_Rule(order, gp, gw)
    DO ig = 1, order
      CALL shape_derivatives(L, gp(ig), dh, ddh)
      CALL element_ab(qe, L, gp(ig), a, b)
      s2 = DOT_PRODUCT(a, a)
      sp = SQRT(s2)
      IF (sp <= GEOM_TOL) THEN
        ErrStat = CD_HTORS_BADINPUT
        ErrMsg = 'element holonomy: zero centreline tangent at a Gauss point'
        RETURN
      END IF
      ta = DOT_PRODUCT(t, a)
      margin = MIN(margin, CD_ONE + ta/sp)
      IF (CD_ONE + ta/sp < CD_HTORS_FOLD_MIN) THEN
        ErrStat = CD_HTORS_FOLD
        ErrMsg = 'element holonomy: tangent turns by more than 120 degrees inside the element (1 + t1.t < 0.5)'
        RETURN
      END IF
      num = DOT_PRODUCT(t, cross(a, b))
      den = s2 + sp*ta
      f = num/den
      w = L*gw(ig)
      h = h + w*f
      gn(1:3) = cross(b, t)
      gn(4:6) = cross(t, a)
      gn(7:9) = cross(a, b)
      gd(1:3) = 2.0_wp*a + (ta/sp)*a + sp*t
      gd(4:6) = CD_ZERO
      gd(7:9) = sp*a
      gf = (gn - f*gd)/den
      DO i = 1, 4
        grad(3*i - 2:3*i) = grad(3*i - 2:3*i) + w*(dh(i)*gf(1:3) + ddh(i)*gf(4:6))
      END DO
      gt = gt + w*gf(7:9)
      IF (.NOT. want_hess) CYCLE
      hn = CD_ZERO
      hn(1:3, 4:6) = -skew(t)
      hn(4:6, 1:3) = skew(t)
      hn(4:6, 7:9) = -skew(a)
      hn(7:9, 4:6) = skew(a)
      hn(7:9, 1:3) = -skew(b)
      hn(1:3, 7:9) = skew(b)
      hd = CD_ZERO
      hd(1:3, 1:3) = 2.0_wp*eye + ta*(eye/sp - outer(a, a)/(sp*s2)) + (outer(a, t) + outer(t, a))/sp
      hd(1:3, 7:9) = outer(a, a)/sp + sp*eye
      hd(7:9, 1:3) = TRANSPOSE(hd(1:3, 7:9))
      hf = (hn - f*hd - outer9(gf, gd) - outer9(gd, gf))/den
      DO j = 1, 4
        DO i = 1, 4
          blk = dh(i)*dh(j)*hf(1:3, 1:3) + dh(i)*ddh(j)*hf(1:3, 4:6) + ddh(i)*dh(j)*hf(4:6, 1:3) &
                + ddh(i)*ddh(j)*hf(4:6, 4:6)
          hess(3*i - 2:3*i, 3*j - 2:3*j) = hess(3*i - 2:3*i, 3*j - 2:3*j) + w*blk
        END DO
        xat(:, :, j) = xat(:, :, j) + w*(dh(j)*hf(1:3, 7:9) + ddh(j)*hf(4:6, 7:9))
      END DO
      htt = htt + w*hf(7:9, 7:9)
    END DO
    ! chain through t1 = m1/|m1| (DOF block 2)
    CALL normalisation_jacobian(m1, pm)
    grad(4:6) = grad(4:6) + MATMUL(pm, gt)
    IF (want_hess) THEN
      DO i = 1, 4
        blk = MATMUL(xat(:, :, i), pm)
        hess(3*i - 2:3*i, 4:6) = hess(3*i - 2:3*i, 4:6) + blk
        hess(4:6, 3*i - 2:3*i) = hess(4:6, 3*i - 2:3*i) + TRANSPOSE(blk)
      END DO
      hess(4:6, 4:6) = hess(4:6, 4:6) + normalised_hessian(m1, gt, htt)
    END IF
    ErrStat = CD_HTORS_OK
    ErrMsg = ''
  END SUBROUTINE element_holonomy

  ! ------------------------------------------------------------------------------------------
  ! chain helpers
  ! ------------------------------------------------------------------------------------------

  PURE SUBROUTINE validate_line(q, Le, ends, order, who, ErrStat, ErrMsg)
    REAL(wp), INTENT(IN) :: q(:), Le(:), ends(3, 4)
    INTEGER, INTENT(IN) :: order
    CHARACTER(*), INTENT(IN) :: who
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: nn, k

    ErrStat = CD_HTORS_BADINPUT
    IF (order < 1 .OR. order > NG_MAX) THEN
      ErrMsg = who//': quadrature order must lie in [1,6]'
      RETURN
    END IF
    nn = SIZE(q)/6
    IF (MOD(SIZE(q), 6) /= 0 .OR. nn < 2 .OR. SIZE(Le) /= nn - 1) THEN
      ErrMsg = who//': q must hold 6 DOFs for each of n >= 2 nodes and Le n - 1 lengths'
      RETURN
    END IF
    IF (.NOT. (CD_All_Finite(q) .AND. CD_All_Finite(Le) .AND. CD_All_Finite(ends))) THEN
      ErrStat = CD_HTORS_NONFINITE
      ErrMsg = who//': non-finite DOFs, lengths or end frames'
      RETURN
    END IF
    IF (ANY(Le <= GEOM_TOL)) THEN
      ErrMsg = who//': element lengths must be positive'
      RETURN
    END IF
    DO k = 1, nn
      IF (NORM2(q(6*k - 2:6*k)) <= GEOM_TOL) THEN
        ErrMsg = who//': zero tangent handle m at a node'
        RETURN
      END IF
    END DO
    DO k = 1, 3, 2
      IF (ABS(NORM2(ends(:, k)) - CD_ONE) > FRAME_TOL .OR. ABS(NORM2(ends(:, k + 1)) - CD_ONE) > FRAME_TOL &
          .OR. ABS(DOT_PRODUCT(ends(:, k), ends(:, k + 1))) > FRAME_TOL) THEN
        IF (k == 1) THEN
          ErrMsg = who//': End-A director and reference normal must be orthonormal'
        ELSE
          ErrMsg = who//': End-B director and reference normal must be orthonormal'
        END IF
        RETURN
      END IF
    END DO
    ErrStat = CD_HTORS_OK
    ErrMsg = ''
  END SUBROUTINE validate_line

  PURE SUBROUTINE chain_vertices(q, ends, u)
    !! u(:, 0) = d_A, u(:, k) = m_k/|m_k| (k = 1..n), u(:, n+1) = d_B.
    REAL(wp), INTENT(IN) :: q(:), ends(3, 4)
    REAL(wp), INTENT(OUT) :: u(:, 0:)
    INTEGER :: nn, k
    nn = SIZE(q)/6
    u(:, 0) = ends(:, 1)
    DO k = 1, nn
      u(:, k) = q(6*k - 2:6*k)/NORM2(q(6*k - 2:6*k))
    END DO
    u(:, nn + 1) = ends(:, 3)
  END SUBROUTINE chain_vertices

  PURE FUNCTION smallest_rotation(a, b, v) RESULT(r)
    !! Smallest rotation taking unit a to unit b, applied to v (requires 1 + a.b > 0).
    REAL(wp), INTENT(IN) :: a(3), b(3), v(3)
    REAL(wp) :: r(3), c, w(3)
    c = DOT_PRODUCT(a, b)
    w = cross(a, b)
    r = c*v + cross(w, v) + w*DOT_PRODUCT(w, v)/(CD_ONE + c)
  END FUNCTION smallest_rotation

  PURE SUBROUTINE solid_angle(x, y, z, want_hess, om, g, h)
    !! Omega = 2 atan2(N, D), N = x.(y x z), D = 1 + x.y + y.z + z.x (signed solid angle of the
    !! geodesic triangle x, y, z), with gradient g and Hessian h in the raw vectors (x, y, z).
    !! Under the fold guard every triangle used here has D >= 1.5.
    REAL(wp), INTENT(IN) :: x(3), y(3), z(3)
    LOGICAL, INTENT(IN) :: want_hess
    REAL(wp), INTENT(OUT) :: om, g(9), h(9, 9)
    REAL(wp) :: nv, dv, r2, gn(9), gd(9), eye(3, 3)

    nv = DOT_PRODUCT(x, cross(y, z))
    dv = CD_ONE + DOT_PRODUCT(x, y) + DOT_PRODUCT(y, z) + DOT_PRODUCT(z, x)
    om = 2.0_wp*ATAN2(nv, dv)
    gn(1:3) = cross(y, z)
    gn(4:6) = cross(z, x)
    gn(7:9) = cross(x, y)
    gd(1:3) = y + z
    gd(4:6) = x + z
    gd(7:9) = x + y
    r2 = nv*nv + dv*dv
    g = 2.0_wp*(dv*gn - nv*gd)/r2
    h = CD_ZERO
    IF (.NOT. want_hess) RETURN
    eye = identity3()
    ! d2N: (x,y) = -S(z), (y,z) = -S(x), (z,x) = -S(y); d2D: identity off-diagonal blocks
    h(1:3, 4:6) = -dv*skew(z) - nv*eye
    h(4:6, 1:3) = dv*skew(z) - nv*eye
    h(4:6, 7:9) = -dv*skew(x) - nv*eye
    h(7:9, 4:6) = dv*skew(x) - nv*eye
    h(7:9, 1:3) = -dv*skew(y) - nv*eye
    h(1:3, 7:9) = dv*skew(y) - nv*eye
    h = 2.0_wp*(h/r2 + ((nv*nv - dv*dv)*(outer9(gn, gd) + outer9(gd, gn)) &
                        - 2.0_wp*dv*nv*(outer9(gn, gn) - outer9(gd, gd)))/(r2*r2))
  END SUBROUTINE solid_angle

  PURE SUBROUTINE normalisation_jacobian(m, p)
    !! d(m/|m|)/dm = (I - t t^T)/|m|.
    REAL(wp), INTENT(IN) :: m(3)
    REAL(wp), INTENT(OUT) :: p(3, 3)
    REAL(wp) :: rho, t(3)
    rho = NORM2(m)
    t = m/rho
    p = (identity3() - outer(t, t))/rho
  END SUBROUTINE normalisation_jacobian

  PURE FUNCTION normalised_hessian(m, g, h) RESULT(hm)
    !! Hessian in m of F(t), t = m/|m|, from g = dF/dt and h = d2F/dt2:
    !! P h P - (g t^T + t g^T + (g.t)(I - 3 t t^T))/|m|**2,  P = (I - t t^T)/|m|.
    REAL(wp), INTENT(IN) :: m(3), g(3), h(3, 3)
    REAL(wp) :: hm(3, 3), p(3, 3), rho, t(3), gt
    rho = NORM2(m)
    t = m/rho
    CALL normalisation_jacobian(m, p)
    gt = DOT_PRODUCT(g, t)
    hm = MATMUL(p, MATMUL(h, p)) - (outer(g, t) + outer(t, g) + gt*(identity3() - 3.0_wp*outer(t, t)))/(rho*rho)
  END FUNCTION normalised_hessian

  PURE FUNCTION rotation_hessian(d, g, h) RESULT(hw)
    !! Hessian in w at w = 0 of F(exp(w) d), from g = dF/dd and h = d2F/dd2:
    !! -S(d) h S(d) + (g d^T + d g^T)/2 - (g.d) I.
    REAL(wp), INTENT(IN) :: d(3), g(3), h(3, 3)
    REAL(wp) :: hw(3, 3), sd(3, 3)
    sd = skew(d)
    hw = -MATMUL(sd, MATMUL(h, sd)) + 0.5_wp*(outer(g, d) + outer(d, g)) - DOT_PRODUCT(g, d)*identity3()
  END FUNCTION rotation_hessian

  PURE SUBROUTINE add_block(band, row0, i0, j0, blk)
    !! Add a 3x3 block at global rows i0+1..i0+3, columns j0+1..j0+3 (local band, diagonal row0).
    REAL(wp), INTENT(INOUT) :: band(:, :)
    INTEGER, INTENT(IN) :: row0, i0, j0
    REAL(wp), INTENT(IN) :: blk(3, 3)
    INTEGER :: i, j
    DO j = 1, 3
      DO i = 1, 3
        band(row0 + i0 + i - j0 - j, j0 + j) = band(row0 + i0 + i - j0 - j, j0 + j) + blk(i, j)
      END DO
    END DO
  END SUBROUTINE add_block

  ! ------------------------------------------------------------------------------------------
  ! small vector helpers
  ! ------------------------------------------------------------------------------------------

  PURE FUNCTION cross(a, b) RESULT(c)
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  END FUNCTION cross

  PURE FUNCTION skew(v) RESULT(s)
    !! S(v) w = v x w.
    REAL(wp), INTENT(IN) :: v(3)
    REAL(wp) :: s(3, 3)
    s = RESHAPE([CD_ZERO, v(3), -v(2), -v(3), CD_ZERO, v(1), v(2), -v(1), CD_ZERO], [3, 3])
  END FUNCTION skew

  PURE FUNCTION outer(a, b) RESULT(m)
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: m(3, 3)
    m = SPREAD(a, 2, 3)*SPREAD(b, 1, 3)
  END FUNCTION outer

  PURE FUNCTION outer9(a, b) RESULT(m)
    REAL(wp), INTENT(IN) :: a(9), b(9)
    REAL(wp) :: m(9, 9)
    m = SPREAD(a, 2, 9)*SPREAD(b, 1, 9)
  END FUNCTION outer9

  PURE FUNCTION identity3() RESULT(m)
    REAL(wp) :: m(3, 3)
    m = RESHAPE([CD_ONE, CD_ZERO, CD_ZERO, CD_ZERO, CD_ONE, CD_ZERO, CD_ZERO, CD_ZERO, CD_ONE], [3, 3])
  END FUNCTION identity3

END MODULE CableDyn_HermiteTorsion
