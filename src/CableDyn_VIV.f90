! File: src/CableDyn_VIV.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_VIV
  !! Cross-flow vortex-induced vibration (VIV) of EI = 0 line nodes: the synchronisation
  !! model of Thorsen, Saevik and Larsen (J. Fluids Struct. 49, 2014, 135-148; 61, 2016,
  !! 1-19) in the form MoorDyn v2 implements it (LINE TYPES columns Cl dF cF; Davies et al.,
  !! J. Fluids Struct. 146, 2026, 104651). Pure node kernels; the model layer
  !! (CableDyn_Model) owns the phase state and the time integration.
  !!
  !! At an interior node with unit tangent q (the chord of its two neighbours), velocity r',
  !! acceleration r'' and ambient fluid velocity u:
  !!
  !!   v_r = u - r'              relative fluid velocity
  !!   v_n = v_r - (q.v_r) q     its component normal to the line, U = |v_n|
  !!   e   = q x v_n / U         cross-flow unit vector
  !!   y'  = r'.e ,  y'' = r''.e cross-flow velocity and acceleration of the node
  !!
  !! The cross-flow velocity phase is phi_y = atan2(-y''/rms(y''), y'/rms(y')) in [0, 2 pi),
  !! with running RMS values. The vortex-shedding (lift) phase phi obeys
  !!
  !!   d phi/dt = 2 pi f_hat U / D ,   f_hat = cF + dF sin(phi_y - phi)
  !!
  !! so the shedding frequency locks onto the cross-flow motion inside the non-dimensional
  !! synchronisation band cF +- dF, and the lift force on the node is
  !!
  !!   F = C cos(phi) U (q x v_n) ,   C = sum over the node's two elements of
  !!                                      (1/2 l0_e w_e) (1/2 rho D_e Cl_e)
  !!
  !! (w_e the submerged fraction of element e): per unit length 1/2 rho D Cl U^2 cos(phi)
  !! along e, exactly MoorDyn's Lf = 1/2 rho_w d |Vp| Cl cos(phi) (q x Vp) l_node.
  !!
  !! In-line VIV (Ulveseter, Thorsen, Saevik and Larsen, Mar. Struct. 60, 2018, 241-260,
  !! Eq. 1, 2, 6 and 7): a second shedding force along the relative normal flow,
  !!
  !!   F_x = C_x cos(phi_x) U v_n ,   per unit length 1/2 rho D C_vx U v_n cos(phi_x)
  !!
  !! whose phase is slaved to twice the cross-flow shedding rate with a free band alpha,
  !!
  !!   d phi_x/dt = 2 (d phi/dt) [1 + alpha sin(phi_xd - phi_x)]
  !!
  !! phi_xd the phase of the in-line velocity x' = r'.v_n/U from x', x'' normalised by their
  !! running RMS values, as phi_y.
  !!
  !! Load convention (CableDyn_Hydro, CD_Cable_Dynamic_Load_Proc): `force` is the external
  !! load; the Jacobians are residual contributions, -dF/d(node velocity) and
  !! -dF/d(neighbour positions) through the node tangent. The wetted-fraction factor in C is
  !! held at its current value in the Jacobians (it moves only for an element that crosses
  !! the waterline); the residual itself is exact.
  USE CableDyn_Precision, ONLY: wp, CD_ZERO, CD_ONE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CD_VIV_Node_Lift
  PUBLIC :: CD_VIV_Lift_Per_Length
  PUBLIC :: CD_VIV_Crossflow_Motion
  PUBLIC :: CD_VIV_Crossflow_Along
  PUBLIC :: CD_VIV_Phase_Rate
  PUBLIC :: CD_VIV_Mean_Square
  PUBLIC :: CD_VIV_Wrap_Phase
  PUBLIC :: CD_VIV_Node_Inline
  PUBLIC :: CD_VIV_Inline_Per_Length
  PUBLIC :: CD_VIV_Inline_Phase_Rate
  PUBLIC :: CD_VIV_Wake_Advance

  !> MoorDyn's running-RMS window in time steps (Line%n_m): each step updates the stored
  !> mean square as ((n_m - 1) old + new) / n_m.
  INTEGER, PARAMETER, PUBLIC :: CD_VIV_NM = 500
  !> MoorDyn's defaults of the synchronisation range dF and centre cF when a LINE TYPES row
  !> gives Cl only.
  REAL(wp), PARAMETER, PUBLIC :: CD_VIV_DEFAULT_DF = 0.08_wp
  REAL(wp), PARAMETER, PUBLIC :: CD_VIV_DEFAULT_CF = 0.18_wp
  REAL(wp), PARAMETER, PUBLIC :: CD_VIV_TWO_PI = 6.283185307179586476925286766559005768394_wp
  !> The in-line frequency band alpha of Ulveseter et al. (2018, Table 3) and Thorsen et al.
  REAL(wp), PARAMETER, PUBLIC :: CD_VIV_DEFAULT_IL_ALPHA = 0.15_wp
  !> The wake oscillator of Facchinetti, de Langre and Biolley (J. Fluids Struct. 19, 2004,
  !> 123-140, Sec. 3): Strouhal number, reference lift coefficient C_L0, van der Pol parameter
  !> epsilon and acceleration coupling A, fixed there from forced-wake data.
  REAL(wp), PARAMETER, PUBLIC :: CD_VIV_WAKE_DEFAULT_ST = 0.2_wp
  REAL(wp), PARAMETER, PUBLIC :: CD_VIV_WAKE_DEFAULT_CL0 = 0.3_wp
  REAL(wp), PARAMETER, PUBLIC :: CD_VIV_WAKE_DEFAULT_EPS = 0.3_wp
  REAL(wp), PARAMETER, PUBLIC :: CD_VIV_WAKE_DEFAULT_A = 12.0_wp

  !> Chord below which a node tangent is undefined (a collapsed pair of neighbours).
  REAL(wp), PARAMETER :: TINY_CHORD = 1.0e-12_wp

CONTAINS

  PURE SUBROUTINE node_frame(r_prev, r_next, rd, u, q, chord, vr, vn, un, ok)
    !! Unit tangent from the neighbour chord and the normal relative velocity of a node.
    REAL(wp), INTENT(IN) :: r_prev(3), r_next(3), rd(3), u(3)
    REAL(wp), INTENT(OUT) :: q(3), chord, vr(3), vn(3), un
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp) :: c(3)

    c = r_next - r_prev
    chord = SQRT(c(1)*c(1) + c(2)*c(2) + c(3)*c(3))
    ok = chord > TINY_CHORD
    IF (.NOT. ok) THEN
      q = CD_ZERO
      vr = CD_ZERO
      vn = CD_ZERO
      un = CD_ZERO
      RETURN
    END IF
    q = c/chord
    vr = u - rd
    vn = vr - (q(1)*vr(1) + q(2)*vr(2) + q(3)*vr(3))*q
    un = SQRT(vn(1)*vn(1) + vn(2)*vn(2) + vn(3)*vn(3))
  END SUBROUTINE node_frame

  PURE FUNCTION cross3(a, b) RESULT(c)
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  END FUNCTION cross3

  PURE SUBROUTINE CD_VIV_Node_Lift(r_prev, r_next, rd, u, coef, phi, force, ok, jac_v, jac_prev, jac_next)
    !! Lift force on one node, F = coef cos(phi) U (q x v_n), and optionally its residual
    !! Jacobians (all three present or none): jac_v = -dF/dr', jac_prev = -dF/dr_prev,
    !! jac_next = -dF/dr_next. ok = .FALSE. (and F = 0) when the two neighbours coincide.
    REAL(wp), INTENT(IN) :: r_prev(3), r_next(3), rd(3), u(3), coef, phi
    REAL(wp), INTENT(OUT) :: force(3)
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp), INTENT(OUT), OPTIONAL :: jac_v(3, 3), jac_prev(3, 3), jac_next(3, 3)
    REAL(wp) :: q(3), chord, vr(3), vn(3), un, w(3), s, alpha
    REAL(wp) :: dfdq(3, 3), dfdv(3, 3), proj(3, 3)
    INTEGER :: i, j

    force = CD_ZERO
    IF (PRESENT(jac_v)) THEN
      jac_v = CD_ZERO
      jac_prev = CD_ZERO
      jac_next = CD_ZERO
    END IF
    CALL node_frame(r_prev, r_next, rd, u, q, chord, vr, vn, un, ok)
    IF (.NOT. ok) RETURN
    s = coef*COS(phi)
    w = cross3(q, vn)
    force = (s*un)*w
    IF (.NOT. PRESENT(jac_v)) RETURN
    IF (.NOT. (un > CD_ZERO)) RETURN   ! F and its derivatives vanish with the normal flow
    ! w = q x v_n = q x v_r. dF/dv_r = s [w (x) v_n/U + U [q]x]; dr' enters v_r with a minus
    ! sign, so jac_v = -dF/dr' = +dF/dv_r.
    DO j = 1, 3
      DO i = 1, 3
        dfdv(i, j) = s*w(i)*vn(j)/un
      END DO
    END DO
    dfdv(1, 2) = dfdv(1, 2) - s*un*q(3)
    dfdv(1, 3) = dfdv(1, 3) + s*un*q(2)
    dfdv(2, 1) = dfdv(2, 1) + s*un*q(3)
    dfdv(2, 3) = dfdv(2, 3) - s*un*q(1)
    dfdv(3, 1) = dfdv(3, 1) - s*un*q(2)
    dfdv(3, 2) = dfdv(3, 2) + s*un*q(1)
    jac_v = dfdv
    ! dF/dq = s [-(alpha/U) w (x) v_r - U [v_r]x], alpha = q.v_r (U^2 = |v_r|^2 - alpha^2 and
    ! d(q x v_r)/dq = -[v_r]x), chained through dq/dr_next = (I - q q^T)/chord = -dq/dr_prev.
    alpha = q(1)*vr(1) + q(2)*vr(2) + q(3)*vr(3)
    DO j = 1, 3
      DO i = 1, 3
        dfdq(i, j) = -s*(alpha/un)*w(i)*vr(j)
      END DO
    END DO
    dfdq(1, 2) = dfdq(1, 2) + s*un*vr(3)
    dfdq(1, 3) = dfdq(1, 3) - s*un*vr(2)
    dfdq(2, 1) = dfdq(2, 1) - s*un*vr(3)
    dfdq(2, 3) = dfdq(2, 3) + s*un*vr(1)
    dfdq(3, 1) = dfdq(3, 1) + s*un*vr(2)
    dfdq(3, 2) = dfdq(3, 2) - s*un*vr(1)
    DO j = 1, 3
      DO i = 1, 3
        proj(i, j) = -q(i)*q(j)/chord
      END DO
      proj(j, j) = proj(j, j) + CD_ONE/chord
    END DO
    jac_next = -MATMUL(dfdq, proj)
    jac_prev = -jac_next
  END SUBROUTINE CD_VIV_Node_Lift

  PURE SUBROUTINE CD_VIV_Crossflow_Motion(r_prev, r_next, rd, rdd, u, yd, ydd, un, ok, xd, xdd)
    !! Cross-flow velocity y' and acceleration y'' of a node along e = q x v_n / U, and the
    !! normal relative flow speed U. Both are zero when there is no normal flow (MoorDyn's
    !! normalisation of a zero vector). Optionally (both or neither) the in-line velocity x' and
    !! acceleration x'' along v_n/U, zero in the same cases.
    REAL(wp), INTENT(IN) :: r_prev(3), r_next(3), rd(3), rdd(3), u(3)
    REAL(wp), INTENT(OUT) :: yd, ydd, un
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp), INTENT(OUT), OPTIONAL :: xd, xdd
    REAL(wp) :: q(3), chord, vr(3), vn(3), e(3)

    yd = CD_ZERO
    ydd = CD_ZERO
    IF (PRESENT(xd)) THEN
      xd = CD_ZERO
      xdd = CD_ZERO
    END IF
    CALL node_frame(r_prev, r_next, rd, u, q, chord, vr, vn, un, ok)
    IF (.NOT. ok) RETURN
    IF (.NOT. (un > CD_ZERO)) RETURN
    e = cross3(q, vn)/un
    yd = rd(1)*e(1) + rd(2)*e(2) + rd(3)*e(3)
    ydd = rdd(1)*e(1) + rdd(2)*e(2) + rdd(3)*e(3)
    IF (PRESENT(xd)) THEN
      xd = (rd(1)*vn(1) + rd(2)*vn(2) + rd(3)*vn(3))/un
      xdd = (rdd(1)*vn(1) + rdd(2)*vn(2) + rdd(3)*vn(3))/un
    END IF
  END SUBROUTINE CD_VIV_Crossflow_Motion

  PURE SUBROUTINE CD_VIV_Lift_Per_Length(urel, rp, k, f, ok, jrel, jtan)
    !! Lift per unit length of a cable centreline point (the finite-EI cable's Gauss points),
    !! f = k U (t x v_r), with k = 1/2 rho D Cl cos(phi), the unit tangent t = rp/|rp| of the raw
    !! centreline derivative rp = dr/dxi, the relative flow v_r = urel and its normal speed U.
    !! Optionally (both or neither) its Jacobians jrel = df/dv_r and jtan = df/drp, for chaining
    !! through the shape functions as the drag primitive is. ok = .FALSE. (and f = 0) for a
    !! degenerate centreline; f and its derivatives vanish with the normal flow.
    REAL(wp), INTENT(IN) :: urel(3), rp(3), k
    REAL(wp), INTENT(OUT) :: f(3)
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp), INTENT(OUT), OPTIONAL :: jrel(3, 3), jtan(3, 3)
    REAL(wp) :: speed, t(3), alpha, vn(3), un, w(3), dfdt(3, 3), proj(3, 3)
    INTEGER :: i, j

    f = CD_ZERO
    IF (PRESENT(jrel)) THEN
      jrel = CD_ZERO
      jtan = CD_ZERO
    END IF
    speed = SQRT(rp(1)*rp(1) + rp(2)*rp(2) + rp(3)*rp(3))
    ok = speed > TINY_CHORD
    IF (.NOT. ok) RETURN
    t = rp/speed
    alpha = t(1)*urel(1) + t(2)*urel(2) + t(3)*urel(3)
    vn = urel - alpha*t
    un = SQRT(vn(1)*vn(1) + vn(2)*vn(2) + vn(3)*vn(3))
    w = cross3(t, urel)                ! = t x v_n
    f = (k*un)*w
    IF (.NOT. PRESENT(jrel)) RETURN
    IF (.NOT. (un > CD_ZERO)) RETURN
    ! df/dv_r = k [w (x) v_n/U + U [t]x]
    DO j = 1, 3
      DO i = 1, 3
        jrel(i, j) = k*w(i)*vn(j)/un
      END DO
    END DO
    jrel(1, 2) = jrel(1, 2) - k*un*t(3)
    jrel(1, 3) = jrel(1, 3) + k*un*t(2)
    jrel(2, 1) = jrel(2, 1) + k*un*t(3)
    jrel(2, 3) = jrel(2, 3) - k*un*t(1)
    jrel(3, 1) = jrel(3, 1) - k*un*t(2)
    jrel(3, 2) = jrel(3, 2) + k*un*t(1)
    ! df/dt = k [-(alpha/U) w (x) v_r - U [v_r]x], chained through dt/drp = (I - t t^T)/|rp|
    DO j = 1, 3
      DO i = 1, 3
        dfdt(i, j) = -k*(alpha/un)*w(i)*urel(j)
      END DO
    END DO
    dfdt(1, 2) = dfdt(1, 2) + k*un*urel(3)
    dfdt(1, 3) = dfdt(1, 3) - k*un*urel(2)
    dfdt(2, 1) = dfdt(2, 1) - k*un*urel(3)
    dfdt(2, 3) = dfdt(2, 3) + k*un*urel(1)
    dfdt(3, 1) = dfdt(3, 1) + k*un*urel(2)
    dfdt(3, 2) = dfdt(3, 2) - k*un*urel(1)
    DO j = 1, 3
      DO i = 1, 3
        proj(i, j) = -t(i)*t(j)/speed
      END DO
      proj(j, j) = proj(j, j) + CD_ONE/speed
    END DO
    jtan = MATMUL(dfdt, proj)
  END SUBROUTINE CD_VIV_Lift_Per_Length

  PURE SUBROUTINE CD_VIV_Crossflow_Along(tdir, rd, rdd, u, yd, ydd, un, ok, xd, xdd)
    !! CD_VIV_Crossflow_Motion with the line tangent given directly (tdir, any length: the
    !! slope DOFs of a finite-EI node) instead of the chord of two neighbours, with the same
    !! optional in-line velocity and acceleration (both or neither).
    REAL(wp), INTENT(IN) :: tdir(3), rd(3), rdd(3), u(3)
    REAL(wp), INTENT(OUT) :: yd, ydd, un
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp), INTENT(OUT), OPTIONAL :: xd, xdd
    REAL(wp) :: speed, t(3), vr(3), vn(3), e(3)

    yd = CD_ZERO
    ydd = CD_ZERO
    un = CD_ZERO
    IF (PRESENT(xd)) THEN
      xd = CD_ZERO
      xdd = CD_ZERO
    END IF
    speed = SQRT(tdir(1)*tdir(1) + tdir(2)*tdir(2) + tdir(3)*tdir(3))
    ok = speed > TINY_CHORD
    IF (.NOT. ok) RETURN
    t = tdir/speed
    vr = u - rd
    vn = vr - (t(1)*vr(1) + t(2)*vr(2) + t(3)*vr(3))*t
    un = SQRT(vn(1)*vn(1) + vn(2)*vn(2) + vn(3)*vn(3))
    IF (.NOT. (un > CD_ZERO)) RETURN
    e = cross3(t, vn)/un
    yd = rd(1)*e(1) + rd(2)*e(2) + rd(3)*e(3)
    ydd = rdd(1)*e(1) + rdd(2)*e(2) + rdd(3)*e(3)
    IF (PRESENT(xd)) THEN
      xd = (rd(1)*vn(1) + rd(2)*vn(2) + rd(3)*vn(3))/un
      xdd = (rdd(1)*vn(1) + rdd(2)*vn(2) + rdd(3)*vn(3))/un
    END IF
  END SUBROUTINE CD_VIV_Crossflow_Along

  ELEMENTAL REAL(wp) FUNCTION CD_VIV_Mean_Square(ms_old, sample, weight) RESULT(ms)
    !! One step of the running mean square, ms_old + weight (sample^2 - ms_old); weight =
    !! 1/n_m is MoorDyn's ((n_m - 1) old^2 + new^2)/n_m.
    REAL(wp), INTENT(IN) :: ms_old, sample, weight
    ms = ms_old + weight*(sample*sample - ms_old)
  END FUNCTION CD_VIV_Mean_Square

  PURE REAL(wp) FUNCTION CD_VIV_Phase_Rate(yd, ydd, ms_yd, ms_ydd, phi, un, diameter, cf, df, phi_rest) RESULT(rate)
    !! Rate of the vortex-shedding phase, d phi/dt = 2 pi (cF + dF sin(phi_y - phi)) U/D, with
    !! the cross-flow velocity phase phi_y from y', y'' normalised by their running RMS values
    !! (unnormalised while either mean square is still zero, as MoorDyn). A node with no
    !! cross-flow motion takes phi_y = phi_rest (default 0, MoorDyn's): pi on a line whose node
    !! tangent runs against MoorDyn's, which turns the cross-flow axis and so phi_y half a turn.
    REAL(wp), INTENT(IN) :: yd, ydd, ms_yd, ms_ydd, phi, un, diameter, cf, df
    REAL(wp), INTENT(IN), OPTIONAL :: phi_rest
    REAL(wp) :: phi_y, rest

    rest = CD_ZERO
    IF (PRESENT(phi_rest)) rest = phi_rest
    phi_y = motion_phase(yd, ydd, ms_yd, ms_ydd, rest)
    rate = CD_VIV_TWO_PI*(cf + df*SIN(phi_y - phi))*un/diameter
  END FUNCTION CD_VIV_Phase_Rate

  PURE REAL(wp) FUNCTION CD_VIV_Inline_Phase_Rate(xd, xdd, ms_xd, ms_xdd, phi_x, rate_y, alpha) RESULT(rate)
    !! Rate of the in-line shedding phase, d phi_x/dt = 2 rate_y [1 + alpha sin(phi_xd - phi_x)]
    !! (Ulveseter et al. 2018, Eq. 6), rate_y the cross-flow shedding rate (CD_VIV_Phase_Rate) and
    !! phi_xd the phase of the in-line velocity from x', x'' as the cross-flow one; a node with no
    !! in-line motion takes phi_xd = 0.
    REAL(wp), INTENT(IN) :: xd, xdd, ms_xd, ms_xdd, phi_x, rate_y, alpha
    rate = 2.0_wp*rate_y*(CD_ONE + alpha*SIN(motion_phase(xd, xdd, ms_xd, ms_xdd, CD_ZERO) - phi_x))
  END FUNCTION CD_VIV_Inline_Phase_Rate

  PURE REAL(wp) FUNCTION motion_phase(sd, sdd, ms_sd, ms_sdd, rest) RESULT(ph)
    !! The phase of a velocity in [0, 2 pi), atan2(-s''/rms(s''), s'/rms(s')), unnormalised while
    !! either mean square is still zero (MoorDyn), and rest for a node with no motion.
    REAL(wp), INTENT(IN) :: sd, sdd, ms_sd, ms_sdd, rest
    IF (.NOT. (ABS(sd) > CD_ZERO .OR. ABS(sdd) > CD_ZERO)) THEN
      ph = rest
    ELSE IF (ms_sd > CD_ZERO .AND. ms_sdd > CD_ZERO) THEN
      ph = ATAN2(-sdd/SQRT(ms_sdd), sd/SQRT(ms_sd))
    ELSE
      ph = ATAN2(-sdd, sd)
    END IF
    IF (ph < CD_ZERO) ph = ph + CD_VIV_TWO_PI
  END FUNCTION motion_phase

  PURE SUBROUTINE CD_VIV_Node_Inline(r_prev, r_next, rd, u, coef, phi_x, force, ok, jac_v, jac_prev, jac_next)
    !! In-line shedding force on one node, F = coef cos(phi_x) U v_n, and optionally its residual
    !! Jacobians (all three present or none) in CD_VIV_Node_Lift's convention: jac_v = -dF/dr',
    !! jac_prev = -dF/dr_prev, jac_next = -dF/dr_next. ok = .FALSE. (and F = 0) when the two
    !! neighbours coincide.
    REAL(wp), INTENT(IN) :: r_prev(3), r_next(3), rd(3), u(3), coef, phi_x
    REAL(wp), INTENT(OUT) :: force(3)
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp), INTENT(OUT), OPTIONAL :: jac_v(3, 3), jac_prev(3, 3), jac_next(3, 3)
    REAL(wp) :: q(3), chord, vr(3), vn(3), un, s, dfdv(3, 3), dfdq(3, 3), proj(3, 3)
    INTEGER :: i, j

    force = CD_ZERO
    IF (PRESENT(jac_v)) THEN
      jac_v = CD_ZERO
      jac_prev = CD_ZERO
      jac_next = CD_ZERO
    END IF
    CALL node_frame(r_prev, r_next, rd, u, q, chord, vr, vn, un, ok)
    IF (.NOT. ok) RETURN
    s = coef*COS(phi_x)
    force = (s*un)*vn
    IF (.NOT. PRESENT(jac_v)) RETURN
    IF (.NOT. (un > CD_ZERO)) RETURN
    ! dF/dv_r and dF/dq of F = s U v_n (inline_jacobians); dr' enters v_r with a minus sign
    CALL inline_jacobians(s, q, vr, vn, un, dfdv, dfdq)
    jac_v = dfdv
    DO j = 1, 3
      DO i = 1, 3
        proj(i, j) = -q(i)*q(j)/chord
      END DO
      proj(j, j) = proj(j, j) + CD_ONE/chord
    END DO
    jac_next = -MATMUL(dfdq, proj)
    jac_prev = -jac_next
  END SUBROUTINE CD_VIV_Node_Inline

  PURE SUBROUTINE CD_VIV_Inline_Per_Length(urel, rp, k, f, ok, jrel, jtan)
    !! In-line shedding force per unit length of a cable centreline point, f = k U v_n with
    !! k = 1/2 rho D C_vx cos(phi_x), the unit tangent t = rp/|rp|, the relative flow v_r = urel,
    !! its normal part v_n and speed U; optionally (both or neither) jrel = df/dv_r and
    !! jtan = df/drp, as CD_VIV_Lift_Per_Length. ok = .FALSE. (and f = 0) for a degenerate
    !! centreline; f and its derivatives vanish with the normal flow.
    REAL(wp), INTENT(IN) :: urel(3), rp(3), k
    REAL(wp), INTENT(OUT) :: f(3)
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp), INTENT(OUT), OPTIONAL :: jrel(3, 3), jtan(3, 3)
    REAL(wp) :: speed, t(3), vn(3), un, dfdt(3, 3), proj(3, 3)
    INTEGER :: i, j

    f = CD_ZERO
    IF (PRESENT(jrel)) THEN
      jrel = CD_ZERO
      jtan = CD_ZERO
    END IF
    speed = SQRT(rp(1)*rp(1) + rp(2)*rp(2) + rp(3)*rp(3))
    ok = speed > TINY_CHORD
    IF (.NOT. ok) RETURN
    t = rp/speed
    vn = urel - (t(1)*urel(1) + t(2)*urel(2) + t(3)*urel(3))*t
    un = SQRT(vn(1)*vn(1) + vn(2)*vn(2) + vn(3)*vn(3))
    f = (k*un)*vn
    IF (.NOT. PRESENT(jrel)) RETURN
    IF (.NOT. (un > CD_ZERO)) RETURN
    CALL inline_jacobians(k, t, urel, vn, un, jrel, dfdt)
    DO j = 1, 3
      DO i = 1, 3
        proj(i, j) = -t(i)*t(j)/speed
      END DO
      proj(j, j) = proj(j, j) + CD_ONE/speed
    END DO
    jtan = MATMUL(dfdt, proj)
  END SUBROUTINE CD_VIV_Inline_Per_Length

  PURE SUBROUTINE inline_jacobians(s, t, vr, vn, un, dfdv, dfdt)
    !! Derivatives of F = s U v_n, v_n = (I - t t^T) v_r, U = |v_n| > 0, with respect to the
    !! relative flow and the unit tangent: dF/dv_r = s [v_n (x) v_n/U + U (I - t t^T)] and, with
    !! alpha = t.v_r (dv_n/dt = -(alpha I + t (x) v_r), dU/dt = -alpha v_n/U as v_n.t = 0),
    !! dF/dt = -s [alpha v_n (x) v_n/U + U (alpha I + t (x) v_r)].
    REAL(wp), INTENT(IN) :: s, t(3), vr(3), vn(3), un
    REAL(wp), INTENT(OUT) :: dfdv(3, 3), dfdt(3, 3)
    REAL(wp) :: alpha
    INTEGER :: i, j
    alpha = t(1)*vr(1) + t(2)*vr(2) + t(3)*vr(3)
    DO j = 1, 3
      DO i = 1, 3
        dfdv(i, j) = s*(vn(i)*vn(j)/un - un*t(i)*t(j))
        dfdt(i, j) = -s*(alpha*vn(i)*vn(j)/un + un*t(i)*vr(j))
      END DO
      dfdv(j, j) = dfdv(j, j) + s*un
      dfdt(j, j) = dfdt(j, j) - s*un*alpha
    END DO
  END SUBROUTINE inline_jacobians

  PURE SUBROUTINE CD_VIV_Wake_Advance(s, om, eps, ad, f0, f1, h)
    !! Advance the wake oscillator of Facchinetti et al. (2004, Eq. 5 with the acceleration
    !! coupling F = A Y''/D; Violette et al., J. Fluids Struct. 26, 2010, Eq. 3) over a time h:
    !!
    !!   q'' + eps om (q^2 - 1) q' + om^2 q = ad y''(t)
    !!
    !! s = (q, q'), om the wake frequency 2 pi St U/D, ad = A/D and the cross-flow acceleration
    !! y'' of the line, taken linear from f0 at the start of the interval to f1 at its end.
    !! Classical Runge-Kutta in equal sub-steps of r dt <= 0.1 (at least one), r = om max(1,
    !! eps q^2) the faster of the oscillation and the van der Pol damping rate eps om (q^2 - 1) at
    !! the start of the interval, so the error is that of a resolved oscillation however long h is
    !! against the wake period and however large eps is. Should the rate at a sub-step exceed
    !! 0.2/dt (q grown on the way) or the result not be finite, the interval is redone from its
    !! start with twice the sub-steps. q = 2 C_L/C_L0: the free oscillator settles on its limit
    !! cycle q = 2 cos(om t), the lift of a fixed cylinder.
    REAL(wp), INTENT(INOUT) :: s(2)
    REAL(wp), INTENT(IN) :: om, eps, ad, f0, f1, h
    INTEGER :: n, i, attempt
    REAL(wp) :: dt, t, k1(2), k2(2), k3(2), k4(2), s0(2)
    LOGICAL :: ok
    IF (.NOT. (h > CD_ZERO)) RETURN
    s0 = s
    n = MAX(1, CEILING(MIN(rate(s0)*h/0.1_wp, 1.0e6_wp)))
    DO attempt = 1, 24
      s = s0
      dt = h/REAL(n, wp)
      ok = .TRUE.
      DO i = 0, n - 1
        IF (rate(s)*dt > 0.2_wp .AND. n < 1000000) THEN
          ok = .FALSE.
          EXIT
        END IF
        t = REAL(i, wp)*dt
        k1 = rhs(s, t)
        k2 = rhs(s + 0.5_wp*dt*k1, t + 0.5_wp*dt)
        k3 = rhs(s + 0.5_wp*dt*k2, t + 0.5_wp*dt)
        k4 = rhs(s + dt*k3, t + dt)
        s = s + (dt/6.0_wp)*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)
      END DO
      IF (ok) ok = ABS(s(1)) <= HUGE(s(1)) .AND. ABS(s(2)) <= HUGE(s(2))
      IF (ok .OR. n >= 1000000) EXIT
      n = MIN(2*n, 1000000)
    END DO
  CONTAINS
    PURE REAL(wp) FUNCTION rate(x) RESULT(r)
      !! The faster of the oscillation and the damping rate of state x (finite or not).
      REAL(wp), INTENT(IN) :: x(2)
      r = ABS(om)*MAX(CD_ONE, eps*x(1)*x(1))
      IF (.NOT. (r <= HUGE(r))) r = HUGE(r)
    END FUNCTION rate
    PURE FUNCTION rhs(x, tau) RESULT(d)
      REAL(wp), INTENT(IN) :: x(2), tau
      REAL(wp) :: d(2)
      d(1) = x(2)
      d(2) = ad*(f0 + (f1 - f0)*tau/h) - eps*om*(x(1)*x(1) - CD_ONE)*x(2) - om*om*x(1)
    END FUNCTION rhs
  END SUBROUTINE CD_VIV_Wake_Advance

  ELEMENTAL REAL(wp) FUNCTION CD_VIV_Wrap_Phase(phi) RESULT(wrapped)
    !! The phase mapped to [0, 2 pi); the model only uses cos(phi) and sin(. - phi).
    REAL(wp), INTENT(IN) :: phi
    wrapped = MODULO(phi, CD_VIV_TWO_PI)
    ! a tiny negative phase can round to exactly 2 pi
    IF (wrapped >= CD_VIV_TWO_PI .OR. wrapped < CD_ZERO) wrapped = CD_ZERO
  END FUNCTION CD_VIV_Wrap_Phase

END MODULE CableDyn_VIV
