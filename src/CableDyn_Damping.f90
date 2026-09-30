! File: src/CableDyn_Damping.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_Damping
  !! MoorDyn-convention axial (BA) damping for the positions-only EI=0 cable path. A
  !! velocity-dependent
  !! internal load that opposes the axial stretch rate of each element, with analytic
  !! position and velocity Jacobians so it enters the implicit generalised-alpha
  !! tangent exactly.
  !!
  !! Per element e with end nodes a, b (current positions r_a, r_b; velocities v_a, v_b):
  !!   chord = r_b - r_a ;  length = ||chord|| ;  t = chord / length  (unit tangent)
  !!   rel   = v_b - v_a ;  ldot = t . rel      (axial stretch rate, m/s)
  !!   k     = BA(e) / L0(e)                     (axial damping coefficient, N.s/m)
  !!   Td    = k * ldot                          (signed damping tension, N)
  !! The nodal force is +Td t on node a and -Td t on node b (internal stress
  !! MoorDyn line convention). It vanishes at rest (ldot = 0 when v = 0), is
  !! dissipative (power f.v <= 0).
  !!
  !! BA resolution (CD_Resolve_Legacy_BA): the BA column is overloaded, as in MoorDyn
  !! -- a non-negative value IS the coefficient (N.s); a NEGATIVE
  !! value is the desired damping ratio -zeta, giving BA = zeta * L_seg * sqrt(EA * w)
  !! (w = mass per unstretched metre), so each segment's axial mode carries that
  !! fraction of critical damping.
  !!
  !! Load convention (matches CD_Cable_Dynamic_Load_Proc / CableDyn_Hydro): the
  !! returned `force` is the external load to SUBTRACT from the residual, and
  !! `jac_q`, `jac_v` are the residual-Jacobian contributions -d(force)/dq and
  !! -d(force)/dv. They are analytic (closed-form, checked against finite differences in
  !! the tests): simple per-element 3x3 blocks.
  USE CableDyn_Precision, ONLY: wp, CD_ZERO, CD_ONE, CD_All_Finite, CD_Is_Finite
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CD_Resolve_Legacy_BA
  PUBLIC :: CD_Cable_Axial_Damping_Force
  PUBLIC :: CD_Cable_Axial_Damping_Load
  PUBLIC :: CD_Axial_Damping_Element_Force
  PUBLIC :: CD_Axial_Damping_Element_Load
  PUBLIC :: CD_Cable_Element_Damping_Tension

  INTEGER, PARAMETER, PUBLIC :: CD_DAMP_OK = 0
  INTEGER, PARAMETER, PUBLIC :: CD_DAMP_BADINPUT = 1

  REAL(wp), PARAMETER :: TINY_LEN = 1.0e-12_wp   ! collapsed-element guard

CONTAINS

  SUBROUTINE CD_Resolve_Legacy_BA(ba_input, segment_length, ea, mass_per_length, ba, ErrStat, ErrMsg)
    !! Resolve a BA/-zeta line-type value to an axial damping coefficient [N.s].
    !! ba_input >= 0 is the coefficient directly (EA/mass ignored); ba_input < 0 is the
    !! -zeta damping ratio, resolved as BA = (-ba_input) * segment_length * sqrt(EA * w).
    REAL(wp), INTENT(IN)  :: ba_input          !! BA value (negative => -zeta ratio)
    REAL(wp), INTENT(IN)  :: segment_length    !! unstretched segment length L_seg [m] (> 0 if -zeta)
    REAL(wp), INTENT(IN)  :: ea                !! axial stiffness EA [N] (> 0 if -zeta)
    REAL(wp), INTENT(IN)  :: mass_per_length   !! dry mass per unstretched metre w [kg/m] (> 0 if -zeta)
    REAL(wp), INTENT(OUT) :: ba                !! resolved damping coefficient [N.s]
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    ba = CD_ZERO
    ErrStat = CD_DAMP_OK
    ErrMsg = ''
    IF (.NOT. CD_Is_Finite(ba_input)) THEN
      CALL fail(ErrStat, ErrMsg, 'BA must be finite'); RETURN
    END IF
    IF (ba_input >= CD_ZERO) THEN
      ba = ba_input                      ! direct coefficient
      RETURN
    END IF
    IF (.NOT. (CD_Is_Finite(segment_length) .AND. segment_length > CD_ZERO) .OR. &
        .NOT. (CD_Is_Finite(ea) .AND. ea > CD_ZERO) .OR. &
        .NOT. (CD_Is_Finite(mass_per_length) .AND. mass_per_length > CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'resolving a -zeta damping ratio needs positive segment_length, EA, '// &
                'and mass_per_length')
      RETURN
    END IF
    ba = (-ba_input)*segment_length*SQRT(ea*mass_per_length)
  END SUBROUTINE CD_Resolve_Legacy_BA

  SUBROUTINE CD_Cable_Element_Damping_Tension(q, v, elem_conn, l0, ba, td, ErrStat, ErrMsg, l0_dot)
    !! Signed per-element axial damping tension Td = (BA/L0)(t . (v_b - v_a)).
    !! Positive when the element is stretching, negative when shortening; zero at rest.
    !! The reported physical tension is the elastic T plus this Td.
    REAL(wp), INTENT(IN)  :: q(:), v(:)        !! (3 n_nodes) positions / velocities
    INTEGER, INTENT(IN)  :: elem_conn(:, :)   !! (2, n_elem), 1-based
    REAL(wp), INTENT(IN)  :: l0(:)             !! (n_elem) unstretched lengths (> 0)
    REAL(wp), INTENT(IN)  :: ba(:)             !! (n_elem) per-element damping coefficient [N.s]
    REAL(wp), INTENT(OUT) :: td(:)             !! (n_elem) signed damping tension [N]
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n_dof, n_elem, e, ia, ib
    ! l0_dot (optional, per element): the unstretched-length RATE (active line
    ! control, MoorDyn's ld): the damped rate is L0 times the strain rate,
    ! L0*d(length/L0)/dt = ldot - (length/L0)*l0_dot (MoorDyn-F), so a segment paid
    ! out at constant strain carries no damping tension.
    REAL(wp), INTENT(IN), OPTIONAL :: l0_dot(:)
    REAL(wp) :: chord(3), tang(3), rel(3), length, ldot

    td = CD_ZERO
    CALL validate(q, v, elem_conn, l0, ba, ErrStat, ErrMsg)
    IF (ErrStat /= CD_DAMP_OK) RETURN
    IF (PRESENT(l0_dot)) THEN
      IF (SIZE(l0_dot) /= SIZE(l0) .OR. .NOT. CD_All_Finite(l0_dot)) THEN
        CALL fail(ErrStat, ErrMsg, 'l0_dot must be finite with one entry per element')
        RETURN
      END IF
    END IF
    n_dof = SIZE(q)
    n_elem = SIZE(elem_conn, 2)
    IF (SIZE(td) /= n_elem) THEN
      CALL fail(ErrStat, ErrMsg, 'td must have shape (n_elem)'); RETURN
    END IF
    DO e = 1, n_elem
      ia = elem_conn(1, e); ib = elem_conn(2, e)
      chord = q(3*ib - 2:3*ib) - q(3*ia - 2:3*ia)
      length = NORM2(chord)
      IF (length <= TINY_LEN) THEN
        CALL fail(ErrStat, ErrMsg, 'axial-damping element collapsed to zero length'); RETURN
      END IF
      tang = chord/length
      rel = v(3*ib - 2:3*ib) - v(3*ia - 2:3*ia)
      ldot = DOT_PRODUCT(tang, rel)
      IF (PRESENT(l0_dot)) ldot = ldot - (length/l0(e))*l0_dot(e)
      td(e) = ba(e)/l0(e)*ldot
    END DO
  END SUBROUTINE CD_Cable_Element_Damping_Tension

  SUBROUTINE CD_Cable_Axial_Damping_Force(q, v, elem_conn, l0, ba, force, ErrStat, ErrMsg, l0_dot)
    !! Assemble the global axial-damping nodal force (force-only; the residual subtracts
    !! it). +Td t on node a, -Td t on node b for every element.
    REAL(wp), INTENT(IN)  :: q(:), v(:)
    INTEGER, INTENT(IN)  :: elem_conn(:, :)
    REAL(wp), INTENT(IN)  :: l0(:), ba(:)
    REAL(wp), INTENT(OUT) :: force(:)          !! (3 n_nodes)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n_dof, n_elem, e, ia, ib
    REAL(wp), INTENT(IN), OPTIONAL :: l0_dot(:)   !! per-element unstretched-length rate (see the tension routine)
    REAL(wp) :: chord(3), tang(3), rel(3), length, ldot, td_vec(3)

    force = CD_ZERO
    CALL validate(q, v, elem_conn, l0, ba, ErrStat, ErrMsg)
    IF (ErrStat /= CD_DAMP_OK) RETURN
    IF (PRESENT(l0_dot)) THEN
      IF (SIZE(l0_dot) /= SIZE(l0) .OR. .NOT. CD_All_Finite(l0_dot)) THEN
        CALL fail(ErrStat, ErrMsg, 'l0_dot must be finite with one entry per element')
        RETURN
      END IF
    END IF
    n_dof = SIZE(q)
    IF (SIZE(force) /= n_dof) THEN
      CALL fail(ErrStat, ErrMsg, 'force must match the state shape'); RETURN
    END IF
    n_elem = SIZE(elem_conn, 2)
    DO e = 1, n_elem
      ia = elem_conn(1, e); ib = elem_conn(2, e)
      chord = q(3*ib - 2:3*ib) - q(3*ia - 2:3*ia)
      length = NORM2(chord)
      IF (length <= TINY_LEN) THEN
        CALL fail(ErrStat, ErrMsg, 'axial-damping element collapsed to zero length'); RETURN
      END IF
      tang = chord/length
      rel = v(3*ib - 2:3*ib) - v(3*ia - 2:3*ia)
      ldot = DOT_PRODUCT(tang, rel)
      IF (PRESENT(l0_dot)) ldot = ldot - (length/l0(e))*l0_dot(e)
      td_vec = (ba(e)/l0(e)*ldot)*tang
      force(3*ia - 2:3*ia) = force(3*ia - 2:3*ia) + td_vec
      force(3*ib - 2:3*ib) = force(3*ib - 2:3*ib) - td_vec
    END DO
  END SUBROUTINE CD_Cable_Axial_Damping_Force

  SUBROUTINE CD_Cable_Axial_Damping_Load(q, v, elem_conn, l0, ba, force, jac_q, jac_v, ErrStat, ErrMsg, l0_dot)
    !! Force + analytic residual Jacobians for the axial damping. `force` is the load
    !! to subtract from the residual; `jac_q = -d(force)/dq` and `jac_v = -d(force)/dv`
    !! (the CD_Cable_Dynamic_Load_Proc convention). Per element, with k = BA/L0:
    !!   d(force_a)/dv = [ -k t t^T (a) , +k t t^T (b) ] ; force_b = -force_a
    !!   d(force_a)/dq = [ -D (a) , +D (b) ] with D = k (t rel^T + ldot I)(I - t t^T)/length
    !! and the returned Jacobians are the NEGATED block patterns of those. With l0_dot,
    !! ldot is the damped rate ldot - (length/L0)*l0_dot and D gains -k (l0_dot/L0) t t^T.
    REAL(wp), INTENT(IN)  :: q(:), v(:)
    INTEGER, INTENT(IN)  :: elem_conn(:, :)
    REAL(wp), INTENT(IN)  :: l0(:), ba(:)
    REAL(wp), INTENT(OUT) :: force(:)          !! (3 n_nodes)
    REAL(wp), INTENT(OUT) :: jac_q(:, :), jac_v(:, :)   !! (3 n_nodes, 3 n_nodes)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n_dof, n_elem, e, ia, ib, i
    REAL(wp), INTENT(IN), OPTIONAL :: l0_dot(:)   !! per-element unstretched-length rate (see the tension routine)
    REAL(wp) :: chord(3), tang(3), rel(3), length, ldot, k, td_vec(3)
    REAL(wp) :: eye(3, 3), ktt(3, 3), m1(3, 3), m2(3, 3), dmat(3, 3)

    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    CALL validate(q, v, elem_conn, l0, ba, ErrStat, ErrMsg)
    IF (ErrStat /= CD_DAMP_OK) RETURN
    IF (PRESENT(l0_dot)) THEN
      IF (SIZE(l0_dot) /= SIZE(l0) .OR. .NOT. CD_All_Finite(l0_dot)) THEN
        CALL fail(ErrStat, ErrMsg, 'l0_dot must be finite with one entry per element')
        RETURN
      END IF
    END IF
    n_dof = SIZE(q)
    IF (SIZE(force) /= n_dof .OR. SIZE(jac_q, 1) /= n_dof .OR. SIZE(jac_q, 2) /= n_dof .OR. &
        SIZE(jac_v, 1) /= n_dof .OR. SIZE(jac_v, 2) /= n_dof) THEN
      CALL fail(ErrStat, ErrMsg, 'force/Jacobian shapes must match a positions-only state'); RETURN
    END IF
    n_elem = SIZE(elem_conn, 2)
    eye = CD_ZERO
    DO i = 1, 3
      eye(i, i) = CD_ONE
    END DO

    DO e = 1, n_elem
      ia = elem_conn(1, e); ib = elem_conn(2, e)
      chord = q(3*ib - 2:3*ib) - q(3*ia - 2:3*ia)
      length = NORM2(chord)
      IF (length <= TINY_LEN) THEN
        CALL fail(ErrStat, ErrMsg, 'axial-damping element collapsed to zero length'); RETURN
      END IF
      tang = chord/length
      rel = v(3*ib - 2:3*ib) - v(3*ia - 2:3*ia)
      ldot = DOT_PRODUCT(tang, rel)
      IF (PRESENT(l0_dot)) ldot = ldot - (length/l0(e))*l0_dot(e)
      k = ba(e)/l0(e)

      ! force
      td_vec = (k*ldot)*tang
      force(3*ia - 2:3*ia) = force(3*ia - 2:3*ia) + td_vec
      force(3*ib - 2:3*ib) = force(3*ib - 2:3*ib) - td_vec

      ! velocity Jacobian block ktt = k t t^T ; d(force_a)/dv = [-ktt(a) +ktt(b)]
      ktt = k*outer(tang, tang)
      ! returned jac_v = -d(force)/dv:  (a,a)=+ktt (a,b)=-ktt (b,a)=-ktt (b,b)=+ktt
      CALL add_block(jac_v, ia, ia, ktt)
      CALL add_block(jac_v, ia, ib, -ktt)
      CALL add_block(jac_v, ib, ia, -ktt)
      CALL add_block(jac_v, ib, ib, ktt)

      ! position Jacobian D = k (t rel^T + ldot I)(I - t t^T)/length ; d(force_a)/dchord = D
      m1 = k*(outer(tang, rel) + ldot*eye)
      m2 = (eye - outer(tang, tang))/length
      dmat = MATMUL(m1, m2)
      ! the paid-out term -(length/L0)*l0_dot stretches with the chord: d/dchord adds -k*(l0_dot/L0) t t^T
      IF (PRESENT(l0_dot)) dmat = dmat - (k*l0_dot(e)/l0(e))*outer(tang, tang)
      ! d(force)/dq blocks: (a,a)=-D (a,b)=+D (b,a)=+D (b,b)=-D ; returned jac_q = -that
      CALL add_block(jac_q, ia, ia, dmat)
      CALL add_block(jac_q, ia, ib, -dmat)
      CALL add_block(jac_q, ib, ia, -dmat)
      CALL add_block(jac_q, ib, ib, dmat)
    END DO
  END SUBROUTINE CD_Cable_Axial_Damping_Load

  SUBROUTINE CD_Axial_Damping_Element_Force(qa, qb, va, vb, l0, ba, force, ErrStat, ErrMsg, l0_dot)
    !! Two-node axial damping force kernel for model hot paths.
    REAL(wp), INTENT(IN) :: qa(3), qb(3), va(3), vb(3), l0, ba
    REAL(wp), INTENT(OUT) :: force(6)
    REAL(wp), INTENT(IN), OPTIONAL :: l0_dot   !! unstretched-length rate (see the cable-path routines)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: chord(3), tang(3), rel(3), length, ldot, td_vec(3)

    CALL validate_element_inputs(qa, qb, va, vb, l0, ba, ErrStat, ErrMsg)
    force = CD_ZERO
    IF (ErrStat /= CD_DAMP_OK) RETURN
    IF (PRESENT(l0_dot)) THEN
      IF (.NOT. CD_Is_Finite(l0_dot)) THEN
        CALL fail(ErrStat, ErrMsg, 'l0_dot must be finite')
        RETURN
      END IF
    END IF
    chord = qb - qa
    length = NORM2(chord)
    IF (length <= TINY_LEN) THEN
      CALL fail(ErrStat, ErrMsg, 'axial-damping element collapsed to zero length')
      RETURN
    END IF
    tang = chord/length
    rel = vb - va
    ldot = DOT_PRODUCT(tang, rel)
    IF (PRESENT(l0_dot)) ldot = ldot - (length/l0)*l0_dot
    td_vec = (ba/l0*ldot)*tang
    force(1:3) = td_vec
    force(4:6) = -td_vec
  END SUBROUTINE CD_Axial_Damping_Element_Force

  SUBROUTINE CD_Axial_Damping_Element_Load(qa, qb, va, vb, l0, ba, force, jac_q, jac_v, ErrStat, ErrMsg, l0_dot)
    !! Two-node axial damping force and residual-Jacobian kernel.
    REAL(wp), INTENT(IN) :: qa(3), qb(3), va(3), vb(3), l0, ba
    REAL(wp), INTENT(OUT) :: force(6), jac_q(6, 6), jac_v(6, 6)
    REAL(wp), INTENT(IN), OPTIONAL :: l0_dot   !! unstretched-length rate: the damped rate is ldot - (length/l0)*l0_dot
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: i
    REAL(wp) :: chord(3), tang(3), rel(3), length, ldot, k, td_vec(3)
    REAL(wp) :: eye(3, 3), ktt(3, 3), m1(3, 3), m2(3, 3), dmat(3, 3)

    CALL validate_element_inputs(qa, qb, va, vb, l0, ba, ErrStat, ErrMsg)
    force = CD_ZERO
    jac_q = CD_ZERO
    jac_v = CD_ZERO
    IF (ErrStat /= CD_DAMP_OK) RETURN
    IF (PRESENT(l0_dot)) THEN
      IF (.NOT. CD_Is_Finite(l0_dot)) THEN
        CALL fail(ErrStat, ErrMsg, 'l0_dot must be finite')
        RETURN
      END IF
    END IF
    chord = qb - qa
    length = NORM2(chord)
    IF (length <= TINY_LEN) THEN
      CALL fail(ErrStat, ErrMsg, 'axial-damping element collapsed to zero length')
      RETURN
    END IF
    tang = chord/length
    rel = vb - va
    ldot = DOT_PRODUCT(tang, rel)
    IF (PRESENT(l0_dot)) ldot = ldot - (length/l0)*l0_dot
    k = ba/l0

    td_vec = (k*ldot)*tang
    force(1:3) = td_vec
    force(4:6) = -td_vec

    eye = CD_ZERO
    DO i = 1, 3
      eye(i, i) = CD_ONE
    END DO
    ktt = k*outer(tang, tang)
    jac_v(1:3, 1:3) = jac_v(1:3, 1:3) + ktt
    jac_v(1:3, 4:6) = jac_v(1:3, 4:6) - ktt
    jac_v(4:6, 1:3) = jac_v(4:6, 1:3) - ktt
    jac_v(4:6, 4:6) = jac_v(4:6, 4:6) + ktt

    m1 = k*(outer(tang, rel) + ldot*eye)
    m2 = (eye - outer(tang, tang))/length
    dmat = MATMUL(m1, m2)
    ! the paid-out term -(length/L0)*l0_dot stretches with the chord: d/dchord adds -k*(l0_dot/L0) t t^T
    IF (PRESENT(l0_dot)) dmat = dmat - (k*l0_dot/l0)*outer(tang, tang)
    jac_q(1:3, 1:3) = jac_q(1:3, 1:3) + dmat
    jac_q(1:3, 4:6) = jac_q(1:3, 4:6) - dmat
    jac_q(4:6, 1:3) = jac_q(4:6, 1:3) - dmat
    jac_q(4:6, 4:6) = jac_q(4:6, 4:6) + dmat
  END SUBROUTINE CD_Axial_Damping_Element_Load

  ! --------------------------------------------------------------------------- !
  ! private helpers                                                             !
  ! --------------------------------------------------------------------------- !

  PURE FUNCTION outer(a, b) RESULT(m)
    !! 3x3 outer product a b^T.
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: m(3, 3)
    INTEGER :: i, j
    DO j = 1, 3
      DO i = 1, 3
        m(i, j) = a(i)*b(j)
      END DO
    END DO
  END FUNCTION outer

  SUBROUTINE add_block(mat, node_i, node_j, blk)
    !! Accumulate a 3x3 block into the (node_i, node_j) DOF positions of a global matrix.
    REAL(wp), INTENT(INOUT) :: mat(:, :)
    INTEGER, INTENT(IN)  :: node_i, node_j
    REAL(wp), INTENT(IN)  :: blk(3, 3)
    INTEGER :: ri, ci
    ri = 3*node_i - 3
    ci = 3*node_j - 3
    mat(ri + 1:ri + 3, ci + 1:ci + 3) = mat(ri + 1:ri + 3, ci + 1:ci + 3) + blk
  END SUBROUTINE add_block

  SUBROUTINE validate(q, v, elem_conn, l0, ba, ErrStat, ErrMsg)
    !! Shared input validation (fail closed) for the damping routines.
    REAL(wp), INTENT(IN)  :: q(:), v(:), l0(:), ba(:)
    INTEGER, INTENT(IN)  :: elem_conn(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n_dof, n_nodes, n_elem
    ErrStat = CD_DAMP_OK
    ErrMsg = ''
    n_dof = SIZE(q)
    IF (MOD(n_dof, 3) /= 0 .OR. n_dof < 6 .OR. SIZE(v) /= n_dof) THEN
      CALL fail(ErrStat, ErrMsg, 'q and v must be a positions-only state (3 n_nodes), n_nodes >= 2'); RETURN
    END IF
    n_nodes = n_dof/3
    n_elem = SIZE(elem_conn, 2)
    IF (SIZE(elem_conn, 1) /= 2 .OR. n_elem < 1) THEN
      CALL fail(ErrStat, ErrMsg, 'elem_conn must be (2, n_elem >= 1)'); RETURN
    END IF
    IF (SIZE(l0) /= n_elem .OR. SIZE(ba) /= n_elem) THEN
      CALL fail(ErrStat, ErrMsg, 'l0 and ba must have shape (n_elem)'); RETURN
    END IF
    IF (MINVAL(elem_conn) < 1 .OR. MAXVAL(elem_conn) > n_nodes) THEN
      CALL fail(ErrStat, ErrMsg, 'elem_conn node indices out of range'); RETURN
    END IF
    IF (.NOT. CD_All_Finite(q) .OR. .NOT. CD_All_Finite(v) .OR. &
        .NOT. CD_All_Finite(l0) .OR. ANY(l0 <= CD_ZERO) .OR. .NOT. CD_All_Finite(ba)) THEN
      CALL fail(ErrStat, ErrMsg, 'q, v, ba must be finite and l0 finite-positive'); RETURN
    END IF
    ! ba is the RESOLVED damping coefficient (CD_Resolve_Legacy_BA returns >= 0). A
    ! negative value -- e.g. an unresolved -zeta input passed by mistake --
    ! would flip the force into anti-damping (energy injection), so reject it.
    IF (ANY(ba < CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'ba must be the resolved coefficient (>= 0); pass it through '// &
                'CD_Resolve_Legacy_BA first'); RETURN
    END IF
  END SUBROUTINE validate

  SUBROUTINE validate_element_inputs(qa, qb, va, vb, l0, ba, ErrStat, ErrMsg)
    REAL(wp), INTENT(IN) :: qa(3), qb(3), va(3), vb(3), l0, ba
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    ErrStat = CD_DAMP_OK
    ErrMsg = ''
    IF (.NOT. (CD_All_Finite(qa) .AND. CD_All_Finite(qb) .AND. &
               CD_All_Finite(va) .AND. CD_All_Finite(vb))) THEN
      CALL fail(ErrStat, ErrMsg, 'axial-damping element state must be finite')
      RETURN
    END IF
    IF (.NOT. (CD_Is_Finite(l0) .AND. l0 > CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'axial-damping element l0 must be finite-positive')
      RETURN
    END IF
    IF (.NOT. (CD_Is_Finite(ba) .AND. ba >= CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'axial-damping element ba must be the resolved coefficient (>= 0)')
    END IF
  END SUBROUTINE validate_element_inputs

  SUBROUTINE fail(ErrStat, ErrMsg, msg)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    CHARACTER(*), INTENT(IN)  :: msg
    ErrStat = CD_DAMP_BADINPUT
    ErrMsg = 'CableDyn_Damping: '//msg
  END SUBROUTINE fail

END MODULE CableDyn_Damping
