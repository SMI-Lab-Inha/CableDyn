! File: src/CableDyn_Modal.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_Modal
  !! Modal analysis of a line about its static equilibrium: the natural frequencies and
  !! mode shapes of the generalised symmetric eigenproblem
  !!   K phi = omega^2 M phi
  !! over the free degrees of freedom (fixed or prescribed ends held), with K the static
  !! tangent stiffness and M the consistent structural mass plus the Morison added mass.
  !! K holds the element stiffness (geometric and material) and the linearised seabed
  !! normal contact; friction, damping and hydrodynamic drag are not part of an undamped
  !! modal analysis.
  !!
  !! K and M are assembled as their symmetric parts in LAPACK symmetric band storage (the
  !! nodes of a line are numbered along it, so the bandwidth is set by one element) and
  !! CD_Modal_Solve_Band computes only the requested lowest modes: the eigenvalues by
  !! LAPACK DSBGVX (band reduction and bisection), each mode shape by inverse iteration
  !! on the banded pencil, M-orthogonalised against the modes before it, and the returned
  !! eigenvalue as the Rayleigh quotient of that shape. Memory and time
  !! grow linearly with the line length for the shapes and quadratically, with a small
  !! constant, for the band reduction. The dense routines (CD_Modal_Bar_Matrices,
  !! CD_Modal_Hermite_Matrices, CD_Modal_Solve with LAPACK DSYGV) are kept as the reference
  !! the banded solver is tested against, for small lines only.
  !!
  !! EI = 0 lines (two-node bars, three translations per node): an element of rest length
  !! l0 and current length l along the unit vector e, tension T = EA (l/l0 - 1) (zero when
  !! a tension-only element is slack), has the tangent
  !!   k = (EA/l0) e e^T + (T/l) (I - e e^T)  on  [[k, -k], [-k, k]],
  !! and the consistent mass (m l0/6) [[2I, I], [I, 2I]].
  !! Finite-EI lines (cubic-Hermite, position and tangent per node): the element Hessian of
  !! CD_HermiteCable_Element and the consistent mass of CD_HermiteCable_Mass.
  USE CableDyn_Precision, ONLY: wp, CD_ZERO, CD_ONE, CD_All_Finite
  USE CableDyn_SeabedContact, ONLY: CD_Seabed_Normal_Law
  USE CableDyn_Hydro, ONLY: CD_Cable_Added_Mass_Matrix
  USE CableDyn_HermiteCable, ONLY: CD_HermiteCable_Element, CD_HermiteCable_Mass
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCable_AddedMass_Element
  USE CableDyn_Linalg, ONLY: CD_Factor_Banded, CD_Solve_Factored_Banded, CD_Blas_Runtime_Check, CD_LINALG_OK, &
                             CD_BLAS_UNAVAILABLE_INFO
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: CD_MODAL_OK = 0, CD_MODAL_BADINPUT = 1, CD_MODAL_FAIL = 2
  INTEGER, PARAMETER :: I8 = SELECTED_INT_KIND(18)
  ! Largest problem of the dense reference solver: three 4000 x 4000 matrices are 400 MB
  INTEGER, PARAMETER :: DENSE_MAX_DOF = 4000
  ! Inverse iterations per mode: at least MIN_ITER, at most MAX_ITER
  INTEGER, PARAMETER :: MIN_ITER = 3, MAX_ITER = 10
  REAL(wp), PARAMETER :: TWO_PI = 6.283185307179586476925286766559005768394_wp

  PUBLIC :: CD_Modal_Bar_Band, CD_Modal_Hermite_Band, CD_Modal_Solve_Band
  PUBLIC :: CD_Modal_Bar_Matrices, CD_Modal_Hermite_Matrices, CD_Modal_Solve
  PUBLIC :: CD_Modal_Write_Header, CD_Modal_Write_Line

CONTAINS

  SUBROUTINE fail(ErrStat, ErrMsg, msg)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    CHARACTER(*), INTENT(IN) :: msg
    ErrStat = CD_MODAL_BADINPUT
    ErrMsg = 'CableDyn_Modal: '//msg
  END SUBROUTINE fail

  SUBROUTINE add_entry(ab, kd, i, j, v)
    !! Add v, the (i, j) entry of a matrix, to the symmetric part held in upper band
    !! storage ab(kd+1+r-c, c), r <= c: a diagonal entry in full, an off-diagonal one
    !! halved (its transpose partner adds the other half).
    REAL(wp), INTENT(INOUT) :: ab(:, :)
    INTEGER, INTENT(IN) :: kd, i, j
    REAL(wp), INTENT(IN) :: v
    IF (i == j) THEN
      ab(kd + 1, j) = ab(kd + 1, j) + v
    ELSE IF (i < j) THEN
      ab(kd + 1 + i - j, j) = ab(kd + 1 + i - j, j) + 0.5_wp*v
    ELSE
      ab(kd + 1 + j - i, i) = ab(kd + 1 + j - i, i) + 0.5_wp*v
    END IF
  END SUBROUTINE add_entry

  SUBROUTINE add_block(ab, kd, gd, blk)
    !! Add every entry of an element matrix blk on the global DOFs gd.
    REAL(wp), INTENT(INOUT) :: ab(:, :)
    INTEGER, INTENT(IN) :: kd, gd(:)
    REAL(wp), INTENT(IN) :: blk(:, :)
    INTEGER :: a, b
    DO b = 1, SIZE(gd)
      DO a = 1, SIZE(gd)
        CALL add_entry(ab, kd, gd(a), gd(b), blk(a, b))
      END DO
    END DO
  END SUBROUTINE add_block

  SUBROUTINE band_to_dense(ab, kd, a)
    !! Expand a symmetric upper band to a full matrix.
    REAL(wp), INTENT(IN) :: ab(:, :)
    INTEGER, INTENT(IN) :: kd
    REAL(wp), INTENT(OUT) :: a(:, :)
    INTEGER :: i, j
    a = CD_ZERO
    DO j = 1, SIZE(ab, 2)
      DO i = MAX(1, j - kd), j
        a(i, j) = ab(kd + 1 + i - j, j)
        a(j, i) = a(i, j)
      END DO
    END DO
  END SUBROUTINE band_to_dense

  SUBROUTINE sym_band_mv(ab, kd, x, y)
    !! y = A x for a symmetric A in upper band storage.
    REAL(wp), INTENT(IN) :: ab(:, :), x(:)
    INTEGER, INTENT(IN) :: kd
    REAL(wp), INTENT(OUT) :: y(:)
    INTEGER :: i, j
    REAL(wp) :: v
    y = CD_ZERO
    DO j = 1, SIZE(x)
      y(j) = y(j) + ab(kd + 1, j)*x(j)
      DO i = MAX(1, j - kd), j - 1
        v = ab(kd + 1 + i - j, j)
        y(i) = y(i) + v*x(j)
        y(j) = y(j) + v*x(i)
      END DO
    END DO
  END SUBROUTINE sym_band_mv

  REAL(wp) FUNCTION sym_band_norm(ab, kd) RESULT(s)
    !! Infinity norm of a symmetric matrix in upper band storage.
    REAL(wp), INTENT(IN) :: ab(:, :)
    INTEGER, INTENT(IN) :: kd
    REAL(wp), ALLOCATABLE :: rows(:)
    INTEGER :: i, j
    ALLOCATE (rows(SIZE(ab, 2)))
    rows = CD_ZERO
    DO j = 1, SIZE(ab, 2)
      rows(j) = rows(j) + ABS(ab(kd + 1, j))
      DO i = MAX(1, j - kd), j - 1
        rows(i) = rows(i) + ABS(ab(kd + 1 + i - j, j))
        rows(j) = rows(j) + ABS(ab(kd + 1 + i - j, j))
      END DO
    END DO
    s = MAXVAL(rows)
  END FUNCTION sym_band_norm

  REAL(wp) FUNCTION rayleigh_quotient(a, b, kd, x) RESULT(rq)
    !! x^T A x / x^T B x for symmetric A and B in upper band storage, accumulated in
    !! quadruple precision (a product of two doubles is exact there): the low eigenvalues
    !! of a line lie many orders below its axial ones, whose terms cancel in the sum.
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL128
    REAL(wp), INTENT(IN) :: a(:, :), b(:, :), x(:)
    INTEGER, INTENT(IN) :: kd
    REAL(REAL128) :: xa, xb, p
    INTEGER :: i, j
    xa = 0.0_REAL128
    xb = 0.0_REAL128
    DO j = 1, SIZE(x)
      p = REAL(x(j), REAL128)*REAL(x(j), REAL128)
      xa = xa + REAL(a(kd + 1, j), REAL128)*p
      xb = xb + REAL(b(kd + 1, j), REAL128)*p
      DO i = MAX(1, j - kd), j - 1
        p = 2.0_REAL128*REAL(x(i), REAL128)*REAL(x(j), REAL128)
        xa = xa + REAL(a(kd + 1 + i - j, j), REAL128)*p
        xb = xb + REAL(b(kd + 1 + i - j, j), REAL128)*p
      END DO
    END DO
    rq = REAL(xa/xb, wp)
  END FUNCTION rayleigh_quotient

  SUBROUTINE CD_Modal_Bar_Band(q, elem_conn, l0, ea, rho_a, tension_only, kd, K, M, ErrStat, ErrMsg, &
                               am_rho, am_waterline, am_diam, am_can, am_cat, seabed_kn, floor_z)
    !! Stiffness and mass of an EI = 0 line (3 DOFs per node) in state q, as their symmetric
    !! parts in upper band storage K(kd+1+i-j, j), i <= j, of half-bandwidth kd (set by the
    !! element connectivity). The optional added mass (rho, nodal waterline, per-element
    !! diameter/Can/Cat) is the Morison added mass of CD_Cable_Added_Mass_Matrix; the
    !! optional seabed (nodal normal stiffness and floor elevation) adds the tangent of the
    !! C1 normal contact law.
    REAL(wp), INTENT(IN) :: q(:), l0(:), ea(:), rho_a(:)
    INTEGER, INTENT(IN) :: elem_conn(:, :)
    LOGICAL, INTENT(IN) :: tension_only
    INTEGER, INTENT(OUT) :: kd
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: K(:, :), M(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: am_rho, am_waterline(:), am_diam(:), am_can(:), am_cat(:)
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_kn(:), floor_z(:)
    INTEGER :: n, nn, ne, e, a, b, ia, ib, i, j, conn1(2, 1), gd(6)
    REAL(wp) :: d(3), l, t, ke(3, 3), mm, ke6(6, 6), me6(6, 6), ma(6, 6), q2(6), wl2(2), l02(1), f, df
    CHARACTER(200) :: em

    ErrStat = CD_MODAL_OK
    ErrMsg = ''
    kd = 0
    n = SIZE(q)
    nn = n/3
    ne = SIZE(l0)
    IF (MOD(n, 3) /= 0 .OR. n < 3 .OR. SIZE(elem_conn, 1) /= 2 .OR. SIZE(elem_conn, 2) /= ne .OR. &
        SIZE(ea) /= ne .OR. SIZE(rho_a) /= ne) THEN
      CALL fail(ErrStat, ErrMsg, 'inconsistent bar-line array sizes')
      RETURN
    END IF
    IF (.NOT. (CD_All_Finite(q) .AND. ALL(l0 > CD_ZERO) .AND. ALL(ea >= CD_ZERO) .AND. ALL(rho_a >= CD_ZERO))) THEN
      CALL fail(ErrStat, ErrMsg, 'bar-line state and properties must be finite with l0 > 0, EA and mass >= 0')
      RETURN
    END IF
    IF (ne > 0) THEN
      IF (ANY(elem_conn < 1) .OR. ANY(elem_conn > nn)) THEN
        CALL fail(ErrStat, ErrMsg, 'bar-line element connectivity out of range')
        RETURN
      END IF
    END IF
    kd = 2
    DO e = 1, ne
      kd = MAX(kd, 3*ABS(elem_conn(2, e) - elem_conn(1, e)) + 2)
    END DO
    kd = MIN(kd, n - 1)
    ALLOCATE (K(kd + 1, n), M(kd + 1, n))
    K = CD_ZERO
    M = CD_ZERO
    conn1(:, 1) = [1, 2]
    DO e = 1, ne
      a = elem_conn(1, e)
      b = elem_conn(2, e)
      ia = 3*a - 3
      ib = 3*b - 3
      gd = [ia + 1, ia + 2, ia + 3, ib + 1, ib + 2, ib + 3]
      d = q(ib + 1:ib + 3) - q(ia + 1:ia + 3)
      l = SQRT(SUM(d*d))
      IF (.NOT. (l > CD_ZERO)) THEN
        CALL fail(ErrStat, ErrMsg, 'a bar element has zero length')
        RETURN
      END IF
      d = d/l
      t = ea(e)*(l/l0(e) - CD_ONE)
      IF (tension_only .AND. t <= CD_ZERO) THEN
        ke = CD_ZERO
      ELSE
        DO j = 1, 3
          DO i = 1, 3
            ke(i, j) = (ea(e)/l0(e) - t/l)*d(i)*d(j)
          END DO
          ke(j, j) = ke(j, j) + t/l
        END DO
      END IF
      ke6(1:3, 1:3) = ke
      ke6(4:6, 4:6) = ke
      ke6(1:3, 4:6) = -ke
      ke6(4:6, 1:3) = -ke
      CALL add_block(K, kd, gd, ke6)
      mm = rho_a(e)*l0(e)/6.0_wp
      me6 = CD_ZERO
      DO i = 1, 3
        me6(i, i) = 2.0_wp*mm
        me6(i + 3, i + 3) = 2.0_wp*mm
        me6(i, i + 3) = mm
        me6(i + 3, i) = mm
      END DO
      IF (PRESENT(am_rho) .AND. PRESENT(am_waterline) .AND. PRESENT(am_diam) .AND. PRESENT(am_can) .AND. &
          PRESENT(am_cat)) THEN
        q2(1:3) = q(ia + 1:ia + 3)
        q2(4:6) = q(ib + 1:ib + 3)
        wl2 = [am_waterline(a), am_waterline(b)]
        l02(1) = l0(e)
        CALL CD_Cable_Added_Mass_Matrix(q2, conn1, l02, wl2, am_rho, am_diam(e), am_can(e), am_cat(e), ma, &
                                        ErrStat, em)
        IF (ErrStat /= 0) THEN
          CALL fail(ErrStat, ErrMsg, 'added mass: '//TRIM(em))
          RETURN
        END IF
        me6 = me6 + ma
      END IF
      CALL add_block(M, kd, gd, me6)
    END DO
    IF (PRESENT(seabed_kn) .AND. PRESENT(floor_z)) THEN
      IF (SIZE(seabed_kn) /= nn .OR. SIZE(floor_z) /= nn) THEN
        CALL fail(ErrStat, ErrMsg, 'seabed arrays must have one entry per node')
        RETURN
      END IF
      DO i = 1, nn
        CALL CD_Seabed_Normal_Law(floor_z(i) - q(3*i), seabed_kn(i), f, df)
        CALL add_entry(K, kd, 3*i, 3*i, df)
      END DO
    END IF
  END SUBROUTINE CD_Modal_Bar_Band

  SUBROUTINE CD_Modal_Hermite_Band(q, l0, ea, ei, rho_a, kd, K, M, ErrStat, ErrMsg, axial_order, bending_order, &
                                   am_rho, am_waterline, am_diam, am_can, am_cat, contact_kn, floor_z)
    !! Stiffness and mass of a cubic-Hermite line (6 DOFs per node: position and tangent)
    !! in state q, as their symmetric parts in upper band storage (half-bandwidth kd = 11):
    !! the element Hessians and consistent masses, the optional Morison added mass
    !! (CD_HermiteCable_AddedMass_Element) and the optional seabed normal contact tangent on
    !! the nodal vertical positions.
    REAL(wp), INTENT(IN) :: q(:), l0(:), ea(:), ei(:), rho_a(:)
    INTEGER, INTENT(OUT) :: kd
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: K(:, :), M(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: axial_order, bending_order
    REAL(wp), INTENT(IN), OPTIONAL :: am_rho, am_waterline, am_diam(:), am_can(:), am_cat(:)
    REAL(wp), INTENT(IN), OPTIONAL :: contact_kn(:), floor_z(:)
    INTEGER :: n, nn, ne, e, i0, i, ao, bo, gd(12)
    REAL(wp) :: qe(12), energy, fint(12), kt(12, 12), me(12, 12), ma(12, 12), f, df
    CHARACTER(200) :: em

    ErrStat = CD_MODAL_OK
    ErrMsg = ''
    kd = 0
    n = SIZE(q)
    nn = n/6
    ne = SIZE(l0)
    ao = 4
    bo = 4
    IF (PRESENT(axial_order)) ao = axial_order
    IF (PRESENT(bending_order)) bo = bending_order
    IF (MOD(n, 6) /= 0 .OR. ne < 1 .OR. ne /= nn - 1 .OR. SIZE(ea) /= ne .OR. SIZE(ei) /= ne .OR. &
        SIZE(rho_a) /= ne) THEN
      CALL fail(ErrStat, ErrMsg, 'inconsistent Hermite-line array sizes')
      RETURN
    END IF
    kd = 11
    ALLOCATE (K(kd + 1, n), M(kd + 1, n))
    K = CD_ZERO
    M = CD_ZERO
    DO e = 1, ne
      i0 = 6*(e - 1)
      DO i = 1, 12
        gd(i) = i0 + i
      END DO
      qe = q(i0 + 1:i0 + 12)
      CALL CD_HermiteCable_Element(qe, l0(e), ea(e), ei(e), energy, fint, kt, ErrStat, em, &
                                   axial_quadrature_order=ao, bending_quadrature_order=bo)
      IF (ErrStat /= 0) THEN
        CALL fail(ErrStat, ErrMsg, 'Hermite element: '//TRIM(em))
        RETURN
      END IF
      CALL add_block(K, kd, gd, kt)
      CALL CD_HermiteCable_Mass(rho_a(e), l0(e), me, ErrStat, em)
      IF (ErrStat /= 0) THEN
        CALL fail(ErrStat, ErrMsg, 'Hermite mass: '//TRIM(em))
        RETURN
      END IF
      IF (PRESENT(am_rho) .AND. PRESENT(am_waterline) .AND. PRESENT(am_diam) .AND. PRESENT(am_can) .AND. &
          PRESENT(am_cat)) THEN
        CALL CD_HermiteCable_AddedMass_Element(qe, l0(e), am_waterline, am_rho, am_diam(e), am_can(e), am_cat(e), &
                                               ma, ErrStat, em)
        IF (ErrStat /= 0) THEN
          CALL fail(ErrStat, ErrMsg, 'Hermite added mass: '//TRIM(em))
          RETURN
        END IF
        me = me + ma
      END IF
      CALL add_block(M, kd, gd, me)
    END DO
    IF (PRESENT(contact_kn) .AND. PRESENT(floor_z)) THEN
      IF (SIZE(contact_kn) /= nn .OR. SIZE(floor_z) /= nn) THEN
        CALL fail(ErrStat, ErrMsg, 'contact arrays must have one entry per node')
        RETURN
      END IF
      DO i = 1, nn
        CALL CD_Seabed_Normal_Law(floor_z(i) - q(6*i - 3), contact_kn(i), f, df)
        CALL add_entry(K, kd, 6*i - 3, 6*i - 3, df)
      END DO
    END IF
  END SUBROUTINE CD_Modal_Hermite_Band

  SUBROUTINE CD_Modal_Bar_Matrices(q, elem_conn, l0, ea, rho_a, tension_only, K, M, ErrStat, ErrMsg, &
                                   am_rho, am_waterline, am_diam, am_can, am_cat, seabed_kn, floor_z)
    !! Dense form of CD_Modal_Bar_Band (the symmetric parts), for the dense reference solver.
    REAL(wp), INTENT(IN) :: q(:), l0(:), ea(:), rho_a(:)
    INTEGER, INTENT(IN) :: elem_conn(:, :)
    LOGICAL, INTENT(IN) :: tension_only
    REAL(wp), INTENT(OUT) :: K(:, :), M(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp), INTENT(IN), OPTIONAL :: am_rho, am_waterline(:), am_diam(:), am_can(:), am_cat(:)
    REAL(wp), INTENT(IN), OPTIONAL :: seabed_kn(:), floor_z(:)
    REAL(wp), ALLOCATABLE :: kb(:, :), mb(:, :)
    INTEGER :: n, kd
    K = CD_ZERO
    M = CD_ZERO
    n = SIZE(q)
    IF (SIZE(K, 1) /= n .OR. SIZE(K, 2) /= n .OR. SIZE(M, 1) /= n .OR. SIZE(M, 2) /= n) THEN
      CALL fail(ErrStat, ErrMsg, 'inconsistent bar-line array sizes')
      RETURN
    END IF
    CALL CD_Modal_Bar_Band(q, elem_conn, l0, ea, rho_a, tension_only, kd, kb, mb, ErrStat, ErrMsg, am_rho, &
                           am_waterline, am_diam, am_can, am_cat, seabed_kn, floor_z)
    IF (ErrStat /= CD_MODAL_OK) RETURN
    CALL band_to_dense(kb, kd, K)
    CALL band_to_dense(mb, kd, M)
  END SUBROUTINE CD_Modal_Bar_Matrices

  SUBROUTINE CD_Modal_Hermite_Matrices(q, l0, ea, ei, rho_a, K, M, ErrStat, ErrMsg, axial_order, bending_order, &
                                       am_rho, am_waterline, am_diam, am_can, am_cat, contact_kn, floor_z)
    !! Dense form of CD_Modal_Hermite_Band (the symmetric parts), for the dense reference solver.
    REAL(wp), INTENT(IN) :: q(:), l0(:), ea(:), ei(:), rho_a(:)
    REAL(wp), INTENT(OUT) :: K(:, :), M(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER, INTENT(IN), OPTIONAL :: axial_order, bending_order
    REAL(wp), INTENT(IN), OPTIONAL :: am_rho, am_waterline, am_diam(:), am_can(:), am_cat(:)
    REAL(wp), INTENT(IN), OPTIONAL :: contact_kn(:), floor_z(:)
    REAL(wp), ALLOCATABLE :: kb(:, :), mb(:, :)
    INTEGER :: n, kd
    K = CD_ZERO
    M = CD_ZERO
    n = SIZE(q)
    IF (SIZE(K, 1) /= n .OR. SIZE(K, 2) /= n .OR. SIZE(M, 1) /= n .OR. SIZE(M, 2) /= n) THEN
      CALL fail(ErrStat, ErrMsg, 'inconsistent Hermite-line array sizes')
      RETURN
    END IF
    CALL CD_Modal_Hermite_Band(q, l0, ea, ei, rho_a, kd, kb, mb, ErrStat, ErrMsg, axial_order, bending_order, &
                               am_rho, am_waterline, am_diam, am_can, am_cat, contact_kn, floor_z)
    IF (ErrStat /= CD_MODAL_OK) RETURN
    CALL band_to_dense(kb, kd, K)
    CALL band_to_dense(mb, kd, M)
  END SUBROUTINE CD_Modal_Hermite_Matrices

  SUBROUTINE CD_Modal_Solve_Band(K, M, kd, free, n_modes, omega2, shapes, n_found, ErrStat, ErrMsg)
    !! The n_modes lowest eigenpairs of K phi = omega^2 M phi on the free DOFs, K and M in
    !! symmetric upper band storage of half-bandwidth kd (CD_Modal_Bar_Band,
    !! CD_Modal_Hermite_Band). The held DOFs are removed, which keeps the band. DSBGVX
    !! returns the lowest n_modes eigenvalues without eigenvectors (its eigenvectors would
    !! need an n x n work matrix); each shape is then found by inverse iteration with the
    !! banded LU factors of K - omega^2 M, M-orthogonalised against the shapes before it
    !! (so a repeated frequency, such as the two transverse planes of a taut string, gets
    !! two orthogonal shapes), and accepted on the residual of its eigen-equation. omega2
    !! is the Rayleigh quotient of the shape (quadruple-precision sums), in ascending order.
    !! shapes(:, j) is the full-length mode vector (zero on the held DOFs), M-orthonormal
    !! on the free DOFs, signed so that its largest component is positive. A negative
    !! omega2 (an unstable or compressed state) is returned as computed.
    REAL(wp), INTENT(IN) :: K(:, :), M(:, :)
    INTEGER, INTENT(IN) :: kd
    LOGICAL, INTENT(IN) :: free(:)
    INTEGER, INTENT(IN) :: n_modes
    REAL(wp), INTENT(OUT) :: omega2(:), shapes(:, :)
    INTEGER, INTENT(OUT) :: n_found, ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n, nf, kc, i, j, ii, jj, info, m_found, mode, it, ntry, prev, pass, es
    INTEGER, ALLOCATABLE :: map(:), iwork(:), ifail(:), ipiv(:)
    INTEGER(I8) :: seed
    REAL(wp), ALLOCATABLE :: ac(:, :), bc(:, :), a(:, :), b(:, :), w(:), work(:), lu(:, :), x(:), y(:), r(:)
    REAL(wp), ALLOCATABLE :: vec(:, :), bvec(:, :)
    REAL(wp) :: qdum(1, 1), zdum(1, 1), anorm, bnorm, shift, scale, nrm, res, rq, noise
    LOGICAL :: converged
    CHARACTER(512) :: em
    EXTERNAL :: dsbgvx

    ErrStat = CD_MODAL_OK
    ErrMsg = ''
    omega2 = CD_ZERO
    shapes = CD_ZERO
    n_found = 0
    n = SIZE(free)
    IF (kd < 0 .OR. SIZE(K, 1) /= kd + 1 .OR. SIZE(M, 1) /= kd + 1 .OR. SIZE(K, 2) /= n .OR. SIZE(M, 2) /= n .OR. &
        SIZE(shapes, 1) /= n .OR. SIZE(shapes, 2) < n_modes .OR. SIZE(omega2) < n_modes .OR. n_modes < 1) THEN
      CALL fail(ErrStat, ErrMsg, 'inconsistent modal array sizes')
      RETURN
    END IF
    nf = COUNT(free)
    IF (nf < 1) THEN
      CALL fail(ErrStat, ErrMsg, 'the line has no free degree of freedom')
      RETURN
    END IF
    IF (.NOT. (CD_All_Finite(K) .AND. CD_All_Finite(M))) THEN
      CALL fail(ErrStat, ErrMsg, 'the stiffness or mass matrix is not finite')
      RETURN
    END IF
    ALLOCATE (map(nf))
    j = 0
    DO i = 1, n
      IF (free(i)) THEN
        j = j + 1
        map(j) = i
      END IF
    END DO
    ! The free-DOF matrices: removing rows and columns never widens the band.
    kc = MIN(kd, nf - 1)
    ALLOCATE (ac(kc + 1, nf), bc(kc + 1, nf))
    ac = CD_ZERO
    bc = CD_ZERO
    DO jj = 1, nf
      j = map(jj)
      DO ii = MAX(1, jj - kc), jj
        i = map(ii)
        IF (j - i > kd) CYCLE
        ac(kc + 1 + ii - jj, jj) = K(kd + 1 + i - j, j)
        bc(kc + 1 + ii - jj, jj) = M(kd + 1 + i - j, j)
      END DO
    END DO
    n_found = MIN(n_modes, nf)

    ! The lowest n_found eigenvalues (DSBGVX overwrites its band copies).
    ALLOCATE (a(kc + 1, nf), b(kc + 1, nf), w(nf), work(7*nf), iwork(5*nf), ifail(nf))
    a = ac
    b = bc
    CALL dsbgvx('N', 'I', 'U', nf, kc, kc, a, kc + 1, b, kc + 1, qdum, 1, CD_ZERO, CD_ZERO, 1, n_found, &
                2.0_wp*TINY(CD_ONE), m_found, w, zdum, 1, work, iwork, ifail, info)
    IF (info /= 0 .OR. m_found /= n_found) THEN
      ErrStat = CD_MODAL_FAIL
      IF (info == CD_BLAS_UNAVAILABLE_INFO) THEN
        CALL CD_Blas_Runtime_Check('CableDyn_Modal (DSBGVX)', es, em)
        ErrMsg = 'CableDyn_Modal: '//TRIM(em)
      ELSE IF (info > nf) THEN
        ErrMsg = 'CableDyn_Modal: the mass matrix is not positive definite (DSBGVX)'
      ELSE
        BLOCK
          ! Local buffer: a record longer than the caller's ErrMsg truncates
          ! instead of aborting the internal write.
          CHARACTER(1024) :: wmsg
          INTEGER :: wios
          wmsg = ''
          WRITE (wmsg, '(A,I0,A,I0,A,I0)', IOSTAT=wios) 'CableDyn_Modal: DSBGVX failed, INFO = ', info, ', ', m_found, &
            ' of the requested eigenvalues found: ', n_found
          ErrMsg = wmsg
        END BLOCK
      END IF
      n_found = 0
      RETURN
    END IF
    DEALLOCATE (a, b, work, iwork, ifail)

    ! The shapes by inverse iteration on the banded pencil.
    anorm = sym_band_norm(ac, kc)
    bnorm = sym_band_norm(bc, kc)
    ! round-off level of an eigenvalue: eps times the largest one (about |K| / min diag M)
    noise = 1.0e3_wp*EPSILON(CD_ONE)*anorm/MINVAL(bc(kc + 1, :))
    ALLOCATE (lu(3*kc + 1, nf), ipiv(nf), x(nf), y(nf), r(nf), vec(nf, n_found), bvec(nf, n_found))
    DO mode = 1, n_found
      ! K - omega^2 M is singular to working precision by construction; an exactly zero
      ! pivot is moved off by a round-off shift of omega^2 (a unit scale when K = 0).
      shift = CD_ZERO
      scale = anorm/bnorm + ABS(w(mode))
      IF (.NOT. (scale > CD_ZERO)) scale = CD_ONE
      DO ntry = 1, 4
        CALL shifted_general_band(w(mode) + shift)
        CALL CD_Factor_Banded(lu, kc, kc, ipiv, es, em, matrix_validated=.TRUE.)
        IF (es == CD_LINALG_OK) EXIT
        shift = 10.0_wp**ntry*EPSILON(CD_ONE)*scale
      END DO
      IF (es /= CD_LINALG_OK) THEN
        CALL mode_fail('the shifted stiffness could not be factored: '//TRIM(em))
        RETURN
      END IF
      ! A deterministic pseudo-random start (a start orthogonal to the mode is improbable).
      seed = 12345_I8 + 7919_I8*INT(mode, I8)
      DO i = 1, nf
        seed = MOD(69069_I8*seed + 1_I8, 2147483647_I8)
        x(i) = REAL(seed, wp)/2147483647.0_wp - 0.5_wp
      END DO
      converged = .FALSE.
      DO it = 1, MAX_ITER
        CALL sym_band_mv(bc, kc, x, y)
        CALL CD_Solve_Factored_Banded(lu, kc, kc, ipiv, y, es, em, factor_validated=.TRUE.)
        IF (es /= CD_LINALG_OK) THEN
          CALL mode_fail('inverse iteration solve failed: '//TRIM(em))
          RETURN
        END IF
        DO pass = 1, 2
          DO prev = 1, mode - 1
            y = y - DOT_PRODUCT(bvec(:, prev), y)*vec(:, prev)
          END DO
        END DO
        CALL sym_band_mv(bc, kc, y, r)
        nrm = SQRT(MAX(DOT_PRODUCT(y, r), CD_ZERO))
        IF (.NOT. (nrm > CD_ZERO) .OR. .NOT. CD_All_Finite(y)) THEN
          CALL mode_fail('inverse iteration lost the mode')
          RETURN
        END IF
        x = y/nrm
        ! residual of K x = omega^2 M x against the scale of the two terms
        CALL sym_band_mv(ac, kc, x, r)
        CALL sym_band_mv(bc, kc, x, y)
        r = r - w(mode)*y
        res = MAXVAL(ABS(r))
        IF (it >= MIN_ITER .AND. res <= 1.0e-9_wp*(anorm + ABS(w(mode))*bnorm)*MAXVAL(ABS(x))) THEN
          converged = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. converged) THEN
        CALL mode_fail('inverse iteration did not converge')
        RETURN
      END IF
      ! The eigenvalue from the converged shape: the Rayleigh quotient is exact to second
      ! order in the shape error, where the band reduction of DSBGVX loses the low end of
      ! the spectrum to the round-off of the axial stiffness. The shape must belong to its
      ! own eigenvalue, not to another computed one distinct from it beyond round-off.
      rq = rayleigh_quotient(ac, bc, kc, x)
      DO j = 1, n_found
        IF (j == mode) CYCLE
        IF (ABS(w(j) - w(mode)) > MAX(1.0e-6_wp*MAX(ABS(w(j)), ABS(w(mode))), noise) .AND. &
            ABS(rq - w(j)) < ABS(rq - w(mode))) THEN
          CALL mode_fail('inverse iteration converged to another mode')
          RETURN
        END IF
      END DO
      i = MAXLOC(ABS(x), 1)
      IF (x(i) < CD_ZERO) x = -x
      vec(:, mode) = x
      CALL sym_band_mv(bc, kc, x, bvec(:, mode))
      omega2(mode) = rq
      shapes(map, mode) = x
    END DO
    ! ascending order (the refined values of a repeated frequency may swap)
    DO mode = 2, n_found
      DO j = mode, 2, -1
        IF (.NOT. (omega2(j) < omega2(j - 1))) EXIT
        rq = omega2(j)
        omega2(j) = omega2(j - 1)
        omega2(j - 1) = rq
        x = shapes(map, j)
        shapes(map, j) = shapes(map, j - 1)
        shapes(map, j - 1) = x
      END DO
    END DO

  CONTAINS

    SUBROUTINE shifted_general_band(sigma)
      !! lu = K - sigma M on the free DOFs in LAPACK general band storage (kl = ku = kc).
      REAL(wp), INTENT(IN) :: sigma
      INTEGER :: ic, jc
      lu = CD_ZERO
      DO jc = 1, nf
        DO ic = MAX(1, jc - kc), MIN(nf, jc + kc)
          IF (ic <= jc) THEN
            lu(2*kc + 1 + ic - jc, jc) = ac(kc + 1 + ic - jc, jc) - sigma*bc(kc + 1 + ic - jc, jc)
          ELSE
            lu(2*kc + 1 + ic - jc, jc) = ac(kc + 1 + jc - ic, ic) - sigma*bc(kc + 1 + jc - ic, ic)
          END IF
        END DO
      END DO
    END SUBROUTINE shifted_general_band

    SUBROUTINE mode_fail(msg)
      CHARACTER(*), INTENT(IN) :: msg
      ErrStat = CD_MODAL_FAIL
      BLOCK
        ! Local buffer: a record longer than the caller's ErrMsg truncates
        ! instead of aborting the internal write.
        CHARACTER(1024) :: wmsg
        INTEGER :: wios
        wmsg = ''
        WRITE (wmsg, '(A,I0,A)', IOSTAT=wios) 'CableDyn_Modal: mode ', mode, ': '//msg
        ErrMsg = wmsg
      END BLOCK
      n_found = 0
      omega2 = CD_ZERO
      shapes = CD_ZERO
    END SUBROUTINE mode_fail

  END SUBROUTINE CD_Modal_Solve_Band

  SUBROUTINE CD_Modal_Solve(K, M, free, n_modes, omega2, shapes, n_found, ErrStat, ErrMsg)
    !! Dense reference solver: the n_modes lowest eigenpairs of K phi = omega^2 M phi on the
    !! free DOFs (DSYGV, symmetric parts of K and M), for lines of at most 4000 free DOFs.
    !! shapes(:, j) is the full-length mode vector (zero on the held DOFs), M-orthonormal on
    !! the free DOFs. A negative omega2 (an unstable or compressed state) is returned as
    !! computed. Production runs use CD_Modal_Solve_Band.
    REAL(wp), INTENT(IN) :: K(:, :), M(:, :)
    LOGICAL, INTENT(IN) :: free(:)
    INTEGER, INTENT(IN) :: n_modes
    REAL(wp), INTENT(OUT) :: omega2(:), shapes(:, :)
    INTEGER, INTENT(OUT) :: n_found, ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n, nf, i, j, ii, jj, info, lwork
    INTEGER, ALLOCATABLE :: map(:)
    REAL(wp), ALLOCATABLE :: a(:, :), b(:, :), w(:), work(:)
    REAL(wp) :: wq(1)
    EXTERNAL :: dsygv

    ErrStat = CD_MODAL_OK
    ErrMsg = ''
    omega2 = CD_ZERO
    shapes = CD_ZERO
    n_found = 0
    n = SIZE(free)
    IF (SIZE(K, 1) /= n .OR. SIZE(K, 2) /= n .OR. SIZE(M, 1) /= n .OR. SIZE(M, 2) /= n .OR. &
        SIZE(shapes, 1) /= n .OR. SIZE(shapes, 2) < n_modes .OR. SIZE(omega2) < n_modes .OR. n_modes < 1) THEN
      CALL fail(ErrStat, ErrMsg, 'inconsistent modal array sizes')
      RETURN
    END IF
    nf = COUNT(free)
    IF (nf < 1) THEN
      CALL fail(ErrStat, ErrMsg, 'the line has no free degree of freedom')
      RETURN
    END IF
    IF (nf > DENSE_MAX_DOF) THEN
      CALL fail(ErrStat, ErrMsg, 'the line has more free DOFs than the dense reference solver takes (4000)')
      RETURN
    END IF
    IF (.NOT. (CD_All_Finite(K) .AND. CD_All_Finite(M))) THEN
      CALL fail(ErrStat, ErrMsg, 'the stiffness or mass matrix is not finite')
      RETURN
    END IF
    ALLOCATE (map(nf), a(nf, nf), b(nf, nf), w(nf))
    j = 0
    DO i = 1, n
      IF (free(i)) THEN
        j = j + 1
        map(j) = i
      END IF
    END DO
    DO jj = 1, nf
      DO ii = 1, nf
        a(ii, jj) = 0.5_wp*(K(map(ii), map(jj)) + K(map(jj), map(ii)))
        b(ii, jj) = 0.5_wp*(M(map(ii), map(jj)) + M(map(jj), map(ii)))
      END DO
    END DO
    CALL dsygv(1, 'V', 'U', nf, a, nf, b, nf, w, wq, -1, info)
    lwork = MAX(1, INT(wq(1)), 3*nf)
    ALLOCATE (work(lwork))
    CALL dsygv(1, 'V', 'U', nf, a, nf, b, nf, w, work, lwork, info)
    IF (info /= 0) THEN
      ErrStat = CD_MODAL_FAIL
      IF (info > nf) THEN
        ErrMsg = 'CableDyn_Modal: the mass matrix is not positive definite (DSYGV)'
      ELSE
        BLOCK
          ! Local buffer: a record longer than the caller's ErrMsg truncates
          ! instead of aborting the internal write.
          CHARACTER(1024) :: wmsg
          INTEGER :: wios
          wmsg = ''
          WRITE (wmsg, '(A,I0)', IOSTAT=wios) 'CableDyn_Modal: DSYGV failed, INFO = ', info
          ErrMsg = wmsg
        END BLOCK
      END IF
      RETURN
    END IF
    n_found = MIN(n_modes, nf)
    DO j = 1, n_found
      omega2(j) = w(j)
      shapes(map, j) = a(:, j)
    END DO
  END SUBROUTINE CD_Modal_Solve

  SUBROUTINE CD_Modal_Write_Header(unit)
    !! The two-table header of <root>.modes.out; the frequency rows of every line follow,
    !! then CD_Modal_Write_Line appends each line's shape rows under a second header.
    INTEGER, INTENT(IN) :: unit
    WRITE (unit, '(A)') '# CableDyn modal analysis about the static equilibrium: K phi = omega^2 M phi '// &
      '(static tangent, structural + added mass, fixed and prescribed ends held)'
  END SUBROUTINE CD_Modal_Write_Header

  SUBROUTINE CD_Modal_Write_Line(unit, line_id, q, dof_per_node, omega2, shapes, n_found, part)
    !! part = 1: the frequency rows "LineID Mode Frequency Period Omega"; part = 2: the
    !! shape rows "LineID Mode Node X Y Z dX dY dZ", the nodal translations of each mode
    !! scaled to a unit largest nodal displacement (the static position X Y Z beside it).
    !! A negative omega^2 is written with a negative frequency and period 0.
    INTEGER, INTENT(IN) :: unit, line_id, dof_per_node, n_found, part
    REAL(wp), INTENT(IN) :: q(:), omega2(:), shapes(:, :)
    INTEGER :: j, i, nn, i0
    REAL(wp) :: om, fr, per, scl, dmax
    CHARACTER(1), PARAMETER :: TB = ACHAR(9)
    nn = SIZE(q)/dof_per_node
    DO j = 1, n_found
      IF (part == 1) THEN
        om = SIGN(SQRT(ABS(omega2(j))), omega2(j))
        fr = om/TWO_PI
        per = CD_ZERO
        IF (om > CD_ZERO) per = CD_ONE/fr
        WRITE (unit, '(I0,A,I0,3(A,ES16.8))') line_id, TB, j, TB, fr, TB, per, TB, om
      ELSE
        dmax = CD_ZERO
        DO i = 1, nn
          i0 = dof_per_node*(i - 1)
          dmax = MAX(dmax, SQRT(SUM(shapes(i0 + 1:i0 + 3, j)**2)))
        END DO
        scl = CD_ONE
        IF (dmax > CD_ZERO) scl = CD_ONE/dmax
        DO i = 1, nn
          i0 = dof_per_node*(i - 1)
          WRITE (unit, '(I0,A,I0,A,I0,6(A,ES16.8))') line_id, TB, j, TB, i, TB, q(i0 + 1), TB, q(i0 + 2), TB, &
            q(i0 + 3), TB, scl*shapes(i0 + 1, j), TB, scl*shapes(i0 + 2, j), TB, scl*shapes(i0 + 3, j)
        END DO
      END IF
    END DO
  END SUBROUTINE CD_Modal_Write_Line

END MODULE CableDyn_Modal
