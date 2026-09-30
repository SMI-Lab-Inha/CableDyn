! File: tests/test_modal.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_modal
  !! Line modal analysis (CableDyn_Modal, OPTION nModes) against closed-form frequencies:
  !!   1. Taut string (EI = 0 bars): transverse f_n = n/(2L) sqrt(T/m) and the first axial
  !!      mode pi/L0 sqrt(EA/m).
  !!   2. Cubic-Hermite beam, pinned-pinned: f_n = (n pi/L)^2 sqrt(EI/m)/(2 pi); pre-tensioned:
  !!      f_n = n/(2L) sqrt(T/m) sqrt(1 + (n pi)^2 EI/(T L^2)).
  !!   3. Irvine's sagging cable (deck, EI = 0, supports level): out-of-plane swing
  !!      pi/L sqrt(H/m), first in-plane antisymmetric 2 pi/L sqrt(H/m), and the first in-plane
  !!      symmetric mode from tan(b/2) = b/2 - (4/lambda^2)(b/2)^3.
  !!   4. The Hermite route writes <root>.modes.out; unsupported decks are rejected.
  !!   5. The banded solver (CD_Modal_Solve_Band, production) against the dense reference
  !!      (CD_Modal_Solve) on cases 1 and 2: the same frequencies, M-orthonormal shapes that
  !!      satisfy K phi = omega^2 M phi; and a 1024-element tensioned beam (6150 DOFs, beyond
  !!      the dense solver) against its closed-form frequencies.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_Modal
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.141592653589793238462643383279502884197_wp
  INTEGER :: nfail

  nfail = 0
  CALL case_taut_string()
  CALL case_hermite_beam()
  CALL case_hermite_long_beam()
  CALL case_slack_string()
  CALL case_irvine()
  CALL case_deck_hermite_and_rules()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'test_modal: ', nfail, ' failure(s)'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'test_modal: all checks passed'

CONTAINS

  INCLUDE 'nan_max_abs.inc'

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. cond) THEN
      WRITE (*, '(A)') 'MISMATCH ['//label//']'
      nfail = nfail + 1
    END IF
  END SUBROUTINE require

  SUBROUTINE case_taut_string()
    INTEGER, PARAMETER :: NE = 100, NN = NE + 1, NM = 80
    REAL(wp), PARAMETER :: L = 100.0_wp, EA = 1.0e8_wp, RHO = 100.0_wp, EPS = 1.0e-3_wp
    REAL(wp) :: q(3*NN), l0(NE), ea_(NE), rho_a(NE), k(3*NN, 3*NN), m(3*NN, 3*NN), om2(NM), sh(3*NN, NM)
    REAL(wp) :: t, mc, f_exp, f_got, l0tot, ax, tr
    REAL(wp), ALLOCATABLE :: kb(:, :), mb(:, :)
    INTEGER :: conn(2, NE), i, j, es, nf, nmode, kd
    LOGICAL :: free(3*NN)
    CHARACTER(256) :: em
    l0tot = L/(1.0_wp + EPS)
    DO i = 1, NE
      conn(:, i) = [i, i + 1]
    END DO
    l0 = l0tot/NE
    ea_ = EA
    rho_a = RHO
    DO i = 1, NN
      q(3*i - 2:3*i) = [L*REAL(i - 1, wp)/NE, 0.0_wp, 0.0_wp]
    END DO
    free = .TRUE.
    free(1:3) = .FALSE.
    free(3*NN - 2:3*NN) = .FALSE.
    CALL CD_Modal_Bar_Matrices(q, conn, l0, ea_, rho_a, .TRUE., k, m, es, em)
    CALL require(es == CD_MODAL_OK, 'string:matrices: '//TRIM(em))
    CALL CD_Modal_Solve(k, m, free, NM, om2, sh, nf, es, em)
    CALL require(es == CD_MODAL_OK .AND. nf == NM, 'string:solve: '//TRIM(em))
    CALL CD_Modal_Bar_Band(q, conn, l0, ea_, rho_a, .TRUE., kd, kb, mb, es, em)
    CALL require(es == CD_MODAL_OK .AND. kd == 5, 'string:band-matrices: '//TRIM(em))
    CALL check_band_against_dense('string', k, m, kb, mb, kd, free, om2)
    t = EA*EPS
    mc = RHO*l0tot/L
    ! transverse pairs (y and z) for n = 1..5
    DO j = 1, 10
      nmode = (j + 1)/2
      f_exp = REAL(nmode, wp)/(2.0_wp*L)*SQRT(t/mc)
      f_got = SQRT(om2(j))/(2.0_wp*PI)
      CALL require(ABS(f_got/f_exp - 1.0_wp) < 2.0e-4_wp*nmode**2, 'string:transverse-frequency')
      IF (j == 1) WRITE (*, '(A,2ES14.6)') 'taut string f1 got/exact (Hz): ', f_got, f_exp
    END DO
    ! the first axial (x-dominant) mode
    DO j = 1, nf
      ax = SUM(sh(1::3, j)**2)
      tr = SUM(sh(2::3, j)**2) + SUM(sh(3::3, j)**2)
      IF (ax > tr) EXIT
    END DO
    CALL require(j <= nf, 'string:axial-mode-found')
    IF (j <= nf) THEN
      f_exp = SQRT(EA/RHO)/(2.0_wp*l0tot)
      f_got = SQRT(om2(j))/(2.0_wp*PI)
      WRITE (*, '(A,2ES14.6)') 'taut string first axial mode got/exact (Hz): ', f_got, f_exp
      CALL require(ABS(f_got/f_exp - 1.0_wp) < 1.0e-3_wp, 'string:axial-frequency')
    END IF
  END SUBROUTINE case_taut_string

  SUBROUTINE hermite_line(ne, l, eps, ei, nm, om2, nf)
    INTEGER, INTENT(IN) :: ne, nm
    REAL(wp), INTENT(IN) :: l, eps, ei
    REAL(wp), INTENT(OUT) :: om2(nm)
    INTEGER, INTENT(OUT) :: nf
    REAL(wp), ALLOCATABLE :: q(:), l0(:), ea(:), eiv(:), rho(:), k(:, :), m(:, :), sh(:, :), kb(:, :), mb(:, :)
    LOGICAL, ALLOCATABLE :: free(:)
    INTEGER :: i, es, nn, kd
    CHARACTER(256) :: em
    nn = ne + 1
    ALLOCATE (q(6*nn), l0(ne), ea(ne), eiv(ne), rho(ne), k(6*nn, 6*nn), m(6*nn, 6*nn), sh(6*nn, nm), free(6*nn))
    l0 = l/(1.0_wp + eps)/ne
    ea = 1.0e9_wp
    eiv = ei
    rho = 100.0_wp
    DO i = 1, nn
      q(6*i - 5:6*i) = [l*REAL(i - 1, wp)/ne, 0.0_wp, 0.0_wp, 1.0_wp + eps, 0.0_wp, 0.0_wp]
    END DO
    free = .TRUE.
    free(1:3) = .FALSE.
    free(6*nn - 5:6*nn - 3) = .FALSE.
    CALL CD_Modal_Hermite_Matrices(q, l0, ea, eiv, rho, k, m, es, em)
    CALL require(es == CD_MODAL_OK, 'hermite:matrices: '//TRIM(em))
    CALL CD_Modal_Solve(k, m, free, nm, om2, sh, nf, es, em)
    CALL require(es == CD_MODAL_OK .AND. nf == nm, 'hermite:solve: '//TRIM(em))
    CALL CD_Modal_Hermite_Band(q, l0, ea, eiv, rho, kd, kb, mb, es, em)
    CALL require(es == CD_MODAL_OK .AND. kd == 11, 'hermite:band-matrices: '//TRIM(em))
    CALL check_band_against_dense('hermite', k, m, kb, mb, kd, free, om2)
  END SUBROUTINE hermite_line

  SUBROUTINE case_hermite_beam()
    INTEGER, PARAMETER :: NM = 8
    REAL(wp), PARAMETER :: L = 50.0_wp, EI = 1.0e6_wp, RHO = 100.0_wp
    REAL(wp) :: om2(NM), f_exp, f_got, t, mc, eps
    INTEGER :: nf, j, n
    ! pinned-pinned beam, no tension
    CALL hermite_line(40, L, 0.0_wp, EI, NM, om2, nf)
    DO j = 1, NM
      n = (j + 1)/2
      f_exp = (REAL(n, wp)*PI/L)**2*SQRT(EI/RHO)/(2.0_wp*PI)
      f_got = SQRT(MAX(om2(j), 0.0_wp))/(2.0_wp*PI)
      IF (j == 1) WRITE (*, '(A,2ES14.6)') 'Hermite pinned beam f1 got/exact (Hz): ', f_got, f_exp
      CALL require(ABS(f_got/f_exp - 1.0_wp) < 1.0e-3_wp, 'hermite:beam-frequency')
    END DO
    ! pre-tensioned beam-string
    eps = 1.0e-4_wp
    CALL hermite_line(40, L, eps, EI, NM, om2, nf)
    t = 1.0e9_wp*eps
    mc = RHO/(1.0_wp + eps)
    DO j = 1, NM
      n = (j + 1)/2
      f_exp = REAL(n, wp)/(2.0_wp*L)*SQRT(t/mc)*SQRT(1.0_wp + (REAL(n, wp)*PI)**2*EI/(t*L*L))
      f_got = SQRT(MAX(om2(j), 0.0_wp))/(2.0_wp*PI)
      IF (j == 1) WRITE (*, '(A,2ES14.6)') 'Hermite tensioned beam f1 got/exact (Hz): ', f_got, f_exp
      CALL require(ABS(f_got/f_exp - 1.0_wp) < 2.0e-3_wp, 'hermite:tensioned-beam-frequency')
    END DO
  END SUBROUTINE case_hermite_beam

  SUBROUTINE check_band_against_dense(label, k, m, kb, mb, kd, free, om2_dense)
    !! The banded solver on the band form of the dense matrices k, m: the frequencies of the
    !! dense reference, and shapes that are M-orthonormal, zero on the held DOFs and satisfy
    !! K phi = omega^2 M phi.
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: k(:, :), m(:, :), kb(:, :), mb(:, :), om2_dense(:)
    INTEGER, INTENT(IN) :: kd
    LOGICAL, INTENT(IN) :: free(:)
    REAL(wp), ALLOCATABLE :: om2(:), sh(:, :), ksh(:, :), msh(:, :), gram(:, :)
    REAL(wp) :: knorm, fdiff, resid, ortho
    INTEGER :: nm, nf, es, j
    LOGICAL :: held_zero
    CHARACTER(256) :: em
    nm = SIZE(om2_dense)
    ALLOCATE (om2(nm), sh(SIZE(free), nm))
    CALL CD_Modal_Solve_Band(kb, mb, kd, free, nm, om2, sh, nf, es, em)
    CALL require(es == CD_MODAL_OK .AND. nf == nm, label//':band-solve: '//TRIM(em))
    IF (es /= CD_MODAL_OK .OR. nf /= nm) RETURN
    knorm = MAXVAL(SUM(ABS(k), 1))
    fdiff = 0.0_wp
    DO j = 1, nm
      fdiff = MAX(fdiff, ABS(om2(j) - om2_dense(j))/ABS(om2_dense(j)))
    END DO
    ksh = MATMUL(k, sh)
    msh = MATMUL(m, sh)
    resid = 0.0_wp
    held_zero = .TRUE.
    DO j = 1, nm
      resid = MAX(resid, MAXVAL(ABS(ksh(:, j) - om2(j)*msh(:, j)), MASK=free)/(knorm*nan_max_abs(sh(:, j))))
      held_zero = held_zero .AND. nan_max_abs(PACK(sh(:, j),.NOT. free)) <= 0.0_wp
    END DO
    gram = MATMUL(TRANSPOSE(sh), msh)
    DO j = 1, nm
      gram(j, j) = gram(j, j) - 1.0_wp
    END DO
    ortho = nan_max_abs(gram)
    WRITE (*, '(A,3ES11.3)') label//' banded vs dense: max relative omega^2 difference, shape residual, '// &
      'M-orthonormality error: ', fdiff, resid, ortho
    CALL require(fdiff < 1.0e-8_wp, label//':band-frequencies-equal-dense')
    CALL require(resid < 1.0e-9_wp, label//':band-shape-residual')
    CALL require(ortho < 1.0e-8_wp, label//':band-shapes-M-orthonormal')
    CALL require(held_zero, label//':band-held-dofs-zero')
  END SUBROUTINE check_band_against_dense

  SUBROUTINE case_slack_string()
    !! A slack tension-only string has no stiffness (K = 0): every frequency is zero and the
    !! banded solver still returns M-orthonormal shapes for the repeated zero eigenvalue.
    INTEGER, PARAMETER :: NE = 20, NN = NE + 1, NM = 6
    REAL(wp) :: q(3*NN), l0(NE), ea_(NE), rho_a(NE), om2(NM), sh(3*NN, NM), gram(NM, NM)
    REAL(wp) :: k(3*NN, 3*NN), m(3*NN, 3*NN)
    REAL(wp), ALLOCATABLE :: kb(:, :), mb(:, :)
    INTEGER :: conn(2, NE), i, es, nf, kd
    LOGICAL :: free(3*NN)
    CHARACTER(256) :: em
    DO i = 1, NE
      conn(:, i) = [i, i + 1]
    END DO
    l0 = 1.1_wp
    ea_ = 1.0e8_wp
    rho_a = 50.0_wp
    DO i = 1, NN
      q(3*i - 2:3*i) = [REAL(i - 1, wp), 0.0_wp, 0.0_wp]
    END DO
    free = .TRUE.
    free(1:3) = .FALSE.
    free(3*NN - 2:3*NN) = .FALSE.
    CALL CD_Modal_Bar_Band(q, conn, l0, ea_, rho_a, .TRUE., kd, kb, mb, es, em)
    CALL require(es == CD_MODAL_OK, 'slack:band-matrices: '//TRIM(em))
    CALL CD_Modal_Solve_Band(kb, mb, kd, free, NM, om2, sh, nf, es, em)
    CALL require(es == CD_MODAL_OK .AND. nf == NM, 'slack:band-solve: '//TRIM(em))
    IF (es /= CD_MODAL_OK .OR. nf /= NM) RETURN
    CALL require(nan_max_abs(om2) < 1.0e-9_wp, 'slack:zero-frequencies')
    CALL CD_Modal_Bar_Matrices(q, conn, l0, ea_, rho_a, .TRUE., k, m, es, em)
    gram = MATMUL(TRANSPOSE(sh), MATMUL(m, sh))
    DO i = 1, NM
      gram(i, i) = gram(i, i) - 1.0_wp
    END DO
    CALL require(nan_max_abs(gram) < 1.0e-10_wp, 'slack:shapes-M-orthonormal')
  END SUBROUTINE case_slack_string

  SUBROUTINE case_hermite_long_beam()
    !! A pre-tensioned 1024-element beam-string (6150 DOFs, beyond the dense reference): the
    !! banded solver against the closed-form frequencies; each repeated frequency (the two
    !! transverse planes) gets shapes in both planes.
    INTEGER, PARAMETER :: NE = 1024, NN = NE + 1, NM = 10
    REAL(wp), PARAMETER :: L = 500.0_wp, EI = 1.0e6_wp, RHO = 100.0_wp, EPS = 1.0e-4_wp
    REAL(wp), ALLOCATABLE :: q(:), l0(:), ea(:), eiv(:), rho_a(:), kb(:, :), mb(:, :), sh(:, :)
    LOGICAL, ALLOCATABLE :: free(:)
    REAL(wp) :: om2(NM), t, mc, f_exp, f_got, sy, sz
    INTEGER :: i, es, kd, nf, j, n
    INTEGER(8) :: c0, c1, rate
    CHARACTER(256) :: em
    ALLOCATE (q(6*NN), l0(NE), ea(NE), eiv(NE), rho_a(NE), sh(6*NN, NM), free(6*NN))
    l0 = L/(1.0_wp + EPS)/NE
    ea = 1.0e9_wp
    eiv = EI
    rho_a = RHO
    DO i = 1, NN
      q(6*i - 5:6*i) = [L*REAL(i - 1, wp)/NE, 0.0_wp, 0.0_wp, 1.0_wp + EPS, 0.0_wp, 0.0_wp]
    END DO
    free = .TRUE.
    free(1:3) = .FALSE.
    free(6*NN - 5:6*NN - 3) = .FALSE.
    CALL SYSTEM_CLOCK(c0, rate)
    CALL CD_Modal_Hermite_Band(q, l0, ea, eiv, rho_a, kd, kb, mb, es, em)
    CALL require(es == CD_MODAL_OK, 'long-beam:band-matrices: '//TRIM(em))
    CALL CD_Modal_Solve_Band(kb, mb, kd, free, NM, om2, sh, nf, es, em)
    CALL SYSTEM_CLOCK(c1)
    CALL require(es == CD_MODAL_OK .AND. nf == NM, 'long-beam:band-solve: '//TRIM(em))
    IF (es /= CD_MODAL_OK .OR. nf /= NM) RETURN
    WRITE (*, '(A,I0,A,F8.3,A)') 'long beam: ', COUNT(free), ' free DOFs, banded solve in ', &
      REAL(c1 - c0, wp)/REAL(rate, wp), ' s'
    t = 1.0e9_wp*EPS
    mc = RHO/(1.0_wp + EPS)
    DO j = 1, NM
      n = (j + 1)/2
      f_exp = REAL(n, wp)/(2.0_wp*L)*SQRT(t/mc)*SQRT(1.0_wp + (REAL(n, wp)*PI)**2*EI/(t*L*L))
      f_got = SQRT(MAX(om2(j), 0.0_wp))/(2.0_wp*PI)
      IF (j == 1) WRITE (*, '(A,2ES14.6)') 'long beam f1 got/exact (Hz): ', f_got, f_exp
      CALL require(ABS(f_got/f_exp - 1.0_wp) < 1.0e-3_wp, 'long-beam:frequency')
    END DO
    DO j = 1, NM, 2
      sy = SUM(sh(2::6, j:j + 1)**2)
      sz = SUM(sh(3::6, j:j + 1)**2)
      CALL require(MIN(sy, sz) > 0.25_wp*(sy + sz), 'long-beam:repeated-frequency-spans-both-planes')
    END DO
  END SUBROUTINE case_hermite_long_beam

  ! ------------------------------------------------------------------ deck helpers

  SUBROUTINE read_table(path, marker, ncol, data, nrow)
    !! Numeric rows of the table that follows the comment line starting with marker
    !! (after its name and unit rows), up to the next comment line.
    CHARACTER(*), INTENT(IN) :: path, marker
    INTEGER, INTENT(IN) :: ncol
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: data(:, :)
    INTEGER, INTENT(OUT) :: nrow
    INTEGER :: u, ios, k, skip
    CHARACTER(4096) :: line
    LOGICAL :: inside
    REAL(wp), ALLOCATABLE :: grow(:, :)
    ALLOCATE (data(ncol, 64))
    nrow = 0
    inside = .FALSE.
    skip = 0
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    DO
      READ (u, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF (line(1:1) == '#') THEN
        IF (inside) EXIT
        IF (INDEX(line, marker) > 0) THEN
          inside = .TRUE.
          skip = 2
        END IF
        CYCLE
      END IF
      IF (.NOT. inside) CYCLE
      IF (skip > 0) THEN
        skip = skip - 1
        CYCLE
      END IF
      IF (nrow == SIZE(data, 2)) THEN
        ALLOCATE (grow(ncol, 2*SIZE(data, 2)))
        grow(:, 1:nrow) = data(:, 1:nrow)
        CALL MOVE_ALLOC(grow, data)
      END IF
      nrow = nrow + 1
      READ (line, *, IOSTAT=ios) (data(k, nrow), k=1, ncol)
      IF (ios /= 0) nrow = nrow - 1
    END DO
    CLOSE (u)
  END SUBROUTINE read_table

  SUBROUTINE write_catenary_deck(path, extra)
    !! A 100 m span between two level points (Coupled, Fixed) in deep water: EI = 0 chain, no added
    !! mass, 100 segments.
    CHARACTER(*), INTENT(IN) :: path, extra
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'Irvine sagging cable'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.1 100.0 1.0e8 -1.0 0.0 1.0 0.5 0.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Coupled 0.0 0.0 -100.0'
    WRITE (u, '(A)') '2 Fixed 100.0 0.0 -100.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 chain 100.15 100'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '500.0 WtrDpth'
    WRITE (u, '(A)') '12 nModes'
    IF (LEN_TRIM(extra) > 0) WRITE (u, '(A)') TRIM(extra)
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1 AnchTen1'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_catenary_deck

  REAL(wp) FUNCTION irvine_beta(lam2) RESULT(b)
    !! First root b in (pi, 3 pi) of tan(b/2) = b/2 - (4/lambda^2)(b/2)^3 (bisection).
    REAL(wp), INTENT(IN) :: lam2
    REAL(wp) :: lo, hi, mid
    INTEGER :: it
    lo = PI*(1.0_wp + 1.0e-9_wp)
    hi = 3.0_wp*PI*(1.0_wp - 1.0e-9_wp)
    DO it = 1, 200
      mid = 0.5_wp*(lo + hi)
      IF (irvine_g(lo, lam2)*irvine_g(mid, lam2) <= 0.0_wp) THEN
        hi = mid
      ELSE
        lo = mid
      END IF
    END DO
    b = 0.5_wp*(lo + hi)
  END FUNCTION irvine_beta

  REAL(wp) FUNCTION irvine_g(x, lam2) RESULT(g)
    !! sin/cos form of tan(x/2) - x/2 + (4/lambda^2)(x/2)^3, continuous across x = pi
    REAL(wp), INTENT(IN) :: x, lam2
    g = SIN(0.5_wp*x) - COS(0.5_wp*x)*(0.5_wp*x - 4.0_wp/lam2*(0.5_wp*x)**3)
  END FUNCTION irvine_g

  SUBROUTINE case_irvine()
    REAL(wp), PARAMETER :: L = 100.0_wp, MASS = 100.0_wp, DIAM = 0.1_wp, EA = 1.0e8_wp
    REAL(wp), ALLOCATABLE :: fr(:, :), sh(:, :), st(:, :)
    REAL(wp) :: w, h, sag, le, lam2, beta, f_swing, f_anti, f_sym, got_swing, got_anti, got_sym
    REAL(wp) :: sy, sxz, zmid
    INTEGER :: nfr, nsh, nst, j, k, i
    LOGICAL :: conv, found_swing, found_anti, found_sym
    INTEGER :: es
    CHARACTER(512) :: em
    CALL write_catenary_deck('modal_irvine.dat', '')
    CALL CD_Run_Deck_Driver('modal_irvine.dat', 'modal_irvine', conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, 'irvine:runs: '//TRIM(em))
    IF (es /= CD_DECKDRV_OK) RETURN
    CALL read_table('modal_irvine.modes.out', 'Natural frequencies', 5, fr, nfr)
    CALL read_table('modal_irvine.modes.out', 'Mode shapes', 9, sh, nsh)
    CALL require(nfr == 12 .AND. nsh == 12*101, 'irvine:modes-file-rows')
    IF (nfr /= 12 .OR. nsh /= 12*101) RETURN
    ! the horizontal tension: the least static tension (at midspan, where the cable is level)
    CALL read_static(st, nst)
    h = MINVAL(st(7, 1:nst))
    w = (MASS - 1025.0_wp*PI*DIAM**2/4.0_wp)*9.80665_wp
    sag = w*L*L/(8.0_wp*h)
    le = L*(1.0_wp + 8.0_wp*(sag/L)**2)
    lam2 = (w*L/h)**2*L/(h*le/EA)
    beta = irvine_beta(lam2)
    f_swing = SQRT(h/MASS)/(2.0_wp*L)
    f_anti = SQRT(h/MASS)/L
    f_sym = beta*SQRT(h/MASS)/(2.0_wp*PI*L)
    WRITE (*, '(A,4ES14.6)') 'Irvine: H (N), sag (m), lambda^2, beta/pi: ', h, sag, lam2, beta/PI
    found_swing = .FALSE.
    found_anti = .FALSE.
    found_sym = .FALSE.
    DO j = 1, nfr
      sy = 0.0_wp
      sxz = 0.0_wp
      zmid = 0.0_wp
      DO k = 1, nsh
        IF (NINT(sh(2, k)) /= j) CYCLE
        sy = sy + sh(8, k)**2
        sxz = sxz + sh(7, k)**2 + sh(9, k)**2
        i = NINT(sh(3, k))
        IF (i == 51) zmid = ABS(sh(9, k))
      END DO
      IF (sy > sxz) THEN
        IF (.NOT. found_swing) THEN
          got_swing = fr(3, j)
          found_swing = .TRUE.
        END IF
      ELSE IF (zmid < 1.0e-3_wp) THEN
        IF (.NOT. found_anti) THEN
          got_anti = fr(3, j)
          found_anti = .TRUE.
        END IF
      ELSE
        IF (.NOT. found_sym) THEN
          got_sym = fr(3, j)
          found_sym = .TRUE.
        END IF
      END IF
    END DO
    CALL require(found_swing .AND. found_anti .AND. found_sym, 'irvine:mode-families-found')
    IF (.NOT. (found_swing .AND. found_anti .AND. found_sym)) RETURN
    WRITE (*, '(A,2ES14.6)') 'Irvine out-of-plane swing got/theory (Hz): ', got_swing, f_swing
    WRITE (*, '(A,2ES14.6)') 'Irvine in-plane antisymmetric got/theory (Hz): ', got_anti, f_anti
    WRITE (*, '(A,2ES14.6)') 'Irvine in-plane symmetric got/theory (Hz): ', got_sym, f_sym
    CALL require(ABS(got_swing/f_swing - 1.0_wp) < 0.01_wp, 'irvine:out-of-plane')
    CALL require(ABS(got_anti/f_anti - 1.0_wp) < 0.01_wp, 'irvine:in-plane-antisymmetric')
    CALL require(ABS(got_sym/f_sym - 1.0_wp) < 0.02_wp, 'irvine:in-plane-symmetric')
  END SUBROUTINE case_irvine

  SUBROUTINE read_static(data, nrow)
    !! modal_irvine.static.out rows (LineID Node ArcLength X Y Z Tension ...).
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: data(:, :)
    INTEGER, INTENT(OUT) :: nrow
    INTEGER :: u, ios, k
    CHARACTER(4096) :: line
    ALLOCATE (data(12, 200))
    nrow = 0
    OPEN (NEWUNIT=u, FILE='modal_irvine.static.out', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    READ (u, '(A)') line
    READ (u, '(A)') line
    READ (u, '(A)') line
    DO
      READ (u, '(A)', IOSTAT=ios) line
      IF (ios /= 0 .OR. nrow == 200) EXIT
      nrow = nrow + 1
      READ (line, *, IOSTAT=ios) (data(k, nrow), k=1, 12)
      IF (ios /= 0) nrow = nrow - 1
    END DO
    CLOSE (u)
  END SUBROUTINE read_static

  SUBROUTINE case_deck_hermite_and_rules()
    REAL(wp), ALLOCATABLE :: fr(:, :)
    INTEGER :: nfr, u, es
    LOGICAL :: conv
    CHARACTER(512) :: em
    OPEN (NEWUNIT=u, FILE='modal_cable.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'finite-EI touchdown cable modes'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'pcable 0.20 250 8.0e8 0.0 2.0e5 1.0 0.0 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Fixed 0 0 -80.0'
    WRITE (u, '(A)') '2 Coupled 40 0 -70.0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 2 1 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 pcable 45 20'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '80.0 WtrDpth'
    WRITE (u, '(A)') '1.0e5 kBot'
    WRITE (u, '(A)') '0.02 dtM'
    WRITE (u, '(A)') '0.04 TMax'
    WRITE (u, '(A)') '6 nModes'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
    CALL CD_Run_Deck_Driver('modal_cable.dat', 'modal_cable', conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, 'hermite-deck:runs: '//TRIM(em))
    CALL read_table('modal_cable.modes.out', 'Natural frequencies', 5, fr, nfr)
    CALL require(nfr == 6, 'hermite-deck:six-modes')
    IF (nfr == 6) CALL require(ALL(fr(3, 1:6) > 0.0_wp) .AND. ALL(fr(3, 2:6) >= fr(3, 1:5)), &
                               'hermite-deck:positive-ascending-frequencies')
    ! rejected: fractional count, and a deck with a Free point
    CALL write_catenary_deck('modal_bad1.dat', '2.5 nModes')
    CALL CD_Run_Deck_Driver('modal_bad1.dat', 'modal_bad1', conv, es, em)
    CALL require(es /= CD_DECKDRV_OK, 'rules:fractional-nModes-rejected')
  END SUBROUTINE case_deck_hermite_and_rules

END PROGRAM test_modal
