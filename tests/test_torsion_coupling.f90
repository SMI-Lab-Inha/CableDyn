! File: tests/test_torsion_coupling.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_torsion_coupling
  !! Torque-tension coupling of a finite-EI cable (LINE TYPES GJ column GJ|K: the section laws
  !! T = EA eps + K tau and M = GJ tau + K eps).
  !!   1. The kernel: with the coupling, CD_HermiteTorsion_Line's effective twist
  !!      Theta - sum (K/GJ)(l_e - l0_e) has the gradient and Hessian of its finite differences
  !!      on a curved line.
  !!   2. A straight neutral cable stretched by eps between its two held ends, both ends rigid in
  !!      torsion: the twist cannot relax, so T = EA eps and M = K eps.
  !!   3. The same cable free in torsion at one end: the stretch twists it, so T = (EA - K^2/GJ) eps
  !!      and M = 0.
  !!   4. Both ends rigid with an imposed twist Phi: tau = Phi/L, so T = EA eps + K Phi/L and
  !!      M = GJ Phi/L + K eps; the twist channel reads the element twists.
  !!   5. K = 0 written as GJ|0 is the uncoupled cable bit for bit; K^2 >= EA GJ fails closed.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_HermiteTorsion, ONLY: CD_HermiteTorsion_Line, CD_HTORS_OK, CD_HTORS_KBAND
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, RHOW = 1025.0_wp
  ! the section: EA 1e7 N, EI 1e6 N m2, GJ 5e4 N m2, K 4e5 N m (K^2 / EA GJ = 0.32)
  REAL(wp), PARAMETER :: EA = 1.0e7_wp, GJ = 5.0e4_wp, KTT = 4.0e5_wp, SPAN = 100.0_wp, L0 = 99.0_wp
  INTEGER :: nfail = 0

  CALL check_kernel()
  CALL check_rigid_ends()
  CALL check_free_end()
  CALL check_pretwist()
  CALL check_zero_and_refusal()
  CALL check_dynamic()
  IF (nfail > 0) THEN
    PRINT '(I0,A)', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  PRINT '(A)', 'all torque-tension coupling checks passed'

CONTAINS

  SUBROUTINE check_kernel()
    !! Finite differences of the effective twist on a curved, stretched 6-node line.
    INTEGER, PARAMETER :: NN = 6, ND = 6*NN, KB = CD_HTORS_KBAND
    REAL(wp) :: q(ND), le(NN - 1), ends(3, 4), cpl(NN - 1), th, thp, thm, g(ND), gp(ND), gm(ND)
    REAL(wp) :: hb(2*KB + 1, ND), hfd(ND, ND), herr, gerr, h, s, r(3), m(3)
    INTEGER :: k, i, j, es
    CHARACTER(256) :: em
    DO k = 1, NN
      s = REAL(k - 1, wp)*2.0_wp
      r = [s, 0.3_wp*SIN(0.4_wp*s), 0.2_wp*s*s/10.0_wp]
      m = [1.0_wp, 0.12_wp*COS(0.4_wp*s), 0.04_wp*s]
      m = 1.01_wp*m/NORM2(m)
      q(6*k - 5:6*k) = [r, m]
    END DO
    le = 2.0_wp
    cpl = [0.3_wp, -0.2_wp, 0.5_wp, 0.1_wp, 0.4_wp]
    ends(:, 1) = q(4:6)/NORM2(q(4:6))
    ends(:, 2) = unit(cross([0.0_wp, 0.0_wp, 1.0_wp], ends(:, 1)))
    ends(:, 3) = q(ND - 2:ND)/NORM2(q(ND - 2:ND))
    ends(:, 4) = unit(cross([0.0_wp, 0.0_wp, 1.0_wp], ends(:, 3)))
    hb = 0.0_wp
    CALL CD_HermiteTorsion_Line(q, le, ends, th, g, es, em, hband=hb, couple=cpl)
    CALL require(es == CD_HTORS_OK, 'kernel with coupling evaluates: '//TRIM(em))
    h = 1.0e-6_wp
    gerr = 0.0_wp
    herr = 0.0_wp
    DO j = 1, ND
      q(j) = q(j) + h
      CALL CD_HermiteTorsion_Line(q, le, ends, thp, gp, es, em, couple=cpl)
      q(j) = q(j) - 2*h
      CALL CD_HermiteTorsion_Line(q, le, ends, thm, gm, es, em, couple=cpl)
      q(j) = q(j) + h
      gerr = MAX(gerr, ABS((thp - thm)/(2*h) - g(j)))
      hfd(:, j) = (gp - gm)/(2*h)
    END DO
    DO j = 1, ND
      DO i = MAX(1, j - KB), MIN(ND, j + KB)
        herr = MAX(herr, ABS(hb(KB + 1 + i - j, j) - hfd(i, j)))
      END DO
    END DO
    PRINT '(A,ES10.3,A,ES10.3)', '  kernel: max |grad - FD| ', gerr, ', max |Hessian - FD| ', herr
    CALL require(gerr <= 1.0e-7_wp .AND. herr <= 1.0e-6_wp, 'kernel: the effective twist''s derivatives')
  END SUBROUTINE check_kernel

  SUBROUTINE check_rigid_ends()
    REAL(wp) :: v(4), eps
    CALL straight_case('tcpl_rigid', 'Rigid 1 0 0 Rigid 0 0 1 0', 'Rigid 1 0 0 Rigid 0 0 1 0', KTT, v)
    eps = SPAN/L0 - 1.0_wp
    PRINT '(A,2ES15.7,A,2ES15.7)', '  rigid ends: T, M ', v(1), v(2), '  expected ', EA*eps, KTT*eps
    CALL require(ABS(v(1) - EA*eps) <= 1.0e-5_wp*EA*eps, 'rigid ends: T = EA eps')
    CALL require(ABS(v(2) - KTT*eps) <= 1.0e-5_wp*KTT*eps, 'rigid ends: M = K eps')
  END SUBROUTINE check_rigid_ends

  SUBROUTINE check_free_end()
    REAL(wp) :: v(4), eps
    CALL straight_case('tcpl_free', 'Rigid 1 0 0 Free 0 0 1', 'Rigid 1 0 0 Rigid 0 0 1 0', KTT, v, torsion=.FALSE.)
    eps = SPAN/L0 - 1.0_wp
    PRINT '(A,ES15.7,A,ES15.7)', '  free end: T ', v(1), '  expected ', (EA - KTT**2/GJ)*eps
    CALL require(ABS(v(1) - (EA - KTT**2/GJ)*eps) <= 1.0e-5_wp*EA*eps, 'free end: T = (EA - K^2/GJ) eps')
  END SUBROUTINE check_free_end

  SUBROUTINE check_pretwist()
    REAL(wp) :: v(4), eps, phi, t_ref, m_ref, tw_ref
    phi = 90.0_wp*PI/180.0_wp
    CALL straight_case('tcpl_twist', 'Rigid 1 0 0 Rigid 0 0 1 0', 'Rigid 1 0 0 Rigid 0 0 1 90', KTT, v)
    eps = SPAN/L0 - 1.0_wp
    t_ref = EA*eps + KTT*phi/L0
    m_ref = GJ*phi/L0 + KTT*eps
    tw_ref = 90.0_wp
    PRINT '(A,3ES15.7,A,3ES15.7)', '  pretwist: T, M, Twist ', v(1:3), '  expected ', t_ref, m_ref, tw_ref
    CALL require(ABS(v(1) - t_ref) <= 1.0e-5_wp*t_ref, 'pretwist: T = EA eps + K Phi/L')
    CALL require(ABS(v(2) - m_ref) <= 1.0e-5_wp*m_ref, 'pretwist: M = GJ Phi/L + K eps')
    CALL require(ABS(v(3) - tw_ref) <= 1.0e-4_wp, 'pretwist: Twist<L> is the imposed twist')
    CALL require(ABS(v(4) - v(1)) <= 1.0e-5_wp*v(1), 'pretwist: the mid-line tension equals the end tension')
  END SUBROUTINE check_pretwist

  SUBROUTINE check_zero_and_refusal()
    REAL(wp) :: v0(4), v1(4)
    LOGICAL :: conv
    INTEGER :: es
    CHARACTER(512) :: em
    CALL straight_case('tcpl_k0', 'Rigid 1 0 0 Rigid 0 0 1 0', 'Rigid 1 0 0 Rigid 0 0 1 30', -1.0_wp, v0)
    CALL straight_case('tcpl_k0b', 'Rigid 1 0 0 Rigid 0 0 1 0', 'Rigid 1 0 0 Rigid 0 0 1 30', 0.0_wp, v1)
    CALL require(ALL(v0 == v1), 'GJ|0 is the uncoupled cable bit for bit')
    CALL write_case('tcpl_bad.dat', 'Rigid 1 0 0 Rigid 0 0 1 0', 'Rigid 1 0 0 Rigid 0 0 1 0', 1.0e6_wp)
    CALL CD_Run_Deck_Driver('tcpl_bad.dat', 'tcpl_bad', conv, es, em)
    CALL require(es /= CD_DECKDRV_OK .AND. INDEX(em, 'K^2 < EA GJ') > 0, 'K^2 >= EA GJ fails closed: '//TRIM(em))
  END SUBROUTINE check_zero_and_refusal

  SUBROUTINE check_dynamic()
    !! The march: End A moves 0.1 m away from End B and rolls -90 deg (Phi + 90 deg) over 20 s on a
    !! C2 quintic, then holds 5 s; the line, far stiffer axially than the ramp is slow, ends at the
    !! static state of the 100.1 m span with Phi = 90 deg: T = EA eps + K Phi/L, M = GJ Phi/L + K eps.
    INTEGER, PARAMETER :: NSTEP = 500
    REAL(wp), PARAMETER :: TR = 20.0_wp, DTM = 0.05_wp, DX = 0.1_wp
    REAL(wp) :: t, tau, s, sd, sdd, v(4), eps, phi, t_ref, m_ref, last(5)
    INTEGER :: u, k, ios
    LOGICAL :: conv
    INTEGER :: es
    CHARACTER(512) :: em
    CHARACTER(4096) :: buf
    OPEN (NEWUNIT=u, FILE='tcpl_motion.txt', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '# time id x y z vx vy vz ax ay az roll'
    DO k = 0, NSTEP
      t = DTM*REAL(k, wp)
      tau = MIN(t/TR, 1.0_wp)
      s = 10.0_wp*tau**3 - 15.0_wp*tau**4 + 6.0_wp*tau**5
      sd = (30.0_wp*tau**2 - 60.0_wp*tau**3 + 30.0_wp*tau**4)/TR
      sdd = (60.0_wp*tau - 180.0_wp*tau**2 + 120.0_wp*tau**3)/TR**2
      IF (t > TR) THEN
        sd = 0.0_wp
        sdd = 0.0_wp
      END IF
      WRITE (u, '(ES16.8,A,ES16.8,A,3ES16.8,A,ES16.8)') t, ' 1 ', -DX*s, ' 0.0 -50.0 ', -DX*sd, 0.0_wp, 0.0_wp, &
        ' 0 0 0 ', -90.0_wp*s
    END DO
    CLOSE (u)
    CALL write_case('tcpl_dyn.dat', 'Rigid 1 0 0 Rigid 0 0 1 0', 'Rigid 1 0 0 Rigid 0 0 1 0', KTT, dynamic=.TRUE.)
    CALL CD_Run_Deck_Driver('tcpl_dyn.dat', 'tcpl_dyn', conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, 'the coupled march runs: '//TRIM(em))
    IF (es /= CD_DECKDRV_OK) RETURN
    OPEN (NEWUNIT=u, FILE='tcpl_dyn.out', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    READ (u, '(A)') buf
    READ (u, '(A)') buf
    last = 0.0_wp
    DO
      READ (u, *, IOSTAT=ios) v(1), v(2:4)
      IF (ios /= 0) EXIT
      last(1:4) = v
    END DO
    CLOSE (u)
    eps = (SPAN + DX)/L0 - 1.0_wp
    phi = 90.0_wp*PI/180.0_wp
    t_ref = EA*eps + KTT*phi/L0
    m_ref = GJ*phi/L0 + KTT*eps
    PRINT '(A,3ES15.7,A,2ES15.7)', '  march end: t, T, M ', last(1:3), '  expected ', t_ref, m_ref
    CALL require(ABS(last(2) - t_ref) <= 1.0e-6_wp*t_ref, 'march: T at the end of the ramp is the static one')
    CALL require(ABS(last(3) - m_ref) <= 1.0e-6_wp*m_ref, 'march: M at the end of the ramp is the static one')
  END SUBROUTINE check_dynamic

  SUBROUTINE straight_case(root, conn_a, conn_b, k, v, torsion)
    !! The straight case at TMax = 0: v = FairTen1, Torq1N1, Twist1, Ten1N21 (Twist and Torq zero
    !! for a line without torsion).
    CHARACTER(*), INTENT(IN) :: root, conn_a, conn_b
    REAL(wp), INTENT(IN) :: k
    REAL(wp), INTENT(OUT) :: v(4)
    LOGICAL, INTENT(IN), OPTIONAL :: torsion
    LOGICAL :: conv, tors
    INTEGER :: es, u, ios
    CHARACTER(512) :: em
    CHARACTER(4096) :: buf
    REAL(wp) :: row(5)
    tors = .TRUE.
    IF (PRESENT(torsion)) tors = torsion
    v = 0.0_wp
    CALL write_case(root//'.dat', conn_a, conn_b, k, tors)
    CALL CD_Run_Deck_Driver(root//'.dat', root, conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, root//' runs: '//TRIM(em))
    IF (es /= CD_DECKDRV_OK) RETURN
    OPEN (NEWUNIT=u, FILE=root//'.out', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    READ (u, '(A)') buf
    READ (u, '(A)') buf
    row = 0.0_wp
    IF (tors) THEN
      READ (u, *, IOSTAT=ios) row(1:5)
      v = row(2:5)
    ELSE
      READ (u, *, IOSTAT=ios) row(1:3)
      v(1) = row(2)
      v(4) = row(3)
    END IF
    CLOSE (u)
  END SUBROUTINE straight_case

  SUBROUTINE write_case(path, conn_a, conn_b, k, torsion, dynamic)
    !! A 99 m neutral cable from a Coupled End A to a Fixed End B 100 m apart along +x, 40 elements, TMax 0; the GJ
    !! column carries GJ|K (k < 0 writes the plain GJ).
    CHARACTER(*), INTENT(IN) :: path, conn_a, conn_b
    REAL(wp), INTENT(IN) :: k
    LOGICAL, INTENT(IN), OPTIONAL :: torsion, dynamic
    INTEGER :: u
    CHARACTER(64) :: gjcol
    LOGICAL :: tors, dyn
    dyn = .FALSE.
    IF (PRESENT(dynamic)) dyn = dynamic
    tors = .TRUE.
    IF (PRESENT(torsion)) tors = torsion
    IF (k < 0.0_wp) THEN
      WRITE (gjcol, '(ES12.5)') GJ
    ELSE
      WRITE (gjcol, '(ES12.5,A,ES12.5)') GJ, '|', k
    END IF
    gjcol = remove_blanks(gjcol)
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'torque-tension coupling gate'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI GAs GJ Irt Irn Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (N) (Nm2) (kgm) (kgm) (-) (-) (-) (-)'
    WRITE (u, '(A,ES24.16,A)') 'cab 0.2 ', RHOW*0.25_wp*PI*0.04_wp, ' 1.0e7 0.0 1.0e6 1.0e8 '//TRIM(gjcol)// &
      ' 1.0 1.0 0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Coupled 0.0 0.0 -50.0 0 0 0 0'
    WRITE (u, '(A)') '2 Fixed 100.0 0.0 -50.0 0 0 0 0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 cab 99.0 40'
    WRITE (u, '(A)') '--- END CONNECTIONS ---'
    WRITE (u, '(A)') 'LineID End Stiffness EzX EzY EzZ TorsStiffness NxX NxY NxZ Pretwist'
    WRITE (u, '(A)') '(-) (-) (N-m/rad) (-) (-) (-) (N-m/rad) (-) (-) (-) (deg)'
    WRITE (u, '(A)') '1 A '//conn_a
    WRITE (u, '(A)') '1 B '//conn_b
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    IF (dyn) THEN
      WRITE (u, '(A)') '0.05 dtM'
      WRITE (u, '(A)') '25.0 TMax'
      WRITE (u, '(A)') 'tcpl_motion.txt motionFile'
    ELSE
      WRITE (u, '(A)') '0.1 dtM'
      WRITE (u, '(A)') '0.0 TMax'
    END IF
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1'
    IF (tors) THEN
      WRITE (u, '(A)') 'Torq1N1'
      WRITE (u, '(A)') 'Twist1'
    END IF
    WRITE (u, '(A)') 'Ten1N21'
    WRITE (u, '(A)') '--- END ---'
    CLOSE (u)
  END SUBROUTINE write_case

  PURE FUNCTION remove_blanks(s) RESULT(r)
    CHARACTER(*), INTENT(IN) :: s
    CHARACTER(LEN(s)) :: r
    INTEGER :: i, k
    r = ''
    k = 0
    DO i = 1, LEN(s)
      IF (s(i:i) /= ' ') THEN
        k = k + 1
        r(k:k) = s(i:i)
      END IF
    END DO
  END FUNCTION remove_blanks

  PURE FUNCTION cross(a, b) RESULT(c)
    REAL(wp), INTENT(IN) :: a(3), b(3)
    REAL(wp) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  END FUNCTION cross

  PURE FUNCTION unit(a) RESULT(b)
    REAL(wp), INTENT(IN) :: a(3)
    REAL(wp) :: b(3)
    b = a/NORM2(a)
  END FUNCTION unit

  SUBROUTINE require(ok, what)
    LOGICAL, INTENT(IN) :: ok
    CHARACTER(*), INTENT(IN) :: what
    IF (ok) THEN
      PRINT '(A)', 'PASS: '//what
    ELSE
      PRINT '(A)', 'FAIL: '//what
      nfail = nfail + 1
    END IF
  END SUBROUTINE require
END PROGRAM test_torsion_coupling
