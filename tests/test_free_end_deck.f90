! File: tests/test_free_end_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_free_end_deck
  !! A finite-EI line ending on a Free point (a free end): a 7 m cable (Diam 0.048 m, 2.8 kg/m,
  !! EA 1e8 N) hanging from a Coupled point at z = -0.5 m with a 200 kg clump (0.02548 m^3, CdA
  !! 0.06 m^2, CA 1) on its Free End B, 70 elements:
  !!   F1 in still water the run holds the static state: the top tension is the submerged weight of
  !!      the line and the clump (to 1e-6 relative), the end hangs one stretched line length below
  !!      the top (the stretch (W_B L + w L^2/2)/EA, to 1e-7 m), and nothing moves over the run;
  !!   F2 in a 0.5 m/s current along x the end rests where the same deck as an EI = 0 line on the
  !!      point-system route puts its Free point (EI 3 N m^2: within 2e-5 m), and the Point
  !!      channels of the free end read its line's End B node;
  !!   F3 a prescribed 0.05 m sway of the top at 0.05 Hz moves the end in phase with an amplitude
  !!      ratio within 1 % of the Rayleigh estimate for a hanging string with a heavy end mass
  !!      (1.104: the end mass with its added mass, a third of the line's, the mean tension);
  !!   F4 a Free point that is not such an end leaves the deck as it was, refused by the routes'
  !!      own messages: a clamped End B, a net-buoyant end, a second Free point.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, G = 9.80665_wp, RHO = 1025.0_wp
  REAL(wp), PARAMETER :: L = 7.0_wp, D = 0.048_wp, M = 2.8_wp, EA = 1.0e8_wp
  REAL(wp), PARAMETER :: MB = 200.0_wp, VB = 0.02548_wp
  INTEGER :: nfail = 0

  CALL check_still()
  CALL check_current()
  CALL check_sway()
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: free cable end deck gates'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE write_profile(path, u)
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), INTENT(IN) :: u
    INTEGER :: k
    OPEN (NEWUNIT=k, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (k, '(A)') 'z ux uy uz'
    WRITE (k, '(A)') '(m) (m/s) (m/s) (m/s)'
    WRITE (k, '(A,F6.3,A)') '-20 ', u, ' 0 0'
    WRITE (k, '(A,F6.3,A)') '0 ', u, ' 0 0'
    CLOSE (k)
  END SUBROUTINE write_profile

  SUBROUTINE write_deck(path, ei, pt2, extra_pts, endconn, opts_extra, outs, dt, tmax)
    !! ei: the type's EI; pt2: the Free point row of End B (after 'ID'); extra_pts: more POINTS rows
    !! (';'-separated, '' for none); endconn: an END CONNECTIONS data row ('' for none);
    !! opts_extra: more OPTIONS rows (';'-separated); outs: OUTPUTS channels (';'-separated).
    CHARACTER(*), INTENT(IN) :: path, ei, pt2, extra_pts, endconn, opts_extra, outs
    REAL(wp), INTENT(IN) :: dt, tmax
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (u, '(A)') 'free cable end gate'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (u, '(A)') 'TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'cable 0.048 2.8 1.0e8 0.0 '//ei//' 1.2 0.0 1.0 0.0'
    WRITE (u, '(A)') '----------------------- POINTS ----------------------------------------------'
    WRITE (u, '(A)') 'ID Attachment X Y Z M V CdA CA'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)'
    WRITE (u, '(A)') '1 Coupled 0.0 0.0 -0.5 0 0 0 0'
    WRITE (u, '(A)') '2 '//pt2
    CALL put_rows(u, extra_pts)
    WRITE (u, '(A)') '----------------------- LINES -----------------------------------------------'
    WRITE (u, '(A)') 'ID LineType AttachA AttachB UnstrLen NumSegs Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-) (m) (-) (-)'
    WRITE (u, '(A)') '1 cable 1 2 7.0 70 -'
    IF (LEN_TRIM(endconn) > 0) THEN
      WRITE (u, '(A)') '--------------------- END CONNECTIONS ----------------------------------'
      WRITE (u, '(A)') 'LineID End Stiffness EzX EzY EzZ'
      WRITE (u, '(A)') '(-) (-) (N-m/rad) (-) (-) (-)'
      WRITE (u, '(A)') endconn
    END IF
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '1025 rhoW'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '20 WtrDpth'
    WRITE (u, '(ES12.5,A)') dt, ' dtM'
    WRITE (u, '(ES12.5,A)') tmax, ' TMax'
    CALL put_rows(u, opts_extra)
    WRITE (u, '(A)') '---------------------- OUTPUTS -----------------------------------------------'
    CALL put_rows(u, outs)
    WRITE (u, '(A)') 'END'
    WRITE (u, '(A)') '------------------------- need this line --------------------------------------'
    CLOSE (u)
  END SUBROUTINE write_deck

  SUBROUTINE put_rows(u, rows)
    INTEGER, INTENT(IN) :: u
    CHARACTER(*), INTENT(IN) :: rows
    INTEGER :: i0, i1
    i0 = 1
    DO WHILE (i0 <= LEN_TRIM(rows))
      i1 = INDEX(rows(i0:), ';')
      IF (i1 == 0) THEN
        WRITE (u, '(A)') TRIM(rows(i0:))
        EXIT
      END IF
      WRITE (u, '(A)') rows(i0:i0 + i1 - 2)
      i0 = i0 + i1
    END DO
  END SUBROUTINE put_rows

  SUBROUTINE run(deck, root, ok, msg)
    CHARACTER(*), INTENT(IN) :: deck, root
    LOGICAL, INTENT(OUT) :: ok
    CHARACTER(*), INTENT(OUT) :: msg
    LOGICAL :: conv
    INTEGER :: es
    CALL CD_Run_Deck_Driver(deck, root, conv, es, msg)
    ok = es == CD_DECKDRV_OK
  END SUBROUTINE run

  SUBROUTINE read_out(path, ncol, tab, n)
    CHARACTER(*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: ncol
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: tab(:, :)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(2048) :: rec
    INTEGER :: u, ios
    ALLOCATE (tab(ncol, 20000))
    n = 0
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    READ (u, '(A)') rec
    READ (u, '(A)') rec
    DO
      READ (u, '(A)', IOSTAT=ios) rec
      IF (ios /= 0) EXIT
      IF (LEN_TRIM(rec) == 0) CYCLE
      IF (n == SIZE(tab, 2)) EXIT
      n = n + 1
      READ (rec, *) tab(:, n)
    END DO
    CLOSE (u)
  END SUBROUTINE read_out

  SUBROUTINE check_still()
    REAL(wp), ALLOCATABLE :: tab(:, :)
    REAL(wp) :: w, wb, stretch
    INTEGER :: n
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('fed_still.dat', '3000', 'Free 0.0 0.0 -7.5 200 0.02548 0.06 1.0', '', '', '', &
                    'FairTen1;Point2px;Point2pz', 0.005_wp, 2.0_wp)
    CALL run('fed_still.dat', 'fed_still', ok, msg)
    CALL require(ok, 'F1 still-water run: '//TRIM(msg))
    IF (.NOT. ok) RETURN
    CALL read_out('fed_still.out', 4, tab, n)
    CALL require(n > 10, 'F1 output rows')
    IF (n <= 10) RETURN
    w = (M - RHO*0.25_wp*PI*D**2)*G
    wb = (MB - RHO*VB)*G
    stretch = (wb*L + 0.5_wp*w*L**2)/EA
    CALL require(ABS(tab(2, 1)/(w*L + wb) - 1.0_wp) < 1.0e-6_wp, 'F1 top tension = submerged weight')
    CALL require(ABS(tab(4, 1) - (-0.5_wp - L - stretch)) < 1.0e-7_wp, 'F1 end one stretched length below')
    CALL require(MAXVAL(ABS(tab(3:4, 1:n) - SPREAD(tab(3:4, 1), 2, n))) < 1.0e-10_wp, 'F1 nothing moves')
  END SUBROUTINE check_still

  SUBROUTINE check_current()
    REAL(wp), ALLOCATABLE :: tf(:, :), t0(:, :)
    INTEGER :: nf, n0
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_profile('current_profile.txt', 0.5_wp)
    CALL write_deck('fed_cur.dat', '3', 'Free 0.0 0.0 -7.5 200 0.02548 0.06 1.0', '', '', '1 Currents', &
                    'Point2px;Con2pz;L1N71px;L1N71pz', 0.005_wp, 0.5_wp)
    CALL run('fed_cur.dat', 'fed_cur', ok, msg)
    CALL require(ok, 'F2 finite-EI free end in current: '//TRIM(msg))
    CALL write_deck('fed_cur0.dat', '0', 'Free 0.0 0.0 -7.5 200 0.02548 0.06 1.0', '', '', '1 Currents', &
                    'Point2px;Point2pz', 0.005_wp, 0.5_wp)
    CALL run('fed_cur0.dat', 'fed_cur0', ok, msg)
    CALL require(ok, 'F2 EI = 0 Free point in current: '//TRIM(msg))
    CALL read_out('fed_cur.out', 5, tf, nf)
    CALL read_out('fed_cur0.out', 3, t0, n0)
    IF (nf < 1 .OR. n0 < 1) THEN
      CALL require(.FALSE., 'F2 output rows')
      RETURN
    END IF
    CALL require(ABS(tf(2, 1) - t0(2, 1)) < 2.0e-5_wp .AND. ABS(tf(3, 1) - t0(3, 1)) < 2.0e-5_wp, &
                 'F2 rest position = the EI = 0 Free point''s')
    CALL require(tf(2, 1) > 0.1_wp, 'F2 the end is carried downstream')
    CALL require(ABS(tf(2, nf) - tf(4, nf)) < 1.0e-12_wp .AND. ABS(tf(3, nf) - tf(5, nf)) < 1.0e-12_wp, &
                 'F2 Point channels read the End B node')
  END SUBROUTINE check_current

  SUBROUTINE check_sway()
    REAL(wp), PARAMETER :: A = 0.05_wp, F = 0.05_wp, DT = 0.01_wp, TMAX = 120.0_wp, TR = 20.0_wp
    REAL(wp), ALLOCATABLE :: tab(:, :)
    REAL(wp) :: t, r, dr, ddr, om, s, c, ata(3, 3), atb(3), coef(3), w, wb, ma, tmean, om2, ratio
    INTEGER :: u, k, n, i
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    om = 2.0_wp*PI*F
    OPEN (NEWUNIT=u, FILE='fed_motion.txt', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '# time point x y z vx vy vz ax ay az'
    DO k = 0, NINT(TMAX/DT)
      t = REAL(k, wp)*DT
      IF (t < TR) THEN
        r = 0.5_wp*(1.0_wp - COS(PI*t/TR))
        dr = 0.5_wp*PI/TR*SIN(PI*t/TR)
        ddr = 0.5_wp*(PI/TR)**2*COS(PI*t/TR)
      ELSE
        r = 1.0_wp; dr = 0.0_wp; ddr = 0.0_wp
      END IF
      s = SIN(om*t); c = COS(om*t)
      WRITE (u, '(F12.6,A,3ES17.9,A,ES17.9,A,ES17.9,A)') t, ' 1', A*r*s, 0.0_wp, -0.5_wp, ' ', &
        A*(dr*s + r*om*c), ' 0 0 ', A*(ddr*s + 2.0_wp*dr*om*c - r*om*om*s), ' 0 0'
    END DO
    CLOSE (u)
    CALL write_deck('fed_sway.dat', '3000', 'Free 0.0 0.0 -7.5 200 0.02548 0.06 1.0', '', '', &
                    'fed_motion.txt motionFile', 'Point2px', DT, TMAX)
    CALL run('fed_sway.dat', 'fed_sway', ok, msg)
    CALL require(ok, 'F3 prescribed sway: '//TRIM(msg))
    IF (.NOT. ok) RETURN
    CALL read_out('fed_sway.out', 2, tab, n)
    ata = 0.0_wp; atb = 0.0_wp
    DO i = 1, n
      IF (tab(1, i) < 40.0_wp) CYCLE
      coef = [SIN(om*tab(1, i)), COS(om*tab(1, i)), 1.0_wp]
      ata = ata + SPREAD(coef, 2, 3)*SPREAD(coef, 1, 3)
      atb = atb + coef*tab(2, i)
    END DO
    CALL solve3(ata, atb, coef)
    ! Rayleigh estimate, shape s/L from the top: omega_p^2 = (mean tension / L)/(end mass with
    ! added mass + a third of the line's mass with added mass)
    w = (M - RHO*0.25_wp*PI*D**2)*G
    wb = (MB - RHO*VB)*G
    tmean = wb + 0.5_wp*w*L
    ma = (MB + RHO*VB) + (M + RHO*0.25_wp*PI*D**2)*L/3.0_wp
    om2 = (tmean/L)/ma
    ratio = 1.0_wp/(1.0_wp - om**2/om2)
    CALL require(ABS(HYPOT(coef(1), coef(2))/A/ratio - 1.0_wp) < 0.01_wp, 'F3 sway amplitude ratio')
    CALL require(ABS(coef(2)) < 0.05_wp*ABS(coef(1)) .AND. coef(1) > 0.0_wp, 'F3 the end moves in phase')
  END SUBROUTINE check_sway

  SUBROUTINE solve3(a, b, x)
    REAL(wp), INTENT(IN) :: a(3, 3), b(3)
    REAL(wp), INTENT(OUT) :: x(3)
    REAL(wp) :: m(3, 4), f
    INTEGER :: i, j
    m(:, 1:3) = a
    m(:, 4) = b
    DO i = 1, 3
      DO j = i + 1, 3
        f = m(j, i)/m(i, i)
        m(j, :) = m(j, :) - f*m(i, :)
      END DO
    END DO
    DO i = 3, 1, -1
      x(i) = (m(i, 4) - SUM(m(i, i + 1:3)*x(i + 1:3)))/m(i, i)
    END DO
  END SUBROUTINE solve3

  SUBROUTINE check_refusals()
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('fed_r1.dat', '3000', 'Free 0.0 0.0 -7.5 200 0.02548 0.06 1.0', '', '1 B Rigid 0 0 -1', '', &
                    'FairTen1', 0.005_wp, 0.1_wp)
    CALL run('fed_r1.dat', 'fed_r1', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'non-pinned END CONNECTIONS need End B') > 0, &
                 'F4 clamped End B refused: '//TRIM(msg))
    CALL write_deck('fed_r2.dat', '3000', 'Free 0.0 0.0 -7.5 10 0.5 0.06 1.0', '', '', '', 'FairTen1', &
                    0.005_wp, 0.1_wp)
    CALL run('fed_r2.dat', 'fed_r2', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'not supported on the multibody march') > 0, &
                 'F4 net-buoyant end refused: '//TRIM(msg))
    CALL write_deck('fed_r3.dat', '3000', 'Free 0.0 0.0 -7.5 200 0.02548 0.06 1.0', &
                    '3 Free 5.0 0.0 -5.0 10 0 0 0', '', '', 'FairTen1', 0.005_wp, 0.1_wp)
    CALL run('fed_r3.dat', 'fed_r3', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'not supported on the multibody march') > 0, &
                 'F4 a second Free point refused: '//TRIM(msg))
  END SUBROUTINE check_refusals
END PROGRAM test_free_end_deck
