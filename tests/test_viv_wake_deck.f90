! File: tests/test_viv_wake_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_wake_deck
  !! Wake-oscillator VIV on a deck (the optional WAKE OSCILLATOR section, TypeName CL0 [St [eps
  !! [A]]]; Facchinetti et al. 2004): the taut finite-EI cable of the in-line VIV gate (100 m,
  !! D 0.1 m, EI 1 N m^2, 50 kN, pinned at 50 m depth) in a uniform 0.15 m/s current along y, its
  !! type with no synchronisation lift (Cl 0) and the row riser 0.3 (St 0.2, eps 0.3, A 12).
  !!   W1 lock-in: with the acceleration coupling A 12 the cross-flow response is captured by the
  !!      cable's first mode, below the wake's own frequency St U/D = 0.30 Hz by more than
  !!      0.01 Hz, at a steady mid-span A/D between 0.4 and 0.8; uncoupled (A 0) the wake drives
  !!      the cable at its own 0.30 Hz (within 0.005 Hz), and the coupled amplitude is at least
  !!      three times the uncoupled one;
  !!   W2 power: over 60-100 s the mean wake-lift power (VIVPow1) is positive and the drag and
  !!      damping powers (DragPow1, DampPow1) take it out within 0.5 %;
  !!   W3 a row with CL0 0 is the deck without the section, to the printed digits;
  !!   W4 refusals by name: a malformed or non-numeric row, an unknown or repeated type, St, eps or
  !!      A out of range, a type that keeps its synchronisation lift (Cl > 0, which in-line VIV
  !!      also needs), a line that mixes the two models, and an EI = 0 line.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, D = 0.1_wp, F_WAKE = 0.2_wp*0.15_wp/0.1_wp
  INTEGER :: nfail = 0

  CALL write_current('wod_current.dat')
  CALL check_runs()
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: wake-oscillator deck gates'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE write_current(path)
    CHARACTER(*), INTENT(IN) :: path
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'MoorDyn v2 water kinematics file'
    WRITE (u, '(A)') 'uniform current along y for the wake-oscillator deck gate'
    WRITE (u, '(A)') '--------------------------- WAVES -------------------------------------'
    WRITE (u, '(A)') '0                    WaveKinMod  - type of wave input'
    WRITE (u, '(A)') '""                   WaveKinFile - file containing wave elevation time series'
    WRITE (u, '(A)') '0.000000E+00  dtWave      - time step to use in setting up wave kinematics grid (s)'
    WRITE (u, '(A)') '0                    WaveDir     - wave heading (deg)'
    WRITE (u, '(A)') '2                                - X wave input type'
    WRITE (u, '(A)') '-10, 110, 2                      - X wave grid point data'
    WRITE (u, '(A)') '2                                - Y wave input type'
    WRITE (u, '(A)') '-10, 10, 2                       - Y wave grid point data'
    WRITE (u, '(A)') '2                                - Z wave input type'
    WRITE (u, '(A)') '-100, 0, 2                       - Z wave grid point data'
    WRITE (u, '(A)') '--------------------------- CURRENT -------------------------------------'
    WRITE (u, '(A)') '1                    CurrentMod  - type of current input'
    WRITE (u, '(A)') 'z-depth     x-current      y-current'
    WRITE (u, '(A)') '(m)           (m/s)         (m/s)'
    WRITE (u, '(A)') '0.0    0.0  0.15'
    WRITE (u, '(A)') '-100.0 0.0  0.15'
    WRITE (u, '(A)') '--------------------- need this line ------------------'
    CLOSE (u)
  END SUBROUTINE write_current

  SUBROUTINE write_deck(path, wake_rows, tmax, cl, ei, viv_row, mixed)
    !! wake_rows: the WAKE OSCILLATOR section's data rows, ';'-separated ('' for no section); cl the
    !! type's LINE TYPES Cl (default 0), ei its EI (default 1), viv_row a VIV section row, mixed:
    !! the line in two sections, the second of a type 'lift' with Cl 0.8.
    CHARACTER(*), INTENT(IN) :: path, wake_rows
    REAL(wp), INTENT(IN) :: tmax
    CHARACTER(*), INTENT(IN), OPTIONAL :: cl, ei, viv_row
    LOGICAL, INTENT(IN), OPTIONAL :: mixed
    CHARACTER(32) :: cl_col, ei_col
    INTEGER :: u, i0, i1
    LOGICAL :: two
    cl_col = '0.0'
    IF (PRESENT(cl)) cl_col = cl
    ei_col = '1.0'
    IF (PRESENT(ei)) ei_col = ei
    two = .FALSE.
    IF (PRESENT(mixed)) two = mixed
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (u, '(A)') 'wake-oscillator deck gate: taut cable in a uniform current'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (u, '(A)') 'TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF'
    WRITE (u, '(A)') 'riser      0.10   10.0     5.0e7     -0.8      '//TRIM(ei_col)//'  1.2   1.0   0.0   0.0   '// &
      TRIM(cl_col)//'  0.08  0.18'
    WRITE (u, '(A)') 'lift       0.10   10.0     5.0e7     -0.8      '//TRIM(ei_col)//'  1.2   1.0   0.0   0.0   '// &
      '0.8  0.08  0.18'
    WRITE (u, '(A)') '----------------------- POINTS ----------------------------------------------'
    WRITE (u, '(A)') 'ID   Attachment   X       Y     Z      M    V     CdA   CA'
    WRITE (u, '(A)') '1    Coupled      0.0     0.0   -50.0  0    0     0     0'
    WRITE (u, '(A)') '2    Fixed        100.0   0.0   -50.0  0    0     0     0'
    WRITE (u, '(A)') '----------------------- LINES -----------------------------------------------'
    WRITE (u, '(A)') 'ID   NodeA  NodeB  Outputs'
    WRITE (u, '(A)') '1    1      2      -'
    WRITE (u, '(A)') '----------------------- SECTIONS --------------------------------------------'
    WRITE (u, '(A)') 'LineID  LineType  Length     NumSegs'
    IF (two) THEN
      WRITE (u, '(A)') '1       riser     49.950050  25'
      WRITE (u, '(A)') '1       lift      49.950050  25'
    ELSE
      WRITE (u, '(A)') '1       riser     99.900100  50'
    END IF
    IF (PRESENT(viv_row)) THEN
      WRITE (u, '(A)') '----------------------- VIV -------------------------------------------------'
      WRITE (u, '(A)') 'TypeName   ClX    alpha'
      WRITE (u, '(A)') viv_row
    END IF
    IF (LEN_TRIM(wake_rows) > 0) THEN
      WRITE (u, '(A)') '----------------------- WAKE OSCILLATOR -------------------------------------'
      WRITE (u, '(A)') 'TypeName   CL0    St     eps    A'
      i0 = 1
      DO
        i1 = INDEX(wake_rows(i0:), ';')
        IF (i1 == 0) THEN
          WRITE (u, '(A)') TRIM(wake_rows(i0:))
          EXIT
        END IF
        WRITE (u, '(A)') wake_rows(i0:i0 + i1 - 2)
        i0 = i0 + i1
      END DO
    END IF
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '0.01     dtM'
    WRITE (u, '(A)') '100.0    WtrDpth'
    WRITE (u, '(A)') '1025.0   rhoW'
    WRITE (u, '(A)') '9.80665  g'
    WRITE (u, '(A)') 'wod_current.dat  WaterKin'
    WRITE (u, '(A,F8.3,A)') ' ', tmax, '    TMax'
    WRITE (u, '(A)') '---------------------- OUTPUTS -----------------------------------------------'
    WRITE (u, '(A)') 'L1N25PZ'
    WRITE (u, '(A)') 'VIVPow1'
    WRITE (u, '(A)') 'DragPow1'
    WRITE (u, '(A)') 'DampPow1'
    WRITE (u, '(A)') 'END'
    WRITE (u, '(A)') '------------------------- need this line --------------------------------------'
    CLOSE (u)
  END SUBROUTINE write_deck

  SUBROUTINE run(deck, root, ok, msg)
    CHARACTER(*), INTENT(IN) :: deck, root
    LOGICAL, INTENT(OUT) :: ok
    CHARACTER(*), INTENT(OUT) :: msg
    LOGICAL :: conv
    INTEGER :: es
    CALL CD_Run_Deck_Driver(deck, root, conv, es, msg)
    ok = es == CD_DECKDRV_OK .AND. conv
  END SUBROUTINE run

  SUBROUTINE read_out(path, tab, n)
    !! The rows (time, L1N25PZ, VIVPow1, DragPow1, DampPow1) of a deck's main output file.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: tab(:, :)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(1024) :: rec
    INTEGER :: u, ios
    ALLOCATE (tab(5, 20000))
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

  REAL(wp) FUNCTION peak_frequency(tt, x) RESULT(fpk)
    !! The frequency of the largest Fourier amplitude of x - mean(x) on 0.05-2 Hz (0.001 Hz grid).
    REAL(wp), INTENT(IN) :: tt(:), x(:)
    REAL(wp) :: xm(SIZE(x)), f, best, re, im
    INTEGER :: j
    xm = x - SUM(x)/REAL(SIZE(x), wp)
    best = -1.0_wp
    fpk = 0.0_wp
    DO j = 0, 1950
      f = 0.05_wp + 0.001_wp*REAL(j, wp)
      re = SUM(xm*COS(2.0_wp*PI*f*tt))
      im = SUM(xm*SIN(2.0_wp*PI*f*tt))
      IF (re*re + im*im > best) THEN
        best = re*re + im*im
        fpk = f
      END IF
    END DO
  END FUNCTION peak_frequency

  SUBROUTINE window_stats(tab, n, t0, t1, ad, f)
    !! Half peak-to-peak mid-span amplitude over D and dominant frequency on [t0, t1].
    REAL(wp), INTENT(IN) :: tab(:, :), t0, t1
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(OUT) :: ad, f
    LOGICAL :: keep(n)
    REAL(wp), ALLOCATABLE :: tt(:), z(:)
    keep = tab(1, 1:n) >= t0 .AND. tab(1, 1:n) <= t1
    tt = PACK(tab(1, 1:n), keep)
    z = PACK(tab(2, 1:n), keep)
    ad = 0.5_wp*(MAXVAL(z) - MINVAL(z))/D
    f = peak_frequency(tt, z)
  END SUBROUTINE window_stats

  SUBROUTINE check_runs()
    REAL(wp), ALLOCATABLE :: tw(:, :), t0(:, :), tn(:, :), tz(:, :)
    INTEGER :: nw, n0, nn, nz
    REAL(wp) :: ad_w, f_w, ad_w2, f_w2, ad_0, f_0, pv, pd, pb
    LOGICAL :: ok, keep(20000)
    CHARACTER(1024) :: msg
    CALL write_deck('wod_wake.dat', 'riser 0.3', 100.0_wp)
    CALL run('wod_wake.dat', 'wod_wake', ok, msg)
    CALL require(ok, 'W1: the wake-oscillator deck runs: '//TRIM(msg))
    CALL write_deck('wod_a0.dat', 'riser 0.3 0.2 0.3 0.0', 100.0_wp)
    CALL run('wod_a0.dat', 'wod_a0', ok, msg)
    CALL require(ok, 'W1: the uncoupled deck (A 0) runs: '//TRIM(msg))
    CALL read_out('wod_wake.out', tw, nw)
    CALL read_out('wod_a0.out', t0, n0)
    CALL require(nw > 1000 .AND. n0 == nw, 'W1: the runs wrote matching rows')
    IF (.NOT. (nw > 1000 .AND. n0 == nw)) RETURN
    CALL window_stats(tw, nw, 50.0_wp, 75.0_wp, ad_w, f_w)
    CALL window_stats(tw, nw, 75.0_wp, 100.0_wp, ad_w2, f_w2)
    CALL window_stats(t0, n0, 50.0_wp, 100.0_wp, ad_0, f_0)
    WRITE (*, '(A,2F8.4,A,2F8.4,A,F7.4)') 'W1: coupled A/D, f [Hz] (50-75 s, 75-100 s) ', ad_w, f_w, ', ', ad_w2, &
      f_w2, '; wake frequency ', F_WAKE
    WRITE (*, '(A,2F8.4)') 'W1: uncoupled (A 0) A/D, f [Hz] ', ad_0, f_0
    CALL require(f_w < F_WAKE - 0.01_wp .AND. ad_w > 0.4_wp .AND. ad_w < 0.8_wp, &
                 'W1: the coupled response locks in below the wake frequency at A/D 0.4-0.8')
    CALL require(ABS(ad_w2/ad_w - 1.0_wp) <= 0.01_wp .AND. ABS(f_w2 - f_w) <= 0.002_wp, 'W1: and is steady')
    CALL require(ABS(f_0 - F_WAKE) <= 0.005_wp, 'W1: uncoupled, the wake drives the cable at its own frequency')
    CALL require(ad_w >= 3.0_wp*ad_0, 'W1: lock-in at least triples the uncoupled amplitude')
    ! W2: the power balance of the steady lock-in
    keep(1:nw) = tw(1, 1:nw) >= 60.0_wp
    pv = SUM(PACK(tw(3, 1:nw), keep(1:nw)))/COUNT(keep(1:nw))
    pd = SUM(PACK(tw(4, 1:nw), keep(1:nw)))/COUNT(keep(1:nw))
    pb = SUM(PACK(tw(5, 1:nw), keep(1:nw)))/COUNT(keep(1:nw))
    WRITE (*, '(A,3ES12.4)') 'W2: mean VIVPow, DragPow, DampPow [W] over 60-100 s: ', pv, pd, pb
    CALL require(pv > 0.0_wp .AND. ABS(pv + pd + pb) <= 0.005_wp*pv, &
                 'W2: the drag and damping take out the wake-lift power within 0.5 %')
    ! W3: CL0 0 is no wake oscillator
    CALL write_deck('wod_none.dat', '', 20.0_wp)
    CALL run('wod_none.dat', 'wod_none', ok, msg)
    CALL require(ok, 'W3: the deck without the section runs: '//TRIM(msg))
    CALL write_deck('wod_zero.dat', 'riser 0.0', 20.0_wp)
    CALL run('wod_zero.dat', 'wod_zero', ok, msg)
    CALL require(ok, 'W3: the CL0 0 deck runs: '//TRIM(msg))
    CALL read_out('wod_none.out', tn, nn)
    CALL read_out('wod_zero.out', tz, nz)
    CALL require(nn > 100 .AND. nz == nn, 'W3: the runs wrote matching rows')
    IF (nn > 100 .AND. nz == nn) CALL require(ALL(tz(:, 1:nz) == tn(:, 1:nn)), &
                                              'W3: a row with CL0 0 is the deck without the section')
  END SUBROUTINE check_runs

  SUBROUTINE check_refusals()
    CALL refused('riser 0.3 0.2 0.3 12 1', 'TypeName CL0 [St [eps [A]]]', 'a row of six columns')
    CALL refused('riser 0.3 x', 'plain numbers', 'a non-numeric St')
    CALL refused('cable 0.3', 'undefined LINE TYPE', 'an unknown type')
    CALL refused('riser 0.3;riser 0.2', 'twice', 'a repeated type')
    CALL refused('riser -0.1', 'CL0 must be non-negative', 'a negative CL0')
    CALL refused('riser 0.3 0.0', 'St > 0, eps > 0 and A >= 0', 'St = 0')
    CALL refused('riser 0.3 0.2 0.0', 'St > 0, eps > 0 and A >= 0', 'eps = 0')
    CALL refused('riser 0.3 0.2 0.3 -1', 'St > 0, eps > 0 and A >= 0', 'a negative A')
    CALL refused('riser 0.3', 'give the type LINE TYPES Cl 0', 'a type that keeps its lift', cl='0.8')
    CALL refused('riser 0.3', 'a line takes one VIV model', 'a line mixing the two models', mixed=.TRUE.)
    CALL refused('riser 0.3', 'finite-EI cables (EI > 0) only', 'an EI = 0 line', ei='0.0')
  END SUBROUTINE check_refusals

  SUBROUTINE refused(rows, needle, what, cl, ei, viv_row, mixed)
    CHARACTER(*), INTENT(IN) :: rows, needle, what
    CHARACTER(*), INTENT(IN), OPTIONAL :: cl, ei, viv_row
    LOGICAL, INTENT(IN), OPTIONAL :: mixed
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    IF (PRESENT(viv_row)) THEN
      CALL write_deck('wod_bad.dat', rows, 1.0_wp, cl=cl, viv_row=viv_row)
    ELSE IF (PRESENT(cl)) THEN
      CALL write_deck('wod_bad.dat', rows, 1.0_wp, cl=cl)
    ELSE IF (PRESENT(ei)) THEN
      CALL write_deck('wod_bad.dat', rows, 1.0_wp, ei=ei)
    ELSE IF (PRESENT(mixed)) THEN
      CALL write_deck('wod_bad.dat', rows, 1.0_wp, mixed=mixed)
    ELSE
      CALL write_deck('wod_bad.dat', rows, 1.0_wp)
    END IF
    CALL run('wod_bad.dat', 'wod_bad', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, needle) > 0, 'W4: '//what//' is refused: '//TRIM(msg))
  END SUBROUTINE refused

END PROGRAM test_viv_wake_deck
