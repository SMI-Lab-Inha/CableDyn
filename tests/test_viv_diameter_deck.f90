! File: tests/test_viv_diameter_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_diameter_deck
  !! The diameter of the VIV model on a deck (the optional VIV DIAMETER section, rows TypeName
  !! Dviv): the taut gate line (100 m, Diam 0.1 m, 50 kN, pinned at 50 m depth, 50 elements) in a
  !! uniform 0.15 m/s current along y, as a finite-EI cable
  !! (EI 1 N m^2) and as an EI = 0 line:
  !!   D1 Dviv equal to the type's Diam is the deck without the section, to the printed digits;
  !!   D2 with the synchronisation lift (Cl 0.8, cF 0.18, dF 0.08) the deck locks in near the first
  !!      mode (cF U/Diam = 0.27 Hz, the lock-in at 0.287 Hz); with Dviv 0.2 m the shedding band
  !!      (cF -+ dF) U/Dviv, 0.075-0.195 Hz, lies below that mode, where the motion lags the lift by
  !!      less than a quarter period and the phase runs at the top of the band, 0.195 Hz (within
  !!      0.005 Hz), on both line kinds;
  !!   D3 a wake oscillator uncoupled from the cable (A 0, St 0.2) drives the finite-EI cable at
  !!      St U/Dviv = 0.15 Hz with Dviv 0.2 m (within 0.005 Hz), half the 0.30 Hz of Diam;
  !!   D4 refusals by name: a malformed row, an unknown or repeated type, a value out of range, and a
  !!      type without VIV (a row with no number is a header row, as in every table section).
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, D = 0.1_wp
  INTEGER :: nfail = 0

  CALL write_current('vdd_current.dat')
  CALL check_kind('1.0', 'finite-EI')
  CALL check_kind('0.0', 'EI = 0')
  CALL check_wake()
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: VIV diameter deck gates'

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
    WRITE (u, '(A)') 'uniform current along y for the VIV diameter deck gate'
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

  SUBROUTINE write_deck(path, dia_rows, tmax, ei, cl, wake_row)
    !! dia_rows: the VIV DIAMETER section's data rows, ';'-separated ('' for no section); ei the
    !! type's EI, cl its LINE TYPES Cl (default 0.8); wake_row a WAKE OSCILLATOR row.
    CHARACTER(*), INTENT(IN) :: path, dia_rows, ei
    REAL(wp), INTENT(IN) :: tmax
    CHARACTER(*), INTENT(IN), OPTIONAL :: cl, wake_row
    CHARACTER(32) :: cl_col
    INTEGER :: u, i0, i1
    cl_col = '0.8'
    IF (PRESENT(cl)) cl_col = cl
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (u, '(A)') 'VIV diameter deck gate: taut cable in a uniform current'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (u, '(A)') 'TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF'
    WRITE (u, '(A)') 'riser      0.10   10.0     5.0e7     -0.8      '//TRIM(ei)//'  1.2   1.0   0.0   0.0   '// &
      TRIM(cl_col)//'  0.08  0.18'
    WRITE (u, '(A)') 'bare       0.10   10.0     5.0e7     -0.8      '//TRIM(ei)//'  1.2   1.0   0.0   0.0   0.0'
    WRITE (u, '(A)') '----------------------- POINTS ----------------------------------------------'
    WRITE (u, '(A)') 'ID   Attachment   X       Y     Z      M    V     CdA   CA'
    WRITE (u, '(A)') '1    Coupled      0.0     0.0   -50.0  0    0     0     0'
    WRITE (u, '(A)') '2    Fixed        100.0   0.0   -50.0  0    0     0     0'
    WRITE (u, '(A)') '----------------------- LINES -----------------------------------------------'
    WRITE (u, '(A)') 'ID   NodeA  NodeB  Outputs'
    WRITE (u, '(A)') '1    1      2      -'
    WRITE (u, '(A)') '----------------------- SECTIONS --------------------------------------------'
    WRITE (u, '(A)') 'LineID  LineType  Length     NumSegs'
    WRITE (u, '(A)') '1       riser     99.900100  50'
    IF (PRESENT(wake_row)) THEN
      WRITE (u, '(A)') '----------------------- WAKE OSCILLATOR -------------------------------------'
      WRITE (u, '(A)') 'TypeName   CL0    St     eps    A'
      WRITE (u, '(A)') wake_row
    END IF
    IF (LEN_TRIM(dia_rows) > 0) THEN
      WRITE (u, '(A)') '----------------------- VIV DIAMETER ----------------------------------------'
      WRITE (u, '(A)') 'TypeName   Dviv'
      WRITE (u, '(A)') '(-)        (m)'
      i0 = 1
      DO
        i1 = INDEX(dia_rows(i0:), ';')
        IF (i1 == 0) THEN
          WRITE (u, '(A)') TRIM(dia_rows(i0:))
          EXIT
        END IF
        WRITE (u, '(A)') dia_rows(i0:i0 + i1 - 2)
        i0 = i0 + i1
      END DO
    END IF
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '0.01     dtM'
    WRITE (u, '(A)') '100.0    WtrDpth'
    WRITE (u, '(A)') '1025.0   rhoW'
    WRITE (u, '(A)') '9.80665  g'
    WRITE (u, '(A)') 'vdd_current.dat  WaterKin'
    WRITE (u, '(A,F8.3,A)') ' ', tmax, '    TMax'
    WRITE (u, '(A)') '---------------------- OUTPUTS -----------------------------------------------'
    WRITE (u, '(A)') 'L1N26PZ'
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
    !! The rows (time, L1N26PZ) of a deck's main output file.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: tab(:, :)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(1024) :: rec
    INTEGER :: u, ios
    ALLOCATE (tab(2, 20000))
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
    !! The frequency of the largest Fourier amplitude of x - mean(x) on 0.05-1 Hz (0.001 Hz grid).
    REAL(wp), INTENT(IN) :: tt(:), x(:)
    REAL(wp) :: xm(SIZE(x)), f, best, re, im
    INTEGER :: j
    xm = x - SUM(x)/REAL(SIZE(x), wp)
    best = -1.0_wp
    fpk = 0.0_wp
    DO j = 0, 950
      f = 0.05_wp + 0.001_wp*REAL(j, wp)
      re = SUM(xm*COS(2.0_wp*PI*f*tt))
      im = SUM(xm*SIN(2.0_wp*PI*f*tt))
      IF (re*re + im*im > best) THEN
        best = re*re + im*im
        fpk = f
      END IF
    END DO
  END FUNCTION peak_frequency

  SUBROUTINE window_stats(tab, n, t0, ad, f)
    !! Half peak-to-peak mid-span cross-flow amplitude over D and dominant frequency from t0 on.
    REAL(wp), INTENT(IN) :: tab(:, :), t0
    INTEGER, INTENT(IN) :: n
    REAL(wp), INTENT(OUT) :: ad, f
    LOGICAL :: keep(n)
    REAL(wp), ALLOCATABLE :: tt(:), z(:)
    keep = tab(1, 1:n) >= t0
    tt = PACK(tab(1, 1:n), keep)
    z = PACK(tab(2, 1:n), keep)
    ad = 0.5_wp*(MAXVAL(z) - MINVAL(z))/D
    f = peak_frequency(tt, z)
  END SUBROUTINE window_stats

  SUBROUTINE check_kind(ei, label)
    CHARACTER(*), INTENT(IN) :: ei, label
    CHARACTER(16) :: tag
    REAL(wp), ALLOCATABLE :: tn(:, :), ts(:, :), tb(:, :)
    INTEGER :: nn, ns, nb
    REAL(wp) :: ad_n, f_n, ad_b, f_b
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    tag = MERGE('fe', 'e0', ei /= '0.0')
    CALL write_deck('vdd_'//TRIM(tag)//'_none.dat', '', 100.0_wp, ei)
    CALL run('vdd_'//TRIM(tag)//'_none.dat', 'vdd_'//TRIM(tag)//'_none', ok, msg)
    CALL require(ok, 'D1 ('//label//'): the deck without the section runs: '//TRIM(msg))
    CALL write_deck('vdd_'//TRIM(tag)//'_same.dat', 'riser 0.10', 100.0_wp, ei)
    CALL run('vdd_'//TRIM(tag)//'_same.dat', 'vdd_'//TRIM(tag)//'_same', ok, msg)
    CALL require(ok, 'D1 ('//label//'): the Dviv = Diam deck runs: '//TRIM(msg))
    CALL write_deck('vdd_'//TRIM(tag)//'_big.dat', 'riser 0.2', 100.0_wp, ei)
    CALL run('vdd_'//TRIM(tag)//'_big.dat', 'vdd_'//TRIM(tag)//'_big', ok, msg)
    CALL require(ok, 'D2 ('//label//'): the Dviv 0.2 deck runs: '//TRIM(msg))
    CALL read_out('vdd_'//TRIM(tag)//'_none.out', tn, nn)
    CALL read_out('vdd_'//TRIM(tag)//'_same.out', ts, ns)
    CALL read_out('vdd_'//TRIM(tag)//'_big.out', tb, nb)
    CALL require(nn > 5000 .AND. ns == nn .AND. nb == nn, 'D1 ('//label//'): the runs wrote matching rows')
    IF (.NOT. (nn > 5000 .AND. ns == nn .AND. nb == nn)) RETURN
    CALL require(ALL(ts(:, 1:ns) == tn(:, 1:nn)), 'D1 ('//label//'): Dviv = Diam is the deck without the section')
    CALL window_stats(tn, nn, 60.0_wp, ad_n, f_n)
    CALL window_stats(tb, nb, 60.0_wp, ad_b, f_b)
    WRITE (*, '(A,2F8.4,A,2F8.4)') 'D2 ('//label//'): A/D, f [Hz] with Diam ', ad_n, f_n, '; with Dviv 0.2 ', &
      ad_b, f_b
    CALL require(ABS(f_n - 0.28_wp) <= 0.02_wp, 'D2 ('//label//'): with Diam the line locks in near its first mode')
    CALL require(ABS(f_b - 0.195_wp) <= 0.005_wp, 'D2 ('//label//'): with Dviv 0.2 the phase runs at the top of '// &
                 'its band, (cF + dF) U/Dviv')
  END SUBROUTINE check_kind

  SUBROUTINE check_wake()
    REAL(wp), ALLOCATABLE :: tw(:, :)
    INTEGER :: nw
    REAL(wp) :: ad, f
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('vdd_wake.dat', 'riser 0.2', 100.0_wp, '1.0', cl='0.0', wake_row='riser 0.3 0.2 0.3 0.0')
    CALL run('vdd_wake.dat', 'vdd_wake', ok, msg)
    CALL require(ok, 'D3: the uncoupled wake-oscillator deck with Dviv 0.2 runs: '//TRIM(msg))
    CALL read_out('vdd_wake.out', tw, nw)
    CALL require(nw > 5000, 'D3: the run wrote its rows')
    IF (nw <= 5000) RETURN
    CALL window_stats(tw, nw, 50.0_wp, ad, f)
    WRITE (*, '(A,2F8.4)') 'D3: uncoupled wake, Dviv 0.2: A/D, f [Hz] ', ad, f
    CALL require(ABS(f - 0.15_wp) <= 0.005_wp, 'D3: the wake drives the cable at St U/Dviv')
  END SUBROUTINE check_wake

  SUBROUTINE check_refusals()
    CALL refused('riser 0.2 0.3', 'TypeName Dviv', 'a row of three columns')
    CALL refused('cable 0.2', 'undefined LINE TYPE', 'an unknown type')
    CALL refused('riser 0.2;riser 0.3', 'twice', 'a repeated type')
    CALL refused('riser 0.0', 'Dviv > 0', 'a zero Dviv')
    CALL refused('riser -0.2', 'Dviv > 0', 'a negative Dviv')
    CALL refused('bare 0.2', 'give the type LINE TYPES Cl > 0 or a WAKE OSCILLATOR row', 'a type without VIV')
  END SUBROUTINE check_refusals

  SUBROUTINE refused(rows, needle, what)
    CHARACTER(*), INTENT(IN) :: rows, needle, what
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('vdd_bad.dat', rows, 1.0_wp, '1.0')
    CALL run('vdd_bad.dat', 'vdd_bad', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, needle) > 0, 'D4: '//what//' is refused: '//TRIM(msg))
  END SUBROUTINE refused

END PROGRAM test_viv_diameter_deck
