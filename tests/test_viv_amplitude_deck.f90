! File: tests/test_viv_amplitude_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_amplitude_deck
  !! Amplitude-dependent VIV coefficients on a deck (the optional VIV AMPLITUDE section, rows
  !! TypeName AD Cv [CdAmp]; the amplitude estimate of Thorsen et al. 2016, Eq. 8 and 9): the taut
  !! cable of the wake-oscillator gate (100 m, D 0.1 m, 50 kN, pinned at 50 m depth, 50 elements)
  !! with the synchronisation lift Cl 0.8 in a uniform 0.15 m/s current along y, as a finite-EI
  !! cable (EI 1 N m^2) and as an EI = 0 line, each run checked on both line kinds:
  !!   A1 a table that gives Cv = Cl and CdAmp = 1 at every amplitude is the deck without the
  !!      section, to the printed digits;
  !!   A2 the amplitude channel VIVAmp1N26 at mid-span follows the half peak-to-peak cross-flow
  !!      displacement of the steady response there within 1 % (it projects the velocity across
  !!      the ambient flow, as Thorsen et al.: across the relative flow, which turns with the
  !!      response velocity, it read 16 % low at this A/D of 0.68);
  !!   A3 a table whose Cv falls from Cl at A/D 0.1 to zero at A/D AD_ZERO, below the free response
  !!      of the constant-Cl deck, holds the steady response between those two amplitudes;
  !!   A4 a table with Cv 0 and CdAmp 2 at every amplitude (no lift, twice the drag) settles at the
  !!      static in-line mid-span offset of the deck with Cd 2.4, within 0.5 %, at rest;
  !!   A5 refusals by name: a malformed or non-numeric row, an unknown type, a type without the
  !!      synchronisation lift, rows of different columns, A/D not increasing, values out of range,
  !!      a one-row table, more than 32 rows, and the VIVAmp channel on a line without a table.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: D = 0.1_wp, AD_ZERO = 0.3_wp
  INTEGER :: nfail = 0

  CALL write_current('vad_current.dat')
  CALL check_kind('1.0', 'finite-EI')
  CALL check_kind('0.0', 'EI = 0')
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: VIV amplitude deck gates'

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
    WRITE (u, '(A)') 'uniform current along y for the VIV amplitude deck gate'
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

  SUBROUTINE write_deck(path, amp_rows, tmax, ei, cl, amp_channel, cd)
    !! amp_rows: the VIV AMPLITUDE section's data rows, ';'-separated ('' for no section); ei the
    !! type's EI, cl its LINE TYPES Cl (default 0.8), cd its Cd (default 1.2); amp_channel adds
    !! VIVAmp1N26 to the outputs (default: with a section).
    CHARACTER(*), INTENT(IN) :: path, amp_rows, ei
    REAL(wp), INTENT(IN) :: tmax
    CHARACTER(*), INTENT(IN), OPTIONAL :: cl, cd
    LOGICAL, INTENT(IN), OPTIONAL :: amp_channel
    CHARACTER(32) :: cl_col, cd_col
    INTEGER :: u, i0, i1
    LOGICAL :: with_amp
    cl_col = '0.8'
    IF (PRESENT(cl)) cl_col = cl
    cd_col = '1.2'
    IF (PRESENT(cd)) cd_col = cd
    with_amp = LEN_TRIM(amp_rows) > 0
    IF (PRESENT(amp_channel)) with_amp = amp_channel
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (u, '(A)') 'VIV amplitude deck gate: taut cable in a uniform current'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (u, '(A)') 'TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF'
    WRITE (u, '(A)') 'riser      0.10   10.0     5.0e7     -0.8      '//TRIM(ei)//'  '//TRIM(cd_col)// &
      '   1.0   0.0   0.0   '//TRIM(cl_col)//'  0.08  0.18'
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
    IF (LEN_TRIM(amp_rows) > 0) THEN
      WRITE (u, '(A)') '----------------------- VIV AMPLITUDE ---------------------------------------'
      WRITE (u, '(A)') 'TypeName   AD     Cv     CdAmp'
      WRITE (u, '(A)') '(-)        (-)    (-)    (-)'
      i0 = 1
      DO
        i1 = INDEX(amp_rows(i0:), ';')
        IF (i1 == 0) THEN
          WRITE (u, '(A)') TRIM(amp_rows(i0:))
          EXIT
        END IF
        WRITE (u, '(A)') amp_rows(i0:i0 + i1 - 2)
        i0 = i0 + i1
      END DO
    END IF
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '0.01     dtM'
    WRITE (u, '(A)') '100.0    WtrDpth'
    WRITE (u, '(A)') '1025.0   rhoW'
    WRITE (u, '(A)') '9.80665  g'
    WRITE (u, '(A)') 'vad_current.dat  WaterKin'
    WRITE (u, '(A,F8.3,A)') ' ', tmax, '    TMax'
    WRITE (u, '(A)') '---------------------- OUTPUTS -----------------------------------------------'
    WRITE (u, '(A)') 'L1N26PZ'
    WRITE (u, '(A)') 'L1N26PY'
    IF (with_amp) WRITE (u, '(A)') 'VIVAmp1N26'
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

  SUBROUTINE read_out(path, ncol, tab, n)
    !! The rows (time and ncol - 1 channels) of a deck's main output file.
    CHARACTER(*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: ncol
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: tab(:, :)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(1024) :: rec
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

  REAL(wp) FUNCTION window_mean(tab, n, row, t0) RESULT(m)
    REAL(wp), INTENT(IN) :: tab(:, :), t0
    INTEGER, INTENT(IN) :: n, row
    LOGICAL :: keep(n)
    keep = tab(1, 1:n) >= t0
    m = SUM(PACK(tab(row, 1:n), keep))/REAL(MAX(1, COUNT(keep)), wp)
  END FUNCTION window_mean

  REAL(wp) FUNCTION window_ad(tab, n, t0) RESULT(ad)
    !! Half peak-to-peak mid-span cross-flow displacement over D from t0 on.
    REAL(wp), INTENT(IN) :: tab(:, :), t0
    INTEGER, INTENT(IN) :: n
    LOGICAL :: keep(n)
    REAL(wp), ALLOCATABLE :: z(:)
    keep = tab(1, 1:n) >= t0
    z = PACK(tab(2, 1:n), keep)
    ad = 0.5_wp*(MAXVAL(z) - MINVAL(z))/D
  END FUNCTION window_ad

  SUBROUTINE check_kind(ei, label)
    CHARACTER(*), INTENT(IN) :: ei, label
    CHARACTER(16) :: tag
    REAL(wp), ALLOCATABLE :: tn(:, :), ti(:, :), tl(:, :), tg(:, :), tc(:, :)
    INTEGER :: nn, ni, nl, ng, nc
    REAL(wp) :: ad_n, ad_i, amp_i, ad_l, amp_l, y_g, y_static, amp_g
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    tag = MERGE('fe', 'e0', ei /= '0.0')
    CALL write_deck('vad_'//TRIM(tag)//'_none.dat', '', 100.0_wp, ei)
    CALL run('vad_'//TRIM(tag)//'_none.dat', 'vad_'//TRIM(tag)//'_none', ok, msg)
    CALL require(ok, 'A1 ('//label//'): the deck without the section runs: '//TRIM(msg))
    CALL write_deck('vad_'//TRIM(tag)//'_ident.dat', 'riser 0.0 0.8 1.0;riser 1.0 0.8 1.0', 100.0_wp, ei)
    CALL run('vad_'//TRIM(tag)//'_ident.dat', 'vad_'//TRIM(tag)//'_ident', ok, msg)
    CALL require(ok, 'A1 ('//label//'): the identity-table deck runs: '//TRIM(msg))
    CALL write_deck('vad_'//TRIM(tag)//'_limit.dat', 'riser 0.0 0.8;riser 0.1 0.8;riser 0.3 0.0', 100.0_wp, ei)
    CALL run('vad_'//TRIM(tag)//'_limit.dat', 'vad_'//TRIM(tag)//'_limit', ok, msg)
    CALL require(ok, 'A3 ('//label//'): the falling-Cv deck runs: '//TRIM(msg))
    CALL write_deck('vad_'//TRIM(tag)//'_drag.dat', 'riser 0.0 0.0 2.0;riser 1.0 0.0 2.0', 100.0_wp, ei)
    CALL run('vad_'//TRIM(tag)//'_drag.dat', 'vad_'//TRIM(tag)//'_drag', ok, msg)
    CALL require(ok, 'A4 ('//label//'): the Cv 0, CdAmp 2 deck runs: '//TRIM(msg))
    CALL write_deck('vad_'//TRIM(tag)//'_cd24.dat', '', 0.1_wp, ei, cd='2.4')
    CALL run('vad_'//TRIM(tag)//'_cd24.dat', 'vad_'//TRIM(tag)//'_cd24', ok, msg)
    CALL require(ok, 'A4 ('//label//'): the Cd 2.4 deck runs: '//TRIM(msg))
    CALL read_out('vad_'//TRIM(tag)//'_cd24.out', 3, tc, nc)
    CALL read_out('vad_'//TRIM(tag)//'_none.out', 3, tn, nn)
    CALL read_out('vad_'//TRIM(tag)//'_ident.out', 4, ti, ni)
    CALL read_out('vad_'//TRIM(tag)//'_limit.out', 4, tl, nl)
    CALL read_out('vad_'//TRIM(tag)//'_drag.out', 4, tg, ng)
    CALL require(nn > 5000 .AND. ni == nn .AND. nl == nn .AND. ng == nn, 'A1 ('//label//'): the runs wrote '// &
                 'matching rows')
    IF (.NOT. (nn > 5000 .AND. ni == nn .AND. nl == nn .AND. ng == nn)) RETURN
    ! A1: the identity table
    CALL require(ALL(ti(1:3, 1:ni) == tn(1:3, 1:nn)), 'A1 ('//label//'): Cv = Cl and CdAmp 1 is the deck '// &
                 'without the section')
    ! A2: the amplitude estimate against the response
    ad_n = window_ad(tn, nn, 60.0_wp)
    ad_i = window_ad(ti, ni, 60.0_wp)
    amp_i = window_mean(ti, ni, 4, 60.0_wp)/D
    WRITE (*, '(A,3F8.4)') 'A2 ('//label//'): A/D constant Cl, identity table, its VIVAmp/D: ', ad_n, ad_i, amp_i
    CALL require(ad_i > 0.05_wp .AND. ABS(amp_i/ad_i - 1.0_wp) <= 0.01_wp, 'A2 ('//label//'): VIVAmp follows '// &
                 'the response amplitude within 1 %')
    ! A3: self-limiting
    ad_l = window_ad(tl, nl, 60.0_wp)
    amp_l = window_mean(tl, nl, 4, 60.0_wp)/D
    WRITE (*, '(A,2F8.4)') 'A3 ('//label//'): falling-Cv A/D and VIVAmp/D: ', ad_l, amp_l
    CALL require(ad_n > AD_ZERO, 'A3 ('//label//'): the constant-Cl response exceeds the zero of the table')
    CALL require(ABS(amp_l/ad_l - 1.0_wp) <= 0.01_wp, 'A3 ('//label//'): VIVAmp follows it within 1 %')
    CALL require(ad_l > 0.1_wp .AND. ad_l < AD_ZERO, 'A3 ('//label//'): the falling Cv holds the response '// &
                 'inside its falling range')
    ! A4: the drag factor, on the mid-span in-line offset from the line between the ends (y = 0)
    CALL require(nc >= 1, 'A4 ('//label//'): the Cd 2.4 deck wrote its static row')
    IF (nc < 1) RETURN
    y_static = tc(3, 1)
    y_g = window_mean(tg, ng, 3, 80.0_wp)
    amp_g = window_mean(tg, ng, 4, 80.0_wp)/D
    WRITE (*, '(A,2F10.6,ES11.3)') 'A4 ('//label//'): Cd 2.4 static, Cv 0 CdAmp 2 final offsets [m], '// &
      'VIVAmp/D: ', y_static, y_g, amp_g
    CALL require(y_static > 0.0_wp .AND. ABS(y_g/y_static - 1.0_wp) <= 0.005_wp, 'A4 ('//label//'): CdAmp 2 '// &
                 'settles at the offset of Cd 2.4')
    CALL require(window_ad(tg, ng, 80.0_wp) < 1.0e-3_wp .AND. amp_g < 1.0e-2_wp, 'A4 ('//label//'): at rest, '// &
                 'without lift')
  END SUBROUTINE check_kind

  SUBROUTINE check_refusals()
    CHARACTER(2048) :: many
    INTEGER :: k
    CALL refused('riser 0.0 0.8 1.0 2.0', 'TypeName AD Cv [CdAmp]', 'a row of five columns')
    CALL refused('riser 0.0;riser 1.0', 'TypeName AD Cv [CdAmp]', 'a row of two columns')
    CALL refused('riser 0.0 x;riser 1.0 0.8', 'plain numbers', 'a non-numeric Cv')
    CALL refused('cable 0.0 0.8;cable 1.0 0.8', 'undefined LINE TYPE', 'an unknown type')
    CALL refused('riser 0.0 0.8;riser 1.0 0.8', 'give the type LINE TYPES Cl > 0', 'a type without lift', &
                 cl='0.0')
    CALL refused('riser 0.0 0.8 1.0;riser 1.0 0.8', 'same columns', 'rows of different columns')
    CALL refused('riser 0.5 0.8;riser 0.5 0.6', 'increase strictly in AD', 'a repeated AD')
    CALL refused('riser -0.1 0.8;riser 1.0 0.8', 'AD >= 0, Cv >= 0 and CdAmp > 0', 'a negative AD')
    CALL refused('riser 0.0 -0.8;riser 1.0 0.8', 'AD >= 0, Cv >= 0 and CdAmp > 0', 'a negative Cv')
    CALL refused('riser 0.0 0.8 0.0;riser 1.0 0.8 1.0', 'AD >= 0, Cv >= 0 and CdAmp > 0', 'a zero CdAmp')
    CALL refused('riser 0.0 0.8', 'at least two', 'a one-row table')
    many = ''
    DO k = 0, 32
      WRITE (many(LEN_TRIM(many) + 1:), '(A,F5.2,A)') 'riser ', 0.01_wp*REAL(k, wp), ' 0.8;'
    END DO
    many(LEN_TRIM(many):LEN_TRIM(many)) = ' '
    CALL refused(TRIM(many), 'more than 32', 'a table of 33 rows')
    CALL refused('', 'has no VIV AMPLITUDE table', 'VIVAmp on a line without a table', amp_channel=.TRUE.)
  END SUBROUTINE check_refusals

  SUBROUTINE refused(rows, needle, what, cl, amp_channel)
    CHARACTER(*), INTENT(IN) :: rows, needle, what
    CHARACTER(*), INTENT(IN), OPTIONAL :: cl
    LOGICAL, INTENT(IN), OPTIONAL :: amp_channel
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    IF (PRESENT(cl)) THEN
      CALL write_deck('vad_bad.dat', rows, 1.0_wp, '1.0', cl=cl)
    ELSE IF (PRESENT(amp_channel)) THEN
      CALL write_deck('vad_bad.dat', rows, 1.0_wp, '1.0', amp_channel=amp_channel)
    ELSE
      CALL write_deck('vad_bad.dat', rows, 1.0_wp, '1.0')
    END IF
    CALL run('vad_bad.dat', 'vad_bad', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, needle) > 0, 'A5: '//what//' is refused: '//TRIM(msg))
  END SUBROUTINE refused

END PROGRAM test_viv_amplitude_deck
