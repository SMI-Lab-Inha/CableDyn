! File: tests/test_viv_cf_damping_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_cf_damping_deck
  !! The split drag of a VIV DAMPING row (TypeName a b [c]: Morison in line, 1/2 rho D [C_dy |w_y| +
  !! c Cdn |w|] w_y across the ambient flow, C_dy = a + b y0/D of Thorsen et al. 2016) on the taut
  !! finite-EI cable of the VIV amplitude gate (100 m, D 0.1 m, 50 kN, Cl 0.8, 0.15 m/s along y):
  !!   D1 `0 0 1` with a table of Cv = Cl and drag factor 1 is the deck without the section (the
  !!      split recombines Morison's drag to round-off, test_viv_cf_damping_kernel K3), to 1e-4 of
  !!      the response: the VIV phase carries the round-off over 40 s, 1.5e-8 in a Release build and
  !!      1.5e-5 in a Debug one;
  !!   D2 `0 0 1` with Cv 0 and drag factor 2 (no lift, the in-line drag doubled, the cross-flow
  !!      drag Morison's) settles at the static in-line offset of the deck with Cd 2.4, within 0.5 %;
  !!   D3 Thorsen's form `0.31 0.89` replaces the cross-flow drag by the quadratic damping, weaker at
  !!      this amplitude: the steady cross-flow response is larger than with Morison's drag;
  !!   D4 the energy ledger closes with the split: over the steady window the mean of VIVPow1 +
  !!      DragPow1 + DampPow1 is within 1 % of the mean VIV power;
  !!   D5 refusals by name: a row without its table, an EI = 0 line, a malformed row, a negative
  !!      coefficient, a type named twice, an undefined type.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: D = 0.1_wp
  CHARACTER(*), PARAMETER :: NEUTRAL = 'riser 0.0 0.8 1.0;riser 1.0 0.8 1.0'
  INTEGER :: nfail = 0

  CALL write_current('vcd_current.dat')
  CALL check_identity()
  CALL check_inline_amplification()
  CALL check_thorsen()
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: VIV cross-flow damping deck gates'

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
    WRITE (u, '(A)') 'uniform current along y for the VIV cross-flow damping deck gate'
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

  SUBROUTINE write_rows(u, rows)
    INTEGER, INTENT(IN) :: u
    CHARACTER(*), INTENT(IN) :: rows
    INTEGER :: i0, i1
    i0 = 1
    DO
      i1 = INDEX(rows(i0:), ';')
      IF (i1 == 0) THEN
        WRITE (u, '(A)') TRIM(rows(i0:))
        EXIT
      END IF
      WRITE (u, '(A)') rows(i0:i0 + i1 - 2)
      i0 = i0 + i1
    END DO
  END SUBROUTINE write_rows

  SUBROUTINE write_deck(path, amp_rows, damp_rows, tmax, ei, cd, power)
    !! amp_rows/damp_rows: the VIV AMPLITUDE and VIV DAMPING data rows, ';'-separated ('' for no
    !! section); ei the type's EI (default 1.0), cd its Cd (default 1.2); power adds the load-power
    !! channels.
    CHARACTER(*), INTENT(IN) :: path, amp_rows, damp_rows
    REAL(wp), INTENT(IN) :: tmax
    CHARACTER(*), INTENT(IN), OPTIONAL :: ei, cd
    LOGICAL, INTENT(IN), OPTIONAL :: power
    CHARACTER(32) :: ei_col, cd_col
    INTEGER :: u
    ei_col = '1.0'
    IF (PRESENT(ei)) ei_col = ei
    cd_col = '1.2'
    IF (PRESENT(cd)) cd_col = cd
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (u, '(A)') 'VIV cross-flow damping deck gate: taut cable in a uniform current'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (u, '(A)') 'TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF'
    WRITE (u, '(A)') 'riser      0.10   10.0     5.0e7     -0.8      '//TRIM(ei_col)//'  '//TRIM(cd_col)// &
      '   1.0   0.0   0.0   0.8  0.08  0.18'
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
      CALL write_rows(u, amp_rows)
    END IF
    IF (LEN_TRIM(damp_rows) > 0) THEN
      WRITE (u, '(A)') '----------------------- VIV DAMPING -----------------------------------------'
      WRITE (u, '(A)') 'TypeName   a      b      c'
      WRITE (u, '(A)') '(-)        (-)    (-)    (-)'
      CALL write_rows(u, damp_rows)
    END IF
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '0.01     dtM'
    WRITE (u, '(A)') '100.0    WtrDpth'
    WRITE (u, '(A)') '1025.0   rhoW'
    WRITE (u, '(A)') '9.80665  g'
    WRITE (u, '(A)') 'vcd_current.dat  WaterKin'
    WRITE (u, '(A,F8.3,A)') ' ', tmax, '    TMax'
    WRITE (u, '(A)') '---------------------- OUTPUTS -----------------------------------------------'
    WRITE (u, '(A)') 'L1N26PZ'
    WRITE (u, '(A)') 'L1N26PY'
    IF (PRESENT(power)) THEN
      IF (power) THEN
        WRITE (u, '(A)') 'VIVPow1'
        WRITE (u, '(A)') 'DragPow1'
        WRITE (u, '(A)') 'DampPow1'
      END IF
    END IF
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

  SUBROUTINE read_out(path, tab, n, ncol)
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: tab(:, :)
    INTEGER, INTENT(OUT) :: n
    INTEGER, INTENT(IN) :: ncol
    CHARACTER(1024) :: rec
    INTEGER :: u, ios
    ALLOCATE (tab(ncol, 40000))
    n = 0
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    READ (u, '(A)') rec
    READ (u, '(A)') rec
    DO
      READ (u, '(A)', IOSTAT=ios) rec
      IF (ios /= 0) EXIT
      IF (LEN_TRIM(rec) == 0) CYCLE
      n = n + 1
      READ (rec, *) tab(:, n)
    END DO
    CLOSE (u)
  END SUBROUTINE read_out

  SUBROUTINE check_identity()
    REAL(wp), ALLOCATABLE :: ta(:, :), tb(:, :)
    INTEGER :: na, nb
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    REAL(wp) :: dz
    CALL write_deck('vcd_morison.dat', NEUTRAL, '', 40.0_wp)
    CALL run('vcd_morison.dat', 'vcd_morison', ok, msg)
    CALL require(ok, 'D1: the Morison deck runs: '//TRIM(msg))
    CALL write_deck('vcd_split.dat', NEUTRAL, 'riser 0 0 1', 40.0_wp)
    CALL run('vcd_split.dat', 'vcd_split', ok, msg)
    CALL require(ok, 'D1: the split deck runs: '//TRIM(msg))
    CALL read_out('vcd_morison.out', ta, na, 3)
    CALL read_out('vcd_split.out', tb, nb, 3)
    IF (na /= nb .OR. na < 100) THEN
      CALL require(.FALSE., 'D1: both runs wrote the same rows'); RETURN
    END IF
    dz = MAXVAL(ABS(ta(2:3, :na) - tb(2:3, :nb)))/MAXVAL(ABS(ta(2, :na) - ta(2, 1)))
    WRITE (*, '(A,ES10.2)') 'D1 0 0 1 against Morison: max response difference / CF range ', dz
    ! the split and recombined drag differ by round-off, which the VIV phase dynamics carry over 40 s
    CALL require(dz <= 1.0e-4_wp, 'D1: the split with c = 1 recombines Morison''s drag')
  END SUBROUTINE check_identity

  SUBROUTINE check_inline_amplification()
    REAL(wp), ALLOCATABLE :: ta(:, :), tb(:, :)
    INTEGER :: na, nb
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    REAL(wp) :: ya, yb
    CALL write_deck('vcd_amp2.dat', 'riser 0.0 0.0 2.0;riser 1.0 0.0 2.0', 'riser 0 0 1', 30.0_wp)
    CALL run('vcd_amp2.dat', 'vcd_amp2', ok, msg)
    CALL require(ok, 'D2: the in-line amplification deck runs: '//TRIM(msg))
    CALL write_deck('vcd_cd24.dat', 'riser 0.0 0.0 1.0;riser 1.0 0.0 1.0', '', 30.0_wp, cd='2.4')
    CALL run('vcd_cd24.dat', 'vcd_cd24', ok, msg)
    CALL require(ok, 'D2: the Cd 2.4 deck runs: '//TRIM(msg))
    CALL read_out('vcd_amp2.out', ta, na, 3)
    CALL read_out('vcd_cd24.out', tb, nb, 3)
    IF (na < 10 .OR. nb < 10) RETURN
    ya = ta(3, na) - 0.0_wp
    yb = tb(3, nb) - 0.0_wp
    WRITE (*, '(A,2ES14.6)') 'D2 in-line offset, split with CdAmp 2 and Cd 2.4: ', ya, yb
    CALL require(ABS(ya - yb) <= 5.0e-3_wp*ABS(yb), 'D2: the in-line drag alone is amplified')
  END SUBROUTINE check_inline_amplification

  SUBROUTINE check_thorsen()
    REAL(wp), ALLOCATABLE :: ta(:, :), tb(:, :)
    INTEGER :: na, nb, i0
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    REAL(wp) :: am, aq, pv, ps
    CALL write_deck('vcd_q.dat', NEUTRAL, 'riser 0.31 0.89', 90.0_wp, power=.TRUE.)
    CALL run('vcd_q.dat', 'vcd_q', ok, msg)
    CALL require(ok, 'D3: the Thorsen deck runs: '//TRIM(msg))
    CALL read_out('vcd_q.out', ta, na, 6)
    CALL read_out('vcd_morison.out', tb, nb, 3)
    IF (na < 100 .OR. nb < 100) RETURN
    i0 = na - 2000
    aq = 0.5_wp*(MAXVAL(ta(2, i0:na)) - MINVAL(ta(2, i0:na)))/D
    am = 0.5_wp*(MAXVAL(tb(2, nb - 1000:nb)) - MINVAL(tb(2, nb - 1000:nb)))/D
    WRITE (*, '(A,F8.4,A,F8.4)') 'D3 steady A/D with Thorsen''s damping ', aq, ', with Morison''s ', am
    CALL require(aq > 1.02_wp*am, 'D3: the quadratic cross-flow damping is weaker at this amplitude')
    pv = SUM(ta(4, i0:na))/REAL(na - i0 + 1, wp)
    ps = SUM(ta(4, i0:na) + ta(5, i0:na) + ta(6, i0:na))/REAL(na - i0 + 1, wp)
    WRITE (*, '(A,2ES12.4)') 'D4 mean VIV power and mean VIV + drag + damping power: ', pv, ps
    CALL require(pv > 0.0_wp .AND. ABS(ps) <= 1.0e-2_wp*pv, 'D4: the load powers balance with the split drag')
  END SUBROUTINE check_thorsen

  SUBROUTINE refuse(label, amp_rows, damp_rows, expect, ei)
    CHARACTER(*), INTENT(IN) :: label, amp_rows, damp_rows, expect
    CHARACTER(*), INTENT(IN), OPTIONAL :: ei
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    IF (PRESENT(ei)) THEN
      CALL write_deck('vcd_refuse.dat', amp_rows, damp_rows, 0.0_wp, ei=ei)
    ELSE
      CALL write_deck('vcd_refuse.dat', amp_rows, damp_rows, 0.0_wp)
    END IF
    CALL run('vcd_refuse.dat', 'vcd_refuse', ok, msg)
    WRITE (*, '(A)') 'D5 '//label//': '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, expect) > 0, 'D5: '//label)
  END SUBROUTINE refuse

  SUBROUTINE check_refusals()
    CALL refuse('no amplitude table', '', 'riser 0.31 0.89', 'give the type Cl > 0 and a VIV AMPLITUDE table')
    CALL refuse('an EI = 0 line', NEUTRAL, 'riser 0.31 0.89', 'acts on finite-EI cables (EI > 0) only', ei='0.0')
    CALL refuse('a malformed row', NEUTRAL, 'riser 0.31', 'rows must be: TypeName a b [c]')
    CALL refuse('a negative coefficient', NEUTRAL, 'riser -0.1 0.89', 'needs a, b and c >= 0')
    CALL refuse('a type named twice', NEUTRAL, 'riser 0.31 0.89;riser 0 0 1', 'twice')
    CALL refuse('an undefined type', NEUTRAL, 'cable 0.31 0.89', 'undefined LINE TYPE')
  END SUBROUTINE check_refusals

END PROGRAM test_viv_cf_damping_deck
