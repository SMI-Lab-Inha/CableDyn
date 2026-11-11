! File: tests/test_load_power_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_load_power_deck
  !! The load power channels VIVPow, DragPow and DampPow (<L> the line total in W, <L>N<J> per unit
  !! length about node J in W/m) on the in-line VIV deck of tests/test_viv_inline_deck.f90: a taut
  !! line, 100 m pinned at 50 m depth in a uniform 0.15 m/s current along y, Cl 0.8, ClX 0.85.
  !! Locked in, the energy the line holds only oscillates, so over 50-100 s
  !!   mean(VIVPow + DragPow + DampPow) = 0
  !! to the window's change of energy.
  !!   P1 the finite-EI cable: the VIV loads put power in, the drag takes it out, the damping power
  !!      is never positive, and the three balance within 1 % of the VIV power;
  !!   P2 the per-length channels of all 51 nodes, times their half-element lengths, sum to the
  !!      line total;
  !!   P3 the same line as an EI = 0 line (the single-line march) balances within 1 % and puts in
  !!      the cable's VIV power within 2 %;
  !!   P4 the EI = 0 line split at a free point into two lines (the multi-line system) balances
  !!      within 1 % over both lines;
  !!   P5 malformed channels, an unknown line, node 0 and a node past the end are refused.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  INTEGER, PARAMETER :: NSEG = 50
  REAL(wp), PARAMETER :: LEN0 = 99.9001_wp
  INTEGER :: nfail = 0

  CALL write_current('lpd_current.dat')
  CALL check_cable()
  CALL check_ei0()
  CALL check_system()
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: load power channel gates'

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
    WRITE (u, '(A)') 'uniform current along y for the load power deck gate'
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

  SUBROUTINE write_deck(path, ei, channels, split, tmax)
    !! channels: the OUTPUTS rows, ';'-separated; split: the EI = 0 line as two lines joined at
    !! a free point at mid-span.
    CHARACTER(*), INTENT(IN) :: path, ei, channels
    LOGICAL, INTENT(IN) :: split
    REAL(wp), INTENT(IN) :: tmax
    INTEGER :: u, i0, i1
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (u, '(A)') 'load power deck gate: taut line with in-line VIV in a uniform current'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (u, '(A)') 'TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF'
    WRITE (u, '(A)') '(name)     (m)    (kg/m)   (N)       (N-s/-)   (-)  (-)   (-)   (-)   (-)   (-)  (-)   (-)'
    WRITE (u, '(A)') 'riser      0.10   10.0     5.0e7     -0.8      '//ei//'  1.2   1.0   0.0   0.0   0.8  0.08  0.18'
    WRITE (u, '(A)') '----------------------- POINTS ----------------------------------------------'
    WRITE (u, '(A)') 'ID   Attachment   X       Y     Z      M    V     CdA   CA'
    WRITE (u, '(A)') '(#)  (-)          (m)     (m)   (m)    (kg) (m^3) (m^2) (-)'
    WRITE (u, '(A)') '1    Coupled      0.0     0.0   -50.0  0    0     0     0'
    WRITE (u, '(A)') '2    Fixed        100.0   0.0   -50.0  0    0     0     0'
    IF (split) WRITE (u, '(A)') '3    Free         50.0    0.0   -50.0  0    0     0     0'
    WRITE (u, '(A)') '----------------------- LINES -----------------------------------------------'
    WRITE (u, '(A)') 'ID   LineType   AttachA  AttachB  UnstrLen  NumSegs  Outputs'
    WRITE (u, '(A)') '(#)  (name)     (#)      (#)      (m)       (-)      (-)'
    IF (split) THEN
      WRITE (u, '(A)') '1    riser      1        3        49.95005   25       -'
      WRITE (u, '(A)') '2    riser      3        2        49.95005   25       -'
    ELSE
      WRITE (u, '(A)') '1    riser      1        2        99.900100  50       -'
    END IF
    WRITE (u, '(A)') '----------------------- VIV -------------------------------------------------'
    WRITE (u, '(A)') 'TypeName   ClX    alpha'
    WRITE (u, '(A)') '(name)     (-)    (-)'
    WRITE (u, '(A)') 'riser      0.85   0.15'
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '0.01     dtM'
    WRITE (u, '(A)') '100.0    WtrDpth'
    WRITE (u, '(A)') '1025.0   rhoW'
    WRITE (u, '(A)') '9.80665  g'
    WRITE (u, '(A)') 'lpd_current.dat  WaterKin'
    WRITE (u, '(A)') '0.25     VIVrmsWindow'
    WRITE (u, '(A,F8.3,A)') ' ', tmax, '    TMax'
    WRITE (u, '(A)') '---------------------- OUTPUTS -----------------------------------------------'
    i0 = 1
    DO
      i1 = INDEX(channels(i0:), ';')
      IF (i1 == 0) THEN
        WRITE (u, '(A)') TRIM(channels(i0:))
        EXIT
      END IF
      WRITE (u, '(A)') channels(i0:i0 + i1 - 2)
      i0 = i0 + i1
    END DO
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
    !! The rows of a deck's main output file (time and ncol channels) from 50 s on.
    CHARACTER(*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: ncol
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: tab(:, :)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(8192) :: rec
    REAL(wp) :: row(ncol + 1)
    INTEGER :: u, ios
    ALLOCATE (tab(ncol + 1, 20000))
    n = 0
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    READ (u, '(A)') rec
    READ (u, '(A)') rec
    DO
      READ (u, '(A)', IOSTAT=ios) rec
      IF (ios /= 0) EXIT
      IF (LEN_TRIM(rec) == 0) CYCLE
      READ (rec, *) row
      IF (row(1) < 50.0_wp - 1.0e-9_wp) CYCLE
      IF (n == SIZE(tab, 2)) EXIT
      n = n + 1
      tab(:, n) = row
    END DO
    CLOSE (u)
  END SUBROUTINE read_out

  SUBROUTINE balance(label, tab, n, lines, p_viv)
    !! The mean of the three powers over the window (columns 2-4 of a line, lines in a row) and
    !! the balance check; p_viv returns the mean VIV power.
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: tab(:, :)
    INTEGER, INTENT(IN) :: n, lines
    REAL(wp), INTENT(OUT) :: p_viv
    REAL(wp) :: m(3)
    INTEGER :: k, j
    m = 0.0_wp
    DO j = 1, lines
      DO k = 1, 3
        m(k) = m(k) + SUM(tab(1 + 3*(j - 1) + k, 1:n))/REAL(n, wp)
      END DO
    END DO
    p_viv = m(1)
    WRITE (*, '(2A,3F10.4,A,ES10.2)') label, ': mean VIVPow, DragPow, DampPow [W] ', m, &
      ', balance/VIV ', SUM(m)/m(1)
    CALL require(m(1) > 0.0_wp .AND. m(2) < 0.0_wp, label//': the VIV loads put power in, the drag takes it out')
    CALL require(ABS(SUM(m)) <= 0.01_wp*m(1), label//': the powers balance within 1 % of the VIV power')
    DO j = 1, lines
      CALL require(ALL(tab(1 + 3*j, 1:n) <= 0.0_wp), label//': the damping power is never positive')
    END DO
  END SUBROUTINE balance

  SUBROUTINE check_cable()
    REAL(wp), ALLOCATABLE :: tab(:, :)
    REAL(wp) :: pv, trib, s, sabs
    INTEGER :: n, j, k
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CHARACTER(4096) :: chans
    CALL write_deck('lpd_cab.dat', '1.0', 'VIVPow1;DragPow1;DampPow1;DampPow1N26', .FALSE., 100.0_wp)
    CALL run('lpd_cab.dat', 'lpd_cab', ok, msg)
    CALL require(ok, 'P1: the finite-EI deck runs: '//TRIM(msg))
    IF (.NOT. ok) RETURN
    CALL read_out('lpd_cab.out', 4, tab, n)
    CALL require(n > 1000, 'P1: the run wrote its rows')
    IF (n <= 1000) RETURN
    CALL balance('P1 finite-EI cable', tab, n, 1, pv)
    CALL require(ALL(tab(5, 1:n) <= 0.0_wp), 'P1: the damping power per length is never positive')
    ! P2: every node's power per length, times its half elements, sums to the line total
    chans = 'VIVPow1'
    DO j = 1, NSEG + 1
      chans = TRIM(chans)//';VIVPow1N'//TRIM(int_text(j))
    END DO
    CALL write_deck('lpd_nodes.dat', '1.0', TRIM(chans), .FALSE., 52.0_wp)
    CALL run('lpd_nodes.dat', 'lpd_nodes', ok, msg)
    CALL require(ok, 'P2: the all-node deck runs: '//TRIM(msg))
    IF (.NOT. ok) RETURN
    CALL read_out('lpd_nodes.out', NSEG + 2, tab, n)
    CALL require(n > 100, 'P2: the run wrote its rows')
    IF (n <= 100) RETURN
    DO k = 1, n
      s = 0.0_wp
      sabs = 0.0_wp
      DO j = 1, NSEG + 1
        trib = LEN0/REAL(NSEG, wp)
        IF (j == 1 .OR. j == NSEG + 1) trib = 0.5_wp*trib
        s = s + tab(2 + j, k)*trib
        sabs = sabs + ABS(tab(2 + j, k))*trib
      END DO
      IF (ABS(s - tab(2, k)) > 1.0e-5_wp*MAX(sabs, 1.0e-12_wp)) THEN
        CALL require(.FALSE., 'P2: the per-length powers sum to the line total')
        EXIT
      END IF
    END DO
  END SUBROUTINE check_cable

  SUBROUTINE check_ei0()
    REAL(wp), ALLOCATABLE :: tab(:, :)
    REAL(wp) :: pv0, pvf
    INTEGER :: n
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('lpd_ei0.dat', '0.0', 'VIVPow1;DragPow1;DampPow1', .FALSE., 100.0_wp)
    CALL run('lpd_ei0.dat', 'lpd_ei0', ok, msg)
    CALL require(ok, 'P3: the EI = 0 deck runs: '//TRIM(msg))
    IF (.NOT. ok) RETURN
    CALL read_out('lpd_ei0.out', 3, tab, n)
    CALL require(n > 1000, 'P3: the run wrote its rows')
    IF (n <= 1000) RETURN
    CALL balance('P3 EI = 0 line', tab, n, 1, pv0)
    CALL read_out('lpd_cab.out', 4, tab, n)
    IF (n <= 1000) RETURN
    pvf = SUM(tab(2, 1:n))/REAL(n, wp)
    CALL require(ABS(pv0/pvf - 1.0_wp) <= 0.02_wp, 'P3: the EI = 0 line takes the cable''s VIV power within 2 %')
  END SUBROUTINE check_ei0

  SUBROUTINE check_system()
    REAL(wp), ALLOCATABLE :: tab(:, :)
    REAL(wp) :: pv
    INTEGER :: n
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('lpd_sys.dat', '0.0', 'VIVPow1;DragPow1;DampPow1;VIVPow2;DragPow2;DampPow2', .TRUE., 100.0_wp)
    CALL run('lpd_sys.dat', 'lpd_sys', ok, msg)
    CALL require(ok, 'P4: the two-line deck runs: '//TRIM(msg))
    IF (.NOT. ok) RETURN
    CALL read_out('lpd_sys.out', 6, tab, n)
    CALL require(n > 1000, 'P4: the run wrote its rows')
    IF (n <= 1000) RETURN
    CALL balance('P4 two EI = 0 lines and a free point', tab, n, 2, pv)
  END SUBROUTINE check_system

  SUBROUTINE check_refusals()
    CALL refused('VIVPow0', 'bad load power channel', 'line 0')
    CALL refused('DragPow7', 'bad load power channel', 'an unknown line')
    CALL refused('DampPowx', 'bad load power channel', 'a channel without a line id')
    CALL refused('VIVPow1N0', 'node numbers start at 1', 'node 0')
    CALL refused('DragPow1N52', 'exceeds line node count', 'a node past End B')
  END SUBROUTINE check_refusals

  SUBROUTINE refused(channel, needle, what)
    CHARACTER(*), INTENT(IN) :: channel, needle, what
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('lpd_bad.dat', '1.0', channel, .FALSE., 1.0_wp)
    CALL run('lpd_bad.dat', 'lpd_bad', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, needle) > 0, 'P5: '//what//' is refused: '//TRIM(msg))
  END SUBROUTINE refused

  FUNCTION int_text(i) RESULT(s)
    INTEGER, INTENT(IN) :: i
    CHARACTER(:), ALLOCATABLE :: s
    CHARACTER(16) :: b
    WRITE (b, '(I0)') i
    s = TRIM(b)
  END FUNCTION int_text

END PROGRAM test_load_power_deck
