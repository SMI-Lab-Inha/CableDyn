! File: tests/test_viv_frame_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_frame_deck
  !! The VIV FRAME section (TypeName Tf) on the synchronisation model: with Tf > 0 the amplitude
  !! estimate (and so the VIV AMPLITUDE coefficients Cv and CdAmp) and the axes and blend of the
  !! VIV DAMPING split take the ambient flow relative to the cable's velocity low-pass filtered
  !! over Tf. The cable is the taut cable of the VIV amplitude and damping gates (100 m, D 0.1 m,
  !! 50 kN, pinned at 50 m depth, 50 elements, Cl 0.8), its type with a VIV AMPLITUDE table (Cv 0.8,
  !! 0.6, 0.0 and CdAmp 1.0, 1.3, 1.6 at A/D 0, 0.5, 1) and, on the finite-EI cable (EI 1 N m^2),
  !! Thorsen's VIV DAMPING (0.31 0.89).
  !!   G1 a Tf 0 row is the deck without the section, bit for bit, on the finite-EI cable and on the
  !!      EI = 0 line;
  !!   G2 the finite-EI cable towed through still water (both ends at 0.15 m/s against y after a
  !!      10 s half-cosine ramp) with Tf 10 s responds as the cable held in a 0.15 m/s current: the
  !!      A/D over 50-80 s within 5 % (measured 1.5 %) and the frequency within 0.005 Hz (the
  !!      30 s window resolves 0.033 Hz; the scan peak moves by 0.002 Hz). The towed cable without
  !!      the frame takes Cv and CdAmp at A/D 0 and the Morison drag (no ambient flow) and responds
  !!      otherwise (A/D more than 10 % off; measured 16 %);
  !!   G3 held in the current, Tf 10 s leaves the response as the ambient flow gives it (A/D within
  !!      2 %, frequency within 0.005 Hz), on both cores (an EI = 0 line cannot run between two
  !!      moving points, so the towed case is finite-EI only);
  !!   G4 refusals by name: a VIV FRAME row on a synchronisation type without a VIV AMPLITUDE table.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, D = 0.1_wp
  CHARACTER(*), PARAMETER :: AMP = 'riser 0.0 0.8 1.0;riser 0.5 0.6 1.3;riser 1.0 0.0 1.6'
  REAL(wp), PARAMETER :: TMAX = 80.0_wp, T0 = 50.0_wp
  INTEGER :: nfail = 0
  REAL(wp) :: fe_cur(2), fe_cur_tf(2), l_cur(2), l_cur_tf(2)

  CALL write_current('vfd_current.dat')
  CALL write_tow('vfd_tow.txt')
  CALL check_fixed('1.0', 'riser 0.31 0.89', 'vfd_fe', fe_cur, fe_cur_tf)
  CALL check_fixed('0.0', '', 'vfd_l', l_cur, l_cur_tf)
  CALL check_tow(fe_cur)
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: VIV FRAME deck gates'

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
    WRITE (u, '(A)') 'uniform current along y for the VIV FRAME deck gate'
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

  SUBROUTINE write_tow(path)
    !! Both ends (points 1 and 2) moved along -y at 0.15 m/s after a 10 s half-cosine velocity ramp,
    !! every 0.01 s to TMAX: y = -U (t/2 - T sin(pi t/T)/(2 pi)) in the ramp (T 10 s).
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), PARAMETER :: U = 0.15_wp, TR = 10.0_wp
    INTEGER :: u_, k
    REAL(wp) :: t, y, v, a
    OPEN (NEWUNIT=u_, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u_, '(A)') '# time point x y z vx vy vz ax ay az'
    DO k = 0, NINT(100.0_wp*TMAX)
      t = 0.01_wp*REAL(k, wp)
      IF (t < TR) THEN
        y = -U*(0.5_wp*t - TR*SIN(PI*t/TR)/(2.0_wp*PI))
        v = -U*0.5_wp*(1.0_wp - COS(PI*t/TR))
        a = -U*0.5_wp*(PI/TR)*SIN(PI*t/TR)
      ELSE
        y = -U*(0.5_wp*TR + (t - TR))
        v = -U
        a = 0.0_wp
      END IF
      WRITE (u_, '(F9.3,A,3ES17.9,A,ES17.9,A,ES17.9,A)') t, ' 1 ', 0.0_wp, y, -50.0_wp, ' 0 ', v, ' 0 0 ', a, ' 0'
      WRITE (u_, '(F9.3,A,3ES17.9,A,ES17.9,A,ES17.9,A)') t, ' 2 ', 100.0_wp, y, -50.0_wp, ' 0 ', v, ' 0 0 ', a, ' 0'
    END DO
    CLOSE (u_)
  END SUBROUTINE write_tow

  SUBROUTINE write_rows(u, title, header, rows)
    INTEGER, INTENT(IN) :: u
    CHARACTER(*), INTENT(IN) :: title, header, rows
    INTEGER :: i0, i1
    WRITE (u, '(A)') '----------------------- '//title//' ---------------------------------------'
    WRITE (u, '(A)') header
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

  SUBROUTINE write_deck(path, ei, amp_rows, damp_rows, tmax_, tow, frame_rows)
    !! ei the type's EI; amp_rows/damp_rows: the VIV AMPLITUDE and VIV DAMPING data rows,
    !! ';'-separated ('' for no section); tow: no current, both ends Coupled and moved by the motion
    !! file vfd_tow.txt; frame_rows: the VIV FRAME data rows (no section without it).
    CHARACTER(*), INTENT(IN) :: path, ei, amp_rows, damp_rows
    REAL(wp), INTENT(IN) :: tmax_
    LOGICAL, INTENT(IN) :: tow
    CHARACTER(*), INTENT(IN), OPTIONAL :: frame_rows
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (u, '(A)') 'VIV FRAME deck gate: taut cable in a uniform current or towed through still water'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (u, '(A)') 'TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF'
    WRITE (u, '(A)') 'riser      0.10   10.0     5.0e7     -0.8      '//ei//'  1.2   1.0   0.0   0.0   0.8  0.08  0.18'
    WRITE (u, '(A)') '----------------------- POINTS ----------------------------------------------'
    WRITE (u, '(A)') 'ID   Attachment   X       Y     Z      M    V     CdA   CA'
    WRITE (u, '(A)') '1    Coupled      0.0     0.0   -50.0  0    0     0     0'
    IF (tow) THEN
      WRITE (u, '(A)') '2    Coupled      100.0   0.0   -50.0  0    0     0     0'
    ELSE
      WRITE (u, '(A)') '2    Fixed        100.0   0.0   -50.0  0    0     0     0'
    END IF
    WRITE (u, '(A)') '----------------------- LINES -----------------------------------------------'
    WRITE (u, '(A)') 'ID   NodeA  NodeB  Outputs'
    WRITE (u, '(A)') '1    1      2      -'
    WRITE (u, '(A)') '----------------------- SECTIONS --------------------------------------------'
    WRITE (u, '(A)') 'LineID  LineType  Length     NumSegs'
    WRITE (u, '(A)') '1       riser     99.900100  50'
    IF (LEN_TRIM(amp_rows) > 0) CALL write_rows(u, 'VIV AMPLITUDE', 'TypeName   AD     Cv     CdAmp', amp_rows)
    IF (LEN_TRIM(damp_rows) > 0) CALL write_rows(u, 'VIV DAMPING', 'TypeName   a      b', damp_rows)
    IF (PRESENT(frame_rows)) CALL write_rows(u, 'VIV FRAME', 'TypeName   Tf', frame_rows)
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '0.01     dtM'
    WRITE (u, '(A)') '100.0    WtrDpth'
    WRITE (u, '(A)') '1025.0   rhoW'
    WRITE (u, '(A)') '9.80665  g'
    IF (tow) THEN
      WRITE (u, '(A)') 'vfd_tow.txt  motionFile'
    ELSE
      WRITE (u, '(A)') 'vfd_current.dat  WaterKin'
    END IF
    WRITE (u, '(A,F8.3,A)') ' ', tmax_, '    TMax'
    WRITE (u, '(A)') '---------------------- OUTPUTS -----------------------------------------------'
    WRITE (u, '(A)') 'L1N25PZ'
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
    !! The rows (time, L1N25PZ) of a deck's main output file.
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

  SUBROUTINE stats(root, res, ok)
    !! Half the peak-to-peak mid-span cross-flow excursion over D and the dominant frequency over
    !! T0-TMAX of the run root.
    CHARACTER(*), INTENT(IN) :: root
    REAL(wp), INTENT(OUT) :: res(2)
    LOGICAL, INTENT(OUT) :: ok
    REAL(wp), ALLOCATABLE :: tab(:, :), tt(:), z(:)
    INTEGER :: n
    LOGICAL, ALLOCATABLE :: keep(:)
    res = 0.0_wp
    CALL read_out(root//'.out', tab, n)
    ok = n > 1000
    IF (.NOT. ok) RETURN
    keep = tab(1, 1:n) >= T0
    tt = PACK(tab(1, 1:n), keep)
    z = PACK(tab(2, 1:n), keep)
    res = [0.5_wp*(MAXVAL(z) - MINVAL(z))/D, peak_frequency(tt, z)]
  END SUBROUTINE stats

  LOGICAL FUNCTION same_file(a, b) RESULT(same)
    !! The two text files are the same, line for line.
    CHARACTER(*), INTENT(IN) :: a, b
    CHARACTER(1024) :: ra, rb
    INTEGER :: ua, ub, ia, ib
    same = .FALSE.
    OPEN (NEWUNIT=ua, FILE=a, STATUS='OLD', ACTION='READ', IOSTAT=ia)
    IF (ia /= 0) RETURN
    OPEN (NEWUNIT=ub, FILE=b, STATUS='OLD', ACTION='READ', IOSTAT=ib)
    IF (ib /= 0) THEN
      CLOSE (ua); RETURN
    END IF
    DO
      READ (ua, '(A)', IOSTAT=ia) ra
      READ (ub, '(A)', IOSTAT=ib) rb
      IF (ia /= 0 .OR. ib /= 0) EXIT
      IF (ra /= rb) EXIT
    END DO
    same = ia /= 0 .AND. ib /= 0
    CLOSE (ua)
    CLOSE (ub)
  END FUNCTION same_file

  SUBROUTINE check_fixed(ei, damp_rows, tag, res, res_tf)
    !! G1 and G3 on the cable of EI ei held in the current: without the section, with a Tf 0 row
    !! and with Tf 10 s; res and res_tf (A/D, f) without the section and with Tf 10 s.
    CHARACTER(*), INTENT(IN) :: ei, damp_rows, tag
    REAL(wp), INTENT(OUT) :: res(2), res_tf(2)
    LOGICAL :: ok, ok2
    CHARACTER(1024) :: msg
    CHARACTER(:), ALLOCATABLE :: lab
    lab = ' (EI '//ei//')'
    CALL write_deck(tag//'_amb.dat', ei, AMP, damp_rows, TMAX, .FALSE.)
    CALL run(tag//'_amb.dat', tag//'_amb', ok, msg)
    CALL require(ok, 'G1: the deck without the section runs'//lab//': '//TRIM(msg))
    CALL write_deck(tag//'_t0.dat', ei, AMP, damp_rows, TMAX, .FALSE., frame_rows='riser 0')
    CALL run(tag//'_t0.dat', tag//'_t0', ok, msg)
    CALL require(ok, 'G1: the deck with Tf 0 runs'//lab//': '//TRIM(msg))
    CALL require(same_file(tag//'_amb.out', tag//'_t0.out'), 'G1: a Tf 0 row is the deck without it, bit for '// &
                 'bit'//lab)
    CALL write_deck(tag//'_tf.dat', ei, AMP, damp_rows, TMAX, .FALSE., frame_rows='riser 10')
    CALL run(tag//'_tf.dat', tag//'_tf', ok, msg)
    CALL require(ok, 'G3: the deck with Tf 10 s runs'//lab//': '//TRIM(msg))
    CALL stats(tag//'_amb', res, ok)
    CALL stats(tag//'_tf', res_tf, ok2)
    CALL require(ok .AND. ok2, 'G3: the runs wrote their rows'//lab)
    WRITE (*, '(A,2F8.4,A,2F8.4)') 'G3'//lab//': in the current A/D, f [Hz] with Tf 10 s ', res_tf, &
      '; without ', res
    CALL require(res(1) > 0.2_wp, 'G3: the cable locks in'//lab)
    CALL require(ABS(res_tf(1)/res(1) - 1.0_wp) <= 0.02_wp .AND. ABS(res_tf(2) - res(2)) <= 0.005_wp, &
                 'G3: held in the current, Tf leaves the response as the ambient flow gives it'//lab)
  END SUBROUTINE check_fixed

  SUBROUTINE check_tow(ref)
    !! G2 on the finite-EI cable: towed with Tf 10 s and without it, against the cable held in the
    !! current (ref, from check_fixed).
    REAL(wp), INTENT(IN) :: ref(2)
    REAL(wp) :: r_tf(2), r_0(2)
    LOGICAL :: ok, ok2
    CHARACTER(1024) :: msg
    CALL write_deck('vfd_tow_tf.dat', '1.0', AMP, 'riser 0.31 0.89', TMAX, .TRUE., frame_rows='riser 10')
    CALL run('vfd_tow_tf.dat', 'vfd_tow_tf', ok, msg)
    CALL require(ok, 'G2: the towed deck with Tf runs: '//TRIM(msg))
    CALL write_deck('vfd_tow_0.dat', '1.0', AMP, 'riser 0.31 0.89', TMAX, .TRUE.)
    CALL run('vfd_tow_0.dat', 'vfd_tow_0', ok, msg)
    CALL require(ok, 'G2: the towed deck without Tf runs: '//TRIM(msg))
    CALL stats('vfd_tow_tf', r_tf, ok)
    CALL stats('vfd_tow_0', r_0, ok2)
    CALL require(ok .AND. ok2, 'G2: the towed runs wrote their rows')
    WRITE (*, '(A,2F8.4,A,2F8.4,A,2F8.4)') 'G2: towed A/D, f [Hz] with Tf 10 s ', r_tf, '; without ', r_0, &
      '; held in the current ', ref
    CALL require(ABS(r_tf(1)/ref(1) - 1.0_wp) <= 0.05_wp .AND. ABS(r_tf(2) - ref(2)) <= 0.005_wp, &
                 'G2: towed with Tf, the cable responds as the cable held in the current')
    CALL require(ABS(r_0(1)/ref(1) - 1.0_wp) > 0.10_wp, 'G2: towed without Tf, it responds otherwise')
  END SUBROUTINE check_tow

  SUBROUTINE check_refusals()
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('vfd_bad.dat', '1.0', '', '', 1.0_wp, .FALSE., frame_rows='riser 10')
    CALL run('vfd_bad.dat', 'vfd_bad', ok, msg)
    WRITE (*, '(A)') 'G4 a synchronisation type without a VIV AMPLITUDE table: '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'give the type a WAKE OSCILLATOR row or a VIV AMPLITUDE table') > 0, &
                 'G4: a VIV FRAME row on a synchronisation type without an amplitude table is refused')
  END SUBROUTINE check_refusals

END PROGRAM test_viv_frame_deck
