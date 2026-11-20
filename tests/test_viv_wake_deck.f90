! File: tests/test_viv_wake_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_wake_deck
  !! Wake-oscillator VIV on a deck (the optional WAKE OSCILLATOR section, TypeName CL0 [St [eps
  !! [A]]], and the VIV FRAME section, TypeName Tf; Facchinetti et al. 2004): the taut cable of the in-line VIV gate (100 m, D 0.1 m,
  !! 50 kN, pinned at 50 m depth) in a uniform 0.15 m/s current along y, its type with no
  !! synchronisation lift (Cl 0) and the row riser 0.3 (St 0.2, eps 0.3, A 12). W1-W3 run on the
  !! finite-EI cable (EI 1 N m^2) and on the EI = 0 line alike.
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
  !!      also needs), and a line that mixes the two models (finite-EI and EI = 0);
  !!   W5 the two cores agree: the EI = 0 line locks in within 3 % of the finite-EI cable's A/D
  !!      and 0.003 Hz of its frequency, and so does its uncoupled response;
  !!   W6 a finite-EI cable towed through still water: both ends move at 0.15 m/s against y (a
  !!      10 s half-cosine ramp), no current. With the filter time Tf 10 s the wake frequency takes
  !!      the flow relative to the cable's filtered motion, and the cable locks in as the fixed
  !!      cable in the 0.15 m/s current does (A/D within 5 %, frequency within 0.003 Hz); without
  !!      it (the ambient flow, zero) the wake does not lock in (A/D below a tenth of that). The
  !!      fixed cable in the current with Tf 10 s locks in as without it (the filter removes the
  !!      VIV velocity: A/D within 2 %, frequency within 0.002 Hz), on both cores, and the two
  !!      cores agree there as in W5 (an EI = 0 line cannot run between two moving points);
  !!   W7 VIV FRAME refusals by name: a malformed row, an unknown Method, Method cycle on a
  !!      wake-oscillator type, an unknown or repeated type, a negative Tf, and a type without VIV
  !!      (no Cl, no WAKE OSCILLATOR row).
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, D = 0.1_wp, F_WAKE = 0.2_wp*0.15_wp/0.1_wp
  INTEGER :: nfail = 0
  REAL(wp) :: ad_fe(4) = 0.0_wp, ad_l(4) = 0.0_wp, af_fe(2) = 0.0_wp, af_l(2) = 0.0_wp

  CALL write_current('wod_current.dat')
  CALL check_runs('1.0', 'wod', ad_fe)
  CALL check_runs('0.0', 'wod0', ad_l)
  WRITE (*, '(A,4F8.4)') 'W5: lock-in A/D, f [Hz] finite-EI, EI = 0: ', ad_fe(1:2), ad_l(1:2)
  WRITE (*, '(A,4F8.4)') 'W5: uncoupled A/D, f [Hz] finite-EI, EI = 0: ', ad_fe(3:4), ad_l(3:4)
  CALL require(ABS(ad_l(1)/ad_fe(1) - 1.0_wp) <= 0.03_wp .AND. ABS(ad_l(2) - ad_fe(2)) <= 0.003_wp, &
               'W5: the EI = 0 line locks in as the finite-EI cable')
  CALL require(ABS(ad_l(3)/ad_fe(3) - 1.0_wp) <= 0.03_wp .AND. ABS(ad_l(4) - ad_fe(4)) <= 0.003_wp, &
               'W5: and responds uncoupled as it')
  CALL write_tow('wod_tow.txt')
  CALL check_tow('wodt', ad_fe)
  CALL check_filter_fixed('1.0', 'wodf', ad_fe, af_fe)
  CALL check_filter_fixed('0.0', 'wodf0', ad_l, af_l)
  CALL require(ABS(af_l(1)/af_fe(1) - 1.0_wp) <= 0.03_wp .AND. ABS(af_l(2) - af_fe(2)) <= 0.003_wp, &
               'W6: with Tf the EI = 0 line locks in as the finite-EI cable')
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

  SUBROUTINE write_deck(path, wake_rows, tmax, cl, ei, viv_row, mixed, tow, frame_rows)
    !! wake_rows: the WAKE OSCILLATOR section's data rows, ';'-separated ('' for no section); cl the
    !! type's LINE TYPES Cl (default 0), ei its EI (default 1), viv_row a VIV section row, mixed:
    !! the line in two sections, the second of a type 'lift' with Cl 0.8; tow: no current, both
    !! ends Coupled and moved by the motion file wod_tow.txt; frame_rows: the VIV FRAME section's
    !! data rows, ';'-separated (none without it).
    CHARACTER(*), INTENT(IN) :: path, wake_rows
    REAL(wp), INTENT(IN) :: tmax
    CHARACTER(*), INTENT(IN), OPTIONAL :: cl, ei, viv_row
    LOGICAL, INTENT(IN), OPTIONAL :: mixed, tow
    CHARACTER(*), INTENT(IN), OPTIONAL :: frame_rows
    CHARACTER(32) :: cl_col, ei_col
    INTEGER :: u, i0, i1
    LOGICAL :: two, towed
    cl_col = '0.0'
    IF (PRESENT(cl)) cl_col = cl
    ei_col = '1.0'
    IF (PRESENT(ei)) ei_col = ei
    two = .FALSE.
    IF (PRESENT(mixed)) two = mixed
    towed = .FALSE.
    IF (PRESENT(tow)) towed = tow
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
    IF (towed) THEN
      WRITE (u, '(A)') '2    Coupled      100.0   0.0   -50.0  0    0     0     0'
    ELSE
      WRITE (u, '(A)') '2    Fixed        100.0   0.0   -50.0  0    0     0     0'
    END IF
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
    IF (PRESENT(frame_rows)) THEN
      WRITE (u, '(A)') '----------------------- VIV FRAME -------------------------------------------'
      WRITE (u, '(A)') 'TypeName   Tf'
      WRITE (u, '(A)') '(-)        (s)'
      i0 = 1
      DO
        i1 = INDEX(frame_rows(i0:), ';')
        IF (i1 == 0) THEN
          WRITE (u, '(A)') TRIM(frame_rows(i0:))
          EXIT
        END IF
        WRITE (u, '(A)') frame_rows(i0:i0 + i1 - 2)
        i0 = i0 + i1
      END DO
    END IF
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '0.01     dtM'
    WRITE (u, '(A)') '100.0    WtrDpth'
    WRITE (u, '(A)') '1025.0   rhoW'
    WRITE (u, '(A)') '9.80665  g'
    IF (towed) THEN
      WRITE (u, '(A)') 'wod_tow.txt  motionFile'
    ELSE
      WRITE (u, '(A)') 'wod_current.dat  WaterKin'
    END IF
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

  SUBROUTINE check_runs(ei, tag, res)
    !! W1-W3 on the cable of EI ei ('1.0' finite-EI, '0.0' the EI = 0 line), the files named
    !! after tag; res = (coupled A/D, f, uncoupled A/D, f) over the steady windows.
    CHARACTER(*), INTENT(IN) :: ei, tag
    REAL(wp), INTENT(OUT) :: res(4)
    REAL(wp), ALLOCATABLE :: tw(:, :), t0(:, :), tn(:, :), tz(:, :)
    INTEGER :: nw, n0, nn, nz
    REAL(wp) :: ad_w, f_w, ad_w2, f_w2, ad_0, f_0, pv, pd, pb
    LOGICAL :: ok, keep(20000)
    CHARACTER(1024) :: msg
    CHARACTER(:), ALLOCATABLE :: lab
    res = 0.0_wp
    lab = ' (EI '//ei//')'
    CALL write_deck(tag//'_wake.dat', 'riser 0.3', 100.0_wp, ei=ei)
    CALL run(tag//'_wake.dat', tag//'_wake', ok, msg)
    CALL require(ok, 'W1: the wake-oscillator deck runs'//lab//': '//TRIM(msg))
    CALL write_deck(tag//'_a0.dat', 'riser 0.3 0.2 0.3 0.0', 100.0_wp, ei=ei)
    CALL run(tag//'_a0.dat', tag//'_a0', ok, msg)
    CALL require(ok, 'W1: the uncoupled deck (A 0) runs'//lab//': '//TRIM(msg))
    CALL read_out(tag//'_wake.out', tw, nw)
    CALL read_out(tag//'_a0.out', t0, n0)
    CALL require(nw > 1000 .AND. n0 == nw, 'W1: the runs wrote matching rows'//lab)
    IF (.NOT. (nw > 1000 .AND. n0 == nw)) RETURN
    CALL window_stats(tw, nw, 50.0_wp, 75.0_wp, ad_w, f_w)
    CALL window_stats(tw, nw, 75.0_wp, 100.0_wp, ad_w2, f_w2)
    CALL window_stats(t0, n0, 50.0_wp, 100.0_wp, ad_0, f_0)
    res = [ad_w, f_w, ad_0, f_0]
    WRITE (*, '(A,2F8.4,A,2F8.4,A,F7.4)') 'W1'//lab//': coupled A/D, f [Hz] (50-75 s, 75-100 s) ', ad_w, f_w, &
      ', ', ad_w2, f_w2, '; wake frequency ', F_WAKE
    WRITE (*, '(A,2F8.4)') 'W1'//lab//': uncoupled (A 0) A/D, f [Hz] ', ad_0, f_0
    CALL require(f_w < F_WAKE - 0.01_wp .AND. ad_w > 0.4_wp .AND. ad_w < 0.8_wp, &
                 'W1: the coupled response locks in below the wake frequency at A/D 0.4-0.8'//lab)
    CALL require(ABS(ad_w2/ad_w - 1.0_wp) <= 0.01_wp .AND. ABS(f_w2 - f_w) <= 0.002_wp, 'W1: and is steady'//lab)
    CALL require(ABS(f_0 - F_WAKE) <= 0.005_wp, 'W1: uncoupled, the wake drives the cable at its own frequency'//lab)
    CALL require(ad_w >= 3.0_wp*ad_0, 'W1: lock-in at least triples the uncoupled amplitude'//lab)
    ! W2: the power balance of the steady lock-in
    keep(1:nw) = tw(1, 1:nw) >= 60.0_wp
    pv = SUM(PACK(tw(3, 1:nw), keep(1:nw)))/COUNT(keep(1:nw))
    pd = SUM(PACK(tw(4, 1:nw), keep(1:nw)))/COUNT(keep(1:nw))
    pb = SUM(PACK(tw(5, 1:nw), keep(1:nw)))/COUNT(keep(1:nw))
    WRITE (*, '(A,3ES12.4)') 'W2'//lab//': mean VIVPow, DragPow, DampPow [W] over 60-100 s: ', pv, pd, pb
    CALL require(pv > 0.0_wp .AND. ABS(pv + pd + pb) <= 0.005_wp*pv, &
                 'W2: the drag and damping take out the wake-lift power within 0.5 %'//lab)
    ! W3: CL0 0 is no wake oscillator
    CALL write_deck(tag//'_none.dat', '', 20.0_wp, ei=ei)
    CALL run(tag//'_none.dat', tag//'_none', ok, msg)
    CALL require(ok, 'W3: the deck without the section runs'//lab//': '//TRIM(msg))
    CALL write_deck(tag//'_zero.dat', 'riser 0.0', 20.0_wp, ei=ei)
    CALL run(tag//'_zero.dat', tag//'_zero', ok, msg)
    CALL require(ok, 'W3: the CL0 0 deck runs'//lab//': '//TRIM(msg))
    CALL read_out(tag//'_none.out', tn, nn)
    CALL read_out(tag//'_zero.out', tz, nz)
    CALL require(nn > 100 .AND. nz == nn, 'W3: the runs wrote matching rows'//lab)
    IF (nn > 100 .AND. nz == nn) CALL require(ALL(tz(:, 1:nz) == tn(:, 1:nn)), &
                                              'W3: a row with CL0 0 is the deck without the section'//lab)
  END SUBROUTINE check_runs

  SUBROUTINE write_tow(path)
    !! Both ends (points 1 and 2) moved along -y at 0.15 m/s after a 10 s half-cosine velocity ramp,
    !! every 0.01 s to 100 s: y = -U (t/2 - T sin(pi t/T)/(2 pi)) in the ramp (T 10 s).
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), PARAMETER :: U = 0.15_wp, TR = 10.0_wp
    INTEGER :: u_, k
    REAL(wp) :: t, y, v, a
    OPEN (NEWUNIT=u_, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u_, '(A)') '# time point x y z vx vy vz ax ay az'
    DO k = 0, 10000
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

  SUBROUTINE check_tow(tag, ref)
    !! W6 on the finite-EI cable: towed with Tf 10 s against the fixed cable in the current (ref,
    !! from check_runs), and without Tf.
    CHARACTER(*), INTENT(IN) :: tag
    REAL(wp), INTENT(IN) :: ref(4)
    CHARACTER(*), PARAMETER :: ei = '1.0'
    REAL(wp), ALLOCATABLE :: tf(:, :), t0(:, :)
    INTEGER :: nf, n0
    REAL(wp) :: ad_f, f_f, ad_0, f_0
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CHARACTER(:), ALLOCATABLE :: lab
    lab = ' (EI '//ei//')'
    CALL write_deck(tag//'_tf.dat', 'riser 0.3', 100.0_wp, ei=ei, tow=.TRUE., frame_rows='riser 10')
    CALL run(tag//'_tf.dat', tag//'_tf', ok, msg)
    CALL require(ok, 'W6: the towed deck with Tf runs'//lab//': '//TRIM(msg))
    CALL write_deck(tag//'_amb.dat', 'riser 0.3', 100.0_wp, ei=ei, tow=.TRUE.)
    CALL run(tag//'_amb.dat', tag//'_amb', ok, msg)
    CALL require(ok, 'W6: the towed deck without Tf runs'//lab//': '//TRIM(msg))
    CALL read_out(tag//'_tf.out', tf, nf)
    CALL read_out(tag//'_amb.out', t0, n0)
    CALL require(nf > 1000 .AND. n0 == nf, 'W6: the towed runs wrote matching rows'//lab)
    IF (.NOT. (nf > 1000 .AND. n0 == nf)) RETURN
    CALL window_stats(tf, nf, 50.0_wp, 75.0_wp, ad_f, f_f)
    CALL window_stats(t0, n0, 50.0_wp, 75.0_wp, ad_0, f_0)
    WRITE (*, '(A,2F8.4,A,2F8.4,A,F8.4)') 'W6'//lab//': towed A/D, f [Hz] with Tf ', ad_f, f_f, &
      '; in the current ', ref(1:2), '; towed without Tf A/D ', ad_0
    CALL require(ABS(ad_f/ref(1) - 1.0_wp) <= 0.05_wp .AND. ABS(f_f - ref(2)) <= 0.003_wp, &
                 'W6: towed with Tf, the cable locks in as the fixed cable in the current'//lab)
    CALL require(ad_0 < 0.1_wp*ref(1), 'W6: towed without Tf (the ambient flow), the wake does not lock in'//lab)
  END SUBROUTINE check_tow

  SUBROUTINE check_filter_fixed(ei, tag, ref, res)
    !! W6 on the fixed cable of EI ei in the current: with Tf 10 s against without it (ref, from
    !! check_runs); res = (A/D, f) with Tf.
    CHARACTER(*), INTENT(IN) :: ei, tag
    REAL(wp), INTENT(IN) :: ref(4)
    REAL(wp), INTENT(OUT) :: res(2)
    REAL(wp), ALLOCATABLE :: tf(:, :)
    INTEGER :: nf
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CHARACTER(:), ALLOCATABLE :: lab
    lab = ' (EI '//ei//')'
    res = 0.0_wp
    CALL write_deck(tag//'.dat', 'riser 0.3', 100.0_wp, ei=ei, frame_rows='riser 10')
    CALL run(tag//'.dat', tag, ok, msg)
    CALL require(ok, 'W6: the fixed cable with Tf runs'//lab//': '//TRIM(msg))
    CALL read_out(tag//'.out', tf, nf)
    CALL require(nf > 1000, 'W6: the fixed run with Tf wrote its rows'//lab)
    IF (nf <= 1000) RETURN
    CALL window_stats(tf, nf, 50.0_wp, 75.0_wp, res(1), res(2))
    WRITE (*, '(A,2F8.4,A,2F8.4)') 'W6'//lab//': fixed in the current A/D, f [Hz] with Tf ', res, '; without ', &
      ref(1:2)
    CALL require(ABS(res(1)/ref(1) - 1.0_wp) <= 0.02_wp .AND. ABS(res(2) - ref(2)) <= 0.002_wp, &
                 'W6: fixed in the current, Tf leaves the lock-in as the ambient flow gives it'//lab)
  END SUBROUTINE check_filter_fixed

  SUBROUTINE check_refusals()
    CALL refused('riser 0.3 0.2 0.3 12 1', 'TypeName CL0 [St [eps [A]]]', 'a row of six columns')
    CALL refused('riser 0.3', 'VIV FRAME rows must be: TypeName Tf [Method]', 'W7: a VIV FRAME row of four columns', &
                 frame_rows='riser 10 cycle 1')
    CALL refused('riser 0.3', 'VIV FRAME Method must be filter or cycle', &
                 'W7: a VIV FRAME row with an unknown Method', frame_rows='riser 10 1')
    CALL refused('riser 0.3', 'acts on the synchronisation model', 'W7: Method cycle on a wake-oscillator type', &
                 frame_rows='riser 10 cycle')
    CALL refused('riser 0.3', 'VIV FRAME row names undefined LINE TYPE', 'W7: an unknown type', &
                 frame_rows='cable 10')
    CALL refused('riser 0.3', 'VIV FRAME names LINE TYPE', 'W7: a repeated type', frame_rows='riser 10;riser 5')
    CALL refused('riser 0.3', 'filter time Tf must be >= 0', 'W7: a negative Tf', frame_rows='riser -1')
    CALL refused('', 'give the type LINE TYPES Cl > 0 or a WAKE OSCILLATOR row', 'W7: a type without VIV', &
                 frame_rows='riser 10')
    CALL refused('riser 0.3 x', 'plain numbers', 'a non-numeric St')
    CALL refused('cable 0.3', 'undefined LINE TYPE', 'an unknown type')
    CALL refused('riser 0.3;riser 0.2', 'twice', 'a repeated type')
    CALL refused('riser -0.1', 'CL0 must be non-negative', 'a negative CL0')
    CALL refused('riser 0.3 0.0', 'St > 0, eps > 0 and A >= 0', 'St = 0')
    CALL refused('riser 0.3 0.2 0.0', 'St > 0, eps > 0 and A >= 0', 'eps = 0')
    CALL refused('riser 0.3 0.2 0.3 -1', 'St > 0, eps > 0 and A >= 0', 'a negative A')
    CALL refused('riser 0.3', 'give the type LINE TYPES Cl 0', 'a type that keeps its lift', cl='0.8')
    CALL refused('riser 0.3', 'a line takes one VIV model', 'a line mixing the two models', mixed=.TRUE.)
    CALL refused('riser 0.3', 'a line takes one VIV model', 'an EI = 0 line mixing the two models', ei='0.0', &
                 mixed=.TRUE.)
    CALL refused('riser 0.3', 'give the type LINE TYPES Cl 0', 'an EI = 0 type that keeps its lift', cl='0.8', &
                 ei='0.0')
  END SUBROUTINE check_refusals

  SUBROUTINE refused(rows, needle, what, cl, ei, viv_row, mixed, frame_rows)
    CHARACTER(*), INTENT(IN) :: rows, needle, what
    CHARACTER(*), INTENT(IN), OPTIONAL :: cl, ei, viv_row, frame_rows
    LOGICAL, INTENT(IN), OPTIONAL :: mixed
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    ! an absent optional passes on absent
    CALL write_deck('wod_bad.dat', rows, 1.0_wp, cl=cl, ei=ei, viv_row=viv_row, mixed=mixed, frame_rows=frame_rows)
    CALL run('wod_bad.dat', 'wod_bad', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, needle) > 0, 'W4: '//what//' is refused: '//TRIM(msg))
  END SUBROUTINE refused

END PROGRAM test_viv_wake_deck
