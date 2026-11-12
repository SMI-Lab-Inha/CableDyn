! File: tests/test_viv_inline_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_inline_deck
  !! In-line VIV on a deck (the optional VIV section, TypeName ClX [alpha]): the taut line of the
  !! cross-flow lock-in benchmark (tests/test_viv_lockin.py) as a finite-EI cable, 100 m pinned
  !! at 50 m depth in a uniform 0.15 m/s current along y, Cl 0.8 with MoorDyn's band.
  !!   D1 with ClX 0.85 the quarter-span in-line (y) response oscillates at twice the mid-span
  !!      cross-flow (z) frequency and is larger than without the VIV section;
  !!   D2 a VIV row with ClX 0 is the deck without the section, to the printed digits;
  !!   D3 malformed rows, an unknown or repeated type, ClX on a type without the cross-flow lift
  !!      and alpha outside [0, 1) are refused by name;
  !!   D4 the same deck as an EI = 0 line (the lumped node force, neighbour-averaged kinematics)
  !!      matches the finite-EI cable: cross-flow and in-line amplitudes within 2 %, the same
  !!      frequencies within one 0.01 Hz bin;
  !!   D5 a higher-harmonic coefficient Chh 1.0 (Ulveseter et al. 2019, the fourth column) at
  !!      least triples the quarter-span in-line response at four times the cross-flow frequency
  !!      (the in-line share of the force at the cross-flow plus in-line phase), on the cable and
  !!      the EI = 0 line alike (within 5 %); Chh 0 is the row without it, to the printed digits;
  !!      Chh on a type without in-line VIV and a negative Chh are refused;
  !!   D6 the cross-flow columns after the 14-column finite-EI row (GAs GJ Irt Irn): a 17-column row
  !!      free in torsion is the 13-column deck, to the printed digits; a 15-column row (Cl alone) is
  !!      the 17-column one with MoorDyn's dF 0.08 and cF 0.18; and on the line restrained in
  !!      torsion at both ends (nodal twist) the lift still locks the cable in, with the cross-flow
  !!      amplitude and frequency of the line free in torsion.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, D = 0.1_wp
  INTEGER :: nfail = 0

  CALL write_current('vid_current.dat')
  CALL check_runs()
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: in-line VIV deck gates'

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
    WRITE (u, '(A)') 'uniform current along y for the in-line VIV deck gate'
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

  SUBROUTINE write_deck(path, viv_rows, tmax, ei, cl, wide, ends)
    !! viv_rows: the VIV section's data rows, ';'-separated ('' for no section). wide: the VIV
    !! columns written after the 14-column finite-EI row (GAs GJ Irt Irn in the CableDyn hydro
    !! order) instead of the 10 stock-order ones; ends: restrain both ends in torsion.
    CHARACTER(*), INTENT(IN) :: path, viv_rows
    REAL(wp), INTENT(IN) :: tmax
    CHARACTER(*), INTENT(IN), OPTIONAL :: ei, cl, wide
    LOGICAL, INTENT(IN), OPTIONAL :: ends
    CHARACTER(32) :: ei_col, cl_col
    INTEGER :: u, i0, i1
    ei_col = '1.0'
    IF (PRESENT(ei)) ei_col = ei
    cl_col = '0.8'
    IF (PRESENT(cl)) cl_col = cl
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (u, '(A)') 'in-line VIV deck gate: taut cable in a uniform current'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    IF (PRESENT(wide)) THEN
      ! GJ = EI/1.3, Irn = m D^2/8, Irt = m D^2/16, GAs = EA/2.6: the circular-section closure
      WRITE (u, '(A)') 'TypeName Diam Mass/m EA BA/-zeta EI GAs GJ Irt Irn Cd_n Cd_t Ca_n Ca_t Cl [dF cF]'
      WRITE (u, '(A)') 'riser 0.10 10.0 5.0e7 -0.8 '//TRIM(ei_col)//' 1.9230769e7 0.76923077 0.00625 0.0125 '// &
        '1.2 0.0 1.0 0.0 '//wide
    ELSE
      WRITE (u, '(A)') 'TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF'
      WRITE (u, '(A)') '(name)     (m)    (kg/m)   (N)       (N-s/-)   (-)  (-)   (-)   (-)   (-)   (-)  (-)   (-)'
      WRITE (u, '(A)') 'riser      0.10   10.0     5.0e7     -0.8      '//TRIM(ei_col)//'  1.2   1.0   0.0   0.0   '// &
        TRIM(cl_col)//'  0.08  0.18'
    END IF
    WRITE (u, '(A)') '----------------------- POINTS ----------------------------------------------'
    WRITE (u, '(A)') 'ID   Attachment   X       Y     Z      M    V     CdA   CA'
    WRITE (u, '(A)') '(#)  (-)          (m)     (m)   (m)    (kg) (m^3) (m^2) (-)'
    WRITE (u, '(A)') '1    Coupled      0.0     0.0   -50.0  0    0     0     0'
    WRITE (u, '(A)') '2    Fixed        100.0   0.0   -50.0  0    0     0     0'
    WRITE (u, '(A)') '----------------------- LINES -----------------------------------------------'
    WRITE (u, '(A)') 'ID   LineType   AttachA  AttachB  UnstrLen  NumSegs  Outputs'
    WRITE (u, '(A)') '(#)  (name)     (#)      (#)      (m)       (-)      (-)'
    WRITE (u, '(A)') '1    riser      1        2        99.900100  50       -'
    IF (PRESENT(ends)) THEN
      IF (ends) THEN
        ! pinned in bending, rigid in torsion, the zero-twist normal along z
        WRITE (u, '(A)') '----------------------- END CONNECTIONS ---------------------------------------'
        WRITE (u, '(A)') 'LineID End Stiffness EzX EzY EzZ TorsStiffness NxX NxY NxZ Pretwist'
        WRITE (u, '(A)') '1      A   0         1.0 0.0 0.0 Rigid         0.0 0.0 1.0 0.0'
        WRITE (u, '(A)') '1      B   0         1.0 0.0 0.0 Rigid         0.0 0.0 1.0 0.0'
      END IF
    END IF
    IF (LEN_TRIM(viv_rows) > 0) THEN
      WRITE (u, '(A)') '----------------------- VIV -------------------------------------------------'
      WRITE (u, '(A)') 'TypeName   ClX    alpha'
      WRITE (u, '(A)') '(name)     (-)    (-)'
      i0 = 1
      DO
        i1 = INDEX(viv_rows(i0:), ';')
        IF (i1 == 0) THEN
          WRITE (u, '(A)') TRIM(viv_rows(i0:))
          EXIT
        END IF
        WRITE (u, '(A)') viv_rows(i0:i0 + i1 - 2)
        i0 = i0 + i1
      END DO
    END IF
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '0.01     dtM'
    WRITE (u, '(A)') '100.0    WtrDpth'
    WRITE (u, '(A)') '1025.0   rhoW'
    WRITE (u, '(A)') '9.80665  g'
    WRITE (u, '(A)') 'vid_current.dat  WaterKin'
    WRITE (u, '(A)') '0.25     VIVrmsWindow'
    WRITE (u, '(A,F8.3,A)') ' ', tmax, '    TMax'
    WRITE (u, '(A)') '---------------------- OUTPUTS -----------------------------------------------'
    WRITE (u, '(A)') 'L1N25PZ'
    WRITE (u, '(A)') 'L1N12PY'
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
    !! The rows (time, L1N25PZ, L1N12PY) of a deck's main output file.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: tab(:, :)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(1024) :: rec
    INTEGER :: u, ios
    ALLOCATE (tab(3, 20000))
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
    !! The frequency of the largest Fourier amplitude of x - mean(x) on 0.05-2 Hz (0.002 Hz grid).
    REAL(wp), INTENT(IN) :: tt(:), x(:)
    REAL(wp) :: xm(SIZE(x)), f, best, re, im
    INTEGER :: i, j
    xm = x - SUM(x)/REAL(SIZE(x), wp)
    best = -1.0_wp
    fpk = 0.0_wp
    DO j = 0, 975
      f = 0.05_wp + 0.002_wp*REAL(j, wp)
      re = 0.0_wp; im = 0.0_wp
      DO i = 1, SIZE(x)
        re = re + xm(i)*COS(2.0_wp*PI*f*tt(i))
        im = im + xm(i)*SIN(2.0_wp*PI*f*tt(i))
      END DO
      IF (re*re + im*im > best) THEN
        best = re*re + im*im
        fpk = f
      END IF
    END DO
  END FUNCTION peak_frequency

  SUBROUTINE window(tab, n, t0, tt, z, y)
    !! The rows of tab from time t0 on.
    REAL(wp), INTENT(IN) :: tab(:, :), t0
    INTEGER, INTENT(IN) :: n
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: tt(:), z(:), y(:)
    LOGICAL :: keep(n)
    keep = tab(1, 1:n) >= t0
    tt = PACK(tab(1, 1:n), keep)
    z = PACK(tab(2, 1:n), keep)
    y = PACK(tab(3, 1:n), keep)
  END SUBROUTINE window

  SUBROUTINE check_runs()
    REAL(wp), ALLOCATABLE :: ti(:, :), tc(:, :), t0(:, :), tt(:), z(:), y(:)
    INTEGER :: ni, nc, n0
    REAL(wp) :: fz, fy, ay_il, ay_cf
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('vid_il.dat', 'riser 0.85 0.15', 100.0_wp)
    CALL run('vid_il.dat', 'vid_il', ok, msg)
    CALL require(ok, 'D1: the in-line VIV deck runs: '//TRIM(msg))
    CALL write_deck('vid_cf.dat', '', 100.0_wp)
    CALL run('vid_cf.dat', 'vid_cf', ok, msg)
    CALL require(ok, 'D1: the cross-flow-only deck runs: '//TRIM(msg))
    CALL write_deck('vid_zero.dat', 'riser 0.0', 100.0_wp)
    CALL run('vid_zero.dat', 'vid_zero', ok, msg)
    CALL require(ok, 'D2: the ClX 0 deck runs: '//TRIM(msg))
    CALL read_out('vid_il.out', ti, ni)
    CALL read_out('vid_cf.out', tc, nc)
    CALL read_out('vid_zero.out', t0, n0)
    CALL require(ni > 100 .AND. ni == nc .AND. ni == n0, 'D1/D2: the runs wrote matching rows')
    IF (.NOT. (ni > 100 .AND. ni == nc .AND. ni == n0)) RETURN
    CALL require(ALL(t0(:, 1:n0) == tc(:, 1:nc)), 'D2: a VIV row with ClX 0 is the deck without the section')
    CALL window(ti, ni, 50.0_wp, tt, z, y)
    fz = peak_frequency(tt, z)
    fy = peak_frequency(tt, y)
    ay_il = 0.5_wp*(MAXVAL(y) - MINVAL(y))
    WRITE (*, '(A,F7.4,A,F7.4,A,F7.3,A,F7.4)') 'D1: cross-flow ', fz, ' Hz, in-line ', fy, ' Hz, ratio ', fy/fz, &
      ', CF amplitude/D ', 0.5_wp*(MAXVAL(z) - MINVAL(z))/D
    CALL window(tc, nc, 50.0_wp, tt, z, y)
    ay_cf = 0.5_wp*(MAXVAL(y) - MINVAL(y))
    WRITE (*, '(A,F8.4,A,F8.4)') 'D1: quarter-span in-line amplitude/D ', ay_il/D, ', without the VIV section ', &
      ay_cf/D
    CALL require(ABS(fy/fz - 2.0_wp) <= 0.15_wp, 'D1: the in-line response is at twice the cross-flow frequency')
    CALL require(ay_il > 2.0_wp*ay_cf, 'D1: the in-line force raises the in-line response')
    ! D4: the EI = 0 line
    BLOCK
      REAL(wp), ALLOCATABLE :: te(:, :)
      INTEGER :: ne0
      REAL(wp) :: az_f, az_0, ay_f, ay_0, fz_f, fy_f
      CALL window(ti, ni, 50.0_wp, tt, z, y)
      az_f = 0.5_wp*(MAXVAL(z) - MINVAL(z))
      ay_f = 0.5_wp*(MAXVAL(y) - MINVAL(y))
      fz_f = peak_frequency(tt, z)
      fy_f = peak_frequency(tt, y)
      CALL write_deck('vid_ei0.dat', 'riser 0.85 0.15', 100.0_wp, ei='0.0')
      CALL run('vid_ei0.dat', 'vid_ei0', ok, msg)
      CALL require(ok, 'D4: the EI = 0 in-line VIV deck runs: '//TRIM(msg))
      CALL read_out('vid_ei0.out', te, ne0)
      CALL require(ne0 == ni, 'D4: the EI = 0 run wrote the same rows')
      IF (ne0 /= ni) RETURN
      CALL window(te, ne0, 50.0_wp, tt, z, y)
      az_0 = 0.5_wp*(MAXVAL(z) - MINVAL(z))
      ay_0 = 0.5_wp*(MAXVAL(y) - MINVAL(y))
      WRITE (*, '(A,2F8.4,A,2F8.4)') 'D4: CF and IL amplitude/D, EI = 0 line ', az_0/D, ay_0/D, &
        ', finite-EI cable ', az_f/D, ay_f/D
      WRITE (*, '(A,2F8.4,A,2F8.4)') 'D4: CF and IL frequency [Hz], EI = 0 line ', peak_frequency(tt, z), &
        peak_frequency(tt, y), ', finite-EI cable ', fz_f, fy_f
      CALL require(ABS(az_0/az_f - 1.0_wp) <= 0.02_wp .AND. ABS(ay_0/ay_f - 1.0_wp) <= 0.02_wp, &
                   'D4: the EI = 0 line matches the cable''s amplitudes within 2 %')
      CALL require(ABS(peak_frequency(tt, z) - fz_f) <= 0.01_wp .AND. ABS(peak_frequency(tt, y) - fy_f) <= 0.01_wp, &
                   'D4: and its frequencies within 0.01 Hz')
    END BLOCK
    ! D5: the higher-harmonic force
    BLOCK
      REAL(wp), ALLOCATABLE :: th(:, :), t0h(:, :), te(:, :)
      INTEGER :: nh, n0h, neh
      REAL(wp) :: f1, a4_il, a4_hh, a4_e
      CALL write_deck('vid_hh.dat', 'riser 0.85 0.15 1.0', 100.0_wp)
      CALL run('vid_hh.dat', 'vid_hh', ok, msg)
      CALL require(ok, 'D5: the Chh deck runs: '//TRIM(msg))
      CALL write_deck('vid_hh0.dat', 'riser 0.85 0.15 0.0', 100.0_wp)
      CALL run('vid_hh0.dat', 'vid_hh0', ok, msg)
      CALL require(ok, 'D5: the Chh 0 deck runs: '//TRIM(msg))
      CALL write_deck('vid_hhe.dat', 'riser 0.85 0.15 1.0', 100.0_wp, ei='0.0')
      CALL run('vid_hhe.dat', 'vid_hhe', ok, msg)
      CALL require(ok, 'D5: the EI = 0 Chh deck runs: '//TRIM(msg))
      CALL read_out('vid_hh.out', th, nh)
      CALL read_out('vid_hh0.out', t0h, n0h)
      CALL read_out('vid_hhe.out', te, neh)
      CALL require(nh == ni .AND. n0h == ni .AND. neh == ni, 'D5: the runs wrote the same rows')
      IF (.NOT. (nh == ni .AND. n0h == ni .AND. neh == ni)) RETURN
      CALL require(ALL(t0h(:, 1:n0h) == ti(:, 1:ni)), 'D5: a VIV row with Chh 0 is the row without it')
      CALL window(ti, ni, 50.0_wp, tt, z, y)
      f1 = peak_frequency(tt, z)
      a4_il = fourier_amplitude(tt, y, 4.0_wp*f1)
      CALL window(th, nh, 50.0_wp, tt, z, y)
      f1 = peak_frequency(tt, z)
      a4_hh = fourier_amplitude(tt, y, 4.0_wp*f1)
      CALL window(te, neh, 50.0_wp, tt, z, y)
      a4_e = fourier_amplitude(tt, y, 4.0_wp*peak_frequency(tt, z))
      WRITE (*, '(A,3F8.4)') 'D5: quarter-span in-line amplitude/D at 4 f, Chh 0, Chh 1, Chh 1 EI = 0: ', &
        a4_il/D, a4_hh/D, a4_e/D
      CALL require(a4_hh >= 3.0_wp*a4_il, 'D5: Chh raises the fourth in-line harmonic at least threefold')
      CALL require(ABS(a4_e/a4_hh - 1.0_wp) <= 0.05_wp, 'D5: the EI = 0 line takes it within 5 %')
    END BLOCK
    CALL check_wide_rows(tc, nc)
  END SUBROUTINE check_runs

  SUBROUTINE check_wide_rows(tc, nc)
    !! D6 against the cross-flow-only deck tc (the 13-column row, no VIV section).
    REAL(wp), INTENT(IN) :: tc(:, :)
    INTEGER, INTENT(IN) :: nc
    REAL(wp), ALLOCATABLE :: tw(:, :), t15(:, :), t17(:, :), tr(:, :), tt(:), z(:), y(:)
    INTEGER :: nw, n15, n17, nr
    REAL(wp) :: az_c, fz_c, az_r, fz_r
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('vid_w13.dat', '', 100.0_wp, wide='0.8 0.08 0.18')
    CALL run('vid_w13.dat', 'vid_w13', ok, msg)
    CALL require(ok, 'D6: the 17-column deck runs: '//TRIM(msg))
    CALL read_out('vid_w13.out', tw, nw)
    CALL require(nw == nc, 'D6: the 17-column run wrote the rows of the 13-column one')
    IF (nw /= nc) RETURN
    CALL require(ALL(tw(:, 1:nw) == tc(:, 1:nc)), 'D6: a 17-column row free in torsion is the 13-column row')
    CALL write_deck('vid_w15.dat', '', 30.0_wp, wide='0.8', ends=.TRUE.)
    CALL run('vid_w15.dat', 'vid_w15', ok, msg)
    CALL require(ok, 'D6: the 15-column deck restrained in torsion runs: '//TRIM(msg))
    CALL write_deck('vid_w17.dat', '', 30.0_wp, wide='0.8 0.08 0.18', ends=.TRUE.)
    CALL run('vid_w17.dat', 'vid_w17', ok, msg)
    CALL require(ok, 'D6: the 17-column deck restrained in torsion runs: '//TRIM(msg))
    CALL read_out('vid_w15.out', t15, n15)
    CALL read_out('vid_w17.out', t17, n17)
    CALL require(n15 > 100 .AND. n15 == n17, 'D6: the 15- and 17-column runs wrote matching rows')
    IF (.NOT. (n15 > 100 .AND. n15 == n17)) RETURN
    CALL require(ALL(t15(:, 1:n15) == t17(:, 1:n17)), 'D6: a 15-column row takes dF 0.08 and cF 0.18')
    CALL write_deck('vid_wr.dat', '', 100.0_wp, wide='0.8 0.08 0.18', ends=.TRUE.)
    CALL run('vid_wr.dat', 'vid_wr', ok, msg)
    CALL require(ok, 'D6: the VIV deck restrained in torsion runs: '//TRIM(msg))
    CALL read_out('vid_wr.out', tr, nr)
    CALL require(nr == nc, 'D6: the run restrained in torsion wrote the same rows')
    IF (nr /= nc) RETURN
    CALL window(tc, nc, 50.0_wp, tt, z, y)
    az_c = 0.5_wp*(MAXVAL(z) - MINVAL(z))
    fz_c = peak_frequency(tt, z)
    CALL window(tr, nr, 50.0_wp, tt, z, y)
    az_r = 0.5_wp*(MAXVAL(z) - MINVAL(z))
    fz_r = peak_frequency(tt, z)
    WRITE (*, '(A,2F8.4,A,2F8.4)') 'D6: mid-span CF amplitude/D and frequency [Hz], restrained in torsion ', &
      az_r/D, fz_r, ', free ', az_c/D, fz_c
    CALL require(az_r > 0.3_wp*D, 'D6: the lift locks in the cable restrained in torsion')
    CALL require(ABS(az_r/az_c - 1.0_wp) <= 0.02_wp .AND. ABS(fz_r - fz_c) <= 0.01_wp, &
                 'D6: with the amplitude (2 %) and frequency (0.01 Hz) of the cable free in torsion')
  END SUBROUTINE check_wide_rows

  REAL(wp) FUNCTION fourier_amplitude(tt, x, f) RESULT(a)
    !! The amplitude of x - mean(x) at frequency f (Hann window, so a sine of amplitude A gives A).
    REAL(wp), INTENT(IN) :: tt(:), x(:), f
    REAL(wp) :: w(SIZE(x)), xm(SIZE(x)), re, im
    INTEGER :: i, n
    n = SIZE(x)
    DO i = 1, n
      w(i) = 0.5_wp - 0.5_wp*COS(2.0_wp*PI*REAL(i - 1, wp)/REAL(n - 1, wp))
    END DO
    xm = x - SUM(x)/REAL(n, wp)
    re = SUM(w*xm*COS(2.0_wp*PI*f*tt))
    im = SUM(w*xm*SIN(2.0_wp*PI*f*tt))
    a = 2.0_wp*SQRT(re*re + im*im)/SUM(w)
  END FUNCTION fourier_amplitude

  SUBROUTINE check_refusals()
    CALL refused('riser 0.85 0.15 7 1', 'TypeName ClX [alpha [Chh]]', 'a row of five columns')
    CALL refused('cable 0.85', 'undefined LINE TYPE', 'an unknown type')
    CALL refused('riser 0.85;riser 0.5', 'twice', 'a repeated type')
    CALL refused('riser -0.1', 'non-negative', 'a negative ClX')
    CALL refused('riser 0.85 1.0', '0 <= alpha < 1', 'alpha = 1')
    CALL refused('riser 0.85 x', 'plain numbers', 'a non-numeric alpha')
    CALL refused('riser 0.85', 'cross-flow lift', 'ClX on a type without the cross-flow lift', cl='0.0')
    CALL refused('riser 0.0 0.15 1.0', 'needs its in-line VIV', 'Chh without in-line VIV')
    CALL refused('riser 0.85 0.15 -1.0', 'Chh must be non-negative', 'a negative Chh')
  END SUBROUTINE check_refusals

  SUBROUTINE refused(rows, needle, what, ei, cl)
    CHARACTER(*), INTENT(IN) :: rows, needle, what
    CHARACTER(*), INTENT(IN), OPTIONAL :: ei, cl
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    IF (PRESENT(ei)) THEN
      CALL write_deck('vid_bad.dat', rows, 1.0_wp, ei=ei)
    ELSE IF (PRESENT(cl)) THEN
      CALL write_deck('vid_bad.dat', rows, 1.0_wp, cl=cl)
    ELSE
      CALL write_deck('vid_bad.dat', rows, 1.0_wp)
    END IF
    CALL run('vid_bad.dat', 'vid_bad', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, needle) > 0, 'D3: '//what//' is refused: '//TRIM(msg))
  END SUBROUTINE refused

END PROGRAM test_viv_inline_deck
