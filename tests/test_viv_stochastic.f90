! File: tests/test_viv_stochastic.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_viv_stochastic
  !! The stochastic band centre of the synchronisation model (VIV STOCHASTIC section; Ulveseter et
  !! al. 2017, Eq. 6, Fig. 4): the centre cF of a node becomes cF + sum_i a sin(p_i + 2 pi f_i
  !! theta), theta = integral U/D dt, of N_f components of a flat spectrum up to c cF.
  !!   S1 the process (CD_VIV_Sto_Node, CD_VIV_Sto_Centre) at their NDP values (cF 0.144, sigma 0.016,
  !!      c 0.07, 100 components) over theta 0-4e6 (about 200 periods of the lowest component)
  !!      sampled every 5: mean within 0.03 sigma of cF, standard deviation within 3 % of sigma (the
  !!      components' amplitudes give exactly sigma), excess kurtosis within 0.3 of a Gaussian's 0;
  !!      no component above c cF; the start phases repeat for the same seed and node and differ
  !!      between nodes and between seeds;
  !!   S2 a sigma 0 row is the deck without the section, bit for bit, on the finite-EI cable and on
  !!      the EI = 0 line (the taut cable of the VIV FRAME gates in a 0.15 m/s current);
  !!   S3 with sigma 0.02 (c 0.07, 100 components) the run repeats bit for bit, differs from the
  !!      constant centre and between seeds, and stays locked in (A/D over 50-80 s within 15 % of the
  !!      constant centre's), on both cores;
  !!   S4 refusals by name: a row on a type without Cl, on a wake-oscillator type, with c out of
  !!      (0, 1], with a non-integer Nf and with a seed out of range.
  !! The checkpoint restart and rewind with the centre are in test_viv_coupled.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_VIV, ONLY: CD_VIV_Sto_Node, CD_VIV_Sto_Centre
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, D = 0.1_wp
  REAL(wp), PARAMETER :: TMAX = 80.0_wp, T0 = 50.0_wp
  INTEGER :: nfail = 0

  CALL check_process()
  CALL write_current('vsd_current.dat')
  CALL check_deck('1.0', 'vsd_fe')
  CALL check_deck('0.0', 'vsd_l')
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: VIV STOCHASTIC gates'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE check_process()
    REAL(wp), PARAMETER :: CF = 0.144_wp, SIG = 0.016_wp, C = 0.07_wp
    INTEGER, PARAMETER :: NF = 100, NS = 800000
    REAL(wp) :: a, f(NF), p(NF), a2, f2(NF), p2(NF), x, s1, s2, s4, mean, sd, kurt
    INTEGER :: k
    CALL CD_VIV_Sto_Node(SIG, C, CF, NF, 7, 12, a, f, p)
    s1 = 0.0_wp
    s2 = 0.0_wp
    s4 = 0.0_wp
    DO k = 0, NS - 1
      x = CD_VIV_Sto_Centre(CF, a, f, p, 5.0_wp*REAL(k, wp)) - CF
      s1 = s1 + x
    END DO
    mean = s1/REAL(NS, wp)
    DO k = 0, NS - 1
      x = CD_VIV_Sto_Centre(CF, a, f, p, 5.0_wp*REAL(k, wp)) - CF - mean
      s2 = s2 + x*x
      s4 = s4 + x**4
    END DO
    sd = SQRT(s2/REAL(NS, wp))
    kurt = (s4/REAL(NS, wp))/(s2/REAL(NS, wp))**2 - 3.0_wp
    WRITE (*, '(A,3ES12.4)') 'S1 mean - cF, sd/sigma, excess kurtosis: ', mean, sd/SIG, kurt
    CALL require(ABS(mean) <= 0.03_wp*SIG, 'S1: the process has the mean cF')
    CALL require(ABS(sd/SIG - 1.0_wp) <= 0.03_wp, 'S1: the process has the standard deviation sigma')
    CALL require(ABS(kurt) <= 0.3_wp, 'S1: the process is close to Gaussian')
    CALL require(ABS(a*a*REAL(NF, wp)/2.0_wp - SIG*SIG) <= 1.0e-15_wp, 'S1: the amplitudes give sigma^2')
    CALL require(MAXVAL(f) < C*CF .AND. MINVAL(f) > 0.0_wp, 'S1: the spectrum lies in (0, c cF]')
    CALL require(ALL(p >= 0.0_wp .AND. p < 2.0_wp*PI), 'S1: the start phases lie in [0, 2 pi)')
    CALL CD_VIV_Sto_Node(SIG, C, CF, NF, 7, 12, a2, f2, p2)
    CALL require(ALL(.NOT. (ABS(p2 - p) > 0.0_wp)), 'S1: the same seed and node repeat the phases')
    CALL CD_VIV_Sto_Node(SIG, C, CF, NF, 7, 13, a2, f2, p2)
    CALL require(COUNT(.NOT. (ABS(p2 - p) > 0.0_wp)) == 0, 'S1: another node draws other phases')
    CALL CD_VIV_Sto_Node(SIG, C, CF, NF, 8, 12, a2, f2, p2)
    CALL require(COUNT(.NOT. (ABS(p2 - p) > 0.0_wp)) == 0, 'S1: another seed draws other phases')
  END SUBROUTINE check_process

  SUBROUTINE write_current(path)
    CHARACTER(*), INTENT(IN) :: path
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'MoorDyn v2 water kinematics file'
    WRITE (u, '(A)') 'uniform current along y for the VIV STOCHASTIC deck gate'
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

  SUBROUTINE write_deck(path, ei, tmax_, sto_row, cl, wake)
    !! The taut cable of EI ei in the current; sto_row: the VIV STOCHASTIC data row (no section
    !! without it); cl the LINE TYPES Cl (default 0.8); wake: a WAKE OSCILLATOR row instead of the lift.
    CHARACTER(*), INTENT(IN) :: path, ei
    REAL(wp), INTENT(IN) :: tmax_
    CHARACTER(*), INTENT(IN), OPTIONAL :: sto_row, cl
    LOGICAL, INTENT(IN), OPTIONAL :: wake
    CHARACTER(8) :: cl_col
    LOGICAL :: wo
    INTEGER :: u
    cl_col = '0.8'
    IF (PRESENT(cl)) cl_col = cl
    wo = .FALSE.
    IF (PRESENT(wake)) wo = wake
    IF (wo) cl_col = '0.0'
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (u, '(A)') 'VIV STOCHASTIC deck gate: taut cable in a uniform current'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (u, '(A)') 'TypeName   Diam   Mass/m   EA        BA/-zeta  EI   Cd    Ca    CdAx  CaAx  Cl   dF    cF'
    WRITE (u, '(A)') 'riser      0.10   10.0     5.0e7     -0.8      '//ei//'  1.2   1.0   0.0   0.0   '// &
      TRIM(cl_col)//'  0.08  0.18'
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
    IF (wo) THEN
      WRITE (u, '(A)') '--------------------- WAKE OSCILLATOR ----------------------------------'
      WRITE (u, '(A)') 'TypeName CL0 St eps A'
      WRITE (u, '(A)') 'riser 0.3 0.2 0.3 12'
    END IF
    IF (PRESENT(sto_row)) THEN
      WRITE (u, '(A)') '--------------------- VIV STOCHASTIC -----------------------------------'
      WRITE (u, '(A)') 'TypeName   sigma   c     Nf   seed'
      WRITE (u, '(A)') sto_row
    END IF
    WRITE (u, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (u, '(A)') '0.01     dtM'
    WRITE (u, '(A)') '100.0    WtrDpth'
    WRITE (u, '(A)') '1025.0   rhoW'
    WRITE (u, '(A)') '9.80665  g'
    WRITE (u, '(A)') 'vsd_current.dat  WaterKin'
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

  REAL(wp) FUNCTION amplitude(root) RESULT(ad)
    !! Half the peak-to-peak mid-span cross-flow excursion over D over T0-TMAX (0 without rows).
    CHARACTER(*), INTENT(IN) :: root
    CHARACTER(1024) :: rec
    REAL(wp) :: row(2), zmin, zmax
    INTEGER :: u, ios
    ad = 0.0_wp
    zmin = HUGE(1.0_wp)
    zmax = -HUGE(1.0_wp)
    OPEN (NEWUNIT=u, FILE=root//'.out', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    READ (u, '(A)') rec
    READ (u, '(A)') rec
    DO
      READ (u, '(A)', IOSTAT=ios) rec
      IF (ios /= 0) EXIT
      IF (LEN_TRIM(rec) == 0) CYCLE
      READ (rec, *) row
      IF (row(1) < T0) CYCLE
      zmin = MIN(zmin, row(2))
      zmax = MAX(zmax, row(2))
    END DO
    CLOSE (u)
    IF (zmax > zmin) ad = 0.5_wp*(zmax - zmin)/D
  END FUNCTION amplitude

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

  SUBROUTINE check_deck(ei, tag)
    !! S2 and S3 on the cable of EI ei.
    CHARACTER(*), INTENT(IN) :: ei, tag
    LOGICAL :: ok, same
    CHARACTER(1024) :: msg
    CHARACTER(:), ALLOCATABLE :: lab
    REAL(wp) :: ad0, ad1
    lab = ' (EI '//ei//')'
    CALL write_deck(tag//'_none.dat', ei, TMAX)
    CALL run(tag//'_none.dat', tag//'_none', ok, msg)
    CALL require(ok, 'S2: the deck without the section runs'//lab//': '//TRIM(msg))
    CALL write_deck(tag//'_s0.dat', ei, TMAX, sto_row='riser 0 0.07 100 3')
    CALL run(tag//'_s0.dat', tag//'_s0', ok, msg)
    CALL require(ok, 'S2: the deck with sigma 0 runs'//lab//': '//TRIM(msg))
    CALL require(same_file(tag//'_none.out', tag//'_s0.out'), 'S2: a sigma 0 row is the deck without it, '// &
                 'bit for bit'//lab)
    CALL write_deck(tag//'_s1.dat', ei, TMAX, sto_row='riser 0.02 0.07 100 3')
    CALL run(tag//'_s1.dat', tag//'_s1', ok, msg)
    CALL require(ok, 'S3: the deck with sigma 0.02 runs'//lab//': '//TRIM(msg))
    CALL run(tag//'_s1.dat', tag//'_s1b', ok, msg)
    same = same_file(tag//'_s1.out', tag//'_s1b.out')
    CALL require(ok .AND. same, 'S3: the run repeats bit for bit'//lab)
    CALL require(.NOT. same_file(tag//'_none.out', tag//'_s1.out'), 'S3: the centre moves the response'//lab)
    CALL write_deck(tag//'_s2.dat', ei, TMAX, sto_row='riser 0.02 0.07 100 4')
    CALL run(tag//'_s2.dat', tag//'_s2', ok, msg)
    same = same_file(tag//'_s1.out', tag//'_s2.out')
    CALL require(ok .AND. .NOT. same, 'S3: another seed, another run'//lab)
    ad0 = amplitude(tag//'_none')
    ad1 = amplitude(tag//'_s1')
    WRITE (*, '(A,2F8.4)') 'S3'//lab//': A/D over 50-80 s with the constant and the stochastic centre ', ad0, ad1
    CALL require(ad0 > 0.2_wp .AND. ABS(ad1/ad0 - 1.0_wp) <= 0.15_wp, 'S3: the cable stays locked in'//lab)
  END SUBROUTINE check_deck

  SUBROUTINE refuse(tag, sto_row, want, label, cl, wake)
    CHARACTER(*), INTENT(IN) :: tag, sto_row, want, label
    CHARACTER(*), INTENT(IN), OPTIONAL :: cl
    LOGICAL, INTENT(IN), OPTIONAL :: wake
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck(tag//'.dat', '1.0', 1.0_wp, sto_row=sto_row, cl=cl, wake=wake)
    CALL run(tag//'.dat', tag, ok, msg)
    WRITE (*, '(A)') 'S4 '//label//': '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, want) > 0, 'S4: '//label//' is refused')
  END SUBROUTINE refuse

  SUBROUTINE check_refusals()
    CALL refuse('vsd_nocl', 'riser 0.02 0.07', 'give the type LINE TYPES Cl > 0 and no WAKE OSCILLATOR row', &
                'a type without Cl', cl='0.0')
    CALL refuse('vsd_wake', 'riser 0.02 0.07', 'give the type LINE TYPES Cl > 0 and no WAKE OSCILLATOR row', &
                'a wake-oscillator type', wake=.TRUE.)
    CALL refuse('vsd_c', 'riser 0.02 1.5', '0 < c <= 1', 'c above 1')
    CALL refuse('vsd_nf', 'riser 0.02 0.07 10.5', 'Nf must be an integer', 'a non-integer Nf')
    CALL refuse('vsd_seed', 'riser 0.02 0.07 100 -1', 'seed must be an integer', 'a negative seed')
  END SUBROUTINE check_refusals

END PROGRAM test_viv_stochastic
