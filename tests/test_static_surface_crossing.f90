! File: tests/test_static_surface_crossing.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_static_surface_crossing
  !! Statics of a light finite-EI cable through the free surface in a cross current: a 4.95 m
  !! laboratory lazy-wave cable (OD 31 mm, 0.919 kg/m, EI 0.217 N m^2, fresh water 2 m deep), its top
  !! 0.753 m above the surface and clamped vertical, its foot clamped on the floor, five buoyancy
  !! modules as attachments, in a uniform current across its plane at 0.390-0.410 m/s every 1 mm/s.
  !! At 0.400 m/s a Gauss point of the element through the surface sits at the waterline in the
  !! equilibrium; with a wet/dry switch per Gauss point the drag jumped there, the residual had no
  !! root and every static route stalled (current load stepping stopped at a drag fraction 0.997).
  !! Gates:
  !!   S1 every speed solves;
  !!   S2 the sideways deflection and the top tension grow monotonically with the speed, and the
  !!      deflection at 0.400 m/s lies within 2 % of a step of the mean of its neighbours.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  INTEGER, PARAMETER :: NS = 21
  REAL(wp) :: speed(NS), ymax(NS), ttop(NS)
  INTEGER :: nfail = 0, k
  LOGICAL :: ok, all_ok
  CHARACTER(1024) :: msg
  CHARACTER(64) :: root

  all_ok = .TRUE.
  DO k = 1, NS
    speed(k) = 0.390_wp + 0.001_wp*REAL(k - 1, wp)
    WRITE (root, '(A,I2.2)') 'ssc_', k
    CALL write_deck(TRIM(root)//'.dat', speed(k))
    CALL run(TRIM(root)//'.dat', TRIM(root), ok, msg)
    CALL require(ok, 'S1: the static solve at the speed of case '//TRIM(root)//' converges: '//TRIM(msg))
    all_ok = all_ok .AND. ok
    IF (ok) CALL read_static(TRIM(root)//'.static.out', ymax(k), ttop(k))
  END DO
  IF (all_ok) THEN
    WRITE (*, '(A,3F10.5)') 'S2: peak sideways deflection [m] at 0.399/0.400/0.401 m/s: ', ymax(10:12)
    CALL require(ALL(ymax(2:NS) > ymax(1:NS - 1)) .AND. ALL(ttop(2:NS) > ttop(1:NS - 1)), &
                 'S2: the deflection and the top tension grow monotonically with the speed')
    CALL require(ABS(ymax(11) - 0.5_wp*(ymax(10) + ymax(12))) <= 0.02_wp*ABS(ymax(12) - ymax(10)), &
                 'S2: the deflection at 0.400 m/s lies between its neighbours (no jump)')
  END IF
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: statics of a cable through the free surface in a cross current'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE write_deck(path, u)
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), INTENT(IN) :: u
    INTEGER :: io
    ! a uniform current across the cable's plane (MoorDyn water-kinematics file, CurrentMod 1)
    OPEN (NEWUNIT=io, FILE='ssc_current.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE (io, '(A)') 'MoorDyn v2 water kinematics file'
    WRITE (io, '(A)') 'uniform current along y for the surface-crossing statics gate'
    WRITE (io, '(A)') '--------------------------- WAVES -------------------------------------'
    WRITE (io, '(A)') '0                    WaveKinMod  - type of wave input'
    WRITE (io, '(A)') '""                   WaveKinFile - file containing wave elevation time series'
    WRITE (io, '(A)') '0.000000E+00  dtWave      - time step to use in setting up wave kinematics grid (s)'
    WRITE (io, '(A)') '0                    WaveDir     - wave heading (deg)'
    WRITE (io, '(A)') '2                                - X wave input type'
    WRITE (io, '(A)') '-5, 5, 2                         - X wave grid point data'
    WRITE (io, '(A)') '2                                - Y wave input type'
    WRITE (io, '(A)') '-5, 5, 2                         - Y wave grid point data'
    WRITE (io, '(A)') '2                                - Z wave input type'
    WRITE (io, '(A)') '-2, 0, 2                         - Z wave grid point data'
    WRITE (io, '(A)') '--------------------------- CURRENT -------------------------------------'
    WRITE (io, '(A)') '1                    CurrentMod  - type of current input'
    WRITE (io, '(A)') 'z-depth     x-current      y-current'
    WRITE (io, '(A)') '(m)           (m/s)         (m/s)'
    WRITE (io, '(A,F8.5)') '0.0    0.0  ', u
    WRITE (io, '(A,F8.5)') '-2.0   0.0  ', u
    WRITE (io, '(A)') '--------------------- need this line ------------------'
    CLOSE (io)
    OPEN (NEWUNIT=io, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (io, '(A)') '--------------------- CableDyn Input File ---------------------------------'
    WRITE (io, '(A)') 'light lazy-wave cable through the free surface in a cross current (statics)'
    WRITE (io, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (io, '(A)') 'TypeName Diam Mass EA BA EI Cd_n Cd_t Ca_n Ca_t Cl dF cF'
    WRITE (io, '(A)') '(-) (m) (kg/m) (N) (N-s) (N-m^2) (-) (-) (-) (-) (-) (-) (-)'
    WRITE (io, '(A)') 'cab 0.031 0.919 1.0e4 0.0 0.217 1.2 0.0 1.0 0.0 0 0.08 0.18'
    WRITE (io, '(A)') '----------------------- POINTS ----------------------------------------------'
    WRITE (io, '(A)') 'ID Attachment X Y Z M V CdA CA'
    WRITE (io, '(A)') '(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)'
    WRITE (io, '(A)') '1 Coupled 0.0 0.0 0.753 0 0 0 0'
    WRITE (io, '(A)') '2 Fixed 2.798 0.0 -1.954 0 0 0 0'
    WRITE (io, '(A)') '----------------------- LINES -----------------------------------------------'
    WRITE (io, '(A)') 'ID LineType AttachA AttachB UnstrLen NumSegs Outputs'
    WRITE (io, '(A)') '(-) (-) (-) (-) (m) (-) (-)'
    WRITE (io, '(A)') '1 cab 1 2 4.95 100 -'
    WRITE (io, '(A)') '--------------------- ATTACHMENTS ----------------------------------'
    WRITE (io, '(A)') 'LineID  ArcLength  Mass  Volume  CdA  Ca'
    WRITE (io, '(A)') '(-)  (m)  (kg)  (m^3)  (m^2)  (-)'
    WRITE (io, '(A)') '1 3.1432 0.012375 1.237473e-04 0.0034 1.0'
    WRITE (io, '(A)') '1 3.2299 0.012375 1.237473e-04 0.0034 1.0'
    WRITE (io, '(A)') '1 3.4845 0.012375 1.237473e-04 0.0034 1.0'
    WRITE (io, '(A)') '1 3.7516 0.012375 1.237473e-04 0.0034 1.0'
    WRITE (io, '(A)') '1 4.0420 0.012375 1.237473e-04 0.0034 1.0'
    WRITE (io, '(A)') '--------------------- END CONNECTIONS -------------------------------'
    WRITE (io, '(A)') 'LineID  End  Stiffness  EzX  EzY  EzZ'
    WRITE (io, '(A)') '(-)     (-)  (N-m/rad)  (-)  (-)  (-)'
    WRITE (io, '(A)') '1       B    Rigid      0.644    0    -0.765'
    WRITE (io, '(A)') '1       A    Rigid      0    0    -1'
    WRITE (io, '(A)') '---------------------- OPTIONS ----------------------------------------------'
    WRITE (io, '(A)') '1000.0 rhoW'
    WRITE (io, '(A)') '9.80665 g'
    WRITE (io, '(A)') '2.0 WtrDpth'
    WRITE (io, '(A)') 'ssc_current.dat WaterKin'
    WRITE (io, '(A)') '0.001 dtM'
    WRITE (io, '(A)') '0.002 TMax'
    WRITE (io, '(A)') '---------------------- OUTPUTS -----------------------------------------------'
    WRITE (io, '(A)') 'FairTen1'
    WRITE (io, '(A)') 'END'
    WRITE (io, '(A)') '------------------------- need this line --------------------------------------'
    CLOSE (io)
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

  SUBROUTINE read_static(path, ym, t1)
    !! The largest sideways (y) node position and the tension at the top node of a static table.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), INTENT(OUT) :: ym, t1
    CHARACTER(512) :: rec
    INTEGER :: io, ios, line_id, node
    REAL(wp) :: v(10)
    ym = -HUGE(1.0_wp)
    t1 = 0.0_wp
    OPEN (NEWUNIT=io, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    READ (io, '(A)') rec
    READ (io, '(A)') rec
    READ (io, '(A)') rec
    DO
      READ (io, '(A)', IOSTAT=ios) rec
      IF (ios /= 0) EXIT
      IF (LEN_TRIM(rec) == 0) CYCLE
      READ (rec, *, IOSTAT=ios) line_id, node, v
      IF (ios /= 0) CYCLE
      ym = MAX(ym, v(3))
      IF (node == 1) t1 = v(5)
    END DO
    CLOSE (io)
  END SUBROUTINE read_static

END PROGRAM test_static_surface_crossing
