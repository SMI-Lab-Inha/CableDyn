! File: tests/test_wavekin_modes.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
!> MoorDyn-C water-kinematics OPTIONS (WaveKin 3 and 7, Currents 1) against MoorDyn-C v2.7.1 on
!> shared inputs (tests/data/wavekin_mdc, reference rows from validation/scripts/
!> wavekin_mdc_probe.cpp): the deck water kinematics CableDyn samples (CD_Deck_Fluid_Probe) at the
!> points and times MoorDyn-C reported.
!>  - WaveKin 7 (component sum at the node, both codes): elevation, velocity and acceleration
!>    identical to 1e-6 of their scale (MoorDyn-C's wave number is Newman's seven-digit
!>    approximation, CableDyn's the dispersion root);
!>  - WaveKin 3 (elevation record reduced to Fourier components): at the grid points and wave
!>    time samples, where MoorDyn-C's grid interpolation is exact, the elevation within 1e-6 m
!>    (MoorDyn-C's own reconstruction departs from the exact record by up to 3e-7 m; CableDyn's
!>    equals the record's two harmonics to 1e-10 m), and velocity and acceleration, which
!>    MoorDyn-C interpolates in the stretched depth on its 0.25 m grid, within 1e-3 of their
!>    scale;
!>  - Currents 1: identical to 1e-12 m/s;
!>  - the MoorDyn-C modes CableDyn does not support fail closed by name.
PROGRAM test_wavekin_modes
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Deck_Fluid_Probe, CD_DECKDRV_OK
  IMPLICIT NONE
  INTEGER :: n_fail = 0

  CALL check_mode('wavekin_mdc/wavekin7.dat', 'wavekin_mdc/mdc_wavekin7.csv', 'WaveKin 7', 1.0e-6_wp, &
                  1.0e-6_wp, 1.0e-6_wp)
  CALL check_mode('wavekin_mdc/wavekin3.dat', 'wavekin_mdc/mdc_wavekin3.csv', 'WaveKin 3', 1.0e-6_wp, &
                  1.0e-3_wp, 1.0e-3_wp)
  CALL check_record()
  CALL check_mode('wavekin_mdc/currents1.dat', 'wavekin_mdc/mdc_currents1.csv', 'Currents 1', 1.0e-12_wp, &
                  1.0e-12_wp, 1.0e-12_wp, absolute=.TRUE.)
  CALL check_rejections()
  CALL check_still_water_ramp()
  IF (n_fail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', n_fail, ' water-kinematics gate(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: MoorDyn-C WaveKin 3/7 and Currents 1 kinematics match MoorDyn-C'

CONTAINS

  INCLUDE 'nan_max_abs.inc'

  SUBROUTINE require(cond, msg)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: msg
    IF (.NOT. cond) THEN
      n_fail = n_fail + 1
      WRITE (*, '(A,A)') 'FAILED: ', msg
    END IF
  END SUBROUTINE require

  SUBROUTINE read_reference(path, ref, ok)
    !! The rows "t,x,y,z,zeta,ux,uy,uz,ax,ay,az" of a probe CSV into ref(11, n).
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: ref(:, :)
    LOGICAL, INTENT(OUT) :: ok
    INTEGER :: u, ios, n, i
    CHARACTER(1024) :: buf
    ok = .FALSE.
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    n = 0
    DO
      READ (u, '(A)', IOSTAT=ios) buf
      IF (ios /= 0) EXIT
      IF (LEN_TRIM(buf) > 0) n = n + 1
    END DO
    REWIND (u)
    READ (u, '(A)') buf
    ALLOCATE (ref(11, n - 1))
    DO i = 1, n - 1
      READ (u, *, IOSTAT=ios) ref(:, i)
      IF (ios /= 0) THEN
        CLOSE (u)
        RETURN
      END IF
    END DO
    CLOSE (u)
    ok = n > 1
  END SUBROUTINE read_reference

  SUBROUTINE check_mode(deck, csv, label, tol_eta, tol_vel, tol_acc, absolute)
    CHARACTER(*), INTENT(IN) :: deck, csv, label
    REAL(wp), INTENT(IN) :: tol_eta, tol_vel, tol_acc
    LOGICAL, INTENT(IN), OPTIONAL :: absolute
    REAL(wp), ALLOCATABLE :: ref(:, :), eta(:), vel(:, :), acc(:, :)
    REAL(wp) :: e_eta, e_vel, e_acc, s_eta, s_vel, s_acc
    INTEGER :: es, n
    LOGICAL :: ok, abs_err
    CHARACTER(512) :: em
    abs_err = .FALSE.
    IF (PRESENT(absolute)) abs_err = absolute
    CALL read_reference(csv, ref, ok)
    CALL require(ok, label//': MoorDyn-C reference '//csv//' readable')
    IF (.NOT. ok) RETURN
    n = SIZE(ref, 2)
    ALLOCATE (eta(n), vel(3, n), acc(3, n))
    CALL CD_Deck_Fluid_Probe(deck, ref(2:4, :), ref(1, :), eta, vel, acc, es, em)
    CALL require(es == CD_DECKDRV_OK, label//': deck water kinematics evaluated: '//TRIM(em))
    IF (es /= CD_DECKDRV_OK) RETURN
    e_eta = nan_max_abs(eta - ref(5, :))
    e_vel = nan_max_abs(vel - ref(6:8, :))
    e_acc = nan_max_abs(acc - ref(9:11, :))
    s_eta = 1.0_wp
    s_vel = 1.0_wp
    s_acc = 1.0_wp
    IF (.NOT. abs_err) THEN
      s_vel = nan_max_abs(ref(6:8, :))
      s_acc = nan_max_abs(ref(9:11, :))
    END IF
    WRITE (*, '(A,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A)') label, ' vs MoorDyn-C (', n, ' samples): elevation ', &
      e_eta, ' m, velocity ', e_vel/s_vel, ', acceleration ', e_acc/s_acc, MERGE(' (absolute)', ' (relative)', &
                                                                                 abs_err)
    CALL require(e_eta <= tol_eta, label//': elevation matches MoorDyn-C')
    CALL require(e_vel <= tol_vel*s_vel, label//': velocity matches MoorDyn-C')
    CALL require(e_acc <= tol_acc*s_acc, label//': acceleration matches MoorDyn-C')
  END SUBROUTINE check_mode

  SUBROUTINE check_record()
    !! WaveKin 3 reproduces the recorded elevation 0.8 cos(2 pi t/8) + 0.4 sin(2 pi t/5 + 0.3) of
    !! wave_elevation.txt (a 40 s periodic record) at the origin, between the samples too.
    REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
    REAL(wp) :: pos(3, 7), t(7), eta(7), vel(3, 7), acc(3, 7), err
    INTEGER :: es, i
    CHARACTER(512) :: em
    pos = 0.0_wp
    pos(3, :) = -10.0_wp
    t = [0.0_wp, 0.37_wp, 1.23_wp, 5.55_wp, 13.01_wp, 27.7_wp, 39.95_wp]
    CALL CD_Deck_Fluid_Probe('wavekin_mdc/wavekin3.dat', pos, t, eta, vel, acc, es, em)
    CALL require(es == CD_DECKDRV_OK, 'WaveKin 3 record check evaluated: '//TRIM(em))
    IF (es /= CD_DECKDRV_OK) RETURN
    err = 0.0_wp
    DO i = 1, SIZE(t)
      err = MAX(err, ABS(eta(i) - (0.8_wp*COS(2.0_wp*PI*t(i)/8.0_wp) + 0.4_wp*SIN(2.0_wp*PI*t(i)/5.0_wp + 0.3_wp))))
    END DO
    WRITE (*, '(A,ES10.3,A)') 'WaveKin 3 vs the recorded harmonics: elevation ', err, ' m'
    CALL require(err <= 1.0e-10_wp, 'WaveKin 3 reproduces the recorded elevation')
  END SUBROUTINE check_record

  SUBROUTINE write_fail_deck(path, wavekin, currents, extra_option)
    !! The shared probe deck with the given WaveKin and Currents rows (and an optional extra
    !! OPTIONS row).
    CHARACTER(*), INTENT(IN) :: path, wavekin, currents
    CHARACTER(*), INTENT(IN), OPTIONAL :: extra_option
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- MoorDyn Input File ------------------------------------'
    WRITE (u, '(A)') 'unsupported water-kinematics mode'
    WRITE (u, '(A)') '----------------------- LINE TYPES ------------------------------------------'
    WRITE (u, '(A)') 'TypeName   Diam    Mass/m     EA         BA/-zeta    EI         Cd     Ca     CdAx    CaAx'
    WRITE (u, '(A)') '(name)     (m)     (kg/m)     (N)        (N-s/-)     (N-m^2)    (-)    (-)    (-)     (-)'
    WRITE (u, '(A)') 'chain      0.1     50.0       1.0e7      0.0         0          1.2    1.0    0.2     0.0'
    WRITE (u, '(A)') '---------------------- POINT PROPERTIES --------------------------------'
    WRITE (u, '(A)') 'ID    Type      X       Y       Z       Mass   Volume  CdA    Ca'
    WRITE (u, '(A)') '(#)   (-)       (m)     (m)     (m)     (kg)   (m^3)   (m^2)  (-)'
    WRITE (u, '(A)') '1     Fixed     -30.0   0.0     -20.0   0      0       0      0'
    WRITE (u, '(A)') '2     Vessel    30.0    0.0     -20.0   0      0       0      0'
    WRITE (u, '(A)') '---------------------- LINES ----------------------------------------'
    WRITE (u, '(A)') 'ID   LineType   AttachA  AttachB  UnstrLen  NumSegs  LineOutputs'
    WRITE (u, '(A)') '(#)   (name)     (#)      (#)       (m)       (-)     (-)'
    WRITE (u, '(A)') '1     chain      2        1         70.0      20      -'
    WRITE (u, '(A)') '---------------------- OPTIONS -----------------------------------------'
    WRITE (u, '(A)') '0.002         dtM'
    WRITE (u, '(A)') '10.0          TMax'
    WRITE (u, '(A)') '50            WtrDpth'
    WRITE (u, '(A)') wavekin//'             WaveKin'
    WRITE (u, '(A)') currents//'             Currents'
    IF (PRESENT(extra_option)) WRITE (u, '(A)') extra_option
    WRITE (u, '(A)') '------------------------- need this line --------------------------------------'
    CLOSE (u)
  END SUBROUTINE write_fail_deck

  SUBROUTINE expect_failure(wavekin, currents, needle)
    CHARACTER(*), INTENT(IN) :: wavekin, currents, needle
    REAL(wp) :: pos(3, 1), t(1), eta(1), vel(3, 1), acc(3, 1)
    INTEGER :: es
    CHARACTER(512) :: em
    pos = 0.0_wp
    pos(3, 1) = -5.0_wp
    t = 0.0_wp
    CALL write_fail_deck('wavekin_fail.dat', wavekin, currents)
    CALL CD_Deck_Fluid_Probe('wavekin_fail.dat', pos, t, eta, vel, acc, es, em)
    CALL require(es /= CD_DECKDRV_OK .AND. INDEX(em, needle) > 0, 'WaveKin '//wavekin//' Currents '//currents// &
                 ' fails closed naming "'//needle//'": '//TRIM(em))
  END SUBROUTINE expect_failure

  SUBROUTINE check_still_water_ramp()
    !! rampTime is valid without waves (it also ramps prescribed motion). The still-water
    !! kinematics must then be exactly zero, including at t = 0 where the ramp is closed
    !! (r = 0, dr/dt = 0), with no 0/0 in the ramped-acceleration term.
    REAL(wp) :: pos(3, 3), t(3), eta(3), vel(3, 3), acc(3, 3)
    INTEGER :: es
    CHARACTER(512) :: em
    pos = 0.0_wp
    pos(3, :) = -5.0_wp
    t = [0.0_wp, 1.0_wp, 5.0_wp]
    CALL write_fail_deck('wavekin_ramp_still.dat', '0', '0', '4.0          rampTime')
    CALL CD_Deck_Fluid_Probe('wavekin_ramp_still.dat', pos, t, eta, vel, acc, es, em)
    CALL require(es == CD_DECKDRV_OK, 'still-water deck with rampTime probes: '//TRIM(em))
    CALL require(ALL(ABS(eta) <= 0.0_wp) .AND. ALL(ABS(vel) <= 0.0_wp) .AND. ALL(ABS(acc) <= 0.0_wp), &
                 'still water with rampTime: zero (and finite) kinematics at t = 0, inside and after the ramp')
  END SUBROUTINE check_still_water_ramp

  SUBROUTINE check_rejections()
    INTEGER :: u
    CALL expect_failure('1', '0', 'MoorDyn-C API')
    CALL expect_failure('2', '0', 'FFT_GRID')
    CALL expect_failure('4', '0', 'no kinematics source in MoorDyn-C')
    CALL expect_failure('6', '0', 'no kinematics source in MoorDyn-C')
    CALL expect_failure('0', '2', 'current_profile_dynamic.txt')
    CALL expect_failure('0', '5', 'current_profile_4d.txt')
    ! short-crested components (rows with different directions)
    OPEN (NEWUNIT=u, FILE='wave_frequencies.txt', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '0.0 0.0 0.0 0.0'
    WRITE (u, '(A)') '0.6 0.7 0.2 0.0'
    WRITE (u, '(A)') '0.9 -0.3 0.5 1.57'
    CLOSE (u)
    CALL expect_failure('7', '0', 'different directions')
  END SUBROUTINE check_rejections
END PROGRAM test_wavekin_modes
