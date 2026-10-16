! File: tests/test_point_seabed_contact.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_point_seabed_contact
  !! Seabed contact of Free/Connect points, and the fold note of EI = 0 lines.
  !!
  !! A 50 t clump weight (Mass only, Vol = 0) joins two polyester lines on the seabed. Its
  !! footprint is A = (Mass/7850)^(2/3) (the volume of solid steel of that mass), so its own
  !! contact has the stiffness kBot*A. Without it only the two line end nodes held the clump
  !! and it sank 2.8 m below the seabed. Checked on every route that steps the point:
  !!   * the static equilibrium (TMax = 0, multibody march): the clump rests within the
  !!     penetration of its own footprint, 0 < pen <= W/(kBot*A), and its contact reaction
  !!     kBot*A*pen balances its weight with the line forces on it (Point<P>Fz);
  !!   * the dynamic march, monolithic and staggered (the point system step), from that
  !!     equilibrium: the clump stays at rest;
  !!   * the OpenFAST aggregate (the coupled point step): the clump stays at rest;
  !!   * a sloped structured seabed: the clump rests on the bathymetry, also at rest in the
  !!     march.
  !! An EI = 0 line pushed past its anchor by a current on a frictionless seabed (a hairpin)
  !! is kept and named by a note; through the OpenFAST aggregate the note is returned with
  !! the line's initialisation report (CD_AGG_GetInitLine).
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  USE CableDyn_OpenFAST_Aggregate, ONLY: CD_AGG_ModuleType, CD_AGG_Init_From_Deck, CD_AGG_GetInitLine, &
                                         CD_AGG_EvalChannel, CD_AGG_GetMovingPointMesh, CD_AGG_Step_Moving, &
                                         CD_AGG_CalcOutput, CD_AGG_End, CD_AGG_OK
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE

  REAL(wp), PARAMETER :: G = 9.80665_wp, KBOT = 1.0e5_wp, MASS = 50000.0_wp, STEEL = 7850.0_wp
  ! columns of the .out rows: Time, FairTen1, AnchTen1, FairTen2, AnchTen2, Point3Fx/Fy/Fz, Point3px/py/pz
  INTEGER, PARAMETER :: NCOL = 11
  INTEGER :: failures
  REAL(wp) :: area, k_clump, w_clump

  failures = 0
  area = (MASS/STEEL)**(2.0_wp/3.0_wp)
  k_clump = KBOT*area
  w_clump = MASS*G

  CALL case_flat_static()
  CALL case_flat_march('monolithic')
  CALL case_flat_march('staggered')
  CALL case_slope()
  CALL case_flat_aggregate()
  CALL case_fold_note_aggregate()

  IF (failures /= 0) THEN
    WRITE (*, '(A,I0,A)') 'test_point_seabed_contact: ', failures, ' failure(s)'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'test_point_seabed_contact: all checks passed'

CONTAINS

  SUBROUTINE case_flat_static()
    REAL(wp), ALLOCATABLE :: rows(:, :)
    REAL(wp) :: pen, residual
    INTEGER :: es
    LOGICAL :: conv
    CHARACTER(512) :: em

    CALL write_clump_deck('point_seabed_flat_static.dat', '0.0', '', .FALSE.)
    CALL CD_Run_Deck_Driver('point_seabed_flat_static.dat', 'point_seabed_flat_static', conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, 'flat static: run: '//TRIM(em))
    IF (es /= CD_DECKDRV_OK) RETURN
    CALL read_rows('point_seabed_flat_static.out', rows)
    IF (SIZE(rows, 2) < 1) THEN
      CALL require(.FALSE., 'flat static: no output row')
      RETURN
    END IF
    pen = -200.0_wp - rows(11, 1)
    WRITE (*, '(A,F8.4,A,F8.4,A)') '  flat seabed: clump ', pen, ' m below the seabed (footprint limit ', &
      w_clump/k_clump, ' m)'
    ! base behaviour: 2.8 m, held by the two line end nodes only
    CALL require(pen > 0.0_wp, 'flat static: the clump rests on the seabed')
    CALL require(pen <= w_clump/k_clump, 'flat static: the clump rests within the penetration of its footprint')
    ! weight = line forces on the point + the footprint's reaction (Point3Fz is printed to 8 digits)
    residual = rows(8, 1) + k_clump*pen - w_clump
    CALL require(ABS(residual) <= 1.0e-4_wp*w_clump, 'flat static: the footprint reaction kBot*A*pen balances '// &
                 'the clump with its lines')
  END SUBROUTINE case_flat_static

  SUBROUTINE case_flat_march(scheme)
    CHARACTER(*), INTENT(IN) :: scheme
    REAL(wp), ALLOCATABLE :: rows(:, :)
    INTEGER :: es
    LOGICAL :: conv
    CHARACTER(512) :: em
    CHARACTER(64) :: root

    root = 'point_seabed_flat_'//scheme
    CALL write_clump_deck(TRIM(root)//'.dat', '2.0', scheme, .FALSE.)
    CALL CD_Run_Deck_Driver(TRIM(root)//'.dat', TRIM(root), conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, 'flat '//scheme//' march: run: '//TRIM(em))
    IF (es /= CD_DECKDRV_OK) RETURN
    CALL read_rows(TRIM(root)//'.out', rows)
    CALL require(SIZE(rows, 2) > 10, 'flat '//scheme//' march: output rows')
    IF (SIZE(rows, 2) < 2) RETURN
    CALL require(ALL(IEEE_IS_FINITE(rows)), 'flat '//scheme//' march: finite outputs')
    CALL require(-200.0_wp - rows(11, 1) > 0.0_wp .AND. -200.0_wp - rows(11, 1) <= w_clump/k_clump, &
                 'flat '//scheme//' march: the clump starts on its footprint')
    CALL require(MAXVAL(ABS(rows(9:11, :) - SPREAD(rows(9:11, 1), 2, SIZE(rows, 2)))) < 1.0e-3_wp, &
                 'flat '//scheme//' march: the clump stays at rest on the seabed')
  END SUBROUTINE case_flat_march

  SUBROUTINE case_slope()
    REAL(wp), ALLOCATABLE :: rows(:, :)
    REAL(wp) :: floor_z, gap, s
    INTEGER :: es, u
    LOGICAL :: conv
    CHARACTER(512) :: em

    ! a plane rising 2 % toward +x under the whole line, depth 200 m at x = 450 m
    OPEN (NEWUNIT=u, FILE='point_seabed_slope_bathy.txt', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '-100.0 -100.0 211.0'
    WRITE (u, '(A)') '-100.0 100.0 211.0'
    WRITE (u, '(A)') '900.0 -100.0 191.0'
    WRITE (u, '(A)') '900.0 100.0 191.0'
    CLOSE (u)
    CALL write_clump_deck('point_seabed_slope.dat', '1.0', 'staggered', .TRUE.)
    CALL CD_Run_Deck_Driver('point_seabed_slope.dat', 'point_seabed_slope', conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, 'slope: run: '//TRIM(em))
    IF (es /= CD_DECKDRV_OK) RETURN
    CALL read_rows('point_seabed_slope.out', rows)
    CALL require(SIZE(rows, 2) > 10, 'slope: output rows')
    IF (SIZE(rows, 2) < 2) RETURN
    floor_z = -(211.0_wp - 0.02_wp*(rows(9, 1) + 100.0_wp))
    gap = floor_z - rows(11, 1)
    s = SQRT(1.0_wp + 0.02_wp**2)
    WRITE (*, '(A,F8.4,A)') '  sloped seabed: clump ', gap, ' m below the seabed'
    ! the normal law pushes along the surface normal with k*gap/s, whose vertical part k*gap/s^2
    ! (with the lines' share) carries the weight
    CALL require(gap > 0.0_wp .AND. gap <= s*s*w_clump/k_clump, 'slope: the clump rests on its footprint')
    CALL require(MAXVAL(ABS(rows(9:11, :) - SPREAD(rows(9:11, 1), 2, SIZE(rows, 2)))) < 1.0e-3_wp, &
                 'slope: the clump stays at rest on the sloped seabed')
  END SUBROUTINE case_slope

  SUBROUTINE case_flat_aggregate()
    !! The OpenFAST aggregate (the coupled point step): the clump starts on its footprint and
    !! stays at rest with the fairlead held. Channel 10 is Point3pz.
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp), PARAMETER :: DT = 0.05_wp
    REAL(wp) :: r0(3, 1), v0(3, 1), a0(3, 1), l0(3, 1), z0, z, zmax_dev
    INTEGER :: es, s, n_iter
    LOGICAL :: cv, st
    CHARACTER(512) :: em

    CALL write_clump_deck('point_seabed_agg.dat', '0.0', '', .FALSE.)
    CALL CD_AGG_Init_From_Deck(agg, 'point_seabed_agg.dat', DT, es, em)
    CALL require(es == CD_AGG_OK, 'aggregate clump: init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_EvalChannel(agg, 10, z0, es, em)
    CALL require(es == CD_AGG_OK, 'aggregate clump: Point3pz: '//TRIM(em))
    CALL require(-200.0_wp - z0 > 0.0_wp .AND. -200.0_wp - z0 <= w_clump/k_clump, &
                 'aggregate clump: the clump starts on its footprint')
    CALL CD_AGG_GetMovingPointMesh(agg, r0, v0, a0, l0, es, em)
    CALL require(es == CD_AGG_OK, 'aggregate clump: moving mesh: '//TRIM(em))
    v0 = 0.0_wp
    a0 = 0.0_wp
    zmax_dev = 0.0_wp
    DO s = 1, 40
      CALL CD_AGG_Step_Moving(agg, DT, r0, v0, a0, cv, st, n_iter, es, em)
      CALL require(es == CD_AGG_OK .AND. cv, 'aggregate clump: step: '//TRIM(em))
      IF (es /= CD_AGG_OK) EXIT
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_EvalChannel(agg, 10, z, es, em)
      zmax_dev = MAX(zmax_dev, ABS(z - z0))
    END DO
    CALL require(zmax_dev < 1.0e-3_wp, 'aggregate clump: the clump stays at rest on the seabed')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_flat_aggregate

  SUBROUTINE case_fold_note_aggregate()
    !! The hairpin of an EI = 0 line through the OpenFAST aggregate: kept, and named in the
    !! line's initialisation report.
    TYPE(CD_AGG_ModuleType) :: agg
    INTEGER :: es, u, line_id
    REAL(wp) :: ten, force(3), incl, decl, azi
    CHARACTER(512) :: em
    CHARACTER(1024) :: note

    OPEN (NEWUNIT=u, FILE='point_seabed_hairpin.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---'
    WRITE (u, '(A)') 'Chain pushed past its anchor by a 2 m/s current on a frictionless seabed'
    WRITE (u, '(A)') '--------------------- LINE TYPES ---'
    WRITE (u, '(A)') 'TypeName Diam MassDenInAir EA BA EI Cd_n Cd_t Ca_n Ca_t'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'chain 0.216 286.56 1229760000 -1 0 1.37 0.64 1 0'
    WRITE (u, '(A)') '--------------------- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Coupled 0 0 -10 0 0 0 0'
    WRITE (u, '(A)') '2 Fixed 450 0 -100 0 0 0 0'
    WRITE (u, '(A)') '--------------------- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '--------------------- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 chain 520 52'
    WRITE (u, '(A)') '--------------------- OPTIONS ---'
    WRITE (u, '(A)') '100.0 WtrDpth'
    WRITE (u, '(A)') 'uniform 2.0 0 0 current'
    WRITE (u, '(A)') '--------------------- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1'
    WRITE (u, '(A)') '--------------------- need this line ---'
    CLOSE (u)
    CALL CD_AGG_Init_From_Deck(agg, 'point_seabed_hairpin.dat', 0.05_wp, es, em)
    CALL require(es == CD_AGG_OK, 'aggregate hairpin: init: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    CALL CD_AGG_GetInitLine(agg, 1, line_id, ten, force, incl, decl, azi, es, em, note=note)
    CALL require(es == CD_AGG_OK .AND. line_id == 1, 'aggregate hairpin: init report')
    CALL require(INDEX(note, 'the static layout carries the grounded run') > 0 .AND. &
                 INDEX(note, 'past its anchor') > 0 .AND. INDEX(note, 'the current drives the run') > 0, &
                 'aggregate hairpin: the line report carries the fold note: '//TRIM(note))
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_fold_note_aggregate

  SUBROUTINE write_clump_deck(path, tmax, scheme, sloped)
    !! The 50 t clump (deck point 3, Mass only) between two polyester lines on a 200 m seabed.
    CHARACTER(*), INTENT(IN) :: path, tmax, scheme
    LOGICAL, INTENT(IN) :: sloped
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ---'
    WRITE (u, '(A)') 'Heavy clump weight resting on the seabed between two polyester lines'
    WRITE (u, '(A)') '--------------------- LINE TYPES ---'
    WRITE (u, '(A)') 'TypeName Diam MassDenInAir EA BA EI Cd_n Cd_t Ca_n Ca_t'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'poly 0.172 31.916 200000000 -1 0 1.2 0.008 1 0'
    WRITE (u, '(A)') '--------------------- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m3) (m2) (-)'
    WRITE (u, '(A)') '1 Coupled 0 0 -10 0 0 0 0'
    IF (sloped) THEN
      WRITE (u, '(A)') '2 Fixed 800 0 -193 0 0 0 0'
    ELSE
      WRITE (u, '(A)') '2 Fixed 800 0 -200 0 0 0 0'
    END IF
    WRITE (u, '(A)') '3 Connect 450 0 -190 50000 0 0 0'
    WRITE (u, '(A)') '--------------------- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 3 -'
    WRITE (u, '(A)') '2 3 2 -'
    WRITE (u, '(A)') '--------------------- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 poly 520 52'
    WRITE (u, '(A)') '2 poly 360 36'
    WRITE (u, '(A)') '--------------------- OPTIONS ---'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    IF (sloped) THEN
      WRITE (u, '(A)') 'point_seabed_slope_bathy.txt bathymetryFile'
    ELSE
      WRITE (u, '(A)') '200.0 WtrDpth'
    END IF
    WRITE (u, '(A)') '100000.0 kBot'
    WRITE (u, '(A)') '1.0e4 cBot'
    WRITE (u, '(A)') '0.05 dtM'
    WRITE (u, '(A)') tmax//' TMax'
    IF (LEN_TRIM(scheme) > 0) WRITE (u, '(A)') scheme//' bodyScheme'
    WRITE (u, '(A)') '--------------------- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1'
    WRITE (u, '(A)') 'AnchTen1'
    WRITE (u, '(A)') 'FairTen2'
    WRITE (u, '(A)') 'AnchTen2'
    WRITE (u, '(A)') 'Point3Fx'
    WRITE (u, '(A)') 'Point3Fy'
    WRITE (u, '(A)') 'Point3Fz'
    WRITE (u, '(A)') 'Point3px'
    WRITE (u, '(A)') 'Point3py'
    WRITE (u, '(A)') 'Point3pz'
    WRITE (u, '(A)') '--------------------- need this line ---'
    CLOSE (u)
  END SUBROUTINE write_clump_deck

  SUBROUTINE read_rows(path, rows)
    !! The numeric rows of a driver .out file (two header lines), one column per row entry.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: rows(:, :)
    REAL(wp) :: row(NCOL)
    REAL(wp), ALLOCATABLE :: tmp(:, :)
    INTEGER :: u, ios, n
    CHARACTER(4096) :: buf
    ALLOCATE (rows(NCOL, 0))
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    CALL require(ios == 0, 'open '//path)
    IF (ios /= 0) RETURN
    READ (u, '(A)', IOSTAT=ios) buf
    READ (u, '(A)', IOSTAT=ios) buf
    n = 0
    DO
      READ (u, '(A)', IOSTAT=ios) buf
      IF (ios /= 0) EXIT
      READ (buf, *, IOSTAT=ios) row
      IF (ios /= 0) EXIT
      n = n + 1
      ALLOCATE (tmp(NCOL, n))
      tmp(:, 1:n - 1) = rows
      tmp(:, n) = row
      CALL MOVE_ALLOC(tmp, rows)
    END DO
    CLOSE (u)
  END SUBROUTINE read_rows

  SUBROUTINE require(ok, label)
    LOGICAL, INTENT(IN) :: ok
    CHARACTER(*), INTENT(IN) :: label
    IF (ok) RETURN
    failures = failures + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

END PROGRAM test_point_seabed_contact
