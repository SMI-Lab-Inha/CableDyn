! File: tests/test_static_dry_weight.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_static_dry_weight
  !! The static start of an EI = 0 line held above the still-water level gives the part above
  !! the surface its dry weight, as the dynamic model's buoyancy recovery does, on the
  !! standalone and the coupled route alike (CD_Static_Line_Equilibrium with the dry-part
  !! arguments). Gates:
  !!   1. A vertical chain hung from 10 m above the surface to a point 90 m below it (99.9 m,
  !!      taut): the vertical loads of its two ends sum to its weight, w_sub times the wet
  !!      length plus the dry weight per metre times the dry length (the dry length 10 m over
  !!      the stretch of the top element), to 1e-6 (measured 1.5e-8). With the buoyancy kept
  !!      above the surface the sum read 1.4e-2 low.
  !!   2. A chain hung 5 m above the surface onto a 200 m seabed (850 m, grounded tail): the
  !!      coupled aggregate's static start gives the standalone driver's fairlead load to the
  !!      driver's 8 output digits, 2e-7 (measured 1.4e-8; before 2.4e-4: the coupled load was
  !!      the dynamic reaction of a start that was not the equilibrium of the dynamics, whose
  !!      buoyancy recovery gives the dry part its dry weight).
  !!   3. The same chain held still on the coupled route for 1 s (40 steps of 0.025 s): the
  !!      fairlead load stays at its static value to 1e-7 (measured 1e-10; before it moved by
  !!      9.7e-4 in the first step and by 1.4e-3 at 1 s, ringing about the heavier
  !!      equilibrium).
  !!   4. A finite-EI cable hung 5 m above the surface (its statics already carry the dry
  !!      weight): coupled and standalone starts agree to 2e-7 (measured 7e-8; a guard).
  !! A deck whose lines are all submerged is solved once, bit-identically (the example A/B).
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_OpenFAST_Aggregate, ONLY: CD_AGG_ModuleType, CD_AGG_Init_From_Deck, CD_AGG_GetMovingPointMesh, &
                                         CD_AGG_NMovingPoints, CD_AGG_CalcOutput, CD_AGG_Step_Moving, CD_AGG_End, &
                                         CD_AGG_OK
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: G = 9.80665_wp, RHOW = 1025.0_wp, PI = 3.14159265358979323846_wp
  REAL(wp), PARAMETER :: D_CHAIN = 0.333_wp, M_CHAIN = 685.0_wp, EA_CHAIN = 3.27e9_wp
  INTEGER :: nfail
  REAL(wp) :: f_drv(3), f_agg(3), f2(3), w_dry, w_sub, l_dry, eps_top, w_true, err, dmax

  nfail = 0
  w_dry = M_CHAIN*G
  w_sub = w_dry - RHOW*G*0.25_wp*PI*D_CHAIN**2

  ! 1. vertical chain across the surface
  CALL write_deck('dry_vertical.dat', 'main 0.333 685.00 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27', &
                  '0.0 0.0 10.0', '0.0 0.0 -90.0', '1 main 99.9 100', 200.0_wp)
  CALL driver_forces('dry_vertical', f_drv, f2)
  ! the dry length: 10 m of the top element's stretch (the tension varies by 0.1 % over it)
  eps_top = -f_drv(3)/EA_CHAIN
  l_dry = 10.0_wp/(1.0_wp + eps_top)
  w_true = w_sub*(99.9_wp - l_dry) + w_dry*l_dry
  err = ABS(f_drv(3) + f2(3) + w_true)/w_true
  WRITE (*, '(A,ES12.4,A,ES12.4,A,ES10.2)') '  vertical chain: sum of the end loads ', -(f_drv(3) + f2(3)), &
    ' N, weight ', w_true, ' N, rel. error ', err
  WRITE (*, '(A,ES10.2)') '    (buoyancy kept above the surface would read it low by ', &
    (w_dry - w_sub)*l_dry/w_true
  CALL check(err <= 1.0e-6_wp, 'vertical chain: the end loads carry the dry weight above the surface')

  ! 2. and 3. chain from 5 m above the surface onto the seabed
  CALL write_deck('dry_chain.dat', 'main 0.333 685.00 3.27e9 -1.0 0.0 2.0 0.4 0.82 0.27', &
                  '-58.0 0.0 5.0', '-837.6 0.0 -200.0', '1 main 850.0 50', 200.0_wp)
  CALL driver_forces('dry_chain', f_drv, f2)
  CALL coupled_force('dry_chain.dat', 40, f_agg, dmax)
  err = MAXVAL(ABS(f_agg - f_drv))/NORM2(f_drv)
  WRITE (*, '(A,3ES16.8)') '  chain, standalone fairlead load ', f_drv
  WRITE (*, '(A,3ES16.8,A,ES10.2)') '  chain, coupled fairlead load    ', f_agg, '  rel. ', err
  CALL check(err <= 2.0e-7_wp, 'chain above the surface: coupled and standalone static starts agree')
  WRITE (*, '(A,ES10.2)') '  chain held still 1 s on the coupled route: largest rel. load change ', dmax
  CALL check(dmax <= 1.0e-7_wp, 'chain above the surface: the static start is the equilibrium of the march')

  ! 4. finite-EI cable from 5 m above the surface
  CALL write_deck('dry_cable.dat', 'bare 0.16 36.70 4.69e8 0.0 1.99e4 1.2 0.1 1.0 0.0', &
                  '30.0 0.0 5.0', '500.0 0.0 -200.0', '1 bare 600.0 240', 200.0_wp)
  CALL driver_forces('dry_cable', f_drv, f2)
  CALL coupled_force('dry_cable.dat', 0, f_agg, dmax)
  err = MAXVAL(ABS(f_agg - f_drv))/NORM2(f_drv)
  WRITE (*, '(A,ES10.2)') '  finite-EI cable above the surface: coupled vs standalone rel. ', err
  CALL check(err <= 2.0e-7_wp, 'finite-EI cable above the surface: coupled and standalone static starts agree')

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAILED: ', nfail, ' check(s)'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: static dry weight above the surface'

CONTAINS

  SUBROUTINE check(ok, what)
    LOGICAL, INTENT(IN) :: ok
    CHARACTER(*), INTENT(IN) :: what
    IF (ok) THEN
      WRITE (*, '(2A)') '  ok: ', what
    ELSE
      WRITE (*, '(2A)') '  FAIL: ', what
      nfail = nfail + 1
    END IF
  END SUBROUTINE check

  SUBROUTINE write_deck(path, type_row, pa, pb, section_row, depth)
    !! One line from a Coupled point pa (End A) to a Fixed point pb; the standalone static start
    !! at TMax = 0, the time step for the coupled route.
    CHARACTER(*), INTENT(IN) :: path, type_row, pa, pb, section_row
    REAL(wp), INTENT(IN) :: depth
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='replace', ACTION='write')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ------------------------------------'
    WRITE (u, '(A)') 'test_static_dry_weight: a line held above the still-water level'
    WRITE (u, '(A)') '--------------------- LINE TYPES ---------------------------------------------'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') type_row
    WRITE (u, '(A)') '--------------------- POINTS ------------------------------------------------'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Coupled '//pa
    WRITE (u, '(A)') '2 Fixed '//pb
    WRITE (u, '(A)') '--------------------- LINES -------------------------------------------------'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '--------------------- SECTIONS ----------------------------------------------'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') section_row
    WRITE (u, '(A)') '--------------------- OPTIONS -----------------------------------------------'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '3.0e6 kBot'
    WRITE (u, '(A)') '3.0e5 cBot'
    WRITE (u, '(A)') '0.025 dtM'
    WRITE (u, '(F8.2,A)') depth, ' WtrDpth'
    WRITE (u, '(A)') '0.0 TMax'
    WRITE (u, '(A)') '--------------------- OUTPUTS -----------------------------------------------'
    WRITE (u, '(A)') '"Point1Fx"'
    WRITE (u, '(A)') '"Point1Fy"'
    WRITE (u, '(A)') '"Point1Fz"'
    WRITE (u, '(A)') '"Point2Fx"'
    WRITE (u, '(A)') '"Point2Fy"'
    WRITE (u, '(A)') '"Point2Fz"'
    WRITE (u, '(A)') 'END'
    WRITE (u, '(A)') '--------------------- need this line ----------------------------------------'
    CLOSE (u)
  END SUBROUTINE write_deck

  SUBROUTINE driver_forces(root, fa, fb)
    !! The standalone static start (TMax = 0): the loads on points 1 and 2 (the t = 0 row).
    CHARACTER(*), INTENT(IN) :: root
    REAL(wp), INTENT(OUT) :: fa(3), fb(3)
    LOGICAL :: conv
    INTEGER :: es, u, ios
    CHARACTER(1024) :: em
    CHARACTER(4096) :: line
    REAL(wp) :: t
    fa = 0.0_wp
    fb = 0.0_wp
    CALL CD_Run_Deck_Driver(root//'.dat', root, conv, es, em)
    CALL check(es == CD_DECKDRV_OK .AND. conv, root//': standalone static start: '//TRIM(em))
    IF (es /= CD_DECKDRV_OK) RETURN
    OPEN (NEWUNIT=u, FILE=root//'.out', STATUS='old', ACTION='read')
    DO
      READ (u, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF (line(1:1) == '#' .OR. INDEX(line, 'Time') > 0 .OR. LEN_TRIM(line) == 0) CYCLE
      READ (line, *) t, fa, fb
      EXIT
    END DO
    CLOSE (u)
  END SUBROUTINE driver_forces

  SUBROUTINE coupled_force(path, nstep, f0, dmax_rel)
    !! The coupled aggregate's static start (OpenFAST host environment): the fairlead load, and
    !! the largest relative change of it over nstep still steps of 0.025 s.
    CHARACTER(*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: nstep
    REAL(wp), INTENT(OUT) :: f0(3), dmax_rel
    TYPE(CD_AGG_ModuleType) :: agg
    REAL(wp), ALLOCATABLE :: pos(:, :), vel(:, :), acc(:, :), ld(:, :), p0(:, :)
    INTEGER :: es, ncp, n, it
    LOGICAL :: conv, stalled
    CHARACTER(1024) :: em
    f0 = 0.0_wp
    dmax_rel = HUGE(1.0_wp)
    CALL CD_AGG_Init_From_Deck(agg, path, 0.025_wp, es, em, env_gravity=G, env_rho_water=RHOW, &
                               env_wtrdpth=200.0_wp)
    CALL check(es == CD_AGG_OK, path//': coupled static start: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    ncp = CD_AGG_NMovingPoints(agg, es, em)
    ALLOCATE (pos(3, ncp), vel(3, ncp), acc(3, ncp), ld(3, ncp), p0(3, ncp))
    CALL CD_AGG_CalcOutput(agg, es, em)
    CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, ld, es, em)
    f0 = ld(:, 1)
    p0 = pos
    vel = 0.0_wp
    acc = 0.0_wp
    dmax_rel = 0.0_wp
    DO n = 1, nstep
      CALL CD_AGG_Step_Moving(agg, 0.025_wp, p0, vel, acc, conv, stalled, it, es, em)
      IF (es /= CD_AGG_OK) THEN
        CALL check(.FALSE., path//': coupled still step: '//TRIM(em))
        EXIT
      END IF
      CALL CD_AGG_CalcOutput(agg, es, em)
      CALL CD_AGG_GetMovingPointMesh(agg, pos, vel, acc, ld, es, em)
      vel = 0.0_wp
      acc = 0.0_wp
      dmax_rel = MAX(dmax_rel, NORM2(ld(:, 1) - f0)/NORM2(f0))
    END DO
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE coupled_force
END PROGRAM test_static_dry_weight
