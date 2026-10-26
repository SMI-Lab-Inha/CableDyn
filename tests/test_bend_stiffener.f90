! File: tests/test_bend_stiffener.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_bend_stiffener
  !! Gates for bend stiffeners (STIFFENERS) on a finite-EI line:
  !!   1. Element properties (CD_Stiffener_Overlay) against closed forms: a cylindrical sleeve
  !!      over half an element gives the compliance average of the bending stiffness, the mean
  !!      sleeve mass, the displaced-area diameter and the mean outer diameter; the same sleeve
  !!      at End B lands on the mirrored element; elements it does not reach are unchanged bit
  !!      for bit.
  !!   2. A cone: the element values follow a 20 000-interval midpoint integration of the same
  !!      averages to 1e-10.
  !!   3. Deck grammar: malformed rows, a stiffener on an EI = 0 line and stiffeners longer than
  !!      the line fail closed by name; on a stock-order LINES row the stiffener sleeves the same
  !!      physical end, and the static profile is the native deck's; StfMom on a line without a
  !!      stiffener and AnchDefl at a pinned end fail closed, FairDefl at the clamped end runs.
  !!   4. Coupled route (the OpenFAST aggregate): the cable carries the stiffened root element,
  !!      and its initial fairlead load is the standalone static solve's.
  !! The static equilibrium against the planar elastica is tests/check_bend_stiffener.py.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK, CD_Stiffener_Overlay
  USE CableDyn_OpenFAST_Aggregate, ONLY: CD_AGG_ModuleType, CD_AGG_Init_From_Deck, CD_AGG_GetMovingPointMesh, &
                                         CD_AGG_End, CD_AGG_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
  INTEGER :: nfail
  nfail = 0

  CALL case_cylinder()
  CALL case_cone()
  CALL case_grammar()
  CALL case_coupled()

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'test_bend_stiffener: ', nfail, ' check(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'test_bend_stiffener: all checks passed'

CONTAINS

  SUBROUTINE require(cond, msg)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: msg
    IF (.NOT. cond) THEN
      nfail = nfail + 1
      WRITE (*, '(A)') 'FAIL: '//TRIM(msg)
    END IF
  END SUBROUTINE require

  LOGICAL FUNCTION close_rel(a, b, tol)
    REAL(wp), INTENT(IN) :: a, b, tol
    close_rel = ABS(a - b) <= tol*MAX(ABS(a), ABS(b))
  END FUNCTION close_rel

  SUBROUTINE case_cylinder()
    !! Four 1 m elements; a 1.5 m cylindrical sleeve (OD 0.5 m, ID 0.2 m) covers element 1 and
    !! half of element 2 from End A (or from End B: elements 4 and 3).
    REAL(wp), PARAMETER :: EI_C = 2.0e4_wp, M_C = 60.0_wp, D_C = 0.2_wp
    REAL(wp), PARAMETER :: OD = 0.5_wp, ID = 0.2_wp, EMOD = 5.0e7_wp, RHO = 1150.0_wp
    REAL(wp) :: stf(8, 1), l0(4), ei(4), mpl(4), dv(4), dh(4), ei_s, ei_half, m_s
    INTEGER :: iend, e_full, e_half, e_far(2)
    l0 = 1.0_wp
    ei_s = EMOD*PI/64.0_wp*(OD**4 - ID**4)
    m_s = RHO*PI/4.0_wp*(OD**2 - ID**2)
    ei_half = 1.0_wp/(0.5_wp/EI_C + 0.5_wp/(EI_C + ei_s))
    DO iend = 1, 2
      stf(:, 1) = [REAL(iend, wp), 0.0_wp, 1.5_wp, OD, OD, ID, EMOD, RHO]
      ei = EI_C; mpl = M_C; dv = D_C; dh = D_C
      CALL CD_Stiffener_Overlay(stf, l0, dv, ei=ei, mpl=mpl, volume_diameter=.TRUE.)
      CALL CD_Stiffener_Overlay(stf, l0, dh, volume_diameter=.FALSE.)
      IF (iend == 1) THEN
        e_full = 1; e_half = 2; e_far = [3, 4]
      ELSE
        e_full = 4; e_half = 3; e_far = [1, 2]
      END IF
      CALL require(close_rel(ei(e_full), EI_C + ei_s, 1.0e-13_wp), 'cylinder: covered element EI')
      CALL require(close_rel(ei(e_half), ei_half, 1.0e-13_wp), 'cylinder: half-covered element takes the '// &
                   'compliance average')
      CALL require(close_rel(mpl(e_full), M_C + m_s, 1.0e-13_wp), 'cylinder: covered element mass')
      CALL require(close_rel(mpl(e_half), M_C + 0.5_wp*m_s, 1.0e-13_wp), 'cylinder: half-covered mass')
      CALL require(close_rel(dv(e_full), SQRT(D_C**2 + OD**2 - ID**2), 1.0e-13_wp), &
                   'cylinder: displaced-area diameter of cable plus solid sleeve')
      CALL require(close_rel(dv(e_half), SQRT(D_C**2 + 0.5_wp*(OD**2 - ID**2)), 1.0e-13_wp), &
                   'cylinder: half-covered displaced-area diameter')
      CALL require(close_rel(dh(e_full), OD, 1.0e-13_wp), 'cylinder: outer diameter for drag')
      CALL require(close_rel(dh(e_half), 0.5_wp*(OD + D_C), 1.0e-13_wp), 'cylinder: half-covered outer diameter')
      CALL require(ALL(ei(e_far) == EI_C) .AND. ALL(mpl(e_far) == M_C) .AND. ALL(dv(e_far) == D_C) .AND. &
                   ALL(dh(e_far) == D_C), 'cylinder: elements the sleeve does not reach are unchanged')
      ! an EI = 0 section: the covered element bends with the sleeve alone; the half-covered one
      ! keeps EI = 0, its bare half having an infinite compliance
      ei = 0.0_wp; mpl = M_C; dv = D_C
      CALL CD_Stiffener_Overlay(stf, l0, dv, ei=ei, mpl=mpl, volume_diameter=.TRUE.)
      CALL require(close_rel(ei(e_full), ei_s, 1.0e-13_wp), 'cylinder: covered EI = 0 element takes the sleeve')
      CALL require(ei(e_half) == 0.0_wp, 'cylinder: half-covered EI = 0 element keeps EI = 0')
    END DO
  END SUBROUTINE case_cylinder

  SUBROUTINE case_cone()
    !! A 2.5 m cone (OD 0.6 -> 0.22 m) from End A over five 0.7 m elements: the fourth element
    !! (2.1 - 2.8 m) holds the cone's tip; the fifth is bare.
    REAL(wp), PARAMETER :: EI_C = 2.0e4_wp, M_C = 60.0_wp, D_C = 0.2_wp
    REAL(wp), PARAMETER :: OD0 = 0.6_wp, OD1 = 0.22_wp, ID = 0.2_wp, EMOD = 5.0e7_wp, RHO = 1150.0_wp
    INTEGER, PARAMETER :: NQ = 20000
    REAL(wp) :: stf(8, 1), l0(5), ei(5), mpl(5), dv(5), dh(5), a, s, od, comp, mass, area, outer, h, tol
    INTEGER :: e, k
    l0 = 0.7_wp
    stf(:, 1) = [1.0_wp, 0.0_wp, 2.5_wp, OD0, OD1, ID, EMOD, RHO]
    ei = EI_C; mpl = M_C; dv = D_C; dh = D_C
    CALL CD_Stiffener_Overlay(stf, l0, dv, ei=ei, mpl=mpl, volume_diameter=.TRUE.)
    CALL CD_Stiffener_Overlay(stf, l0, dh, volume_diameter=.FALSE.)
    DO e = 1, 4
      a = 0.7_wp*REAL(e - 1, wp)
      h = 0.7_wp/REAL(NQ, wp)
      comp = 0.0_wp; mass = 0.0_wp; area = 0.0_wp; outer = 0.0_wp
      DO k = 1, NQ
        s = a + (REAL(k, wp) - 0.5_wp)*h
        IF (s < 2.5_wp) THEN
          od = OD0 + s/2.5_wp*(OD1 - OD0)
          comp = comp + h/(EI_C + EMOD*PI/64.0_wp*(od**4 - ID**4))
          mass = mass + h*RHO*PI/4.0_wp*(od**2 - ID**2)
          area = area + h*(od**2 - ID**2)
          outer = outer + h*MAX(od, D_C)
        ELSE
          comp = comp + h/EI_C
          outer = outer + h*D_C
        END IF
      END DO
      ! the cone tip inside element 4 limits the midpoint rule there to O(h/l) = 5e-5
      tol = MERGE(2.0e-4_wp, 1.0e-6_wp, e == 4)
      CALL require(close_rel(ei(e), 0.7_wp/comp, tol), 'cone: compliance-average EI of element')
      CALL require(close_rel(mpl(e), M_C + mass/0.7_wp, tol), 'cone: mean mass of element')
      CALL require(close_rel(dv(e), SQRT(D_C**2 + area/0.7_wp), tol), 'cone: displaced-area diameter')
      CALL require(close_rel(dh(e), outer/0.7_wp, tol), 'cone: mean outer diameter')
      IF (e < 4) CALL require(close_rel(ei(e), 0.7_wp/comp, 1.0e-10_wp), &
                              'cone: Gauss average matches the fine integration inside the cone')
    END DO
    CALL require(ei(5) == EI_C .AND. mpl(5) == M_C .AND. dv(5) == D_C .AND. dh(5) == D_C, &
                 'cone: bare elements unchanged')
  END SUBROUTINE case_cone

  SUBROUTINE write_deck(path, stiffener_rows, ei, stock, channel)
    !! A 60 m cable clamped horizontal at a 20 m deep hang-off and fixed 30 m lower; stock: the
    !! 7-column LINES row written anchor first, so the parser swaps the ends.
    CHARACTER(*), INTENT(IN) :: path
    CHARACTER(*), INTENT(IN) :: stiffener_rows(:)
    REAL(wp), INTENT(IN) :: ei
    LOGICAL, INTENT(IN) :: stock
    CHARACTER(*), INTENT(IN), OPTIONAL :: channel
    INTEGER :: u, k
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ------------------------------------'
    WRITE (u, '(A)') 'Bend stiffener test deck'
    WRITE (u, '(A)') '--------------------- LINE TYPES ---------------------------------------------'
    WRITE (u, '(A)') 'TypeName Diam MassDenInAir EA BA/-zeta EI Cd_n Cd_t Ca_n Ca_t'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (N-s/-) (N-m^2) (-) (-) (-) (-)'
    WRITE (u, '(A,ES12.4,A)') 'cable 0.20 60.0 5.0e8 -1.0 ', ei, ' 1.2 0.0 1.0 0.0'
    WRITE (u, '(A)') '--------------------- POINTS -------------------------------------------------'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)'
    WRITE (u, '(A)') '1 Coupled 0.0 0.0 -20.0 0 0 0 0'
    WRITE (u, '(A)') '2 Fixed 51.90 0.0 -50.0 0 0 0 0'
    WRITE (u, '(A)') '--------------------- LINES --------------------------------------------------'
    IF (stock) THEN
      WRITE (u, '(A)') 'ID LineType AttachA AttachB UnstrLen NumSegs Outputs'
      WRITE (u, '(A)') '(-) (-) (-) (-) (m) (-) (-)'
      WRITE (u, '(A)') '1 cable 2 1 60.0 120 -'
    ELSE
      WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
      WRITE (u, '(A)') '(-) (-) (-) (-)'
      WRITE (u, '(A)') '1 1 2 -'
      WRITE (u, '(A)') '--------------------- SECTIONS -----------------------------------------------'
      WRITE (u, '(A)') 'LineID LineType Length NumSegs'
      WRITE (u, '(A)') '(-) (-) (m) (-)'
      WRITE (u, '(A)') '1 cable 60.0 120'
    END IF
    IF (SIZE(stiffener_rows) > 0) THEN
      WRITE (u, '(A)') '--------------------- STIFFENERS ---------------------------------------------'
      WRITE (u, '(A)') 'LineID End Length RootOD TipOD InnerD E Density'
      WRITE (u, '(A)') '(-) (-) (m) (m) (m) (m) (Pa) (kg/m^3)'
      DO k = 1, SIZE(stiffener_rows)
        WRITE (u, '(A)') TRIM(stiffener_rows(k))
      END DO
    END IF
    IF (ei > 0.0_wp) THEN
      WRITE (u, '(A)') '--------------------- END CONNECTIONS ----------------------------------------'
      WRITE (u, '(A)') 'LineID End Stiffness EzX EzY EzZ'
      WRITE (u, '(A)') '(-) (-) (N-m/rad) (-) (-) (-)'
      IF (stock) THEN
        WRITE (u, '(A)') '1 B Rigid -1.0 0.0 0.0'
      ELSE
        WRITE (u, '(A)') '1 A Rigid 1.0 0.0 0.0'
      END IF
    END IF
    WRITE (u, '(A)') '--------------------- OPTIONS ------------------------------------------------'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '0.05 dtM'
    WRITE (u, '(A)') '0.0 TMax'
    WRITE (u, '(A)') '--------------------- OUTPUTS ------------------------------------------------'
    WRITE (u, '(A)') '"FairTen1"'
    IF (PRESENT(channel)) WRITE (u, '(A)') '"'//channel//'"'
    WRITE (u, '(A)') 'END'
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

  SUBROUTINE read_static(path, prof, n)
    !! Columns 3-9 of a .static.out (arc, X, Y, Z, tension, curvature, moment), node by node.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: prof(:, :)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(1024) :: rec
    INTEGER :: u, ios, idum(2)
    ALLOCATE (prof(7, 2000))
    n = 0
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    READ (u, '(A)') rec
    READ (u, '(A)') rec
    READ (u, '(A)') rec
    DO
      READ (u, '(A)', IOSTAT=ios) rec
      IF (ios /= 0) EXIT
      IF (LEN_TRIM(rec) == 0) CYCLE
      n = n + 1
      READ (rec, *) idum, prof(:, n)
    END DO
    CLOSE (u)
  END SUBROUTINE read_static

  SUBROUTINE case_grammar()
    CHARACTER(80) :: rows(2), none(0)
    CHARACTER(512) :: msg
    LOGICAL :: ok
    REAL(wp), ALLOCATABLE :: pa(:, :), pb(:, :)
    INTEGER :: na, nb
    rows(1) = '1 A 0.5 0.60 0.60 0.20 5.0e7 1150'
    rows(2) = '1 A 4.5 0.60 0.22 0.20 5.0e7 1150'
    CALL write_deck('stf_native.dat', rows, 2.0e4_wp, .FALSE.)
    CALL run('stf_native.dat', 'stf_native', ok, msg)
    CALL require(ok, 'grammar: native deck solves: '//TRIM(msg))
    ! stock order: the deck's End B is the hang-off; the stiffener sleeves the same physical end
    rows(1) = '1 B 0.5 0.60 0.60 0.20 5.0e7 1150'
    rows(2) = '1 B 4.5 0.60 0.22 0.20 5.0e7 1150'
    CALL write_deck('stf_stock.dat', rows, 2.0e4_wp, .TRUE.)
    CALL run('stf_stock.dat', 'stf_stock', ok, msg)
    CALL require(ok, 'grammar: stock-order deck solves: '//TRIM(msg))
    CALL read_static('stf_native.static.out', pa, na)
    CALL read_static('stf_stock.static.out', pb, nb)
    CALL require(na == nb .AND. na > 0, 'grammar: both decks write a static profile')
    IF (na == nb .AND. na > 0) THEN
      CALL require(MAXVAL(ABS(pa(2:4, :na) - pb(2:4, :nb))) < 1.0e-6_wp .AND. &
                   MAXVAL(ABS(pa(7, :na) - pb(7, :nb))) < 1.0e-6_wp*MAXVAL(ABS(pa(7, :na))), &
                   'grammar: on a stock-order row the stiffener sleeves the same physical end')
    END IF
    rows(1) = '1 A 0.5 0.60 0.60 0.20 5.0e7'
    CALL write_deck('stf_bad0.dat', rows(1:1), 2.0e4_wp, .FALSE.)
    CALL run('stf_bad0.dat', 'stf_bad0', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'LineID End Length') > 0, 'grammar: a seven-column row: '//TRIM(msg))
    rows(1) = '1 C 0.5 0.60 0.60 0.20 5.0e7 1150'
    CALL write_deck('stf_bad1.dat', rows(1:1), 2.0e4_wp, .FALSE.)
    CALL run('stf_bad1.dat', 'stf_bad1', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'End must be A or B') > 0, 'grammar: an unknown end: '//TRIM(msg))
    rows(1) = '1 A 0.5 0.60 0.60 0.70 5.0e7 1150'
    CALL write_deck('stf_bad2.dat', rows(1:1), 2.0e4_wp, .FALSE.)
    CALL run('stf_bad2.dat', 'stf_bad2', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'InnerD must be smaller') > 0, 'grammar: InnerD above OD: '//TRIM(msg))
    rows(1) = '1 A 0.5 0.60 0.60 0.20 0.0 1150'
    CALL write_deck('stf_bad3.dat', rows(1:1), 2.0e4_wp, .FALSE.)
    CALL run('stf_bad3.dat', 'stf_bad3', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'must be > 0') > 0, 'grammar: a zero modulus: '//TRIM(msg))
    rows(1) = '2 A 0.5 0.60 0.60 0.20 5.0e7 1150'
    CALL write_deck('stf_bad4.dat', rows(1:1), 2.0e4_wp, .FALSE.)
    CALL run('stf_bad4.dat', 'stf_bad4', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'undefined LINE') > 0, 'grammar: an unknown line: '//TRIM(msg))
    rows(1) = '1 A 40.0 0.60 0.60 0.20 5.0e7 1150'
    rows(2) = '1 B 25.0 0.60 0.60 0.20 5.0e7 1150'
    CALL write_deck('stf_bad5.dat', rows, 2.0e4_wp, .FALSE.)
    CALL run('stf_bad5.dat', 'stf_bad5', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'longer than the line') > 0, 'grammar: stiffeners overlap: '//TRIM(msg))
    rows(1) = '1 A 0.5 0.60 0.60 0.20 5.0e7 1150'
    CALL write_deck('stf_bad6.dat', rows(1:1), 0.0_wp, .FALSE.)
    CALL run('stf_bad6.dat', 'stf_bad6', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'not a finite-EI line') > 0, 'grammar: an EI = 0 line: '//TRIM(msg))
    ! no STIFFENERS table: the same deck still solves (the table is optional)
    CALL write_deck('stf_none.dat', none, 2.0e4_wp, .FALSE.)
    CALL run('stf_none.dat', 'stf_none', ok, msg)
    CALL require(ok, 'grammar: the deck without stiffeners solves: '//TRIM(msg))
    ! the stiffener moment needs a stiffener; the deflection angle a held end
    CALL write_deck('stf_ch1.dat', none, 2.0e4_wp, .FALSE., channel='StfMom1N3')
    CALL run('stf_ch1.dat', 'stf_ch1', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'has no STIFFENERS') > 0, 'grammar: StfMom without a stiffener: '//TRIM(msg))
    rows(1) = '1 A 0.5 0.60 0.60 0.20 5.0e7 1150'
    CALL write_deck('stf_ch2.dat', rows(1:1), 2.0e4_wp, .FALSE., channel='AnchDefl1')
    CALL run('stf_ch2.dat', 'stf_ch2', ok, msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'END CONNECTIONS row at End B') > 0, &
                 'grammar: AnchDefl at a pinned End B: '//TRIM(msg))
    CALL write_deck('stf_ch3.dat', rows(1:1), 2.0e4_wp, .FALSE., channel='FairDefl1')
    CALL run('stf_ch3.dat', 'stf_ch3', ok, msg)
    CALL require(ok, 'grammar: FairDefl at the clamped End A: '//TRIM(msg))
  END SUBROUTINE case_grammar

  SUBROUTINE case_coupled()
    !! The native deck of case_grammar on the coupled route: its root element, wholly inside the
    !! 0.5 m cylindrical root (OD 0.6 m, ID 0.2 m), has EI_c + E pi/64 (OD^4 - ID^4), and the
    !! initial fairlead load matches the standalone static fairlead tension.
    TYPE(CD_AGG_ModuleType), SAVE :: agg
    REAL(wp), PARAMETER :: EI_ROOT = 2.0e4_wp + 5.0e7_wp*PI/64.0_wp*(0.6_wp**4 - 0.2_wp**4)
    REAL(wp) :: r0(3, 1), vel(3, 1), acc(3, 1), load(3, 1)
    REAL(wp), ALLOCATABLE :: prof(:, :)
    CHARACTER(512) :: em
    INTEGER :: es, ne, n
    CALL CD_AGG_Init_From_Deck(agg, 'stf_native.dat', 0.05_wp, es, em)
    CALL require(es == CD_AGG_OK, 'coupled: the stiffener deck initializes: '//TRIM(em))
    IF (es /= CD_AGG_OK) RETURN
    ! internal element order runs from End B: the root element at End A is the last
    ne = SIZE(agg%cables(1)%line%EI)
    CALL require(close_rel(agg%cables(1)%line%EI(ne), EI_ROOT, 1.0e-12_wp), &
                 'coupled: the root element carries the stiffener')
    CALL CD_AGG_GetMovingPointMesh(agg, r0, vel, acc, load, es, em)
    CALL read_static('stf_native.static.out', prof, n)
    IF (n > 0) CALL require(close_rel(NORM2(load(:, 1)), prof(5, 1), 1.0e-6_wp), &
                            'coupled: the initial fairlead load is the standalone static tension')
    CALL CD_AGG_End(agg, es, em)
  END SUBROUTINE case_coupled

END PROGRAM test_bend_stiffener
