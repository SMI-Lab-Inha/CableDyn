! File: tests/test_bending_damping_deck.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_bending_damping_deck
  !! Bending damping on a deck (LINE TYPES EI|Bb, EI|-zeta|f and EI|-zeta|f1|f2): an 80 m power cable hanging from
  !! a clamped hang-off to a seabed anchor, in a regular sea.
  !!   E1 EI|-zeta|f and the EI|Bb it resolves to (Bb = zeta EI/(pi f)) give the same run, to the
  !!      printed digits;
  !!   E2 the damping lowers the hang-off curvature range against the undamped deck, and EI|0 is
  !!      the undamped deck exactly;
  !!   E3 bending damping on an EI = 0 type (as Bb or as a ratio) and malformed EI columns are
  !!      refused by name;
  !!   E4 EI|-zeta|f1|f2 runs; its mass-proportional half acts (the run differs from the EI|Bb deck
  !!      of the same Bb = zeta EI/(pi (f1 + f2))) and lowers the hang-off curvature range further;
  !!   E5 the two-frequency form with a positive ratio, a zero frequency or on an EI = 0 type is
  !!      refused by name.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp, EI = 2.0e4_wp, ZETA = 0.05_wp, FREQ = 0.5_wp
  INTEGER :: nfail = 0

  CALL check_runs()
  CALL check_pair()
  CALL check_refusals()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: bending-damping deck gates'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE write_deck(path, ei_col, tmax, rope)
    CHARACTER(*), INTENT(IN) :: path, ei_col
    REAL(wp), INTENT(IN) :: tmax
    LOGICAL, INTENT(IN), OPTIONAL :: rope
    INTEGER :: u
    LOGICAL :: on_rope
    on_rope = .FALSE.
    IF (PRESENT(rope)) on_rope = rope
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '--------------------- CableDyn Input File ------------------------------------'
    WRITE (u, '(A)') 'hang-off cable with bending damping'
    WRITE (u, '(A)') '--------------------- LINE TYPES ---------------------------------------------'
    WRITE (u, '(A)') 'TypeName Diam MassDenInAir EA BA EI Cd_n Cd_t Ca_n Ca_t'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (N-m2) (-) (-) (-) (-)'
    IF (on_rope) THEN
      WRITE (u, '(A)') 'pc 0.15 40.0 5.0e8 0.0 2.0e4 1.2 0.1 1.0 0.0'
      WRITE (u, '(A)') 'rope 0.10 20.0 5.0e7 0.0 '//TRIM(ei_col)//' 1.2 0.1 1.0 0.0'
    ELSE
      WRITE (u, '(A)') 'pc 0.15 40.0 5.0e8 0.0 '//TRIM(ei_col)//' 1.2 0.1 1.0 0.0'
    END IF
    WRITE (u, '(A)') '--------------------- POINTS -------------------------------------------------'
    WRITE (u, '(A)') 'ID Type X Y Z'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m)'
    WRITE (u, '(A)') '1 Coupled 0.0 0.0 -5.0'
    WRITE (u, '(A)') '2 Fixed 60.0 0.0 -50.0'
    WRITE (u, '(A)') '--------------------- LINES --------------------------------------------------'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 -'
    WRITE (u, '(A)') '--------------------- SECTIONS -----------------------------------------------'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 pc 80.0 320'
    WRITE (u, '(A)') '--------------------- END CONNECTIONS ----------------------------------------'
    WRITE (u, '(A)') 'LineID End Stiffness EzX EzY EzZ'
    WRITE (u, '(A)') '(-) (-) (N-m/rad) (-) (-) (-)'
    WRITE (u, '(A)') '1 A Rigid 0.3746 0.0 -0.9272'
    WRITE (u, '(A)') '--------------------- OPTIONS ------------------------------------------------'
    WRITE (u, '(A)') '9.80665 g'
    WRITE (u, '(A)') '1025.0 rhoW'
    WRITE (u, '(A)') '50.0 WtrDpth'
    WRITE (u, '(A)') '0.02 dtM'
    WRITE (u, '(A,F8.3,A)') ' ', tmax, ' TMax'
    WRITE (u, '(A)') 'airy 3.0 8.0 0.0 waves'
    WRITE (u, '(A)') '5.0 rampTime'
    WRITE (u, '(A)') '--------------------- OUTPUTS ------------------------------------------------'
    WRITE (u, '(A)') '"FairTen1"'
    WRITE (u, '(A)') '"BendMom1N1"'
    WRITE (u, '(A)') '"Curv1N1"'
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

  SUBROUTINE read_out(path, tab, n)
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: tab(:, :)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(1024) :: rec
    INTEGER :: u, ios
    ALLOCATE (tab(4, 20000))
    n = 0
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    READ (u, '(A)') rec
    READ (u, '(A)') rec
    DO
      READ (u, '(A)', IOSTAT=ios) rec
      IF (ios /= 0) EXIT
      IF (LEN_TRIM(rec) == 0) CYCLE
      n = n + 1
      READ (rec, *) tab(:, n)
    END DO
    CLOSE (u)
  END SUBROUTINE read_out

  SUBROUTINE check_runs()
    REAL(wp), ALLOCATABLE :: tz(:, :), tb(:, :), tp(:, :), t0(:, :)
    INTEGER :: nz, nb, np, n0
    REAL(wp) :: bb, rz, rp, r0
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CHARACTER(64) :: col
    bb = ZETA*EI/(PI*FREQ)
    WRITE (col, '(ES22.15)') bb
    CALL write_deck('bdd_zeta.dat', '2.0e4|-0.05|0.5', 20.0_wp)
    CALL run('bdd_zeta.dat', 'bdd_zeta', ok, msg)
    CALL require(ok, 'E1: the EI|-zeta|f deck runs: '//TRIM(msg))
    CALL write_deck('bdd_bb.dat', '2.0e4|'//TRIM(ADJUSTL(col)), 20.0_wp)
    CALL run('bdd_bb.dat', 'bdd_bb', ok, msg)
    CALL require(ok, 'E1: the EI|Bb deck runs: '//TRIM(msg))
    CALL write_deck('bdd_plain.dat', '2.0e4', 20.0_wp)
    CALL run('bdd_plain.dat', 'bdd_plain', ok, msg)
    CALL require(ok, 'E2: the undamped deck runs: '//TRIM(msg))
    CALL write_deck('bdd_zero.dat', '2.0e4|0', 20.0_wp)
    CALL run('bdd_zero.dat', 'bdd_zero', ok, msg)
    CALL require(ok, 'E2: the EI|0 deck runs: '//TRIM(msg))
    CALL read_out('bdd_zeta.out', tz, nz)
    CALL read_out('bdd_bb.out', tb, nb)
    CALL read_out('bdd_plain.out', tp, np)
    CALL read_out('bdd_zero.out', t0, n0)
    CALL require(nz > 10 .AND. nz == nb .AND. nz == np .AND. nz == n0, 'E1/E2: the runs wrote matching rows')
    IF (.NOT. (nz > 10 .AND. nz == nb .AND. nz == np .AND. nz == n0)) RETURN
    WRITE (*, '(A,ES11.4)') 'E1: max |EI|-zeta|f - EI|Bb| over the channels ', MAXVAL(ABS(tz(:, :nz) - tb(:, :nb)))
    CALL require(MAXVAL(ABS(tz(:, :nz) - tb(:, :nb))) <= 1.0e-6_wp*MAXVAL(ABS(tz(:, :nz))), &
                 'E1: EI|-zeta|f is EI|Bb with Bb = zeta EI/(pi f)')
    CALL require(MAXVAL(ABS(t0(:, :n0) - tp(:, :np))) == 0.0_wp, 'E2: EI|0 is the undamped deck exactly')
    rz = MAXVAL(tz(4, nz/2:nz)) - MINVAL(tz(4, nz/2:nz))
    rp = MAXVAL(tp(4, np/2:np)) - MINVAL(tp(4, np/2:np))
    r0 = MAXVAL(tz(3, nz/2:nz)) - MINVAL(tz(3, nz/2:nz))
    WRITE (*, '(A,ES11.4,A,ES11.4,A)') 'E2: hang-off curvature range, damped ', rz, ' 1/m, undamped ', rp, ' 1/m'
    CALL require(rz < rp, 'E2: the bending damping lowers the hang-off curvature range')
    CALL require(r0 > 0.0_wp, 'E2: the hang-off moment varies in the sea')
  END SUBROUTINE check_runs

  SUBROUTINE check_refusals()
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CALL write_deck('bdd_rope.dat', '0.0|250.0', 0.0_wp, rope=.TRUE.)
    CALL run('bdd_rope.dat', 'bdd_rope', ok, msg)
    WRITE (*, '(A)') 'E3 EI = 0: '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'bending damping (EI|Bb) needs a finite-EI cable') > 0, &
                 'E3: bending damping on an EI = 0 type')
    CALL write_deck('bdd_rope_ratio.dat', '0.0|-0.05|0.5', 0.0_wp, rope=.TRUE.)
    CALL run('bdd_rope_ratio.dat', 'bdd_rope_ratio', ok, msg)
    WRITE (*, '(A)') 'E3 EI = 0, ratio form: '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'bending damping (EI|Bb) needs a finite-EI cable') > 0, &
                 'E3: a damping ratio on an EI = 0 type')
    CALL write_deck('bdd_neg.dat', '2.0e4|-250.0', 0.0_wp)
    CALL run('bdd_neg.dat', 'bdd_neg', ok, msg)
    WRITE (*, '(A)') 'E3 negative without f: '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'needs the frequency it holds at') > 0, 'E3: negative Bb without f')
    CALL write_deck('bdd_five.dat', '2.0e4|-0.05|1|2|3', 0.0_wp)
    CALL run('bdd_five.dat', 'bdd_five', ok, msg)
    WRITE (*, '(A)') 'E3 five parts: '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'takes EI, EI|Bb, EI|-zeta|f or EI|-zeta|f1|f2') > 0, &
                 'E3: five-part EI column')
    CALL write_deck('bdd_pair_pos.dat', '2.0e4|0.05|0.5|2.0', 0.0_wp)
    CALL run('bdd_pair_pos.dat', 'bdd_pair_pos', ok, msg)
    WRITE (*, '(A)') 'E5 positive ratio: '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'EI|-zeta|f1|f2 takes a damping ratio 0 < zeta < 1') > 0, &
                 'E5: a positive ratio in the two-frequency form')
    CALL write_deck('bdd_pair_f0.dat', '2.0e4|-0.05|0.5|0', 0.0_wp)
    CALL run('bdd_pair_f0.dat', 'bdd_pair_f0', ok, msg)
    WRITE (*, '(A)') 'E5 zero frequency: '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'frequencies f1 and f2 of EI|-zeta|f1|f2 must be positive') > 0, &
                 'E5: a zero frequency in the two-frequency form')
    CALL write_deck('bdd_rope_pair.dat', '0.0|-0.05|0.5|2.0', 0.0_wp, rope=.TRUE.)
    CALL run('bdd_rope_pair.dat', 'bdd_rope_pair', ok, msg)
    WRITE (*, '(A)') 'E5 EI = 0, two-frequency form: '//TRIM(msg)
    CALL require(.NOT. ok .AND. INDEX(msg, 'bending damping (EI|Bb) needs a finite-EI cable') > 0, &
                 'E5: the two-frequency form on an EI = 0 type')
  END SUBROUTINE check_refusals

  SUBROUTINE check_pair()
    REAL(wp), PARAMETER :: F1 = 0.5_wp, F2 = 2.0_wp
    REAL(wp), ALLOCATABLE :: tr(:, :), tb(:, :), tp(:, :)
    INTEGER :: nr, nb, np
    REAL(wp) :: bb, rr, rb, rp
    LOGICAL :: ok
    CHARACTER(1024) :: msg
    CHARACTER(64) :: col
    bb = ZETA*EI/(PI*(F1 + F2))
    WRITE (col, '(ES22.15)') bb
    CALL write_deck('bdd_pair.dat', '2.0e4|-0.05|0.5|2.0', 20.0_wp)
    CALL run('bdd_pair.dat', 'bdd_pair', ok, msg)
    CALL require(ok, 'E4: the EI|-zeta|f1|f2 deck runs: '//TRIM(msg))
    CALL write_deck('bdd_pair_bb.dat', '2.0e4|'//TRIM(ADJUSTL(col)), 20.0_wp)
    CALL run('bdd_pair_bb.dat', 'bdd_pair_bb', ok, msg)
    CALL require(ok, 'E4: the EI|Bb deck of its Bb runs: '//TRIM(msg))
    CALL read_out('bdd_pair.out', tr, nr)
    CALL read_out('bdd_pair_bb.out', tb, nb)
    CALL read_out('bdd_plain.out', tp, np)
    CALL require(nr > 10 .AND. nr == nb .AND. nr == np, 'E4: the runs wrote matching rows')
    IF (.NOT. (nr > 10 .AND. nr == nb .AND. nr == np)) RETURN
    rr = MAXVAL(tr(4, nr/2:nr)) - MINVAL(tr(4, nr/2:nr))
    rb = MAXVAL(tb(4, nb/2:nb)) - MINVAL(tb(4, nb/2:nb))
    rp = MAXVAL(tp(4, np/2:np)) - MINVAL(tp(4, np/2:np))
    WRITE (*, '(A,ES11.4,A)') 'E4: max |EI|-zeta|f1|f2 - EI|Bb| over the channels ', &
      MAXVAL(ABS(tr(:, :nr) - tb(:, :nb))), ' (the alpha share)'
    WRITE (*, '(A,3ES11.4,A)') 'E4: hang-off curvature range, Rayleigh pair, its Bb alone, undamped ', rr, rb, rp, &
      ' 1/m'
    CALL require(MAXVAL(ABS(tr(:, :nr) - tb(:, :nb))) > 1.0e-6_wp*MAXVAL(ABS(tb(:, :nb))), &
                 'E4: the mass-proportional half acts')
    CALL require(rr < rb .AND. rb < rp, 'E4: the mass-proportional half lowers the curvature range further')
  END SUBROUTINE check_pair

END PROGRAM test_bending_damping_deck
