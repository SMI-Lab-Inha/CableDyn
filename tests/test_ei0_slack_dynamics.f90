! File: tests/test_ei0_slack_dynamics.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_ei0_slack_dynamics
  !! Tension-only (EI = 0) lines that go slack and snap taut, at the default rhoInf = 0.4:
  !!   1. A light, axially stiff cable (80 kg/m, EA 8e8 N, 115 x 2 m segments) under 3 m / 10 s
  !!      heave at dtM 0.1 s. The configuration-blended step reported 0.61 MN at the fairlead and
  !!      20 MN at the anchor against about 0.1 MN (rhoInf 0 and small-step references); the
  !!      force blend stays below 0.2 MN.
  !!   2. The same cable at dtM 0.01 s, where the force blend alone still fed a rattle of
  !!      single-step segment spikes (2.3 MN by 24 s) through the energy a release creates;
  !!      with the release cap every segment stays below 0.2 MN.
  !!   3. The shipped viscoelastic polyester (Es|Ed, Bs|Bd) under 15 m surge at dtM 0.1 s:
  !!      1.80 MN at the fairlead (2.69 MN configuration-blended), as at dtM 0.0025 s.
  !!   4. The run audit of the driver (stdout of the executable given as the first argument):
  !!      the configuration blend (False alpha_force_blend) of case 1 warns of a non-physical
  !!      tension spike, the force blend does not; rhoInf above 0.95 prints its note; a 10 m
  !!      sway at dtM 0.2 s prints the segment-turn advisory; the unresolved slack-to-taut
  !!      snaps of case 1 print the snap note.
  !!   5. The Newton tangent of the force-blended step against a central finite difference of its
  !!      residual (the step's probe_q hook), on a three-segment tension-only chain with one fixed
  !!      node: all segments taut at t_n and at the trial state (release cap inactive), a segment
  !!      taut at t_n and deeply slack at the trial state (cap active, its energy-gradient force
  !!      W_n/(e_n - e_{n+1}) and tangent), and a segment slack at both states; each on the
  !!      banded and the dense path. The cap's residual is also checked against the closed form.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  USE CableDyn_Dynamic, ONLY: GenAlphaConfig, CD_CableGenAlphaWorkspace, CD_Cable_Gen_Alpha_Step, CD_DYN_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.141592653589793238462643383279502884197_wp
  CHARACTER(*), PARAMETER :: CABLE = 'cable 0.2 80.0 8.0e8 0.0 0.0 1.2 0.1 1.0 0.0'
  CHARACTER(*), PARAMETER :: POLY_VE = 'poly_ve 0.1438 22.42 1.424e8|2.50e8 4.0e9|1.1e7 0.0 1.2 0.2 1.0 0.0'
  INTEGER :: nfail
  CHARACTER(1024) :: driver

  nfail = 0
  driver = ''
  IF (COMMAND_ARGUMENT_COUNT() >= 1) CALL GET_COMMAND_ARGUMENT(1, driver)
  CALL case_tangent_consistency()
  CALL case_light_cable_coarse_step()
  CALL case_light_cable_release_cap()
  CALL case_viscoelastic_surge()
  IF (LEN_TRIM(driver) > 0) THEN
    CALL case_run_audit()
  ELSE
    WRITE (*, '(A)') 'MISMATCH [audit: the driver executable is the first argument]'
    nfail = nfail + 1
  END IF

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: EI = 0 slack-taut dynamics (force blend, release cap, run audit)'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. cond) THEN
      WRITE (*, '(A)') 'MISMATCH ['//label//']'
      nfail = nfail + 1
    END IF
  END SUBROUTINE require

  SUBROUTINE write_motion(path, pid, x0, axis, amp, period, ramp, dt, tmax)
    !! Harmonic translation of one point, amp*r(t)*sin(w t), with the C2 quintic start-up
    !! ramp r over `ramp` seconds and its exact velocity and acceleration.
    CHARACTER(*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: pid, axis
    REAL(wp), INTENT(IN) :: x0(3), amp, period, ramp, dt, tmax
    INTEGER :: u, i, n
    REAL(wp) :: t, s, r, rd, rdd, w, p(3), v(3), a(3)
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') '# time id x y z vx vy vz ax ay az'
    n = NINT(tmax/dt)
    w = 2.0_wp*PI/period
    DO i = 0, n
      t = REAL(i, wp)*dt
      s = MIN(t/ramp, 1.0_wp)
      r = 10.0_wp*s**3 - 15.0_wp*s**4 + 6.0_wp*s**5
      rd = 0.0_wp
      rdd = 0.0_wp
      IF (t < ramp) THEN
        rd = (30.0_wp*s**2 - 60.0_wp*s**3 + 30.0_wp*s**4)/ramp
        rdd = (60.0_wp*s - 180.0_wp*s**2 + 120.0_wp*s**3)/ramp**2
      END IF
      p = x0
      v = 0.0_wp
      a = 0.0_wp
      p(axis) = p(axis) + amp*r*SIN(w*t)
      v(axis) = amp*(rd*SIN(w*t) + r*w*COS(w*t))
      a(axis) = amp*(rdd*SIN(w*t) + 2.0_wp*rd*w*COS(w*t) - r*w*w*SIN(w*t))
      WRITE (u, '(ES17.10,1X,I0,9(1X,ES17.10))') t, pid, p, v, a
    END DO
    CLOSE (u)
  END SUBROUTINE write_motion

  SUBROUTINE write_deck(path, linetype, anchor, fairlead, length, nseg, depth, outputs, options)
    !! One line from a Coupled fairlead (point 1, End A) to a Fixed anchor (point 2).
    CHARACTER(*), INTENT(IN) :: path, linetype, outputs, options(:)
    REAL(wp), INTENT(IN) :: anchor(3), fairlead(3), length, depth
    INTEGER, INTENT(IN) :: nseg
    INTEGER :: u, i
    CHARACTER(32) :: name
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'EI = 0 slack-taut dynamics regression deck'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') TRIM(linetype)
    READ (linetype, *) name
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)'
    WRITE (u, '(A,3(1X,ES16.9),A)') '1 Coupled', fairlead, ' 0 0 0 0'
    WRITE (u, '(A,3(1X,ES16.9),A)') '2 Fixed', anchor, ' 0 0 0 0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A,A)') '1 1 2 ', outputs
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A,A,1X,ES16.9,1X,I0)') '1 ', TRIM(name), length, nseg
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(ES16.9,A)') depth, ' WtrDpth'
    DO i = 1, SIZE(options)
      IF (LEN_TRIM(options(i)) > 0) WRITE (u, '(A)') TRIM(options(i))
    END DO
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1 AnchTen1'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_deck

  SUBROUTINE light_cable_deck(path, dt, tmax, extra, outputs)
    !! The light, axially stiff cable from the fairlead (5, 0, -14) to the anchor (205, 0, -80),
    !! 230 m in 115 segments, heaving 3 m at 10 s after a 10 s ramp.
    CHARACTER(*), INTENT(IN) :: path, extra, outputs
    REAL(wp), INTENT(IN) :: dt, tmax
    CHARACTER(64) :: opts(4)
    CALL write_motion('ei0slack_heave.txt', 1, [5.0_wp, 0.0_wp, -14.0_wp], 3, 3.0_wp, 10.0_wp, 10.0_wp, dt, tmax)
    WRITE (opts(1), '(ES14.7,A)') dt, ' dtM'
    WRITE (opts(2), '(ES14.7,A)') tmax, ' TMax'
    opts(3) = 'ei0slack_heave.txt motionFile'
    opts(4) = extra
    CALL write_deck(path, CABLE, [205.0_wp, 0.0_wp, -80.0_wp], [5.0_wp, 0.0_wp, -14.0_wp], 230.0_wp, 115, 80.0_wp, &
                    outputs, opts)
  END SUBROUTINE light_cable_deck

  SUBROUTINE read_out(path, data, nrow)
    !! Numeric rows of a driver output file (comment and header lines skipped).
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), ALLOCATABLE, INTENT(OUT) :: data(:, :)
    INTEGER, INTENT(OUT) :: nrow
    INTEGER :: u, ios, ncol, k
    CHARACTER(8192) :: line
    LOGICAL :: header_seen
    REAL(wp), ALLOCATABLE :: grow(:, :)
    nrow = 0
    ncol = 0
    header_seen = .FALSE.
    ALLOCATE (data(0, 0))
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    DO
      READ (u, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF (line(1:1) == '#' .OR. LEN_TRIM(line) == 0) CYCLE
      IF (.NOT. header_seen) THEN
        header_seen = .TRUE.
        ncol = count_tokens(line)
        DEALLOCATE (data)
        ALLOCATE (data(ncol, 64))
        CYCLE
      END IF
      IF (nrow == SIZE(data, 2)) THEN
        ALLOCATE (grow(ncol, 2*SIZE(data, 2)))
        grow(:, 1:nrow) = data(:, 1:nrow)
        CALL MOVE_ALLOC(grow, data)
      END IF
      nrow = nrow + 1
      READ (line, *) (data(k, nrow), k=1, ncol)
    END DO
    CLOSE (u)
  END SUBROUTINE read_out

  INTEGER FUNCTION count_tokens(line) RESULT(n)
    CHARACTER(*), INTENT(IN) :: line
    INTEGER :: i
    LOGICAL :: inside
    n = 0
    inside = .FALSE.
    DO i = 1, LEN_TRIM(line)
      IF (line(i:i) == ' ' .OR. line(i:i) == CHAR(9)) THEN
        inside = .FALSE.
      ELSE IF (.NOT. inside) THEN
        inside = .TRUE.
        n = n + 1
      END IF
    END DO
  END FUNCTION count_tokens

  REAL(wp) FUNCTION column_max(path, first_col, last_col) RESULT(mx)
    !! Largest value of the columns first_col..last_col (last_col 0: through the last column).
    CHARACTER(*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: first_col, last_col
    REAL(wp), ALLOCATABLE :: d(:, :)
    INTEGER :: nrow, lc
    mx = -HUGE(1.0_wp)
    CALL read_out(path, d, nrow)
    IF (nrow < 1) RETURN
    lc = last_col
    IF (lc == 0) lc = SIZE(d, 1)
    mx = MAXVAL(d(first_col:lc, 1:nrow))
  END FUNCTION column_max

  LOGICAL FUNCTION file_contains(path, needle) RESULT(found)
    CHARACTER(*), INTENT(IN) :: path, needle
    INTEGER :: u, ios
    CHARACTER(4096) :: line
    found = .FALSE.
    OPEN (NEWUNIT=u, FILE=path, STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) RETURN
    DO
      READ (u, '(A)', IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF (INDEX(line, needle) > 0) THEN
        found = .TRUE.
        EXIT
      END IF
    END DO
    CLOSE (u)
  END FUNCTION file_contains

  SUBROUTINE run_driver_logged(deck, root, log, ok)
    !! Run the driver executable on a deck, its stdout and stderr to log.
    CHARACTER(*), INTENT(IN) :: deck, root, log
    LOGICAL, INTENT(OUT) :: ok
    INTEGER :: ex, cs
    CALL EXECUTE_COMMAND_LINE('"'//TRIM(driver)//'" '//deck//' '//root//' > '//log//' 2>&1', &
                              WAIT=.TRUE., EXITSTAT=ex, CMDSTAT=cs)
    ok = cs == 0 .AND. ex == 0
  END SUBROUTINE run_driver_logged

  SUBROUTINE zero_load(q, v, force, jac_q, jac_v, ErrStat, ErrMsg)
    !! No dynamic load; selects the dense (load_proc) path of the step.
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    force = 0.0_wp
    jac_q = 0.0_wp
    jac_v = 0.0_wp
    ErrStat = 0
    ErrMsg = ''
  END SUBROUTINE zero_load

  SUBROUTINE chain_nodes(dirs, strain, l0, q)
    !! Nodes of the chain from the fixed node at the origin: segment e runs along dirs(:, e)
    !! with the length l0(e)*(1 + strain(e)).
    REAL(wp), INTENT(IN) :: dirs(3, 3), strain(3), l0(3)
    REAL(wp), INTENT(OUT) :: q(12)
    INTEGER :: e
    q(1:3) = 0.0_wp
    DO e = 1, 3
      q(3*e + 1:3*e + 3) = q(3*e - 2:3*e) + dirs(:, e)/NORM2(dirs(:, e))*l0(e)*(1.0_wp + strain(e))
    END DO
  END SUBROUTINE chain_nodes

  SUBROUTINE blended_force(qn, q1, l0, ea, alpha_f, capped, f)
    !! Closed form of the blended internal force (1 - alpha_f) f(q1) + alpha_f f(qn) of the
    !! tension-only chain, with the release cap (capped) or without it.
    REAL(wp), INTENT(IN) :: qn(12), q1(12), l0(3), ea(3), alpha_f
    LOGICAL, INTENT(IN) :: capped
    REAL(wp), INTENT(OUT) :: f(12)
    INTEGER :: e
    REAL(wp) :: dn(3), d1(3), ln, l1, t0, t1, e0, e1, dg, share(3)
    f = 0.0_wp
    DO e = 1, 3
      dn = qn(3*e + 1:3*e + 3) - qn(3*e - 2:3*e)
      d1 = q1(3*e + 1:3*e + 3) - q1(3*e - 2:3*e)
      ln = NORM2(dn)
      l1 = NORM2(d1)
      t0 = MAX(0.0_wp, ea(e)*(ln/l0(e) - 1.0_wp))
      t1 = MAX(0.0_wp, ea(e)*(l1/l0(e) - 1.0_wp))
      e0 = ln - l0(e)
      e1 = l1 - l0(e)
      dg = alpha_f*t0
      IF (capped .AND. t0 > 0.0_wp .AND. e1 < 0.0_wp) dg = MIN(dg, 0.5_wp*t0*e0/(e0 - e1))
      share = (1.0_wp - alpha_f)*t1*d1/l1 + dg*dn/ln
      f(3*e - 2:3*e) = f(3*e - 2:3*e) - share
      f(3*e + 1:3*e + 3) = f(3*e + 1:3*e + 3) + share
    END DO
  END SUBROUTINE blended_force

  SUBROUTINE probe_step(qn, vn, an, q1, l0, ea, rho_a, dense, R, K)
    !! Residual and effective tangent of the force-blended step at q1 (the probe_q hook).
    REAL(wp), INTENT(IN) :: qn(12), vn(12), an(12), q1(12), l0(3), ea(3), rho_a(3)
    LOGICAL, INTENT(IN) :: dense
    REAL(wp), INTENT(OUT) :: R(12), K(9, 9)
    INTEGER :: conn(2, 3), n_iter, es
    REAL(wp) :: qo(12), vo(12), ao(12)
    LOGICAL :: converged, stalled
    TYPE(GenAlphaConfig) :: cfg
    TYPE(CD_CableGenAlphaWorkspace) :: work
    CHARACTER(120) :: em
    conn = RESHAPE([1, 2, 2, 3, 3, 4], [2, 3])
    cfg%rho_inf = 0.4_wp
    IF (dense) THEN
      CALL CD_Cable_Gen_Alpha_Step(qn, vn, an, conn, l0, ea, rho_a, .TRUE., 0.0_wp*qn, [1, 2, 3], 0.05_wp, cfg, &
                                   qo, vo, ao, converged, stalled, n_iter, es, em, load_proc=zero_load, &
                                   workspace=work, probe_q=q1, probe_R=R, probe_K=K)
    ELSE
      CALL CD_Cable_Gen_Alpha_Step(qn, vn, an, conn, l0, ea, rho_a, .TRUE., 0.0_wp*qn, [1, 2, 3], 0.05_wp, cfg, &
                                   qo, vo, ao, converged, stalled, n_iter, es, em, &
                                   workspace=work, probe_q=q1, probe_R=R, probe_K=K)
    END IF
    IF (es /= CD_DYN_OK) THEN
      WRITE (*, '(A)') 'MISMATCH [tangent:probe '//TRIM(em)//']'
      nfail = nfail + 1
      R = 0.0_wp
      K = 0.0_wp
    END IF
  END SUBROUTINE probe_step

  SUBROUTINE tangent_config(label, strain_n, strain_1, expect_cap)
    !! One configuration of the three-segment chain (node 1 fixed): segment strains at t_n and
    !! at the trial state of the step, t_n at rest but for a small velocity and acceleration.
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: strain_n(3), strain_1(3)
    LOGICAL, INTENT(IN) :: expect_cap
    REAL(wp), PARAMETER :: L0(3) = [1.0_wp, 1.2_wp, 0.8_wp], EA(3) = [1.0e6_wp, 1.5e6_wp, 0.8e6_wp]
    REAL(wp), PARAMETER :: DIRS_N(3, 3) = RESHAPE([1.0_wp, 0.2_wp, 0.1_wp, 0.3_wp, 1.0_wp, -0.2_wp, &
                                                   -0.2_wp, 0.4_wp, 1.0_wp], [3, 3])
    REAL(wp), PARAMETER :: DIRS_1(3, 3) = RESHAPE([1.0_wp, 0.27_wp, 0.04_wp, 0.36_wp, 1.0_wp, -0.12_wp, &
                                                   -0.12_wp, 0.46_wp, 1.0_wp], [3, 3])
    REAL(wp), PARAMETER :: H = 1.0e-6_wp, ALPHA_F = 0.4_wp/1.4_wp
    INTEGER :: p, j, free(9), i
    REAL(wp) :: qn(12), vn(12), an(12), q1(12), qp(12), rho_a(3), R(12), Rp(12), Rm(12), K(9, 9), Kfd(9, 9)
    REAL(wp) :: Kdump(9, 9), f_cap(12), f_raw(12), scale, err
    LOGICAL :: dense
    free = [(i + 3, i=1, 9)]
    CALL chain_nodes(DIRS_N, strain_n, L0, qn)
    CALL chain_nodes(DIRS_1, strain_1, L0, q1)
    vn = 0.0_wp
    an = 0.0_wp
    DO i = 4, 12
      vn(i) = 0.01_wp*SIN(1.7_wp*REAL(i, wp))
      an(i) = 0.02_wp*COS(2.3_wp*REAL(i, wp))
    END DO
    ! the closed-form blended force: the residual of a nearly massless chain is the force
    rho_a = 1.0e-9_wp
    CALL probe_step(qn, 0.0_wp*vn, 0.0_wp*an, q1, L0, EA, rho_a, .FALSE., R, Kdump)
    CALL blended_force(qn, q1, L0, EA, ALPHA_F, .TRUE., f_cap)
    CALL blended_force(qn, q1, L0, EA, ALPHA_F, .FALSE., f_raw)
    scale = MAXVAL(ABS(f_raw))
    CALL require(MAXVAL(ABS(R(4:12) - f_cap(4:12))) < 1.0e-6_wp*scale, 'tangent:'//label//':residual-closed-form')
    IF (expect_cap) THEN
      CALL require(MAXVAL(ABS(f_cap - f_raw)) > 1.0e-2_wp*scale, 'tangent:'//label//':cap-active')
    ELSE
      CALL require(MAXVAL(ABS(f_cap - f_raw)) < 1.0e-12_wp*scale, 'tangent:'//label//':cap-inactive')
    END IF
    ! the tangent against the central difference of the residual, with inertia
    rho_a = [2.0_wp, 3.0_wp, 1.5_wp]
    DO p = 1, 2
      dense = p == 2
      CALL probe_step(qn, vn, an, q1, L0, EA, rho_a, dense, R, K)
      DO j = 1, 9
        qp = q1
        qp(free(j)) = q1(free(j)) + H
        CALL probe_step(qn, vn, an, qp, L0, EA, rho_a, dense, Rp, Kdump)
        qp(free(j)) = q1(free(j)) - H
        CALL probe_step(qn, vn, an, qp, L0, EA, rho_a, dense, Rm, Kdump)
        Kfd(:, j) = (Rp(4:12) - Rm(4:12))/(2.0_wp*H)
      END DO
      scale = MAXVAL(ABS(K))
      err = MAXVAL(ABS(K - Kfd))/scale
      WRITE (*, '(A,A,A,L1,A,ES9.2)') '  tangent ', label, ' dense=', dense, ': max |K - K_fd|/max|K| = ', err
      CALL require(err < 1.0e-6_wp, 'tangent:'//label//TRIM(MERGE(':dense ', ':banded', dense))//':fd-tangent')
    END DO
  END SUBROUTINE tangent_config

  SUBROUTINE case_tangent_consistency()
    ! (a) all segments taut at t_n and at the trial state: the cap is inactive
    CALL tangent_config('taut-taut', [0.010_wp, 0.020_wp, 0.010_wp], [0.020_wp, 0.015_wp, 0.020_wp], .FALSE.)
    ! (b) segment 1 taut at t_n and deeply slack at the trial state (release deeper than
    !     0.75 e_n): it carries the energy gradient W_n/(e_n - e_{n+1}) instead of alpha_f T_n
    CALL tangent_config('release-cap', [0.010_wp, 0.020_wp, 0.010_wp], [-0.050_wp, 0.015_wp, 0.020_wp], .TRUE.)
    ! (c) segment 2 slack at both states
    CALL tangent_config('slack-slack', [0.010_wp, -0.020_wp, 0.010_wp], [0.020_wp, -0.030_wp, 0.015_wp], .FALSE.)
  END SUBROUTINE case_tangent_consistency

  SUBROUTINE case_light_cable_coarse_step()
    LOGICAL :: conv
    INTEGER :: es
    CHARACTER(512) :: em
    REAL(wp) :: fair, anch
    CALL light_cable_deck('ei0slack_coarse.dat', 0.1_wp, 40.0_wp, '', '-')
    CALL CD_Run_Deck_Driver('ei0slack_coarse.dat', 'ei0slack_coarse', conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, 'coarse:run '//TRIM(em))
    fair = column_max('ei0slack_coarse.out', 2, 2)
    anch = column_max('ei0slack_coarse.out', 3, 3)
    WRITE (*, '(A,ES10.3,A,ES10.3,A)') '  light cable, dtM 0.1 s: FairTen1 max ', fair, ' N, AnchTen1 max ', anch, ' N'
    CALL require(fair > 5.0e4_wp .AND. fair < 2.0e5_wp, 'coarse:fairlead-peak')
    CALL require(anch > 2.0e4_wp .AND. anch < 2.0e5_wp, 'coarse:anchor-peak')
  END SUBROUTINE case_light_cable_coarse_step

  SUBROUTINE case_light_cable_release_cap()
    LOGICAL :: conv
    INTEGER :: es
    CHARACTER(512) :: em
    REAL(wp) :: seg
    CALL light_cable_deck('ei0slack_release.dat', 0.01_wp, 24.0_wp, '', 't')
    CALL CD_Run_Deck_Driver('ei0slack_release.dat', 'ei0slack_release', conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, 'release:run '//TRIM(em))
    seg = column_max('ei0slack_release.Line1.t.out', 2, 0)
    WRITE (*, '(A,ES10.3,A)') '  light cable, dtM 0.01 s: largest segment tension ', seg, ' N'
    CALL require(seg > 5.0e4_wp .AND. seg < 2.0e5_wp, 'release:segment-peak')
  END SUBROUTINE case_light_cable_release_cap

  SUBROUTINE case_viscoelastic_surge()
    LOGICAL :: conv
    INTEGER :: es
    CHARACTER(512) :: em
    CHARACTER(64) :: opts(4)
    REAL(wp) :: fair
    CALL write_motion('ei0slack_surge.txt', 1, [58.0_wp, 0.0_wp, -14.0_wp], 1, 15.0_wp, 12.0_wp, 12.0_wp, &
                      0.1_wp, 48.0_wp)
    opts(1) = '0.1 dtM'
    opts(2) = '48.0 TMax'
    opts(3) = 'ei0slack_surge.txt motionFile'
    opts(4) = '3.0e6 kBot'
    CALL write_deck('ei0slack_ve.dat', POLY_VE, [700.0_wp, 0.0_wp, -200.0_wp], [58.0_wp, 0.0_wp, -14.0_wp], &
                    690.0_wp, 50, 200.0_wp, '-', opts)
    CALL CD_Run_Deck_Driver('ei0slack_ve.dat', 'ei0slack_ve', conv, es, em)
    CALL require(es == CD_DECKDRV_OK .AND. conv, 've:run '//TRIM(em))
    fair = column_max('ei0slack_ve.out', 2, 2)
    WRITE (*, '(A,ES10.3,A)') '  viscoelastic polyester, 15 m surge, dtM 0.1 s: FairTen1 max ', fair, ' N'
    CALL require(fair > 1.6e6_wp .AND. fair < 1.9e6_wp, 've:fairlead-peak')
  END SUBROUTINE case_viscoelastic_surge

  SUBROUTINE case_run_audit()
    LOGICAL :: ok
    CHARACTER(64) :: opts(4)
    ! the configuration blend of case 1: its spikes are reported once, as a warning
    CALL light_cable_deck('ei0slack_audit_cfg.dat', 0.1_wp, 40.0_wp, 'False alpha_force_blend', '-')
    CALL run_driver_logged('ei0slack_audit_cfg.dat', 'ei0slack_audit_cfg', 'ei0slack_audit_cfg.log', ok)
    CALL require(ok, 'audit:configuration-blend-run')
    CALL require(file_contains('ei0slack_audit_cfg.log', 'Warning: line 1: segment '), &
                 'audit:configuration-blend-spike-warned')
    CALL require(file_contains('ei0slack_audit_cfg.log', 'for a single step'), 'audit:spike-message')
    ! the snaps of the light cable (axial wave crossing a segment ~160 times per step) are noted
    CALL require(file_contains('ei0slack_audit_cfg.log', 'went from slack to taut above the initial tension'), &
                 'audit:snap-note')
    CALL require(file_contains('ei0slack_audit_cfg.log', 'snap-load peaks are not resolved in time'), &
                 'audit:snap-note-text')
    ! the default force blend: no spike warning
    CALL light_cable_deck('ei0slack_audit_fb.dat', 0.1_wp, 40.0_wp, '', '-')
    CALL run_driver_logged('ei0slack_audit_fb.dat', 'ei0slack_audit_fb', 'ei0slack_audit_fb.log', ok)
    CALL require(ok, 'audit:force-blend-run')
    CALL require(.NOT. file_contains('ei0slack_audit_fb.log', 'Warning: line'), 'audit:force-blend-no-warning')
    CALL require(.NOT. file_contains('ei0slack_audit_fb.log', 'Note: rhoInf'), 'audit:default-rhoinf-no-note')
    ! rhoInf above 0.95
    CALL light_cable_deck('ei0slack_audit_rho.dat', 0.1_wp, 12.0_wp, '0.99 rhoInf', '-')
    CALL run_driver_logged('ei0slack_audit_rho.dat', 'ei0slack_audit_rho', 'ei0slack_audit_rho.log', ok)
    CALL require(ok, 'audit:rhoinf-run')
    CALL require(file_contains('ei0slack_audit_rho.log', 'Note: rhoInf = 0.990 leaves almost no numerical'), &
                 'audit:rhoinf-note')
    ! a 10 m sway at dtM 0.2 s turns the segments by several degrees per step
    CALL write_motion('ei0slack_sway.txt', 1, [5.0_wp, 0.0_wp, -14.0_wp], 2, 10.0_wp, 10.0_wp, 10.0_wp, &
                      0.2_wp, 20.0_wp)
    opts(1) = '0.2 dtM'
    opts(2) = '20.0 TMax'
    opts(3) = 'ei0slack_sway.txt motionFile'
    opts(4) = ''
    CALL write_deck('ei0slack_audit_sway.dat', CABLE, [205.0_wp, 0.0_wp, -80.0_wp], [5.0_wp, 0.0_wp, -14.0_wp], &
                    230.0_wp, 115, 80.0_wp, '-', opts)
    CALL run_driver_logged('ei0slack_audit_sway.dat', 'ei0slack_audit_sway', 'ei0slack_audit_sway.log', ok)
    CALL require(ok, 'audit:sway-run')
    CALL require(file_contains('ei0slack_audit_sway.log', 'Note: line 1: a taut segment turned up to'), &
                 'audit:segment-turn-note')
    CALL require(.NOT. file_contains('ei0slack_audit_sway.log', 'Warning: line'), 'audit:sway-no-warning')
  END SUBROUTINE case_run_audit

END PROGRAM test_ei0_slack_dynamics
