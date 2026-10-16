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
  !!   6. The same chain with a segment whose tension-only law a load callback carries (ea = 0,
  !!      as the viscoelastic and Syrope segments of a line model): with its step-start tension
  !!      (start_tension) its deep release is capped like an elastic one, closed form and finite
  !!      difference tangent; without it the segment keeps the plain force blend.
  !!   7. The finite-difference tangent of the step with nonzero dynamic loads (a quadratic drag
  !!      and a cubic position load, plus a position-dependent added mass on the dense path), and
  !!      the closed form of the blended load (1 - alpha_f) f(q_{n+1}, v_{n+1}) + alpha_f f(q_n, v_n).
  !!   8. A warm re-solve of the same step (guess_armed, the junction coupling) reuses the
  !!      step-start force, load and release-cap data of the previous call: its solution equals a
  !!      cold solve bit for bit while the start data it is given are wrong (so they are not
  !!      read); after an unarmed call that fails before computing them, an armed call recomputes
  !!      them instead of reusing those of an older state.
  !!   9. A lightly damped viscoelastic version of the light cable at dtM 0.1 s (the cap on its
  !!      law segments) and a long taut Syrope line in 4 m surge that goes slack: no spike
  !!      warning, and peaks within the band of small-step references.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_DeckDriver, ONLY: CD_Run_Deck_Driver, CD_DECKDRV_OK
  USE CableDyn_Dynamic, ONLY: GenAlphaConfig, CD_CableGenAlphaWorkspace, CD_Cable_Gen_Alpha_Step, CD_DYN_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.141592653589793238462643383279502884197_wp
  CHARACTER(*), PARAMETER :: CABLE = 'cable 0.2 80.0 8.0e8 0.0 0.0 1.2 0.1 1.0 0.0'
  CHARACTER(*), PARAMETER :: POLY_VE = 'poly_ve 0.1438 22.42 1.424e8|2.50e8 4.0e9|1.1e7 0.0 1.2 0.2 1.0 0.0'
  ! the light cable as a lightly damped viscoelastic line (Es|Ed, Bs|Bd)
  CHARACTER(*), PARAMETER :: CABLE_VE = 'cable 0.2 80.0 8.0e8|1.6e9 1.0e8|1.0e6 0.0 1.2 0.1 1.0 0.0'
  INTEGER :: nfail
  CHARACTER(1024) :: driver
  ! dynamic load of the chain cases 6 to 8 (host-associated by the load callbacks): nodal
  ! quadratic drag -cb_drag |v| v, nodal cubic load -cb_spring |q - cb_q0|^2 (q - cb_q0), the
  ! tension-only law (EA cb_ea_law, length cb_l0_law) of segment cb_law_elem (0: none), and the
  ! added mass cb_madd (1 + 0.3 z^2) per DOF of the dense path
  REAL(wp) :: cb_drag = 0.0_wp, cb_spring = 0.0_wp, cb_ea_law = 0.0_wp, cb_l0_law = 1.0_wp, cb_madd = 0.0_wp
  REAL(wp) :: cb_q0(12) = 0.0_wp
  INTEGER :: cb_law_elem = 0

  nfail = 0
  driver = ''
  IF (COMMAND_ARGUMENT_COUNT() >= 1) CALL GET_COMMAND_ARGUMENT(1, driver)
  CALL case_tangent_consistency()
  CALL case_law_release_cap()
  CALL case_load_share_tangent()
  CALL case_warm_resolve_reuse()
  CALL case_light_cable_coarse_step()
  CALL case_light_cable_release_cap()
  CALL case_viscoelastic_surge()
  IF (LEN_TRIM(driver) > 0) THEN
    CALL case_run_audit()
    CALL case_law_slack_runs()
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

  ! ------------------------------------------------------------------------------------------
  ! Chain cases 6 to 8: a load callback with drag, a position load and a segment law
  ! ------------------------------------------------------------------------------------------

  FUNCTION outer3(x, y) RESULT(m)
    REAL(wp), INTENT(IN) :: x(3), y(3)
    REAL(wp) :: m(3, 3)
    INTEGER :: i
    DO i = 1, 3
      m(:, i) = x*y(i)
    END DO
  END FUNCTION outer3

  SUBROUTINE chain_load(q, v, force, jac_q, jac_v)
    !! The chain load (see cb_drag) and its residual Jacobians -d(force)/dq, -d(force)/dv.
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER :: n, i, a, b
    REAL(wp) :: d(3), w(3), dd, sp, c(3), ell, t(3), ten, kb(3, 3), eye(3, 3)
    force = 0.0_wp
    jac_q = 0.0_wp
    jac_v = 0.0_wp
    eye = 0.0_wp
    DO i = 1, 3
      eye(i, i) = 1.0_wp
    END DO
    DO n = 1, SIZE(q)/3
      d = q(3*n - 2:3*n) - cb_q0(3*n - 2:3*n)
      w = v(3*n - 2:3*n)
      dd = DOT_PRODUCT(d, d)
      sp = NORM2(w)
      force(3*n - 2:3*n) = -cb_drag*sp*w - cb_spring*dd*d
      jac_q(3*n - 2:3*n, 3*n - 2:3*n) = cb_spring*(dd*eye + 2.0_wp*outer3(d, d))
      IF (sp > 0.0_wp) jac_v(3*n - 2:3*n, 3*n - 2:3*n) = cb_drag*(sp*eye + outer3(w, w)/sp)
    END DO
    IF (cb_law_elem > 0) THEN
      ! the tension-only segment law: force T t on node a, -T t on node b
      a = cb_law_elem
      b = a + 1
      c = q(3*b - 2:3*b) - q(3*a - 2:3*a)
      ell = NORM2(c)
      t = c/ell
      ten = cb_ea_law*(ell/cb_l0_law - 1.0_wp)
      IF (ten > 0.0_wp) THEN
        force(3*a - 2:3*a) = force(3*a - 2:3*a) + ten*t
        force(3*b - 2:3*b) = force(3*b - 2:3*b) - ten*t
        kb = (cb_ea_law/cb_l0_law)*outer3(t, t) + ten*(eye - outer3(t, t))/ell
        jac_q(3*a - 2:3*a, 3*a - 2:3*a) = jac_q(3*a - 2:3*a, 3*a - 2:3*a) + kb
        jac_q(3*a - 2:3*a, 3*b - 2:3*b) = jac_q(3*a - 2:3*a, 3*b - 2:3*b) - kb
        jac_q(3*b - 2:3*b, 3*a - 2:3*a) = jac_q(3*b - 2:3*b, 3*a - 2:3*a) - kb
        jac_q(3*b - 2:3*b, 3*b - 2:3*b) = jac_q(3*b - 2:3*b, 3*b - 2:3*b) + kb
      END IF
    END IF
  END SUBROUTINE chain_load

  SUBROUTINE cb_load(q, v, force, jac_q, jac_v, ErrStat, ErrMsg)
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:), jac_q(:, :), jac_v(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    CALL chain_load(q, v, force, jac_q, jac_v)
    ErrStat = 0
    ErrMsg = ''
  END SUBROUTINE cb_load

  SUBROUTINE cb_force(q, v, force, ErrStat, ErrMsg)
    REAL(wp), INTENT(IN) :: q(:), v(:)
    REAL(wp), INTENT(OUT) :: force(:)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: jq(SIZE(q), SIZE(q)), jv(SIZE(q), SIZE(q))
    CALL chain_load(q, v, force, jq, jv)
    ErrStat = 0
    ErrMsg = ''
  END SUBROUTINE cb_force

  SUBROUTINE cb_load_band(q, v, free, kl, ku, force, jac_q_band, jac_v_band, ErrStat, ErrMsg)
    !! The chain load with its Jacobians in the reduced free/free band storage of the step.
    REAL(wp), INTENT(IN) :: q(:), v(:)
    INTEGER, INTENT(IN) :: free(:), kl, ku
    REAL(wp), INTENT(OUT) :: force(:), jac_q_band(:, :), jac_v_band(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    REAL(wp) :: jq(SIZE(q), SIZE(q)), jv(SIZE(q), SIZE(q))
    INTEGER :: i, j
    CALL chain_load(q, v, force, jq, jv)
    jac_q_band = 0.0_wp
    jac_v_band = 0.0_wp
    DO j = 1, SIZE(free)
      DO i = 1, SIZE(free)
        IF (i - j > kl .OR. j - i > ku) CYCLE
        jac_q_band(kl + ku + 1 + i - j, j) = jq(free(i), free(j))
        jac_v_band(kl + ku + 1 + i - j, j) = jv(free(i), free(j))
      END DO
    END DO
    ErrStat = 0
    ErrMsg = ''
  END SUBROUTINE cb_load_band

  SUBROUTINE cb_added_mass(q, accel, M_add, dMa_a_dq, ErrStat, ErrMsg)
    !! Position-dependent nodal added mass cb_madd (1 + 0.3 z^2) and d(M_add a)/dq.
    REAL(wp), INTENT(IN) :: q(:), accel(:)
    REAL(wp), INTENT(OUT) :: M_add(:, :), dMa_a_dq(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    INTEGER :: n, i
    M_add = 0.0_wp
    dMa_a_dq = 0.0_wp
    DO n = 1, SIZE(q)/3
      DO i = 3*n - 2, 3*n
        M_add(i, i) = cb_madd*(1.0_wp + 0.3_wp*q(3*n)**2)
        dMa_a_dq(i, 3*n) = cb_madd*0.6_wp*q(3*n)*accel(i)
      END DO
    END DO
    ErrStat = 0
    ErrMsg = ''
  END SUBROUTINE cb_added_mass

  SUBROUTINE probe_chain(qn, vn, an, q1, l0, ea, rho_a, banded, R, K, start_tension)
    !! Residual and effective tangent at q1 of the force-blended step of the chain with the
    !! callback load: banded (load_band_proc + load_force_proc) or dense (load_proc, with the
    !! added mass when cb_madd > 0).
    REAL(wp), INTENT(IN) :: qn(12), vn(12), an(12), q1(12), l0(3), ea(3), rho_a(3)
    LOGICAL, INTENT(IN) :: banded
    REAL(wp), INTENT(OUT) :: R(12), K(9, 9)
    REAL(wp), INTENT(IN), OPTIONAL :: start_tension(:)
    INTEGER :: conn(2, 3), n_iter, es
    REAL(wp) :: qo(12), vo(12), ao(12)
    LOGICAL :: converged, stalled
    TYPE(GenAlphaConfig) :: cfg
    TYPE(CD_CableGenAlphaWorkspace) :: work
    CHARACTER(120) :: em
    conn = RESHAPE([1, 2, 2, 3, 3, 4], [2, 3])
    cfg%rho_inf = 0.4_wp
    IF (banded) THEN
      CALL CD_Cable_Gen_Alpha_Step(qn, vn, an, conn, l0, ea, rho_a, .TRUE., 0.0_wp*qn, [1, 2, 3], 0.05_wp, cfg, &
                                   qo, vo, ao, converged, stalled, n_iter, es, em, load_band_proc=cb_load_band, &
                                   load_force_proc=cb_force, workspace=work, start_tension=start_tension, &
                                   probe_q=q1, probe_R=R, probe_K=K)
    ELSE IF (cb_madd > 0.0_wp) THEN
      CALL CD_Cable_Gen_Alpha_Step(qn, vn, an, conn, l0, ea, rho_a, .TRUE., 0.0_wp*qn, [1, 2, 3], 0.05_wp, cfg, &
                                   qo, vo, ao, converged, stalled, n_iter, es, em, load_proc=cb_load, &
                                   added_mass_proc=cb_added_mass, workspace=work, start_tension=start_tension, &
                                   probe_q=q1, probe_R=R, probe_K=K)
    ELSE
      CALL CD_Cable_Gen_Alpha_Step(qn, vn, an, conn, l0, ea, rho_a, .TRUE., 0.0_wp*qn, [1, 2, 3], 0.05_wp, cfg, &
                                   qo, vo, ao, converged, stalled, n_iter, es, em, load_proc=cb_load, &
                                   workspace=work, start_tension=start_tension, probe_q=q1, probe_R=R, probe_K=K)
    END IF
    IF (es /= CD_DYN_OK) THEN
      WRITE (*, '(A)') 'MISMATCH [chain:probe '//TRIM(em)//']'
      nfail = nfail + 1
      R = 0.0_wp
      K = 0.0_wp
    END IF
  END SUBROUTINE probe_chain

  SUBROUTINE fd_tangent_check(label, qn, vn, an, q1, l0, ea, rho_a, banded, start_tension)
    !! The step's tangent against the central difference of its residual.
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), INTENT(IN) :: qn(12), vn(12), an(12), q1(12), l0(3), ea(3), rho_a(3)
    LOGICAL, INTENT(IN) :: banded
    REAL(wp), INTENT(IN), OPTIONAL :: start_tension(:)
    REAL(wp), PARAMETER :: H = 1.0e-6_wp
    INTEGER :: j
    REAL(wp) :: qp(12), R(12), Rp(12), Rm(12), K(9, 9), Kfd(9, 9), Kdump(9, 9), err
    CALL probe_chain(qn, vn, an, q1, l0, ea, rho_a, banded, R, K, start_tension)
    DO j = 1, 9
      qp = q1
      qp(3 + j) = q1(3 + j) + H
      CALL probe_chain(qn, vn, an, qp, l0, ea, rho_a, banded, Rp, Kdump, start_tension)
      qp(3 + j) = q1(3 + j) - H
      CALL probe_chain(qn, vn, an, qp, l0, ea, rho_a, banded, Rm, Kdump, start_tension)
      Kfd(:, j) = (Rp(4:12) - Rm(4:12))/(2.0_wp*H)
    END DO
    err = MAXVAL(ABS(K - Kfd))/MAXVAL(ABS(K))
    WRITE (*, '(A,A,A,L1,A,ES9.2)') '  tangent ', label, ' banded=', banded, ': max |K - K_fd|/max|K| = ', err
    CALL require(err < 1.0e-6_wp, label//TRIM(MERGE(':banded', ':dense ', banded))//':fd-tangent')
  END SUBROUTINE fd_tangent_check

  SUBROUTINE chain_states(strain_n, strain_1, l0, qn, vn, an, q1)
    !! The chain of case 5 at t_n and at the trial state, with a small velocity and
    !! acceleration at t_n.
    REAL(wp), INTENT(IN) :: strain_n(3), strain_1(3), l0(3)
    REAL(wp), INTENT(OUT) :: qn(12), vn(12), an(12), q1(12)
    REAL(wp), PARAMETER :: DIRS_N(3, 3) = RESHAPE([1.0_wp, 0.2_wp, 0.1_wp, 0.3_wp, 1.0_wp, -0.2_wp, &
                                                   -0.2_wp, 0.4_wp, 1.0_wp], [3, 3])
    REAL(wp), PARAMETER :: DIRS_1(3, 3) = RESHAPE([1.0_wp, 0.27_wp, 0.04_wp, 0.36_wp, 1.0_wp, -0.12_wp, &
                                                   -0.12_wp, 0.46_wp, 1.0_wp], [3, 3])
    INTEGER :: i
    CALL chain_nodes(DIRS_N, strain_n, l0, qn)
    CALL chain_nodes(DIRS_1, strain_1, l0, q1)
    vn = 0.0_wp
    an = 0.0_wp
    DO i = 4, 12
      vn(i) = 0.01_wp*SIN(1.7_wp*REAL(i, wp))
      an(i) = 0.02_wp*COS(2.3_wp*REAL(i, wp))
    END DO
  END SUBROUTINE chain_states

  SUBROUTINE reset_chain_load()
    cb_drag = 0.0_wp
    cb_spring = 0.0_wp
    cb_ea_law = 0.0_wp
    cb_l0_law = 1.0_wp
    cb_madd = 0.0_wp
    cb_q0 = 0.0_wp
    cb_law_elem = 0
  END SUBROUTINE reset_chain_load

  SUBROUTINE case_law_release_cap()
    !! Segment 3 carries a tension-only law through the load callback (ea(3) = 0) and
    !! releases deeply within the step.
    REAL(wp), PARAMETER :: L0(3) = [1.0_wp, 1.2_wp, 0.8_wp], EA(3) = [1.0e6_wp, 1.5e6_wp, 0.0_wp]
    REAL(wp), PARAMETER :: EA_LAW = 0.8e6_wp, ALPHA_F = 0.4_wp/1.4_wp
    REAL(wp) :: qn(12), vn(12), an(12), q1(12), R(12), Rnone(12), Rel(12), K(9, 9), f_cap(12), f_raw(12)
    REAL(wp) :: t0, scale
    INTEGER, PARAMETER :: I8 = SELECTED_INT_KIND(18)
    CALL reset_chain_load()
    cb_law_elem = 3
    cb_ea_law = EA_LAW
    cb_l0_law = L0(3)
    CALL chain_states([0.010_wp, 0.020_wp, 0.010_wp], [0.020_wp, 0.015_wp, -0.050_wp], L0, qn, vn, an, q1)
    t0 = EA_LAW*0.010_wp
    ! nearly massless chain at rest at t_n: the residual is the blended force
    CALL probe_chain(qn, 0.0_wp*vn, 0.0_wp*an, q1, L0, EA, [1.0e-9_wp, 1.0e-9_wp, 1.0e-9_wp], .TRUE., R, K, &
                     start_tension=[0.0_wp, 0.0_wp, t0])
    CALL probe_chain(qn, 0.0_wp*vn, 0.0_wp*an, q1, L0, EA, [1.0e-9_wp, 1.0e-9_wp, 1.0e-9_wp], .TRUE., Rnone, K)
    ! the entries of the elastic segments are ignored
    CALL probe_chain(qn, 0.0_wp*vn, 0.0_wp*an, q1, L0, EA, [1.0e-9_wp, 1.0e-9_wp, 1.0e-9_wp], .TRUE., Rel, K, &
                     start_tension=[99.0_wp, -5.0_wp, t0])
    CALL blended_force(qn, q1, L0, [EA(1), EA(2), EA_LAW], ALPHA_F, .TRUE., f_cap)
    CALL blended_force(qn, q1, L0, [EA(1), EA(2), EA_LAW], ALPHA_F, .FALSE., f_raw)
    scale = MAXVAL(ABS(f_raw))
    CALL require(MAXVAL(ABS(f_cap - f_raw)) > 1.0e-2_wp*scale, 'law-cap:cap-active')
    CALL require(MAXVAL(ABS(R(4:12) - f_cap(4:12))) < 1.0e-6_wp*scale, 'law-cap:residual-closed-form')
    CALL require(MAXVAL(ABS(Rnone(4:12) - f_raw(4:12))) < 1.0e-6_wp*scale, 'law-cap:no-start-tension-uncapped')
    CALL require(ALL(TRANSFER(Rel, [0_I8]) == TRANSFER(R, [0_I8])), 'law-cap:elastic-entries-ignored')
    ! the tangent with inertia, on both paths
    CALL fd_tangent_check('law-cap', qn, vn, an, q1, L0, EA, [2.0_wp, 3.0_wp, 1.5_wp], .TRUE., &
                          start_tension=[0.0_wp, 0.0_wp, t0])
    CALL fd_tangent_check('law-cap', qn, vn, an, q1, L0, EA, [2.0_wp, 3.0_wp, 1.5_wp], .FALSE., &
                          start_tension=[0.0_wp, 0.0_wp, t0])
    CALL reset_chain_load()
  END SUBROUTINE case_law_release_cap

  SUBROUTINE case_load_share_tangent()
    !! All segments elastic and taut; the callback load carries a quadratic drag and a cubic
    !! position load (and an added mass on the dense path).
    REAL(wp), PARAMETER :: L0(3) = [1.0_wp, 1.2_wp, 0.8_wp], EA(3) = [1.0e6_wp, 1.5e6_wp, 0.8e6_wp]
    REAL(wp), PARAMETER :: RHO = 0.4_wp, DT = 0.05_wp
    REAL(wp) :: qn(12), vn(12), an(12), q1(12), v1(12), R(12), K(9, 9), f_int(12), f1(12), fn0(12)
    REAL(wp) :: jq(12, 12), jv(12, 12), alpha_m, alpha_f, beta, gamma, scale, load_scale
    INTEGER :: i
    CALL reset_chain_load()
    CALL chain_states([0.010_wp, 0.020_wp, 0.010_wp], [0.020_wp, 0.015_wp, 0.020_wp], L0, qn, vn, an, q1)
    cb_drag = 2.0e4_wp
    cb_spring = 1.0e8_wp
    DO i = 1, 12
      cb_q0(i) = qn(i) + 0.05_wp*SIN(0.9_wp*REAL(i, wp))
    END DO
    alpha_m = (2.0_wp*RHO - 1.0_wp)/(RHO + 1.0_wp)
    alpha_f = RHO/(RHO + 1.0_wp)
    beta = 0.25_wp*(1.0_wp - alpha_m + alpha_f)**2
    gamma = 0.5_wp - alpha_m + alpha_f
    ! closed form, nearly massless chain at rest at t_n: v_{n+1} = gamma/(beta dt) (q1 - qn) at
    ! the free nodes, the fixed node held
    v1 = gamma/(beta*DT)*(q1 - qn)
    v1(1:3) = 0.0_wp
    CALL probe_chain(qn, 0.0_wp*vn, 0.0_wp*an, q1, L0, EA, [1.0e-9_wp, 1.0e-9_wp, 1.0e-9_wp], .TRUE., R, K)
    CALL blended_force(qn, q1, L0, EA, alpha_f, .TRUE., f_int)
    CALL chain_load(q1, v1, f1, jq, jv)
    CALL chain_load(qn, 0.0_wp*vn, fn0, jq, jv)
    scale = MAXVAL(ABS(f_int))
    load_scale = MAXVAL(ABS(f1(4:12)))
    CALL require(load_scale > 0.1_wp*scale, 'load-share:load-significant')
    CALL require(MAXVAL(ABS(R(4:12) - (f_int(4:12) - (1.0_wp - alpha_f)*f1(4:12) - alpha_f*fn0(4:12)))) < &
                 1.0e-6_wp*scale, 'load-share:residual-closed-form')
    ! the tangent with inertia: banded, and dense with the added mass
    CALL fd_tangent_check('load-share', qn, vn, an, q1, L0, EA, [2.0_wp, 3.0_wp, 1.5_wp], .TRUE.)
    cb_madd = 50.0_wp
    CALL fd_tangent_check('load-share+added-mass', qn, vn, an, q1, L0, EA, [2.0_wp, 3.0_wp, 1.5_wp], .FALSE.)
    CALL reset_chain_load()
  END SUBROUTINE case_load_share_tangent

  SUBROUTINE solve_chain(work, qn, vn, an, ends, dt, start_load, start_tension, q_out, v_out, a_out, n_iter, es)
    !! One banded step of the chain with node 1 prescribed at ends (at rest).
    TYPE(CD_CableGenAlphaWorkspace), INTENT(INOUT) :: work
    REAL(wp), INTENT(IN) :: qn(12), vn(12), an(12), ends(3), dt, start_load(12), start_tension(3)
    REAL(wp), INTENT(OUT) :: q_out(12), v_out(12), a_out(12)
    INTEGER, INTENT(OUT) :: n_iter, es
    REAL(wp), PARAMETER :: L0(3) = [1.0_wp, 1.2_wp, 0.8_wp], EA(3) = [1.0e6_wp, 1.5e6_wp, 0.0_wp]
    INTEGER :: conn(2, 3)
    REAL(wp) :: pq(12)
    LOGICAL :: converged, stalled
    TYPE(GenAlphaConfig) :: cfg
    CHARACTER(120) :: em
    conn = RESHAPE([1, 2, 2, 3, 3, 4], [2, 3])
    cfg%rho_inf = 0.4_wp
    pq = qn
    pq(1:3) = ends
    CALL CD_Cable_Gen_Alpha_Step(qn, vn, an, conn, L0, EA, [2.0_wp, 3.0_wp, 1.5_wp], .TRUE., 0.0_wp*qn, [1, 2, 3], &
                                 dt, cfg, q_out, v_out, a_out, converged, stalled, n_iter, es, em, &
                                 load_band_proc=cb_load_band, load_force_proc=cb_force, prescribed_q=pq, &
                                 prescribed_v=0.0_wp*qn, prescribed_a=0.0_wp*qn, workspace=work, &
                                 start_load=start_load, start_tension=start_tension)
    IF (es == CD_DYN_OK .AND. .NOT. converged) es = -1
  END SUBROUTINE solve_chain

  LOGICAL FUNCTION same_bits(x, y) RESULT(same)
    REAL(wp), INTENT(IN) :: x(:), y(:)
    INTEGER, PARAMETER :: I8 = SELECTED_INT_KIND(18)
    same = ALL(TRANSFER(x, [0_I8]) == TRANSFER(y, [0_I8]))
  END FUNCTION same_bits

  SUBROUTINE case_warm_resolve_reuse()
    !! A junction-style re-solve of the same step (guess_armed) from the same committed state
    !! with moved ends: the step-start data of the first call are reused.
    REAL(wp), PARAMETER :: L0(3) = [1.0_wp, 1.2_wp, 0.8_wp], EA_LAW = 0.8e6_wp, ALPHA_F = 0.4_wp/1.4_wp
    REAL(wp), PARAMETER :: DIRS(3, 3) = RESHAPE([1.0_wp, 0.2_wp, 0.1_wp, 0.3_wp, 1.0_wp, -0.2_wp, &
                                                 -0.2_wp, 0.4_wp, 1.0_wp], [3, 3])
    TYPE(CD_CableGenAlphaWorkspace) :: warm, cold, wrong, cold2
    REAL(wp) :: qn(12), vn(12), an(12), qn2(12), sl(12), sl2(12), ends_a(3), ends_b(3), u3(3)
    REAL(wp) :: qw(12), vw(12), aw(12), qc(12), vc(12), ac(12), qx(12), vx(12), ax(12), t0, t02, e0, e1
    REAL(wp) :: jq(12, 12), jv(12, 12)
    INTEGER :: itw, itc, itx, es
    CALL reset_chain_load()
    cb_law_elem = 3
    cb_ea_law = EA_LAW
    cb_l0_law = L0(3)
    cb_drag = 500.0_wp
    cb_spring = 1.0e6_wp
    CALL chain_nodes(DIRS, [0.010_wp, 0.020_wp, 0.010_wp], L0, qn)
    cb_q0 = qn
    ! node 4 runs toward node 3, so segment 3 (the law) releases deeply within the step
    u3 = (qn(10:12) - qn(7:9))/NORM2(qn(10:12) - qn(7:9))
    vn = 0.0_wp
    vn(10:12) = -1.5_wp*u3
    vn(4:6) = [0.05_wp, -0.02_wp, 0.01_wp]
    an = 0.0_wp
    ends_a = qn(1:3) + [0.01_wp, 0.0_wp, 0.0_wp]
    ends_b = qn(1:3) + [0.0_wp, 0.012_wp, -0.004_wp]
    CALL chain_load(qn, vn, sl, jq, jv)
    t0 = EA_LAW*0.010_wp
    ! first call (cold), then the armed re-solve given wrong start data, against a cold solve
    CALL solve_chain(warm, qn, vn, an, ends_a, 0.05_wp, sl, [0.0_wp, 0.0_wp, t0], qw, vw, aw, itw, es)
    CALL require(es == CD_DYN_OK, 'reuse:first-call')
    warm%guess_armed = .TRUE.
    CALL solve_chain(warm, qn, vn, an, ends_b, 0.05_wp, sl + 1.0e3_wp, [0.0_wp, 0.0_wp, 5.0_wp*t0], qw, vw, aw, &
                     itw, es)
    CALL require(es == CD_DYN_OK, 'reuse:warm-call')
    CALL solve_chain(cold, qn, vn, an, ends_b, 0.05_wp, sl, [0.0_wp, 0.0_wp, t0], qc, vc, ac, itc, es)
    CALL require(es == CD_DYN_OK, 'reuse:cold-call')
    CALL require(itw == itc .AND. same_bits(qw, qc) .AND. same_bits(vw, vc) .AND. same_bits(aw, ac), &
                 'reuse:warm-equals-cold-bitwise')
    ! the cap is active at the solution, and the wrong start data would change it
    e0 = NORM2(qn(10:12) - qn(7:9)) - L0(3)
    e1 = NORM2(qc(10:12) - qc(7:9)) - L0(3)
    CALL require(e0 - e1 > e0/(2.0_wp*ALPHA_F), 'reuse:cap-active')
    CALL solve_chain(wrong, qn, vn, an, ends_b, 0.05_wp, sl + 1.0e3_wp, [0.0_wp, 0.0_wp, 5.0_wp*t0], qx, vx, ax, &
                     itx, es)
    CALL require(es == CD_DYN_OK .AND. .NOT. same_bits(qx, qc), 'reuse:start-data-matter')
    ! an unarmed call that fails before computing its start data (dt < 0), then an armed call
    ! from another state: it must not reuse the start data of the first state
    qn2 = qn
    qn2(7:9) = qn2(7:9) + [0.003_wp, -0.002_wp, 0.001_wp]
    CALL chain_load(qn2, vn, sl2, jq, jv)
    t02 = EA_LAW*(NORM2(qn2(10:12) - qn2(7:9))/L0(3) - 1.0_wp)
    CALL solve_chain(warm, qn2, vn, an, ends_b, -1.0_wp, sl2, [0.0_wp, 0.0_wp, t02], qx, vx, ax, itx, es)
    CALL require(es /= CD_DYN_OK, 'reuse:early-failure')
    warm%guess_armed = .TRUE.
    CALL solve_chain(warm, qn2, vn, an, ends_b, 0.05_wp, sl2, [0.0_wp, 0.0_wp, t02], qw, vw, aw, itw, es)
    CALL require(es == CD_DYN_OK, 'reuse:armed-after-failure')
    CALL solve_chain(cold2, qn2, vn, an, ends_b, 0.05_wp, sl2, [0.0_wp, 0.0_wp, t02], qc, vc, ac, itc, es)
    CALL require(es == CD_DYN_OK .AND. itw == itc .AND. same_bits(qw, qc) .AND. same_bits(vw, vc) .AND. &
                 same_bits(aw, ac), 'reuse:no-stale-start-data')
    CALL reset_chain_load()
  END SUBROUTINE case_warm_resolve_reuse

  ! ------------------------------------------------------------------------------------------
  ! Case 9: viscoelastic and Syrope lines that go slack (driver runs)
  ! ------------------------------------------------------------------------------------------

  SUBROUTINE write_syrope_files()
    !! The settings and OWC table of the shipped Syrope example (examples/data/syrope).
    REAL(wp), PARAMETER :: OWC(2, 30) = RESHAPE([ &
                                                0.00000e+00_wp, 0.00000e+00_wp, 2.06897e-03_wp, 1.71768e+05_wp, &
                                                4.13793e-03_wp, 3.30952e+05_wp, 6.20690e-03_wp, 4.78788e+05_wp, &
                                                8.27586e-03_wp, 6.16510e+05_wp, 1.03448e-02_wp, 7.45355e+05_wp, &
                                                1.24138e-02_wp, 8.66556e+05_wp, 1.44828e-02_wp, 9.81351e+05_wp, &
                                                1.65517e-02_wp, 1.09097e+06_wp, 1.86207e-02_wp, 1.19666e+06_wp, &
                                                2.06897e-02_wp, 1.29964e+06_wp, 2.27586e-02_wp, 1.40116e+06_wp, &
                                                2.48276e-02_wp, 1.50245e+06_wp, 2.68966e-02_wp, 1.60474e+06_wp, &
                                                2.89655e-02_wp, 1.70927e+06_wp, 3.10345e-02_wp, 1.81728e+06_wp, &
                                                3.31034e-02_wp, 1.93000e+06_wp, 3.51724e-02_wp, 2.04866e+06_wp, &
                                                3.72414e-02_wp, 2.17450e+06_wp, 3.93103e-02_wp, 2.30876e+06_wp, &
                                                4.13793e-02_wp, 2.45267e+06_wp, 4.34483e-02_wp, 2.60747e+06_wp, &
                                                4.55172e-02_wp, 2.77439e+06_wp, 4.75862e-02_wp, 2.95467e+06_wp, &
                                                4.96552e-02_wp, 3.14954e+06_wp, 5.17241e-02_wp, 3.36024e+06_wp, &
                                                5.37931e-02_wp, 3.58800e+06_wp, 5.58621e-02_wp, 3.83406e+06_wp, &
                                               5.79310e-02_wp, 4.09965e+06_wp, 6.00000e-02_wp, 4.38601e+06_wp], [2, 30])
    INTEGER :: u, i
    OPEN (NEWUNIT=u, FILE='ei0slack_syrope_settings.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'ei0slack_syrope_owc.dat  OWC     - Original working-curve table file'
    WRITE (u, '(A)') 'EXP             WCType  - Working-curve formulation'
    WRITE (u, '(A)') '0.20            k1      - First working-curve shape parameter'
    WRITE (u, '(A)') '1.50            k2      - Second working-curve shape parameter'
    CLOSE (u)
    OPEN (NEWUNIT=u, FILE='ei0slack_syrope_owc.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'Strain      Tension'
    WRITE (u, '(A)') '(-)         (N)'
    DO i = 1, 30
      WRITE (u, '(ES13.5,1X,ES13.5)') OWC(1, i), OWC(2, i)
    END DO
    CLOSE (u)
  END SUBROUTINE write_syrope_files

  SUBROUTINE write_syrope_deck(path, dt, tmax)
    !! A 230 m Syrope polyester line (115 segments, 80 kg/m) taut from the fairlead
    !! (5, 0, -14) to the anchor (227, 0, -80), working curve of a 1 MN running maximum.
    CHARACTER(*), INTENT(IN) :: path
    REAL(wp), INTENT(IN) :: dt, tmax
    INTEGER :: u
    OPEN (NEWUNIT=u, FILE=path, STATUS='REPLACE', ACTION='WRITE')
    WRITE (u, '(A)') 'EI = 0 slack-taut dynamics regression deck (Syrope)'
    WRITE (u, '(A)') '--- LINE TYPES ---'
    WRITE (u, '(A)') 'Name Diam Mass EA BA EI Cdn Cdt Can Cat'
    WRITE (u, '(A)') '(-) (m) (kg/m) (N) (-) (Nm2) (-) (-) (-) (-)'
    WRITE (u, '(A)') 'rope 0.2 80.0 "SYROPE:ei0slack_syrope_settings.dat|1.53e8|23.12" 1.0e8|1.0e5 0.0 1.2 0.1 1.0 0.0'
    WRITE (u, '(A)') '--- POINTS ---'
    WRITE (u, '(A)') 'ID Type X Y Z Mass Vol CdA Ca'
    WRITE (u, '(A)') '(-) (-) (m) (m) (m) (kg) (m^3) (m^2) (-)'
    WRITE (u, '(A)') '1 Coupled 5.0 0.0 -14.0 0 0 0 0'
    WRITE (u, '(A)') '2 Fixed 227.0 0.0 -80.0 0 0 0 0'
    WRITE (u, '(A)') '--- LINES ---'
    WRITE (u, '(A)') 'ID NodeA NodeB Outputs'
    WRITE (u, '(A)') '(-) (-) (-) (-)'
    WRITE (u, '(A)') '1 1 2 t'
    WRITE (u, '(A)') '--- SYROPE IC ---'
    WRITE (u, '(A)') 'Line(s) Tmax0 Tmean0'
    WRITE (u, '(A)') '(-) (N) (N)'
    WRITE (u, '(A)') '1 1.0e6 5.0e5'
    WRITE (u, '(A)') '--- SECTIONS ---'
    WRITE (u, '(A)') 'LineID LineType Length NumSegs'
    WRITE (u, '(A)') '(-) (-) (m) (-)'
    WRITE (u, '(A)') '1 rope 230.0 115'
    WRITE (u, '(A)') '--- OPTIONS ---'
    WRITE (u, '(A)') '80.0 WtrDpth'
    WRITE (u, '(ES14.7,A)') dt, ' dtM'
    WRITE (u, '(ES14.7,A)') tmax, ' TMax'
    WRITE (u, '(A)') 'ei0slack_surge_syr.txt motionFile'
    WRITE (u, '(A)') '--- OUTPUTS ---'
    WRITE (u, '(A)') 'FairTen1 AnchTen1'
    WRITE (u, '(A)') '--- end ---'
    CLOSE (u)
  END SUBROUTINE write_syrope_deck

  REAL(wp) FUNCTION column_min(path, first_col, last_col) RESULT(mn)
    !! Smallest value of the columns first_col..last_col (last_col 0: through the last column).
    CHARACTER(*), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: first_col, last_col
    REAL(wp), ALLOCATABLE :: d(:, :)
    INTEGER :: nrow, lc
    mn = HUGE(1.0_wp)
    CALL read_out(path, d, nrow)
    IF (nrow < 1) RETURN
    lc = last_col
    IF (lc == 0) lc = SIZE(d, 1)
    mn = MINVAL(d(first_col:lc, 1:nrow))
  END FUNCTION column_min

  SUBROUTINE case_law_slack_runs()
    LOGICAL :: ok
    CHARACTER(64) :: opts(4)
    REAL(wp) :: fair, anch, seg_min
    ! the light cable as a lightly damped viscoelastic line, 3 m heave at dtM 0.1 s: anchor
    ! 7.6e4 to 9.4e4 N and fairlead 0.97e5 to 1.08e5 N at dtM 0.5 to 2.5 ms; without the cap
    ! on its law segments 1.45e5 N at the anchor
    CALL write_motion('ei0slack_heave_ve.txt', 1, [5.0_wp, 0.0_wp, -14.0_wp], 3, 3.0_wp, 10.0_wp, 10.0_wp, &
                      0.1_wp, 40.0_wp)
    opts(1) = '0.1 dtM'
    opts(2) = '40.0 TMax'
    opts(3) = 'ei0slack_heave_ve.txt motionFile'
    opts(4) = ''
    CALL write_deck('ei0slack_ve_light.dat', CABLE_VE, [205.0_wp, 0.0_wp, -80.0_wp], [5.0_wp, 0.0_wp, -14.0_wp], &
                    230.0_wp, 115, 80.0_wp, 't', opts)
    CALL run_driver_logged('ei0slack_ve_light.dat', 'ei0slack_ve_light', 'ei0slack_ve_light.log', ok)
    CALL require(ok, 've-slack:run')
    fair = column_max('ei0slack_ve_light.out', 2, 2)
    anch = column_max('ei0slack_ve_light.out', 3, 3)
    seg_min = column_min('ei0slack_ve_light.Line1.t.out', 2, 0)
    WRITE (*, '(A,ES10.3,A,ES10.3,A,ES10.3,A)') '  viscoelastic light cable, dtM 0.1 s: FairTen1 max ', fair, &
      ' N, AnchTen1 max ', anch, ' N, smallest segment tension ', seg_min, ' N'
    CALL require(seg_min <= 0.0_wp, 've-slack:goes-slack')
    CALL require(anch > 6.5e4_wp .AND. anch < 1.0e5_wp, 've-slack:anchor-peak')
    CALL require(fair > 8.0e4_wp .AND. fair < 1.2e5_wp, 've-slack:fairlead-peak')
    CALL require(.NOT. file_contains('ei0slack_ve_light.log', 'Warning: line'), 've-slack:no-spike-warning')
    ! a taut Syrope line in 4 m surge that goes slack, dtM 0.1 s: fairlead 1.78 MN at dtM 1 and
    ! 2.5 ms
    CALL write_syrope_files()
    CALL write_motion('ei0slack_surge_syr.txt', 1, [5.0_wp, 0.0_wp, -14.0_wp], 1, 4.0_wp, 10.0_wp, 10.0_wp, &
                      0.1_wp, 40.0_wp)
    CALL write_syrope_deck('ei0slack_syrope.dat', 0.1_wp, 40.0_wp)
    CALL run_driver_logged('ei0slack_syrope.dat', 'ei0slack_syrope', 'ei0slack_syrope.log', ok)
    CALL require(ok, 'syrope-slack:run')
    fair = column_max('ei0slack_syrope.out', 2, 2)
    seg_min = column_min('ei0slack_syrope.Line1.t.out', 2, 0)
    WRITE (*, '(A,ES10.3,A,ES10.3,A)') '  Syrope line, 4 m surge, dtM 0.1 s: FairTen1 max ', fair, &
      ' N, smallest segment tension ', seg_min, ' N'
    CALL require(seg_min <= 0.0_wp, 'syrope-slack:goes-slack')
    CALL require(fair > 1.5e6_wp .AND. fair < 1.95e6_wp, 'syrope-slack:fairlead-peak')
    CALL require(.NOT. file_contains('ei0slack_syrope.log', 'Warning: line'), 'syrope-slack:no-spike-warning')
  END SUBROUTINE case_law_slack_runs

END PROGRAM test_ei0_slack_dynamics
