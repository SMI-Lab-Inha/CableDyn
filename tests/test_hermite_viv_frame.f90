! File: tests/test_hermite_viv_frame.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_hermite_viv_frame
  !! The flow frame of the synchronisation model on the cubic-Hermite cable
  !! (CD_HermiteCable_Dyn_Set_VIV tf): with a filter time Tf > 0 the amplitude estimate and the
  !! cross-flow damping split take the ambient flow relative to the node velocity low-pass
  !! filtered over Tf, so these quantities are the same in the cylinder's frame and in the water's.
  !!
  !! The pipe: 20 m, 10 elements, D 0.1 m, EI 1000 N m^2, pinned, about 21 kN, its mass that of the
  !! water it displaces (8.05 kg/m, no net weight), Cd 1.2, Ca 1, Cl 0.8 with an amplitude table
  !! (Cv 0.8, 0.6, 0.0 and CdAmp 1.0, 1.3, 1.6 at A/D 0, 0.5, 1) and Thorsen's cross-flow damping
  !! (a, b) = (0.31, 0.89). The flow along y oscillates with amplitude 0.6 m/s and period 40 s
  !! (KC 240), its amplitude ramped up over the first period; the cross-flow axis is z.
  !! * moving: still water, both ends moved along y with velocity V(t);
  !! * fixed: the ends held, the water moving at -V(t) (a held fluid field with its acceleration,
  !!   so the Froude-Krylov and added-mass loads act).
  !! With the pipe's mass equal to the displaced mass, the in-line inertia of the two problems is
  !! the same ((m + Ca rho A) dV/dt against (1 + Ca) rho A dV/dt); for any other mass they differ by
  !! the Froude-Krylov term (m - rho A) dV/dt, the inertial force of the moving frame. So the two
  !! runs differ only where the model depends on the frame.
  !!
  !! Gates:
  !! * F1 the setter refuses tf of the wrong size, a negative or non-finite tf, tf > 0 on an element
  !!   without an amplitude table, and tf with the wake oscillator; an all-zero tf installs no
  !!   frame and reproduces the run without tf bit for bit;
  !! * F2 with Tf 0.6 s the moving and the fixed pipe agree: the peak cross-flow response of the
  !!   mid node over the second period within 5 % (measured 1.0-2.2 % over the Release and Debug
  !!   builds and the two OpenBLAS kernels), its RMS within 12 % (measured 7 %). Two differences
  !!   are left. The filter lags the imposed velocity by about omega Tf = 0.094 rad (5 degrees).
  !!   The phase rate keeps the published form: its cross-flow velocity is the node's velocity
  !!   across the relative flow, r'.e = u.e, which is zero in still water, so the moving pipe's
  !!   phase follows y'' alone (the RMS gap, and the moving run carries round-off further; with y'
  !!   taken relative to the frame the two runs agree within 0.1 %);
  !! * F3 without the frame they differ: the moving pipe in still water takes Cv and CdAmp at A/D
  !!   0 and the Morison drag (no ambient flow); its peak response is more than 10 % (measured
  !!   16 %) and its RMS more than 15 % (measured 23 %) from the fixed pipe's;
  !! * F4 a snapshot taken mid-run and restored replays the following steps bit for bit, the
  !!   filtered velocity included.
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_VIV, ONLY: CD_VIV_AMP_MAXP, CD_VIV_NMS_FRAME
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Set_Drag, CD_HermiteCable_Dyn_Set_AddedMass, &
                                          CD_HermiteCable_Dyn_Set_VIV, CD_HermiteCable_Dyn_Set_Held_Fluid, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_Snapshot, &
                                          CD_HermiteCable_Dyn_Restore, CD_HermiteCable_Dyn_End, CD_HCDYN_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.141592653589793_wp
  REAL(wp), PARAMETER :: RHOW = 1025.0_wp, D = 0.1_wp, EA = 5.0e7_wp, EI = 1.0e3_wp, TENS = 2.1e4_wp
  REAL(wp), PARAMETER :: SPAN = 20.0_wp, DT = 0.01_wp, U0 = 0.6_wp, TP = 40.0_wp, TF = 0.6_wp
  REAL(wp), PARAMETER :: MLIN = RHOW*0.25_wp*PI*D*D
  INTEGER, PARAMETER :: NE = 10, NN = NE + 1, MID = 6*(NE/2), NSTEP = NINT(2.0_wp*TP/DT)

  TYPE(CD_HermiteCableDynType) :: model, twin
  REAL(wp) :: cl(NE), cf(NE), df(NE), amp(3, CD_VIV_AMP_MAXP, NE), cdy(2, NE), tf0(NE), tfv(NE)
  INTEGER :: ampn(NE)
  REAL(wp) :: pk_mf, rms_mf, pk_ff, rms_ff, pk_m0, rms_m0, pk_f0, rms_f0
  INTEGER :: es, nfail, k
  CHARACTER(1024) :: em

  nfail = 0
  cl = 0.8_wp; cf = 0.18_wp; df = 0.08_wp
  amp = 0.0_wp
  amp(:, 1, :) = SPREAD([0.0_wp, 0.8_wp, 1.0_wp], 2, NE)
  amp(:, 2, :) = SPREAD([0.5_wp, 0.6_wp, 1.3_wp], 2, NE)
  amp(:, 3, :) = SPREAD([1.0_wp, 0.0_wp, 1.6_wp], 2, NE)
  ampn = 3
  cdy(1, :) = 0.31_wp; cdy(2, :) = 0.89_wp
  tf0 = 0.0_wp
  tfv = TF

  ! ---- F1 refusals; an all-zero tf is no frame ----
  CALL build(model)
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 2.0_wp, es, em, amp=amp, amp_n=ampn, tf=tfv(1:NE - 1))
  CALL require(es /= CD_HCDYN_OK, 'F1: tf of the wrong size is refused')
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 2.0_wp, es, em, amp=amp, amp_n=ampn, tf=-tfv)
  CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'must be >= 0') > 0, 'F1: a negative tf is refused: '//TRIM(em))
  BLOCK
    REAL(wp) :: bad(NE)
    INTEGER :: an(NE)
    bad = tfv
    bad(4) = IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN)
    CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 2.0_wp, es, em, amp=amp, amp_n=ampn, tf=bad)
    CALL require(es /= CD_HCDYN_OK, 'F1: a non-finite tf is refused')
    CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 2.0_wp, es, em, tf=tfv)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'needs an amplitude table') > 0, &
                 'F1: tf > 0 without amplitude tables is refused: '//TRIM(em))
    an = ampn
    an(3) = 0
    CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 2.0_wp, es, em, amp=amp, amp_n=an, tf=tfv)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'needs an amplitude table') > 0, &
                 'F1: tf > 0 on an element without its table is refused: '//TRIM(em))
    CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 2.0_wp, es, em, &
                                     wake=RESHAPE([(0.2_wp, 0.3_wp, 12.0_wp, k=1, NE)], [3, NE]), tf=tf0)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'wake(4, :)') > 0, 'F1: tf with the wake oscillator is refused: '// &
                 TRIM(em))
  END BLOCK
  CALL require(.NOT. model%has_viv .AND. .NOT. model%viv_frame, 'F1: a refused call installs nothing')
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 2.0_wp, es, em, amp=amp, amp_n=ampn, cdy=cdy, tf=tf0)
  CALL require(es == CD_HCDYN_OK .AND. .NOT. model%viv_frame .AND. SIZE(model%viv_ms, 1) == 6, &
               'F1: an all-zero tf installs no frame: '//TRIM(em))
  CALL build(twin)
  CALL CD_HermiteCable_Dyn_Set_VIV(twin, cl, cf, df, 2.0_wp, es, em, amp=amp, amp_n=ampn, cdy=cdy)
  CALL require(es == CD_HCDYN_OK, 'F1: the twin without tf installs: '//TRIM(em))
  DO k = 1, 600
    CALL advance(model, k, .TRUE.)
    CALL advance(twin, k, .TRUE.)
    IF (es /= CD_HCDYN_OK) EXIT
    IF (ANY(model%q /= twin%q) .OR. ANY(model%viv_phi /= twin%viv_phi) .OR. ANY(model%viv_ms /= twin%viv_ms)) EXIT
  END DO
  CALL require(k > 600, 'F1: an all-zero tf reproduces the run without tf bit for bit')
  CALL CD_HermiteCable_Dyn_End(model)
  CALL CD_HermiteCable_Dyn_End(twin)

  ! ---- F2 and F3 the two frames, with and without the flow frame ----
  CALL run_case(.TRUE., tfv, pk_mf, rms_mf)
  CALL run_case(.FALSE., tfv, pk_ff, rms_ff)
  CALL run_case(.TRUE., tf0, pk_m0, rms_m0)
  CALL run_case(.FALSE., tf0, pk_f0, rms_f0)
  WRITE (*, '(A,2F9.4,A,2F9.4)') '  F2 with Tf: moving peak, RMS A/D ', pk_mf/D, rms_mf/D, '; fixed ', pk_ff/D, rms_ff/D
  WRITE (*, '(A,2F9.4,A,2F9.4)') '  F3 without: moving peak, RMS A/D ', pk_m0/D, rms_m0/D, '; fixed ', pk_f0/D, rms_f0/D
  CALL require(pk_ff > 0.2_wp*D, 'F2: the fixed pipe locks in')
  CALL require(ABS(pk_mf/pk_ff - 1.0_wp) <= 0.05_wp .AND. ABS(rms_mf/rms_ff - 1.0_wp) <= 0.12_wp, &
               'F2: with the flow frame the moving and the fixed pipe agree')
  CALL require(ABS(pk_m0/pk_f0 - 1.0_wp) > 0.10_wp .AND. ABS(rms_m0/rms_f0 - 1.0_wp) > 0.15_wp, &
               'F3: without it the moving pipe responds otherwise')

  ! ---- F4 snapshot replay with the filtered velocity ----
  BLOCK
    REAL(wp), ALLOCATABLE :: q_ref(:, :), ms_ref(:, :, :)
    INTEGER :: s
    CALL build(model)
    CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 2.0_wp, es, em, amp=amp, amp_n=ampn, cdy=cdy, tf=tfv)
    CALL require(es == CD_HCDYN_OK .AND. model%viv_frame .AND. SIZE(model%viv_ms, 1) == CD_VIV_NMS_FRAME, &
                 'F4: the frame installs: '//TRIM(em))
    DO s = 1, 1500
      CALL advance(model, s, .TRUE.)
    END DO
    CALL require(MAXVAL(ABS(model%viv_ms(8, :))) > 0.1_wp*U0, 'F4: the filtered velocity follows the motion')
    CALL CD_HermiteCable_Dyn_Snapshot(model, es, em)
    CALL require(es == CD_HCDYN_OK, 'F4: snapshot: '//TRIM(em))
    ALLOCATE (q_ref(SIZE(model%q), 200), ms_ref(CD_VIV_NMS_FRAME, NN, 200))
    DO s = 1, 200
      CALL advance(model, 1500 + s, .TRUE.)
      q_ref(:, s) = model%q
      ms_ref(:, :, s) = model%viv_ms
    END DO
    CALL CD_HermiteCable_Dyn_Restore(model, es, em)
    CALL require(es == CD_HCDYN_OK, 'F4: restore: '//TRIM(em))
    DO s = 1, 200
      CALL advance(model, 1500 + s, .TRUE.)
      IF (ANY(model%q /= q_ref(:, s)) .OR. ANY(model%viv_ms /= ms_ref(:, :, s))) EXIT
    END DO
    CALL require(s > 200, 'F4: the restored snapshot replays the state and the filtered velocity bit for bit')
    CALL CD_HermiteCable_Dyn_End(model)
  END BLOCK

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: the flow frame of the synchronisation model on the cubic-Hermite cable'

CONTAINS

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (cond) RETURN
    nfail = nfail + 1
    WRITE (*, '(A)') 'FAIL: '//TRIM(label)
  END SUBROUTINE require

  SUBROUTINE build(m)
    !! The pinned taut pipe along x at 10 m depth in still water, pretensioned through its
    !! unstretched length, with drag and the normal added mass.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: m
    REAL(wp) :: seed(6*NN), l0(NE), eav(NE), eiv(NE), wv(NE), rhoa(NE), hd(NE), hn(NE), ht(NE)
    INTEGER :: i
    seed = 0.0_wp
    DO i = 1, NN
      seed(6*(i - 1) + 1) = SPAN*REAL(i - 1, wp)/REAL(NE, wp)
      seed(6*(i - 1) + 3) = -10.0_wp
      seed(6*(i - 1) + 4) = 1.0_wp
    END DO
    l0 = SPAN/REAL(NE, wp)/(1.0_wp + TENS/EA)
    eav = EA; eiv = EI; rhoa = MLIN; wv = 0.0_wp
    CALL CD_HermiteCable_Dyn_Init(m, l0, eav, eiv, rhoa, wv, seed, end_dofs(), -200.0_wp, 0.0_wp, 0.5_wp, es, em)
    CALL require(es == CD_HCDYN_OK, 'init: '//TRIM(em))
    hd = D; hn = 1.2_wp; ht = 0.0_wp
    CALL CD_HermiteCable_Dyn_Set_Drag(m, RHOW, hd, hn, ht, 0.0_wp, [0.0_wp, 0.0_wp, 0.0_wp], es, em)
    CALL require(es == CD_HCDYN_OK, 'drag: '//TRIM(em))
    hn = 1.0_wp
    CALL CD_HermiteCable_Dyn_Set_AddedMass(m, RHOW, hd, hn, ht, 0.0_wp, es, em)
    CALL require(es == CD_HCDYN_OK, 'added mass: '//TRIM(em))
  END SUBROUTINE build

  PURE FUNCTION end_dofs() RESULT(fixed)
    INTEGER :: fixed(6)
    fixed = [1, 2, 3, 6*NE + 1, 6*NE + 2, 6*NE + 3]
  END FUNCTION end_dofs

  SUBROUTINE imposed(t, y, v, a)
    !! The imposed motion along y: V = U0 r(t) sin(omega t), r = (1 - cos(pi t/TP))/2 over the
    !! first period and 1 after it; y its integral (exact), a its derivative.
    REAL(wp), INTENT(IN) :: t
    REAL(wp), INTENT(OUT) :: y, v, a
    REAL(wp) :: om, p, r, dr
    om = 2.0_wp*PI/TP
    p = PI/TP
    IF (t < TP) THEN
      r = 0.5_wp*(1.0_wp - COS(p*t))
      dr = 0.5_wp*p*SIN(p*t)
      ! int_0^t (1 - cos(p s)) sin(om s)/2 ds with om = 2 p
      y = U0*0.5_wp*((1.0_wp - COS(om*t))/om - (0.5_wp*(1.0_wp - COS(3.0_wp*p*t))/(3.0_wp*p) + &
                                                0.5_wp*(1.0_wp - COS(p*t))/p))
    ELSE
      r = 1.0_wp
      dr = 0.0_wp
      y = U0*0.5_wp*(-(0.5_wp*(1.0_wp - COS(3.0_wp*PI))/(3.0_wp*p) + 0.5_wp*(1.0_wp - COS(PI))/p)) + &
          U0*(COS(om*TP) - COS(om*t))/om
    END IF
    v = U0*r*SIN(om*t)
    a = U0*(dr*SIN(om*t) + r*om*COS(om*t))
  END SUBROUTINE imposed

  SUBROUTINE advance(m, s, moving)
    !! Step s (to t = s DT) of the moving pipe (ends prescribed) or of the fixed pipe (the held
    !! field at the step's mid time, -V and -dV/dt).
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: m
    INTEGER, INTENT(IN) :: s
    LOGICAL, INTENT(IN) :: moving
    REAL(wp) :: y, v, a, pq(6), pv(6), pa(6), fu(3, NN), fa(3, NN), wl(NN)
    INTEGER :: i
    IF (moving) THEN
      CALL imposed(REAL(s, wp)*DT, y, v, a)
      pq = [0.0_wp, y, -10.0_wp, SPAN, y, -10.0_wp]
      pv = [0.0_wp, v, 0.0_wp, 0.0_wp, v, 0.0_wp]
      pa = [0.0_wp, a, 0.0_wp, 0.0_wp, a, 0.0_wp]
      CALL CD_HermiteCable_Dyn_Step(m, DT, 50, 1.0e-8_wp, es, em, pres_dofs=end_dofs(), pres_q=pq, pres_v=pv, &
                                    pres_a=pa)
    ELSE
      CALL imposed((REAL(s, wp) - 0.5_wp)*DT, y, v, a)
      DO i = 1, NN
        fu(:, i) = [0.0_wp, -v, 0.0_wp]
        fa(:, i) = [0.0_wp, -a, 0.0_wp]
      END DO
      wl = 0.0_wp
      CALL CD_HermiteCable_Dyn_Set_Held_Fluid(m, fu, fa, wl, es, em)
      IF (es /= CD_HCDYN_OK) THEN
        CALL require(.FALSE., 'held field: '//TRIM(em)); RETURN
      END IF
      CALL CD_HermiteCable_Dyn_Step(m, DT, 50, 1.0e-8_wp, es, em)
    END IF
    IF (es /= CD_HCDYN_OK) CALL require(.FALSE., 'step: '//TRIM(em))
  END SUBROUTINE advance

  SUBROUTINE run_case(moving, tfe, pk, rms)
    !! Two periods of the moving or the fixed pipe with the filter times tfe; over the second
    !! period the peak and RMS cross-flow (z) displacement of the mid node from its rest depth.
    LOGICAL, INTENT(IN) :: moving
    REAL(wp), INTENT(IN) :: tfe(NE)
    REAL(wp), INTENT(OUT) :: pk, rms
    TYPE(CD_HermiteCableDynType), ALLOCATABLE :: m
    REAL(wp) :: z, s2
    INTEGER :: s, n
    ALLOCATE (m)
    CALL build(m)
    CALL CD_HermiteCable_Dyn_Set_VIV(m, cl, cf, df, 2.0_wp, es, em, amp=amp, amp_n=ampn, cdy=cdy, tf=tfe)
    CALL require(es == CD_HCDYN_OK, 'VIV: '//TRIM(em))
    pk = 0.0_wp
    s2 = 0.0_wp
    n = 0
    DO s = 1, NSTEP
      CALL advance(m, s, moving)
      IF (es /= CD_HCDYN_OK) EXIT
      IF (s <= NSTEP/2) CYCLE
      z = m%q(MID + 3) + 10.0_wp
      pk = MAX(pk, ABS(z))
      s2 = s2 + z*z
      n = n + 1
    END DO
    rms = SQRT(s2/REAL(MAX(n, 1), wp))
    CALL CD_HermiteCable_Dyn_End(m)
  END SUBROUTINE run_case

END PROGRAM test_hermite_viv_frame
