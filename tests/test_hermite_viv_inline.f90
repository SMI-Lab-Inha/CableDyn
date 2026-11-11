! File: tests/test_hermite_viv_inline.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_hermite_viv_inline
  !! In-line VIV on the cubic-Hermite cable (CD_HermiteCable_Dyn_Set_VIV clx, il_alpha; Ulveseter,
  !! Thorsen, Saevik and Larsen, Mar. Struct. 60, 2018): a second shedding force
  !! 1/2 rho D C_vx cos(phi_x) U v_n along the relative normal flow, its phase slaved to twice the
  !! cross-flow shedding rate, d phi_x/dt = 2 (d phi/dt) [1 + alpha sin(phi_xd - phi_x)].
  !!
  !! The line is test_hermite_viv's: 100 m taut, pinned, 20 elements, D 0.1 m, in a 0.15 m/s
  !! current along y, so the cross-flow lift acts along z and the in-line force along y.
  !!
  !! Gates:
  !! * I1 the setter refuses clx without il_alpha, mismatched sizes, a negative clx, alpha >= 1 and
  !!   an in-line coefficient on an element without the cross-flow lift; an all-zero clx installs
  !!   no in-line state; the in-line phase starts at twice the cross-flow phase;
  !! * I2 the element load: a straight element at rest across the current carries
  !!   U^2 l0 (k_x1 + k_x2)/2 along the flow; the in-line part adds to the cross-flow one; its
  !!   Jacobians match central differences on a curved, moving element;
  !! * I3 an all-zero clx reproduces the cross-flow-only run bit for bit;
  !! * I4 at rest the committed node load carries the in-line force along y, absent without it;
  !! * I5 in lock-in the in-line response at the quarter span oscillates at twice the cross-flow
  !!   frequency of the mid span, and is larger than the cross-flow-only run's;
  !! * I6 a snapshot taken mid-run and restored replays the following steps bit for bit (state and
  !!   both phases), and the recovering step equals the plain step.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCableDynType, CD_HermiteCable_Dyn_Init, &
                                          CD_HermiteCable_Dyn_Set_Drag, CD_HermiteCable_Dyn_Set_AddedMass, &
                                          CD_HermiteCable_Dyn_Set_VIV, CD_HermiteCable_Dyn_VIV_Node_Lift, &
                                          CD_HermiteCable_Dyn_Step, CD_HermiteCable_Dyn_Step_Recovering, &
                                          CD_HermiteCable_Dyn_Snapshot, CD_HermiteCable_Dyn_Restore, &
                                          CD_HermiteCable_Dyn_End, CD_HermiteCable_VIV_Lift_Element, CD_HCDYN_OK
  IMPLICIT NONE

  REAL(wp), PARAMETER :: PI = 3.141592653589793_wp
  REAL(wp), PARAMETER :: RHOW = 1025.0_wp, GACC = 9.80665_wp
  REAL(wp), PARAMETER :: D = 0.1_wp, MLIN = 8.10_wp, EA = 5.0e7_wp, EI = 1.0_wp, TENS = 5.0e4_wp
  REAL(wp), PARAMETER :: SPAN = 100.0_wp, U = 0.15_wp, DT = 0.02_wp, CVX = 0.85_wp, ALPHA = 0.15_wp
  INTEGER, PARAMETER :: NE = 20

  TYPE(CD_HermiteCableDynType) :: model, twin
  REAL(wp), ALLOCATABLE :: cl(:), cf(:), df(:), clx(:), al(:), q_ref(:, :), phi_ref(:, :), phx_ref(:, :)
  REAL(wp) :: lift(3, NE + 1)
  INTEGER :: es, nfail, k, mid, quarter
  CHARACTER(1024) :: em

  nfail = 0
  mid = 6*(NE/2)          ! base index of the mid-span node
  quarter = 6*(NE/4)      ! and of the quarter-span node
  ALLOCATE (cl(NE), cf(NE), df(NE), clx(NE), al(NE))
  cl = 0.8_wp; cf = 0.18_wp; df = 0.08_wp; clx = CVX; al = ALPHA

  ! ---- I1 refusals and the start ----
  CALL build(model)
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=clx)
  CALL require(es /= CD_HCDYN_OK, 'I1: clx without il_alpha is refused')
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=clx(1:NE - 1), il_alpha=al(1:NE - 1))
  CALL require(es /= CD_HCDYN_OK, 'I1: mismatched in-line arrays are refused')
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=-clx, il_alpha=al)
  CALL require(es /= CD_HCDYN_OK, 'I1: a negative clx is refused')
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=clx, il_alpha=al + 1.0_wp)
  CALL require(es /= CD_HCDYN_OK, 'I1: alpha >= 1 is refused')
  BLOCK
    REAL(wp) :: cl2(NE)
    cl2 = cl
    cl2(3) = 0.0_wp
    CALL CD_HermiteCable_Dyn_Set_VIV(model, cl2, cf, df, 0.25_wp, es, em, clx=clx, il_alpha=al)
    CALL require(es /= CD_HCDYN_OK .AND. INDEX(em, 'cross-flow') > 0, &
                 'I1: an in-line coefficient without the cross-flow lift is refused: '//TRIM(em))
  END BLOCK
  CALL require(.NOT. model%has_viv .AND. .NOT. model%has_viv_il, 'I1: a refused call installs nothing')
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=0.0_wp*clx, il_alpha=al)
  CALL require(es == CD_HCDYN_OK .AND. model%has_viv .AND. .NOT. model%has_viv_il, &
               'I1: an all-zero clx installs the cross-flow VIV only')
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=clx, il_alpha=al)
  CALL require(es == CD_HCDYN_OK .AND. model%has_viv_il, 'I1: in-line VIV installs: '//TRIM(em))
  CALL require(ALL(ABS(model%vil_phi - MODULO(2.0_wp*model%viv_phi, 2.0_wp*PI)) <= 1.0e-12_wp), &
               'I1: the in-line phase starts at twice the cross-flow phase')
  CALL require(ALL(ABS(model%vil_k - 0.5_wp*RHOW*D*CVX) <= 1.0e-12_wp), 'I1: the in-line scale is 1/2 rho D C_vx')

  ! ---- I4 the committed node load at rest: in-line along y, cross-flow along z ----
  CALL CD_HermiteCable_Dyn_VIV_Node_Lift(model, lift, es, em)
  CALL require(es == CD_HCDYN_OK, 'I4: node lift: '//TRIM(em))
  BLOCK
    REAL(wp) :: lift0(3, NE + 1), ref_y, term, scale
    INTEGER :: e
    ! the in-line force of each (stretched) element at rest: U^2 (SPAN/NE) (k_x1 + k_x2)/2 along y;
    ! the spread start phases nearly cancel in the sum, so the tolerance scales with its terms
    ref_y = 0.0_wp
    scale = 0.0_wp
    DO e = 1, NE
      term = U*U*(SPAN/REAL(NE, wp))*0.5_wp*model%vil_k(e)*(COS(model%vil_phi(e)) + COS(model%vil_phi(e + 1)))
      ref_y = ref_y + term
      scale = scale + ABS(term)
    END DO
    CALL require(ABS(SUM(lift(2, :)) - ref_y) <= 1.0e-12_wp*scale .AND. MAXVAL(ABS(lift(2, :))) > 1.0e-2_wp, &
                 'I4: the node load carries the in-line force along the flow')
    CALL twin_without_il(twin)
    CALL CD_HermiteCable_Dyn_VIV_Node_Lift(twin, lift0, es, em)
    CALL require(MAXVAL(ABS(lift0(2, :))) <= 1.0e-12_wp .AND. &
                 MAXVAL(ABS(lift0(3, :) - lift(3, :))) <= 1.0e-12_wp*MAXVAL(ABS(lift0(3, :))), &
                 'I4: without in-line VIV no force along y, the same cross-flow lift')
    CALL CD_HermiteCable_Dyn_End(twin)
  END BLOCK

  ! ---- I2 element load and Jacobians ----
  CALL check_element()

  ! ---- I3 an all-zero clx is the cross-flow-only run ----
  CALL CD_HermiteCable_Dyn_End(model)
  CALL build(model)
  CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=0.0_wp*clx, il_alpha=al)
  CALL twin_without_il(twin)
  DO k = 1, 300
    CALL advance(model, 1, .FALSE.)
    CALL advance(twin, 1, .FALSE.)
    IF (ANY(model%q /= twin%q) .OR. ANY(model%viv_phi /= twin%viv_phi)) EXIT
  END DO
  CALL require(k > 300, 'I3: an all-zero clx reproduces the cross-flow-only run bit for bit')
  CALL CD_HermiteCable_Dyn_End(model)
  CALL CD_HermiteCable_Dyn_End(twin)

  ! ---- I5 lock-in: the in-line response at twice the cross-flow frequency ----
  BLOCK
    REAL(wp) :: fy, fx, ax_il, ax_cf, ay
    CALL build(model)
    CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=clx, il_alpha=al)
    CALL run_response(model, fy, fx, ay, ax_il)
    CALL CD_HermiteCable_Dyn_End(model)
    CALL twin_without_il(twin)
    CALL run_response(twin, fy, fx, ay, ax_cf)
    CALL CD_HermiteCable_Dyn_End(twin)
    CALL build(model)
    CALL CD_HermiteCable_Dyn_Set_VIV(model, cl, cf, df, 0.25_wp, es, em, clx=clx, il_alpha=al)
    CALL run_response(model, fy, fx, ay, ax_il)
    WRITE (*, '(A,F7.4,A,F7.4,A,F7.3)') '  I5: cross-flow ', fy, ' Hz, in-line ', fx, ' Hz, ratio ', fx/fy
    WRITE (*, '(A,F8.4,A,F8.4,A,F8.4)') '  I5: CF amplitude/D ', ay/D, ', IL amplitude/D ', ax_il/D, &
      ', IL amplitude/D without in-line VIV ', ax_cf/D
    CALL require(ay > 0.3_wp*D, 'I5: the line locks in cross-flow')
    CALL require(ABS(fx/fy - 2.0_wp) <= 0.15_wp, 'I5: the in-line response is at twice the cross-flow frequency')
    CALL require(ax_il > 2.0_wp*ax_cf, 'I5: the in-line force raises the in-line response')
  END BLOCK

  ! ---- I6 snapshot replay and the recovering step ----
  CALL CD_HermiteCable_Dyn_Snapshot(model, es, em)
  CALL require(es == CD_HCDYN_OK, 'I6: snapshot: '//TRIM(em))
  ALLOCATE (q_ref(SIZE(model%q), 100), phi_ref(NE + 1, 100), phx_ref(NE + 1, 100))
  DO k = 1, 100
    CALL advance(model, 1, .FALSE.)
    q_ref(:, k) = model%q
    phi_ref(:, k) = model%viv_phi
    phx_ref(:, k) = model%vil_phi
  END DO
  CALL CD_HermiteCable_Dyn_Restore(model, es, em)
  CALL require(es == CD_HCDYN_OK, 'I6: restore: '//TRIM(em))
  DO k = 1, 100
    CALL advance(model, 1, .FALSE.)
    IF (ANY(model%q /= q_ref(:, k)) .OR. ANY(model%viv_phi /= phi_ref(:, k)) .OR. &
        ANY(model%vil_phi /= phx_ref(:, k))) EXIT
  END DO
  CALL require(k > 100, 'I6: the restored snapshot replays the state and both phases bit for bit')
  CALL CD_HermiteCable_Dyn_Restore(model, es, em)
  DO k = 1, 100
    CALL advance(model, 1, .TRUE.)
    IF (ANY(model%q /= q_ref(:, k)) .OR. ANY(model%viv_phi /= phi_ref(:, k)) .OR. &
        ANY(model%vil_phi /= phx_ref(:, k))) EXIT
  END DO
  CALL require(k > 100, 'I6: the recovering step equals the plain step')
  CALL CD_HermiteCable_Dyn_End(model)

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: in-line VIV on the cubic-Hermite cable'

CONTAINS

  SUBROUTINE build(m)
    !! The pinned taut line along x, pretensioned through its unstretched length, with the current
    !! along y (drag) and the normal added mass.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: m
    REAL(wp) :: seed(6*(NE + 1)), l0(NE), eav(NE), eiv(NE), wv(NE), rhoa(NE), hd(NE), hn(NE), ht(NE)
    INTEGER :: fixed(6), i
    seed = 0.0_wp
    DO i = 1, NE + 1
      seed(6*(i - 1) + 1) = SPAN*REAL(i - 1, wp)/REAL(NE, wp)
      seed(6*(i - 1) + 3) = -50.0_wp
      seed(6*(i - 1) + 4) = 1.0_wp
    END DO
    l0 = SPAN/REAL(NE, wp)/(1.0_wp + TENS/EA)
    eav = EA; eiv = EI; rhoa = MLIN
    wv = (MLIN - RHOW*0.25_wp*PI*D**2)*GACC
    fixed = [1, 2, 3, 6*NE + 1, 6*NE + 2, 6*NE + 3]
    CALL CD_HermiteCable_Dyn_Init(m, l0, eav, eiv, rhoa, wv, seed, fixed, -200.0_wp, 0.0_wp, 0.5_wp, es, em)
    CALL require(es == CD_HCDYN_OK, 'init: '//TRIM(em))
    hd = D; hn = 1.2_wp; ht = 0.0_wp
    CALL CD_HermiteCable_Dyn_Set_Drag(m, RHOW, hd, hn, ht, 0.0_wp, [0.0_wp, U, 0.0_wp], es, em)
    CALL require(es == CD_HCDYN_OK, 'drag: '//TRIM(em))
    hn = 1.0_wp
    CALL CD_HermiteCable_Dyn_Set_AddedMass(m, RHOW, hd, hn, ht, 0.0_wp, es, em)
    CALL require(es == CD_HCDYN_OK, 'added mass: '//TRIM(em))
  END SUBROUTINE build

  SUBROUTINE twin_without_il(m)
    !! The same line with the cross-flow VIV only.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: m
    CALL build(m)
    CALL CD_HermiteCable_Dyn_Set_VIV(m, cl, cf, df, 0.25_wp, es, em)
    CALL require(es == CD_HCDYN_OK, 'cross-flow VIV: '//TRIM(em))
  END SUBROUTINE twin_without_il

  SUBROUTINE advance(m, n, recovering)
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: m
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(IN) :: recovering
    INTEGER :: s
    DO s = 1, n
      IF (recovering) THEN
        CALL CD_HermiteCable_Dyn_Step_Recovering(m, DT, 50, 1.0e-8_wp, es, em)
      ELSE
        CALL CD_HermiteCable_Dyn_Step(m, DT, 50, 1.0e-8_wp, es, em)
      END IF
      IF (es /= CD_HCDYN_OK) THEN
        CALL require(.FALSE., 'step: '//TRIM(em)); RETURN
      END IF
    END DO
  END SUBROUTINE advance

  SUBROUTINE run_response(m, fy, fx, ay, ax)
    !! March 100 s and, over 50-100 s, the dominant frequencies of the mid-span cross-flow (z)
    !! and quarter-span in-line (y) motions (the peak of a direct Fourier scan of the detrended
    !! records) and the amplitudes, half the peak-to-peak excursions.
    TYPE(CD_HermiteCableDynType), INTENT(INOUT) :: m
    REAL(wp), INTENT(OUT) :: fy, fx, ay, ax
    INTEGER, PARAMETER :: NS = NINT(100.0_wp/DT), N0 = NINT(50.0_wp/DT)
    REAL(wp) :: zrec(NS - N0), yrec(NS - N0)
    INTEGER :: s
    DO s = 1, NS
      CALL advance(m, 1, .FALSE.)
      IF (es /= CD_HCDYN_OK) EXIT
      IF (s <= N0) CYCLE
      zrec(s - N0) = m%q(mid + 3)
      yrec(s - N0) = m%q(quarter + 2)
    END DO
    ay = 0.5_wp*(MAXVAL(zrec) - MINVAL(zrec))
    ax = 0.5_wp*(MAXVAL(yrec) - MINVAL(yrec))
    fy = peak_frequency(zrec)
    fx = peak_frequency(yrec)
  END SUBROUTINE run_response

  REAL(wp) FUNCTION peak_frequency(x) RESULT(fpk)
    !! The frequency of the largest Fourier amplitude of x - mean(x) on 0.05-2 Hz (0.002 Hz grid).
    REAL(wp), INTENT(IN) :: x(:)
    REAL(wp) :: xm(SIZE(x)), f, best, re, im, tt
    INTEGER :: i, j
    xm = x - SUM(x)/REAL(SIZE(x), wp)
    best = -1.0_wp
    fpk = 0.0_wp
    DO j = 0, 975
      f = 0.05_wp + 0.002_wp*REAL(j, wp)
      re = 0.0_wp; im = 0.0_wp
      DO i = 1, SIZE(x)
        tt = REAL(i, wp)*DT
        re = re + xm(i)*COS(2.0_wp*PI*f*tt)
        im = im + xm(i)*SIN(2.0_wp*PI*f*tt)
      END DO
      IF (re*re + im*im > best) THEN
        best = re*re + im*im
        fpk = f
      END IF
    END DO
  END FUNCTION peak_frequency

  SUBROUTINE check_element()
    !! I2: the element load at rest across the current, additivity, and the Jacobians against
    !! central differences on a curved, moving element.
    REAL(wp) :: qe(12), ve(12), fl(12), f_cf(12), f_il(12), djq(12, 12), djv(12, 12), fp(12), fm(12)
    REAL(wp) :: fdq(12, 12), fdv(12, 12), kn(2), kx(2), le, h, uu(3)
    INTEGER :: j
    le = 5.0_wp
    uu = [0.0_wp, U, 0.0_wp]
    kn = [0.0_wp, 0.0_wp]
    kx = [30.0_wp, 50.0_wp]
    ! straight along x at rest, unit tangents
    qe = 0.0_wp
    qe(4) = 1.0_wp; qe(7) = le; qe(10) = 1.0_wp
    ve = 0.0_wp
    CALL CD_HermiteCable_VIV_Lift_Element(qe, ve, le, uu, 0.0_wp, kn, fl, ErrStat=es, ErrMsg=em, kx_node=kx)
    CALL require(es == CD_HCDYN_OK, 'I2: element: '//TRIM(em))
    CALL require(ABS(fl(2) + fl(8) - U*U*le*0.5_wp*(kx(1) + kx(2))) <= 1.0e-12_wp*U*U*le*kx(2) .AND. &
                 ABS(fl(1) + fl(7)) + ABS(fl(3) + fl(9)) <= 1.0e-12_wp*U*U*le*kx(2), &
                 'I2: a straight element at rest carries U^2 l0 (k_x1 + k_x2)/2 along the flow')
    ! curved and moving, off the flow plane
    qe = [0.1_wp, -0.2_wp, -10.0_wp, 0.9_wp, 0.3_wp, 0.2_wp, 4.8_wp, 0.4_wp, -9.6_wp, 1.0_wp, -0.2_wp, 0.15_wp]
    ve = [0.02_wp, -0.05_wp, 0.08_wp, 0.01_wp, 0.02_wp, -0.03_wp, -0.04_wp, 0.06_wp, 0.03_wp, 0.0_wp, -0.01_wp, &
          0.02_wp]
    uu = [0.05_wp, U, -0.02_wp]
    kn = [20.0_wp, -15.0_wp]
    CALL CD_HermiteCable_VIV_Lift_Element(qe, ve, le, uu, 0.0_wp, kn, f_cf, ErrStat=es, ErrMsg=em)
    CALL CD_HermiteCable_VIV_Lift_Element(qe, ve, le, uu, 0.0_wp, [0.0_wp, 0.0_wp], f_il, ErrStat=es, ErrMsg=em, &
                                          kx_node=kx)
    CALL CD_HermiteCable_VIV_Lift_Element(qe, ve, le, uu, 0.0_wp, kn, fl, djq, djv, es, em, kx_node=kx)
    CALL require(MAXVAL(ABS(fl - f_cf - f_il)) <= 1.0e-12_wp*MAXVAL(ABS(fl)), &
                 'I2: the in-line load adds to the cross-flow one')
    h = 1.0e-6_wp
    DO j = 1, 12
      qe(j) = qe(j) + h
      CALL CD_HermiteCable_VIV_Lift_Element(qe, ve, le, uu, 0.0_wp, kn, fp, ErrStat=es, ErrMsg=em, kx_node=kx)
      qe(j) = qe(j) - 2.0_wp*h
      CALL CD_HermiteCable_VIV_Lift_Element(qe, ve, le, uu, 0.0_wp, kn, fm, ErrStat=es, ErrMsg=em, kx_node=kx)
      qe(j) = qe(j) + h
      fdq(:, j) = (fp - fm)/(2.0_wp*h)
      ve(j) = ve(j) + h
      CALL CD_HermiteCable_VIV_Lift_Element(qe, ve, le, uu, 0.0_wp, kn, fp, ErrStat=es, ErrMsg=em, kx_node=kx)
      ve(j) = ve(j) - 2.0_wp*h
      CALL CD_HermiteCable_VIV_Lift_Element(qe, ve, le, uu, 0.0_wp, kn, fm, ErrStat=es, ErrMsg=em, kx_node=kx)
      ve(j) = ve(j) + h
      fdv(:, j) = (fp - fm)/(2.0_wp*h)
    END DO
    WRITE (*, '(A,ES10.2,A,ES10.2)') '  I2: Jacobian errors, position ', MAXVAL(ABS(djq - fdq))/MAXVAL(ABS(fdq)), &
      ', velocity ', MAXVAL(ABS(djv - fdv))/MAXVAL(ABS(fdv))
    CALL require(MAXVAL(ABS(djq - fdq)) <= 1.0e-6_wp*MAXVAL(ABS(fdq)), 'I2: position Jacobian matches differences')
    CALL require(MAXVAL(ABS(djv - fdv)) <= 1.0e-6_wp*MAXVAL(ABS(fdv)), 'I2: velocity Jacobian matches differences')
  END SUBROUTINE check_element

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. cond) THEN
      WRITE (*, '(A,A)') 'MISMATCH: ', label
      nfail = nfail + 1
    END IF
  END SUBROUTINE require

END PROGRAM test_hermite_viv_inline
