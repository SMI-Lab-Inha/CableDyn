! File: tests/test_bathymetry.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_bathymetry
  !! Unit tests for the structured bathymetry service and variable-elevation
  !! seabed contact tangent. These tests are self-contained analytical checks:
  !! bilinear interpolation on a plane, edge clamping, and fail-closed validation.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_Bathymetry, ONLY: CD_BathymetryType, CD_Init_Bathymetry, CD_End_Bathymetry, &
                                 CD_Bathymetry_Depth, CD_Bathymetry_Floor, &
                                 CD_Bathymetry_Floor_Gradient, CD_Bathymetry_Seabed_Load, &
                                 CD_Bathymetry_Is_Initialized, CD_BATHY_OK, CD_BATHY_BADINPUT
  USE CableDyn_SeabedContact, ONLY: CD_SEABED_CONTACT_BLEND
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
  IMPLICIT NONE

  INTEGER :: nfail
  nfail = 0

  CALL case_bilinear_depth_and_floor()
  CALL case_variable_seabed_contact_tangent()
  CALL case_variable_seabed_tangent_fd()
  CALL case_fail_closed_inputs()

  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAIL: ', nfail, ' assertion(s) failed'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: CableDyn_Bathymetry'

CONTAINS

  INCLUDE 'nan_max_abs.inc'

  SUBROUTINE expect(got, want, label)
    REAL(wp), INTENT(IN) :: got, want
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), PARAMETER :: atol = 1.0e-10_wp, rtol = 1.0e-12_wp
    IF (.NOT. (ABS(got - want) <= atol + rtol*ABS(want))) THEN
      WRITE (*, '(A,A,A,ES23.15,A,ES23.15)') 'MISMATCH [', label, ']: got ', got, ' want ', want
      nfail = nfail + 1
    END IF
  END SUBROUTINE expect

  SUBROUTINE expect_es(es, want, label)
    INTEGER, INTENT(IN) :: es, want
    CHARACTER(*), INTENT(IN) :: label
    IF (es /= want) THEN
      WRITE (*, '(A,A,A,I0,A,I0)') 'MISMATCH [', label, ']: ErrStat ', es, ' want ', want
      nfail = nfail + 1
    END IF
  END SUBROUTINE expect_es

  SUBROUTINE expect_relaxed(got, want, label)
    REAL(wp), INTENT(IN) :: got, want
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), PARAMETER :: atol = 1.0e-5_wp, rtol = 1.0e-10_wp
    IF (.NOT. (ABS(got - want) <= atol + rtol*ABS(want))) THEN
      WRITE (*, '(A,A,A,ES23.15,A,ES23.15)') 'MISMATCH [', label, ']: got ', got, ' want ', want
      nfail = nfail + 1
    END IF
  END SUBROUTINE expect_relaxed

  SUBROUTINE require(cond, label)
    LOGICAL, INTENT(IN) :: cond
    CHARACTER(*), INTENT(IN) :: label
    IF (.NOT. cond) THEN
      WRITE (*, '(A,A)') 'MISMATCH [', label//']'
      nfail = nfail + 1
    END IF
  END SUBROUTINE require

  SUBROUTINE init_plane_bathy(bathy, es, em)
    TYPE(CD_BathymetryType), INTENT(OUT) :: bathy
    INTEGER, INTENT(OUT) :: es
    CHARACTER(*), INTENT(OUT) :: em

    REAL(wp) :: x(3), y(3), depth(3, 3)
    INTEGER :: i, j

    x = [0.0_wp, 10.0_wp, 20.0_wp]
    y = [0.0_wp, 5.0_wp, 15.0_wp]
    DO j = 1, SIZE(y)
      DO i = 1, SIZE(x)
        depth(i, j) = 100.0_wp + 0.2_wp*x(i) + 0.4_wp*y(j)
      END DO
    END DO
    CALL CD_Init_Bathymetry(bathy, x, y, depth, es, em)
  END SUBROUTINE init_plane_bathy

  SUBROUTINE case_bilinear_depth_and_floor()
    TYPE(CD_BathymetryType) :: bathy
    REAL(wp) :: depth, floor, dfdx, dfdy
    INTEGER :: es
    CHARACTER(160) :: em

    CALL init_plane_bathy(bathy, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'init-plane')
    CALL require(CD_Bathymetry_Is_Initialized(bathy), 'is-initialized')
    CALL expect(bathy%average_depth, 104.66666666666666667_wp, 'average-depth')
    CALL expect(bathy%minimum_depth, 100.0_wp, 'minimum-depth')

    CALL CD_Bathymetry_Depth(bathy, 7.5_wp, 2.5_wp, depth, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'depth-interior-es')
    CALL expect(depth, 102.5_wp, 'depth-interior')

    CALL CD_Bathymetry_Floor_Gradient(bathy, 7.5_wp, 2.5_wp, floor, dfdx, dfdy, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'floor-gradient-es')
    CALL expect(floor, -102.5_wp, 'floor-interior')
    CALL expect(dfdx, -0.2_wp, 'floor-dfdx')
    CALL expect(dfdy, -0.4_wp, 'floor-dfdy')

    CALL CD_Bathymetry_Floor(bathy, -2.0_wp, 20.0_wp, floor, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'floor-clamped-es')
    CALL expect(floor, -106.0_wp, 'floor-clamped-edge')
    CALL CD_Bathymetry_Floor_Gradient(bathy, -2.0_wp, 7.5_wp, floor, dfdx, dfdy, es, em)
    CALL expect(dfdx, 0.0_wp, 'floor-clamped-x-gradient')
    CALL expect(dfdy, -0.4_wp, 'floor-clamped-y-gradient')

    CALL CD_End_Bathymetry(bathy)
    CALL require(.NOT. CD_Bathymetry_Is_Initialized(bathy), 'end-resets')
  END SUBROUTINE case_bilinear_depth_and_floor

  SUBROUTINE case_variable_seabed_contact_tangent()
    TYPE(CD_BathymetryType) :: bathy
    REAL(wp) :: nodes(3, 3), kn(3), f(9), kres(9, 9)
    INTEGER :: es
    CHARACTER(160) :: em

    CALL init_plane_bathy(bathy, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'contact:init')

    nodes(:, 1) = [7.5_wp, 2.5_wp, -103.0_wp]   ! 0.5 m penetration
    nodes(:, 2) = [7.5_wp, 2.5_wp, -102.5_wp]   ! exactly on floor
    nodes(:, 3) = [7.5_wp, 2.5_wp, -102.0_wp]   ! above floor
    kn = [1000.0_wp, 2000.0_wp, 3000.0_wp]

    ! Frictionless contact along the surface normal. The test plane has gradient
    ! (dz/dx, dz/dy) = (-0.2, -0.4), so the upward unit normal is n = (0.2, 0.4, 1)/s,
    ! s = sqrt(1.2); a vertical penetration g is a normal penetration g/s, the force is
    ! k (g/s) n and the residual tangent k n n^T.
    CALL CD_Bathymetry_Seabed_Load(bathy, nodes, kn, f, kres, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'contact:load-es')
    CALL expect(f(3), 500.0_wp/1.2_wp, 'contact:penetrating-force')
    CALL expect(f(1), 0.2_wp*500.0_wp/1.2_wp, 'contact:penetrating-force-x')
    CALL expect(f(2), 0.4_wp*500.0_wp/1.2_wp, 'contact:penetrating-force-y')
    CALL expect(f(6), 500.0_wp*CD_SEABED_CONTACT_BLEND/SQRT(1.2_wp), 'contact:on-floor-force')
    CALL expect(f(9), 0.0_wp, 'contact:above-force')
    CALL expect(kres(3, 1), 200.0_wp/1.2_wp, 'contact:node1-dx-tangent')
    CALL expect(kres(3, 2), 400.0_wp/1.2_wp, 'contact:node1-dy-tangent')
    CALL expect(kres(3, 3), 1000.0_wp/1.2_wp, 'contact:node1-dz-tangent')
    CALL expect(kres(1, 1), 40.0_wp/1.2_wp, 'contact:node1-xx-tangent')
    CALL expect(kres(1, 2), 80.0_wp/1.2_wp, 'contact:node1-xy-tangent')
    CALL expect(kres(1, 3), 200.0_wp/1.2_wp, 'contact:node1-symmetric-tangent')
    CALL expect_relaxed(kres(6, 4), 200.0_wp/1.2_wp, 'contact:node2-dx-blend')
    CALL expect_relaxed(kres(6, 5), 400.0_wp/1.2_wp, 'contact:node2-dy-blend')
    CALL expect_relaxed(kres(6, 6), 1000.0_wp/1.2_wp, 'contact:node2-dz-blend')
    CALL expect(kres(9, 7), 0.0_wp, 'contact:node3-inactive-x')
    CALL expect(kres(9, 9), 0.0_wp, 'contact:node3-inactive-z')
    CALL expect(nan_max_abs(f(7:9)), 0.0_wp, 'contact:no-force-above-floor')

    ! half-way into the touchdown blend along the NORMAL (vertical offset s*blend/2)
    nodes(:, 1) = [7.5_wp, 2.5_wp, -102.5_wp - 0.5_wp*SQRT(1.2_wp)*CD_SEABED_CONTACT_BLEND]
    CALL CD_Bathymetry_Seabed_Load(bathy, nodes(:, 1:1), kn(1:1), f(1:3), kres(1:3, 1:3), es, em)
    CALL expect_es(es, CD_BATHY_OK, 'contact:blend-load-es')
    CALL expect_relaxed(f(3), 562.5_wp*CD_SEABED_CONTACT_BLEND/SQRT(1.2_wp), 'contact:blend-force')
    CALL expect_relaxed(kres(3, 1), 150.0_wp/1.2_wp, 'contact:blend-dx-tangent')
    CALL expect_relaxed(kres(3, 2), 300.0_wp/1.2_wp, 'contact:blend-dy-tangent')
    CALL expect_relaxed(kres(3, 3), 750.0_wp/1.2_wp, 'contact:blend-dz-tangent')

    CALL CD_End_Bathymetry(bathy)
  END SUBROUTINE case_variable_seabed_contact_tangent

  SUBROUTINE case_variable_seabed_tangent_fd()
    !! Independent FD gate for the variable-floor residual tangent:
    !! k_resid(:,j) = -d f_contact / d q_j on an active penetrating node.
    TYPE(CD_BathymetryType) :: bathy
    REAL(wp) :: nodes(3, 1), kn(1), f0(3), fp(3), fm(3), kres(3, 3), kbase(3, 3), fd_col(3), trial(3, 1)
    REAL(wp) :: h
    INTEGER :: es, j
    CHARACTER(160) :: em

    CALL init_plane_bathy(bathy, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'fd:init')
    nodes(:, 1) = [7.5_wp, 2.5_wp, -103.0_wp]
    kn = [1000.0_wp]
    h = 1.0e-6_wp

    CALL CD_Bathymetry_Seabed_Load(bathy, nodes, kn, f0, kres, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'fd:base-load')
    kbase = kres
    DO j = 1, 3
      trial = nodes
      trial(j, 1) = trial(j, 1) + h
      CALL CD_Bathymetry_Seabed_Load(bathy, trial, kn, fp, kres, es, em)
      CALL expect_es(es, CD_BATHY_OK, 'fd:plus-load')
      trial = nodes
      trial(j, 1) = trial(j, 1) - h
      CALL CD_Bathymetry_Seabed_Load(bathy, trial, kn, fm, kres, es, em)
      CALL expect_es(es, CD_BATHY_OK, 'fd:minus-load')
      fd_col = -(fp - fm)/(2.0_wp*h)
      CALL expect_vec_fd(kbase(:, j), fd_col, 'fd:variable-floor-kcol')
    END DO
    CALL expect(f0(3), 500.0_wp/1.2_wp, 'fd:base-force')
    CALL CD_End_Bathymetry(bathy)
  END SUBROUTINE case_variable_seabed_tangent_fd

  SUBROUTINE case_fail_closed_inputs()
    TYPE(CD_BathymetryType) :: bathy
    REAL(wp) :: x(3), y(2), depth(3, 2), d, f(6), k(6, 6), nodes(3, 2)
    INTEGER :: es
    CHARACTER(160) :: em

    x = [0.0_wp, 10.0_wp, 20.0_wp]
    y = [0.0_wp, 5.0_wp]
    depth = RESHAPE([100.0_wp, 102.0_wp, 104.0_wp, 101.0_wp, 103.0_wp, 105.0_wp], [3, 2])

    CALL CD_Bathymetry_Depth(bathy, 0.0_wp, 0.0_wp, d, es, em)
    CALL expect_es(es, CD_BATHY_BADINPUT, 'uninitialized-depth')

    CALL CD_Init_Bathymetry(bathy, [0.0_wp, 0.0_wp, 10.0_wp], y, depth, es, em)
    CALL expect_es(es, CD_BATHY_BADINPUT, 'duplicate-x')

    CALL CD_Init_Bathymetry(bathy, x, y, -depth, es, em)
    CALL expect_es(es, CD_BATHY_BADINPUT, 'negative-depth')

    depth(2, 1) = IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN)
    CALL CD_Init_Bathymetry(bathy, x, y, depth, es, em)
    CALL expect_es(es, CD_BATHY_BADINPUT, 'nan-depth')

    depth = RESHAPE([100.0_wp, 102.0_wp, 104.0_wp, 101.0_wp, 103.0_wp, 105.0_wp], [3, 2])
    CALL CD_Init_Bathymetry(bathy, x, y, depth, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'valid-after-failures')
    CALL CD_Bathymetry_Depth(bathy, 10.0_wp, 5.0_wp, d, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'valid-before-bad-reinit-depth-es')
    CALL expect(d, 103.0_wp, 'valid-before-bad-reinit-depth')
    CALL CD_Init_Bathymetry(bathy, [0.0_wp, 0.0_wp, 10.0_wp], y, depth, es, em)
    CALL expect_es(es, CD_BATHY_BADINPUT, 'bad-reinit-preserves-grid-es')
    CALL require(CD_Bathymetry_Is_Initialized(bathy), 'bad-reinit-preserves-grid-init')
    CALL CD_Bathymetry_Depth(bathy, 10.0_wp, 5.0_wp, d, es, em)
    CALL expect_es(es, CD_BATHY_OK, 'bad-reinit-preserves-depth-es')
    CALL expect(d, 103.0_wp, 'bad-reinit-preserves-depth')

    nodes(:, 1) = [0.0_wp, 0.0_wp, -101.0_wp]
    nodes(:, 2) = [1.0_wp, 0.0_wp, -101.0_wp]
    CALL CD_Bathymetry_Seabed_Load(bathy, nodes, [100.0_wp], f, k, es, em)
    CALL expect_es(es, CD_BATHY_BADINPUT, 'bad-kn-size')
    CALL CD_Bathymetry_Seabed_Load(bathy, nodes, [100.0_wp, -1.0_wp], f, k, es, em)
    CALL expect_es(es, CD_BATHY_BADINPUT, 'bad-kn-value')

    nodes(1, 2) = IEEE_VALUE(1.0_wp, IEEE_QUIET_NAN)
    CALL CD_Bathymetry_Seabed_Load(bathy, nodes, [100.0_wp, 100.0_wp], f, k, es, em)
    CALL expect_es(es, CD_BATHY_BADINPUT, 'nan-node')

    CALL CD_End_Bathymetry(bathy)
  END SUBROUTINE case_fail_closed_inputs

  SUBROUTINE expect_vec_fd(got, want, label)
    REAL(wp), INTENT(IN) :: got(:), want(:)
    CHARACTER(*), INTENT(IN) :: label
    REAL(wp), PARAMETER :: atol = 2.0e-5_wp, rtol = 1.0e-9_wp
    INTEGER :: i

    DO i = 1, SIZE(got)
      IF (.NOT. (ABS(got(i) - want(i)) <= atol + rtol*ABS(want(i)))) THEN
        WRITE (*, '(A,A,A,ES23.15,A,ES23.15)') 'MISMATCH [', label, ']: got ', got(i), ' want ', want(i)
        nfail = nfail + 1
      END IF
    END DO
  END SUBROUTINE expect_vec_fd

END PROGRAM test_bathymetry
