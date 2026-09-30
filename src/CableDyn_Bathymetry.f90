! File: src/CableDyn_Bathymetry.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_Bathymetry
  !! Structured 3D bathymetry service for production seabed contact and local
  !! water-depth queries. The grid stores positive water depth below still-water
  !! level at monotone x/y coordinates; the seabed elevation is z_floor = -depth.
  !! Queries clamp to the nearest grid edge outside the supplied domain and use
  !! bilinear interpolation inside each cell. Contact is a frictionless penalty force
  !! along the local upward surface normal at the normal penetration; its tangent holds
  !! the surface gradient fixed over the node (the bilinear curvature is omitted).
  USE CableDyn_Precision, ONLY: wp, CD_ZERO, CD_All_Finite, CD_Is_Finite
  USE CableDyn_SeabedContact, ONLY: CD_Seabed_Normal_Contact
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CD_BathymetryType
  PUBLIC :: CD_Init_Bathymetry
  PUBLIC :: CD_End_Bathymetry
  PUBLIC :: CD_Bathymetry_Is_Initialized
  PUBLIC :: CD_Bathymetry_Depth
  PUBLIC :: CD_Bathymetry_Floor
  PUBLIC :: CD_Bathymetry_Floor_Gradient
  PUBLIC :: CD_Bathymetry_Seabed_Load
  INTEGER, PARAMETER, PUBLIC :: CD_BATHY_OK = 0, CD_BATHY_BADINPUT = 1, CD_BATHY_ALLOCFAIL = 2

  TYPE :: CD_BathymetryType
    INTEGER :: nx = 0, ny = 0
    REAL(wp), ALLOCATABLE :: x(:)
    REAL(wp), ALLOCATABLE :: y(:)
    REAL(wp), ALLOCATABLE :: depth(:, :)
    REAL(wp) :: average_depth = CD_ZERO
    REAL(wp) :: minimum_depth = CD_ZERO
  END TYPE CD_BathymetryType

CONTAINS

  SUBROUTINE CD_Init_Bathymetry(bathy, x, y, depth, ErrStat, ErrMsg)
    !! Initialize a validated structured bathymetry grid.
    TYPE(CD_BathymetryType), INTENT(INOUT) :: bathy
    REAL(wp), INTENT(IN) :: x(:)
    REAL(wp), INTENT(IN) :: y(:)
    REAL(wp), INTENT(IN) :: depth(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    TYPE(CD_BathymetryType) :: next_bathy
    INTEGER :: i, istat

    ErrStat = CD_BATHY_OK
    ErrMsg = ''

    IF (SIZE(x) < 2 .OR. SIZE(y) < 2) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry x/y grids need at least two points')
      RETURN
    END IF
    IF (SIZE(depth, 1) /= SIZE(x) .OR. SIZE(depth, 2) /= SIZE(y)) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry depth must have shape (nx, ny)')
      RETURN
    END IF
    IF (.NOT. (CD_All_Finite(x) .AND. CD_All_Finite(y) .AND. CD_All_Finite(depth))) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry coordinates and depths must be finite')
      RETURN
    END IF
    DO i = 2, SIZE(x)
      IF (x(i) <= x(i - 1)) THEN
        CALL fail(ErrStat, ErrMsg, 'bathymetry x grid must be strictly increasing')
        RETURN
      END IF
    END DO
    DO i = 2, SIZE(y)
      IF (y(i) <= y(i - 1)) THEN
        CALL fail(ErrStat, ErrMsg, 'bathymetry y grid must be strictly increasing')
        RETURN
      END IF
    END DO
    IF (ANY(depth <= CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry depth entries must be positive')
      RETURN
    END IF

    next_bathy%nx = SIZE(x)
    next_bathy%ny = SIZE(y)
    ALLOCATE (next_bathy%x(next_bathy%nx), next_bathy%y(next_bathy%ny), &
              next_bathy%depth(next_bathy%nx, next_bathy%ny), STAT=istat)
    IF (istat /= 0) THEN
      ErrStat = CD_BATHY_ALLOCFAIL
      ErrMsg = 'CableDyn_Bathymetry: bathymetry grid allocation failed'
      RETURN
    END IF
    next_bathy%x = x
    next_bathy%y = y
    next_bathy%depth = depth
    next_bathy%average_depth = SUM(depth)/REAL(SIZE(depth), wp)
    next_bathy%minimum_depth = MINVAL(depth)

    CALL CD_End_Bathymetry(bathy)
    bathy%nx = next_bathy%nx
    bathy%ny = next_bathy%ny
    bathy%average_depth = next_bathy%average_depth
    bathy%minimum_depth = next_bathy%minimum_depth
    CALL MOVE_ALLOC(next_bathy%x, bathy%x)
    CALL MOVE_ALLOC(next_bathy%y, bathy%y)
    CALL MOVE_ALLOC(next_bathy%depth, bathy%depth)
  END SUBROUTINE CD_Init_Bathymetry

  SUBROUTINE CD_End_Bathymetry(bathy)
    !! Release grid storage and reset metadata.
    TYPE(CD_BathymetryType), INTENT(INOUT) :: bathy

    IF (ALLOCATED(bathy%x)) DEALLOCATE (bathy%x)
    IF (ALLOCATED(bathy%y)) DEALLOCATE (bathy%y)
    IF (ALLOCATED(bathy%depth)) DEALLOCATE (bathy%depth)
    bathy%nx = 0
    bathy%ny = 0
    bathy%average_depth = CD_ZERO
    bathy%minimum_depth = CD_ZERO
  END SUBROUTINE CD_End_Bathymetry

  LOGICAL FUNCTION CD_Bathymetry_Is_Initialized(bathy) RESULT(ok)
    !! True when the bathymetry object owns a complete validated grid.
    TYPE(CD_BathymetryType), INTENT(IN) :: bathy

    ok = bathy%nx >= 2 .AND. bathy%ny >= 2 .AND. ALLOCATED(bathy%x) .AND. &
         ALLOCATED(bathy%y) .AND. ALLOCATED(bathy%depth)
  END FUNCTION CD_Bathymetry_Is_Initialized

  SUBROUTINE CD_Bathymetry_Depth(bathy, xq, yq, depth, ErrStat, ErrMsg)
    !! Return positive water depth at a horizontal location.
    TYPE(CD_BathymetryType), INTENT(IN) :: bathy
    REAL(wp), INTENT(IN) :: xq, yq
    REAL(wp), INTENT(OUT) :: depth
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp) :: floor, dfdx, dfdy

    CALL floor_eval(bathy, xq, yq, floor, dfdx, dfdy, ErrStat, ErrMsg)
    IF (ErrStat /= CD_BATHY_OK) THEN
      depth = CD_ZERO
      RETURN
    END IF
    depth = -floor
  END SUBROUTINE CD_Bathymetry_Depth

  SUBROUTINE CD_Bathymetry_Floor(bathy, xq, yq, z_floor, ErrStat, ErrMsg)
    !! Return seabed elevation z_floor at a horizontal location.
    TYPE(CD_BathymetryType), INTENT(IN) :: bathy
    REAL(wp), INTENT(IN) :: xq, yq
    REAL(wp), INTENT(OUT) :: z_floor
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    REAL(wp) :: dfdx, dfdy

    CALL floor_eval(bathy, xq, yq, z_floor, dfdx, dfdy, ErrStat, ErrMsg)
  END SUBROUTINE CD_Bathymetry_Floor

  SUBROUTINE CD_Bathymetry_Floor_Gradient(bathy, xq, yq, z_floor, dfdx, dfdy, ErrStat, ErrMsg)
    !! Return seabed elevation and horizontal gradients dz_floor/dx, dz_floor/dy.
    TYPE(CD_BathymetryType), INTENT(IN) :: bathy
    REAL(wp), INTENT(IN) :: xq, yq
    REAL(wp), INTENT(OUT) :: z_floor, dfdx, dfdy
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    CALL floor_eval(bathy, xq, yq, z_floor, dfdx, dfdy, ErrStat, ErrMsg)
  END SUBROUTINE CD_Bathymetry_Floor_Gradient

  SUBROUTINE CD_Bathymetry_Seabed_Load(bathy, nodes, k_n, f, k_resid, ErrStat, ErrMsg)
    !! Variable-elevation penalty contact against the bathymetry surface.
    !! For node i, g_i = z_floor(x_i,y_i) - z_i is the vertical penetration. The contact
    !! is frictionless: the force acts along the local upward surface normal with the
    !! shared C1 touchdown law evaluated at the normal penetration (CD_Seabed_Normal_Contact;
    !! exactly the vertical law on a level patch). The residual tangent is the node's
    !! -d(force)/dq block, F' n n^T.
    TYPE(CD_BathymetryType), INTENT(IN) :: bathy
    REAL(wp), INTENT(IN) :: nodes(:, :)
    REAL(wp), INTENT(IN) :: k_n(:)
    REAL(wp), INTENT(OUT) :: f(:)
    REAL(wp), INTENT(OUT) :: k_resid(:, :)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: n_nodes, i, ix, iz
    REAL(wp) :: z_floor, dfdx, dfdy, gap, fvec(3), jac(3, 3)

    f = CD_ZERO
    k_resid = CD_ZERO
    ErrStat = CD_BATHY_OK
    ErrMsg = ''
    IF (.NOT. CD_Bathymetry_Is_Initialized(bathy)) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry object is not initialized')
      RETURN
    END IF
    IF (SIZE(nodes, 1) /= 3) THEN
      CALL fail(ErrStat, ErrMsg, 'nodes must have shape (3, n_nodes)')
      RETURN
    END IF
    n_nodes = SIZE(nodes, 2)
    IF (SIZE(k_n) /= n_nodes) THEN
      CALL fail(ErrStat, ErrMsg, 'k_n must have shape (n_nodes)')
      RETURN
    END IF
    IF (SIZE(f) /= 3*n_nodes .OR. SIZE(k_resid, 1) /= 3*n_nodes .OR. SIZE(k_resid, 2) /= 3*n_nodes) THEN
      CALL fail(ErrStat, ErrMsg, 'f / k_resid shapes inconsistent with nodes')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(nodes)) THEN
      CALL fail(ErrStat, ErrMsg, 'node coordinates must be finite')
      RETURN
    END IF
    IF (.NOT. CD_All_Finite(k_n) .OR. ANY(k_n <= CD_ZERO)) THEN
      CALL fail(ErrStat, ErrMsg, 'k_n must be finite and positive')
      RETURN
    END IF

    DO i = 1, n_nodes
      CALL floor_eval(bathy, nodes(1, i), nodes(2, i), z_floor, dfdx, dfdy, ErrStat, ErrMsg)
      IF (ErrStat /= CD_BATHY_OK) RETURN
      gap = z_floor - nodes(3, i)
      ix = 3*i - 2
      iz = ix + 2
      CALL CD_Seabed_Normal_Contact(gap, dfdx, dfdy, k_n(i), fvec, jac)
      f(ix:iz) = fvec
      k_resid(ix:iz, ix:iz) = jac
    END DO
  END SUBROUTINE CD_Bathymetry_Seabed_Load

  SUBROUTINE floor_eval(bathy, xq, yq, z_floor, dfdx, dfdy, ErrStat, ErrMsg)
    TYPE(CD_BathymetryType), INTENT(IN) :: bathy
    REAL(wp), INTENT(IN) :: xq, yq
    REAL(wp), INTENT(OUT) :: z_floor, dfdx, dfdy
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg

    INTEGER :: ix, iy
    REAL(wp) :: tx, ty
    REAL(wp) :: f00, f10, f01, f11, dx, dy
    LOGICAL :: clamp_x, clamp_y

    z_floor = CD_ZERO
    dfdx = CD_ZERO
    dfdy = CD_ZERO
    ErrStat = CD_BATHY_OK
    ErrMsg = ''
    IF (.NOT. CD_Bathymetry_Is_Initialized(bathy)) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry object is not initialized')
      RETURN
    END IF
    IF (.NOT. (CD_Is_Finite(xq) .AND. CD_Is_Finite(yq))) THEN
      CALL fail(ErrStat, ErrMsg, 'bathymetry query coordinates must be finite')
      RETURN
    END IF

    CALL locate_axis(bathy%x, xq, ix, tx, clamp_x)
    CALL locate_axis(bathy%y, yq, iy, ty, clamp_y)
    dx = bathy%x(ix + 1) - bathy%x(ix)
    dy = bathy%y(iy + 1) - bathy%y(iy)
    f00 = -bathy%depth(ix, iy)
    f10 = -bathy%depth(ix + 1, iy)
    f01 = -bathy%depth(ix, iy + 1)
    f11 = -bathy%depth(ix + 1, iy + 1)

    z_floor = (1.0_wp - tx)*(1.0_wp - ty)*f00 + tx*(1.0_wp - ty)*f10 + &
              (1.0_wp - tx)*ty*f01 + tx*ty*f11
    IF (.NOT. clamp_x) dfdx = ((1.0_wp - ty)*(f10 - f00) + ty*(f11 - f01))/dx
    IF (.NOT. clamp_y) dfdy = ((1.0_wp - tx)*(f01 - f00) + tx*(f11 - f10))/dy
  END SUBROUTINE floor_eval

  SUBROUTINE locate_axis(axis, q, i0, t, clamped)
    REAL(wp), INTENT(IN) :: axis(:), q
    INTEGER, INTENT(OUT) :: i0
    REAL(wp), INTENT(OUT) :: t
    LOGICAL, INTENT(OUT) :: clamped

    INTEGER :: i, n

    n = SIZE(axis)
    clamped = .FALSE.
    IF (q <= axis(1)) THEN
      i0 = 1
      t = CD_ZERO
      clamped = q < axis(1)
      RETURN
    END IF
    IF (q >= axis(n)) THEN
      i0 = n - 1
      t = 1.0_wp
      clamped = q > axis(n)
      RETURN
    END IF
    DO i = 1, n - 1
      IF (q >= axis(i) .AND. q <= axis(i + 1)) THEN
        i0 = i
        t = (q - axis(i))/(axis(i + 1) - axis(i))
        RETURN
      END IF
    END DO

    i0 = n - 1
    t = 1.0_wp
    clamped = .TRUE.
  END SUBROUTINE locate_axis

  SUBROUTINE fail(ErrStat, ErrMsg, msg)
    INTEGER, INTENT(OUT) :: ErrStat
    CHARACTER(*), INTENT(OUT) :: ErrMsg
    CHARACTER(*), INTENT(IN) :: msg

    ErrStat = CD_BATHY_BADINPUT
    ErrMsg = 'CableDyn_Bathymetry: '//msg
  END SUBROUTINE fail

END MODULE CableDyn_Bathymetry
