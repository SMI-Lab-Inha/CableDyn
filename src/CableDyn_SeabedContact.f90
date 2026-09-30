! File: src/CableDyn_SeabedContact.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
MODULE CableDyn_SeabedContact
  !! Shared normal-contact law for flat and structured-bathymetry seabed contact.
  !! This module implements the C1 touchdown smoothing required by the production
  !! mooring program: the normal force and residual tangent are continuous through
  !! first contact, while the law is exactly linear beyond a small blend length.
  !! It also holds the seabed friction laws shared by the EI = 0 and finite-EI lines:
  !! the static friction spring and the stick-slip spring of the dynamics.
  USE CableDyn_Precision, ONLY: wp, CD_ZERO, CD_ONE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: CD_SEABED_CONTACT_BLEND
  PUBLIC :: CD_Seabed_Normal_Law
  PUBLIC :: CD_Seabed_Normal_Contact
  PUBLIC :: CD_Seabed_Contact_Gate
  PUBLIC :: CD_Seabed_Friction_Spring
  PUBLIC :: CD_Seabed_Friction_Stick_Slip
  PUBLIC :: CD_Seabed_Friction_Aniso, CD_Seabed_Friction_Mu_Dir
  PUBLIC :: CD_FRICTION_SPRING, CD_FRICTION_STICK_SLIP

  INTEGER, PARAMETER :: CD_FRICTION_SPRING = 1, CD_FRICTION_STICK_SLIP = 2
    !! Law selectors of CD_Seabed_Friction_Aniso: the static spring or the stick-slip spring.

  REAL(wp), PARAMETER :: CD_SEABED_CONTACT_BLEND = 1.0e-6_wp
    !! Penetration length [m] over which the normal penalty ramps C1 from zero
    !! tangent to the full penalty stiffness.

CONTAINS

  PURE SUBROUTINE CD_Seabed_Normal_Law(gap, stiffness, force, dforce_dgap)
    !! C1 normal contact law for signed penetration ``gap = z_floor - z``.
    !! For gap <= -blend, force and tangent are zero. For -blend < gap < blend,
    !! force = k*(gap + blend)^2/(4*blend); for gap >= blend, force = k*gap.
    !! The returned tangent is d(force)/d(gap), which is also the residual
    !! tangent contribution for a flat seabed z DOF.
    REAL(wp), INTENT(IN) :: gap
    REAL(wp), INTENT(IN) :: stiffness
    REAL(wp), INTENT(OUT) :: force
    REAL(wp), INTENT(OUT) :: dforce_dgap

    IF (gap <= -CD_SEABED_CONTACT_BLEND) THEN
      force = CD_ZERO
      dforce_dgap = CD_ZERO
    ELSE IF (gap < CD_SEABED_CONTACT_BLEND) THEN
      force = stiffness*(gap + CD_SEABED_CONTACT_BLEND)**2/(4.0_wp*CD_SEABED_CONTACT_BLEND)
      dforce_dgap = stiffness*(gap + CD_SEABED_CONTACT_BLEND)/(2.0_wp*CD_SEABED_CONTACT_BLEND)
    ELSE
      force = stiffness*gap
      dforce_dgap = stiffness
    END IF
  END SUBROUTINE CD_Seabed_Normal_Law

  PURE SUBROUTINE CD_Seabed_Normal_Contact(gap, dfdx, dfdy, stiffness, force, jac)
    !! Frictionless penalty contact with a locally planar seabed z = f(x, y) of gradient
    !! (dfdx, dfdy). ``gap = f(x, y) - z`` is the vertical penetration; the penetration
    !! normal to the surface is gap/s, s = sqrt(1 + dfdx^2 + dfdy^2), and the contact force
    !! acts along the upward unit normal n = (-dfdx, -dfdy, 1)/s with the magnitude of
    !! CD_Seabed_Normal_Law at that normal penetration. On a level floor this is exactly
    !! the vertical law. ``jac`` is the residual tangent -d(force)/dq = F'(gap/s) n n^T
    !! (the surface gradient is held fixed over the node, i.e. the bathymetry curvature
    !! term is omitted).
    REAL(wp), INTENT(IN) :: gap, dfdx, dfdy, stiffness
    REAL(wp), INTENT(OUT) :: force(3), jac(3, 3)
    REAL(wp) :: s, normal(3), magnitude, dmag
    INTEGER :: a, b

    s = SQRT(CD_ONE + dfdx*dfdx + dfdy*dfdy)
    normal = [-dfdx, -dfdy, CD_ONE]/s
    CALL CD_Seabed_Normal_Law(gap/s, stiffness, magnitude, dmag)
    force = magnitude*normal
    DO b = 1, 3
      DO a = 1, 3
        jac(a, b) = dmag*normal(a)*normal(b)
      END DO
    END DO
  END SUBROUTINE CD_Seabed_Normal_Contact

  PURE SUBROUTINE CD_Seabed_Contact_Gate(gap, gate, dgate_dgap)
    !! C0 contact activation and piecewise-constant derivative consistent with
    !! CD_Seabed_Normal_Law. The gate is zero above the blend, one below it, and
    !! ramps linearly through the touchdown blend. It is intended for dissipative
    !! terms whose force should be continuous through first contact.
    REAL(wp), INTENT(IN) :: gap
    REAL(wp), INTENT(OUT) :: gate
    REAL(wp), INTENT(OUT) :: dgate_dgap

    IF (gap <= -CD_SEABED_CONTACT_BLEND) THEN
      gate = CD_ZERO
      dgate_dgap = CD_ZERO
    ELSE IF (gap < CD_SEABED_CONTACT_BLEND) THEN
      gate = (gap + CD_SEABED_CONTACT_BLEND)/(2.0_wp*CD_SEABED_CONTACT_BLEND)
      dgate_dgap = CD_ONE/(2.0_wp*CD_SEABED_CONTACT_BLEND)
    ELSE
      gate = CD_ONE
      dgate_dgap = CD_ZERO
    END IF
  END SUBROUTINE CD_Seabed_Contact_Gate

  PURE SUBROUTINE CD_Seabed_Friction_Spring(d, stiffness, capacity, force, dforce_dd, dforce_dcap)
    !! Static seabed friction spring: the resisting force of a horizontal spring of
    !! stiffness k stretched by d from its reference, capped smoothly at the capacity C,
    !!   f = k d / s,   s = sqrt(1 + (k |d| / C)^2),
    !! which is k d for small d and tends to C along d once the node slides. Returns f (the
    !! seabed acts on the node with -f), df/dd = (k/s) (I - (k/(C s))^2 d d^T) and
    !! df/dC = k d (k |d|)^2 / (C^3 s^3). C <= 0 (no contact) gives no force.
    REAL(wp), INTENT(IN) :: d(2), stiffness, capacity
    REAL(wp), INTENT(OUT) :: force(2), dforce_dd(2, 2), dforce_dcap(2)
    REAL(wp) :: s, kr2, g
    force = CD_ZERO
    dforce_dd = CD_ZERO
    dforce_dcap = CD_ZERO
    IF (.NOT. (capacity > CD_ZERO) .OR. .NOT. (stiffness > CD_ZERO)) RETURN
    kr2 = stiffness*stiffness*(d(1)*d(1) + d(2)*d(2))
    s = SQRT(CD_ONE + kr2/(capacity*capacity))
    force = stiffness*d/s
    g = stiffness*stiffness/(capacity*capacity*s*s)
    dforce_dd(1, 1) = (stiffness/s)*(CD_ONE - g*d(1)*d(1))
    dforce_dd(2, 2) = (stiffness/s)*(CD_ONE - g*d(2)*d(2))
    dforce_dd(1, 2) = -(stiffness/s)*g*d(1)*d(2)
    dforce_dd(2, 1) = dforce_dd(1, 2)
    dforce_dcap = stiffness*d*kr2/(capacity*capacity*capacity*s*s*s)
  END SUBROUTINE CD_Seabed_Friction_Spring

  PURE SUBROUTINE CD_Seabed_Friction_Stick_Slip(d, stiffness, capacity, force, dforce_dd, dforce_dcap)
    !! Stick-slip seabed friction of the dynamics: a horizontal spring of stiffness k from
    !! the node's anchor, stretched by d. It sticks, f = k d, while k |d| <= C, and slides
    !! at the capacity, f = C d/|d|, beyond (the anchor then follows the node, see the
    !! callers' return mapping). Returns f (the seabed acts on the node with -f), df/dd
    !! (k I sticking, (C/|d|) (I - d d^T/|d|^2) sliding) and df/dC (0 sticking, d/|d|
    !! sliding). C <= 0 (no contact) gives no force.
    REAL(wp), INTENT(IN) :: d(2), stiffness, capacity
    REAL(wp), INTENT(OUT) :: force(2), dforce_dd(2, 2), dforce_dcap(2)
    REAL(wp) :: dlen, dhat(2)
    force = CD_ZERO
    dforce_dd = CD_ZERO
    dforce_dcap = CD_ZERO
    IF (.NOT. (capacity > CD_ZERO) .OR. .NOT. (stiffness > CD_ZERO)) RETURN
    dlen = SQRT(d(1)*d(1) + d(2)*d(2))
    IF (stiffness*dlen <= capacity) THEN
      force = stiffness*d
      dforce_dd(1, 1) = stiffness
      dforce_dd(2, 2) = stiffness
      RETURN
    END IF
    dhat = d/dlen
    force = capacity*dhat
    dforce_dd(1, 1) = (capacity/dlen)*(CD_ONE - dhat(1)*dhat(1))
    dforce_dd(2, 2) = (capacity/dlen)*(CD_ONE - dhat(2)*dhat(2))
    dforce_dd(1, 2) = -(capacity/dlen)*dhat(1)*dhat(2)
    dforce_dd(2, 1) = dforce_dd(1, 2)
    dforce_dcap = dhat
  END SUBROUTINE CD_Seabed_Friction_Stick_Slip

  PURE SUBROUTINE CD_Seabed_Friction_Mu_Dir(d, mu_axial, mu_lateral, axis, mu, q, gradq)
    !! Direction-dependent Coulomb coefficient of anisotropic seabed friction (the OrcaFlex
    !! axial/normal pair): for a horizontal slip direction d with components d_a along the
    !! unit horizontal line axis a and d_n across it,
    !!   mu(d) = Q(d)/|d|,   Q(d) = sqrt((mu_a d_a)^2 + (mu_n d_n)^2),
    !! so a slip along the line resists with mu_a and one across it with mu_n, and the
    !! force stays collinear with d. axis is the horizontal projection of the local line
    !! tangent (any length); a vertical tangent (|axis| ~ 0) has no axial direction and
    !! takes mu_n in every direction. Returns mu, Q and dQ/dd; at d = 0, mu = mu_n and
    !! dQ/dd = 0.
    REAL(wp), INTENT(IN) :: d(2), mu_axial, mu_lateral, axis(2)
    REAL(wp), INTENT(OUT) :: mu, q, gradq(2)
    REAL(wp) :: alen, a(2), nrm(2), da, dn, dlen

    alen = SQRT(axis(1)*axis(1) + axis(2)*axis(2))
    dlen = SQRT(d(1)*d(1) + d(2)*d(2))
    IF (.NOT. (alen > 1.0e-12_wp) .OR. .NOT. (dlen > CD_ZERO)) THEN
      mu = mu_lateral
      q = mu_lateral*dlen
      gradq = CD_ZERO
      IF (dlen > CD_ZERO) gradq = mu_lateral*d/dlen
      RETURN
    END IF
    a = axis/alen
    nrm = [-a(2), a(1)]
    da = d(1)*a(1) + d(2)*a(2)
    dn = d(1)*nrm(1) + d(2)*nrm(2)
    q = SQRT((mu_axial*da)**2 + (mu_lateral*dn)**2)
    mu = q/dlen
    IF (q > CD_ZERO) THEN
      gradq = (mu_axial*mu_axial*da*a + mu_lateral*mu_lateral*dn*nrm)/q
    ELSE
      gradq = CD_ZERO
    END IF
  END SUBROUTINE CD_Seabed_Friction_Mu_Dir

  PURE SUBROUTINE CD_Seabed_Friction_Aniso(law, d, stiffness, normal, mu_axial, mu_lateral, axis, &
                                           force, dforce_dd, dforce_dnormal)
    !! Anisotropic seabed friction: the static spring (law = CD_FRICTION_SPRING) or the
    !! stick-slip spring (CD_FRICTION_STICK_SLIP) of the isotropic laws above with the
    !! capacity C = mu(d) N of the direction-dependent coefficient CD_Seabed_Friction_Mu_Dir
    !! (N = normal reaction). Returns f (the seabed acts on the node with -f), df/dd
    !! (including the dependence of the capacity on the slip direction; the line axis is
    !! held fixed) and df/dN. N <= 0 or k <= 0 gives no force.
    !!   spring:     f = k d / s,  s = sqrt(1 + r^2),  r = k |d|^2 / (N Q(d))
    !!   stick-slip: f = k d while k |d| <= mu(d) N, else f = N Q(d) d / |d|^2
    INTEGER, INTENT(IN) :: law
    REAL(wp), INTENT(IN) :: d(2), stiffness, normal, mu_axial, mu_lateral, axis(2)
    REAL(wp), INTENT(OUT) :: force(2), dforce_dd(2, 2), dforce_dnormal(2)
    REAL(wp) :: mu, q, gq(2), d2, s, r, drdd(2), c
    INTEGER :: i, j

    force = CD_ZERO
    dforce_dd = CD_ZERO
    dforce_dnormal = CD_ZERO
    IF (.NOT. (normal > CD_ZERO) .OR. .NOT. (stiffness > CD_ZERO)) RETURN
    d2 = d(1)*d(1) + d(2)*d(2)
    IF (.NOT. (d2 > CD_ZERO)) THEN
      dforce_dd(1, 1) = stiffness
      dforce_dd(2, 2) = stiffness
      RETURN
    END IF
    CALL CD_Seabed_Friction_Mu_Dir(d, mu_axial, mu_lateral, axis, mu, q, gq)
    IF (.NOT. (q > CD_ZERO)) RETURN
    IF (law == CD_FRICTION_SPRING) THEN
      r = stiffness*d2/(normal*q)
      s = SQRT(CD_ONE + r*r)
      force = stiffness*d/s
      drdd = r*(2.0_wp*d/d2 - gq/q)
      c = stiffness*r/(s*s*s)
      DO j = 1, 2
        DO i = 1, 2
          dforce_dd(i, j) = -c*d(i)*drdd(j)
        END DO
        dforce_dd(j, j) = dforce_dd(j, j) + stiffness/s
      END DO
      dforce_dnormal = stiffness*d*r*r/(normal*s*s*s)
    ELSE
      IF (stiffness*SQRT(d2) <= mu*normal) THEN
        force = stiffness*d
        dforce_dd(1, 1) = stiffness
        dforce_dd(2, 2) = stiffness
        RETURN
      END IF
      force = normal*q*d/d2
      DO j = 1, 2
        DO i = 1, 2
          dforce_dd(i, j) = normal*(d(i)*gq(j)/d2 - 2.0_wp*q*d(i)*d(j)/(d2*d2))
        END DO
        dforce_dd(j, j) = dforce_dd(j, j) + normal*q/d2
      END DO
      dforce_dnormal = q*d/d2
    END IF
  END SUBROUTINE CD_Seabed_Friction_Aniso

END MODULE CableDyn_SeabedContact
