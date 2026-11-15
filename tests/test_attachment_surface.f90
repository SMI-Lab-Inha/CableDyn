! File: tests/test_attachment_surface.f90
! SPDX-License-Identifier: Apache-2.0
! Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
PROGRAM test_attachment_surface
  !! The immersed fraction of a point attachment at the free surface (the sphere of its volume,
  !! CD_HermiteCable_Point_Immersion) and the buoyancy it loses there in the finite-EI static solve.
  !! Gates:
  !!   K1 the fraction: 1 from r below the surface, 0 from r above, 1/2 at it; its height
  !!      derivative matches central FD in the band; the potential psi has the derivative wf; wf
  !!      and dwf are continuous at the band edges; r = 0 is the wet/dry switch; the sphere radius
  !!      of a volume and of a drag area.
  !!   K2 the static solve (CD_HermiteCable_Static_Solve, att_node/att_buoyancy/att_radius): with
  !!      the attachment node in the band the residual changes by rho g V (1 - wf) at its height
  !!      and the tangent's column matches central FD; deeper than r the residual and tangent are
  !!      those without the attachment term, bit for bit; above the band the residual carries
  !!      rho g V in full.
  USE CableDyn_Precision, ONLY: wp
  USE CableDyn_HermiteCableDynamic, ONLY: CD_HermiteCable_Point_Immersion, CD_HermiteCable_Point_Radius
  USE CableDyn_HermiteCableStatic, ONLY: CD_HermiteCable_Static_Solve, CD_HCSTAT_OK
  IMPLICIT NONE
  REAL(wp), PARAMETER :: PI = 3.14159265358979323846_wp
  INTEGER, PARAMETER :: NE = 10, NN = NE + 1, NDOF = 6*NN, LDAB = 2*11 + 11 + 1
  INTEGER :: nfail = 0

  CALL gate_k1()
  CALL gate_k2()
  IF (nfail > 0) THEN
    WRITE (*, '(A,I0,A)') 'FAILED: ', nfail, ' check(s)'
    ERROR STOP 1
  END IF
  WRITE (*, '(A)') 'PASS: attachment immersion at the free surface'

CONTAINS

  SUBROUTINE gate_k1()
    REAL(wp), PARAMETER :: R = 0.4_wp, H = 1.0e-6_wp
    REAL(wp) :: wf, dwf, psi, wp_, wm, dum, pp, pm, err_d, err_p, d
    INTEGER :: k
    CALL CD_HermiteCable_Point_Immersion(1.0_wp, R, wf, dwf, psi)
    CALL check(wf == 1.0_wp .AND. dwf == 0.0_wp .AND. psi == -1.0_wp, 'K1: fully submerged below r')
    CALL CD_HermiteCable_Point_Immersion(-1.0_wp, R, wf, dwf, psi)
    CALL check(wf == 0.0_wp .AND. dwf == 0.0_wp .AND. psi == 0.0_wp, 'K1: dry above r')
    CALL CD_HermiteCable_Point_Immersion(0.0_wp, R, wf, dwf, psi)
    CALL check(ABS(wf - 0.5_wp) < 1.0e-15_wp .AND. ABS(dwf + 3.0_wp/(4.0_wp*R)) < 1.0e-14_wp, &
               'K1: half immersed at the surface, slope -3/(4 r)')
    err_d = 0.0_wp
    err_p = 0.0_wp
    DO k = -19, 19
      d = R*REAL(k, wp)/20.0_wp
      CALL CD_HermiteCable_Point_Immersion(d, R, wf, dwf, psi)
      ! z = -depth: d wf/d z = -(d wf/d depth)
      CALL CD_HermiteCable_Point_Immersion(d - H, R, wp_, dum, pp)
      CALL CD_HermiteCable_Point_Immersion(d + H, R, wm, dum, pm)
      err_d = MAX(err_d, ABS((wp_ - wm)/(2.0_wp*H) - dwf))
      err_p = MAX(err_p, ABS((pp - pm)/(2.0_wp*H) - wf))
    END DO
    WRITE (*, '(A,2ES10.2)') '  K1 FD error of dwf/dz and dpsi/dz in the band: ', err_d, err_p
    CALL check(err_d < 1.0e-8_wp, 'K1: dwf/dz matches central FD')
    CALL check(err_p < 1.0e-8_wp, 'K1: dpsi/dz = wf')
    CALL CD_HermiteCable_Point_Immersion(R*(1.0_wp - 1.0e-9_wp), R, wf, dwf, psi)
    CALL check(ABS(wf - 1.0_wp) < 1.0e-12_wp .AND. ABS(dwf) < 1.0e-8_wp .AND. ABS(psi + R) < 1.0e-9_wp, &
               'K1: continuous with its slope and potential at the lower band edge')
    CALL CD_HermiteCable_Point_Immersion(-R*(1.0_wp - 1.0e-9_wp), R, wf, dwf, psi)
    CALL check(ABS(wf) < 1.0e-12_wp .AND. ABS(dwf) < 1.0e-8_wp .AND. ABS(psi) < 1.0e-12_wp, &
               'K1: continuous with its slope and potential at the upper band edge')
    CALL CD_HermiteCable_Point_Immersion(0.0_wp, 0.0_wp, wf, dwf, psi)
    CALL check(wf == 1.0_wp .AND. dwf == 0.0_wp, 'K1: a point (r = 0) is wet at the surface')
    CALL CD_HermiteCable_Point_Immersion(-1.0e-12_wp, 0.0_wp, wf, dwf, psi)
    CALL check(wf == 0.0_wp, 'K1: a point (r = 0) is dry above the surface')
    CALL check(ABS(CD_HermiteCable_Point_Radius(4.0_wp*PI/3.0_wp, 0.5_wp) - 1.0_wp) < 1.0e-14_wp, &
               'K1: the sphere of the volume')
    CALL check(ABS(CD_HermiteCable_Point_Radius(0.0_wp, PI) - 1.0_wp) < 1.0e-14_wp, &
               'K1: without a volume, the sphere of the drag area')
    CALL check(CD_HermiteCable_Point_Radius(0.0_wp, 0.0_wp) == 0.0_wp, 'K1: neither: a point')
  END SUBROUTINE gate_k1

  SUBROUTINE gate_k2()
    !! A 10 m line of ten elements along x at z0, attachment at node 6 (rho g V 2000 N, r 0.3 m).
    REAL(wp), PARAMETER :: B = 2000.0_wp, R = 0.3_wp, HH = 1.0e-6_wp
    REAL(wp) :: q(NDOF), r0(NDOF), r1(NDOF), ab0(LDAB, NDOF), ab1(LDAB, NDOF), rp(NDOF), rm(NDOF), dum(LDAB, NDOF)
    REAL(wp) :: wf, dwf, col_fd(NDOF), col_an(NDOF), err, z0
    INTEGER :: jz, i
    jz = 6*5 + 3
    ! in the band
    z0 = 0.1_wp
    CALL line_state(z0, q)
    CALL evaluate(q, .FALSE., r0, ab0)
    CALL evaluate(q, .TRUE., r1, ab1)
    CALL CD_HermiteCable_Point_Immersion(-z0, R, wf, dwf)
    CALL check(ABS(r1(jz) - r0(jz) - B*(1.0_wp - wf)) < 1.0e-9_wp*B, &
               'K2: the residual carries the lost buoyancy rho g V (1 - wf) at the node height')
    CALL check(MAXVAL(ABS(r1(1:jz - 1) - r0(1:jz - 1))) == 0.0_wp .AND. &
               MAXVAL(ABS(r1(jz + 1:NDOF) - r0(jz + 1:NDOF))) == 0.0_wp, 'K2: only the node''s vertical DOF changes')
    q(jz) = z0 + HH
    CALL evaluate(q, .TRUE., rp, dum)
    q(jz) = z0 - HH
    CALL evaluate(q, .TRUE., rm, dum)
    col_fd = (rp - rm)/(2.0_wp*HH)
    col_an = 0.0_wp
    DO i = MAX(1, jz - 11), MIN(NDOF, jz + 11)
      col_an(i) = ab1(2*11 + 1 + i - jz, jz)
    END DO
    err = MAXVAL(ABS(col_fd - col_an))/MAX(1.0_wp, MAXVAL(ABS(col_an)))
    WRITE (*, '(A,ES10.2)') '  K2 static tangent column vs central FD (node in the band): ', err
    CALL check(err < 1.0e-6_wp, 'K2: the static tangent column of the node matches central FD')
    CALL check(ABS((ab1(2*11 + 1, jz) - ab0(2*11 + 1, jz)) + B*dwf) < 1.0e-9_wp*B, &
               'K2: the tangent adds the buoyancy stiffness -rho g V dwf/dz')
    ! deeper than r: unchanged, bit for bit
    CALL line_state(-0.31_wp, q)
    CALL evaluate(q, .FALSE., r0, ab0)
    CALL evaluate(q, .TRUE., r1, ab1)
    CALL check(ALL(r1 == r0) .AND. ALL(ab1 == ab0), 'K2: an attachment deeper than r leaves the solve bit-identical')
    ! above the band: the whole buoyancy is lost
    CALL line_state(0.31_wp, q)
    CALL evaluate(q, .FALSE., r0, ab0)
    CALL evaluate(q, .TRUE., r1, ab1)
    CALL check(ABS(r1(jz) - r0(jz) - B) < 1.0e-12_wp*B .AND. ALL(ab1 == ab0), &
               'K2: above the band the attachment loses its whole buoyancy, with no stiffness')
  END SUBROUTINE gate_k2

  SUBROUTINE line_state(z0, q)
    !! A straight line along x at height z0, stretched 1e-4 (a non-trivial elastic residual).
    REAL(wp), INTENT(IN) :: z0
    REAL(wp), INTENT(OUT) :: q(NDOF)
    INTEGER :: i
    q = 0.0_wp
    DO i = 1, NN
      q(6*(i - 1) + 1) = 1.0001_wp*REAL(i - 1, wp)
      q(6*(i - 1) + 3) = z0 + 0.002_wp*SIN(REAL(i, wp))
      q(6*(i - 1) + 4) = 1.0001_wp
    END DO
    q(6*5 + 3) = z0
  END SUBROUTINE line_state

  SUBROUTINE evaluate(qx, with_att, r, ab)
    !! The solver's residual and banded tangent at qx (no Newton update: tol = HUGE).
    REAL(wp), INTENT(IN) :: qx(NDOF)
    LOGICAL, INTENT(IN) :: with_att
    REAL(wp), INTENT(OUT) :: r(NDOF), ab(LDAB, NDOF)
    REAL(wp) :: l0(NE), ea(NE), ei(NE), w(NE), qd(NDOF), cd(NN), rr
    INTEGER :: itx, esx
    CHARACTER(512) :: emx
    l0 = 1.0_wp; ea = 1.0e6_wp; ei = 50.0_wp; w = 40.0_wp
    IF (with_att) THEN
      CALL CD_HermiteCable_Static_Solve(l0, ea, ei, w, qx, [1, 2, 3, NDOF - 5, NDOF - 4, NDOF - 3], -1.0e4_wp, &
                                        0.0_wp, 1, 1, HUGE(1.0_wp), 1.0_wp, qd, cd, rr, itx, esx, emx, &
                                        residual_out=r, tangent_out=ab, att_node=[6], att_buoyancy=[2000.0_wp], &
                                        att_radius=[0.3_wp])
    ELSE
      CALL CD_HermiteCable_Static_Solve(l0, ea, ei, w, qx, [1, 2, 3, NDOF - 5, NDOF - 4, NDOF - 3], -1.0e4_wp, &
                                        0.0_wp, 1, 1, HUGE(1.0_wp), 1.0_wp, qd, cd, rr, itx, esx, emx, &
                                        residual_out=r, tangent_out=ab)
    END IF
    CALL check(esx == CD_HCSTAT_OK, 'residual and tangent evaluation: '//TRIM(emx))
  END SUBROUTINE evaluate

  SUBROUTINE check(ok, msg)
    LOGICAL, INTENT(IN) :: ok
    CHARACTER(*), INTENT(IN) :: msg
    IF (ok) THEN
      WRITE (*, '(2A)') 'ok    ', msg
    ELSE
      WRITE (*, '(2A)') 'FAIL  ', msg
      nfail = nfail + 1
    END IF
  END SUBROUTINE check
END PROGRAM test_attachment_surface
