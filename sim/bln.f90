     MODULE BLN 

     IMPLICIT NONE

     CONTAINS
!
!> @brief Calculate the energy and gradient for a given configuration of a BLN polymer chain.
!>
!> @param[out] ENERGY
!> @param[in] GRADT
!

      SUBROUTINE BLN(QO,GRAD,ENERGY,GRADT)
C{{{
      USE COMMONS

      IMPLICIT NONE
      DOUBLE PRECISION QO(3*NATOMS), GRAD(3*NATOMS)
      LOGICAL GRADT
      DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XR(NATOMS,NATOMS), YR(NATOMS,NATOMS), 
     &                 ZR(NATOMS,NATOMS), DOT_PROD(NATOMS,3),
     &                 X_PROD(NATOMS), BOND_ANGLE(NATOMS), TOR_ANGLE(NATOMS), RADII(NATOMS,NATOMS), 
     &                 ENERGY, COSTOR(NATOMS), DFAC(NATOMS), SINBOND(NATOMS)

!
! Without these initialisations the NAG compiler fills in random numbers for
! unassigned elements with optimisation turned on.
!
      BOND_ANGLE(1:NATOMS)=0.0D0
      TOR_ANGLE(1:NATOMS)=0.0D0
      COSTOR(1:NATOMS)=0.0D0
      DFAC(1:NATOMS)=0.0D0
      SINBOND(1:NATOMS)=0.0D0
      DOT_PROD(1:NATOMS,1:3)=0.0D0
      X_PROD(1:NATOMS)=0.0D0
      RADII(1:NATOMS,1:NATOMS)=0.0D0
      XR(1:NATOMS,1:NATOMS)=0.0D0
      YR(1:NATOMS,1:NATOMS)=0.0D0
      ZR(1:NATOMS,1:NATOMS)=0.0D0

      CALL CALC_INT_COORDS_BLN(QO,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,&
                BOND_ANGLE,TOR_ANGLE,RADII,&
                NATOMS,COSTOR,SINBOND,A_BLN,B_BLN,C_BLN,D_BLN,DFAC)

      CALL CALC_ENERGY_BLN(QO,ENERGY,LJREP_BLN,LJATT_BLN,&
                A_BLN,B_BLN,C_BLN,D_BLN,&
                X,Y,Z,& 
                XR,YR,ZR,DOT_PROD,X_PROD,BOND_ANGLE,TOR_ANGLE,RADII,&
                NATOMS,RK_R,RK_THETA,COSTOR)

      IF (.NOT.GRADT) RETURN
 
      CALL CALC_GRADIENT_BLN(QO,GRAD,LJREP_BLN,LJATT_BLN,&
         A_BLN,B_BLN,C_BLN,D_BLN,&
         X,Y,Z,&
         XR,YR,ZR,&
         DOT_PROD,X_PROD,&
         BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,&
         RK_R,RK_THETA,COSTOR,DFAC,SINBOND)

      RETURN
      END
C }}}

!> @brief Calculate the internal coordinates of a BLN chain

      SUBROUTINE CALC_INT_COORDS_BLN(QO,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,BOND_ANGLE,TOR_ANGLE,RADII,NATOMS, 
     &                              COSTOR,SINBOND,A_BLN,B_BLN,C_BLN,D_BLN,DFAC)
! {{{
! Declarations {{{
      IMPLICIT NONE

      INTEGER I, J, N

      DOUBLE PRECISION COS_PHI, COS_THETA, DUMMY, DUMMY2
      DOUBLE PRECISION QO(3*NATOMS)
      DOUBLE PRECISION, PARAMETER :: TWOPI=6.283185307179586477D0
      DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XR(NATOMS,NATOMS), DFAC(NATOMS),
     &  YR(NATOMS,NATOMS), ZR(NATOMS,NATOMS), DOT_PROD(NATOMS,3), COSTOR(NATOMS), D_BLN(NATOMS), 
     &  A_BLN(NATOMS), B_BLN(NATOMS), X_PROD(NATOMS), BOND_ANGLE(NATOMS), TOR_ANGLE(NATOMS), 
     &  RADII(NATOMS,NATOMS), SINBOND(NATOMS), C_BLN(NATOMS)
! }}}

      N=NATOMS

      DO I = 1, N
         J = (I-1)*3
         X(I) = QO((I-1)*3+1)
         Y(I) = QO((I-1)*3+2)
         Z(I) = QO((I-1)*3+3)
      ENDDO
!
! Inter-particle distances {{{

      DO I = 1, N-1
         DO J = I+1, N
            XR(I,J) = X(J) - X(I)
            YR(I,J) = Y(J) - Y(I)
            ZR(I,J) = Z(J) - Z(I)
            RADII(I,J) = SQRT(XR(I,J)*XR(I,J) + YR(I,J)*YR(I,J) + ZR(I,J)*ZR(I,J))
            RADII(J,I) = RADII(I,J)
         ENDDO
      ENDDO
C }}}

! Dot products between bond vectors {{{

      DO I = 1, N-3
         DOT_PROD(I,1) = XR(I,I+1)*XR(I,I+1) + YR(I,I+1)*YR(I,I+1) + ZR(I,I+1)*ZR(I,I+1)
         DOT_PROD(I,2) = XR(I,I+1)*XR(I+1,I+2)+YR(I,I+1)*YR(I+1,I+2)+ ZR(I,I+1)*ZR(I+1,I+2)
         DOT_PROD(I,3) = XR(I,I+1)*XR(I+2,I+3)+YR(I,I+1)*YR(I+2,I+3)+ ZR(I,I+1)*ZR(I+2,I+3)
      ENDDO

!     I = N-2
      DOT_PROD(N-2,1) = XR(N-2,N-2+1)*XR(N-2,N-2+1) + YR(N-2,N-2+1)*YR(N-2,N-2+1) + ZR(N-2,N-2+1)*ZR(N-2,N-2+1)
      DOT_PROD(N-2,2) = XR(N-2,N-2+1)*XR(N-2+1,N-2+2)+YR(N-2,N-2+1)*YR(N-2+1,N-2+2)+ ZR(N-2,N-2+1)*ZR(N-2+1,N-2+2)
!     I = N-1
      DOT_PROD(N-1,1) = XR(N-1,N-1+1)*XR(N-1,N-1+1) + YR(N-1,N-1+1)*YR(N-1,N-1+1) + ZR(N-1,N-1+1)*ZR(N-1,N-1+1)
! }}}

! Cross-products between adjacent bond vectors {{{

      do i = 1, n-2
         x_prod(i) = dot_prod(i,1)*dot_prod(i+1,1) - dot_prod(i,2)*dot_prod(i,2)   
      enddo
! }}}
! Bond angles {{{

      do i = 1, n-2
         cos_theta=-dot_prod(i,2)/(sqrt(dot_prod(i,1)*dot_prod(i+1,1)))
         bond_angle(i+1) = dacos(cos_theta)
         SINBOND(i+1)=SIN(bond_angle(i+1))*SQRT(dot_prod(i,1)*dot_prod(i+1,1))
      enddo
! }}}
! Torsional angles {{{

      do i = 1, n-3
         cos_phi = (dot_prod(i,2)*dot_prod(i+1,2) - dot_prod(i,3)*dot_prod(i+1,1))/sqrt(x_prod(i)*x_prod(i+1))
         IF (ABS(cos_phi).GT.1.0D0) cos_phi=SIGN(1.0D0,cos_phi)
         tor_angle(i+1) = dacos(cos_phi)
! }}}
!
! tor_angle is returned in the range 0 to Pi.
! dummy should take the opposite sign from the dihedral angle. Negative
! values of the dihedral should be subtracted from 2*pi.
! This is only necessary when the potential contains cos(phi+pi/4) terms because
! the gradient is discontinuous when phi goes through pi for such terms if phi
! is restricted to 0 < phi < pi.
! {{{
!        dummy=xr(i+2,i+3)*(yr(i+1,i)*zr(i+1,i+2)-yr(i+1,i+2)*zr(i+1,i))+
!    &         yr(i+2,i+3)*(xr(i+1,i+2)*zr(i+1,i)-xr(i+1,i)*zr(i+1,i+2))+
!    &         zr(i+2,i+3)*(xr(i+1,i)*yr(i+1,i+2)-xr(i+1,i+2)*yr(i+1,i))
         dummy=xr(i+2,i+3)*(-yr(i,i+1)*zr(i+1,i+2)+yr(i+1,i+2)*zr(i,i+1))+
     &         yr(i+2,i+3)*(-xr(i+1,i+2)*zr(i,i+1)+xr(i,i+1)*zr(i+1,i+2))+
     &         zr(i+2,i+3)*(-xr(i,i+1)*yr(i+1,i+2)+xr(i+1,i+2)*yr(i,i+1))
         IF (DUMMY.GT.0.0D0) tor_angle(i+1)=TWOPI-tor_angle(i+1)
         COSTOR(i+1)=COS(tor_angle(i+1))
! }}}
!  This is an ugly hack to prevent division by zero. There will be a loss of precision
!  if a dihedral is 0, PI, 2*PI, if D_BLN is non-zero.
! {{{
         IF (TAN(tor_angle(i+1)).EQ.0.0D0) THEN
            PRINT '(A,I8,A,G20.10)','WARNING in BLN, dihedral angle ',i+1,' is ',tor_angle(i+1)
            tor_angle(i+1)=tor_angle(i+1)+1.0D-10
            PRINT '(A,G20.10)','WARNING in BLN, TAN perturbed angle=',tor_angle(i+1)
         ENDIF
         DUMMY2=TAN(tor_angle(i+1))
         DFAC(i+1)=(A_BLN(i+1)+D_BLN(i+1)*( 1.0D0+1.0D0/DUMMY2 )*0.7071067811865475244D0-B_BLN(i+1)
     &              +C_BLN(i+1)*(12.0*costor(i+1)**2-3.0))/sqrt(x_prod(i+1)*x_prod(i))
! }}}
      ENDDO

      RETURN
      END
C }}}

!> @brief Calculate the energy of a BLN chain

      SUBROUTINE CALC_ENERGY_BLN(QO,ENERGY,N,LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     &  BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,RK_R,RK_THETA,COSTOR)
C {{{
      IMPLICIT NONE
      INTEGER NATOMS, N, I, J
      DOUBLE PRECISION, parameter :: theta_0 = 1.8326D0, PI4=0.7853981633974483096D0 ! 1.8326 radians is 105 degrees
      DOUBLE PRECISION QO(3*NATOMS), ENERGY
      DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XR(NATOMS,NATOMS), 
     &  YR(NATOMS,NATOMS), ZR(NATOMS,NATOMS), DOT_PROD(NATOMS,3), COSTOR(NATOMS),
     &  X_PROD(NATOMS), BOND_ANGLE(NATOMS), TOR_ANGLE(NATOMS), RADII(NATOMS,NATOMS)
      DOUBLE PRECISION RK_R, RK_THETA, E_NBOND, E_BOND, E_BANGLE, E_TANGLE, RAD6
      DOUBLE PRECISION LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),C_BLN(NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),D_BLN(NATOMS)

      e_nbond=0.0D0
      e_bond=0.0D0
      e_bangle=0.0D0
      e_tangle=0.0D0

      do i = 1, n-2
         do j = i+2, n
            rad6 = radii(i,j)**6
            e_nbond = e_nbond + (LJREP_BLN(i,j)/rad6 + LJATT_BLN(i,j))/rad6
         enddo
      enddo
      e_nbond=e_nbond*4.0D0

      do i = 1, n-1
         e_bond = e_bond + (radii(i,i+1)-1.0D0)**2
      enddo
      e_bond=e_bond*rk_r/2.0D0

      do i = 2, n-1
         e_bangle = e_bangle + (bond_angle(i)-theta_0)**2
      enddo
      e_bangle=e_bangle*rk_theta/2.0D0

      do i = 2, n-2
         e_tangle = e_tangle + A_BLN(i)*(1.0D0 + costor(i)) + C_BLN(i)*(1.0 + cos(3.0*tor_angle(i)))
     &                       + B_BLN(i)*(1.0D0 - costor(i)) + D_BLN(i)*(1.0 + cos(tor_angle(i)+PI4))
      enddo

      energy = e_nbond + e_bond + e_bangle + e_tangle

C     write(*,'(A,4F20.10)') 'nbond,bond,bangle,tangle=',e_nbond,e_bond,e_bangle,e_tangle

      return
      end
C }}}

!> @brief Calculate the gradients

      SUBROUTINE CALC_GRADIENT_BLN(QO,FQ,N,LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,
     &                            X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     &                            BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,RK_R,RK_THETA,COSTOR,DFAC,SINBOND)
C {{{
C Declarations {{{
      USE COMMONS,ONLY : MYUNIT
      IMPLICIT NONE
      INTEGER NATOMS, N, I, J
      DOUBLE PRECISION, parameter :: theta_0 = 1.8326D0
      DOUBLE PRECISION qo(3*NATOMS),fq(3*NATOMS),fx(NATOMS),fy(NATOMS), fz(NATOMS), DFAC(NATOMS), SINBOND(NATOMS)
      DOUBLE PRECISION fnb_x(NATOMS),fnb_y(NATOMS),fnb_z(NATOMS), fb_x(NATOMS),fb_y(NATOMS)
      DOUBLE PRECISION fb_z(NATOMS),fba_x(NATOMS),fba_y(NATOMS),fba_z(NATOMS), COSTOR(NATOMS)
      DOUBLE PRECISION fta_x(NATOMS),fta_y(NATOMS),fta_z(NATOMS), a3, coef, coef1, coef2, coef3, a4
      DOUBLE PRECISION RK_R, RK_THETA, rad7, rad14, df, fxx, fzz, fyy, rvar, den, rnum, den1, a1, a2, den2

      DOUBLE PRECISION x(NATOMS), y(NATOMS), z(NATOMS), xr(NATOMS,NATOMS), 
     &  yr(NATOMS,NATOMS), zr(NATOMS,NATOMS), dot_prod(NATOMS,3), 
     &  x_prod(NATOMS), bond_angle(NATOMS), tor_angle(NATOMS), radii(NATOMS,NATOMS)

      DOUBLE PRECISION LJREP_BLN(NATOMS,NATOMS), LJATT_BLN(NATOMS,NATOMS),C_BLN(NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),D_BLN(NATOMS)
C }}}
C
C Gradients of potential
C {{{
      do i = 1,n

         fnb_x(i) = 0.0  
         fnb_y(i) = 0.0 
         fnb_z(i) = 0.0 
   
         fb_x(i)  = 0.0 
         fb_y(i)  = 0.0 
         fb_z(i)  = 0.0 
   
         fba_x(i) = 0.0 
         fba_y(i) = 0.0 
         fba_z(i) = 0.0 
   
         fta_x(i) = 0.0 
         fta_y(i) = 0.0 
         fta_z(i) = 0.0 
   
         fx(i)= 0.0 
         fy(i)= 0.0 
         fz(i)= 0.0 

      enddo
C }}}
C ..... Non-bonded interaction forces ..... 
C {{{
      do i = 1, n-2
         do j = i+2, n

            rad7 = radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)   
            rad14 = rad7*rad7 

            df = -24.0*((2.0*LJREP_BLN(i,j)/rad14) + (LJATT_BLN(i,j)/(rad7*radii(i,j))))

            fxx = df*xr(i,j) 
            fyy = df*yr(i,j) 
            fzz = df*zr(i,j) 

            fnb_x(i) = fxx + fnb_x(i)
            fnb_y(i) = fyy + fnb_y(i)
            fnb_z(i) = fzz + fnb_z(i)

            fnb_x(j) = -fxx + fnb_x(j)
            fnb_y(j) = -fyy + fnb_y(j)
            fnb_z(j) = -fzz + fnb_z(j)

         enddo
      enddo
C }}}
C ... Bond interaction forces ... 
C {{{
      do i = 1, n-1

         rvar = 1.0D0/radii(i,i+1) 

         df = rk_r*(1.0 - rvar) 
         fxx = df*xr(i,i+1) 
         fyy = df*yr(i,i+1) 
         fzz = df*zr(i,i+1) 

         fb_x(i) = fxx + fb_x(i)
         fb_y(i) = fyy + fb_y(i)
         fb_z(i) = fzz + fb_z(i)

         fb_x(i+1) = -fxx + fb_x(i+1)
         fb_y(i+1) = -fyy + fb_y(i+1)
         fb_z(i+1) = -fzz + fb_z(i+1)

      enddo
C }}}
C bond angle forces  particle 1
C particles 1,2,n-1, and n done outside of the loop
C {{{
!     i = 1
      den = sinbond(1+1)
      rnum = rk_theta*(bond_angle(1+1) - theta_0)

      fba_x(1) = -rnum*((dot_prod(1,2)/dot_prod(1,1))*xr(1,1+1) - xr(1+1,1+2))/den

      fba_y(1) = -rnum*((dot_prod(1,2)/dot_prod(1,1))*yr(1,1+1) - yr(1+1,1+2))/den

      fba_z(1) = -rnum*((dot_prod(1,2)/dot_prod(1,1))*zr(1,1+1) - zr(1+1,1+2))/den

C }}}
C particle 2
C {{{
!     i = 2
      den = sinbond(2)
      den1 = sinbond(3)

      a1 = -rk_theta*(bond_angle(2) - theta_0)*( (dot_prod(2-1,2)/
     1  dot_prod(2,1))*xr(2,2+1) - (dot_prod(2-1,2)/dot_prod(2-1,1))
     1      *xr(2-1,2) + xr(2,2+1) - xr(2-1,2))/den

      a2 = -rk_theta*(bond_angle(2+1) - theta_0)*((dot_prod(2,2)/dot_prod(2,1))*xr(2,2+1) - xr(2+1,2+2))/den1

      fba_x(2) = a1 + a2 

      a1 = -rk_theta*(bond_angle(2) - theta_0)*( (dot_prod(2-1,2)/
     1  dot_prod(2,1))*yr(2,2+1) - (dot_prod(2-1,2)/dot_prod(2-1,1))
     1      *yr(2-1,2) + yr(2,2+1) - yr(2-1,2))/den

      a2 = -rk_theta*(bond_angle(2+1) - theta_0)*((dot_prod(2,2)/dot_prod(2,1))*yr(2,2+1) - yr(2+1,2+2))/den1

      fba_y(2) = a1 + a2 

      a1 = -rk_theta*(bond_angle(2) - theta_0)*( (dot_prod(2-1,2)/
     1  dot_prod(2,1))*zr(2,2+1) - (dot_prod(2-1,2)/dot_prod(2-1,1))
     1      *zr(2-1,2) + zr(2,2+1) - zr(2-1,2))/den

      a2 = -rk_theta*(bond_angle(2+1) - theta_0)*((dot_prod(2,2)/dot_prod(2,1))*zr(2,2+1) - zr(2+1,2+2))/den1

      fba_z(2) = a1 + a2 

C }}}
C particles 3 thru n-2 
C {{{
      do i = 3, n-2

         den = sinbond(i)
         den1 = sinbond(i+1)
         den2 = sinbond(i-1)

         a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1      *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

         a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

         a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den2

         fba_x(i) = a1 + a2 + a3 

         a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1   dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1      *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

         a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1      dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

         a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1      dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den2

         fba_y(i) = a1 + a2 + a3 

         a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1   dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1      *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

         a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1      dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

         a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1      dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den2

         fba_z(i) = a1 + a2 + a3 

      enddo
C }}}
C particle n-1 
C {{{
!     i = n-1
      den = sinbond(n-1)
      den1 = sinbond(n-1-1)

      a1 = -rk_theta*(bond_angle(n-1) - theta_0)*( (dot_prod(n-1-1,2)/
     1   dot_prod(n-1,1))*xr(n-1,n-1+1) - (dot_prod(n-1-1,2)/dot_prod(n-1-1,1))
     1      *xr(n-1-1,n-1) + xr(n-1,n-1+1) - xr(n-1-1,n-1))/den

      a2 = rk_theta*(bond_angle(n-1-1) - theta_0)*((dot_prod(n-1-2,2)/
     1      dot_prod(n-1-1,1))*xr(n-1-1,n-1) - xr(n-1-2,n-1-1))/den1

      fba_x(n-1) = a1 + a2

      a1 = -rk_theta*(bond_angle(n-1) - theta_0)*( (dot_prod(n-1-1,2)/
     1   dot_prod(n-1,1))*yr(n-1,n-1+1) - (dot_prod(n-1-1,2)/dot_prod(n-1-1,1))
     1      *yr(n-1-1,n-1) + yr(n-1,n-1+1) - yr(n-1-1,n-1))/den
   
      a2 = rk_theta*(bond_angle(n-1-1) - theta_0)*((dot_prod(n-1-2,2)/
     1      dot_prod(n-1-1,1))*yr(n-1-1,n-1) - yr(n-1-2,n-1-1))/den1

      fba_y(n-1) = a1 + a2

      a1 = -rk_theta*(bond_angle(n-1) - theta_0)*( (dot_prod(n-1-1,2)/
     1  dot_prod(n-1,1))*zr(n-1,n-1+1) - (dot_prod(n-1-1,2)/dot_prod(n-1-1,1))
     1      *zr(n-1-1,n-1) + zr(n-1,n-1+1) - zr(n-1-1,n-1))/den

      a2 = rk_theta*(bond_angle(n-1-1) - theta_0)*((dot_prod(n-1-2,2)/
     1      dot_prod(n-1-1,1))*zr(n-1-1,n-1) - zr(n-1-2,n-1-1))/den1

      fba_z(n-1) = a1 + a2
C }}}
C particle n
C {{{
!     i = n
      den = sinbond(n-1)

      fba_x(n) = rk_theta*(bond_angle(n-1) - theta_0)*
     1      ((dot_prod(n-2,2)/dot_prod(n-1,1))*xr(n-1,n) 
     1      - xr(n-2,n-1))/den

      fba_y(n) = rk_theta*(bond_angle(n-1) - theta_0)*
     1      ((dot_prod(n-2,2)/dot_prod(n-1,1))*yr(n-1,n) 
     1      - yr(n-2,n-1))/den

      fba_z(n) = rk_theta*(bond_angle(n-1) - theta_0)*
     1      ((dot_prod(n-2,2)/dot_prod(n-1,1))*zr(n-1,n) 
     1      - zr(n-2,n-1))/den

C }}}
C Torsional angle forces
C particles 1, 2, 3, n-2, n-1, and n are done outside of the loop
C particle 1
C {{{
!     i = 1
      coef =DFAC(1+1)

      fta_x(1) = -coef*(-dot_prod(1+1,2)*xr(1+1,1+2) +
     1       dot_prod(1+1,1)*xr(1+2,1+3) -
     1      (1.0/x_prod(1))*(dot_prod(1+1,2)*dot_prod(1,2) -
     1      dot_prod(1,3)*dot_prod(1+1,1))*(-dot_prod(1+1,1)*xr(1,1+1) +
     1      dot_prod(1,2)*xr(1+1,1+2))) 

      fta_y(1) = -coef*(-dot_prod(1+1,2)*yr(1+1,1+2) +
     1      dot_prod(1+1,1)*yr(1+2,1+3) -
     1      (1.0/x_prod(1))*(dot_prod(1+1,2)*dot_prod(1,2) -
     1      dot_prod(1,3)*dot_prod(1+1,1))*(-dot_prod(1+1,1)*yr(1,1+1) +
     1      dot_prod(1,2)*yr(1+1,1+2))) 

      fta_z(1) = -coef*(-dot_prod(1+1,2)*zr(1+1,1+2) +
     1      dot_prod(1+1,1)*zr(1+2,1+3) -
     1      (1.0/x_prod(1))*(dot_prod(1+1,2)*dot_prod(1,2) -
     1      dot_prod(1,3)*dot_prod(1+1,1))*(-dot_prod(1+1,1)*zr(1,1+1) +
     1      dot_prod(1,2)*zr(1+1,1+2))) 

C }}}
C particle 2
C {{{
!     i = 2
      coef =DFAC(2+1)

      coef1 = DFAC(2)

      a1 =  -coef*(-dot_prod(2+1,2)*xr(2+1,2+2) +
     1      dot_prod(2+1,1)*xr(2+2,2+3) -
     1      (1.0/x_prod(2))*(dot_prod(2+1,2)*dot_prod(2,2) -
     1      dot_prod(2,3)*dot_prod(2+1,1))*(-dot_prod(2+1,1)*xr(2,2+1) +
     1      dot_prod(2,2)*xr(2+1,2+2))) 

      a2 = -coef1*(-dot_prod(2-1,2)*xr(2+1,2+2) +
     1      dot_prod(2,2)*xr(2,2+1) - dot_prod(2,2)*xr(2-1,2) -
     1      dot_prod(2,1)*xr(2+1,2+2) + 2.0*dot_prod(2-1,3)*xr(2,2+1) -
     1      (1.0/x_prod(2-1))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(dot_prod(2,1)*xr(2-1,2) -
     1      dot_prod(2-1,1)*xr(2,2+1) - dot_prod(2-1,2)*xr(2,2+1) +
     1      dot_prod(2-1,2)*xr(2-1,2)) -
     1      (1.0/x_prod(2))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(-dot_prod(2+1,1)*xr(2,2+1) +
     1      dot_prod(2,2)*xr(2+1,2+2))) 

      fta_x(2) = a1 + a2 

      a1 = -dot_prod(3,2)*yr(3,4) + dot_prod(3,1)*yr(4,5) 
      a1=a1-(dot_prod(3,2)*dot_prod(2,2) - dot_prod(2,3)*dot_prod(3,1))*
     &     (-dot_prod(3,1)*yr(2,3) + dot_prod(2,2)*yr(3,4))/x_prod(2)
      a1=-coef*a1

      a2 = -coef1*(-dot_prod(2-1,2)*yr(2+1,2+2) +
     1      dot_prod(2,2)*yr(2,2+1) - dot_prod(2,2)*yr(2-1,2) -
     1      dot_prod(2,1)*yr(2+1,2+2) + 2.0*dot_prod(2-1,3)*yr(2,2+1) -
     1      (1.0/x_prod(2-1))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(dot_prod(2,1)*yr(2-1,2) -
     1      dot_prod(2-1,1)*yr(2,2+1) - dot_prod(2-1,2)*yr(2,2+1) +
     1      dot_prod(2-1,2)*yr(2-1,2)) -
     1      (1.0/x_prod(2))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(-dot_prod(2+1,1)*yr(2,2+1) +
     1      dot_prod(2,2)*yr(2+1,2+2)))

      fta_y(2) = a1 + a2 
!     WRITE(MYUNIT,'(A,I6,5G15.5)') 'y i,a1,a2,coef,coef1,fta_y=',2,a1,a2,coef,coef1,fta_y(2)
      
      a1 = -coef*(-dot_prod(3,2)*zr(3,4) +
     1      dot_prod(3,1)*zr(4,5) -
     1      (1.0/x_prod(2))*(dot_prod(3,2)*dot_prod(2,2) -
     1      dot_prod(2,3)*dot_prod(3,1))*(-dot_prod(3,1)*zr(2,3) +
     1      dot_prod(2,2)*zr(3,4))) 

      a2 = -coef1*(-dot_prod(2-1,2)*zr(2+1,2+2) +
     1      dot_prod(2,2)*zr(2,2+1) - dot_prod(2,2)*zr(2-1,2) -
     1      dot_prod(2,1)*zr(2+1,2+2) + 2.0*dot_prod(2-1,3)*zr(2,2+1) -
     1      (1.0/x_prod(2-1))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(dot_prod(2,1)*zr(2-1,2) -
     1      dot_prod(2-1,1)*zr(2,2+1) - dot_prod(2-1,2)*zr(2,2+1) +
     1      dot_prod(2-1,2)*zr(2-1,2)) -
     1      (1.0/x_prod(2))*(dot_prod(2,2)*dot_prod(2-1,2) -
     1      dot_prod(2-1,3)*dot_prod(2,1))*(-dot_prod(2+1,1)*zr(2,2+1) +
     1      dot_prod(2,2)*zr(2+1,2+2))) 

      fta_z(2) = a1 + a2 
!     WRITE(MYUNIT,'(A,I6,4G15.5)') 'z i,a1,a2,coef,coef1=',2,a1,a2,coef,coef1
C }}}
C particle 3
C {{{
!     i = 3
      coef=DFAC(3+1)

      coef1=DFAC(3)

      coef2=DFAC(3-1)

      a1 = -coef*(-dot_prod(3+1,2)*xr(3+1,3+2) +
     1      dot_prod(3+1,1)*xr(3+2,3+3) -
     1      (1.0/x_prod(3))*(dot_prod(3+1,2)*dot_prod(3,2) -
     1      dot_prod(3,3)*dot_prod(3+1,1))*(-dot_prod(3+1,1)*xr(3,3+1) +
     1      dot_prod(3,2)*xr(3+1,3+2))) 

      a2 = -coef1*(-dot_prod(3-1,2)*xr(3+1,3+2) +
     1      dot_prod(3,2)*xr(3,3+1) - dot_prod(3,2)*xr(3-1,3) -
     1      dot_prod(3,1)*xr(3+1,3+2) + 2.0*dot_prod(3-1,3)*xr(3,3+1) -
     1      (1.0/x_prod(3-1))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(dot_prod(3,1)*xr(3-1,3) -
     1      dot_prod(3-1,1)*xr(3,3+1) - dot_prod(3-1,2)*xr(3,3+1) +
     1      dot_prod(3-1,2)*xr(3-1,3)) -
     1      (1.0/x_prod(3))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(-dot_prod(3+1,1)*xr(3,3+1) +
     1      dot_prod(3,2)*xr(3+1,3+2))) 

      a3 = -coef2*(dot_prod(3-2,2)*xr(3,3+1) -
     1      dot_prod(3-2,2)*xr(3-1,3) + dot_prod(3-1,2)*xr(3-2,3-1) +
     1      dot_prod(3-1,1)*xr(3-2,3-1) - 2.0*dot_prod(3-2,3)*xr(3-1,3) -
     1      (1.0/x_prod(3-2))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3-2,1)*xr(3-1,3) -
     1      dot_prod(3-2,2)*xr(3-2,3-1)) -
     1      (1.0/x_prod(3-1))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3,1)*xr(3-1,3) -
     1      dot_prod(3-1,1)*xr(3,3+1) - dot_prod(3-1,2)*xr(3,3+1) +
     1      dot_prod(3-1,2)*xr(3-1,3))) 

      fta_x(3) = a1 + a2 + a3 
 
      a1 = -coef*(-dot_prod(3+1,2)*yr(3+1,3+2) +
     1      dot_prod(3+1,1)*yr(3+2,3+3) -
     1      (1.0/x_prod(3))*(dot_prod(3+1,2)*dot_prod(3,2) -
     1      dot_prod(3,3)*dot_prod(3+1,1))*(-dot_prod(3+1,1)*yr(3,3+1) +
     1      dot_prod(3,2)*yr(3+1,3+2))) 
      
      a2 = -coef1*(-dot_prod(3-1,2)*yr(3+1,3+2) +
     1      dot_prod(3,2)*yr(3,3+1) - dot_prod(3,2)*yr(3-1,3) -
     1      dot_prod(3,1)*yr(3+1,3+2) + 2.0*dot_prod(3-1,3)*yr(3,3+1) -
     1      (1.0/x_prod(3-1))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(dot_prod(3,1)*yr(3-1,3) -
     1      dot_prod(3-1,1)*yr(3,3+1) - dot_prod(3-1,2)*yr(3,3+1) +
     1      dot_prod(3-1,2)*yr(3-1,3)) -
     1      (1.0/x_prod(3))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(-dot_prod(3+1,1)*yr(3,3+1) +
     1      dot_prod(3,2)*yr(3+1,3+2))) 

      a3 = -coef2*(dot_prod(3-2,2)*yr(3,3+1) -
     1      dot_prod(3-2,2)*yr(3-1,3) + dot_prod(3-1,2)*yr(3-2,3-1) +
     1      dot_prod(3-1,1)*yr(3-2,3-1) - 2.0*dot_prod(3-2,3)*yr(3-1,3) -
     1      (1.0/x_prod(3-2))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3-2,1)*yr(3-1,3) -
     1      dot_prod(3-2,2)*yr(3-2,3-1)) -
     1      (1.0/x_prod(3-1))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3,1)*yr(3-1,3) -
     1      dot_prod(3-1,1)*yr(3,3+1) - dot_prod(3-1,2)*yr(3,3+1) +
     1      dot_prod(3-1,2)*yr(3-1,3)))

      fta_y(3) = a1 + a2 + a3 
 
      a1 = -coef*(-dot_prod(3+1,2)*zr(3+1,3+2) +
     1      dot_prod(3+1,1)*zr(3+2,3+3) -
     1      (1.0/x_prod(3))*(dot_prod(3+1,2)*dot_prod(3,2) -
     1      dot_prod(3,3)*dot_prod(3+1,1))*(-dot_prod(3+1,1)*zr(3,3+1) +
     1      dot_prod(3,2)*zr(3+1,3+2))) 

      a2 =  -coef1*(-dot_prod(3-1,2)*zr(3+1,3+2) +
     1      dot_prod(3,2)*zr(3,3+1) - dot_prod(3,2)*zr(3-1,3) -
     1      dot_prod(3,1)*zr(3+1,3+2) + 2.0*dot_prod(3-1,3)*zr(3,3+1) -
     1      (1.0/x_prod(3-1))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(dot_prod(3,1)*zr(3-1,3) -
     1      dot_prod(3-1,1)*zr(3,3+1) - dot_prod(3-1,2)*zr(3,3+1) +
     1      dot_prod(3-1,2)*zr(3-1,3)) -
     1      (1.0/x_prod(3))*(dot_prod(3,2)*dot_prod(3-1,2) -
     1      dot_prod(3-1,3)*dot_prod(3,1))*(-dot_prod(3+1,1)*zr(3,3+1) +
     1      dot_prod(3,2)*zr(3+1,3+2))) 

      a3 = -coef2*(dot_prod(3-2,2)*zr(3,3+1) -
     1      dot_prod(3-2,2)*zr(3-1,3) + dot_prod(3-1,2)*zr(3-2,3-1) +
     1      dot_prod(3-1,1)*zr(3-2,3-1) - 2.0*dot_prod(3-2,3)*zr(3-1,3) -
     1      (1.0/x_prod(3-2))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3-2,1)*zr(3-1,3) -
     1      dot_prod(3-2,2)*zr(3-2,3-1)) -
     1      (1.0/x_prod(3-1))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3,1)*zr(3-1,3) -
     1      dot_prod(3-1,1)*zr(3,3+1) - dot_prod(3-1,2)*zr(3,3+1) +
     1      dot_prod(3-1,2)*zr(3-1,3))) 

      fta_z(3) = a1 + a2 + a3 
C }}}
C particles 4 to n-3
C {{{
      do i = 4, n-3

         coef=DFAC(i+1)

         coef1=DFAC(i)

         coef2=DFAC(i-1)

         coef3=DFAC(i-2)

         a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1      dot_prod(i+1,1)*xr(i+2,i+3) -
     1      (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1      dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1      dot_prod(i,2)*xr(i+1,i+2))) 

         a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1      dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1      dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1      (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1      dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1      dot_prod(i-1,2)*xr(i-1,i)) -
     1      (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1      dot_prod(i,2)*xr(i+1,i+2))) 

         a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1      dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1      dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1      dot_prod(i-2,2)*xr(i-2,i-1)) -
     1      (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1  dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1      dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1      dot_prod(i-1,2)*xr(i-1,i))) 

         a4 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) -
     1      dot_prod(i-2,1)*xr(i-3,i-2) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1      dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1      dot_prod(i-2,2)*xr(i-2,i-1))) 

         fta_x(i) = a1 + a2 + a3 + a4 

         a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1      dot_prod(i+1,1)*yr(i+2,i+3) -
     1      (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1      dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1      dot_prod(i,2)*yr(i+1,i+2))) 

         a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1      dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1      dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1      (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1      dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1      dot_prod(i-1,2)*yr(i-1,i)) -
     1      (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1      dot_prod(i,2)*yr(i+1,i+2))) 

         a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1      dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1      dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1      dot_prod(i-2,2)*yr(i-2,i-1)) -
     1      (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1      dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1      dot_prod(i-1,2)*yr(i-1,i))) 

         a4 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) -
     1      dot_prod(i-2,1)*yr(i-3,i-2) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1      dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1      dot_prod(i-2,2)*yr(i-2,i-1))) 

         fta_y(i) = a1 + a2 + a3 + a4 

         a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1      dot_prod(i+1,1)*zr(i+2,i+3) -
     1      (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1      dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1      dot_prod(i,2)*zr(i+1,i+2))) 

         a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1      dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1      dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1      (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1      dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1      dot_prod(i-1,2)*zr(i-1,i)) -
     1      (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1      dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1      dot_prod(i,2)*zr(i+1,i+2))) 

         a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1      dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1      dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1      dot_prod(i-2,2)*zr(i-2,i-1)) -
     1      (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1      dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1      dot_prod(i-1,2)*zr(i-1,i))) 
      
         a4 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) -
     1      dot_prod(i-2,1)*zr(i-3,i-2) -
     1      (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1      dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1      dot_prod(i-2,2)*zr(i-2,i-1))) 

         fta_z(i) = a1 + a2 + a3 + a4 

      enddo
C }}}
C particle n-2
C {{{
!     i = n-2
      coef1=DFAC(n-2)

      coef2=DFAC(n-2-1)

      coef3=DFAC(n-2-2)

      a1 = -coef1*(-dot_prod(n-2-1,2)*xr(n-2+1,n-2+2) + 
     1      dot_prod(n-2,2)*xr(n-2,n-2+1) - dot_prod(n-2,2)*xr(n-2-1,n-2) -
     1      dot_prod(n-2,1)*xr(n-2+1,n-2+2) + 2.0*dot_prod(n-2-1,3)*xr(n-2,n-2+1) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*xr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*xr(n-2,n-2+1) - dot_prod(n-2-1,2)*xr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*xr(n-2-1,n-2)) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(-dot_prod(n-2+1,1)*xr(n-2,n-2+1) +
     1      dot_prod(n-2,2)*xr(n-2+1,n-2+2)))


      a2 = -coef2*(dot_prod(n-2-2,2)*xr(n-2,n-2+1) -
     1      dot_prod(n-2-2,2)*xr(n-2-1,n-2) + dot_prod(n-2-1,2)*xr(n-2-2,n-2-1) +
     1      dot_prod(n-2-1,1)*xr(n-2-2,n-2-1) - 2.0*dot_prod(n-2-2,3)*xr(n-2-1,n-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2-2,1)*xr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*xr(n-2-2,n-2-1)) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2,1)*xr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*xr(n-2,n-2+1) - dot_prod(n-2-1,2)*xr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*xr(n-2-1,n-2))) 
      
      a3 = -coef3*(dot_prod(n-2-3,2)*xr(n-2-2,n-2-1) -
     1      dot_prod(n-2-2,1)*xr(n-2-3,n-2-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-2,2)*dot_prod(n-2-3,2) -
     1      dot_prod(n-2-3,3)*dot_prod(n-2-2,1))*(dot_prod(n-2-2,1)*xr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*xr(n-2-2,n-2-1))) 

      fta_x(n-2) = a1 + a2 + a3 

      a1 =  -coef1*(-dot_prod(n-2-1,2)*yr(n-2+1,n-2+2) +  
     1      dot_prod(n-2,2)*yr(n-2,n-2+1) - dot_prod(n-2,2)*yr(n-2-1,n-2) -
     1      dot_prod(n-2,1)*yr(n-2+1,n-2+2) + 2.0*dot_prod(n-2-1,3)*yr(n-2,n-2+1) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*yr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*yr(n-2,n-2+1) - dot_prod(n-2-1,2)*yr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*yr(n-2-1,n-2)) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(-dot_prod(n-2+1,1)*yr(n-2,n-2+1) +
     1      dot_prod(n-2,2)*yr(n-2+1,n-2+2))) 

      a2 =  -coef2*(dot_prod(n-2-2,2)*yr(n-2,n-2+1) -
     1      dot_prod(n-2-2,2)*yr(n-2-1,n-2) + dot_prod(n-2-1,2)*yr(n-2-2,n-2-1) +
     1      dot_prod(n-2-1,1)*yr(n-2-2,n-2-1) - 2.0*dot_prod(n-2-2,3)*yr(n-2-1,n-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2-2,1)*yr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*yr(n-2-2,n-2-1)) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2,1)*yr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*yr(n-2,n-2+1) - dot_prod(n-2-1,2)*yr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*yr(n-2-1,n-2)))

      a3 = -coef3*(dot_prod(n-2-3,2)*yr(n-2-2,n-2-1) -
     1      dot_prod(n-2-2,1)*yr(n-2-3,n-2-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-2,2)*dot_prod(n-2-3,2) -
     1      dot_prod(n-2-3,3)*dot_prod(n-2-2,1))*(dot_prod(n-2-2,1)*yr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*yr(n-2-2,n-2-1))) 

      fta_y(n-2) = a1 + a2 + a3 
 
      a1 = -coef1*(-dot_prod(n-2-1,2)*zr(n-2+1,n-2+2) +  
     1      dot_prod(n-2,2)*zr(n-2,n-2+1) - dot_prod(n-2,2)*zr(n-2-1,n-2) -
     1      dot_prod(n-2,1)*zr(n-2+1,n-2+2) + 2.0*dot_prod(n-2-1,3)*zr(n-2,n-2+1) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*zr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*zr(n-2,n-2+1) - dot_prod(n-2-1,2)*zr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*zr(n-2-1,n-2)) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(-dot_prod(n-2+1,1)*zr(n-2,n-2+1) +
     1      dot_prod(n-2,2)*zr(n-2+1,n-2+2))) 

      a2 = -coef2*(dot_prod(n-2-2,2)*zr(n-2,n-2+1) -
     1      dot_prod(n-2-2,2)*zr(n-2-1,n-2) + dot_prod(n-2-1,2)*zr(n-2-2,n-2-1) +
     1      dot_prod(n-2-1,1)*zr(n-2-2,n-2-1) - 2.0*dot_prod(n-2-2,3)*zr(n-2-1,n-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2-2,1)*zr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*zr(n-2-2,n-2-1)) -
     1      (1.0/x_prod(n-2-1))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2,1)*zr(n-2-1,n-2) -
     1      dot_prod(n-2-1,1)*zr(n-2,n-2+1) - dot_prod(n-2-1,2)*zr(n-2,n-2+1) +
     1      dot_prod(n-2-1,2)*zr(n-2-1,n-2))) 

      a3 = -coef3*(dot_prod(n-2-3,2)*zr(n-2-2,n-2-1) -
     1      dot_prod(n-2-2,1)*zr(n-2-3,n-2-2) -
     1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-2,2)*dot_prod(n-2-3,2) -
     1      dot_prod(n-2-3,3)*dot_prod(n-2-2,1))*(dot_prod(n-2-2,1)*zr(n-2-1,n-2) -
     1      dot_prod(n-2-2,2)*zr(n-2-2,n-2-1))) 

      fta_z(n-2) = a1 + a2 + a3 
C }}}
C particle n-1
C {{{
!     i = n-1
      coef2=DFAC(n-1-1)

      coef3=DFAC(n-1-2)

      a1 = -coef2*(dot_prod(n-1-2,2)*xr(n-1,n-1+1) - 
     1      dot_prod(n-1-2,2)*xr(n-1-1,n-1) +
     1      dot_prod(n-1-1,2)*xr(n-1-2,n-1-1) +  dot_prod(n-1-1,1)*xr(n-1-2,n-1-1) -
     1      2.0*dot_prod(n-1-2,3)*xr(n-1-1,n-1) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1-2,1)*xr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*xr(n-1-2,n-1-1)) -
     1      (1.0/x_prod(n-1-1))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1,1)*xr(n-1-1,n-1) -
     1      dot_prod(n-1-1,1)*xr(n-1,n-1+1) - dot_prod(n-1-1,2)*xr(n-1,n-1+1) +
     1      dot_prod(n-1-1,2)*xr(n-1-1,n-1))) 

      a2 = -coef3*(dot_prod(n-1-3,2)*xr(n-1-2,n-1-1) - 
     1      dot_prod(n-1-2,1)*xr(n-1-3,n-1-2) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-2,2)*dot_prod(n-1-3,2) -
     1      dot_prod(n-1-3,3)*dot_prod(n-1-2,1))*(dot_prod(n-1-2,1)*xr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*xr(n-1-2,n-1-1))) 

      fta_x(n-1) = a1 + a2  

      a1 = -coef2*(dot_prod(n-1-2,2)*yr(n-1,n-1+1) - 
     1      dot_prod(n-1-2,2)*yr(n-1-1,n-1) +
     1      dot_prod(n-1-1,2)*yr(n-1-2,n-1-1) +  dot_prod(n-1-1,1)*yr(n-1-2,n-1-1) -
     1      2.0*dot_prod(n-1-2,3)*yr(n-1-1,n-1) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1-2,1)*yr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*yr(n-1-2,n-1-1)) -
     1      (1.0/x_prod(n-1-1))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1,1)*yr(n-1-1,n-1) -
     1  dot_prod(n-1-1,1)*yr(n-1,n-1+1) - dot_prod(n-1-1,2)*yr(n-1,n-1+1) +
     1      dot_prod(n-1-1,2)*yr(n-1-1,n-1))) 

      a2 = -coef3*(dot_prod(n-1-3,2)*yr(n-1-2,n-1-1) - 
     1      dot_prod(n-1-2,1)*yr(n-1-3,n-1-2) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-2,2)*dot_prod(n-1-3,2) -
     1      dot_prod(n-1-3,3)*dot_prod(n-1-2,1))*(dot_prod(n-1-2,1)*yr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*yr(n-1-2,n-1-1))) 

      fta_y(n-1) = a1 + a2  

      a1 = -coef2*(dot_prod(n-1-2,2)*zr(n-1,n-1+1) - 
     1      dot_prod(n-1-2,2)*zr(n-1-1,n-1) +
     1      dot_prod(n-1-1,2)*zr(n-1-2,n-1-1) +  dot_prod(n-1-1,1)*zr(n-1-2,n-1-1) -
     1      2.0*dot_prod(n-1-2,3)*zr(n-1-1,n-1) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1-2,1)*zr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*zr(n-1-2,n-1-1)) -
     1      (1.0/x_prod(n-1-1))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1,1)*zr(n-1-1,n-1) -
     1      dot_prod(n-1-1,1)*zr(n-1,n-1+1) - dot_prod(n-1-1,2)*zr(n-1,n-1+1) +
     1      dot_prod(n-1-1,2)*zr(n-1-1,n-1))) 

      a2 = -coef3*(dot_prod(n-1-3,2)*zr(n-1-2,n-1-1) - 
     1      dot_prod(n-1-2,1)*zr(n-1-3,n-1-2) -
     1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-2,2)*dot_prod(n-1-3,2) -
     1      dot_prod(n-1-3,3)*dot_prod(n-1-2,1))*(dot_prod(n-1-2,1)*zr(n-1-1,n-1) -
     1      dot_prod(n-1-2,2)*zr(n-1-2,n-1-1))) 

      fta_z(n-1) = a1 + a2 
C }}} 
C particle n
C {{{
!     i = n
      coef3=DFAC(n-2)

      fta_x(n) = -coef3*(dot_prod(n-3,2)*xr(n-2,n-1) 
     1      - dot_prod(n-2,1)*xr(n-3,n-2) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-3,2) -
     1      dot_prod(n-3,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*xr(n-1,n) -
     1      dot_prod(n-2,2)*xr(n-2,n-1))) 

      fta_y(n) = -coef3*(dot_prod(n-3,2)*yr(n-2,n-1) - 
     1      dot_prod(n-2,1)*yr(n-3,n-2) 
     1      - (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-3,2) 
     1      - dot_prod(n-3,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*yr(n-1,n) -
     1      dot_prod(n-2,2)*yr(n-2,n-1))) 

      fta_z(n) = -coef3*(dot_prod(n-3,2)*zr(n-2,n-1) - 
     1      dot_prod(n-2,1)*zr(n-3,n-2) -
     1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-3,2) -
     1      dot_prod(n-3,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*zr(n-1,n) -
     1      dot_prod(n-2,2)*zr(n-2,n-1))) 
C }}}
C Total up the gradients
C {{{
      do i = 1, n
!        IF (I.EQ.2) THEN
!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'i,fnbx,fbx,fbax,ftax=',i,fnb_x(i),fb_x(i),fba_x(i),fta_x(i)
!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'i,fnby,fby,fbay,ftay=',i,fnb_y(i),fb_y(i),fba_y(i),fta_y(i)
!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'i,fnbz,fbz,fbaz,ftaz=',i,fnb_z(i),fb_z(i),fba_z(i),fta_z(i)
!        ENDIF
         fx(i) = fnb_x(i) + fb_x(i) + fba_x(i) + fta_x(i) 
         fy(i) = fnb_y(i) + fb_y(i) + fba_y(i) + fta_y(i) 
         fz(i) = fnb_z(i) + fb_z(i) + fba_z(i) + fta_z(i) 
      enddo

      do i = 1, n
         j = (i-1)*3
         fq(j+1) = -fx(i)
         fq(j+2) = -fy(i)
         fq(j+3) = -fz(i)
      enddo
C }}}
      return
      end
C }}}

!> @brief  Fill the parameter arrays

      SUBROUTINE PARAM_ARRAY_BLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     &   LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, 
     &   HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, N)
C {{{
C DECLARATIONS {{{
      IMPLICIT NONE
      INTEGER N
      INTEGER NTYPE(N), J1, I, J, ICOUNT
      DOUBLE PRECISION LJREP_BLN(N,N), LJATT_BLN(N,N), A_BLN(N), B_BLN(N), C_BLN(N), D_BLN(N),
     &                 LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, 
     &                 HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN
      CHARACTER(LEN=1) BEADLETTER(N), BLNSSTRUCT(N)
C }}}
C AMINO ACID TYPES {{{
C B=1 L=2 N=3
C
      DO J1=1,N
         IF (BEADLETTER(J1).EQ.'B') THEN
            NTYPE(J1)=1
         ELSEIF (BEADLETTER(J1).EQ.'L') THEN
            NTYPE(J1)=2
         ELSEIF (BEADLETTER(J1).EQ.'N') THEN
            NTYPE(J1)=3
         ELSE
            PRINT '(A,A1)','ERROR IN PARAM_ARRAYBLN, UNRECOGNISED BEAD TYPE: ',BEADLETTER(J1)
            STOP
         ENDIF
      ENDDO

C }}}

C PARAMETERS FOR THE DIHEDRAL ANGLE POTENTIAL 
C THE END BONDS HAVE NO DIHEDRAL TERM, SO THE TOTAL NUMBER OF TERMS
C IS N-3 (E.G. 4 ATOMS, 1 DIHEDRAL). NON-ZERO TERMS FOR
C 2 TO 3, 3 TO 4, 4 TO 5, ... , N-3 TO N-2, N-2 TO N-1. THE
C H, E AND T PARAMETERS ARE DEFINED FOR THE FIRST BEAD OF EACH EDGE,
C I.E. FOR 2, 3, 4, ..., N-2.
C {{{

      A_BLN(1:N)=0.0D0
      B_BLN(1:N)=0.0D0
      C_BLN(1:N)=0.0D0
      D_BLN(1:N)=0.0D0
      DO I=1,N-3
         IF (BLNSSTRUCT(I).EQ.'H') THEN
            A_BLN(I+1)=HABLN
            B_BLN(I+1)=HBBLN
            C_BLN(I+1)=HCBLN
            D_BLN(I+1)=HDBLN
         ELSE IF (BLNSSTRUCT(I).EQ.'E') THEN
            A_BLN(I+1)=EABLN
            B_BLN(I+1)=EBBLN
            C_BLN(I+1)=ECBLN
            D_BLN(I+1)=EDBLN
         ELSE IF (BLNSSTRUCT(I).EQ.'T') THEN
            A_BLN(I+1)=TABLN
            B_BLN(I+1)=TBBLN
            C_BLN(I+1)=TCBLN
            D_BLN(I+1)=TDBLN
         ELSE
            PRINT '(A,A1)','ERROR IN PARAM_ARRAYBLN, UNRECOGNISED SS TYPE: ',BLNSSTRUCT(J1)
            STOP
         ENDIF
C        PRINT '(A,I6,A,A1,A,4F12.4)','I+1=',I+1,' SYMBOL=',BLNSSTRUCT(I),' A,B,C,D=',A_BLN(I+1),B_BLN(I+1),C_BLN(I+1),D_BLN(I+1)
      ENDDO
C }}}
C  PARAMETERS FOR THE L-J INTERACTION BETWEEN NON-BONDED PARTICLES {{{

      DO I = 1, N-1
         DO J = I+1, N
            IF (NTYPE(I) .EQ. 3 .OR. NTYPE(J) .EQ. 3) THEN
               LJREP_BLN(I,J) = LJREPNN
               LJATT_BLN(I,J) = LJATTNN
               LJREP_BLN(J,I) = LJREPNN
               LJATT_BLN(J,I) = LJATTNN
            ELSE IF (NTYPE(I) .EQ. 1 .AND. NTYPE(J) .EQ. 1) THEN
               LJREP_BLN(I,J) = LJREPBB
               LJATT_BLN(I,J) = LJATTBB
               LJREP_BLN(J,I) = LJREPBB
               LJATT_BLN(J,I) = LJATTBB
            ELSE
               LJREP_BLN(I,J) = LJREPLL
               LJATT_BLN(I,J) = LJATTLL
               LJREP_BLN(J,I) = LJREPLL
               LJATT_BLN(J,I) = LJATTLL
            ENDIF
         ENDDO
      ENDDO
C }}}

      RETURN
      END
C }}}
!------------------------------------------------------------
! Doxygen: G46MERDIFF {{{
!
!> @name G46MERDIFF
!
!> @brief Calculate the energy, gradient, and second derivatives matrix for the Go-like BLN model \n
!> @author John Rose
!>
!> A particle configuration and number of particles is passed to the subroutine and
!> the energy, gradient, and matrix of second derivatives is returned.
!>
!> \param N          number of particles
!> \param QO         array of cartesian particle coordinates
!> \param GRAD       array of gradients
!> \param ENERGY     energy
! }}}
!------------------------------------------------------------
        SUBROUTINE G46MERDIFF(QO, N, GRAD, ENERGY, GTEST)
C {{{ 
C declarations {{{
        USE MODHESS

        IMPLICIT NONE

        LOGICAL GTEST, STEST
        INTEGER NTYPE(46), N
        DOUBLE PRECISION QO(3*N), GRAD(3*N), ENERGY
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N),D_PARAM(N),
     1                   C_PARAM(N), RK_THETA, RK_R, EPSILON, SIGMA, THETA_0, DELTA, RMASS
        PARAMETER (RMASS = 40.0, EPSILON = 0.0100570)
        PARAMETER (SIGMA=3.4, DELTA=1.0D-6, THETA_0 = 1.8326)
        PARAMETER (RK_R = 20.0*0.0100570, RK_THETA = 20.0*0.0100570)
        DOUBLE PRECISION X(N), Y(N), Z(N), XR(N,N), YR(N,N), ZR(N,N),
     2                  DOT_PROD(N,3), X_PROD(N), BOND_ANGLE(N), TOR_ANGLE(N), RADII(N,N)
C }}}
C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(46),
C    1  d_param(n),c_param(n)

        STEST=.FALSE.

        call gparam_array(a_param,b_param,c_param,d_param,n)
        call calc_int_coords(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        call calc_energy(qo,energy,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)
        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
        call calc_gradient(qo,grad,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)

C commented section {{{
C        DIF=1.0D-4
C        DO J1=1,3*N
C           TEMP1=QO(J1)
C           QO(J1)=QO(J1)+DIF
C           call calc_int_coords(qo,n)
C           call calc_energy(qo,V1,n)
C           QO(J1)=QO(J1)-2.0D0*DIF
C           call calc_int_coords(qo,n)
C           call calc_energy(qo,V2,n)
C           tgrad(J1)=(V1-V2)/(2.0D0*DIF)
C           QO(J1)=TEMP1
C        ENDDO
C        call calc_int_coords(qo,n)

C       PRINT*,'Analytical/Numerical first derivatives:'
C       WRITE(*,'(3G20.10)') (GRAD(J1)/TGRAD(J1),J1=1,3*N)
C }}}

        IF (.NOT.STEST) RETURN
        call calc_dyn(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, bond_angle,tor_angle,
     1                            radii,ntype)

        return
        end
C }}}
C
C Doxygen: gparam_array {{{
C>
C> \brief Fill the parameter arrays which specify interaction potentials
C> \param N INTEGER  - number of particles
C> \param a_param \param b_param - LJ interaction between non-bonded particles 
C> \param c_param \param d_param - dihedral angle potential
C>
C}}}
        SUBROUTINE GPARAM_ARRAY(A_PARAM,B_PARAM,C_PARAM,D_PARAM,N)
C {{{
C Declarations {{{
        IMPLICIT NONE
        LOGICAL CONNECT(46,46)
        INTEGER J, ICOUNT, I, J2, J1, N
        DOUBLE PRECISION NTYPE(46), A_PARAM(N,N), B_PARAM(N,N)
        DOUBLE PRECISION C_PARAM(N), D_PARAM(N), EPSILON
        PARAMETER (EPSILON = 0.0100570)
C }}}
C Specify amino acid types by filling in the array ntype(:) {{{

        ntype(1) = 1
        ntype(2) = 1
        ntype(3) = 1
        ntype(4) = 1
        ntype(5) = 1
        ntype(6) = 1
        ntype(7) = 1
        ntype(8) = 1
        ntype(9) = 1
        ntype(10) = 3
        ntype(11) = 3
        ntype(12) = 3
        ntype(13) = 2
        ntype(14) = 1
        ntype(15) = 2
        ntype(16) = 1
        ntype(17) = 2
        ntype(18) = 1
        ntype(19) = 2
        ntype(20) = 1
        ntype(21) = 3
        ntype(22) = 3
        ntype(23) = 3
        ntype(24) = 1
        ntype(25) = 1
        ntype(26) = 1
        ntype(27) = 1
        ntype(28) = 1
        ntype(29) = 1
        ntype(30) = 1
        ntype(31) = 1
        ntype(32) = 1
        ntype(33) = 3
        ntype(34) = 3
        ntype(35) = 3
        ntype(36) = 2
        ntype(37) = 1
        ntype(38) = 2
        ntype(39) = 1
        ntype(40) = 2
        ntype(41) = 1
        ntype(42) = 2
        ntype(43) = 1
        ntype(44) = 2
        ntype(45) = 1
        ntype(46) = 2
     
C }}}
C Go-like model connectivities: fill in array CONNECT(:,:) {{{
C
        DO J1=1,46
           DO J2=J1,46
              CONNECT(J2,J1)=.FALSE.
           ENDDO
        ENDDO
        CONNECT(20, 1)=.TRUE.
        CONNECT(24, 1)=.TRUE.
        CONNECT(45, 1)=.TRUE.
        CONNECT(24, 2)=.TRUE.
        CONNECT(43, 2)=.TRUE.
        CONNECT(45, 2)=.TRUE.
        CONNECT(18, 3)=.TRUE.
        CONNECT(20, 3)=.TRUE.
        CONNECT(25, 3)=.TRUE.
        CONNECT(26, 3)=.TRUE.
        CONNECT(43, 3)=.TRUE.
        CONNECT(26, 4)=.TRUE.
        CONNECT(41, 4)=.TRUE.
        CONNECT(16, 5)=.TRUE.
        CONNECT(18, 5)=.TRUE.
        CONNECT(26, 5)=.TRUE.
        CONNECT(27, 5)=.TRUE.
        CONNECT(28, 5)=.TRUE.
        CONNECT(41, 5)=.TRUE.
        CONNECT(28, 6)=.TRUE.
        CONNECT(39, 6)=.TRUE.
        CONNECT(16, 7)=.TRUE.
        CONNECT(28, 7)=.TRUE.
        CONNECT(29, 7)=.TRUE.
        CONNECT(30, 7)=.TRUE.
        CONNECT(39, 7)=.TRUE.
        CONNECT(30, 8)=.TRUE.
        CONNECT(37, 8)=.TRUE.
        CONNECT(14, 9)=.TRUE.
        CONNECT(30, 9)=.TRUE.
        CONNECT(31, 9)=.TRUE.
        CONNECT(32, 9)=.TRUE.
        CONNECT(37, 9)=.TRUE.
        CONNECT(30, 14)=.TRUE.
        CONNECT(31, 14)=.TRUE.
        CONNECT(28, 16)=.TRUE.
        CONNECT(29, 16)=.TRUE.
        CONNECT(26, 18)=.TRUE.
        CONNECT(24, 20)=.TRUE.
        CONNECT(25, 20)=.TRUE.
        CONNECT(45, 24)=.TRUE.
        CONNECT(41, 26)=.TRUE.
        CONNECT(43, 26)=.TRUE.
        CONNECT(39, 28)=.TRUE.
        CONNECT(41, 28)=.TRUE.
        CONNECT(39, 30)=.TRUE.
        CONNECT(37, 32)=.TRUE.
        CONNECT(1, 20)=.TRUE.
        CONNECT(1, 24)=.TRUE.
        CONNECT(1, 45)=.TRUE.
        CONNECT(2, 24)=.TRUE.
        CONNECT(2, 43)=.TRUE.
        CONNECT(2, 45)=.TRUE.
        CONNECT(3, 18)=.TRUE.
        CONNECT(3, 20)=.TRUE.
        CONNECT(3, 25)=.TRUE.
        CONNECT(3, 26)=.TRUE.
        CONNECT(3, 43)=.TRUE.
        CONNECT(4, 26)=.TRUE.
        CONNECT(4, 41)=.TRUE.
        CONNECT(5, 16)=.TRUE.
        CONNECT(5, 18)=.TRUE.
        CONNECT(5, 26)=.TRUE.
        CONNECT(5, 27)=.TRUE.
        CONNECT(5, 28)=.TRUE.
        CONNECT(5, 41)=.TRUE.
        CONNECT(6, 28)=.TRUE.
        CONNECT(6, 39)=.TRUE.
        CONNECT(7, 16)=.TRUE.
        CONNECT(7, 28)=.TRUE.
        CONNECT(7, 29)=.TRUE.
        CONNECT(7, 30)=.TRUE.
        CONNECT(7, 39)=.TRUE.
        CONNECT(8, 30)=.TRUE.
        CONNECT(8, 37)=.TRUE.
        CONNECT(9, 14)=.TRUE.
        CONNECT(9, 30)=.TRUE.
        CONNECT(9, 31)=.TRUE.
        CONNECT(9, 32)=.TRUE.
        CONNECT(9, 37)=.TRUE.
        CONNECT(14, 30)=.TRUE.
        CONNECT(14, 31)=.TRUE.
        CONNECT(16, 28)=.TRUE.
        CONNECT(16, 29)=.TRUE.
        CONNECT(18, 26)=.TRUE.
        CONNECT(20, 24)=.TRUE.
        CONNECT(20, 25)=.TRUE.
        CONNECT(24, 45)=.TRUE.
        CONNECT(26, 41)=.TRUE.
        CONNECT(26, 43)=.TRUE.
        CONNECT(28, 39)=.TRUE.
        CONNECT(28, 41)=.TRUE.
        CONNECT(30, 39)=.TRUE.
        CONNECT(32, 37)=.TRUE.

C }}}
C Parameters for the dihedral angle potential: fill in arrays c_param(:,:), d_param(:,:) {{{

        do i = 1, n-3
        icount = 0

        do j = 0,3
        if(ntype(i+j) .eq. 3)then
        icount = icount + 1
        endif
        enddo

        if(icount .ge. 2)then
        c_param(i+1) = 0.0
        d_param(i+1) = 0.2*epsilon
        else
        c_param(i+1) = 1.2*epsilon
        d_param(i+1) = 1.2*epsilon
        endif

        icount = 0

        enddo
C }}}
C Parameters for the L-J interaction between non-bonded particles:
C arrays a_param(:,:), b_param(:,:)
C {{{

        do i = 1, n-1
           do j = i+1, n

           if (ntype(i) .eq. 3 .or. ntype(j) .eq. 3) then
             a_param(i,j) = 1.0*epsilon 
             b_param(i,j) = 0.0 
             a_param(j,i) = 1.0*epsilon 
             b_param(j,i) = 0.0
           elseif (ntype(i) .eq. 1 .and. ntype(j) .eq. 1)then
             a_param(i,j) =  epsilon
             a_param(j,i) =  epsilon
             IF (CONNECT(I,J)) THEN
                b_param(i,j) = -epsilon 
                b_param(j,i) = -epsilon
             ELSE
                b_param(i,j) = 0.0D0
                b_param(j,i) = 0.0D0
             ENDIF
           else
             a_param(i,j) = epsilon*2.0/3.0 
             b_param(i,j) = epsilon*2.0/3.0 
             a_param(j,i) = epsilon*2.0/3.0 
             b_param(j,i) = epsilon*2.0/3.0 
           endif
   
           enddo
        enddo
C }}}
        return
        end
C }}}
!------------------------------------------------------------
! Doxygen - P46MERDIFF {{{
!
!> @name P46MERDIFF
!
!> @brief
!> Calculate the energy, gradient, and second
!> derivatives for a given configuration of the 46 particle polymer chain.
!> A configuration and number of particles is passed to the subroutine and
!> the energy (ENERGY), gradient (GRAD), and 
!> the matrix of second derivatives are returned.
!
!> @author John Rose
!
! }}}
!------------------------------------------------------------

        SUBROUTINE P46MERDIFF(QO, N, GRAD, ENERGY, GTEST)
C {{{
        IMPLICIT NONE
        INTEGER NTYPE(46),N
        DOUBLE PRECISION RMASS, EPSILON,SIGMA,DELTA,THETA_0,RK_R,RK_THETA,QO(3*N),GRAD(3*N),ENERGY
        PARAMETER (RMASS = 40.0, EPSILON = 0.0100570)
        PARAMETER (SIGMA=3.4, DELTA=1.0D-6, THETA_0 = 1.8326)
        PARAMETER (RK_R = 20.0*0.0100570, RK_THETA = 20.0*0.0100570)
        LOGICAL GTEST, STEST
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1                  X(N), Y(N), Z(N), XR(N,N), YR(N,N), ZR(N,N),
     2                  DOT_PROD(N,3), X_PROD(N), BOND_ANGLE(N), TOR_ANGLE(N), RADII(N,N)

C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(46),
C    2  d_param(n),c_param(n),
C    3  x(n), y(n), z(n), 
C    4  xr(n,n), yr(n,n), zr(n,n), 
C    5  dot_prod(n,3), x_prod(n), 
C    6  bond_angle(n), stest, tor_angle(n), radii(n,n)

        STEST=.FALSE.

        CALL PARAM_ARRAY(A_PARAM,B_PARAM,C_PARAM,D_PARAM,N)
        CALL CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                       BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
        CALL CALC_ENERGY(QO,ENERGY,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                   BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
        CALL CALC_GRADIENT(QO,GRAD,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                     BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)

        IF (.NOT.STEST) RETURN
        CALL CALC_DYN(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)

        RETURN
        END
C }}}
C> Calculate the internal coordinates

        SUBROUTINE CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                             BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
C {{{
        IMPLICIT NONE
        INTEGER ntype(46),I,J,N
        DOUBLE PRECISION QO(3*N)
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1       x(n), y(n), z(n), xr(n,n), yr(n,n), zr(n,n),
     2       dot_prod(n,3), x_prod(n), bond_angle(n), tor_angle(n), radii(n,n),
     3       COS_THETA, COS_PHI

C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(46),
C    2  d_param(n),c_param(n),
C    3  x(n), y(n), z(n),
C    4  xr(n,n), yr(n,n), zr(n,n),
C    5  dot_prod(n,3), x_prod(n),
C    6  bond_angle(n), tor_angle(n), radii(n,n)

        do i = 1, n
        j = (i-1)*3
        x(i) = qo(j+1)
        y(i) = qo(j+2)
        z(i) = qo(j+3)
        enddo

C Inter-particle distances

        do i = 1, n-1
        do j = i+1, n
        xr(i,j) = x(j) - x(i)
        yr(i,j) = y(j) - y(i)
        zr(i,j) = z(j) - z(i)
        radii(i,j) = dsqrt(xr(i,j)*xr(i,j) + yr(i,j)*yr(i,j) 
     1  + zr(i,j)*zr(i,j))
        radii(j,i) = radii(i,j)
        enddo
        enddo

C Dot products between bond vectors

        do i = 1, n-3
        dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) +
     1  zr(i,i+1)*zr(i,i+1)

        dot_prod(i,2) = xr(i,i+1)*xr(i+1,i+2)+yr(i,i+1)*yr(i+1,i+2)+ 
     1  zr(i,i+1)*zr(i+1,i+2)

        dot_prod(i,3) = xr(i,i+1)*xr(i+2,i+3)+yr(i,i+1)*yr(i+2,i+3)+
     1  zr(i,i+1)*zr(i+2,i+3)
        enddo

        i = n-2
        dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) +
     1  zr(i,i+1)*zr(i,i+1)

        dot_prod(i,2) = xr(i,i+1)*xr(i+1,i+2)+yr(i,i+1)*yr(i+1,i+2)+
     1  zr(i,i+1)*zr(i+1,i+2)

        i = n-1
        dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) +
     1  zr(i,i+1)*zr(i,i+1)

C Cross-products between adjacent bond vectors

        do i = 1, n-2
        x_prod(i) = dot_prod(i,1)*dot_prod(i+1,1) -
     1               dot_prod(i,2)*dot_prod(i,2)   
        enddo

C Bond angles

        do i = 1, n-2
        cos_theta=-dot_prod(i,2)/(dsqrt(dot_prod(i,1)
     1  *dot_prod(i+1,1)))
        bond_angle(i+1) = dacos(cos_theta)
        enddo

C Torsional angles

        do i = 1, n-3
        cos_phi = (dot_prod(i,2)*dot_prod(i+1,2) -
     1  dot_prod(i,3)*dot_prod(i+1,1))/dsqrt(x_prod(i)*x_prod(i+1))
        IF (ABS(cos_phi).GT.1.0D0) cos_phi=cos_phi/abs(cos_phi)
        tor_angle(i+1) = dacos(cos_phi)
C       WRITE(*,'(A,I4,4F20.10)') 'i,tor_angle,cos_phi,dacos=',i,tor_angle(i+1),cos_phi,dacos(cos_phi)
        enddo

        return
        end

C }}} 
C> Calculate the energy
        SUBROUTINE CALC_ENERGY(QO,ENERGY,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                         BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
C {{{
        IMPLICIT NONE
        INTEGER I,J,N
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, THETA_0, RK_R, RK_THETA, RAD6, E_TANGLE,
     1                   QO(*), ENERGY, S6, E_NBOND, E_BOND, E_BANGLE
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        INTEGER ntype(46)
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1                  x(n), y(n), z(n), xr(n,n), yr(n,n), zr(n,n),
     2                  dot_prod(n,3), x_prod(n), bond_angle(n), tor_angle(n), radii(n,n)

C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(46),
C    2  d_param(n),c_param(n),
C    3  x(n), y(n), z(n),
C    4  xr(n,n), yr(n,n), zr(n,n),
C    5  dot_prod(n,3), x_prod(n),
C    6  bond_angle(n), tor_angle(n), radii(n,n)

        s6 = sigma*sigma*sigma*sigma*sigma*sigma
        e_nbond=0.0D0
        e_bond=0.0D0
        e_bangle=0.0D0
        e_tangle=0.0D0

        do i = 1, n-2
        do j = i+2, n

        rad6 = radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*
     1  radii(i,j)

        e_nbond = e_nbond + 4.0*((a_param(i,j)*s6*s6/(rad6*rad6)) + 
     1  (b_param(i,j)*s6/rad6))

        enddo
        enddo

        do i = 1, n-1

        e_bond = e_bond + 0.5*rk_r*(radii(i,i+1)-sigma)*
     1  (radii(i,i+1)-sigma)

        enddo

        
        do i = 2, n-1

        e_bangle = e_bangle + 0.5*rk_theta*(bond_angle(i)-theta_0)
     1  *(bond_angle(i)-theta_0)

        enddo

        do i = 2, n-2

        e_tangle = e_tangle + c_param(i)*(1.0 + cos(tor_angle(i))) 
     1  + d_param(i)*(1.0 + cos(3.0*tor_angle(i)))

        enddo

        energy = e_nbond + e_bond + e_bangle + e_tangle
C       WRITE(*,'(A,4F20.10)') 'nbond,bond,bangle,tangle=',e_nbond,e_bond,e_bangle,e_tangle

        return
        end
C }}}
C Calculate the gradients

        SUBROUTINE CALC_GRADIENT(QO,FQ,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, 
     1                           BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
C {{{
C Declarations {{{
        IMPLICIT NONE
        INTEGER I, J, N
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, RK_R, RK_THETA,
     1                   A4, COEF, COEF1, COEF2, COEF3, A3, DEN2, A2, A1, DEN1, RNUM, DEN,
     2                   RVAR, FZZ, FYY, FXX, DF, RAD14, RAD7, S6, QO(3*n), theta_0
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6,theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        INTEGER ntype(46)
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1                  x(n),y(n),z(n),xr(n,n), yr(n,n), zr(n,n),
     2     dot_prod(n,3), x_prod(n),bond_angle(n),tor_angle(n), radii(n,n)
        DOUBLE PRECISION FNB_X(N),FNB_Y(N),FNB_Z(N),
     1  fb_x(n),fb_y(n),fb_z(n)
        DOUBLE PRECISION FBA_X(N),FBA_Y(N),FBA_Z(N), FX(N),FY(N),FZ(N)
        DOUBLE PRECISION FTA_X(N),FTA_Y(N),FTA_Z(N), FQ(3*N)
C }}}
C       common/work/a_param(n,n),
C    1  b_param(n,n),ntype(46),
C    2  d_param(n),c_param(n),
C    3  x(n), y(n), z(n),
C    4  xr(n,n), yr(n,n), zr(n,n),
C    5  dot_prod(n,3), x_prod(n),
C    6  bond_angle(n), tor_angle(n), radii(n,n)

        s6 = sigma*sigma*sigma*sigma*sigma*sigma

C Gradients of potential
C {{{

        do i = 1,n

        fnb_x(i) = 0.0  
        fnb_y(i) = 0.0 
        fnb_z(i) = 0.0 

        fb_x(i)  = 0.0 
        fb_y(i)  = 0.0 
        fb_z(i)  = 0.0 

        fba_x(i) = 0.0 
        fba_y(i) = 0.0 
        fba_z(i) = 0.0 

        fta_x(i) = 0.0 
        fta_y(i) = 0.0 
        fta_z(i) = 0.0 

        fx(i)= 0.0 
        fy(i)= 0.0 
        fz(i)= 0.0 

        enddo
C }}}
C ..... Non-bonded interaction forces ..... {{{

        do i = 1, n-2
        do j = i+2, n

        rad7 = radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*
     1        radii(i,j)*radii(i,j)*radii(i,j)   
        rad14 = rad7*rad7 

        df = -24.0*((2.0*a_param(i,j)*s6*s6/rad14) + 
     1              (b_param(i,j)*s6/(rad7*radii(i,j))))

        fxx = df*xr(i,j) 
        fyy = df*yr(i,j) 
        fzz = df*zr(i,j) 

        fnb_x(i) = fxx + fnb_x(i)
        fnb_y(i) = fyy + fnb_y(i)
        fnb_z(i) = fzz + fnb_z(i)

        fnb_x(j) = -fxx + fnb_x(j)
        fnb_y(j) = -fyy + fnb_y(j)
        fnb_z(j) = -fzz + fnb_z(j)

        enddo
        enddo
C }}}
C ... Bond interaction forces ... {{{

        do i = 1, n-1

        rvar = sigma/radii(i,i+1) 

        df = rk_r*(1.0 - rvar) 
        fxx = df*xr(i,i+1) 
        fyy = df*yr(i,i+1) 
        fzz = df*zr(i,i+1) 

        fb_x(i) = fxx + fb_x(i)
        fb_y(i) = fyy + fb_y(i)
        fb_z(i) = fzz + fb_z(i)

        fb_x(i+1) = -fxx + fb_x(i+1)
        fb_y(i+1) = -fyy + fb_y(i+1)
        fb_z(i+1) = -fzz + fb_z(i+1)

        enddo
C }}}
C bond angle forces  particle 1
C particles 1,2,n-1, and n done outside of the loop
C {{{
        i = 1
        den = dsin(bond_angle(i+1))
     1        *dsqrt(dot_prod(i+1,1)*dot_prod(i,1))
        rnum = rk_theta*(bond_angle(i+1) - theta_0)

        fba_x(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*xr(i,i+1) -
     1        xr(i+1,i+2))/den

        fba_y(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*yr(i,i+1) -
     1        yr(i+1,i+2))/den

        fba_z(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*zr(i,i+1) -
     1        zr(i+1,i+2))/den

C particle 2

        i = 2
        den = dsin(bond_angle(i))
     1        *dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i+1))*dsqrt(dot_prod(i+1,1)
     1         *dot_prod(i,1))

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

        fba_x(i) = a1 + a2 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

        fba_y(i) = a1 + a2 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

        fba_z(i) = a1 + a2 

C particles 3 thru n-2 

        do i = 3, n-2

        den = dsin(bond_angle(i))*
     1              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i+1))*
     1               dsqrt(dot_prod(i+1,1)*dot_prod(i,1))
        den2 = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     1         *dot_prod(i-1,1))

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den2

        fba_x(i) = a1 + a2 + a3 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den2

        fba_y(i) = a1 + a2 + a3 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den2

        fba_z(i) = a1 + a2 + a3 

        enddo

C particle n-1 

        i = n-1
        den = dsin(bond_angle(i))*
     1              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     1         *dot_prod(i-1,1))

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den1

        fba_x(i) = a1 + a2

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den1

        fba_y(i) = a1 + a2

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den1

        fba_z(i) = a1 + a2

C particle n

        i = n
        den = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     1        *dot_prod(i-1,1))

        fba_x(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*xr(i-1,i) 
     1        - xr(i-2,i-1))/den

        fba_y(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*yr(i-1,i) 
     1        - yr(i-2,i-1))/den

        fba_z(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*zr(i-1,i) 
     1        - zr(i-2,i-1))/den
C }}}
C Torsional angle forces
C particles 1, 2, 3, n-2, n-1, and n are done outside of the loop
C particle 1

        i = 1
             coef =(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1         *dcos(tor_angle(i+1))-3.0))
     1  *(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))  

        fta_x(i) = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1         dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        fta_y(i) = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        fta_z(i) = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 


C particle 2

        i = 2
        coef =(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1  *dcos(tor_angle(i+1)) - 3.0))
     1        *(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))  

             coef1 = (c_param(i) + d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))  

        a1 =  -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1        dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        fta_x(i) = a1 + a2 

        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2)))

        fta_y(i) = a1 + a2 
        
        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        fta_z(i) = a1 + a2 

C particle 3

        i = 3
        coef=(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1  *dcos(tor_angle(i+1)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))  

        coef1=(c_param(i)+d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))  

        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1        dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 

        fta_x(i) = a1 + a2 + a3 
 
        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 
        
        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)))

        fta_y(i) = a1 + a2 + a3 
 
        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 =  -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 

        fta_z(i) = a1 + a2 + a3 

C particles 4 to n-3

        do i = 4, n-3

        coef = (c_param(i+1) + d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1  *dcos(tor_angle(i+1)) - 3.0))*(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))

        coef1 = (c_param(i) + d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) -3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))  

        coef2 = (c_param(i-1) + d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1  3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        coef3 = (c_param(i-2) + d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1  3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1        dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1  dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 

        a4 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) -
     1        dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_x(i) = a1 + a2 + a3 + a4 

        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i))) 

        a4 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) -
     1        dot_prod(i-2,1)*yr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_y(i) = a1 + a2 + a3 + a4 

        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 
        
        a4 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) -
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

        fta_z(i) = a1 + a2 + a3 + a4 

        enddo

C particle n-2

        i = n-2
        coef1=(c_param(i)+d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))  

        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1  3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        a1 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) + 
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2)))


        a2 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 
        
        a3 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) -
     1        dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_x(i) = a1 + a2 + a3 

        a1 =  -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +  
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a2 =  -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)))

        a3 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) -
     1        dot_prod(i-2,1)*yr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_y(i) = a1 + a2 + a3 
 
        a1 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +  
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 

        a3 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) -
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

        fta_z(i) = a1 + a2 + a3 

C particle n-1

        i = n-1
        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        a1 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) - 
     1        dot_prod(i-2,2)*xr(i-1,i) +
     1        dot_prod(i-1,2)*xr(i-2,i-1) +  dot_prod(i-1,1)*xr(i-2,i-1) -
     1        2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 

        a2 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) - 
     1        dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_x(i) = a1 + a2  

        a1 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) - 
     1        dot_prod(i-2,2)*yr(i-1,i) +
     1        dot_prod(i-1,2)*yr(i-2,i-1) +  dot_prod(i-1,1)*yr(i-2,i-1) -
     1        2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1  dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i))) 

        a2 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - 
     1        dot_prod(i-2,1)*yr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_y(i) = a1 + a2  

        a1 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) - 
     1        dot_prod(i-2,2)*zr(i-1,i) +
     1        dot_prod(i-1,2)*zr(i-2,i-1) +  dot_prod(i-1,1)*zr(i-2,i-1) -
     1        2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 

        a2 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - 
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

        fta_z(i) = a1 + a2 
 
C particle n

        i = n
        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        fta_x(i) = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) 
     1        - dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_y(i) = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - 
     1        dot_prod(i-2,1)*yr(i-3,i-2) 
     1        - (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) 
     1        - dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_z(i) = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - 
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

C Total up the gradients

        do i = 1, n
        fx(i) = fnb_x(i) + fb_x(i) + fba_x(i) + fta_x(i) 
        fy(i) = fnb_y(i) + fb_y(i) + fba_y(i) + fta_y(i) 
        fz(i) = fnb_z(i) + fb_z(i) + fba_z(i) + fta_z(i) 
        enddo

        do i = 1, n
        j = (i-1)*3
        fq(j+1) = -fx(i)
        fq(j+2) = -fy(i)
        fq(j+3) = -fz(i)
        enddo

        return
        end
C }}}
C> Calculate the second derivative matrix (two-sided numerical approach)

        SUBROUTINE CALC_DYN(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, 
     1                      BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
C {{{
C Declarations {{{
        USE MODHESS
        IMPLICIT NONE
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, THETA_0, RK_R, RK_THETA
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4 ,delta=1.0d-4, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        INTEGER ntype(46), N, I, J
        DOUBLE PRECISION QO(3*N), FQ1(3*N), 
     1  fq2(3*n)
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1                  x(n), y(n),z(n),xr(n,n),yr(n,n),zr(n,n),
     2                  dot_prod(n,3), x_prod(n), bond_angle(n), tor_angle(n), 
     3                  radii(n,n)
C }}}
C Fill in the Hessian matrix
C {{{

        do j = 1, 3*n

        qo(j) = qo(j) + delta
        call calc_int_coords(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod,
     1                       bond_angle,tor_angle,radii,ntype)
        call calc_gradient(qo,fq2,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod,
     1                     bond_angle,tor_angle,radii,ntype)
        qo(j) = qo(j) - 2.0*delta
        call calc_int_coords(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod,
     1                       bond_angle,tor_angle,radii,ntype)
        call calc_gradient(qo,fq1,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod,
     1                     bond_angle,tor_angle,radii,ntype)
        qo(j) = qo(j) + delta

        do i = j, 3*n

        HESS(i,j) = (fq2(i) -  fq1(i))/(2.0*delta)
        HESS(j,i) = HESS(i,j)

        enddo
        enddo

C }}}
        return
        end
C }}}
C> Fill the parameter arrays

        SUBROUTINE PARAM_ARRAY(A_PARAM,B_PARAM,C_PARAM,D_PARAM,N)
C {{{
C Declarations {{{
        implicit NONE
        INTEGER ntype(46), N, ICOUNT, J, I
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), EPSILON
        DOUBLE PRECISION C_PARAM(N), D_PARAM(N)
        parameter (epsilon = 0.0100570D0)
C }}}
C Firstly, specify amino acid types by filling in array ntype(:)
C 1 -> Hydrophobic (B)
C 2 -> Hydrophilic (L)
C 3 -> Neutral (N)
C {{{
        ntype(1) = 1
        ntype(2) = 1
        ntype(3) = 1
        ntype(4) = 1
        ntype(5) = 1
        ntype(6) = 1
        ntype(7) = 1
        ntype(8) = 1
        ntype(9) = 1
        ntype(10) = 3
        ntype(11) = 3
        ntype(12) = 3
        ntype(13) = 2
        ntype(14) = 1
        ntype(15) = 2
        ntype(16) = 1
        ntype(17) = 2
        ntype(18) = 1
        ntype(19) = 2
        ntype(20) = 1
        ntype(21) = 3
        ntype(22) = 3
        ntype(23) = 3
        ntype(24) = 1
        ntype(25) = 1
        ntype(26) = 1
        ntype(27) = 1
        ntype(28) = 1
        ntype(29) = 1
        ntype(30) = 1
        ntype(31) = 1
        ntype(32) = 1
        ntype(33) = 3
        ntype(34) = 3
        ntype(35) = 3
        ntype(36) = 2
        ntype(37) = 1
        ntype(38) = 2
        ntype(39) = 1
        ntype(40) = 2
        ntype(41) = 1
        ntype(42) = 2
        ntype(43) = 1
        ntype(44) = 2
        ntype(45) = 1
        ntype(46) = 2
C }}}
C Specify parameters for the dihedral angle potential.
C       These are stored in arrays 
C       c_param(:), d_param(:) 
C {{{

        do i = 1, n-3
        icount = 0

        do j = 0,3
        if(ntype(i+j) .eq. 3)then
        icount = icount + 1
        endif
        enddo

        if(icount .ge. 2)then
        c_param(i+1) = 0.0
        d_param(i+1) = 0.2*epsilon
        else
        c_param(i+1) = 1.2*epsilon
        d_param(i+1) = 1.2*epsilon
        endif

        icount = 0

        enddo
!}}}
C Specify parameters for the L-J interaction between non-bonded particles.
C       These are stored in arrays 
C       a_param(:,:), b_param(:,:) 
C {{{

        do i = 1, n-1
        do j = i+1, n

        if (ntype(i) .eq. 3 .or. ntype(j) .eq. 3)then
        a_param(i,j) = 1.0*epsilon 
        b_param(i,j) = 0.0 
        a_param(j,i) = 1.0*epsilon 
        b_param(j,i) = 0.0

        elseif (ntype(i) .eq. 1 .and. ntype(j) .eq. 1)then
        a_param(i,j) =  epsilon
        b_param(i,j) = -epsilon 
        a_param(j,i) =  epsilon
        b_param(j,i) = -epsilon
        
        else

        a_param(i,j) = epsilon*2.0/3.0 
        b_param(i,j) = epsilon*2.0/3.0 
        a_param(j,i) = epsilon*2.0/3.0 
        b_param(j,i) = epsilon*2.0/3.0 

        endif

        enddo
        enddo
!}}}
        return
        end
C }}}

        ENDMODULE
