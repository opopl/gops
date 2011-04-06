     MODULE BLN 

     IMPLICIT NONE

     CONTAINS

!! generic BLN - commented {{{
!!
!!> @brief Calculate the energy and gradient for a given configuration of a BLN polymer chain.
!!>
!!> @param[out] ENERGY
!!> @param[in] GRADT
!!

      !SUBROUTINE BLN(N,QO,GRAD,ENERGY,GRADT)
!!{{{
      !USE COMMONS

!include "blnvars.inc.f90"

!!
!! Without these initialisations the NAG compiler fills in random numbers for
!! unassigned elements with optimisation turned on.
!!
      !ANG=0.0D0
      !COSTOR=0.0D0
      !DFAC=0.0D0
      !SINBOND=0.0D0
      !DOT_PROD=0.0D0
      !X_PROD=0.0D0
      !RADII=0.0D0
      !DR=0.0D0

      !CALL CALC_INT_COORDS_BLN(N,QO,R,DR,DOT_PROD,X_PROD,ANG,RADII,COSTOR,SINBOND,A,DFAC)

      !CALL CALC_ENERGY_BLN(N,QO,ENERGY,LJREP,LJATT,A,R,DR,DOT_PROD,X_PROD,ANG,RADII,&
                !RK_R,RK_THETA,COSTOR)

      !IF (.NOT.GRADT) RETURN
 
      !CALL CALC_GRADIENT_BLN(N,QO,GRAD,LJREP,LJATT,A,R,DR,&
         !DOT_PROD,X_PROD,ANG,RADII,RK_R,RK_THETA,COSTOR,DFAC,SINBOND)

      !RETURN
      !END
!! }}}

!!> @brief Calculate the internal coordinates of a BLN chain

      !SUBROUTINE CALC_INT_COORDS_BLN(N,QO,R,DR,DOT_PROD,X_PROD,ANG,RADII, 
     !&                              COSTOR,SINBOND,A,DFAC)
!! {{{
!include "blnvars.inc.f90"

      !DO I = 1, N
         !J = (I-1)*3
         !DO K=1,3
                !R(I,1) = QO((I-1)*3+K)
         !ENDDO
      !ENDDO
!!
!! Inter-particle distances {{{

      !DO I = 1, N-1
         !DO J = I+1, N
            !DR(I,J,1:3) = R(J,1:3) - R(I,1:3)
            !RADII(I,J) = SQRT(SUM(DR(I,J,1:3)**2))
            !RADII(J,I) = RADII(I,J)
         !ENDDO
      !ENDDO
!! }}}

!! Dot products between bond vectors {{{

      !DO I = 1, N-1
        !IF ( I .LE. N-3 ) THEN JMAX=I+2
        !IF ( I .EQ. N-2 ) THEN JMAX=I+1
        !IF ( I .LE. N-1 ) THEN JMAX=I
        !DO J=I,JMAX
         !! K=1..3   FOR I=1..N-3
         !! K=1..2   FOR I=N-2
         !! K=1      FOR I=N-1
         !K=J-I+1 
         !DOT_PROD(I,K) = SUM(DR(I,I+1,1:3)*DR(J,J+1,1:3))
        !ENDDO
      !ENDDO

!! }}}

!! Cross-products between adjacent bond vectors {{{

      !DO I = 1, N-2
         !X_PROD(I) = DOT_PROD(I,1)*DOT_PROD(I+1,1) - DOT_PROD(I,2)*DOT_PROD(I,2)   
      !ENDDO
!! }}}
!! Bond angles {{{

      !DO I = 1, N-2
         !COS_THETA=-DOT_PROD(I,2)/(SQRT(DOT_PROD(I,1)*DOT_PROD(I+1,1)))
         !ANG(I+1,1) = DACOS(COS_THETA)
         !SINBOND(I+1)=SIN(ANG(I+1,1))*SQRT(DOT_PROD(I,1)*DOT_PROD(I+1,1))
      !ENDDO
!! }}}
!! Torsional angles {{{

      !DO I = 1, N-3
         !COS_PHI = DOT_PROD(I,2)*DOT_PROD(I+1,2) - DOT_PROD(I,3)*DOT_PROD(I+1,1)
         !COS_PHI = COS_PHI/SQRT(X_PROD(I)*X_PROD(I+1))
         !IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=SIGN(1.0D0,COS_PHI)
         !ANG(I+1,2) = DACOS(COS_PHI)
!! }}}
!!
!! ANG is returned in the range 0 to Pi.
!! dummy should take the opposite sign from the dihedral angle. Negative
!! values of the dihedral should be subtracted from 2*pi.
!! This is only necessary when the potential contains cos(phi+pi/4) terms because
!! the gradient is discontinuous when phi goes through pi for such terms if phi
!! is restricted to 0 < phi < pi.
!! {{{
         !DUMMY=DR(I+2,I+3,1)*(-DR(I,I+1,2)*DR(I+1,I+2,3)+DR(I+1,I+2,2)*DR(I,I+1,3))+
     !&         DR(I+2,I+3,2)*(-DR(I+1,I+2,1)*DR(I,I+1,3)+DR(I,I+1,1)*DR(I+1,I+2,3))+
     !&         DR(I+2,I+3,3)*(-DR(I,I+1,1)*DR(I+1,I+2,2)+DR(I+1,I+2,1)*DR(I,I+1,2))

         !IF (DUMMY.GT.0.0D0) ANG(I+1,2)=TWOPI-ANG(I+1,2)
         !COSTOR(I+1)=COS(ANG(I+1,2))
!! }}}
!!  This is an ugly hack to prevent division by zero. There will be a loss of precision
!!  if a dihedral is 0, PI, 2*PI, if D_BLN is non-zero.
!! {{{
         !IF (TAN(ANG(I+1,2)).EQ.0.0D0) THEN
            !PRINT '(A,I8,A,G20.10)','WARNING in BLN, dihedral angle ',i+1,' is ',ANG(I+1,2)
            !ANG(I+1,2)=ANG(I+1,2)+1.0D-10
            !PRINT '(A,G20.10)','WARNING in BLN, TAN perturbed angle=',ANG(i+1,2)
         !ENDIF
         !DUMMY2=TAN(ANG(i+1,2))
         !DFAC(I+1)=(A(I+1,1)+D_BLN(I+1)*( 1.0D0+1.0D0/DUMMY2 )*0.7071067811865475244D0-B_BLN(I+1)
     !&              +C_BLN(I+1)*(12.0*COSTOR(I+1)**2-3.0))/SQRT(X_PROD(I+1)*X_PROD(I))
!! }}}
      !ENDDO

      !RETURN
      !END
!! }}}

!!> @brief Calculate the energy of a BLN chain

      !SUBROUTINE CALC_ENERGY_BLN(N,QO,ENERGY,LJREP,LJATT,A,ANG,RADII,RK_R,RK_THETA,COSTOR)
!! {{{
!include "blnvars.inc.f90"

      !E= 0.0D0  
      !DO I = 1, N-2
         !DO J = I+2, N
            !RAD(6) = RADII(I,J)**6; RAD(12)=RAD(6)**2
            !E(1)=E(1)+LJREP(I,J)/RAD(12) + LJATT(I,J)/RAD(6)
         !ENDDO
      !ENDDO
      !E(1)=E(1)*4.0D0

      !DO I = 1, N-1
         !E(2) = E(2) + (RADII(I,I+1)-1.0D0)**2
      !ENDDO
      !E(2)=E(2)*RK_R/2.0D0

      !DO I = 2, N-1
         !E(3) = E(3) + (ANG(I,1)-THETA_0)**2
      !ENDDO
      !E(3)=E(3)*RK_THETA/2.0D0

      !DO I = 2, N-2
         !E(4) = E(4) + A(I,1)*(1.0D0 + COSTOR(I)) + A(I,3)*(1.0 + COS(3.0*ANG(I,2))) &
                             !+ A(I,2)*(1.0D0 - COSTOR(I)) + A(I,4)*(1.0 + COS(ANG(I,2)+PI4))
      !ENDDO

      !ENERGY = SUM(E)

      !RETURN
      !END
!! }}}

!!> @brief Calculate the gradients

      !SUBROUTINE CALC_GRADIENT_BLN(N,QO,FQ,LJREP,LJATT,A,R,DR,DOT_PROD,X_PROD,
     !&                            ANG,RADII,RK_R,RK_THETA,COSTOR,DFAC,SINBOND)
!! {{{
!include "blnvars.inc.f90" 
!!
!! Gradients of potential
!! {{{

!FNB= 0.0D0  
!FB= 0.0D0  
!FBA= 0.0D0  
!FTA= 0.0D0  
!F= 0.0D0  

!! }}}
!! ..... Non-bonded interaction forces ..... 
!! {{{
      !DO I = 1, N-2
         !DO J = I+2, N
            !RAD(7)=RADII(I,J)**7
            !RAD(8)=RADII(I,J)*RAD(7)
            !RAD(14) = RAD(7)**2

            !DF = 2.0*LJREP(I,J)/RAD(14) + LJATT(I,J)/RAD(8)
            !DF = -24.0*DF

            !F(1:3) = DF*DR(I,J,1:3) 

            !FNB(I,1:3) = F(1:3) + FNB(I,1:3)
            !FNB(J,1:3) = -F(1:3) + FNB(J,1:3)
         !ENDDO
      !ENDDO
!! }}}
!! ... Bond interaction forces ... 
!! {{{
      !DO I = 1, N-1
         !RVAR = 1.0D0/RADII(I,I+1) 
         !DF = RK_R*(1.0 - RVAR) 
         !F(1:3) = DF*DR(I,I+1,1:3) 
         !FB(I,1:3) = F(1:3) + FB(I,1:3)
         !FB(I+1,1:3) = -F(1:3) + FB(I+1,1:3)
      !ENDDO
!! }}}
!! bond angle forces  particle 1
!! particles 1,2,n-1, and n done outside of the loop
!! {{{
!!     i = 1
      !DEN = SINBOND(1+1)
      !RNUM = RK_THETA*(ANG(1+1,1) - THETA_0)
      !FBA(1,1:3) = -RNUM*((DOT_PROD(1,2)/DOT_PROD(1,1))*DR(1,1+1,1:3) - DR(1+1,1+2,1:3))/DEN
      
!! }}}
!! particle 2
!! {{{
!!     i = 2
      !DEN = SINBOND(2)
      !DEN1 = SINBOND(3)

      !A1 = -RK_THETA*(ANG(2,1) - THETA_0)*( (DOT_PROD(2-1,2)/
     !1  DOT_PROD(2,1))*DR(2,2+1,1) - (DOT_PROD(2-1,2)/DOT_PROD(2-1,1))
     !1      *DR(2-1,2,1) + DR(2,2+1,1) - DR(2-1,2,1))/DEN

      !A2 = -RK_THETA*(ANG(2+1,1) - THETA_0)*((DOT_PROD(2,2)/DOT_PROD(2,1))*DR(2,2+1,1) - DR(2+1,2+2,1))/DEN1

      !FBA_X(2) = A1 + A2 

      !A1 = -RK_THETA*(ANG(2,1) - THETA_0)*( (DOT_PROD(2-1,2)/
     !1  DOT_PROD(2,1))*DR(2,2+1,2) - (DOT_PROD(2-1,2)/DOT_PROD(2-1,1))
     !1      *DR(2-1,2,2) + DR(2,2+1,2) - DR(2-1,2,2))/DEN

      !A2 = -RK_THETA*(ANG(2+1,1) - THETA_0)*((DOT_PROD(2,2)/DOT_PROD(2,1))*DR(2,2+1,2) - DR(2+1,2+2,2))/DEN1

      !FBA_Y(2) = A1 + A2 

      !A1 = -RK_THETA*(ANG(2,1) - THETA_0)*( (DOT_PROD(2-1,2)/
     !1  DOT_PROD(2,1))*DR(2,2+1,3) - (DOT_PROD(2-1,2)/DOT_PROD(2-1,1))
     !1      *DR(2-1,2,3) + DR(2,2+1,3) - DR(2-1,2,3))/DEN

      !A2 = -RK_THETA*(ANG(2+1,1) - THETA_0)*((DOT_PROD(2,2)/DOT_PROD(2,1))*DR(2,2+1,3) - DR(2+1,2+2,3))/DEN1

      !FBA_Z(2) = A1 + A2 

!! }}}
!! particles 3 thru n-2 
!! {{{
      !do i = 3, n-2

         !den = sinbond(i)
         !den1 = sinbond(i+1)
         !den2 = sinbond(i-1)

         !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1      *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

         !a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

         !a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den2

         !fba_x(i) = a1 + a2 + a3 

         !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1   dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1      *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

         !a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     !1      dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

         !a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     !1      dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den2

         !fba_y(i) = a1 + a2 + a3 

         !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1   dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1      *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

         !a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     !1      dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

         !a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     !1      dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den2

         !fba_z(i) = a1 + a2 + a3 

      !enddo
!! }}}
!! particle n-1 
!! {{{
!!     i = n-1
      !den = sinbond(n-1)
      !den1 = sinbond(n-1-1)

      !a1 = -rk_theta*(bond_angle(n-1) - theta_0)*( (dot_prod(n-1-1,2)/
     !1   dot_prod(n-1,1))*xr(n-1,n-1+1) - (dot_prod(n-1-1,2)/dot_prod(n-1-1,1))
     !1      *xr(n-1-1,n-1) + xr(n-1,n-1+1) - xr(n-1-1,n-1))/den

      !a2 = rk_theta*(bond_angle(n-1-1) - theta_0)*((dot_prod(n-1-2,2)/
     !1      dot_prod(n-1-1,1))*xr(n-1-1,n-1) - xr(n-1-2,n-1-1))/den1

      !fba_x(n-1) = a1 + a2

      !a1 = -rk_theta*(bond_angle(n-1) - theta_0)*( (dot_prod(n-1-1,2)/
     !1   dot_prod(n-1,1))*yr(n-1,n-1+1) - (dot_prod(n-1-1,2)/dot_prod(n-1-1,1))
     !1      *yr(n-1-1,n-1) + yr(n-1,n-1+1) - yr(n-1-1,n-1))/den
   
      !a2 = rk_theta*(bond_angle(n-1-1) - theta_0)*((dot_prod(n-1-2,2)/
     !1      dot_prod(n-1-1,1))*yr(n-1-1,n-1) - yr(n-1-2,n-1-1))/den1

      !fba_y(n-1) = a1 + a2

      !a1 = -rk_theta*(bond_angle(n-1) - theta_0)*( (dot_prod(n-1-1,2)/
     !1  dot_prod(n-1,1))*zr(n-1,n-1+1) - (dot_prod(n-1-1,2)/dot_prod(n-1-1,1))
     !1      *zr(n-1-1,n-1) + zr(n-1,n-1+1) - zr(n-1-1,n-1))/den

      !a2 = rk_theta*(bond_angle(n-1-1) - theta_0)*((dot_prod(n-1-2,2)/
     !1      dot_prod(n-1-1,1))*zr(n-1-1,n-1) - zr(n-1-2,n-1-1))/den1

      !fba_z(n-1) = a1 + a2
!! }}}
!! particle n
!! {{{
!!     i = n
      !den = sinbond(n-1)

      !fba_x(n) = rk_theta*(bond_angle(n-1) - theta_0)*
     !1      ((dot_prod(n-2,2)/dot_prod(n-1,1))*xr(n-1,n) 
     !1      - xr(n-2,n-1))/den

      !fba_y(n) = rk_theta*(bond_angle(n-1) - theta_0)*
     !1      ((dot_prod(n-2,2)/dot_prod(n-1,1))*yr(n-1,n) 
     !1      - yr(n-2,n-1))/den

      !fba_z(n) = rk_theta*(bond_angle(n-1) - theta_0)*
     !1      ((dot_prod(n-2,2)/dot_prod(n-1,1))*zr(n-1,n) 
     !1      - zr(n-2,n-1))/den

!! }}}
!! Torsional angle forces
!! particles 1, 2, 3, n-2, n-1, and n are done outside of the loop
!! particle 1
!! {{{
!!     i = 1
      !coef =DFAC(1+1)

      !fta_x(1) = -coef*(-dot_prod(1+1,2)*xr(1+1,1+2) +
     !1       dot_prod(1+1,1)*xr(1+2,1+3) -
     !1      (1.0/x_prod(1))*(dot_prod(1+1,2)*dot_prod(1,2) -
     !1      dot_prod(1,3)*dot_prod(1+1,1))*(-dot_prod(1+1,1)*xr(1,1+1) +
     !1      dot_prod(1,2)*xr(1+1,1+2))) 

      !fta_y(1) = -coef*(-dot_prod(1+1,2)*yr(1+1,1+2) +
     !1      dot_prod(1+1,1)*yr(1+2,1+3) -
     !1      (1.0/x_prod(1))*(dot_prod(1+1,2)*dot_prod(1,2) -
     !1      dot_prod(1,3)*dot_prod(1+1,1))*(-dot_prod(1+1,1)*yr(1,1+1) +
     !1      dot_prod(1,2)*yr(1+1,1+2))) 

      !fta_z(1) = -coef*(-dot_prod(1+1,2)*zr(1+1,1+2) +
     !1      dot_prod(1+1,1)*zr(1+2,1+3) -
     !1      (1.0/x_prod(1))*(dot_prod(1+1,2)*dot_prod(1,2) -
     !1      dot_prod(1,3)*dot_prod(1+1,1))*(-dot_prod(1+1,1)*zr(1,1+1) +
     !1      dot_prod(1,2)*zr(1+1,1+2))) 

!! }}}
!! particle 2
!! {{{
!!     i = 2
      !coef =DFAC(2+1)

      !coef1 = DFAC(2)

      !a1 =  -coef*(-dot_prod(2+1,2)*xr(2+1,2+2) +
     !1      dot_prod(2+1,1)*xr(2+2,2+3) -
     !1      (1.0/x_prod(2))*(dot_prod(2+1,2)*dot_prod(2,2) -
     !1      dot_prod(2,3)*dot_prod(2+1,1))*(-dot_prod(2+1,1)*xr(2,2+1) +
     !1      dot_prod(2,2)*xr(2+1,2+2))) 

      !a2 = -coef1*(-dot_prod(2-1,2)*xr(2+1,2+2) +
     !1      dot_prod(2,2)*xr(2,2+1) - dot_prod(2,2)*xr(2-1,2) -
     !1      dot_prod(2,1)*xr(2+1,2+2) + 2.0*dot_prod(2-1,3)*xr(2,2+1) -
     !1      (1.0/x_prod(2-1))*(dot_prod(2,2)*dot_prod(2-1,2) -
     !1      dot_prod(2-1,3)*dot_prod(2,1))*(dot_prod(2,1)*xr(2-1,2) -
     !1      dot_prod(2-1,1)*xr(2,2+1) - dot_prod(2-1,2)*xr(2,2+1) +
     !1      dot_prod(2-1,2)*xr(2-1,2)) -
     !1      (1.0/x_prod(2))*(dot_prod(2,2)*dot_prod(2-1,2) -
     !1      dot_prod(2-1,3)*dot_prod(2,1))*(-dot_prod(2+1,1)*xr(2,2+1) +
     !1      dot_prod(2,2)*xr(2+1,2+2))) 

      !fta_x(2) = a1 + a2 

      !a1 = -dot_prod(3,2)*yr(3,4) + dot_prod(3,1)*yr(4,5) 
      !a1=a1-(dot_prod(3,2)*dot_prod(2,2) - dot_prod(2,3)*dot_prod(3,1))*
     !&     (-dot_prod(3,1)*yr(2,3) + dot_prod(2,2)*yr(3,4))/x_prod(2)
      !a1=-coef*a1

      !a2 = -coef1*(-dot_prod(2-1,2)*yr(2+1,2+2) +
     !1      dot_prod(2,2)*yr(2,2+1) - dot_prod(2,2)*yr(2-1,2) -
     !1      dot_prod(2,1)*yr(2+1,2+2) + 2.0*dot_prod(2-1,3)*yr(2,2+1) -
     !1      (1.0/x_prod(2-1))*(dot_prod(2,2)*dot_prod(2-1,2) -
     !1      dot_prod(2-1,3)*dot_prod(2,1))*(dot_prod(2,1)*yr(2-1,2) -
     !1      dot_prod(2-1,1)*yr(2,2+1) - dot_prod(2-1,2)*yr(2,2+1) +
     !1      dot_prod(2-1,2)*yr(2-1,2)) -
     !1      (1.0/x_prod(2))*(dot_prod(2,2)*dot_prod(2-1,2) -
     !1      dot_prod(2-1,3)*dot_prod(2,1))*(-dot_prod(2+1,1)*yr(2,2+1) +
     !1      dot_prod(2,2)*yr(2+1,2+2)))

      !fta_y(2) = a1 + a2 
!!     WRITE(MYUNIT,'(A,I6,5G15.5)') 'y i,a1,a2,coef,coef1,fta_y=',2,a1,a2,coef,coef1,fta_y(2)
      
      !a1 = -coef*(-dot_prod(3,2)*zr(3,4) +
     !1      dot_prod(3,1)*zr(4,5) -
     !1      (1.0/x_prod(2))*(dot_prod(3,2)*dot_prod(2,2) -
     !1      dot_prod(2,3)*dot_prod(3,1))*(-dot_prod(3,1)*zr(2,3) +
     !1      dot_prod(2,2)*zr(3,4))) 

      !a2 = -coef1*(-dot_prod(2-1,2)*zr(2+1,2+2) +
     !1      dot_prod(2,2)*zr(2,2+1) - dot_prod(2,2)*zr(2-1,2) -
     !1      dot_prod(2,1)*zr(2+1,2+2) + 2.0*dot_prod(2-1,3)*zr(2,2+1) -
     !1      (1.0/x_prod(2-1))*(dot_prod(2,2)*dot_prod(2-1,2) -
     !1      dot_prod(2-1,3)*dot_prod(2,1))*(dot_prod(2,1)*zr(2-1,2) -
     !1      dot_prod(2-1,1)*zr(2,2+1) - dot_prod(2-1,2)*zr(2,2+1) +
     !1      dot_prod(2-1,2)*zr(2-1,2)) -
     !1      (1.0/x_prod(2))*(dot_prod(2,2)*dot_prod(2-1,2) -
     !1      dot_prod(2-1,3)*dot_prod(2,1))*(-dot_prod(2+1,1)*zr(2,2+1) +
     !1      dot_prod(2,2)*zr(2+1,2+2))) 

      !fta_z(2) = a1 + a2 
!!     WRITE(MYUNIT,'(A,I6,4G15.5)') 'z i,a1,a2,coef,coef1=',2,a1,a2,coef,coef1
!! }}}
!! particle 3
!! {{{
!!     i = 3
      !coef=DFAC(3+1)

      !coef1=DFAC(3)

      !coef2=DFAC(3-1)

      !a1 = -coef*(-dot_prod(3+1,2)*xr(3+1,3+2) +
     !1      dot_prod(3+1,1)*xr(3+2,3+3) -
     !1      (1.0/x_prod(3))*(dot_prod(3+1,2)*dot_prod(3,2) -
     !1      dot_prod(3,3)*dot_prod(3+1,1))*(-dot_prod(3+1,1)*xr(3,3+1) +
     !1      dot_prod(3,2)*xr(3+1,3+2))) 

      !a2 = -coef1*(-dot_prod(3-1,2)*xr(3+1,3+2) +
     !1      dot_prod(3,2)*xr(3,3+1) - dot_prod(3,2)*xr(3-1,3) -
     !1      dot_prod(3,1)*xr(3+1,3+2) + 2.0*dot_prod(3-1,3)*xr(3,3+1) -
     !1      (1.0/x_prod(3-1))*(dot_prod(3,2)*dot_prod(3-1,2) -
     !1      dot_prod(3-1,3)*dot_prod(3,1))*(dot_prod(3,1)*xr(3-1,3) -
     !1      dot_prod(3-1,1)*xr(3,3+1) - dot_prod(3-1,2)*xr(3,3+1) +
     !1      dot_prod(3-1,2)*xr(3-1,3)) -
     !1      (1.0/x_prod(3))*(dot_prod(3,2)*dot_prod(3-1,2) -
     !1      dot_prod(3-1,3)*dot_prod(3,1))*(-dot_prod(3+1,1)*xr(3,3+1) +
     !1      dot_prod(3,2)*xr(3+1,3+2))) 

      !a3 = -coef2*(dot_prod(3-2,2)*xr(3,3+1) -
     !1      dot_prod(3-2,2)*xr(3-1,3) + dot_prod(3-1,2)*xr(3-2,3-1) +
     !1      dot_prod(3-1,1)*xr(3-2,3-1) - 2.0*dot_prod(3-2,3)*xr(3-1,3) -
     !1      (1.0/x_prod(3-2))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     !1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3-2,1)*xr(3-1,3) -
     !1      dot_prod(3-2,2)*xr(3-2,3-1)) -
     !1      (1.0/x_prod(3-1))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     !1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3,1)*xr(3-1,3) -
     !1      dot_prod(3-1,1)*xr(3,3+1) - dot_prod(3-1,2)*xr(3,3+1) +
     !1      dot_prod(3-1,2)*xr(3-1,3))) 

      !fta_x(3) = a1 + a2 + a3 
 
      !a1 = -coef*(-dot_prod(3+1,2)*yr(3+1,3+2) +
     !1      dot_prod(3+1,1)*yr(3+2,3+3) -
     !1      (1.0/x_prod(3))*(dot_prod(3+1,2)*dot_prod(3,2) -
     !1      dot_prod(3,3)*dot_prod(3+1,1))*(-dot_prod(3+1,1)*yr(3,3+1) +
     !1      dot_prod(3,2)*yr(3+1,3+2))) 
      
      !a2 = -coef1*(-dot_prod(3-1,2)*yr(3+1,3+2) +
     !1      dot_prod(3,2)*yr(3,3+1) - dot_prod(3,2)*yr(3-1,3) -
     !1      dot_prod(3,1)*yr(3+1,3+2) + 2.0*dot_prod(3-1,3)*yr(3,3+1) -
     !1      (1.0/x_prod(3-1))*(dot_prod(3,2)*dot_prod(3-1,2) -
     !1      dot_prod(3-1,3)*dot_prod(3,1))*(dot_prod(3,1)*yr(3-1,3) -
     !1      dot_prod(3-1,1)*yr(3,3+1) - dot_prod(3-1,2)*yr(3,3+1) +
     !1      dot_prod(3-1,2)*yr(3-1,3)) -
     !1      (1.0/x_prod(3))*(dot_prod(3,2)*dot_prod(3-1,2) -
     !1      dot_prod(3-1,3)*dot_prod(3,1))*(-dot_prod(3+1,1)*yr(3,3+1) +
     !1      dot_prod(3,2)*yr(3+1,3+2))) 

      !a3 = -coef2*(dot_prod(3-2,2)*yr(3,3+1) -
     !1      dot_prod(3-2,2)*yr(3-1,3) + dot_prod(3-1,2)*yr(3-2,3-1) +
     !1      dot_prod(3-1,1)*yr(3-2,3-1) - 2.0*dot_prod(3-2,3)*yr(3-1,3) -
     !1      (1.0/x_prod(3-2))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     !1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3-2,1)*yr(3-1,3) -
     !1      dot_prod(3-2,2)*yr(3-2,3-1)) -
     !1      (1.0/x_prod(3-1))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     !1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3,1)*yr(3-1,3) -
     !1      dot_prod(3-1,1)*yr(3,3+1) - dot_prod(3-1,2)*yr(3,3+1) +
     !1      dot_prod(3-1,2)*yr(3-1,3)))

      !fta_y(3) = a1 + a2 + a3 
 
      !a1 = -coef*(-dot_prod(3+1,2)*zr(3+1,3+2) +
     !1      dot_prod(3+1,1)*zr(3+2,3+3) -
     !1      (1.0/x_prod(3))*(dot_prod(3+1,2)*dot_prod(3,2) -
     !1      dot_prod(3,3)*dot_prod(3+1,1))*(-dot_prod(3+1,1)*zr(3,3+1) +
     !1      dot_prod(3,2)*zr(3+1,3+2))) 

      !a2 =  -coef1*(-dot_prod(3-1,2)*zr(3+1,3+2) +
     !1      dot_prod(3,2)*zr(3,3+1) - dot_prod(3,2)*zr(3-1,3) -
     !1      dot_prod(3,1)*zr(3+1,3+2) + 2.0*dot_prod(3-1,3)*zr(3,3+1) -
     !1      (1.0/x_prod(3-1))*(dot_prod(3,2)*dot_prod(3-1,2) -
     !1      dot_prod(3-1,3)*dot_prod(3,1))*(dot_prod(3,1)*zr(3-1,3) -
     !1      dot_prod(3-1,1)*zr(3,3+1) - dot_prod(3-1,2)*zr(3,3+1) +
     !1      dot_prod(3-1,2)*zr(3-1,3)) -
     !1      (1.0/x_prod(3))*(dot_prod(3,2)*dot_prod(3-1,2) -
     !1      dot_prod(3-1,3)*dot_prod(3,1))*(-dot_prod(3+1,1)*zr(3,3+1) +
     !1      dot_prod(3,2)*zr(3+1,3+2))) 

      !a3 = -coef2*(dot_prod(3-2,2)*zr(3,3+1) -
     !1      dot_prod(3-2,2)*zr(3-1,3) + dot_prod(3-1,2)*zr(3-2,3-1) +
     !1      dot_prod(3-1,1)*zr(3-2,3-1) - 2.0*dot_prod(3-2,3)*zr(3-1,3) -
     !1      (1.0/x_prod(3-2))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     !1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3-2,1)*zr(3-1,3) -
     !1      dot_prod(3-2,2)*zr(3-2,3-1)) -
     !1      (1.0/x_prod(3-1))*(dot_prod(3-1,2)*dot_prod(3-2,2) -
     !1      dot_prod(3-2,3)*dot_prod(3-1,1))*(dot_prod(3,1)*zr(3-1,3) -
     !1      dot_prod(3-1,1)*zr(3,3+1) - dot_prod(3-1,2)*zr(3,3+1) +
     !1      dot_prod(3-1,2)*zr(3-1,3))) 

      !fta_z(3) = a1 + a2 + a3 
!! }}}
!! particles 4 to n-3
!! {{{
      !do i = 4, n-3

         !coef=DFAC(i+1)

         !coef1=DFAC(i)

         !coef2=DFAC(i-1)

         !coef3=DFAC(i-2)

         !a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     !1      dot_prod(i+1,1)*xr(i+2,i+3) -
     !1      (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     !1      dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     !1      dot_prod(i,2)*xr(i+1,i+2))) 

         !a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     !1      dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     !1      dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     !1      (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     !1      dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     !1      dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     !1      dot_prod(i-1,2)*xr(i-1,i)) -
     !1      (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     !1      dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     !1      dot_prod(i,2)*xr(i+1,i+2))) 

         !a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     !1      dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     !1      dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     !1      (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     !1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     !1      dot_prod(i-2,2)*xr(i-2,i-1)) -
     !1      (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     !1  dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     !1      dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     !1      dot_prod(i-1,2)*xr(i-1,i))) 

         !a4 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) -
     !1      dot_prod(i-2,1)*xr(i-3,i-2) -
     !1      (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     !1      dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     !1      dot_prod(i-2,2)*xr(i-2,i-1))) 

         !fta_x(i) = a1 + a2 + a3 + a4 

         !a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     !1      dot_prod(i+1,1)*yr(i+2,i+3) -
     !1      (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     !1      dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     !1      dot_prod(i,2)*yr(i+1,i+2))) 

         !a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     !1      dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     !1      dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     !1      (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     !1      dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     !1      dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     !1      dot_prod(i-1,2)*yr(i-1,i)) -
     !1      (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     !1      dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     !1      dot_prod(i,2)*yr(i+1,i+2))) 

         !a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     !1      dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     !1      dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     !1      (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     !1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     !1      dot_prod(i-2,2)*yr(i-2,i-1)) -
     !1      (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     !1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     !1      dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     !1      dot_prod(i-1,2)*yr(i-1,i))) 

         !a4 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) -
     !1      dot_prod(i-2,1)*yr(i-3,i-2) -
     !1      (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     !1      dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     !1      dot_prod(i-2,2)*yr(i-2,i-1))) 

         !fta_y(i) = a1 + a2 + a3 + a4 

         !a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     !1      dot_prod(i+1,1)*zr(i+2,i+3) -
     !1      (1.0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     !1      dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     !1      dot_prod(i,2)*zr(i+1,i+2))) 

         !a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     !1      dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     !1      dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     !1      (1.0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     !1      dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     !1      dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     !1      dot_prod(i-1,2)*zr(i-1,i)) -
     !1      (1.0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     !1      dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     !1      dot_prod(i,2)*zr(i+1,i+2))) 

         !a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     !1      dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     !1      dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     !1      (1.0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     !1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     !1      dot_prod(i-2,2)*zr(i-2,i-1)) -
     !1      (1.0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     !1      dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     !1      dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     !1      dot_prod(i-1,2)*zr(i-1,i))) 
      
         !a4 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) -
     !1      dot_prod(i-2,1)*zr(i-3,i-2) -
     !1      (1.0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     !1      dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     !1      dot_prod(i-2,2)*zr(i-2,i-1))) 

         !fta_z(i) = a1 + a2 + a3 + a4 

      !enddo
!! }}}
!! particle n-2
!! {{{
!!     i = n-2
      !coef1=DFAC(n-2)

      !coef2=DFAC(n-2-1)

      !coef3=DFAC(n-2-2)

      !a1 = -coef1*(-dot_prod(n-2-1,2)*xr(n-2+1,n-2+2) + 
     !1      dot_prod(n-2,2)*xr(n-2,n-2+1) - dot_prod(n-2,2)*xr(n-2-1,n-2) -
     !1      dot_prod(n-2,1)*xr(n-2+1,n-2+2) + 2.0*dot_prod(n-2-1,3)*xr(n-2,n-2+1) -
     !1      (1.0/x_prod(n-2-1))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     !1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*xr(n-2-1,n-2) -
     !1      dot_prod(n-2-1,1)*xr(n-2,n-2+1) - dot_prod(n-2-1,2)*xr(n-2,n-2+1) +
     !1      dot_prod(n-2-1,2)*xr(n-2-1,n-2)) -
     !1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     !1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(-dot_prod(n-2+1,1)*xr(n-2,n-2+1) +
     !1      dot_prod(n-2,2)*xr(n-2+1,n-2+2)))


      !a2 = -coef2*(dot_prod(n-2-2,2)*xr(n-2,n-2+1) -
     !1      dot_prod(n-2-2,2)*xr(n-2-1,n-2) + dot_prod(n-2-1,2)*xr(n-2-2,n-2-1) +
     !1      dot_prod(n-2-1,1)*xr(n-2-2,n-2-1) - 2.0*dot_prod(n-2-2,3)*xr(n-2-1,n-2) -
     !1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     !1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2-2,1)*xr(n-2-1,n-2) -
     !1      dot_prod(n-2-2,2)*xr(n-2-2,n-2-1)) -
     !1      (1.0/x_prod(n-2-1))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     !1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2,1)*xr(n-2-1,n-2) -
     !1      dot_prod(n-2-1,1)*xr(n-2,n-2+1) - dot_prod(n-2-1,2)*xr(n-2,n-2+1) +
     !1      dot_prod(n-2-1,2)*xr(n-2-1,n-2))) 
      
      !a3 = -coef3*(dot_prod(n-2-3,2)*xr(n-2-2,n-2-1) -
     !1      dot_prod(n-2-2,1)*xr(n-2-3,n-2-2) -
     !1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-2,2)*dot_prod(n-2-3,2) -
     !1      dot_prod(n-2-3,3)*dot_prod(n-2-2,1))*(dot_prod(n-2-2,1)*xr(n-2-1,n-2) -
     !1      dot_prod(n-2-2,2)*xr(n-2-2,n-2-1))) 

      !fta_x(n-2) = a1 + a2 + a3 

      !a1 =  -coef1*(-dot_prod(n-2-1,2)*yr(n-2+1,n-2+2) +  
     !1      dot_prod(n-2,2)*yr(n-2,n-2+1) - dot_prod(n-2,2)*yr(n-2-1,n-2) -
     !1      dot_prod(n-2,1)*yr(n-2+1,n-2+2) + 2.0*dot_prod(n-2-1,3)*yr(n-2,n-2+1) -
     !1      (1.0/x_prod(n-2-1))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     !1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*yr(n-2-1,n-2) -
     !1      dot_prod(n-2-1,1)*yr(n-2,n-2+1) - dot_prod(n-2-1,2)*yr(n-2,n-2+1) +
     !1      dot_prod(n-2-1,2)*yr(n-2-1,n-2)) -
     !1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     !1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(-dot_prod(n-2+1,1)*yr(n-2,n-2+1) +
     !1      dot_prod(n-2,2)*yr(n-2+1,n-2+2))) 

      !a2 =  -coef2*(dot_prod(n-2-2,2)*yr(n-2,n-2+1) -
     !1      dot_prod(n-2-2,2)*yr(n-2-1,n-2) + dot_prod(n-2-1,2)*yr(n-2-2,n-2-1) +
     !1      dot_prod(n-2-1,1)*yr(n-2-2,n-2-1) - 2.0*dot_prod(n-2-2,3)*yr(n-2-1,n-2) -
     !1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     !1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2-2,1)*yr(n-2-1,n-2) -
     !1      dot_prod(n-2-2,2)*yr(n-2-2,n-2-1)) -
     !1      (1.0/x_prod(n-2-1))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     !1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2,1)*yr(n-2-1,n-2) -
     !1      dot_prod(n-2-1,1)*yr(n-2,n-2+1) - dot_prod(n-2-1,2)*yr(n-2,n-2+1) +
     !1      dot_prod(n-2-1,2)*yr(n-2-1,n-2)))

      !a3 = -coef3*(dot_prod(n-2-3,2)*yr(n-2-2,n-2-1) -
     !1      dot_prod(n-2-2,1)*yr(n-2-3,n-2-2) -
     !1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-2,2)*dot_prod(n-2-3,2) -
     !1      dot_prod(n-2-3,3)*dot_prod(n-2-2,1))*(dot_prod(n-2-2,1)*yr(n-2-1,n-2) -
     !1      dot_prod(n-2-2,2)*yr(n-2-2,n-2-1))) 

      !fta_y(n-2) = a1 + a2 + a3 
 
      !a1 = -coef1*(-dot_prod(n-2-1,2)*zr(n-2+1,n-2+2) +  
     !1      dot_prod(n-2,2)*zr(n-2,n-2+1) - dot_prod(n-2,2)*zr(n-2-1,n-2) -
     !1      dot_prod(n-2,1)*zr(n-2+1,n-2+2) + 2.0*dot_prod(n-2-1,3)*zr(n-2,n-2+1) -
     !1      (1.0/x_prod(n-2-1))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     !1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*zr(n-2-1,n-2) -
     !1      dot_prod(n-2-1,1)*zr(n-2,n-2+1) - dot_prod(n-2-1,2)*zr(n-2,n-2+1) +
     !1      dot_prod(n-2-1,2)*zr(n-2-1,n-2)) -
     !1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-2-1,2) -
     !1      dot_prod(n-2-1,3)*dot_prod(n-2,1))*(-dot_prod(n-2+1,1)*zr(n-2,n-2+1) +
     !1      dot_prod(n-2,2)*zr(n-2+1,n-2+2))) 

      !a2 = -coef2*(dot_prod(n-2-2,2)*zr(n-2,n-2+1) -
     !1      dot_prod(n-2-2,2)*zr(n-2-1,n-2) + dot_prod(n-2-1,2)*zr(n-2-2,n-2-1) +
     !1      dot_prod(n-2-1,1)*zr(n-2-2,n-2-1) - 2.0*dot_prod(n-2-2,3)*zr(n-2-1,n-2) -
     !1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     !1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2-2,1)*zr(n-2-1,n-2) -
     !1      dot_prod(n-2-2,2)*zr(n-2-2,n-2-1)) -
     !1      (1.0/x_prod(n-2-1))*(dot_prod(n-2-1,2)*dot_prod(n-2-2,2) -
     !1      dot_prod(n-2-2,3)*dot_prod(n-2-1,1))*(dot_prod(n-2,1)*zr(n-2-1,n-2) -
     !1      dot_prod(n-2-1,1)*zr(n-2,n-2+1) - dot_prod(n-2-1,2)*zr(n-2,n-2+1) +
     !1      dot_prod(n-2-1,2)*zr(n-2-1,n-2))) 

      !a3 = -coef3*(dot_prod(n-2-3,2)*zr(n-2-2,n-2-1) -
     !1      dot_prod(n-2-2,1)*zr(n-2-3,n-2-2) -
     !1      (1.0/x_prod(n-2-2))*(dot_prod(n-2-2,2)*dot_prod(n-2-3,2) -
     !1      dot_prod(n-2-3,3)*dot_prod(n-2-2,1))*(dot_prod(n-2-2,1)*zr(n-2-1,n-2) -
     !1      dot_prod(n-2-2,2)*zr(n-2-2,n-2-1))) 

      !fta_z(n-2) = a1 + a2 + a3 
!! }}}
!! particle n-1
!! {{{
!!     i = n-1
      !coef2=DFAC(n-1-1)

      !coef3=DFAC(n-1-2)

      !a1 = -coef2*(dot_prod(n-1-2,2)*xr(n-1,n-1+1) - 
     !1      dot_prod(n-1-2,2)*xr(n-1-1,n-1) +
     !1      dot_prod(n-1-1,2)*xr(n-1-2,n-1-1) +  dot_prod(n-1-1,1)*xr(n-1-2,n-1-1) -
     !1      2.0*dot_prod(n-1-2,3)*xr(n-1-1,n-1) -
     !1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     !1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1-2,1)*xr(n-1-1,n-1) -
     !1      dot_prod(n-1-2,2)*xr(n-1-2,n-1-1)) -
     !1      (1.0/x_prod(n-1-1))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     !1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1,1)*xr(n-1-1,n-1) -
     !1      dot_prod(n-1-1,1)*xr(n-1,n-1+1) - dot_prod(n-1-1,2)*xr(n-1,n-1+1) +
     !1      dot_prod(n-1-1,2)*xr(n-1-1,n-1))) 

      !a2 = -coef3*(dot_prod(n-1-3,2)*xr(n-1-2,n-1-1) - 
     !1      dot_prod(n-1-2,1)*xr(n-1-3,n-1-2) -
     !1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-2,2)*dot_prod(n-1-3,2) -
     !1      dot_prod(n-1-3,3)*dot_prod(n-1-2,1))*(dot_prod(n-1-2,1)*xr(n-1-1,n-1) -
     !1      dot_prod(n-1-2,2)*xr(n-1-2,n-1-1))) 

      !fta_x(n-1) = a1 + a2  

      !a1 = -coef2*(dot_prod(n-1-2,2)*yr(n-1,n-1+1) - 
     !1      dot_prod(n-1-2,2)*yr(n-1-1,n-1) +
     !1      dot_prod(n-1-1,2)*yr(n-1-2,n-1-1) +  dot_prod(n-1-1,1)*yr(n-1-2,n-1-1) -
     !1      2.0*dot_prod(n-1-2,3)*yr(n-1-1,n-1) -
     !1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     !1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1-2,1)*yr(n-1-1,n-1) -
     !1      dot_prod(n-1-2,2)*yr(n-1-2,n-1-1)) -
     !1      (1.0/x_prod(n-1-1))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     !1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1,1)*yr(n-1-1,n-1) -
     !1  dot_prod(n-1-1,1)*yr(n-1,n-1+1) - dot_prod(n-1-1,2)*yr(n-1,n-1+1) +
     !1      dot_prod(n-1-1,2)*yr(n-1-1,n-1))) 

      !a2 = -coef3*(dot_prod(n-1-3,2)*yr(n-1-2,n-1-1) - 
     !1      dot_prod(n-1-2,1)*yr(n-1-3,n-1-2) -
     !1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-2,2)*dot_prod(n-1-3,2) -
     !1      dot_prod(n-1-3,3)*dot_prod(n-1-2,1))*(dot_prod(n-1-2,1)*yr(n-1-1,n-1) -
     !1      dot_prod(n-1-2,2)*yr(n-1-2,n-1-1))) 

      !fta_y(n-1) = a1 + a2  

      !a1 = -coef2*(dot_prod(n-1-2,2)*zr(n-1,n-1+1) - 
     !1      dot_prod(n-1-2,2)*zr(n-1-1,n-1) +
     !1      dot_prod(n-1-1,2)*zr(n-1-2,n-1-1) +  dot_prod(n-1-1,1)*zr(n-1-2,n-1-1) -
     !1      2.0*dot_prod(n-1-2,3)*zr(n-1-1,n-1) -
     !1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     !1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1-2,1)*zr(n-1-1,n-1) -
     !1      dot_prod(n-1-2,2)*zr(n-1-2,n-1-1)) -
     !1      (1.0/x_prod(n-1-1))*(dot_prod(n-1-1,2)*dot_prod(n-1-2,2) -
     !1      dot_prod(n-1-2,3)*dot_prod(n-1-1,1))*(dot_prod(n-1,1)*zr(n-1-1,n-1) -
     !1      dot_prod(n-1-1,1)*zr(n-1,n-1+1) - dot_prod(n-1-1,2)*zr(n-1,n-1+1) +
     !1      dot_prod(n-1-1,2)*zr(n-1-1,n-1))) 

      !a2 = -coef3*(dot_prod(n-1-3,2)*zr(n-1-2,n-1-1) - 
     !1      dot_prod(n-1-2,1)*zr(n-1-3,n-1-2) -
     !1      (1.0/x_prod(n-1-2))*(dot_prod(n-1-2,2)*dot_prod(n-1-3,2) -
     !1      dot_prod(n-1-3,3)*dot_prod(n-1-2,1))*(dot_prod(n-1-2,1)*zr(n-1-1,n-1) -
     !1      dot_prod(n-1-2,2)*zr(n-1-2,n-1-1))) 

      !fta_z(n-1) = a1 + a2 
!! }}} 
!! particle n
!! {{{
!!     i = n
      !coef3=DFAC(n-2)

      !fta_x(n) = -coef3*(dot_prod(n-3,2)*xr(n-2,n-1) 
     !1      - dot_prod(n-2,1)*xr(n-3,n-2) -
     !1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-3,2) -
     !1      dot_prod(n-3,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*xr(n-1,n) -
     !1      dot_prod(n-2,2)*xr(n-2,n-1))) 

      !fta_y(n) = -coef3*(dot_prod(n-3,2)*yr(n-2,n-1) - 
     !1      dot_prod(n-2,1)*yr(n-3,n-2) 
     !1      - (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-3,2) 
     !1      - dot_prod(n-3,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*yr(n-1,n) -
     !1      dot_prod(n-2,2)*yr(n-2,n-1))) 

      !fta_z(n) = -coef3*(dot_prod(n-3,2)*zr(n-2,n-1) - 
     !1      dot_prod(n-2,1)*zr(n-3,n-2) -
     !1      (1.0/x_prod(n-2))*(dot_prod(n-2,2)*dot_prod(n-3,2) -
     !1      dot_prod(n-3,3)*dot_prod(n-2,1))*(dot_prod(n-2,1)*zr(n-1,n) -
     !1      dot_prod(n-2,2)*zr(n-2,n-1))) 
!! }}}
!! Total up the gradients
!! {{{
      !do i = 1, n
!!        IF (I.EQ.2) THEN
!!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'i,fnbx,fbx,fbax,ftax=',i,fnb_x(i),fb_x(i),fba_x(i),fta_x(i)
!!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'i,fnby,fby,fbay,ftay=',i,fnb_y(i),fb_y(i),fba_y(i),fta_y(i)
!!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'i,fnbz,fbz,fbaz,ftaz=',i,fnb_z(i),fb_z(i),fba_z(i),fta_z(i)
!!        ENDIF
         !fx(i) = fnb_x(i) + fb_x(i) + fba_x(i) + fta_x(i) 
         !fy(i) = fnb_y(i) + fb_y(i) + fba_y(i) + fta_y(i) 
         !fz(i) = fnb_z(i) + fb_z(i) + fba_z(i) + fta_z(i) 
      !enddo

      !do i = 1, n
         !j = (i-1)*3
         !fq(j+1) = -fx(i)
         !fq(j+2) = -fy(i)
         !fq(j+3) = -fz(i)
      !enddo
!! }}}
      !return
      !end
!! }}}

!!> @brief  Fill the parameter arrays

      !SUBROUTINE PARAM_ARRAY_BLN(LJREP,LJATT,A,BEADLETTER,BLNSSTRUCT,
     !&   LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, 
     !&   HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, N)
!! {{{
!! DECLARATIONS {{{

!include  "blnvars.inc.f90"

      !INTEGER J1, ICOUNT

      !DOUBLE PRECISION, DIMENSION(N,N) :: LJREP,LJATT
      !DOUBLE PRECISION, DIMENSION(N,4) :: A

      !DOUBLE PRECISION LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN
      !DOUBLE PRECISION HABLN, HBBLN, HCBLN, HDBLN
      !DOUBLE PRECISION EABLN, EBBLN, ECBLN, EDBLN
      !DOUBLE PRECISION TABLN, TBBLN, TCBLN, TDBLN

      !CHARACTER(LEN=1) BEADLETTER(N), BLNSSTRUCT(N)
!! }}}
!! AMINO ACID TYPES {{{
!! B=1 L=2 N=3
!!
      !DO J1=1,N
         !IF (BEADLETTER(J1).EQ.'B') THEN
            !NTYPE(J1)=1
         !ELSEIF (BEADLETTER(J1).EQ.'L') THEN
            !NTYPE(J1)=2
         !ELSEIF (BEADLETTER(J1).EQ.'N') THEN
            !NTYPE(J1)=3
         !ELSE
            !PRINT '(A,A1)','ERROR IN PARAM_ARRAY_BLN, UNRECOGNISED BEAD TYPE: ',BEADLETTER(J1)
            !STOP
         !ENDIF
      !ENDDO

!! }}}

!! PARAMETERS FOR THE DIHEDRAL ANGLE POTENTIAL 
!! THE END BONDS HAVE NO DIHEDRAL TERM, SO THE TOTAL NUMBER OF TERMS
!! IS N-3 (E.G. 4 ATOMS, 1 DIHEDRAL). NON-ZERO TERMS FOR
!! 2 TO 3, 3 TO 4, 4 TO 5, ... , N-3 TO N-2, N-2 TO N-1. THE
!! H, E AND T PARAMETERS ARE DEFINED FOR THE FIRST BEAD OF EACH EDGE,
!! I.E. FOR 2, 3, 4, ..., N-2.
!! {{{

      !A_BLN(1:N)=0.0D0
      !B_BLN(1:N)=0.0D0
      !C_BLN(1:N)=0.0D0
      !D_BLN(1:N)=0.0D0
      !DO I=1,N-3
         !IF (BLNSSTRUCT(I).EQ.'H') THEN
            !A_BLN(I+1)=HABLN
            !B_BLN(I+1)=HBBLN
            !C_BLN(I+1)=HCBLN
            !D_BLN(I+1)=HDBLN
         !ELSE IF (BLNSSTRUCT(I).EQ.'E') THEN
            !A_BLN(I+1)=EABLN
            !B_BLN(I+1)=EBBLN
            !C_BLN(I+1)=ECBLN
            !D_BLN(I+1)=EDBLN
         !ELSE IF (BLNSSTRUCT(I).EQ.'T') THEN
            !A_BLN(I+1)=TABLN
            !B_BLN(I+1)=TBBLN
            !C_BLN(I+1)=TCBLN
            !D_BLN(I+1)=TDBLN
         !ELSE
            !PRINT '(A,A1)','ERROR IN PARAM_ARRAYBLN, UNRECOGNISED SS TYPE: ',BLNSSTRUCT(J1)
            !STOP
         !ENDIF
!!        PRINT '(A,I6,A,A1,A,4F12.4)','I+1=',I+1,' SYMBOL=',BLNSSTRUCT(I),' A,B,C,D=',A_BLN(I+1),B_BLN(I+1),C_BLN(I+1),D_BLN(I+1)
      !ENDDO
!! }}}
!!  PARAMETERS FOR THE L-J INTERACTION BETWEEN NON-BONDED PARTICLES {{{

      !DO I = 1, N-1
         !DO J = I+1, N
            !IF (NTYPE(I) .EQ. 3 .OR. NTYPE(J) .EQ. 3) THEN
               !LJREP(I,J) = LJREPNN
               !LJATT(I,J) = LJATTNN
               !LJREP(J,I) = LJREPNN
               !LJATT(J,I) = LJATTNN
            !ELSE IF (NTYPE(I) .EQ. 1 .AND. NTYPE(J) .EQ. 1) THEN
               !LJREP(I,J) = LJREPBB
               !LJATT(I,J) = LJATTBB
               !LJREP(J,I) = LJREPBB
               !LJATT(J,I) = LJATTBB
            !ELSE
               !LJREP(I,J) = LJREPLL
               !LJATT(I,J) = LJATTLL
               !LJREP(J,I) = LJREPLL
               !LJATT(J,I) = LJATTLL
            !ENDIF
         !ENDDO
      !ENDDO
!! }}}

      !RETURN
      !END
!! }}}
!! }}}
! end generic BLN
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
        SUBROUTINE G46MERDIFF(N,QO,GRAD,ENERGY,GTEST)
! {{{ 

include "blnvars.inc.f90"

        STEST=.FALSE.

        CALL GPARAM_ARRAY(N,AB,CD)
        CALL CALC_INT_COORDS(N,QO,PARAM,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
        CALL CALC_ENERGY(N,QO,ENERGY,PARAM,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
        CALL CALC_GRADIENT(N,QO,GRAD,PARAM,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)

        IF (.NOT.STEST) RETURN
        CALL CALC_DYN(N,QO,PARAM,R,DR,DOT_PROD,X_PROD, ANG,RADII,NTYPE)

        RETURN
        END
! }}}
!
! Doxygen: gparam_array {{{
!>
!> \brief Fill the parameter arrays which specify interaction potentials
!> \param N INTEGER  - number of particles
!> \param a_param \param b_param - LJ interaction between non-bonded particles 
!> \param c_param \param d_param - dihedral angle potential
!>
!}}}
        SUBROUTINE GPARAM_ARRAY(N,AB,CD)
! {{{
include "blnvars.inc.f90" 

! Specify amino acid types by filling in the array ntype(:) {{{
     
! }}}
! Go-like model connectivities: fill in array CONNECT(:,:) {{{
!
        DO J1=1,N
           DO J2=J1,N
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

! }}}
! Parameters for the dihedral angle potential: fill in arrays c_param(:,:), d_param(:,:) {{{

        DO I = 1, N-3
          ICOUNT = 0

          DO J = 0,3
            IF(NTYPE(I+J) .EQ. 3)THEN
              ICOUNT = ICOUNT + 1
            ENDIF
          ENDDO

          IF (ICOUNT .GE. 2) THEN
            CD(I+1,1:2) = (/ 0.0, 0.2*EPSILON /)
          ELSE
            CD(I+1,1:2) = 1.2*EPSILON
        ENDIF

        ICOUNT = 0

        ENDDO
! }}}
! Parameters for the L-J interaction between non-bonded particles:
! arrays a_param(:,:), b_param(:,:)
! {{{

        DO I = 1, N
           DO J = 1, N

           IF (NTYPE(I) .EQ. 3 .OR. NTYPE(J) .EQ. 3) THEN
             AB(I,J,1) = 1.0*EPSILON 
             AB(I,J,2) = 0.0 
           ELSEIF (NTYPE(I) .EQ. 1 .AND. NTYPE(J) .EQ. 1)THEN
             AB(I,J,1) =  EPSILON
             IF (CONNECT(I,J)) THEN
                AB(I,J,2) = -EPSILON 
             ELSE
                AB(I,J,2) = 0.0D0
             ENDIF
           ELSE
             AB(I,J,1:2) = EPSILON*2.0/3.0 
           ENDIF
   
           ENDDO
        ENDDO
! }}}
        RETURN
        END
! }}}
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

        SUBROUTINE P46MERDIFF(N,QO,GRAD,ENERGY,GTEST)
! {{{

include "blnvars.inc.f90"

        STEST=.FALSE.

        CALL PARAM_ARRAY(N,AB,CD)
        CALL CALC_INT_COORDS(N,QO,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
        CALL CALC_ENERGY(N,QO,ENERGY,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)

        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN

        CALL CALC_GRADIENT(N,QO,GRAD,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)

        IF (.NOT.STEST) RETURN

        CALL CALC_DYN(N,QO,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)

        RETURN
        END
! }}}
!> Calculate the internal coordinates

        SUBROUTINE CALC_INT_COORDS(N,QO,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
! {{{
include "blnvars.inc.f90" 

        DO I = 1, N
          J = (I-1)*3
          DO K=1,3
            R(I,K) = QO(J+K)
          ENDDO
        ENDDO

! INTER-PARTICLE DISTANCES

        DO I = 1, N-1
          DO J = I+1, N
            DR(I,J,1:3) = R(J,1:3) - R(I,1:3)
            RADII(I,J) = DSQRT(SUM(DR(I,J,1:3)**2))
            RADII(J,I) = RADII(I,J)
          ENDDO
        ENDDO

! DOT PRODUCTS BETWEEN BOND VECTORS

      DO I = 1, N-1
        IF ( I .LE. N-3 ) THEN JMAX=I+2
        IF ( I .EQ. N-2 ) THEN JMAX=I+1
        IF ( I .LE. N-1 ) THEN JMAX=I
        DO J=I,JMAX
         ! K=1..3   FOR I=1..N-3
         ! K=1..2   FOR I=N-2
         ! K=1      FOR I=N-1
         K=J-I+1 
         DOT_PROD(I,K) = SUM(DR(I,I+1,1:3)*DR(J,J+1,1:3))
        ENDDO
      ENDDO

! CROSS-PRODUCTS BETWEEN ADJACENT BOND VECTORS

        DO I = 1, N-2
           X_PROD(I) = DOT_PROD(I,1)*DOT_PROD(I+1,1)-DOT_PROD(I,2)*DOT_PROD(I,2)   
        ENDDO

! BOND ANGLES

        DO I = 1, N-2
            COS_THETA=-DOT_PROD(I,2)/(DSQRT(DOT_PROD(I,1)*DOT_PROD(I+1,1)))
            ANG(I+1,1) = DACOS(COS_THETA)
        ENDDO

! TORSIONAL ANGLES

        DO I = 1, N-3
            COS_PHI = (DOT_PROD(I,2)*DOT_PROD(I+1,2)-DOT_PROD(I,3)*DOT_PROD(I+1,1))
            COS_PHI = COS_PHI/DSQRT(X_PROD(I)*X_PROD(I+1))
            IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=COS_PHI/ABS(COS_PHI)
            ANG(I+1,2) = DACOS(COS_PHI)
        ENDDO

        RETURN
        END

! }}} 
!> Calculate the energy
        SUBROUTINE CALC_ENERGY(N,QO,ENERGY,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
! {{{
       include "blnvars.inc.f90"   

       S(6)=SIGMA**6; E=0.0D0 
        
! non-bonded
        DO I = 1, N-2
          DO J = I+2, N
            RAD(6)=RADII(I,J)**6
            E(1) = E(1) + 4.0*((AB(I,J,1)*S(12)/(RAD(12))) + (AB(I,J,2)*S(6)/RAD(6)))
          ENDDO
        ENDDO

! bonded
        DO I = 1, N-1
          E(2) = E(2) + 0.5*RK_R*(RADII(I,I+1)-SIGMA)*(RADII(I,I+1)-SIGMA)
        ENDDO

! bond angles 
        DO I = 2, N-1
          E(3) = E(3) + 0.5*RK_THETA*(ANG(I,1)-THETA_0)*(ANG(I,1)-THETA_0)
        ENDDO

! torsional angles 
        DO I = 2, N-2
          E(4) = E(4) + CD(I,1)*(1.0 + COS(ANG(I,2)))+CD(I,2)*(1.0+COS(3.0*ANG(I,2)))
        ENDDO

        ENERGY=SUM(E)

        RETURN
        END
! }}}
! Calculate the gradients

        SUBROUTINE CALC_GRADIENT(N,QO,FQ,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
! {{{
include "blnvars.inc.f90" 

      S(6)=S(1)**6 ; S(12)=S(6)**2

      FNB= 0.0D0; FB= 0.0D0 ; FBA= 0.0D0 ; FTA= 0.0D0 ; F= 0.0D0 

! ..... Non-bonded interaction forces ..... {{{

        DO I = 1, N-2
          DO J = I+2, N
  
            RAD(1)=RADII(I,J) ; RAD(7)=RAD(1)**7 ; RAD(14)=RAD(7)**2 ; RAD(8)=RAD(7)*RAD(1)
        
            DF = 2.0*AB(I,J,1)*S(12)/RAD(14) + AB(I,J,2)*S(6)/RAD(8)
            DF=-24.0*DF

            FRR(1:3) = DF*DR(I,J,1:3) 
            FNB(I,1:3) = FRR(1:3) + FNB(I,1:3)
            FNB(J,1:3) = -FRR(1:3) + FNB(J,1:3)

          ENDDO
        ENDDO
! }}}
! ... Bond interaction forces ... {{{

        DO I = 1, N-1
           RVAR = SIGMA/RADII(I,I+1) 

           DF = RK_R*(1.0 - RVAR) 
           FRR(1:3) = DF*DR(I,I+1,1:3) 
           FB(I,1:3) = FRR(1:3) + FB(I,1:3)
           FB(I+1,1:3) = -FRR(1:3) + FB(I+1,1:3)

        ENDDO
! }}}
! bond angle forces: particle 1
! particles 1, 2, N-1, and N are done outside of the loop
! {{{
        I = 1
        DEN = DSIN(ANG(I+1,1))*DSQRT(DOT_PROD(I+1,1)*DOT_PROD(I,1))
        RNUM = RK_THETA*(ANG(I+1,1) - THETA_0)
        FBA(I,1:3)=-RNUM*((DOT_PROD(I,2)/DOT_PROD(I,1))*DR(I,I+1,1:3)-DR(I+1,I+2,1:3))/DEN

     ! particle 2

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

! particles 3 thru n-2 

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

! particle n-1 

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

! particle n

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
! }}}
! Torsional angle forces
! particles 1, 2, 3, n-2, n-1, and n are done outside of the loop
! particle 1

        I = 1
             COEF =(CD(I+1,1)+CD(I+1,2)*(12.0*DCOS(ANG(I+1,2))
     1         *DCOS(ANG(I+1,2))-3.0))
     1  *(1.0D0/DSQRT(X_PROD(I+1)*X_PROD(I)))  

        FTA_X(I) = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,1) +
     1         DOT_PROD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,1) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,1))) 

        FTA(I,2) = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,2) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,2) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,2))) 

        FTA_Z(I) = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,3) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,3) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,3))) 


! PARTICLE 2

        I = 2
        COEF =(CD(I+1,1)+CD(I+1,2)*(12.0*DCOS(ANG(I+1,2))
     1  *DCOS(ANG(I+1,2)) - 3.0))
     1        *(1.0D0/DSQRT(X_PROD(I+1)*X_PROD(I)))  

             COEF1 = (CD(I,1) + CD(I,2)*(12.0*DCOS(ANG(I,2))
     1        *DCOS(ANG(I,2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I)*X_PROD(I-1)))  

        A1 =  -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,1) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,1) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,1))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,1) +
     1        DOT_PROD(I,2)*DR(I,I+1,1) - DOT_PROD(I,2)*DR(I-1,I,1) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,1) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,1) - DOT_PROD(I-1,2)*DR(I,I+1,1) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,1) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,1))) 

        FTA_X(I) = A1 + A2 

        A1 = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,2) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,2) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,2))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,2) +
     1        DOT_PROD(I,2)*DR(I,I+1,2) - DOT_PROD(I,2)*DR(I-1,I,2) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,2) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,2) - DOT_PROD(I-1,2)*DR(I,I+1,2) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,2) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,2)))

        FTA(I,2) = A1 + A2 
        
        A1 = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,3) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,3) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,3))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,3) +
     1        DOT_PROD(I,2)*DR(I,I+1,3) - DOT_PROD(I,2)*DR(I-1,I,3) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,3) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,3) - DOT_PROD(I-1,2)*DR(I,I+1,3) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,3) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,3))) 

        FTA_Z(I) = A1 + A2 

! PARTICLE 3

        I = 3
        COEF=(CD(I+1,1)+CD(I+1,2)*(12.0*DCOS(ANG(I+1,2))
     1  *DCOS(ANG(I+1,2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I+1)*X_PROD(I)))  

        COEF1=(CD(I,1)+CD(I,2)*(12.0*DCOS(ANG(I,2))
     1        *DCOS(ANG(I,2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I)*X_PROD(I-1)))  

        COEF2=(CD(I-1,1)+CD(I-1,2)*(12.0*DCOS(ANG(I-1,2))
     1        *DCOS(ANG(I-1,2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-1)*X_PROD(I-2)))  

        A1 = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,1) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,1) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,1))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,1) +
     1        DOT_PROD(I,2)*DR(I,I+1,1) - DOT_PROD(I,2)*DR(I-1,I,1) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,1) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,1) - DOT_PROD(I-1,2)*DR(I,I+1,1) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,1) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,1))) 

        A3 = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,1) -
     1        DOT_PROD(I-2,2)*DR(I-1,I,1) + DOT_PROD(I-1,2)*DR(I-2,I-1,1) +
     1        DOT_PROD(I-1,1)*DR(I-2,I-1,1) - 2.0*DOT_PROD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,1) - DOT_PROD(I-1,2)*DR(I,I+1,1) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,1))) 

        FTA(I,1) = sum(AA(1:3))
 
        A1 = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,2) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,2) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,2))) 
        
        A2 = -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,2) +
     1        DOT_PROD(I,2)*DR(I,I+1,2) - DOT_PROD(I,2)*DR(I-1,I,2) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,2) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,2) - DOT_PROD(I-1,2)*DR(I,I+1,2) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,2) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,2))) 

        A3 = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,2) -
     1        DOT_PROD(I-2,2)*DR(I-1,I,2) + DOT_PROD(I-1,2)*DR(I-2,I-1,2) +
     1        DOT_PROD(I-1,1)*DR(I-2,I-1,2) - 2.0*DOT_PROD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,2) - DOT_PROD(I-1,2)*DR(I,I+1,2) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,2)))

        FTA(I,2) = sum(AA(1:3))
 
        A1 = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,3) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,3) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,3))) 

        A2 =  -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,3) +
     1        DOT_PROD(I,2)*DR(I,I+1,3) - DOT_PROD(I,2)*DR(I-1,I,3) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,3) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,3) - DOT_PROD(I-1,2)*DR(I,I+1,3) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,3) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,3))) 

        A3 = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,3) -
     1        DOT_PROD(I-2,2)*DR(I-1,I,3) + DOT_PROD(I-1,2)*DR(I-2,I-1,3) +
     1        DOT_PROD(I-1,1)*DR(I-2,I-1,3) - 2.0*DOT_PROD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,3) - DOT_PROD(I-1,2)*DR(I,I+1,3) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,3))) 

        FTA(I,3) = SUM(AA(1:3))

! PARTICLES 4 TO N-3

        DO I=1,N
           XPD(I)=1.0D0/DSQRT(X_PROD(I+1)*X_PROD(I)) 
        ENDDO

        DO I = 4, N-3
          DO K=1,4
            COEF(K) = CD(I-K+2,1) + CD(I-K+2,2)*(12.0*DCOS(ANG(I-K+2,2))**2-3.0)
            COEF(K) = COEF(K)*XPD(I+1-K)
          ENDDO

        AA(1) = -COEF(1)*(-DOT_PROD(I+1,2)*DR(I+1,I+2,1) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,1) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,1))) 

        AA(2) = -COEF(2)*(-DOT_PROD(I-1,2)*DR(I+1,I+2,1) +
     1        DOT_PROD(I,2)*DR(I,I+1,1) - DOT_PROD(I,2)*DR(I-1,I,1) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,1) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,1) - DOT_PROD(I-1,2)*DR(I,I+1,1) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,1) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,1))) 

        AA(3) = -COEF(3)*(DOT_PROD(I-2,2)*DR(I,I+1,1) -
     1        DOT_PROD(I-2,2)*DR(I-1,I,1) + DOT_PROD(I-1,2)*DR(I-2,I-1,1) +
     1        DOT_PROD(I-1,1)*DR(I-2,I-1,1) - 2.0*DOT_PROD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1  DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,1) - DOT_PROD(I-1,2)*DR(I,I+1,1) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,1))) 

        AA(4) = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,1) -
     1        DOT_PROD(I-2,1)*DR(I-3,I-2,1) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,1))) 

        FTA(I,1) =sum(AA)

        AA(1) = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,2) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,2) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,2))) 

        AA(2) = -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,2) +
     1        DOT_PROD(I,2)*DR(I,I+1,2) - DOT_PROD(I,2)*DR(I-1,I,2) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,2) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,2) - DOT_PROD(I-1,2)*DR(I,I+1,2) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,2) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,2))) 

        AA(3) = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,2) -
     1        DOT_PROD(I-2,2)*DR(I-1,I,2) + DOT_PROD(I-1,2)*DR(I-2,I-1,2) +
     1        DOT_PROD(I-1,1)*DR(I-2,I-1,2) - 2.0*DOT_PROD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,2) - DOT_PROD(I-1,2)*DR(I,I+1,2) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,2))) 

        AA(4) = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,2) -
     1        DOT_PROD(I-2,1)*DR(I-3,I-2,2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,2))) 

        FTA(I,2) = sum(AA) 

        AA(1) = -COEF*(-DOT_PROD(I+1,2)*DR(I+1,I+2,3) +
     1        DOT_PROD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,3) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,3))) 

        AA(2) = -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,3) +
     1        DOT_PROD(I,2)*DR(I,I+1,3) - DOT_PROD(I,2)*DR(I-1,I,3) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,3) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,3) - DOT_PROD(I-1,2)*DR(I,I+1,3) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,3) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,3))) 

        AA(3) = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,3) -
     1        DOT_PROD(I-2,2)*DR(I-1,I,3) + DOT_PROD(I-1,2)*DR(I-2,I-1,3) +
     1        DOT_PROD(I-1,1)*DR(I-2,I-1,3) - 2.0*DOT_PROD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,3) - DOT_PROD(I-1,2)*DR(I,I+1,3) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,3))) 
        
        AA(4) = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,3) -
     1        DOT_PROD(I-2,1)*DR(I-3,I-2,3) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,3))) 

        FTA(I,3) = SUM(AA)  

        ENDDO

! PARTICLE N-2

        I = N-2
        COEF(1)=(CD(I,1)+CD(I,2)*(12.0*DCOS(ANG(I,2))**2-3.0)*XPD(I-1)

        COEF(2)=(CD(I-1,1)+CD(I-1,2)*(12.0*DCOS(ANG(I-1,2))
     1        *DCOS(ANG(I-1,2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-1)*X_PROD(I-2)))  

        COEF3=(CD(I-2,1)+CD(I-2,2)*(12.0*DCOS(ANG(I-2,2))
     1        *DCOS(ANG(I-2,2)) - 
     1  3.0))*(1.0D0/DSQRT(X_PROD(I-2)*X_PROD(I-3)))  

        A1 = -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,1) + 
     1        DOT_PROD(I,2)*DR(I,I+1,1) - DOT_PROD(I,2)*DR(I-1,I,1) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,1) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,1) - DOT_PROD(I-1,2)*DR(I,I+1,1) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,1) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,1)))


        A2 = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,1) -
     1        DOT_PROD(I-2,2)*DR(I-1,I,1) + DOT_PROD(I-1,2)*DR(I-2,I-1,1) +
     1        DOT_PROD(I-1,1)*DR(I-2,I-1,1) - 2.0*DOT_PROD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,1) - DOT_PROD(I-1,2)*DR(I,I+1,1) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,1))) 
        
        A3 = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,1) -
     1        DOT_PROD(I-2,1)*DR(I-3,I-2,1) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,1))) 

        FTA(I,1) = sum(AA(1:3)) 

        A1 =  -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,2) +  
     1        DOT_PROD(I,2)*DR(I,I+1,2) - DOT_PROD(I,2)*DR(I-1,I,2) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,2) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,2) - DOT_PROD(I-1,2)*DR(I,I+1,2) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,2) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,2))) 

        A2 =  -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,2) -
     1        DOT_PROD(I-2,2)*DR(I-1,I,2) + DOT_PROD(I-1,2)*DR(I-2,I-1,2) +
     1        DOT_PROD(I-1,1)*DR(I-2,I-1,2) - 2.0*DOT_PROD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,2) - DOT_PROD(I-1,2)*DR(I,I+1,2) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,2)))

        A3 = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,2) -
     1        DOT_PROD(I-2,1)*DR(I-3,I-2,2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,2))) 

        FTA(I,2) = sum(AA(1:3)) 
 
        AA(1) = -COEF1*(-DOT_PROD(I-1,2)*DR(I+1,I+2,3) +  
     1        DOT_PROD(I,2)*DR(I,I+1,3) - DOT_PROD(I,2)*DR(I-1,I,3) -
     1        DOT_PROD(I,1)*DR(I+1,I+2,3) + 2.0*DOT_PROD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,3) - DOT_PROD(I-1,2)*DR(I,I+1,3) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*DR(I,I+1,3) +
     1        DOT_PROD(I,2)*DR(I+1,I+2,3))) 

        AA(2) = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,3) -
     1        DOT_PROD(I-2,2)*DR(I-1,I,3) + DOT_PROD(I-1,2)*DR(I-2,I-1,3) +
     1        DOT_PROD(I-1,1)*DR(I-2,I-1,3) - 2.0*DOT_PROD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,3) - DOT_PROD(I-1,2)*DR(I,I+1,3) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,3))) 

        AA(3) = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,3) -
     1        DOT_PROD(I-2,1)*DR(I-3,I-2,3) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,3))) 

        FTA(I,3) = sum(AA(1:3)) 

! PARTICLE N-1

        I = N-1
        COEF2=(CD(I-1,1)+CD(I-1,2)*(12.0*DCOS(ANG(I-1,2))
     1        *DCOS(ANG(I-1,2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-1)*X_PROD(I-2)))  

        COEF3=(CD(I-2,1)+CD(I-2,2)*(12.0*DCOS(ANG(I-2,2))
     1        *DCOS(ANG(I-2,2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-2)*X_PROD(I-3)))  

        A1 = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,1) - 
     1        DOT_PROD(I-2,2)*DR(I-1,I,1) +
     1        DOT_PROD(I-1,2)*DR(I-2,I-1,1) +  DOT_PROD(I-1,1)*DR(I-2,I-1,1) -
     1        2.0*DOT_PROD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,1) - DOT_PROD(I-1,2)*DR(I,I+1,1) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,1))) 

        A2 = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,1) - 
     1        DOT_PROD(I-2,1)*DR(I-3,I-2,1) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,1))) 

        FTA(I,1) = sum(AA(1:2))  

        A1 = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,2) - 
     1        DOT_PROD(I-2,2)*DR(I-1,I,2) +
     1        DOT_PROD(I-1,2)*DR(I-2,I-1,2) +  DOT_PROD(I-1,1)*DR(I-2,I-1,2) -
     1        2.0*DOT_PROD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,2) -
     1  DOT_PROD(I-1,1)*DR(I,I+1,2) - DOT_PROD(I-1,2)*DR(I,I+1,2) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,2))) 

        A2 = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,2) - 
     1        DOT_PROD(I-2,1)*DR(I-3,I-2,2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,2) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,2))) 

        FTA(I,2) = sum(AA(1:2))  

        A1 = -COEF2*(DOT_PROD(I-2,2)*DR(I,I+1,3) - 
     1        DOT_PROD(I-2,2)*DR(I-1,I,3) +
     1        DOT_PROD(I-1,2)*DR(I-2,I-1,3) +  DOT_PROD(I-1,1)*DR(I-2,I-1,3) -
     1        2.0*DOT_PROD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-1,1)*DR(I,I+1,3) - DOT_PROD(I-1,2)*DR(I,I+1,3) +
     1        DOT_PROD(I-1,2)*DR(I-1,I,3))) 

        A2 = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,3) - 
     1        DOT_PROD(I-2,1)*DR(I-3,I-2,3) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,3) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,3))) 

        FTA(I,3) = sum(AA(1:2))  
 
! PARTICLE N

        I = N
        COEF3=(CD(I-2,1)+CD(I-2,2)*(12.0*DCOS(ANG(I-2,2))
     1        *DCOS(ANG(I-2,2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-2)*X_PROD(I-3)))  

        FTA(I,1:3) = -COEF3*(DOT_PROD(I-3,2)*DR(I-2,I-1,1:3) 
     1        - DOT_PROD(I-2,1)*DR(I-3,I-2,1:3) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*DR(I-1,I,1) -
     1        DOT_PROD(I-2,2)*DR(I-2,I-1,1:3))) 

 ! TOTAL UP THE GRADIENTS

             F=FNB+FB+FBA+FTA 

        DO I = 1, N
          J = (I-1)*3
          FQ(J+1) = -F(I,1)
          FQ(J+2) = -F(I,2)
          FQ(J+3) = -F(I,3)
        ENDDO

        RETURN
        END
! }}}

!> @brief Calculate the second derivative matrix (two-sided numerical approach)

        SUBROUTINE CALC_DYN(N,QO,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
! {{{

include "blnvars.inc.f90"

! Fill in the Hessian matrix
! {{{

        DO J = 1, 3*N
            QO(J) = QO(J) + DELTA

            CALL CALC_INT_COORDS(N,QO,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
            CALL CALC_GRADIENT(N,QO,FQ2,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)

            QO(J) = QO(J) - 2.0*DELTA

            CALL CALC_INT_COORDS(N,QO,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)
            CALL CALC_GRADIENT(N,QO,FQ1,AB,CD,R,DR,DOT_PROD,X_PROD,ANG,RADII,NTYPE)

            QO(J) = QO(J) + DELTA
    
            DO I = J, 3*N
                HESS(I,J) = (FQ2(I) -  FQ1(I))/(2.0*DELTA)
                HESS(J,I) = HESS(I,J)
            ENDDO
        ENDDO

! }}}
        RETURN
        END
! }}}
!> Fill the parameter arrays

        SUBROUTINE PARAM_ARRAY(N,AB,CD,PTYPE)
! {{{
include blnvars.inc.f90
! Firstly, specify amino acid types by filling in array NTYPE(:)
! 1 => Hydrophobic (B);  2 => Hydrophilic (L); 3 => Neutral (N)
! {{{
        NTYPE(1:9) = 1
        NTYPE(10:12) = 3
        NTYPE(13:19:2) = 2
        NTYPE(14:20:2) = 1
        NTYPE(21:23) = 3
        NTYPE(24:32) = 1
        NTYPE(33:35) = 3
        NTYPE(36:46:2) = 2
        NTYPE(37:45:2) = 1
             ! }}}
! Specify parameters for the dihedral angle potential.
!       These are stored in arrays 
!       c_param(:), d_param(:) 
! {{{

        DO I = 1, N-3
          ICOUNT = 0
          DO J = 0,3
            IF(NTYPE(I+J) .EQ. 3) THEN
              ICOUNT = ICOUNT + 1
            ENDIF
          ENDDO

          IF(ICOUNT .GE. 2) THEN
            CD(I+1,1:2) = (/ 0.0 0.2*EPSILON /)
          ELSE
            CD(I+1,1:2) = 1.2*EPSILON
          ENDIF

          ICOUNT = 0
        ENDDO
!}}}
! Specify parameters for the L-J interaction between non-bonded particles.
!       These are stored in arrays 
!       a_param(:,:), b_param(:,:) 
! {{{

        DO I = 1, N
          DO J = 1, N
            IF (NTYPE(I) .EQ. 3 .OR. NTYPE(J) .EQ. 3)THEN
              AB(I,J,1:2) = (/ 1.0*EPSILON, 0.0D0  /) 
            ELSEIF (NTYPE(I) .EQ. 1 .AND. NTYPE(J) .EQ. 1)THEN
              AB(I,J,1:2)= EPSILON*(/ 1.0D0, -1.0D0 /) 
            ELSE
              AB(I,J,1:2) = EPSILON*2.0/3.0 
            ENDIF
          ENDDO
        ENDDO
!}}}
        RETURN
        END
! }}}

        ENDMODULE
