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

!include "bln.vars.inc.f90"

!!
!! Without these initialisations the NAG compiler fills in random numbers for
!! unassigned elements with optimisation turned on.
!!
      !ANG=0.0D0
      !COSTOR=0.0D0
      !DFAC=0.0D0
      !SINBOND=0.0D0
      !DPD=0.0D0
      !XPD=0.0D0
      !LEN_DR=0.0D0
      !DR=0.0D0

      !CALL CALC_INT_COORDS_BLN(N,QO,R,DR,DPD,XPD,ANG,LEN_DR,COSTOR,SINBOND,A,DFAC)

      !CALL CALC_ENERGY_BLN(N,QO,ENERGY,LJREP,LJATT,A,R,DR,DPD,XPD,ANG,LEN_DR,&
                !RK_R,RK_THETA,COSTOR)

      !IF (.NOT.GRADT) RETURN
 
      !CALL CALC_GRADIENT_BLN(N,QO,GRAD,LJREP,LJATT,A,R,DR,&
         !DPD,XPD,ANG,LEN_DR,RK_R,RK_THETA,COSTOR,DFAC,SINBOND)

      !RETURN
      !END
!! }}}

!!> @brief Calculate the internal coordinates of a BLN chain

      !SUBROUTINE CALC_INT_COORDS_BLN(N,QO,R,DR,DPD,XPD,ANG,LEN_DR, 
     !&                              COSTOR,SINBOND,A,DFAC)
!! {{{
!include "bln.vars.inc.f90"

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
            !LEN_DR(I,J) = SQRT(SUM(DR(I,J,1:3)**2))
            !LEN_DR(J,I) = LEN_DR(I,J)
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
         !DPD(I,K) = SUM(DR(I,I+1,1:3)*DR(J,J+1,1:3))
        !ENDDO
      !ENDDO

!! }}}

!! Cross-products between adjacent bond vectors {{{

      !DO I = 1, N-2
         !XPD(I) = DPD(I,1)*DPD(I+1,1) - DPD(I,2)*DPD(I,2)   
      !ENDDO
!! }}}
!! Bond angles {{{

      !DO I = 1, N-2
         !COS_THETA=-DPD(I,2)/(SQRT(DPD(I,1)*DPD(I+1,1)))
         !ANG(I+1,1) = ACOS(COS_THETA)
         !SINBOND(I+1)=SIN(ANG(I+1,1))*SQRT(DPD(I,1)*DPD(I+1,1))
      !ENDDO
!! }}}
!! Torsional angles {{{

      !DO I = 1, N-3
         !COS_PHI = DPD(I,2)*DPD(I+1,2) - DPD(I,3)*DPD(I+1,1)
         !COS_PHI = COS_PHI/SQRT(XPD(I)*XPD(I+1))
         !IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=SIGN(1.0D0,COS_PHI)
         !ANG(I+1,2) = ACOS(COS_PHI)
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
     !&              +C_BLN(I+1)*(12.0*COSTOR(I+1)**2-3.0))/SQRT(XPD(I+1)*XPD(I))
!! }}}
      !ENDDO

      !RETURN
      !END
!! }}}

!!> @brief Calculate the energy of a BLN chain

      !SUBROUTINE CALC_ENERGY_BLN(N,QO,ENERGY,LJREP,LJATT,A,ANG,LEN_DR,RK_R,RK_THETA,COSTOR)
!! {{{
!include "bln.vars.inc.f90"

      !E= 0.0D0  
      !DO I = 1, N-2
         !DO J = I+2, N
            !RAD(6) = LEN_DR(I,J)**6; RAD(12)=RAD(6)**2
            !E(1)=E(1)+LJREP(I,J)/RAD(12) + LJATT(I,J)/RAD(6)
         !ENDDO
      !ENDDO
      !E(1)=E(1)*4.0D0

      !DO I = 1, N-1
         !E(2) = E(2) + (LEN_DR(I,I+1)-1.0D0)**2
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

      !SUBROUTINE CALC_GRADIENT_BLN(N,QO,FQ,LJREP,LJATT,A,R,DR,DPD,XPD,
     !&                            ANG,LEN_DR,RK_R,RK_THETA,COSTOR,DFAC,SINBOND)
!! {{{
!include "bln.vars.inc.f90" 
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
            !RAD(7)=LEN_DR(I,J)**7
            !RAD(8)=LEN_DR(I,J)*RAD(7)
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
         !RVAR = 1.0D0/LEN_DR(I,I+1) 
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
      !FBA(1,1:3) = -RNUM*((DPD(1,2)/DPD(1,1))*DR(1,1+1,1:3) - DR(1+1,1+2,1:3))/DEN
      
!! }}}
!! particle 2
!! {{{
!!     i = 2
      !DEN = SINBOND(2)
      !DEN1 = SINBOND(3)

      !A1 = -RK_THETA*(ANG(2,1) - THETA_0)*( (DPD(2-1,2)/
     !1  DPD(2,1))*DR(2,2+1,1) - (DPD(2-1,2)/DPD(2-1,1))
     !1      *DR(2-1,2,1) + DR(2,2+1,1) - DR(2-1,2,1))/DEN

      !A2 = -RK_THETA*(ANG(2+1,1) - THETA_0)*((DPD(2,2)/DPD(2,1))*DR(2,2+1,1) - DR(2+1,2+2,1))/DEN1

      !FBA_X(2) = A1 + A2 

      !A1 = -RK_THETA*(ANG(2,1) - THETA_0)*( (DPD(2-1,2)/
     !1  DPD(2,1))*DR(2,2+1,2) - (DPD(2-1,2)/DPD(2-1,1))
     !1      *DR(2-1,2,2) + DR(2,2+1,2) - DR(2-1,2,2))/DEN

      !A2 = -RK_THETA*(ANG(2+1,1) - THETA_0)*((DPD(2,2)/DPD(2,1))*DR(2,2+1,2) - DR(2+1,2+2,2))/DEN1

      !FBA_Y(2) = A1 + A2 

      !A1 = -RK_THETA*(ANG(2,1) - THETA_0)*( (DPD(2-1,2)/
     !1  DPD(2,1))*DR(2,2+1,3) - (DPD(2-1,2)/DPD(2-1,1))
     !1      *DR(2-1,2,3) + DR(2,2+1,3) - DR(2-1,2,3))/DEN

      !A2 = -RK_THETA*(ANG(2+1,1) - THETA_0)*((DPD(2,2)/DPD(2,1))*DR(2,2+1,3) - DR(2+1,2+2,3))/DEN1

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

!include  "bln.vars.inc.f90"

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

