     MODULE BLN 

     IMPLICIT NONE

     CONTAINS

include bln.gen.f90
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
! }}}
!
!------------------------------------------------------------
! Doxygen - EBLN {{{
!
!> @name EBLN
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
        SUBROUTINE EBLN(N,QO,GRAD,ENERGY,GTEST,PTYPE)
! {{{
include bln.vars.inc.f90

        STEST=.FALSE.

        CALL PARAM_ARRAY(N,AB,CD,PTYPE)
        CALL CALC_INT_COORDS(N,QO,R,DR,BVR,LEN_DR,LEN_BVR,DPD,XPD_2,ANG)
        CALL CALC_ENERGY(N,ENERGY,AB,CD,LEN_DR,LEN_BVR,ANG)

        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN

        CALL CALC_GRADIENT(N,QO,GRAD,AB,CD,R,DR,DPD,XPD,ANG,LEN_DR,NTYPE)

        IF (.NOT.STEST) RETURN

        CALL CALC_DYN(N,QO,AB,CD,R,DR,DPD,XPD,ANG,LEN_DR,NTYPE)

        RETURN
        END
! }}}
!> Calculate the internal coordinates

        SUBROUTINE CALC_INT_COORDS(N,QO,R,DR,BVR,LEN_DR,LEN_BVR,DPD,XPD_2,ANG)
! {{{
include bln.vars.inc.f90

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
            LEN_DR(I,J) = SQRT(SUM(DR(I,J,1:3)**2))
            LEN_DR(J,I) = LEN_DR(I,J)
          ENDDO
            BVR(I,1:3)=DR(I,I+1,1:3)
        ENDDO

! DOT PRODUCTS BETWEEN BOND VECTORS

      DO I = 1, N-1
        IF ( I .LE. N-3 ) THEN KMAX=3
        IF ( I .EQ. N-2 ) THEN KMAX=2
        IF ( I .LE. N-1 ) THEN KMAX=1
        DO J=1,KMAX
         J=I+K-1
         DPD(I,K) = SUM(BVR(I,1:3)*BVR(J,1:3))
        ENDDO
        LEN_BV(I)=DPD(I,1)
      ENDDO

! Squared cross-products between adjacent bond vectors i and i+1 

        DO I = 1, N-2
           XPD_2(I) = DPD(I,1)*DPD(I+1,1)-DPD(I,2)**2 
        ENDDO

! BOND ANGLES

        DO I = 1, N-2
            COS_THETA=-DPD(I,2)/(LEN_BV(I)*LEN_BV(I+1))
            ANG(I+1,1) = ACOS(COS_THETA)
        ENDDO

! TORSIONAL ANGLES

        DO I = 1, N-3
            COS_PHI = (DPD(I,2)*DPD(I+1,2)-DPD(I,3)*DPD(I+1,1))
            COS_PHI = COS_PHI/SQRT(XPD_2(I)*XPD_2(I+1))
            IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=COS_PHI/ABS(COS_PHI)
            ANG(I+1,2) = ACOS(COS_PHI)
        ENDDO

        RETURN
        END

! }}} 
!> Calculate the energy
        SUBROUTINE CALC_ENERGY(N,ENERGY,AB,CD,LEN_DR,LEN_BVR,ANG)
! {{{
include bln.vars.inc.f90  

       S(6)=SIGMA**6; E=0.0D0 
        
! non-bonded
        DO I = 1, N-2
          DO J = I+2, N
            RAD(6)=LEN_DR(I,J)**6
            E(1) = E(1) + 4.0*((AB(I,J,1)*S(12)/(RAD(12))) + (AB(I,J,2)*S(6)/RAD(6)))
          ENDDO
        ENDDO

! bonded
        DO I = 1, N-1
          E(2) = E(2) + 0.5*RK_R*(LEN_BVR(I)-SIGMA)**2
        ENDDO

! bond angles 
        DO I = 2, N-1
          E(3) = E(3) + 0.5*RK_THETA*(ANG(I,1)-THETA_0)**2
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

        SUBROUTINE CALC_GRADIENT(N,QO,FQ,AB,CD,R,DR,DPD,XPD,ANG,LEN_DR,NTYPE)
! {{{
include bln.vars.inc.f90 

      S(6)=S(1)**6 ; S(12)=S(6)**2

      FNB= 0.0D0; FB= 0.0D0 ; FBA= 0.0D0 ; FTA= 0.0D0 ; F= 0.0D0 

! ..... Non-bonded interaction forces ..... {{{

        DO I = 1, N-2
          DO J = I+2, N
  
            RAD(1)=LEN_DR(I,J) ; RAD(7)=RAD(1)**7 ; RAD(14)=RAD(7)**2 ; RAD(8)=RAD(7)*RAD(1)
        
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
           RVAR = SIGMA/LEN_DR(I,I+1) 

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
        DEN = DSIN(ANG(I+1,1))*SQRT(DPD(I+1,1)*DPD(I,1))
        RNUM = RK_THETA*(ANG(I+1,1) - THETA_0)
        FBA(I,1:3)=-RNUM*((DPD(I,2)/DPD(I,1))*DR(I,I+1,1:3)-DR(I+1,I+2,1:3))/DEN

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
     1  *(1.0D0/SQRT(XPD(I+1)*XPD(I)))  

        FTA_X(I) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,1) +
     1         DPD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        FTA(I,2) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,2) +
     1        DPD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        FTA_Z(I) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,3) +
     1        DPD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 


! PARTICLE 2

        I = 2
        COEF =(CD(I+1,1)+CD(I+1,2)*(12.0*DCOS(ANG(I+1,2))
     1  *DCOS(ANG(I+1,2)) - 3.0))
     1        *(1.0D0/SQRT(XPD(I+1)*XPD(I)))  

             COEF1 = (CD(I,1) + CD(I,2)*(12.0*DCOS(ANG(I,2))
     1        *DCOS(ANG(I,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I)*XPD(I-1)))  

        A1 =  -COEF*(-DPD(I+1,2)*DR(I+1,I+2,1) +
     1        DPD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,1) +
     1        DPD(I,2)*DR(I,I+1,1) - DPD(I,2)*DR(I-1,I,1) -
     1        DPD(I,1)*DR(I+1,I+2,1) + 2.0*DPD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        FTA_X(I) = A1 + A2 

        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,2) +
     1        DPD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,2) +
     1        DPD(I,2)*DR(I,I+1,2) - DPD(I,2)*DR(I-1,I,2) -
     1        DPD(I,1)*DR(I+1,I+2,2) + 2.0*DPD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2)))

        FTA(I,2) = A1 + A2 
        
        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,3) +
     1        DPD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,3) +
     1        DPD(I,2)*DR(I,I+1,3) - DPD(I,2)*DR(I-1,I,3) -
     1        DPD(I,1)*DR(I+1,I+2,3) + 2.0*DPD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        FTA_Z(I) = A1 + A2 

! PARTICLE 3

        I = 3
        COEF=(CD(I+1,1)+CD(I+1,2)*(12.0*DCOS(ANG(I+1,2))
     1  *DCOS(ANG(I+1,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I+1)*XPD(I)))  

        COEF1=(CD(I,1)+CD(I,2)*(12.0*DCOS(ANG(I,2))
     1        *DCOS(ANG(I,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I)*XPD(I-1)))  

        COEF2=(CD(I-1,1)+CD(I-1,2)*(12.0*DCOS(ANG(I-1,2))
     1        *DCOS(ANG(I-1,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-1)*XPD(I-2)))  

        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,1) +
     1        DPD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,1) +
     1        DPD(I,2)*DR(I,I+1,1) - DPD(I,2)*DR(I-1,I,1) -
     1        DPD(I,1)*DR(I+1,I+2,1) + 2.0*DPD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        A3 = -COEF2*(DPD(I-2,2)*DR(I,I+1,1) -
     1        DPD(I-2,2)*DR(I-1,I,1) + DPD(I-1,2)*DR(I-2,I-1,1) +
     1        DPD(I-1,1)*DR(I-2,I-1,1) - 2.0*DPD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1))) 

        FTA(I,1) = sum(AA(1:3))
 
        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,2) +
     1        DPD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 
        
        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,2) +
     1        DPD(I,2)*DR(I,I+1,2) - DPD(I,2)*DR(I-1,I,2) -
     1        DPD(I,1)*DR(I+1,I+2,2) + 2.0*DPD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        A3 = -COEF2*(DPD(I-2,2)*DR(I,I+1,2) -
     1        DPD(I-2,2)*DR(I-1,I,2) + DPD(I-1,2)*DR(I-2,I-1,2) +
     1        DPD(I-1,1)*DR(I-2,I-1,2) - 2.0*DPD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)))

        FTA(I,2) = sum(AA(1:3))
 
        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,3) +
     1        DPD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        A2 =  -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,3) +
     1        DPD(I,2)*DR(I,I+1,3) - DPD(I,2)*DR(I-1,I,3) -
     1        DPD(I,1)*DR(I+1,I+2,3) + 2.0*DPD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        A3 = -COEF2*(DPD(I-2,2)*DR(I,I+1,3) -
     1        DPD(I-2,2)*DR(I-1,I,3) + DPD(I-1,2)*DR(I-2,I-1,3) +
     1        DPD(I-1,1)*DR(I-2,I-1,3) - 2.0*DPD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3))) 

        FTA(I,3) = SUM(AA(1:3))

! PARTICLES 4 TO N-3

        DO I=1,N
           IXPD(I)=1.0D0/SQRT(XPD(I+1)*XPD(I)) 
        ENDDO

        DO I = 4, N-3
          DO K=1,4
            COEF(K) = CD(I-K+2,1) + CD(I-K+2,2)*(12.0*DCOS(ANG(I-K+2,2))**2-3.0)
            COEF(K) = COEF(K)*IXPD(I+1-K)
          ENDDO

        AA(1) = -COEF(1)*(-DPD(I+1,2)*DR(I+1,I+2,1) +
     1        DPD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        AA(2) = -COEF(2)*(-DPD(I-1,2)*DR(I+1,I+2,1) +
     1        DPD(I,2)*DR(I,I+1,1) - DPD(I,2)*DR(I-1,I,1) -
     1        DPD(I,1)*DR(I+1,I+2,1) + 2.0*DPD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        AA(3) = -COEF(3)*(DPD(I-2,2)*DR(I,I+1,1) -
     1        DPD(I-2,2)*DR(I-1,I,1) + DPD(I-1,2)*DR(I-2,I-1,1) +
     1        DPD(I-1,1)*DR(I-2,I-1,1) - 2.0*DPD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1  DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1))) 

        AA(4) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,1) -
     1        DPD(I-2,1)*DR(I-3,I-2,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1))) 

        FTA(I,1) =sum(AA)

        AA(1) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,2) +
     1        DPD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        AA(2) = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,2) +
     1        DPD(I,2)*DR(I,I+1,2) - DPD(I,2)*DR(I-1,I,2) -
     1        DPD(I,1)*DR(I+1,I+2,2) + 2.0*DPD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        AA(3) = -COEF2*(DPD(I-2,2)*DR(I,I+1,2) -
     1        DPD(I-2,2)*DR(I-1,I,2) + DPD(I-1,2)*DR(I-2,I-1,2) +
     1        DPD(I-1,1)*DR(I-2,I-1,2) - 2.0*DPD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2))) 

        AA(4) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,2) -
     1        DPD(I-2,1)*DR(I-3,I-2,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2))) 

        FTA(I,2) = sum(AA) 

        AA(1) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,3) +
     1        DPD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        AA(2) = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,3) +
     1        DPD(I,2)*DR(I,I+1,3) - DPD(I,2)*DR(I-1,I,3) -
     1        DPD(I,1)*DR(I+1,I+2,3) + 2.0*DPD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        AA(3) = -COEF2*(DPD(I-2,2)*DR(I,I+1,3) -
     1        DPD(I-2,2)*DR(I-1,I,3) + DPD(I-1,2)*DR(I-2,I-1,3) +
     1        DPD(I-1,1)*DR(I-2,I-1,3) - 2.0*DPD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3))) 
        
        AA(4) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,3) -
     1        DPD(I-2,1)*DR(I-3,I-2,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3))) 

        FTA(I,3) = SUM(AA)  

        ENDDO

! PARTICLE N-2

        I = N-2
        COEF(1)=(CD(I,1)+CD(I,2)*(12.0*DCOS(ANG(I,2))**2-3.0)*IXPD(I-1)

        COEF(2)=(CD(I-1,1)+CD(I-1,2)*(12.0*DCOS(ANG(I-1,2))
     1        *DCOS(ANG(I-1,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-1)*XPD(I-2)))  

        COEF3=(CD(I-2,1)+CD(I-2,2)*(12.0*DCOS(ANG(I-2,2))
     1        *DCOS(ANG(I-2,2)) - 
     1  3.0))*(1.0D0/SQRT(XPD(I-2)*XPD(I-3)))  

        A1 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,1) + 
     1        DPD(I,2)*DR(I,I+1,1) - DPD(I,2)*DR(I-1,I,1) -
     1        DPD(I,1)*DR(I+1,I+2,1) + 2.0*DPD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1)))


        A2 = -COEF2*(DPD(I-2,2)*DR(I,I+1,1) -
     1        DPD(I-2,2)*DR(I-1,I,1) + DPD(I-1,2)*DR(I-2,I-1,1) +
     1        DPD(I-1,1)*DR(I-2,I-1,1) - 2.0*DPD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1))) 
        
        A3 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,1) -
     1        DPD(I-2,1)*DR(I-3,I-2,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1))) 

        FTA(I,1) = sum(AA(1:3)) 

        A1 =  -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,2) +  
     1        DPD(I,2)*DR(I,I+1,2) - DPD(I,2)*DR(I-1,I,2) -
     1        DPD(I,1)*DR(I+1,I+2,2) + 2.0*DPD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        A2 =  -COEF2*(DPD(I-2,2)*DR(I,I+1,2) -
     1        DPD(I-2,2)*DR(I-1,I,2) + DPD(I-1,2)*DR(I-2,I-1,2) +
     1        DPD(I-1,1)*DR(I-2,I-1,2) - 2.0*DPD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)))

        A3 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,2) -
     1        DPD(I-2,1)*DR(I-3,I-2,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2))) 

        FTA(I,2) = sum(AA(1:3)) 
 
        AA(1) = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,3) +  
     1        DPD(I,2)*DR(I,I+1,3) - DPD(I,2)*DR(I-1,I,3) -
     1        DPD(I,1)*DR(I+1,I+2,3) + 2.0*DPD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        AA(2) = -COEF2*(DPD(I-2,2)*DR(I,I+1,3) -
     1        DPD(I-2,2)*DR(I-1,I,3) + DPD(I-1,2)*DR(I-2,I-1,3) +
     1        DPD(I-1,1)*DR(I-2,I-1,3) - 2.0*DPD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3))) 

        AA(3) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,3) -
     1        DPD(I-2,1)*DR(I-3,I-2,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3))) 

        FTA(I,3) = sum(AA(1:3)) 

! PARTICLE N-1

        I = N-1
        COEF2=(CD(I-1,1)+CD(I-1,2)*(12.0*DCOS(ANG(I-1,2))
     1        *DCOS(ANG(I-1,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-1)*XPD(I-2)))  

        COEF3=(CD(I-2,1)+CD(I-2,2)*(12.0*DCOS(ANG(I-2,2))
     1        *DCOS(ANG(I-2,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-2)*XPD(I-3)))  

        A1 = -COEF2*(DPD(I-2,2)*DR(I,I+1,1) - 
     1        DPD(I-2,2)*DR(I-1,I,1) +
     1        DPD(I-1,2)*DR(I-2,I-1,1) +  DPD(I-1,1)*DR(I-2,I-1,1) -
     1        2.0*DPD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1))) 

        A2 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,1) - 
     1        DPD(I-2,1)*DR(I-3,I-2,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1))) 

        FTA(I,1) = sum(AA(1:2))  

        A1 = -COEF2*(DPD(I-2,2)*DR(I,I+1,2) - 
     1        DPD(I-2,2)*DR(I-1,I,2) +
     1        DPD(I-1,2)*DR(I-2,I-1,2) +  DPD(I-1,1)*DR(I-2,I-1,2) -
     1        2.0*DPD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,2) -
     1  DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2))) 

        A2 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,2) - 
     1        DPD(I-2,1)*DR(I-3,I-2,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2))) 

        FTA(I,2) = sum(AA(1:2))  

        A1 = -COEF2*(DPD(I-2,2)*DR(I,I+1,3) - 
     1        DPD(I-2,2)*DR(I-1,I,3) +
     1        DPD(I-1,2)*DR(I-2,I-1,3) +  DPD(I-1,1)*DR(I-2,I-1,3) -
     1        2.0*DPD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3))) 

        A2 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,3) - 
     1        DPD(I-2,1)*DR(I-3,I-2,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3))) 

        FTA(I,3) = sum(AA(1:2))  
 
! PARTICLE N

        I = N
        COEF3=(CD(I-2,1)+CD(I-2,2)*(12.0*DCOS(ANG(I-2,2))
     1        *DCOS(ANG(I-2,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-2)*XPD(I-3)))  

        FTA(I,1:3) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,1:3) 
     1        - DPD(I-2,1)*DR(I-3,I-2,1:3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1:3))) 

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

! Doxygen - CALC_DYN {{{
!> @name CALC_DYN
!> @brief Calculate the second derivative matrix (two-sided numerical approach)
!
!> @param[in] N        INTEGER         Number of particles
!> @param[in] QO
! }}}

        SUBROUTINE CALC_DYN(N,QO,AB,CD,R,DR,DPD,XPD,ANG,LEN_DR,NTYPE)
include bln.vars.inc.f90
include bln.hess.inc.f90
        RETURN
        END

! Doxygen - PARAM_ARRAY {{{
!> @name PARAM_ARRAY
!>
!> @brief Fill the parameter arrays for the wild-type (WT) BLN model 
!
!> @param[in]   N  INTEGER              Number of particles 
!> @param[out] AB  DP, DIMENSION(N,N,2) Parameters for the L-J interaction between non-bonded particles. 
!> @param[out] CD  DP, DIMENSION(N,2)   Parameters for the dihedral angle potential. 
!
! }}}

        SUBROUTINE PARAM_ARRAY(N,AB,CD,PTYPE)

        include bln.vars.inc.f90
        include bln.ntype.inc.f90

        SELECTCASE(PTYPE)
                CASE("GO")
			include bln.go.connect.inc.f90
			include bln.go.ab.inc.f90
			include bln.go.cd.inc.f90
                CASE("WT")
                        include bln.wt.ab.inc.f90       ! L-J interaction between non-bonded particles.
                        include bln.go.cd.inc.f90       ! Dihedral angle potential
        ENDSELECT

        RETURN
        END

        ENDMODULE
