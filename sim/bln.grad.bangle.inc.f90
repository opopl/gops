
! bond angle forces: particle 1
! particles 1, 2, N-1, and N are done outside of the loop
! {{{
        I = 1
        DEN = SIN(ANG(I+1,1))*SQRT(DPD(I+1,1)*DPD(I,1))
        RNUM = RK_THETA*(ANG(I+1,1) - THETA_0)
        FBA(I,1:3)=-RNUM*((DPD(I,2)/DPD(I,1))*DR(I,I+1,1:3)-DR(I+1,I+2,1:3))/DEN

     ! PARTICLE 2

        I = 2
        DEN = DSIN(BOND_ANGLE(I))
     1        *DSQRT(DOT_PROD(I,1)*DOT_PROD(I-1,1))
        DEN1 = DSIN(BOND_ANGLE(I+1))*DSQRT(DOT_PROD(I+1,1)
     1         *DOT_PROD(I,1))

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*XR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *XR(I-1,I) + XR(I,I+1) - XR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*XR(I,I+1) - XR(I+1,I+2))/DEN1

        FBA_X(I) = A1 + A2 

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*YR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *YR(I-1,I) + YR(I,I+1) - YR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*YR(I,I+1) - YR(I+1,I+2))/DEN1

        FBA_Y(I) = A1 + A2 

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*ZR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *ZR(I-1,I) + ZR(I,I+1) - ZR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*ZR(I,I+1) - ZR(I+1,I+2))/DEN1

        FBA_Z(I) = A1 + A2 

! PARTICLES 3 THRU N-2 

        DO I = 3, N-2

        DEN = DSIN(BOND_ANGLE(I))*
     1              DSQRT(DOT_PROD(I,1)*DOT_PROD(I-1,1))
        DEN1 = DSIN(BOND_ANGLE(I+1))*
     1               DSQRT(DOT_PROD(I+1,1)*DOT_PROD(I,1))
        DEN2 = DSIN(BOND_ANGLE(I-1))*DSQRT(DOT_PROD(I-2,1)
     1         *DOT_PROD(I-1,1))

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*XR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *XR(I-1,I) + XR(I,I+1) - XR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*XR(I,I+1) - XR(I+1,I+2))/DEN1

        A3 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*XR(I-1,I) - XR(I-2,I-1))/DEN2

        FBA_X(I) = A1 + A2 + A3 

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*YR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *YR(I-1,I) + YR(I,I+1) - YR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*YR(I,I+1) - YR(I+1,I+2))/DEN1

        A3 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*YR(I-1,I) - YR(I-2,I-1))/DEN2

        FBA_Y(I) = A1 + A2 + A3 

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*ZR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *ZR(I-1,I) + ZR(I,I+1) - ZR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*ZR(I,I+1) - ZR(I+1,I+2))/DEN1

        A3 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*ZR(I-1,I) - ZR(I-2,I-1))/DEN2

        FBA_Z(I) = A1 + A2 + A3 

        ENDDO

! PARTICLE N-1 

        I = N-1
        DEN = DSIN(BOND_ANGLE(I))*
     1              DSQRT(DOT_PROD(I,1)*DOT_PROD(I-1,1))
        DEN1 = DSIN(BOND_ANGLE(I-1))*DSQRT(DOT_PROD(I-2,1)
     1         *DOT_PROD(I-1,1))

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*XR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *XR(I-1,I) + XR(I,I+1) - XR(I-1,I))/DEN

        A2 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*XR(I-1,I) - XR(I-2,I-1))/DEN1

        FBA_X(I) = A1 + A2

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*YR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *YR(I-1,I) + YR(I,I+1) - YR(I-1,I))/DEN

        A2 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*YR(I-1,I) - YR(I-2,I-1))/DEN1

        FBA_Y(I) = A1 + A2

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*ZR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *ZR(I-1,I) + ZR(I,I+1) - ZR(I-1,I))/DEN

        A2 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*ZR(I-1,I) - ZR(I-2,I-1))/DEN1

        FBA_Z(I) = A1 + A2

! PARTICLE N

        I = N
        DEN = DSIN(BOND_ANGLE(I-1))*DSQRT(DOT_PROD(I-2,1)
     1        *DOT_PROD(I-1,1))

        FBA_X(I) = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*
     1        ((DOT_PROD(I-2,2)/DOT_PROD(I-1,1))*XR(I-1,I) 
     1        - XR(I-2,I-1))/DEN

        FBA_Y(I) = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*
     1        ((DOT_PROD(I-2,2)/DOT_PROD(I-1,1))*YR(I-1,I) 
     1        - YR(I-2,I-1))/DEN

        FBA_Z(I) = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*
     1        ((DOT_PROD(I-2,2)/DOT_PROD(I-1,1))*ZR(I-1,I) 
     1        - ZR(I-2,I-1))/DEN
! }}}

