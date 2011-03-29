!
!   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
!   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
!   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
!   (AT YOUR OPTION) ANY LATER VERSION.
!
!   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
!   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
!   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
!   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
!
!   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
!   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
!   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
!
!
!  ENERGY AND GRADIENT FOR THE STOCKMAYER POTENTIAL USING TWO POLAR
!  COORDINATES FOR THE DIPOLE DIRECTION
!
      SUBROUTINE STOCK(NATOMS,X,V,ESTOCK,GTEST,STEST)
      USE KEY,ONLY : STOCKMU, STOCKLAMBDA, STOCKSPIN, STOCKZTOL, STOCKMAXSPIN
      USE MODHESS
      IMPLICIT NONE
      LOGICAL GTEST,STEST
      LOGICAL ZALIGNTEST
      INTEGER J1, J2, K1, K2, NAT2, J3, J4, REALNATOMS, NATOMS, OFFSET, SPINITER
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), ESTOCK
      DOUBLE PRECISION DPRAND, S, SIGMA1, SIGMA2, THETA1, THETA2, Q(4), PI, M(3,3)
      DOUBLE PRECISION MU2, DUMMY
      DOUBLE PRECISION X1, X2, Y1, Y2, Z1, Z2, T1, P1, T2, P2
      DOUBLE PRECISION R12, R122, R125, R126, R127, R128, R129, R1211, R1212, R1214, R1216
      DOUBLE PRECISION N1DOTN2, N1DOTR12, N2DOTR12, CP1P2, SP1PP2, SP1MP2
      DOUBLE PRECISION CT1, CT2, ST1, ST2, CP1, CP2, SP1, SP2, DX, DY, DZ, DX2, DY2, DZ2

      REALNATOMS=NATOMS/2
      OFFSET = 3*REALNATOMS
      MU2=STOCKMU**2
      ESTOCK=0.0D0
      V(1:6*REALNATOMS)=0.0D0

      IF (STOCKSPIN) THEN
         PI = ATAN2(1.0D0, 1.0D0) * 4.0D0
!        CHECK THAT NO SPINS ARE TOO CLOSELY ALIGNED WITH Z, SINCE THAT WOULD
!        INTRODUCE REDUNDANT PHI ANGLES.
         SPINITER = 1
         DO
            IF (.NOT.ZALIGNTEST(NATOMS, X)) THEN
               EXIT
            END IF
            IF (SPINITER.GT.STOCKMAXSPIN) THEN
               PRINT*,'WARNING: RANDOMISATION OF ORIENTATION FAILED TO REMOVE Z-ALIGNED DIPOLES'
               EXIT
            END IF

!           RANDOM QUATERNION
            S = DPRAND();
            SIGMA1 = SQRT(1.0D0-S);
            SIGMA2 = SQRT(S);
            THETA1 = 2.0D0*PI*DPRAND();
            THETA2 = 2.0*PI*DPRAND();
            Q(1) = COS(THETA2) * SIGMA2;
            Q(2) = SIN(THETA1) * SIGMA1;
            Q(3) = COS(THETA1) * SIGMA1;
            Q(4) = SIN(THETA2) * SIGMA2;

!           ROTATION MATRIX CORRESPONDING TO THE QUATERNION
            M(1,1) = Q(1)*Q(1) + Q(2)*Q(2) - Q(3)*Q(3) - Q(4)*Q(4)
            M(2,1) = 2.0*(Q(2)*Q(3) + Q(1)*Q(4))
            M(3,1) = 2.0*(Q(2)*Q(4) - Q(1)*Q(3))
            M(1,2) = 2.0*(Q(2)*Q(3) - Q(1)*Q(4))
            M(2,2) = Q(1)*Q(1) - Q(2)*Q(2) + Q(3)*Q(3) - Q(4)*Q(4)
            M(3,2) = 2.0*(Q(3)*Q(4) + Q(1)*Q(2))
            M(1,3) = 2.0*(Q(2)*Q(4) + Q(1)*Q(3))
            M(2,3) = 2.0*(Q(3)*Q(4) - Q(1)*Q(2))
            M(3,3) = Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) + Q(4)*Q(4)

            CALL NEWROTGEOMSTOCK(NATOMS, X, M, 0.0D0, 0.0D0, 0.0D0)

            SPINITER = SPINITER + 1
         END DO

         IF ( (SPINITER.GT.1).AND.(SPINITER.LE.STOCKMAXSPIN) ) THEN
            PRINT*,'WARNING: ORIENTATION HAS BEEN RANDOMISED TO REMOVE Z-ALIGNED DIPOLES'
         END IF
      END IF
      
      DO J1=1,REALNATOMS
         J3=3*J1
         X1=X(J3-2)
         Y1=X(J3-1)
         Z1=X(J3)
         T1=X(OFFSET+J3-2)
         CT1=COS(T1)
         ST1=SIN(T1)
         P1=X(OFFSET+J3-1)
         CP1=COS(P1)
         SP1=SIN(P1)
         DO J2=J1+1,REALNATOMS
            J4=3*J2
            X2=X(J4-2)
            Y2=X(J4-1)
            Z2=X(J4)
            T2=X(OFFSET+J4-2)
            CT2=COS(T2)
            ST2=SIN(T2)
            P2=X(OFFSET+J4-1)
            CP2=COS(P2)
            SP2=SIN(P2)
            N1DOTR12=ST1*CP1*(X1-X2)+ST1*SP1*(Y1-Y2)+CT1*(Z1-Z2)
            N2DOTR12=ST2*CP2*(X1-X2)+ST2*SP2*(Y1-Y2)+CT2*(Z1-Z2)
            N1DOTN2=ST1*CP1*ST2*CP2 + ST1*SP1*ST2*SP2 + CT1*CT2
            R122=(X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2
            R12=SQRT(R122)
            R126=R122**3
            R127=R126*R12
            R129=R127*R122
            R1212=R126**2
            R125=R126/R12
            R1214=R1212*R122
           
            DUMMY= (MU2*R127*(-3*N1DOTR12*N2DOTR12 + N1DOTN2*R122) +
     &             (4 - 4*STOCKLAMBDA*R126))/R1212
         
            ESTOCK=ESTOCK+DUMMY

! DERIVATIVES FOR POSITIONS

      DUMMY = (-3*(-((-16 + 8*STOCKLAMBDA*R126 + 
     &   MU2*(5*N1DOTR12*N2DOTR12*R127 - N1DOTN2*R129))*(X1 - X2))
     &   + MU2*N2DOTR12*R129*CP1*ST1 + MU2*N1DOTR12*R129*CP2*ST2))/R1214
      V(J3-2)=V(J3-2)+DUMMY
      V(J4-2)=V(J4-2)-DUMMY

      DUMMY = (-3*(-((-16 + 8*STOCKLAMBDA*R126 + 
     &   MU2*(5*N1DOTR12*N2DOTR12*R127 - N1DOTN2*R129))*(Y1 - Y2))
     &   + MU2*N2DOTR12*R129*SP1*ST1 + MU2*N1DOTR12*R129*SP2*ST2))/R1214
      V(J3-1)=V(J3-1)+DUMMY
      V(J4-1)=V(J4-1)-DUMMY

      DUMMY = (-3*(-((-16 + 8*STOCKLAMBDA*R126 + 
     &   MU2*(5*N1DOTR12*N2DOTR12*R127 - N1DOTN2*R129))*(Z1 - Z2))
     &   + MU2*N2DOTR12*R129*CT1 + MU2*N1DOTR12*R129*CT2))/R1214
      V(J3)=V(J3)+DUMMY
      V(J4)=V(J4)-DUMMY

! DERIVATIVES FOR ANGULAR VARIABLES OF ATOM J1

      V(OFFSET+J3-2) = V(OFFSET+J3-2) +
     &   (MU2*((3*N2DOTR12*Z1 - 3*N2DOTR12*Z2 - R122*CT2)*ST1 + 
     &   CP1*CT1*(3*N2DOTR12*(-X1 + X2) + R122*CP2*ST2) + 
     &   CT1*SP1*(-3*N2DOTR12*Y1 + 3*N2DOTR12*Y2 + R122*SP2*ST2)))/R125

      V(OFFSET+J3-1) = V(OFFSET+J3-1) +
     &   (MU2*ST1*(SP1*(3*N2DOTR12*(X1 - X2) - R122*CP2*ST2) + 
     &   CP1*(3*N2DOTR12*(-Y1 + Y2) + R122*SP2*ST2)))/R125

! DERIVATIVES FOR ANGULAR VARIABLES OF ATOM J2

      V(OFFSET+J4-2) = V(OFFSET+J4-2) +
     &   (MU2*(CP2*CT2*(3*N1DOTR12*(-X1 + X2) + R122*CP1*ST1) + 
     &   CT2*SP2*(-3*N1DOTR12*Y1 + 3*N1DOTR12*Y2 + R122*SP1*ST1) + 
     &   (3*N1DOTR12*Z1 - 3*N1DOTR12*Z2 - R122*CT1)*ST2))/R125

      V(OFFSET+J4-1) = V(OFFSET+J4-1) +
     &   (MU2*(SP2*(3*N1DOTR12*(X1 - X2) - R122*CP1*ST1) + 
     &   CP2*(3*N1DOTR12*(-Y1 + Y2) + R122*SP1*ST1))*ST2)/R125

         ENDDO
      ENDDO

!     ***********************
!     ANALYTIC HESSIAN MATRIX
!     ***********************
      IF (STEST) THEN
         HESS(1:3*NATOMS, 1:3*NATOMS) = 0.0D0
         DO J1=1, REALNATOMS
            J3 = J1*3
!           COORDINATES AND FUNCTIONS OF FIRST ATOM
            X1 = X(J3 - 2)
            Y1 = X(J3 - 1)
            Z1 = X(J3)
            T1 = X(OFFSET + J3 - 2)
            P1 = X(OFFSET + J3 - 1)
            ST1 = SIN(T1)
            CT1 = COS(T1)
            SP1 = SIN(P1)
            CP1 = COS(P1)

            DO J2=1, REALNATOMS
               IF (J1 == J2) CYCLE
               J4 = J2*3

!              COORDINATES AND FUNCTIONS OF SECOND ATOM
               X2 = X(J4 - 2)
               Y2 = X(J4 - 1)
               Z2 = X(J4)
               T2 = X(OFFSET + J4 - 2)
               P2 = X(OFFSET + J4 - 1)
               ST2 = SIN(T2)
               CT2 = COS(T2)
               SP2 = SIN(P2)
               CP2 = COS(P2)

!              COORDINATES AND FUNCTIONS RELATED TO BOTH ATOMS
               DX = X1 - X2
               DY = Y1 - Y2
               DZ = Z1 - Z2
               DX2 = DX * DX
               DY2 = DY * DY
               DZ2 = DZ * DZ
               N1DOTR12 = ST1*CP1*DX + ST1*SP1*DY + CT1*DZ
               N2DOTR12 = ST2*CP2*DX + ST2*SP2*DY + CT2*DZ
               N1DOTN2 = ST1*CP1*ST2*CP2 + ST1*SP1*ST2*SP2 + CT1*CT2
               R122 = DX2 + DY2 + DZ2
               R12 = SQRT(R122)
               R125 = R122*R122*R12
               R126 = R125*R12
               R127 = R126*R12
               R128 = R126*R122
               R129 = R127*R122
               R1211 = R126*R125
               R1212 = R126*R126
               R1214 = R1212*R122
               R1216 = R128*R128
               SP1PP2 = SIN(P1 + P2)
               SP1MP2 = SIN(P1 - P2)
               CP1P2 = COS(P1 - P2)

!              [1] THE FIVE COMPLETELY DIAGONAL TERMS: SAME ATOM, SAME COORDINATE
!              X1,X1
               HESS(J3-2, J3-2) = HESS(J3-2, J3-2) +
     &            (-3*(MU2*N1DOTN2*R1211 + 16*R122 - 8*STOCKLAMBDA*R128 +
     &            DX2*(-224 + 64*STOCKLAMBDA*R126 + 35*MU2*N1DOTR12*N2DOTR12*R127 -
     &            5*MU2*N1DOTN2*R129) -
     &            MU2*R129*(5*N1DOTR12*(N2DOTR12 + 2*CP2*DX*ST2) +
     &            2*CP1*ST1*(5*DX*N2DOTR12 - CP2*R122*ST2))))/R1216
!              Y1,Y1
               HESS(J3-1, J3-1) = HESS(J3-1, J3-1) +
     &            (3*(-(MU2*N1DOTN2*R1211) - 16*R122 + 8*STOCKLAMBDA*R128 +
     &            DY2*(224 - 64*STOCKLAMBDA*R126 +
     &            5*MU2*(-7*N1DOTR12*N2DOTR12*R127 + N1DOTN2*R129)) +
     &            MU2*R129*(5*N1DOTR12*(N2DOTR12 + 2*DY*SP2*ST2) +
     &            2*SP1*ST1*(5*DY*N2DOTR12 - R122*SP2*ST2))))/R1216
!              Z1,Z1
               HESS(J3, J3) = HESS(J3, J3) +
     &            (-3*(MU2*N1DOTN2*R1211 + 16*R122 +
     &            DZ2*(-224 + 64*STOCKLAMBDA*R126 + 35*MU2*N1DOTR12*N2DOTR12*R127) -
     &            8*STOCKLAMBDA*R128) + 3*MU2*(5*
     &            (DZ2*N1DOTN2 + 2*CT2*DZ*N1DOTR12 + (2*CT1*DZ + N1DOTR12)*N2DOTR12)
     &            - 2*CT1*CT2*R122)*R129)/R1216
!              T1,T1
               HESS(OFFSET+J3-2, OFFSET+J3-2) = HESS(OFFSET+J3-2, OFFSET+J3-2) +
     &            (MU2*(3*CT1*DZ*N2DOTR12 - CT1*CT2*R122 + 3*CP1*DX*N2DOTR12*ST1 +
     &            3*DY*N2DOTR12*SP1*ST1 - CP1P2*R122*ST1*ST2))/R125
!              P1,P1
               HESS(OFFSET+J3-1, OFFSET+J3-1) = HESS(OFFSET+J3-1, OFFSET+J3-1) +
     &            (MU2*ST1*(3*CP1*DX*N2DOTR12 + 3*DY*N2DOTR12*SP1 - CP1P2*R122*ST2))/R125

!              [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCKS: SAME ATOM, DIFFERENT COORDINATES
!              X1,Y1
               DUMMY =
     &            (3*(MU2*R129*(5*CP1*DY*N2DOTR12*ST1 + 5*CP2*DY*N1DOTR12*ST2 -
     &            R122*SP1PP2*ST1*ST2) +
     &            DX*(DY*(224 - 64*STOCKLAMBDA*R126 + 5*MU2*(-7*N1DOTR12*N2DOTR12*R127 +
     &            N1DOTN2*R129)) + 5*MU2*R129*(N2DOTR12*SP1*ST1 + N1DOTR12*SP2*ST2))))/
     &            R1216
               HESS(J3-2, J3-1) = HESS(J3-2, J3-1) + DUMMY
               HESS(J3-1, J3-2) = HESS(J3-1, J3-2) + DUMMY
!              X1,Z1
               DUMMY =
     &            (-3*(DX*DZ*(-224 + 64*STOCKLAMBDA*R126 + 35*MU2*N1DOTR12*N2DOTR12*R127) -
     &            5*DX*MU2*(DZ*N1DOTN2 + CT2*N1DOTR12 + CT1*N2DOTR12)*R129 +
     &            MU2*R129*(CP1*(-5*DZ*N2DOTR12 + CT2*R122)*ST1 +
     &            CP2*(-5*DZ*N1DOTR12 + CT1*R122)*ST2)))/R1216
               HESS(J3-2, J3) = HESS(J3-2, J3) + DUMMY
               HESS(J3, J3-2) = HESS(J3, J3-2) + DUMMY
!              Y1,Z1
               DUMMY =
     &            (3*(DY*DZ*(224 - 64*STOCKLAMBDA*R126 - 35*MU2*N1DOTR12*N2DOTR12*R127) +
     &            5*DY*MU2*(DZ*N1DOTN2 + CT2*N1DOTR12 + CT1*N2DOTR12)*R129 +
     &            MU2*R129*(-(R122*(CT2*SP1*ST1 + CT1*SP2*ST2)) +
     &            5*DZ*(N2DOTR12*SP1*ST1 + N1DOTR12*SP2*ST2))))/R1216
               HESS(J3-1, J3) = HESS(J3-1, J3) + DUMMY
               HESS(J3, J3-1) = HESS(J3, J3-1) + DUMMY
!              X1,T1
               DUMMY =
     &            (-3*MU2*(ST1*(5*DX*DZ*N2DOTR12 - CT2*DX*R122 - CP2*DZ*R122*ST2) +
     &            CT1*SP1*(-5*DX*DY*N2DOTR12 + CP2*DY*R122*ST2 + DX*R122*SP2*ST2) +
     &            CP1*CT1*(-5*DX2*N2DOTR12 + R122*(N2DOTR12 + 2*CP2*DX*ST2))))/R127
               HESS(J3-2, OFFSET+J3-2) = HESS(J3-2, OFFSET+J3-2) + DUMMY
               HESS(OFFSET+J3-2, J3-2) = HESS(OFFSET+J3-2, J3-2) + DUMMY
!              X1,P1
               DUMMY =
     &            (3*MU2*ST1*(5*CP1*DX*DY*N2DOTR12 - CP1*R122*(CP2*DY + DX*SP2)*ST2 +
     &            SP1*(-5*DX2*N2DOTR12 + R122*(N2DOTR12 + 2*CP2*DX*ST2))))/R127
               HESS(J3-2, OFFSET+J3-1) = HESS(J3-2, OFFSET+J3-1) + DUMMY
               HESS(OFFSET+J3-1, J3-2) = HESS(OFFSET+J3-1, J3-2) + DUMMY
!              Y1,T1
               DUMMY =
     &            (-3*MU2*(CP1*CT1*(-5*DX*DY*N2DOTR12 + CP2*DY*R122*ST2 + DX*R122*SP2*ST2) +
     &            ST1*(5*DY*DZ*N2DOTR12 - CT2*DY*R122 - DZ*R122*SP2*ST2) +
     &            CT1*SP1*(-5*DY2*N2DOTR12 + R122*(N2DOTR12 + 2*DY*SP2*ST2))))/R127
               HESS(J3-1, OFFSET+J3-2) = HESS(J3-1, OFFSET+J3-2) + DUMMY
               HESS(OFFSET+J3-2, J3-1) = HESS(OFFSET+J3-2, J3-1) + DUMMY
!              Y1,P1
               DUMMY =
     &            (3*MU2*ST1*(N2DOTR12*(5*CP1*DY2 - CP1*R122 - 5*DX*DY*SP1) +
     &            R122*(CP2*DY*SP1 - 2*CP1*DY*SP2 + DX*SP1*SP2)*ST2))/R127
               HESS(J3-1, OFFSET+J3-1) = HESS(J3-1, OFFSET+J3-1) + DUMMY
               HESS(OFFSET+J3-1, J3-1) = HESS(OFFSET+J3-1, J3-1) + DUMMY
!              Z1,T1
               DUMMY =
     &            (-3*MU2*((5*DZ2*N2DOTR12 - (2*CT2*DZ + N2DOTR12)*R122)*ST1 +
     &            CP1*CT1*(-5*DX*DZ*N2DOTR12 + CT2*DX*R122 + CP2*DZ*R122*ST2) +
     &            CT1*SP1*(-5*DY*DZ*N2DOTR12 + CT2*DY*R122 + DZ*R122*SP2*ST2)))/R127
               HESS(J3, OFFSET+J3-2) = HESS(J3, OFFSET+J3-2) + DUMMY
               HESS(OFFSET+J3-2, J3) = HESS(OFFSET+J3-2, J3) + DUMMY
!              Z1,P1
               DUMMY =
     &            (3*MU2*ST1*((5*DZ*N2DOTR12 - CT2*R122)*(CP1*DY - DX*SP1) +
     &            DZ*R122*SP1MP2*ST2))/R127
               HESS(J3, OFFSET+J3-1) = HESS(J3, OFFSET+J3-1) + DUMMY
               HESS(OFFSET+J3-1, J3) = HESS(OFFSET+J3-1, J3) + DUMMY
!              T1,P1
               DUMMY =
     &            -((CT1*MU2*(3*CP1*DY*N2DOTR12 - 3*DX*N2DOTR12*SP1 + R122*SP1MP2*ST2))/R125)
               HESS(OFFSET+J3-2, OFFSET+J3-1) = HESS(OFFSET+J3-2, OFFSET+J3-1) + DUMMY
               HESS(OFFSET+J3-1, OFFSET+J3-2) = HESS(OFFSET+J3-1, OFFSET+J3-2) + DUMMY

!              [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT PARTICLE, SAME COORDINATE
!              X1,X2
               HESS(J3-2, J4-2) =
     &            (3*(MU2*N1DOTN2*R1211 + 16*R122 - 8*STOCKLAMBDA*R128 +
     &            DX2*(-224 + 64*STOCKLAMBDA*R126 + 35*MU2*N1DOTR12*N2DOTR12*R127 -
     &            5*MU2*N1DOTN2*R129) - MU2*R129*(5*N1DOTR12*(N2DOTR12 + 2*CP2*DX*ST2) +
     &            2*CP1*ST1*(5*DX*N2DOTR12 - CP2*R122*ST2))))/R1216
!              Y1,Y2
               HESS(J3-1, J4-1) =
     &            (3*(-224*DY2 + MU2*N1DOTN2*R1211 + 16*R122 + 64*DY2*STOCKLAMBDA*R126 +
     &            35*DY2*MU2*N1DOTR12*N2DOTR12*R127 - 8*STOCKLAMBDA*R128 -
     &            5*MU2*(DY2*N1DOTN2 + N1DOTR12*N2DOTR12)*R129 +
     &            2*MU2*R129*(R122*SP1*SP2*ST1*ST2 - 5*DY*(N2DOTR12*SP1*ST1 +
     &            N1DOTR12*SP2*ST2))))/R1216
!              Z1,Z2
               HESS(J3, J4) =
     &            (3*(-224*DZ2 + MU2*N1DOTN2*R1211 + 16*R122 + 64*DZ2*STOCKLAMBDA*R126 +
     &            35*DZ2*MU2*N1DOTR12*N2DOTR12*R127 - 8*STOCKLAMBDA*R128 -
     &            MU2*(5*(DZ2*N1DOTN2 + 2*CT2*DZ*N1DOTR12 + (2*CT1*DZ + N1DOTR12)*
     &            N2DOTR12) - 2*CT1*CT2*R122)*R129))/R1216
!              T1,T2
               HESS(OFFSET+J3-2, OFFSET+J4-2) =
     &            (MU2*(-3*(CP1*CT1*DX + CT1*DY*SP1 - DZ*ST1)*(CP2*CT2*DX + CT2*DY*SP2 -
     &            DZ*ST2) + R122*(CP1P2*CT1*CT2 + ST1*ST2)))/R125
!              P1,P2
               HESS(OFFSET+J3-1, OFFSET+J4-1) =
     &            (MU2*(CP1P2*R122 - 3*(CP1*DY - DX*SP1)*(CP2*DY - DX*SP2))*ST1*ST2)/R125

!              [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT PARTICLE, DIFFERENT COORDINATE
!              X1,Y2 AND Y1,X2
               HESS(J3-2, J4-1) =
     &            (3*(DX*DY*(-224 + 64*STOCKLAMBDA*R126 + 35*MU2*N1DOTR12*N2DOTR12*R127 -
     &            5*MU2*N1DOTN2*R129) - MU2*R129*(5*N2DOTR12*(CP1*DY + DX*SP1)*ST1 +
     &            (5*CP2*DY*N1DOTR12 + 5*DX*N1DOTR12*SP2 - R122*SP1PP2*ST1)*ST2)))/R1216
               HESS(J3-1, J4-2) = HESS(J3-2, J4-1)
!              X1,Z2 AND Z1,X2
               HESS(J3-2, J4) =
     &            (3*(DX*DZ*(-224 + 64*STOCKLAMBDA*R126 + 35*MU2*N1DOTR12*N2DOTR12*R127) -
     &            5*DX*MU2*(DZ*N1DOTN2 + CT2*N1DOTR12 + CT1*N2DOTR12)*R129 +
     &            MU2*R129*(CP1*(-5*DZ*N2DOTR12 + CT2*R122)*ST1 +
     &            CP2*(-5*DZ*N1DOTR12 + CT1*R122)*ST2)))/R1216
               HESS(J3, J4-2) = HESS(J3-2, J4)
!              Y1,Z2 AND Z1,Y2
               HESS(J3-1, J4) =
     &            (3*(DY*DZ*(-224 + 64*STOCKLAMBDA*R126 + 35*MU2*N1DOTR12*N2DOTR12*R127) -
     &            5*DY*MU2*(DZ*N1DOTN2 + CT2*N1DOTR12 + CT1*N2DOTR12)*R129 +
     &            MU2*R129*(R122*(CT2*SP1*ST1 + CT1*SP2*ST2) -
     &            5*DZ*(N2DOTR12*SP1*ST1 + N1DOTR12*SP2*ST2))))/R1216
               HESS(J3, J4-1) = HESS(J3-1, J4)
!              X1,T2
               HESS(J3-2, OFFSET+J4-2) =
     &            (-3*MU2*(CT2*SP2*(-5*DX*DY*N1DOTR12 + CP1*DY*R122*ST1 + DX*R122*SP1*ST1) +
     &            CP2*CT2*(-5*DX2*N1DOTR12 + R122*(N1DOTR12 + 2*CP1*DX*ST1)) +
     &            (5*DX*DZ*N1DOTR12 - CT1*DX*R122 - CP1*DZ*R122*ST1)*ST2))/R127
!              X1,P2
               HESS(J3-2, OFFSET+J4-1) =
     &            (-3*MU2*(-(N1DOTR12*(5*CP2*DX*DY + (-5*DX2 + R122)*SP2)) +
     &            R122*(CP1*CP2*DY + CP2*DX*SP1 - 2*CP1*DX*SP2)*ST1)*ST2)/R127
!              Y1,T2
               HESS(J3-1, OFFSET+J4-2) =
     &            (-3*MU2*(CP2*CT2*(-5*DX*DY*N1DOTR12 + CP1*DY*R122*ST1 + DX*R122*SP1*ST1) +
     &            CT2*SP2*(-5*DY2*N1DOTR12 + R122*(N1DOTR12 + 2*DY*SP1*ST1)) +
     &            (5*DY*DZ*N1DOTR12 - CT1*DY*R122 - DZ*R122*SP1*ST1)*ST2))/R127
!              Y1,P2
               HESS(J3-1, OFFSET+J4-1) =
     &            (3*MU2*(N1DOTR12*(5*CP2*DY2 - CP2*R122 - 5*DX*DY*SP2) +
     &            R122*(-2*CP2*DY*SP1 + CP1*DY*SP2 + DX*SP1*SP2)*ST1)*ST2)/R127
!              Z1,T2
               HESS(J3, OFFSET+J4-2) =
     &            (-3*MU2*(CP2*CT2*(-5*DX*DZ*N1DOTR12 + CT1*DX*R122 + CP1*DZ*R122*ST1) +
     &            CT2*SP2*(-5*DY*DZ*N1DOTR12 + CT1*DY*R122 + DZ*R122*SP1*ST1) +
     &            (5*DZ2*N1DOTR12 - (2*CT1*DZ + N1DOTR12)*R122)*ST2))/R127
!              Z1,P2
               HESS(J3, OFFSET+J4-1) =
     &            (3*MU2*((5*DZ*N1DOTR12 - CT1*R122)*(CP2*DY - DX*SP2) - DZ*R122*SP1MP2*ST1)*
     &            ST2)/R127
!              T1,X2
               HESS(OFFSET+J3-2, J4-2) =
     &            (3*MU2*(-(ST1*(-5*DX*DZ*N2DOTR12 + CT2*DX*R122 + CP2*DZ*R122*ST2)) +
     &            CT1*SP1*(-5*DX*DY*N2DOTR12 + CP2*DY*R122*ST2 + DX*R122*SP2*ST2) +
     &            CP1*CT1*(-5*DX2*N2DOTR12 + R122*(N2DOTR12 + 2*CP2*DX*ST2))))/R127
!              P1,X2
               HESS(OFFSET+J3-1, J4-2) =
     &            (3*MU2*ST1*(-(N2DOTR12*(5*CP1*DX*DY + (-5*DX2 + R122)*SP1)) +
     &            R122*(CP1*CP2*DY - 2*CP2*DX*SP1 + CP1*DX*SP2)*ST2))/R127
!              T1,Y2
               HESS(OFFSET+J3-2, J4-1) =
     &            (3*MU2*(CP1*CT1*(-5*DX*DY*N2DOTR12 + CP2*DY*R122*ST2 + DX*R122*SP2*ST2) +
     &            ST1*(5*DY*DZ*N2DOTR12 - CT2*DY*R122 - DZ*R122*SP2*ST2) +
     &            CT1*SP1*(-5*DY2*N2DOTR12 + R122*(N2DOTR12 + 2*DY*SP2*ST2))))/R127
!              P1,Y2
               HESS(OFFSET+J3-1, J4-1) =
     &            (3*MU2*ST1*(5*DX*DY*N2DOTR12*SP1 - R122*SP1*(CP2*DY + DX*SP2)*ST2 +
     &            CP1*(-5*DY2*N2DOTR12 + R122*(N2DOTR12 + 2*DY*SP2*ST2))))/R127
!              T1,Z2
               HESS(OFFSET+J3-2, J4) =
     &            (3*MU2*((5*DZ2*N2DOTR12 - (2*CT2*DZ + N2DOTR12)*R122)*ST1 +
     &            CP1*CT1*(-5*DX*DZ*N2DOTR12 + CT2*DX*R122 + CP2*DZ*R122*ST2) +
     &            CT1*SP1*(-5*DY*DZ*N2DOTR12 + CT2*DY*R122 + DZ*R122*SP2*ST2)))/R127
!              P1,Z2
               HESS(OFFSET+J3-1, J4) =
     &            (3*MU2*ST1*(-((5*DZ*N2DOTR12 - CT2*R122)*(CP1*DY - DX*SP1)) -
     &            DZ*R122*SP1MP2*ST2))/R127
!              T1,P2
               HESS(OFFSET+J3-2, OFFSET+J4-1) =
     &            (MU2*(CT1*R122*SP1MP2 + 3*(-(CP2*DY) + DX*SP2)*(CP1*CT1*DX + CT1*DY*SP1 -
     &            DZ*ST1))*ST2)/R125
!              P1,T2
               HESS(OFFSET+J3-1, OFFSET+J4-2) =
     &            (MU2*ST1*(-(CT2*R122*SP1MP2) + 3*(-(CP1*DY) + DX*SP1)*
     &            (CP2*CT2*DX + CT2*DY*SP2 - DZ*ST2)))/R125

            END DO

         END DO
      ENDIF

      RETURN
      END
C
C  ORTHOGONALISE VEC1 TO OVERALL TRANSLATIONS AND ROTATIONS.
C  STILL NOT WORKING PROPERLY FOR NO APPARENT REASON!
C
      SUBROUTINE ORTHOGSTOCK(VEC1,COORDS,OTEST) 
      USE COMMONS
      USE KEY
      USE VECCK
      IMPLICIT NONE
      INTEGER J2, J3, NCHECK
      DOUBLE PRECISION COORDS(*), VEC1(*), DUMMY1, DUMMY2, ROOTN, VDOT, DUMMYX, DUMMYY, DUMMYZ, DDOT,
     1                 CMX, CMY, CMZ, AMASS(NATOMS), TMASS, RMASS(NATOMS), THETA, PHI, DUMMY,
     &                 ETX(3*NATOMS), ETY(3*NATOMS), ETZ(3*NATOMS), ERX(3*NATOMS), ERY(3*NATOMS), ERZ(3*NATOMS)
      LOGICAL OTEST

      IF (.NOT.ALLOCATED(VECCHK)) ALLOCATE(VECCHK(3*NATOMS,30))
      ROOTN=SQRT(1.0D0*(NATOMS/2)) ! FACTOR OF 2 FOR STOCKMAYER
      NCHECK=0
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      TMASS=0.0D0
      DO J2=1,(NATOMS/2) ! FACTOR OF 2 FOR STOCKMAYER
         AMASS(J2)=1.0D0
         RMASS(J2)=1.0D0
         IF (MASST) AMASS(J2)=ATMASS(J2)
         IF (MASST) RMASS(J2)=SQRT(ATMASS(J2))
         TMASS=TMASS+AMASS(J2)
      ENDDO
C
C  IF MASST THEN THE COORDINATES ALREADY HAVE A SQUARE ROOT OF THE MASS IN THEM
C
      DO J2=1,(NATOMS/2) ! FACTOR OF 2 FOR STOCKMAYER
         CMX=CMX+COORDS(3*(J2-1)+1)*RMASS(J2)
         CMY=CMY+COORDS(3*(J2-1)+2)*RMASS(J2)
         CMZ=CMZ+COORDS(3*(J2-1)+3)*RMASS(J2)
      ENDDO
      CMX=CMX/TMASS
      CMY=CMY/TMASS
      CMZ=CMZ/TMASS
C
C  ORTHOGONALIZE TO KNOWN EIGENVECTORS CORRESPONDING TO NEGATIVE EIGENVALUES
C  FOR CHECKINDEX RUN WITH NOHESS.
C
      IF (CHECKINDEX) THEN
         DO J2=1,NMDONE
            DUMMY1=0.0D0
            DO J3=1,NOPT
               DUMMY1=DUMMY1+VECCHK(J3,J2)*VEC1(J3)
            ENDDO
            VDOT=MAX(VDOT,ABS(DUMMY1))
            DUMMY2=0.0D0
            DO J3=1,NOPT
               VEC1(J3)=VEC1(J3)-DUMMY1*VECCHK(J3,J2)
               DUMMY2=DUMMY2+VEC1(J3)**2
            ENDDO
            DUMMY2=1.0D0/DSQRT(DUMMY2)
            DO J3=1,NOPT
               VEC1(J3)=VEC1(J3)*DUMMY2
            ENDDO
         ENDDO
         IF (OTEST) CALL VECNORM(VEC1,NOPT)
      ENDIF
!
!  GRAM-SCMIDT ORTHOGONALISATION OF OVERALL TRANS/ROT EIGENVECTORS.
!
      ETX(1:3*NATOMS)=0.0D0; ETY(1:3*NATOMS)=0.0D0; ETZ(1:3*NATOMS)=0.0D0
      ERX(1:3*NATOMS)=0.0D0; ERY(1:3*NATOMS)=0.0D0; ERZ(1:3*NATOMS)=0.0D0

      DO J2=1,(NATOMS/2) 
         ETX(3*(J2-1)+1)=1.0D0/ROOTN
         ETY(3*(J2-1)+2)=1.0D0/ROOTN
         ETZ(3*(J2-1)+3)=1.0D0/ROOTN
      ENDDO

      DUMMYX=0.0D0
      DUMMYY=0.0D0
      DUMMYZ=0.0D0
      DO J2=1,(NATOMS/2) ! ROTATION 
         J3=3*J2
         ERX(3*(J2-1)+2)=  COORDS(J3)  -CMZ
         ERX(3*(J2-1)+3)=-(COORDS(J3-1)-CMY)
         DUMMYX=DUMMYX+ERX(3*(J2-1)+2)**2+ERX(3*(J2-1)+3)**2
         ERY(3*(J2-1)+1)=  COORDS(J3)  -CMZ
         ERY(3*(J2-1)+3)=-(COORDS(J3-2)-CMX)
         DUMMYY=DUMMYY+ERY(3*(J2-1)+1)**2+ERY(3*(J2-1)+3)**2
         ERZ(3*(J2-1)+1)=  COORDS(J3-1)-CMY
         ERZ(3*(J2-1)+2)=-(COORDS(J3-2)-CMX)
         DUMMYZ=DUMMYZ+ERZ(3*(J2-1)+1)**2+ERZ(3*(J2-1)+2)**2
      ENDDO
      DO J2=(NATOMS/2)+1,NATOMS ! XROTATION FOR THETA, PHI COORDINATES
         THETA=COORDS(3*(J2-1)+1)
         PHI=COORDS(3*(J2-1)+2)
         DUMMY=TAN(THETA)
         ERX(3*(J2-1)+1)=SIN(PHI)
         IF (DUMMY.NE.0.0D0) ERX(3*(J2-1)+2)=COS(PHI)/DUMMY
         DUMMYX=DUMMYX+ERX(3*(J2-1)+1)**2+ERX(3*(J2-1)+2)**2
         ERY(3*(J2-1)+1)=COS(PHI)
         IF (DUMMY.NE.0.0D0) ERY(3*(J2-1)+2)=-SIN(PHI)/DUMMY
         DUMMYY=DUMMYY+ERY(3*(J2-1)+1)**2+ERY(3*(J2-1)+2)**2
         ERZ(3*(J2-1)+2)=-1.0D0
         DUMMYZ=DUMMYZ+ERZ(3*(J2-1)+2)**2
       ENDDO 
       DUMMYX=SQRT(DUMMYX); DUMMYY=SQRT(DUMMYY); DUMMYZ=SQRT(DUMMYZ)
       DO J2=1,3*NATOMS
          ERX(J2)=ERX(J2)/DUMMYX
          ERY(J2)=ERY(J2)/DUMMYY
          ERZ(J2)=ERZ(J2)/DUMMYZ
       ENDDO
!
! ORTHOGONALISE ERY TO ERX
!
       DUMMY=0.0D0
       DO J2=1,3*NATOMS
          DUMMY=DUMMY+ERX(J2)*ERY(J2)
       ENDDO
       DUMMY2=0.0D0
       DO J2=1,3*NATOMS
          ERY(J2)=ERY(J2)-DUMMY*ERX(J2)
          DUMMY2=DUMMY2+ERY(J2)**2
       ENDDO
       DUMMY2=SQRT(DUMMY2)
       DO J2=1,3*NATOMS
          ERY(J2)=ERY(J2)/DUMMY2
       ENDDO
!
! ORTHOGONALISE ERZ TO ERX
!
       DUMMY=0.0D0
       DO J2=1,3*NATOMS
          DUMMY=DUMMY+ERX(J2)*ERZ(J2)
       ENDDO
       DUMMY2=0.0D0
       DO J2=1,3*NATOMS
          ERZ(J2)=ERZ(J2)-DUMMY*ERX(J2)
          DUMMY2=DUMMY2+ERZ(J2)**2
       ENDDO
       DUMMY2=SQRT(DUMMY2)
       DO J2=1,3*NATOMS
          ERZ(J2)=ERZ(J2)/DUMMY2
       ENDDO
!
! ORTHOGONALISE ERZ TO ERY
!
       DUMMY=0.0D0
       DO J2=1,3*NATOMS
          DUMMY=DUMMY+ERY(J2)*ERZ(J2)
       ENDDO
       DUMMY2=0.0D0
       DO J2=1,3*NATOMS
          ERZ(J2)=ERZ(J2)-DUMMY*ERY(J2)
          DUMMY2=DUMMY2+ERZ(J2)**2
       ENDDO
       DUMMY2=SQRT(DUMMY2)
       DO J2=1,3*NATOMS
          ERZ(J2)=ERZ(J2)/DUMMY2
       ENDDO
!
!  CHECK ORTHOGONALISATION OF TRANS/ROT
!
!       DUMMY=DDOT(3*NATOMS,ETX,1,ETX,1)
!       PRINT *,'ETX,ETX=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETY,1,ETY,1)
!       PRINT *,'ETY,ETY=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETZ,1,ETZ,1)
!       PRINT *,'ETZ,ETZ=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ERX,1,ERX,1)
!       PRINT *,'ERX,ERX=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ERY,1,ERY,1)
!       PRINT *,'ERY,ERY=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ERZ,1,ERZ,1)
!       PRINT *,'ERZ,ERZ=',DUMMY
!
!       DUMMY=DDOT(3*NATOMS,ETX,1,ETY,1)
!       PRINT *,'ETX,ETY=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETX,1,ETZ,1)
!       PRINT *,'ETX,ETZ=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETX,1,ERX,1)
!       PRINT *,'ETX,ERX=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETX,1,ERY,1)
!       PRINT *,'ETX,ERY=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETX,1,ERZ,1)
!       PRINT *,'ETX,ERZ=',DUMMY
!
!       DUMMY=DDOT(3*NATOMS,ETY,1,ETZ,1)
!       PRINT *,'ETY,ETZ=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETY,1,ERX,1)
!       PRINT *,'ETY,ERX=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETY,1,ERY,1)
!       PRINT *,'ETY,ERY=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETY,1,ERZ,1)
!       PRINT *,'ETY,ERZ=',DUMMY
!    
!       DUMMY=DDOT(3*NATOMS,ETZ,1,ERX,1)
!       PRINT *,'ETZ,ERX=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETZ,1,ERY,1)
!       PRINT *,'ETZ,ERY=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETZ,1,ERZ,1)
!       PRINT *,'ETZ,ERZ=',DUMMY
!
!       DUMMY=DDOT(3*NATOMS,ERX,1,ERY,1)
!       PRINT *,'ERX,ERY=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ERX,1,ERZ,1)
!       PRINT *,'ERX,ERZ=',DUMMY
!
!       DUMMY=DDOT(3*NATOMS,ERY,1,ERZ,1)
!       PRINT *,'ERY,ERZ=',DUMMY
!
!  PROJECT TRANS/ROT SET OUT OF VEC1
!
       DUMMY=DDOT(3*NATOMS,ETX,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ETX(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ETY,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ETY(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ETZ,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ETZ(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ERX,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ERX(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ERY,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ERY(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ERZ,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ERZ(J2)
       ENDDO
!
! THE SIXTH COORDINATE FOR EACH PARTICLE IS A DUMMY.
!
      DO J2=(NATOMS/2)+1,NATOMS
         VEC1(3*(J2-1)+3)=0.0D0
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      RETURN
      END

C
      SUBROUTINE SHIFTSTOCK(COORDS,NATOMS)
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, NOPT, NATOMS, NIADD
      DOUBLE PRECISION CMX, CMY, CMZ, COORDS(3*NATOMS), DUMMYX, DUMMYY, DUMMYZ, DDOT, THETA, PHI, DUMMY, VEC6(3*NATOMS),
     &                 ETX(3*NATOMS), ETY(3*NATOMS), ETZ(3*NATOMS), ERX(3*NATOMS), ERY(3*NATOMS), ERZ(3*NATOMS), ROOTN

      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO J2=1,(NATOMS/2)
         CMX=CMX+COORDS(3*(J2-1)+1)
         CMY=CMY+COORDS(3*(J2-1)+2)
         CMZ=CMZ+COORDS(3*(J2-1)+3)
      ENDDO
      CMX=CMX/(NATOMS/2)
      CMY=CMY/(NATOMS/2)
      CMZ=CMZ/(NATOMS/2)
      ETX(1:3*NATOMS)=0.0D0; ETY(1:3*NATOMS)=0.0D0; ETZ(1:3*NATOMS)=0.0D0
      ERX(1:3*NATOMS)=0.0D0; ERY(1:3*NATOMS)=0.0D0; ERZ(1:3*NATOMS)=0.0D0

      ROOTN=SQRT(1.0D0*(NATOMS/2)) ! FACTOR OF 2 FOR STOCKMAYER
      DO J2=1,(NATOMS/2) 
         ETX(3*(J2-1)+1)=1.0D0/ROOTN
         ETY(3*(J2-1)+2)=1.0D0/ROOTN
         ETZ(3*(J2-1)+3)=1.0D0/ROOTN
      ENDDO

      DUMMYX=0.0D0
      DUMMYY=0.0D0
      DUMMYZ=0.0D0
      DO J2=1,(NATOMS/2) ! ROTATION 
         J3=3*J2
         ERX(3*(J2-1)+2)=  COORDS(J3)  -CMZ
         ERX(3*(J2-1)+3)=-(COORDS(J3-1)-CMY)
         DUMMYX=DUMMYX+ERX(3*(J2-1)+2)**2+ERX(3*(J2-1)+3)**2
         ERY(3*(J2-1)+1)=  COORDS(J3)  -CMZ
         ERY(3*(J2-1)+3)=-(COORDS(J3-2)-CMX)
         DUMMYY=DUMMYY+ERY(3*(J2-1)+1)**2+ERY(3*(J2-1)+3)**2
         ERZ(3*(J2-1)+1)=  COORDS(J3-1)-CMY
         ERZ(3*(J2-1)+2)=-(COORDS(J3-2)-CMX)
         DUMMYZ=DUMMYZ+ERZ(3*(J2-1)+1)**2+ERZ(3*(J2-1)+2)**2
      ENDDO
      DO J2=(NATOMS/2)+1,NATOMS ! XROTATION FOR THETA, PHI COORDINATES
         THETA=COORDS(3*(J2-1)+1)
         PHI=COORDS(3*(J2-1)+2)
         DUMMY=TAN(THETA)
         ERX(3*(J2-1)+1)=SIN(PHI)
         IF (DUMMY.NE.0.0D0) ERX(3*(J2-1)+2)=COS(PHI)/DUMMY
         DUMMYX=DUMMYX+ERX(3*(J2-1)+1)**2+ERX(3*(J2-1)+2)**2
         ERY(3*(J2-1)+1)=COS(PHI)
         IF (DUMMY.NE.0.0D0) ERY(3*(J2-1)+2)=-SIN(PHI)/DUMMY
         DUMMYY=DUMMYY+ERY(3*(J2-1)+1)**2+ERY(3*(J2-1)+2)**2
         ERZ(3*(J2-1)+2)=-1.0D0
         DUMMYZ=DUMMYZ+ERZ(3*(J2-1)+2)**2
      ENDDO 
      DUMMYX=SQRT(DUMMYX); DUMMYY=SQRT(DUMMYY); DUMMYZ=SQRT(DUMMYZ)
      DO J2=1,3*NATOMS
         ERX(J2)=ERX(J2)/DUMMYX
         ERY(J2)=ERY(J2)/DUMMYY
         ERZ(J2)=ERZ(J2)/DUMMYZ
      ENDDO
!
! THE SIXTH COORDINATE FOR EACH PARTICLE IS A DUMMY.
!
      DO J2=(NATOMS/2)+1,NATOMS
         HESS(3*(J2-1)+3,3*(J2-1)+3)=HESS(3*(J2-1)+3,3*(J2-1)+3)+SHIFTL(1)
      ENDDO

C     PRINT*,'UNSHIFTED HESS:'
C     WRITE(*,'(6F15.5)') ((HESS(J1,J2),J1=1,NOPT),J2=1,NOPT)

      DO J1=1,3*NATOMS
         DO J2=1,3*NATOMS
            HESS(J2,J1)=HESS(J2,J1)
     &                             +SHIFTL(1)*ETX(J2)*ETX(J1)
     &                             +SHIFTL(2)*ETY(J2)*ETY(J1)
     &                             +SHIFTL(3)*ETZ(J2)*ETZ(J1)
     &                             +SHIFTL(4)*ERX(J2)*ERX(J1)
     &                             +SHIFTL(5)*ERY(J2)*ERY(J1)
     &                             +SHIFTL(6)*ERZ(J2)*ERZ(J1)
         ENDDO
      ENDDO

C     PRINT*,'SHIFTED HESS:'
C     WRITE(*,'(6F15.5)') ((HESS(J1,J2),J1=1,NOPT),J2=1,NOPT)
C     PRINT*,'COORDINATES:'
C     WRITE(*,'(6F15.5)') (Q(J1),J1=1,NOPT)

      RETURN
      END


C     TEST FOR Z-ALIGNED DIPOLES, RETURNING TRUE IF FOUND AND FALSE OTHERWISE.
      LOGICAL FUNCTION ZALIGNTEST(NATOMS, X)
         USE KEY,ONLY : STOCKZTOL, STOCKMAXSPIN

         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)

         INTEGER J1

         DO J1 = NATOMS/2+1, NATOMS
            IF ( (1.0D0 - ABS(COS(X(J1*3-2)))) .LT. STOCKZTOL) THEN
               ZALIGNTEST = .TRUE.
               RETURN
            END IF
         END DO

         ZALIGNTEST = .FALSE.

      END FUNCTION ZALIGNTEST
