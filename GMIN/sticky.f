C   GMIN: A PROGRAM FOR FINDING GLOBAL MINIMA
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF GMIN.
C
C   GMIN IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   GMIN IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C
C  ENERGY AND GRADIENT FOR STICKY PATCHES MODEL.
C
      SUBROUTINE STICKY(X,V,ESTICKY,GTEST,SECT)
      USE COMMONS

      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, NAT2, J3, J4, J3MIN, J4MIN
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), X1, X2, Y1, Y2, Z1, Z2, 
     1                 ESTICKY, DUMMY, RALPHA12, RALPHA22, SA1, SA2,
     2                 M1, L1, N1, M2, L2, N2, ALPHA1, CA1, CA2, ALPHA2, S1, S2, C3A1, C3A2,
     3                 L12, M12, N12, L22, M22, N22, C3A1DOT1, C3A2DOT2,
     4                 GX1, GY1, GZ1, C2A2, C2A1, DOT1, DOT2, DIST, 
     6                 P1X, P1Y, P1Z, P2X, P2Y, P2Z, RSQ, R123, R1212, R126, 
     7           SIGMA22, R128, R1214, R12,
     8           OMG1,OMG2,CSOMG1,CSOMG2,PX,PY,PZ,CSOMGM1,CSOMGM2,
     9           FOMG1,FOMG2,
     +           OMG1B,OMG2B,CSOMG1B,CSOMG2B,DUMMYB


      SIGMA22=STICKYSIG**2
      NAT2=NATOMS/2
      ESTICKY=0.0D0
      DO J1=1,NATOMS/2
         VT(J1)=0.0D0
      ENDDO
      DO J1=1,3*NATOMS
         V(J1)=0.0D0
      ENDDO

C     DO J1=1,NRBSITES
C        PRINT '(A,I5,3G20.10)','PATCH ',J1, SITE(J1,1),SITE(J1,2),SITE(J1,3)
C     ENDDO

C  POTENTIAL ENERGY FIRST.
C
      DO J1=1,NAT2-1
         X1=X(3*(J1-1)+1)
         Y1=X(3*(J1-1)+2)
         Z1=X(3*(J1-1)+3)
         L1=X(3*(NAT2+J1-1)+1)
         M1=X(3*(NAT2+J1-1)+2)
         N1=X(3*(NAT2+J1-1)+3)
         L12=L1**2
         M12=M1**2
         N12=N1**2
         ALPHA1=SQRT(L12+M12+N12)
!        RALPHA12=1.0D0/MAX(ALPHA1**2,1.0D-20)
         CA1=COS(ALPHA1)
         SA1=SIN(ALPHA1)
!        C2A1=CA1
!	ANOTHER CHANGE HERE! THERE WAS A MISTAKE IN THE LIMIT CASE OF C3A2
         IF (ALPHA1.LT.0.0010D0) THEN
C           C3A1=-1/2+ALPHA1**2/24  ! BUG !! DJW
            C3A1=-0.5D0+ALPHA1**2/24
            S1=1.0D0-ALPHA1**2/6
C            C3A1=-1/2+ALPHA1**2/24-ALPHA1**4/720+ALPHA1**6/40320 ! THE 1/2 IS WRONG HERE AS WELL!
C            S1=1.0D0-ALPHA1**2/6+ALPHA1**4/120-ALPHA1**6/5040
         ELSE
            C3A1=(CA1-1.0D0)/ALPHA1**2
            S1=SIN(ALPHA1)/ALPHA1
         ENDIF
C        WRITE(*,'(A,6F15.5)') 'ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1=',ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1
         DO J2=J1+1,NAT2
            X2=X(3*(J2-1)+1)
            Y2=X(3*(J2-1)+2)
            Z2=X(3*(J2-1)+3)
            L2=X(3*(NAT2+J2-1)+1)
            M2=X(3*(NAT2+J2-1)+2)
            N2=X(3*(NAT2+J2-1)+3)
            L22=L2**2
            M22=M2**2
            N22=N2**2
            ALPHA2=SQRT(L22+M22+N22)
!           RALPHA22=1.0D0/MAX(ALPHA2**2,1.0D-20)
            CA2=COS(ALPHA2)
            SA2=SIN(ALPHA2)
!           C2A2=CA2
!	ANOTHER CHANGE HERE! THERE WAS A MISTAKE IN THE LIMIT CASE OF C3A2
            IF (ALPHA2.LT.0.00001D0) THEN
C              C3A2=-1/2+ALPHA2**2/24  ! BUG !!
               C3A2=-0.5D0+ALPHA2**2/24
               S2=1.0D0-ALPHA2**2/6
C               C3A2=-1/2+ALPHA2**2/24-ALPHA2**4/720+ALPHA2**6/40320 ! 1/2 IS A BUG HERE AS WELL
C               S2=1.0D0-ALPHA2**2/6+ALPHA2**4/120-ALPHA2**6/5040
            ELSE
               C3A2=(CA2-1.0D0)/ALPHA2**2
               S2=SIN(ALPHA2)/ALPHA2
            ENDIF
            RSQ=(X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2
            R12=SQRT(RSQ)
            R123=R12*RSQ
            R126=RSQ**3
            R1212=R126*R126
            R128=R126*RSQ
            R1214=R1212*RSQ
C
C  IDENTIFY CLOSEST PATCHES AND LOAD INTO P1X, P1Y, P1Z, P2X, P2Y, P2Z
C
CC
C            DMIN=1.0D100
C            DO J3=1,NRBSITES
C               P1X=SITE(J3,1)
C               P1Y=SITE(J3,2)
C               P1Z=SITE(J3,3)
C               DO J4=1,NRBSITES
C                  P2X=SITE(J4,1)
C                  P2Y=SITE(J4,2)
C                  P2Z=SITE(J4,3)
C                  DIST=
C     1  SQRT((C2A1*P1X-C3A1*L1*(L1*P1X+M1*P1Y+N1*P1Z)-C2A2*P2X+C3A2*L2*(L2*P2X + M2*P2Y + N2*P2Z)+
C     1      (N1*P1Y - M1*P1Z)*S1 - (N2*P2Y - M2*P2Z)*S2 + X1 - X2)**2 +
C     1   (C2A1*P1Y - C3A1*M1*(L1*P1X + M1*P1Y + N1*P1Z) - C2A2*P2Y+C3A2*M2*(L2*P2X+M2*P2Y+N2*P2Z)+
C     1      (-(N1*P1X) + L1*P1Z)*S1 - (-(N2*P2X) + L2*P2Z)*S2 + Y1-Y2)**2 +
C     1   (C2A1*P1Z - C3A1*N1*(L1*P1X + M1*P1Y + N1*P1Z) - C2A2*P2Z+C3A2*N2*(L2*P2X+M2*P2Y+N2*P2Z)+
C     1      (M1*P1X - L1*P1Y)*S1 - (M2*P2X - L2*P2Y)*S2 + Z1 - Z2)**2)
C                  IF (DIST.LT.DMIN) THEN
C                     J3MIN=J3
C                     J4MIN=J4
C                     DMIN=DIST
C                  ENDIF
CC                 PRINT '(A,2I5,2G20.10)','J3,J4,DIST,DMIN=',J3,J4,DIST,DMIN
C               ENDDO
C            ENDDO
            CSOMGM1=-1.0D100
            CSOMGM2=-1.0D100
            DO J3=1,NRBSITES
                  PX=SITE(J3,1)
                  PY=SITE(J3,2)
                  PZ=SITE(J3,3)
                  DOT1=L1*PX + M1*PY + N1*PZ
                  C3A1DOT1=C3A1*DOT1
                  DOT2=L2*PX + M2*PY + N2*PZ
                  C3A2DOT2=C3A2*DOT2
                  CSOMG1=((C3A1DOT1*L1 - CA1*PX - N1*PY*S1 +
     @                 M1*PZ*S1)*(X1 - X2) +
     @              (C3A1DOT1*M1 - CA1*PY + N1*PX*S1 - L1*PZ*S1)*
     @               (Y1 - Y2) + (C3A1DOT1*N1 - CA1*PZ - M1*PX*S1 +
     @                 L1*PY*S1)*(Z1 - Z2))/R12
                  CSOMG2=((-(C3A2DOT2*L2) + CA2*PX + N2*PY*S2 - M2*PZ*S2)*
     @                (X1 - X2) + (-(C3A2DOT2*M2) + CA2*PY - N2*PX*S2 +
     @                  L2*PZ*S2)*(Y1 - Y2) +
     @               (-(C3A2DOT2*N2) + CA2*PZ + M2*PX*S2 - L2*PY*S2)*
     @                (Z1 - Z2))/R12
                  IF (CSOMG1.GT.CSOMGM1) THEN
                     J3MIN=J3
                     CSOMGM1=CSOMG1
                  ENDIF
                  IF (CSOMG2.GT.CSOMGM2) THEN
                     J4MIN=J3
                     CSOMGM2=CSOMG2
                  ENDIF
            ENDDO
C           PRINT*,'J3MIN,J4MIN=',J3MIN,J4MIN
           P1X=SITE(J3MIN,1)
           P1Y=SITE(J3MIN,2)
           P1Z=SITE(J3MIN,3)
           P2X=SITE(J4MIN,1)
           P2Y=SITE(J4MIN,2)
            P2Z=SITE(J4MIN,3)

            DOT1=L1*P1X + M1*P1Y + N1*P1Z 
            C3A1DOT1=C3A1*DOT1
            DOT2=L2*P2X + M2*P2Y + N2*P2Z 
            C3A2DOT2=C3A2*DOT2
C            OJO, ESTO LO HE COMENTADO YO
C            P1X=SITE(1,1)
C            P1Y=SITE(1,2)
C            P1Z=SITE(1,3)
C            P2X=SITE(1,1)
C            P2Y=SITE(1,2)
C            P2Z=SITE(1,3)
C           PRINT*,'COORDINATES OF MOLECULE PAIR:'
C           WRITE(*,'(3F20.10)') X1,Y1,Z1
C           WRITE(*,'(3F20.10)') L1,M1,N1
C           WRITE(*,'(3F20.10)') X2,Y2,Z2
C           WRITE(*,'(3F20.10)') L2,M2,N2
            CSOMG1=((C3A1DOT1*L1 - CA1*P1X - N1*P1Y*S1 + 
     @           M1*P1Z*S1)*(X1 - X2) +
     @        (C3A1DOT1*M1 - CA1*P1Y + N1*P1X*S1 - L1*P1Z*S1)*
     @         (Y1 - Y2) + (C3A1DOT1*N1 - CA1*P1Z - M1*P1X*S1 +
     @           L1*P1Y*S1)*(Z1 - Z2))/R12
            CSOMG1B= (-P1X*(X1-X2)-P1Y*(Y1-Y2)-P1Z*(Z1-Z2))/R12
            OMG1B=ACOS(CSOMG1B)
            OMG1=ACOS(CSOMG1)
            IF(OMG1.LT.0.001D00) THEN
                FOMG1=1+OMG1**2/6
            ELSE
                FOMG1=OMG1/SIN(OMG1)
            ENDIF
            CSOMG2=((-C3A2DOT2*L2 + CA2*P2X + N2*P2Y*S2 - M2*P2Z*S2)*
     @         (X1 - X2) + (-(C3A2DOT2*M2) + CA2*P2Y - N2*P2X*S2 +
     @           L2*P2Z*S2)*(Y1 - Y2) +
     @        (-(C3A2DOT2*N2) + CA2*P2Z + M2*P2X*S2 - L2*P2Y*S2)*
     @         (Z1 - Z2))/R12
            CSOMG2B= (P2X*(X1-X2)+P2Y*(Y1-Y2)+P2Z*(Z1-Z2))/R12
            OMG2B=ACOS(CSOMG2B)
            OMG2=ACOS(CSOMG2)
            IF(OMG2.LT.0.001D00) THEN
                FOMG2=1+OMG2**2/6
            ELSE
                FOMG2=OMG2/SIN(OMG2)
            ENDIF
            DUMMY= (-4*(R1212 - R126)*EXP(-(OMG1**2+OMG2**2)/(2.*SIGMA22)))/(R1212*R126)
            DUMMYB= (-4*(R1212 - R126)*EXP(-(OMG1B**2+OMG2B**2)/(2.*SIGMA22)))/(R1212*R126)
C	PRINT*,'DUMMY=',DUMMY,'DUMMYB=',DUMMYB
            IF (GTEST.OR.SECT) THEN
               GX1=
     @ (4*(-12/R1214 + 6/R128)*(X1 - X2) - 
     @ (4*(1/R1212 - 1/R126)*( (-((-C3A2DOT2*L2 + CA2*P2X + 
     @              (N2*P2Y - M2*P2Z)*S2)/R12) + 
     @         ((X1 - X2)*CSOMG2)/RSQ)*FOMG2 +
     @    ( -((C3A1DOT1*L1 - CA1*P1X - (N1*P1Y - M1*P1Z)*S1)/R12) + 
     @         ((X1 - X2)*CSOMG1)/RSQ)*FOMG1))/SIGMA22)*
     @    EXP((-OMG1**2-OMG2**2)/(2.*SIGMA22))

C	GX1=(4*(-12/R1214 + 6/R128)*(X1 - X2) - 
C     @ (4*(1/R1212 - 1/R126)*(( -P2X/R12+ ((X1 - X2)*CSOMG2)/RSQ)*OMG2B/SIN(OMG2B) +
C     @   (P1X/R12+((X1-X2)*CSOMG1B)/RSQ)*OMG1B/SIN(OMG1B)))/SIGMA22)*
C     @    EXP((-OMG1B**2-OMG2B**2)/(2.*SIGMA22))


               V(3*(J1-1)+1)=V(3*(J1-1)+1)+GX1

               GY1=
     @  (4*(-12/R1214 + 6/R128)*(Y1 - Y2) - 
     @  (4*(1/R1212 - 1/R126)*( (-((-C3A2DOT2*M2 + CA2*P2Y + 
     @              (-(N2*P2X) + L2*P2Z)*S2)/R12) + 
     @         ((Y1 - Y2)*CSOMG2)/RSQ)*FOMG2+
     @    ( -((C3A1DOT1*M1 - CA1*P1Y - (-(N1*P1X) + L1*P1Z)*S1)/R12) + 
     @         ((Y1 - Y2)*CSOMG1)/RSQ)*FOMG1))/SIGMA22)*
     @   EXP((-OMG1**2-OMG2**2)/(2.*SIGMA22))

               V(3*(J1-1)+2)=V(3*(J1-1)+2)+GY1

         GZ1=
     @   (4*(-12/R1214 + 6/R128)*(Z1 - Z2) - 
     @   (4*(1/R1212 - 1/R126)*( (-((-C3A2DOT2*N2 + CA2*P2Z + 
     @              (M2*P2X - L2*P2Y)*S2)/R12) + 
     @         ((Z1-Z2)*CSOMG2)/RSQ)*FOMG2+
     @    ( -((C3A1DOT1*N1 - CA1*P1Z - (M1*P1X - L1*P1Y)*S1)/R12) + 
     @         ((Z1 - Z2)*CSOMG1)/RSQ)*FOMG1))/SIGMA22)*
     @  EXP((-OMG1**2-OMG2**2)/(2.*SIGMA22))
               V(3*(J1-1)+3)=V(3*(J1-1)+3)+GZ1

               V(3*(NAT2+J1-1)+1)=V(3*(NAT2+J1-1)+1)+
     @   (4*(1/R1212 - 1/R126)*
     @   ((-C3A1DOT1 - (2*(1 - CA1)*DOT1*L1**2)/ALPHA1**4 - C3A1*L1*P1X + 
     @    (CA1*L1*(N1*P1Y - M1*P1Z))/ALPHA1**2 - L1*P1X*S1 + 
     @    (DOT1*L1**2*SA1)/ALPHA1**3 - 
     @    (L1*(N1*P1Y - M1*P1Z)*SA1)/(ALPHA1**3))*(-X1 + X2) + 
     @ (P1Z*S1 + L1*(-(P1Y*S1) + 
     @       (-(N1*P1X) + L1*P1Z)*(CA1/ALPHA1**2 - SA1/ALPHA1**3)) 
     @  + M1*(-(C3A1*P1X) + DOT1*L1*((-2*(1 - CA1))/ALPHA1**4 + 
     @          SA1/(ALPHA1**3))))*(-Y1 + Y2) + 
     @ (-(P1Y*S1) + L1*((CA1*(M1*P1X - L1*P1Y))/ALPHA1**2 - P1Z*S1) - 
     @    (L1*(M1*P1X - L1*P1Y)*SA1)/(ALPHA1**3) + 
     @    N1*((-2*(1 - CA1)*DOT1*L1)/ALPHA1**4 - C3A1*P1X + 
     @       (DOT1*L1*SA1)/(ALPHA1**3)))*(-Z1 + Z2))*FOMG1*
     @   EXP((-OMG1**2-OMG2**2)/(2.*SIGMA22)))/
     @   (R12*SIGMA22)

               V(3*(NAT2+J1-1)+2)=V(3*(NAT2+J1-1)+2)+
     @   (4*(1/R1212 - 1/R126)*
     @   ((-P1Z*S1 + M1*((CA1*(N1*P1Y - M1*P1Z))/ALPHA1**2 - P1X*S1) - 
     @    (M1*(N1*P1Y - M1*P1Z)*SA1)/ALPHA1**3 + 
     @    L1*((-2*(1 - CA1)*DOT1*M1)/ALPHA1**4 - C3A1*P1Y + 
     @       (DOT1*M1*SA1)/ALPHA1**3))*(-X1 + X2) + 
     @ (-C3A1DOT1 - (2*(1 - CA1)*DOT1*M1**2)/ALPHA1**4 - C3A1*M1*P1Y + 
     @    (CA1*M1*(-N1*P1X + L1*P1Z))/ALPHA1**2 - M1*P1Y*S1 + 
     @    (DOT1*M1**2*SA1)/ALPHA1**3 - 
     @    (M1*(-N1*P1X + L1*P1Z)*SA1)/ALPHA1**3)*(-Y1 + Y2) + 
     @ (P1X*S1 + M1*(-(P1Z*S1) + 
     @       (M1*P1X - L1*P1Y)*(CA1/ALPHA1**2 - SA1/ALPHA1**3)) + 
     @    N1*(-(C3A1*P1Y) + DOT1*M1*
     @        ((-2*(1 - CA1))/ALPHA1**4 + SA1/ALPHA1**3)))*(-Z1 + Z2)
     @  )*FOMG1*
     @   EXP((-OMG1**2-OMG2**2)/(2.*SIGMA22)))/
     @  (R12*SIGMA22)

               V(3*(NAT2+J1-1)+3)=V(3*(NAT2+J1-1)+3)+
     @ (4*(1/R1212 - 1/R126)* ((P1Y*S1 + N1*(-(P1X*S1) + 
     @       (N1*P1Y - M1*P1Z)*(CA1/ALPHA1**2 - SA1/ALPHA1**3)) + 
     @    L1*(-(C3A1*P1Z) + DOT1*N1*
     @        ((-2*(1 - CA1))/ALPHA1**4 + SA1/ALPHA1**3)))*
     @  (-X1 + X2) + (-(P1X*S1) + 
     @    N1*((CA1*(-(N1*P1X) + L1*P1Z))/ALPHA1**2 - P1Y*S1) - 
     @    (N1*(-(N1*P1X) + L1*P1Z)*SA1)/ALPHA1**3 + 
     @    M1*((-2*(1 - CA1)*DOT1*N1)/ALPHA1**4 - C3A1*P1Z + 
     @       (DOT1*N1*SA1)/ALPHA1**3))*(-Y1 + Y2) + 
     @ (-C3A1DOT1 - (2*(1 - CA1)*DOT1*N1**2)/ALPHA1**4 + 
     @    (CA1*N1*(M1*P1X - L1*P1Y))/ALPHA1**2 - C3A1*N1*P1Z - N1*P1Z*S1 + 
     @    (DOT1*N1**2*SA1)/ALPHA1**3 - 
     @    (N1*(M1*P1X - L1*P1Y)*SA1)/ALPHA1**3)*(-Z1 + Z2))*
     @  FOMG1*EXP((-OMG1**2-OMG2**2)/(2.*SIGMA22)))/
     @  (R12*SIGMA22)

               V(3*(J2-1)+1)=V(3*(J2-1)+1)-GX1

               V(3*(J2-1)+2)=V(3*(J2-1)+2)-GY1

               V(3*(J2-1)+3)=V(3*(J2-1)+3)-GZ1

               V(3*(NAT2+J2-1)+1)=V(3*(NAT2+J2-1)+1)+
     @  (4*(1/R1212 - 1/R126)*
     @  ((-C3A2DOT2 - (2*(1 - CA2)*DOT2*L2**2)/ALPHA2**4 - C3A2*L2*P2X + 
     @    (CA2*L2*(N2*P2Y - M2*P2Z))/ALPHA2**2 - L2*P2X*S2 + (DOT2*L2**2*SA2)/ALPHA2**3 - 
     @    (L2*(N2*P2Y - M2*P2Z)*SA2)/ALPHA2**3)*(X1 - X2) + (P2Z*S2 + L2*(-(P2Y*S2) + 
     @       (-(N2*P2X) + L2*P2Z)*(CA2/ALPHA2**2 - SA2/ALPHA2**3)) 
     @    + M2*(-(C3A2*P2X) + DOT2*L2*((-2*(1 - CA2))/ALPHA2**4 + 
     @          SA2/ALPHA2**3)))*(Y1 - Y2) + 
     @ (-(P2Y*S2) + L2*((CA2*(M2*P2X - L2*P2Y))/ALPHA2**2 - P2Z*S2) - 
     @    (L2*(M2*P2X - L2*P2Y)*SA2)/ALPHA2**3 + 
     @    N2*((-2*(1 - CA2)*DOT2*L2)/ALPHA2**4 - C3A2*P2X + 
     @       (DOT2*L2*SA2)/ALPHA2**3))*(Z1 - Z2))*FOMG2*
     @    EXP((-OMG1**2-OMG2**2)/(2.*SIGMA22)))/(R12*SIGMA22)

               V(3*(NAT2+J2-1)+2)=V(3*(NAT2+J2-1)+2)+
     @     (4*(1/R1212 - 1/R126)* ((-(P2Z*S2) + M2*((CA2*(N2*P2Y - M2*P2Z))/ALPHA2**2 - P2X*S2) - 
     @    (M2*(N2*P2Y - M2*P2Z)*SA2)/ALPHA2**3 + 
     @    L2*((-2*(1 - CA2)*DOT2*M2)/ALPHA2**4 - C3A2*P2Y + 
     @       (DOT2*M2*SA2)/ALPHA2**3))*(X1 - X2) + 
     @ (-C3A2DOT2 - (2*(1 - CA2)*DOT2*M2**2)/ALPHA2**4 - C3A2*M2*P2Y + 
     @    (CA2*M2*(-(N2*P2X) + L2*P2Z))/ALPHA2**2 - M2*P2Y*S2 + 
     @    (DOT2*M2**2*SA2)/ALPHA2**3 - 
     @    (M2*(-(N2*P2X) + L2*P2Z)*SA2)/ALPHA2**3)*(Y1 - Y2) + 
     @ (P2X*S2 + M2*(-(P2Z*S2) + 
     @       (M2*P2X - L2*P2Y)*(CA2/ALPHA2**2 - SA2/ALPHA2**3)) + 
     @    N2*(-(C3A2*P2Y) + DOT2*M2*
     @        ((-2*(1 - CA2))/ALPHA2**4 + SA2/ALPHA2**3)))*(Z1 - Z2))*
     @  FOMG2*EXP((-OMG1**2-OMG2**2)/(2.*SIGMA22)))/ (R12*SIGMA22)

               V(3*(NAT2+J2-1)+3)=V(3*(NAT2+J2-1)+3)+
     @  (4*(1/R1212 - 1/R126)* ((P2Y*S2 + N2*(-(P2X*S2) + 
     @       (N2*P2Y - M2*P2Z)*(CA2/ALPHA2**2 - SA2/ALPHA2**3)) + L2*(-(C3A2*P2Z) + DOT2*N2*
     @        ((-2*(1 - CA2))/ALPHA2**4 + SA2/ALPHA2**3)))*(X1 - X2) 
     @   + (-(P2X*S2) + N2*((CA2*(-(N2*P2X) + L2*P2Z))/ALPHA2**2 - P2Y*S2) - 
     @    (N2*(-(N2*P2X) + L2*P2Z)*SA2)/ALPHA2**3 + 
     @    M2*((-2*(1 - CA2)*DOT2*N2)/ALPHA2**4 - C3A2*P2Z + 
     @       (DOT2*N2*SA2)/ALPHA2**3))*(Y1 - Y2) + 
     @ (-C3A2DOT2 - (2*(1 - CA2)*DOT2*N2**2)/ALPHA2**4 + 
     @    (CA2*N2*(M2*P2X - L2*P2Y))/ALPHA2**2 - C3A2*N2*P2Z - N2*P2Z*S2 + 
     @    (DOT2*N2**2*SA2)/ALPHA2**3 - 
     @    (N2*(M2*P2X - L2*P2Y)*SA2)/ALPHA2**3)*(Z1 - Z2))*FOMG2*
     @   EXP((-OMG1**2-OMG2**2)/(2.*SIGMA22)))/ (R12*SIGMA22)

            ENDIF
   
            ESTICKY=ESTICKY+DUMMY
            VT(J1)=VT(J1)+DUMMY
            VT(J2)=VT(J2)+DUMMY

         ENDDO
      ENDDO

C      WRITE(*,'(A,G20.10)') 'ENERGY=',ESTICKY
C      PRINT*,'COORDS:'
C      WRITE(*,'(I6,G20.10)') (J1,X(J1),J1=1,3*NATOMS)
C      PRINT*,'GRADIENT:'
C      WRITE(*,'(I6,G20.10)') (J1,V(J1),J1=1,3*NATOMS)

      RETURN
      END
