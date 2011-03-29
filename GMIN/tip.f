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
C  ENERGY AND GRADIENT FOR RIGID BODY TIP4P USING NEW RIGID BODY
C  DERIVATIVE FUNCTIONS ETC.
C
      SUBROUTINE TIP4P(X,V,ETIP,GTEST,SECT)
      USE COMMONS

      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, K1, K2, NAT2, NSITE
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), X1, X2, Y1, Y2, Z1, Z2, FLJ, DFLJ, FC, DFC,
     1                 ETIP, DUMMY, RALPHA12, RALPHA22, RDIST, 
     2                 M1, L1, N1, M2, L2, N2, ALPHA1, CA1, CA2, ALPHA2, S1, S2, C3A1, C3A2,
     3                 L12, M12, N12, L22, M22, N22, C6, C12,
     4                 GX1, GY1, GZ1, GL1, GM1, GN1, GL2, GM2, GN2, C2A2, C2A1,
     5                 D1, D2, D3, D4, D5, D6, DA, DB, DC, DIST, DUMMY2,
     6                 P1X, P1Y, P1Z, P2X, P2Y, P2Z, XDUMM, Q1, Q2,
     7           D1S, D2S, D3S, D4S, D5S, D6S, DAS, DBS, DCS, CHARGE(5), HTOKJ
      PARAMETER (HTOKJ=1389.354848D0) ! CONVERSION FACTOR FOR COULOMB ENERGY AND GRADIENTS
C
C  STATEMENT FUNCTIONS.
C  SITE-SITE ENERGY TERMS - COULOMBIC AND LJ. XDUMM IS THE SITE-SITE DISTANCE.
C
      FLJ(XDUMM)=(C12/XDUMM**6-C6)/XDUMM**6
      DFLJ(XDUMM)=-6.0D0*(-C6+2.0D0*C12/XDUMM**6)/XDUMM**7
      FC(Q1,Q2,XDUMM)=Q1*Q2/XDUMM
      DFC(Q1,Q2,XDUMM)=-Q1*Q2/XDUMM**2

C     PRINT*,'COORDS IN TIP:'
C     WRITE(*,'(3F20.10)') (X(J1),J1=1,3*NATOMS)
C
C DISTINGUISH TIP1P, TIP2P, TIP3P, TIP4P
C
      SITE(1,1)=0.0D0
      SITE(1,2)=0.0D0
      SITE(1,3)=0.0D0
C     SITE(1,3)=0.065098D0
      SITE(2,1)=0.756950327D0
      SITE(2,2)=0.0D0
      SITE(2,3)=-0.5858822760D0
C     SITE(2,3)=-0.5207842760D0
      SITE(3,1)=-0.756950327D0
      SITE(3,2)=0.0D0
      SITE(3,3)=-0.5858822760D0
C     SITE(3,3)=-0.5207842760D0
      SITE(4,1)=0.0D0
      SITE(4,2)=0.0D0
      CHARGE(1)=0.0D0
      NSITE=4

      IF (TIPID.EQ.1) THEN
         C6=2510.4D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
         C12=2426720.0D0
         SITE(4,3)=0.0D0
C        SITE(4,3)=0.065098D0
         CHARGE(2)=0.4D0
         CHARGE(3)=0.4D0
         CHARGE(4)=-0.8D0
      ELSE IF (TIPID.EQ.2) THEN
         C6=2510.4D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
         C12=2907880.0D0
         SITE(4,3)=-0.15D0
C        SITE(4,3)=-0.084902D0
         CHARGE(2)=0.535D0
         CHARGE(3)=0.535D0
         CHARGE(4)=-1.07D0
      ELSE IF (TIPID.EQ.3) THEN
         C6=2489.48D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
         C12=2435088.0D0
         SITE(4,3)=0.0D0
C        SITE(4,3)=0.065098D0
         CHARGE(2)=0.417D0
         CHARGE(3)=0.417D0
         CHARGE(4)=-0.834D0
      ELSE IF (TIPID.EQ.4) THEN
         C6=2552.24D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
         C12=2510.4D3
         SITE(4,3)=-0.15D0
C        SITE(4,3)=-0.084902D0
         CHARGE(2)=0.52D0
         CHARGE(3)=0.52D0
         CHARGE(4)=-1.04D0
      ELSE IF (TIPID.EQ.5) THEN
! C6  = 590.3472412D0 KCAL/MOL ANGSTROM**6
! C12 = 544546.6644D0 KCAL/MOL ANGSTROM**12 
         C6=2470.012857D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
         C12=2278383.244D0
         NSITE=5
         SITE(4,1)=0.0D0; SITE(4,2)= 0.571543301D0; SITE(4,3)=0.404151276D0
         SITE(5,1)=0.0D0; SITE(5,2)=-0.571543301D0; SITE(5,3)=0.404151276D0
         CHARGE(2)=0.241D0
         CHARGE(3)=0.241D0
         CHARGE(4)=-0.241D0
         CHARGE(5)=-0.241D0
      ENDIF

      NAT2=NATOMS/2
      ETIP=0.0D0
      DO J1=1,NATOMS/2
         VT(J1)=0.0D0
      ENDDO
      DO J1=1,3*NATOMS
         V(J1)=0.0D0
      ENDDO

C     PRINT*,'SITE 1', SITE(1,1),SITE(1,2),SITE(1,3)
C     PRINT*,'SITE 2', SITE(2,1),SITE(2,2),SITE(2,3)
C     PRINT*,'SITE 3', SITE(3,1),SITE(3,2),SITE(3,3)
C     PRINT*,'SITE 4', SITE(4,1),SITE(4,2),SITE(4,3)

C
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
         RALPHA12=1.0D0/MAX(ALPHA1**2,1.0D-20)
         CA1=COS(ALPHA1)
         C2A1=CA1
         IF (ALPHA1.LT.0.0001D0) THEN
C           C3A1=-ALPHA1/2+ALPHA1**3/24 ! BUG SPOTTED BY TIM!
            C3A1=-0.5D0+ALPHA1**2/24.0D0
            S1=1.0D0-ALPHA1**2/6
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
            RALPHA22=1.0D0/MAX(ALPHA2**2,1.0D-20)
            CA2=COS(ALPHA2)
            C2A2=CA2
            IF (ALPHA2.LT.0.00001D0) THEN
C              C3A2=-ALPHA2/2+ALPHA2**3/24 ! BUG SPOTTED BY TIM!
               C3A2=-0.5D0+ALPHA2**2/24.0D0
               S2=1.0D0-ALPHA2**2/6
            ELSE
               C3A2=(CA2-1.0D0)/ALPHA2**2
               S2=SIN(ALPHA2)/ALPHA2
            ENDIF
C
C  O-O LJ CONTRIBUTION.
C
!
!  SINCE SITE 1 (O) IS AT THE ORIGIN WE DON'T NEED TO ROTATE IT, SO THE LJ TERM
!  CAN BE DONE MUCH FASTER.
!
            DIST=SQRT((X1-X(3*(J2-1)+1))**2+(Y1-X(3*(J2-1)+2))**2+(Z1-X(3*(J2-1)+3))**2)
            DUMMY=FLJ(DIST)
            DUMMY2=DFLJ(DIST)
            GX1=DUMMY2*(X1-X(3*(J2-1)+1))/DIST
            GY1=DUMMY2*(Y1-X(3*(J2-1)+2))/DIST
            GZ1=DUMMY2*(Z1-X(3*(J2-1)+3))/DIST
            GL1=0.0D0
            GM1=0.0D0
            GN1=0.0D0
            GL2=0.0D0
            GM2=0.0D0
            GN2=0.0D0

!             P1X=SITE(1,1)
!             P1Y=SITE(1,2)
!             P1Z=SITE(1,3)
!             P2X=SITE(1,1)
!             P2Y=SITE(1,2)
!             P2Z=SITE(1,3)
!             DIST= 
!      1  SQRT((C2A1*P1X-C3A1*L1*(L1*P1X+M1*P1Y+N1*P1Z)-C2A2*P2X+C3A2*L2*(L2*P2X + M2*P2Y + N2*P2Z)+ 
!      1      (N1*P1Y - M1*P1Z)*S1 - (N2*P2Y - M2*P2Z)*S2 + X1 - X2)**2 + 
!      1   (C2A1*P1Y - C3A1*M1*(L1*P1X + M1*P1Y + N1*P1Z) - C2A2*P2Y+C3A2*M2*(L2*P2X+M2*P2Y+N2*P2Z)+ 
!      1      (-(N1*P1X) + L1*P1Z)*S1 - (-(N2*P2X) + L2*P2Z)*S2 + Y1-Y2)**2 + 
!      1   (C2A1*P1Z - C3A1*N1*(L1*P1X + M1*P1Y + N1*P1Z) - C2A2*P2Z+C3A2*N2*(L2*P2X+M2*P2Y+N2*P2Z)+ 
!      1      (M1*P1X - L1*P1Y)*S1 - (M2*P2X - L2*P2Y)*S2 + Z1 - Z2)**2)
! C           PRINT*,'COORDINATES OF MOLECULE PAIR:'
! C           WRITE(*,'(3F20.10)') X1,Y1,Z1
! C           WRITE(*,'(3F20.10)') L1,M1,N1
! C           WRITE(*,'(3F20.10)') X2,Y2,Z2
! C           WRITE(*,'(3F20.10)') L2,M2,N2
! C           PRINT*,'DISTANCE=',DIST
!             DUMMY=FLJ(DIST)
!             IF (GTEST.OR.SECT) THEN
!                DUMMY2=DFLJ(DIST)
!                RDIST=1.0D0/DIST
!                D1S=D1(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
!      1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
!                D2S=D2(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
!      1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
!                D3S=D3(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
!      1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
!                D4S=D4(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
!      1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
!                D5S=D5(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
!      1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
!                D6S=D6(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
!      1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
!                DAS=DA(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
!      1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
!                DBS=DB(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
!      1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
!                DCS=DC(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
!      1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
!                GX1=DUMMY2*D1S
!                GY1=DUMMY2*D2S
!                GZ1=DUMMY2*D3S
!                GL1=DUMMY2*D4S
!                GM1=DUMMY2*D5S
!                GN1=DUMMY2*D6S
!                GL2=DUMMY2*DAS
!                GM2=DUMMY2*DBS
!                GN2=DUMMY2*DCS
!             ENDIF
C
C  SUM OVER CHARGED SITES. THIS COULD ALSO BE FASTER IF WE DIDN'T USE THE GENERAL
C  RIGID BODY FORMULATION.
C
            DO K1=2,NSITE
               P1X=SITE(K1,1)
               P1Y=SITE(K1,2)
               P1Z=SITE(K1,3)
               DO K2=2,NSITE
                  P2X=SITE(K2,1)
                  P2Y=SITE(K2,2)
                  P2Z=SITE(K2,3)
                  DIST=
     1  SQRT((C2A1*P1X-C3A1*L1*(L1*P1X+M1*P1Y+N1*P1Z)-C2A2*P2X+C3A2*L2*(L2*P2X + M2*P2Y + N2*P2Z)+ 
     1      (N1*P1Y - M1*P1Z)*S1 - (N2*P2Y - M2*P2Z)*S2 + X1 - X2)**2 + 
     1   (C2A1*P1Y - C3A1*M1*(L1*P1X + M1*P1Y + N1*P1Z) - C2A2*P2Y+C3A2*M2*(L2*P2X+M2*P2Y+N2*P2Z)+ 
     1      (-(N1*P1X) + L1*P1Z)*S1 - (-(N2*P2X) + L2*P2Z)*S2 + Y1-Y2)**2 + 
     1   (C2A1*P1Z - C3A1*N1*(L1*P1X + M1*P1Y + N1*P1Z) - C2A2*P2Z+C3A2*N2*(L2*P2X+M2*P2Y+N2*P2Z)+ 
     1      (M1*P1X - L1*P1Y)*S1 - (M2*P2X - L2*P2Y)*S2 + Z1 - Z2)**2)
                  DUMMY=DUMMY+HTOKJ*FC(CHARGE(K1),CHARGE(K2),DIST)
                  DUMMY2=HTOKJ*DFC(CHARGE(K1),CHARGE(K2),DIST)
                  IF (GTEST) THEN
                     RDIST=1.0D0/DIST
                     GX1=GX1+DUMMY2*D1(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GY1=GY1+DUMMY2*D2(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GZ1=GZ1+DUMMY2*D3(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GL1=GL1+DUMMY2*D4(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GM1=GM1+DUMMY2*D5(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GN1=GN1+DUMMY2*D6(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GL2=GL2+DUMMY2*DA(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GM2=GM2+DUMMY2*DB(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GN2=GN2+DUMMY2*DC(P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,RDIST,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                  ENDIF
               ENDDO
            ENDDO
   
            ETIP=ETIP+DUMMY
            VT(J1)=VT(J1)+DUMMY
            VT(J2)=VT(J2)+DUMMY

            IF (GTEST) THEN
               V(3*(J1-1)+1)=V(3*(J1-1)+1)+GX1
               V(3*(J1-1)+2)=V(3*(J1-1)+2)+GY1
               V(3*(J1-1)+3)=V(3*(J1-1)+3)+GZ1
               V(3*(NAT2+J1-1)+1)=V(3*(NAT2+J1-1)+1)+GL1
               V(3*(NAT2+J1-1)+2)=V(3*(NAT2+J1-1)+2)+GM1
               V(3*(NAT2+J1-1)+3)=V(3*(NAT2+J1-1)+3)+GN1
               V(3*(J2-1)+1)=V(3*(J2-1)+1)-GX1
               V(3*(J2-1)+2)=V(3*(J2-1)+2)-GY1
               V(3*(J2-1)+3)=V(3*(J2-1)+3)-GZ1
               V(3*(NAT2+J2-1)+1)=V(3*(NAT2+J2-1)+1)+GL2
               V(3*(NAT2+J2-1)+2)=V(3*(NAT2+J2-1)+2)+GM2
               V(3*(NAT2+J2-1)+3)=V(3*(NAT2+J2-1)+3)+GN2
            ENDIF
         ENDDO
      ENDDO

C     WRITE(*,'(A,G20.10)') 'ENERGY=',ETIP
C     PRINT*,'COORDS:'
C     WRITE(*,'(I6,G20.10)') (J1,X(J1),J1=1,3*NATOMS)
C     PRINT*,'GRADIENT:'
C     WRITE(*,'(I6,G20.10)') (J1,V(J1),J1=1,3*NATOMS)

      RETURN
      END
