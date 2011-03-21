C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C  Energy and gradient for LJ plus quadrupole consisting of four charged sites.
C  Site one is for LJ, sites 2 to 5 are charges. Use C6=C12=1.0D0 Use C6=C12=1.0D0.
C
      SUBROUTINE QUAD(X,V,EQUAD,GTEST,SECT)
      USE commons

      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, K1, K2, NAT2, NSITE
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), X1, X2, Y1, Y2, Z1, Z2, FLJ, DFLJ, FC, DFC,
     1                 EQUAD, DUMMY, RALPHA12, RALPHA22, RDIST, 
     2                 M1, L1, N1, M2, L2, N2, ALPHA1, CA1, CA2, ALPHA2, S1, S2, C3A1, C3A2,
     3                 L12, M12, N12, L22, M22, N22, C6, C12,
     4                 GX1, GY1, GZ1, GL1, GM1, GN1, GL2, GM2, GN2, C2A2, C2A1,
     5                 D1, D2, D3, D4, D5, D6, DA, DB, DC, DIST, DUMMY2,
     6                 P1X, P1Y, P1Z, P2X, P2Y, P2Z, XDUMM, Q1, Q2,
     7           D1S, D2S, D3S, D4S, D5S, D6S, DAS, DBS, DCS, CHARGE(5), HTOKJ
!     PARAMETER (HTOKJ=1389.354848D0) ! conversion factor for coulomb energy and gradients
      PARAMETER (HTOKJ=1.0D00) ! conversion factor for coulomb energy and gradients
      LOGICAL :: CALLED=.FALSE.
      SAVE CHARGE
C
C  Statement functions.
C  Site-site energy terms - Coulombic and LJ. XDUMM is the site-site distance.
C
      FLJ(XDUMM)=(1.0D0/XDUMM**6-1.0D0)/XDUMM**6
      DFLJ(XDUMM)=-6.0D0*(-1.0D0+2.0D0/XDUMM**6)/XDUMM**7
      FC(Q1,Q2,XDUMM)=Q1*Q2/XDUMM
      DFC(Q1,Q2,XDUMM)=-Q1*Q2/XDUMM**2

      NSITE=5
      IF (.NOT.CALLED) THEN
         OPEN(UNIT=1,FILE='quad.data',STATUS='OLD')
         CALLED=.TRUE.
         DO J1=1,NSITE
            READ(1,*) CHARGE(J1),SITE(J1,1),SITE(J1,2),SITE(J1,3)
         ENDDO
         CLOSE(1)
         PRINT '(A,3G20.10)','quad> LJ site is at: ',SITE(1,1),SITE(2,2),SITE(3,3)
         PRINT '(A)','quad> The four charges and their positions are:'
         PRINT '(F15.5,3G20.10)', (CHARGE(J1),SITE(J1,1),SITE(J1,2),SITE(J1,3),J1=2,5)
      ENDIF

      NAT2=NATOMS/2
      EQUAD=0.0D0
      DO J1=1,NATOMS/2
         VT(J1)=0.0D0
      ENDDO
      DO J1=1,3*NATOMS
         V(J1)=0.0D0
      ENDDO
C
C  Potential energy first.
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
C           C3A1=-ALPHA1/2+ALPHA1**3/24 ! bug spotted by Tim!
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
               C3A2=-0.5D0+ALPHA2**2/24.0D0
               S2=1.0D0-ALPHA2**2/6
            ELSE
               C3A2=(CA2-1.0D0)/ALPHA2**2
               S2=SIN(ALPHA2)/ALPHA2
            ENDIF
C
C  LJ contribution.
C
            P1X=SITE(1,1)
            P1Y=SITE(1,2)
            P1Z=SITE(1,3)
            P2X=SITE(1,1)
            P2Y=SITE(1,2)
            P2Z=SITE(1,3)
            DIST= 
     1  Sqrt((c2a1*p1x-c3a1*l1*(l1*p1x+m1*p1y+n1*p1z)-c2a2*p2x+c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)+ 
     1      (n1*p1y - m1*p1z)*s1 - (n2*p2y - m2*p2z)*s2 + x1 - x2)**2 + 
     1   (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y+c3a2*m2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (-(n1*p1x) + l1*p1z)*s1 - (-(n2*p2x) + l2*p2z)*s2 + y1-y2)**2 + 
     1   (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z+c3a2*n2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (m1*p1x - l1*p1y)*s1 - (m2*p2x - l2*p2y)*s2 + z1 - z2)**2)
C           PRINT*,'coordinates of molecule pair:'
C           WRITE(*,'(3F20.10)') X1,Y1,Z1
C           WRITE(*,'(3F20.10)') L1,M1,N1
C           WRITE(*,'(3F20.10)') X2,Y2,Z2
C           WRITE(*,'(3F20.10)') L2,M2,N2
C           PRINT*,'distance=',DIST
            DUMMY=FLJ(DIST)
            IF (GTEST.OR.SECT) THEN
               DUMMY2=DFLJ(DIST)
               RDIST=1.0D0/DIST
               D1S=D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D2S=D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D3S=D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D4S=D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D5S=D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D6S=D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               DAS=DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               DBS=DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               DCS=DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1             C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               GX1=DUMMY2*D1S
               GY1=DUMMY2*D2S
               GZ1=DUMMY2*D3S
               GL1=DUMMY2*D4S
               GM1=DUMMY2*D5S
               GN1=DUMMY2*D6S
               GL2=DUMMY2*DAS
               GM2=DUMMY2*DBS
               GN2=DUMMY2*DCS
            ENDIF
C
C  Sum over charged sites.
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
     1  Sqrt((c2a1*p1x-c3a1*l1*(l1*p1x+m1*p1y+n1*p1z)-c2a2*p2x+c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)+ 
     1      (n1*p1y - m1*p1z)*s1 - (n2*p2y - m2*p2z)*s2 + x1 - x2)**2 + 
     1   (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y+c3a2*m2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (-(n1*p1x) + l1*p1z)*s1 - (-(n2*p2x) + l2*p2z)*s2 + y1-y2)**2 + 
     1   (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z+c3a2*n2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (m1*p1x - l1*p1y)*s1 - (m2*p2x - l2*p2y)*s2 + z1 - z2)**2)
                  DUMMY=DUMMY+HTOKJ*FC(CHARGE(K1),CHARGE(K2),DIST)
                  DUMMY2=HTOKJ*DFC(CHARGE(K1),CHARGE(K2),DIST)
                  IF (GTEST) THEN
                     RDIST=1.0D0/DIST
                     GX1=GX1+DUMMY2*D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GY1=GY1+DUMMY2*D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GZ1=GZ1+DUMMY2*D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GL1=GL1+DUMMY2*D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GM1=GM1+DUMMY2*D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GN1=GN1+DUMMY2*D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GL2=GL2+DUMMY2*DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GM2=GM2+DUMMY2*DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GN2=GN2+DUMMY2*DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                  ENDIF
               ENDDO
            ENDDO
   
            EQUAD=EQUAD+DUMMY
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

C     WRITE(*,'(A,G20.10)') 'energy=',EQUAD
C     PRINT*,'coords:'
C     WRITE(*,'(I6,G20.10)') (J1,X(J1),J1=1,3*NATOMS)
C     PRINT*,'gradient:'
C     WRITE(*,'(I6,G20.10)') (J1,V(J1),J1=1,3*NATOMS)

      RETURN
      END
