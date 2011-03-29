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
C-----------------------------------------------------------------------*
C
C  ENERGY AND GRADIENT FOR THE DZUGUTOV POTENTIAL.
C
C  16.02.2000
C  SERGEI SIMDYANKIN <SSIM@NADA.KTH.SE>
C
C CORRESPONDENCE WITH THE PARAMETERS IN THE FILE "POTENTIAL.F"
C V     = GRAD
C EDZ   = EREAL
C GTEST = GRADT
C STEST 
C
C MEANING OF THE VARIABLES
C DIST - DISTANCE
C EDZ  - TOTAL INTERACTION ENERGY
C VT   - (3*N) VECTOR OF INTERACTION ENERGIES OF EACH PARTICLE
C G    - (N,N) STORAGE FOR GRADIENT INTERMEDIATE BITS
C V    - (3*N) VECTOR - GRADIENT (FORCE ACTING ON A PARTICLE)
! 
C-----------------------------------------------------------------------*
      SUBROUTINE DZPOT(X,V,EDZ,GTEST,STEST)
      USE COMMONS
      IMPLICIT NONE
      LOGICAL GTEST, STEST, DZT
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
     1                 EDZ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY,
     2                 DDUMMY, DDDUMMY, NEARD(NATOMS), DZP72
      INTEGER NEAREST(NATOMS)
      COMMON /DZ/ DZT

      DZP72=DZP7**2
      DZT=.FALSE.
      EDZ=0.0D0
      DO J1=1,NATOMS 
         VT(J1)=0.0D0
         NEARD(J1)=1.0D100
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            G(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.LT.NEARD(J2)) THEN
                  NEARD(J2)=DIST
                  NEAREST(J2)=J1
               ENDIF
               IF (DIST.LT.NEARD(J1)) THEN
                  NEARD(J1)=DIST
                  NEAREST(J1)=J2
               ENDIF
               DDUMMY=0.0D0
               IF (DIST.LT.DZP72) THEN
                  CALL DERPHI(DIST,DUMMY,DDUMMY,DDDUMMY,STEST,DZP1,DZP2,DZP3,DZP4,DZP5,DZP6,DZP7)
                  VT(J1)=VT(J1)+DUMMY
                  VT(J2)=VT(J2)+DUMMY
                  EDZ=EDZ+DUMMY
               ENDIF
               G(J2,J1)=DDUMMY
               G(J1,J2)=G(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.LT.NEARD(J2)) THEN
                  NEARD(J2)=DIST
                  NEAREST(J2)=J1
               ENDIF
               IF (DIST.LT.NEARD(J1)) THEN
                  NEARD(J1)=DIST
                  NEAREST(J1)=J2
               ENDIF
               IF (DIST.LT.DZP72) THEN
                  CALL DERPHI(DIST,DUMMY,DDUMMY,DDDUMMY,STEST,DZP1,DZP2,DZP3,DZP4,DZP5,DZP6,DZP7)
                  VT(J1)=VT(J1)+DUMMY
                  VT(J2)=VT(J2)+DUMMY
                  EDZ=EDZ+DUMMY
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.GTEST) GOTO 10

      DO J1=1,NATOMS
         J3=3*J1
         IF (SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE) THEN
            V(J3-2)=0.0D0
            V(J3-1)=0.0D0
            V(J3)=0.0D0
         ELSE
            DUMMYX=0.0D0
            DUMMYY=0.0D0
            DUMMYZ=0.0D0
            DO J4=1,NATOMS
               J2=3*J4
               XMUL2=G(J4,J1)
               DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
               DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
               DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
            ENDDO
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDIF
      ENDDO

10    DO J1=1,NATOMS 
C        IF (VT(J1).EQ.0.0D0) THEN
         IF (NEARD(J1).GT.2.00D0) THEN
            DZT=.TRUE.
C           DIST=SQRT(X(3*(J1-1)+1)**2+X(3*(J1-1)+2)**2+X(3*(J1-1)+3)**2)
C           X(3*(J1-1)+1)=X(3*(J1-1)+1)*(DIST-1.0D0)/DIST
C           X(3*(J1-1)+2)=X(3*(J1-1)+2)*(DIST-1.0D0)/DIST
C           X(3*(J1-1)+3)=X(3*(J1-1)+3)*(DIST-1.0D0)/DIST
C           PRINT*,'J1,VT,DIST=',J1,VT(J1),DIST
C           PRINT*,'MOVING ATOM ',J1,' NEARER TO ATOM ',NEAREST(J1)
            X(3*(J1-1)+1)=X(3*(NEAREST(J1)-1)+1)+(X(3*(J1-1)+1)-X(3*(NEAREST(J1)-1)+1))/NEARD(J1)
            X(3*(J1-1)+2)=X(3*(NEAREST(J1)-1)+2)+(X(3*(J1-1)+2)-X(3*(NEAREST(J1)-1)+2))/NEARD(J1)
            X(3*(J1-1)+3)=X(3*(NEAREST(J1)-1)+3)+(X(3*(J1-1)+3)-X(3*(NEAREST(J1)-1)+3))/NEARD(J1)
            NEARD(NEAREST(J1))=1.0D0
         ENDIF
      ENDDO

      RETURN
      END

C%%%%%%%
C VALUES OF THE POTENTIAL, FIRST, AND SECOND DERIVATIVES
C IF SECDER = .FALSE. DDPHI = 0 ON OUTPUT
C%%%%%%%
C-----------------------------------------------------------------------*
      SUBROUTINE DERPHI(R2,PHI,DPHI,DDPHI,SECDER,M,A,C,AA,B,D,BB)
      IMPLICIT NONE
      DOUBLE PRECISION R2, PHI, DPHI, DDPHI, M
      LOGICAL SECDER
      INTEGER M2

      DOUBLE PRECISION A, AA, B, BB, C, D 
C     PARAMETER(M=16, A=5.82D0, C=1.10D0, AA=1.87D0, B=1.28D0,  D=0.27D0, BB=1.94D0, M2=8)
C     PARAMETER(M=4,  A=3.00D0, C=0.52D0, AA=1.65D0, B=2.109D0, D=0.55D0, BB=1.94D0, M2=2)
      DOUBLE PRECISION R, V1, DV1, DDV1, V2, DV2, DDV2, DUMMY
      DOUBLE PRECISION EXPO, EXPARG

      M2=M/2
      R = SQRT(R2)
C
C  INITIALIZATION NEEDED FOR SUN!
C
      V1=0.0D0
      V2=0.0D0
      DV1=0.0D0
      DDV1=0.0D0
      DV2=0.0D0
      DDV2=0.0D0

      IF (R .LT. BB) THEN
         EXPARG = D/(R-BB)
C        IF (DABS(EXPARG).GT.708.0D0) THEN 
C           EXPO = 0.0D0
C        ELSE
            EXPO = DEXP(EXPARG)
C        ENDIF
         V2 =  B*EXPO
         DV2 = -V2*D/(R-BB)**2
         IF (SECDER) DDV2 = -DV2*D/(R-BB)**2 + 2*V2*D/(R-BB)**3
         IF (R.LT.AA) THEN 
            EXPARG = C/(R-AA)
C           IF (DABS(EXPARG).GT.708.0D0) THEN 
C              EXPO=0.0D0
C           ELSE
               EXPO = DEXP(EXPARG)
C           ENDIF
            DUMMY=1.0D0/R2**M2
            V1=A*(DUMMY-B)*EXPO
            DV1=-A*M*(DUMMY/R)*EXPO - V1*C/(R-AA)**2
            IF (SECDER) THEN
               DDV1 = A*M*(M+1)*(R**(-M-2))*EXPO 
     1             +  A*M*(R**(-M-1))*EXPO*C/((R-AA)**2) 
     2             -  DV1*C/((R-AA)**2) + 2*V1*C/((R-AA)**3) 
            ENDIF
         ENDIF
      ENDIF
  
      PHI=V1+V2
      DPHI=(DV1+DV2)/R
      IF (SECDER) DDPHI = DDV1 + DDV2

      RETURN
      END                       
