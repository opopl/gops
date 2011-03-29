C   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF OPTIM.
C
C   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
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
C G    - (N,N) STORAGE FOR GRADIENT INTERMEDIATE BITS
C V    - (3*N) VECTOR - GRADIENT (FORCE ACTING ON A PARTICLE)
! 
C-----------------------------------------------------------------------*
      SUBROUTINE DZUGUTOV(NATOMS,X,V,EDZ,GTEST,STEST,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,PARAM7)
      USE MODHESS
      IMPLICIT NONE
      LOGICAL GTEST, STEST
      INTEGER J1, J2, J3, J4, NATOMS, J5, J6
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
     1                 EDZ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY, F(NATOMS,NATOMS), 
     2                 DDUMMY, DDDUMMY,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6, R2(NATOMS,NATOMS),PARAM7
      
      EDZ=0.0D0
      IF (GTEST.OR.STEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            G(J1,J1)=0.0D0
            F(J1,J1)=0.0D0
            R2(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
                     DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               CALL DERPHI(DIST,DUMMY,DDUMMY,DDDUMMY,STEST,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,PARAM7)
               EDZ=EDZ+DUMMY
               R2(J2,J1)=1.0D0/DIST
               R2(J1,J2)=R2(J2,J1)
               G(J2,J1)=DDUMMY
               G(J1,J2)=G(J2,J1)
               IF (STEST) THEN
                  F(J2,J1)=DDDUMMY
                  F(J1,J2)=DDDUMMY
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               CALL DERPHI(DIST,DUMMY,DDUMMY,DDDUMMY,STEST,PARAM1,PARAM2,PARAM3,PARAM4,PARAM5,PARAM6,PARAM7)
               EDZ=EDZ+DUMMY
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.(GTEST.OR.STEST)) RETURN

      DO J1=1,NATOMS
         J3=3*J1
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
      ENDDO

      IF (.NOT.STEST) RETURN
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO J1=1,NATOMS
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,NATOMS
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  NEXT ARE THE TERMS WHERE X_I AND X_J ARE ON THE SAME ATOM
C  BUT ARE DIFFERENT, E.G. Y AND Z.
C
      DO J1=1,NATOMS
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,NATOMS
                  DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)* 
     1           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  CASE III, DIFFERENT ATOMS, SAME CARTESIAN COORDINATE.
C
      DO J1=1,NATOMS
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,NATOMS
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
C
C  CASE IV: DIFFERENT ATOMS AND DIFFERENT CARTESIAN COORDINATES.
C
      DO J1=1,NATOMS
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,NATOMS
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  SYMMETRISE HESSIAN
C
      DO J1=1,3*NATOMS
         DO J2=J1+1,3*NATOMS
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
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
      DOUBLE PRECISION R2, PHI, DPHI, DDPHI
      LOGICAL SECDER
      INTEGER M2

      DOUBLE PRECISION A, AA, B, BB, C, D, M
C     PARAMETER(M=16, A=5.82D0, C=1.10D0, AA=1.87D0, B=1.28D0,  D=0.27D0, BB=1.94D0)
C     PARAMETER(M=4,  A=3.00D0, C=0.52D0, AA=1.65D0, B=2.109D0, D=0.55D0, BB=1.94D0)
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
