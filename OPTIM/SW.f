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
C
C*************************************************************************
C
C  SUBROUTINE SW2 CALCULATES THE ENERGY, CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY FOR THE SW SI POTENTIAL.
C
C*************************************************************************
C
      SUBROUTINE SWTWO(N, X, V, ENERGY, E3, BOXLX, BOXLY, BOXLZ, GTEST, STEST, LAMBDA)
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, ANV(N,N,3)
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, 
     1                 V(3*N), R(N,N), R4T, VEC(N,N,3),
     2                 G(N,N), BOXLX, BOXLY, BOXLZ, LAMBDA,
     3                 ER(N,N), F(N,N), R2(N,N), E3
      DOUBLE PRECISION BIGA, BIGB, SMALLA
      PARAMETER (BIGA=7.049556277D0, BIGB=0.6022245584D0, SMALLA=1.8D0)
C     COMMON /ANVCOMMON/ ANV

C
C  DEAL WITH ANY ATOMS THAT HAVE LEFT THE BOX.
C
      IF (BULKT.AND.(.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF
C 
C  STORE DISTANCE MATRICES.
C
      ENERGY=0.0D0
C     IF (BULKT.AND.(.NOT.FIXIMAGE)) PRINT*,'RECALCULATING ANV'
C     IF (GTEST) THEN
         DO J1=1, N
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(J2,J1,1)=X(J3+1)-X(J4+1)
               VEC(J2,J1,2)=X(J3+2)-X(J4+2)
               VEC(J2,J1,3)=X(J3+3)-X(J4+3)
               IF (BULKT.AND.(.NOT.FIXIMAGE)) THEN
                  ANV(J2,J1,1)=NINT(VEC(J2,J1,1)/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC(J2,J1,2)/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC(J2,J1,3)/BOXLZ)
               ENDIF
               VEC(J2,J1,1)=VEC(J2,J1,1)-BOXLX*ANV(J2,J1,1)
               VEC(J2,J1,2)=VEC(J2,J1,2)-BOXLY*ANV(J2,J1,2)
               VEC(J2,J1,3)=VEC(J2,J1,3)-BOXLZ*ANV(J2,J1,3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
            ENDDO
         ENDDO 
         DO J1=1,N
            R2(J1,J1)=0.0D0
            R(J1,J1)=0.0D0
            ER(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2(J2,J1)=VEC(J2,J1,1)**2+VEC(J2,J1,2)**2+VEC(J2,J1,3)**2
               R(J2,J1)=SQRT(R2(J2,J1)) 
               IF (R(J2,J1).LT.SMALLA) THEN
                  R2(J2,J1)=1.0D0/R2(J2,J1)
                  R4T=R2(J2,J1)**2 
                  ER(J2,J1)=EXP(1.0D0/(-SMALLA + R(J2,J1)))
                  ENERGY=ENERGY+ER(J2,J1)*(-1.0D0 + BIGB*R4T) 
               ENDIF
               R(J1,J2)=R(J2,J1)
               R2(J1,J2)=R2(J2,J1)
               ER(J1,J2)=ER(J2,J1)
            ENDDO 
         ENDDO
C
C  SWTHREE NEEDS TO KNOW VEC
C
C     ELSE
C        DO J1=1,N
C           J3=3*(J1-1)
C           DO J2=J1+1,N
C              J4=3*(J2-1)
C              VEC1=X(J3+1)-X(J4+1)
C              VEC2=X(J3+2)-X(J4+2)
C              VEC3=X(J3+3)-X(J4+3)
C              IF (BULKT.AND.(.NOT.FIXIMAGE)) THEN
C                 ANV(J2,J1,1)=NINT(VEC1/BOXLX)
C                 ANV(J2,J1,2)=NINT(VEC2/BOXLY)
C                 ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
C              ENDIF
C              VEC1=VEC1-BOXLX*ANV(J2,J1,1)
C              VEC2=VEC2-BOXLY*ANV(J2,J1,2)
C              VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
C              R2T=VEC1**2+VEC2**2+VEC3**2
C              DIST=SQRT(R2T)
C              R4T=R2T**2
C              IF (DIST.LT.SMALLA) ENERGY=ENERGY+EXP(1.0D0/(-SMALLA + DIST))*(-1.0D0 + BIGB/R4T)
C           ENDDO
C        ENDDO
C     ENDIF

      ENERGY=BIGA*ENERGY

      IF (GTEST) CALL SWTWOG(G,R,R2,ER,V,X,N,VEC)
      
      IF (STEST) CALL SWTWOS(G,F,R,R2,ER,X,N,VEC)

      CALL SWTHREE(N, X, V, E3, GTEST, STEST, LAMBDA, R, R2, G, F, VEC)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE SWTWOG(G,R,R2,ER,V,X,N,VEC)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION G(N,N), R(N,N), R2(N,N),
     1                 V(3*N), X(3*N), DUMMY, ER(N,N), VEC(N,N,3)
      DOUBLE PRECISION BIGA, BIGB, SMALLA
      PARAMETER (BIGA=7.049556277D0, BIGB=0.6022245584D0, SMALLA=1.8D0)
C
C  CALCULATE THE G TENSOR.
C
      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            IF (R(J2,J1).LT.SMALLA) THEN
               G(J2,J1)=(BIGA*ER(J2,J1)*(-4*SMALLA**2*BIGB + 
     1                    (-1.0D0 + 8.0D0*SMALLA)*BIGB*R(J2,J1) - 4*BIGB/R2(J2,J1) + R(J2,J1)**5))*
     2                    R2(J2,J1)**3/(SMALLA - R(J2,J1))**2
            ELSE
               G(J2,J1)=0.0D0
            ENDIF
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
C  FROM HERE ON DOWN THE CODE IS SYSTEM-INDEPENDENT!
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               IF (R(J4,J1).LT.SMALLA) DUMMY=DUMMY+G(J4,J1)*VEC(J4,J1,J2)
            ENDDO
            V(J3)=DUMMY
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE SWTWOS(G,F,R,R2,ER,X,N,VEC)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(N,N), R(N,N), VEC(N,N,3),
     1                 F(N,N), 
     2                 X(3*N),DUMMY, R2(N,N), ER(N,N)
      DOUBLE PRECISION BIGA, BIGB, SMALLA
      PARAMETER (BIGA=7.049556277D0, BIGB=0.6022245584D0, SMALLA=1.8D0)

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            IF (R(J2,J1).LT.SMALLA) THEN
               F(J2,J1)=BIGA*ER(J2,J1)*((-1.0D0 + 2*SMALLA - 2*R(J2,J1))/(SMALLA - R(J2,J1))**4 + 
     1                  R2(J2,J1)**3*(BIGB*(20.0D0 + (R(J2,J1)*(2.0D0*(4.0D0*SMALLA - 5.0D0*R(J2,J1))*
     1                  (SMALLA - R(J2,J1)) + R(J2,J1)))/(SMALLA - R(J2,J1))**4))) - G(J2,J1)
            ELSE
               F(J2,J1)=0.0D0
            ENDIF
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               IF (R(J4,J1).LT.SMALLA) DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  NEXT ARE THE TERMS WHERE X_I AND X_J ARE ON THE SAME ATOM
C  BUT ARE DIFFERENT, E.G. Y AND Z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,N
                  IF (R(J5,J1).LT.SMALLA) DUMMY=DUMMY + 
     1                F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  CASE III, DIFFERENT ATOMS, SAME CARTESIAN COORDINATE.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               IF (R(J4,J1).LT.SMALLA) THEN
                  HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1) 
               ELSE
                  HESS(3*(J4-1)+J2,J3)=0.0D0
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C  CASE IV: DIFFERENT ATOMS AND DIFFERENT CARTESIAN COORDINATES.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  IF (R(J4,J1).LT.SMALLA) THEN
                     HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                  ELSE
                     HESS(J6,J3)=0.0D0
                  ENDIF
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  SYMMETRISE HESSIAN
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE SWLATMIN(N,XC,BOXLX,BOXLY,BOXLZ,V,LAMBDA)
      IMPLICIT NONE
      INTEGER N, NCOUNT, J1
      DOUBLE PRECISION XC(3*N),F1,F2,F3,GRAD,SECOND,XSAVE(3*N),
     1                 P2, P3,V(3*N),LAMBDA,BOXLX,BOXLY,BOXLZ,DIF,TEMP1
C
C  VALUE OF DIF IS THE ORDER OF MAGNITUDE TO WHICH THE LATTICE
C  CONSTANT CAN BE OPTIMISED. SETTING IT SMALLER THAN 10^(-7)
C  CAUSES NUMERICAL PROBLEMS ON THE DEC.
C
      DIF=1.0D-7
      NCOUNT=1
      DO 20 J1=1,3*N
         XSAVE(J1)=XC(J1)
20    CONTINUE

10    TEMP1=BOXLX
      BOXLX=TEMP1+DIF
      BOXLY=TEMP1+DIF
      BOXLZ=TEMP1+DIF

      CALL SWTWO  (N, XC, V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC, V, P3, .FALSE., .FALSE.,LAMBDA)
      F1=P2+P3

      BOXLX=TEMP1-DIF
      BOXLY=TEMP1-DIF
      BOXLZ=TEMP1-DIF
      DO 30 J1=1,3*N
         XC(J1)=XSAVE(J1)
30    CONTINUE
      CALL SWTWO  (N, XC,  V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC,  V, P3, .FALSE., .FALSE.,LAMBDA)
      F2=P2+P3

      GRAD=(F1-F2)/(2.0D0*DIF)

      BOXLX=TEMP1
      BOXLY=TEMP1
      BOXLZ=TEMP1
      DO 40 J1=1,3*N
         XC(J1)=XSAVE(J1)
40    CONTINUE
      CALL SWTWO  (N, XC, V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC, V, P3, .FALSE., .FALSE.,LAMBDA)
      F1=P2+P3

      BOXLX=TEMP1+2.0D0*DIF
      BOXLY=TEMP1+2.0D0*DIF
      BOXLZ=TEMP1+2.0D0*DIF
      DO 50 J1=1,3*N
         XC(J1)=XSAVE(J1)
50    CONTINUE
      CALL SWTWO  (N, XC, V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC, V, P3, .FALSE., .FALSE.,LAMBDA)
      F2=P2+P3

      BOXLX=TEMP1-2.0D0*DIF
      BOXLY=TEMP1-2.0D0*DIF
      BOXLZ=TEMP1-2.0D0*DIF
      DO 60 J1=1,3*N
         XC(J1)=XSAVE(J1)
60    CONTINUE
      CALL SWTWO  (N, XC, V, P2, P3, BOXLX, BOXLY, BOXLZ, .FALSE., .FALSE., LAMBDA)
C     CALL SWTHREE(N, XC, V, P3, .FALSE., .FALSE.,LAMBDA)
      F3=P2+P3

      SECOND=(F3+F2-2.0D0*F1)/(4.0D0*DIF*DIF)
      PRINT*,'ENERGY FOR LATTICE CYCLE ',NCOUNT,' IS ',F1
      PRINT*,'GRADIENT WRT BOX LENGTH=',GRAD
      PRINT*,'SECOND DERIVATIVE WRT BOX LENGTH=',SECOND
      NCOUNT=NCOUNT+1
      IF (DABS(GRAD/SECOND).GT.1.0D-4) THEN
         BOXLX=BOXLX-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         BOXLY=BOXLY-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         BOXLZ=BOXLZ-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         PRINT*,'STEP=',-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         GOTO 10
      ELSE
         BOXLX=BOXLX-GRAD/SECOND
         BOXLY=BOXLY-GRAD/SECOND
         BOXLZ=BOXLZ-GRAD/SECOND
         PRINT*,'STEP=',-GRAD/SECOND
         IF (DABS(GRAD/SECOND).GT.1.0D-6) GOTO 10
      ENDIF
      RETURN
      END
C
C*************************************************************************
C
C  SUBROUTINE SW2 CALCULATES THE ENERGY, CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY FOR THE SW SI POTENTIAL.
C
C*************************************************************************
C
      SUBROUTINE SWTHREE(N, X, V, ENERGY, GTEST, STEST, LAMBDA, R, R2, G, F, VEC)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, K1, K2, K3
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, GAMMA2, ER2131,
     1                 V(3*N), R(N,N), R122, R123, R132, R133, R12S, R13S,
     2                 G(N,N), R12DOTR13PT, R12DOTR13, R134, R124, ERLDOT,
     3                 ER(N,N), F(N,N), R2(N,N), THIRD,
     4                 ER21, ER31, R12, R13, R23, R13S3, R12S3, DOT1,
     5                 R13S4, R12S4, DOT2, ERL, 
     6                 VEC(N,N,3), V211, V212, V213, V311, V312, V313, V2112, V2122, V2132,
     7                 V3112, V3122, V3132
      DOUBLE PRECISION LAMBDA, GAMMA, SMALLA
      PARAMETER (GAMMA=1.2D0, SMALLA=1.8D0, THIRD=0.33333333333333333333D0)
      
C 
C  STORE THINGS
C
      DO J1=1,N
         ER(J1,J1)=0.0D0
         DO J2=J1+1,N
            IF (R(J2,J1).LT.SMALLA) THEN
               ER(J2,J1)=EXP(GAMMA/(-SMALLA + R(J2,J1)))
            ELSE
               ER(J2,J1)=0.0D0
            ENDIF
            ER(J1,J2)=ER(J2,J1)
C           WRITE(*,'(A,2I4,2E20.10)') 'J1,J2,R,ER=',J1,J2,R(J2,J1),ER(J2,J1)
         ENDDO
      ENDDO

      ENERGY=0.0D0
      DO J1=1,N
         K1=3*(J1-1)
         DO J2=1,N
            IF ((R(J2,J1).LT.SMALLA).AND.(J2.NE.J1)) THEN
               K2=3*(J2-1)
               ER21=ER(J2,J1)
               R12=R(J2,J1)
               V211=VEC(J2,J1,1)
               V212=VEC(J2,J1,2)
               V213=VEC(J2,J1,3)
               V2112=V211**2
               V2122=V212**2
               V2132=V213**2
               DO J3=1,N
                  IF (((R(J3,J1).LT.SMALLA).AND.(J3.NE.J1)).AND.(J3.NE.J2)) THEN
                     K3=3*(J3-1)
                     ER31=ER(J3,J1)
                     ER2131=ER21*ER31
                     R13=R(J3,J1)
                     R23=R(J3,J2)
                     V311=VEC(J3,J1,1)
                     V312=VEC(J3,J1,2)
                     V313=VEC(J3,J1,3)
                     V3112=V311**2
                     V3122=V312**2
                     V3132=V313**2
                     R12DOTR13=(V211*V311+V212*V312+V213*V313)/(R12*R13)
                     DOT1=R12DOTR13**2
                     R12DOTR13PT=(R12DOTR13+THIRD)
                     DOT2=R12DOTR13PT**2

                     ENERGY=ENERGY+ER2131*DOT2

                     IF (GTEST) THEN
                        R122=R12**2
                        R132=R13**2
                        R12S=(R12 - SMALLA)**2
                        R13S=(R13 - SMALLA)**2
                        ERLDOT=ER2131*LAMBDA*R12DOTR13PT

                       V(K1+1)=V(K1+1)+ (ERLDOT*
     1    ((2*R13*(R12 - R12DOTR13*R13)*V211 + 
     1         2*R12*(-(R12*R12DOTR13) + R13)*V311)/(R122*R132) + 
     1      GAMMA*R12DOTR13PT*(-(V211/(R12*R12S)) - 
     1         V311/(R13*R13S))))/2.0D0
                        V(K1+2)=V(K1+2)+ (ERLDOT*
     1    ((2*R13*(R12 - R12DOTR13*R13)*V212 + 
     1         2*R12*(-(R12*R12DOTR13) + R13)*V312)/(R122*R132) + 
     1      GAMMA*R12DOTR13PT*(-(V212/(R12*R12S)) - 
     1         V312/(R13*R13S))))/2.0D0
                        V(K1+3)=V(K1+3)+ (ERLDOT*
     1    ((2*R13*(R12 - R12DOTR13*R13)*V213 + 
     1         2*R12*(-(R12*R12DOTR13) + R13)*V313)/(R122*R132) + 
     1      GAMMA*R12DOTR13PT*(-(V213/(R12*R12S)) - 
     1         V313/(R13*R13S))))/2.0D0
                        V(K2+1)=V(K2+1)+ (ERLDOT*
     1    ((2*R12DOTR13 + (GAMMA*R12*R12DOTR13PT)/R12S)*
     1       V211 - (2*R12*V311)/R13))/(2.*R122)
                        V(K2+2)=V(K2+2)+ (ERLDOT*
     1    ((2*R12DOTR13 + (GAMMA*R12*R12DOTR13PT)/R12S)*
     1       V212 - (2*R12*V312)/R13))/(2.*R122)
                        V(K2+3)=V(K2+3)+ (ERLDOT*
     1    ((2*R12DOTR13 + (GAMMA*R12*R12DOTR13PT)/R12S)*
     1       V213 - (2*R12*V313)/R13))/(2.*R122)
                        V(K3+1)=V(K3+1)+ (ERLDOT*
     1    ((-2*R13*V211)/R12 + 
     1      (2*R12DOTR13 + (GAMMA*R12DOTR13PT*R13)/R13S)*
     1       V311))/(2.*R132)
                       V(K3+2)=V(K3+2)+ (ERLDOT*
     1    ((-2*R13*V212)/R12 + 
     1      (2*R12DOTR13 + (GAMMA*R12DOTR13PT*R13)/R13S)*
     1       V312))/(2.*R132)
                        V(K3+3)=V(K3+3)+ (ERLDOT*
     1    ((-2*R13*V213)/R12 + 
     1      (2*R12DOTR13 + (GAMMA*R12DOTR13PT*R13)/R13S)*
     1       V313))/(2.*R132)

                     ENDIF
                     IF (STEST) THEN
                        ERL=ER2131*LAMBDA
                        GAMMA2=GAMMA**2
                        R122=R12**2
                        R123=R12*R122
                        R124=R12*R123
                        R132=R13**2
                        R133=R13*R132
                        R134=R13*R133
                        R12S=(R12 - SMALLA)**2
                        R13S=(R13 - SMALLA)**2
                        R12S3=R12S*(R12 - SMALLA)
                        R13S3=R13S*(R13 - SMALLA)
                        R12S4=R12S3*(R12 - SMALLA)
                        R13S4=R13S3*(R13 - SMALLA)

                        HESS(K1+1,K1+1)=HESS(K1+1,K1+1)+ (ERL*((2*(R13*(-R12 + R12DOTR13*R13)*V211 + 
     1            R12*(R12*R12DOTR13 - R13)*V311)**2)/(R124*R134)
     1       + (4*GAMMA*R12DOTR13PT*
     1         (R13*(R12 - R12DOTR13*R13)*V211 + 
     1           R12*(-(R12*R12DOTR13) + R13)*V311)*
     1         (-(V211/(R12*R12S)) - 
     1           V311/(R13*R13S)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (V211/(R12*R12S) + 
     1          V311/(R13*R13S))**2 + 
     1      GAMMA*DOT2*
     1       (-(1/(R12*R12S)) - 1/(R13*R13S) + 
     1         ((3*R12 - SMALLA)*V2112)/
     1          (R123*R12S3) + 
     1         ((3*R13 - SMALLA)*V3112)/(R133*R13S3))
     1       - (2*R12DOTR13PT*(R133*(2*R12 - 3*R12DOTR13*R13)*
     1            V2112 + 
     1           2*R12*R13*(R122 - R12*R12DOTR13*R13 + R132)*V211*
     1            V311 + 
     1           R122*(R132*
     1               (R122*R12DOTR13 - 2*R12*R13 + R12DOTR13*R132) + 
     1              R12*(-3*R12*R12DOTR13 + 2*R13)*V3112)))/
     1       (R124*R134)))/2.
                        HESS(K1+1,K1+2)=HESS(K1+1,K1+2)+ (ERL*((2*(R13*(-R12 + R12DOTR13*R13)*V211 + 
     1           R12*(R12*R12DOTR13 - R13)*V311)*
     1         (R13*(-R12 + R12DOTR13*R13)*V212 + 
     1           R12*(R12*R12DOTR13 - R13)*V312))/(R124*R134) + 
     1      (2*GAMMA*R12DOTR13PT*
     1         (-(V211/(R12*R12S)) - 
     1           V311/(R13*R13S))*
     1         (R13*(R12 - R12DOTR13*R13)*V212 + 
     1           R12*(-(R12*R12DOTR13) + R13)*V312))/(R122*R132)
     1       + (2*GAMMA*R12DOTR13PT*
     1         (R13*(R12 - R12DOTR13*R13)*V211 + 
     1           R12*(-(R12*R12DOTR13) + R13)*V311)*
     1         (-(V212/(R12*R12S)) - 
     1           V312/(R13*R13S)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (-(V211/(R12*R12S)) - 
     1         V311/(R13*R13S))*
     1       (-(V212/(R12*R12S)) - 
     1         V312/(R13*R13S)) + 
     1      GAMMA*DOT2*
     1       (((3*R12 - SMALLA)*V211*V212)/
     1          (R123*R12S3) + 
     1         ((3*R13 - SMALLA)*V311*V312)/
     1          (R133*R13S3)) + 
     1      (2*R12DOTR13PT*(R12*V311*
     1            (-(R13*(R122 - R12*R12DOTR13*R13 + R132)*
     1                 V212) + 
     1              R122*(3*R12*R12DOTR13 - 2*R13)*V312) - 
     1           R13*V211*
     1            (R132*(2*R12 - 3*R12DOTR13*R13)*V212 + 
     1              R12*(R122 - R12*R12DOTR13*R13 + R132)*V312)))/
     1       (R124*R134)))/2.
                        HESS(K1+1,K1+3)=HESS(K1+1,K1+3)+ (ERL*((2*(R13*(-R12 + R12DOTR13*R13)*V211 + 
     1           R12*(R12*R12DOTR13 - R13)*V311)*
     1         (R13*(-R12 + R12DOTR13*R13)*V213 + 
     1           R12*(R12*R12DOTR13 - R13)*V313))/(R124*R134) + 
     1      (2*GAMMA*R12DOTR13PT*
     1         (-(V211/(R12*R12S)) - 
     1           V311/(R13*R13S))*
     1         (R13*(R12 - R12DOTR13*R13)*V213 + 
     1           R12*(-(R12*R12DOTR13) + R13)*V313))/(R122*R132)
     1       + (2*GAMMA*R12DOTR13PT*
     1         (R13*(R12 - R12DOTR13*R13)*V211 + 
     1           R12*(-(R12*R12DOTR13) + R13)*V311)*
     1         (-(V213/(R12*R12S)) - 
     1           V313/(R13*R13S)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (-(V211/(R12*R12S)) - 
     1         V311/(R13*R13S))*
     1       (-(V213/(R12*R12S)) - 
     1         V313/(R13*R13S)) + 
     1      GAMMA*DOT2*
     1       (((3*R12 - SMALLA)*V211*V213)/
     1          (R123*R12S3) + 
     1         ((3*R13 - SMALLA)*V311*V313)/
     1          (R133*R13S3)) + 
     1      (2*R12DOTR13PT*(R12*V311*
     1            (-(R13*(R122 - R12*R12DOTR13*R13 + R132)*
     1                 V213) + 
     1              R122*(3*R12*R12DOTR13 - 2*R13)*V313) - 
     1           R13*V211*
     1            (R132*(2*R12 - 3*R12DOTR13*R13)*V213 + 
     1              R12*(R122 - R12*R12DOTR13*R13 + R132)*V313)))/
     1       (R124*R134)))/2.
                        HESS(K1+1,K2+1)=HESS(K1+1,K2+1)+ (ERL*((-2*DOT1 + 
     1         R12DOTR13*(-6*R12DOTR13PT + (2*R12)/R13 - 
     1            (4*GAMMA*R12*R12DOTR13PT)/R12S) + 
     1         R12*R12DOTR13PT*(2/R13 - 
     1            (GAMMA2*R12*R12DOTR13PT)/R12S4 - 
     1            (3*GAMMA*R12*R12DOTR13PT)/R12S3 + 
     1            (2*GAMMA*R12)/(R13*R12S) + 
     1            (GAMMA*R12DOTR13PT*SMALLA)/R12S3))*
     1       V2112 + (R12*
     1         (4*(R12DOTR13 + R12DOTR13PT)*R13 + 
     1           R12*(-2 - 2*DOT1 - 2*R12DOTR13*R12DOTR13PT + 
     1              (4*GAMMA*R12DOTR13PT*R13)/R12S - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S) + 
     1           (GAMMA*R122*R12DOTR13PT*
     1              (-2*R12DOTR13 - (GAMMA*R12DOTR13PT*R13)/R13S)
     1              )/R12S)*V211*V311)/R132 + 
     1      R122*(R12DOTR13PT*
     1          (2*R12DOTR13 - (2*R12)/R13 + 
     1            (GAMMA*R122*R12DOTR13PT)/R12S3 - 
     1            (GAMMA*R12*R12DOTR13PT*SMALLA)/R12S3) + 
     1         (2*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V3112)/R133)))/(2.*R124)
                        HESS(K1+1,K2+2)=HESS(K1+1,K2+2)+ (ERL*(V211*
     1       ((-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12*R12DOTR13)/R13 + (2*R12*R12DOTR13PT)/R13 - 
     1            (GAMMA2*R122*DOT2)/R12S4 - 
     1            (3*GAMMA*R122*DOT2)/R12S3 - 
     1            (4*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R122*R12DOTR13PT)/(R13*R12S) + 
     1            (GAMMA*R12*DOT2*SMALLA)/R12S3)*
     1          V212 + (2*R12*
     1            (-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V312)/R132) + 
     1      (R12*V311*(R13*
     1            (-2*R12*DOT1 - 2*R12*R12DOTR13*R12DOTR13PT + 
     1              2*R12DOTR13*R13 + 2*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R122*R12DOTR13*R12DOTR13PT)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R12S - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R13S - 
     1              (GAMMA2*R122*DOT2*R13)/
     1               (R12S*R13S))*V212 + 
     1           2*R12*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V312))/R133))/(2.*R124)
                        HESS(K1+1,K2+3)=HESS(K1+1,K2+3)+ (ERL*(V211*
     1       ((-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12*R12DOTR13)/R13 + (2*R12*R12DOTR13PT)/R13 - 
     1            (GAMMA2*R122*DOT2)/R12S4 - 
     1            (3*GAMMA*R122*DOT2)/R12S3 - 
     1            (4*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R122*R12DOTR13PT)/(R13*R12S) + 
     1            (GAMMA*R12*DOT2*SMALLA)/R12S3)*
     1          V213 + (2*R12*
     1            (-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V313)/R132) + 
     1      (R12*V311*(R13*
     1            (-2*R12*DOT1 - 2*R12*R12DOTR13*R12DOTR13PT + 
     1              2*R12DOTR13*R13 + 2*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R122*R12DOTR13*R12DOTR13PT)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R12S - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R13S - 
     1              (GAMMA2*R122*DOT2*R13)/
     1               (R12S*R13S))*V213 + 
     1           2*R12*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V313))/R133))/(2.*R124)
                        HESS(K1+1,K3+1)=HESS(K1+1,K3+1)+ (ERL*(R12DOTR13PT*R132*
     1       (2*R12DOTR13 - (2*R13)/R12 + 
     1         (GAMMA*R12DOTR13PT*R132)/R13S3 - 
     1         (GAMMA*R12DOTR13PT*R13*SMALLA)/R13S3) + 
     1      (2*R132*((R12DOTR13 + R12DOTR13PT)*R13 + 
     1           R12*(-1 + (GAMMA*R12DOTR13PT*R13)/R12S))*
     1         V2112)/R123 + 
     1      (R13*(2*R13*(-1 - DOT1 - R12DOTR13*R12DOTR13PT - 
     1              (GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S) + 
     1           R12*(4*R12DOTR13 + 4*R12DOTR13PT - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R12S + 
     1              (4*GAMMA*R12DOTR13PT*R13)/R13S - 
     1              (GAMMA2*DOT2*R132)/
     1               (R12S*R13S)))*V211*
     1         V311)/R122 + 
     1      (-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1         (2*R12DOTR13*R13)/R12 + (2*R12DOTR13PT*R13)/R12 - 
     1         (GAMMA2*DOT2*R132)/R13S4 - 
     1         (3*GAMMA*DOT2*R132)/R13S3 - 
     1         (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1         (2*GAMMA*R12DOTR13PT*R132)/(R12*R13S) + 
     1         (GAMMA*DOT2*R13*SMALLA)/R13S3)*
     1       V3112))/(2.*R134)
                        HESS(K1+1,K3+2)=HESS(K1+1,K3+2)+ (ERL*((R13*V211*
     1         (2*R13*(-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V212 + 
     1           R12*(2*R12*R12DOTR13 + 2*R12*R12DOTR13PT - 
     1              2*DOT1*R13 - 2*R12DOTR13*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R13S - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R132)/
     1               R13S - 
     1              (GAMMA2*R12*DOT2*R132)/
     1               (R12S*R13S))*V312))/
     1       R123 + V311*
     1       ((2*R13*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V212)/R122 + 
     1         (-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12DOTR13*R13)/R12 + (2*R12DOTR13PT*R13)/R12 - 
     1            (GAMMA2*DOT2*R132)/R13S4 - 
     1            (3*GAMMA*DOT2*R132)/R13S3 - 
     1            (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (2*GAMMA*R12DOTR13PT*R132)/(R12*R13S) + 
     1            (GAMMA*DOT2*R13*SMALLA)/R13S3)*
     1          V312)))/(2.*R134)
                        HESS(K1+1,K3+3)=HESS(K1+1,K3+3)+ (ERL*((R13*V211*
     1         (2*R13*(-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V213 + 
     1           R12*(2*R12*R12DOTR13 + 2*R12*R12DOTR13PT - 
     1              2*DOT1*R13 - 2*R12DOTR13*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R13S - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R132)/
     1               R13S - 
     1              (GAMMA2*R12*DOT2*R132)/
     1               (R12S*R13S))*V313))/
     1       R123 + V311*
     1       ((2*R13*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V213)/R122 + 
     1         (-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12DOTR13*R13)/R12 + (2*R12DOTR13PT*R13)/R12 - 
     1            (GAMMA2*DOT2*R132)/R13S4 - 
     1            (3*GAMMA*DOT2*R132)/R13S3 - 
     1            (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (2*GAMMA*R12DOTR13PT*R132)/(R12*R13S) + 
     1            (GAMMA*DOT2*R13*SMALLA)/R13S3)*
     1          V313)))/(2.*R134)
                        HESS(K1+2,K1+2)=HESS(K1+2,K1+2)+ (ERL*((2*(R13*(-R12 + R12DOTR13*R13)*V212 + 
     1            R12*(R12*R12DOTR13 - R13)*V312)**2)/(R124*R134)
     1       + (4*GAMMA*R12DOTR13PT*
     1         (R13*(R12 - R12DOTR13*R13)*V212 + 
     1           R12*(-(R12*R12DOTR13) + R13)*V312)*
     1         (-(V212/(R12*R12S)) - 
     1           V312/(R13*R13S)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (V212/(R12*R12S) + 
     1          V312/(R13*R13S))**2 + 
     1      GAMMA*DOT2*
     1       (-(1/(R12*R12S)) - 1/(R13*R13S) + 
     1         ((3*R12 - SMALLA)*V2122)/
     1          (R123*R12S3) + 
     1         ((3*R13 - SMALLA)*V3122)/(R133*R13S3))
     1       - (2*R12DOTR13PT*(R133*(2*R12 - 3*R12DOTR13*R13)*
     1            V2122 + 
     1           2*R12*R13*(R122 - R12*R12DOTR13*R13 + R132)*V212*
     1            V312 + 
     1           R122*(R132*
     1               (R122*R12DOTR13 - 2*R12*R13 + R12DOTR13*R132) + 
     1              R12*(-3*R12*R12DOTR13 + 2*R13)*V3122)))/
     1       (R124*R134)))/2.
                        HESS(K1+2,K1+3)=HESS(K1+2,K1+3)+ (ERL*((2*(R13*(-R12 + R12DOTR13*R13)*V212 + 
     1           R12*(R12*R12DOTR13 - R13)*V312)*
     1         (R13*(-R12 + R12DOTR13*R13)*V213 + 
     1           R12*(R12*R12DOTR13 - R13)*V313))/(R124*R134) + 
     1      (2*GAMMA*R12DOTR13PT*
     1         (-(V212/(R12*R12S)) - 
     1           V312/(R13*R13S))*
     1         (R13*(R12 - R12DOTR13*R13)*V213 + 
     1           R12*(-(R12*R12DOTR13) + R13)*V313))/(R122*R132)
     1       + (2*GAMMA*R12DOTR13PT*
     1         (R13*(R12 - R12DOTR13*R13)*V212 + 
     1           R12*(-(R12*R12DOTR13) + R13)*V312)*
     1         (-(V213/(R12*R12S)) - 
     1           V313/(R13*R13S)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (-(V212/(R12*R12S)) - 
     1         V312/(R13*R13S))*
     1       (-(V213/(R12*R12S)) - 
     1         V313/(R13*R13S)) + 
     1      GAMMA*DOT2*
     1       (((3*R12 - SMALLA)*V212*V213)/
     1          (R123*R12S3) + 
     1         ((3*R13 - SMALLA)*V312*V313)/
     1          (R133*R13S3)) + 
     1      (2*R12DOTR13PT*(R12*V312*
     1            (-(R13*(R122 - R12*R12DOTR13*R13 + R132)*
     1                 V213) + 
     1              R122*(3*R12*R12DOTR13 - 2*R13)*V313) - 
     1           R13*V212*
     1            (R132*(2*R12 - 3*R12DOTR13*R13)*V213 + 
     1              R12*(R122 - R12*R12DOTR13*R13 + R132)*V313)))/
     1       (R124*R134)))/2.
                        HESS(K1+2,K2+1)=HESS(K1+2,K2+1)+ (ERL*((2*R12*V311*
     1         (R13*(-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V212 + 
     1           R12*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V312))/R133 + 
     1      V211*((-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12*R12DOTR13)/R13 + (2*R12*R12DOTR13PT)/R13 - 
     1            (GAMMA2*R122*DOT2)/R12S4 - 
     1            (3*GAMMA*R122*DOT2)/R12S3 - 
     1            (4*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R122*R12DOTR13PT)/(R13*R12S) + 
     1            (GAMMA*R12*DOT2*SMALLA)/R12S3)*
     1          V212 + (R12*
     1            (-2*R12*DOT1 - 2*R12*R12DOTR13*R12DOTR13PT + 
     1              2*R12DOTR13*R13 + 2*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R122*R12DOTR13*R12DOTR13PT)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R12S - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R13S - 
     1              (GAMMA2*R122*DOT2*R13)/
     1               (R12S*R13S))*V312)/
     1          R132)))/(2.*R124)
                        HESS(K1+2,K2+2)=HESS(K1+2,K2+2)+ (ERL*((-2*DOT1 + 
     1         R12DOTR13*(-6*R12DOTR13PT + (2*R12)/R13 - 
     1            (4*GAMMA*R12*R12DOTR13PT)/R12S) + 
     1         R12*R12DOTR13PT*(2/R13 - 
     1            (GAMMA2*R12*R12DOTR13PT)/R12S4 - 
     1            (3*GAMMA*R12*R12DOTR13PT)/R12S3 + 
     1            (2*GAMMA*R12)/(R13*R12S) + 
     1            (GAMMA*R12DOTR13PT*SMALLA)/R12S3))*
     1       V2122 + (R12*
     1         (4*(R12DOTR13 + R12DOTR13PT)*R13 + 
     1           R12*(-2 - 2*DOT1 - 2*R12DOTR13*R12DOTR13PT + 
     1              (4*GAMMA*R12DOTR13PT*R13)/R12S - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S) + 
     1           (GAMMA*R122*R12DOTR13PT*
     1              (-2*R12DOTR13 - (GAMMA*R12DOTR13PT*R13)/R13S)
     1              )/R12S)*V212*V312)/R132 + 
     1      R122*(R12DOTR13PT*
     1          (2*R12DOTR13 - (2*R12)/R13 + 
     1            (GAMMA*R122*R12DOTR13PT)/R12S3 - 
     1            (GAMMA*R12*R12DOTR13PT*SMALLA)/R12S3) + 
     1         (2*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V3122)/R133)))/(2.*R124)
                        HESS(K1+2,K2+3)=HESS(K1+2,K2+3)+ (ERL*(V212*
     1       ((-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12*R12DOTR13)/R13 + (2*R12*R12DOTR13PT)/R13 - 
     1            (GAMMA2*R122*DOT2)/R12S4 - 
     1            (3*GAMMA*R122*DOT2)/R12S3 - 
     1            (4*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R122*R12DOTR13PT)/(R13*R12S) + 
     1            (GAMMA*R12*DOT2*SMALLA)/R12S3)*
     1          V213 + (2*R12*
     1            (-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V313)/R132) + 
     1      (R12*V312*(R13*
     1            (-2*R12*DOT1 - 2*R12*R12DOTR13*R12DOTR13PT + 
     1              2*R12DOTR13*R13 + 2*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R122*R12DOTR13*R12DOTR13PT)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R12S - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R13S - 
     1              (GAMMA2*R122*DOT2*R13)/
     1               (R12S*R13S))*V213 + 
     1           2*R12*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V313))/R133))/(2.*R124)
                        HESS(K1+2,K3+1)=HESS(K1+2,K3+1)+ (ERL*((2*R13*V211*
     1         (R13*(-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V212 + 
     1           R12*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V312))/R123 + 
     1      V311*((R13*(2*R12*R12DOTR13 + 2*R12*R12DOTR13PT - 
     1              2*DOT1*R13 - 2*R12DOTR13*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R13S - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R132)/
     1               R13S - 
     1              (GAMMA2*R12*DOT2*R132)/
     1               (R12S*R13S))*V212)/
     1          R122 + (-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12DOTR13*R13)/R12 + (2*R12DOTR13PT*R13)/R12 - 
     1            (GAMMA2*DOT2*R132)/R13S4 - 
     1            (3*GAMMA*DOT2*R132)/R13S3 - 
     1            (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (2*GAMMA*R12DOTR13PT*R132)/(R12*R13S) + 
     1            (GAMMA*DOT2*R13*SMALLA)/R13S3)*
     1          V312)))/(2.*R134)
                        HESS(K1+2,K3+2)=HESS(K1+2,K3+2)+ (ERL*(R12DOTR13PT*R132*
     1       (2*R12DOTR13 - (2*R13)/R12 + 
     1         (GAMMA*R12DOTR13PT*R132)/R13S3 - 
     1         (GAMMA*R12DOTR13PT*R13*SMALLA)/R13S3) + 
     1      (2*R132*((R12DOTR13 + R12DOTR13PT)*R13 + 
     1           R12*(-1 + (GAMMA*R12DOTR13PT*R13)/R12S))*
     1         V2122)/R123 + 
     1      (R13*(2*R13*(-1 - DOT1 - R12DOTR13*R12DOTR13PT - 
     1              (GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S) + 
     1           R12*(4*R12DOTR13 + 4*R12DOTR13PT - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R12S + 
     1              (4*GAMMA*R12DOTR13PT*R13)/R13S - 
     1              (GAMMA2*DOT2*R132)/
     1               (R12S*R13S)))*V212*
     1         V312)/R122 + 
     1      (-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1         (2*R12DOTR13*R13)/R12 + (2*R12DOTR13PT*R13)/R12 - 
     1         (GAMMA2*DOT2*R132)/R13S4 - 
     1         (3*GAMMA*DOT2*R132)/R13S3 - 
     1         (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1         (2*GAMMA*R12DOTR13PT*R132)/(R12*R13S) + 
     1         (GAMMA*DOT2*R13*SMALLA)/R13S3)*
     1       V3122))/(2.*R134)
                        HESS(K1+2,K3+3)=HESS(K1+2,K3+3)+ (ERL*((R13*V212*
     1         (2*R13*(-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V213 + 
     1           R12*(2*R12*R12DOTR13 + 2*R12*R12DOTR13PT - 
     1              2*DOT1*R13 - 2*R12DOTR13*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R13S - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R132)/
     1               R13S - 
     1              (GAMMA2*R12*DOT2*R132)/
     1               (R12S*R13S))*V313))/
     1       R123 + V312*
     1       ((2*R13*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V213)/R122 + 
     1         (-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12DOTR13*R13)/R12 + (2*R12DOTR13PT*R13)/R12 - 
     1            (GAMMA2*DOT2*R132)/R13S4 - 
     1            (3*GAMMA*DOT2*R132)/R13S3 - 
     1            (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (2*GAMMA*R12DOTR13PT*R132)/(R12*R13S) + 
     1            (GAMMA*DOT2*R13*SMALLA)/R13S3)*
     1          V313)))/(2.*R134)
                        HESS(K1+3,K1+3)=HESS(K1+3,K1+3)+ (ERL*((2*(R13*(-R12 + R12DOTR13*R13)*V213 + 
     1            R12*(R12*R12DOTR13 - R13)*V313)**2)/(R124*R134)
     1       + (4*GAMMA*R12DOTR13PT*
     1         (R13*(R12 - R12DOTR13*R13)*V213 + 
     1           R12*(-(R12*R12DOTR13) + R13)*V313)*
     1         (-(V213/(R12*R12S)) - 
     1           V313/(R13*R13S)))/(R122*R132) + 
     1      GAMMA2*DOT2*
     1       (V213/(R12*R12S) + 
     1          V313/(R13*R13S))**2 + 
     1      GAMMA*DOT2*
     1       (-(1/(R12*R12S)) - 1/(R13*R13S) + 
     1         ((3*R12 - SMALLA)*V2132)/
     1          (R123*R12S3) + 
     1         ((3*R13 - SMALLA)*V3132)/(R133*R13S3))
     1       - (2*R12DOTR13PT*(R133*(2*R12 - 3*R12DOTR13*R13)*
     1            V2132 + 
     1           2*R12*R13*(R122 - R12*R12DOTR13*R13 + R132)*V213*
     1            V313 + 
     1           R122*(R132*
     1               (R122*R12DOTR13 - 2*R12*R13 + R12DOTR13*R132) + 
     1              R12*(-3*R12*R12DOTR13 + 2*R13)*V3132)))/
     1       (R124*R134)))/2.
                        HESS(K1+3,K2+1)=HESS(K1+3,K2+1)+ (ERL*((2*R12*V311*
     1         (R13*(-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V213 + 
     1           R12*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V313))/R133 + 
     1      V211*((-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12*R12DOTR13)/R13 + (2*R12*R12DOTR13PT)/R13 - 
     1            (GAMMA2*R122*DOT2)/R12S4 - 
     1            (3*GAMMA*R122*DOT2)/R12S3 - 
     1            (4*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R122*R12DOTR13PT)/(R13*R12S) + 
     1            (GAMMA*R12*DOT2*SMALLA)/R12S3)*
     1          V213 + (R12*
     1            (-2*R12*DOT1 - 2*R12*R12DOTR13*R12DOTR13PT + 
     1              2*R12DOTR13*R13 + 2*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R122*R12DOTR13*R12DOTR13PT)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R12S - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R13S - 
     1              (GAMMA2*R122*DOT2*R13)/
     1               (R12S*R13S))*V313)/
     1          R132)))/(2.*R124)
                        HESS(K1+3,K2+2)=HESS(K1+3,K2+2)+ (ERL*((2*R12*V312*
     1         (R13*(-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V213 + 
     1           R12*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V313))/R133 + 
     1      V212*((-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12*R12DOTR13)/R13 + (2*R12*R12DOTR13PT)/R13 - 
     1            (GAMMA2*R122*DOT2)/R12S4 - 
     1            (3*GAMMA*R122*DOT2)/R12S3 - 
     1            (4*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R122*R12DOTR13PT)/(R13*R12S) + 
     1            (GAMMA*R12*DOT2*SMALLA)/R12S3)*
     1          V213 + (R12*
     1            (-2*R12*DOT1 - 2*R12*R12DOTR13*R12DOTR13PT + 
     1              2*R12DOTR13*R13 + 2*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R122*R12DOTR13*R12DOTR13PT)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R12S - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R13S - 
     1              (GAMMA2*R122*DOT2*R13)/
     1               (R12S*R13S))*V313)/
     1          R132)))/(2.*R124)
                        HESS(K1+3,K2+3)=HESS(K1+3,K2+3)+ (ERL*((-2*DOT1 + 
     1         R12DOTR13*(-6*R12DOTR13PT + (2*R12)/R13 - 
     1            (4*GAMMA*R12*R12DOTR13PT)/R12S) + 
     1         R12*R12DOTR13PT*(2/R13 - 
     1            (GAMMA2*R12*R12DOTR13PT)/R12S4 - 
     1            (3*GAMMA*R12*R12DOTR13PT)/R12S3 + 
     1            (2*GAMMA*R12)/(R13*R12S) + 
     1            (GAMMA*R12DOTR13PT*SMALLA)/R12S3))*
     1       V2132 + (R12*
     1         (4*(R12DOTR13 + R12DOTR13PT)*R13 + 
     1           R12*(-2 - 2*DOT1 - 2*R12DOTR13*R12DOTR13PT + 
     1              (4*GAMMA*R12DOTR13PT*R13)/R12S - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S) + 
     1           (GAMMA*R122*R12DOTR13PT*
     1              (-2*R12DOTR13 - (GAMMA*R12DOTR13PT*R13)/R13S)
     1              )/R12S)*V213*V313)/R132 + 
     1      R122*(R12DOTR13PT*
     1          (2*R12DOTR13 - (2*R12)/R13 + 
     1            (GAMMA*R122*R12DOTR13PT)/R12S3 - 
     1            (GAMMA*R12*R12DOTR13PT*SMALLA)/R12S3) + 
     1         (2*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V3132)/R133)))/(2.*R124)
                        HESS(K1+3,K3+1)=HESS(K1+3,K3+1)+ (ERL*((2*R13*V211*
     1         (R13*(-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V213 + 
     1           R12*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V313))/R123 + 
     1      V311*((R13*(2*R12*R12DOTR13 + 2*R12*R12DOTR13PT - 
     1              2*DOT1*R13 - 2*R12DOTR13*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R13S - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R132)/
     1               R13S - 
     1              (GAMMA2*R12*DOT2*R132)/
     1               (R12S*R13S))*V213)/
     1          R122 + (-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12DOTR13*R13)/R12 + (2*R12DOTR13PT*R13)/R12 - 
     1            (GAMMA2*DOT2*R132)/R13S4 - 
     1            (3*GAMMA*DOT2*R132)/R13S3 - 
     1            (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (2*GAMMA*R12DOTR13PT*R132)/(R12*R13S) + 
     1            (GAMMA*DOT2*R13*SMALLA)/R13S3)*
     1          V313)))/(2.*R134)
                        HESS(K1+3,K3+2)=HESS(K1+3,K3+2)+ (ERL*((2*R13*V212*
     1         (R13*(-R12 + R12DOTR13*R13 + R12DOTR13PT*R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1            V213 + 
     1           R12*(R12*R12DOTR13 + R12*R12DOTR13PT - R13 + 
     1              (GAMMA*R12*R12DOTR13PT*R13)/R13S)*
     1            V313))/R123 + 
     1      V312*((R13*(2*R12*R12DOTR13 + 2*R12*R12DOTR13PT - 
     1              2*DOT1*R13 - 2*R12DOTR13*R12DOTR13PT*R13 - 
     1              (2*GAMMA*R12*R12DOTR13*R12DOTR13PT*R13)/
     1               R12S + 
     1              (2*GAMMA*R12*R12DOTR13PT*R13)/R13S - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R132)/
     1               R13S - 
     1              (GAMMA2*R12*DOT2*R132)/
     1               (R12S*R13S))*V213)/
     1          R122 + (-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1            (2*R12DOTR13*R13)/R12 + (2*R12DOTR13PT*R13)/R12 - 
     1            (GAMMA2*DOT2*R132)/R13S4 - 
     1            (3*GAMMA*DOT2*R132)/R13S3 - 
     1            (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (2*GAMMA*R12DOTR13PT*R132)/(R12*R13S) + 
     1            (GAMMA*DOT2*R13*SMALLA)/R13S3)*
     1          V313)))/(2.*R134)
                        HESS(K1+3,K3+3)=HESS(K1+3,K3+3)+ (ERL*(R12DOTR13PT*R132*
     1       (2*R12DOTR13 - (2*R13)/R12 + 
     1         (GAMMA*R12DOTR13PT*R132)/R13S3 - 
     1         (GAMMA*R12DOTR13PT*R13*SMALLA)/R13S3) + 
     1      (2*R132*((R12DOTR13 + R12DOTR13PT)*R13 + 
     1           R12*(-1 + (GAMMA*R12DOTR13PT*R13)/R12S))*
     1         V2132)/R123 + 
     1      (R13*(2*R13*(-1 - DOT1 - R12DOTR13*R12DOTR13PT - 
     1              (GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S) + 
     1           R12*(4*R12DOTR13 + 4*R12DOTR13PT - 
     1              (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R12S + 
     1              (4*GAMMA*R12DOTR13PT*R13)/R13S - 
     1              (GAMMA2*DOT2*R132)/
     1               (R12S*R13S)))*V213*
     1         V313)/R122 + 
     1      (-2*DOT1 - 6*R12DOTR13*R12DOTR13PT + 
     1         (2*R12DOTR13*R13)/R12 + (2*R12DOTR13PT*R13)/R12 - 
     1         (GAMMA2*DOT2*R132)/R13S4 - 
     1         (3*GAMMA*DOT2*R132)/R13S3 - 
     1         (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1         (2*GAMMA*R12DOTR13PT*R132)/(R12*R13S) + 
     1         (GAMMA*DOT2*R13*SMALLA)/R13S3)*
     1       V3132))/(2.*R134)
                        HESS(K2+1,K2+1)=HESS(K2+1,K2+1)+ (ERL*((2*DOT1 + 
     1         R12DOTR13*(6*R12DOTR13PT + 
     1            (4*GAMMA*R12*R12DOTR13PT)/R12S) + 
     1         (GAMMA*R12*DOT2*
     1            (GAMMA*R12 + 3*R122 - 4*R12*SMALLA + SMALLA**2))/
     1          R12S4)*V2112 - 
     1      (4*R12*(R122*(R12DOTR13 + R12DOTR13PT) + 
     1           (R12DOTR13 + R12DOTR13PT)*SMALLA**2 + 
     1           R12*(GAMMA*R12DOTR13PT - 2*R12DOTR13*SMALLA - 
     1              2*R12DOTR13PT*SMALLA))*V211*V311)/
     1       (R13*R12S) + 
     1      R122*(R12DOTR13PT*
     1          (-2*R12DOTR13 - (GAMMA*R12*R12DOTR13PT)/R12S) + 
     1         (2*V3112)/R132)))/(2.*R124)
                        HESS(K2+1,K2+2)=HESS(K2+1,K2+2)+ (ERL*((2*R12*V311*
     1         (R13*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12*R12DOTR13PT)/R12S)*V212 + 
     1           R12*V312))/R132 + 
     1      V211*((2*DOT1 + 6*R12DOTR13*R12DOTR13PT + 
     1            (GAMMA2*R122*DOT2)/R12S4 + 
     1            (2*GAMMA*R122*DOT2)/R12S3 + 
     1            (4*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (GAMMA*R12*DOT2)/R12S)*V212
     1          + (2*R12*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12*R12DOTR13PT)/R12S)*V312)/
     1          R13)))/(2.*R124)
                        HESS(K2+1,K2+3)=HESS(K2+1,K2+3)+ (ERL*((2*R12*V311*
     1         (R13*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12*R12DOTR13PT)/R12S)*V213 + 
     1           R12*V313))/R132 + 
     1      V211*((2*DOT1 + 6*R12DOTR13*R12DOTR13PT + 
     1            (GAMMA2*R122*DOT2)/R12S4 + 
     1            (2*GAMMA*R122*DOT2)/R12S3 + 
     1            (4*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (GAMMA*R12*DOT2)/R12S)*V213
     1          + (2*R12*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12*R12DOTR13PT)/R12S)*V313)/
     1          R13)))/(2.*R124)
                        HESS(K2+1,K3+1)=HESS(K2+1,K3+1)+ (ERL*(2*R132*
     1       (-R12DOTR13 + R12DOTR13PT*(-1 - (GAMMA*R12)/R12S))*
     1       V2112 + R12*R13*
     1       (2 + 2*DOT1 + 
     1         2*R12DOTR13*R12DOTR13PT*
     1          (1 + (GAMMA*R12)/R12S + 
     1            (GAMMA*R13)/R13S) + 
     1         (GAMMA2*R12*DOT2*R13)/
     1          (R12S*R13S))*V211*
     1       V311 + 2*R122*
     1       (R12DOTR13PT*R132 + 
     1         (-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12DOTR13PT*R13)/R13S)*V3112))
     1    )/(2.*R123*R133)
                        HESS(K2+1,K3+2)=HESS(K2+1,K3+2)+ (ERL*(2*R12*V311*
     1       (R13*V212 + 
     1         R12*(-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12DOTR13PT*R13)/R13S)*V312) + 
     1      R13*V211*((-2*R12DOTR13*R13 - 2*R12DOTR13PT*R13 - 
     1            (2*GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1          V212 + R12*
     1          (2*DOT1 + 2*R12DOTR13*R12DOTR13PT + 
     1            (2*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (GAMMA2*R12*DOT2*R13)/
     1             (R12S*R13S))*V312)))/
     1  (2.*R123*R133)
                        HESS(K2+1,K3+3)=HESS(K2+1,K3+3)+ (ERL*(2*R12*V311*
     1       (R13*V213 + 
     1         R12*(-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12DOTR13PT*R13)/R13S)*V313) + 
     1      R13*V211*((-2*R12DOTR13*R13 - 2*R12DOTR13PT*R13 - 
     1            (2*GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1          V213 + R12*
     1          (2*DOT1 + 2*R12DOTR13*R12DOTR13PT + 
     1            (2*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (GAMMA2*R12*DOT2*R13)/
     1             (R12S*R13S))*V313)))/
     1  (2.*R123*R133)
                        HESS(K2+2,K2+2)=HESS(K2+2,K2+2)+ (ERL*((2*DOT1 + 
     1         R12DOTR13*(6*R12DOTR13PT + 
     1            (4*GAMMA*R12*R12DOTR13PT)/R12S) + 
     1         (GAMMA*R12*DOT2*
     1            (GAMMA*R12 + 3*R122 - 4*R12*SMALLA + SMALLA**2))/
     1          R12S4)*V2122 - 
     1      (4*R12*(R122*(R12DOTR13 + R12DOTR13PT) + 
     1           (R12DOTR13 + R12DOTR13PT)*SMALLA**2 + 
     1           R12*(GAMMA*R12DOTR13PT - 2*R12DOTR13*SMALLA - 
     1              2*R12DOTR13PT*SMALLA))*V212*V312)/
     1       (R13*R12S) + 
     1      R122*(R12DOTR13PT*
     1          (-2*R12DOTR13 - (GAMMA*R12*R12DOTR13PT)/R12S) + 
     1         (2*V3122)/R132)))/(2.*R124)
                        HESS(K2+2,K2+3)=HESS(K2+2,K2+3)+ (ERL*((2*R12*V312*
     1         (R13*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12*R12DOTR13PT)/R12S)*V213 + 
     1           R12*V313))/R132 + 
     1      V212*((2*DOT1 + 6*R12DOTR13*R12DOTR13PT + 
     1            (GAMMA2*R122*DOT2)/R12S4 + 
     1            (2*GAMMA*R122*DOT2)/R12S3 + 
     1            (4*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (GAMMA*R12*DOT2)/R12S)*V213
     1          + (2*R12*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12*R12DOTR13PT)/R12S)*V313)/
     1          R13)))/(2.*R124)
                        HESS(K2+2,K3+1)=HESS(K2+2,K3+1)+ (ERL*(2*R13*V211*
     1       (R13*(-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12*R12DOTR13PT)/R12S)*V212 + 
     1         R12*V312) + 
     1      R12*V311*(R13*
     1          (2*DOT1 + 2*R12DOTR13*R12DOTR13PT + 
     1            (2*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (GAMMA2*R12*DOT2*R13)/
     1             (R12S*R13S))*V212 + 
     1         2*R12*(-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12DOTR13PT*R13)/R13S)*V312)))/
     1  (2.*R123*R133)
                        HESS(K2+2,K3+2)=HESS(K2+2,K3+2)+ (ERL*(2*R132*
     1       (-R12DOTR13 + R12DOTR13PT*(-1 - (GAMMA*R12)/R12S))*
     1       V2122 + R12*R13*
     1       (2 + 2*DOT1 + 
     1         2*R12DOTR13*R12DOTR13PT*
     1          (1 + (GAMMA*R12)/R12S + 
     1            (GAMMA*R13)/R13S) + 
     1         (GAMMA2*R12*DOT2*R13)/
     1          (R12S*R13S))*V212*
     1       V312 + 2*R122*
     1       (R12DOTR13PT*R132 + 
     1         (-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12DOTR13PT*R13)/R13S)*V3122))
     1    )/(2.*R123*R133)
                        HESS(K2+2,K3+3)=HESS(K2+2,K3+3)+ (ERL*(2*R12*V312*
     1       (R13*V213 + 
     1         R12*(-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12DOTR13PT*R13)/R13S)*V313) + 
     1      R13*V212*((-2*R12DOTR13*R13 - 2*R12DOTR13PT*R13 - 
     1            (2*GAMMA*R12*R12DOTR13PT*R13)/R12S)*
     1          V213 + R12*
     1          (2*DOT1 + 2*R12DOTR13*R12DOTR13PT + 
     1            (2*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (GAMMA2*R12*DOT2*R13)/
     1             (R12S*R13S))*V313)))/
     1  (2.*R123*R133)
                        HESS(K2+3,K2+3)=HESS(K2+3,K2+3)+ (ERL*((2*DOT1 + 
     1         R12DOTR13*(6*R12DOTR13PT + 
     1            (4*GAMMA*R12*R12DOTR13PT)/R12S) + 
     1         (GAMMA*R12*DOT2*
     1            (GAMMA*R12 + 3*R122 - 4*R12*SMALLA + SMALLA**2))/
     1          R12S4)*V2132 - 
     1      (4*R12*(R122*(R12DOTR13 + R12DOTR13PT) + 
     1           (R12DOTR13 + R12DOTR13PT)*SMALLA**2 + 
     1           R12*(GAMMA*R12DOTR13PT - 2*R12DOTR13*SMALLA - 
     1              2*R12DOTR13PT*SMALLA))*V213*V313)/
     1       (R13*R12S) + 
     1      R122*(R12DOTR13PT*
     1          (-2*R12DOTR13 - (GAMMA*R12*R12DOTR13PT)/R12S) + 
     1         (2*V3132)/R132)))/(2.*R124)
                        HESS(K2+3,K3+1)=HESS(K2+3,K3+1)+ (ERL*(2*R13*V211*
     1       (R13*(-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12*R12DOTR13PT)/R12S)*V213 + 
     1         R12*V313) + 
     1      R12*V311*(R13*
     1          (2*DOT1 + 2*R12DOTR13*R12DOTR13PT + 
     1            (2*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (GAMMA2*R12*DOT2*R13)/
     1             (R12S*R13S))*V213 + 
     1         2*R12*(-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12DOTR13PT*R13)/R13S)*V313)))/
     1  (2.*R123*R133)
                        HESS(K2+3,K3+2)=HESS(K2+3,K3+2)+ (ERL*(2*R13*V212*
     1       (R13*(-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12*R12DOTR13PT)/R12S)*V213 + 
     1         R12*V313) + 
     1      R12*V312*(R13*
     1          (2*DOT1 + 2*R12DOTR13*R12DOTR13PT + 
     1            (2*GAMMA*R12*R12DOTR13*R12DOTR13PT)/R12S + 
     1            (2*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (GAMMA2*R12*DOT2*R13)/
     1             (R12S*R13S))*V213 + 
     1         2*R12*(-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12DOTR13PT*R13)/R13S)*V313)))/
     1  (2.*R123*R133)
                        HESS(K2+3,K3+3)=HESS(K2+3,K3+3)+ (ERL*(2*R132*
     1       (-R12DOTR13 + R12DOTR13PT*(-1 - (GAMMA*R12)/R12S))*
     1       V2132 + R12*R13*
     1       (2 + 2*DOT1 + 
     1         2*R12DOTR13*R12DOTR13PT*
     1          (1 + (GAMMA*R12)/R12S + 
     1            (GAMMA*R13)/R13S) + 
     1         (GAMMA2*R12*DOT2*R13)/
     1          (R12S*R13S))*V213*
     1       V313 + 2*R122*
     1       (R12DOTR13PT*R132 + 
     1         (-R12DOTR13 - R12DOTR13PT - 
     1            (GAMMA*R12DOTR13PT*R13)/R13S)*V3132))
     1    )/(2.*R123*R133)
                        HESS(K3+1,K3+1)=HESS(K3+1,K3+1)+ (ERL*(R12DOTR13PT*R132*
     1       (-2*R12DOTR13 - (GAMMA*R12DOTR13PT*R13)/R13S) + 
     1      (2*R132*V2112)/R122 - 
     1      (4*R13*(GAMMA*R12DOTR13PT*R13 + 
     1           (R12DOTR13 + R12DOTR13PT)*R13S)*V211*
     1         V311)/(R12*R13S) + 
     1      (2*DOT1 + 6*R12DOTR13*R12DOTR13PT + 
     1         (GAMMA2*DOT2*R132)/R13S4 + 
     1         (2*GAMMA*DOT2*R132)/R13S3 + 
     1         (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1         (GAMMA*DOT2*R13)/R13S)*V3112))
     1   /(2.*R134)
                        HESS(K3+1,K3+2)=HESS(K3+1,K3+2)+ (ERL*((2*R13*V211*
     1         (R13*V212 + 
     1           R12*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12DOTR13PT*R13)/R13S)*V312))/
     1       R122 + V311*
     1       ((2*R13*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12DOTR13PT*R13)/R13S)*V212)/
     1          R12 + (2*DOT1 + 6*R12DOTR13*R12DOTR13PT + 
     1            (GAMMA2*DOT2*R132)/R13S4 + 
     1            (2*GAMMA*DOT2*R132)/R13S3 + 
     1            (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (GAMMA*DOT2*R13)/R13S)*V312))
     1    )/(2.*R134)
                        HESS(K3+1,K3+3)=HESS(K3+1,K3+3)+ (ERL*((2*R13*V211*
     1         (R13*V213 + 
     1           R12*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12DOTR13PT*R13)/R13S)*V313))/
     1       R122 + V311*
     1       ((2*R13*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12DOTR13PT*R13)/R13S)*V213)/
     1          R12 + (2*DOT1 + 6*R12DOTR13*R12DOTR13PT + 
     1            (GAMMA2*DOT2*R132)/R13S4 + 
     1            (2*GAMMA*DOT2*R132)/R13S3 + 
     1            (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (GAMMA*DOT2*R13)/R13S)*V313))
     1    )/(2.*R134)
                        HESS(K3+2,K3+2)=HESS(K3+2,K3+2)+ (ERL*(R12DOTR13PT*R132*
     1       (-2*R12DOTR13 - (GAMMA*R12DOTR13PT*R13)/R13S) + 
     1      (2*R132*V2122)/R122 - 
     1      (4*R13*(GAMMA*R12DOTR13PT*R13 + 
     1           (R12DOTR13 + R12DOTR13PT)*R13S)*V212*
     1         V312)/(R12*R13S) + 
     1      (2*DOT1 + 6*R12DOTR13*R12DOTR13PT + 
     1         (GAMMA2*DOT2*R132)/R13S4 + 
     1         (2*GAMMA*DOT2*R132)/R13S3 + 
     1         (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1         (GAMMA*DOT2*R13)/R13S)*V3122))
     1   /(2.*R134)
                        HESS(K3+2,K3+3)=HESS(K3+2,K3+3)+ (ERL*((2*R13*V212*
     1         (R13*V213 + 
     1           R12*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12DOTR13PT*R13)/R13S)*V313))/
     1       R122 + V312*
     1       ((2*R13*(-R12DOTR13 - R12DOTR13PT - 
     1              (GAMMA*R12DOTR13PT*R13)/R13S)*V213)/
     1          R12 + (2*DOT1 + 6*R12DOTR13*R12DOTR13PT + 
     1            (GAMMA2*DOT2*R132)/R13S4 + 
     1            (2*GAMMA*DOT2*R132)/R13S3 + 
     1            (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1            (GAMMA*DOT2*R13)/R13S)*V313))
     1    )/(2.*R134)
                        HESS(K3+3,K3+3)=HESS(K3+3,K3+3)+ (ERL*(R12DOTR13PT*R132*
     1       (-2*R12DOTR13 - (GAMMA*R12DOTR13PT*R13)/R13S) + 
     1      (2*R132*V2132)/R122 - 
     1      (4*R13*(GAMMA*R12DOTR13PT*R13 + 
     1           (R12DOTR13 + R12DOTR13PT)*R13S)*V213*
     1         V313)/(R12*R13S) + 
     1      (2*DOT1 + 6*R12DOTR13*R12DOTR13PT + 
     1         (GAMMA2*DOT2*R132)/R13S4 + 
     1         (2*GAMMA*DOT2*R132)/R13S3 + 
     1         (4*GAMMA*R12DOTR13*R12DOTR13PT*R13)/R13S + 
     1         (GAMMA*DOT2*R13)/R13S)*V3132))
     1   /(2.*R134)

                        HESS(K1+2,K1+1)=HESS(K1+1,K1+2) 
                        HESS(K1+3,K1+1)=HESS(K1+1,K1+3) 
                        HESS(K2+1,K1+1)=HESS(K1+1,K2+1)
                        HESS(K2+2,K1+1)=HESS(K1+1,K2+2)
                        HESS(K2+3,K1+1)=HESS(K1+1,K2+3)
                        HESS(K3+1,K1+1)=HESS(K1+1,K3+1)
                        HESS(K3+2,K1+1)=HESS(K1+1,K3+2)
                        HESS(K3+3,K1+1)=HESS(K1+1,K3+3)

                        HESS(K1+3,K1+2)=HESS(K1+2,K1+3)
                        HESS(K2+1,K1+2)=HESS(K1+2,K2+1)
                        HESS(K2+2,K1+2)=HESS(K1+2,K2+2)
                        HESS(K2+3,K1+2)=HESS(K1+2,K2+3)
                        HESS(K3+1,K1+2)=HESS(K1+2,K3+1)
                        HESS(K3+2,K1+2)=HESS(K1+2,K3+2)
                        HESS(K3+3,K1+2)=HESS(K1+2,K3+3)

                        HESS(K2+1,K1+3)=HESS(K1+3,K2+1)
                        HESS(K2+2,K1+3)=HESS(K1+3,K2+2)
                        HESS(K2+3,K1+3)=HESS(K1+3,K2+3)
                        HESS(K3+1,K1+3)=HESS(K1+3,K3+1)
                        HESS(K3+2,K1+3)=HESS(K1+3,K3+2)
                        HESS(K3+3,K1+3)=HESS(K1+3,K3+3)

                        HESS(K2+2,K2+1)=HESS(K2+1,K2+2)
                        HESS(K2+3,K2+1)=HESS(K2+1,K2+3)
                        HESS(K3+1,K2+1)=HESS(K2+1,K3+1)
                        HESS(K3+2,K2+1)=HESS(K2+1,K3+2)
                        HESS(K3+3,K2+1)=HESS(K2+1,K3+3)

                        HESS(K2+3,K2+2)=HESS(K2+2,K2+3)
                        HESS(K3+1,K2+2)=HESS(K2+2,K3+1)
                        HESS(K3+2,K2+2)=HESS(K2+2,K3+2)
                        HESS(K3+3,K2+2)=HESS(K2+2,K3+3)

                        HESS(K3+1,K2+3)=HESS(K2+3,K3+1)
                        HESS(K3+2,K2+3)=HESS(K2+3,K3+2)
                        HESS(K3+3,K2+3)=HESS(K2+3,K3+3)

                        HESS(K3+2,K3+1)=HESS(K3+1,K3+2)
                        HESS(K3+3,K3+1)=HESS(K3+1,K3+3)

                        HESS(K3+3,K3+2)=HESS(K3+2,K3+3)

                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      ENERGY=ENERGY*LAMBDA/2.0D0

      RETURN
      END
