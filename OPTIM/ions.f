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
C  SHIELDED BORN-MAYER POTENTIAL. IONS MUST BE ENERTED ALTERNATING
C  PLUS, MINUS, PLUS AND THERE MUST BE EQUAL NUMBERS OF EACH SIGN.
C
C*************************************************************************
C
      SUBROUTINE IONS(N, X, V, ENERGY, GAMMA, CHARGE, RHO, APP, AMM, APM, NTEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NTEST
!
! SHOULD REALLY MAKE THE 2D ARRAYS ALLOCATBLE
!
      DOUBLE PRECISION X(3*N), ENERGY, GAMMA, CHARGE,
     1                 V(3*N), R(N,N), EXG(N,N), RHO, APP, 
     2                 R3(N,N), G(N,N), AMM, APM, RIJ, 
     3                 R2(N,N), EXR(N,N), F(N,N)

      IF (CHARGE.EQ.0.0D0) CHARGE=1.0D0
      IF (RHO.EQ.0.0D0) RHO=0.636848D0
      IF (APP.EQ.0.0D0) APP=57.153D0
      IF (AMM.EQ.0.0D0) AMM=70.735D0
      IF (APM.EQ.0.0D0) APM=65.667D0
C 
C  STORE DISTANCE MATRICES.
C
      DO 20 J1=1,N
         R(J1,J1)=0.0D0
         R2(J1,J1)=0.0D0
         R3(J1,J1)=0.0D0
         EXR(J1,J1)=1.0D0
         EXG(J1,J1)=1.0D0
         DO 10 J2=J1+1,N
            R(J2,J1)=DSQRT( (X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                     +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                     +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2 )
            EXR(J2,J1)=DEXP(-R(J2,J1)/RHO)
            EXG(J2,J1)=DEXP(-R(J2,J1)*GAMMA)
            R2(J2,J1)=1.0D0/R(J2,J1)**2
            R3(J2,J1)=R2(J2,J1)/R(J2,J1)
            EXR(J1,J2)=EXR(J2,J1)
            EXG(J1,J2)=EXG(J2,J1)
            R(J1,J2)=R(J2,J1)
            R2(J1,J2)=R2(J2,J1)
            R3(J1,J2)=R3(J2,J1)
10       CONTINUE
20    CONTINUE
C
C  CALCULATE THE G AND F TENSORS AND THE ENERGY.
C
C  PLUS - PLUS
C
      ENERGY=0.0D0

      DO J1=1,N-1,2
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO J2=J1+2,N-1,2 
            RIJ=R(J2,J1)
            G(J2,J1)=-APP*EXR(J2,J1)/(RHO*RIJ) - CHARGE**2*(1.0D0 + 
     1                GAMMA*RIJ)*EXG(J2,J1)*R3(J2,J1)
            F(J2,J1)=APP*EXR(J2,J1)/RHO**2 + CHARGE**2*(2.0D0 + 
     1                2.0D0*GAMMA*RIJ + (GAMMA*RIJ)**2)*EXG(J2,J1)
     2               *R3(J2,J1)-G(J2,J1)            
            G(J1,J2)=G(J2,J1)
            F(J1,J2)=F(J2,J1)
            ENERGY=ENERGY+CHARGE**2*EXG(J2,J1)/RIJ + APP*EXR(J2,J1)
         ENDDO
      ENDDO
C
C  MINUS - MINUS
C
      DO J1=2,N,2
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO J2=J1+2,N,2 
            RIJ=R(J2,J1)
            G(J2,J1)=-(AMM*EXR(J2,J1)/(RHO*RIJ)) - CHARGE**2*(1.0D0 
     1               + GAMMA*RIJ)*EXG(J2,J1)*R3(J2,J1)
            F(J2,J1)=AMM*EXR(J2,J1)/RHO**2 + CHARGE**2*(2.0D0 
     1 + 2.0D0*GAMMA*RIJ + (GAMMA*RIJ)**2)*EXG(J2,J1)*R3(J2,J1)-G(J2,J1)
            G(J1,J2)=G(J2,J1)
            F(J1,J2)=F(J2,J1)
            ENERGY=ENERGY+CHARGE**2*EXG(J2,J1)/RIJ + AMM*EXR(J2,J1)
         ENDDO
      ENDDO
C
C  PLUS - MINUS
C
      DO J1=1,N-1,2
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO J2=2,N,2 
            RIJ=R(J2,J1)
            G(J2,J1)=-(APM*EXR(J2,J1)/(RHO*RIJ)) + CHARGE**2*(1.0D0 
     1               + GAMMA*RIJ)*EXG(J2,J1)*R3(J2,J1)
            F(J2,J1)=APM*EXR(J2,J1)/RHO**2 - CHARGE**2*(2.0D0 + 
     1   2.0D0*GAMMA*RIJ + (GAMMA*RIJ)**2)*EXG(J2,J1)*R3(J2,J1)-G(J2,J1)
            G(J1,J2)=G(J2,J1)
            F(J1,J2)=F(J2,J1)
            ENERGY=ENERGY-CHARGE**2*EXG(J2,J1)/RIJ + APM*EXR(J2,J1)
         ENDDO
      ENDDO

      IF (NTEST.EQ.0) RETURN
C
C  FROM HERE ON DOWN THE CODE IS SYSTEM-INDEPENDENT!
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
               V(J3)=V(J3)+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
30          CONTINUE
C           PRINT*,'J3,V(J3)=',J3,V(J3)
40       CONTINUE
50    CONTINUE
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO 60 J4=1,N
               HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
60          CONTINUE
70       CONTINUE
80    CONTINUE
C
C  NEXT ARE THE TERMS WHERE X_I AND X_J ARE ON THE SAME ATOM
C  BUT ARE DIFFERENT, E.G. Y AND Z.
C
      DO 120 J1=1,N
         DO 110 J2=1,3
            J3=3*(J1-1)+J2
            DO 100 J4=J2+1,3
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO 90 J5=1,N
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3) + 
     1            F(J5,J1)*R2(J5,J1)* 
     2           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
90             CONTINUE
100         CONTINUE
110      CONTINUE
120   CONTINUE
C
C  CASE III, DIFFERENT ATOMS, SAME CARTESIAN COORDINATE.
C
      DO 150 J1=1,N
         DO 140 J2=1,3
            J3=3*(J1-1)+J2
            DO 130 J4=J1+1,N
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
130         CONTINUE
140      CONTINUE
150   CONTINUE
C
C  CASE IV: DIFFERENT ATOMS AND DIFFERENT CARTESIAN COORDINATES.
C
      DO 180 J1=1,N
         DO 170 J2=1,3
            J3=3*(J1-1)+J2
            DO 160 J4=J1+1,N
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
155            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
C
C  SYMMETRISE HESSIAN
C
      DO 200 J1=1,3*N
         DO 190 J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
190      CONTINUE
200   CONTINUE
      RETURN
      END
