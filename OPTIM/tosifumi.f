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
C  ENERGY AND DERIVATIVES OF A GENERAL BORN-MAYER POTENTIAL WITH
C  TOSI-FUMI PARAMETERS SPECIFIED AFTER THE KEYWORD TOSI IN THE
C  ODATA FILE. ENERGY IS IN HARTREE, LENGTH IN BOHR.
C
      SUBROUTINE TOSIFUMI(N, X, V, POTEL, APP, AMM, APM, RHO, ZSYM, GTEST, STEST)
      USE PORFUNCS
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, I, J, J1, J2, K, L, J3
      DOUBLE PRECISION X(3*N), V(3*N), POTEL, 
     1                 APP, AMM, APM, RHO, DUM
      DOUBLE PRECISION AC(N,N),Q(N),RD(N,N)
      CHARACTER(LEN=5) ZSYM(N)
      LOGICAL GTEST, STEST

      DO I=1,N
         IF (ZSYM(I).EQ.'MI') THEN
            Q(I)=-1.0D0
         ELSE IF (ZSYM(I).EQ.'PL') THEN
            Q(I)=1.0D0
         ELSE
            WRITE(*,'(A)') ' ALL ATOMS MUST BE TYPE PL OR MI FOR TOSI-FUMI'
            STOP
         ENDIF
      ENDDO

       DO I=1,N
          DO J=I,N
             IF (ZSYM(I).EQ.'PL') THEN 
                IF (ZSYM(J).EQ.'PL') THEN 
                   AC(I,J)=APP
                   AC(J,I)=APP
                ELSE IF (ZSYM(J).EQ.'MI') THEN 
                   AC(I,J)=APM
                   AC(J,I)=APM
                ENDIF
             ELSE IF (ZSYM(I).EQ.'MI') THEN
                IF (ZSYM(J).EQ.'MI') THEN 
                   AC(I,J)=AMM
                   AC(J,I)=AMM
                ELSE IF (ZSYM(J).EQ.'PL') THEN 
                   AC(I,J)=APM
                   AC(J,I)=APM
                ENDIF
             ENDIF 
         ENDDO
      ENDDO
C
C  CALCULATE THE INTERPARTICLE DISTANCES.
C 
      DO I=1,N
         K=1+3*(I-1)         
         DO J=I+1,N
            L=K+3*(J-I)
            RD(I,J)= DSQRT((X(K)-X(L))*(X(K)-X(L))+
     1                     (X(K+1)-X(L+1))*(X(K+1)-X(L+1))+
     2                     (X(K+2)-X(L+2))*(X(K+2)-X(L+2))) 
            RD(J,I)=RD(I,J)
         ENDDO
      ENDDO

      POTEL=0.0D0
      DO I=1,N
         DO J=I+1,N
            POTEL=POTEL+Q(I)*Q(J)/RD(J,I)+AC(J,I)*DEXP(-RD(J,I)/RHO)
         ENDDO
      ENDDO

      IF (.NOT.GTEST) RETURN
C
C  GRADIENT.
C
      DO J1=1,N
         DO J2=1,3
            DUM=0.0D0
            DO J3=1,N
               IF (J1.NE.J3) THEN
                  DUM=DUM-(Q(J3)*Q(J1)/RD(J3,J1)**3+AC(J3,J1)*DEXP(-RD(J3,J1)/RHO)/(RD(J3,J1)*RHO))
     1                       *(X(3*(J1-1)+J2)-X(3*(J3-1)+J2))
               ENDIF
            ENDDO
            V(3*(J1-1)+J2)=DUM
         ENDDO
      ENDDO

      IF (.NOT.STEST) RETURN
C
C  HESSIAN.
C
      DO I=1,3*N
         DO J=1,3*N
            HESS(I,J)=0.0D0
         ENDDO
      ENDDO
      DO I=1,3*N
         DO J=I,3*N
C 
C  DETERMINE THE CARTESIAN COORDINATE OF X(I) AND X(J).
C
            K=MOD((I-1),3)
            L=MOD((J-1),3)
C 
C  THIS IF STATEMENT CALLS SHESS IF I AND J ARE THE SAME
C  CARTESIAN COORDINATE AND CALLS DHESS IF THE ARE DIFFERENT
C  CARTESIAN COORDINATES
C
            IF(K .EQ. L)THEN
               CALL NSHESS(I,J,K,L,X,N,RHO,AC,Q,RD)
            ELSE
               CALL NDHESS(I,J,K,L,X,N,RHO,AC,Q,RD)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
C....................................................
C.....THIS SUBROUTINE IF THE COORDINATES ARE THE SAME
 
      SUBROUTINE NSHESS(I,J,K,L,CFG,NUM,RHO,AC,Q,RD)
      USE MODHESS
      IMPLICIT NONE
      INTEGER I,J,K,L,I1,J1,M,N,NUM
      DOUBLE PRECISION CFG(3*NUM), RW, RHO, R, QQ, AB, 
     1                 EP, DR, DF
      DOUBLE PRECISION AC(NUM,NUM),Q(NUM),RD(NUM,NUM)

      RW=RHO
 
C.....I1 AND J1 ARE THE PARTICLE NUMBERS FOR X(I) AND X(J)

      I1 = ((I-1)/3) + 1
      J1 = ((J-1)/3) + 1
 
C.....IF THE PARTICLES ARE THE SAME DO THE NEXT TWO LOOPS
C.....THERE ARE TWO LOOPS SO THAT WE CAN SKIP PARTICLE I1
C.....WITHOUT RESORTING TO AN IF STATEMENT

      IF(I1 .EQ. J1)THEN
      
      DO 10 M = 1, I1-1
        
      N = 3*(M-1) + K + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N))*(CFG(I)-CFG(N))
      
      DF = -(QQ/(R*R*R)) + (3.0D0*QQ*DR/(R*R*R*R*R))
     1     + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
     2     - (AB*EP/(RW*R))
 
      HESS(I,I) = HESS(I,I) + DF
      
   10 CONTINUE
 
      DO 20 M = I1+1, NUM
         
      N = 3*(M-1) + K + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N))*(CFG(I)-CFG(N))
      
      DF = -(QQ/(R*R*R)) + (3.0D0*QQ*DR/(R*R*R*R*R))
     1     + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
     2     - (AB*EP/(RW*R))
 
      HESS(I,I) = HESS(I,I) + DF
      
   20 CONTINUE
 
      ELSE
 
C.....DO THIS IF THE PARTICLES ARE NOT THE SAME
      R = RD(I1,J1)
      QQ = Q(I1)*Q(J1)
      AB = AC(I1,J1)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(J))*(CFG(I)-CFG(J))
      
      HESS(I,J) = (QQ/(R*R*R)) - (3.0D0*QQ*DR/(R*R*R*R*R))
     1     - (AB*DR*EP/(RW*R*R*R)) - (AB*DR*EP/(RW*RW*R*R))
     2     + (AB*EP/(RW*R))
 
      HESS(J,I) = HESS(I,J)
 
      ENDIF
 
      RETURN
      END
C..........................................................
C.....THIS SUBROUTINE IF THE COORDINATES ARE NOT THE SAME
 
      SUBROUTINE NDHESS(I,J,K,L,CFG,NUM,RHO,AC,Q,RD)
      USE MODHESS
      IMPLICIT NONE
      INTEGER I,J,K,L,I1,J1,M,N,NUM,N1
      DOUBLE PRECISION CFG(3*NUM),RHO,RW,R,QQ,AB,
     1                 EP,DR,DF
      DOUBLE PRECISION AC(NUM,NUM),Q(NUM),RD(NUM,NUM)
      RW=RHO
      
      I1 = ((I-1)/3) + 1
      J1 = ((J-1)/3) + 1
 
C.....DO THE NEXT TO LOOPS IF THE PARTICLES ARE THE SAME

      IF(I1 .EQ. J1)THEN
 
      DO 10 M = 1, I1-1
 
      N = 3*(M-1) + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N+K))*(CFG(J)-CFG(N+L))
      
      DF = (3.0D0*QQ*DR/(R*R*R*R*R))
     1    + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
 
      HESS(I,J) = HESS(I,J) + DF
      
   10 CONTINUE
 
      DO 20 M = I1+1 ,NUM
         
      N = 3*(M-1) + 1
 
      R = RD(I1,M)
      QQ = Q(I1)*Q(M)
      AB = AC(I1,M)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N+K))*(CFG(J)-CFG(N+L))
      
      DF = (3.0D0*QQ*DR/(R*R*R*R*R))
     1    + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
 
      HESS(I,J) = HESS(I,J) + DF
      
   20 CONTINUE
      HESS(J,I) = HESS(I,J)
 
      ELSE
 
C.....DO THIS IF THE PARTICLES ARE NOT THE SAME
      N = K + 3*(J1-1) + 1
      N1 = L + 3*(I1-1) + 1
 
      R = RD(I1,J1)
      QQ = Q(I1)*Q(J1)
      AB = AC(I1,J1)
      EP = DEXP(-R/RW)
      DR = (CFG(I)-CFG(N))*(CFG(N1)-CFG(J))
      
      HESS(I,J) = -(3.0D0*QQ*DR/(R*R*R*R*R))
     1         - (AB*DR*EP/(RW*R*R*R)) - (AB*DR*EP/(RW*RW*R*R))
 
      HESS(J,I) = HESS(I,J)
 
      ENDIF
 
      RETURN
      END
C
C  ENERGY AND ANALYTIC DERIVATIVES OF POSSIBLE DISPERSION TERMS.
C
      SUBROUTINE TOSIFUMIC6(N, X, V, EDISP, C6PP, C6MM, C6PM, ZSYM, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J6, J4, J5
      DOUBLE PRECISION X(3*N), R6,
     1                 V(3*N), R2(N,N), DUMMY3, DUMMY4,
     2                 R8(N,N), G(N,N), 
     3                 F(N,N), DUMMY, EDISP, C6PP, C6MM, C6PM
      CHARACTER(LEN=5) ZSYM(N)
      LOGICAL GTEST, STEST
C 
C  STORE DISTANCE MATRICES.
C
      EDISP=0.0D0
      DO J1=1,N
         R2(J1,J1)=0.0D0
         R8(J1,J1)=0.0D0
         DO J2=J1+1,N
            R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1               +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2               +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R2(J2,J1)=1.0D0/R2(J2,J1)
            R6=R2(J2,J1)**3
            IF (ZSYM(J2).EQ.'PL') THEN
               IF (ZSYM(J1).EQ.'PL') THEN
                  EDISP=EDISP-C6PP*R6
               ELSE
                  EDISP=EDISP-C6PM*R6
               ENDIF
            ELSE
               IF (ZSYM(J1).EQ.'PL') THEN
                  EDISP=EDISP-C6PM*R6
               ELSE
                  EDISP=EDISP-C6MM*R6
               ENDIF
            ENDIF
            R8(J2,J1)=R2(J2,J1)**4
            R2(J1,J2)=R2(J2,J1)
         ENDDO
      ENDDO

      IF (.NOT.GTEST) RETURN
C
C  CALCULATE THE G TENSOR.
C
      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            IF (ZSYM(J2).EQ.'PL') THEN
               IF (ZSYM(J1).EQ.'PL') THEN
                  G(J2,J1)=6.0D0*C6PP*R8(J2,J1)
               ELSE
                  G(J2,J1)=6.0D0*C6PM*R8(J2,J1)
               ENDIF
            ELSE
               IF (ZSYM(J1).EQ.'PL') THEN
                  G(J2,J1)=6.0D0*C6PM*R8(J2,J1)
               ELSE
                  G(J2,J1)=6.0D0*C6MM*R8(J2,J1)
               ENDIF
            ENDIF
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
            ENDDO
            V(J3)=V(J3)+DUMMY
         ENDDO
      ENDDO

      IF (.NOT.STEST) RETURN
C
C  F TENSOR
C
      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            IF (ZSYM(J2).EQ.'PL') THEN
               IF (ZSYM(J1).EQ.'PL') THEN
                  F(J2,J1)=-42.0D0*C6PP*R8(J2,J1)-G(J2,J1)
               ELSE
                  F(J2,J1)=-42.0D0*C6PM*R8(J2,J1)-G(J2,J1)
               ENDIF
            ELSE
               IF (ZSYM(J1).EQ.'PL') THEN
                  F(J2,J1)=-42.0D0*C6PM*R8(J2,J1)-G(J2,J1)
               ELSE
                  F(J2,J1)=-42.0D0*C6MM*R8(J2,J1)-G(J2,J1)
               ENDIF
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
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=HESS(J3,J3)+DUMMY
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
                  DUMMY3=F(J5,J1)*R2(J5,J1)
                  DUMMY4=(X(J3)-X(3*J5-3+J2))*(X(3*J1-3+J4)-X(3*J5-3+J4))
                  DUMMY=DUMMY + DUMMY3*DUMMY4
               ENDDO
               HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+DUMMY
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
               HESS(3*(J4-1)+J2,J3)=HESS(3*(J4-1)+J2,J3)-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
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
                  HESS(J6,J3)=HESS(J6,J3)-F(J4,J1)*R2(J4,J1)
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
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END

C**************************************************************************
C
C  ENERGY AND ANALYTIC DERIVATIVES OF POSSIBLE FIRST ORDER INDUCTION.
C
      SUBROUTINE TOSIFUMIPOL(N, X, V, EIND, ALPHAP, ALPHAM, ZSYM, DAMP, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2
      DOUBLE PRECISION X(3*N), EIND, DAMP,
     1                 V(3*N), 
     2                 ALPHAP, ALPHAM, VEC1(3*N),
     3                 VEC2(3*N), DIF, TEMP1
      CHARACTER(LEN=5) ZSYM(N)
      LOGICAL GTEST, STEST

      CALL TFIND(N, X, EIND, ALPHAP, ALPHAM, ZSYM, DAMP)

      IF (.NOT.GTEST) RETURN
C
C  ANALYTIC GRADIENT.
C
      CALL TFGRAD(N, X, VEC1, ALPHAP, ALPHAM, ZSYM, DAMP)
      DO J1=1,3*N
         V(J1)=V(J1)+VEC1(J1)
      ENDDO

      IF (.NOT.STEST) RETURN
C
C  ANALYTIC SECOND DERIVATIVES. NOT DONE FOR DAMPING!
C
C     PRINT*,'ANALYTIC SECOND DERIVATIVES:'
C     CALL TFSEC(N, X, ALPHAP, ALPHAM, ZSYM, DAMP, NPROD)
C     DO J1=1,3*N
C        DO J2=1,3*N
C           WRITE(*,'(A,2I4,F20.10)') 'J1,J2,A=',J1,J2,HESS(J2,J1)
C        ENDDO
C     ENDDO
C
C  NUMERICAL SECOND DERIVATIVES.
C
C     PRINT*,'NUMERICAL SECOND DERIVATIVES:'
      DIF=1.0D-4
      DO J1=1,3*N
         TEMP1=X(J1)
         X(J1)=X(J1)+DIF
         CALL TFGRAD(N, X, VEC1, ALPHAP, ALPHAM, ZSYM, DAMP)
         X(J1)=X(J1)-2.0D0*DIF
         CALL TFGRAD(N, X, VEC2, ALPHAP, ALPHAM, ZSYM, DAMP)
         DO J2=1,3*N
            HESS(J2,J1)=HESS(J2,J1)+(VEC1(J2)-VEC2(J2))/(2.0D0*DIF)
C           WRITE(*,'(A,2I3,F20.10)') 'J1,J2,A=',J1,J2,(VEC1(J2)-VEC2(J2))/(2.0D0*DIF)
         ENDDO
         X(J1)=TEMP1
      ENDDO

      RETURN
      END

C*******************************************************************************
C
C  FIRST ORDER INDUCTION ENERGY
C
      SUBROUTINE TFIND(N, X, EIND, ALPHAP, ALPHAM, ZSYM, DAMP)
      IMPLICIT NONE
      INTEGER N, J1, J2, NPROD(N,N)
      DOUBLE PRECISION X(3*N), EIND, DAMP, DUMMY3,
     1                 R(N,N), 
     2                 DUMMY, ALPHAP, ALPHAM, VEC(3), 
     3                 F2, RR(N,N)
      CHARACTER(LEN=5) ZSYM(N)
C
C STATEMENT FUNCTIONS
C
      F2(DUMMY)=1.0D0-DEXP(-DUMMY/DAMP)*(1.0D0+DUMMY/DAMP+DUMMY**2/(2.0D0*DAMP**2))
C     F2(DUMMY)=1.0D0
C 
C  STORE DISTANCE MATRICES.
C
      DO J1=1,N
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         NPROD(J1,J1)=1
         DO J2=J1+1,N
            R(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1              +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2              +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R(J2,J1)=1.0D0/DSQRT(R(J2,J1))
            RR(J2,J1)=1.0D0/R(J2,J1)
            R(J1,J2)=R(J2,J1)
            RR(J1,J2)=RR(J2,J1)
            IF (ZSYM(J2).EQ.ZSYM(J1)) THEN
               NPROD(J2,J1)=1
               NPROD(J1,J2)=1
            ELSE
               NPROD(J2,J1)=-1
               NPROD(J1,J2)=-1
            ENDIF
         ENDDO
      ENDDO
C
C  FIRST ORDER INDUCTION ENERGY.
C
      EIND=0.0D0
      DO J1=1,N
         DUMMY=0.0D0
         VEC(1)=0.0D0
         VEC(2)=0.0D0
         VEC(3)=0.0D0
         DO J2=1,N
            DUMMY3=F2(RR(J2,J1))*R(J2,J1)**3
            IF (ZSYM(J2).EQ.'PL') THEN
               VEC(1)=VEC(1)+(X(3*J1-2)-X(3*J2-2))*DUMMY3
               VEC(2)=VEC(2)+(X(3*J1-1)-X(3*J2-1))*DUMMY3
               VEC(3)=VEC(3)+(X(3*J1)-X(3*J2))*DUMMY3
            ELSE
               VEC(1)=VEC(1)-(X(3*J1-2)-X(3*J2-2))*DUMMY3
               VEC(2)=VEC(2)-(X(3*J1-1)-X(3*J2-1))*DUMMY3
               VEC(3)=VEC(3)-(X(3*J1)-X(3*J2))*DUMMY3
            ENDIF
         ENDDO
         IF (ZSYM(J1).EQ.'PL') THEN
            EIND=EIND-ALPHAP*(VEC(1)**2+VEC(2)**2+VEC(3)**2)/2.0D0
         ELSE
            EIND=EIND-ALPHAM*(VEC(1)**2+VEC(2)**2+VEC(3)**2)/2.0D0
         ENDIF
      ENDDO

      RETURN
      END

C*************************************************************************
C
C  ANALYTIC FIRST DERIVATIVES OF POSSIBLE FIRST ORDER INDUCTION.
C
      SUBROUTINE TFGRAD(N, X, V, ALPHAP, ALPHAM, ZSYM, DAMP)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, NPROD(N,N)
      DOUBLE PRECISION X(3*N), DAMP,
     1                 V(3*N), R(N,N), DUMMY3, DUMMY4,
     2                 DUMMY, ALPHAP, ALPHAM, DUMMY5, DUMMY6,
     3                 DOT(N,N), R2(N,N), DUMMY7, DUMMY8,
     4                 R3(N,N), XJ3, RR(N,N),
     5                 DAMP2, DAMP3, RJ4J1, DSAVE(N,N)
      CHARACTER(LEN=5) ZSYM(N)
C
C STATEMENT FUNCTIONS
C
C     F2(DUMMY)=1.0D0-DEXP(-DUMMY/DAMP)*(1.0D0+DUMMY/DAMP+DUMMY**2/(2.0D0*DAMP**2))
C     F2P(DUMMY)=DEXP(-DUMMY/DAMP)*DUMMY**2/(2.0D0*DAMP**3)

      DAMP2=DAMP**2
      DAMP3=DAMP**3
C 
C  STORE DISTANCE MATRICES.
C
      DO J1=1,N
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         R3(J1,J1)=0.0D0
         DOT(J1,J1)=X(3*J1-2)**2 + X(3*J1-1)**2 + X(3*J1)**2
         NPROD(J1,J1)=1
         DO J2=J1+1,N
            DUMMY=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1           +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2           +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            RR(J2,J1)=DSQRT(DUMMY)
            DSAVE(J2,J1)=DEXP(-RR(J2,J1)/DAMP)
            DUMMY=1.0D0/DUMMY
            R2(J2,J1)=DUMMY
            R(J2,J1)=DSQRT(DUMMY)
            R3(J2,J1)=R(J2,J1)*DUMMY
            R(J1,J2)=R(J2,J1)
            RR(J1,J2)=RR(J2,J1)
            R2(J1,J2)=R2(J2,J1)
            R3(J1,J2)=R3(J2,J1)
            DSAVE(J1,J2)=DSAVE(J2,J1)
            DOT(J2,J1)=X(3*J2-2)*X(3*J1-2) + X(3*J2-1)*X(3*J1-1) + X(3*J2)*X(3*J1)
            DOT(J1,J2)=DOT(J2,J1)
            IF (ZSYM(J2).EQ.ZSYM(J1)) THEN
               NPROD(J2,J1)=1
               NPROD(J1,J2)=1
            ELSE
               NPROD(J2,J1)=-1
               NPROD(J1,J2)=-1
            ENDIF
         ENDDO
      ENDDO

      DO J1=1,3*N
         V(J1)=0.0D0
      ENDDO
C
C  GRADIENT.
C
      DO J1=1,N           !   J1 = L
         DO J2=1,3        !   J2 = ALPHA
            J3=3*J1-3+J2 
            XJ3=X(J3)
            DO J4=1,N        !  J4 = I
               DUMMY=0.0D0
               DUMMY3=DOT(J4,J4)-DOT(J4,J1)
               DUMMY4=(XJ3-X(3*J4-3+J2))*R2(J4,J1)
               DUMMY5=1.0D0-DSAVE(J4,J1)*(1.0D0+RR(J4,J1)/DAMP+RR(J4,J1)**2/(2.0D0*DAMP2))
               DUMMY6=DSAVE(J4,J1)*RR(J4,J1)**2/(2.0D0*DAMP**3)
               DUMMY7=X(3*J4-3+J2)
               DUMMY8=R(J4,J1)
               DO J5=1,N     !  J5 = K
                  DUMMY=DUMMY-(NPROD(J5,J1)*R3(J5,J4)*(1.0D0 + 
     1      DSAVE(J5,J4)*(-1.0D0 - RR(J5,J4)/DAMP - (0.5D0*RR(J5,J4)**2)/DAMP2))*
     2    (DUMMY6*DUMMY8*(DUMMY7 - XJ3)*(DUMMY3 + DOT(J5,J1) - DOT(J5,J4)) + 
     3      DUMMY5*(DUMMY7 + 3.0D0*DUMMY4*(DUMMY3 + DOT(J5,J1) - DOT(J5,J4)) - X(-3 + J2 + 3*J5))) )
               ENDDO
               IF (ZSYM(J4).EQ.'PL') THEN
                  V(J3)=V(J3)-ALPHAP*R3(J4,J1)*DUMMY
               ELSE
                  V(J3)=V(J3)-ALPHAM*R3(J4,J1)*DUMMY
               ENDIF
            ENDDO

            DUMMY=0.0D0
            DO J4=1,N            ! J4 = J
               DUMMY3=DOT(J1,J1)-DOT(J4,J1)
               DUMMY4=(XJ3-X(3*J4-3+J2))*R2(J4,J1)
               DUMMY5=R3(J4,J1)
               DUMMY6=X(3*J4-3+J2)
               DUMMY7=1.0D0-DSAVE(J4,J1)*(1.0D0+RR(J4,J1)/DAMP+RR(J4,J1)**2/(2.0D0*DAMP2))
               DUMMY8=DSAVE(J4,J1)*RR(J4,J1)**2/(2.0D0*DAMP**3)
               RJ4J1=R(J4,J1)
               DO J5=1,N         ! J5 = K
                  DUMMY=DUMMY+DUMMY5*NPROD(J5,J4)*R3(J5,J1)*
     1  ( (1.0D0-DSAVE(J5,J1)*(1.0D0+RR(J5,J1)/DAMP+RR(J5,J1)**2/(2.0D0*DAMP2)))
     2       *(DUMMY8*(-DUMMY6 + XJ3)*(DUMMY3 - DOT(J5,J1) + DOT(J5,J4))*
     2        RJ4J1 + DUMMY7*(-DUMMY6 + 2.0D0*XJ3 - 
     3          3.0D0*(DUMMY3 - DOT(J5,J1) + DOT(J5,J4))*
     4           (DUMMY4 + R2(J5,J1)*(XJ3 - X(-3 + J2 + 3*J5))) - X(-3 + J2 + 3*J5))
     5       ) + DUMMY7*(DUMMY3 - DOT(J5,J1) + DOT(J5,J4))*
     6           DSAVE(J5,J1)*RR(J5,J1)/(2.0D0*DAMP3)*(XJ3 - X(-3 + J2 + 3*J5)))
               ENDDO
            ENDDO
            IF (ZSYM(J1).EQ.'PL') THEN
               V(J3)=V(J3)-ALPHAP*DUMMY/2.0D0
            ELSE
               V(J3)=V(J3)-ALPHAM*DUMMY/2.0D0
            ENDIF
         ENDDO
      ENDDO

      RETURN 
      END
C
C  ANALYTIC SECOND DERIVATIVES OF POSSIBLE FIRST ORDER INDUCTION.
C
      SUBROUTINE TFSEC(N, X, ALPHAP, ALPHAM, ZSYM, DAMP, NPROD)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NPROD(N,N), J7
      DOUBLE PRECISION X(3*N), DAMP,
     1                 R(N,N), 
     2                 DUMMY, ALPHAP, ALPHAM, 
     3                 DOT(N,N), R2(N,N),
     4                 R3(N,N), XJ3, XJ6, DOT11, DOT44, DOT41,
     5                 XJ1J5, XJ4J2, R2J4J1, R3J4J1, EX, AL(N), XJ5, R3J6J1,
     6                 XJ6J2, XJ6J4, DOT66, DOT61, R2J6J1, RR(N,N)
      CHARACTER(LEN=5) ZSYM(N)
C
C  STATEMENT FUNCTIONS.
C
C     F2(DUMMY)=1.0D0-DEXP(-DUMMY/DAMP)*(1.0D0+DUMMY/DAMP+DUMMY**2/(2.0D0*DAMP**2))
C     F2P(DUMMY)=DEXP(-DUMMY/DAMP)*DUMMY**2/(2.0D0*DAMP**3)
C     F2PP(DUMMY)=DEXP(-DUMMY/DAMP)*DUMMY*(1.0D0-DUMMY/(2.0D0*DAMP))/DAMP**3
C 
C  STORE DISTANCE MATRICES.
C
      DO J1=1,N
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         DOT(J1,J1)=X(3*J1-2)**2 + X(3*J1-1)**2 + X(3*J1)**2
         IF (ZSYM(J1).EQ.'PL') THEN
            AL(J1)=ALPHAP
         ELSE
            AL(J1)=ALPHAM
         ENDIF
         DO J2=J1+1,N
            DUMMY=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1           +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2           +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            RR(J2,J1)=DSQRT(DUMMY)
            DUMMY=1.0D0/DUMMY
            R2(J2,J1)=DUMMY
            R(J2,J1)=DSQRT(DUMMY)
            R3(J2,J1)=R(J2,J1)*DUMMY
            RR(J1,J2)=RR(J2,J1)
            R(J1,J2)=R(J2,J1)
            R2(J1,J2)=R2(J2,J1)
            R3(J1,J2)=R3(J2,J1)
            DOT(J2,J1)=X(3*J2-2)*X(3*J1-2) + X(3*J2-1)*X(3*J1-1) + X(3*J2)*X(3*J1)
            DOT(J1,J2)=DOT(J2,J1)
         ENDDO
      ENDDO

      DO J1=1,3*N
         DO J2=1,3*N
            HESS(J2,J1)=0.0D0
         ENDDO
      ENDDO
C
C  SECOND DERIVATIVES - COORDINATES ON THE SAME ATOM.
C
      DO J1=1,N                         ! J1 = L     
         DOT11=DOT(J1,J1)
         DO J2=1,3                      ! J2 = ALPHA
            J3=3*J1-3+J2
            XJ3=X(J3)
            DO J4=J2,3                  ! J4 = BETA
               J5=3*J1-3+J4
               XJ5=X(J5)
               EX=0.0D0
               IF (J4.EQ.J2) EX=1.0D0
               DO J6=1,N                ! J6 = I
                  DOT66=DOT(J6,J6)
                  DOT61=DOT(J6,J1)
                  R2J6J1=R2(J6,J1)
                  R3J6J1=R3(J6,J1)
                  XJ6J2=X(3*J6-3+J2)
                  XJ6J4=X(3*J6-3+J4)
                  DUMMY=2.0D0*EX*R3J6J1**2 - 12.0D0*R2J6J1**4*(XJ3 - XJ6J2)*(XJ5 - XJ6J4) + 
     1         18.0D0*(DOT11 - 2.0D0*DOT61 + DOT66)*R2J6J1**5*(XJ3 - XJ6J2)*(XJ5 - XJ6J4)
                  DO J7=1,N             ! J7 = K
                     DUMMY=DUMMY+2.0D0*R2J6J1*R3J6J1*NPROD(J7,J1)*R3(J7,J6)*
     1                    (3.0D0*EX*(DOT61 - DOT66 - DOT(J7,J1) + DOT(J7,J6)) - 
     2                  15.0D0*R2J6J1*(XJ3 - XJ6J2)*(XJ5 - XJ6J4)*(DOT61 - DOT66 - DOT(J7,J1) + DOT(J7,J6)) + 
     3                    3.0D0*(XJ5 - XJ6J4)*(XJ6J2 - X(-3 + J2 + 3*J7)) + 
     4                    3.0D0*(XJ3 - XJ6J2)*(XJ6J4 - X(-3 + J4 + 3*J7)))
                  ENDDO
                  HESS(J5,J3)=HESS(J5,J3)-AL(J6)*DUMMY/2.0D0
               ENDDO
C
C  CASE L = I.
C
               DUMMY=0.0D0
               DO J6=1,N              ! J6 = J
                  R3J6J1=R3(J6,J1)
                  R2J6J1=R2(J6,J1)
                  XJ6J2=X(3*J6-3+J2)
                  XJ6J4=X(3*J6-3+J4)
                  DOT61=DOT(J6,J1)
                  DO J7=1,N           ! J7 = K
                     DUMMY=DUMMY+R3J6J1*NPROD(J7,J6)*R3(J7,J1)*(EX*
     1     (2.0D0 - 3.0D0*(DOT11 - DOT61 - DOT(J7,J1) + DOT(J7,J6))*(R2J6J1 + R2(J7,J1))) - 
     2    3.0D0*(2.0D0*XJ3 - XJ6J2 - 3.0D0*(DOT11 - DOT61 - DOT(J7,J1) + DOT(J7,J6))*
     3        (R2J6J1*(XJ3 - XJ6J2) + R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7))) - X(-3 + J2 + 3*J7))*
     4     (R2J6J1*(XJ5 - XJ6J4) + R2(J7,J1)*(XJ5 - X(-3 + J4 + 3*J7))) + 
     5    6.0D0*(DOT11 - DOT61 - DOT(J7,J1) + DOT(J7,J6))*
     6     (R2J6J1**2*(XJ3 - XJ6J2)*(XJ5 - XJ6J4) + 
     7       R2(J7,J1)**2*(XJ3 - X(-3 + J2 + 3*J7))*(XJ5 - X(-3 + J4 + 3*J7))) - 
     8    3.0D0*(R2J6J1*(XJ3 - XJ6J2) + R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7)))*
     9     (2.0D0*XJ5 - XJ6J4 - X(-3 + J4 + 3*J7)))
                  ENDDO 
               ENDDO
               HESS(J5,J3)=HESS(J5,J3)-AL(J1)*DUMMY/2.0D0
               HESS(J3,J5)=HESS(J5,J3)
            ENDDO
         ENDDO
      ENDDO
C
C  SECOND DERIVATIVES - DIFFERENT ATOMS.
C
      DO J1=1,N             ! J1 = L     
         DOT11=DOT(J1,J1)
         DO J2=1,3          ! J2 = ALPHA
            J3=3*J1-3+J2 
            XJ3=X(J3)
            DO J4=J1+1,N    ! J4 = M
               DOT44=DOT(J4,J4)
               DOT41=DOT(J4,J1)
               XJ4J2=X(3*J4-3+J2)
               R2J4J1=R2(J4,J1)
               R3J4J1=R3(J4,J1)
               DO J5=1,3    ! J5 = BETA
                  J6=3*J4-3+J5
                  XJ1J5=X(3*J1-3+J5)
                  XJ6=X(J6)
                  EX=0.0D0
                  IF (J5.EQ.J2) EX=1.0D0
                  DO J7=1,N   ! J7 = I
                     DUMMY=2.0D0*NPROD(J4,J1)*R3(J7,J1)*R3(J7,J4)*
     1  (EX - 3.0D0*(R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7))*(XJ1J5 - X(-3 + J5 + 3*J7)) + 
     2       R2(J7,J4)*(XJ4J2 - X(-3 + J2 + 3*J7))*(XJ6 - X(-3 + J5 + 3*J7))) + 
     3    9.0D0*(DOT41 - DOT(J7,J1) - DOT(J7,J4) + DOT(J7,J7))*R2(J7,J1)*R2(J7,J4)*
     -     (XJ3 - X(-3 + J2 + 3*J7))*(XJ6 - X(-3 + J5 + 3*J7)))
                     HESS(J6,J3)=HESS(J6,J3)-AL(J7)*DUMMY/2.0D0
                  ENDDO
C
C  I = L TERM
C
                  DUMMY=0.0D0
                  DO J7=1,N    ! J7 = J
                     DUMMY=DUMMY+2.0D0*R3J4J1*NPROD(J7,J4)*R3(J7,J1)*
     1  (-EX + 3.0D0*EX*R2J4J1*(DOT11 - DOT41 - DOT(J7,J1) + DOT(J7,J4)) - 
     2    3.0D0*R2J4J1*(-XJ1J5 + XJ6)*(2.0D0*XJ3 - XJ4J2 - 
     3       3.0D0*(DOT11 - DOT41 - DOT(J7,J1) + DOT(J7,J4))*
     4        (R2J4J1*(XJ3 - XJ4J2) + R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7))) - X(-3 + J2 + 3*J7))
     5     - 3.0D0*(2.0D0*R2J4J1**2*(XJ3 - XJ4J2)*(XJ1J5 - XJ6)*
     6        (DOT11 - DOT41 - DOT(J7,J1) + DOT(J7,J4)) + 
     7       (R2J4J1*(XJ3 - XJ4J2) + R2(J7,J1)*(XJ3 - X(-3 + J2 + 3*J7)))*
     8        (-XJ1J5 + X(-3 + J5 + 3*J7))))
                  ENDDO
                  HESS(J6,J3)=HESS(J6,J3)-AL(J1)*DUMMY/2.0D0
C
C  I = M TERM
C
                  DUMMY=0.0D0
                  DO J7=1,N    ! J7 = J
                     DUMMY=DUMMY+2.0D0*R3J4J1*NPROD(J7,J1)*R3(J7,J4)*
     1  (-EX - 3.0D0*EX*R2J4J1*(DOT41 - DOT44 - DOT(J7,J1) + DOT(J7,J4)) - 
     2    3.0D0*(-2.0D0*R2J4J1**2*(XJ3 - XJ4J2)*(XJ1J5 - XJ6)*
     3        (DOT41 - DOT44 - DOT(J7,J1) + DOT(J7,J4)) + 
     4       (-XJ4J2 + X(-3 + J2 + 3*J7))*
     5        (R2J4J1*(-XJ1J5 + XJ6) + R2(J7,J4)*(XJ6 - X(-3 + J5 + 3*J7)))) - 
     6    3.0D0*R2J4J1*(XJ3 - XJ4J2)*(-XJ1J5 + 2.0D0*XJ6 - 
     7       3.0D0*(-DOT41 + DOT44 + DOT(J7,J1) - DOT(J7,J4))*
     8        (R2J4J1*(-XJ1J5 + XJ6) + R2(J7,J4)*(XJ6 - X(-3 + J5 + 3*J7))) - X(-3 + J5 + 3*J7)))
                  ENDDO
                  HESS(J6,J3)=HESS(J6,J3)-AL(J4)*DUMMY/2.0D0
                  HESS(J3,J6)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN 
      END

