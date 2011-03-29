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

C  ENERGY AND DERIVATIVES OF A GENERAL BORN-MAYER POTENTIAL WITH
C  TOSI-FUMI PARAMETERS SPECIFIED AFTER THE KEYWORD TOSI IN THE
C  ODATA FILE. ENERGY IS IN HARTREE, LENGTH IN BOHR.
C
      SUBROUTINE TOSIFUMI(X, V, POTEL, GTEST, STEST)
      USE COMMONS
      IMPLICIT NONE 
      INTEGER N, I, J, J1, J2, K, L, J3, J4
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), POTEL, 
     1                 DUM, AC(NATOMS,NATOMS), Q(NATOMS), RD(NATOMS,NATOMS), RRD(NATOMS,NATOMS)
      LOGICAL GTEST, STEST

      N=NATOMS
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
         RD(I,I)=0.0D0
         RRD(I,I)=0.0D0
         DO J=I+1,N
            L=K+3*(J-I)
            RD(I,J)= DSQRT((X(K)-X(L))*(X(K)-X(L))+
     1                     (X(K+1)-X(L+1))*(X(K+1)-X(L+1))+
     2                     (X(K+2)-X(L+2))*(X(K+2)-X(L+2))) 
            RRD(I,J)=1.0D0/RD(I,J)
            RD(J,I)=RD(I,J)
            RRD(J,I)=RRD(I,J)
         ENDDO
      ENDDO

      POTEL=0.0D0
      DO I=1,N
         DO J=I+1,N
            POTEL=POTEL+Q(I)*Q(J)*RRD(J,I)+AC(J,I)*DEXP(-RD(J,I)/RHO)
         ENDDO
      ENDDO

      IF (.NOT.GTEST) RETURN
C
C  GRADIENT.
C
      DO J1=1,N
         DO J2=1,3
            DUM=0.0D0
            J4=3*(J1-1)+J2
            DO J3=1,N
               DUM=DUM-(Q(J3)*Q(J1)*RRD(J3,J1)**2+AC(J3,J1)*DEXP(-RD(J3,J1)/RHO)/RHO)*RRD(J3,J1)*(X(J4)-X(3*(J3-1)+J2))
            ENDDO
            V(J4)=DUM
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
               CALL NSHESS(I,J,K,L,X,N,AC,Q,RD)
            ELSE
               CALL NDHESS(I,J,K,L,X,N,AC,Q,RD)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
C....................................................
C.....THIS SUBROUTINE IF THE COORDINATES ARE THE SAME
 
      SUBROUTINE NSHESS(I,J,K,L,CFG,NUM,AC,Q,RD)
      USE COMMONS
      IMPLICIT NONE
      INTEGER I,J,K,L,I1,J1,M,N,NUM
      DOUBLE PRECISION CFG(3*NATOMS), RW, R, QQ, AB,
     1                 EP, DR, DF,AC(NATOMS,NATOMS),Q(NATOMS),RD(NATOMS,NATOMS)

      RW=RHO
 
C.....I1 AND J1 ARE THE PARTICLE NUMBERS FOR X(I) AND X(J)

      I1 = ((I-1)/3) + 1
      J1 = ((J-1)/3) + 1
 
C.....IF THE PARTICLES ARE THE SAME DO THE NEXT TWO LOOPS
C.....THERE ARE TWO LOOPS SO THAT WE CAN SKIP PARTICLE I1
C.....WITHOUT RESORTING TO AN IF STATEMENT

      IF(I1 .EQ. J1)THEN
      
      DO M = 1, I1-1
        
         N = 3*(M-1) + K + 1
 
         R = RD(I1,M)
         QQ = Q(I1)*Q(M)
         AB = AC(I1,M)
         EP = DEXP(-R/RW)
         DR = (CFG(I)-CFG(N))*(CFG(I)-CFG(N))
      
         DF = -(QQ/(R*R*R)) + (3.0D0*QQ*DR/(R*R*R*R*R))
     1        + (AB*DR*EP/(RW*R*R*R)) + (AB*DR*EP/(RW*RW*R*R))
     2        - (AB*EP/(RW*R))
 
         HESS(I,I) = HESS(I,I) + DF
      
      ENDDO
 
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
 
      SUBROUTINE NDHESS(I,J,K,L,CFG,NUM,AC,Q,RD)
      USE COMMONS
      IMPLICIT NONE
      INTEGER I,J,K,L,I1,J1,M,N,NUM,N1
      DOUBLE PRECISION CFG(3*NATOMS),RW,R,QQ,AB,
     1                 EP,DR,DF,AC(NATOMS,NATOMS),Q(NATOMS),RD(NATOMS,NATOMS)
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
