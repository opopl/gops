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
C  TWO BODY TERM OF MURRELL POTENTIAL
C  SUBROUTINE V2DIFF CALCULATES THE CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY. 
C
C*************************************************************************
C
      SUBROUTINE JM2(N, X, V)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION RE, D, AN
      PARAMETER(D=2.918D0)
      PARAMETER(AN=6.50D0)
      PARAMETER(RE=2.389D0)
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), ABY, ABX, 
     2                 RAB(N,N), RRAB(N,N),
     3                 RHOAB(N,N)
C
C  STORE DISTANCE MATRICES.
C
      DO 20 J1=1,N
         RRAB(J1,J1)=0.0D0
         RAB(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            RAB(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            RAB(J2,J1)=DSQRT(RAB(J2,J1))
            RAB(J1,J2)=RAB(J2,J1)
            RRAB(J2,J1)=1.0D0/RAB(J2,J1)
            RRAB(J1,J2)=RRAB(J2,J1)
            RHOAB(J2,J1)=(RAB(J2,J1)-RE)/RE
            RHOAB(J1,J2)=RHOAB(J2,J1)
10       CONTINUE
20    CONTINUE
C
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
               ABX=X(3*(J4-1)+J2)-X(J3)
               V(J3)=V(J3)
     1          +AN*ABX*D*(RRAB(J4,J1)-(1+AN*RHOAB(J4,J1))*RRAB(J4,J1)) 
     2           /(DEXP(AN*RHOAB(J4,J1))*RE) 
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
               ABX=X(3*(J4-1)+J2)-X(J3)
               HESS(J3,J3)=HESS(J3,J3)+
     1   D*(AN**2*ABX**2*(1-AN*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+ 
     2   AN*(AN*RHOAB(J4,J1)*RRAB(J4,J1)+ABX**2*(RRAB(J4,J1)**3-
     3   (1+AN*RHOAB(J4,J1))*RRAB(J4,J1)**3))/RE)/DEXP(AN*RHOAB(J4,J1)) 
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
                  ABX=-X(J3)+X(3*(J5-1)+J2)
                  ABY=-X(3*(J1-1)+J4)+X(3*(J5-1)+J4)
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+
     1    ABX*ABY*D*(AN**2*(1-AN*RHOAB(J5,J1))*RRAB(J5,J1)**2/RE**2+  
     2    AN*(RRAB(J5,J1)**3-(1+AN*RHOAB(J5,J1))*RRAB(J5,J1)**3)/RE)/ 
     3    DEXP((AN*RHOAB(J5,J1)))
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
              ABX=-X(J3)+X(3*(J4-1)+J2)
                HESS(3*(J4-1)+J2,J3)=
     1     D*(AN**2*ABX**2*(-1+AN*RHOAB(J4,J1))/(RAB(J4,J1)**2*RE**2)+   
     2     AN*(RRAB(J4,J1)-(1+AN*RHOAB(J4,J1))*RRAB(J4,J1)+
     3     AN*ABX**2*RHOAB(J4,J1)*RRAB(J4,J1)**3)/RE)/
     4     DEXP((AN*RHOAB(J4,J1)))
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
                  ABX=-X(J3)+X(3*(J4-1)+J2)
                  ABY=-X(3*(J1-1)+J5)+X(J6)
                  HESS(J6,J3)=
     1          AN**2*ABX*ABY*D*((-1+AN*RHOAB(J4,J1))/
     2          (RAB(J4,J1)**2*RE**2)+
     3          RHOAB(J4,J1)*RRAB(J4,J1)**3/RE)/DEXP(AN*RHOAB(J4,J1)) 
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
C
C*************************************************************************
C
C  TWO BODY TERM OF MURRELL POTENTIAL - ANALYTIC DERIVATIVES WITH
C  PERIODIC BOUNDARY CONDITIONS
C
C*************************************************************************
C
      SUBROUTINE JM2C (N, X, V, CUTOFF, RAB, RRAB, VEC, RHOAB, M, NDUM, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      INTEGER M(N,N,3)
      LOGICAL STEST
      DOUBLE PRECISION RAB(N,N), RRAB(N,N), 
     1              VEC(N,N,3), RHOAB(N,N), NDUM(N,N)
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF, SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), ABY, ABX
C 
C
C  STORE DISTANCE MATRICES.
C
      DO 20 J1=1,N
         DO 10 J2=J1+1,N
            IF (RAB(J2,J1).GE.CUTOFF) THEN
               RRAB(J2,J1)=0.0D0
               RRAB(J1,J2)=0.0D0
            ENDIF
            RHOAB(J2,J1)=(RAB(J2,J1)-RE)/RE
            RHOAB(J1,J2)=RHOAB(J2,J1)
10       CONTINUE
20    CONTINUE
C
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
                  ABX=VEC(J4,J1,J2)
                  V(J3)=V(J3)
     1              +AN2*ABX*D*(RRAB(J4,J1)-
     2              (1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)) 
     3              /(DEXP(AN2*RHOAB(J4,J1))*RE) 
30          CONTINUE
C           PRINT*,'J3,V(J3)=',J3,V(J3)
40       CONTINUE
50    CONTINUE
      IF (.NOT.STEST) RETURN
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO 60 J4=1,N
               ABX=VEC(J4,J1,J2)
               HESS(J3,J3)=HESS(J3,J3)+
     1 D*(AN2**2*ABX**2*(1-AN2*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+ 
     2 AN2*(AN2*RHOAB(J4,J1)*RRAB(J4,J1)+ABX**2*(RRAB(J4,J1)**3-
     3 (1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)**3))/RE)
     4 /DEXP(AN2*RHOAB(J4,J1)) 
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
                  ABX=VEC(J5,J1,J2)
                  ABY=VEC(J5,J1,J4)
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+
     1    ABX*ABY*D*(AN2**2*(1-AN2*RHOAB(J5,J1))*RRAB(J5,J1)**2/RE**2+  
     2    AN2*(RRAB(J5,J1)**3-(1+AN2*RHOAB(J5,J1))*RRAB(J5,J1)**3)/RE)/ 
     3    DEXP((AN2*RHOAB(J5,J1)))
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
              ABX=VEC(J4,J1,J2)
                HESS(3*(J4-1)+J2,J3)=
     1     D*(AN2**2*ABX**2*(-1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+ 
     2     AN2*(RRAB(J4,J1)-(1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)+
     3     AN2*ABX**2*RHOAB(J4,J1)*RRAB(J4,J1)**3)/RE)/
     4     DEXP((AN2*RHOAB(J4,J1)))
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
                  ABX=VEC(J4,J1,J2)
                  ABY=VEC(J4,J1,J5)
                  HESS(J6,J3)=AN2**2*ABX*ABY*D*((-1+AN2*RHOAB(J4,J1))*
     1                     RRAB(J4,J1)**2/RE**2+RHOAB(J4,J1)
     2                    *RRAB(J4,J1)**3/RE)/DEXP(AN2*RHOAB(J4,J1)) 

                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=
     1          + AN2**2*ABX*ABY*D*((-1+AN2*RHOAB(J4,J1))*
     2          RRAB(J4,J1)**2/RE**2+RHOAB(J4,J1)*RRAB(J4,J1)**3/RE)/
     3          DEXP(AN2*RHOAB(J4,J1)) 
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
C
C*************************************************************************
C
C  TWO BODY TERM OF MURRELL POTENTIAL - ANALYTIC DERIVATIVES
C
C*************************************************************************
C
      SUBROUTINE JM2CC(N, X, V,RAB,RRAB,VEC,RHOAB,M,NDUM)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      INTEGER M(N,N,3)
      DOUBLE PRECISION RAB(N,N), RRAB(N,N),
     1              VEC(N,N,3), RHOAB(N,N), NDUM(N,N)
      DOUBLE PRECISION CUTOFF
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), ABY, ABX, ABIG, B, C
      PARAMETER (ABIG=13.29D0, B=0.574D0, C=1.0D0, CUTOFF=2.8D0)
C   
C      INITIALISE THE FIRST AND SECOND DERIVATIVES
C
      DO 21 J1=1,3*N
         V(J1)=0.0D0
         DO 22 J2=1,3*N
            HESS(J2,J1)=0.0D0
22       CONTINUE
21    CONTINUE
C
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            DO 30 J4=1,N
            IF ((RAB(J4,J1).LT.CUTOFF).AND.(J4.NE.J1)) THEN
               ABX=VEC(J4,J1,J2)
               V(J3)=V(J3)+
     1         ABIG*ABX*DEXP(C/RHOAB(J4,J1))*(4.0D0*B/RAB(J4,J1)**6 
     2      + C*(-1.0D0 + B/RAB(J4,J1)**4)*RRAB(J4,J1)/RHOAB(J4,J1)**2)
            ENDIF
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
            DO 60 J4=1,N
            IF ((RAB(J4,J1).GE.CUTOFF).OR.(J4.EQ.J1)) GOTO 60
               ABX=VEC(J4,J1,J2)
               HESS(J3,J3)=HESS(J3,J3)+
     1   ABIG*DEXP(C/RHOAB(J4,J1))*(B*(-4/RAB(J4,J1)**6 + 
     2   ABX**2*(24/RAB(J4,J1)**8+8*C/(RAB(J4,J1)**7*RHOAB(J4,J1)**2))) 
     3   + (-1 + B/RAB(J4,J1)**4)*(-(C*RRAB(J4,J1)/RHOAB(J4,J1)**2) + 
     4   ABX**2*((C**2/RHOAB(J4,J1)**4 +
     5   2*C/RHOAB(J4,J1)**3)/RAB(J4,J1)**2 + 
     6   C*RRAB(J4,J1)**3/RHOAB(J4,J1)**2)))
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
                  IF ((RAB(J5,J1).GE.CUTOFF).OR.(J5.EQ.J1)) GOTO 90 
                  ABX=VEC(J5,J1,J2)
                  ABY=VEC(J5,J1,J4)
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+
     1    ABIG*ABX*ABY*DEXP(C/RHOAB(J5,J1))*(B*(24/RAB(J5,J1)**8 +
     2    8*C/(RAB(J5,J1)**7*RHOAB(J5,J1)**2)) +
     3    (-1.0D0+B/RAB(J5,J1)**4)*((C**2/RHOAB(J5,J1)**4 + 
     4    2*C/RHOAB(J5,J1)**3)/RAB(J5,J1)**2 + 
     5    C*RRAB(J5,J1)**3/RHOAB(J5,J1)**2))
90             CONTINUE
            HESS(J3,3*(J1-1)+J4)=HESS(3*(J1-1)+J4,J3)
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
                IF (RAB(J4,J1).GE.CUTOFF) GOTO 130 
                ABX=VEC(J4,J1,J2)
                HESS(3*(J4-1)+J2,J3)=
     1    ABIG*DEXP(C/RHOAB(J4,J1))*(B*(4/RAB(J4,J1)**6 +
     2    ABX**2*(-24/RAB(J4,J1)**8 
     3    - 8*C/(RAB(J4,J1)**7*RHOAB(J4,J1)**2))) + 
     4    (-1.0D0 + B/RAB(J4,J1)**4)*(C*RRAB(J4,J1)/RHOAB(J4,J1)**2 + 
     5    ABX**2*((-C**2/RHOAB(J4,J1)**4 -
     6    2*C/RHOAB(J4,J1)**3)/RAB(J4,J1)**2 - 
     7    C*RRAB(J4,J1)**3/RHOAB(J4,J1)**2)))
            HESS(J3,3*(J4-1)+J2)=HESS(3*(J4-1)+J2,J3)
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
               IF (RAB(J4,J1).GE.CUTOFF) GOTO 160 
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                  ABX=VEC(J4,J1,J2)
                  ABY=VEC(J4,J1,J5)
                  HESS(J6,J3)=
     1    ABIG*ABX*ABY*DEXP(C/RHOAB(J4,J1))*(B*(-24/RAB(J4,J1)**8 -  
     2    8*C/(RAB(J4,J1)**7*RHOAB(J4,J1)**2)) + 
     3    (-1.0D0+B/RAB(J4,J1)**4)*((-C**2/RHOAB(J4,J1)**4 -
     4    2.0D0*C/RHOAB(J4,J1)**3)/RAB(J4,J1)**2 - 
     5    C*RRAB(J4,J1)**3/RHOAB(J4,J1)**2))
               HESS(J3,J6)=HESS(J6,J3)
               HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               HESS(3*(J1-1)+J5,3*(J4-1)+J2)=HESS(J6,J3)
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
C
C*************************************************************************
C
C  TWO BODY TERM OF MURRELL POTENTIAL - ANALYTIC DERIVATIVES WITH
C  PERIODIC BOUNDARY CONDITIONS
C
C*************************************************************************
C
      SUBROUTINE JM2P (N, X, V, BOXLX, BOXLY, BOXLZ, CUTOFF, RAB, VEC)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION RAB(N,N), RRAB(N,N),
     1              VEC(N,N,3), RHOAB(N,N)
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF
      DOUBLE PRECISION X(3*N), V(3*N), ABY, ABX,
     1                 SR2, SR3, SR6, C0, C1, C2, C3, C4, BOXLX, BOXLY, BOXLZ,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
C 
C
C  STORE DISTANCE MATRICES.
C
      DO 20 J1=1,N
         DO 10 J2=J1+1,N
            IF (RAB(J2,J1).GE.CUTOFF) THEN
               RRAB(J2,J1)=0.0D0
               RRAB(J1,J2)=0.0D0
            ENDIF
            RHOAB(J2,J1)=(RAB(J2,J1)-RE)/RE
            RHOAB(J1,J2)=RHOAB(J2,J1)
10       CONTINUE
20    CONTINUE
C
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
                  ABX=VEC(J4,J1,J2)
                  V(J3)=V(J3)
     1              +AN2*ABX*D*(RRAB(J4,J1)-
     2              (1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)) 
     3              /(DEXP(AN2*RHOAB(J4,J1))*RE) 
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
C           HESS(J3,J3)=0.0D0
            DO 60 J4=1,N
               ABX=VEC(J4,J1,J2)
               HESS(J3,J3)=HESS(J3,J3)+
     1 D*(AN2**2*ABX**2*(1-AN2*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+ 
     2 AN2*(AN2*RHOAB(J4,J1)*RRAB(J4,J1)+ABX**2*(RRAB(J4,J1)**3-
     3 (1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)**3))/RE)
     4 /DEXP(AN2*RHOAB(J4,J1)) 
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
C               HESS(3*(J1-1)+J4,J3)=0.0D0
                DO 90 J5=1,N
                  ABX=VEC(J5,J1,J2)
                  ABY=VEC(J5,J1,J4)
                  HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+
     1    ABX*ABY*D*(AN2**2*(1-AN2*RHOAB(J5,J1))*RRAB(J5,J1)**2/RE**2+  
     2    AN2*(RRAB(J5,J1)**3-(1+AN2*RHOAB(J5,J1))*RRAB(J5,J1)**3)/RE)/ 
     3    DEXP((AN2*RHOAB(J5,J1)))
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
              ABX=VEC(J4,J1,J2)
                HESS(3*(J4-1)+J2,J3)=HESS(3*(J4-1)+J2,J3) +
     1     D*(AN2**2*ABX**2*(-1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)**2/RE**2+   
     2     AN2*(RRAB(J4,J1)-(1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)+
     3     AN2*ABX**2*RHOAB(J4,J1)*RRAB(J4,J1)**3)/RE)/
     4     DEXP((AN2*RHOAB(J4,J1)))
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
                  ABX=VEC(J4,J1,J2)
                  ABY=VEC(J4,J1,J5)
                  HESS(J6,J3)=
     1          AN2**2*ABX*ABY*D*((-1+AN2*RHOAB(J4,J1))*
     2          RRAB(J4,J1)**2/RE**2+
     3          RHOAB(J4,J1)*RRAB(J4,J1)**3/RE)/DEXP(AN2*RHOAB(J4,J1)) 
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=
     1          + AN2**2*ABX*ABY*D*((-1+AN2*RHOAB(J4,J1))*
     2          RRAB(J4,J1)**2/RE**2+
     3          RHOAB(J4,J1)*RRAB(J4,J1)**3/RE)/DEXP(AN2*RHOAB(J4,J1)) 
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
C
C*************************************************************************
C
C  DERIVATIVES OF THE JM THREE-BODY TERMS - TO SAVE SPACE THESE ARE
C  ADDED DIRECTLY TO THE SUPPLIED GRADIENT AND SECOND DERIVATIVE
C  MATRIX.
C                                        
C*************************************************************************
C
      SUBROUTINE JM3(N, X, V)
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5
      DOUBLE PRECISION SR2, SR3, SR6, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, RE, D, AN
      PARAMETER(SR2=1.414213562D0)
      PARAMETER(SR3=1.732050808D0)
      PARAMETER(SR6=2.449489743D0)
      PARAMETER(C0=3.598D0)
      PARAMETER(C1=-11.609D0)
      PARAMETER(C2=13.486D0)
      PARAMETER(C3=-18.174D0)
      PARAMETER(C4=-5.570D0)
      PARAMETER(C5=79.210D0)
      PARAMETER(C6=-6.458D0)
      PARAMETER(C7=23.383D0)
      PARAMETER(C8=-111.890D0)
      PARAMETER(C9=9.705D0)
      PARAMETER(C10=38.297D0)
      PARAMETER(RE=2.389D0)
      PARAMETER(D=2.918D0)
      PARAMETER(AN=6.50D0)
      DOUBLE PRECISION X(3*N), ABX, ACX, BCX,
     1                 R2(N,N), V(3*N), 
     2                 VEC(N,N,3), TEMP, RR2(N,N), 
     3                 RAB, RRAB, RAC, RRAC, RBC, RRBC, 
     4                 ABY, ACY, BCY, RHOAB, RHOBC ,RHOAC ,QA, QQ2,
     5                 QQ3, QB, QC
      DO 20 J1=1,N
         R2(J1,J1)=0.0D0
         RR2(J1,J1)=0.0D0
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 10 J2=J1+1,N
            R2(J2,J1)=
     1                ( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     2                ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     3                ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2 
            R2(J2,J1)=DSQRT(R2(J2,J1))
            RR2(J2,J1)=1.0/R2(J2,J1)
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            RR2(J1,J2)=RR2(J2,J1) 
            R2(J1,J2)=R2(J2,J1) 
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
10       CONTINUE 
20    CONTINUE 
C
C  FIRST THE GRADIENT.
 
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               IF (J3.NE.J1) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 95 J4=J3+1,N
                  IF (J4.NE.J1) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1) 
                  RHOAC=(RAC-RE)/RE
                  RHOBC=(RBC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QQ2=(RHOBC-RHOAC)/SR2
                  QB=(QQ2**2)+(QQ3**2)
                  QC=QQ3**3-3*QQ3*QQ2**2
                  TEMP=TEMP+
     1           D*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     2           C10*QC-AN*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3           QA**2*(C2+C8*QB)+C6*QC+QA*( C1+C5*QB+C10*QC)))*
     4           (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     5           (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     6           (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     7           (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     8           (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/RE)/
     9           DEXP(AN*QA)
               ENDIF
95             CONTINUE
            ENDIF
100         CONTINUE
            V(3*(J1-1)+J2)=TEMP
     1                     + V(3*(J1-1)+J2)
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
C
C  DIAGONAL BITS OF THE HESSIAN.
C
      DO 160 J1=1,N
         DO 150 J2=1,3
            TEMP=0.0D0
            DO 140 J3=1,N
               IF (J3.NE.J1) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 130 J4=J3+1,N
                  IF (J4.NE.J1) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
                  RHOBC=(RBC-RE)/RE
                  RHOAC=(RAC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ2=(RHOBC-RHOAC)/SR2
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QB=(QQ2**2)+(QQ3**2)
                  QC=(QQ3**3)-3*QQ3*(QQ2**2)
                  TEMP=TEMP+
     1      D*((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     2      AN**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))*
     4      (-(ABX*RRAB/RE)-ACX*RRAC/RE)**2/SR3**2+
     5      ((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     6      C10*QC-AN*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     7      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))*
     8      (RRAB-ABX**2*RRAB**3+RRAC-ACX**2*RRAC**3)/RE+
     9      (-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     A      (-2*AN*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     B      2*QA*(C2+C8*QB)+C10*QC)*
     C      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     D      (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     E      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     F      +(C6+C10*QA)*
     G      (-6*ACX*QQ2*QQ3*RRAC/SR2+
     H      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     I      )/RE)+2*
     J      ((2*C5+4*C8*QA)*
     K      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     L      C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     M      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/
     N      RE))/SR3+2*C9*
     O      (QB*(-2*QQ2*RRAC/(RE*SR2)+
     P      2*(ACX**2*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     Q      (-2*ABX*RRAB+ACX*RRAC)**2/(RE**2*SR6**2)+
     R      QQ3*(2*RRAB-2*ABX**2*RRAB**3-RRAC+ACX**2*RRAC**3)/
     S      (RE*SR6)))+
     T      4*(ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)**2/
     U      RE**2)+(C3+C5*QA+C8*QA**2)*
     V      (-2*QQ2*RRAC/(RE*SR2)+
     W      2*(ACX**2*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+ 
     X      (-2*ABX*RRAB+ACX*RRAC)**2/(RE**2*SR6**2)+
     Y      QQ3*(2*RRAB-2*ABX**2*RRAB**3-RRAC+ACX**2*RRAC**3)/
     Z      (RE*SR6)))+(C6+C10*QA)*
     A      (QQ3*(-6*ACX**2*(1/(RAC**2*RE**2*SR2**2)+
     B      QQ2*RRAC**3/(RE*SR2))+ 
     C      6*(QQ2*RRAC/(RE*SR2)+ 
     D      (-2*ABX*RRAB+ACX*RRAC)**2/(RE**2*SR6**2)))+
     E      ((-3*QQ2**2+3*QQ3**2)*
     F      (2*RRAB-2*ABX**2*RRAB**3-RRAC+ACX**2*RRAC**3)/RE-
     G      12*ACX*QQ2*RRAC*(-2*ABX*RRAB+ACX*RRAC)/(RE**2*SR2))/SR6))/ 
     H      DEXP(AN*QA)
               ENDIF
130            CONTINUE
            ENDIF
140         CONTINUE
            HESS(3*(J1-1)+J2,3*(J1-1)+J2)=TEMP
     1                                 + HESS(3*(J1-1)+J2,3*(J1-1)+J2)
C           PRINT*,'K2,A=',3*(J1-1)+J2,HESS(3*(J1-1)+J2,3*(J1-1)+J2)
150      CONTINUE
160   CONTINUE
C
C  SAME ATOM, DIFFERENT COMPONENT.
C
      DO 210 J1=1,N
        DO 200 J2=1,3
          DO 190 J5=J2+1,3 
            TEMP=0.0D0
            DO 180 J3=1,N
              IF (J3.NE.J1) THEN
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 170 J4=J3+1,N
                IF (J4.NE.J1) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2) 
                QC=QQ3**3-3*QQ3*QQ2**2
               TEMP=TEMP+ 
     1 D*(((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+C10*QC- 
     2           AN*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3     QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))*
     4         (-ABX*ABY*RRAB**3/RE-ACX*ACY*RRAC**3/RE)+
     5        (-ABY*RRAB/RE-ACY*RRAC/RE)*
     6         (-(AN*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     7        2*QA*(C2+C8*QB)+C10*QC)*
     8      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     9    (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     A     (ACX*QQ2*RRAC/SR2+
     B       QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     C       (C6+C10*QA)*
     D     (-6*ACX*QQ2*QQ3*RRAC/SR2+
     E       (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/
     F        SR6))/RE))+
     G     (2*(C5+2*C8*QA)*
     H      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     I     C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     J     (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/
     K      RE))/SR3+(-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     L   ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     M      AN**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     N      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))* 
     O    (-(ABY*RRAB/RE)-ACY*RRAC/RE)/SR3**2+
     P     (-(AN*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     Q     2*QA*(C2+C8*QB)+C10*QC)*
     R   (-(ABY*RRAB/RE)-ACY*RRAC/RE)/SR3+
     S    (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     T     (ACY*QQ2*RRAC/SR2+
     U    QQ3*(-2*ABY*RRAB+ACY*RRAC)/SR6)+
     V    (C6+C10*QA)*
     W     (-6*ACY*QQ2*QQ3*RRAC/SR2+
     X    (-3*QQ2**2+3*QQ3**2)*(-2*ABY*RRAB+ACY*RRAC)/
     Y     SR6))/RE))+
     Z     (2*(C5+2*C8*QA)*
     A   (ACY*QQ2*RRAC/SR2+QQ3*(-2*ABY*RRAB+ACY*RRAC)/SR6)+
     B  C10*(-6*ACY*QQ2*QQ3*RRAC/SR2+
     C  (-3*QQ2**2+3*QQ3**2)*(-2*ABY*RRAB+ACY*RRAC)/SR6))/
     D      RE)/SR3)+2*(C9*
     E   (4*(ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     F    (ACY*QQ2*RRAC/SR2+QQ3*(-2*ABY*RRAB+ACY*RRAC)/SR6)/RE**2 
     G   +2*QB*(ACX*ACY*
     H   (1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     I  (-2*ABX*RRAB+ACX*RRAC)*(-2*ABY*RRAB+ACY*RRAC)/
     J   (RE**2*SR6**2)+
     K  QQ3*(-2*ABX*ABY*RRAB**3+ACX*ACY*RRAC**3)/(RE*SR6)))+
     L  (C3+C5*QA+C8*QA**2)*
     M   (ACX*ACY*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     N     (-2*ABX*RRAB+ACX*RRAC)*(-2*ABY*RRAB+ACY*RRAC)/
     O   (RE**2*SR6**2)+
     P  QQ3*(-2*ABX*ABY*RRAB**3+ACX*ACY*RRAC**3)/(RE*SR6)))+
     Q     (C6+C10*QA)*(-6*(ACX*
     R   (ACY*QQ3*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     S  QQ2*RRAC*(-2*ABY*RRAB+ACY*RRAC)/(RE**2*SR2*SR6))+
     T  ACY*QQ2*RRAC*(-2*ABX*RRAB+ACX*RRAC)/(RE**2*SR2*SR6))+
     U  6*QQ3*(-2*ABX*RRAB+ACX*RRAC)*(-2*ABY*RRAB+ACY*RRAC)/
     V   (RE**2*SR6**2)+
     W  (-3*QQ2**2+3*QQ3**2)*(-2*ABX*ABY*RRAB**3+ACX*ACY*RRAC**3)/ 
     X   (RE*SR6)))*DEXP(-AN*QA)
               ENDIF
170            CONTINUE
            ENDIF
180         CONTINUE
            HESS(3*(J1-1)+J5,3*(J1-1)+J2)=TEMP
     1                                 + HESS(3*(J1-1)+J5,3*(J1-1)+J2)
190      CONTINUE
200   CONTINUE
210   CONTINUE
C
C  DIFFERENT ATOMS, SAME COMPONENT.
C
      DO 260 J1=1,N
        DO 250 J2=1,3
          DO 230 J3=J1+1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            ABX=VEC(J3,J1,J2)
            TEMP=0.0D0
            DO 220 J4=1,N
              IF ((J4.NE.J3).AND.(J4.NE.J1)) THEN
              BCX=VEC(J4,J3,J2)
              ACX=VEC(J4,J1,J2)
              RBC=R2(J4,J3)
              RRBC=RR2(J4,J3)
              RAC=R2(J4,J1)
              RRAC=RR2(J4,J1)
              RHOAB=(RAB-RE)/RE
              RHOBC=(RBC-RE)/RE
              RHOAC=(RAC-RE)/RE
              QA=(RHOAB+RHOBC+RHOAC)/SR3
              QQ2=(RHOBC-RHOAC)/SR2
              QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
              QB=(QQ2**2)+(QQ3**2)
              QC=QQ3**3-3*QQ3*QQ2**2
              TEMP=TEMP+
     1        D*((C3+C5*QA+C8*QA**2+2*C9*QB)*
     2        (-2*ACX*BCX*RRAC*RRBC/(RE**2*SR2**2)+
     3        2*((-2*ABX*RRAB+ACX*RRAC)*(2*ABX*RRAB+BCX*RRBC)/
     4        (RE**2*SR6**2)+QQ3*(-2*RRAB+2*ABX**2*RRAB**3)/(RE*SR6)))
     5        +(((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     6        C10*QC-AN*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     7        QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))
     8        *(-RRAB+ABX**2*RRAB**3)+
     9        (ABX*RRAB-BCX*RRBC)*
     A        (-(AN*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     B        2*QA*(C2+C8*QB)+C10*QC)*
     C        (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     D        (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     E        (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     F        (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     G        (-3*QQ2**2+3*QQ3**2)*
     H        (-2*ABX*RRAB+ACX*RRAC)/SR6))/RE))+
     I        (2*(C5+2*C8*QA)*
     J        (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     K        C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     L        (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     M        )/RE))/SR3+(-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     N        ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     O        AN**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     P        QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC))
     Q        )*(ABX*RRAB-BCX*RRBC)/SR3**2+
     R        (-(AN*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     S        2*QA*(C2+C8*QB)+C10*QC)*(ABX*RRAB-BCX*RRBC)
     T        /SR3+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     U        (-2*BCX*QQ2*RRBC/SR2+
     V        2*QQ3*(2*ABX*RRAB+BCX*RRBC)/SR6)+
     W        (C6+C10*QA)*(6*BCX*QQ2*QQ3*RRBC/SR2+
     X        (-3*QQ2**2+3*QQ3**2)*(2*ABX*RRAB+BCX*RRBC)/SR6
     Y        )))+(C5+2*C8*QA)*
     Z        (-2*BCX*QQ2*RRBC/SR2+2*QQ3*(2*ABX*RRAB+BCX*RRBC)/SR6)
     A        +C10*(6*BCX*QQ2*QQ3*RRBC/SR2+
     B        (-3*QQ2**2+3*QQ3**2)*(2*ABX*RRAB+BCX*RRBC)/SR6))/
     C        SR3)+4*C9*(ACX*QQ2*RRAC/SR2+
     D        QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     E        (-2*BCX*QQ2*RRBC/SR2+2*QQ3*(2*ABX*RRAB+BCX*RRBC)/SR6)/RE)/
     F        RE+(C6+C10*QA)*(6*
     G        (QQ3*(ACX*BCX*RRAC*RRBC/SR2**2+
     H        (-2*ABX*RRAB+ACX*RRAC)*(2*ABX*RRAB+BCX*RRBC)/SR6**2)
     I        +BCX*QQ2*(-2*ABX*RRAB+ACX*RRAC)*RRBC/(SR2*SR6))/RE**2
     J        +((-3*QQ2**2+3*QQ3**2)*(-2*RRAB+2*ABX**2*RRAB**3)/RE-
     K        6*ACX*QQ2*RRAC*(2*ABX*RRAB+BCX*RRBC)/(RE**2*SR2))/SR6))/
     L        DEXP(AN*QA)

            ENDIF
220         CONTINUE
            HESS(3*(J3-1)+J2,3*(J1-1)+J2)=TEMP
     1                                 + HESS(3*(J3-1)+J2,3*(J1-1)+J2)
230       CONTINUE
250     CONTINUE
260   CONTINUE
C
C  DIFFERENT ATOMS AND DIFFERENT COMPONENTS
C
      DO 310 J1=1,N
        DO 300 J2=1,3
          DO 280 J3=J1+1,N
            DO 290 J5=1,J2-1 
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 270 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((ABX*ABY*(C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*
     2      (C2+C8*QB)+
     3      C10*QC-AN*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     4      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))
     5      *RRAB**3+(ABY*RRAB-BCY*RRBC)*
     6      (-(AN*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     7      2*QA*(C2+C8*QB)+C10*QC)*
     8      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     9      (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     A      (ACX*QQ2*RRAC/SR2+
     B      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     C      (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     D      (-3*QQ2**2+3*QQ3**2)*
     E      (-2*ABX*RRAB+ACX*RRAC)/SR6))/RE))+
     F      (2*(C5+2*C8*QA)*
     G      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     H      +C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     I      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     J      )/RE))/SR3+
     K      (-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     L      ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     M      AN**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     N      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC))
     O      )*(ABY*RRAB-BCY*RRBC)/SR3**2+
     P      (-(AN*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     Q      2*QA*(C2+C8*QB)+C10*QC)*(ABY*RRAB-BCY*RRBC)
     R      /SR3+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     S      (-2*BCY*QQ2*RRBC/SR2+
     T      2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)+
     U      (C6+C10*QA)*(6*BCY*QQ2*QQ3*RRBC/SR2+
     V      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6
     W      )))+(C5+2*C8*QA)*
     X      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)
     Y      +C10*(6*BCY*QQ2*QQ3*RRBC/SR2+
     Z      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6))/
     A      SR3)+4*C9*(ACX*QQ2*RRAC/SR2+
     B      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     C      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)/RE)/ 
     D      RE+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     E      ((-2*ACX*BCY*RRAC*RRBC/SR2**2+
     F      2*(-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     G      RE**2+4*ABX*ABY*QQ3*RRAB**3/(RE*SR6))+
     H      (C6+C10*QA)*(6*(QQ3*
     I      (ACX*BCY*RRAC*RRBC/SR2**2+
     J      (-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     K      RE**2+(ABX*ABY*QQ3**2*RRAB**3+
     L      BCY*QQ2*(-2*ABX*RRAB+ACX*RRAC)*RRBC/(RE*SR2))/(RE*SR6))
     M      -6*(ABX*ABY*QQ2**2*RRAB**3+
     N      ACX*QQ2*RRAC*(2*ABY*RRAB+BCY*RRBC)/(RE*SR2))/(RE*SR6)))/ 
     O      DEXP(AN*QA)
              ENDIF
270           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=TEMP
     1                                   +HESS(3*(J3-1)+J5,3*(J1-1)+J2)
290         CONTINUE
            DO 295 J5=J2+1,3
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 275 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOAC+RHOBC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOAC-RHOBC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((ABX*ABY*(C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*
     2      (C2+C8*QB)+
     3      C10*QC-AN*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     4      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))
     5      *RRAB**3+(ABY*RRAB-BCY*RRBC)*
     6      (-(AN*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     7      2*QA*(C2+C8*QB)+C10*QC)*
     8      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     9      (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     A      (ACX*QQ2*RRAC/SR2+
     B      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     C      (C6+C10*QA)*
     D      (-6*ACX*QQ2*QQ3*RRAC/SR2+
     E      (-3*QQ2**2+3*QQ3**2)*
     F      (-2*ABX*RRAB+ACX*RRAC)/SR6))/RE))+
     G      (2*(C5+2*C8*QA)*
     H      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     I      +C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     J      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     K      )/RE))/SR3+
     L      (-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     M      ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     O      AN**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     P      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC))
     Q      )*(ABY*RRAB-BCY*RRBC)/SR3**2+
     R      (-(AN*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     S      2*QA*(C2+C8*QB)+C10*QC)*(ABY*RRAB-BCY*RRBC)
     T      /SR3+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     U      (-2*BCY*QQ2*RRBC/SR2+
     V      2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)+
     W      (C6+C10*QA)*(6*BCY*QQ2*QQ3*RRBC/SR2+
     X      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6
     Y      )))+(C5+2*C8*QA)*
     Z      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)
     A      +C10*(6*BCY*QQ2*QQ3*RRBC/SR2+
     B      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6))/
     C      SR3)+4*C9*(ACX*QQ2*RRAC/SR2+
     D      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     E      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)/RE)/
     F      RE+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     G      ((-2*ACX*BCY*RRAC*RRBC/SR2**2+
     H      2*(-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     I      RE**2+4*ABX*ABY*QQ3*RRAB**3/(RE*SR6))+
     J      (C6+C10*QA)*(6*(QQ3*
     K      (ACX*BCY*RRAC*RRBC/SR2**2+
     L      (-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     M      RE**2+(ABX*ABY*QQ3**2*RRAB**3+
     N      BCY*QQ2*(-2*ABX*RRAB+ACX*RRAC)*RRBC/(RE*SR2))/(RE*SR6))
     O      -6*(ABX*ABY*QQ2**2*RRAB**3+
     P      ACX*QQ2*RRAC*(2*ABY*RRAB+BCY*RRBC)/(RE*SR2))/(RE*SR6)))/
     Q      DEXP(AN*QA)
              ENDIF
275           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=TEMP
     1                                   + HESS(3*(J3-1)+J5,3*(J1-1)+J2)
295         CONTINUE
280       CONTINUE
300     CONTINUE
310   CONTINUE
C
C  SYMMETRISE
C
      DO 1000 J1=1,3*N
         DO 1010 J2=J1+1,3*N
C           PRINT*,'J1,J2,A=',J1,J2,HESS(J2,J1)
            HESS(J1,J2)=HESS(J2,J1)
1010     CONTINUE
1000  CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  HERE WE CALCULATE THE ANALYTIC GRADIENT AND SECOND DERIVATIVES
C  FOR THE THREE-BODY TERM WITH PERIODIC BOUNDARY CONDITIONS.
C                                        
C*************************************************************************
C
      SUBROUTINE JM3C (N, X, V, CUTOFF, R2, RR2, VEC, RH, M, NDUM, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5
      INTEGER M(N,N,3)
      LOGICAL STEST
      DOUBLE PRECISION R2(N,N), RR2(N,N), 
     1              VEC(N,N,3), RH(N,N), NDUM(N,N)
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION X(3*N), ABX, ACX, BCX, V(3*N),
     1                 TEMP, RRBC,
     2                 RAB, RRAB, RAC, RRAC, RBC,
     3                 ABY, ACY, BCY, RHOAB, RHOBC ,RHOAC ,QA, QQ2,
     4                 QQ3, QB, QC, 
     5                 SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF
C
C
C  FIRST THE GRADIENT.
 
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 95 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF)) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1) 
                  RHOAC=(RAC-RE)/RE
                  RHOBC=(RBC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QQ2=(RHOBC-RHOAC)/SR2
                  QB=(QQ2**2)+(QQ3**2)
                  QC=QQ3**3-3*QQ3*QQ2**2
                  TEMP=TEMP+
     1           D*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     2           C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3           QA**2*(C2+C8*QB)+C6*QC+QA*( C1+C5*QB+C10*QC)))*
     4           (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     5           (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     6           (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     7           (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     8           (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/RE)/
     9           DEXP(AN3*QA)
               ENDIF
95             CONTINUE
            ENDIF
100         CONTINUE
            V(3*(J1-1)+J2)=V(3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
      IF (.NOT.STEST) RETURN
C
C  DIAGONAL BITS OF THE HESSIAN.
C
      DO 160 J1=1,N
         DO 150 J2=1,3
            TEMP=0.0D0
            DO 140 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 130 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF)) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
                  RHOBC=(RBC-RE)/RE
                  RHOAC=(RAC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ2=(RHOBC-RHOAC)/SR2
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QB=(QQ2**2)+(QQ3**2)
                  QC=(QQ3**3)-3*QQ3*(QQ2**2)
                  TEMP=TEMP+
     1      D*((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     2      AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))*
     4      (-(ABX*RRAB/RE)-ACX*RRAC/RE)**2/SR3**2+
     5      ((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     6      C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     7      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))*
     8      (RRAB-ABX**2*RRAB**3+RRAC-ACX**2*RRAC**3)/RE+
     9      (-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     A      (-2*AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     B      2*QA*(C2+C8*QB)+C10*QC)*
     C      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     D      (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     E      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     F      +(C6+C10*QA)*
     G      (-6*ACX*QQ2*QQ3*RRAC/SR2+
     H      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     I      )/RE)+2*
     J      ((2*C5+4*C8*QA)*
     K      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     L      C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     M      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/
     N      RE))/SR3+2*C9*
     O      (QB*(-2*QQ2*RRAC/(RE*SR2)+
     P      2*(ACX**2*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     Q      (-2*ABX*RRAB+ACX*RRAC)**2/(RE**2*SR6**2)+
     R      QQ3*(2*RRAB-2*ABX**2*RRAB**3-RRAC+ACX**2*RRAC**3)/
     S      (RE*SR6)))+
     T      4*(ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)**2/
     U      RE**2)+(C3+C5*QA+C8*QA**2)*
     V      (-2*QQ2*RRAC/(RE*SR2)+
     W      2*(ACX**2*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+ 
     X      (-2*ABX*RRAB+ACX*RRAC)**2/(RE**2*SR6**2)+
     Y      QQ3*(2*RRAB-2*ABX**2*RRAB**3-RRAC+ACX**2*RRAC**3)/
     Z      (RE*SR6)))+(C6+C10*QA)*
     A      (QQ3*(-6*ACX**2*(1/(RAC**2*RE**2*SR2**2)+
     B      QQ2*RRAC**3/(RE*SR2))+ 
     C      6*(QQ2*RRAC/(RE*SR2)+ 
     D      (-2*ABX*RRAB+ACX*RRAC)**2/(RE**2*SR6**2)))+
     E      ((-3*QQ2**2+3*QQ3**2)*
     F      (2*RRAB-2*ABX**2*RRAB**3-RRAC+ACX**2*RRAC**3)/RE-
     G      12*ACX*QQ2*RRAC*(-2*ABX*RRAB+ACX*RRAC)/(RE**2*SR2))/SR6))/ 
     H      DEXP(AN3*QA)
               ENDIF
130            CONTINUE
            ENDIF
140         CONTINUE
            HESS(3*(J1-1)+J2,3*(J1-1)+J2)=HESS(3*(J1-1)+J2,3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,A=',3*(J1-1)+J2,HESS(3*(J1-1)+J2,3*(J1-1)+J2)
150      CONTINUE
160   CONTINUE
C
C  SAME ATOM, DIFFERENT COMPONENT.
C
      DO 210 J1=1,N
        DO 200 J2=1,3
          DO 190 J5=J2+1,3 
            TEMP=0.0D0
            DO 180 J3=1,N
              IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 170 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2) 
                QC=QQ3**3-3*QQ3*QQ2**2
               TEMP=TEMP+ 
     1 D*(((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+C10*QC- 
     2           AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3     QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))*
     4         (-ABX*ABY*RRAB**3/RE-ACX*ACY*RRAC**3/RE)+
     5        (-ABY*RRAB/RE-ACY*RRAC/RE)*
     6         (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     7        2*QA*(C2+C8*QB)+C10*QC)*
     8      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     9    (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     A     (ACX*QQ2*RRAC/SR2+
     B       QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     C       (C6+C10*QA)*
     D     (-6*ACX*QQ2*QQ3*RRAC/SR2+
     E       (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/
     F        SR6))/RE))+
     G     (2*(C5+2*C8*QA)*
     H      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     I     C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     J     (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/
     K      RE))/SR3+(-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     L   ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     M      AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     N      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))* 
     O    (-(ABY*RRAB/RE)-ACY*RRAC/RE)/SR3**2+
     P     (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     Q     2*QA*(C2+C8*QB)+C10*QC)*
     R   (-(ABY*RRAB/RE)-ACY*RRAC/RE)/SR3+
     S    (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     T     (ACY*QQ2*RRAC/SR2+
     U    QQ3*(-2*ABY*RRAB+ACY*RRAC)/SR6)+
     V    (C6+C10*QA)*
     W     (-6*ACY*QQ2*QQ3*RRAC/SR2+
     X    (-3*QQ2**2+3*QQ3**2)*(-2*ABY*RRAB+ACY*RRAC)/
     Y     SR6))/RE))+
     Z     (2*(C5+2*C8*QA)*
     A   (ACY*QQ2*RRAC/SR2+QQ3*(-2*ABY*RRAB+ACY*RRAC)/SR6)+
     B  C10*(-6*ACY*QQ2*QQ3*RRAC/SR2+
     C  (-3*QQ2**2+3*QQ3**2)*(-2*ABY*RRAB+ACY*RRAC)/SR6))/
     D      RE)/SR3)+2*(C9*
     E   (4*(ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     F    (ACY*QQ2*RRAC/SR2+QQ3*(-2*ABY*RRAB+ACY*RRAC)/SR6)/RE**2 
     G   +2*QB*(ACX*ACY*
     H   (1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     I  (-2*ABX*RRAB+ACX*RRAC)*(-2*ABY*RRAB+ACY*RRAC)/
     J   (RE**2*SR6**2)+
     K  QQ3*(-2*ABX*ABY*RRAB**3+ACX*ACY*RRAC**3)/(RE*SR6)))+
     L  (C3+C5*QA+C8*QA**2)*
     M   (ACX*ACY*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     N     (-2*ABX*RRAB+ACX*RRAC)*(-2*ABY*RRAB+ACY*RRAC)/
     O   (RE**2*SR6**2)+
     P  QQ3*(-2*ABX*ABY*RRAB**3+ACX*ACY*RRAC**3)/(RE*SR6)))+
     Q     (C6+C10*QA)*(-6*(ACX*
     R   (ACY*QQ3*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     S  QQ2*RRAC*(-2*ABY*RRAB+ACY*RRAC)/(RE**2*SR2*SR6))+
     T  ACY*QQ2*RRAC*(-2*ABX*RRAB+ACX*RRAC)/(RE**2*SR2*SR6))+
     U  6*QQ3*(-2*ABX*RRAB+ACX*RRAC)*(-2*ABY*RRAB+ACY*RRAC)/
     V   (RE**2*SR6**2)+
     W  (-3*QQ2**2+3*QQ3**2)*(-2*ABX*ABY*RRAB**3+ACX*ACY*RRAC**3)/ 
     X   (RE*SR6)))*DEXP(-AN3*QA)
               ENDIF
170            CONTINUE
            ENDIF
180         CONTINUE
            HESS(3*(J1-1)+J5,3*(J1-1)+J2)=HESS(3*(J1-1)+J5,3*(J1-1)+J2)+TEMP
190      CONTINUE
200   CONTINUE
210   CONTINUE
C
C  DIFFERENT ATOMS, SAME COMPONENT.
C
      DO 260 J1=1,N
        DO 250 J2=1,3
          DO 230 J3=J1+1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            ABX=VEC(J3,J1,J2)
            TEMP=0.0D0
            DO 220 J4=1,N
              IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1            R2(J4,J3).LT.CUTOFF.AND.
     2            R2(J4,J1).LT.CUTOFF.AND.
     3            R2(J3,J1).LT.CUTOFF) THEN
              BCX=VEC(J4,J3,J2)
              ACX=VEC(J4,J1,J2)
              RBC=R2(J4,J3)
              RRBC=RR2(J4,J3)
              RAC=R2(J4,J1)
              RRAC=RR2(J4,J1)
              RHOAB=(RAB-RE)/RE
              RHOBC=(RBC-RE)/RE
              RHOAC=(RAC-RE)/RE
              QA=(RHOAB+RHOBC+RHOAC)/SR3
              QQ2=(RHOBC-RHOAC)/SR2
              QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
              QB=(QQ2**2)+(QQ3**2)
              QC=QQ3**3-3*QQ3*QQ2**2
              TEMP=TEMP+
     1        D*((C3+C5*QA+C8*QA**2+2*C9*QB)*
     2        (-2*ACX*BCX*RRAC*RRBC/(RE**2*SR2**2)+
     3        2*((-2*ABX*RRAB+ACX*RRAC)*(2*ABX*RRAB+BCX*RRBC)/
     4        (RE**2*SR6**2)+QQ3*(-2*RRAB+2*ABX**2*RRAB**3)/(RE*SR6)))
     5        +(((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     6        C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     7        QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))
     8        *(-RRAB+ABX**2*RRAB**3)+
     9        (ABX*RRAB-BCX*RRBC)*
     A        (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     B        2*QA*(C2+C8*QB)+C10*QC)*
     C        (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     D        (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     E        (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     F        (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     G        (-3*QQ2**2+3*QQ3**2)*
     H        (-2*ABX*RRAB+ACX*RRAC)/SR6))/RE))+
     I        (2*(C5+2*C8*QA)*
     J        (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     K        C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     L        (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     M        )/RE))/SR3+(-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     N        ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     O        AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     P        QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC))
     Q        )*(ABX*RRAB-BCX*RRBC)/SR3**2+
     R        (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     S        2*QA*(C2+C8*QB)+C10*QC)*(ABX*RRAB-BCX*RRBC)
     T        /SR3+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     U        (-2*BCX*QQ2*RRBC/SR2+
     V        2*QQ3*(2*ABX*RRAB+BCX*RRBC)/SR6)+
     W        (C6+C10*QA)*(6*BCX*QQ2*QQ3*RRBC/SR2+
     X        (-3*QQ2**2+3*QQ3**2)*(2*ABX*RRAB+BCX*RRBC)/SR6
     Y        )))+(C5+2*C8*QA)*
     Z        (-2*BCX*QQ2*RRBC/SR2+2*QQ3*(2*ABX*RRAB+BCX*RRBC)/SR6)
     A        +C10*(6*BCX*QQ2*QQ3*RRBC/SR2+
     B        (-3*QQ2**2+3*QQ3**2)*(2*ABX*RRAB+BCX*RRBC)/SR6))/
     C        SR3)+4*C9*(ACX*QQ2*RRAC/SR2+
     D        QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     E        (-2*BCX*QQ2*RRBC/SR2+2*QQ3*(2*ABX*RRAB+BCX*RRBC)/SR6)/RE)/
     F        RE+(C6+C10*QA)*(6*
     G        (QQ3*(ACX*BCX*RRAC*RRBC/SR2**2+
     H        (-2*ABX*RRAB+ACX*RRAC)*(2*ABX*RRAB+BCX*RRBC)/SR6**2)
     I        +BCX*QQ2*(-2*ABX*RRAB+ACX*RRAC)*RRBC/(SR2*SR6))/RE**2
     J        +((-3*QQ2**2+3*QQ3**2)*(-2*RRAB+2*ABX**2*RRAB**3)/RE-
     K        6*ACX*QQ2*RRAC*(2*ABX*RRAB+BCX*RRBC)/(RE**2*SR2))/SR6))/
     L        DEXP(AN3*QA)

            ENDIF
220         CONTINUE
            HESS(3*(J3-1)+J2,3*(J1-1)+J2)=HESS(3*(J3-1)+J2,3*(J1-1)+J2)+TEMP
230       CONTINUE
250     CONTINUE
260   CONTINUE
C
C  DIFFERENT ATOMS AND DIFFERENT COMPONENTS
C
      DO 310 J1=1,N
        DO 300 J2=1,3
          DO 280 J3=J1+1,N
            DO 290 J5=1,J2-1 
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 270 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             R2(J4,J3).LT.CUTOFF.AND.
     2             R2(J4,J1).LT.CUTOFF.AND.
     3             R2(J3,J1).LT.CUTOFF) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((ABX*ABY*(C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*
     2      (C2+C8*QB)+
     3      C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     4      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))
     5      *RRAB**3+(ABY*RRAB-BCY*RRBC)*
     6      (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     7      2*QA*(C2+C8*QB)+C10*QC)*
     8      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     9      (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     A      (ACX*QQ2*RRAC/SR2+
     B      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     C      (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     D      (-3*QQ2**2+3*QQ3**2)*
     E      (-2*ABX*RRAB+ACX*RRAC)/SR6))/RE))+
     F      (2*(C5+2*C8*QA)*
     G      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     H      +C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     I      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     J      )/RE))/SR3+
     K      (-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     L      ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     M      AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     N      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC))
     O      )*(ABY*RRAB-BCY*RRBC)/SR3**2+
     P      (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     Q      2*QA*(C2+C8*QB)+C10*QC)*(ABY*RRAB-BCY*RRBC)
     R      /SR3+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     S      (-2*BCY*QQ2*RRBC/SR2+
     T      2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)+
     U      (C6+C10*QA)*(6*BCY*QQ2*QQ3*RRBC/SR2+
     V      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6
     W      )))+(C5+2*C8*QA)*
     X      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)
     Y      +C10*(6*BCY*QQ2*QQ3*RRBC/SR2+
     Z      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6))/
     A      SR3)+4*C9*(ACX*QQ2*RRAC/SR2+
     B      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     C      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)/RE)/ 
     D      RE+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     E      ((-2*ACX*BCY*RRAC*RRBC/SR2**2+
     F      2*(-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     G      RE**2+4*ABX*ABY*QQ3*RRAB**3/(RE*SR6))+
     H      (C6+C10*QA)*(6*(QQ3*
     I      (ACX*BCY*RRAC*RRBC/SR2**2+
     J      (-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     K      RE**2+(ABX*ABY*QQ3**2*RRAB**3+
     L      BCY*QQ2*(-2*ABX*RRAB+ACX*RRAC)*RRBC/(RE*SR2))/(RE*SR6))
     M      -6*(ABX*ABY*QQ2**2*RRAB**3+
     N      ACX*QQ2*RRAC*(2*ABY*RRAB+BCY*RRBC)/(RE*SR2))/(RE*SR6)))/ 
     O      DEXP(AN3*QA)
              ENDIF
270           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
290         CONTINUE
            DO 295 J5=J2+1,3
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 275 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             R2(J4,J3).LT.CUTOFF.AND.
     2             R2(J4,J1).LT.CUTOFF.AND.
     3             R2(J3,J1).LT.CUTOFF) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOAC+RHOBC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOAC-RHOBC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((ABX*ABY*(C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*
     2      (C2+C8*QB)+
     3      C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     4      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))
     5      *RRAB**3+(ABY*RRAB-BCY*RRBC)*
     6      (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     7      2*QA*(C2+C8*QB)+C10*QC)*
     8      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     9      (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     A      (ACX*QQ2*RRAC/SR2+
     B      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     C      (C6+C10*QA)*
     D      (-6*ACX*QQ2*QQ3*RRAC/SR2+
     E      (-3*QQ2**2+3*QQ3**2)*
     F      (-2*ABX*RRAB+ACX*RRAC)/SR6))/RE))+
     G      (2*(C5+2*C8*QA)*
     H      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     I      +C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     J      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     K      )/RE))/SR3+
     L      (-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     M      ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     O      AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     P      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC))
     Q      )*(ABY*RRAB-BCY*RRBC)/SR3**2+
     R      (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     S      2*QA*(C2+C8*QB)+C10*QC)*(ABY*RRAB-BCY*RRBC)
     T      /SR3+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     U      (-2*BCY*QQ2*RRBC/SR2+
     V      2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)+
     W      (C6+C10*QA)*(6*BCY*QQ2*QQ3*RRBC/SR2+
     X      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6
     Y      )))+(C5+2*C8*QA)*
     Z      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)
     A      +C10*(6*BCY*QQ2*QQ3*RRBC/SR2+
     B      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6))/
     C      SR3)+4*C9*(ACX*QQ2*RRAC/SR2+
     D      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     E      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)/RE)/
     F      RE+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     G      ((-2*ACX*BCY*RRAC*RRBC/SR2**2+
     H      2*(-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     I      RE**2+4*ABX*ABY*QQ3*RRAB**3/(RE*SR6))+
     J      (C6+C10*QA)*(6*(QQ3*
     K      (ACX*BCY*RRAC*RRBC/SR2**2+
     L      (-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     M      RE**2+(ABX*ABY*QQ3**2*RRAB**3+
     N      BCY*QQ2*(-2*ABX*RRAB+ACX*RRAC)*RRBC/(RE*SR2))/(RE*SR6))
     O      -6*(ABX*ABY*QQ2**2*RRAB**3+
     P      ACX*QQ2*RRAC*(2*ABY*RRAB+BCY*RRBC)/(RE*SR2))/(RE*SR6)))/
     Q      DEXP(AN3*QA)
              ENDIF
275           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
295         CONTINUE
280       CONTINUE
300     CONTINUE
310   CONTINUE
C
C  SYMMETRISE
C
      DO 1000 J1=1,3*N
         DO 1010 J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
1010     CONTINUE
1000  CONTINUE
      RETURN
      END

C*************************************************************************
C
C  HERE WE CALCULATE THE ANALYTIC GRADIENT AND SECOND DERIVATIVES
C  FOR THE THREE-BODY TERM 
C                                        
C*************************************************************************
C
      SUBROUTINE JM3CC(N, X, V,R2,RR2,VEC,RH,M,NDUM)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5
      INTEGER M(N,N,3)
      DOUBLE PRECISION R2(N,N), RR2(N,N), C,
     1              VEC(N,N,3), RH(N,N), 
     2              NDUM(N,N)
      DOUBLE PRECISION X(3*N), ABX, ACX, BCX, V(3*N),
     1                 TEMP, RRBC, V0, LAMBDA, SR2, SR3, SR6,
     2                 RAB, RRAB, RAC, RRAC, RBC,
     3                 ABY, ACY, BCY, RHOAB, RHOBC ,RHOAC, QQ2,
     4                 QQ3, QC
      DOUBLE PRECISION CUTOFF
      PARAMETER (SR2=1.414213562D0, SR3=1.732050808D0,SR6=2.449489743D0,
     1           C=1.0D0, 
     2           V0=60.0D0, LAMBDA=-4.0D0, CUTOFF=2.8D0)
C
C
C  FIRST THE GRADIENT.
C 
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=RH(J3,J1)
               DO 95 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(J4.NE.J3).AND.(R2(J4,J3).LT.CUTOFF)) THEN 
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1) 
                  RHOAC=RH(J4,J1)
                  RHOBC=RH(J4,J3)
                  QQ3=(2*RAB-RBC-RAC)/SR6
                  QQ2=(RBC-RAC)/SR2
                  QC=QQ3**3-3.0D0*QQ3*QQ2**2
                  TEMP=TEMP+
     1   DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2   (C*(1.0D0+LAMBDA*QC)
     3    *(ABX*RRAB/RHOAB**2+ACX*RRAC/RHOAC**2) + 
     4   LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     5   (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6))*V0
               ENDIF
95             CONTINUE
            ENDIF
100         CONTINUE
            V(3*(J1-1)+J2)=V(3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
C
C  DIAGONAL BITS OF THE HESSIAN.
C
      DO 160 J1=1,N
         DO 150 J2=1,3
            TEMP=0.0D0
            DO 140 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=RH(J3,J1)
               DO 130 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(J4.NE.J3).AND.(R2(J4,J3).LT.CUTOFF)) THEN 
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
                  RHOBC=RH(J4,J3)
                  RHOAC=RH(J4,J1)
                  QQ2=(RBC-RAC)/SR2
                  QQ3=(2*RAB-RBC-RAC)/SR6
                  QC=(QQ3**3)-3.0D0*QQ3*(QQ2**2)
                  TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (C**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)**2 + 
     4  C*(2*(ABX**2/(RAB**2*RHOAB**3) + ACX**2/(RAC**2*RHOAC**3))+  
     5  (-RRAB + ABX**2*RRAB**3)/RHOAB**2 + 
     6  (-RRAC + ACX**2*RRAC**3)/RHOAC**2)) + 
     7  2*C*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     8  (LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     A  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     F  LAMBDA*(QQ3*(-6*ACX**2*(1/(RAC**2*SR2**2) + QQ2*RRAC**3/SR2) +   
     G  6*(QQ2*RRAC/SR2 + (-2*ABX*RRAB + ACX*RRAC)**2/SR6**2)) + 
     H  ((-3*QQ2**2 + 3*QQ3**2)*
     I  (2*RRAB - 2*ABX**2*RRAB**3 - RRAC + ACX**2*RRAC**3) - 
     J  12*ACX*QQ2*RRAC*(-2*ABX*RRAB + ACX*RRAC)/SR2)/SR6))*V0
               ENDIF
130            CONTINUE
            ENDIF
140         CONTINUE
            HESS(3*(J1-1)+J2,3*(J1-1)+J2)=HESS(3*(J1-1)+J2,3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,A=',3*(J1-1)+J2,HESS(3*(J1-1)+J2,3*(J1-1)+J2)
150      CONTINUE
160   CONTINUE
C
C  SAME ATOM, DIFFERENT COMPONENT.
C
      DO 210 J1=1,N
        DO 200 J2=1,3
          DO 190 J5=J2+1,3 
            TEMP=0.0D0
            DO 180 J3=1,N
              IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 170 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(J4.NE.J3).AND.(R2(J4,J3).LT.CUTOFF)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=RH(J3,J1)
                RHOBC=RH(J4,J3)
                RHOAC=RH(J4,J1)
                QQ2=(RBC-RAC)/SR2
                QQ3=(2*RAB-RBC-RAC)/SR6
                QC=QQ3**3-3.0D0*QQ3*QQ2**2
               TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (C**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     4  (ABY*RRAB/RHOAB**2 + ACY*RRAC/RHOAC**2) + 
     5  C*(2*(ABX*ABY/(RAB**2*RHOAB**3) + ACX*ACY/(RAC**2*RHOAC**3)) + 
     6  ABX*ABY*RRAB**3/RHOAB**2 + ACX*ACY*RRAC**3/RHOAC**2)) + 
     7  C*((ABY*RRAB/RHOAB**2 + ACY*RRAC/RHOAC**2)*
     8  (  LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     A  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     B  (ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     D  (LAMBDA*(-6*ACY*QQ2*QQ3*RRAC/SR2 + 
     E  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABY*RRAB + ACY*RRAC)/SR6))) + 
     I  LAMBDA*(-6*(ACX*(ACY*QQ3*(1/(RAC**2*SR2**2) + QQ2*RRAC**3/SR2)+ 
     J  QQ2*RRAC*(-2*ABY*RRAB + ACY*RRAC)/(SR2*SR6)) + 
     K  ACY*QQ2*RRAC*(-2*ABX*RRAB + ACX*RRAC)/(SR2*SR6)) + 
     L  6*QQ3*(-2*ABX*RRAB + ACX*RRAC)*(-2*ABY*RRAB + ACY*RRAC)/SR6**2 + 
     M  (-3*QQ2**2+3*QQ3**2)*(-2*ABX*ABY*RRAB**3 + ACX*ACY*RRAC**3)/SR6 
     N  ))*V0
               ENDIF
170            CONTINUE
            ENDIF
180         CONTINUE
            HESS(3*(J1-1)+J5,3*(J1-1)+J2)=HESS(3*(J1-1)+J5,3*(J1-1)+J2)+TEMP
190      CONTINUE
200   CONTINUE
210   CONTINUE
C
C  DIFFERENT ATOMS, SAME COMPONENT.
C
      DO 260 J1=1,N
        DO 250 J2=1,3
          DO 230 J3=J1+1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            ABX=VEC(J3,J1,J2)
            TEMP=0.0D0
            DO 220 J4=1,N
              IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1            R2(J4,J3).LT.CUTOFF.AND.
     2            R2(J4,J1).LT.CUTOFF.AND.
     3            R2(J3,J1).LT.CUTOFF) THEN
              BCX=VEC(J4,J3,J2)
              ACX=VEC(J4,J1,J2)
              RBC=R2(J4,J3)
              RRBC=RR2(J4,J3)
              RAC=R2(J4,J1)
              RRAC=RR2(J4,J1)
              RHOAB=RH(J3,J1)
              RHOBC=RH(J4,J3)
              RHOAC=RH(J4,J1)
              QQ2=(RBC-RAC)/SR2
              QQ3=(2*RAB-RBC-RAC)/SR6
              QC=QQ3**3-3.0D0*QQ3*QQ2**2
              TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (C*(RRAB/RHOAB**2 + ABX**2*
     4  (-2/(RAB**2*RHOAB**3) - RRAB**3/RHOAB**2)) + 
     5  C**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     6  (-(ABX*RRAB/RHOAB**2) + BCX*RRBC/RHOBC**2)) + 
     7  C*((-(ABX*RRAB/RHOAB**2) + BCX*RRBC/RHOBC**2)*
     8  (  LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     A  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     B  (ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     D  (LAMBDA*(6*BCX*QQ2*QQ3*RRBC/SR2 + 
     E  (-3*QQ2**2 + 3*QQ3**2)*(2*ABX*RRAB + BCX*RRBC)/SR6))) + 
     I  LAMBDA*(6*(QQ3*(ACX*BCX*RRAC*RRBC/SR2**2 + 
     J  (-2*ABX*RRAB + ACX*RRAC)*(2*ABX*RRAB + BCX*RRBC)/SR6**2) +  
     K  BCX*QQ2*(-2*ABX*RRAB + ACX*RRAC)*RRBC/(SR2*SR6)) + 
     L  ((-3*QQ2**2 + 3*QQ3**2)*(-2*RRAB + 2*ABX**2*RRAB**3) - 
     M  6*ACX*QQ2*RRAC*(2*ABX*RRAB + BCX*RRBC)/SR2)/SR6))*V0
            ENDIF
220         CONTINUE
            HESS(3*(J3-1)+J2,3*(J1-1)+J2)=HESS(3*(J3-1)+J2,3*(J1-1)+J2)+TEMP
230       CONTINUE
250     CONTINUE
260   CONTINUE
C
C  DIFFERENT ATOMS AND DIFFERENT COMPONENTS
C
      DO 310 J1=1,N
        DO 300 J2=1,3
          DO 280 J3=J1+1,N
            DO 290 J5=1,J2-1
              TEMP=0.0D0
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 270 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             (R2(J3,J1).LT.CUTOFF).AND.
     2             (R2(J4,J3).LT.CUTOFF).AND.
     3             (R2(J4,J1).LT.CUTOFF)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=RH(J3,J1)
                RHOBC=RH(J4,J3)
                RHOAC=RH(J4,J1)
                QQ2=(RBC-RAC)/SR2
                QQ3=(2*RAB-RBC-RAC)/SR6
                QC=QQ3**3-3.0D0*QQ3*QQ2**2
                TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (ABX*ABY*C*(-2/(RAB**2*RHOAB**3) - RRAB**3/RHOAB**2) + 
     4  C**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     5  (-(ABY*RRAB/RHOAB**2) + BCY*RRBC/RHOBC**2)) + 
     6  C*((-(ABY*RRAB/RHOAB**2) + BCY*RRBC/RHOBC**2)*
     8  (LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     9  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     A  (ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     C  (LAMBDA*(6*BCY*QQ2*QQ3*RRBC/SR2 + 
     D  (-3*QQ2**2 + 3*QQ3**2)*(2*ABY*RRAB + BCY*RRBC)/SR6))) + 
     H  LAMBDA*(6*(QQ3*(ACX*BCY*RRAC*RRBC/SR2**2 + 
     I  (-2*ABX*RRAB + ACX*RRAC)*(2*ABY*RRAB + BCY*RRBC)/SR6**2) + 
     J  (ABX*ABY*QQ3**2*RRAB**3 + 
     K  BCY*QQ2*(-2*ABX*RRAB + ACX*RRAC)*RRBC/SR2)/SR6) - 
     L  6*(ABX*ABY*QQ2**2*RRAB**3 + 
     M  ACX*QQ2*RRAC*(2*ABY*RRAB + BCY*RRBC)/SR2)/SR6))*V0
              ENDIF
270           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP 
290         CONTINUE
            DO 295 J5=J2+1,3
              TEMP=0.0D0
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 275 J4=1,N
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             (R2(J3,J1).LT.CUTOFF).AND.
     2             (R2(J4,J3).LT.CUTOFF).AND.
     3             (R2(J4,J1).LT.CUTOFF)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=RH(J3,J1)
                RHOBC=RH(J4,J3)
                RHOAC=RH(J4,J1)
                QQ2=(RBC-RAC)/SR2
                QQ3=(2*RAB-RAC-RBC)/SR6
                QC=QQ3**3-3.0D0*QQ3*QQ2**2
                TEMP=TEMP+
     1  DEXP(C*(1/RHOAB + 1/RHOAC + 1/RHOBC))*
     2  ((1 + LAMBDA*QC)*
     3  (ABX*ABY*C*(-2/(RAB**2*RHOAB**3) - RRAB**3/RHOAB**2) + 
     4  C**2*(ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     5  (-(ABY*RRAB/RHOAB**2) + BCY*RRBC/RHOBC**2)) + 
     6  C*((-(ABY*RRAB/RHOAB**2) + BCY*RRBC/RHOBC**2)*
     8  (LAMBDA*(-6*ACX*QQ2*QQ3*RRAC/SR2 + 
     9  (-3*QQ2**2 + 3*QQ3**2)*(-2*ABX*RRAB + ACX*RRAC)/SR6)) + 
     A  (ABX*RRAB/RHOAB**2 + ACX*RRAC/RHOAC**2)*
     C  (LAMBDA*(6*BCY*QQ2*QQ3*RRBC/SR2 + 
     D  (-3*QQ2**2 + 3*QQ3**2)*(2*ABY*RRAB + BCY*RRBC)/SR6))) + 
     H  LAMBDA*(6*(QQ3*(ACX*BCY*RRAC*RRBC/SR2**2 + 
     I  (-2*ABX*RRAB + ACX*RRAC)*(2*ABY*RRAB + BCY*RRBC)/SR6**2) + 
     J  (ABX*ABY*QQ3**2*RRAB**3 + 
     K  BCY*QQ2*(-2*ABX*RRAB + ACX*RRAC)*RRBC/SR2)/SR6) - 
     L  6*(ABX*ABY*QQ2**2*RRAB**3 + 
     M  ACX*QQ2*RRAC*(2*ABY*RRAB + BCY*RRBC)/SR2)/SR6))*V0
              ENDIF
275           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
295         CONTINUE
280       CONTINUE
300     CONTINUE
310   CONTINUE
C
C  SYMMETRISE
C
      DO 1000 J1=1,3*N
         DO 1010 J2=J1+1,3*N
C           PRINT*,'J1,J2,A=',J1,J2,HESS(J2,J1)
            HESS(J1,J2)=HESS(J2,J1)
1010     CONTINUE
1000  CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  HERE WE CALCULATE THE ANALYTIC GRADIENT AND SECOND DERIVATIVES
C  FOR THE THREE-BODY TERM WITH PERIODIC BOUNDARY CONDITIONS.
C                                        
C*************************************************************************
C
      SUBROUTINE JM3P (N, X, V, BOXLX, BOXLY, BOXLZ, CUTOFF, R2, RR2, VEC, M)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, NMX, NMY, NMZ, J3, J4, J5
      INTEGER M(N,N,3)
      DOUBLE PRECISION R2(N,N), RR2(N,N),
     1              VEC(N,N,3)
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION X(3*N), ABX, ACX, BCX, V(3*N),
     1                 TEMP, RRBC,
     2                 RAB, RRAB, RAC, RRAC, RBC,
     3                 ABY, ACY, BCY, RHOAB, RHOBC ,RHOAC ,QA, QQ2,
     4                 QQ3, QB, QC, 
     5                 SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3, BOXLX, BOXLY, BOXLZ
      DOUBLE PRECISION CUTOFF
C
C
C  FIRST THE GRADIENT.
 
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 95 J4=J3+1,N
                  NMX=M(J4,J3,1)+M(J3,J1,1)
                  NMY=M(J4,J3,2)+M(J3,J1,2)
                  NMZ=M(J4,J3,3)+M(J3,J1,3)
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF).AND.
     3                NMZ.EQ.M(J4,J1,3).AND.
     2                NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1) 
                  RHOAC=(RAC-RE)/RE
                  RHOBC=(RBC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QQ2=(RHOBC-RHOAC)/SR2
                  QB=(QQ2**2)+(QQ3**2)
                  QC=QQ3**3-3*QQ3*QQ2**2
                  TEMP=TEMP+
     1           D*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     2           C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3           QA**2*(C2+C8*QB)+C6*QC+QA*( C1+C5*QB+C10*QC)))*
     4           (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     5           (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     6           (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     7           (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     8           (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/RE)/
     9           DEXP(AN3*QA)
               ENDIF
95             CONTINUE
            ENDIF
100         CONTINUE
            V(3*(J1-1)+J2)=V(3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
C
C  DIAGONAL BITS OF THE HESSIAN.
C
      DO 160 J1=1,N
         DO 150 J2=1,3
            TEMP=0.0D0
            DO 140 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 130 J4=J3+1,N
                  NMX=M(J4,J3,1)+M(J3,J1,1)
                  NMY=M(J4,J3,2)+M(J3,J1,2)
                  NMZ=M(J4,J3,3)+M(J3,J1,3)
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF).AND.
     3                NMZ.EQ.M(J4,J1,3).AND.
     2                NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
                  RHOBC=(RBC-RE)/RE
                  RHOAC=(RAC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ2=(RHOBC-RHOAC)/SR2
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QB=(QQ2**2)+(QQ3**2)
                  QC=(QQ3**3)-3*QQ3*(QQ2**2)
                  TEMP=TEMP+
     1      D*((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     2      AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))*
     4      (-(ABX*RRAB/RE)-ACX*RRAC/RE)**2/SR3**2+
     5      ((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     6      C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     7      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))*
     8      (RRAB-ABX**2*RRAB**3+RRAC-ACX**2*RRAC**3)/RE+
     9      (-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     A      (-2*AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     B      2*QA*(C2+C8*QB)+C10*QC)*
     C      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     D      (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     E      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     F      +(C6+C10*QA)*
     G      (-6*ACX*QQ2*QQ3*RRAC/SR2+
     H      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     I      )/RE)+2*
     J      ((2*C5+4*C8*QA)*
     K      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     L      C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     M      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/
     N      RE))/SR3+2*C9*
     O      (QB*(-2*QQ2*RRAC/(RE*SR2)+
     P      2*(ACX**2*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     Q      (-2*ABX*RRAB+ACX*RRAC)**2/(RE**2*SR6**2)+
     R      QQ3*(2*RRAB-2*ABX**2*RRAB**3-RRAC+ACX**2*RRAC**3)/
     S      (RE*SR6)))+
     T      4*(ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)**2/
     U      RE**2)+(C3+C5*QA+C8*QA**2)*
     V      (-2*QQ2*RRAC/(RE*SR2)+
     W      2*(ACX**2*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+ 
     X      (-2*ABX*RRAB+ACX*RRAC)**2/(RE**2*SR6**2)+
     Y      QQ3*(2*RRAB-2*ABX**2*RRAB**3-RRAC+ACX**2*RRAC**3)/
     Z      (RE*SR6)))+(C6+C10*QA)*
     A      (QQ3*(-6*ACX**2*(1/(RAC**2*RE**2*SR2**2)+
     B      QQ2*RRAC**3/(RE*SR2))+ 
     C      6*(QQ2*RRAC/(RE*SR2)+ 
     D      (-2*ABX*RRAB+ACX*RRAC)**2/(RE**2*SR6**2)))+
     E      ((-3*QQ2**2+3*QQ3**2)*
     F      (2*RRAB-2*ABX**2*RRAB**3-RRAC+ACX**2*RRAC**3)/RE-
     G      12*ACX*QQ2*RRAC*(-2*ABX*RRAB+ACX*RRAC)/(RE**2*SR2))/SR6))/ 
     H      DEXP(AN3*QA)
               ENDIF
130            CONTINUE
            ENDIF
140         CONTINUE
            HESS(3*(J1-1)+J2,3*(J1-1)+J2)=HESS(3*(J1-1)+J2,3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,A=',3*(J1-1)+J2,HESS(3*(J1-1)+J2,3*(J1-1)+J2)
150      CONTINUE
160   CONTINUE
C
C  SAME ATOM, DIFFERENT COMPONENT.
C
      DO 210 J1=1,N
        DO 200 J2=1,3
          DO 190 J5=J2+1,3 
            TEMP=0.0D0
            DO 180 J3=1,N
              IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              DO 170 J4=J3+1,N
                  NMX=M(J4,J3,1)+M(J3,J1,1)
                  NMY=M(J4,J3,2)+M(J3,J1,2)
                  NMZ=M(J4,J3,3)+M(J3,J1,3)
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF).AND.
     3                NMZ.EQ.M(J4,J1,3).AND.
     2                NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2) 
                QC=QQ3**3-3*QQ3*QQ2**2
               TEMP=TEMP+ 
     1 D*(((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+C10*QC- 
     2           AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3     QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))*
     4         (-ABX*ABY*RRAB**3/RE-ACX*ACY*RRAC**3/RE)+
     5        (-ABY*RRAB/RE-ACY*RRAC/RE)*
     6         (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     7        2*QA*(C2+C8*QB)+C10*QC)*
     8      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     9    (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     A     (ACX*QQ2*RRAC/SR2+
     B       QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     C       (C6+C10*QA)*
     D     (-6*ACX*QQ2*QQ3*RRAC/SR2+
     E       (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/
     F        SR6))/RE))+
     G     (2*(C5+2*C8*QA)*
     H      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     I     C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     J     (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/
     K      RE))/SR3+(-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     L   ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     M      AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     N      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))* 
     O    (-(ABY*RRAB/RE)-ACY*RRAC/RE)/SR3**2+
     P     (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     Q     2*QA*(C2+C8*QB)+C10*QC)*
     R   (-(ABY*RRAB/RE)-ACY*RRAC/RE)/SR3+
     S    (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     T     (ACY*QQ2*RRAC/SR2+
     U    QQ3*(-2*ABY*RRAB+ACY*RRAC)/SR6)+
     V    (C6+C10*QA)*
     W     (-6*ACY*QQ2*QQ3*RRAC/SR2+
     X    (-3*QQ2**2+3*QQ3**2)*(-2*ABY*RRAB+ACY*RRAC)/
     Y     SR6))/RE))+
     Z     (2*(C5+2*C8*QA)*
     A   (ACY*QQ2*RRAC/SR2+QQ3*(-2*ABY*RRAB+ACY*RRAC)/SR6)+
     B  C10*(-6*ACY*QQ2*QQ3*RRAC/SR2+
     C  (-3*QQ2**2+3*QQ3**2)*(-2*ABY*RRAB+ACY*RRAC)/SR6))/
     D      RE)/SR3)+2*(C9*
     E   (4*(ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     F    (ACY*QQ2*RRAC/SR2+QQ3*(-2*ABY*RRAB+ACY*RRAC)/SR6)/RE**2 
     G   +2*QB*(ACX*ACY*
     H   (1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     I  (-2*ABX*RRAB+ACX*RRAC)*(-2*ABY*RRAB+ACY*RRAC)/
     J   (RE**2*SR6**2)+
     K  QQ3*(-2*ABX*ABY*RRAB**3+ACX*ACY*RRAC**3)/(RE*SR6)))+
     L  (C3+C5*QA+C8*QA**2)*
     M   (ACX*ACY*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     N     (-2*ABX*RRAB+ACX*RRAC)*(-2*ABY*RRAB+ACY*RRAC)/
     O   (RE**2*SR6**2)+
     P  QQ3*(-2*ABX*ABY*RRAB**3+ACX*ACY*RRAC**3)/(RE*SR6)))+
     Q     (C6+C10*QA)*(-6*(ACX*
     R   (ACY*QQ3*(1/(RAC**2*RE**2*SR2**2)+QQ2*RRAC**3/(RE*SR2))+
     S  QQ2*RRAC*(-2*ABY*RRAB+ACY*RRAC)/(RE**2*SR2*SR6))+
     T  ACY*QQ2*RRAC*(-2*ABX*RRAB+ACX*RRAC)/(RE**2*SR2*SR6))+
     U  6*QQ3*(-2*ABX*RRAB+ACX*RRAC)*(-2*ABY*RRAB+ACY*RRAC)/
     V   (RE**2*SR6**2)+
     W  (-3*QQ2**2+3*QQ3**2)*(-2*ABX*ABY*RRAB**3+ACX*ACY*RRAC**3)/ 
     X   (RE*SR6)))*DEXP(-AN3*QA)
               ENDIF
170            CONTINUE
            ENDIF
180         CONTINUE
            HESS(3*(J1-1)+J5,3*(J1-1)+J2)=HESS(3*(J1-1)+J5,3*(J1-1)+J2)+TEMP
190      CONTINUE
200   CONTINUE
210   CONTINUE
C
C  DIFFERENT ATOMS, SAME COMPONENT.
C
      DO 260 J1=1,N
        DO 250 J2=1,3
          DO 230 J3=J1+1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            ABX=VEC(J3,J1,J2)
            TEMP=0.0D0
            DO 220 J4=1,N
               NMX=M(J4,J3,1)+M(J3,J1,1)
               NMY=M(J4,J3,2)+M(J3,J1,2)
               NMZ=M(J4,J3,3)+M(J3,J1,3)
              IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1            R2(J4,J3).LT.CUTOFF.AND.
     2            R2(J4,J1).LT.CUTOFF.AND.
     3            R2(J3,J1).LT.CUTOFF.AND.
     3            NMZ.EQ.M(J4,J1,3).AND.
     4            NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
              BCX=VEC(J4,J3,J2)
              ACX=VEC(J4,J1,J2)
              RBC=R2(J4,J3)
              RRBC=RR2(J4,J3)
              RAC=R2(J4,J1)
              RRAC=RR2(J4,J1)
              RHOAB=(RAB-RE)/RE
              RHOBC=(RBC-RE)/RE
              RHOAC=(RAC-RE)/RE
              QA=(RHOAB+RHOBC+RHOAC)/SR3
              QQ2=(RHOBC-RHOAC)/SR2
              QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
              QB=(QQ2**2)+(QQ3**2)
              QC=QQ3**3-3*QQ3*QQ2**2
              TEMP=TEMP+
     1        D*((C3+C5*QA+C8*QA**2+2*C9*QB)*
     2        (-2*ACX*BCX*RRAC*RRBC/(RE**2*SR2**2)+
     3        2*((-2*ABX*RRAB+ACX*RRAC)*(2*ABX*RRAB+BCX*RRBC)/
     4        (RE**2*SR6**2)+QQ3*(-2*RRAB+2*ABX**2*RRAB**3)/(RE*SR6)))
     5        +(((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     6        C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     7        QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))
     8        *(-RRAB+ABX**2*RRAB**3)+
     9        (ABX*RRAB-BCX*RRBC)*
     A        (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     B        2*QA*(C2+C8*QB)+C10*QC)*
     C        (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     D        (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     E        (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     F        (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     G        (-3*QQ2**2+3*QQ3**2)*
     H        (-2*ABX*RRAB+ACX*RRAC)/SR6))/RE))+
     I        (2*(C5+2*C8*QA)*
     J        (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     K        C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     L        (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     M        )/RE))/SR3+(-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     N        ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     O        AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     P        QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC))
     Q        )*(ABX*RRAB-BCX*RRBC)/SR3**2+
     R        (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     S        2*QA*(C2+C8*QB)+C10*QC)*(ABX*RRAB-BCX*RRBC)
     T        /SR3+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     U        (-2*BCX*QQ2*RRBC/SR2+
     V        2*QQ3*(2*ABX*RRAB+BCX*RRBC)/SR6)+
     W        (C6+C10*QA)*(6*BCX*QQ2*QQ3*RRBC/SR2+
     X        (-3*QQ2**2+3*QQ3**2)*(2*ABX*RRAB+BCX*RRBC)/SR6
     Y        )))+(C5+2*C8*QA)*
     Z        (-2*BCX*QQ2*RRBC/SR2+2*QQ3*(2*ABX*RRAB+BCX*RRBC)/SR6)
     A        +C10*(6*BCX*QQ2*QQ3*RRBC/SR2+
     B        (-3*QQ2**2+3*QQ3**2)*(2*ABX*RRAB+BCX*RRBC)/SR6))/
     C        SR3)+4*C9*(ACX*QQ2*RRAC/SR2+
     D        QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     E        (-2*BCX*QQ2*RRBC/SR2+2*QQ3*(2*ABX*RRAB+BCX*RRBC)/SR6)/RE)/
     F        RE+(C6+C10*QA)*(6*
     G        (QQ3*(ACX*BCX*RRAC*RRBC/SR2**2+
     H        (-2*ABX*RRAB+ACX*RRAC)*(2*ABX*RRAB+BCX*RRBC)/SR6**2)
     I        +BCX*QQ2*(-2*ABX*RRAB+ACX*RRAC)*RRBC/(SR2*SR6))/RE**2
     J        +((-3*QQ2**2+3*QQ3**2)*(-2*RRAB+2*ABX**2*RRAB**3)/RE-
     K        6*ACX*QQ2*RRAC*(2*ABX*RRAB+BCX*RRBC)/(RE**2*SR2))/SR6))/
     L        DEXP(AN3*QA)

            ENDIF
220         CONTINUE
            HESS(3*(J3-1)+J2,3*(J1-1)+J2)=HESS(3*(J3-1)+J2,3*(J1-1)+J2)+TEMP
230       CONTINUE
250     CONTINUE
260   CONTINUE
C
C  DIFFERENT ATOMS AND DIFFERENT COMPONENTS
C
      DO 310 J1=1,N
        DO 300 J2=1,3
          DO 280 J3=J1+1,N
            DO 290 J5=1,J2-1 
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 270 J4=1,N
                NMX=M(J4,J3,1)+M(J3,J1,1)
                NMY=M(J4,J3,2)+M(J3,J1,2)
                NMZ=M(J4,J3,3)+M(J3,J1,3)
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             R2(J4,J3).LT.CUTOFF.AND.
     2             R2(J4,J1).LT.CUTOFF.AND.
     3             R2(J3,J1).LT.CUTOFF.AND.
     3             NMZ.EQ.M(J4,J1,3).AND.
     4             NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOBC+RHOAC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((ABX*ABY*(C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*
     2      (C2+C8*QB)+
     3      C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     4      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))
     5      *RRAB**3+(ABY*RRAB-BCY*RRBC)*
     6      (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     7      2*QA*(C2+C8*QB)+C10*QC)*
     8      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     9      (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     A      (ACX*QQ2*RRAC/SR2+
     B      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     C      (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     D      (-3*QQ2**2+3*QQ3**2)*
     E      (-2*ABX*RRAB+ACX*RRAC)/SR6))/RE))+
     F      (2*(C5+2*C8*QA)*
     G      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     H      +C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     I      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     J      )/RE))/SR3+
     K      (-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     L      ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     M      AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     N      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC))
     O      )*(ABY*RRAB-BCY*RRBC)/SR3**2+
     P      (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     Q      2*QA*(C2+C8*QB)+C10*QC)*(ABY*RRAB-BCY*RRBC)
     R      /SR3+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     S      (-2*BCY*QQ2*RRBC/SR2+
     T      2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)+
     U      (C6+C10*QA)*(6*BCY*QQ2*QQ3*RRBC/SR2+
     V      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6
     W      )))+(C5+2*C8*QA)*
     X      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)
     Y      +C10*(6*BCY*QQ2*QQ3*RRBC/SR2+
     Z      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6))/
     A      SR3)+4*C9*(ACX*QQ2*RRAC/SR2+
     B      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     C      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)/RE)/ 
     D      RE+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     E      ((-2*ACX*BCY*RRAC*RRBC/SR2**2+
     F      2*(-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     G      RE**2+4*ABX*ABY*QQ3*RRAB**3/(RE*SR6))+
     H      (C6+C10*QA)*(6*(QQ3*
     I      (ACX*BCY*RRAC*RRBC/SR2**2+
     J      (-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     K      RE**2+(ABX*ABY*QQ3**2*RRAB**3+
     L      BCY*QQ2*(-2*ABX*RRAB+ACX*RRAC)*RRBC/(RE*SR2))/(RE*SR6))
     M      -6*(ABX*ABY*QQ2**2*RRAB**3+
     N      ACX*QQ2*RRAC*(2*ABY*RRAB+BCY*RRBC)/(RE*SR2))/(RE*SR6)))/ 
     O      DEXP(AN3*QA)
              ENDIF
270           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
290         CONTINUE
            DO 295 J5=J2+1,3
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              ABX=VEC(J3,J1,J2)
              ABY=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 275 J4=1,N
                NMX=M(J4,J3,1)+M(J3,J1,1)
                NMY=M(J4,J3,2)+M(J3,J1,2)
                NMZ=M(J4,J3,3)+M(J3,J1,3)
                IF ((J4.NE.J3).AND.(J4.NE.J1).AND.
     1             R2(J4,J3).LT.CUTOFF.AND.
     2             R2(J4,J1).LT.CUTOFF.AND.
     3             R2(J3,J1).LT.CUTOFF.AND.
     3             NMZ.EQ.M(J4,J1,3).AND.
     4             NMX.EQ.M(J4,J1,1).AND.NMY.EQ.M(J4,J1,2)) THEN
                BCX=VEC(J4,J3,J2)
                ACX=VEC(J4,J1,J2)
                BCY=VEC(J4,J3,J5)
                ACY=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                RHOAB=(RAB-RE)/RE
                RHOBC=(RBC-RE)/RE
                RHOAC=(RAC-RE)/RE
                QA=(RHOAB+RHOAC+RHOBC)/SR3
                QQ2=(RHOBC-RHOAC)/SR2
                QQ3=(2*RHOAB-RHOAC-RHOBC)/SR6
                QB=(QQ2**2)+(QQ3**2)
                QC=QQ3**3-3*QQ3*QQ2**2
                TEMP=TEMP+
     1      D*(((ABX*ABY*(C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*
     2      (C2+C8*QB)+
     3      C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     4      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC)))
     5      *RRAB**3+(ABY*RRAB-BCY*RRBC)*
     6      (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     7      2*QA*(C2+C8*QB)+C10*QC)*
     8      (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     9      (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     A      (ACX*QQ2*RRAC/SR2+
     B      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     C      (C6+C10*QA)*
     D      (-6*ACX*QQ2*QQ3*RRAC/SR2+
     E      (-3*QQ2**2+3*QQ3**2)*
     F      (-2*ABX*RRAB+ACX*RRAC)/SR6))/RE))+
     G      (2*(C5+2*C8*QA)*
     H      (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     I      +C10*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     J      (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6)
     K      )/RE))/SR3+
     L      (-(ABX*RRAB/RE)-ACX*RRAC/RE)*
     M      ((6*C4*QA+12*C7*QA**2+2*(C2+C8*QB)+
     O      AN3**2*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     P      QA**2*(C2+C8*QB)+C6*QC+QA*(C1+C5*QB+C10*QC))
     Q      )*(ABY*RRAB-BCY*RRBC)/SR3**2+
     R      (-(AN3*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+
     S      2*QA*(C2+C8*QB)+C10*QC)*(ABY*RRAB-BCY*RRBC)
     T      /SR3+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     U      (-2*BCY*QQ2*RRBC/SR2+
     V      2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)+
     W      (C6+C10*QA)*(6*BCY*QQ2*QQ3*RRBC/SR2+
     X      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6
     Y      )))+(C5+2*C8*QA)*
     Z      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)
     A      +C10*(6*BCY*QQ2*QQ3*RRBC/SR2+
     B      (-3*QQ2**2+3*QQ3**2)*(2*ABY*RRAB+BCY*RRBC)/SR6))/
     C      SR3)+4*C9*(ACX*QQ2*RRAC/SR2+
     D      QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)*
     E      (-2*BCY*QQ2*RRBC/SR2+2*QQ3*(2*ABY*RRAB+BCY*RRBC)/SR6)/RE)/
     F      RE+(C3+C5*QA+C8*QA**2+2*C9*QB)*
     G      ((-2*ACX*BCY*RRAC*RRBC/SR2**2+
     H      2*(-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     I      RE**2+4*ABX*ABY*QQ3*RRAB**3/(RE*SR6))+
     J      (C6+C10*QA)*(6*(QQ3*
     K      (ACX*BCY*RRAC*RRBC/SR2**2+
     L      (-2*ABX*RRAB+ACX*RRAC)*(2*ABY*RRAB+BCY*RRBC)/SR6**2)/
     M      RE**2+(ABX*ABY*QQ3**2*RRAB**3+
     N      BCY*QQ2*(-2*ABX*RRAB+ACX*RRAC)*RRBC/(RE*SR2))/(RE*SR6))
     O      -6*(ABX*ABY*QQ2**2*RRAB**3+
     P      ACX*QQ2*RRAC*(2*ABY*RRAB+BCY*RRBC)/(RE*SR2))/(RE*SR6)))/
     Q      DEXP(AN3*QA)
              ENDIF
275           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=HESS(3*(J3-1)+J5,3*(J1-1)+J2)+TEMP
295         CONTINUE
280       CONTINUE
300     CONTINUE
310   CONTINUE
C
C  SYMMETRISE
C
      DO 1000 J1=1,3*N
         DO 1010 J2=J1+1,3*N
C           PRINT*,'J1,J2,A=',J1,J2,HESS(J2,J1)
            HESS(J1,J2)=HESS(J2,J1)
1010     CONTINUE
1000  CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  HERE WE CALCULATE THE POTENTIAL USING THE TWO AND THREE-BODY
C  TERMS OF THE MURRELL POTENTIAL
C                                        
C*************************************************************************
C
      SUBROUTINE JME(N, X, P2, P3, POTEL)
      IMPLICIT NONE 
      INTEGER N, J1, J2, I, J, IJ
      DOUBLE PRECISION SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN
      PARAMETER(SR2=1.414213562D0)
      PARAMETER(SR3=1.732050808D0)
      PARAMETER(SR6=2.449489743D0)
      PARAMETER(C0=3.598D0)
      PARAMETER(C1=-11.609D0)
      PARAMETER(C2=13.486D0)
      PARAMETER(C3=-18.174D0)
      PARAMETER(C4=-5.570D0)
      PARAMETER(C5=79.210D0)
      PARAMETER(C6=-6.458D0)
      PARAMETER(C7=23.383D0)
      PARAMETER(C8=-111.890D0)
      PARAMETER(C9=9.705D0)
      PARAMETER(C10=38.297D0)
      PARAMETER(RE=2.389D0)
      PARAMETER(D=2.918D0)
      PARAMETER(AN=6.50D0)
      DOUBLE PRECISION X(3*N), POTEL, P2, P3, RAB, RHOAB,
     1                 DIST(N,N), SDIST(N,N), RBC, 
     2                 RAC, QQ2, RHOBC, RHOAC, QA, QB, QC, QQ3
C
      P2=0.0D0
      P3=0.0D0
C
C  IT IS VITAL TO TAKE AS MUCH AS POSSIBLE OUT OF THE LOOPS,
C  ESPECIALLY THE INNER LOOP.
C
      DO 20 J1=1,N
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            SDIST(J1,J2)=( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     1                   ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     2                   ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2 
            DIST(J1,J2)=DSQRT(SDIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  CALCULATE THE ENERGY
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J) THEN
               RAB=DIST(I,J)
               RHOAB=(RAB-RE)/RE
               P2=P2-D*(1+AN*RHOAB)*DEXP(-AN*RHOAB)
               DO 24 IJ=1,N
                  IF ((IJ.NE.I).AND.(IJ.NE.J)) THEN
                     RBC=DIST(IJ,J)
                     RAC=DIST(IJ,I)
                     RHOBC=(RBC-RE)/RE
                     RHOAC=(RAC-RE)/RE
                     QA=(RHOAB+RHOBC+RHOAC)/SR3
                     QQ2=(RHOBC-RHOAC)/SR2
                     QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                     QB=QQ2**2+QQ3**2
                     QC=QQ3**3-3*QQ3*QQ2**2
                     P3=P3+
     1                 (C0+C1*QA+C2*(QA**2)+C3*(QB)+C4*(QA**3)+
     2                  C5*QA*QB+C6*QC+C7*(QA**4) +
     3                  C8*(QA**2)*QB+C9*(QB**2)+C10*QA*QC)
     4                *D*DEXP(-AN*QA)
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3/6.0D0
      P2=P2/2.0D0
      POTEL=P2+P3
      RETURN
      END
C
C*************************************************************************
C
C  TWO- AND THREE-BODY TERMS IN THE ENERGY OF THE JM POTENTIAL WITH
C  PERIODIC BOUNDARY CONDITIONS. THIS MUST BE CALLED BEFORE JM2P
C  OR JM3P AS IT DOES SOME SETTING UP FOR THEM.
C                                        
C*************************************************************************
C
      SUBROUTINE JMEC(N,X,P2,P3,VNEW,POTEL,CUTOFF,GTEST,STEST)
      USE PORFUNCS
      IMPLICIT NONE
      INTEGER N, J1, J2, I, J, IJ
      LOGICAL YESNO, GTEST,STEST
      INTEGER M(N,N,3)
      DOUBLE PRECISION DIST(N,N), SDIST(N,N), VNEW(3*N), 
     1              VEC(N,N,3), RH(N,N), NDUM(N,N),
     2              SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF, X(3*N), POTEL, P2, P3, RAB, RHOAB, RBC, 
     1                 RAC, QQ2, RHOBC, RHOAC, QA, QB, QC, QQ3
C
      SR2=1.414213562D0
      SR3=1.732050808D0
      SR6=2.449489743D0
      P2=0.0D0
      P3=0.0D0
      CUTOFF=1.0D10
      INQUIRE(FILE='JMPARAMS',EXIST=YESNO)
      IF (.NOT.YESNO) THEN
         PRINT*,'DATA FILE JMPARAMS NOT FOUND - QUIT'
         STOP
      ELSE
         OPEN(UNIT=33,FILE='JMPARAMS',STATUS='OLD')
         READ(33,*) C0
         READ(33,*) C1
         READ(33,*) C2
         READ(33,*) C3
         READ(33,*) C4
         READ(33,*) C5
         READ(33,*) C6
         READ(33,*) C7
         READ(33,*) C8
         READ(33,*) C9
         READ(33,*) C10
         READ(33,*) RE
         READ(33,*) D
         READ(33,*) AN2
         READ(33,*) AN3
         READ(33,*,END=666) CUTOFF
666      CLOSE(33)
      ENDIF
C
C  CALCULATION OF CONNECTING VECTORS
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  CALCULATION OF DISTANCES:
C
      DO 20 J1=1,N
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            DIST(J1,J2)= VEC(J2,J1,1)**2 +
     1                   VEC(J2,J1,2)**2 +
     2                   VEC(J2,J1,3)**2 
            DIST(J1,J2)=DSQRT(DIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J1,J2)=1/DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  CALCULATE THE ENERGY
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J.AND.DIST(I,J).LT.CUTOFF) THEN 
               RAB=DIST(I,J)
               RHOAB=(RAB-RE)/RE
               P2=P2-D*(1+AN2*RHOAB)*DEXP(-AN2*RHOAB)
               DO 24 IJ=1,N
                  IF ((IJ.NE.I).AND.(IJ.NE.J).AND.
     1                DIST(IJ,I).LT.CUTOFF.AND.
     2                DIST(IJ,J).LT.CUTOFF) THEN
                     RBC=DIST(IJ,J)
                     RAC=DIST(IJ,I)
                     RHOBC=(RBC-RE)/RE
                     RHOAC=(RAC-RE)/RE
                     QA=(RHOAB+RHOBC+RHOAC)/SR3
                     QQ2=(RHOBC-RHOAC)/SR2
                     QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                     QB=QQ2**2+QQ3**2
                     QC=QQ3**3-3*QQ3*QQ2**2
                     P3=P3+
     1                 (C0+C1*QA+C2*(QA**2)+C3*(QB)+C4*(QA**3)+
     2                  C5*QA*QB+C6*QC+C7*(QA**4) +
     3                  C8*(QA**2)*QB+C9*(QB**2)+C10*QA*QC)
     4                *D*DEXP(-AN3*QA)
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3/6.0D0
      P2=P2/2.0D0
      POTEL=P2+P3
      IF (GTEST.OR.STEST) THEN
         CALL JM2C(N, X, VNEW, CUTOFF, DIST, SDIST, VEC, RH, M, NDUM,STEST)
         CALL JM3C(N, X, VNEW, CUTOFF, DIST, SDIST, VEC, RH, M, NDUM,STEST)
      ENDIF
      RETURN
      END
C
C*************************************************************************
C
C  TWO- AND THREE-BODY TERMS IN THE ENERGY OF THE JM POTENTIAL WITH
C  PERIODIC BOUNDARY CONDITIONS. THIS MUST BE CALLED BEFORE JM2P
C  OR JM3P AS IT DOES SOME SETTING UP FOR THEM.
C                                        
C*************************************************************************
C
      SUBROUTINE JMECC(N,X,P2,P3,VNEW,POTEL)
      IMPLICIT NONE
      INTEGER N, J1, J2, I, J, IJ
      INTEGER M(N,N,3)
      DOUBLE PRECISION DIST(N,N), SDIST(N,N), VNEW(3*N), 
     1              VEC(N,N,3), RH(N,N), NDUM(N,N)
      DOUBLE PRECISION CUTOFF, X(3*N), POTEL, P2, P3, RAB, RHOAB, RBC, 
     1                 RAC, QQ2, RHOBC, RHOAC, QC, QQ3, V0, 
     2                 LAMBDA, SR2, SR6, ABIG, B, C, ASMALL
      PARAMETER (SR2=1.414213562D0, SR6=2.449489743D0,  
     1           ABIG=13.29D0, B=0.574D0, C=1.0D0, ASMALL=2.8D0,
     2           V0=60.0D0, LAMBDA=-4.0D0, CUTOFF=ASMALL)
C
      P2=0.0D0
      P3=0.0D0
C
C  CALCULATION OF CONNECTING VECTORS
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  CALCULATION OF DISTANCES:
C
      DO 20 J1=1,N
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            DIST(J1,J2)= VEC(J2,J1,1)**2 +
     1                   VEC(J2,J1,2)**2 +
     2                   VEC(J2,J1,3)**2 
            DIST(J1,J2)=DSQRT(DIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J1,J2)=1/DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
            RH(J1,J2)=(DIST(J1,J2)-ASMALL)
            RH(J2,J1)=RH(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  CALCULATE THE ENERGY
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J.AND.DIST(I,J).LT.CUTOFF) THEN 
               RAB=DIST(I,J)
               RHOAB=RH(I,J)
               P2=P2+ABIG*((B/RAB**4)-1.0D0)*DEXP(C/RHOAB)
               DO 24 IJ=1,N
                  IF ((IJ.NE.I).AND.(IJ.NE.J).AND.
     1                DIST(IJ,I).LT.CUTOFF.AND.
     2                DIST(IJ,J).LT.CUTOFF) THEN
                     RBC=DIST(IJ,J)
                     RAC=DIST(IJ,I)
                     RHOBC=RH(IJ,J)
                     RHOAC=RH(IJ,I)
                     QQ2=(RBC-RAC)/SR2
                     QQ3=(2*RAB-RBC-RAC)/SR6
                     QC=QQ3**3-3.0D0*QQ3*QQ2**2
                     P3=P3+(1.0D0+LAMBDA*QC)*V0
     1                 *DEXP(C*(1.0D0/RHOAB+1.0D0/RHOAC+1.0D0/RHOBC))
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3/6.0D0
      P2=P2/2.0D0
      POTEL=P2+P3

      CALL JM2CC(N, X, VNEW,DIST,SDIST,VEC,RH,M,NDUM)
      CALL JM3CC(N, X, VNEW,DIST,SDIST,VEC,RH,M,NDUM)
      RETURN
      END
C
C*************************************************************************
C
C  TWO- AND THREE-BODY TERMS IN THE ENERGY OF THE JM POTENTIAL WITH
C  PERIODIC BOUNDARY CONDITIONS. THIS MUST BE CALLED BEFORE JM2P
C  OR JM3P AS IT DOES SOME SETTING UP FOR THEM.
C                                        
C*************************************************************************
C

      SUBROUTINE JMEP(N,X,P2,P3,VNEW,POTEL,BOXLX,BOXLY,BOXLZ,CUTOFF)
      USE PORFUNCS
      IMPLICIT NONE
      INTEGER N, J1, J2, I, J, IJ, NMX, NMY, NMZ
      LOGICAL YESNO
      INTEGER M(N,N,3)
      DOUBLE PRECISION DIST(N,N), SDIST(N,N),
     1              VEC(N,N,3), R2(N,N), RR2(N,N),
     2              SR2, SR3, SR6, C0, C1, C2, C3, C4, VNEW(3*N),
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF, X(3*N), POTEL, P2, P3, RAB(N,N), RHOAB, RBC, 
     1                 RAC, QQ2, RHOBC, RHOAC, QA, QB, QC, QQ3, BOXLX, BOXLY, BOXLZ
C
      P2=0.0D0
      P3=0.0D0
      SR2=1.414213562D0
      SR3=1.732050808D0
      SR6=2.449489743D0
      PRINT*,'WARNING - R2, RAB AND RR2 WERE NOT SET - NOT TESTED!'
      INQUIRE(FILE='JMPARAMS',EXIST=YESNO)
      IF (.NOT.YESNO) THEN
         PRINT*,'DATA FILE JMPARAMS NOT FOUND - QUIT'
         STOP
      ELSE
         OPEN(UNIT=33,FILE='JMPARAMS',STATUS='OLD')
         READ(33,*) C0
         READ(33,*) C1
         READ(33,*) C2
         READ(33,*) C3
         READ(33,*) C4
         READ(33,*) C5
         READ(33,*) C6
         READ(33,*) C7
         READ(33,*) C8
         READ(33,*) C9
         READ(33,*) C10
         READ(33,*) RE
         READ(33,*) D
         READ(33,*) AN2
         READ(33,*) AN3
         CLOSE(33)
      ENDIF
C
C  ATOMS LEAVING THE BOX IN X OR Y DIRECTION ENTER AGAIN ON THE
C  OTHER SIDE.
C
      DO 21 J1=1,N
         X(3*(J1-1)+1)=X(3*(J1-1)+1) - BOXLX *
     1                  DNINT(X(3*(J1-1)+1)/BOXLX)
         X(3*(J1-1)+2)=X(3*(J1-1)+2) -  BOXLY *
     1                  DNINT(X(3*(J1-1)+2)/BOXLY)
         X(3*(J1-1)+3)=X(3*(J1-1)+3) -  BOXLZ *
     1                  DNINT(X(3*(J1-1)+3)/BOXLZ)
21    CONTINUE
C
C  CALCULATION OF CONNECTING VECTORS; TO IMPLEMENT THE PERIODIC
C  BOUNDARY CONDITIONS, THE SHORTEST VECTOR BETWEEN TWO ATOMS IS
C  USED:
C
      DO 25 J1=1,N
         R2(J1,J1)=0.0D0
         RR2(J1,J1)=0.0D0
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         M(J1,J1,1)=0
         M(J1,J1,2)=0
         M(J1,J1,3)=0
         DO 15 J2=J1+1,N
            R2(J2,J1)=
     1                ( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     2                ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     3                ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2 
            R2(J2,J1)=DSQRT(R2(J2,J1))
            RR2(J2,J1)=1.0/R2(J2,J1)
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            M(J2,J1,1)=NINT(VEC(J2,J1,1)/BOXLX)
            M(J2,J1,2)=NINT(VEC(J2,J1,2)/BOXLY)
            M(J2,J1,3)=NINT(VEC(J2,J1,3)/BOXLZ)
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * FLOAT(M(J2,J1,1))
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * FLOAT(M(J2,J1,2))
            VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * FLOAT(M(J2,J1,3))
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
            M(J1,J2,1)=-M(J2,J1,1)
            M(J1,J2,2)=-M(J2,J1,2)
            M(J1,J2,3)=-M(J2,J1,3)
            RR2(J1,J2)=RR2(J2,J1) 
            R2(J1,J2)=R2(J2,J1) 
15       CONTINUE
25    CONTINUE
C
C  CALCULATION OF DISTANCES:
C
      DO 20 J1=1,N
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            DIST(J1,J2)= VEC(J2,J1,1)**2 +
     1                   VEC(J2,J1,2)**2 +
     2                   VEC(J2,J1,3)**2 
            DIST(J1,J2)=DSQRT(DIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J1,J2)=1/DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  CALCULATE THE ENERGY
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J.AND.DIST(I,J).LT.CUTOFF) THEN 
               RAB(I,J)=DIST(I,J)
               RHOAB=(RAB(I,J)-RE)/RE
               P2=P2-D*(1+AN2*RHOAB)*DEXP(-AN2*RHOAB)
               DO 24 IJ=1,N
                  NMX=M(I,J,1)+M(J,IJ,1)
                  NMY=M(I,J,2)+M(J,IJ,2)
                  NMZ=M(I,J,3)+M(J,IJ,3)
                  IF ((IJ.NE.I).AND.(IJ.NE.J).AND.
     1                DIST(IJ,I).LT.CUTOFF.AND.
     2                DIST(IJ,J).LT.CUTOFF.AND.
     3                NMZ.EQ.M(I,IJ,3).AND.
     3                NMX.EQ.M(I,IJ,1).AND.NMY.EQ.M(I,IJ,2)) THEN
                     RBC=DIST(IJ,J)
                     RAC=DIST(IJ,I)
                     RHOBC=(RBC-RE)/RE
                     RHOAC=(RAC-RE)/RE
                     QA=(RHOAB+RHOBC+RHOAC)/SR3
                     QQ2=(RHOBC-RHOAC)/SR2
                     QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                     QB=QQ2**2+QQ3**2
                     QC=QQ3**3-3*QQ3*QQ2**2
                     P3=P3+
     1                 (C0+C1*QA+C2*(QA**2)+C3*(QB)+C4*(QA**3)+
     2                  C5*QA*QB+C6*QC+C7*(QA**4) +
     3                  C8*(QA**2)*QB+C9*(QB**2)+C10*QA*QC)
     4                *D*DEXP(-AN3*QA)
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3/6.0D0
      P2=P2/2.0D0
      POTEL=P2+P3
      CALL JM2P(N,X,VNEW,BOXLX,BOXLY,BOXLZ,CUTOFF,RAB,VEC)
      CALL JM3P(N,X,VNEW,BOXLX,BOXLY,BOXLZ,CUTOFF, R2, RR2, VEC, M)
      RETURN
      END
