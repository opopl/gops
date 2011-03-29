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
C  SUBROUTINE MDIFF CALCULATES THE CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY FOR THE MORSE POTENTIAL.
C
C*************************************************************************
C
      SUBROUTINE MDIFF(N, X, V, RHO) 
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5, J6 
      DOUBLE PRECISION X(3*N), RHO,
     1                 V(3*N), DIST, R(N,N), RR(N,N)
C 
C  STORE DISTANCE MATRICES.  
C 
      DO 20 J1=1,N 
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            DIST=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1          +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2          +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R(J2,J1)=RHO*DSQRT(DIST)
            RR(J2,J1)=1.0D0/R(J2,J1)
            R(J1,J2)=R(J2,J1)
            RR(J1,J2)=RR(J2,J1)
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
               V(J3)=V(J3)-2.0*RHO*RHO*(X(J3)-X(3*(J4-1)+J2))*
     1               DEXP(RHO-R(J4,J1))* 
     2              (DEXP(RHO-R(J4,J1))-1.0D0)*RR(J4,J1)
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
               HESS(J3,J3)=HESS(J3,J3)+2.0*RHO*RHO*DEXP(RHO-R(J4,J1))
     1         *(  (DEXP(RHO-R(J4,J1))-1.0D0) 
     2         * ( (RHO * (X(J3)-X(3*(J4-1)+J2) )*RR(J4,J1) )**2-1.0D0)  
     4         + ( RHO*(X(J3)-X(3*(J4-1)+J2)) )**2*RR(J4,J1)*
     5             (2.0*DEXP(RHO-R(J4,J1))-1.0D0) ) * RR(J4,J1)
60          CONTINUE
C           PRINT*,'HESS(',J3,',',J3,')=',HESS(J3,J3)
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
               HESS(J3,3*(J1-1)+J4)=0.0D0
               DO 90 J5=1,N
                 HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           +2.0*RHO**4*DEXP(RHO-R(J5,J1))
     2           *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
     3           *( (DEXP(RHO-R(J5,J1))-1.0D0)*RR(J5,J1)+
     4              2.0*DEXP(RHO-R(J5,J1))-1.0D0) * RR(J5,J1)**2
90             CONTINUE
               HESS(3*(J1-1)+J4,J3)=HESS(J3,3*(J1-1)+J4)
C              PRINT*,'HESS(',J3,',',3*(J1-1)+J4,')=',HESS(J3,3*(J1-1)+J4)
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
               HESS(J3,3*(J4-1)+J2)=2.0*RHO**2*DEXP(RHO-R(J4,J1))
     1      *(  (-RHO**2*(X(J3)-X(3*(J4-1)+J2))**2*RR(J4,J1)**2+1.0) 
     2           *   (DEXP(RHO-R(J4,J1))-1.0D0)
     3           -   (X(J3)-X(3*(J4-1)+J2))**2 * RR(J4,J1) * RHO**2
     4           *   (2.0*DEXP(RHO-R(J4,J1))-1.0D0)   )  * RR(J4,J1)
C              PRINT*,'HESS(',J3,',',3*(J4-1)+J2,')=',HESS(J3,3*(J4-1)+J2)
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
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
                  HESS(J3,J6)=
     1           -2.0*RHO**4*DEXP(RHO-R(J4,J1))
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3           *( (DEXP(RHO-R(J4,J1))-1.0D0)*RR(J4,J1)+
     4              2.0*DEXP(RHO-R(J4,J1))-1.0D0) * RR(J4,J1)**2

                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=
     1           -2.0*RHO**4*DEXP(RHO-R(J4,J1))
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3           *( (DEXP(RHO-R(J4,J1))-1.0D0)*RR(J4,J1)+
     4              2.0*DEXP(RHO-R(J4,J1))-1.0D0) * RR(J4,J1)**2
                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
