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
      SUBROUTINE MSDIFF(N,X,V,P2,RHO,BOXLX,BOXLY,CUTOFF)
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5, J6, I, J, MX, MY
      DOUBLE PRECISION X(3*N), RHO, 
     1                 V(3*N), DIST, R(N,N), RR(N,N), P2, MEXP(N,N),
     2                 BOXLX, BOXLY, VEC(N,N,3), 
     3                 CUTOFF, RCUT, RRCUT
C    1              VEC(N,N,3)

      IF (BOXLX.LT.BOXLY) THEN
        RCUT=BOXLX*CUTOFF
      ELSE 
        RCUT=BOXLY*CUTOFF
      ENDIF
      PRINT*,'CUTOFF USED = ',RCUT

C
C  DEAL WITH ATOMS LEAVING THE BOX:
C
      DO 41 J1=1,N
         X(3*(J1-1)+1)=X(3*(J1-1)+1) - BOXLX*DNINT(X(3*(J1-1)+1)/BOXLX)
         X(3*(J1-1)+2)=X(3*(J1-1)+2) - BOXLY*DNINT(X(3*(J1-1)+2)/BOXLY)
41    CONTINUE
C
C  CALCULATION OF CONNECTING VECTORS; TO IMPLEMENT THE PERIODIC
C  BOUNDARY CONDITIONS, THE SHORTEST VECTOR BETWEEN TWO ATOMS IS
C  USED:
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            MX=NINT(VEC(J2,J1,1)/BOXLX)
            MY=NINT(VEC(J2,J1,2)/BOXLY)
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * MX
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * MY
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  STORE DISTANCE MATRICES.
C
      P2=0.0D0
      DO 20 J1=1,N
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         DO 10 J2=1,J1-1
            DIST=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 +
     1           VEC(J1,J2,3)**2
            R(J2,J1)=RHO*DSQRT(DIST)
            R(J1,J2)=R(J2,J1)
            RR(J2,J1)=1.0D0/R(J2,J1)
            RR(J1,J2)=RR(J2,J1)
            MEXP(J1,J2)=DEXP(RHO-R(J1,J2))
            MEXP(J2,J1)=MEXP(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  HERE WE CALCULATE THE ENERGY
C                                        
      RRCUT=RHO*RCUT
      DO 22 I=1,N
         DO 23 J=I+1,N
            IF ((I.NE.J).AND.(R(I,J).LT.RRCUT)) THEN
               P2=P2+MEXP(I,J)*(MEXP(I,J)-2.0D0)
            ENDIF
23       CONTINUE
22    CONTINUE
C
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
              IF (R(J4,J1).LT.RRCUT) THEN 
                V(J3)=V(J3)-2.0*RHO*RHO*VEC(J1,J4,J2)*
     1               MEXP(J4,J1)* 
     2              (MEXP(J4,J1)-1.0D0)*RR(J4,J1)
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
            HESS(J3,J3)=0.0D0
            DO 60 J4=1,N
              IF (R(J4,J1).LT.RRCUT) THEN 
               HESS(J3,J3)=HESS(J3,J3)+2.0*RHO*RHO*MEXP(J4,J1)
     1         *(  (MEXP(J4,J1)-1.0D0) 
     2         * ( (RHO * VEC(J1,J4,J2) *RR(J4,J1) )**2-1.0D0)  
     4         + ( RHO*VEC(J1,J4,J2) )**2*RR(J4,J1)*
     5             (2.0*MEXP(J4,J1)-1.0D0) ) * RR(J4,J1)
              ENDIF 
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
                IF (R(J5,J1).LT.RRCUT) THEN 
                 HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           +2.0*RHO**4*MEXP(J5,J1)
     2           *VEC(J1,J5,J2)*VEC(J1,J5,J4)
     3           *( (MEXP(J5,J1)-1.0D0)*RR(J5,J1)+
     4              2.0*MEXP(J5,J1)-1.0D0) * RR(J5,J1)**2
                ENDIF 
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
              IF (R(J4,J1).LT.RRCUT) THEN 
               HESS(J3,3*(J4-1)+J2)=2.0*RHO**2*MEXP(J4,J1)
     1      *(  (-RHO**2*VEC(J1,J4,J2)**2*RR(J4,J1)**2+1.0) 
     2           *   (MEXP(J4,J1)-1.0D0)
     3           -   VEC(J1,J4,J2)**2 * RR(J4,J1) * RHO**2
     4           *   (2.0*MEXP(J4,J1)-1.0D0)   )  * RR(J4,J1)
              ENDIF 
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
                IF (R(J4,J1).LT.RRCUT) THEN
                  HESS(J3,J6)=
     1           -2.0*RHO**4*MEXP(J4,J1)
     2           *VEC(J1,J4,J2)*VEC(J1,J4,J5)
     3           *( (MEXP(J4,J1)-1.0D0)*RR(J4,J1)+
     4              2.0*MEXP(J4,J1)-1.0D0) * RR(J4,J1)**2
                ENDIF
                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                 J6=3*(J4-1)+J5
                IF (R(J4,J1).LT.RRCUT) THEN
                  HESS(J3,J6)=
     1           -2.0*RHO**4*MEXP(J4,J1)
     2           *VEC(J1,J4,J2)*VEC(J1,J4,J5)
     3           *( (MEXP(J4,J1)-1.0D0)*RR(J4,J1)+
     4              2.0*MEXP(J4,J1)-1.0D0) * RR(J4,J1)**2
                ENDIF
                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
