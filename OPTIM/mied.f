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
C  CALCULATE THE ANALYTIC FIRST AND SECOND DERIVATIVES FOR THE BULK
C  PERIODIC MIE POTENTIAL.
C
C*************************************************************************
C
      SUBROUTINE MIED(N,X,V,RHO,RE2,ALPHA,XN,AWELL,RE1) 
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), RNF, RMF,RN,RM
      DOUBLE PRECISION VEC(N,N,3),RSQ1(N,N),
     1              R8(N,N),R10(N,N),R14(N,N),
     2              R16(N,N)
C 
      DOUBLE PRECISION RHO, RE2, ALPHA, XN, AWELL, RE1, BOXLX,BOXLY,BOXLZ,CUTOFF
      INTEGER NEXP,MEXP
      NEXP=DINT(RHO)
      MEXP=DINT(RE2)
      BOXLX=ALPHA
      BOXLY=XN
      BOXLZ=AWELL
      CUTOFF=RE1
      
      RN=FLOAT(NEXP)
      RM=FLOAT(MEXP)
      RNF=RM/(RN-RM)
      RMF=RN/(RN-RM)
      PRINT*, 'N, M, N/(N-M), M/(N-M)=', IDNINT(RN), IDNINT(RM), RNF, RMF
C     PRINT*,'INCOORDS'
C     PRINT*, (X(I),X(I+1),X(I+2),I=1,3*N,3)
C
C  STORE DISTANCE MATRICES.
C
      PRINT*,'WARNING - VEC WAS NOT SET! NOT TESTED!'

      DO 20 J1=1,N
         R8(J1,J1)=0.0D0
         R10(J1,J1)=0.0D0
         R14(J1,J1)=0.0D0
         R16(J1,J1)=0.0D0
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 10 J2=1,J1-1
            IF (DSQRT(RSQ1(J2,J1)).LT.CUTOFF) THEN
               R8(J2,J1)=RSQ1(J2,J1)**(-(RM/2.0D0+1.0D0))
               R10(J2,J1)=RSQ1(J2,J1)**(-(RM/2.0D0+2.0D0))
               R14(J2,J1)=RSQ1(J2,J1)**(-(RN/2.0D0+1.0D0))
               R16(J2,J1)=RSQ1(J2,J1)**(-(RN/2.0D0+2.0D0))
               VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
               VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
               VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            ELSE
               R8(J2,J1)=0.0D0
               R10(J2,J1)=0.0D0
               R14(J2,J1)=0.0D0
               R16(J2,J1)=0.0D0
               VEC(J2,J1,1)=0.0D0
               VEC(J2,J1,2)=0.0D0
               VEC(J2,J1,3)=0.0D0
            ENDIF
            R8(J1,J2)=R8(J2,J1)
            R10(J1,J2)=R10(J2,J1)
            R14(J1,J2)=R14(J2,J1)
            R16(J1,J2)=R16(J2,J1)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
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
               V(J3)=V(J3)-(RN*RNF*R14(J1,J4)-RM*RMF*R8(J1,J4))
     1              *VEC(J1,J4,J2)
30          CONTINUE
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
               HESS(J3,J3)=HESS(J3,J3)+(RN*(RN+2.0D0)*RNF*R16(J4,J1)
     1                  -RM*(RM+2.0D0)*RMF*R10(J4,J1))
     2                  *VEC(J1,J4,J2)**2
     3                  -RN*RNF*R14(J4,J1)+RM*RMF*R8(J4,J1)
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
            DO 100 J4=1,J2-1
               HESS(J3,3*(J1-1)+J4)=0.0D0
               DO 90 J5=1,N
                  HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           +(RN*(RN+2.0D0)*RNF*R16(J1,J5)
     2           -RM*(RM+2.0D0)*RMF*R10(J1,J5))
     3           *VEC(J1,J5,J4)*VEC(J1,J5,J2)
90             CONTINUE
               HESS(3*(J1-1)+J4,J3)=HESS(J3,3*(J1-1)+J4)
100         CONTINUE
110      CONTINUE
120   CONTINUE
C
C  CASE III, DIFFERENT ATOMS, SAME CARTESIAN COORDINATE.
C
      DO 150 J1=1,N
         DO 140 J2=1,3
            J3=3*(J1-1)+J2
            DO 130 J4=1,J1-1
               HESS(J3,3*(J4-1)+J2)=-(RN*(RN+2.0D0)*RNF*R16(J1,J4)
     1                            -RM*(RM+2.0D0)*RMF*R10(J1,J4)) 
     2                            *VEC(J1,J4,J2)**2
     3                            +RN*RNF*R14(J1,J4)
     4                            -RM*RMF*R8(J1,J4)
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
            DO 160 J4=1,J1-1
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=-(RN*(RN+2.0D0)*RNF*R16(J1,J4)
     1                     -RM*(RM+2.0D0)*RMF*R10(J1,J4))
     2                     *VEC(J1,J4,J2)*VEC(J1,J4,J5)
                  HESS(J6,J3)=HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=-(RN*(RN+2)*RNF*R16(J1,J4)
     1                     -RM*(RM+2)*RMF*R10(J1,J4))
     2                     *VEC(J1,J4,J2)*VEC(J1,J4,J5)
                  HESS(J6,J3)=HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
