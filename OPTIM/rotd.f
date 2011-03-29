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
C**************************************************************************
C
C  SUBROUTINE ROTD CALCULATES ANALYTICALLY THE CARTESIAN GRADIENT AND
C  SECOND DERIVATIVE MATRIX DUE TO THE ROTATIONAL TERM ONLY AND INCREMENTS
C  THE EXISTING MATRICES CONTAINING THE POTENTIAL CONTRIBUTIONS 
C  ACCORDINGLY.
C
C**************************************************************************
C
      SUBROUTINE ROTD(N, X, V, MASS, ZJ, FLAG, ROTE)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION X(3*N), ROTE,
     1                 V(3*N), MASS, ZJ, IZINV, XA, XB, YA, YB, ZA,ZB,  
     2                 A, B, C, D, E, F, TOP, BOTTOM, XM, YM, ZM, XOM
      LOGICAL FLAG
      XOM=1.0D0-1.0D0/FLOAT(N)
C
C  CALCULATE CENTRE OF MASS COMPONENTS XM, YM, ZM
C
      XM=0.0D0
      YM=0.0D0
      ZM=0.0D0
      DO J1=1,N
         XM=XM+X(3*(J1-1)+1)
         YM=YM+X(3*(J1-1)+2)
         ZM=ZM+X(3*(J1-1)+3)
      ENDDO
      XM=XM/N
      YM=YM/N
      ZM=ZM/N
C
C  CALCULATE CONSTANTS A, B, C, D, E, F
C
      A=0.0D0
      B=0.0D0
      C=0.0D0
      D=0.0D0
      E=0.0D0
      F=0.0D0
      DO J1=1,N
         J3=3*(J1-1)
         A=A+(X(J3+2)-YM)**2
         B=B+(X(J3+3)-ZM)**2
         E=E+(X(J3+1)-XM)**2
         C=C+(X(J3+2)-YM)*(X(J3+1)-XM)
         D=D+(X(J3+3)-ZM)*(X(J3+1)-XM)
         F=F+(X(J3+3)-ZM)*(X(J3+2)-YM)
      ENDDO
      TOP=A*B+B*B-C*C+A*E+B*E
      BOTTOM=A*A*B+A*B*B-A*C*C-B*D*D+A*A*E+2.0D0*A*B*E+B*B*E-C*C*E-D*D*E
     1      +A*E*E+B*E*E-2.0D0*C*D*F-A*F*F-B*F*F
      IF (DABS(BOTTOM).LT.1.0D-6) THEN
         PRINT*,'PANIC CALLED, BOTTOM=',BOTTOM
         CALL EPANIC(X,N,ZJ,ROTE,FLAG,V,MASS)
         RETURN
      ELSE
         IZINV=TOP/(MASS*BOTTOM)
         ROTE=ZJ**2*IZINV/2.0D0
      ENDIF

      IF (FLAG) RETURN
C
C  ANALYTIC GRADIENT:
C
      DO J1=1,N
         J3=3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         V(J3+1)=V(J3+1)+ZJ**2*(2*(-(A*XM)-B*XM+A*XA+B*XA+C*YM-C*YA)
     &    /BOTTOM +
     &    2*TOP*(A**2*XM+2*A*B*XM+B**2*XM-C**2*XM-D**2*XM+2*A*E*XM+  
     &    2*B*E*XM-A**2*XA-2*A*B*XA-B**2*XA+C**2*XA+D**2*XA -
     &    2*A*E*XA-2*B*E*XA-A*C*YM-C*E*YM-D*F*YM+A*C*YA+C*E*YA +
     &    D*F*YA-B*D*ZM-D*E*ZM-C*F*ZM+B*D*ZA+D*E*ZA+C*F*ZA)/
     &    BOTTOM**2) / (2.0D0*MASS) 

         V(J3+2)=V(J3+2)+ZJ**2*(
     &    2*(C*XM-C*XA-B*YM-E*YM+B*YA+E*YA)/BOTTOM -
     &    2*TOP*(A*C*XM+C*E*XM+D*F*XM-A*C*XA-C*E*XA-D*F*XA-2*A*B*YM- 
     &    B**2*YM+C**2*YM-2*A*E*YM-2*B*E*YM-E**2*YM+F**2*YM +
     &    2*A*B*YA+B**2*YA-C**2*YA+2*A*E*YA+2*B*E*YA+E**2*YA -
     &    F**2*YA+C*D*ZM+A*F*ZM+B*F*ZM-C*D*ZA-A*F*ZA-B*F*ZA)/
     &    BOTTOM**2 ) / (2.0D0*MASS)

         V(J3+3)=V(J3+3)+ZJ**2*(
     &    2*(A+2*B+E)*(-ZM+ZA)/BOTTOM -
     &    2*TOP*(B*D*XM+D*E*XM+C*F*XM-B*D*XA-D*E*XA-C*F*XA+C*D*YM +
     &    A*F*YM+B*F*YM-C*D*YA-A*F*YA-B*F*YA-A**2*ZM-2*A*B*ZM +
     &    D**2*ZM-2*A*E*ZM-2*B*E*ZM-E**2*ZM+F**2*ZM+A**2*ZA +
     &    2*A*B*ZA-D**2*ZA+2*A*E*ZA+2*B*E*ZA+E**2*ZA-F**2*ZA)/
     &    BOTTOM**2) / (2.0D0*MASS)

      ENDDO
C
C  CASE I: DIAGONAL TERMS.
C
      DO J1 = 1, N
         J3 = 3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         HESS(J3+1, J3+1) = HESS(J3+1, J3+1)+ZJ**2 * (
     &  (2*A+2*B-2*A/N-2*B/N-2*YA**2+4*YA*YM-2*YM**2)/BOTTOM- 
     &  TOP*(2*A**2*XOM+4*A*B*XOM+
     &      2*B**2*XOM-2*C**2*XOM-
     &      2*D**2*XOM+4*A*E*XOM+
     &      4*B*E*XOM+8*A*(XA-XM)**2+8*B*(XA-XM)**2-
     &      8*C*(XA-XM)*(YA-YM)-2*A*(YA-YM)**2-2*E*(YA-YM)**2-
     &      8*D*(XA-XM)*(ZA-ZM)-4*F*(YA-YM)*(ZA-ZM)-
     &      2*B*(ZA-ZM)**2-2*E*(ZA-ZM)**2)/BOTTOM**2-
     &  8*(A*XA+B*XA-A*XM-B*XM-C*YA+C*YM)*
     &    (A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+2*A*E*XA+
     &      2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+D**2*XM-
     &      2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+A*C*YM+
     &      C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+D*E*ZM+
     &      C*F*ZM)/BOTTOM**2+
     &  8*TOP*(A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+2*A*E*XA+
     &       2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+D**2*XM-
     &       2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+A*C*YM+
     &       C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+D*E*ZM+
     &       C*F*ZM)**2/BOTTOM**3)/(2.0D0*MASS)

         HESS(J3+2, J3+2) = HESS(J3+2, J3+2)+ZJ**2 * (
     &   (2*B+2*E-2*B/N-2*E/N-2*XA**2+4*XA*XM-2*XM**2)/BOTTOM-   
     &  TOP*(4*A*B*XOM+2*B**2*XOM-
     &      2*C**2*XOM+4*A*E*XOM+
     &      4*B*E*XOM+2*E**2*XOM-
     &      2*F**2*XOM-2*A*(XA-XM)**2-2*E*(XA-XM)**2-
     &      8*C*(XA-XM)*(YA-YM)+8*B*(YA-YM)**2+8*E*(YA-YM)**2-
     &      4*D*(XA-XM)*(ZA-ZM)-8*F*(YA-YM)*(ZA-ZM)-
     &      2*A*(ZA-ZM)**2-2*B*(ZA-ZM)**2)/BOTTOM**2+
     &  8*(C*XA-C*XM-B*YA-E*YA+B*YM+E*YM)*
     &    (-(A*C*XA)-C*E*XA-D*F*XA+A*C*XM+C*E*XM+D*F*XM+
     &      2*A*B*YA+B**2*YA-C**2*YA+2*A*E*YA+2*B*E*YA+E**2*YA-
     &      F**2*YA-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM-2*B*E*YM-
     &      E**2*YM+F**2*YM-C*D*ZA-A*F*ZA-B*F*ZA+C*D*ZM+A*F*ZM+
     &      B*F*ZM)/BOTTOM**2+
     &  8*TOP*(-(A*C*XA)-C*E*XA-D*F*XA+A*C*XM+C*E*XM+D*F*XM+
     &       2*A*B*YA+B**2*YA-C**2*YA+2*A*E*YA+2*B*E*YA+E**2*YA-
     &       F**2*YA-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM-2*B*E*YM-
     &       E**2*YM+F**2*YM-C*D*ZA-A*F*ZA-B*F*ZA+C*D*ZM+A*F*ZM+
     &       B*F*ZM)**2/BOTTOM**3)/(2.0D0*MASS)

        HESS(J3+3, J3+3) = HESS(J3+3, J3+3)+ZJ**2 * (
     &   -(TOP*(2*A**2*XOM+4*A*B*XOM-2*D**2*XOM+
     &        4*A*E*XOM+4*B*E*XOM+
     &        2*E**2*XOM-2*F**2*XOM-
     &        2*B*(XA-XM)**2-2*E*(XA-XM)**2-4*C*(XA-XM)*(YA-YM)-
     &        2*A*(YA-YM)**2-2*B*(YA-YM)**2-8*D*(XA-XM)*(ZA-ZM)-
     &        8*F*(YA-YM)*(ZA-ZM)+8*A*(ZA-ZM)**2+8*E*(ZA-ZM)**2)/
     &     BOTTOM**2)+8*TOP*(B*D*XA+D*E*XA+C*F*XA-B*D*XM-D*E*XM-
     &       C*F*XM+C*D*YA+A*F*YA+B*F*YA-C*D*YM-A*F*YM-B*F*YM-
     &       A**2*ZA-2*A*B*ZA+D**2*ZA-2*A*E*ZA-2*B*E*ZA-E**2*ZA+
     &       F**2*ZA+A**2*ZM+2*A*B*ZM-D**2*ZM+2*A*E*ZM+2*B*E*ZM+
     &       E**2*ZM-F**2*ZM)**2/BOTTOM**3-
     &  8*(A+2*B+E)*(ZA-ZM)*
     &    (-(B*D*XA)-D*E*XA-C*F*XA+B*D*XM+D*E*XM+C*F*XM-C*D*YA-
     &      A*F*YA-B*F*YA+C*D*YM+A*F*YM+B*F*YM+A**2*ZA+2*A*B*ZA-
     &      D**2*ZA+2*A*E*ZA+2*B*E*ZA+E**2*ZA-F**2*ZA-A**2*ZM-
     &      2*A*B*ZM+D**2*ZM-2*A*E*ZM-2*B*E*ZM-E**2*ZM+F**2*ZM)/
     &   BOTTOM**2+(2*A+4*B+2*E-2*A/N-4*B/N-2*E/N+
     &     8*ZA**2-16*ZA*ZM+8*ZM**2)/BOTTOM)/(2.0D0*MASS)

      ENDDO
C
C  CASE II: SAME CARTESIAN COORDINATE, DIFFERENT ATOM.
C
      DO J1=1,N
         J3=3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         DO J2=J1+1,N
            J4=3*(J2-1)
            XB=X(J4+1)
            YB=X(J4+2)
            ZB=X(J4+3)
            HESS(J3+1,J4+1)=HESS(J3+1,J4+1)+ZJ**2*(
     1      (-2*A/N-2*B/N-2*(YA-YM)*(YB-YM))/BOTTOM-
     1     TOP*(-2*A**2/N-4*A*B/N-2*B**2/N+2*C**2/N+
     1     2*D**2/N-4*A*E/N-4*B*E/N+
     1     8*A*(XA-XM)*(XB-XM)+8*B*(XA-XM)*(XB-XM)-
     1     4*C*(XB-XM)*(YA-YM)-4*C*(XA-XM)*(YB-YM)-
     1     2*A*(YA-YM)*(YB-YM)-2*E*(YA-YM)*(YB-YM)-
     1     4*D*(XB-XM)*(ZA-ZM)-2*F*(YB-YM)*(ZA-ZM)-
     1     4*D*(XA-XM)*(ZB-ZM)-2*F*(YA-YM)*(ZB-ZM)-
     1     2*B*(ZA-ZM)*(ZB-ZM)-2*E*(ZA-ZM)*(ZB-ZM))/BOTTOM**2-
     1     4*(A*XB+B*XB-A*XM-B*XM-C*YB+C*YM)*
     1     (A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+2*A*E*XA+
     1     2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+D**2*XM-
     1     2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+A*C*YM+
     1     C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+D*E*ZM+
     1     C*F*ZM)/BOTTOM**2-
     1     4*(A*XA+B*XA-A*XM-B*XM-C*YA+C*YM)*
     1     (A**2*XB+2*A*B*XB+B**2*XB-C**2*XB-D**2*XB+2*A*E*XB+
     1     2*B*E*XB-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+D**2*XM-
     1     2*A*E*XM-2*B*E*XM-A*C*YB-C*E*YB-D*F*YB+A*C*YM+
     1     C*E*YM+D*F*YM-B*D*ZB-D*E*ZB-C*F*ZB+B*D*ZM+D*E*ZM+
     1     C*F*ZM)/BOTTOM**2+
     1     8*TOP*(A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+2*A*E*XA+
     1     2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+D**2*XM-
     1     2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+A*C*YM+
     1     C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+D*E*ZM+
     1     C*F*ZM)*(A**2*XB+2*A*B*XB+B**2*XB-C**2*XB-D**2*XB+
     1     2*A*E*XB+2*B*E*XB-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+
     1     D**2*XM-2*A*E*XM-2*B*E*XM-A*C*YB-C*E*YB-D*F*YB+
     1     A*C*YM+C*E*YM+D*F*YM-B*D*ZB-D*E*ZB-C*F*ZB+B*D*ZM+
     1     D*E*ZM+C*F*ZM)/BOTTOM**3)/(2.0D0*MASS)
            HESS(J4+1,J3+1)=HESS(J3+1,J4+1)

            HESS(J3+2,J4+2)=HESS(J3+2,J4+2)+ZJ**2*(
     &      (-2*B/N-2*E/N-2*(XA-XM)*(XB-XM))/BOTTOM-
     &  TOP*(-4*A*B/N-2*B**2/N+2*C**2/N-4*A*E/N-
     &      4*B*E/N-2*E**2/N+2*F**2/N-
     &      2*A*(XA-XM)*(XB-XM)-2*E*(XA-XM)*(XB-XM)-
     &      4*C*(XB-XM)*(YA-YM)-4*C*(XA-XM)*(YB-YM)+
     &      8*B*(YA-YM)*(YB-YM)+8*E*(YA-YM)*(YB-YM)-
     &      2*D*(XB-XM)*(ZA-ZM)-4*F*(YB-YM)*(ZA-ZM)-
     &      2*D*(XA-XM)*(ZB-ZM)-4*F*(YA-YM)*(ZB-ZM)-
     &      2*A*(ZA-ZM)*(ZB-ZM)-2*B*(ZA-ZM)*(ZB-ZM))/BOTTOM**2+
     &  8*TOP*(A*C*XA+C*E*XA+D*F*XA-A*C*XM-C*E*XM-D*F*XM-
     &      2*A*B*YA-B**2*YA+C**2*YA-2*A*E*YA-2*B*E*YA-E**2*YA+
     &      F**2*YA+2*A*B*YM+B**2*YM-C**2*YM+2*A*E*YM+2*B*E*YM+
     &      E**2*YM-F**2*YM+C*D*ZA+A*F*ZA+B*F*ZA-C*D*ZM-A*F*ZM-
     &      B*F*ZM)*(A*C*XB+C*E*XB+D*F*XB-A*C*XM-C*E*XM-D*F*XM-
     &      2*A*B*YB-B**2*YB+C**2*YB-2*A*E*YB-2*B*E*YB-E**2*YB+
     &      F**2*YB+2*A*B*YM+B**2*YM-C**2*YM+2*A*E*YM+2*B*E*YM+
     &      E**2*YM-F**2*YM+C*D*ZB+A*F*ZB+B*F*ZB-C*D*ZM-A*F*ZM-
     &      B*F*ZM)/BOTTOM**3+
     &  4*(C*XB-C*XM-B*YB-E*YB+B*YM+E*YM)*
     &    (-(A*C*XA)-C*E*XA-D*F*XA+A*C*XM+C*E*XM+D*F*XM+
     &      2*A*B*YA+B**2*YA-C**2*YA+2*A*E*YA+2*B*E*YA+E**2*YA-
     &      F**2*YA-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM-2*B*E*YM-
     &      E**2*YM+F**2*YM-C*D*ZA-A*F*ZA-B*F*ZA+C*D*ZM+A*F*ZM+
     &      B*F*ZM)/BOTTOM**2+
     &  4*(C*XA-C*XM-B*YA-E*YA+B*YM+E*YM)*
     &    (-(A*C*XB)-C*E*XB-D*F*XB+A*C*XM+C*E*XM+D*F*XM+
     &      2*A*B*YB+B**2*YB-C**2*YB+2*A*E*YB+2*B*E*YB+E**2*YB-
     &      F**2*YB-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM-2*B*E*YM-
     &      E**2*YM+F**2*YM-C*D*ZB-A*F*ZB-B*F*ZB+C*D*ZM+A*F*ZM+
     &      B*F*ZM)/BOTTOM**2)/(2.0D0*MASS)
            HESS(J4+2,J3+2)=HESS(J3+2,J4+2)

            HESS(J3+3,J4+3)=HESS(J3+3,J4+3)+ZJ**2*(
     &   (-2*A/N-4*B/N-2*E/N+8*(ZA-ZM)*(ZB-ZM))/BOTTOM-  
     &  TOP*(-2*A**2/N-4*A*B/N+2*D**2/N-4*A*E/N-
     &      4*B*E/N-2*E**2/N+2*F**2/N-
     &      2*B*(XA-XM)*(XB-XM)-2*E*(XA-XM)*(XB-XM)-
     &      2*C*(XB-XM)*(YA-YM)-2*C*(XA-XM)*(YB-YM)-
     &      2*A*(YA-YM)*(YB-YM)-2*B*(YA-YM)*(YB-YM)-
     &      4*D*(XB-XM)*(ZA-ZM)-4*F*(YB-YM)*(ZA-ZM)-
     &      4*D*(XA-XM)*(ZB-ZM)-4*F*(YA-YM)*(ZB-ZM)+
     &      8*A*(ZA-ZM)*(ZB-ZM)+8*E*(ZA-ZM)*(ZB-ZM))/BOTTOM**2-
     &  4*(A+2*B+E)*(ZB-ZM)*
     &    (-(B*D*XA)-D*E*XA-C*F*XA+B*D*XM+D*E*XM+C*F*XM-C*D*YA-
     &      A*F*YA-B*F*YA+C*D*YM+A*F*YM+B*F*YM+A**2*ZA+2*A*B*ZA-
     &      D**2*ZA+2*A*E*ZA+2*B*E*ZA+E**2*ZA-F**2*ZA-A**2*ZM-
     &      2*A*B*ZM+D**2*ZM-2*A*E*ZM-2*B*E*ZM-E**2*ZM+F**2*ZM)/
     &   BOTTOM**2-4*(A+2*B+E)*(ZA-ZM)*
     &    (-(B*D*XB)-D*E*XB-C*F*XB+B*D*XM+D*E*XM+C*F*XM-C*D*YB-
     &      A*F*YB-B*F*YB+C*D*YM+A*F*YM+B*F*YM+A**2*ZB+2*A*B*ZB-
     &      D**2*ZB+2*A*E*ZB+2*B*E*ZB+E**2*ZB-F**2*ZB-A**2*ZM-
     &      2*A*B*ZM+D**2*ZM-2*A*E*ZM-2*B*E*ZM-E**2*ZM+F**2*ZM)/
     &   BOTTOM**2+8*TOP*(-(B*D*XA)-D*E*XA-C*F*XA+B*D*XM+D*E*XM+
     &      C*F*XM-C*D*YA-A*F*YA-B*F*YA+C*D*YM+A*F*YM+B*F*YM+
     &      A**2*ZA+2*A*B*ZA-D**2*ZA+2*A*E*ZA+2*B*E*ZA+E**2*ZA-
     &      F**2*ZA-A**2*ZM-2*A*B*ZM+D**2*ZM-2*A*E*ZM-2*B*E*ZM-
     &      E**2*ZM+F**2*ZM)*(-(B*D*XB)-D*E*XB-C*F*XB+B*D*XM+
     &      D*E*XM+C*F*XM-C*D*YB-A*F*YB-B*F*YB+C*D*YM+A*F*YM+
     &      B*F*YM+A**2*ZB+2*A*B*ZB-D**2*ZB+2*A*E*ZB+2*B*E*ZB+
     &      E**2*ZB-F**2*ZB-A**2*ZM-2*A*B*ZM+D**2*ZM-2*A*E*ZM-
     &      2*B*E*ZM-E**2*ZM+F**2*ZM)/BOTTOM**3)/(2.0D0*MASS)
            HESS(J4+3,J3+3)=HESS(J3+3,J4+3)
         ENDDO
      ENDDO
C
C  CASE III: SAME ATOM, DIFFERENT CARTESIAN COORDINATE
C
      DO J1=1,N
         J3=3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         HESS(J3+1,J3+2)=HESS(J3+1,J3+2)+ZJ**2*(
     &  (-2*C+2*C/N+2*XA*YA-2*XM*YA-2*XA*YM+2*XM*YM)/BOTTOM-
     &  TOP*(-2*A*C*XOM-2*C*E*XOM-
     &      2*D*F*XOM-4*C*(XA-XM)**2+
     &      6*A*(XA-XM)*(YA-YM)+8*B*(XA-XM)*(YA-YM)+
     &      6*E*(XA-XM)*(YA-YM)-4*C*(YA-YM)**2-
     &      2*F*(XA-XM)*(ZA-ZM)-2*D*(YA-YM)*(ZA-ZM)-
     &      2*C*(ZA-ZM)**2)/BOTTOM**2-
     &  4*(A*XA+B*XA-A*XM-B*XM-C*YA+C*YM)*
     &    (-(A*C*XA)-C*E*XA-D*F*XA+A*C*XM+C*E*XM+D*F*XM+
     &      2*A*B*YA+B**2*YA-C**2*YA+2*A*E*YA+2*B*E*YA+E**2*YA-
     &      F**2*YA-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM-2*B*E*YM-
     &      E**2*YM+F**2*YM-C*D*ZA-A*F*ZA-B*F*ZA+C*D*ZM+A*F*ZM+
     &      B*F*ZM)/BOTTOM**2+
     &  4*(C*XA-C*XM-B*YA-E*YA+B*YM+E*YM)*
     &    (A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+2*A*E*XA+
     &      2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+D**2*XM-
     &      2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+A*C*YM+
     &      C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+D*E*ZM+
     &      C*F*ZM)/BOTTOM**2+
     &  8*TOP*(-(A*C*XA)-C*E*XA-D*F*XA+A*C*XM+C*E*XM+D*F*XM+
     &      2*A*B*YA+B**2*YA-C**2*YA+2*A*E*YA+2*B*E*YA+E**2*YA-
     &      F**2*YA-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM-2*B*E*YM-
     &      E**2*YM+F**2*YM-C*D*ZA-A*F*ZA-B*F*ZA+C*D*ZM+A*F*ZM+
     &      B*F*ZM)*(A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+
     &      2*A*E*XA+2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+
     &      D**2*XM-2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+
     &      A*C*YM+C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+
     &      D*E*ZM+C*F*ZM)/BOTTOM**3)/(2.0D0*MASS)
         HESS(J3+2,J3+1)=HESS(J3+1,J3+2)

         HESS(J3+1,J3+3)=HESS(J3+1,J3+3)+ZJ**2*(
     & -(TOP*(-2*B*D*XOM-2*D*E*XOM-2*C*F*XOM-
     &        4*D*(XA-XM)**2-2*F*(XA-XM)*(YA-YM)-2*D*(YA-YM)**2+
     &        8*A*(XA-XM)*(ZA-ZM)+6*B*(XA-XM)*(ZA-ZM)+
     &        6*E*(XA-XM)*(ZA-ZM)-2*C*(YA-YM)*(ZA-ZM)-
     &        4*D*(ZA-ZM)**2)/BOTTOM**2)+4*(XA-XM)*(ZA-ZM)/BOTTOM-
     &  4*(A+2*B+E)*(ZA-ZM)*
     &    (A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+2*A*E*XA+
     &      2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+D**2*XM-
     &      2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+A*C*YM+
     &      C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+D*E*ZM+
     &      C*F*ZM)/BOTTOM**2+
     &  4*(A*XA+B*XA-A*XM-B*XM-C*YA+C*YM)*
     &    (B*D*XA+D*E*XA+C*F*XA-B*D*XM-D*E*XM-C*F*XM+C*D*YA+
     &      A*F*YA+B*F*YA-C*D*YM-A*F*YM-B*F*YM-A**2*ZA-2*A*B*ZA+
     &      D**2*ZA-2*A*E*ZA-2*B*E*ZA-E**2*ZA+F**2*ZA+A**2*ZM+
     &      2*A*B*ZM-D**2*ZM+2*A*E*ZM+2*B*E*ZM+E**2*ZM-F**2*ZM)/
     &   BOTTOM**2+8*TOP*(A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-
     &      D**2*XA+2*A*E*XA+2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+
     &      C**2*XM+D**2*XM-2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-
     &      D*F*YA+A*C*YM+C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+
     &      B*D*ZM+D*E*ZM+C*F*ZM)*
     &    (-(B*D*XA)-D*E*XA-C*F*XA+B*D*XM+D*E*XM+C*F*XM-C*D*YA-
     &      A*F*YA-B*F*YA+C*D*YM+A*F*YM+B*F*YM+A**2*ZA+2*A*B*ZA-
     &      D**2*ZA+2*A*E*ZA+2*B*E*ZA+E**2*ZA-F**2*ZA-A**2*ZM-
     &      2*A*B*ZM+D**2*ZM-2*A*E*ZM-2*B*E*ZM-E**2*ZM+F**2*ZM)/
     &   BOTTOM**3 )/(2.0D0*MASS)
         HESS(J3+3,J3+1)=HESS(J3+1,J3+3)

         HESS(J3+2,J3+3)=HESS(J3+2,J3+3)+ZJ**2*(
     &  -(TOP*(-2*C*D*XOM-2*A*F*XOM-2*B*F*XOM-
     &        2*F*(XA-XM)**2-2*D*(XA-XM)*(YA-YM)-4*F*(YA-YM)**2-
     &        2*C*(XA-XM)*(ZA-ZM)+6*A*(YA-YM)*(ZA-ZM)+
     &        6*B*(YA-YM)*(ZA-ZM)+8*E*(YA-YM)*(ZA-ZM)-
     &        4*F*(ZA-ZM)**2)/BOTTOM**2)+4*(YA-YM)*(ZA-ZM)/BOTTOM+   
     &  4*(A+2*B+E)*(ZA-ZM)*
     &    (A*C*XA+C*E*XA+D*F*XA-A*C*XM-C*E*XM-D*F*XM-2*A*B*YA-
     &      B**2*YA+C**2*YA-2*A*E*YA-2*B*E*YA-E**2*YA+F**2*YA+
     &      2*A*B*YM+B**2*YM-C**2*YM+2*A*E*YM+2*B*E*YM+E**2*YM-
     &      F**2*YM+C*D*ZA+A*F*ZA+B*F*ZA-C*D*ZM-A*F*ZM-B*F*ZM)/
     &   BOTTOM**2+4*(-(C*XA)+C*XM+B*YA+E*YA-B*YM-E*YM)*
     &    (B*D*XA+D*E*XA+C*F*XA-B*D*XM-D*E*XM-C*F*XM+C*D*YA+
     &      A*F*YA+B*F*YA-C*D*YM-A*F*YM-B*F*YM-A**2*ZA-2*A*B*ZA+
     &      D**2*ZA-2*A*E*ZA-2*B*E*ZA-E**2*ZA+F**2*ZA+A**2*ZM+
     &      2*A*B*ZM-D**2*ZM+2*A*E*ZM+2*B*E*ZM+E**2*ZM-F**2*ZM)/
     &   BOTTOM**2+8*TOP*(-(A*C*XA)-C*E*XA-D*F*XA+A*C*XM+C*E*XM+
     &      D*F*XM+2*A*B*YA+B**2*YA-C**2*YA+2*A*E*YA+2*B*E*YA+
     &      E**2*YA-F**2*YA-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM-
     &      2*B*E*YM-E**2*YM+F**2*YM-C*D*ZA-A*F*ZA-B*F*ZA+
     &      C*D*ZM+A*F*ZM+B*F*ZM)*
     &    (-(B*D*XA)-D*E*XA-C*F*XA+B*D*XM+D*E*XM+C*F*XM-C*D*YA-
     &      A*F*YA-B*F*YA+C*D*YM+A*F*YM+B*F*YM+A**2*ZA+2*A*B*ZA-
     &      D**2*ZA+2*A*E*ZA+2*B*E*ZA+E**2*ZA-F**2*ZA-A**2*ZM-
     &      2*A*B*ZM+D**2*ZM-2*A*E*ZM-2*B*E*ZM-E**2*ZM+F**2*ZM)/
     &   BOTTOM**3 )/(2.0D0*MASS)
         HESS(J3+3,J3+2)=HESS(J3+2,J3+3)
      ENDDO
C
C  CASE IV: DIFFERENT ATOMS AND DIFFERENT CARTESIAN COORDINATES.
C
      DO J1=1,N
         J3=3*(J1-1)
         XA=X(J3+1)
         YA=X(J3+2)
         ZA=X(J3+3)
         DO J2=1,N
            IF (J2.NE.J1) THEN
            J4=3*(J2-1)
            XB=X(J4+1)
            YB=X(J4+2)
            ZB=X(J4+3)
            HESS(J3+1,J4+2)=HESS(J3+1,J4+2)+ZJ**2*(
     & (2*C/N-2*(XB-XM)*(YA-YM)+4*(XA-XM)*(YB-YM))/BOTTOM- 
     &  TOP*(2*A*C/N+2*C*E/N+2*D*F/N-
     &      4*C*(XA-XM)*(XB-XM)-2*A*(XB-XM)*(YA-YM)-
     &      2*E*(XB-XM)*(YA-YM)+8*A*(XA-XM)*(YB-YM)+
     &      8*B*(XA-XM)*(YB-YM)+8*E*(XA-XM)*(YB-YM)-
     &      4*C*(YA-YM)*(YB-YM)-2*F*(XB-XM)*(ZA-ZM)-
     &      2*D*(YA-YM)*(ZB-ZM)-2*C*(ZA-ZM)*(ZB-ZM))/BOTTOM**2-
     &  4*(A*XA+B*XA-A*XM-B*XM-C*YA+C*YM)*
     &    (-(A*C*XB)-C*E*XB-D*F*XB+A*C*XM+C*E*XM+D*F*XM+
     &      2*A*B*YB+B**2*YB-C**2*YB+2*A*E*YB+2*B*E*YB+E**2*YB-
     &      F**2*YB-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM-2*B*E*YM-
     &      E**2*YM+F**2*YM-C*D*ZB-A*F*ZB-B*F*ZB+C*D*ZM+A*F*ZM+
     &      B*F*ZM)/BOTTOM**2+
     &  4*(C*XB-C*XM-B*YB-E*YB+B*YM+E*YM)*
     &    (A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+2*A*E*XA+
     &      2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+D**2*XM-
     &      2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+A*C*YM+
     &      C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+D*E*ZM+
     &      C*F*ZM)/BOTTOM**2+
     &  8*TOP*(-(A*C*XB)-C*E*XB-D*F*XB+A*C*XM+C*E*XM+D*F*XM+
     &      2*A*B*YB+B**2*YB-C**2*YB+2*A*E*YB+2*B*E*YB+E**2*YB-
     &      F**2*YB-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM-2*B*E*YM-
     &      E**2*YM+F**2*YM-C*D*ZB-A*F*ZB-B*F*ZB+C*D*ZM+A*F*ZM+
     &      B*F*ZM)*(A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+
     &      2*A*E*XA+2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+
     &      D**2*XM-2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+
     &      A*C*YM+C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+
     &      D*E*ZM+C*F*ZM)/BOTTOM**3 )/(2.0D0*MASS)
            HESS(J4+2,J3+1)=HESS(J3+1,J4+2)

            HESS(J3+1,J4+3)=HESS(J3+1,J4+3)+ZJ**2*(
     &  -(TOP*(2*B*D/N+2*D*E/N+2*C*F/N-4*D*(XA-XM)*(XB-XM) -
     &        2*F*(XB-XM)*(YA-YM)-2*D*(YA-YM)*(YB-YM) -
     &        2*B*(XB-XM)*(ZA-ZM)-2*E*(XB-XM)*(ZA-ZM) -
     &        2*C*(YB-YM)*(ZA-ZM)+8*A*(XA-XM)*(ZB-ZM) +
     &        8*B*(XA-XM)*(ZB-ZM)+8*E*(XA-XM)*(ZB-ZM) -
     &        4*D*(ZA-ZM)*(ZB-ZM))/BOTTOM**2) +
     &  4*(XA-XM)*(ZB-ZM)/BOTTOM -
     &  4*(A+2*B+E)*(ZB-ZM)*
     &    (A**2*XA+2*A*B*XA+B**2*XA-C**2*XA-D**2*XA+2*A*E*XA +
     &      2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM+C**2*XM+D**2*XM -
     &      2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA-D*F*YA+A*C*YM +
     &      C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA+B*D*ZM+D*E*ZM +
     &      C*F*ZM)/BOTTOM**2 +
     &  4*(A*XA+B*XA-A*XM-B*XM-C*YA+C*YM)*
     &    (B*D*XB+D*E*XB+C*F*XB-B*D*XM-D*E*XM-C*F*XM+C*D*YB +
     &      A*F*YB+B*F*YB-C*D*YM-A*F*YM-B*F*YM-A**2*ZB-2*A*B*ZB +
     &      D**2*ZB-2*A*E*ZB-2*B*E*ZB-E**2*ZB+F**2*ZB+A**2*ZM +
     &      2*A*B*ZM-D**2*ZM+2*A*E*ZM+2*B*E*ZM+E**2*ZM-F**2*ZM)/
     &   BOTTOM**2+8*TOP*(A**2*XA+2*A*B*XA+B**2*XA-C**2*XA -
     &      D**2*XA+2*A*E*XA+2*B*E*XA-A**2*XM-2*A*B*XM-B**2*XM +
     &      C**2*XM+D**2*XM-2*A*E*XM-2*B*E*XM-A*C*YA-C*E*YA -
     &      D*F*YA+A*C*YM+C*E*YM+D*F*YM-B*D*ZA-D*E*ZA-C*F*ZA +
     &      B*D*ZM+D*E*ZM+C*F*ZM)*
     &    (-(B*D*XB)-D*E*XB-C*F*XB+B*D*XM+D*E*XM+C*F*XM-C*D*YB -
     &      A*F*YB-B*F*YB+C*D*YM+A*F*YM+B*F*YM+A**2*ZB+2*A*B*ZB -
     &      D**2*ZB+2*A*E*ZB+2*B*E*ZB+E**2*ZB-F**2*ZB-A**2*ZM -
     &      2*A*B*ZM+D**2*ZM-2*A*E*ZM-2*B*E*ZM-E**2*ZM+F**2*ZM)/
     &   BOTTOM**3 )/(2.0D0*MASS)
            HESS(J4+3,J3+1)=HESS(J3+1,J4+3)

            HESS(J3+2,J4+3)=HESS(J3+2,J4+3)+ZJ**2*(
     &   -(TOP*(2*C*D/N+2*A*F/N+2*B*F/N-2*F*(XA-XM)*(XB-XM) -
     &        2*D*(XA-XM)*(YB-YM)-4*F*(YA-YM)*(YB-YM) -
     &        2*C*(XB-XM)*(ZA-ZM)-2*A*(YB-YM)*(ZA-ZM) -
     &        2*B*(YB-YM)*(ZA-ZM)+8*A*(YA-YM)*(ZB-ZM) +
     &        8*B*(YA-YM)*(ZB-ZM)+8*E*(YA-YM)*(ZB-ZM) -
     &        4*F*(ZA-ZM)*(ZB-ZM))/BOTTOM**2) +
     &  4*(YA-YM)*(ZB-ZM)/BOTTOM +
     &  4*(A+2*B+E)*(ZB-ZM)*
     &    (A*C*XA+C*E*XA+D*F*XA-A*C*XM-C*E*XM-D*F*XM-2*A*B*YA -
     &      B**2*YA+C**2*YA-2*A*E*YA-2*B*E*YA-E**2*YA+F**2*YA +
     &      2*A*B*YM+B**2*YM-C**2*YM+2*A*E*YM+2*B*E*YM+E**2*YM -
     &      F**2*YM+C*D*ZA+A*F*ZA+B*F*ZA-C*D*ZM-A*F*ZM-B*F*ZM)/
     &   BOTTOM**2+4*(-(C*XA)+C*XM+B*YA+E*YA-B*YM-E*YM)*
     &    (B*D*XB+D*E*XB+C*F*XB-B*D*XM-D*E*XM-C*F*XM+C*D*YB +
     &      A*F*YB+B*F*YB-C*D*YM-A*F*YM-B*F*YM-A**2*ZB-2*A*B*ZB +
     &      D**2*ZB-2*A*E*ZB-2*B*E*ZB-E**2*ZB+F**2*ZB+A**2*ZM +
     &      2*A*B*ZM-D**2*ZM+2*A*E*ZM+2*B*E*ZM+E**2*ZM-F**2*ZM)/
     &   BOTTOM**2+8*TOP*(-(A*C*XA)-C*E*XA-D*F*XA+A*C*XM+C*E*XM +
     &      D*F*XM+2*A*B*YA+B**2*YA-C**2*YA+2*A*E*YA+2*B*E*YA +
     &      E**2*YA-F**2*YA-2*A*B*YM-B**2*YM+C**2*YM-2*A*E*YM -
     &      2*B*E*YM-E**2*YM+F**2*YM-C*D*ZA-A*F*ZA-B*F*ZA +
     &      C*D*ZM+A*F*ZM+B*F*ZM)*
     &    (-(B*D*XB)-D*E*XB-C*F*XB+B*D*XM+D*E*XM+C*F*XM-C*D*YB -
     &      A*F*YB-B*F*YB+C*D*YM+A*F*YM+B*F*YM+A**2*ZB+2*A*B*ZB -
     &      D**2*ZB+2*A*E*ZB+2*B*E*ZB+E**2*ZB-F**2*ZB-A**2*ZM -
     &      2*A*B*ZM+D**2*ZM-2*A*E*ZM-2*B*E*ZM-E**2*ZM+F**2*ZM)/
     &   BOTTOM**3 )/(2.0D0*MASS)
            HESS(J4+3,J3+2)=HESS(J3+2,J4+3)
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END
C
C  ENERGY FOR LINEAR SYSTEM
C
      SUBROUTINE EPANIC(X,N,ZJ,ROTE,FLAG,V,MASS)
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, N, J2
      LOGICAL FLAG
      DOUBLE PRECISION X(3*N), ROTE, ZJ, 
     1                 V(3*N), TEMP2, TEMP1,DIF,
     4                 V1, V2, V3, V4, MASS
      CALL NROTE(X,N,ROTE,MASS,ZJ)
      IF (FLAG) RETURN
      DIF=1.0D-4
      DO 20 J1=1,3*N
         TEMP1=X(J1)
         X(J1)=X(J1)+DIF
C        PRINT*,'J1,X=',J1,X(J1)
         CALL NROTE(X,N,V1,MASS,ZJ)
         X(J1)=X(J1)-2.0D0*DIF
C        PRINT*,'J1,X=',J1,X(J1)
         CALL NROTE(X,N,V2,MASS,ZJ)
         V(J1)=(V1-V2)/(2.0D0*DIF)
         X(J1)=TEMP1
C        PRINT*,'J1,GRAD=',J1,V(J1)
         DO 10 J2=J1,3*N
            TEMP1=X(J1)
            TEMP2=X(J2)
            X(J1)=X(J1)+DIF
            X(J2)=X(J2)+DIF
C           PRINT*,'J1,J2,X=',J1,J2,X(J1),X(J2)
            CALL NROTE(X,N,V1,MASS,ZJ)
            X(J1)=X(J1)-2.0D0*DIF
C           PRINT*,'J1,J2,X=',J1,J2,X(J1),X(J2)
            CALL NROTE(X,N,V2,MASS,ZJ)
            X(J2)=X(J2)-2.0D0*DIF
C           PRINT*,'J1,J2,X=',J1,J2,X(J1),X(J2)
            CALL NROTE(X,N,V3,MASS,ZJ)
            X(J1)=X(J1)+2.0D0*DIF
C           PRINT*,'J1,J2,X=',J1,J2,X(J1),X(J2)
            CALL NROTE(X,N,V4,MASS,ZJ)
            X(J1)=TEMP1
            X(J2)=TEMP2
            HESS(J1,J2)=(V1-V2-V4+V3)/(4.0D0*DIF*DIF)
            HESS(J2,J1)=HESS(J1,J2)
C           PRINT*,'V1,V2,V3,V4=',V1,V2,V3,V4
C           PRINT*,'J1,J2,HESS=',J1,J2,HESS(J1,J2)
10       CONTINUE
20    CONTINUE

      RETURN
      END

      SUBROUTINE NROTE(X,N,ROTE,MASS,ZJ)
      IMPLICIT NONE
      INTEGER J1, N, J3
      DOUBLE PRECISION MASS, X(3*N), ROTE, ZJ,
     1                 XM, YM, ZM, TOP, BOTTOM, 
     2                 A, B, C, D, E, F, XT(N)
C
C  CALCULATE CENTRE OF MASS COMPONENTS XM, YM, ZM
C
      XM=0.0D0
      YM=0.0D0
      ZM=0.0D0
      DO J1=1,N
         XM=XM+X(3*(J1-1)+1)
         YM=YM+X(3*(J1-1)+2)
         ZM=ZM+X(3*(J1-1)+3)
      ENDDO
      XM=XM/N
      YM=YM/N
      ZM=ZM/N
C
C  CALCULATE CONSTANTS A, B, C, D, E, F
C
      A=0.0D0
      B=0.0D0
      C=0.0D0
      D=0.0D0
      E=0.0D0
      F=0.0D0
      DO J1=1,N
         J3=3*(J1-1)
         A=A+(X(J3+2)-YM)**2
         B=B+(X(J3+3)-ZM)**2
         E=E+(X(J3+1)-XM)**2
         C=C+(X(J3+2)-YM)*(X(J3+1)-XM)
         D=D+(X(J3+3)-ZM)*(X(J3+1)-XM)
         F=F+(X(J3+3)-ZM)*(X(J3+2)-YM)
      ENDDO
      TOP=A*B+B*B-C*C+A*E+B*E
      BOTTOM=A*A*B+A*B*B-A*C*C-B*D*D+A*A*E+2.0D0*A*B*E+B*B*E-C*C*E-D*D*E
     1      +A*E*E+B*E*E-2.0D0*C*D*F-A*F*F-B*F*F
C     IF (DABS(BOTTOM).LT.1.0D-10) THEN
         DO J1=1,N
            IF (X(3*(J1-1)+1).NE.0.0D0) THEN
               XT(J1)=DSQRT(X(3*(J1-1)+1)**2
     1         +X(3*(J1-1)+2)**2)*X(3*(J1-1)+1)/DABS(X(3*(J1-1)+1))
            ELSE
               XT(J1)=0.0D0
            ENDIF
C           PRINT*,'J1,XT=',J1,XT(J1)
         ENDDO
         XM=0.0D0
         DO J1=1,N
            XM=XM+XT(J1)
         ENDDO
         XM=XM/N
         E=0.0D0
         DO J1=1,N
            E=E+(XT(J1)-XM)**2
         ENDDO
         ROTE=ZJ**2/(2.0D0*E*MASS)
C        PRINT*,'PANIC ROTE, E, BOTTOM=',ROTE,E, BOTTOM
C     ELSE
C        IZINV=TOP/(MASS*BOTTOM)
C        ROTE=ZJ**2*IZINV/2.0D0
C        PRINT*,'NORMAL ROTE, BOTTOM=',ROTE, BOTTOM
C     ENDIF
      RETURN
      END
C
C**************************************************************************
C
C  SUBROUTINE ROTDERIV CALCULATES ANALYTICALLY THE CARTESIAN GRADIENT AND
C  SECOND DERIVATIVE MATRIX DUE TO THE ROTATIONAL TERM ONLY AND INCREMENTS
C  THE EXISTING MATRICES CONTAINING THE POTENTIAL CONTRIBUTIONS
C  ACCORDINGLY.
C
C**************************************************************************
C
      SUBROUTINE ROTDERIV(N, X, V, MASS, AVEL, IZ)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), MASS, AVEL, IZ
C
C  ANALYTIC GRADIENT:
C
      DO J1 = 1, N
         DO J2 = 1, 2
            J3 = 3*(J1-1)+J2
            V(J3) = V(J3) + (AVEL**2) * X(J3) * MASS
         ENDDO
      ENDDO
C
C  CASE I: ENTIRELY DIAGONAL TERMS.
C
      DO J1 = 1, N
         DO J2 = 1, 2
            J3 = 3*(J1-1)+J2
            HESS(J3, J3) = HESS(J3, J3) + (AVEL**2) * MASS
         ENDDO
      ENDDO
      RETURN
      END
C
C*************************************************************************
C
C  THIS SUBROUTINE CALCULATES THE MOMENT OF INERTIA IZ ABOUT THE Z AXIS
C  AND HENCE THE ROTATIONAL ENERGY
C
C  ENERGY = OMEGA^2 * I_{ZZ} / 2 
C
C*************************************************************************
C
      SUBROUTINE ROTENERGY(N, X, AVEL, MASS, IZ, ROT)
      IMPLICIT NONE
      INTEGER N, I
      DOUBLE PRECISION X(3*N), AVEL, ROT, IZ, MASS
      IZ = 0.0D0
      DO 10 I = 0, N-1
         IZ = IZ + MASS * ( X(I*3+1)**2 + X(I*3+2)**2 )
10    CONTINUE
      ROT = 0.5D0 * IZ * AVEL**2
      RETURN
      END
