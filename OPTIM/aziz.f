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
C***********************************************************************
C
C  SUBROUTINE AZIZ CALCULATES THE CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY FOR THE ACCURATE AR-AR PAIR POTENTIAL
C  OF AZIZ (JCP, 99, 4518, 1993). EPSILON AND R_E ARE SET TO UNITY.
C  INPUT COORDINATES ARE ASSUMED TO BE IN SIGMA=R_E/2**(1/6).
C
C***********************************************************************
C
      SUBROUTINE AZIZ(N, X, V, ENERGY, NTEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NTEST
      DOUBLE PRECISION X(3*N), CONV, SR2, SR3, SR6, 
     1                 SR10, SR14,
     1              V(3*N), R2(N,N), ENERGY, R, R6(N,N),
     2              R8(N,N), G(N,N), R1(N,N), R10(N,N),
     3              R14(N,N), F(N,N), R12(N,N)
C     PARAMETER (CONV=2.0D0**(1.0D0/6.0D0), SR14=DSQRT(14.0D0), 
C    1           SR3=DSQRT(3.0D0))
C     PARAMETER (SR10=DSQRT(10.0D0), SR2=DSQRT(2.0D0), SR6=DSQRT(6.0D0))
      PARAMETER (SR10=3.16227766D0, SR2=1.414213562D0, SR6=2.44948743D0)
      PARAMETER (CONV=1.122462048D0, SR14=3.741657387D0, 
     1           SR3=1.732050808D0)
      DOUBLE PRECISION AA, ALPHA, BETA, C6, C8, C10, C12, C14, RHO, RM,
     1                 ATOB, Z1, Z2, Z3, RC1, RC2, RC3, RC4, RC5
      PARAMETER (AA= 8.73933927D4, ALPHA=9.0328328D0, 
     1           BETA=-2.37132823D0,
     2           C6=1.09309955D0, C8=0.51568309D0, C10=0.32521242D0, 
     3           C12=0.27818156D0, C14=0.31111959D0,
     4           RHO=1.107D0, RM=3.757D0, ATOB=1.889726164D0, 
     5           Z1=2.1D0, Z2=0.109D0, Z3=0.78D0)
C
C  CONVERT DISTANCES TO R_E RATHER THAN SIGMA.
C
      DO J1=1,3*N
         X(J1)=X(J1)/CONV
      ENDDO
C 
C  STORE DISTANCE MATRICES.
C
      DO 20 J1=1,N
         R2(J1,J1)=0.0D0
         R8(J1,J1)=0.0D0
         R14(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
         R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1               +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2               +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R1(J2,J1)=DSQRT(R2(J2,J1))
            R2(J2,J1)=1.0D0/R2(J2,J1)
            R8(J2,J1)=R2(J2,J1)**4
            R6(J2,J1)=R8(J2,J1)/R2(J2,J1)
            R10(J2,J1)=R8(J2,J1)*R2(J2,J1)
            R12(J2,J1)=R10(J2,J1)*R2(J2,J1)
            R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)

            R1(J1,J2)=R1(J2,J1)
            R2(J1,J2)=R2(J2,J1)
            R6(J1,J2)=R6(J2,J1)
            R8(J1,J2)=R8(J2,J1)
            R10(J1,J2)=R10(J2,J1)
            R12(J1,J2)=R12(J2,J1)
            R14(J1,J2)=R14(J2,J1)
C
10       CONTINUE
20    CONTINUE
C
C  CALCULATE THE ENERGY IN EPSILON.
C
      ENERGY=0.0D0
      DO J1=1,N
         DO J2=J1+1,N
            R=R1(J2,J1)
            RC1=ATOB*R*RHO*RM*Z1
            RC2=ATOB**2*R**2*RHO**2*RM**2*Z2
            RC3=ATOB*R*RHO*RM
            ENERGY=ENERGY + AA*DEXP(-(ALPHA*R) + BETA*R**2) - 
     1             (1 - RC3**1.68/DEXP(RC3*Z3))*
     2   (C10*(1 - DEXP(-RC1/10 - RC2/SR10))**10*R10(J2,J1) + 
     3   C12*(1 - DEXP(-RC1/12 - RC2/(2*SR3)))**12*R12(J2,J1) + 
     4    C14*(1 - DEXP(-RC1/14 - RC2/SR14))**14*R14(J2,J1) 
     5   + C6*(1 - DEXP(-RC1/6 - RC2/SR6))**6*R6(J2,J1) + 
     6     C8*(1 - DEXP(-RC1/8 - RC2/(2*SR2)))**8*R8(J2,J1))
         ENDDO
      ENDDO

      IF (NTEST.EQ.0) GOTO 11
C
C  CALCULATE THE G AND F TENSORS.
C
      DO 21 J1=1,N
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO 22 J2=J1+1,N 
            R=R1(J2,J1)
            RC1=ATOB*R*RHO*RM*Z1
            RC2=ATOB**2*R**2*RHO**2*RM**2*Z2
            RC3=ATOB*R*RHO*RM
            RC4=ATOB*RHO*RM*Z1
            RC5=ATOB**2*R*RHO**2*RM**2

            G(J2,J1)= (AA*DEXP(-(ALPHA*R) + BETA*R**2)*(-ALPHA 
     %       + 2*BETA*R) - 
     %                (1 - RC3**1.68/DEXP(RC3*Z3))*
     %     (-14*C14*(1 - DEXP(-RC1/14 - RC2/SR14))**14/R**15 - 
     %       12*C12*(1 - DEXP(-RC1/12 - RC2/(2*SR3)))**12/R**13 - 
     %       10*C10*(1 - DEXP(-RC1/10 - RC2/SR10))**10/R**11 - 
     %       8*C8*(1 - DEXP(-RC1/8 - RC2/(2*SR2)))**8/R**9 - 
     %       6*C6*(1 - DEXP(-RC1/6 - RC2/SR6))**6/R**7 - 
     %       14*C14*DEXP(-RC1/14 - RC2/SR14)*(1 - DEXP(-RC1/14 
     %       - RC2/SR14))**13*
     %        R14(J2,J1)*(-RC4/14 - DSQRT(2.0D0/7.0D0)*RC5*Z2) - 
     %       10*C10*DEXP(-RC1/10 - RC2/SR10)*(1 - DEXP(-RC1/10 
     %       - RC2/SR10))**9*
     %        R10(J2,J1)*(-RC4/10 - DSQRT(2.0D0/5.0D0)*RC5*Z2) - 
     %       6*C6*DEXP(-RC1/6 - RC2/SR6)*(1 - DEXP(-RC1/6 
     %       - RC2/SR6))**5*R6(J2,J1)*
     %        (-RC4/6 - DSQRT(2.0D0/3.0D0)*RC5*Z2) - 
     %       8*C8*DEXP(-RC1/8 - RC2/(2*SR2))*
     %        (1 - DEXP(-RC1/8 - RC2/(2*SR2)))**7*R8(J2,J1)*(-RC4/8 
     %        - RC5*Z2/SR2) - 
     %       12*C12*DEXP(-RC1/12 - RC2/(2*SR3))*
     %        (1 - DEXP(-RC1/12 - RC2/(2*SR3)))**11*R12(J2,J1)
     %       *(-RC4/12 - RC5*Z2/SR3))
     %      - (C10*(1 - DEXP(-RC1/10 - RC2/SR10))**10*R10(J2,J1) + 
     %       C12*(1 - DEXP(-RC1/12 - RC2/(2*SR3)))**12*R12(J2,J1) + 
     %       C14*(1 - DEXP(-RC1/14 - RC2/SR14))**14*R14(J2,J1) + 
     %       C6*(1 - DEXP(-RC1/6 - RC2/SR6))**6*R6(J2,J1) + 
     %       C8*(1 - DEXP(-RC1/8 - RC2/(2*SR2)))**8*R8(J2,J1))*
     %     (-1.68*ATOB*RC3**0.68*RHO*RM/DEXP(RC3*Z3) + 
     %       ATOB*RC3**1.68*RHO*RM*Z3/DEXP(RC3*Z3)))/R

            F(J2,J1)= 2*AA*BETA*DEXP(-(ALPHA*R) + BETA*R**2) + 
     -  AA*DEXP(-(ALPHA*R) + BETA*R**2)*(-ALPHA + 2*BETA*R)**2 - 
     -  (1.0D0 - RC3**1.68/DEXP(RC3*Z3))*
     -   (210*C14*(1.0D0 - DEXP(-RC1/14.0D0 - RC2/SR14))**14/R**16 + 
     -   72*C8*(1.0D0 - DEXP(-RC1/8.0D0 - RC2/(2*SR2)))**8*R10(J2,J1) + 
     -   110*C10*(1.0D0 - DEXP(-RC1/10.0D0 - RC2/SR10))**10*R12(J2,J1) +
     -  156*C12*(1.0D0-DEXP(-RC1/12.0D0 - RC2/(2*SR3)))**12*R14(J2,J1) +
     -     42*C6*(1.0D0 - DEXP(-RC1/6.0D0 - RC2/SR6))**6*R8(J2,J1) + 
     -     2*SR10*ATOB**2*C10*DEXP(-RC1/10.0D0 - RC2/SR10)*
     -      (1.0D0 - DEXP(-RC1/10.0D0 - RC2/SR10))**9
     -    *RHO**2*RM**2*R10(J2,J1)*Z2 + 
     -     4*SR3*ATOB**2*C12*DEXP(-RC1/12.0D0 - RC2/(2*SR3))*
     -      (1.0D0 - DEXP(-RC1/12.0D0 - RC2/(2.0D0*SR3)))**11
     -   *RHO**2*RM**2*R12(J2,J1)*Z2 + 
     -     2*SR14*ATOB**2*C14*DEXP(-RC1/14.0D0 - RC2/SR14)*
     -      (1.0D0 - DEXP(-RC1/14.0D0 - RC2/SR14))**13*RHO**2
     -    *RM**2*R14(J2,J1)*Z2 + 
     -     2*SR6*ATOB**2*C6*DEXP(-RC1/6.0D0 - RC2/SR6)*
     -      (1.0D0 - DEXP(-RC1/6.0D0 - RC2/SR6))**5
     -    *RHO**2*RM**2*R6(J2,J1)*Z2 + 
     -     4*SR2*ATOB**2*C8*DEXP(-RC1/8.0D0 - RC2/(2*SR2))*
     -      (1.0D0 - DEXP(-RC1/8.0D0 - RC2/(2*SR2)))**7*RHO**2
     -     *RM**2*R8(J2,J1)*Z2 + 
     -     392*C14*DEXP(-RC1/14 - RC2/SR14)*(1 - DEXP(-RC1/14 - 
     -     RC2/SR14))**13*
     -       (-RC4/14.0D0 - DSQRT(2.0D0/7.0D0)*RC5*Z2)/R**15 + 
     -     182*C14*DEXP(-RC1/7.0D0 - DSQRT(2.0D0/7.0D0)*RC2)*
     -      (1.0D0 - DEXP(-RC1/14.0D0 - RC2/SR14))**12*R14(J2,J1)*
     -      (-RC4/14.0D0 - DSQRT(2.0D0/7.0D0)*RC5*Z2)**2 - 
     -     14*C14*DEXP(-RC1/14.0D0 - RC2/SR14)*(1.0D0 - DEXP(-RC1/14 
     -     - RC2/SR14))**13*
     -      R14(J2,J1)*(-RC4/14.0D0 - DSQRT(2.0D0/7.0D0)*RC5*Z2)**2 + 
     -     200*C10*DEXP(-RC1/10.0D0 - RC2/SR10)*(1.0D0 - 
     -      DEXP(-RC1/10.0D0 - RC2/SR10))**9*
     -       (-RC4/10.0D0 - DSQRT(2.0D0/5.0D0)*RC5*Z2)/R**11 + 
     -     90*C10*DEXP(-RC1/5.0D0 - DSQRT(2.0D0/5.0D0)*RC2)*
     -      (1.0D0 - DEXP(-RC1/10.0D0 - RC2/SR10))**8*R10(J2,J1)*
     -      (-RC4/10.0D0 - DSQRT(2.0D0/5.0D0)*RC5*Z2)**2 - 
     -     10*C10*DEXP(-RC1/10.0D0 - RC2/SR10)*(1.0D0 - 
     -     DEXP(-RC1/10.0D0 - RC2/SR10))**9*
     -      R10(J2,J1)*(-RC4/10.0D0 - DSQRT(2.0D0/5.0D0)*RC5*Z2)**2 + 
     -     72*C6*DEXP(-RC1/6.0D0 - RC2/SR6)*(1.0D0 - DEXP(-RC1/6.0D0 
     -      - RC2/SR6))**5*
     -       (-RC4/6.0D0 - DSQRT(2.0D0/3.0D0)*RC5*Z2)/R**7 + 
     -     30*C6*DEXP(-RC1/3.0D0 - DSQRT(2.0D0/3.0D0)*RC2)*(1.0D0 
     -     - DEXP(-RC1/6.0D0 - RC2/SR6))**4*
     -      R6(J2,J1)*(-RC4/6.0D0 - DSQRT(2.0D0/3.0D0)*RC5*Z2)**2 - 
     -     6*C6*DEXP(-RC1/6.0D0 - RC2/SR6)*(1.0D0 - DEXP(-RC1/6 
     -     - RC2/SR6))**5*R6(J2,J1)*
     -      (-RC4/6.0D0 - DSQRT(2.0D0/3.0D0)*RC5*Z2)**2 + 
     -     128*C8*DEXP(-RC1/8.0D0 - RC2/(2*SR2))*
     -       (1.0D0 - DEXP(-RC1/8.0D0 - RC2/(2*SR2)))**7*(-RC4/8 
     -      - RC5*Z2/SR2)/R**9 + 
     -     56*C8*DEXP(-RC1/4.0D0 - RC2/SR2)*(1.0D0 - DEXP(-RC1/8 
     -    - RC2/(2*SR2)))**6*
     -      R8(J2,J1)*(-RC4/8.0D0 - RC5*Z2/SR2)**2 - 
     -     8*C8*DEXP(-RC1/8.0D0 - RC2/(2*SR2))*(1.0D0 - DEXP(-RC1/8 
     -     - RC2/(2*SR2)))**7*
     -      R8(J2,J1)*(-RC4/8.0D0 - RC5*Z2/SR2)**2 + 
     -     288*C12*DEXP(-RC1/12.0D0 - RC2/(2*SR3))*
     -       (1.0D0 - DEXP(-RC1/12.0D0 - RC2/(2*SR3)))**11*(-RC4/12 
     -     - RC5*Z2/SR3)/R**13
     -       + 132*C12*DEXP(-RC1/6.0D0 - RC2/SR3)*
     -      (1.0D0 - DEXP(-RC1/12.0D0 - RC2/(2*SR3)))**10*R12(J2,J1)*
     -      (-RC4/12.0D0 - RC5*Z2/SR3)**2 - 
     -     12*C12*DEXP(-RC1/12.0D0 - RC2/(2*SR3))*
     -      (1.0D0 - DEXP(-RC1/12.0D0 - RC2/(2*SR3)))**11*R12(J2,J1)*
     -      (-RC4/12.0D0 - RC5*Z2/SR3)**2) - 
     -  2*(-14*C14*(1.0D0 - DEXP(-RC1/14.0D0 - RC2/SR14))**14/R**15 - 
     -     12*C12*(1.0D0 - DEXP(-RC1/12.0D0 - RC2/(2*SR3)))**12/R**13 - 
     -     10*C10*(1.0D0 - DEXP(-RC1/10.0D0 - RC2/SR10))**10/R**11 - 
     -     8*C8*(1.0D0 - DEXP(-RC1/8.0D0 - RC2/(2*SR2)))**8/R**9 - 
     -     6*C6*(1.0D0 - DEXP(-RC1/6.0D0 - RC2/SR6))**6/R**7 - 
     -     14*C14*DEXP(-RC1/14.0D0 - RC2/SR14)*(1.0D0 - DEXP(-RC1/14 
     -      - RC2/SR14))**13*
     -      R14(J2,J1)*(-RC4/14.0D0 - DSQRT(2.0D0/7.0D0)*RC5*Z2) - 
     -     10*C10*DEXP(-RC1/10.0D0 - RC2/SR10)*(1.0D0 - DEXP(-RC1/10 
     -     - RC2/SR10))**9*
     -      R10(J2,J1)*(-RC4/10.0D0 - DSQRT(2.0D0/5.0D0)*RC5*Z2) - 
     -     6*C6*DEXP(-RC1/6.0D0 - RC2/SR6)*(1.0D0 - DEXP(-RC1/6 - 
     -     RC2/SR6))**5*R6(J2,J1)*
     -      (-RC4/6.0D0 - DSQRT(2.0D0/3.0D0)*RC5*Z2) - 
     -     8*C8*DEXP(-RC1/8.0D0 - RC2/(2*SR2))*(1.0D0 - DEXP(-RC1/8 
     -     - RC2/(2*SR2)))**7*
     -      R8(J2,J1)*(-RC4/8.0D0 - RC5*Z2/SR2) - 
     -     12*C12*DEXP(-RC1/12.0D0 - RC2/(2*SR3))*
     -      (1.0D0 - DEXP(-RC1/12.0D0 - RC2/(2*SR3)))**11*R12(J2,J1)*
     -      (-RC4/12.0D0 - RC5*Z2/SR3))*
     -   (-1.68*ATOB*RC3**0.68*RHO*RM/DEXP(RC3*Z3) + 
     -     ATOB*RC3**1.68*RHO*RM*Z3/DEXP(RC3*Z3)) - 
     -  (C10*(1.0D0 - DEXP(-RC1/10.0D0 - RC2/SR10))**10*R10(J2,J1) + 
     -     C12*(1.0D0-DEXP(-RC1/12.0D0 - RC2/(2*SR3)))**12*R12(J2,J1) + 
     -     C14*(1.0D0 - DEXP(-RC1/14.0D0 - RC2/SR14))**14*R14(J2,J1) + 
     -     C6*(1.0D0 - DEXP(-RC1/6.0D0 - RC2/SR6))**6*R6(J2,J1) + 
     -     C8*(1.0D0 - DEXP(-RC1/8.0D0 - RC2/(2*SR2)))**8*R8(J2,J1))*
     -   (-1.1424*ATOB**2*RHO**2*RM**2/
     -      (DEXP(RC3*Z3)*RC3**0.32) + 
     -     3.36*ATOB**2*RC3**0.68*RHO**2*RM**2*Z3/DEXP(RC3*Z3) - 
     -     ATOB**2*RC3**1.68*RHO**2*RM**2*Z3**2/DEXP(RC3*Z3)) 

            F(J2,J1)=F(J2,J1)-G(J2,J1)

            G(J1,J2)=G(J2,J1)
            F(J1,J2)=F(J2,J1)
22       CONTINUE
21    CONTINUE
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
         V(J1)=V(J1)/CONV
         HESS(J1,J1)=HESS(J1,J1)/(CONV**2)
         DO 190 J2=J1+1,3*N
            HESS(J2,J1)=HESS(J2,J1)/(CONV**2)
            HESS(J1,J2)=HESS(J2,J1)
190      CONTINUE
200   CONTINUE
C
C  CONVERT DISTANCES BACK TO SIGMA.
C
11    DO J1=1,3*N
          X(J1)=X(J1)*CONV
      ENDDO
      RETURN
      END
