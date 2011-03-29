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
      SUBROUTINE CTEST(N, VAR, V, ENERGY, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION VAR(3*N), V(3*N), ENERGY, X, Y, R, PHI, E
      PARAMETER (E=2.718281828459045D0)
      LOGICAL CART, METRIC, EFFHESS, GTEST, STEST
      PARAMETER (CART=.TRUE.,METRIC=.FALSE.,EFFHESS=.FALSE.)

      IF (CART) THEN
         X=VAR(1)
         Y=VAR(2)
         ENERGY=(-2*X*Y)/(E**(-E + SQRT(X**2 + Y**2))**2*(X**2 + Y**2))
         IF (GTEST) THEN
            V(1)=(2*Y*((-2*E*X**2)/(X**2 + Y**2)**1.5 + (2*X**4 - Y**2 + X**2*(1 + 2*Y**2))/(X**2 + Y**2)**2))/
     1  E**(E - SQRT(X**2 + Y**2))**2
            V(2)=(2*X*((-2*E*Y**2)/(X**2 + Y**2)**1.5 + (Y**2 + 2*Y**4 + X**2*(-1 + 2*Y**2))/(X**2 + Y**2)**2))/
     1  E**(E - SQRT(X**2 + Y**2))**2
         ENDIF
         IF (STEST) THEN
            HESS(1,1)=(4*X*Y*((-4*X**2 + (3 - 2*(2 + E**2)*X**2)*(X**2 + Y**2) + (3 - 2*X**2)*(X**2 + Y**2)**2)/
     1       (X**2 + Y**2)**3 + (E*(4*X**4 - 3*Y**2 + X**2*(2 + 4*Y**2)))/(X**2 + Y**2)**2.5))/
     1       E**(E - SQRT(X**2 + Y**2))**2
            HESS(1,2)=(2*((2*E*(-Y**4 + X**4*(-1 + 4*Y**2) + X**2*Y**2*(3 + 4*Y**2)))/(X**2 + Y**2)**2.5 + 
     1      (Y**4 + 2*Y**6 + X**6*(2 - 4*Y**2) - 2*X**2*Y**2*(3 + (1 + 2*E**2)*Y**2 + 2*Y**4) - 
     1         X**4*(-1 + (2 + 4*E**2)*Y**2 + 8*Y**4))/(X**2 + Y**2)**3))/E**(E - SQRT(X**2 + Y**2))**2
            HESS(2,1)=HESS(1,2)
            HESS(2,2)=(4*X*Y*(-((-3*(X**2 + X**4) + (1 + 2*X**2*(-1 + E**2 + X**2))*Y**2 + 
     1           (1 + 2*E**2 + 4*X**2)*Y**4 + 2*Y**6)/(X**2 + Y**2)**3) + 
     1      (E*(2*Y**2 + 4*Y**4 + X**2*(-3 + 4*Y**2)))/(X**2 + Y**2)**2.5))/
     1  E**(E - SQRT(X**2 + Y**2))**2
         ENDIF
      ELSE
         R=VAR(1)
         PHI=VAR(2)
         ENERGY=-(SIN(2*PHI)/E**(-E + R)**2)
         IF (GTEST) THEN
            V(1)=(2*(-E + R)*SIN(2*PHI))/EXP((E - R)**2)
            V(2)=(-2*COS(2*PHI))/EXP((E - R)**2)
         ENDIF
         IF (STEST) THEN
            IF (METRIC) THEN
               IF (EFFHESS) THEN
                  HESS(1,1)=(-2*(-1 + 2*E**2 + 2*R*(-2*E + R))*SIN(2*PHI))/EXP((E - R)**2)
                  HESS(1,2)=(2*(1 + 2*R*(-E + R))*COS(2*PHI))/(EXP((E - R)**2)*R)
                  HESS(2,1)=(2*(1 + 2*R*(-E + R))*COS(2*PHI))/(EXP((E - R)**2)*R**3)
                  HESS(2,2)=(2*(2 - E*R + R**2)*SIN(2*PHI))/(EXP((E - R)**2)*R**2)
C                 WRITE(*,'(A,F20.10)') 'A11=',HESS(1,1)
C                 WRITE(*,'(A,F20.10)') 'A12=',HESS(1,2)
C                 WRITE(*,'(A,F20.10)') 'A21=',HESS(2,1)
C                 WRITE(*,'(A,F20.10)') 'A21=',HESS(2,2)
               ELSE
                  HESS(1,1)=(-2*(-1 + 2*E**2 + 2*R*(-2*E + R))*SIN(2*PHI))/E**(E - R)**2
                  HESS(1,2)=(4*(-E + R)*COS(2*PHI))/E**(E - R)**2
                  HESS(2,1)=(4*(-E + R)*COS(2*PHI))/(E**(E - R)**2*R**2)
                  HESS(2,2)=(4*SIN(2*PHI))/(E**(E - R)**2*R**2)
               ENDIF
            ELSE
               IF (EFFHESS) THEN
                  HESS(1,1)=(-2*(-1 + 2*E**2 + 2*R*(-2*E + R))*SIN(2*PHI))/E**(E - R)**2
                  HESS(1,2)=(4*(-E + R)*COS(2*PHI))/E**(E - R)**2-V(2)/R
                  HESS(2,1)=HESS(1,2)
                  HESS(2,2)=(4*SIN(2*PHI))/E**(E - R)**2+R*V(1)
               ELSE
                  HESS(1,1)=(-2*(-1 + 2*E**2 + 2*R*(-2*E + R))*SIN(2*PHI))/E**(E - R)**2
                  HESS(1,2)=(4*(-E + R)*COS(2*PHI))/E**(E - R)**2
                  HESS(2,1)=HESS(1,2)
                  HESS(2,2)=(4*SIN(2*PHI))/E**(E - R)**2
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END
