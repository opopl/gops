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
C  SUBROUTINE P4DIFF CALCULATES THE CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY. REDUCED UNITS.
C
C*************************************************************************
C
      SUBROUTINE P4DIFF(N, X, V, ENERGY, K, GTEST, STEST)
      IMPLICIT NONE
      INTEGER N, L, J
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(N), ENERGY, K, DUM, V(N)

C     K=0.0 ! NOW ASSIGNED VIA PARAM1 ALLOWING US TO SET IT IN ODATA
      L=SQRT(N*1.0D0)
C 
C  CALCULATE ENERGY.
C
      ENERGY=0.0D0
      DO J=1,N
         DUM=0.5*(X(J))**2
         IF (J+L.GT.N) THEN
            IF (J+1.GT.N) THEN
               ENERGY=ENERGY+DUM*(DUM-1.0)-K*X(J)*(X(J+1-N)+X(J+L-N))
            ELSE 
               ENERGY=ENERGY+DUM*(DUM-1.0)-K*X(J)*(X(J+1)+X(J+L-N))
            ENDIF
         ELSE
            ENERGY=ENERGY+DUM*(DUM-1.0)-K*X(J)*(X(J+1)+X(J+L))
         ENDIF
      ENDDO

      IF (.NOT.GTEST) RETURN
      CALL P4G(V,X,N,L,K)
      
      IF (.NOT.STEST) RETURN
      CALL P4S(X,N,L,K)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE P4G(V,X,N,L,K)
      IMPLICIT NONE
      INTEGER N, L, J
      DOUBLE PRECISION V(N), X(N), K

      DO J=1,N
         V(J)=X(J)*(X(J)*X(J)-1.0)
         IF (J+L .GT. N) THEN
            IF (J+1 .GT. N) THEN
               V(J)=V(J)-K*(X(J+1-N)+X(J+L-N)+X(J-1)+X(J-L))
            ELSE
               V(J)=V(J)-K*(X(J+1)+X(J+L-N)+X(J-1)+X(J-L))
            ENDIF
         ELSE IF (J-L .LT. 1) THEN
            IF (J-1 .LT. 1) THEN
               V(J)=V(J)-K*(X(J+1)+X(J+L)+X(J-1+N)+X(J-L+N))
            ELSE
               V(J)=V(J)-K*(X(J+1)+X(J+L)+X(J-1)+X(J-L+N))
            ENDIF
         ELSE
            V(J)=V(J)-K*(X(J+1)+X(J+L)+X(J-1)+X(J-L))
         ENDIF
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE P4S(X,N,L,K)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, L, J1, J2
      DOUBLE PRECISION X(N), K
C
C  NOW DO THE HESSIAN.
C
      HESS(1:N,1:N)=0.0D0
      DO J1=1,N
         HESS(J1,J1)=3.0*X(J1)*X(J1)-1.0
      ENDDO

      HESS(N,1)=-K
      HESS(1,N)=-K
      DO J1=1,N-1
         HESS(J1,J1+1)=-K
         HESS(J1+1,J1)=-K
      ENDDO

      DO J1=1,N-L
         HESS(J1,J1+L)=-K
         HESS(J1+L,J1)=-K
      ENDDO

      DO J1=N-L+1,N
         HESS(J1,J1-N+L)=-K
         HESS(J1-N+L,J1)=-K
      ENDDO

      RETURN
      END
