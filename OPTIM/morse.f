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
C  HERE WE CALCULATE THE TWO-BODY MORSE POTENTIAL 
C                                        
C*************************************************************************
C
      SUBROUTINE MORSE (N,X,P2,VNEW,RHO,GTEST,STEST)
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4
      LOGICAL GTEST,STEST
      DOUBLE PRECISION X(3*N), P2, RHO, TEMP, DIST, 
     1                 VNEW(3*N), R, RR(N,N), ERMR(N,N), ERMRM(N,N), ERMRT(N,N)

      P2=0.0D0
      DO J1=1,N
         RR(J1,J1)=0.0D0
         ERMR(J1,J1)=0.0D0
         ERMRM(J1,J1)=0.0D0
         ERMRT(J1,J1)=0.0D0
         J3=3*(J1-1)
         DO J2=J1+1,N
            J4=3*(J2-1)
            DIST=(X(J4+1)-X(J3+1))**2+(X(J4+2)-X(J3+2))**2+(X(J4+3)-X(J3+3))**2
            DIST=DSQRT(DIST)
            TEMP=DEXP(-RHO*(DIST-1.0D0))
            P2=P2+TEMP*(TEMP-2.0)
            IF (GTEST.OR.STEST) THEN
               R=RHO*DIST
               RR(J2,J1)=1.0D0/R
               ERMR(J2,J1)=TEMP
               ERMRM(J2,J1)=TEMP-1.0D0
               ERMRT(J2,J1)=2.0D0*TEMP-1.0D0
               RR(J1,J2)=RR(J2,J1)
               ERMR(J1,J2)=ERMR(J2,J1)
               ERMRM(J1,J2)=ERMRM(J2,J1)
               ERMRT(J1,J2)=ERMRT(J2,J1)
            ENDIF
         ENDDO
      ENDDO

      IF (.NOT.GTEST) RETURN
      CALL MG(N, X, VNEW, RHO, RR, ERMR, ERMRM) 

      IF (.NOT.STEST) RETURN
      CALL MS(N, X, RHO, RR, ERMR, ERMRM, ERMRT) 

      RETURN
      END
C
C*************************************************************************
C
C  SUBROUTINE MG CALCULATES THE CARTESIAN GRADIENT ANALYTICALLY FOR THE MORSE POTENTIAL.
C
C*************************************************************************
C
      SUBROUTINE MG(N, X, V, RHO, RR, ERMR, ERMRM) 
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION X(3*N), RHO, V(3*N), RR(N,N),
     1                 DUMMY1, ERMR(N,N), ERMRM(N,N)
C
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY1=0.0D0
            DO J4=1,N
               DUMMY1=DUMMY1-2.0*RHO*RHO*(X(J3)-X(3*(J4-1)+J2))*
     1               ERMR(J4,J1)*ERMRM(J4,J1)*RR(J4,J1)
            ENDDO
            V(J3)=DUMMY1
         ENDDO
      ENDDO

      RETURN
      END
C
C*************************************************************************
C
C  SUBROUTINE MS CALCULATES THE CARTESIAN SECOND DERIVATIVE MATRIX 
C  ANALYTICALLY FOR THE MORSE POTENTIAL.
C
C*************************************************************************
C
      SUBROUTINE MS(N, X, RHO, RR, ERMR, ERMRM, ERMRT) 
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5, J6 
      DOUBLE PRECISION X(3*N), RHO, ERMR(N,N),
     1                 RR(N,N), DUMMY1, ERMRM(N,N),
     2                 ERMRT(N,N)
C
C  FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY1=0.0D0
            DO J4=1,N
               DUMMY1=DUMMY1+2.0*RHO*RHO*ERMR(J4,J1)
     1         *(  ERMRM(J4,J1) 
     2         * ( (RHO * (X(J3)-X(3*(J4-1)+J2) )*RR(J4,J1) )**2-1.0D0)  
     4         + (  RHO*(X(J3)-X(3*(J4-1)+J2)) )**2*RR(J4,J1)*
     5              ERMRT(J4,J1) ) * RR(J4,J1)
            ENDDO
            HESS(J3,J3)=DUMMY1
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
               DUMMY1=0.0D0
               DO J5=1,N
                 DUMMY1=DUMMY1
     1           +2.0*RHO**4*ERMR(J5,J1)
     2           *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
     3           *( ERMRM(J5,J1)*RR(J5,J1)+
     4           ERMRT(J5,J1)) * RR(J5,J1)**2
               ENDDO
               HESS(J3,3*(J1-1)+J4)=DUMMY1
               HESS(3*(J1-1)+J4,J3)=DUMMY1
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
               HESS(J3,3*(J4-1)+J2)=2.0*RHO**2*ERMR(J4,J1)
     1      *(  (-RHO**2*(X(J3)-X(3*(J4-1)+J2))**2*RR(J4,J1)**2+1.0) 
     2           *   ERMRM(J4,J1)
     3           -   (X(J3)-X(3*(J4-1)+J2))**2 * RR(J4,J1) * RHO**2
     4           *   ERMRT(J4,J1)   )  * RR(J4,J1)
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
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
                  HESS(J3,J6)=
     1           -2.0*RHO**4*ERMR(J4,J1)
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3           *( ERMRM(J4,J1)*RR(J4,J1)+
     4           ERMRT(J4,J1)) * RR(J4,J1)**2

                  HESS(J6,J3)=HESS(J3,J6)
               ENDDO
               DO J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=
     1           -2.0*RHO**4*ERMR(J4,J1)
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3           *( ERMRM(J4,J1)*RR(J4,J1)+
     4           ERMRT(J4,J1)) * RR(J4,J1)**2
                  HESS(J6,J3)=HESS(J3,J6)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
