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
C  SUBROUTINE LJMS (LJ-MS POTENTIALS) CALCULATES THE CARTESIAN 
C  GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY. REDUCED UNITS, EPSILON AND RM.
C
C*************************************************************************
C
      SUBROUTINE LJMS (N, EPSILON, RM, GAMMA, X, V, ENERGY, GTEST, STEST)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, 
     1                 V(3*N), R2(N,N), R2T,
     2                 R8(N,N), G(N,N), ENERGY1,
     3                 ENERGY2, R14(N,N), F(N,N), 
     4                 GAMMA, R6(N,N), R4T, R6T,
     5                 EPSILON,RM,RM4,RM6,RM12
C 
C  STORE DISTANCE MATRICES.
C
      ENERGY1=0.0D0
      ENERGY2=0.0D0
C     EPSILON=12.D0
C     RM=0.78D0
      RM12=RM**12
      RM6=RM**6
      RM4=RM**4

      IF (GTEST) THEN
         DO J1=1,N-1
            R2(J1,J1)=0.0D0
            R6(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N-1
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6(J2,J1)=R2(J2,J1)**3
               ENERGY1=ENERGY1+R6(J2,J1)*(R6(J2,J1)-2.0D0)
               R8(J2,J1)=R2(J2,J1)**4
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
            ENDDO
         ENDDO
         J1=N
            R2(J1,J1)=0.0D0
            R6(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=1,N-1
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R2(J2,J1)=R2(J2,J1)
               R4T=R2(J2,J1)**2
               R6(J2,J1)=R2(J2,J1)**3
               ENERGY2=ENERGY2+(1.D0+GAMMA)*RM12*R6(J2,J1)*R6(J2,J1)
     1                -4.D0*GAMMA*RM6*R6(J2,J1)
     2                -3.D0*(1.D0-GAMMA)*RM4*R4T
               R8(J2,J1)=R2(J2,J1)**4
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
            ENDDO
      ELSE
         DO J1=1,N-1
            J3=3*(J1-1)
            DO J2=J1+1,N-1
               J4=3*(J2-1)
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               R2T=1.0D0/R2T
               R6T=R2T**3
               ENERGY1=ENERGY1+R6T*(R6T-2.0D0)
            ENDDO
         ENDDO
         J1=N
            J3=3*(J1-1)
            DO J2=1,N-1
               J4=3*(J2-1)
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               R2T=1.0D0/R2T
               R2T=R2T
               R4T=R2T**2
               R6T=R2T**3
               ENERGY2=ENERGY2+(1.D0+GAMMA)*RM12*R6T*R6T
     1                -4.D0*GAMMA*RM6*R6T
     2                -3.D0*(1.D0-GAMMA)*RM4*R4T
            ENDDO

      ENDIF
      ENERGY=ENERGY1+EPSILON*0.5D0*ENERGY2

      IF (.NOT.GTEST) RETURN
      CALL LJMSG(G,R14,R8,R6,V,X,N,GAMMA,RM12,RM6,RM4,EPSILON)
      
      IF (.NOT.STEST) RETURN
      CALL LJMSS(G,F,R2,R14,R8,R6,X,N,GAMMA,RM12,RM6,RM4,EPSILON)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJMSG(G,R14,R8,R6,V,X,N,GAMMA,RM12,RM6,RM4,EPSILON)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 V(3*N), X(3*N), DUMMY, R6(N,N), 
     2                 EPSILON,GAMMA,RM12,RM6,RM4
C
C  CALCULATE THE G TENSOR.
C
      DO J1=1,N-1
         G(J1,J1)=0.0D0
         DO J2=J1+1,N-1
            G(J2,J1)=-12.0D0*(R14(J2,J1)-R8(J2,J1))
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO

      J1=N
         G(J1,J1)=0.0D0
         DO J2=1,N-1
            G(J2,J1)=-6.D0*((1.D0+GAMMA)*RM12*R14(J2,J1)
     1               -2.D0*GAMMA*RM6*R8(J2,J1)
     2               -(1.D0-GAMMA)*RM4*R6(J2,J1))       
            G(J2,J1)=EPSILON*G(J2,J1)
            G(J1,J2)=G(J2,J1)
         ENDDO
C
C  FROM HERE ON DOWN THE CODE IS SYSTEM-INDEPENDENT!
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
            ENDDO
            V(J3)=DUMMY
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJMSS(G,F,R2,R14,R8,R6,X,N,GAMMA,RM12,RM6,RM4,EPSILON)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 R2(N,N), F(N,N), X(3*N),
     3                 DUMMY,GAMMA,RM12,RM6,RM4,EPSILON,R6(N,N)

      DO J1=1,N-1
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=168.0D0*R14(J2,J1)-96.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO

      J1=N
         F(J1,J1)=0.0D0
         DO J2=1,N-1
           F(J2,J1)=84.0D0*(1.D0+GAMMA)*RM12*R14(J2,J1)-16.0D0*GAMMA*RM6*
     1              R8(J2,J1)-36.D0*(1.D0-GAMMA)*RM4*R6(J2,J1)
           F(J2,J1)=EPSILON*F(J2,J1)
           F(J1,J2)=F(J2,J1)
         END DO

C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=DUMMY
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
               DUMMY=0.0D0
               DO J5=1,N
                  DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)* 
     1           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
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
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
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
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  SYMMETRISE HESSIAN
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END
