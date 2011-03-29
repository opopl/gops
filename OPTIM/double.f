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
C  SUBROUTINE DOUBLE ADDS A DOUBLE-WELL POTENTIAL BETWEEN THE
C  FIRST TWO ATOMS.
C
      SUBROUTINE DOUBLE(N, X, V, EDOUBLE, GTEST, STEST, CUTOFF, H, W)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), VD(6), R2, EDOUBLE,
     1                 V(3*N), R, G(2,2), H, W, CUTOFF, AD(6,6)

      R2=(X(1)-X(4))**2+(X(2)-X(5))**2+(X(3)-X(6))**2
      R=SQRT(R2)
      R2=1.0D0/R2
      EDOUBLE=H*(1.0D0-((R-CUTOFF-W)/W)**2)**2

      IF (.NOT.GTEST) RETURN
      CALL DG(G,R,VD,X,H,W,CUTOFF)
      DO J1=1,6
         V(J1)=V(J1)+VD(J1)
      ENDDO
      
      IF (.NOT.STEST) RETURN
      CALL DS(G,R,R2,AD,X,H,W,CUTOFF)
      DO J1=1,6
         DO J2=1,6
            HESS(J2,J1)=HESS(J2,J1)+AD(J2,J1)
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE DG(G,R,V,X,H,W,RWCA)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION G(2,2), V(6), X(*), DUMMY, R, H, W, RWCA
C
C  CALCULATE THE G TENSOR.
C
      G(1,1)=0.0D0
      G(2,2)=0.0D0
      G(1,2)=4*H*(R - RWCA)*(R - RWCA - 2*W)*(R - RWCA - W)/(R*W**4)
      G(2,1)=G(1,2)
C     WRITE(*,'(A,F20.10)') 'G(1,2)=',G(1,2)

      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,2
               DUMMY=DUMMY+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
            ENDDO
            V(J3)=DUMMY
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE DS(G,R,R2,A,X,H,W,RWCA)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(2,2), R, A(6,6), X(*), DUMMY, H, W, RWCA, R2, F(2,2)

      F(1,1)=0.0D0
      F(2,2)=0.0D0
      F(2,1)=4*H*(2*R**3 - 3*R**2*(RWCA + W) + RWCA*(RWCA + W)*(RWCA + 2*W))/(R*W**4)
      F(1,2)=F(2,1)
C     WRITE(*,'(A,F20.10)') 'F(1,2)=',F(1,2)
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,2
               DUMMY=DUMMY+F(J4,J1)*R2*(X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            A(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  NEXT ARE THE TERMS WHERE X_I AND X_J ARE ON THE SAME ATOM
C  BUT ARE DIFFERENT, E.G. Y AND Z.
C
      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,2
                  DUMMY=DUMMY + F(J5,J1)*R2*(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               A(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  CASE III, DIFFERENT ATOMS, SAME CARTESIAN COORDINATE.
C
      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,2
               A(3*(J4-1)+J2,J3)=-F(J4,J1)*R2*(X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
C
C  CASE IV: DIFFERENT ATOMS AND DIFFERENT CARTESIAN COORDINATES.
C
      DO J1=1,2
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,2
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  A(J6,J3)=-F(J4,J1)*R2*(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
                  A(3*(J4-1)+J2,3*(J1-1)+J5)=A(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  SYMMETRISE HESSIAN
C
      DO J1=1,6
         DO J2=J1+1,6
            A(J1,J2)=A(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END
