
      SUBROUTINE MATMULV(A,B,C,NA,NB,NC)
      IMPLICIT NONE
      INTEGER NA, NB, NC, I, K, J
      DOUBLE PRECISION A(NC,NA),B(NB,NA),C(NB,NC),Z

      DO I=1,NA
         DO K=1,NC
            Z=0.D0
            DO J=1,NB
               Z=Z+B(J,I)*C(J,K)
            ENDDO
            A(K,I)=Z
         ENDDO
      ENDDO
      RETURN 
      END

