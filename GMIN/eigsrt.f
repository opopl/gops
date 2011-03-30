      SUBROUTINE EIGSRT(D,V,N,NP)
      DOUBLE PRECISION D(N), V(NP,N), P
      DO 30 I=1,N-1
         K=I
         P=D(I)
         DO 10 J=I+1,N
            IF (D(J).GE.P) THEN
               K=J
               P=D(J)
            ENDIF
10       CONTINUE
         IF (K.NE.I) THEN
            D(K)=D(I)
            D(I)=P
C
C  Interchanges columns of V
C
            DO 20 J=1,N
               P=V(J,I)
               V(J,I)=V(J,K)
               V(J,K)=P
20          CONTINUE
         ENDIF
30    CONTINUE
      RETURN
      END
