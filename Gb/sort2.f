
      SUBROUTINE GSORT2(N,NATOMS)

      USE QMODULE

      IMPLICIT NONE

      INTEGER NATOMS
      INTEGER J1, L, N, J3, J2, NTEMP
      DOUBLE PRECISION TEMP, C
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (QMIN(L).GT.QMIN(J2)) L=J2
10       CONTINUE
         TEMP=QMIN(L)
         QMIN(L)=QMIN(J1)
         QMIN(J1)=TEMP
         NTEMP=FF(L)
         FF(L)=FF(J1)
         FF(J1)=NTEMP
         DO J2=1,3*NATOMS
            C=QMINP(L,J2)
            QMINP(L,J2)=QMINP(J1,J2)
            QMINP(J1,J2)=C
         ENDDO
20    CONTINUE
      RETURN
      END
