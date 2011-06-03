
!> 
!> @name FINALQ
!> @brief Make sure the lowest minima are tightly converged and then sort them just to be on the safe side.

      SUBROUTINE FINALQ

      USE COMMONS

      IMPLICIT NONE

      INTEGER J1, J2, ITERATIONS, BRUN,QDONE,J3
      DOUBLE PRECISION POTEL, TIME

      COMMON /MYPOT/ POTEL
 
      NQ=0
      MAXIT=MAXIT2

      DO J1=1,NSAVE
         IF (QMIN(J1).LT.1.0D10) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,1)=QMINP(J1,J2)
            ENDDO
            NQ=NQ+1
            CALL QUENCH(.TRUE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            WRITE(MYUNIT,'(A,I6,A,F20.10,A,I5,A,F15.7,A,F12.2)') 'Final Quench ',NQ,' energy=',
     1                POTEL,' steps=',ITERATIONS,' RMS force=',RMS,' time=',TIME-TSTART

            QMIN(J1)=POTEL
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=COORDS(J2,1)
            ENDDO
         ENDIF
      ENDDO

      NSAVE=NQ

      CALL GSORT(NSAVE,NATOMS)

      RETURN
      END
