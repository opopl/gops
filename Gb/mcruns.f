
      SUBROUTINE MCRUNS(SCREENC)

      USE COMMONS
      USE MCFUNC

      IMPLICIT NONE

      DOUBLE PRECISION SCREENC(3*NATOMS)
      LOGICAL LOPEN

      INTEGER J1

      INQUIRE(UNIT=1,OPENED=LOPEN)

      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'mcruns> ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      DO J1=1,NRUNS
         CALL MC(MCSTEPS(J1),TFAC(J1),SCREENC)
      ENDDO

      CALL FINALQ
      CALL FINALIO

      RETURN
      END
