      SUBROUTINE MCRUNS(SCREENC)
      USE COMMONS
      IMPLICIT NONE
      LOGICAL LOPEN
      DOUBLE PRECISION SCREENC(3*NATOMS)

      INTEGER J1

      IF (PTMC.OR.BSPT) THEN
         RETURN
      ENDIF

      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'mcruns> ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      IF(CHECKDT) THEN
      ENDIF

      DO J1=1,NRUNS
      ENDDO


      RETURN
      END
