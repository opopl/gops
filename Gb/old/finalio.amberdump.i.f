
SUBROUTINE amberdump(J1,QMINP)
      USE commons
      USE modamber
      IMPLICIT NONE


      CHARACTER(LEN=25) coordfile
      CHARACTER(LEN=2) FNAME
      INTEGER J1
      DOUBLE PRECISION QMINP(NSAVE,3*NATOMS)

      IF (J1.LT.10) THEN
         WRITE (FNAME,'(I1)') J1
      ELSE
         WRITE (FNAME,'(I2)') J1
      ENDIF

      DO a=1,atoms
        x(a)=QMINP(J1,3*a-2)
        y(a)=QMINP(J1,3*a-1)
        z(a)=QMINP(J1,3*a)
      END DO

      coordfile='acoords.dump.'//FNAME

      OPEN (UNIT=4,IOSTAT=ios,FILE=coordfile,STATUS='UNKNOWN')

      DO a=1,atoms
        WRITE (UNIT=4,FMT='(A1,2X,A2,2X,I3,2X,I3,2X,F7.3,3X,F7.3,3X,F7.3)') label(a),typech(a),
     1        a,bondedto(a),x(a),y(a),z(a)
      ENDDO

      WRITE (UNIT=4,FMT='(A3)') 'end'
      WRITE (UNIT=4,FMT='(A)') ' '
      WRITE (UNIT=4,FMT='(A4,7X,I2)') 'loop',rings

      DO a=1,rings
        WRITE (UNIT=4,FMT='(I3,4X,I3)') loopatom(2*a-1),loopatom(2*a)
      END DO

      WRITE (UNIT=4,FMT='(A)') ' '
      WRITE (UNIT=4,FMT='(A7)') 'charges'

      DO a=1,atoms
        q(a)=q(a)/18.2223
        WRITE (UNIT=4,FMT='(I3,2X,F7.4)') a,q(a)
      END DO

      WRITE (UNIT=4,FMT='(A3)') 'end'

      RETURN

      END

