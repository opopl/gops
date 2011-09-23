      SUBROUTINE CENTRECOM(X)
      USE commons
      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION XMASS, YMASS, ZMASS, X(3*NATOMS), TOTMASS

      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      TOTMASS=0.0D0

      DO I=1,NATOMS
         XMASS=XMASS+MASSES(I)*X(3*(I-1)+1)
         YMASS=YMASS+MASSES(I)*X(3*(I-1)+2)
         ZMASS=ZMASS+MASSES(I)*X(3*(I-1)+3)
         TOTMASS = TOTMASS + MASSES(I)
      ENDDO

      XMASS=XMASS/TOTMASS
      YMASS=YMASS/TOTMASS
      ZMASS=ZMASS/TOTMASS

      DO I=1,NATOMS
         X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
         X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
         X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
      ENDDO

      IF (DEBUG) WRITE(LFH,'(A,3F15.10)') 'centre of mass reset to the origin from ',XMASS,YMASS,ZMASS
C     PRINT*,'final coordinates in centrecom:'
C     WRITE(*,'(I5,3F15.5)') (I,X(3*(I-1)+1),X(3*(I-1)+2),X(3*(I-1)+3),I=1,NATOMS)

      RETURN
      END
