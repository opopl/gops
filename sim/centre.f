
      SUBROUTINE SETCENTRE(X)

      USE COMMONS

      IMPLICIT NONE
      DOUBLE PRECISION RMASS(3), R(NATOMS,3)
      INTEGER I, J1

C csw34> XMASS, YMASS and ZMASS are the components of the COM position
C vector
      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      IF (RIGID) THEN
C csw34> First need to calculate the position of the COM
         DO I=1,NATOMS/2
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         XMASS=2*XMASS/NATOMS
         YMASS=2*YMASS/NATOMS
         ZMASS=2*ZMASS/NATOMS
C csw34> Then need to move the COM to the new centre via the origin
         DO I=1,NATOMS/2
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS+CENTX
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS+CENTY
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS+CENTZ
         ENDDO
      ELSE
         DO I=1,NATOMS
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         XMASS=XMASS/NATOMS
         YMASS=YMASS/NATOMS
         ZMASS=ZMASS/NATOMS
         DO I=1,NATOMS
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS+CENTX
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS+CENTY
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS+CENTZ
         ENDDO
      ENDIF
      IF (DEBUG) WRITE(MYUNIT,'(A,3F12.4)') 'centre of mass moved to ',CENTX,CENTY,CENTZ
      RETURN
      END

