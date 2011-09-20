!op226>=================================== 

      SUBROUTINE FINALIO

      USE COMMONS
      USE QMODULE

      IMPLICIT NONE
      
      LOGICAL :: GTEST
      INTEGER MYUNIT2,J1,J2,J3

! subroutine body {{{
      PI = 4.D0*DATAN(1.D0)

      MYUNIT2=25 
      OPEN(MYUNIT2,FILE='lowest',STATUS='UNKNOWN')

      DO J1=1,NSAVE
         WRITE(MYUNIT2,*) NATOMS
         WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)
         IF (BLNT.AND.(.NOT.P46).AND.(.NOT.G46)) THEN
            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(2A1,1X,3F20.10)') BEADLETTER(J2),'L',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE
            WRITE(MYUNIT2,30) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
30          FORMAT('LA ',3F20.10)
         ENDIF
      CLOSE(MYUNIT2)

      ! }}}
      RETURN
      END

      !include "finalio.amberdump.i.f"
      !include "finalio.capsidio.i.f"
      !include "finalio.rbio.i.f"
      !include "finalio.tipio.i.f"
