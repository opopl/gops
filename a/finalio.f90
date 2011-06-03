
!> @brief Produce final quenches  
      SUBROUTINE FINALIO
      USE COMMONS

      IMPLICIT NONE

      INTEGER MYUNIT2

      MYUNIT2=25 

      CALL OPENF(MYUNIT2,">",FILE='lowest')
      
      DO J1=1,NSAVE
         WRITE(MYUNIT2,*) NATOMS
         WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)
         WRITE(MYUNIT2,30) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
30       FORMAT('LA ',3F20.10)
      ENDDO

      CLOSE(MYUNIT2)

      RETURN
      END

      
