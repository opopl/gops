
      SUBROUTINE FINALIO
!op226> Declarations {{{ 
      USE COMMONS
      USE V

      IMPLICIT NONE

      INTEGER J1, J2, J3
      ! }}}
! subroutine body {{{

      LE_FH=25 
      OPEN(LE_FH,FILE=LE_FILE,STATUS='UNKNOWN')
      OPEN(EA_FH,FILE=EA_FILE,STATUS='UNKNOWN',ACCESS='APPEND')

      WRITE(EA_FH,'(A1,A)') "#",CMDLINE
      WRITE(EA_FH,'()') PFORCE,EAMIN(1,1:6)
      DO J1=1,NSAVE
      ! {{{

         WRITE(LE_FH,'(I)') NATOMS
         WRITE(LE_FH,10) J1,QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)
            DO J2=1,NATOMS
               WRITE(LE_FH,'(1A1,1X,3F20.10)') BEADLETTER(J2),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
       ENDDO

      CLOSE(LE_FH)
      CLOSE(EA_FH)

      ! }}}
      RETURN
      END
     
