
      SUBROUTINE FINALIO
!op226> Declarations {{{ 
      USE COMMONS
      USE V

      IMPLICIT NONE

      INTEGER MYUNIT2, J1, J2, J3
      ! }}}
! subroutine body {{{

      MYUNIT2=25 
      OPEN(MYUNIT2,FILE=LE_FILE,STATUS='UNKNOWN')

      DO J1=1,NSAVE
      ! {{{

         WRITE(MYUNIT2,*) NATOMS
         WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)
!
         IF (BLNT.AND.(.NOT.P46).AND.(.NOT.G46)) THEN
!
! this writes 'lowest' in xyz (Xmakemol) format
!
            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(2A1,1X,3F20.10)') BEADLETTER(J2),'L',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE
             !WRITE(MYUNIT2,50) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
             !WRITE(MYUNIT2,50) (QMINP(J1,1:NR),J2=1,NR)
             WRITE(MYUNIT2,50) QMINP(J1,1:NR)
50           FORMAT('LB',3F20.10)
         ENDIF

         !IF(NS.GT.0) THEN
             !WRITE(MYUNIT2,50) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
!50           FORMAT('LB',3F20.10)
            !ENDIF
         !ENDIF
         ! }}}
      ENDDO

!  End of loop over dump to file lowest or equivalent.

      CLOSE(MYUNIT2)

      ! }}}
      RETURN
      END
     
