
      SUBROUTINE FINALIO
!op226> Declarations {{{ 
      USE COMMONS
      USE QMODULE

      IMPLICIT NONE
      ! }}}
! subroutine body {{{
      PI = 4.D0*DATAN(1.D0)

         MYUNIT2=25 
         OPEN(MYUNIT2,FILE='lowest',STATUS='UNKNOWN')

      DO J1=1,NSAVE
      ! {{{

         WRITE(MYUNIT2,*) NATOMS
         WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)
C        
C   OUTPUT CORDS FOR LOWEST IN X,Y,Z FORMAT
C   OUTPUT CORDS FOR MOVIE_GMIN IN MOVIESEG FORMAT
C
         IF (BLNT.AND.(.NOT.P46).AND.(.NOT.G46)) THEN
C
C this writes 'lowest' in xyz (Xmakemol) format
C
            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(2A1,1X,3F20.10)') BEADLETTER(J2),'L',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE

         if(ns.gt.0) then
             WRITE(MYUNIT2,50) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
50           FORMAT('LB',3F20.10)
            ENDIF
         ENDIF
         ! }}}
      ENDDO

C  End of loop over dump to file lowest or equivalent.

      CLOSE(MYUNIT2)

      ! }}}
      RETURN
      END
     
