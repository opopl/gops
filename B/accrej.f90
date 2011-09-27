      SUBROUTINE ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
      ! {{{
      USE COMMONS
      USE V
      IMPLICIT NONE

      ! sub 
      INTEGER,DIMENSION(:),INTENT(INOUT) ::  NSUCCESS, NFAIL, NFAILT, NSUCCESST
      INTEGER,INTENT(IN) :: JP

      ! loc 
      INTEGER J1, J2 
      LOGICAL EVAP, EVAPREJECT
      DOUBLE PRECISION P0,FAC
      COMMON /EV/ EVAP, EVAPREJECT

      P0=1.D0*NSUCCESS(JP)/(1.D0*(NSUCCESS(JP)+NFAIL(JP)))
      
      IF (P0.GT.ACCRAT(JP)) THEN
        FAC=1.05D0
      ELSE
        FAC=1.D0/1.05D0
      ENDIF
        STEP(JP)=FAC*STEP(JP)
        ASTEP(JP)=FAC*ASTEP(JP)

      WRITE(LFH,'(A,I6,A,F8.4,A,F8.4)') &
		&            'Acceptance ratio for previous ',NACCEPT,&
		&            ' steps=',P0,&
		&            '  FAC=',FAC


      IF (FIXBOTH(JP)) THEN
      ELSE IF (FIXSTEP(JP)) THEN
         IF(.NOT.FIXTEMP(JP)) WRITE(LFH,'(A,F12.4)') 'Temperature is now:',TEMP(JP)
      ELSE
         WRITE(LFH,'(A)',ADVANCE='NO') 'Steps are now:'
         WRITE(LFH,'(A,F10.4)',ADVANCE='NO') '  STEP=',STEP(JP)    
         IF(ASTEP(JP).GT.0.D0) WRITE(LFH,'(A,F10.4)',ADVANCE='NO')'  ASTEP=',ASTEP(JP) 
         IF(.NOT.FIXTEMP(JP)) WRITE(LFH,'(A,F10.4)') ' Temperature is now:',TEMP(JP)
         IF (RIGID) WRITE(LFH,'(A,F12.6,A,F12.6)') 'Maximum rigid body rotational move is now ',OSTEP(JP)
      ENDIF

      IF (FIXD) WRITE(LFH,'(A,I4)') 'hard sphere collision moves=',NHSMOVE
!
      NSUCCESST(JP)=NSUCCESST(JP)+NSUCCESS(JP)
      NFAILT(JP)=NFAILT(JP)+NFAIL(JP)
      NSUCCESS(JP)=0
      NFAIL(JP)=0 
!
      RETURN
      ! }}}
      END

