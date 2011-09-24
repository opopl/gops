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
        ! com {{{
 !!        IF(ARMT) THEN
           !FAC=LOG(ARMA*ACCRAT(JP)+ARMB)/LOG(ARMA*P0+ARMB)
         !ELSE
           !FAC=1.05D0
         !ENDIF
!         IF (FIXBOTH(JP)) THEN
         !ELSE IF (FIXSTEP(JP)) THEN
            !IF (.NOT.FIXTEMP(JP)) TEMP(JP)=TEMP(JP)/1.05D0
         !ELSE
            !IF (FIXD) THEN
               !NHSMOVE=NHSMOVE+1 
            !ELSE
               !IF (RIGID) THEN
                  !IF (TMOVE(JP)) STEP(JP)=STEP(JP)*1.05D0 
                  !IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)*1.05D0
               !ELSE
                  !STEP(JP)=FAC*STEP(JP)
               !ENDIF
            !ENDIF
            !ASTEP(JP)=ASTEP(JP)*1.05D0
         !ENDIF
         ! }}}
      ELSE
           FAC=1.D0/1.05D0
           STEP(JP)=FAC*STEP(JP)
        ! com {{{
!         IF(ARMT) THEN
           !FAC=LOG(ARMA*ACCRAT(JP)+ARMB)/LOG(ARMA*P0+ARMB)
         !ELSE
           !FAC=1.D0/1.05D0
         !ENDIF
!         IF (FIXBOTH(JP)) THEN
         !ELSE IF (FIXSTEP(JP)) THEN
            !IF (.NOT.FIXTEMP(JP)) TEMP(JP)=TEMP(JP)*1.05D0
         !ELSE
            !STEP(JP)=FAC*STEP(JP)
            !!IF (FIXD) THEN
               !!NHSMOVE=MAX(1,NHSMOVE-1)
            !!ELSE
               !!IF (RIGID) THEN
                  !!IF (TMOVE(JP)) STEP(JP)=STEP(JP)/1.05D0
                  !!IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)/1.05D0
               !!ELSE
                  !!STEP(JP)=FAC*STEP(JP)
               !!ENDIF
            !!ENDIF
            !ASTEP(JP)=ASTEP(JP)/1.05D0
         !ENDIF
         ! }}}
      ENDIF

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

