      SUBROUTINE ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
      ! {{{
      USE COMMONS
      USE V
      IMPLICIT NONE

      ! sub 
      INTEGER,DIMENSION(:) ::  NSUCCESS, NFAIL, NFAILT, NSUCCESST
      INTEGER JP

      ! loc 
      INTEGER J1, J2, NDUMMY
      LOGICAL evap, evapreject
      DOUBLE PRECISION P0,FAC
!     COMMON /IG/ IGNOREBIN, FIXBIN
     !COMMON /MOVE/ TMOVE, OMOVE
      common /ev/ evap, evapreject

      P0=1.D0*NSUCCESS(JP)/(1.D0*(NSUCCESS(JP)+NFAIL(JP)))
      
      IF (P0.GT.ACCRAT(JP)) THEN
        ! {{{
         IF(ARMT) THEN
           FAC=LOG(ARMA*ACCRAT(JP)+ARMB)/LOG(ARMA*P0+ARMB)
         ELSE
           FAC=1.05D0
         ENDIF
         IF (FIXBOTH(JP)) THEN
         ELSE IF (FIXSTEP(JP)) THEN
            IF (.NOT.FIXTEMP(JP)) TEMP(JP)=TEMP(JP)/1.05D0
         ELSE
            IF (FIXD) THEN
               NHSMOVE=NHSMOVE+1 
            ELSE
               IF (RIGID) THEN
                  IF (TMOVE(JP)) STEP(JP)=STEP(JP)*1.05D0 
                  IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)*1.05D0
               ELSE
                  STEP(JP)=FAC*STEP(JP)
               ENDIF
            ENDIF
            ASTEP(JP)=ASTEP(JP)*1.05D0
         ENDIF
         ! }}}
      ELSE
        ! {{{
         IF(ARMT) THEN
           FAC=LOG(ARMA*ACCRAT(JP)+ARMB)/LOG(ARMA*P0+ARMB)
         ELSE
           FAC=1.D0/1.05D0
         ENDIF
         IF (FIXBOTH(JP)) THEN
         ELSE IF (FIXSTEP(JP)) THEN
            IF (.NOT.FIXTEMP(JP)) TEMP(JP)=TEMP(JP)*1.05D0
         ELSE
            STEP(JP)=FAC*STEP(JP)
            !IF (FIXD) THEN
               !NHSMOVE=MAX(1,NHSMOVE-1)
            !ELSE
               !IF (RIGID) THEN
                  !IF (TMOVE(JP)) STEP(JP)=STEP(JP)/1.05D0
                  !IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)/1.05D0
               !ELSE
                  !STEP(JP)=FAC*STEP(JP)
               !ENDIF
            !ENDIF
            ASTEP(JP)=ASTEP(JP)/1.05D0
         ENDIF
         ! }}}
      ENDIF
!
! Prevent steps from growing out of bounds. The value of 1000 seems sensible, until
! we do something with such huge dimensions?!
!
      STEP(JP)=MIN(STEP(JP),1.0D3)
      OSTEP(JP)=MIN(OSTEP(JP),1.0D3)
      ASTEP(JP)=MIN(ASTEP(JP),1.0D3)
!
      IF (NPAR.GT.1) THEN
         WRITE(LFH,'(A,I2,A,I6,A,F8.4,A,F8.4)') '[',JP,']Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
      ELSE
         WRITE(LFH,'(A,I6,A,F8.4,A,F8.4)') 'Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
      ENDIF
      IF (FIXBOTH(JP)) THEN
      ELSE IF (FIXSTEP(JP)) THEN
         IF(.NOT.FIXTEMP(JP)) WRITE(LFH,'(A,F12.4)') 'Temperature is now:',TEMP(JP)
      ELSE
         IF (NPAR.GT.1) THEN
            WRITE(LFH,'(A,I2,A)',ADVANCE='NO') '[',JP,']Steps are now:'
         ELSE
            WRITE(LFH,'(A)',ADVANCE='NO') 'Steps are now:'
         ENDIF
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

