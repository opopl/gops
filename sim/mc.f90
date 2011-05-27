
      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)

      USE COMMONS
      USE PORFUNCS

      IMPLICIT NONE
! inits {{{

      INTEGER NFAIL, NFAILT, NSUCCESS, NSUCCESST, NQTOT

      NRMS=0
      STAY=.FALSE.
      JACCPREV=0
      NQTOT=0
      RMINO=1.0D100
      RMIN=1.0D100
      NREN=NRENORM

      NSUCCESS=0
      NFAIL=0
      NSUCCESST=0
      NFAILT=0
! }}}

include fmt.inc.f90

!  Calculate the initial energy and save in EPREV
!op226>{{{ 
WRITE(LFH,'(A)') 'Calculating initial energy'
EPSSAVE=EPSSPHERE
EPSSPHERE=0.0D0
CALL QUENCH(.FALSE.,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
NQTOT=NQTOT+1
WRITE(LFH,111)  'Qu ',          NQ,                        	&
                ' E=',          POTEL,	&
                ' steps=',      ITERATIONS,	&
                ' RMS=',RMS,	&
                ' Markov E=',POTEL,	&
                ' t=',TIME-TSTART

EPREV=POTEL
EPPREV=0.0D0
ELASTSYM=0.0D0
EBEST=POTEL
BESTCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS)
JBEST=0
RMIN=POTEL
RCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,1)
COORDSO(1:3*NATOMS)=COORDS(1:3*NATOMS)
VATO(1:NATOMS)=VAT(1:NATOMS)
EPSSPHERE=EPSSAVE
!op226>}}} 

WRITE(LFH,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
WRITE(LFH,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

include mc.bh.inc.f90 ! Main basin-hopping loop

37    CONTINUE

ARATIO=NSUCCESST*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST+NFAILT))

WRITE(LFH,111)  'Qu ',          NQ,	&
                ' E=',          POTEL,	&
                ' steps=',      ITERATIONS,	&
                ' RMS=',        RMS,	&
                ' Markov E=',   POTEL,	&
                ' t=',          TIME-TSTART	

WRITE(LFH,21)   ' Acceptance ratio for run=',  ARATIO,	&
                ' Step=',                      STEP,     	&
                ' Angular step factor=',       ASTEP,	&
                ' T=',                         TEMP	&
RETURN

END

      SUBROUTINE ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)

      USE COMMONS
      IMPLICIT NONE

      INTEGER NSUCCESS, NFAIL, JP, NFAILT, NSUCCESST, J1, J2, NDUMMY
      LOGICAL evap, evapreject
      DOUBLE PRECISION DUMMY, DUMMY2, DUMMY3, DUMMY4, HWMAX,P0,FAC
      COMMON /EV/ EVAP, EVAPREJECT

      P0=1.D0*NSUCCESS/(1.D0*(NSUCCESS+NFAIL))
      
      IF (P0.GT.ACCRAT) THEN
           FAC=1.05D0
         IF (FIXBOTH) THEN
         ELSE IF (FIXSTEP) THEN
            IF (.NOT.FIXTEMP) TEMP(JP)=TEMP(JP)/1.05D0
         ELSE
            IF (FIXD) THEN
               NHSMOVE=NHSMOVE+1 
            ELSE
               IF (RIGID) THEN
                  IF (TMOVE) STEP=STEP*1.05D0
                  IF (OMOVE) OSTEP=OSTEP*1.05D0
               ELSE
                  STEP=FAC*STEP
                  IF(CHRIGIDTRANST.AND.CHRMMT) TRANSMAX=FAC*TRANSMAX
                  IF(CHRIGIDROTT.AND.CHRMMT) ROTMAX=FAC*ROTMAX  
               ENDIF
            ENDIF
            ASTEP=ASTEP*1.05D0
         ENDIF
      ELSE
           FAC=LOG(ARMA*ACCRAT+ARMB)/LOG(ARMA*P0+ARMB)
         IF (FIXBOTH) THEN
         ELSE IF (FIXSTEP) THEN
            IF (.NOT.FIXTEMP) TEMP(JP)=TEMP(JP)*1.05D0
         ELSE
            IF (FIXD) THEN
               NHSMOVE=MAX(1,NHSMOVE-1)
            ELSE
               IF (RIGID) THEN
                  IF (TMOVE) STEP=STEP/1.05D0
                  IF (OMOVE) OSTEP=OSTEP/1.05D0
               ELSE
                  STEP=FAC*STEP
               ENDIF
            ENDIF
            ASTEP=ASTEP/1.05D0
         ENDIF
      ENDIF
!
! Prevent steps from growing out of bounds. The value of 1000 seems sensible, until
! we do something with such huge dimensions?!
c
      STEP=MIN(STEP,1.0D3)
      OSTEP=MIN(OSTEP,1.0D3)
      ASTEP=MIN(ASTEP,1.0D3)
!
         WRITE(MYUNIT,'(A,I6,A,F8.4,A,F8.4)') 'Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
      IF (FIXBOTH) THEN
      ELSE IF (FIXSTEP) THEN
         IF(.NOT.FIXTEMP) WRITE(MYUNIT,'(A,F12.4)') 'Temperature is now:',TEMP(JP)
      ELSE
            WRITE(MYUNIT,'(A)',ADVANCE='NO') 'Steps are now:'
         WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO') '  STEP=',STEP    
         IF(ASTEP.GT.0.D0) WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO')'  ASTEP=',ASTEP 
         IF(.NOT.FIXTEMP) WRITE(MYUNIT,'(A,F10.4)') ' Temperature is now:',TEMP(JP)
         IF (RIGID) WRITE(MYUNIT,'(A,F12.6,A,F12.6)') 'Maximum rigid body rotational move is now ',OSTEP
      ENDIF
      IF (FIXD) WRITE(MYUNIT,'(A,I4)') 'hard sphere collision moves=',NHSMOVE
!
      NSUCCESST=NSUCCESST+NSUCCESS
      NFAILT=NFAILT+NFAIL
      NSUCCESS=0
      NFAIL=0 
!
      RETURN
      END

      SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,RANDOM,MCTEMP)

      IMPLICIT NONE

      DOUBLE PRECISION ENEW, EOLD, DPRAND, RANDOM, MCTEMP
      LOGICAL ATEST

         IF (ENEW.LT.EOLD) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DEXP(-(ENEW-EOLD)/MAX(MCTEMP,1.0D-100)).GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF

      RETURN 
      END 

