
      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)

      USE COMMONS
      USE PORFUNCS

      IMPLICIT NONE

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

C  Calculate the initial energy and save in EPREV
!op226>{{{ 
      WRITE(LOG_FH,'(A)') 'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
         CALL QUENCH(.FALSE.,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NQTOT=NQTOT+1
            WRITE(LOG_FH,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ,' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART

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

      WRITE(LOG_FH,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
      WRITE(LOG_FH,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

include mc.bh.inc.f90 ! Main basin-hopping loop

37    CONTINUE
           WRITE(LOG_FH,21) NSUCCESST*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST+NFAILT)),
     1               STEP,ASTEP,TEMP
21         FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
     
      RETURN

      END

      SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,NP,RANDOM,MCTEMP)

      USE COMMONS
      USE QMODULE

      IMPLICIT NONE

      DOUBLE PRECISION ENEW, EOLD, DUMMY, DPRAND, RANDOM, EREF, TEOLD, TENEW, RATIO,MCTEMP
      DOUBLE PRECISION TRANS, DISTMIN, DISTMINOLD
      LOGICAL ATEST, FLAT, evap, evapreject
      INTEGER NP,INDEXOLD, INDEXNEW, J1, NDUMMY
      DATA DISTMINOLD /0.0D0/
      COMMON /DMIN/ DISTMIN
      common /ev/ evap, evapreject

      IF (DISTMINOLD.EQ.0.0D0) DISTMINOLD=DISTMIN  ! this should allow for the first step
      
      TEOLD=EOLD
      TENEW=ENEW

C  Standard canonical sampling.
C
         IF (TENEW.LT.TEOLD) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DEXP(-(TENEW-TEOLD)/MAX(MCTEMP,1.0D-100)).GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF

      RETURN 
      END 

