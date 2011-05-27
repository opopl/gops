
      DO J1=NDONE+1,NSTEPS 
         ISTEP = J1
         CALL FLUSH(LOG_FH)
         MCTEMP = TEMP

23       CONTINUE

         SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS)
 
         CALL TAKESTEP
              
         NQ=NQ+1
CALL QUENCH(.FALSE.,ITERATIONS,TIME,BRUN,QDONE,SCREENC)  
NQTOT=NQTOT+1
 
WRITE(LOG_FH,111) 'Qu ',NQ,' E=',POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV,' t=',TIME-TSTART
111 format(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)

!     mp466>  writes structure and energetic data at regular increments
!             to *plot and movie files for AMH potential
               CALL FLUSH(LOG_FH)

!
include mc.bh.pairdist.inc.f90
include mc.bh.trackdata.inc.f90

            ELSE
               ATEST=.TRUE. 

                      COLDFUSION=.FALSE.
               IF ((QDONE.EQ.0).AND.TIP) THEN
                  ATEST=.FALSE.
               ELSEIF (ATEST) THEN
                  CALL TRANSITION(POTEL,EPREV,ATEST,RANDOM,MCTEMP)
               ENDIF
!
!  Sanity check to make sure the Markov energy agrees with COORDSO. 
!  Stop if not true.
!
               IF (DEBUG.OR.CHECKMARKOVT) THEN
                  CALL POTENTIAL(COORDSO(:),GRAD,OPOTEL,.FALSE.,.FALSE.)
                  IF (ABS(OPOTEL-EPREV).GT.ECONV) THEN
                     IF (EVAP) THEN
                        WRITE(LOG_FH,'(3(A,G20.10))') 'mc> WARNING - energy for saved coordinates ',OPOTEL,
     &                     ' differs from Markov energy ',EPREV,' because an atom moved outside the container'
                     ELSE
                        WRITE(LOG_FH,'(2(A,G20.10))') 'mc> ERROR - energy for coordinates in COORDSO=',OPOTEL,
     &                                                 ' but Markov energy=',EPREV
                        STOP
                     ENDIF
                  ENDIF
               ENDIF
! Accept or reject step. If the quench did not converge then allow a
! potenial move, but count it as a rejection in terms of NSUCCESS and
! NFAIL. This way we will accept a lower minimum if found, but the steps won;t become so big.
! However, for TIP5P some cold fusion events that had not actually reached the threshold for
! rejection were actually accepted. Must prevent this!
               IF (ATEST) THEN
                  IF (DEBUG) THEN
                     WRITE(LOG_FH,34) JP,RANDOM,POTEL,EPREV,NSUCCESS,NFAIL
34                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF
                  IF ((J1-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV).GT.ECONV) THEN
!                    NRELAX=J1-JACCPREV
!                    IF (RESTART) WRITE(LOG_FH,'(A,I6,A)') ' relaxation time set to ',NRELAX,' steps'
                     JACCPREV=J1
                  ENDIF
                  IF (QDONE.EQ.1) THEN
                     NSUCCESS=NSUCCESS+1
                  ELSE
                     NFAIL=NFAIL+1
                  ENDIF
                  EPPREV=EPREV
                  EPREV=POTEL
                  COORDSO(1:3*NATOMS)=COORDS(1:3*NATOMS)
                  VATO(1:NATOMS)=VAT(1:NATOMS)
               ELSE
                  NFAIL=NFAIL+1
                  CALL MYRESET(JP,NATOMS,NPAR,NSEED)
                  IF (DEBUG) THEN
                     WRITE(LOG_FH,36) JP,RANDOM,POTEL,EPREV,NSUCCESS,NFAIL
36                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
               ENDIF
            ENDIF
!
!  If RESTART then reseed if we haven t accepted a step in twice the relaxation time.
!
            IF (RESTART.AND.(J1-JACCPREV.GT.1.1D0*NRELAX)) CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
!
!  Check the acceptance ratio.
! 
            IF ((MOD(J1,NACCEPT).EQ.0).AND.(NSEED.EQ.0).AND.(.NOT.STAY)) CALL ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)

            TEMP=TEMP*SCALEFAC
            IF (HIT) GOTO 37
            IF (DUMPINT.GT.0) THEN
               IF (MOD(J1,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(J1,EBEST,BESTCOORDS,JBEST)
               ENDIF
            ENDIF
           IF (NQ.GT.NSTEPS) GOTO 37

         ENDDO
 
         CALL FLUSH(LOG_FH)
      ENDDO

