
      DO ISTEP=NDONE+1,NSTEPS 
         ISTEP = ISTEP
         CALL FLUSH(LFH)
         MCTEMP = TEMP

23       CONTINUE

         SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS)
 
         CALL TAKESTEP
              
         NQ=NQ+1
CALL QUENCH(.FALSE.,ITERATIONS,TIME,BRUN,QDONE,SCREENC)  
NQTOT=NQTOT+1
 
WRITE(LFH,111)  'Qu ',          NQ,	&
                ' E=',          POTEL,	&
                ' steps=',      ITERATIONS,	&
                ' RMS=',        RMS,	&
                ' Markov E=',   EPREV,	&
                ' t=',TIME-TSTART	

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
                        WRITE(LFH,22) 
                        'mc> WARNING - energy for saved coordinates ',OPOTEL,&
                        ' differs from Markov energy ',
                        EPREV,' because an atom moved outside the container'
                     ELSE
                        WRITE(LFH,'(2(A,G20.10))') 'mc> ERROR - energy for coordinates in COORDSO=',OPOTEL,
     &                                                 ' but Markov energy=',EPREV
                        STOP
                     ENDIF
                  ENDIF
               ENDIF
! Accept or reject step. If the quench did not converge then allow a
! potenial move, but count it as a rejection in terms of NSUCCESS and
! NFAIL. This way we will accept a lower minimum if found, 
! but the steps won't become so big.
               IF (ATEST) THEN
                  IF (DEBUG) THEN
                     WRITE(LFH,34) JP,RANDOM,POTEL,EPREV,NSUCCESS,NFAIL
                  ENDIF
                  IF ((ISTEP-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV).GT.ECONV) THEN
                     JACCPREV=ISTEP
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
                     WRITE(LFH,36) JP,RANDOM,POTEL,EPREV,NSUCCESS,NFAIL
36                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
               ENDIF
            ENDIF
!
!  If RESTART then reseed if we haven t accepted a step in twice the relaxation time.
!
            IF (RESTART.AND.(ISTEP-JACCPREV.GT.1.1D0*NRELAX)) CALL REST(ITERATIONS,TIME,ISTEP,RCOORDS,RMIN,RVAT,JACCPREV)
!
!  Check the acceptance ratio.
! 
            IF ( (MOD(ISTEP,NACCEPT).EQ.0) &
                .AND. (NSEED.EQ.0) &
                .AND. (.NOT.STAY)) CALL ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)

            TEMP=TEMP*SCALEFAC
            IF (HIT) GOTO 37
            IF (DUMPINT.GT.0) THEN
               IF (MOD(ISTEP,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(ISTEP,EBEST,BESTCOORDS,JBEST)
               ENDIF
            ENDIF
           IF (NQ.GT.NSTEPS) GOTO 37

         ENDDO
 
         CALL FLUSH(LFH)
      ENDDO

