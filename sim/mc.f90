
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

C  Main basin-hopping loop 
! {{{
      DO J1=NDONE+1,NSTEPS 
         ISTEP = J1
         CALL FLUSH(LOG_FH)
         MCTEMP = TEMP

23          CONTINUE

                  SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS)
 
                  CALL TAKESTEP
              
               NQ=NQ+1
               CALL QUENCH(.FALSE.,ITERATIONS,TIME,BRUN,QDONE,SCREENC)  
               NQTOT=NQTOT+1
 
C
C  Output
C
                        WRITE(LOG_FH,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ,' E=',
     1                 POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV,' t=',TIME-TSTART

!     mp466>  writes structure and energetic data at regular increments
!             to *plot and movie files for AMH potential
               CALL FLUSH(LOG_FH)

!     csw34> PAIRDIST prints the distance between specified atom pairs
!     after each quench. 
!
          IF (PAIRDISTT) THEN
!     Write end of previous line as using ADVANCE="NO"
                WRITE(MYPUNIT,*) " "
!     Write current quench number
                WRITE(MYPUNIT,'(I10)',ADVANCE="NO") NQ
!     For each pair, assign ATOM1 and ATOM2 arrays containing coordinates
             DO PAIRCOUNTER=1,NPAIRS
                ATOM1(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,1)-2:3*PAIRDIST(PAIRCOUNTER,1))
                ATOM2(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,2)-2:3*PAIRDIST(PAIRCOUNTER,2))
!     Call PAIRDISTANCE with (x,y,z) for each atom
                WRITE(MYPUNIT,'(F10.4)',ADVANCE="NO") PAIRDISTANCE(ATOM1,ATOM2) 
             ENDDO
             FLUSH(MYPUNIT)
          ENDIF
          
!
!     csw34> TRACKDATA keyword prints the quench energy, markov energy
!     and energy of the current lowest minimum to files for viewing during a run. 
!
          IF (TRACKDATAT) THEN
             WRITE(MYEUNIT,'(I10,F20.10)') J1,POTEL
             WRITE(MYMUNIT,'(I10,G20.10)') J1,EPREV
             WRITE(MYBUNIT,'(I10,G20.10)') J1,QMIN(1)
             CALL FLUSH(MYEUNIT)
             CALL FLUSH(MYMUNIT)
             CALL FLUSH(MYBUNIT)

             IF (RMST.AND.CHRMMT) THEN
                WRITE(4428,'(I10,F15.5)') NQ,RMSD
                CALL FLUSH(4428)
             ENDIF
          ENDIF
 
            ELSE
               ATEST=.TRUE. 

                      COLDFUSION=.FALSE.
C     csw34> Check to see if LOCALSAMPLE constraints have been violated - if they have, reject the step
               IF (LOCALSAMPLET) THEN
                  DISTOK=.FALSE.
                  CALL A9DISTCHECK(COORDS(:),DISTOK)
                  IF (.NOT.DISTOK) ATEST=.FALSE.
               ENDIF
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
C Accept or reject step. If the quench did not converge then allow a
C potenial move, but count it as a rejection in terms of NSUCCESS and
C NFAIL. This way we will accept a lower minimum if found, but the steps won;t become so big.
C However, for TIP5P some cold fusion events that had not actually reached the threshold for
C rejection were actually accepted. Must prevent this!
               IF (ATEST) THEN
                  IF (DEBUG) THEN
                     WRITE(LOG_FH,34) JP,RANDOM,POTEL,EPREV,NSUCCESS,NFAIL
34                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF
                  IF ((J1-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV).GT.ECONV) THEN
C                    NRELAX=J1-JACCPREV
C                    IF (RESTART) WRITE(LOG_FH,'(A,I6,A)') ' relaxation time set to ',NRELAX,' steps'
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
C
C  If RESTART then reseed if we haven t accepted a step in twice the relaxation time.
C
            IF (RESTART.AND.(J1-JACCPREV.GT.1.1D0*NRELAX)) CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
C
C  Check the acceptance ratio.
C 
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
! }}}

37    CONTINUE
           WRITE(LOG_FH,21) NSUCCESST*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST+NFAILT)),
     1               STEP,ASTEP,TEMP
21         FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
     
      RETURN

      END

