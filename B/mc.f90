
! Doxygen {{{
!> \name MC
!> \brief This subroutine ...
!> \param NSTEPS
!> \param SCALEFAC
!> \param SCREENC
! }}}
!op226>=================================== 
      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)
!op226> Declarations {{{ 
      ! modules {{{
      USE COMMONS
      USE F
      USE V
      USE PORFUNCS

      ! }}}

      IMPLICIT NONE

      ! subroutine  {{{
      INTEGER ::    NSTEPS
      DOUBLE PRECISION ::   SCALEFAC
      DOUBLE PRECISION, DIMENSION(:) ::   SCREENC
      ! }}}
      ! local  {{{
      INTEGER J1, J2, JP, J5
      DOUBLE PRECISION :: MCTEMP, OPOTEL
      LOGICAL EVAP, ATEST, EVAPREJECT, STAY
      INTEGER ITERATIONS,NQTOT,JACCPREV,BRUN,JBEST(NPAR)
      INTEGER NDONE,QDONE
      DOUBLE PRECISION ::   POTEL,RANDOM,EPSSAVE
      DOUBLE PRECISION,DIMENSION(3*NATOMS) ::   SAVECOORDS,TEMPCOORDS,RCOORDS
      DOUBLE PRECISION,DIMENSION(NPAR) ::  EPPREV,EBEST
      INTEGER,DIMENSION(NPAR) :: NSUCCESS,NFAIL,NFAILT,NSUCCESST
      DOUBLE PRECISION,DIMENSION(NATOMS) ::   RVAT,RVATO
      DOUBLE PRECISION,DIMENSION(3*NATOMS,NPAR) ::   BESTCOORDS
      DOUBLE PRECISION,DIMENSION(3*NATOMS) :: RCOORDSO, GRAD
      DOUBLE PRECISION ::   TIME
      DOUBLE PRECISION ::   ELASTSYM(NPAR)
      DOUBLE PRECISION ::   RMIN,RMINO
      INTEGER NRMS,NLAST,NREN,NSUPERCOUNT
      DOUBLE PRECISION ::   RANNJ
       ! common {{{
      COMMON /EV/ EVAP, EVAPREJECT
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT
      ! }}}

      !INTEGER UNT
      !integer ITERATIONS, NSUPERCOUNT, NQTOT, JACCPREV, NREN, NLAST, NSTEPREN, BRUN,QDONE,JBEST(NPAR)
      !integer NRMS, NDONE, I, RNDSEED, J, NTOT, IMESG, ITRAJ, ITRAJO, NEACCEPT
      !integer J3, J4, ISTAT, LOCALCOUNT
      !INTEGER :: NSYMCALL=0
      !DOUBLE PRECISION POTEL, RANDOM, DPRAND, SAVECOORDS(3*NATOMS), TEMPCOORDS(3*NATOMS)
      !DOUBLE PRECISION :: TIME, SPOTEL(NSUPER), SCOORDS(3*NATOMS,NSUPER)
      !DOUBLE PRECISION :: EPPREV(NPAR), QSTART, QFINISH, RANNJ, RMIN
      !DOUBLE PRECISION :: RMINO, RCOORDS(3*NATOMS),ELASTSYM(NPAR)
      !DOUBLE PRECISION ::  RCOORDSO(3*NATOMS), RVAT(NATOMS), RVATO(NATOMS), EPSSAVE, EBEST(NPAR)
      !DOUBLE PRECISION ::  BESTCOORDS(3*NATOMS,NPAR), endtime, RMSD, VINIT, CTE, TEMPTRAJ(0:NPAR-1)
      !DOUBLE PRECISION ::  T, BETA(0:NPAR-1), GRAD(3*NATOMS), E, ER, W, DELTA, DBETA, A9ANGLE 
      !DOUBLE PRECISION ::  DUMMY1, DUMMY2, DUMMY3, INTE, OPOTEL, DUMGRAD(3*NATOMS), DJWPOTEL
      !LOGICAL CHANGEDE, EXCHANGEACCEPT, EXCHANGE, FLAG 
      !LOGICAL CHIRALFAIL,AMIDEFAIL, LOGDUMMY, DISTOK, ATOMINGROUP(NATOMS)
      !CHARACTER FNAME*9
      !CHARACTER (LEN= 3)  ISTR
      !CHARACTER (LEN=20) QUENCHNUM, QUNAME, DUMMYCHAR
      !CHARACTER (LEN=20) BESTNAME, CURRENTBESTNAME
!!  AMH 
      !INTEGER :: gly_count,iii,i2,i500,snapcount, DUMMYINT
      !DOUBLE PRECISION prcord(NATOMS,3,3,3)
!!  csw34> PAIRDIST variables
      !INTEGER :: PAIRCOUNTER
      !DOUBLE PRECISION :: ATOM1(3),ATOM2(3)

      !LOGICAL EVAP, ATEST, STAY, evapreject, LOPEN

      !character(len=10)       :: datechar,timechar,zonechar
      !integer                 :: values(8),itime1
      !double precision :: DISTGROUPX2,DISTGROUPY2,DISTGROUPZ2,DISTGROUPCENTRE,TESTANGLE
      !integer :: J6
      ! }}}
!op226>}}} 
!op226> Subroutine body {{{ 

      EVAPREJECT=.FALSE.

      NDONE=0
      NQ(:)=NDONE

      IF (NACCEPT.EQ.0) NACCEPT=NSTEPS+1
      JACCPREV=0
      NQTOT=0
      RMINO=1.0D100
      RMIN=1.0D100
      DO JP=1,NPAR 
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         NSUCCESST(JP)=0
         NFAILT(JP)=0
      ENDDO

!  Calculate the initial energy and save in EPREV
      ! {{{
      WRITE(LFH,'(A)') 'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      DO JP=1,NPAR
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NQTOT=NQTOT+1
         WRITE(LFH,101) &
             &   'Qu ',NQ(JP),	&
             &   ' E=',POTEL,	&
             &   ' steps=',ITERATIONS,	&
             & ' RMS=',RMS,	&
             & ' Markov E=',POTEL,	&
             & ' t=',TIME-TSTART	

101      format(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)

!  EPREV saves the previous energy in the Markov chain.
!  EBEST and JBEST record the lowest energy since the last reseeding and the
!  step it was attained at. BESTCOORDS contains the corresponding coordinates.
 
         EPREV(JP)=POTEL
         EPPREV(JP)=0.0D0
         ELASTSYM(JP)=0.0D0
         EBEST(JP)=POTEL
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=0
         RMIN=POTEL
         RCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,1)
         COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
      ENDDO
      EPSSPHERE=EPSSAVE
      ! }}}

      WRITE(LFH,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
      WRITE(LFH,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

!  Main basin-hopping loop 
      JP=1
      DO J1=NDONE+1,NSTEPS 
! {{{
         ISTEP = J1

         CALL FLUSH(LFH)
         MCTEMP = TEMP(JP)
23       CONTINUE
         SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
         CALL TAKESTEP(JP)

            NQ(JP)=NQ(JP)+1
            CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)  
            NQTOT=NQTOT+1
            WRITE(LFH,100) &
                & 'Qu ',NQ(JP),&
                & ' E=',POTEL,&
                & ' steps=',ITERATIONS,&
                & ' RMS=',RMS,&
                & ' Markov E=',EPREV(JP),&
                & ' t=',TIME-TSTART
        
        100 format(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)

        CALL FLUSH(LFH)

        ! trackdata {{{
        IF (TRACKDATAT) THEN
             WRITE(MYEUNIT,'(I10,F20.10)') J1,POTEL
             WRITE(MYMUNIT,'(I10,G20.10)') J1,EPREV(JP)
             WRITE(MYBUNIT,'(I10,G20.10)') J1,QMIN(1)
             CALL FLUSH(MYEUNIT)
             CALL FLUSH(MYMUNIT)
             CALL FLUSH(MYBUNIT)
        ENDIF
        ! }}}
            IF (EVAP .AND. .NOT.EVAPREJECT) THEN
               NFAIL(JP)=NFAIL(JP)+1
               CALL MYRESET(JP,NATOMS,NPAR,NSEED)
               IF (DEBUG) THEN
                  WRITE(LFH,33) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
33                FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' EVAP,REJ')
               ENDIF
            ELSE
               CALL TRANSITION(POTEL,EPREV(JP),ATEST,JP,RANDOM,MCTEMP)

!  check: Markov energy agrees with COORDSO.{{{
!  Stop if not true.
!
               IF (DEBUG.OR.CHECKMARKOVT) THEN
                  CALL POTENTIAL(COORDSO(1:3*NATOMS,JP),GRAD,OPOTEL,.FALSE.,.FALSE.)
                  IF (ABS(OPOTEL-EPREV(JP)).GT.ECONV) THEN
                     IF (EVAP) THEN
                        WRITE(LFH,'(3(A,G20.10))') 'mc> WARNING - energy for saved coordinates ',OPOTEL,&
     &                     ' differs from Markov energy ',EPREV(JP),' because an atom moved outside the container'
                     ELSE
                        WRITE(LFH,'(2(A,G20.10))') 'mc> ERROR - energy for coordinates in COORDSO=',OPOTEL,&
     &                                                 ' but Markov energy=',EPREV(JP) 
                        STOP
                     ENDIF
                  ENDIF
               ENDIF

! }}}
!Accept or reject step {{{
!  If the quench did not converge then allow a
! potenial move, but count it as a rejection in terms of NSUCCESS and
! NFAIL. This way we will accept a lower minimum if found, but the steps won;t become so big.
! However, for TIP5P some cold fusion events that had not actually reached the threshold for
! rejection were actually accepted. Must prevent this!
               IF (ATEST) THEN
                  IF (DEBUG) THEN
                     WRITE(LFH,34) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
34                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF
                  IF ((J1-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV(JP)).GT.ECONV) THEN
!                    NRELAX=J1-JACCPREV
!                    IF (RESTART) WRITE(LFH,'(A,I6,A)') ' relaxation time set to ',NRELAX,' steps'
                     JACCPREV=J1
                  ENDIF
                  IF (QDONE.EQ.1) THEN
                     NSUCCESS(JP)=NSUCCESS(JP)+1
                  ELSE
                     NFAIL(JP)=NFAIL(JP)+1
                  ENDIF
                  EPPREV(JP)=EPREV(JP)
                  EPREV(JP)=POTEL
                  COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
                  VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
               ELSE
                  NFAIL(JP)=NFAIL(JP)+1
                  ! myreset:
                  !      coords=coordso
                  !      vat=vato
                  CALL MYRESET(JP,NATOMS,NPAR,NSEED)
                  IF (DEBUG) THEN
                     WRITE(LFH,36) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
36                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
               ENDIF
            ENDIF
            ! }}}

            IF ((MOD(J1,NACCEPT).EQ.0)) CALL ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
            TEMP(JP)=TEMP(JP)*SCALEFAC
            IF (HIT) GOTO 37
            IF (DUMPINT.GT.0) THEN
               IF (MOD(J1,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(J1,EBEST,BESTCOORDS,JBEST,JP)
               ENDIF
            ENDIF
           IF (NQ(JP).GT.NSTEPS) GOTO 37

         CALL FLUSH(LFH)
         ! }}}
         ENDDO
37    CONTINUE

      WRITE(LFH,21) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),&
     &               STEP(JP),ASTEP(JP),TEMP(JP)
21    FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)

      RETURN
!op226>}}} 
      END SUBROUTINE MC 



