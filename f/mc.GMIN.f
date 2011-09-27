
      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)
      ! declarations {{{
      USE COMMONS
      USE QMODULE , ONLY : QMIN, INTEQMIN

      USE PORFUNCS

      IMPLICIT NONE

      INTEGER J1, NSUCCESS(NPAR), NFAIL(NPAR), NFAILT(NPAR), NSUCCESST(NPAR), J2, NSTEPS, JP, J5, 
     1        UNT, ITERATIONS, NSUPERCOUNT, NQTOT, JACCPREV, NREN, NLAST, NSTEPREN, BRUN,QDONE,JBEST(NPAR),
     2        NRMS, NDONE, I, RNDSEED, J, NTOT, IMESG, ITRAJ, ITRAJO, NEACCEPT, J3, J4, ISTAT, LOCALCOUNT
      INTEGER :: NSYMCALL=0
      DOUBLE PRECISION POTEL, SCALEFAC, RANDOM, DPRAND, SAVECOORDS(3*NATOMS), TEMPCOORDS(3*NATOMS),
     1                 TIME, SPOTEL(NSUPER), SCOORDS(3*NATOMS,NSUPER), SCREENC(3*NATOMS),
     2                 EPPREV(NPAR), QSTART, QFINISH, RANNJ, RMIN, RMINO, RCOORDS(3*NATOMS),ELASTSYM(NPAR),
     3                 RCOORDSO(3*NATOMS), RVAT(NATOMS), RVATO(NATOMS), EPSSAVE, EBEST(NPAR),
     4                 BESTCOORDS(3*NATOMS,NPAR), endtime, RMSD, VINIT, CTE, TEMPTRAJ(0:NPAR-1),
     5                 T, BETA(0:NPAR-1), GRAD(3*NATOMS), E, ER, W, DELTA, DBETA, A9ANGLE, 
     &                 DUMMY1, DUMMY2, DUMMY3, INTE, OPOTEL
      LOGICAL CHANGEDE, EXCHANGEACCEPT, EXCHANGE, FLAG 
      LOGICAL CHIRALFAIL,AMIDEFAIL,GOODSTRUCTURE, LOGDUMMY, DISTOK, ATOMINGROUP(NATOMS)
      CHARACTER FNAME*9
      CHARACTER (LEN= 3)  ISTR
      CHARACTER (LEN=20) QUENCHNUM, QUNAME, DUMMYCHAR
c  AMH 
      INTEGER :: gly_count,iii,i2,i500,snapcount, DUMMYINT
      DOUBLE PRECISION prcord(NATOMS,3,3,3)
      DOUBLE PRECISION :: MCTEMP
!  csw34> PAIRDIST variables
      INTEGER :: PAIRCOUNTER
      DOUBLE PRECISION, EXTERNAL :: PAIRDISTANCE
      DOUBLE PRECISION :: ATOM1(3),ATOM2(3)

      LOGICAL EVAP, ATEST, STAY, evapreject, LOPEN
      COMMON /EV/ EVAP, evapreject
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT
      COMMON /Q4C/ QSTART, QFINISH

      character(len=10)       :: datechar,timechar,zonechar
      integer                 :: values(8),itime1
      double precision :: DISTGROUPX2,DISTGROUPY2,DISTGROUPZ2,DISTGROUPCENTRE,TESTANGLE
      integer :: J6
!op226>}}} 
!op226> Subroutine body {{{ 


      ALLOCATE(TMOVE(NPAR), OMOVE(NPAR))
      snapcount=0
      NSTEPREN=0
      EVAPREJECT=.FALSE.

      NDONE=0
            NQ(:)=NDONE

      IF (NACCEPT.EQ.0) NACCEPT=NSTEPS+1
      NRMS=0
      NLAST=0
      STAY=.FALSE.
      JACCPREV=0
      NQTOT=0
      RMINO=1.0D100
      RMIN=1.0D100
      NREN=NRENORM
      DO JP=1,NPAR 
         TMOVE(JP)=.TRUE.
         OMOVE(JP)=.TRUE.
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         NSUCCESST(JP)=0
         NFAILT(JP)=0
       ENDDO

C  Calculate the initial energy and save in EPREV
! {{{
      WRITE(MYUNIT,'(A)') 'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      DO JP=1,NPAR
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NQTOT=NQTOT+1
             WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART

C  EPREV saves the previous energy in the Markov chain.
C  EBEST and JBEST record the lowest energy since the last reseeding and the
C  step it was attained at. BESTCOORDS contains the corresponding coordinates.
 
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

!op226> GMIN_out: Starting MC run ...; Temperature will ... {{{ 
      WRITE(MYUNIT,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
      WRITE(MYUNIT,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'
!op226>}}} 

      NSUPERCOUNT=NSUPER

C  Main basin-hopping loop 
! {{{
      DO J1=NDONE+1,NSTEPS 
         ISTEP = J1

         CALL FLUSH(MYUNIT)
C
C  ********************************* Loop over NPAR parallel runs ******************************
C
         DO JP=1,NPAR 
C
            MCTEMP = TEMP(JP)
23          CONTINUE
 
            IF (EVAP .AND. .NOT.EVAPREJECT) THEN
               NFAIL(JP)=NFAIL(JP)+1
               CALL MYRESET(JP,NATOMS,NPAR,NSEED)
               IF (DEBUG) THEN
                  WRITE(MYUNIT,33) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
33                FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' EVAP,REJ')
               ENDIF
            ELSE
                  CALL TRANSITION(POTEL,EPREV(JP),ATEST,JP,RANDOM,MCTEMP)
!
!  Sanity check to make sure the Markov energy agrees with COORDSO. 
!  Stop if not true.
!
               IF (DEBUG.OR.CHECKMARKOVT) THEN
                  CALL POTENTIAL(COORDSO(:,JP),GRAD,OPOTEL,.FALSE.,.FALSE.)
                  IF (ABS(OPOTEL-EPREV(JP)).GT.ECONV) THEN
                     IF (EVAP) THEN
                        WRITE(MYUNIT,'(3(A,G20.10))') 'mc> WARNING - energy for saved coordinates ',OPOTEL,
     &                     ' differs from Markov energy ',EPREV(JP),' because an atom moved outside the container'
                     ELSE
                        WRITE(MYUNIT,'(2(A,G20.10))') 'mc> ERROR - energy for coordinates in COORDSO=',OPOTEL,
     &                                                 ' but Markov energy=',EPREV(JP)
                        STOP
                     ENDIF
                  ENDIF
               ENDIF
C Accept or reject step. 
               IF (ATEST) THEN
                  IF (DEBUG) THEN
                     WRITE(MYUNIT,34) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
34                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF
                  IF ((J1-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV(JP)).GT.ECONV) THEN
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
                  CALL MYRESET(JP,NATOMS,NPAR,NSEED)
                  IF (DEBUG) THEN
                     WRITE(MYUNIT,36) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
36                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
               ENDIF
            ENDIF
C  Check the acceptance ratio.
            IF ((MOD(J1,NACCEPT).EQ.0).AND.(NSEED.EQ.0)) CALL ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)

            TEMP(JP)=TEMP(JP)*SCALEFAC

            IF (HIT) GOTO 37
            IF (DUMPINT.GT.0) THEN
               IF (MOD(J1,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(J1,EBEST,BESTCOORDS,JBEST,JP)
               ENDIF
            ENDIF
           IF (NQ(JP).GT.NSTEPS) GOTO 37

         ENDDO
         CALL FLUSH(MYUNIT)
      ENDDO
! }}}

37    CONTINUE
      WRITE(MYUNIT,21) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),
     1               STEP(JP),ASTEP(JP),TEMP(JP)
21          FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
     
      RETURN
!op226>}}} 
      END
 

