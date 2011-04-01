
      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)

      USE COMMONS
      USE qmodule , only : qmin, INTEQMIN
      USE modcharmm
      USE modamber9, only : mdstept,cisarray1,cisarray2,chiarray1,chiarray2,nocistransdna,nocistransrna,
     &                      setchiral,amchpmax,doligmove,ligmovefreq, grouprott, grouprotfreq, dogrouprot,
     &                      amchnmax
      USE porfuncs
      USE AMHGLOBALS, ONLY: NMRES,OMOVI,AVEP,NUMPRO,IRES
      USE AMH_INTERFACES, ONLY:E_WRITE

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
      DOUBLE PRECISION :: mctemp
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

! Write a list of FROZEN atoms for use in an (o)data file
!op226> IF (FREEZEGROUPT) THEN {{{
      IF (FREEZEGROUPT) THEN
         OPEN(UNIT=4431,FILE='frozen.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         DO J6=1,NATOMS
!
! Work out the distance from GROUPCENTRE to the current atom J1
! 
            DISTGROUPX2=(COORDS(3*GROUPCENTRE-2,1)-COORDS(3*J6-2,1))**2
            DISTGROUPY2=(COORDS(3*GROUPCENTRE-1,1)-COORDS(3*J6-1,1))**2
            DISTGROUPZ2=(COORDS(3*GROUPCENTRE  ,1)-COORDS(3*J6  ,1))**2
            DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
! If working in GT mode (default), FREEZE all atoms >GROUPRADIUS from the GROUPCENTRE atom
            IF((FREEZEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
               NFREEZE=NFREEZE+1
               FROZEN(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'FREEZE ',J6
! If working in LT mode, FREEZE all atoms <GROUPRADIUS from the GROUPCENTRE atom
            ELSE IF((FREEZEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
               NFREEZE=NFREEZE+1
               FROZEN(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'FREEZE ',J6
            END IF
         END DO
         CLOSE(4431)
! Prevent it doing this again
         FREEZEGROUPT=.FALSE.     
      ENDIF
!op226>}}} 

! Write a list of DONTMOVE atoms for use in an (o)data file
!op226> IF (DONTMOVEGROUPT) THEN {{{
      IF (DONTMOVEGROUPT) THEN
              OPEN(UNIT=4431,FILE='dontmove.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         DO J6=1,NATOMS
!
! Work out the distance from DONTMOVECENTRE to the current atom J1
! 
            DISTGROUPX2=(COORDS(3*DONTMOVECENTRE-2,1)-COORDS(3*J6-2,1))**2
            DISTGROUPY2=(COORDS(3*DONTMOVECENTRE-1,1)-COORDS(3*J6-1,1))**2
            DISTGROUPZ2=(COORDS(3*DONTMOVECENTRE  ,1)-COORDS(3*J6  ,1))**2
            DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
! If working in GT mode (default), DONTMOVE all atoms >GROUPRADIUS from the DONTMOVECENTRE atom
            IF((DONTMOVEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
               NDONTMOVE=NDONTMOVE+1
               DONTMOVE(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'DONTMOVE ',J6
! IF working in LT mode, DONTMOVE all atoms <GROUPRADIUS from the DONTMOVECENTRE atom
       ELSE IF((DONTMOVEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
               NDONTMOVE=NDONTMOVE+1
               DONTMOVE(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'DONTMOVE ',J6
            END IF
         END DO
         CLOSE(4431)
! Prevent it doing this again
         DONTMOVEGROUPT=.FALSE.     
      ENDIF
!op226>}}} 
      
!     csw34> Added defaults to prevent accidentaly discarding
!     structures for AMH      

      CHIRALFAIL=.FALSE.
      AMIDEFAIL=.FALSE.
      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'mc> A ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      ALLOCATE(TMOVE(NPAR), OMOVE(NPAR))
      snapcount=0
      NSTEPREN=0
      EVAPREJECT=.FALSE.
      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'mc> B ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

#ifdef MPI
      IF (MPIT) THEN 
         IF (DEBUG) WRITE(MYUNIT,'(A,I6)') 'MPIERR=',MPIERR
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPAR,MPIERR)
         IF (DEBUG) WRITE(MYUNIT, '(A,2I6)') 'NPAR,MPIERR=',NPAR,MPIERR
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
         JP=MYNODE+1
         ITRAJ=MYNODE
         IF (DEBUG) WRITE(MYUNIT, '(A,3I6)') 'In mc after MPI_MPI_COMM_RANK MPIERR,MYNODE,JP=',MPIERR,MYNODE,JP
      ENDIF
      NEACCEPT=0
#endif

      NDONE=0
      IF (RESTORET) THEN
#ifdef MPI
         CALL RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
#else
         DO JP=1,NPAR
            CALL RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
         ENDDO
#endif
         WRITE(MYUNIT, '(A,I10)') 'MC> restore NDONE=',NDONE
!     csw34> Sets the quench counter so that the GMIN_out file makes sense after using RESTORE!
!         NQ(:)=NDONE
      ENDIF
      NQ(:)=NDONE

C tvb requesting a basin-sampling MC run: {{{
      
      IF (BSWL.and.(.not.TETHER)) then
         CALL BasinSampling
         RETURN
      ELSEIF (TETHER) THEN
         CALL TetheredWL
         RETURN
      ENDIF
C }}}

#ifdef MPI
      WRITE(MYUNIT, '(A,I10,A,I10,A)') "Processor", mynode+1, " of", NPAR, " speaking:"
      WRITE(MYUNIT, '(A,I10)') 'Number of atoms', natoms
#endif

      IF (NACCEPT.EQ.0) NACCEPT=NSTEPS+1
      NRMS=0
      NLAST=0
      STAY=.FALSE.
      JACCPREV=0
      NQTOT=0
      RMINO=1.0D100
      RMIN=1.0D100
      NREN=NRENORM
#ifdef MPI
#else
      DO JP=1,NPAR 
#endif
         TMOVE(JP)=.TRUE.
         OMOVE(JP)=.TRUE.
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         NSUCCESST(JP)=0
         NFAILT(JP)=0
         IF (JDUMP(JP).AND.(.NOT.NEWJUMP)) THEN
            WRITE(FNAME,'(A,I1)') 'ebuffer.',JP
            UNT=70+JP
            OPEN(UNIT=UNT,FILE=FNAME,STATUS='UNKNOWN')
            WRITE(FNAME,'(A,I1)') 'cbuffer.',JP
            UNT=70+NPAR+JP
            OPEN(UNIT=UNT,FILE=FNAME,STATUS='UNKNOWN')
         ENDIF
#ifdef MPI
#else
      ENDDO
#endif

      IF (AMHT) THEN
         write(omovi,1334)nmres,3,1,INT(real(NSTEPS)/real(NINT_AMH))
1334     format(4(i8,1x),' nmres nmcrd numpro nmsnap')
      ENDIF
    
      IF (.NOT.RESTORET) THEN
!       csw34> Set the centre of mass to be at the specified location
!       contained in the SETCENTRE X Y Z keyword
         IF (SETCENT) CALL SETCENTRE(COORDS)
C
C For MAKEOLIGOT and MAKEOLIGOSTART=TRUE: generate oligomers by placing new segments.
         IF (CHRMMT.AND.MAKEOLIGOT.AND.MAKEOLIGOSTART) THEN
#ifdef MPI
C        seed the random number generator with system time  + MYNODE (for MPI runs)
             IF (RANDOMSEEDT) THEN
                CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
                itime1= values(6)*60 + values(7)
                CALL SDPRND(itime1+MYNODE)
             ENDIF
             CALL CHMAKEOLIGOMER(JP)
             call flush(6)
#else
             DO JP=1,NPAR
                CALL CHMAKEOLIGOMER(JP)
             ENDDO
#endif
         ENDIF
      ENDIF

C  Calculate the initial energy and save in EPREV
!op226>{{{ 
#ifdef MPI
!op226> MPI section {{{ 
      WRITE(MYUNIT, '(A)')  'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
      NQTOT=NQTOT+1
      WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
      CALL FLUSH(MYUNIT,ISTAT)

 
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
      EPSSPHERE=EPSSAVE

! Initialisation 

      IF (PTTMIN < 0.000001D0 ) PTTMIN = 0.000001D0 ! to avoid devision by zero
      CTE=(LOG(PTTMAX/PTTMIN))/(NPAR-1)
      CTE=EXP(CTE)

      DO I=0, NPAR-1
         TEMPTRAJ(I)=PTTMIN*CTE**I
         T=TEMPTRAJ(I)
         BETA(I)=1.0D0/T
      ENDDO
      DO I=1, NPAR
         TEMP(I)=TEMPTRAJ(I-1)
      ENDDO
      CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE.)
      VINIT=POTEL

      WRITE(MYUNIT, '(A, 2G20.10)') 'Temperature range', TEMPTRAJ(0), TEMPTRAJ(NPAR-1)
      WRITE(MYUNIT, '(A, G20.10)') 'For this replica T=', TEMPTRAJ(MYNODE)
      WRITE(MYUNIT, '(A, G20.10)') 'Starting potential energy=', VINIT

      RNDSEED=2002+MYNODE
      CALL SDPRND(RNDSEED)
      RANDOM=DPRAND()
      WRITE(MYUNIT, '(A, G20.10)') 'Starting random number=', RANDOM
      WRITE(ISTR,'(I3)') JP
      TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
      IF ((CHRMMT).AND.(.NOT.RESTORET)) CALL CHARMMDUMP(TEMPCOORDS,'initialmin.'//TRIM(ADJUSTL(ISTR)))

!op226>}}} 
#else
      WRITE(MYUNIT,'(A)') 'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      DO JP=1,NPAR
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NQTOT=NQTOT+1
         IF (NPAR.GT.1) THEN
            WRITE(MYUNIT,'(A,I2,A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,']Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ELSE
            WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ENDIF


C Added dump of the initial structure, this is very useful for flu! csw34
      IF (CHRMMT) CALL CHARMMDUMP(COORDS,'initialmin')

C     csw34> Added initial call to check_cistrans_protein to store cis/trans info for initial structure
         IF (AMBERT.AND.NOCISTRANS.AND.(.NOT.NOCISTRANSDNA).AND.(.NOT.NOCISTRANSRNA)) THEN
            WRITE(MYUNIT,'(A)') ' mc> Storing cis/trans information for initial structure'
            CALL check_cistrans_protein(COORDS(:,1),NATOMS,LOGDUMMY,MINOMEGA,cisarray1)
         ENDIF
C     abc> Added initial call to set_check_chiral for L/D info for initial structure
         IF (AMBERT.AND.SETCHIRAL.AND.NOCISTRANS.AND.(.NOT.NOCISTRANSDNA).AND.(.NOT.NOCISTRANSRNA)) THEN
            WRITE(MYUNIT,'(A)') ' mc> Storing chiral information for initial structure'
            CALL set_check_chiral(COORDS(:,1),NATOMS,LOGDUMMY,chiarray1)
         ENDIF 
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
#endif
!op226>}}} 

!op226> IF (THOMSONT) THEN - Thomson problem {{{
      IF (THOMSONT) THEN
!
! Scale maximum step size for the Thomson problem according to the mean nearest-neighbour 
! distance after the first quench.
!
         DUMMY1=0.0D0
         DO J1=1,NATOMS
            DUMMY2=1.0D100
            DO J2=1,NATOMS
               IF (J2.EQ.J1) CYCLE
               DUMMY3=(COORDS(3*(J1-1)+1,1)-COORDS(3*(J2-1)+1,1))**2 +  
     &                (COORDS(3*(J1-1)+2,1)-COORDS(3*(J2-1)+2,1))**2 + 
     &                (COORDS(3*(J1-1)+3,1)-COORDS(3*(J2-1)+3,1))**2 
               IF (DUMMY3.LT.DUMMY2) DUMMY2=DUMMY3
            ENDDO
            DUMMY1=DUMMY1+DUMMY2
         ENDDO
         DUMMY1=SQRT(DUMMY1/NATOMS)
         DO J1=1,NPAR
            STEP(J1)=STEP(J1)*DUMMY1
         ENDDO
         WRITE(MYUNIT, '(2(A,G20.10))') 'Maximum step size scaled by mean nearest neighbour distance of ',DUMMY1,' to ',STEP(1)
      ENDIF
!op226>}}} 

!op226> GMIN_out: Starting MC run ...; Temperature will ... {{{ 
      IF (NPAR.EQ.1) THEN
         WRITE(MYUNIT,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
      ELSE
         WRITE(MYUNIT,'(A,I3,A,I10,A)') 'Starting ',NPAR,' parallel MC runs of ',NSTEPS,' steps'
      ENDIF
      WRITE(MYUNIT,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'
!op226>}}} 

!op226> csw34> IF (LOCALSAMPLET.AND.AMBERT) THEN  {{{
! csw34> Before we start BH steps, 
! check the structure satisfies the conditions in LOCALSAMPLE if specified
!        This is important to prevent an infinite loop occuring!
!        Also, make sure that the ligand has been specified in movableatoms
      IF (LOCALSAMPLET.AND.AMBERT) THEN
         IF(NMOVABLEATOMS==0) THEN
            WRITE(MYUNIT,*) 'must have MOVABLEATOMS specified when using LOCALSAMPLE' 
            STOP
         ENDIF
         DISTOK=.FALSE.
         CALL A9DISTCHECK(COORDS(:,JP),DISTOK)
         IF (.NOT.DISTOK) THEN
            WRITE(MYUNIT,*) 'initial structure violates LOCALSAMPLE conditions - exiting!'
            STOP
         ENDIF
         DISTOK=.FALSE.
      ENDIF   
!op226>}}} 
      NSUPERCOUNT=NSUPER

C  Main basin-hopping loop 
! {{{
      DO J1=NDONE+1,NSTEPS 
         ISTEP = J1

         CALL FLUSH(MYUNIT)
         IF (NEWJUMP) RANNJ=DPRAND()
C
C  ********************************* Loop over NPAR parallel runs ******************************
C
#ifdef MPI
#else
         DO JP=1,NPAR 
#endif
!       csw34> If QUCENTRE is specified, move the centre of coordinates
!       to (0,0,0) before taking the next step (improve this so that you
!       can specify where to move the centre like SETCENTRE?
            IF (QUCENTRE) THEN 
!       csw34> David mentioned a possible compiler bug causing problems
!       when you just pass a part of the COORDS array so reading the
!       right bit into TEMPCOORDS first and then back out.
               TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
               CALL CENTRE2(TEMPCOORDS)
               COORDS(1:3*NATOMS,JP)=TEMPCOORDS(1:3*NATOMS)
            ENDIF
            IF (RIGID.AND.(BLOCK(JP).GT.0)) THEN
               IF (MOD(J1-1,BLOCK(JP)).EQ.0) THEN
                  IF (MOD((J1-1)/BLOCK(JP),2).EQ.0) THEN
                     WRITE(MYUNIT,'(A,I6,A)') 'Starting a block of ',BLOCK(JP),' rigid body translational moves'
                     TMOVE(JP)=.TRUE.
                     OMOVE(JP)=.FALSE.
                  ELSE 
                     WRITE(MYUNIT,'(A,I6,A)') 'Starting a block of ',BLOCK(JP),' rigid body angular moves'
                     OMOVE(JP)=.TRUE.
                     TMOVE(JP)=.FALSE.
                  ENDIF
               ENDIF
            ENDIF
C
C  MAM1000> The default temperature used for the MC acceptance criterion is the one derived in
C           the initialisation section above.  MCTEMP is passed to the subroutine TRANSITION, where
C           the acceptance/rejection decision is made.  However, individual MC moves can override
C           this temperature by setting MCTEMP to a different value for the current step.
            MCTEMP = TEMP(JP)
C 
C  Jump-moves.
C 
            IF (JUMPMOVE(JP).AND.(MOD(J1,JUMPINT(JP)).EQ.0)) CALL JUMPM(RANNJ,J1,JP,EPPREV)
C 
C  Ordinary steps.
C 
23          CONTINUE
C
C Don;t call symmetry unless the minimum in the Markov chain has changed.
C We should really check if the minimum has changed since the last call to SYMMETRY,
C which can be done with ABS(ELASTSYM(JP)-EPREV(JP)) if NSYMINTERVAL=1.
C
C           IF (SYMMETRIZE.AND.(MOD(J1-1,NSYMINTERVAL).EQ.0).AND.(ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV)) THEN
C           WRITE(MYUNIT,'(A,3G20.10)') 'ELASTSYM,EPREV,diff=',ELASTSYM(JP),EPREV(JP),ABS(ELASTSYM(JP)-EPREV(JP))
            IF (SYMMETRIZE.AND.(MOD(NQ(JP),NSYMINTERVAL).EQ.0).AND.
     &               ((SYMMETRIZECSM).OR.(ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV))) THEN
               IF ((ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV)) NSYMCALL=0
               ELASTSYM(JP)=EPREV(JP)
C              IF ((.NOT.MOVESHELLT).OR.(NSURFMOVES(JP).LT.0)) THEN
                  IF (SYMMETRIZECSM) THEN
                     CALL SYMMETRYCSM(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL,  
     &                                J1,NSUCCESS,NFAIL,EBEST,BESTCOORDS,JBEST,EPPREV)
                  ELSE
                     CALL SYMMETRY(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL,  
     &                                J1,EBEST,BESTCOORDS,JBEST,EPPREV)
                  ENDIF
                  IF (HIT) GOTO 37 
!              ELSE
!                 CALL SYMMETRY2(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL)
!              ENDIF
!              WRITE(MYUNIT,'(A,I2,A,2I6)') '[',JP,']mc> NCORE: ',NCORE(1:NPAR)
!              IF (HIT) GOTO 37 ! hit cannot change in symmetry2 
C
C  Check for reseeding.
C 
               POTEL=EPREV(JP) ! NEWRES assumes POTEL is the energy of the current structure in COORDS
               IF (CHANGEDE.AND.NEWRESTART) THEN
                  CALL NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,
     1                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST)
               ENDIF
            ELSEIF (ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV) THEN ! Markov minimum has changed, but SYMMETRY not called
               NSYMREM=0                                       ! Should therefore reset NSYMREM.
            ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! START OF STEP TAKING CALLS!                                                                            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! csw34> Coordinates are saved so that moves can be undone
!
                  SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
! csw34> If you want to look at the effect of moves, you can dump out
! the structure BEFORE the move here.
!                 CALL A9DUMPPDB(COORDS(:,JP),"beforemove")
!                 CALL CHARMMDUMP(COORDS(:,JP),'beforemove')

!
! CHARMM STEP TAKING
!
               IF (CHRMMT) THEN
                  IF (CHMDT.AND.MOD(J1,CHMDFREQ).EQ.0) THEN
                     CALL CHMD(JP)
                  ELSE
!                    CALL CHARMMDUMP(COORDS(1:3*NATOMS,JP),'beforestep')
                     IF (CHRIGIDTRANST.AND.MOD(J1,FTRANS).EQ.0) CALL MKRIGIDTRANS(JP)
                     IF (CHRIGIDROTT.AND.MOD(J1,FROT).EQ.0) CALL MKRIGIDROT(JP)
                     IF (MOD(J1,CHFREQ).EQ.0) CALL TAKESTEPCH(JP)
! bs360: check the perturbed structure prior minimization (for debugging reasons)
!#ifdef MPI
!                     CALL CHARMMDUMP(COORDS(1:3*NATOMS,JP),'afterstep.'//TRIM(ADJUSTL(ISTR)))
!#else
!                     CALL CHARMMDUMP(COORDS(1:3*NATOMS,JP),'afterstep')
!#endif
C               seed the random number generator with system time  + MYNODE (for MPI runs)
                     IF (RANDOMSEEDT) THEN
                         CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
                         itime1= values(6)*60 + values(7)
                         CALL SDPRND(itime1+MYNODE)
                     END IF
                              
! jmc49> only do this undoing if the steps are in Cartesians (i.e. CHNMAX <= 0.0D0), not internals, 
!        in order to allow rigid-body movements of the frozen parts relative to each other when the moves 
!        are in internals
                  ENDIF
!
! ANISOTROPIC STEP TAKING
!
               ELSE IF (ELLIPSOIDT.AND..NOT.PYGPERIODICT.AND..NOT.PYBINARYT) THEN
                CALL TAKESTEPELLIPSOIDS(JP) 

               ELSE IF (GAYBERNEDCT) THEN
                CALL TAKESTEPGB(JP)

               ELSE IF (MULTISITEPYT) THEN
C               seed the random number generator with system time + MYNODE (for MPI runs)
                  IF(RANDOMSEEDT) THEN
                     CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
                     itime1= values(6)*60 + values(7)
                     CALL SDPRND(itime1+MYNODE)
                  END IF
                CALL TAKESTEPMULTISITEPY(JP)

               ELSE IF (PYGPERIODICT.OR.PYBINARYT.OR.LJCAPSIDT) THEN
C               seed the random number generator with system time + MYNODE (for MPI runs)
                  IF(RANDOMSEEDT) THEN
                     CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
                     itime1= values(6)*60 + values(7)
                     CALL SDPRND(itime1+MYNODE)
                  END IF
                  IF(SWAPMOVEST) THEN
                     CALL TAKESTEPSWAPMOVES(JP)
                     CALL TAKESTEPELPSD(JP)
                  ELSE
                     CALL TAKESTEPELPSD(JP)
                  END IF

!               ELSE IF (DBPT) THEN
!                  CALL TAKESTEPDB(JP)

!
! MARK'S MIXED CLUSTER STEP TAKING
!
C     MAM1000 >  Charged/neutral particle swaps for LJ+Coulomb mixed clusters
               ELSE IF (LJCOULT) THEN
                  RANDOM=DPRAND()
                  IF (RANDOM < COULSWAP) THEN
                     CALL TAKESTEPLJC(JP)
                     MCTEMP = COULTEMP
                  ELSE
                     CALL TAKESTEP(JP)
                  END IF

!               ELSE IF (LWOTPT) THEN
!                  CALL TAKESTEPLWOTP(JP)

               ELSE IF (GBT .OR. GBDT .OR. GBDPT .OR. PYGT .OR. PYGDPT) THEN
                  CALL TKSTDCELPSD(JP)
      
               ELSE IF (MSGBT) THEN
                  CALL TAKESTEPMSGB(JP) 

               ELSE IF (MSPYGT) THEN
                  CALL TAKESTEPMSPY(JP)
!
! AMBER STEP TAKING
!               
               ELSE IF (AMBERT) THEN
                  IF (.NOT.ALLOCATED(MOVABLEATOMLIST)) THEN 
                     ALLOCATE(MOVABLEATOMLIST(NATOMS))
                     NMOVABLEATOMS=NATOMS
                  ENDIF
                  DISTOK=.FALSE.
                  LOCALCOUNT=0
                  DO WHILE (.NOT.DISTOK)
                     LOCALCOUNT=LOCALCOUNT+1
                     IF (LIGMOVET.AND.MOD(J1,ligmovefreq).EQ.0) THEN
                        doligmove=.TRUE.
                     ENDIF
                     IF(GROUPROTT.AND.MOD(J1,grouprotfreq).EQ.0) THEN
                        dogrouprot=.TRUE.
                     ENDIF
                     CALL TAKESTEPAMBER(JP,COORDS(:,JP),MOVABLEATOMLIST,NMOVABLEATOMS,LIGMOVET,MDSTEPT,RANDOMSEEDT)
                     doligmove=.FALSE.
                     dogrouprot=.FALSE.
C
C  seed the random number generator with system time + MYNODE (for MPI runs)
C
                     IF (RANDOMSEEDT) THEN
                        CALL DATE_AND_TIME(DATECHAR,TIMECHAR,ZONECHAR,VALUES)
                        ITIME1= VALUES(6)*60 + VALUES(7)
                        CALL SDPRND(ITIME1+MYNODE)
                     ENDIF
                     IF (AMCHPMAX.EQ.0) THEN
                        CALL TAKESTEP(JP)
                     ELSE
C     msb50> New AMBER dihedral move routine
                        CALL TAKESTEPAMM(COORDS(:,JP), DEBUG, STEP(JP))
                     ENDIF
C     csw34> Check distances for groups defined in movableatoms file
C            If A->B > ABTHRESH or A->C > ACTHRESH, the step is discarded and we try again 
                     IF (LOCALSAMPLET) THEN
                        CALL A9DISTCHECK(COORDS(:,JP),DISTOK)
                        IF (.NOT.DISTOK) COORDS(:,JP)=SAVECOORDS(:)
                     ELSE
                        DISTOK=.TRUE.
                     ENDIF
                  ENDDO
!
! ALL OTHER STEP TAKING
!
               ELSE
 
C  These coordinates are overwritten if we try to call MPI_IPROBE
C
C                 WRITE(MYUNIT,'(I6)') NATOMS
C                 WRITE(MYUNIT,'(A,G20.10)') 'eprev=',EPREV
C                 WRITE(MYUNIT,'(A,3G20.10)') ('LA ',COORDS(3*(J3-1)+1:3*(J3-1)+3,JP),J3=1,NATOMS)
                  CALL TAKESTEP(JP)
               ENDIF
! Restore atom coordinates if atom is FROZEN or DONTMOVE as long as
! we're not taking internal coordinate moves in CHARMM or AMBER
               IF ((CHNMAX.LE.0.0D0).AND.(AMCHNMAX.LE.0.0D0)) THEN
                  IF(FREEZE) THEN
                     DO J2=1,NATOMS
                        IF (FROZEN(J2)) THEN
                           COORDS(3*(J2-1)+1:3*(J2-1)+3,JP)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                        ENDIF
                     ENDDO
                  ENDIF
!
! Same for DONTMOVE atoms
! 
                  IF(DONTMOVET) THEN
                     DO J2=1,NATOMS
                        IF (DONTMOVE(J2)) THEN
                           COORDS(3*(J2-1)+1:3*(J2-1)+3,JP)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
               
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END OF STEP TAKING CALLS!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! csw34> If you want to look at the effect of moves, you can dump out
! the structure AFTER the move here.
!              CALL A9DUMPPDB(COORDS(:,JP),"aftermove")
!              CALL CHARMMDUMP(COORDS(:,JP),"aftermove")
            
! KEYWORD <DUMPSTEPS> BLOCK            
! csw34> Dump the coordinates after every step in AMBER pdb and rst format 
                  IF (DUMPSTEPST) THEN
                      WRITE(QUENCHNUM,*) NQ(JP)
                      QUNAME='afterstep'//TRIM(ADJUSTL(QUENCHNUM))//'.rst'
                      OPEN(UNIT=20,FILE=QUNAME,STATUS='UNKNOWN')  
                      WRITE(20,'(a20)') QUNAME
                      WRITE(20,'(i5)') NATOMS
                      WRITE(20,'(6f12.7)') COORDS(:,JP) 
                      CLOSE(20)
! csw34> Dump to PDB using routine in amberinterface.f
                      QUNAME='afterstep'//TRIM(ADJUSTL(QUENCHNUM))
                      CALL A9DUMPPDB(COORDS(:,JP),QUNAME)
                  ENDIF
! END OF KEYWORD <DUMPSTEPS> BLOCK

               NQ(JP)=NQ(JP)+1
               IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
               CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)  
               NQTOT=NQTOT+1
 
C  Check for results of taboo list. SELFT is set in taboo.
 
               IF (SELFT) THEN
                  CALL MYRESET(JP,NATOMS,NPAR,NSEED)
                  IF (DEBUG) THEN
                     WRITE(MYUNIT,'(A)') 'Taboo list:'
                     WRITE(MYUNIT,'(6F20.10)') (ESAVE(J2,1),J2=1,NT(1))
                     WRITE(MYUNIT,73) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
73                   FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' TABOO')
                  ENDIF
                  GOTO 23
               ENDIF
   
C
C  Output
C
               IF (NPAR.GT.1) THEN
                  WRITE(MYUNIT,'(A,I2,A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,']Qu ',NQ(JP),' E=',
     1                 POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART
               ELSE
                  WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1                 POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART

!     mp466>  writes structure and energetic data at regular increments
!             to *plot and movie files for AMH potential
!
                  IF ((MOD(J1,NINT_AMH).EQ.0).AND.AMHT) THEN

                     GLY_COUNT = 0
                     SNAPCOUNT = SNAPCOUNT + 1

                     DO III = 1,NMRES
                        IF (IRES(III).EQ.8) THEN
                           PRCORD(III,1,1,1) = COORDS(9*(III-1)+1- GLY_COUNT*3,JP) !  CA X
                           PRCORD(III,2,1,1) = COORDS(9*(III-1)+2- GLY_COUNT*3,JP) !  CA Y
                           PRCORD(III,3,1,1) = COORDS(9*(III-1)+3- GLY_COUNT*3,JP) !  CA Z
!    SWAP  CA for CB
                           PRCORD(III,1,1,2) = COORDS(9*(III-1)+1- GLY_COUNT*3,JP) !  CB X
                           PRCORD(III,2,1,2) = COORDS(9*(III-1)+2- GLY_COUNT*3,JP) !  CB Y
                           PRCORD(III,3,1,2) = COORDS(9*(III-1)+3- GLY_COUNT*3,JP) !  CB Z
                           PRCORD(III,1,1,3) = COORDS(9*(III-1)+4- GLY_COUNT*3,JP) !  O X
                           PRCORD(III,2,1,3) = COORDS(9*(III-1)+5- GLY_COUNT*3,JP) !  O Y
                           PRCORD(III,3,1,3) = COORDS(9*(III-1)+6- GLY_COUNT*3,JP) !  O Z
                           GLY_COUNT = GLY_COUNT +1
                        ELSE
                           PRCORD(III,1,1,1) = COORDS(9*(III-1)+1 - GLY_COUNT*3,JP) !  CA X
                           PRCORD(III,2,1,1) = COORDS(9*(III-1)+2 - GLY_COUNT*3,JP) !  CA Y
                           PRCORD(III,3,1,1) = COORDS(9*(III-1)+3 - GLY_COUNT*3,JP) !  CA Z
                           PRCORD(III,1,1,2) = COORDS(9*(III-1)+4 - GLY_COUNT*3,JP) !  CB X
                           PRCORD(III,2,1,2) = COORDS(9*(III-1)+5 - GLY_COUNT*3,JP) !  CB Y
                           PRCORD(III,3,1,2) = COORDS(9*(III-1)+6 - GLY_COUNT*3,JP) !  CB Z
                           PRCORD(III,1,1,3) = COORDS(9*(III-1)+7 - GLY_COUNT*3,JP) !  O X
                           PRCORD(III,2,1,3) = COORDS(9*(III-1)+8 - GLY_COUNT*3,JP) !  O Y
                           PRCORD(III,3,1,3) = COORDS(9*(III-1)+9 - GLY_COUNT*3,JP) !  O Z
                        ENDIF
                     ENDDO

                     WRITE(OMOVI,683)1,SNAPCOUNT,1,TEMP(JP),1
683                  FORMAT(3(i6,1x),f8.4,1x,i5,' stuct snap t T Tid')

                     DO  I500=1,NMRES
                         WRITE(OMOVI,332)(PRCORD(I500,I2,1,1),I2=1,3),(PRCORD(I500,I2,1,2),I2=1,3),(PRCORD(I500,I2,1,3),I2=1,3)
                     ENDDO
332                  FORMAT('CA: ',3(F8.3,1X),'CB: ',3(F8.3,1X),'OX: ', 3(F8.3,1x))

                     !CALL E_WRITE(AVEP(:,:,:),TEMP(JP),NUMPRO,SNAPCOUNT)
                  ENDIF  !                  IF (MOD(J1,NINT_AMH).EQ.0)
               ENDIF     !                  IF (NPAR.GT.1)
               CALL FLUSH(MYUNIT)

C  RMS compared to reference structure 'compare'
          IF (RMST.AND.CHRMMT) THEN
             CALL CHRMS(JP,RMSD)
             IF (DEBUG) WRITE(MYUNIT,'(A,F15.5)')'RMSD = ',RMSD
             IF (RMSD.LE.RMSLIMIT) CALL SAVERMS(JP,POTEL,RMSD)
C                NRMS=NRMS+1
C                WRITE(CNRMS,'(I6)') NRMS
C                CALL CHARMMDUMP(COORDS(1:3*NATOMS,JP),'rms.'//TRIM(ADJUSTL(CNRMS)))
C                OPEN(UNIT=20,FILE='rms.'//TRIM(ADJUSTL(CNRMS)),POSITION='APPEND',STATUS='OLD')
C                WRITE(20,'(A,I6,A,F15.5,A,F15.5)') '*   Qu ',NQ(JP),' E=',POTEL,' RMSD=',RMSD
C                CLOSE(20)
C            ENDIF
          ENDIF

!
!     csw34> PAIRDIST prints the distance between specified atom pairs
!     after each quench. 
!
          IF (PAIRDISTT) THEN
!     Write end of previous line as using ADVANCE="NO"
                WRITE(MYPUNIT,*) " "
!     Write current quench number
                WRITE(MYPUNIT,'(I10)',ADVANCE="NO") NQ(JP)
!     For each pair, assign ATOM1 and ATOM2 arrays containing coordinates
             DO PAIRCOUNTER=1,NPAIRS
                ATOM1(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,1)-2:3*PAIRDIST(PAIRCOUNTER,1),JP)
                ATOM2(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,2)-2:3*PAIRDIST(PAIRCOUNTER,2),JP)
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
             WRITE(MYMUNIT,'(I10,G20.10)') J1,EPREV(JP)
             WRITE(MYBUNIT,'(I10,G20.10)') J1,QMIN(1)
             CALL FLUSH(MYEUNIT)
             CALL FLUSH(MYMUNIT)
             CALL FLUSH(MYBUNIT)
             IF (A9INTET.AND.(NQ(JP).GT.2)) THEN
                WRITE(3999,'(I10,G20.10)') NQ(JP),INTEQMIN(1)
                CALL FLUSH(3999)
             ENDIF
!
!     csw34> if RMS is also specified, prints the RMSD from the
!     comparison structure to the file 'rmsd'. Which RMSD depends on the
!     final arguement of the RMS keyword - see documentation!
!
             IF (RMST.AND.CHRMMT) THEN
                WRITE(4428,'(I10,F15.5)') NQ(JP),RMSD
                CALL FLUSH(4428)
             ENDIF
          ENDIF

C DAESTAT keyword prints two files :
C stat.all which contains all quenches and their energies
C stat.acc which contains only accepted quenches and their energies
C used to analyse how many new minima are being found
 
            IF (DAESTAT) THEN
               PRINT*,'DAESTAT block in mc.f not implemented'
               STOP
C              CALL CALCMIND(JP,MIND)
C              CALL CALCDIHE(DIHE)
C              WRITE(MYUNIT,'(A,I6,3F20.10)') 'NQALL POTEL',NQ(JP),POTEL,MIND,DIHE
C              WRITE(36,'(I6,3F20.10)') NQ(JP),POTEL,MIND,DIHE
            ENDIF
 
C  Check for reseeding.
 
            IF (NEWRESTART.AND.(.NOT.SEEDT)) THEN 
              CALL NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,RMIN,RVAT,BRUN,SCREENC,QDONE,
     &                    JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DJW
!             IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,JP))
!             COORDSO(1:3*(NATOMS-NSEED),JP)=COORDS(1:3*(NATOMS-NSEED),JP)
!             WRITE(MYUNIT,'(A,2G20.10)'),'mc> coordso changed: ',COORDSO(1,JP),COORDS(1,JP)     
!             VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DJW
            ENDIF
            IF (STAY) THEN
            ELSE IF (EVAP .and. .not.evapreject) THEN
               NFAIL(JP)=NFAIL(JP)+1
               CALL MYRESET(JP,NATOMS,NPAR,NSEED)
               IF (DEBUG) THEN
                  WRITE(MYUNIT,33) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
33                FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' EVAP,REJ')
               ENDIF
            ELSE
               ATEST=.TRUE. 

C       csw34> Test for cold fusion in mylbfgs.f - if so, set ATEST to .FALSE.
               IF (COLDFUSION) THEN
                  ATEST=.FALSE.
               ENDIF
               COLDFUSION=.FALSE.
C       csw34> New chirality and peptide bond cis/trans checking block
C              These checks are turned ON my default for both AMBER and
C              CHARMM and can be disabled with ALLOWCISTRANS and 
C              NOCHIRALCHECKS in data. IF using DNA or RNA in AMBER,
C              you should use NOCISTRANSDNA or NOCISTRANSRNA.

C       csw34> First check is for inverted CA or sidechain (CB) centres
               IF (CHECKCHIRALITY.AND.ATEST) THEN
                  IF (CHRMMT) THEN
                     TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
                     CALL CHECKCHIRAL(TEMPCOORDS,CHIRALFAIL)
                  ENDIF
                  IF (AMBERT) THEN 
C       csw34> CHIRALFAIL must be reset, it is not in the AMBER routine
                     CHIRALFAIL=.FALSE.
C       csw34> Call Szilard's routine in amberinterface.f
C       abc> Edit for set_check_chiral, check current chiral indicies against initial
                     IF(SETCHIRAL) THEN
                       CALL set_check_chiral(COORDS(:,JP),NATOMS,GOODSTRUCTURE,chiarray2)               
                       DO J5=1,NATOMS
                         IF ((CHIARRAY1(J5)-CHIARRAY2(J5))/=0) THEN
                           WRITE(MYUNIT,*) ' mc> inverted chirality with respect to initial structure'
                           CHIRALFAIL=.TRUE.
                         ENDIF
                       ENDDO
                     ELSE
                       CALL check_chirality(COORDS(:,JP),NATOMS,GOODSTRUCTURE) 
                       IF (.NOT.(GOODSTRUCTURE)) CHIRALFAIL=.TRUE.
                     ENDIF
                  ENDIF
C       csw34> ATEST is the acceptance test variable. If a centre has
C              been inverted, the step is rejected.
                  IF (CHIRALFAIL) THEN
                     ATEST=.FALSE.
                     WRITE(MYUNIT,*) ' mc> CHIRALITY CHECK FAILED - discarding structure'
                  ENDIF
               ENDIF
C       csw34> As before, but this time checking for trans->cis peptide
C              bonds. The deformation threshold defaults to 150 degrees
C              but can be set with NOCISTRANS
               IF (NOCISTRANS.AND.ATEST) THEN
                  IF (CHRMMT) THEN
                     TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
                     CALL CHECKOMEGA(TEMPCOORDS,AMIDEFAIL)
                  ENDIF
                  IF (AMBERT) THEN
C       csw34> As with CHIRALFAIL, AMIDEFAIL needs resetting
                     AMIDEFAIL=.FALSE.
                     IF (NOCISTRANSDNA) THEN
                        CALL CHECK_CISTRANS_DNA(COORDS(:,JP),NATOMS,ZSYM,GOODSTRUCTURE)
                        IF (.NOT.(GOODSTRUCTURE)) AMIDEFAIL=.TRUE.
                     ELSE IF (NOCISTRANSRNA) THEN
                        CALL CHECK_CISTRANS_RNA(COORDS(:,JP),NATOMS,ZSYM,GOODSTRUCTURE)
                        IF (.NOT.(GOODSTRUCTURE)) AMIDEFAIL=.TRUE.
                     ELSE
                        CALL check_cistrans_protein(COORDS(:,JP),NATOMS,GOODSTRUCTURE,MINOMEGA,cisarray2)
                        DO J5=1,NATOMS
                          IF ((CISARRAY1(J5)-CISARRAY2(J5))/=0) THEN
                             WRITE(*,'(A,I6)') 'cis-trans isomerisation of a peptide bond detected involving atom ', J5
                             AMIDEFAIL=.TRUE.
                          ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
                  IF (AMIDEFAIL) THEN
                     ATEST=.FALSE.
                     WRITE(MYUNIT,*) 'AMIDE BOND CHECK FAILED - discarding structure'
                  ENDIF
               ENDIF
C     csw34> Dumping every quench in AMBER9 format - only dumps if no chiral problems found 
C  Only for AMBER! DJW Otherwise we get the cis-trans message for every other potential!
               IF (AMBERT) THEN
C     csw34> Calculate interaction energy to residue specified in the AMBGMINintE.sh script
                  IF (A9INTET.AND.ATEST) THEN
                     CALL A9INTE(COORDS(:,JP),INTE)
                     CALL A9INTESAVEIT(INTE,COORDS(:,JP),JP)
                     CALL A9INTEDUMPLOWEST()
                     IF (DEBUG) WRITE(MYUNIT,*) 'intE=',INTE
                     IF (TRACKDATAT) THEN
                        WRITE(3998,'(I10,G20.10)') NQ(JP),INTE
                        CALL FLUSH(3998)
                     ENDIF
                  ENDIF
C     csw34> Dump the quench coordinates in AMBER pdb and rst format if DUMPQ is present and chirality checks satisfied
                  IF (DUMPQUT.AND.ATEST) THEN
                      WRITE(QUENCHNUM,*) NQ(JP)
                      QUNAME='quench'//TRIM(ADJUSTL(QUENCHNUM))//'.rst'
                      OPEN(UNIT=20,FILE=QUNAME,STATUS='UNKNOWN')  
                      WRITE(20,'(a20)') QUNAME
                      WRITE(20,'(i5)') NATOMS
                      WRITE(20,'(6f12.7)') COORDS(:,JP) 
                      CLOSE(20)
C     csw34> Dump to PDB using routine in amberinterface.f
                     QUNAME='quench'//TRIM(ADJUSTL(QUENCHNUM))
                     CALL A9DUMPPDB(COORDS(:,JP),QUNAME)
                  ELSEIF (DUMPQUT) THEN
                      WRITE(*,'(A,I6)') 'DUMPQU> chirality/cis-trans problem detected, not dumping quench ',NQ(JP)
                  ENDIF
               ENDIF
C     csw34> Check to see if LOCALSAMPLE constraints have been violated - if they have, reject the step
               IF (LOCALSAMPLET) THEN
                  DISTOK=.FALSE.
                  CALL A9DISTCHECK(COORDS(:,JP),DISTOK)
                  IF (.NOT.DISTOK) ATEST=.FALSE.
               ENDIF
               IF ((QDONE.EQ.0).AND.TIP) THEN
                  ATEST=.FALSE.
               ELSEIF (ATEST) THEN
                  CALL TRANSITION(POTEL,EPREV(JP),ATEST,JP,RANDOM,MCTEMP)
               ENDIF
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
C Accept or reject step. If the quench did not converge then allow a
C potenial move, but count it as a rejection in terms of NSUCCESS and
C NFAIL. This way we will accept a lower minimum if found, but the steps won;t become so big.
C However, for TIP5P some cold fusion events that had not actually reached the threshold for
C rejection were actually accepted. Must prevent this!
               IF (ATEST) THEN
#ifdef MPI
		  IF (DEBUG) THEN
                     WRITE(MYUNIT,334) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
334                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF
#else
                  IF (DEBUG) THEN
                     WRITE(MYUNIT,34) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
34                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF
#endif
                  IF ((J1-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV(JP)).GT.ECONV) THEN
C                    NRELAX=J1-JACCPREV
C                    IF (RESTART) WRITE(MYUNIT,'(A,I6,A)') ' relaxation time set to ',NRELAX,' steps'
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
C           WRITE(MYUNIT,'(A,4F20.10)') 'Q4 values ',QS,QF,QSTART,QFINISH
C
C  Dump coordinates and energy if another run is attempting to jump to this one.
C
            IF ((NPAR.GT.1).AND.(.NOT.NEWJUMP)) CALL DUMPJ(JP,JUMPTO,NPAR,COORDS(1:3*NATOMS,1:NPAR),NATOMS,EPREV)
C
C  If RESTART then reseed if we haven t accepted a step in twice the relaxation time.
C
            IF (RESTART.AND.(J1-JACCPREV.GT.1.1D0*NRELAX)) CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
C
C  Check the acceptance ratio.
C 
            IF ((MOD(J1,NACCEPT).EQ.0).AND.(NSEED.EQ.0).AND.(.NOT.STAY)) CALL ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)

            TEMP(JP)=TEMP(JP)*SCALEFAC
#ifdef MPI
!           CALL MPI_REDUCE(HIT,HITANY,1,MPI_LOGICAL,0,MPI_COMM_WORLD, MPIERR)
!           WRITE(MYUNIT,'(A,I6,2L5)') 'mc> JP,HIT,HITANY=',JP,HIT,HITANY
!           HIT=HITANY
!           IF (HIT) THEN
!              DO J3=0,NPAR-1
!                 IF (J3.EQ.MYNODE) CYCLE
!                 CALL MPI_SEND(HIT,1,MPI_LOGICAL,J3,999,MPI_COMM_WORLD,MPIERR)
!              ENDDO
!              CALL SYSTEM_SUBR('echo hit > hit',ISTAT)
!              GOTO 37
!           ELSE
!
! Calling this MPI routine overwrites memory in both LAM and MPICH. No idea why.
!
!              CALL MPI_IPROBE(MPI_ANY_SOURCE,999,MPI_COMM_WORLD,FLAG,MPIERR)
!              IF (FLAG) THEN
!                 WRITE(MYUNIT,'(A)') 'mc> Target located in some other run'
!                 HIT=.TRUE.
!                 GOTO 37
!              ENDIF
!              INQUIRE(FILE='hit',EXIST=HIT)
!              IF (HIT) THEN
!                 WRITE(MYUNIT,'(A)') 'mc> Target located in some other run'
!                 CALL FLUSH(MYUNIT)
!                 GOTO 37
!              ENDIF
!           ENDIF
#else
            IF (HIT) GOTO 37
#endif
            IF (DUMPINT.GT.0) THEN
               IF (MOD(J1,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(J1,EBEST,BESTCOORDS,JBEST,JP)
               ENDIF
            ENDIF
#ifdef MPI
         E=POTEL
         IF(MYNODE.EQ.0) THEN
            RANDOM=DPRAND()
            J=(NPAR-1)*RANDOM
            RANDOM=DPRAND()
            IF (RANDOM.GT.EXCHPROB) THEN ! 0.1 probability of exchange 
               J=-2
               EXCHANGE=.FALSE.
            ELSE
               EXCHANGE=.TRUE.
            ENDIF
          ENDIF
          CALL MPI_BCAST(J,1,MPI_INTEGER,0,MPI_COMM_WORLD, MPIERR)
          IF (MYNODE.EQ.(J+1)) THEN
             CALL MPI_SEND(E,1,MPI_DOUBLE_PRECISION,J,0,MPI_COMM_WORLD,MPIERR)
          ENDIF
          IF (MYNODE.EQ.J) THEN
             CALL MPI_RECV(ER,1,MPI_DOUBLE_PRECISION,J+1,0,MPI_COMM_WORLD,IS,MPIERR)

             DBETA=BETA(J)-BETA(J+1)
             DELTA=E-ER
             W=MIN(1.0D0,DEXP(DELTA*DBETA))
             NTOT=NTOT+1 
             RANDOM=DPRAND()
             IF (W.GT.RANDOM) THEN
                EXCHANGEACCEPT=.TRUE.
                IMESG=1
                CALL MPI_SEND(IMESG,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(ITRAJ,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,MPIERR)

                CALL MPI_SEND(COORDS(1:3*NATOMS,JP), 3*NATOMS,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(E,1,MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,MPIERR)
                CALL MPI_RECV(ITRAJ,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(COORDS(1:3*NATOMS,JP), 3*NATOMS, MPI_DOUBLE_PRECISION,J+1,1,MPI_COMM_WORLD,IS,MPIERR)
                E=ER
                NEACCEPT=NEACCEPT+1
              ELSE
                EXCHANGEACCEPT=.FALSE.
                IMESG=0
                CALL MPI_SEND(IMESG,1,MPI_INTEGER,J+1,0,MPI_COMM_WORLD,MPIERR)
              ENDIF
          ENDIF
          IF (MYNODE.EQ.(J+1)) THEN
             CALL MPI_RECV(IMESG,1,MPI_INTEGER,J,0,MPI_COMM_WORLD,IS,MPIERR)
             NTOT=NTOT+1
             IF (IMESG.EQ.1) THEN
                CALL MPI_RECV(ITRAJO,1,MPI_INTEGER,J,0,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(COORDSO(1:3*NATOMS,JP),3*NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_RECV(E,1,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,IS,MPIERR)
                CALL MPI_SEND(ITRAJ,1,MPI_INTEGER,J,0,MPI_COMM_WORLD,MPIERR)
                CALL MPI_SEND(COORDS(1:3*NATOMS,JP),3*NATOMS,MPI_DOUBLE_PRECISION,J,1,MPI_COMM_WORLD,MPIERR)
                COORDS(1:3*NATOMS,JP)=COORDSO(1:3*NATOMS,JP)
                ITRAJ=ITRAJO
                NEACCEPT=NEACCEPT+1
              ENDIF
           ENDIF
C
C  When SYMMETRISE is turned on we can have multiple quenches on each pass through the main
C  basin-hopping loop. If we really want to limit the number of quenches to NSTEPS then we
C  need to test NQ(JP) as well.
C
  
           IF (NQ(1).GT.NSTEPS) EXIT
#else
           IF (NQ(JP).GT.NSTEPS) GOTO 37

         ENDDO
#endif
C
C  ****************************** End of loop over NPAR parallel runs *****************************
C
C  Supersteps
C
C        IF (SUPERSTEP) CALL SUPERMC(SPOTEL,SCOORDS,NSUPERCOUNT,POTEL)
C
C  Renormalised type steps
C
C        IF (RENORM) CALL REN(J1,RMIN,RCOORDS,RVAT,NREN,RMINO,RCOORDSO,RVATO,ITERATIONS,TIME,NLAST,JACCPREV,NSTEPREN)
C        IF (STAY) CALL MYRESET(1,NATOMS,NPAR,NSEED)
  
         CALL FLUSH(MYUNIT)
      ENDDO
#ifdef MPI
      WRITE(MYUNIT, '(G20.10,A,G20.10)') NEACCEPT, ' exchanges accepted out of ', NTOT
      IF (MPIT) THEN
! sf344> Apparently MPI_FINALIZE can only be called once during a run, so we call it at the end of main.F
!         CALL MPI_FINALIZE(MPIERR)
!         IF (DEBUG) WRITE(MYUNIT, '(A,I6)') 'In mc after MPI_FINALIZE MPIERR=',MPIERR
      ENDIF
#endif
! }}}

37    CONTINUE
#ifdef MPI
      WRITE(MYUNIT,10) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),STEP(JP),ASTEP(JP),TEMP(JP)
10    FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
      IF (RIGID) WRITE(MYUNIT,'(A,F12.5)') 'Rigid body orientational step=',OSTEP(JP)
#else
      DO JP=1,NPAR
         IF (NPAR.GT.1) THEN
            WRITE(MYUNIT,20) JP,NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),
     1               STEP(JP),ASTEP(JP),TEMP(JP)
20          FORMAT('[',I2,']Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
         ELSE
            WRITE(MYUNIT,21) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),
     1               STEP(JP),ASTEP(JP),TEMP(JP)
21          FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
         ENDIF
         IF (RIGID) WRITE(MYUNIT,'(A,F12.5)') 'Rigid body orientational step=',OSTEP(JP)
      ENDDO
#endif
     
      RETURN
!op226>}}} 
      END
C
C
! csw34> Subroutine to dump the current lowest interaction energy structure 
!> \brief Subroutine to dump the current lowest interaction energy structure 
!> \author Chris Whittleston, csw34@cam.ac.uk
!
      SUBROUTINE A9INTEDUMPLOWEST()
      USE COMMONS
      USE QMODULE
      IMPLICIT NONE
      OPEN(UNIT=20,FILE="bestintE.rst",STATUS='UNKNOWN')  
      WRITE(20,'(g20.10,i5)') INTEQMIN(1), INTEFF(1)
      WRITE(20,'(i5)') NATOMS
      WRITE(20,'(6f12.7)') INTEQMINP(1,:) 
      CLOSE(20)
C     csw34> Dump to PDB using routine in amberinterface.f
      CALL A9DUMPPDB(INTEQMINP(1,:),"bestintE")
      END SUBROUTINE A9INTEDUMPLOWEST
! csw34> Subroutine to work out the interaction energy INTE between residue RESLIG for geometry INTECOORDS
      SUBROUTINE A9INTE(INTECOORDS,INTE)
      USE PORFUNCS
      USE COMMONS, ONLY : NATOMS
      IMPLICIT NONE
      DOUBLE PRECISION :: INTECOORDS(3*NATOMS), INTE
      INTEGER ISTAT
! Dump current coordinates to file for SANDER to read in
      OPEN(UNIT=20,FILE='coords.intres',STATUS='UNKNOWN')
      WRITE(20,'(a20)') 'intE coordinates'
      WRITE(20,'(i5)') NATOMS
      WRITE(20,'(6f12.7)') INTECOORDS(:)
      CLOSE(20)
! Run the script that does the interaction energy calculation (downloadable from the group website)
      CALL SYSTEM_SUBR('bash AMBGMINintE.sh',ISTAT)
! Read interaction energy in from the file called intE
      OPEN(UNIT=20,FILE="intE",STATUS='OLD')
      READ(20,*) INTE
      CLOSE(20)
      END SUBROUTINE A9INTE

!
!> \brief Subroutine to check the  distance 
!> \brief between groups of atoms defined in the movableatoms file.
!> Checks the A->B and A->C distances,
!> and if they are greater than the A(B/C)THRESH values defined
!> in the data file, the routine returns DISTOK=.FALSE.
!> \author Chris Whittleston, csw34@cam.ac.uk
!
      SUBROUTINE A9DISTCHECK(COORDS,DISTOK)
      USE MODAMBER9, ONLY : NATOMSINA,NATOMSINB,NATOMSINC,ATOMSINALIST,ATOMSINBLIST,ATOMSINCLIST
      USE COMMONS, ONLY : NATOMS,ABTHRESH,ACTHRESH,DEBUG
      IMPLICIT NONE
      LOGICAL :: DISTOK
      INTEGER :: I,J
      DOUBLE PRECISION :: COORDS(3*NATOMS)
      DOUBLE PRECISION :: CENTREOFA(3),CENTREOFB(3),CENTREOFC(3)  
      DOUBLE PRECISION :: ABDIST,ACDIST 
! initialise variables
      CENTREOFA(:)=0.0D0
      CENTREOFB(:)=0.0D0
      CENTREOFC(:)=0.0D0
      ABDIST=0.0D0
      ACDIST=0.0D0
! find centre of ligand (group A)
      DO I=1,NATOMSINA
         J=ATOMSINALIST(I)
         CENTREOFA(1)=CENTREOFA(1)+COORDS(3*J-2)
         CENTREOFA(2)=CENTREOFA(2)+COORDS(3*J-1)
         CENTREOFA(3)=CENTREOFA(3)+COORDS(3*J  )
      END DO
      CENTREOFA(:) = CENTREOFA(:)/NATOMSINA
! find centre of group B 
      DO I=1,NATOMSINB
         J=ATOMSINBLIST(I)
         CENTREOFB(1)=CENTREOFB(1)+COORDS(3*J-2)
         CENTREOFB(2)=CENTREOFB(2)+COORDS(3*J-1)
         CENTREOFB(3)=CENTREOFB(3)+COORDS(3*J  )
      END DO
      CENTREOFB(:) = CENTREOFB(:)/NATOMSINB
! find centre of group B 
      DO I=1,NATOMSINC
         J=ATOMSINCLIST(I)
         CENTREOFC(1)=CENTREOFC(1)+COORDS(3*J-2)
         CENTREOFC(2)=CENTREOFC(2)+COORDS(3*J-1)
         CENTREOFC(3)=CENTREOFC(3)+COORDS(3*J  )
      END DO
      CENTREOFC(:) = CENTREOFC(:)/NATOMSINC
! calculate A->B distance
      ABDIST=SQRT((CENTREOFA(1)-CENTREOFB(1))**2+(CENTREOFA(2)-CENTREOFB(2))**2+(CENTREOFA(3)-CENTREOFB(3))**2)
! calculate A->C distance
      ACDIST=SQRT((CENTREOFA(1)-CENTREOFC(1))**2+(CENTREOFA(2)-CENTREOFC(2))**2+(CENTREOFA(3)-CENTREOFC(3))**2)
! some DEBUG printing
      IF (DEBUG) THEN
         WRITE(*,*) 'AB distance=',ABDIST
         WRITE(*,*) 'AC distance=',ACDIST
         IF (ABDIST.LT.ABTHRESH) THEN
            WRITE(*,*) 'A->B condition met! :)'
         ELSE 
            WRITE(*,*) 'A->B condition broken :('
         ENDIF   
         IF (ACDIST.LT.ACTHRESH) THEN
            WRITE(*,*) 'A->C condition met! :)'
         ELSE 
            WRITE(*,*) 'A->C condition broken :('
         ENDIF   
      ENDIF
! do the check for both conditions
      IF ((ABDIST.LT.ABTHRESH).AND.(ACDIST.LT.ACTHRESH)) DISTOK=.TRUE.
! more debug printing
      IF (DEBUG) WRITE(*,*) 'DISTOK=',DISTOK
 
      END SUBROUTINE A9DISTCHECK


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
C     COMMON /IG/ IGNOREBIN, FIXBIN
      common /ev/ evap, evapreject

      IF (DISTMINOLD.EQ.0.0D0) DISTMINOLD=DISTMIN  ! this should allow for the first step
      IF (TUNNELT) THEN
         TEOLD=TRANS(EOLD,QMIN(1),GAMMA)
         TENEW=TRANS(ENEW,QMIN(1),GAMMA)
C        WRITE(MYUNIT,'(A,4F20.10)') 'TEOLD,TENEW,QMIN(1),GAMMA=',TEOLD,TENEW,QMIN(1),GAMMA
      ELSE
         TEOLD=EOLD
         TENEW=ENEW
      ENDIF

      IF (TSALLIST) THEN
         EREF=QMIN(NP)*1.1D0
         DUMMY=(1.0D0-(1.0D0-QTSALLIS)*(TENEW-EREF)/MCTEMP)/(1.0D0-(1.0D0-QTSALLIS)*(TEOLD-EREF)/MCTEMP)
         DUMMY=DUMMY**(QTSALLIS/(1.0D0-QTSALLIS))
C        WRITE(MYUNIT,'(A,4F20.10)') 'TENEW,TEOLD,EREF,DUMMY=',TENEW,TEOLD,EREF,DUMMY
         IF (DUMMY.GE.1.0D0) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DUMMY.GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF
      ELSE
C
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
      ENDIF 

      RETURN 
      END 

      SUBROUTINE JUMPM(RANNJ,J1,JP,EPPREV)
      USE commons
      IMPLICIT NONE
      INTEGER J1, JP, J2, UNT, NDUM, ITERATIONS, BRUN, QDONE
      DOUBLE PRECISION RANNJ, RANDOM, DPRAND, EPPREV(NPAR), DUMMY, TIME, EJ, SCREENC(3*NATOMS)

      IF (NEWJUMP) THEN
         IF (RANNJ.LT.PNEWJUMP) THEN
            RANDOM=DPRAND()
            IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
C           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))/TEMP(JP)).GT.RANDOM) THEN
C           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
               WRITE(MYUNIT,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1                ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' accepted before quench ',J1
               DUMMY=EPREV(JP)
               EPREV(JP)=EPREV(JUMPTO(JP))
               EPREV(JUMPTO(JP))=DUMMY
               DUMMY=EPPREV(JP)
               EPPREV(JP)=EPPREV(JUMPTO(JP))
               EPPREV(JUMPTO(JP))=DUMMY
               DO J2=1,NATOMS
                  DUMMY=VATO(J2,JP)
                  VATO(J2,JP)=VATO(J2,JUMPTO(JP))
                  VATO(J2,JUMPTO(JP))=DUMMY
                  DUMMY=VAT(J2,JP)
                  VAT(J2,JP)=VAT(J2,JUMPTO(JP))
                  VAT(J2,JUMPTO(JP))=DUMMY
               ENDDO
               DO J2=1,3*NATOMS
                  DUMMY=COORDS(J2,JP)
                  COORDS(J2,JP)=COORDS(J2,JUMPTO(JP))
                  COORDS(J2,JUMPTO(JP))=DUMMY
                  DUMMY=COORDSO(J2,JP)
                  COORDSO(J2,JP)=COORDSO(J2,JUMPTO(JP))
                  COORDSO(J2,JUMPTO(JP))=DUMMY
               ENDDO
            ELSE
               WRITE(MYUNIT,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1                  ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' rejected before quench ',J1
            ENDIF
         ENDIF
      ELSE
         UNT=70+JUMPTO(JP)
         REWIND(UNT)           
         NDUM=INT(DPRAND()*(NQ(JUMPTO(JP))-1))
         IF (DEBUG) WRITE(MYUNIT,'(A, G20.10)') 'Should be choosing buffer energy number ',NDUM
         DO J2=1,NDUM
            READ(UNT,*) EJ
         ENDDO
         DO J2=1,NQ(JUMPTO(JP))-1-NDUM
            READ(UNT,*) DUMMY
         ENDDO
         RANDOM=DPRAND()
C
C  Coordinates are only read if the jump is successful.
C
         IF (DEXP((EPREV(JP)-EJ)*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
            WRITE(MYUNIT,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' accepted before quench ',NQ(JP)
            EPREV(JP)=EJ
            UNT=70+NPAR+JUMPTO(JP)
            REWIND(UNT)           
            DO J2=1,(NDUM-1)*NATOMS
               READ(UNT,*) DUMMY
            ENDDO
            READ(UNT,*) (COORDS(J2,JP),J2=1,3*NATOMS)
C
C  Coordinates should be converged already, but we need to reset VAT and VATO.
C
C  next line should be uncommented if routine is made availabe to use with CHARMM
C            IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
            CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            DO J2=1,NATOMS
               VATO(J2,JP)=VAT(J2,JP)
            ENDDO
            DO J2=1,3*NATOMS
               COORDSO(J2,JP)=COORDS(J2,JP)
            ENDDO
            IF (DEBUG) THEN
               WRITE(MYUNIT,'(A)') 'Jump coordinates:'
               WRITE(MYUNIT,'(3F20.10)') (COORDS(J2,JP),J2=1,3*NATOMS)
            ENDIF
            DO J2=1,NATOMS*(NQ(JUMPTO(JP))-1)-NDUM*NATOMS
               READ(UNT,*) DUMMY
            ENDDO
         ELSE
            WRITE(MYUNIT,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,
     1              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' rejected before quench ',NQ(JP)
         ENDIF
      ENDIF

      RETURN
      END

      SUBROUTINE DUMPJ(JP,JUMPTO,NPAR,COORDS,NATOMS,EPREV)
      IMPLICIT NONE
      LOGICAL TEST
      INTEGER NPAR, J2, JP, JUMPTO(NPAR), NATOMS, UNT
      DOUBLE PRECISION COORDS(3*NATOMS,NPAR), EPREV(NPAR)

      DO J2=1,NPAR
         IF (JUMPTO(J2).EQ.JP) TEST=.TRUE.
      ENDDO
C     TEST=.FALSE.
C     IF (TEST) THEN
C        UNT=70+JP
C        WRITE(UNT,'(F20.10)') EPREV(JP)
C        UNT=70+NPAR+JP
C        WRITE(UNT,'(3F20.10)') (COORDS(J2,JP),J2=1,3*NATOMS)
C     ENDIF

      RETURN
      END

      SUBROUTINE MYRESET(JP,NATOMS,NPAR,NSEED)
      USE COMMONS,ONLY : MYUNIT,COORDS,COORDSO,VAT,VATO
      IMPLICIT NONE
      INTEGER JP, NSEED, J2, NATOMS, NPAR

      DO J2=1,3*(NATOMS-NSEED)
         COORDS(J2,JP)=COORDSO(J2,JP)
      ENDDO
      DO J2=1,NATOMS
         VAT(J2,JP)=VATO(J2,JP)
      ENDDO

      RETURN
      END
 
 
 
      SUBROUTINE REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
      USE commons
      IMPLICIT NONE
      INTEGER ITERATIONS, J2, JACCPREV, J1, NQTOT, BRUN, QDONE
      DOUBLE PRECISION TIME, POTEL, RCOORDS(3*NATOMS), RMIN, RVAT(NATOMS), SCREENC(3*NATOMS)
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT

10    CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),1,NHSRESTART)
C  next line should be uncommented if routine is made availabe to use with CHARMM
C      IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
      CALL QUENCH(.FALSE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
C
C  Bad idea to accept this quench configuration unconditionally - it could be unconvergeable.
C
      IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
         DO J2=1,3*NATOMS
            COORDS(J2,1)=COORDSO(J2,1)
         ENDDO
         GOTO 10
      ENDIF
      JACCPREV=J1
      NQTOT=NQTOT+1
      WRITE(MYUNIT,'(A,I6,A)') ' Restarting using ',NHSRESTART,' hard sphere moves'
      WRITE(MYUNIT,'(A,I7,A,F20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Restart Qu ',NQ(1),' E=',
     1              POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
      DO J2=1,3*NATOMS
         COORDSO(J2,1)=COORDS(J2,1)
         RCOORDS(J2)=COORDS(J2,1)
      ENDDO
      DO J2=1,NATOMS
         VATO(J2,1)=VAT(J2,1)
         RVAT(J2)=VAT(J2,1)
      ENDDO
      EPREV(1)=POTEL
      RMIN=POTEL

      RETURN
      END
 
 
 
      SUBROUTINE ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
      USE commons
      USE modcharmm
      IMPLICIT NONE
      INTEGER NSUCCESS(NPAR), NFAIL(NPAR), JP, NFAILT(NPAR), NSUCCESST(NPAR), J1, J2, NDUMMY
      LOGICAL evap, evapreject
      DOUBLE PRECISION DUMMY, DUMMY2, DUMMY3, DUMMY4, HWMAX,P0,FAC
C     COMMON /IG/ IGNOREBIN, FIXBIN
C     COMMON /MOVE/ TMOVE, OMOVE
      common /ev/ evap, evapreject

      P0=1.D0*NSUCCESS(JP)/(1.D0*(NSUCCESS(JP)+NFAIL(JP)))
      
      IF (P0.GT.ACCRAT(JP)) THEN
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
                  IF(CHRIGIDTRANST.AND.CHRMMT) TRANSMAX=FAC*TRANSMAX
                  IF(CHRIGIDROTT.AND.CHRMMT) ROTMAX=FAC*ROTMAX  
               ENDIF
            ENDIF
            ASTEP(JP)=ASTEP(JP)*1.05D0
         ENDIF
      ELSE
         IF(ARMT) THEN
           FAC=LOG(ARMA*ACCRAT(JP)+ARMB)/LOG(ARMA*P0+ARMB)
         ELSE
           FAC=1.D0/1.05D0
         ENDIF
         IF (FIXBOTH(JP)) THEN
         ELSE IF (FIXSTEP(JP)) THEN
            IF (.NOT.FIXTEMP(JP)) TEMP(JP)=TEMP(JP)*1.05D0
         ELSE
            IF (FIXD) THEN
               NHSMOVE=MAX(1,NHSMOVE-1)
            ELSE
               IF (RIGID) THEN
                  IF (TMOVE(JP)) STEP(JP)=STEP(JP)/1.05D0
                  IF (OMOVE(JP)) OSTEP(JP)=OSTEP(JP)/1.05D0
               ELSE
                  STEP(JP)=FAC*STEP(JP)
                  IF(CHRIGIDTRANST.AND.CHRMMT) TRANSMAX=FAC*TRANSMAX
                  IF(CHRIGIDROTT.AND.CHRMMT) ROTMAX=FAC*ROTMAX
               ENDIF
            ENDIF
            ASTEP(JP)=ASTEP(JP)/1.05D0
         ENDIF
      ENDIF
C
C Prevent steps from growing out of bounds. The value of 1000 seems sensible, until
C we do something with such huge dimensions?!
c
      STEP(JP)=MIN(STEP(JP),1.0D3)
      OSTEP(JP)=MIN(OSTEP(JP),1.0D3)
      ASTEP(JP)=MIN(ASTEP(JP),1.0D3)
C
      IF (NPAR.GT.1) THEN
         WRITE(MYUNIT,'(A,I2,A,I6,A,F8.4,A,F8.4)') '[',JP,']Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
      ELSE
         WRITE(MYUNIT,'(A,I6,A,F8.4,A,F8.4)') 'Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
      ENDIF
      IF (FIXBOTH(JP)) THEN
      ELSE IF (FIXSTEP(JP)) THEN
         IF(.NOT.FIXTEMP(JP)) WRITE(MYUNIT,'(A,F12.4)') 'Temperature is now:',TEMP(JP)
      ELSE
         IF (NPAR.GT.1) THEN
            WRITE(MYUNIT,'(A,I2,A)',ADVANCE='NO') '[',JP,']Steps are now:'
         ELSE
            WRITE(MYUNIT,'(A)',ADVANCE='NO') 'Steps are now:'
         ENDIF
         WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO') '  STEP=',STEP(JP)    
         IF(ASTEP(JP).GT.0.D0) WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO')'  ASTEP=',ASTEP(JP) 
         IF(CHRIGIDTRANST.AND.CHRMMT) WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO')'  TRANSMAX=',TRANSMAX
         IF(CHRIGIDROTT.AND.CHRMMT) WRITE(MYUNIT,'(A,F10.4)')'  ROTMAX=',ROTMAX
         IF(.NOT.FIXTEMP(JP)) WRITE(MYUNIT,'(A,F10.4)') ' Temperature is now:',TEMP(JP)
         IF (RIGID) WRITE(MYUNIT,'(A,F12.6,A,F12.6)') 'Maximum rigid body rotational move is now ',OSTEP(JP)
      ENDIF
      IF (FIXD) WRITE(MYUNIT,'(A,I4)') 'hard sphere collision moves=',NHSMOVE
C
      NSUCCESST(JP)=NSUCCESST(JP)+NSUCCESS(JP)
      NFAILT(JP)=NFAILT(JP)+NFAIL(JP)
      NSUCCESS(JP)=0
      NFAIL(JP)=0 
C
      RETURN
      END

C
C
      SUBROUTINE REN(J1,RMIN,RCOORDS,RVAT,NREN,RMINO,RCOORDSO,RVATO,ITERATIONS,TIME,NLAST,JACCPREV,NSTEPREN)
      USE commons
      IMPLICIT NONE
      LOGICAL STAY, REJECT, METROPOLIS
      INTEGER J1, J2, NREN, ITERATIONS, NQTOT, NLAST, JACCPREV, NSTEPREN, J3, BRUN, QDONE
      DOUBLE PRECISION POTEL, RMIN, RCOORDS(3*NATOMS), RVAT(NATOMS), RANDOM, DPRAND, RMINO, RCOORDSO(3*NATOMS), RVATO(NATOMS),
     1                 TIME, XIP, DUMMY, SCREENC(3*NATOMS)
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT

      STAY=.FALSE.
      IF (POTEL.LT.RMIN) THEN
         RMIN=POTEL          
         DO J2=1,3*NATOMS
            RCOORDS(J2)=COORDS(J2,1)
         ENDDO
         DO J2=1,NATOMS
            RVAT(J2)=VAT(J2,1)
         ENDDO
      ENDIF
      IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') 'RMIN,POTEL=',RMIN,POTEL
C     PRINT*,'J1,JACCPREV+NRENSTUCK,NLAST+NREN=',J1,JACCPREV+NRENSTUCK,NLAST+NREN
      IF ((J1.GE.JACCPREV+NRENSTUCK).OR.(J1.GE.NLAST+NREN)) THEN
         JACCPREV=J1
         NLAST=J1
         RANDOM=DPRAND()
         METROPOLIS=DEXP(-(RMIN-RMINO)/TRENORM).GT.RANDOM
         REJECT=.FALSE.
C
C  Taboo list for renormalised energies. Skip if the step is going to be rejected by
C  the Metropolis condition.
C
         IF (TABOOT.AND.METROPOLIS) THEN
            IF (NSTEPREN.EQ.0) NT(1)=0
            CALL NEWINERTIA(RCOORDS,NATOMS,NATOMS,XIP)
            DO J1=1,NT(1)
               IF (DABS(RMIN-ESAVE(J1,1)).LT.ECONV) THEN
                  IF (DABS(XIP-XINSAVE(J1,1))/(XIP+XINSAVE(J1,1)).LT.1.0D-2) THEN
                     REJECT=.TRUE.
                     GOTO 20
                  ELSE
                     WRITE(MYUNIT,'(A, 2G20.10)') 'Energies nearly degenerate:',RMIN,ESAVE(J1,1)
                     WRITE(MYUNIT,'(A, 2G20.10)') 'But  different  structures:',XIP,XINSAVE(J1,1)
                  ENDIF
               ENDIF
               IF (RMIN.LT.ESAVE(J1,1)) THEN
                  NT(1)=MIN(NT(1)+1,NTAB)
                  DO J3=NT(1),J1+1,-1
                     ESAVE(J3,1)=ESAVE(J3-1,1)
                     XINSAVE(J3,1)=XINSAVE(J3-1,1)
                  ENDDO
                  ESAVE(J1,1)=RMIN
                  XINSAVE(J1,1)=XIP
                  GOTO 20
               ENDIF
            ENDDO
            
            NT(1)=MIN(NT(1)+1,NTAB)
            ESAVE(NT(1),1)=RMIN
            XINSAVE(NT(1),1)=XIP

20          CONTINUE

            WRITE(MYUNIT,'(A,I10)') ' Number of entries in taboo list=',NT(1)
            IF (DEBUG) THEN
               WRITE(MYUNIT,'(6F20.10)') (ESAVE(J2,1),J2=1,NT(1))
            ENDIF
         ENDIF
C
C  Accept/reject for renormalised step
C
         IF (METROPOLIS.AND.(.NOT.REJECT)) THEN
            IF (NSTEPREN.GT.0) WRITE(MYUNIT,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' accepted'
            NREN=MAX(NREN/1.1D0,NRENORM/2.0D0)
            RMINO=RMIN
            IF (.NOT.STAY) THEN
               DO J2=1,3*NATOMS
                  RCOORDSO(J2)=RCOORDS(J2)
                  COORDS(J2,1)=RCOORDS(J2)
               ENDDO
               DO J2=1,NATOMS
                  RVATO(J2)=RVAT(J2)
               ENDDO
            ELSE
               DO J2=1,3*NATOMS
                  RCOORDSO(J2)=COORDSO(J2,1)
                  COORDS(J2,1)=COORDSO(J2,1)
               ENDDO
            ENDIF
         ELSE
            IF (REJECT) THEN
               WRITE(MYUNIT,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected by taboo criterion'
            ELSE
               WRITE(MYUNIT,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected'
            ENDIF
            DO J2=1,3*NATOMS
               COORDS(J2,1)=RCOORDSO(J2)
            ENDDO
            DO J2=1,NATOMS
               RVAT(J2)=RVATO(J2)
            ENDDO
            NREN=NREN*1.1D0
         ENDIF
         NSTEPREN=NSTEPREN+1
         IF (NSTEPREN.EQ.1) WRITE(MYUNIT,'(A,G20.10)') ' first renorm energy is ',RMIN
         WRITE(MYUNIT,'(A,I6)') ' renormalisation interval is now ',NREN
C
C  Longer renorm step
C
10       IF ((XMOVERENORM.GT.3.0D0).OR.FIXD) THEN
            CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),1,INT(XMOVERENORM))
         ELSE
            DUMMY=STEP(1)
            STEP(1)=XMOVERENORM
            CALL TAKESTEP(1)
            STEP(1)=DUMMY
         ENDIF
C  next line should be uncommented if routine is made availabe to use with CHARMM
C         IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
         CALL QUENCH(.FALSE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
C
C  Bad idea to accept this quench configuration unconditionally - it could be unconvergeable.
C
         IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,1)=COORDSO(J2,1)
            ENDDO
            GOTO 10
         ENDIF
         NQTOT=NQTOT+1
         WRITE(MYUNIT,'(A,I7,A,F20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Renorm Qu ',NQ(1),' E=',
     1        POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
         DO J2=1,3*NATOMS
            COORDSO(J2,1)=COORDS(J2,1)
            RCOORDS(J2)=COORDS(J2,1)
         ENDDO
         DO J2=1,NATOMS
            VATO(J2,1)=VAT(J2,1)
            RVAT(J2)=VAT(J2,1)
         ENDDO
         EPREV(1)=POTEL
         RMIN=POTEL
      ENDIF

      RETURN
      END
C
C  Reseed if the energy has not improved by more than ECONV over the
C  last NRELAX mc steps.
C  If AVOID is true then save the energy and coordinates of the lowest
C  minimum achieved before each restart and restart if we get too close
C  to any one of them using bipartite matching and mind. Note that bipartite
C  matching and mind can give a local minimum of distance if the optimal 
C  permutation-inversion isn;t found. Using ORIENT to effect a standard
C  orientation first seems to help. It should at least ensure that permutation-inversion
C  isomers are always found. 
C  Would it be possible to just use the energy and inertia components to identify
C  permutation-inversion isomers instead of MINPERM and MIND?
C  This would be like using a small value of AVOIDDIST, which doesn;t seem to be
C  as good.
C
      SUBROUTINE NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,
     1                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST)
      USE commons
      USE modamber9, only : mdstept 
      
      IMPLICIT NONE
      INTEGER J1, JP, JBEST(NPAR), ITERATIONS, J2, JACCPREV, BRUN, QDONE, J3, PERM(NATOMS), NPERM, NTRIES
      INTEGER, PARAMETER :: MAXIMUMTRIES=20
      DOUBLE PRECISION EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR), EPPREV(NPAR), POTEL, TIME, RCOORDS(3*NATOMS), DIST2, DUMMY2,
     1                 RVAT(NATOMS), RMIN, RANDOM, SR3, SCREENC(3*NATOMS), DPRAND, FCOORDS(3*NATOMS), XMSBSAVE(3*NATOMS),
     2                 DUMMY(3*NATOMS), DISTANCE, XMSB(3*NATOMS), EBESTP, BESTCOORDSP(3*NATOMS), WORSTRAD, RMAT(3,3), QENERGY,
     &                 ROTA(3,3), ROTINVA(3,3)
      INTEGER NSUCCESS(NPAR), NFAIL(NPAR), NFAILT(NPAR), NSUCCESST(NPAR), NORBIT1, NORBIT2, INVERT, NPOSITION
      LOGICAL RES1, RES2, HIGHEST
      COMMON /MYPOT/ QENERGY

      SR3=DSQRT(3.0D0)
      IF (POTEL.LT.EBEST(JP)) THEN
         IF (EBEST(JP)-POTEL.GT.ECONV) JBEST(JP)=J1
         EBESTP=EBEST(JP)
         BESTCOORDSP(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP) ! save previous BESTCOORDS for possible use
         EBEST(JP)=POTEL ! reset ebest, but not necessarily jbest
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      ENDIF
C
C  Reseed if the energy has not improved in the last NRELAX mc cycles,
C  or if the current minimum is too close to one of the NMSBSAVE structures
C  saved in MSBCOORDS.
C
C  Instead of using the current minimum, employ the current minimum in 
C  BESTCOORDS for the AVOID check. Then we need only do the AVOID check when we have
C  a new minimum in BESTCOORDS, i.e. J1.EQ.JBEST(JP).
C
      RES1=.FALSE.
      IF (J1-JBEST(JP).GT.NRELAX) RES1=.TRUE.
!     WRITE(MYUNIT,'(A,I5,2G17.7,3I5,L8)') 'J1,POTEL,EBEST,JBEST,J1-JBEST,NRELAX,RES1=',
!    1                                J1,POTEL,EBEST(JP),JBEST(JP),J1-JBEST(JP),NRELAX,RES1
      RES2=.FALSE.
C     IF ((.NOT.RES1).AND.AVOID) THEN
!     PRINT*,'RES1,AVOID,J1,JBEST(JP)=',RES1,AVOID,J1,JBEST(JP)
      IF ((.NOT.RES1).AND.AVOID.AND.(J1.EQ.JBEST(JP)).AND.(NMSBSAVE.GT.0)) THEN ! best minimum has just changed.
         FCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
         savedloop: DO J2=1,MIN(NMSBSAVE,MAXSAVE)
C
C  Bipartite matching routine for permutations. Coordinates in FCOORDS do not change
C  but the coordinates in XMSB do. DISTANCE is the distance in this case.
C
! If the energy is lower no reseeding regardless of separation ? DJW
!           PRINT '(A,2G20.10,L8)','POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV=',
!    1                              POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV
            IF (POTEL.LT.MSBE(J2)-ECONV) CYCLE savedloop
            XMSB(1:3*NATOMS)=MSBCOORDS(1:3*NATOMS,J2)
            CALL MINPERMDIST(FCOORDS,XMSB,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DISTANCE,DIST2,RIGID,RMAT)
            IF (DISTANCE.LT.AVOIDDIST) THEN
               RES2=.TRUE.
               WRITE(MYUNIT,'(A,I6,A,G20.10,A,F10.3)') 'newres> Current minimum is too close to saved structure ',
     1                                 J2,' with energy ',MSBE(J2),' dist=',DISTANCE
               GOTO 20
            ENDIF
         ENDDO savedloop
      ENDIF

20    CONTINUE 


      IF (RES1.OR.RES2) THEN
          !  Does not seem necessary to save MSB data for restarts.
          !  If we are reseeding because RES2 is true then:
          !  (1) If the AVOID test is based on the distance for BESTCOORDS, add
          !      these to the  MSB list if the distance is > 0.01D0, i.e. a different structure.
          !      Otherwise, add the coordinates of the previous lowest energy minimum to
          !      the MSB list; these have been saved in BESTCOORDSP and the energy in EBESTP.
          !  (2) Alternatively, if the AVOID test is based on the coordinates of the current
          !      structure, we should add these to the MSB list instead. Testing
          !      every step involves an overhead, which grows with the number of saved taboo
          !      structures. We could make the taboo list cyclic, so that new structures to be
          !      avoided replace old ones at the beginning of the list when the length of the
          !      list is exceeded. We are currently only doing the bipartite matching test if
          !      the energy of the lowest minimum since the last reseeding has changed, which
          !      saves time.
          !
          !  Change to a cyclic AVOID list 30/12/04.

         IF (RES1.OR.(RES2.AND.(DISTANCE.GT.0.01D0))) THEN ! new condition
            NMSBSAVE=NMSBSAVE+1
            NPOSITION=MOD(NMSBSAVE,MAXSAVE)
            IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
            MSBE(NPOSITION)=EBEST(JP)
            FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP)

            CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)

            MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)

            WRITE(MYUNIT,'(A,I6,A,G20.10)') 'newres> Moving current best minimum to position ',NPOSITION,
     &                         ' in the AVOID list E=',EBEST(JP)
C           OPEN(UNIT=34,FILE='MSBdata',POSITION='APPEND')
C           WRITE(34,'(G20.10)') MSBE(NPOSITION) 
C           WRITE(34,'(3G20.10)') BESTCOORDS(1:3*NATOMS,JP)
C           CLOSE(34)
         ELSEIF (RES2.AND.(DISTANCE.LE.0.01D0)) THEN
C           IF (NMSBSAVE.LT.MAXSAVE) THEN
               NMSBSAVE=NMSBSAVE+1
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
               MSBE(NPOSITION)=EBESTP
               FCOORDS(1:3*NATOMS)=BESTCOORDSP(1:3*NATOMS)
               CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)
               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
               WRITE(MYUNIT,'(A,I6,A,G20.10)') 'newres> Moving previous best minimum to position ',NPOSITION,
     1                             ' in the AVOID list E=',EBESTP
!           ENDIF
         ENDIF ! end new condition

         IF (NEWRESTART_MD) THEN ! lb415
            IF (RES1) WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
            WRITE(MYUNIT,'(A,I8,A)') 'newres> Reseeding due to a short high temperature MD run'
            CHANGE_TEMP = .true.

! th368: 20-10-2009 Extending MD-Reseeding to the CHARMM interface
! terminate if neither AMBER nor CHARMM was requested in the data file

             IF (AMBERT) THEN
               CALL TAKESTEPAMBER(JP,COORDS(:,JP),movableatomlist,nmovableatoms,ligmovet,mdstept,randomseedt)
             ELSEIF (CHRMMT) THEN
               CALL CHMD(JP)
             ELSE
               WRITE(MYUNIT,'(A,I8,A)') 'newres> Molecular Dynamics reseeding is available for AMBER or CHARMM runs only.'
               STOP
             ENDIF
! end th368: 20-10-2009

            CHANGE_TEMP = .false.
         ELSE
            IF (NHSRESTART.GT.0) THEN
               IF (RES1) WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
               IF (RES2) WRITE(MYUNIT,'(A,I8,A)') 'newres> Reseeding due to taboo condition'
               IF (SHELLMOVES(JP)) WRITE(MYUNIT,'(A)') 'newres> Turning off shell moves'
               SHELLMOVES(JP)=.FALSE.
               CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
            ELSE
               IF (RES1) WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - reseeding'
               IF (RES2) WRITE(MYUNIT,'(A,I8,A)') 'newres> Reseeding due to taboo condition'
               HIGHEST=.TRUE.
               DO J2=1,NPAR
                  IF (J2.EQ.JP) CYCLE
                  IF (EBEST(J2).GT.EBEST(JP)) THEN
                     HIGHEST=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF (HIGHEST) THEN
                  WRITE(MYUNIT,'(A,I6,A,F20.10)') 'newres> Parallel run ',JP,' has the highest energy ',EBEST(JP)
                  IF (NPAR.GT.1) WRITE(MYUNIT,'(6F20.10)') EBEST(1:NPAR)
                  WRITE(MYUNIT,'(A)') 'newres> Full reseeding'
                  DO J2=1,3*NATOMS
                     RANDOM=(DPRAND()-0.5D0)*2.0D0
                     COORDS(J2,JP)=RANDOM*DSQRT(RADIUS)/SR3
                  ENDDO
                  NCORE(JP)=0
                  PTGROUP(JP)='   '
               ELSEIF (NCORE(JP).GT.0) THEN
C         
C  Reseed everything if this is the lowest of a set of parallel runs, otherwise
C  just reseed the surface.
C         
                  WRITE(MYUNIT,'(A)') 'newres> Accepting an enforced surface reseeding'
                  DUMMY2=-1.0D0
                  DO J2=NATOMS-NCORE(JP)+1, NATOMS
                     DISTANCE=COORDS(3*(J2-1)+1,JP)**2+COORDS(3*(J2-1)+2,JP)**2+COORDS(3*(J2-1)+3,JP)**2
                     IF (DISTANCE.GT.DUMMY2) DUMMY2=DISTANCE
                  ENDDO
                  DISTANCE=SQRT(DISTANCE)
                  WRITE(MYUNIT,'(A,F15.5)') 'newres> largest core radius=',DISTANCE
                  DO J2=1,NATOMS-NCORE(JP)
                     COORDS(3*(J2-1)+1,JP)=(DPRAND()-0.5D0)*2.0D0
                     COORDS(3*(J2-1)+2,JP)=(DPRAND()-0.5D0)*2.0D0
                     COORDS(3*(J2-1)+3,JP)=(DPRAND()-0.5D0)*2.0D0
                     RANDOM=DISTANCE+1.0D0+DPRAND()*(SQRT(RADIUS)-DISTANCE-1.0D0)
                     DUMMY2=SQRT(COORDS(3*(J2-1)+1,JP)**2+COORDS(3*(J2-1)+2,JP)**2+COORDS(3*(J2-1)+3,JP)**2)
                     COORDS(3*(J2-1)+1,JP)=COORDS(3*(J2-1)+1,JP)*RANDOM/DUMMY2
                     COORDS(3*(J2-1)+2,JP)=COORDS(3*(J2-1)+2,JP)*RANDOM/DUMMY2
                     COORDS(3*(J2-1)+3,JP)=COORDS(3*(J2-1)+3,JP)*RANDOM/DUMMY2
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
C
C  This step will be accepted because EPREV(JP)=POTEL. However, there
C  could be additional steps, so we need
C  to change COORDSO and VATO. Should reset EBEST and BESTCOORDS as well.
C
C  next line should be uncommented if routine is made availabe to use with CHARMM
C         IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
654      CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         IF (RMS.GT.BQMAX) THEN
            WRITE(MYUNIT,'(A)') 'newres> Quench from reseeded geometry failed - try again'
            DO J2=1,3*NATOMS
               RANDOM=(DPRAND()-0.5D0)*2.0D0
               COORDS(J2,JP)=RANDOM*DSQRT(RADIUS)/SR3
            ENDDO
            GOTO 654
         ENDIF
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         POTEL=QENERGY
         EBEST(JP)=POTEL ! this is communicated via common block MYPOT, which is in quench.f
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=J1
         EPREV(JP)=POTEL !
         EPPREV(JP)=0.0D0
         NSYMREM=0
! th368: 20.10.2009 updates should only take place in case
! reseeding took place - transfered from subroutine MC lines 363 to 367

         IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,JP))
         COORDSO(1:3*(NATOMS-NSEED),JP)=COORDS(1:3*(NATOMS-NSEED),JP)
!        WRITE(MYUNIT,'(A,2G20.10)') 'newres> coordso changed: ',COORDSO(1,JP),COORDS(1,JP)
         VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)

! end th368 20.10.2009
      ENDIF

      RETURN
      END
