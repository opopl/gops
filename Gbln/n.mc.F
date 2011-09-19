      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)
      USE COMMONS
      USE qmodule , only : qmin, QMINP, INTEQMIN
      USE modcharmm
      USE modamber9, only : mdstept,cisarray1,cisarray2,chiarray1,chiarray2,nocistransdna,nocistransrna,
     &                      setchiral,amchpmax,doligmove,ligmovefreq, amchnmax
      USE porfuncs
      USE AMHGLOBALS, ONLY: NMRES,OMOVI,AVEP,NUMPRO,IRES
      USE AMH_INTERFACES, ONLY:E_WRITE

      IMPLICIT NONE
#ifdef MPI
      INTEGER  MPIERR
      INCLUDE 'mpif.h'
      INTEGER IS(MPI_STATUS_SIZE)
#endif

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
     &                 DUMMY1, DUMMY2, DUMMY3, INTE, OPOTEL, DUMGRAD(3*NATOMS), DJWPOTEL
      LOGICAL CHANGEDE, EXCHANGEACCEPT, EXCHANGE, FLAG 
      LOGICAL CHIRALFAIL,AMIDEFAIL, LOGDUMMY, DISTOK, ATOMINGROUP(NATOMS)
      INTEGER :: gly_count,iii,i2,i500,snapcount, DUMMYINT
      DOUBLE PRECISION prcord(NATOMS,3,3,3)
      DOUBLE PRECISION :: mctemp
      INTEGER :: PAIRCOUNTER
      DOUBLE PRECISION, EXTERNAL :: PAIRDISTANCE
      DOUBLE PRECISION :: ATOM1(3),ATOM2(3)

      LOGICAL EVAP, ATEST, STAY, evapreject, LOPEN

      integer                 :: values(8),itime1
      double precision :: DISTGROUPX2,DISTGROUPY2,DISTGROUPZ2,DISTGROUPCENTRE,TESTANGLE
      integer :: J6

      IF (FREEZEGROUPT) THEN
         OPEN(UNIT=4431,FILE='frozen.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         DO J6=1,NATOMS
            DISTGROUPX2=(COORDS(3*GROUPCENTRE-2,1)-COORDS(3*J6-2,1))**2
            DISTGROUPY2=(COORDS(3*GROUPCENTRE-1,1)-COORDS(3*J6-1,1))**2
            DISTGROUPZ2=(COORDS(3*GROUPCENTRE  ,1)-COORDS(3*J6  ,1))**2
            DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
            IF((FREEZEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
               NFREEZE=NFREEZE+1
               FROZEN(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'FREEZE ',J6
            ELSE IF((FREEZEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
               NFREEZE=NFREEZE+1
               FROZEN(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'FREEZE ',J6
            END IF
         END DO
         FREEZEGROUPT=.FALSE.     
      ENDIF

      IF (DONTMOVEGROUPT) THEN
              OPEN(UNIT=4431,FILE='dontmove.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         DO J6=1,NATOMS
            DISTGROUPX2=(COORDS(3*DONTMOVECENTRE-2,1)-COORDS(3*J6-2,1))**2
            DISTGROUPY2=(COORDS(3*DONTMOVECENTRE-1,1)-COORDS(3*J6-1,1))**2
            DISTGROUPZ2=(COORDS(3*DONTMOVECENTRE  ,1)-COORDS(3*J6  ,1))**2
            DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
            IF((DONTMOVEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
               NDONTMOVE=NDONTMOVE+1
               DONTMOVE(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'DONTMOVE ',J6
       ELSE IF((DONTMOVEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
               NDONTMOVE=NDONTMOVE+1
               DONTMOVE(J6)=.TRUE.
               WRITE(4431,'(A,I6)') 'DONTMOVE ',J6
            END IF
         END DO
         DONTMOVEGROUPT=.FALSE.     
      ENDIF
      

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
         IF (DEBUG) WRITE(MYUNIT, '(A,2I6)') 'NPAR,MPIERR=',NPAR,MPIERR
         JP=MYNODE+1
         ITRAJ=MYNODE
         IF (DEBUG) WRITE(MYUNIT, '(A,3I6)') 'In mc after MPI_MPI_COMM_RANK MPIERR,MYNODE,JP=',MPIERR,MYNODE,JP
      ENDIF
      NEACCEPT=0
#endif

      NDONE=0
      IF (RESTORET) THEN
#ifdef MPI
#else
         DO JP=1,NPAR
         ENDDO
#endif
         WRITE(MYUNIT, '(A,I10)') 'MC> restore NDONE=',NDONE
      ENDIF
      NQ(:)=NDONE

      
      IF (BSWL.and.(.not.TETHER)) then
         RETURN
      ELSEIF (TETHER) THEN
         RETURN
      ENDIF

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
         IF (SETCENT) CALL SETCENTRE(COORDS)
         IF (CHRMMT.AND.MAKEOLIGOT.AND.MAKEOLIGOSTART) THEN
#ifdef MPI
             IF (RANDOMSEEDT) THEN
                itime1= values(6)*60 + values(7)
             ENDIF
#else
             DO JP=1,NPAR
             ENDDO
#endif
         ENDIF
      ENDIF

#ifdef MPI
      WRITE(MYUNIT, '(A)')  'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      NQTOT=NQTOT+1
      WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART

 
 
      EPREV(JP)=POTEL
      EPPREV(JP)=0.0D0
      ELASTSYM(JP)=0.0D0
      EBEST(JP)=POTEL
      BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      JBEST(JP)=0
      RMIN=POTEL
      RCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,1)
      VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
      EPSSPHERE=EPSSAVE


      IF (PTTMIN < 0.000001D0 ) PTTMIN = 0.000001D0 ! to avoid devision by zero

      DO I=0, NPAR-1
         TEMPTRAJ(I)=PTTMIN*CTE**I
         T=TEMPTRAJ(I)
         BETA(I)=1.0D0/T
      ENDDO
      DO I=1, NPAR
         TEMP(I)=TEMPTRAJ(I-1)
      ENDDO
      VINIT=POTEL

      WRITE(MYUNIT, '(A, 2G20.10)') 'Temperature range', TEMPTRAJ(0), TEMPTRAJ(NPAR-1)
      WRITE(MYUNIT, '(A, G20.10)') 'For this replica T=', TEMPTRAJ(MYNODE)
      WRITE(MYUNIT, '(A, G20.10)') 'Starting potential energy=', VINIT

      RNDSEED=2002+MYNODE
      RANDOM=DPRAND()
      WRITE(MYUNIT, '(A, G20.10)') 'Starting random number=', RANDOM
      WRITE(ISTR,'(I3)') JP
      TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
      IF ((CHRMMT).AND.(.NOT.RESTORET)) CALL CHARMMDUMP(TEMPCOORDS,'initialmin.'//TRIM(ADJUSTL(ISTR)))

#else
      WRITE(MYUNIT,'(A)') 'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      DO JP=1,NPAR
         NQTOT=NQTOT+1
         IF (NPAR.GT.1) THEN
            WRITE(MYUNIT,'(A,I2,A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,']Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ELSE
            WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ENDIF


      IF (CHRMMT) CALL CHARMMDUMP(COORDS,'initialmin')

         IF (AMBERT.AND.NOCISTRANS.AND.(.NOT.NOCISTRANSDNA).AND.(.NOT.NOCISTRANSRNA)) THEN
            WRITE(MYUNIT,'(A)') ' mc> Storing cis/trans information for initial structure'
         ENDIF
         IF (AMBERT.AND.SETCHIRAL.AND.NOCISTRANS.AND.(.NOT.NOCISTRANSDNA).AND.(.NOT.NOCISTRANSRNA)) THEN
            WRITE(MYUNIT,'(A)') ' mc> Storing chiral information for initial structure'
         ENDIF 
 
         EPREV(JP)=POTEL
         EPPREV(JP)=0.0D0
         ELASTSYM(JP)=0.0D0
         EBEST(JP)=POTEL
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=0
         RMIN=POTEL
         RCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,1)
         VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
      ENDDO
      EPSSPHERE=EPSSAVE
#endif

      IF (THOMSONT) THEN
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

      IF (NPAR.EQ.1) THEN
         WRITE(MYUNIT,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
      ELSE
         WRITE(MYUNIT,'(A,I3,A,I10,A)') 'Starting ',NPAR,' parallel MC runs of ',NSTEPS,' steps'
      ENDIF
      WRITE(MYUNIT,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

      IF (LOCALSAMPLET.AND.AMBERT) THEN
         IF(NMOVABLEATOMS==0) THEN
            WRITE(MYUNIT,*) 'must have MOVABLEATOMS specified when using LOCALSAMPLE' 
            STOP
         ENDIF
         DISTOK=.FALSE.
         IF (.NOT.DISTOK) THEN
            WRITE(MYUNIT,*) 'initial structure violates LOCALSAMPLE conditions - exiting!'
            STOP
         ENDIF
         DISTOK=.FALSE.
      ENDIF   
      NSUPERCOUNT=NSUPER

      DO J1=NDONE+1,NSTEPS 
         ISTEP = J1

         IF (NEWJUMP) RANNJ=DPRAND()
#ifdef MPI
#else
         DO JP=1,NPAR 
#endif
            IF (QUCENTRE) THEN 
               TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
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
            MCTEMP = TEMP(JP)
            IF (JUMPMOVE(JP).AND.(MOD(J1,JUMPINT(JP)).EQ.0)) CALL JUMPM(RANNJ,J1,JP,EPPREV)
23          CONTINUE
            IF (SYMMETRIZE.AND.(MOD(NQ(JP),NSYMINTERVAL).EQ.0).AND.
     &               ((SYMMETRIZECSM).OR.(ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV))) THEN
               IF ((ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV)) NSYMCALL=0
               ELASTSYM(JP)=EPREV(JP)
                  IF (SYMMETRIZECSM) THEN
     &                                J1,NSUCCESS,NFAIL,EBEST,BESTCOORDS,JBEST,EPPREV)
                  ELSE
     &                                J1,EBEST,BESTCOORDS,JBEST,EPPREV)
                  ENDIF
                  IF (HIT) GOTO 37 
               POTEL=EPREV(JP) ! NEWRES assumes POTEL is the energy of the current structure in COORDS
               IF (CHANGEDE.AND.NEWRESTART) THEN
     1                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST)
               ENDIF
            ELSEIF (ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV) THEN ! Markov minimum has changed, but SYMMETRY not called
               NSYMREM=0                                       ! Should therefore reset NSYMREM.
            ENDIF
               IF(GROUPROTT.AND.MOD(J1,GROUPROTFREQ).EQ.0) THEN
                       DOGROUPROT=.TRUE.
               ENDIF
                  SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)

               IF (CHRMMT) THEN
                  IF (CHMDT.AND.MOD(J1,CHMDFREQ).EQ.0) THEN
                  ELSE
                     IF (CHRIGIDTRANST.AND.MOD(J1,FTRANS).EQ.0) CALL MKRIGIDTRANS(JP)
                     IF (CHRIGIDROTT.AND.MOD(J1,FROT).EQ.0) CALL MKRIGIDROT(JP)
                     IF (MOD(J1,CHFREQ).EQ.0) CALL TAKESTEPCH(JP)

                     IF (RANDOMSEEDT) THEN
                         itime1= values(6)*60 + values(7)
                     END IF
                              
                  ENDIF
                  IF(GROUPROTT.AND.DOGROUPROT) CALL GROUPROTSTEP(JP) 
               ELSE IF (ELLIPSOIDT.AND..NOT.PYGPERIODICT.AND..NOT.PYBINARYT) THEN

               ELSE IF (GAYBERNEDCT) THEN

               ELSE IF (MULTISITEPYT) THEN
                  IF(RANDOMSEEDT) THEN
                     itime1= values(6)*60 + values(7)
                  END IF

               ELSE IF (PYGPERIODICT.OR.PYBINARYT.OR.LJCAPSIDT) THEN
                  IF(RANDOMSEEDT) THEN
                     itime1= values(6)*60 + values(7)
                  END IF
                  IF(SWAPMOVEST) THEN
                  ELSE
                  END IF


               ELSE IF (LJCOULT) THEN
                  RANDOM=DPRAND()
                  IF (RANDOM < COULSWAP) THEN
                     MCTEMP = COULTEMP
                  ELSE
                  END IF


               ELSE IF (GBT .OR. GBDT .OR. GBDPT .OR. PYGT .OR. PYGDPT) THEN
      
               ELSE IF (MSGBT) THEN

               ELSE IF (MSPYGT) THEN
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
                     
                     IF (DEBUG) THEN
                        WRITE(MYUNIT, '(A)') '=== PRE-TAKESTEP COORDS ENERGY ==='
                        WRITE(MYUNIT, '(A, F20.10)') 'Energy = ', OPOTEL
                        WRITE(MYUNIT, '(A)') '=== PRE-TAKESTEP COORDSO ENERGY ==='
                        WRITE(MYUNIT, '(A, F20.10)') 'Energy = ', OPOTEL

                        OPEN(1473, FILE="prestepcoordsfile", POSITION="APPEND", STATUS="unknown", form="formatted")

                        WRITE(1473, '(A)') '====================='
                        WRITE(1473, '(A, I6)') 'Step ', NQ(JP)
                        WRITE(1473, '(A)') '====================='
                        WRITE(1473, '(A)') '====================='
                        WRITE(1473, '(A)') 'coordso'
                        WRITE(1473, '(A)') '====================='
                        WRITE(1473, '(3F20.10)') COORDSO(:,JP)
                     END IF


                     IF (DEBUG) THEN
                        OPEN(1483, FILE="poststepcoordsfile", POSITION="APPEND", STATUS="unknown", form="formatted")
                        WRITE(1483, '(A)') '====================='
                        WRITE(1483, '(A, I6)') 'Step ', NQ(JP)
                        WRITE(1483, '(A)') '====================='
                        WRITE(1483, '(A)') '====================='
                        WRITE(1483, '(A)') 'coordso'
                        WRITE(1483, '(A)') '====================='
                        WRITE(1483, '(3F20.10)') COORDSO(:,JP)
 
                     
                        WRITE(MYUNIT, '(A)') '=== POST-TAKESTEP COORDS ENERGY ==='
                        WRITE(MYUNIT, '(A, F20.10)') 'Energy = ', OPOTEL
                        WRITE(MYUNIT, '(A)') '=== POST-TAKESTEP COORDSO ENERGY ==='
                        WRITE(MYUNIT, '(A, F20.10)') 'Energy = ', OPOTEL
                     END IF
                     
                     doligmove=.FALSE.
                     IF (RANDOMSEEDT) THEN
                        ITIME1= VALUES(6)*60 + VALUES(7)
                     ENDIF
                     IF (AMCHPMAX.EQ.0) THEN
                     ELSE
                     ENDIF
                     IF(GROUPROTT.AND.DOGROUPROT) CALL GROUPROTSTEP(JP) 
                     IF (LOCALSAMPLET) THEN
                        IF (.NOT.DISTOK) COORDS(:,JP)=SAVECOORDS(:)
                     ELSE
                        DISTOK=.TRUE.
                     ENDIF
                  ENDDO
               ELSE
 
               ENDIF
               IF ((CHNMAX.LE.0.0D0).AND.(AMCHNMAX.LE.0.0D0)) THEN
                  IF(FREEZE) THEN
                     DO J2=1,NATOMS
                        IF (FROZEN(J2)) THEN
                        ENDIF
                     ENDDO
                  ENDIF
                  IF(DONTMOVET) THEN
                     DO J2=1,NATOMS
                        IF (DONTMOVE(J2)) THEN
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
 
                     dogrouprot=.FALSE.
               
            
                  IF (DUMPSTEPST) THEN
                      WRITE(QUENCHNUM,*) NQ(JP)
                      QUNAME='afterstep'//TRIM(ADJUSTL(QUENCHNUM))//'.rst'
                      OPEN(UNIT=20,FILE=QUNAME,STATUS='UNKNOWN')  
                      WRITE(20,'(a20)') QUNAME
                      WRITE(20,'(i5)') NATOMS
                      WRITE(20,'(6f12.7)') COORDS(:,JP) 
                      QUNAME='afterstep'//TRIM(ADJUSTL(QUENCHNUM))
                  ENDIF

               NQ(JP)=NQ(JP)+1
               IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
               NQTOT=NQTOT+1
 
 
               IF (SELFT) THEN
                  IF (DEBUG) THEN
                     WRITE(MYUNIT,'(A)') 'Taboo list:'
                     WRITE(MYUNIT,'(6F20.10)') (ESAVE(J2,1),J2=1,NT(1))
                     WRITE(MYUNIT,73) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
73                   FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' TABOO')
                  ENDIF
                  GOTO 23
               ENDIF
   
               IF (NPAR.GT.1) THEN
                  WRITE(MYUNIT,'(A,I2,A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,']Qu ',NQ(JP),' E=',
     1                 POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART
               ELSE
                  WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1                 POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART



                  IF ((MOD(J1,NINT_AMH).EQ.0).AND.AMHT) THEN

                     GLY_COUNT = 0
                     SNAPCOUNT = SNAPCOUNT + 1

                     DO III = 1,NMRES
                        IF (IRES(III).EQ.8) THEN
                           PRCORD(III,1,1,1) = COORDS(9*(III-1)+1- GLY_COUNT*3,JP) !  CA X
                           PRCORD(III,2,1,1) = COORDS(9*(III-1)+2- GLY_COUNT*3,JP) !  CA Y
                           PRCORD(III,3,1,1) = COORDS(9*(III-1)+3- GLY_COUNT*3,JP) !  CA Z
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

                  ENDIF  !                  IF (MOD(J1,NINT_AMH).EQ.0)
               ENDIF     !                  IF (NPAR.GT.1)

          IF (RMST.AND.CHRMMT) THEN
             IF (DEBUG) WRITE(MYUNIT,'(A,F15.5)')'RMSD = ',RMSD
             IF (RMSD.LE.RMSLIMIT) CALL SAVERMS(JP,POTEL,RMSD)
          ENDIF

          IF (PAIRDISTT) THEN
                WRITE(MYPUNIT,*) " "
                WRITE(MYPUNIT,'(I10)',ADVANCE="NO") NQ(JP)
             DO PAIRCOUNTER=1,NPAIRS
                ATOM1(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,1)-2:3*PAIRDIST(PAIRCOUNTER,1),JP)
                ATOM2(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,2)-2:3*PAIRDIST(PAIRCOUNTER,2),JP)
                WRITE(MYPUNIT,'(F10.4)',ADVANCE="NO") PAIRDISTANCE(ATOM1,ATOM2) 
             ENDDO
             FLUSH(MYPUNIT)
          ENDIF
          
          IF (TRACKDATAT) THEN
             WRITE(MYEUNIT,'(I10,F20.10)') J1,POTEL
             WRITE(MYMUNIT,'(I10,G20.10)') J1,EPREV(JP)
             WRITE(MYBUNIT,'(I10,G20.10)') J1,QMIN(1)
             IF (A9INTET.AND.(NQ(JP).GT.2)) THEN
                WRITE(3999,'(I10,G20.10)') NQ(JP),INTEQMIN(1)
             ENDIF
             IF (RMST.AND.CHRMMT) THEN
                WRITE(4428,'(I10,F15.5)') NQ(JP),RMSD
             ENDIF
          ENDIF

 
            IF (DAESTAT) THEN
               PRINT*,'DAESTAT block in mc.f not implemented'
               STOP
            ENDIF
 
 
            IF (NEWRESTART.AND.(.NOT.SEEDT)) THEN 
     &                    JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST)
            ENDIF
            IF (STAY) THEN
            ELSE IF (EVAP .and. .not.evapreject) THEN
               NFAIL(JP)=NFAIL(JP)+1
               IF (DEBUG) THEN
                  WRITE(MYUNIT,33) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
33                FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' EVAP,REJ')
               ENDIF
            ELSE
               ATEST=.TRUE. 
               IF (COLDFUSION) THEN
                  ATEST=.FALSE.
               ENDIF

               IF (PERCOLATET) THEN
                 IF(PERCACCEPTED .AND. (.NOT. PERCT)) THEN
                   WRITE(MYUNIT,*) 'perc> Disconnected atoms, rejecting step'
                   ATEST=.FALSE.
                 ELSE IF((.NOT. PERCACCEPTED) .AND. PERCT) THEN
                   PERCACCEPTED=.TRUE.
                   WRITE(MYUNIT,*) 'perc> Found first connected configuration'
                 ENDIF
               ENDIF

               IF (.NOT.GOODSTRUCTURE) ATEST=.FALSE.

               IF (AMBERT) THEN
                  IF (A9INTET.AND.ATEST) THEN
                     IF (DEBUG) WRITE(MYUNIT,*) 'intE=',INTE
                     IF (TRACKDATAT) THEN
                        WRITE(3998,'(I10,G20.10)') NQ(JP),INTE
                     ENDIF
                  ENDIF
                  IF (DUMPQUT.AND.ATEST) THEN
                      WRITE(QUENCHNUM,*) NQ(JP)
                      QUNAME='quench'//TRIM(ADJUSTL(QUENCHNUM))//'.rst'
                      OPEN(UNIT=20,FILE=QUNAME,STATUS='UNKNOWN')  
                      WRITE(20,'(a20)') QUNAME
                      WRITE(20,'(i5)') NATOMS
                      WRITE(20,'(6f12.7)') COORDS(:,JP) 
                     QUNAME='quench'//TRIM(ADJUSTL(QUENCHNUM))
                  ELSEIF (DUMPQUT) THEN
                      WRITE(*,'(A,I6)') 'DUMPQU> chirality/cis-trans problem detected, not dumping quench ',NQ(JP)
                  ENDIF
                  IF (DUMPBESTT) THEN
                     WRITE(QUENCHNUM,*) NQ(JP)
                     BESTNAME='best_'//TRIM(ADJUSTL(QUENCHNUM))
                  ENDIF
               ENDIF
               IF (LOCALSAMPLET) THEN
                  DISTOK=.FALSE.
                  IF (.NOT.DISTOK) ATEST=.FALSE.
               ENDIF

               IF ((QDONE.EQ.0).AND.TIP) THEN
                  ATEST=.FALSE.
               ELSEIF (ATEST) THEN
               ENDIF
               IF (DEBUG.OR.CHECKMARKOVT) THEN
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
                     JACCPREV=J1
                  ENDIF
                  IF (QDONE.EQ.1) THEN
                     NSUCCESS(JP)=NSUCCESS(JP)+1
                  ELSE
                     NFAIL(JP)=NFAIL(JP)+1
                  ENDIF
                  EPPREV(JP)=EPREV(JP)
                  EPREV(JP)=POTEL
                  VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
               ELSE
                  NFAIL(JP)=NFAIL(JP)+1
                  IF (DEBUG) THEN
                     WRITE(MYUNIT,36) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
36                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
               ENDIF
            ENDIF
            IF ((NPAR.GT.1).AND.(.NOT.NEWJUMP)) CALL DUMPJ(JP,JUMPTO,NPAR,COORDS(1:3*NATOMS,1:NPAR),NATOMS,EPREV)
            IF (RESTART.AND.(J1-JACCPREV.GT.1.1D0*NRELAX)) CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
            IF ((MOD(J1,NACCEPT).EQ.0).AND.(NSEED.EQ.0).AND.(.NOT.STAY)) CALL ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)

            TEMP(JP)=TEMP(JP)*SCALEFAC
#ifdef MPI
#else
            IF (HIT) GOTO 37
#endif
            IF (DUMPINT.GT.0) THEN
               IF (MOD(J1,DUMPINT).EQ.0) THEN
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
          IF (MYNODE.EQ.(J+1)) THEN
          ENDIF
          IF (MYNODE.EQ.J) THEN

             DBETA=BETA(J)-BETA(J+1)
             DELTA=E-ER
             W=MIN(1.0D0,DEXP(DELTA*DBETA))
             NTOT=NTOT+1 
             RANDOM=DPRAND()
             IF (W.GT.RANDOM) THEN
                EXCHANGEACCEPT=.TRUE.
                IMESG=1

                E=ER
                NEACCEPT=NEACCEPT+1
              ELSE
                EXCHANGEACCEPT=.FALSE.
                IMESG=0
              ENDIF
          ENDIF
          IF (MYNODE.EQ.(J+1)) THEN
             NTOT=NTOT+1
             IF (IMESG.EQ.1) THEN
                ITRAJ=ITRAJO
                NEACCEPT=NEACCEPT+1
              ENDIF
           ENDIF
  
           IF (NQ(1).GT.NSTEPS) EXIT
#else
           IF (NQ(JP).GT.NSTEPS) GOTO 37

         ENDDO
#endif
  
      ENDDO
#ifdef MPI
      WRITE(MYUNIT, '(G20.10,A,G20.10)') NEACCEPT, ' exchanges accepted out of ', NTOT
      IF (MPIT) THEN
      ENDIF
#endif

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
      DEALLOCATE(TMOVE)
      DEALLOCATE(OMOVE)
      RETURN
      END
      SUBROUTINE A9INTEDUMPLOWEST()
      USE COMMONS
      USE QMODULE
      IMPLICIT NONE
      OPEN(UNIT=20,FILE="bestintE.rst",STATUS='UNKNOWN')  
      WRITE(20,'(g20.10,i5)') INTEQMIN(1), INTEFF(1)
      WRITE(20,'(i5)') NATOMS
      WRITE(20,'(6f12.7)') INTEQMINP(1,:) 
      END SUBROUTINE A9INTEDUMPLOWEST
      SUBROUTINE A9INTE(INTECOORDS,INTE)
      USE PORFUNCS
      USE COMMONS, ONLY : NATOMS
      IMPLICIT NONE
      DOUBLE PRECISION :: INTECOORDS(3*NATOMS), INTE
      INTEGER ISTAT
      OPEN(UNIT=20,FILE='coords.intres',STATUS='UNKNOWN')
      WRITE(20,'(a20)') 'intE coordinates'
      WRITE(20,'(i5)') NATOMS
      WRITE(20,'(6f12.7)') INTECOORDS(:)
      OPEN(UNIT=20,FILE="intE",STATUS='OLD')
      READ(20,*) INTE
      END SUBROUTINE A9INTE

      SUBROUTINE A9DISTCHECK(COORDS,DISTOK)
      USE MODAMBER9, ONLY : NATOMSINA,NATOMSINB,NATOMSINC,ATOMSINALIST,ATOMSINBLIST,ATOMSINCLIST
      USE COMMONS, ONLY : NATOMS,ABTHRESH,ACTHRESH,DEBUG
      IMPLICIT NONE
      LOGICAL :: DISTOK
      INTEGER :: I,J
      DOUBLE PRECISION :: COORDS(3*NATOMS)
      DOUBLE PRECISION :: CENTREOFA(3),CENTREOFB(3),CENTREOFC(3)  
      DOUBLE PRECISION :: ABDIST,ACDIST 
      ABDIST=0.0D0
      ACDIST=0.0D0
      DO I=1,NATOMSINA
         J=ATOMSINALIST(I)
      END DO
      DO I=1,NATOMSINB
         J=ATOMSINBLIST(I)
      END DO
      DO I=1,NATOMSINC
         J=ATOMSINCLIST(I)
      END DO
      ABDIST=SQRT((CENTREOFA(1)-CENTREOFB(1))**2+(CENTREOFA(2)-CENTREOFB(2))**2+(CENTREOFA(3)-CENTREOFB(3))**2)
      ACDIST=SQRT((CENTREOFA(1)-CENTREOFC(1))**2+(CENTREOFA(2)-CENTREOFC(2))**2+(CENTREOFA(3)-CENTREOFC(3))**2)
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
      IF ((ABDIST.LT.ABTHRESH).AND.(ACDIST.LT.ACTHRESH)) DISTOK=.TRUE.
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

      IF (DISTMINOLD.EQ.0.0D0) DISTMINOLD=DISTMIN  ! this should allow for the first step
      IF (TUNNELT) THEN
         TEOLD=TRANS(EOLD,QMIN(1),GAMMA)
         TENEW=TRANS(ENEW,QMIN(1),GAMMA)
      ELSE
         TEOLD=EOLD
         TENEW=ENEW
      ENDIF

      IF (TSALLIST) THEN
         EREF=QMIN(NP)*1.1D0
         DUMMY=(1.0D0-(1.0D0-QTSALLIS)*(TENEW-EREF)/MCTEMP)/(1.0D0-(1.0D0-QTSALLIS)*(TEOLD-EREF)/MCTEMP)
         DUMMY=DUMMY**(QTSALLIS/(1.0D0-QTSALLIS))
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
                  DUMMY=COORDSO(J2,JP)
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
            DO J2=1,NATOMS
               VATO(J2,JP)=VAT(J2,JP)
            ENDDO
            DO J2=1,3*NATOMS
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

      RETURN
      END

      SUBROUTINE MYRESET(JP,NATOMS,NPAR,NSEED)
      USE COMMONS,ONLY : MYUNIT,COORDS,COORDSO,VAT,VATO
      IMPLICIT NONE
      INTEGER JP, NSEED, J2, NATOMS, NPAR

      DO J2=1,3*(NATOMS-NSEED)
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

10    CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),1,NHSRESTART)
      IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
         DO J2=1,3*NATOMS
         ENDDO
         GOTO 10
      ENDIF
      JACCPREV=J1
      NQTOT=NQTOT+1
      WRITE(MYUNIT,'(A,I6,A)') ' Restarting using ',NHSRESTART,' hard sphere moves'
      WRITE(MYUNIT,'(A,I7,A,F20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Restart Qu ',NQ(1),' E=',
     1              POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
      DO J2=1,3*NATOMS
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
            IF (PERCOLATET .AND. ( STEP(JP) .GT. PERCCUT))  STEP(JP) = PERCCUT
            IF (PERCOLATET .AND. (ASTEP(JP) .GT. PERCCUT)) ASTEP(JP) = PERCCUT
            IF (PERCOLATET .AND. (OSTEP(JP) .GT. PERCCUT)) OSTEP(JP) = PERCCUT
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
      STEP(JP)=MIN(STEP(JP),1.0D3)
      OSTEP(JP)=MIN(OSTEP(JP),1.0D3)
      ASTEP(JP)=MIN(ASTEP(JP),1.0D3)
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
      NSUCCESST(JP)=NSUCCESST(JP)+NSUCCESS(JP)
      NFAILT(JP)=NFAILT(JP)+NFAIL(JP)
      NSUCCESS(JP)=0
      NFAIL(JP)=0 
      RETURN
      END

      SUBROUTINE REN(J1,RMIN,RCOORDS,RVAT,NREN,RMINO,RCOORDSO,RVATO,ITERATIONS,TIME,NLAST,JACCPREV,NSTEPREN)
      USE commons
      IMPLICIT NONE
      LOGICAL STAY, REJECT, METROPOLIS
      INTEGER J1, J2, NREN, ITERATIONS, NQTOT, NLAST, JACCPREV, NSTEPREN, J3, BRUN, QDONE
      DOUBLE PRECISION POTEL, RMIN, RCOORDS(3*NATOMS), RVAT(NATOMS), RANDOM, DPRAND, RMINO, RCOORDSO(3*NATOMS), RVATO(NATOMS),
     1                 TIME, XIP, DUMMY, SCREENC(3*NATOMS)

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
      IF ((J1.GE.JACCPREV+NRENSTUCK).OR.(J1.GE.NLAST+NREN)) THEN
         JACCPREV=J1
         NLAST=J1
         RANDOM=DPRAND()
         METROPOLIS=DEXP(-(RMIN-RMINO)/TRENORM).GT.RANDOM
         REJECT=.FALSE.
         IF (TABOOT.AND.METROPOLIS) THEN
            IF (NSTEPREN.EQ.0) NT(1)=0
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
         IF (METROPOLIS.AND.(.NOT.REJECT)) THEN
            IF (NSTEPREN.GT.0) WRITE(MYUNIT,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' accepted'
            NREN=MAX(NREN/1.1D0,NRENORM/2.0D0)
            RMINO=RMIN
            IF (.NOT.STAY) THEN
               DO J2=1,3*NATOMS
                  RCOORDSO(J2)=RCOORDS(J2)
               ENDDO
               DO J2=1,NATOMS
                  RVATO(J2)=RVAT(J2)
               ENDDO
            ELSE
               DO J2=1,3*NATOMS
                  RCOORDSO(J2)=COORDSO(J2,1)
               ENDDO
            ENDIF
         ELSE
            IF (REJECT) THEN
               WRITE(MYUNIT,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected by taboo criterion'
            ELSE
               WRITE(MYUNIT,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected'
            ENDIF
            DO J2=1,3*NATOMS
            ENDDO
            DO J2=1,NATOMS
               RVAT(J2)=RVATO(J2)
            ENDDO
            NREN=NREN*1.1D0
         ENDIF
         NSTEPREN=NSTEPREN+1
         IF (NSTEPREN.EQ.1) WRITE(MYUNIT,'(A,G20.10)') ' first renorm energy is ',RMIN
         WRITE(MYUNIT,'(A,I6)') ' renormalisation interval is now ',NREN
10       IF ((XMOVERENORM.GT.3.0D0).OR.FIXD) THEN
         ELSE
            DUMMY=STEP(1)
            STEP(1)=XMOVERENORM
            STEP(1)=DUMMY
         ENDIF
         IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
            DO J2=1,3*NATOMS
            ENDDO
            GOTO 10
         ENDIF
         NQTOT=NQTOT+1
         WRITE(MYUNIT,'(A,I7,A,F20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Renorm Qu ',NQ(1),' E=',
     1        POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
         DO J2=1,3*NATOMS
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

      SR3=DSQRT(3.0D0)
      IF (POTEL.LT.EBEST(JP)) THEN
         IF (EBEST(JP)-POTEL.GT.ECONV) JBEST(JP)=J1
         EBESTP=EBEST(JP)
         BESTCOORDSP(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP) ! save previous BESTCOORDS for possible use
         EBEST(JP)=POTEL ! reset ebest, but not necessarily jbest
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      ENDIF
      RES1=.FALSE.
      IF (J1-JBEST(JP).GT.NRELAX) RES1=.TRUE.
      RES2=.FALSE.
      IF ((.NOT.RES1).AND.AVOID.AND.(J1.EQ.JBEST(JP)).AND.(NMSBSAVE.GT.0)) THEN ! best minimum has just changed.
         FCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
         savedloop: DO J2=1,MIN(NMSBSAVE,MAXSAVE)
            IF (POTEL.LT.MSBE(J2)-ECONV) CYCLE savedloop
            XMSB(1:3*NATOMS)=MSBCOORDS(1:3*NATOMS,J2)
            IF (DISTANCE.LT.AVOIDDIST) THEN
               RES2=.TRUE.
               WRITE(MYUNIT,'(A,G20.10,A,I6,A,G20.10,A,F10.3)') 'newres> Minimum energy ',POTEL,
     &                                    ' is too close to saved structure ',
     &                                 J2,' with energy ',MSBE(J2),' dist=',DISTANCE
               GOTO 20
            ENDIF
         ENDDO savedloop
      ENDIF

20    CONTINUE 


      IF (RES1.OR.RES2) THEN

         IF (RES2.AND.(.NOT.AVOIDRESEEDT)) THEN 
            POTEL=MAX(1.0D10,10.0D0*POTEL) ! This should be enough to reject the step until we do a massive Thomson problem!
            WRITE(MYUNIT,'(A,I8,A)') 'newres> Resetting energy to a large value due to taboo condition'
            RETURN
         ELSEIF (RES1.OR.(RES2.AND.(DISTANCE.GT.0.1D0))) THEN ! new condition
            NMSBSAVE=NMSBSAVE+1
            NPOSITION=MOD(NMSBSAVE,MAXSAVE)
            IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
            MSBE(NPOSITION)=EBEST(JP)
            FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP)


            MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)

            WRITE(MYUNIT,'(A,I6,A,G20.10)') 'newres> Moving current best minimum to position ',NPOSITION,
     &                         ' in the AVOID list E=',EBEST(JP)
         ELSEIF (RES2.AND.(DISTANCE.LE.0.01D0)) THEN
               NMSBSAVE=NMSBSAVE+1
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
               MSBE(NPOSITION)=EBESTP
               FCOORDS(1:3*NATOMS)=BESTCOORDSP(1:3*NATOMS)
               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
               WRITE(MYUNIT,'(A,I6,A,G20.10)') 'newres> Moving previous best minimum to position ',NPOSITION,
     1                             ' in the AVOID list E=',EBESTP
         ENDIF ! end new condition

         IF (NEWRESTART_MD) THEN ! lb415
            IF (RES1) WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
            WRITE(MYUNIT,'(A,I8,A)') 'newres> Reseeding via a short high temperature MD run'


             IF (AMBERT) THEN
             ELSEIF (CHRMMT) THEN
             ELSE
               WRITE(MYUNIT,'(A,I8,A)') 'newres> Molecular Dynamics reseeding is available for AMBER or CHARMM runs only.'
               STOP
             ENDIF

         ELSE
            IF (NHSRESTART.GT.0) THEN
               IF (RES1) WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
               IF (RES2) WRITE(MYUNIT,'(A,I8,A)') 'newres> Reseeding due to taboo condition'
               IF (SHELLMOVES(JP)) WRITE(MYUNIT,'(A)') 'newres> Turning off shell moves'
               SHELLMOVES(JP)=.FALSE.
            ELSE
               IF (RES1) THEN
                  IF (AVOIDRESEEDT) THEN
                     WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - reseeding'
                  ELSE
                     WRITE(MYUNIT,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP)
                  ENDIF
               ENDIF
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
                  IF (.NOT.AVOIDRESEEDT) THEN
                     WRITE(MYUNIT,'(A,I8,A)') 'newres> Resetting energy to a large value to reject this step (not reseeding)'
                     JBEST(JP)=J1
                     EBEST(JP)=POTEL ! this is communicated via common block MYPOT, which is in quench.f
                     POTEL=MAX(1.0D10,10.0D0*POTEL) ! This should be enough to reject the step until we do a massive Thomson problem!
                     RETURN
                  ELSE
                     WRITE(MYUNIT,'(A)') 'newres> Full reseeding'
                     DO J2=1,3*NATOMS
                        RANDOM=(DPRAND()-0.5D0)*2.0D0
                     ENDDO
                     NCORE(JP)=0
                     PTGROUP(JP)='   '
                  ENDIF
               ELSEIF (NCORE(JP).GT.0) THEN
                  WRITE(MYUNIT,'(A)') 'newres> Accepting an enforced surface reseeding'
                  DUMMY2=-1.0D0
                  DO J2=NATOMS-NCORE(JP)+1, NATOMS
                     DISTANCE=COORDS(3*(J2-1)+1,JP)**2+COORDS(3*(J2-1)+2,JP)**2+COORDS(3*(J2-1)+3,JP)**2
                     IF (DISTANCE.GT.DUMMY2) DUMMY2=DISTANCE
                  ENDDO
                  DISTANCE=SQRT(DISTANCE)
                  WRITE(MYUNIT,'(A,F15.5)') 'newres> largest core radius=',DISTANCE
                  DO J2=1,NATOMS-NCORE(JP)
                     RANDOM=DISTANCE+1.0D0+DPRAND()*(SQRT(RADIUS)-DISTANCE-1.0D0)
                     DUMMY2=SQRT(COORDS(3*(J2-1)+1,JP)**2+COORDS(3*(J2-1)+2,JP)**2+COORDS(3*(J2-1)+3,JP)**2)
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
654      CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         IF (RMS.GT.BQMAX) THEN
            WRITE(MYUNIT,'(A)') 'newres> Quench from reseeded geometry failed - try again'
            DO J2=1,3*NATOMS
               RANDOM=(DPRAND()-0.5D0)*2.0D0
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

         IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,JP))
         VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)

      ENDIF

      RETURN
      END
