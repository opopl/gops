
      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)

      USE COMMONS
      USE QMODULE , ONLY : QMIN, INTEQMIN

      USE PORFUNCS

      IMPLICIT NONE

      SNAPCOUNT=0
      NSTEPREN=0
      EVAPREJECT=.FALSE.

      NRMS=0
      NLAST=0
      STAY=.FALSE.
      JACCPREV=0
      NQTOT=0
      RMINO=1.0D100
      RMIN=1.0D100
      NREN=NRENORM

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

C  Calculate the initial energy and save in EPREV
!op226>{{{ 
      WRITE(MYUNIT,'(A)') 'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      DO JP=1,NPAR
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NQTOT=NQTOT+1
            WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',
     1           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART



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
!op226>}}} 

      WRITE(MYUNIT,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
      WRITE(MYUNIT,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

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

                  SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
 
                  CALL TAKESTEP
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
 
               ENDIF
               
               NQ=NQ+1
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
           WRITE(MYUNIT,21) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),
     1               STEP(JP),ASTEP(JP),TEMP(JP)
21         FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
     
      RETURN
!op226>}}} 
      END
C

