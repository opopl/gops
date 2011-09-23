

!op226>=================================== 
!op226> GPL License Info {{{ 
!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!op226>}}} 
!op226>=================================== 
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
      INTEGER J1, NSUCCESS(NPAR), NFAIL(NPAR), NFAILT(NPAR), NSUCCESST(NPAR), J2, JP, J5

      INTEGER UNT
      integer ITERATIONS, NSUPERCOUNT, NQTOT, JACCPREV, NREN, NLAST, NSTEPREN, BRUN,QDONE,JBEST(NPAR)
      integer NRMS, NDONE, I, RNDSEED, J, NTOT, IMESG, ITRAJ, ITRAJO, NEACCEPT
      integer J3, J4, ISTAT, LOCALCOUNT
      INTEGER :: NSYMCALL=0
      DOUBLE PRECISION POTEL, RANDOM, DPRAND, SAVECOORDS(3*NATOMS), TEMPCOORDS(3*NATOMS)
      DOUBLE PRECISION :: TIME, SPOTEL(NSUPER), SCOORDS(3*NATOMS,NSUPER)
      DOUBLE PRECISION :: EPPREV(NPAR), QSTART, QFINISH, RANNJ, RMIN
      DOUBLE PRECISION :: RMINO, RCOORDS(3*NATOMS),ELASTSYM(NPAR)
      DOUBLE PRECISION ::  RCOORDSO(3*NATOMS), RVAT(NATOMS), RVATO(NATOMS), EPSSAVE, EBEST(NPAR)
      DOUBLE PRECISION ::  BESTCOORDS(3*NATOMS,NPAR), endtime, RMSD, VINIT, CTE, TEMPTRAJ(0:NPAR-1)
      DOUBLE PRECISION ::  T, BETA(0:NPAR-1), GRAD(3*NATOMS), E, ER, W, DELTA, DBETA, A9ANGLE 
      DOUBLE PRECISION ::  DUMMY1, DUMMY2, DUMMY3, INTE, OPOTEL, DUMGRAD(3*NATOMS), DJWPOTEL
      LOGICAL CHANGEDE, EXCHANGEACCEPT, EXCHANGE, FLAG 
      LOGICAL CHIRALFAIL,AMIDEFAIL, LOGDUMMY, DISTOK, ATOMINGROUP(NATOMS)
      CHARACTER FNAME*9
      CHARACTER (LEN= 3)  ISTR
      CHARACTER (LEN=20) QUENCHNUM, QUNAME, DUMMYCHAR
      CHARACTER (LEN=20) BESTNAME, CURRENTBESTNAME
!  AMH 
      INTEGER :: gly_count,iii,i2,i500,snapcount, DUMMYINT
      DOUBLE PRECISION prcord(NATOMS,3,3,3)
      DOUBLE PRECISION :: mctemp
!  csw34> PAIRDIST variables
      INTEGER :: PAIRCOUNTER
      DOUBLE PRECISION, EXTERNAL :: PAIRDISTANCE
      DOUBLE PRECISION :: ATOM1(3),ATOM2(3)

      LOGICAL EVAP, ATEST, STAY, evapreject, LOPEN

      character(len=10)       :: datechar,timechar,zonechar
      integer                 :: values(8),itime1
      double precision :: DISTGROUPX2,DISTGROUPY2,DISTGROUPZ2,DISTGROUPCENTRE,TESTANGLE
      integer :: J6
      ! }}}
      ! common {{{
      COMMON /EV/ EVAP, evapreject
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT
      COMMON /Q4C/ QSTART, QFINISH
      ! }}}
!op226>}}} 
!op226> Subroutine body {{{ 

! Write a list of FROZEN atoms for use in an (o)data file
!op226> IF (FREEZEGROUPT) THEN {{{
      !IF (FREEZEGROUPT) THEN
         !OPEN(UNIT=4431,FILE='frozen.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         !DO J6=1,NATOMS
!!
!! Work out the distance from GROUPCENTRE to the current atom J1
!! 
            !DISTGROUPX2=(COORDS(3*GROUPCENTRE-2,1)-COORDS(3*J6-2,1))**2
            !DISTGROUPY2=(COORDS(3*GROUPCENTRE-1,1)-COORDS(3*J6-1,1))**2
            !DISTGROUPZ2=(COORDS(3*GROUPCENTRE  ,1)-COORDS(3*J6  ,1))**2
            !DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
!! If working in GT mode (default), FREEZE all atoms >GROUPRADIUS from the GROUPCENTRE atom
            !IF((FREEZEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
               !NFREEZE=NFREEZE+1
               !FROZEN(J6)=.TRUE.
               !WRITE(4431,'(A,I6)') 'FREEZE ',J6
!! If working in LT mode, FREEZE all atoms <GROUPRADIUS from the GROUPCENTRE atom
            !ELSE IF((FREEZEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
               !NFREEZE=NFREEZE+1
               !FROZEN(J6)=.TRUE.
               !WRITE(4431,'(A,I6)') 'FREEZE ',J6
            !END IF
         !END DO
         !CLOSE(4431)
!! Prevent it doing this again
         !FREEZEGROUPT=.FALSE.     
      !ENDIF
!op226>}}} 

! Write a list of DONTMOVE atoms for use in an (o)data file
!!op226> IF (DONTMOVEGROUPT) THEN {{{
      !IF (DONTMOVEGROUPT) THEN
              !OPEN(UNIT=4431,FILE='dontmove.dat',STATUS='UNKNOWN',FORM='FORMATTED')
         !DO J6=1,NATOMS
!!
!! Work out the distance from DONTMOVECENTRE to the current atom J1
!! 
            !DISTGROUPX2=(COORDS(3*DONTMOVECENTRE-2,1)-COORDS(3*J6-2,1))**2
            !DISTGROUPY2=(COORDS(3*DONTMOVECENTRE-1,1)-COORDS(3*J6-1,1))**2
            !DISTGROUPZ2=(COORDS(3*DONTMOVECENTRE  ,1)-COORDS(3*J6  ,1))**2
            !DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
!! If working in GT mode (default), DONTMOVE all atoms >GROUPRADIUS from the DONTMOVECENTRE atom
            !IF((DONTMOVEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
               !NDONTMOVE=NDONTMOVE+1
               !DONTMOVE(J6)=.TRUE.
               !WRITE(4431,'(A,I6)') 'DONTMOVE ',J6
!! IF working in LT mode, DONTMOVE all atoms <GROUPRADIUS from the DONTMOVECENTRE atom
       !ELSE IF((DONTMOVEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
               !NDONTMOVE=NDONTMOVE+1
               !DONTMOVE(J6)=.TRUE.
               !WRITE(4431,'(A,I6)') 'DONTMOVE ',J6
            !END IF
         !END DO
         !CLOSE(4431)
!! Prevent it doing this again
         !DONTMOVEGROUPT=.FALSE.     
      !ENDIF
!!op226>}}} 
      
!     csw34> Added defaults to prevent accidentaly discarding
!     structures for AMH      

      !CHIRALFAIL=.FALSE.
      !AMIDEFAIL=.FALSE.
      !INQUIRE(UNIT=1,OPENED=LOPEN)
      !IF (LOPEN) THEN
         !WRITE(*,'(A,I2,A)') 'mc> A ERROR *** Unit ', 1, ' is not free '
         !STOP
      !ENDIF

      ALLOCATE(TMOVE(NPAR), OMOVE(NPAR))
      !snapcount=0
      !NSTEPREN=0
      EVAPREJECT=.FALSE.
      !INQUIRE(UNIT=1,OPENED=LOPEN)
      !IF (LOPEN) THEN
         !WRITE(*,'(A,I2,A)') 'mc> B ERROR *** Unit ', 1, ' is not free '
         !STOP
      !ENDIF

      NDONE=0
      IF (RESTORET) THEN
         DO JP=1,NPAR
            CALL RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
         ENDDO
         WRITE(LFH, '(A,I10)') 'MC> restore NDONE=',NDONE
!     csw34> Sets the quench counter so that the GMIN_out file makes sense after using RESTORE!
!         NQ(:)=NDONE
      ENDIF
      NQ(:)=NDONE

!! tvb requesting a basin-sampling MC run: {{{
      
      !IF (BSWL.and.(.not.TETHER)) then
         !CALL BasinSampling
         !RETURN
      !ELSEIF (TETHER) THEN
         !CALL TetheredWL
         !RETURN
      !ENDIF
!! }}}

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
         IF (JDUMP(JP).AND.(.NOT.NEWJUMP)) THEN
            WRITE(FNAME,'(A,I1)') 'ebuffer.',JP
            UNT=70+JP
            OPEN(UNIT=UNT,FILE=FNAME,STATUS='UNKNOWN')
            WRITE(FNAME,'(A,I1)') 'cbuffer.',JP
            UNT=70+NPAR+JP
            OPEN(UNIT=UNT,FILE=FNAME,STATUS='UNKNOWN')
         ENDIF
      ENDDO

      !IF (AMHT) THEN
         !write(omovi,1334)nmres,3,1,INT(real(NSTEPS)/real(NINT_AMH))
!1334     format(4(i8,1x),' nmres nmcrd numpro nmsnap')
      !ENDIF
    
      IF (.NOT.RESTORET) THEN
!       csw34> Set the centre of mass to be at the specified location
!       contained in the SETCENTRE X Y Z keyword
         IF (SETCENT) CALL SETCENTRE(COORDS)
!
! For MAKEOLIGOT and MAKEOLIGOSTART=TRUE: generate oligomers by placing new segments.
      ENDIF

!  Calculate the initial energy and save in EPREV
!op226>{{{ 
      ! {{{
      WRITE(LFH,'(A)') 'Calculating initial energy'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      DO JP=1,NPAR
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NQTOT=NQTOT+1
         IF (NPAR.GT.1) THEN
            WRITE(LFH,'(A,I2,A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,']Qu ',NQ(JP),' E=',&
     &           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ELSE
            WRITE(LFH,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',&
     &           POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ENDIF

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
!op226>}}} 

      WRITE(LFH,'(A,I10,A)') 'Starting MC run of ',NSTEPS,' steps'
      WRITE(LFH,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

      NSUPERCOUNT=NSUPER

!  Main basin-hopping loop 
! {{{
      DO J1=NDONE+1,NSTEPS 
         ISTEP = J1

         CALL FLUSH(LFH)
         IF (NEWJUMP) RANNJ=DPRAND()
!
!  ********************************* Loop over NPAR parallel runs ******************************
!
         DO JP=1,NPAR 
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
                     WRITE(LFH,'(A,I6,A)') 'Starting a block of ',BLOCK(JP),' rigid body translational moves'
                     TMOVE(JP)=.TRUE.
                     OMOVE(JP)=.FALSE.
                  ELSE 
                     WRITE(LFH,'(A,I6,A)') 'Starting a block of ',BLOCK(JP),' rigid body angular moves'
                     OMOVE(JP)=.TRUE.
                     TMOVE(JP)=.FALSE.
                  ENDIF
               ENDIF
            ENDIF
!
!  MAM1000> The default temperature used for the MC acceptance criterion is the one derived in
!           the initialisation section above.  MCTEMP is passed to the subroutine TRANSITION, where
!           the acceptance/rejection decision is made.  However, individual MC moves can override
!           this temperature by setting MCTEMP to a different value for the current step.
            MCTEMP = TEMP(JP)
! 
!  Jump-moves.
! 
            !IF (JUMPMOVE(JP).AND.(MOD(J1,JUMPINT(JP)).EQ.0)) CALL JUMPM(RANNJ,J1,JP,EPPREV)
! 
!  Ordinary steps.
! 
23          CONTINUE
!
! Don;t call symmetry unless the minimum in the Markov chain has changed.
! We should really check if the minimum has changed since the last call to SYMMETRY,
! which can be done with ABS(ELASTSYM(JP)-EPREV(JP)) if NSYMINTERVAL=1.
!
!           IF (SYMMETRIZE.AND.(MOD(J1-1,NSYMINTERVAL).EQ.0).AND.(ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV)) THEN
!           WRITE(LFH,'(A,3G20.10)') 'ELASTSYM,EPREV,diff=',ELASTSYM(JP),EPREV(JP),ABS(ELASTSYM(JP)-EPREV(JP))
            !IF (SYMMETRIZE.AND.(MOD(NQ(JP),NSYMINTERVAL).EQ.0).AND.&
     !&               ((SYMMETRIZECSM).OR.(ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV))) THEN
               !IF ((ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV)) NSYMCALL=0
               !ELASTSYM(JP)=EPREV(JP)
!!              IF ((.NOT.MOVESHELLT).OR.(NSURFMOVES(JP).LT.0)) THEN
                  !IF (SYMMETRIZECSM) THEN
                     !CALL SYMMETRYCSM(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL,  &
     !&                                J1,NSUCCESS,NFAIL,EBEST,BESTCOORDS,JBEST,EPPREV)
                  !ELSE
                     !CALL SYMMETRY(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL, &
     !&                                J1,EBEST,BESTCOORDS,JBEST,EPPREV)
                  !ENDIF
                  !IF (HIT) GOTO 37 
!!              ELSE
!!                 CALL SYMMETRY2(JP,SCREENC,QDONE,BRUN,ITERATIONS,TIME,CHANGEDE,NSYMCALL)
!!              ENDIF
!!              WRITE(LFH,'(A,I2,A,2I6)') '[',JP,']mc> NCORE: ',NCORE(1:NPAR)
!!              IF (HIT) GOTO 37 ! hit cannot change in symmetry2 
!!
!!  Check for reseeding.
!! 
               !POTEL=EPREV(JP) ! NEWRES assumes POTEL is the energy of the current structure in COORDS
               !IF (CHANGEDE.AND.NEWRESTART) THEN
                  !CALL NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,&
     !&                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST)
               !ENDIF
            !ELSEIF (ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV) THEN ! Markov minimum has changed, but SYMMETRY not called
            IF (ABS(ELASTSYM(JP)-EPREV(JP)).GT.ECONV) THEN ! Markov minimum has changed, but SYMMETRY not called
               NSYMREM=0                                       ! Should therefore reset NSYMREM.
            ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! START OF STEP TAKING CALLS!                                                                            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set switch variables for moves that are not always performed
               !IF(GROUPROTT.AND.MOD(J1,GROUPROTFREQ).EQ.0) THEN
                       !DOGROUPROT=.TRUE.
               !ENDIF
!
! csw34> Coordinates are saved so that moves can be undone
!
                  SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
! csw34> If you want to look at the effect of moves, you can dump out
! the structure BEFORE the move here.
!                 CALL A9DUMPPDB(COORDS(:,JP),"beforemove")
!                 CALL CHARMMDUMP(COORDS(:,JP),'beforemove')

!
!
! ALL OTHER STEP TAKING
!
!  These coordinates are overwritten if we try to call MPI_IPROBE
!
!                 WRITE(LFH,'(I6)') NATOMS
!                 WRITE(LFH,'(A,G20.10)') 'eprev=',EPREV
!                 WRITE(LFH,'(A,3G20.10)') ('LA ',COORDS(3*(J3-1)+1:3*(J3-1)+3,JP),J3=1,NATOMS)
                  CALL TAKESTEP(JP)
!! Restore atom coordinates if atom is FROZEN or DONTMOVE as long as
!! we're not taking internal coordinate moves in CHARMM or AMBER
               !IF ((CHNMAX.LE.0.0D0).AND.(AMCHNMAX.LE.0.0D0)) THEN
                  !IF(FREEZE) THEN
                     !DO J2=1,NATOMS
                        !IF (FROZEN(J2)) THEN
                           !COORDS(3*(J2-1)+1:3*(J2-1)+3,JP)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                        !ENDIF
                     !ENDDO
                  !ENDIF
!!
!! Same for DONTMOVE atoms
!! 
                  !IF(DONTMOVET) THEN
                     !DO J2=1,NATOMS
                        !IF (DONTMOVE(J2)) THEN
                           !COORDS(3*(J2-1)+1:3*(J2-1)+3,JP)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                        !ENDIF
                     !ENDDO
                  !ENDIF
               !ENDIF
 
!
! Reset switch variables for steps not done every time
!
                     dogrouprot=.FALSE.
               
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
                      !CALL A9DUMPPDB(COORDS(:,JP),QUNAME)
                  ENDIF
! END OF KEYWORD <DUMPSTEPS> BLOCK

               NQ(JP)=NQ(JP)+1
               CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)  
               NQTOT=NQTOT+1
 
!  Check for results of taboo list. SELFT is set in taboo.
 
               IF (SELFT) THEN
                  CALL MYRESET(JP,NATOMS,NPAR,NSEED)
                  IF (DEBUG) THEN
                     WRITE(LFH,'(A)') 'Taboo list:'
                     WRITE(LFH,'(6F20.10)') (ESAVE(J2,1),J2=1,NT(1))
                     WRITE(LFH,73) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
73                   FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' TABOO')
                  ENDIF
                  GOTO 23
               ENDIF
   
!
!  Output
!
                  WRITE(LFH,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=',&
     &                 POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART

!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!                 CALL POTENTIAL(COORDSO(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> A energy for coordinates in COORDSO=',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!                 CALL POTENTIAL(COORDS(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> A energy for coordinates in COORDS= ',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!1


!     mp466>  writes structure and energetic data at regular increments
!             to *plot and movie files for AMH potential
               CALL FLUSH(LFH)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!1
!                 CALL POTENTIAL(COORDSO(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> B energy for coordinates in COORDSO=',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!                 CALL POTENTIAL(COORDS(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> B energy for coordinates in COORDS= ',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!1

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
             !IF (RMST.AND.CHRMMT) THEN
                !WRITE(4428,'(I10,F15.5)') NQ(JP),RMSD
                !CALL FLUSH(4428)
             !ENDIF
          ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!                 CALL POTENTIAL(COORDSO(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> C energy for coordinates in COORDSO=',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!                 CALL POTENTIAL(COORDS(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> C energy for coordinates in COORDS= ',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!1

! DAESTAT keyword prints two files :
! stat.all which contains all quenches and their energies
! stat.acc which contains only accepted quenches and their energies
! used to analyse how many new minima are being found
 
            !IF (DAESTAT) THEN
               !PRINT*,'DAESTAT block in mc.f not implemented'
               !STOP
!!              CALL CALCMIND(JP,MIND)
!!              CALL CALCDIHE(DIHE)
!!              WRITE(LFH,'(A,I6,3F20.10)') 'NQALL POTEL',NQ(JP),POTEL,MIND,DIHE
!!              WRITE(36,'(I6,3F20.10)') NQ(JP),POTEL,MIND,DIHE
            !ENDIF
 
!  Check for reseeding.
 
            !IF (NEWRESTART.AND.(.NOT.SEEDT)) THEN 
              !CALL NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,RMIN,RVAT,BRUN,SCREENC,QDONE,&
     !&                    JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DJW
!!             IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,JP))
!!             COORDSO(1:3*(NATOMS-NSEED),JP)=COORDS(1:3*(NATOMS-NSEED),JP)
!!             WRITE(LFH,'(A,2G20.10)'),'mc> coordso changed: ',COORDSO(1,JP),COORDS(1,JP)     
!!             VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DJW
            !ENDIF
            IF (STAY) THEN
            ELSE IF (EVAP .and. .not.evapreject) THEN
               NFAIL(JP)=NFAIL(JP)+1
               CALL MYRESET(JP,NATOMS,NPAR,NSEED)
               IF (DEBUG) THEN
                  WRITE(LFH,33) JP,J1,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
33                FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' EVAP,REJ')
               ENDIF
            ELSE
!     csw34> A series of tests start here to check if a structure should
!     be allowed into the markov chain. If it fails a test, the ATEST
!     variable will be set to .FALSE. 
               ATEST=.TRUE. 
!     csw34> Test for cold fusion in mylbfgs.f - if so, set ATEST to .FALSE.
               IF (COLDFUSION) THEN
                  ATEST=.FALSE.
               ENDIF
               COLDFUSION=.FALSE.

!     jwrm2> If using percolation, and it's not a connected network, reject the step, unless no
!            step has yet been connected.
!            If this is the first step that's fully connected, set PERCACCEPTED to true.
               IF (PERCOLATET) THEN
                 IF(PERCACCEPTED .AND. (.NOT. PERCT)) THEN
                   WRITE(LFH,*) 'perc> Disconnected atoms, rejecting step'
                   ATEST=.FALSE.
                 ELSE IF((.NOT. PERCACCEPTED) .AND. PERCT) THEN
                   PERCACCEPTED=.TRUE.
                   WRITE(LFH,*) 'perc> Found first connected configuration'
                 ENDIF
               ENDIF

!     csw34> If the chirality or peptide bond checks in quench.f have
!            failed, GOODSTRUCTURE should be .FALSE. and we should set ATEST to
!            .FALSE. to prevent these bad structures entering the markov chain!
               IF (.NOT.GOODSTRUCTURE) ATEST=.FALSE.

!     csw34> Dumping every quench in AMBER9 format - only dumps if no chiral problems found 
!     Only for AMBER! DJW Otherwise we get the cis-trans message for every other potential!
               !IF (AMBERT) THEN
!!     csw34> Calculate interaction energy to residue specified in the AMBGMINintE.sh script
!!                  IF (A9INTET.AND.ATEST) THEN
                     !!CALL A9INTE(COORDS(:,JP),INTE)
                     !!CALL A9INTESAVEIT(INTE,COORDS(:,JP),JP)
                     !!CALL A9INTEDUMPLOWEST()
                     !!IF (DEBUG) WRITE(LFH,*) 'intE=',INTE
                     !!IF (TRACKDATAT) THEN
                        !!WRITE(3998,'(I10,G20.10)') NQ(JP),INTE
                        !!CALL FLUSH(3998)
                     !!ENDIF
                  !!ENDIF
!!     csw34> Dump the quench coordinates in AMBER pdb and rst format if DUMPQ is present and chirality checks satisfied
                  !IF (DUMPQUT.AND.ATEST) THEN
                      !WRITE(QUENCHNUM,*) NQ(JP)
                      !QUNAME='quench'//TRIM(ADJUSTL(QUENCHNUM))//'.rst'
                      !OPEN(UNIT=20,FILE=QUNAME,STATUS='UNKNOWN')  
                      !WRITE(20,'(a20)') QUNAME
                      !WRITE(20,'(i5)') NATOMS
                      !WRITE(20,'(6f12.7)') COORDS(:,JP) 
                      !CLOSE(20)
!!     csw34> Dump to PDB using routine in amberinterface.f
                     !QUNAME='quench'//TRIM(ADJUSTL(QUENCHNUM))
                     !CALL A9DUMPPDB(COORDS(:,JP),QUNAME)
                  !ELSEIF (DUMPQUT) THEN
                      !WRITE(*,'(A,I6)') 'DUMPQU> chirality/cis-trans problem detected, not dumping quench ',NQ(JP)
                  !ENDIF
!!     khs26> Dump current best structure to PDB using routine in amberinterface.f
                  !IF (DUMPBESTT) THEN
                     !WRITE(QUENCHNUM,*) NQ(JP)
                     !BESTNAME='best_'//TRIM(ADJUSTL(QUENCHNUM))
                     !CURRENTBESTNAME='current_best'
                     !CALL A9DUMPPDB(QMINP(1,:), BESTNAME)
                     !CALL A9DUMPPDB(QMINP(1,:), CURRENTBESTNAME)
                  !ENDIF
               !ENDIF
!     csw34> Check to see if LOCALSAMPLE constraints have been violated - if they have, reject the step
              ! IF (LOCALSAMPLET) THEN
                  !DISTOK=.FALSE.
                  !CALL A9DISTCHECK(COORDS(:,JP),DISTOK)
                  !IF (.NOT.DISTOK) ATEST=.FALSE.
               !ENDIF

               IF ((QDONE.EQ.0).AND.TIP) THEN
                  ATEST=.FALSE.
               ELSEIF (ATEST) THEN
                  CALL TRANSITION(POTEL,EPREV(JP),ATEST,JP,RANDOM,MCTEMP)
               ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!                 CALL POTENTIAL(COORDSO(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> D energy for coordinates in COORDSO=',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!                 CALL POTENTIAL(COORDS(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> D energy for coordinates in COORDS= ',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  Sanity check to make sure the Markov energy agrees with COORDSO. 
!  Stop if not true.
!
               IF (DEBUG.OR.CHECKMARKOVT) THEN
                  CALL POTENTIAL(COORDSO(:,JP),GRAD,OPOTEL,.FALSE.,.FALSE.)
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
! Accept or reject step. If the quench did not converge then allow a
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
                  CALL MYRESET(JP,NATOMS,NPAR,NSEED)
                  IF (DEBUG) THEN
                     WRITE(LFH,36) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
36                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
               ENDIF
            ENDIF
!           WRITE(LFH,'(A,4F20.10)') 'Q4 values ',QS,QF,QSTART,QFINISH
!
!  Dump coordinates and energy if another run is attempting to jump to this one.
!
            IF ((NPAR.GT.1).AND.(.NOT.NEWJUMP)) CALL DUMPJ(JP,JUMPTO,NPAR,COORDS(1:3*NATOMS,1:NPAR),NATOMS,EPREV)
!
!  If RESTART then reseed if we haven t accepted a step in twice the relaxation time.
!
            IF (RESTART.AND.(J1-JACCPREV.GT.1.1D0*NRELAX)) CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
!
!  Check the acceptance ratio.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!                 CALL POTENTIAL(COORDSO(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> E energy for coordinates in COORDSO=',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!                 CALL POTENTIAL(COORDS(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> E energy for coordinates in COORDS= ',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            IF ((MOD(J1,NACCEPT).EQ.0).AND.(NSEED.EQ.0).AND.(.NOT.STAY)) CALL ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!                 CALL POTENTIAL(COORDSO(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> F energy for coordinates in COORDSO=',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!                 CALL POTENTIAL(COORDS(1:3*NATOMS,JP),DUMGRAD,DJWPOTEL,.FALSE.,.FALSE.)
!                 WRITE(LFH,'(2(A,G20.10))') 'mc> F energy for coordinates in COORDS= ',DJWPOTEL, 
!    &                                                 ' Markov energy=',EPREV(JP) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!1

            TEMP(JP)=TEMP(JP)*SCALEFAC
            IF (HIT) GOTO 37
            IF (DUMPINT.GT.0) THEN
               IF (MOD(J1,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(J1,EBEST,BESTCOORDS,JBEST,JP)
               ENDIF
            ENDIF
           IF (NQ(JP).GT.NSTEPS) GOTO 37

         ENDDO
!
!  ****************************** End of loop over NPAR parallel runs *****************************
!
!  Supersteps
!
!        IF (SUPERSTEP) CALL SUPERMC(SPOTEL,SCOORDS,NSUPERCOUNT,POTEL)
!
!  Renormalised type steps
!
!        IF (RENORM) CALL REN(J1,RMIN,RCOORDS,RVAT,NREN,RMINO,RCOORDSO,RVATO,ITERATIONS,TIME,NLAST,JACCPREV,NSTEPREN)
!        IF (STAY) CALL MYRESET(1,NATOMS,NPAR,NSEED)
  
         CALL FLUSH(LFH)
      ENDDO
! }}}

37    CONTINUE
      DO JP=1,NPAR
         IF (NPAR.GT.1) THEN
            WRITE(LFH,20) JP,NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),&
     &               STEP(JP),ASTEP(JP),TEMP(JP)
20          FORMAT('[',I2,']Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
         ELSE
            WRITE(LFH,21) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),&
     &               STEP(JP),ASTEP(JP),TEMP(JP)
21          FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
         ENDIF
         IF (RIGID) WRITE(LFH,'(A,F12.5)') 'Rigid body orientational step=',OSTEP(JP)
      ENDDO
!mo361>Deallocating these arrays to cope with multiple runs of this subroutine in GA
      DEALLOCATE(TMOVE)
      DEALLOCATE(OMOVE)
      RETURN
!op226>}}} 
      END SUBROUTINE MC 


