
! A9INTEDUMPLOWEST {{{
!
!
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
!    csw34> Dump to PDB using routine in amberinterface.f
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
      ! }}}
! A9DISTCHECK {{{
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


      SUBROUTINE JUMPM(RANNJ,J1,JP,EPPREV)
      USE commons
      IMPLICIT NONE
      INTEGER J1, JP, J2, UNT, NDUM, ITERATIONS, BRUN, QDONE
      DOUBLE PRECISION RANNJ, RANDOM, DPRAND, EPPREV(NPAR), DUMMY, TIME, EJ, SCREENC(3*NATOMS)

      IF (NEWJUMP) THEN
         IF (RANNJ.LT.PNEWJUMP) THEN
            RANDOM=DPRAND()
            IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
!           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))/TEMP(JP)).GT.RANDOM) THEN
!           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
               WRITE(LFH,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,&
     &                ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' accepted before quench ',J1
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
               WRITE(LFH,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,&
     &                  ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' rejected before quench ',J1
            ENDIF
         ENDIF
      ELSE
         UNT=70+JUMPTO(JP)
         REWIND(UNT)           
         NDUM=INT(DPRAND()*(NQ(JUMPTO(JP))-1))
         IF (DEBUG) WRITE(LFH,'(A, G20.10)') 'Should be choosing buffer energy number ',NDUM
         DO J2=1,NDUM
            READ(UNT,*) EJ
         ENDDO
         DO J2=1,NQ(JUMPTO(JP))-1-NDUM
            READ(UNT,*) DUMMY
         ENDDO
         RANDOM=DPRAND()
!
!  Coordinates are only read if the jump is successful.
!
         IF (DEXP((EPREV(JP)-EJ)*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
            WRITE(LFH,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP, &
     &              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' accepted before quench ',NQ(JP)
            EPREV(JP)=EJ
            UNT=70+NPAR+JUMPTO(JP)
            REWIND(UNT)           
            DO J2=1,(NDUM-1)*NATOMS
               READ(UNT,*) DUMMY
            ENDDO
            READ(UNT,*) (COORDS(J2,JP),J2=1,3*NATOMS)
!
!  Coordinates should be converged already, but we need to reset VAT and VATO.
!
!  next line should be uncommented if routine is made availabe to use with CHARMM
!            IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
            CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            DO J2=1,NATOMS
               VATO(J2,JP)=VAT(J2,JP)
            ENDDO
            DO J2=1,3*NATOMS
               COORDSO(J2,JP)=COORDS(J2,JP)
            ENDDO
            IF (DEBUG) THEN
               WRITE(LFH,'(A)') 'Jump coordinates:'
               WRITE(LFH,'(3F20.10)') (COORDS(J2,JP),J2=1,3*NATOMS)
            ENDIF
            DO J2=1,NATOMS*(NQ(JUMPTO(JP))-1)-NDUM*NATOMS
               READ(UNT,*) DUMMY
            ENDDO
         ELSE
            WRITE(LFH,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,&
     &              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' rejected before quench ',NQ(JP)
         ENDIF
      ENDIF

      RETURN
      END
      ! }}}
      ! DUMPJ {{{
      SUBROUTINE DUMPJ(JP,JUMPTO,NPAR,COORDS,NATOMS,EPREV)
      IMPLICIT NONE
      LOGICAL TEST
      INTEGER NPAR, J2, JP, JUMPTO(NPAR), NATOMS, UNT
      DOUBLE PRECISION COORDS(3*NATOMS,NPAR), EPREV(NPAR)

      DO J2=1,NPAR
         IF (JUMPTO(J2).EQ.JP) TEST=.TRUE.
      ENDDO
!     TEST=.FALSE.
!     IF (TEST) THEN
!        UNT=70+JP
!        WRITE(UNT,'(F20.10)') EPREV(JP)
!        UNT=70+NPAR+JP
!        WRITE(UNT,'(3F20.10)') (COORDS(J2,JP),J2=1,3*NATOMS)
!     ENDIF

      RETURN
      END
      ! }}}
      ! MYRESET {{{

      SUBROUTINE MYRESET(JP,NATOMS,NPAR,NSEED)
      USE COMMONS,ONLY : LFH,COORDS,COORDSO,VAT,VATO
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
      ! }}}
      ! REST {{{
      SUBROUTINE REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)

      USE COMMONS

      IMPLICIT NONE
      ! sub 
      INTEGER ::    ITERATIONS
      DOUBLE PRECISION ::   TIME
      INTEGER :: J1
      DOUBLE PRECISION,DIMENSION(3*NATOMS) ::   RCOORDS
      DOUBLE PRECISION ::   RMIN
      DOUBLE PRECISION,DIMENSION(NATOMS) ::   RVAT
      INTEGER JACCPREV
      ! local 
      INTEGER ITERATIONS, J2, JACCPREV, J1, NQTOT, BRUN, QDONE
      DOUBLE PRECISION POTEL, SCREENC(3*NATOMS)
      ! common
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT

10    CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),1,NHSRESTART)
!  next line should be uncommented if routine is made availabe to use with CHARMM
!      IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
      CALL QUENCH(.FALSE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
!
!  Bad idea to accept this quench configuration unconditionally - it could be unconvergeable.
!
      IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
         DO J2=1,3*NATOMS
            COORDS(J2,1)=COORDSO(J2,1)
         ENDDO
         GOTO 10
      ENDIF
      JACCPREV=J1
      NQTOT=NQTOT+1
      WRITE(LFH,'(A,I6,A)') ' Restarting using ',NHSRESTART,' hard sphere moves'
      WRITE(LFH,'(A,I7,A,F20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Restart Qu ',NQ(1),' E=',&
     &              POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
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
      ! }}}
      ! REN {{{
      SUBROUTINE REN(J1,RMIN,RCOORDS,RVAT,NREN,RMINO,RCOORDSO,RVATO,ITERATIONS,TIME,NLAST,JACCPREV,NSTEPREN)
      USE COMMONS

      IMPLICIT NONE
      ! sub
      INTEGER ::    J1,JACCPREV, NSTEPREN
      INTEGER NREN,ITERATIONS,NLAST
      DOUBLE PRECISION ::   RMIN,RMINO,TIME
      DOUBLE PRECISION, DIMENSION(3*NATOMS) ::   RCOORDS,RCOORDSO
      DOUBLE PRECISION, DIMENSION(NATOMS) ::   RVAT,RVATO
      ! loc
      LOGICAL STAY, REJECT, METROPOLIS
      INTEGER J2,  NQTOT,  J3, BRUN, QDONE
      DOUBLE PRECISION POTEL,  RANDOM, DPRAND 
      DOUBLE PRECISION ::  XIP, DUMMY, SCREENC(3*NATOMS)
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
      IF (DEBUG) WRITE(LFH,'(A,2G20.10)') 'RMIN,POTEL=',RMIN,POTEL
!     PRINT*,'J1,JACCPREV+NRENSTUCK,NLAST+NREN=',J1,JACCPREV+NRENSTUCK,NLAST+NREN
      IF ((J1.GE.JACCPREV+NRENSTUCK).OR.(J1.GE.NLAST+NREN)) THEN
         JACCPREV=J1
         NLAST=J1
         RANDOM=DPRAND()
         METROPOLIS=DEXP(-(RMIN-RMINO)/TRENORM).GT.RANDOM
         REJECT=.FALSE.
!
!  Taboo list for renormalised energies. Skip if the step is going to be rejected by
!  the Metropolis condition.
!
         IF (TABOOT.AND.METROPOLIS) THEN
            IF (NSTEPREN.EQ.0) NT(1)=0
            CALL NEWINERTIA(RCOORDS,NATOMS,NATOMS,XIP)
            DO J1=1,NT(1)
               IF (DABS(RMIN-ESAVE(J1,1)).LT.ECONV) THEN
                  IF (DABS(XIP-XINSAVE(J1,1))/(XIP+XINSAVE(J1,1)).LT.1.0D-2) THEN
                     REJECT=.TRUE.
                     GOTO 20
                  ELSE
                     WRITE(LFH,'(A, 2G20.10)') 'Energies nearly degenerate:',RMIN,ESAVE(J1,1)
                     WRITE(LFH,'(A, 2G20.10)') 'But  different  structures:',XIP,XINSAVE(J1,1)
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

            WRITE(LFH,'(A,I10)') ' Number of entries in taboo list=',NT(1)
            IF (DEBUG) THEN
               WRITE(LFH,'(6F20.10)') (ESAVE(J2,1),J2=1,NT(1))
            ENDIF
         ENDIF
!
!  Accept/reject for renormalised step
!
         IF (METROPOLIS.AND.(.NOT.REJECT)) THEN
            IF (NSTEPREN.GT.0) WRITE(LFH,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' accepted'
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
               WRITE(LFH,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected by taboo criterion'
            ELSE
               WRITE(LFH,'(A,G20.10,A,G20.10,A)') ' renorm step from ',RMINO,' to ',RMIN,' rejected'
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
         IF (NSTEPREN.EQ.1) WRITE(LFH,'(A,G20.10)') ' first renorm energy is ',RMIN
         WRITE(LFH,'(A,I6)') ' renormalisation interval is now ',NREN
!
!  Longer renorm step
!
10       IF ((XMOVERENORM.GT.3.0D0).OR.FIXD) THEN
            CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),1,INT(XMOVERENORM))
         ELSE
            DUMMY=STEP(1)
            STEP(1)=XMOVERENORM
            CALL TAKESTEP(1)
            STEP(1)=DUMMY
         ENDIF
!  next line should be uncommented if routine is made availabe to use with CHARMM
!         IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
         CALL QUENCH(.FALSE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
!
!  Bad idea to accept this quench configuration unconditionally - it could be unconvergeable.
!
         IF (POTEL-EPREV(1).GT.10.0D0*ABS(EPREV(1))) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,1)=COORDSO(J2,1)
            ENDDO
            GOTO 10
         ENDIF
         NQTOT=NQTOT+1
         WRITE(LFH,'(A,I7,A,F20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Renorm Qu ',NQ(1),' E=',&
     &        POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
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
      ! }}}
!  NEWRES {{{ 
!  Reseed if the energy has not improved by more than ECONV over the
!  last NRELAX mc steps.
!  If AVOID is true then save the energy and coordinates of the lowest
!  minimum achieved before each restart and restart if we get too close
!  to any one of them using bipartite matching and mind. Note that bipartite
!  matching and mind can give a local minimum of distance if the optimal 
!  permutation-inversion isn;t found. Using ORIENT to effect a standard
!  orientation first seems to help. It should at least ensure that permutation-inversion
!  isomers are always found. 
!  Would it be possible to just use the energy and inertia components to identify
!  permutation-inversion isomers instead of MINPERM and MIND?
!  This would be like using a small value of AVOIDDIST, which doesn;t seem to be
!  as good.
!
      SUBROUTINE NEWRES(J1,JP,JBEST,EBEST,BESTCOORDS,EPPREV,POTEL,ITERATIONS,TIME,RCOORDS,&
     &                  RMIN,RVAT,BRUN,SCREENC,QDONE,JACCPREV,NSUCCESS,NFAIL,NFAILT,NSUCCESST)
      USE commons
      USE modamber9, only : mdstept 
      
      IMPLICIT NONE
      INTEGER J1, JP, JBEST(NPAR), ITERATIONS, J2, JACCPREV, BRUN, QDONE, J3, PERM(NATOMS), NPERM, NTRIES
      INTEGER, PARAMETER :: MAXIMUMTRIES=20
      DOUBLE PRECISION EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR)
      DOUBLE PRECISION :: EPPREV(NPAR), POTEL, TIME, RCOORDS(3*NATOMS), DIST2, DUMMY2
      DOUBLE PRECISION :: RVAT(NATOMS), RMIN, RANDOM, SR3, SCREENC(3*NATOMS)
      DOUBLE PRECISION :: DPRAND, FCOORDS(3*NATOMS), XMSBSAVE(3*NATOMS)
      DOUBLE PRECISION :: DUMMY(3*NATOMS), DISTANCE, XMSB(3*NATOMS)
      DOUBLE PRECISION :: EBESTP, BESTCOORDSP(3*NATOMS), WORSTRAD, RMAT(3,3), QENERGY
      DOUBLE PRECISION :: ROTA(3,3), ROTINVA(3,3)
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
!
!  Reseed if the energy has not improved in the last NRELAX mc cycles,
!  or if the current minimum is too close to one of the NMSBSAVE structures
!  saved in MSBCOORDS.
!
!  Instead of using the current minimum, employ the current minimum in 
!  BESTCOORDS for the AVOID check. Then we need only do the AVOID check when we have
!  a new minimum in BESTCOORDS, i.e. J1.EQ.JBEST(JP).
!
      RES1=.FALSE.
      IF (J1-JBEST(JP).GT.NRELAX) RES1=.TRUE.
!     WRITE(LFH,'(A,I5,2G17.7,3I5,L8)') 'J1,POTEL,EBEST,JBEST,J1-JBEST,NRELAX,RES1=',
!    1                                J1,POTEL,EBEST(JP),JBEST(JP),J1-JBEST(JP),NRELAX,RES1
      RES2=.FALSE.
!     IF ((.NOT.RES1).AND.AVOID) THEN
!     PRINT*,'RES1,AVOID,J1,JBEST(JP)=',RES1,AVOID,J1,JBEST(JP)
      IF ((.NOT.RES1).AND.AVOID.AND.(J1.EQ.JBEST(JP)).AND.(NMSBSAVE.GT.0)) THEN ! best minimum has just changed.
         FCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
         savedloop: DO J2=1,MIN(NMSBSAVE,MAXSAVE)
!
!  Bipartite matching routine for permutations. Coordinates in FCOORDS do not change
!  but the coordinates in XMSB do. DISTANCE is the distance in this case.
!
! If the energy is lower no reseeding regardless of separation ? DJW
!           PRINT '(A,2G20.10,L8)','POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV=',
!    1                              POTEL,MSBE(J2)-ECONV,POTEL.LT.MSBE(J2)-ECONV
            IF (POTEL.LT.MSBE(J2)-ECONV) CYCLE savedloop
            XMSB(1:3*NATOMS)=MSBCOORDS(1:3*NATOMS,J2)
            CALL MINPERMDIST(FCOORDS,XMSB,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DISTANCE,DIST2,RIGID,RMAT)
            IF (DISTANCE.LT.AVOIDDIST) THEN
               RES2=.TRUE.
               WRITE(LFH,'(A,G20.10,A,I6,A,G20.10,A,F10.3)') 'newres> Minimum energy ',POTEL,&
     &                                    ' is too close to saved structure ',&
     &                                 J2,' with energy ',MSBE(J2),' dist=',DISTANCE
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

         IF (RES2.AND.(.NOT.AVOIDRESEEDT)) THEN 
            POTEL=MAX(1.0D10,10.0D0*POTEL) ! This should be enough to reject the step until we do a massive Thomson problem!
            WRITE(LFH,'(A,I8,A)') 'newres> Resetting energy to a large value due to taboo condition'
            RETURN
         ELSEIF (RES1.OR.(RES2.AND.(DISTANCE.GT.0.1D0))) THEN ! new condition
            NMSBSAVE=NMSBSAVE+1
            NPOSITION=MOD(NMSBSAVE,MAXSAVE)
            IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
            MSBE(NPOSITION)=EBEST(JP)
            FCOORDS(1:3*NATOMS)=BESTCOORDS(1:3*NATOMS,JP)

            CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)

            MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)

            WRITE(LFH,'(A,I6,A,G20.10)') 'newres> Moving current best minimum to position ',NPOSITION,&
     &                         ' in the AVOID list E=',EBEST(JP)
!           OPEN(UNIT=34,FILE='MSBdata',POSITION='APPEND')
!           WRITE(34,'(G20.10)') MSBE(NPOSITION) 
!           WRITE(34,'(3G20.10)') BESTCOORDS(1:3*NATOMS,JP)
!           CLOSE(34)
         ELSEIF (RES2.AND.(DISTANCE.LE.0.01D0)) THEN
!           IF (NMSBSAVE.LT.MAXSAVE) THEN
               NMSBSAVE=NMSBSAVE+1
               NPOSITION=MOD(NMSBSAVE,MAXSAVE)
               IF (NPOSITION.EQ.0) NPOSITION=MAXSAVE
               MSBE(NPOSITION)=EBESTP
               FCOORDS(1:3*NATOMS)=BESTCOORDSP(1:3*NATOMS)
               CALL MYORIENT(FCOORDS,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)
               MSBCOORDS(1:3*NATOMS,NPOSITION)=DUMMY(1:3*NATOMS)
               WRITE(LFH,'(A,I6,A,G20.10)') 'newres> Moving previous best minimum to position ',NPOSITION,&
     &                             ' in the AVOID list E=',EBESTP
!           ENDIF
         ENDIF ! end new condition

         IF (NEWRESTART_MD) THEN ! lb415
            IF (RES1) WRITE(LFH,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
            WRITE(LFH,'(A,I8,A)') 'newres> Reseeding via a short high temperature MD run'
            CHANGE_TEMP = .true.

! th368: 20-10-2009 Extending MD-Reseeding to the CHARMM interface
! terminate if neither AMBER nor CHARMM was requested in the data file

             IF (AMBERT) THEN
               CALL TAKESTEPAMBER(JP,COORDS(:,JP),movableatomlist,nmovableatoms,ligmovet,mdstept,randomseedt)
             ELSEIF (CHRMMT) THEN
               CALL CHMD(JP)
             ELSE
               WRITE(LFH,'(A,I8,A)') 'newres> Molecular Dynamics reseeding is available for AMBER or CHARMM runs only.'
               STOP
             ENDIF
! end th368: 20-10-2009

            CHANGE_TEMP = .false.
         ELSE
            IF (NHSRESTART.GT.0) THEN
               IF (RES1) WRITE(LFH,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - perturbing'
               IF (RES2) WRITE(LFH,'(A,I8,A)') 'newres> Reseeding due to taboo condition'
               IF (SHELLMOVES(JP)) WRITE(LFH,'(A)') 'newres> Turning off shell moves'
               SHELLMOVES(JP)=.FALSE.
               CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
            ELSE
               IF (RES1) THEN
                  IF (AVOIDRESEEDT) THEN
                     WRITE(LFH,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP),' - reseeding'
                  ELSE
                     WRITE(LFH,'(A,I8,A)') 'newres> Energy has not improved since step ',JBEST(JP)
                  ENDIF
               ENDIF
               IF (RES2) WRITE(LFH,'(A,I8,A)') 'newres> Reseeding due to taboo condition'
               HIGHEST=.TRUE.
               DO J2=1,NPAR
                  IF (J2.EQ.JP) CYCLE
                  IF (EBEST(J2).GT.EBEST(JP)) THEN
                     HIGHEST=.FALSE.
                     EXIT
                  ENDIF
               ENDDO
               IF (HIGHEST) THEN
                  WRITE(LFH,'(A,I6,A,F20.10)') 'newres> Parallel run ',JP,' has the highest energy ',EBEST(JP)
                  IF (NPAR.GT.1) WRITE(LFH,'(6F20.10)') EBEST(1:NPAR)
                  IF (.NOT.AVOIDRESEEDT) THEN
                     WRITE(LFH,'(A,I8,A)') 'newres> Resetting energy to a large value to reject this step (not reseeding)'
                     JBEST(JP)=J1
                     EBEST(JP)=POTEL ! this is communicated via common block MYPOT, which is in quench.f
                     POTEL=MAX(1.0D10,10.0D0*POTEL) ! This should be enough to reject the step until we do a massive Thomson problem!
                     RETURN
                  ELSE
                     WRITE(LFH,'(A)') 'newres> Full reseeding'
                     DO J2=1,3*NATOMS
                        RANDOM=(DPRAND()-0.5D0)*2.0D0
                        COORDS(J2,JP)=RANDOM*DSQRT(RADIUS)/SR3
                     ENDDO
                     NCORE(JP)=0
                     PTGROUP(JP)='   '
                  ENDIF
               ELSEIF (NCORE(JP).GT.0) THEN
!         
!  Reseed everything if this is the lowest of a set of parallel runs, otherwise
!  just reseed the surface.
!         
                  WRITE(LFH,'(A)') 'newres> Accepting an enforced surface reseeding'
                  DUMMY2=-1.0D0
                  DO J2=NATOMS-NCORE(JP)+1, NATOMS
                     DISTANCE=COORDS(3*(J2-1)+1,JP)**2+COORDS(3*(J2-1)+2,JP)**2+COORDS(3*(J2-1)+3,JP)**2
                     IF (DISTANCE.GT.DUMMY2) DUMMY2=DISTANCE
                  ENDDO
                  DISTANCE=SQRT(DISTANCE)
                  WRITE(LFH,'(A,F15.5)') 'newres> largest core radius=',DISTANCE
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
!
!  This step will be accepted because EPREV(JP)=POTEL. However, there
!  could be additional steps, so we need
!  to change COORDSO and VATO. Should reset EBEST and BESTCOORDS as well.
!
!  next line should be uncommented if routine is made available to use with CHARMM
!         IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
654      CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         IF (RMS.GT.BQMAX) THEN
            WRITE(LFH,'(A)') 'newres> Quench from reseeded geometry failed - try again'
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
!        WRITE(LFH,'(A,2G20.10)') 'newres> coordso changed: ',COORDSO(1,JP),COORDS(1,JP)
         VATO(1:NATOMS,JP)=VAT(1:NATOMS,JP)

! end th368 20.10.2009
      ENDIF

      RETURN
      END
      ! }}}
