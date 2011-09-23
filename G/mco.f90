
!! A9INTEDUMPLOWEST {{{
!!
!!
!! csw34> Subroutine to dump the current lowest interaction energy structure 
!!> \brief Subroutine to dump the current lowest interaction energy structure 
!!> \author Chris Whittleston, csw34@cam.ac.uk
!!
      !SUBROUTINE A9INTEDUMPLOWEST()
      !USE COMMONS
      !USE QMODULE
      !IMPLICIT NONE
      !OPEN(UNIT=20,FILE="bestintE.rst",STATUS='UNKNOWN')  
      !WRITE(20,'(g20.10,i5)') INTEQMIN(1), INTEFF(1)
      !WRITE(20,'(i5)') NATOMS
      !WRITE(20,'(6f12.7)') INTEQMINP(1,:) 
      !CLOSE(20)
!!    csw34> Dump to PDB using routine in amberinterface.f
      !CALL A9DUMPPDB(INTEQMINP(1,:),"bestintE")
      !END SUBROUTINE A9INTEDUMPLOWEST
!! csw34> Subroutine to work out the interaction energy INTE between residue RESLIG for geometry INTECOORDS
      !SUBROUTINE A9INTE(INTECOORDS,INTE)
      !USE PORFUNCS
      !USE COMMONS, ONLY : NATOMS
      !IMPLICIT NONE
      !DOUBLE PRECISION :: INTECOORDS(3*NATOMS), INTE
      !INTEGER ISTAT
!! Dump current coordinates to file for SANDER to read in
      !OPEN(UNIT=20,FILE='coords.intres',STATUS='UNKNOWN')
      !WRITE(20,'(a20)') 'intE coordinates'
      !WRITE(20,'(i5)') NATOMS
      !WRITE(20,'(6f12.7)') INTECOORDS(:)
      !CLOSE(20)
!! Run the script that does the interaction energy calculation (downloadable from the group website)
      !CALL SYSTEM_SUBR('bash AMBGMINintE.sh',ISTAT)
!! Read interaction energy in from the file called intE
      !OPEN(UNIT=20,FILE="intE",STATUS='OLD')
      !READ(20,*) INTE
      !CLOSE(20)
      !END SUBROUTINE A9INTE
      !! }}}
!! A9DISTCHECK {{{
!!> \brief Subroutine to check the  distance 
!!> \brief between groups of atoms defined in the movableatoms file.
!!> Checks the A->B and A->C distances,
!!> and if they are greater than the A(B/C)THRESH values defined
!!> in the data file, the routine returns DISTOK=.FALSE.
!!> \author Chris Whittleston, csw34@cam.ac.uk
!!
      !SUBROUTINE A9DISTCHECK(COORDS,DISTOK)
      !USE MODAMBER9, ONLY : NATOMSINA,NATOMSINB,NATOMSINC,ATOMSINALIST,ATOMSINBLIST,ATOMSINCLIST
      !USE COMMONS, ONLY : NATOMS,ABTHRESH,ACTHRESH,DEBUG
      !IMPLICIT NONE
      !LOGICAL :: DISTOK
      !INTEGER :: I,J
      !DOUBLE PRECISION :: COORDS(3*NATOMS)
      !DOUBLE PRECISION :: CENTREOFA(3),CENTREOFB(3),CENTREOFC(3)  
      !DOUBLE PRECISION :: ABDIST,ACDIST 
!! initialise variables
      !CENTREOFA(:)=0.0D0
      !CENTREOFB(:)=0.0D0
      !CENTREOFC(:)=0.0D0
      !ABDIST=0.0D0
      !ACDIST=0.0D0
!! find centre of ligand (group A)
      !DO I=1,NATOMSINA
         !J=ATOMSINALIST(I)
         !CENTREOFA(1)=CENTREOFA(1)+COORDS(3*J-2)
         !CENTREOFA(2)=CENTREOFA(2)+COORDS(3*J-1)
         !CENTREOFA(3)=CENTREOFA(3)+COORDS(3*J  )
      !END DO
      !CENTREOFA(:) = CENTREOFA(:)/NATOMSINA
!! find centre of group B 
      !DO I=1,NATOMSINB
         !J=ATOMSINBLIST(I)
         !CENTREOFB(1)=CENTREOFB(1)+COORDS(3*J-2)
         !CENTREOFB(2)=CENTREOFB(2)+COORDS(3*J-1)
         !CENTREOFB(3)=CENTREOFB(3)+COORDS(3*J  )
      !END DO
      !CENTREOFB(:) = CENTREOFB(:)/NATOMSINB
!! find centre of group B 
      !DO I=1,NATOMSINC
         !J=ATOMSINCLIST(I)
         !CENTREOFC(1)=CENTREOFC(1)+COORDS(3*J-2)
         !CENTREOFC(2)=CENTREOFC(2)+COORDS(3*J-1)
         !CENTREOFC(3)=CENTREOFC(3)+COORDS(3*J  )
      !END DO
      !CENTREOFC(:) = CENTREOFC(:)/NATOMSINC
!! calculate A->B distance
      !ABDIST=SQRT((CENTREOFA(1)-CENTREOFB(1))**2+(CENTREOFA(2)-CENTREOFB(2))**2+(CENTREOFA(3)-CENTREOFB(3))**2)
!! calculate A->C distance
      !ACDIST=SQRT((CENTREOFA(1)-CENTREOFC(1))**2+(CENTREOFA(2)-CENTREOFC(2))**2+(CENTREOFA(3)-CENTREOFC(3))**2)
!! some DEBUG printing
      !IF (DEBUG) THEN
         !WRITE(*,*) 'AB distance=',ABDIST
         !WRITE(*,*) 'AC distance=',ACDIST
         !IF (ABDIST.LT.ABTHRESH) THEN
            !WRITE(*,*) 'A->B condition met! :)'
         !ELSE 
            !WRITE(*,*) 'A->B condition broken :('
         !ENDIF   
         !IF (ACDIST.LT.ACTHRESH) THEN
            !WRITE(*,*) 'A->C condition met! :)'
         !ELSE 
            !WRITE(*,*) 'A->C condition broken :('
         !ENDIF   
      !ENDIF
!! do the check for both conditions
      !IF ((ABDIST.LT.ABTHRESH).AND.(ACDIST.LT.ACTHRESH)) DISTOK=.TRUE.
!! more debug printing
      !IF (DEBUG) WRITE(*,*) 'DISTOK=',DISTOK
 
      !END SUBROUTINE A9DISTCHECK


      !SUBROUTINE JUMPM(RANNJ,J1,JP,EPPREV)
      !USE commons
      !IMPLICIT NONE
      !INTEGER J1, JP, J2, UNT, NDUM, ITERATIONS, BRUN, QDONE
      !DOUBLE PRECISION RANNJ, RANDOM, DPRAND, EPPREV(NPAR), DUMMY, TIME, EJ, SCREENC(3*NATOMS)

      !IF (NEWJUMP) THEN
         !IF (RANNJ.LT.PNEWJUMP) THEN
            !RANDOM=DPRAND()
            !IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
!!           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))/TEMP(JP)).GT.RANDOM) THEN
!!           IF (DEXP((EPREV(JP)-EPREV(JUMPTO(JP)))*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
               !WRITE(LFH,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,&
     !&                ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' accepted before quench ',J1
               !DUMMY=EPREV(JP)
               !EPREV(JP)=EPREV(JUMPTO(JP))
               !EPREV(JUMPTO(JP))=DUMMY
               !DUMMY=EPPREV(JP)
               !EPPREV(JP)=EPPREV(JUMPTO(JP))
               !EPPREV(JUMPTO(JP))=DUMMY
               !DO J2=1,NATOMS
                  !DUMMY=VATO(J2,JP)
                  !VATO(J2,JP)=VATO(J2,JUMPTO(JP))
                  !VATO(J2,JUMPTO(JP))=DUMMY
                  !DUMMY=VAT(J2,JP)
                  !VAT(J2,JP)=VAT(J2,JUMPTO(JP))
                  !VAT(J2,JUMPTO(JP))=DUMMY
               !ENDDO
               !DO J2=1,3*NATOMS
                  !DUMMY=COORDS(J2,JP)
                  !COORDS(J2,JP)=COORDS(J2,JUMPTO(JP))
                  !COORDS(J2,JUMPTO(JP))=DUMMY
                  !DUMMY=COORDSO(J2,JP)
                  !COORDSO(J2,JP)=COORDSO(J2,JUMPTO(JP))
                  !COORDSO(J2,JUMPTO(JP))=DUMMY
               !ENDDO
            !ELSE
               !WRITE(LFH,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,&
     !&                  ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EPREV(JUMPTO(JP)),' rejected before quench ',J1
            !ENDIF
         !ENDIF
      !ELSE
         !UNT=70+JUMPTO(JP)
         !REWIND(UNT)           
         !NDUM=INT(DPRAND()*(NQ(JUMPTO(JP))-1))
         !IF (DEBUG) WRITE(LFH,'(A, G20.10)') 'Should be choosing buffer energy number ',NDUM
         !DO J2=1,NDUM
            !READ(UNT,*) EJ
         !ENDDO
         !DO J2=1,NQ(JUMPTO(JP))-1-NDUM
            !READ(UNT,*) DUMMY
         !ENDDO
         !RANDOM=DPRAND()
!!
!!  Coordinates are only read if the jump is successful.
!!
         !IF (DEXP((EPREV(JP)-EJ)*(1.0D0/TEMP(JP)-1.0D0/TEMP(JUMPTO(JP)))).GT.RANDOM) THEN
            !WRITE(LFH,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP, &
     !&              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' accepted before quench ',NQ(JP)
            !EPREV(JP)=EJ
            !UNT=70+NPAR+JUMPTO(JP)
            !REWIND(UNT)           
            !DO J2=1,(NDUM-1)*NATOMS
               !READ(UNT,*) DUMMY
            !ENDDO
            !READ(UNT,*) (COORDS(J2,JP),J2=1,3*NATOMS)
!!
!!  Coordinates should be converged already, but we need to reset VAT and VATO.
!!
!!  next line should be uncommented if routine is made availabe to use with CHARMM
!!            IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
            !CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            !DO J2=1,NATOMS
               !VATO(J2,JP)=VAT(J2,JP)
            !ENDDO
            !DO J2=1,3*NATOMS
               !COORDSO(J2,JP)=COORDS(J2,JP)
            !ENDDO
            !IF (DEBUG) THEN
               !WRITE(LFH,'(A)') 'Jump coordinates:'
               !WRITE(LFH,'(3F20.10)') (COORDS(J2,JP),J2=1,3*NATOMS)
            !ENDIF
            !DO J2=1,NATOMS*(NQ(JUMPTO(JP))-1)-NDUM*NATOMS
               !READ(UNT,*) DUMMY
            !ENDDO
         !ELSE
            !WRITE(LFH,'(A,I2,A,F20.10,A,I2,A,F20.10,A,I6)') 'Jump move from parallel run ',JP,&
     !&              ' energy ',EPREV(JP),' to run ',JUMPTO(JP),' energy ',EJ,' rejected before quench ',NQ(JP)
         !ENDIF
      !ENDIF

      !RETURN
      !END
      !! }}}
      ! DUMPJ {{{
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

