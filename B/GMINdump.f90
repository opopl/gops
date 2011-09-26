
SUBROUTINE DUMPSTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)

USE COMMONS
USE V
USE F

IMPLICIT NONE

INTEGER,INTENT(IN) :: NDONE, JBEST(NPAR), JP
DOUBLE PRECISION,INTENT(IN) ::  EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR)

! loc
INTEGER J1, MYUNIT2
DOUBLE PRECISION  DUMGRAD(3*NATOMS), OPOTEL, SAVEP
CHARACTER(LEN=20) :: ISTR
LOGICAL LOPEN
DOUBLE PRECISION POTEL
COMMON /MYPOT/ POTEL

MYUNIT2=1
CALL SYSTEM('cp GMIN.dump GMIN.dump.save')
OPEN(MYUNIT2,FILE='GMIN.dump',STATUS='UNKNOWN')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  SAVEP=POTEL
!                  CALL POTENTIAL(COORDSO(1:3*NATOMS,JP),DUMGRAD,OPOTEL,.FALSE.,.FALSE.)
!                  WRITE(LFH,'(3(A,G20.10))') 'dumpstate> energy for coordinates in COORDSO=',OPOTEL, &
!     &                                                 ' Markov energy=',EPREV(JP),' potel=',SAVEP
!                  CALL POTENTIAL(COORDS(1:3*NATOMS,JP),DUMGRAD,OPOTEL,.FALSE.,.FALSE.)
!                  WRITE(LFH,'(3(A,G20.10))') 'dumpstate> energy for coordinates in COORDS= ',OPOTEL, &
!     &                                                 ' Markov energy=',EPREV(JP),' potel=',SAVEP
!                  POTEL=SAVEP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  WRITE(MYUNIT2, '(A)') 'steps completed J1 in mc'
  WRITE(MYUNIT2, '(I8)') NDONE
  WRITE(MYUNIT2, '(I8)') NSAVE
  WRITE(MYUNIT2, '(A)') 'COORDS'
  WRITE(MYUNIT2, '(A,I8)') 'run number ',JP
  CALL WRITE_MARKOV_COORDS(MYUNIT2, '(3F25.15)', JP)
!  WRITE(MYUNIT2, '(3F25.15)') COORDSO(1:3*NATOMS,JP)
  WRITE(MYUNIT2, '(A)') 'STEP, ASTEP, OSTEP, TEMP:'
  WRITE(MYUNIT2, '(4F25.15)') STEP(JP),ASTEP(JP),OSTEP(JP),TEMP(JP)
  WRITE(MYUNIT2, '(A)') 'QMIN and QMINP'
! best NSAVE minima
  DO J1=1,NSAVE
     WRITE(MYUNIT2, '(A,I8)') 'saved minimum ',J1
     WRITE(MYUNIT2, '(G25.15)') QMIN(J1)
     WRITE(MYUNIT2, '(3F25.15)') QMINP(J1,1:3*NATOMS)
  ENDDO
! bs360: what is it good for? EBEST and BESTCOORDS never get updated. (July 2010)
  WRITE(MYUNIT2,'(A)') 'new restart procedure - JBEST, EBEST, BESTCOORDS'
  WRITE(MYUNIT2, '(A,I8)') 'run number ',JP
  WRITE(MYUNIT2, '(I8,F25.15)') JBEST(JP), EBEST(JP)
  WRITE(MYUNIT2, '(3F25.15)') BESTCOORDS(1:3*NATOMS,JP)
! bs360: the following is not tested for MPI (July 2010)
  WRITE(MYUNIT2, '(A)') 'new restart procedure - saved energies and coordinates in MSB list'
  WRITE(MYUNIT2, '(I8)') NMSBSAVE
  DO J1=1,MIN(NMSBSAVE,MAXSAVE)
     WRITE(MYUNIT2, '(A,I8)') 'structure number ',J1
     WRITE(MYUNIT2, '(G25.15)') MSBE(J1)
     WRITE(MYUNIT2, '(3F25.15)') MSBCOORDS(1:3*NATOMS,J1)
  ENDDO
! WRITE(MYUNIT2, '(A)','old taboo procedure: energy and inertia determinant'
! DO J1=1,NPAR
   ! WRITE(MYUNIT2, '(A,I8)') 'taboo structure ',J1
   ! WRITE(MYUNIT2, '(2G25.15)') ESAVE(J1), XINSAVE(J1)
! ENDDO
!
  CLOSE(MYUNIT2)


!   csw34> Add in dumping of the lowest interaction energies if A9INTE is specified
IF (NPAR.EQ.1) THEN
  IF (A9INTET) THEN
     CALL SYSTEM('cp GMINinte.dump GMINinte.dump.save')

     OPEN(UNIT=1,FILE='GMINinte.dump',STATUS='UNKNOWN')
     WRITE(1, '(A)') 'SAVEINTE'
     WRITE(1, '(I8)') NSAVEINTE
     WRITE(1, '(A)') 'steps completed J1 in mc'
     WRITE(1, '(I8)') NDONE
     WRITE(1, '(A)') 'INTEQMIN and INTEQMINP'
     DO J1=1,NSAVEINTE
        WRITE(1, '(A,I8)') 'saved interaction enthalpy minimum ',J1
        WRITE(1, '(G25.15)') INTEQMIN(J1)
        WRITE(1, '(3F25.15)') INTEQMINP(J1,1:3*NATOMS)
     ENDDO
     CLOSE(UNIT=1)
  ENDIF
ENDIF

RETURN
END SUBROUTINE DUMPSTATE

SUBROUTINE RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
USE COMMONS,ONLY : NATOMS, COORDS, NPAR, STEP, ASTEP, OSTEP, TEMP, NMSBSAVE, MSBE, MSBCOORDS, DUMPFILE, &
  &                NSAVE, MAXSAVE, MPIT, AMHT, SEEDT, A9INTET, NSAVEINTE, INTEDUMPFILE, INTERESTORE
!USE QMODULE
USE V 

IMPLICIT NONE
INTEGER NDONE, JBEST(NPAR), JP, J1, OLDSAVE, MYUNIT2
DOUBLE PRECISION EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR), DUMMYARRAY(3*NATOMS)
CHARACTER(LEN=1) DUMMYS
CHARACTER(LEN=20) :: ISTR

IF (NPAR.GT.1) THEN
   WRITE (ISTR, '(i10)') JP
   MYUNIT2=JP+30000
   OPEN(MYUNIT2,FILE="GMIN.dump."//trim(adjustl(istr)),STATUS='OLD')
ELSE
   MYUNIT2=1
   OPEN(MYUNIT2,FILE='GMIN.dump',STATUS='OLD')
ENDIF

   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,'(I8)') NDONE
   READ(MYUNIT2,'(I8)') OLDSAVE
   READ(MYUNIT2,*) DUMMYS
! Read in last minimum, which is also the last minimum in the Markov chain.
   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,*) COORDS(1:3*NATOMS,JP)
   IF ((.NOT.SEEDT).AND.(.NOT.AMHT)) THEN
      WRITE(LFH,'(A,I4)') 'Initial coordinates: process',JP
      WRITE(LFH,'(3F20.10)') (COORDS(J1,JP),J1=1,3*NATOMS)
   ENDIF
! step and temperature information
   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,*) STEP(JP),ASTEP(JP),OSTEP(JP),TEMP(JP)
! best OLDSAVE minima
  READ(MYUNIT2,*) DUMMYS
!       csw34> this is the loop that has been edited to allow SAVE to vary between
!       runs. There are three cases which must be considered. 1) NSAVE=OLDSAVE.
!       In this case, nothing new needs to be done. 2) NSAVE<OLDSAVE. Now, the
!       extra minima that are saved in the dump file need to be cycled over
!       (ignored). 3) NSAVE>OLDSAVE. The array must be allocated with
!       dimensions suitable for the new number of minima, and the dump file read
!       in as far as possible. The rest of the array must then be padded with
!       zeros. The QMIN and QMINP arrays are already allocated in main.F with
!       the correct dimension so that need not worry us here :)
        IF (NSAVE.EQ.OLDSAVE) THEN
           DO J1=1,NSAVE
              READ(MYUNIT2,*) DUMMYS
              READ(MYUNIT2,*) QMIN(J1)
              READ(MYUNIT2,*) QMINP(J1,1:3*NATOMS)
           ENDDO
        ELSEIF (NSAVE.LT.OLDSAVE) THEN
           PRINT *,'NSAVE<OLDSAVE - truncating read in'
           DO J1=1,OLDSAVE
              IF (J1.LE.NSAVE) THEN
                 READ(MYUNIT2,*) DUMMYS
                 READ(MYUNIT2,*) QMIN(J1)
                 READ(MYUNIT2,*) QMINP(J1,1:3*NATOMS)
              ELSE 
                 READ(MYUNIT2,*) DUMMYS
                 READ(MYUNIT2,*) DUMMYS
                 READ(MYUNIT2,*) DUMMYARRAY(1:3*NATOMS)
              ENDIF
           ENDDO
        ELSEIF (NSAVE.GT.OLDSAVE) THEN
           PRINT *,'NSAVE>OLDSAVE - padding QMIN and QMINP with zeros'
           DO J1=1,NSAVE
              IF (J1.LE.OLDSAVE) THEN
                 READ(MYUNIT2,*) DUMMYS
                 READ(MYUNIT2,*) QMIN(J1)
                 READ(MYUNIT2,*) QMINP(J1,1:3*NATOMS)
              ELSE 
                 QMIN(J1)=0.0D0
                 QMINP(J1,1:3*NATOMS)=0.0D0
              ENDIF
           ENDDO
        ENDIF
! read in EBEST and BESTCOORDS
! bs360: what is it good for? EBEST and BESTCOORDS never get updated. (July 2010)
   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,*) DUMMYS
   READ(MYUNIT2,*) JBEST(JP), EBEST(JP)
   READ(MYUNIT2,*) BESTCOORDS(1:3*NATOMS,JP)
! bs360: the following is not tested for MPI (July 2010)
  READ(MYUNIT2,*) DUMMYS
  READ(MYUNIT2,*) NMSBSAVE
  DO J1=1,MIN(NMSBSAVE,MAXSAVE)
     READ(MYUNIT2,*) DUMMYS
     READ(MYUNIT2,*) MSBE(J1)
     READ(MYUNIT2,*) MSBCOORDS(1:3*NATOMS,J1)
  ENDDO
! READ(MYUNIT2,*) DUMMYS
! DO J1=1,NPAR
   ! READ(MYUNIT2,*) DUMMYS
   ! READ(MYUNIT2,*) ESAVE(J1), XINSAVE(J1)
! ENDDO

  CLOSE(MYUNIT2)


!   csw34> Add in restoring of the lowest interaction energies if A9INTE is specified
!          INTERESTORE is .TRUE. if the RESTORE keyword has two arguements. The second is
!          the name of the interaction energy dump file!
IF (NPAR.EQ.1) THEN
IF (A9INTET.AND.INTERESTORE) THEN
   OPEN(UNIT=1,FILE=TRIM(ADJUSTL(INTEDUMPFILE)),STATUS='OLD')
   READ(1,*) DUMMYS
   READ(1,'(I8)') OLDSAVE
   READ(1,*) DUMMYS
   READ(1,'(I8)') NDONE
   READ(1,*) DUMMYS
!       csw34> loop as above to allow SAVE to vary between runs
   IF (NSAVEINTE.EQ.OLDSAVE) THEN
           DO J1=1,NSAVEINTE
                   READ(1,*) DUMMYS
                   READ(1,*) INTEQMIN(J1)
                   READ(1,*) INTEQMINP(J1,1:3*NATOMS)
           ENDDO
!       csw34> NSAVEINTE<OLDSAVE
   ELSE IF (NSAVEINTE.LT.OLDSAVE) THEN
           PRINT *,'NSAVEINTE<OLDSAVE - truncating read in'
           DO J1=1,OLDSAVE
                   IF (J1.LE.NSAVEINTE) THEN
                           READ(1,*) DUMMYS
                           READ(1,*) INTEQMIN(J1)
                           READ(1,*) INTEQMINP(J1,1:3*NATOMS)
                   ELSE
                           READ(1,*) DUMMYS
                           READ(1,*) DUMMYS
!       csw34> This is the key bit, you've got to skip the correct number of
!       lines and that includes those for the coordinates. The new array
!       DUMMYARRAY allows this in a very simple way.
                           READ(1,*) DUMMYARRAY(1:3*NATOMS)
                  ENDIF
           ENDDO
!       csw34> NSAVEINTE>OLDSAVE
   ELSE IF (NSAVEINTE.GT.OLDSAVE) THEN
           PRINT *,'NSAVEINTE>OLDSAVE - padding INTEQMIN and INTEQMINP with zeros'
           DO J1=1,NSAVEINTE
                   IF (J1.LE.OLDSAVE) THEN
                           READ(1,*) DUMMYS
                           READ(1,*) INTEQMIN(J1)
                           READ(1,*) INTEQMINP(J1,1:3*NATOMS)
                   ELSE
                           INTEQMIN(J1)=0.0D0
                           INTEQMINP(J1,1:3*NATOMS)=0.0D0
                   ENDIF
           ENDDO
   ENDIF
CLOSE(UNIT=1)
ENDIF
ENDIF

RETURN
END SUBROUTINE RESTORESTATE
