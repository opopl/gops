
      PROGRAM GMIN

      USE KW
      USE COMMONS
      USE MCFUNC
      USE FUNC
      USE QMODULE
      USE PORFUNCS

      IMPLICIT NONE

      ! variables {{{

      INTEGER J1,J2, JP, MPIERR, NDUMMY3,NPTOTAL
      DOUBLE PRECISION, ALLOCATABLE :: SCREENC(:)
      DOUBLE PRECISION POTEL
      DOUBLE PRECISION, ALLOCATABLE :: TMPCOORDS(:)
      INTEGER, ALLOCATABLE :: NDUMMY(:), NDUMMY2(:,:)

      LOGICAL LOPEN

      CHARACTER(LEN=130) ISTR
      CHARACTER(LEN=40) :: ATOM1,ATOM2,ATOMPAIR

      COMMON /MYPOT/ POTEL

      ! }}}

      CALL CPU_TIME(TSTART)

      CALL INITVARS
      CALL RCA

      MYNODE=0
      OPEN(LFH,FILE=LF_FILE, STATUS="unknown", form="formatted")
      WRITE(LFH, '(A,I10,A,I10,A)') "Starting serial execution"

!op226> Allocate memory; open files; initialize different things  {{{ 

      CALL COUNTATOMS
      CALL AM
 
      ALLOCATE(SCREENC(3*NATOMS))
      INQUIRE(UNIT=1,OPENED=LOPEN)

      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'main> A ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      CALL KEYWORD(1)

      INQUIRE(UNIT=1,OPENED=LOPEN)

      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'main> B ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      ALLOCATE(FF(NSAVE),QMIN(NSAVE))
      ALLOCATE(QMINP(NSAVE,3*NATOMS))

      QMINP(1:NSAVE,1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
      COORDSO(1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
      FF(1:NSAVE)=0 ! to prevent reading from uninitialised memorY
      VATO(1:NATOMS)=0.0D0 ! to prevent reading from uninitialised memory

!op226> DUMPT {{{ 
!      IF (DUMPT) THEN
         !DUMPXYZUNIT=40
         !DUMPVUNIT=39
         !DO J1=1,NPAR
            !WRITE (ISTR,'(A5,I1,A4)') 'dump.',J1,'.xyz'
            !J2=DUMPXYZUNIT+J1
            !OPEN(UNIT=J2,FILE=TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
            !WRITE (ISTR,'(A5,I1,A2)') 'dump.',J1,'.V'
            !J2=DUMPVUNIT-J1
            !OPEN(UNIT=J2,FILE=TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
         !ENDDO
      !ENDIF
!op226>}}} 

      ! PAIRDISTT {{{
      IF (PAIRDISTT) THEN
         MYPUNIT=3000+MYNODE
         OPEN(MYPUNIT,FILE="pairdists",STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
         WRITE(MYPUNIT,'(A10)',ADVANCE="NO") "Quench  "
         DO J1=1,NPAIRS
            WRITE(ATOM1,*) PAIRDIST(J1,1)
            WRITE(ATOM2,*) PAIRDIST(J1,2)
            WRITE(ATOMPAIR,*) TRIM(ADJUSTL(ATOM1))//"-"//TRIM(ADJUSTL(ATOM2))
            WRITE(MYPUNIT,'(A10)',ADVANCE="NO") TRIM(ADJUSTL(ATOMPAIR))//"  " 
         ENDDO
         WRITE(LFH,'(A)') ""
      ENDIF
      ! }}}

!op226> TRACKDATAT {{{ 
!
!     csw34> TRACKDATA keyword prints the energy and markov energy 
!     to files for viewing during a run. If RMS is also specified it
!     prints the rmsd from the comparison structure into a file.
!
      IF (TRACKDATAT) THEN
         MYEUNIT=4000+MYNODE
         MYMUNIT=6000+MYNODE
         MYRUNIT=8000+MYNODE
         MYBUNIT=10000+MYNODE
            OPEN(MYEUNIT,FILE='energy',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYMUNIT,FILE='markov',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYBUNIT,FILE='best',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
      ENDIF
!op226>}}} 

      CALL FLUSH(6)
      CALL IO1
      CALL FLUSH(6)

!op226> SEEDT {{{ 

!      IF (SEEDT) THEN
         !CALL GSEED
      !ELSE
         !IF ((.NOT.FIELDT).AND.CENT) THEN
            !DO J1=1,NPAR
               !IF (.NOT.SEEDT) CALL CENTRE2(COORDS(1:3*NATOMS,J1))
            !ENDDO
         !ELSEIF ((.NOT.FIELDT).AND.FIXCOM) THEN
            !DO J1=1,NPAR
               !IF (.NOT.SEEDT) CALL CENTRECOM(COORDS(1:3*NATOMS,J1))
            !ENDDO
         !ENDIF
      !ENDIF
!op226>}}} 

         NQ=1
      DO J1=1,NSAVE
         QMIN(J1)=1.0D10
      ENDDO
!op226> End initializations and allocations }}} 

      CALL MCRUNS(SCREENC)

!op226> Deallocate memory; close files {{{ 
      !IF (ALLOCATED(XICOM)) DEALLOCATE(XICOM)
      !IF (ALLOCATED(PCOM)) DEALLOCATE(PCOM)
      CALL FLUSH(LFH)

      CLOSE(LFH)
! csw34> close pairdists.* files
      IF (PAIRDISTT) CLOSE(MYPUNIT)
!op226> closing files for: TRACKDATAT RMST A9INTET {{{ 
      IF (TRACKDATAT) THEN
         CLOSE(MYEUNIT)
         CLOSE(MYMUNIT)
         CLOSE(MYBUNIT)
      ENDIF
!op226>}}} 

!op226>}}} 
      STOP

      END
