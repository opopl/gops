
      PROGRAM GMIN
!op226> Declarations {{{ 
      !op226> Modules {{{ 
      USE COMMONS
      USE V
      USE F 
      USE PORFUNCS
      !op226>  }}}
      !op226> Other {{{
      IMPLICIT NONE
      !EXTERNAL READ_CMD_ARGS

      INTEGER J1,J2, JP, MPIERR, NDUMMY3,NPTOTAL,VERSIONTEMP
      DOUBLE PRECISION, ALLOCATABLE :: SCREENC(:)
      DOUBLE PRECISION POTEL
      DOUBLE PRECISION, ALLOCATABLE :: TMPCOORDS(:)
      INTEGER, ALLOCATABLE :: NDUMMY(:), NDUMMY2(:,:)
      LOGICAL LOPEN

      CHARACTER(LEN=130) ISTR
      CHARACTER(LEN=40) :: atom1,atom2,atompair

      COMMON /MYPOT/ POTEL
!op226>  }}}
!op226> End declarations }}}
! intro {{{
      CALL CPU_TIME(TSTART)
      CALL RCA

      LFH=10
      OPEN(LFH,FILE=O_FILE, STATUS="unknown", form="formatted")
      WRITE(LFH, '(A,I10,A,I10,A)') "Starting serial execution"

      write(*,'(2A20)') "Input data file:      ",D_FILE
      write(*,'(2A20)') "Input coords file:    ",C_FILE
      write(*,'(2A20)') "Output runtime file:  ",O_FILE
      NPAR=1
      MYNODE=0
      ! }}}
!op226> Allocate memory; open files; initialize different things  {{{ 
      CALL COUNTATOMS
      ALLOCATE(SCREENC(3*NATOMS),VT(NATOMS))

      CALL KEYWORD
      !CALL KW
      CALL SETVARS

!!op226> RMST {{{ 
      !IF (RMST) THEN
         !ALLOCATE(RMSBEST(RMSSAVE,2),RMSCOOR(RMSSAVE,3*NATOMS))
         !RMSBEST(1:RMSSAVE,1)=RMSLIMIT+RMSTOL
         !RMSBEST(1:RMSSAVE,2)=0.D0
         !RMSCOOR(1:RMSSAVE,1:3*NATOMS)=0.D0
         !ALLOCATE(COORCOMP(1:3*NATOMS))
!!
!!        csw34> Need to add the read in for other file formats
!!        (crd,pdb). Still not sure how GMIN reads in the input.crd file
!!        for CHARMM. It appears to be reading into unit 7 in io1.f
!!        though - but only from xyz format! Is there a hidden
!!        conversion?
!!         
         !OPEN(UNIT=1,FILE='compare',STATUS='OLD')
         !READ(1,*) (COORCOMP(J1),J1=1,3*NATOMS)
         !CLOSE(1)
      !ENDIF
!!op226>}}} 

      ALLOCATE(FF(NSAVE),QMIN(NSAVE))
      ALLOCATE(QMINP(NSAVE,3*NATOMS))

!        csw34> ALLOCATE the interaction energy tracking arrays if A9INTE in data

      QMINP(1:NSAVE,1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
      COORDSO(1:3*NATOMS,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      FF(1:NSAVE)=0 ! to prevent reading from uninitialised memorY
      VATO(1:NATOMS,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory

!op226> DUMPT {{{ 
      IF (DUMPT) THEN
         DUMPXYZUNIT=40
         DUMPVUNIT=39
         DO J1=1,NPAR
            WRITE (ISTR,'(A5,I1,A4)') 'dump.',J1,'.xyz'
            J2=DUMPXYZUNIT+J1
            OPEN(UNIT=J2,FILE=TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
            WRITE (ISTR,'(A5,I1,A2)') 'dump.',J1,'.V'
            J2=DUMPVUNIT-J1
            OPEN(UNIT=J2,FILE=TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
         ENDDO
      ENDIF
!op226>}}} 

!        ! PAIRDIST {{{
      !IF (PAIRDISTT) THEN
         !MYPUNIT=3000+MYNODE
         !IF (NPAR.GT.1) THEN
            !OPEN(MYPUNIT,FILE="pairdists."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
         !ELSE
            !OPEN(MYPUNIT,FILE="pairdists",STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
         !ENDIF
         !WRITE(MYPUNIT,'(A10)',ADVANCE="NO") "Quench  "
         !DO J1=1,NPAIRS
            !WRITE(atom1,*) PAIRDIST(J1,1)
            !WRITE(atom2,*) PAIRDIST(J1,2)
            !WRITE(atompair,*) TRIM(ADJUSTL(atom1))//"-"//TRIM(ADJUSTL(atom2))
            !WRITE(MYPUNIT,'(A10)',ADVANCE="NO") TRIM(ADJUSTL(atompair))//"  " 
         !ENDDO
         !WRITE(LFH,'(A)') ""
      !ENDIF
         !! }}}

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
         IF (NPAR.GT.1) THEN
            OPEN(MYEUNIT,FILE="energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYMUNIT,FILE="markov."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYBUNIT,FILE="best."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
         ELSE
            OPEN(MYEUNIT,FILE='energy',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYMUNIT,FILE='markov',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            IF (RMST) OPEN(MYRUNIT,FILE='rmsd',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYBUNIT,FILE='best',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            IF (A9INTET) THEN
               OPEN(UNIT=3998,FILE='intE.dat',STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(UNIT=3999,FILE='bestintE.dat',STATUS='UNKNOWN',FORM='FORMATTED')
            ENDIF            
         ENDIF
      ENDIF
!op226>}}} 
      CALL FLUSH(6)
      CALL IOM
      CALL FLUSH(6)

!op226> SEEDT {{{ 

      IF (SEEDT) THEN
         CALL GSEED
      ELSE
         IF ((.NOT.FIELDT).AND.CENT) THEN
            DO J1=1,NPAR
               IF (.NOT.SEEDT) CALL CENTRE2(COORDS(1:3*NATOMS,J1))
            ENDDO
         ELSEIF ((.NOT.FIELDT).AND.FIXCOM) THEN
            DO J1=1,NPAR
               IF (.NOT.SEEDT) CALL CENTRECOM(COORDS(1:3*NATOMS,J1))
            ENDDO
         ENDIF
      ENDIF
!op226>}}} 
      
      DO JP=1,NPAR
         NQ(JP)=1
      ENDDO
      DO J1=1,NSAVE
         QMIN(J1)=1.0D10
      ENDDO
!op226> End initializations and allocations }}} 

!op226> Main program run
      !IF ((NRUNS.GT.0).OR.PTMC.OR.BSPT) CALL MCRUNS(SCREENC)

      CALL MC(MCSTEPS(1),TFAC(1),SCREENC)

      CALL FINALQ
      CALL FINALIO

!op226> Deallocate memory; close files {{{ 
      IF (ALLOCATED(FIN)) DEALLOCATE(FIN)
      CALL FLUSH(LFH)

      CLOSE(LFH)
! csw34> close pairdists.* files
      IF (PAIRDISTT) CLOSE(MYPUNIT)
!op226> closing files for: TRACKDATAT RMST A9INTET {{{ 
      IF (TRACKDATAT) THEN
         CLOSE(MYEUNIT)
         CLOSE(MYMUNIT)
         IF (RMST) CLOSE(MYRUNIT)
         CLOSE(MYBUNIT)
      ENDIF
!op226>}}} 

!op226>}}} 
      STOP

      END
