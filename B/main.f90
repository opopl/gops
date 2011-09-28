
      PROGRAM GMIN
!op226> Declarations {{{ 
      !op226> Modules {{{ 
      USE COMMONS
      USE V
      USE F 
      USE KW
      USE DV
      USE PORFUNCS
      !op226>  }}}
      IMPLICIT NONE

      integer j1

! }}}
! intro {{{
      CALL CPU_TIME(TSTART)
      CALL RCA
      CALL INITVARS("FILES")
      CALL INITVARS("PREFFILES")

      ! write initial output {{{
  
      OPEN(LFH,FILE=LO_FILE, STATUS="unknown", form="formatted")
      OPEN(EA_FH,FILE=EA_FILE, STATUS="unknown", form="formatted")

      WRITE(LFH, '(A,I10,A,I10,A)') "Starting serial execution"

      CALL ED(LFH) 
      CALL PRINTTIME(LFH)
      CALL DISPLAY_VERSION(LFH)
      IF (USERCA) then 
        WRITE(LFH,'(A)') 'Command-line:  '
        WRITE(LFH,'(A)') ADJUSTL(CMDLINE)
        CALL ED(LFH) 
      ENDIF

      WRITE(LFH,'(A)') "Input files:"
      WRITE(LFH,'(A)') ''
      WRITE(LFH,'(2A50)') "Input data file:      ",D_FILE
      WRITE(LFH,'(2A50)') "Input coords file:    ",C_FILE
      WRITE(LFH,'(A)') ''
      WRITE(LFH,'(A)') "Output files:"
      WRITE(LFH,'(A)') ''
      WRITE(LFH,'(2A50)') "Runtime info file:  ",LO_FILE
      WRITE(LFH,'(2A50)') "Lowest energy geometries:  ",LE_FILE
      WRITE(LFH,'(2A50)') "Lowest energies:  ",EA_FILE
      CALL ED(LFH) 
      ! }}}

      ! }}}
!op226> Allocate memory; open files; initialize different things  {{{ 
      CALL AM("INIT")
      CALL INITVARS("LOG")
      CALL INITVARS("VARS")
      CALL COUNTATOMS
      CALL INITVARS("ARR")
      CALL KEYWORD(1)
      CALL RCA
      CALL AM("MAIN")
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

      QMINP(1:NSAVE,1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
      COORDSO(1:3*NATOMS,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      FF(1:NSAVE)=0 ! to prevent reading from uninitialised memorY
      VATO(1:NATOMS,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory

!!op226> DUMPT {{{ 
      !IF (DUMPT) THEN
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
!!op226>}}} 
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
      IF (TRACKENERGY) THEN
            OPEN(ENERGY_FH,FILE='energy',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
      ENDIF
      IF (TRACKMARKOV) THEN
            OPEN(MARKOV_FH,FILE='markov',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
      ENDIF
      IF (TRACKBEST) THEN
            OPEN(BEST_FH,FILE='best',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
      ENDIF
!op226>}}} 
      CALL FLUSH(6)
      CALL IOM
      CALL FLUSH(6)

!op226> SEEDT {{{ 

         !CALL GSEED
         !IF ((.NOT.FIELDT).AND.CENT) THEN
               !IF (.NOT.SEEDT) CALL CENTRE2(COORDS(1:3*NATOMS,1))
         !ELSEIF ((.NOT.FIELDT).AND.FIXCOM) THEN
            !DO J1=1,NPAR
               !IF (.NOT.SEEDT) CALL CENTRECOM(COORDS(1:3*NATOMS,J1))
            !ENDDO
         !ENDIF
      !ENDIF
!op226>}}} 
      
         NQ=1
         QMIN=1.0D10
!op226> End initializations and allocations }}} 

      CALL MC(MCSTEPS,TFAC,MSCREENC)

      CALL FINALQ
      CALL FINALIO
      CALL DEAM

!op226> Deallocate memory; close files {{{ 
      CALL FLUSH(LFH)

      CLOSE(LFH)
      CLOSE(EA_FH)
! csw34> close pairdists.* files
      IF (PAIRDISTT) CLOSE(PAIRDIST_FH)
!op226> closing files: energy best markov {{{ 
      IF (TRACKENERGY) THEN
         CLOSE(ENERGY_FH)
      ENDIF
      IF (TRACKMARKOV) THEN
         CLOSE(MARKOV_FH)
      ENDIF
      IF (TRACKBEST) THEN
         CLOSE(BEST_FH)
      ENDIF
!op226>}}} 

!op226>}}} 
      STOP

      END
