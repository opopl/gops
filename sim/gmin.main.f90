      
!> @name GMIN
!
!> @brief Main program file for GMIN
!

      PROGRAM GMIN
! declarations {{{
      USE V
      USE KW
      USE PORFUNCS
      USE FUNC
      USE MC
      
      IMPLICIT NONE
      
      DOUBLE PRECISION, ALLOCATABLE :: SCREENC(:,:)
! }}}

! init {{{
      ! initialize variables
      CALL INITVARS

      CALL CPU_TIME(TSTART)

      CALL READ_CMD_ARGS
      CALL COUNTATOMS
      CALL KEYWORD(1)

      CALL OPENF(LFH,">",LFN)

      WRITE(LFH,'(A)') "Starting serial execution" 

      ALLOCATE(SCREENC(NATOMS,3))
      
      !include "gmin.main.pairdist.inc.f90"
      
      IF (TRACKDATAT) THEN
         CALL OPENF(ENERGY_FH,">>","energy.dat")
         CALL OPENF(MARKOV_FH,">>","markov.dat")
         CALL OPENF(BEST_FH,">>","best.dat")
         !IF (RMST) 
        CALL OPENF(RMSD_FH,">>","rmsd.dat")
      ENDIF
      
      CALL IO
      CALL CENTRE2(COORDS)

      NQ=1

! }}}

      CALL MCRUN(MCSTEPS,TFAC,SCREENC)

! end {{{
      CALL FLUSH(LFH)
      CLOSE(LFH)

      IF (PAIRDISTT) CLOSE(PAIRDIST_FH)
      IF (TRACKDATAT) THEN
         CLOSE(ENERGY_FH)
         CLOSE(MARKOV_FH)
         IF (RMST) CLOSE(RMSD_FH)
         CLOSE(BEST_FH)
      ENDIF
! }}}

      STOP

      END
