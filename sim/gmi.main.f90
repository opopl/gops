      
!> @name GMIN
!
!> @brief Main program file for GMIN
!

      PROGRAM GMI
! declarations {{{
      USE V
      USE KW
      USE PORFUNCS
      USE FUNC
      USE MC
      
      IMPLICIT NONE
      
! }}}

! init {{{
      ! initialize variables
      CALL INITVARS

      CALL CPU_TIME(TSTART)

      CALL READ_CMD_ARGS
      ! read the input 'data' file
      CALL KEYWORD(1)
      ! set vars depending on the input provided
      ! through KEYWORS and READ_CMD_ARGS
      CALL SETVARS

      ! open the output log file for writing
      CALL OPENF(LFH,">",LFN)

      WRITE(LFH,'(A)') "Starting serial execution" 

      ALLOCATE(SCREENC(3*NATOMS))
      
      CALL IO
      CALL CENTRE(SCREENC)

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
