      
!> @name GMI
!
!> @brief Main program file for GMI
!

! conventions: SC

      PROGRAM GMI
! declarations {{{
      USE V
      USE KW
      USE PORFUNCS
      USE FUNC
      USE MC
      
      IMPLICIT NONE
      
! }}}
 
!  init {{{
      ! initialize variables
      CALL INITVARS

      CALL CPU_TIME(TSTART)

      ! read command-line arguments
      CALL RCA
      ! read the input 'data' file
      CALL KEYWORD(1)
      ! SETVARS():
      ! 
      ! set variables depending on the input provided
      ! after having called KEYWORD() and RCA()
      ! the variables set are: 
      !         NATOMS RADIUS
      CALL SETVARS

      ! open the output log file for writing
      CALL OPENF(LFH,">",LFN)

      WRITE(LFH,'(A)') "Starting serial execution" 
      WRITE(*,'(A)') "Starting serial execution" 

      ALLOCATE(SCREENC(NATOMS,3))

      CALL RCOORDS(NATOMS,RADIUS,SCREENC)
      WRITE(*,'(A)') "RCOORDS done" 
      
      CALL IO
      WRITE(*,'(A)') "IO done" 
      CALL CENTRE(SCREENC)
      WRITE(*,'(A)') "CENTRE done" 


! }}}

      CALL MCRUN
      WRITE(*,'(A)') "MCRUN done" 

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
