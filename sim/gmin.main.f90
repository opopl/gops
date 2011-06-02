      
!> @name GMIN
!
!> @brief Main program file for GMIN
!

      PROGRAM GMIN
! declarations {{{
      USE COMMONS
      USE PORFUNCS
      USE FUNC
      
      IMPLICIT NONE
      
      DOUBLE PRECISION, ALLOCATABLE :: SCREENC(:,:)
      DOUBLE PRECISION E

      ! USED IN WRITING TO PAIRDIST      
      CHARACTER(LEN=40) :: ATOM1,ATOM2,ATOMPAIR

      COMMON /MYPOT/ E
! }}}

! init {{{
      CALL CPU_TIME(TSTART)
      CALL COUNTATOMS
      CALL KEYWORD(1)
      CALL READ_CMD_ARGS

      CALL OPENF(LFH,">",LFN)

      WRITE(LFH,'(A)') "Starting serial execution" 

      ALLOCATE(SCREENC(NATOMS,3))
      
      include gmin.main.pairdist.inc.f90
      
      IF (TRACKDATAT) THEN
         CALL OPENF(ENERGY_FH,">>","energy.dat")
         CALL OPENF(MARKOV_FH,">>","markov.dat")
         CALL OPENF(BEST_FH,">>","best.dat")
         IF (RMST) CALL OPENF(RMSD_FH,">>","rmsd.dat")
      ENDIF
      
      CALL IO
      CALL CENTRE2(COORDS)

      NQ=1

! }}}
      CALL MCRUNS(SCREENC)
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
