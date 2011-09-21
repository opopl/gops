MODULE OUTPUT
      
   USE FILE_MANAGER, ONLY : CHECK_FILE
      
   PUBLIC :: WRITE_COORDS, WRITE_MARKOV_COORDS
      
   CONTAINS
      
      SUBROUTINE WRITE_COORDS(FILE_UNIT, FORMAT_SPEC, RUN_NUMBER)
         ! Does a sanity check and then writes COORDS to the specified unit with a 
         ! specified format.  Optional check for run number if we're doing parallel stuff.
         ! Commons
         USE COMMONS, ONLY : COORDS, NPAR
         IMPLICIT NONE
         ! Arguments
         INTEGER, INTENT(IN)           :: FILE_UNIT
         CHARACTER (LEN=*), INTENT(IN) :: FORMAT_SPEC
         ! Optional arguments
         INTEGER, INTENT(IN), OPTIONAL :: RUN_NUMBER
         ! Local variables
         INTEGER                       :: COUNTER
         
         ! Some quick sanity checks to make sure that input makes sense
         IF (PRESENT(RUN_NUMBER)) THEN
            IF (RUN_NUMBER .GT. NPAR) THEN
               STOP 'The run number is larger than the number of parallel runs.  Cannot write the output file.'
            END IF
         END IF
         
         ! Sanity checks on the file that we're writing to
         IF (.NOT. CHECK_FILE(FILE_UNIT, FOR_READ=.FALSE., FOR_WRITE=.TRUE.)) THEN
            STOP 'File did not open correctly.'
         END IF
         
         ! Write coords
          WRITE(FILE_UNIT, FORMAT_SPEC) COORDS(:)
         
      END SUBROUTINE WRITE_COORDS
      
      SUBROUTINE WRITE_MARKOV_COORDS(FILE_UNIT, FORMAT_SPEC, RUN_NUMBER)
         ! Does a sanity check and then writes COORDSO (the Markov coords) to the specified 
         ! unit with a specified format.  Optional check for run number if we're doing 
         ! parallel stuff.
         ! Commons
         USE COMMONS, ONLY : COORDSO, NPAR
         IMPLICIT NONE
         ! Arguments
         INTEGER, INTENT(IN)           :: FILE_UNIT
         CHARACTER (LEN=*), INTENT(IN) :: FORMAT_SPEC
         ! Optional arguments
         INTEGER, INTENT(IN), OPTIONAL :: RUN_NUMBER
         ! Local variables
         INTEGER                       :: COUNTER
      
         ! Some quick sanity checks to make sure that input makes sense
         IF (PRESENT(RUN_NUMBER)) THEN
            IF (RUN_NUMBER .GT. NPAR) THEN
               STOP 'The run number is larger than the number of parallel runs.  Cannot write the output file.'
            END IF
         END IF
         
         ! Sanity checks on the file that we're writing to
         IF (.NOT. CHECK_FILE(FILE_UNIT, FOR_READ=.FALSE., FOR_WRITE=.TRUE.)) THEN
            STOP 'File did not open correctly.'
         END IF
         
         ! Write coordso
         WRITE(FILE_UNIT, FORMAT_SPEC) COORDSO(:)
         
      END SUBROUTINE WRITE_MARKOV_COORDS
      
END MODULE OUTPUT
