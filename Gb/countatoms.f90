
MODULE NOA

      IMPLICIT NONE
      SAVE
      INTEGER :: NUMBER_OF_ATOMS

      CONTAINS

      SUBROUTINE COUNTATOMS(MYUNIT)

      INTEGER EOF,MYUNIT
      LOGICAL :: YESNO

      YESNO=.FALSE.
      
      INQUIRE(FILE='coords',EXIST=YESNO)

      NUMBER_OF_ATOMS=0

      IF (YESNO) THEN
         OPEN(UNIT=7,FILE='coords',STATUS='OLD')
         DO
            READ(7,*,IOSTAT=EOF)
            IF (EOF==0) THEN
               NUMBER_OF_ATOMS = NUMBER_OF_ATOMS + 1
            ELSE
               EXIT
            ENDIF
         ENDDO
        ELSEIF (YESNOAMH) THEN
      ELSE
         PRINT '(A)','ERROR - no coords, input.crd, coords.inpcrd or coords.amber file'
         STOP
      ENDIF

      CLOSE(7)
      
      END SUBROUTINE COUNTATOMS

END MODULE NOA
