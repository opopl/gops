
      !include "kwde.i.f90"
      !include "kwi.i.f90"
      
      MODULE KW
! {{{
USE V
USE F
USE S

IMPLICIT NONE

! }}}
CONTAINS
      SUBROUTINE KEYWORD(PROG)
      ! declarations {{{
      INTEGER PROG
      LOGICAL YESNO

      ! INPUT RELATED VARIABLES

      CHARACTER(LEN=100) :: BUFFER, WORD
      CHARACTER(LEN=20),DIMENSION(20) :: ARGS
      INTEGER :: POS, LINE
      INTEGER, PARAMETER :: DATA_FH = 15
      INTEGER :: IOS=0
      INTEGER NARGS

      COMMON /BUFINF/ BUFFER,POS,WORD
      include "kwde.i.f90"
!,IOS
      ! }}}
! body {{{
      !IOS=0

      SELECTCASE(PROG)
        CASE(1) 
          ! "data" file {{{
      CALL INQF(D_FILE,YESNO)

      IF (YESNO) THEN 
         USEKW=.true.
         CALL OPENF(DATA_FH,"<",D_FILE)
      ELSE
         USEKW=.FALSE.
         PRINT '(A)','No data file found. Using build-in values + command-line parameters...'
      ENDIF

! ios<0 end of file;
! ios=0 
! ios>0 error 
      ! }}}
      DO WHILE (IOS == 0)
      ! {{{
100     READ(DATA_FH, '(A)', IOSTAT=IOS) BUFFER
        IF (IOS == 0) THEN

           POS = SCAN(BUFFER, '    ')
           WORD = BUFFER(1:POS)
           CALL PARSE(BUFFER,' ',ARGS,NARGS)
           BUFFER = BUFFER(POS+1:)
           WORD=ARGS(1)
           include "kwif.i.f90"
       END IF
       ! }}}
      END DO

      ENDSELECT

      IF (USEKW) THEN
        CLOSE(DATA_FH)
      ENDIF

      RETURN
      ! }}}
      END SUBROUTINE KEYWORD
! }}}

      ENDMODULE
