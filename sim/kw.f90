
MODULE KW
! {{{
USE V
uSE FUNC
uSE STRINGS

IMPLICIT NONE

CONTAINS

      SUBROUTINE KEYWORD(PROG)
      ! declarations {{{
      INTEGER PROG
      LOGICAL YESNO

      ! INPUT RELATED VARIABLES

      CHARACTER(LEN=100) :: BUFFER, LABEL
      INTEGER :: POS, LINE
      INTEGER, PARAMETER :: DATA_FH = 15
      INTEGER :: IOS=0

      COMMON /BUFINF/ BUFFER,POS,LABEL
!,IOS
      ! }}}
! body {{{
      !IOS=0

      SELECTCASE(PROG)
        CASE(1) ! GMI {{{
          ! "data" file {{{
      CALL INQF('data',YESNO)

      IF (YESNO) THEN 
         USEKW=.true.
         CALL OPENF(DATA_FH,"<","data")
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
        READ(DATA_FH, '(A)', IOSTAT=IOS) BUFFER
        IF (IOS == 0) THEN

           POS = SCAN(BUFFER, '    ')
           LABEL = BUFFER(1:POS)
           CALL PARSE(BUFFER,' ',ARGS,NARGS)
           BUFFER = BUFFER(POS+1:)
           LABEL=ARGS(1)
   
          SELECTCASE(LABEL)
              ! {{{
    	      CASE('  ','NOTE','COMMENT','\\','#')
    	      CASE('STOP')
    	         RETURN
                 ! A-L {{{
    	      CASE('ACCRAT')
                     READ(BUFFER, *) ACCRAT
    	      CASE('FQMAX','TIGHTCONV')
    	             READ(BUFFER, *) FQMAX
              CASE('EDIFF')
    	         READ(BUFFER, *) EDIFF
              CASE('G46')
                 G46=.TRUE.
                 BLNT=.TRUE.
                 WRITE(*,*) "P46 BLN model - Go-like model"
                 ! }}}
                 ! M-S {{{
              CASE('MAXIT')
    	         READ(BUFFER, *) MAXIT
              CASE('M_LBFGS')
    	         READ(BUFFER, *) M_LBFGS
              CASE('NACCEPT')
    	         READ(BUFFER, *) NACCEPT
              CASE('NORESET')
    	         NORESET=.TRUE.
              CASE('NSAVE')
    	         READ(BUFFER, *) NSAVE
              CASE('P46')
                 P46=.TRUE.
                 BLNT=.TRUE.
              CASE('PULL')
                 PULLT=.TRUE.
    	         READ(BUFFER, *) PATOM1,PATOM2,PFORCE
              CASE('RADIUS')
    	         READ(BUFFER, *) RADIUS
    	      CASE('SQMAX','BASIN','SLOPPYCONV')
    	         READ(BUFFER, *) SQMAX
              CASE('STEPS')
    	         READ(BUFFER, *) MCSTEPS,TFAC
                 ! }}} 
              CASE('TARGET')
                TARGET=.TRUE.
                NTARGETS=NARGS-1
                ALLOCATE(TARGETS(NTARGETS))
    	        READ(BUFFER, *) TARGETS(1:NTARGETS)
              CASE('TRACKDATA')
                TRACKDATAT=.TRUE.
              CASE DEFAULT
                 !CALL REPORT('Unrecognized command '//BUFFER,.TRUE.)
                 write(*,*) "Unrecognized command"
    	         STOP
                 ! }}}
          ENDSELECT
       END IF
       ! }}}
      END DO
      ! }}}
      ENDSELECT

      IF (USEKW) THEN
        CLOSE(DATA_FH)
      ENDIF

      RETURN
      ! }}}
      END SUBROUTINE
! }}}
ENDMODULE
