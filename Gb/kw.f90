
MODULE KW
! {{{
USE COMMONS
USE FUNC
USE STRINGS

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
        CASE(1) 
          ! "data" file {{{
      CALL INQF(D_FILE,YESNO)

      IF (YESNO) THEN 
         USEKW=.true.
         CALL OPENF(DATA_FH,"<",D_FILE)
      ELSE
         USEKW=.FALSE.
         PRINT '(A)','No input ',D_FILE,' file found. Using build-in values + command-line parameters...'
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
    	      CASE('  ','NOTE','COMMENT','\\','#','!')
    	      CASE('STOP')
    	         RETURN
                 ! A-L {{{
    	      CASE('ACCRAT') ; READ(BUFFER, *) ACCRAT
    	      CASE('CENTRE') ; CENT=.TRUE.
              CASE('DEBUG') ; DEBUG=.TRUE.
              CASE('DGUESS') ; READ(BUFFER, *) DGUESS
              CASE('EDIFF') ;  READ(BUFFER, *) EDIFF
    	      CASE('FQMAX','TIGHTCONV')
    	             READ(BUFFER, *) FQMAX
              CASE('G46')
                 G46=.TRUE.
                 BLNT=.TRUE.
                 WRITE(*,*) "P46 BLN model - Go-like model"
                 ! }}}
                 ! M-S {{{
              CASE('MAXIT')
    	         READ(BUFFER, *) MAXIT
              CASE('UPDATES','MUPDATE') ; READ(BUFFER, *) MUPDATE
              CASE('MAXBFGS');  READ(BUFFER, *) MAXBFGS
              CASE('NACCEPT','CHANGEACCEPT') ; READ(BUFFER, *) NACCEPT
              CASE('NORESET') ; NORESET=.TRUE.
              CASE('NSAVE','SAVE') ; READ(BUFFER, *) NSAVE
              CASE('P46')
                 P46=.TRUE.
                 BLNT=.TRUE.
                 WRITE(*,*) "P46 BLN model - Wild-type"
              CASE('PULL')
                 PULLT=.TRUE.
    	         READ(BUFFER, *) PATOM1,PATOM2,PFORCE
              CASE('RADIUS');  READ(BUFFER, *) RADIUS
              CASE('STEP');  READ(BUFFER, *) STEP, ASTEP
    	      CASE('SQMAX','BASIN','SLOPPYCONV') ; READ(BUFFER, *) SQMAX
              CASE('STEPS') ; READ(BUFFER, *) MCSTEPS,TFAC
                 ! }}} 
              CASE('TARGET')
                TARGET=.TRUE.
                NTARGETS=NARGS-1
                ALLOCATE(TARGETS(NTARGETS))
    	        READ(BUFFER, *) TARGETS(1:NTARGETS)
              CASE('TEMPERATURE'); READ(BUFFER, *) TEMP
              CASE('TRACKDATA') ; TRACKDATAT=.TRUE.
              CASE DEFAULT
                 write(*,'(A10,A1,A40,A20)') adjustr(PNAME),">","Unrecognized Keyword: ",LABEL
    	         STOP
                 ! }}}
          ENDSELECT
       END IF
       ! }}}
      END DO

      ENDSELECT

      IF (USEKW) THEN
        CLOSE(DATA_FH)
      ENDIF

      RETURN
      ! }}}
      END SUBROUTINE
! }}}
ENDMODULE
