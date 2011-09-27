SUBROUTINE RCA
! declarations {{{
USE PORFUNCS
USE DV
USE V
USE F
USE COMMONS

IMPLICIT NONE

INTEGER NARGS, I
CHARACTER(LEN=80) BFF,VAR,S
! }}}
! subroutine body {{{

CALL IARGC_SUBR(NARGS)
I=0
CALL GETARG_SUBR(I,PROGNAME)

CMDLINE=trim(PROGNAME)

IF (NARGS.GT.0) THEN
  USERCA=.TRUE.
  DO WHILE(I.LT.NARGS) ; I=I+1 ; CALL GETARG_SUBR(I,BFF)
     CMDLINE=trim(CMDLINE)//' '//BFF
     SELECTCASE(BFF)
       CASE('-v')  ; CALL DISPLAY_VERSION(0) ; STOP
       CASE('-pv') ; CALL INITVARS("ALL") ; CALL PRINTVARS(0) ; STOP
       CASE('-h')  ; CALL PRINTHELP(0) ; STOP
       CASE DEFAULT
         CALL GETARG_SUBR(I+1,VAR)
         SELECTCASE(BFF)
                CASE('-f')
                        PULLT=.TRUE.
                        READ(VAR,*) PFORCE
                        !WRITE(*,*) 'FORCE:',PFORCE
				CASE('-ca') ; READ(VAR,*) NACCEPT 
				CASE('-ediff') ; READ(VAR,*) ECONV
				CASE('-g') ;  DEBUG=.TRUE.
				CASE('-g46') ;  BLNT=.TRUE.; G46=.TRUE.; P46=.FALSE.
				CASE('-p46') ;  BLNT=.TRUE.; P46=.TRUE.; G46=.FALSE.
				CASE('-mbln') ; MYBLNT=.TRUE.
				CASE('-nrg') ; READ(VAR,*) NRG
				CASE('-steps') ; READ(VAR,*) MCSTEPS
				CASE('-p') ; READ(VAR,*) S; USEPREF=.TRUE. ; PREF=trim(S)//"."
				CASE('-crd') ; READ(VAR,*) S; C_FILE=TRIM(S)
				CASE('-cmarkov') ; CHECKMARKOVT=.TRUE.
                ! tracking
				CASE('-track') ; TRACKDATAT=.TRUE.
				CASE('-tenergy') ; TRACKENERGY=.TRUE.
				CASE('-tbest') ; TRACKBEST=.TRUE.
				CASE('-tmarkov') ; TRACKMARKOV=.TRUE.
                CASE DEFAULT
         ENDSELECT
     ENDSELECT
  ENDDO
ENDIF
! }}}
END SUBROUTINE RCA

