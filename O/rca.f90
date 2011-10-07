SUBROUTINE RCA
! declarations {{{
USE PORFUNCS
USE DV
USE COMMONS

IMPLICIT NONE

INTEGER NARGS, I
CHARACTER(LEN=80) BFF,VAR,S
! }}}
! subroutine body {{{

CALL IARGC_SUBR(NARGS)
I=0
<<<<<<< HEAD

IF (NARGS.GT.0) THEN
  DO WHILE(I.LT.NARGS) ; I=I+1 ; CALL GETARG_SUBR(I,BFF)
     SELECTCASE(BFF)
       CASE('-v')  ; CALL DISPLAY_VERSION(0) ; STOP
       !CASE('-pv') ; CALL INITVARS ; CALL PRINTVARS ; STOP
       !CASE('-h')  ; CALL PRINTHELP ; STOP
       CASE DEFAULT
!         CALL GETARG_SUBR(I+1,VAR)
!         SELECTCASE(BFF)
!                CASE('-f')
!                        PULLT=.TRUE.
!                        READ(VAR,*) PFORCE
!                        !WRITE(*,*) 'FORCE:',PFORCE
!				CASE('-ca') ; READ(VAR,*) NACCEPT 
!				CASE('-ediff') ; READ(VAR,*) ECONV
!				CASE('-g') ;  DEBUG=.TRUE.
!				CASE('-g46') ;  BLNT=.TRUE.; G46=.TRUE.
!				CASE('-p46') ;  BLNT=.TRUE.; P46=.TRUE.
!				CASE('-mbln') ; MYBLNT=.TRUE.
!				CASE('-nrg') ; READ(VAR,*) NRG
!                CASE DEFAULT
!         ENDSELECT
         !CALL MAKETEST(BFF)
       !CASE DEFAULT
         !INFILE=TRIM(ADJUSTL(BFF))
     ENDSELECT
         I=I+1
  ENDDO
ENDIF
! }}}
END SUBROUTINE RCA

