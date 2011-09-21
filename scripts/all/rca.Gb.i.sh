# {{{
cat << EOF
SUBROUTINE RCA
! declarations {{{
USE PORFUNCS
USE DV
USE COMMONS
USE FUNC

IMPLICIT NONE

INTEGER CMD_NARGS, I
CHARACTER(LEN=80) BFF,VAR
! }}}
! subroutine body {{{

CALL IARGC_SUBR(CMD_NARGS)
I=0

IF (CMD_NARGS.GT.0) THEN
  DO WHILE(I.LT.CMD_NARGS) ; I=I+1 ; CALL GETARG_SUBR(I,BFF)
     SELECTCASE(BFF)
       CASE('-v')  ; CALL DISPLAY_VERSION(0) ; STOP
       CASE('-pv') ; CALL INITVARS ; CALL PRINTVARS ; STOP
       !CASE('-h')  ; CALL PRINTHELP ; STOP
       CASE DEFAULT
         CALL GETARG_SUBR(i+1,var)
         SELECTCASE(BFF)
                CASE('-f')
                        READ(VAR,*) PFORCE
                        WRITE(*,*) 'FORCE:',PFORCE
				CASE('-g')
				  		DEBUG=.TRUE.
                CASE DEFAULT
         ENDSELECT
         !CALL MAKETEST(BFF)
       !CASE DEFAULT
         !INFILE=TRIM(ADJUSTL(BFF))
     ENDSELECT
         I=I+1
  ENDDO
ENDIF
! }}}
END SUBROUTINE RCA
EOF
# }}}
