#!/bin/bash

export shd="`dirname $(readlink -f $0)`"
export this_script=` basename $0 `

prog="$1"

rca_lbfgs(){
# {{{
cat << EOF
SUBROUTINE READ_CMD_ARGS

USE PORFUNCS
USE DV
!USE COMMONS,ONLY : INFILE
use commons, only: infile,pforce

IMPLICIT NONE

INTEGER NARGS, I
CHARACTER(LEN=80) BFF,VAR

CALL IARGC_SUBR(NARGS)
I=0

IF (NARGS.GT.0) THEN
  DO WHILE(I.LT.NARGS) ; I=I+1 ; CALL GETARG_SUBR(I,BFF)
     SELECTCASE(BFF)
       CASE('-v')
         CALL DISPLAY_VERSION(0)
         STOP
       case default
         CALL GETARG_SUBR(i+1,var)
         SELECTCASE(BFF)
                CASE('-f')
                        READ(VAR,*) PFORCE
                        WRITE(*,*) 'FORCE:',PFORCE
                CASE DEFAULT
         ENDSELECT
         !CALL MAKETEST(BFF)
       !CASE DEFAULT
         !INFILE=TRIM(ADJUSTL(BFF))
     ENDSELECT
         I=I+1
  ENDDO
ENDIF

END SUBROUTINE READ_CMD_ARGS
EOF
# }}}
}

rca_gmi(){
rca_gmin
}

rca_G(){
rca_GMIN
}

rca_GMIN(){
# {{{
cat << EOF
SUBROUTINE RCA
! declarations {{{
USE PORFUNCS
USE DV
USE COMMONS

IMPLICIT NONE

INTEGER NARGS, I
CHARACTER(LEN=80) BFF,VAR
! }}}
! subroutine body {{{

CALL IARGC_SUBR(NARGS)
I=0

IF (NARGS.GT.0) THEN
  DO WHILE(I.LT.NARGS) ; I=I+1 ; CALL GETARG_SUBR(I,BFF)
     SELECTCASE(BFF)
       CASE('-v')  ; CALL DISPLAY_VERSION(0) ; STOP
       !CASE('-pv') ; CALL INITVARS ; CALL PRINTVARS ; STOP
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
}

rca_gmin(){
# {{{
cat << EOF
SUBROUTINE RCA
! declarations {{{
USE PORFUNCS
USE FUNC, ONLY: PRINTVARS, INITVARS, PRINTHELP
USE DV
!USE COMMONS,ONLY : INFILE
!USE V, ONLY: INFILE,PFORCE
USE V, ONLY: PFORCE,DEBUG

IMPLICIT NONE

INTEGER NARGS, I
CHARACTER(LEN=80) BFF,VAR
! }}}
! subroutine body {{{

CALL IARGC_SUBR(NARGS)
I=0

IF (NARGS.GT.0) THEN
  DO WHILE(I.LT.NARGS) ; I=I+1 ; CALL GETARG_SUBR(I,BFF)
     SELECTCASE(BFF)
       CASE('-v')  ; CALL DISPLAY_VERSION(0) ; STOP
       CASE('-pv') ; CALL INITVARS ; CALL PRINTVARS ; STOP
       CASE('-h')  ; CALL PRINTHELP ; STOP
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
}

rca_OPTIM(){
rca_optim
}

rca_PATHSAMPLE(){
rca_pathsample
}

rca_optim(){
# {{{
cat << EOF
SUBROUTINE READ_CMD_ARGS

USE PORFUNCS
USE DV
USE COMMONS,ONLY: RINDEX

IMPLICIT NONE

INTEGER NARGS, I
CHARACTER(LEN=80) BFF,VAR

CALL IARGC_SUBR(NARGS)
I=0

IF (NARGS.GT.0) THEN
  DO WHILE(I.LT.NARGS) ; I=I+1 ; CALL GETARG_SUBR(I,BFF)
     SELECTCASE(BFF)
       CASE('-v')
         CALL DISPLAY_VERSION(0)
         STOP
       case default
         CALL GETARG_SUBR(i+1,var)
         SELECTCASE(BFF)
                CASE('-r')
                        READ(VAR,*) RINDEX 
                        WRITE(*,*) 'Run Index:',RINDEX
                CASE DEFAULT
         ENDSELECT
     ENDSELECT
         I=I+1
  ENDDO
ENDIF

END SUBROUTINE READ_CMD_ARGS
EOF
# }}}
}

rca_pathsample(){
# {{{
cat << EOF
SUBROUTINE RCA
! declarations {{{
USE PORFUNCS
USE DV
!USE COMMONS,ONLY : INFILE
use commons, only: infile,pforce

IMPLICIT NONE

INTEGER NARGS, I
CHARACTER(LEN=80) BFF,VAR
! }}}
CALL IARGC_SUBR(NARGS)
I=0

IF (NARGS.GT.0) THEN
  DO WHILE(I.LT.NARGS) ; I=I+1 ; CALL GETARG_SUBR(I,BFF)
     SELECTCASE(BFF)
       CASE('-v')
         CALL DISPLAY_VERSION(0)
         STOP
       case default
         CALL GETARG_SUBR(i+1,var)
         SELECTCASE(BFF)
                CASE('-f')
                        READ(VAR,*) PFORCE
                        WRITE(*,*) 'FORCE:',PFORCE
                CASE DEFAULT
         ENDSELECT
         !CALL MAKETEST(BFF)
       !CASE DEFAULT
         !INFILE=TRIM(ADJUSTL(BFF))
     ENDSELECT
         I=I+1
  ENDDO
ENDIF
END SUBROUTINE RCA

EOF
# }}}
}

rca_"$prog" 
