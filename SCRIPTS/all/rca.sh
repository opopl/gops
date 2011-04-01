#!/bin/bash

export shd="`dirname $(readlink -f $0)`"
export this_script=` basename $0 `

prog="$1"

rca_GMIN(){
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

rca_OPTIM(){
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

rca_PATHSAMPLE(){
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

rca_"$prog" 
