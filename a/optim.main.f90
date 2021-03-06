      PROGRAM OPTIM4
! declarations {{{
      USE COMMONS
      USE PORFUNCS
      USE KEY, ONLY: FILTHSTR,SEQ,NUMGLY,TARFL,CASTEPJOB,CP2KJOB,ONETEPJOB
      USE MODAMBER9, ONLY: AMBERSTR,AMBERSTR1,INPCRD,ATMASS1

      IMPLICIT NONE

      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST, J1
      LOGICAL VARIABLES,CASTEP,ONETEP,CP2K,DFTP,CPMD,END,CAT,SKIPBL,CLEAR,ECHO,AMBER,AMBERT,NABT,RINGPOLYMERT
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,NERROR, IR, ECHO, LAST, CAT
      DOUBLE PRECISION DUMMY
      CHARACTER ZDUM*5
      INTEGER J2, NELEMENTS, LSYS, NTYPE(105), IOS, NARG, FILTH, FILTH2
      CHARACTER FNAME*80, TSTRING*80
      CHARACTER(LEN=80) :: SYS
      CHARACTER WORD*16
      CHARACTER(LEN=10)  check
      CHARACTER(LEN=20) OTEMP, OSTRING, CSTRING
      CHARACTER(LEN=21) DSTRING1, DSTRING2
      CHARACTER(LEN=80) ARGSTRING, MYLINE
      LOGICAL AMH
      INTEGER :: NRES,I_RES,NOGLY
! }}}

      CASTEP=.FALSE.
      ONETEP=.FALSE.
      CP2K=.FALSE. 
      CPMD=.FALSE.
      VARIABLES=.FALSE.
      RINGPOLYMERT=.FALSE.
      AMBER=.FALSE.
      AMBERT=.FALSE.
      NABT=.FALSE.

      CALL READ_CMD_ARGS

include optim.read.odata.inc.f90
include optim.get.natoms.inc.f90 

      WRITE(*,'(A,I6)') ' getparams> Number of atoms (or variables)  determined as ',NATOMS

      CALL OPTIM(FILTH,FILTH2,ARGSTRING)

      STOP
      END
