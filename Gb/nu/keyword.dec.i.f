Cop226> Declarations {{{
      USE V
      USE MODMXATMS   
      USE PORFUNCS

      IMPLICIT NONE

      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST, IX, J1, JP, NPCOUNT, NTYPEA, NPCALL, NDUMMY, INDEX, J2, J3
      INTEGER MOVABLEATOMINDEX
      LOGICAL CAT, YESNO, PERMFILE
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      DOUBLE PRECISION XX, EPSAB, EPSBB, SIGAB, SIGBB
      LOGICAL END, SKIPBL, CLEAR, ECHO
      CHARACTER WORD*16,PBC*3,WORD2*10
      DOUBLE PRECISION EAMLJA0, EAMLJBETA, EAMLJZ0, DUMMY
      DOUBLE PRECISION SLENGTH, EPS
      INTEGER NOK, NBAD
      DOUBLE PRECISION EPS2, RAD, HEIGHT

      DOUBLE PRECISION    PMAX,PMIN,NMAX,NMIN,SIDESTEP

      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /EAMLJCOMM/ EAMLJA0, EAMLJBETA, EAMLJZ0
      COMMON /BSNEW/ SLENGTH, NOK, NBAD, EPS
      COMMON /CAPS/ EPS2, RAD, HEIGHT
      COMMON /AMBWORD/    PMAX,PMIN,NMAX,NMIN,SIDESTEP
      COMMON /PCALL/ NPCALL

      INTEGER NATOM, DMODE, NDUM
C
C These arrays should have dimension MXATMS
C
      DOUBLE PRECISION, ALLOCATABLE :: CHX(:), CHY(:), CHZ(:), CHMASS(:)
      CHARACTER(LEN=1) DUMMYCH
      CHARACTER(LEN=100) TOPFILE,PARFILE
      CHARACTER(LEN=9) UNSTRING
      DOUBLE PRECISION LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     &                 HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN
      DOUBLE PRECISION LJREPBL, LJATTBL, LJREPBN, LJATTBN, LJREPLN, LJATTLN

!     DC430 >
      DOUBLE PRECISION :: LPL, LPR

C
C       sf344> added stuff
C
      CHARACTER(LEN=10) check1
      CHARACTER(LEN=1) readswitch
      INTEGER iostatus, groupsize, groupatom,groupoffset,axis1,axis2
Cop226> End declarations </begin> }}}

