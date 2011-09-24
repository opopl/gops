!op226> Declarations {{{
      !USE commons
!      use COMMONS
      !USE V
      !USE PORFUNCS
      ! comm {{{
!      use MODMXATMS   ! NEEDED FOR charmm
      !USE modcharmm
!C       sf344> AMBER additions
      !USE modamber9, only : coords1,amberstr,amberstr1,mdstept,inpcrd,amberenergiest, nocistransdna, nocistransrna,
     !&                      uachiral, ligrotscale, setchiral, STEEREDMINT, SMINATOMA, SMINATOMB, SMINK, SMINKINC,
     !&                      SMINDISTSTART, SMINDISTFINISH, natomsina, natomsinb, natomsinc, atomsinalist, atomsinblist,
     !&                      atomsinclist, atomsinalistlogical, atomsinblistlogical, atomsinclistlogical, ligcartstep,
     !&                      ligtransstep, ligmovefreq, amchnmax, amchnmin, amchpmax, amchpmin, rotamert, rotmaxchange, 
     !&                      rotpselect, rotoccuw, rotcentre, rotcutoff 
      !USE modamber
      !USE PORFUNCS
      !USE MYGA_PARAMS
      ! }}}
      
      !COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,&
     !&                NERROR, IR, ECHO, LAST, CAT

      INTEGER ITEM, NITEMS, LOC, NCR, NERROR, IR, LAST, IX, J1, JP, NPCOUNT, NTYPEA, NPCALL, NDUMMY, INDEX, J2, J3
      LOGICAL END, SKIPBL, CLEAR, ECHO
      LOGICAL ::    CAT
      !CHARACTER WORD*16
      DOUBLE PRECISION ::   XX

      !INTEGER MOVABLEATOMINDEX
      !LOGICAL CAT, YESNO, PERMFILE
      !DOUBLE PRECISION XX, EPSAB, EPSBB, SIGAB, SIGBB
      !CHARACTER WORD*16,PBC*3,WORD2*10
      !COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      !DOUBLE PRECISION EAMLJA0, EAMLJBETA, EAMLJZ0, DUMMY
      !COMMON /EAMLJCOMM/ EAMLJA0, EAMLJBETA, EAMLJZ0
      !DOUBLE PRECISION SLENGTH, EPS
      !INTEGER NOK, NBAD
      !COMMON /BSNEW/ SLENGTH, NOK, NBAD, EPS
      !DOUBLE PRECISION EPS2, RAD, HEIGHT
      !COMMON /CAPS/ EPS2, RAD, HEIGHT

!     LOGICAL IGNOREBIN(HISTBINMAX), FIXBIN
!     COMMON /IG/ IGNOREBIN, FIXBIN
      !DOUBLE PRECISION    PMAX,PMIN,NMAX,NMIN,SIDESTEP
      !COMMON /AMBWORD/    PMAX,PMIN,NMAX,NMIN,SIDESTEP
      COMMON /PCALL/ NPCALL
      INTEGER NATOM, DMODE, NDUM
!
! These arrays should have dimension MXATMS
!!
      !DOUBLE PRECISION, ALLOCATABLE :: CHX(:), CHY(:), CHZ(:), CHMASS(:)
      !CHARACTER(LEN=1) DUMMYCH
      !CHARACTER(LEN=100) TOPFILE,PARFILE
      !CHARACTER(LEN=9) UNSTRING
      !DOUBLE PRECISION LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     !&                 HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN
      !DOUBLE PRECISION LJREPBL, LJATTBL, LJREPBN, LJATTBN, LJREPLN, LJATTLN

!!     DC430 >
      !DOUBLE PRECISION :: LPL, LPR

!!
!!       sf344> added stuff
!!
      !CHARACTER(LEN=10) check1
      !CHARACTER(LEN=1) readswitch
      !INTEGER iostatus, groupsize, groupatom,groupoffset,axis1,axis2
!op226> End declarations </begin> }}}

