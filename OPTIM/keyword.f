C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C  All the keywords possible for the odata file are contained here in 
C  alphabetical order. Initialisation statements procede the big IF block.
C
      !op226 <begin>
      SUBROUTINE KEYWORD(Q)
      ! Declarations {{{
      ! Modules  {{{
      USE COMMONS 
      USE KEY
      USE MODMEC
      USE MODTWOEND
      USE MODAMBER
      USE MODAMBER9, ONLY : COORDS1,IH,M04,PRMTOP,SALTCON,IGB,CUT,RGBMAX,
     & NOCISTRANSRNA,NOCISTRANSDNA,ATMASS1,CHECKCISTRANSALWAYS,CHECKCISTRANSALWAYSRNA,CHECKCISTRANSALWAYSDNA,
     & AMBERICT,AMBSTEPT,AMBIT,AMBPERTT, PERTHRESH, AMBOLDPERTT,AMBICDNEBT, 
     & AMBPDB_UNIT, AMBRST_UNIT, MDCRD_UNIT, MDINFO_UNIT,
     & KTWN, KTWNT, DUMPMODEN, UACHIRAL, NOPERMPROCHIRAL, FROZENAMBER
      USE MODNEB
      USE modmxatms   ! needed for CHARMM
      USE modcharmm   
      USE MODUNRES
      USE KEYNEB, NNNIMAGE=>NIMAGE
      USE KEYCONNECT
      USE MODMEC
      USE MODGUESS
      USE PORFUNCS
      USE PYMODULE, only : BOXLX,BOXLY,BOXLZ
      USE msevb_common, ONLY: shellsToCount, maxHbondLength, minHbondAngle, OOclash_sq, printCoefficients
      use wc
      use binaryio
      USE GSDATA, ONLY : CUBSPLT, GSUPDATE,
     $     GSGROWTOL, GSMXSTP,GSCONV, REPARAMTOL, EVOLVESTRINGT,
     $     GSITERDENSITY, FIXATMS, MAXLENPERIM,
     $     HESSGRAD, GSMAXTOTITD, MAXGROWSTEPS, GSDGUESS,
     $     NOLBFGS, PREROTATE, GSTANTYPE=>TANTYPE
      USE CUBSPLSTRING, ONLY : ARCTOL, DQAGKEY
      USE INTCOMMONS, ONLY : NATINT, INTNEWT, BBCART, INTINTERPT, INTERPSIMPLE,
     $     INTMINPERMT, INTERPCHOICE, NINTIM, CARTRESSTART, INTPARFILE,
     $     MINBACKTCUT, INTERPBACKTCUT, PRINTCOORDS, DESMINT, NURINGS, URINGS,
     $     NUBONDS, UBONDS, USEPARFILE, CHICDNEB, OLDINTMINPERMT, 
     $     GLYCART, INTDISTANCET
!     USE BENCHMARKS, ONLY : MINBMT, MINBMNSAMP
C     MCP
      USE AMHGLOBALS!}}}


      IMPLICIT NONE
      
      DOUBLE PRECISION ::  Q(3*NATOMS)

      ! local parameters {{{

      INTEGER NDUM
      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST, NTYPEA, J1, J2, J3
      COMMON /BUFINF/ ITEM, NITEMS, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      DOUBLE PRECISION :: XX, EPSAB, EPSBB, SIGAB, SIGBB, RANDOM, DPRAND
      LOGICAL END, SKIPBL, CLEAR, ECHO, CAT, CISTRANS, RBSYMTEST
      CHARACTER WORD*25, WW*20, PBC*3
      CHARACTER WORD2*25
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA

      INTEGER NATOM, DMODE
      DOUBLE PRECISION CHX(MXATMS), CHY(MXATMS), CHZ(MXATMS), CHMASS(MXATMS)
      DOUBLE PRECISION DPERT
      DOUBLE PRECISION CHPMIN, CHPMAX, CHNMIN, CHNMAX
      INTEGER ISEED

      DOUBLE PRECISION UNRX(NATOMS), UNRY(NATOMS), UNRZ(NATOMS) ! UNRES
      DOUBLE PRECISION DUMMY1(NATOMS)

      DOUBLE PRECISION SLENGTH, EPS
      INTEGER NOK, NBAD
      COMMON /BSNEW/ SLENGTH, NOK, NBAD, EPS
      DOUBLE PRECISION GSQSCALE, GSTHRESH
      INTEGER NSPECIAL, NALLOW, NINFO
      COMMON /G2/ GSTHRESH, GSQSCALE, NSPECIAL, NALLOW, NINFO
      LOGICAL CUBIC
      COMMON /CUB/ CUBIC
      LOGICAL PATHT, DRAGT
      INTEGER NPATHFRAME
      COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
      LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST
      INTEGER NCONNECT
      DOUBLE PRECISION TEMPERATURE, HRED
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
      INTEGER NMOVE
      COMMON /HS/ NMOVE
C     DOUBLE PRECISION REPELTS(3*NATOMS,100), REPELPUSH 
C     INTEGER NREPELTS, REPELFROM
C     LOGICAL REPELTST, REPEL
C     COMMON /OTS/ NREPELTS, REPELTST, REPELPUSH, REPEL, REPELFROM
      INTEGER ISTAT, NDUMMY
      DOUBLE PRECISION STOPDISP
      LOGICAL STOPDISPT, PERMFILE
      COMMON /STOPD/ STOPDISP, STOPDISPT
      DOUBLE PRECISION CAPSRHO, CAPSEPS2, CAPSRAD, HEIGHT
      COMMON /CAPS/ CAPSRHO, CAPSEPS2, CAPSRAD, HEIGHT
      CHARACTER(LEN=20) OSTRING, OTEMP
      CHARACTER(LEN=20) :: PINFOSTRING
      CHARACTER(LEN=5) :: TEMPSTRING
      CHARACTER(LEN=9) UNSTRING
      CHARACTER(LEN=1) DUMMYCH
      CHARACTER(LEN=100) TOPFILE,PARFILE
      DOUBLE PRECISION LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     &               HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN
      INTEGER IGBNAB     ! sf344
C    LOCAL AMH VARIABLES
      INTEGER NRES_AMH, I_RES, GLY_COUNT
C      CHARACTER(LEN=5) TARFL
C      DOUBLE PRECISION X, Y, Z
      INTEGER :: GROUPCENTRE
      DOUBLE PRECISION :: GROUPRADIUS,DISTGROUPX2,DISTGROUPY2,DISTGROUPZ2,DISTGROUPCENTRE
      CHARACTER (LEN=2) :: FREEZEGROUPTYPE
      LOGICAL :: FREEZEGROUPT
      ! }}}
! }}}
    
      ! </begin>
      ! Initialize variables/keywords {{{
      ! <init>
      DESMAXEJUMP = HUGE(1.0D0)
      DESMAXAVGE = HUGE(1.0D0)
      UNSTRING='UNDEFINED'
      WELCH=.FALSE.
      TOSI=.FALSE.
      TOSIC6=.FALSE.
      SIO2T=.FALSE.
      SIO2C6T=.FALSE.
      TOSIPOL=.FALSE.
      ALLOCATE(TAGFAC(NATOMS),TAGNUM(NATOMS))
      TAGFAC(1:NATOMS)=1.0D0
      TAGNUM(1:NATOMS)=0
      NTAG=0
      TAGT=.FALSE.
  
      REPEL=.FALSE.

      INR=-1

      INTMINPERMT =.FALSE. !msb50 internal permutation
      NOPERMPROCHIRAL = .FALSE.
      BFGSMINT=.FALSE.
      GMAX=0.001D0
      BFGSTST=.FALSE.
      HYBRIDMINT=.FALSE.
      REOPT=.FALSE.
      NOHESS=.FALSE.
      NOFRQS=.FALSE.
      NOIT=.FALSE.
      NEVL=100
      NEVS=500
      NINTS=0
      NBFGSMAX1=10
      NBFGSMAX2=100
      CEIG=1.0D-10  ! changed to a small default to make people change it!
      CHECKINDEX=.FALSE.
      CHECKCONT=.FALSE.
      BFGSSTEP=.FALSE.
      EXTRASTEPS=0.0D0

      DCHECK=.TRUE.

      PRESSURE=.FALSE.
      PV=.FALSE.
      PVTS=.FALSE.
      FRACTIONAL=.FALSE.
      PRESS=0.0D0
      PVCONV=1.0D-3
      PVTOL=1.0D60
      PVSTEPS=100
      NBOXTS=1

      VARIABLES=.FALSE.
      NZERO=0
      EVCUT=0.0D0

      GAUSSIAN=.FALSE.
      CADPAC=.FALSE.
      GAMESSUS=.FALSE.
      GAMESSUK=.FALSE.
      CASTEP=.FALSE.
      CASTEPJOB=''
      ONETEP=.FALSE.
      ONETEPJOB=''
      CP2K=.FALSE. 
      CP2KJOB=''
      DFTP=.FALSE.
      CPMD=.FALSE.
      CPMDC=.FALSE.
      PARALLEL=.FALSE.
      NPROC='1'
      DFTBT=.FALSE.
      CPMD_COMMAND='/home/trj25/bin/cpmd.x'
      SCORE_QUEUE=.FALSE.

!     DC430 >

      DBPT     = .FALSE.
      DBPTDT   = .FALSE.
      LWOTPT   = .FALSE.
      GBT      = .FALSE.
      GBDT     = .FALSE.
      MSSTOCKT = .FALSE.
      NCAPT    = .FALSE.
      NTIPT    = .FALSE.
      PAHAT    = .FALSE.
      PATCHYDT = .FALSE.
      PYGT     = .FALSE.
      RADIFT   = .FALSE.
      RBAAT    = .FALSE.
      RBSYMT   = .FALSE.
      STOCKAAT = .FALSE.
      UNIAXT   = .FALSE.

!     -----------------------

      ISTCRT=10
      IPRNT=0
      IVEC=0
      IVEC2=0

      MXSTP=0.2D0
      MINMAX=0.01D0
      MAXMAX=0.5D0
      MAXBFGS=0.2D0
      MAXXBFGS=0.2D0
      MAXMBFGS=0.2D0
      MAXNEBBFGS=0.2D0
      MAXINTBFGS=0.2D0

      DTEST=.FALSE.

      MASST=.FALSE.

      VALUEST=.TRUE.
      EFSTEPST=.FALSE.
      EFSTEPS=1
      NVALUES=20
      NSTEPS=1
      BFGSSTEPS=1
      DUMPV=.FALSE.
      ALLSTEPS=.FALSE.
      ALLVECTORS=.FALSE.
      MWVECTORS=.FALSE.
      READV=.FALSE.

      PGRAD=.FALSE.
      NGRADIENTS=1

      VECTORST=.FALSE.
      NVECTORS=1

      SUMMARYT=.TRUE.
      NSUMMARY=20

      ADMT=.FALSE.
      NADM=20

      CONVU=1.0D-5
      CONVR=1.0D-5
      INDEXT=.TRUE.

      SYMCUT=0.001D0
      TOLD=0.0001D0
      TOLE=0.0001D0
      NHCHECK=6

      TRAD=2.0
      RESIZE=1.0D0

      RTEST=.FALSE.
      JZ=0.0D0
      OMEGA=0.0D0

      PUSHOFF=0.01D0
      PUSHCUT=1.0D-5

      BINARY=.FALSE.
      NSTEPMIN=0

      HUPDATE=.FALSE.
      NSTHUP=0
      INTHUP=0
      PHIG=0.0D0
      READHESS=.FALSE.

      SHIFTV=1.0D6

      NORESET=.FALSE.

      MUPDATE=4
      XMUPDATE=4
      MMUPDATE=4
      NEBMUPDATE=4
      INTMUPDATE=100
      GSUPDATE = 4
      GCUPDATE=4
      DGUESS=0.1D0
      XDGUESS=0.1D0
      NEBDGUESS=0.001D0
      INTDGUESS=0.001D0
      GSDGUESS = 0.001D0
      AMBER=.FALSE.
      AMBERT=.FALSE.
      NABT=.FALSE.
      NOCISTRANSRNA=.FALSE.
      NOCISTRANSDNA=.FALSE.
      CHECKCISTRANSALWAYS=.FALSE.
      CHECKCISTRANSALWAYSRNA=.FALSE.
      CHECKCISTRANSALWAYSDNA=.FALSE.
      UACHIRAL=.FALSE.
      
      FAKEWATER=.FALSE.

      CHRMMT=.FALSE.
      REDUCEDBONDLENGTHT=.FALSE.
      BLFACTOR=1.D0
      ACESOLV=.FALSE.
      ACEUPSTEP=50
      TWISTDIHET=.FALSE.
      PERTDIHET=.FALSE.
      CHPMAX=0.5d0
      CHPMIN=0.25d0
      CHNMAX=1.0d0
      CHNMIN=0.d0
      ISEED=0
      TOMEGAC=.FALSE.
      TSIDECHAIN=.FALSE.
      INTMINT=.FALSE.
      IMINCUT=0.0D0
      GUESSTST=.False.
      CALCDIHE=.False.
      TRYNEB=.FALSE.
      NOCISTRANS=.FALSE.
      CISTRANS=.FALSE.
      CHECKOMEGAT=.FALSE.
      MINOMEGA=150.D0
      CHECKCHIRALT=.FALSE.
      NORANDOM=.FALSE.
      RANDOMCUTOFF=0.d0
!     GUESSTHRESH=1.0D100
      ENDHESS=.FALSE.
      NENDHESS=0
      ENDNUMHESS=.FALSE.
      NPERMDIHE=0
      TWISTTYPE=0
      NGUESS=3
      FAILT=.FALSE.
      OSASAT=.FALSE.
      RPRO=1.4D0
      ODIHET=.FALSE.

C unres stuff
      UNRST=.FALSE.
      CONSECT=.FALSE.
      STARTRES=0
      ENDRES=0
      NUMSEC=0
C 
C AMH  stuff
      AMHT=.FALSE.

      FREEZE=.FALSE.
      FREEZEGROUPT=.FALSE.
      FREEZEGROUPTYPE='GT'
      FREEZERES=.FALSE.
      NFREEZE=0
      DO J1=1,NATOMS
         FROZENRES(J1)=.FALSE.
         FROZEN(J1)=.FALSE.
      ENDDO
      ALLOCATE(DUMPMODEN(3*NATOMS))
      DO J1=1,3*NATOMS
         DUMPMODEN(J1)=.FALSE.
      ENDDO
      KEEPINDEX=.FALSE.
      BSMIN=.FALSE.
      RKMIN=.FALSE.
      SLENGTH=0.0D0
      FIXAFTER=-1
      HINDEX=1
      NOK=0
      NBAD=0
      EPS=1.0D-3
      CONTAINER=.FALSE.
      FIXD=.FALSE.
      T12FAC=1.1D0
      PRINTPTS=.FALSE.
      GRADSQ=.FALSE.
      GSQSCALE=1.0D0
      NSPECIAL=-1
      NALLOW=100
      NINFO=0
      GSTHRESH=0.0D0
      TWOD=.FALSE.
      DOUBLET=.FALSE.
      TWOENDS=.FALSE.
      FSTART=1.0D0
      FINC=1.0D0
      RMSTWO=0.001D0
      NTWO=100
      NTWOITER=25
      TWOEVAL=0.0D0
      PATHT=.FALSE.
      STOPFIRST=.FALSE.
      CONNECTT=.FALSE.
      DUMPPATH=.FALSE.
      DUMPALLPATHS=.FALSE.
      READPATH=.FALSE.
      CALCRATES=.FALSE.
      TEMPERATURE=1.0D0
      KTWN=207.11
      KTWNT=.FALSE.
      HRED=1.0D0
      NCONNECT=100
      NEWNEBT=.FALSE.
      NEBT=.FALSE.
      NEWCONNECTT=.False.
      SQVVGuess=.FALSE.
      SQVVGuessRMSTol=2.0D0
      NIterSQVVGuessMax=300
      DEBUG=.FALSE.
      CHDEBUG=.FALSE.
      EDEBUG=.FALSE.
      NIMAGE=10
      RMSNEB=0.1
      DTHRESH=2.0D0
      NSTEPNEB=1
      NEBMAG=0
      NMOVE=1
      NPATHFRAME=0
      FRAMEEDIFF=0.0D0
      CUBIC=.FALSE.
      BULKT=.FALSE.
      SDT=.FALSE.
      SDOXYGEN=0
      SDHYDROGEN=0
      SDCHARGE=0
      BOWMANT=.FALSE.
      BOWMANPES=2
      BOWMANDIR='~/svn/OPTIM/source/Bowman/coef-3b/'
      RATIOS=.FALSE.
      QSPCFWT=.FALSE.
      QTIP4PFT=.FALSE.


      ! EFK: growing strings and freezing nodes
      GROWSTRINGT = .FALSE.
      NOLBFGS = .FALSE.
      HESSGRAD = .FALSE.
      ARCTOL = 1.0D-4
      DQAGKEY = 6
      DESMDEBUG = .FALSE.
      GSMAXTOTITD = -1
      MAXGROWSTEPS = 1.0D3
      EVOLVESTRINGT = .FALSE.
      FREEZENODEST = .FALSE.
      FIXATMS = .FALSE.
      PREROTATE = .FALSE.
      CUBSPLT = .FALSE.
      MAXLENPERIM = 100.0D0
      GSTANTYPE = 1
      REPARAMTOL = 0.75
      GSGROWTOL = 0.25
      GSCONV = 1.0D-3
      GSMXSTP = 0.1
      STOCKT=.FALSE.
      STOCKSPIN = .FALSE.
      STOCKZTOL = 1.0D-4
      STOCKMAXSPIN = 20
      GEOMDIFFTOL=1.0D-1
      EDIFFTOL=1.0D-6
      NSECDIAG=1

! MSEVB parameters

      shellsToCount = 3
      maxHbondLength = 2.5d0
      minHbondAngle = 130.0d0
      OOclash_sq = 4.41d0 ! 2.1^2
      printCoefficients = .FALSE.

      NEBRESEEDT=.FALSE.
      NEBRESEEDINT=100
      NEBRESEEDEMAX=1.0D100
      NEBRESEEDBMAX=1.0D100
      NEBRESEEDDEL1=1.0D5
      NEBRESEEDDEL2=1.0D5
      NEBRESEEDPOW1=2
      NEBRESEEDPOW2=10
      ADDREPT=.FALSE.

      INTLJT=.FALSE.
      INTLJSTEPS=1000
      INTLJTOL=1.0D-3
      INTLJDEL=0.1D0
      INTLJEPS=1.0D0

      INTCONSTRAINTT=.FALSE.
      INTCONSTRAINTTOL=0.1D0
      INTCONSTRAINTDEL=1.0D5
      INTCONSTRAINTREP=1.0D0
      INTCONSTRAINREPCUT=20.0D0
      INTCONSTEPS=1000
      INTRELSTEPS=100
      INTCONFRAC=1.0D-4
      INTREPSEP=0
      INTCONSEP=10000
      MAXINTIMAGE=75
      MAXCONUSE=3
      MAXCONE=0.1D0
      INTRMSTOL=1.0D-3
      INTIMAGE=30
      CHECKCONINT=.FALSE.
      DUMPINTXYZ=.FALSE.
      DUMPINTEOS=.FALSE.
      DUMPINTXYZFREQ=100
      DUMPINTEOSFREQ=100

      CONPOTT=.FALSE.
      CPCONSTRAINTTOL=0.1D0
      CPCONSTRAINTDEL=1.0D5
      CPCONSTRAINTREP=1.0D0
      CPCONSTRAINREPCUT=20.0D0
      CPCONFRAC=1.0D-4
      CPREPSEP=0
      CPCONSEP=10000
C
C  UNDOCUMENTED keywords/parameters
C
      MORPHT=.FALSE.
      GREATCIRCLET=.FALSE.
      MAXTSENERGY=1.0D100
      MAXBARRIER=1.0D100
      MAXMAXBARRIER=1.0D100
      ReoptimiseEndpoints=.False.
      ANGLEAXIS=.FALSE.
      NFAILMAX=2
      NATBT=.FALSE.
      READSP=.FALSE.
      DUMPSP=.FALSE.
      TIMELIMIT=HUGE(TIMELIMIT)
      RIGIDBODY=.FALSE.
      STOPDISPT=.FALSE.
      NCHENCALLS=0

      REPELTST=.FALSE.
      NREPELTS=0
      REPELPUSH=0.1D0

      DRAGT=.FALSE.

      LANCZOST=.FALSE.
      ACCLAN=1.0D-8
      SHIFTLAN=1.0D-2
      CUTLAN=-1.0D0

      GFRACTION=0.0D0
      MFRACTION1=0.0D0
      MFRACTION2=0.0D0
      FTEST=.FALSE.
      GALPHA=6.0D0
      MALPHA1=6.0D0
      MALPHA2=6.0D0

      FIELDT=.FALSE.
      OHT=.FALSE.
      IHT=.FALSE.
      TDT=.FALSE.
      D5HT=.FALSE.
      FOH=0.0D0
      FIH=0.0D0
      FTD=0.0D0
      FD5H=0.0D0
      MAXERISE=1.0D-10
      XMAXERISE=1.0D-3
      INTEPSILON=1.0D-6

      EFIELD=0.0D0
      COLDFUSIONLIMIT=-1.0D6
      BLNT=.FALSE.
      DUMPDATAT=.FALSE.
      LOWESTFRQT=.FALSE.
      REDOPATH=.FALSE.
      REDOFRAC=0.5D0
      REDOK=0.0D0
      REDOKADD=.FALSE.
      REDOPATHNEB=.FALSE.
      REDOBFGSSTEPS=100
      REDOPATHXYZ=.FALSE.
      REALIGNXYZ=.FALSE.
      PERMDIST=.FALSE.
      LOCALPERMDIST=.FALSE.
      LPDGEOMDIFFTOL=0.3D0
      RBCUTOFF=4.0D0
      NRBTRIES=1
      ALLOCATE(BESTPERM(NATOMS))
      PERMDISTINIT=.FALSE.
      NEBK=1.0D0 ! changed DJW 14/5/08
      NEBKINITIAL=1.0D0
      NEBKFINAL=1.0D0
      NEBFACTOR=1.01D0
      KADJUSTFRQ=-1
      KADJUSTTOL=0.1D0
      KADJUSTFRAC=0.1D0  ! ten percent

      BHDEBUG=.FALSE.
      BHDISTTHRESH=1.0D0
      BHMAXENERGY=1.D100
      BHINTERPT=.FALSE.
      BHACCREJ=0.5D0
      BHSTEPSIZE=0.4D0
      BHCONV=0.01D0
      BHSTEPS=1000
      BHTEMP=1.0D0
      BHINTERPUSELOWEST=.FALSE.
      BHCHECKENERGYT=.FALSE.
      BHSTEPSMIN=0
      BHK=1.0D0
      ICINTERPT=.FALSE.
      CHBIT=.FALSE.
      BISECTT=.FALSE.
      BISECTMAXENERGY=1.D100
      BISECTDEBUG=.FALSE.
      BISECTSTEPS=1
      BISECTMINDIST=1.0
      BISECTMAXATTEMPTS=5
      CHRIGIDT=.FALSE.
      PTRANS=0.D0
      TRANSMAX=0.D0
      PROT=0.D0
      ROTMAX=0.D0
      BBRSDMT=.FALSE.
      AMBERICT=.FALSE.
      AMBSTEPT=.FALSE.
      AMBICDNEBT = .FALSE.
      AMBPERTT=.FALSE.
      AMBOLDPERTT=.FALSE.
      AMBIT = .FALSE. 

      DIJKSTRALOCAL=1.0D0
      DNEBSWITCH=-1.0D0
      PATHSDSTEPS=-1
      NUSEEV=-1
      ACKLANDID=5

      RINGPOLYMERT=.FALSE.
      RPSYSTEM='     '
      RPIMAGES=1
      RPBETA=1.0D0
      GRAD4T=.FALSE.
      RPCYCLICT=.TRUE.
      RPFIXT=.FALSE.
      EYTRAPT=.FALSE.
      BFGSTSTOL=0.0001D0
! sf344
      PYGPERIODICT=.FALSE.
      PYBINARYT=.FALSE.
      LJSITE=.FALSE.
      BLJSITE=.FALSE.
      LJSITECOORDST=.FALSE.
      LJSITEATTR=.FALSE.
      PCUTOFF=999.0D0
      CLOSESTALIGNMENT=.FALSE.
      DF1T=.FALSE.
      PULLT=.FALSE.
      IMSEPMIN=0.0D0
      IMSEPMAX=HUGE(1.0D0)
! </init>
! }}}
 
      ! open odata; read WORD {{{
      ! <read>
      IF (FILTH2.EQ.0) THEN
         OPEN (5,FILE='odata',STATUS='OLD')
      ELSE
         WRITE(OTEMP,*) FILTH2
         WRITE(OSTRING,'(A)') 'odata.' // TRIM(ADJUSTL(OTEMP))
         OPEN (5,FILE=OSTRING,STATUS='OLD')
      ENDIF

190   CALL INPUT(END)
      IF (.NOT. END) THEN
        CALL READU(WORD)
      ENDIF
C
C  POINTS - keyword at the end of the list of options after which
C           the Cartesian coordinates follow. Must be present unless VARIABLES or RINGPOLYMER
C           is present instead. MACHINE keyword overrides POINTS. If MACHINE is
C           true coordinates that were read from odata file will be overwritten
C           with coordinates from a direct access file, in which case section of
C           odata file after POINTS keyword is used only to read in the labels. (SAT)
C
      IF (END.OR.WORD.EQ.'STOP'.OR.WORD.EQ.'POINTS') THEN
        RETURN
      ENDIF

      IF (WORD.EQ.'    ' .OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'
     &                          .OR. WORD .EQ. '\\') THEN 
         GOTO 190
         ! </read>
         ! }}}
         
        ! loop over the rest of keywords {{{
        ! <kwd>
C
C  Enforce flatland.
C
      ELSE IF (WORD .EQ. '2D') THEN
         TWOD=.TRUE.
C
C  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C
C  bs360: ACE is to be used together with CHARMM and the ACE solvent model, 
C  it makes sure that the Born radii are regularly updated
C
      ELSE IF (WORD.EQ.'ACE') THEN
          ACESOLV=.TRUE.
          IF (NITEMS.GT.1) CALL READI(ACEUPSTEP)
C
C  Adjust NEB force constant values between different images on-the-fly in order
C  to try and equispace them.
C
      ELSE IF (WORD.EQ.'ADJUSTK') THEN
          CALL READI(KADJUSTFRQ)
          CALL READF(KADJUSTTOL)
          CALL READF(KADJUSTFRAC)
C
C  ADM [OFF | ON n] prints the atomic distance matrix every n 
C                   if switched on                 cycles       - default n=20      
C
      ELSE IF (WORD .EQ. 'ADM') THEN
         ADMT=.TRUE.
         CALL READI(NADM)
C
C  Ackland embedded atom metal potentials.
C
      ELSE IF (WORD.EQ.'ACKLAND') THEN
         CALL READI(ACKLANDID) ! default is 5 = Fe

C  Keywork ALPHA enables exponent values to be set for the averaged
C  Gaussian and Morse potentials. All defaults = 6.
C
      ELSE IF (WORD.EQ.'ALPHA') THEN
         CALL READF(GALPHA)
         IF (NITEMS.GT.2) THEN
            CALL READF(MALPHA1)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(MALPHA2)
         ENDIF
C
C SAT: ALLPOINTS turns on printing of coordinates to file points for intermediate steps.
C This is the default.
C
      ELSE IF (WORD.EQ.'ALLPOINTS') THEN
         PRINTPTS=.TRUE.
C
C  AMBER stuff
C
      ELSE IF (WORD.EQ.'AMBER') THEN
         AMBER=.TRUE.
         CALL APARAMS
         CALL AREAD
         NATOMS=ATOMS
         DO J1=1,NATOMS
            Q(3*(J1-1)+1)=x(J1)
            Q(3*(J1-1)+2)=y(J1)
            Q(3*(J1-1)+3)=z(J1)
         ENDDO
         t=0
         ang=0
         imp=0
         count=0
C MCP
      ELSE IF (WORD.EQ.'AMH') THEN
         WRITE(6,*)'USING AMH ENERGIES FORCES' 
         WRITE(6,*)'CALCULATE ENERGY AND FORCE TABLES  '
         AMHT=.TRUE.
         WRITE(6,*)'AMH FLAG ', AMHT
         WRITE(6,*)'AMH NATOMS ',  NATOMS
         IF (DEBUG) WRITE(6,*)'Entering WALESAMH_INITIAL'

         CALL WALESAMH_INITIAL

         IF (DEBUG)WRITE(6,*)'Leaving WALESAMH_INITIAL'
         IF (DEBUG)WRITE(6,*)'TARFL ',TARFL

           OPEN(30,FILE='proteins/'//TARFL,STATUS='OLD')
           READ(30,*)
           READ(30,*)NRES_AMH
           IF (NRES_AMH.GT.500) THEN
              WRITE(6,*) 'FAILURE NRES_AMH GR THAN 500 CONNECTODATA'
              STOP
           ENDIF
           READ (30,25)(SEQ(I_RES),I_RES=1,NRES_AMH)
25         FORMAT(25(I2,1X))
           CLOSE(30)

           WRITE(6,*)'NRES ',NRES_AMH
           NRES_AMH_TEMP=NRES_AMH
 
         DO J1=1,NRES_AMH
             Q(9*(J1-1)+1)=X_MCP(9*(J1-1)+1)
             Q(9*(J1-1)+2)=X_MCP(9*(J1-1)+2)
             Q(9*(J1-1)+3)=X_MCP(9*(J1-1)+3)
             Q(9*(J1-1)+4)=X_MCP(9*(J1-1)+4)
             Q(9*(J1-1)+5)=X_MCP(9*(J1-1)+5)
             Q(9*(J1-1)+6)=X_MCP(9*(J1-1)+6)
             Q(9*(J1-1)+7)=X_MCP(9*(J1-1)+7)
             Q(9*(J1-1)+8)=X_MCP(9*(J1-1)+8)
             Q(9*(J1-1)+9)=X_MCP(9*(J1-1)+9)
         ENDDO

         t=0
         ang=0
         imp=0
         count=0
!
! sf344> start of AMBER-related keywords
!
      ELSE IF (WORD.EQ.'AMBER9') THEN
        AMBERT=.TRUE.
!
! csw34> if FREEZERES specified, populate the FROZEN array with A9RESTOATOM
!
        IF (FREEZERES) CALL A9RESTOATOM(FROZENRES,FROZEN,NFREEZE)
        IF ((PERMDIST.OR.LOCALPERMDIST).AND.(NPERMSIZE(1).EQ.NATOMS)) THEN
           PRINT '(A)','keyword> ERROR - PERMDIST or LOCALPERMDIST is specfied for AMBER, but there is no perm.allow file present'
           STOP
        ENDIF
        IF (FREEZEGROUPT) THEN
! Write a list of FROZEN atoms for use in an (o)data file
           OPEN(UNIT=4431,FILE='frozen.dat',STATUS='UNKNOWN',FORM='FORMATTED')
           DO J1=1,NATOMS
!
! Work out the distance from GROUPCENTRE to the current atom J1
! 
              DISTGROUPX2=(COORDS1(3*GROUPCENTRE-2)-COORDS1(3*J1-2))**2
              DISTGROUPY2=(COORDS1(3*GROUPCENTRE-1)-COORDS1(3*J1-1))**2
              DISTGROUPZ2=(COORDS1(3*GROUPCENTRE  )-COORDS1(3*J1  ))**2
              DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
! If working in GT mode (default), FREEZE all atoms >GROUPRADIUS from the GROUPCENTRE atom
              IF((FREEZEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
                 NFREEZE=NFREEZE+1
                 FROZEN(J1)=.TRUE.
                 WRITE(4431,'(A,I6)') 'FREEZE ',J1
! IF working in LT mode, FREEZE all atoms <GROUPRADIUS from the GROUPCENTRE atom
              ELSE IF((FREEZEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
                 NFREEZE=NFREEZE+1
                 FROZEN(J1)=.TRUE.
                 WRITE(4431,'(A,I6)') 'FREEZE ',J1
              END IF
           END DO
           CLOSE(4431)     
        ENDIF
!
! csw34> A copy of the FROZEN array called FROZENAMBER is created to be passed through to AMBERINTERFACE
!
        ALLOCATE(FROZENAMBER(NATOMS))
        FROZENAMBER(:)=FROZEN(:)
        IF(.NOT.ALLOCATED(ATMASS)) ALLOCATE(ATMASS(NATOMS))
        ATMASS(1:NATOMS) = ATMASS1(1:NATOMS)
        DO i=1,3*NATOMS
                Q(i) = COORDS1(i)
        END DO
! save atom names in array zsym
        do i=1,natoms
                zsym(i) = ih(m04+i-1)
        end do
        RETURN
! initialise unit numbers
        ambpdb_unit=1110
        ambrst_unit=1111
        mdinfo_unit=1112
        mdcrd_unit =1113
 
      ELSE IF (WORD.EQ.'AMBERIC') THEN
        PRINT*, "amberic"
        AMBERICT = .TRUE.
        IF (NITEMS .GT. 1) THEN
           CALL READA(WORD2)
           IF (WORD2.EQ.'BACKBONE')  THEN
              PRINT*, "backbone interpolated"
              AMBIT = .TRUE.
           ELSE 
              PRINT*, "keyword error in amberic"
              RETURN
           ENDIF
        ENDIF
 
      ELSE IF (WORD.eq.'AMBERSTEP') THEN
        PRINT*, "amberstept"
        AMBSTEPT = .TRUE. 
 
      ELSE IF (WORD.eq.'AMBPERTOLD') THEN
        PRINT*, "original perturbation scheme"
        AMBOLDPERTT = .TRUE.

      ELSE IF (WORD.eq. 'AMBPERTONLY') THEN
        AMBPERTT = .TRUE.
        CALL READF(PERTHRESH)        
        PRINT*, "amber pertonly, perthresh", perthresh

      ELSE IF (WORD.eq. 'AMBICDNEB') THEN
        AMBICDNEBT = .TRUE.

      ELSE IF (WORD.eq.'NAB') THEN
        IF (FREEZERES) CALL A9RESTOATOM(FROZENRES,FROZEN,NFREEZE)
        NABT=.TRUE.
        DO i=1,3*NATOMS
                Q(i) = COORDS1(i)
        END DO
! save atom names in array zsym
        do i=1,natoms
                zsym(i) = ih(m04+i-1)
        end do
        IF(.NOT.ALLOCATED(ATMASS)) ALLOCATE(ATMASS(NATOMS))
! for the NAB interface, ATMASS is also set up in mme2wrapper, and that setting 
! overrides the one from below. However, both originate from the same prmtop file, 
! so they should be the same. ATMASS is being assigned here so that it's somewhat consistent
! with the AMBER interface.
        ATMASS(1:NATOMS) = ATMASS1(1:NATOMS)
        WRITE(prmtop,'(A)') 'coords.prmtop'
        igbnab=igb
        if(igb==6) igbnab=0     ! this is also in vacuo, but NAB doesn't understand igb=6!
        CALL MMEINITWRAPPER(trim(adjustl(prmtop)),igbnab,saltcon,rgbmax,sqrt(cut))
        RETURN

      ELSE IF (WORD.eq.'DF1') THEN
        DF1T=.TRUE.

      ELSE IF (WORD.eq.'DUMPSTRUCTURES') THEN
        DUMPSTRUCTURES=.TRUE.
        WRITE(*,'(A)') ' keywords> Final structures will be dumped in different formats (.rst, .xyz, .pdb)'
!
! Distinguish between old C of M/Euler and new angle/axis coordinates for
! rigid body TIP potentials
!
      ELSE IF (WORD.EQ.'ANGLEAXIS') THEN
         ANGLEAXIS=.TRUE.
C
C  Growing string arc tolerance.
C
      ELSE IF (WORD.EQ.'ARCTOL') THEN
         CALL READF(ARCTOL)
C
C  Specifies the highest symmetry axis to search for in routine {\bf symmetry}; default is six.
C
      ELSE IF (WORD .EQ. 'AXIS') THEN
         CALL READI(NHCHECK)
C
C  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
      ELSE IF (WORD.EQ.'BBCART') THEN
         BBCART = .TRUE. ! use cartesians for backbone

      ELSE IF (WORD.EQ.'BBRSDM') THEN
C
C  BBSDM minimiser.
C
         BBRSDMT = .TRUE. 
         CALL READF(BBRGAM)
         CALL READF(BBREPS)
         CALL READF(BBRSIGMA1)
         CALL READF(BBRSIGMA2)
         CALL READI(BBRM)
         CALL READF(BBRALPHA)
         CALL READF(BBRCONV)
         CALL READI(BBRSTEPS)

      ELSE IF (WORD.EQ.'BFGSCONV') THEN
C
C  Turn on LBFGS gradient minimization. GMAX is the convergence
C  criterion for the RMS gradient, default 0.001.
C  For BFGSTS NEVL and NEVS are the maximum iterations allowed in the searches for 
C  the largest and smallest eigenvectors, respectively and NBFGSMAX1 is the largest
C  number of BFGS steps allowed in the subsequent restricted minimization.
C  If the negative eigenvalue appears to have converged then NBFGSMAX2 steps
C  are allowed in the tangent space.
C  CONVU is used to determine convergence in such runs and BFGSCONV can be used
C  to set GMAX, the convergence criteria for the subspace optimization.
C
C  IF REOPT is true the smallest Hessian eigenvector is redetermined after the
C  EF step before the tangent space minimisation.
C
         IF (NITEMS.GT.1) THEN
            CALL READF(GMAX)
         ENDIF

      ELSE IF (WORD.EQ.'BFGSMIN') THEN
C
C  instructs the program to perform an LBFGS minimisation.
C  {\it gmax\/} is the convergence criterion 
C  for the root-mean-square gradient, default $0.001$.
C
         BFGSMINT=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(GMAX)
         ENDIF

      ELSE IF (WORD.EQ.'BFGSSTEP') THEN
C
C  If starting from a transition state we just want to take one EF step using
C  BFGSTS before calling MYLBFGS (or something else). 
C
         BFGSSTEP=.TRUE.
         BFGSTST=.TRUE.
         IF (NITEMS.GT.1) CALL READF(PUSHOFF)

      ELSE IF (WORD .EQ. 'BFGSSTEPS') THEN
C
C  BFGSSTEPS n sets the number of BFGS optimisation steps to perform
C          per call to OPTIM                                    - default n=1     
C  If BFGSSTEPS is not specified then it is set to the same value as NSTEPS
C
        CALL READI(BFGSSTEPS)
        IF (NSTEPS.EQ.1) NSTEPS=BFGSSTEPS

      ELSE IF (WORD.EQ.'BFGSTS') THEN
C
C  Hybrid BFGS/eigenvector-following transition state search.
C
         BFGSTST=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(NEVS)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READI(NBFGSMAX1)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READI(NBFGSMAX2)
         ENDIF
         IF (NITEMS.GT.4) THEN
            CALL READF(CEIG)
         ENDIF
         IF (NITEMS.GT.5) THEN
            CALL READI(NEVL)
         ENDIF
         BFGSTST=.TRUE.
C
C  Tolerance for eigenvector overlap in BFGSTS where the number of tangent space
C  steps switches from small to large. 0.0001 was the traditional value (default). 
C
      ELSE IF (WORD.EQ.'BFGSTSTOL') THEN
         CALL READF(BFGSTSTOL)
C
C  Debug for basin-hopping interpolation
C
      ELSE IF (WORD.EQ.'BHDEBUG') THEN
         BHDEBUG=.TRUE.
C
C  Parameters for basin-hopping interpolation
C
      ELSE IF (WORD.EQ.'BHINTERP') THEN
         BHINTERPT=.TRUE.
         CALL READF(BHDISTTHRESH)
         CALL READF(BHMAXENERGY)
         CALL READI(BHSTEPS)
         CALL READF(BHCONV)
         CALL READF(BHTEMP)
         CALL READF(BHSTEPSIZE)
         CALL READF(BHACCREJ)
         CALL READF(BHK)
         CALL READF(BHSFRAC)
C
C  Additional parameter for basin-hopping interpolation.
C  Save the lowest energy minimum, rather than the lowest with the true PE plus spring energy.
C
C
      ELSE IF (WORD.EQ.'BHINTERPUSELOWEST') THEN
         BHINTERPUSELOWEST=.TRUE.
         IF (NITEMS.GT.1) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'CHECKENER') BHCHECKENERGYT=.TRUE.
         IF (NITEMS.GT.2) CALL READI(BHSTEPSMIN)
C
C  Binary LJ parameters for use with the LP or LS atom types.
C
      ELSE IF (WORD.EQ.'BINARY') THEN
         BINARY=.TRUE.
         CALL READI(NTYPEA)
         CALL READF(EPSAB)
         CALL READF(EPSBB)
         CALL READF(SIGAB)
         CALL READF(SIGBB)
C
C  Parameters for bisection runs
C
      ELSE IF (WORD.EQ.'BISECT') THEN
         BISECTT=.TRUE.
         CALL READF(BISECTMINDIST)
         CALL READF(BISECTMAXENERGY)
         CALL READI(BISECTSTEPS)
         CALL READI(BISECTMAXATTEMPTS)
         IF (NITEMS.GT.5) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ICINTERP') ICINTERPT=.TRUE.
C
C  Debug printing for BISECT runs.
C
      ELSE IF (WORD.EQ.'BISECTDEBUG') THEN
         BISECTDEBUG=.TRUE.
      ELSE IF (WORD.EQ.'BOND') THEN
         NUBONDS = NUBONDS + 1
         CALL READI(UBONDS(NUBONDS,1))
         CALL READI(UBONDS(NUBONDS,2))
C
C  General BLN model.
C
      ELSE IF (WORD.EQ.'BLN') THEN
         BLNT=.TRUE.
         CALL READF(RK_R)
         CALL READF(RK_THETA)
         ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     &            LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS))
         OPEN(UNIT=100,FILE='BLNsequence',STATUS='OLD')
         READ(100,*) DUMMYCH
         READ(100,*) LJREPBB, LJATTBB
         READ(100,*) LJREPLL, LJATTLL
         READ(100,*) LJREPNN, LJATTNN
         READ(100,*) DUMMYCH
         READ(100,*) DUMMYCH
         READ(100,*) HABLN, HBBLN, HCBLN, HDBLN
         READ(100,*) EABLN, EBBLN, ECBLN, EDBLN
         READ(100,*) TABLN, TBBLN, TCBLN, TDBLN
         DO J1=1,NATOMS-1
            READ(100,'(A1)',ADVANCE='NO') BEADLETTER(J1)
         ENDDO
         READ(100,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
         DO J1=1,NATOMS-3
            READ(100,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
         ENDDO
         CLOSE(100)
         PRINT '(A,I8,A)','BLN sequence of ',NATOMS,' beads read:'
         WRITE(*,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
         PRINT '(A)',' '
         PRINT '(A,I8,A)','BLN dihedral types:'
         WRITE(*,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
         PRINT '(A)',' '
         PRINT '(A,2F15.5)','B-B LJ coefficients: ',LJREPBB, LJATTBB
         PRINT '(A,2F15.5)','L-L LJ coefficients: ',LJREPLL, LJATTLL
         PRINT '(A,2F15.5)','N-N LJ coefficients: ',LJREPNN, LJATTNN
         PRINT '(A,4F15.5)','Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
         PRINT '(A,4F15.5)','Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
         PRINT '(A,4F15.5)','Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
         call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     &                       HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, NATOMS)         
C
C  Yimin Wang and Joel Bowman's water potential (2010)
C
      ELSE IF (WORD.EQ.'BOWMAN') THEN
         BOWMANT=.TRUE.
         CALL READI(BOWMANPES)
         CALL READA(BOWMANDIR)
C
C  BSMIN calculates a steepest-descent path using gradient only information
C  with convergence criterion GMAX for the RMS force and initial precision
C  EPS. The Bulirsch-Stoer algorithm is used.
C
      ELSE IF (WORD.EQ.'BSMIN') THEN
         BSMIN=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GMAX)
         IF (NITEMS.GT.2) CALL READF(EPS)

      ELSE IF (WORD.EQ.'BULK') THEN
         BULKT=.TRUE.
C
C  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CADPAC tells the program to read derivative information in 
C         CADPAC format.                                        - default FALSE
C
      ELSE IF (WORD.EQ.'CADPAC') THEN
         CADPAC=.TRUE.
         CALL READA(SYS)
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 10
            ENDIF
         ENDDO
10       IF (NITEMS.GT.2) THEN
            CALL READA(EDITIT)
         ELSE
            EDITIT='editit.' // SYS(1:LSYS)
         ENDIF
      ELSE IF (WORD.EQ.'CALCDIHE') THEN
         CALCDIHE=.TRUE.
C
C  If READPATH is specified with CALCRATES then the rates are calculated from the 
C  information in an existing path.info file without any stationary point searches.
C  A CONNECT or PATH run must be performed first unless READPATH is specified.
C
      ELSE IF (WORD.EQ.'CALCRATES') THEN
         CALCRATES=.TRUE.
         IF (NITEMS.GT.1) CALL READF(TEMPERATURE)
         IF (NITEMS.GT.2) CALL READF(HRED)
C
C  Double-ended connection keyword for ts candidates.
C
      ELSE IF (WORD == 'CANDIDATES') THEN
          CALL READA(CANDIDATES)
C
C  Virus capsid specification.
C
      ELSE IF (WORD.EQ.'CAPSID') THEN
         RIGIDBODY=.TRUE.
         ANGLEAXIS=.TRUE.
         HEIGHT=0.5D0
         CALL READF(CAPSRHO)
         CALL READF(CAPSEPS2)
         CALL READF(CAPSRAD)
         IF (NITEMS.GT.4) CALL READF(HEIGHT)
      ELSE IF (WORD.EQ.'CAPSID2') THEN
!         RIGIDBODY=.TRUE.
         ANGLEAXIS2=.TRUE.
         HEIGHT=0.5D0
         CALL READF(CAPSRHO)
         CALL READF(CAPSEPS2)
         CALL READF(CAPSRAD)
         IF (NITEMS.GT.4) CALL READF(HEIGHT)

C starting from a given residue, use cartesians for everything
      ELSE IF (WORD.EQ.'CARTRESSTART') THEN
         CALL READI(CARTRESSTART)
C
C  CASTEP tells the program to read derivative information in 
C         CASTEP format.                                        - default FALSE
C
      ELSE IF ((WORD.EQ.'CASTEP').OR.(WORD.EQ.'CASTEPC')) THEN
         CASTEP=.TRUE.
         IF (WORD.EQ.'CASTEP') DFTP=.TRUE.
         IF (NITEMS.GT.2) THEN
            CALL READA(CASTEPJOB)
            CALL READA(SYS)
            CASTEPJOB=TRIM(ADJUSTL(CASTEPJOB)) // ' ' // TRIM(ADJUSTL(SYS))
         ELSE
            WRITE(*,'(A)') 'keywords> ERROR - CASTEP job or system unspecified'
            CALL FLUSH(6,ISTAT)
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 22
            ENDIF
         ENDDO
22       CONTINUE
C
C charmm stuff (DAE)
C
      ELSE IF (WORD.EQ.'CHARMM') THEN
         CHRMMT=.TRUE.
         IF (.NOT.CISTRANS) THEN
            NOCISTRANS=.TRUE.
            CHECKOMEGAT=.TRUE.
         ENDIF
         CHECKCHIRALT=.TRUE. 

         IF ((PERMDIST.OR.LOCALPERMDIST).AND.(NPERMSIZE(1).EQ.NATOMS)) THEN
            PRINT '(A)','keyword> ERROR - PERMDIST or LOCALPERMDIST is specfied for CHARMM, but there is no perm.allow file present'
            STOP
         ENDIF
         CALL CHALLOCATE(NATOMS)
         ALLOCATE(ATMASS(NATOMS))
         IF (MACHINE) THEN 
              ! SAT: we will read in the coords ourselves and pass them to CHARMM {{{

              ! --- start ---
              ! read in the coords
              INQUIRE(IOLENGTH=J1) (Q(J),J=1,3*NATOMS)
              IF (FILTH2==0) THEN
                   OTEMP='points1.inp'
              ELSE
                   WRITE(OTEMP,*) FILTH2
                   OTEMP='points1.inp.'//TRIM(ADJUSTL(OTEMP))
              ENDIF
              OPEN(113,FILE=OTEMP,ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=J1)
              READ(113,REC=1) (Q(J),J=1,3*NATOMS)
              IF (MAXVAL(Q)==0.0D0) THEN
                  PRINT *, 'Zero coordinates - stop'
                  CALL FLUSH(6,ISTAT)
                  STOP
              ENDIF
              CLOSE(113)
              ! --- end ---
              ! SAT: line below was intended to replace the block of code above
              ! (marked); unfortunately, due to the miscompilation with pgi this
              ! does not work. The compiler does not really want to reuse the
              ! code. Sigh... 
              ! call ReadInpFile(Q)

              ! save them into CH. arrays and pass to CHARMM
              DO J1=1,NATOMS
                 CHX(J1)=Q(3*(J1-1)+1)
                 CHY(J1)=Q(3*(J1-1)+2)
                 CHZ(J1)=Q(3*(J1-1)+3)
              ENDDO
              CALL CHSETUP(CHX,CHY,CHZ,CHMASS,NATOM,TOPFILE,PARFILE)
C              CALL FILLICT(CHX,CHY,CHZ,DUMMY1,.TRUE.)
              CALL FILLICTABLE(Q)
              ! }}}
         ELSE 
              ! charmm will read the coords and will return them to OPTIM via CH. vecs {{{
              CHX(1)=13.13d13 ! this way we will tell CHARMM to save its coords into CH. arrays; otherwise it will
              CALL CHSETUP(CHX,CHY,CHZ,CHMASS,NATOM,TOPFILE,PARFILE)
              ! }}}
         ENDIF ! SAT
         CALL CHSETZSYMATMASS
         IF (FILTH.NE.0) THEN
            OPEN(UNIT=20,FILE='coords.read',STATUS='REPLACE')
            CLOSE(20)
         ENDIF
C        NATOMS=NATOM  ! should already know NATOMS from getparams
         IF (NATOM /= NATOMS) THEN
            WRITE(*,'(A)') 'No. of atoms in "input.crd" and file specified in CHARMM part of odata conflict'
            PRINT *, 'NATOM,NATOMS=',NATOM, NATOMS
            CALL FLUSH(6,ISTAT)
            STOP
         ENDIF
         CALL CHALLOCATE(NATOMS)
         CALL CHSETDIHE
!       csw34> If FREEZERES specified, call CHRESTOATOM to populate the
!       FROZEN array (from ocharmm.src)
         IF (FREEZERES) CALL CHRESTOATOM(FROZENRES,FROZEN)

         IF (CONNECTT) CALL SETSEED
C         IF (CALCDIHE) CALL READREF(NATOMS)
         DO J1=1,NATOMS 
            Q(3*(J1-1)+1)=CHX(J1)
            Q(3*(J1-1)+2)=CHY(J1)
            Q(3*(J1-1)+3)=CHZ(J1)
            ATMASS(J1) = CHMASS(J1)
C           PRINT *,'ATMASS',ATMASS(J1)
         ENDDO
         IF (TWISTDIHET) THEN
            CALL TWISTDIHE(Q,DMODE,DPERT)
         ENDIF
         IF (PERTDIHET) THEN
            CALL PERTDIHE(Q,CHPMIN,CHPMAX,CHNMIN,CHNMAX,ISEED)
         ENDIF
         IF (INTMINT) CALL GETNINT(NINTS)  ! DJW - this is OK because CHARMM is the last keyword!
C
      ELSE IF (WORD.EQ.'CHARMMTYPE') THEN
         IF (NITEMS.GT.1) THEN
            CALL READA(TOPFILE)
            TOPFILE=TRIM(ADJUSTL(TOPFILE))
         ENDIF   
         IF (NITEMS.GT.2) THEN
            CALL READA(PARFILE)
            PARFILE=TRIM(ADJUSTL(PARFILE))
         ELSE
            WRITE(*,*) 'keywords> TOPFILE and PARFILE have to be defined for CHARMMTYPE'
            STOP
         ENDIF 
         IF (TOPFILE(1:6).EQ."toph19") THEN
            CHARMMTYPE=2
         ELSEIF (TOPFILE(1:9).EQ."top_all22") THEN 
            CHARMMTYPE = 1
         ELSE
             WRITE(*,*) 'keywords> TOPFILE ', TRIM(ADJUSTL(TOPFILE)),' is not recognised by OPTIM'
             STOP
         ENDIF
         WRITE(*,'(A,I2)') 'CHARMMTYPE set to ',CHARMMTYPE
C
C If CHDEBUG is on, CHARMM related debug messages are printed
C
      ELSE IF (WORD.EQ.'CHDEBUG') THEN
         CHDEBUG=.TRUE.
         IF (NITEMS.GT.1) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'EDEBUG') EDEBUG=.TRUE.

C
C CHARMM related keyword, avoids inversion around C_alpha
C -- also implemented to AMBER (sf344) 
      ELSE IF (WORD.EQ.'CHECKCHIRALITY') THEN
         CHECKCHIRALT=.TRUE.
C
C  If CHECKINDEX is .TRUE. and the BFGSTS routine converges an attempt is
C  made to count the number of negative Hessian eigenvalues using projection,
C  orthogonalization and iteration. We also need the opportunity to change the
C  parameters NEVL and NEVS within BFGSTS if BFGSTS isn t true.
C  CHECKINDEX can also be used with BFGSMIN and should understand NOHESS too.
C
      ELSE IF (WORD.EQ.'CHECKINDEX') THEN
         CHECKINDEX=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(NEVS)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READF(CEIG)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READI(NEVL)
         ENDIF
C
C  If the index found by checkindex does not correspond to BFGSMIN or BFGSTS then
C  CHECKCONT causes a pushoff along the eigenvector correpsonding to the softest
C  undesired negative eigenvalue. 
C
      ELSE IF (WORD.EQ.'CHECKCONT') THEN
         CHECKCONT=.TRUE.
C
C  Check for internal minimum in constraint terms for INTBCONSTRAINT
C
      ELSE IF (WORD.EQ.'CONINT') THEN
         CHECKCONINT=.TRUE.
C
C  CHINTERPOLATE controls the interpolation for BHINTERP using CHARMM's primitive 
C  internal coordinates. The 1st argument has to be either BC or BI for the backbone
C  interpolation with Cartesians and Internals, respectively. The 2nd argument
C  has to be either SC or SI for the sidechain interpolation with Cartesians and 
C  Internals, respectively. If DNEB is given as 3rd argument, this interpolation scheme
C  will be used for DNEB. If CHINTERPOLATE is not defined in the odata file the default is
C  that DNEB and BHINTERP are done in Cartesians
C
      ELSE IF (WORD.EQ.'CHINTERPOLATE') THEN
         CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BI') CHBIT=.TRUE.
         CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SI') ICINTERPT=.TRUE.
         IF (NITEMS.GT.3) THEN
            CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'DNEB') CHICDNEB=.TRUE.
         ENDIF
C
C  If BHINTERPolation, and CHRIGID is set for the CHARMM potential, rigid body
C  translation and rotation is applied to the peptides/proteins if more
C  than one peptide/protein is prsent.
C
      ELSE IF (WORD.EQ.'CHRIGID') THEN
         CHRIGIDT=.TRUE.
         CALL READF(PTRANS)
         CALL READF(TRANSMAX)
         CALL READF(PROT)
         CALL READF(ROTMAX)
C
C CISTRANS is a CHARMM related keyword, which allows cis-trans isomerisation of the peptide bond .
C
      ELSE IF (WORD.EQ.'CISTRANS') THEN
         CISTRANS=.TRUE.
C
C  Sometimes have to modify the cold fusion limit when using high electric fields
C
      ELSE IF (WORD.EQ.'COLDFUSION') THEN
         IF (NITEMS.GT.1) call READF(COLDFUSIONLIMIT)
C
C  Connect initial minimum in odata to final minimum in file finish - maximum 
C  number of transiiton states=NCONNECT. Obsolete - use NEWCONNECT instead.
C
      ELSE IF (WORD.EQ.'CONNECT') THEN
         CONNECTT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNECT)
C
C  Constraint potential for interpolation between minima.
C
      ELSE IF (WORD.EQ.'CONPOT') THEN
         CONPOTT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(CPCONSTRAINTTOL)
         IF (NITEMS.GT.2) CALL READF(CPCONSTRAINTDEL)
         IF (NITEMS.GT.3) CALL READF(CPCONSTRAINTREP)
         IF (NITEMS.GT.4) CALL READF(CPCONSTRAINREPCUT)
         IF (NITEMS.GT.5) CALL READF(CPCONFRAC)
         IF (NITEMS.GT.6) CALL READI(CPCONSEP)
         IF (NITEMS.GT.7) CALL READI(CPREPSEP)
C
C jmc unres 
C Note also use some of the non-specific charmm keywords like INTMIN, NGUESS, TWISTTYPE etc...
C
      ELSE IF (WORD.EQ.'CONSEC') THEN
         CONSECT=.TRUE.
         DO J1=1,(NITEMS-1)/2
            CALL READI(STARTRES(J1))
            CALL READI(ENDRES(J1))
         END DO
         IF (NITEMS.GT.21) WRITE(*,'(A)') 'Too many sections requested - please adapt code!'
         NUMSEC=(NITEMS-1)/2
         PRINT *,'CONSEC ',(STARTRES(J1),J1=1,10),(ENDRES(J1),J1=1,10), NUMSEC
C
C  CONVERGE n m INDEX/NOINDEX sets the convergence criteria for the maximum 
C               unscaled step and RMS force                     - default n=0.0001, m=0.000001
C                                                           or m < 0.00001 .AND. n < m*100000  
C               If NOINDEX is set the Hessian index isn t checked - the default is
C               INDEX.
C
      ELSE IF (WORD .EQ. 'CONVERGE') THEN
        CALL READF(CONVU)
        IF (NITEMS.GT.2) THEN
           CALL READF(CONVR)
        ENDIF
        IF (NITEMS.GT.3) THEN
           CALL READU(WORD)
           IF (WORD.EQ.'NOINDEX') INDEXT=.FALSE.
        ENDIF
C
C  Probably prints the copyright info?
C
      ELSE IF (WORD == 'COPYRIGHT') THEN
          CALL COPYRIGHT
C     CP2K tells the program to read derivative information in 
C         CP2K format.                                        - default FALSE
C
      ELSE IF ((WORD.EQ.'CP2K').OR.(WORD.EQ.'CP2KC')) THEN
         CP2K=.TRUE.
         IF (WORD.EQ.'CP2K') DFTP=.TRUE.
         IF (NITEMS.GT.2) THEN
            CALL READA(CP2KJOB)
            CALL READA(SYS)
            CP2KJOB=TRIM(ADJUSTL(CP2KJOB)) // ' ' // TRIM(ADJUSTL(SYS))
         ELSE
            WRITE(*,'(A)') 'keywords> ERROR - no CP2K system specified'
            CALL FLUSH(6,ISTAT)
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 281
            ENDIF
         ENDDO
281      CONTINUE 
C
C  CPMD tells the program to read derivative information in 
C         CPMD format.                                        - default FALSE
C
      ELSE IF ((WORD.EQ.'CPMD').OR.(WORD.EQ.'CPMDC')) THEN
         CPMD=.TRUE.
         IF (WORD.EQ.'CPMDC') CPMDC=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') ' ERROR - no CPMD system specified'
            CALL FLUSH(6,ISTAT)
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 12
            ENDIF
         ENDDO
12       CONTINUE
         CALL SYSTEM(' grep -c DUMMY ' // SYS(1:LSYS) // ' > temp ')
         OPEN(UNIT=7,FILE='temp',STATUS='OLD')
         READ(7,*) J1
         IF (J1.NE.1) THEN
            WRITE(*,'(A)') 'ERROR, no dummy line in CPMD input file'
            CALL FLUSH(6,ISTAT)
            STOP
         ENDIF
C
C   Option to specify a different CPMD executible
C
         ELSE IF (WORD.EQ.'CPMD_COMMAND') THEN
            IF (NITEMS.GT.1) CALL READA(CPMD_COMMAND)
C
C  CUBIC: maintains cubic supercell for PV calculations
C
      ELSE IF (WORD.EQ.'CUBIC') THEN
         CUBIC=.TRUE.
C
C  For the growing string or evolving string double-ended
C  transition state search methods, use a cubic spline interpolation between
C  the image points.
C
      ELSE IF (WORD.EQ.'CUBSPL') THEN
         CUBSPLT = .TRUE.
C
C  DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
C
C  Add a decahedral field to the potential of magnitude FTD.
C
      ELSE IF (WORD.EQ.'D5H') THEN
         FIELDT=.TRUE.
         D5HT=.TRUE.
         CALL READF(FD5H)

      ELSE IF (WORD.EQ.'DB') THEN

         DBPT   = .TRUE.
         RBAAT  = .TRUE.
         CALL READF(DBEPSBB)
         CALL READF(DBEPSAB)
         CALL READF(DBSIGBB)
         CALL READF(DBSIGAB)
         CALL READF(DBPMU)
         IF (NITEMS > 6) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF

         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))

         NTSITES = NATOMS*NRBSITES/2

      ELSE IF (WORD.EQ.'DBTD') THEN

         DBPTDT = .TRUE.
         RBAAT  = .TRUE.
         CALL READF(DBEPSBB)
         CALL READF(DBEPSAB)
         CALL READF(DBSIGBB)
         CALL READF(DBSIGAB)
         CALL READF(DBPMU)
         IF (NITEMS > 6) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF 

         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))

         NTSITES = (NATOMS/2-1)*NRBSITES + 4

C
C  DCHECK  turns ON/OFF warnings about short interatomic distances
C                                                     - default ON
C
      ELSE IF (WORD.EQ.'DCHECK') THEN
         CALL READU(WW)
         IF (WW .EQ. 'ON' .OR. WW .EQ. ' ') THEN
            DCHECK=.TRUE.
         ELSE IF (WW .EQ. 'OFF') THEN
            DCHECK=.FALSE.
         ENDIF
C
C  DEBUG ON/OFF sets n=1 for EFSTEPS, VALUES, SUMMARY above     - default OFF 
C
      ELSE IF (WORD .EQ. 'DEBUG') THEN
         BHDEBUG=.TRUE.
        CALL READU(WW)
        IF (WW .EQ. 'ON' .OR. WW .EQ. ' ') THEN
          EFSTEPST=.TRUE.
          PGRAD=.TRUE.
          NGRADIENTS=1
          EFSTEPS=1
          NSUMMARY=1 
          NVALUES=1
          DEBUG=.TRUE.
          PRINTOPTIMIZETS=.TRUE.
          DUMPNEBXYZ=.TRUE.
          DUMPINTXYZ=.TRUE.
          DUMPNEBPTS=.TRUE.
          DUMPNEBEOS=.TRUE.
          DUMPINTEOS=.TRUE.
        ENDIF

      ELSE IF (WORD.EQ.'DESMAXAVGE') THEN
C maximum average energy before double ended search method can stop
         CALL READF(DESMAXAVGE)

      
C         ! maximum energy jump in one step
C         ! for an image in a double-ended search method
         CALL READF(DESMAXEJUMP)
C
C  Produces extra printing for the double-ended
C  transition state search method runs (DNEB, GS or ES).
C
      ELSE IF (WORD.EQ.'DESMDEBUG') THEN
         DESMDEBUG = .TRUE.

      ELSE IF (WORD.EQ.'DESMINT') THEN
         DESMINT = .TRUE.
         INTINTERPT = .FALSE. ! desmint and intinterp are mutually exclusive
         NATINT = .TRUE. ! must use natural internals for double ended search
C
C  DFTBT tells the program to call dftb for Tiffany s tight-binding.
C                                                  - default FALSE
      ELSE IF (WORD.EQ.'DFTB') THEN
         DFTBT=.TRUE.
C
C  Initial diagonal elements for LBFGS
C
       ELSE IF (WORD.EQ.'DGUESS') THEN
          CALL READF(DGUESS)
          IF (NITEMS.GT.2) CALL READF(XDGUESS)
          IF (NITEMS.GT.3) CALL READF(NEBDGUESS)
          IF (NITEMS.GT.4) CALL READF(INTDGUESS)
          IF (NITEMS.GT.5) CALL READF(GSDGUESS)
C
C  If DIJKSTRA is true then decide in newconnect uses Dijkstra;s algorithm in
C  deciding which connections to try next.
C  First argument on DIJKSTRA line controls the cost function. SAT
C
      ELSE IF (WORD.EQ.'DIJKSTRA') THEN
         IF (NITEMS.GT.1) THEN
            CALL READU(WW)
            IF (TRIM(ADJUSTL(WW))=='EXP') THEN
                 EXPCOSTFUNCTION = .TRUE.
            ELSEIF (TRIM(ADJUSTL(WW))=='INDEX') THEN
                 INDEXCOSTFUNCTION = .TRUE.
            ELSEIF (trim(adjustl(WW))=='INTERP') THEN
                 INTERPCOSTFUNCTION = .TRUE.
                 CALL READF(INTERPDIFF)
                 CALL READU(WW)
                 IF (TRIM(ADJUSTL(WW))=='EXP') THEN
                    EXPCOSTFUNCTION = .TRUE.
                 ELSE
C                   CALL READI(COSTFUNCTIONPOWER)
                    READ(WW,'(I20)') COSTFUNCTIONPOWER
                 ENDIF
            ELSE IF (WW(1:1) /= ' ') THEN
                 READ(WW,'(I20)') COSTFUNCTIONPOWER
            ENDIF
            IF (NITEMS.GT.2) THEN
               CALL READU(WW)
               IF (trim(adjustl(WW))=='INTDISTANCE') THEN
                  IF (.NOT.INTINTERPT) THEN
                     PRINT*, "INTDISTANCE doesn,t work without INTINTERP"
                     PRINT*, "specify the latter before DIJKSRA in odata"
                  ELSE
                     INTDISTANCET = .TRUE.
                  ENDIF
               ELSE
                  READ(WW,*) DIJKSTRADMAX
               ENDIF
            ENDIF
         ENDIF
C
C  DIJKSTRALOCAL specifies an adjustable factor used to multiply the
C  distances between minima found within one DNEB cycle. Decreasing
C  this metric will encourage attempts to complete the connection, which
C  might otherwise never be tried if shorter distances exist. We are
C  trying to correct for the imperfect nature of the distance criterion
C  used for the DIJKSTRA metric in choosing new connection pairs.
C
      ELSE IF (WORD.EQ.'DIJKSTRALOCAL') THEN
         CALL READF(DIJKSTRALOCAL)
C
C  Double well potential between first two atoms
C
      ELSE IF (WORD.EQ.'DOUBLE') THEN
         DOUBLET=.TRUE.
C
C  DNEB RMS threshold for switching to NEB
C
      ELSE IF (WORD.EQ.'DNEBSWITCH') THEN
         CALL READF(DNEBSWITCH)
C
C  Strings keyword.
C
      ELSE IF (WORD.EQ.'DQAGKEY') THEN
         CALL READI(DQAGKEY)
C
C  Obsolete: create a trajectory between the endpoints by increasing
C  a spring constant. Sounds like MD steering, but it doesn`t actually
C  work very well!
C
      ELSE IF (WORD.EQ.'DRAG') THEN
         DRAGT=.TRUE.
C
C  DUMPALLPATHS prints a summary of all min-sad-min triples produced by NEWCONNECT to
C  file path.info. For each stationary point the energy, point group order and symbol,
C  Hessian eigenvalues and coordinates are given. Hessian eigenvalues are computed
C  if not yet calculated, otherwise they are saved during the CONNECT process.
C
      ELSE IF (WORD.EQ.'DUMPALLPATHS') THEN
         DUMPALLPATHS=.TRUE.
         IF (FILTH.EQ.0) THEN
            WRITE(PINFOSTRING,'(A9)') 'path.info'
         ELSE
            WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
         ENDIF
         IF (MACHINE) THEN
              OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED')
         ELSE
             OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN')
         ENDIF
C
C  Creates a file in pathsample min.data format for the minimum found
C  following a minimisation. Useful for a DPS initial path run in
C  creating entries for the two endpoints.
C  Can also be used with BHINTERP alone to generate a list of entries 
C  for interpolated minima.
C
      ELSE IF (WORD.EQ.'DUMPDATA') THEN
         DUMPDATAT=.TRUE.
         IF (FILTH.EQ.0) THEN
            WRITE(PINFOSTRING,'(A13)') 'min.data.info'
         ELSE
            WRITE(PINFOSTRING,'(A)') 'min.data.info.'//TRIM(ADJUSTL(FILTHSTR))
         ENDIF
         IF (MACHINE) THEN
              OPEN(UNIT=881,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED')
         ELSE
              OPEN(UNIT=881,FILE=PINFOSTRING,STATUS='UNKNOWN')
         ENDIF
C
C  Explicit dump of interpolation EofS for intlbfgs. Should be set .TRUE. if DEBUG is set.
C
      ELSE IF (WORD == 'DUMPINTEOS') THEN
          DUMPINTEOS=.TRUE.
          IF (NITEMS>1) CALL READI(DUMPINTEOSFREQ)
C
C  Explicit dump of EofS.neb for DNEB. Should be set .TRUE. if DEBUG is set.
C
      ELSE IF (WORD == 'DUMPNEBEOS') THEN
          DUMPNEBEOS=.TRUE.
          IF (NITEMS>1) CALL READI(DUMPNEBEOSFREQ)
C
C  Explicit dump of something for DNEB. Should be set .TRUE. if DEBUG is set.
C
      ELSE IF (WORD == 'DUMPNEBPTS') THEN
          DUMPNEBPTS=.TRUE.
          IF (NITEMS>1) CALL READI(DUMPNEBPTSFREQ)
C
C  Explicit dump of image coordinates in xyz format for intlbfgs. Should
C  be set .TRUE. if DEBUG is set.
C
      ELSE IF (WORD == 'DUMPINTXYZ') THEN
          DUMPINTXYZ=.TRUE.
          IF (NITEMS>1) CALL READI(DUMPINTXYZFREQ)
C
C  Explicit dump of image coordinates in xyz format for DNEB. Should
C  be set .TRUE. if DEBUG is set.
C
      ELSE IF (WORD == 'DUMPNEBXYZ') THEN
          DUMPNEBXYZ=.TRUE.
          IF (NITEMS>1) CALL READI(DUMPNEBXYZFREQ)
C
C  DUMPPATH prints a summary of a min-sad-min-...-min path produced by CONNECT to
C  file path.info. For each stationary point the energy, point group order and symbol, 
C  Hessian eigenvalues and coordinates are given. Hessian eigenvalues are computed
C  if not yet calculated, otherwise they are saved during the CONNECT process.
C
      ELSE IF (WORD.EQ.'DUMPPATH') THEN
         DUMPPATH=.TRUE.
         IF (FILTH.EQ.0) THEN
            WRITE(PINFOSTRING,'(A9)') 'path.info'
         ELSE
            WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
         ENDIF
         IF (MACHINE) THEN
              OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED')
         ELSE
              OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN')
         ENDIF
C
C  If DUMPSP is true then OPTIM will dump minima and ts data in the pathsample format
C
      ELSE IF (WORD.EQ.'DUMPSP') THEN
         DUMPSP=.TRUE.
C
C  DUMPVECTOR switches on dumping of eigenvectors to file 
C              vectors.dump                                     - default OFF
C  ALLSTEPS dumps the vector(s) at each step. ALLVECTORS dumps all the vectors.
C  The defaults are for only the vector corresponding to the softest non-zero
C  eigenvalue to be dumped for the last step.
C
      ELSE IF (WORD .EQ. 'DUMPVECTOR') THEN
        DUMPV=.TRUE.
        IF (NITEMS.GT.1) THEN
           CALL READU(WORD)
           IF (WORD.EQ.'ALLSTEPS') ALLSTEPS=.TRUE.
           IF (WORD.EQ.'ALLVECTORS') ALLVECTORS=.TRUE.
           IF (WORD.EQ.'MWVECTORS') MWVECTORS=.TRUE.
        ENDIF
        IF (NITEMS.GT.2) THEN
           CALL READU(WORD)
           IF (WORD.EQ.'ALLSTEPS') ALLSTEPS=.TRUE.
           IF (WORD.EQ.'ALLVECTORS') ALLVECTORS=.TRUE.
           IF (WORD.EQ.'MWVECTORS') MWVECTORS=.TRUE.
        ENDIF
        IF (NITEMS.GT.3) THEN
           CALL READU(WORD)
           IF (WORD.EQ.'ALLSTEPS') ALLSTEPS=.TRUE.
           IF (WORD.EQ.'ALLVECTORS') ALLVECTORS=.TRUE.
           IF (WORD.EQ.'MWVECTORS') MWVECTORS=.TRUE.
        ENDIF
C
C  EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C
C  EDIFFTOL specifies the maximum energy difference between permutational isomers in connect.
C
      ELSE IF (WORD.EQ.'EDIFFTOL') THEN
         CALL READF(EDIFFTOL)
C
C  Specify an electric field in the z-direction, units are V/A
C  So far only implemented for use with TIPnP potentials
C
      ELSE IF (WORD.EQ.'EFIELD') THEN
         IF (NITEMS.GT.1) CALL READF(EFIELD)
C
C  EFSTEPS n print the unscaled steps calculated for each mode
C          every n cycles                                       - default OFF
C
      ELSE IF (WORD .EQ. 'EFSTEPS') THEN
        EFSTEPST=.TRUE.
        CALL READI(EFSTEPS)
C
C  Calculate analytical Hessian and normal mode frequencies at end of run.
C  ENDHESS is only intended for use in single geometry optimisations, and
C  should not be needed for CONNECT or PATH runs if DUMPPATH is specified.
C  If the argument NENDHESS is omitted then all the eigenvalues are
C  calculated - otherwise just the lowest NENDHESS.
C
      ELSE IF (WORD.EQ.'ENDHESS') THEN
         ENDHESS=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NENDHESS)
C
C  Calculate numerical Hessian and normal mode frequencies at end of run.
C  Required if DUMPPATH or ENDHESS is specified for an UNRES run,
C  in which case it`s an internal coordinate Hessian, or for other potentials
C  that lack analytic second derivatives.
C 
      ELSE IF (WORD.EQ.'ENDNUMHESS') THEN
         ENDNUMHESS=.TRUE.
C
C
C
      ELSE IF (WORD == 'ERROREXIT') THEN
         PRINT *, 'ERROR EXIT'
         CALL FLUSH(6,ISTAT)
         STOP
C
C  Cutoff below which Hessian eigenvalues are considered to be zero.
C
      ELSE IF (WORD.EQ.'EVCUT') THEN
         CALL READF(EVCUT)
C
C  Specify evolving strings.
C
      ELSE IF (WORD.EQ.'EVOLVESTRING') THEN
         EVOLVESTRINGT = .TRUE.      
C
C  sf344> extra repulsive LJ site for PY ellipsoids
C
      ELSE IF (WORD.EQ.'EXTRALJSITE') THEN
         LJSITE=.TRUE.
         CALL READF(PEPSILON1(1))
         CALL READF(PSCALEFAC1(1))
          MAXINTERACTIONS=1
         IF(NITEMS.GT.3) THEN
          CALL READF(PSCALEFAC2(1))
          WRITE(MYUNIT,'(A,3F8.3)') ' keyword> primary and secondary apex sites will be used, epsilon and heights: ', 
     &                              PEPSILON1(1), PSCALEFAC1(1), PSCALEFAC2(1)
          IF(.NOT.LJSITEATTR) THEN
                MAXINTERACTIONS=3
          ELSE
                MAXINTERACTIONS=4
          END IF
         ELSE
          WRITE(MYUNIT,'(A,2F8.3)') ' keyword> primary apex sites will be used, epsilon and height: ', PEPSILON1(1), PSCALEFAC1(1)
         END IF
         IF(NITEMS.GT.4) THEN           ! binary ellipsoidal clusters will be set up only for two apex sites, not one
           BLJSITE=.TRUE.               ! we also won't use the sigma parameter from now on, epsilon is enough for repulsive sites
           CALL READF(PEPSILON1(2))
           CALL READF(PSCALEFAC1(2))
           CALL READF(PSCALEFAC2(2))
           CALL READF(PEPSILON1(3))     ! this is epsilon for the interaction between A and B type ellipsoids
           MAXINTERACTIONS=3 ! attractive secondary apex sites not incorporated for binary systems
           WRITE(MYUNIT,'(A,3F8.3)') ' keyword> binary system with primary and secondary apex sites, ' //  
     &  'epsilon and heights for 2nd type particle: ', PEPSILON1(2), PSCALEFAC1(2), PSCALEFAC2(2)

         END IF
      ELSE IF (WORD.EQ.'EXTRALJSITEATTR') THEN
         LJSITE=.TRUE.
         LJSITEATTR=.TRUE.
         CALL READF(PSIGMAATTR(1))
         CALL READF(PEPSILONATTR(1))
         CALL READF(PSIGMAATTR(2))
         CALL READF(PEPSILONATTR(2))

         WRITE(MYUNIT,'(A,4F8.3)') 'keyword> primary and secondary apex sites '//
     &                             'with normal LJ attraction, sigmas and epsilons: ', 
     &                             PSIGMAATTR(1), PEPSILONATTR(1), PSIGMAATTR(2), PEPSILONATTR(2)
         MAXINTERACTIONS=4
      ELSE IF (WORD.EQ.'LJSITECOORDS') THEN
           LJSITECOORDST=.TRUE.
           CALL READF(LJSITECOORDS(1))
           CALL READF(LJSITECOORDS(2))
           CALL READF(LJSITECOORDS(3))
      ELSE IF (WORD.EQ.'PYBINARY') THEN
         PYBINARYT=.TRUE.
         ANGLEAXIS2=.TRUE.
         RBAAT=.TRUE.
         ELLIPSOIDT=.TRUE.
         RADIFT=.TRUE.
         NRBSITES = 1
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL READI(PYBINARYTYPE1)
         CALL READF(PYA11(1))
         CALL READF(PYA11(2))
         CALL READF(PYA11(3))
         CALL READF(PYA21(1))
         CALL READF(PYA21(2))
         CALL READF(PYA21(3))
         CALL READF(PYA12(1))
         CALL READF(PYA12(2))
         CALL READF(PYA12(3))
         CALL READF(PYA22(1))
         CALL READF(PYA22(2))
         CALL READF(PYA22(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)
         IF(NITEMS.GT.16) THEN
            CALL READF(PCUTOFF)
            PARAMONOVCUTOFF=.TRUE.
            PCUTOFF=PCUTOFF*PYSIGNOT
            write (MYUNIT,*) "PY Potential. PCutoff ON:",PCUTOFF
         END IF
         IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMS/2,3))
         IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMS/2,3))
         DO J1=1,NATOMS/2
          IF(J1<=PYBINARYTYPE1) THEN
           PYA1bin(J1,:)=PYA11(:)
           PYA2bin(J1,:)=PYA21(:)
          ELSE
           PYA1bin(J1,:)=PYA12(:)
           PYA2bin(J1,:)=PYA22(:)
          END IF
         END DO
C
C  Obsolete. Allows for extra steps in LBFGS minimisations for CHARMM.
C
      ELSE IF (WORD.EQ.'EXTRASTEPS') THEN
         CALL READF(EXTRASTEPS)
C
C  FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
C
C
C  Distance dependent dielectric for Paul Mortenson`s amber
C
      ELSE IF (WORD.EQ.'FAKEWATER') THEN
         FAKEWATER=.TRUE.
         WRITE (*,'(A)') ' SETTINGS Distance dependent dielectric will be used'
C
C  Integer variable to distinguish output files from parallel maiden jobs
C
      ELSE IF (WORD.EQ.'FILTH') THEN
         IF (FILTH.EQ.0) THEN
            CALL READI(FILTH)
         ELSE
            WRITE(*,'(A)') 'WARNING **** FILTH keyword in odata was overridden by command line argument'
         ENDIF
C
C  Specifies that FIXIMAGE should be set permanently after step
C  FIXAFTER. This effectively freezes the interacting images in different supercells
C  for calculations with periodic boundary conditions.
C
      ELSE IF (WORD.EQ.'FIXAFTER') THEN
         CALL READI(FIXAFTER)
C
C  Strings keyword.
C
      ELSE IF (WORD.EQ.'FIXATMS') THEN
         FIXATMS = .TRUE.
C
C  Fix uphill direction until force changes sign.
C  T12FAC is the fraction of the first collision time to be used in HSMOVE
C
      ELSE IF (WORD.EQ.'FIXD') THEN
         FIXD=.TRUE.
         IF (NITEMS.GT.1) CALL READF(T12FAC)
         NMOVE=1
         IF (NITEMS.GT.2) CALL READF(DTHRESH)
C
C  FRACTIONAL: constant pressure calculation using fractional coordinates
C
      ELSE IF (WORD.EQ.'FRACTIONAL') THEN
         FRACTIONAL=.TRUE.
C
C  Frozen atoms.
C
      ELSE IF (WORD.EQ.'FREEZE') THEN
         FREEZE=.TRUE.
         DO J1=1,NITEMS-1
            NFREEZE=NFREEZE+1
            CALL READI(NDUM)
            FROZEN(NDUM)=.TRUE.
         ENDDO

C csw34
C Frozen residues (to be converted to frozen atoms)
C
      ELSE IF (WORD.EQ.'FREEZERES') THEN
         FREEZE=.TRUE.
         FREEZERES=.TRUE.
C The FROZENRES array is then filled with the residue number from the
C data file
         DO J1=1,NITEMS-1
            CALL READI(NDUM)
            FROZENRES(NDUM)=.TRUE.
         ENDDO
C Finally, the frozen residue numbers are converted into frozen atom
C numbers. This is also forcefield dependant and must be done when we
C know which forcefield to use (i.e. in the CHARMM block above)
!
! csw34> FREEZEGROUP centreatom radius
! FREEZEs all atoms within radius angstroms of centreatom (labelled by index)
!
      ELSE IF (WORD.EQ.'FREEZEGROUP') THEN
         FREEZE=.TRUE.
         FREEZEGROUPT=.TRUE.
         CALL READI(GROUPCENTRE)
         CALL READF(GROUPRADIUS)
         IF(NITEMS.GT.3) CALL READA(FREEZEGROUPTYPE)
     
      ELSE IF (WORD.EQ.'DUMPMODE') THEN
                IF(DUMPV.AND.MWVECTORS) THEN
                        DO J1=1,NITEMS-1
                                CALL READI(NDUM)
                                DUMPMODEN(NDUM)=.TRUE.
                        ENDDO
                ENDIF

         
         IF ((PERMDIST.OR.LOCALPERMDIST).OR.PERMDISTINIT) THEN
            NDUMMY=0
            DO J1=1,NPERMGROUP
               DO J2=1,NPERMSIZE(J1)
                  IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                     PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                     STOP
                  ENDIF
               ENDDO
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
         ENDIF
C
C  Strings keyword.
C
      ELSE IF (WORD.EQ.'FREEZENODES') THEN
         FREEZENODEST=.TRUE.
         CALL READF(FREEZETOL)         
C
C  GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
C
C
C  GAMESS-UK tells the program to read derivative information in 
C         GAMESS-UK format.                                        - default FALSE
      ELSE IF (WORD.EQ.'GAMESS-UK') THEN
         GAMESSUK=.TRUE.
         CALL READA(SYS)
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 112
            ENDIF
         ENDDO
112       IF (NITEMS.GT.2) THEN
            CALL READA(EDITIT)
         ELSE
            EDITIT='editit.' // SYS(1:LSYS)
         ENDIF
C
C  GAMESS-US tells the program to read derivative information in 
C         GAMESS-US format.                                        - default FALSE
      ELSE IF (WORD.EQ.'GAMESS-US') THEN
         GAMESSUS=.TRUE.
         CALL READA(SYS)
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 111
            ENDIF
         ENDDO
111       IF (NITEMS.GT.2) THEN
            CALL READA(EDITIT)
         ELSE
            EDITIT='editit.' // SYS(1:LSYS)
         ENDIF
C
C  GAUSSIAN tells the program to read derivative information in 
C           Gaussian92 format.                                  - default FALSE
      ELSE IF (WORD.EQ.'GAUSSIAN') THEN
         GAUSSIAN=.TRUE.

!     DC430 >

      ELSE IF (WORD .EQ. 'GB') THEN
         GBT   = .TRUE.
         RBAAT = .TRUE.
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(GBSIGNOT)
         CALL READF(GBEPSNOT)

         GBCHI    = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
         GBCHIPRM = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)

      ELSE IF (WORD .EQ. 'GBD') THEN

         GBDT  = .TRUE.
         RBAAT = .TRUE.
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(GBSIGNOT)
         CALL READF(GBEPSNOT)

         GBCHI    = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
         GBCHIPRM = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)

!     -----------------------------
C  GDIIS x y z  x=cutoff on previous RMS force below which GDIIS
C               may be applied, y=NDIIA the dimension of the DIIS
C               problem to solve, z=NINTV the interval between
C               GDIIS steps                                     - default OFF       
C
      ELSE IF (WORD .EQ. 'GDIIS') THEN
        PRINT '(A)','keyword> GDIIS keyword not available'
        STOP
!       IF (NITEMS.LT.4) THEN
!          DTEST=.FALSE.
!          PRINT*,'Error in GDIIS input - insufficient items'
!       ELSE
!          DTEST=.TRUE.
!          CALL READF(PCUT)
!          CALL READI(NDIIA)
!          CALL READI(NINTV)
!       ENDIF
!       IF (NDIIA.GT.NDIIS) THEN
!          WRITE(*,'(A,I6)') ' NDIIA too large=',NDIIA
!          STOP
!       ENDIF
C
C  GEOMDIFFTOL specifies the maximum displacement between identical permutational isomers in connect.
C
      ELSE IF (WORD.EQ.'GEOMDIFFTOL') THEN
         CALL READF(GEOMDIFFTOL)
C
C  Paul Whitford Go model
C
      ELSE IF (WORD.EQ.'GOT') THEN
         GOT=.TRUE.
C
C  GRADIENT n prints the gradients along the Hessian eigendirections
C             every n cycles                                    - default OFF
C
      ELSE IF (WORD .EQ. 'GRAD4') THEN ! 4-point gradient in finite_differences.f90
         GRAD4T=.TRUE.
         print *, 'use 4-point gradient'
      ELSE IF (WORD .EQ. 'GRADIENTS') THEN
        PGRAD=.TRUE.
        CALL READI(NGRADIENTS)
C
C  GRADSQ specifies optimisation of the modulus gradient. This is a really bad idea!
C
      ELSE IF (WORD.EQ.'GRADSQ') THEN
         GRADSQ=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GSTHRESH)
         IF (NITEMS.GT.2) CALL READI(NSPECIAL)
         IF (NITEMS.GT.3) CALL READI(NALLOW)
C
C  Approximation to use for the gradient in DNEB routine NEB/grad.f90 - default is "dneb"
C
      ELSE IF (WORD == 'GRADTYPE') THEN
          CALL READA(GRADTYPE)
C
C  Attempt to interpolate between endpoints using a great circle. Not a huge success.
C
      ELSE IF (WORD.EQ.'GREATCIRCLE') THEN
         GREATCIRCLET=.TRUE.
         CALL READF(GCMXSTP)
         CALL READI(GCIMAGE)
         CALL READI(GCSTEPS)
         CALL READF(GCCONV)
C
C EF: growing strings
C number of images and iterations for first iteration; 
C reparametrization tolerance, growth tolerance, convergence tolerance
C maximum LBFGS step; LBFGS memory
C
      ELSE IF (WORD.EQ.'GROWSTRING') THEN
         GROWSTRINGT = .TRUE.
         FCD = .TRUE.         
         IF (NITEMS.GT.1) CALL READI(nnNIMAGE)
         IF (NITEMS.GT.2) CALL READF(GSITERDENSITY)
         IF (NITEMS.GT.3) CALL READF(REPARAMTOL)
         IF (NITEMS.GT.4) CALL READF(GSGROWTOL)
         IF (NITEMS.GT.5)CALL READF(GSCONV)
         IF (NITEMS.GT.6) CALL READF(GSMXSTP)
C
C  Set the maximum total iteration density for the
C  growing string method. This specifies the maximum evolution iterations allowed per
C  total image number, including the iterations while the string is still 
C  growing. If {\it itd\/} is less than 0, then this parameter is turned off
C  and there is no limit on the total iterations (this is the default).
C
      ELSE IF (WORD.EQ.'GSMAXTOTITD') THEN
         CALL READF(GSMAXTOTITD)
C
C  Try to guess an interpolated path for sequential coordinate changes.
C  GSTEPS is the number of step to be tried for each sequential coordinate.
C  MAXGCYCLES is the number of sweeps through pairwise exchanges
C  GTHRESHOLD is the coordinate change above which sequential changes are considered.
C  MAXINTE is the convergence criterion for the maximum allowed interpolated energy.
C
      ELSE IF (WORD.EQ.'GUESSPATH') THEN
         GUESSPATHT=.TRUE.
         CALL READI(GSTEPS)
         CALL READI(MAXGCYCLES)
         CALL READF(GTHRESHOLD)
C        CALL READF(MAXINTE)
C
C  Use dihedral twisting in place of DNEB for
C  transition state guesses with CONNECT for CHARMM and UNRES.
C
      ELSE IF (WORD.EQ.'GUESSTS') THEN
         GUESSTST=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GUESSTHRESH)
      ELSE IF (WORD.EQ.'GUPTA') THEN
         CALL READI(GUPTATYPE)
C
C  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
C
C  For the growing string or evolving string double-ended
C  transition state search method, use the method described in the appendix of
C  Peters et al\cite{PetersHBC04} to calculate the Newton-Raphston search
C  direction. Namely, the Hessian is approximated based on changes in the
C  gradient, and the tangential component of $-\mathbf{Hf^\perp}$ is projected
C  out. By default, the Hessian used is actually an approximation to the
C  derivative matrix of $f^\perp$ rather than the gradient.
C
      ELSE IF (WORD.EQ.'HESSGRAD') THEN 
         HESSGRAD = .TRUE.
C
C  HIGHESTIMAGE - only use the highest non-endpoint image in a double-ended
C  MECCANO-type run.
C
      ELSE IF (WORD.EQ.'HIGHESTIMAGE') THEN
         HIGHESTIMAGE=.TRUE.
C
C  HUPDATE specifies that a Hessian updating procedure should be used.
C
      ELSE IF (WORD .EQ. 'HUPDATE') THEN
        HUPDATE=.TRUE.
        NHUP=0
        IF (NITEMS.GT.1) THEN
           CALL READI(NSTHUP)
        ENDIF
        IF (NITEMS.GT.2) THEN
           CALL READI(INTHUP)
        ENDIF
        IF (NITEMS.GT.3) THEN
           CALL READF(PHIG)
        ENDIF
C
C  Hybrid BFGS/eigenvector-following minimisation.
C
      ELSE IF (WORD.EQ.'HYBRIDMIN') THEN
         HYBRIDMINT=.TRUE.
         CALL READI(HMNEVS)      ! maximum steps to converge smallest eigenvalue in Rayleigh-Ritz
         CALL READI(HMNBFGSMAX1) ! maximum tangent space LBFGS steps if eigenvalue unconverged
         CALL READI(HMNBFGSMAX2) ! maximum tangent space LBFGS steps if eigenvalue converged
         CALL READF(HMCEIG)      ! convegence criterion for eigenvalue
         CALL READF(HMMXSTP)     ! maximum step size for EF steps
         CALL READI(HMNSTEPS)    ! maximum number of hybrid minimisation steps
         CALL READF(HMEVMAX)     ! If the lowest eigenvalue goes above HMEVMAX then exit
         CALL READA(HMMETHOD)    ! Choose between EF and Page-McIver steepest-descent steps
         IF (NITEMS.GT.9) THEN
            CALL READI(HMNEVL)   ! maximum steps for iterative calculation of largest eigenvalue if applicable
         ENDIF
C
C  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
C
C  Add an icosahedral field to the potential of magnitude FIH.
C
      ELSE IF (WORD.EQ.'IH') THEN
         FIELDT=.TRUE.
         IHT=.TRUE.
         CALL READF(FIH)
C
C
C
      ELSE IF (WORD.EQ.'IMSEP') THEN
         IF (NITEMS.GT.1) CALL READF(IMSEPMIN)
         IF (NITEMS.GT.2) CALL READF(IMSEPMAX)
C
C  Search for a saddle of index INDEX if
C  SEARCH 2 is specified. See also KEEPINDEX. Also works with BFGSTS
C  up to a maximum of index 50, but NOIT must be set and a Hessian is needed.  
C
      ELSE IF (WORD.EQ.'INDEX') THEN
         CALL READI(HINDEX)
C
C  Use constraint potential for initial interpolation in each cycle.
C
      ELSE IF (WORD.EQ.'INTCONSTRAINT') THEN
         INTCONSTRAINTT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(INTCONSTRAINTTOL)
         IF (NITEMS.GT.2) CALL READF(INTCONSTRAINTDEL)
         IF (NITEMS.GT.3) CALL READF(INTCONSTRAINTREP)
         IF (NITEMS.GT.4) CALL READF(INTCONSTRAINREPCUT)
         IF (NITEMS.GT.5) CALL READF(INTCONFRAC)
         IF (NITEMS.GT.6) CALL READI(INTCONSEP)
         IF (NITEMS.GT.7) CALL READI(INTREPSEP)
         IF (NITEMS.GT.8) CALL READI(INTSTEPS1)
         IF (NITEMS.GT.9) CALL READI(INTCONSTEPS)
         IF (NITEMS.GT.10) CALL READI(INTRELSTEPS)
         IF (NITEMS.GT.11) CALL READF(MAXCONE)
         IF (NITEMS.GT.12) CALL READF(INTRMSTOL)
         IF (NITEMS.GT.13) CALL READI(INTIMAGE)
         IF (NITEMS.GT.14) CALL READI(MAXINTIMAGE)
!
! Use the quasi-continuous metric for connection attempts, instead of distance.
!
         INTERPCOSTFUNCTION=.TRUE.
C
C  Use interpolation potential for LJ.
C
      ELSE IF (WORD.EQ.'INTLJ') THEN
         INTLJT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(INTLJSTEPS)
         IF (NITEMS.GT.2) CALL READF(INTLJDEL)
         IF (NITEMS.GT.3) CALL READF(INTLJTOL)
         IF (NITEMS.GT.4) CALL READI(INTIMAGE)
         IF (NITEMS.GT.5) CALL READF(INTLJEPS)
!
! Use the quasi-continuous metric for connection attempts, instead of distance.
!
         INTERPCOSTFUNCTION=.TRUE.
C
C  Epsilon value in internal coordinate optimisation.
C
      ELSE IF (WORD == 'INTEPSILON') THEN
         IF (NITEMS.GT.1) CALL READF(INTEPSILON)

C     back transformation cutoff for interpolation in internals         
      ELSE IF (WORD.EQ.'INTERPBACKTCUT') THEN
         CALL READF(INTERPBACKTCUT)

C must set INTINTERP as well to use INTERPCHOICE
C use internals or cartesian interpolation depending on which gives
C the lower max energy
      ELSE IF (WORD.EQ.'INTERPCHOICE') THEN
         INTERPCHOICE = .TRUE.

      ELSE IF (WORD.EQ.'INTINTERP') THEN
         INTINTERPT = .TRUE. ! interpolate with internals
         NATINT = .TRUE. ! if interpolating, assume natural internal coords
         DESMINT = .FALSE. ! intinterp and desmint are mutually exclusive
         IF (NITEMS.GT.1) CALL READI(NINTIM)
        
C when interpolating with internals, keep actual interpolation points.
C don't distribute images between them to make them equidistant in cartesians
      ELSE IF (WORD.EQ.'INTERPSIMPLE') THEN
         INTERPSIMPLE = .TRUE.
C
C  Internal coordinate minimisation - do not use.
C   IMINCUT is the RMSG below which we take steps in internal coordinates
C
C
      ELSE IF (WORD.EQ.'INTMIN') THEN
         INTMINT=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(IMINCUT)
         ENDIF
         
C align permutations of starting structures to match up internals
      ELSE IF (WORD.EQ.'INTMINPERM') THEN
         INTMINPERMT = .TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READA(WORD2)
            IF (WORD2.EQ."GLYCART") THEN
               GLYCART = .TRUE.
            ELSE 
               PRINT*, "keyword error intminperm"
            ENDIF
         ENDIF

      ELSE IF (WORD.EQ.'INTPARFILE') THEN
         USEPARFILE = .TRUE.
         CALL READA(INTPARFILE) ! file with internals parameters
C
C  JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
C
C
C  KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
C
C
C  KEEPINDEX: specifies that INDEX is set with 
C  the number of negative Hessian eigenvalues at the initial point.
C
      ELSE IF (WORD.EQ.'KEEPINDEX') THEN
         KEEPINDEX=.TRUE.

!       csw34> Specify kT in wavenumbers, below which a normal mode is 
!              determined to be thermally accessible. KTWN defaults to
!              room temperature (207.11cm-1). This is used by the 
!              CHARMMDUMPMODES subroutine
      ELSE IF (WORD.EQ.'KTWN') THEN
         KTWNT=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(KTWN)
         ENDIF

C
C  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
C  Use Lanczos to diagonalize the Hamiltonian. Defaults for the three
C  associated parameters are ACCLAN=1.0D-8 SHIFTLAN=1.0D-2 CUTLAN=-1.0D0.
C
      ELSE IF (WORD.EQ.'LANCZOS') THEN
         LANCZOST=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(ACCLAN)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READF(SHIFTLAN)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(CUTLAN)
         ENDIF
      ELSE IF (WORD.EQ.'LWOTP') THEN
         LWOTPT = .TRUE.
         RBAAT  = .TRUE.    
         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))
         NTSITES = NATOMS*NRBSITES/2

         ELSE IF (WORD.EQ.'LOWESTFRQ') THEN
            LOWESTFRQT=.TRUE.
C
C  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
      ELSE IF (WORD.EQ.'MACHINE') THEN
          MACHINE=.TRUE.
C
C  MASS ON/OFF takes steps with a fictitious "kinetic" metric   - default OFF       
C
      ELSE IF (WORD .EQ. 'MASS') THEN
        CALL READU(WW)
        IF (WW .EQ. 'ON' .OR. WW .EQ. ' ') THEN
          MASST=.TRUE.
        ELSE IF (WW .EQ. 'OFF') THEN
          MASST=.FALSE.
        ENDIF
C
C  Maximum value for the smaller barrier height that is allowed to constitute a connection during the
C  Dijkstra connection procedure.
C  MAXMAXBARRIER specifies a maximum for the maximum barrier. 
C  MAXBARRIER requires both sides to be greater than MAXBARRIER to discard.
C
      ELSE IF (WORD.EQ.'MAXBARRIER') THEN
         CALL READF(MAXBARRIER)
      ELSE IF (WORD.EQ.'MAXMAXBARRIER') THEN
         CALL READF(MAXMAXBARRIER)
C
C  MAXBFGS x1 x2 x3 x4\/}: {\it x\/} specifies the maximum allowed step length in LBFGS
C  minimisations, {\it x1\/} for  normal minimisations, {\it x2\/} for Rayleigh-Ritz ratio
C  minimisation, {\it x3\/} for putting structures in closest coincidence with
C  {\bf mind} (NO LONGER USED!!), and {\it x4\/} for NEB minimisations. Default values all 0.2.
C
      ELSE IF (WORD.EQ.'MAXBFGS') THEN
         CALL READF(MAXBFGS)
         IF (NITEMS.GT.2) CALL READF(MAXXBFGS)
         IF (NITEMS.GT.3) CALL READF(MAXMBFGS)
         IF (NITEMS.GT.4) CALL READF(MAXNEBBFGS)
         IF (NITEMS.GT.5) CALL READF(MAXINTBFGS)
C
C  The maximum number of constraints to use in the constrained potential.
C  The deafult is 3.
C
      ELSE IF (WORD.EQ.'MAXCON') THEN
         CALL READI(MAXCONUSE)
C
C  The maximum energy increase above which mylbfgs will reject a proposed step.
C
      ELSE IF (WORD.EQ.'MAXERISE') THEN
         CALL READF(MAXERISE)
         IF (NITEMS.GT.1) CALL READF(XMAXERISE)
C
C  Maximum number of failures allowed in a minimisation before giving up.
C
      ELSE IF (WORD.EQ.'MAXFAIL') THEN
         IF (NITEMS.GT.1) CALL READI(NFAILMAX)
C
C  For the growing string double-ended connection
C  method, specify a maximum number of steps allowed before another image is
C  added to the growing string. Default is 1000.
C
      ELSE IF (WORD.EQ.'MAXGROWSTEPS') THEN
         CALL READI(MAXGROWSTEPS)
C
C  Will stop the entire job if the total string
C  length for the growing strings or evolving strings method goes above {\it x}
C  times the total number of images. This usually means that something is going
C  wrong with the string. Default is 1000.
C
      ELSE IF (WORD.EQ.'MAXLENPERIM') THEN
         CALL READF(MAXLENPERIM)
C
C  Specifies the maximum value that the maximum step size
C  is allowed to rise to. The default value is $0.5$.
C
      ELSE IF (WORD.EQ.'MAXMAX') THEN
         CALL READF(MAXMAX)
C
C  MAXSTEP n specifies the maximum step size in real units      - default n=0.2
C  Applies to eigenvector-following and steepest-descent calculations.
C
      ELSE IF (WORD.EQ.'MAXSTEP') THEN
         CALL READF(MXSTP)
C
C  Maximum ts energy that is allowed to constitute a connection during the
C  Dijkstra connection procedure.
C
      ELSE IF (WORD.EQ.'MAXTSENERGY') THEN
         CALL READF(MAXTSENERGY)
C
C  MECCANO - an interpolation via rigid rods of variable length
C
      ELSE IF (WORD.EQ.'MECCANO') THEN
         MECCANOT=.TRUE.
         CALL READF(MECIMDENS) ! now an image density
         CALL READI(MECMAXIMAGES)  ! maximum number of images
         CALL READF(MECITDENS) ! iteration density
         CALL READI(MECMAXIT)  ! maximum number of iterations
         CALL READF(MECLAMBDA)
         CALL READF(MECDIST)
         CALL READF(MECRMSTOL)
         CALL READF(MECSTEP)
         CALL READF(MECDGUESS)
         CALL READI(MECUPDATE)
      ELSE IF (WORD.EQ.'MINMAX') THEN
         CALL READF(MINMAX)
!     ELSE IF (WORD.EQ.'MINBM') THEN
!        MINBMT = .TRUE.
!        CALL READI(MINBMNSAMP)
      ELSE IF (WORD.EQ.'MINBACKTCUT') THEN
         CALL READF(MINBACKTCUT)
C
C  MODE n  specifies the eigenvector to follow                  - default n=0
C
      ELSE IF (WORD.EQ.'MODE') THEN
         CALL READI(IVEC)
         IF (NITEMS.GT.2) THEN
            CALL READI(IVEC2)
         ELSE
C           IVEC2=IVEC
         ENDIF
C
C  Attempt to morph between endpoints by taking steps towards or
C  away from the endpoint finish.
C
      ELSE IF (WORD.EQ.'MORPH') THEN
         MORPHT=.TRUE.
         CALL READF(MORPHMXSTP)
         CALL READI(MNBFGSMAX1)
         CALL READI(MNBFGSMAX2)
         CALL READF(MORPHEMAX)
         CALL READF(MORPHERISE)
         CALL READI(MSTEPS)
         IF (MAXTSENERGY.EQ.1.0D100) MAXTSENERGY=MORPHEMAX
C
C  Movie dump for Paul Mortenson`s amber
C
      ELSE IF (WORD.EQ.'MOVIE') THEN
         MOVIE=.TRUE.
         OPEN (UNIT=27, FILE='amber.movie', STATUS='UNKNOWN')
C
C  MSEVB parameters - probably shouldn`t be changed on a regular basis
C
      ELSE IF (WORD.EQ.'MSEVBPARAMS') THEN
         IF (NITEMS.GT.1) CALL READI(shellsToCount)
         IF (NITEMS.GT.2) CALL READF(maxHbondLength)
         IF (NITEMS.GT.3) CALL READF(minHbondAngle)
         IF (NITEMS.GT.4) CALL READF(OOclash_sq)

      ELSE IF (WORD.EQ.'MSSTOCK') THEN

         MSSTOCKT = .TRUE.
         RBAAT    = .TRUE.
         CALL READI(NRBSITES)
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(RBSTLA(NRBSITES,3))
         ALLOCATE(DPMU(NRBSITES))

         NTSITES = NATOMS*NRBSITES/2

         DO J1 = 1, NRBSITES
            CALL READF(DPMU(J1))
         ENDDO

         IF (NITEMS > (NRBSITES+2)) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF

         CALL DEFMULTSTOCK()
         IF (PERMDIST) THEN ! correct all permutations allowed if perm.allow is not given explicitly
            IF (NPERMSIZE(1).EQ.NATOMS) NPERMSIZE(1)=NATOMS/2
         ENDIF

C
C  NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
C
C  Specifies a tight-binding potential for sodium, silver and lithium
C
      ELSE IF (WORD.EQ.'NATB') THEN
         NATBT=.TRUE.
      ELSE IF (WORD.EQ.'NATINT') THEN
         NATINT = .TRUE.
C
      ELSE IF (WORD.EQ.'NCAP') THEN
         NCAPT    = .TRUE.
         RBAAT    = .TRUE.
         HEIGHT   = 0.5D0
         CALL READF(CAPSRHO)
         CALL READF(CAPSEPS2)
         CALL READF(CAPSRAD)
         IF (NITEMS.GT.4) CALL READF(HEIGHT)
         NRBSITES = 6
         ALLOCATE(RBSITE(NRBSITES,3))
         NTSITES = NATOMS*NRBSITES/2
C
C  Nudged elastic band calculation using a maximum of NSTEPNEB steps with
C  NIMAGE images and RMS convergence criterion RMSNEB.
C
      ELSE IF (WORD.EQ.'NEB') THEN
         NEBT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NSTEPNEB)
         IF (NITEMS.GT.2) CALL READI(NIMAGE)
         IF (NITEMS.GT.3) CALL READF(RMSNEB)
         IF (NITEMS.GT.3) CALL READI(NEBMAG)
      ELSE IF (WORD == 'NEBK') THEN
         CALL READF(NEBK)
         NEBKFINAL=NEBK
         NEBKINITIAL=NEBK
         NEBFACTOR=1.01D0
         IF (NITEMS.GT.2) CALL READF(NEBKFINAL)
         IF (NITEMS.GT.3) CALL READF(NEBFACTOR)
C
C  Read dneb guess images from file GUESSFILE, default name guess.xyz
C
      ELSE IF (WORD.EQ.'NEBREADGUESS') THEN
         READGUESS=.TRUE.
         IF (NITEMS.GT.1) CALL READA(GUESSFILE)
C
C  Reseed DNEB images if they exceed a certain energy.
C
      ELSE IF (WORD.EQ.'NEBRESEED') THEN
         NEBRESEEDT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NEBRESEEDINT)
         IF (NITEMS.GT.2) CALL READF(NEBRESEEDEMAX)
         IF (NITEMS.GT.3) CALL READF(NEBRESEEDBMAX)
         IF (NITEMS.GT.4) CALL READF(NEBRESEEDDEL1)
         IF (NITEMS.GT.5) CALL READI(NEBRESEEDPOW1)
         IF (NITEMS.GT.6) CALL READF(NEBRESEEDDEL2)
         IF (NITEMS.GT.7) CALL READI(NEBRESEEDPOW2)
      ELSE IF (WORD == 'NEWCONNECT') THEN
         NEWCONNECTT = .TRUE.
         CONNECTT = .TRUE.
         OPTIMIZETS = .TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONMAX)
         IF (NITEMS.GT.2) CALL READI(NTRIESMAX)
         IF (NITEMS.GT.3) CALL READF(IMAGEDENSITY)
         IF (NITEMS.GT.4) CALL READF(ITERDENSITY)
         IF (NITEMS.GT.5) CALL READI(IMAGEMAX)
         IF (NITEMS.GT.6) CALL READF(IMAGEINCR)
         IF (NITEMS.GT.7) CALL READF(RMSTOL)
C
C  If NEWCONNECT is specified the values read below are only used for the first cycle.
C  If NEWNEB is used with OLDCONNECT then the values read on the NEWNEB line are
C  used in every cycle. If NEWCONNECT is used then a NEWNEB line isn;t necessary.
C
      ELSE IF (WORD == 'NEWNEB') THEN
         NEWNEBT=.TRUE.
         FCD=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NNNIMAGE)
         IF (NITEMS.GT.2) CALL READI(NITERMAX)
         IF (NITEMS.GT.3) CALL READF(RMSTOL)
C
C  NGUESS specifies the number of transition state guesses tried in GUESSTS for CHARMM
C  before switching back to NEB or NEWNEB.
C
      ELSE IF (WORD.EQ.'NGUESS') THEN
         CALL READI(NGUESS)
C
C  CHARMM related keyword to reject transition states
C  that connect two minima with different omega angles, i.e. to prevent cis-trans peptide
C  isomerisation.
C
      ELSE IF (WORD.EQ.'UACHIRAL') THEN
         UACHIRAL=.TRUE.

      ELSE IF (WORD.EQ.'NOCISTRANS') THEN
         NOCISTRANS=.TRUE.   ! is used in connect.f
         CHECKOMEGAT=.TRUE.  ! is used in NEWNEB
         IF (NITEMS.GT.1) CALL READF(MINOMEGA)
        
         IF (NITEMS.GT.2) THEN
            CALL READU(WW)
            IF (TRIM(ADJUSTL(WW))=='RNA') THEN
                NOCISTRANSRNA = .TRUE.
                write(*,*) ' keywords> NOCISTRANSRNA set to .TRUE.'
            ELSE IF (TRIM(ADJUSTL(WW))=='DNA') THEN
                NOCISTRANSDNA = .TRUE.
                write(*,*) ' keywords> NOCISTRANSDNA set to .TRUE.'
            ELSE IF (TRIM(ADJUSTL(WW))=='ALWAYS') THEN
                CHECKCISTRANSALWAYS = .TRUE.
            ELSE IF (TRIM(ADJUSTL(WW))=='ALWAYSRNA') THEN
                CHECKCISTRANSALWAYSRNA = .TRUE.
            ELSE IF (TRIM(ADJUSTL(WW))=='ALWAYSDNA') THEN
                CHECKCISTRANSALWAYSDNA = .TRUE.
            ELSE 
                WRITE(*,*) ' keywords> ERROR - currently no other nocistrans options implemented than for RNA and DNA'
            ENDIF
         ENDIF
C
C  No frequencies should be evaluated or placed in the path.info file.
C
      ELSE IF (WORD.EQ.'NOFRQS') THEN
         NOFRQS=.TRUE.
C
C  No Hessian should be calculated during geometry optimisation.
C
      ELSE IF (WORD.EQ.'NOHESS') THEN
         NOHESS=.TRUE.
C
C  If NOIT is true and we have a Hessian then use DSYEVR to calculate eigenvectors
C
      ELSE IF (WORD.EQ.'NOINTNEWT') THEN
         ! dont use newtons method to converge internals back transform
         INTNEWT = .FALSE.

      ELSE IF (WORD.EQ.'NOIT') THEN
         NOIT=.TRUE.
C
C  For the growing string or evolving string double-ended
C  transition state search methods, instead of using L-BFGS optimization to
C  evolve the strings, simply take steps in the direction of the perpendicular force. 
C

      ELSE IF (WORD.EQ.'NOLBFGS') THEN
         NOLBFGS = .TRUE.
      ELSE IF (WORD == 'NONEBMIND') THEN
          NEBMIND=.FALSE.
          PRINT *, 'keywords> Structures supplied to NEB will NOT be put in the closest coincidence'
C
C  NONLOCAL x y z factors for averaged Gaussian, Morse type 1 and Morse
C                 type 2 potentials to include                  - default 0 0 0  
C
      ELSE IF (WORD.EQ.'NONLOCAL') THEN
         CALL READF(GFRACTION)
         IF (NITEMS.GT.2) THEN
            CALL READF(MFRACTION1)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(MFRACTION2)
         ENDIF
         FTEST=.TRUE.
         IF ((GFRACTION.EQ.0.0D0).AND.(MFRACTION1.EQ.0.0D0).AND.
     1       (MFRACTION2.EQ.0.0D0)) FTEST=.FALSE.

      ELSE IF (WORD.EQ.'NOPERMPROCHIRAL') THEN
         NOPERMPROCHIRAL = .TRUE.

C
C  Reduce printing of coordinates.
C
      ELSE IF (WORD.EQ.'NOPOINTS') THEN
         PRINTPTS=.FALSE.
C
C  Used in CHARMM transition state guessing procedure
C  together with TWISTTYPE. Setting randomcutoff very large prevents random
C  steps, and is recommended. 

C
      ELSE IF (WORD.EQ.'NORANDOM') THEN
         NORANDOM=.TRUE.
         IF (NITEMS.GT.1) CALL READF(RANDOMCUTOFF)
C
C  Whether to put periodic images back in the primary supercell.
C
      ELSE IF (WORD .EQ. 'NORESET') THEN
         NORESET=.TRUE.

      ELSE IF (WORD.EQ.'NTIP') THEN
         CALL READI(TIPID)
         IF (TIPID /= 4) THEN
            PRINT *, 'NOT YET INCLUDED'
            STOP
         ENDIF
         NTIPT = .TRUE.
         RBAAT = .TRUE.
         NRBSITES = 4
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(STCHRG(NRBSITES))
         NTSITES = NATOMS*NRBSITES/2
         IF (PERMDIST) THEN ! correct all permutations allowed if perm.allow is not given explicitly
            IF (NPERMSIZE(1).EQ.NATOMS) NPERMSIZE(1)=NATOMS/2
         ENDIF
C
C  OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
      ELSE IF (WORD.EQ.'ODIHE') THEN
         ODIHET=.TRUE.
         WRITE(*,'(A)') 'ODIHE set: dihedral-angle order parameter will be calculated'
         WRITE(*,'(A)') 'using the reference structure supplied in ref.crd'
C
C  Add an octahedral field to the potential of magnitude FOH.
C
      ELSE IF (WORD.EQ.'OH') THEN
         FIELDT=.TRUE.
         OHT=.TRUE.
         CALL READF(FOH)

      ELSE IF (WORD.EQ.'OLDINTMINPERM') THEN
         INTMINPERMT = .TRUE.
         OLDINTMINPERMT=.TRUE.

C
C  ONETEP tells the program to read derivative information in 
C         ONETEP format.                                        - default FALSE
C
      ELSE IF ((WORD.EQ.'ONETEP').OR.(WORD.EQ.'ONETEPC')) THEN
         ONETEP=.TRUE.
         IF (WORD.EQ.'ONETEP') DFTP=.TRUE.
         IF (NITEMS.GT.2) THEN
            CALL READA(ONETEPJOB)
            CALL READA(SYS)
            ONETEPJOB=TRIM(ADJUSTL(ONETEPJOB)) // ' ' // TRIM(ADJUSTL(SYS)) // ' >& ' // TRIM(ADJUSTL(SYS)) // '.onetep'
         ELSE
            WRITE(*,'(A)') 'keywords> ERROR - ONETEP job or system unspecified'
            CALL FLUSH(6,ISTAT)
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 24
            ENDIF
         ENDDO
24       CONTINUE
C
C  Optimise TS with SQVV
C
      ELSE IF (WORD == 'OPTIMIZETS') THEN
          OPTIMIZETS=.TRUE.
C
C  Calculates order parameters and theire derivatives wrt normal modes at the end of a geometry optimisation.
C  The 1st argument is the number of order parameters to be calculated. The next arguments then specify
C  the order parameters (defined by a 4 letters) and, if necessary, further information regarding this
C  order parameter can be given. If such details are not required, set them to -9999.
C  Following order parameters are currently supported: DIHEdral angles for CHARMM.
C
      ELSE IF (WORD == 'ORDERPARAM') THEN
          ORDERPARAMT=.TRUE.
          CALL READI(NORDER)
          ALLOCATE(WHICHORDER(NORDER),ORDERNUM(NORDER))
          DO J1=1,NORDER
             CALL READA(WHICHORDER(J1))
             ORDERNUM(J1)=-9999
             CALL READI(ORDERNUM(J1))
          ENDDO
C
C  Remove overall trans/rot with SQVV
C
      ELSE IF (WORD == 'ORT') THEN
         ORT = .TRUE.
      ELSE IF (WORD.EQ.'OSASA') THEN
         OSASAT=.TRUE.
         CALL READF(RPRO)
         WRITE(*,'(A)') 'OSASA set: solvent accessible surface area order parameter will be calculated'
         WRITE(*,'(A,F3.1)') 'using probe radius ',RPRO
C
C  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C
      ELSE IF (WORD .EQ. 'PAHA') THEN

         CALL READI(PAHID)

         IF (PAHID == 1) THEN
            NRBSITES = 12
         ELSEIF (PAHID == 2) THEN
            NRBSITES = 18
         ELSEIF (PAHID == 3) THEN
            NRBSITES = 24
         ELSEIF (PAHID == 4) THEN
            NRBSITES = 26
         ENDIF

         PAHAT    = .TRUE.
         RBAAT    = .TRUE.
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(RBSTLA(NRBSITES,3))
         ALLOCATE(STCHRG(NRBSITES))

         NTSITES = (NATOMS/2)*NRBSITES

         CALL DEFPAHA()

         IF (PAHID == 1) THEN
            NCARBON  = 6
            CALL DEFBENZENE()
         ELSEIF (PAHID == 2) THEN
            NCARBON  = 10
            CALL DEFNAPHTHALENE()
         ELSEIF (PAHID == 3) THEN
            NCARBON  = 14
            CALL DEFANTHRACENE()
         ELSEIF (PAHID == 4) THEN
            NCARBON  = 16
            CALL DEFPYRENE()
         ENDIF

      ELSE IF (WORD.EQ.'PARALLEL') THEN
         PARALLEL=.TRUE.
         CALL READA(NPROC)
C
C  PARAMS n1 n2 ... up to seven real input parameters used for the
C                   following atom types:
C  AX: Z*
C  M:  rho
C  MV: rho, delta
C  ME: N, M, BOXLENGTHS X, Y, Z AND CUTOFF (N, M ARE READ DOUBLE PRECISION)
C  JM: box lengths x, y, z and cutoff
C  SC: box lengths x, y, z and cutoff (epsilon, c, sigma are read from SCparams)
C  P6: box lengths x, y, z and cutoff
C  AU: epsilon, c, sigma
C  AG: epsilon, c, sigma
C  NI: epsilon, c, sigma
C
      ELSE IF (WORD.EQ.'PARAMS') THEN
         CALL READF(PARAM1)
         GALPHA=PARAM1
         MALPHA1=PARAM1
         MALPHA2=PARAM1
         IF (NITEMS.GT.2) THEN
            CALL READF(PARAM2)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(PARAM3)
         ENDIF
         IF (NITEMS.GT.4) THEN
            CALL READF(PARAM4)
         ENDIF
         IF (NITEMS.GT.5) THEN
            CALL READF(PARAM5)
         ENDIF
         IF (NITEMS.GT.6) THEN
            CALL READF(PARAM6)
         ENDIF
         IF (NITEMS.GT.7) THEN
            CALL READF(PARAM7)
         ENDIF

      ELSE IF (WORD.EQ.'PATCHYD') THEN

         PATCHYDT   = .TRUE.
         RBAAT  = .TRUE.
         CALL READI(NRBSITES)

         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(RBSTLA(NRBSITES,3))

         CALL DEFPATCHES

         NTSITES = NATOMS*NRBSITES/2

C
C  PATH specifies calculation of the pathway connecting two minima from the transition
C  state specified in odata. NPATHFRAME is the number of points files to save on either
C  side. A complete xyz file is printed to path.xyz and the energy as a function of
C  path length is printed to file EofS.
C  Movies generated in this way tend to move too fast for the interesting bits, and too
C  slow around stationary points. Specify FRAMEEDIFF to give a lower bound to the energy difference
C  between frames for which the structure is considered different.
C
      ELSE IF (WORD.EQ.'PATH') THEN   
         PATHT=.TRUE.
         IF (NITEMS.GT.1) THEN
              CALL READI(NPATHFRAME)
!             if (NPATHFRAME<3) THEN
!                  PRINT *, 'Number of path frames cannot be less than 3 - stop'
!                  stop
!             ELSE IF (NPATHFRAME>3) THEN
!                  IF (.NOT.PRINTPTS) THEN
!                       PRINT *, 'Number of path frames is more than 3 - dumping all points!'
!                       PRINTPTS=.TRUE.
!                  ENDIF
!             ENDIF
         ENDIF
         IF (NITEMS.GT.2) CALL READF(FRAMEEDIFF)
         IF (NPATHFRAME.LE.0) PRINTPTS=.FALSE.
         IF (NPATHFRAME.GT.0) PRINTPTS=.TRUE.
      ELSE IF (WORD.EQ.'PERMDIHE') THEN
!
!  PATHSDSTEPS sets the number of SD steps allowed at the beginning of a path
!  calculation. We switch to LBFGS from RKMIN, BSMIN and SEARCH INR methods if
!  they don't converge in PATHSDSTEPS steps. If not set then the default is NSTEPS.
!
      ELSE IF (WORD.EQ.'PATHSDSTEPS') THEN
         CALL READI(PATHSDSTEPS)
      ELSE IF (WORD.EQ.'PERMDIHE') THEN
         PERMDIHET=.TRUE.
         DO J1=1,NITEMS-1
            CALL READI(NDUM)
            PERMDIHE(J1)=NDUM
         ENDDO
         NPERMDIHE=NITEMS-1
         DO J1=1,NITEMS-1
            PRINT *,'PERMDIHE',PERMDIHE(J1)
         ENDDO
C
C  Whether to optimise the permutational isomers in assessing optimal
C  alignment.
C
      ELSE IF ((WORD.EQ.'PERMDIST').OR.(WORD.EQ.'PERMDISTINIT').OR.(WORD.EQ.'LOCALPERMDIST')) THEN
         PERMDIST=.TRUE.
         IF (WORD.EQ.'PERMDISTINIT') PERMDISTINIT=.TRUE.
         IF (WORD.EQ.'LOCALPERMDIST') THEN
            LOCALPERMDIST=.TRUE.
            IF (NITEMS.GT.1) CALL READF(LPDGEOMDIFFTOL)
            IF (NITEMS.GT.2) CALL READF(RBCUTOFF)
            IF (NITEMS.GT.3) CALL READI(NRBTRIES)
            PRINT '(A)',' keyword> Local rigid body permutational alignment:'
            PRINT '(2(A,F12.4),A,I6)','          distance tolerance=',LPDGEOMDIFFTOL,' cutoff=',RBCUTOFF, 
     &                  ' number of passes through alignment phase=',NRBTRIES
         ENDIF

         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
!           ALLOCATE(NPERMSIZE(NATOMS),PERMGROUP(NATOMS),NSWAP(NATOMS),SWAP1(NATOMS,3),SWAP2(NATOMS,3))
            ALLOCATE(NPERMSIZE(3*NATOMS),PERMGROUP(3*NATOMS),NSETS(3*NATOMS),SETS(NATOMS,3))
!
!  The above dimensions were fixed at NATOMS because:
!  (a) Atoms were not allowed to appear in more than one group.
!  (b) The maximum number of pair exchanges associated with a group is three.
!
! However, for flexible water models we need to exchange all waters,
! and we can exchange H's attached to the same O. The dimension required
! becomes 3*NATOMS
!
            NDUMMY=1
            DO J1=1,NPERMGROUP
               READ(1,*) NPERMSIZE(J1),NSETS(J1)
!
!  Sanity checks!
!
               IF (NSETS(J1).GT.3) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1),' is > 3'
                  STOP
               ENDIF
!              IF (NDUMMY+NPERMSIZE(J1)-1.GT.NATOMS) THEN
               IF (NDUMMY+NPERMSIZE(J1)-1.GT.3*NATOMS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J2=1,NSETS(J1)),
!    &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
!
!  And another sanity check! This condition is now allowed.
!
!           DO J1=1,NDUMMY
!              DO J2=J1+1,NDUMMY
!                 IF (PERMGROUP(J2).EQ.PERMGROUP(J1)) THEN
!                    PRINT '(2(A,I8))','keyword> ERROR - atom ',PERMGROUP(J1),' appears more than once'
!                    STOP
!                 ENDIF
!              ENDDO
!           ENDDO
            CLOSE(1)
!
!  And yet another!
!
            IF (NFREEZE.GT.0) THEN
               NDUMMY=0
               DO J1=1,NPERMGROUP
                  DO J2=1,NPERMSIZE(J1)
                     IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                        PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                        STOP
                     ENDIF
                  ENDDO
                  NDUMMY=NDUMMY+NPERMSIZE(J1)
               ENDDO
            ENDIF
         ELSE
            ALLOCATE(NPERMSIZE(NATOMS),PERMGROUP(NATOMS),NSETS(NATOMS),SETS(NATOMS,2))
            NSETS(1:NATOMS)=0
            NPERMGROUP=1 ! all atoms can be permuted - default
            NPERMSIZE(1)=NATOMS ! all atoms can be permuted - default
            IF (RBAAT) NPERMSIZE(1)=NATOMS/2 ! for rigid bodies
            DO J1=1,NPERMSIZE(1)
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         PRINT '(A,I6)',' keyword> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         IF (DEBUG) THEN
            DO J1=1,NPERMGROUP
               PRINT '(A,3(I6,A))',' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                    NSETS(J1),' additional atom sets:'
               WRITE(*,'(22I6)',ADVANCE='NO') PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1) 
               IF (NSETS(J1).GT.0) THEN
                  WRITE(*,'(A)',ADVANCE='NO') ' with '
                  DO J2=1,NSETS(J1)
                     DO J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1
                        WRITE(*,'(I6)',ADVANCE='NO') SETS(PERMGROUP(J3),J2)
                        IF (J3.LT.NDUMMY+NPERMSIZE(J1)-1) WRITE(*,'(A3)',ADVANCE='NO') ' / '
                     ENDDO
                     IF (J2.LT.NSETS(J1)) WRITE(*,'(A3)',ADVANCE='NO') ' ; '
                  ENDDO
               ENDIF
               PRINT *,' '
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
         ENDIF
C
C  CHARMM and UNRES dihedral angle perturbation specification.
C  Performs random, {\tt GMIN}-style twists before starting optimisation.
C

      ELSE IF (WORD.EQ.'PERTDIHE') THEN
         PERTDIHET=.TRUE.
         CALL READF(CHPMAX)
         CHPMIN=CHPMAX
         CALL READF(CHNMIN)
         CALL READF(CHNMAX)
         CALL READI(ISEED)
C        PRINT *,'CHPMIN,CHPMAX,CHNMIN,CHNMAX',CHPMIN,CHPMAX,CHNMIN,CHNMAX

C
C  SQVV keyword.
C
      ELSE IF (WORD == 'PRINTOPTIMIZETS') THEN
          PRINTOPTIMIZETS=.TRUE.
C
C  For the GS and ES double-ended transition state
C  search methods, if using {\it FIXATMS\/} to zero some coordinates of the
C  forces to avoid overall translation and rotation, this keyword will rotate
C  the start and end points so that those coordinates are zero in both.
C
      ELSE IF (WORD.EQ.'PREROTATE') THEN
         PREROTATE = .TRUE.
C
C
C  PRESSURE tells the program to perform a constant pressure optimisation
C           for SC, ME and P6 with periodic boundary conditions - default off
C
      ELSE IF (WORD.EQ.'PRESSURE') THEN
         PRESSURE=.TRUE.
C
C  PRINT n sets the value of IPRNT                              - default n=0
C
      ELSE IF (WORD.EQ.'PRINT') THEN
         CALL READI(IPRNT)
C
C  Print ground state coefficients - only valid for MSEVB potential
C
      ELSE IF (WORD.EQ.'PRINTCOEFFICIENTS') THEN
         printCoefficients=.TRUE.

C print out info on coordinates and stop; for debugging internals
      ELSE IF (WORD.EQ.'PRINTCOORDS') THEN
         PRINTCOORDS = .TRUE.
C
C
C
C
C  Keyword for applied static force.
C
      ELSE IF (WORD.EQ.'PULL') THEN
         PULLT=.TRUE.
         CALL READI(PATOM1)
         CALL READI(PATOM2)
         CALL READF(PFORCE)
         IF (PFORCE.EQ.0.0D0) THEN
            WRITE(*,'(A,I6,A,I6,A,G20.10)') ' keyword> WARNING *** Pulling force is zero, turning off pulling directive'
            PULLT=.FALSE.
         ELSE
            WRITE(*,'(A,I6,A,I6,A,G20.10)') ' keyword> Pulling atoms ',PATOM1,' and ',PATOM2,' force=',PFORCE
         ENDIF
C
C  PUSHCUT sets the threshold for when a PUSHOFF will be applied, i.e.
C  the RMS force must be less than PUSHCUT.
C
      ELSE IF (WORD .EQ. 'PUSHCUT') THEN
         CALL READF(PUSHCUT)
C
C  PUSHOFF x sets the magnitude of the step away from a converged 
C            transition state if detected on the first cycle of 
C            a minimisation                                     - default x=0.01
C
      ELSE IF (WORD .EQ. 'PUSHOFF') THEN
         CALL READF(PUSHOFF)
C
C  PV
C
      ELSE IF (WORD.EQ.'PV') THEN
         PV=.TRUE.
         IF (NITEMS.GT.1) CALL READF(PRESS)
         IF (NITEMS.GT.2) CALL READF(PVCONV)
         IF (NITEMS.GT.3) CALL READF(PVTOL)
         IF (NITEMS.GT.4) CALL READI(PVSTEPS)
      ELSE IF (WORD.EQ.'PVTS') THEN
         PV=.TRUE.
         PVTS=.TRUE.
         NBOXTS=1
         IF (NITEMS.GT.1) CALL READF(PRESS)
         IF (NITEMS.GT.2) CALL READF(PVCONV)
         IF (NITEMS.GT.3) CALL READF(PVTOL)
         IF (NITEMS.GT.4) CALL READI(PVSTEPS)
         IF (NITEMS.GT.5) CALL READI(NBOXTS)
         WRITE(*,'(A,I5)') ' Searching uphill for a transition state in box length coordinate ',NBOXTS

      ELSE IF (WORD.EQ.'PYG') THEN
         NRBSITES = 1
         ALLOCATE(RBSITE(NRBSITES,3))
         PYGT  = .TRUE.
         RBAAT = .TRUE.
         CALL READF(PYA1(1))
         CALL READF(PYA1(2))
         CALL READF(PYA1(3))
         CALL READF(PYA2(1))
         CALL READF(PYA2(2))
         CALL READF(PYA2(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)

         IF (PYA1(1) == PYA2(1) .AND. PYA1(2) == PYA2(2) .AND. PYA1(3) == PYA2(3)) THEN
            RADIFT = .FALSE.
         ELSE
            RADIFT = .TRUE.
         ENDIF
! sf344> PY potential and extra LJ site
      ELSE IF (WORD.EQ.'PYBINARY') THEN
         NRBSITES = 1
         ALLOCATE(RBSITE(NRBSITES,3))
         PYBINARYT=.TRUE.
         ANGLEAXIS2=.TRUE.
         RBAAT=.TRUE.
         RADIFT=.TRUE.
         CALL READI(PYBINARYTYPE1)
         CALL READF(PYA11(1))
         CALL READF(PYA11(2))
         CALL READF(PYA11(3))
         CALL READF(PYA21(1))
         CALL READF(PYA21(2))
         CALL READF(PYA21(3))
         CALL READF(PYA12(1))
         CALL READF(PYA12(2))
         CALL READF(PYA12(3))
         CALL READF(PYA22(1))
         CALL READF(PYA22(2))
         CALL READF(PYA22(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)

         IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMS/2,3))
         IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMS/2,3))
         DO J1=1,NATOMS/2
          IF(J1<=PYBINARYTYPE1) THEN
           PYA1bin(J1,:)=PYA11(:)
           PYA2bin(J1,:)=PYA21(:)
          ELSE
           PYA1bin(J1,:)=PYA12(:)
           PYA2bin(J1,:)=PYA22(:)
          END IF
         END DO
      ELSE IF (WORD.EQ.'CLOSESTALIGNMENT') THEN
         CLOSESTALIGNMENT=.TRUE.
         WRITE(*,*) 'Putting structures into closest alignment, then exiting'
      ELSE IF (WORD.EQ.'PYGPERIODIC') THEN
         PYGPERIODICT = .TRUE.
         ANGLEAXIS2=.TRUE.
         RBAAT=.TRUE.
         CALL READF(PYA1(1))
         CALL READF(PYA1(2))
         CALL READF(PYA1(3))
         CALL READF(PYA2(1))
         CALL READF(PYA2(2))
         CALL READF(PYA2(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)

         IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMS/2,3))
         IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMS/2,3))
         DO J1=1,NATOMS/2
           PYA1bin(J1,:)=PYA1(:)
           PYA2bin(J1,:)=PYA2(:)
         END DO
         IF (PYA1(1) == PYA2(1) .AND. PYA1(2) == PYA2(2) .AND. PYA1(3) == PYA2(3)) THEN
            RADIFT = .FALSE.
         ELSE
            RADIFT = .TRUE.
         ENDIF

         IF (NITEMS.GT.9) THEN
            CALL READF(PCUTOFF)
            PARAMONOVCUTOFF=.TRUE.
            PCUTOFF=PCUTOFF*PYSIGNOT
            write (MYUNIT,*) "PY Potential. PCutoff ON:",PCUTOFF
         ENDIF
         IF (NITEMS.GT.10) THEN
! control which dimensions have periodic boundaries with a string 'XYZ', always put x before y before z.
! eg ...  Xz 20 30  specifies PBC on X and Z directions.  The X box size will be 20, the Z box size 30
            CALL READA(PBC)
            write (*,*) "PBCs are: ",PBC
            BOXLX=0
            BOXLY=0
            BOXLZ=0
            IF (SCAN(PBC,'Xx').NE.0) THEN
                PARAMONOVPBCX=.TRUE.
                CALL READF(BOXLX)       ! BOXLX is a scaling factor, not the actual box length!
                BOXLX=BOXLX*PCUTOFF     ! now BOXLX is the actual box length
                write(*,*) "PY Periodic Boundary Condition X active. BOXLX:",BOXLX
            ENDIF
            IF (SCAN(PBC,'Yy').NE.0) THEN
                PARAMONOVPBCY=.TRUE.
                CALL READF(BOXLY)
                BOXLY=BOXLY*PCUTOFF
                write(*,*) "PY Periodic Boundary Condition Y active. BOXLY:",BOXLY
            ENDIF
            IF (SCAN(PBC,'Zz').NE.0) THEN
                PARAMONOVPBCZ=.TRUE.
                CALL READF(BOXLZ)
                BOXLZ=BOXLZ*PCUTOFF
                write(*,*) "PY Periodic Boundary Condition Z active. BOXLZ",BOXLZ
            ENDIF
         ENDIF
         ALLOCATE(RBSITE(NRBSITES,3))
C
C  QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
C
C
C  qSPCFw  flexible water model introduced by Paesani et al. (JCP 125, 184507 (2006))
C  Coded by Javier.
C
      ELSE IF (WORD.EQ.'QSPCFW') THEN
         QSPCFWT=.TRUE.
C
C  qTIP4PF flexible water model introduced by Habershon et al. (JCP 131, 024501 (2009))
C  Coded by Javier.
C
      ELSE IF (WORD.EQ.'QTIP4PF') THEN
         QTIP4PFT=.TRUE.
C
C  RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C
C
C  Spherical container
C
      ELSE IF (WORD.EQ.'RADIUS') THEN
         CONTAINER=.TRUE.
         CALL READF(RADIUS)
         RADIUS=RADIUS**2
C
C  integer seed for random number generator.
C
      ELSE IF (WORD.EQ.'RANSEED') THEN
         CALL READI(NDUM)
         CALL SDPRND(NDUM)
         IF ((NDUM.LT.0).OR.(NDUM.GT.9999)) THEN
C
C  if we ever need more than 20000 searches from the same minimum
C  then this could be a problem
C
            DO J1=1,3*NATOMS
                RANDOM=DPRAND()
            ENDDO
         ENDIF
         WRITE(*,'(A,I6)') ' SETTINGS Random number generator seed=',NDUM
C
C  TVB: Requests to print out pathway parameters necessary to calculate catastrophe
C  ratios. Affects path routine only.
C
      ELSE IF (WORD.EQ.'RATIOS') THEN
               RATIOS=.TRUE.
C
C  RBSYM defines the internal symmetry operations for each sort of rigid body
C  coded via RBAAT.
C
      ELSE IF (WORD.EQ.'RBSYM') THEN
         RBSYMT=.TRUE.
         INQUIRE(FILE='rbsymops',EXIST=RBSYMTEST)
         IF (RBSYMTEST) THEN
            OPEN(UNIT=1,FILE='rbsymops',STATUS='OLD')
            READ(1,*) NRBGROUP
            ALLOCATE(RBOPS(4,NRBGROUP))
            READ(1,*) ((RBOPS(J1,J2),J1=1,4),J2=1,NRBGROUP)
            PRINT '(A,I6)',' keywords> number of symmetry operations for rigid body=',NRBGROUP
            DO J1=1,NRBGROUP
               PRINT '(A,I6)',' keywords> rigid-body symmetry operation', J1
               RBOPS(4,J1) = RBOPS(4,J1)*ATAN(1.D0)/45.D0
               PRINT '(3F20.10)',RBOPS(1:4,J1)
            ENDDO
         ELSE
            PRINT '(A)',' keywords> ERROR *** missing file rbsymops'
            STOP
         ENDIF
C
C  If READPATH is specified with CALCRATES then the rates are calculated from the 
C  information in an existing path.info file without any stationary point searches.
C
      ELSE IF (WORD.EQ.'READPATH') THEN
         READPATH=.TRUE.
C
C  If READSP is true then OPTIM will read minima and ts data in the pathsample format
C
      ELSE IF (WORD.EQ.'READSP') THEN
         READSP=.TRUE.
C
C  READHESS tells the program to read a Hessian at the first step.
C
      ELSE IF (WORD .EQ. 'READHESS') THEN
        READHESS=.TRUE.
C
C  READVEC "file" reads the eigenvalue and associated eigenvector corresponding 
C  to the reaction coordinate for use in a pathway calculation. The format
C  is the same as that used for vector.dump. If there is more than one vector
C  in the file the program reads down to the last entry.
C
      ELSE IF (WORD(1:7) .EQ. 'READVEC') THEN
         READV=.TRUE.
C     ELSE IF (WORD.EQ.'REBUILDSC') THEN
C        CALL READF(REBUILDSC)
C
C  sf344> read in coordinates from path.xyz files for rigid bodies, and
C         bring the frames in the best alignment
C
      ELSE IF (WORD.EQ.'REALIGNXYZ') THEN
         REALIGNXYZ=.TRUE.
C
C  Whether to use a redopoints file if it exists.
C
      ELSE IF (WORD.EQ.'REDOPATH') THEN
         REDOPATH=.TRUE.
         IF (NITEMS.GT.1) CALL READF(REDOK)
         IF (NITEMS.GT.2) CALL READF(REDOFRAC)
C
C  Whether to use a redopoints file if it exists.
C
      ELSE IF (WORD.EQ.'REDOPATHNEB') THEN
         REDOPATHNEB=.TRUE.
         REDOPATH=.TRUE.
         FREEZENODEST=.TRUE.
         FREEZETOL=-1.0D0
         IF (NITEMS.GT.1) CALL READI(REDOBFGSSTEPS)
C
C  Whether to use path.<n>.xyz files in the current directory
C
      ELSE IF (WORD.EQ.'REDOPATHXYZ') THEN
         REDOPATHXYZ=.TRUE.
         REDOPATH=.TRUE.
         IF (NITEMS.GT.1) CALL READF(REDOK)
         IF (NITEMS.GT.2) CALL READF(REDOFRAC)
C
C Whether to reduce the bond lengths for side chains during the connection runs.
C To be used together with CHARMM (and AMBER not yet).
C
      ELSE IF (WORD.EQ.'REDUCEDBONDLENGTH') THEN
         REDUCEDBONDLENGTHT=.TRUE.
         CALL READF(BLFACTOR)
         IF (NITEMS.GT.2) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'CB') CBT=.TRUE.
C
C  Specifies that the eigenvector to be followed should be reoptimised
C  in a BFGSTS search after the EF step and before the tangent space minimisation.
C  This is probably not a good idea.
C
      ELSE IF (WORD.EQ.'REOPT') THEN
         REOPT=.TRUE.
      ELSE IF (WORD.EQ.'REOPTIMISEENDPOINTS') THEN
         REOPTIMISEENDPOINTS=.TRUE.
C
C  coordinates to orthogonalise search directions to are to be found in 
C  points.repel
C
      ELSE IF (WORD.EQ.'REPELTS') THEN
         REPELTST=.TRUE.
         IF (NITEMS.GT.1) CALL READF(REPELPUSH)
C
C  RESIZE x scales the radial distances by x on the first 
C           step only                                           - default n=1    
C
      ELSE IF (WORD .EQ. 'RESIZE') THEN
        CALL READF(RESIZE)

C specifies additional rings other than the usual ones in 
C PHE, PRO, TYR, HIS, and TRP residues
      ELSE IF (WORD.EQ.'RING') THEN
         NURINGS = NURINGS + 1
         IF (NITEMS.EQ.6) THEN            
            URINGS(NURINGS,0) = 5
            DO J1 = 1,5
               CALL READI(URINGS(NURINGS,J1))
            ENDDO
         ELSE IF (NITEMS.EQ.7) THEN
            URINGS(NURINGS,0) = 6
            DO J1 = 1,6
               CALL READI(URINGS(NURINGS,J1))
            ENDDO
         ENDIF
C
C  RINGPOLYMER specifies a ring polymer system with harmonic springs between
C  NRP images of the same system that generally have different geometries.
C  RPSYSTEM is a string specifying the system, e.g. LJ.
C  RPIMAGES is the number of RP images.
C  RPBETA is 1/kT in reduced units.
C  RINGPOLYMER keyword takes the place of POINTS and must be the last
C  keyword in the odata file before the points.
C
      ELSE IF (WORD.EQ.'RINGPOLYMER') THEN
         RINGPOLYMERT=.TRUE.
         CALL READA(RPSYSTEM)
         CALL READI(RPIMAGES)
         CALL READF(RPBETA)
C
C  Sanity checks.
C
         TEMPSTRING=TRIM(ADJUSTL(RPSYSTEM))
         IF (TEMPSTRING(1:2).EQ.' ') THEN
            PRINT '(A)','keyword> ERROR *** Ring polymer potential type is not set'
         ENDIF
         IF (RPIMAGES.LT.1) THEN
            PRINT '(A)','keyword> ERROR *** Ring polymer images too small, value is ',RPIMAGES
         ENDIF
         RETURN
C
C  RKMIN calculates a steepest-descent path using gradient only information
C  with convergence criterion GMAX for the RMS force and initial precision
C  EPS. A fifth order Runga-Kutta algorithm is used.
C
      ELSE IF (WORD.EQ.'RKMIN') THEN
         RKMIN=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GMAX)
         IF (NITEMS.GT.2) CALL READF(EPS)
C
C  ROT [JZ n or OMEGA n] sets the value of J_z, the angular 
C                          momentum about the z axis or
C                          OMEGA, the corresponding angular velocity
C
      ELSE IF (WORD .EQ. 'ROT') THEN
        RTEST=.TRUE.
        CALL READU(WORD)
        CALL READF(XX)
        IF (WORD.EQ.'JZ') THEN
           JZ=XX
        ELSE
           OMEGA=XX
        ENDIF
C
C fix linear polymer at its ends
C
      ELSE IF (WORD .EQ. 'RPFIX') THEN
         RPFIXT=.TRUE.
         print *, 'fixed ends'
C
C make ring polymer system into linear polymer
C
      ELSE IF (WORD .EQ. 'RPLINEAR') THEN
         RPCYCLICT=.FALSE.
         print *, 'use linear polymer'
C
C  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
C
C
C  Save candidate TS`s in SQVV run.
C
      ELSE IF (WORD == 'SAVECANDIDATES') THEN
          SAVECANDIDATES=.TRUE.
C
C  SCALE n sets the value of ISTCRT                             - default n=10
C
      ELSE IF (WORD.EQ.'SCALE') THEN
         CALL READI(ISTCRT)
C
C  Specify that we are running in a SCore environment. Currently never used.
C
         ELSE IF (WORD.EQ.'SCORE_QUEUE') THEN
            SCORE_QUEUE = .TRUE.
C
C  SEARCH specifies the value of INR, i.e. the search type.     - default n=0
C
      ELSE IF (WORD.EQ.'SEARCH') THEN
         CALL READI(INR)
C
C  Eigenvalue shift parameter.
C
      ELSE IF (WORD .EQ. 'SHIFT') THEN
         CALL READF(SHIFTV)
C
C  Parameters for Edwin;s SiO2 model
C
      ELSE IF (WORD.EQ.'SIO2') THEN
         SIO2T=.TRUE.
         CALL READF(PARAM1)
         IF (NITEMS.GT.2) THEN
            CALL READF(PARAM2)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(PARAM3)
         ENDIF
         IF (NITEMS.GT.4) THEN
            CALL READF(PARAM4)
         ENDIF
      ELSE IF (WORD.EQ.'SIO2C6') THEN
         SIO2C6T=.TRUE.
         CALL READF(C6PP)
         CALL READF(C6MM)
         CALL READF(C6PM)
C
C  SQVV allows the first NIterSQVVGuessMax DNEB iterations to be done using
C  SQVV - switches to LBFGS minimisation after NIterSQVVGuessMax iterations
C         or if the RMS force goes below SQVVGuessRMSTol.
C
      ELSE IF (WORD == 'SQVV') THEN
         SQVVGUESS=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NITERSQVVGUESSMAX)
         IF (NITEMS.GT.2) CALL READF(SQVVGUESSRMSTOL)
C
C  NSTEPMIN sets the minimum number of steps allowed before convergence.
C 
      ELSE IF (WORD .EQ. 'STEPMIN') THEN
         CALL READI(NSTEPMIN)
C
C  STEPS n sets the number of optimisation steps to perform
C          per call to OPTIM                                    - default n=1     
C  If BFGSSTEPS is not specified then it is set to the same value as NSTEPS
C
      ELSE IF (WORD .EQ. 'STEPS') THEN
        CALL READI(NSTEPS)
        IF (BFGSSTEPS.EQ.1) BFGSSTEPS=NSTEPS
C
C  Stillinger-David water potential - coded by Jeremy Richardson
C
      ELSE IF (WORD.EQ.'SD') THEN
         SDT=.TRUE.
         CALL READI(SDOXYGEN)
         CALL READI(SDHYDROGEN)
         CALL READI(SDCHARGE)
         IF (SDOXYGEN*SDHYDROGEN.EQ.0) THEN
            PRINT '(A,2I6)', ' keyword> ERROR *** number of SD oxygens and hydrogens=',SDOXYGEN,SDHYDROGEN
            STOP
         ENDIF
      ELSE IF (WORD.EQ.'STOCK') THEN
         STOCKT=.TRUE.
!        RIGIDBODY=.TRUE.
!        NRBSITES=1 ! used in current GMIN
         CALL READF(STOCKMU)
         CALL READF(STOCKLAMBDA)
!        ALLOCATE(SITE(NRBSITES,3))
C
C    STOCKSPIN randomises the orientation of a Stockmayer cluster at any point in
C    an optimisation where a dipole vector becomes aligned with the z axis (which
C    make the phi angle for that dipole redundant).  STOCKZTOL is the amount by
C    which cos(theta) may differ from 1.0 for alignment to be recognised.
C    STOCKMAXSPIN is the maximum number of random orientations that will be attempted.
C
      ELSE IF (WORD.EQ.'STOCKSPIN') THEN
         STOCKSPIN = .TRUE.
         CALL READF(STOCKZTOL)
         CALL READI(STOCKMAXSPIN)
C
C  STOPDIST specifies an alternative stopping criterion based on displacement
C  between the first or last minimum and the furthest connected minimum.
C
      ELSE IF (WORD.EQ.'ST') THEN

         STOCKAAT =.TRUE.
         RBAAT    = .TRUE.
         CALL READF(STOCKMU)
         IF (NITEMS .GT. 2) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         NRBSITES = 1
         ALLOCATE(RBSITE(NRBSITES,3))
         NTSITES = NATOMS*NRBSITES/2
      ELSE IF (WORD.EQ.'STOPDISP') THEN
         CALL READF(STOPDISP)
         STOPDISPT=.TRUE.
C
C  In a CONNECT run, stop as soon as the initial minimum has a transition state
C  connection.
C
      ELSE IF (WORD.EQ.'STOPFIRST') THEN
         STOPFIRST=.TRUE.
C
C  SUMMARY n print a summary of the steps taken every n cycles  - default n=20   
C
      ELSE IF (WORD .EQ. 'SUMMARY') THEN
        IF (NITEMS.GT.1) CALL READI(NSUMMARY)
C
C  SYMCUT n RMS force below which symmetry subroutine is called - default 0.001
C
      ELSE IF (WORD .EQ. 'SYMCUT') THEN
        CALL READF(SYMCUT)
C
C  TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
C
C  Tagged particle - atom in question has mass increased by TAGFAC in symmetry.f and inertia.f
C
      ELSE IF (WORD.EQ.'TAG') THEN
         TAGT=.TRUE.
         NTAG=NTAG+1
         CALL READI(TAGNUM(NTAG))
         CALL READF(TAGFAC(TAGNUM(NTAG)))
      ELSE IF (WORD.EQ.'TANTYPE') THEN
         CALL READI(TANTYPE)
         GSTANTYPE = TANTYPE
C
C  Add a tetrahedral field to the potential of magnitude FTD.
C
      ELSE IF (WORD.EQ.'TD') THEN
         FIELDT=.TRUE.
         TDT=.TRUE.
         CALL READF(FTD)
C
C  TIMELIMIT - in seconds - OPTIM will stop if this limit is exceeded.
C
      ELSE IF (WORD.EQ.'TIMELIMIT') THEN
         CALL READF(TIMELIMIT)
C
C  TOLD n initial distance tolerance in symmetry subroutine     - default 0.0001
C
      ELSE IF (WORD .EQ. 'TOLD') THEN
        CALL READF(TOLD)
C
C  TOLE n initial tolerance for the difference in principal moments 
C         of inertia divided by the sum of the principal moments 
C         in symmetry subroutine                                - default 0.0001
C
      ELSE IF (WORD .EQ. 'TOLE') THEN
        CALL READF(TOLE)
C
C  Includes omega angles in the TWISTDIHE list.
C
      ELSE IF (WORD.EQ.'TOMEGA') THEN
         TOMEGAC=.TRUE.
      ELSE IF (WORD.EQ.'TOSIPOL') THEN
         TOSIPOL=.TRUE.
         CALL READF(ALPHAP)
         CALL READF(ALPHAM)
         CALL READF(DAMP)
         WRITE(*,'(A)') ' Polarizabilities:'
         WRITE(*,'(A,F12.8,A,F12.8)') ' alpha+=',ALPHAP,' alpha-=',ALPHAM
         WRITE(*,'(A,F12.8,A)') ' damping coefficent=',DAMP,' per bohr'
C
C  TRAD n sets the trust radius to n                            - default n=4       
C
      ELSE IF (WORD .EQ. 'TRAD') THEN
        CALL READF(TRAD)
C
C  TRAP is used for the trap potential in EYtrap coded by Ersin Yurtsever.
C
      ELSE IF (WORD .EQ. 'TRAP') THEN
        EYTRAPT=.TRUE.
        CALL READF(TRAPK)
        CALL READI(NTRAPPOW)
C
C  Includes sidechain angles in the TWISTDIHE list.
C
      ELSE IF (WORD.EQ.'TSIDECHAIN') THEN
         TSIDECHAIN=.TRUE.
C
C  Twist phi/psi dihedral angle nmode by xpert degrees before starting optimisation.
C
      ELSE IF (WORD.EQ.'TWISTDIHE') THEN
         TWISTDIHET=.TRUE.
         CALL READI(DMODE)
         CALL READF(DPERT)
C
C  TWISTTYPE specifies the type of twisting done to guess transition states in GUESSTS for CHARMM
C
      ELSE IF (WORD.EQ.'TWISTTYPE') THEN
         CALL READI(TWISTTYPE)
C
C  Double ended ts search.
C
      ELSE IF (WORD.EQ.'TWOENDS') THEN
         TWOENDS=.TRUE.
         IF (NITEMS.GT.1) CALL READF(FSTART)
         IF (NITEMS.GT.2) CALL READF(FINC)
         IF (NITEMS.GT.3) CALL READI(NTWO)
         IF (NITEMS.GT.4) CALL READF(RMSTWO)
         IF (NITEMS.GT.5) CALL READI(NTWOITER)
         IF (NITEMS.GT.6) CALL READF(TWOEVAL)
C
C  UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
C

      ELSE IF (WORD.EQ.'UNIAX') THEN
           UNIAXT = .TRUE.
           
      ELSE IF (WORD.EQ.'UNRES') THEN
         UNRST=.TRUE.
!CALL UNRESINIT
C CALPHAS AND THE SIDE CHAIN CENTROIDS ARE COUNTED AS ATOMS, BUT NOT THE PEPTIDE BOND CENTRES.
         NATOM=2*NRES
         NATOM=2*nres
         IF (NATOM /= NATOMS) THEN
            WRITE(*,'(A)') 'No. of atoms in "coords" conflicts with that deduced from unres part of odata'
            CALL FLUSH(6,ISTAT)
            STOP
         ENDIF
         NINTS=2*nres-5+2*nside ! jmc change this depending on how want to deal with non-capping glycines!
                                ! jmc NINTS was previously set in fetchz, but need either it or nvaru earlier (i.e. here)
                                ! so may as well set it when we first know nres and nside.

         IF (ENDHESS.AND.(.NOT.ENDNUMHESS)) THEN
            PRINT *,'**ERROR - to calculate normal mode frequencies for UNRES, please specify ENDNUMHESS keyword'
            CALL FLUSH(6,ISTAT)
            STOP
         ELSEIF ((DUMPPATH.OR.DUMPALLPATHS).AND.(.NOT.ENDHESS)) THEN
            PRINT *,'**ERROR - to calculate normal mode frequencies for UNRES, please specify ENDHESS and ENDNUMHESS keywords'
            CALL FLUSH(6,ISTAT)
            STOP
         ENDIF

c        DO J1=1,nres
C jmc c contains x,y,z for all the Calphas
c           UNRX(2*J1-1)=c(1,J1)
c           UNRY(2*J1-1)=c(2,J1)
c           UNRZ(2*J1-1)=c(3,J1)
C jmc then x,y,z for the side chain centroids
c           UNRX(2*J1)=c(1,J1+nres)
c           UNRY(2*J1)=c(2,J1+nres)
c           UNRZ(2*J1)=c(3,J1+nres)
c        ENDDO

C new read replaces random configuration coordinates with alternative from file coords
         CALL UNEWREAD(UNRX,UNRY,UNRZ,NATOMS,FILTH,FILTHSTR)
         DO J1=1,nres
            c(1,J1)=UNRX(2*J1-1)
            c(2,J1)=UNRY(2*J1-1)
            c(3,J1)=UNRZ(2*J1-1)
            c(1,J1+nres)=UNRX(2*J1)
            c(2,J1+nres)=UNRY(2*J1)
            c(3,J1+nres)=UNRZ(2*J1)
         ENDDO
         CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL CHAINBUILD
C JMC PUT COORDS IN STANDARD ORIENTATION (1ST ATOM AT 0,0,0 ETC...) INTO UNR ARRAY.  FIXES PROBLEM IN PATH,
C FOR CALCULATING THE STEP OFF THE TS
         DO J1=1,NRES
            UNRX(2*J1-1)=C(1,J1)
            UNRY(2*J1-1)=C(2,J1)
            UNRZ(2*J1-1)=C(3,J1)
            UNRX(2*J1)=C(1,J1+NRES)
            UNRY(2*J1)=C(2,J1+NRES)
            UNRZ(2*J1)=C(3,J1+NRES)
         ENDDO
         CALL UNRSETZSYMATMASS
         IF (FILTH.NE.0) THEN
            OPEN(UNIT=20,FILE='coords.read',STATUS='REPLACE')
            CLOSE(20)
         ENDIF
         ALLOCATE(UREFCOORD(3*NATOMS),UREFPPSANGLE(3*NATOMS))
         IF (TWISTDIHET.OR.PERTDIHET.OR.GUESSTST.OR.CALCDIHE) THEN
            CALL UNRSETDIHE
         ENDIF
         IF (TWISTDIHET) THEN
            CALL UNRSTWISTDIHE(UNRX,UNRY,UNRZ,DMODE,DPERT)
         ENDIF
         IF (PERTDIHET) THEN
            CALL UNRSPERTDIHE(UNRX,UNRY,UNRZ,CHPMIN,CHPMAX,CHNMIN,CHNMAX,ISEED)
         ENDIF
         IF (CALCDIHE) THEN
            CALL UNREADREF(NATOMS)
C jmc readref2 leaves reference coords in unres c and internal coord arrays, so replace with UNR{X,Y,Z} here.
            DO J1=1,nres
               c(1,J1)=UNRX(2*J1-1)
               c(2,J1)=UNRY(2*J1-1)
               c(3,J1)=UNRZ(2*J1-1)
               c(1,J1+nres)=UNRX(2*J1)
               c(2,J1+nres)=UNRY(2*J1)
               c(3,J1+nres)=UNRZ(2*J1)
            ENDDO
            CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
         END IF

         DO J1=1,NATOMS
            Q(3*(J1-1)+1)=UNRX(J1)
            Q(3*(J1-1)+2)=UNRY(J1)
            Q(3*(J1-1)+3)=UNRZ(J1)
         ENDDO
!
! USEDIAG enables the user to select DIAG or DIAG2 as the eigenvalue estimate in
! Rayleigh-Ritz routine secdiag. Default is currently one, but two may be better!
!
      ELSE IF (WORD.EQ.'USEDIAG') THEN
         CALL READI(NSECDIAG)

!
! USEEV allows the lowest NUSEEV eigenvalues and associated eigenvectors to be
! used in second-order searches with efol.f90.
!
      ELSE IF (WORD.EQ.'USEEV') THEN
         CALL READI(NUSEEV)
C
C  Number of BFGS updates before resetting, default=4
C
      ELSE IF (WORD.EQ.'UPDATES') THEN
         CALL READI(MUPDATE)
         IF (NITEMS.GT.2) CALL READI(XMUPDATE)
         IF (NITEMS.GT.3) CALL READI(MMUPDATE)
         IF (NITEMS.GT.4) CALL READI(NEBMUPDATE)
         IF (NITEMS.GT.5) CALL READI(INTMUPDATE)
         IF (NITEMS.GT.6) CALL READI(GSUPDATE)     
         IF (NITEMS.GT.7) CALL READI(GCUPDATE)     
C
C  VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
C
C
C  VALUES n print the Hessian eigenvalues every n cycles        - default n=20     
C
      ELSE IF (WORD .EQ. 'VALUES') THEN
        CALL READI(NVALUES)
C
C  VARIABLES - keyword at the end of the list of options after which
C           the general variables follow. NZERO is the number of zero
C  eigenvalues, default 0.
C
      ELSE IF (WORD.EQ.'VARIABLES') THEN
         VARIABLES=.TRUE.
         RETURN
C
C  VECTORS n prints the eigenvectors every n cycles             - default OFF
C
      ELSE IF (WORD .EQ. 'VECTORS') THEN
        VECTORST=.TRUE.
        CALL READI(NVECTORS)
C
C  WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
C
      ELSE IF (WORD == 'WARRANTY') THEN
          CALL WARRANTY
C
C  Welch parameters for Born-Meyer binary salt potentials.
C  These are A++, A--, A+- and rho, in order, followed by 
C
      ELSE IF (WORD.EQ.'WELCH') THEN
         WELCH=.TRUE.
         CALL READF(APP)  
         CALL READF(AMM) 
         CALL READF(APM)
         CALL READF(RHO) 
         CALL READF(XQP)
         CALL READF(XQM)
         CALL READF(ALPHAP)
         CALL READF(ALPHAM)
C
C  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C
C  YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
C
C
C  ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
C
      ELSE IF (WORD.EQ.'ZEROS') THEN
         CALL READI(NZERO)
         ! </kwd>
         ! }}}
         ! <end>
      ELSE
        CALL REPORT(' keyword> Unrecognized command '//WORD,.TRUE.)
        STOP
      ENDIF

      CALL FLUSH(6,ISTAT)
      GOTO 190

      RETURN
      END
      ! </end>
