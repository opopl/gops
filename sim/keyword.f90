      SUBROUTINE KEYWORD(PROG)

      INTEGER PROG

      SELECTCASE(PROG)
        CASE(1) ! GMIN {{{
!op226>=================================== 
Cop226> Declarations {{{

      USE COMMONS
      USE PORFUNCS

      IMPLICIT NONE

      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST, IX, J1, JP, NPCOUNT, NTYPEA, NPCALL, NDUMMY, INDEX, J2, J3
      LOGICAL CAT, YESNO, PERMFILE
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      DOUBLE PRECISION XX, EPSAB, EPSBB, SIGAB, SIGBB
      LOGICAL END, SKIPBL, CLEAR, ECHO
      CHARACTER WORD*16,PBC*3,WORD2*10
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      DOUBLE PRECISION EAMLJA0, EAMLJBETA, EAMLJZ0, DUMMY
      COMMON /EAMLJCOMM/ EAMLJA0, EAMLJBETA, EAMLJZ0
      DOUBLE PRECISION SLENGTH, EPS
      INTEGER NOK, NBAD
      COMMON /BSNEW/ SLENGTH, NOK, NBAD, EPS
      DOUBLE PRECISION EPS2, RAD, HEIGHT
      COMMON /CAPS/ EPS2, RAD, HEIGHT

      DOUBLE PRECISION    PMAX,PMIN,NMAX,NMIN,SIDESTEP
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

      DOUBLE PRECISION :: LPL, LPR

Cop226> End declarations </begin> }}}
!op226>=================================== 
Cop226> Initializations <init> {{{
      NPCOUNT=0
      NPCALL=0
      NSEED=0
      NS=0
      NSSTOP=0
      HIT=.FALSE.
      SAVEQ=.TRUE.
      NSAVE=5
      NSAVEINTE=0
      TFAC=1.0D0
      RESIZET=.FALSE.
      STEPOUT=.FALSE.
      SUPERSTEP=.FALSE.
      NSUPER=10
      SUPSTEP=1.1D0
      SACCRAT=0.5D0
      NSACCEPT=100
      EVSTEPT=.FALSE.
      NEVS=100
      CEIG=0.1D0
      NEVL=100
      NVECTORS=2
      TEMPS=0.8
      NRBSITES=0
      CHFREQ=1
      FTRANS=1
      FROT=1
      ALLOCATE(FIXSTEP(1),FIXTEMP(1),FIXBOTH(1),TEMP(1),ACCRAT(1),STEP(1),ASTEP(1),OSTEP(1),BLOCK(1),NT(1),NQ(1),EPREV(1),
     @         JUMPMOVE(1),JUMPINT(1),JDUMP(1),COORDS(3*NATOMS,1),COORDSO(3*NATOMS,1),VAT(NATOMS,1),VATO(NATOMS,1),
     @         JUMPTO(1),SHELLMOVES(1),PTGROUP(1),NSURFMOVES(1),NCORE(1))
      DO JP=1,1
         FIXSTEP(JP)=.FALSE.
         FIXTEMP(JP)=.FALSE.
         FIXBOTH(JP)=.FALSE.
         TEMP(JP)=0.3D0
         ACCRAT(JP)=0.5D0
         STEP(JP)=0.3D0
         ASTEP(JP)=0.3D0
         OSTEP(JP)=0.3D0
         BLOCK(JP)=0
         NT(JP)=0
         JUMPMOVE(JP)=.FALSE.
         JUMPINT(JP)=100
         JDUMP(JP)=.FALSE.
         SHELLMOVES(JP)=.FALSE.
         PTGROUP(JP)='    '
         NSURFMOVES(JP)=0
         NCORE(JP)=0
      ENDDO
      NEWJUMP=.FALSE.
      PNEWJUMP=0.2D0
      ECONV=0.02D0
      TABOOT=.FALSE.
      NTAB=10
      CUTOFF=1.0D6
      PCUTOFF=1.0D6
      FINALCUTOFF=1.0D6
      MYPOWER=5
      NEON=.FALSE.
      RGCL2=.FALSE.
      AXTELL=.FALSE.
      ZSTAR=0.0D0
      GROUND=.FALSE.
      ARGON=.FALSE.
      ARNO=.FALSE.
      STAR=.FALSE.
      PLUS=.FALSE.
      TWOPLUS=.FALSE.
      DIPOLE=.FALSE.
      DUMPT=.FALSE.
      TARGET=.FALSE.
      SORTT=.FALSE.
      NTARGETS=0
      MSORIGT=.FALSE.
      MSTRANST=.FALSE.
      FRAUSIT=.FALSE.
      ANGST=.FALSE.
      MORSET=.FALSE.
      LB2T=.FALSE.
      DZTEST=.FALSE.
      ZETT1=.FALSE.
      ZETT2=.FALSE.
      P46=.FALSE.
      G46=.FALSE.
      BLNT=.FALSE.
      DFTBT=.FALSE.
      SW=.FALSE.
      XMUL=1
      SCT=.FALSE.
      SQUEEZET=.FALSE.
      NVEC=0
C     SQUEEZER=5.0D0
C     SQUEEZED=0.95D0
      DEBUG=.FALSE.
      SEEDT=.FALSE.
      FREEZECORE=.TRUE.
      FREEZE=.FALSE.
      FREEZERES=.FALSE.
      FREEZEALL=.FALSE.
      UNFREEZERES =.FALSE.
      NFREEZE=0
      ALLOCATE(FROZEN(NATOMS))
C csw34> The FROZENRES array is bigger than needed
      ALLOCATE(FROZENRES(NATOMS))
      DO J1=1,NATOMS
         FROZEN(J1)=.FALSE.
         FROZENRES(J1)=.FALSE.
      ENDDO
      FREEZEGROUPTYPE='GT'
      FREEZEGROUPT=.FALSE.
C csw34> DONTMOVE defaults
      DONTMOVET=.FALSE.
      NDONTMOVE=0
      DONTMOVEREST=.FALSE.
      DONTMOVEALL=.FALSE.
      DOMOVEREST=.FALSE.
      DONTMOVEGROUPT=.FALSE.
      DONTMOVEGROUPTYPE='GT'
      ALLOCATE(DONTMOVE(NATOMS))
      ALLOCATE(DONTMOVERES(NATOMS))
      DO J1=1,NATOMS
         DONTMOVE(J1)=.FALSE.
         DONTMOVERES(J1)=.FALSE.
      ENDDO
C END of DONTMOVE defaults
      CHECKCHIRALITY=.TRUE.
      NOCISTRANS=.TRUE.
      NOCISTRANSRNA=.FALSE.
      NOCISTRANSDNA=.FALSE.
      MINOMEGA=150.D0
      UACHIRAL=.FALSE.
      SETCHIRAL=.FALSE.
      FIELDT=.FALSE.
      OHT=.FALSE.
      IHT=.FALSE.
      TDT=.FALSE.
      D5HT=.FALSE.
      CENT=.FALSE.
      CENTXY=.FALSE.
      SETCENT=.FALSE.
      CENTX=0.0D0
      CENTY=0.0D0
      CENTZ=0.0D0
      QUCENTRE=.FALSE.
      FIXCOM=.FALSE.
      FIH=0.0D0
      FTD=0.0D0
      FD5H=0.0D0
      TOLD=0.0001D0
      TOLE=0.0001D0
      CUTT=.FALSE.
      PERIODIC=.FALSE.
      PARAMONOVPBCX=.FALSE.
      PARAMONOVPBCY=.FALSE.
      PARAMONOVPBCZ=.FALSE.
      PARAMONOVCUTOFF=.FALSE.
      LJSITE=.FALSE.
      BLJSITE=.FALSE.
      LJSITECOORDST=.FALSE.
      LJSITEATTR=.FALSE.
      NRUNS=0
      PCUT=1.0D0
      RADIUS=0.0D0
      MAXIT=500
      MAXIT2=500
      EXPFAC=10.0D0
      EXPD=1.0D0
      CQMAX=1.0D-10
      BQMAX=1.0D-3
      RHO=6.0D0
      NACCEPT=50
      NORESET=.FALSE.
      TSALLIST=.FALSE.
      QTSALLIS=0.0D0
      PARALLELT=.FALSE.
      TOSI=.FALSE.
      WELCH=.FALSE.
      BINARY=.FALSE.
      SHIFTCUT=.FALSE.
      FAL=.FALSE.
      FNI=.FALSE.
!     AMBER=.FALSE.
      AMHT=.FALSE.
      NINT_AMH=1
      DPARAM=1.0D0
      FAKEWATER=.FALSE.
      AMCUT= .FALSE.
      MGBWATER=.FALSE.
      BIN=.FALSE.
      AMBERSEED= .FALSE.
      FIXT= .FALSE.
      FIX= .FALSE.
      CAP= .TRUE.
      WATERSTEP= .FALSE.
      QCUTOFF= 1.0D6
      RCUTOFF= 1.0D6
      REALQCUTOFF= 1.0D6
      REALRCUTOFF= 1.0D6
      RINGROTSCALE=0.0D0
      TRACKDATAT=.FALSE.
      PROGRESS=.FALSE.
      listupdate=20

      BLJCLUSTER=.FALSE.

      CHRMMT=.FALSE.
      ACESOLV=.FALSE.
      ACEUPSTEP=50
      CHRIGIDTRANST=.FALSE.
      CHRIGIDROTT=.FALSE.
      CHNMAX=0.0D0
      CHMDT=.FALSE.
      CHMDFREQ=HUGE(1)
      CURRENTIMP=0
      BOXT=.FALSE.
      SPHERET=.FALSE.
      RMST=.FALSE.
      NEWCONFT=.FALSE.
      INTMINT=.FALSE.
      DAESTAT=.FALSE.
      MAKEOLIGOT=.FALSE.
      MAKEOLIGOSTART=.FALSE.
      TRANSXYT=.FALSE.
      ROTZT=.FALSE.
      NREPEAT=0
      NFIXSEG=0

C  sf344> AMBER stuff
      AMBERT=.FALSE.
      AMCHNMAX=0.0D0
      AMCHPMAX=0.0D0
      MDSTEPT=.FALSE.
      DUMPSTRUCTURES=.FALSE.
C csw34> RANDOMSEED now works for CHARMM also!
      RANDOMSEEDT=.FALSE.
C csw34> Dumping structures after every quench
      DUMPQUT=.FALSE.
C csw34> Dumping structures after every step (before quenching)
      DUMPSTEPST=.FALSE.
C csw34> Local sampling within distance constraints 
      LOCALSAMPLET=.FALSE. 
      ABTHRESH=999.99
      ACTHRESH=999.99
C csw34> AMBER interaction energy logical
      A9INTET=.FALSE.
      INTERESTORE=.FALSE.
C csw34> set COLDFUSION flag to .FALSE.
      COLDFUSION=.FALSE.
C
C  sf344> for specifying movable atoms and ligand rotating steptaking moves
C
      MOVABLEATOMST=.FALSE.
      LIGMOVET=.FALSE.
      LIGROTSCALE=0.0D0
      LIGCARTSTEP=0.0D0
      LIGTRANSSTEP=0.0D0
      LIGMOVEFREQ=1

C
C  csw34> rotamer move stuff
C
      ROTAMERT=.FALSE.
C
C  csw34> some defaults (just in case)
C
      ROTMAXCHANGE=1
      ROTPSELECT=0.2
      ROTOCCUW=0.004
      ROTCENTRE=1
      ROTCUTOFF=999.99
C
C  csw34> atom group rotation moves
C
      GROUPROTT=.FALSE.
      NGROUPS=0
      GROUPROTFREQ=1
      GROUPOFFSET=0
C
      NOPHIPSIT=.FALSE.
      OMEGAT=.FALSE.

      OSASAT=.FALSE.
      RPRO=1.4D0
      ODIHET=.FALSE.
      ORGYT=.FALSE.
      OEINTT=.FALSE.
      MON1(1:2)=1
      MON2(1:2)=1

      BSMIN=.FALSE.
      RKMIN=.FALSE.
      PERMDIST=.FALSE.
      PERMOPT=.FALSE.

      GAMMA=1.0D0
      TUNNELT=.FALSE.
      
      TWOD=.FALSE.
      COMPRESST=.FALSE.

      MUPDATE=4
      DGUESS=0.1D0
      BFGS=.FALSE.
      LBFGST=.TRUE.
      CONJG=.FALSE.
      TNT=.FALSE.
      TOLB=0.1D0
      DBRENTT=.FALSE.
      GUIDECUT=0.0001D0
      CPMD=.FALSE.
      DL_POLY=.FALSE.
      EFAC=0.0D0
      EAMP=0.01D0
      FIXD=.FALSE.
      NHSMOVE=1
      T12FAC=1.1D0
      RENORM=.FALSE.
      NRENORM=10
      NRENSTUCK=20
      XMOVERENORM=6.0
      TRENORM=1.0D0
      PACHECO=.FALSE.
      EAMLJT=.FALSE.
      PBGLUET=.FALSE.
      EAMALT=.FALSE.
      ALGLUET=.FALSE.
      MGGLUET=.FALSE.
      GUPTAT=.FALSE.
      FST=.FALSE.
      WENZEL=.FALSE.
      RESTART=.FALSE.
      NEWRESTART=.FALSE.
      NRELAX=0
      NMSBSAVE=0
      AVOID=.FALSE.
      AVOIDDIST=1.0D0
      MAXSAVE=10
      NHSRESTART=0
      MAXBFGS=0.4D0

      CAPSID=.FALSE.
      STRANDT=.FALSE.
      PAHT=.FALSE.
      TIP=.FALSE.
      QUADT=.FALSE.
      STOCKT=.FALSE.
      LJCOULT=.FALSE.
      COULN=0
      COULQ=0.0D0
      COULSWAP = 0.0D0
      COULTEMP = 0.0D0
      GAYBERNET=.FALSE.
      PARAMONOVT=.FALSE.
      ELLIPSOIDT=.FALSE.
      PYGPERIODICT=.FALSE.
      LJCAPSIDT=.FALSE.
      PYBINARYT=.FALSE.
      MULTISITEPYT=.FALSE.
      PYOVERLAPTHRESH=1.0D0
      LJSITE=.FALSE.
      SWAPMOVEST=.FALSE.
      STICKYT=.FALSE.
      RIGID=.FALSE.
      TIPID=4
      HEIGHT=0.5D0
      OTPT=.FALSE.
      LJMFT=.FALSE.
      Q4T=.FALSE.

!     DC430 >
      DBPT        = .FALSE.
      DBPTDT      = .FALSE.
      DBLPYT      = .FALSE.
      DMBLMT      = .FALSE.
      EFIELDT     = .FALSE.
      GAYBERNEDCT = .FALSE.
      GBDT        = .FALSE.
      GBDPT       = .FALSE.
      GEMT        = .FALSE.
      LINRODT     = .FALSE.
      LWOTPT      = .FALSE.
      MMRSDPT     = .FALSE.
      MSGBT       = .FALSE. 
      MSPYGT      = .FALSE.
      MSTBINT     = .FALSE.
      MSSTOCKT    = .FALSE.
      MULTPAHAT   = .FALSE.
      NCAPT       = .FALSE.
      NPAHT       = .FALSE.
      NTIPT       = .FALSE.
      PAHAT       = .FALSE.
      PYGT        = .FALSE.
      PYGDPT      = .FALSE.
      STOCKAAT    = .FALSE.
      TDHDT       = .FALSE.
      WATERDCT    = .FALSE.
      WATERKZT    = .FALSE.
!|gd351>
      PATCHY = .FALSE.
      ASAOOS = .FALSE.
!<gd351|


       
      THRESHOLDT=.FALSE.
      BSWL=.FALSE.
      BSPT=.FALSE.
      BSPTQMAX=1.0D100
      BSPTQMIN=-1.0D100
      BSPTRESTART=.FALSE.
      HISTSMOOTH=.FALSE.
      NSpline=1
      EPSSPHERE=0.0D0
      FIXBIN=.FALSE.
      QUENCHFRQ=1
      NQUENCH=0

      DECAY=.FALSE.
      DECAYPARAM=0.0D0
      COOP=.FALSE.
      NCOOP=5
      COOPCUT=1.0D0

      UNSTRING='UNDEFINED'
      BOXSIZE=20.D0
      SPHERERAD=20.D0
      NCHENCALLS=0
      NATBT=.FALSE.
      MAXERISE=1.0D-10
      MAXEFALL=-HUGE(1.0D0)
      SYMMETRIZE=.FALSE.
      SYMMETRIZECSM=.FALSE.
      NSYMINTERVAL=10
      SYMTOL1=0.1D0
      SYMTOL2=0.1D0
      SYMTOL3=0.1D0
      SYMTOL4=0.1D0
      SYMTOL5=0.1D0
      NSYMQMAX=20
      MATDIFF=0.1D0
      DISTFAC=0.0D0
      ARMA=0.4D0
      ARMB=0.4D0

      BINSTRUCTURES=.FALSE.
      TETHER=.FALSE.
      EQUIL=0
      PTMC=.FALSE.
      VISITPROP=.FALSE.
      HWINDOWS=1

      FixedEndMoveT=.FALSE.
      PIVOTP=0.0D0
      SIDECHAINP=0.0D0

      DIFFRACTT=.FALSE.
      THOMSONT=.FALSE.
      GAUSST=.FALSE.
      MPIT=.FALSE.
      DUMPINT=1000 ! default is to dump a restart file every 1000 cycles of mc.f
      RESTORET=.FALSE.
      DUMPFILE=''
      INTEDUMPFILE=''
      MOVESHELLT=.FALSE.
      SHELLMOVEMAX=0
      SHELLPROB=0.0D0
      COREFRAC=0.0D0
      MYSDT=.FALSE.
      TSTAR=-1.0D0
      PRTFRQ=1
      BSPTDUMPFRQ=100000
      QUENCHDOS=.FALSE.
      QDT=.FALSE.
      QD2T=.FALSE.
      MULLERBROWNT=.FALSE.
      QDLIMIT=-1
      CHARMMENERGIES=.FALSE.
      FIXDIHEFLAG=.FALSE.
      JMT=.FALSE.
      PROJIT=.FALSE.
      PROJIHT=.FALSE.
      COLDFUSIONLIMIT=-1.0D6
      MODEL1T=.FALSE.

      VGW=.FALSE.              ! VGW PARAMETERS
      LJSIGMA=1.0D0
      LJEPSILON=1.0D0
      TAUMAX=5.0D0
      TAUMAXFULL=7.0D0
      CPFACTORSG=3.0D0
      CPFACTORFG=1.0D0
      CPS=1
      CPF=1
      VGWTOL=0.0001D0
      ACKLANDT=.FALSE.
      ACKLANDID=5
      STEEREDMINT=.FALSE.
      DF1T=.FALSE.
      PULLT=.FALSE.
      CSMT=.FALSE.
      CSMGUIDET=.FALSE.
      CSMEPS=1.0D-6
      CSMSTEPS=1
      CSMQUENCHES=1
      CSMMAXIT=0
      CHECKMARKOVT=.FALSE.
Cop226> End initializations </init>  }}}
!op226>=================================== 
!op226> 1) Open file 'data' => unit number 5
!op226> 2) Read in a single WORD from this file
!op226> 3) Check whether WORD equals any of keywords in the below ELSEIF
!op226>         loop
!op226>=================================== 
!op226> <read> {{{ 

      CALL OPENF(DATA_FH,"<","data")

190   CALL INPUT(END)
      IF (.NOT. END) THEN
        CALL READU(WORD)
      ENDIF

!op226> IF (END .OR. WORD .EQ. 'STOP') THEN {{{
      IF (END .OR. WORD .EQ. 'STOP') THEN
        RETURN
      ENDIF
!op226>}}} 

      IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'
     &                          .OR. WORD .EQ. '\\') THEN 
         GOTO 190

!op226> </read> }}}
!op226>=================================== 
C
C  Remaining documented keywords should be in alphabetical order.
C
Cop226> ============================================== 
!op226> Starting to check whether WORD fits any keyword from a 
!op226> really looong list below {{{ 
Cop226> ============================================== 
Cop226> Keywords: 2D ACCEPTRATIO ACE ACKLAND ALGLUE ...
Cop226> AMBER (Commented)
Cop226> ============================================== 
Cop226> {{{
!op226> <kwd>
      ELSE IF (WORD.EQ.'2D') THEN
         TWOD=.TRUE.

      ELSE IF (WORD.EQ.'ACCEPTRATIO') THEN
         IF (NITEMS-1.GT.NPAR) THEN
            WRITE(MYUNIT,'(A)') 'Number of acceptance ratios exceeds NPAR - quit'
            STOP
         ENDIF
         DO J1=1,NITEMS-1
            CALL READF(ACCRAT(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            DO J1=NITEMS,NPAR
               ACCRAT(J1)=ACCRAT(1)
            ENDDO
         ENDIF
C
C  bs360: ACE is to be used together with CHARMM and the ACE solvent model,
C  it makes sure that the Born radii are regularly updated
C
      ELSE IF (WORD.EQ.'ACE') THEN
          ACESOLV=.TRUE.
          IF (NITEMS.GT.1) CALL READI(ACEUPSTEP)
C
C  Ackland embedded atom metal potentials.
C
      ELSE IF (WORD.EQ.'ACKLAND') THEN
         ACKLANDT=.TRUE.
         CALL READI(ACKLANDID) ! default is 5 = Fe

      ELSE IF (WORD.EQ.'ALGLUE') THEN
          ALGLUET=.TRUE.

      ELSE IF (WORD.EQ.'CISTRANS') THEN
          NOCISTRANS=.FALSE.

!     ELSE IF (WORD.EQ.'AMBER') THEN
!        AMBER=.TRUE.
!        CALL APARAMS
!        CALL AREAD
!        NATOMS=ATOMS
!        DO J1=1,NATOMS
!           COORDS(3*(J1-1)+1,1)=x(J1)
!           COORDS(3*(J1-1)+2,1)=y(J1)
!           COORDS(3*(J1-1)+3,1)=z(J1)
!        ENDDO
!        t=0
!        ang=0
!        imp=0
!        count=0
C
C Dump Amber9 energy components at every step
C
      ELSE IF (WORD.EQ.'AMBERENERGIES') THEN
         AMBERENERGIEST=.TRUE.

      ELSE IF (WORD.EQ.'AMCHPMAX') THEN
         CALL READF(AMCHPMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'AMCHPMAX=  ',AMCHPMAX

      ELSE IF (WORD.EQ.'AMCHPMIN') THEN
         CALL READF(AMCHPMIN)
         WRITE(MYUNIT,'(A,F14.10)') 'AMCHPMIN=  ',AMCHPMIN

      ELSE IF (WORD.EQ.'AMCHNMAX') THEN
         CALL READF(AMCHNMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'AMCHNMAX=  ',AMCHNMAX

      ELSE IF (WORD.EQ.'AMCHNMIN') THEN
         CALL READF(AMCHNMIN)
         WRITE(MYUNIT,'(A,F14.10)') 'AMCHNMIN=  ',AMCHNMIN

      ELSE IF (WORD.EQ.'AMH') THEN
         AMHT=.TRUE.
         WRITE(MYUNIT,'(A)')'USING AMH ENERGIES FORCES'
         WRITE(MYUNIT,'(A)')'CALCULATE ENERGY AND FORCE TABLES  '
         CALL WALESAMH_INITIAL 
      ELSE IF (WORD.EQ.'NINT_AMH') THEN
         CALL READI(NINT_AMH)
         WRITE(MYUNIT,*)'NINT_AMH' , NINT_AMH
      ELSE IF (WORD.EQ.'HARM_AMH') THEN
         CALL READI(HARM_AMH)
         WRITE(MYUNIT,*)'HARM_AMH' , HARM_AMH
      ELSE IF (WORD.EQ.'ANGSTROM') THEN
         ANGST=.TRUE.

      ELSE IF (WORD.EQ.'ARGON') THEN
         ARGON=.TRUE.
  
      ELSE IF (WORD.EQ.'ARM') THEN
         ARMT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(ARMA)
         IF (NITEMS.GT.2) CALL READF(ARMB)

      ELSE IF (WORD.EQ.'ARNO') THEN
         ARNO=.TRUE.
C
C  Specify resetting if the latest structure gets too close to minima saved
C  in MSBCOORDS. Use bipartite matching and closest approach distance AVOIDDIST.
C  Maximum number of saved structures is specified by MAXSAVE.
C 
      ELSE IF (WORD.EQ.'AVOID') THEN
         AVOID=.TRUE.
         IF (NITEMS.GT.1) CALL READF(AVOIDDIST)
         IF (NITEMS.GT.2) CALL READI(MAXSAVE)
         IF (.NOT.ALLOCATED(MSBCOORDS)) THEN
            ALLOCATE(MSBCOORDS(3*NATOMS,MAXSAVE))
         ELSE
            WRITE(MYUNIT,*) 'reallocating MSBCOORDS'
            DEALLOCATE(MSBCOORDS)
            ALLOCATE(MSBCOORDS(3*NATOMS,MAXSAVE))
         ENDIF
         IF (.NOT.ALLOCATED(MSBE)) THEN
            ALLOCATE(MSBE(MAXSAVE))
         ELSE
            WRITE(MYUNIT,*) 'reallocating MSBE'
            DEALLOCATE(MSBE)
            ALLOCATE(MSBE(MAXSAVE))
         ENDIF

      ELSE IF (WORD.EQ.'AXTELL') THEN
         AXTELL=.TRUE.
         CALL READF(ZSTAR)
C }}}
Cop226 ============================================== 
Cop226 BASIN, BFGS, BHPT, BINARY, BINSTRUCTURES
Cop226 ==============================================
Cop226> {{{
      ELSE IF ((WORD.EQ.'BASIN').OR.(WORD.EQ.'SLOPPYCONV')) THEN
         IF (NITEMS.GT.1) CALL READF(BQMAX)

      ELSE IF (WORD.EQ.'BFGS') THEN
         BFGS=.TRUE.
C
C PT basin-hopping.
C
      ELSE IF (WORD.EQ.'BHPT') THEN
         CALL READF(PTTMIN)
         CALL READF(PTTMAX)
         CALL READF(EXCHPROB)

      ELSE IF (WORD.EQ.'BINARY') THEN
         BINARY=.TRUE.
         CALL READI(NTYPEA)
         CALL READF(EPSAB)
         CALL READF(EPSBB)
         CALL READF(SIGAB)
         CALL READF(SIGBB)
C
C Saves every n'th structure to the file with corresponding bin label
C
      ELSE IF (WORD.EQ.'BINSTRUCTURES') THEN
         BINSTRUCTURES=.TRUE.
         CALL READI(SAVENTH)
C BLJCLUSTER {{{
      ELSE IF (WORD.EQ.'BLJCLUSTER') THEN
         BLJCLUSTER=.TRUE.
         CALL READI(NTYPEA)
         CALL READF(EPSAB)
         CALL READF(EPSBB)
         CALL READF(SIGAB)
         CALL READF(SIGBB)
         CALL READF(CUTOFF)
C }}}
C BLN {{{
      ELSE IF (WORD.EQ.'BLN') THEN
         BLNT=.TRUE.
         CALL READF(RK_R)
         CALL READF(RK_THETA)
         ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     &            LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS))
         OPEN(UNIT=1,FILE='BLNsequence',STATUS='OLD')
         READ(1,*) DUMMYCH
         READ(1,*) LJREPBB, LJATTBB
         READ(1,*) LJREPLL, LJATTLL
         READ(1,*) LJREPNN, LJATTNN
         READ(1,*) DUMMYCH
         READ(1,*) DUMMYCH
         READ(1,*) HABLN, HBBLN, HCBLN, HDBLN
         READ(1,*) EABLN, EBBLN, ECBLN, EDBLN
         READ(1,*) TABLN, TBBLN, TCBLN, TDBLN
         DO J1=1,NATOMS-1
            READ(1,'(A1)',ADVANCE='NO') BEADLETTER(J1)
         ENDDO
         READ(1,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
         DO J1=1,NATOMS-3
            READ(1,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
         ENDDO
         CLOSE(1)
         WRITE(MYUNIT,'(A,I8,A)') 'BLN sequence of ',NATOMS,' beads read:'
         WRITE(MYUNIT,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A,I8,A)') 'BLN dihedral types:'
         WRITE(MYUNIT,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
         WRITE(MYUNIT,'(A)') ' '
         WRITE(MYUNIT,'(A,2F15.5)') 'B-B LJ coefficients: ',LJREPBB, LJATTBB
         WRITE(MYUNIT,'(A,2F15.5)') 'L-L LJ coefficients: ',LJREPLL, LJATTLL
         WRITE(MYUNIT,'(A,2F15.5)') 'N-N LJ coefficients: ',LJREPNN, LJATTNN
         WRITE(MYUNIT,'(A,4F15.5)') 'Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
         WRITE(MYUNIT,'(A,4F15.5)') 'Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
         WRITE(MYUNIT,'(A,4F15.5)') 'Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
         call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, 
     &                       HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, NATOMS)
C        call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
C    &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, NATOMS) 
C End BLN }}}
      ELSE IF (WORD.EQ.'BSMIN') THEN
         BSMIN=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GMAX)
         IF (NITEMS.GT.2) CALL READF(EPS)
         WRITE(MYUNIT,'(A,2L5)') 'BSMIN branch RKMIN,BSMIN=',RKMIN,BSMIN
C
C Basin-sampling. {{{
C
      ELSE IF (WORD.EQ.'BSPT') THEN
         BSPT=.TRUE.
         CALL READF(HISTMIN)
         CALL READF(HISTMAX)
         CALL READF(PTEMIN)
         CALL READF(PTEMAX)
         CALL READF(PTTMIN)
         CALL READF(PTTMAX)
         CALL READF(EXCHPROB)
         CALL READF(NEQUIL)  ! equilibration only
         CALL READF(PTSTEPS) ! PT following equilibration
         CALL READF(NQUENCH) ! combined PT and quench steps
         CALL READI(NENRPER)
         CALL READI(HBINS)
         CALL READI(QUENCHFRQ)
C }}}
C
C Frequency of dumping the Visits.his.n files
C 
      ELSE IF (WORD.EQ.'BSPTDUMPFRQ') THEN
         CALL READI(BSPTDUMPFRQ)
C
C BSPT restriction on quench energy.
C
      ELSE IF (WORD.EQ.'BSPTQRANGE') THEN
         CALL READF(BSPTQMIN)
         CALL READF(BSPTQMAX)
C
C Restart from Visits and bsptrestart files.
C
      ELSE IF (WORD.EQ.'BSPTRESTART') THEN
         BSPTRESTART=.TRUE.
C
C WL Basin-sampling. {{{
C
      ELSE IF (WORD.EQ.'BSWL') THEN
         BSWL=.TRUE.
         CALL READF(HISTMIN)
         CALL READF(HISTMAX)
         CALL READF(HISTFAC)
         CALL READI(HBINS)
         CALL READF(HISTFACMUL)
         CALL READI(TargetWL)
         CALL READF(HPERCENT)
         ALLOCATE(HDIST(HBINS),HWEIGHT(HBINS),HISTVALS(HBINS),LHISTVALS(HBINS),IGNOREBIN(HBINS))
         DO J1=1,HBINS
            HISTVALS(J1)=0
            LHISTVALS(J1)=0
            HWEIGHT(J1)=1.0D0
            HDIST(J1)=0.0D0
C           DO J2=1,HBINS
C              HTRANS(J1,J2)=1.0D0 ! transition matrix
C           ENDDO
         ENDDO
C
C  During the BSWL run HDIST contains the sum of the distances found for minima in each bin. The
C  average is saved in hist.new.
C
         DO J1=1,HBINS
            HDIST(J1)=HDIST(J1)*HISTVALS(J1)
         ENDDO
C }}}
      ELSE IF (WORD.EQ.'CAPSID') THEN
         CAPSID=.TRUE.
         RIGID=.TRUE.
         CALL READF(RHO)
         CALL READF(EPS2)
         CALL READF(RAD)
         IF (NITEMS.GT.4) CALL READF(HEIGHT)
C
C  The six reference site positions per capped pentagon. These need to
C  be multiplied by RAD, including the repulsive site!
C
         NRBSITES=6
         ALLOCATE(SITE(NRBSITES,3))
         SITE(1,1)=1.0D0
         SITE(1,2)=0.0D0
         SITE(1,3)=0.0D0
         SITE(2,1)=((-1.0D0 + Sqrt(5.0D0)))/4.0D0
         SITE(2,2)=(Sqrt((5.0D0 + Sqrt(5.0D0))/2.0D0))/2.0D0
         SITE(2,3)=0.0D0
         SITE(3,1)=((-1.0D0 - Sqrt(5.0D0)))/4.0D0
         SITE(3,2)=(Sqrt((5.0D0 - Sqrt(5.0D0))/2.0D0))/2.0D0
         SITE(3,3)=0.0D0
         SITE(4,1)=((-1 - Sqrt(5.0D0)))/4.0D0
         SITE(4,2)=-(Sqrt((5.0D0 - Sqrt(5.0D0))/2.))/2.0D0
         SITE(4,3)=0.0D0
         SITE(5,1)=((-1 + Sqrt(5.0D0)))/4.0D0
         SITE(5,2)=-(Sqrt((5.0D0 + Sqrt(5.0D0))/2.))/2.0D0
         SITE(5,3)=0.0D0
         SITE(6,1)=0.0D0
         SITE(6,2)=0.0D0
         SITE(6,3)=HEIGHT
Cop226> }}}
Cop226> ============================================== 
Cop226> CENTRE CENTREXY SETCENTRE QUCENTRE CG
Cop226> ============================================== 
Cop226> {{{
      ELSE IF (WORD.EQ.'CENTRE') THEN
         CENT=.TRUE.
       
C csw34> When using the implicit membrane potential IMM1, we want to let
C the molecule move in and out of the membrane (z direction), but fix it
C in x and y.
       
      ELSE IF (WORD.EQ.'CENTREXY') THEN
         CENT=.TRUE.
         CENTXY=.TRUE.

C csw34> SETCENTRE moves the centre of coordinates to the specified
C location before the initial quench is done.

      ELSE IF (WORD.EQ.'SETCENTRE') THEN
         SETCENT=.TRUE.
         IF (NITEMS.EQ.2) THEN
            CALL READF(CENTX) 
         ELSE IF (NITEMS.EQ.3) THEN
            CALL READF(CENTX)
            CALL READF(CENTY)
         ELSE IF (NITEMS.EQ.4) THEN
            CALL READF(CENTX)
            CALL READF(CENTY)
            CALL READF(CENTZ)
         ENDIF

C       csw34> QUCENTRE moves the centre of coordinates to (0,0,0)
C       before each step is taken. 
      ELSE IF (WORD.EQ.'QUCENTRE') THEN
         QUCENTRE=.TRUE.
C
C  Conjugate gradient optimisation instead of default LBFGS
C
      ELSE IF (WORD.EQ.'CG') THEN
         LBFGST=.FALSE.
         CONJG=.TRUE.
Cop226> }}}
Cop226> ============================================== 
C 
C sf344> Start of AMBER-related keywords 
C
Cop226> ============================================== 
Cop226> AMBERMDSTEPS AMBER9 MODEL1 MOVABLEATOMS
Cop226> ============================================== 
Cop226> {{{ 
      ELSE IF (WORD.EQ.'AMBERMDSTEPS') THEN
        MDSTEPT = .TRUE.
!op226>=================================== 
!op226> AMBER 9
!op226>=================================== 
!op226>{{{ 
      ELSE IF (WORD.EQ.'AMBER9') THEN
        AMBERT=.TRUE.
C
C csw34> if residues are frozen with FREEZERES, call the amber routine
C to fill the FROZEN array correctly (in amberinterface.f) 
C
        IF (PERMDIST) THEN
          IF(NPERMSIZE(1).EQ.NATOMS) THEN
           PRINT '(A)','keyword> ERROR - PERMDIST is specfied for AMBER, but there is no perm.allow file present'
           STOP
          ENDIF
        ENDIF

! sf344> file open unit used to conflict with AMBER's IO units (mdin opened with unit = 5),

!               call amberinterface(natom,1,trim(adjustl(inpcrd)),MYUNIT)
        IF(NITEMS==2) then
         inpcrd='coords.inpcrd'
         CALL READA(amberstr)
         WRITE(MYUNIT,'(A)') 'keywords> input coordinates for AMBER9 system will be read from ', trim(adjustl(amberstr))
               call amberinterface(natom,2,inpcrd,MYUNIT)
               CALL amber_readcoords(amberstr)

        ELSE IF(NITEMS==3) then
         CALL READA(amberstr)
         CALL READA(amberstr1)
         WRITE(MYUNIT,'(A)') 'keywords> input coordinates for AMBER9 system will be read from ', trim(adjustl(amberstr)),
     &                              'type: ', trim(adjustl(amberstr1))
          IF(trim(adjustl(amberstr1)).EQ.'inpcrd') then
               inpcrd=amberstr
               call amberinterface(natom,2,inpcrd,MYUNIT)
           WRITE(MYUNIT,'(A)') 'keywords> reading AMBER inpcrd coordinate format'
          ELSE
           WRITE(MYUNIT,'(A)') 'keywords> ERROR - no other types defined currently than inpcrd'
           STOP
          END IF
        END IF
               IF(.NOT.ALLOCATED(COORDS1)) ALLOCATE(COORDS1(3*NATOM))
               IF(.NOT.ALLOCATED(MOVABLEATOMLIST)) ALLOCATE(MOVABLEATOMLIST(NATOMS))
               IF(ALLOCATED(COORDS)) DEALLOCATE(COORDS)
               ALLOCATE(COORDS(3*NATOM,NPAR))
               NATOMS = NATOM
             DO J1=1,NPAR
               COORDS(:,J1) = COORDS1(:)
             END DO
!                natoms = natom

!                WRITE(*,*) 'sf344> keywords.f, natoms = ', natoms
!        call amopen(9,AMBERPRMTOP,'O','F','R')
        IF (FREEZERES) THEN 
           IF (UNFREEZERES) THEN
              CALL A9RESTOATOM(FROZENRES,FROZEN,NFREEZE,.TRUE.)
           ELSE
              CALL A9RESTOATOM(FROZENRES,FROZEN,NFREEZE,.FALSE.)
           ENDIF
        ENDIF

        IF(DONTMOVEREST) THEN
           IF (DOMOVEREST) THEN
              CALL A9RESTOATOM(DONTMOVERES,DONTMOVE,NDONTMOVE,.TRUE.)
           ELSE
              CALL A9RESTOATOM(DONTMOVERES,DONTMOVE,NDONTMOVE,.FALSE.)
           ENDIF
        ENDIF

Cop226 }}}
!op226>=================================== 
!op226> MODEL1
!op226>=================================== 
!op226>{{{ 
      ELSE IF(WORD.EQ.'MODEL1') THEN
         MODEL1T=.TRUE.
         CALL READF(ME1)
         CALL READF(ME2)
         CALL READF(ME3)
         CALL READF(MSTART)
         CALL READF(MFINISH)
         CALL READF(MBSTART1)
         CALL READF(MBFINISH1)
         CALL READF(MBSTART2)
         CALL READF(MBFINISH2)
         CALL READF(MBHEIGHT1)
         CALL READF(MBHEIGHT2)

! sf344> keyword for taking constrained random steps of just one part of the molecule between quenches
!        should be useful for systems with ligands docked in a protein, this would be correspond to some
!        partially flexible docking scheme
!op226>}}}
Cop226> ============================================== 
Cop226> MOVABLEATOMS
Cop226> ============================================== 
!op226>{{{ 
      ELSE IF(WORD.EQ.'MOVABLEATOMS') THEN
         MOVABLEATOMST=.TRUE.
         nmovableatoms=0
! csw34> for steered minimisation, need to read in the atoms used to define the points force acts on/from
!        group A is the ligand, group B is the protein. They start at -1 as groups delimited by A and B
! csw34> also, for LOCALSAMPLEing, need to define a group of atoms in the ligand (A) and two groups of atoms (B and C) 
!        that are used to calculate distances to the ligand
         natomsina=-1
         natomsinb=-1
         natomsinc=-1
!        readswitch controls which part of movableatoms we're reading in i.e. 
!        M=whole ligand, A=group A, B=group B, C=group C
         readswitch='M'
           WRITE(MYUNIT,'(A)') ' keyword> list of movable atoms will be read from file <movableatoms>'
           OPEN(unit=222,file='movableatoms',status='old')
                 do
                   read (unit=222,iostat=iostatus,fmt='(A6)') check1
                     if (iostatus<0) then
                        close(222)
                        exit
                     else if (TRIM(ADJUSTL(check1)).EQ.'A') then
                        readswitch='A'
                     else if (TRIM(ADJUSTL(check1)).EQ.'B') then
                        readswitch='B'
                     else if (TRIM(ADJUSTL(check1)).EQ.'C') then
                        readswitch='C'
                     end if
                     if (readswitch.EQ.'M') then 
                        nmovableatoms = nmovableatoms + 1
                     else if (readswitch.EQ.'A') then
                        natomsina = natomsina + 1
                     else if (readswitch.EQ.'B') then
                        natomsinb = natomsinb + 1
                     else if (readswitch.EQ.'C') then
                        natomsinc = natomsinc + 1
                     endif
                 end do
! setup arrays for movableatoms
           if(.not.allocated(movableatomlist)) ALLOCATE(movableatomlist(nmovableatoms))
           if(.not.allocated(movableatomlistlogical)) ALLOCATE(movableatomlistlogical(natoms))
           movableatomlistlogical(:)=.false.
! setup arrays for steered minimisation groups
           if (natomsina.gt.0) then 
              if(.not.allocated(atomsinalist)) ALLOCATE(atomsinalist(natomsina))
              if(.not.allocated(atomsinalistlogical)) ALLOCATE(atomsinalistlogical(natoms))
              atomsinalistlogical(:)=.false.
           endif
           if (natomsinb.gt.0) then 
              if(.not.allocated(atomsinblist)) ALLOCATE(atomsinblist(natomsinb))
              if(.not.allocated(atomsinblistlogical)) ALLOCATE(atomsinblistlogical(natoms))
              atomsinblistlogical(:)=.false.
           endif
           if (natomsinc.gt.0) then 
              if(.not.allocated(atomsinclist)) ALLOCATE(atomsinclist(natomsinc))
              if(.not.allocated(atomsinclistlogical)) ALLOCATE(atomsinclistlogical(natoms))
              atomsinclistlogical(:)=.false.
           endif
! now open movableatoms for the second time to actually read in the atom indicies
           OPEN(unit=222,file='movableatoms',status='old')
              do i=1,nmovableatoms
                 read (unit=222,iostat=iostatus,fmt='(I6)') movableatomlist(i) 
              end do
! if groups for steered minimisation are specified, read them in
           if (natomsina.gt.0) then 
              read (unit=222,iostat=iostatus,fmt='(A6)') check1
                 do i=1,natomsina
                    read (unit=222,iostat=iostatus,fmt='(I6)') atomsinalist(i) 
                 end do
           endif
           if (natomsinb.gt.0) then 
              read (unit=222,iostat=iostatus,fmt='(A6)') check1
                 do i=1,natomsinb
                    read (unit=222,iostat=iostatus,fmt='(I6)') atomsinblist(i) 
                 end do
           endif
           if (natomsinc.gt.0) then 
              read (unit=222,iostat=iostatus,fmt='(A6)') check1
                 do i=1,natomsinc
                    read (unit=222,iostat=iostatus,fmt='(I6)') atomsinclist(i) 
                 end do
           endif
! now need to set the logicals to .true. for all atoms in movableatom groups
              do i=1,nmovableatoms
                 movableatomlistlogical(movableatomlist(i))=.true.
              end do
           if (natomsina.gt.0) then
              do i=1,natomsina
              end do
           endif
           if (natomsinb.gt.0) then
              do i=1,natomsinb
                 atomsinblistlogical(atomsinblist(i))=.true.
              end do
           endif
           if (natomsinb.gt.0) then
              do i=1,natomsinc
                 atomsinclistlogical(atomsinclist(i))=.true.
              end do
           endif
! write some output for the user
           WRITE(MYUNIT,'(A15,I6,A46)') ' keyword> ligand defined with ',nmovableatoms,' atoms'
           IF(natomsina.GT.0) WRITE(MYUNIT,'(A10,I6,A64)') ' keyword> ',natomsina,' atoms in group A'
           IF(natomsinb.GT.0) WRITE(MYUNIT,'(A10,I6,A63)') ' keyword> ',natomsinb,' atoms in group B'
           IF(natomsinc.GT.0) WRITE(MYUNIT,'(A10,I6,A63)') ' keyword> ',natomsinc,' atoms in group C'
! sf344> keyword for taking moves between quenches in which one part of the molecule is moved
! csw34> updated on 29th September 2009 to allow local cartesian moves 
!op226>}}} 
!op226>=================================== 
!op226> LIGMOVE
!op226>=================================== 
!op226>{{{ 
      ELSE IF(WORD.EQ.'LIGMOVE') THEN
        LIGMOVET=.TRUE.
        IF (NITEMS.GT.1) THEN
           CALL READF(ligrotscale)
        ENDIF
        IF (NITEMS.GT.2) THEN
           CALL READF(ligcartstep)
        ENDIF
        IF (NITEMS.GT.3) THEN
           CALL READF(ligtransstep)
        ENDIF
        IF (NITEMS.GT.4) THEN
           CALL READI(ligmovefreq)
! csw34> to prevent an arithmetic exception (divide by 0), we need to
! check that the frequency of ligand moves is > 0. Otherwise, disable
! the moves.
           IF(ligmovefreq.EQ.0) THEN 
              LIGMOVET=.FALSE.
              WRITE(MYUNIT,'(A)') ' keyword> WARNING: frequency of LIGMOVE moves set to 0 - moves DISABLED!'
           ENDIF
        ENDIF

! csw34> some user info about ligand moves 
        IF (ligrotscale.gt.0) THEN
           WRITE(MYUNIT,'(A)') ' keyword> one part of the molecule (ligand) will be randomly rotated during MC steptaking moves'
           WRITE(MYUNIT,'(A,G8.3)') ' keyword> ligand rotations scaled will be scaled by ',ligrotscale 
        ENDIF
        IF (ligcartstep.gt.0) THEN 
           WRITE(MYUNIT,'(A,G8.3,A)') ' keyword> cartesian perturbations of up to ',ligcartstep,' will be applied to the ligand' 
        ENDIF 

      ELSE IF (WORD.eq.'DUMPSTRUCTURES') THEN
        DUMPSTRUCTURES=.TRUE.
        WRITE(MYUNIT,'(A)') ' keywords> Final structures will be dumped in different formats (.rst, .xyz, .pdb)'
      ELSE IF (WORD.eq.'RANDOMSEED') THEN
        WRITE(MYUNIT,'(A)') ' keywords> The random number generator for the random steptaking moves will be seeded with system time'
        RANDOMSEEDT=.TRUE.
      ELSE IF (WORD.EQ.'CHANGEACCEPT') THEN
         CALL READI(IX)
         NACCEPT=IX

! csw34> Dumping after every quench in AMBER rst and pdb format
      ELSE IF (WORD.eq.'DUMPQU') THEN
         DUMPQUT=.TRUE.
! csw34> Dumping after every step (before quenching) in AMBER rst and pdb format
      ELSE IF (WORD.eq.'DUMPSTEPS') THEN
         DUMPSTEPST=.TRUE.
! csw34> AMBER interaction energy using script AMBGMINintE.sh
      ELSE IF (WORD.eq.'A9INTE') THEN
         YESNO=.FALSE.
         INQUIRE(FILE='AMBGMINintE.sh',EXIST=YESNO)
         IF (YESNO) THEN
            A9INTET=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> The interaction enthalpy to the specified residue will be calculated after each quench'
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: NEED AMBGMINintE.sh SCRIPT TO USE A9INTE - see www-wales.ch.cam.ac.uk/software'
            STOP
         ENDIF

! csw34> FREEZEGROUP 
C FREEZEGROUP lets you freeze all atoms within or beyond a radius
      ELSE IF (WORD.EQ.'FREEZEGROUP') THEN
         FREEZE=.TRUE.
         FREEZEGROUPT=.TRUE.
         CALL READI(GROUPCENTRE)
         CALL READF(GROUPRADIUS)
         IF(NITEMS.GT.3) CALL READA(FREEZEGROUPTYPE)

! csw34> Rotamer moves 
      ELSE IF (WORD.EQ.'ROTAMER') THEN
         YESNO=.FALSE.
         INQUIRE(FILE='PdbRotamerSearch',EXIST=YESNO)
         IF (YESNO) THEN
            ROTAMERT=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> AMBER rotamer moves enabled'
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: NEED PdbRotamerSearch TO USE ROTAMER - see www-wales.ch.cam.ac.uk/software'
            STOP
         ENDIF
         CALL READI(ROTMAXCHANGE)
         CALL READF(ROTPSELECT)
         CALL READF(ROTOCCUW)
         IF (NITEMS.GT.4) CALL READI(ROTCENTRE)
         IF (NITEMS.GT.5) CALL READF(ROTCUTOFF)

! csw34> Group rotation moves 
      ELSE IF (WORD.EQ.'GROUPROTATION') THEN
! csw34> Check the group file is present
         YESNO=.FALSE.
         INQUIRE(FILE='atomgroups',EXIST=YESNO)
         IF (YESNO) THEN
            GROUPROTT=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> AMBER group rotation moves enabled'
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: atom groups must be defined in atomgroups file'
            STOP
         ENDIF
         IF (NITEMS.GT.1) CALL READI(GROUPROTFREQ)
! csw34> if the frequency is 0, we need to disable the moves to present
! a divide by 0!
         IF(GROUPROTFREQ.EQ.0) THEN 
            GROUPROTT=.FALSE.
            WRITE(MYUNIT,'(A)') ' keyword> WARNING: frequency of GROUPROTATION moves set to 0 - moves DISABLED!'
         ENDIF
         IF (NITEMS.GT.2) CALL READI(GROUPOFFSET)
! csw34> Figure out how many atom groups have been defined
         NGROUPS=0
         OPEN(UNIT=222,FILE='atomgroups',status='old')
         DO
            READ(222,*,IOSTAT=iostatus) CHECK1
            IF (iostatus<0) THEN
            CLOSE(222)
            EXIT
            ELSE IF (TRIM(ADJUSTL(check1)).EQ.'GROUP') then
               NGROUPS=NGROUPS+1
            ENDIF
         END DO        
         CLOSE(222)
! csw34> Allocate atom group info arrays appropriately
         ALLOCATE(ATOMGROUPNAMES(NGROUPS))
         ALLOCATE(ATOMGROUPAXIS(NGROUPS,2))
         ALLOCATE(ATOMGROUPPSELECT(NGROUPS))
         ALLOCATE(ATOMGROUPSCALING(NGROUPS))
         ALLOCATE(ATOMGROUPS(NGROUPS,NATOMS))
! csw34> Set safe defaults
         ATOMGROUPS(:,:)=.FALSE.
         ATOMGROUPNAMES(:)='EMPTY'
         ATOMGROUPAXIS(:,:)=0
         ATOMGROUPSCALING(:)=1.0D0
         ATOMGROUPPSELECT(:)=1.0D0
! csw34> Read in group info
! Here is an example entry:
! GROUP OME 6 5 4 1.0
! 1
! 2
! 3
! 4
! This says that group OME is to be rotated about the bond from atom 6->5.
! There are 4 atoms in the OME group. Rotations of -pi->+pi are to be scaled by 1.0. 
! Finally, the group members are specified one per line
         OPEN(UNIT=222,FILE='atomgroups',status='unknown')
         WRITE(MYUNIT,*) 'keyword> Reading in atom groups for GROUPROTATION'
         IF(GROUPOFFSET.NE.0) WRITE(MYUNIT,*) 'keyword> Group atom numbering offset by ',GROUPOFFSET
         DO J1=1,NGROUPS
            READ(222,*) CHECK1,ATOMGROUPNAMES(J1),AXIS1,AXIS2,GROUPSIZE,ATOMGROUPSCALING(J1),
     &                       ATOMGROUPPSELECT(J1) 
            ATOMGROUPAXIS(J1,1)=AXIS1+GROUPOFFSET
            ATOMGROUPAXIS(J1,2)=AXIS2+GROUPOFFSET
            CALL FLUSH(MYUNIT)
            IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') THEN
               DO J2=1,GROUPSIZE
                  READ(222,*) GROUPATOM
                  IF(GROUPOFFSET.GT.0) GROUPATOM=GROUPATOM+GROUPOFFSET
                  ATOMGROUPS(J1,GROUPATOM)=.TRUE. 
               END DO 
            ELSE
               WRITE(MYUNIT,'(A)') ' keyword: ERROR! Group file not formatted correctly!'
               STOP
            ENDIF
            WRITE(MYUNIT,'(3A)') '<GROUP ',TRIM(ADJUSTL(ATOMGROUPNAMES(J1))),'>'
            WRITE(MYUNIT,'(A,I3)') 'Index: ',J1
            WRITE(MYUNIT,'(A,I4)') 'Size: ',GROUPSIZE
            WRITE(MYUNIT,'(A,2I5)') 'Atoms defining axis: ',ATOMGROUPAXIS(J1,1),ATOMGROUPAXIS(J1,2)
            WRITE(MYUNIT,'(A,F4.2)') 'Rotation scaling: ',ATOMGROUPSCALING(J1)
            WRITE(MYUNIT,'(A,F4.2)') 'Selection probablity: ',ATOMGROUPPSELECT(J1)
            WRITE(MYUNIT,'(A)') 'Members:'
            DO J2=1,NATOMS
               IF(ATOMGROUPS(J1,J2)) WRITE(MYUNIT,*) J2
            ENDDO
         ENDDO
         CLOSE(222)
C }}}
!op226>}}} 
!op226>=================================== 
C
C Start of CHARMM-related keywords, including options for MC moves and order parameter specifications.
C
!op226>=================================== 
!op226>{{{ 
!op226>=================================== 
!op226> CHARMM 
!op226>=================================== 
!op226>{{{ 
      ELSE IF (WORD.EQ.'CHARMM') THEN
         CHRMMT=.TRUE.
         IF (MXATMS.EQ.0) THEN
            WRITE(MYUNIT,'(A)') 'keyword> ERROR *** MXATMS is zero'
            STOP
         ENDIF

         IF (PERMDIST) THEN
            IF(NPERMSIZE(1).EQ.NATOMS) THEN
            WRITE(MYUNIT,'(A)') 'keyword> ERROR - PERMDIST is specfied for CHARMM, but there is no perm.allow file present'
            STOP
            ENDIF
         ENDIF

         ALLOCATE(CHX(MXATMS),CHY(MXATMS),CHZ(MXATMS),CHMASS(MXATMS))

         CHX(1)=13.13d13 ! this way we will tell CHARMM to save it's coords into CH. arrays; otherwise it will
                         ! use input.crd only which is the default now
         CALL CHSETUP(CHX,CHY,CHZ,CHMASS,NATOM,TOPFILE,PARFILE)
         CALL CHSETZSYMATMASS
         CALL CHALLOCATE(NATOMS)
C jmc49>         CALL CHSETDIHE
C Moved the above call to CHSETDIHE to below the "IF (FREEZERES) CALL CHRESTOATOM(FROZENRES,FROZEN)" below.
C Necessary for the case of freezing residues when the MC moves are made in internal coordinates, and 
C inconsequential otherwise.

C csw34> If we're freezing RESIDUES, fill the FROZEN array with atoms that are
C within the frozen residues. This requires info from a CHARMM routine
C and so will be done in charmmgmin.src

         IF (FREEZERES) THEN 
            IF (UNFREEZERES) THEN
               CALL CHRESTOATOM(FROZENRES,FROZEN,NFREEZE,.TRUE.)
            ELSE
               CALL A9RESTOATOM(FROZENRES,FROZEN,NFREEZE,.FALSE.)
            ENDIF
         ENDIF

         IF(DONTMOVEREST) THEN
            IF (DOMOVEREST) THEN
               CALL CHRESTOATOM(DONTMOVERES,DONTMOVE,NDONTMOVE,.TRUE.)
            ELSE
               CALL CHRESTOATOM(DONTMOVERES,DONTMOVE,NDONTMOVE,.FALSE.)
            ENDIF
         ENDIF
         CALL CHSETDIHE
         
         IF (NATOM /= NATOMS) THEN
            WRITE(MYUNIT,'(A)') 'No. of atoms in "input.crd" and file specified in CHARMM part of data conflict'
            WRITE(MYUNIT,'(A,2I8)') 'NATOM,NATOMS=',NATOM, NATOMS
            CALL EXIT(10)
         ENDIF
         IF (MPIT) THEN
            DO J1=1,NATOMS
               COORDS(3*(J1-1)+1,MYNODE+1)=CHX(J1)
               COORDS(3*(J1-1)+2,MYNODE+1)=CHY(J1)
               COORDS(3*(J1-1)+3,MYNODE+1)=CHZ(J1)
            ENDDO
         ELSE
            DO J1=1,NATOMS
               COORDS(3*(J1-1)+1,1)=CHX(J1)
               COORDS(3*(J1-1)+2,1)=CHY(J1)
               COORDS(3*(J1-1)+3,1)=CHZ(J1)
            ENDDO
         ENDIF
         DEALLOCATE(CHX,CHY,CHZ,CHMASS)
         IF (INTMINT) CALL GETNINT(NINTS)  ! DJW - this is OK because CHARMM is the last keyword!
         IF (ODIHET) CALL READREF ! likewise, this is OK and in fact essential because CHARMM is the last keyword!
C }}}
!op226>=================================== 
!op226> CHARMMENERGIES
!op226>=================================== 
!op226>{{{ 
C
C Dump charmm energy components at every step
C
      ELSE IF (WORD.EQ.'CHARMMENERGIES') THEN
         CHARMMENERGIES=.TRUE.
!op226>}}}
!op226>=================================== 
!op226> CHARMMTYPE
!op226>=================================== 
!op226>{{{ 
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
C Sanity check for the energy in COORDSO.
C
!op226>}}} 
!op226>=================================== 
!op226> CHECKMARKOV
!op226>=================================== 
!op226>{{{ 
      ELSE IF (WORD.EQ.'CHECKMARKOV') THEN
         CHECKMARKOVT=.TRUE.
C
C MD for CHARMM for the generation of new structures.
C
!op226>}}} 
!op226>=================================== 
!op226> CHMD
!op226>=================================== 
!op226>{{{ 
      ELSE IF (WORD.EQ.'CHMD') THEN
         CHMDT=.TRUE.
         INQUIRE(FILE='chmd.par',EXIST=YESNO)
         IF (YESNO) THEN
            OPEN(99,file='chmd.par')
            READ(99,'(A)') CHMDPAR
            CLOSE(99)
         ENDIF
         CALL READI(CHMDFREQ) 
C
!op226>}}} 
!op226>=================================== 
!op226> CHPMAX CHPMIN CHNMAX CHNMIN NOPHIPSI TOMEGA CHFREQ
!op226>=================================== 
!op226>{{{ 
      ELSE IF (WORD.EQ.'CHPMAX') THEN
         CALL READF(CHPMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'CHPMAX=  ',CHPMAX

      ELSE IF (WORD.EQ.'CHPMIN') THEN
         CALL READF(CHPMIN)
         WRITE(MYUNIT,'(A,F14.10)') 'CHPMIN=  ',CHPMIN

      ELSE IF (WORD.EQ.'CHNMAX') THEN
         CALL READF(CHNMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'CHNMAX=  ',CHNMAX

      ELSE IF (WORD.EQ.'CHNMIN') THEN
         CALL READF(CHNMIN)
         WRITE(MYUNIT,'(A,F14.10)') 'CHNMIN=  ',CHNMIN

      ELSE IF (WORD.EQ.'NOPHIPSI') THEN
         NOPHIPSIT=.TRUE.
         WRITE(MYUNIT,'(A)') 'NOPHIPSIT set : only sidechain dihedrals will be twisted'

      ELSE IF (WORD.EQ.'TOMEGA') THEN
         OMEGAT=.TRUE.
         WRITE(MYUNIT,'(A)') 'TOMEGA set : peptide bonds will be twisted along with all other dihedrals'

      ELSE IF (WORD.EQ.'CHFREQ') THEN
         CALL READI(CHFREQ)
         WRITE(MYUNIT,'(A,I4,A)') 'Every ',CHFREQ,' steps TAKESTEPCH is called'
!op226>}}} 
!op226>=================================== 
!op226> CHRIGIDTRANS CHRIGIDROT FIXEDEND OSASA ODIHE OEINT
!op226>=================================== 
!op226>{{{ 

      ELSE IF (WORD.EQ.'CHRIGIDTRANS') THEN
         CHRIGIDTRANST=.TRUE.
         CALL READF(PTRANS)
         CALL READF(TRANSMAX)
         CALL READI(FTRANS)
         IF(NITEMS.GT.4) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BOX') BOXT=.TRUE.
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SPHERE') SPHERET=.TRUE.
         IF (BOXT.AND.NITEMS.GT.5) CALL READF(BOXSIZE)
         IF (SPHERET.AND.NITEMS.GT.5) CALL READF(SPHERERAD)

      ELSE IF (WORD.EQ.'CHRIGIDROT') THEN
         CHRIGIDROTT=.TRUE.
         CALL READF(PROT)
         CALL READF(ROTMAX)
         CALL READI(FROT)

      ELSE IF (WORD.EQ.'FIXEDEND') THEN
         FixedEndMoveT = .TRUE.
         IF (NITEMS.GT.1) CALL READF(PIVOTP)
         IF (NITEMS.GT.2) CALL READF(SIDECHAINP)

      ELSE IF (WORD.EQ.'OSASA') THEN
         OSASAT=.TRUE.
         CALL READF(RPRO)
         WRITE(MYUNIT,'(A)') 'OSASA set: solvent accessible surface area order parameter will be calculated'
         WRITE(MYUNIT,'(A,F3.1)') 'using probe radius ',RPRO

      ELSE IF (WORD.EQ.'ODIHE') THEN
         ODIHET=.TRUE.
         WRITE(MYUNIT,'(A)') 'ODIHE set: dihedral-angle order parameter will be calculated'
         WRITE(MYUNIT,'(A)') 'using the reference structure supplied in ref.crd'
 
      ELSE IF (WORD.EQ.'OEINT') THEN
         OEINTT=.TRUE.
         CALL READI(MON1(1))
         CALL READI(MON1(2))
         CALL READI(MON2(1))
         CALL READI(MON2(2))
         WRITE(MYUNIT,'(A)') 'OEINTT set: interaction energy between 2 peptides will be used as an order parameter'
!op226>}}} 
!op226>=================================== 
!op226> ORGYR NORANDOM PERMDIHE 
!op226>=================================== 
!op226>{{{ 

      ELSE IF (WORD.EQ.'ORGYR') THEN
         ORGYT=.TRUE.
         WRITE(MYUNIT,'(A)') 'ORGYT set: radius of gyration will be calculated as an order parameter'

C     ELSE IF (WORD.EQ.'NORANDOM') THEN
C        NORANDOM=.TRUE.
C        IF (NITEMS.GT.1) CALL READF(RANDOMCUTOFF)

C     ELSE IF (WORD.EQ.'PERMDIHE') THEN
C        PERMDIHET=.TRUE.
C        DO J1=1,NITEMS-1
C           CALL READI(NDUM)
C           PERMDIHE(J1)=NDUM
C        ENDDO
C        NPERMDIHE=NITEMS-1
C        DO J1=1,NITEMS-1
C           print *,'PERMDIHE',PERMDIHE(J1)
C        ENDDO
C
!op226>}}} 
C  End of CHARMM block.  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!op226>}}} 
C
C
C  Sometimes have to modify the cold fusion limit when using high electric fields
C
      ELSE IF (WORD.EQ.'COLDFUSION') THEN
         CALL READF(COLDFUSIONLIMIT)

      ELSE IF (WORD.EQ.'COMPRESS') THEN
         COMPRESST=.TRUE.
         CALL READF(COMP)
C
C  Alternative correlated moves: NCOOP nearest neighbours of a randomly selected atom
C  all move by the same amount. NCOOP default is 5.
C
      ELSE IF (WORD.EQ.'COOPMOVE') THEN
         COOP=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCOOP)
         IF (NITEMS.GT.2) CALL READF(COOPCUT)

      ELSE IF (WORD.EQ.'CPMD') THEN
         CPMD=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READA(SYS)
         ELSE
            WRITE(MYUNIT,'(A)') ' ERROR - no CPMD system specified'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 12
            ENDIF
         ENDDO
12       CONTINUE
C
C  Calculate the continuous symmetry measure from Pinsky, Dryzun, Casanova, Alemany and Avnir,
C  J. Comp. Chem., 29, 2712, 2008.
C
      ELSE IF (WORD.EQ.'CSM') THEN
         CSMT=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READA(CSMGP)
            CALL MYUPCASE(CSMGP)
         ELSE
            PRINT '(A)','keyword> ERROR - point group must be specified for CMS keyword'
            STOP
         ENDIF
         IF (NITEMS.GT.2) CALL READF(CSMEPS)
         IF (NITEMS.GT.3) THEN
            CSMGUIDET=.TRUE.
            CALL READA(CSMGUIDEGP)
            CALL MYUPCASE(CSMGUIDEGP)
         ENDIF

         IF (.NOT.PERMDIST) THEN ! set permdist if not done already
         PERMDIST=.TRUE.
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
!        ALLOCATE(NPERMSIZE(NATOMS),PERMGROUP(NATOMS),NSWAP(NATOMS),SWAP1(NATOMS,2),SWAP2(NATOMS,2))
         ALLOCATE(NPERMSIZE(3*NATOMS),PERMGROUP(3*NATOMS),NSETS(3*NATOMS),SETS(NATOMS,3))
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
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
               IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),
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
            ENDIF
         ELSE
            NSETS(1:NATOMS)=0
            NPERMGROUP=1 ! all atoms can be permuted - default
            NPERMSIZE(1)=NATOMS ! all atoms can be permuted - default
            DO J1=1,NATOMS
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         WRITE(MYUNIT,'(A,I6)') ' keyword> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            PRINT '(A,3(I6,A))',' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
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
      ELSE IF (WORD.EQ.'CUTOFF') THEN
         CUTT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(CUTOFF)
         FINALCUTOFF=CUTOFF
         IF (NITEMS.GT.2) CALL READF(FINALCUTOFF)

C
C NOT DOCUMENTED - INTENTIONAL
C
      ELSE IF (WORD.EQ.'D5H') THEN
         FIELDT=.TRUE.
         D5HT=.TRUE.
         CALL READF(XX)
         FD5H=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            EXPFAC=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(EXPD)
         ENDIF

      ELSE IF (WORD.EQ.'DBRENT') THEN
         DBRENTT=.TRUE.

      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.
C
C  Correlated random moves, in the sense that the magnitude of the step
C  decays exponentially as DECAYPARAM from a randomly chosen atom.
C
      ELSE IF (WORD.EQ.'DECAY') THEN
         DECAY=.TRUE.
         CALL READF(DECAYPARAM)
C
C  2D binary potential - Daan Frenkel #1
C

      ELSE IF (WORD.eq.'DF1') THEN
         DF1T=.TRUE.

      ELSE IF (WORD.EQ.'DFTB') THEN
         DFTBT=.TRUE.
C
C  Initial guess for diagonal matrix elements in LBFGS.
C
      ELSE IF (WORD.EQ.'DGUESS') THEN
         CALL READF(DGUESS)

      ELSE IF (WORD.EQ.'DIELEC') THEN
         CALL READF(XX)
         DPARAM=XX
         WRITE(MYUNIT,'(A,F9.5)') ' Dielectric constant = ',DPARAM
C
C  NOT DOCUMENTED - INTENTIONAL
C
      ELSE IF (WORD.EQ.'DIFFRACT') THEN
         DIFFRACTT=.TRUE.

      ELSE IF (WORD.EQ.'DIPOLES') THEN
         DIPOLE=.TRUE.
C
C  NOT DOCUMENTED - INTENTIONAL
C
      ELSE IF (WORD.EQ.'DL_POLY') THEN
         DL_POLY=.TRUE.
         CALL READI(NATOMS)
!
! csw34> DONTMOVE is an extension of the FREEZE idea. Atoms are replaced
! after the STEP, but their gradients are NOT set to 0, hence they may
! move during minimisation.
!
      ELSE IF (WORD.EQ.'DONTMOVE') THEN
         DONTMOVET=.TRUE.
         DO J1=1,NITEMS-1
            NDONTMOVE=NDONTMOVE+1
            CALL READI(NDUMMY)
            DONTMOVE(NDUMMY)=.TRUE.
         ENDDO
         
! csw34> DONTMOVEGROUP is analagous to FREEZEGROUP 
      ELSE IF (WORD.EQ.'DONTMOVEGROUP') THEN
         DONTMOVET=.TRUE.
         DONTMOVEGROUPT=.TRUE.
         CALL READI(DONTMOVECENTRE)
         CALL READF(DONTMOVERADIUS)
         IF(NITEMS.GT.3) CALL READA(DONTMOVEGROUPTYPE)
         
      ELSE IF (WORD.EQ.'DONTMOVEALL') THEN
         DONTMOVET=.TRUE.
         DONTMOVEALL=.TRUE.
         NDONTMOVE=NATOMS
         DO J1=1,NATOMS
            DONTMOVE(J1)=.TRUE.
            DONTMOVERES(J1)=.TRUE.
         ENDDO
C csw34
C Things are then moved using the DOMOVE and DOMOVERES keywords
C This is only a valid keyword if DONTMOVEALL is also specified
      ELSE IF ((WORD.EQ.'DOMOVE').AND.DONTMOVEALL) THEN
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            DONTMOVE(NDUMMY)=.FALSE.
            NDONTMOVE=NDONTMOVE-1
         ENDDO

      ELSE IF (WORD.EQ.'DONTMOVERES') THEN
         DONTMOVET=.TRUE.
         DONTMOVEREST=.TRUE.
C The FROZENRES array is then filled with the residue number from the
C data file
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            DONTMOVERES(NDUMMY)=.TRUE.
         ENDDO
         
         ELSEIF ((WORD.EQ.'DOMOVERES').AND.DONTMOVEALL) THEN
         DOMOVEREST=.TRUE.
C Set the right parts of the DONTMOVERES array to FALSE
         DO J1=1,NITEMS-1 
            CALL READI(NDUMMY)
            DONTMOVERES(NDUMMY)=.FALSE.
         ENDDO
         DONTMOVEREST=.TRUE.
         
      ELSE IF (WORD.EQ.'DUMP') THEN
         DUMPT=.TRUE.

      ELSE IF (WORD.EQ.'DUMPINT') THEN
         CALL READI(DUMPINT)

      ELSE IF (WORD.EQ.'DZUGUTOV') THEN
         DZTEST=.TRUE.
         CALL READF(DZP1)
         CALL READF(DZP2)
         CALL READF(DZP3)
         CALL READF(DZP4)
         CALL READF(DZP5)
         CALL READF(DZP6)
         CALL READF(DZP7)

      ELSE IF (WORD.EQ.'EAMAL') THEN
         EAMALT=.TRUE.

      ELSE IF (WORD.EQ.'EAMLJ') THEN
         EAMLJT=.TRUE.
         CALL READF(EAMLJA0)
         CALL READF(EAMLJBETA)
         CALL READF(EAMLJZ0)

      ELSE IF (WORD.EQ.'EDIFF') THEN
         CALL READF(ECONV)
C
C Accumulation of thermodynamic statistics starting after Equil steps, 
C calculated thermodynamic properties is dumped every DumpEveryNthQuench quench.
C
      ELSE IF (WORD.EQ.'EQUILIBRATION') THEN
         CALL READI(EQUIL)
         CALL READI(DumpEveryNthQuench)
C
C  Steps using transition state search-type moves. Obsolete.
C  NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'EVSTEP') THEN
         EVSTEPT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NEVS)
         IF (NITEMS.GT.2) CALL READF(CEIG)
         IF (NITEMS.GT.3) CALL READI(NEVL)
         IF (NITEMS.GT.4) CALL READI(NVECTORS)
C
C  NOT DOCUMENTED - INTENTIONAL
C
      ELSE IF (WORD.EQ.'EXPFAC') THEN
         CALL READF(EFAC)
         IF (NITEMS.GT.2) CALL READF(EAMP)

      ELSE IF (WORD.EQ.'FAKEWATER') THEN
         FAKEWATER=.TRUE.
         WRITE (MYUNIT,'(A)') '**********************************************************'
         WRITE (MYUNIT,'(A)') '* DISTANCE DEPENDENT DIELECTRIC BEING USED - FAKE WATER! *'
         WRITE (MYUNIT,'(A)') '**********************************************************'

      ELSE IF (WORD.EQ.'FAL') THEN
         FAL=.TRUE.
      
      ELSE IF (WORD.EQ.'FIXBOTH') THEN
         IF (NITEMS.EQ.1) THEN
            FIXBOTH(1)=.TRUE.
            IF (NPAR.GT.1) THEN
               DO J1=2,NPAR
                  FIXBOTH(J1)=.TRUE.
               ENDDO
            ENDIF
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX) 
               FIXBOTH(IX)=.TRUE.
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'FIXCOM') THEN
         FIXCOM=.TRUE.
         ALLOCATE(MASSES(NATOMS))
C
C  Take hard sphere type moves.
C  T12FAC is the fraction of the first collision time to be used in HSMOVE
C
      ELSE IF (WORD.EQ.'FIXD') THEN
         FIXD=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(NHSMOVE)
         ENDIF
         IF (NITEMS.GT.2) THEN
            CALL READF(T12FAC)
         ENDIF

      ELSE IF (WORD.EQ.'FIXSTEP') THEN
         IF (NITEMS.EQ.1) THEN
            FIXSTEP(1)=.TRUE.
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX) 
               FIXSTEP(IX)=.TRUE.
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'FIXTEMP') THEN
         IF (NITEMS.EQ.1) THEN
            FIXTEMP(1)=.TRUE.
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX) 
               FIXTEMP(IX)=.TRUE.
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'FNI') THEN
         FNI=.TRUE.

      ELSE IF (WORD.EQ.'FRAUSI') THEN
         FRAUSIT=.TRUE.
C
C  Frozen atoms.
C
      ELSE IF (WORD.EQ.'FREEZE') THEN
         FREEZE=.TRUE.
         DO J1=1,NITEMS-1
            NFREEZE=NFREEZE+1
            CALL READI(NDUMMY)
            FROZEN(NDUMMY)=.TRUE.
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
            CALL READI(NDUMMY)
            FROZENRES(NDUMMY)=.TRUE.
         ENDDO
C Finally, the frozen residue numbers are converted into frozen atom
C numbers. This is also forcefield dependant and must be done when we
C know which forcefield to use (i.e. in the CHARMM block above)

C csw34
C Freezing EVERYTHING and then permitting small parts to move
C This is useful for large system to prevent the data file getting silly
      ELSEIF (WORD.EQ.'FREEZEALL') THEN
         FREEZE=.TRUE.
         FREEZEALL=.TRUE.
         NFREEZE=NATOMS
         DO J1=1,NATOMS
            FROZEN(J1)=.TRUE.
            FROZENRES(J1)=.TRUE.
         ENDDO
C csw34
C Things are then UNFROZEN using the UNFREEZE and UNFREEZERES keywords
C This is only a valid keyword if FREEZEALL is also specified
      ELSEIF ((WORD.EQ.'UNFREEZE').AND.FREEZEALL) THEN
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            FROZEN(NDUMMY)=.FALSE.
            NFREEZE=NFREEZE-1
         ENDDO

      ELSEIF ((WORD.EQ.'UNFREEZERES').AND.FREEZEALL) THEN
         UNFREEZERES=.TRUE.
C Set the right parts of the FROZENRES array to FALSE
         DO J1=1,NITEMS-1 
            CALL READI(NDUMMY)
            FROZENRES(NDUMMY)=.FALSE.
         ENDDO
         FREEZERES=.TRUE.
C The FREEZERES routines for AMBER and CHARMM do the rest :)
      
      
C Finnis-Sinclair potential coded by James Elliott
C
      ELSE IF (WORD.EQ.'FS') THEN
         FST=.TRUE.
         CALL READI(GATOM)

      ELSE IF (WORD.EQ.'G46') THEN
         G46=.TRUE.
         BLNT=.TRUE.
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'GAMMA') THEN
         CALL READF(GAMMA)
         TUNNELT=.TRUE.
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'GAUSS') THEN
         GAUSST=.TRUE.
         CALL READI(GMODES) ! number of nodes

      ELSE IF (WORD.EQ.'GROUND') THEN
         GROUND=.TRUE.

      ELSE IF (WORD.EQ.'GUIDE') THEN
         CALL READF(GUIDECUT)

      ELSE IF (WORD.EQ.'GUPTA') THEN
         GUPTAT=.TRUE.
         CALL READI(GATOM)
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'HISTSMOOTH') THEN
         CALL READI(NSpline)
C
C Parameters of the temperature range on which to calculate thermodynamic properties in Basin Sampling
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'HISTTEMP') THEN
         CALL READF(MinimalTemperature)
         CALL READF(MaximalTemperature)
         CALL READI(NTempPoints)
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'IH') THEN
         FIELDT=.TRUE.
         IHT=.TRUE.
         CALL READF(XX)
         FIH=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            EXPFAC=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(XX)
            EXPD=XX
         ENDIF

      ELSE IF (WORD.EQ.'INTMIN') THEN
         INTMINT=.TRUE.
C        IF (NITEMS.GT.1) THEN
C           CALL READF(IMINCUT)
C        ENDIF
 
      ELSE IF (WORD.EQ.'JM') THEN
         JMT=.TRUE.

      ELSE IF (WORD.EQ.'JUMPMOVE') THEN
         CALL READI(IX)
         JUMPMOVE(IX)=.TRUE.
         CALL READI(JUMPTO(IX))
         JDUMP(JUMPTO(IX))=.TRUE.
         IF (NITEMS.GT.3) CALL READI(JUMPINT(IX))

      ELSE IF (WORD.EQ.'LB2') THEN
         LB2T=.TRUE.
      ELSE IF (WORD.EQ.'LOCALSAMPLE') THEN
         LOCALSAMPLET=.TRUE.
         IF (NITEMS.EQ.2) THEN
            CALL READF(ABTHRESH)
         ELSEIF (NITEMS.EQ.3) THEN
            CALL READF(ABTHRESH)
            CALL READF(ACTHRESH)
         ENDIF

C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'LJMF') THEN
         LJMFT=.TRUE.
C        CALL LJPARAMMF
         WRITE(MYUNIT,'(A)') 'LJMF not currently maintained'
         STOP

C start = an N-oligomer is constructed by relocating NREPEAT units and placing them
C at random distance R with Rmin <= R <= Rmax and angle alpha in the xy-plane
C transxy: rigid body translation only in the xy-plane
C rotz:  rigid body rotation only around the z-axis
C dihesc: only perturbation to the side chains
C
      ELSE IF (WORD.EQ.'MAKEOLIGO') THEN
         MAKEOLIGOT=.TRUE.
         CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'START') THEN
            MAKEOLIGOSTART=.TRUE.
            CALL READI(NFIXSEG)
            CALL READI(NREPEAT)
            ALLOCATE(REPATF(NREPEAT),REPATL(NREPEAT),REPPHIF(NREPEAT),REPPHIL(NREPEAT))
            REPATF(1:NREPEAT)=0
            REPATL(1:NREPEAT)=0
            REPPHIF(1:NREPEAT)=0.D0
            REPPHIL(1:NREPEAT)=0.D0
            DO J1=1,NREPEAT
               CALL READI(REPATF(J1))
               CALL READI(REPATL(J1))
               CALL READF(REPPHIF(J1))
               CALL READF(REPPHIL(J1))
            ENDDO
            CALL READF(PLACERMIN)
            CALL READF(PLACERMAX)
            TRANSXYT=.TRUE.
            ROTZT=.TRUE.
            NOPHIPSIT=.TRUE.
         ELSEIF (TRIM(ADJUSTL(UNSTRING)).EQ.'INITROT') THEN
            MAKEOLIGOSTART=.TRUE.
            INITROT=.TRUE.
            CALL READI(NFIXSEG)
            CALL READI(NREPEAT)
            IF (NREPEAT.GT.0) THEN
               ALLOCATE(REPATF(NREPEAT),REPATL(NREPEAT),REPPHIF(NREPEAT),REPPHIL(NREPEAT))
               REPATF(1:NREPEAT)=0
               REPATL(1:NREPEAT)=0
               REPPHIF(1:NREPEAT)=0.D0
               REPPHIL(1:NREPEAT)=0.D0
               DO J1=1,NREPEAT
                  CALL READI(REPATF(J1))
                  CALL READI(REPATL(J1))
                  CALL READF(REPPHIF(J1))
                  CALL READF(REPPHIL(J1))
               ENDDO
               CALL READF(PLACERMIN)
               CALL READF(PLACERMAX)
            ENDIF
            TRANSXYT=.TRUE.
            IF (NREPEAT.GT.0) THEN
               IF (NITEMS.GT.(6+NREPEAT*4)) THEN
                  CALL READA(UNSTRING)
                  IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') NOPHIPSIT=.TRUE.
               ENDIF
            ELSE
               IF (NITEMS.GT.4) THEN
                  CALL READA(UNSTRING)
                  IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') NOPHIPSIT=.TRUE.
               ENDIF
            ENDIF
         ELSEIF (TRIM(ADJUSTL(UNSTRING)).EQ.'REFINE') THEN
            MAKEOLIGOSTART=.FALSE.
            IF (NITEMS.GT.2) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'TRANSXY') TRANSXYT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ROTZ') ROTZT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') NOPHIPSIT=.TRUE.
            IF (NITEMS.GT.3) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'TRANSXY') TRANSXYT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ROTZ') ROTZT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') NOPHIPSIT=.TRUE.
            IF (NITEMS.GT.4) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'TRANSXY') TRANSXYT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ROTZ') ROTZT=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SCONLY') NOPHIPSIT=.TRUE.
         ELSE
            WRITE(MYUNIT,'(A)') 'The first argument to MAKEOLIGO has to be START, INITROT or REFINE - quit.'
            STOP
         ENDIF

      ELSE IF (WORD.EQ.'MAXBFGS') THEN
         CALL READF(MAXBFGS)

      ELSE IF (WORD.EQ.'MAXERISE') THEN
         CALL READF(MAXERISE)

      ELSE IF (WORD.EQ.'MAXIT') THEN
         CALL READI(IX)
         MAXIT=IX
         IF (NITEMS.GT.2) THEN
            CALL READI(IX)
            MAXIT2=IX
         ENDIF

      ELSE IF (WORD.EQ.'MGGLUE') THEN
         MGGLUET=.TRUE.

      ELSE IF (WORD.EQ.'MORSE') THEN
         MORSET=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(XX)
            RHO=XX
         ENDIF
C
C  MPI keyword
C
      ELSE IF (WORD.EQ.'MPI') THEN
         MPIT=.TRUE.
         DEALLOCATE(FIXSTEP,FIXTEMP,FIXBOTH,TEMP,ACCRAT,STEP,ASTEP,OSTEP,BLOCK,NT,JUMPMOVE,JUMPINT,JDUMP,COORDS,NQ,
     @              JUMPTO,EPREV,COORDSO,VAT,VATO,SHELLMOVES,PTGROUP,NSURFMOVES,NCORE)
         ALLOCATE(FIXSTEP(NPAR),FIXTEMP(NPAR),FIXBOTH(NPAR),TEMP(NPAR),ACCRAT(NPAR),STEP(NPAR),ASTEP(NPAR),OSTEP(NPAR),
     @         BLOCK(NPAR),NT(NPAR),JUMPMOVE(NPAR),JUMPINT(NPAR),JDUMP(NPAR),COORDS(3*NATOMS,NPAR),NQ(NPAR),
     @         JUMPTO(NPAR),EPREV(NPAR),
     @         COORDSO(3*NATOMS,NPAR),VAT(NATOMS,NPAR),VATO(natoms,NPAR))
         ALLOCATE(SHELLMOVES(NPAR))
         ALLOCATE(PTGROUP(NPAR))
         ALLOCATE(NSURFMOVES(NPAR))
         ALLOCATE(NCORE(NPAR))
         DO JP=1,NPAR
            FIXSTEP(JP)=.FALSE.
            FIXTEMP(JP)=.FALSE.
            FIXBOTH(JP)=.FALSE.
            TEMP(JP)=0.3D0
            ACCRAT(JP)=0.5D0
            STEP(JP)=0.3D0
            ASTEP(JP)=0.3D0
            OSTEP(JP)=0.3D0
            BLOCK(JP)=0
            NT(JP)=0
            JUMPMOVE(JP)=.FALSE.
            JUMPINT(JP)=100
            JDUMP(JP)=.FALSE.
            SHELLMOVES(JP)=.FALSE.
            PTGROUP(JP)='    '
            NSURFMOVES(JP)=0
            NCORE(JP)=0
         ENDDO

      ELSE IF (WORD.EQ.'MSORIG') THEN
         MSORIGT=.TRUE.

      ELSE IF (WORD.EQ.'MSTRANS') THEN
         MSTRANST=.TRUE.

      ELSE IF (WORD.EQ.'MULLERBROWN') THEN
         MULLERBROWNT=.TRUE.

      ELSE IF (WORD.EQ.'MULTIPLICITY') THEN
         CALL READI(XMUL)
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'MYSD') THEN
         MYSDT=.TRUE.
         CALL READF(SDTOL)

      ELSE IF (WORD.EQ.'NATB') THEN
         NATBT=.TRUE.

      ELSE IF (WORD.EQ.'NEON') THEN
         NEON=.TRUE.
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'NEWCONF') THEN
         NEWCONFT=.TRUE.
         CALL READI(NEWCONFST)
         CALL READF(NCWALL)

      ELSE IF (WORD.EQ.'NEWJUMP') THEN
         NEWJUMP=.TRUE.
         IF (NITEMS.GT.1) CALL READF(PNEWJUMP)
C
C  Reseed runs if the energy does not decrease within NRELAX mc steps.
C  NHSRESTART defines the number of hard sphere moves used to produce the new starting
C  configuration. If NHSRESTART=0 then the geometry is changed using RESEED.
C
      ELSE IF (WORD.EQ.'NEWRESTART') THEN
         NEWRESTART=.TRUE.
         NEWRESTART_MD = .FALSE.                   ! lb415
         NEWRES_TEMP = 0.0D0
         IF (NITEMS.GT.1) CALL READI(NRELAX)
         IF (NITEMS.GT.2) CALL READI(NHSRESTART)
         IF (NITEMS.GT.3) CALL READA(WORD2)        ! lb415
         IF (WORD2.EQ.'MD') NEWRESTART_MD = .TRUE. ! lb415
         IF (NITEMS.GT.4) THEN                     ! lb415
            CALL READF(NEWRES_TEMP)
         ELSE
            NEWRES_TEMP = 1000.0D0
            WRITE(MYUNIT,'(A)') 'keyword> WARNING - temperature unspecified for NEWRESTART. Default for high T MD is 1000K'
         ENDIF                                     ! lb415
         IF (.NOT.ALLOCATED(MSBE)) ALLOCATE(MSBE(MAXSAVE))
         IF (.NOT.ALLOCATED(MSBCOORDS)) ALLOCATE(MSBCOORDS(3*NATOMS,MAXSAVE))

      ELSE IF (WORD.EQ.'NMAX') THEN
         CALL READF(NMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'NMAX=  ',NMAX

      ELSE IF (WORD.EQ.'NMIN') THEN
         CALL READF(NMIN)
         WRITE(MYUNIT,'(A,F14.10)') 'NMIN=  ',NMIN
 
      ELSE IF (WORD.EQ.'NOCHIRALCHECKS') THEN
         CHECKCHIRALITY=.FALSE.

      ELSE IF (WORD.EQ.'UACHIRAL') THEN
         UACHIRAL=.TRUE.

      ELSE IF (WORD.EQ.'NOCISTRANS') THEN
         IF (NITEMS.GT.1) CALL READF(MINOMEGA)

      ELSE IF (WORD.EQ.'NOCISTRANSDNA') THEN
         NOCISTRANSDNA=.TRUE.
         IF (NITEMS.GT.1) CALL READF(MINOMEGA)

      ELSE IF (WORD.EQ.'NOCISTRANSRNA') THEN
         NOCISTRANSRNA=.TRUE.
         IF (NITEMS.GT.1) CALL READF(MINOMEGA)
      
      ELSE IF (WORD.EQ.'NOFREEZE') THEN
         FREEZECORE=.FALSE.

      ELSE IF (WORD.EQ.'NORESET') THEN
         NORESET=.TRUE.
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'OH') THEN
         FIELDT=.TRUE.
         OHT=.TRUE.
         CALL READF(XX)
         FOH=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(EXPFAC)
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(EXPD)
         ENDIF
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'OTP') THEN
         OTPT=.TRUE.
         RIGID=.TRUE.
C        CALL OTPPARAMMF
      ELSE IF (WORD.EQ.'P46') THEN
         P46=.TRUE.
         BLNT=.TRUE.

      ELSE IF (WORD.EQ.'PACHECO') THEN
         PACHECO=.TRUE.

      ELSE IF (WORD.EQ.'PAH') THEN
         PAHT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=36
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD.EQ.'PAIRDIST') THEN
! csw34> PAIRDIST allows the monitoring of the distance between pairs of
! atoms during a BH run. The atom pairs are read from either as
! arguements to the PAIRDIST keyword, or from the pairdist file
         NPAIRS=0
         IF (NITEMS.GT.1) THEN
! If arguements are specified, assume NITEMS/2 pairs on the line 
            PAIRDISTT=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> Pairwise atom distances will be output to pairdists*'
            NPAIRS=(NITEMS-1)/2
            ALLOCATE(PAIRDIST(NPAIRS,2))
            DO J1=1,NPAIRS
               CALL READI(PAIRDIST(J1,1))
               CALL READI(PAIRDIST(J1,2))
               IF (PAIRDIST(J1,1).GT.NATOMS) THEN
                  WRITE(MYUNIT,'(A)') ' keyword> ERROR: PAIRDIST atom index larger than system specified!'
                  STOP
               ELSEIF (PAIRDIST(J1,2).GT.NATOMS) THEN
                  WRITE(MYUNIT,'(A)') ' keyword> ERROR: PAIRDIST atom index larger than system specified!'
                  STOP
               ENDIF
            ENDDO
         ELSE
! If there are no atoms specified on the PAIRDIST line, assume reading them from the 'pairdist' file
! First step - check the pairdist file is present
            YESNO=.FALSE.
            INQUIRE(FILE='pairdist',EXIST=YESNO)
            IF (YESNO) THEN
               PAIRDISTT=.TRUE.
               WRITE(MYUNIT,'(A)') ' keyword> Pairwise atom distances will be output to pairdists*'
            ELSE
               WRITE(MYUNIT,'(A)') ' keyword> ERROR: pairdist input file missing for PAIRDIST'
               FLUSH(MYUNIT)
               STOP
            ENDIF
! Determine NPAIRS to allow allocation of the PAIRDIST array 
            OPEN(UNIT=222,FILE='pairdist',status='old')
            DO
               READ(222,*,IOSTAT=iostatus) CHECK1
               IF (iostatus<0) THEN
                  CLOSE(222)
                  EXIT
               ELSE 
                  NPAIRS=NPAIRS+1
               ENDIF
            END DO        
            CLOSE(222)
! Allocate the PAIRDIST array and read the pairs in
            ALLOCATE(PAIRDIST(NPAIRS,2))
            OPEN(UNIT=222,FILE='pairdist',status='old')
            DO J1=1,NPAIRS
               READ(222,*) PAIRDIST(J1,1),PAIRDIST(J1,2)
            ENDDO
            CLOSE(222)
         ENDIF
! Print list of pairs to GMIN output for checking
         WRITE(MYUNIT,'(A)') ' keyword> Atom pair distances to monitor:'
         DO J1=1,NPAIRS
            WRITE(MYUNIT,*) PAIRDIST(J1,:)
         ENDDO

C
C  PARALLEL must come before STEP and ACCRAT
C  This keyword is for the serial parallel implementation - now obsolete.
C
      ELSE IF (WORD.EQ.'PARALLEL') THEN
         PARALLELT=.TRUE.
         CALL READI(NPAR)
         DEALLOCATE(FIXSTEP,FIXTEMP,FIXBOTH,TEMP,ACCRAT,STEP,ASTEP,OSTEP,BLOCK,NT,JUMPMOVE,JUMPINT,JDUMP,COORDS,NQ,
     @              JUMPTO,EPREV,COORDSO,VAT,VATO,SHELLMOVES,PTGROUP,NSURFMOVES,NCORE) 
         ALLOCATE(FIXSTEP(NPAR),FIXTEMP(NPAR),FIXBOTH(NPAR),TEMP(NPAR),ACCRAT(NPAR),STEP(NPAR),ASTEP(NPAR),OSTEP(NPAR), 
     @         BLOCK(NPAR),NT(NPAR),JUMPMOVE(NPAR),JUMPINT(NPAR),JDUMP(NPAR),NQ(NPAR),JUMPTO(NPAR),COORDS(3*NATOMS,NPAR),
     @         COORDSO(3*NATOMS,NPAR),VAT(NATOMS,NPAR),VATO(NATOMS,NPAR),EPREV(NPAR),SHELLMOVES(NPAR),PTGROUP(NPAR),
     @         NSURFMOVES(NPAR),NCORE(NPAR))
         DO JP=1,NPAR
            FIXSTEP(JP)=.FALSE.
            FIXTEMP(JP)=.FALSE.
            FIXBOTH(JP)=.FALSE.
            TEMP(JP)=0.3D0
            ACCRAT(JP)=0.5D0
            STEP(JP)=0.3D0
            ASTEP(JP)=0.3D0
            OSTEP(JP)=0.3D0
            BLOCK(JP)=0
            NT(JP)=0
            JUMPMOVE(JP)=.FALSE.
            JUMPINT(JP)=100
            JDUMP(JP)=.FALSE.
            SHELLMOVES(JP)=.FALSE.
            PTGROUP(JP)='    '
            NSURFMOVES(JP)=0
            NCORE(JP)=0
         ENDDO

      ELSE IF (WORD.EQ.'PBGLUE') THEN
         PBGLUET=.TRUE.

      ELSE IF (WORD.EQ.'PERIODIC') THEN
         PERIODIC=.TRUE.
         CALL READF(XX)
         BOXLX=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            BOXLY=XX
            IF (NITEMS.GT.3) THEN
               CALL READF(XX)
               BOXLZ=XX
            ENDIF
         ELSE
            BOXLY=BOXLX
            BOXLZ=BOXLX
         ENDIF
!
! If permdist is set then distance calculations are performed with minperdist instead
! of newmindist in procedures such as AVOID and CSM. This keyword is now independent
! from PERMOPT
!
      ELSE IF (WORD.EQ.'PERMDIST') THEN
         PERMDIST=.TRUE.
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
         ALLOCATE(NPERMSIZE(NATOMS),PERMGROUP(NATOMS),NSETS(NATOMS),SETS(NATOMS,2))
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
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
!              IF (NDUMMY+NPERMSIZE(J1).GT.NATOMS) THEN
               IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),
!    &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
!
!  And another sanity check!
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
            NSETS(1:NATOMS)=0
            NPERMGROUP=1 ! all atoms can be permuted - default
            NPERMSIZE(1)=NATOMS ! all atoms can be permuted - default
            DO J1=1,NATOMS
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         WRITE(MYUNIT,'(A,I6)') ' keyword> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            PRINT '(A,3(I6,A))',' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
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
!
!  This keyword is for optimising the distance between permutational isomers.
!
      ELSE IF (WORD.EQ.'PERMOPT') THEN
         PERMOPT=.TRUE.
         PERMDIST=.TRUE.
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
         ALLOCATE(NPERMSIZE(NATOMS),PERMGROUP(NATOMS),NSETS(NATOMS),SETS(NATOMS,2))
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
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
               IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),
!    &                                 ((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))

               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
!
!  And another sanity check!
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
            NSETS(1:NATOMS)=0
            NPERMGROUP=1 ! all atoms can be permuted - default
            NPERMSIZE(1)=NATOMS ! all atoms can be permuted - default
            DO J1=1,NATOMS
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         WRITE(MYUNIT,'(A,I6)') ' keyword> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            PRINT '(A,3(I6,A))',' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
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
      ELSE IF (WORD.EQ.'PLUS') THEN
         PLUS=.TRUE.

      ELSE IF (WORD.EQ.'PMAX') THEN
         CALL READF(PMAX)
         WRITE(MYUNIT,'(A,F14.10)') 'PMAX=  ',PMAX

      ELSE IF (WORD.EQ.'PMIN') THEN
         CALL READF(PMIN)
C
C  POWER provides a means to set the initial premultiplication factor for the
C  gradient in MYLINMIN
C
      ELSE IF (WORD.EQ.'POWER') THEN
         CALL READI(IX)
         MYPOWER=IX
C
C  Purify the geometry in mylbfgs to preserve icosahedral (I)
C  symmetry.
C
      ELSE IF (WORD.EQ.'PROJI') THEN
         PROJIT=.TRUE.
C
C  Purify the geometry in mylbfgs to preserve icosahedral (Ih)
C  symmetry.
C
      ELSE IF (WORD.EQ.'PROJIH') THEN
         PROJIHT=.TRUE.
C
C  Frequency of printing in lbfgs to reduce size of output files
C  NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'PRTFRQ') THEN
         CALL READI(PRTFRQ)
C
C  Plain parallel tempering. Same as BSPT but without quenches.
C
      ELSE IF (WORD.EQ.'PTMC') THEN
         PTMC=.TRUE.
         CALL READF(HISTMIN)
         CALL READF(HISTMAX)
         CALL READF(PTEMIN)
         CALL READF(PTEMAX)
         CALL READF(PTTMIN)
         CALL READF(PTTMAX)
         CALL READF(EXCHPROB)
         CALL READF(NEQUIL)
         CALL READF(PTSTEPS)
         NQUENCH=0.0D0
         CALL READI(NENRPER)
         CALL READI(HBINS)
         QUENCHFRQ=1 ! must be set to avoid division by zero in bspt.F
!        CALL READI(QUENCHFRQ)
C
C  Keyword for applied static force.
C
      ELSE IF (WORD.EQ.'PULL') THEN
         PULLT=.TRUE.
         CALL READI(PATOM1)
         CALL READI(PATOM2)
         CALL READF(PFORCE)
         WRITE(MYUNIT,'(A,I6,A,I6,A,G20.10)') ' keyword> Pulling atoms ',PATOM1,' and ',PATOM2,' force=',PFORCE
C
C Request calculation of structural order parameter Q4 on the fly 
C NOT DOCUMENTED.
C
      ELSE IF (WORD.EQ.'Q4') THEN
         Q4T=.TRUE.
C
C Distance cut-off for Coulomb interactions in AMBER.
C
      ELSE IF (WORD.EQ.'QCUTOFF') THEN
         AMCUT=.TRUE.
         CALL READF(REALQCUTOFF)
         QCUTOFF=1.1D0*REALQCUTOFF
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'QDTEST2') THEN
         QD2T=.TRUE.
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'QDTEST') THEN
         QDT=.TRUE.

      ELSE IF ((WORD.EQ.'QMAX').OR.(WORD.EQ.'TIGHTCONV')) THEN
         CALL READF(CQMAX)

      ELSE IF (WORD.EQ.'QUAD') THEN
         QUADT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=5
         ALLOCATE(SITE(NRBSITES,3))
C
C Collect data from quenches for hopeful conversion into relative density of states.
C Requires saved lowestdirect and firstfit files to be present.
C NOT DOCUMENTED.
C
      ELSE IF (WORD.EQ.'QUENCHDOS') THEN
         QUENCHDOS=.TRUE.
         IF (NITEMS.GT.1) CALL READI(QDLIMIT)

      ELSE IF (WORD.EQ.'RADIUS') THEN
         CALL READF(XX)
         RADIUS=XX
C
C  integer seed for random number generator.
C
      ELSE IF (WORD.EQ.'RANSEED') THEN
         CALL READI(NDUMMY)
         CALL SDPRND(NDUMMY)

      ELSE IF (WORD.EQ.'RCUTOFF') THEN
         AMCUT=.TRUE.
         CALL READF(REALRCUTOFF)
         RCUTOFF=1.1D0*REALRCUTOFF
C
C  Read data for previous geometries that lead to reseeding, which
C  probably approximate MSB bottoms.
C  NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'READMSB') THEN
         INQUIRE(FILE='MSBdata',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            WRITE(MYUNIT,'(A)') 'ERROR - READMSB specified, but no MSBdata data file found'
            STOP
         ELSE
            IF (.NOT.ALLOCATED(MSBCOORDS)) ALLOCATE(MSBCOORDS(3*NATOMS,MAXSAVE))
            IF (.NOT.ALLOCATED(MSBE)) ALLOCATE(MSBE(MAXSAVE))
            OPEN(UNIT=34,FILE='MSBdata',STATUS='OLD')
57          READ(34,*,END=56) DUMMY
            NMSBSAVE=NMSBSAVE+1
            MSBE(NMSBSAVE)=DUMMY
            READ(34,*) (MSBCOORDS(J1,NMSBSAVE),J1=1,3*NATOMS)
            IF (NMSBSAVE.LT.MAXSAVE) GOTO 57
56          WRITE(MYUNIT,'(A,I6,A)') 
     1         'Energies and coordinates read for ',NMSBSAVE,' previous structures from MSBdata'
            CLOSE(34)
         ENDIF
C
C  Renormalisation attempt
C
C  TRENORM is the temperature for the Metropolis accept/reject comparison
C          of lowest energies calculated over NRENORM steps having moved
C          XMOVERENORM atoms. NRENORM is dynamically adjusted with a
C          minimum value equal to half the original NRENORM.
C  NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'RENORM') THEN
         RENORM=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NRENORM)
         IF (NITEMS.GT.2) CALL READF(XMOVERENORM)
         IF (NITEMS.GT.3) CALL READF(TRENORM)
         IF (NITEMS.GT.4) CALL READI(NRENSTUCK)

      ELSE IF (WORD.EQ.'RESIZE') THEN
         RESIZET=.TRUE.
         CALL READF(XX)
         RESIZE=XX
C
C  Reseed runs if a step is not accepted in twice the relaxation time,
C  defined in terms of a number of mc steps NRELAX. NHSRESTART defines
C  the number of hard sphere moves used to produce the new starting
C  configuration. 
C
      ELSE IF (WORD.EQ.'RESTART') THEN
         RESTART=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NRELAX)
         IF (NITEMS.GT.2) CALL READI(NHSRESTART)
C
C  Restore the state of a previous GMIN run from dumpfile.
C
      ELSE IF (WORD.EQ.'RESTORE') THEN
         RESTORET=.TRUE.
         CALL READA(DUMPFILE)
         IF (NITEMS.GT.2) THEN
            CALL READA(INTEDUMPFILE)
            INTERESTORE=.TRUE.
         ENDIF

      ELSE IF (WORD.EQ.'RGCL2') THEN
         RGCL2=.TRUE.

      ELSE IF (WORD.EQ.'RINGROTSCALE') THEN
         CALL READF(RINGROTSCALE)

      ELSE IF (WORD.EQ.'RKMIN') THEN
         RKMIN=.TRUE.
         IF (NITEMS.GT.1) CALL READF(GMAX)
         IF (NITEMS.GT.2) CALL READF(EPS)
         WRITE(MYUNIT,'(A,2L5)') 'RKMIN branch RKMIN,BSMIN=',RKMIN,BSMIN

      ELSE IF (WORD.EQ.'RMS') THEN
         RMST=.TRUE.
         CALL READF(RMSLIMIT)
         CALL READF(RMSTOL)
         CALL READI(RMSSAVE)
         CALL READI(J1)
         CALL READI(J2)
         IF(J1.EQ.1) THEN
           SELECTT=.TRUE.
         ELSE
           SELECTT=.FALSE.
         ENDIF
         IF (J2.EQ.1) PROGRESS=.TRUE.
         WRITE(MYUNIT,'(A)') 'RMST set'

      ELSE IF (WORD.EQ.'SAVE') THEN
         CALL READI(NSAVE)
         IF (A9INTET.AND.(NSAVEINTE.EQ.0)) NSAVEINTE=NSAVE

      ELSE IF (WORD.EQ.'SAVEINTE') THEN
         CALL READI(NSAVEINTE)

      ELSE IF (WORD.EQ.'SC') THEN
         SCT=.TRUE.
         CALL READI(IX)
         NN=IX
         CALL READI(IX)
         MM=IX
         CALL READF(XX)
         SIG=XX
         CALL READF(XX)
         SCEPS=XX
         CALL READF(XX)
         SCC=XX

      ELSE IF (WORD.EQ.'SEED') THEN
         SEEDT=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(IX)
            NSSTOP=IX
         ENDIF
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'SHELLMOVE') THEN
         MOVESHELLT=.TRUE.
         CALL READI(SHELLMOVEMAX)
         CALL READF(SHELLPROB)
         CALL READF(COREFRAC)

      ELSE IF (WORD.EQ.'SHIFTCUT') THEN
         SHIFTCUT=.TRUE.
         CALL READF(XX)
         CUTOFF=XX

      ELSE IF (WORD.EQ.'SIDESTEP') THEN
         CALL READF(SIDESTEP)
         WRITE(MYUNIT,'(A,F14.10)') 'SIDESTEP=  ',SIDESTEP

      ELSE IF (WORD.EQ.'SORT') THEN
         SORTT=.TRUE.

C     ELSE IF (WORD.EQ.'SQUEEZE') THEN
C        CALL READI(NVEC)
C        SQUEEZET=.TRUE.
C        IF (NITEMS.GT.2) THEN
C           CALL READF(XX)
C           SQUEEZER=XX
C        ENDIF
C        IF (NITEMS.GT.3) THEN
C           CALL READF(XX)
C           SQUEEZED=XX
C        ENDIF

      ELSE IF (WORD.EQ.'STAR') THEN
         STAR=.TRUE.
C
C Read in the maximum initial step size, factor for determining angular
C moves, and for rigid bodies the angular step size and the size of the
C blocks for Cartesian and angular moves.
C
C For parallel runs different values can be used for different runs by
C adding additional "STEP" lines to the data file. Otherwise the
C parameters for subsequent parallel runs are set to the values for the
C first one.
C
      ELSE IF (WORD.EQ.'STEP') THEN
         NPCOUNT=NPCOUNT+1
         IF (NPCOUNT.GT.NPAR) THEN
            WRITE(MYUNIT,'(A)') 'Number of STEP lines exceeds NPAR - quit'
            STOP
         ENDIF
         CALL READF(STEP(NPCOUNT))
         CALL READF(ASTEP(NPCOUNT))
         IF (NITEMS.GT.3) CALL READF(OSTEP(NPCOUNT))
         IF (NITEMS.GT.4) CALL READI(BLOCK(NPCOUNT))
C
C  Steered minimisation. This is for basin-hopping steps involving two well-defined
C  objects, e.g. a ligand and a protein.
C
      ELSE IF (WORD.EQ.'STEEREDMIN') THEN
         STEEREDMINT=.TRUE.
         CALL READF(SMINK)          ! final value of force constant
         CALL READF(SMINKINC)        ! increment of force constant per LBFGS step
         CALL READF(SMINDISTSTART)  ! initial distance for atoms SMINATOMA and SMINATOMB
         CALL READF(SMINDISTFINISH) ! final distance for atoms SMINATOMA and SMINATOMB (force turned off)
         CALL READI(SMINATOMA)      ! Atom A in the body to be rotated
         CALL READI(SMINATOMB)      ! Atom B in the other body (fixed for step)

      ELSE IF (WORD.EQ.'TRACKDATA') THEN
         TRACKDATAT=.TRUE.     
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'STEPOUT') THEN
         STEPOUT=.TRUE.

      ELSE IF (WORD.EQ.'STEPS') THEN
         NRUNS=1
         CALL READI(IX)
         MCSTEPS(1)=IX
         IF (NITEMS.GT.2) THEN
         CALL READF(XX)
         TFAC(1)=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            NRUNS=2
            CALL READI(IX)
            MCSTEPS(2)=IX
            CALL READF(XX)
            TFAC(2)=XX
         ENDIF
         IF (NITEMS.GT.5) THEN
            NRUNS=3
            CALL READI(IX)
            MCSTEPS(3)=IX
            CALL READF(XX)
            TFAC(3)=XX
         ENDIF

      ELSE IF (WORD.EQ.'STICKY') THEN
         STICKYT=.TRUE.
         RIGID=.TRUE.
         CALL READI(NRBSITES)
         CALL READF(STICKYSIG)
!        WRITE(MYUNIT,*) 'NRBSITES=',NRBSITES 
!        WRITE(MYUNIT,*) 'STICKYSIG=',STICKYSIG 
         ALLOCATE(SITE(NRBSITES,3))
         DO J1=1,NRBSITES
            READ(5,*) SITE(J1,1:3)
C           CALL READF(SITE(J1,1))
C           CALL READF(SITE(J1,2))
C           CALL READF(SITE(J1,3))
!           WRITE(MYUNIT,'(A,I5,3G20.10)') 'J1,site: ',J1,SITE(J1,1:3)
         ENDDO

      ELSE IF (WORD.EQ.'LJCOUL') THEN
         LJCOULT=.TRUE.
         CALL READI(COULN)
         CALL READF(COULQ)
         CALL READF(COULSWAP)
         CALL READF(COULTEMP)
         NRBSITES=1
         ALLOCATE(SITE(NRBSITES,3))
C        Maybe the above two lines are not necessary!

      ELSE IF (WORD.EQ.'STOCK') THEN
         STOCKT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=1
         CALL READF(STOCKMU)
         CALL READF(STOCKLAMBDA)
         ALLOCATE(SITE(NRBSITES,3))

C       Anisotropic potentials:

!     DC430 >

      ELSE IF (WORD.EQ.'CHECKD') THEN
         CHECKDT = .TRUE.

      ELSE IF (WORD.EQ.'DBP') THEN

         DBPT   = .TRUE.
         RIGID  = .TRUE.
         CALL READF(DBEPSBB)
         CALL READF(DBSIGBB)
         CALL READF(DBPMU)
         IF (NITEMS > 4) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF 
         NRBSITES = 3
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD.EQ.'DBPTD') THEN

         DBPTDT = .TRUE.
         RIGID  = .TRUE.
         CALL READF(DBEPSBB)
         CALL READF(DBSIGBB)
         CALL READF(DBPMU)
         IF (NITEMS > 4) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         NRBSITES = 3 
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD .EQ.'DBLPY') THEN

         DBLPYT = .TRUE.
         CALL READF(PYA1(1))
         CALL READF(PYA1(2))
         CALL READF(PYA1(3))
         CALL READF(PYA2(1))
         CALL READF(PYA2(2))
         CALL READF(PYA2(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)

         RIGID     = .TRUE.

         IF (PYA1(1) == PYA2(1) .AND. PYA1(2) == PYA2(2) .AND. PYA1(3) == PYA2(3)) THEN
            RADIFT = .FALSE.
         ELSE
            RADIFT = .TRUE.
         ENDIF

         ESA = PYA1

      ELSE IF (WORD.EQ.'DMBLM') THEN

         DMBLMT   = .TRUE.
         RIGID    = .TRUE.
         CALL READF(EPS11)
         CALL READF(EPS22)
         CALL READF(MRHO11)
         CALL READF(MRHO22)
         CALL READF(REQ11)
         CALL READF(REQ22)
         CALL READF(DBPMU)
         IF (NITEMS > 8) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         NRBSITES = 2
         ALLOCATE(SITE(NRBSITES,3))
         ALLOCATE(RBUV(NRBSITES,3))
         ALLOCATE(DPMU(NRBSITES))

         CALL DEFDMBL()

      ELSE IF (WORD.EQ.'LINROD') THEN

         LINRODT = .TRUE.
         RIGID   = .TRUE.

         NRBSITES = 4
         ALLOCATE(SITE(NRBSITES,3))

         CALL DEFLINROD()

      ELSE IF (WORD.EQ.'LWOTP') THEN

         LWOTPT = .TRUE.
         RIGID  = .TRUE.
         IF (NITEMS > 1) CALL READF(LWRCUT)
         NRBSITES = 3
         ALLOCATE(SITE(NRBSITES,3)) 
 
         CALL DEFLWOTP()

      ELSE IF (WORD.EQ.'NCAP') THEN

         NCAPT  = .TRUE.
         RIGID  = .TRUE.
         CALL READF(RHO)
         CALL READF(EPS2)
         CALL READF(RAD)
         CALL READF(HEIGHT)
         NRBSITES = 6
         ALLOCATE(SITE(NRBSITES,3))

         CALL DEFCAPSID(RAD,HEIGHT)

      ELSE IF (WORD .EQ. 'NPAH') THEN

         CALL READI(PAHID)         
         NPAHT    = .TRUE.
         RIGID    = .TRUE.
         IF (PAHID == 1) THEN
            NRBSITES = 36
         ELSEIF (PAHID == 2) THEN
            NRBSITES = 26
         ELSEIF (PAHID == 3) THEN
            NRBSITES = 12
         ENDIF
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD .EQ. 'NTIP') THEN
 
         NTIPT =.TRUE.
         RIGID =.TRUE.
         IF (NITEMS > 1) THEN 
            CALL READI(TIPID)
         ELSE
            PRINT *, 'ERROR, TIPID is missing'
            STOP
         ENDIF 
         IF (TIPID == 1) NRBSITES = 3
         IF (TIPID == 2) NRBSITES = 4
         IF (TIPID == 3) NRBSITES = 3
         IF (TIPID == 4) NRBSITES = 4
         IF (TIPID == 5) NRBSITES = 5

         ALLOCATE(SITE(NRBSITES,3))

!|gd351>

      ELSE IF (WORD .EQ. 'PATCHY') THEN
 
         PATCHY =.TRUE.
         RIGID =.TRUE.
         IF (NITEMS > 1) THEN 
            CALL READI(NRBSITES)
         ELSE
            PRINT *, 'ERROR, NRBSITES is missing'
            STOP
         ENDIF 

         ALLOCATE(SITE(NRBSITES,3))

         SIGMASQ = (1.D0)**2
         RANGESQ = (1.9D0)**2
         FACTOR =  (2*3.14159265358979D0*0.05)**2

         CALL DEFINE_PATCHES(7.298D0)

      ELSE IF (WORD .EQ. 'ASAOOS') THEN
 
         ASAOOS =.TRUE.
 
!<gd351|

      ELSE IF (WORD .EQ. 'STOCKAA') THEN

         STOCKAAT = .TRUE.
         RIGID    = .TRUE.
!         CALL READF(MRHO11)
         CALL READF(STOCKMU)
         IF (NITEMS > 2) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         IF (.NOT. EFIELDT) EFIELD = 0.D0
         NRBSITES = 1

      ELSE IF (WORD .EQ. 'MSSTK') THEN

         CALL READI(NRBSITES)
         MSSTOCKT = .TRUE.
         RIGID    = .TRUE.
         ALLOCATE(SITE(NRBSITES,3))
         ALLOCATE(RBUV(NRBSITES,3))
         ALLOCATE(DPMU(NRBSITES))
         DO J1 = 1, NRBSITES
            CALL READF(DPMU(J1))
         ENDDO    
         IF (NRBSITES == 2) THEN
            CALL DEFMULT2STOCK()
         ELSE IF (NRBSITES == 3) THEN
            CALL DEFMULT3STOCK()
         ELSE IF (NRBSITES == 4) THEN
            CALL DEFMULT4STOCK()
         ELSE IF (NRBSITES == 5) THEN
            CALL DEFMULT5STOCK()
         ELSE
            WRITE(*,*) 'NOT ALLOWED NRBSITES=',NRBSITES
            STOP 
         ENDIF
         IF (NITEMS > (2 + NRBSITES)) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         IF (.NOT. EFIELDT) EFIELD = 0.D0

      ELSE IF (WORD .EQ. 'MSBIN') THEN

!         CALL READI(NRBSITES)
!         CALL READI(NRBSITES1)
         CALL READI(NPS)
         NRBSITES  = 11
         NRBSITES1 = 5        
         MSTBINT = .TRUE.
         RIGID    = .TRUE.
         CALL READF(STOCKMU)
         IF (NITEMS > 3) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF
         ALLOCATE(SITE(NRBSITES,3))
         ALLOCATE(RBUV(NRBSITES,3))
         CALL DEFMSTBIN()
         IF (.NOT. EFIELDT) EFIELD = 0.D0

      ELSE IF (WORD .EQ. 'MULTPAHA') THEN

         TPAHA     = 4
         ALLOCATE(NCMP(TPAHA))
         CALL READI (NCMP(1))
         CALL READI (NCMP(2))
         CALL READI (NCMP(3))
         CALL READI (NCMP(4))
         MULTPAHAT = .TRUE.
         RIGID     = .TRUE.

      ELSE IF (WORD .EQ. 'TDHD') THEN

         CALL READF(RHO)
         CALL READF(MREQ)
         CALL READF(EPSR)
         TDHDT    = .TRUE.
         RIGID    = .TRUE.
         NRBSITES = 4
         ALLOCATE(SITE(NRBSITES,3))

         CALL DEFTDHD()

      ELSE IF (WORD .EQ. 'WATERDC') THEN

         WATERDCT = .TRUE.
         RIGID    = .TRUE.
         NRBSITES = 4
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD .EQ. 'WATERKZ') THEN

         WATERKZT = .TRUE.
         RIGID    = .TRUE.
         NRBSITES = 4
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD.EQ.'GB') THEN

         GBT = .TRUE.
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(GBSIGNOT)
         CALL READF(GBEPSNOT)

         RIGID     = .TRUE.
         ESA       = 0.5D0*GBSIGNOT*(/1.D0, 1.D0, GBKAPPA/)
         GBCHI     = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
         GBCHIPRM  = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)

      ELSE IF (WORD.EQ.'GBD') THEN

         GBDT = .TRUE.
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(GBSIGNOT)
         CALL READF(GBEPSNOT)

         RIGID     = .TRUE.
         SIGMAF     = GBSIGNOT * GBKAPPA
         ESA        = 0.5D0*(/GBSIGNOT, GBSIGNOT, SIGMAF/)
         INVKAP     = 1.D0/GBKAPPA
         GBCHI      = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
         GBCHIPRM   = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)

      ELSE IF (WORD.EQ.'GBDP') THEN

         GBDPT = .TRUE.
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(GBSIGNOT)
         CALL READF(GBEPSNOT)
         CALL READF(GBDPMU)
         CALL READF(GBDPEPS)

         RIGID     = .TRUE.
         SIGMAF     = GBSIGNOT * GBKAPPA
         ESA        = 0.5D0*(/GBSIGNOT, GBSIGNOT, SIGMAF/)
         INVKAP     = 1.D0/GBKAPPA
         GBCHI      = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
         IF (GBMU == 0.D0) THEN
            GBCHIPRM = -1.D0
         ELSE
            GBCHIPRM   = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)
         ENDIF
         GBDPFCT    = 3.D0*GBDPEPS*GBDPMU*GBDPMU*GBSIGNOT**3.D0

      ELSE IF (WORD.EQ.'GEM') THEN

         GEMT   = .TRUE.
         CALL READF(GEMRC)

!     ----------------------------------------------------------------------------------------------

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
         RIGID    = .TRUE.
         ALLOCATE(SITE(NRBSITES,3))
         ALLOCATE(RBSTLA(NRBSITES,3))
         ALLOCATE(STCHRG(NRBSITES))

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

!     ----------------------------------------------------------------------------------------------

      ELSE IF (WORD .EQ. 'PAHW99') THEN

         CALL READI(PAHID)

         IF (PAHID == 1) THEN
            NRBSITES = 18
         ENDIF

         PAHW99T  = .TRUE.
         RIGID    = .TRUE.
         ALLOCATE(SITE(NRBSITES,3))

!     ----------------------------------------------------------------------------------------------

      ELSE IF (WORD .EQ.'PYG') THEN

         PYGT = .TRUE.
         CALL READF(PYA1(1))
         CALL READF(PYA1(2))
         CALL READF(PYA1(3))
         CALL READF(PYA2(1))
         CALL READF(PYA2(2))
         CALL READF(PYA2(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)

         RIGID     = .TRUE.

         IF (PYA1(1) == PYA2(1) .AND. PYA1(2) == PYA2(2) .AND. PYA1(3) == PYA2(3)) THEN
            RADIFT = .FALSE.
         ELSE
            RADIFT = .TRUE.
         ENDIF
     
         ESA = PYA1

         NRBSITES = 1 ! required for finalio.f

       ELSE IF (WORD .EQ.'PYGDP') THEN

         PYGDPT = .TRUE.
         CALL READF(PYA1(1))
         CALL READF(PYA1(2))
         CALL READF(PYA1(3))
         CALL READF(PYA2(1))
         CALL READF(PYA2(2))
         CALL READF(PYA2(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)
         CALL READF(PYDPMU)
         CALL READF(PYDPEPS)
         IF (NITEMS > 11) THEN
            CALL READF(EFIELD)
            EFIELDT = .TRUE.
         ENDIF

         RIGID     = .TRUE.

         IF (PYA1(1) == PYA2(1) .AND. PYA1(2) == PYA2(2) .AND. PYA1(3) == PYA2(3)) THEN
            RADIFT = .FALSE.
         ELSE
            RADIFT = .TRUE.
         ENDIF

         PYDPFCT  = 3.D0*PYDPEPS*PYDPMU*PYDPMU*PYSIGNOT**3.D0

         ESA = PYA1

         NRBSITES = 1 ! required for finalio.f

       ELSE IF (WORD .EQ.'MSGB') THEN

         MSGBT = .TRUE.
!         RIGID = .TRUE.
         CALL READI(NGBSITE)
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(SIGNOT)
         CALL READF(EPSNOT)
!         ALLOCATE(SITE(NRBSITES,3))

         LPL    = 0.5D0 * SIGNOT * GBKAPPA
         LPR    = 0.5D0 * SIGNOT
         LPRSQ  = LPR * LPR
         LSQDFR = LPL * LPL - LPRSQ

 
      ELSE IF (WORD .EQ. 'MSPYG') THEN

         MSPYGT = .TRUE.
         CALL READI(NPYSITE)
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

      ELSE IF (WORD .EQ. 'MULTISITEPY') THEN

         MULTISITEPYT = .TRUE.
         RIGID = .TRUE.
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)

      ELSE IF (WORD.EQ.'GAYBERNE') THEN

         GAYBERNET=.TRUE.
         ELLIPSOIDT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=1
         CALL READF(GBANISOTROPYR)
         CALL READF(GBWELLDEPTHR)
         CALL READF(PSIGMA0(1))
         CALL READF(PEPSILON0)
         CALL READF(GBMU)
         CALL READF(GBNU)
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD.EQ.'PARAMONOV') THEN
         PARAMONOVT=.TRUE.
         ELLIPSOIDT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=1
         CALL READF(PARAMa1)
         CALL READF(PARAMb1)
         CALL READF(PARAMc1)
         CALL READF(PARAMa2)
         CALL READF(PARAMb2)
         CALL READF(PARAMc2)
         CALL READF(PSIGMA0(1))
         CALL READF(PSIGMA0(2))
         CALL READF(PEPSILON0)
         IF (NITEMS.GT.10) THEN
            CALL READF(PCUTOFF)
            PARAMONOVCUTOFF=.TRUE.
            PCUTOFF=PCUTOFF*MAX(PSIGMA0(1),PSIGMA0(2))
            write (MYUNIT,*) "PY Potential. PCutoff ON:",PCUTOFF
         ENDIF
         IF (NITEMS.GT.11) THEN
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
                write(*,*) "Paramonov Periodic Boundary Condition X active. BOXLX:",BOXLX
            ENDIF
            IF (SCAN(PBC,'Yy').NE.0) THEN
                PARAMONOVPBCY=.TRUE.
                CALL READF(BOXLY)
                BOXLY=BOXLY*PCUTOFF
                write(*,*) "Paramonov Periodic Boundary Condition Y active. BOXLY:",BOXLY
            ENDIF
            IF (SCAN(PBC,'Zz').NE.0) THEN
                PARAMONOVPBCZ=.TRUE.
                CALL READF(BOXLZ)
                BOXLZ=BOXLZ*PCUTOFF
                write(*,*) "Paramonov Periodic Boundary Condition Z active. BOXLZ",BOXLZ
            ENDIF
         ENDIF 
         ALLOCATE(SITE(NRBSITES,3))
      ELSE IF (WORD.EQ.'PYOVERLAPTHRESH') THEN
         CALL READF(PYOVERLAPTHRESH)
         WRITE(MYUNIT,'(A,F8.3)') 'keyword> ellipsoids considered to overlap for an ECF value of ', PYOVERLAPTHRESH
      ELSE IF (WORD .EQ.'PYGPERIODIC') THEN

         PYGPERIODICT = .TRUE.
         ELLIPSOIDT = .TRUE.
         RIGID = .TRUE.
         CALL READF(PYA1(1))
         CALL READF(PYA1(2))
         CALL READF(PYA1(3))
         CALL READF(PYA2(1))
         CALL READF(PYA2(2))
         CALL READF(PYA2(3))
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)
         PARAMa1=PYA1(1)
         PARAMb1=PYA1(2)
         PARAMc1=PYA1(3)

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
                write(*,*) "Paramonov Periodic Boundary Condition X active. BOXLX:",BOXLX
            ENDIF
            IF (SCAN(PBC,'Yy').NE.0) THEN
                PARAMONOVPBCY=.TRUE.
                CALL READF(BOXLY)
                BOXLY=BOXLY*PCUTOFF
                write(*,*) "Paramonov Periodic Boundary Condition Y active. BOXLY:",BOXLY
            ENDIF
            IF (SCAN(PBC,'Zz').NE.0) THEN
                PARAMONOVPBCZ=.TRUE.
                CALL READF(BOXLZ)
                BOXLZ=BOXLZ*PCUTOFF
                write(*,*) "Paramonov Periodic Boundary Condition Z active. BOXLZ",BOXLZ
            ENDIF
         ENDIF
         ALLOCATE(SITE(NRBSITES,3))
      ELSE IF (WORD .EQ.'LJCAPSID') THEN
!         Three-site Lennard-Jones based capsid model. Sites at the origin are standard LJ sites, the two apex sites 
!         are repulsive LJ sites, polarised. The site in the middle only interacts with sites in the middle of other 
!         molecules.
         LJCAPSIDT = .TRUE.
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)
         CALL READF(PEPSILON1(1))
         CALL READF(PSCALEFAC1(1))
         CALL READF(PSCALEFAC2(1))
         MAXINTERACTIONS=4
        
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD.EQ.'EXTRALJSITE') THEN
         LJSITE=.TRUE.
         CALL READF(PEPSILON1(1))
         CALL READF(PSCALEFAC1(1))
          MAXINTERACTIONS=1
         IF(NITEMS.GT.3) THEN
          CALL READF(PSCALEFAC2(1))
          WRITE(MYUNIT,'(A,3F8.3)') 'keyword> primary and secondary apex sites will be used, epsilon and heights: ', 
     &                              PEPSILON1(1), PSCALEFAC1(1), PSCALEFAC2(1)
          IF(.NOT.LJSITEATTR) THEN
                MAXINTERACTIONS=3
          ELSE
                MAXINTERACTIONS=4
          END IF
         ELSE
          WRITE(MYUNIT,'(A,2F8.3)') 'keyword> primary apex sites will be used, epsilon and height: ', PEPSILON1(1), PSCALEFAC1(1)
         END IF
         IF(NITEMS.GT.4) THEN           ! binary ellipsoidal clusters will be set up only for two apex sites, not one
           BLJSITE=.TRUE.               ! we also won't use the sigma parameter from now on, epsilon is enough for repulsive sites
           WRITE(MYUNIT,'(A,3F8.3)') 'keyword> binary system with primary and secondary apex sites, ' //  
     &  'epsilon and heights for 2nd type particle: ', PEPSILON1(1), PSCALEFAC1(1), PSCALEFAC2(1)

           CALL READF(PEPSILON1(2))
           CALL READF(PSCALEFAC1(2))
           CALL READF(PSCALEFAC2(2))
           CALL READF(PEPSILON1(3))     ! this is epsilon for the interaction between A and B type ellipsoids
           MAXINTERACTIONS=3 ! attractive secondary apex sites not incorporated for binary systems
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
      ELSE IF (WORD.EQ.'SWAPMOVES') THEN
         SWAPMOVEST=.TRUE.
         IF(PYBINARYT) THEN
                PYSWAP(1) = 1
                PYSWAP(2) = PYBINARYTYPE1 + 1
                PYSWAP(3) = 1
                IF(NITEMS.GT.1) CALL READI(PYSWAP(3))
                WRITE(MYUNIT,'(A,I5,A)') 'keyword> ',PYSWAP(3), ' pairs of atoms will be swapped at once'
         END IF
      ELSE IF (WORD.EQ.'PYBINARY') THEN
         PYBINARYT=.TRUE.
         ELLIPSOIDT=.TRUE.
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
         IF(NITEMS.GT.16) THEN
            CALL READF(PCUTOFF)
            PARAMONOVCUTOFF=.TRUE.
            PCUTOFF=PCUTOFF*PYSIGNOT
            write (MYUNIT,*) "PY Potential. PCutoff ON:",PCUTOFF
         END IF
         IF(SWAPMOVEST) THEN
                PYSWAP(1) = 1
                PYSWAP(2) = PYBINARYTYPE1 + 1
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

      ELSE IF (WORD.EQ.'GAYBERNEDC') THEN
         GAYBERNEDCT=.TRUE.
C         ELLIPSOIDT=.TRUE.
         RIGID=.TRUE.
         NRBSITES=1
         CALL READF(GBKAPPA)
         CALL READF(GBKAPPRM)
         CALL READF(GBMU)
         CALL READF(GBNU)
         CALL READF(SIGNOT)
         CALL READF(EPSNOT)
         ALLOCATE(SITE(NRBSITES,3))

      ELSE IF (WORD.EQ.'STRAND') THEN
         STRANDT=.TRUE.
         RIGID=.TRUE.
C
C  The nine reference site positions per strand.
C
         NRBSITES=9
         ALLOCATE(SITE(NRBSITES,3))
         SITE(1,1)=-2.7298862082
         SITE(1,2)=2.3622865625 
         SITE(1,3)=0.6475151629
         SITE(2,1)=-1.7492122114
         SITE(2,2)=2.3331194664 
         SITE(2,3)=0.5887015133
         SITE(3,1)=-1.5963638586
         SITE(3,2)=1.4304320585 
         SITE(3,3)=0.2442792479
         SITE(4,1)=-0.6166461313
         SITE(4,2)=1.4301805389 
         SITE(4,3)=0.1327546571
         SITE(5,1)=-0.4460267836
         SITE(5,2)=0.5254809645  
         SITE(5,3)=-0.2196837962
         SITE(6,1)=0.5313983749 
         SITE(6,2)=0.5210707739  
         SITE(6,3)=-0.3409645197
         SITE(7,1)=0.7065341613  
         SITE(7,2)=-0.3914277962 
         SITE(7,3)=-0.6719579835 
         SITE(8,1)=1.6776397940  
         SITE(8,2)=-0.3830053500 
         SITE(8,3)=-0.8355266604
         SITE(9,1)=1.8162689403  
         SITE(9,2)=-1.3093381947 
         SITE(9,3)=-1.1427874015
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'SUPERSTEP') THEN
         SUPERSTEP=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NSUPER)
         IF (NITEMS.GT.2) CALL READF(SUPSTEP)
         IF (NITEMS.GT.3) CALL READF(TEMPS)
         IF (NITEMS.GT.4) CALL READF(SACCRAT)
         IF (NITEMS.GT.5) CALL READI(NSACCEPT)

      ELSE IF (WORD.EQ.'SW') THEN
         SW=.TRUE.
      ELSE IF (WORD.EQ.'SETCHIRAL') THEN
         SETCHIRAL=.TRUE.
C
C  Keyword and parameters for symmetrisation.
C
      ELSE IF (WORD.EQ.'SYMMETRISE') THEN
         SYMMETRIZE=.TRUE.
         NCORE=0
         IF (NITEMS.GT.1) CALL READI(NSYMINTERVAL)
         IF (NITEMS.GT.2) CALL READF(SYMTOL1)
         IF (NITEMS.GT.3) CALL READF(SYMTOL2)
         IF (NITEMS.GT.4) CALL READF(SYMTOL3)
         IF (NITEMS.GT.5) CALL READF(SYMTOL4)
         IF (NITEMS.GT.6) CALL READF(SYMTOL5)
         IF (NITEMS.GT.7) CALL READI(NSYMQMAX)
         IF (NITEMS.GT.8) CALL READF(MATDIFF) ! appears to have little effect now
         IF (NITEMS.GT.9) CALL READF(DISTFAC)
C
C  Keyword and parameters for symmetrisation according to a continuous symmetry measure.
C
      ELSE IF (WORD.EQ.'SYMMETRISECSM') THEN
         SYMMETRIZE=.TRUE.
         SYMMETRIZECSM=.TRUE.
         CSMT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NSYMINTERVAL)
         IF (NITEMS.GT.2) THEN
            CALL READA(CSMGP)
            CALL MYUPCASE(CSMGP)
         ELSE
            PRINT '(A)','keyword> ERROR - point group must be specified for SYMMETRIZECMS keyword'
            STOP
         ENDIF
         IF (NITEMS.GT.3) CALL READF(CSMEPS)
         IF (NITEMS.GT.4) CALL READI(CSMSTEPS)
         IF (NITEMS.GT.5) CALL READI(CSMQUENCHES)
         IF (NITEMS.GT.6) CALL READI(CSMMAXIT)
         IF (.NOT.PERMDIST) THEN
         PERMDIST=.TRUE.
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
!        ALLOCATE(NPERMSIZE(NATOMS),PERMGROUP(NATOMS),NSWAP(NATOMS),SWAP1(NATOMS,2),SWAP2(NATOMS,2))
         ALLOCATE(NPERMSIZE(3*NATOMS),PERMGROUP(3*NATOMS),NSETS(3*NATOMS),SETS(NATOMS,3))
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
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
!              IF (NDUMMY+NPERMSIZE(J1).GT.NATOMS) THEN
               IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),
!    &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))

               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
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
            NSETS(1:NATOMS)=0
            NPERMGROUP=1 ! all atoms can be permuted - default
            NPERMSIZE(1)=NATOMS ! all atoms can be permuted - default
            DO J1=1,NATOMS
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         WRITE(MYUNIT,'(A,I6)') ' keyword> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            PRINT '(A,3(I6,A))',' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
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
      ELSE IF (WORD.EQ.'TABOO') THEN
         TABOOT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NTAB)

      ELSE IF (WORD.EQ.'TARGET') THEN
         TARGET=.TRUE.
         NTARGETS=NITEMS-1
         ALLOCATE(TARGETS(NTARGETS))
         INQUIRE(FILE='coords.target',EXIST=YESNO)
         IF (YESNO) THEN
            ALLOCATE(TCOORDS(NTARGETS,3*NATOMS))
            OPEN(UNIT=1,FILE='coords.target',STATUS='OLD')
            READ(1,*) ((TCOORDS(J1,J2),J2=1,3*NATOMS),J1=1,NTARGETS)
            CLOSE(1)
         ENDIF
         DO J1=2,NITEMS
            CALL READF(XX)
            TARGETS(J1-1)=XX
         ENDDO
C
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'TD') THEN
         FIELDT=.TRUE.
         TDT=.TRUE.
         CALL READF(XX)
         FTD=XX
         IF (NITEMS.GT.2) THEN
            CALL READF(XX)
            EXPFAC=XX
         ENDIF
         IF (NITEMS.GT.3) THEN
            CALL READF(XX)
            EXPD=XX
         ENDIF

      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         DO J1=1,NITEMS-1
            CALL READF(TEMP(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            DO J1=NITEMS,NPAR
               TEMP(J1)=TEMP(1)
            ENDDO
         ENDIF
C
C Tethered WL walk to determine anharmonic vibrational density of states
C
      ELSE IF (WORD.EQ.'TETHER') THEN
         TETHER=.TRUE.
         CALL READF(hdistconstraint)
         CALL READI(hwindows)
         lhbins=int(hbins/hwindows)
         CALL READF(ExtrapolationPercent)
         lhbins=int(hbins/hwindows)
         sampledbins=int((1.0d0-ExtrapolationPercent)*hbins/hwindows)
         CALL READF(lnHarmFreq)
      ELSE IF (WORD.EQ.'THOMSON') THEN
         THOMSONT=.TRUE.
         ODDCHARGE=1.0D0
         IF (NITEMS.GT.1) CALL READF(ODDCHARGE)
C
C  Threshold acceptance rather than Metropolis, i.e. the energy change
C  can;t increase by more than a certain amount.
C  NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'THRESHOLD') THEN
         THRESHOLDT=.TRUE.
!        WRITE(MYUNIT,*) 'keyword THRESHOLD doesnt appear to do anything at the moment'
         STOP

      ELSE IF (WORD.EQ.'TIP') THEN
         TIP=.TRUE.
         RIGID=.TRUE.
         IF (NITEMS.GT.1) CALL READI(TIPID)
         IF (TIPID.EQ.5) NRBSITES=5
         IF (TIPID.EQ.4) NRBSITES=4
         IF (TIPID.EQ.3) NRBSITES=3
         IF (TIPID.EQ.2) NRBSITES=4
         IF (TIPID.EQ.1) NRBSITES=3
         ALLOCATE(SITE(NRBSITES,3))
C     ELSE IF (WORD.EQ.'TN') THEN
C        TNT=.TRUE.
C        WRITE(MYUNIT,'(A)') 'optimisation with tn no longer supported'
C        STOP

      ELSE IF (WORD.EQ.'TOLBRENT') THEN
         CALL READF(TOLB)
C
C NOT DOCUMENTED
C
      ELSE IF (WORD .EQ. 'TOLD') THEN
        CALL READF(XX)
        TOLD=XX
C
C NOT DOCUMENTED
C
      ELSE IF (WORD .EQ. 'TOLE') THEN
        CALL READF(XX)
        TOLE=XX

      ELSE IF (WORD.EQ.'TOSI') THEN
         TOSI=.TRUE.
         CALL READF(APP)
         CALL READF(AMM)
         CALL READF(APM)
         CALL READF(RHO)
C
C  Set Tsallis statistics with some q value.
C
      ELSE IF (WORD.EQ.'TSALLIS') THEN
         TSALLIST=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READF(QTSALLIS)
         ENDIF

      ELSE IF (WORD.EQ.'TWOPLUS') THEN
         TWOPLUS=.TRUE.
         WRITE(MYUNIT,'(A,F14.10)') 'PMIN=  ',PMIN
C
C  Number of BFGS updates before resetting, default=4
C
      ELSE IF (WORD.EQ.'UPDATES') THEN
         CALL READI(MUPDATE)

C
C  Use VGW (Variational Gaussian Wavepacket) Minimization (Quantum Quenching)
C
      ELSE IF (WORD.EQ.'VGW') THEN
         VGW=.TRUE.
         LBFGST=.FALSE.
         CALL READF(LJSIGMA)
         CALL READF(LJEPSILON)
         CALL READF(TAUMAX)
         CALL READF(TAUMAXFULL)
                 
      ELSE IF (WORD.EQ.'VGWCPS') THEN
         CALL READI(CPS)
         CALL READF(CPFACTORSG)

      ELSE IF (WORD.EQ.'VGWCPF') THEN
         CALL READI(CPF)
         CALL READF(CPFACTORFG)

      ELSE IF (WORD.EQ.'VGWTOL') THEN
         CALL READF(VGWTOL)

C
C Choice of convergence regime for Wang-Landau runs: histogram flatness (default), 
C VisitProp - minimal number of visits proportional to 1/sqrt(ln(f))
C
      ELSE IF (WORD.EQ.'VISITPROP') THEN
         VISITPROP=.TRUE.
C
C Maximum PE for an instantaneous configuration above a basin bottom in BSPT
C
      ELSE IF (WORD.EQ.'TSTAR') THEN
         CALL READF(TSTAR)
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
C NOT DOCUMENTED
C
      ELSE IF (WORD.EQ.'WENZEL') THEN
         WENZEL=.TRUE.

      ELSE IF (WORD.EQ.'ZETT1') THEN
         ZETT1=.TRUE.

      ELSE IF (WORD.EQ.'ZETT2') THEN
         ZETT2=.TRUE.
!op226> </kwd>
         !op226>}}} 
         !op226> <end>
         !op226> Uhh! Went through all the available keywords; now the final
         !op226> ELSE... 
      ELSE
         CALL REPORT('Unrecognized command '//WORD,.TRUE.)
         STOP
      ENDIF
      CALL FLUSH(MYUNIT)
      GOTO 190
      ! }}}
        CASE(2) ! OPTIM {{{
          ! }}}
        CASE(3) ! PATHSAMPLE {{{
!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales                      {{{
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!                                                                }}}
!op226 <begin>
      ! Declarations 
      ! {{{
      ! Modules {{{
      USE PORFUNCS
      USE NODES, ONLY: JPN, GETNODES, NNODES
      USE KEY
      USE COMMON
      USE RIGIDBODYMOD, ONLY: CAPSOMER, NUMRBTYPES, CAPSOMERDEFS
        ! }}}

      ! Variables {{{
      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, ISTAT, NDUMMY2, LAST
      LOGICAL CAT
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR, NERROR, IR, ECHO, LAST, CAT
      INTEGER J1, NDUMMY, J2, J3
      DOUBLE PRECISION DUMMY, DBSIGBB
      LOGICAL END, SKIPBL, CLEAR, ECHO, PERMFILE, RBSYMTEST
      CHARACTER(LEN=20) WW
      CHARACTER(LEN=80) RWENERGYFILE
      CHARACTER(LEN=20) WORD
      CHARACTER(LEN=1) LOYNO
      CHARACTER(LEN=9) UNSTRING

      TYPE (CAPSOMER), ALLOCATABLE :: TEMPCAPSOMERDEFS(:)
      INTEGER :: MAXRBTYPES
      ! }}}
!}}}
!op226 </begin>
!op226 Initializations {{{ 
!op226 <init>
      DBPT     = .FALSE.
      DBPTDT   = .FALSE.
      MSSTOCKT = .FALSE.
      PAHAT    = .FALSE.
      PATCHYDT = .FALSE.
      STOCKAAT = .FALSE.
      EFIELDT  = .FALSE.
      DEBUG=.FALSE.
      EXTRACTMINT=.FALSE.
      EXTRACTTST=.FALSE.
      ADDPT=.FALSE.
      DUMPGROUPST=.FALSE.
      DSCALE=3.0D0
      PSCALE=0.5D0
      CONNECTIONS=0
      MAXTSATTEMPTS=10
      ISEED=1
      NATTEMPT=0
      PERTVALUE=0.9D0
      PERTMAX=2.0D0*PERTVALUE
      PERTMIN=0.5D0*PERTVALUE
      EDIFFTOL=1.0D-8
      GEOMDIFFTOL=1.0D-1
      IDIFFTOL=1.0D-3
      DIRECTION='AB'
      ENSEMBLE='T'
      TWOD=.FALSE.
      BULKT=.FALSE.
      ANGLEAXIS=.FALSE.
      NUMRBTYPES = 0
      MAXRBTYPES = 0

      AMBERT=.FALSE.
      AMHT=.FALSE.
      AMHQT=.FALSE.
      AMHQCONTT=.FALSE.
      AMHRMSDT=.FALSE.
      AMHRELQT=.FALSE.
      AMH_RELCOT=.FALSE.
      AMHALLATOMMINT=.FALSE.
      AMHALLATOMTST=.FALSE.
      NOFRQS=.FALSE.

      JPN=1
      NNODES=1
      NATOMS=-1
      CHARMMT=.FALSE.
      STARTFROMPATH=.FALSE.
      READMINT=.FALSE.
      ADDPATH=.FALSE.
      KMCCOMMITT=.FALSE.
      GTT=.FALSE.
      NGTT=.FALSE.
      NGTDISCONNECTALL=.FALSE.
      NGTSWITCH=0.3D0
      NGTSIZE=11000
      NGTCRSWITCH=2.0D0
C     Sem: begin GT2 controls
      GT2T=.FALSE.
      GT2SPARSE=.TRUE.
      GT2SWITCH=.TRUE.
      GT2DisconnectSources=.TRUE.
      GT2ALTPBB=.TRUE.
      GT2RESCALE=.FALSE.
      GT2NORMALISE=.FALSE.
      GT2RSWITCH=0.08D0
      GT2PTOL=1.0D-5
C     Sem: end GT2 controls
      GTINT=0
      REGROUPT=.FALSE.
      REGROUPTHRESH=-1.0D100
      REGROUPRATET=.FALSE.
      REGROUPRATETHRESH=1.0D100
      REGROUPPET=.FALSE.
      REGROUPPETHRESH=-1.0D100
      REGROUPFREET=.FALSE.
      REGROUPFREEABT=.FALSE.
      REGROUPFREETHRESH=-1.0D100
      PABCONV=1.0D-8
      OMEGA=1.0D0 ! GAUSS-SEIDEL ITERATION
      KMCT=.FALSE.
      NCONNMIN=0
      MAXBREAK=1.0D6
      NKMCCYCLES=100.0D0
      PAIRTHRESH=1.0D0
      UNRST=.FALSE.
      NINTS=-1
      NCPU=1
      UNSTRING='UNDEFINED'
      ZSYM='UNKNOWN'
      EXEC='UNDEFINED'
      EXECGMIN='UNDEFINED'
      PATHNAME='UNDEFINED'
      NOPOINTS=.FALSE.
      TRIPLES=.FALSE.
      STARTTRIPLES=.FALSE.
      ADDTRIPLES=.FALSE.
      DIJKSTRAT=.FALSE.
      DIJKSTRAWAITT=.FALSE.
      DIJPAIRT=.FALSE.
      BARRIERSORT=.FALSE.
      BARRIERSHORT=.FALSE.
      RATESHORT=.FALSE.
      DIJINITT=.FALSE.
      DIJINITFLYT=.FALSE.
      DIJINITSTARTT=.FALSE.
      DIJINITCONTT=.FALSE.
      TSTHRESH=HUGE(1.0D0)
      MAXBARRIER=HUGE(1.0D0)
      MAXDOWNBARRIER=HUGE(1.0D0)
      COSTFUNCTIONPOWER=1
      EXPCOSTFUNCTION=.FALSE.
      INDEXCOSTFUNCTION=.FALSE.
      COPYFILES=''
      COPYOPTIMT=.FALSE.
      NPFOLD=0
      TFOLDT=.FALSE.
      NTFOLD=1.0D4 ! real not integer !
      TOMEGA=1.0D0
      TFOLDTHRESH=1.0D-3
      MAXATTEMPT=1
      CALCORDERT=.FALSE.
      CONNECTREGIONT=.FALSE.
      CONNECTMIN1=1
      CONNECTMIN2=1
      CONNECTDIST=1.0D10
      SHORTCUTT=.FALSE.
!     PTAUT=.FALSE.
      NPAIRFRQ=0
      MERGEDBT=.FALSE.
      UNTRAPT=.FALSE.
      EUNTRAPTHRESH=1.0D100
      FREEPAIRT=.FALSE.
      PERMDIST=.FALSE.
      PERMISOMER=.FALSE.
      NTAG=0
      TAGT=.FALSE.
      FREEZE=.FALSE.
      NFREEZE=0
      PLANCK=1.0D0
      DUMMYRUNT=.FALSE.
      DUMMYTST=.FALSE.
      REWEIGHTT=.FALSE.
      KSHORTESTPATHST = .FALSE.
      KSHORT_FULL_PRINTT = .FALSE.
      NPATHS = 0
      BHINTERPT=.FALSE.
      BHACCREJ=0.5D0
      BHSTEPSIZE=0.4D0
      BHCONV=0.01D0
      BHSTEPS=1000
      BHTEMP=1.0D0
      BHK=1.0D0
      ICINTERPT=.FALSE.

      USEPAIRST=.FALSE.
      LOWESTFRQT=.FALSE.
      IMFRQT=.FALSE.
      EVCUT=2.0D-6

      BISECTT=.FALSE.
      DIAGT=.FALSE.
      ARNOLDIT=.FALSE.
      SLURMT=.FALSE.
      CVT=.FALSE.
      DOST=.FALSE.
      CHECKCONNECTIONST=.FALSE.
      CLOSEFILEST=.FALSE.
      PULLT=.FALSE.
!     DC430 >
      RBAAT  = .FALSE.
      RBSYMT = .FALSE.
      FRICTIONT=.FALSE.
      GAMMAFRICTION=0.0D0
      REMOVEUNCONNECTEDT=.FALSE.
      UNCONNECTEDS='AB'
!
! Constraint potential
!
      INTERPCOSTFUNCTION=.FALSE.
      INTCONSTRAINTT=.FALSE.
      INTCONSTRAINTTOL=0.1D0
      INTCONSTRAINTDEL=1.0D5
      INTCONSTRAINTREP=1.0D0
      INTCONSTRAINREPCUT=20.0D0
      INTREPSEP=0
      INTCONSEP=10000
      CHECKCONINT=.FALSE.
      MAXCONUSE=3
!op226 </init>
!op226 }}}
!op226 Open and read pathdata {{{
!op226 <read>
C
C SAT: Now, there it is! Unit 5 is a standard input unit, and we must not use it!
C 5--->959
C
      OPEN (959,FILE='pathdata',STATUS='OLD')
C
C  Read from unit 959
C
      IR=959
190   CALL INPUT(END)
      IF (.NOT. END) THEN
        CALL READU(WORD)
      ENDIF

      IF (END .OR. WORD .EQ. 'STOP') THEN
        CLOSE(IR)
        RETURN
      ENDIF

      IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'
     &                          .OR. WORD .EQ. '\\') THEN 
         GOTO 190

! </read> 
!op226 }}}

!op226 keyword loop {{{
! <kwd>
C
C  Add the minima and transition= states in path.info.<PATHNAME>
C  to an existing database. The end points are NOT assumed to
C  belong to the A and B sets.
C

include pathsample.keyword.addpath.inc.f90
      ELSE IF (WORD.EQ.'ADDPATH') THEN
         ADDPATH=.TRUE.
         CALL READA(PATHNAME)
C
C  ADDPT determines whether we call ADDPERM to add permutational isomers of every stationary point to the
C  min.data and ts.data databases. Speeds up 2DLJ7! For other tagged atom situations the number of isomers
C  will equal the number of atoms for C1 symmetry, so this is probably a bad idea.
C
      ELSE IF (WORD.EQ.'ADDPERM') THEN
         ADDPT=.TRUE.
C
C  STARTTRIPLES and ADDTRIPLES specify the triples format for
C  path.info files read using the STARTFROM and ADDPATH keywords.
C
      ELSE IF (WORD.EQ.'ADDTRIPLES') THEN
         ADDTRIPLES=.TRUE.
C
C  Use the angle-axis system for rigid bodies

      ELSE IF (WORD.EQ.'ANGLEAXIS') THEN
         ANGLEAXIS=.TRUE.
      ELSE IF (WORD.EQ.'ANGLEAXIS2') THEN
         ANGLEAXIS2=.TRUE.
C
C  sf344> Set AMBER9 potential.
C
       ELSE IF (WORD.EQ.'AMBER9') THEN
          AMBERT=.TRUE.
       ELSE IF (WORD.EQ.'NAB') THEN
          AMBERT=.TRUE.
C
C  Set AMH potential.
C
       ELSE IF (WORD.EQ.'AMH') THEN
          AMHT=.TRUE.

C
C  Calculated  Q between minima and native .
C
      ELSE IF (WORD.EQ.'AMHQ') THEN
         AMHQT=.TRUE.
         CALL READI(WHICHMIN)
          PRINT '(A,I6,A,I6)','keywords> Calculate AMH Q  ',WHICHMIN
         IF (NITEMS.GT.2) THEN
          PRINT '(A)','keywords> ERROR - AMHQ '
          STOP
         ENDIF
C
C  Calculated  Q_contact  between minima and native .
C
      ELSE IF (WORD.EQ.'AMHQCONT') THEN
         AMHQCONTT=.TRUE.
         CALL READI(WHICHMIN)
         CALL READF(QCONTCUT)
          PRINT '(A,I6)','keywords> Calculate AMH Q CONTACT ',WHICHMIN
          PRINT '(A,G3.2)','keywords> Within a distance cutoff ',QCONTCUT
         IF (NITEMS.GT.3) THEN
          PRINT '(A)','keywords> ERROR - AMH Q CONTACT '
          STOP
         ENDIF
C
C  Calculated  RMSD between minima and native .
C
      ELSE IF (WORD.EQ.'AMHRMSD') THEN
         AMHRMSDT=.TRUE.
         CALL READI(WHICHMIN)
          PRINT '(A,I6,A,I6)','keywords> Calculate AMH RMSD  ',WHICHMIN
         IF (NITEMS.GT.2) THEN
          PRINT '(A)','keywords> ERROR - AMHRMSD '
          STOP
         ENDIF
C
C  Calculated Relative Qs between minima.
C
      ELSE IF (WORD.EQ.'AMHRELQ') THEN
         AMHRELQT=.TRUE.
         CALL READI(QRELONE)
         CALL READI(QRELTWO)
          PRINT '(A,I6,A,I6)','keywords> AMHRELQ min 1 ',QRELONE,' min 2',QRELTWO
         IF (NITEMS.GT.3) THEN
          PRINT '(A)','keywords> ERROR - AMHRELQ GT 2'
          STOP
         ENDIF

C
C  Calculated Relative Contact Order
C
       ELSE IF (WORD.EQ.'AMH_RELCO') THEN
         AMH_RELCOT=.TRUE.
         CALL READI(WHICHMIN)
         CALL READF(RELCOCUT)
          PRINT '(A,I6)','keywords> Calculate AMH Relative Contact Order ',WHICHMIN
          PRINT '(A,G20.10)','keywords> Within a distance cutoff ',RELCOCUT

         IF (NITEMS.GT.3) THEN
          PRINT '(A)','keywords> ERROR - RELCO GT 2'
          STOP
         ENDIF

C
C  Output AMH all-atom MIN .
C
      ELSE IF (WORD.EQ.'AMHALLATOMMIN') THEN
         AMHALLATOMMINT=.TRUE.
         PRINT '(A,I6,A,I6)','keywords> output AMH all atom MIN structures'
C
C  Output AMH all-atom TS .
C
      ELSE IF (WORD.EQ.'AMHALLATOMTS') THEN
         AMHALLATOMTST=.TRUE.
         PRINT '(A,I6,A,I6)','keywords> output AMH all atom TS structures'

C
C  Rates from eigenvalues obtained with weighted Arnoldi subspace method.
C
      ELSE IF (WORD.EQ.'ARNOLDI') THEN
         ARNOLDIT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
C
C  Parameters for OPTIM basin-hopping interpolation jobs.
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
         IF (NITEMS.GT.10) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ICINTERP') ICINTERPT=.TRUE.
C
C  Parameters for OPTIM bisection interpolation jobs.
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
C  Bulk with periodic boundary conditions.
C
      ELSE IF (WORD.EQ.'BULK') THEN
         BULKT=.TRUE.
         CALL READF(BOXLX)
         CALL READF(BOXLY)
         CALL READF(BOXLZ)
C
C  Calculate order parameter for all stationary points and exit.
C  The order parameter routine will be system specific and must
C  replace the dummy routine provided.
C  Minima and transition states will be listed in separate files
C  for points with order parameters above/below the specified threshold,
C  with default value 0.0.
C
      ELSE IF (WORD.EQ.'CALCORDER') THEN
         CALCORDERT=.TRUE.
C
C Capsomere definition
C
      ELSE IF (WORD.EQ.'CAPSID') THEN
         ANGLEAXIS=.TRUE.
         RIGIDBODY=.TRUE.
         IF (.NOT.ALLOCATED(CAPSOMERDEFS)) THEN
            MAXRBTYPES = 10
            ALLOCATE(CAPSOMERDEFS(MAXRBTYPES))

         ELSE IF (NUMRBTYPES.EQ.MAXRBTYPES) THEN
            ALLOCATE(TEMPCAPSOMERDEFS(MAXRBTYPES))
            TEMPCAPSOMERDEFS = CAPSOMERDEFS

            MAXRBTYPES = MAXRBTYPES + 10
            DEALLOCATE(CAPSOMERDEFS)
            ALLOCATE(CAPSOMERDEFS(MAXRBTYPES))

            CAPSOMERDEFS(1:NUMRBTYPES) = TEMPCAPSOMERDEFS
            DEALLOCATE(TEMPCAPSOMERDEFS)
         ENDIF

         NUMRBTYPES = NUMRBTYPES + 1

         CAPSOMERDEFS(NUMRBTYPES)%NBASALSITES = 5

         CALL READF(CAPSOMERDEFS(NUMRBTYPES)%RHO)
         CALL READF(CAPSOMERDEFS(NUMRBTYPES)%EPSILONREP)
         CALL READF(CAPSOMERDEFS(NUMRBTYPES)%RADIUS)
         IF (NITEMS.GT.4) THEN
            CALL READF(CAPSOMERDEFS(NUMRBTYPES)%HEIGHT)
            IF (NITEMS.GT.5) CALL READI(CAPSOMERDEFS(NUMRBTYPES)%NBASALSITES)
         ELSE
            CAPSOMERDEFS(NUMRBTYPES)%HEIGHT = 0.5D0*CAPSOMERDEFS(NUMRBTYPES)%RADIUS
         ENDIF
C
C  Set CHARMM potential. NDIHE is the number of dihedral angles in the system, which
C is used in perturb and tssearch to perturb the system randomly
C
      ELSE IF (WORD.EQ.'CHARMM') THEN
         CHARMMT=.TRUE.
         CALL READI(NDIHE)
C
C  Check minimum number of connections for each minimum at startup.
C
      ELSE IF (WORD.EQ.'CHECKCONNECTIONS') THEN
         CHECKCONNECTIONST=.TRUE.
C
C  Check for internal minimum in constraint terms for INTCONSTRAINT
C
      ELSE IF (WORD.EQ.'CONINT') THEN
         CHECKCONINT=.TRUE.
C
C  Close and open all files as needed. The objective is to work around
C  the misinteraction of nfs with the Linux kernel, which causes random
C  control characters to be written to open nfs mounted files.
C
      ELSE IF (WORD.EQ.'CLOSEFILES') THEN
         CLOSEFILEST=.TRUE.
C
C  Minimum number of connections for each minimum.
C
      ELSE IF (WORD.EQ.'CONNECTIONS') THEN
         CALL READI(CONNECTIONS)
         IF (NITEMS.GT.2) CALL READI(MAXTSATTEMPTS)
C
C  Specify that the database is grown by attempting all connections between
C  known minima within a distance cutoff.
C
      ELSE IF (WORD.EQ.'CONNECTREGION') THEN
         CONNECTREGIONT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(CONNECTMIN1)
         IF (NITEMS.GT.2) CALL READI(CONNECTMIN2)
         IF (NITEMS.GT.3) CALL READF(CONNECTDIST)
C
C  Specify files to be copied to distributed nodes in addition to the default
C  odata.<pid> and finish.<pid> files.
C
      ELSE IF (WORD.EQ.'COPYFILES') THEN
C        PRINT '(A,I8)','NITEMS=',NITEMS
         DO J1=2,NITEMS
            CALL READA(WW)
            WRITE(COPYFILES,'(3A)') TRIM(ADJUSTL(COPYFILES)),' ',TRIM(ADJUSTL(WW))
C           PRINT '(3A)','WW,COPYFILES: ',WW,COPYFILES
         ENDDO
C
C  All OPTIM output files will be copied from distributed nodes and not deleted. 
C  This also occurs if DEBUG is .TRUE.
C
      ELSE IF (WORD.EQ.'COPYOPTIM') THEN
         COPYOPTIMT=.TRUE.
C
C  Use NCPU's by starting up to NCPU's OPTIM jobs.
C
      ELSE IF (WORD.EQ.'CPUS') THEN
         CALL READI(NCPU)
      ELSE IF (WORD.EQ.'CV') THEN
         CVT=.TRUE.
         CVTMIN=1.0D0
         CVTMAX=2.0D0
         CVTINC=0.1D0
         CALL READF(CVTMIN)
         CALL READF(CVTMAX)
         CALL READF(CVTINC)
C
C  "JOBSPERNODE JPN" specifies parallel run with JPN jobs per node. Mutually exclusive with CPUS keyword.
C  This requests to read node file. Number of nodes and nodefile must be
C  specified on the command line for this to work. If no command line args are
C  given "JOBSPERNODE JPN" acts same as "CPUS NCPUS".
C
C
C  Maximum number of connect cycles in CYCLE.
C
      ELSE IF (WORD.EQ.'CYCLES') THEN
         CALL READI(NATTEMPT)

      ELSE IF (WORD.EQ.'DB') THEN
! Set the EFIELDT separately, if required
         DBPT = .TRUE.
         CALL READF(DBSIGBB)
         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFDUM(DBSIGBB)

      ELSE IF (WORD.EQ.'DBTD') THEN
! Set the EFIELDT separately, if required
         DBPTDT = .TRUE.
         CALL READF(DBSIGBB)
         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFDUM(DBSIGBB)
C
C  Turn on debug printing.
C
      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.
C
C  Rates from direct diagonalisation.
C
      ELSE IF (WORD.EQ.'DIAG') THEN
         DIAGT=.TRUE.
         NDIAGEIG=10
         DIAGSCALE=1.0D0
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) CALL READI(NDIAGEIG)
         IF (NITEMS.GT.3) CALL READF(DIAGSCALE)
C
C  DIJINIT specifies a Dijkstra analysis to try and construct an initial path.
C  We only need files start and finish containing the end point coordinates for
C  DIJINITSTART - this keyword creates new min.A, min.B, points.min, points.ts,
C  min.data and ts.data files, so it is potentially dangerous!
C  DIJINITCONT specifies an initial path run where the above files already exist.
C  In a DIJINITSTART run C  file min.A is created exist with entries 1 1 and file min.B
C  with entries 1 2. min.data is set up with the two entries for
C  these minima. The coordinates are obtained via odata.start and odata.finish and put in records
C  1 and 2 in the points.min file. 
C
      ELSE IF (WORD.EQ.'DIJINITCONT') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(I20)') COSTFUNCTIONPOWER
               ENDIF
         ENDIF
         DIJINITCONTT=.TRUE.
         DIJINITT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'DIJINITSTART') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(i20)') CostFunctionPower
               ENDIF
         ENDIF
         DIJINITSTARTT=.TRUE.
         DIJINITT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'DIJINITCONTFLY') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(I20)') COSTFUNCTIONPOWER
               ENDIF
         ENDIF
         DIJINITCONTT=.TRUE.
         DIJINITFLYT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'DIJINITSTARTFLY') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(i20)') CostFunctionPower
               ENDIF
         ENDIF
         DIJINITSTARTT=.TRUE.
         DIJINITFLYT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
C
C  Dijkstra analysis to find the largest contribution to the SS rate constant
C  from each B (or A) minimum.
C
      ELSE IF (WORD.EQ.'DIJKSTRA') THEN
         DIJKSTRAT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
C
C  Dijkstra analysis to find the path with the lowest sum of waiting times
C  times initial conditional probability.
C
      ELSE IF (WORD.EQ.'DIJKSTRAWAIT') THEN
         DIJKSTRAWAITT=.TRUE.
         DIJKSTRAT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
C
C  DIJPAIR specifies a Dijkstra analysis to propose the pairs of minima to
C  try and connect during the database construction.
C
      ELSE IF (WORD.EQ.'DIJPAIR') THEN
         CALL READI(MAXATTEMPT)
         IF (NITEMS.GT.2) CALL READI(NCONNMIN)
         IF (NITEMS.GT.3) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BARRIER') BARRIERSORT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
C
C  Direction of sampling.
C
      ELSE IF (WORD.EQ.'DIRECTION') THEN
         CALL READA(DIRECTION)
C
C  Docking routine to calculate the binding free energy of a ligand to a protein
C
      ELSE IF (WORD.EQ.'DOCK') THEN
         DOCKT=.TRUE.
         CALL READI(PARALLEL)
         IF (NITEMS.GT.2) THEN
          DO J1=1,5
            CALL READI(NDUMMY)
            IF(NDUMMY==0) DSTAGE(J1)=.FALSE.
          END DO
         END IF
      ELSE IF (WORD.EQ.'DOS') THEN
         DOST=.TRUE.
         DOSEMIN=1.0D0
         DOSEMAX=2.0D0
         DOSEINC=0.1D0
         CALL READF(DOSEMIN)
         CALL READF(DOSEMAX)
         CALL READF(DOSEINC)
C
C  DSCALE is the distance scaling for decisions about whether to connect a given
C  pair of local minima. PSCALE is the scaling for the difference in GPFOLD.
C  IF (EXP(-(DISTANCE-DSCALE)/DSCALE).LT.RANDOM) is used in cycle.
C
      ELSE IF (WORD.EQ.'DSCALE') THEN
         CALL READF(DSCALE)
C
C  Setting DUMMYRUN to true means that no OPTIM jobs are actually submitted.
C
      ELSE IF (WORD.EQ.'DUMMYRUN') THEN
         DUMMYRUNT=.TRUE.
C
C  Setting DUMMYTST to true means that we create phoney entries in ts.data
C  corresponding to nearest-neighbour minima. For use with BHINTERP and BISECTT.
C
      ELSE IF (WORD.EQ.'DUMMYTS') THEN
         DUMMYTST=.TRUE.
C
C  DUMPGROUPS specifies that the groups of potential energy minima and transition
C  states obtained from subroutine regroupfree2 should be saved in files
C  minima_groups and ts_groups
C
      ELSE IF (WORD.EQ.'DUMPGROUPS') THEN
         DUMPGROUPST=.TRUE.
C
C  Energy difference  criterion for distinguishing stationary points
C  Can also be specified with ETOL
C
      ELSE IF (WORD.EQ.'EDIFFTOL') THEN
         CALL READF(EDIFFTOL)

!  Turns on the electric field   
      ELSE IF (WORD.EQ.'EFLD') THEN
         EFIELDT = .TRUE.
C
      ELSE IF (WORD.EQ.'ENERGY') THEN
         CALL READF(TOTALE)
         ENSEMBLE='E'
C
C  Energy difference  criterion for distinguishing stationary points
C
      ELSE IF (WORD.EQ.'ETOL') THEN
         CALL READF(EDIFFTOL)
C
C  Threshold for zero eigenvalues used to skip triples in GETALLPATHS
C
      ELSE IF (WORD.EQ.'EVCUT') THEN
         CALL READF(EVCUT)
C
C  OPTIM executable.
C
      ELSE IF (WORD.EQ.'EXEC') THEN
         CALL READA(EXEC)
C
C  GMIN executable.
C
      ELSE IF (WORD.EQ.'EXECGMIN') THEN
         CALL READA(EXECGMIN)
C
C  Write the coordinates of minimum WHICHMIN to file extractedmin and stop.
C
      ELSE IF (WORD.EQ.'EXTRACTMIN') THEN
         CALL READI(WHICHMIN)
         EXTRACTMINT=.TRUE.
C
C  Write the coordinates of minimum WHICHTS to file extractedts and stop.
C
      ELSE IF (WORD.EQ.'EXTRACTTS') THEN
         EXTRACTTST=.TRUE.
         CALL READI(WHICHTS)
C
C Choose connection pairs based on free energy barriers between regrouped minima
C and product minima.
C
      ELSE IF (WORD.EQ.'FREEPAIRS') THEN
         FREEPAIRT=.TRUE.
         CALL READF(REGROUPFREETHRESH)
         CALL READF(FREETHRESH)
         CALL READF(EINC)
         IF (EINC.LE.0.0D0) THEN
            PRINT '(2(A,G20.10))','keywords> WARNING EINC for FREEPAIRS reset from ',EINC,' to ',1.0D0
            EINC=1.0D0
         ENDIF
C
C  Frozen atoms (to adjust zeros when reading frequencies)
C  and needed for freezing atoms is tssearch.
C
      ELSE IF (WORD.EQ.'FREEZE') THEN
         IF (NATOMS.LE.0) THEN
            PRINT '(A,I6,A)','keywords> ERROR - NATOMS=',NATOMS,' NATOMS keyword must preceed FREEZE'
            STOP
         ENDIF
         IF (.NOT.ALLOCATED(FROZEN)) THEN
            ALLOCATE(FROZEN(NATOMS))
            DO J1=1,NATOMS
               FROZEN(J1)=.FALSE.
            ENDDO
         ENDIF
         FREEZE=.TRUE.
         DO J1=1,NITEMS-1
            NFREEZE=NFREEZE+1
            CALL READI(NDUMMY)
            FROZEN(NDUMMY)=.TRUE.
         ENDDO
         IF (PERMDIST) THEN
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
C Friction coefficient for Berezhkovsii, Pollak, Zitserman formulation
C JCP, 97, 2422, 1992
C
      ELSE IF (WORD.EQ.'FRICTION') THEN
         FRICTIONT=.TRUE.
         IMFRQT=.TRUE.
         CALL READF(GAMMAFRICTION)
C
C  Distance criterion for distinguishing stationary points
C
      ELSE IF (WORD.EQ.'GEOMDIFFTOL') THEN
         CALL READF(GEOMDIFFTOL)
C
C  Specify graph transformation rate calculation. NCONNMIN is the connectivity
C  at which (and below) minima are removed. GTINT is the interval for performing
C  the analysis during a DPS calculation. IF GTINT=1 we do the analysis every cycle,
C  consisting of NCPU jobs, if GTINT=2 we do the analysis every other cycle. 
C  GTINT is zero by default, in which case we only do the analysis if NATTEMPTS=0.
C  GT is the original DJW implementation where we renormalise out all intervening
C  miinima. The resulting branching probablities and waiting times do not allow
C  for return to the starting minimum - they are like the rejection-free BKL approach.
C
      ELSE IF (WORD.EQ.'GT') THEN
         GTT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) CALL READI(GTINT)
C
C  Semen Trygubenko, last modification Thu Mar 15 12:52:36 GMT 2007
C
C  A different implementation of GT method, as described in our GT JCP paper.
C  Input format: {SGT|DGT|SDGT} [DisconnectSources [AltPbb [Rescale [Normalise]]]]
C  Example input line: SDGT F T F F
C 
C  If set to true DisconnectSources instructs to disconnect sources, yielding rate
C  that theoretically corresponds to the one obtained stochastically with KMC. Disconnection
C  of the sources is currently implemented in DGT part of GT. For SGT runs with DisconnectSources=T
C  a call is made to DGT once all SGT work is done. DisconnectSources=T can be memory hungry.
C
C  If set to true, GT2Sparse instructs to perform GT analysis using
C  sparse-optimised algorithms and data structures.
C  
C  GT2Switch=T instructs the program to change the data structures and algorithms
C  from sparse to dense when the density requirement (can be adjusted by GT2RSWITCH) is
C  met.
C  
C  GT2AltPbb, when set to true, triggers the evaluation of Pbb sums using the algorithm
C  suggested by David to maintain precision. Alternative remedy is provided with
C  GT2Rescale which I have devised to stop the errors propagation. Both of these
C  should double the execution time, and are recognized by both
C  SGTDetachNode and DGTDetachNode routines. NB: Disconnect routine, another place where
C  roundoff errors can breed like rabbits, does not recognize either of these as of
C  Thu Mar 15 12:45:55 GMT 2007.
C  
C  GT2Normalise=T will instruct GT2input.f90 to normalize the branching probabilities when
C  obtained from PATHSAMPLE. This is purely for debugging purposes, as they always should be.
C 
      ELSE IF (word=='SGT'.or.word=='DGT'.or.word=='SDGT'.or.word=='GT2') then
         GT2T=.TRUE.
         if (word=='SGT') then
              GT2Sparse=.True.
              GT2Switch=.False.
         else if (word=='DGT') then
              GT2Sparse=.False.
              GT2Switch=.False.
         else if (word=='SDGT') then
              GT2Sparse=.True.
              GT2Switch=.True.
         else
              print '(A)','Please specify SGT, DGT or SDGT. GT2 keyword is now obsolete';stop
         endif
         PRINT *,'keywords> GT2Sparse=',GT2Sparse
         PRINT *,'keywords> GT2Switch=',GT2Switch
    
         IF (NITEMS>1) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2DisconnectSources
                   PRINT '(1x,A,L5)','keywords> GT2DisconnectSources=',GT2DisconnectSources
              ELSE
                   PRINT '(1x,A,L5)','keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>2) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2AltPbb
                   PRINT '(1x,A,L5)','keywords> GT2AltPbb=',GT2AltPbb
              ELSE
                   PRINT '(1x,A,L5)','keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>3) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2Rescale
                   PRINT *,'keywords> GT2Rescale=',GT2Rescale
              ELSE
                   PRINT *,'keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>4) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2Normalise
                   PRINT *,'keywords> GT2Normalise=',GT2Normalise
              ELSE
                   PRINT *,'keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>5) THEN
              print '(1x,a)','Input error: more than 4 entities following SGT, DGT or SDGT keyword!';stop
         ENDIF
C
C Switching ratio. Specifies when to switch from SGT to DGT.
C See GT2.f90 for details.
C
      ELSE IF (WORD=='GT2RSWITCH') then
         CALL READF(GT2RSWITCH)
C
C Is used to establish whether a node is a dead end or not.
C See GT2.f90 for details.
C
      ELSE IF (WORD=='GT2PTOL') then
         CALL READF(GT2PTOL)
C
C Prints the negative eigenvalue of each TS in ts.data as final (ninth) column
C
      ELSE IF (WORD=='IMFRQ') THEN
         IMFRQT=.TRUE.
C
C  Use constraint potential for interpolation as a connection metric (instead of distance).
C
      ELSE IF (WORD.EQ.'INTCONSTRAINT') THEN
         INTCONSTRAINTT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(INTCONSTRAINTTOL)
         IF (NITEMS.GT.2) CALL READF(INTCONSTRAINTDEL)
         IF (NITEMS.GT.3) CALL READF(INTCONSTRAINTREP)
         IF (NITEMS.GT.4) CALL READF(INTCONSTRAINREPCUT)
         IF (NITEMS.GT.5) CALL READI(INTCONSEP)
         IF (NITEMS.GT.6) CALL READI(INTREPSEP)
         INTERPCOSTFUNCTION=.TRUE.
C
C  Use interpolation potential for LJ.
C
      ELSE IF (WORD.EQ.'INTLJ') THEN
         INTLJT=.TRUE. 
         IF (NITEMS.GT.1) CALL READF(INTLJDEL)
         IF (NITEMS.GT.2) CALL READF(INTLJEPS)
         INTERPCOSTFUNCTION=.TRUE.
C
C  Inertia difference criterion - no longer used for distinguishing stationary points!
C
      ELSE IF (WORD.EQ.'ITOL') THEN
         CALL READF(IDIFFTOL)
C
C  Number of jobs to run per node for a distributed memory architecture.
C
      ELSE IF (WORD.EQ.'JOBSPERNODE') THEN
         CALL READI(JPN)
         CALL GETNODES(NCPU)
C
C  Semen Trygubenko, Thu Mar 15 16:21:41 GMT 2007
C  Minimum number of connections a minimum ought to have to be included in the
C  rate calculation or regrouping schemes. Can be set directly for some
C  keywords, but not others, so this keyword is provided just to make sure
C  it casn be set for all methods.
C
      ELSE IF (WORD.EQ.'NCONNMIN') THEN
         if (NITEMS.GT.1) then
              CALL READI(NCONNMIN)
         else
              print *, 'keywords> Usage: RateNConnMin <integer>'; stop
         endif
      ELSE IF (WORD.EQ.'KMC') THEN
         KMCT=.TRUE.
         NOPOINTS=.TRUE.
         IF (NITEMS.GT.1) CALL READF(NKMCCYCLES)
         IF (NITEMS.GT.2) CALL READF(PAIRTHRESH)
         IF (NITEMS.GT.3) CALL READI(NCONNMIN)
C
C  Specify new graph transformation rate calculation. NCONNMIN is the connectivity
C  at which (and below) minima are removed. 
C  NGT is the new implementation where we renormalise out all intervening
C  miinima but allow return to the starting state. This should give committor
C  probabilities if we remove all the I minima, unlike GT and GT2.
C
      ELSE IF (WORD.EQ.'NGT') THEN
         NGTT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.OR.LOYNO=='F') then
                 READ(LOYNO,'(l1)') NGTDISCONNECTALL
                 IF (NGTDISCONNECTALL) PRINT '(A)','keywords> NGT will calculate kSS, kNSS and kKMC'
                 IF (.NOT.NGTDISCONNECTALL) PRINT '(A)','keywords> NGT will calculate kSS and kNSS but not kKMC'
              ELSE
                 PRINT '(1X,A,L5)','keywords> Invalid after "NGT" keyword'; STOP
              ENDIF
         ENDIF
         IF (NITEMS.GT.3) CALL READF(NGTSWITCH)
         IF (NITEMS.GT.4) CALL READI(NGTSIZE)
         PRINT '(A,F12.4)','keywords> NGT will switch to dense renormalisation scheme at threshold ',NGTSWITCH
         PRINT '(A,I6)','keywords> NGT maximum square matrix size in dense phase=',NGTSIZE
         IF (NITEMS.GT.5) CALL READF(NGTCRSWITCH) ! N.B. threshold in giga-bytes
         PRINT '(A,F12.4)','keywords> NGT will use compressed row storage scheme beyond threshold ',NGTCRSWITCH
C
C  Set KMC parameters: number of KMC runs for averages and value of PAIRTHRESH.
C  If the product of branching probabilities p12*p21 > PAIRTHRESH then we
C  renormalise this pair.
C
      ELSE IF (WORD.EQ.'KMCCOMMIT') THEN
         NOPOINTS=.TRUE.
         KMCCOMMITT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(NKMCCYCLES)
         IF (NITEMS.GT.2) CALL READF(MAXBREAK)
         IF (NITEMS.GT.3) CALL READF(PABCONV)
         IF (NITEMS.GT.4) CALL READF(OMEGA)
         IF (NITEMS.GT.5) CALL READF(PAIRTHRESH)
         IF (NITEMS.GT.6) CALL READI(NCONNMIN)
C
C  k-th shortest paths analysis for each B (or A) minimum.
C
      ELSE IF (WORD.EQ.'KSHORTESTPATHS') THEN
         KSHORTESTPATHST=.TRUE.
         CALL READI(NPATHS)
         CALL READI(NCONNMIN)
         IF (NITEMS.GT.3) THEN 
            CALL READA(LOYNO)
            IF (LOYNO == 'T') KSHORT_FULL_PRINTT = .TRUE.
         ENDIF
C
C  Whether to read extra curvatures from min.data.info files in DUMMYTS runs
C
      ELSE IF (WORD.EQ.'LOWESTFRQ') THEN
         LOWESTFRQT=.TRUE.

      ELSE IF (WORD.EQ.'MACHINE') THEN
         MACHINE=.TRUE.
C
C  MAXBARRIER requires both sides to be greater than MAXBARRIER to discard.
C
      ELSE IF (WORD.EQ.'MAXBARRIER') THEN
         CALL READF(MAXBARRIER)
C
C  The maximum number of constraints to use in the constrained potential.
C  The default is 3.
C
      ELSE IF (WORD.EQ.'MAXCON') THEN
         CALL READI(MAXCONUSE)
C
C  MAXDOWNBARRIER checks the downhill barrier in Dijkstra.
C
      ELSE IF (WORD.EQ.'MAXDOWNBARRIER') THEN
         CALL READF(MAXDOWNBARRIER)
C
C  TSTHRESH discards transition states above the specfied threshold. May be useful
C  for producing a better initial path and excluding CHARMM transition states with
C  silly energies from the database. MAXTSENERGY does the same as TSTHRESH keyword 
C  for compatability with OPTIM.
C
      ELSE IF (WORD.EQ.'MAXTSENERGY') THEN
         CALL READF(TSTHRESH)
C
C  Add the minima and transition= states in path.info.<PATHNAME> and
C  output.<PATHNAME> to an existing database. The end points are NOT assumed to
C  belong to the A and B sets.
C
      ELSE IF (WORD.EQ.'MERGEDB') THEN
         MERGEDBT=.TRUE.
         CALL READA(PATHNAME)

      ELSE IF (WORD.EQ.'MSSTOCK') THEN

         MSSTOCKT = .TRUE.
         CALL READI(NRBSITES)
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFMULTSTOCK()

C
C  Number of atoms - essential unless RIGIDBODIES/RBAA keyword is present.
C
      ELSE IF (WORD.EQ.'NATOMS') THEN
         CALL READI(NATOMS)
      ELSE IF (WORD.EQ.'NINTS') THEN
         CALL READI(NINTS)
C
C  If NOFRQS is specified then frequencies are assumed to be absent from path.info files.
C
      ELSE IF (WORD.EQ.'NOFRQS') THEN
         NOFRQS=.TRUE.
C
C  If NOPOINTS is specified then setup should not try to read the min.points or ts.points
C  files. Should be the default for post-DPS database kinetics analysis such as KMC.
C
      ELSE IF (WORD.EQ.'NOPOINTS') THEN
         NOPOINTS=.TRUE.
C
C  Read in an order parameter threshold, ORDERPARAM, that tells us when we are in the other 
C  phase for a DOALLMIN run.
C
      ELSE IF (WORD.EQ.'ORDER') THEN
         CALL READF(ORDERPARAM)
C
C  Set the frequency with which the list of pairs for future searches is recreated.
C  Currently works for default search type based on Pfold difference.
C

      ELSE IF (WORD.EQ.'PAHA') THEN

         CALL READI(PAHID)

         IF (PAHID == 1) THEN
            NRBSITES = 12
         ENDIF

         PAHAT    = .TRUE.
         ALLOCATE(RBSITE(NRBSITES,3))

         IF (PAHID == 1) THEN
            CALL DEFBENZENE()
         ENDIF


      ELSE IF (WORD.EQ.'PAIRLIST') THEN
         CALL READI(NPAIRFRQ)

      ELSE IF (WORD.EQ.'PATCHYD') THEN
         PATCHYDT = .TRUE.
         CALL READI(NRBSITES)
         IF (NRBSITES .NE. 4) THEN
            PRINT *, 'NRBSITES has to be 4'
            STOP 
         ENDIF
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFPATCHYD()
!
!  Whether to optimise the permutational isomers in assessing optimal alignment.
!  For PERMDIST all minimum distances will be minimsed with respect to the
!  specified permutations. For PERMISOMER we only check for permutational isomers
!  if the energy difference is below EDIFFTOL. This should save a great deal
!  of time for large systems containing many equivalent atoms, although the
!  distance won't be minimised with repect to permutations for inequivalent minima.
!
      ELSE IF ((WORD.EQ.'PERMDIST').OR.(WORD.EQ.'PERMISOMER')) THEN
         IF (NATOMS.LE.0) THEN
            PRINT '(A,I6,A)','keywords> ERROR - NATOMS=',NATOMS,' NATOMS keyword must preceed PERMDIST'
            STOP
         ENDIF
         IF (WORD.EQ.'PERMDIST') PERMDIST=.TRUE.
         IF (WORD.EQ.'PERMISOMER') PERMISOMER=.TRUE.
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
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
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
            NDUMMY=1
            DO J1=1,NPERMGROUP
!
!  Sanity checks!
! 
               READ(1,*) NPERMSIZE(J1),NSETS(J1)
               IF (NSETS(J1).GT.3) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1),' is > 3'
                  STOP
               ENDIF
               IF (NDUMMY+NPERMSIZE(J1)-1.GT.3*NATOMS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),
!    &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
            CLOSE(1)
         ELSE
            NPERMGROUP=1 ! ALL ATOMS CAN BE PERMUTED - DEFAULT
            NPERMSIZE(1)=NATOMS ! ALL ATOMS CAN BE PERMUTED - DEFAULT
            DO J1=1,NATOMS
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         PRINT '(A,I6)','keywords> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            PRINT '(A,3(I6,A))',' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
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
!
!  Another check.
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
C
C  Initial PERT parameter for geometry perturbations used in single-ended ts searches
C
      ELSE IF (WORD.EQ.'PERTURB') THEN
         CALL READF(PERTVALUE)
         PERTMIN=PERTVALUE/2.0D0
         PERTMAX=PERTVALUE*2.0D0
         IF (NITEMS.GT.2) CALL READF(PERTMIN)
         IF (NITEMS.GT.3) CALL READF(PERTMAX)
C
C  NPFOLD is the number of iterations per call to the global Pfold subroutine
C  and PFOLDINT is the frequency at which we call this GPFOLD calculation, 
C  i.e. PFOLDINT=1 means for every cycle, PFOLDINT=2 means every other cycle etc.
C
      ELSE IF (WORD.EQ.'PFOLD') THEN
         CALL READI(NPFOLD)
         CALL READI(PFOLDINT)
         CALL READF(OMEGA)
C
C  The value of the Planck constant in units of prevailing energy * seconds. 
C  This is needed for regrouped ts free energies. 
C  For example, for CHARMM we enter the temperature in kcal/mol (i.e. k_B * T),
C  so h needs to be in kcal/mol * seconds i.e. 9.546 * 10^(-14).
C
      ELSE IF (WORD.EQ.'PLANCK') THEN
         CALL READF(PLANCK)
      ELSE IF (WORD.EQ.'PSCALE') THEN
         CALL READF(PSCALE)
C
C  Try connecting minima that have the largest values of equilibrium occupation probability times
C  waiting time to minima on the fastest path.
C  Don't bother - UNTRAP with free energy sorting should do the job.
C
!     ELSE IF (WORD.EQ.'PTAU') THEN
!        PTAUT=.TRUE.
!        DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'PULL') THEN
         PULLT=.TRUE.
         PRINT '(A)','keywords> Constant fulling force with 4 zero eigenvalues'
      ELSE IF (WORD.EQ.'RBAA') THEN
         CALL READI(NATOMS)
         NATOMS=NATOMS*2
         NTSITES  = (NATOMS/2)*NRBSITES
         IF (DBPTDT) NTSITES = (NATOMS/2-1)*NRBSITES + 4
         IF (NRBSITES == 0) THEN
            PRINT *, 'NRBSITES not yet defined'
            STOP
         ENDIF
         RBAAT = .TRUE.
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

      ELSE IF (WORD.EQ.'READMIN') THEN
         READMINT=.TRUE.
         CALL READA(MINNAME)
C
C  Threshold for redefining the A/B/I sets on the basis of a PE superbasin analysis
C  at energy REGROUPTHRESH.
C
      ELSE IF (WORD.EQ.'REGROUP') THEN
         NOPOINTS=.TRUE.
         REGROUPT=.TRUE.
         CALL READF(REGROUPTHRESH)
C
C  Regroup on the basis of free energies (equivalent to intergroup rates).
C  Subsequent rate calculations use the new groups and free energies.
C
      ELSE IF (WORD.EQ.'REGROUPFREE') THEN
         NOPOINTS=.TRUE.
         REGROUPFREET=.TRUE.
         CALL READF(REGROUPFREETHRESH)
C
C  Regroup on the basis of free energies (equivalent to intergroup rates).
C  Subsequent rate calculations use pe stationary points, but the
C  A and B groups are expanded based on the free energy regrouping.
C
      ELSE IF (WORD.EQ.'REGROUPFREEAB') THEN
         NOPOINTS=.TRUE.
         REGROUPFREEABT=.TRUE.
         CALL READF(REGROUPFREETHRESH)
C
C  Threshold for regrouping and calculating free energies
C  on the basis of a superbasin analysis at PE REGROUPPETHRESH
C
      ELSE IF (WORD.EQ.'REGROUPPE') THEN
         NOPOINTS=.TRUE.
         REGROUPPET=.TRUE.
         CALL READF(REGROUPPETHRESH)
C
C  Threshold for regrouping and calculating free energies
C  on the basis of rates.
C
      ELSE IF (WORD.EQ.'REGROUPRATE') THEN
         NOPOINTS=.TRUE.
         REGROUPRATET=.TRUE.
         CALL READF(REGROUPRATETHRESH)
C
C  Remove stationary points specified in file min.remove and ts.remove
C
      ELSE IF (WORD.EQ.'REMOVESP') THEN
         REMOVESP=.TRUE.
C
C  Remove stationary points specified in file min.remove and ts.remove
C
      ELSE IF (WORD.EQ.'REMOVEUNCONNECTED') THEN
         REMOVEUNCONNECTEDT=.TRUE.
          IF (NITEMS.GT.1) CALL READA(UNCONNECTEDS)
C
C  Retain only stationary points specified in file min.retain and ts.retain
C
      ELSE IF (WORD.EQ.'RETAINSP') THEN
         RETAINSP=.TRUE.
C
C  Reweighting for reactant minima to allow stochastic sampling of reactant
C  in a GT calculation.
C
      ELSE IF (WORD.EQ.'REWEIGHT') THEN
         REWEIGHTT=.TRUE.
         CALL READI(NRWBINS)      ! number of bins in probability distribution
         CALL READI(NRWREACTANT)  ! number of reactant minima to choose
         CALL READA(RWENERGYFILE) ! name of file containing quench data
         IF (ALLOCATED(RWPROB)) DEALLOCATE(RWPROB)
         ALLOCATE(RWPROB(NRWBINS))
         OPEN(UNIT=1,FILE=TRIM(ADJUSTL(RWENERGYFILE)),STATUS='OLD')
         RWEMAX=-1.0D100
         RWEMIN=1.0D100
         NDUMMY=0
         DO
            READ(1,*,END=222) DUMMY
            IF (DUMMY.GT.RWEMAX) RWEMAX=DUMMY
            IF (DUMMY.LT.RWEMIN) RWEMIN=DUMMY
            NDUMMY=NDUMMY+1
         ENDDO
222      CONTINUE
         PRINT '(A,I8,2A)','keyword> ',NDUMMY,' energies read from file ',TRIM(ADJUSTL(RWENERGYFILE))
         RWBINWIDTH=(RWEMAX-RWEMIN)/NRWBINS
         REWIND(1)
         RWPROB(1:NRWBINS)=0.0D0
         DO J1=1,NDUMMY
            READ(1,*) DUMMY
            NDUMMY2=INT((DUMMY-RWEMIN-1.0D-10)/RWBINWIDTH) + 1
!           PRINT '(A,3G20.10,I6)','RWEMIN,RWEMAX,DUMMY,NDUMMY2=',RWEMIN,RWEMAX,DUMMY,NDUMMY2
            RWPROB(NDUMMY2)=RWPROB(NDUMMY2)+1
         ENDDO
         CLOSE(1)
         RWPROB(1:NRWBINS)=RWPROB(1:NRWBINS)/NDUMMY
         IF (DEBUG) THEN
            PRINT '(A)','keyword> bin energy ranges and probabilities:'
            DUMMY=0.0D0
            DO J1=1,NRWBINS
               PRINT '(I6,2F20.10,G20.10)',J1,RWEMIN+RWBINWIDTH*(J1-1),RWEMIN+RWBINWIDTH*J1,RWPROB(J1)
               DUMMY=DUMMY+RWPROB(J1)
            ENDDO
            PRINT '(A,G20.10)','keyword> sum of probabilities=',DUMMY
         ENDIF
C
C  Number of rigid bodies - essential unless NATOMS keyword is present.
C  All we should have to do is then set NATOMS equal to twice the number
C  of rigid bodies to get all the dimensions right.
C
      ELSE IF (WORD.EQ.'RIGIDBODIES') THEN
         CALL READI(NATOMS)
         NATOMS=NATOMS*2
C
C  Random number seed.
C
      ELSE IF (WORD.EQ.'SEED') THEN
         CALL READI(ISEED)
C
C  Try connecting closest minima that are further apart than a given
C  separation in terms of steps in the best path.
C
      ELSE IF (WORD.EQ.'SHORTCUT') THEN
         SHORTCUTT=.TRUE.
         DIJKSTRAT=.TRUE.
         CALL READI(MINSEP)
         IF (NITEMS.GT.2) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BARRIER') BARRIERSHORT=.TRUE.
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'RATE') RATESHORT=.TRUE.

      ELSE IF (WORD.EQ.'SIS')THEN
!
! SIS epidemiological model
!
         SIST=.TRUE.
         CALL READI(SMAX)
         CALL READI(IMAX)
         CALL READI(POPSS)
         CALL READF(SISMU)
         CALL READF(SISKAPPA)
         CALL READF(SISBETA)
         ZSYM='C' ! just setting this randomly in case it shouldn't be undefined...
         NATOMS=1 ! Just setting this randomly in case it shouldn't be undefined...
         PRINT '(A,3I6,1X,3G15.10)','keywords> SIS parameters ',SMAX,IMAX,POPSS,SISMU,SISKAPPA,SISBETA

      ELSE IF (WORD.EQ.'ST') THEN

         STOCKAAT = .TRUE.
         NRBSITES = 1
         ALLOCATE(RBSITE(NRBSITES,3))
         RBSITE(1,:) = 0.D0
C
C  Node specification in nodes.info file in slurm format.
C
      ELSE IF (WORD.EQ.'SLURM') THEN
         SLURMT=.TRUE.
         CALL READI(JPN)
         CALL GETNODES(NCPU)
C
C  Make the initial min.A and min.B files using information in path.info.<PATHNAME> and
C  output.<PATHNAME>. This sets up a single A and a single B minimum, which are specified
C  by STARTMINA and STARTMINB. For path.info files in the DUMPALLPATHS format, these
C  will usually not be the first and last entries!
C
      ELSE IF (WORD.EQ.'STARTFROMPATH') THEN
         STARTFROMPATH=.TRUE.
         CALL READA(PATHNAME)
         CALL READI(STARTMINA)
         CALL READI(STARTMINB)
C
C  Similarly, STARTTRIPLES and ADDTRIPLES specify the triples format for
C  path.info files read using the STARTFROM and ADDPATH keywords.
C
      ELSE IF (WORD.EQ.'STARTTRIPLES') THEN
         STARTTRIPLES=.TRUE.
C
C  OPTIM system symbol, e.g. AX, LS, etc.
C
      ELSE IF (WORD.EQ.'SYSTEM') THEN
         CALL READA(ZSYM)
C
C  Number and mass of a tagged atom. 
C
      ELSE IF (WORD.EQ.'TAG') THEN
         IF (NATOMS.LT.1) THEN
            PRINT '(A)','keywords> ERROR - number of atoms must be set before TAG keyword in pathdata'
            STOP
         ENDIF
         IF (.NOT.ALLOCATED(TAGFAC)) THEN
            TAGT=.TRUE.
            ALLOCATE(TAGFAC(NATOMS),TAGNUM(NATOMS))
            TAGFAC(1:NATOMS)=1.0D0
            TAGNUM(1:NATOMS)=0
         ENDIF
         NTAG=NTAG+1
         CALL READI(TAGNUM(NTAG))
         CALL READF(TAGFAC(TAGNUM(NTAG)))
C
C  Canonical temperature.
C
      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         CALL READF(TEMPERATURE)
         ENSEMBLE='T'
C
C  NTFOLD is the maximum number of iterations for the Tfold subroutine.
C  Make it real to prevent overflow.
C
      ELSE IF (WORD.EQ.'TFOLD') THEN
         TFOLDT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) CALL READF(NTFOLD)
         IF (NITEMS.GT.3) CALL READF(TFOLDTHRESH)
         IF (NITEMS.GT.4) CALL READF(TOMEGA)
C
C  If TRIPLES is .TRUE. we read path.info as min-sad-min triples, rather
C  than the traditional min-sad-min-sad-...
C  This keyword must be used if odata.connect uses DUMPALLPATHS rather than
C  DUMPPATH.
C
      ELSE IF (WORD.EQ.'TRIPLES') THEN
         TRIPLES=.TRUE.
C
C  TSTHRESH discards transition states above the specfied threshold. May be useful
C  for producing a better initial path and excluding CHARMM transition states with
C  silly energies from the database.
C
      ELSE IF (WORD.EQ.'TSTHRESH') THEN
         CALL READF(TSTHRESH)
C
C  Two-dimensional flatland.
C
      ELSE IF (WORD.EQ.'TWOD') THEN
         TWOD=.TRUE.
C
C jmc Set unres potential. NDIHE is the number of dihedral angles in the system, which
C is used in perturb and tssearch to perturb the system randomly
C
      ELSE IF (WORD.EQ.'UNRES') THEN
         UNRST=.TRUE.
         CALL READI(NDIHE)
         ZSYM='C'
C
C Choose connection pairs based on pe barriers to minima in the product set
C
      ELSE IF (WORD.EQ.'UNTRAP') THEN
         UNTRAPT=.TRUE.
         CALL READF(EINC)
         CALL READF(EUNTRAPTHRESH)
C
C Choose connection pairs based on the sequence read from file USEPAIRSFILE
C
      ELSE IF (WORD.EQ.'USEPAIRS') THEN
         USEPAIRST=.TRUE.
         CALL READA(USEPAIRSFILE)
         ! </kwd>
         ! }}}
!op226 End section {{{
!op226 <end>
      ELSE

         CALL REPORT('Unrecognized command '//WORD,.TRUE.)
      ENDIF

      CALL FLUSH(6,ISTAT)
      GOTO 190

      RETURN
      END
      !op226 </end>
      !op226 }}}
          ! }}}
      endselect

      RETURN
      END
