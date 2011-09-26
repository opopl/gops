!op226> Initializations <init> {{{
      NPCALL=0
      NSEED=0
      NS=0
      NSSTOP=0
      !NSAVEINTE=0
!      RESIZET=.FALSE.
      !STEPOUT=.FALSE.
      !SUPERSTEP=.FALSE.
      !NSUPER=10
      !SUPSTEP=1.1D0
      !SACCRAT=0.5D0
      !NSACCEPT=100
      !EVSTEPT=.FALSE.
      DO JP=1,1
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
      MYBLNT=.FALSE.
      CHAPERONINT=.FALSE.
      DFTBT=.FALSE.
      SW=.FALSE.
      XMUL=1
      SCT=.FALSE.
      SQUEEZET=.FALSE.
      NVEC=0
!     SQUEEZER=5.0D0
!     SQUEEZED=0.95D0
      DEBUG=.FALSE.
      SEEDT=.FALSE.
      FREEZECORE=.TRUE.
      FREEZE=.FALSE.
      FREEZERES=.FALSE.
      FREEZEALL=.FALSE.
      UNFREEZERES =.FALSE.
! sf344> unfreeze structures at the final quench
      UNFREEZEFINALQ=.FALSE.
      NFREEZE=0
      ALLOCATE(FROZEN(NATOMS))
! csw34> The FROZENRES array is bigger than needed
      ALLOCATE(FROZENRES(NATOMS))
      DO J1=1,NATOMS
         FROZEN(J1)=.FALSE.
         FROZENRES(J1)=.FALSE.
      ENDDO
      FREEZEGROUPTYPE='GT'
      FREEZEGROUPT=.FALSE.
! csw34> DONTMOVE defaults
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
! END of DONTMOVE defaults
      !CHECKCHIRALITY=.TRUE.
      !NOCISTRANS=.TRUE.
      !NOCISTRANSRNA=.FALSE.
      !NOCISTRANSDNA=.FALSE.
      !MINOMEGA=150.D0
      !UACHIRAL=.FALSE.
      !SETCHIRAL=.FALSE.
      !FIELDT=.FALSE.
      !OHT=.FALSE.
      !IHT=.FALSE.
      !TDT=.FALSE.
      !D5HT=.FALSE.
      !CENT=.FALSE.
      !CENTXY=.FALSE.
      !SETCENT=.FALSE.
      !CENTX=0.0D0
      !CENTY=0.0D0
      !CENTZ=0.0D0
      !QUCENTRE=.FALSE.
      !FIXCOM=.FALSE.
      !FIH=0.0D0
      !FTD=0.0D0
      !FD5H=0.0D0
      !TOLD=0.0001D0
      !TOLE=0.0001D0
      !CUTT=.FALSE.
      !PERIODIC=.FALSE.
      !PARAMONOVPBCX=.FALSE.
      !PARAMONOVPBCY=.FALSE.
      !PARAMONOVPBCZ=.FALSE.
      !PARAMONOVCUTOFF=.FALSE.
      !LJSITE=.FALSE.
      !BLJSITE=.FALSE.
      !LJSITECOORDST=.FALSE.
      !LJSITEATTR=.FALSE.
      !NRUNS=0
      !PCUT=1.0D0
      !RADIUS=0.0D0
      !MAXIT=500
      !MAXIT2=500
      !EXPFAC=10.0D0
      !EXPD=1.0D0
      !CQMAX=1.0D-10
      !BQMAX=1.0D-3
      !RHO=6.0D0
      !NACCEPT=50
      !NORESET=.FALSE.
      !TSALLIST=.FALSE.
      !QTSALLIS=0.0D0
      !PARALLELT=.FALSE.
      !TOSI=.FALSE.
      !WELCH=.FALSE.
      !BINARY=.FALSE.
      !SHIFTCUT=.FALSE.
      !FAL=.FALSE.
      !FNI=.FALSE.
!!     AMBER=.FALSE.
      !AMHT=.FALSE.
      !NINT_AMH=1
      !DPARAM=1.0D0
      !FAKEWATER=.FALSE.
      !AMCUT= .FALSE.
      !MGBWATER=.FALSE.
      !BIN=.FALSE.
      !AMBERSEED= .FALSE.
      !FIXT= .FALSE.
      !FIX= .FALSE.
      !CAP= .TRUE.
      !WATERSTEP= .FALSE.
      !QCUTOFF= 1.0D6
      !RCUTOFF= 1.0D6
      !REALQCUTOFF= 1.0D6
      !REALRCUTOFF= 1.0D6
      !RINGROTSCALE=0.0D0
      !TRACKDATAT=.FALSE.
      !PROGRESS=.FALSE.
      !listupdate=20

      !BLJCLUSTER=.FALSE.

      !CHRMMT=.FALSE.
      !ACESOLV=.FALSE.
      !ACEUPSTEP=50
      !CHRIGIDTRANST=.FALSE.
      !CHRIGIDROTT=.FALSE.
      !CHNMAX=0.0D0
      !CHMDT=.FALSE.
      !CHMDFREQ=HUGE(1)
      !CURRENTIMP=0
      !BOXT=.FALSE.
      !SPHERET=.FALSE.
      !RMST=.FALSE.
      !NEWCONFT=.FALSE.
      !INTMINT=.FALSE.
      !DAESTAT=.FALSE.
      !MAKEOLIGOT=.FALSE.
      !MAKEOLIGOSTART=.FALSE.
      !TRANSXYT=.FALSE.
      !ROTZT=.FALSE.
      !NREPEAT=0
      !NFIXSEG=0
      !OHCELLT=.FALSE.

!!  sf344> AMBER stuff
      !AMBERT=.FALSE.
      !AMCHNMAX=0.0D0
      !AMCHPMAX=0.0D0
      !MDSTEPT=.FALSE.
      !DUMPSTRUCTURES=.FALSE.
!! csw34> RANDOMSEED now works for CHARMM also!
      !RANDOMSEEDT=.FALSE.
!! csw34> Dumping structures after every quench
      !DUMPQUT=.FALSE.
!! csw34> Dumping structures after every step (before quenching)
      !DUMPSTEPST=.FALSE.
!! khs26> Dump best structures after every step
      !DUMPBESTT=.FALSE.
!! csw34> Local sampling within distance constraints 
      !LOCALSAMPLET=.FALSE. 
      !ABTHRESH=999.99
      !ACTHRESH=999.99
!! csw34> AMBER interaction energy logical
      !A9INTET=.FALSE.
      !INTERESTORE=.FALSE.
!! csw34> set COLDFUSION flag to .FALSE.
      !COLDFUSION=.FALSE.
!!
!!  sf344> for specifying movable atoms and ligand rotating steptaking moves
!!
      !MOVABLEATOMST=.FALSE.
      !LIGMOVET=.FALSE.
      !LIGROTSCALE=0.0D0
      !LIGCARTSTEP=0.0D0
      !LIGTRANSSTEP=0.0D0
      !LIGMOVEFREQ=1

!!
!!  csw34> rotamer move stuff
!!
      !ROTAMERT=.FALSE.
!!
!!  csw34> some defaults (just in case)
!!
      !ROTMAXCHANGE=1
      !ROTPSELECT=0.2
      !ROTOCCUW=0.004
      !ROTCENTRE=1
      !ROTCUTOFF=999.99
!!
!!  csw34> atom group rotation moves
!!
      !GROUPROTT=.FALSE.
      !NGROUPS=0
      !GROUPROTFREQ=1
      !GROUPOFFSET=0
!!
      !NOPHIPSIT=.FALSE.
      !OMEGAT=.FALSE.

      !OSASAT=.FALSE.
      !RPRO=1.4D0
      !ODIHET=.FALSE.
      !ORGYT=.FALSE.
      !OEINTT=.FALSE.
      !!MON1(1:2)=1
      !!MON2(1:2)=1

      !BSMIN=.FALSE.
      !RKMIN=.FALSE.
      !PERMDIST=.FALSE.
      !PERMOPT=.FALSE.

      !GAMMA=1.0D0
      !TUNNELT=.FALSE.
      
      !TWOD=.FALSE.
      !COMPRESST=.FALSE.

      !MUPDATE=4
      !DGUESS=0.1D0
      !BFGS=.FALSE.
      !LBFGST=.TRUE.
      !CONJG=.FALSE.
      !TNT=.FALSE.
      !TOLB=0.1D0
      !DBRENTT=.FALSE.
      !GUIDECUT=0.0001D0
      !CPMD=.FALSE.
      !DL_POLY=.FALSE.
      !EFAC=0.0D0
      !EAMP=0.01D0
      !FIXD=.FALSE.
      !NHSMOVE=1
      !T12FAC=1.1D0
      !RENORM=.FALSE.
      !NRENORM=10
      !NRENSTUCK=20
      !XMOVERENORM=6.0
      !TRENORM=1.0D0
      !PACHECO=.FALSE.
      !EAMLJT=.FALSE.
      !PBGLUET=.FALSE.
      !EAMALT=.FALSE.
      !ALGLUET=.FALSE.
      !MGGLUET=.FALSE.
      !GUPTAT=.FALSE.
      !FST=.FALSE.
      !WENZEL=.FALSE.
      !RESTART=.FALSE.
      !NEWRESTART=.FALSE.
      !NRELAX=0
      !NMSBSAVE=0
      !AVOID=.FALSE.
      !AVOIDDIST=1.0D0
      !AVOIDRESEEDT=.TRUE.
      !MAXSAVE=10
      !NHSRESTART=0
      !MAXBFGS=0.4D0

      !CAPSID=.FALSE.
      !STRANDT=.FALSE.
      !PAHT=.FALSE.
      !TIP=.FALSE.
      !QUADT=.FALSE.
      !STOCKT=.FALSE.
      !LJCOULT=.FALSE.
      !COULN=0
      !COULQ=0.0D0
      !COULSWAP = 0.0D0
      !COULTEMP = 0.0D0
      !GAYBERNET=.FALSE.
      !PARAMONOVT=.FALSE.
      !ELLIPSOIDT=.FALSE.
      !PYGPERIODICT=.FALSE.
      !LJCAPSIDT=.FALSE.
      !PYBINARYT=.FALSE.
      !MULTISITEPYT=.FALSE.
      !LJGSITET=.FALSE.
      !PYOVERLAPTHRESH=1.0D0
      !LJSITE=.FALSE.
      !SWAPMOVEST=.FALSE.
      !STICKYT=.FALSE.
      !RIGID=.FALSE.
      !TIPID=4
      !HEIGHT=0.5D0
      !OTPT=.FALSE.
      !LJMFT=.FALSE.
      !Q4T=.FALSE.

!!     DC430 >
      !DBPT        = .FALSE.
      !DBPTDT      = .FALSE.
      !DBLPYT      = .FALSE.
      !DMBLMT      = .FALSE.
      !EFIELDT     = .FALSE.
      !GAYBERNEDCT = .FALSE.
      !GBDT        = .FALSE.
      !GBDPT       = .FALSE.
      !GEMT        = .FALSE.
      !LINRODT     = .FALSE.
      !LWOTPT      = .FALSE.
      !MMRSDPT     = .FALSE.
      !MSGBT       = .FALSE. 
      !MSPYGT      = .FALSE.
      !MSTBINT     = .FALSE.
      !MSSTOCKT    = .FALSE.
      !MULTPAHAT   = .FALSE.
      !NCAPT       = .FALSE.
      !NPAHT       = .FALSE.
      !NTIPT       = .FALSE.
      !PAHAT       = .FALSE.
      !PAPT        = .FALSE.
      !PYGT        = .FALSE.
      !PYGDPT      = .FALSE.
      !SILANET     = .FALSE.
      !STOCKAAT    = .FALSE.
      !TDHDT       = .FALSE.
      !WATERDCT    = .FALSE.
      !WATERKZT    = .FALSE.
!|gd351>
      PATCHY = .FALSE.
      ASAOOS = .FALSE.
!<gd351|

      !THRESHOLDT=.FALSE.
      !BSWL=.FALSE.
      !BSPT=.FALSE.
      !MINDENSITYT=.FALSE.
      !BSPTQMAX=1.0D100
      !BSPTQMIN=-1.0D100
      !BSPTRESTART=.FALSE.
      !HISTSMOOTH=.FALSE.
      !NSpline=1
      !FIXBIN=.FALSE.
      !QUENCHFRQ=1
      !NQUENCH=0

      !DECAY=.FALSE.
      !DECAYPARAM=0.0D0
      !COOP=.FALSE.
      !NCOOP=5
      !COOPCUT=1.0D0

      !!UNSTRING='UNDEFINED'
      !BOXSIZE=20.D0
      !SPHERERAD=20.D0
      !NCHENCALLS=0
      !NATBT=.FALSE.
      !SYMMETRIZE=.FALSE.
      !SYMMETRIZECSM=.FALSE.
      !NSYMINTERVAL=10
      !SYMTOL1=0.1D0
      !SYMTOL2=0.1D0
      !SYMTOL3=0.1D0
      !SYMTOL4=0.1D0
      !SYMTOL5=0.1D0
      !NSYMQMAX=20
      !MATDIFF=0.1D0
      !DISTFAC=0.0D0
      MAXERISE=1.0D-10
      MAXEFALL=-HUGE(1.0D0)
      ARMA=0.4D0
      ARMB=0.4D0
      EPSSPHERE=0.0D0

      !BINSTRUCTURES=.FALSE.
      !TETHER=.FALSE.
      !EQUIL=0
      !PTMC=.FALSE.
      !VISITPROP=.FALSE.
      !HWINDOWS=1

      !FixedEndMoveT=.FALSE.
      !PIVOTP=0.0D0
      !SIDECHAINP=0.0D0

      !DIFFRACTT=.FALSE.
      !THOMSONT=.FALSE.
      !GAUSST=.FALSE.
      !MPIT=.FALSE.
      !DUMPINT=1000 ! default is to dump a restart file every 1000 cycles of mc.f
      !RESTORET=.FALSE.
      !DUMPFILE=''
      !INTEDUMPFILE=''
      !MOVESHELLT=.FALSE.
      !SHELLMOVEMAX=0
      !SHELLPROB=0.0D0
      !COREFRAC=0.0D0
      !MYSDT=.FALSE.
      !TSTAR=-1.0D0
      !PRTFRQ=1
      !BSPTDUMPFRQ=100000
      !QUENCHDOS=.FALSE.
      !QDT=.FALSE.
      !QD2T=.FALSE.
      !MULLERBROWNT=.FALSE.
      !QDLIMIT=-1
      !CHARMMENERGIES=.FALSE.
      !FIXDIHEFLAG=.FALSE.
      !JMT=.FALSE.
      !PROJIT=.FALSE.
      !PROJIHT=.FALSE.
      !COLDFUSIONLIMIT=-1.0D6
      !MODEL1T=.FALSE.

      !VGW=.FALSE.              ! VGW PARAMETERS
      !LJSIGMA=1.0D0
      !LJEPSILON=1.0D0
      !TAUMAX=5.0D0
      !TAUMAXFULL=7.0D0
      !CPFACTORSG=3.0D0
      !CPFACTORFG=1.0D0
      !CPS=1
      !CPF=1
      !VGWTOL=0.0001D0
      !ACKLANDT=.FALSE.
      !ACKLANDID=5
      !STEEREDMINT=.FALSE.
      !DF1T=.FALSE.
      !PULLT=.FALSE.
      !CSMT=.FALSE.
      !CSMGUIDET=.FALSE.
      !CSMEPS=1.0D-6
      !CSMSTEPS=1
      !CSMQUENCHES=1
      !CSMMAXIT=0
      !CHECKMARKOVT=.FALSE.
      !PERCOLATET=.FALSE.
      !PERCCUT=1.0D100
      !GENALT=.FALSE.
!op226> End initializations </init>  }}}

