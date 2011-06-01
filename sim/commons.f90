
!> @name COMMONS
!
!> @brief Declarations for common variables 
!
MODULE COMMONS

IMPLICIT NONE
SAVE

! Doxygen {{{
!
!> @param RMASS(3)      (dp)    center-of-mass coordinates
!> @param MCSTEPS       (i)     length of a basin-hopping (BH) run   
!> @param NATOMS        (i)     number of particles in the system
!> @param NCORE         (i)     
!> @param NQ            (i)     
!> @param NSAVE         (i)     number of lowest energy geometries to be saved
!> @param NSTEPS        (i)     number of basin-hopping steps
!> @param STEP          (dp)    maximum step size in BH calculations
!> @param ISTEP         (i)     run index during a BH calculation, 1+NDONE ... NSTEPS 
!> @param TFAC          (dp)    specifies the annealing protocol - temperature TEMP is multiplied by TFAC after every MC step
!> @param M_LBFGS       (i)     used in LBFGS ( used there as M )
!>                              ...
!>             is an INTEGER variable that must be set by the user to
!>             the number of corrections used in the BFGS update. It
!>             is not altered by the routine. Values of M less than 3 are
!>             not recommended; large values of M will result in excessive
!>             computing time. 3<= M <=7 is recommended. Restriction: M_LBFGS>0.
!                              ...
! @param NQ             (i)     quench number
! @param ARATIO         (dp)    acceptance ratio
! @param COORDS         dp(N,3) coordinates
!}}}
!
! declarations {{{


! LOGICALS
LOGICAL PULLT, P46, G46, BLNT

! pairdist vars
INTEGER :: NPAIRS
INTEGER, ALLOCATABLE :: PAIRDIST(:,:)
LOGICAL :: PAIRDISTT

DOUBLE PRECISION :: TFAC=1.0D0
DOUBLE PRECISION :: ECONV=0.02D0

INTEGER :: NACCEPT=50
INTEGER :: NATOMS
INTEGER :: MCSTEPS=10000
INTEGER :: NSAVE=10
INTEGER :: ISTEP
INTEGER :: NCORE
INTEGER :: NQ
INTEGER :: M_LBFGS=4
INTEGER :: MAXBFGS

! File handles
INTEGER :: FH=20
INTEGER :: LFH=FH+1
INTEGER :: ENERGY_FH=FH+2
INTEGER :: MARKOV_FH=FH+3
INTEGER :: BEST_FH=FH+4
INTEGER :: PAIRDIST_FH=FH+5

CHARACTER(LEN=100) LFN
PARAMETER(LFH=10,LFN="GMIN.log")

! Index variables
! 1\le IA\le NATOMS
!       JA=3*IA
INTEGER IA, JA
DOUBLE PRECISION ::     RND(3)  
! =========== 
! FILES {{{
! ===========
! --- Filehandle variables 
!       All filehandle variables end with _FH or FH
!
!> @param LFH                 (i)     the output log file
!> @param ENERGY_FH           (i)     energy.dat
!> @param MARKOV_FH           (i)     markov.dat
!> @param BEST_FH             (i)     best.dat
!> @param PAIRDIST_FH         (i)     pairdist.dat
!
!
! --- File names
!
!> @param LFN           (s)     name for the output log file 
!
! --- Values
!
! }}}
! =========== 

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COORDS, COORDSO, VAT, VATO
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HESS
DOUBLE PRECISION ::     RMASS(3)
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GRAD

INTEGER,ALLOCATABLE :: FF(:),INTEFF(:) ! NSAVE
DOUBLE PRECISION, ALLOCATABLE :: QMIN(:), INTEQMIN(:) ! NSAVE
DOUBLE PRECISION, ALLOCATABLE :: QMINP(:,:), INTEQMINP(:,:)

DOUBLE PRECISION TEMP, STEP, OSTEP, ASTEP, ACCRAT, EPREV, ARATIO

INTEGER :: INFIX
! for pulling 
INTEGER :: PATOM1
INTEGER :: PATOM2
! }}}
PARAMETER (PI=3.141592654D0)

SUBROUTINE INIT_PARS
MAXBFGS=0.4D0
END

! optim {{{
MODULE COMMONS
      IMPLICIT NONE
      ! \param RINDEX   run index - integer suffix to odata, i.e., odata.1, odata.2
      INTEGER :: RINDEX 
      ! \param NATOMS   number of atoms
      INTEGER :: NATOMS
      ! \param ISTATUS
      INTEGER :: ISTATUS
      INTEGER :: NUNIQUE, NOPT, IPRNT, INR, IVEC, IVEC2, ISTCRT, NMOL, NINTS, NRBSITES
      DOUBLE PRECISION :: ZSTAR, RHO, PARAM1=0.0D0, PARAM2, PARAM3, PARAM4, PARAM5, PARAM6, PARAM7, EVDISTTHRESH, NEWRES_TEMP
      LOGICAL REDOPATH, REDOPATHXYZ, RATIOS, REDOPATHNEB

      INTEGER, ALLOCATABLE, DIMENSION(:) :: NR, IATNUM      !   MXATMS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ATMASS !   MXATMS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: STPMAX ! 3*MXATMS
      CHARACTER(LEN=5), ALLOCATABLE, DIMENSION(:) :: ZSYM   !   MXATMS
      DOUBLE PRECISION, ALLOCATABLE :: TAGFAC(:)
      DOUBLE PRECISION, ALLOCATABLE :: RBSITE(:,:), RBSTLA(:,:), STCHRG(:), DPMU(:)
      INTEGER, ALLOCATABLE :: TAGNUM(:)
      INTEGER NTAG
      LOGICAL TAGT, CHANGE_TEMP
! }}}
! pathsample {{{
      INTEGER NMINA, NMINB, NMIN, NTS, SAVELENGTH, NPFOLD, PFOLDINT, NRWBINS, NUSEPAIRS, &
     &        NATOMS, MAXLENGTH, CONNECTIONS, HORDER, ERROR, MINSEP, NRWREACTANT, &
     &        NATTEMPT, NNEW, NTOTAL, NEXCLUDE, NPERMGROUP, BHSTEPS, NGTSIZE,  &
     &        MAXCONN, KAPPA, ISEED, NTAG, NDIHE, NCPU, GTINT, NCONNMAX, BESTPATHLENGTH, NGLY, &
     &        STARTMINA, STARTMINB, WHICHMIN, SECONDMIN, WHICHTS, MAXATTEMPT, COSTFUNCTIONPOWER, NPAIRDONE, NMINDONE, NPAIRFRQ, &
     &        BISECTSTEPS, BISECTMAXATTEMPTS, NDIAGEIG, QRELONE, QRELTWO, MAXTSATTEMPTS, &
     &        INTCONSEP, INTREPSEP 
     
      INTEGER, ALLOCATABLE :: LOCATIONA(:), LOCATIONB(:)
      INTEGER, ALLOCATABLE :: NCONN(:)     ! reallocate MAXMIN when used
      INTEGER, ALLOCATABLE :: BESTPATH(:)  
      INTEGER, ALLOCATABLE :: USEPAIRSMIN(:)
!
! I/O unit numbers
!
      INTEGER, PARAMETER :: UMINDATA=11, UTSDATA=12, UMIN=13, UTS=15

      INTEGER :: MAXMIN=26
      INTEGER :: DMINMAX=10
      INTEGER :: MAXTS=10
      INTEGER :: MAXPAIRS=10
      INTEGER :: MAXDONE=10
      INTEGER :: PAIRDISTMAX=1000
      DOUBLE PRECISION, ALLOCATABLE :: FRQS(:)
      DOUBLE PRECISION, POINTER :: MASS(:)
! MAXMIN
      DOUBLE PRECISION, ALLOCATABLE :: EMIN(:), FVIBMIN(:), PFMIN(:), IXMIN(:),  IYMIN(:), IZMIN(:), &
     &                                 GPFOLD(:), MINDISTMIN(:), MINCURVE(:), MINFRQ2(:)
!
! Changed PAIRDIST to linear scaling rather than quadratic with MAXMIN. DJW 14/4/08
!
      DOUBLE PRECISION, ALLOCATABLE :: PAIRDIST(:) ! dimension MAXMIN*PAIRDISTMAX for DIJINITT only runs
!     REAL*4, ALLOCATABLE :: PAIRDIST(:) ! dimension MAXMIN*(MAXMIN-1)/2 for DIJINITT only runs
! MAXTS
      DOUBLE PRECISION, ALLOCATABLE :: ETS(:), FVIBTS(:), KPLUS(:), KMINUS(:), IXTS(:),  IYTS(:), IZTS(:), NEGEIG(:)
! NATOMS
      DOUBLE PRECISION, ALLOCATABLE :: TAGFAC(:)
! NRWBINS
      DOUBLE PRECISION, ALLOCATABLE :: RWPROB(:)

      DOUBLE PRECISION EDIFFTOL, IDIFFTOL, GEOMDIFFTOL, PFMEAN, TOTALE, TEMPERATURE, PFTOTALA, PFTOTALB, PERTMAX, PERTMIN, &
     &                 PERTVALUE, TAGMASS, PABCONV, REGROUPTHRESH, REGROUPRATETHRESH, CONNECTDIST, &
     &                 ORDERPARAM, BOXLX, BOXLY, BOXLZ, DSCALE, PSCALE, TSTHRESH, MAXBARRIER, MAXDOWNBARRIER, REGROUPPETHRESH, &
     &                 PAIRTHRESH, MAXBREAK, PRODTHRESH, PBTHRESH, OMEGA, EINC, RWBINWIDTH, RWEMAX, RWEMIN, &
     &                 GT2RSwitch, GT2Ptol, EUNTRAPTHRESH, PLANCK, REGROUPFREETHRESH, FREETHRESH, &
     &                 BHACCREJ, BHSTEPSIZE, BHCONV, BHTEMP, BHDISTTHRESH, BHK, BHMAXENERGY, BHSFRAC, &
     &                 BISECTMINDIST, BISECTMAXENERGY, NKMCCYCLES, NGTSWITCH, NTFOLD, TOMEGA, TFOLDTHRESH, DIAGSCALE, &
     &                 CVTMIN, CVTMAX, CVTINC, DOSEMIN, DOSEMAX, DOSEINC, EVCUT, GAMMAFRICTION, &
     &                 INTEPSILON, INTCONSTRAINTDEL, INTCONSTRAINTREP, INTCONSTRAINREPCUT, INTLJDEL, INTLJEPS, NGTCRSWITCH

! AMH
      DOUBLE PRECISION QCONTCUT, RELCOCUT

      DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
      DOUBLE PRECISION TTSSEARCH, TPFOLD, TTFOLD, TGT, TDIJKSTRA, TCONNECTDIST, TKSHORTESTPATHS ! timers

! MAXMIN
      INTEGER, ALLOCATABLE :: HORDERMIN(:), TOPPOINTER(:), MINGROUP(:)
! MAXTS
      INTEGER, ALLOCATABLE :: HORDERTS(:), PLUS(:), MINUS(:), POINTERM(:), POINTERP(:), TSATTEMPT(:)
      INTEGER, ALLOCATABLE :: DMIN1(:), DMIN2(:) ! dimension MINA*MINB for DIJPAIRT runs
! MAXPAIRS
      INTEGER, ALLOCATABLE :: PAIR1(:), PAIR2(:)
! MAXDONE
      INTEGER, ALLOCATABLE :: MINDONE(:)
! NATOMS
      INTEGER, ALLOCATABLE :: NPERMSIZE(:), PERMGROUP(:)
      INTEGER, ALLOCATABLE :: NSETS(:)
      INTEGER, ALLOCATABLE :: SETS(:,:)
      INTEGER, ALLOCATABLE :: TAGNUM(:)
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: FROZEN  

      INTEGER NFSTART, NFFINISH, NINTS, NCONNMIN, CONNECTMIN1, CONNECTMIN2, NFREEZE, NPATHS, PARALLEL

      LOGICAL YESNO, TEST1, TEST2, DEBUG, PRINTT, ADDPT, TWOD, BULKT, ANGLEAXIS, TAGT, CHARMMT, AMBERT, STARTFROMPATH, &
     &        KMCT, UNRST, KMCCOMMITT, REGROUPT, REGROUPRATET, REGROUPPET, NOPOINTS, ADDPATH, NGTT, GTT, GT2T, TRIPLES, &
     &        STARTTRIPLES, ADDTRIPLES, DIJKSTRAT, DIJPAIRT, DIJINITT, EXTRACTMINT, EXTRACTTST, DIJKSTRAWAITT, &
     &        EXPCOSTFUNCTION, COPYOPTIMT, CALCORDERT, CONNECTREGIONT, SHORTCUTT, MERGEDBT, UNTRAPT, AMHT,  AMHALLATOMMINT, &
     &        CHECKCONNECTIONST, AMHALLATOMTST, AMHQT, AMHQCONTT ,AMHRMSDT, AMHRELQT, AMH_RELCOT, DIAGT, ARNOLDIT, &
     &        GT2Sparse, GT2Switch, GT2AltPbb, GT2Rescale, GT2Normalise, GT2DisconnectSources, BARRIERSORT, &
     &        PERMDIST, PERMISOMER, RIGIDBODY, DIJINITSTARTT, DIJINITCONTT, RETAINSP, REMOVESP, NOFRQS, &
     &        BARRIERSHORT, FREEZE, RATESHORT, DUMMYRUNT, REWEIGHTT, REGROUPFREET, REGROUPFREEABT, READMINT, &
     &        DUMPGROUPST, FREEPAIRT, KSHORTESTPATHST, KSHORT_FULL_PRINTT, DIJINITFLYT, BHINTERPT, ICINTERPT, &
     &        DUMMYTST, DOCKT, DSTAGE(6), USEPAIRST, LOWESTFRQT, BISECTT, NGTDISCONNECTALL, ANGLEAXIS2, TFOLDT, &
     &        SLURMT, INDEXCOSTFUNCTION, CVT, DOST, IMFRQT, CLOSEFILEST, PULLT, FRICTIONT, &
     &        INTCONSTRAINTT, CHECKCONINT, INTLJT, INTERPCOSTFUNCTION, REMOVEUNCONNECTEDT, &
     &        DBPT, DBPTDT, EFIELDT, MSSTOCKT, PAHAT, PATCHYDT, STOCKAAT, RBAAT, RBSYMT 
      LOGICAL, ALLOCATABLE :: SHIFTABLE(:)
      CHARACTER(LEN=80) COORDSLIGANDSTR, COORDSCOMPLEXSTR, COORDSPROTEINSTR
      CHARACTER(LEN=80) EXEC,EXECGMIN
      CHARACTER(LEN=80) PATHNAME, MINNAME
      CHARACTER(LEN=80) COPYFILES
      CHARACTER(LEN=80) USEPAIRSFILE
      CHARACTER(LEN=2) DIRECTION
      CHARACTER(LEN=2) UNCONNECTEDS
      CHARACTER(LEN=5) ZSYM
      CHARACTER(LEN=2), ALLOCATABLE ::  ZSYMBOL(:)
      CHARACTER(LEN=1) ENSEMBLE

      CHARACTER(LEN=4), ALLOCATABLE :: RESLABEL(:), ATOMLABEL(:)
      INTEGER, ALLOCATABLE :: RESNUMBER(:)
      INTEGER, ALLOCATABLE :: NCONNGROUP(:) ! jmc

!     DC430 >
      INTEGER :: NRBSITES, NTSITES, NRBGROUP, PAHID
      DOUBLE PRECISION, ALLOCATABLE :: RBSITE(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: RBOPS(:,:)

!
! Constraint potential stuff.
!
      INTEGER, ALLOCATABLE :: CONI(:), CONJ(:), CONION(:), CONJON(:)
      INTEGER, ALLOCATABLE :: REPI(:), REPJ(:)
      DOUBLE PRECISION, ALLOCATABLE :: REPCUT(:),  NREPCUT(:), CONDISTREF(:), CONDISTREFLOCAL(:)
      DOUBLE PRECISION INTCONSTRAINTTOL, REPCON
      INTEGER, ALLOCATABLE :: NREPI(:), NREPJ(:)
      INTEGER :: NNREPULSIVE, NCONSTRAINT, NREPULSIVE, MAXCONUSE
      INTEGER :: NREPMAX=10
      INTEGER :: INTIMAGE=0
      INTEGER :: INTCONMAX=10
      INTEGER :: NTRYING=0
      LOGICAL, ALLOCATABLE :: CONACTIVE(:)
      LOGICAL, ALLOCATABLE :: ATOMACTIVE(:)
!
! SIS epidemiological model stuff
!
      INTEGER :: SMAX, IMAX, POPSS
      DOUBLE PRECISION :: SISMU, SISKAPPA, SISBETA
      LOGICAL :: SIST
! }}}

END MODULE COMMONS
