
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
!> @param TEMP          (dp)    temperature
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

! LOGICALS {{{
LOGICAL PULLT=.TRUE.
LOGICAL P46=.FALSE.
LOGICAL G46=.FALSE.
LOGICAL BLNT=.FALSE.
LOGICAL TARGET=.FALSE.
LOGICAL TRACKDATAT=.FALSE.
LOGICAL DEBUG=.FALSE.
! whether we are doing a final quench
LOGICAL :: FQFLAG=.FALSE.
! }}}

! reading the data file/command-line

INTEGER, PARAMETER :: SLEN=100
INTEGER, PARAMETER :: MAXNARGS=20

INTEGER NARGS
CHARACTER(LEN=SLEN),DIMENSION(MAXNARGS) :: ARGS

! pairdist vars
INTEGER :: NPAIRS
INTEGER, ALLOCATABLE :: PAIRDIST(:,:)
LOGICAL :: PAIRDISTT

DOUBLE PRECISION, DIMENSION(:) :: TARGETS
DOUBLE PRECISION :: TFAC=1.0D0
DOUBLE PRECISION :: ECONV=0.02D0
DOUBLE PRECISION :: ACCRAT=0.5D0
DOUBLE PRECISION :: TEMP=0.035D0

! initial time
DOUBLE PRECISION :: TSTART

! sloppy quenches
DOUBLE PRECISION :: SQMAX
! final quenches
DOUBLE PRECISION :: FQMAX

! integers {{{
INTEGER :: NACCEPT=50
INTEGER :: NATOMS=46
INTEGER :: MCSTEPS=10000
INTEGER :: NSAVE=10
INTEGER :: ISTEP
INTEGER :: NCORE
INTEGER :: NQ
INTEGER :: M_LBFGS=4
INTEGER :: MAXBFGS=0.4D0
! }}}

! File handling {{{
INTEGER :: FH=20
INTEGER :: LFH=FH+1
INTEGER :: ENERGY_FH=FH+2
INTEGER :: MARKOV_FH=FH+3
INTEGER :: BEST_FH=FH+4
INTEGER :: PAIRDIST_FH=FH+5
INTEGER :: COORDS_FH=FH+6

CHARACTER(LEN=100) LFN
PARAMETER(LFH=10,LFN="GMIN.log")
! }}}

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

DOUBLE PRECISION :: STEP=0.3D0
DOUBLE PRECISION :: OSTEP=0.3D0
DOUBLE PRECISION :: ASTEP=0.3D0

EPREV, ARATIO

INTEGER :: INFIX
! for pulling 

INTEGER :: PATOM1
INTEGER :: PATOM2
DOUBLE PRECISION :: PFORCE

PARAMETER (PI=3.141592654D0)
! }}}

END MODULE COMMONS
