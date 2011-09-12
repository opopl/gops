
!> @name V
!> @brief Variable/Constants declarations 

MODULE V

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
! CONSTANTS {{{
DOUBLE PRECISION :: PI,ZERO,ONE
INTEGER ::      SLEN
! reading the data file/command-line

INTEGER :: MAXNARGS

PARAMETER (PI=3.141592654D0)
PARAMETER (ZERO=1.0D-100)
PARAMETER (ONE=1.0D0)
PARAMETER (SLEN=100,MAXNARGS=100)
! }}}
     ! LOGICALS {{{
    
     LOGICAL :: PULLT
     LOGICAL :: P46
     LOGICAL :: G46
     LOGICAL :: BLNT
     LOGICAL :: TARGET
     LOGICAL :: TRACKDATAT
     LOGICAL :: DEBUG
     LOGICAL :: LBFGST
     LOGICAL :: BFGST
     LOGICAL :: RMST
     ! use the 'data' file?
     LOGICAL :: USEKW
     ! whether we are doing a final quench
     LOGICAL :: FQFLAG

     ! }}}
! LBFGS {{{
     ! DGUESS: Guess for initial diagonal elements in LBFGS
     DOUBLE PRECISION :: DGUESS
     
     ! Number of LBFGS updates
     INTEGER :: M_LBFGS

     ! Maximum BFGS step size
     DOUBLE PRECISION :: MAXBFGS

     ! maximal number of iterations (for sloppy and final quenches)
     INTEGER :: MAXIT
     
     ! sloppy quenches
     DOUBLE PRECISION :: SQMAX
     ! final quenches
     DOUBLE PRECISION :: FQMAX
     ! GMAX 
     DOUBLE PRECISION :: GMAX
     ! }}}
! FILES {{{
! comments {{{
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
! }}}
! {{{
INTEGER, PARAMETER :: FH=20
INTEGER :: LFH
INTEGER :: ENERGY_FH
INTEGER :: MARKOV_FH
INTEGER :: BEST_FH
INTEGER :: PAIRDIST_FH
INTEGER :: COORDS_FH
INTEGER :: DATA_FH
INTEGER :: RMSD_FH

CHARACTER(LEN=100) LFN
PARAMETER(LFN="gmi.log")
! }}}
! }}}
! MC {{{

     DOUBLE PRECISION :: TFAC
     DOUBLE PRECISION :: EDIFF
     DOUBLE PRECISION :: ACCRAT
     DOUBLE PRECISION :: TEMP
     DOUBLE PRECISION :: RADIUS
     ! Maximum allowed energy rise during a minimisation
     DOUBLE PRECISION ::     MAXERISE
     ! Maximum allowed energy fall during a minimisation
     DOUBLE PRECISION ::     MAXEFALL
     ! Used in ACCREJ
     DOUBLE PRECISION :: FAC0
     ! "Fixing" option (regarding STEP, TEMP and accept ratio for quenches) 
     CHARACTER(LEN=100) :: FIXOPT

     DOUBLE PRECISION :: STEP
     DOUBLE PRECISION :: OSTEP
     DOUBLE PRECISION :: ASTEP
     DOUBLE PRECISION :: QEPREV
     DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:) :: TARGETS
     INTEGER NTARGETS

     INTEGER :: MCSTEPS
          
     INTEGER :: NQ
     INTEGER :: NACCEPT
     INTEGER :: NRELAX

     
! }}}
! pull {{{

     INTEGER :: PATOM1
     INTEGER :: PATOM2
     DOUBLE PRECISION :: PFORCE

! }}}
! general {{{

! main coordinates variable
DOUBLE PRECISION, ALLOCATABLE :: SCREENC(:)

! reading the data file/command-line

INTEGER NARGS
CHARACTER(LEN=SLEN), DIMENSION(MAXNARGS) :: ARGS

! number of lowest energy configurations to save
INTEGER :: NSAVE
INTEGER :: ISTEP
INTEGER :: NCORE

     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COORDS, COORDSO
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HESS
     DOUBLE PRECISION ::     RMASS(3)
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: GRAD

     ! initial time

     DOUBLE PRECISION :: TSTART, TFINISH
     INTEGER :: NATOMS

     ! potential type

     CHARACTER(LEN=SLEN) PTYPE

     ! name for the model system
     CHARACTER(LEN=100) MODEL

     INTEGER :: NSEED

     DOUBLE PRECISION ::       RMS
! }}}
! pairdist {{{

     ! pairdist vars
     
     INTEGER :: NPAIRS
     INTEGER, ALLOCATABLE :: PAIRDIST(:,:)
     LOGICAL :: PAIRDISTT

! }}}

END MODULE V

