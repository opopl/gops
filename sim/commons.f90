
!> @name COMMONS
!
!> @brief Declarations for common variables 
!
MODULE GMIN_COMMONS

IMPLICIT NONE
SAVE

! Doxygen {{{
!> @param RMASS(3)      (dp)    center-of-mass coordinates
!> @param MCSTEPS       (i)     length of a basin-hopping (BH) run   
!> @param NATOMS        (i)     number of particles in the system
!> @param NCORE         (i)     
!> @param NQ            (i)     
!> @param NSAVE         (i)     number of lowest energy geometries to be saved
!> @param NSTEPS        (i)     
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
INTEGER :: NATOMS
INTEGER :: NSAVE
INTEGER :: ISTEP
INTEGER :: NCORE
INTEGER :: NQ
INTEGER :: M_LBFGS
INTEGER :: MAXBFGS
INTEGER :: LFH
INTEGER :: ENERGY_FH, MARKOV_FH, BEST_FH, PAIRDIST_FH
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

END MODULE COMMONS
