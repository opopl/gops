
!> @name COMMONS
!
!> @brief Declarations for common variables 
!
MODULE COMMONS

IMPLICIT NONE
SAVE

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
!>                              ...
INTEGER :: NATOMS 
INTEGER :: NSAVE
INTEGER :: ISTEP
INTEGER :: NCORE
INTEGER :: NQ
INTEGER :: M_LBFGS

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COORDS
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HESS
!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GRAD

INTEGER,ALLOCATABLE :: FF(:),INTEFF(:) ! NSAVE
DOUBLE PRECISION, ALLOCATABLE :: QMIN(:), INTEQMIN(:) ! NSAVE
DOUBLE PRECISION, ALLOCATABLE :: QMINP(:,:), INTEQMINP(:,:)

DOUBLE PRECISION TEMP, STEP, OSTEP, ASTEP, ACCRAT, EPREV

INTEGER :: INFIX
! for pulling 
INTEGER :: PATOM1
INTEGER :: PATOM2

PARAMETER (PI=3.141592654D0)

END MODULE COMMONS
