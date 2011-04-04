
!> @name COMMONS
!
!> @brief Declarations for common variables 
!
MODULE COMMONS

IMPLICIT NONE
SAVE

!> @param MCSTEPS       (i)     length of a BH run   
!> @param NATOMS        (i)     number of particles in the system
!> @param NCORE         (i)     
!> @param NQ            (i)     
!> @param NSAVE         (i)     number of lowest energy geometries to be saved
!> @param NSTEPS        (i)     
!> @param STEP          (dp)    maximum step size in BH calculations
!> @param ISTEP         (i)     run index during a BH calculation, 1+NDONE ... NSTEPS 
!> @param TFAC          (dp)    specifies the annealing protocol - temperature TEMP is multiplied by TFAC after every MC step
INTEGER :: NATOMS 
INTEGER :: NSAVE
INTEGER :: ISTEP
INTEGER :: NCORE
INTEGER :: NQ

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COORDS

INTEGER,ALLOCATABLE :: FF(:),INTEFF(:) ! NSAVE
DOUBLE PRECISION, ALLOCATABLE :: QMIN(:), INTEQMIN(:) ! NSAVE
DOUBLE PRECISION, ALLOCATABLE :: QMINP(:,:), INTEQMINP(:,:)

DOUBLE PRECISION TEMP, STEP, OSTEP, ASTEP, ACCRAT, EPREV

INTEGER :: INFIX

PARAMETER (PI=3.141592654D0)

END MODULE COMMONS
