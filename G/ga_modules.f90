!mo361
!Modules for genetic algorithms
!
! Parameters for genetic algorithm searches
!
MODULE MYGA_PARAMS
   IMPLICIT NONE
   INTEGER MYGA_NSTRUC, MYGA_NOFF! Population size
   INTEGER :: MYGA_NMUT=0 !Number of mutants. Re-sets every generation
   INTEGER MYGA_GENS ! Number of generations
   INTEGER :: MYGA_TOURN_SIZE = 3
   INTEGER :: MYGA_EPOCH_SAVE = 0
   DOUBLE PRECISION :: MYGA_MUT_RATE=0.1 ! Mutation rate
   DOUBLE PRECISION MYGA_BQMAX,MYGA_CQMAX !Loose & tight convergence thresholds
   DOUBLE PRECISION :: MYGA_DUPLICATE_ETHRESH=0.00001 !Energy threshold for duplicate predator
!  DOUBLE PRECISION :: MYGA_DUPLICATE_GTHRESH=10*4.d0*DATAN(1.D0)/180d0 !Geometry threshold for duplicate predator
   DOUBLE PRECISION :: MYGA_DUPLICATE_GTHRESH=0.174532925D0
   DOUBLE PRECISION :: MYGA_EPOCH_THRESH=0.01
   DOUBLE PRECISION :: MYGA_LAST_ENERGY=1D10
   LOGICAL :: MYGA_L_ROUL=.FALSE.
   LOGICAL :: MYGA_L_EPOCH=.FALSE.
   LOGICAL :: MYGA_L_CHAIN=.FALSE.!Generate random structures as a continuous chain
   LOGICAL :: MYGA_L_SPHERE=.FALSE.!
   DOUBLE PRECISION :: MYGA_SPHERE_RADIUS=3.0D0
   INTEGER :: MYGA_COUNT_EPOCH = 1
   INTEGER :: MYGA_COUNT_MIN = 0
   INTEGER :: CURR_GEN=0
END MODULE MYGA_PARAMS
!
! Module to hold whole population of genetic algorithm
!
MODULE MYGA_POPULATION
   USE MYGA_PARAMS
   IMPLICIT NONE
   INTEGER, ALLOCATABLE :: MYGA_POP_FOUND(:)
   DOUBLE PRECISION, ALLOCATABLE :: MYGA_POP_ENERGY(:)
   DOUBLE PRECISION, ALLOCATABLE :: MYGA_POP_GENOME(:,:) !3*natoms,popsize
   DOUBLE PRECISION, ALLOCATABLE :: MYGA_POP_COORDS(:,:) !3*natoms,popsize
   DOUBLE PRECISION, ALLOCATABLE :: MYGA_POP_FITNESS(:) !Fitness for roulette selection
END MODULE MYGA_POPULATION
