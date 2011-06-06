     
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

     
