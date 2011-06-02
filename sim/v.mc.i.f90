     
     DOUBLE PRECISION, DIMENSION(:) :: TARGETS
     DOUBLE PRECISION :: TFAC=1.0D0
     DOUBLE PRECISION :: EDIFF=0.02D0
     DOUBLE PRECISION :: ACCRAT=0.5D0
     DOUBLE PRECISION :: TEMP=0.035D0
     DOUBLE PRECISION :: RADIUS=0.0D0
     ! Maximum allowed energy rise during a minimisation
     DOUBLE PRECISION ::     MAXERISE
     ! Maximum allowed energy fall during a minimisation
     DOUBLE PRECISION ::     MAXEFALL
     ! Used in ACCREJ
     DOUBLE PRECISION :: FAC0=1.05D0
     ! "Fixing" option (regarding STEP, TEMP and accept ratio for quenches) 
     CHARACTER(LEN=100) :: FIXOPT='T'
p
     DOUBLE PRECISION :: STEP=0.3D0
     DOUBLE PRECISION :: OSTEP=0.3D0
     DOUBLE PRECISION :: ASTEP=0.3D0
     DOUBLE PRECISION :: EPREV, ARATIO

     INTEGER :: MCSTEPS=10000
          
     INTEGER :: NQ
     INTEGER :: NACCEPT=50

     