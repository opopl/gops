
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
