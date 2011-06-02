
     ! DGUESS: Guess for initial diagonal elements in LBFGS
     DOUBLE PRECISION :: DGUESS=0.0D0
     
     ! Number of LBFGS updates
     INTEGER :: M_LBFGS=4

     ! Maximum BFGS step size
     DOUBLE PRECISION :: MAXBFGS=0.4D0

     ! maximal number of iterations (for sloppy and final quenches)
     INTEGER MAXIT=500
     
     ! sloppy quenches
     DOUBLE PRECISION :: SQMAX
     ! final quenches
     DOUBLE PRECISION :: FQMAX
     
