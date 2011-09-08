
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: COORDS, COORDSO
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HESS
     DOUBLE PRECISION ::     RMASS(3)
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GRAD

     ! initial time

     DOUBLE PRECISION :: TSTART, TFINISH
     INTEGER :: NATOMS

     ! potential type

     CHARACTER(LEN=SLEN) PTYPE

     INTEGER :: NSEED

     DOUBLE PRECISION ::       RMS
