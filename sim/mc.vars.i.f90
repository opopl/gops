      ! local {{{

      COMMON /MYPOT/ E
    
      ! }}}
      ! subroutine {{{

      INTEGER, INTENT(IN) :: NSTEPS
      DOUBLE PRECISION, INTENT(IN) :: SCALEFAC 
      DOUBLE PRECISION, INTENT(INOUT) :: SCREENC(NATOMS,3)

      ! local

      INTEGER NFAIL, NFAILT, NSUCCESS, NSUCCESST
      INTEGER ITERATIONS,QDONE,NDONE
      DOUBLE PRECISION, DIMENSION(NATOMS,3) ::  GRAD
      DOUBLE PRECISION ::       TIME,E,EPREV

      NSUCCESS=0
      NFAIL=0
      NSUCCESST=0
      NFAILT=0
      ! }}}

