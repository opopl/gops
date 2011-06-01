      ! local {{{

      COMMON /MYPOT/ E
    
      ! }}}
      ! subroutine {{{

      INTEGER,INTENT(IN) :: NSTEPS
      DOUBLE PRECISION,INTENT(IN) :: SCALEFAC 
      DOUBLE PRECISION, INTENT(INOUT) :: SCREENC(NATOMS,3)

      ! local

      INTEGER NFAIL, NFAILT, NSUCCESS, NSUCCESST
      INTEGER ITERATIONS,JACCPREV,QDONE,NDONE
      DOUBLE PRECISION, DIMENSION(NATOMS,3) ::  GRAD
      DOUBLE PRECISION, DIMENSION(NATOMS) :: RVAT, RVATO
!  EPREV saves the previous energy in the Markov chain.
!  EBEST and JBEST record the lowest energy since the last reseeding and the
!  step it was attained at. BESTCOORDS contains the corresponding coordinates.
      DOUBLE PRECISION :: EBEST,EPREV,TIME,E
      INTEGER JBEST

      NSUCCESS=0
      NFAIL=0
      NSUCCESST=0
      NFAILT=0
      ! }}}


