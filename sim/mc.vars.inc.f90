      ! local {{{

      COMMON /MYPOT/ E
    
!  EPREV saves the previous energy in the Markov chain.
!  EBEST and JBEST record the lowest energy since the last reseeding and the
!  step it was attained at. BESTCOORDS contains the corresponding coordinates.
      ! }}}
      ! subroutine {{{

      INTEGER NSTEPS
      DOUBLE PRECISION :: SCALEFAC, SCREENC(NATOMS,3)

      ! local

      INTEGER NFAIL, NFAILT, NSUCCESS, NSUCCESST, NQTOT
      INTEGER ITERATIONS,JACCPREV,QDONE,NDONE
      DOUBLE PRECISION, DIMENSION(NATOMS,3) ::  GRAD,SAVECOORDS,BESTCOORDS,TEMPCOORDS,RCOORDS,RCOORDSO
      DOUBLE PRECISION, DIMENSION(NATOMS) :: RVAT, RVATO
      DOUBLE PRECISION :: EBEST,EPREV,TIME,E

      NSUCCESS=0
      NFAIL=0
      NSUCCESST=0
      NFAILT=0
      ! }}}


