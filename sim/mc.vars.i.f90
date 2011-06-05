      ! subroutine {{{

      INTEGER, INTENT(IN) :: NSTEPS
      DOUBLE PRECISION, INTENT(IN) :: SCALEFAC 
      DOUBLE PRECISION, INTENT(INOUT) :: SCREENC(NATOMS,3)

      ! local

      ! ATEST: .TRUE. if the move is accepted; .FALSE. otherwise
      LOGICAL :: ATEST
      INTEGER NFAIL, NFAILT, NSUCCESS, NSUCCESST
      INTEGER ITERATIONS,QDONE,NDONE
      !
      ! JACCPREV: the index of the last step when the move 
      !           was accepted
      INTEGER JACCPREV
      DOUBLE PRECISION ::       TIME
      ! QE       - energy after each quench
      ! QEPREV   - energy from the previous quench
      ! DQE=QE-QEPREV
      DOUBLE PRECISION :: QE, QEPPREV, DQE
      DOUBLE PRECISION :: MCTEMP

      DOUBLE PRECISION :: EBEST
      INTEGER :: JBEST

      COMMON /MYPOT/ QE

      ! }}}


