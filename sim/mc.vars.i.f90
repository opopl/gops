
      ! local {{{
      LOGICAL CONVG

      ! ATEST: .TRUE. if the move is accepted; .FALSE. otherwise
      LOGICAL :: ATEST
      INTEGER NFAIL, NFAILT, NSUCCESS, NSUCCESST
      INTEGER ITERATIONS,NDONE
      !
      ! JACCPREV: the index of the last step when the move 
      !           was accepted
      INTEGER JACCPREV
      DOUBLE PRECISION ::       TIME
      DOUBLE PRECISION ::       ARATIO
      ! QE       - energy after each quench
      ! QEPREV   - energy from the previous quench
      ! DQE=QE-QEPREV
      DOUBLE PRECISION :: QE, QEPPREV, DQE
      DOUBLE PRECISION :: MCTEMP

      DOUBLE PRECISION :: EBEST
      INTEGER :: JBEST

      COMMON /MYPOT/ QE

      ! }}}


