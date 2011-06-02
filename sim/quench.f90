
      SUBROUTINE QUENCH(P,ITER,TIME,QDONE)
! declarations {{{
      USE V
      USE PORFUNCS

      IMPLICIT NONE

      ! SUBROUTINE PARAMETERS {{{

      ! number of iterations spent in the computation
      INTEGER,INTENT(OUT) :: ITER
      ! time spent in computation 
      DOUBLE PRECISION, INTENT(OUT) :: TIME
      ! exit code:
      !         QDONE=0 if convergence achieved
      !         QDONE=1 otherwise
      INTEGER,INTENT(OUT) :: QDONE
      ! coordinate array used in minimizations 
      DOUBLE PRECISION,INTENT(INOUT) :: P(3*NATOMS)

      ! }}}
      ! LOCAL PARAMETERS  {{{

      DOUBLE PRECISION SSAVE
      INTEGER NOPT

      ! CFLAG - test convergence
      LOGICAL CFLAG
      ! E - energy 
      DOUBLE PRECISION E
      COMMON /MYPOT/ E

      ! }}}
     ! }}} 
      ! body {{{
!  FQFLAG is set for the final quenches with tighter convergence criteria.

      IF (FQFLAG) THEN
         GMAX=FQMAX
      ELSE
         GMAX=SQMAX
      ENDIF

      QDONE=1

      ! Invoking LBFGS:
      !         P                       input coordinates
      !         DIAGCO=.FALSE           don't provide Hk0 at each iteration
      !         GMAX                    = EPS in LBFGS
      !         CFLAG                   convergence
      !         EREAL                   returns energy
      !         MAXIT                   = ITMAX in LBFGS
      !         ITER                    = ITDONE in LBFGS
      !                                 (number of iterations needed to obtain convergence)
      !         RESET=.TRUE.            Reset ITER=0 in LBFGS
      !
      CALL MYLBFGS(P,.FALSE.,GMAX,CFLAG,E,MAXIT,ITER,.TRUE.)

      IF (CFLAG) QDONE=0

      IF (.NOT.CFLAG) THEN
            WRITE(MYUNIT,'(A,I6,A)') 'WARNING - Final Quench ',NQ,'  did not converge'
      ENDIF

      CALL MYCPU_TIME(TIME)

      RETURN
      ! }}}

      END
