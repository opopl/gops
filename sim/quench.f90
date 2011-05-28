
      SUBROUTINE QUENCH(FINALQUENCH,ITER,TIME,BRUN,QDONE,P)

      USE COMMONS
      USE PORFUNCS

      IMPLICIT NONE

      ! SUBROUTINE PARAMETERS {{{

      LOGICAL,INTENT(IN) :: FINALQUENCH
      INTEGER ITER
      DOUBLE PRECISION TIME
      INTEGER BRUN,QDONE
      DOUBLE PRECISION P(3*NATOMS)

      ! }}}

      ! LOCAL PARAMETERS  {{{

      DOUBLE PRECISION SSAVE
      INTEGER NOPT

      ! CFLAG - test convergence
      LOGICAL CFLAG

      ! }}}
      
      DOUBLE PRECISION POTEL
      
      COMMON /MYPOT/ POTEL

!  Turn on guiding potentials. These get turned off in potential.F when the RMS force is small enough.

      SSAVE=STEP
      NFIX=0
      NOPT=3*NATOMS

!  FINALQUENCH is set for the final quenches with tighter convergence criteria.

      IF (FINALQUENCH) THEN
         GMAX=CQMAX
      ELSE
         GMAX=BQMAX
      ENDIF

      QDONE=0
      P=COORDS

      COMPON=.FALSE.

      CALL MYLBFGS(3*NATOMS,M_LBFGS,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.)
      POTEL=EREAL

      IF (CFLAG) QDONE=1

      IF (.NOT.CFLAG) THEN
            WRITE(MYUNIT,'(A,I6,A)') 'WARNING - Final Quench ',NQ,'  did not converge'
      ENDIF

      CALL MYCPU_TIME(TIME)

      RES=.FALSE.

      IF (.NOT.NORESET) THEN
         DO J1=1,3*(NATOMS-NSEED)
            COORDS(J1)=P(J1)
         ENDDO
         DO J1=1,NATOMS
            VAT(J1)=VT(J1)
         ENDDO
      ENDIF

      RETURN

      END
