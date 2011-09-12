
MODULE MC

USE V
USE FUNC
USE PORFUNCS

IMPLICIT NONE

CONTAINS

SUBROUTINE QUENCH(R,ITERS,TIME,CONVG)
! declarations {{{

      IMPLICIT NONE

      ! SUBROUTINE PARAMETERS {{{

      ! number of iterations spent in the computation
      INTEGER,INTENT(OUT) :: ITERS
      ! time spent in computation 
      DOUBLE PRECISION, INTENT(OUT) :: TIME
      ! exit code:
      !         QDONE=0 if convergence achieved
      !         QDONE=1 otherwise
      LOGICAL,INTENT(OUT) :: CONVG
      ! input/output coordinate array 
      DOUBLE PRECISION,INTENT(INOUT) :: R(NATOMS,3)
      ! coordinates used in the minimizations
      DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: X

      ! }}}
      ! LOCAL PARAMETERS  {{{

      DOUBLE PRECISION SSAVE
      INTEGER NOPT

      ! E - energy 
      DOUBLE PRECISION QE
      COMMON /MYPOT/ QE

      ! }}}
     ! }}} 
      ! body {{{
!  FQFLAG is set for the final quenches with tighter convergence criteria.

      if (.not.allocated(X)) ALLOCATE(X(3*NATOMS))
      !write(*,*) 'in QUENCH: allocate(X)'

      X=PACK(R,.true.)          ! R(NATOMS,3) => X(3*NATOMS)

      IF (FQFLAG) THEN
         GMAX=FQMAX
      ELSE
         GMAX=SQMAX
      ENDIF

      CONVG=.FALSE.

      ! Invoking LBFGS:
      !         X                       input coordinates
      !         DIAGCO=.FALSE           don't provide Hk0 at each iteration
      !         GMAX                    = EPS in LBFGS
      !         CONVG                   convergence, MFLAG in LBFGS
      !         EREAL                   returns energy
      !         MAXIT                   = ITMAX in LBFGS
      !         ITERS                    = ITDONE in LBFGS
      !                                 (number of iterations needed to obtain convergence)
      !         RESET=.TRUE.            Reset ITER=0 in LBFGS
      !
      write(*,*) 'invoking LBFGS'
      CALL MYLBFGS(X,.FALSE.,GMAX,CONVG,QE,MAXIT,ITERS,.TRUE.)
      write(*,*) 'MYLBFGS done'

      IF (.NOT.CONVG) THEN
            WRITE(LFH,'(A,I6,A)') 'WARNING - Final Quench ',NQ,'  did not converge'
      ENDIF

      CALL CPU_TIME(TIME)
      IF(ALLOCATED(X)) DEALLOCATE(X)

      RETURN
      ! }}}
END SUBROUTINE

SUBROUTINE MCRUN(NSTEPS,SCALEFAC,SCREENC)
! {{{

IMPLICIT NONE

include "mc.vars.i.f90" ! Variable declarations 
include "fmt.i.f90"     ! Formats
include "mc.init.i.f90" ! Initializations 
include "mc.cien.i.f90" ! Calculate the initial energy and save in EPREV
include "mc.bh.i.f90"   ! Main basin-hopping loop
37 CONTINUE
include "mc.w.111.i.f90" ! write to LFH: Qu, E, Steps, RMS, Markov E, t (elapsed time)

WRITE(LFH,21)   ' Acceptance ratio for run=',  ARATIO,	&
                ' Step=',                      STEP, 	&
                ' Angular step factor=',       ASTEP,	&
                ' Temperature T=',             TEMP	
RETURN
! }}}
END SUBROUTINE

SUBROUTINE ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)
! {{{

IMPLICIT NONE
include "fmt.i.f90"     ! Formats

INTEGER,INTENT(INOUT) :: NSUCCESS, NFAIL, NFAILT, NSUCCESST
DOUBLE PRECISION :: ARATIO,FAC

ARATIO=ONE*NSUCCESS/(ONE*(NSUCCESS+NFAIL))

IF (ARATIO.GT.ACCRAT) THEN
	FAC=FAC0
ELSE
	FAC=ONE/FAC0
ENDIF
        SELECTCASE(FIXOPT)
               CASE('T') ! fix temperature, but adjust step(s)
	                STEP=FAC*STEP
	                ASTEP=ASTEP*FAC
               CASE('S') ! fix step(s), but adjust temperature
                        TEMP=TEMP/FAC
               CASE('TS') ! fix both steps and temperature
        ENDSELECT

WRITE(LFH,23) 'Acceptance ratio for previous ',        NACCEPT,	&
                ' steps=',                             ARATIO,	&
                '  FAC=',                                FAC

WRITE(LFH,'(A)',ADVANCE='NO') 'Steps are now:'
WRITE(LFH,'(A,F10.4)',ADVANCE='NO') '  STEP=',STEP    
IF(ASTEP.GT.0.D0) WRITE(LFH,'(A,F10.4)',ADVANCE='NO')   '  ASTEP=',ASTEP 
WRITE(LFH,'(A,F10.4)') ' Temperature is now:',TEMP
!
NSUCCESST=NSUCCESST+NSUCCESS
NFAILT=NFAILT+NFAIL
NSUCCESS=0
NFAIL=0 
!
RETURN
! }}}
END SUBROUTINE

SUBROUTINE TRAN(ENEW,EOLD,ATEST,MCTEMP)
! {{{
      IMPLICIT NONE

      DOUBLE PRECISION,INTENT(IN) :: ENEW, EOLD, MCTEMP
      LOGICAL,INTENT(OUT) :: ATEST
      DOUBLE PRECISION :: RN

      DOUBLE PRECISION :: DE,T 

      DE=ENEW-EOLD
      T=MAX(MCTEMP,ZERO)
      ATEST=.TRUE.

         IF (DE.GE.ZERO) THEN
            RN=DPRAND()
            IF (DEXP(-DE/T).LE.RN) THEN
               ATEST=.FALSE.
            ENDIF
         ENDIF

      RETURN 
      ! }}}
END  SUBROUTINE

END MODULE

