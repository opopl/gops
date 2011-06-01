
      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)
! {{{
      USE COMMONS
      USE PORFUNCS

      IMPLICIT NONE

include mc.vars.inc.f90
include fmt.inc.f90
include mc.cien.inc.f90 ! Calculate the initial energy and save in EPREV
include mc.bh.inc.f90   ! Main basin-hopping loop

37    CONTINUE

ARATIO=NSUCCESST*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST+NFAILT))

WRITE(LFH,111)  ' Qu ',         NQ,	&
                ' E=',          POTEL,	&
                ' steps=',      ITERATIONS,	&
                ' RMS=',        RMS,	&
                ' Markov E=',   POTEL,	&
                ' t=',          TIME-TSTART	

WRITE(LFH,21)   ' Acceptance ratio for run=',  ARATIO,	&
                ' Step=',                      STEP,     	&
                ' Angular step factor=',       ASTEP,	&
                ' T=',                         TEMP	&
RETURN
! }}}
END

SUBROUTINE ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)
! {{{
USE COMMONS

IMPLICIT NONE

INTEGER NSUCCESS, NFAIL, NFAILT, NSUCCESST, J1, J2 
DOUBLE PRECISION HWMAX,P0,FAC
COMMON /EV/ EVAP, EVAPREJECT

P0=1.D0*NSUCCESS/(1.D0*(NSUCCESS+NFAIL))

IF (P0.GT.ACCRAT) THEN
	FAC=1.05D0
ELSE
	FAC=1.D0/1.05D0
ENDIF
	STEP=FAC*STEP
	ASTEP=ASTEP*FAC
!
! Prevent steps from growing out of bounds. The value of 1000 seems sensible, until
! we do something with such huge dimensions?!

STEP=MIN(STEP,1.0D3)
OSTEP=MIN(OSTEP,1.0D3)
ASTEP=MIN(ASTEP,1.0D3)
!
WRITE(LFH,23) 'Acceptance ratio for previous ',        NACCEPT,	&
                ' steps=',                               P0,	&
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
END

      SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,RANDOM,MCTEMP)
! {{{
      IMPLICIT NONE

      DOUBLE PRECISION ENEW, EOLD, DPRAND, RANDOM, MCTEMP

      DOUBLE PRECISION ::       DE,T 
      LOGICAL ATEST

      DE=ENEW-EOLD
      T=MAX(MCTEMP,1.0D-100)

         IF (ENEW.LT.EOLD) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DEXP(-DE/T).GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF

      RETURN 
      ! }}}
      END 

