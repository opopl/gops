
SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)
! {{{
      USE COMMONS
      USE PORFUNCS

      IMPLICIT NONE

include mc.vars.inc.f90 ! Variable declarations 
include fmt.inc.f90     ! Formats
include mc.cien.inc.f90 ! Calculate the initial energy and save in EPREV
include mc.bh.inc.f90   ! Main basin-hopping loop
include write.111.inc.f90 ! write to LFH: Qu, E, Steps, RMS, Markov E, t (elapsed time)

WRITE(LFH,21)   ' Acceptance ratio for run=',  ARATIO,	&
                ' Step=',                      STEP, 	&
                ' Angular step factor=',       ASTEP,	&
                ' Temperature T=',             TEMP	
RETURN
! }}}
END

SUBROUTINE ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)
! {{{
USE COMMONS

IMPLICIT NONE

INTEGER,INTENT(INOUT) : NSUCCESS, NFAIL, NFAILT, NSUCCESST
DOUBLE PRECISION ARATIO,FAC
COMMON /EV/ EVAP, EVAPREJECT

ARATIO=1.D0*NSUCCESS/(1.D0*(NSUCCESS+NFAIL))

IF (ARATIO.GT.ACCRAT) THEN
	FAC=FAC0
ELSE
	FAC=1.D0/FAC0
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
END

SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,MCTEMP)
! {{{
      IMPLICIT NONE

      DOUBLE PRECISION,INTENT(IN) :: ENEW, EOLD, MCTEMP
      LOGICAL,INTENT(OUT) :: ATEST
      DOUBLE PRECISION, INTENT(OUT) :: RANDOM

      DOUBLE PRECISION :: DE,T 

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

