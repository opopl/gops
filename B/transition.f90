
      SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,NP,RANDOM,MCTEMP)
      ! {{{
      USE COMMONS
      USE V
      USE F

      IMPLICIT NONE
      ! sub 
      DOUBLE PRECISION,INTENT(IN) :: ENEW, EOLD, RANDOM, MCTEMP
      LOGICAL,INTENT(OUT) :: ATEST
      INTEGER,INTENT(IN) :: NP

      ! loc
      DOUBLE PRECISION :: DUMMY, DPRAND, EREF
      DOUBLE PRECISION DISTMIN, DISTMINOLD
      LOGICAL EVAP, EVAPREJECT
      COMMON /EV/ EVAP, EVAPREJECT

!  Standard canonical sampling.
!
         IF (ENEW.LT.EOLD) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DEXP(-(ENEW-EOLD)/MAX(MCTEMP,1.0D-100)).GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF

      RETURN 
      ! }}}
      END 

