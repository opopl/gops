
      SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,NP,RANDOM,MCTEMP)
      ! {{{
      USE COMMONS
      USE V
      USE F

      IMPLICIT NONE
      ! sub 
      DOUBLE PRECISION ENEW, EOLD, RANDOM, MCTEMP
      LOGICAL ATEST
      INTEGER NP

      ! loc
      DOUBLE PRECISION :: DUMMY, DPRAND, EREF, TEOLD, TENEW, RATIO
      DOUBLE PRECISION TRANS, DISTMIN, DISTMINOLD
      LOGICAL ATEST, FLAT, evap, evapreject
      INTEGER NP,INDEXOLD, INDEXNEW, J1, NDUMMY
      DATA DISTMINOLD /0.0D0/
      COMMON /DMIN/ DISTMIN
      common /ev/ evap, evapreject

      IF (DISTMINOLD.EQ.0.0D0) DISTMINOLD=DISTMIN  ! this should allow for the first step
         TEOLD=EOLD
         TENEW=ENEW

!  Standard canonical sampling.
!
         IF (TENEW.LT.TEOLD) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DEXP(-(TENEW-TEOLD)/MAX(MCTEMP,1.0D-100)).GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF

      RETURN 
      ! }}}
      END 

