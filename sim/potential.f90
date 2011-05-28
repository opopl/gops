
!> @param[in]  X        dp(N,3)    input coordinates 
!> @param[out] GRAD     dp(N,3)    gradient 
!> @param[out] EREAL    dp         energy 

      SUBROUTINE POTENTIAL(R,GRAD,EREAL,GRADT,SECT,PTYPE)

      USE COMMONS
      USE PORFUNCS
      USE BLN

      IMPLICIT NONE

      ! subroutine parameters 

      DOUBLE PRECISION,INTENT(OUT) :: EREAL, GRAD(:,3)
      DOUBLE PRECISION,INTENT(IN) :: R(:,3)

      LOGICAL GRADT,SECT

      ! potential type

      CHARACTER(LEN=*) PTYPE
    
      CALL EBLN(N,R,GRAD,ENERGY,HESS,PTYPE,.TRUE.,.TRUE.)

      IF (PULLT) THEN
         EREAL=EREAL-PFORCE*R(PATOM1,3)+PFORCE*R(PATOM2,3)
         GRAD(PATOM1,3)=GRAD(PATOM1,3)-PFORCE
         GRAD(PATOM2,3)=GRAD(PATOM2,3)+PFORCE
      ENDIF

      RMS=MAX(DSQRT(SUM(GRAD**2)/(3*NATOMS)),1.0D-100)

      RETURN

      END
