
!> @param[in]  X        dp(3*N)    input coordinates 
!> @param[out] GRAD     dp(3*N)    gradient 
!> @param[out] EREAL    dp(3*N)    energy 

      SUBROUTINE POTENTIAL(X,GRAD,EREAL,GRADT,SECT,PTYPE)

      USE COMMONS
      USE PORFUNCS
      USE BLN

      IMPLICIT NONE

      ! subroutine parameters 

      DOUBLE PRECISION,INTENT(OUT) :: EREAL, GRAD(*)
      DOUBLE PRECISION,INTENT(IN) :: X(*)
      LOGICAL GRADT,SECT

      ! potential type

      CHARACTER(LEN=*) PTYPE
    
      CALL EBLN(N,QO,GRAD,ENERGY,HESS,PTYPE)

C  --------------- End of possible potentials - now add fields if required ------------------------------

      IF (PULLT) THEN
         EREAL=EREAL-PFORCE*X(3*PATOM1)-PFORCE*X(3*PATOM2)
         GRAD(3*PATOM1)=GRAD(3*PATOM1)-PFORCE
         GRAD(3*PATOM2)=GRAD(3*PATOM2)+PFORCE
      ENDIF

      RETURN

      END
