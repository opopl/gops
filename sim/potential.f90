
      SUBROUTINE POTENTIAL(X,GRAD,EREAL,GRADT,SECT,PTYPE)

      USE COMMONS
      USE PORFUNCS
      USE BLN

      IMPLICIT NONE

      ! subroutine parameters 

      DOUBLE PRECISION,INTENT(OUT) :: EREAL, GRAD(*)
      DOUBLE PRECISION, INTENT(IN) :: X(*)
      LOGICAL GRADT,SECT

      ! potential type

      CHARACTER(LEN=*) PTYPE
     
      SELECTCASE(PTYPE)
        CASE("P46") ; CALL P46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
        CASE("G46") ; CALL G46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
        CASE("BLN") ; CALL BLN(X,GRAD,EREAL,GRADT)
      ENDSELECT

C  --------------- End of possible potentials - now add fields if required ------------------------------

      IF (PULLT) THEN
         EREAL=EREAL-PFORCE*(X(3*(PATOM1-1)+3)-X(3*(PATOM2-1)+3))
         GRAD(3*(PATOM1-1)+3)=GRAD(3*(PATOM1-1)+3)-PFORCE
         GRAD(3*(PATOM2-1)+3)=GRAD(3*(PATOM2-1)+3)+PFORCE
      ENDIF

      RETURN

      END
