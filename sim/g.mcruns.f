
      SUBROUTINE MCRUNS(SCREENC)

      USE COMMONS

      IMPLICIT NONE
      LOGICAL LOPEN
      DOUBLE PRECISION SCREENC(3*NATOMS)

      INTEGER J1

      CALL MC(MCSTEPS(1),TFAC(1),SCREENC)

      CALL FINALQ
      CALL FINALIO

      RETURN
      END