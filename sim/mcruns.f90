
      SUBROUTINE MCRUNS(SCREENC)

      DOUBLE PRECISION SCREENC(:,:)

      INTEGER J1

      CALL MC(MCSTEPS,TFAC,SCREENC)

      !CALL FINALQ
      !CALL FINALIO

      RETURN
      END SUBROUTINE
