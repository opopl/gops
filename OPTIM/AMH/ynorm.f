
C     ---------------------  YNORM ----------------------

      SUBROUTINE YNORM(AMHMAXSIZ,JSTRT,JFINS,OBSERV,RNORM)

C     ---------------------------------------------------

C     YNORM  NORMALIZES THE ARRAY OBSERV BY RNORM

C     ARGUMENTS:

C        AMHMAXSIZ- MAXIMUM ARRAY LENGTH
C        JSTRT - FIRST SITE TO BE INCLUDED
C        JFINS - LAST SITE TO BE INCLUDED
C        OBSERV- SET OF OBSERVABLES TO BE NORMALIZED
C        RNORM - NORMALIZING FACTOR

C     ---------------------------------------------------

      IMPLICIT NONE

C      INCLUDE 'UTILITY'  ! UTILITY FILE

C     ARGUMENT DECLARATIONS:

         INTEGER AMHMAXSIZ,JSTRT,JFINS
     
         DOUBLE PRECISION OBSERV(AMHMAXSIZ),RNORM

C     INTERNAL VARIABLES:
          INTEGER I_RES


C     --------------------- BEGIN -----------------------

C     NORMALIZE OBSERV BY RNORM

      DO 515 I_RES=JSTRT,JFINS
         OBSERV(I_RES)=RNORM*OBSERV(I_RES)
  515 CONTINUE

C     --------------------- DONE  -----------------------
      
      RETURN
      END
