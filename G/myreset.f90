      ! MYRESET {{{

      SUBROUTINE MYRESET(JP,NATOMS,NPAR,NSEED)
      USE COMMONS,ONLY : LFH,COORDS,COORDSO,VAT,VATO
      IMPLICIT NONE
      INTEGER JP, NSEED, J2, NATOMS, NPAR

      ! sub 
      INTEGER JP,NATOMS,NPAR,NSEED

      DO J2=1,3*(NATOMS-NSEED)
         COORDS(J2,JP)=COORDSO(J2,JP)
      ENDDO
      DO J2=1,NATOMS
         VAT(J2,JP)=VATO(J2,JP)
      ENDDO

      RETURN
      END
      ! }}}

