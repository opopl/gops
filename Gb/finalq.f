
      SUBROUTINE FINALQ
!op226> Declarations {{{ 

      USE COMMONS
      USE QMODULE

      IMPLICIT NONE

      INTEGER J1, J2, ITERATIONS, BRUN,QDONE, J3
      DOUBLE PRECISION POTEL, SCREENC(3*NATOMS), X(3*NATOMS), ENERGY, GRAD(3*NATOMS), TIME
      DOUBLE PRECISION SAVECSMNORM, DUMMY2, DIST2, RMAT(3,3), AVVAL, CSMGRAD(3), XTEMP(1:3*NATOMS)
      DOUBLE PRECISION DUMMY(3*NATOMS), AA(3)

      COMMON /MYPOT/ POTEL
!op226> End declarations }}} 
 
      NQ(1)=0
      IF (SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
      ENDIF
      DO J1=1,NSAVE
      ! {{{
         IF (QMIN(J1).LT.1.0D10) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,1)=QMINP(J1,J2)
            ENDDO
            NQ(1)=NQ(1)+1
            CALL QUENCH(.TRUE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            WRITE(LFH,'(A,I6,A,F20.10,A,I5,A,F15.7,A,F12.2)') 'Final Quench ',NQ(1),' energy=',
     1                POTEL,' steps=',ITERATIONS,' RMS force=',RMS,' time=',TIME-TSTART

            QMIN(J1)=POTEL
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=COORDS(J2,1)
            ENDDO
         ENDIF
         ! }}}
      ENDDO

      NSAVE=NQ(1)
      CALL GSORT2(NSAVE,NATOMS)

C  Optionally sort the atoms from most bound to least bound according to VAT. {{{
C
       
!op226>}}} 

      RETURN
      END
