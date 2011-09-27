
      SUBROUTINE FINALQ
!op226>=================================== 
!op226> Declarations {{{ 
      USE COMMONS
      USE V
      USE F

      IMPLICIT NONE

      INTEGER J1, J2, ITERATIONS, BRUN,QDONE, J3
      DOUBLE PRECISION ::   POTEL,TIME,ENERGY
      DOUBLE PRECISION,DIMENSION(NR):: SCREENC, X, GRAD, xtemp

      COMMON /MYPOT/ POTEL
!op226> End declarations }}} 
!op226>=================================== 
 
!
!  Make sure the lowest minima are tightly converged and then sort
!  them just to be on the safe side.
!
      SAVEQ=.FALSE.
      NQ(1)=0
      IF (SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
      ENDIF
      MAXIT=MAXIT2

      DO J1=1,NSAVE
         IF (QMIN(J1).LT.1.0D10) THEN
            COORDS(1:NR,1)=QMINP(J1,1:NR)
            NQ(1)=NQ(1)+1
            CALL QUENCH(.TRUE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            WRITE(LFH,'(A,I6,A,F20.10,A,I5,A,F15.7,A,F12.2)') &
			&                'Final Quench ',NQ(1),&
			&                ' energy=', POTEL,&
			&                ' steps=',ITERATIONS,&
			&                ' RMS force=',RMS,&
			&                ' time=',TIME-TSTART

            QMIN(J1)=POTEL
            EAMIN(J1,1:6)=EA(1:6)
            QMINP(J1,1:NR)=COORDS(1:NR,1)
         ENDIF
      ENDDO

      NSAVE=NQ(1)
      CALL GSORT2(NSAVE,NATOMS)
!
!  Optionally sort the atoms from most bound to least bound according to VAT. {{{
!
       
      IF (SORTT) THEN
         DO J1=1,NSAVE
            IF (QMIN(J1).LT.0.0D0) THEN
               DO J2=1,3*NATOMS
                  COORDS(J2,1)=QMINP(J1,J2)
                  X(J2)=QMINP(J1,J2)
               ENDDO
               CALL POTENTIAL(X,GRAD,ENERGY,.FALSE.,.FALSE.)
               DO J2=1,NATOMS
                  VAT(J2,1)=VT(J2)
               ENDDO
               CALL SORT3(NATOMS,NATOMS,VAT(1:NATOMS,1),COORDS(1:3*NATOMS,1))
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=COORDS(J2,1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

!op226>}}} 

!     IF (DEBUG) THEN
         IF (TABOOT) THEN
            IF (NPAR.GT.1) THEN
               WRITE(LFH,'(A)') 'Taboo lists:'
               DO J1=1,NPAR
                  WRITE(LFH,'(A,G20.10)') 'Parallel run ',J1
                  WRITE(LFH,'(6F15.7)') (ESAVE(J2,J1),J2=1,NT(J1))
               ENDDO
            ELSE
               WRITE(LFH,'(A)') 'Taboo list:'
               WRITE(LFH,'(6F15.7)') (ESAVE(J2,1),J2=1,NT(1))
            ENDIF
         ENDIF
!     ENDIF
      RETURN
      END
