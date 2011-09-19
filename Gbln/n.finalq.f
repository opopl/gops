      SUBROUTINE FINALQ
      USE COMMONS
      use QMODULE
      USE MODCHARMM, ONLY: ACESOLV,NCHENCALLS,ACEUPSTEP
      IMPLICIT NONE


      INTEGER J1, J2, ITERATIONS, BRUN,QDONE, J3
      DOUBLE PRECISION POTEL, SCREENC(3*NATOMS), X(3*NATOMS), ENERGY, GRAD(3*NATOMS), TIME
      DOUBLE PRECISION SAVECSMNORM, DUMMY2, DIST2, RMAT(3,3), AVVAL, CSMGRAD(3), XTEMP(1:3*NATOMS)
      DOUBLE PRECISION DUMMY(3*NATOMS), AA(3)

 
      SHELLMOVES(1:NPAR)=.FALSE.
      IF (CUTT) CUTOFF=FINALCUTOFF
      SAVEQ=.FALSE.
      NQ(1)=0
      IF (FIELDT) FIELDT=.FALSE.
      IF (SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
      ENDIF
      IF (SQUEEZET) SQUEEZET=.FALSE.
      MAXIT=MAXIT2
      DO J1=1,NSAVE
         IF (QMIN(J1).LT.1.0D10) THEN
            DO J2=1,3*NATOMS
            ENDDO
            NQ(1)=NQ(1)+1
            IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
            IF(UNFREEZEFINALQ) FROZEN(:)=.FALSE.  ! unfreeze everything before the final quenches
            WRITE(MYUNIT,'(A,I6,A,F20.10,A,I5,A,F15.7,A,F12.2)') 'Final Quench ',NQ(1),' energy=',
     1                POTEL,' steps=',ITERATIONS,' RMS force=',RMS,' time=',TIME-TSTART

            QMIN(J1)=POTEL
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=COORDS(J2,1)
            ENDDO

            IF (CSMT) THEN
               DO J2=1,CSMGPINDEX
                  DO J3=1,NATOMS
                     XTEMP(3*(J3-1)+1)=CSMPMAT(1,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1)
     &                                +CSMPMAT(1,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2)
     &                                +CSMPMAT(1,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                     XTEMP(3*(J3-1)+2)=CSMPMAT(2,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1)
     &                                +CSMPMAT(2,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2)
     &                                +CSMPMAT(2,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                     XTEMP(3*(J3-1)+3)=CSMPMAT(3,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1)
     &                                +CSMPMAT(3,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2)
     &                                +CSMPMAT(3,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                  ENDDO
               ENDDO
               IF (PERMDIST) THEN
                  DO J2=1,CSMGPINDEX
                     XTEMP(1:3*NATOMS)=CSMAV(1:3*NATOMS)
                  ENDDO
               ELSE
                  DO J2=1,CSMGPINDEX
                  ENDDO
               ENDIF
               AA(1)=0.0D0; AA(2)=0.0D0; AA(3)=6.283185307D0 ! should give an identity matrix
               SAVECSMNORM=CSMNORM
               DO J2=1,NATOMS
               ENDDO
               IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'finalq> CSM for averaged structure=',AVVAL
               QMINAV(J1)=AVVAL
               IF (DEBUG) WRITE(MYUNIT,'(A,I6,2G20.10)') 'finalq> J1,QMIN,QMINAV=',J1,QMIN(J1),QMINAV(J1)
               QMINPCSMAV(J1,1:3*NATOMS)=CSMAV(1:3*NATOMS)
            ENDIF
         ENDIF
      ENDDO
      NSAVE=NQ(1)
       
      IF (SORTT) THEN
         DO J1=1,NSAVE
            IF (QMIN(J1).LT.0.0D0) THEN
               DO J2=1,3*NATOMS
                  X(J2)=QMINP(J1,J2)
               ENDDO
               DO J2=1,NATOMS
                  VAT(J2,1)=VT(J2)
               ENDDO
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=COORDS(J2,1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF


         IF (TABOOT) THEN
            IF (NPAR.GT.1) THEN
               WRITE(MYUNIT,'(A)') 'Taboo lists:'
               DO J1=1,NPAR
                  WRITE(MYUNIT,'(A,G20.10)') 'Parallel run ',J1
                  WRITE(MYUNIT,'(6F15.7)') (ESAVE(J2,J1),J2=1,NT(J1))
               ENDDO
            ELSE
               WRITE(MYUNIT,'(A)') 'Taboo list:'
               WRITE(MYUNIT,'(6F15.7)') (ESAVE(J2,1),J2=1,NT(1))
            ENDIF
         ENDIF
      RETURN
      END
