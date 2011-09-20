
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
 
      SAVEQ=.FALSE.
      NQ(1)=0
      IF (SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
      ENDIF
      MAXIT=MAXIT2
      DO J1=1,NSAVE
      ! {{{
         IF (QMIN(J1).LT.1.0D10) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,1)=QMINP(J1,J2)
            ENDDO
            NQ(1)=NQ(1)+1
            IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
            IF(UNFREEZEFINALQ) FROZEN(:)=.FALSE.  ! unfreeze everything before the final quenches
            CALL QUENCH(.TRUE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            WRITE(MYUNIT,'(A,I6,A,F20.10,A,I5,A,F15.7,A,F12.2)') 'Final Quench ',NQ(1),' energy=',
     1                POTEL,' steps=',ITERATIONS,' RMS force=',RMS,' time=',TIME-TSTART

            QMIN(J1)=POTEL
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=COORDS(J2,1)
            ENDDO

!op226> IF (CSMT) THEN {{{ 
            IF (CSMT) THEN
               CSMAV(1:3*NATOMS)=0.0D0
               DO J2=1,CSMGPINDEX
!
! rotate permuted image to best orientation with CSMPMAT
! apply point group operation J2
! 
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
                  CALL CSMROT(XTEMP,DUMMY,1,J2)
                  CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)+DUMMY(1:3*NATOMS)
               ENDDO
               CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)/CSMGPINDEX
!
!  Check the CSM for the averaged structure. It should be zero if this structure has the
!  right point group. Need to reset CSMIMAGES and CSMNORM temporarily. 
!  
               IF (PERMDIST) THEN
                  DO J2=1,CSMGPINDEX
                     XTEMP(1:3*NATOMS)=CSMAV(1:3*NATOMS)
                     CALL CSMROT(XTEMP,DUMMY,1,J2)
                     CALL MINPERMDIST(XTEMP,DUMMY,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY2,DIST2,RIGID,RMAT)
                     CALL CSMROT(DUMMY,XTEMP,-1,J2) ! need to rotate the permuted rotated images back to the reference orientation
                     CSMIMAGES(1+3*NATOMS*(J2-1):3*NATOMS*J2)=XTEMP(1:3*NATOMS)
                  ENDDO
               ELSE
                  DO J2=1,CSMGPINDEX
                     CSMIMAGES(1+3*NATOMS*(J2-1):3*NATOMS*J2)=CSMAV(1:3*NATOMS)
                  ENDDO
               ENDIF
               AA(1)=0.0D0; AA(2)=0.0D0; AA(3)=6.283185307D0 ! should give an identity matrix
               SAVECSMNORM=CSMNORM
               CSMNORM=0.0D0
               DO J2=1,NATOMS
                  CSMNORM=CSMNORM+CSMAV(3*(J2-1)+1)**2+CSMAV(3*(J2-1)+2)**2+CSMAV(3*(J2-1)+3)**2
               ENDDO
               CSMNORM=2*CSMGPINDEX*CSMNORM
               CALL CSMPOTGRAD(CSMAV,AA,AVVAL,.TRUE.,CSMGRAD)
               CSMNORM=SAVECSMNORM
               IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'finalq> CSM for averaged structure=',AVVAL
               QMINAV(J1)=AVVAL
               IF (DEBUG) WRITE(MYUNIT,'(A,I6,2G20.10)') 'finalq> J1,QMIN,QMINAV=',J1,QMIN(J1),QMINAV(J1)
               QMINPCSMAV(J1,1:3*NATOMS)=CSMAV(1:3*NATOMS)
            ENDIF
!op226>}}} 
         ENDIF
         ! }}}
      ENDDO

      NSAVE=NQ(1)
      CALL GSORT2(NSAVE,NATOMS)

C  Optionally sort the atoms from most bound to least bound according to VAT. {{{
C
       
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

      RETURN
      END
