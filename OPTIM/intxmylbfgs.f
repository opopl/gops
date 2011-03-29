C   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF OPTIM.
C
C   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** JULY 1990 ***
C
C        LINE SEARCH REMOVED PLUS SMALL MODIFICATIONS, DJW 2001
C        THIS ROUTINE IS FOR MINIMISING A RAYLEIGH-RITZ RATIO
C        IN BFGS/BFGS TRANSITION STATE SEARCHES.
C        CALLED IF INTMINT IS TRUE, COORDS PASSED FROM INTBFGSTS 
C        VIA INTBEIG CONTAINS INTERNALS, N IS NINTS. JMC
C
      SUBROUTINE INTXMYLBFGS(N,M,X,DIAGCO,DIAG,EPS,W,EENERGY,COORDS,ITDONE,ITMAX,ENERGY,PTEST,CONVERGED)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE PORFUNCS
      IMPLICIT NONE

      INTEGER N,M,J1,J2,ITMAX,ITDONE,NFAIL
C JMC      DOUBLE PRECISION X(N),G(3*NATOMS),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT,DUMMY2
      DOUBLE PRECISION X(3*NATOMS),G(3*NATOMS),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT,DUMMY2
      DOUBLE PRECISION EPS,DUMMY1,ENERGY,ENEW,GNEW(3*NATOMS),RMS,EENERGY
      LOGICAL DIAGCO, PTEST
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,COORDS(3*NATOMS),OVERLAP,DOT1,DOT2
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE, ISTAT
      LOGICAL MFLAG, FAILED, CONVERGED

      ITER=0
      ITDONE=0
      NFAIL=0
      FAILED=.FALSE.

C JMC X IS DISPLACEMENT VECTOR (CONVERGING TO EIGENVECTOR), COORDS IS COORDS, EENERGY IS REAL ENERGY, 
C G IS GRADIENT OF LAMBDA WRT X, ENERGY IS (ALMOST) LAMBDA, RMS IS RMS G

      CALL INTSECDIAG(X,COORDS,EENERGY,G,ENERGY,.TRUE.,RMS)

      IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A)') 'INTXMYLBFGS> EIGENVALUE AND RMS FORCE=',ENERGY,RMS,' AFTER ',ITDONE,' STEPS'

10    CALL FLUSH(6,ISTAT)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
         MFLAG=.TRUE.
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            IF (PTEST) WRITE(*,'(A,I4,A,F15.7,A,F15.7)') 'INTXMYLBFGS> SMALLEST EIGENVALUE CONVERGED IN ',ITDONE,
     1              ' STEPS. EIGENVALUE=',ENERGY,' RMS FORCE=',RMS
            CALL FLUSH(6,ISTAT)
            IF (PTEST.AND.(ITER.GT.0))
     1          WRITE(*,'(A,F20.10)') 'INTXMYLBFGS> DIAGONAL INVERSE HESSIAN ELEMENTS ARE NOW ',DIAG(1)
            CONVERGED=.TRUE.
            RETURN
         ENDIF
      ENDIF

      IF ((ITDONE.EQ.ITMAX).OR.FAILED) THEN
         FIXIMAGE=.FALSE.
         IF (PTEST) WRITE(*,'(A,G15.7,A,G15.7)') 'INTXMYLBFGS> **WARNING - SMALLEST EIGENVALUE DID NOT CONVERGE, VALUE=',ENERGY,
     1           ' RMS FORCE=',RMS
         IF (PTEST) WRITE(*,'(A,F20.10)') 'INTXMYLBFGS> DIAGONAL INVERSE HESSIAN ELEMENTS ARE NOW ',DIAG(1)
         CONVERGED=.FALSE.
         RETURN
      ENDIF

      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(*,240)
 240        FORMAT('INTXMYLBFGS> IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            PRINT*,'USING ESTIMATE OF THE INVERSE DIAGONAL ELEMENTS'
            DO I=1,N
               IF (DIAG(I).LE.0.0D0) THEN
                  WRITE(*,235) I
 235              FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
     1                   ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
                  STOP
               ENDIF
            ENDDO
         ELSE
C           INQUIRE(FILE='DIAG',EXIST=YESNO)
C           IF (YESNO) THEN
C              OPEN(UNIT=34,FILE='DIAG',STATUS='OLD')
C              READ(34,*) (DIAG(I),I=1,N)
C              PRINT*,'INTXMYLBFGS> DIAG READ IN LBFGS'
C              WRITE(*,'(6F15.5)') (DIAG(I),I=1,N)
C           ELSE
            DO I=1,N
               DIAG(I)=XDGUESS
            ENDDO
         ENDIF
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
         ISPT= N+2*M
         IYPT= ISPT+N*M
C
C  NR STEP FOR DIAGONAL INVERSE HESSIAN
C
         DO I=1,N
            W(ISPT+I)= -G(I)*DIAG(I)
            W(I)= -G(I)*DIAG(I)
         ENDDO
         GNORM= DSQRT(DDOT(N,G,1,G,1))
C
C  MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS.
C
         IF (GNORM.EQ.0.0D0) THEN
            GNORM=1.0D0 ! EXACT ZERO IS PRESUMABLY WRONG!
            PRINT '(A)','WARNING - GNORM WAS ZERO IN XMYLBFGS, RESETTING TO ONE'
         ENDIF
         STP=MIN(GNORM,1.0D0/GNORM)
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
         IF (YS.EQ.0.0D0) YS=1.0D0
C
C  UPDATE ESTIMATE OF DIAGONAL INVERSE HESSIAN ELEMENTS
C  WE DIVIDE BY BOTH YS AND YY AT DIFFERENT POINTS, SO
C  THEY HAD BETTER NOT BE ZERO!
C
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            IF (YY.EQ.0.0D0) YY=1.0D0
            DO I=1,N
               DIAG(I)= YS/YY
C              DIAG(I)= ABS(YS/YY)
            ENDDO
         ELSE
            PRINT*,'USING ESTIMATE OF THE INVERSE DIAGONAL ELEMENTS'
            DO I=1,N
               IF (DIAG(I).LE.0.0D0) THEN
                  WRITE(*,235) I
                  STOP
               ENDIF
            ENDDO
         ENDIF
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: NOCEDAL, J. 1980,
C     "UPDATING QUASI-NEWTON MATRICES WITH LIMITED STORAGE",
C     MATHEMATICS OF COMPUTATION, VOL.24, NO.151, PP. 773-782.
C     ---------------------------------------------------------
C
         CP= POINT
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS
         DO I=1,N
            W(I)= -G(I)
         ENDDO
         CP= POINT
         DO I= 1,BOUND
            CP=CP-1
            IF (CP.EQ. -1)CP=M-1
            SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)= W(N+CP+1)*SQ
            CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
         ENDDO
        
         DO I=1,N
            W(I)=DIAG(I)*W(I)
         ENDDO

         DO I=1,BOUND
            YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
            BETA= W(N+CP+1)*YR
            INMC=N+M+CP+1
            BETA= W(INMC)-BETA
            ISCN=ISPT+CP*N
            CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
            CP=CP+1
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0
      ENDIF
C
C  ORTHOGONALISATION OF THE SEARCH DIRECTION TO TRANSLATIONS AND ROTATIONS
C  APPEARS TO MAKE NO DIFFERENCE.
C
C      CALL ORTHOGOPT(W,COORDS,.TRUE.)

C
C  STORE THE NEW SEARCH DIRECTION
C
C     IF (ITER.NE.0) THEN
         DO I=1,N
            W(ISPT+POINT*N+I)= W(I)
         ENDDO
C     ENDIF

C     OVERLAP=DDOT(N,G,1,W,1)/SQRT(DDOT(N,G,1,G,1)*DDOT(N,W,1,W,1))
      DOT1=SQRT(DDOT(N,G,1,G,1))
      DOT2=SQRT(DDOT(N,W,1,W,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

      IF (OVERLAP.GT.0.0D0) THEN
         IF (PTEST) PRINT*,'SEARCH DIRECTION HAS POSITIVE PROJECTION ONTO GRADIENT - REVERSING STEP'
         DO I=1,N
            W(ISPT+POINT*N+I)= -W(I)
         ENDDO
C        ITER=0
C        GOTO 10
      ENDIF
      
      DO I=1,N
         W(I)=G(I)
      ENDDO
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
      IF (STP*SLENGTH.GT.MAXXBFGS) STP=MAXXBFGS/SLENGTH
C
C  WE NOW HAVE THE PROPOSED STEP.
C
      DO J1=1,N
         X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
      ENDDO 
      NDECREASE=0
20    CONTINUE
      CALL INTSECDIAG(X,COORDS,EENERGY,GNEW,ENEW,.TRUE.,RMS)

      IF ((ENEW-ENERGY)/MAX(ABS(ENEW),1.0D-100).LE.1.0D-3) THEN
C     IF ((ENEW-ENERGY)/MAX(ABS(ENEW),1.0D-100).LE.1.0D-8) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         DO J1=1,N
            G(J1)=GNEW(J1)
         ENDDO
         IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A,F13.10)') 
     1          ' EIGENVALUE AND RMS FORCE=',ENERGY,RMS,' AFTER ',ITDONE,' STEPS, STEP:',STP*SLENGTH
      ELSE
C
C  EIGENVALUE INCREASED - TRY AGAIN WITH A SMALLER STEP SIZE?
C
         IF (NDECREASE.GT.5) THEN
C JMC         IF (NDECREASE.GT.10) THEN
C           PRINT*,' LBFGS STEP CANNOT FIND A LOWER NON-ZERO EIGENVALUE'
            NFAIL=NFAIL+1
            ITER=0  !  TRY RESETTING
            DO J1=1,N
               X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
            ENDDO
            IF (NFAIL.GT.2) FAILED=.TRUE. ! ALLOW ONE RESTART DJW
            GOTO 10
         ENDIF
         DO J1=1,N
            X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
         ENDDO 
         NDECREASE=NDECREASE+1
         STP=STP/10.0D0
         IF (PTEST) WRITE(*,'(A,F19.10,A,G16.10,A,F15.8)') 
     1                  ' EIGENVALUE INCREASED FROM ',ENERGY,' TO ',ENEW,' DECREASING STEP TO ',STP*SLENGTH
         FIXIMAGE=.TRUE.
         GOTO 20
      ENDIF
C
C     COMPUTE THE NEW STEP AND GRADIENT CHANGE
C
      NPT=POINT*N
      DO I=1,N
         W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
         W(IYPT+NPT+I)= G(I)-W(I)
      ENDDO
C
C  DOES IT HELP TO KEEP X NORMALISED?
C
      CALL VECNORM(X,N)

      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      FIXIMAGE=.FALSE.
      IF ((FIXAFTER.GT.0).AND.(ITDONE.GE.FIXAFTER)) FIXIMAGE=.TRUE.
      GOTO 10

      RETURN
      END
