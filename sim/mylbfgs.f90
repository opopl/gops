! ------------------------------------------------------------------
! Doxygen {{{
!>
!> @name MYLBFGS
!
!> @brief Adaptation of the original LBFGS routine
!
!> @param[in] N       
!>             is an INTEGER variable that must be set by the user to the
!>             number of variables. It is not altered by the routine.
!>             Restriction: N>0.
!
!> @param[in] DIAGCO 
!>             is a LOGICAL variable that must be set to .TRUE. if the
!>             user wishes to provide the diagonal matrix Hk0 at each
!>             iteration. Otherwise it should be set to .FALSE., in which
!>             case  LBFGS will use a default value described below. If
!>             DIAGCO is set to .TRUE. the routine will return at each
!>             iteration of the algorithm with IFLAG=2, and the diagonal
!>             matrix Hk0  must be provided in the array DIAG.
!
!> @param[in] EPS 
!>             is a positive DOUBLE PRECISION variable that must be set by
!>             the user, and determines the accuracy with which the solution
!>             is to be found. The subroutine terminates when
!>
!>                         ||G|| < EPS max(1,||X||),
!>
!>             where ||.|| denotes the Euclidean norm.
!
!> @param[inout] XCOORDS is a DOUBLE PRECISION array of length N. On initial entry
!>             it must be set by the user to the values of the initial
!>             estimate of the solution vector. On exit with IFLAG=0, it
!>             contains the values of the variables at the best point
!>             found (usually a solution).
!
!
! }}}
! ------------------------------------------------------------------

      SUBROUTINE MYLBFGS(N,M,XCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET)
      
      USE COMMONS
      USE PORFUNCS
      
      IMPLICIT NONE
      ! ====================================
      ! LBFGS comparisons
      ! ====================================
      ! {{{
      !
      ! LBFGS           MY
      !
      ! dp X(N)         dp XCOORDS(3*NATOMS)
      ! dp DIAG(N)      
      ! i IPRINT
      !                 i ITMAX, ITDONE
      ! dp XTOL 
      ! dp EPS          dp EPS
      !                 dp ENERGY 
      !                 L RESET 
      !                 L MFLAG
      ! L DIAGCO        L DIAGCO 
      ! dp W(N)
      !       
      ! }}}
      ! ====================================
      ! original LBFGS declaration
      ! ====================================
      ! {{{
      !
      ! SUBROUTINE LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
      ! SUBROUTINE MY...(N,M,XCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET)
      !
      !

      ! INTEGER N,M,IPRINT(2),IFLAG
      ! DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
      ! DOUBLE PRECISION F,EPS,XTOL
      ! LOGICAL DIAGCO
      !
      ! }}}
      ! ====================================
      ! subroutine parameters 
      INTEGER, INTENT(IN) :: N,M
      DOUBLE PRECISION, INTENT(INOUT) XCOORDS(3*NATOMS)
      LOGICAL,INTENT(IN) :: DIAGCO
      DOUBLE PRECISION,INTENT(IN) :: EPS
      LOGICAL,INTENT(OUT) :: MFLAG
      DOUBLE PRECISION,INTENT(OUT) :: ENERGY
      INTEGER,INTENT(IN) :: ITMAX
      INTEGER,INTENT(OUT) :: ITDONE
      LOGICAL,INTENT(IN) RESET

      ! local parameters 
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W
      INTEGER ITER
      SAVE W,DIAG,ITER
      !
      integer   ISPT ! index for storage of search steps
      integer   IYPT ! index for storage of gradient differences
      !
      INTEGER J1,J2,J3,NFAIL,NDECREASE,NGUESS,NDUMMY
      DOUBLE PRECISION GRAD(3*NATOMS),SLENGTH,DDOT,EPLUS,EMINUS,DIFF,DUMMY,WTEMP(3*NATOMS)
      DOUBLE PRECISION TMPANG(3*NATOMS), TMPCOORDS(3*NATOMS)
      DOUBLE PRECISION ENEW,GNEW(3*NATOMS),OVERLAP,OLDX(3*NATOMS),OLDOLDX(3*NATOMS),VGUESS(3),
      1                 X1, Y1, Z1, X2, Y2, Z2, TRY(3*NATOMS), D1, D2, RBCOORDS(18), DUMMY2, DIST, DIST1
      LOGICAL YESNO, NOTCALLED, CTEST
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,POTEL,QSTART,QFINISH
      DOUBLE PRECISION OLDCART(3*NATOMS), DELTAQ(N),DELTACART(3*NATOMS),LEPSILON,DOT1,DOT2
      DOUBLE PRECISION LCART(3*NATOMS),OLDQ(N),NEWQ(N),OLDGINT(N),GINT(N),XINT(N),XSAVE(N),SMINKCURRENTP
      DOUBLE PRECISION, ALLOCATABLE :: FRAMES(:,:), PE(:), MODGRAD(:)
      LOGICAL NOCOOR, FAILED, COREDONE
      INTEGER ITER,POINT,BOUND,NPT,CP,INMC,IYCN,ISCN
      INTEGER KD, NNZ
      COMMON /MYPOT/ POTEL
      COMMON /Q4C/ QSTART, QFINISH
      LOGICAL EVAP, GUIDECHANGET, GUIDET, EVAPREJECT, SMINKCHANGET, CSMDOGUIDET
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /EV/ EVAP, EVAPREJECT
      SAVE POINT, ISPT, IYPT, NPT

      IF (.NOT.ALLOCATED(DIAG)) ALLOCATE(DIAG(N))       ! SAVE doesn't work otherwise for Sun
      IF (.NOT.ALLOCATED(W)) ALLOCATE(W(N*(2*M+1)+2*M)) ! SAVE doesn't work otherwise for Sun

      IF (SIZE(W,1).NE.N*(2*M+1)+2*M) THEN ! mustn't call mylbfgs with changing number of variables!!!
         WRITE(LOG_FH, '(A,I10,A,I10,A)') 'ERROR, dimension of W=',SIZE(W,1),' but N*(2*M+1)+2*M=',N*(2*M+1)+2*M,' in mylbfgs'
         call exit(10)
      ENDIF

      COREDONE=.FALSE.
      SMINKCHANGET=.FALSE.
      SMINKCURRENT=0.0D0
      SMINKCURRENTP=0.0D0
      LOCALSTEEREDMINT=.FALSE.
      IF (STEEREDMINT) LOCALSTEEREDMINT=.TRUE.

      NFAIL=0
      IF (GUIDECHANGET) ITER=0
      IF (RESET) ITER=0
      ITDONE=0
      FIXIMAGE=.FALSE.

      CALL POTENTIAL(XCOORDS,GRAD,ENERGY,.TRUE.,.FALSE.)

         IF (EVAPREJECT) RETURN
      POTEL=ENERGY

      IF (DEBUG) WRITE(LOG_FH,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'

C  Termination test. 
C
10    CALL FLUSH(LOG_FH)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
                     MFLAG=.TRUE.
         IF (EVAP) MFLAG=.FALSE. ! do not allow convergence if we happen to have a small RMS and EVAP is true'
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            IF (DEBUG) WRITE(LOG_FH,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'
            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
         IF (DEBUG) FIXIMAGE=.FALSE.
         IF (DEBUG) WRITE(LOG_FH,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
         RETURN
      ENDIF

      IF (ITER.EQ.0) THEN
         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            WRITE(LOG_FH,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(LOG_FH,235) J1
 235              FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
     1                   ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
                  STOP
               ENDIF
            ENDDO
         ELSE
            DO J1=1,N
               DIAG(J1)=DGUESS
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
         ISPT= N+2*M    ! index for storage of search steps
         IYPT= ISPT+N*M ! index for storage of gradient differences
C
C  NR step for diagonal inverse Hessian
C

C
C  Make the first guess for the step length cautious.
C
         STP=MIN(1.0D0/GNORM,GNORM)
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
C
C  Update estimate of diagonal inverse Hessian elements
C
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            IF (YY.EQ.0.0D0) THEN
               WRITE(LOG_FH,'(A)') 'WARNING, resetting YY to one in mylbfgs'
               YY=1.0D0
            ENDIF
            IF (YS.EQ.0.0D0) THEN
               WRITE(LOG_FH,'(A)') 'WARNING, resetting YS to one in mylbfgs'
               YS=1.0D0
            ENDIF
C           WRITE(*,'(A,2F20.10)') 'YS/YY,STP=',YS/YY,STP
            DO J1=1,N
C              DIAG(J1)= ABS(YS/YY) ! messes up after step reversals!
               DIAG(J1)= YS/YY
            ENDDO
         ELSE
            WRITE(LOG_FH,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(LOG_FH,235) J1
                  STOP
               ENDIF
            ENDDO
         ENDIF
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C
         CP= POINT
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS
            DO J1=1,N
               W(J1)= -GRAD(J1)
            ENDDO
         CP= POINT
         DO J1= 1,BOUND
            CP=CP-1
            IF (CP.EQ.-1) CP=M-1
            SQ=DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)=W(N+CP+1)*SQ
            CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
         ENDDO
        
         DO J1=1,N
            W(J1)=DIAG(J1)*W(J1)
         ENDDO

         DO J1=1,BOUND
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
C  Store the new search direction
C
      IF (ITER.GT.0) THEN
         DO J1=1,N
            W(ISPT+POINT*N+J1)= W(J1)
         ENDDO
      ENDIF

            DOT1=SQRT(DDOT(N,GRAD,1,GRAD,1))
C
C  Overflow has occasionally occurred here.
C  We only need the sign of the overlap, so use a temporary array with
C  reduced elements.
C
      DUMMY=1.0D0
      DO J1=1,N
         IF (ABS(W(J1)).GT.DUMMY) DUMMY=ABS(W(J1))
      ENDDO
      DO J1=1,N
         WTEMP(J1)=W(J1)/DUMMY
      ENDDO
      DOT2=SQRT(DDOT(N,WTEMP,1,WTEMP,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) THEN
                   OVERLAP=DDOT(N,GRAD,1,WTEMP,1)/(DOT1*DOT2)
      ENDIF
      IF (OVERLAP.GT.0.0D0) THEN
         IF (DEBUG) WRITE(LOG_FH,'(A)') 'Search direction has positive projection onto gradient - reversing step'
         DO J1=1,N
            W(ISPT+POINT*N+J1)= -W(J1)  !!! DJW, reverses step
         ENDDO
      ENDIF

         DO J1=1,N
            W(J1)=GRAD(J1)
         ENDDO
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
      IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH
C
C  We now have the proposed step.
C
! Save XCOORDS here so that we can undo the step reliably including the
! non-linear projection for Thomson for the angular coordinates.
!
         XSAVE(1:N)=XCOORDS(1:N) 
         DO J1=1,N
            XCOORDS(J1)=XCOORDS(J1)+STP*W(ISPT+POINT*N+J1)
         ENDDO 
C
C  For charmm internals must transform and back-transform!
C
      NDECREASE=0
      LEPSILON=1.0D-6


! csw34> INCREMENT THE FORCE CONSTANT FOR STEERED MINIMISATION
      SMINKCHANGET=.FALSE.
      IF (LOCALSTEEREDMINT) THEN
         SMINKCURRENT=MIN(SMINKCURRENT+SMINKINC,SMINK)
         IF (SMINKCURRENT.NE.SMINKCURRENTP) SMINKCHANGET=.TRUE.
! a bit of useful debug printing         
        IF (DEBUG) WRITE(LOG_FH,'(A,2F20.10,L5)') 'SMINKCURRENT,SMINKCURRENTP,SMINKCHANGET=',SMINKCURRENT,SMINKCURRENTP,SMINKCHANGET
      ENDIF

      CALL POTENTIAL(XCOORDS,GNEW,ENEW,.TRUE.,.FALSE.)

      IF (((ENEW-ENERGY.LE.MAXERISE).OR.EVAP.OR.GUIDECHANGET.OR.SMINKCHANGET).AND.(ENEW-ENERGY.GT.MAXEFALL)) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         DO J1=1,3*NATOMS
            GRAD(J1)=GNEW(J1)
         ENDDO
         IF (DEBUG) WRITE(LOG_FH,'(A,F20.10,G20.10,A,I6,A,F13.10)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,
     1           ' LBFGS steps, step:',STP*SLENGTH

C  May want to prevent the PE from falling too much if we are trying to visit all the
C  PE bins. Halve the step size until the energy decrease is in range.
C
      ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
C
C  Energy decreased too much - try again with a smaller step size
C
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            WRITE(LOG_FH,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find an energy in the required range, NFAIL=',NFAIL
      !
! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
               XCOORDS(1:N)=XSAVE(1:N)
               GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point

            ITER=0   !  try resetting
            IF (NFAIL.GT.20) THEN
               WRITE(LOG_FH,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.
               RETURN
            ENDIF
            GOTO 30
         ENDIF
!
! Resetting to XSAVE and adding half the step should be the same as subtracting 
! half the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
            XCOORDS(1:N)=XSAVE(1:N)
            DO J1=1,N
               XCOORDS(J1)=XCOORDS(J1)+0.5*STP*W(ISPT+POINT*N+J1)
            ENDDO 
         ENDIF
         STP=STP/2.0D0
         NDECREASE=NDECREASE+1
         IF (DEBUG) WRITE(LOG_FH,'(A,F19.10,A,F16.10,A,F15.8)') 
     1                      ' energy decreased too much from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         
         FIXIMAGE=.TRUE.
         GOTO 20
      ELSE
C
C  Energy increased - try again with a smaller step size
C
         IF (NDECREASE.GT.10) THEN ! DJW
            NFAIL=NFAIL+1
            WRITE(LOG_FH,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL
            !
! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
               XCOORDS(1:N)=XSAVE(1:N)
               GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point
!              DO J1=1,N
!                 GRAD(J1)=GNEW(J1) ! GRAD contains the gradient at the lowest energy point
!                 XCOORDS(J1)=XCOORDS(J1)-STP*W(ISPT+POINT*N+J1)
!              ENDDO
            ITER=0   !  try resetting
             IF (NFAIL.GT.5) THEN         
               WRITE(LOG_FH,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.
!              STOP
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
               RETURN
            ENDIF
            GOTO 30
         ENDIF
      !
! Resetting to XSAVE and adding 0.1 of the step should be the same as subtracting 
! 0.9 of the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
            XCOORDS(1:N)=XSAVE(1:N)
            DO J1=1,N
               XCOORDS(J1)=XCOORDS(J1)+0.1D0*STP*W(ISPT+POINT*N+J1)
            ENDDO 
         ENDIF
         STP=STP/1.0D1
         NDECREASE=NDECREASE+1
         IF (DEBUG) WRITE(LOG_FH,'(A,F20.10,A,F20.10,A,F20.10)') 
     1                      ' energy increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         FIXIMAGE=.TRUE.
         GOTO 20
      ENDIF
C
C     Compute the new step and gradient change
C
30    NPT=POINT*N

         DO J1=1,N
            W(ISPT+NPT+J1)= STP*W(ISPT+NPT+J1) ! save the step taken
            W(IYPT+NPT+J1)= GRAD(J1)-W(J1)     ! save gradient difference: W(1:N) contains the old gradient
         ENDDO
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      FIXIMAGE=.FALSE.
            IF (CENT) CALL CENTRE2(XCOORDS)
      SMINKCURRENTP=SMINKCURRENT
      GOTO 10

      RETURN
      END
