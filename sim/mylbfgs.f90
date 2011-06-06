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
!> @param[in] M   
!>             is an INTEGER variable that must be set by the user to
!>             the number of corrections used in the BFGS update. It
!>             is not altered by the routine. Values of M less than 3 are
!>             not recommended; large values of M will result in excessive
!>             computing time. 3<= M <=7 is recommended. Restriction: M>0.
!
!> @param[out] MFLAG 
!>             returns .TRUE. is the algorithm has converged; .FALSE. otherwise
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
!> @param[inout] R is a vector of coordinates. On initial entry
!>             it must be set by the user to the values of the initial
!>             estimate of the solution vector. On exit with IFLAG=0, it
!>             contains the values of the variables at the best point
!>             found (usually a solution).
!
!
! }}}
! ------------------------------------------------------------------
      SUBROUTINE MYLBFGS(X,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET)
      ! ====================================
     ! declarations {{{ 
      USE V
      USE PORFUNCS
      
      IMPLICIT NONE
! subroutine parameters  {{{
      DOUBLE PRECISION, INTENT(INOUT) :: X(:)
      LOGICAL,INTENT(IN) :: DIAGCO
      DOUBLE PRECISION,INTENT(IN) :: EPS
      LOGICAL,INTENT(OUT) :: MFLAG
      DOUBLE PRECISION,INTENT(OUT) :: ENERGY
      INTEGER,INTENT(IN) :: ITMAX
      INTEGER,INTENT(OUT) :: ITDONE
      LOGICAL,INTENT(IN) :: RESET
      ! }}}
       ! local parameters {{{
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W, WTEMP, GRADX, GNEW, XSAVE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ALPHA, RHO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: WSS,WDG
      INTEGER ITER, NFAIL, CP, NDECREASE
      DOUBLE PRECISION :: GNORM,SLENGTH,E,ENEW,YS,YY,SQ,YR,BETA,POTEL, DDOT, STP
      DOUBLE PRECISION ::       DOT1,DOT2,OVERLAP, DUMMY

      INTEGER :: N,M,K,IX,IXMIN,IXMAX
      COMMON /MYPOT/ QE

      SAVE W,DIAG,ITER
      ! 
      M=M_LBFGS
      !
      !   }}}
      ! ====================================
      ! Labels:  {{{
      !         10  - termination test
      !         30  - compute the new step and gradient change
      ! }}}
      ! ====================================
      ! subroutine body  {{{

      ! initializations {{{

      ALLOCATE(RHO(M),ALPHA(M),W(N),WSS(N,M),WDG(N,M))

      NFAIL=0 ; ITDONE=0 ; IF (RESET) ITER=0 

      IF (DEBUG) THEN
         IF (RESET) WRITE(LFH,'(A)')            'mylbfgs> Resetting LBFGS minimiser'
         IF (.NOT.RESET) WRITE(LFH,'(A)')       'mylbfgs> Not resetting LBFGS minimiser'
      ENDIF

      !        evaluate gradient (.TRUE.) but not Hessian (.FALSE.)
      CALL POTENTIAL(X,GRADX,E,RMS,.TRUE.,.FALSE.)

      IF (DEBUG) WRITE(LFH,101) ' Energy and RMS force=', E,RMS, &
                                ' after ',ITDONE,' LBFGS steps'
101   FORMAT(A,F20.10,G20.10,A,I6,A)
102   FORMAT(A,I10,A,I10,A)

! }}}

!  10 Termination test {{{

10    CALL FLUSH(LFH)

      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
        MFLAG=.TRUE.
         IF (MFLAG) THEN
            IF (DEBUG) WRITE(LFH,101) ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'
            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
         IF (DEBUG) WRITE(LFH,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
         RETURN
      ENDIF
      ! }}}

      IF (ITER.EQ.0) THEN
         ! {{{
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(LFH,240)
 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF

         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            WRITE(LFH,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(LFH,235) J1
 235              FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,&
                ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
                  STOP
               ENDIF
            ENDDO
         ELSE
            DIAG=DGUESS
         ENDIF
!
!     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
!     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
!

! W     =>      W 
! RHO   =>      N,...,N+M      
! ALPHA =>      N+M+1,...,N+2M
! WSS   =>      storage of the last M steps 
! WDG   =>      storage of the last M gradient differences 

!
!  NR step for diagonal inverse Hessian
!
     DO IX=1,N
            DUMMY=-GRAD(IX)*DIAG(IX)
            WSS(IX,1)=DUMMY
            W(IX)=DUMMY
     ENDDO
     GNORM=DSQRT(DDOT(N,GRAD,1,GRAD,1))

!
!  Make the first guess for the step length cautious.
!
         STP=MIN(1.0D0/GNORM,GNORM)
         ! }}}
      ELSE 
         ! {{{
         BOUND=MIN(ITER,M) 
         YS= DDOT(N,WDG(1,1+POINT),1,WSS(1,1+POINT),1)
!
!  Update estimate of diagonal inverse Hessian elements {{{
!
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,WDG(1,1+POINT),1,WDG(1,1+POINT),1)
            IF (YY.EQ.0.0D0) THEN
               WRITE(LFH,'(A)') 'WARNING, resetting YY to one in mylbfgs'
               YY=1.0D0
            ENDIF
            IF (YS.EQ.0.0D0) THEN
               WRITE(LFH,'(A)') 'WARNING, resetting YS to one in mylbfgs'
               YS=1.0D0
            ENDIF
               DIAG= YS/YY
         ELSE
            WRITE(LFH,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(LFH,235) J1
                  STOP
               ENDIF
            ENDDO
         ENDIF
! }}}
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!
         CP= POINT
         IF (POINT.EQ.0) CP=M
               RHO(CP)= 1.0D0/YS
               W= -GRAD
               CP= POINT
         DO J1= 1,BOUND
            CP=CP-1
            IF (CP.EQ.-1) CP=M-1
            SQ=DDOT(N,WSS(1,CP+1),1,W,1)
            ALPHA(CP+1)=RHO(CP+1)*SQ
            CALL DAXPY(N,-ALPHA(CP+1),WDG(1,CP+1),1,W,1)
         ENDDO
        
         W=DIAG*W

         DO J1=1,BOUND
            YR=DDOT(N,W(1,CP+1),1,W,1)
            BETA=RHO(CP+1)*YR
            BETA=ALPHA(CP+1)-BETA
            CALL DAXPY(N,BETA,WSS(1,CP+1),1,W,1)
            CP=CP+1
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0  
         ! }}}
      ENDIF
!
!  Store the new search direction (160)
!
      IF (ITER.GT.0) THEN
            WSS(:,1+POINT)=W
      ENDIF

!  test for overlap {{{
!
!  Overflow has occasionally occurred here.
!  We only need the sign of the overlap, so use a temporary array with
!  reduced elements.
!
      DUMMY=1.0D0
      DO IX=1,N
         IF (ABS(W(IX)).GT.DUMMY) DUMMY=ABS(W(IX))
      ENDDO
      
      WTEMP=W/DUMMY

      DOT1=SQRT(DDOT(N,GRAD,1,GRAD,1))
      DOT2=SQRT(DDOT(N,WTEMP,1,WTEMP,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) THEN
                   OVERLAP=DDOT(N,GRAD,1,WTEMP,1)/(DOT1*DOT2)
      ENDIF
      IF (OVERLAP.GT.0.0D0) THEN
         IF (DEBUG) WRITE(LFH,'(A)') 'Search direction has positive projection onto gradient - reversing step'
         WSS(1:N,1+POINT)=-W  !!! DJW, reverses step
      ENDIF
! }}}

      W=GRAD
      SLENGTH=SQRT(SUM(WSS(:,1+POINT)**2))
      IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH
!
!  We now have the proposed step.
!
! Save X here so that we can undo the step reliably including the
! non-linear projection for Thomson for the angular coordinates.
!
         XSAVE=X
         X=X+STP*WSS(:,1+POINT)
      CALL POTENTIAL(X,GNEW,ENEW,RMS,.TRUE.,.FALSE.)

      IF ((ENEW-ENERGY.LE.MAXERISE).AND.(ENEW-ENERGY.GT.MAXEFALL)) THEN
        ! {{{;
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         GRAD=GNEW

      IF (DEBUG) WRITE(LFH,103) ' Energy and RMS force=', ENERGY,RMS, &
                                ' after ',ITDONE,' LBFGS steps, ', &
                                ' step:',STP*SLENGTH

103   FORMAT(A,F20.10,G20.10,A,I6,A,A,F13.10)

!  May want to prevent the PE from falling too much if we are trying to visit all the
!  PE bins. Halve the step size until the energy decrease is in range.
!
! }}}
      ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
      ! {{{
!
!  Energy decreased too much - try again with a smaller step size
!
         IF (NDECREASE.GT.5) THEN
           ! {{{
            NFAIL=NFAIL+1
            WRITE(LFH,'(A,A,G20.10)')   ' in mylbfgs LBFGS step cannot find an energy in the required range, ',&
                                        '    NFAIL=',NFAIL
      !
! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
               X=XSAVE
               GRAD=GNEW ! GRAD contains the gradient at the lowest energy point

            ITER=0   !  try resetting
            IF (NFAIL.GT.20) THEN
               WRITE(LFH,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.
               RETURN
            ENDIF
            GOTO 30
            ! }}}
         ENDIF
!
! Resetting to XSAVE and adding half the step should be the same as subtracting 
! half the step. 
!
         X=XSAVE
         X=X+0.5*STP*WSS(:,1+POINT)
         STP=STP/2.0D0
         NDECREASE=NDECREASE+1
         IF (DEBUG) WRITE(LFH,105) &
                ' energy decreased too much from ',ENERGY,&
                ' to ',ENEW,&
                ' decreasing step to ', STP*SLENGTH
        
105 FORMAT(A,F19.10,A,F16.10,A,F15.8) 

         FIXIMAGE=.TRUE.
         GOTO 20
         ! }}}
      ELSE
        ! {{{
!
!  Energy increased - try again with a smaller step size
!
         IF (NDECREASE.GT.10) THEN ! DJW
           ! {{{
            NFAIL=NFAIL+1
            WRITE(LFH,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL
            !
! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
            X=XSAVE
            GRAD=GNEW ! GRAD contains the gradient at the lowest energy point
            ITER=0   !  try resetting
             IF (NFAIL.GT.5) THEN         
               WRITE(LFH,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.
               RETURN
            ENDIF
            GOTO 30
            ! }}}
         ENDIF
      !
! Resetting to XSAVE and adding 0.1 of the step should be the same as subtracting 
! 0.9 of the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
            X=XSAVE
            X=X+0.1D0*STP*SUM(WSS(:,1+POINT))
            ! }}}
      ENDIF
         STP=STP/1.0D1
         NDECREASE=NDECREASE+1
         IF (DEBUG) WRITE(LFH,104) ' energy increased from ',ENERGY,&
                                  ' to ',ENEW,&
                                  ' decreasing step to ',STP*SLENGTH
104 FORMAT(A,F20.10,A,F20.10,A,F20.10) 
         FIXIMAGE=.TRUE.

         GOTO 20
      ENDIF
!
!     Compute the new step and gradient change
!
30    NPT=POINT*N

      WSS(:,1+POINT)=STP*WSS(:,1+POINT)     ! save the step taken
      WDG(:,1+POINT)=GRAD-W                 ! save gradient difference: W contains the old gradient
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      FIXIMAGE=.FALSE.
      IF (CENT) CALL CENTRE2(R)

      GOTO 10           ! Go to the termination test
! }}}
      RETURN
      END
