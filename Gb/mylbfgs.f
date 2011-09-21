
      SUBROUTINE MYLBFGS(N,M,XCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET)
      ! declarations {{{
      USE COMMONS
      USE PORFUNCS

      IMPLICIT NONE

      ! subroutine {{{
      INTEGER N,M,ITMAX,ITDONE,NP,J2,J3,NFAIL,NDECREASE,NGUESS,NDUMMY
      DOUBLE PRECISION, DIMENSION(3*NATOMS) :: XCOORDS   
      LOGICAL RESET
      ! }}}
      ! local {{{
      INTEGER :: J1
      DOUBLE PRECISION :: GRAD(3*NATOMS),SLENGTH,DDOT,EPLUS,EMINUS,DIFF,DUMMY,WTEMP(3*NATOMS)
      DOUBLE PRECISION TMPANG(3*NATOMS), TMPCOORDS(3*NATOMS)
      DOUBLE PRECISION EPS,ENERGY,ENEW,GNEW(3*NATOMS),OVERLAP,OLDX(3*NATOMS),OLDOLDX(3*NATOMS),VGUESS(3),
     1                 X1, Y1, Z1, X2, Y2, Z2, TRY(3*NATOMS), D1, D2, RBCOORDS(18), DUMMY2, DIST, DIST1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W

      LOGICAL DIAGCO, YESNO NOTCALLED, CTEST, MFLAG

      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,POTEL,QSTART,QFINISH
      DOUBLE PRECISION OLDCART(3*NATOMS), DELTAQ(N),DELTACART(3*NATOMS),LEPSILON,DOT1,DOT2
      DOUBLE PRECISION LCART(3*NATOMS),OLDQ(N),NEWQ(N),OLDGINT(N),GINT(N),XINT(N),XSAVE(N)
      DOUBLE PRECISION, ALLOCATABLE :: FRAMES(:,:), PE(:), MODGRAD(:)
      LOGICAL NOCOOR, FAILED, COREDONE
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,INMC,IYCN,ISCN
      INTEGER KD, NNZ
      LOGICAL EVAP, GUIDECHANGET, GUIDET, EVAPREJECT, CSMDOGUIDET
      ! }}}
      ! common, save {{{
      COMMON /MYPOT/ POTEL
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /EV/ EVAP, evapreject

      SAVE W, DIAG, ITER, POINT, ISPT, IYPT, NPT
      ! }}}
      ! }}}
      ! subroutine body {{{

      ! intro {{{

      IF (.NOT.ALLOCATED(DIAG)) ALLOCATE(DIAG(N))       ! SAVE doesn't work otherwise for Sun
      IF (.NOT.ALLOCATED(W)) ALLOCATE(W(N*(2*M+1)+2*M)) ! SAVE doesn't work otherwise for Sun
      IF (SIZE(W,1).NE.N*(2*M+1)+2*M) THEN ! mustn't call mylbfgs with changing number of variables!!!
         WRITE(LFH, '(A,I10,A,I10,A)') 'ERROR, dimension of W=',SIZE(W,1),' but N*(2*M+1)+2*M=',N*(2*M+1)+2*M,' in mylbfgs'
         call exit(10)
      ENDIF

      COREDONE=.FALSE.

      NFAIL=0
      IF (GUIDECHANGET) ITER=0
      IF (RESET) ITER=0
      ITDONE=0
      FIXIMAGE=.FALSE.

      IF (DEBUG) THEN
         IF (RESET.OR.GUIDECHANGET) WRITE(LFH,'(A)') 'mylbfgs> Resetting LBFGS minimiser'
         IF (.NOT.(RESET.OR.GUIDECHANGET)) WRITE(LFH,'(A)') 'mylbfgs> Not resetting LBFGS minimiser'
      ENDIF

      IF (DUMPT) THEN
            WRITE(DUMPXYZUNIT,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT,11) NQ
11          FORMAT(1X,'QUENCH NUMBER ',I6,' initial points in mylbfgs')
            WRITE(DUMPXYZUNIT,'(A2,3F20.10)') ('LA ',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=1,NATOMS-NS)
            IF (NS.GT.0) WRITE(DUMPXYZUNIT,'(A2,3F20.10)') 
     1          ('LB',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
      ENDIF
      CALL POTENTIAL(XCOORDS,GRAD,ENERGY,.TRUE.,.FALSE.)

      IF (EVAPREJECT) RETURN
      POTEL=ENERGY

      IF (DEBUG) WRITE(LFH,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'
! }}}
C  10 Termination test.  {{{

10    CALL FLUSH(LFH)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
            MFLAG=.TRUE.
         IF (EVAP) MFLAG=.FALSE. ! do not allow convergence if we happen to have a small RMS and EVAP is true'
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            IF (DEBUG) WRITE(LFH,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'

            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
         IF (DEBUG) FIXIMAGE=.FALSE.
         IF (DEBUG) WRITE(LFH,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
         RETURN
      ENDIF
!}}}

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
            DO J1=1,N
               DUMMY=-GRAD(J1)*DIAG(J1)
               W(ISPT+J1)=DUMMY
               W(J1)=DUMMY
            ENDDO
            GNORM=DSQRT(DDOT(N,GRAD,1,GRAD,1))

C  Make the first guess for the step length cautious.

         STP=MIN(1.0D0/GNORM,GNORM)
         ! }}}
      ELSE 
         ! {{{
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
C
C  Update estimate of diagonal inverse Hessian elements
C
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            IF (YY.EQ.0.0D0) THEN
               WRITE(LFH,'(A)') 'WARNING, resetting YY to one in mylbfgs'
               YY=1.0D0
            ENDIF
            IF (YS.EQ.0.0D0) THEN
               WRITE(LFH,'(A)') 'WARNING, resetting YS to one in mylbfgs'
               YS=1.0D0
            ENDIF
            DO J1=1,N
               DIAG(J1)= YS/YY
            ENDDO
         ELSE
            WRITE(LFH,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(LFH,235) J1
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
         ! }}}
      ENDIF

C  Store the new search direction

      IF (ITER.GT.0) THEN
         DO J1=1,N
            W(ISPT+POINT*N+J1)= W(J1)
         ENDDO
      ENDIF

      DOT1=SQRT(DDOT(N,GRAD,1,GRAD,1))

C  Overflow has occasionally occurred here.
C  We only need the sign of the overlap, so use a temporary array with
C  reduced elements.

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
C        IF (DEBUG) PRINT*,'Search direction has positive projection onto gradient - resetting'
C        ITER=0
C        GOTO 10
         IF (DEBUG) WRITE(LFH,'(A)') 'Search direction has positive projection onto gradient - reversing step'
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
!
! Save XCOORDS here so that we can undo the step reliably including the
! non-linear projection for Thomson for the angular coordinates.
!
         XSAVE(1:N)=XCOORDS(1:N) 
         DO J1=1,N
            XCOORDS(J1)=XCOORDS(J1)+STP*W(ISPT+POINT*N+J1)
         ENDDO 
!

      NDECREASE=0
      LEPSILON=1.0D-6

! csw34> INCREMENT THE FORCE CONSTANT FOR STEERED MINIMISATION

      CALL POTENTIAL(XCOORDS,GNEW,ENEW,.TRUE.,.FALSE.)

      IF (EVAPREJECT) RETURN
C
C  Catch cold fusion for ionic potentials and discard.
C
C  Changed EREAL for cold fusion to 1.0D6 rather than 0.0D0, which could result in steps being accepted
C  for systems with positive energies. - khs26 26/11/09

C  We need to transform the newly obtained Cartesian gradient for CHARMM and internals.
C  NOCOOR is true because we dont need to transform the coordinates.
C

      IF (((ENEW-ENERGY.LE.MAXERISE).OR.EVAP.OR.GUIDECHANGET).AND.(ENEW-ENERGY.GT.MAXEFALL)) THEN
        ! {{{
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         DO J1=1,3*NATOMS
            GRAD(J1)=GNEW(J1)
         ENDDO
         IF (DEBUG) WRITE(LFH,'(A,F20.10,G20.10,A,I6,A,F13.10)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,
     1           ' LBFGS steps, step:',STP*SLENGTH
C
C  Step finished so can reset OLDQ to new XINT, OLDCART to new LCART,
C  as well as the Cartesian and internal gradients.
C
C
C  Try to take an extra step using the two previous geometries.
C 
C          GOTO 112
C          IF (MOD(ITDONE,3).EQ.0) THEN
C             NGUESS=0
C 111         CONTINUE
C             DO J1=1,NATOMS
C                X1=OLDX(3*(J1-1)+1)-OLDOLDX(3*(J1-1)+1)
C                Y1=OLDX(3*(J1-1)+2)-OLDOLDX(3*(J1-1)+2)
C                Z1=OLDX(3*(J1-1)+3)-OLDOLDX(3*(J1-1)+3)
C                X2=XCOORDS(3*(J1-1)+1)-OLDX(3*(J1-1)+1)
C                Y2=XCOORDS(3*(J1-1)+2)-OLDX(3*(J1-1)+2)
C                Z2=XCOORDS(3*(J1-1)+3)-OLDX(3*(J1-1)+3)
C                VGUESS(1)=(x2*(x1*x2 + y1*y2 + z1*z2))/(Sqrt(x1**2 + y1**2 + z1**2)*Sqrt(x2**2 + y2**2 + z2**2)) + 
C      -  ((x2*(y1*y2 + z1*z2) - x1*(y2**2 + z2**2))*
C      -     Sqrt(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))))/
C      -   Sqrt((x2*y1 - x1*y2)**2 + (x2*z1 - x1*z2)**2 + (y2*z1 - y1*z2)**2)
C                VGUESS(2)=(y2*(x1*x2 + y1*y2 + z1*z2))/(Sqrt(x1**2 + y1**2 + z1**2)*Sqrt(x2**2 + y2**2 + z2**2)) + 
C      -  ((-(x2**2*y1) + x1*x2*y2 + z2*(y2*z1 - y1*z2))*
C      -     Sqrt(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))))/
C      -   Sqrt((x2*y1 - x1*y2)**2 + (x2*z1 - x1*z2)**2 + (y2*z1 - y1*z2)**2)
C                VGUESS(3)=(z2*(x1*x2 + y1*y2 + z1*z2))/(Sqrt(x1**2 + y1**2 + z1**2)*Sqrt(x2**2 + y2**2 + z2**2)) + 
C      -  ((-(x2**2*z1) + x1*x2*z2 + y2*(-(y2*z1) + y1*z2))*
C      -     Sqrt(1 - (x1*x2 + y1*y2 + z1*z2)**2/((x1**2 + y1**2 + z1**2)*(x2**2 + y2**2 + z2**2))))/
C      -   Sqrt((x2*y1 - x1*y2)**2 + (x2*z1 - x1*z2)**2 + (y2*z1 - y1*z2)**2)
C                D1=SQRT(VGUESS(1)**2+VGUESS(2)**2+VGUESS(3)**2)
C                IF (D1.LT.0.1) THEN
C                   TRY(3*(J1-1)+1)=XCOORDS(3*(J1-1)+1)+VGUESS(1)*1.0D0
C                   TRY(3*(J1-1)+2)=XCOORDS(3*(J1-1)+2)+VGUESS(2)*1.0D0
C                   TRY(3*(J1-1)+3)=XCOORDS(3*(J1-1)+3)+VGUESS(3)*1.0D0
C                ENDIF
C             ENDDO
C             CALL POTENTIAL(TRY,GNEW,EGUESS,.FALSE.,.FALSE.)
C             WRITE(*,'(A,3G20.10)') 'ENEW,EGUESS,change=',ENEW,EGUESS,EGUESS-ENEW
C             IF (EGUESS-ENEW.LT.0.0D0) THEN
C                NGUESS=NGUESS+1
C                ENEW=EGUESS
C                DO J1=1,N
C                   OLDOLDX(J1)=OLDX(J1)
C                   OLDX(J1)=XCOORDS(J1)
C                   XCOORDS(J1)=TRY(J1)
C                ENDDO
C                IF (NGUESS.LT.6) GOTO 111
C             ENDIF
C          ENDIF
C          DO J1=1,N
C             OLDOLDX(J1)=OLDX(J1)
C             OLDX(J1)=XCOORDS(J1)
C          ENDDO
C 
C 112      CONTINUE
C
C  May want to prevent the PE from falling too much if we are trying to visit all the
C  PE bins. Halve the step size until the energy decrease is in range.
C
         ! }}}
      ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
      ! {{{
C
C  Energy decreased too much - try again with a smaller step size
C
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            WRITE(LFH,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find an energy in the required range, NFAIL=',NFAIL
!
! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
            XCOORDS(1:N)=XSAVE(1:N)
            GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point

            ITER=0   !  try resetting
            IF (NFAIL.GT.20) THEN
               WRITE(LFH,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.
!              STOP
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
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
!           DO J1=1,N
!              XCOORDS(J1)=XCOORDS(J1)-0.5*STP*W(ISPT+POINT*N+J1)
!           ENDDO 
            STP=STP/2.0D0
         NDECREASE=NDECREASE+1
         IF (DEBUG) WRITE(LFH,'(A,F19.10,A,F16.10,A,F15.8)') 
     1                      ' energy decreased too much from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         
         FIXIMAGE=.TRUE.
         ! }}}
      ELSE
        ! {{{

C  Energy increased - try again with a smaller step size

         IF (NDECREASE.GT.10) THEN ! DJW
            NFAIL=NFAIL+1
            WRITE(LFH,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL
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
!            IF (NFAIL.GT.20) THEN
! bs360: smaller NFAIL 
             IF (NFAIL.GT.5) THEN         
               WRITE(LFH,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.

! jwrm2> For testing: dumps coords and exits when failing
!        Allows use of CHECKD on failed coordinates
!               OPEN (UNIT = 26, FILE = 'failcoords')
!               DO J1=1,NATOMS
!                 WRITE (26, '(3F20.10)') XCOORDS(J1), XCOORDS(J1+1), XCOORDS(J1+2) 
!               ENDDO
!               CLOSE(26)
!               STOP
! jwrm2> End coordinate dumping

!              STOP
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
               RETURN
            ENDIF
            GOTO 30
         ENDIF

            XCOORDS(1:N)=XSAVE(1:N)
            DO J1=1,N
               XCOORDS(J1)=XCOORDS(J1)+0.1D0*STP*W(ISPT+POINT*N+J1)
            ENDDO 

         STP=STP/1.0D1
         NDECREASE=NDECREASE+1
         IF (DEBUG) WRITE(LFH,'(A,F20.10,A,F20.10,A,F20.10)') 
     1                      ' energy increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         FIXIMAGE=.TRUE.
         ! }}}
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
      IF (DUMPT.AND.DEBUG) THEN
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,'(A,I8,A,G20.10)') 'at step ',ITER,' energy=',ENERGY
            WRITE(DUMPXYZUNIT+NP,'(A2,3F20.10)') ('LA ',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=1,NATOMS-NS)
            IF (NS.GT.0) WRITE(DUMPXYZUNIT+NP,'(A2,3F20.10)') 
     1          ('LB',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
      ENDIF
      IF (CENT) CALL CENTRE2(XCOORDS)
      GOTO 10

      RETURN
      ! }}}
      END
