      SUBROUTINE MYLBFGS(N,M,XCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET,NP)
      USE COMMONS
      USE MODAMBER
      USE MODAMBER9, ONLY : STEEREDMINT, LOCALSTEEREDMINT, SMINK, SMINKINC, SMINKCURRENT
      USE MODCHARMM
      USE PORFUNCS
      IMPLICIT NONE
      INTEGER N,M,J1,ITMAX,ITDONE,NP,J2,J3,NFAIL,NDECREASE,NGUESS,NDUMMY
      DOUBLE PRECISION XCOORDS(3*NATOMS),GRAD(3*NATOMS),SLENGTH,DDOT,EPLUS,EMINUS,DIFF,DUMMY,WTEMP(3*NATOMS)
      DOUBLE PRECISION TMPANG(3*NATOMS), TMPCOORDS(3*NATOMS)
      DOUBLE PRECISION EPS,ENERGY,ENEW,GNEW(3*NATOMS),OVERLAP,OLDX(3*NATOMS),OLDOLDX(3*NATOMS),VGUESS(3),
     1                 X1, Y1, Z1, X2, Y2, Z2, TRY(3*NATOMS), D1, D2, RBCOORDS(18), DUMMY2, DIST, DIST1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W
      LOGICAL DIAGCO, YESNO, RESET, NOTCALLED, CTEST, MFLAG
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,POTEL,QSTART,QFINISH
      DOUBLE PRECISION OLDCART(3*NATOMS), DELTAQ(N),DELTACART(3*NATOMS),LEPSILON,DOT1,DOT2
      DOUBLE PRECISION LCART(3*NATOMS),OLDQ(N),NEWQ(N),OLDGINT(N),GINT(N),XINT(N),XSAVE(N),SMINKCURRENTP
      DOUBLE PRECISION, ALLOCATABLE :: FRAMES(:,:), PE(:), MODGRAD(:)
      LOGICAL NOCOOR, FAILED, COREDONE
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,INMC,IYCN,ISCN
      INTEGER KD, NNZ
      LOGICAL EVAP, GUIDECHANGET, GUIDET, EVAPREJECT, SMINKCHANGET, CSMDOGUIDET
      SAVE W, DIAG, ITER, POINT, ISPT, IYPT, NPT

      IF (.NOT.ALLOCATED(DIAG)) ALLOCATE(DIAG(N))       ! SAVE doesn't work otherwise for Sun
      IF (.NOT.ALLOCATED(W)) ALLOCATE(W(N*(2*M+1)+2*M)) ! SAVE doesn't work otherwise for Sun
      IF (SIZE(W,1).NE.N*(2*M+1)+2*M) THEN ! mustn't call mylbfgs with changing number of variables!!!
         WRITE(MYUNIT, '(A,I10,A,I10,A)') 'ERROR, dimension of W=',SIZE(W,1),' but N*(2*M+1)+2*M=',N*(2*M+1)+2*M,' in mylbfgs'
      ENDIF
      SMINKCHANGET=.FALSE.
      SMINKCURRENT=0.0D0
      SMINKCURRENTP=0.0D0
      LOCALSTEEREDMINT=.FALSE.
      IF (STEEREDMINT) LOCALSTEEREDMINT=.TRUE.

      IF(CHRMMT) CALL UPDATENBONDS(XCOORDS)

      NFAIL=0
      IF (GUIDECHANGET) ITER=0
      IF (RESET) ITER=0
      ITDONE=0
      FIXIMAGE=.FALSE.
      IF (DEBUG) THEN
         IF (RESET.OR.GUIDECHANGET) WRITE(MYUNIT,'(A)') 'mylbfgs> Resetting LBFGS minimiser'
         IF (.NOT.(RESET.OR.GUIDECHANGET)) WRITE(MYUNIT,'(A)') 'mylbfgs> Not resetting LBFGS minimiser'
      ENDIF

      IF (Q4T) CALL ORDERQ4(NATOMS,XCOORDS,QSTART)

      IF (DUMPT) THEN
         IF (ARNO) THEN
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS+2
            WRITE(DUMPXYZUNIT+NP,11) NP,NQ(NP)
            WRITE(DUMPXYZUNIT+NP,'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(DUMPXYZUNIT+NP,'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            WRITE(DUMPXYZUNIT+NP,65) (XCOORDS(J1),J1=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
         ELSE IF (TIP) THEN
            WRITE(DUMPXYZUNIT+NP,'(I6)') (NATOMS/2)*3
            WRITE(DUMPXYZUNIT+NP,'(A,I5)') 'LBFGS iteration ',ITER
            DO J2=1,NATOMS/2
     1              XCOORDS(3*(NATOMS/2+J2-1)+1),XCOORDS(3*(NATOMS/2+J2-1)+2),XCOORDS(3*(NATOMS/2+J2-1)+3),RBCOORDS)
               WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
               WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
               WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            ENDDO
         ELSE IF (AMBER) THEN
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,11) NP,NQ(NP)
            DO J2=1,NATOMS
               WRITE(DUMPXYZUNIT+NP,'(A,3F20.10)') TYPECH(J2)(1:1),(XCOORDS(3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSEIF (NCORE(NP).GT.0) THEN
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,11) NQ(NP)
            WRITE(DUMPXYZUNIT+NP,'(A2,3F20.10)') ('LB',XCOORDS(3*(I-1)+1),XCOORDS(3*(I-1)+2),XCOORDS(3*(I-1)+3),
     &                                               I=1,NATOMS-NCORE(NP))
            IF (NCORE(NP).GT.0) WRITE(DUMPXYZUNIT+NP,'(A2,3F20.10)') 
     &            ('LA ',XCOORDS(3*(I-1)+1),XCOORDS(3*(I-1)+2),XCOORDS(3*(I-1)+3),I=NATOMS-NCORE(NP)+1,NATOMS)
         ELSE
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,11) NQ(NP)
11          FORMAT(1X,'QUENCH NUMBER ',I6,' initial points in mylbfgs')
            WRITE(DUMPXYZUNIT+NP,'(A2,3F20.10)') ('LA ',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=1,NATOMS-NS)
            IF (NS.GT.0) WRITE(DUMPXYZUNIT+NP,'(A2,3F20.10)') 
     1          ('LB',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
         ENDIF
      ENDIF
      IF ((TOSI.OR.WELCH.OR.RGCL2.OR.AMBER.OR.ARNO.OR.PACHECO.OR.TIP.OR.CHRMMT.OR.AMBERT 
     &   .OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT.OR.JMT) 
     &   .AND.(ENERGY.LT.COLDFUSIONLIMIT)) THEN
         WRITE(MYUNIT,'(A,G20.10)') 'ENERGY=',ENERGY
         ENERGY=1.0D6
         POTEL=1.0D6
         RMS=1.0D0
         WRITE(MYUNIT,'(A)') ' Cold fusion diagnosed - step discarded'
         RETURN
      ENDIF


      IF (CHRMMT .AND. GCHARMMFAIL) THEN
          WRITE(MYUNIT,'(A)') 'Failure in CHARMM energy/gradient evaluation - geometry discarded.'
          RETURN
      ENDIF

      IF (INTMINT) THEN
         OLDCART(1:3*NATOMS)=XCOORDS(1:3*NATOMS) ! store cartesians in OLDCART for both CHARMM and UNRES
            NOCOOR=.FALSE. ! calculate internals therefore NOCOOR is false
            GINT(1:N)=0.0D0 ! to prevent NaN's for Sun!
            XINT(1:N)=0.0D0 ! to prevent NaN's for Sun!
            OLDQ(1:N)=XINT(1:N)    ! store internals
            OLDGINT(1:N)=GINT(1:N) ! store gradient in internals
      ENDIF
      IF (EVAPREJECT) RETURN
      POTEL=ENERGY

      IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'

      IF ((DBPT.OR.DBPTDT.OR.MSTBINT.OR.MSSTOCKT.OR.MULTPAHAT.OR.NPAHT.OR.PAHW99T.OR.PYGT.OR.TDHDT.OR.SILANET) 
     &   .AND.(ENERGY.LT.-5.0D4)) THEN
         WRITE(MYUNIT,'(A,G20.10)') 'ENERGY=',ENERGY
         ENERGY=0.0D0
         POTEL=0.0D0
         RMS=1.0D0
         WRITE(MYUNIT,'(A)') ' Cold fusion diagnosed - step discarded'
         RETURN
      ENDIF
10    CALL FLUSH(MYUNIT)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
         IF (CHRMMT.AND.ACESOLV) THEN
            NCHENCALLS=ACEUPSTEP-1
            IF (DEBUG) WRITE(*,'(A,2G20.10,A)') ' mylbfgs> Energy and RMS force=',ENERGY,RMS,' after ACE update'
            IF (RMS.LE.EPS) MFLAG=.TRUE.
         ELSE
            MFLAG=.TRUE.
         ENDIF
         IF (EVAP) MFLAG=.FALSE. ! do not allow convergence if we happen to have a small RMS and EVAP is true'
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'


            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
         IF (DEBUG) FIXIMAGE=.FALSE.
         IF (DEBUG) WRITE(MYUNIT,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
         RETURN
      ENDIF

      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(MYUNIT,240)
 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            WRITE(MYUNIT,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(MYUNIT,235) J1
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
         ISPT= N+2*M    ! index for storage of search steps
         IYPT= ISPT+N*M ! index for storage of gradient differences
         IF (CHRMMT.AND.INTMINT) THEN
            DO I=1,N
               W(ISPT+I)= -GINT(I)*DIAG(I)
               W(I)= -GINT(I)*DIAG(I)
            ENDDO
            GNORM= DSQRT(DDOT(N,GINT,1,GINT,1))
         ELSE
            DO J1=1,N
               DUMMY=-GRAD(J1)*DIAG(J1)
               W(ISPT+J1)=DUMMY
               W(J1)=DUMMY
            ENDDO
            GNORM=DSQRT(DDOT(N,GRAD,1,GRAD,1))
         ENDIF
         STP=MIN(1.0D0/GNORM,GNORM)
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            IF (YY.EQ.0.0D0) THEN
               WRITE(MYUNIT,'(A)') 'WARNING, resetting YY to one in mylbfgs'
               YY=1.0D0
            ENDIF
            IF (YS.EQ.0.0D0) THEN
               WRITE(MYUNIT,'(A)') 'WARNING, resetting YS to one in mylbfgs'
               YS=1.0D0
            ENDIF
            DO J1=1,N
               DIAG(J1)= YS/YY
            ENDDO
         ELSE
            WRITE(MYUNIT,'(A)') 'using estimate of the inverse diagonal elements'
            DO J1=1,N
               IF (DIAG(J1).LE.0.0D0) THEN
                  WRITE(MYUNIT,235) J1
                  STOP
               ENDIF
            ENDDO
         ENDIF
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS
         IF (CHRMMT.AND.INTMINT) THEN
            DO I=1,N
               W(I)= -GINT(I)
            ENDDO
         ELSE
            DO J1=1,N
               W(J1)= -GRAD(J1)
            ENDDO
         ENDIF
         DO J1= 1,BOUND
            IF (CP.EQ.-1) CP=M-1
            SQ=DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)=W(N+CP+1)*SQ
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
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0  
      ENDIF
      IF (ITER.GT.0) THEN
         DO J1=1,N
            W(ISPT+POINT*N+J1)= W(J1)
         ENDDO
      ENDIF

      IF (CHRMMT.AND.INTMINT) THEN
         DOT1=SQRT(DDOT(N,GINT,1,GINT,1))
      ELSE
         DOT1=SQRT(DDOT(N,GRAD,1,GRAD,1))
      ENDIF
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
         IF (CHRMMT.AND.INTMINT) THEN
            OVERLAP=DDOT(N,GINT,1,WTEMP,1)/(DOT1*DOT2)
         ELSE
            OVERLAP=DDOT(N,GRAD,1,WTEMP,1)/(DOT1*DOT2)
        ENDIF
      ENDIF
      IF (OVERLAP.GT.0.0D0) THEN
         IF (DEBUG) WRITE(MYUNIT,'(A)') 'Search direction has positive projection onto gradient - reversing step'
         DO J1=1,N
            W(ISPT+POINT*N+J1)= -W(J1)  !!! DJW, reverses step
         ENDDO
      ENDIF

      IF (CHRMMT.AND.INTMINT) THEN
         DO I=1,N
            W(I)=GINT(I)
         ENDDO
      ELSE
         DO J1=1,N
            W(J1)=GRAD(J1)
         ENDDO
      ENDIF
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
      IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH
      IF (CHRMMT.AND.INTMINT) THEN
         DO J1=1,N
            XINT(J1)=XINT(J1)+STP*W(ISPT+POINT*N+J1)
            DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)
         ENDDO
      ELSE
         XSAVE(1:N)=XCOORDS(1:N) 
         DO J1=1,N
            XCOORDS(J1)=XCOORDS(J1)+STP*W(ISPT+POINT*N+J1)
         ENDDO 
         IF (PROJIT) THEN
            IF (THOMSONT) THEN
               TMPCOORDS(1:N)=XCOORDS(1:N)
               XCOORDS(1:N)=TMPANG(1:N)
            ELSE
               TMPCOORDS(1:N)=XCOORDS(1:N)
               XCOORDS(1:N)=TMPCOORDS(1:N)
            ENDIF
         ENDIF
         IF (PROJIHT) THEN
            IF (THOMSONT) THEN
               TMPCOORDS(1:N)=XCOORDS(1:N)
               XCOORDS(1:N)=TMPANG(1:N)
            ELSE
               TMPCOORDS(1:N)=XCOORDS(1:N)
               XCOORDS(1:N)=TMPCOORDS(1:N)
            ENDIF
         ENDIF
      ENDIF
      NDECREASE=0
      LEPSILON=1.0D-6

20    IF (INTMINT) THEN
         IF (CHRMMT) THEN
            NEWQ(1:N)=OLDQ(1:N)
            LCART(1:3*NATOMS)=OLDCART(1:3*NATOMS)
            IF (FAILED) THEN
               MFLAG=.FALSE.
               RETURN
            ENDIF
            LCART(1:3*NATOMS)=OLDCART(1:3*NATOMS)+DELTACART(1:3*NATOMS)
            XCOORDS(1:3*NATOMS)=OLDCART(1:3*NATOMS)+DELTACART(1:3*NATOMS)
         ENDIF
      ENDIF

      SMINKCHANGET=.FALSE.
      IF (LOCALSTEEREDMINT) THEN
         SMINKCURRENT=MIN(SMINKCURRENT+SMINKINC,SMINK)
         IF (SMINKCURRENT.NE.SMINKCURRENTP) SMINKCHANGET=.TRUE.
        IF (DEBUG) WRITE(MYUNIT,'(A,2F20.10,L5)') 'SMINKCURRENT,SMINKCURRENTP,SMINKCHANGET=',SMINKCURRENT,SMINKCURRENTP,SMINKCHANGET
      ENDIF



      IF (EVAPREJECT) return
      IF (CHRMMT .AND. GCHARMMFAIL) THEN
          WRITE(MYUNIT,'(A)') 'Failure in CHARMM energy/gradient evaluation - step discarded.'
          RETURN
      ENDIF
      IF ((TOSI.OR.WELCH.OR.RGCL2.OR.AMBER.OR.ARNO.OR.PACHECO.OR.TIP.OR.CHRMMT.OR.AMBERT 
     &   .OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT.OR.JMT)
     &   .AND.(ENEW.LT.COLDFUSIONLIMIT)) THEN
         ENERGY=1.0D6
         ENEW=1.0D6
         POTEL=1.0D6
         RMS=1.0D0
         WRITE(MYUNIT,'(A)') ' Cold fusion diagnosed - step discarded'
         RETURN
      ENDIF
      IF ((DBPT.OR.DBPTDT.OR.MSTBINT.OR.MSSTOCKT.OR.MULTPAHAT.OR.NPAHT.OR.PAHW99T.OR.PYGT.OR.TDHDT) .AND.(ENEW.LT.-5.0D4)) THEN
         ENERGY=0.0D0
         ENEW=0.0D0
         POTEL=0.0D0
         RMS=1.0D0
         WRITE(MYUNIT,'(A)') ' Cold fusion diagnosed - step discarded'
         RETURN
      ENDIF


      IF (CHRMMT.AND.INTMINT) THEN
         NOCOOR=.TRUE.
      ENDIF



      IF (FIXDIHEFLAG) ENERGY=ENEW

      IF (((ENEW-ENERGY.LE.MAXERISE).OR.EVAP.OR.GUIDECHANGET.OR.SMINKCHANGET).AND.(ENEW-ENERGY.GT.MAXEFALL)) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         DO J1=1,3*NATOMS
            GRAD(J1)=GNEW(J1)
         ENDDO
         IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,G20.10,A,I6,A,F13.10)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,
     1           ' LBFGS steps, step:',STP*SLENGTH
         IF (CHRMMT.AND.INTMINT) THEN
            OLDGINT(1:N)=GINT(1:N)
            OLDCART(1:3*NATOMS)=LCART(1:3*NATOMS)
            OLDQ(1:N)=XINT(1:N)
         ENDIF
      ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find an energy in the required range, NFAIL=',NFAIL
            IF (CHRMMT.AND.INTMINT) THEN ! need to reset X, XINT, G, GINT to original values
               XINT(1:N)=XINT(1:N)-STP*W(ISPT+POINT*N+1:ISPT+POINT*N+N)
               GINT(1:N)=OLDGINT(1:N)
               GRAD(1:3*NATOMS)=GNEW(1:3*NATOMS) ! here OPTIM uses GLAST ! DJW
               XCOORDS(1:3*NATOMS)=OLDCART(1:3*NATOMS)
            ELSE
               XCOORDS(1:N)=XSAVE(1:N)
               GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point


            ENDIF
            ITER=0   !  try resetting
            IF (NFAIL.GT.20) THEN
               WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.
               RETURN
            ENDIF
            GOTO 30
         ENDIF
         IF (CHRMMT.AND.INTMINT) THEN
            DO J1=1,N
               XINT(J1)=XINT(J1)-0.5*STP*W(ISPT+POINT*N+J1)
               DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)*0.5D0
            ENDDO
         ELSE
            XCOORDS(1:N)=XSAVE(1:N)
            DO J1=1,N
               XCOORDS(J1)=XCOORDS(J1)+0.5*STP*W(ISPT+POINT*N+J1)
            ENDDO 
            IF (PROJIT) THEN
               IF (THOMSONT) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF
            IF (PROJIHT) THEN
               IF (THOMSONT) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF

         ENDIF
         STP=STP/2.0D0
         NDECREASE=NDECREASE+1
         IF (DEBUG) WRITE(MYUNIT,'(A,F19.10,A,F16.10,A,F15.8)') 
     1                      ' energy decreased too much from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         
         FIXIMAGE=.TRUE.
         GOTO 20
      ELSE
         IF (NDECREASE.GT.10) THEN ! DJW
            NFAIL=NFAIL+1
            WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL
            IF (CHRMMT.AND.INTMINT) THEN ! need to reset X, XINT, G, GINT to original values
               XINT(1:N)=XINT(1:N)-STP*W(ISPT+POINT*N+1:ISPT+POINT*N+N)
               GINT(1:N)=OLDGINT(1:N)
               GRAD(1:3*NATOMS)=GNEW(1:3*NATOMS) ! here OPTIM uses GLAST ! DJW
               XCOORDS(1:3*NATOMS)=OLDCART(1:3*NATOMS)
            ELSE
               XCOORDS(1:N)=XSAVE(1:N)
               GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point
            ENDIF
            ITER=0   !  try resetting
             IF (NFAIL.GT.5) THEN         
               WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.


               RETURN
            ENDIF
            GOTO 30
         ENDIF
         IF (CHRMMT.AND.INTMINT) THEN
            DO J1=1,N
               XINT(J1)=XINT(J1)-0.9*STP*W(ISPT+POINT*N+J1)
               DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)*0.1D0
            ENDDO
         ELSE
            XCOORDS(1:N)=XSAVE(1:N)
            DO J1=1,N
               XCOORDS(J1)=XCOORDS(J1)+0.1D0*STP*W(ISPT+POINT*N+J1)
            ENDDO 

            IF (PROJIT) THEN
               IF (THOMSONT) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF
            IF (PROJIHT) THEN
               IF (THOMSONT) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF
         ENDIF
         STP=STP/1.0D1
         NDECREASE=NDECREASE+1
         IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,A,F20.10,A,F20.10)') 
     1                      ' energy increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         FIXIMAGE=.TRUE.
         GOTO 20
      ENDIF
30    NPT=POINT*N

      IF (CHRMMT.AND.INTMINT) THEN
         DO I=1,N
            W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
            W(IYPT+NPT+I)= GINT(I)-W(I)
         ENDDO
      ELSE
         DO J1=1,N
            W(ISPT+NPT+J1)= STP*W(ISPT+NPT+J1) ! save the step taken
            W(IYPT+NPT+J1)= GRAD(J1)-W(J1)     ! save gradient difference: W(1:N) contains the old gradient
         ENDDO
      ENDIF
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      FIXIMAGE=.FALSE.
      IF (DUMPT.AND.DEBUG) THEN
         IF (AMBER) THEN
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,'(A,I4,A,F15.5)') 'At step number ',ITER,' energy=',ENERGY
            DO J2=1,NATOMS
               WRITE(DUMPXYZUNIT+NP,'(A,3F20.10)') typech(J2)(1:1),(XCOORDS(3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,'(A,I8,A,G20.10)') 'at step ',ITER,' energy=',ENERGY
            WRITE(DUMPXYZUNIT+NP,'(A2,3F20.10)') ('LA ',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=1,NATOMS-NS)
            IF (NS.GT.0) WRITE(DUMPXYZUNIT+NP,'(A2,3F20.10)') 
     1          ('LB',XCOORDS(3*(J1-1)+1),XCOORDS(3*(J1-1)+2),XCOORDS(3*(J1-1)+3),J1=NATOMS-NS+1,NATOMS)
         ENDIF
      ENDIF
      IF (CENT) CALL CENTRE2(XCOORDS)
      SMINKCURRENTP=SMINKCURRENT
      GOTO 10

      RETURN
      END
