C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
C        Line search removed plus small modifications, DJW 2001
C
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
      COMMON /MYPOT/ POTEL
      COMMON /Q4C/ QSTART, QFINISH
      LOGICAL EVAP, GUIDECHANGET, GUIDET, EVAPREJECT, SMINKCHANGET, CSMDOGUIDET
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /EV/ EVAP, evapreject
      SAVE W, DIAG, ITER, POINT, ISPT, IYPT, NPT

      IF (.NOT.ALLOCATED(DIAG)) ALLOCATE(DIAG(N))       ! SAVE doesn't work otherwise for Sun
      IF (.NOT.ALLOCATED(W)) ALLOCATE(W(N*(2*M+1)+2*M)) ! SAVE doesn't work otherwise for Sun
!     IF (QUENCHDOS) ALLOCATE(FRAMES(N,ITMAX), PE(ITMAX), MODGRAD(ITMAX))
      IF (SIZE(W,1).NE.N*(2*M+1)+2*M) THEN ! mustn't call mylbfgs with changing number of variables!!!
         WRITE(MYUNIT, '(A,I10,A,I10,A)') 'ERROR, dimension of W=',SIZE(W,1),' but N*(2*M+1)+2*M=',N*(2*M+1)+2*M,' in mylbfgs'
         call exit(10)
      ENDIF
      COREDONE=.FALSE.
      SMINKCHANGET=.FALSE.
      SMINKCURRENT=0.0D0
      SMINKCURRENTP=0.0D0
      LOCALSTEEREDMINT=.FALSE.
      IF (STEEREDMINT) LOCALSTEEREDMINT=.TRUE.

! for CHARMM: update nonbonded list at the start of each minimization
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
               CALL TIPIO(XCOORDS(3*(J2-1)+1),XCOORDS(3*(J2-1)+2),XCOORDS(3*(J2-1)+3),
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
      CALL POTENTIAL(XCOORDS,GRAD,ENERGY,.TRUE.,.FALSE.)
C
C  Catch cold fusion for ionic potentials and discard.
C
C  Changed EREAL for cold fusion to 1.0D6 rather than 0.0D0, which could result in steps being accepted
C  for systems with positive energies. - khs26 26/11/09
C
      IF ((TOSI.OR.WELCH.OR.RGCL2.OR.AMBER.OR.ARNO.OR.PACHECO.OR.TIP.OR.CHRMMT.OR.AMBERT 
     &   .OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT.OR.JMT) 
     &   .AND.(ENERGY.LT.COLDFUSIONLIMIT)) THEN
         WRITE(MYUNIT,'(A,G20.10)') 'ENERGY=',ENERGY
         ENERGY=1.0D6
         POTEL=1.0D6
         RMS=1.0D0
         WRITE(MYUNIT,'(A)') ' Cold fusion diagnosed - step discarded'
!     csw34> set COLDFUSION=.TRUE. so that ATEST=.FALSE. in MC
         COLDFUSION=.TRUE.
!        IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
         RETURN
      ENDIF

!     IF (QUENCHDOS) THEN
!        MODGRAD(1)=DSQRT(DDOT(N,GRAD,1,GRAD,1))
!        PE(1)=ENERGY
!        FRAMES(1:N,1)=XCOORDS(1:N)
!     ENDIF
!
! Stop the core from changing morphology easily, but allow it to relax
!
!     IF ((NCORE(NP).GT.0).AND.(.NOT.COREDONE)) THEN
!        DUMMY2=0.0D0
!        DO J1=1,3*(NATOMS-NCORE(NP))
!           DUMMY2=DUMMY2+GRAD(J1)**2
!        ENDDO
!        DUMMY2=DSQRT(DUMMY2/(3*NATOMS))
!        IF (DUMMY2.GT.EPS*2.0D0) THEN
!           GRAD(3*(NATOMS-NCORE(NP))+1:3*NATOMS)=0.0D0
!        ELSE
!           COREDONE=.TRUE.
!        ENDIF
!     ENDIF

      IF (CHRMMT .AND. GCHARMMFAIL) THEN
          WRITE(MYUNIT,'(A)') 'Failure in CHARMM energy/gradient evaluation - geometry discarded.'
!         IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
          RETURN
      ENDIF

C
C  If INTMINT and CHRMMT need to transform to internal coordinates
C  See COPTIM.2.3 for switching to internals from Cartesians using LIMINCUT.
C
      IF (INTMINT) THEN
         OLDCART(1:3*NATOMS)=XCOORDS(1:3*NATOMS) ! store cartesians in OLDCART for both CHARMM and UNRES
C         IF (UNRST) THEN
CC
CC store internals (in OLDQ) and update X to contain internals
CC
C            CALL geom_to_var(N,OLDQ)
C            XCOORDS(1:N)=OLDQ(1:N)
C         ELSE IF (CHRMMT) THEN
            CALL GETKD(KD) ! get width of sparse band in G matrix KD
            CALL GETNNZ(NNZ) ! get number of non-zero elements in B-matrix
            NOCOOR=.FALSE. ! calculate internals therefore NOCOOR is false
            GINT(1:N)=0.0D0 ! to prevent NaN's for Sun!
            XINT(1:N)=0.0D0 ! to prevent NaN's for Sun!
            CALL TRANSFORM(XCOORDS,GRAD,XINT,GINT,N,3*NATOMS,NNZ,NOCOOR,KD)
            OLDQ(1:N)=XINT(1:N)    ! store internals
            OLDGINT(1:N)=GINT(1:N) ! store gradient in internals
C         ENDIF
      ENDIF
C
C  for CHRMMT:
C  XCOORDS contains current Cartesians
C  GRAD    contains current gradient
C  XINT    contains current internals
C  GINT    contains current gradient in internals
C  OLDQ    contains internals for initial geometry
C  OLDGINT contains gradient in internals for initial geometry
C  OLDCART contains Cartesian coordinates for initial geometry
C
      IF (EVAPREJECT) RETURN
      POTEL=ENERGY

      IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'

C
C  Catch cold fusion for ionic potentials and discard.
C
      IF ((DBPT.OR.DBPTDT.OR.MSTBINT.OR.MSSTOCKT.OR.MULTPAHAT.OR.NPAHT.OR.PAHW99T.OR.PYGT.OR.TDHDT.OR.SILANET) 
     &   .AND.(ENERGY.LT.-5.0D4)) THEN
         WRITE(MYUNIT,'(A,G20.10)') 'ENERGY=',ENERGY
         ENERGY=0.0D0
         POTEL=0.0D0
         RMS=1.0D0
         WRITE(MYUNIT,'(A)') ' Cold fusion diagnosed - step discarded'
         RETURN
      ENDIF
C
C  Termination test. 
C
10    CALL FLUSH(MYUNIT)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
         IF (CHRMMT.AND.ACESOLV) THEN
            NCHENCALLS=ACEUPSTEP-1
            CALL POTENTIAL(XCOORDS,GRAD,ENERGY,.TRUE.,.FALSE.)
            IF (DEBUG) WRITE(*,'(A,2G20.10,A)') ' mylbfgs> Energy and RMS force=',ENERGY,RMS,' after ACE update'
            IF (RMS.LE.EPS) MFLAG=.TRUE.
         ELSE
            MFLAG=.TRUE.
         ENDIF
         IF (EVAP) MFLAG=.FALSE. ! do not allow convergence if we happen to have a small RMS and EVAP is true'
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            IF (DEBUG) WRITE(MYUNIT,'(A,F20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'

!             IF (QUENCHDOS) THEN
!                DO J1=1,ITDONE+1
!                   DIST=0.0D0
!                   DO J2=1,N
!                      DIST=DIST+(FRAMES(J2,J1)-FRAMES(J2,ITDONE+1))**2
!                   ENDDO
!                   DIST=SQRT(DIST)
!                   IF (J1.EQ.1) DIST1=DIST
! !                   DOSSTATS(J1,1)=PE(J1)
!                   IF ((MODGRAD(J1).GT.0.0D0).AND.(DIST.GT.0.0D0)) THEN
! !                    DOSSTATS(J1,2)=(MODGRAD(1)/MODGRAD(J1))*(DIST/DIST1)**(N-1)
! !                    DOSSTATS(J1,2)=(N-1)*LOG(DIST)-LOG(MODGRAD(J1))
!                      DOSSTATS(J1,2)=DIST**(N-1)/MODGRAD(J1)
!                   ELSE
!                      DOSSTATS(J1,2)=0.0D0
!                   ENDIF
! !                 WRITE(MYUNIT,'(A,I6,4G18.8)') 'lbfgs> J1,MODGRAD,DIST,DOSSTATS(J1,2),DOSSTATS(J1,1)=',
! !    &                      J1,MODGRAD(J1),DIST,DOSSTATS(J1,2),DOSSTATS(J1,1)
!                ENDDO
!                DEALLOCATE(FRAMES, PE, MODGRAD)
!             ENDIF

            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
         IF (DEBUG) FIXIMAGE=.FALSE.
         IF (DEBUG) WRITE(MYUNIT,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
!        IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
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
C           INQUIRE(FILE='diag',EXIST=YESNO)
C           IF (YESNO) THEN
C              OPEN(UNIT=34,FILE='diag',STATUS='OLD')
C              READ(34,*) (DIAG(J1),J1=1,N)
C              PRINT*,'diag read in LBFGS'
C              WRITE(*,'(6F15.5)') (DIAG(J1),J1=1,N)
C           ELSE
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
               WRITE(MYUNIT,'(A)') 'WARNING, resetting YY to one in mylbfgs'
               YY=1.0D0
            ENDIF
            IF (YS.EQ.0.0D0) THEN
               WRITE(MYUNIT,'(A)') 'WARNING, resetting YS to one in mylbfgs'
               YS=1.0D0
            ENDIF
            IF (DEBUG) WRITE(MYUNIT,'(A20,F20.5)') 'YY= ',YY
            IF (DEBUG) WRITE(MYUNIT,'(A20,F20.5)') 'YS= ',YS
C           WRITE(*,'(A,2F20.10)') 'YS/YY,STP=',YS/YY,STP
            DO J1=1,N
C              DIAG(J1)= ABS(YS/YY) ! messes up after step reversals!
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
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C
         CP= POINT
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

      IF (CHRMMT.AND.INTMINT) THEN
         DOT1=SQRT(DDOT(N,GINT,1,GINT,1))
      ELSE
         DOT1=SQRT(DDOT(N,GRAD,1,GRAD,1))
      ENDIF
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
         IF (CHRMMT.AND.INTMINT) THEN
            OVERLAP=DDOT(N,GINT,1,WTEMP,1)/(DOT1*DOT2)
         ELSE
            OVERLAP=DDOT(N,GRAD,1,WTEMP,1)/(DOT1*DOT2)
        ENDIF
      ENDIF
C     PRINT*,'OVERLAP,DIAG(1)=',OVERLAP,DIAG(1)
C     PRINT*,'GRAD . GRAD=',DDOT(N,GRAD,1,GRAD,1)
C     PRINT*,'W . W=',DDOT(N,W,1,W,1)
      IF (OVERLAP.GT.0.0D0) THEN
C        IF (DEBUG) PRINT*,'Search direction has positive projection onto gradient - resetting'
C        ITER=0
C        GOTO 10
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
C
C  We now have the proposed step.
C
      IF (CHRMMT.AND.INTMINT) THEN
         DO J1=1,N
            XINT(J1)=XINT(J1)+STP*W(ISPT+POINT*N+J1)
            DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)
         ENDDO
      ELSE
!
! Save XCOORDS here so that we can undo the step reliably including the
! non-linear projection for Thomson for the angular coordinates.
!
         XSAVE(1:N)=XCOORDS(1:N) 
         DO J1=1,N
            XCOORDS(J1)=XCOORDS(J1)+STP*W(ISPT+POINT*N+J1)
         ENDDO 
!
! For Thomson try projection for the geometry after the step.
!
         IF (PROJIT) THEN
            IF (THOMSONT) THEN
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
               CALL PROJI(TMPCOORDS,NATOMS)
               CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
               XCOORDS(1:N)=TMPANG(1:N)
            ELSE
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL PROJI(TMPCOORDS,NATOMS)
               XCOORDS(1:N)=TMPCOORDS(1:N)
            ENDIF
         ENDIF
         IF (PROJIHT) THEN
            IF (THOMSONT) THEN
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
               CALL PROJIH(TMPCOORDS,NATOMS)
               CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
               XCOORDS(1:N)=TMPANG(1:N)
            ELSE
               TMPCOORDS(1:N)=XCOORDS(1:N)
               CALL PROJIH(TMPCOORDS,NATOMS)
               XCOORDS(1:N)=TMPCOORDS(1:N)
            ENDIF
         ENDIF
      ENDIF
C
C  For charmm internals must transform and back-transform!
C
      NDECREASE=0
      LEPSILON=1.0D-6

20    IF (INTMINT) THEN
         IF (CHRMMT) THEN
            NEWQ(1:N)=OLDQ(1:N)
            LCART(1:3*NATOMS)=OLDCART(1:3*NATOMS)
C
C Need to keep OLDQ constant for repeated back-transformations if first step size fails.
C Therefore pass dummy array newq that can change.
C Similarly with LCART and OLDCART.
C
C           CALL TRANSBACK(XINT,NEWQ,LCART,N,3*NATOMS,NNZ,KD)
            CALL TRANSBACKDELTA(DELTAQ,DELTACART,LCART,N,3*NATOMS,NNZ,KD,FAILED,.FALSE.,LEPSILON) ! transform step to Cartesians
            IF (FAILED) THEN
C              NCOUNT=NCOUNT+1
C              IF (NCOUNT.GT.1) STOP
C              LEPSILON=1.0D-5*DPRAND()
C              GOTO 21
C or
               MFLAG=.FALSE.
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
               RETURN
            ENDIF
C
C now add DELTACART to LCART to get new cartesians. Put these in X.
C
            LCART(1:3*NATOMS)=OLDCART(1:3*NATOMS)+DELTACART(1:3*NATOMS)
            XCOORDS(1:3*NATOMS)=OLDCART(1:3*NATOMS)+DELTACART(1:3*NATOMS)
C
C  for CHRMMT:
C  LCART    contains new Cartesians (after step)
C  XCOORDS contains new Cartesians (after step)
C  XINT    contains new internals (after step)
C  GRAD    contains old gradient
C  GINT    contains old gradient in internals
C  OLDQ    contains old internals
C  OLDGINT contains old gradient in internals for the last successful geometry
C  NEWQ    contains old internals for the last successful geometry
C  OLDCART contains old Cartesians for the last successful geometry
C
C         ELSEIF (UNRST) THEN
C            NEWQ(1:N)=X(1:N) ! store new internals in NEWQ
CC
CC need a temporary array NEWQ here as argument to var_to_geom to keep X unchanged in case we need to
CC modify the step below.
CC
C            CALL var_to_geom(N,NEWQ) ! update internals
C            CALL chainbuild ! get cartesians
         ENDIF
      ENDIF

! csw34> INCREMENT THE FORCE CONSTANT FOR STEERED MINIMISATION
      SMINKCHANGET=.FALSE.
      IF (LOCALSTEEREDMINT) THEN
         SMINKCURRENT=MIN(SMINKCURRENT+SMINKINC,SMINK)
         IF (SMINKCURRENT.NE.SMINKCURRENTP) SMINKCHANGET=.TRUE.
! a bit of useful debug printing         
        IF (DEBUG) WRITE(MYUNIT,'(A,2F20.10,L5)') 'SMINKCURRENT,SMINKCURRENTP,SMINKCHANGET=',SMINKCURRENT,SMINKCURRENTP,SMINKCHANGET
      ENDIF

      CALL POTENTIAL(XCOORDS,GNEW,ENEW,.TRUE.,.FALSE.)

!     IF (QUENCHDOS) THEN
!        MODGRAD(ITDONE+2)=DSQRT(DDOT(N,GNEW,1,GNEW,1))
!        PE(ITDONE+2)=ENEW
!        FRAMES(1:N,ITDONE+2)=XCOORDS(1:N)
!     ENDIF
!
! Stop the core from changing morphology easily, but allow it to relax
!
!     IF ((NCORE(NP).GT.0).AND.(.NOT.COREDONE)) THEN
!        DUMMY2=0.0D0
!        DO J1=1,3*(NATOMS-NCORE(NP))
!           DUMMY2=DUMMY2+GNEW(J1)**2
!        ENDDO
!        DUMMY2=DSQRT(DUMMY2/(3*NATOMS))
!        IF (DUMMY2.GT.EPS*2.0D0) THEN
!           GNEW(3*(NATOMS-NCORE(NP))+1:3*NATOMS)=0.0D0
!        ELSE
!           COREDONE=.TRUE.
!        ENDIF
!     ENDIF

      IF (EVAPREJECT) return
      IF (CHRMMT .AND. GCHARMMFAIL) THEN
          WRITE(MYUNIT,'(A)') 'Failure in CHARMM energy/gradient evaluation - step discarded.'
!         IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
          RETURN
      ENDIF
C
C  Catch cold fusion for ionic potentials and discard.
C
C  Changed EREAL for cold fusion to 1.0D6 rather than 0.0D0, which could result in steps being accepted
C  for systems with positive energies. - khs26 26/11/09
C
      IF ((TOSI.OR.WELCH.OR.RGCL2.OR.AMBER.OR.ARNO.OR.PACHECO.OR.TIP.OR.CHRMMT.OR.AMBERT 
     &   .OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT.OR.JMT)
     &   .AND.(ENEW.LT.COLDFUSIONLIMIT)) THEN
         ENERGY=1.0D6
         ENEW=1.0D6
         POTEL=1.0D6
         RMS=1.0D0
         WRITE(MYUNIT,'(A)') ' Cold fusion diagnosed - step discarded'
!     csw34> set COLDFUSION=.TRUE. so that ATEST=.FALSE. in MC
         COLDFUSION=.TRUE.
!        IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
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


C
C  We need to transform the newly obtained Cartesian gradient for CHARMM and internals.
C  NOCOOR is true because we dont need to transform the coordinates.
C
      IF (CHRMMT.AND.INTMINT) THEN
         NOCOOR=.TRUE.
         CALL TRANSFORM(XCOORDS,GNEW,XINT,GINT,N,3*NATOMS,NNZ,NOCOOR,KD)
      ENDIF

C     IF (TIP) THEN
C           WRITE(DUMPXYZUNIT+NP,'(I6)') (NATOMS/2)*3
C           WRITE(DUMPXYZUNIT+NP,'(A,I5,A,F20.10)') 'LBFGS iteration ',ITER,' energy =',ENEW
C           DO J2=1,NATOMS/2
C              CALL TIPIO(XCOORDS(3*(J2-1)+1),XCOORDS(3*(J2-1)+2),XCOORDS(3*(J2-1)+3),
C    1              XCOORDS(3*(NATOMS/2+J2-1)+1),XCOORDS(3*(NATOMS/2+J2-1)+2),XCOORDS(3*(NATOMS/2+J2-1)+3),RBCOORDS)
C              WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
C              WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
C              WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
C           ENDDO
C     ENDIF

C     WRITE(*,'(A,F20.10)') 'ENEW=',ENEW
C     WRITE(*,'(I6,F20.10)') (J1,GNEW(J1),J1=1,N)

C
C csw34 Force acceptance of step if FIXDIHEFLAG is TRUE
C
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
C
C  Step finished so can reset OLDQ to new XINT, OLDCART to new LCART,
C  as well as the Cartesian and internal gradients.
C
         IF (CHRMMT.AND.INTMINT) THEN
            OLDGINT(1:N)=GINT(1:N)
            OLDCART(1:3*NATOMS)=LCART(1:3*NATOMS)
C
C  Need to remake XINT because step was only projected in Cartesians?
C  Actually, just setting OLDQ=XINT without this correction seems to
C  be OK. Due to numerical imprecision, it might still be possible
C  for X and XINT to get out of register. Perhaps this doesn't matter
C  because the energy and gradient are always calculated in Cartesians.
C
C           IF (BFGSTST) CALL TRANSDELTA(DELTACART,DELTAQ,LCART,N,3*NATOMS,NNZ,KD)
C           OLDQ(1:N)=OLDQ(1:N)+DELTAQ(1:N)
            OLDQ(1:N)=XINT(1:N)
C         ELSEIF (UNRST) THEN
C!           TEST1(1:N)=X(1:N)
C            CALL geom_to_var(N,X(1:N)) ! testing!!! - to put X back into register with the common block internals (and g)
C!           CALL geom_to_var(N,TEST1(1:N))
C!           do j1=1,N
C!           if (abs((TEST1(j1)-x(j1))/x(j1))*100.0d0.gt.1.0D-6) print *,'hello coords ',J1
C!           enddo
         ENDIF
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
      ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
C
C  Energy decreased too much - try again with a smaller step size
C
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find an energy in the required range, NFAIL=',NFAIL
            IF (CHRMMT.AND.INTMINT) THEN ! need to reset X, XINT, G, GINT to original values
               XINT(1:N)=XINT(1:N)-STP*W(ISPT+POINT*N+1:ISPT+POINT*N+N)
C              XINT=OLDQ ! should be the same as subtracting the step
               GINT(1:N)=OLDGINT(1:N)
               GRAD(1:3*NATOMS)=GNEW(1:3*NATOMS) ! here OPTIM uses GLAST ! DJW
               XCOORDS(1:3*NATOMS)=OLDCART(1:3*NATOMS)
            ELSE
!
! Resetting to XSAVE should be the same as subtracting the step. 
! If we have tried PROJI with Thomson then the projection is non-linear
! and we need to reset to XSAVE. This should always be reliable!
!
               XCOORDS(1:N)=XSAVE(1:N)
               GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point

!              DO J1=1,N
!                 XCOORDS(J1)=XCOORDS(J1)-STP*W(ISPT+POINT*N+J1)
!              ENDDO

            ENDIF
            ITER=0   !  try resetting
            IF (NFAIL.GT.20) THEN
               WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
               FIXIMAGE=.FALSE.
!              STOP
!              IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)
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
!
! For Thomson try projection for the geometry after the step.
!        
            IF (PROJIT) THEN
               IF (THOMSONT) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF
            IF (PROJIHT) THEN
               IF (THOMSONT) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJIH(TMPCOORDS,NATOMS)
                  CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL PROJIH(TMPCOORDS,NATOMS)
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
C
C  Energy increased - try again with a smaller step size
C
         IF (NDECREASE.GT.10) THEN ! DJW
            NFAIL=NFAIL+1
            WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL
            IF (CHRMMT.AND.INTMINT) THEN ! need to reset X, XINT, G, GINT to original values
               XINT(1:N)=XINT(1:N)-STP*W(ISPT+POINT*N+1:ISPT+POINT*N+N)
C              XINT=OLDQ ! should be the same as subtracting the step
               GINT(1:N)=OLDGINT(1:N)
               GRAD(1:3*NATOMS)=GNEW(1:3*NATOMS) ! here OPTIM uses GLAST ! DJW
               XCOORDS(1:3*NATOMS)=OLDCART(1:3*NATOMS)
            ELSE
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
            ENDIF
            ITER=0   !  try resetting
!            IF (NFAIL.GT.20) THEN
! bs360: smaller NFAIL 
             IF (NFAIL.GT.5) THEN         
               WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
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
         IF (CHRMMT.AND.INTMINT) THEN
            DO J1=1,N
               XINT(J1)=XINT(J1)-0.9*STP*W(ISPT+POINT*N+J1)
               DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)*0.1D0
            ENDDO
         ELSE
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

!           DO J1=1,N
!              XCOORDS(J1)=XCOORDS(J1)-0.9*STP*W(ISPT+POINT*N+J1)
!           ENDDO 
!
! For Thomson try projection for the geometry after the step.
!        
            IF (PROJIT) THEN
               IF (THOMSONT) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL PROJI(TMPCOORDS,NATOMS)
                  XCOORDS(1:N)=TMPCOORDS(1:N)
               ENDIF
            ENDIF
            IF (PROJIHT) THEN
               IF (THOMSONT) THEN
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL THOMSONANGTOC(TMPCOORDS,NATOMS)
                  CALL PROJIH(TMPCOORDS,NATOMS)
                  CALL THOMSONCTOANG(TMPCOORDS,TMPANG,NATOMS)
                  XCOORDS(1:N)=TMPANG(1:N)
               ELSE
                  TMPCOORDS(1:N)=XCOORDS(1:N)
                  CALL PROJIH(TMPCOORDS,NATOMS)
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
C
C     Compute the new step and gradient change
C
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
!     IF (QUENCHDOS) DEALLOCATE(FRAMES, PE, MODGRAD)

      RETURN
      END
