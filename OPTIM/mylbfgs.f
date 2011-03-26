C   GPL License Info {{{      
C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C}}}
C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
C        Line search removed plus small modifications, DJW 2001
C        If INTMINT is true, then N is NINTS, 
C        otherwise N = NOPT JMC
C        changed declaration of X(N) to X(3*NATOMS) 30/4/04
C        X is passed in and out in Cartesians.
C
C  Cartesian coordinate and gradient vectors are declared 3*NATOMS - the complication
C  is that for internal coordinate optimisations we can have a number of degrees of
C  freedom that is more or less than 3*NATOMS. N should specify this dimension.
C
      SUBROUTINE MYLBFGS(N,M,X,DIAGCO,MFLAG,ENERGY,RMS,EREAL,REALRMS,ITMAX,
     1                   RESET,ITDONE,PTEST,GSAVE,NODUMP,PROJECT)
      ! Declarations {{{
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE MODHESS
      USE ZWK
      USE MODUNRES
      USE MODCHARMM
      USE MODAMBER9, ONLY : NRESPA2, IRESPA2
      use PORFUNCS
      USE SPFUNCTS, ONLY : DUMPCOORDS
      
      IMPLICIT NONE
      INTEGER N,M,J1,J2,ITMAX,ITDONE,NFAIL,NCOUNT
C     DOUBLE PRECISION X(*),G(3*NATOMS),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT,VEC2(3*NATOMS),OVERLAP
      DOUBLE PRECISION X(3*NATOMS),G(3*NATOMS),SLENGTH,DDOT,VEC2(3*NATOMS),OVERLAP,DISTF,DISTS,GAMMA
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W
      DOUBLE PRECISION DUMMY1,ENERGY,ENEW,RMS,EREAL,REALRMS,RVEC(3*NATOMS),ALPHA,GSAVE(3*NATOMS),DUMMY2
      LOGICAL DIAGCO, RESET, PTEST, NODUMP, PROJECT
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA
      LOGICAL NOCOOR, NODERV ! internals stuff
      DOUBLE PRECISION DELTAQ(N),DELTACART(3*NATOMS)
      DOUBLE PRECISION CART(3*NATOMS),OLDQ(N),NEWQ(N),OLDGINT(N),GINT(N),XINT(N) ! JMC
      DOUBLE PRECISION OLDCART(3*NATOMS),TMPINT(NINTS) ! JMC
      DOUBLE PRECISION GLAST(3*NATOMS),XSAVE(N)
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,I1,NREV
      INTEGER KD, NNZ
      LOGICAL MFLAG,FAILED
      DOUBLE PRECISION GSQSCALE, GSTHRESH, DOT1, DOT2
      INTEGER NSPECIAL, NALLOW, NINFO, FRAME
      LOGICAL PVFLAG
      COMMON /PVF/ PVFLAG
      COMMON /G2/ GSTHRESH, GSQSCALE, NSPECIAL, NALLOW, NINFO
      LOGICAL PATHT, DRAGT
      INTEGER NPATHFRAME, NDECREASE, ISTAT
      COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      CHARACTER ESTRING*87, GPSTRING*80, NSTRING*80, FSTRING*80
      COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING
      LOGICAL RADMOVED, OVERLP
      COMMON /DISCON/ RADMOVED, OVERLP
      LOGICAL PUSH, PULL
      COMMON /MORPHDATA/ PUSH, PULL
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE

C
C  numerical derivative for g^2
C
C     DOUBLE PRECISION SHIT, EPLUS, EMINUS, GDUM(3*NATOMS), RMSDUM, HDUM(3*NATOMS,3*NATOMS)
C     LOGICAL FIXSAVE
C
C  SGI appears to need this SAVE statement!
C
      SAVE W, DIAG, ITER, POINT, ISPT, IYPT, NPT
!}}}

      IF (.NOT.ALLOCATED(DIAG)) ALLOCATE(DIAG(N))       ! SAVE doesn't work otherwise for Sun
      IF (.NOT.ALLOCATED(W)) ALLOCATE(W(N*(2*M+1)+2*M)) ! SAVE doesn't work otherwise for Sun
      IF (SIZE(W,1).NE.N*(2*M+1)+2*M) THEN ! mustn't call mylbfgs with changing number of variables!!!
         PRINT '(A,I10,A,I10,A)', ' mylbfgs> ERROR, dimension of W=',SIZE(W,1),' but N*(2*M+1)+2*M=',N*(2*M+1)+2*M
         STOP
      ENDIF

      IF (N.NE.3*NATOMS) THEN
         IF ((.NOT.UNRST).AND.(.NOT.(CHRMMT.AND.INTMINT)).AND.(.NOT.VARIABLES).AND.(.NOT.RINGPOLYMERT)) THEN
            PRINT*,'ERROR - N and 3*NATOMS are different in mylbfgs: ',N,3*NATOMS
            STOP
         ENDIF
      ENDIF
! for CHARMM: update nonbonded list at the start of each minimization
      IF(CHRMMT) CALL UPDATENBONDS(X)

      ALPHA=1.0D0
      NFAIL=0
      FRAME=1
      CFUSIONT=.FALSE.
      IF (RESET) ITER=0
      ITDONE=0
      IF (RESET.AND.PTEST) WRITE(*,'(A)') ' mylbfgs> Resetting LBFGS minimiser'
      IF ((.NOT.RESET).AND.PTEST) WRITE(*,'(A)') ' mylbfgs> Not resetting LBFGS minimiser'
1     FIXIMAGE=.FALSE.
      IF (PV.AND.(.NOT.PROJECT)) THEN
         IF (.NOT.KNOWE) THEN 
            CALL POTENTIAL(X,ENERGY,GSAVE,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)

! Check for cold fusion
            PRINT*, "msb50 in mylbfgs, energy", energy
            if (ENERGY.LT.coldFusionLimit) then
               WRITE(*,'(A,2G20.10)') ' mylbfgs> Cold fusion diagnosed - step discarded, energy, limit=',ENERGY,coldFusionLimit
               ENERGY=1.0d60
               EREAL=1.0d60
               RMS=1.0d1
               CFUSIONT=.TRUE.
               RETURN
            endif
         ENDIF
         PVFLAG=.FALSE.
         CALL PVOPT(X,ENERGY,GSAVE)
      ENDIF

      IF (UNRST) THEN
         DO J2=1,nres
            c(1,J2)=X(6*(J2-1)+1)
            c(2,J2)=X(6*(J2-1)+2)
            c(3,J2)=X(6*(J2-1)+3)
            c(1,J2+nres)=X(6*(J2-1)+4)
            c(2,J2+nres)=X(6*(J2-1)+5)
            c(3,J2+nres)=X(6*(J2-1)+6)
         END DO
         CALL UPDATEDC
         CALL int_from_cart(.true.,.false.)
C
C jmc 9/1/03 Need this call to chainbuild to calculate xrot,xloc matrices which are used
C in the calculation of d(bond vectors) by d(internals).
C
         CALL chainbuild
      ENDIF

      IF ((.NOT.KNOWE).OR.(.NOT.KNOWG)) THEN
         IF(AMBERT.OR.NABT) irespa2=ITDONE+1
         RMS=0.0D0
         CALL POTENTIAL(X,ENERGY,GSAVE,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)
         REALRMS=RMS
         IF (ENERGY.LT.COLDFUSIONLIMIT) then
            WRITE(*,'(A,2G20.10)') ' mylbfgs> Cold fusion diagnosed - step discarded, energy, limit=',ENERGY,coldFusionLimit
            ENERGY=1.0d60
            EREAL=1.0d60
            RMS=1.0d1
            REALRMS=RMS
            CFUSIONT=.TRUE.
            RETURN
         endif
      ELSE
         RMS=REALRMS
      ENDIF

      G(1:3*NATOMS)=GSAVE(1:3*NATOMS)
      GLAST(1:3*NATOMS)=GSAVE(1:3*NATOMS)
      IF (.NOT.(NODUMP)) CALL DUMPP(X,ENERGY)


      IF (TWOENDS.AND.(.NOT.TTDONE).AND.(FORCE.NE.0.0D0)) THEN
         IF (CHRMMT.AND.INTMINT) THEN
            PRINT*,'ERROR - TWOENDS not available for CHARMM internal coordinates'
            STOP
         ENDIF
         DUMMY1=0.0D0
         DO J1=1,N
            RVEC(J1)=FIN(J1)-X(J1)
            DUMMY1=DUMMY1+RVEC(J1)**2
         ENDDO
         DUMMY1=1.0D0/SQRT(DUMMY1)
         IF (1.0D0/DUMMY1.GT.0.05D0) THEN
            RMS=0.0D0
            DO J1=1,N
               G(J1)=G(J1)+FORCE*RVEC(J1)*DUMMY1
               RMS=RMS+G(J1)**2
            ENDDO
            RMS=SQRT(RMS/N)
         ENDIF
      ELSE IF (DRAGT) THEN
         IF (CHRMMT.AND.INTMINT) THEN
            PRINT*,'ERROR - DRAG not available for CHARMM internal coordinates'
            STOP
         ENDIF
         IF (N.GT.3*NATOMS) THEN
            PRINT*,'ERROR N > 3*NATOMS in mylbfgs'
            STOP
         ENDIF
         DUMMY1=0.0D0
         DO J1=1,N
            RVEC(J1)=FIN(J1)-X(J1)
            DUMMY1=DUMMY1+RVEC(J1)**2
         ENDDO
         DUMMY1=1.0D0/SQRT(DUMMY1)
         WRITE(*,'(A,F20.10)') 'DIST=',1.0D0/DUMMY1
         DO J1=1,N
            RVEC(J1)=RVEC(J1)*DUMMY1
         ENDDO
         DUMMY1=0.0D0
         DO J1=1,N
            DUMMY1=DUMMY1+RVEC(J1)*G(J1)
         ENDDO
         PRINT*,'Projection of gradient=',DUMMY1
C        IF (DUMMY1.GT.0.0D0) THEN
            RMS=0.0D0
            DO J1=1,N
               G(J1)=G(J1)-DUMMY1*RVEC(J1)-0.1D0*RVEC(J1)
               RMS=RMS+G(J1)**2
            ENDDO
            RMS=SQRT(RMS/N)
C        ENDIF
      ELSE IF (PROJECT) THEN
C
C  For CHARMM internal coordinate minimisation we project the uphill direction
C  out of the Cartesian gradient every time. If MYLBFGS is reset on each call then
C  the total step should be a linear combination of gradients that all have this
C  component removed. However - if we don't reset then on previous calls the
C  uphill direction will be different! Suggests that we should reset in bfgsts.
C
         DO J2=1,NUP
            IF (FREEZE) THEN
               DO J1=1,NATOMS
                  IF (.NOT.FROZEN(J1)) CYCLE
                  ZWORK(3*(J1-1)+1,J2)=0.0D0
                  ZWORK(3*(J1-1)+2,J2)=0.0D0
                  ZWORK(3*(J1-1)+3,J2)=0.0D0
               ENDDO
            ENDIF
            DUMMY2=0.0D0
            DO J1=1,N
               DUMMY2=DUMMY2+ZWORK(J1,J2)**2
            ENDDO
!           PRINT '(A,G20.10)','  mylbfgs> uphill vector mod squared=',DUMMY2
            IF (ABS(DUMMY2-1.0D0).GT.1.0D-10) THEN
               DUMMY2=1.0D0/SQRT(DUMMY2)
               PRINT '(A,G20.10)',' mylbfgs> renormalising uphill vector by factor of ',DUMMY2
               DO J1=1,N
                  ZWORK(J1,J2)=ZWORK(J1,J2)*DUMMY2
               ENDDO
            ENDIF
            DUMMY1=0.0D0
            IF (UNRST) THEN ! for UNRST N < 3*NATOMS
               DO J1=1,N
                  DUMMY1=DUMMY1+ZWORK(J1,J2)*G(J1)
               ENDDO
               RMS=0.0D0
               DO J1=1,N
                  G(J1)=G(J1)-ZWORK(J1,J2)*DUMMY1
                  RMS=RMS+G(J1)**2
               ENDDO
               RMS=SQRT(RMS/N)
            ELSE
               DUMMY2=0.0D0
               DO J1=1,N
                  DUMMY1=DUMMY1+ZWORK(J1,J2)*G(J1)
               ENDDO
               RMS=0.0D0
               DO J1=1,N
                  G(J1)=G(J1)-ZWORK(J1,J2)*DUMMY1
                  RMS=RMS+G(J1)**2
               ENDDO
               RMS=SQRT(RMS/N)
            ENDIF
         ENDDO
C        CALL ORTHOGOPT(W,X,.FALSE.)
      ENDIF
C
C  If INTMINT and CHRMMT need to transform to internal coordinates
C  See COPTIM.2.3 for switching to internals from Cartesians using LIMINCUT.
C
      IF (INTMINT) THEN
         OLDCART(1:3*NATOMS)=X(1:3*NATOMS) ! store cartesians in OLDCART for both CHARMM and UNRES
         IF (UNRST) THEN
C
C store internals (in OLDQ) and update X to contain internals
C
            CALL geom_to_var(N,OLDQ)
            X(1:N)=OLDQ(1:N)
         ELSE IF (CHRMMT) THEN 
            CALL GETKD(KD) ! get width of sparse band in G matrix KD
            CALL GETNNZ(NNZ) ! get number of non-zero elements in B-matrix
            NOCOOR=.FALSE. ! calculate internals therefore NOCOOR is false
            NODERV = .FALSE.
            GINT(1:N)=0.0D0 ! to prevent NaN's for Sun!
            XINT(1:N)=0.0D0 ! to prevent NaN's for Sun!
            CALL TRANSFORM(X,G,XINT,GINT,N,3*NATOMS,NNZ,NOCOOR,NODERV,KD,INTEPSILON)
            OLDQ(1:N)=XINT(1:N)    ! store internals
            OLDGINT(1:N)=GINT(1:N) ! store gradient in internals
         ELSE IF (AMBERT.OR.NABT) THEN
            PRINT*, "not yet set up for amber"
            STOP
         ENDIF
      ENDIF
C
C  for CHRMMT:
C  X       contains current Cartesians
C  G       contains current gradient
C  XINT    contains current internals
C  GINT    contains current gradient in internals
C  OLDQ    contains internals for initial geometry
C  OLDGINT contains gradient in internals for initial geometry
C  OLDCART contains Cartesian coordinates for initial geometry
C
      IF (GRADSQ) THEN
         IF (CHRMMT.AND.INTMINT) THEN
            PRINT*,'ERROR - GRADSQ minimisation incompatible with internal coordinates'
            STOP
         ENDIF
         EREAL=ENERGY
         REALRMS=RMS
         ENERGY=DDOT(3*NATOMS,G,1,G,1)
C        ENERGY=SQRT(DDOT(3*NATOMS,G,1,G,1))
         CALL DSYMV('U',3*NATOMS,2.0D0,HESS,SIZE(HESS,1),G,1,0.0D0,VEC2,1)
         RMS=DSQRT(DDOT(3*NATOMS,VEC2,1,VEC2,1)/(3*NATOMS))
         IF (RMS.LT.GSTHRESH) THEN
            IF (NSPECIAL.GT.0) THEN
               IF (MOD(ITER+1,NSPECIAL).EQ.0) THEN
                  CALL G2SPECIAL(NATOMS,X,G,VEC2,ENERGY,RMS,EREAL,REALRMS,MFLAG)
               ENDIF
            ELSE IF (FIXAFTER.EQ.0) THEN
               FIXAFTER=ITER+1
               PRINT*,'mylbfgs> changing to atom type LC and setting FIXAFTER=',ITER+1
               DO J1=1,NATOMS
                  ZSYM(J1)='LC'
               ENDDO
               GOTO 1
            ENDIF
         ENDIF
         DO J1=1,3*NATOMS
            G(J1)=VEC2(J1)
C           G(J1)=VEC2(J1)/ENERGY
         ENDDO
      ELSE
         EREAL=ENERGY
         REALRMS=RMS
      ENDIF

      IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A)') 
     1             ' mylbfgs> Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps'
      WRITE(ESTRING,16) ' mylbfgs> Energy for last cycle=',ENERGY,' '
16    FORMAT(A,27X,F20.10,A)

10    CALL FLUSH(6,ISTAT)

      MFLAG=.FALSE.
      IF (RMS.LE.GMAX) THEN ! GMAX is in key module, so can be changed by changep call
         IF (CHRMMT.AND.ACESOLV) THEN
            NCHENCALLS=ACEUPSTEP-1
            CALL POTENTIAL(X,ENERGY,GSAVE,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            IF(DEBUG) WRITE(*,'(A,2G20.10,A)')
     $           ' mylbfgs> Energy and RMS force=',ENERGY,RMS,
     1           ' after ACE update'
            IF (RMS.LE.GMAX) MFLAG=.TRUE.
!         ELSE IF (AMBERT.OR.NABT) THEN
!            Born radii update disabled for the moment, until we figure out what the hell is going wrong
!            IRESPA2=NRESPA2
!            RMS=0.0D0
!            CALL POTENTIAL(X,ENERGY,GSAVE,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!            IF(DEBUG) WRITE(*,'(A,2G20.10,A)')
!     $           ' mylbfgs> Energy and RMS force=',ENERGY,RMS,
!     1           ' after Born radii update'
!            IF(RMS.LE.GMAX) MFLAG=.TRUE.
         ELSE
            MFLAG=.TRUE.
         ENDIF
         IF (ITDONE.LT.NSTEPMIN) MFLAG=.FALSE.
         IF (PV.AND.(.NOT.PVFLAG)) MFLAG=.FALSE.
         IF (MFLAG) THEN
            IF (GRADSQ.AND.PTEST) WRITE(*,'(A,4F20.10)')
     $           ' mylbfgs> g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS
            FIXIMAGE=.FALSE.
            IF (INTMINT) THEN
C
C jmc put Cartesians in X for return
C
               IF (UNRST) THEN
                  TMPINT(1:N)=X(1:N)
                  CALL var_to_geom(N,TMPINT) ! jmc update internals
                  CALL chainbuild ! get cartesians
                  DO I1=1,nres
                     X(6*(I1-1)+1)=c(1,I1)
                     X(6*(I1-1)+2)=c(2,I1)
                     X(6*(I1-1)+3)=c(3,I1)
                     X(6*(I1-1)+4)=c(1,I1+nres)
                     X(6*(I1-1)+5)=c(2,I1+nres)
                     X(6*(I1-1)+6)=c(3,I1+nres)
                  ENDDO
               ENDIF
            ENDIF
C           WRITE(*,'(A,F20.10)') ' mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
C           DGUESS=DIAG(1) ! saved for subsequent calls - should be OK for the same system?
C                          ! May make redopath runs unreproducible?
            IF (PTEST) WRITE(*,'(A,G25.17)') ' mylbfgs> Final energy is ',ENERGY
            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
         FIXIMAGE=.FALSE.
         IF (INTMINT) THEN
            IF (UNRST) THEN
               TMPINT(1:N)=X(1:N) ! update internals
               CALL var_to_geom(N,TMPINT)
               CALL chainbuild ! get cartesians
               DO I1=1,nres
                  X(6*(I1-1)+1)=c(1,I1)
                  X(6*(I1-1)+2)=c(2,I1)
                  X(6*(I1-1)+3)=c(3,I1)
                  X(6*(I1-1)+4)=c(1,I1+nres)
                  X(6*(I1-1)+5)=c(2,I1+nres)
                  X(6*(I1-1)+6)=c(3,I1+nres)
               ENDDO
            ENDIF
         ENDIF
C        WRITE(*,'(A,F20.10)') ' mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
C        DGUESS=DIAG(1) ! saved for subsequent calls - should be OK for the same system?
C                          ! May make redopath runs unreproducible?
         RETURN
      ENDIF

      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(*,240)
 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            PRINT*,'using estimate of the inverse diagonal elements'
            DO I=1,N
               IF (DIAG(I).LE.0.0D0) THEN
                  WRITE(*,235) I
 235              FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
     1                   ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
                  STOP
               ENDIF
            ENDDO
         ELSE
C           INQUIRE(FILE='diag',EXIST=YESNO)
C           IF (YESNO) THEN
C              OPEN(UNIT=34,FILE='diag',STATUS='OLD')
C              READ(34,*) (DIAG(I),I=1,N)
C              PRINT*,'diag read in LBFGS'
C              WRITE(*,'(6F15.5)') (DIAG(I),I=1,N)
C           ELSE
            DO I=1,N
               DIAG(I)=DGUESS
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
C  NR step for diagonal inverse Hessian
C
         IF (CHRMMT.AND.INTMINT) THEN
            DO I=1,N
               W(ISPT+I)= -GINT(I)*DIAG(I)
               W(I)= -GINT(I)*DIAG(I)
            ENDDO
            GNORM= DSQRT(DDOT(N,GINT,1,GINT,1))
         ELSE
            DO I=1,N
               W(ISPT+I)= -G(I)*DIAG(I)
               W(I)= -G(I)*DIAG(I)
            ENDDO
            GNORM= DSQRT(DDOT(N,G,1,G,1))
         ENDIF
C
C  Make the first guess for the step length cautious.
C
         IF (GNORM.EQ.0.0D0) THEN
            GNORM=1.0D0 ! exact zero is presumably wrong!
            PRINT '(A)','WARNING - GNORM was zero in mylbfgs, resetting to one'
         ENDIF
         STP=MIN(1.0D0/GNORM,GNORM)
C        STP=1.0D0
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
C        PRINT*,'before overlap W, W: ITER,M,ISPT,IYPT,NPT=',ITER,M,ISPT,IYPT,NPT
C        WRITE(*,'(I5,2E20.10)') (J1,W(ISPT+NPT+J1),W(IYPT+NPT+J1),J1=1,10)
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
C        WRITE(*,'(A,E20.10)') 'YS=',YS
         IF (YS.EQ.0.0D0) YS=1.0D0
C
C  Update estimate of diagonal inverse Hessian elements
C  We divide by both YS and YY at different points, so
C  they had better not be zero!
C  W(ISPT+NPT+1:ISPT+NPT+N) stores the previous step vector
C  W(IYPT+NPT+1:IYPT+NPT+N) stores the previous difference of gradient vectors
C
         IF (.NOT.DIAGCO) THEN
C
C  Scaling individual components
C
            IF (.FALSE.) THEN
               YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
               IF (YY.EQ.0.0D0) YY=1.0D0
               GAMMA=YS/YY
               DO I=1,N
                  DUMMY1=0.0D0
                  DUMMY2=0.0D0
                  DO J1=1,M
                     DUMMY1=DUMMY1+W(ISPT+N*(J1-1)+I)*W(IYPT+N*(J1-1)+I)
                     DUMMY2=DUMMY2+W(IYPT+N*(J1-1)+I)**2
                  ENDDO
                  IF (DUMMY2.LT.1.0D-5) THEN
                     DIAG(I)=GAMMA
                  ELSE 
                     DIAG(I)=DUMMY1/DUMMY2
                  ENDIF
               ENDDO
            ELSE
               YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
C              WRITE(*,'(A,E20.10)') 'YY=',YY
               IF (YY.EQ.0.0D0) YY=1.0D0
               DUMMY1=YS/YY
C              DUMMY1=ABS(YS/YY) ! if the previous step was reversed this ABS seems a bad idea!  DJW
                              ! An alternative would be to reverse the gradient difference as well.
                              ! However, YS is the dot product of the gradient change with the step
                              ! so if we don't take the ABS value the signs should probably be right.
C              WRITE(*,'(A,E20.10)') 'DUMMY1=',DUMMY1
               DO I=1,N
                  DIAG(I)=DUMMY1
               ENDDO
            ENDIF
         ELSE
            PRINT*,'using estimate of the inverse diagonal elements'
            DO I=1,N
               IF (DIAG(I).LE.0.0D0) THEN
                  WRITE(*,235) I
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
C        PRINT*,'W(I) gets set to -G(I):'
C        WRITE(*,'(I5,2E20.10)') (J1,W(J1),G(J1),J1=1,10)
         IF (CHRMMT.AND.INTMINT) THEN
            DO I=1,N
               W(I)= -GINT(I)
            ENDDO
         ELSE
            DO I=1,N
               W(I)= -G(I)
            ENDDO
         ENDIF
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
C  If this is a BFGSTST or MORPHT  run project out the uphill direction.
C  For CHARMM internals we project the Cartesian step below.
C
!     PRINT '(A,6G20.10)','DIAG=',DIAG(1:3)
!     PRINT '(A,6G20.10)','COORDS=',X(1:3)
!     PRINT '(A,6G20.10)','G=',G(1:3)
!     PRINT '(A,6G20.10)','W=',W(1:3)
!     PRINT '(A,6G20.10)','ZWORK=',ZWORK(1:3,1)
      IF ((PROJECT).AND.(.NOT.TWOENDS).AND.(.NOT.(CHRMMT.AND.INTMINT))) THEN
         DO J2=1,NUP

            IF (FREEZE) THEN ! this may duplicate the block in the projection of the gradient
               DO J1=1,NATOMS
                  IF (.NOT.FROZEN(J1)) CYCLE
                  ZWORK(3*(J1-1)+1,J2)=0.0D0
                  ZWORK(3*(J1-1)+2,J2)=0.0D0
                  ZWORK(3*(J1-1)+3,J2)=0.0D0
               ENDDO
            ENDIF

            DUMMY1=0.0D0
            DO J1=1,N
               DUMMY1=DUMMY1+ZWORK(J1,J2)*W(J1)
            ENDDO
            DO J1=1,N
               W(J1)=W(J1)-ZWORK(J1,J2)*DUMMY1
            ENDDO
         ENDDO
      ENDIF
C
C  Store the new search direction
C
C     PRINT*,'W(I):'
C     WRITE(*,'(I5,E20.10)') (J1,W(J1),J1=1,10)
      IF (ITER.GT.0) THEN
         DO I=1,N
            W(ISPT+POINT*N+I)= W(I)
         ENDDO
      ENDIF

!     PRINT*,'before overlap test ITER, DIAG(1)=',ITER, DIAG(1)
!     PRINT*,'before overlap test  X, G, W:'
!     WRITE(*,'(I5,3E20.10)') (J1,X(J1),G(J1),W(ISPT+POINT*N+J1),J1=1,N)

      IF (CHRMMT.AND.INTMINT) THEN
         DOT1=SQRT(DDOT(N,GINT,1,GINT,1))
      ELSE
         DOT1=SQRT(DDOT(N,G,1,G,1))
      ENDIF
      DOT2=SQRT(DDOT(N,W,1,W,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) THEN
         IF (CHRMMT.AND.INTMINT) THEN
            OVERLAP=DDOT(N,GINT,1,W,1)/(DOT1*DOT2)
         ELSE
            OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)
        ENDIF
      ENDIF
C     PRINT*,'OVERLAP,DIAG(1)=',OVERLAP,DIAG(1)
C     PRINT*,'G . G=',DDOT(N,G,1,G,1)
C     PRINT*,'W . W=',DDOT(N,W,1,W,1)
C
C  The step is saved in W(ISPT+POINT*N+1:ISPT+POINT*N+N). 
C  W(1:N) is overwritten by the gradient.
C
      IF (OVERLAP.GT.0.0D0) THEN
         IF (PTEST) PRINT '(A,G20.10,A)','Search direction has positive projection onto gradient ',OVERLAP,' reversing step'
         DO I=1,N
            W(ISPT+POINT*N+I)= -W(I)  ! if we reverse the step it is important not to take the ABS value of YS/YY!
         ENDDO
C        ITER=0
C        GOTO 10
      ENDIF
C
C  Is it better to reverse individual components? No!
C
!     NREV=0
!     DO I=1,N
!        IF (W(I)*G(I).GT.0.0D0) THEN
!           W(ISPT+POINT*N+I)= -W(I)
!           NREV=NREV+1
!        ENDIF
!     ENDDO
!     IF (PTEST.AND.(NREV.GT.0)) PRINT '(A,I6,A)',' mylbfgs> ',NREV,' search direction components were reversed'

      IF (CHRMMT.AND.INTMINT) THEN
         DO I=1,N
            W(I)=GINT(I)
         ENDDO
      ELSE
         DO I=1,N
            W(I)=G(I)
         ENDDO
      ENDIF
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
      IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH
      NCOUNT = 0
 19   CONTINUE

      IF (CHRMMT.AND.INTMINT) THEN
         GNORM= DSQRT(DDOT(N,GINT,1,GINT,1))
         DO J1=1,N
            XINT(J1)=XINT(J1)+STP*W(ISPT+POINT*N+J1)
            DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)
         ENDDO 
      ELSE
         GNORM= DSQRT(DDOT(N,G,1,G,1))
!
! Save X here so that we can undo the step reliably.
!
         XSAVE(1:N)=X(1:N)
         DO J1=1,N
            X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
         ENDDO 
      ENDIF
!     PRINT '(A,L5)',' mylbfgs> X,grad,proposed step: PROJECT=',PROJECT
!     DO J1=1,N
!        PRINT '(3F20.10)', X(J1),G(J1),STP*W(ISPT+POINT*N+J1)
!     ENDDO
C
C  To play with step sizes for different sorts of internal coordinates
C  see COPTIM.2.3. Would use MAKESTPVEC.
C
      KNOWE=.FALSE.
      KNOWG=.FALSE.
      KNOWH=.FALSE.
C
C At this point we have new Cartesian or internal coordinates after taking a full
C or decreased step. The gradient is not known at this geometry.
C If INTMIN must transform to Cartesians here.
C

      NDECREASE=0

20    IF (INTMINT) THEN
         IF (CHRMMT) THEN

            NEWQ(1:N)=OLDQ(1:N)
            CART(1:3*NATOMS)=OLDCART(1:3*NATOMS)
C
C Need to keep OLDQ constant for repeated back-transformations if first step size fails.
C Therefore pass dummy array newq that can change.
C Similarly with CART and OLDCART.
C
            CALL TRANSBACKDELTA(DELTAQ,DELTACART,CART,N,3*NATOMS,NNZ,KD,FAILED,.FALSE.,INTEPSILON) ! transform step to Cartesians
            IF (FAILED) THEN
              NCOUNT=NCOUNT+1
              IF (NCOUNT.GT.10) THEN
                 print*, 'Failed to successfully transform
     $                from internal to cartesian step'
                 STOP
              ENDIF

              DO J1=1,N
                 XINT(J1)=XINT(J1)-STP*W(ISPT+POINT*N+J1)
              ENDDO
              STP = STP*0.1
              print*, ' mylbfgs>> backtransform failed; decreasing step', STP
              GOTO 19

            ENDIF
C
C If this is a BFGSTST or MORPHT run project out the uphill direction in Cartesians.
C X is updated from the OLDCARTS value with the projected step. However,
C XINT contains the internal coordinates without projection.
C
            IF ((PROJECT).AND.(.NOT.TWOENDS)) THEN
               DO J2=1,NUP
                  DUMMY1=0.0D0
                  DO J1=1,NOPT
                     DUMMY1=DUMMY1+ZWORK(J1,J2)*DELTACART(J1)
                  ENDDO
                  DO J1=1,NOPT
                     DELTACART(J1)=DELTACART(J1)-ZWORK(J1,J2)*DUMMY1
                  ENDDO
               ENDDO
            ENDIF
C
C now add DELTACART to CART to get new cartesians. Put these in X.
C
            CART(1:3*NATOMS)=OLDCART(1:3*NATOMS)+DELTACART(1:3*NATOMS)
            X(1:3*NATOMS)=OLDCART(1:3*NATOMS)+DELTACART(1:3*NATOMS)
C
C  for CHRMMT:
C  CART    contains new Cartesians (after step) 
C  X       contains new Cartesians (after step)
C  XINT    contains new internals (after step)
C  G       contains old gradient
C  GINT    contains old gradient in internals
C  OLDQ    contains old internals
C  OLDGINT contains old gradient in internals for the last successful geometry
C  NEWQ    contains old internals for the last successful geometry
C  OLDCART contains old Cartesians for the last successful geometry
C
         ELSEIF (UNRST) THEN
            NEWQ(1:N)=X(1:N) ! store new internals in NEWQ
C
C need a temporary array NEWQ here as argument to var_to_geom to keep X unchanged in case we need to
C modify the step below.
C
            CALL var_to_geom(N,NEWQ) ! update internals
            CALL chainbuild ! get cartesians
         ENDIF
      ENDIF

      IF (PV) THEN
         CALL POTENTIAL(X,ENEW,GSAVE,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         PVFLAG=.FALSE.
         CALL PVOPT(X,ENEW,GSAVE)
      ENDIF
      IF(AMBERT.OR.NABT) irespa2=ITDONE+1
      CALL POTENTIAL(X,ENEW,GSAVE,.TRUE.,GRADSQ,RMS,.FALSE.,.FALSE.)

      IF (ENEW.LT.coldFusionLimit) then
         WRITE(*,'(A,2G20.10)') ' mylbfgs> Cold fusion diagnosed - step discarded, energy, limit=',ENEW,coldFusionLimit
         ENERGY=1.0d60
         EREAL=1.0D60
         RMS=1.0d1
         CFUSIONT=.TRUE.
         RETURN
      ENDIF

C     WRITE(*,'(A3,6F20.10)') ('LA ',X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3),
C    1                     GSAVE(3*(J1-1)+1),GSAVE(3*(J1-1)+2),GSAVE(3*(J1-1)+3),J1=1,N/3)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     FIXSAVE=FIXIMAGE
C     FIXIMAGE=.TRUE.
C     SHIT=1.0D-6
C     DO J1=1,3*NATOMS
C        VEC2(J1)=0.0D0
C     ENDDO
C     CALL DSYMV('U',3*NATOMS,1.0D0,HESS,3*NATOMS,GSAVE,1,0.0D0,VEC2,1)
C     DO J1=1,3*NATOMS
C        VEC2(J1)=2*VEC2(J1)
C     ENDDO
C     DO J1=150,160
C        DUMMY1=X(J1)
C        X(J1)=X(J1)+SHIT
C        CALL POTENTIAL(X,EPLUS,GDUM,.TRUE.,.FALSE.,RMSDUM,.FALSE.,.FALSE.)
C        EPLUS=DDOT(3*NATOMS,GDUM,1,GDUM,1)
C        X(J1)=X(J1)-2.0D0*SHIT
C        CALL POTENTIAL(X,EMINUS,GDUM,.TRUE.,.FALSE.,RMSDUM,.FALSE.,.FALSE.)
C        EMINUS=DDOT(3*NATOMS,GDUM,1,GDUM,1)
C        X(J1)=DUMMY1
C        IF (100.0D0*ABS((VEC2(J1)-(EPLUS-EMINUS)/(2.0D0*SHIT))/VEC2(J1)).GT.1.0D0) 
C    1     WRITE(*,'(A,I5,5G20.10)') 'J1,anal,num,%,+,-=',J1,VEC2(J1),(EPLUS-EMINUS)/(2.0D0*SHIT),
C    2                  100.0D0*ABS((VEC2(J1)-(EPLUS-EMINUS)/(2.0D0*SHIT))/VEC2(J1)),EPLUS,EMINUS
C     ENDDO
C     FIXIMAGE=FIXSAVE
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      G(1:N)=GSAVE(1:N)

      IF (TWOENDS.AND.(.NOT.TTDONE)) THEN
         IF (CHRMMT.AND.INTMINT) THEN
            PRINT*,'ERROR - cant use twoends with internals'
            STOP
         ENDIF
         IF (FORCE.NE.0.0D0) THEN
            DUMMY1=0.0D0
            DO J1=1,N
               RVEC(J1)=FIN(J1)-X(J1)
               DUMMY1=DUMMY1+RVEC(J1)**2
            ENDDO
            DUMMY1=1.0D0/SQRT(DUMMY1)
            IF (1.0D0/DUMMY1.GT.0.05D0) THEN
               RMS=0.0D0
                  DO J1=1,N
                  G(J1)=G(J1)+FORCE*RVEC(J1)*DUMMY1
                  RMS=RMS+G(J1)**2
               ENDDO
               RMS=SQRT(RMS/N)
            ENDIF
         ENDIF
      ELSE IF (DRAGT) THEN
         IF (CHRMMT.AND.INTMINT) THEN
            PRINT*,'ERROR - cant use dragt with internals'
            STOP
         ENDIF
         DUMMY1=0.0D0
         DO J1=1,N
            RVEC(J1)=FIN(J1)-X(J1)
            DUMMY1=DUMMY1+RVEC(J1)**2
         ENDDO
         DUMMY1=1.0D0/SQRT(DUMMY1)
         WRITE(*,'(A,F20.10)') 'DIST=',1.0D0/DUMMY1
         IF (1.0D0/DUMMY1.GT.0.3D0) THEN
            DO J1=1,N
               RVEC(J1)=RVEC(J1)*DUMMY1
            ENDDO
            DUMMY1=0.0D0
            DO J1=1,N
               DUMMY1=DUMMY1+RVEC(J1)*G(J1)
            ENDDO
            PRINT*,'Projection of gradient=',DUMMY1
            RMS=0.0D0
            DO J1=1,N
               G(J1)=G(J1)-DUMMY1*RVEC(J1)-ALPHA*RVEC(J1)
               RMS=RMS+G(J1)**2
            ENDDO
            RMS=SQRT(RMS/N)
         ENDIF
      ELSE IF (PROJECT) THEN
         DO J2=1,NUP
            DUMMY1=0.0D0
            IF (UNRST) THEN
               DO J1=1,N
                  DUMMY1=DUMMY1+ZWORK(J1,J2)*G(J1)
               ENDDO
               RMS=0.0D0
               DO J1=1,N
                  G(J1)=G(J1)-ZWORK(J1,J2)*DUMMY1
                  RMS=RMS+G(J1)**2
               ENDDO
               RMS=SQRT(RMS/N)
            ELSE
               DO J1=1,N
                  DUMMY1=DUMMY1+ZWORK(J1,J2)*G(J1)
               ENDDO
               RMS=0.0D0
               DO J1=1,N
                  G(J1)=G(J1)-ZWORK(J1,J2)*DUMMY1
                  RMS=RMS+G(J1)**2
               ENDDO
               RMS=SQRT(RMS/N)
            ENDIF
         ENDDO
      ENDIF
C
C  We need to transform the newly obtained Cartesian gradient for CHARMM and internals. 
C  NOCOOR is true because we dont need to transform the coordinates.
C
      IF (CHRMMT.AND.INTMINT) THEN
         NOCOOR=.TRUE.; NODERV = .FALSE.
         CALL TRANSFORM(X,G,XINT,GINT,N,3*NATOMS,NNZ,NOCOOR,NODERV,KD,INTEPSILON)
      ENDIF

      IF (GRADSQ) THEN
         EREAL=ENEW
         REALRMS=RMS
         ENEW=DDOT(3*NATOMS,G,1,G,1)
C        ENEW=SQRT(DDOT(3*NATOMS,G,1,G,1))
         CALL DSYMV('U',3*NATOMS,1.0D0,HESS,SIZE(HESS,1),G,1,0.0D0,VEC2,1)
         DO J1=1,3*NATOMS
            G(J1)=2*VEC2(J1)
C           G(J1)=VEC2(J1)/ENERGY
         ENDDO
         RMS=DSQRT(DDOT(3*NATOMS,G,1,G,1)/(3*NATOMS))
C
C  NINFO was designed to find GRADSQ discontinuities. Not set at present - we can
C  deal with this using the "Too many failures" route below.
C
         IF (NINFO.GT.100) THEN
            OPEN(UNIT=96,FILE='disconn',STATUS='UNKNOWN')
            PRINT*,' mylbfgs> intractable discontinuity - quit '
            WRITE(96,'(A)') ' mylbfgs> intractable discontinuity'
            CLOSE(96)
            WRITE(*,'(A,4F20.10)')
     $           ' mylbfgs>  g^2, RMS force and real energy and RMS=',ENERGY,RMS,EREAL,REALRMS
            CALL DUMPIT(X,'points.final')
            STOP
         ENDIF
      ELSE
         EREAL=ENEW
         REALRMS=RMS
      ENDIF

      IF ((ENEW-ENERGY.LE.MAXERISE).AND.DRAGT.AND.(ITDONE.GT.NSTEPMIN)) THEN
         WRITE(*,'(A,F20.10,A,F20.10,A)')
     $        ' mylbfgs> Energy falls from ',ENERGY,' to ',ENEW,
     $        ' try ts search from previous geometry'
         DO J1=1,N
            X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
         ENDDO 
         KNOWE=.FALSE.
         KNOWG=.FALSE.
         KNOWH=.FALSE.
C        DGUESS=DIAG(1) ! saved for subsequent calls - should be OK for the same system?
C                          ! May make redopath runs unreproducible?
         RETURN
      ENDIF
C
C  Must allow the energy to rise during a minimisation to allow for numerical noise or
C  systematic errors due to discontinuities or SCF convergence problems.
C
      IF ((ENEW-ENERGY.LE.MAXERISE).OR.PVTS.OR.DRAGT.OR.TWOENDS.OR.RADMOVED) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A,G13.5)')
     $        ' mylbfgs> Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,
     1           ' steps, step:',STP*SLENGTH
         WRITE(ESTRING,16) ' mylbfgs> Energy for last cycle=',ENERGY,' '
C
C  Step finished so can reset OLDQ to new XINT, OLDCART to new CART,
C  as well as the Cartesian and internal gradients.
C
         IF (CHRMMT.AND.INTMINT) THEN
            OLDGINT(1:N)=GINT(1:N)
            OLDCART(1:3*NATOMS)=CART(1:3*NATOMS)
C
C  Need to remake XINT because step was only projected in Cartesians?
C  Actually, just setting OLDQ=XINT without this correction seems to
C  be OK. Due to numerical imprecision, it might still be possible
C  for X and XINT to get out of register. Perhaps this doesn't matter
C  because the energy and gradient are always calculated in Cartesians.
C
C           IF (PROJECT) CALL TRANSDELTA(DELTACART,DELTAQ,CART,N,3*NATOMS,NNZ,KD,INTEPSILON)
C           OLDQ(1:N)=OLDQ(1:N)+DELTAQ(1:N)
            OLDQ(1:N)=XINT(1:N)
         ELSEIF (UNRST) THEN
!           TEST1(1:N)=X(1:N)
            CALL geom_to_var(N,X(1:N)) ! testing!!! - to put X back into register with the common block internals (and g)
!           CALL geom_to_var(N,TEST1(1:N))
!           do j1=1,N
!           if (abs((TEST1(j1)-x(j1))/x(j1))*100.0d0.gt.1.0D-6) print *,'hello coords ',J1
!           enddo
         ENDIF
         GLAST(1:N)=GSAVE(1:N) 
         IF (.NOT.(NODUMP)) THEN ! jmc dumpp dumps x but x is in internals...
            IF (INTMINT) THEN
               IF (UNRST) THEN
!                 TMPINT(1:N)=X(1:N) ! the Cartesians should not have changed from when they were set before the 
                                     ! call to potential above, so commented out these three lines
!                 CALL var_to_geom(N,TMPINT)
!                 CALL chainbuild
                  DO I1=1,nres
                     CART(6*(I1-1)+1)=c(1,I1)
                     CART(6*(I1-1)+2)=c(2,I1)
                     CART(6*(I1-1)+3)=c(3,I1)
                     CART(6*(I1-1)+4)=c(1,I1+nres)
                     CART(6*(I1-1)+5)=c(2,I1+nres)
                     CART(6*(I1-1)+6)=c(3,I1+nres)
                  ENDDO
                  CALL DUMPP(CART,ENERGY)
               ELSEIF (CHRMMT) THEN
                  CALL DUMPP(X,ENERGY)
               ENDIF
            ELSE
               CALL DUMPP(X,ENERGY)
            ENDIF
         ENDIF
      ELSE 
C
C  Energy increased - try again with a smaller step size. Must cater for possible enormous
C  values of SLENGTH. Decreasing the step size doesn;t seem to help for CASTEP.
C
         IF (((ITER.GT.1).AND.(NDECREASE.GT.10)).OR.((ITER.LE.1).AND.(NDECREASE.GT.10)).OR.
     1              ((CASTEP.OR.ONETEP.OR.CP2K).AND.(NDECREASE.GT.10))) THEN 
            NFAIL=NFAIL+1
            IF (PTEST) WRITE(*,'(2(A,I6))') ' mylbfgs> in mylbfgs step ',ITER,
     $           ' cannot find a lower energy, NFAIL=',NFAIL
C
C  try resetting - go back to previous coordinates, ENERGY is not set to ENEW
C  we need to save the gradient corresponding to the last successful step
C              
            ITER=0
            IF (CHRMMT.AND.INTMINT) THEN ! need to reset X, XINT, G, GINT to original values 
               XINT(1:N)=XINT(1:N)-STP*W(ISPT+POINT*N+1:ISPT+POINT*N+N)
C              XINT=OLDQ ! should be the same as subtracting the step
               GINT(1:N)=OLDGINT(1:N)
               G(1:3*NATOMS)=GLAST(1:3*NATOMS)
               X(1:3*NATOMS)=OLDCART(1:3*NATOMS)
            ELSE
!
! Resetting to XSAVE should be the same as subtracting the step.
!
               X(1:N)=XSAVE(1:N)
               DO J1=1,N
C                 X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
                  G(J1)=GLAST(J1)
C                 G(J1)=GSAVE(J1) ! DJW 6/5/04
               ENDDO
            ENDIF
            IF (NFAIL.GT.NFAILMAX) THEN
               WRITE(*,'(A)') ' mylbfgs> Too many failures - give up'
               CALL DUMPCOORDS(X, 'lbfgsfailed.xyz', .FALSE.)
               FIXIMAGE=.FALSE.
C              DGUESS=DIAG(1) ! saved for subsequent calls - should be OK for the same system?
C                          ! May make redopath runs unreproducible?
               RETURN
            ENDIF
            GOTO 30
         ENDIF
C
C  Try a smaller step.
C
         IF (CHRMMT.AND.INTMINT) THEN
            DO J1=1,N
               XINT(J1)=XINT(J1)-0.9*STP*W(ISPT+POINT*N+J1)
               DELTAQ(J1)=STP*W(ISPT+POINT*N+J1)*0.1D0
            ENDDO 
         ELSE
!
! Resetting to XSAVE and adding 0.1 of the step should be the same as subtracting
! 0.9 of the step.
!
!        PRINT*,'X should match XSAVE:'
!        WRITE(*,'(I5,2E20.10)') (J1,X(J1)-STP*W(ISPT+POINT*N+J1),XSAVE(J1),J1=1,N)
            DO J1=1,N
               X(J1)=X(J1)-0.9D0*STP*W(ISPT+POINT*N+J1)
            ENDDO 
         ENDIF
         KNOWE=.FALSE.
         KNOWG=.FALSE.
         KNOWH=.FALSE.
         STP=STP/10.0D0
         NDECREASE=NDECREASE+1
         IF (PTEST) 
     1    WRITE(*,'(A,G25.15,A,G25.15,A,G15.8)') ' energy increase:',ENERGY,' to ',ENEW,
     2            ' decreasing step to ',STP*SLENGTH
!        PRINT*,'X:'
!        WRITE(*,'(I5,2E20.10)') (J1,X(J1),XSAVE(J1)+STP*W(ISPT+POINT*N+J1),J1=1,N)
C        FIXIMAGE=.TRUE. ! BLJ seems to work better without this
         GOTO 20
      ENDIF
C
C     Compute the new step and gradient change. Note that the step
C     length is accounted for when the step taken is saved.
C
30    NPT=POINT*N

      IF (CHRMMT.AND.INTMINT) THEN
         DO I=1,N
            W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
            W(IYPT+NPT+I)= GINT(I)-W(I)
         ENDDO
      ELSE
         DO I=1,N
            W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
            W(IYPT+NPT+I)= G(I)-W(I)
         ENDDO
      ENDIF
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      IF (.NOT.(NODUMP)) THEN
         IF ((AMBER).AND.(MOVIE)) CALL AMOVIEDUMP(FRAME)
      ENDIF
      FIXIMAGE=.FALSE.
      IF ((FIXAFTER.GT.0).AND.(ITER.GE.FIXAFTER)) FIXIMAGE=.TRUE.
      GOTO 10

      RETURN
      END
