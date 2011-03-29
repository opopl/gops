!   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
!   COPYRIGHT (C) 1999-2006 DAVID J. WALES
!   THIS FILE IS PART OF OPTIM.
!   
!   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
!   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
!   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
!   (AT YOUR OPTION) ANY LATER VERSION.
!   
!   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
!   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
!   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
!   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
!   
!   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
!   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
!   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
!
SUBROUTINE UNMECCANO(ITEST,PTEST,EMAX,PERM,Q,FIN,E1,E2,RMSINITIAL,RMSFINAL)
USE COMMONS
USE KEY
USE MODNEB
USE MODUNRES
USE MODMEC
USE KEYNEB, ONLY: NITERMAX
IMPLICIT NONE
INTEGER :: J1, J2, J3, ITDONE, JMAX
DOUBLE PRECISION :: RMS, EMAX, TIME, TIME0, Q(3*NATOMS), FIN(3*NATOMS), &
     &              RMSVEC(NIMAGE), EVEC(NIMAGE), E1, E2, &
     &              POINTS(3*NATOMS*NIMAGE+1), SEPARATION, EINITIAL, EFINAL, RMSINITIAL, RMSFINAL, &
     &              DELTA(3*NATOMS), &
     &              GRAD(3*NATOMS), &
     &              EFORIG, EIORIG, FINORIG(3*NATOMS)
DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793D0
DOUBLE PRECISION :: PEPCOORDS(3*NATOMS)
LOGICAL :: MFLAG, ITEST, PTEST, PERM
COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
CHARACTER FNAME*17
INTEGER IADD, LNOPT, LNATOMS
COMMON /INTS/ IADD, LNOPT, LNATOMS
CHARACTER(LEN=80) LFNAME

IF (.NOT.ITEST) THEN
   EINITIAL=E1
   EFINAL=E2
ELSE
!CALL VAR_TO_GEOM(NINTS,Q)
!CALL CHAINBUILD
   CALL POTENTIAL(Q,EINITIAL,GRAD,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
!CALL VAR_TO_GEOM(NINTS,FIN)
!CALL CHAINBUILD
   CALL POTENTIAL(FIN,EFINAL,GRAD,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
ENDIF
IF (ITEST.AND.PTEST) WRITE(*,'(A,F20.10,A,F15.6,A,F20.10,A,F15.6)') &
     &           ' INITIAL ENERGY=',EINITIAL,' RMS FORCE=',RMSINITIAL,' FINAL ENERGY=',EFINAL,' RMS=',RMSFINAL


IADD=0 ; LNOPT=NINTS ; LNATOMS=NATOMS

CALL MYCPU_TIME(TIME0,.TRUE.)
EIORIG=EINITIAL ; EFORIG=EFINAL
FINORIG(1:LNOPT)=FIN(1:LNOPT)
DO J1=1,LNOPT
   IF (ABS(FIN(J1)-Q(J1)).LE.PI) THEN
      DELTA(J1)=FIN(J1)-Q(J1) ! BETWEEN 0 AND PI
   ELSE
      IF (FIN(J1).GT.Q(J1)) THEN
         DELTA(J1)=FIN(J1)-Q(J1)-2*PI ! BETWEEN 0 AND -PI
      ELSE
         DELTA(J1)=FIN(J1)-Q(J1)+2*PI ! BETWEEN 0 AND PI
      ENDIF
   ENDIF
ENDDO
DO J1=1,NIMAGE
   DO J2=1,LNOPT
      POINTS(LNOPT*(J1-1)+J2)=Q(J2)+DELTA(J2)*J1/(NIMAGE+1)
   ENDDO
!CALL VAR_TO_GEOM(LNOPT,POINTS(LNOPT*(J1-1)+1:LNOPT*(J1-1)+LNOPT))
!CALL GEOM_TO_VAR(LNOPT,POINTS(LNOPT*(J1-1)+1:LNOPT*(J1-1)+LNOPT))
ENDDO
POINTS(LNOPT*NIMAGE+1)=MECDIST
!
!  MINIMISATION OF LAGRANGIAN
!
IF (FILTH.EQ.0) THEN
   WRITE(FNAME,'(A8)') 'EOFS.NEB'
ELSE 
   WRITE(FNAME,'(A)') 'EOFS.NEB.'//TRIM(ADJUSTL(FILTHSTR))
ENDIF
OPEN(UNIT=4,FILE=FNAME,STATUS='UNKNOWN')
CALL UNMECBFGS(LNOPT*NIMAGE+1,MECUPDATE,POINTS,.FALSE.,MECRMSTOL,MFLAG,EVEC,RMSVEC, &
     &             NITERMAX,ITDONE,PTEST,Q,FIN)

         WRITE(4,'(3F20.10)') 0.0D0,EINITIAL,RMSINITIAL
         WRITE(4,'(2F20.10)') (EVEC(J1),RMSVEC(J1),J1=1,NIMAGE)
         WRITE(4,'(2F20.10)') EFINAL,RMSFINAL
      CLOSE(4)
!
!  DUMP PATHWAY AS AN XYZ FILE
!
! IF (FILTH.EQ.0) THEN
   ! WRITE(FNAME,'(A12)') 'NEB.PATH.XYZ'
! ELSE 
   ! WRITE(FNAME,'(A)') 'NEB.PATH.XYZ.'//TRIM(ADJUSTL(FILTHSTR))
! ENDIF
! OPEN(UNIT=3,FILE=FNAME,STATUS='UNKNOWN')
! WRITE(3,'(I6)') LNATOMS
! WRITE(3,'(F20.10)') EINITIAL
! WRITE(3,'(A2,4X,3F20.10)') (ZSYM(J1),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),J1=1,LNATOMS)
! DO J1=1,NIMAGE
   ! WRITE(3,'(I6)') LNATOMS
   ! WRITE(3,'(F20.10)') EVEC(J1) 
   ! WRITE(3,'(A2,4X,3F20.10)') &
   ! &  (ZSYM(J2),POINTS(LNOPT*(J1-1)+3*(J2-1)+1),POINTS(LNOPT*(J1-1)+3*(J2-1)+2),POINTS(LNOPT*(J1-1)+3*(J2-1)+3),J2=1,LNATOMS)
! ENDDO
! WRITE(3,'(I6)') LNATOMS
! WRITE(3,'(F20.10)') EFINAL
! WRITE(3,'(A2,4X,3F20.10)') (ZSYM(J1),FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3),J1=1,LNATOMS)
! CLOSE(3)
!
!  FIND THE HIGHEST ENERGY IMAGE.
!
JMAX=1
EMAX=-1.0D100
DO J1=1,NIMAGE
   IF (EVEC(J1).GT.EMAX) THEN
      EMAX=EVEC(J1)
      JMAX=J1
   ENDIF
ENDDO
CALL MYCPU_TIME(TIME,.FALSE.)
WRITE(*,'(A,I4,A,G15.5,A,I6,A,F12.4,A,I4,A,F11.2)') ' IMAGE ',JMAX,' HAS HIGHEST ENERGY=',EMAX, &
     &                   ' STEPS=',ITDONE,' RMS=',RMS,' IMAGES=',NIMAGE,'    TIME=',TIME-TIME0
TIME0=TIME
!
!  WRITE COORDINATES TO GUESS.XYZ AS A CANDIDATE INTERPOLATION
!
IF (FILTH.EQ.0) THEN
   OPEN(UNIT=7,FILE='GUESS.XYZ',STATUS='UNKNOWN')
ELSE
   LFNAME='GUESS.XYZ.'//TRIM(ADJUSTL(FILTHSTR))
   OPEN(UNIT=7,FILE=TRIM(ADJUSTL(LFNAME)),STATUS='UNKNOWN')
ENDIF
IF (FILTH.EQ.0) THEN
   OPEN(UNIT=8,FILE='GUESS.UNRES.XYZ',STATUS='UNKNOWN')
ELSE
   LFNAME='GUESS.UNRES.XYZ.'//TRIM(ADJUSTL(FILTHSTR))
   OPEN(UNIT=8,FILE=TRIM(ADJUSTL(LFNAME)),STATUS='UNKNOWN')
ENDIF

DO J1=0,NIMAGE+1
   IF (J1.EQ.0) THEN 
!CALL VAR_TO_GEOM(LNOPT,Q(1:LNOPT))
   ELSE IF (J1.EQ.NIMAGE+1) THEN
!CALL VAR_TO_GEOM(LNOPT,FIN(1:LNOPT))
   ELSE
!CALL VAR_TO_GEOM(LNOPT,POINTS(LNOPT*(J1-1)+1:LNOPT*J1))
   ENDIF
!CALL CHAINBUILD
   DO J2=1,NRES
      WRITE(7,'(3G20.10)') C(1,J2),C(2,J2),C(3,J2) ! BACKBONE
      WRITE(7,'(3G20.10)') C(1,J2+NRES),C(2,J2+NRES),C(3,J2+NRES) ! SIDE CHAINS
   ENDDO
   DO J2=1,(NATOMS/2)-1 ! JMC ADD PEPTIDE ATOMS...
      DO J3=1,3
         PEPCOORDS(6*(J2-1)+J3)=(2.0D0*C(J3,J2)+C(J3,J2+1))/3.0D0
         PEPCOORDS(6*(J2-1)+J3+3)=(C(J3,J2)+2.0D0*C(J3,J2+1))/3.0D0
      ENDDO
   ENDDO
   WRITE(8,'(I6)') 2*NATOMS-2
   WRITE(8,'(A)') ' '
   WRITE(8,'(A2,3G20.10)') ('C  ',C(1,J2),C(2,J2),C(3,J2),J2=1,NATOMS)
   WRITE(8,'(A2,4X,3F20.10)') ('O ',PEPCOORDS(6*(J2-1)+1),PEPCOORDS(6*(J2-1)+2),PEPCOORDS(6*(J2-1)+3) &
                             ,J2=1,(NATOMS/2)-1)
   WRITE(8,'(A2,4X,3F20.10)') ('N ',PEPCOORDS(6*(J2-1)+4),PEPCOORDS(6*(J2-1)+5),PEPCOORDS(6*(J2-1)+6) &
                             ,J2=1,(NATOMS/2)-1)
ENDDO
CLOSE(7) ; CLOSE(8) 
! OPEN(UNIT=7,FILE='INTERNALS',STATUS='UNKNOWN')
! WRITE(7,'(3G20.10)') POINTS(1:LNOPT*NIMAGE)
! CLOSE(7)

RETURN
END SUBROUTINE UNMECCANO
!
!        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!                          JORGE NOCEDAL
!                        *** JULY 1990 ***
!
!        LINE SEARCH REMOVED PLUS SMALL MODIFICATIONS, DJW 2001
!        IF INTMINT IS TRUE, THEN N IS NINTS, 
!        OTHERWISE N = NOPT JMC
!        CHANGED DECLARATION OF X(N) TO X(3*NATOMS) 30/4/04
!        X IS PASSED IN AND OUT IN CARTESIANS.
!
!  CARTESIAN COORDINATE AND GRADIENT VECTORS ARE DECLARED 3*NATOMS - THE COMPLICATION
!  IS THAT FOR INTERNAL COORDINATE OPTIMISATIONS WE CAN HAVE A NUMBER OF DEGREES OF
!  FREEDOM THAT IS MORE OR LESS THAN 3*NATOMS. N SHOULD SPECIFY THIS DIMENSION.
!
SUBROUTINE UNMECBFGS(N,M,X,DIAGCO,EPS,MFLAG,EVEC,RMSVEC,ITMAX,ITDONE,PTEST, &
     &    QSTART,QFINISH)
USE COMMONS
USE KEY
USE MODTWOEND
USE ZWK
USE MODUNRES
USE MODCHARMM
USE MODNEB
USE MODMEC
USE PORFUNCS
IMPLICIT NONE
INTEGER N,M,J1,ITMAX,ITDONE,NFAIL,NCOUNT
DOUBLE PRECISION X(N),G(N),SLENGTH,DDOT,OVERLAP,QSTART(*),QFINISH(*),XLAST(N)
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W
DOUBLE PRECISION EPS,DUMMY1,ENERGY,ENEW,RMS,EREAL,REALRMS,ALPHA,GSAVE(N)
LOGICAL DIAGCO, PTEST
DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA
DOUBLE PRECISION EPSILON,EVEC(NIMAGE),RMSVEC(NIMAGE)
DOUBLE PRECISION GINT(N) ! JMC
DOUBLE PRECISION GLAST(N)
INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN
LOGICAL MFLAG
DOUBLE PRECISION DOT1, DOT2
INTEGER FRAME
DOUBLE PRECISION :: EINITIAL, EFINAL, SEPARATION, EPLUS, EMINUS
COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
LOGICAL PATHT, DRAGT
INTEGER NPATHFRAME, NDECREASE, ISTAT
COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
LOGICAL KNOWE, KNOWG, KNOWH
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
CHARACTER ESTRING*87, GPSTRING*80, NSTRING*80, FSTRING*80
COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING
DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793D0

ALLOCATE(DIAG(N)) 
ALLOCATE(W(N*(2*M+1)+2*M)) 

ALPHA=1.0D0
FRAME=1
NFAIL=0
ITER=0
ITDONE=0

FIXIMAGE=.FALSE.

! PUT COORDINATES INTO THE +PI TO -PI RANGE - SHOULD HAVE BEEN DONE IN CALLING ROUTINE!
DO J1=1,NIMAGE
!CALL VAR_TO_GEOM(NINTS,X(NINTS*(J1-1)+1:NINTS*(J1-1)+NINTS))
!CALL GEOM_TO_VAR(NINTS,X(NINTS*(J1-1)+1:NINTS*(J1-1)+NINTS))
ENDDO
XLAST(1:N)=X(1:N)

CALL UNMAKEGRADMEC(X,GSAVE,ENERGY,EVEC,RMSVEC,QSTART,QFINISH)
G(1:N)=GSAVE(1:N)
GLAST(1:N)=GSAVE(1:N)

! CHECK NUMERICAL DERIVATIVES
IF (.FALSE.) THEN
   DO J1=1,N
      X(J1)=X(J1)+1.0D-5
      CALL UNMAKEGRADMEC(X,G,EPLUS,EVEC,RMSVEC,QSTART,QFINISH)
      X(J1)=X(J1)-2.0D-5
      CALL UNMAKEGRADMEC(X,G,EMINUS,EVEC,RMSVEC,QSTART,QFINISH)
      X(J1)=X(J1)+1.0D-5
      IF (ABS(100*(GSAVE(J1)-(EPLUS-EMINUS)/(2.0D-5))/GSAVE(J1)).GT.1.0D0) &
 &    PRINT '(A,I5,4G20.10)', 'J1,X,ANAL,NUM,%ERROR=', &
       &   J1,X(J1),GSAVE(J1),(EPLUS-EMINUS)/(2.0D-5),ABS(100*(GSAVE(J1)-(EPLUS-EMINUS)/(2.0D-5))/GSAVE(J1))
   ENDDO
   STOP
ENDIF

EREAL=ENERGY
REALRMS=RMS

IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A)') &
     &             ' ENERGY AND RMS FORCE=',ENERGY,RMS,' AFTER ',ITDONE,' LBFGS STEPS'
      WRITE(ESTRING,16) ' ENERGY FOR LAST CYCLE=',ENERGY,' '
16 FORMAT(A,27X,F20.10,A)

10 CALL FLUSH(6,ISTAT)
DO J1=1,NIMAGE
!CALL VAR_TO_GEOM(NINTS,X(NINTS*(J1-1)+1:NINTS*(J1-1)+NINTS))
!CALL GEOM_TO_VAR(NINTS,X(NINTS*(J1-1)+1:NINTS*(J1-1)+NINTS))
ENDDO
XLAST(1:N)=X(1:N)

MFLAG=.FALSE.
IF (RMS.LE.EPS) THEN
   MFLAG=.TRUE.
!  PRINT*,'RMS,EPS,ITDONE,NSTEPMIN=',RMS,EPS,ITDONE,NSTEPMIN
   IF (ITDONE.LT.NSTEPMIN) MFLAG=.FALSE.
   IF (MFLAG) THEN
      FIXIMAGE=.FALSE.
      DEALLOCATE(DIAG,W)
      RETURN
   ENDIF
ENDIF

IF (ITDONE.EQ.ITMAX) THEN
   FIXIMAGE=.FALSE.
!  WRITE(*,'(A,F20.10)') ' DIAGONAL INVERSE HESSIAN ELEMENTS ARE NOW ',DIAG(1)
   MECDGUESS=DIAG(1) ! SAVED FOR SUBSEQUENT CALLS - SHOULD BE OK FOR THE SAME SYSTEM?
   DEALLOCATE(DIAG,W)
   RETURN
ENDIF

IF (ITER.EQ.0) THEN
   IF (N.LE.0.OR.M.LE.0) THEN
      WRITE(*,240)
240  FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
      PRINT*,'N,M=',N,M
      STOP
   ENDIF
   POINT=0
   MFLAG=.FALSE.
   IF (DIAGCO) THEN
      PRINT*,'USING ESTIMATE OF THE INVERSE DIAGONAL ELEMENTS'
      DO I=1,N
         IF (DIAG(I).LE.0.0D0) THEN
            WRITE(*,235) I
235         FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/, &
     &                   ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
            STOP
         ENDIF
      ENDDO
   ELSE
      DO I=1,N
         DIAG(I)=MECDGUESS
      ENDDO
   ENDIF
   ISPT= N+2*M
   IYPT= ISPT+N*M
!
!  NR STEP FOR DIAGONAL INVERSE HESSIAN
!
   DO I=1,N
      W(ISPT+I)= -G(I)*DIAG(I)
      W(I)= -G(I)*DIAG(I)
   ENDDO
!  PRINT*,'I,W,W=',I,W(1),W(ISPT+1)
   GNORM= DSQRT(DDOT(N,G,1,G,1))
!
!  MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS.
!
   STP=MIN(1.0D0/GNORM,GNORM)
!  STP=1.0D0
ELSE 
   BOUND=ITER
   IF (ITER.GT.M) BOUND=M
!  PRINT*,'BEFORE OVERLAP W, W: ITER,M,ISPT,IYPT,NPT=',ITER,M,ISPT,IYPT,NPT
!  WRITE(*,'(I5,2E20.10)') (J1,W(ISPT+NPT+J1),W(IYPT+NPT+J1),J1=1,10)
   YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
!  WRITE(*,'(A,E20.10)') 'YS=',YS
   IF (YS.EQ.0.0D0) YS=1.0D0
!
!  UPDATE ESTIMATE OF DIAGONAL INVERSE HESSIAN ELEMENTS
!  WE DIVIDE BY BOTH YS AND YY AT DIFFERENT POINTS, SO
!  THEY HAD BETTER NOT BE ZERO!
!
   IF (.NOT.DIAGCO) THEN
      YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
!     WRITE(*,'(A,E20.10)') 'YY=',YY
      IF (YY.EQ.0.0D0) YY=1.0D0
      DUMMY1=YS/YY
!     DUMMY1=ABS(YS/YY)
!     WRITE(*,'(A,E20.10)') 'DUMMY1=',DUMMY1
      DO I=1,N
         DIAG(I)=DUMMY1
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
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: NOCEDAL, J. 1980,
!     "UPDATING QUASI-NEWTON MATRICES WITH LIMITED STORAGE",
!     MATHEMATICS OF COMPUTATION, VOL.24, NO.151, PP. 773-782.
!     ---------------------------------------------------------
!
   CP= POINT
   IF (POINT.EQ.0) CP=M
   W(N+CP)= 1.0D0/YS
!  PRINT*,'W(I) GETS SET TO -G(I):'
!  WRITE(*,'(I5,2E20.10)') (J1,W(J1),G(J1),J1=1,10)
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
!
!  STORE THE NEW SEARCH DIRECTION
!
!     PRINT*,'W(I):'
!     WRITE(*,'(I5,E20.10)') (J1,W(J1),J1=1,10)
IF (ITER.GT.0) THEN
   DO I=1,N
      W(ISPT+POINT*N+I)= W(I)
   ENDDO
ENDIF

DOT1=SQRT(DDOT(N,G,1,G,1))
DOT2=SQRT(DDOT(N,W,1,W,1))

OVERLAP=0.0D0
IF (DOT1*DOT2.NE.0.0D0) THEN
   IF (CHRMMT.AND.INTMINT) THEN
      OVERLAP=DDOT(N,GINT,1,W,1)/(DOT1*DOT2)
   ELSE
      OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)
  ENDIF
ENDIF
!     PRINT*,'OVERLAP,DIAG(1)=',OVERLAP,DIAG(1)
!     PRINT*,'G . G=',DDOT(N,G,1,G,1)
!     PRINT*,'W . W=',DDOT(N,W,1,W,1)
IF (OVERLAP.GT.0.0D0) THEN
   IF (PTEST) PRINT*,'SEARCH DIRECTION HAS POSITIVE PROJECTION ONTO GRADIENT - REVERSING STEP'
   DO I=1,N
      W(ISPT+POINT*N+I)= -W(I)
   ENDDO
!        ITER=0
!        GOTO 10
ENDIF

DO I=1,N
   W(I)=G(I)
ENDDO
SLENGTH=0.0D0
DO J1=1,N
   SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
ENDDO
SLENGTH=SQRT(SLENGTH)
IF (STP*SLENGTH.GT.MECSTEP) STP=MECSTEP/SLENGTH
!
!  WE NOW HAVE THE PROPOSED STEP.
!
DO J1=1,N
   X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
!  PRINT '(A,I6,2G20.10)','J1,X,W,STP=',J1,X(J1),W(ISPT+POINT*N+J1)
ENDDO 
KNOWE=.FALSE.
KNOWG=.FALSE.
KNOWH=.FALSE.
!
! AT THIS POINT WE HAVE NEW CARTESIAN OR INTERNAL COORDINATES AFTER TAKING A FULL
! OR DECREASED STEP. THE GRADIENT IS NOT KNOWN AT THIS GEOMETRY.
! IF INTMIN MUST TRANSFORM TO CARTESIANS HERE.
!
NDECREASE=0
EPSILON=1.0D-6
NCOUNT=0
20    CONTINUE
! DO J1=1,NIMAGE ! NOT HERE - DO THIS ONLY AT THE TOP OF THE LOOP
   ! CALL VAR_TO_GEOM(NINTS,X(NINTS*(J1-1)+1:NINTS*(J1-1)+NINTS))
   ! CALL GEOM_TO_VAR(NINTS,X(NINTS*(J1-1)+1:NINTS*(J1-1)+NINTS))
! ENDDO

CALL UNMAKEGRADMEC(X,GSAVE,ENEW,EVEC,RMSVEC,QSTART,QFINISH)

G(1:N)=GSAVE(1:N)

EREAL=ENEW
REALRMS=RMS
!
!  MUST ALLOW THE ENERGY TO RISE DURING A MINIMISATION TO ALLOW FOR NUMERICAL NOISE OR
!  SYSTEMATIC ERRORS DUE TO DISCONTINUITIES OR SCF CONVERGENCE PROBLEMS.
!
! IF (ENEW-ENERGY.LE.1.0D-5) THEN
IF (ENEW-ENERGY.LE.MAXERISE) THEN
   ITER=ITER+1
   ITDONE=ITDONE+1
   ENERGY=ENEW
   IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A,G13.5)') ' ENERGY AND RMS FORCE=',ENERGY,RMS,' AFTER ',ITDONE, &
     &     ' LBFGS STEPS, STEP:',STP*SLENGTH
   WRITE(ESTRING,16) ' ENERGY FOR LAST CYCLE=',ENERGY,' '
!
!  STEP FINISHED SO CAN RESET OLDQ TO NEW XINT, OLDCART TO NEW CART,
!  AS WELL AS THE CARTESIAN AND INTERNAL GRADIENTS.
!
   GLAST(1:N)=GSAVE(1:N)
!  IF (.NOT.BFGSTST) CALL DUMPP(X,ENERGY)
ELSE 
!
!  ENERGY INCREASED - TRY AGAIN WITH A SMALLER STEP SIZE. MUST CATER FOR POSSIBLE ENORMOUS
!  VALUES OF SLENGTH. DECREASING THE STEP SIZE DOESN;T SEEM TO HELP FOR CASTEP.
!
   IF (((ITDONE.GT.1).AND.(NDECREASE.GT.9)).OR.((ITDONE.LE.1).AND.(NDECREASE.GT.9)).OR. &
     &        ((CASTEP.OR.ONETEP.OR.CP2K).AND.(NDECREASE.GT.1))) THEN 
      NFAIL=NFAIL+1
      PRINT*,' IN MECBFGS LBFGS STEP CANNOT FIND A LOWER ENERGY, NFAIL=',NFAIL
      IF (PTEST) PRINT*,' IN MECBFGS LBFGS STEP CANNOT FIND A LOWER ENERGY, NFAIL=',NFAIL
!
!  TRY RESETTING - GO BACK TO PREVIOUS COORDINATES, ENERGY IS NOT SET TO ENEW
!  WE NEED TO SAVE THE GRADIENT CORRESPONDING TO THE LAST SUCCESSFUL STEP
!              
      ITER=0
      DO J1=1,N
         X(J1)=XLAST(J1)
         G(J1)=GLAST(J1)
!        G(J1)=GSAVE(J1) ! DJW 6/5/04
      ENDDO
! WITH ALL THE POSSIBLE DISCONTINUITIES WE PROBABLY NEED TO RESET EVERYTHING

      CALL UNMAKEGRADMEC(X,GSAVE,ENEW,EVEC,RMSVEC,QSTART,QFINISH)
      ENERGY=ENEW
      G(1:N)=GSAVE(1:N)
      PRINT*,'NDECREASE,NFAIL=',NDECREASE,NFAIL
      PRINT*,'ABOUT TO TEST THE VALUE OF NFAIL'
      IF (NFAIL.GT.5) THEN
         PRINT*,' TOO MANY FAILURES - GIVE UP'
! CHECK NUMERICAL DERIVATIVES
!        DO J1=1,N
!           X(J1)=X(J1)+1.0D-5
!           PRINTT=.FALSE.
!           CALL UNMAKEGRADMEC(X,G,EPLUS,EVEC,RMSVEC,QSTART,QFINISH)
!           X(J1)=X(J1)-2.0D-5
!           CALL UNMAKEGRADMEC(X,G,EMINUS,EVEC,RMSVEC,QSTART,QFINISH)
!           X(J1)=X(J1)+1.0D-5
!           IF (ABS(100*(GSAVE(J1)-(EPLUS-EMINUS)/(2.0D-5))/GSAVE(J1)).GT.0.01D0) THEN
!              PRINT '(A,I5,4G20.10)', 'J1,X,ANAL,NUM,%ERROR=', &
!         &      J1,X(J1),GSAVE(J1),(EPLUS-EMINUS)/(2.0D-5),ABS(100*(GSAVE(J1)-(EPLUS-EMINUS)/(2.0D-5))/GSAVE(J1))
!              PRINT '(A,2G20.10)', 'EPLUS,EMINUS=',EPLUS,EMINUS
!              PRINT*,'VALUE OF THIS COORDINATE IN ALL IMAGES:'
!              PRINT '(I5,G20.10)', 0,QSTART(J1-(J1/NINTS)*NINTS)
!              DO J2=1,NIMAGE
!                 PRINT '(I5,G20.10)', J2,X(J1-(J1/NINTS)*NINTS+(J2-1)*NINTS)
!              ENDDO
!              PRINT '(I5,G20.10)', NIMAGE+1,QFINISH(J1-(J1/NINTS)*NINTS)
!           ENDIF
!        ENDDO
!        STOP
         FIXIMAGE=.FALSE.
         MECDGUESS=DIAG(1) ! SAVED FOR SUBSEQUENT CALLS - SHOULD BE OK FOR THE SAME SYSTEM?
         PRINT*,'FIXIMAGE=',FIXIMAGE
         PRINT*,'DIAG(1)=',DIAG(1)
         PRINT*,'MECDGUESS=',MECDGUESS
         DEALLOCATE(DIAG,W)
         RETURN
      ENDIF
      GOTO 30
   ENDIF
!
!  TRY A SMALLER STEP.
!
   DO J1=1,N
      X(J1)=XLAST(J1)-0.9*STP*W(ISPT+POINT*N+J1)
   ENDDO 
   KNOWE=.FALSE.
   KNOWG=.FALSE.
   KNOWH=.FALSE.
   STP=STP/10.0D0
   NDECREASE=NDECREASE+1
   IF (PTEST) &
     &    WRITE(*,'(A,G19.10,A,G16.10,A,G15.8)') ' ENERGY INCREASED FROM ',ENERGY,' TO ',ENEW, &
     &            ' DECREASING STEP TO ',STP*SLENGTH
   FIXIMAGE=.TRUE.
   GOTO 20
ENDIF
!
!     COMPUTE THE NEW STEP AND GRADIENT CHANGE. NOTE THAT THE STEP
!     LENGTH IS ACCOUNTED FOR WHEN THE STEP TAKEN IS SAVED.
!
30    NPT=POINT*N

DO I=1,N
   W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
   W(IYPT+NPT+I)= G(I)-W(I)
ENDDO
POINT=POINT+1
IF (POINT.EQ.M) POINT=0
FIXIMAGE=.FALSE.
IF ((FIXAFTER.GT.0).AND.(ITER.GE.FIXAFTER)) FIXIMAGE=.TRUE.
GOTO 10

RETURN
END SUBROUTINE UNMECBFGS
!
!****************************************************************************************
!
SUBROUTINE UNMAKEGRADMEC(X,G,SUM,EVEC,RMSVEC,QSTART,FIN)
USE COMMONS
USE MODHESS
USE MODMEC
USE KEY
USE MODNEB
USE PORFUNCS
IMPLICIT NONE
INTEGER J1,J2,J3
DOUBLE PRECISION DUMMY1,EVEC(NIMAGE),ETOTAL,FIN(3*NATOMS), &
     &                 XX(3*NATOMS),GG(3*NATOMS),RMS,G(*),X(*), &
     &                 QSTART(3*NATOMS),RMSVEC(NIMAGE),SUM,EINITIAL,EFINAL,SEPARATION, &
     &                 LARG(0:NIMAGE+1,0:NIMAGE+1), LDIST, DIFF, SDIFF, NORM1, NORM2
COMMON /NEBRMS/ RMS,EINITIAL,EFINAL,SEPARATION
INTEGER IADD, LNOPT, LNATOMS, ISTAT
COMMON /INTS/ IADD, LNOPT, LNATOMS
DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793D0

NORM1=1.0D0/NIMAGE
NORM2=2.0D0/((NIMAGE+1)*(NIMAGE+2))
LDIST=X(LNOPT*NIMAGE+1)
ETOTAL=0.0D0
DO J1=1,NIMAGE
   DO J2=1,LNOPT
      XX(J2)=X(J2+LNOPT*(J1-1))
   ENDDO
!CALL VAR_TO_GEOM(LNOPT,XX(1:LNOPT))
!CALL CHAINBUILD
   CALL POTENTIAL(XX,EVEC(J1),GG,.TRUE.,.FALSE.,RMSVEC(J1),.FALSE.,.FALSE.)
   DO J2=1,LNOPT
      G(J2+LNOPT*(J1-1))=NORM1*GG(J2)
   ENDDO
   CALL FLUSH(6,ISTAT)
   ETOTAL=ETOTAL+EVEC(J1)
!  PRINT*,'IMAGE,EVEC,ETOTAL=',J1,EVEC(J1),ETOTAL
ENDDO

SUM=ETOTAL*NORM1

DO J2=1,NIMAGE
   DUMMY1=0.0D0
   DO J1=1,LNOPT
      DIFF=X(J1+(J2-1)*LNOPT)-QSTART(J1)
      SDIFF=DIFF-INT(DIFF/PI)*2*PI
      DUMMY1=DUMMY1+SIN(SDIFF/(2.0D0*J2))**2 
   ENDDO
   LARG(J2,0)=DUMMY1-LDIST**2 
   SUM=SUM+NORM2*MECLAMBDA*LARG(J2,0)**2  
ENDDO

DO J2=1,NIMAGE
   DO J3=J2+1,NIMAGE
      DUMMY1=0.0D0
      DO J1=1,LNOPT
         DIFF=X(J1+LNOPT*(J2-1))-X(J1+LNOPT*(J3-1))
         SDIFF=DIFF-INT(DIFF/PI)*2*PI
         DUMMY1=DUMMY1+SIN(SDIFF/(2.0D0*(J2-J3)))**2 
      ENDDO
      LARG(J3,J2)=DUMMY1-LDIST**2
      SUM=SUM+NORM2*MECLAMBDA*LARG(J3,J2)**2
   ENDDO
ENDDO

DO J2=1,NIMAGE
   DUMMY1=0.0D0
   DO J1=1,LNOPT
      DIFF=X(J1+LNOPT*(J2-1))-FIN(J1)
      SDIFF=DIFF-INT(DIFF/PI)*2*PI
      DUMMY1=DUMMY1+SIN(SDIFF/(2.0D0*(J2-NIMAGE-1)))**2 
   ENDDO
   LARG(NIMAGE+1,J2)=DUMMY1-LDIST**2
   SUM=SUM+NORM2*MECLAMBDA*LARG(NIMAGE+1,J2)**2
ENDDO

RMS=0.0D0
DO J1=1,LNOPT ! K
   DO J2=1,NIMAGE ! GAMMA
      DIFF=X(J1+LNOPT*(J2-1))-FIN(J1)
      SDIFF=DIFF-INT(DIFF/PI)*2*PI
      DUMMY1=LARG(NIMAGE+1,J2)*SIN(SDIFF/(J2-NIMAGE-1))/(J2-NIMAGE-1)
      DIFF=X(J1+LNOPT*(J2-1))-QSTART(J1)
      SDIFF=DIFF-INT(DIFF/PI)*2*PI
      DUMMY1=DUMMY1+LARG(J2,0)*SIN(SDIFF/J2)/J2  
      DO J3=1,J2-1
         DIFF=X(J1+LNOPT*(J2-1))-X(J1+LNOPT*(J3-1))
         SDIFF=DIFF-INT(DIFF/PI)*2*PI
         DUMMY1=DUMMY1+LARG(J2,J3)*SIN(SDIFF/(J2-J3))/(J2-J3)
      ENDDO
      DO J3=J2+1,NIMAGE
         DIFF=X(J1+LNOPT*(J2-1))-X(J1+LNOPT*(J3-1))
         SDIFF=DIFF-INT(DIFF/PI)*2*PI
         DUMMY1=DUMMY1+LARG(J3,J2)*SIN(SDIFF/(J2-J3))/(J2-J3)
      ENDDO
      G(J1+LNOPT*(J2-1))=G(J1+LNOPT*(J2-1))+NORM2*MECLAMBDA*DUMMY1
      RMS=RMS+G(J1)**2
   ENDDO
ENDDO

G(LNOPT*NIMAGE+1)=0.0D0
LARG(NIMAGE+1,0)=0.0D0
DO J1=0,NIMAGE
   DO J2=J1+1,NIMAGE+1
      G(LNOPT*NIMAGE+1)=G(LNOPT*NIMAGE+1)+LARG(J2,J1)
   ENDDO
ENDDO
G(LNOPT*NIMAGE+1)=-4*NORM2*MECLAMBDA*LDIST*G(LNOPT*NIMAGE+1)
RMS=RMS+G(LNOPT*NIMAGE+1)**2
RMS=SQRT(RMS/(LNOPT*NIMAGE+1))

IF (.FALSE.) THEN ! POTENTIALLY USEFUL DEBUG PRINTING
   DUMMY1=0.0D0
   DO J1=1,LNOPT
      DUMMY1=DUMMY1+SIN((X(J1)-QSTART(J1))/2.0D0)**2
   ENDDO
   PRINT '(A,I5,4G20.10)', 'IMAGE,E,RMS,DIST^2,D^2=',0,0.0D0,0.0D0,DUMMY1,LDIST**2
   DO J1=1,NIMAGE-1
      DUMMY1=0.0D0
      DO J2=1,LNOPT
         DUMMY1=DUMMY1+SIN((X(J2+LNOPT*(J1-1))-X(J2+LNOPT*J1))/2.0D0)**2
      ENDDO
      PRINT '(A,I5,4G20.10)', 'IMAGE,E,RMS,DIST^2,D^2=',J1,EVEC(J1),RMSVEC(J1),DUMMY1,LDIST**2
   ENDDO
   DUMMY1=0.0D0
   DO J1=1,LNOPT
      DUMMY1=DUMMY1+SIN((X(J1+LNOPT*(NIMAGE-1))-FIN(J1))/2.0D0)**2
   ENDDO
   PRINT '(A,I5,4G20.10)', 'IMAGE,E,RMS,DIST^2,D^2=',NIMAGE,EVEC(NIMAGE),RMSVEC(NIMAGE),DUMMY1,LDIST**2
ENDIF

RETURN
END SUBROUTINE UNMAKEGRADMEC
