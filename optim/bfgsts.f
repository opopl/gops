C
! GPL License Info {{{
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
C
C   }}}
C
!BFGSTS: hybrid eigenvector-following/gradient minimization {{{
C***********************************************************************
C
C  This subroutine is designed to perform a hybrid eigenvector-following/gradient
C  minimization optimization for large systems where we do not want to 
C  diagonalize the whole Hessian.
C
C***********************************************************************
C
      SUBROUTINE BFGSTS(ITMAX,COORDS,ENERGY,GRAD,MFLAG,RMS,EVALMIN,EVALMAX,VECS,ITER,POTCALL,PTEST)
! Doxygen comments {{{
!> \name BFGSTS
!> \brief perform hybrid eigenvector-following/gradient minimization for large systems 
!> \param ITMAX
!> \param COORDS
!> \param ENERGY
!> \param GRAD
!> \param MFLAG
!> \param RMS 
!> \param EVALMIN
!> \param EVALMAX
!> \param VECS
!> \param ITER
!> \param POTCALL
!> \param PTEST
! }}}
! Modules {{{
      USE COMMONS
      USE KEY
      USE VECCK
      USE ZWK
      USE MODCHARMM
      USE MODHESS
      USE PORFUNCS
! }}}
      IMPLICIT NONE
! subroutine parameters {{{

      INTEGER ITMAX, ITER 
      LOGICAL POTCALL, PTEST, MFLAG
      DOUBLE PRECISION, DIMENSION(3*NATOMS) :: COORDS, GRAD, VECS
      DOUBLE PRECISION  ENERGY,RMS, EVALMAX, EVALMIN 
!}}}
! local parameters  {{{

      INTEGER J1, J2, INEG, J, NS, I, K1, NBFGS, FRAME, ITDONE, J3, NSTEPMINSAVE
      DOUBLE PRECISION SCRATCH(6*NATOMS),SUM,AVG,FOBNEW,DUMMY,
     1                 FOB,PSTEP,DPRAND,FIXDIR(3*NATOMS),TEMPA(9*NATOMS),
     2                 STPMAG,COORDSN(3*NATOMS),VECSP(3*NATOMS),
     3                 XP1,XP2,E1,E2,DELE,SCALE,EPER,EOLD,RMS2,SOVER,EREAL,XRAT(3*NATOMS),
     4                 XFOB(3*NATOMS), STEP(3*NATOMS), XPSTEP(3*NATOMS), CSTEP(3*NATOMS), PFOB(3*NATOMS),
     5                 AV(6),SSTPMAG,RAT1,RAT2,TEMP,LP1,LP2,LP,VECL(3*NATOMS),DOTOPT,PROJ1,PROJ2
C     EXTERNAL OP, IOVECT
      CHARACTER(LEN=80) :: FNAME
      LOGICAL AWAY, LINE, PVFLAG, RESET, STEST,  FIXDSAVE
      COMMON /PVF/ PVFLAG
      LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST, CONVERGED
      INTEGER NCONNECT
      DOUBLE PRECISION TEMPERATURE, HRED
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
C
C  Assign enough memory to WORK for a blocksize of 32 to be possible.
C  This is for DSYEVR.
C
      INTEGER IWORK(33*3*NATOMS)
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NATOMS)
      INTEGER INFO, ISTAT
      DOUBLE PRECISION WORK(33*3*NATOMS), ABSTOL, DIAG(3*NATOMS), DLAMCH
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
! }}}
      SAVE
! Subroutine body  {{{

      LWORK=33*3*NATOMS
      ILWORK=33*3*NATOMS

      IF ((ZSYM(1)(1:1).EQ.'W').AND.(.NOT.BFGSSTEP)) THEN
         PRINT*,'BFGSTS procedures have not been programmed for TIP potentials'
         CALL FLUSH(6,ISTAT)
         STOP
      ENDIF
      FIXDSAVE=FIXD
      IF ((HINDEX.GT.1).AND.(.NOT.NOIT)) THEN
         WRITE(*,'(A)') 'For HINDEX > 1 you must use NOIT with BFGSTS'
         CALL FLUSH(6,ISTAT)
         STOP
      ENDIF
!
! Check for consistent convergence criteria on RMS gradient in MYLBFGS and EF part.
!     
      IF (GMAX.NE.CONVR) THEN
         IF (DEBUG) PRINT '(2(A,G20.10),A)','bfgsts> WARNING - GMAX ',GMAX,' is different from CONVR ',CONVR,' - resetting'
         GMAX=MIN(CONVR,GMAX)
         CONVR=MIN(CONVR,GMAX)
      ENDIF
      IF (MINMAX.GE.CONVU) THEN
         IF (DEBUG) PRINT '(2(A,G20.10),A)','bfgsts> WARNING - CONVU <= MINMAX - resetting CONVU'
         CONVU=2*MINMAX
      ENDIF
!
!  Reset maximum step sizes in case this isn't the first call to EFOL.
!
      IF (DEBUG) PRINT '(A,G20.10)',' bfgsts> resetting maximum step sizes to ',MXSTP
      DO J1=1,NOPT
         STPMAX(J1)=MXSTP
      ENDDO

      ITER=1
      frame=1
      CONVERGED=.TRUE.
C
C  If DUMPV is .TRUE. then vector.dump is already open and attached to unit 44.
C  SGI compiler won;t allow us to attach it to another unit.
C  vector.dump could contain multiple dumps for more than the last step, so we 
C  have to make sure we get the results for the last step.
C
      IF (READV.AND.(ITER.EQ.1)) THEN   
         IF (DUMPV) THEN
            REWIND(44)
211         READ(44,*,END=111) EVALMIN
            READ(44,*) (VECS(J1),J1=1,NOPT)
            GOTO 211
111         CONTINUE
         ELSE
            IF (FILTH.EQ.0) THEN
               FNAME='vector.dump'
            ELSE
               WRITE(FNAME,'(A)') 'vector.dump.'//TRIM(ADJUSTL(FILTHSTR))
            ENDIF
   
            OPEN(UNIT=45,FILE=FNAME,STATUS='OLD')
20          READ(45,*,END=10) EVALMIN
            READ(45,*) (VECS(J1),J1=1,NOPT)
            GOTO 20
10          CLOSE(45)
         ENDIF
         WRITE(*,'(A,F20.10)') ' Reaction vector read successfully. Eigenvalue=   ',EVALMIN
         NS=100
      ENDIF
      IF (FREEZE) THEN
         DO J1=1,NATOMS
            IF (FROZEN(J1)) THEN
               VECS(3*(J1-1)+1)=0.0D0
               VECS(3*(J1-1)+2)=0.0D0
               VECS(3*(J1-1)+3)=0.0D0
            ENDIF
         ENDDO
      ENDIF

      STEST=.TRUE.
      IF (NOHESS) STEST=.FALSE.
      IF (READV.AND.(ITER.EQ.1)) STEST=.FALSE.
90    IF (PTEST) PRINT*
      NUP=HINDEX
      IF ((.NOT.BFGSSTEP).AND.PTEST) WRITE(*,11) ITER
11          FORMAT (' BFGSTS> Beginning of optimization cycle ', I4,'.',/
     1              ' ---------------------------------------')
      FIXIMAGE=.FALSE.
      IF ((FIXAFTER.GT.0).AND.(ITER.GT.FIXAFTER)) FIXIMAGE=.TRUE.
      IF (CHRMMT) NCHENCALLS = 999 ! make sure non-bond list is updated at the start of each TS search.
      IF (POTCALL) THEN
         IF (PV.AND.(.NOT.BFGSSTEP)) THEN
            IF (.NOT.KNOWE) CALL POTENTIAL(COORDS,ENERGY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PVFLAG=.FALSE.
            CALL PVOPT(COORDS,ENERGY,GRAD)
         ENDIF
         IF ((.NOT.KNOWG).OR.(STEST.AND.(.NOT.KNOWH))) CALL POTENTIAL(COORDS,ENERGY,GRAD,.TRUE.,STEST,RMS,PTEST,.FALSE.)
      ENDIF
      CALL DUMPP(COORDS,ENERGY)
      IF ((.NOT.VARIABLES).AND.NOIT.AND.STEST) THEN
         IF (RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')) THEN
         ELSE IF (ZSYM(NATOMS).EQ.'SY') THEN
            CALL SHIFTSTOCK(COORDS,NATOMS)
         ELSEIF (RBAAT) THEN
            CALL SHIFTRIGID(COORDS,NATOMS)
         ELSEIF (ZSYM(NATOMS).EQ.'TH') THEN
            CALL SHIFTHTH(COORDS,NATOMS)
         ELSE
            CALL SHIFTH(COORDS,.TRUE.,NOPT,NATOMS,ATMASS)
            IF (TWOD) THEN
               DO J1=1,NATOMS
                  J2=3*J1
                  DO J3=1,NATOMS
                     HESS(J2,3*(J3-1)+1)=0.0D0
                     HESS(J2,3*(J3-1)+2)=0.0D0
                     HESS(J2,3*(J3-1)+3)=0.0D0
                     HESS(3*(J3-1)+1,J2)=0.0D0
                     HESS(3*(J3-1)+2,J2)=0.0D0
                     HESS(3*(J3-1)+3,J2)=0.0D0
                  ENDDO
                  HESS(J2,J2)=SHIFTV
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      IF (.NOT.((READV.AND.BFGSSTEP).OR.(.NOT.POTCALL))) THEN
C        IF (DEBUG) PRINT*,'This line is needed to prevent miscompilation by Sun Forte Developer 7' ! SIGH !!!
         IF (FIXD.AND.(ITER.EQ.1)) THEN
            IF (CONNECTT) THEN
               DO J1=1,NOPT
                  FIXDIR(J1)=VECS(J1)
               ENDDO
            ELSE
               DO J1=1,NOPT
                  DUMMY=2*(DPRAND()-0.5D0)
                  FIXDIR(J1)=DUMMY
               ENDDO
            ENDIF
            CALL ORTHOGOPT(FIXDIR,COORDS,.TRUE.)
            IF (.TRUE.) CALL HSMOVE(COORDS,COORDSN,FIXDIR,PSTEP,PTEST)
            EVALMIN=1.0D0
         ELSE IF (FIXD) THEN
            CALL ORTHOGOPT(FIXDIR,COORDS,.TRUE.)
            IF (.TRUE.) CALL HSMOVE(COORDS,COORDSN,FIXDIR,PSTEP,PTEST)
            EVALMIN=1.0D0
         ENDIF
         IF (.NOT.NOHESS) THEN 
            IF (.NOT.NOIT) THEN
               CALL ITEIG(ITER,COORDS,VECS,EVALMIN,EVALMAX,NS,SOVER,PTEST,VECL,CONVERGED)
            ELSE
               ABSTOL=DLAMCH('Safe  minimum')
               IF (ITER.GT.1) THEN
                  DO J1=1,NOPT
                     VECSP(J1)=VECS(J1)
                  ENDDO
               ENDIF
               CALL DSYEVR('V','I','U',NOPT,HESS,SIZE(HESS,1),0.0D0,1.0D0,1,HINDEX,ABSTOL,NFOUND,DIAG,
     1                        ZWORK,3*NATOMS,ISUPPZ,WORK,
     2                        LWORK, IWORK, ILWORK, INFO )
               IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEVR'
C              PRINT*,'Optimal and actual values of LWORK=',WORK(1),LWORK
C              PRINT*,'Optimal and actual values of ILWORK=',IWORK(1),ILWORK
               EVALMIN=DIAG(1)
               SOVER=0.0D0
               DO J1=1,NOPT
                  VECS(J1)=ZWORK(J1,1)
                  IF (ITER.GT.1) SOVER=SOVER+VECS(J1)*VECSP(J1)
               ENDDO
               IF (PTEST) THEN
                  IF (ITER.GT.1) THEN
                     WRITE(*,'(A,F15.7,A,F15.7)') ' Smallest eigenvalue=',EVALMIN,' overlap with previous vector=',SOVER
                  ELSE
                     WRITE(*,'(A,F15.7,A,F15.7)') ' Smallest eigenvalue=',EVALMIN
                  ENDIF
               ENDIF
               DO I=1,NOPT
                  SCRATCH(I) = COORDS(I)
               ENDDO
               DO I=1,HINDEX
                  XFOB(I)=0.0D0
                  DO J = 1,NOPT
                     XFOB(I)=XFOB(I)+GRAD(J)*ZWORK(J,I)
                  ENDDO
               ENDDO
C
C  Find the vector of STPMAX values by comparing predicted and
C  actual second derivatives for each eigenvector.
C
               IF ((ITER.GT.1).AND.(ISTCRT.EQ.10)) THEN
                  DO J1=1,NOPT
                     TEMPA(J1)=STPMAX(J1)
                  ENDDO
                  K1=0
                  DO J1=1,HINDEX
                     XRAT(J1)=0.0D0
                     K1=K1+1
                     IF (DABS(XPSTEP(K1)).GT.1.0D-40) THEN
C
C  Allow for possible phase change in the eigenvector. Just take the smaller value.
C
                        RAT1=DABS(( XFOB(J1)-PFOB(K1))/(DIAG(J1)*XPSTEP(K1))-1.0D0)
                        RAT2=DABS((-XFOB(J1)-PFOB(K1))/(DIAG(J1)*XPSTEP(K1))-1.0D0)
                        XRAT(J1)=MIN(RAT1,RAT2)
C                       WRITE(*,'(A,2I4,5E15.7)') 'J1,K1,FOB,PFOB,PSTEP,RAT1,DIAG=', J1,K1,XFOB(J1),PFOB(K1),PSTEP(K1),RAT1,DIAG(J1)
                        IF (XRAT(J1).GT.TRAD) THEN
                           STPMAX(J1)=MAX(TEMPA(K1)/1.11D0,MINMAX)
                        ELSE
                           IF (MASST) THEN
                              SUM=0
                              DO J2=1,NATOMS
                                 SUM=SUM+ATMASS(J2)
                              ENDDO
                              AVG=SQRT(SUM/NATOMS)
C                             PRINT *,'the average is',AVG
                              STPMAX(J1)=MIN(MAX(TEMPA(K1)*1.09D0,MINMAX),AVG*MAXMAX)
                           ELSE
                              STPMAX(J1)=MIN(MAX(TEMPA(K1)*1.09D0,MINMAX),MAXMAX)
                           ENDIF
                        ENDIF
                     ELSE
                        STPMAX(J1)=MAX(TEMPA(K1),MINMAX)
                     ENDIF
                  ENDDO
               ENDIF
               TEMP=-1.0D0
            ENDIF
C
C  Do not use.
C
C           CALL ITEIGN(ITER,COORDS,VECS,EVALMIN,EVALMAX,PTEST)
C
C  Lanczos routine - does work but seems to be slower? For NFIG=3 it is only
C  a bit slower than ITEIG with CEIG=0.01.
C
C           NVAL=-10
C           ANV=ABS(NVAL)
C           NFIG=3
C           NPERM=0
C           IF (ITER.GT.1) NPERM=1
C           MAXOP=NEVS
C           NBLOCK=1
C           IF ((ANV.GT.MANV).OR.(NBLOCK.GT.MAXBLOCK)) THEN
C              WRITE(*,'(A)') ' Too many eigenvectors or too large a Lanczos blocksize requested'
c              STOP
C           ENDIF
C
C           CALL DNLASO(HESS, Q, NATOMS, OP, IOVECT, 3*NATOMS, NVAL, NFIG, NPERM,
C    *                  NMVAL, VAL, NMVEC, LVEC, NBLOCK, MAXOP, MAXJ, LWORK,
C    *                  IND, IERR)
C           WRITE(*,'(I3,A,100F20.10)') NPERM,' eigenvectors determined: ',(VAL(J1,1),J1=1,NPERM)
C           IF (IERR.NE.0) THEN
C              WRITE(*,'(A,I4)') 'Lanczos call completed with error code ',IERR
C           ELSE
C              WRITE(*,'(A,F20.10)') ' Smallest eigenvalue=',VAL(1,1)
C           ENDIF

C           EVALMIN=VAL(1,1)
C           DO J1=1,NOPT
C              VECS(J1)=LVEC(J1,1)
C           ENDDO
         ELSE
            CALL BEIG(ITER,COORDS,ENERGY,VECS,EVALMIN,NS,SOVER,PTEST,CONVERGED)

C
C  The following two routines are also obsolete.
C
C           CALL POWEIG(ITER,COORDS,ENERGY,VECS,EVALMIN)
C           CALL MCEIG(ITER,COORDS,ENERGY,VECS,EVALMIN)
         ENDIF
      ENDIF

      IF ((EVALMIN.LT.0.0D0).AND.FIXD) THEN
         IF (PTEST) WRITE(*,'(A)') ' Negative eigenvalue, changing to hybrid EF'
         FIXD=.FALSE.
         CALL VECNORM(VECS,NOPT)
         DO J1=1,NOPT
            ZWORK(J1,1)=VECS(J1)
         ENDDO
      ELSE IF (FIXD) THEN
C     IF (FIXD) THEN
         DO J1=1,NOPT
            ZWORK(J1,1)=FIXDIR(J1)
            VECS(J1)=FIXDIR(J1)
         ENDDO
      ELSE
         CALL VECNORM(VECS,NOPT)
         DUMMY=0.0D0
         DO J1=1,NOPT
            ZWORK(J1,1)=VECS(J1)
            DUMMY=DUMMY+VECS(J1)**2
         ENDDO
      ENDIF
C
C  Dump the smallest non-zero eigenvalue and eigenvector 
C  in file vector.dump, if required.
C
      IF (DUMPV) THEN
         IF (.NOT.ALLSTEPS) REWIND(44)
         WRITE(44,'(E20.10)') EVALMIN
         WRITE(44,'(3F20.10)') (ZWORK(J1,1),J1=1,NOPT)
         CALL FLUSH(44,ISTAT)
      ENDIF
C
C  Take an eigenvector-following step uphill along the direction corresponding
C  to the smallest non-zero eigenvalue. Then do a line minimization along the
C  gradient vector with component along the uphill direction projected out. 
C  Do we need a new gradient vector or will the old one from the point before
C  stepping uphill do?
C
      FOB=0.0D0
      DO J1=1,NOPT
         FOB=FOB+ZWORK(J1,1)*GRAD(J1)
      ENDDO
      XFOB(1)=FOB
      REDOPATH1=.FALSE.; REDOPATH2=.FALSE.
      IF (REDOPATH.AND.BFGSSTEP) THEN
         PROJ1=0.0D0
         PROJ2=0.0D0
         DO J1=1,NOPT
            PROJ1=PROJ1+ZWORK(J1,1)*(MIN1REDO(J1)-COORDS(J1))/D1INIT
            PROJ2=PROJ2+ZWORK(J1,1)*(MIN2REDO(J1)-COORDS(J1))/D2INIT
         ENDDO
         IF (IVEC.LT.0) THEN
            PROJ1=-PROJ1
            PROJ2=-PROJ2
         ENDIF
         IF (PROJ1.GT.PROJ2) REDOPATH1=.TRUE.
         IF (PROJ2.GT.PROJ1) REDOPATH2=.TRUE.
         PRINT '(A,2G20.10)',' bfgsts> projections of step onto directions of minima 1 and 2 are: ',PROJ1,PROJ2 
      ENDIF
      IF (FIXD) GOTO 666
      IF (HINDEX.LE.1) THEN
         IF (EVALMIN.LT.0.0D0) THEN
C           PRINT*,'There is at least one negative eigenvalue'
            INEG=1
         ELSE
C           PRINT*,'There are no negative eigenvalues'
            INEG=0
         ENDIF
      ELSE
         INEG=0
         DO J1=1,HINDEX
            IF (DIAG(J1).LT.0.0D0) INEG=INEG+1
         ENDDO
         WRITE(*,'(A,I5,A)') ' There are at least ',INEG,' negative Hessian eigenvalues'
      ENDIF
C
C  Take a step away from a stationary point along the appropriate
C  Hessian eigenvector. This enables us to start from converged minima.
C  Distinguish the case where we want to take a very small step away
C  from a transition state from others where we want a big displacement
C  to get unstuck. 
C
      AWAY=.FALSE.
      IF ((RMS.LT.PUSHCUT).AND.(.NOT.FIXD)) THEN
         IF ((INEG.NE.HINDEX).AND.CONVERGED) THEN ! don;t try pushoff if CONVERGED is .FALSE.
            IF (HINDEX.LE.1) THEN
               IF ((ITER.EQ.1).AND.(INEG.EQ.0)) THEN
                  IF (PTEST) THEN
                     IF (IVEC.GE.0) PRINT*,'Stepping away from minimum along softest mode + direction'
                     IF (IVEC.LT.0) PRINT*,'Stepping away from minimum along softest mode - direction'
                  ENDIF
                  AWAY=.TRUE.
               ELSE IF (MOD(ITER-1,4).EQ.0) THEN
                  PRINT*,'Stepping away from solution of wrong index'
                  IF (PUSHOFF.EQ.0.0D0) THEN
                     FOB=STPMAX(1)
                  ELSE 
                     FOB=PUSHOFF
                  ENDIF
               ENDIF
            ELSE
               DO J1=1,HINDEX
                  IF (DIAG(J1).GT.0.0D0) THEN
                     WRITE(*,'(A,I6)') 'Stepping away from solution of wrong index along eigenvector ',J1
                     IF (PUSHOFF.EQ.0.0D0) THEN
                        XFOB(J1)=STPMAX(J1)
                     ELSE
                        XFOB(J1)=PUSHOFF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      IF (BFGSSTEP) THEN
         IF (INEG.EQ.0) THEN
            PRINT*,'****WARNING - BFGSSTEP set for a minimum'
         ELSE
            IF (PTEST) THEN
               IF (IVEC.GE.0) PRINT*,'Stepping away from saddle along softest mode + direction'
               IF (IVEC.LT.0) PRINT*,'Stepping away from saddle along softest mode - direction'
            ENDIF
            AWAY=.TRUE.
         ENDIF
      ENDIF
C
C  EF determination of steps
C
      IF (HINDEX.LE.1) THEN
         XP1=DABS(EVALMIN)/2.0D0
         IF (MASST) THEN
            SUM=0
            DO J2=1,NATOMS
               SUM=SUM+ATMASS(J2)
            ENDDO
            AVG=SUM/NATOMS
            XP2=1.0D0 + 4.0D0*(FOB/EVALMIN)**2/AVG
         ELSE
            XP2=1.0D0 + 4.0D0*(FOB/EVALMIN)**2
         ENDIF
         XP1=-XP1*(1.0D0+DSQRT(XP2))
         PSTEP=-FOB/XP1
   
         IF (AWAY) THEN
            IF (PUSHOFF.NE.0.0D0) THEN
               PSTEP=PUSHOFF
            ELSE
               PSTEP=STPMAX(1)/1.0D1
            ENDIF
            IF (IVEC.LT.0) PSTEP=-PSTEP
         ENDIF
C        WRITE(*,'(A,F20.10)') ' Unscaled step=',PSTEP
C
C  Scale according to step size in ev basis:
C
         STPMAG=ABS(PSTEP)
         IF (.NOT.AWAY) THEN
            SCALE=MIN(STPMAX(1)/MAX(DABS(PSTEP),1D-10),1.0D0)
            PSTEP=SCALE*PSTEP
         ENDIF
         STEP(1)=PSTEP ! because we save this value in XPSTEP later from STEP array
      ELSE
         DO I=HINDEX,1,-1
            STEP(I)=0.0D0
            NZERO=0
            IF (DABS(DIAG(NZERO+1)/MAX(DABS(DIAG(I)),1.0D-10)).GT.TEMP) TEMP=DABS(DIAG(NZERO+1)/MAX(DABS(DIAG(I)),1.0D-10))
            LP1=DABS(DIAG(I))/2.0D0
            IF (MASST) THEN
               SUM=0
               DO J2=1,NATOMS
                  SUM=SUM+ATMASS(J2)
               ENDDO
               AVG=SUM/NATOMS
               LP2=1.0D0 + 4.0D0*((XFOB(I)/DIAG(I))**2)/AVG
            ELSE
               LP2=1.0D0 + 4.0D0*(XFOB(I)/DIAG(I))**2
            ENDIF
            LP=LP1*(1.0D0+DSQRT(LP2))
            WRITE(*,'(A,I4,A,4X,F19.10)') ' Mode ',I,' will be searched uphill. Eigenvalue=',DIAG(I)
            LP=-LP
            STEP(I)=-XFOB(I)/LP
         ENDDO
         DO J1=1,NOPT
            SCRATCH(NOPT+J1)=STEP(J1)
         ENDDO
C
C  Convert the steps to the Cartesian rather than the EV basis
C
         DO J=1,NOPT
            SCRATCH(J+NOPT)=0.0D0
            DO I=1,HINDEX
               SCRATCH(J+NOPT)=SCRATCH(J+NOPT)+STEP(I)*ZWORK(J,I)
            ENDDO
         ENDDO
C
C  Scale according to step size in ev basis:
C
         CALL VSTAT(STEP(1),AV,NOPT,3*NATOMS)
         STPMAG=MAX(AV(1),1D-10)
         DO J1=1,NOPT
            SCALE=MIN(STPMAX(J1)/MAX(DABS(STEP(J1)),1D-10),1.0D0)
            STEP(J1)=SCALE*STEP(J1)
         ENDDO
         CALL VSTAT(STEP(1),AV,NOPT,3*NATOMS)
         SSTPMAG=MAX(AV(1),1D-10)
         SCALE=1.0D0
         E1=0.0D0
         DO J1=1,NOPT
            E1=E1+STEP(J1)*XFOB(J1)
         ENDDO
         E2=0.0D0
         DO J1=1,NOPT
            E2=E2+DIAG(J1)*STEP(J1)**2
         ENDDO
         E2=E2/2.0D0
         DELE=E1*SCALE+E2*SCALE**2
         CALL VADD(SCRATCH(1),SCRATCH(1),SCRATCH(NOPT+1),NOPT,1)

         IF (EFSTEPST.AND.(MOD(ITER-1,EFSTEPS).EQ.0)) THEN
            DO I=NZERO+1,HINDEX
                WRITE(*,360) I, STEP(I)
360             FORMAT(' Unscaled step for mode ',I3,'=',F20.10)
            ENDDO
         ENDIF
      ENDIF
C
C  Use MAXMAX until we have a negative eigenvalue. Be sure to
C  step uphill!
C
      LINE=.FALSE.
      IF (CONVERGED) THEN
         IF (((EVALMIN.GT.0.0D0).OR.FIXD).AND.(HINDEX.LE.1)) THEN
            PSTEP=FOB*MAXMAX/ABS(FOB)
            STEP(1)=PSTEP
         ELSE IF (HINDEX.GT.1) THEN
            DO J1=1,HINDEX
               IF (DIAG(J1).GT.0.0D0) STEP(J1)=XFOB(J1)*MAXMAX/ABS(XFOB(J1))
            ENDDO
         ENDIF
      ENDIF

      IF (HINDEX.LE.1) THEN
         E1=PSTEP*FOB
         E2=EVALMIN*PSTEP**2/2.0D0
         DELE=E1+E2
      ENDIF
C
C  Regenerate full Q vector. 
C
666   IF (.NOT.FIXD) THEN
         IF (HINDEX.LE.1) THEN
            DO J=1,NOPT
               COORDS(J)=COORDS(J)+PSTEP*ZWORK(J,1)
            ENDDO
         ELSE
C
C  Unpack SCRATCH(1:NOPT) to regenerate full Q vector.
C  CSTEP contains the step in the Cartesian basis.
C
            DO J=1,NOPT
               COORDS(J)=SCRATCH(J)
               CSTEP(J)=SCRATCH(NOPT+J)
            ENDDO
         ENDIF
      ELSE
         DO J=1,NOPT
            COORDS(J)=COORDSN(J)
         ENDDO
      ENDIF

!     WRITE(667,'(I8)') NATOMS
!     WRITE(667,'(F20.10)') EVALMIN
!     DO J1=1,NATOMS
!        WRITE(667,'(A2,1X,6F20.10)') ZSYM(J1),COORDS(3*(J1-1)+1),COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3),
!    &                            ZWORK(3*(J1-1)+1,1),ZWORK(3*(J1-1)+2,1),ZWORK(3*(J1-1)+3,1)
!     ENDDO
!     CALL FLUSH(667,ISTAT)
      KNOWE=.FALSE.
      KNOWG=.FALSE.
      KNOWH=.FALSE.

      IF (.NOT.LINE) THEN
         EOLD=ENERGY
         IF (.NOT.BFGSSTEP) THEN
C
C  Optimising the box lengths here would change the critical eigenvalue and
C  eigenvector.
C
C           IF (PV) THEN
C              CALL POTENTIAL(COORDS,ENERGY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
C              PVFLAG=.FALSE.
C              CALL PVOPT(COORDS,ENERGY,GRAD)
C           ENDIF
C
            CALL POTENTIAL(COORDS,ENERGY,GRAD,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            KNOWE=.TRUE.
            KNOWG=.TRUE.
            CALL DUMPP(COORDS,ENERGY)
         ENDIF
      ENDIF
C
C  If REOPT is true then reoptimise the eigenvector corresponding to the uphill direction
C  before doing the tangent space minimisation. Does not make sense with FIXD.
C
      IF (REOPT) THEN
         WRITE(*,'(A)') 'Reoptimising uphill direction after EF step'
         IF (.NOT.NOHESS) THEN
            IF (NOIT) THEN
               PRINT*,'reoptimisation with NOIT not implemented'
               CALL FLUSH(6,ISTAT)
               STOP
            ENDIF
            CALL ITEIG(2,COORDS,VECS,EVALMIN,EVALMAX,NS,SOVER,PTEST,VECL,CONVERGED)
         ELSE
            CALL BEIG(2,COORDS,ENERGY,VECS,EVALMIN,NS,SOVER,PTEST,CONVERGED)
         ENDIF
         DO J1=1,NOPT
            ZWORK(J1,1)=VECS(J1)
         ENDDO
      ENDIF
      CALL FLUSH(6,ISTAT)

      FOBNEW=0.0D0
      DO J1=1,NOPT
         FOBNEW=FOBNEW+GRAD(J1)*ZWORK(J1,1)
      ENDDO
C
C  Only scale if we have a -ve eigenvalue
C
      IF ((EOLD.NE.0.0D0).AND.(EVALMIN.LT.0.0D0).AND.(.NOT.FIXD)) THEN
         EPER=MIN(DABS(1.0D0-(FOBNEW-FOB)/(PSTEP*EVALMIN)),DABS(1.0D0-(-FOBNEW-FOB)/(PSTEP*EVALMIN)))
C        WRITE(*,'(A,4F20.10)') 'FOB,FOBNEW,PSTEP,EVALMIN=',FOB,FOBNEW,PSTEP,EVALMIN
C        WRITE(*,'(A,3F20.10)') 'EPER,EPER1,EPER2=',
C    1          EPER,DABS(1.0D0-(FOBNEW-FOB)/(PSTEP*EVALMIN)),DABS(1.0D0-(-FOBNEW-FOB)/(PSTEP*EVALMIN))
         IF (EPER.GT.TRAD) THEN
            STPMAX(1)=MAX(STPMAX(1)/1.1D0,MINMAX)
C           WRITE(*,'(A,E12.4,A,E12.4)') ' Decreasing allowed EF step; trust radius=',TRAD,' calculated ratio=',EPER
         ELSE
            STPMAX(1)=MIN(STPMAX(1)*1.1D0,MAXMAX)
C           WRITE(*,'(A,E12.4,A,E12.4)') ' Increasing allowed EF step; trust radius=',TRAD,' calculated ratio=',EPER
         ENDIF
      ENDIF
C
C Summarize
C
C     IF (SUMMARYT.AND.(MOD(ITER-1,NSUMMARY).EQ.0)) THEN
         IF (PTEST) THEN
            IF (HINDEX.LE.1) THEN
               WRITE(*,30)
30             FORMAT(1X,79('-'))
               WRITE(*,40)
40             FORMAT(' Vector      Gradient        Secder       Step          Max step    Trust ratio')
               WRITE(*,30)
               WRITE(*,50) 1,FOB,EVALMIN,PSTEP,STPMAX(1),EPER
50             FORMAT(1X,I4,1X,E15.6,1X,E15.6,1X,E13.6,1X,E12.6,1X,E15.6)
               WRITE(*,30)
            ELSE
               WRITE(*,30)
               WRITE(*,40)
               WRITE(*,30)
               DO I=HINDEX,1,-1
                  WRITE(*,50) I,XFOB(I),DIAG(I),STEP(I),STPMAX(I),XRAT(I)
               ENDDO
               WRITE(*,30)
            ENDIF
         ELSE IF (.NOT.BFGSSTEP) THEN
C           WRITE(*,'(A,I6,A,F20.10,A,F15.7,A,F12.4,A,F12.4)') 
C    1          ' Cycle ',ITER,' E=',ENERGY,' RMS=',RMS,' eigenvalue=',EVALMIN,' overlap=',SOVER
         ENDIF
C     ENDIF
C
      IF (BFGSSTEP) THEN
         FIXIMAGE=.FALSE.
         FIXD=FIXDSAVE
         RETURN
      ENDIF
C
C  Tangent space minimization next.
C  Uphill direction is projected out of the step in mylbfgs
C  The next IF block allows for zero tangent space steps in the initial phase
C
      IF ((NBFGSMAX1.EQ.0).AND.((1.0D0-DABS(SOVER).GT.BFGSTSTOL).OR.(INEG.EQ.0))) THEN
         FIXIMAGE=.FALSE.
         ITER=ITER+1
         FIXD=FIXDSAVE
         MFLAG=.FALSE.
         IF (ITER.GT.ITMAX) RETURN
         GOTO 90
      ENDIF

      IF (((1.0D0-DABS(SOVER).GT.BFGSTSTOL).OR.(HINDEX.NE.INEG)).OR.(INEG.EQ.0).OR.(TWOENDS.AND.(ITER.EQ.1)).OR.(ITER.EQ.1)) THEN
          NBFGS=NBFGSMAX1
      ELSE
          NBFGS=NBFGSMAX2
      ENDIF

!
! If we don't reset then the previous steps and gradients saved in mylbfgs
! will not necessarily be perpendicular to the new uphill direction.
!
      RESET=.TRUE.
      IF (ITER.EQ.1) RESET=.TRUE.
!     RMS2=RMS
      NSTEPMINSAVE=NSTEPMIN
      NSTEPMIN=0
      IF (CHRMMT.AND.INTMINT) THEN
         CALL MYLBFGS(NINTS,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,NBFGS,
     1                RESET,ITDONE,PTEST,GRAD,.TRUE.,.TRUE.)
      ELSE
         CALL MYLBFGS(NOPT,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,NBFGS,
     1                RESET,ITDONE,PTEST,GRAD,.TRUE.,.TRUE.)
      ENDIF
      NSTEPMIN=NSTEPMINSAVE
      IF (CFUSIONT) THEN
         MFLAG=.FALSE.
         RETURN
      ENDIF

      IF (REPELTST) CALL REPELSP(COORDS,RMS,INEG,MFLAG)

      RMS2=RMS ! save subspace RMS
      RMS=DSQRT(DOTOPT(GRAD(1),GRAD(1),NOPT))/SQRT(1.0D0*NOPT) ! true RMS
      IF (PTEST) WRITE(*,'(A,F15.7,A,F15.7,A,F15.7)') ' bfgsts> RMS grad=',RMS,
     &         ' subspace grad=',RMS2,' unscaled EF step length=',STPMAG
      IF (MFLAG) THEN
!        IF (NFREEZE.GT.0) THEN
!           PRINT '(A,I8,A)',' bfgsts> unfreezing ',NFREEZE,' atoms'
!           MFLAG=.FALSE.
!           FROZEN(1:NATOMS)=.FALSE.
!           NFREEZE=0
!        ENDIF
         IF (((RMS.GT.CONVR).OR.(INEG.EQ.0)).OR.(STPMAG.GT.CONVU)) MFLAG=.FALSE.
         IF (MFLAG) THEN
            FIXD=FIXDSAVE
            RETURN
         ENDIF
      ENDIF
!
! Bug fix 11/4/08 DJW - would have converged on subspace gradient, not total after subspace minimisation
! without projection. Probably caused when projection was changed in MYLBFGS.
!
      ITER=ITER+1
      IF (ITER.GT.ITMAX) THEN
         FIXD=FIXDSAVE
         RETURN
      ENDIF
C     IF ((AMBER).AND.(MOVIE)) CALL AMOVIEDUMP(FRAME)
      EOLD=ENERGY
      DO J1=1,NOPT
         XPSTEP(J1)=STEP(J1)
         PFOB(J1)=XFOB(J1)
      ENDDO
      GOTO 90
! }}}
      RETURN
      END
! }}}
C
!ITEG: Orthogonalise VEC1 to overall translations and rotations.{{{
C
      SUBROUTINE ORTHOGOPT(VEC1,COORDS,OTEST) 
! {{{
      USE COMMONS
      USE KEY
      USE VECCK
      IMPLICIT NONE
      INTEGER J2, J3, NCHECK
      DOUBLE PRECISION COORDS(*), VEC1(*), DUMMY1, DUMMY2, ROOTN, VDOT,
     1                 CMX, CMY, CMZ, AMASS(NATOMS), TMASS, RMASS(NATOMS)
      LOGICAL OTEST

      IF (RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')) THEN
         IF (OTEST) CALL VECNORM(VEC1,NOPT)
         RETURN
      ENDIF
      IF (ZSYM(NATOMS).EQ.'TH') THEN
         CALL ORTHOGTH(VEC1,COORDS,OTEST)
         RETURN
      ENDIF
      IF (STOCKT) THEN
         CALL ORTHOGSTOCK(VEC1,COORDS,OTEST)
         RETURN
      ENDIF
      IF (RBAAT) THEN
         CALL ORTHOGRIGID(VEC1,COORDS,OTEST)
         RETURN
      ENDIF 
      IF (.NOT.ALLOCATED(VECCHK)) ALLOCATE(VECCHK(3*NATOMS,30))
      ROOTN=SQRT(1.0D0*NATOMS)
      NCHECK=0
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      TMASS=0.0D0
      DO J2=1,NATOMS
         AMASS(J2)=1.0D0
         RMASS(J2)=1.0D0
         IF (MASST) AMASS(J2)=ATMASS(J2)
         IF (MASST) RMASS(J2)=SQRT(ATMASS(J2))
         TMASS=TMASS+AMASS(J2)
      ENDDO
C
C  If MASST then the coordinates already have a square root of the mass in them
C
      DO J2=1,NATOMS
         CMX=CMX+COORDS(3*(J2-1)+1)*RMASS(J2)
         CMY=CMY+COORDS(3*(J2-1)+2)*RMASS(J2)
         CMZ=CMZ+COORDS(3*(J2-1)+3)*RMASS(J2)
      ENDDO
      CMX=CMX/TMASS
      CMY=CMY/TMASS
      CMZ=CMZ/TMASS
C     WRITE(*,'(A,3F20.10)') 'centre at ',CMX,CMY,CMZ
1     VDOT=0.0D0
      NCHECK=NCHECK+1
C
C  Orthogonalize to known eigenvectors corresponding to negative eigenvalues
C  for CHECKINDEX run with NOHESS.
C
      IF (CHECKINDEX) THEN
         DO J2=1,NMDONE
            DUMMY1=0.0D0
            DO J3=1,NOPT
               DUMMY1=DUMMY1+VECCHK(J3,J2)*VEC1(J3)
            ENDDO
            VDOT=MAX(VDOT,ABS(DUMMY1))
            DUMMY2=0.0D0
            DO J3=1,NOPT
               VEC1(J3)=VEC1(J3)-DUMMY1*VECCHK(J3,J2)
               DUMMY2=DUMMY2+VEC1(J3)**2
            ENDDO
            DUMMY2=1.0D0/DSQRT(DUMMY2)
            DO J3=1,NOPT
               VEC1(J3)=VEC1(J3)*DUMMY2
            ENDDO
         ENDDO
         IF (OTEST) CALL VECNORM(VEC1,NOPT)
      ENDIF
C
      IF ((VDOT.GT.1.0D-6).AND.(SHIFTED).AND.(NCHECK.LT.100)) GOTO 1
      IF (NCHECK.GE.100) THEN
         PRINT*,'*** WARNING, cannot orthogonalise to known eigenvectors in ORTHOGOPT'
      ENDIF
      IF ((SHIFTED).AND.(.NOT.FIXD)) RETURN  ! zeros are shifted - lowest eigenvector should be uncontaminated

      IF (ZSYM(1).EQ.'BE') GOTO 10
      IF (EYTRAPT) GOTO 10
      IF (RTEST.AND.(JZ.EQ.0.0D0)) GOTO 20
      IF (RTEST.AND.(JZ.NE.0.0D0)) GOTO 30

      DUMMY1=0.0D0
      DO J2=1,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/ROOTN
      VDOT=MAX(VDOT,ABS(DUMMY1))
C     WRITE(*,'(A,F20.10)') 'X dot:',DUMMY1
      DO J2=1,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1/ROOTN
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      DUMMY1=0.0D0
      DO J2=2,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/ROOTN
      VDOT=MAX(VDOT,ABS(DUMMY1))
C     WRITE(*,'(A,F20.10)') 'Y dot:',DUMMY1
      DO J2=2,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1/ROOTN
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

30    DUMMY1=0.0D0
      DO J2=3,NOPT,3
         DUMMY1=DUMMY1+VEC1(J2)
      ENDDO
      DUMMY1=DUMMY1/ROOTN
      VDOT=MAX(VDOT,ABS(DUMMY1))
C     WRITE(*,'(A,F20.10)') 'Z dot:',DUMMY1
      DO J2=3,NOPT,3
         VEC1(J2)=VEC1(J2)-DUMMY1/ROOTN
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

31    CONTINUE
      IF ((VDOT.GT.1.0D-6).AND.BULKT.AND.(NCHECK.LT.100)) GOTO 1
!     PRINT*,'after next part VDOT=',VDOT
      IF (NCHECK.GE.100) THEN
         PRINT*,'*** WARNING, cannot orthogonalise to known eigenvectors in ORTHOGOPT'
      ENDIF
      IF (BULKT.AND.TWOD) GOTO 20
      IF (BULKT) RETURN
      IF (DFTP) RETURN
      IF (RTEST.AND.(JZ.NE.0.0D0)) GOTO 20
      IF (PULLT.OR.EFIELDT) GOTO 20

10    DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         J3=3*J2
         DUMMY1=DUMMY1+RMASS(J2)*(VEC1(J3-1)*(COORDS(J3)-CMZ)-VEC1(J3)*(COORDS(J3-1)-CMY))
         DUMMY2=DUMMY2+AMASS(J2)*((COORDS(J3)-CMZ)**2+(COORDS(J3-1)-CMY)**2)
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
!        WRITE(*,'(A,F20.10)') 'XRdot:',DUMMY1/SQRT(DUMMY2)
         VDOT=MAX(VDOT,ABS(DUMMY1)/SQRT(DUMMY2))
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NATOMS
            J3=3*J2
            VEC1(J3-1)=VEC1(J3-1)-DUMMY2*RMASS(J2)*(COORDS(J3)-CMZ)
            VEC1(J3)=VEC1(J3)+DUMMY2*RMASS(J2)*(COORDS(J3-1)-CMY)
         ENDDO
      ENDIF
      IF (OTEST) CALL VECNORM(VEC1,NOPT)
!     PRINT*,'rot x component=',DUMMY2

      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         J3=3*J2
         DUMMY1=DUMMY1-VEC1(J3-2)*RMASS(J2)*(COORDS(J3)-CMZ)+VEC1(J3)*RMASS(J2)*(COORDS(J3-2)-CMX)
         DUMMY2=DUMMY2+AMASS(J2)*((COORDS(J3)-CMZ)**2+(COORDS(J3-2)-CMX)**2)
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
!        WRITE(*,'(A,F20.10)') 'YRdot:',DUMMY1/SQRT(DUMMY2)
         VDOT=MAX(VDOT,ABS(DUMMY1)/SQRT(DUMMY2))
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NATOMS
            J3=3*J2
            VEC1(J3-2)=VEC1(J3-2)+DUMMY2*RMASS(J2)*(COORDS(J3)-CMZ)
            VEC1(J3)=VEC1(J3)-DUMMY2*RMASS(J2)*(COORDS(J3-2)-CMX)
         ENDDO
      ENDIF
      IF (OTEST) CALL VECNORM(VEC1,NOPT)
!     PRINT*,'rot y component=',DUMMY2

20    DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         J3=3*J2
         DUMMY1=DUMMY1+VEC1(J3-2)*RMASS(J2)*(COORDS(J3-1)-CMY)-VEC1(J3-1)*RMASS(J2)*(COORDS(J3-2)-CMX)
         DUMMY2=DUMMY2+AMASS(J2)*((COORDS(J3-1)-CMY)**2+(COORDS(J3-2)-CMX)**2)
      ENDDO
      IF (DUMMY2.GT.0.0D0) THEN
!        WRITE(*,'(A,F20.10)') 'ZRdot:',DUMMY1/SQRT(DUMMY2)
         VDOT=MAX(VDOT,ABS(DUMMY1)/SQRT(DUMMY2))
         DUMMY2=DUMMY1/DUMMY2
         DO J2=1,NATOMS
            J3=3*J2
            VEC1(J3-2)=VEC1(J3-2)-DUMMY2*RMASS(J2)*(COORDS(J3-1)-CMY)
            VEC1(J3-1)=VEC1(J3-1)+DUMMY2*RMASS(J2)*(COORDS(J3-2)-CMX)
         ENDDO
      ENDIF
      IF (OTEST) CALL VECNORM(VEC1,NOPT)
!     PRINT*,'rot z component=',DUMMY2

C     WRITE(*,'(A,F20.10)') 'Largest remaining component in ORTHOGOPT=',VDOT
      IF ((VDOT.GT.1.0D-6).AND.(NCHECK.LT.100)) GOTO 1
C     PRINT*,'after next part VDOT=',VDOT
      IF (NCHECK.GE.100) THEN
         PRINT*,'*** WARNING, cannot orthogonalise to known eigenvectors in ORTHOGOPT'
      ENDIF
! }}}
      RETURN
      END
! }}}

!ITEIG: Est. smallest non-zero eigenvalue & eigenvector by repeated iteration {{{
C
C  Estimate the smallest non-zero eigenvalue and the associated
C  eigenvector by repeated iteration. Modified to find more than
C  one eigenvector. The eigenvalue of the chosen eigenvector is
C  returned as EVALMIN, though it may not be the minimum value.
C  VECS contains the chosen eigenvector, and is not altered by the
C  calling routine.
C
      SUBROUTINE ITEIG(ITER,COORDS,VECS,EVALMIN,EVALMAX,NS,SOVER,PTEST,VECL,CONVERGED)
! Doxygen {{{
!>
!> \name ITEIG
!> \brief Estimate the smallest non-zero eigenvalue and the associated
!> \brief eigenvector by repeated iteration. 
!>  Modified to find more than
!>  one eigenvector. The eigenvalue of the chosen eigenvector is
!>  returned as EVALMIN, though it may not be the minimum value.
!>  VECS contains the chosen eigenvector, and is not altered by the
!>  calling routine.
!> \param[in]
!> \param[in]
!> \param[out]
!> \param[out]
!>
! }}}
! Modules {{{

        USE COMMONS
        USE KEY
        USE MODTWOEND
        USE VECCK
        USE MODHESS
! }}}

      IMPLICIT NONE

! subroutine parameters {{{

      INTEGER ITER 
      DOUBLE PRECISION, DIMENSION(3*NATOMS) :: COORDS, VECS 
      DOUBLE PRECISION EVALMIN, EVALMAX
      INTEGER NS  
      DOUBLE PRECISION SOVER
      LOGICAL PTEST
      DOUBLE PRECISION,DIMENSION(3*NATOMS) :: VECL
      LOGICAL CONVERGED

! }}}
! local parameters {{{

      INTEGER J1, J2, ISEED
      DOUBLE PRECISION DUMMY1,VEC1(3*NATOMS),XSIGN,
     1                 VEC2(3*NATOMS),DUMMY2,DUMMY, DUMMY4,
     2                 PERCENT,OVLAP(30), OVMAX, DPRAND, EVALS(30)
      INTEGER JVEC, NDUM, NMPREV, JVECP
      COMMON /IS/ ISEED
      SAVE

! }}}
! subroutine body {{{

      NMPREV=NMDONE
      IF (ITER.EQ.1) JVECP=30
C
C  First estimate the largest eigenvalue of HESS. 
C
      IF ((ITER.EQ.1).AND.(.NOT.TTDONE)) THEN
         DO J1=1,NOPT
            VEC1(J1)=2*(DPRAND()-0.5D0)
         ENDDO
         CALL VECNORM(VEC1,NOPT)

C        DO J1=1,NOPT
C           VEC1(J1)=0.0D0
C        ENDDO
C        VEC1(1)=1.0D0

      ELSE
         DO J1=1,NOPT
            VEC1(J1)=VECL(J1)
         ENDDO
      ENDIF
      IF (FREEZE) THEN
         DO J1=1,NATOMS
            IF (FROZEN(J1)) THEN
               VEC1(3*(J1-1)+1)=0.0D0
               VEC1(3*(J1-1)+2)=0.0D0
               VEC1(3*(J1-1)+3)=0.0D0
            ENDIF
         ENDDO
         CALL VECNORM(VEC1,NOPT)
      ENDIF

      DUMMY1=-1.0D10
      DO J1=1,NEVL
         IF (.NOT.FREEZE) THEN
            CALL ORTHOGOPT(VEC1,COORDS,.TRUE.)
         ELSE
            CALL VECNORM(VEC1,NOPT)
         ENDIF
         CALL DSYMV('U',NOPT,1.0D0,HESS,SIZE(HESS,1),VEC1,1,0.0D0,VEC2,1)
C        CALL MATMUL(VEC2,HESS,VEC1,3*NATOMS,3*NATOMS,1,NOPT,NOPT,1)
         DUMMY2=0.0D0
         DO J2=1,NOPT
            DUMMY2=DUMMY2+VEC2(J2)**2
         ENDDO
         EVALMAX=DSQRT(DUMMY2)
         DO J2=1,NOPT
            VEC2(J2)=VEC2(J2)/EVALMAX
         ENDDO
         PERCENT=100.0D0*DABS((DUMMY1-EVALMAX)/DUMMY1)
C        CALL EXTRAP(J1,NOPT,VEC2,PERCENT)
C        WRITE(*,'(A,4F20.10)') 'DUMMY1,EVALMAX,percent,CEIG=',DUMMY1,EVALMAX,PERCENT,CEIG
         IF (100.0D0*DABS((DUMMY1-EVALMAX)/DUMMY1).LT.CEIG) THEN
            IF (PTEST) WRITE(*,'(A,I4,A,F15.7)') ' Largest  eigenvalue converged in ',J1,' steps. Eigenvalue=',EVALMAX
            GOTO 10
         ENDIF
C        WRITE(*,'(A,I4,A,F15.7,A,F15.7)') ' At step ',J1,' % change=',PERCENT,' eigenvalue=',EVALMAX
         DUMMY1=EVALMAX

         DO J2=1,NOPT
            VEC1(J2)=VEC2(J2)
         ENDDO
      ENDDO

      IF (PTEST) WRITE(*,'(A,F15.7)') ' ****WARNING - Largest  eigenvalue did not converge, value=',EVALMAX
10    IF (ITER.GT.1) THEN
         DUMMY1=0.0D0
         DO J2=1,NOPT
            DUMMY1=DUMMY1+VECL(J2)*VEC2(J2)
         ENDDO
         IF (PTEST) WRITE(*,'(A,G12.4)') 
     1          ' Overlap between previous and new eigenvectors for largest  eigenvalue=',DUMMY1
      ENDIF
      DO J2=1,NOPT
         VECL(J2)=VEC2(J2)
      ENDDO
C
C  Shift all the eigenvalues according to the size of EVALMAX.
C
      DO J1=1,NATOMS
         J2=3*(J1-1)
C        IF (.NOT.FROZEN(J1)) THEN
            HESS(J2+1,J2+1)=HESS(J2+1,J2+1)-EVALMAX*0.55D0
            HESS(J2+2,J2+2)=HESS(J2+2,J2+2)-EVALMAX*0.55D0
            HESS(J2+3,J2+3)=HESS(J2+3,J2+3)-EVALMAX*0.55D0
C        ENDIF
      ENDDO
C
C  Now estimate the new largest magnitude eigenvalue of HESS.
C
      IF ((ITER.EQ.1).AND.(.NOT.TWOENDS).AND.(.NOT.READV)) THEN
         DO J1=1,NOPT
            VEC1(J1)=2*(DPRAND()-0.5D0)
         ENDDO
         CALL VECNORM(VEC1,NOPT)
C        DO J1=1,NOPT
C           VEC1(J1)=0.0D0
C        ENDDO
C        VEC1(1)=1.0D0
      ELSE
         DO J1=1,NOPT
            VEC1(J1)=VECS(J1)
         ENDDO
      ENDIF
      IF (FREEZE) THEN
         DO J1=1,NATOMS
            IF (FROZEN(J1)) THEN
               VEC1(3*(J1-1)+1)=0.0D0
               VEC1(3*(J1-1)+2)=0.0D0
               VEC1(3*(J1-1)+3)=0.0D0
            ENDIF
         ENDDO
         CALL VECNORM(VEC1,NOPT)
      ENDIF
      DUMMY1=-1.0D10
      DO J1=1,NEVS
         IF (.NOT.FREEZE) THEN
            CALL ORTHOGOPT(VEC1,COORDS,.TRUE.)
         ELSE
            CALL VECNORM(VEC1,NOPT)
         ENDIF
         CALL DSYMV('U',NOPT,1.0D0,HESS,SIZE(HESS,1),VEC1,1,0.0D0,VEC2,1)

         DUMMY2=0.0D0
         DO J2=1,NOPT
            DUMMY2=DUMMY2+VEC2(J2)**2
         ENDDO
         EVALMIN=DSQRT(DUMMY2)
C
C  Normalize. Changing the phase may be necessary for extrapolation to work.
C
         DO J2=1,NOPT
            VEC2(J2)=-VEC2(J2)/EVALMIN
         ENDDO
         EVALMIN=-EVALMIN+EVALMAX*0.55D0
         PERCENT=100.0D0*DABS((DUMMY1-EVALMIN)/DUMMY1)
C        CALL EXTRAP(J1,NOPT,VEC2,PERCENT)
         IF (PERCENT.LT.CEIG) THEN
            IF (PTEST) WRITE(*,'(A,I4,A,F15.7)') ' Smallest eigenvalue converged in ',J1,' steps. Eigenvalue=',EVALMIN
            CONVERGED=.TRUE.
            GOTO 20
         ENDIF
C        WRITE(*,'(A,I4,A,F15.7,A,F15.7)') ' At step ',J1,' % change=',PERCENT,' eigenvalue=',EVALMIN
         DUMMY1=EVALMIN

         DO J2=1,NOPT
            VEC1(J2)=VEC2(J2)
         ENDDO
      ENDDO
      IF (PTEST) WRITE(*,'(A,F15.7)') ' **WARNING - Smallest eigenvalue did not converge, value=',EVALMIN
      CONVERGED=.FALSE.
     
C 20    EVALMIN=-EVALMIN+EVALMAX*0.55D0
20    CONTINUE
      NS=J1
      IF (ITER.GT.1) THEN
         SOVER=0.0D0
         DO J1=1,NOPT
            SOVER=SOVER+VECS(J1)*VEC2(J1)
         ENDDO
         IF (PTEST) WRITE(*,'(A,G12.4)') ' Overlap between previous and new eigenvectors for smallest eigenvalue=',SOVER
      ENDIF
C
C  The following code should allow us to search along other eigenvectors.
C
      IF ((IVEC.EQ.0).OR.((ITER.EQ.1).AND.(ABS(IVEC).EQ.1))) THEN
         DO J1=1,NOPT
            VECS(J1)=VEC2(J1)
         ENDDO
         JVECP=1
         RETURN
      ENDIF

      EVALS(1)=EVALMIN
      DO J1=1,NOPT
         VECCHK(J1,1)=VEC2(J1)
      ENDDO
      IF (ITER.GT.1) THEN
         DUMMY=0.0D0
         DO J1=1,NOPT
            DUMMY=DUMMY+VECS(J1)*VEC2(J1)
         ENDDO
         OVLAP(1)=ABS(DUMMY)
         IF (ABS(DUMMY).GT.0.8D0) THEN
            DO J1=1,NOPT
               VECS(J1)=VEC2(J1)
            ENDDO
            IF (PTEST) WRITE(*,'(A,F12.5,A,F12.5)') ' Mode    1 will be searched uphill, overlap=',OVLAP(1),' eigenvalue=',EVALMIN
            JVECP=1
            RETURN
         ENDIF
      ENDIF
      NMDONE=1
      JVEC=1
      OVMAX=OVLAP(1)
C
C  The largest eigenvalue in magnitude of the shifted Hessian is 
C  guaranteed to be negative. However, if we need to find other
C  (shifted) eigenvalues we must allow for the possibility that they could
C  be positive.
C
      CHECKINDEX=.TRUE.
80    CONTINUE
      IF ((ITER.GT.1).AND.(NMDONE+1.LE.NMPREV)) THEN
         DO J1=1,NOPT
            VEC1(J1)=VECCHK(J1,NMDONE+1)
         ENDDO
      ELSE
         DO J1=1,NOPT
            VEC1(J1)=2*(DPRAND()-0.5D0)
         ENDDO
         CALL VECNORM(VEC1,NOPT)
      ENDIF
      IF (FREEZE) THEN
         DO J1=1,NATOMS
            IF (FROZEN(J1)) THEN
               VEC1(3*(J1-1)+1)=0.0D0
               VEC1(3*(J1-1)+2)=0.0D0
               VEC1(3*(J1-1)+3)=0.0D0
            ENDIF
         ENDDO
         CALL VECNORM(VEC1,NOPT)
      ENDIF

      DUMMY4=-1.0D10
      DO J1=1,NEVS
         IF (.NOT.FREEZE) THEN
            CALL ORTHOGOPT(VEC1,COORDS,.TRUE.)
         ELSE
            CALL VECNORM(VEC1,NOPT)
         ENDIF

         CALL DSYMV('U',NOPT,1.0D0,HESS,SIZE(HESS,1),VEC1,1,0.0D0,VEC2,1)
C        CALL MATMUL(VEC2,HESS,VEC1,3*NATOMS,3*NATOMS,1,NOPT,NOPT,1)

         DUMMY2=0.0D0
         DUMMY1=0.0D0
         NDUM=0
         DO J2=1,NOPT
            DUMMY2=DUMMY2+VEC2(J2)**2
            IF (VEC1(J2).NE.0.0D0) THEN
               DUMMY1=DUMMY1+VEC2(J2)/VEC1(J2)
               NDUM=NDUM+1
            ENDIF
         ENDDO
         EVALMIN=DSQRT(DUMMY2)
         DUMMY1=DUMMY1/MAX(NDUM,1)
         XSIGN=-1.0D0
         IF (DABS((DABS(DUMMY1)-DABS(EVALMIN))/EVALMIN).LT.0.1D0) XSIGN=DUMMY1/DABS(DUMMY1)
         EVALMIN=EVALMIN*XSIGN
C
C  Normalize.
C
         DO J2=1,NOPT
            VEC2(J2)=VEC2(J2)/EVALMIN
         ENDDO

         EVALMIN=EVALMIN+EVALMAX*0.55D0
         PERCENT=100.0D0*DABS((DUMMY4-EVALMIN)/DUMMY4)
         IF (PERCENT.LT.CEIG) THEN
            IF (PTEST) WRITE(*,'(A,I4,A,F15.7)') ' New eigenvalue converged in      ',J1,' steps. Eigenvalue=',EVALMIN
            GOTO 90
         ENDIF
         DUMMY4=EVALMIN

         DO J2=1,NOPT
            VEC1(J2)=VEC2(J2)
         ENDDO
      ENDDO
      IF (PTEST) WRITE(*,'(A,F15.7)') ' ****WARNING - New eigenvalue  did not converge, value=  ',EVALMIN

90    CONTINUE

      NMDONE=NMDONE+1

      EVALS(NMDONE)=EVALMIN
      DO J2=1,NOPT
         VECCHK(J2,NMDONE)=VEC2(J2)
      ENDDO

      IF (ITER.GT.1) THEN
         DUMMY=0.0D0
         DO J1=1,NOPT
            DUMMY=DUMMY+VECS(J1)*VEC2(J1)
         ENDDO
         OVLAP(NMDONE)=ABS(DUMMY)
C        WRITE(*,'(A,I4,F15.5)') 'NMDONE,OVLAP=',NMDONE,OVLAP(NMDONE)
         IF (OVLAP(NMDONE).GT.OVMAX) THEN
            JVEC=NMDONE
            OVMAX=OVLAP(NMDONE)
         ENDIF
         IF (ABS(DUMMY).GT.0.8D0) THEN
            DO J1=1,NOPT
               VECS(J1)=VEC2(J1)
            ENDDO
            IF (PTEST) WRITE(*,'(A,I4,A,F12.5,A,F12.5)') ' Mode ',NMDONE,' will be searched uphill, overlap=',OVLAP(NMDONE),
     1                                ' eigenvalue=',EVALS(NMDONE)
            CHECKINDEX=.FALSE.
            JVECP=NMDONE
            RETURN
         ENDIF
      ELSE
         IF (NMDONE.EQ.ABS(IVEC)) THEN
            DO J1=1,NOPT
               VECS(J1)=VEC2(J1)
            ENDDO
            WRITE(*,'(A,I4,A,F12.5)') 'Mode ',NMDONE,' will be searched uphill, eigenvalue=',EVALS(NMDONE)
            CHECKINDEX=.FALSE.
            JVECP=NMDONE
            RETURN
         ENDIF
      ENDIF

C     IF (NMDONE.EQ.30) THEN
      IF ((NMDONE.EQ.30).OR.(NMDONE.EQ.JVECP+3)) THEN
         DO J1=1,NOPT
            VECS(J1)=VECCHK(J1,JVEC)
         ENDDO
         EVALMIN=EVALS(JVEC)
         IF (PTEST) WRITE(*,'(I3,A,I3,A,F12.5,A,F12.5)') NMDONE,' eigenvectors found; using number ',JVEC,
     1                                               ' eigenvalue=',EVALS(JVEC),' overlap=',OVMAX
         CHECKINDEX=.FALSE.
         JVECP=JVEC
         RETURN
      ENDIF

      GOTO 80
! }}}
      RETURN
      END
! }}}
C
C  Estimate the smallest non-zero eigenvalue and the associated
C  eigenvector by repeated iteration using Nesbet's algorithm.
C  This seems to give the largest not the smallest eigenvector?
C     NOT USED
C
      SUBROUTINE ITEIGN(ITER,COORDS,VECS,EVALMIN,EVALMAX,PTEST)
      USE COMMONS
      USE KEY
      USE VECCK
      USE MODHESS
      IMPLICIT NONE

      INTEGER J1, J2, J3
      DOUBLE PRECISION COORDS(3*NATOMS),DUMMY1,VEC1(3*NATOMS),
     1                 VEC2(3*NATOMS),DUMMY2,EVALMAX,EVALMIN,
     2                 VECS(3*NATOMS),DELC,VECL(3*NATOMS)
      INTEGER ITER
      LOGICAL PTEST
      SAVE

      IF (ITER.EQ.1) THEN
         DO J1=1,NOPT
            VEC1(J1)=0.0D0
         ENDDO
         VEC1(1)=1.0D0
      ELSE
         DO J1=1,NOPT
            VEC1(J1)=VECL(J1)
         ENDDO
      ENDIF

      CALL DSYMV('U',NOPT,1.0D0,HESS,SIZE(HESS,1),VEC1,1,0.0D0,VEC2,1)

      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J1=1,NOPT
         DUMMY1=DUMMY1+VEC1(J1)*VEC2(J1)
         DUMMY2=DUMMY2+VEC1(J1)**2
      ENDDO
      EVALMAX=DUMMY1/DUMMY2
      DUMMY1=EVALMAX
      
      DO J1=1,NEVL
         DO J2=2,NOPT
            DELC=(VEC2(J2)-DUMMY1*VEC1(J2))/(DUMMY1-HESS(J2,J2))
            VEC1(J2)=VEC1(J2)+DELC

            CALL DSYMV('U',NOPT,1.0D0,HESS,SIZE(HESS,1),VEC1,1,0.0D0,VEC2,1)
   
            DUMMY1=0.0D0
            DUMMY2=0.0D0
            DO J3=1,NOPT
               DUMMY1=DUMMY1+VEC1(J3)*VEC2(J3)
               DUMMY2=DUMMY2+VEC1(J3)**2
            ENDDO
            DUMMY1=DUMMY1/DUMMY2
            DUMMY2=1.0D0/DSQRT(DUMMY2)
         ENDDO
         IF ((100.0D0*DABS((DUMMY1-EVALMAX)/DUMMY1).LT.CEIG).AND.(J1.GT.1)) THEN
            IF (PTEST) WRITE(*,'(A,I4,A,F15.7)') ' Largest eigenvalue converged in ',J1,' steps. Eigenvalue=',DUMMY1
            GOTO 20
         ENDIF
         IF (PTEST) WRITE(*,'(A,I4,A,F15.7,A,F15.7)') ' At step ',J1,' convergence=',DABS(DUMMY1-EVALMAX),' eigenvalue=',DUMMY1
         EVALMAX=DUMMY1
      ENDDO
      IF (PTEST) WRITE(*,'(A,F15.7)') ' ****WARNING - Largest eigenvalue did not converge, value=',EVALMAX
     
20    CONTINUE
      DO J2=1,NOPT
         VECL(J2)=VEC1(J2)
      ENDDO
C
C  Shift all the eigenvalues according to the size of EVALMAX.
C
      DO J1=1,NOPT
         HESS(J1,J1)=HESS(J1,J1)-EVALMAX*0.55D0
      ENDDO
C
C  Now estimate the new largest magnitude eigenvalue of HESS.
C
      IF (ITER.EQ.1) THEN
         DO J1=1,NOPT
            VEC1(J1)=0.0D0
         ENDDO
         VEC1(1)=1.0D0
      ELSE
         DO J1=1,NOPT
            VEC1(J1)=VECS(J1)
         ENDDO
      ENDIF

      CALL DSYMV('U',NOPT,1.0D0,HESS,SIZE(HESS,1),VEC1,1,0.0D0,VEC2,1)

      DUMMY1=0.0D0
      DUMMY2=0.0D0
      DO J1=1,NOPT
         DUMMY1=DUMMY1+VEC1(J1)*VEC2(J1)
         DUMMY2=DUMMY2+VEC1(J1)**2
      ENDDO
      EVALMIN=DUMMY1/DUMMY2
      DUMMY1=EVALMIN

      DO J1=1,NEVS
         DO J2=1,NOPT-1
            DELC=(VEC2(J2)-DUMMY1*VEC1(J2))/(DUMMY1-HESS(J2,J2))
            VEC1(J2)=VEC1(J2)+DELC

            CALL DSYMV('U',NOPT,1.0D0,HESS,SIZE(HESS,1),VEC1,1,0.0D0,VEC2,1)
   
            DUMMY1=0.0D0
            DUMMY2=0.0D0
            DO J3=1,NOPT
               DUMMY1=DUMMY1+VEC1(J3)*VEC2(J3)
               DUMMY2=DUMMY2+VEC1(J3)**2
            ENDDO

            DUMMY1=DUMMY1/DUMMY2

            DUMMY2=1.0D0/DSQRT(DUMMY2)
            DO J3=1,NOPT
               VEC1(J3)=VEC1(J3)*DUMMY2
            ENDDO
            CALL ORTHOGOPT(VEC1,COORDS,.TRUE.)
         ENDDO
         IF (100.0D0*DABS((DUMMY1-EVALMIN)/DUMMY1).LT.CEIG) THEN
            IF (PTEST) WRITE(*,'(A,I4,A,F15.7)') ' Smallest eigenvalue converged in ',J1,' steps. Eigenvalue=', EVALMIN+EVALMAX*0.55D0
            GOTO 30
         ENDIF
         IF (PTEST) WRITE(*,'(A,I4,A,F15.7,A,F15.7)') ' At step ',J1,' convergence=',DABS(DUMMY1-EVALMIN),
     1           ' eigenvalue=', EVALMIN+EVALMAX*0.55D0
         EVALMIN=DUMMY1
      ENDDO
      IF (PTEST) WRITE(*,'(A,F15.7)') ' **WARNING - Smallest eigenvalue did not converge, value=', EVALMIN+EVALMAX*0.55D0
     
30    EVALMIN= EVALMIN+EVALMAX*0.55D0
      DUMMY1=0.0D0
      DO J2=1,NOPT
         DUMMY1=DUMMY1+VEC1(J2)**2
      ENDDO
      DUMMY1=1.0D0/SQRT(DUMMY1)
      DO J2=1,NOPT
         VECS(J2)=VEC1(J2)*DUMMY1
      ENDDO

      RETURN
      END
C
C  Count the number of negative Hessian eigenvalues by further iteration and
C  orthogonalization if required.
C
      SUBROUTINE CHECKIND(COORDS,MFLAG,INEG,ENERGY,EVALMIN,EVALMAX,PTEST)
      USE COMMONS
      USE KEY
      USE VECCK
      USE ZWK
      USE MODHESS
      IMPLICIT NONE

      INTEGER J1, J2, INEG, ISEED, NDUM, NS
      DOUBLE PRECISION COORDS(3*NATOMS),DUMMY1,VEC1(3*NATOMS),
     1                 VEC2(3*NATOMS),DUMMY2,EVALMAX,EVALMIN,
     2                 ENERGY,DPRAND,XSIGN,VECS(3*NATOMS),
     3                 DUMMY4,SOVER,VECL(3*NATOMS)
      LOGICAL MFLAG,PTEST,CONVERGED
      COMMON /IS/ ISEED
C
C  Assign enough memory to WORK for a blocksize of 32 to be possible.
C  This is for DSYEVR.
C
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NATOMS)
      INTEGER IWORK(33*3*NATOMS), INFO
      DOUBLE PRECISION WORK(33*3*NATOMS),  ABSTOL, DIAG(3*NATOMS), DLAMCH

      FIXIMAGE=.TRUE.
      LWORK=33*3*NATOMS
      ILWORK=33*3*NATOMS

      IF (NOIT) THEN
         CALL POTENTIAL(COORDS,ENERGY,VEC1,.TRUE.,.TRUE.,DUMMY1,PTEST,.FALSE.) 
         IF (.NOT.VARIABLES) THEN
            IF (ZSYM(NATOMS).EQ.'TH') THEN
               CALL SHIFTHTH(COORDS,NATOMS)
            ELSEIF (RBAAT) THEN
               CALL SHIFTRIGID(COORDS,NATOMS)
            ELSE
               CALL SHIFTH(COORDS,.TRUE.,NOPT,NATOMS,ATMASS)
            ENDIF
         ENDIF
         ABSTOL=DLAMCH('Safe  minimum')
         CALL DSYEVR('V','I','U',NOPT,HESS,SIZE(HESS,1),0.0D0,1.0D0,1,10,ABSTOL,NFOUND,DIAG,ZWORK,3*NATOMS,ISUPPZ,WORK,
     1                        LWORK, IWORK, ILWORK, INFO )
         IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEVR'
C        PRINT*,'Optimal and actual values of LWORK=',WORK(1),LWORK
C        PRINT*,'Optimal and actual values of ILWORK=',IWORK(1),ILWORK
         EVALMIN=DIAG(1)
         SOVER=0.0D0
         WRITE(*,'(A)') 'Lowest 10 non-zero Hessian eigenvalues:'
         WRITE(*,'(5F15.5)') (DIAG(J1),J1=1,10)
         DO J1=1,10
            DO J2=1,NOPT
               VECCHK(J2,J1)=ZWORK(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NOPT
            VECS(J1)=ZWORK(J1,1)
         ENDDO
C
C  Get the largest and smallest eigenvalues if we don't already know them.
C
         IF (BFGSMINT.OR.RKMIN.OR.BSMIN) THEN
            CALL ITEIG(1,COORDS,VECS,EVALMIN,EVALMAX,NS,SOVER,PTEST,VECL,CONVERGED)
            IF (.NOT.CONVERGED) WRITE(*,'(A,F15.7)') ' **WARNING - Smallest eigenvalue did not converge, value=',EVALMIN
            IF (EVALMIN.LT.0.0D0) THEN
               INEG=1
            ELSE
               INEG=0
            ENDIF
         ENDIF
         NMDONE=INEG
         IF (EVALMIN.GT.0.0D0) THEN
            WRITE(*,'(A,F20.10)') ' Hessian index is 0, smallest eigenvector=',EVALMIN
            GOTO 666
         ENDIF
   
         IF (INEG.EQ.1) THEN
            DO J1=1,NOPT
               VECCHK(J1,1)=VECS(J1)
            ENDDO
         ENDIF
C
C  The largest eigenvalue in magnitude of the shifted Hessian is 
C  guaranteed to be negative. However, if we need to find other
C  (shifted) eigenvalues we must allow for the possibility that they could
C  be positive.
C

80       NMDONE=INEG
         DO J1=1,NOPT
            VEC1(J1)=2*(DPRAND()-0.5D0)
         ENDDO
         CALL VECNORM(VEC1,NOPT)
   
         DUMMY4=-1.0D10
         DO J1=1,NEVS
            CALL ORTHOGOPT(VEC1,COORDS,.TRUE.)

            CALL DSYMV('U',NOPT,1.0D0,HESS,SIZE(HESS,1),VEC1,1,0.0D0,VEC2,1)
C           CALL MATMUL(VEC2,HESS,VEC1,3*NATOMS,3*NATOMS,1,NOPT,NOPT,1)

            DUMMY2=0.0D0
            DUMMY1=0.0D0
            NDUM=0
            DO J2=1,NOPT
               DUMMY2=DUMMY2+VEC2(J2)**2
               IF (VEC1(J2).NE.0.0D0) THEN
                  DUMMY1=DUMMY1+VEC2(J2)/VEC1(J2)
                  NDUM=NDUM+1
               ENDIF
            ENDDO
            EVALMIN=DSQRT(DUMMY2)
            DUMMY1=DUMMY1/MAX(NDUM,1)
            XSIGN=-1.0D0
            IF (DABS((DABS(DUMMY1)-DABS(EVALMIN))/EVALMIN).LT.0.1D0) XSIGN=DUMMY1/DABS(DUMMY1)
            EVALMIN=EVALMIN*XSIGN
C
C  Normalize.
C
            DO J2=1,NOPT
               VEC2(J2)=VEC2(J2)/EVALMIN
            ENDDO
            IF (100.0D0*DABS((DUMMY4-EVALMIN)/DUMMY4).LT.CEIG) THEN
               IF (PTEST) WRITE(*,'(A,I4,A,F15.7)') ' New eigenvector converged in ',J1,' steps. Eigenvalue=',EVALMIN+EVALMAX*0.55D0
               GOTO 90
            ENDIF
C           WRITE(*,'(A,I4,A,F15.7,A,F15.7)') ' At step ',J1,' % change=',100.0D0*DABS((DUMMY4-EVALMIN)/DUMMY4),
C    1              ' eigenvalue=',EVALMIN+EVALMAX*0.55D0
            DUMMY4=EVALMIN

            DO J2=1,NOPT
               VEC1(J2)=VEC2(J2)
            ENDDO
         ENDDO
         IF (PTEST) WRITE(*,'(A,F15.7)') ' ****WARNING - New eigenvector did not converge, value=',EVALMIN+EVALMAX*0.55D0

90       EVALMIN=EVALMIN+EVALMAX*0.55D0
         IF (EVALMIN.GT.0.0D0) THEN
            PRINT*,'Hessian index=',INEG
         ELSE
            INEG=INEG+1
            IF (INEG.GT.10) THEN
               PRINT*,'Hessian index is 10 or greater'
               FIXIMAGE=.FALSE.
               RETURN
            ENDIF
            DO J2=1,NOPT
               VECCHK(J2,INEG)=VEC2(J2)
            ENDDO
            GOTO 80
         ENDIF
      ENDIF

666   CONTINUE

      IF (((BFGSMINT).AND.(INEG.NE.0)).AND.CHECKCONT) THEN
         MFLAG=.FALSE.
         IF (IVEC.GE.0) PRINT*,'Stepping away from saddle along mode ',1,' + direction'
         IF (IVEC.LT.0) PRINT*,'Stepping away from saddle along mode ',1,' - direction'
         DUMMY1=PUSHOFF
         IF (IVEC.LT.0) DUMMY1=-DUMMY1
         DO J1=1,NOPT
            COORDS(J1)=COORDS(J1)+DUMMY1*VECCHK(J1,1)
         ENDDO
         INEG=0
         NMDONE=0
      ENDIF
      IF (((BFGSTST).AND.(INEG.NE.1)).AND.CHECKCONT) THEN
         MFLAG=.FALSE.
         IF (IVEC.GE.0) PRINT*,'Stepping away from saddle along mode ',1,' + direction'
         IF (IVEC.LT.0) PRINT*,'Stepping away from saddle along mode ',1,' - direction'
         DUMMY1=PUSHOFF
         IF (IVEC.LT.0) DUMMY1=-DUMMY1
         DO J1=1,NOPT
            COORDS(J1)=COORDS(J1)+DUMMY1*VECCHK(J1,2)
         ENDDO
         INEG=0
         NMDONE=0
      ENDIF

      FIXIMAGE=.FALSE.
      RETURN
      END
C
C*************************************************************************
C
C  Check the Hessian index - assumes no Hessian is available.
C
      SUBROUTINE CHECKIND2(COORDS,MFLAG,INEG,ENERGY)
      USE COMMONS
      USE KEY
      USE VECCK
      USE ZWK
      IMPLICIT NONE

      INTEGER J1, J2, INEG, NS, ISEED
      DOUBLE PRECISION COORDS(3*NATOMS),ENERGY,DUMMY1, VEC2(3*NATOMS),EVALMIN, VECS(3*NATOMS),DPRAND, SOVER
      LOGICAL MFLAG, CONVERGED
      COMMON /IS/ ISEED

      DO J1=1,NOPT
         VECS(J1)=ZWORK(J1,1)
      ENDDO
      FIXIMAGE=.TRUE.

      IF (INEG.EQ.1) THEN
         DO J1=1,NOPT
            VECCHK(J1,1)=VECS(J1)
         ENDDO
      ENDIF

80    NMDONE=INEG
      DO J1=1,NOPT
         VEC2(J1)=2*(DPRAND()-0.5D0)
      ENDDO
      CALL VECNORM(VEC2,NOPT)

      CALL BEIG(1,COORDS,ENERGY,VEC2,EVALMIN,NS,SOVER,.TRUE.,CONVERGED)

      IF (EVALMIN.GT.0.0D0) THEN
         PRINT*,'Hessian index=',INEG
      ELSE
         INEG=INEG+1
         IF (INEG.GT.10) THEN
            PRINT*,' Hessian index is 10 or greater'
            FIXIMAGE=.FALSE.
            RETURN
         ENDIF
         DO J2=1,NOPT
            VECCHK(J2,INEG)=VEC2(J2)
         ENDDO
         GOTO 80
      ENDIF

      IF (((BFGSMINT).AND.(INEG.NE.0)).AND.CHECKCONT) THEN
         MFLAG=.FALSE.
         IF (IVEC.GE.0) PRINT*,'Stepping away from saddle along mode ',1,' + direction'
         IF (IVEC.LT.0) PRINT*,'Stepping away from saddle along mode ',1,' - direction'
         DUMMY1=PUSHOFF
         IF (IVEC.LT.0) DUMMY1=-DUMMY1
         DO J1=1,NOPT
            COORDS(J1)=COORDS(J1)+DUMMY1*VECCHK(J1,1)
         ENDDO
         INEG=0
         NMDONE=0
      ENDIF
      IF (((BFGSTST).AND.(INEG.NE.1)).AND.CHECKCONT) THEN
         MFLAG=.FALSE.
         IF (IVEC.GE.0) PRINT*,'Stepping away from saddle along mode ',1,' + direction'
         IF (IVEC.LT.0) PRINT*,'Stepping away from saddle along mode ',1,' - direction'
         DUMMY1=PUSHOFF
         IF (IVEC.LT.0) DUMMY1=-DUMMY1
         DO J1=1,NOPT
            COORDS(J1)=COORDS(J1)+DUMMY1*VECCHK(J1,2)
         ENDDO
         INEG=0
         NMDONE=0
      ENDIF

      FIXIMAGE=.FALSE.
      RETURN
      END

C
C  Aitken's extrapolation method. See "Numerical Methods that Work"  p. 216.
C
      SUBROUTINE EXTRAP(NSTEP,NOPT,VEC,PERCENT)
      IMPLICIT NONE
      INTEGER I,NOPT,NSTEP
      DOUBLE PRECISION VECT(NOPT,3),VEC(NOPT),Z,B,D,DUMMY,PERCENT

      DO I=1,NOPT
         VECT(I,1)=VECT(I,2)
         VECT(I,2)=VECT(I,3)
         VECT(I,3)=VEC(I)   !  is this right? 
      ENDDO
      IF ((NSTEP.LT.3).OR.(MOD(NSTEP,3).NE.0)) RETURN
      IF (PERCENT.GT.10.0D0) RETURN

      DO I=1,NOPT
         D=VECT(I,1)-2.0*VECT(I,2)+VECT(I,3)
         B=VECT(I,2)-VECT(I,3)
         IF (D.EQ.0.0D0) THEN
            Z=0.0D0
         ELSE
            Z=B*B/D
         ENDIF
         IF (ABS(Z).GT.ABS(VEC(I))/20.0D0) Z=0.0D0
         VECT(I,1)=Z
      ENDDO
      PRINT*,'Extrapolating'
      DUMMY=0.0D0
      DO I=1,NOPT
         VEC(I)=VECT(I,3)-VECT(I,1)
         DUMMY=DUMMY+VEC(I)**2
      ENDDO
      DUMMY=1.0D0/SQRT(DUMMY)
      PRINT*,'DUMMY=',DUMMY
      DO I=1,NOPT
         VEC(I)=VEC(I)*DUMMY
      ENDDO

      RETURN
      END

      DOUBLE PRECISION FUNCTION MINIM(a,b,boxl)
      IMPLICIT NONE
      DOUBLE PRECISION A, B, BOXL

      MINIM=A-B-BOXL*ANINT((A-B)/BOXL)

      RETURN
      END

      SUBROUTINE HSMOVE(COORDS,COORDSN,FIXDIR,STEP,PTEST)
      USE COMMONS
      USE KEY
      use porfuncs
      IMPLICIT NONE
      INTEGER J1, J2, PARTNR(NATOMS), J3, NTYPEA, NMOVE, NUMBERS(NATOMS), NPAIR, ISTAT, NCOLL
      DOUBLE PRECISION COORDS(3*NATOMS), FIXDIR(3*NATOMS), COLTIM(NATOMS), TIMBIG, RX12, RY12, RZ12, MINIM, STEP,
     1                 VX12, VY12, VZ12, B12, R12SQ, V12SQ, DISCR, SIGSQ, T12, COLTIM2(NATOMS), T122, DNEW, COORDSN(3*NATOMS)
      DOUBLE PRECISION EPSAB, EPSBB, SIGAB, SIGBB
      LOGICAL PTEST
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /HS/ NMOVE

      TIMBIG=1.0D100
      IF (.NOT.BINARY) SIGSQ=1.0D0
      IF (ZSYM(1)(1:2).EQ.'C6') SIGSQ=7.0D0
      IF (ZSYM(1)(1:2).EQ.'CA') SIGSQ=45.0D0

      DO J1=1,NATOMS
         COLTIM(J1)=TIMBIG
      ENDDO
      NCOLL=0

      DO J1=1,3*NATOMS
         COORDSN(J1)=COORDS(J1)
      ENDDO
     
C     DO J1=1,NATOMS
C        J2=3*(J1-1)
C        WRITE(*,'(A4,6G17.10)') 'AX  ',COORDS(J2+1),COORDS(J2+2),COORDS(J2+3),FIXDIR(J2+1),FIXDIR(J2+2),FIXDIR(J2+3)
C     ENDDO

      IF (TWOD) THEN
         DO J1=1,NATOMS
            FIXDIR(3*(J1-1)+3)=0.0D0
         ENDDO
      ENDIF

      DO J1=1,NATOMS
         DO J2=J1+1,NATOMS

            IF (BULKT) THEN
               RX12=MINIM(COORDSN(3*(J1-1)+1),COORDSN(3*(J2-1)+1),PARAM1)
               RY12=MINIM(COORDSN(3*(J1-1)+2),COORDSN(3*(J2-1)+2),PARAM2)
               RZ12=MINIM(COORDSN(3*(J1-1)+3),COORDSN(3*(J2-1)+3),PARAM3)
            ELSE
               RX12=COORDSN(3*(J1-1)+1)-COORDSN(3*(J2-1)+1)
               RY12=COORDSN(3*(J1-1)+2)-COORDSN(3*(J2-1)+2)
               RZ12=COORDSN(3*(J1-1)+3)-COORDSN(3*(J2-1)+3)
            ENDIF

            VX12=FIXDIR(3*(J1-1)+1)-FIXDIR(3*(J2-1)+1)
            VY12=FIXDIR(3*(J1-1)+2)-FIXDIR(3*(J2-1)+2)
            VZ12=FIXDIR(3*(J1-1)+3)-FIXDIR(3*(J2-1)+3)

            B12=RX12*VX12+RY12*VY12+RZ12*VZ12
C           PRINT*,'J1,J2,B12=',J1,J2,B12

            IF (B12.LT.0.0D0) THEN
               R12SQ=RX12**2+RY12**2+RZ12**2
               V12SQ=VX12**2+VY12**2+VZ12**2
               IF (BINARY) THEN
                  IF (J1.LE.NTYPEA) THEN
                     IF (J2.LE.NTYPEA) THEN
                        DISCR=B12**2-V12SQ*(R12SQ-1.0D0)
                     ELSE
                        DISCR=B12**2-V12SQ*(R12SQ-SIGAB**2*0.81D0)
                     ENDIF
                  ELSE
                     DISCR=B12**2-V12SQ*(R12SQ-SIGBB**2*0.81D0)
                  ENDIF
               ELSE
                  DISCR=B12**2-V12SQ*(R12SQ-SIGSQ)
               ENDIF
C              WRITE(*,'(A,5G20.10)') 'B12,V12SQ,R12SQ,SIGSQ,DISCR=',B12,V12SQ,R12SQ,SIGSQ,DISCR
               IF (DISCR.GT.0.0D0) THEN
                  T12=(-B12-SQRT(DISCR))/V12SQ
C                 WRITE(*,'(A,2I4,4G20.10)') 'J1,J2,T12,B12,SQRT(DISCR),V12SQ=',J1,J2,T12,B12,SQRT(DISCR),V12SQ
   
                  IF (ABS(T12).LT.COLTIM(J1)) THEN
                     NCOLL=NCOLL+1
                     COLTIM(J1)=ABS(T12)
                     COLTIM2(J1)=(-B12+SQRT(DISCR))/V12SQ
                     PARTNR(J1)=J2
                  ENDIF
   
                  IF (ABS(T12).LT.COLTIM(J2)) THEN
                     COLTIM(J2)=ABS(T12)
                     COLTIM2(J2)=(-B12+SQRT(DISCR))/V12SQ
                     PARTNR(J2)=J1
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      IF (NCOLL.LT.NMOVE) THEN
         WRITE(*,'(A)') ' WARNING - insufficient  collisions for input velocities'
         COORDSN(1:3*NATOMS)=COORDSN(1:3*NATOMS)+FIXDIR(1:3*NATOMS)/2.0D0
         FIXDIR(1:3*NATOMS)=COORDSN(1:3*NATOMS)-COORDS(1:3*NATOMS)
         STEP=0.0D0
         DO J3=1,3*NATOMS
            STEP=STEP+(COORDSN(J3)-COORDS(J3))**2
         ENDDO
         STEP=SQRT(STEP)
         CALL ORTHOGOPT(FIXDIR,COORDS,.TRUE.)
         RETURN
      ENDIF

      DO J1=1,NATOMS
         NUMBERS(J1)=J1
      ENDDO
      CALL SORTN(NATOMS,COLTIM,NUMBERS)

C     T12=TIMBIG
C     DO J1=1,NATOMS
C        IF (COLTIM(J1).LT.T12) THEN
C           T12=COLTIM(J1)
C           T122=COLTIM2(J1)
C           J2=J1
C        ENDIF
C     ENDDO

      DO NPAIR=1,NMOVE

      J2=NUMBERS(NPAIR)
      J1=PARTNR(J2)
C     PRINT*,'NPAIR,J2,J1=',NPAIR,J2,J1
      T12=COLTIM(NPAIR)
      T122=COLTIM2(J2)

      IF (PTEST) WRITE(*,'(A,I4,A,I4,A,I4,A,G15.6)') ' collision ',NPAIR,' is between atoms ',J1,' and ',J2,' at time ',T12
C
C  Advance all positions to time T12*T12FAC. 
C
      IF (T12FAC.LE.1.0D0) THEN
         DO J3=1,NATOMS
            COLTIM(J3)=COLTIM(J3)-T12*T12FAC
            COORDSN(3*(J3-1)+1)=COORDSN(3*(J3-1)+1)+FIXDIR(3*(J3-1)+1)*T12*T12FAC
            COORDSN(3*(J3-1)+2)=COORDSN(3*(J3-1)+2)+FIXDIR(3*(J3-1)+2)*T12*T12FAC
            COORDSN(3*(J3-1)+3)=COORDSN(3*(J3-1)+3)+FIXDIR(3*(J3-1)+3)*T12*T12FAC
C
C  Could put atoms leaving the primary supercell back in the box here
C  if desired.
C
         ENDDO
C
C  Advance positions of the colliding pair only to time T12. 
C
C        COORDSN(3*(J1-1)+1)=COORDSN(3*(J1-1)+1)+FIXDIR(3*(J1-1)+1)*T12
C        COORDSN(3*(J1-1)+2)=COORDSN(3*(J1-1)+2)+FIXDIR(3*(J1-1)+2)*T12
C        COORDSN(3*(J1-1)+3)=COORDSN(3*(J1-1)+3)+FIXDIR(3*(J1-1)+3)*T12
C        COORDSN(3*(J2-1)+1)=COORDSN(3*(J2-1)+1)+FIXDIR(3*(J2-1)+1)*T12
C        COORDSN(3*(J2-1)+2)=COORDSN(3*(J2-1)+2)+FIXDIR(3*(J2-1)+2)*T12
C        COORDSN(3*(J2-1)+3)=COORDSN(3*(J2-1)+3)+FIXDIR(3*(J2-1)+3)*T12

         IF (BULKT) THEN
            IF (PTEST) PRINT*,'separation between these atoms is now ',
     1                                  SQRT(MINIM(COORDSN(3*(J1-1)+1),COORDSN(3*(J2-1)+1),PARAM1)**2
     2                                      +MINIM(COORDSN(3*(J1-1)+2),COORDSN(3*(J2-1)+2),PARAM2)**2
     3                                      +MINIM(COORDSN(3*(J1-1)+3),COORDSN(3*(J2-1)+3),PARAM3)**2)
         ELSE
            IF (PTEST) PRINT*,'separation between these atoms is now ',SQRT((COORDSN(3*(J1-1)+1)-COORDSN(3*(J2-1)+1))**2
     1                                      +(COORDSN(3*(J1-1)+2)-COORDSN(3*(J2-1)+2))**2
     2                                      +(COORDSN(3*(J1-1)+3)-COORDSN(3*(J2-1)+3))**2)
         ENDIF
      ELSE 
C
C  Put the first colliding pair half way between their entrance/exit positions
C

         COORDSN(3*(J1-1)+1)=COORDSN(3*(J1-1)+1)+FIXDIR(3*(J1-1)+1)*(T122-T12)/2
         COORDSN(3*(J1-1)+2)=COORDSN(3*(J1-1)+2)+FIXDIR(3*(J1-1)+2)*(T122-T12)/2
         COORDSN(3*(J1-1)+3)=COORDSN(3*(J1-1)+3)+FIXDIR(3*(J1-1)+3)*(T122-T12)/2
         COORDSN(3*(J2-1)+1)=COORDSN(3*(J2-1)+1)+FIXDIR(3*(J2-1)+1)*(T122-T12)/2
         COORDSN(3*(J2-1)+2)=COORDSN(3*(J2-1)+2)+FIXDIR(3*(J2-1)+2)*(T122-T12)/2
         COORDSN(3*(J2-1)+3)=COORDSN(3*(J2-1)+3)+FIXDIR(3*(J2-1)+3)*(T122-T12)/2

         IF (BULKT) THEN
            DNEW=SQRT(MINIM(COORDSN(3*(J1-1)+1),COORDSN(3*(J2-1)+1),PARAM1)**2
     1            +MINIM(COORDSN(3*(J1-1)+2),COORDSN(3*(J2-1)+2),PARAM2)**2
     2            +MINIM(COORDSN(3*(J1-1)+3),COORDSN(3*(J2-1)+3),PARAM3)**2)
         ELSE
            DNEW=SQRT((COORDSN(3*(J1-1)+1)-COORDSN(3*(J2-1)+1))**2
     1            +(COORDSN(3*(J1-1)+2)-COORDSN(3*(J2-1)+2))**2
     2            +(COORDSN(3*(J1-1)+3)-COORDSN(3*(J2-1)+3))**2)
         ENDIF
         IF (PTEST) PRINT*,'separation half way between entrance and exit=',DNEW

         IF (BINARY) THEN
            IF (J1.LE.NTYPEA) THEN
               IF (J2.LE.NTYPEA) THEN
                  SIGSQ=1.0D0
               ELSE
                  SIGSQ=SIGAB**2*0.81D0
               ENDIF
            ELSE
               SIGSQ=SIGBB**2*0.81D0
            ENDIF
         ENDIF
C
C  Rescale the distance between the first colliding pair
C
         IF (BULKT) THEN
            RX12=MINIM(COORDSN(3*(J1-1)+1),COORDSN(3*(J2-1)+1),PARAM1)*SQRT(SIGSQ)/DNEW/2
            RY12=MINIM(COORDSN(3*(J1-1)+2),COORDSN(3*(J2-1)+2),PARAM2)*SQRT(SIGSQ)/DNEW/2
            RZ12=MINIM(COORDSN(3*(J1-1)+3),COORDSN(3*(J2-1)+3),PARAM3)*SQRT(SIGSQ)/DNEW/2

            VX12=MINIM(COORDSN(3*(J1-1)+1),COORDSN(3*(J2-1)+1),PARAM1)
            VY12=MINIM(COORDSN(3*(J1-1)+2),COORDSN(3*(J2-1)+2),PARAM2)
            VZ12=MINIM(COORDSN(3*(J1-1)+3),COORDSN(3*(J2-1)+3),PARAM3)

            COORDSN(3*(J1-1)+1)=COORDSN(3*(J1-1)+1)-VX12/2+RX12
            COORDSN(3*(J1-1)+2)=COORDSN(3*(J1-1)+2)-VY12/2+RY12
            COORDSN(3*(J1-1)+3)=COORDSN(3*(J1-1)+3)-VZ12/2+RZ12
            COORDSN(3*(J2-1)+1)=COORDSN(3*(J2-1)+1)+VX12/2-RX12
            COORDSN(3*(J2-1)+2)=COORDSN(3*(J2-1)+2)+VY12/2-RY12
            COORDSN(3*(J2-1)+3)=COORDSN(3*(J2-1)+3)+VZ12/2-RZ12

            DNEW=SQRT(MINIM(COORDSN(3*(J1-1)+1),COORDSN(3*(J2-1)+1),PARAM1)**2
     1            +MINIM(COORDSN(3*(J1-1)+2),COORDSN(3*(J2-1)+2),PARAM2)**2
     2            +MINIM(COORDSN(3*(J1-1)+3),COORDSN(3*(J2-1)+3),PARAM3)**2)
         ELSE
            RX12=(COORDSN(3*(J1-1)+1)-COORDSN(3*(J2-1)+1))*SQRT(SIGSQ)/DNEW/2
            RY12=(COORDSN(3*(J1-1)+2)-COORDSN(3*(J2-1)+2))*SQRT(SIGSQ)/DNEW/2
            RZ12=(COORDSN(3*(J1-1)+3)-COORDSN(3*(J2-1)+3))*SQRT(SIGSQ)/DNEW/2
            VX12=(COORDSN(3*(J1-1)+1)+COORDSN(3*(J2-1)+1))/2
            VY12=(COORDSN(3*(J1-1)+2)+COORDSN(3*(J2-1)+2))/2
            VZ12=(COORDSN(3*(J1-1)+3)+COORDSN(3*(J2-1)+3))/2
            COORDSN(3*(J1-1)+1)=VX12+RX12
            COORDSN(3*(J1-1)+2)=VY12+RY12
            COORDSN(3*(J1-1)+3)=VZ12+RZ12

            COORDSN(3*(J2-1)+1)=VX12-RX12
            COORDSN(3*(J2-1)+2)=VY12-RY12
            COORDSN(3*(J2-1)+3)=VZ12-RZ12

            DNEW=SQRT((COORDSN(3*(J1-1)+1)-COORDSN(3*(J2-1)+1))**2
     1               +(COORDSN(3*(J1-1)+2)-COORDSN(3*(J2-1)+2))**2
     2               +(COORDSN(3*(J1-1)+3)-COORDSN(3*(J2-1)+3))**2)
         ENDIF
         IF (PTEST) PRINT*,'separation reset to ',DNEW
      ENDIF

      STEP=0.0D0
      DO J3=1,3*NATOMS
         FIXDIR(J3)=COORDSN(J3)-COORDS(J3)
         STEP=STEP+(COORDSN(J3)-COORDS(J3))**2
      ENDDO
      STEP=SQRT(STEP)
      CALL ORTHOGOPT(FIXDIR,COORDS,.TRUE.)

      ENDDO

      CALL FLUSH(6,ISTAT)

      RETURN
      END

      SUBROUTINE SORTN(N,A,NA)
      IMPLICIT NONE
      INTEGER J1, L, N, J2, NA(N), NTEMP
      DOUBLE PRECISION TEMP, A(N)
C
      DO J1=1,N-1
         L=J1
         DO J2=J1+1,N
            IF (A(L).GT.A(J2)) L=J2
         ENDDO
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         NTEMP=NA(L)
         NA(L)=NA(J1)
         NA(J1)=NTEMP
      ENDDO
      RETURN
      END
