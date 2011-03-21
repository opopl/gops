C
C GPL License Info {{{
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
C }}}
C
C jmc Note nothing has been done here to fix the unres pathlength coordinate resetting problem...
C
      SUBROUTINE PATH(Q,ENERGY,VNEW,RMS,EVTS,VECS,POTCALL,QPLUS,QMINUS,PTEST,ETS,EPLUS,EMINUS,SLENGTH,DISP,GAMMA,NTILDE,
     1                FRQSTS,FRQSPLUS,FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)
! Declarations {{{
      ! Module declarations {{{
      USE COMMONS
      USE KEY
      USE SYMINF
      USE modcharmm
      USE MODUNRES
      USE MODHESS
      USE MODEFOL
      USE PORFUNCS
      USE AMHGLOBALS, ONLY : NMRES
      ! }}}

      IMPLICIT NONE

! subroutine parameters {{{

      DOUBLE PRECISION ENERGY, VNEW(3*NATOMS), RMS, EVTS, VECS(3*NATOMS)
      DOUBLE PRECISION Q(3*NATOMS), QPLUS(3*NATOMS), QMINUS(3*NATOMS)
      DOUBLE PRECISION ETS, EPLUS, EMINUS, SLENGTH, DISP, GAMMA, NTILDE
      DOUBLE PRECISION FRQSTS(3*NATOMS), FRQSPLUS(3*NATOMS), FRQSMINUS(3*NATOMS)
       
      LOGICAL POTCALL, PTEST
      LOGICAL PATHFAILT
      CHARACTER(LEN=*) ITSTRING,EOFSSTRING

! }}}

! local parameters {{{

      INTEGER NSTEPPLUS, ITDONE, NSTEPMINUS, J1, J2, NPATHFRAME, NATOMSSAVE, NEWINR, 
     1        INEG, HORDER, INFO, IVECSAVE, IVEC2SAVE, IPOT, J3, NFPLUS, NFMINUS, RECLEN, ISTAT, NUSE
      INTEGER, PARAMETER :: NSTEPMAX=50000, NFMAX=20000

      DOUBLE PRECISION EVALMAX, RANDOM,
     1                 DIAG(3*NATOMS), EREAL, RMS2, STEP(3*NATOMS), QINIT(3*NATOMS),
     2                 EOFS(NSTEPMAX), PATHLENGTH(NSTEPMAX), TEMP, MINIM, EOFSFRAMEP(NSTEPMAX),
     3                 SUM2, SUM4, EVPLUS, EVMINUS, EOFSFRAMEM(NSTEPMAX),
     4                 SPLUS, SMINUS, STS, STEMP, ETEMP, DUMMY, TIME, TIME0,
     5                 OVEC(3), H1VEC(3), H2VEC(3),
     6                 TEMPA(9*NATOMS), CAPSCOORDS1(18), CAPSCOORDS2(18), 
     7                 DPRAND, PPLUS, PMINUS, lambdats, lambdap, lambdam, distp, distm, RMAT(3,3) !, P(3)
      LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST, ETEST, NOSHIFTSAVE, NOHESSSAVE
      DOUBLE PRECISION TEMPERATURE, HRED, DIHE, ALLANG, LASTE
      INTEGER NCONNECT
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
      DOUBLE PRECISION CAPSRHO, EPS2, RAD, HEIGHT
      COMMON /CAPS/ CAPSRHO, EPS2, RAD, HEIGHT

C jmc
      CHARACTER(LEN=80) ITSTRING2
      DOUBLE PRECISION NEWINT(NINTS),TSINT(NINTS)
      DOUBLE PRECISION PEPCOORDS(3*NATOMS), INERTIA(3,3)
      INTEGER K1,K2,KD,NNZ,NINTB
C jmc

C    LOCAL AMH VARIABLES
      INTEGER GLY_COUNT
C      CHARACTER(LEN=5) TARFL

      LOGICAL MFLAG, BFGSTSTSAVE
      LOGICAL PATHT, DRAGT
      COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      INTEGER NATOMSIMUL
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Q1, Q2, QW
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: QFRAMEP, QFRAMEM
! }}}
! }}}
! Subroutine body {{{ 
      PATHFAILT=.FALSE.

      IF (BULKT.AND.(PARAM1*PARAM2*PARAM3.EQ.0.0D0)) THEN
        ! account for zero box size {{{
         IF (TWOD) THEN
            IF (PARAM1*PARAM2.EQ.0.0D0) THEN
               PRINT '(A)',' path> ERROR - BULKT is true but a box parameter is zero' 
               STOP
            ENDIF
         ELSE 
            PRINT '(A)',' path> ERROR - BULKT is true but a box parameter is zero' 
            STOP
         ENDIF
         ! }}}
      ENDIF

C  Calls to dumpp dump energies and coordinates on units 1 and 2. Must rewind for path to work.
C  If PRINTPTS is .FALSE. then reduce the I/O to a bare minimum. Use
C  QPLUS, QINIT and QMINUS for the stationary points rather than reading through
C  points saved on disk.
C
      IF (PRINTPTS) THEN 
         REWIND(1)  ! points file
         REWIND(2)  ! energies file
      ENDIF
      BFGSTSTSAVE=BFGSTST
      IVECSAVE=IVEC
      IVEC2SAVE=IVEC2
C
C  Plus side first. {{{
C
      IVEC=1
      DO J1=1,NOPT
         QINIT(J1)=Q(J1)
      ENDDO
!     IF (DEBUG) PRINT*,'ts points in path:'
!     IF (DEBUG) WRITE(*,'(3G20.10)') (Q(J1),J1=1,NOPT)

      CALL MYCPU_TIME(TIME0,.FALSE.)
      IF (HYBRIDMINT) THEN
        ! use hybrid minimization, keywords: HYBRIDMIN {{{
!
! We must have the ts energy somewhere already. Should;t really need the call to potential.
!
         CALL POTENTIAL(Q,ETS,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         ENERGY=ETS
         POTCALL=.TRUE.
         CALL DUMPP(Q,ENERGY)
         CALL HYBRIDMIN(HMNSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVTS,EVALMAX,VECS,ITDONE,POTCALL,PTEST)
         NSTEPPLUS=ITDONE
         IF (.NOT.MFLAG) THEN
            IF (PTEST) PRINT '(A,I8,A)','efol> switching to LBFGS minimisation after ',NSTEPPLUS,' hybrid minimisation steps'
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            END IF
            NSTEPPLUS=NSTEPPLUS+ITDONE
         ENDIF
         ! }}}
      ELSE IF (BFGSTST) THEN
        ! keywords: BFGSTS BFGSSTEP 
        !       BFGSTS: hybrid BFGS/eigenvector-following transition state search {{{
         BFGSSTEP=.TRUE.
         IF (UNRST) THEN
            CALL INTBFGSTS(NSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVTS,EVALMAX,VECS,ITDONE,POTCALL,PTEST)
         ELSE
            CALL BFGSTS(NSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVTS,EVALMAX,VECS,ITDONE,POTCALL,PTEST)
         ENDIF
         ETS=ENERGY
!        IF (DEBUG) PRINT*,'ts step off plus points in path:'
!        IF (DEBUG) WRITE(*,'(3G20.10)') (Q(J1),J1=1,NOPT)
         DO J1=1,NOPT
            STEP(J1)=Q(J1)-QINIT(J1)
         ENDDO
         
         IF (RKMIN) RMS=1.0D0
         MFLAG=.FALSE.
         BFGSSTEP=.FALSE.
         BFGSTST=.FALSE.
         IF (BFGSMINT) THEN
            IF (UNRST.OR.(CHRMMT.AND.INTMINT)) THEN
                CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                       .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
                CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                       .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEPPLUS=ITDONE+1
         ELSE IF (BBRSDMT) THEN
            CALL BBRSDM(Q,MFLAG,ITDONE,ENERGY,RMS,.FALSE.,VNEW,PTEST)
            NSTEPPLUS=ITDONE+1
C
C DAE to switch to BFGS after NSTEPS sd
C
            IF (.NOT.MFLAG) THEN
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ENDIF
         ELSE IF (BSMIN.OR.RKMIN) THEN
            NUSE=NSTEPS
            IF (PATHSDSTEPS.GT.0) NUSE=PATHSDSTEPS
            CALL ODESD(NUSE,Q,MFLAG,ITDONE,PTEST)
            NSTEPPLUS=ITDONE+1
C
C DAE to switch to BFGS after NSTEPS sd
C
            IF (.NOT.MFLAG) THEN
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ELSE
!
!  ODESD does not return the energy!
!
               CALL POTENTIAL(Q,ENERGY,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ENDIF
C end DAE
         ELSE
            IF (INR.LT.6) THEN
               NEWINR=0 ! so we can use Page-McIver
            ELSE
               NEWINR=INR
            ENDIF
C bs360 (29/07/08): ACE does not converge with SD. I will investigate it later.
            IF (ACESOLV) THEN
               NUSE=20
            ELSE
               NUSE=NSTEPS
            ENDIF
            IF ((PATHSDSTEPS.GT.0).AND.(INR.LE.6)) NUSE=PATHSDSTEPS
            NOSHIFTSAVE=NOSHIFT; NOHESSSAVE=NOHESS
            NOSHIFT=.FALSE.; NOHESS=.FALSE.
            CALL EFOL(Q,MFLAG,NUSE,ENERGY,ITDONE,EVPLUS,PTEST,FRQSPLUS,NEWINR)
            NOSHIFT=NOSHIFTSAVE; NOHESS=NOHESSSAVE
            NSTEPPLUS=ITDONE
!
!  Switch to LBFGS if SD did not finish.
!
            IF (.NOT.MFLAG) THEN
               IF (PTEST) PRINT '(A,I8,A)','efol> switching to LBFGS minimisation after ',NUSE,' steepest-descent steps'
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ENDIF
         ENDIF
         ! }}}
      ELSE
         ! {{{
         CALL EFOL(Q,MFLAG,1,ENERGY,ITDONE,EVTS,PTEST,FRQSTS,0)
         ETS=ENERGY
         DO J1=1,NOPT
            STEP(J1)=Q(J1)-QINIT(J1)
         ENDDO
         KNOWE=.FALSE.
         KNOWG=.FALSE. ! we don`t know the gradient at the point in Q because VNEW isn`t passed from efol?
         IF (BFGSMINT) THEN
C           NOSHIFT=.FALSE.
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            END IF

            NSTEPPLUS=ITDONE+1
         ELSE IF (BBRSDMT) THEN
            CALL BBRSDM(Q,MFLAG,ITDONE,ENERGY,RMS,.FALSE.,VNEW,PTEST)
            NSTEPPLUS=ITDONE+1
C
C DAE to switch to BFGS after NSTEPS sd
C
            IF (.NOT.MFLAG) THEN
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ENDIF
         ELSE IF (BSMIN.OR.RKMIN) THEN
C           NOSHIFT=.FALSE.
            NUSE=NSTEPS
            IF (PATHSDSTEPS.GT.0) NUSE=PATHSDSTEPS
            CALL ODESD(NUSE,Q,MFLAG,ITDONE,PTEST)
            NSTEPPLUS=ITDONE+1
C DAE to switch to BFGS after NSTEPS sd
            IF (.NOT.MFLAG) THEN
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ELSE
!
!  ODESD does not return the energy!
!
               CALL POTENTIAL(Q,ENERGY,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ENDIF
C end DAE
         ELSE
            IF (INR.LT.6) THEN
               NEWINR=0 ! so we can use Page-McIver
            ELSE
               NEWINR=INR
            ENDIF
            NUSE=NSTEPS
            IF ((PATHSDSTEPS.GT.0).AND.(INR.LE.6)) NUSE=PATHSDSTEPS
            NOSHIFTSAVE=NOSHIFT; NOHESSSAVE=NOHESS
            NOSHIFT=.FALSE.; NOHESS=.FALSE.
            CALL EFOL(Q,MFLAG,NUSE,ENERGY,ITDONE,EVPLUS,PTEST,FRQSPLUS,NEWINR)
            NOSHIFT=NOSHIFTSAVE; NOHESS=NOHESSSAVE
            NSTEPPLUS=ITDONE ! bug fix DJW 7/10/08
!
!  Switch to LBFGS if SD did not finish.
!
            IF (.NOT.MFLAG) THEN
               IF (PTEST) PRINT '(A,I8,A)','efol> switching to LBFGS minimisation after ',NUSE,' steepest-descent steps'
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ENDIF
         ENDIF
         ! }}}
      ENDIF
      CALL MYCPU_TIME(TIME,.FALSE.)
      IF (MFLAG) THEN
         ! {{{
         WRITE(*,'(A,I20,A,G20.10,2X,A,F11.2)') ' Plus  side of path:    ',NSTEPPLUS,' steps. Energy=',ENERGY,' time=',TIME-TIME0
         CALL FLUSH(6,ISTAT)
         ! }}}
      ELSE
         ! {{{
         WRITE(*,'(A,I20,A)') ' Plus  side of path failed to converge in ',NSTEPPLUS,' steps'
c        STOP
         PATHFAILT=.TRUE.
         BFGSTST=BFGSTSTSAVE
         IVEC=IVECSAVE
         IVEC2=IVEC2SAVE
C
C jmc Note that before the plus-side minimization above, BFGSTST is set to false, so if we`re doing a connect run (and BFGSTS
C was initially true), the next ts search will mess up, calling efol not intbfgsts.  Reset to BFGSTSTSAVE and also reset IVEC
C and IVEC2 here (as at the end of this subroutine).
C
         IF (ALLOCATED(Q1)) DEALLOCATE(Q1)
         IF (ALLOCATED(Q2)) DEALLOCATE(Q2)
         IF (ALLOCATED(QW)) DEALLOCATE(QW)
         IF (ALLOCATED(QFRAMEP)) DEALLOCATE(QFRAMEP)
         IF (ALLOCATED(QFRAMEM)) DEALLOCATE(QFRAMEM)
         ! }}}
         RETURN
      ENDIF

C Check Hessian index  

      IF ((MFLAG.AND.CHECKINDEX).AND.(BFGSMINT.OR.BFGSTST.OR.BSMIN.OR.RKMIN)) THEN
C  {{{
         IF (NOHESS) THEN
            CALL CHECKIND2(Q,MFLAG,INEG,ENERGY)
C }}}
         ELSE
C  We need the Hessian in CHECKIND. {{{
C
            IF (BFGSMINT.OR.BSMIN.OR.RKMIN) CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
            CALL CHECKIND(Q,MFLAG,INEG,ENERGY,EVPLUS,EVALMAX,.FALSE.)
         ENDIF
C }}}
      ENDIF
C }}}
C  Minus side. {{{
C
      IVEC=-1
      IF (REDOPATH1) THEN
         REDOPATH1=.FALSE.
         REDOPATH2=.TRUE.
      ELSEIF (REDOPATH2) THEN
         REDOPATH1=.TRUE.
         REDOPATH2=.FALSE.
      ENDIF
      DO J1=1,NOPT
         QPLUS(J1)=Q(J1)
         IF (.NOT.UNRST) Q(J1)=QINIT(J1)-STEP(J1)
         IF (HYBRIDMINT) Q(J1)=QINIT(J1)
      ENDDO
      EPLUS=ENERGY
! }}}
!     IF (DEBUG) PRINT*,'ts step off minus points in path:'
!     IF (DEBUG) WRITE(*,'(3G20.10)') (Q(J1),J1=1,NOPT)

      IF (UNRST) THEN ! jmc new intstep stuff 
         ! {{{
         DO J1=1,nres
            c(1,J1)=QINIT(6*(J1-1)+1)
            c(2,J1)=QINIT(6*(J1-1)+2)
            c(3,J1)=QINIT(6*(J1-1)+3)
            c(1,J1+nres)=QINIT(6*(J1-1)+4)
            c(2,J1+nres)=QINIT(6*(J1-1)+5)
            c(3,J1+nres)=QINIT(6*(J1-1)+6)
         END DO
         CALL UPDATEDC
         CALL int_from_cart(.true.,.false.)
         CALL geom_to_var(NINTS,TSINT)
         NEWINT=TSINT-INTSTEP
         CALL var_to_geom(NINTS,NEWINT)
         CALL chainbuild
         DO J1=1,nres
            Q(6*(J1-1)+1)=c(1,J1)
            Q(6*(J1-1)+2)=c(2,J1)
            Q(6*(J1-1)+3)=c(3,J1)
            Q(6*(J1-1)+4)=c(1,J1+nres)
            Q(6*(J1-1)+5)=c(2,J1+nres)
            Q(6*(J1-1)+6)=c(3,J1+nres)
         END DO
         ! }}}
      ENDIF ! jmc end new stuff

      KNOWE=.FALSE.
      KNOWG=.FALSE.

      CALL MYCPU_TIME(TIME0,.FALSE.)
      
      ! HYBRIDMINT BFGSMIT BBRSDMT BSMIN RKMIN;  ELSE
! {{{
      IF (HYBRIDMINT) THEN
        ! {{{
         POTCALL=.TRUE.
         ENERGY=ETS
         CALL HYBRIDMIN(HMNSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVTS,EVALMAX,VECS,ITDONE,POTCALL,PTEST)
         NSTEPMINUS=ITDONE
         IF (.NOT.MFLAG) THEN
            IF (PTEST) PRINT '(A,I8,A)','efol> switching to LBFGS minimisation after ',NSTEPMINUS,' hybrid minimisation steps'
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEPMINUS=NSTEPMINUS+ITDONE
         ENDIF
         ! }}}
      ELSE IF (BFGSMINT) THEN
         ! {{{
         IF (UNRST.OR.(CHRMMT.AND.INTMINT)) THEN
            CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                   .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
         ELSE
            CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                   .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
         ENDIF
         NSTEPMINUS=ITDONE+1
         ! }}}
      ELSE IF (BBRSDMT) THEN
         ! {{{
         CALL BBRSDM(Q,MFLAG,ITDONE,ENERGY,RMS,.FALSE.,VNEW,PTEST)
         NSTEPMINUS=ITDONE+1
C
C DAE to switch to BFGS after NSTEPS sd
C
         IF (.NOT.MFLAG) THEN
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEPMINUS=NSTEPMINUS+ITDONE
         ENDIF
         ! }}}
      ELSE IF (BSMIN.OR.RKMIN) THEN
         ! {{{
         NUSE=NSTEPS
         IF (PATHSDSTEPS.GT.0) NUSE=PATHSDSTEPS
         CALL ODESD(NUSE,Q,MFLAG,ITDONE,PTEST)
         NSTEPMINUS=ITDONE+1
C DAE to switch to BFGS after NSTEPS sd
         IF (.NOT.MFLAG) THEN
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEPMINUS=NSTEPMINUS+ITDONE
         ELSE
!
!  ODESD does not return the energy!
!
            CALL POTENTIAL(Q,ENERGY,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         ENDIF
C end DAE
         ! }}}
      ELSE
        ! {{{
         IF (INR.LT.6) THEN
            NEWINR=0 ! so we can use Page-McIver
         ELSE
            NEWINR=INR
         ENDIF
C bs360 (29/07/08): ACE does not converge with SD. I will investigate it later.
         IF (ACESOLV) THEN
            NUSE=20
         ELSE
            NUSE=NSTEPS
         ENDIF
         IF ((PATHSDSTEPS.GT.0).AND.(INR.LE.6)) NUSE=PATHSDSTEPS
         NOSHIFTSAVE=NOSHIFT; NOHESSSAVE=NOHESS
         NOSHIFT=.FALSE.; NOHESS=.FALSE.
         CALL EFOL(Q,MFLAG,NUSE,ENERGY,ITDONE,EVMINUS,PTEST,FRQSMINUS,NEWINR)
         NOSHIFT=NOSHIFTSAVE; NOHESS=NOHESSSAVE
         NSTEPMINUS=ITDONE
!
!  Switch to LBFGS if SD did not finish.
!
         IF (.NOT.MFLAG) THEN
            IF (PTEST) PRINT '(A,I8,A)','efol> switching to LBFGS minimisation after ',NUSE,' steepest-descent steps'
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEPMINUS=NSTEPMINUS+ITDONE
         ENDIF
         ! }}}
      ENDIF
      ! }}}

      CALL MYCPU_TIME(TIME,.FALSE.)

      ! write info about minus side of energy {{{
      IF (MFLAG) THEN
        ! {{{
         WRITE(*,'(A,I20,A,G20.10,2X,A,F11.2)') ' Minus side of path:    ',NSTEPMINUS,' steps. Energy=',ENERGY,' time=',TIME-TIME0
         CALL FLUSH(6,ISTAT)
         ! }}}
      ELSE
        ! {{{
         WRITE(*,'(A,I20,A)') ' Minus side of path failed to converge in ',NSTEPMINUS,' steps'
c        STOP
         PATHFAILT=.TRUE.
         BFGSTST=BFGSTSTSAVE
         IVEC=IVECSAVE
         IVEC2=IVEC2SAVE
C jmc Note that before the plus-side minimization above, BFGSTST is set to false, so if we`re doing a connect run (and BFGSTS 
C was initially true), the next ts search will mess up, calling efol not intbfgsts.  Reset to BFGSTSTSAVE and also reset IVEC
C and IVEC2 here (as at the end of this subroutine).
         IF (ALLOCATED(Q1)) DEALLOCATE(Q1)
         IF (ALLOCATED(Q2)) DEALLOCATE(Q2)
         IF (ALLOCATED(QW)) DEALLOCATE(QW)
         IF (ALLOCATED(QFRAMEP)) DEALLOCATE(QFRAMEP)
         IF (ALLOCATED(QFRAMEM)) DEALLOCATE(QFRAMEM)
         ! }}}
         RETURN
      ENDIF

      ! }}}

      DO J1=1,NOPT
         QMINUS(J1)=Q(J1)
      ENDDO
      EMINUS=ENERGY
C
C  Check Hessian index {{{
C
      IF ((MFLAG.AND.CHECKINDEX).AND.(BFGSMINT.OR.BFGSTST.OR.BSMIN.OR.RKMIN)) THEN
         IF (NOHESS) THEN
            CALL CHECKIND2(Q,MFLAG,INEG,ENERGY)
         ELSE
C }}}
C  We need the Hessian in CHECKIND. {{{
C
            IF (BFGSMINT.OR.BSMIN.OR.RKMIN) CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
            CALL CHECKIND(Q,MFLAG,INEG,ENERGY,EVMINUS,EVALMAX,.FALSE.)
         ENDIF
      ENDIF
C }}}
C  The total number of energies and coordinates is NSTEPPLUS + NSTEPMINUS + 1 for the transition state.
C  The rest of this subroutine is post-processing {{{
C  Minimise IO if PRINTPTS is .FALSE.
C
      IF (ZSYMSAVE(1:1).EQ.'W') THEN  !  WCOMMENT
         ALLOCATE(QFRAMEP(9*(NATOMS/2),NFMAX),QFRAMEM(9*(NATOMS/2),NFMAX))
      ELSE
         ALLOCATE(QFRAMEP(3*NATOMS,NFMAX),QFRAMEM(3*NATOMS,NFMAX))
      ENDIF
      IF (.NOT.PRINTPTS) THEN 
         ! {{{
         NSTEPPLUS=1
         NSTEPMINUS=1
         NFPLUS=1
         NFMINUS=1
         EOFS(1)=EPLUS
         EOFS(2)=ETS
         EOFS(3)=EMINUS
         EOFSFRAMEP(1)=EPLUS
         EOFSFRAMEM(1)=EMINUS
         PATHLENGTH(1)=0.0D0
         DUMMY=0.0D0
C        PRINT*,'WARNING - S,N,gamma are not calculated from Cartesian coordinates here'
         DO J2=1,NOPT
            DUMMY=DUMMY+(QINIT(J2)-QPLUS(J2))**2
         ENDDO
         PATHLENGTH(2)=SQRT(DUMMY)
         DUMMY=0.0D0
         DO J2=1,NOPT
            DUMMY=DUMMY+(QINIT(J2)-QMINUS(J2))**2
         ENDDO
         PATHLENGTH(3)=SQRT(DUMMY)+PATHLENGTH(2)
         DO J2=1,NOPT
            QFRAMEP(J2,1)=QPLUS(J2)
            QFRAMEM(J2,1)=QMINUS(J2)
         ENDDO
         GOTO 555
         ! }}}
      ENDIF

      REWIND(1)
      REWIND(2)
      READ(2,*) EOFS(NSTEPPLUS+1)
      IF (ZSYMSAVE(1:1).EQ.'W') THEN  !  WCOMMENT
         ALLOCATE(Q1(9*(NATOMS/2)),Q2(9*(NATOMS/2)),QW(9*(NATOMS/2)))
C        READ(1,*) (Q1(J1),J1=1,9*NATOMS)
         READ(1,*) (Q1(J1),J1=1,9*(NATOMS/2))
      ELSE
         ALLOCATE(Q1(3*NATOMS),Q2(3*NATOMS))
         READ(1,*) (Q1(J1),J1=1,NOPT)
      ENDIF
      PATHLENGTH(NSTEPPLUS+1)=0.0D0
C     
C  The number of frames for each side of the path, NPATHFRAME, is now treated 
C  in an average way. We must dump the three stationary points at the very least,
C  and other frames on the two paths are then dumped with a probability proportional
C  to NPATHFRAME. If NPATHFRAME is > NFMAX then use NFMAX in its place, otherwise
C  we don;t have enough storage declared for the frames.
C
C  The new variable FRAMEDIST can be used to exclude frames with configurations that
C  are separated by less than FRAMEDIST distance units.
C
!     PPLUS=MIN(MIN(MIN(NPATHFRAME,NSTEPPLUS),NFMAX)*1.0D0/(1.0D0*(NSTEPPLUS-1)),1.0D0)
!     PMINUS=MIN(MIN(MIN(NPATHFRAME,NSTEPMINUS),NFMAX)*1.0D0/(1.0D0*(NSTEPMINUS-1)),1.0D0)

      if (NSTEPPLUS.lt.NPATHFRAME .or. NSTEPPLUS.eq.1) then
         PPLUS = 1.0d0
      else
         PPLUS = MIN(MIN(NPATHFRAME,NFMAX)*1.0D0/(1.0D0*(NSTEPPLUS-1)),1.0D0)
      endif

      if (NSTEPMINUS.lt.NPATHFRAME .or. NSTEPMINUS.eq.1) then
         PMINUS = 1.0d0
      else
         PMINUS=MIN(MIN(NPATHFRAME,NFMAX)*1.0D0/(1.0D0*(NSTEPMINUS-1)),1.0D0)
      endif

      IF (PTEST) WRITE(*,'(A,F10.4,A,F10.4,A)') ' Frames will be dumped to points.path.xyz with probability ',
     1                          PPLUS,'/',PMINUS,' steps on the plus/minus sides'
      NFPLUS=0
      DO J1=1,NSTEPPLUS
         ! {{{
         READ(2,*) EOFS(NSTEPPLUS+1-J1)
         IF (ZSYMSAVE(1:1).EQ.'W') THEN  
C           READ(1,*) (Q2(J2),J2=1,9*NATOMS)  !  WCOMMENT
            READ(1,*) (Q2(J2),J2=1,9*(NATOMS/2))
         ELSE
            READ(1,*) (Q2(J2),J2=1,NOPT)
         ENDIF

         ETEST=.FALSE.
         IF ((J1.EQ.1).OR.(J1.EQ.NSTEPPLUS)) THEN ! always take the end points
            ETEST=.TRUE.
         ELSE
            IF (ABS(EOFS(NSTEPPLUS+1-J1)-LASTE).GE.FRAMEEDIFF) ETEST=.TRUE. ! energy must have changed enough
         ENDIF
         IF (ETEST) LASTE=EOFS(NSTEPPLUS+1-J1)
         IF (((DPRAND().LE.PPLUS).OR.(J1.EQ.NSTEPPLUS)).AND.(NFPLUS.LT.NFMAX).AND.ETEST) THEN
            IF ((J1.EQ.NSTEPPLUS).OR.(NFPLUS.LT.NFMAX-1)) THEN ! save room for the last endpoint
               NFPLUS=NFPLUS+1
               EOFSFRAMEP(NFPLUS)=EOFS(NSTEPPLUS+1-J1)
               LASTE=EOFS(NSTEPPLUS+1-J1)
               IF (ZSYMSAVE(1:1).EQ.'W') THEN
C                 DO J2=1,9*NATOMS  !  WCOMMENT
                  DO J2=1,9*(NATOMS/2)
                     QFRAMEP(J2,NFPLUS)=Q2(J2)
                  ENDDO
               ELSE IF (ZSYMSAVE(1:2).EQ.'CD') THEN
                   DO J2=1,NATOMS/2
                      CALL CAPSIDIO(Q2(3*(J2-1)+1),Q2(3*(J2-1)+2),Q2(3*(J2-1)+3),
     1                              Q2(3*(NATOMS/2+J2-1)+1),Q2(3*(NATOMS/2+J2-1)+2),Q2(3*(NATOMS/2+J2-1)+3),
     2                              CAPSCOORDS2,RAD,HEIGHT)
                      DO J3=1,18
                         QFRAMEP(18*(J2-1)+J3,NFPLUS)=CAPSCOORDS2(J3)
                      ENDDO
                  ENDDO
               ELSE
                  DO J2=1,NOPT
                     QFRAMEP(J2,NFPLUS)=Q2(J2)
                  ENDDO
               ENDIF
!              PRINT '(A,I6,A,I6,A,G20.10)','dumping plus frame ',J1,' NFPLUS=',NFPLUS,' energy=',EOFSFRAMEP(NFPLUS)
            ENDIF
         ENDIF
         TEMP=0.0D0
         IF (BULKT) THEN
            DO J2=1,NATOMS
               TEMP=TEMP+MINIM(Q2(3*(J2-1)+1),Q1(3*(J2-1)+1),PARAM1)**2
     1                  +MINIM(Q2(3*(J2-1)+2),Q1(3*(J2-1)+2),PARAM2)**2
               IF (.NOT.TWOD) TEMP=TEMP+MINIM(Q2(3*(J2-1)+3),Q1(3*(J2-1)+3),PARAM3)**2
            ENDDO
         ELSE
            IF (ZSYMSAVE(1:1).EQ.'W') THEN
C              DO J2=1,9*NATOMS  !  WCOMMENT
               DO J2=1,9*(NATOMS/2)
                  TEMP=TEMP+(Q2(J2)-Q1(J2))**2
               ENDDO
            ELSE IF (ZSYMSAVE.EQ.'CD') THEN
               DO J2=1,NATOMS/2
                  CALL CAPSIDIO(Q2(3*(J2-1)+1),Q2(3*(J2-1)+2),Q2(3*(J2-1)+3),
     1                          Q2(3*(NATOMS/2+J2-1)+1),Q2(3*(NATOMS/2+J2-1)+2),Q2(3*(NATOMS/2+J2-1)+3),
     2                          CAPSCOORDS2,RAD,HEIGHT)
                  CALL CAPSIDIO(Q1(3*(J2-1)+1),Q1(3*(J2-1)+2),Q1(3*(J2-1)+3),
     1                          Q1(3*(NATOMS/2+J2-1)+1),Q1(3*(NATOMS/2+J2-1)+2),Q1(3*(NATOMS/2+J2-1)+3),
     2                          CAPSCOORDS1,RAD,HEIGHT)
                  DO J3=1,18
                      TEMP=TEMP+(CAPSCOORDS1(J3)-CAPSCOORDS2(J3))**2
                  ENDDO
               ENDDO
            ELSE IF((PYGPERIODICT.OR.PYBINARYT).AND.UNIAXT) THEN
C           calculate path lengths wrt. Cartesian coordinate of the centre and two coordinates along the symmetry axis only
               CALL UNIAXGETPATHLENGTH(Q1,Q2,TEMP)
            ELSE
               DO J2=1,NOPT
                  TEMP=TEMP+(Q2(J2)-Q1(J2))**2
               ENDDO
            ENDIF
         ENDIF
         PATHLENGTH(NSTEPPLUS+1-J1)=PATHLENGTH(NSTEPPLUS+2-J1)-SQRT(TEMP)
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
C           DO J2=1,9*NATOMS  !  WCOMMENT
            DO J2=1,9*(NATOMS/2)
               Q1(J2)=Q2(J2)
            ENDDO
         ELSE
            DO J2=1,NOPT
               Q1(J2)=Q2(J2)
            ENDDO
         ENDIF
         ! }}}
      ENDDO
      IF (PTEST) WRITE(*,'(A,I6,A,I6,A)') ' Transition state will be frame number ',NFPLUS+1

      IF (ZSYMSAVE(1:1).EQ.'W') THEN
         ! {{{
C        DO J1=1,NATOMS  !  WCOMMENT
         DO J1=1,NATOMS/2
            CALL CONVERT(QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3),
C    1                   QINIT(3*(NATOMS+J1-1)+1),QINIT(3*(NATOMS+J1-1)+2),QINIT(3*(NATOMS+J1-1)+3),
     1                   QINIT(3*(NATOMS/2+J1-1)+1),QINIT(3*(NATOMS/2+J1-1)+2),QINIT(3*(NATOMS/2+J1-1)+3),
     2                   OVEC,H1VEC,H2VEC)
            Q1(9*(J1-1)+1)=OVEC(1)
            Q1(9*(J1-1)+2)=OVEC(2)
            Q1(9*(J1-1)+3)=OVEC(3)
            Q1(9*(J1-1)+4)=H1VEC(1)
            Q1(9*(J1-1)+5)=H1VEC(2)
            Q1(9*(J1-1)+6)=H1VEC(3)
            Q1(9*(J1-1)+7)=H2VEC(1)
            Q1(9*(J1-1)+8)=H2VEC(2)
            Q1(9*(J1-1)+9)=H2VEC(3)
         ENDDO
         ! }}}
      ELSE
         ! {{{
         DO J1=1,NOPT
            Q1(J1)=QINIT(J1)
         ENDDO
         ! }}}
      ENDIF

      NFMINUS=0
      DO J1=1,NSTEPMINUS
      ! {{{
         READ(2,*) EOFS(NSTEPPLUS+1+J1)
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
C           READ(1,*) (Q2(J2),J2=1,9*NATOMS)  !  WCOMMENT
            READ(1,*) (Q2(J2),J2=1,9*(NATOMS/2))
         ELSE
            READ(1,*) (Q2(J2),J2=1,NOPT)
         ENDIF

         ETEST=.FALSE.
         IF ((J1.EQ.1).OR.(J1.EQ.NSTEPMINUS)) THEN ! always take the end points
            ETEST=.TRUE.
         ELSE
            IF (ABS(EOFS(NSTEPPLUS+1+J1)-LASTE).GE.FRAMEEDIFF) ETEST=.TRUE. ! energy must have changed enough
         ENDIF
         IF (ETEST) LASTE=EOFS(NSTEPPLUS+1+J1)
         RANDOM=DPRAND()
         IF (((RANDOM.LE.PMINUS).OR.(J1.EQ.NSTEPMINUS)).AND.(NFMINUS.LE.NFMAX).AND.ETEST) THEN
            IF ((J1.EQ.NSTEPMINUS).OR.(NFMINUS.LT.NFMAX-1)) THEN ! save a space for the stationary point at the end
               NFMINUS=NFMINUS+1
               EOFSFRAMEM(NFMINUS)=EOFS(NSTEPPLUS+1+J1)
               LASTE=EOFS(NSTEPPLUS+1+J1)
               IF (ZSYMSAVE(1:1).EQ.'W') THEN
C                 DO J2=1,9*NATOMS  !  WCOMMENT
                     DO J2=1,9*(NATOMS/2)
                     QFRAMEM(J2,NFMINUS)=Q2(J2)
                  ENDDO
               ELSE IF (ZSYMSAVE.EQ.'CD') THEN
                   DO J2=1,NATOMS/2
                      CALL CAPSIDIO(Q2(3*(J2-1)+1),Q2(3*(J2-1)+2),Q2(3*(J2-1)+3),
     1                              Q2(3*(NATOMS/2+J2-1)+1),Q2(3*(NATOMS/2+J2-1)+2),Q2(3*(NATOMS/2+J2-1)+3),
     2                              CAPSCOORDS2,RAD,HEIGHT)
                      DO J3=1,18
                         QFRAMEM(18*(J2-1)+J3,NFMINUS)=CAPSCOORDS2(J3)
                      ENDDO
                  ENDDO
               ELSE
                  DO J2=1,NOPT
                     QFRAMEM(J2,NFMINUS)=Q2(J2)
                  ENDDO
               ENDIF
!              PRINT '(A,I6,A,I6,A,G20.10)','dumping minus frame ',J1,' NFMINUS=',NFMINUS,' energy=',EOFSFRAMEM(NFMINUS)

            ENDIF
         ENDIF
         TEMP=0.0D0
         IF (BULKT) THEN
            DO J2=1,NATOMS
               TEMP=TEMP+MINIM(Q2(3*(J2-1)+1),Q1(3*(J2-1)+1),PARAM1)**2
     1                  +MINIM(Q2(3*(J2-1)+2),Q1(3*(J2-1)+2),PARAM2)**2
               IF (.NOT.TWOD) TEMP=TEMP+MINIM(Q2(3*(J2-1)+3),Q1(3*(J2-1)+3),PARAM3)**2
            ENDDO
         ELSE
            IF (ZSYMSAVE(1:1).EQ.'W') THEN
C              DO J2=1,9*NATOMS  !  WCOMMENT
               DO J2=1,9*(NATOMS/2)
                  TEMP=TEMP+(Q2(J2)-Q1(J2))**2
               ENDDO
            ELSE IF (ZSYMSAVE.EQ.'CD') THEN
               DO J2=1,NATOMS/2
                  CALL CAPSIDIO(Q2(3*(J2-1)+1),Q2(3*(J2-1)+2),Q2(3*(J2-1)+3),
     1                          Q2(3*(NATOMS/2+J2-1)+1),Q2(3*(NATOMS/2+J2-1)+2),Q2(3*(NATOMS/2+J2-1)+3),
     2                          CAPSCOORDS2,RAD,HEIGHT)
                  CALL CAPSIDIO(Q1(3*(J2-1)+1),Q1(3*(J2-1)+2),Q1(3*(J2-1)+3),
     1                          Q1(3*(NATOMS/2+J2-1)+1),Q1(3*(NATOMS/2+J2-1)+2),Q1(3*(NATOMS/2+J2-1)+3),
     2                          CAPSCOORDS1,RAD,HEIGHT)
                  DO J3=1,18
                      TEMP=TEMP+(CAPSCOORDS1(J3)-CAPSCOORDS2(J3))**2
                  ENDDO
               ENDDO
            ELSE IF((PYGPERIODICT.OR.PYBINARYT).AND.UNIAXT) THEN
C           calculate path lengths wrt. Cartesian coordinate of the centre and two coordinates along the symmetry axis only
               CALL UNIAXGETPATHLENGTH(Q1,Q2,TEMP)
            ELSE
               DO J2=1,NOPT
                  TEMP=TEMP+(Q2(J2)-Q1(J2))**2
               ENDDO
            ENDIF
         ENDIF
         PATHLENGTH(NSTEPPLUS+1+J1)=PATHLENGTH(NSTEPPLUS+J1)+SQRT(TEMP)
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
C           DO J2=1,9*NATOMS  !  WCOMMENT
            DO J2=1,9*(NATOMS/2)
               Q1(J2)=Q2(J2)
            ENDDO
         ELSE
            DO J2=1,NOPT
               Q1(J2)=Q2(J2)
            ENDDO
         ENDIF
         ! }}}
      ENDDO

555   CONTINUE
      OPEN(UNIT=3,FILE=EOFSSTRING,STATUS='UNKNOWN')
      WRITE(3,'(2G20.10,I6)') (PATHLENGTH(J1),EOFS(J1),J1,J1=1,NSTEPPLUS+NSTEPMINUS+1)
      CLOSE(3)
      
      SUM2=0.0D0
      SUM4=0.0D0
C     NDUMMY=1 !  WCOMMENT
C     IF (ZSYMSAVE(1:1).EQ.'W') NDUMMY=9  !  WCOMMENT 9 should have been 3?
      NATOMSIMUL=NATOMS
      IF (ZSYMSAVE(1:1).EQ.'W') NATOMSIMUL=3*(NATOMS/2)
      IF (ZSYMSAVE(1:1).EQ.'W') THEN
        ! {{{
         DO J1=1,NATOMSIMUL  !  WCOMMENT
            SUM2=SUM2+
     1           (QFRAMEP(3*(J1-1)+1,NFPLUS)-QFRAMEM(3*(J1-1)+1,NFMINUS))**2
     2          +(QFRAMEP(3*(J1-1)+2,NFPLUS)-QFRAMEM(3*(J1-1)+2,NFMINUS))**2
     3          +(QFRAMEP(3*(J1-1)+3,NFPLUS)-QFRAMEM(3*(J1-1)+3,NFMINUS))**2
            SUM4=SUM4+
     1          ((QFRAMEP(3*(J1-1)+1,NFPLUS)-QFRAMEM(3*(J1-1)+1,NFMINUS))**2
     2          +(QFRAMEP(3*(J1-1)+2,NFPLUS)-QFRAMEM(3*(J1-1)+2,NFMINUS))**2
     3          +(QFRAMEP(3*(J1-1)+3,NFPLUS)-QFRAMEM(3*(J1-1)+3,NFMINUS))**2)**2
         ENDDO
         ! }}}
      ELSE IF (BULKT) THEN
        ! {{{
         IF (TWOD) THEN
            DO J1=1,NATOMSIMUL  !  WCOMMENT
               SUM2=SUM2+
     1              MINIM(QPLUS(3*(J1-1)+1),Q(3*(J1-1)+1),PARAM1)**2
     2             +MINIM(QPLUS(3*(J1-1)+2),Q(3*(J1-1)+2),PARAM2)**2
               SUM4=SUM4+
     1             (MINIM(QPLUS(3*(J1-1)+1),Q(3*(J1-1)+1),PARAM1)**2
     2             +MINIM(QPLUS(3*(J1-1)+2),Q(3*(J1-1)+2),PARAM2)**2)**2
            ENDDO
         ELSE
C           DO J1=1,NATOMS*NDUMMY
            DO J1=1,NATOMSIMUL  !  WCOMMENT
               SUM2=SUM2+
     1              MINIM(QPLUS(3*(J1-1)+1),Q(3*(J1-1)+1),PARAM1)**2
     2             +MINIM(QPLUS(3*(J1-1)+2),Q(3*(J1-1)+2),PARAM2)**2
     3             +MINIM(QPLUS(3*(J1-1)+3),Q(3*(J1-1)+3),PARAM3)**2
               SUM4=SUM4+
     1             (MINIM(QPLUS(3*(J1-1)+1),Q(3*(J1-1)+1),PARAM1)**2
     2             +MINIM(QPLUS(3*(J1-1)+2),Q(3*(J1-1)+2),PARAM2)**2
     3             +MINIM(QPLUS(3*(J1-1)+3),Q(3*(J1-1)+3),PARAM3)**2)**2
            ENDDO
         ENDIF
         ! }}}
      ELSE IF (ZSYMSAVE.EQ.'CD') THEN
        ! {{{
         DO J2=1,NATOMS/2
            CALL CAPSIDIO(QPLUS(3*(J2-1)+1),QPLUS(3*(J2-1)+2),QPLUS(3*(J2-1)+3),
     1                    QPLUS(3*(NATOMS/2+J2-1)+1),QPLUS(3*(NATOMS/2+J2-1)+2),QPLUS(3*(NATOMS/2+J2-1)+3),
     2                    CAPSCOORDS2,RAD,HEIGHT)
            CALL CAPSIDIO(Q(3*(J2-1)+1),Q(3*(J2-1)+2),Q(3*(J2-1)+3),
     1                    Q(3*(NATOMS/2+J2-1)+1),Q(3*(NATOMS/2+J2-1)+2),Q(3*(NATOMS/2+J2-1)+3),
     2                    CAPSCOORDS1,RAD,HEIGHT)
            DO J1=1,18
               SUM2=SUM2+
     1              (CAPSCOORDS1(3*(J1-1)+1)-CAPSCOORDS2(3*(J1-1)+1))**2
     2             +(CAPSCOORDS1(3*(J1-1)+2)-CAPSCOORDS2(3*(J1-1)+2))**2
     3             +(CAPSCOORDS1(3*(J1-1)+3)-CAPSCOORDS2(3*(J1-1)+3))**2
               SUM4=SUM4+
     1             ((CAPSCOORDS1(3*(J1-1)+1)-CAPSCOORDS2(3*(J1-1)+1))**2
     2             +(CAPSCOORDS1(3*(J1-1)+2)-CAPSCOORDS2(3*(J1-1)+2))**2
     3             +(CAPSCOORDS1(3*(J1-1)+3)-CAPSCOORDS2(3*(J1-1)+3))**2)**2
            ENDDO
         ENDDO
         ! }}}
      ELSEIF (RINGPOLYMERT) THEN
      ! {{{
         SUM2=1.0D0
         SUM4=1.0D0
         ! }}}
      ELSE
        ! {{{
C        DO J1=1,NATOMS*NDUMMY  !  WCOMMENT
         DO J1=1,NATOMSIMUL
            SUM2=SUM2+
     1           (QPLUS(3*(J1-1)+1)-Q(3*(J1-1)+1))**2
     2          +(QPLUS(3*(J1-1)+2)-Q(3*(J1-1)+2))**2
     3          +(QPLUS(3*(J1-1)+3)-Q(3*(J1-1)+3))**2
            SUM4=SUM4+
     1          ((QPLUS(3*(J1-1)+1)-Q(3*(J1-1)+1))**2
     2          +(QPLUS(3*(J1-1)+2)-Q(3*(J1-1)+2))**2
     3          +(QPLUS(3*(J1-1)+3)-Q(3*(J1-1)+3))**2)**2
         ENDDO
         ! }}}
      ENDIF

      ETS=EOFS(NSTEPPLUS+1)
      EPLUS=EOFS(1)
      EMINUS=EOFS(NSTEPPLUS+NSTEPMINUS+1)
      PRINT*

      IF (RINGPOLYMERT) THEN
        ! {{{
         WRITE(*,'(A)')
     1'         E+        Ets - E+           Ets       Ets - E-           E-'
         WRITE(*,60) EPLUS,ETS-EPLUS,ETS,ETS-EMINUS,EMINUS
      ELSE
         WRITE(*,'(A)')
     1'         E+        Ets - E+           Ets       Ets - E-           E-          S       D' //
     2 '      gamma   ~N'
         WRITE(*,60) EPLUS,ETS-EPLUS,ETS,ETS-EMINUS,EMINUS,
     1         PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)-PATHLENGTH(1),SQRT(SUM2),SUM4*NATOMSIMUL/SUM2**2,SUM2**2/SUM4
60       FORMAT(F17.7,G12.5,F17.7,G12.5,F17.7,F8.3,F8.3,F8.3,F8.3)
         ! }}}
      ENDIF
C
C tvb Calculation of catastrophe ratios
      IF (RATIOS) THEN
            ! {{{
            CALL POTENTIAL(QPLUS,EPLUS,VNEW,.true.,.TRUE.,RMS,.FALSE.,.FALSE.)
            CALL DSYEV('V','U',NOPT,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO) 
            call eigensort_val_asc(diag,hess,nopt,3*natoms)
            lambdap=diag(3*natoms-6)
            CALL POTENTIAL(QMINUS,EMINUS,VNEW,.true.,.TRUE.,RMS,.FALSE.,.FALSE.)
            CALL DSYEV('N','U',NOPT,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO) 
            call eigensort_val_asc(diag,hess,nopt,3*natoms)
            lambdam=diag(3*natoms-6)
            CALL POTENTIAL(QINIT,EOFS(NSTEPPLUS+1),VNEW,.true.,.TRUE.,RMS,.FALSE.,.FALSE.)
            CALL DSYEV('N','U',NOPT,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO) 
            call eigensort_val_asc(diag,hess,nopt,3*natoms)
            lambdats=diag(3*natoms)
            CALL NEWMINDIST(QPLUS,QINIT,NATOMS,DISTP,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
            CALL NEWMINDIST(QMINUS,QINIT,NATOMS,DISTM,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
            OPEN(UNIT=24,FILE='folddata',STATUS='UNKNOWN')
C           Folddata: Em1, Em2, Ets, Sm1, Sm2, Evts, Ev1, Ev2, Sm1min, Sm2min
            WRITE(24,'(10F13.7)') Eofs(1), Eofs(nstepplus+nstepminus+1), Eofs(nstepplus+1), PATHLENGTH(1),
     1      PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1), lambdats, lambdap, lambdam, distp, distm
            ! }}}
      ENDIF
C end tvb
C
      SLENGTH=PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)-PATHLENGTH(1)
      DISP=SQRT(SUM2)
C     GAMMA=SUM4*NDUMMY*NATOMS/SUM2**2  !  WCOMMENT
      GAMMA=SUM4*NATOMSIMUL/SUM2**2
      NTILDE=SUM2**2/SUM4
      IF (CHECKINDEX) THEN
        ! {{{
         SMINUS=PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)
         SPLUS=PATHLENGTH(1)
         STS=PATHLENGTH(NSTEPPLUS+1)
         IF (EPLUS.LT.EMINUS) THEN
            ETEMP=EPLUS
            STEMP=SPLUS
            DUMMY=EVPLUS
            EPLUS=EMINUS
            SPLUS=SMINUS
            EVPLUS=EVMINUS
            EMINUS=ETEMP
            SMINUS=STEMP
            EVMINUS=DUMMY
         ENDIF
         PRINT*
         WRITE(*,'(A)')
     1      '         evts         evplus/ts    evminus/ts  del e plus      del e minus      s plus         s minus          frat +'
     2        // '    frat -  '

         WRITE(*,'(A5,3F13.7,4G16.8,2F10.4)') 'fold ',EVTS,EVPLUS/ABS(EVTS),EVMINUS/ABS(EVTS),ETS-EPLUS,ETS-EMINUS,ABS(STS-SPLUS),
     1              ABS(SMINUS-STS),6*(ETS-EPLUS)/(ABS(EVTS)*(STS-SPLUS)**2),6*(ETS-EMINUS)/(ABS(EVTS)*(STS-SMINUS)**2)

         WRITE(*,'(A)') '        del e plus      del e minus     s plus           s minus       '
     1        // '   evrat     dev1      dev2      dev3      rat+      rat-'

         WRITE(*,'(A5,4G16.8,6F10.4)') 'cusp ',ETS-EPLUS,ETS-EMINUS,ABS(STS-SPLUS),ABS(SMINUS-STS),
     1                           EVPLUS*EVMINUS/((EVPLUS+EVMINUS)*EVTS),
     2                           (ETS-EPLUS)*EVMINUS**3*(EVMINUS+2.0D0*EVPLUS)/((ETS-EMINUS)*EVPLUS**3*(EVPLUS+2.0D0*EVMINUS)),
     3                           (ETS-EPLUS)*EVTS**3*(EVTS+2.0D0*EVPLUS)/((ETS-EMINUS)*(EVTS-EVPLUS)*(EVPLUS+EVTS)**3),
     4                           (ETS-EPLUS)*(EVTS-EVMINUS)*(EVMINUS+EVTS)**3/((ETS-EMINUS)*EVTS**3*(2.0D0*EVMINUS+EVTS)),
     5                           12*(ETS-EPLUS) /((EVPLUS -EVTS)*(STS-SPLUS )**2),
     6                           12*(ETS-EMINUS)/((EVMINUS-EVTS)*(STS-SMINUS)**2)
         ! }}}
      ENDIF

      OPEN(UNIT=3,FILE=ITSTRING,STATUS='UNKNOWN')

      ! ZSYMSAVE == W CD; ELSE  {{{
C     IF (ZSYMSAVE(1:1).EQ.'W') THEN
         ! {{{
C        WRITE(3,'(I6)') 3*NATOMS
C        WRITE(3,'(g20.10)') EOFS(1)
C        DO J1=1,NATOMS
C           CALL CONVERT(QPLUS(3*(J1-1)+1),QPLUS(3*(J1-1)+2),QPLUS(3*(J1-1)+3),
C    1                   QPLUS(3*(NATOMS+J1-1)+1),QPLUS(3*(NATOMS+J1-1)+2),QPLUS(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
C           WRITE(3,'(A2,4X,3g20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
C        ENDDO
         ! }}}
C     ELSE IF (ZSYMSAVE.EQ.'CD') THEN
         ! {{{
C        WRITE(3,'(I6)') NATOMS*6/2
C        WRITE(3,'(g20.10)') EOFS(1)
C        DO J2=1,NATOMS/2
C           CALL CAPSIDIO(QPLUS(3*(J2-1)+1),QPLUS(3*(J2-1)+2),QPLUS(3*(J2-1)+3),
C    1                    QPLUS(3*(NATOMS/2+J2-1)+1),QPLUS(3*(NATOMS/2+J2-1)+2),QPLUS(3*(NATOMS/2+J2-1)+3),CAPSCOORDS2,RAD,HEIGHT)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(1),CAPSCOORDS2(2),CAPSCOORDS2(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(4),CAPSCOORDS2(5),CAPSCOORDS2(6)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(7),CAPSCOORDS2(8),CAPSCOORDS2(9)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(10),CAPSCOORDS2(11),CAPSCOORDS2(12)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(13),CAPSCOORDS2(14),CAPSCOORDS2(15)
C           WRITE(3,'(A2,4X,3g20.10)') 'C4  ',CAPSCOORDS2(16),CAPSCOORDS2(17),CAPSCOORDS2(18)
C        ENDDO
         ! }}}
C     ELSE
         ! {{{
C        WRITE(3,'(I6)') NATOMS
C        WRITE(3,'(g20.10)') EOFS(1)
C        WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J1),QPLUS(3*(J1-1)+1),QPLUS(3*(J1-1)+2),QPLUS(3*(J1-1)+3),J1=1,NATOMS)
         ! }}}
C     ENDIF
         ! }}}

      DO J1=NFPLUS,1,-1
         ! {{{
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
           ! {{{
C           WRITE(3,'(I6)') 3*NATOMS  !  WCOMMENT
            WRITE(3,'(I6)') 3*(NATOMS/2)
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEP(J1)
C           IF (J1.EQ.NFPLUS) THEN 
C              WRITE(3,'(g20.10)') EOFS(1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A2,4X,3g20.10)') 
     1   ('O  ',QFRAMEP(9*(J2-1)+1,J1),QFRAMEP(9*(J2-1)+2,J1),QFRAMEP(9*(J2-1)+3,J1),
     1    'H  ',QFRAMEP(9*(J2-1)+4,J1),QFRAMEP(9*(J2-1)+5,J1),QFRAMEP(9*(J2-1)+6,J1),
C    2    'H  ',QFRAMEP(9*(J2-1)+7,J1),QFRAMEP(9*(J2-1)+8,J1),QFRAMEP(9*(J2-1)+9,J1),J2=1,NATOMS)
     2    'H  ',QFRAMEP(9*(J2-1)+7,J1),QFRAMEP(9*(J2-1)+8,J1),QFRAMEP(9*(J2-1)+9,J1),J2=1,NATOMS/2)
            ! }}}
         ELSE IF (ZSYMSAVE.EQ.'CD') THEN
           ! {{{
            WRITE(3,'(I6)') NATOMS*6/2
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEP(J1)
C           IF (J1.EQ.NFPLUS) THEN
C              WRITE(3,'(g20.10)') EOFS(1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A)') ' '
            DO J2=1,NATOMS/2
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+1,J1),QFRAMEP(18*(J2-1)+2,J1),QFRAMEP(18*(J2-1)+3,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+4,J1),QFRAMEP(18*(J2-1)+5,J1),QFRAMEP(18*(J2-1)+6,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+7,J1),QFRAMEP(18*(J2-1)+8,J1),QFRAMEP(18*(J2-1)+9,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+10,J1),QFRAMEP(18*(J2-1)+11,J1),QFRAMEP(18*(J2-1)+12,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+13,J1),QFRAMEP(18*(J2-1)+14,J1),QFRAMEP(18*(J2-1)+15,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C4  ',QFRAMEP(18*(J2-1)+16,J1),QFRAMEP(18*(J2-1)+17,J1),QFRAMEP(18*(J2-1)+18,J1)
            ENDDO
            ! }}}
         ELSEIF (STOCKT) THEN
            ! {{{
            WRITE(3,'(I6)') (NATOMS/2)
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
            DO J2=1,(NATOMS/2)
               WRITE(3,'(A2,4X,3G20.10,A13,3G20.10)')
     &         ZSYM(J2),QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1),
     &                        ' atom_vector ',
     &                        SIN(QFRAMEP(3*((NATOMS/2)+J2-1)+1,J1))*COS(QFRAMEP(3*((NATOMS/2)+J2-1)+2,J1)),
     &                        SIN(QFRAMEP(3*((NATOMS/2)+J2-1)+1,J1))*SIN(QFRAMEP(3*((NATOMS/2)+J2-1)+2,J1)),
     &                        COS(QFRAMEP(3*((NATOMS/2)+J2-1)+1,J1))
            ENDDO

            ! }}}
!         ELSEIF (STOCKAAT) THEN
            ! {{{
!            WRITE(3,'(I6)') (NATOMS/2)
!            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
!            DO J2=1,(NATOMS/2)
               
!               P(:) = QFRAMEP(3*((NATOMS/2)+J2-1)+1:3*((NATOMS/2)+J2-1)+3,J1)
!               CALL ROTMAT(P(:), RMAT(:,:))
!               P(:) = RMAT(:,3)     ! THE DIPOLE-VECTOR HAS TO BE ALONG THE Z-AXIS IN THE BODY-FRAME

!               WRITE(3,'(A1,4X,3G20.10,A13,3G20.10)')
!     &         'O', QFRAMEP(3*(J2-1)+1,J1), QFRAMEP(3*(J2-1)+2,J1), QFRAMEP(3*(J2-1)+3,J1),
!     &         ' atom_vector ', P(1), P(2), P(3)
!            ENDDO
        ! }}}
         ELSEIF (RBAAT) THEN
         ! {{{
            WRITE(3,'(I6)') (NATOMS/2)
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
            DO J2=1,(NATOMS/2)
               WRITE(3,'(A1,4X,6G20.10)')
     &         'O', QFRAMEP(3*(J2-1)+1,J1), QFRAMEP(3*(J2-1)+2,J1), QFRAMEP(3*(J2-1)+3,J1),
     &         QFRAMEP(3*((NATOMS/2)+J2-1)+1,J1), QFRAMEP(3*((NATOMS/2)+J2-1)+2,J1),
     &         QFRAMEP(3*((NATOMS/2)+J2-1)+3,J1) 
            ENDDO
            ! }}}
         ELSEIF (AMHT) THEN
         ! {{{
            GLY_COUNT = 0
            DO J2=1,NMRES
               IF (SEQ(J2).EQ.8) GLY_COUNT = GLY_COUNT +1
            ENDDO
            WRITE(3,'(I6)') NATOMS+GLY_COUNT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)

            GLY_COUNT = 0
            DO J2=1,NMRES
              IF (SEQ(J2).EQ.8) THEN
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEP(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEP(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'O    ',QFRAMEP(9*(J2-1)+4-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+5-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+6-GLY_COUNT*3,J1)
                GLY_COUNT = GLY_COUNT +1
              ELSE
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEP(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'C2   ',QFRAMEP(9*(J2-1)+4-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+5-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+6-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'O    ',QFRAMEP(9*(J2-1)+7-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+8-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+9-GLY_COUNT*3,J1)
              ENDIF
          ENDDO
          ! }}}
         ELSEIF (RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')) THEN
         ! {{{
            WRITE(3,'(I6)') NOPT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
            WRITE(3,'(A2,4X,3G20.10)') (ZSYM(J2),QFRAMEP(J2,J1),0.0D0,0.0D0,J2=1,NOPT)
            ! }}}
         ELSE
           ! {{{
            WRITE(3,'(I6)') NATOMS
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
C           IF (J1.EQ.NFPLUS) THEN
C              WRITE(3,'(g20.10)') EOFS(1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            if (ZSYM(NATOMS).eq.'SV') then
               do j2=1, NATOMS
                  if (MOD(j2,3).eq.0) then
                     WRITE(3,'(A,4X,3g20.10)') 'O',QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1)
                  else
                     WRITE(3,'(A,4X,3g20.10)') 'H',QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1)
                  endif
               enddo
            ELSE
               WRITE(3,'(A2,4X,3G20.10)') 
     &         (ZSYM(J2),QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1),J2=1,NATOMS)
            ENDIF
            ! }}}
         ENDIF
         ! }}}
      ENDDO
!
!  ZSYMSAVE=W,CD  
!  STOCKT RBAAT AMHT RINGPOLYMERT ELSE
!  {{{
      IF (ZSYMSAVE(1:1).EQ.'W') THEN
        ! {{{
C        WRITE(3,'(I6)') 3*NATOMS ! WCOMMENT
         WRITE(3,'(I6)') 3*(NATOMS/2)
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
C        DO J1=1,NATOMS ! WCOMMENT
         DO J1=1,NATOMS/2
            CALL CONVERT(QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3),
C    1                   QINIT(3*(NATOMS+J1-1)+1),QINIT(3*(NATOMS+J1-1)+2),QINIT(3*(NATOMS+J1-1)+3),
     1                   QINIT(3*(NATOMS/2+J1-1)+1),QINIT(3*(NATOMS/2+J1-1)+2),QINIT(3*(NATOMS/2+J1-1)+3),
     2                   OVEC,H1VEC,H2VEC)
            WRITE(3,'(A2,4X,3g20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
            WRITE(3,'(A2,4X,3g20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
            WRITE(3,'(A2,4X,3g20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
         ENDDO
         ! }}}
      ELSE IF (ZSYMSAVE.EQ.'CD') THEN
        ! {{{
         WRITE(3,'(I6)') NATOMS*6/2
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
         DO J2=1,NATOMS/2
            CALL CAPSIDIO(QINIT(3*(J2-1)+1),QINIT(3*(J2-1)+2),QINIT(3*(J2-1)+3),
     1                    QINIT(3*(NATOMS/2+J2-1)+1),QINIT(3*(NATOMS/2+J2-1)+2),QINIT(3*(NATOMS/2+J2-1)+3),
     2                    CAPSCOORDS2,RAD,HEIGHT)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(1),CAPSCOORDS2(2),CAPSCOORDS2(3)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(4),CAPSCOORDS2(5),CAPSCOORDS2(6)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(7),CAPSCOORDS2(8),CAPSCOORDS2(9)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(10),CAPSCOORDS2(11),CAPSCOORDS2(12)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(13),CAPSCOORDS2(14),CAPSCOORDS2(15)
            WRITE(3,'(A2,4X,3g20.10)') 'C4  ',CAPSCOORDS2(16),CAPSCOORDS2(17),CAPSCOORDS2(18)
         ENDDO
         ! }}}
      ELSEIF (STOCKT) THEN
         ! {{{
         WRITE(3,'(I6)') (NATOMS/2)
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
         DO J2=1,(NATOMS/2)
            WRITE(3,'(A2,4X,3G20.10,A13,3G20.10)')
     &         ZSYM(J2),QINIT(3*(J2-1)+1),QINIT(3*(J2-1)+2),QINIT(3*(J2-1)+3),
     &                                  ' atom_vector ',
     &                                  SIN(QINIT(3*((NATOMS/2)+J2-1)+1))*COS(QINIT(3*((NATOMS/2)+J2-1)+2)),
     &                                  SIN(QINIT(3*((NATOMS/2)+J2-1)+1))*SIN(QINIT(3*((NATOMS/2)+J2-1)+2)),
     &                                  COS(QINIT(3*((NATOMS/2)+J2-1)+1))
         ENDDO
         ! }}}
      ELSEIF (RBAAT) THEN
         ! {{{
         WRITE(3,'(I6)') (NATOMS/2)
         WRITE(3,'(A,G25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
         DO J2=1,(NATOMS/2)
            WRITE(3,'(A1,4X,6G20.10)')
     &      'O', QINIT(3*(J2-1)+1), QINIT(3*(J2-1)+2), QINIT(3*(J2-1)+3),
     &       QINIT(3*((NATOMS/2)+J2-1)+1), QINIT(3*((NATOMS/2)+J2-1)+2),
     &       QINIT(3*((NATOMS/2)+J2-1)+3)
         ENDDO
         ! }}}
      ELSEIF (AMHT) THEN
         ! {{{
            GLY_COUNT = 0
            DO J2=1,NMRES
               IF (SEQ(J2).EQ.8) GLY_COUNT = GLY_COUNT +1
            ENDDO
            WRITE(3,'(I6)') NATOMS+GLY_COUNT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
            GLY_COUNT = 0
            DO J2=1,NMRES
              IF (SEQ(J2).EQ.8) THEN
        WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QINIT(9*(J2-1)+1-GLY_COUNT*3),QINIT(9*(J2-1)+2-GLY_COUNT*3),QINIT(9*(J2-1)+3-GLY_COUNT*3)
        WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QINIT(9*(J2-1)+1-GLY_COUNT*3),QINIT(9*(J2-1)+2-GLY_COUNT*3),QINIT(9*(J2-1)+3-GLY_COUNT*3)
        WRITE(3,'(a5,1x,3f20.10)') 'O    ',QINIT(9*(J2-1)+4-GLY_COUNT*3),QINIT(9*(J2-1)+5-GLY_COUNT*3),QINIT(9*(J2-1)+6-GLY_COUNT*3)
                GLY_COUNT = GLY_COUNT +1
              ELSE
        WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QINIT(9*(J2-1)+1-GLY_COUNT*3),QINIT(9*(J2-1)+2-GLY_COUNT*3),QINIT(9*(J2-1)+3-GLY_COUNT*3)
        WRITE(3,'(a5,1x,3f20.10)') 'C2   ',QINIT(9*(J2-1)+4-GLY_COUNT*3),QINIT(9*(J2-1)+5-GLY_COUNT*3),QINIT(9*(J2-1)+6-GLY_COUNT*3)
        WRITE(3,'(a5,1x,3f20.10)') 'O    ',QINIT(9*(J2-1)+7-GLY_COUNT*3),QINIT(9*(J2-1)+8-GLY_COUNT*3),QINIT(9*(J2-1)+9-GLY_COUNT*3)
              ENDIF
          ENDDO

          ! }}}
      ELSEIF (RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')) THEN
         ! {{{
         WRITE(3,'(I6)') NOPT
         WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(NSTEPPLUS+1)
         WRITE(3,'(A2,4X,3G20.10)') (ZSYM(J2),QINIT(J2),0.0D0,0.0D0,J2=1,NOPT)

         ! }}}
      ELSE
         ! {{{
         WRITE(3,'(I6)') NATOMS
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)

         if (ZSYM(NATOMS).eq.'SV') then
            do j1=1, NATOMS
               if (MOD(j1,3).eq.0) then
                  WRITE(3,'(A,4X,3g20.10)') 'O',QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3)
               else
                  WRITE(3,'(A,4X,3g20.10)') 'H',QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3)
               endif
            enddo
         else
            WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J1),QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3),J1=1,NATOMS)
         endif
         ! }}}
      ENDIF
! }}}

      DO J1=1,NFMINUS
         ! {{{
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
           ! {{{
C           WRITE(3,'(I6)') 3*NATOMS  !  WCOMMENT
            WRITE(3,'(I6)') 3*(NATOMS/2)
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
C           IF (J1.EQ.NFMINUS) THEN
C              WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A2,4X,3g20.10)') 
     1            ('O  ',QFRAMEM(9*(J2-1)+1,J1),QFRAMEM(9*(J2-1)+2,J1),QFRAMEM(9*(J2-1)+3,J1),
     1             'H  ',QFRAMEM(9*(J2-1)+4,J1),QFRAMEM(9*(J2-1)+5,J1),QFRAMEM(9*(J2-1)+6,J1),
C    2             'H  ',QFRAMEM(9*(J2-1)+7,J1),QFRAMEM(9*(J2-1)+8,J1),QFRAMEM(9*(J2-1)+9,J1),J2=1,NATOMS)
     2             'H  ',QFRAMEM(9*(J2-1)+7,J1),QFRAMEM(9*(J2-1)+8,J1),QFRAMEM(9*(J2-1)+9,J1),J2=1,NATOMS/2)
            ! }}}
         ELSE IF (ZSYMSAVE.EQ.'CD') THEN
           ! {{{
            WRITE(3,'(I6)') NATOMS*6/2
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
C           IF (J1.EQ.NFMINUS) THEN
C              WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            DO J2=1,NATOMS/2
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+1,J1),QFRAMEM(18*(J2-1)+2,J1),QFRAMEM(18*(J2-1)+3,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+4,J1),QFRAMEM(18*(J2-1)+5,J1),QFRAMEM(18*(J2-1)+6,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+7,J1),QFRAMEM(18*(J2-1)+8,J1),QFRAMEM(18*(J2-1)+9,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+10,J1),QFRAMEM(18*(J2-1)+11,J1),QFRAMEM(18*(J2-1)+12,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+13,J1),QFRAMEM(18*(J2-1)+14,J1),QFRAMEM(18*(J2-1)+15,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C4  ',QFRAMEM(18*(J2-1)+16,J1),QFRAMEM(18*(J2-1)+17,J1),QFRAMEM(18*(J2-1)+18,J1)
            ENDDO
            ! }}}
         ELSEIF (STOCKT) THEN
            ! {{{
            WRITE(3,'(I6)') (NATOMS/2)
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
            DO J2=1,(NATOMS/2)
               WRITE(3,'(A2,4X,3G20.10,A13,3G20.10)')
     &         ZSYM(J2),QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1),
     &                  ' atom_vector ',
     &                  SIN(QFRAMEM(3*((NATOMS/2)+J2-1)+1,J1))*COS(QFRAMEM(3*((NATOMS/2)+J2-1)+2,J1)),
     &                  SIN(QFRAMEM(3*((NATOMS/2)+J2-1)+1,J1))*SIN(QFRAMEM(3*((NATOMS/2)+J2-1)+2,J1)),
     &                  COS(QFRAMEM(3*((NATOMS/2)+J2-1)+1,J1))
            ENDDO
            ! }}}
         ELSEIF (RBAAT) THEN
            ! {{{
            WRITE(3,'(I6)') (NATOMS/2)
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEM(J1)
            DO J2=1,(NATOMS/2)
               WRITE(3,'(A1,4X,6G20.10)')
     &         'O', QFRAMEM(3*(J2-1)+1,J1), QFRAMEM(3*(J2-1)+2,J1), QFRAMEM(3*(J2-1)+3,J1),
     &         QFRAMEM(3*((NATOMS/2)+J2-1)+1,J1), QFRAMEM(3*((NATOMS/2)+J2-1)+2,J1),
     &         QFRAMEM(3*((NATOMS/2)+J2-1)+3,J1)
            ENDDO
            ! }}}
         ELSEIF (AMHT) THEN
           ! {{{
            GLY_COUNT = 0
            DO J2=1,NMRES
               IF (SEQ(J2).EQ.8) GLY_COUNT = GLY_COUNT +1
            ENDDO
            WRITE(3,'(I6)') NATOMS+GLY_COUNT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEM(J1)

            GLY_COUNT = 0
            DO J2=1,NMRES
              IF (SEQ(J2).EQ.8) THEN
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEM(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEM(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'O    ',QFRAMEM(9*(J2-1)+4-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+5-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+6-GLY_COUNT*3,J1)
                GLY_COUNT = GLY_COUNT +1
              ELSE
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEM(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'C2   ',QFRAMEM(9*(J2-1)+4-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+5-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+6-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'O    ',QFRAMEM(9*(J2-1)+7-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+8-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+9-GLY_COUNT*3,J1)
              ENDIF
          ENDDO
        ! }}}
         ELSEIF (RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')) THEN
         ! {{{
            WRITE(3,'(I6)') NOPT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEM(J1)
            WRITE(3,'(A2,4X,3G20.10)') (ZSYM(J2),QFRAMEM(J2,J1),0.0D0,0.0D0,J2=1,NOPT)
        ! }}}
         ELSE
           ! {{{
            WRITE(3,'(I6)') NATOMS
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
C           IF (J1.EQ.NFMINUS) THEN
C              WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF

            if (ZSYM(NATOMS).eq.'SV') then
               do j2=1, NATOMS
                  if (MOD(j2,3).eq.0) then
                     WRITE(3,'(A,4X,3g20.10)') 'O',QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1)
                  else
                     WRITE(3,'(A,4X,3g20.10)') 'H',QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1)
                  endif
               enddo
            else
               WRITE(3,'(A2,4X,3g20.10)') 
     &         (ZSYM(J2),QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1),J2=1,NATOMS)
            endif
            ! }}}
         ENDIF
         ! }}}
      ENDDO
      ! commented {{{
C     IF (ZSYMSAVE(1:1).EQ.'W') THEN
C        WRITE(3,'(I6)') 3*NATOMS
C        WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C        DO J1=1,NATOMS
C           CALL CONVERT(Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),
C    1                   Q(3*(NATOMS+J1-1)+1),Q(3*(NATOMS+J1-1)+2),Q(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
C           WRITE(3,'(A2,4X,3g20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
C        ENDDO
C     ELSE IF (ZSYMSAVE.EQ.'CD') THEN
C        WRITE(3,'(I6)') NATOMS*6/2
C        WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C        DO J2=1,NATOMS/2
C           CALL CAPSIDIO(Q(3*(J2-1)+1),Q(3*(J2-1)+2),Q(3*(J2-1)+3),
C    1                    Q(3*(NATOMS/2+J2-1)+1),Q(3*(NATOMS/2+J2-1)+2),Q(3*(NATOMS/2+J2-1)+3),CAPSCOORDS2,RAD,HEIGHT)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(1),CAPSCOORDS2(2),CAPSCOORDS2(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(4),CAPSCOORDS2(5),CAPSCOORDS2(6)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(7),CAPSCOORDS2(8),CAPSCOORDS2(9)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(10),CAPSCOORDS2(11),CAPSCOORDS2(12)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(13),CAPSCOORDS2(14),CAPSCOORDS2(15)
C           WRITE(3,'(A2,4X,3g20.10)') 'C4  ',CAPSCOORDS2(16),CAPSCOORDS2(17),CAPSCOORDS2(18)
C        ENDDO
C     ELSE
C        WRITE(3,'(I6)') NATOMS
C        WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C        WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J1),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),J1=1,NATOMS)
C     ENDIF
      ! }}}

      CLOSE(3)

C jmc for unres, to put in the dummy peptide groups
      IF (UNRST) THEN
        ! {{{
         WRITE(ITSTRING2,'(A)') 'unr.'//TRIM(ADJUSTL(ITSTRING))
         OPEN(UNIT=3,FILE=ITSTRING2,STATUS='UNKNOWN')
         DO J1=NFPLUS,1,-1
            DO K1=1,(NATOMS/2)-1
               DO K2=1,3
                  PEPCOORDS(6*(K1-1)+K2)=(2.0D0*QFRAMEP(6*(K1-1)+K2,J1)+QFRAMEP(6*K1+K2,J1))/3.0D0
                  PEPCOORDS(6*(K1-1)+K2+3)=(QFRAMEP(6*(K1-1)+K2,J1)+2.0D0*QFRAMEP(6*K1+K2,J1))/3.0D0
               END DO
            END DO
            WRITE(3,'(I6)') 2*NATOMS-2
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEP(J1)
C           IF (J1.EQ.NFPLUS) THEN
C              WRITE(3,'(g20.10)') EOFS(1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J2),QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1),J2=1,NATOMS)
            WRITE(3,'(A2,4X,3g20.10)') ('O ',PEPCOORDS(6*(K1-1)+1),PEPCOORDS(6*(K1-1)+2),PEPCOORDS(6*(K1-1)+3),K1=1,(NATOMS/2)-1)
            WRITE(3,'(A2,4X,3g20.10)') ('N ',PEPCOORDS(6*(K1-1)+4),PEPCOORDS(6*(K1-1)+5),PEPCOORDS(6*(K1-1)+6),K1=1,(NATOMS/2)-1)
         ENDDO
         DO K1=1,(NATOMS/2)-1
            DO K2=1,3
               PEPCOORDS(6*(K1-1)+K2)=(2.0D0*QINIT(6*(K1-1)+K2)+QINIT(6*K1+K2))/3.0D0
               PEPCOORDS(6*(K1-1)+K2+3)=(QINIT(6*(K1-1)+K2)+2.0D0*QINIT(6*K1+K2))/3.0D0
            END DO
         END DO
         WRITE(3,'(I6)') 2*NATOMS-2
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
         WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J1),QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3),J1=1,NATOMS)
         WRITE(3,'(A2,4X,3g20.10)') ('O ',PEPCOORDS(6*(K1-1)+1),PEPCOORDS(6*(K1-1)+2),PEPCOORDS(6*(K1-1)+3),K1=1,(NATOMS/2)-1)
         WRITE(3,'(A2,4X,3g20.10)') ('N ',PEPCOORDS(6*(K1-1)+4),PEPCOORDS(6*(K1-1)+5),PEPCOORDS(6*(K1-1)+6),K1=1,(NATOMS/2)-1)
         DO J1=1,NFMINUS
            DO K1=1,(NATOMS/2)-1
               DO K2=1,3
                  PEPCOORDS(6*(K1-1)+K2)=(2.0D0*QFRAMEM(6*(K1-1)+K2,J1)+QFRAMEM(6*K1+K2,J1))/3.0D0
                  PEPCOORDS(6*(K1-1)+K2+3)=(QFRAMEM(6*(K1-1)+K2,J1)+2.0D0*QFRAMEM(6*K1+K2,J1))/3.0D0
               END DO
            END DO
            WRITE(3,'(I6)') 2*NATOMS-2
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
C           IF (J1.EQ.NFMINUS) THEN
C              WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J2),QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1),J2=1,NATOMS)
            WRITE(3,'(A2,4X,3g20.10)') ('O ',PEPCOORDS(6*(K1-1)+1),PEPCOORDS(6*(K1-1)+2),PEPCOORDS(6*(K1-1)+3),K1=1,(NATOMS/2)-1)
            WRITE(3,'(A2,4X,3g20.10)') ('N ',PEPCOORDS(6*(K1-1)+4),PEPCOORDS(6*(K1-1)+5),PEPCOORDS(6*(K1-1)+6),K1=1,(NATOMS/2)-1)
         ENDDO
         CLOSE(3)
         ! }}}
      END IF ! unrst

      IF ((DUMPPATH.OR.DUMPALLPATHS).AND.(.NOT.CONNECTT)) THEN
        ! {{{
         IF (UNRST) WRITE(*,'(A)') '*** NOTE - pathlengths calculated from saved Cartesian coords will be rubbish
     & as they have been placed in the standard unres orientation.'
         IF (ZSYMSAVE.EQ.'CD') WRITE(*,'(A)') 'WARNING, symmetry and normal modes not implemented properly for CAPSID'
         IF (UNRST) THEN
           ! {{{
            IF (CALCDIHE) THEN
                CALL UNRESCALCDIHEREF(DIHE,ALLANG,QPLUS)
            ELSE
                DIHE=0.5D0 ! dummy order param for pathsample related purposes
            ENDIF
C jmc         WRITE(88,'(3g20.10)') EPLUS, DIHE, ALLANG
            WRITE(88,'(2g20.10)') EPLUS, DIHE
            ! }}}
         ELSE
            WRITE(88,'(g20.10)') EPLUS
         ENDIF
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
           ! {{{
            IF (ZSYMSAVE.EQ.'W4') IPOT=4
            IF (ZSYMSAVE.EQ.'W3') IPOT=3
            IF (ZSYMSAVE.EQ.'W2') IPOT=2
            IF (ZSYMSAVE.EQ.'W1') IPOT=1
C           DO J2=1,NATOMS
            DO J2=1,NATOMS/2 ! WCOMMENT
               CALL CONVERT(QPLUS(3*(J2-1)+1),QPLUS(3*(J2-1)+2),QPLUS(3*(J2-1)+3),
C    1                      QPLUS(3*(NATOMS+J2-1)+1),QPLUS(3*(NATOMS+J2-1)+2),QPLUS(3*(NATOMS+J2-1)+3),
     1                      QPLUS(3*(NATOMS/2+J2-1)+1),QPLUS(3*(NATOMS/2+J2-1)+2),QPLUS(3*(NATOMS/2+J2-1)+3),
     2                      OVEC,H1VEC,H2VEC)
               QW(9*(J2-1)+1)=OVEC(1) ! WCOMMENT
               QW(9*(J2-1)+2)=OVEC(2)
               QW(9*(J2-1)+3)=OVEC(3)
               QW(9*(J2-1)+4)=H1VEC(1)
               QW(9*(J2-1)+5)=H1VEC(2)
               QW(9*(J2-1)+6)=H1VEC(3)
               QW(9*(J2-1)+7)=H2VEC(1)
               QW(9*(J2-1)+8)=H2VEC(2)
               QW(9*(J2-1)+9)=H2VEC(3)
            ENDDO 
C           NATOMS=NATOMS*3 ! WCOMMENT
            NATOMSSAVE=NATOMS
            NATOMS=(NATOMS/2)*3
            CALL SYMMETRY(HORDER,.FALSE.,QW,INERTIA) ! WCOMMENT
C           NATOMS=NATOMS/3 ! WCOMMENT
            NATOMS=2*(NATOMS/3)
            NATOMS=NATOMSSAVE
            WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
C           DO J2=1,6*NATOMS ! WCOMMENT
            DO J2=1,3*NATOMS
               Q(J2)=QPLUS(J2)
            ENDDO
C           CALL H2OMODES(NATOMS,IPOT,Q,DIAG) ! WCOMMENT
            CALL H2OMODES(NATOMS/2,IPOT,Q,DIAG)
C           WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,6*NATOMS) ! WCOMMENT
            IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
            ! }}}
         ELSE IF ((FRQSPLUS(NOPT).EQ.0.0D0).OR.CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT) THEN
           ! {{{
            DO J2=1,3*NATOMS
               Q(J2)=QPLUS(J2)
            ENDDO
!           IF (.NOT.UNRST) CALL POTENTIAL(Q,EPLUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
            IF (CHRMMT) THEN
              ! {{{
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,EPLUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               IF (.NOT.NOFRQS) CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ! }}}
            ELSE IF (AMBERT) THEN
               ! {{{
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN
                   CALL MAKENUMHESS(Q,NATOMS)
               ENDIF
               IF (.NOT.NOFRQS) CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ! }}}
            ELSE IF (NABT) THEN
               ! {{{
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,EPLUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               IF (.NOT.NOFRQS) CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ! }}}
            ELSE IF (UNRST) THEN
               ! {{{
               DO J2=1,nres
                  c(1,J2)=Q(6*(J2-1)+1)
                  c(2,J2)=Q(6*(J2-1)+2)
                  c(3,J2)=Q(6*(J2-1)+3)
                  c(1,J2+nres)=Q(6*(J2-1)+4)
                  c(2,J2+nres)=Q(6*(J2-1)+5)
                  c(3,J2+nres)=Q(6*(J2-1)+6)
               ENDDO
               CALL UPDATEDC
               CALL int_from_cart(.true.,.false.)
               CALL chainbuild
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMINTHESS(NINTS,NATOMS)
                  CALL GETSTUFF(KD,NNZ,NINTB)
                  CALL INTSECDET(Q,3*NATOMS,KD,NNZ,NINTB,DIAG)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,EPLUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               ! }}}
            ELSE
               ! {{{
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,EPLUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
               IF (.NOT.NOFRQS) CALL MASSWT(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ! }}}
            ENDIF
            if (machine) then
                 WRITE(88) HORDER,FPGRP
            else
                 WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
            endif
            IF (.NOT.(UNRST.OR.NOFRQS)) THEN
               CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
               if (diag(1).lt.diag(3*natoms)) call eigensort_val_asc(diag,hess,3*natoms,3*natoms)
            ENDIF
            IF (CHRMMT.OR.AMBERT.OR.NABT) THEN
               if (machine) then
                    IF (.NOT.NOFRQS) WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
               else
                    IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
               endif
            ELSE
               IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
            ENDIF
            ! }}}
         ELSE
           ! {{{
            DO J2=1,3*NATOMS
               Q(J2)=QPLUS(J2)
            ENDDO
            CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
            WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
            WRITE(88,'(3g20.10)') (FRQSPLUS(J2),J2=1,3*NATOMS)
            ! }}}
         ENDIF

         if (machine) then
              WRITE(88) (QPLUS(J2),J2=1,NOPT)
              WRITE(88) ETS

         else
              WRITE(88,'(3F25.15)') (QPLUS(J2),J2=1,NOPT)
              WRITE(88,'(F20.10)') ETS
         endif

         IF (ZSYMSAVE(1:1).EQ.'W') THEN
           ! {{{
            IF (ZSYMSAVE.EQ.'W4') IPOT=4
            IF (ZSYMSAVE.EQ.'W3') IPOT=3
            IF (ZSYMSAVE.EQ.'W2') IPOT=2
            IF (ZSYMSAVE.EQ.'W1') IPOT=1
C           DO J2=1,NATOMS ! WCOMMENT
            DO J2=1,NATOMS/2
               CALL CONVERT(QINIT(3*(J2-1)+1),QINIT(3*(J2-1)+2),QINIT(3*(J2-1)+3),
C    1                      QINIT(3*(NATOMS+J2-1)+1),QINIT(3*(NATOMS+J2-1)+2),QINIT(3*(NATOMS+J2-1)+3),
     1                      QINIT(3*(NATOMS/2+J2-1)+1),QINIT(3*(NATOMS/2+J2-1)+2),QINIT(3*(NATOMS/2+J2-1)+3),
     2                      OVEC,H1VEC,H2VEC)
               QW(9*(J2-1)+1)=OVEC(1)
               QW(9*(J2-1)+2)=OVEC(2)
               QW(9*(J2-1)+3)=OVEC(3)
               QW(9*(J2-1)+4)=H1VEC(1)
               QW(9*(J2-1)+5)=H1VEC(2)
               QW(9*(J2-1)+6)=H1VEC(3)
               QW(9*(J2-1)+7)=H2VEC(1)
               QW(9*(J2-1)+8)=H2VEC(2)
               QW(9*(J2-1)+9)=H2VEC(3)
            ENDDO 
C           NATOMS=NATOMS*3 ! WCOMMENT
            NATOMSSAVE=NATOMS
            NATOMS=(NATOMS/2)*3
            CALL SYMMETRY(HORDER,.FALSE.,QW,INERTIA) ! WCOMMENT
C           NATOMS=NATOMS/3 ! WCOMMENT
            NATOMS=2*(NATOMS/3)
            NATOMS=NATOMSSAVE
            WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
C           DO J2=1,6*NATOMS ! WCOMMENT
            DO J2=1,3*NATOMS
               Q(J2)=QINIT(J2)
            ENDDO
C           CALL H2OMODES(NATOMS,IPOT,Q,DIAG) ! WCOMMENT
            CALL H2OMODES(NATOMS/2,IPOT,Q,DIAG)
C           WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,6*NATOMS) ! WCOMMENT
            IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
            ! }}}
         ELSE IF ((FRQSTS(NOPT).EQ.0.0D0).OR.CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT) THEN
           ! {{{
            DO J2=1,3*NATOMS
               Q(J2)=QINIT(J2)
            ENDDO
!           IF (.NOT.UNRST) CALL POTENTIAL(Q,ETS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
            IF (CHRMMT) THEN
               ! {{{
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,ETS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               IF (.NOT.NOFRQS) CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.) ! ?should this be MASSWT2 for CHARMM?
               ! }}}
            ELSE IF (AMBERT) THEN
               ! {{{
               HORDER=1
               FPGRP='C1'
               IF (.NOT.NOFRQS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
                  CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ENDIF
               ! }}}
            ELSE IF (NABT) THEN
               ! {{{
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN
                CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                CALL POTENTIAL(Q,ETS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               END IF
               IF (.NOT.NOFRQS) CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ! }}}
            ELSEIF (UNRST) THEN
               ! {{{
               DO J2=1,nres
                  c(1,J2)=Q(6*(J2-1)+1)
                  c(2,J2)=Q(6*(J2-1)+2)
                  c(3,J2)=Q(6*(J2-1)+3)
                  c(1,J2+nres)=Q(6*(J2-1)+4)
                  c(2,J2+nres)=Q(6*(J2-1)+5)
                  c(3,J2+nres)=Q(6*(J2-1)+6)
               ENDDO
               CALL UPDATEDC
               CALL int_from_cart(.true.,.false.)
               CALL chainbuild
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMINTHESS(NINTS,NATOMS)
                  CALL GETSTUFF(KD,NNZ,NINTB)
                  CALL INTSECDET(Q,3*NATOMS,KD,NNZ,NINTB,DIAG)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,ETS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               ! }}}
            ELSE
               ! {{{
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,ETS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
               IF (.NOT.NOFRQS) CALL MASSWT(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ! }}}
            ENDIF
            if (machine) then
                 WRITE(88) HORDER,FPGRP
            else
                 WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
            endif
            IF (.NOT.(UNRST.OR.NOFRQS)) THEN
               CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
               if (diag(1).lt.diag(3*natoms)) call eigensort_val_asc(diag,hess,3*natoms,3*natoms)
            ENDIF
            IF (CHRMMT.OR.AMBERT.OR.NABT) THEN
               if (machine) then
                    IF (.NOT.NOFRQS) WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
               else
                    IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
               endif
            ELSE
               IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
            ENDIF
            ! }}}
         ELSE
           ! {{{
            DO J2=1,3*NATOMS
               Q(J2)=QINIT(J2)
            ENDDO
            CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
            WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
            IF (.NOT.NOFRQS) WRITE(88,'(3g20.10)') (FRQSTS(J2),J2=1,3*NATOMS)
            ! }}}
         ENDIF
         if (machine) then
              WRITE(88) (QINIT(J2),J2=1,NOPT)
         else
              WRITE(88,'(3F25.15)') (QINIT(J2),J2=1,NOPT)
         endif

         IF (UNRST) THEN
            ! {{{
            IF (CALCDIHE) THEN
                CALL UNRESCALCDIHEREF(DIHE,ALLANG,QMINUS)
            ELSE
                DIHE=0.5D0 ! dummy order param for pathsample related purposes
            ENDIF
C jmc         WRITE(88,'(3g20.10)') EPLUS, DIHE, ALLANG
            WRITE(88,'(2g20.10)') EMINUS, DIHE
            ! }}}
         ELSE
            WRITE(88,'(g20.10)') EMINUS
         ENDIF

         IF (ZSYMSAVE(1:1).EQ.'W') THEN
            ! {{{
            IF (ZSYMSAVE.EQ.'W4') IPOT=4
            IF (ZSYMSAVE.EQ.'W3') IPOT=3
            IF (ZSYMSAVE.EQ.'W2') IPOT=2
            IF (ZSYMSAVE.EQ.'W1') IPOT=1
C           DO J2=1,NATOMS ! WCOMMENT
            DO J2=1,NATOMS/2
            ! {{{
               CALL CONVERT(QMINUS(3*(J2-1)+1),QMINUS(3*(J2-1)+2),QMINUS(3*(J2-1)+3),
C    1                   QMINUS(3*(NATOMS+J2-1)+1),QMINUS(3*(NATOMS+J2-1)+2),QMINUS(3*(NATOMS+J2-1)+3),
     1                   QMINUS(3*(NATOMS/2+J2-1)+1),QMINUS(3*(NATOMS/2+J2-1)+2),QMINUS(3*(NATOMS/2+J2-1)+3),
     2                   OVEC,H1VEC,H2VEC)
               QW(9*(J2-1)+1)=OVEC(1) ! WCOMMENT
               QW(9*(J2-1)+2)=OVEC(2)
               QW(9*(J2-1)+3)=OVEC(3)
               QW(9*(J2-1)+4)=H1VEC(1)
               QW(9*(J2-1)+5)=H1VEC(2)
               QW(9*(J2-1)+6)=H1VEC(3)
               QW(9*(J2-1)+7)=H2VEC(1)
               QW(9*(J2-1)+8)=H2VEC(2)
               QW(9*(J2-1)+9)=H2VEC(3)
               ! }}}
            ENDDO 
C           NATOMS=NATOMS*3 ! WCOMMENT
            NATOMSSAVE=NATOMS
            NATOMS=(NATOMS/2)*3
            CALL SYMMETRY(HORDER,.FALSE.,QW,INERTIA) ! WCOMMENT
C           NATOMS=NATOMS/3 ! WCOMMENT
            NATOMS=2*(NATOMS/3)
            NATOMS=NATOMSSAVE
            WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
C           DO J2=1,6*NATOMS ! WCOMMENT
            DO J2=1,3*NATOMS
               Q(J2)=QMINUS(J2)
            ENDDO
C           CALL H2OMODES(NATOMS,IPOT,Q,DIAG)
            CALL H2OMODES(NATOMS/2,IPOT,Q,DIAG)
C           WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,6*NATOMS) ! WCOMMENT
            IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
            ! }}}
         ELSE IF ((FRQSMINUS(NOPT).EQ.0.0D0).OR.CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT) THEN
           ! {{{
            DO J2=1,3*NATOMS
               Q(J2)=QMINUS(J2)
            ENDDO
!           IF (.NOT.UNRST) CALL POTENTIAL(Q,EMINUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
            IF (CHRMMT) THEN
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,EMINUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               IF (.NOT.NOFRQS) CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.) ! ?should this be MASSWT2 for CHARMM?
            ELSE IF (AMBERT) THEN
               HORDER=1
               FPGRP='C1'
               IF (.NOT.NOFRQS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
                  CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ENDIF
            ELSE IF (NABT) THEN
               ! {{{
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN              
                CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                CALL POTENTIAL(Q,EMINUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               END IF
               IF (.NOT.NOFRQS) CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ! }}}
            ELSE IF (UNRST) THEN
               ! {{{
               DO J2=1,nres
                  c(1,J2)=Q(6*(J2-1)+1)
                  c(2,J2)=Q(6*(J2-1)+2)
                  c(3,J2)=Q(6*(J2-1)+3)
                  c(1,J2+nres)=Q(6*(J2-1)+4)
                  c(2,J2+nres)=Q(6*(J2-1)+5)
                  c(3,J2+nres)=Q(6*(J2-1)+6)
               ENDDO
               CALL UPDATEDC
               CALL int_from_cart(.true.,.false.)
               CALL chainbuild
               HORDER=1
               FPGRP='C1'
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMINTHESS(NINTS,NATOMS)
                  CALL GETSTUFF(KD,NNZ,NINTB)
                  CALL INTSECDET(Q,3*NATOMS,KD,NNZ,NINTB,DIAG)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,EMINUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               ! }}}
            ELSE
               ! {{{
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                  CALL POTENTIAL(Q,EMINUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
               IF (.NOT.NOFRQS) CALL MASSWT(NATOMS,ATMASS,Q,VNEW,.TRUE.)
               ! }}}
            ENDIF
            if (machine) then
                 WRITE(88) HORDER,FPGRP
            else
                 WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
            endif
            IF (.NOT.(UNRST.OR.NOFRQS)) THEN
               CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
               if (diag(1).lt.diag(3*natoms)) call eigensort_val_asc(diag,hess,3*natoms,3*natoms)
            ENDIF
            IF (CHRMMT.OR.AMBERT.OR.NABT) THEN
               if (machine) then
                    IF (.NOT.NOFRQS) WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
               else
                    IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
               endif
            ELSE
               IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
            ENDIF
            ! }}}
         ELSE
           ! {{{
            DO J2=1,3*NATOMS
               Q(J2)=QMINUS(J2)
            ENDDO
            CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
            WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
            IF (.NOT.NOFRQS) WRITE(88,'(3g20.10)') (FRQSMINUS(J2),J2=1,3*NATOMS)
            ! }}}
         ENDIF
         if (machine) then
              WRITE(88) (QMINUS(J2),J2=1,NOPT)

         else
              WRITE(88,'(3F25.15)') (QMINUS(J2),J2=1,NOPT)
         endif
         CLOSE(88)
         ! }}}
      else if (machine.and..not.connectt) then
        ! {{{
C SAT this is for the case when we need points for minima to be output in binary format, but do not want expensive Hessian
C diagonalization, which is required to produce "path.info" file
         inquire(iolength=reclen) (diag(J1),J1=1,3*Natoms)
         open(unit=38,file="points1.out",status='unknown',form='unformatted',access='direct',recl=reclen)
         write(38,rec=1) (QPLUS(J2),J2=1,NOPT)
         close(38)
         open(unit=38,file="points2.out",status='unknown',form='unformatted',access='direct',recl=reclen)
         write(38,rec=1) (QMINUS(J2),J2=1,NOPT)
         close(38)
         ! }}}
      endif

      BFGSTST=BFGSTSTSAVE
      IVEC=IVECSAVE
      IVEC2=IVEC2SAVE

      IF (ALLOCATED(Q1)) DEALLOCATE(Q1)
      IF (ALLOCATED(Q2)) DEALLOCATE(Q2)
      IF (ALLOCATED(QW)) DEALLOCATE(QW)
      IF (ALLOCATED(QFRAMEP)) DEALLOCATE(QFRAMEP)
      IF (ALLOCATED(QFRAMEM)) DEALLOCATE(QFRAMEM)
! }}}
      ! }}}
      RETURN
      END

