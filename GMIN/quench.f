!op226>=================================== 
!op226> GPL License Info {{{ 
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
!op226>}}} 
!op226>=================================== 
C  Conjugate gradient driver. 
C  CFLAG convergence test
C  CTEST checks for changes in chirality for AMBER runs
C
!> \param QTEST LOGICAL 
!> \param NP INTEGER
!> \param ITER INTEGER
!> \param TIME DOUBLE PRECISION
!> \param BRUN INTEGER
!> \param QDONE INTEGER
!> \param P DOUBLE PRECISION(3*NATOMS)
      SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,BRUN,QDONE,P)
!op226> Declarations {{{ 
      use COMMONS
      USE QMODULE
      use porfuncs
      IMPLICIT NONE

      INTEGER I, J1, NSQSTEPS, NP, IFLAG, ITER, NOPT, J2, NDUMMY, J3, CSMIT
      DOUBLE PRECISION P(3*NATOMS),POTEL,TIME,EREAL,RBCOORDS(18),TMPCOORDS(3*NATOMS), DIST, QE, QX, AVVAL, CSMRMS
      LOGICAL QTEST, CFLAG, RES, COMPON, evapreject
      DOUBLE PRECISION  GRAD(3*NATOMS), DUMMY, DUM(3*NATOMS), DISTMIN, SSAVE, DIST2, RMAT(3,3)
C     DOUBLE PRECISION  WORK(60*NATOMS)
      DOUBLE PRECISION, PARAMETER :: HALFPI=1.570796327D0

      CHARACTER(LEN=80) DSTRING
      COMMON /MYPOT/ POTEL
      COMMON /CO/ COMPON
      COMMON /DMIN/ DISTMIN
      LOGICAL GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      common /ev/ evapreject
      DOUBLE PRECISION QSTART, QFINISH
      COMMON /Q4C/ QSTART, QFINISH
      COMMON /CSMAVVAL/ AVVAL, CSMRMS, CSMIT

C
C   sf344> gradually changing parameters to prevent dissociation of PY ellipsoids with repulsive sites 
C
      DOUBLE PRECISION epssave(3)

C
C  Data for the screen saver.
C
      INTEGER BRUN, QDONE,ii
!op226>}}} 
C
C  Turn on guiding potentials. These get turned off in potential.F when
C  the RMS force is small enough.
C
      SSAVE=STEP(NP)
C
C csw34 Reset the NFIX counter
C
      NFIX=0

11    IF (WELCH) TOSI=.TRUE.
      IF (PACHECO) AXTELL=.FALSE.
      IF (CPMD) SCT=.TRUE.
      IF (ZETT1.OR.ZETT2) THEN
         MORSET=.TRUE.
         RHO=6.0D0
      ENDIF
      IF (NATBT.AND.GUPTAT) GUIDET=.TRUE.
      IF (NATBT.AND.GUIDET) GUPTAT=.TRUE.
      IF (CSMGUIDET) CSMDOGUIDET=.TRUE.
      NOPT=3*NATOMS
      IF (WENZEL) NOPT=2
      IF (MULLERBROWNT) NOPT=2
C
C  QTEST is set for the final quenches with tighter convergence criteria.
C
      IF (QTEST) THEN
         GMAX=CQMAX
      ELSE
         GMAX=BQMAX
      ENDIF

      QDONE=0
      DO I=1,3*NATOMS
         P(I)=COORDS(I,NP)
      ENDDO
C
C     IF (TIP) THEN
C        WRITE(DUMPXYZUNIT+NP,'(I6)') (NATOMS/2)*3
C        WRITE(DUMPXYZUNIT+NP,70) NP,NQ(NP), EREAL, RMS
C        DO J2=1,NATOMS/2
C           CALL TIPIO(P(3*(J2-1)+1),P(3*(J2-1)+2),P(3*(J2-1)+3),
C    1           P(3*(NATOMS/2+J2-1)+1),P(3*(NATOMS/2+J2-1)+2),P(3*(NATOMS/2+J2-1)+3),RBCOORDS)
C           WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
C           WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
C           WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
C        ENDDO
C     ENDIF


      IF (COMPRESST.AND.(.NOT.QTEST)) THEN
         COMPON=.TRUE.
         IF (PATCHY) THEN
           CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,1.D1*GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         ELSE
           CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         END IF
         POTEL=EREAL
         IF (.NOT.CFLAG) WRITE(MYUNIT,'(A,I7,A)') ' WARNING - compressed quench ',NQ(NP),'  did not converge'
         WRITE(MYUNIT,'(A,I7,A,F20.10,A,I5,A,F15.7,A,I4,A,F12.2)') 'Comp Q ',NQ(NP),' energy=',
     1              POTEL,' steps=',ITER,' RMS force=',RMS
      ENDIF
      COMPON=.FALSE.


10    IF (PERMOPT) THEN ! lb415
         !IF ( NQ(NP) .eq. 1) THEN
         IF (DUMPT) THEN
            WRITE(DUMPXYZUNIT+NP,'(I6)') NDUMMY
            WRITE(DUMPXYZUNIT+NP,'(A,I6)') 'quench> initial points before quench ',NQ(NP)
            WRITE(DUMPXYZUNIT+NP,'(A,3G20.10)') ('LA ',P(3*(J2-1)+1),P(3*(J2-1)+2),P(3*(J2-1)+3),J2=1,NATOMS)
         ENDIF
         CALL POTENTIAL(P,GRAD,EREAL,.FALSE.,.FALSE.)
         CFLAG=.TRUE.
!        ITER=1
!        RMS=0.0D0
         RMS=CSMRMS
         ITER=CSMIT
         !ELSE
         !   CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP) ! minimize structure
         !   write(*,*) 'permdist mylbfgs', EREAL, ITER, RMS
         !   POTEL=EREAL
         !   IF (.NOT.CFLAG) WRITE(MYUNIT,'(A,I7,A)') 'WARNING - Quench ',NQ(NP),'  did not converge'
         !   DO II=1,NSAVE
         !      IF ( II .GE. NQ(NP) ) EXIT ! There's no need to check further, there's nothing
         !      CALL MINPERMDIST(P,QMINP(II,:),NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY,DIST2,RIGID,RMAT)
         !      write(*,*) DUMMY, 'dummy',ii
         !      IF (DUMMY .LT. 0.5D0) THEN
         !         !DO NOT ACCEPT THIS QUENCH
         !         WRITE(MYUNIT,*) 'This quench ended in a known minimum. It won`t be counted.'
         !         RETURN
         !      ENDIF
         !   ENDDO
         !ENDIF 
      ELSEIF (MODEL1T) THEN
         CALL MODEL1(P,GRAD,EREAL,QE,QX)
         EREAL=QE
         CFLAG=.TRUE.
         ITER=1
         RMS=0.0D0
         P(1)=QX
      ELSE IF (DL_POLY) THEN
C
C  Need to make DL_POLY input file for current coordinates.
C
         OPEN (UNIT=91,FILE='CONFIG',STATUS='OLD')
         OPEN (UNIT=92,FILE='config',STATUS='UNKNOWN')
         READ(91,'(A80)') DSTRING
         WRITE(92,'(A80)') DSTRING
         READ(91,'(A80)') DSTRING
         WRITE(92,'(A80)') DSTRING
         DO J1=1,NATOMS
            READ(91,'(A80)') DSTRING
            WRITE(92,'(A80)') DSTRING
            READ(91,'(A80)') DSTRING
            WRITE(92,'(3G20.10)') P(3*(J1-1)+1),P(3*(J1-1)+2),P(3*(J1-1)+3)
            READ(91,'(A80)') DSTRING
            WRITE(92,'(A80)') DSTRING
            READ(91,'(A80)') DSTRING
            WRITE(92,'(A80)') DSTRING
         ENDDO
         CLOSE(91)
         CLOSE(92)
         CALL SYSTEM('cp CONFIG CONFIG.old; cp config CONFIG')
         CALL SYSTEM('DLPOLY.X > output.DL_POLY ; tail -9 STATIS > energy')
         OPEN (UNIT=91,FILE='energy',STATUS='OLD')
         READ(91,*) EREAL
         WRITE(MYUNIT,'(A,G20.10)') 'energy=',EREAL
         CLOSE(91)
         OPEN(UNIT=91,FILE='REVCON',STATUS='OLD')
         READ(91,'(A1)') DUMMY
         READ(91,'(A1)') DUMMY
         NATOMS=0
13       READ(91,'(A1)',END=14) DUMMY
         NATOMS=NATOMS+1
         READ(91,*) P(3*(NATOMS-1)+1),P(3*(NATOMS-1)+2),P(3*(NATOMS-1)+3)
         READ(91,'(A1)') DUMMY
         READ(91,'(A1)') DUMMY
C        WRITE(MYUNIT,'(3G20.10)') P(3*(NATOMS-1)+1),P(3*(NATOMS-1)+2),P(3*(NATOMS-1)+3)
         GOTO 13
14       CONTINUE
         CLOSE(91)
         CFLAG=.TRUE.
C
C  Read the coordinates of the minimised geometry into vector P.
C
C     ELSE IF (BFGS .AND.(.NOT.QTEST)) THEN
      ELSE IF (BFGS) THEN
C        CALL CGMIN(100,P,CFLAG,ITER,EREAL,NP)
         CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,100,ITER,.TRUE.,NP)
         CALL DFPMIN(MAXIT,P,3*NATOMS,GMAX,ITER,EREAL,CFLAG)
      ELSEIF (TNT) THEN
C        CALL CGMIN(100,P,CFLAG,ITER,EREAL,NP)
         CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,100,ITER,.TRUE.,NP)
          WRITE(MYUNIT, '(A)') 'subroutine tn does not compile with NAG/PG'
         STOP
C        CALL TN(IFLAG,3*NATOMS,P,EREAL,GRAD,WORK,60*NATOMS,GMAX,ITER,MAXIT,CFLAG,DEBUG)
      ELSEIF (CONJG) THEN
         CALL CGMIN(MAXIT,P,CFLAG,ITER,EREAL,NP)
    ! 
! Compute quantum energy with Variation Gaussian Wavepacket.
! Coords are scaled by VGW LJ sigma (LJSIGMA) inputed with VGW params.
! Coords are then scaled back to unit sigma.
! 
      ELSEIF (VGW) THEN    
        IF(QTEST) THEN              
          CALL VGWQUENCH(P,EREAL,CFLAG)
          ELSE
            CALL VGWQUENCHSP(P,EREAL,CFLAG)
        ENDIF 

      ELSEIF (MYSDT) THEN
         CALL MYSD(MAXIT,P,CFLAG,ITER,EREAL)
      ELSEIF (RKMIN) THEN
         CALL ODESD(MAXIT,P,CFLAG,ITER,EREAL,NP)
      ELSEIF (BSMIN) THEN
         CALL ODESD(MAXIT,P,CFLAG,ITER,EREAL,NP)
      ELSE
C        CALL CGMIN(5,P,CFLAG,ITER,EREAL,NP)
         IF (CHRMMT.AND.INTMINT) THEN
            CALL MYLBFGS(NINTS,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         ELSE IF (THOMSONT) THEN
            TMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,NP)
            CALL THOMSONCTOANG(TMPCOORDS,P,NATOMS)
            CALL MYLBFGS(2*NATOMS,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
            CALL THOMSONANGTOC(P,NATOMS)

!         ELSE IF(PYBINARYT) THEN
!! sf344> trying out some sort of systematic parameter change to prevent particles from dissociating:
!! first decrease repulsive epsilon values, converge, then gradually increase them
!           epssave(:)=PEPSILON1(:)
!          IF(.NOT.QTEST) THEN
!           WRITE(MYUNIT,*) 'first iteration: decreasing epsilon_rep values by a factor of 10000' 
!           PEPSILON1(:)=PEPSILON1(:)/10000.0D0
!            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
!           WRITE(MYUNIT,*) 'second iteration: increasing epsilon_rep values by a factor of 100' 
!           PEPSILON1(:)=PEPSILON1(:)*100.0D0
!            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
!            WRITE(MYUNIT,*) 'third iteration: increasing epsilon_rep values by a factor of 100' 
!           PEPSILON1(:)=PEPSILON1(:)*100.0D0
!            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
!          END IF
         ELSE
            CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
         ENDIF
         IF (EVAPREJECT) RETURN
      ENDIF
      POTEL=EREAL

      IF (CFLAG) QDONE=1
      IF (.NOT.CFLAG) THEN
         IF (QTEST) THEN
            WRITE(MYUNIT,'(A,I6,A)') 'WARNING - Final Quench ',NQ(NP),'  did not converge'
         ELSE
            IF (NPAR.GT.1) THEN
               WRITE(MYUNIT,'(A,I7,A)') 'WARNING - Quench ',NQ(NP),'  did not converge'
            ELSE
               WRITE(MYUNIT,'(A,I7,A)') 'WARNING - Quench ',NQ(NP),'  did not converge'
            ENDIF
         ENDIF
      ENDIF

      CALL MYCPU_TIME(TIME)

      RES=.FALSE.
      IF (TABOOT.AND.(.NOT.QTEST).AND.(.NOT.RENORM)) THEN
         CALL TABOO(EREAL,POTEL,P,NP,RES)
         IF (RES) GOTO 10
      ENDIF

C     PRINT*,'Taboo lists:'
C     DO J1=1,NPAR
C        PRINT*,'Parallel run ',J1
C        WRITE(*,'(6F15.7)') (ESAVE(J2,J1),J2=1,NT(J1))
C     ENDDO
C     PRINT*,'Inertia lists:'
C     DO J1=1,NPAR
C        PRINT*,'Parallel run ',J1
C        WRITE(MYUNIT,'(6F15.7)') (XINSAVE(J2,J1),J2=1,NT(J1))
C     ENDDO

      IF (SAVEQ) CALL GSAVEIT(EREAL,P,NP)
!     IF (QDONE.EQ.0) THEN
!        PRINT '(A)','WARNING quench did not converge from starting coodinates:'
!        WRITE(MYUNIT,'(3G20.10)') (COORDS(J1,NP),J1=1,3*NATOMS)
!     ENDIF
C
C  If EPSSPHERE is non-zero we are presumably doing a calculation of the 
C  energy density of local minima. We need to know the minimum distance
C  between the starting point and the quenched minima.
C
      IF ((EPSSPHERE.NE.0.0D0).OR.BSWL) THEN
         DO J1=1,3*NATOMS
            DUM(J1)=COORDS(J1,NP)
         ENDDO
C
C  DUM is returned in the closest orientation to P; P should not change.
C  This is nearly the same mind as OPTIM! To execute a random walk we must take 
C  another step and minimise until the distance between the starting point
C  and the quench minimum is less than EPSSPHERE.
C
!        CALL MINDGMIN(P,DUM,NATOMS,DISTMIN,PERIODIC,TWOD)
         CALL NEWMINDIST(P,DUM,NATOMS,DISTMIN,PERIODIC,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
      ENDIF
C
C  Deal with EPSSPHERE sampling.
C
      IF (EPSSPHERE.NE.0.0D0) THEN
         IF ((DISTMIN.GT.EPSSPHERE).OR.(ABS(EREAL-EPREV(NP)).LE.ECONV)) THEN
            WRITE(MYUNIT,'(A,F12.5,A,4F14.5)') 'step ',STEP(NP),' EREAL, EPREV, DISTMIN, EPSSPHERE=',
     1                                     EREAL, EPREV(NP), DISTMIN, EPSSPHERE
            DO J1=1,3*NATOMS
               COORDS(J1,NP)=COORDSO(J1,NP)
            ENDDO
            CALL TAKESTEP(NP)
             WRITE(MYUNIT,'(A,G20.10)' ) 'reseeding step, maximum displacement reset to ',STEP(NP)
            GOTO 11
         ELSE
            WRITE(MYUNIT,'(A,2F20.10)') 'valid step, DISTMIN, EPSSPHERE=',DISTMIN, EPSSPHERE
         ENDIF
      ENDIF
C
C  If we are provided with target minimum coordinates in file coords.target then
C  calculate the minimum distances. May be useful for algorithm development.
C  If we get close, we don;t want to escape without a hit!
C
!     IF (ALLOCATED(TCOORDS)) THEN
!        DO J1=1,NTARGETS
!           TMPCOORDS(1:3*NATOMS)=TCOORDS(J1,1:3*NATOMS)
!           CALL MINPERMDIST(P,TMPCOORDS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY,DIST2,RIGID)
!           WRITE(MYUNIT, '(A,I5,A,F15.3,A,F15.3,A,F20.10)') 'for target structure ',J1,' dist=',DUMMY,' dist2=',DIST2,' V=',POTEL
!        ENDDO
!        DO J1=1,MIN(NMSBSAVE,MAXSAVE)
!           TMPCOORDS(1:3*NATOMS)=MSBCOORDS(1:3*NATOMS,J1)
!           CALL MINPERMDIST(P,TMPCOORDS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY,DIST2,RIGID)
!           PRINT '(A,I5,A,F15.3,A,F15.3,A,F20.10)','for taboo  structure ',J1,' dist=',DUMMY,' dist2=',DIST2,' V=',POTEL
!        ENDDO
!     ENDIF
C
C  NORESET true does not set the configuration point to the quench geometry
C  A relaxed frozen core does not get saved, but the lowest minima are saved
C  by GSAVEIT.
C
      IF (.NOT.NORESET) THEN
         DO J1=1,3*(NATOMS-NSEED)
            COORDS(J1,NP)=P(J1)
         ENDDO
         DO J1=1,NATOMS
            VAT(J1,NP)=VT(J1)
         ENDDO
      ENDIF

      IF (Q4T) CALL ORDERQ4(NATOMS,P,QFINISH)
C
C  Calling CENTRE here without an evaporation check can put particles
C  outside the container, and make a valid step in takestep impossible.
C
C     PRINT*,'Calling centre from quench'
C     IF ((.NOT.FIELDT).AND.(.NOT.SEEDT).AND.CENT) CALL CENTRE2(COORDS(1:3*NATOMS,NP))

      IF (DUMPT) THEN
         IF (ARNO) THEN
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS+2
            WRITE(DUMPXYZUNIT+NP,70) NP,NQ(NP),EREAL,RMS
            WRITE(DUMPXYZUNIT+NP,'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(DUMPXYZUNIT+NP,'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            IF (NS.NE.0) WRITE(DUMPXYZUNIT+NP,65) (P(I),I=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
         ELSE IF (TIP) THEN
            WRITE(DUMPVUNIT-NP,'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT+NP,'(I6)') (NATOMS/2)*3
            WRITE(DUMPXYZUNIT+NP,70) NP,NQ(NP), EREAL, RMS
            DO J2=1,NATOMS/2
               CALL TIPIO(P(3*(J2-1)+1),P(3*(J2-1)+2),P(3*(J2-1)+3),
     1              P(3*(NATOMS/2+J2-1)+1),P(3*(NATOMS/2+J2-1)+2),P(3*(NATOMS/2+J2-1)+3),RBCOORDS)
               WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
               WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
               WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            ENDDO
         ELSE IF (CHRMMT) THEN
            CALL CHARMMDUMP3(P)
            CALL CHARMMDUMP2(P,DUMPXYZUNIT+NP) ! xyz
         ELSEIF (NCORE(NP).GT.0) THEN
            WRITE(DUMPVUNIT-NP,'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,70) NQ(NP), EREAL, RMS
C           WRITE(DUMPXYZUNIT+NP,80) ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NCORE(NP))
C           WRITE(DUMPXYZUNIT+NP,80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NCORE(NP)+1,NATOMS)
            WRITE(DUMPXYZUNIT+NP,80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NCORE(NP))
            IF (NCORE(NP).GT.0) WRITE(DUMPXYZUNIT+NP,80) 
     &                     ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NCORE(NP)+1,NATOMS)
         ELSE
            WRITE(DUMPVUNIT-NP,'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,70) NQ(NP), EREAL, RMS
70          FORMAT(1X,'QUENCH NUMBER ',I6,' final energy=',F20.10,' RMS force=',E20.10)
            WRITE(DUMPXYZUNIT+NP,80) ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NS)
            IF (NS.NE.0) WRITE(DUMPXYZUNIT+NP,80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NS+1,NATOMS)
80          FORMAT(A2,3F20.10)
         ENDIF
      ENDIF

      IF (SQUEEZET) THEN
         IF ((EREAL.GT.0.0D0).AND.(SQUEEZED.LT.1.0D0)) THEN
            SQUEEZED=2.0D0-SQUEEZED
            NSQSTEPS=NQ(NP)
         ELSE
            NSQSTEPS=100000
         ENDIF
         DO J1=1,3*NVEC
            VEC(J1)=VEC(J1)*SQUEEZED
         ENDDO
         IF (NQ(NP).GT.2*NSQSTEPS) SQUEEZET=.FALSE.
      ENDIF
    
      IF ((NQ(NP).GE.NSSTOP).AND.SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
         WRITE(MYUNIT,'(I6,A,G20.10)') NSSTOP,' quenches completed, setting coordinates to the lowest minimum, E=',QMIN(1)
         DO J1=1,3*NATOMS
            COORDS(J1,NP)=QMINP(1,J1)
         ENDDO
         POTEL=QMIN(1)
         EREAL=POTEL
      ENDIF

      RETURN
      END
