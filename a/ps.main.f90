      ! info {{{
!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

C  New version of PATHSAMPLE using the committor probability formulation
C  and a new sampling strategy for stationary points. Paths themselves
C  are not generally referenced directly.
C  Initially we need to know
C  which minima are in which funnel, i.e. members of the sets A and B for
C  which the rate constant is required, and there may be pre-existing databases
C  of minima and transition states containing energies and products of non-zero
C  eigenvalues. We need the ability to start from a single
C  path.info file.
C
C  We are allowing connections between any A and B permutational isomers.
C  Use TAG to distinguish permutational isomers.
C
C SAT: BOUNDARYT variable is referenced but never set! Assume default value of .False. and ignore relevant bits of code.
C SAT: SAVEALL variable is set to .True. in setup if 'paths' file exists and 'pathpointer' does not!
C
      ! }}}
      PROGRAM PATHSAMPLE
      ! declarations {{{
      USE COMMONS
      USE NODES, ONLY: SSHPARALLEL
      USE DOCKMODULE
      USE RIGIDBODYMOD, ONLY: INITIALISERIGIDBODY, CLEANRIGIDBODIES
      IMPLICIT NONE
      INTEGER J1, NSIZE, NWORST, NAVAIL, NMINSAVE
      INTEGER, ALLOCATABLE :: MINMAP(:)
!     INTEGER, ALLOCATABLE :: NSEED(:)
      LOGICAL MASSFILE
      DOUBLE PRECISION TINIT, TNEW, DUMMY1(1), DUMMY2(1)

!     AMH LOCAL VARIABLES
      INTEGER :: NRES,I_RES,J2
      CHARACTER(LEN=5) :: TARFL
      CHARACTER(LEN=2) :: SDUMMY
      INTEGER :: SEQ(500)
! }}}
      ! body {{{
!
!  CONNECTIONS is the minimum number of connected minima.
!  TEMPERATURE is the temperature in reduced units.
!  Use EMICRO for a microcanonical total energy. Don't try to access inaccessible
!  stationary points above EMICRO!
!
      CALL CPU_TIME(TINIT)
      CALL SDPRND(ISEED)

      PRINTT=.TRUE.
      TTSSEARCH=0.0D0
      TPFOLD=0.0D0
      TTFOLD=0.0D0
      TDIJKSTRA=0.0D0
      TKSHORTESTPATHS=0.0D0
      TCONNECTDIST=0.0D0
      TGT=0.0D0

      CALL KEYWORD(3)

! memory allocations
      include ps.main.am.inc.f90
! set MASS, ZSYMBOL      
      include ps.main.mass.inc.f90
! set NFSTART, NFFINISH, KAPPA
include ps.main.nfstart.inc.f90

include ps.main.write.inc.f90

      IF (REGROUPT.AND.(NATTEMPT.GT.0)) THEN
         PRINT '(A)','ERROR - regroup can only be called once, set CYCLES=0'
!        STOP
      ENDIF
      IF (DOCKT) CALL DOCK
      IF (SIST) THEN
         CALL SETUP_SIS
      ELSE
         CALL SETUP
      END IF
      CALL CPU_TIME(TNEW)
      WRITE(*,'(A,F15.5)') 'main> time spent in setup=',TNEW-TINIT

      IF (NATTEMPT.GT.0) THEN
         ! {{{
         IF (CONNECTIONS.GT.2) THEN
            CALL CYCLE
         ELSE
            CALL CYCLE2
         ENDIF
         OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
         WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
         CLOSE(1)
         IF (NPAIRDONE.GT.0) THEN
            OPEN(UNIT=1,FILE='pairs.data',STATUS='UNKNOWN')
            WRITE(1,'(2I8)') (PAIR1(J1),PAIR2(J1),J1=1,NPAIRDONE)
            CLOSE(1)
         ENDIF
         IF (NMINDONE.GT.0) THEN
            OPEN(UNIT=1,FILE='min.done',STATUS='UNKNOWN')
            WRITE(1,'(I8)') (MINDONE(J1),J1=1,NMINDONE)
            CLOSE(1)
         ENDIF
         ! }}}
      ELSE IF (TFOLDT) THEN          ! just run TFOLD subroutine
         CALL TFOLD
      ELSE IF (ARNOLDIT) THEN           ! just run diagonalise2 subroutine
         CALL REWEIGHT
      ELSE IF (DIAGT) THEN           ! just run diagonalise2 subroutine
         CALL DIAGONALISE2
      ELSE IF (NGTT) THEN            ! just run NGT subroutine
         CALL NGT
      ELSE IF (GTT.OR.GT2T) THEN     ! just run GT subroutine
         CALL GT
      ELSE IF (NPFOLD.GT.0) THEN     ! just run PFOLD subroutine
         CALL PFOLD
         OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
         WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
         CLOSE(1)
      ELSE IF (KMCT) THEN            ! just run KMC 
        ! {{{
!
!  REGROUP and REGROUPFREE change NMINA, NMINB, LOCATIONA, LOCATIONB, so save the values and reset
!  to call GT more than once.
!  In fact, REGROUPFREE changes NMIN, NTS, etc. so we must stop after such a run or
!  do a complete reset somehow. This is done explicitly in routines like getppair
!  using the SAVESTATE module. Not needed here because we assume that NGT cannot be
!  called more than once!
!
         IF (REGROUPFREET.OR.REGROUPFREEABT) THEN
            CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
            CALL REGROUPFREE2(.FALSE.,1,NAVAIL)
            TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                                 ! before calling GETNCONN
            MAXBARRIER=HUGE(1.0D0)
         ENDIF
   
         IF (REGROUPRATET.OR.REGROUPPET) THEN
            CALL REGROUPFREE
            TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                                 ! before calling GETNCONN
            MAXBARRIER=HUGE(1.0D0)
         ENDIF
!
!  Must not call GETNCONN here or MAXCONN may be reset to a lower value, which messes things up.
!
         CALL GETNCONN
!
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!
         ALLOCATE(MINMAP(NMIN))
         CALL REGROUP(MINMAP)
         CALL KMC                    ! standard stochastic KMC
         ! }}}
      ELSE IF (DIJKSTRAT .OR. KSHORTESTPATHST) THEN ! run shortest path analysis
         NMINSAVE=NMIN
         IF (REGROUPFREET.OR.REGROUPFREEABT) THEN
            CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
            CALL REGROUPFREE2(.FALSE.,1,NAVAIL)
            TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                                 ! before calling GETNCONN
            MAXBARRIER=HUGE(1.0D0)
            CALL GETNCONN
!  
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!  
            ALLOCATE(MINMAP(NMIN))
            CALL REGROUP(MINMAP)
         ELSEIF (REGROUPRATET.OR.REGROUPPET) THEN
            CALL REGROUPFREE
            TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                                 ! before calling GETNCONN
            MAXBARRIER=HUGE(1.0D0)
            CALL GETNCONN
!  
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!  REGROUP also reorders the stationary points so that the most connected minima appear
!  at the end in order to accelerate GT.
!  
            ALLOCATE(MINMAP(NMIN))
            CALL REGROUP(MINMAP)
         ELSE
            CALL GETNCONN
         ENDIF
         IF (.NOT.ALLOCATED(MINMAP)) THEN ! allow for the case of no regrouping
            ALLOCATE(MINMAP(NMIN))
            DO J1=1,NMIN
               MINMAP(J1)=J1
            ENDDO
         ENDIF
         IF (DIJKSTRAT) THEN
            CALL DIJKSTRA(NWORST,.FALSE.,0,NMINSAVE,MINMAP)
         ELSEIF (KSHORTESTPATHST) THEN
            CALL KSHORTESTPATHS(NWORST,.FALSE.,0,NMINSAVE,MINMAP)
         ENDIF
         IF (ALLOCATED(MINMAP))  DEALLOCATE(MINMAP)
         IF (ALLOCATED(MINGROUP))  DEALLOCATE(MINGROUP)
      ELSE IF (REGROUPFREET.OR.REGROUPFREEABT) THEN    ! just run new REGROUPFREE2 subroutine
         CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
         CALL REGROUPFREE2(.FALSE.,1,NAVAIL)
         DEALLOCATE(MINGROUP)
      ELSE IF (REGROUPRATET) THEN    ! just run new REGROUPFREE subroutine
         CALL REGROUPFREE
         DEALLOCATE(MINGROUP)
      ELSE IF (REGROUPPET) THEN      ! just run new REGROUPFREE subroutine
         CALL REGROUPFREE
         DEALLOCATE(MINGROUP)
      ELSE IF (KMCCOMMITT) THEN      ! just run new KMC assuming equilibrium in B region
         CALL GETNCONN
         CALL KMCCOMMIT
      ENDIF
      IF (DIJPAIRT) THEN
         OPEN(UNIT=1,FILE='ts.attempts',STATUS='UNKNOWN')
         WRITE(1,'(I8)') TSATTEMPT(1:NTS)
         CLOSE(1)
      ENDIF

! memory deallocations
      include ps.main.deam.inc.f90

      include ps.main.write.cputime.inc.f90

      ! }}}
      STOP
      END
!
!  Subroutine GETNCONN sets up array NCONN containing the number of 
!  connections for each minimum after pruning according to the value of
!  NCONNMIN.
!
      SUBROUTINE GETNCONN
      ! declarations {{{
      USE COMMON, ONLY: NMIN,NTS,NCONN,PLUS,MINUS,NCONNMIN,NCONNMAX,MAXMIN,DEBUG,ETS, 
     &                  KPLUS, KMINUS, EMIN
      IMPLICIT NONE
      INTEGER J1, PNCONNECTED, NCONNECTED, NZERO, JMAX
      LOGICAL, ALLOCATABLE :: CONNECTED(:) ! reallocate MAXMIN when used
      LOGICAL :: DEADTS
      DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0
      ! }}}
      ! subroutine body {{{
!
!  Record the number of connections for each minimum in NCONN.
!  NCONN is the number of connections to minima with more
!  than NCONNMIN connections.
!
      NCONNECTED=0
      IF (ALLOCATED(NCONN)) DEALLOCATE(NCONN)
      ALLOCATE(CONNECTED(MAXMIN),NCONN(MAXMIN))
      DO J1=1,NMIN
         CONNECTED(J1)=.TRUE.
      ENDDO
11    DO J1=1,NMIN
         NCONN(J1)=0
      ENDDO
      PNCONNECTED=NCONNECTED
      DO J1=1,NTS
! JMC n.b. don't apply the nconnmin criteria at this point, hence the huge(1) 's in place of NCONN() for the plus and minus minima.        
         CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),HUGE(1),HUGE(1), 
     &                PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,DEADTS)
         IF ((.NOT.DEADTS).AND.(PLUS(J1).NE.MINUS(J1))) THEN
            IF (CONNECTED(MINUS(J1))) NCONN(PLUS(J1))=NCONN(PLUS(J1))+1
            IF (CONNECTED(PLUS(J1)))  NCONN(MINUS(J1))=NCONN(MINUS(J1))+1
         ENDIF
      ENDDO
      NCONNECTED=0
      DO J1=1,NMIN
         CONNECTED(J1)=.FALSE.
         IF (NCONN(J1).GT.NCONNMIN) THEN
            CONNECTED(J1)=.TRUE.
            NCONNECTED=NCONNECTED+1
         ENDIF
      ENDDO
      IF (NCONNECTED.NE.PNCONNECTED) GOTO 11

!     DO J1=1,NMIN
!        IF (DEBUG) WRITE(*,'(A,I6,A,I6)') 'getnconn> number of connections for minimum ',J1,' is ',NCONN(J1)
!     ENDDO

      NCONNMAX=NCONN(1)
      NZERO=0
      IF (NCONN(1).EQ.0) NZERO=1
      JMAX=1
      DO J1=2,NMIN
         IF (NCONN(J1).EQ.0) NZERO=NZERO+1
         IF (NCONN(J1).GT.NCONNMAX) THEN
            NCONNMAX=NCONN(J1)
            JMAX=J1
         ENDIF
      ENDDO
      WRITE(*,'(4(A,I6))') 'getnconn> max connections: ',NCONNMAX,' for min ',JMAX,' # of zeros=',NZERO,
     &                     ' after removing minima with < ',NCONNMIN+1
      DEALLOCATE(CONNECTED)

      RETURN
      END
