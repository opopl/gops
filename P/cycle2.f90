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

SUBROUTINE CYCLE2
USE COMMONS
USE PORFUNCS
IMPLICIT NONE
DOUBLE PRECISION SPOINTS(3*NATOMS), FPOINTS(3*NATOMS), DISTANCE
INTEGER J1, J2, J3, PID(NCPU), MINS(NCPU), MINF(NCPU), STATUS, NCYCLES, NMINOLD, NTSOLD, ISTAT, NEWJOB, PIDDONE, NREJ, NINITIAL
INTEGER NAVAIL, NUSED, CHILDPID
LOGICAL KILLED, STOPCYCLE, CHANGEPD
CHARACTER(LEN=10) CONNSTR

NAVAIL=0
NMINOLD=NMIN; NTSOLD=NTS
IF ((NPFOLD.GT.0).AND.(.NOT.DIJPAIRT)) CALL PFOLD
OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
CLOSE(1)

IF (DEBUG) PRINT '(A)','cycle2> removing previous OPTIM files'
! CALL MYSYSTEM(STATUS,DEBUG,'rm -f *[0-9] >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f odata.[0-9] >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f path.[0-9]*.xyz.[0-9]* >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f EofS.[0-9]* >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f OPTIM.connect.[0-9]* >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f path.info.[0-9]* >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f min.data.info.[0-9]* >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f input.crd.[0-9] >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f finish.[0-9] >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f points1.inp.[0-9] >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f points2.inp.[0-9] >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f vector.dump.[0-9] >& /dev/null ')
CALL MYSYSTEM(STATUS,DEBUG,'rm -f coords.[0-9] >& /dev/null ')

! NUSED is otherwise uninitialised
NUSED=0 
IF ((NMIN.EQ.2).AND.(DIJINITT.OR.DIJINITFLYT)) THEN ! no point running more than one initial search in this case
   CALL GETDPAIR(NAVAIL,NUSED,MINS(1),MINF(1),SPOINTS,FPOINTS)
   CALL CONNECTODATA(1,SPOINTS,FPOINTS)
   CALL FLUSH(6,ISTAT) ! the child process may duplicate output without this line!
   CALL FORK_SUBR(PID(1)) ! PID is zero in the child, non-zero in the parent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  IF (PID(1).NE.0) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(1)
!  PRINT *,'PID(1)=',PID(1)
!  IF (PID(1).EQ.0) THEN
!     WRITE(*,'(A,I8)') 'cycle2> I am the child! PID=',PID(1)
!     CALL GETPID_SUBR(CHILDPID)
!     PRINT *,'CHILDPID=',CHILDPID
!  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (PID(1).EQ.0) CALL SUBMITOPTIMJOB(1,CHARMMT,UNRST,1,EXEC,DEBUG,'OPTIM.connect.')
   CALL WAIT_SUBR(PIDDONE,STATUS)
!
! Failure to return the right pid here is not fatal if we have the path.info file.
! It shouldn;t happen, though! Calling wait a second time seems to allow the system
! to catch up!
!
   IF (PIDDONE.NE.PID(1)) THEN
      PRINT '(2(A,I10))','cycle2> ERROR - initial path WAIT returned process id',PIDDONE,' instead of ',PID(1)
      CALL WAIT_SUBR(PIDDONE,STATUS)
      PRINT '(2(A,I10))','cycle2> on calling wait again pid=',PIDDONE,' status=',STATUS
      IF (PIDDONE.NE.PID(1)) STOP
   ENDIF
   WRITE(CONNSTR,'(I10)') PID(1)
   WRITE(*,'(A,I8)') 'cycle2> analysing result of initial search for process id ',PID(1)
   IF (BHINTERPT.OR.BISECTT) THEN
      CALL MYSYSTEM(STATUS,DEBUG,'cp min.data.info.'//TRIM(ADJUSTL(CONNSTR))//' min.data.info')
   ELSE
      CALL MYSYSTEM(STATUS,DEBUG,'cp path.info.'//TRIM(ADJUSTL(CONNSTR))//' path.info')
   ENDIF
   IF (STATUS.EQ.0) THEN ! The file exists, so try to analyse it
      IF (BHINTERPT.OR.BISECTT) THEN
         CALL GETALLMIN(MINS(1),MINF(1))
      ELSE ! IF (TRIPLES) THEN
         CALL GETALLPATHS
!     ELSE
!        CALL GETNEWPATH(MINS(1),MINF(1))
      ENDIF
   ELSE
      PRINT '(A)','cycle2> ERROR - no path.info file generated by initial search'
      STOP
   ENDIF
ENDIF

DO J3=1,NCPU
   CALL SYSTEM('sleep 1') ! to prevent us running out of source ports. Needed for M13.
   IF (DIJPAIRT) THEN
      CALL GETDPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (CONNECTREGIONT) THEN
      CALL GETRPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
!  ELSEIF (PTAUT) THEN
!     CALL GETPPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (SHORTCUTT) THEN
      CALL GETSPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (UNTRAPT) THEN
      CALL GETUPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (FREEPAIRT) THEN
      CALL GETFREEPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSEIF (USEPAIRST) THEN
      CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in GETUSEPAIR
      CALL GETUSEPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ELSE
      CALL GETPAIR(NAVAIL,NUSED,MINS(J3),MINF(J3),SPOINTS,FPOINTS)
   ENDIF
   CALL CONNECTODATA(J3,SPOINTS,FPOINTS)
   CALL FLUSH(6,ISTAT) ! the child process may duplicate output without this line!
   CALL FORK_SUBR(PID(J3)) ! PID is zero in the child, non-zero in the parent
   IF (DEBUG.AND.(PID(J3).NE.0)) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(J3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  IF (PID(J3).NE.0) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(J3)
!  IF (PID(J3).EQ.0) THEN
!     WRITE(*,'(A,I8)') 'cycle2> I am the child! PID=',PID(J3)
!     CALL GETPID_SUBR(CHILDPID)
!     PRINT *,'CHILDPID=',CHILDPID
!  ENDIF
!  PRINT '(A,2I8)','cycle2> J3,PID=',J3,PID(J3)
   IF (PID(J3).EQ.0) CALL SUBMITOPTIMJOB(J3,CHARMMT,UNRST,J3,EXEC,DEBUG,'OPTIM.connect.')
ENDDO

cycles: DO NCYCLES=1,NATTEMPT*NCPU ! the last NCPU steps do not submit new jobs

   KILLED=.FALSE.
   CALL FLUSH(6,ISTAT) ! the child process may duplicate output without this line!
   CALL FLUSH(6,ISTAT)
   CALL WAIT_SUBR(PIDDONE,STATUS)
11 CONTINUE
!  PRINT '(A,5I8)','cycle2> PIDDONE,STATUS,PID=',PIDDONE,STATUS,PID(1:NCPU)
   CALL FLUSH(6,ISTAT) ! the child process may duplicate output without this line!
   IF (PIDDONE.GT.0) THEN
      IF (DEBUG) PRINT '(A,I8,A,I6)','cycle2> PID ',PIDDONE,' has finished with exit status ',STATUS
      DO J2=1,NCPU
         IF (PIDDONE.EQ.PID(J2)) THEN
            IF (STATUS.NE.0) KILLED=.TRUE. ! INCOMPLETE OPTIM JOBS WOULD IDEALLY RETURN A NON-ZERO EXIT CODE 
            NEWJOB=J2
            GOTO 10
         ENDIF
      ENDDO
      PRINT*,'ERROR - PID of completed child process not recognised: ',PIDDONE
      STOP
   ELSE
      CALL FLUSH(6,ISTAT) ! the child process may duplicate output without this line!
      PRINT '(A,I20)','cycle2> WARNING - WAIT returned system error code ',-PIDDONE
!
! Try calling wait again to see if this fixes things. 
! For very short OPTIM jobs WAIT may have trouble keeping up!
!
      CALL MYSYSTEM(STATUS,DEBUG,' sleep 1')
      CALL WAIT_SUBR(PIDDONE,STATUS)
      PRINT '(2(A,I8))','cycle2> on calling wait again pid=',PIDDONE,' status=',STATUS
      IF (PIDDONE.GT.0) GOTO 11
      STOP
   ENDIF
10 WRITE(*,'(3(A,I8))') 'cycle2> forked connect run ',NCYCLES,' on CPU ',NEWJOB,' completed or killed process id ',PID(NEWJOB)
!
!  It is important to identify OPTIM jobs that did not terminate with exit code 0 
!  In such cases KILLED should be .TRUE.
!
   CALL FLUSH(6,ISTAT) ! the child process may duplicate output without this line!
   WRITE(*,'(3(A,I8))') 'cycle2> analysing result of search ',NCYCLES,' on CPU ',NEWJOB,' for process id ',PID(NEWJOB)
   WRITE(CONNSTR,'(I10)') PID(NEWJOB)
   IF (KILLED) WRITE(*,'(3(A,I8))') 'cycle2> connection ',NCYCLES,' on CPU ',NEWJOB,' was unsuccessful'
!
!  If KILLED is .TRUE. there could be a viable path.info file if DUMPALLPATHS is set in
!  OPTIM. 
!  It would be nice if we had the OPTIM exit status when running in a distributed
!  environment - but we don;t! Instead we have the exit status of the attempt to copy
!  back the path.info file for distributed environments, and the OPTIM exit status
!  for SMP. We can go ahead and try to analyse a path.info file so long as it exists!
!
   IF (BHINTERPT.OR.BISECTT) THEN
      CALL MYSYSTEM(STATUS,DEBUG,'cp min.data.info.'//TRIM(ADJUSTL(CONNSTR))//' min.data.info')
   ELSE
      CALL MYSYSTEM(STATUS,DEBUG,'cp path.info.'//TRIM(ADJUSTL(CONNSTR))//' path.info')
   ENDIF
   IF (STATUS.EQ.0) THEN ! the file exists, so try to analyse it!
      IF (BHINTERPT.OR.BISECTT) THEN
         CALL GETALLMIN(MINS(NEWJOB),MINF(NEWJOB))
      ELSE ! IF (TRIPLES) THEN
         CALL GETALLPATHS
!     ELSE
!        CALL GETNEWPATH(MINS(NEWJOB),MINF(NEWJOB))
      ENDIF
   ENDIF

   IF (NCYCLES.LE.NATTEMPT*NCPU-NCPU) THEN ! submit replacement job
      IF ((.NOT.DEBUG).AND.(.NOT.COPYOPTIMT)) &
   &  CALL MYSYSTEM(STATUS,DEBUG,'rm -f *' // TRIM(ADJUSTL(CONNSTR)) // ' >& /dev/null' ) ! remove old output
      CALL SYSTEM('sleep 1') ! to prevent us running out of source ports. Needed for M13.
      IF (DIJPAIRT) THEN
         CALL GETDPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (CONNECTREGIONT) THEN
         CALL GETRPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
!     ELSEIF (PTAUT) THEN
!        CALL GETSPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (SHORTCUTT) THEN
         CALL GETSPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (UNTRAPT) THEN
         CALL GETUPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (FREEPAIRT) THEN
         CALL GETFREEPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSEIF (USEPAIRST) THEN
         CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in GETUSEPAIR
         CALL GETUSEPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ELSE
         CALL GETPAIR(NAVAIL,NUSED,MINS(NEWJOB),MINF(NEWJOB),SPOINTS,FPOINTS)
      ENDIF
!     WRITE(*,'(2(A,I6),A,F12.1,A,F12.3,A,I8)') 'cycle2> connecting minima ',MINS(NEWJOB),' and ', &
! &      MINF(NEWJOB), ' distance=',DISTANCE,' |Pfold diff|=',ABS(GPFOLD(MINS(NEWJOB))-GPFOLD(MINF(NEWJOB))),' rejects=',NREJ
      CALL CONNECTODATA(NEWJOB,SPOINTS,FPOINTS)
      CALL FLUSH(6,ISTAT) ! the child process may duplicate output without this line!
      CALL FORK_SUBR(PID(NEWJOB))
!     IF (DEBUG.AND.(PID(NEWJOB).NE.0)) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(NEWJOB)
!     IF (PID(NEWJOB).NE.0) WRITE(*,'(A,I8)') 'cycle2> forked connect run process id=',PID(NEWJOB)
!     IF (PID(NEWJOB).EQ.0) WRITE(*,'(A,I8)') 'cycle2> I am the child! PID=',PID(NEWJOB)
      IF (PID(NEWJOB).EQ.0) CALL SUBMITOPTIMJOB(NEWJOB,CHARMMT,UNRST,NEWJOB,EXEC,DEBUG,'OPTIM.connect.')
   ENDIF

   IF (MOD(NCYCLES,NCPU).EQ.0) THEN 
      WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                    '--------------------------------------------------'
      WRITE(*,'(5(A,I8))') 'cycle2> end of cycle ',NCYCLES/NCPU,' new min=',NMIN-NMINOLD,' new ts=',NTS-NTSOLD, &
  &                          ' total min=',NMIN,' total ts=',NTS
      WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                    '--------------------------------------------------'
      NMINOLD=NMIN; NTSOLD=NTS
      IF ((NPAIRDONE.GT.0).AND.(.NOT.DUMMYRUNT)) THEN
         OPEN(UNIT=1,FILE='pairs.data',STATUS='UNKNOWN')
         WRITE(1,'(2I8)') (PAIR1(J1),PAIR2(J1),J1=1,NPAIRDONE)
         CLOSE(1)
      ENDIF
      IF (NMINDONE.GT.0) THEN
         OPEN(UNIT=1,FILE='min.done',STATUS='UNKNOWN')
         WRITE(1,'(I8)') (MINDONE(J1),J1=1,NMINDONE)
         CLOSE(1)
      ENDIF
      IF (PFOLDINT.NE.0) THEN
         IF (MOD(NCYCLES,PFOLDINT*NCPU).EQ.0) THEN
            IF ((NPFOLD.GT.0).AND.(.NOT.DIJPAIRT)) CALL PFOLD
            OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
            WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
            CLOSE(1)
         ENDIF
      ENDIF
   ENDIF

   INQUIRE(FILE='STOP',EXIST=STOPCYCLE)
   IF (STOPCYCLE) THEN
      PRINT '(A)','cycle2> File STOP detected - exit'
      EXIT
   ENDIF
   INQUIRE(FILE='pathdata.change',EXIST=CHANGEPD)
   IF (CHANGEPD) THEN
      PRINT '(A)','cycle2> rereading parameter file'
      CALL MYSYSTEM(STATUS,DEBUG,'mv pathdata pathdata.orig')
      CALL MYSYSTEM(STATUS,DEBUG,'mv pathdata.change pathdata')
      CALL KEYWORD
   ENDIF

ENDDO CYCLES

RETURN

END SUBROUTINE CYCLE2
