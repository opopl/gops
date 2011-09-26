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

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Produce CONNECTIONS connected minima for minimum MINDEX.
C  New transition states must be added to ts.data and new minima must be added to min.data.
C  If we find a connection to an existing minimum we should update its sum of rates. 
C
      SUBROUTINE TSSEARCH(MINDEX,NADD)
      USE PORFUNCS
      USE KEY
      USE COMMONS
      IMPLICIT NONE
      INTEGER NOFFSET, ISTAT
      INTEGER MINDEX, L1, L2, HORDERPLUS, HTS, HORDERMINUS, ECON, NADD, NINIT, J1, DMODE, J3,
     1        NDONE, PID(NCPU+NCPU+1), J2, LNTS, STATUS, J4, DUMMYI
! file numbers are offset by an additional NCPU to avoid overwriting output from connect runs.
      LOGICAL LTEST1, LTEST2, T1, T2, FINISHED(NCPU+NCPU+1), NOPATH(NCPU+NCPU+1), KILLED(NCPU+NCPU+1), LDEBUG
      DOUBLE PRECISION POINTS(3*NATOMS), POINTSPLUSLOCAL(3*NATOMS), POINTSMINUSLOCAL(3*NATOMS), POINTSTSLOCAL(3*NATOMS),
     1                 EPLUS, ETSLOCAL, EMINUS, FRQSPLUS(3*NATOMS), FRQSTS(3*NATOMS), FRQSMINUS(3*NATOMS),
     2                 IXPLUS, IYPLUS, IZPLUS, IXMINUS, IYMINUS, IZMINUS, IXM ,IYM ,IZM, DPERT, 
     3                 DUMMY, RANDOM, RANARRAY(3*NATOMS), DPRAND, LOCALPOINTS2(3*NATOMS),
     4                 DISTANCE, DIST2, RMAT(3,3), NEWIXMIN,NEWIYMIN,NEWIZMIN,FRICTIONFAC
      CHARACTER(LEN=10) J3STR, PIDSTR
      CHARACTER(LEN=80) FPOO
C     INTEGER CONNECTEDMIN(MAXCONN) ! would need to make this allocatable
      DOUBLE PRECISION TINIT, TNEW

      NOFFSET=NCPU+1
      NNEW=0
      IF (.NOT.ALLOCATED(FROZEN)) ALLOCATE(FROZEN(NATOMS))

      CALL CPU_TIME(TINIT)
C     CALL MYSYSTEM(STATUS,DEBUG,'rm points.repel')
C
C  How many connections do we already have?
C
      NTOTAL=0
C     OPEN(UNIT=47,FILE='points.repel',STATUS='UNKNOWN')
      DO L1=1,NTS
         LTEST1=PLUS(L1).EQ.MINDEX
         LTEST2=MINUS(L1).EQ.MINDEX
         IF ((LTEST1.OR.LTEST2).AND.(PLUS(L1).NE.MINUS(L1))) THEN
            NTOTAL=NTOTAL+1
C           IF (NTOTAL.GT.MAXCONN) THEN
C              WRITE(*,'(A)') 'tssearch> too many connected minima - increase MAXCONN'
C              STOP
C           ENDIF
C           IF (LTEST1) CONNECTEDMIN(NTOTAL)=MINUS(L1)
C           IF (LTEST2) CONNECTEDMIN(NTOTAL)=PLUS(L1)
C           CONNECTEDBYTS(NTOTAL)=L1
C           IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,I6)') 'minimum ',MINDEX,' is connected to minimum ',CONNECTEDMIN(NTOTAL),' by ts ',L1
            READ(UTS,REC=L1) (POINTS(L2),L2=1,3*NATOMS)
C           WRITE(47,'(3F20.10)') (POINTS(L2),L2=1,3*NATOMS)
         ENDIF
      ENDDO
C     CLOSE(47)
      IF (DEBUG) WRITE(*,'(A,I6,A,I6,A)') 'tssearch> minimum ',MINDEX,' has ',NTOTAL,' connections to different structures' 
      NINIT=NTOTAL
      IF (NTOTAL.GE.MAX(CONNECTIONS,NINIT+NADD)) RETURN
      WRITE(*,'(A,I6,A,I6,A)') 'tssearch> minimum ',MINDEX,' has ',NTOTAL,' connections to different structures - looking for more'

      READ(UMIN,REC=MINDEX) (POINTS(L2),L2=1,3*NATOMS)
      DO L1=1,MAXTSATTEMPTS
C
C  To use NCPU cpu.s first set up NCPU odata.n files and run
C  NCPU OPTIM jobs. Need to offset the indices or the path.info.n
C  files will overwrite the ones generated by connect runs.
C
         DO J3=1+NOFFSET,NCPU+NOFFSET
            WRITE(J3STR,'(I10)') J3
            IF (CHARMMT) THEN
               FPOO='odata.'//TRIM(ADJUSTL(J3STR)) ! workaround for Sun compiler bug
               OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
               DPERT=60.0D0 * PERTVALUE
C              CALL RANDOM_NUMBER(RANDOM)
               RANDOM=DPRAND()
               DMODE=NINT(RANDOM*NDIHE*2)
               DMODE=DMODE-NDIHE
               IF (DMODE.EQ.0) DMODE=1
               WRITE(1,'(A,I6,F15.5)') 'TWISTDIHE  ',DMODE,DPERT
               CLOSE(1)
               CALL MYSYSTEM(STATUS,DEBUG,'cat odata.tssearch >> odata.'//TRIM(ADJUSTL(J3STR)))
               if (machine) then
                    CALL CHARMMDUMP(POINTS,'points1.inp.'//TRIM(ADJUSTL(J3STR)))
               else
                    CALL CHARMMDUMP(POINTS,'input.crd.'//TRIM(ADJUSTL(J3STR)))
               endif
            ELSE IF (UNRST) THEN 
               FPOO='odata.'//TRIM(ADJUSTL(J3STR)) ! workaround for Sun compiler bug
               OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
C jmc do something slightly smarter to determine which dihedral to twist
C first determine whether to twist backbone or side chain dihe:
C              CALL RANDOM_NUMBER(RANDOM)
               RANDOM=DPRAND()
               IF (RANDOM.LT.0.75D0) THEN ! bias towards twisting backbone
                  DPERT=30.0D0 * PERTVALUE
C                 CALL RANDOM_NUMBER(RANDOM)
                  RANDOM=DPRAND()
                  DMODE=NINT(RANDOM*((NATOMS/2)-3)*2)
                  DMODE=DMODE-(NATOMS/2)+3
                  IF (DMODE.EQ.0) DMODE=1
               ELSE
C side chain: there are ndihe - natoms/2 + 3 of these.  Choose any with equal probability.
                  DPERT=60.0D0 * PERTVALUE
C                 CALL RANDOM_NUMBER(RANDOM)
                  RANDOM=DPRAND()
                  DMODE=(NATOMS/2)-3+NINT(RANDOM*(NDIHE-(NATOMS/2)+3))
                  IF (DMODE.EQ.(NATOMS/2-3)) DMODE=1+DMODE
                  IF (RANDOM.GT.0.5D0) DMODE=-DMODE
               ENDIF
C jmc          DPERT=60.0D0 * PERTVALUE
C C jmc          CALL RANDOM_NUMBER(RANDOM)
C jmc          RANDOM=DPRAND()
C jmc          DMODE=NINT(RANDOM*NDIHE*2) 
C jmc          DMODE=DMODE-NDIHE 
C jmc          IF (DMODE.EQ.0) DMODE=1
               WRITE(1,'(A,I6,F15.5)') 'TWISTDIHE  ',DMODE,DPERT
               CLOSE(1)
               CALL MYSYSTEM(STATUS,DEBUG,'cat odata.tssearch >> odata.'//TRIM(ADJUSTL(J3STR)))
               CALL MYUNRESDUMP(POINTS,'coords.'//TRIM(ADJUSTL(J3STR)))
            ELSE IF (AMBERT) THEN
               CALL MYSYSTEM(STATUS,DEBUG,'cp odata.tssearch odata.'//TRIM(ADJUSTL(J3STR)))
               FPOO='start.'//TRIM(ADJUSTL(J3STR)) ! workaround for Sun compiler bug
               OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
               DO L2=1,NATOMS
                  IF (FROZEN(L2)) THEN
                     RANARRAY(3*(L2-1)+1)=0.0D0
                     RANARRAY(3*(L2-1)+2)=0.0D0
                     RANARRAY(3*(L2-1)+3)=0.0D0
                  ELSE
                     RANDOM=DPRAND()
                     RANARRAY(3*(L2-1)+1)=(RANDOM-0.5D0)*2.0D0*PERTVALUE
                     RANDOM=DPRAND()
                     RANARRAY(3*(L2-1)+2)=(RANDOM-0.5D0)*2.0D0*PERTVALUE
                     RANDOM=DPRAND()
                     RANARRAY(3*(L2-1)+3)=(RANDOM-0.5D0)*2.0D0*PERTVALUE
                  ENDIF
               ENDDO
               WRITE(1,'(3F20.10)') (POINTS(3*(L2-1)+1)+RANARRAY(3*(L2-1)+1),
     &                               POINTS(3*(L2-1)+2)+RANARRAY(3*(L2-1)+2),
     &                               POINTS(3*(L2-1)+3)+RANARRAY(3*(L2-1)+3),L2=1,NATOMS)
               CLOSE(1)
            ELSE
               CALL MYSYSTEM(STATUS,DEBUG,'cp odata.tssearch odata.'//TRIM(ADJUSTL(J3STR)))
               FPOO='odata.'//TRIM(ADJUSTL(J3STR)) ! this line works around a Sun compiler bug
               OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='OLD',POSITION='APPEND')
C              WRITE(1,'(A,F20.10)') 'FIXD  ',PERTVALUE
               WRITE(1,'(A)') 'POINTS'
C              WRITE(1,'(A2,2X,3F20.10)') (ZSYMBOL(L2),POINTS(3*(L2-1)+1),POINTS(3*(L2-1)+2),POINTS(3*(L2-1)+3),L2=1,NATOMS)
               DO L2=1,NATOMS
                  IF (FROZEN(L2)) THEN
                     RANARRAY(3*(L2-1)+1)=0.0D0
                     RANARRAY(3*(L2-1)+2)=0.0D0
                     RANARRAY(3*(L2-1)+3)=0.0D0
                  ELSE
                     RANDOM=DPRAND()
                     RANARRAY(3*(L2-1)+1)=(RANDOM-0.5D0)*2.0D0*PERTVALUE
                     RANDOM=DPRAND()
                     RANARRAY(3*(L2-1)+2)=(RANDOM-0.5D0)*2.0D0*PERTVALUE
                     RANDOM=DPRAND()
                     RANARRAY(3*(L2-1)+3)=(RANDOM-0.5D0)*2.0D0*PERTVALUE
                     IF (TWOD) RANARRAY(3*(L2-1)+3)=0.0D0
                  ENDIF
               ENDDO
               WRITE(1,'(A2,2X,3F20.10)') (ZSYMBOL(L2),POINTS(3*(L2-1)+1)+RANARRAY(3*(L2-1)+1),
     &                                                 POINTS(3*(L2-1)+2)+RANARRAY(3*(L2-1)+2),
     &                                                 POINTS(3*(L2-1)+3)+RANARRAY(3*(L2-1)+3),L2=1,NATOMS)
               CLOSE(1)
            ENDIF
            CALL FLUSH(6,ISTAT) ! The child process may duplicate output without this line
            call fork_subr(PID(J3))
!
!  We have to set DEBUG to true here so that all the necessary files are copied back !
!
            LDEBUG=.TRUE.
            IF (PID(J3).EQ.0) CALL SUBMITOPTIMJOB(J3-NOFFSET,CHARMMT,UNRST,J3,EXEC,LDEBUG,'OPTIM.tssearch.')
         ENDDO
         NOPATH(1:NCPU+NCPU+1)=.FALSE.
         CALL MYWAIT(NCPU,NOFFSET,PID,NOPATH,KILLED,DEBUG) ! manage the forked processes
         WRITE(*,'(A)') 'tssearch> all forked ts searches completed or killed'
C
C  We now need the paths for any new transition states.
C  Must allow for the fact that the same new transition state could be found more than
C  once in the parallel batch of runs. Use a temporary local counter for this, LNTS.
C
         LNTS=0
         DO J2=1+NOFFSET,NCPU+NOFFSET
            FINISHED(J2)=.FALSE.
         ENDDO
         NDONE=0
         analyse_ts: DO J3=1+NOFFSET,NCPU+NOFFSET
!
! It isn't possible to use the exit status from OPTIM this way
!
            IF (KILLED(J3)) THEN
               WRITE(*,'(3(A,I5))') 'tssearch> transition state for search ',L1,' CPU ',J3-NOFFSET,' was unsuccessful'
               NOPATH(J3)=.TRUE.
               CYCLE analyse_ts
            ENDIF
            NOPATH(J3)=.FALSE.
            WRITE(J3STR,'(I10)') J3
            WRITE(PIDSTR,'(I10)') PID(J3)
C
C  Not much point calculating the path for an old transition state.
C
C           CALL MYSYSTEM(STATUS,DEBUG,'tail -1 energies.'//TRIM(ADJUSTL(PIDSTR))//' > energy')
C DAE don.t necessarily have energies file
C
            CALL MYSYSTEM(STATUS,DEBUG, 
     &            'grep -c "**** CONVERGED ****" OPTIM.tssearch.'//TRIM(ADJUSTL(PIDSTR))//' > energy')
            OPEN(1,FILE='energy',STATUS='OLD')
            READ(1,*) DUMMYI
            CLOSE(1)
            IF (DUMMYI.NE.1) THEN
               NOPATH(J3)=.TRUE.
               CYCLE analyse_ts
            END IF
            CALL MYSYSTEM(STATUS,DEBUG, 
     &            'grep Energy OPTIM.tssearch.'//TRIM(ADJUSTL(PIDSTR))//' | tail -1 | sed -e "s/.*=//g" > energy')
            OPEN(1,FILE='energy',STATUS='OLD')
            READ(1,*) ETSLOCAL
            CLOSE(1)
            IF (CHARMMT) THEN
               CALL CHARMMREAD(POINTSTSLOCAL,'points.final.'//TRIM(ADJUSTL(PIDSTR)))
            ELSE IF (UNRST) THEN
               CALL UNRESREAD(POINTSTSLOCAL,'points.final.'//TRIM(ADJUSTL(PIDSTR)))
            ELSE
               FPOO='points.final.'//TRIM(ADJUSTL(PIDSTR)) ! work around for Sun compiler bug
               OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='OLD')
               READ(1,*) (POINTSTSLOCAL(J1),J1=1,3*NATOMS)
               CLOSE(1)
            ENDIF
            CALL INERTIAWRAPPER(POINTSTSLOCAL,NATOMS,ANGLEAXIS,IXM,IYM,IZM)
            DO J1=1,NTS+LNTS
               T1=(ABS(ETSLOCAL-ETS(J1)).LT.EDIFFTOL)
               T2=.FALSE.
               IF (T1) THEN
                  READ(UTS,REC=J1) (LOCALPOINTS2(J4),J4=1,3*NATOMS)
                  CALL MINPERMDIST(POINTSTSLOCAL,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,
     &                             DIST2,RIGIDBODY,RMAT,.FALSE.)
                  IF (DISTANCE.LT.GEOMDIFFTOL) T2=.TRUE.
               ENDIF

               IF (T1.AND.T2) THEN
                  IF (J1.LE.NTS) THEN
                     WRITE(*,'(3(A,I5))') 'tssearch> transition state for search ',L1,' CPU ',J3-NOFFSET,' is old ts ',J1
                  ELSE
                     WRITE(*,'(3(A,I5))') 'tssearch> transition state for search ',L1,' CPU ',
     1               J3-NOFFSET,' already found in this batch'
                  ENDIF
                  FINISHED(J3)=.TRUE.
                  NDONE=NDONE+1
                  NOPATH(J3)=.TRUE.
                  CYCLE analyse_ts ! skip path calculation
               ENDIF
            ENDDO
C
C  Local update of ts information to avoid duplication within the same batch of parallel runs.
C  These values get set properly after the path has been calculated.
C
            LNTS=LNTS+1
111         CONTINUE
            IF (NTS+LNTS.GT.MAXTS) THEN
               CALL TSDOUBLE
               GOTO 111 ! might need to double more than once
            ENDIF
            ETS(NTS+LNTS)=ETSLOCAL
            IXTS(NTS+LNTS)=IXM
            IYTS(NTS+LNTS)=IYM
            IZTS(NTS+LNTS)=IZM
!
! this is a pain - the output files from the ts search runs have old PID's appended, so some
! of them must be renamed twice if we are going to append the PID of the path search
!
            IF (CHARMMT) THEN
            if (machine) then
               CALL MYSYSTEM(STATUS,DEBUG,'cp odata.path odata.'//TRIM(ADJUSTL(J3STR))//
     1                           '; cp points.final.'//TRIM(ADJUSTL(PIDSTR))//' points1.inp.'//TRIM(ADJUSTL(J3STR)))
            else
               CALL MYSYSTEM(STATUS,DEBUG,'cp odata.path odata.'//TRIM(ADJUSTL(J3STR))//
     1                           ' ; cp points.final.'//TRIM(ADJUSTL(PIDSTR))//' input.crd.'//TRIM(ADJUSTL(J3STR)))
            endif
            ELSE IF (UNRST) THEN
               CALL MYSYSTEM(STATUS,DEBUG,'cp odata.path odata.'//TRIM(ADJUSTL(J3STR))//
     1                           ' ; cp points.final.'//TRIM(ADJUSTL(PIDSTR))//' coords.'//TRIM(ADJUSTL(J3STR)))
            ELSE IF (AMBERT) THEN
               CALL MYSYSTEM(STATUS,DEBUG,'cp odata.path odata.'//TRIM(ADJUSTL(J3STR))//
     1                           ' ; cp points.final.'//TRIM(ADJUSTL(PIDSTR))//' start.'//TRIM(ADJUSTL(J3STR)))
            ELSE
               CALL MYSYSTEM(STATUS,DEBUG,'cp odata.path odata.'//TRIM(ADJUSTL(J3STR))//
     &                           ' ; echo "POINTS " >> odata.'//TRIM(ADJUSTL(J3STR))//
     1                           ' ; sed -e "1,/POINTS/d" odata.new.'//TRIM(ADJUSTL(PIDSTR))//' >> odata.'//TRIM(ADJUSTL(J3STR)))
            ENDIF
            CALL MYSYSTEM(STATUS,DEBUG,'cp vector.dump.'//TRIM(ADJUSTL(PIDSTR))//' vector.dump.'//TRIM(ADJUSTL(J3STR)))
            
            call fork_subr(PID(J3))
            IF (PID(J3).EQ.0) CALL SUBMITOPTIMJOB(J3-NOFFSET,CHARMMT,UNRST,J3,EXEC,DEBUG,'OPTIM.path.')
         ENDDO analyse_ts
 
         CALL MYWAIT(NCPU,NOFFSET,PID,NOPATH,KILLED,DEBUG) ! manage the forked processes
         WRITE(*,'(A)') 'tssearch> all forked pathway calculations completed or killed'

         analyse_paths: DO J3=1+NOFFSET,NCPU+NOFFSET
            IF (NOPATH(J3)) CYCLE analyse_paths
            IF (KILLED(J3)) THEN
               WRITE(*,'(3(A,I6),A,F10.3)') 'tssearch> pathway for search ',L1,' CPU ',J3,' was unsuccessful'
               CYCLE analyse_paths
            ENDIF
            WRITE(PIDSTR,'(I10)') PID(J3)
            FPOO='path.info.'//TRIM(ADJUSTL(PIDSTR)) ! work around for Sun compiler bug
            if (machine) then
                 OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='old',form='unformatted')
            else
                 OPEN(1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='OLD')
            endif
C
C  The path.info file may contain point group strings - avoid reading these.
C
            if (machine) then
                 READ(1) EPLUS
                 READ(1) HORDERPLUS
            else
                 READ(1,*) EPLUS
                 READ(1,*) HORDERPLUS
            endif
            IF (MACHINE) THEN
               IF (.NOT.NOFRQS) READ(1) (FRQSPLUS(L2),L2=1,3*NATOMS)
               READ(1) (POINTSPLUSLOCAL(L2),L2=1,3*NATOMS)
               READ(1) ETSLOCAL
               READ(1) HTS
            ELSE
               IF (.NOT.NOFRQS) READ(1,*) (FRQSPLUS(L2),L2=1,3*NATOMS)
               READ(1,*) (POINTSPLUSLOCAL(L2),L2=1,3*NATOMS), ETSLOCAL, HTS
            endif
            if (machine) then
                 IF (.NOT.NOFRQS) READ(1) (FRQSTS(L2),L2=1,3*NATOMS)
                 READ(1) (POINTSTSLOCAL(L2),L2=1,3*NATOMS)
                 READ(1) EMINUS
                 READ(1) HORDERMINUS
            else
                 IF (.NOT.NOFRQS) READ(1,*) (FRQSTS(L2),L2=1,3*NATOMS)
                 READ(1,*) (POINTSTSLOCAL(L2),L2=1,3*NATOMS),EMINUS,HORDERMINUS
            endif
            if (machine) then
                 IF (.NOT.NOFRQS) READ(1) (FRQSMINUS(L2),L2=1,3*NATOMS)
                 READ(1) (POINTSMINUSLOCAL(L2),L2=1,3*NATOMS)
            else
                 IF (.NOT.NOFRQS) READ(1,*) (FRQSMINUS(L2),L2=1,3*NATOMS)
                 READ(1,*) (POINTSMINUSLOCAL(L2),L2=1,3*NATOMS)
            endif
            CLOSE(1)
            CALL INERTIAWRAPPER(POINTSTSLOCAL,NATOMS,angleAxis,IXM,IYM,IZM) ! since the values obtained previously may have been overwritten
            IF (UNRST.AND.(.NOT.NOFRQS)) THEN
               DO L2=1,NINTS-1
                  IF (FRQSTS(L2).LT.0.0D0) THEN
                     WRITE(*,'(3(A,I5))') 'tssearch> transition state search ',L1,' CPU ',J3,
     1                     ' converged to a higher order saddle point - rejecting'
                     WRITE(*,'(A,F20.10)') 'tssearch> eigenvalue ',FRQSTS(L2)
                     CYCLE analyse_paths
                  ENDIF
               ENDDO
            ENDIF
   
            IF (DEBUG) THEN 
               WRITE(*,*) '         E+              nu+      h+         Ets             nu' //
     1                 ' ts   h ts        E-              nu-      h-'
               WRITE(*,'(3(F17.7,F15.5,I4))') EPLUS, 0.0D0, HORDERPLUS,ETSLOCAL, DUMMY, HTS, EMINUS, DUMMY, HORDERMINUS
            ENDIF

            CALL INERTIAWRAPPER(POINTSPLUSLOCAL,NATOMS,angleAxis,IXPLUS,IYPLUS,IZPLUS)
            LTEST1=ABS(EPLUS-EMIN(MINDEX)).LT.EDIFFTOL
            IF (LTEST1) THEN
               READ(UMIN,REC=MINDEX) (LOCALPOINTS2(J4),J4=1,3*NATOMS)
               CALL MINPERMDIST(POINTSPLUSLOCAL,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, 
     &                          DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
               LTEST1=.FALSE.
               IF (DISTANCE.LT.GEOMDIFFTOL) LTEST1=.TRUE.
            ENDIF

            CALL INERTIAWRAPPER(POINTSMINUSLOCAL,NATOMS,ANGLEAXIS,IXMINUS,IYMINUS,IZMINUS)
            LTEST2=ABS(EMINUS-EMIN(MINDEX)).LT.EDIFFTOL
            IF (LTEST2) THEN
               READ(UMIN,REC=MINDEX) (LOCALPOINTS2(J4),J4=1,3*NATOMS)
               CALL MINPERMDIST(POINTSMINUSLOCAL,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, 
     &                          DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
               LTEST2=.FALSE.
               IF (DISTANCE.LT.GEOMDIFFTOL) LTEST2=.TRUE.
            ENDIF 

            IF (.NOT.(LTEST1.OR.LTEST2)) THEN
               PERTVALUE=MAX(PERTMIN,PERTVALUE/1.05D0)
               WRITE(*,'(3(A,I6),A,F10.3)') 
     1              'tssearch> transition state for search ',L1,' CPU ',J3-NOFFSET,
     1              ' does not involve minimum ',MINDEX,'             PERT is now ',PERTVALUE
C              CALL MYSYSTEM(STATUS,DEBUG,'cat points.final.'//TRIM(ADJUSTL(PIDSTR))//' >> points.repel')
            ELSE IF (ABS(EPLUS-EMINUS).LT.EDIFFTOL) THEN
               WRITE(*,'(3(A,I6))') 'tssearch> transition state for search ',L1,' CPU ',J3-NOFFSET,
     1         ' corresponds to a degenerate rearrangement according to the EDIFF criterion'
C              CALL MYSYSTEM(STATUS,DEBUG,'cat points.final.'//TRIM(ADJUSTL(PIDSTR))//' >> points.repel')
            ELSE
               PERTVALUE=MIN(PERTMAX,PERTVALUE*1.05D0)
               WRITE(*,'(3(A,I6),A,F10.3)') 
     1              'tssearch> transition state for search ',L1,' CPU ',J3-NOFFSET,
     2              ' connection for minimum   ',MINDEX,' is new,     PERT is now ',PERTVALUE
CC          CALL MYSYSTEM(STATUS,DEBUG,'cat points.final >> points.repel')
               NTOTAL=NTOTAL+1
               NNEW=NNEW+1
C
C  The other minimum need not be new.
C
               DO L2=1,NMIN
                  IF (LTEST1.AND.(ABS(EMINUS-EMIN(L2)).LT.EDIFFTOL)) THEN
                     READ(UMIN,REC=L2) (LOCALPOINTS2(J4),J4=1,3*NATOMS)
                     CALL MINPERMDIST(POINTSMINUSLOCAL,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, 
     &                                   DIST2,RIGIDBODY,RMAT,.FALSE.)
                     IF (DISTANCE.LT.GEOMDIFFTOL) THEN
                        WRITE(*,'(A,I6)') 'tssearch> the connected minimum is intermediate ',L2
                        ECON=L2
C                       CONNECTEDMIN(NTOTAL)=L2
                        GOTO 11
                     ENDIF
                  ELSE IF (LTEST2.AND.(ABS(EPLUS-EMIN(L2)).LT.EDIFFTOL)) THEN
                     READ(UMIN,REC=L2) (LOCALPOINTS2(J4),J4=1,3*NATOMS)
                     CALL MINPERMDIST(POINTSPLUSLOCAL,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, 
     &                                   DIST2,RIGIDBODY,RMAT,.FALSE.)
                     IF (DISTANCE.LT.GEOMDIFFTOL) THEN
                        WRITE(*,'(A,I6)') 'tssearch> the connected minimum is intermediate ',L2
C                       CONNECTEDMIN(NTOTAL)=L2
                        ECON=L2
                        GOTO 11
                     ENDIF
                  ENDIF
               ENDDO
C
C  The other minimum is new.
C
               NMIN=NMIN+1
               IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
               WRITE(*,'(A,I6,A)') 'tssearch> new connected minimum ',NMIN,' -  writing parameters to min.data'
               ECON=NMIN
C              CONNECTEDMIN(NTOTAL)=NMIN
               GPFOLD(NMIN)=0.0D0
               IF (LTEST1) THEN
                  EMIN(NMIN)=EMINUS  
                  IF (NOFRQS) THEN
                     DUMMY=4.675754133D0 ! 2 ln(2pi) +1
                  ELSE
                     DUMMY=0.0D0
                     DO L2=NFSTART,NFFINISH
                        IF (FRQSMINUS(L2).LE.0.0D0) THEN 
                           PRINT '(A,I8,A,G15.5,A)','tssearch> error, this frequency ',J2,' is not positive: ',
     &                                           FRQSMINUS(L2),' skipping'
                        ELSE
                           DUMMY=DUMMY+LOG(FRQSMINUS(L2))
                        ENDIF
                     ENDDO
                  ENDIF
                  FVIBMIN(NMIN)=DUMMY
                  HORDERMIN(NMIN)=HORDERMINUS
                  IF (ENSEMBLE.EQ.'T') THEN
                     PFMIN(NMIN) = -EMIN(NMIN)/TEMPERATURE - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
                  ELSEIF (ENSEMBLE.EQ.'E') THEN
                     IF (TOTALE.GT.EMIN(NMIN)) THEN
                        PFMIN(NMIN) = (KAPPA-1)*LOG(TOTALE-EMIN(NMIN)) - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
                     ELSE
                        PFMIN(NMIN) = -1.0D250
                     ENDIF
                  ENDIF
                  PFMIN(NMIN)=PFMIN(NMIN)-PFMEAN

                  CALL INERTIAWRAPPER(POINTSMINUSLOCAL,NATOMS,ANGLEAXIS,NEWIXMIN,NEWIYMIN,NEWIZMIN)
                  IXMIN(NMIN)=NEWIXMIN; IYMIN(NMIN)=NEWIYMIN; IZMIN(NMIN)=NEWIZMIN
                  IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
                  WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(NMIN),FVIBMIN(NMIN),
     1                                             HORDERMIN(NMIN),IXMIN(NMIN),IYMIN(NMIN),IZMIN(NMIN)
                  WRITE(*,'(2F20.10,I6,3F20.10)') EMIN(NMIN),FVIBMIN(NMIN),
     1                                             HORDERMIN(NMIN),IXMIN(NMIN),IYMIN(NMIN),IZMIN(NMIN)
                  CALL FLUSH(UMINDATA,ISTAT)
                  IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
                  WRITE(UMIN,REC=NMIN) (POINTSMINUSLOCAL(L2),L2=1,3*NATOMS)
                  CALL FLUSH(UMIN,ISTAT)
               
               ELSE IF (LTEST2) THEN
                  EMIN(NMIN)=EPLUS
                  IF (NOFRQS) THEN
                     DUMMY=4.675754133D0 ! 2 ln(2pi) +1
                  ELSE
                     DUMMY=0.0D0
                     DO L2=NFSTART,NFFINISH
                        IF (FRQSPLUS(L2).LE.0.0D0) THEN 
                           PRINT '(A,I8,A,G15.5,A)','tssearch> error, this frequency ',J2,' is not positive: ',
     &                                           FRQSPLUS(L2),' skipping'
                        ELSE
                           DUMMY=DUMMY+LOG(FRQSPLUS(L2))
                        ENDIF
                     ENDDO
                  ENDIF
                  FVIBMIN(NMIN)=DUMMY
                  HORDERMIN(NMIN)=HORDERPLUS
                  IF (ENSEMBLE.EQ.'T') THEN
                     PFMIN(NMIN) = -EMIN(NMIN)/TEMPERATURE - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
                  ELSEIF (ENSEMBLE.EQ.'E') THEN
                     IF (TOTALE.GT.EMIN(NMIN)) THEN
                        PFMIN(NMIN) = (KAPPA-1)*LOG(TOTALE-EMIN(NMIN)) - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
                     ELSE
                        PFMIN(NMIN) = -1.0D250
                     ENDIF
                  ENDIF
                  PFMIN(NMIN)=PFMIN(NMIN)-PFMEAN

                  CALL INERTIAWRAPPER(POINTSPLUSLOCAL,NATOMS,ANGLEAXIS,NEWIXMIN,NEWIYMIN,NEWIZMIN)
                  IXMIN(NMIN)=NEWIXMIN; IYMIN(NMIN)=NEWIYMIN; IZMIN(NMIN)=NEWIZMIN
                  IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
                  WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(NMIN),FVIBMIN(NMIN),HORDERMIN(NMIN),IXMIN(NMIN),IYMIN(NMIN),
     &                                                   IZMIN(NMIN)
                  WRITE(*,'(2F20.10,I6,3F20.10)') EMIN(NMIN),FVIBMIN(NMIN),HORDERMIN(NMIN),IXMIN(NMIN),IYMIN(NMIN),
     &                                                   IZMIN(NMIN)
                  CALL FLUSH(UMINDATA,ISTAT)
                  IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
                  WRITE(UMIN,REC=NMIN) (POINTSPLUSLOCAL(L2),L2=1,3*NATOMS)
                  CALL FLUSH(UMIN,ISTAT)
               ENDIF
           
11             CONTINUE
C
C  It must be a new transition state.
C
               NTS=NTS+1
               IF (NTS.GT.MAXTS) CALL TSDOUBLE
C              CONNECTEDBYTS(NTOTAL)=NTS
               ETS(NTS)=ETSLOCAL
               HORDERTS(NTS)=HTS
               IF (NOFRQS) THEN
                  DUMMY=1.0D0
                  NEGEIG(NTS)=-1.0D0
               ELSE
                  DUMMY=0.0D0
                  DO L2=NFSTART,NFFINISH-1
                     IF (FRQSTS(J2).LT.2.0D-6) THEN 
                        WRITE(*,'(A,I6,A,G20.10)') 'getallpaths> WARNING - eigenvalue ',J2,' of this transition state is only ',
     &                                           FRQSTS(J2)
                        DUMMY=DUMMY+LOG(FRQSTS(J2))
                     ELSE
                        DUMMY=DUMMY+LOG(FRQSTS(L2))
                     ENDIF
                  ENDDO
                  NEGEIG(NTS)=FRQSTS(3*NATOMS)
               ENDIF
               FVIBTS(NTS)=DUMMY
               IF (LTEST1) THEN
                  PLUS(NTS)=MINDEX
                  MINUS(NTS)=ECON
               ELSE IF (LTEST2) THEN
                  PLUS(NTS)=ECON
                  MINUS(NTS)=MINDEX
               ENDIF
               IXTS(NTS)=IXM
               IYTS(NTS)=IYM
               IZTS(NTS)=IZM
               IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0 ! this array isn;t being allocated unless 
                                                                    ! DIJKSTRAT or KSHORTESTPATHST is true
               WRITE(*,'(A,I6,A)') 'tssearch> new ts ',NTS,' writing parameters to file ts.data'
               IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
               IF (IMFRQT) THEN
                  WRITE(UTSDATA,'(2F20.10,3I10,4F20.10)') 
     1                     ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),IXTS(NTS),IYTS(NTS),IZTS(NTS),NEGEIG(NTS)
               ELSE
                  WRITE(UTSDATA,'(2F20.10,3I10,3F20.10)') 
     1                     ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),IXTS(NTS),IYTS(NTS),IZTS(NTS)
               ENDIF
               CALL FLUSH(UTSDATA,ISTAT)
               IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
               WRITE(UTS,REC=NTS) (POINTSTSLOCAL(L2),L2=1,3*NATOMS)
               CALL FLUSH(UTS,ISTAT)
C
C  Update ts pointers.
C
               POINTERP(NTS)=-1
               POINTERM(NTS)=-1
               TOPPOINTER(PLUS(NTS))=NTS
               TOPPOINTER(MINUS(NTS))=NTS

               DO J1=NTS-1,1,-1
                  IF (PLUS(J1).EQ.PLUS(NTS)) THEN
                     POINTERP(NTS)=J1
                     GOTO 41
                  ELSE IF (MINUS(J1).EQ.PLUS(NTS)) THEN
                     POINTERP(NTS)=J1
                     GOTO 41
                  ENDIF
               ENDDO
41             CONTINUE

               DO J1=NTS-1,1,-1
                  IF (PLUS(J1).EQ.MINUS(NTS)) THEN
                     POINTERM(NTS)=J1
                     GOTO 42
                  ELSE IF (MINUS(J1).EQ.MINUS(NTS)) THEN
                     POINTERM(NTS)=J1
                     GOTO 42
                  ENDIF
               ENDDO
42             CONTINUE
C
C  ts pointers have been updated.
C
C  Must calculate rates. Don;t bother with KSUM, as this is calculated now when needed.
C
               IF (ENSEMBLE.EQ.'T') THEN
                  KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(PLUS(NTS))  /
     1               (2.0D0 * PI*HORDERTS(NTS))) +
     2             (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS)) / 2.0D0
     3           - (ETS(NTS) - EMIN(PLUS(NTS)))/TEMPERATURE
                  IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))

                  KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(MINUS(NTS)) /
     1                (2.0D0 * PI*HORDERTS(NTS))) +
     2             (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS)) / 2.0D0
     3           - (ETS(NTS) - EMIN(MINUS(NTS)))/TEMPERATURE
                  IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
               ELSE
                  IF (TEMPERATURE.GT.ETS(NTS)) THEN
                     KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(PLUS(NTS))  / (2*PI*HORDERTS(NTS))) +
     1                      (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS))/2 + 
     2                      (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(PLUS(NTS))))
                     KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2*PI*HORDERTS(NTS))) +
     1                      (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS))/2 + 
     2                      (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(MINUS(NTS))))
                  ELSE
                     KPLUS(NTS)=-1.0D250
                     KMINUS(NTS)=-1.0D250
                  ENDIF
               ENDIF

               IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
               IF (ZSYM(1:2).EQ.'CA') KMINUS(NTS)=KMINUS(NTS)+30.66356D0
               IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)
               IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)

!              PRINT '(3I8,4G20.10)',NTS,PLUS(NTS),MINUS(NTS),KPLUS(NTS),KMINUS(NTS)
C
C  Update sum of rates out of the connected minima.
C
!              IF (KSUM(PLUS(NTS)).EQ.0.0D0) THEN
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(PLUS(NTS))=KPLUS(NTS)
!              ELSE
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(PLUS(NTS)) =LOG(EXP(KSUM(PLUS(NTS))-KMEAN) + EXP(KPLUS(NTS) -KMEAN)) + KMEAN
!              ENDIF
!              IF (KSUM(MINUS(NTS)).EQ.0.0D0) THEN
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(MINUS(NTS))=KMINUS(NTS)
!              ELSE
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(MINUS(NTS))=LOG(EXP(KSUM(MINUS(NTS))-KMEAN) + EXP(KMINUS(NTS)-KMEAN)) + KMEAN
!              ENDIF
            ENDIF
            IF (ADDPT) CALL ADDPERM(POINTSTSLOCAL,POINTSPLUSLOCAL,POINTSMINUSLOCAL)
         ENDDO analyse_paths
!         IF (NTOTAL.GE.MAX(CONNECTIONS,NINIT+NADD)) GOTO 20
         IF (NTOTAL.GE.MAX(CONNECTIONS,NINIT+NADD)) CALL CPU_TIME(TNEW) 

         IF (NTOTAL.GE.MAX(CONNECTIONS,NINIT+NADD)) GOTO 20
      ENDDO
 
!20    CALL CPU_TIME(TNEW)
20      TTSSEARCH=TTSSEARCH+TNEW-TINIT

      RETURN
      END
