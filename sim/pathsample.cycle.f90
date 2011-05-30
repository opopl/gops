
SUBROUTINE CYCLE
! declarations {{{
USE COMMONS
USE PORFUNCS

IMPLICIT NONE

DOUBLE PRECISION SPOINTS(3*NATOMS), FPOINTS(3*NATOMS), RANDOM, DISTANCE, DPRAND, RMAT(3,3), DIST2
INTEGER ICPU, PID(NCPU), MINS(NCPU), MINF(NCPU), STATUS, NCYCLES, NMINOLD, NTSOLD, ISTAT, NREJ, J1
LOGICAL KILLED(NCPU), NOJOB(NCPU), STOPCYCLE
CHARACTER(LEN=10) CONNSTR
! }}}

NMINOLD=NMIN; NTSOLD=NTS
CALL PFOLD
OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
CLOSE(1)

cycles: DO NCYCLES=1,NATTEMPT
! {{{
   IF (DEBUG) PRINT '(A)','cycle> removing previous OPTIM files'
!  CALL MYSYSTEM(STATUS,DEBUG, &
!&   'rm odata.[0-9]* *.xyz.[0-9]* finish.[0-9]* EofS.*[0-9]* OPTIM.*.[0-9]* energies.[0-9]* points.[0-9]*') 
   CALL MYSYSTEM(STATUS,DEBUG,'rm *[0-9]')
!
!  Create NCPU start/finish pairs, which require NCPU OPTIM connect runs.
!  Make the NCONNRUNS odata.connect files, then run them, then process the results.
!
   DO ICPU=1,NCPU
   ! {{{
      NREJ=-1
! 20    CALL RANDOM_NUMBER(RANDOM)
20    RANDOM=DPRAND()
      NREJ=NREJ+1
      MINS(ICPU)=NINT(0.5D0+NMIN*RANDOM)     ! starting min - change variable name
! 10    CALL RANDOM_NUMBER(RANDOM)
10    RANDOM=DPRAND()
      MINF(ICPU)=NINT(0.5D0+NMIN*RANDOM)     ! finishing min - change variable name
      IF (MINF(ICPU).EQ.MINS(ICPU)) GOTO 10
!
!  Do the check based on GPFOLD difference first so that we save the MIND
!  call if it fails.
!
!  Any difference > PSCALE will be accepted.
!
!     CALL RANDOM_NUMBER(RANDOM)
      RANDOM=DPRAND()
      IF (0.0D0+ABS(GPFOLD(MINS(ICPU))-GPFOLD(MINF(ICPU)))/PSCALE.LT.RANDOM) THEN
         IF (NREJ.LT.1000000) THEN
!           IF (DEBUG) PRINT '(A,2I6,A,G20.10)','rejected pair ',MINS(ICPU),MINF(ICPU),' difference=', &
! &                     ABS(GPFOLD(MINS(ICPU))-GPFOLD(MINF(ICPU)))/PSCALE
            GOTO 20
         ENDIF
      ENDIF
      READ(UMIN,REC=MINS(ICPU)) SPOINTS(1:3*NATOMS)
      READ(UMIN,REC=MINF(ICPU)) FPOINTS(1:3*NATOMS)
      CALL MINPERMDIST(SPOINTS,FPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
!     CALL RANDOM_NUMBER(RANDOM)
      RANDOM=DPRAND()
!
!  Calculate metric for connection attempt probability. We probably want this
!  to depend upon the minimised distance and the difference in committor probabilities.
!
!  Any DISTANCE < DSCALE is accepted; then the probability decreases exponentially.
!
      IF (EXP(-(DISTANCE-DSCALE)/DSCALE).LT.RANDOM) THEN
         IF (NREJ.LT.1000000) THEN
!           IF (DEBUG) PRINT '(A,2I6,A,G20.10,A,G20.10)','rejected pair ',MINS(ICPU),MINF(ICPU),' distance=', &
! &                     DISTANCE,' exponent=',EXP(-(DISTANCE-DSCALE)/DSCALE)
            GOTO 20
         ENDIF
      ENDIF

      CALL CONNECTODATA(ICPU,SPOINTS,FPOINTS)
      WRITE(*,'(2(A,I6),A,F12.1,A,F12.3,A,I8)') 'cycle> connecting between minima ',MINS(ICPU),' and ', &
  &      MINF(ICPU), ' distance=',DISTANCE,' |Pfold diff|=',ABS(GPFOLD(MINS(ICPU))-GPFOLD(MINF(ICPU))),' rejects=',NREJ
      NOJOB(ICPU)=.FALSE.
!
!  Don;t bother to check whether we already have a single ts that does this connection.
!
      CALL FLUSH(6,ISTAT) ! the child process may duplicate output without this line!
      call fork_subr(PID(ICPU))
!     IF (DEBUG.AND.(PID(ICPU).NE.0)) WRITE(*,'(A,I8)') 'forked connect run process id=',PID(ICPU)
!     IF (PID(ICPU).NE.0) WRITE(*,'(A,I8)') 'forked connect run process id=',PID(ICPU)
!     IF (PID(ICPU).EQ.0) WRITE(*,'(A,I8)') 'I am the child! PID=',PID(ICPU)
      IF (PID(ICPU).EQ.0) CALL SUBMITOPTIMJOB(ICPU,CHARMMT,UNRST,ICPU,EXEC,DEBUG,'OPTIM.connect.')
      ! }}}
   ENDDO

   CALL MYWAIT(NCPU,0,PID,NOJOB,KILLED,DEBUG) ! manage the forked processes
   WRITE(*,'(A)') 'cycle> all forked connect runs are completed or killed'

!
!  It is important to identify OPTIM jobs that did not terminate with exit code 0.
!  In such cases KILLED should be .TRUE.
!
   analyse_connect: DO ICPU=1,NCPU
   ! {{{
      WRITE(*,'(A)') ' '
      WRITE(*,'(A,I6,A,I8)') 'cycle> analysing result of search ',ICPU,' for process id ',PID(ICPU)
      IF (KILLED(ICPU)) THEN
         WRITE(*,'(A,I6,A)') 'cycle> connection ',ICPU,' was unsuccessful'
!
!  Nevertheless, there could be a viable path.info file if DUMPALLPATHS is set in
!  OPTIM and TRIPLES in pathsample. Give it a try?
!
         IF (TRIPLES) THEN
            PRINT '(A)','cycle> attempting to analyse the path.info file nevertheless'
         ELSE
            CYCLE analyse_connect
         ENDIF
      ENDIF
      WRITE(CONNSTR,'(I10)') PID(ICPU)
!     CALL MYSYSTEM(STATUS,DEBUG,'cp OPTIM.connect.'//TRIM(ADJUSTL(CONNSTR))//' OPTIM.connect') ! not needed
      CALL MYSYSTEM(STATUS,DEBUG,'cp path.info.'//TRIM(ADJUSTL(CONNSTR))//' path.info')
      IF (TRIPLES) THEN
         CALL GETALLPATHS
      ELSE
         CALL GETNEWPATH(MINS(ICPU),MINF(ICPU))
      ENDIF
      ! }}}
   ENDDO analyse_connect

   WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                 '--------------------------------------------------'
   WRITE(*,'(5(A,I8))') 'cycle> end of cycle ',NCYCLES,' new min=',NMIN-NMINOLD,' new ts=',NTS-NTSOLD, &
  &                       ' total min=',NMIN,' total ts=',NTS
   WRITE(*,'(A)')   '------------------------------------------------------------' // &
  &                 '--------------------------------------------------'
   NMINOLD=NMIN; NTSOLD=NTS

   IF (NPAIRDONE.GT.0) THEN
      OPEN(UNIT=1,FILE='pairs.data',STATUS='UNKNOWN')
      WRITE(1,'(2I8)') (PAIR1(J1),PAIR2(J1),J1=1,NPAIRDONE)
      CLOSE(1)
   ENDIF
   IF (NMINDONE.GT.0) THEN
      OPEN(UNIT=1,FILE='min.done',STATUS='UNKNOWN')
      WRITE(1,'(2I8)') (MINDONE(J1),J1=1,NMINDONE)
      CLOSE(1)
   ENDIF
   IF (PFOLDINT.NE.0) THEN
      IF (MOD(NCYCLES,PFOLDINT).EQ.0) THEN
         CALL PFOLD
         OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
         WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
         CLOSE(1)
      ENDIF
   ENDIF
   INQUIRE(FILE='STOP',EXIST=STOPCYCLE)
   IF (STOPCYCLE) THEN
      PRINT '(A)','File STOP detected - exit'
      EXIT
   ENDIF
! }}}
ENDDO cycles

RETURN

END SUBROUTINE CYCLE
