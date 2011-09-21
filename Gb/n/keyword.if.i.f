      IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'.OR.WORD.EQ.'!'
     &                          .OR. WORD .EQ. '\\') THEN 
         GOTO 190
      ELSE IF ((WORD.EQ.'BASIN').OR.(WORD.EQ.'SLOPPYCONV')) THEN
         IF (NITEMS.GT.1) CALL READF(SQMAX)
      ELSE IF (WORD.EQ.'CENTRE') THEN
         CENT=.TRUE.
      ELSE IF (WORD.EQ.'CHANGEACCEPT') THEN
         CALL READI(IX)
         NACCEPT=IX
      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.
      ELSE IF (WORD.EQ.'DGUESS') THEN
         CALL READF(DGUESS)
      ELSE IF (WORD.EQ.'EDIFF') THEN
         CALL READF(EDIFF)
      ELSE IF (WORD.EQ.'G46') THEN
         G46=.TRUE.
         BLNT=.TRUE.
      ELSE IF (WORD.EQ.'MAXBFGS') THEN
         CALL READF(MAXBFGS)
      ELSE IF (WORD.EQ.'P46') THEN
         P46=.TRUE.
         BLNT=.TRUE.
      ELSE IF (WORD.EQ.'PAIRDIST') THEN
!        ! {{{
         !NPAIRS=0
         !IF (NITEMS.GT.1) THEN
            !PAIRDISTT=.TRUE.
            !WRITE(LFH,'(A)') ' keyword> Pairwise atom distances will be output to pairdists*'
            !NPAIRS=(NITEMS-1)/2
            !ALLOCATE(PAIRDIST(NPAIRS,2))
            !DO J1=1,NPAIRS
               !CALL READI(PAIRDIST(J1,1))
               !CALL READI(PAIRDIST(J1,2))
               !IF (PAIRDIST(J1,1).GT.NATOMS) THEN
                  !WRITE(LFH,'(A)') ' keyword> ERROR: PAIRDIST atom index larger than system specified!'
                  !STOP
               !ELSEIF (PAIRDIST(J1,2).GT.NATOMS) THEN
                  !WRITE(LFH,'(A)') ' keyword> ERROR: PAIRDIST atom index larger than system specified!'
                  !STOP
               !ENDIF
            !ENDDO
         !ELSE
            !YESNO=.FALSE.
            !INQUIRE(FILE='pairdist',EXIST=YESNO)
            !IF (YESNO) THEN
               !PAIRDISTT=.TRUE.
               !WRITE(LFH,'(A)') ' keyword> Pairwise atom distances will be output to pairdists*'
            !ELSE
               !WRITE(LFH,'(A)') ' keyword> ERROR: pairdist input file missing for PAIRDIST'
               !FLUSH(LFH)
               !STOP
            !ENDIF
            !OPEN(UNIT=222,FILE='pairdist',status='old')
            !DO
               !READ(222,*,IOSTAT=iostatus) CHECK1
               !IF (iostatus<0) THEN
                  !CLOSE(222)
                  !EXIT
               !ELSE 
                  !NPAIRS=NPAIRS+1
               !ENDIF
            !END DO        
            !CLOSE(222)
            !ALLOCATE(PAIRDIST(NPAIRS,2))
            !OPEN(UNIT=222,FILE='pairdist',status='old')
            !DO J1=1,NPAIRS
               !READ(222,*) PAIRDIST(J1,1),PAIRDIST(J1,2)
            !ENDDO
            !CLOSE(222)
         !ENDIF
         !WRITE(LFH,'(A)') ' keyword> Atom pair distances to monitor:'
         !DO J1=1,NPAIRS
            !WRITE(LFH,*) PAIRDIST(J1,:)
         !ENDDO
         !! }}}
      ELSE IF (WORD.EQ.'PULL') THEN
         PULLT=.TRUE.
         CALL READI(PATOM1)
         CALL READI(PATOM2)
         CALL READF(PFORCE)
         WRITE(LFH,'(A,I6,A,I6,A,G20.10)') ' keyword> Pulling atoms ',PATOM1,' and ',PATOM2,' force=',PFORCE
      ELSE IF ((WORD.EQ.'QMAX').OR.(WORD.EQ.'TIGHTCONV')) THEN
         CALL READF(FQMAX)
      ELSE IF (WORD.EQ.'RESTART') THEN
         RESTART=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NRELAX)
         IF (NITEMS.GT.2) CALL READI(NHSRESTART)
      ELSE IF (WORD.EQ.'RESTORE') THEN
         RESTORET=.TRUE.
         CALL READA(DUMPFILE)
      ELSE IF (WORD.EQ.'SAVE') THEN
         CALL READI(NSAVE)
      ELSE IF (WORD.EQ.'SEED') THEN
         SEEDT=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(IX)
            NSSTOP=IX
         ENDIF
      ELSE IF (WORD.EQ.'STEP') THEN
         NPCOUNT=NPCOUNT+1
         IF (NPCOUNT.GT.NPAR) THEN
            WRITE(LFH,'(A)') 'Number of STEP lines exceeds NPAR - quit'
            STOP
         ENDIF
         CALL READF(STEP(NPCOUNT))
         CALL READF(ASTEP(NPCOUNT))
         IF (NITEMS.GT.3) CALL READF(OSTEP(NPCOUNT))
      ELSE IF (WORD.EQ.'TRACKDATA') THEN
         TRACKDATAT=.TRUE.     
      ELSE IF (WORD.EQ.'STEPS') THEN
         NRUNS=1
         CALL READI(IX)
         MCSTEPS(1)=IX
      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         DO J1=1,NITEMS-1
            CALL READF(TEMP(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            DO J1=NITEMS,NPAR
               TEMP(J1)=TEMP(1)
            ENDDO
         ENDIF
      ELSE IF (WORD.EQ.'UPDATES') THEN
         CALL READI(MUPDATE)
      ELSE
         CALL REPORT('Unrecognized command '//WORD,.TRUE.)
         STOP
      ENDIF
