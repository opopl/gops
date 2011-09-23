
      IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'.OR.WORD.EQ.'!'&
     &                          .OR. WORD .EQ. '\\') THEN 
         GOTO 190

         ! P46 G46 MYBLN BLN BLNGO {{{
      ELSE IF (WORD.EQ.'P46') THEN
         P46=.TRUE.
         BLNT=.TRUE.
      ELSE IF (WORD.EQ.'G46') THEN
         G46=.TRUE.
         BLNT=.TRUE.
      ELSE IF (WORD.EQ.'MYBLN') THEN
         MYBLNT=.TRUE.
         BLNT=.TRUE.
      ELSE IF (WORD.EQ.'BLN') THEN
!        ! BLN {{{
         !BLNT=.TRUE.
         !CALL READF(RK_R)
         !CALL READF(RK_THETA)
         !ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     !&            LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS))
         !OPEN(UNIT=1,FILE='BLNsequence',STATUS='OLD')
         !READ(1,*) DUMMYCH
         !READ(1,*) LJREPBB, LJATTBB
         !READ(1,*) LJREPLL, LJATTLL
         !READ(1,*) LJREPNN, LJATTNN
         !READ(1,*) DUMMYCH
         !READ(1,*) DUMMYCH
         !READ(1,*) HABLN, HBBLN, HCBLN, HDBLN
         !READ(1,*) EABLN, EBBLN, ECBLN, EDBLN
         !READ(1,*) TABLN, TBBLN, TCBLN, TDBLN
         !DO J1=1,NATOMS-1
            !READ(1,'(A1)',ADVANCE='NO') BEADLETTER(J1)
         !ENDDO
         !READ(1,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
         !DO J1=1,NATOMS-3
            !READ(1,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
         !ENDDO
         !CLOSE(1)
         !WRITE(LFH,'(A,I8,A)') 'BLN sequence of ',NATOMS,' beads read:'
         !WRITE(LFH,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
         !WRITE(LFH,'(A)') ' '
         !WRITE(LFH,'(A,I8,A)') 'BLN dihedral types:'
         !WRITE(LFH,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
         !WRITE(LFH,'(A)') ' '
         !WRITE(LFH,'(A,2F15.5)') 'B-B LJ coefficients: ',LJREPBB, LJATTBB
         !WRITE(LFH,'(A,2F15.5)') 'L-L LJ coefficients: ',LJREPLL, LJATTLL
         !WRITE(LFH,'(A,2F15.5)') 'N-N LJ coefficients: ',LJREPNN, LJATTNN
         !WRITE(LFH,'(A,4F15.5)') 'Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
         !WRITE(LFH,'(A,4F15.5)') 'Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
         !WRITE(LFH,'(A,4F15.5)') 'Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
         !call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     !&                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, 
     !&                       HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, NATOMS)
!C        call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
!C    &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, NATOMS) 
!C End BLN }}}
      ELSE IF (WORD.EQ.'BLNGO') THEN
!C BLN-Go Model {{{
         !GOTYPE=.TRUE.
         !BLNT=.TRUE.
         !CALL READF(RK_R)
         !CALL READF(RK_THETA)
         !IF (NITEMS.GT.3) THEN
            !CALL READF(GOFACTOR)
         !ENDIF
         !ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     !&            LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS))
         !LJREP_BLN=0
         !LJATT_BLN=0
         !OPEN(UNIT=1,FILE='BLNsequence',STATUS='OLD')
         !READ(1,*) DUMMYCH
         !READ(1,*) LJREPBB, LJATTBB
         !READ(1,*) LJREPLL, LJATTLL
         !READ(1,*) LJREPNN, LJATTNN
         !READ(1,*) DUMMYCH
         !READ(1,*) DUMMYCH
         !READ(1,*) HABLN, HBBLN, HCBLN, HDBLN
         !READ(1,*) EABLN, EBBLN, ECBLN, EDBLN
         !READ(1,*) TABLN, TBBLN, TCBLN, TDBLN
         !DO J1=1,NATOMS-1
            !READ(1,'(A1)',ADVANCE='NO') BEADLETTER(J1)
         !ENDDO
         !READ(1,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
         !DO J1=1,NATOMS-3
            !READ(1,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
         !ENDDO
         !CLOSE(1)
         !WRITE(LFH,'(A,I8,A)') 'BLN sequence of ',NATOMS,' beads read:'
         !WRITE(LFH,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
         !WRITE(LFH,'(A)') ' '
         !WRITE(LFH,'(A,I8,A)') 'BLN dihedral types:'
         !WRITE(LFH,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
         !WRITE(LFH,'(A)') ' '
         !WRITE(LFH,'(A,2F15.5)') 'B-B LJ coefficients: ',LJREPBB, LJATTBB
         !WRITE(LFH,'(A,2F15.5)') 'L-L LJ coefficients: ',LJREPLL, LJATTLL
         !WRITE(LFH,'(A,2F15.5)') 'N-N LJ coefficients: ',LJREPNN, LJATTNN
         !WRITE(LFH,'(A,4F15.5)') 'Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
         !WRITE(LFH,'(A,4F15.5)') 'Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
         !WRITE(LFH,'(A,4F15.5)') 'Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
         !call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     !&                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     !&                       HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, NATOMS)
!C        call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
!C    &                       LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, NATOMS)
!C End BLN }}}
         ! }}}
         ! EDIFF STEPS UPDATES {{{
      ELSE IF (WORD.EQ.'EDIFF') THEN
         CALL READF(ECONV)
      ELSE IF (WORD.EQ.'STEPS') THEN
         NRUNS=1
         CALL READI(IX)
         MCSTEPS(1)=IX
         IF (NITEMS.GT.2) THEN
         CALL READF(XX)
         TFAC(1)=XX
         ENDIF
      ELSE IF (WORD.EQ.'UPDATES') THEN
         CALL READI(MUPDATE)
         ! }}}
         ! SLOPPYCONV TIGHTCONV DGUESS  {{{
      ELSE IF ((WORD.EQ.'BASIN').OR.(WORD.EQ.'SLOPPYCONV')) THEN
         IF (NITEMS.GT.1) CALL READF(BQMAX)
      ELSE IF ((WORD.EQ.'QMAX').OR.(WORD.EQ.'TIGHTCONV')) THEN
         CALL READF(CQMAX)
      ELSE IF (WORD.EQ.'DGUESS') THEN
         CALL READF(DGUESS)
         ! }}}
         ! SAVE DEBUG CENTRE CHANGEACCEPT MAXBFGS  MAXIT  {{{
      ELSE IF (WORD.EQ.'SAVE') THEN
         CALL READI(NSAVE)
      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.
      ELSE IF (WORD.EQ.'CENTRE') THEN
         CENT=.TRUE.
      ELSE IF (WORD.EQ.'CHANGEACCEPT') THEN
         CALL READI(IX)
         NACCEPT=IX
      ELSE IF (WORD.EQ.'MAXBFGS') THEN
         CALL READF(MAXBFGS)
      ELSE IF (WORD.EQ.'MAXIT') THEN
         CALL READI(IX)
         MAXIT=IX
         IF (NITEMS.GT.2) THEN
            CALL READI(IX)
            MAXIT2=IX
         ENDIF
         ! }}}
         ! STEP TEMPERATURE  {{{
      ELSE IF (WORD.EQ.'STEP') THEN
         NPCOUNT=NPCOUNT+1
         IF (NPCOUNT.GT.NPAR) THEN
            WRITE(LFH,'(A)') 'Number of STEP lines exceeds NPAR - quit'
            STOP
         ENDIF
         CALL READF(STEP(NPCOUNT))
         CALL READF(ASTEP(NPCOUNT))
         IF (NITEMS.GT.3) CALL READF(OSTEP(NPCOUNT))
         IF (NITEMS.GT.4) CALL READI(BLOCK(NPCOUNT))
      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         DO J1=1,NITEMS-1
            CALL READF(TEMP(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            DO J1=NITEMS,NPAR
               TEMP(J1)=TEMP(1)
            ENDDO
         ENDIF
         ! }}}
      ELSE IF (WORD.EQ.'TRACKDATA') THEN
         TRACKDATAT=.TRUE.     
      ELSE
         CALL REPORT('Unrecognized command '//WORD,.TRUE.)
         STOP
      ENDIF

