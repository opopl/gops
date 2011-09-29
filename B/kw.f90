      
      MODULE KW
! {{{
USE V
USE F
USE S

IMPLICIT NONE

! }}}
CONTAINS
      SUBROUTINE KEYWORD(PROG)
      ! declarations {{{
      INTEGER PROG
      LOGICAL YESNO

      ! INPUT RELATED VARIABLES

      CHARACTER(LEN=100) :: BUFFER, WORD
      CHARACTER(LEN=20),DIMENSION(20) :: ARGS
      INTEGER :: POS, LINE
      INTEGER, PARAMETER :: DATA_FH = 15
      INTEGER :: IOS=0
      INTEGER NARGS
      integer i
      COMMON /BUFINF/ BUFFER,POS,WORD
!,IOS
      ! }}}
! body {{{
      !IOS=0

      SELECTCASE(PROG)
        CASE(1) 
          ! "data" file {{{
      CALL INQF(D_FILE,YESNO)

      IF (YESNO) THEN 
         USEKW=.true.
         CALL OPENF(DATA_FH,"<",D_FILE)
      ELSE
         USEKW=.FALSE.
         PRINT '(A)','No data file found. Using build-in values + command-line parameters...'
      ENDIF

! ios<0 end of file;
! ios=0 
! ios>0 error 
      ! }}}
      DO WHILE (IOS == 0)
      ! {{{
100     READ(DATA_FH, '(A)', IOSTAT=IOS) BUFFER
        IF (IOS == 0) THEN

           POS = SCAN(BUFFER, '    ')
           WORD = BUFFER(1:POS)
           CALL PARSE(BUFFER,' ',ARGS,NARGS)
           BUFFER = BUFFER(POS+1:)
           WORD=ARGS(1)
           !if..else loop{{{
      IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'.OR.WORD.EQ.'!'&
     &                          .OR. WORD .EQ. '\\') THEN 
         GOTO 100

         ! P46 G46 MYBLN BLN BLNGO {{{
      ELSE IF (WORD.EQ.'P46'.or.'P69') THEN
         BLNT=.TRUE.
         P46=.TRUE.
         G46=.FALSE.
         BLNTYPE="WT"
      ELSE IF (WORD.EQ.'G46') THEN
         BLNT=.TRUE.
         G46=.TRUE.
         P46=.FALSE.
         BLNTYPE="GO"
      ELSE IF (WORD.EQ.'MYBLN') THEN
         MYBLNT=.TRUE.
         BLNT=.TRUE.
      ELSE IF (WORD.EQ.'BLN') THEN
!        ! BLN {{{
         !BLNT=.TRUE.
         !READ(ARGS(2),*) RK_R
         !READ(ARGS(2),*) RK_THETA
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
         !READ(ARGS(2),*) RK_R
         !READ(ARGS(2),*) RK_THETA
         !IF (NITEMS.GT.3) THEN
            !READ(ARGS(2),*) GOFACTOR
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
          READ(ARGS(2),*) ECONV
      ELSE IF (WORD.EQ.'STEPS') THEN
          READ(ARGS(2),*) MCSTEPS
          IF (NARGS .GT. 2) READ(ARGS(3),*) TFAC
      ELSE IF (WORD.EQ.'UPDATES') THEN
          READ(ARGS(2),*) MUPDATE
         ! }}}
         ! SLOPPYCONV TIGHTCONV DGUESS  {{{
      ELSE IF ((WORD.EQ.'BASIN').OR.(WORD.EQ.'SLOPPYCONV')) THEN
          READ(ARGS(2),*) SQMAX
      ELSE IF ((WORD.EQ.'QMAX').OR.(WORD.EQ.'TIGHTCONV')) THEN
          READ(ARGS(2),*) FQMAX
      ELSE IF (WORD.EQ.'DGUESS') THEN
         READ(ARGS(2),*) DGUESS
         ! }}}
         ! SAVE DEBUG CENTRE CHANGEACCEPT MAXBFGS  MAXIT  {{{
      ELSE IF (WORD.EQ.'SAVE') THEN
         READ(ARGS(2),*) NSAVE
      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.
      ELSE IF (WORD.EQ.'CENTRE') THEN
         CENT=.TRUE.
      ELSE IF (WORD.EQ.'CHANGEACCEPT') THEN
         READ(ARGS(2),*) NACCEPT
      ELSE IF (WORD.EQ.'MAXBFGS') THEN
         READ(ARGS(2),*) MAXBFGS
      ELSE IF (WORD.EQ.'MAXIT') THEN
         READ(ARGS(2),*) MAXIT
         IF (NARGS.GT.2) THEN
            READ(ARGS(3),*) MAXIT2
         ENDIF
         ! }}}
         ! STEP TEMPERATURE  {{{
      ELSE IF (WORD.EQ.'STEP') THEN
         READ(ARGS(2),*) STEP(1)
         IF (NARGS.GT.2) READ(ARGS(3),*) ASTEP(1)
         IF (NARGS.GT.3) READ(ARGS(4),*) OSTEP(1)
         IF (NARGS.GT.4) READ(ARGS(5),*) BLOCK(1)
      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         READ(ARGS(2),*) TEMP(1)
                  ! }}}
                  !NRG {{{
      ELSE IF (WORD.EQ.'NRG') THEN
         READ(ARGS(2),*) NRG
                  ! }}}
                  !PULL {{{
      ELSE IF (WORD.EQ.'PULL') THEN
        PULLT=.TRUE.
        READ(ARGS(2),*) PATOM1
        READ(ARGS(3),*) PATOM2
        READ(ARGS(4),*) PFORCE
        ! }}}
        ! track {{{
      ELSE IF (WORD.EQ.'TRACKDATA') THEN
         TRACKDATAT=.TRUE.     
      ELSE IF (WORD.EQ.'TENERGY') THEN
         TRACKENERGY=.TRUE.     
      ELSE IF (WORD.EQ.'TBEST') THEN
         TRACKBEST=.TRUE.     
      ELSE IF (WORD.EQ.'TMARKOV') THEN
         TRACKMARKOV=.TRUE.     
      ELSE IF (WORD.EQ.'TXYZ') THEN
         TXYZ=.TRUE.     
         NSQ=NARGS-1
         DO I=1,NARGS-1
            READ(ARGS(I+1),*) SQU(I)
         ENDDO
         ! }}}
      ELSE
        GOTO 100
      ENDIF

           ! }}}
       END IF
       ! }}}
      END DO

      ENDSELECT

      IF (USEKW) THEN
        CLOSE(DATA_FH)
      ENDIF

      RETURN
      ! }}}
      END SUBROUTINE KEYWORD
! }}}

      ENDMODULE
