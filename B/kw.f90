      
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
           !select..case{{{
           selectcase(word)
         case('  ','NOTE','COMMENT','!','#','\\')
            GOTO 100
         ! P46 G46 MYBLN BLN BLNGO {{{
      CASE('P46','P69')
         BLNT=.TRUE.
         P46=.TRUE.
         G46=.FALSE.
         BLNTYPE="WT"
      CASE('G46')
         BLNT=.TRUE.
         G46=.TRUE.
         P46=.FALSE.
         BLNTYPE="GO"
      CASE('MYBLN')
         MYBLNT=.TRUE.
         BLNT=.TRUE.
      CASE('BLN')
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
      CASE('BLNGO')
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
      CASE('EDIFF')
          READ(ARGS(2),*) ECONV
      CASE('STEPS')
          READ(ARGS(2),*) MCSTEPS
          IF (NARGS .GT. 2) READ(ARGS(3),*) TFAC
      CASE('UPDATES')
          READ(ARGS(2),*) MUPDATE
         ! }}}
         ! SLOPPYCONV TIGHTCONV DGUESS  {{{
      CASE('BASIN','SLOPPYCONV')
          READ(ARGS(2),*) SQMAX
      CASE('QMAX','TIGHTCONV')
          READ(ARGS(2),*) FQMAX
      CASE('DGUESS')
         READ(ARGS(2),*) DGUESS
         ! }}}
         ! SAVE DEBUG CENTRE CHANGEACCEPT MAXBFGS  MAXIT  {{{
      CASE('SAVE')
         READ(ARGS(2),*) NSAVE
      CASE('DEBUG')
         DEBUG=.TRUE.
      CASE('CENTRE')
         CENT=.TRUE.
      CASE('CHANGEACCEPT')
         READ(ARGS(2),*) NACCEPT
      CASE('MAXBFGS')
         READ(ARGS(2),*) MAXBFGS
      CASE('MAXIT')
         READ(ARGS(2),*) MAXIT
         IF (NARGS.GT.2) THEN
            READ(ARGS(3),*) MAXIT2
         ENDIF
         ! }}}
         ! STEP TEMPERATURE  {{{
      CASE('STEP')
         READ(ARGS(2),*) STEP(1)
         IF (NARGS.GT.2) READ(ARGS(3),*) ASTEP(1)
         IF (NARGS.GT.3) READ(ARGS(4),*) OSTEP(1)
         IF (NARGS.GT.4) READ(ARGS(5),*) BLOCK(1)
      CASE('TEMPERATURE')
         READ(ARGS(2),*) TEMP(1)
                  ! }}}
                  !NRG {{{
      CASE('NRG')
         READ(ARGS(2),*) NRG
                  ! }}}
                  !PULL {{{
      CASE('PULL')
        PULLT=.TRUE.
        READ(ARGS(2),*) PATOM1
        READ(ARGS(3),*) PATOM2
        READ(ARGS(4),*) PFORCE
        ! }}}
        ! track {{{
      CASE('TRACKDATA')
         TRACKDATAT=.TRUE.     
      CASE('TENERGY')
         TRACKENERGY=.TRUE.     
      CASE('TBEST')
         TRACKBEST=.TRUE.     
      CASE('TMARKOV')
         TRACKMARKOV=.TRUE.     
      CASE('TXYZ')
         TXYZ=.TRUE.     
         NSQ=NARGS-1
         DO I=1,NARGS-1
            READ(ARGS(I+1),*) SQU(I)
         ENDDO
         ! }}}
       case default
        GOTO 100
      endselect

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
