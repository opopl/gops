
MODULE FUNC

USE V

IMPLICIT NONE

SAVE

CONTAINS

SUBROUTINE IO
! {{{
      INTEGER IA,J,K

      CALL OPENF(COORDS_FH,'O','coords')
      REWIND(COORDS_FH)

      DO IA=1,NATOMS
              READ(COORDS_FH,*) COORDS(IA,1:3)
      ENDDO

      CLOSE(COORDS_FH)

      WRITE(LFH,20) 
20    FORMAT('Initial coordinates:')
30    FORMAT(3F20.10)

      DO IA=1,NATOMS
              WRITE(LFH,30) COORDS(IA,1:3)
      ENDDO

      IF (P46) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' 3-COLOUR, 46 BEAD MODEL POLYPEPTIDE'
      ELSE IF (BLNT) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' BEAD BLN MODEL'
      ENDIF
      
      IF (RADIUS.EQ.0.0D0) THEN
        RADIUS=2.0D0+(3.0D0*NATOMS/17.77153175D0)**(1.0D0/3.0D0)
        IF (P46) THEN
           RADIUS=RADIUS*3.0D0
        ENDIF
      ENDIF

      if (LBFGST) then
         WRITE(LFH,'(A)') 'Nocedal LBFGS minimization'
         WRITE(LFH,'(A,I6)') 'Number of updates before reset in LBFGS=',M_LBFGS
         WRITE(LFH,'(A,F20.10)') 'Maximum step size=',MAXBFGS
         WRITE(LFH,'(A,G12.4)') 'Guess for initial diagonal elements in LBFGS=',DGUESS
      ENDIF

      WRITE(LFH,'(A,F15.10)') 'Final quench tolerance for RMS gradient ',FQMAX
      WRITE(LFH,'(A,F15.10)') 'Energy difference criterion for minima=',EDIFF
      WRITE(LFH,'(A,I5)') 'Maximum number of iterations (sloppy and final quenches)  ',MAXIT

      IF (DEBUG) THEN
         WRITE(LFH,160) 
160      FORMAT('Debug printing is on')
      ENDIF
       
      WRITE(LFH, '(A,G20.10)') 'Maximum allowed energy rise during a minimisation=',MAXERISE

      IF (TARGET) THEN
         WRITE(LFH,'(A)',ADVANCE='NO') 'Target energies: '
         WRITE(LFH,'(F20.10)',ADVANCE='NO') (TARGETS(J),J=1,NTARGETS)
         WRITE(LFH,'(A)') ' '
      ENDIF
! }}}
END SUBROUTINE  

SUBROUTINE INITVARS
! logicals {{{
PULLT=.TRUE.
P46=.FALSE.
G46=.FALSE.
BLNT=.FALSE.
TARGET=.FALSE.
TRACKDATAT=.FALSE.
DEBUG=.FALSE.
LBFGST=.TRUE.
RMST=.FALSE.
! whether we are doing a final quench
FQFLAG=.FALSE.
! }}}
! files {{{
LFH=FH+1
ENERGY_FH=FH+2
MARKOV_FH=FH+3
BEST_FH=FH+4
PAIRDIST_FH=FH+5
COORDS_FH=FH+6
DATA_FH=FH+7
RMSD_FH=FH+8
! }}}
! {{{
NSAVE=10
MAXNARGS=20
! Number of LBFGS updates
M_LBFGS=4
! Maximum BFGS step size
MAXBFGS=0.4D0
DGUESS=0.1D0
BFGS=.FALSE.
LBFGST=.TRUE.

! DGUESS: Guess for initial diagonal elements in LBFGS
DGUESS=0.0D0

TFAC=1.0D0
EDIFF=0.02D0
ACCRAT=0.5D0
TEMP=0.035D0
RADIUS=0.0D0
! Maximum allowed energy rise during a minimisation
MAXERISE=1.0D-10
! Maximum allowed energy fall during a minimisation
MAXEFALL=-HUGE(ONE)
! Used in ACCREJ
FAC0=1.05D0
! "Fixing" option (regarding STEP, TEMP and accept ratio for quenches) 
FIXOPT='T'

STEP=0.3D0
OSTEP=0.3D0
ASTEP=0.3D0

MCSTEPS=10000
          
NACCEPT=50
NRELAX=0
! }}}
END SUBROUTINE

SUBROUTINE COUNTATOMS
! Declarations {{{ 
      IMPLICIT NONE

      INTEGER :: EOF
      LOGICAL :: YESNO
! }}} 
! {{{
      CALL INQF('coords',YESNO)

      IF (YESNO) THEN
         CALL OPENF(COORDS_FH,'O','coords')
         DO
            READ(COORDS_FH,*,IOSTAT=EOF)
            IF (EOF==0) THEN
               NATOMS=NATOMS+1 
            ELSE
               EXIT
            ENDIF
         ENDDO
      ELSE
         PRINT '(A)','ERROR - no coords file'
         STOP
      ENDIF

      CLOSE(COORDS_FH)
! }}}      
END SUBROUTINE COUNTATOMS

SUBROUTINE MYSYSTEM(STATUS,DEBUG,JOBSTRING)
! {{{
USE PORFUNCS

IMPLICIT NONE

LOGICAL DEBUG
INTEGER STATUS
CHARACTER(LEN=*) JOBSTRING

IF (DEBUG) WRITE(*,'(A)') 'mysystem> '//trim(adjustl(jobstring)) 
CALL SYSTEM_SUBR(JOBSTRING,STATUS)

! IF (DEBUG) PRINT '(A,I6)','command '//JOBSTRING//' exit status=',STATUS
! IF (STATUS.NE.0) PRINT '(A,I8)','mysystem> WARNING - '//JOBSTRING//' exit status=',STATUS

RETURN
! }}}
END SUBROUTINE MYSYSTEM

SUBROUTINE INQF(FILENAME,YESNO)
! {{{
LOGICAL,INTENT(OUT) :: YESNO
CHARACTER(LEN=*),INTENT(IN) :: FILENAME
CHARACTER(LEN=100) :: FLN

FLN=TRIM(ADJUSTL(FILENAME))
INQUIRE(FILE=FLN,EXIST=YESNO)
! }}}
END SUBROUTINE INQF

include "mylbfgs.i.f90"

!> @name OPENF
!! @brief open files 

SUBROUTINE OPENF(FILEHANDLE,MODE,FILENAME)
! {{{

INTEGER, INTENT(IN) :: FILEHANDLE
CHARACTER (LEN=*), INTENT(IN) :: FILENAME
CHARACTER (LEN=*), INTENT(IN) :: MODE
CHARACTER (LEN=100) :: FLN

FLN=TRIM(ADJUSTL(FILENAME))

SELECTCASE(MODE)
        CASE(">")
                OPEN(UNIT=FILEHANDLE,FILE=FILENAME,STATUS="UNKNOWN",FORM="FORMATTED")
        CASE("O")
                OPEN(UNIT=FILEHANDLE,FILE=FILENAME,STATUS="OLD")
        CASE("<")
                OPEN(UNIT=FILEHANDLE,FILE=FILENAME,STATUS="OLD",ACTION="READ")
        CASE(">>")
                OPEN(FILEHANDLE,FILE=FILENAME,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
        CASE("RW>>")
                OPEN(FILEHANDLE,FILE=FILENAME,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND',ACTION="READWRITE")
        case("DA")
                OPEN(FILEHANDLE,FILE=FILENAME,ACCESS="DIRECT",STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=8*3*NATOMS)
                
ENDSELECT
! }}}
ENDSUBROUTINE OPENF

SUBROUTINE TAKESTEP
      ! {{{

      USE V

      IMPLICIT NONE

      INTEGER IA
      DOUBLE PRECISION ::  RND(3)

      DO IA=1,NATOMS
         CALL GETRND(RND,3,-1.0D0,1.0D0)
         COORDS(IA,1:3)=COORDS(IA,1:3)+STEP*RND(1:3)
      ENDDO
      
      RETURN
      ! }}}
END SUBROUTINE 

SUBROUTINE CENTRE2(R)
! {{{
USE V

IMPLICIT NONE

DOUBLE PRECISION, INTENT(INOUT) :: R(:,:)
INTEGER I,K

RMASS=SUM(R,DIM=1)/SIZE(R,DIM=1)

do K=1,3
        R(:,K)=R(:,K)-RMASS(K)
ENDDO

IF (DEBUG) WRITE(LFH,'(A,3G20.10)') 'centre2> centre of mass reset to the origin from ',RMASS

! }}}
END SUBROUTINE

! doxygen - GETRND {{{
!> @name         GETRND
! 
!> @brief        Get an array of random numbers inside the interval [XMIN,XMAX]
!
!> @param[in]    N              - dimension of the array RND
!> @param[out]   RND            - the generated array of random numbers 
!> @param[in]    XMIN,XMAX      
!
! }}}
SUBROUTINE GETRND(RND,N,XMIN,XMAX)
! {{{
IMPLICIT NONE

! random number vector
DOUBLE PRECISION :: RND(:),XMIN,XMAX,DX
! dimension of RND(:)
INTEGER N,I

DX=XMAX-XMIN

DO I=1,N 
        RND(I)=XMAX-DX*DPRAND()
ENDDO

RETURN
! }}}
END SUBROUTINE 

SUBROUTINE SDPRND (ISEED)
        ! declarations {{{
        DOUBLE PRECISION XMOD, YMOD, POLY(101), OTHER, OFFSET, X
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0)
        INTEGER ISEED, INDEX, IX, IY, IZ, I
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
        ! }}}
        ! subroutine body {{{
!
!   ISEED should be set to an integer between 0 and 9999 inclusive;
!   a value of 0 will initialise the generator only if it has not
!   already been done.
!
        IF (INITAL .OR. ISEED .NE. 0) THEN
            INITAL = .FALSE.
        ELSE
            RETURN
        END IF
!
!   INDEX must be initialised to an integer between 1 and 101
!   inclusive, POLY(1...N) to integers between 0 and 1000009710
!   inclusive (not all 0), and OTHER to a non-negative proper fraction
!   with denominator 33554432.  It uses the Wichmann-Hill generator to
!   do this.
!
        IX = MOD(ABS(ISEED),10000)+1
        IY = 2*IX+1
        IZ = 3*IX+1
        DO 10 I = -10,101
            IF (I .GE. 1) POLY(I) = AINT(XMOD*X)
            IX = MOD(171*IX,30269)
            IY = MOD(172*IY,30307)
            IZ = MOD(170*IZ,30323)
            X = MOD(DBLE(IX)/30269.0D0+DBLE(IY)/30307.0D0+DBLE(IZ)/30323.0D0,1.0D0)
  10    CONTINUE
        OTHER = AINT(YMOD*X)/YMOD
        OFFSET = 1.0D0/YMOD
        INDEX = 1
        ! }}}
END SUBROUTINE

FUNCTION DPRAND()
        ! DECLARATIONS {{{
        DOUBLE PRECISION XMOD, YMOD, XMOD2, XMOD4, TINY, POLY(101), DPRAND, &
         OTHER, OFFSET, X, Y
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0, &
        XMOD2 = 2000019422.0D0, XMOD4 = 4000038844.0D0,&
        TINY = 1.0D-17)
        INTEGER INDEX, N
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
        ! }}}
        ! SUBROUTINE BODY {{{
!
!   THIS RETURNS A UNIFORM (0,1) RANDOM NUMBER, WITH EXTREMELY GOOD
!   UNIFORMITY PROPERTIES.  IT ASSUMES THAT DOUBLE PRECISION PROVIDES
!   AT LEAST 33 BITS OF ACCURACY, AND USES A POWER OF TWO BASE.
!
        IF (INITAL) THEN
            CALL SDPRND (0)
            INITAL = .FALSE.
        END IF
!
!   SEE [KNUTH] FOR WHY THIS IMPLEMENTS THE ALGORITHM DESCRIBED IN
!   THE PAPER.  NOTE THAT THIS CODE IS TUNED FOR MACHINES WITH FAST
!   DOUBLE PRECISION, BUT SLOW MULTIPLY AND DIVIDE; MANY, MANY OTHER
!   OPTIONS ARE POSSIBLE.
!
        N = INDEX-64
        IF (N .LE. 0) N = N+101
        X = POLY(INDEX)+POLY(INDEX)
        X = XMOD4-POLY(N)-POLY(N)-X-X-POLY(INDEX)
        IF (X .LT. 0.0D0) THEN
            IF (X .LT. -XMOD) X = X+XMOD2
            IF (X .LT. 0.0D0) X = X+XMOD
        ELSE
            IF (X .GE. XMOD2) THEN
                X = X-XMOD2
                IF (X .GE. XMOD) X = X-XMOD
            END IF
            IF (X .GE. XMOD) X = X-XMOD
        END IF
        POLY(INDEX) = X
        INDEX = INDEX+1
        IF (INDEX .GT. 101) INDEX = INDEX-101
!
!   ADD IN THE SECOND GENERATOR MODULO 1, AND FORCE TO BE NON-ZERO.
!   THE RESTRICTED RANGES LARGELY CANCEL THEMSELVES OUT.
!
   10   Y = 37.0D0*OTHER+OFFSET
        OTHER = Y-AINT(Y)
        IF (OTHER .EQ. 0.0D0) GO TO 10
        X = X/XMOD+OTHER
        IF (X .GE. 1.0D0) X = X-1.0D0
        DPRAND = X+TINY
        ! }}}
END FUNCTION
 
ENDMODULE





