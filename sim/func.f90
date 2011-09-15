
MODULE FUNC

USE V
USE PORFUNCS
USE BLN

IMPLICIT NONE

SAVE

CONTAINS

! doxygen - POTENTIAL   {{{
!>
!> @name POTENTIAL
!> @brief Given the input coordinates X, calculate the energy (EREAL), the
!> gradient (GRADX), the root-mean-square force (RMS). The logical flags DOGRAD
!> and DOHESS control whether one needs to calculate the gradient and the Hessian
!
!> @param[in]  X        dp(N)    input coordinates 
!> @param[out] GRADX    dp(N)    gradient 
!> @param[out] EREAL    dp         energy 
!> @param[out] RMS      dp         RMS 
! }}}
      SUBROUTINE POTENTIAL(X,EREAL,GRADX,RMS,DOGRAD,DOHESS)
! declarations {{{

      IMPLICIT NONE

      ! subroutine parameters  {{{

      ! conventions for energies:
      !     EREAL  - always total energy; 
      !     EO - other energies

      DOUBLE PRECISION, INTENT(OUT) :: EREAL, GRADX(:)
      DOUBLE PRECISION, INTENT(IN) :: X(:)
      LOGICAL DOGRAD, DOHESS
      DOUBLE PRECISION, INTENT(OUT) :: RMS
      ! }}}
      ! local parameters  {{{

      INTEGER NX,NR

      DOUBLE PRECISION,ALLOCATABLE :: R(:,:)
      DOUBLE PRECISION,ALLOCATABLE :: GRAD(:,:),HESS(:,:)
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: E   
      ! }}}
      ! }}}
! subroutine body {{{

      ! NX is the number of atoms multiplied by the dimension, i.e., NATOMS*3
      ! NR is the number of atoms, i.e., NATOMS
      NX=SIZE(X)
      NR=NX/3

      ALLOCATE(R(NR,3),GRAD(NR,3),HESS(NX,NX))

      R=RESHAPE(X,(/ NR,3 /))

      IF (BLNT .OR. PULLT) THEN 
            CALL EBLN(NR,R,EA,GRAD,HESS,PTYPE,DOGRAD,DOHESS)
            EREAL=EA(1)
      ENDIF

      IF (PULLT) THEN
         EA(6)=-PFORCE*R(PATOM1,3)+PFORCE*R(PATOM2,3)
         EREAL=EREAL+EA(6)
         GRAD(PATOM1,3)=GRAD(PATOM1,3)-PFORCE
         GRAD(PATOM2,3)=GRAD(PATOM2,3)+PFORCE
      ENDIF

      GRADX=PACK(GRAD,.TRUE.)
      RMS=DSQRT(SUM(GRADX**2)/(1.0D0*NX))
      CALL ECHO_S
      WRITE(*,*) "potential> NX=  ",NX,"RMS= ",RMS 
      write(*,*) "potential> Total energy: ",EA(1)
      write(*,*) "potential> Non-bonded ",EA(2)
      write(*,*) "potential> Bonded ",EA(3)
      write(*,*) "potential> Bond angles ",EA(4)
      write(*,*) "potential> Torsional angles ",EA(5)
      write(*,*) "potential> Force term: -F*(z1-z2)=  ",EA(6)
      CALL ECHO_S

      DEALLOCATE(R,GRAD,HESS)

      RETURN
! }}}
      END SUBROUTINE POTENTIAL

SUBROUTINE GETMODEL
! {{{

IF (P46) THEN 
  MODEL="P46 THREE-COLOUR OFF-LATTICE PROTEIN MODEL, WILD-TYPE"
ELSEIF(G46) THEN 
  MODEL="P46 THREE-COLOUR OFF-LATTICE PROTEIN MODEL, GO-LIKE"
ELSEIF(BLNT) THEN
  MODEL="GENERAL BLN MODEL"
ENDIF
! }}}
ENDSUBROUTINE GETMODEL

! RCOORDS IO PRINTVARS SETVARS INITVARS {{{

SUBROUTINE AM

allocate(SCREENC(NATOMS,3))
allocate(EA(NE))

ENDSUBROUTINE AM

! initialize variables
SUBROUTINE INITVARS
! subroutine body {{{
! logicals {{{
BFGST=.FALSE.
LBFGST=.TRUE.
PULLT=.TRUE.
P46=.FALSE.
G46=.TRUE.
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
! other {{{
NSAVE=10            ! number of saved lowest-energy geometries

NATOMS=46

PATOM1=1
PATOM2=NATOMS

M_LBFGS=4           ! Number of LBFGS updates
MAXBFGS=0.4D0       ! Maximum BFGS step size
DGUESS=0.1D0        ! DGUESS: Guess for initial diagonal elements in LBFGS

FQMAX=1.0D-5        ! FQMAX: same meaning as for SQMAX, but for final quenches only.
SQMAX=1.0D-3        ! SQMAX: convergence criterion for the RMS force in the basin-hopping quenches.
                    ! note: used in QUENCH() 

TFAC=1.0D0
EDIFF=0.02D0
ACCRAT=0.5D0
TEMP=0.035D0        ! Temperature
RADIUS=0.0D0
! maximum number of iterations allowed in conjugate gradient searches
MAXIT=500
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
CALL SETVARS
! }}}
END SUBROUTINE INITVARS

! print variables 
SUBROUTINE PRINTVARS
! {{{
! vars                                                                       {{{
character(100) fmt(10)
character(40) s
integer i
! }}}
include '../include/fmt.i.f90'

CALL ECHO_S
write(*,10) "PARAMETER VALUES" !                                             {{{
write(*,1)  "PARAMETER DESCRIPTION",         "NAME",                  "VALUE"
! pd:general                                                                 {{{
write(*,11) s_stars
write(*,10) "GENERAL" 
write(*,11) s_stars

call getmodel
write(s,*) "Model" 
!s=adjustl(s)
!write(*,1)  "Model"                                         ,"Model",   trim(model)
write(*,*)  trim(model)
write(*,11) s_stars
write(*,3)  "Number of particles",                          "NATOMS",     NATOMS
write(*,2)  "Container radius",                          "RADIUS",     RADIUS
write(*,3)  "Number of saved lowest energy geometries",     "NSAVE",      NSAVE
write(*,3)  "Number of basin-hopping steps",                "MCSTEPS",    MCSTEPS
write(*,2)  "Temperature",                                  "TEMP",       TEMP
write(*,2)  "Acceptance ratio",                             "ACCRAT",     ACCRAT
write(*,2)  "Energy difference criterion for minima",       "EDIFF",      EDIFF
write(*,2)  "Final quench tolerance for RMS gradient ",     "FQMAX",      FQMAX
write(*,2)  "Quench convergence criterion for RMS gradient ", "SQMAX",      SQMAX
write(*,3)  "Maximum number of iterations",   "MAXIT", MAXIT
write(*,11) "(sloppy and final quenches)"
write(*,2)  "",                                             "TFAC",       TFAC
Write(*,3)  "",                                             "NACCEPT",    NACCEPT
write(*,3)  "",                                             "NRELAX",     NRELAX
write(*,11) s_stars !                                                        }}}
! pd:lbfgs                                                                     {{{
write(*,10) "LBFGS parameters"
write(*,11) s_stars
write(*,3) "Number of LBFGS updates",   "M_LBFGS",  M_LBFGS
write(*,2) "Maximum BFGS step size",   "MAXBFGS",   MAXBFGS
write(*,2) "Guess for initial diagonal elements in BFGS", "DGUESS", DGUESS
write(*,11) s_stars
! }}} 
! }}}
CALL ECHO_S
! }}}
END SUBROUTINE PRINTVARS

! set variables 
SUBROUTINE SETVARS
! {{{
IF (P46 .OR. G46) THEN
  NATOMS=46
ENDIF
IF (P46) PTYPE="WT"
IF (G46) PTYPE="GO"

NE=20

! radius {{{
!IF (RADIUS.EQ.0.0D0) THEN
    RADIUS=2.0D0+(3.0D0*NATOMS/17.77153175D0)**(1.0D0/3.0D0)
    IF (P46) THEN
        RADIUS=RADIUS*3.0D0
    ENDIF
    RADIUS=10.0D0
!ENDIF
! }}}


! }}}
END SUBROUTINE SETVARS


! doxygen RCOORDS {{{
!> @name RCOORDS
!> @brief Generate a random set of coordinates, based on:
!> @param[in] NATOMS number of particles
!> @param[in] RADIUS container radius
!> @param[out] COORDS  randomly generated coordinates
! }}}
SUBROUTINE RCOORDS(NATOMS,RADIUS,COORDS)
! declarations {{{
! subroutine parameters 
INTEGER, INTENT(IN) :: NATOMS
DOUBLE PRECISION,INTENT(IN) :: RADIUS 
DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: COORDS

! local parameters 
DOUBLE PRECISION :: SR3
DOUBLE PRECISION, ALLOCATABLE :: RND(:)
! }}}
! {{{
SR3=DSQRT(3.0D0)
ALLOCATE(RND(3*NATOMS))
CALL GETRND(RND,3*NATOMS,-1.0D0,1.0D0)
RND=RADIUS*RND/SR3
COORDS=RESHAPE(RND,(/ NATOMS, 3 /))
DEALLOCATE(RND)
! }}}
END SUBROUTINE RCOORDS

SUBROUTINE IO
! {{{
      INTEGER IA,J,K

      IF (TRACKDATAT) THEN
         CALL OPENF(ENERGY_FH,">>","energy.dat")
         CALL OPENF(MARKOV_FH,">>","markov.dat")
         CALL OPENF(BEST_FH,">>","best.dat")
         !IF (RMST) 
        CALL OPENF(RMSD_FH,">>","rmsd.dat")
      ENDIF

      WRITE(LFH,20) 
20    FORMAT('Initial coordinates:')
30    FORMAT(3F20.10)

      WRITE(LFH,30) SCREENC

      IF (P46) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' 3-COLOUR, 46 BEAD MODEL POLYPEPTIDE'
      ELSE IF (BLNT) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' BEAD BLN MODEL'
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
END SUBROUTINE IO



! }}}
! INQF OPENF {{{
SUBROUTINE INQF(FILENAME,YESNO)
! {{{
LOGICAL,INTENT(OUT) :: YESNO
CHARACTER(LEN=*),INTENT(IN) :: FILENAME
CHARACTER(LEN=100) :: FLN

FLN=TRIM(ADJUSTL(FILENAME))
INQUIRE(FILE=FLN,EXIST=YESNO)
! }}}
END SUBROUTINE INQF

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
! }}}
! GETRND SDPRND DPRAND {{{

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
! }}}
! LBFGS {{{

! Doxygen - MYLBFGS {{{
!>
!> @name MYLBFGS
!
!> @brief Adaptation of the original LBFGS routine
!
!> @param[in] N       
!>             is an INTEGER variable that must be set by the user to the
!>             number of variables. It is not altered by the routine.
!>             Restriction: N>0.
!
!> @param[in] M   
!>             is an INTEGER variable that must be set by the user to
!>             the number of corrections used in the BFGS update. It
!>             is not altered by the routine. Values of M less than 3 are
!>             not recommended; large values of M will result in excessive
!>             computing time. 3<= M <=7 is recommended. Restriction: M>0.
!
!> @param[out] MFLAG 
!>             returns .TRUE. is the algorithm has converged; .FALSE. otherwise
!
!> @param[in] DIAGCO 
!>             is a LOGICAL variable that must be set to .TRUE. if the
!>             user wishes to provide the diagonal matrix Hk0 at each
!>             iteration. Otherwise it should be set to .FALSE., in which
!>             case  LBFGS will use a default value described below. If
!>             DIAGCO is set to .TRUE. the routine will return at each
!>             iteration of the algorithm with IFLAG=2, and the diagonal
!>             matrix Hk0  must be provided in the array DIAG.
!
!> @param[in] EPS 
!>             is a positive DOUBLE PRECISION variable that must be set by
!>             the user, and determines the accuracy with which the solution
!>             is to be found. The subroutine terminates when
!>
!>                         ||G|| < EPS max(1,||X||),
!>
!>             where ||.|| denotes the Euclidean norm.
!
!> @param[inout] R is a vector of coordinates. On initial entry
!>             it must be set by the user to the values of the initial
!>             estimate of the solution vector. On exit with IFLAG=0, it
!>             contains the values of the variables at the best point
!>             found (usually a solution).
!
!
! }}}
SUBROUTINE MYLBFGS(X,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET)
! ====================================
! subroutine parameters  {{{

IMPLICIT NONE

DOUBLE PRECISION, INTENT(INOUT) :: X(:)
LOGICAL,INTENT(IN) :: DIAGCO
DOUBLE PRECISION,INTENT(IN) :: EPS
LOGICAL,INTENT(OUT) :: MFLAG
DOUBLE PRECISION,INTENT(OUT) :: ENERGY
INTEGER,INTENT(IN) :: ITMAX
INTEGER,INTENT(OUT) :: ITDONE
LOGICAL,INTENT(IN) :: RESET
! }}}
! local parameters {{{
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DIAG, W, WTEMP, GRADX, GNEW, XSAVE
! ALPHA, RHO: M coefficients in the LBFGS algorithm
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ALPHA, RHO
! WSS - for storing step directions
! WDG - for storing gradient differences
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: WSS,WDG
INTEGER ITER, NFAIL, CP, NDECREASE
DOUBLE PRECISION :: GNORM,SLENGTH,E,ENEW,YS,YY,SQ,YR,BETA,POTEL, DDOT, STP
DOUBLE PRECISION :: DOT1,DOT2,OVERLAP, DUMMY, QE

INTEGER :: BOUND, POINT

INTEGER :: N,M,K,IX,IXMIN,IXMAX
COMMON /MYPOT/ QE

SAVE W,DIAG,ITER
! 
N=3*NATOMS
M=M_LBFGS
!
!   }}}
! ====================================
! Labels:  {{{
!         10  - termination test
!         30  - compute the new step and gradient change
! }}}
! ====================================
! subroutine body  {{{

! Initializations {{{

IF(.NOT.ALLOCATED(W)) ALLOCATE(W(N))
IF(.NOT.ALLOCATED(WTEMP)) ALLOCATE(WTEMP(N))
IF(.NOT.ALLOCATED(GRADX)) ALLOCATE(GRADX(N))
IF(.NOT.ALLOCATED(GNEW)) ALLOCATE(GNEW(N))
IF(.NOT.ALLOCATED(DIAG)) ALLOCATE(DIAG(N))
IF(.NOT.ALLOCATED(XSAVE)) ALLOCATE(XSAVE(N))

ALLOCATE(RHO(M),ALPHA(M),WSS(N,M),WDG(N,M))

NFAIL=0 ; ITDONE=0 ; IF (RESET) ITER=0 

IF (DEBUG) THEN
    IF (RESET) WRITE(LFH,'(A)')            'mylbfgs> Resetting LBFGS minimiser'
    IF (.NOT.RESET) WRITE(LFH,'(A)')       'mylbfgs> Not resetting LBFGS minimiser'
ENDIF

!        evaluate gradient (.TRUE.) but not Hessian (.FALSE.)
CALL POTENTIAL(X,QE,GRADX,RMS,.TRUE.,.FALSE.)
write(*,*) 'mylbfgs: init, after potential(), QE=',QE

IF (DEBUG) WRITE(LFH,101) ' Energy and RMS force=', E,RMS, &
    ' after ',ITDONE,' LBFGS steps'
101   FORMAT(A,F20.10,G20.10,A,I6,A)
102   FORMAT(A,I10,A,I10,A)

! }}}
!  Termination test (label 10){{{

10    CALL FLUSH(LFH)

MFLAG=.FALSE.
IF (RMS.LE.EPS) THEN
    MFLAG=.TRUE.
    IF (MFLAG) THEN
        IF (DEBUG) WRITE(LFH,101) ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'
        RETURN
    ENDIF
ENDIF

IF (ITDONE.EQ.ITMAX) THEN
    IF (DEBUG) WRITE(LFH,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
    RETURN
ENDIF
! }}}

IF (ITER.EQ.0) THEN
	! {{{
	IF (N.LE.0.OR.M.LE.0) THEN
    	WRITE(LFH,240)
	    240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
	    STOP
	ENDIF
	
	POINT=0
	MFLAG=.FALSE.
    ! Set the DIAG matrix depending on DIAGCO {{{
	IF (DIAGCO) THEN
	    WRITE(LFH,'(A)') 'using estimate of the inverse diagonal elements'
	    DO IX=1,N
		    IF (DIAG(IX).LE.0.0D0) THEN
	    	    WRITE(LFH,235) IX
		        235 FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
		        STOP
	        ENDIF
	    ENDDO
	ELSE
	    DIAG=DGUESS
	ENDIF
    ! }}}
    ! Rules for storage: {{{
	!
	!     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
	!     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
	!
	
	! W     =>      W 
	! RHO   =>      N,...,N+M      
	! ALPHA =>      N+M+1,...,N+2M
	! WSS   =>      storage of the last M steps 
	! WDG   =>      storage of the last M gradient differences 
    ! }}}
	
	!  NR step for diagonal inverse Hessian
	
	DO IX=1,N
		DUMMY=-GRADX(IX)*DIAG(IX)
		WSS(IX,1)=DUMMY
		W(IX)=DUMMY
	ENDDO

	GNORM=DSQRT(DDOT(N,GRADX,1,GRADX,1))
    !  Guess for the step length
	STP=MIN(1.0D0/GNORM,GNORM)              
	! }}}
ELSE 
	! {{{
	BOUND=MIN(ITER,M) 
	YS= DDOT(N,WDG(1,1+POINT),1,WSS(1,1+POINT),1)
	!
	!  Update estimate of diagonal inverse Hessian elements {{{
	!
	IF (.NOT.DIAGCO) THEN
		YY= DDOT(N,WDG(1,1+POINT),1,WDG(1,1+POINT),1)
		IF (YY.EQ.0.0D0) THEN
			WRITE(LFH,'(A)') 'WARNING, resetting YY to one in mylbfgs'
			YY=1.0D0
	    ENDIF
		IF (YS.EQ.0.0D0) THEN
		    WRITE(LFH,'(A)') 'WARNING, resetting YS to one in mylbfgs'
		    YS=1.0D0
		ENDIF
	    DIAG= YS/YY
	ELSE
		WRITE(LFH,'(A)') 'using estimate of the inverse diagonal elements'
		DO IX=1,N
			IF (DIAG(IX).LE.0.0D0) THEN
			    WRITE(LFH,235) IX
			    STOP
	    	ENDIF
	    ENDDO
	ENDIF
	! }}}
    ! -HG computation {{{
	!
	!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
	!     "Updating quasi-Newton matrices with limited storage",
	!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
	!     ---------------------------------------------------------
	!
	CP= POINT
	IF (POINT.EQ.0) CP=M
	RHO(CP)= 1.0D0/YS
	W= -GRADX

	CP= POINT
	DO IX= 1,BOUND
		CP=CP-1
		IF (CP.EQ.-1) CP=M-1
        IF (DEBUG_LBFGS) WRITE(*,*) 'IX= ',IX,'CP= ',CP
        IF (DEBUG_LBFGS) WRITE(*,*) "SQ=DDOT(N,WSS(1,CP+1),1,W,1)"
		SQ=DDOT(N,WSS(1,CP+1),1,W,1)
		ALPHA(CP+1)=RHO(CP+1)*SQ
		CALL DAXPY(N,-ALPHA(CP+1),WDG(1,CP+1),1,W,1)
	ENDDO
	
	W=DIAG*W
	
	DO IX=1,BOUND
		YR=DDOT(N,WSS(1,CP+1),1,W,1)
		BETA=RHO(CP+1)*YR
		BETA=ALPHA(CP+1)-BETA
		CALL DAXPY(N,BETA,WSS(1,CP+1),1,W,1)
		CP=CP+1
		IF (CP.EQ.M) CP=0
	ENDDO
	STP=1.0D0  
    ! }}}
	! }}}
ENDIF

!  Store the new search direction (160)

IF (ITER.GT.0) THEN
    WSS(:,1+POINT)=W
ENDIF

!  test for overlap {{{
!
!  Overflow has occasionally occurred here.
!  We only need the sign of the overlap, so use a temporary array with
!  reduced elements.
!
DUMMY=1.0D0
DO IX=1,N
    IF (ABS(W(IX)).GT.DUMMY) DUMMY=ABS(W(IX))
ENDDO

WTEMP=W/DUMMY

DOT1=SQRT(DDOT(N,GRADX,1,GRADX,1))
DOT2=SQRT(DDOT(N,WTEMP,1,WTEMP,1))
OVERLAP=0.0D0
IF (DOT1*DOT2.NE.0.0D0) THEN
    OVERLAP=DDOT(N,GRADX,1,WTEMP,1)/(DOT1*DOT2)
ENDIF
IF (OVERLAP.GT.0.0D0) THEN
    IF (DEBUG) WRITE(LFH,'(A)') 'Search direction has positive projection onto gradient - reversing step'
    WSS(1:N,1+POINT)=-W  !!! DJW, reverses step
ENDIF
! }}}

W=GRADX
SLENGTH=SQRT(SUM(WSS(:,1+POINT)**2))
IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH

! We now have the proposed step; save X here so that we can undo the step reliably.

XSAVE=X; X=X+STP*WSS(:,1+POINT)

CALL POTENTIAL(X,ENEW,GNEW,RMS,.TRUE.,.FALSE.)
write(*,*) 'ENEW= ',ENEW

IF ((ENEW-ENERGY.LE.MAXERISE).AND.(ENEW-ENERGY.GT.MAXEFALL)) THEN
	! {{{
	ITER=ITER+1
	ITDONE=ITDONE+1
	ENERGY=ENEW
	GRADX=GNEW
	
	IF (DEBUG) WRITE(LFH,103) ' Energy and RMS force=', ENERGY,RMS, &
	' after ',ITDONE,' LBFGS steps, ', &
	' step:',STP*SLENGTH
	
	103   FORMAT(A,F20.10,G20.10,A,I6,A,A,F13.10)
	
	!  May want to prevent the PE from falling too much if we are trying to visit all the
	!  PE bins. Halve the step size until the energy decrease is in range.
	!
	! }}}
ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
	! {{{
	!
	!  Energy decreased too much - try again with a smaller step size
	!
	IF (NDECREASE.GT.5) THEN
	! {{{
	NFAIL=NFAIL+1
	WRITE(LFH,'(A,A,G20.10)')   ' in mylbfgs LBFGS step cannot find an energy in the required range, ',&
	'    NFAIL=',NFAIL
	!
	! Resetting to XSAVE should be the same as subtracting the step. 
	! If we have tried PROJI with Thomson then the projection is non-linear
	! and we need to reset to XSAVE. This should always be reliable!
	!
	X=XSAVE
	GRADX=GNEW ! GRADX contains the gradient at the lowest energy point
	
	ITER=0   !  try resetting
	IF (NFAIL.GT.20) THEN
	WRITE(LFH,'(A)') ' Too many failures - giving up '
	RETURN
	ENDIF
	GOTO 30
	! }}}
	ENDIF
	!
	! Resetting to XSAVE and adding half the step should be the same as subtracting 
	! half the step. 
	!
	X=XSAVE
	X=X+0.5*STP*WSS(:,1+POINT)
	STP=STP/2.0D0
	NDECREASE=NDECREASE+1
	IF (DEBUG) WRITE(LFH,105) &
	' energy decreased too much from ',ENERGY,&
	' to ',ENEW,&
	' decreasing step to ', STP*SLENGTH
	
	105 FORMAT(A,F19.10,A,F16.10,A,F15.8) 
	
	!GOTO 20
	! }}}
ELSE
	! {{{
	!
	!  Energy increased - try again with a smaller step size
	!
	IF (NDECREASE.GT.10) THEN ! DJW
	! {{{
		NFAIL=NFAIL+1
		WRITE(LFH,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL
		!
		! Resetting to XSAVE should be the same as subtracting the step. 
		! If we have tried PROJI with Thomson then the projection is non-linear
		! and we need to reset to XSAVE. This should always be reliable!
		!
		X=XSAVE
		GRADX=GNEW ! GRAD contains the gradient at the lowest energy point
		ITER=0   !  try resetting
		IF (NFAIL.GT.5) THEN         
			WRITE(LFH,'(A)') ' Too many failures - giving up '
			RETURN
		ENDIF
		GOTO 30
		! }}}
	ENDIF
	!
	! Resetting to XSAVE and adding 0.1 of the step should be the same as subtracting 
	! 0.9 of the step. 
	! If we have tried PROJI with Thomson then the projection is non-linear
	! and we need to reset to XSAVE. This should always be reliable!
	!
	X=XSAVE
	X=X+0.1D0*STP*SUM(WSS(:,1+POINT))
	! }}}
ENDIF

STP=STP/1.0D1
NDECREASE=NDECREASE+1
IF (DEBUG) WRITE(LFH,104) ' energy increased from ',ENERGY,&
    ' to ',ENEW,&
    ' decreasing step to ',STP*SLENGTH
104 FORMAT(A,F20.10,A,F20.10,A,F20.10) 

!
!     Compute the new step and gradient change

30 WSS(:,1+POINT)=STP*WSS(:,1+POINT)     ! save the step taken
WDG(:,1+POINT)=GRADX-W                 ! save gradient difference: W contains the old gradient
POINT=POINT+1
IF (POINT.EQ.M) POINT=0

DEALLOCATE(RHO,ALPHA,WSS,WDG)

GOTO 10           ! Go to the termination test
! }}}
RETURN

END SUBROUTINE
! }}}

SUBROUTINE ECHO_S
WRITE(*,'(A)') "**********************************************************************************************" 
ENDSUBROUTINE ECHO_S

SUBROUTINE PRINTHELP
! {{{
write(*,*) '======================='
write(*,*) 'gmi - A program for finding global minima'
write(*,*) ''
write(*,*) '-h      display help' 
write(*,*) '-g      turn on debugging info' 
write(*,*) ''
write(*,*) '======================='
! }}}
END SUBROUTINE PRINTHELP

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

SUBROUTINE TAKESTEP
      ! {{{

      IMPLICIT NONE

      INTEGER IA
      DOUBLE PRECISION ::  RND(3)

      DO IA=1,NATOMS
         CALL GETRND(RND,3,-1.0D0,1.0D0)
         SCREENC(IA,1:3)=SCREENC(IA,1:3)+STEP*RND(1:3)
      ENDDO
      
      RETURN
      ! }}}
END SUBROUTINE TAKESTEP

SUBROUTINE CENTRE(R)
! {{{

IMPLICIT NONE

! subroutine parameters 
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:,:) :: R
! local parameters 
INTEGER I,K

RMASS=SUM(R,DIM=1)/SIZE(R,DIM=1)

do K=1,3
        R(:,K)=R(:,K)-RMASS(K)
ENDDO

IF (DEBUG) WRITE(LFH,'(A,3G20.10)') 'centre> centre of mass reset to the origin from ',RMASS
! }}}
END SUBROUTINE CENTRE

ENDMODULE FUNC
