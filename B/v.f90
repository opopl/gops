      MODULE V

      IMPLICIT NONE
      SAVE
        
      LOGICAL ::  MYBLNT, RMST, USEKW,USEPREF, USERCA
      LOGICAL :: TRACKENERGY,TRACKBEST,TRACKMARKOV, TXYZ
      ! number of specific quenches, at which 
      ! lowest minimum geometries are saved 
      INTEGER :: NSQ
      ! 
      INTEGER, DIMENSION(:), ALLOCATABLE :: SQU
      LOGICAL ::    DEB_BLN=.FALSE.
      DOUBLE PRECISION ::   dE_fz
      INTEGER PATOM1,PATOM2
      INTEGER LFH, IFH, EA_FH, FOFH
      INTEGER PAIRDIST_FH, COORDS_FH, DATA_FH, RMSD_FH, LE_FH
      INTEGER ENERGY_FH, BEST_FH, MARKOV_FH, SXYZ_FH
      DOUBLE PRECISION :: PI,ZERO,ONE
      INTEGER ::      SLEN
      INTEGER THE_SEED

      INTEGER MCSTEPS
      DOUBLE PRECISION ::    TFAC

      PARAMETER (PI=3.141592654D0)
      PARAMETER (ZERO=1.0D-100)
      PARAMETER (ONE=1.0D0)
      PARAMETER (SLEN=100)

      ! name for the model system
      CHARACTER(LEN=SLEN) MODEL

      !LOGICAL ::  GUIDECHANGET,  CSMDOGUIDET, GUIDET
      !COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET

      ! files {{{
      CHARACTER(LEN=130) C_FILE,D_FILE,FO_FILE,LO_FILE,LE_FILE,E_FILE,SEED_FILE,EA_FILE
      CHARACTER(LEN=130) BLNTYPE, SXYZ_FILE

      ! suffix to files
      CHARACTER(LEN=30) SUF, PREF

      !> @name CMDLINE 
      !> @brief Contains command-line arguments
      !> @name PROGNAME 
      !> @brief program name, e.g. B 
      CHARACTER(LEN=SLEN) CMDLINE, PROGNAME

      ! list of energies, EA(1) is the total one 
      DOUBLE PRECISION, DIMENSION(10) :: EA
      DOUBLE PRECISION ::   rgyr

      INTEGER NR
      INTEGER NRG
      ! }}}

      DOUBLE PRECISION, dimension(:), ALLOCATABLE :: MSCREENC

      ! from module PERMU
      DOUBLE PRECISION, ALLOCATABLE :: FIN(:)
      ! from f1com
      INTEGER :: NCOM
      DOUBLE PRECISION, ALLOCATABLE :: XICOM(:),PCOM(:)
      ! from modhess
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HESS  !  3*MXATMS,3*MXATMS

      ! from qmodule
      integer,allocatable :: FF(:),INTEFF(:) ! NSAVE
      DOUBLE PRECISION, ALLOCATABLE :: QMIN(:), INTEQMIN(:) ! NSAVE
      DOUBLE PRECISION, ALLOCATABLE :: EAMIN(:,:)! NSAVE
      DOUBLE PRECISION, ALLOCATABLE :: RGMIN(:)! NSAVE, radius of gyration
      DOUBLE PRECISION, ALLOCATABLE :: QMINP(:,:), INTEQMINP(:,:)

      INTEGER :: MXATMS=0

      ENDMODULE V
