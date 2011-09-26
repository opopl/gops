      MODULE V

      IMPLICIT NONE
      SAVE
        
      LOGICAL ::  MYBLNT, RMST, USEKW, USESUF
      LOGICAL ::    DEB_BLN=.FALSE.
      DOUBLE PRECISION ::   dE_fz
      INTEGER PATOM1,PATOM2
      INTEGER LFH, IFH, EA_FH, ENERGY_FH, MARKOV_FH, BEST_FH
      INTEGER PAIRDIST_FH, COORDS_FH, DATA_FH, RMSD_FH, LE_FH
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
      CHARACTER(LEN=130) C_FILE,D_FILE,O_FILE,LE_FILE,E_FILE,SEED_FILE,EA_FILE
      CHARACTER(LEN=130) BLNTYPE

      ! suffix to files
      CHARACTER(LEN=30) SUF

      ! list of energies, EA(1) is the total one 
      DOUBLE PRECISION, DIMENSION(10) :: EA

      INTEGER NR
      INTEGER NRG
      ! }}}

      DOUBLE PRECISION, dimension(:), ALLOCATABLE :: MSCREENC


      ! from module PERMU
      DOUBLE PRECISION, ALLOCATABLE :: FIN(:)
      ! from f1com
      integer :: NCOM
      double precision, allocatable :: xicom(:),pcom(:)
      ! from modhess
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: HESS  !  3*MXATMS,3*MXATMS

      ! from qmodule
      integer,allocatable :: FF(:),INTEFF(:) ! NSAVE
      DOUBLE PRECISION, ALLOCATABLE :: QMIN(:), INTEQMIN(:) ! NSAVE
      DOUBLE PRECISION, ALLOCATABLE :: EAMIN(:)! NSAVE
      DOUBLE PRECISION, ALLOCATABLE :: QMINP(:,:), INTEQMINP(:,:)

      INTEGER :: MXATMS=0

      ENDMODULE V
