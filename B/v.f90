      MODULE V

      IMPLICIT NONE
      SAVE
        
      LOGICAL ::  MYBLNT, RMST, USEKW
      LOGICAL ::    DEB_BLN=.FALSE.
      DOUBLE PRECISION ::   dE_fz
      INTEGER PATOM1,PATOM2
      INTEGER LFH, IFH, EA_FH, ENERGY_FH, MARKOV_FH, BEST_FH
      INTEGER PAIRDIST_FH, COORDS_FH, DATA_FH, RMSD_FH
      DOUBLE PRECISION :: PI,ZERO,ONE
      INTEGER ::      SLEN

      PARAMETER (PI=3.141592654D0)
      PARAMETER (ZERO=1.0D-100)
      PARAMETER (ONE=1.0D0)
      PARAMETER (SLEN=100)

      !LOGICAL ::  GUIDECHANGET,  CSMDOGUIDET, GUIDET
      !COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET

      ! files {{{
      CHARACTER(LEN=130) C_FILE,D_FILE,O_FILE,LE_FILE,E_FILE,SEED_FILE,EA_FILE
      CHARACTER(LEN=130) BLNTYPE

      ! list of energies, EA(1) is the total one 
      DOUBLE PRECISION, DIMENSION(10) :: EA

      INTEGER NR
      INTEGER NRG
      ! }}}

      DOUBLE PRECISION, dimension(:), ALLOCATABLE :: MSCREENC

      parameter(C_FILE="coords",D_FILE="data")
      parameter(LE_FILE="lowest",E_FILE="e.tex")
      parameter(O_FILE="out",SEED_FILE="seed")
      parameter(EA_FILE="ea")

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
      DOUBLE PRECISION, ALLOCATABLE :: QMINP(:,:), INTEQMINP(:,:)

      INTEGER :: MXATMS=0

      ENDMODULE V
