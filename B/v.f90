      MODULE V

      IMPLICIT NONE
      SAVE
        
      LOGICAL ::  MYBLNT, RMST, USEKW
      LOGICAL ::    deb_bln=.false.
      DOUBLE PRECISION ::   dE_fz
      !LOGICAL ::  GUIDECHANGET,  CSMDOGUIDET, GUIDET
      !COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET

      ! files {{{
      CHARACTER(LEN=130) C_FILE,D_FILE,O_FILE,LE_FILE,E_FILE,SEED_FILE,EA_FILE
      CHARACTER(LEN=130) BLNTYPE

      parameter(C_FILE="coords",D_FILE="data")
      parameter(LE_FILE="lowest",E_FILE="e.tex")
      parameter(O_FILE="out",SEED_FILE="seed")
      parameter(EA_FILE="ea")

      INTEGER NR
      INTEGER NRG
      ! }}}

      INTEGER,PARAMETER :: DATA_FH=5
      INTEGER,parameter :: EAFH=30
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
      DOUBLE PRECISION, ALLOCATABLE :: QMINP(:,:), INTEQMINP(:,:)

      INTEGER :: MXATMS=0

      ENDMODULE V
