
!> @name V
!> @brief Variable/Constants declarations 

MODULE V

IMPLICIT NONE

SAVE

include "v.dox.i.f90"
include "v.const.i.f90"
include "v.lg.i.f90"
include "v.lbfgs.i.f90"
include "v.files.i.f90"
include "v.mc.i.f90"
include "v.pull.i.f90"
include "v.gen.i.f90"
include "v.pairdist.i.f90"

! reading the data file/command-line

INTEGER, PARAMETER :: MAXNARGS=20

INTEGER NARGS
CHARACTER(LEN=SLEN),DIMENSION(MAXNARGS) :: ARGS

INTEGER :: NSAVE=10
INTEGER :: ISTEP
INTEGER :: NCORE

END MODULE V
