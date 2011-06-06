
! =========== 
! FILES {{{
! ===========
! --- Filehandle variables 
!       All filehandle variables end with _FH or FH
!
!> @param LFH                 (i)     the output log file
!> @param ENERGY_FH           (i)     energy.dat
!> @param MARKOV_FH           (i)     markov.dat
!> @param BEST_FH             (i)     best.dat
!> @param PAIRDIST_FH         (i)     pairdist.dat
!
!
! --- File names
!
!> @param LFN           (s)     name for the output log file 
!
! --- Values
!
! }}}
! =========== 
! File handling {{{
INTEGER, PARAMETER :: FH=20
INTEGER :: LFH
INTEGER :: ENERGY_FH
INTEGER :: MARKOV_FH
INTEGER :: BEST_FH
INTEGER :: PAIRDIST_FH
INTEGER :: COORDS_FH
INTEGER :: DATA_FH
INTEGER :: RMSD_FH

CHARACTER(LEN=100) LFN
PARAMETER(LFN="GMIN.log")
! }}}

