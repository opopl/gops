
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
INTEGER :: FH=20
INTEGER :: LFH=FH+1
INTEGER :: ENERGY_FH=FH+2
INTEGER :: MARKOV_FH=FH+3
INTEGER :: BEST_FH=FH+4
INTEGER :: PAIRDIST_FH=FH+5
INTEGER :: COORDS_FH=FH+6

CHARACTER(LEN=100) LFN
PARAMETER(LFH=10,LFN="GMIN.log")
! }}}

