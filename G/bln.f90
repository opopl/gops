
     MODULE MODBLN 

     IMPLICIT NONE

     CONTAINS

!------------------------------------------------------------
! Doxygen - EBLN {{{
!
!> @name EBLN
!
!> @brief
!>      Calculate the energy, gradient, and second
!>          derivatives for a given configuration of the 46 particle polymer chain.
!>          A configuration and number of particles is passed to the subroutine and
!>          the energy (ENERGY), gradient (GRADIENT), and 
!>          the matrix of second derivatives (HESS) are returned.
!
!> @param[in]    integer N      number of particles
!> @param[in]    dp R(:,:)      input coordinates
!> @param[out]   dp E(:)        array of calculated energies:
!>                      E(1)    total energy
!>                      E(2)    non-bonded contribution
!>                      E(3)    bonded
!>                      E(4)    bond angles
!>                      E(5)    torsional angles
!> @param[in]    ch(*) PTYPE    model type. Values:
!>                      GO      Go-like 
!>                      WT      Wild-type 
!> @param[in]    logical GRADT  Do we need to calculate the gradient?
!> @param[in]    logical HESST  Do we need to calculate the Hessian?
! }}}
!------------------------------------------------------------
        SUBROUTINE EBLN(N,X,E,GRADX,HESS,PTYPE,GRADT,HESST)
! {{{

include "bln.vars.inc.f90"        ! variables
include "bln.am.inc.f90"          ! allocate memory
include "bln.param.inc.f90"       ! calculate BLN model parameters 
include "bln.cic.inc.f90"         ! calculate internal coordinates
include "bln.ce.inc.f90"          ! calculate energies
include "bln.grad.inc.f90"        ! calculate gradient G
include "bln.hess.inc.f90"        ! calculate Hessian

        RETURN
        END SUBROUTINE
! }}} 

        ENDMODULE MODBLN
