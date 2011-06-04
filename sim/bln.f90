     MODULE BLN 

     IMPLICIT NONE

     CONTAINS

!------------------------------------------------------------
! Doxygen - EBLN {{{
!
!> @name EBLN
!
!> @brief
!> Calculate the energy, gradient, and second
!> derivatives for a given configuration of the 46 particle polymer chain.
!> A configuration and number of particles is passed to the subroutine and
!> the energy (ENERGY), gradient (GRADIENT), and 
!> the matrix of second derivatives (HESS) are returned.
!
! }}}
!------------------------------------------------------------
        SUBROUTINE EBLN(N,R,GRAD,ENERGY,HESS,PTYPE,GRADT,HESST)
! {{{

include "bln.vars.inc.f90"        ! variables
include "bln.param.inc.f90"       ! calculate BLN model parameters 
include "bln.cic.inc.f90"         ! calculate internal coordinates
include "bln.ce.inc.f90"          ! calculate energies
IF (.NOT. GRADT) RETURN
include "bln.grad.inc.f90"        ! calculate gradient G
GRAD=G
IF (.NOT. HESST) RETURN
include "bln.hess.inc.f90"        ! calculate Hessian

        RETURN
        END SUBROUTINE
! }}} 

        ENDMODULE
