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
!> the energy (ENERGY), gradient (GRAD), and 
!> the matrix of second derivatives are returned.
!
!> @author John Rose
!
! }}}
!------------------------------------------------------------
        SUBROUTINE EBLN(N,QO,GRAD,ENERGY,GTEST,PTYPE)
! {{{
include bln.vars.inc.f90

        STEST=.FALSE.

include bln.param.inc.f90       ! calculate BLN model parameters 
include bln.cic.inc.f90         ! calculate internal coordinates
include bln.ce.inc.f90          ! calculate energies

        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN

include bln.grad.inc.f90        ! calculate gradients

        IF (.NOT.STEST) RETURN

include bln.hess.inc.f90        ! calculate Hessian

        RETURN
        END
! }}}
! }}} 


 !subroutine cosba()
!subroutine cosba(i,j,ii,k,l,rn,nn, mut, m.dml.dm2,ddm) E


        ENDMODULE
