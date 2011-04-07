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

        CALL PARAM_ARRAY(N,AB,CD,PTYPE)
        CALL CALC_INT_COORDS(N,QO,R,DR,BVR,LEN_DR,LEN_BVR,DPD,XPD_2,ANG)
        CALL CALC_ENERGY(N,ENERGY,AB,CD,LEN_DR,LEN_BVR,ANG)

        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN

        CALL CALC_GRADIENT(N,GRAD,AB,CD,DR,DPD,XPD,ANG,LEN_DR)

        IF (.NOT.STEST) RETURN

        CALL CALC_HESS(N,QO,AB,CD,R,DR,DPD,XPD,ANG,LEN_DR)

        RETURN
        END
! }}}
!> Calculate the internal coordinates

        SUBROUTINE CALC_INT_COORDS(N,QO,R,DR,BVR,LEN_DR,LEN_BVR,DPD,XPD_2,ANG)
! {{{
include bln.vars.inc.f90

        DO I = 1, N
          J = (I-1)*3
          DO K=1,3
            R(I,K) = QO(J+K)
          ENDDO
        ENDDO

! INTER-PARTICLE DISTANCES

        DO I = 1, N-1
          DO J = I+1, N
            DR(I,J,1:3) = R(J,1:3) - R(I,1:3)
            LEN_DR(I,J) = SQRT(SUM(DR(I,J,1:3)**2))
            LEN_DR(J,I) = LEN_DR(I,J)
          ENDDO
            BVR(I,1:3)=DR(I,I+1,1:3)
        ENDDO

! DOT PRODUCTS BETWEEN BOND VECTORS

      DO I = 1, N-1
        IF ( I .LE. N-3 ) THEN KMAX=3
        IF ( I .EQ. N-2 ) THEN KMAX=2
        IF ( I .LE. N-1 ) THEN KMAX=1
        DO J=1,KMAX
         J=I+K-1
         DPD(I,K) = SUM(BVR(I,1:3)*BVR(J,1:3))
        ENDDO
        LEN_BV(I)=DPD(I,1)
      ENDDO

! Squared cross-products between adjacent bond vectors i and i+1 

        DO I = 1, N-2
           XPD_2(I) = DPD(I,1)*DPD(I+1,1)-DPD(I,2)**2 
        ENDDO

! BOND ANGLES

        DO I = 1, N-2
            COS_THETA=-DPD(I,2)/(LEN_BV(I)*LEN_BV(I+1))
            ANG(I+1,1) = ACOS(COS_THETA)
        ENDDO

! TORSIONAL ANGLES

        DO I = 1, N-3
            COS_PHI = (DPD(I,2)*DPD(I+1,2)-DPD(I,3)*DPD(I+1,1))
            COS_PHI = COS_PHI/SQRT(XPD_2(I)*XPD_2(I+1))
            IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=COS_PHI/ABS(COS_PHI)
            ANG(I+1,2) = ACOS(COS_PHI)
        ENDDO

        RETURN
        END

! }}} 
!> Calculate the energy
        SUBROUTINE CALC_ENERGY(N,ENERGY,AB,CD,LEN_DR,LEN_BVR,ANG)
! {{{
include bln.vars.inc.f90  

       S(6)=SIGMA**6; E=0.0D0 
        
! non-bonded
        DO I = 1, N-2
          DO J = I+2, N
            RAD(6)=LEN_DR(I,J)**6
            E(1) = E(1) + 4.0*((AB(I,J,1)*S(12)/(RAD(12))) + (AB(I,J,2)*S(6)/RAD(6)))
          ENDDO
        ENDDO

! bonded
        DO I = 1, N-1
          E(2) = E(2) + 0.5*RK_R*(LEN_BVR(I)-SIGMA)**2
        ENDDO

! bond angles 
        DO I = 2, N-1
          E(3) = E(3) + 0.5*RK_THETA*(ANG(I,1)-THETA_0)**2
        ENDDO

! torsional angles 
        DO I = 2, N-2
          E(4) = E(4) + CD(I,1)*(1.0 + COS(ANG(I,2)))+CD(I,2)*(1.0+COS(3.0*ANG(I,2)))
        ENDDO

        ENERGY=SUM(E)

        RETURN
        END
! }}}
! Calculate the gradients

        SUBROUTINE CALC_GRADIENT(N,GRAD,AB,CD,DR,DPD,XPD,ANG,LEN_DR)
! {{{
include bln.vars.inc.f90 
include bln.grad.start.inc.f90
include bln.grad.nbond.inc.f90
include bln.grad.bond.inc.f90
include bln.grad.bangle.inc.f90
include bln.grad.tangle.inc.f90
include bln.grad.end.inc.f90

        RETURN
        END
! }}}

! Doxygen - CALC_HESS {{{
!> @name CALC_HESS
!> @brief Calculate the second derivative matrix (two-sided numerical approach)
!
!> @param[in] N        INTEGER         Number of particles
!> @param[in] QO
! }}}

        SUBROUTINE CALC_HESS(N,QO,AB,CD,R,DR,DPD,XPD,ANG,LEN_DR,NTYPE)
include bln.vars.inc.f90
include bln.hess.inc.f90
        RETURN
        END

! Doxygen - PARAM_ARRAY {{{
!> @name PARAM_ARRAY
!>
!> @brief Fill the parameter arrays for the wild-type (WT) BLN model 
!
!> @param[in]   N  INTEGER              Number of particles 
!> @param[out] AB  DP, DIMENSION(N,N,2) Parameters for the L-J interaction between non-bonded particles. 
!> @param[out] CD  DP, DIMENSION(N,2)   Parameters for the dihedral angle potential. 
!
! }}}

        SUBROUTINE PARAM_ARRAY(N,AB,CD,PTYPE)

        include bln.vars.inc.f90
        include bln.ntype.inc.f90

        SELECTCASE(PTYPE)
                CASE("GO")
			include bln.go.connect.inc.f90
			include bln.go.ab.inc.f90
			include bln.go.cd.inc.f90
                CASE("WT")
                        include bln.wt.ab.inc.f90       ! L-J interaction between non-bonded particles.
                        include bln.go.cd.inc.f90       ! Dihedral angle potential
        ENDSELECT

        RETURN
        END

        ENDMODULE
