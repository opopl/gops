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

        CALL CALC_HESS(N,QO,HESS,AB,CD)

        RETURN
        END
! }}}
! Doxygen: CALC_INT_COORDS {{{
!
!> @name CALC_INT_COORDS
!>
!> @brief Calculate the internal coordinates
!>
!> @param[in] N         integer                 - number of coordinates
!> @param[in] QO        double precision (3N)   - vector input coordinates
!
!> @param[out]     R    double precision (N,3)    - output vector with Cartesian coordinates
!> @param[out]     DR   DOUBLE PRECISION (N,N,3)  - vector of coordinate differences
!> @param[out]     BVR  DOUBLE PRECISION (N-1,3)  - bond vectors, BVR(i)=R(i)-R(i+1), i=1,N
!
!> @param[out]     LEN_DR    DOUBLE PRECISION(N,N,3)    - lengths of coordinate differences  
!> @param[out]     LEN_BVR   DOUBLE PRECISION(N-1,3)    - lengths of the bond vectors 
!>
!}}}
        SUBROUTINE CALC_INT_COORDS(N,QO,R,DR,BVR,LEN_DR,LEN_BVR,DPD,XPD,ANG)
! {{{
include bln.vars.inc.f90

! QO => R{{{
        DO I = 1, N
          J = (I-1)*3
          DO K=1,3
            R(I,K) = QO(J+K)
          ENDDO
        ENDDO
! }}}
! INTER-PARTICLE DISTANCES: DR, LEN_DR; BOND VECTORS: BVR {{{

        DO I = 1, N-1
          DO J = I+1, N
            DR(I,J,1:3) = R(J,1:3) - R(I,1:3)
            DR(J,I,1:3) = -DR(I,J,1:3)
            LEN_DR(I,J) = SQRT(SUM(DR(I,J,1:3)**2))
            LEN_DR(J,I) = LEN_DR(I,J)
          ENDDO
            BVR(I,1:3)=DR(I,I+1,1:3)
        ENDDO
! }}}
! DOT PRODUCTS BETWEEN BOND VECTORS: DPD, LEN_BVR  {{{

      DO I = 1, N-1
        IF ( I .LE. N-3 ) THEN KMAX=3
        IF ( I .EQ. N-2 ) THEN KMAX=2
        IF ( I .LE. N-1 ) THEN KMAX=1
        DO K=1,KMAX
         J=I+K-1
         DPD(I,K) = SUM(BVR(I,1:3)*BVR(J,1:3))
        ENDDO
        LEN_BVR(I)=SQRT(DPD(I,1))
      ENDDO
! }}}
! Squared cross-products between adjacent bond vectors i and i+1: XPD_2 {{{

        DO I = 1, N-2
           XPD_2(I) = DPD(I,1)*DPD(I+1,1)-DPD(I,2)**2 
           LEN_XPD(I)=SQRT(XPD_2(I))
        ENDDO
! }}}
! BOND ANGLES: ANG(I,1), I=2,...,N-1 {{{

        DO I = 1, N-2
            COS_THETA=-DPD(I,2)/(LEN_BVR(I)*LEN_BVR(I+1))
            ANG(I+1,1) = ACOS(COS_THETA)
        ENDDO
! }}}
! TORSIONAL ANGLES: ANG(I,2), I=2,...,N-2 {{{

        DO I = 1, N-3
            COS_PHI = (DPD(I,2)*DPD(I+1,2)-DPD(I,3)*DPD(I+1,1))
            COS_PHI = COS_PHI/(XPD(I)*XPD(I+1))
            IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=COS_PHI/ABS(COS_PHI)
            ANG(I+1,2) = ACOS(COS_PHI)
        ENDDO
! }}}
        RETURN
        END

! }}} 
!> @brief Calculate the energy
        SUBROUTINE CALC_ENERGY(N,ENERGY,AB,CD,LEN_DR,LEN_BVR,ANG)
! {{{
include bln.vars.inc.f90  

            E=0.0D0 
        
! non-bonded: E(1) {{{

        DO I = 1, N-2
          DO J = I+2, N
            RAD(6)=LEN_DR(I,J)**6; 
            RAD(12)=RAD(6)**2
            E(1) = E(1) + 4.0*AB(I,J,1)*S(12)/RAD(12)
            E(1) = E(1) + 4.0*AB(I,J,2)*S(6)/RAD(6)
          ENDDO
        ENDDO
! }}}
! bonded: E(2) {{{

        DO I = 1, N-1
          E(2) = E(2) + 0.5*RK_R*(LEN_BVR(I)-SIGMA)**2
        ENDDO
! }}}
! bond angles: E(3) {{{

        DO I = 2, N-1
          E(3) = E(3) + 0.5*RK_THETA*(ANG(I,1)-THETA_0)**2
        ENDDO
! }}}
! torsional angles: E(4) {{{

        DO I = 2, N-2
          E(4) = E(4) + CD(I,1)*(1.0 + COS(ANG(I,2)))
          E(4) = E(4) + CD(I,2)*(1.0 + COS(3.0*ANG(I,2)))
        ENDDO
! }}}
        ENERGY=SUM(E)

        RETURN
        END
! }}}
!> @brief Calculate the gradients

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

        SUBROUTINE CALC_HESS(N,QO,HESS,AB,CD)
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
! {{{
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
! }}}
        RETURN
        END

!subroutine cosba()
!subroutine cosba(i,j,ii,k,l,rn,nn, mut, m.dml.dm2,ddm) E


        ENDMODULE
