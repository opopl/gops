
        IMPLICIT NONE
        
        ! subroutine parameters  {{{

        ! N - number of particles
        INTEGER, INTENT(IN) :: N
        ! PTYPE: type of BLN potential {{{
        !
        !       GO   Go-like
        !       WT   Wild-type
        !
        CHARACTER(LEN=*), INTENT(IN) :: PTYPE
        ! }}}

        ! calculate the gradient if GRADT=.TRUE.
        ! calculate the Hessian if HESST=.TRUE. 
        LOGICAL, INTENT(IN) :: GRADT, HESST
        ! X - input coordinates
        ! GRADX - output gradient
        !DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: X
        !DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: GRADX
        DOUBLE PRECISION, DIMENSION(3*N), INTENT(IN) :: X
        DOUBLE PRECISION, DIMENSION(3*N), INTENT(OUT) :: GRADX

        ! Output Energies are specified by array E(:) {{{
        !
        ! E(1) is the total energy.
        !       Other energies:
        !           E(2) - non-bonded
        !           E(3) - bonded
        !           E(4) - bond angles
        !           E(5) - torsional angles
        DOUBLE PRECISION, DIMENSION(:),INTENT(OUT) :: E
        ! }}}

        ! }}}

        ! local parameters {{{
  
        LOGICAL CONNECT(N,N)

        ! R GRAD HESS 
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: R, GRAD
        DOUBLE PRECISION, DIMENSION(:,:) :: HESS

        INTEGER J1,J2

        ! AB, CD - BLN model coefficients {{{
        ! AB(,1)  =>  A_PARAM(N,N)
        ! AB(,2)  =>  B_PARAM(N,N)
        ! CD(,1)  =>  C_PARAM(N)
        ! CD(,2)  =>  D_PARAM(N)
        DOUBLE PRECISION, DIMENSION(N,N,2) :: AB
        DOUBLE PRECISION, DIMENSION(N,2)   :: CD
        ! }}}
        ! R DR LEN_DR - particle positions and distances {{{
        ! particle positions R => (X,Y,Z)
        ! particle relative positions DR_{ij} =>  R_i-R_j 
        ! distances between particles LEN_DR_ij => || DR_ij ||
        DOUBLE PRECISION, DIMENSION(N,N,3) :: DR
        DOUBLE PRECISION, DIMENSION(N,N) :: LEN_DR
        ! }}}
        ! Angles {{{
        ! 1 => bond angles
        ! 2 => torsion (dihedral) angles
        DOUBLE PRECISION, DIMENSION(N,2) :: ANG 
        ! F => d(potential)/d(angle)
        DOUBLE PRECISION, DIMENSION(-1:N+1,2) :: F
        ! ==============================
        ! for torsional angles:
        ! 
        ! FB(:,1) => F*B
        ! FB(:,2) => F/B
        ! 
        ! ==============================
        DOUBLE PRECISION, DIMENSION(N,2) :: FB
        ! AN temporary angle variable
        DOUBLE PRECISION   AN
        ! }}}
        ! cross products {{{
        !
        !> @param XPD_2 - squared cross product 
        !> @param XPD - length of cross product 
        !> @param VXPD - cross product vector
        !> @param HVXPD - cross product direction
        DOUBLE PRECISION, DIMENSION(N-1) :: XPD_2, XPD
        DOUBLE PRECISION, DIMENSION(N-1,3) :: VXPD, HVXPD, PP
        ! }}}
        ! Bond vectors {{{
        ! B(:) =>  bond vectors lengths
        DOUBLE PRECISION, DIMENSION(N) :: B
        ! bond vectors, BVR_i => DR(i,i+1) => R_{i+1}-R_i 
        DOUBLE PRECISION, DIMENSION(N-1,3) :: BVR, EB

        ! DPD(1:N-2) array of dot-products between the bond vectors BVR(i)
        DOUBLE PRECISION, DIMENSION(N-2) :: DPD
        ! }}}
        ! Gradients {{{
        ! G, GNB, GB, GBA, GTA: vectors representing gradients of different kinds
        ! 
        !       G       => total gradient
        !       GNB     => non-bonded
        !       GB      => bonded
        !       GBA     => bond angles
        !       GTA     => torsional angles
        !
        DOUBLE PRECISION GRAD_MIN(N,3), GRAD_PLUS(N,3)
        DOUBLE PRECISION, DIMENSION(N,3) :: GBA, GNB, GTA, GB 
        DOUBLE PRECISION, DIMENSION(N,3) :: GTA_I, GTA_J, GTA_K, GTA_L
        DOUBLE PRECISION, DIMENSION(N,3) :: GBA_I, GBA_J, GBA_K
        DOUBLE PRECISION ::     DF, FRR(3)
        DOUBLE PRECISION, DIMENSION(N,3) :: G
        ! }}}
        INTEGER NTYPE(N), I, J, K, ICOUNT

        DOUBLE PRECISION COS_PHI, COS_THETA

        ! Hessian {{{
        ! Hessian - (N,N) matrix of second-order derivatives

        ! IK,JK,KJ,KI - used for Hessian calculation
        INTEGER ::     IK,JK,KI,KJ,IH,JH,KH
        ! }}}
      
        ! LJREP => repulsion
        ! LJATT => attraction
        !DOUBLE PRECISION, DIMENSION(N,N) :: LJREP, LJATT

        ! Other constants {{{
        ! 1.8326 RADIANS IS 105 DEGREES

        DOUBLE PRECISION RMASS, SIGMA, THETA_0, EPSILON, DELTA, PI4, TWOPI, RK_R, RK_THETA
        DOUBLE PRECISION ::     RVAR
        PARAMETER(TWOPI=6.283185307179586477D0)
        PARAMETER(THETA_0 = 1.8326D0, PI4=0.7853981633974483096D0 )
        PARAMETER (RMASS = 40.0, EPSILON = 0.0100570)
        PARAMETER (SIGMA=3.4, DELTA=1.0D-6)
        PARAMETER (RK_R = 20.0*0.0100570, RK_THETA = 20.0*0.0100570)

        DOUBLE PRECISION, DIMENSION(15) :: RAD, S   

                ! 1 => non-bonded  
                ! 2 => bonded 
                ! 3 => bond angles 
                ! 4 => torsional angles
        ! }}}

        S(1)=SIGMA
        S(6)=S(1)**6 
        S(12)=S(6)**2
       
        ! }}}
