
        IMPLICIT NONE
        
        ! number of particles
        INTEGER, INTENT(IN) :: N

        ! input vector of coordinates
        !DOUBLE PRECISION, DIMENSION(N,3) :: QO, FQ
        
        ! calculate the gradient if GRADT=.TRUE.
        ! calculate the Hessian if HESST=.TRUE. 
        LOGICAL GRADT, HESST

        DOUBLE PRECISION ENERGY

        ! AB(,1)  =>  A_PARAM(N,N)
        ! AB(,2)  =>  B_PARAM(N,N)
        ! CD(,1)  =>  C_PARAM(N)
        ! CD(,2)  =>  D_PARAM(N)
        DOUBLE PRECISION, DIMENSION(N,N,2) AB
        DOUBLE PRECISION, DIMENSION(N,2)   CD

        ! particle positions R => (X,Y,Z)
        ! particle relative positions DR_{ij} =>  R_i-R_j 
        ! distances between particles LEN_DR_ij => || DR_ij ||
        DOUBLE PRECISION, DIMENSION(N,3) :: R
        DOUBLE PRECISION, DIMENSION(N,N,3) :: DR
        DOUBLE PRECISION, DIMENSION(N,N) :: LEN_DR
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
        ! 
        ! cross products:
        !
        ! XPD_2 - squared cross product 
        ! XPD - length of cross product 
        ! VXPD - cross product vector
        ! HVXPD - cross product direction
        DOUBLE PRECISION, DIMENSION(N-1) :: XPD_2, XPD
        DOUBLE PRECISION, DIMENSION(N-1,3) :: VXPD, HVXPD, pp
        ! B(:) =>  bond vectors lengths
        DOUBLE PRECISION, DIMENSION(N) :: B
        ! bond vectors, BVR_i => DR(i,i+1) => R_{i+1}-R_i 
        DOUBLE PRECISION, DIMENSION(N-1,3) :: BVR, EB

        DOUBLE PRECISION, DIMENSION(N-1) :: DPD
        ! G, GNB, GB, GBA, GTA: vectors representing gradients of different kinds
        ! 
        !       GRAD    => total gradient
        !       GNB     => non-bonded
        !       GB      => bonded
        !       GBA     => bond angles
        !       GTA     => torsional angles
        !
        DOUBLE PRECISION, DIMENSION(N,3) :: GBA, GNB, GTA, GB 
        DOUBLE PRECISION, DIMENSION(N,3) :: GTA_I, GTA_J, GTA_K, GTA_L
        DOUBLE PRECISION, DIMENSION(N,3) :: GBA_I, GBA_J, GBA_K
        DOUBLE PRECISION ::     DF, FRR(3)
        ! type of BLN potential
        !
        !       GO   Go-like
        !       WT   Wild-type
        !

        CHARACTER(LEN=*) :: PTYPE
        INTEGER NTYPE(N), I, J, JMAX, K, KMAX, ICOUNT
        DOUBLE PRECISION RK_R, RK_THETA

        DOUBLE PRECISION COS_PHI, COS_THETA

        ! Hessian - (N,N) matrix of second-order derivatives

        DOUBLE PRECISION, DIMENSION(N,N) :: HESS
        ! IK,JK,KJ,KI - used for Hessian calculation
        DOUBLE PRECISION ::     IK,JK,KI,KJ,IH,JH,KH

        DOUBLE PRECISION RMASS, SIGMA, EPSILON, DELTA, THETA_0

        DOUBLE PRECISION, DIMENSION(N,3) :: GRAD, GRADIENT
      
        ! LJREP => repulsion
        ! LJATT => attraction

        DOUBLE PRECISION, DIMENSION(N,N) :: LJREP, LJATT

        ! 1.8326 RADIANS IS 105 DEGREES

        DOUBLE PRECISION, PARAMETER :: TWOPI=6.283185307179586477D0
        DOUBLE PRECISION, PARAMETER :: THETA_0 = 1.8326D0, PI4=0.7853981633974483096D0 
        PARAMETER (RMASS = 40.0, EPSILON = 0.0100570)
        PARAMETER (SIGMA=3.4, DELTA=1.0D-6, THETA_0 = 1.8326)
        PARAMETER (RK_R = 20.0*0.0100570, RK_THETA = 20.0*0.0100570)

        DOUBLE PRECISION, DIMENSION(15) :: RAD, S   

                ! 1 => non-bonded  
                ! 2 => bonded 
                ! 3 => bond angles 
                ! 4 => torsional angles

        DOUBLE PRECISION, DIMENSION(4) :: E
        
        LOGICAL GTEST, STEST
        LOGICAL CONNECT(N,N)

        INTEGER J1,J2

        DOUBLE PRECISION GRAD_MIN(N,3), GRAD_PLUS(N,3)

        S(1)=SIGMA
        S(6)=S(1)**6 
        S(12)=S(6)**2
        
