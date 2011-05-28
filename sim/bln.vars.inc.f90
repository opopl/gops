
        IMPLICIT NONE
        
        ! number of particles
        INTEGER, INTENT(IN) :: N

        ! input vector of coordinates
        DOUBLE PRECISION, DIMENSION(3*N) :: QO, FQ
        
        LOGICAL GRADT

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
        !
        DOUBLE PRECISION, DIMENSION(N) :: COSTOR, SINBOND, DFAC
        ! 1 => bond angles
        ! 2 => torsion (dihedral) angles
        DOUBLE PRECISION, DIMENSION(N,2) :: ANG 
        ! 
        ! cross products:
        !
        ! XPD_2 - squared cross product 
        ! XPD - length of cross product 
        ! VXPD - cross product vector
        ! HVXPD - cross product direction
        DOUBLE PRECISION, DIMENSION(N-1) :: XPD_2, XPD
        DOUBLE PRECISION, DIMENSION(N-1,3) :: VXPD, HVXPD        
        ! inverse cross products
        DOUBLE PRECISION, DIMENSION(N) :: IXPD
        ! bond vectors lengths
        DOUBLE PRECISION, DIMENSION(N) :: LEN_BVR
        ! bond vectors, BVR_i => DR(i,i+1) => R_{i+1}-R_i 
        DOUBLE PRECISION, DIMENSION(N-1) :: BVR

        DOUBLE PRECISION, DIMENSION(N-1,3) :: DPD
        DOUBLE PRECISION, DIMENSION(N,3) :: FBA, FNB, FTA, F
        
        ! type of BLN potential
        !
        !       GO   Go-like
        !       WT   Wild-type
        !

        CHARACTER(LEN=*) :: PTYPE
        INTEGER NTYPE(N), I, J, JMAX, K, KMAX, ICOUNT
        DOUBLE PRECISION RK_R, RK_THETA

        DOUBLE PRECISION COS_PHI, COS_THETA, DUMMY, DUMMY2

        ! Hessian - (N,N) matrix of second-order derivatives

        DOUBLE PRECISION, DIMENSION(N,N) :: HESS

        DOUBLE PRECISION E_NBOND, E_BOND, E_BANGLE, E_TANGLE, RAD6
        DOUBLE PRECISION RMASS, SIGMA, EPSILON, DELTA
        DOUBLE PRECISION THETA_0

        DOUBLE PRECISION GRAD(3*N),GRADIENT(3*N)
      
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

        DOUBLE PRECISION  A4, COEF, COEF1, COEF2, COEF3, A3, DEN2, A2, A1, &
                DEN1, RNUM, DEN, RVAR, FRR(3), DF, RAD(15), S(15)

                ! 1 => non-bonded  
                ! 2 => bonded 
                ! 3 => bond angles 
                ! 4 => torsional angles

        DOUBLE PRECISION, DIMENSION(4) :: E
        
        LOGICAL GTEST, STEST
        LOGICAL CONNECT(N,N)

        INTEGER J1,J2

        DOUBLE PRECISION FQ_PLUS(3*N), FQ_MINUS(3*N)
        ! F, FNB, FB, FBA, FTA: vectors used in the gradient calculations
        ! 
        !       F       => total gradient
        !       FNB     => non-bonded
        !       FB      => bonded
        !       FBA     => bond angles
        !       FTA     => torsional angles
        !
        DOUBLE PRECISION, DIMENSION(N,3) :: F, FNB, FB, FBA, FTA
        DOUBLE PRECISION RAD7, RAD14, DF, RVAR, DEN, RNUM, DEN1, A1, A2, DEN2
        DOUBLE PRECISION A3, COEF, COEF1, COEF2, COEF3, A4
        ! old vars
        double precision, dimension(n): fba_x,fba_y,fba_z 

        S(1)=SIGMA
        S(6)=S(1)**6 
        S(12)=S(6)**2
       
