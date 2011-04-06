
        IMPLICIT NONE
        
        ! number of particles
        INTEGER, INTENT(IN) :: N

        ! input vector of coordinates
        DOUBLE PRECISION, DIMENSION(3*N) :: QO,FQ
        
        LOGICAL GRADT

        DOUBLE PRECISION ENERGY

        ! A => array A, B, C, D coefficients for the generic BLN model
        DOUBLE PRECISION, DIMENSION(N,4) :: A

        ! AB(,1)=>A_PARAM(N,N)
        ! AB(,2)=>B_PARAM(N,N)
        ! CD(,1)=>C_PARAM(N)
        ! CD(,2)=>D_PARAM(N)
        DOUBLE PRECISION, DIMENSION(N,N,2) AB
        DOUBLE PRECISION, DIMENSION(N,2)   CD

        DOUBLE PRECISION, DIMENSION(N,3) :: R
        DOUBLE PRECISION, DIMENSION(N,N,3) :: DR
        DOUBLE PRECISION, DIMENSION(N,N) :: RADII
        DOUBLE PRECISION, DIMENSION(N) :: X_PROD, COSTOR, BOND_ANGLE, TOR_ANGLE, SINBOND, DFAC
        ! 1=> bond angles
        ! 2=> torsion (dihedral) angles
        DOUBLE PRECISION, DIMENSION(N,2) :: ANG 

        DOUBLE PRECISION, DIMENSION(N,3) :: DOT_PROD, FBA, FNB, FTA, F
        
        CHARACTER(LEN=*) :: PTYPE
        INTEGER NTYPE(N), I, J, JMAX, K, ICOUNT
        DOUBLE PRECISION RK_R, RK_THETA

        DOUBLE PRECISION COS_PHI, COS_THETA, DUMMY, DUMMY2

        DOUBLE PRECISION E_NBOND, E_BOND, E_BANGLE, E_TANGLE, RAD6
        DOUBLE PRECISION RMASS, SIGMA, EPSILON, DELTA
        DOUBLE PRECISION THETA_0

        DOUBLE PRECISION GRAD(3*N)
      
        ! LJREP => repulsion
        ! LJATT => attraction
        DOUBLE PRECISION, DIMENSION(N,N) :: LJREP,LJATT

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

        DOUBLE PRECISION FQ1(3*N), FQ2(3*N)
        DOUBLE PRECISION, DIMENSION(N,3) :: F, FNB, FB, FBA, FTA
      DOUBLE PRECISION RAD7, RAD14, DF, FXX, FZZ, FYY, RVAR, DEN, RNUM, DEN1, A1, A2, DEN2
      DOUBLE PRECISION A3, COEF, COEF1, COEF2, COEF3, A4

