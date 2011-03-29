MODULE SDWATER

IMPLICIT NONE

REAL*8 :: ALPHA = 1.444 !ANGSTROMS**3 !POLARIZABILITY OF O. (H+ CAN'T BE POLARIZED AT ALL)
REAL*8 :: RE = 0.9584 !ANGSTROMS
INTEGER :: NO, NA
INTEGER :: LWORK
REAL*8, DIMENSION(:), ALLOCATABLE :: WORK, Q

INTERFACE DIST
   MODULE PROCEDURE DIST_X, DIST_D
END INTERFACE

CONTAINS

!#######################################################################

!VECTOR BETWEEN ITH PARTICLE AND JTH PARTICLE
PURE FUNCTION DIFF(X, J, I)
REAL*8, DIMENSION(3) :: DIFF
REAL*8,  INTENT(IN), DIMENSION(:) :: X
INTEGER, INTENT(IN)               :: I, J
DIFF = X(I*3-2:I*3) - X(J*3-2:J*3)
END FUNCTION DIFF

!DISTANCE BETWEEN ITH PARTICLE AND JTH PARTICLE IN ANGSTROM
PURE REAL*8 FUNCTION DIST_X(X, I, J)
REAL*8,  INTENT(IN), DIMENSION(:) :: X
INTEGER, INTENT(IN)               :: I, J
DIST_X = SQRT(SUM(DIFF(X, I, J)**2))
END FUNCTION DIST_X

PURE REAL*8 FUNCTION DIST_D(DIFF)
REAL*8,  INTENT(IN), DIMENSION(3) :: DIFF
DIST_D = SQRT(SUM(DIFF**2))
END FUNCTION DIST_D

!#######################################################################
! FUNCTIONS
!#######################################################################

PURE REAL*8 FUNCTION PHIHH(R)
REAL*8, INTENT(IN) :: R ! ANGSTROM
PHIHH = 332.1669 / R
END FUNCTION PHIHH

!#######################################################################

PURE REAL*8 FUNCTION PHIOH(R)
REAL*8, INTENT(IN) :: R
PHIOH = 332.1669 / R * (10*EXP(-3.699392820*R)-2) &
        + (-184.6966743*(R-RE) + 123.9762188*(R-RE)**2) * EXP(-8*(R-RE)**2)
END FUNCTION PHIOH

!#######################################################################

PURE REAL*8 FUNCTION PHIOO(R)
REAL*8, INTENT(IN) :: R
PHIOO = 1328.6676/R               &
        + 24/(1+EXP(2.5*(R-2.9))) &
        + 90/(1+EXP(8*(R-2.45)))  &
        + EXP(-6*(R-2.7))
END FUNCTION PHIOO

!#######################################################################

PURE REAL*8 FUNCTION OMKR3(R) ! (1 - K) / R**3
REAL*8, INTENT(IN) :: R
OMKR3 = 1 / ( R**3                                        &
               + 1.855785223*(R-RE)**2 * EXP(-8*(R-RE)**2) &
               + 16.95145727 * EXP(-2.702563425*R) )
END FUNCTION OMKR3

!#######################################################################

PURE REAL*8 FUNCTION OMLR3(R) ! (1 - L) / R**3
REAL*8, INTENT(IN) :: R
OMLR3 = 1 - EXP(-3.169888166*R) * (1 + (3.169888166 + (5.024095492 + (-17.99599078 + 23.92285*R)*R)*R)*R )
OMLR3 = OMLR3 / R**3
END FUNCTION OMLR3

!#######################################################################
! FUNCTION DERIVATIVES
!#######################################################################

PURE REAL*8 FUNCTION DPHIHHDROR(R) ! MINUS DERIVATIVE OF PHIHH W.R.T R DIVIDED BY R
REAL*8, INTENT(IN) :: R ! ANGSTROM
DPHIHHDROR = 332.1669 / R**3
END FUNCTION DPHIHHDROR

!#######################################################################

PURE REAL*8 FUNCTION DPHIOHDROR(R) ! MINUS DERIVATIVE OF PHIOH W.R.T R DIVIDED BY R
REAL*8, INTENT(IN) :: R
REAL*8 :: TMP
TMP = 3.699392820*R
DPHIOHDROR = 332.1669 / R**2 * (10*EXP(-TMP)*(1+TMP) - 2)
TMP = R - RE
DPHIOHDROR = DPHIOHDROR + EXP(-8*TMP**2) * (184.6966743 + TMP*(-2*123.9762188 + 16*TMP*(-184.6966743 + 123.9762188*TMP)))
DPHIOHDROR = DPHIOHDROR / R
END FUNCTION DPHIOHDROR

!#######################################################################

PURE REAL*8 FUNCTION DPHIOODROR(R) ! MINUS DERIVATIVE OF PHIOO W.R.T R DIVIDED BY R
REAL*8, INTENT(IN) :: R
REAL*8 :: TMP
DPHIOODROR = 1328.6676 / R**2
TMP = EXP(2.5*(R-2.9))
DPHIOODROR = DPHIOODROR + 24 / (1+TMP)**2 * 2.5*TMP
TMP = EXP(8*(R-2.45))
DPHIOODROR = DPHIOODROR + 90 / (1+TMP)**2 * 8 * TMP
DPHIOODROR = DPHIOODROR + 6 * EXP(-6*(R-2.7))
DPHIOODROR = DPHIOODROR / R
END FUNCTION DPHIOODROR

!#######################################################################

PURE REAL*8 FUNCTION DOMKR3DR(R) ! D/DR (1 - K)/R**3
REAL*8, INTENT(IN) :: R
REAL*8 :: RRE2
RRE2 = (R - RE)**2
DOMKR3DR = R**3 + 1.855785223*RRE2 * EXP(-8*RRE2) + 16.95145727 * EXP(-2.702563425*R)
DOMKR3DR = - ( 3*R**2 + 1.855785223*(R-RE)*(2 - 16.95145727*RRE2)*EXP(-8*RRE2) - 2.702563425 * 16.95145727 * EXP(-2.702563425*R) ) &
  &      / DOMKR3DR**2
END FUNCTION DOMKR3DR

!#######################################################################

PURE REAL*8 FUNCTION DOMLR3DR(R) ! D/DR (1 - L)/R**3
REAL*8, INTENT(IN) :: R
DOMLR3DR = R**2
DOMLR3DR = - (3 - EXP(-3.169888166*R) * (3 + (2*3.169888166 + (5.024095492 - 23.92285*DOMLR3DR)*R)*R +  &
  &           3.169888166*(1 + (3.169888166 + (5.024095492 + (-17.99599078 + 23.92285*R)*R)*R)*R)*R ) ) / DOMLR3DR**2
END FUNCTION DOMLR3DR

!#######################################################################
! DIPOLE
!#######################################################################

! CALCULATES DIPOLE MU(I,:) ON THE ITH OXYGEN ATOM

! MU = ALPHA * G
! WHERE ALPHA IS POLARIZABILITY AND G IS MODIFIED FIELD
! G = - SUM( |R> * OMK * Q / R**3) - SUM( <T|MU> * OMK / R**3)
! T = 1 - 3 |R><R| / R**2
! PHI = 0.5 * SUM( <MU|R> * OML * Q / R**3)

! LINEAR COUPLED EQUATIONS: MU=ALPHA*G, G=B+C.MU
! G = B + ALPHA*C.G
! I.G = B + D.G
! (I-D).G = B
! A.G = B !SOLVE THIS WITH ?GETRS

!######################################################################

FUNCTION FIELD1(X) ! AT OXYGENS = - SUM( |R> * OMK * Q / R**3 )
REAL*8, INTENT(IN) :: X(3*NA)
REAL*8 :: FIELD1(3*NO)
REAL*8 :: RV(3), R
INTEGER :: I, J
FIELD1 = 0
DO I = 1, NO !OXYGEN ATOMS
	DO J = 1, NA
      IF (J==I) CYCLE
      RV = DIFF(X,I,J)
      R = DIST(RV)
		FIELD1(3*I-2:3*I) = FIELD1(3*I-2:3*I) - RV*Q(J)*OMKR3(R)
   END DO
END DO
END FUNCTION FIELD1

!######################################################################

FUNCTION FIELD2(X) ! AT OXYGENS = - SUM( <T| * OMK * Q / R**3 )
REAL*8, INTENT(IN) :: X(3*NA)
REAL*8 :: FIELD2(3*NO,3*NO)
REAL*8 :: RV(3), R
INTEGER :: I, J
FIELD2 = 0
DO I = 1, NO !OXYGEN ATOMS
	DO J = 1, NO
      IF (J==I) CYCLE
		RV = DIFF(X, I, J)
		R = DIST(RV)
      FIELD2(3*I-2:3*I,3*J-2:3*J) = FIELD2(3*I-2:3*I,3*J-2:3*J) - T(R, RV)*OMKR3(R)
   END DO
END DO
END FUNCTION FIELD2

!######################################################################

!FIELD FROM AN ELECTRIC DIPOLE
!POT(R) = MU.\HAT{R} / R**2 / 4PIE_0
!WHERE \HAT{R} IS UNIT VECTOR IN DIRECTION OF R
!FIELD = (3*(MU.\HAT{R})*\HAT{R} - MU)/R**3 / 4PIE_0
!T_IJ = \DELTA_IJ - 3*R_I*R_J
PURE FUNCTION T(R, RV)
REAL*8, DIMENSION(3,3) :: T
REAL*8, INTENT(IN) :: R, RV(3)
REAL*8 :: TMP
INTEGER :: K, L
TMP = - 3 / R**2
DO K = 1, 3
   DO L = K, 3
      T(K,L) = TMP * RV(K) * RV(L)
   END DO
END DO
DO K = 2, 3
   DO L = 1, K-1
      T(K,L) = T(L,K)
   END DO
END DO
FORALL(K=1:3) T(K,K) = T(K,K) + 1
END FUNCTION T

!#######################################################################

SUBROUTINE FACTORIZE(A, IPIV)
REAL*8, INTENT(INOUT) :: A(3*NO,3*NO)
INTEGER, INTENT(OUT) :: IPIV(3*NO)
INTEGER :: INFO

CALL DSYTRF('U', 3*NO, A, 3*NO, IPIV, WORK, LWORK, INFO)

IF (INT(WORK(1)) /= LWORK) PRINT *, WORK(1), "IS BETTER THAN", LWORK
IF (INFO < 0) THEN
   PRINT '("THE ", I0, "TH PARAMETER HAD AN ILLEGAL VALUE")', -INFO
   STOP
ELSE IF (INFO > 0) THEN
   PRINT '("ERROR: U(",I0,",",I0,") = 0")', INFO, INFO
   STOP
END IF
END SUBROUTINE FACTORIZE

!#######################################################################

SUBROUTINE SOLVE(A, IPIV, MU)
REAL*8, INTENT(IN) :: A(3*NO,3*NO)
INTEGER, INTENT(IN) :: IPIV(3*NO)
REAL*8, INTENT(INOUT) :: MU(NO,3)
INTEGER :: INFO

CALL DSYTRS('U', 3*NO, 1, A, 3*NO, IPIV, MU, 3*NO, INFO)
IF (INFO < 0) THEN
   PRINT '("THE ", I0, "TH PARAMETER HAD AN ILLEGAL VALUE")', -INFO
   STOP
ELSE IF (INFO > 0) THEN
   PRINT '("ERROR: U(",I0,",",I0,") = 0")', INFO, INFO
   STOP
END IF
END SUBROUTINE SOLVE

!#######################################################################
! DIPOLE DERIVATIVES
!######################################################################

FUNCTION DTDX(R, RV, NDIM, S)
REAL*8, DIMENSION(3,3) :: DTDX
REAL*8, INTENT(IN) :: R, RV(3)
INTEGER, INTENT(IN) :: NDIM, S
REAL*8 :: TMP
INTEGER :: K, L
FORALL(K=1:3,L=1:3) DTDX(K,L) = 2 * RV(K) * RV(L) * RV(NDIM)
TMP = R**2
FORALL(L=1:3) DTDX(NDIM,L) = DTDX(NDIM,L) - TMP*RV(L)
FORALL(K=1:3) DTDX(K,NDIM) = DTDX(K,NDIM) - TMP*RV(K)
DTDX = DTDX * S * 3 / TMP**2
END FUNCTION DTDX

!#######################################################################

FUNCTION DBDX(X, N, MU)
REAL*8, INTENT(IN) :: X(3*NA), MU(3*NO)
INTEGER, INTENT(IN) :: N
REAL*8 :: DBDX(3*NO)
INTEGER :: I, J
REAL*8 :: R, RV(3), TMP
INTEGER :: NAT, NDIM
NAT = (N-1)/3 + 1
NDIM = N - 3*(NAT-1)

IF (NAT > NO) THEN ! .: NAT != I
	DO I = 1, NO
		RV = DIFF(X, I, NAT)
		R = DIST(RV)
		DBDX(3*I-2:3*I) = Q(NAT) * DOMKR3DR(R) * RV(NDIM) / R * RV
		DBDX(3*(I-1)+NDIM) = DBDX(3*(I-1)+NDIM) + OMKR3(R) * Q(NAT)
	END DO
ELSE ! NAT <= NO
	! I != NAT
	DO I = 1, NO
		IF (I==NAT) CYCLE
		RV = DIFF(X, I, NAT)
		R = DIST(RV)
		DBDX(3*I-2:3*I) = DOMKR3DR(R)*RV(NDIM)/R * (RV*Q(NAT) + MATMUL(T(R,RV),MU(3*NAT-2:3*NAT)))

		TMP = OMKR3(R)
		DBDX(3*(I-1)+NDIM) = DBDX(3*(I-1)+NDIM) + TMP * Q(NAT)
		DBDX(3*I-2:3*I) = DBDX(3*I-2:3*I) + TMP * MATMUL(DTDX(R,RV,NDIM,1), MU(3*NAT-2:3*NAT))
	END DO

	! I = NAT
	DBDX(3*NAT-2:3*NAT) = 0
	DO J = 1, NA
		IF (J==NAT) CYCLE
		RV = DIFF(X, NAT, J)
		R = DIST(RV)
                DBDX(3*NAT-2:3*NAT) = DBDX(3*NAT-2:3*NAT) - DOMKR3DR(R) * RV(NDIM)/R * Q(J) * RV
		DBDX(N) = DBDX(N) - OMKR3(R) * Q(J)
	END DO

	DO J = 1, NO
		IF (J==NAT) CYCLE
		RV = DIFF(X, NAT, J)
		R = DIST(RV)
		DBDX(3*NAT-2:3*NAT) = DBDX(3*NAT-2:3*NAT) - DOMKR3DR(R) * RV(NDIM)/R * MATMUL(T(R,RV),MU(3*J-2:3*J))
		DBDX(3*NAT-2:3*NAT) = DBDX(3*NAT-2:3*NAT) + OMKR3(R) * MATMUL(DTDX(R,RV,NDIM,-1),MU(3*J-2:3*J))
	END DO
END IF
DBDX = - ALPHA * DBDX
END FUNCTION DBDX

!#######################################################################
! MAIN SUBROUTINES
!#######################################################################

SUBROUTINE SDINIT(NO_, NP_)
INTEGER, INTENT(IN) :: NO_, NP_
INTEGER :: INFO
REAL*8, DIMENSION(3*NO_, 3*NO_) :: A
REAL*8, DIMENSION(3*NO_) :: B, IPIV

NO = NO_
NA = 3*NO_ + NP_
ALLOCATE(Q(NA))
Q(1:NO) = -2 * SQRT(332.1669)
Q(NO+1:) = SQRT(332.1669)

ALLOCATE(WORK(1))
CALL DSYSV('U', SIZE(A,DIM=2), 1, A, SIZE(A,DIM=1), IPIV, B, SIZE(B), WORK, -1, INFO)
IF (INFO < 0) THEN
   PRINT '("THE ", I0, "TH PARAMETER HAD AN ILLEGAL VALUE")', -INFO
   STOP
ELSE IF (INFO > 0) THEN
   PRINT '("ERROR: U(",I0,",",I0,") = 0")', INFO, INFO
   STOP
END IF
LWORK = WORK(1)
!PRINT *, 'BEST LWORK =', LWORK
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
END SUBROUTINE SDINIT

!#######################################################################

SUBROUTINE DESTRUCT()
DEALLOCATE(WORK, Q)
END SUBROUTINE DESTRUCT

!#######################################################################

REAL*8 FUNCTION SDPOTENTIAL(X) RESULT(V)
REAL*8, INTENT(IN) :: X(:)
INTEGER :: I, J
REAL*8 :: MU(3*NO), RV(3), R, PHI2
REAL*8, DIMENSION(3*NO,3*NO) :: A
INTEGER :: IPIV(3*NO)

!SUM OVER HH PAIRS
V = 0
DO J = NO+2, NA
   DO I = NO+1, J-1
      R = DIST(X,I,J)
      V = V + PHIHH(R)
   END DO
END DO

!SUM OVER OH PAIRS
DO J = 1, NO !O ATOMS
   DO I = NO+1, NA !H ATOMS
      R = DIST(X,I,J)
      V = V + PHIOH(R)
   END DO
END DO 

!SUM OVER OO PAIRS
DO J = 2, NO
   DO I = 1, J-1
      R = DIST(X, I, J)
      V = V + PHIOO(R)
   END DO
END DO

MU = FIELD1(X) * ALPHA

A = FIELD2(X)
A = A * ALPHA ! IS D
A = - A
FORALL(I=1:3*NO) A(I,I) = 1 + A(I,I)

CALL FACTORIZE(A, IPIV)
CALL SOLVE(A, IPIV, MU)

PHI2 = 0
DO I = 1, NA !ALL ATOMS
   DO J = 1, NO !OXYGEN ATOMS
      IF (I==J) CYCLE
      RV = DIFF(X,J,I)
      R = DIST(RV)
      PHI2 = PHI2 + DOT_PRODUCT(MU(3*J-2:3*J), RV) * Q(I) * OMLR3(R)
   END DO
END DO
PHI2 = PHI2 / 2
V = V + PHI2

END FUNCTION SDPOTENTIAL

!######################################################################

FUNCTION SDGRAD(X) RESULT(G)
REAL*8, INTENT(IN) :: X(:)
REAL*8 :: G(SIZE(X))
INTEGER :: I, J, NAT, NDIM, N
REAL*8 :: MU(3*NO), RV(3), R, DMUDX(3*NO)
REAL*8, DIMENSION(3*NO,3*NO) :: A
INTEGER :: IPIV(3*NO)

G = 0

! HH PAIRS
DO J = NO+2, NA
   DO I = NO+1, J-1
		RV = DIFF(X,I,J)
		R = DIST(RV)
		RV = RV * DPHIHHDROR(R)
		G(3*J-2:3*J) = G(3*J-2:3*J) - RV
		G(3*I-2:3*I) = G(3*I-2:3*I) + RV
   END DO
END DO

! OH PAIRS
DO J = 1, NO !O ATOMS
   DO I = NO+1, NA !H ATOMS
  		RV = DIFF(X,I,J)
		R = DIST(RV)
		RV = RV * DPHIOHDROR(R)
		G(3*J-2:3*J) = G(3*J-2:3*J) - RV
		G(3*I-2:3*I) = G(3*I-2:3*I) + RV
   END DO
END DO 

! OO PAIRS
DO J = 2, NO
   DO I = 1, J-1
  		RV = DIFF(X,I,J)
		R = DIST(RV)
		RV = RV * DPHIOODROR(R)
		G(3*J-2:3*J) = G(3*J-2:3*J) - RV
		G(3*I-2:3*I) = G(3*I-2:3*I) + RV
   END DO
END DO

MU = FIELD1(X) * ALPHA

A = FIELD2(X)
A = A * ALPHA ! IS D
A = - A
FORALL(I=1:3*NO) A(I,I) = 1 + A(I,I)

CALL FACTORIZE(A, IPIV)
CALL SOLVE(A, IPIV, MU)

DO N = 1, 3*NA
	NAT = (N-1)/3 + 1
	NDIM = N - 3*(NAT-1)
	DMUDX = DBDX(X, N, MU)
	CALL SOLVE(A, IPIV, DMUDX)
DO I = 1, NA
		DO J = 1, NO
			IF (I==J) CYCLE
			RV = DIFF(X,J,I)
			R = DIST(RV)
			G(N) = G(N) + Q(I) * OMLR3(R) * DOT_PRODUCT(DMUDX(3*J-2:3*J),RV) / 2
		END DO
	END DO

	DO J = 1, NO
		IF (J==NAT) CYCLE
		RV = DIFF(X,J,NAT)
		R = DIST(RV)
		G(N) = G(N) + Q(NAT) / 2 * (MU(3*(J-1)+NDIM) * OMLR3(R) + DOT_PRODUCT(MU(3*J-2:3*J),RV) * DOMLR3DR(R) * RV(NDIM)/R)
	END DO

	IF (NAT <= NO) THEN
		DO I = 1, NA
			IF (I==NAT) CYCLE
			RV = DIFF(X,NAT,I)
			R = DIST(RV)
			G(N) = G(N) - Q(I) / 2 * (MU(N) * OMLR3(R) + DOMLR3DR(R) * RV(NDIM)/R * DOT_PRODUCT(MU(3*NAT-2:3*NAT),RV))
		END DO
	END IF
END DO
END FUNCTION SDGRAD

FUNCTION SDHESS(X, G) RESULT(HESS)
REAL*8, INTENT(IN) :: X(:), G(:)
REAL*8 :: HESS(SIZE(X),SIZE(X)), X0(SIZE(X))
INTEGER :: I, J
REAL*8 :: EPS = 1.0D-8, TMP

X0 = X
DO I = 1, 3*NA
	X0(I) = X(I) + EPS
	HESS(I,:) = (SDGRAD(X0) - G) / EPS
	X0(I) = X(I)
END DO
! SYMMETRIZE
DO I = 1, 3*NA
	DO J = 1, 3*NA
		TMP = 0.5 * (HESS(I,J) + HESS(J,I))
		HESS(I,J) = TMP
		HESS(J,I) = TMP
	END DO
END DO
END FUNCTION SDHESS

!######################################################################

END MODULE SDWATER
