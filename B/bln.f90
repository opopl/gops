
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
        SUBROUTINE EBLN(N,X,E,GRADX,HESS,PTYPE,GRADT,HESST,FH,DEB)
! vars {{{
!commented  {{{
!include "bln.vars.inc.f90"        ! variables
!include "bln.am.inc.f90"          ! allocate memory
!include "bln.param.inc.f90"       ! calculate BLN model parameters 
!include "bln.cic.inc.f90"         ! calculate internal coordinates
!include "bln.ce.inc.f90"          ! calculate energies
!include "bln.grad.inc.f90"        ! calculate gradient G
!include "bln.hess.inc.f90"        ! calculate Hessian

        IMPLICIT NONE
        ! }}}
        ! subroutine parameters  {{{

        INTEGER,INTENT(IN) :: FH
        LOGICAL,INTENT(IN) :: DEB
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
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: X
        DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: GRADX
        DOUBLE PRECISION, DIMENSION(:,:),INTENT(OUT) :: HESS

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
        ! CONNECT GRAD R HESS J1 J2 {{{ 
        ! N, N
        !LOGICAL, DIMENSION(:,:),ALLOCATABLE :: CONNECT
        LOGICAL, DIMENSION(N,N) :: CONNECT

        ! R GRAD HESS 
        DOUBLE PRECISION, DIMENSION(N,3) :: R, GRAD
        ! 
        INTEGER NCALL
        SAVE NCALL

        INTEGER J1,J2
        ! }}}
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
        DOUBLE PRECISION, DIMENSION(N-1,2) :: ANG 
        ! F => d(potential)/d(angle)
        DOUBLE PRECISION, DIMENSION(N-1,2) :: F
        ! ==============================
        ! for torsional angles:
        ! 
        ! FB(:,1) => F*B
        ! FB(:,2) => F/B
        ! 
        ! ==============================
        DOUBLE PRECISION, DIMENSION(N-2,2) :: FB
        ! AN temporary angle variable
        DOUBLE PRECISION   AN
        ! }}}
        ! cross products {{{
        !
        !> @param XPD_2 - squared cross product 
        !> @param XPD - length of cross product 
        !> @param VXPD - cross product vector
        !> @param HVXPD - cross product direction
        DOUBLE PRECISION, DIMENSION(N) :: XPD_2, XPD
        DOUBLE PRECISION, DIMENSION(N,3) :: VXPD, HVXPD, PP
        ! }}}
        ! Bond vectors {{{
        ! B(:) =>  bond vectors lengths
        DOUBLE PRECISION ,DIMENSION(N-1) :: B
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
        DOUBLE PRECISION, DIMENSION(N,3) :: GBA, GNB, GTA, GB, G
        DOUBLE PRECISION, DIMENSION(N,3) :: GTA_I, GTA_J, GTA_K, GTA_L
        DOUBLE PRECISION, DIMENSION(N,3) :: GBA_I, GBA_J, GBA_K
        DOUBLE PRECISION ::     DF, FRR(3)
        ! }}}
        ! Other {{{
        INTEGER NTYPE(N), I, J, K, ICOUNT

        DOUBLE PRECISION COS_PHI, COS_THETA
        INTEGER NR,NCALLMAX
        DOUBLE PRECISION, dimension(10) :: RMS

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

        ! }}}
       
        ! }}}
        ! }}}
        !intro {{{
        S(1)=SIGMA
        S(6)=S(1)**6 
        S(12)=S(6)**2

        NR=3*N
        include "ncallmax.i.f90"
        include "fmt.i.f90"
        ! }}}

!! am {{{
        !ALLOCATE(R(N,3),GRAD(N,3))
        !r=0.0d0; grad=0.0d0
        !ALLOCATE(CONNECT(N,N))

        !ALLOCATE(AB(N,N,2))
        !!ALLOCATE(CD(2:N-2,2))
        !ALLOCATE(CD(N-2,2))
        !AB=0.0D0; CD=0.0D0

        !ALLOCATE(DR(N-1,N-1,3))
        !ALLOCATE(LEN_DR(N-1,N-1),B(N-1))
        !dr=0.0D0; len_dr=0.0d0

        !ALLOCATE(ANG(N,2))
        !ALLOCATE(F(-1:N+1,2))
        !ALLOCATE(FB(N,2))
        !ALLOCATE(XPD_2(N-1), XPD(N-1))
        !ALLOCATE(VXPD(N-1,3), HVXPD(N-1,3), PP(N-1,3))
        !allocate(BVR(N-1,3), EB(N-1,3))
        !allocate(DPD(N-2))
        !allocate(GBA(N,3))
        !allocate(GB(N,3))
        !allocate(GNB(N,3))
        !allocate(GTA(N,3))
        !allocate(GTA_I(N,3))
        !allocate(GTA_J(N,3))
        !allocate(GTA_K(N,3))
        !allocate(GTA_L(N,3))
        !allocate(GBA_I(N,3))
        !allocate(GBA_J(N,3))
        !allocate(GBA_K(N,3))
        !allocate(G(N,3))

        !! }}}
        ! param (Parameters ) {{{
      include "bln.ntype.inc.f90"         ! specify bead types 

        SELECTCASE(trim(PTYPE))
                CASE("GO")                ! Go-like model
      ! ==================================================
      include "bln.go.connect.inc.f90"    ! specify native contacts 
      include "bln.go.ab.inc.f90"         ! parameters for the non-bonded LJ interaction 
      include "bln.go.cd.inc.f90"         ! parameters for the torsion angle interaction
      ! ==================================================
                CASE("WT")                ! Wild-type (WT) original frustrated BLN model 
      ! ==================================================
      include "bln.wt.ab.inc.f90"         ! L-J interaction between non-bonded particles
      include "bln.go.cd.inc.f90"         ! Torsion angle potential
      ! ==================================================
        ENDSELECT
        ! }}}
        ! cic (Internal Coordinates ) {{{
        R=RESHAPE(X,(/ N,3 /))

! INTER-PARTICLE DISTANCES: DR, LEN_DR; BOND VECTORS: BVR {{{

        DO I = 1, N-1
          DO J = I+1, N
            DR(I,J,1:3) = R(I,1:3) - R(J,1:3)
            DR(J,I,1:3) = -DR(I,J,1:3)
            LEN_DR(I,J) = SQRT(SUM(DR(I,J,1:3)**2))
            LEN_DR(J,I) = LEN_DR(I,J)
          ENDDO
            BVR(I,1:3)=DR(I+1,I,1:3)
        ENDDO
! }}}
! BOND VECTORS: DPD, B, EB  {{{

      DO I = 1, N-1
        B(I)= SQRT(SUM(BVR(I,1:3)**2))
        EB(I,1:3)=BVR(I,1:3)/B(I)
        IF (I.LT.N-1) DPD(I)= SUM(BVR(I,1:3)*BVR(I+1,1:3))
      ENDDO
! }}}
! Cross-products between adjacent bond vectors i and i+1 {{{

        DO I = 1, N-2
           XPD_2(I) = (B(I)*B(I+1))**2-DPD(I)**2 
           XPD(I)=SQRT(XPD_2(I))
           VXPD(I,1)=BVR(I,2)*BVR(I+1,3)-BVR(I,3)*BVR(I+1,2)
           VXPD(I,2)=BVR(I,3)*BVR(I+1,1)-BVR(I,1)*BVR(I+1,3)
           VXPD(I,3)=BVR(I,1)*BVR(I+1,2)-BVR(I,2)*BVR(I+1,1)
           HVXPD(I,1:3)=VXPD(I,1:3)/XPD(I) 
           PP(I,1:3)=VXPD(I,1:3)/XPD_2(I)
        ENDDO

! }}}
! BOND ANGLES: ANG(I,1), I=2,...,N-1 {{{

        DO I = 2, N-1
            COS_THETA=-DPD(I-1)/(B(I-1)*B(I))
            ANG(I,1) = ACOS(COS_THETA)
            F(I,1)=2.0D0*RK_THETA*(ANG(I,1)-THETA_0)
        ENDDO
! }}}
! TORSIONAL ANGLES: ANG(I,2), I=2,...,N-2 {{{

        DO I = 2, N-2
            COS_PHI=-SUM(HVXPD(I-1,1:3)*HVXPD(I,1:3))
            IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=COS_PHI/ABS(COS_PHI)
            ANG(I,2) = ACOS(COS_PHI)
            AN=ANG(I,2)
            F(I,2)=-CD(I,1)*SIN(AN)-3.0*CD(I,2)*SIN(3.0*AN)
            FB(I,1)=F(I,2)*B(I)
            FB(I,2)=F(I,2)/B(I)
        ENDDO
! }}}
! }}}
! ce  ( Calculate Energies ) {{{
            E(1:5)=0.0D0 
! non-bonded: E(2) {{{

        DO I = 1, N-2
          DO J = I+2, N
            RAD(6)=LEN_DR(I,J)**6; 
            RAD(12)=RAD(6)**2
            E(2) = E(2) + 4.0D0*AB(I,J,1)*S(12)/RAD(12)
            E(2) = E(2) + 4.0D0*AB(I,J,2)*S(6)/RAD(6)
          ENDDO
        ENDDO
! }}}
! bonded: E(3) {{{

        DO I = 1, N-1
          E(3) = E(3) + 0.5*RK_R*(B(I)-SIGMA)**2
        ENDDO
! }}}
! bond angles: E(4) {{{

        DO I = 2, N-1
          E(4) = E(4) + 0.5*RK_THETA*(ANG(I,1)-THETA_0)**2
        ENDDO
! }}}
! torsional angles: E(5) {{{

        DO I = 2, N-2
          E(5) = E(5) + CD(I,1)*(1.0 + COS(ANG(I,2)))
          E(5) = E(5) + CD(I,2)*(1.0 + COS(3.0*ANG(I,2)))
        ENDDO
! }}}

        ! Now the total energy
        E(1)=SUM(E(2:5))
        include "deb.ebln_e.i.f90"
        ! }}}
! grad - Gradients {{{
IF (.NOT. GRADT) THEN 
    RETURN
ENDIF

S(6)=S(1)**6 ; S(12)=S(6)**2
GNB= 0.0D0; GB= 0.0D0 ; GBA= 0.0D0 ; GTA= 0.0D0 ; GRAD= 0.0D0 

!include "bln.grad.nbond.inc.f90"  ! GNB         non-bonded 
!include "bln.grad.bond.inc.f90"   ! GB          bonded  
!include "bln.grad.bangle.inc.f90" ! GBA         bond angles
!include "bln.grad.tangle.inc.f90" ! GTA         torsional angles

!  Non-bonded interaction forces ..... {{{

        DO I = 1, N-2
          DO J = I+2, N
  
            RAD(1)=LEN_DR(I,J) ; 
            RAD(2)=RAD(1)**2 ; 
            RAD(7)=RAD(1)**7 ; 
            RAD(14)=RAD(7)**2 ; 
            RAD(8)=RAD(7)*RAD(1)
        
            DF = 2.0*AB(I,J,1)*S(12)/RAD(12) + AB(I,J,2)*S(6)/RAD(6)
            DF=-24.0*DF/RAD(2)

            FRR(1:3) = DF*DR(I,J,1:3) 
            GNB(I,1:3) = FRR(1:3) + GNB(I,1:3)
            GNB(J,1:3) = -FRR(1:3) + GNB(J,1:3)

          ENDDO
        ENDDO
! }}}
! Bond interaction forces ... {{{

        DO I = 1, N-1
           RVAR = SIGMA/B(I) ; DF = RK_R*(1.0 - RVAR) 
           FRR(1:3)=DF*BVR(I,1:3)
           GB(I,1:3) = FRR(1:3)+GB(I,1:3)             
           GB(I+1,1:3) = -FRR(1:3)+GB(I+1,1:3)
        ENDDO
! }}}
! Bond angles {{{
        do i=2,N-1
                AN=ANG(I,1)
                GBA_I(i,1:3)=-(F(i,1)/(B(i-1)*sin(AN)))*(EB(i,1:3)+EB(i-1,1:3)*cos(AN))
                GBA_K(i,1:3)=(F(i,1)/(B(i)*sin(AN)))*(EB(i,1:3)*cos(AN)+EB(i-1,1:3))
                GBA_J(i,1:3)=-GBA_I(i,1:3)-GBA_K(i,1:3)
        enddo

        GBA(1,1:3)=GBA_I(2,1:3)
        GBA(2,1:3)=GBA_I(3,1:3)+GBA_J(2,1:3)

        do i=3,N-2
                GBA(i,1:3)=GBA_I(i+1,1:3)+GBA_J(i,1:3)+GBA_K(i-1,1:3)
        enddo 

        GBA(N-1,1:3)=GBA_J(N-1,1:3)+GBA_K(N-2,1:3)
        GBA(N,1:3)=GBA_K(N-1,1:3)
        ! }}}
        ! Torsional angles {{{
        do i=2,N-2
                GTA_I(i,1:3)=-PP(i-1,1:3)*FB(i,1)
                GTA_L(i,1:3)=PP(i,1:3)*FB(i,1)

                ! use Bekker formulas
                GTA_J(i,1:3)=-GTA_I(i,1:3)*(1.0D0+DPD(i-1)/(B(i)**2)) & 
                        !-GTA_L(i,1:3)*DPD(i+1)/(B(i)**2)
                        -GTA_L(i,1:3)*DPD(i)/(B(i)**2)
                        
                GTA_K(i,1:3)=-GTA_I(i,1:3)-GTA_J(i,1:3)-GTA_L(i,1:3)
        enddo

        GTA(1,1:3)=GTA_I(2,1:3)
        GTA(2,1:3)=GTA_I(3,1:3)+GTA_J(2,1:3)
        GTA(3,1:3)=GTA_I(4,1:3)+GTA_J(3,1:3)+GTA_K(2,1:3)

        do i=4,N-3
                GTA(i,1:3)=GTA_I(i+1,1:3)+GTA_J(i,1:3)+GTA_K(i-1,1:3)+GTA_L(i-2,1:3)
        enddo 

        GTA(N-2,1:3)=GTA_J(N-2,1:3)+GTA_K(N-3,1:3)+GTA_L(N-4,1:3)
        GTA(N-1,1:3)=GTA_K(N-2,1:3)+GTA_L(N-3,1:3)
        GTA(N,1:3)=GTA_L(N-2,1:3)
        ! }}}

G=GNB+GB+GBA+GTA 
GRADX=PACK(G,.TRUE.)


! }}}
        ! debug {{{
        RMS(1)=SQRT(SUM(GRADX**2)/NR)
        RMS(2)=SQRT(SUM(GNB**2)/NR)
        RMS(3)=SQRT(SUM(GB**2)/NR)
        RMS(4)=SQRT(SUM(GBA**2)/NR)
        RMS(5)=SQRT(SUM(GTA**2)/NR)

        include "deb.ebln_g.i.f90"
                ! }}}

!!deam {{{
!DEALLOCATE(R,GRAD)
!deallocate(CONNECT)
!deallocate(AB)
!!if (allocated(CD)) deallocate(CD)
!!deallocate(DR,LEN_DR)
!!deallocate(ANG,F,FB)
!!deallocate(XPD_2,XPD)
!!DEALLOCATE(VXPD,HVXPD,PP,BVR,EB,B,ANG,DPD)
!!DEALLOCATE(G,GB,GBA_I,GBA_J,GBA_K,GTA,GTA_I,GTA_J,GTA_K,GTA_L,GBA,GNB)
!! }}}

        RETURN
        END SUBROUTINE EBLN

        ENDMODULE MODBLN
