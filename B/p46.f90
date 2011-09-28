
        SUBROUTINE EP46(FH,DEB,QO, N, GRAD, E, GTEST,GOTYPE)
! dec {{{
        IMPLICIT NONE
        ! sub
        INTEGER,INTENT(IN) :: N
        LOGICAL,INTENT(IN) :: DEB
        INTEGER,INTENT(IN) :: FH
        DOUBLE PRECISION,DIMENSION(3*N),INTENT(IN) :: QO
        DOUBLE PRECISION,DIMENSION(3*N),INTENT(OUT) :: GRAD
        DOUBLE PRECISION,DIMENSION(10) :: E
        LOGICAL,INTENT(IN) :: GTEST
        ! loc
        INTEGER NCALLMAX
        INTEGER,dimension(N) :: NTYPE
        LOGICAL STEST
        LOGICAL GOTYPE
        character(len=10) ptype
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N), &
     &                  x(n), y(n), z(n), xr(n,n), yr(n,n), zr(n,n), &
     &                  dot_prod(n,3), x_prod(n), bond_angle(n), tor_angle(n), radii(n,n)
        !
        DOUBLE PRECISION RMASS, EPSILON,SIGMA,DELTA,THETA_0,RK_R,RK_THETA
        integer,save :: ncall
        DOUBLE PRECISION :: rms
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        ! }}}
        ! body {{{
        ptype="GO"
        STEST=.FALSE.
        CALL PARAM_ARRAY(N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,NTYPE,GOTYPE)
        CALL CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
     &                       BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)
        CALL CALC_ENERGY(FH,DEB,QO,E,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
     &                   BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)
        !IF (NCALL .LE. 0) THEN
            !NCALL=1
          !ELSE
            !NCALL=1+NCALL
        !ENDIF
        !if (ncall .gt. 2000) stop
        !rms=sqrt(sum(grad**2)/(3*N))
        !WRITE(*,*) PTYPE,E(1),NCALL,RMS

        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
        CALL CALC_GRADIENT(FH,DEB,QO,GRAD,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
     &                     BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)
        IF (.NOT.STEST) RETURN
        CALL CALC_DYN(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
     &                BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
        RETURN
        ! }}}
        END SUBROUTINE EP46

        SUBROUTINE PARAM_ARRAY(N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,NTYPE,GOTYPE)
! {{{
! Declarations {{{
        implicit NONE
        
        ! sub  
        INTEGER,INTENT(IN) :: N
        CHARACTER(LEN=10) PTYPE
        INTEGER,DIMENSION(N),INTENT(OUT) :: NTYPE
        DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: A_PARAM,B_PARAM
        DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: C_PARAM,D_PARAM
        LOGICAL GOTYPE
        ! loc
        INTEGER ICOUNT, J, I, j1, j2
        DOUBLE PRECISION EPSILON
        LOGICAL CONNECT(N,N)
        DOUBLE PRECISION :: CON(N,N)
        parameter (epsilon = 0.0100570D0)
! }}}
! Firstly, specify amino acid types by filling in array ntype(:) {{{
! 1 -> Hydrophobic (B)
! 2 -> Hydrophilic (L)
! 3 -> Neutral (N)
        ntype(1:9) = 1
        ntype(10:12) = 3
        ntype(13:19:2) = 2
        ntype(14:20:2) = 1
        ntype(21:23) = 3
        ntype(24:32) = 1
        ntype(33:35) = 3
        ntype(36:46:2) = 2
        ntype(37:45:2) = 1
! }}}
! Specify parameters for the dihedral angle potential. {{{
!       These are stored in arrays
!       c_param(:), d_param(:)
        do i = 1, n-3
        icount = 0
        do j = 0,3
        if(ntype(i+j) .eq. 3)then
        icount = icount + 1
        endif
        enddo
        if(icount .ge. 2)then
        c_param(i+1) = 0.0
        d_param(i+1) = 0.2*epsilon
        else
        c_param(i+1) = 1.2*epsilon
        d_param(i+1) = 1.2*epsilon
        endif
        icount = 0
        enddo
!}}}
! Go-like model connectivities: fill in array CONNECT(:,:) {{{
!
        con=0.0D0
        IF (GOTYPE) THEN 
          ! {{{
        CON(20, 1)=1.0D0
        CON(24, 1)=1.0D0
        CON(45, 1)=1.0D0
        CON(24, 2)=1.0D0
        CON(43, 2)=1.0D0
        CON(45, 2)=1.0D0
        CON(18, 3)=1.0D0
        CON(20, 3)=1.0D0
        CON(25, 3)=1.0D0
        CON(26, 3)=1.0D0
        CON(43, 3)=1.0D0
        CON(26, 4)=1.0D0
        CON(41, 4)=1.0D0
        CON(16, 5)=1.0D0
        CON(18, 5)=1.0D0
        CON(26, 5)=1.0D0
        CON(27, 5)=1.0D0
        CON(28, 5)=1.0D0
        CON(41, 5)=1.0D0
        CON(28, 6)=1.0D0
        CON(39, 6)=1.0D0
        CON(16, 7)=1.0D0
        CON(28, 7)=1.0D0
        CON(29, 7)=1.0D0
        CON(30, 7)=1.0D0
        CON(39, 7)=1.0D0
        CON(30, 8)=1.0D0
        CON(37, 8)=1.0D0
        CON(14, 9)=1.0D0
        CON(30, 9)=1.0D0
        CON(31, 9)=1.0D0
        CON(32, 9)=1.0D0
        CON(37, 9)=1.0D0
        CON(30, 14)=1.0D0
        CON(31, 14)=1.0D0
        CON(28, 16)=1.0D0
        CON(29, 16)=1.0D0
        CON(26, 18)=1.0D0
        CON(24, 20)=1.0D0
        CON(25, 20)=1.0D0
        CON(45, 24)=1.0D0
        CON(41, 26)=1.0D0
        CON(43, 26)=1.0D0
        CON(39, 28)=1.0D0
        CON(41, 28)=1.0D0
        CON(39, 30)=1.0D0
        CON(37, 32)=1.0D0
        CON(1, 20)=1.0D0
        CON(1, 24)=1.0D0
        CON(1, 45)=1.0D0
        CON(2, 24)=1.0D0
        CON(2, 43)=1.0D0
        CON(2, 45)=1.0D0
        CON(3, 18)=1.0D0
        CON(3, 20)=1.0D0
        CON(3, 25)=1.0D0
        CON(3, 26)=1.0D0
        CON(3, 43)=1.0D0
        CON(4, 26)=1.0D0
        CON(4, 41)=1.0D0
        CON(5, 16)=1.0D0
        CON(5, 18)=1.0D0
        CON(5, 26)=1.0D0
        CON(5, 27)=1.0D0
        CON(5, 28)=1.0D0
        CON(5, 41)=1.0D0
        CON(6, 28)=1.0D0
        CON(6, 39)=1.0D0
        CON(7, 16)=1.0D0
        CON(7, 28)=1.0D0
        CON(7, 29)=1.0D0
        CON(7, 30)=1.0D0
        CON(7, 39)=1.0D0
        CON(8, 30)=1.0D0
        CON(8, 37)=1.0D0
        CON(9, 14)=1.0D0
        CON(9, 30)=1.0D0
        CON(9, 31)=1.0D0
        CON(9, 32)=1.0D0
        CON(9, 37)=1.0D0
        CON(14, 30)=1.0D0
        CON(14, 31)=1.0D0
        CON(16, 28)=1.0D0
        CON(16, 29)=1.0D0
        CON(18, 26)=1.0D0
        CON(20, 24)=1.0D0
        CON(20, 25)=1.0D0
        CON(24, 45)=1.0D0
        CON(26, 41)=1.0D0
        CON(26, 43)=1.0D0
        CON(28, 39)=1.0D0
        CON(28, 41)=1.0D0
        CON(30, 39)=1.0D0
        CON(32, 37)=1.0D0
        !}}}
      else
        con=1.0D0
      endif

! }}}
! Specify parameters for the L-J interaction between non-bonded particles. {{{
!       These are stored in arrays
!       a_param(:,:), b_param(:,:)
        do i = 1, n-1
	        do j = i+1, n
			        if (ntype(i) .eq. 3 .or. ntype(j) .eq. 3)then
				        a_param(i,j) = 1.0*epsilon
				        b_param(i,j) = 0.0
				        a_param(j,i) = 1.0*epsilon
				        b_param(j,i) = 0.0
			        elseif (ntype(i) .eq. 1 .and. ntype(j) .eq. 1)then
				        a_param(i,j) =  epsilon
				        a_param(j,i) =  epsilon
				        b_param(i,j) = -epsilon*con(i,j)
				        b_param(j,i) = -epsilon*con(i,j)
   			        else
				        a_param(i,j) = epsilon*2.0/3.0
				        b_param(i,j) = epsilon*2.0/3.0
				        a_param(j,i) = epsilon*2.0/3.0
				        b_param(j,i) = epsilon*2.0/3.0
			        endif
            enddo
        enddo
!  !! }}}

!}}}
!}}}
        return
        end

! }}}
!> Calculate the internal coordinates
        SUBROUTINE CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
     &                             BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,ptype)
! {{{
        IMPLICIT NONE

        ! sub
        INTEGER,INTENT(IN) :: N 
        DOUBLE PRECISION,DIMENSION(3*N),INTENT(IN) :: QO
        DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A_PARAM,B_PARAM
        DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: C_PARAM,D_PARAM
        DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: X,Y,Z
        DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: XR,YR,ZR
        DOUBLE PRECISION,DIMENSION(N,3),INTENT(OUT) :: DOT_PROD
        DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: X_PROD, BOND_ANGLE, TOR_ANGLE
        DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: RADII 
        INTEGER,dimension(N),INTENT(IN) :: Ntype
        CHARACTER(LEN=*),INTENT(IN) :: PTYPE
        ! local 
        DOUBLE PRECISION :: COS_THETA,COS_PHI
        INTEGER I,J

        do i = 1, n
        j = (i-1)*3
        x(i) = qo(j+1)
        y(i) = qo(j+2)
        z(i) = qo(j+3)
        enddo
! Inter-particle distances
        do i = 1, n-1
        do j = i+1, n
        xr(i,j) = x(j) - x(i)
        yr(i,j) = y(j) - y(i)
        zr(i,j) = z(j) - z(i)
        radii(i,j) = dsqrt(xr(i,j)*xr(i,j) + yr(i,j)*yr(i,j) &
     &  + zr(i,j)*zr(i,j))
        radii(j,i) = radii(i,j)
        enddo
        enddo
! Dot products between bond vectors
        do i = 1, n-3
        dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) + &
     &  zr(i,i+1)*zr(i,i+1)
        dot_prod(i,2) = xr(i,i+1)*xr(i+1,i+2)+yr(i,i+1)*yr(i+1,i+2)+ &
     &  zr(i,i+1)*zr(i+1,i+2)
        dot_prod(i,3) = xr(i,i+1)*xr(i+2,i+3)+yr(i,i+1)*yr(i+2,i+3)+ &
     &  zr(i,i+1)*zr(i+2,i+3)
        enddo
        i = n-2
        dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) + &
     &  zr(i,i+1)*zr(i,i+1)
        dot_prod(i,2) = xr(i,i+1)*xr(i+1,i+2)+yr(i,i+1)*yr(i+1,i+2)+ &
     &  zr(i,i+1)*zr(i+1,i+2)
        i = n-1
        dot_prod(i,1) = xr(i,i+1)*xr(i,i+1) + yr(i,i+1)*yr(i,i+1) + &
     &  zr(i,i+1)*zr(i,i+1)
! Cross-products between adjacent bond vectors
        do i = 1, n-2
        x_prod(i) = dot_prod(i,1)*dot_prod(i+1,1) - &
     &               dot_prod(i,2)*dot_prod(i,2)
        enddo
! Bond angles
        do i = 1, n-2
        cos_theta=-dot_prod(i,2)/(dsqrt(dot_prod(i,1) &
     &  *dot_prod(i+1,1)))
        bond_angle(i+1) = dacos(cos_theta)
        enddo
! Torsional angles
        do i = 1, n-3
        cos_phi = (dot_prod(i,2)*dot_prod(i+1,2) - &
     &  dot_prod(i,3)*dot_prod(i+1,1))/dsqrt(x_prod(i)*x_prod(i+1))
        IF (ABS(cos_phi).GT.1.0D0) cos_phi=cos_phi/abs(cos_phi)
        tor_angle(i+1) = dacos(cos_phi)
!       WRITE(*,'(A,I4,4F20.10)') 'i,tor_angle,cos_phi,dacos=',i,tor_angle(i+1),cos_phi,dacos(cos_phi)
        enddo
        return
        ENDSUBROUTINE CALC_INT_COORDS
! }}}
!> Calculate the energy
        SUBROUTINE CALC_ENERGY(FH,DEB,QO,E,N,&
            & A_PARAM,B_PARAM,C_PARAM,D_PARAM,&
            & X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
            & BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)
! dec {{{

        ! sub {{{
            IMPLICIT NONE
            LOGICAL DEB
            INTEGER,INTENT(IN) :: FH
            DOUBLE PRECISION,DIMENSION(10) :: E 
            DOUBLE PRECISION,DIMENSION(3*N),INTENT(IN) :: QO
            INTEGER,INTENT(IN) :: N
            DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A_PARAM,B_PARAM
            DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: C_PARAM,D_PARAM
            DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X,Y,Z
            DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: XR,YR,ZR
            DOUBLE PRECISION,DIMENSION(N,3),INTENT(IN) :: DOT_PROD
            DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X_PROD, BOND_ANGLE, TOR_ANGLE
            DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: RADII 
            INTEGER,DIMENSION(N),INTENT(IN) :: NTYPE
            CHARACTER(LEN=*),INTENT(IN) :: PTYPE
            ! }}}
        ! local {{{
        INTEGER I,J
        DOUBLE PRECISION :: ENERGY
        INTEGER, SAVE :: NCALL
        INTEGER NCALLMAX
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, THETA_0, RK_R, RK_THETA, RAD6, E_TANGLE, &
     &                    S6, E_NBOND, E_BOND, E_BANGLE
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        ! }}}
        ! }}}
        ! body {{{

        include "ncallmax.i.f90"
        
        s6=sigma**6
        e_nbond=0.0D0
        e_bond=0.0D0
        e_bangle=0.0D0
        e_tangle=0.0D0
        do i = 1, n-2
        do j = i+2, n
        rad6 = radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)* &
     &  radii(i,j)
        e_nbond = e_nbond + 4.0*((a_param(i,j)*s6*s6/(rad6*rad6)) + &
     &  (b_param(i,j)*s6/rad6))
        enddo
        enddo
        do i = 1, n-1
        e_bond = e_bond + 0.5*rk_r*(radii(i,i+1)-sigma)* &
     &  (radii(i,i+1)-sigma)
        enddo
        do i = 2, n-1
        e_bangle = e_bangle + 0.5*rk_theta*(bond_angle(i)-theta_0) &
     &  *(bond_angle(i)-theta_0)
        enddo
        do i = 2, n-2
        e_tangle = e_tangle + c_param(i)*(1.0 + cos(tor_angle(i))) &
     &  + d_param(i)*(1.0 + cos(3.0*tor_angle(i)))
        enddo
        energy = e_nbond + e_bond + e_bangle + e_tangle
        E(1:5)=(/ ENERGY,E_NBOND,E_BOND,E_BANGLE,E_TANGLE /)

        include "deb.calc_energy.i.f90"
        !       WRITE(*,'(A,4F20.10)') 'nbond,bond,bangle,tangle=',e_nbond,e_bond,e_bangle,e_tangle
        return
        end
! }}}

! Calculate the gradients
        SUBROUTINE CALC_GRADIENT(FH,DEB,QO,FQ,N,&
                & A_PARAM,B_PARAM,C_PARAM,D_PARAM,&
                & X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
                & BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)
        ! {{{
! Declarations {{{
           IMPLICIT NONE
    
            ! sub {{{
            LOGICAL DEB
            DOUBLE PRECISION,DIMENSION(3*N),INTENT(IN) :: QO
            INTEGER,INTENT(IN) :: N,FH
            DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A_PARAM,B_PARAM
            DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: C_PARAM,D_PARAM
            DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X,Y,Z
            DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: XR,YR,ZR
            DOUBLE PRECISION,DIMENSION(N,3),INTENT(IN) :: DOT_PROD
            DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X_PROD, BOND_ANGLE, TOR_ANGLE
            DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: RADII 
            INTEGER,DIMENSION(N),INTENT(IN) :: NTYPE
            CHARACTER(LEN=*),INTENT(IN) :: PTYPE
            ! FQ - gradient
            DOUBLE PRECISION,DIMENSION(3*N),INTENT(OUT) :: FQ
            ! }}}
        ! local {{{
        INTEGER I, J, NR
        INTEGER NCALLMAX
        DOUBLE PRECISION :: RMS(10)
        INTEGER, SAVE :: NCALL
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, RK_R, RK_THETA, &
     &                   A4, COEF, COEF1, COEF2, COEF3, A3, DEN2, A2, A1, DEN1, RNUM, DEN, &
     &                   RVAR, FZZ, FYY, FXX, DF, RAD14, RAD7, S6, theta_0
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6,theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)

        DOUBLE PRECISION,DIMENSION(N) :: FNB_X,FNB_Y,FNB_Z,FB_X,FB_Y,FB_Z
        DOUBLE PRECISION,DIMENSION(N) :: FBA_X,FBA_Y,FBA_Z, FX,FY,FZ
        DOUBLE PRECISION,DIMENSION(N) :: FTA_X,FTA_Y,FTA_Z
        ! }}}
! }}}
                
        include "ncallmax.i.f90"
                
        s6 = sigma**6
        NR=3*N
! Gradients of potential :=0 {{{
        do i = 1,n
        fnb_x(i) = 0.0
        fnb_y(i) = 0.0
        fnb_z(i) = 0.0
        fb_x(i)  = 0.0
        fb_y(i)  = 0.0
        fb_z(i)  = 0.0
        fba_x(i) = 0.0
        fba_y(i) = 0.0
        fba_z(i) = 0.0
        fta_x(i) = 0.0
        fta_y(i) = 0.0
        fta_z(i) = 0.0
        fx(i)= 0.0
        fy(i)= 0.0
        fz(i)= 0.0
        enddo
! }}}
        ! ..... Non-bonded interaction forces ..... {{{
        do i = 1, n-2
        do j = i+2, n
        rad7 = radii(i,j)*radii(i,j)*radii(i,j)*radii(i,j)* &
     &        radii(i,j)*radii(i,j)*radii(i,j)
        rad14 = rad7*rad7
        df = -24.0*((2.0*a_param(i,j)*s6*s6/rad14) + &
     &              (b_param(i,j)*s6/(rad7*radii(i,j))))
        fxx = df*xr(i,j)
        fyy = df*yr(i,j)
        fzz = df*zr(i,j)
        fnb_x(i) = fxx + fnb_x(i)
        fnb_y(i) = fyy + fnb_y(i)
        fnb_z(i) = fzz + fnb_z(i)
        fnb_x(j) = -fxx + fnb_x(j)
        fnb_y(j) = -fyy + fnb_y(j)
        fnb_z(j) = -fzz + fnb_z(j)
        enddo
        enddo
! }}}
! ... Bond interaction forces ... {{{
        do i = 1, n-1
        rvar = sigma/radii(i,i+1)
        df = rk_r*(1.0 - rvar)
        fxx = df*xr(i,i+1)
        fyy = df*yr(i,i+1)
        fzz = df*zr(i,i+1)
        fb_x(i) = fxx + fb_x(i)
        fb_y(i) = fyy + fb_y(i)
        fb_z(i) = fzz + fb_z(i)
        fb_x(i+1) = -fxx + fb_x(i+1)
        fb_y(i+1) = -fyy + fb_y(i+1)
        fb_z(i+1) = -fzz + fb_z(i+1)
        enddo
! }}}
! bond angle forces {{{ 
! particle 1
! particles 1,2,n-1, and n done outside of the loop
        i=1
        den = dsin(bond_angle(i+1)) &
     &        *dsqrt(dot_prod(i+1,1)*dot_prod(i,1))
        rnum = rk_theta*(bond_angle(i+1) - theta_0)
        fba_x(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*xr(i,i+1) - &
     &        xr(i+1,i+2))/den
        fba_y(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*yr(i,i+1) - &
     &        yr(i+1,i+2))/den
        fba_z(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*zr(i,i+1) - &
     &        zr(i+1,i+2))/den
! particle 2
        i = 2
        den = dsin(bond_angle(i)) &
     &        *dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i+1))*dsqrt(dot_prod(i+1,1) &
     &         *dot_prod(i,1))
        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/ &
     &  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1)) &
     &        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den
        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/ &
     &        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1
        fba_x(i) = a1 + a2
        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/ &
     &  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1)) &
     &        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den
        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/ &
     &        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1
        fba_y(i) = a1 + a2
        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/ &
     &  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1)) &
     &        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den
        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/ &
     &        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1
        fba_z(i) = a1 + a2
! particles 3 thru n-2
        do i = 3, n-2
        den = dsin(bond_angle(i))* &
     &              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i+1))* &
     &               dsqrt(dot_prod(i+1,1)*dot_prod(i,1))
        den2 = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1) &
     &         *dot_prod(i-1,1))
        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/ &
     &  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1)) &
     &        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den
        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/ &
     &        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1
        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/ &
     &        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den2
        fba_x(i) = a1 + a2 + a3
        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/ &
     &  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1)) &
     &        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den
        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/ &
     &        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1
        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/ &
     &        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den2
        fba_y(i) = a1 + a2 + a3
        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/ &
     &  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1)) &
     &        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den
        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/ &
     &        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1
        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/ &
     &        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den2
        fba_z(i) = a1 + a2 + a3
        enddo
! particle n-1
        i = n-1
        den = dsin(bond_angle(i))* &
     &              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1) &
     &         *dot_prod(i-1,1))
        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/ &
     &  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1)) &
     &        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den
        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/ &
     &        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den1
        fba_x(i) = a1 + a2
        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/ &
     &  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1)) &
     &        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den
        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/ &
     &        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den1
        fba_y(i) = a1 + a2
        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/ &
     &  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1)) &
     &        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den
        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/ &
     &        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den1
        fba_z(i) = a1 + a2
! particle n
        i = n
        den = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1) &
     &        *dot_prod(i-1,1))
        fba_x(i) = rk_theta*(bond_angle(i-1) - theta_0)* &
     &        ((dot_prod(i-2,2)/dot_prod(i-1,1))*xr(i-1,i) &
     &        - xr(i-2,i-1))/den
        fba_y(i) = rk_theta*(bond_angle(i-1) - theta_0)* &
     &        ((dot_prod(i-2,2)/dot_prod(i-1,1))*yr(i-1,i) &
     &        - yr(i-2,i-1))/den
        fba_z(i) = rk_theta*(bond_angle(i-1) - theta_0)* &
     &        ((dot_prod(i-2,2)/dot_prod(i-1,1))*zr(i-1,i) &
     &        - zr(i-2,i-1))/den
! }}}
! Torsional angle forces ! {{{
! particles 1, 2, 3, n-2, n-1, and n are done outside of the loop! particle 1
        i = 1
             coef =(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1)) &
     &         *dcos(tor_angle(i+1))-3.0)) &
     &  *(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))
        fta_x(i) = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) + &
     &         dot_prod(i+1,1)*xr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) + &
     &        dot_prod(i,2)*xr(i+1,i+2)))
        fta_y(i) = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) + &
     &        dot_prod(i+1,1)*yr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) + &
     &        dot_prod(i,2)*yr(i+1,i+2)))
        fta_z(i) = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) + &
     &        dot_prod(i+1,1)*zr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) + &
     &        dot_prod(i,2)*zr(i+1,i+2)))
! particle 2
        i = 2
        coef =(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1)) &
     &  *dcos(tor_angle(i+1)) - 3.0)) &
     &        *(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))
             coef1 = (c_param(i) + d_param(i)*(12.0*dcos(tor_angle(i)) &
     &        *dcos(tor_angle(i)) - &
     &        3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))
        a1 =  -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) + &
     &        dot_prod(i+1,1)*xr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) + &
     &        dot_prod(i,2)*xr(i+1,i+2)))
        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) + &
     &        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) - &
     &        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) - &
     &        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) + &
     &        dot_prod(i-1,2)*xr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) + &
     &        dot_prod(i,2)*xr(i+1,i+2)))
        fta_x(i) = a1 + a2
        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) + &
     &        dot_prod(i+1,1)*yr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) + &
     &        dot_prod(i,2)*yr(i+1,i+2)))
        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) + &
     &        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) - &
     &        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) - &
     &        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) + &
     &        dot_prod(i-1,2)*yr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) + &
     &        dot_prod(i,2)*yr(i+1,i+2)))
        fta_y(i) = a1 + a2
        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) + &
     &        dot_prod(i+1,1)*zr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) + &
     &        dot_prod(i,2)*zr(i+1,i+2)))
        a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) + &
     &        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) - &
     &        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) - &
     &        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) + &
     &        dot_prod(i-1,2)*zr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) + &
     &        dot_prod(i,2)*zr(i+1,i+2)))
        fta_z(i) = a1 + a2
! particle 3
        i = 3
        coef=(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1)) &
     &  *dcos(tor_angle(i+1)) - &
     &        3.0))*(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))
        coef1=(c_param(i)+d_param(i)*(12.0*dcos(tor_angle(i)) &
     &        *dcos(tor_angle(i)) - &
     &        3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))
        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1)) &
     &        *dcos(tor_angle(i-1)) - &
     &        3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))
        a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) + &
     &        dot_prod(i+1,1)*xr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) + &
     &        dot_prod(i,2)*xr(i+1,i+2)))
        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) + &
     &        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) - &
     &        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) - &
     &        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) + &
     &        dot_prod(i-1,2)*xr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) + &
     &        dot_prod(i,2)*xr(i+1,i+2)))
        a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) - &
     &        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) + &
     &        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) - &
     &        dot_prod(i-2,2)*xr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) - &
     &        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) + &
     &        dot_prod(i-1,2)*xr(i-1,i)))
        fta_x(i) = a1 + a2 + a3
        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) + &
     &        dot_prod(i+1,1)*yr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) + &
     &        dot_prod(i,2)*yr(i+1,i+2)))
        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) + &
     &        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) - &
     &        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) - &
     &        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) + &
     &        dot_prod(i-1,2)*yr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) + &
     &        dot_prod(i,2)*yr(i+1,i+2)))
        a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) - &
     &        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) + &
     &        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) - &
     &        dot_prod(i-2,2)*yr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) - &
     &        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) + &
     &        dot_prod(i-1,2)*yr(i-1,i)))
        fta_y(i) = a1 + a2 + a3
        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) + &
     &        dot_prod(i+1,1)*zr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) + &
     &        dot_prod(i,2)*zr(i+1,i+2)))
        a2 =  -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) + &
     &        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) - &
     &        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) - &
     &        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) + &
     &        dot_prod(i-1,2)*zr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) + &
     &        dot_prod(i,2)*zr(i+1,i+2)))
        a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) - &
     &        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) + &
     &        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) - &
     &        dot_prod(i-2,2)*zr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) - &
     &        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) + &
     &        dot_prod(i-1,2)*zr(i-1,i)))
        fta_z(i) = a1 + a2 + a3
! particles 4 to n-3
        do i = 4, n-3
        coef = (c_param(i+1) + d_param(i+1)*(12.0*dcos(tor_angle(i+1)) &
     &  *dcos(tor_angle(i+1)) - 3.0))*(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))
        coef1 = (c_param(i) + d_param(i)*(12.0*dcos(tor_angle(i)) &
     &        *dcos(tor_angle(i)) -3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))
        coef2 = (c_param(i-1) + d_param(i-1)*(12.0*dcos(tor_angle(i-1)) &
     &        *dcos(tor_angle(i-1)) - &
     &  3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))
        coef3 = (c_param(i-2) + d_param(i-2)*(12.0*dcos(tor_angle(i-2)) &
     &        *dcos(tor_angle(i-2)) - &
     &  3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))
        a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) + &
     &        dot_prod(i+1,1)*xr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) + &
     &        dot_prod(i,2)*xr(i+1,i+2)))
        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) + &
     &        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) - &
     &        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) - &
     &        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) + &
     &        dot_prod(i-1,2)*xr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) + &
     &        dot_prod(i,2)*xr(i+1,i+2)))
        a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) - &
     &        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) + &
     &        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) - &
     &        dot_prod(i-2,2)*xr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &  dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) - &
     &        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) + &
     &        dot_prod(i-1,2)*xr(i-1,i)))
        a4 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) - &
     &        dot_prod(i-2,1)*xr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) - &
     &        dot_prod(i-2,2)*xr(i-2,i-1)))
        fta_x(i) = a1 + a2 + a3 + a4
        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) + &
     &        dot_prod(i+1,1)*yr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) + &
     &        dot_prod(i,2)*yr(i+1,i+2)))
        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) + &
     &        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) - &
     &        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) - &
     &        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) + &
     &        dot_prod(i-1,2)*yr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) + &
     &        dot_prod(i,2)*yr(i+1,i+2)))
        a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) - &
     &        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) + &
     &        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) - &
     &        dot_prod(i-2,2)*yr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) - &
     &        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) + &
     &        dot_prod(i-1,2)*yr(i-1,i)))
        a4 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - &
     &        dot_prod(i-2,1)*yr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) - &
     &        dot_prod(i-2,2)*yr(i-2,i-1)))
        fta_y(i) = a1 + a2 + a3 + a4
        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) + &
     &        dot_prod(i+1,1)*zr(i+2,i+3) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) - &
     &        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) + &
     &        dot_prod(i,2)*zr(i+1,i+2)))
        a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) + &
     &        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) - &
     &        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) - &
     &        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) + &
     &        dot_prod(i-1,2)*zr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) + &
     &        dot_prod(i,2)*zr(i+1,i+2)))
        a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) - &
     &        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) + &
     &        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) - &
     &        dot_prod(i-2,2)*zr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) - &
     &        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) + &
     &        dot_prod(i-1,2)*zr(i-1,i)))
        a4 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - &
     &        dot_prod(i-2,1)*zr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) - &
     &        dot_prod(i-2,2)*zr(i-2,i-1)))
        fta_z(i) = a1 + a2 + a3 + a4
        enddo
! particle n-2
        i = n-2
        coef1=(c_param(i)+d_param(i)*(12.0*dcos(tor_angle(i)) &
     &        *dcos(tor_angle(i)) - &
     &        3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))
        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1)) &
     &        *dcos(tor_angle(i-1)) - &
     &        3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))
        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2)) &
     &        *dcos(tor_angle(i-2)) - &
     &  3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))
        a1 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) + &
     &        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) - &
     &        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) - &
     &        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) + &
     &        dot_prod(i-1,2)*xr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) + &
     &        dot_prod(i,2)*xr(i+1,i+2)))
        a2 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) - &
     &        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) + &
     &        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) - &
     &        dot_prod(i-2,2)*xr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) - &
     &        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) + &
     &        dot_prod(i-1,2)*xr(i-1,i)))
        a3 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) - &
     &        dot_prod(i-2,1)*xr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) - &
     &        dot_prod(i-2,2)*xr(i-2,i-1)))
        fta_x(i) = a1 + a2 + a3
        a1 =  -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) + &
     &        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) - &
     &        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) - &
     &        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) + &
     &        dot_prod(i-1,2)*yr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) + &
     &        dot_prod(i,2)*yr(i+1,i+2)))
        a2 =  -coef2*(dot_prod(i-2,2)*yr(i,i+1) - &
     &        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) + &
     &        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) - &
     &        dot_prod(i-2,2)*yr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) - &
     &        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) + &
     &        dot_prod(i-1,2)*yr(i-1,i)))
        a3 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - &
     &        dot_prod(i-2,1)*yr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) - &
     &        dot_prod(i-2,2)*yr(i-2,i-1)))
        fta_y(i) = a1 + a2 + a3
        a1 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) + &
     &        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) - &
     &        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) - &
     &        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) + &
     &        dot_prod(i-1,2)*zr(i-1,i)) - &
     &        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) - &
     &        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) + &
     &        dot_prod(i,2)*zr(i+1,i+2)))
        a2 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) - &
     &        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) + &
     &        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) - &
     &        dot_prod(i-2,2)*zr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) - &
     &        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) + &
     &        dot_prod(i-1,2)*zr(i-1,i)))
        a3 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - &
     &        dot_prod(i-2,1)*zr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) - &
     &        dot_prod(i-2,2)*zr(i-2,i-1)))
        fta_z(i) = a1 + a2 + a3
! particle n-1
        i = n-1
        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1)) &
     &        *dcos(tor_angle(i-1)) - &
     &        3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))
        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2)) &
     &        *dcos(tor_angle(i-2)) - &
     &        3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))
        a1 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) - &
     &        dot_prod(i-2,2)*xr(i-1,i) + &
     &        dot_prod(i-1,2)*xr(i-2,i-1) +  dot_prod(i-1,1)*xr(i-2,i-1) - &
     &        2.0*dot_prod(i-2,3)*xr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) - &
     &        dot_prod(i-2,2)*xr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) - &
     &        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) + &
     &        dot_prod(i-1,2)*xr(i-1,i)))
        a2 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) - &
     &        dot_prod(i-2,1)*xr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) - &
     &        dot_prod(i-2,2)*xr(i-2,i-1)))
        fta_x(i) = a1 + a2
        a1 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) - &
     &        dot_prod(i-2,2)*yr(i-1,i) + &
     &        dot_prod(i-1,2)*yr(i-2,i-1) +  dot_prod(i-1,1)*yr(i-2,i-1) - &
     &        2.0*dot_prod(i-2,3)*yr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) - &
     &        dot_prod(i-2,2)*yr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) - &
     &  dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) + &
     &        dot_prod(i-1,2)*yr(i-1,i)))
        a2 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - &
     &        dot_prod(i-2,1)*yr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) - &
     &        dot_prod(i-2,2)*yr(i-2,i-1)))
        fta_y(i) = a1 + a2
        a1 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) - &
     &        dot_prod(i-2,2)*zr(i-1,i) + &
     &        dot_prod(i-1,2)*zr(i-2,i-1) +  dot_prod(i-1,1)*zr(i-2,i-1) - &
     &        2.0*dot_prod(i-2,3)*zr(i-1,i) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) - &
     &        dot_prod(i-2,2)*zr(i-2,i-1)) - &
     &        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) - &
     &        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) - &
     &        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) + &
     &        dot_prod(i-1,2)*zr(i-1,i)))
        a2 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - &
     &        dot_prod(i-2,1)*zr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) - &
     &        dot_prod(i-2,2)*zr(i-2,i-1)))
        fta_z(i) = a1 + a2
! particle n
        i = n
        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2)) &
     &        *dcos(tor_angle(i-2)) - &
     &        3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))
        fta_x(i) = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) &
     &        - dot_prod(i-2,1)*xr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) - &
     &        dot_prod(i-2,2)*xr(i-2,i-1)))
        fta_y(i) = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - &
     &        dot_prod(i-2,1)*yr(i-3,i-2) &
     &        - (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) &
     &        - dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) - &
     &        dot_prod(i-2,2)*yr(i-2,i-1)))
        fta_z(i) = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - &
     &        dot_prod(i-2,1)*zr(i-3,i-2) - &
     &        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) - &
     &        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) - &
     &        dot_prod(i-2,2)*zr(i-2,i-1)))
        ! }}}
! Total up the gradients {{{
        do i = 1, n
        fx(i) = fnb_x(i) + fb_x(i) + fba_x(i) + fta_x(i)
        fy(i) = fnb_y(i) + fb_y(i) + fba_y(i) + fta_y(i)
        fz(i) = fnb_z(i) + fb_z(i) + fba_z(i) + fta_z(i)
        enddo
        do i = 1, n
        j = (i-1)*3
        fq(j+1) = -fx(i)
        fq(j+2) = -fy(i)
        fq(j+3) = -fz(i)
        enddo
        RMS(1)=SQRT(SUM(FX**2+FY**2+FZ**2)/NR)
        RMS(2)=SQRT(SUM(FNB_X**2+FNB_Y**2+FNB_Z**2)/NR)
        RMS(3)=SQRT(SUM(FB_X**2+FB_Y**2+FB_Z**2)/NR)
        RMS(4)=SQRT(SUM(FBA_X**2+FBA_Y**2+FBA_Z**2)/NR)
        RMS(5)=SQRT(SUM(FTA_X**2+FTA_Y**2+FTA_Z**2)/NR)

        include "deb.calc_gradient.i.f90"

        return
        end
! }}}
        ! }}}

!> Calculate the second derivative matrix (two-sided numerical approach)
        SUBROUTINE CALC_DYN(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
     &                      BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)
! {{{! Declarations {{{
        USE V,ONLY: HESS

        IMPLICIT NONE

        ! sub {{{
            !LOGICAL DEB
            !INTEGER,INTENT(IN) :: FH
            DOUBLE PRECISION,DIMENSION(3*N),INTENT(INOUT) :: QO
            INTEGER,INTENT(IN) :: N
            DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A_PARAM,B_PARAM
            DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: C_PARAM,D_PARAM
            DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X,Y,Z
            DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: XR,YR,ZR
            DOUBLE PRECISION,DIMENSION(N,3),INTENT(IN) :: DOT_PROD
            DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X_PROD, BOND_ANGLE, TOR_ANGLE
            DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: RADII 
            INTEGER,DIMENSION(N),INTENT(IN) :: NTYPE
            CHARACTER(LEN=*),INTENT(IN) :: PTYPE
        ! }}}
        ! loc  {{{
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, THETA_0, RK_R, RK_THETA
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4 ,delta=1.0d-4, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        INTEGER I, J
        DOUBLE PRECISION FQ1(3*N), FQ2(3*N)
        ! }}}

! }}} 
!Fill in the Hessian matrix {{{
        do j = 1, 3*n
        qo(j) = qo(j) + delta
        call calc_int_coords(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, &
     &                       bond_angle,tor_angle,radii,ntype)
        call calc_gradient(qo,fq2,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, &
     &                     bond_angle,tor_angle,radii,ntype)
        qo(j) = qo(j) - 2.0*delta
        call calc_int_coords(qo,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, &
     &                       bond_angle,tor_angle,radii,ntype)
        call calc_gradient(qo,fq1,n,a_param,b_param,c_param,d_param,x,y,z,xr,yr,zr,dot_prod,x_prod, &
     &                     bond_angle,tor_angle,radii,ntype)
        qo(j) = qo(j) + delta
        do i = j, 3*n
        HESS(i,j) = (fq2(i) -  fq1(i))/(2.0*delta)
        HESS(j,i) = HESS(i,j)
        enddo
        enddo
! }}}
        return
        end
! }}}
!> Fill the parameter arrays

        ! }}}


