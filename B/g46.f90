

        SUBROUTINE EG46(FH,DEB,QO, N, GRAD, E, GTEST)
! {{{
! declarations {{{
        USE V, ONLY : HESS
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
        logical stest
        INTEGER ntype(46)
        CHARACTER(LEN=10) PTYPE
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N),D_PARAM(N)
        DOUBLE PRECISION ::  c_param(n), rk_theta, rk_r, epsilon, sigma, theta_0, delta, rmass
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        DOUBLE PRECISION X(N), Y(N), Z(N), XR(N,N), YR(N,N), ZR(N,N), &
     &                  dot_prod(n,3), x_prod(n), bond_angle(n), tor_angle(n), radii(n,n)
! }}}
        ! body {{{
        STEST=.FALSE.
        PTYPE="GO"
        CALL GPARAM_ARRAY(A_PARAM,B_PARAM,C_PARAM,D_PARAM,N,NTYPE)
        CALL CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, BOND_ANGLE,TOR_ANGLE, &
     &                            RADII,NTYPE,PTYPE)
        CALL CALC_ENERGY(FH,DEB,QO,E,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, BOND_ANGLE,TOR_ANGLE, &
     &                            RADII,NTYPE,PTYPE)
        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
        CALL CALC_GRADIENT(FH,DEB,QO,GRAD,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, BOND_ANGLE,TOR_ANGLE, &
     &                            RADII,NTYPE,PTYPE)
        IF (.NOT.STEST) RETURN
        CALL CALC_DYN(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, BOND_ANGLE,TOR_ANGLE, &
     &                            RADII,NTYPE)
        return
        ! }}}
        END SUBROUTINE EG46
! }}}

!        SUBROUTINE GPARAM_ARRAY(A_PARAM,B_PARAM,C_PARAM,D_PARAM,N,NTYPE)
!! {{{
!! Declarations {{{
            !! sub
            !IMPLICIT NONE
			!INTEGER,INTENT(IN) :: N
			!INTEGER,DIMENSION(N),INTENT(OUT) :: NTYPE
			!DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: A_PARAM,B_PARAM
			!DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: C_PARAM,D_PARAM
            !! loc
        !LOGICAL CONNECT(N,N)
        !INTEGER J, ICOUNT, I, J2, J1
        !DOUBLE PRECISION EPSILON
        !parameter (epsilon = 0.0100570)
!! }}}
!! Specify amino acid types by filling in the array ntype(:) {{{
        !ntype(1:9) = 1
        !ntype(10:12) = 3
        !ntype(13:19:2) = 2
        !ntype(14:20:2) = 1
        !ntype(21:23) = 3
        !ntype(24:32) = 1
        !ntype(33:35) = 3
        !ntype(36:46:2) = 2
        !ntype(37:45:2) = 1
        !! }}}
!! Go-like model connectivities: fill in array CONNECT(:,:) {{{
!!
        !DO J1=1,46
           !DO J2=J1,46
              !CONNECT(J2,J1)=.FALSE.
           !ENDDO
        !ENDDO
        !CONNECT(20, 1)=.TRUE.
        !CONNECT(24, 1)=.TRUE.
        !CONNECT(45, 1)=.TRUE.
        !CONNECT(24, 2)=.TRUE.
        !CONNECT(43, 2)=.TRUE.
        !CONNECT(45, 2)=.TRUE.
        !CONNECT(18, 3)=.TRUE.
        !CONNECT(20, 3)=.TRUE.
        !CONNECT(25, 3)=.TRUE.
        !CONNECT(26, 3)=.TRUE.
        !CONNECT(43, 3)=.TRUE.
        !CONNECT(26, 4)=.TRUE.
        !CONNECT(41, 4)=.TRUE.
        !CONNECT(16, 5)=.TRUE.
        !CONNECT(18, 5)=.TRUE.
        !CONNECT(26, 5)=.TRUE.
        !CONNECT(27, 5)=.TRUE.
        !CONNECT(28, 5)=.TRUE.
        !CONNECT(41, 5)=.TRUE.
        !CONNECT(28, 6)=.TRUE.
        !CONNECT(39, 6)=.TRUE.
        !CONNECT(16, 7)=.TRUE.
        !CONNECT(28, 7)=.TRUE.
        !CONNECT(29, 7)=.TRUE.
        !CONNECT(30, 7)=.TRUE.
        !CONNECT(39, 7)=.TRUE.
        !CONNECT(30, 8)=.TRUE.
        !CONNECT(37, 8)=.TRUE.
        !CONNECT(14, 9)=.TRUE.
        !CONNECT(30, 9)=.TRUE.
        !CONNECT(31, 9)=.TRUE.
        !CONNECT(32, 9)=.TRUE.
        !CONNECT(37, 9)=.TRUE.
        !CONNECT(30, 14)=.TRUE.
        !CONNECT(31, 14)=.TRUE.
        !CONNECT(28, 16)=.TRUE.
        !CONNECT(29, 16)=.TRUE.
        !CONNECT(26, 18)=.TRUE.
        !CONNECT(24, 20)=.TRUE.
        !CONNECT(25, 20)=.TRUE.
        !CONNECT(45, 24)=.TRUE.
        !CONNECT(41, 26)=.TRUE.
        !CONNECT(43, 26)=.TRUE.
        !CONNECT(39, 28)=.TRUE.
        !CONNECT(41, 28)=.TRUE.
        !CONNECT(39, 30)=.TRUE.
        !CONNECT(37, 32)=.TRUE.
        !CONNECT(1, 20)=.TRUE.
        !CONNECT(1, 24)=.TRUE.
        !CONNECT(1, 45)=.TRUE.
        !CONNECT(2, 24)=.TRUE.
        !CONNECT(2, 43)=.TRUE.
        !CONNECT(2, 45)=.TRUE.
        !CONNECT(3, 18)=.TRUE.
        !CONNECT(3, 20)=.TRUE.
        !CONNECT(3, 25)=.TRUE.
        !CONNECT(3, 26)=.TRUE.
        !CONNECT(3, 43)=.TRUE.
        !CONNECT(4, 26)=.TRUE.
        !CONNECT(4, 41)=.TRUE.
        !CONNECT(5, 16)=.TRUE.
        !CONNECT(5, 18)=.TRUE.
        !CONNECT(5, 26)=.TRUE.
        !CONNECT(5, 27)=.TRUE.
        !CONNECT(5, 28)=.TRUE.
        !CONNECT(5, 41)=.TRUE.
        !CONNECT(6, 28)=.TRUE.
        !CONNECT(6, 39)=.TRUE.
        !CONNECT(7, 16)=.TRUE.
        !CONNECT(7, 28)=.TRUE.
        !CONNECT(7, 29)=.TRUE.
        !CONNECT(7, 30)=.TRUE.
        !CONNECT(7, 39)=.TRUE.
        !CONNECT(8, 30)=.TRUE.
        !CONNECT(8, 37)=.TRUE.
        !CONNECT(9, 14)=.TRUE.
        !CONNECT(9, 30)=.TRUE.
        !CONNECT(9, 31)=.TRUE.
        !CONNECT(9, 32)=.TRUE.
        !CONNECT(9, 37)=.TRUE.
        !CONNECT(14, 30)=.TRUE.
        !CONNECT(14, 31)=.TRUE.
        !CONNECT(16, 28)=.TRUE.
        !CONNECT(16, 29)=.TRUE.
        !CONNECT(18, 26)=.TRUE.
        !CONNECT(20, 24)=.TRUE.
        !CONNECT(20, 25)=.TRUE.
        !CONNECT(24, 45)=.TRUE.
        !CONNECT(26, 41)=.TRUE.
        !CONNECT(26, 43)=.TRUE.
        !CONNECT(28, 39)=.TRUE.
        !CONNECT(28, 41)=.TRUE.
        !CONNECT(30, 39)=.TRUE.
        !CONNECT(32, 37)=.TRUE.
!! }}}
!! Parameters for the dihedral angle potential: fill in arrays c_param(:,:), d_param(:,:) {{{
        !do i = 1, n-3
        !icount = 0
        !do j = 0,3
        !if(ntype(i+j) .eq. 3)then
        !icount = icount + 1
        !endif
        !enddo
        !if(icount .ge. 2)then
        !c_param(i+1) = 0.0
        !d_param(i+1) = 0.2*epsilon
        !else
        !c_param(i+1) = 1.2*epsilon
        !d_param(i+1) = 1.2*epsilon
        !endif
        !icount = 0
        !enddo
!! }}}
!! Parameters for the L-J interaction between non-bonded particles:! arrays a_param(:,:), b_param(:,:)
!! {{{
        !do i = 1, n-1
           !do j = i+1, n
           !if (ntype(i) .eq. 3 .or. ntype(j) .eq. 3) then
             !a_param(i,j) = 1.0*epsilon
             !b_param(i,j) = 0.0
             !a_param(j,i) = 1.0*epsilon
             !b_param(j,i) = 0.0
           !elseif (ntype(i) .eq. 1 .and. ntype(j) .eq. 1)then
             !a_param(i,j) =  epsilon
             !a_param(j,i) =  epsilon
             !IF (CONNECT(I,J)) THEN
                !b_param(i,j) = -epsilon
                !b_param(j,i) = -epsilon
             !ELSE
                !b_param(i,j) = 0.0D0
                !b_param(j,i) = 0.0D0
             !ENDIF
           !else
             !a_param(i,j) = epsilon*2.0/3.0
             !b_param(i,j) = epsilon*2.0/3.0
             !a_param(j,i) = epsilon*2.0/3.0
             !b_param(j,i) = epsilon*2.0/3.0
           !endif
           !enddo
        !enddo
!! }}}
!! }}}
        !return
        !end
!! }}}


        subroutine gparam_array(a_param,b_param,c_param,d_param,n,ntype)
! Declarations {{{
        IMPLICIT NONE
        logical connect(46,46)
        INTEGER J, ICOUNT, I, J2, J1, N
        DOUBLE PRECISION NTYPE(46), A_PARAM(N,N), B_PARAM(N,N)
        DOUBLE PRECISION C_PARAM(N), D_PARAM(N), EPSILON
        parameter (epsilon = 0.0100570)
! }}}
! Specify amino acid types by filling in the array ntype(:) {{{
        ntype(1) = 1
        ntype(2) = 1
        ntype(3) = 1
        ntype(4) = 1
        ntype(5) = 1
        ntype(6) = 1
        ntype(7) = 1
        ntype(8) = 1
        ntype(9) = 1
        ntype(10) = 3
        ntype(11) = 3
        ntype(12) = 3
        ntype(13) = 2
        ntype(14) = 1
        ntype(15) = 2
        ntype(16) = 1
        ntype(17) = 2
        ntype(18) = 1
        ntype(19) = 2
        ntype(20) = 1
        ntype(21) = 3
        ntype(22) = 3
        ntype(23) = 3
        ntype(24) = 1
        ntype(25) = 1
        ntype(26) = 1
        ntype(27) = 1
        ntype(28) = 1
        ntype(29) = 1
        ntype(30) = 1
        ntype(31) = 1
        ntype(32) = 1
        ntype(33) = 3
        ntype(34) = 3
        ntype(35) = 3
        ntype(36) = 2
        ntype(37) = 1
        ntype(38) = 2
        ntype(39) = 1
        ntype(40) = 2
        ntype(41) = 1
        ntype(42) = 2
        ntype(43) = 1
        ntype(44) = 2
        ntype(45) = 1
        ntype(46) = 2
! }}}
! Go-like model connectivities: fill in array CONNECT(:,:) {{{
        DO J1=1,46
           DO J2=J1,46
              CONNECT(J2,J1)=.FALSE.
           ENDDO
        ENDDO
        CONNECT(20, 1)=.TRUE.
        CONNECT(24, 1)=.TRUE.
        CONNECT(45, 1)=.TRUE.
        CONNECT(24, 2)=.TRUE.
        CONNECT(43, 2)=.TRUE.
        CONNECT(45, 2)=.TRUE.
        CONNECT(18, 3)=.TRUE.
        CONNECT(20, 3)=.TRUE.
        CONNECT(25, 3)=.TRUE.
        CONNECT(26, 3)=.TRUE.
        CONNECT(43, 3)=.TRUE.
        CONNECT(26, 4)=.TRUE.
        CONNECT(41, 4)=.TRUE.
        CONNECT(16, 5)=.TRUE.
        CONNECT(18, 5)=.TRUE.
        CONNECT(26, 5)=.TRUE.
        CONNECT(27, 5)=.TRUE.
        CONNECT(28, 5)=.TRUE.
        CONNECT(41, 5)=.TRUE.
        CONNECT(28, 6)=.TRUE.
        CONNECT(39, 6)=.TRUE.
        CONNECT(16, 7)=.TRUE.
        CONNECT(28, 7)=.TRUE.
        CONNECT(29, 7)=.TRUE.
        CONNECT(30, 7)=.TRUE.
        CONNECT(39, 7)=.TRUE.
        CONNECT(30, 8)=.TRUE.
        CONNECT(37, 8)=.TRUE.
        CONNECT(14, 9)=.TRUE.
        CONNECT(30, 9)=.TRUE.
        CONNECT(31, 9)=.TRUE.
        CONNECT(32, 9)=.TRUE.
        CONNECT(37, 9)=.TRUE.
        CONNECT(30, 14)=.TRUE.
        CONNECT(31, 14)=.TRUE.
        CONNECT(28, 16)=.TRUE.
        CONNECT(29, 16)=.TRUE.
        CONNECT(26, 18)=.TRUE.
        CONNECT(24, 20)=.TRUE.
        CONNECT(25, 20)=.TRUE.
        CONNECT(45, 24)=.TRUE.
        CONNECT(41, 26)=.TRUE.
        CONNECT(43, 26)=.TRUE.
        CONNECT(39, 28)=.TRUE.
        CONNECT(41, 28)=.TRUE.
        CONNECT(39, 30)=.TRUE.
        CONNECT(37, 32)=.TRUE.
        CONNECT(1, 20)=.TRUE.
        CONNECT(1, 24)=.TRUE.
        CONNECT(1, 45)=.TRUE.
        CONNECT(2, 24)=.TRUE.
        CONNECT(2, 43)=.TRUE.
        CONNECT(2, 45)=.TRUE.
        CONNECT(3, 18)=.TRUE.
        CONNECT(3, 20)=.TRUE.
        CONNECT(3, 25)=.TRUE.
        CONNECT(3, 26)=.TRUE.
        CONNECT(3, 43)=.TRUE.
        CONNECT(4, 26)=.TRUE.
        CONNECT(4, 41)=.TRUE.
        CONNECT(5, 16)=.TRUE.
        CONNECT(5, 18)=.TRUE.
        CONNECT(5, 26)=.TRUE.
        CONNECT(5, 27)=.TRUE.
        CONNECT(5, 28)=.TRUE.
        CONNECT(5, 41)=.TRUE.
        CONNECT(6, 28)=.TRUE.
        CONNECT(6, 39)=.TRUE.
        CONNECT(7, 16)=.TRUE.
        CONNECT(7, 28)=.TRUE.
        CONNECT(7, 29)=.TRUE.
        CONNECT(7, 30)=.TRUE.
        CONNECT(7, 39)=.TRUE.
        CONNECT(8, 30)=.TRUE.
        CONNECT(8, 37)=.TRUE.
        CONNECT(9, 14)=.TRUE.
        CONNECT(9, 30)=.TRUE.
        CONNECT(9, 31)=.TRUE.
        CONNECT(9, 32)=.TRUE.
        CONNECT(9, 37)=.TRUE.
        CONNECT(14, 30)=.TRUE.
        CONNECT(14, 31)=.TRUE.
        CONNECT(16, 28)=.TRUE.
        CONNECT(16, 29)=.TRUE.
        CONNECT(18, 26)=.TRUE.
        CONNECT(20, 24)=.TRUE.
        CONNECT(20, 25)=.TRUE.
        CONNECT(24, 45)=.TRUE.
        CONNECT(26, 41)=.TRUE.
        CONNECT(26, 43)=.TRUE.
        CONNECT(28, 39)=.TRUE.
        CONNECT(28, 41)=.TRUE.
        CONNECT(30, 39)=.TRUE.
        CONNECT(32, 37)=.TRUE.
! }}}
! Parameters for the dihedral angle potential: fill in arrays c_param(:,:), d_param(:,:) {{{
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
! }}}
! Parameters for the L-J interaction between non-bonded particles:! arrays a_param(:,:), b_param(:,:)! {{{
        do i = 1, n-1
           do j = i+1, n
           if (ntype(i) .eq. 3 .or. ntype(j) .eq. 3) then
             a_param(i,j) = 1.0*epsilon
             b_param(i,j) = 0.0
             a_param(j,i) = 1.0*epsilon
             b_param(j,i) = 0.0
           elseif (ntype(i) .eq. 1 .and. ntype(j) .eq. 1)then
             a_param(i,j) =  epsilon
             a_param(j,i) =  epsilon
             IF (CONNECT(I,J)) THEN
                b_param(i,j) = -epsilon
                b_param(j,i) = -epsilon
             ELSE
                b_param(i,j) = 0.0D0
                b_param(j,i) = 0.0D0
             ENDIF
           else
             a_param(i,j) = epsilon*2.0/3.0
             b_param(i,j) = epsilon*2.0/3.0
             a_param(j,i) = epsilon*2.0/3.0
             b_param(j,i) = epsilon*2.0/3.0
           endif
           enddo
        enddo
! }}}
        return
        end

