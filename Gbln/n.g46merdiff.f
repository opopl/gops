
        subroutine g46merdiff(qo, n, grad, energy, gtest)
        USE MODHESS
        IMPLICIT NONE
        logical gtest, stest
        INTEGER ntype(46), N
        DOUBLE PRECISION QO(3*N), GRAD(3*N), ENERGY
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N),D_PARAM(N),
     1                   c_param(n), rk_theta, rk_r, epsilon, sigma, theta_0, delta, rmass
        parameter (rmass = 40.0, epsilon = 0.0100570)
        parameter (sigma=3.4, delta=1.0d-6, theta_0 = 1.8326)
        parameter (rk_r = 20.0*0.0100570, rk_theta = 20.0*0.0100570)
        DOUBLE PRECISION X(N), Y(N), Z(N), XR(N,N), YR(N,N), ZR(N,N),
     2                  dot_prod(n,3), x_prod(n), bond_angle(n), tor_angle(n), radii(n,n)

        STEST=.FALSE.

     1                            radii,ntype)
     1                            radii,ntype)
        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
     1                            radii,ntype)



        IF (.NOT.STEST) RETURN
     1                            radii,ntype)

        return
        end
        subroutine gparam_array(a_param,b_param,c_param,d_param,n)
        IMPLICIT NONE
        logical connect(46,46)
        INTEGER J, ICOUNT, I, J2, J1, N
        DOUBLE PRECISION NTYPE(46), A_PARAM(N,N), B_PARAM(N,N)
        DOUBLE PRECISION C_PARAM(N), D_PARAM(N), EPSILON
        parameter (epsilon = 0.0100570)

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
     
        DO J1=1,46
           DO J2=J1,46
           ENDDO
        ENDDO


        do i = 1, n-3
        icount = 0

        do j = 0,3
        if(ntype(i+j) .eq. 3)then
        icount = icount + 1
        endif
        enddo

        if(icount .ge. 2)then
        d_param(i+1) = 0.2*epsilon
        else
        d_param(i+1) = 1.2*epsilon
        endif

        icount = 0

        enddo

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
        return
        end
