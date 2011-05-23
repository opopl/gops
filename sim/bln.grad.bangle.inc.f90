
! sim part {{{
!! bond angle forces: particle 1
!! particles 1, 2, N-1, and N are done outside of the loop
!! {{{
        !i = 1
        !den = sin(ang(i+1,1))*sqrt(dpd(i+1,1)*dpd(i,1))
        !rnum = rk_theta*(ang(i+1,1) - theta_0)
        !fba(i,1:3)=-rnum*((dpd(i,2)/dpd(i,1))*dr(i,i+1,1:3)-dr(i+1,i+2,1:3))/den
!! }}}
     !! particle 2
!! {{{

        !i = 2
        !den = sin(bond_angle(i))
     !1        *dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        !den1 = sin(bond_angle(i+1))*dsqrt(dot_prod(i+1,1)
     !1         *dot_prod(i,1))

        !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        !a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     !1        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

        !fba_x(i) = a1 + a2 

        !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        !a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     !1        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

        !fba_y(i) = a1 + a2 

        !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        !a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     !1        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

        !fba_z(i) = a1 + a2 

!! }}}
!! particles 3 thru n-2 !{{{

        !do i = 3, n-2

        !den = sin(bond_angle(i))*
     !1              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        !den1 = sin(bond_angle(i+1))*
     !1               dsqrt(dot_prod(i+1,1)*dot_prod(i,1))
        !den2 = sin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     !1         *dot_prod(i-1,1))

        !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        !a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     !1        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

        !a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     !1        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den2

        !fba_x(i) = a1 + a2 + a3 

        !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        !a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     !1        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

        !a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     !1        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den2

        !fba_y(i) = a1 + a2 + a3 

        !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        !a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     !1        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

        !a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     !1        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den2

        !fba_z(i) = a1 + a2 + a3 

        !enddo
!!}}}
!! particle n-1 {{{
        !i = n-1
        !den = sin(bond_angle(i))*
     !1              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        !den1 = sin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     !1         *dot_prod(i-1,1))

        !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        !a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     !1        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den1

        !fba_x(i) = a1 + a2

        !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        !a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     !1        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den1

        !fba_y(i) = a1 + a2

        !a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     !1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     !1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        !a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     !1        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den1

        !fba_z(i) = a1 + a2
!!}}}
!! particle n!{{{

        !i = n
        !den = sin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     !1        *dot_prod(i-1,1))

        !fba_x(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     !1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*xr(i-1,i) 
     !1        - xr(i-2,i-1))/den

        !fba_y(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     !1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*yr(i-1,i) 
     !1        - yr(i-2,i-1))/den

        !fba_z(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     !1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*zr(i-1,i) 
     !1        - zr(i-2,i-1))/den!}}}
!! }}}


        do i=1,N
                dot_prod=
        enddo
C bond angle forces  particle 1
C particles 1,2,n-1, and n done outside of the loop
C {{{
        i = 1
        den = dsin(bond_angle(i+1))
     1        *dsqrt(dot_prod(i+1,1)*dot_prod(i,1))
        rnum = rk_theta*(bond_angle(i+1) - theta_0)

        fba_x(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*xr(i,i+1) -
     1        xr(i+1,i+2))/den

        fba_y(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*yr(i,i+1) -
     1        yr(i+1,i+2))/den

        fba_z(i) = -rnum*((dot_prod(i,2)/dot_prod(i,1))*zr(i,i+1) -
     1        zr(i+1,i+2))/den

C particle 2

        i = 2
        den = dsin(bond_angle(i))
     1        *dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i+1))*dsqrt(dot_prod(i+1,1)
     1         *dot_prod(i,1))

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

        fba_x(i) = a1 + a2 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

        fba_y(i) = a1 + a2 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

        fba_z(i) = a1 + a2 

C particles 3 thru n-2 

        do i = 3, n-2

        den = dsin(bond_angle(i))*
     1              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i+1))*
     1               dsqrt(dot_prod(i+1,1)*dot_prod(i,1))
        den2 = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     1         *dot_prod(i-1,1))

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*xr(i,i+1) - xr(i+1,i+2))/den1

        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den2

        fba_x(i) = a1 + a2 + a3 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*yr(i,i+1) - yr(i+1,i+2))/den1

        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den2

        fba_y(i) = a1 + a2 + a3 

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        a2 = -rk_theta*(bond_angle(i+1) - theta_0)*((dot_prod(i,2)/
     1        dot_prod(i,1))*zr(i,i+1) - zr(i+1,i+2))/den1

        a3 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den2

        fba_z(i) = a1 + a2 + a3 

        enddo

C particle n-1 

        i = n-1
        den = dsin(bond_angle(i))*
     1              dsqrt(dot_prod(i,1)*dot_prod(i-1,1))
        den1 = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     1         *dot_prod(i-1,1))

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*xr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *xr(i-1,i) + xr(i,i+1) - xr(i-1,i))/den

        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*xr(i-1,i) - xr(i-2,i-1))/den1

        fba_x(i) = a1 + a2

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*yr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *yr(i-1,i) + yr(i,i+1) - yr(i-1,i))/den

        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*yr(i-1,i) - yr(i-2,i-1))/den1

        fba_y(i) = a1 + a2

        a1 = -rk_theta*(bond_angle(i) - theta_0)*( (dot_prod(i-1,2)/
     1  dot_prod(i,1))*zr(i,i+1) - (dot_prod(i-1,2)/dot_prod(i-1,1))
     1        *zr(i-1,i) + zr(i,i+1) - zr(i-1,i))/den

        a2 = rk_theta*(bond_angle(i-1) - theta_0)*((dot_prod(i-2,2)/
     1        dot_prod(i-1,1))*zr(i-1,i) - zr(i-2,i-1))/den1

        fba_z(i) = a1 + a2

C particle n

        i = n
        den = dsin(bond_angle(i-1))*dsqrt(dot_prod(i-2,1)
     1        *dot_prod(i-1,1))

        fba_x(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*xr(i-1,i) 
     1        - xr(i-2,i-1))/den

        fba_y(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*yr(i-1,i) 
     1        - yr(i-2,i-1))/den

        fba_z(i) = rk_theta*(bond_angle(i-1) - theta_0)*
     1        ((dot_prod(i-2,2)/dot_prod(i-1,1))*zr(i-1,i) 
     1        - zr(i-2,i-1))/den
C }}}


! torsional angle forces
! particles 1, 2, 3, n-2, n-1, and n are done outside of the loop

        DO I=1,N
           IXPD(I)=1.0D0/SQRT(XPD_2(I+1)*XPD_2(I)) 
        ENDDO

! extra {{{
! particle 1
! {{{

        i = 1
             coef =(cd(i+1,1)+cd(i+1,2)*(12.0*cos(ang(i+1,2))
     1         *cos(ang(i+1,2))-3.0))
     1  *(1.0d0/sqrt(xpd(i+1)*xpd(i)))  

        fta_x(i) = -coef*(-dpd(i+1,2)*dr(i+1,i+2,1) +
     1         dpd(i+1,1)*dr(i+2,i+3,1) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,1) +
     1        dpd(i,2)*dr(i+1,i+2,1))) 

        fta(i,2) = -coef*(-dpd(i+1,2)*dr(i+1,i+2,2) +
     1        dpd(i+1,1)*dr(i+2,i+3,2) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,2) +
     1        dpd(i,2)*dr(i+1,i+2,2))) 

        fta_z(i) = -coef*(-dpd(i+1,2)*dr(i+1,i+2,3) +
     1        dpd(i+1,1)*dr(i+2,i+3,3) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,3) +
     1        dpd(i,2)*dr(i+1,i+2,3))) 

! }}}
! particle 2
! {{{
        i = 2
        coef =(cd(i+1,1)+cd(i+1,2)*(12.0*cos(ang(i+1,2))
     1  *cos(ang(i+1,2)) - 3.0))
     1        *(1.0d0/sqrt(xpd(i+1)*xpd(i)))  

             coef1 = (cd(i,1) + cd(i,2)*(12.0*cos(ang(i,2))
     1        *cos(ang(i,2)) - 
     1        3.0))*(1.0d0/sqrt(xpd(i)*xpd(i-1)))  

        a1 =  -coef*(-dpd(i+1,2)*dr(i+1,i+2,1) +
     1        dpd(i+1,1)*dr(i+2,i+3,1) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,1) +
     1        dpd(i,2)*dr(i+1,i+2,1))) 

        a2 = -coef1*(-dpd(i-1,2)*dr(i+1,i+2,1) +
     1        dpd(i,2)*dr(i,i+1,1) - dpd(i,2)*dr(i-1,i,1) -
     1        dpd(i,1)*dr(i+1,i+2,1) + 2.0*dpd(i-1,3)*dr(i,i+1,1) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,1) -
     1        dpd(i-1,1)*dr(i,i+1,1) - dpd(i-1,2)*dr(i,i+1,1) +
     1        dpd(i-1,2)*dr(i-1,i,1)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,1) +
     1        dpd(i,2)*dr(i+1,i+2,1))) 

        fta_x(i) = a1 + a2 

        a1 = -coef*(-dpd(i+1,2)*dr(i+1,i+2,2) +
     1        dpd(i+1,1)*dr(i+2,i+3,2) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,2) +
     1        dpd(i,2)*dr(i+1,i+2,2))) 

        a2 = -coef1*(-dpd(i-1,2)*dr(i+1,i+2,2) +
     1        dpd(i,2)*dr(i,i+1,2) - dpd(i,2)*dr(i-1,i,2) -
     1        dpd(i,1)*dr(i+1,i+2,2) + 2.0*dpd(i-1,3)*dr(i,i+1,2) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,2) -
     1        dpd(i-1,1)*dr(i,i+1,2) - dpd(i-1,2)*dr(i,i+1,2) +
     1        dpd(i-1,2)*dr(i-1,i,2)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,2) +
     1        dpd(i,2)*dr(i+1,i+2,2)))

        fta(i,2) = a1 + a2 
        
        a1 = -coef*(-dpd(i+1,2)*dr(i+1,i+2,3) +
     1        dpd(i+1,1)*dr(i+2,i+3,3) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,3) +
     1        dpd(i,2)*dr(i+1,i+2,3))) 

        a2 = -coef1*(-dpd(i-1,2)*dr(i+1,i+2,3) +
     1        dpd(i,2)*dr(i,i+1,3) - dpd(i,2)*dr(i-1,i,3) -
     1        dpd(i,1)*dr(i+1,i+2,3) + 2.0*dpd(i-1,3)*dr(i,i+1,3) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,3) -
     1        dpd(i-1,1)*dr(i,i+1,3) - dpd(i-1,2)*dr(i,i+1,3) +
     1        dpd(i-1,2)*dr(i-1,i,3)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,3) +
     1        dpd(i,2)*dr(i+1,i+2,3))) 

        fta_z(i) = a1 + a2 
! }}}
! particle 3
! {{{
        i = 3
        coef=(cd(i+1,1)+cd(i+1,2)*(12.0*cos(ang(i+1,2))
     1  *cos(ang(i+1,2)) - 
     1        3.0))*(1.0d0/sqrt(xpd(i+1)*xpd(i)))  

        coef1=(cd(i,1)+cd(i,2)*(12.0*cos(ang(i,2))
     1        *cos(ang(i,2)) - 
     1        3.0))*(1.0d0/sqrt(xpd(i)*xpd(i-1)))  

        coef2=(cd(i-1,1)+cd(i-1,2)*(12.0*cos(ang(i-1,2))
     1        *cos(ang(i-1,2)) - 
     1        3.0))*(1.0d0/sqrt(xpd(i-1)*xpd(i-2)))  

        a1 = -coef*(-dpd(i+1,2)*dr(i+1,i+2,1) +
     1        dpd(i+1,1)*dr(i+2,i+3,1) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,1) +
     1        dpd(i,2)*dr(i+1,i+2,1))) 

        a2 = -coef1*(-dpd(i-1,2)*dr(i+1,i+2,1) +
     1        dpd(i,2)*dr(i,i+1,1) - dpd(i,2)*dr(i-1,i,1) -
     1        dpd(i,1)*dr(i+1,i+2,1) + 2.0*dpd(i-1,3)*dr(i,i+1,1) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,1) -
     1        dpd(i-1,1)*dr(i,i+1,1) - dpd(i-1,2)*dr(i,i+1,1) +
     1        dpd(i-1,2)*dr(i-1,i,1)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,1) +
     1        dpd(i,2)*dr(i+1,i+2,1))) 

        a3 = -coef2*(dpd(i-2,2)*dr(i,i+1,1) -
     1        dpd(i-2,2)*dr(i-1,i,1) + dpd(i-1,2)*dr(i-2,i-1,1) +
     1        dpd(i-1,1)*dr(i-2,i-1,1) - 2.0*dpd(i-2,3)*dr(i-1,i,1) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,1) -
     1        dpd(i-2,2)*dr(i-2,i-1,1)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,1) -
     1        dpd(i-1,1)*dr(i,i+1,1) - dpd(i-1,2)*dr(i,i+1,1) +
     1        dpd(i-1,2)*dr(i-1,i,1))) 

        fta(i,1) = sum(aa(1:3))
 
        a1 = -coef*(-dpd(i+1,2)*dr(i+1,i+2,2) +
     1        dpd(i+1,1)*dr(i+2,i+3,2) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,2) +
     1        dpd(i,2)*dr(i+1,i+2,2))) 
        
        a2 = -coef1*(-dpd(i-1,2)*dr(i+1,i+2,2) +
     1        dpd(i,2)*dr(i,i+1,2) - dpd(i,2)*dr(i-1,i,2) -
     1        dpd(i,1)*dr(i+1,i+2,2) + 2.0*dpd(i-1,3)*dr(i,i+1,2) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,2) -
     1        dpd(i-1,1)*dr(i,i+1,2) - dpd(i-1,2)*dr(i,i+1,2) +
     1        dpd(i-1,2)*dr(i-1,i,2)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,2) +
     1        dpd(i,2)*dr(i+1,i+2,2))) 

        a3 = -coef2*(dpd(i-2,2)*dr(i,i+1,2) -
     1        dpd(i-2,2)*dr(i-1,i,2) + dpd(i-1,2)*dr(i-2,i-1,2) +
     1        dpd(i-1,1)*dr(i-2,i-1,2) - 2.0*dpd(i-2,3)*dr(i-1,i,2) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,2) -
     1        dpd(i-2,2)*dr(i-2,i-1,2)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,2) -
     1        dpd(i-1,1)*dr(i,i+1,2) - dpd(i-1,2)*dr(i,i+1,2) +
     1        dpd(i-1,2)*dr(i-1,i,2)))

        fta(i,2) = sum(aa(1:3))
 
        a1 = -coef*(-dpd(i+1,2)*dr(i+1,i+2,3) +
     1        dpd(i+1,1)*dr(i+2,i+3,3) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,3) +
     1        dpd(i,2)*dr(i+1,i+2,3))) 

        a2 =  -coef1*(-dpd(i-1,2)*dr(i+1,i+2,3) +
     1        dpd(i,2)*dr(i,i+1,3) - dpd(i,2)*dr(i-1,i,3) -
     1        dpd(i,1)*dr(i+1,i+2,3) + 2.0*dpd(i-1,3)*dr(i,i+1,3) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,3) -
     1        dpd(i-1,1)*dr(i,i+1,3) - dpd(i-1,2)*dr(i,i+1,3) +
     1        dpd(i-1,2)*dr(i-1,i,3)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,3) +
     1        dpd(i,2)*dr(i+1,i+2,3))) 

        a3 = -coef2*(dpd(i-2,2)*dr(i,i+1,3) -
     1        dpd(i-2,2)*dr(i-1,i,3) + dpd(i-1,2)*dr(i-2,i-1,3) +
     1        dpd(i-1,1)*dr(i-2,i-1,3) - 2.0*dpd(i-2,3)*dr(i-1,i,3) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,3) -
     1        dpd(i-2,2)*dr(i-2,i-1,3)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,3) -
     1        dpd(i-1,1)*dr(i,i+1,3) - dpd(i-1,2)*dr(i,i+1,3) +
     1        dpd(i-1,2)*dr(i-1,i,3))) 

        fta(i,3) = sum(aa(1:3))
! }}}
! particle n-2
! {{{
        i = n-2
        coef(1)=(cd(i,1)+cd(i,2)*(12.0*cos(ang(i,2))**2-3.0)*ixpd(i-1)

        coef(2)=(cd(i-1,1)+cd(i-1,2)*(12.0*cos(ang(i-1,2))
     1        *cos(ang(i-1,2)) - 
     1        3.0))*(1.0d0/sqrt(xpd(i-1)*xpd(i-2)))  

        coef3=(cd(i-2,1)+cd(i-2,2)*(12.0*cos(ang(i-2,2))
     1        *cos(ang(i-2,2)) - 
     1  3.0))*(1.0d0/sqrt(xpd(i-2)*xpd(i-3)))  

        a1 = -coef1*(-dpd(i-1,2)*dr(i+1,i+2,1) + 
     1        dpd(i,2)*dr(i,i+1,1) - dpd(i,2)*dr(i-1,i,1) -
     1        dpd(i,1)*dr(i+1,i+2,1) + 2.0*dpd(i-1,3)*dr(i,i+1,1) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,1) -
     1        dpd(i-1,1)*dr(i,i+1,1) - dpd(i-1,2)*dr(i,i+1,1) +
     1        dpd(i-1,2)*dr(i-1,i,1)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,1) +
     1        dpd(i,2)*dr(i+1,i+2,1)))


        a2 = -coef2*(dpd(i-2,2)*dr(i,i+1,1) -
     1        dpd(i-2,2)*dr(i-1,i,1) + dpd(i-1,2)*dr(i-2,i-1,1) +
     1        dpd(i-1,1)*dr(i-2,i-1,1) - 2.0*dpd(i-2,3)*dr(i-1,i,1) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,1) -
     1        dpd(i-2,2)*dr(i-2,i-1,1)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,1) -
     1        dpd(i-1,1)*dr(i,i+1,1) - dpd(i-1,2)*dr(i,i+1,1) +
     1        dpd(i-1,2)*dr(i-1,i,1))) 
        
        a3 = -coef3*(dpd(i-3,2)*dr(i-2,i-1,1) -
     1        dpd(i-2,1)*dr(i-3,i-2,1) -
     1        (1.0d0/xpd(i-2))*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,1) -
     1        dpd(i-2,2)*dr(i-2,i-1,1))) 

        fta(i,1) = sum(aa(1:3)) 

        a1 =  -coef1*(-dpd(i-1,2)*dr(i+1,i+2,2) +  
     1        dpd(i,2)*dr(i,i+1,2) - dpd(i,2)*dr(i-1,i,2) -
     1        dpd(i,1)*dr(i+1,i+2,2) + 2.0*dpd(i-1,3)*dr(i,i+1,2) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,2) -
     1        dpd(i-1,1)*dr(i,i+1,2) - dpd(i-1,2)*dr(i,i+1,2) +
     1        dpd(i-1,2)*dr(i-1,i,2)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,2) +
     1        dpd(i,2)*dr(i+1,i+2,2))) 

        a2 =  -coef2*(dpd(i-2,2)*dr(i,i+1,2) -
     1        dpd(i-2,2)*dr(i-1,i,2) + dpd(i-1,2)*dr(i-2,i-1,2) +
     1        dpd(i-1,1)*dr(i-2,i-1,2) - 2.0*dpd(i-2,3)*dr(i-1,i,2) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,2) -
     1        dpd(i-2,2)*dr(i-2,i-1,2)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,2) -
     1        dpd(i-1,1)*dr(i,i+1,2) - dpd(i-1,2)*dr(i,i+1,2) +
     1        dpd(i-1,2)*dr(i-1,i,2)))

        a3 = -coef3*(dpd(i-3,2)*dr(i-2,i-1,2) -
     1        dpd(i-2,1)*dr(i-3,i-2,2) -
     1        (1.0d0/xpd(i-2))*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,2) -
     1        dpd(i-2,2)*dr(i-2,i-1,2))) 

        fta(i,2) = sum(aa(1:3)) 
 
        aa(1) = -coef1*(-dpd(i-1,2)*dr(i+1,i+2,3) +  
     1        dpd(i,2)*dr(i,i+1,3) - dpd(i,2)*dr(i-1,i,3) -
     1        dpd(i,1)*dr(i+1,i+2,3) + 2.0*dpd(i-1,3)*dr(i,i+1,3) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,3) -
     1        dpd(i-1,1)*dr(i,i+1,3) - dpd(i-1,2)*dr(i,i+1,3) +
     1        dpd(i-1,2)*dr(i-1,i,3)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,3) +
     1        dpd(i,2)*dr(i+1,i+2,3))) 

        aa(2) = -coef2*(dpd(i-2,2)*dr(i,i+1,3) -
     1        dpd(i-2,2)*dr(i-1,i,3) + dpd(i-1,2)*dr(i-2,i-1,3) +
     1        dpd(i-1,1)*dr(i-2,i-1,3) - 2.0*dpd(i-2,3)*dr(i-1,i,3) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,3) -
     1        dpd(i-2,2)*dr(i-2,i-1,3)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,3) -
     1        dpd(i-1,1)*dr(i,i+1,3) - dpd(i-1,2)*dr(i,i+1,3) +
     1        dpd(i-1,2)*dr(i-1,i,3))) 

        aa(3) = -coef3*(dpd(i-3,2)*dr(i-2,i-1,3) -
     1        dpd(i-2,1)*dr(i-3,i-2,3) -
     1        (1.0d0/xpd(i-2))*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,3) -
     1        dpd(i-2,2)*dr(i-2,i-1,3))) 

        fta(i,3) = sum(aa(1:3)) 
! }}}
! particle n-1
! {{{

        i = n-1
        coef2=(cd(i-1,1)+cd(i-1,2)*(12.0*cos(ang(i-1,2))
     1        *cos(ang(i-1,2)) - 
     1        3.0))*(1.0d0/sqrt(xpd(i-1)*xpd(i-2)))  

        coef3=(cd(i-2,1)+cd(i-2,2)*(12.0*cos(ang(i-2,2))
     1        *cos(ang(i-2,2)) - 
     1        3.0))*(1.0d0/sqrt(xpd(i-2)*xpd(i-3)))  

        a1 = -coef2*(dpd(i-2,2)*dr(i,i+1,1) - 
     1        dpd(i-2,2)*dr(i-1,i,1) +
     1        dpd(i-1,2)*dr(i-2,i-1,1) +  dpd(i-1,1)*dr(i-2,i-1,1) -
     1        2.0*dpd(i-2,3)*dr(i-1,i,1) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,1) -
     1        dpd(i-2,2)*dr(i-2,i-1,1)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,1) -
     1        dpd(i-1,1)*dr(i,i+1,1) - dpd(i-1,2)*dr(i,i+1,1) +
     1        dpd(i-1,2)*dr(i-1,i,1))) 

        a2 = -coef3*(dpd(i-3,2)*dr(i-2,i-1,1) - 
     1        dpd(i-2,1)*dr(i-3,i-2,1) -
     1        (1.0d0/xpd(i-2))*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,1) -
     1        dpd(i-2,2)*dr(i-2,i-1,1))) 

        fta(i,1) = sum(aa(1:2))  

        a1 = -coef2*(dpd(i-2,2)*dr(i,i+1,2) - 
     1        dpd(i-2,2)*dr(i-1,i,2) +
     1        dpd(i-1,2)*dr(i-2,i-1,2) +  dpd(i-1,1)*dr(i-2,i-1,2) -
     1        2.0*dpd(i-2,3)*dr(i-1,i,2) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,2) -
     1        dpd(i-2,2)*dr(i-2,i-1,2)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,2) -
     1  dpd(i-1,1)*dr(i,i+1,2) - dpd(i-1,2)*dr(i,i+1,2) +
     1        dpd(i-1,2)*dr(i-1,i,2))) 

        a2 = -coef3*(dpd(i-3,2)*dr(i-2,i-1,2) - 
     1        dpd(i-2,1)*dr(i-3,i-2,2) -
     1        (1.0d0/xpd(i-2))*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,2) -
     1        dpd(i-2,2)*dr(i-2,i-1,2))) 

        fta(i,2) = sum(aa(1:2))  

        a1 = -coef2*(dpd(i-2,2)*dr(i,i+1,3) - 
     1        dpd(i-2,2)*dr(i-1,i,3) +
     1        dpd(i-1,2)*dr(i-2,i-1,3) +  dpd(i-1,1)*dr(i-2,i-1,3) -
     1        2.0*dpd(i-2,3)*dr(i-1,i,3) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,3) -
     1        dpd(i-2,2)*dr(i-2,i-1,3)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,3) -
     1        dpd(i-1,1)*dr(i,i+1,3) - dpd(i-1,2)*dr(i,i+1,3) +
     1        dpd(i-1,2)*dr(i-1,i,3))) 

        a2 = -coef3*(dpd(i-3,2)*dr(i-2,i-1,3) - 
     1        dpd(i-2,1)*dr(i-3,i-2,3) -
     1        (1.0d0/xpd(i-2))*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,3) -
     1        dpd(i-2,2)*dr(i-2,i-1,3))) 

        fta(i,3) = sum(aa(1:2))  
! }}} 
! particle n
! {{{

        i = n
        coef3=(cd(i-2,1)+cd(i-2,2)*(12.0*cos(ang(i-2,2))
     1        *cos(ang(i-2,2)) - 
     1        3.0))*(1.0d0/sqrt(xpd(i-2)*xpd(i-3)))  

        fta(i,1:3) = -coef3*(dpd(i-3,2)*dr(i-2,i-1,1:3) 
     1        - dpd(i-2,1)*dr(i-3,i-2,1:3) -
     1        (1.0d0/xpd(i-2))*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,1) -
     1        dpd(i-2,2)*dr(i-2,i-1,1:3))) 
! }}}
! }}}
! particles 4 to n-3
!
! define bb0 coef {{{
!
        do i = 4, n-3
          bb0(i)=1.0/xpd_2(i)
          do k=1,4
            coef(k) = cd(i-k+2,1) + cd(i-k+2,2)*(12.0*cos(ang(i-k+2,2))**2-3.0)
            coef(k) = coef(k)*ixpd(i+1-k)
          enddo

! }}}
! aa(1) {{{
        ! 2 terms summed 
        aa(1) = -dpd(i+1,2)*bvr(i+1,1) + dpd(i+1,1)*bvr(i+2,1) 

        B1= -dpd(i+1,1)*bvr(i,1) + dpd(i,2)*bvr(i+1,1)
        B2=  dpd(i+1,2)*dpd(i,2) - dpd(i,3)*dpd(i+1,1) 

        aa(1) = aa(1) - bb0(i)*B1*B2
! }}}
! aa(2) {{{
        ! 5 terms summed 
        aa(2) = -dpd(i-1,2)*bvr(i+1,1)  &
                +dpd(i,2)*bvr(i,1)      &
                +2.0*dpd(i-1,3)*bvr(i,1) &
                -dpd(i,2)*bvr(i-1,1)    &
                -dpd(i,1)*bvr(i+1,1)    

        aa(2) = aa(2) &
                - bb0(i-1)&     ! 2  {{{
                 *(                     &
                   dpd(i,2)*dpd(i-1,2)  &
                   - dpd(i-1,3)*dpd(i,1)    &
                  ) &
                 *( 
                    dpd(i,1)*dr(i-1,i,1)  &
                    - dpd(i-1,1)*dr(i,i+1,1) &
                    - dpd(i-1,2)*dr(i,i+1,1) & 
                    + dpd(i-1,2)*dr(i-1,i,1) &
                  ) &           ! }}}
                - bb0(i)&       ! ( 2 )( 2 ) {{{
                 *(                         &
                   dpd(i,2)*dpd(i-1,2)      &
                   - dpd(i-1,3)*dpd(i,1)    &
                  )                         &
                 *(                         &
                   -dpd(i+1,1)*dr(i,i+1,1)  &
                   + dpd(i,2)*dr(i+1,i+2,1) &
                  )
                                ! }}}

! aa(3) {{{
        aa(3) = -dpd(i-2,2)*bvr(i,1) -
     1        dpd(i-2,2)*bvr(i-1,1) + dpd(i-1,2)*bvr(i-2,1) +
     1        dpd(i-1,1)*bvr(i-2,1) - 2.0*dpd(i-2,3)*bvr(i-1,1) -
     1        bb0(i-2)*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*bvr(i-1,1) -
     1        dpd(i-2,2)*bvr(i-2,1)) -
     1        bb0(i-1)*(dpd(i-1,2)*dpd(i-2,2) -
     1  dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*bvr(i-1,1) -
     1        dpd(i-1,1)*bvr(i,1) - dpd(i-1,2)*bvr(i,1) +
     1        dpd(i-1,2)*bvr(i-1,1))
! }}}
! aa(4) {{{
        aa(4) = -dpd(i-3,2)*dr(i-2,i-1,1) -
     1        dpd(i-2,1)*dr(i-3,i-2,1) -
     1        bb0(i-2)*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,1) -
     1        dpd(i-2,2)*dr(i-2,i-1,1))
! }}}

aa=-coef*aa

! }}}
        fta(i,1) = sum(aa)
! {{{
        aa(1) = -coef*(-dpd(i+1,2)*dr(i+1,i+2,2) +
     1        dpd(i+1,1)*dr(i+2,i+3,2) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,2) +
     1        dpd(i,2)*dr(i+1,i+2,2))) 

        aa(2) = -coef1*(-dpd(i-1,2)*dr(i+1,i+2,2) +
     1        dpd(i,2)*dr(i,i+1,2) - dpd(i,2)*dr(i-1,i,2) -
     1        dpd(i,1)*dr(i+1,i+2,2) + 2.0*dpd(i-1,3)*dr(i,i+1,2) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,2) -
     1        dpd(i-1,1)*dr(i,i+1,2) - dpd(i-1,2)*dr(i,i+1,2) +
     1        dpd(i-1,2)*dr(i-1,i,2)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,2) +
     1        dpd(i,2)*dr(i+1,i+2,2))) 

        aa(3) = -coef2*(dpd(i-2,2)*dr(i,i+1,2) -
     1        dpd(i-2,2)*dr(i-1,i,2) + dpd(i-1,2)*dr(i-2,i-1,2) +
     1        dpd(i-1,1)*dr(i-2,i-1,2) - 2.0*dpd(i-2,3)*dr(i-1,i,2) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,2) -
     1        dpd(i-2,2)*dr(i-2,i-1,2)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,2) -
     1        dpd(i-1,1)*dr(i,i+1,2) - dpd(i-1,2)*dr(i,i+1,2) +
     1        dpd(i-1,2)*dr(i-1,i,2))) 

        aa(4) = -coef3*(dpd(i-3,2)*dr(i-2,i-1,2) -
     1        dpd(i-2,1)*dr(i-3,i-2,2) -
     1        (1.0d0/xpd(i-2))*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,2) -
     1        dpd(i-2,2)*dr(i-2,i-1,2))) 
! }}}
        fta(i,2) = sum(aa) 
! {{{

        aa(1) = -coef*(-dpd(i+1,2)*dr(i+1,i+2,3) +
     1        dpd(i+1,1)*dr(i+2,i+3,3) -
     1        (1.0d0/xpd(i))*(dpd(i+1,2)*dpd(i,2) -
     1        dpd(i,3)*dpd(i+1,1))*(-dpd(i+1,1)*dr(i,i+1,3) +
     1        dpd(i,2)*dr(i+1,i+2,3))) 

        aa(2) = -coef1*(-dpd(i-1,2)*dr(i+1,i+2,3) +
     1        dpd(i,2)*dr(i,i+1,3) - dpd(i,2)*dr(i-1,i,3) -
     1        dpd(i,1)*dr(i+1,i+2,3) + 2.0*dpd(i-1,3)*dr(i,i+1,3) -
     1        (1.0d0/xpd(i-1))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(dpd(i,1)*dr(i-1,i,3) -
     1        dpd(i-1,1)*dr(i,i+1,3) - dpd(i-1,2)*dr(i,i+1,3) +
     1        dpd(i-1,2)*dr(i-1,i,3)) -
     1        (1.0d0/xpd(i))*(dpd(i,2)*dpd(i-1,2) -
     1        dpd(i-1,3)*dpd(i,1))*(-dpd(i+1,1)*dr(i,i+1,3) +
     1        dpd(i,2)*dr(i+1,i+2,3))) 

        aa(3) = -coef2*(dpd(i-2,2)*dr(i,i+1,3) -
     1        dpd(i-2,2)*dr(i-1,i,3) + dpd(i-1,2)*dr(i-2,i-1,3) +
     1        dpd(i-1,1)*dr(i-2,i-1,3) - 2.0*dpd(i-2,3)*dr(i-1,i,3) -
     1        (1.0d0/xpd(i-2))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i-2,1)*dr(i-1,i,3) -
     1        dpd(i-2,2)*dr(i-2,i-1,3)) -
     1        (1.0d0/xpd(i-1))*(dpd(i-1,2)*dpd(i-2,2) -
     1        dpd(i-2,3)*dpd(i-1,1))*(dpd(i,1)*dr(i-1,i,3) -
     1        dpd(i-1,1)*dr(i,i+1,3) - dpd(i-1,2)*dr(i,i+1,3) +
     1        dpd(i-1,2)*dr(i-1,i,3))) 
        
        aa(4) = -coef3*(dpd(i-3,2)*dr(i-2,i-1,3) -
     1        dpd(i-2,1)*dr(i-3,i-2,3) -
     1        (1.0d0/xpd(i-2))*(dpd(i-2,2)*dpd(i-3,2) -
     1        dpd(i-3,3)*dpd(i-2,1))*(dpd(i-2,1)*dr(i-1,i,3) -
     1        dpd(i-2,2)*dr(i-2,i-1,3))) 
! }}}
        fta(i,3) = sum(aa)  

        enddo
! }}}


C Torsional angle forces 
C particles 1, 2, 3, n-2, n-1, and n are done outside of the loop
C particle 1
C {{{

        i = 1
             coef =(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1         *dcos(tor_angle(i+1))-3.0))
     1  *(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))  

        fta_x(i) = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1         dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        fta_y(i) = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        fta_z(i) = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 


C particle 2

        i = 2
        coef =(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1  *dcos(tor_angle(i+1)) - 3.0))
     1        *(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))  

             coef1 = (c_param(i) + d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))  

        a1 =  -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1        dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        fta_x(i) = a1 + a2 

        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2)))

        fta_y(i) = a1 + a2 
        
        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        fta_z(i) = a1 + a2 

C particle 3

        i = 3
        coef=(c_param(i+1)+d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1  *dcos(tor_angle(i+1)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))  

        coef1=(c_param(i)+d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))  

        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1        dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 

        fta_x(i) = a1 + a2 + a3 
 
        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 
        
        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)))

        fta_y(i) = a1 + a2 + a3 
 
        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 =  -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 

        fta_z(i) = a1 + a2 + a3 

C particles 4 to n-3

        do i = 4, n-3

        coef = (c_param(i+1) + d_param(i+1)*(12.0*dcos(tor_angle(i+1))
     1  *dcos(tor_angle(i+1)) - 3.0))*(1.0D0/dsqrt(x_prod(i+1)*x_prod(i)))

        coef1 = (c_param(i) + d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) -3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))  

        coef2 = (c_param(i-1) + d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1  3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        coef3 = (c_param(i-2) + d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1  3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        a1 = -coef*(-dot_prod(i+1,2)*xr(i+1,i+2) +
     1        dot_prod(i+1,1)*xr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) +
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1  dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 

        a4 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) -
     1        dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_x(i) = a1 + a2 + a3 + a4 

        a1 = -coef*(-dot_prod(i+1,2)*yr(i+1,i+2) +
     1        dot_prod(i+1,1)*yr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i))) 

        a4 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) -
     1        dot_prod(i-2,1)*yr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_y(i) = a1 + a2 + a3 + a4 

        a1 = -coef*(-dot_prod(i+1,2)*zr(i+1,i+2) +
     1        dot_prod(i+1,1)*zr(i+2,i+3) -
     1        (1.0D0/x_prod(i))*(dot_prod(i+1,2)*dot_prod(i,2) -
     1        dot_prod(i,3)*dot_prod(i+1,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a3 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 
        
        a4 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) -
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

        fta_z(i) = a1 + a2 + a3 + a4 

        enddo

C particle n-2

        i = n-2
        coef1=(c_param(i)+d_param(i)*(12.0*dcos(tor_angle(i))
     1        *dcos(tor_angle(i)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i)*x_prod(i-1)))  

        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1  3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        a1 = -coef1*(-dot_prod(i-1,2)*xr(i+1,i+2) + 
     1        dot_prod(i,2)*xr(i,i+1) - dot_prod(i,2)*xr(i-1,i) -
     1        dot_prod(i,1)*xr(i+1,i+2) + 2.0*dot_prod(i-1,3)*xr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*xr(i,i+1) +
     1        dot_prod(i,2)*xr(i+1,i+2)))


        a2 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) -
     1        dot_prod(i-2,2)*xr(i-1,i) + dot_prod(i-1,2)*xr(i-2,i-1) +
     1        dot_prod(i-1,1)*xr(i-2,i-1) - 2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 
        
        a3 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) -
     1        dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_x(i) = a1 + a2 + a3 

        a1 =  -coef1*(-dot_prod(i-1,2)*yr(i+1,i+2) +  
     1        dot_prod(i,2)*yr(i,i+1) - dot_prod(i,2)*yr(i-1,i) -
     1        dot_prod(i,1)*yr(i+1,i+2) + 2.0*dot_prod(i-1,3)*yr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*yr(i,i+1) +
     1        dot_prod(i,2)*yr(i+1,i+2))) 

        a2 =  -coef2*(dot_prod(i-2,2)*yr(i,i+1) -
     1        dot_prod(i-2,2)*yr(i-1,i) + dot_prod(i-1,2)*yr(i-2,i-1) +
     1        dot_prod(i-1,1)*yr(i-2,i-1) - 2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1        dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i)))

        a3 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) -
     1        dot_prod(i-2,1)*yr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_y(i) = a1 + a2 + a3 
 
        a1 = -coef1*(-dot_prod(i-1,2)*zr(i+1,i+2) +  
     1        dot_prod(i,2)*zr(i,i+1) - dot_prod(i,2)*zr(i-1,i) -
     1        dot_prod(i,1)*zr(i+1,i+2) + 2.0*dot_prod(i-1,3)*zr(i,i+1) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i)) -
     1        (1.0D0/x_prod(i))*(dot_prod(i,2)*dot_prod(i-1,2) -
     1        dot_prod(i-1,3)*dot_prod(i,1))*(-dot_prod(i+1,1)*zr(i,i+1) +
     1        dot_prod(i,2)*zr(i+1,i+2))) 

        a2 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) -
     1        dot_prod(i-2,2)*zr(i-1,i) + dot_prod(i-1,2)*zr(i-2,i-1) +
     1        dot_prod(i-1,1)*zr(i-2,i-1) - 2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 

        a3 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) -
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

        fta_z(i) = a1 + a2 + a3 

C particle n-1

        i = n-1
        coef2=(c_param(i-1)+d_param(i-1)*(12.0*dcos(tor_angle(i-1))
     1        *dcos(tor_angle(i-1)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-1)*x_prod(i-2)))  

        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        a1 = -coef2*(dot_prod(i-2,2)*xr(i,i+1) - 
     1        dot_prod(i-2,2)*xr(i-1,i) +
     1        dot_prod(i-1,2)*xr(i-2,i-1) +  dot_prod(i-1,1)*xr(i-2,i-1) -
     1        2.0*dot_prod(i-2,3)*xr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*xr(i-1,i) -
     1        dot_prod(i-1,1)*xr(i,i+1) - dot_prod(i-1,2)*xr(i,i+1) +
     1        dot_prod(i-1,2)*xr(i-1,i))) 

        a2 = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) - 
     1        dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_x(i) = a1 + a2  

        a1 = -coef2*(dot_prod(i-2,2)*yr(i,i+1) - 
     1        dot_prod(i-2,2)*yr(i-1,i) +
     1        dot_prod(i-1,2)*yr(i-2,i-1) +  dot_prod(i-1,1)*yr(i-2,i-1) -
     1        2.0*dot_prod(i-2,3)*yr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*yr(i-1,i) -
     1  dot_prod(i-1,1)*yr(i,i+1) - dot_prod(i-1,2)*yr(i,i+1) +
     1        dot_prod(i-1,2)*yr(i-1,i))) 

        a2 = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - 
     1        dot_prod(i-2,1)*yr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_y(i) = a1 + a2  

        a1 = -coef2*(dot_prod(i-2,2)*zr(i,i+1) - 
     1        dot_prod(i-2,2)*zr(i-1,i) +
     1        dot_prod(i-1,2)*zr(i-2,i-1) +  dot_prod(i-1,1)*zr(i-2,i-1) -
     1        2.0*dot_prod(i-2,3)*zr(i-1,i) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1)) -
     1        (1.0D0/x_prod(i-1))*(dot_prod(i-1,2)*dot_prod(i-2,2) -
     1        dot_prod(i-2,3)*dot_prod(i-1,1))*(dot_prod(i,1)*zr(i-1,i) -
     1        dot_prod(i-1,1)*zr(i,i+1) - dot_prod(i-1,2)*zr(i,i+1) +
     1        dot_prod(i-1,2)*zr(i-1,i))) 

        a2 = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - 
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

        fta_z(i) = a1 + a2 
 
C particle n

        i = n
        coef3=(c_param(i-2)+d_param(i-2)*(12.0*dcos(tor_angle(i-2))
     1        *dcos(tor_angle(i-2)) - 
     1        3.0))*(1.0D0/dsqrt(x_prod(i-2)*x_prod(i-3)))  

        fta_x(i) = -coef3*(dot_prod(i-3,2)*xr(i-2,i-1) 
     1        - dot_prod(i-2,1)*xr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*xr(i-1,i) -
     1        dot_prod(i-2,2)*xr(i-2,i-1))) 

        fta_y(i) = -coef3*(dot_prod(i-3,2)*yr(i-2,i-1) - 
     1        dot_prod(i-2,1)*yr(i-3,i-2) 
     1        - (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) 
     1        - dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*yr(i-1,i) -
     1        dot_prod(i-2,2)*yr(i-2,i-1))) 

        fta_z(i) = -coef3*(dot_prod(i-3,2)*zr(i-2,i-1) - 
     1        dot_prod(i-2,1)*zr(i-3,i-2) -
     1        (1.0D0/x_prod(i-2))*(dot_prod(i-2,2)*dot_prod(i-3,2) -
     1        dot_prod(i-3,3)*dot_prod(i-2,1))*(dot_prod(i-2,1)*zr(i-1,i) -
     1        dot_prod(i-2,2)*zr(i-2,i-1))) 

C}}}


