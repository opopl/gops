
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

