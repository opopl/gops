
! Torsional angle forces
! particles 1, 2, 3, n-2, n-1, and n are done outside of the loop

        DO I=1,N
           IXPD(I)=1.0D0/SQRT(XPD_2(I+1)*XPD_2(I)) 
        ENDDO

! particle 1
! {{{

        I = 1
             COEF =(CD(I+1,1)+CD(I+1,2)*(12.0*COS(ANG(I+1,2))
     1         *COS(ANG(I+1,2))-3.0))
     1  *(1.0D0/SQRT(XPD(I+1)*XPD(I)))  

        FTA_X(I) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,1) +
     1         DPD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        FTA(I,2) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,2) +
     1        DPD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        FTA_Z(I) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,3) +
     1        DPD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

! }}}
! PARTICLE 2
! {{{
        I = 2
        COEF =(CD(I+1,1)+CD(I+1,2)*(12.0*COS(ANG(I+1,2))
     1  *COS(ANG(I+1,2)) - 3.0))
     1        *(1.0D0/SQRT(XPD(I+1)*XPD(I)))  

             COEF1 = (CD(I,1) + CD(I,2)*(12.0*COS(ANG(I,2))
     1        *COS(ANG(I,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I)*XPD(I-1)))  

        A1 =  -COEF*(-DPD(I+1,2)*DR(I+1,I+2,1) +
     1        DPD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,1) +
     1        DPD(I,2)*DR(I,I+1,1) - DPD(I,2)*DR(I-1,I,1) -
     1        DPD(I,1)*DR(I+1,I+2,1) + 2.0*DPD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        FTA_X(I) = A1 + A2 

        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,2) +
     1        DPD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,2) +
     1        DPD(I,2)*DR(I,I+1,2) - DPD(I,2)*DR(I-1,I,2) -
     1        DPD(I,1)*DR(I+1,I+2,2) + 2.0*DPD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2)))

        FTA(I,2) = A1 + A2 
        
        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,3) +
     1        DPD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,3) +
     1        DPD(I,2)*DR(I,I+1,3) - DPD(I,2)*DR(I-1,I,3) -
     1        DPD(I,1)*DR(I+1,I+2,3) + 2.0*DPD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        FTA_Z(I) = A1 + A2 
! }}}
! PARTICLE 3
! {{{
        I = 3
        COEF=(CD(I+1,1)+CD(I+1,2)*(12.0*COS(ANG(I+1,2))
     1  *COS(ANG(I+1,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I+1)*XPD(I)))  

        COEF1=(CD(I,1)+CD(I,2)*(12.0*COS(ANG(I,2))
     1        *COS(ANG(I,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I)*XPD(I-1)))  

        COEF2=(CD(I-1,1)+CD(I-1,2)*(12.0*COS(ANG(I-1,2))
     1        *COS(ANG(I-1,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-1)*XPD(I-2)))  

        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,1) +
     1        DPD(I+1,1)*DR(I+2,I+3,1) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,1) +
     1        DPD(I,2)*DR(I,I+1,1) - DPD(I,2)*DR(I-1,I,1) -
     1        DPD(I,1)*DR(I+1,I+2,1) + 2.0*DPD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        A3 = -COEF2*(DPD(I-2,2)*DR(I,I+1,1) -
     1        DPD(I-2,2)*DR(I-1,I,1) + DPD(I-1,2)*DR(I-2,I-1,1) +
     1        DPD(I-1,1)*DR(I-2,I-1,1) - 2.0*DPD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1))) 

        FTA(I,1) = SUM(AA(1:3))
 
        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,2) +
     1        DPD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 
        
        A2 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,2) +
     1        DPD(I,2)*DR(I,I+1,2) - DPD(I,2)*DR(I-1,I,2) -
     1        DPD(I,1)*DR(I+1,I+2,2) + 2.0*DPD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        A3 = -COEF2*(DPD(I-2,2)*DR(I,I+1,2) -
     1        DPD(I-2,2)*DR(I-1,I,2) + DPD(I-1,2)*DR(I-2,I-1,2) +
     1        DPD(I-1,1)*DR(I-2,I-1,2) - 2.0*DPD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)))

        FTA(I,2) = SUM(AA(1:3))
 
        A1 = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,3) +
     1        DPD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        A2 =  -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,3) +
     1        DPD(I,2)*DR(I,I+1,3) - DPD(I,2)*DR(I-1,I,3) -
     1        DPD(I,1)*DR(I+1,I+2,3) + 2.0*DPD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        A3 = -COEF2*(DPD(I-2,2)*DR(I,I+1,3) -
     1        DPD(I-2,2)*DR(I-1,I,3) + DPD(I-1,2)*DR(I-2,I-1,3) +
     1        DPD(I-1,1)*DR(I-2,I-1,3) - 2.0*DPD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3))) 

        FTA(I,3) = SUM(AA(1:3))
! }}}
! PARTICLES 4 TO N-3
! {{{

        DO I = 4, N-3
          DO K=1,4
            COEF(K) = CD(I-K+2,1) + CD(I-K+2,2)*(12.0*COS(ANG(I-K+2,2))**2-3.0)
            COEF(K) = COEF(K)*IXPD(I+1-K)
          ENDDO

        AA(1) = -DPD(I+1,2)*BVR(I+1,1)+DPD(I+1,1)*BVR(I+2,1) 

        aa0=1.0D0/xpd_2(i)
        aa1= -dpd(i+1,1)*bvr(i,1) + dpd(i,2)*bvr(i+1,1)

        AA(1) = AA(1) - (1.0D0/XPD_2(I))*(DPD(I+1,2)*DPD(I,2)-DPD(I,3)*DPD(I+1,1))*AA0*AA1

AA(1)=-COEF(1)*AA(1)

        AA(2) = -COEF(2)*(-DPD(I-1,2)*DR(I+1,I+2,1) +
     1        DPD(I,2)*DR(I,I+1,1) - DPD(I,2)*DR(I-1,I,1) -
     1        DPD(I,1)*DR(I+1,I+2,1) + 2.0*DPD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1))) 

        AA(3) = -COEF(3)*(DPD(I-2,2)*DR(I,I+1,1) -
     1        DPD(I-2,2)*DR(I-1,I,1) + DPD(I-1,2)*DR(I-2,I-1,1) +
     1        DPD(I-1,1)*DR(I-2,I-1,1) - 2.0*DPD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1  DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1))) 

        AA(4) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,1) -
     1        DPD(I-2,1)*DR(I-3,I-2,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1))) 

        FTA(I,1) =sum(AA)

        AA(1) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,2) +
     1        DPD(I+1,1)*DR(I+2,I+3,2) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        AA(2) = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,2) +
     1        DPD(I,2)*DR(I,I+1,2) - DPD(I,2)*DR(I-1,I,2) -
     1        DPD(I,1)*DR(I+1,I+2,2) + 2.0*DPD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        AA(3) = -COEF2*(DPD(I-2,2)*DR(I,I+1,2) -
     1        DPD(I-2,2)*DR(I-1,I,2) + DPD(I-1,2)*DR(I-2,I-1,2) +
     1        DPD(I-1,1)*DR(I-2,I-1,2) - 2.0*DPD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2))) 

        AA(4) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,2) -
     1        DPD(I-2,1)*DR(I-3,I-2,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2))) 

        FTA(I,2) = sum(AA) 

        AA(1) = -COEF*(-DPD(I+1,2)*DR(I+1,I+2,3) +
     1        DPD(I+1,1)*DR(I+2,I+3,3) -
     1        (1.0D0/XPD(I))*(DPD(I+1,2)*DPD(I,2) -
     1        DPD(I,3)*DPD(I+1,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        AA(2) = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,3) +
     1        DPD(I,2)*DR(I,I+1,3) - DPD(I,2)*DR(I-1,I,3) -
     1        DPD(I,1)*DR(I+1,I+2,3) + 2.0*DPD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        AA(3) = -COEF2*(DPD(I-2,2)*DR(I,I+1,3) -
     1        DPD(I-2,2)*DR(I-1,I,3) + DPD(I-1,2)*DR(I-2,I-1,3) +
     1        DPD(I-1,1)*DR(I-2,I-1,3) - 2.0*DPD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3))) 
        
        AA(4) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,3) -
     1        DPD(I-2,1)*DR(I-3,I-2,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3))) 

        FTA(I,3) = SUM(AA)  

        ENDDO
! }}}
! PARTICLE N-2
! {{{
        I = N-2
        COEF(1)=(CD(I,1)+CD(I,2)*(12.0*COS(ANG(I,2))**2-3.0)*IXPD(I-1)

        COEF(2)=(CD(I-1,1)+CD(I-1,2)*(12.0*COS(ANG(I-1,2))
     1        *COS(ANG(I-1,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-1)*XPD(I-2)))  

        COEF3=(CD(I-2,1)+CD(I-2,2)*(12.0*COS(ANG(I-2,2))
     1        *COS(ANG(I-2,2)) - 
     1  3.0))*(1.0D0/SQRT(XPD(I-2)*XPD(I-3)))  

        A1 = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,1) + 
     1        DPD(I,2)*DR(I,I+1,1) - DPD(I,2)*DR(I-1,I,1) -
     1        DPD(I,1)*DR(I+1,I+2,1) + 2.0*DPD(I-1,3)*DR(I,I+1,1) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,1) +
     1        DPD(I,2)*DR(I+1,I+2,1)))


        A2 = -COEF2*(DPD(I-2,2)*DR(I,I+1,1) -
     1        DPD(I-2,2)*DR(I-1,I,1) + DPD(I-1,2)*DR(I-2,I-1,1) +
     1        DPD(I-1,1)*DR(I-2,I-1,1) - 2.0*DPD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1))) 
        
        A3 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,1) -
     1        DPD(I-2,1)*DR(I-3,I-2,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1))) 

        FTA(I,1) = sum(AA(1:3)) 

        A1 =  -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,2) +  
     1        DPD(I,2)*DR(I,I+1,2) - DPD(I,2)*DR(I-1,I,2) -
     1        DPD(I,1)*DR(I+1,I+2,2) + 2.0*DPD(I-1,3)*DR(I,I+1,2) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,2) +
     1        DPD(I,2)*DR(I+1,I+2,2))) 

        A2 =  -COEF2*(DPD(I-2,2)*DR(I,I+1,2) -
     1        DPD(I-2,2)*DR(I-1,I,2) + DPD(I-1,2)*DR(I-2,I-1,2) +
     1        DPD(I-1,1)*DR(I-2,I-1,2) - 2.0*DPD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,2) -
     1        DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2)))

        A3 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,2) -
     1        DPD(I-2,1)*DR(I-3,I-2,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2))) 

        FTA(I,2) = sum(AA(1:3)) 
 
        AA(1) = -COEF1*(-DPD(I-1,2)*DR(I+1,I+2,3) +  
     1        DPD(I,2)*DR(I,I+1,3) - DPD(I,2)*DR(I-1,I,3) -
     1        DPD(I,1)*DR(I+1,I+2,3) + 2.0*DPD(I-1,3)*DR(I,I+1,3) -
     1        (1.0D0/XPD(I-1))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3)) -
     1        (1.0D0/XPD(I))*(DPD(I,2)*DPD(I-1,2) -
     1        DPD(I-1,3)*DPD(I,1))*(-DPD(I+1,1)*DR(I,I+1,3) +
     1        DPD(I,2)*DR(I+1,I+2,3))) 

        AA(2) = -COEF2*(DPD(I-2,2)*DR(I,I+1,3) -
     1        DPD(I-2,2)*DR(I-1,I,3) + DPD(I-1,2)*DR(I-2,I-1,3) +
     1        DPD(I-1,1)*DR(I-2,I-1,3) - 2.0*DPD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3))) 

        AA(3) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,3) -
     1        DPD(I-2,1)*DR(I-3,I-2,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3))) 

        FTA(I,3) = sum(AA(1:3)) 
! }}}
! PARTICLE N-1
! {{{

        I = N-1
        COEF2=(CD(I-1,1)+CD(I-1,2)*(12.0*COS(ANG(I-1,2))
     1        *COS(ANG(I-1,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-1)*XPD(I-2)))  

        COEF3=(CD(I-2,1)+CD(I-2,2)*(12.0*COS(ANG(I-2,2))
     1        *COS(ANG(I-2,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-2)*XPD(I-3)))  

        A1 = -COEF2*(DPD(I-2,2)*DR(I,I+1,1) - 
     1        DPD(I-2,2)*DR(I-1,I,1) +
     1        DPD(I-1,2)*DR(I-2,I-1,1) +  DPD(I-1,1)*DR(I-2,I-1,1) -
     1        2.0*DPD(I-2,3)*DR(I-1,I,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,1) -
     1        DPD(I-1,1)*DR(I,I+1,1) - DPD(I-1,2)*DR(I,I+1,1) +
     1        DPD(I-1,2)*DR(I-1,I,1))) 

        A2 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,1) - 
     1        DPD(I-2,1)*DR(I-3,I-2,1) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1))) 

        FTA(I,1) = sum(AA(1:2))  

        A1 = -COEF2*(DPD(I-2,2)*DR(I,I+1,2) - 
     1        DPD(I-2,2)*DR(I-1,I,2) +
     1        DPD(I-1,2)*DR(I-2,I-1,2) +  DPD(I-1,1)*DR(I-2,I-1,2) -
     1        2.0*DPD(I-2,3)*DR(I-1,I,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,2) -
     1  DPD(I-1,1)*DR(I,I+1,2) - DPD(I-1,2)*DR(I,I+1,2) +
     1        DPD(I-1,2)*DR(I-1,I,2))) 

        A2 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,2) - 
     1        DPD(I-2,1)*DR(I-3,I-2,2) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,2) -
     1        DPD(I-2,2)*DR(I-2,I-1,2))) 

        FTA(I,2) = sum(AA(1:2))  

        A1 = -COEF2*(DPD(I-2,2)*DR(I,I+1,3) - 
     1        DPD(I-2,2)*DR(I-1,I,3) +
     1        DPD(I-1,2)*DR(I-2,I-1,3) +  DPD(I-1,1)*DR(I-2,I-1,3) -
     1        2.0*DPD(I-2,3)*DR(I-1,I,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3)) -
     1        (1.0D0/XPD(I-1))*(DPD(I-1,2)*DPD(I-2,2) -
     1        DPD(I-2,3)*DPD(I-1,1))*(DPD(I,1)*DR(I-1,I,3) -
     1        DPD(I-1,1)*DR(I,I+1,3) - DPD(I-1,2)*DR(I,I+1,3) +
     1        DPD(I-1,2)*DR(I-1,I,3))) 

        A2 = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,3) - 
     1        DPD(I-2,1)*DR(I-3,I-2,3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,3) -
     1        DPD(I-2,2)*DR(I-2,I-1,3))) 

        FTA(I,3) = SUM(AA(1:2))  
! }}} 
! PARTICLE N
! {{{

        I = N
        COEF3=(CD(I-2,1)+CD(I-2,2)*(12.0*COS(ANG(I-2,2))
     1        *COS(ANG(I-2,2)) - 
     1        3.0))*(1.0D0/SQRT(XPD(I-2)*XPD(I-3)))  

        FTA(I,1:3) = -COEF3*(DPD(I-3,2)*DR(I-2,I-1,1:3) 
     1        - DPD(I-2,1)*DR(I-3,I-2,1:3) -
     1        (1.0D0/XPD(I-2))*(DPD(I-2,2)*DPD(I-3,2) -
     1        DPD(I-3,3)*DPD(I-2,1))*(DPD(I-2,1)*DR(I-1,I,1) -
     1        DPD(I-2,2)*DR(I-2,I-1,1:3))) 
! }}}

