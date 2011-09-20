
C  SUBROUTINE to convert capsid CofM and DV coordinates to penatgons.
C
      SUBROUTINE CAPSIDIO(X1, Y1, Z1, L1, M1, N1,COORDS,RAD,HEIGHT)
      ! {{{
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), HEIGHT, C2A1,
     2                 M1, L1, N1, ALPHA1, RAD, CA1, S1, C3A1,
     3                 NUM1, NUM2, NUM3, NUM4, NUM5, L12, M12, N12

      NUM1=-(1.0D0+SQRT(5.0D0))/4.0D0
      NUM2=SQRT((5.0D0-SQRT(5.0D0))/2.0D0)/2.0D0
      NUM3=SQRT((5.0D0+SQRT(5.0D0))/2.0D0)/2.0D0
      NUM4=(SQRT(5.0D0)-1.0D0)/4.0D0
      NUM5=-(1.0D0+SQRT(5.0D0))/4.0D0

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=RAD*CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=RAD*(-ALPHA1/2+ALPHA1**3/24)
         C3A1=RAD*(-0.5D0+ALPHA1**2/24.0D0)
         S1=RAD*(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=RAD*(CA1-1.0D0)/ALPHA1**2
         S1=RAD*SIN(ALPHA1)/ALPHA1
      ENDIF

      COORDS(1) =     c2a1 - c3a1*l12 + x1
      COORDS(2) =     -(c3a1*l1*m1) - n1*s1 + y1
      COORDS(3) =     -(c3a1*l1*n1) + m1*s1 + z1
      COORDS(4) =     c2a1*num4 - c3a1*l1*(m1*num3 + l1*num4) + n1*num3*s1 + x1
      COORDS(5) =     c2a1*num3 - c3a1*m1*(m1*num3 + l1*num4) - n1*num4*s1 + y1
      COORDS(6) =     -(c3a1*n1*(m1*num3 + l1*num4)) - l1*num3*s1 + m1*num4*s1 + z1
      COORDS(7) =     c2a1*num1 - c3a1*l1*(l1*num1 + m1*num2) + n1*num2*s1 + x1
      COORDS(8) = c2a1*num2 - c3a1*m1*(l1*num1 + m1*num2) - n1*num5*s1 + y1
      COORDS(9) = -(c3a1*n1*(l1*num1 + m1*num2)) + m1*num1*s1 - l1*num2*s1 + z1
      COORDS(10) = c2a1*num1 + c3a1*l1*(-(l1*num1) + m1*num2) - n1*num2*s1 + x1
      COORDS(11) = -(c2a1*num2) + c3a1*m1*(-(l1*num1) + m1*num2) - n1*num5*s1 + y1
      COORDS(12) = -(c3a1*l1*n1*num1) + c3a1*m1*n1*num2 + m1*num1*s1 + l1*num2*s1 + z1
      COORDS(13) = c2a1*num4 + c3a1*l1*(m1*num3 - l1*num4) - n1*num3*s1 + x1
      COORDS(14) = -(c2a1*num3) + c3a1*m1*(m1*num3 - l1*num4) - n1*num4*s1 + y1
      COORDS(15) = c3a1*n1*(m1*num3 - l1*num4) + l1*num3*s1 + m1*num4*s1 + z1
C     COORDS(16)= (-(c3a1*l1*n1) - m1*s1 + 2*x1)/2.
C     COORDS(17)= -(c3a1*m1*n1)/2. + (l1*s1)/2. + y1
C     COORDS(18)= (c2a1 - c3a1*n12 + 2*z1)/2.
      COORDS(16)= -(c3a1*height*l1*n1) - height*m1*s1 + x1
      COORDS(17)= -(c3a1*height*m1*n1) + height*l1*s1 + y1
      COORDS(18)= c2a1*height - c3a1*height*n12 + z1

      RETURN
      ! }}}
      END

