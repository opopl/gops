C  SUBROUTINE to convert TIP oxygen and DV coordinates to Cartesians.
C
      SUBROUTINE TIPIO(X1, Y1, Z1, L1, M1, N1, COORDS)
      ! {{{
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=(-ALPHA1/2+ALPHA1**3/24)
         C3A1=(-0.5D0+ALPHA1**2/24.0D0)
         S1=(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=(CA1-1.0D0)/ALPHA1**2
         S1=SIN(ALPHA1)/ALPHA1
      ENDIF

      COORDS(1) = X1
      COORDS(2) = Y1
      COORDS(3) = Z1    
      COORDS(4) = 0.756950327*c2a1 - c3a1*l1*(0.756950327*l1 - 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(5) = -(c3a1*m1*(0.756950327*l1 - 0.585882276*n1)) + (-0.585882276*l1 - 0.756950327*n1)*s1 + Y1
      COORDS(6) = -0.585882276*c2a1 - c3a1*(0.756950327*l1 - 0.585882276*n1)*n1 + 0.756950327*m1*s1 + Z1
      COORDS(7) = -0.756950327*c2a1 + c3a1*l1*(0.756950327*l1 + 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(8) = c3a1*m1*(0.756950327*l1 + 0.585882276*n1) + (-0.585882276*l1 + 0.756950327*n1)*s1 + Y1
      COORDS(9) = -0.585882276*c2a1 + c3a1*(0.756950327*l1 + 0.585882276*n1)*n1 - 0.756950327*m1*s1 + Z1

      RETURN
      ! }}}
      END

