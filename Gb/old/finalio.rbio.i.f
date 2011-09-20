C
C  Subroutine to convert rigid body CofM and DV coordinates to molecular sites.
C
      SUBROUTINE RBIO(X1, Y1, Z1, L1, M1, N1, COORDS, NRBSITES, SITE)
      ! {{{
      IMPLICIT NONE
      INTEGER NRBSITES
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, SITE(NRBSITES,3),
     2                 M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12
      INTEGER J1

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
   
      DO J1=1,NRBSITES
         COORDS(3*(J1-1)+1)=c2a1*SITE(J1,1) + s1*(n1*SITE(J1,2) - m1*SITE(J1,3)) - 
     1                      c3a1*l1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + X1
         COORDS(3*(J1-1)+2)=c2a1*SITE(J1,2) + s1*(-(n1*SITE(J1,1)) + l1*SITE(J1,3)) 
     1                    - c3a1*m1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + Y1
         COORDS(3*(J1-1)+3)=s1*(m1*SITE(J1,1) - l1*SITE(J1,2)) + c2a1*SITE(J1,3) 
     1                    - c3a1*n1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + Z1
      ENDDO 

      RETURN
      ! }}}
      END

