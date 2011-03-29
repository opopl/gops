C     ROUTINE TO CALCULATION FORCES AND ENERGY.
C     X AND Y ARE PARTICLE COORDINATES.
C     D1 IS THE DIAMETER OF SMALL PARTICLES.
C     D IS THE DIAMETER OF PARTICLES IN THE UNIT OF SMALL PARTICLE DIAMETER,
C     SO D=1 FOR PARTICLES 1 TO N / 2 AND D=1.4 FOR PARTICLES N/2+1 TO N.
C     N IS THE NUMBER OF PARTICLES.  THE PARTICLES ARE BIDISPERSE WITH DIAMETER
C     RATIO OF 1.4.
C     FX AND FY ARE FORCES IN X AND Y DIRECTION ACTING ON PARTICLES.
C     V IS THE TOTAL POTENTIAL ENERGY.


      SUBROUTINE DF1GRAD(COORDS,N,VNEW,V,GTEST,BOXLX,BOXLY)
      USE COMMONS,ONLY : FIXIMAGE
      IMPLICIT NONE
      
      INTEGER N, J1, I, J, J2
      LOGICAL GTEST
      DOUBLE PRECISION X(N), Y(N), COORDS(3*N), VNEW(3*N), BOXLX, BOXLY
      DOUBLE PRECISION D(N), D1
      DOUBLE PRECISION FX(N), FY(N), V, FR
      DOUBLE PRECISION XIJ, YIJ, RIJ, DIJ

      D1=1.4D0
      DO J1=1,N/2
         D(J1)=1.0D0
      ENDDO
      DO J1=N/2+1,N
         D(J1)=D1
      ENDDO
C
C  DEAL WITH ANY ATOMS THAT HAVE LEFT THE BOX.
C
C     IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
      IF (.NOT.FIXIMAGE) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
            COORDS(J2+1)=COORDS(J2+1) - BOXLX*ANINT(COORDS(J2+1)/BOXLX)
            COORDS(J2+2)=COORDS(J2+2) - BOXLY*ANINT(COORDS(J2+2)/BOXLY)
         ENDDO
      ENDIF

      DO J1=1,N
         X(J1)=COORDS(3*(J1-1)+1)
         Y(J1)=COORDS(3*(J1-1)+2)
      ENDDO

      V = 0.D0

      DO I = 1, N
         FX(I) = 0.D0
         FY(I) = 0.D0
      ENDDO
      
      DO I = 1, N - 1
         DO J = I + 1, N
            XIJ = X(I) - X(J)
            XIJ = XIJ - BOXLX*NINT(XIJ/BOXLX)
            YIJ = Y(I) - Y(J)
            YIJ = YIJ - BOXLY*NINT(YIJ/BOXLY)
            RIJ = DSQRT(XIJ * XIJ + YIJ * YIJ)
            DIJ = (D(I) + D(J)) / 2.D0
            IF (RIJ .LT. DIJ) THEN
               FR = (1.D0 - RIJ / DIJ) / DIJ
               FX(I) = FX(I) + FR * XIJ / RIJ
               FX(J) = FX(J) - FR * XIJ / RIJ
               FY(I) = FY(I) + FR * YIJ / RIJ
               FY(J) = FY(J) - FR * YIJ / RIJ
              
               V = V + (1.D0 - RIJ / DIJ) ** 2 / 2.D0
            ENDIF
         ENDDO
      ENDDO
      DO J1=1,N
         VNEW(3*(J1-1)+1)=-FX(J1)
         VNEW(3*(J1-1)+2)=-FY(J1)
      ENDDO

      RETURN
      END
