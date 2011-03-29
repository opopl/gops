C     ROUTINE TO CALCULATION FORCES AND ENERGY.
C     X AND Y ARE PARTICLE COORDINATES.
C     D1 IS THE DIAMETER OF SMALL PARTICLES.
C     D IS THE DIAMETER OF PARTICLES IN THE UNIT OF SMALL PARTICLE DIAMETER,
C     SO D=1 FOR PARTICLES 1 TO N / 2 AND D=1.4 FOR PARTICLES N/2+1 TO N.
C     N IS THE NUMBER OF PARTICLES.  THE PARTICLES ARE BIDISPERSE WITH DIAMETER
C     RATIO OF 1.4.
C     FX AND FY ARE FORCES IN X AND Y DIRECTION ACTING ON PARTICLES.
C     V IS THE TOTAL POTENTIAL ENERGY.


      SUBROUTINE DF1GRAD(COORDS,N,VNEW,V,GTEST,SSTEST,BOXLX,BOXLY)
      USE KEY
      IMPLICIT NONE
      
      INTEGER N, J1, I, J, J2
      LOGICAL GTEST, SSTEST
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
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
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

      IF (.NOT.SSTEST) RETURN
      CALL DF1HESS(X,Y,D1,D,BOXLX,BOXLY,N)
 
      RETURN
      END

C     ROUTINE TO CALCULATION HESSIAN MATRIX.
C     X AND Y ARE PARTICLE COORDINATES.
C     D1 IS THE DIAMETER OF SMALL PARTICLES.
C     D IS THE DIAMETER OF PARTICLES IN THE UNIT OF SMALL PARTICLE DIAMETER,
C     SO D=1 FOR PARTICLES 1 TO N / 2 AND D=1.4 FOR PARTICLES N/2+1 TO N.
C     N IS THE NUMBER OF PARTICLES.  THE PARTICLES ARE BIDISPERSE WITH DIAMETER
C     RATIO OF 1.4.
C     M IS THE RETURNED HESSIAN MATRIX.

      SUBROUTINE DF1HESS(X, Y, D1, D, BOXLX,BOXLY,N)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(N), Y(N), BOXLX, BOXLY
      DOUBLE PRECISION D(N), D1
      DOUBLE PRECISION XIJ, YIJ, RIJ, DIJ, TIJ, CIJ, COE
      INTEGER I1, I2, J1, J2, I, J

      HESS(1:3*N,1:3*N)=0.0D0
      
      DO I = 1, N - 1
         I1 = 3 * I - 2  ! X
         I2 = I1 + 1     ! Y
         DO J = I + 1, N
            XIJ = X(I) - X(J)
            XIJ = XIJ - BOXLX*NINT(XIJ/BOXLX)
            YIJ = Y(I) - Y(J)
            YIJ = YIJ - BOXLY*NINT(YIJ/BOXLY)
            RIJ = DSQRT(XIJ * XIJ + YIJ * YIJ)
            DIJ = (D(I) + D(J)) / 2.D0
            IF(RIJ .LT. DIJ) THEN
               TIJ = - (1.D0 - RIJ / DIJ) / DIJ
               CIJ = 1.D0 / DIJ ** 2

               COE = - CIJ / RIJ ** 2 + TIJ / RIJ ** 3
               
               J1 = 3 * J - 2 ! X
               J2 = J1 + 1    ! Y
               
               HESS(I1, J1) = COE * XIJ * XIJ - TIJ / RIJ
               HESS(I1, J2) = COE * XIJ * YIJ
               HESS(I2, J1) = HESS(I1, J2)
               HESS(I2, J2) = COE * YIJ * YIJ - TIJ / RIJ
               
               HESS(I1, I1) = HESS(I1, I1) - HESS(I1, J1)
               HESS(J1, J1) = HESS(J1, J1) - HESS(I1, J1)
               HESS(I1, I2) = HESS(I1, I2) - HESS(I1, J2)
               HESS(J1, J2) = HESS(J1, J2) - HESS(I1, J2)
               HESS(I2, I2) = HESS(I2, I2) - HESS(I2, J2)
               HESS(J2, J2) = HESS(J2, J2) - HESS(I2, J2)

            END IF
         ENDDO
      ENDDO

      DO I = 2, 3*N
         DO J = 1, I - 1
            HESS(I, J) = HESS(J, I)
         ENDDO
      ENDDO
!     WRITE(*,'(2I6,G15.5)') ((I,J,HESS(I,J),I=3,3),J=3,3)
      
      RETURN
      END
