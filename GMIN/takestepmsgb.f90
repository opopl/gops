      SUBROUTINE TAKESTEPMSGB (NP)

!     THIS ROUTINE TAKES STEP FOR MULTI-SITE RIGID BODIES ENSURING NO OVERLAP

      USE COMMONS  

      IMPLICIT NONE

      INTEGER          :: NP, JMAX, JMAX2, REALNATOMS, OFFSET, PTINDX
      INTEGER          :: J1, J2, J3, J4, J5, J6, I, J
      LOGICAL          :: OVRLPT

      DOUBLE PRECISION :: PI, DUMMY, DUMMY2
      DOUBLE PRECISION :: DIST(3*NATOMS/2), XMASS, YMASS, ZMASS, DMAX, VMAX, VMAX2
      DOUBLE PRECISION :: VMIN, CMMAX, CMDIST(NATOMS/2), LOCALSTEP
      DOUBLE PRECISION :: DPRAND, RANDOM, THETA, PHI
      DOUBLE PRECISION :: PST(NGBSITE,3), OST(NGBSITE,3), ROTMAT(3,3) 
      DOUBLE PRECISION :: RI(3), RJ(3), RISITE(3), RJSITE(3), RIJ(3)
      DOUBLE PRECISION :: P(3), EI(3), EJ(3), ABSRIJ, RCUT, ECFVAL


      PI         = 4.D0*ATAN(1.0D0)
      RCUT       = GBKAPPA + 1.D0
      REALNATOMS = NATOMS/2
      OFFSET     = 3 * REALNATOMS

      CALL DEFMOL(PST,OST)

      DO J1 = 1,REALNATOMS

         J2     = 3*J1
         DUMMY2 = COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2
         IF (DUMMY2 .GT. RADIUS) THEN
            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,X,Y,Z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), &
                                       COORDS(J2-1,NP), COORDS(J2,NP)
            PRINT*, 'INITIAL COORDINATE OUTSIDE CONTAINER -- INCREASE CONTAINER RADIUS'
            STOP
         END IF

      END DO

      DO J1 = 1,3*NATOMS
         COORDSO(J1,NP) = COORDS(J1,NP)
      END DO

      DO J1 = 1,NATOMS/2
         VATO(J1,NP) = VAT(J1,NP)
      END DO

!     FIND THE CENTRE OF MASS

      XMASS = 0.0D0
      YMASS = 0.0D0
      ZMASS = 0.0D0

      DO J1 = 1,NATOMS/2

         XMASS = XMASS + COORDS(3*(J1-1)+1,NP)
         YMASS = YMASS + COORDS(3*(J1-1)+2,NP)
         ZMASS = ZMASS + COORDS(3*(J1-1)+3,NP)

      ENDDO

      XMASS = XMASS/(REALNATOMS)
      YMASS = YMASS/(REALNATOMS)
      ZMASS = ZMASS/(REALNATOMS)

!     FIND THE MOST WEAKLY BOUND ATOM, JMAX, THE SECOND MOST WEAKLY BOUND ATOM, JMAX2,
!     AND THE PAIR ENERGY OF THE MOST TIGHTLY BOUND ATOM, VMIN. AN ANGULAR STEP IS
!     TAKEN FOR JMAX IF ITS PAIR ENERGY IS > ASTEP*VMIN PUTTING THE ATOM AT A RADIUS OF
!     DMAX (OR CMMAX FROM CM OF THE CLUSTER).

      DMAX  =  1.0D0
      VMAX  = -1.0D3
      VMAX2 = -1.0D3
      VMIN  =  0.0D0
      CMMAX =  1.0D0

      DO J1 = 1, REALNATOMS

         J2 = 3*J1
         DIST(J1)   = DSQRT( COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2)
         CMDIST(J1) = DSQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF (CMDIST(J1) .GT. CMMAX) CMMAX = CMDIST(J1)
         IF (DIST(J1) .GT. DMAX) DMAX = DIST(J1)
         IF (VAT(J1,NP) .GT. VMAX) THEN
            VMAX = VAT(J1,NP)
            JMAX = J1
         ELSE IF ((VAT(J1,NP).LT. VMAX) .AND. (VAT(J1,NP) .GT. VMAX2)) THEN
              VMAX2 = VAT(J1,NP)
              JMAX2 = J1
         ENDIF
         IF (VAT(J1,NP) .LT. VMIN) VMIN = VAT(J1,NP)

      ENDDO

      IF (VAT(JMAX,NP) > (ASTEP(NP)*VMIN) .AND. (.NOT.NORESET)) THEN

         J2 = 3*JMAX
         THETA           = DPRAND()*PI
         PHI             = DPRAND()*PI*2.0D0
         COORDS(J2-2,NP) = XMASS + (CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
         COORDS(J2-1,NP) = YMASS + (CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
         COORDS(J2,NP)   = ZMASS + (CMMAX+1.0D0)*DCOS(THETA)
         DUMMY           = COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2

         IF (DUMMY > RADIUS) THEN

            DUMMY           = DSQRT(RADIUS*0.99D0/DUMMY)
            COORDS(J2-2,NP) = COORDS(J2-2,NP)*DUMMY
            COORDS(J2-1,NP) = COORDS(J2-1,NP)*DUMMY
            COORDS(J2,NP)   = COORDS(J2,NP)*DUMMY

         END IF

      ENDIF

      DO J1 = 1, NATOMS

         J3 = 3*J1

         LOCALSTEP = STEP(NP)

         IF (J1 > REALNATOMS) THEN

            LOCALSTEP = 0.0D0
            IF (OMOVE(NP)) LOCALSTEP = OSTEP(NP)

         ELSE IF (J1 <= REALNATOMS) THEN

            LOCALSTEP=0.0D0
            IF (TMOVE(NP)) LOCALSTEP = STEP(NP)

         END IF

!     CHECK FOR OVERLAP

         OVRLPT = .TRUE.

         DO WHILE (OVRLPT)
        
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-2,NP) = COORDS(J3-2,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-1,NP) = COORDS(J3-1,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3,NP)   = COORDS(J3,NP) + LOCALSTEP*RANDOM
          
            IF (J1 <= REALNATOMS) THEN

               PTINDX = J1
               J5    = OFFSET + J3
               RI(:) = COORDS(J3-2:J3,NP)
               P(:)  = COORDS(J5-2:J5,NP)

            ELSE 

               PTINDX = J1 - REALNATOMS
               J5    = J3 - OFFSET
               RI(:) = COORDS(J5-2:J5,NP)
               P(:)  = COORDS(J3-2:J3,NP)

            ENDIF  

!     ROTATION MATRIX

            CALL ROTATION (P, ROTMAT)

            DO I = 1, NGBSITE

!     OBTAIN THE SITE POSITION IN THE SPACE-FIXED FRAME  

                  RISITE = RI +  MATMUL(ROTMAT,PST(I,:))
                  
!     OBTAIN THE SITE ORIENTATION IN THE SPACE-FIXED FRAME

                  EI = MATMUL(ROTMAT,OST(I,:))

                  DO J2 = 1, REALNATOMS

                     IF (J2 == PTINDX) CYCLE

                     J4    = 3*J2
                     J6    = OFFSET + J4
                     RJ(:) = COORDS(J4-2:J4,NP)
                     P(:)  = COORDS(J6-2:J6,NP)

!     ROTATION MATRIX

                     CALL ROTATION (P, ROTMAT)

                     DO J = 1, NGBSITE

!     OBTAIN THE SITE POSITION IN THE SPACE-FIXED FRAME   

                        RJSITE = RJ + MATMUL(ROTMAT,PST(J,:))

!     OBTAIN THE SITE ORIENTATION IN THE SPACE-FIXED FRAME
                   
                        EJ = MATMUL(ROTMAT,OST(J,:))

                        RIJ    = RISITE - RJSITE
                        ABSRIJ = DSQRT(DOT_PRODUCT(RIJ,RIJ))
                       
                        IF (ABSRIJ < RCUT) THEN 

!     DETERMINE ELLIPTIC CONTACT FUNCTION

                           CALL ECFDC (EI, EJ, RIJ, ECFVAL)
 
                           IF (ECFVAL >= 1.D0) THEN
                              OVRLPT = .FALSE.
                           ENDIF

                        ELSE

                           OVRLPT = .FALSE.

                        ENDIF

                  ENDDO  ! END LOOP OVER J

               ENDDO  ! END LOOP OVER J2

            ENDDO  ! END LOOP OVER I

         ENDDO  ! END WHILE
 
      ENDDO  ! END LOOP OVER J1   

      END SUBROUTINE TAKESTEPMSGB

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ROTATION (PI, ROTMAT)

      IMPLICIT NONE

      INTEGER          :: K1, K2
      DOUBLE PRECISION :: E(3,3), I3(3,3), ROT2(3,3), ROTMAT(3,3), PI(3), THETA, THETA2 

      I3(:,:)  = 0.D0
      I3(1,1)  = 1.D0
      I3(2,2)  = 1.D0
      I3(3,3)  = 1.D0
       
      THETA   = DSQRT(DOT_PRODUCT(PI,PI))
      THETA2  = THETA * THETA
      E(:,:)  = 0.D0
      E(1,2)  = -PI(3)
      E(1,3)  =  PI(2)
      E(2,3)  = -PI(1)
      E(2,1)  = -E(1,2)
      E(3,1)  = -E(1,3)
      E(3,2)  = -E(2,3)
      E       = E/THETA

      DO K1 = 1, 3
         DO K2 = 1, 3
            ROT2(K1,K2) = PI(K1)*PI(K2)/THETA2
         ENDDO
      ENDDO

      ROTMAT = I3*COS(THETA) + (1.D0-COS(THETA))*ROT2 - E*SIN(THETA)

      END SUBROUTINE ROTATION

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ECFDC (EI, EJ, RIJ, ECFVAL)

      USE COMMONS, ONLY : LPRSQ, LSQDFR

      IMPLICIT NONE

      INTEGER          :: K1, K2
      DOUBLE PRECISION :: RIJ(3), EI(3), EJ(3)
      DOUBLE PRECISION :: AE(3,3), BE(3,3), FMIN, ECFVAL

      DO K1 = 1, 3

         DO K2 = K1, 3

            IF (K1 == K2) THEN

               AE(K1,K2) = LPRSQ + LSQDFR * EI(K1) * EI(K2)
               BE(K1,K2) = LPRSQ + LSQDFR * EJ(K1) * EJ(K2)

            ELSE

               AE(K1,K2) = LSQDFR * EI(K1)*EI(K2)
               BE(K1,K2) = LSQDFR * EJ(K1)*EJ(K2)

            ENDIF

         ENDDO

      ENDDO

      CALL BRENTMINGB(AE, BE, RIJ, 0.D0, 0.51D0, 1.D0, FMIN)

      ECFVAL = - FMIN

      RETURN
      END SUBROUTINE ECFDC

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE OBJCTFGB(AE, BE, RIJ, LAMDA, SLMD)

      IMPLICIT NONE

      INTEGER          :: I, J
      DOUBLE PRECISION :: AE(3,3), BE(3,3), RIJ(3)
      DOUBLE PRECISION :: G(3,3), GINVR(3), QFGINV
      DOUBLE PRECISION :: LAMDA, SLMD

      DO I = 1, 3

         DO J = I, 3

            G(I,J) = (1.D0 - LAMDA) * AE(I,J) + LAMDA * BE(I,J)

         ENDDO

      ENDDO

!     DETERMINE THE INVERSE MATRIX

      CALL MTRXINGB(G)

      GINVR =  MATMUL(G, RIJ)

      SLMD  =  - LAMDA * (1.D0 - LAMDA) * DOT_PRODUCT(RIJ,GINVR)

      RETURN
      END SUBROUTINE OBJCTFGB

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE MTRXINGB(S)

      IMPLICIT NONE

      INTEGER            :: I, J
      INTEGER, PARAMETER :: M = 3
      DOUBLE PRECISION   :: S(M,M), A(M,M), DET, INVDET

!     ADJOINT

      A(1,1) = S(2,2) * S(3,3) - S(2,3) * S(2,3)
      A(2,2) = S(1,1) * S(3,3) - S(1,3) * S(1,3)
      A(3,3) = S(1,1) * S(2,2) - S(1,2) * S(1,2)
      A(1,2) = S(1,3) * S(2,3) - S(1,2) * S(3,3)
      A(1,3) = S(1,2) * S(2,3) - S(1,3) * S(2,2)
      A(2,3) = S(1,2) * S(1,3) - S(1,1) * S(2,3)

!     DETERMINANT

      DET    =  S(1,1)* A(1,1) + S(1,2) * A(1,2) + S(1,3) * A(1,3)
      INVDET = 1.D0 / DET

      DO I = 1, 3

         DO J = 1, 3

            IF (I .GT. J) A(I,J) = A(J,I)
            S(I,J) = A(I,J) * INVDET

         ENDDO

      ENDDO

      END SUBROUTINE MTRXINGB

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE BRENTMINGB (AE, BE, RIJ, AX, BX, CX, FMIN)

      IMPLICIT NONE

      
      INTEGER            :: ITR
      INTEGER, PARAMETER :: ITRMX = 200
      DOUBLE PRECISION   :: MA(3,3), MB(3,3)
      DOUBLE PRECISION   :: AX, BX, CX, A, B, D, E, P, Q, R, U, V, W, X, XM
      DOUBLE PRECISION   :: XMIN, FX, FU, FV, FW, FMIN, F, ETMP
      DOUBLE PRECISION   :: AE(3,3), BE(3,3), RIJ(3)
      DOUBLE PRECISION   :: TOL1, TOL2, CGOLD
      DOUBLE PRECISION, PARAMETER :: TOL = 1.D-08, ZEPS = 1.D-10

      CGOLD = 0.5D0 * (3 - DSQRT(5.D0))
      A     = MIN(AX,CX)
      B     = MAX(AX,CX)
      V     = BX
      W     = V
      X     = V
      E     = 0.D0
!     FX    = F(X)

      CALL OBJCTFGB (AE, BE, RIJ, X, FX)

      FV    = FX
      FW    = FX

      DO 10 ITR = 1, ITRMX

         XM   = 0.5D0 * (A + B)
         TOL1 = TOL * ABS(X) + ZEPS
         TOL2 = 2.D0 * TOL1
         IF (ABS(X - XM) .LE. (TOL2 - 0.5D0 * (B - A))) GOTO 3
         IF (ABS(E) .GT. TOL1) THEN
            R    = (X - W) * (FX - FV)
            Q    = (X - V) * (FX - FW)
            P    = (X - V) * Q - (X - W) * R
            Q    = 2.D0 * (Q - R)
            IF (Q .GT. 0.D0) P = -P
            Q    = ABS(Q)
            ETMP = E
            E    = D

            IF (ABS(P) .GE. ABS(0.5D0*Q*ETMP) .OR. P .LE. Q*(A-X) .OR. P .GE. Q*(B-X)) GOTO 1
!     %         .OR. P .GE. Q*(B-X)) GOTO 1
!     THE ABOVE CONDITIONS DETERMINE THE ACCEPTABILITY OF THE PARABOLIC FIT. HERE IT IS O.K.
            D = P / Q ! TAKE THE PARABOLIC STEP.
            U = X + D
            IF(U-A .LT. TOL2 .OR. B-U .LT. TOL2) D=SIGN(TOL1,XM-X)
            GOTO 2 !SKIP OVER THE GOLDEN SECTION STEP.

         ENDIF

1        IF (X .GE. XM) THEN
!      WE ARRIVE HERE FOR A GOLDEN SECTION STEP, WHICH WE TAKE
            E = A - X !INTO THE LARGER OF THE TWO SEGMENTS.
         ELSE
            E = B - X
         ENDIF

         D = CGOLD * E !TAKE THE GOLDEN SECTION STEP.
2        IF (ABS(D) .GE. TOL1) THEN ! ARRIVE HERE WITH D COMPUTED EITHER FROM PARABOLIC FIT, OR
            U = X + D !ELSE FROM GOLDEN SECTION.
         ELSE
            U = X + SIGN(TOL1,D)
         ENDIF

!         FU = F(U) !THIS IS THE ONE FUNCTION EVALUATION PER ITERATION,

         CALL OBJCTFGB (AE, BE, RIJ, U, FU)

         IF (FU .LE. FX) THEN !AND NOW WE HAVE TO DECIDE WHAT TO DO WITH OUR FUNCTION
            IF (U .GE. X) THEN !EVALUATION. HOUSEKEEPING FOLLOWS:
               A = X
            ELSE
               B = X
            ENDIF
            V  = W
            FV = FW
            W  = X
            FW = FX
            X  = U
            FX = FU
         ELSE
            IF (U .LT. X) THEN
               A = U
            ELSE
               B = U
            ENDIF

            IF (FU .LE. FW .OR. W .EQ. X) THEN
               V  = W
               FV = FW
               W  = U
               FW = FU
            ELSE IF (FU .LE. FV .OR. V .EQ. X .OR. V .EQ. W) THEN
               V = U
               FV = FU
            ENDIF

         ENDIF ! DONE WITH HOUSEKEEPING. BACK FOR ANOTHER ITERATION.

10    CONTINUE
      PAUSE 'BRENT EXCEED MAXIMUM ITERATIONS'
3     XMIN = X ! ARRIVE HERE READY TO EXIT WITH BEST VALUES.
      FMIN = FX
     
      RETURN
      END SUBROUTINE BRENTMINGB
