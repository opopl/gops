      SUBROUTINE NEWTIP (X, G, ENERGY, GTEST, STEST)

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE, STCHRG
      USE KEY, ONLY: TIPID  

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, REALNATOMS, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), FRQN(3*NATOMS)
      DOUBLE PRECISION :: ENERGY, R2, R6, R12, RSS2, ABSR, DUMMY
      DOUBLE PRECISION :: RI(3), RSS(3), P(3), R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DVDR(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), D2VDR2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: D2RMI1(3,3), D2RMI2(3,3), D2RMI3(3,3), D2RMI12(3,3), D2RMI23(3,3), D2RMI31(3,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: D2R1(NRBSITES*NATOMS/2,3), D2R2(NRBSITES*NATOMS/2,3), D2R3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: D2R12(NRBSITES*NATOMS/2,3), D2R23(NRBSITES*NATOMS/2,3), D2R31(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: DOTI1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTI2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTI3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2) 
      DOUBLE PRECISION :: DOTJ1(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2), DOTJ2(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2)
      DOUBLE PRECISION :: DOTJ3(NRBSITES*NATOMS/2,NRBSITES*NATOMS/2) 
      DOUBLE PRECISION :: C12, C6, CH2O
      LOGICAL          :: GTEST, STEST

      D2VDR2(:,:) = 0.D0; DVDR(:,:) = 0.D0
      DOTI1(:,:) = 0.D0; DOTI2(:,:) = 0.D0; DOTI3(:,:) = 0.D0
      DOTJ1(:,:) = 0.D0; DOTJ2(:,:) = 0.D0; DOTJ3(:,:) = 0.D0

      CALL DEFTIP4(C12, C6, CH2O)
 
      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0
      IF (STEST) HESS(:,:) = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
  
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         P  = X(J5-2:J5)

         CALL RMDFAS(P, RMI, DRMI1, DRMI2, DRMI3, D2RMI1, D2RMI2, D2RMI3, D2RMI12, D2RMI23, D2RMI31, GTEST, STEST)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,RBSITE(J2,:))

            IF (GTEST .OR. STEST) THEN
 
               DR1(J4,:) = MATMUL(DRMI1,RBSITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,RBSITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,RBSITE(J2,:))

            ENDIF

            IF (STEST) THEN

               D2R1(J4,:) = MATMUL(D2RMI1,RBSITE(J2,:))
               D2R2(J4,:) = MATMUL(D2RMI2,RBSITE(J2,:))
               D2R3(J4,:) = MATMUL(D2RMI3,RBSITE(J2,:))

               D2R12(J4,:) = MATMUL(D2RMI12,RBSITE(J2,:))
               D2R23(J4,:) = MATMUL(D2RMI23,RBSITE(J2,:))
               D2R31(J4,:) = MATMUL(D2RMI31,RBSITE(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, REALNATOMS 

         J3 = 3*J1
         J5 = OFFSET + J3

         DO J2 = J1 + 1, REALNATOMS

            J4 = 3*J2
            J6 = OFFSET + J4

!     O-O LJ CONTRIBUTION

            J7 = NRBSITES*(J1-1) + 1
            J8 = NRBSITES*(J2-1) + 1
            RSS(:) = R(J7,:) - R(J8,:)
            RSS2   = DOT_PRODUCT(RSS(:),RSS(:))
            ABSR   = DSQRT(RSS2)
            R2     = 1.D0/RSS2
            R6     = R2*R2*R2
            R12    = R6*R6
            ENERGY = ENERGY + C12*R12 - C6*R6

            IF (GTEST .OR. STEST) THEN
!     DVDR = DVDR/R
               DVDR(J7,J8) =-6.D0*(2.D0*C12*R12 - C6*R6)*R2
               DVDR(J8,J7) = DVDR(J7,J8)

               DOTI1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J7,:))
               DOTI2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J7,:))
               DOTI3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J7,:))

               DOTJ1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J8,:))
               DOTJ2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J8,:))
               DOTJ3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J8,:))

               G(J3-2:J3) = G(J3-2:J3) + DVDR(J7,J8)*RSS(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR(J7,J8)*RSS(:)

               G(J5-2) = G(J5-2) + DVDR(J7,J8)*DOTI1(J7,J8)
               G(J5-1) = G(J5-1) + DVDR(J7,J8)*DOTI2(J7,J8)
               G(J5)   = G(J5)   + DVDR(J7,J8)*DOTI3(J7,J8)

               G(J6-2) = G(J6-2) - DVDR(J7,J8)*DOTJ1(J7,J8)
               G(J6-1) = G(J6-1) - DVDR(J7,J8)*DOTJ2(J7,J8)
               G(J6)   = G(J6)   - DVDR(J7,J8)*DOTJ3(J7,J8)

               D2VDR2(J7,J8) = (168.D0*C12*R12 - 48.D0*C6*R6)*R2*R2

            ENDIF

            IF (STEST) THEN

                D2VDR2(J8,J7) = D2VDR2(J7,J8)
                DOTI1(J8,J7)  =-DOTJ1(J7,J8)
                DOTI2(J8,J7)  =-DOTJ2(J7,J8)
                DOTI3(J8,J7)  =-DOTJ3(J7,J8)
                DOTJ1(J8,J7)  =-DOTI1(J7,J8)
                DOTJ2(J8,J7)  =-DOTI2(J7,J8)
                DOTJ3(J8,J7)  =-DOTI3(J7,J8)

            ENDIF

            DO I = 2, NRBSITES 

               J7 = NRBSITES*(J1-1) + I

               DO J = 2, NRBSITES 

                  J8     = NRBSITES*(J2-1) + J
                  RSS(:) = R(J7,:) - R(J8,:)
                  RSS2   = DOT_PRODUCT(RSS(:),RSS(:))
                  R2     = 1.D0/RSS2
                  ABSR   = DSQRT(RSS2)
                  ENERGY = ENERGY + CH2O*STCHRG(I)*STCHRG(J)/ABSR

                  IF (GTEST .OR. STEST) THEN
!     DVDR = DVDR/R
                     DVDR(J7,J8) =-CH2O*STCHRG(I)*STCHRG(J)*R2/ABSR
                     DVDR(J8,J7) = DVDR(J7,J8)

                     DOTI1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J7,:))
                     DOTI2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J7,:))
                     DOTI3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J7,:))

                     DOTJ1(J7,J8) = DOT_PRODUCT(RSS(:),DR1(J8,:))
                     DOTJ2(J7,J8) = DOT_PRODUCT(RSS(:),DR2(J8,:))
                     DOTJ3(J7,J8) = DOT_PRODUCT(RSS(:),DR3(J8,:))

                     G(J3-2:J3)  = G(J3-2:J3) + DVDR(J7,J8)*RSS(:)
                     G(J4-2:J4)  = G(J4-2:J4) - DVDR(J7,J8)*RSS(:)

                     G(J5-2)     = G(J5-2) + DVDR(J7,J8)*DOTI1(J7,J8)
                     G(J5-1)     = G(J5-1) + DVDR(J7,J8)*DOTI2(J7,J8)
                     G(J5)       = G(J5)   + DVDR(J7,J8)*DOTI3(J7,J8)

                     G(J6-2)     = G(J6-2) - DVDR(J7,J8)*DOTJ1(J7,J8)
                     G(J6-1)     = G(J6-1) - DVDR(J7,J8)*DOTJ2(J7,J8)
                     G(J6)       = G(J6)   - DVDR(J7,J8)*DOTJ3(J7,J8)

                     D2VDR2(J7,J8) = 3.D0*CH2O*STCHRG(I)*STCHRG(J)*R2*R2/ABSR
 
                  ENDIF

                  IF (STEST) THEN

                     D2VDR2(J8,J7) = D2VDR2(J7,J8)
                     DOTI1(J8,J7)  =-DOTJ1(J7,J8)
                     DOTI2(J8,J7)  =-DOTJ2(J7,J8)
                     DOTI3(J8,J7)  =-DOTJ3(J7,J8)
                     DOTJ1(J8,J7)  =-DOTI1(J7,J8)
                     DOTJ2(J8,J7)  =-DOTI2(J7,J8)
                     DOTJ3(J8,J7)  =-DOTI3(J7,J8)

                  ENDIF

               ENDDO

            ENDDO

         ENDDO

      ENDDO

      IF (STEST) THEN

         DO J1 = 1, REALNATOMS

            J3 = 3*J1
            J5 = OFFSET + J3

            DO J2 = 1, REALNATOMS

               IF (J1 == J2) CYCLE

               J4 = 3*J2
               J6 = OFFSET + J4

               DO I = 1, NRBSITES

                  J7 = NRBSITES*(J1 - 1) + I

                  DO J = 1, NRBSITES 

                     J8 = NRBSITES*(J2 - 1) + J

                     RSS(:) = R(J7,:) - R(J8,:) 
                 
!     [1] SIX COMPLETELY DIAGONAL TERMS: SAME MOLECULE, SAME COORDINATES

!     xi,xi
                     HESS(J3-2,J3-2) = HESS(J3-2,J3-2) + D2VDR2(J7,J8)*RSS(1)*RSS(1) + DVDR(J7,J8)
!     yi,yi             
                     HESS(J3-1,J3-1) = HESS(J3-1,J3-1) + D2VDR2(J7,J8)*RSS(2)*RSS(2) + DVDR(J7,J8)
!     zi,zi
                     HESS(J3,J3)     = HESS(J3,J3)     + D2VDR2(J7,J8)*RSS(3)*RSS(3) + DVDR(J7,J8)
!     pi1,pi1
                     HESS(J5-2,J5-2) = HESS(J5-2,J5-2) + D2VDR2(J7,J8)*DOTI1(J7,J8)*DOTI1(J7,J8) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(DR1(J7,:),DR1(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R1(J7,:))
!     pi2,pi2
                     HESS(J5-1,J5-1) = HESS(J5-1,J5-1) + D2VDR2(J7,J8)*DOTI2(J7,J8)*DOTI2(J7,J8) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(DR2(J7,:),DR2(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R2(J7,:))
!     pi3,pi3
                     HESS(J5,J5)     = HESS(J5,J5) + D2VDR2(J7,J8)*DOTI3(J7,J8)*DOTI3(J7,J8) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(DR3(J7,:),DR3(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R3(J7,:))

!     [2] OFF-DIAGONAL TERMS ON THE DIAGONAL BLOCKS: SAME MOLECULE, DIFFERENT COORDINATES

!     xi,yi
                     DUMMY           = D2VDR2(J7,J8)*RSS(1)*RSS(2)
                     HESS(J3-2,J3-1) = HESS(J3-2,J3-1) + DUMMY
                     HESS(J3-1,J3-2) = HESS(J3-1,J3-2) + DUMMY
!     yi,zi
                     DUMMY           = D2VDR2(J7,J8)*RSS(2)*RSS(3)
                     HESS(J3-1,J3)   = HESS(J3-1,J3) + DUMMY
                     HESS(J3,J3-1)   = HESS(J3,J3-1) + DUMMY
!     zi,xi
                     DUMMY           = D2VDR2(J7,J8)*RSS(3)*RSS(1)
                     HESS(J3,J3-2)   = HESS(J3,J3-2) + DUMMY
                     HESS(J3-2,J3)   = HESS(J3-2,J3) + DUMMY
!     xi,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*RSS(1) + DVDR(J7,J8)*DR1(J7,1)
                     HESS(J3-2,J5-2) = HESS(J3-2,J5-2) + DUMMY
                     HESS(J5-2,J3-2) = HESS(J5-2,J3-2) + DUMMY
!     yi,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*RSS(2) + DVDR(J7,J8)*DR1(J7,2)
                     HESS(J3-1,J5-2) = HESS(J3-1,J5-2) + DUMMY
                     HESS(J5-2,J3-1) = HESS(J5-2,J3-1) + DUMMY
!     zi,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*RSS(3) + DVDR(J7,J8)*DR1(J7,3)
                     HESS(J3,J5-2)   = HESS(J3,J5-2) + DUMMY
                     HESS(J5-2,J3)   = HESS(J5-2,J3) + DUMMY
!     xi,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*RSS(1) + DVDR(J7,J8)*DR2(J7,1)
                     HESS(J3-2,J5-1) = HESS(J3-2,J5-1) + DUMMY
                     HESS(J5-1,J3-2) = HESS(J5-1,J3-2) + DUMMY
!     yi,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*RSS(2) + DVDR(J7,J8)*DR2(J7,2)
                     HESS(J3-1,J5-1) = HESS(J3-1,J5-1) + DUMMY
                     HESS(J5-1,J3-1) = HESS(J5-1,J3-1) + DUMMY
!     zi,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*RSS(3) + DVDR(J7,J8)*DR2(J7,3)
                     HESS(J3,J5-1)   = HESS(J3,J5-1) + DUMMY
                     HESS(J5-1,J3)   = HESS(J5-1,J3) + DUMMY
!     xi,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*RSS(1) + DVDR(J7,J8)*DR3(J7,1)
                     HESS(J3-2,J5)   = HESS(J3-2,J5) + DUMMY
                     HESS(J5,J3-2)   = HESS(J5,J3-2) + DUMMY
!     yi,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*RSS(2) + DVDR(J7,J8)*DR3(J7,2)
                     HESS(J3-1,J5)   = HESS(J3-1,J5) + DUMMY
                     HESS(J5,J3-1)   = HESS(J5,J3-1) + DUMMY
!     zi,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*RSS(3) + DVDR(J7,J8)*DR3(J7,3)
                     HESS(J3,J5)     = HESS(J3,J5) + DUMMY
                     HESS(J5,J3)     = HESS(J5,J3) + DUMMY
!     pi1,pi2
                     DUMMY           = D2VDR2(J7,J8)*DOTI1(J7,J8)*DOTI2(J7,J8) + DVDR(J7,J8)*DOT_PRODUCT(DR2(J7,:),DR1(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R12(J7,:))
                     HESS(J5-2,J5-1) = HESS(J5-2,J5-1) + DUMMY
                     HESS(J5-1,J5-2) = HESS(J5-1,J5-2) + DUMMY
!     pi2,pi3
                     DUMMY           = D2VDR2(J7,J8)*DOTI2(J7,J8)*DOTI3(J7,J8) + DVDR(J7,J8)*DOT_PRODUCT(DR3(J7,:),DR2(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R23(J7,:))
                     HESS(J5-1,J5)   = HESS(J5-1,J5) + DUMMY
                     HESS(J5,J5-1)   = HESS(J5,J5-1) + DUMMY
!     pi3,pi1
                     DUMMY           = D2VDR2(J7,J8)*DOTI3(J7,J8)*DOTI1(J7,J8) + DVDR(J7,J8)*DOT_PRODUCT(DR1(J7,:),DR3(J7,:)) &
                                     + DVDR(J7,J8)*DOT_PRODUCT(RSS,D2R31(J7,:))
                     HESS(J5,J5-2)   = HESS(J5,J5-2) + DUMMY
                     HESS(J5-2,J5)   = HESS(J5-2,J5) + DUMMY

                 
!     [3] DIAGONAL ELEMENTS ON OFF-DIAGONAL BLOCKS: DIFFERENT MOLECULES, SAME COORDINATE

!     xi,xj
                     HESS(J3-2,J4-2) = HESS(J3-2,J4-2) - D2VDR2(J7,J8)*RSS(1)*RSS(1) - DVDR(J7,J8)
!     yi,yj
                     HESS(J3-1,J4-1) = HESS(J3-1,J4-1) - D2VDR2(J7,J8)*RSS(2)*RSS(2) - DVDR(J7,J8)
!     zi,zj
                     HESS(J3,J4)     = HESS(J3,J4)     - D2VDR2(J7,J8)*RSS(3)*RSS(3) - DVDR(J7,J8)
!     pi1,pj1
                     HESS(J5-2,J6-2) = HESS(J5-2,J6-2) - D2VDR2(J7,J8)*DOTJ1(J7,J8)*DOTI1(J7,J8) &
                                     - DVDR(J7,J8)*DOT_PRODUCT(DR1(J8,:),DR1(J7,:))
!     pi2,pj2
                     HESS(J5-1,J6-1) = HESS(J5-1,J6-1) - D2VDR2(J7,J8)*DOTJ2(J7,J8)*DOTI2(J7,J8) &
                                     - DVDR(J7,J8)*DOT_PRODUCT(DR2(J8,:),DR2(J7,:))
!     pi3,pj3
                     HESS(J5,J6)     = HESS(J5,J6)     - D2VDR2(J7,J8)*DOTJ3(J7,J8)*DOTI3(J7,J8) &
                                    - DVDR(J7,J8)*DOT_PRODUCT(DR3(J8,:),DR3(J7,:))

!     [4] COMPLETELY OFF-DIAGONAL TERMS: DIFFERENT MOLECULES, DIFFERENT COORDINATES

!     xi,yj
                     DUMMY           = - D2VDR2(J7,J8)*RSS(1)*RSS(2)
                     HESS(J3-2,J4-1) = HESS(J3-2,J4-1) + DUMMY
                     HESS(J4-1,J3-2) = HESS(J4-1,J3-2) + DUMMY
!     yi,zj
                     DUMMY           = - D2VDR2(J7,J8)*RSS(2)*RSS(3)
                     HESS(J3-1,J4)   = HESS(J3-1,J4) + DUMMY
                     HESS(J4,J3-1)   = HESS(J4,J3-1) + DUMMY
!     zi,xj
                     DUMMY           = - D2VDR2(J7,J8)*RSS(3)*RSS(1)
                     HESS(J3,J4-2)   = HESS(J3,J4-2) + DUMMY
                     HESS(J4-2,J3)   = HESS(J4-2,J3) + DUMMY


!     xi,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ1(J7,J8)*RSS(1) - DVDR(J7,J8)*DR1(J8,1)
                     HESS(J3-2,J6-2) = HESS(J3-2,J6-2) + DUMMY
                     HESS(J6-2,J3-2) = HESS(J6-2,J3-2) + DUMMY
!     yi,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ1(J7,J8)*RSS(2) - DVDR(J7,J8)*DR1(J8,2)
                     HESS(J3-1,J6-2) = HESS(J3-1,J6-2) + DUMMY
                     HESS(J6-2,J3-1) = HESS(J6-2,J3-1) + DUMMY
!     zi,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ1(J7,J8)*RSS(3) - DVDR(J7,J8)*DR1(J8,3)
                     HESS(J3,J6-2)   = HESS(J3,J6-2) + DUMMY
                     HESS(J6-2,J3)   = HESS(J6-2,J3) + DUMMY
!     xi,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ2(J7,J8)*RSS(1) - DVDR(J7,J8)*DR2(J8,1)
                     HESS(J3-2,J6-1) = HESS(J3-2,J6-1) + DUMMY
                     HESS(J6-1,J3-2) = HESS(J6-1,J3-2) + DUMMY
!     yi,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ2(J7,J8)*RSS(2) - DVDR(J7,J8)*DR2(J8,2)
                     HESS(J3-1,J6-1) = HESS(J3-1,J6-1) + DUMMY
                     HESS(J6-1,J3-1) = HESS(J6-1,J3-1) + DUMMY
!     zi,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ2(J7,J8)*RSS(3) - DVDR(J7,J8)*DR2(J8,3)
                     HESS(J3,J6-1)   = HESS(J3,J6-1) + DUMMY
                     HESS(J6-1,J3)   = HESS(J6-1,J3) + DUMMY
!     xi,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ3(J7,J8)*RSS(1) - DVDR(J7,J8)*DR3(J8,1)
                     HESS(J3-2,J6)   = HESS(J3-2,J6) + DUMMY
                     HESS(J6,J3-2)   = HESS(J6,J3-2) + DUMMY
!     yi,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ3(J7,J8)*RSS(2) - DVDR(J7,J8)*DR3(J8,2)
                     HESS(J3-1,J6)   = HESS(J3-1,J6) + DUMMY
                     HESS(J6,J3-1)   = HESS(J6,J3-1) + DUMMY
!     zi,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTJ3(J7,J8)*RSS(3) - DVDR(J7,J8)*DR3(J8,3)
                     HESS(J3,J6)     = HESS(J3,J6) + DUMMY
                     HESS(J6,J3)     = HESS(J6,J3) + DUMMY
!     pi1,pj2
                     DUMMY           = - D2VDR2(J7,J8)*DOTI1(J7,J8)*DOTJ2(J7,J8) - DVDR(J7,J8)*DOT_PRODUCT(DR2(J8,:),DR1(J7,:))
                     HESS(J5-2,J6-1) = HESS(J5-2,J6-1) + DUMMY
                     HESS(J6-1,J5-2) = HESS(J6-1,J5-2) + DUMMY
!     pi2,pj3
                     DUMMY           = - D2VDR2(J7,J8)*DOTI2(J7,J8)*DOTJ3(J7,J8) - DVDR(J7,J8)*DOT_PRODUCT(DR3(J8,:),DR2(J7,:))
                     HESS(J5-1,J6)   = HESS(J5-1,J6) + DUMMY
                     HESS(J6,J5-1)   = HESS(J6,J5-1) + DUMMY
!     pi3,pj1
                     DUMMY           = - D2VDR2(J7,J8)*DOTI3(J7,J8)*DOTJ1(J7,J8) - DVDR(J7,J8)*DOT_PRODUCT(DR1(J8,:),DR3(J7,:))
                     HESS(J5,J6-2)   = HESS(J5,J6-2) + DUMMY
                     HESS(J6-2,J5)   = HESS(J6-2,J5) + DUMMY

                  ENDDO

               ENDDO

            ENDDO

         ENDDO

!      CALL NRMLMD(X,FRQN)

!      STOP

      ENDIF

      END SUBROUTINE NEWTIP 

!     ---------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP4(C12, C6, CH2O)
!     TIP4P water

      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE, STCHRG

      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETA, ROH, ROM, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      ROM   = 0.15D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      RBSITE(1,1) = 0.D0
      RBSITE(1,2) = 0.D0
      RBSITE(1,3) = 0.D0

      RBSITE(2,1) = 0.D0
      RBSITE(2,2) = SIN(0.5D0*THETA)*ROH
      RBSITE(2,3) = COS(0.5D0*THETA)*ROH

      RBSITE(3,1) = 0.D0
      RBSITE(3,2) = -SIN(0.5D0*THETA)*ROH
      RBSITE(3,3) = COS(0.5D0*THETA)*ROH

      RBSITE(4,1) = 0.D0
      RBSITE(4,2) = 0.D0
      RBSITE(4,3) = ROM

      STCHRG(:) = (/0.D0, 0.52D0, 0.52D0, -1.04D0/)
      C6        = 2552.24D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
      C12       = 2510.4D3
      CH2O      = 1389.354848D0 ! Conversion factor for coulomb energy

      M(:)  = (/16.D0, 1.D0, 1.D0, 0.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*RBSITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         RBSITE(I,:) = RBSITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP4

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE NRMLMD (X, FRQN)

      USE COMMONS
      USE MODHESS

      IMPLICIT NONE

      INTEGER          :: I, J, K, J1, J2, J3, J5, OFFSET, NDIM, IR, IC, K1, K2 
      DOUBLE PRECISION :: X(3*NATOMS), FRQN(3*NATOMS)
      DOUBLE PRECISION :: MS(3), TMASS, FRQCNV
      DOUBLE PRECISION :: KBLOCK(3,3), KBEGNV(3)
      DOUBLE PRECISION :: P(3), RMI(3,3), DRMI(3,3), DR(3)
      DOUBLE PRECISION :: KD(3*NATOMS), U(3*NATOMS,3*NATOMS)
      DOUBLE PRECISION :: HUK(3*NATOMS,3*NATOMS), AP(3*NATOMS,3*NATOMS)
      LOGICAL          :: GTEST, STEST
!     the following required to call the LAPACK routine DSYEV
      INTEGER          :: INFO
      INTEGER, PARAMETER :: LWORK = 10000 ! the dimension is set arbitrarily
      DOUBLE PRECISION :: WORK(LWORK)

!     Computes the normal modes and frequencies
!     following Pohorille et al. JCP 87, 6070 (1987)
!     INPUT:
!     ndim   : number of degrees of freedom
!     second : second derivatives (ORIENT convention) from derivs
!     common "sites" : coordinates of sites in global coordinates
!     OUTPUT:
!     freq   : eigenfrequencies of normal modes

!     We adopt Pohorille's notation for clarity:
!     K matrix : block diagonal kinetic energy matrix (6N x 6N)
!     kblock : for each rigid body we have one 6x6 matrix which consists
!     of two blocks : left upper or "mass" diagonal submatrix
!     and right lower or inertia tensor matrix. Here kblock
!     is the 3 x 3 inertia tensor.
!     KD     : using the diagonalized kinetic energy tensor
!     we keep track of the diagonal of KD only
!     U      : eigenvector matrix (Pohorille's S)

!     Frequency conversion factor: Energy in KJ/mol and length in angstrom
!     to get frequencies in cm^{-1}
      
      FRQCNV = 1.D03/(2.D0*4.D0*DATAN(1.D0)*2.998D0)
      MS(:)  = (/16.D0, 1.D0, 1.D0/)
      TMASS  = 18.D0

!     Initialize

      U(:,:) = 0.D0
      IR     = 0
      IC     = 0


!     Get the site positions

      OFFSET = 3*NATOMS/2
      GTEST = .FALSE.; STEST = .FALSE.

      DO J1 = 1, NATOMS/2

         J3 = 3*J1
         J5 = OFFSET + J3
         P  = X(J5-2:J5)
         KBLOCK(:,:) = 0.D0

         CALL RMDFAS(P, RMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, GTEST, STEST)

         DO J2 = 1, NRBSITES - 1

            DR(:)  = MATMUL(RMI(:,:),RBSITE(J2,:))

            DO I = 1, 3

               KBLOCK(I,I) = KBLOCK(I,I) + MS(J2)*(DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))

               DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix

                  KBLOCK(I,J) = KBLOCK(I,J) - MS(J2)*DR(I)*DR(J)

               ENDDO

            ENDDO

         ENDDO

!     Diagonalise KBLOCK using LAPACK routine DSYEV
!     DSYEV computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix KBLOCK 

         CALL DSYEV('V','L',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)

!     The character 'V' instructs to return the eigenvector as well
!     On exit, if INFO = 0, KBLOCK contains the orthonormal eigenvectors of the matrix KBLOCK in columns
!     The character 'L' tells that the Lower triangle of KBLOCK is stored 
!     Next is the order of the matrix KBLOCK
!     The integer after KBLOCK is the leading dimension of the array KBLOCK
!     KEGNV holds the eigenvalues in ascending order if INFO = 0

         IF (INFO /= 0) THEN
            WRITE(*,*) 'NRMLMD > Error in DSYEV with KBLOCK, INFO =', INFO
            STOP
         ENDIF

!     Construction of the matrix U
!     First: translation coordinates

!         U(IR+1,IC+1) = 1.D0; U(IR+2,IC+2) = 1.D0; U(IR+3,IC+3) = 1.D0
         U(IR+1,IC+1) = 1.D0/SQRT(3.D0); U(IR+2,IC+2) = 1.D0/SQRT(3.D0); U(IR+3,IC+3) = 1.D0/SQRT(3.D0)
         KD(IC+1:IC+3) = 1.D0/SQRT(TMASS)
            
!     Now rotational coordinates

         U(OFFSET+IR+1:OFFSET+IR+3,OFFSET+IC+1:OFFSET+IC+3) = KBLOCK(:,:) 
         KD(OFFSET+IC+1:OFFSET+IC+3) = 1.D0/SQRT(KBEGNV(:))

         IR = IR + 3
         IC = IC + 3 

      ENDDO

      NDIM = 3*NATOMS

!      DO I = 1, NDIM
!         DO J = 1, NDIM
!            HUK(I,J) = 0.D0
!            DO K = 1, NDIM
!               HUK(I,J) = HUK(I,J) + HESS(I,K)*U(K,J)*KD(J)
!            ENDDO
!         ENDDO
!      ENDDO
!      DO I = 1, NDIM
!         DO J = 1, NDIM
!            AP(I,J) = 0.D0
!            DO K  = 1, NDIM
!               AP(I,J) = AP(I,J) + KD(I)*U(K,I)*HUK(K,J) ! times U transpose
!            ENDDO
!         ENDDO 
!      ENDDO

      AP(:,:) = 0.D0
      DO I = 1, NDIM
         DO J = 1, I
            DO K1 = 1, NDIM
               DO K2 = 1, NDIM
                  AP(I,J) = AP(I,J) + U(K1,I)*HESS(K1,K2)*U(K2,J)
               ENDDO
            ENDDO
            AP(I,J) = KD(I)*AP(I,J)*KD(J)
         ENDDO
      ENDDO

      CALL DSYEV('V','L',NDIM,AP,NDIM,FRQN,WORK,LWORK,INFO)

      FRQN(:) = FRQCNV*SQRT(ABS(FRQN(:)))

      PRINT *, 'FRQCNV=', FRQCNV
 
      PRINT *, FRQN(:)

      END SUBROUTINE NRMLMD 
