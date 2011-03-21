      SUBROUTINE NEWCAPSID (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, RHO
 
      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: RISITE(3), RJSITE(3), RI(3), RJ(3), DSS(3), PI(3)
      DOUBLE PRECISION :: DSS2, RSS, DVDR, ENERGY, SIGMASQ, VSS
      DOUBLE PRECISION :: R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: CAPEPS2, CAPRAD, CAPHEIGHT
      LOGICAL          :: GTEST
      COMMON / CAPS /  CAPEPS2, CAPRAD, CAPHEIGHT

      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0
      
      OFFSET  = 3*NATOMS/2
      SIGMASQ = (1.D0 + CAPRAD*SQRT(0.5D0*(5.D0 + SQRT(5.D0))))**2
    
      DO J1 = 1, NATOMS/2

         J3 = 3*J1
         J5 = OFFSET + J3
         RI = X(J3-2:J3)
         PI = X(J5-2:J5)

         CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, GTEST)

         DO J2 = 1, NRBSITES

            J4        = NRBSITES*(J1-1) + J2
            R(J4,:)   = RI(:) + MATMUL(RMI,SITE(J2,:))

            IF (GTEST) THEN

               DR1(J4,:) = MATMUL(DRMI1,SITE(J2,:))
               DR2(J4,:) = MATMUL(DRMI2,SITE(J2,:))
               DR3(J4,:) = MATMUL(DRMI3,SITE(J2,:))

            ENDIF

         ENDDO

      ENDDO

      DO J1 = 1, NATOMS/2 - 1 

         J3 = 3*J1
         J5 = OFFSET + J3

         DO J2 = J1 + 1, NATOMS/2

            J4 = 3*J2
            J6 = OFFSET + J4

!     REPULSIVE SITE CONTRIBUTION

            J7     = NRBSITES*(J1-1) + 6
            J8     = NRBSITES*(J2-1) + 6
            DSS(:) = R(J7,:) - R(J8,:)
            DSS2   = DOT_PRODUCT(DSS,DSS)
              
            VSS    = CAPEPS2*(SIGMASQ/DSS2)**6
            ENERGY = ENERGY + VSS

            IF (GTEST) THEN

               DVDR   = -12.D0*VSS/DSS2
               G(J3-2:J3) = G(J3-2:J3) + DVDR*DSS(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR*DSS(:)

               G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(DSS,DR1(J7,:))
               G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(DSS,DR2(J7,:))
               G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(DSS,DR3(J7,:))

               G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(DSS,DR1(J8,:))
               G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(DSS,DR2(J8,:))
               G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(DSS,DR3(J8,:))

            ENDIF

!     SUM OVER MORSE SITES

            DO I = 1, NRBSITES - 1

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES - 1

                  J8 = NRBSITES*(J2-1) + J
                  DSS(:) = R(J7,:) - R(J8,:)
                  DSS2   = DOT_PRODUCT(DSS(:),DSS(:))
                  RSS    = SQRT(DSS2)
                  VSS    = EXP(RHO*(1.D0 - RSS))
                  ENERGY = ENERGY + (1.D0 - VSS)*(1.D0 - VSS) - 1.D0 

                  IF (GTEST) THEN
 
                     DVDR       = 2.D0*RHO*(-VSS*VSS + VSS)/RSS
                     G(J3-2:J3) = G(J3-2:J3) + DVDR*DSS(:)
                     G(J4-2:J4) = G(J4-2:J4) - DVDR*DSS(:)

                     G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(DSS,DR1(J7,:))
                     G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(DSS,DR2(J7,:))
                     G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(DSS,DR3(J7,:))

                     G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(DSS,DR1(J8,:))
                     G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(DSS,DR2(J8,:))
                     G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(DSS,DR3(J8,:))

                  ENDIF

               ENDDO

            ENDDO

         ENDDO

      ENDDO

      END SUBROUTINE NEWCAPSID 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFCAPSID(CAPRAD,CAPHEIGHT)
    
      USE COMMONS, ONLY: SITE  
      
      IMPLICIT NONE
  
      DOUBLE PRECISION :: CAPRAD, CAPHEIGHT

      SITE(1,1) = 1.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0
      SITE(2,1) = ((-1.D0 + SQRT(5.D0)))/4.D0
      SITE(2,2) = (SQRT((5.D0 + SQRT(5.D0))/2.D0))/2.D0
      SITE(2,3) = 0.D0
      SITE(3,1) = ((-1.D0 - SQRT(5.D0)))/4.0D0
      SITE(3,2) = (SQRT((5.D0 - SQRT(5.D0))/2.D0))/2.D0
      SITE(3,3) = 0.0D0
      SITE(4,1) = ((-1 - SQRT(5.D0)))/4.D0
      SITE(4,2) = -(SQRT((5.D0 - SQRT(5.D0))/2.D0))/2.D0
      SITE(4,3) = 0.D0
      SITE(5,1) = ((-1 + SQRT(5.D0)))/4.D0
      SITE(5,2) = -(SQRT((5.D0 + SQRT(5.D0))/2.D0))/2.D0
      SITE(5,3) = 0.0D0
      SITE(6,1) = 0.0D0
      SITE(6,2) = 0.0D0
      SITE(6,3) = CAPHEIGHT

!      SITE(1:5,:) = CAPRAD*SITE(1:5,:)
      SITE(:,:)   = CAPRAD*SITE(:,:)

      END SUBROUTINE DEFCAPSID
