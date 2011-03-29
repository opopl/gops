      SUBROUTINE NEWTIP(X, G, ENERGY, GTEST) 

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, TIPID

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, OFFSET 
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), CHARGE(NRBSITES)
      DOUBLE PRECISION :: RISITE(3), RJSITE(3), RI(3), RJ(3), DR(3), DSS(3), PI(3), PJ(3)
      DOUBLE PRECISION :: DSS2, R2, R6, R12, ABSR, DVDR, ENERGY, C12, C6, CH2O
      DOUBLE PRECISION :: R(NRBSITES*NATOMS/2,3)
      DOUBLE PRECISION :: DR1(NRBSITES*NATOMS/2,3), DR2(NRBSITES*NATOMS/2,3), DR3(NRBSITES*NATOMS/2,3) 
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
      DOUBLE PRECISION :: DOTI1, DOTI2, DOTI3, DOTJ1, DOTJ2, DOTJ3
      LOGICAL          :: GTEST

      IF (TIPID == 1) THEN
         CALL DEFTIP1(CHARGE, C12, C6, CH2O)
      ELSE IF (TIPID == 2) THEN
         CALL DEFTIP2(CHARGE, C12, C6, CH2O)
      ELSE IF (TIPID == 3) THEN
         CALL DEFTIP3(CHARGE, C12, C6, CH2O)
      ELSE IF (TIPID == 4) THEN
         CALL DEFTIP4(CHARGE, C12, C6, CH2O)
      ELSE IF (TIPID == 5) THEN
         CALL DEFTIP5(CHARGE, C12, C6, CH2O)
      ENDIF

      ENERGY          = 0.D0
      IF (GTEST) G(:) = 0.D0
      OFFSET          = 3*NATOMS/2

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

!     O-O LJ CONTRIBUTION

            J7 = NRBSITES*(J1-1) + 1
            J8 = NRBSITES*(J2-1) + 1   
            DSS(:) = R(J7,:) - R(J8,:)
            DSS2   = DOT_PRODUCT(DSS,DSS)
            R2     = 1.D0 / DSS2
            R6     = R2 * R2 * R2
            R12    = R6 ** 2
            ENERGY = ENERGY + C12*R12 - C6*R6

            IF (GTEST) THEN

!     DVDR = DVDR/R
               DVDR       = -6.D0*(2.D0*C12*R12 - C6*R6)*R2
               G(J3-2:J3) = G(J3-2:J3) + DVDR*DSS(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR*DSS(:)

               G(J5-2) = G(J5-2) + DVDR*DOT_PRODUCT(DSS,DR1(J7,:))
               G(J5-1) = G(J5-1) + DVDR*DOT_PRODUCT(DSS,DR2(J7,:))
               G(J5)   = G(J5)   + DVDR*DOT_PRODUCT(DSS,DR3(J7,:))

               G(J6-2) = G(J6-2) - DVDR*DOT_PRODUCT(DSS,DR1(J8,:))
               G(J6-1) = G(J6-1) - DVDR*DOT_PRODUCT(DSS,DR2(J8,:))
               G(J6)   = G(J6)   - DVDR*DOT_PRODUCT(DSS,DR3(J8,:))

            ENDIF

!     SUM OVER CHARGED SITES

            DO I = 1, NRBSITES

               J7 = NRBSITES*(J1-1) + I

               DO J = 1, NRBSITES

                  J8 = NRBSITES*(J2-1) + J 
                  DSS(:) = R(J7,:) - R(J8,:)
                  DSS2   = DSS(1)*DSS(1) + DSS(2)*DSS(2) + DSS(3)*DSS(3)
                  R2     = 1.D0 / DSS2
                  ABSR   = SQRT(DSS2)
                  ENERGY = ENERGY + CH2O*CHARGE(I)*CHARGE(J)/ABSR

                  IF (GTEST) THEN

!     DVDR = DVDR/R
                     DVDR       = -CH2O*CHARGE(I)*CHARGE(J)*R2/ABSR 
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

      END SUBROUTINE NEWTIP 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP1(CHARGE, C12, C6, CH2O)
!     TIPS WATER    

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE  
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETA, ROH, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETA)*ROH
      SITE(2,3) = COS(0.5D0*THETA)*ROH

      SITE(3,1) = 0.D0
      SITE(3,2) = -SIN(0.5D0*THETA)*ROH
      SITE(3,3) = COS(0.5D0*THETA)*ROH
  
      CHARGE(:) = (/-0.8D0, 0.4D0, 0.4D0/)
      C6        = 2510.4D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
      C12       = 2426720.D0
      CH2O      = 1389.354848D0 ! CONVERSION FACTOR FOR COULOMB ENERGY

      M(:)  = (/16.D0, 1.D0, 1.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP1

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP2(CHARGE, C12, C6, CH2O)
!     TIPS2 WATER    
    
      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE  
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETA, ROH, ROM, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      ROM   = 0.15D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETA)*ROH
      SITE(2,3) = COS(0.5D0*THETA)*ROH

      SITE(3,1) = 0.D0
      SITE(3,2) = -SIN(0.5D0*THETA)*ROH
      SITE(3,3) = COS(0.5D0*THETA)*ROH
  
      SITE(4,1) = 0.D0
      SITE(4,2) = 0.D0
      SITE(4,3) = ROM

      CHARGE(:) = (/0.D0, 0.535D0, 0.535D0, -1.07D0/)
      C6        = 2510.4D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
      C12       = 2907880.D0
      CH2O      = 1389.354848D0 ! CONVERSION FACTOR FOR COULOMB ENERGY

      M(:)  = (/16.D0, 1.D0, 1.D0, 0.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP2

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP3(CHARGE, C12, C6, CH2O)
!     TIP3P WATER    
    
      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE  
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETA, ROH, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETA)*ROH
      SITE(2,3) = COS(0.5D0*THETA)*ROH

      SITE(3,1) = 0.D0
      SITE(3,2) = -SIN(0.5D0*THETA)*ROH
      SITE(3,3) = COS(0.5D0*THETA)*ROH
  
      CHARGE(:) = (/-0.834D0, 0.417D0, 0.417D0/)
      C6        = 2489.48D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
      C12       = 2435088.D0
      CH2O      = 1389.354848D0 ! CONVERSION FACTOR FOR COULOMB ENERGY

      M(:)  = (/16.D0, 1.D0, 1.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP3

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP4(CHARGE, C12, C6, CH2O)
!     TIP4P WATER    
    
      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE  
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETA, ROH, ROM, PI

      PI    = 4.D0*DATAN(1.D0)
      ROH   = 0.9572D0
      ROM   = 0.15D0
      THETA = 104.52D0
      THETA = PI*THETA/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETA)*ROH
      SITE(2,3) = COS(0.5D0*THETA)*ROH

      SITE(3,1) = 0.D0
      SITE(3,2) = -SIN(0.5D0*THETA)*ROH
      SITE(3,3) = COS(0.5D0*THETA)*ROH
  
      SITE(4,1) = 0.D0
      SITE(4,2) = 0.D0
      SITE(4,3) = ROM

      CHARGE(:) = (/0.D0, 0.52D0, 0.52D0, -1.04D0/)
      C6        = 2552.24D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
      C12       = 2510.4D3
      CH2O      = 1389.354848D0 ! CONVERSION FACTOR FOR COULOMB ENERGY

      M(:)  = (/16.D0, 1.D0, 1.D0, 0.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP4

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE DEFTIP5(CHARGE, C12, C6, CH2O)
!     TIP5P WATER    
    
      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE  
      
      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: CHARGE(NRBSITES), M(NRBSITES), MASS, CM(3), C12, C6, CH2O
      DOUBLE PRECISION :: THETAH, THETAM, ROH, ROM, PI

      PI     = 4.D0*DATAN(1.D0)
      ROH    = 0.9572D0
      ROM    = 0.70D0
      THETAH = 104.52D0
      THETAM = 109.47D0
      THETAH = PI*THETAH/180.D0
      THETAM = PI*THETAM/180.D0

!     THE REFERENCE GEOMETRY IS ON THE Y-Z PLANE

      SITE(1,1) = 0.D0
      SITE(1,2) = 0.D0
      SITE(1,3) = 0.D0

      SITE(2,1) = 0.D0
      SITE(2,2) = SIN(0.5D0*THETAH)*ROH
      SITE(2,3) = COS(0.5D0*THETAH)*ROH

      SITE(3,1) = 0.D0
      SITE(3,2) =-SIN(0.5D0*THETAH)*ROH
      SITE(3,3) = COS(0.5D0*THETAH)*ROH
  
      SITE(4,1) = SIN(0.5D0*THETAM)*ROM 
      SITE(4,2) = 0.D0
      SITE(4,3) =-COS(0.5D0*THETAM)*ROM 

      SITE(5,1) =-SIN(0.5D0*THETAM)*ROM 
      SITE(5,2) = 0.D0
      SITE(5,3) =-COS(0.5D0*THETAM)*ROM
 
     CHARGE(:) = (/0.D0, 0.241D0, 0.241D0, -0.241D0, -0.241D0/)
      C6        = 2470.012857D0 ! LJ COEFFICIENTS IN KJ/MOL ANGSTROM**6 OR ANGSTROM**12
      C12       = 2278383.244D0
      CH2O      = 1389.354848D0 ! CONVERSION FACTOR FOR COULOMB ENERGY

      M(:)  = (/16.D0, 1.D0, 1.D0, 0.D0, 0.D0/)
      CM(:) = 0.D0; MASS = 0.D0
      DO I = 1, NRBSITES
         CM(:) = CM(:) + M(I)*SITE(I,:)
         MASS = MASS + M(I)
      ENDDO
      CM(:) = CM(:)/MASS
      DO I = 1, NRBSITES
         SITE(I,:) = SITE(I,:) - CM(:)
      ENDDO

      END SUBROUTINE DEFTIP5

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE VIEWNEWTIP()

      USE COMMONS, ONLY: NATOMS, NRBSITES, SITE, NSAVE
      USE QMODULE

      IMPLICIT NONE

      INTEGER          :: I, J1, J2, J3, J5, J7
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3)
      LOGICAL          :: GTEST

      OPEN(UNIT=26, FILE='NEWTIP.XYZ', STATUS='UNKNOWN')

      GTEST = .FALSE. 

      DO J1 = 1, NSAVE

         WRITE(26,'(I6)') (NATOMS/2)*3
         WRITE(26,10) J1, QMIN(J1), FF(J1)
10       FORMAT('ENERGY OF MINIMUM ',I6,'=',F20.10,' FIRST FOUND AT STEP ',I8)

         DO J3 = 1, NATOMS/2

            J5   = 3*J3
            J7   = 3*NATOMS/2 + J5
            P(:) = QMINP(J1,J7-2:J7)

            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

            DO J2 = 1, 3

               RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + MATMUL(RMI(:,:),SITE(J2,:))
               IF (J2 == 1) THEN
                  WRITE(26,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSE
                  WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ENDIF

            ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWNEWTIP
