! Calculates Energy and gradients for the patch-antipatch potential, using the rigid body angle axis framework described in
! 'Simulations of rigid bodies in an angle axis framework', Dwaipayan Chakrabarti and David Wales, 
! Phys. Chem. Chem. Phys., 2009, 11, 1970-1976
!
! X: the positions and orientations of the bodies
! G: the gradients of the potential energy surface for changes in positions and orientations of the bodies
! ENERGY: the potential energy for the configuration stored in X
! GTEST: logical, true if gradients need to be calculated
      SUBROUTINE PAP (X, G, ENERGY, GTEST)

! NATOMS: twice the number of bodies, one for position and one for orientation
! NRBSITES: the number of patches and antipatches per body
! RBSTLA: the directions of the patch and antipatch vectors
! PAPALP: the pap potential 'alpha' parameter, controls the range of the isotropic LJ potential
! PAPS: the pap potential 's' parameter, controls the smoothness of the patch-antipatch attraction
! PAPCD: the pap potential 'cos(delta)' parameter, controls the width of the patches
! PAPEPS: the pap potential 'epsilon' parameter, controls the depth of the patch-antipatch attraction
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSTLA, PAPALP, PAPS, PAPCD, PAPEPS

      IMPLICIT NONE

! I,J J1-8: indices and counters to iterate over
! NMOL: the number of bodies
! OFFSET: the indicial offset to the start of orientations in X
      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, J7, J8, NMOL, OFFSET

      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS) 

! DVDR: derivative of LJ potential with respect to the distance between bodies I and J, divided by that distance
! DSS2: square modulus of RIJ, so the distance between bodies I and J squared
! EXP6: R2D**6, denominator in LJ potential
! EXP3: R2D**3, denomiator in LJ potential
! R2D: 1/(DSS2-1), relevant to LJ potential
! PREFAC: factor in the LJ potential, 4/PAPALP**2
! ABSR: square root of DSS2, so the absolute distance between bodies I and J
! PI: the mathematical constant pi
      DOUBLE PRECISION :: ENERGY, DVDR, DSS2, EXP6, EXP3, R2D, PREFAC, ABSR, PI 

! LJSIGMASQ: the sigma value in the LJ potential squared
! LJN: the exponent in the LJ potential
! R2: 1/DSS2, relevant to LJ potential
! R2LJN: R2**LJN, relevant to LJ potential
! RLJN: R2**(LJN/2), relevant to LJ potential 
      DOUBLE PRECISION :: LJSIGMASQ, LJN, R2, R2LJN, RLJN

! RIJ: the vector displacement between bodies I and J
! NR: normalised RIJ
! P: the rotation vector for some body
! WP: the patch-antipatch attraction potential
! DELR: RIJ - LAMBDA, argument for patch-antipatch attraction
! ANGFAC: factor in the angular potential phi, PI/(1-PAPCD)
! DWPDR: derivative of WP with respect to the ditance between bodies I and J
      DOUBLE PRECISION :: RIJ(3), NR(3), P(3), WP, DELR, ANGFAC, DWPDR 

! DOTI: angle between orientation of patch I and the vector displacement between two bodies
! DOTJ: angle between orientation of patch J and the vector displacement between two bodies
! ARGI: argument for PHII in the angular potential
! ARGJ: argument for PHIJ in the angular potential
! PHII: function in the angular potential
! PHIJ: function in the angular potential
      DOUBLE PRECISION :: DOTI, DOTJ, ARGI, ARGJ, PHII, PHIJ

! E: the orientations of patches and antipatches in the lab frame
      DOUBLE PRECISION :: E(NATOMS*NRBSITES/2,3)

! DEk: matrix product of the derivative of the rotation matrix with respect to kth component of the rotation vector and the
!      position of each site in the reference geometry, as in equation (12) of the paper
      DOUBLE PRECISION :: DE1(NATOMS*NRBSITES/2,3), DE2(NATOMS*NRBSITES/2,3), DE3(NATOMS*NRBSITES/2,3) 

! RMI: the rotation matrix for body I
! DRMIk: the derivative RMI with respect to the kth component of the rotation vector
      DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)

! DDOTIDR:  derivative of DOTI wrt a change in position of body I
! DDOTJDR:  derivative of DOTJ wrt a change in position of body I
! DPHIIDR:  derivative of PHII wrt a change in position of body I
! DPHIJDR:  derivative of PHIJ wrt a change in position of body I
! DPHIIDPI: derivative of PHII wrt a change in orientation of body I (note DPHIIDPJ would be zero)
! PPHIJDPJ: derivative of PHIJ wrt a change in orientation of body J (note DPHIJDPI would be zero)
      DOUBLE PRECISION :: DDOTIDR(3), DDOTJDR(3), DPHIIDR(3), DPHIJDR(3), DPHIIDPI(3), DPHIJDPJ(3)

! LAMBDA: the range at which patch-antipatch attraction starts to decrease, and the minimum of the LJ potential
! INVS: PAPS**-1, LAMBDA+INVS is the range at which the patch-antipatch attraction becomes zero
      DOUBLE PRECISION :: LAMBDA, INVS 

      LOGICAL          :: GTEST

! Set up pap potential, subroutine at the bottom of this file
      CALL DEFPAP()

! Define useful factors for calculating the potential
      PI      = 4.D0*DATAN(1.D0)
      ANGFAC  = PI/(1.D0 - PAPCD)
      INVS    = 1.D0/PAPS
      PREFAC  = 4.D0/(PAPALP*PAPALP)
      LAMBDA  = SQRT(1.D0 + (2.D0/PAPALP)**(1.D0/3.D0))

! Initialise energy and gradients
      ENERGY  = 0.D0
      IF (GTEST) G(:) = 0.D0

! Find NMOL as the actual number of bodies, and OFFSET as the start of rotational coordinates
      NMOL    = NATOMS/2
      OFFSET  = 3*NMOL

! J1 Loop over all bodies  
      DO J1 = 1, NMOL

! Set J3 as the index of the third translational coordinate of the J1th body, and J5 as the index of the third rotational
! coordinate of the J1th body
! Set P to be the rotation vector of body J1
         J3 = 3*J1
         J5 = OFFSET + J3
         P  = X(J5-2:J5)

! Calculate the rotation matrix and derivatives thereof
         CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

! J2 Loop over all patches and antipatches for body J1
         DO J2 = 1, NRBSITES

! Set J4 as the index of the J2th patch of body J1
            J4        = NRBSITES*(J1-1) + J2
! Calculate E from the rotation matrix acting on the patch orientation in the body frame
            E(J4,:)   = MATMUL(RMI(:,:),RBSTLA(J2,:))

            IF (GTEST) THEN

! If derivatives are required, calculate the kth derivative of the rotation matrix for body J1 acting on the orientation of
! the J2th patch in the body frame
               DE1(J4,:) = MATMUL(DRMI1(:,:),RBSTLA(J2,:))
               DE2(J4,:) = MATMUL(DRMI2(:,:),RBSTLA(J2,:))
               DE3(J4,:) = MATMUL(DRMI3(:,:),RBSTLA(J2,:))

            ENDIF

         ENDDO

      ENDDO

! J1 Loop over all bodies, except the last, to perform the double sum over bodies without double counting
      DO J1 = 1, NMOL - 1  

! Set J3 as the index of the third translational coordinate of body J1, and J5 as the index of the third rotational
! coordinate of body J1
         J3 = 3*J1
         J5 = OFFSET + J3

! J2 Loop over bodies index greater than J1, to perform the double sum over bodies without double counting
         DO J2 = J1 + 1, NMOL

! Set J4 as the index of the third translational coordinate of body J1, and J6 as the index of the third rotational
! coordinate of body J1
            J4 = 3*J2
            J6 = OFFSET + J4

! Calculate RIJ as the vector displacement from body J2 to J1, then calculate the absolute distance bewteen I and J
            RIJ(:) = X(J3-2:J3) - X(J4-2:J4)
            DSS2   = DOT_PRODUCT(RIJ(:),RIJ(:))
            ABSR   = SQRT(DSS2)

! Calculate denomiators relevant to the LJ potential
! Add the LJ contribution to the potential energy from the interaction of bodies J1 and J2
! LJ potential from proposed PAP potential description, generates structures with overlapping cores
!            R2D    = 1.D0/(DSS2 - 1.D0)
!            EXP6   = R2D**6
!            EXP3   = R2D**3
!            ENERGY = ENERGY + PREFAC*(EXP6 - PAPALP*EXP3)

! LJ 2n-n potential, with epsilon equal to one
! Sigma must be adjusted to match the proposed potential above
! LJN chosen to give a reasonable match to the above potential
            R2 = 1.D0/DSS2
            LJN = 23
            LJSIGMASQ = 1.D0 + 1.D0/PAPALP**(1.D0/3.D0)
            RLJN = (LJSIGMASQ*R2)**(LJN/2.D0)
            R2LJN = RLJN**2
            ENERGY = ENERGY + 4*(R2LJN - RLJN)

            IF (GTEST) THEN

! Add the LJ contribution to the gradient from the change in position of J1 due to J2, and vice versa
! Simple due to isotropic LJ potential, and single site bodies
! Note that DVDR is the derivative of potential wrt distance divided by the distance
! LJ potential from proposed PAP potential description, generates structures with overlapping cores
!               DVDR = 2.D0*PREFAC*(-6*EXP6 + 3.D0*PAPALP*EXP3)*R2D

! LJ 2n-n potential, with epsilon equal to one
               DVDR = 4.D0*LJN*R2*(RLJN-2*R2LJN)

               G(J3-2:J3) = G(J3-2:J3) + DVDR*RIJ(:)
               G(J4-2:J4) = G(J4-2:J4) - DVDR*RIJ(:)

            ENDIF

! Now need to calculate angular part of the potential
! Find DELR, argument for patch-antipatch attraction
            DELR = ABSR - LAMBDA

! Set specific cases for the patch-antipatch potential and its gradient
            IF (DELR < 0.D0) THEN
! Distance less than LAMBDA, so full attraction
               WP =-1.D0 
               DWPDR = 0.D0           
            ELSE IF (DELR <= INVS) THEN
! Distance is between LAMBDA and LAMBDA+INVS, so potential varies from -1 to 0 
               WP =-0.5D0*(1.D0 + COS(PI*DELR*PAPS))
               DWPDR = 0.5D0*PI*PAPS*SIN(PI*DELR*PAPS) 
            ELSE
! Distance is greater than LAMBDA+INVS, so no attraction
               WP = 0.D0
               DWPDR = 0.D0
            ENDIF

! Find the normalised vector parallel to the displacement between bodies J1 and J2
            NR(:) = RIJ(:)/ABSR

! I Loop over patches in body J1
            DO I = 1, NRBSITES/2    ! for patch-antipatch, I only loops over patches
!            DO I = 1, NRBSITES      ! for only one kind of patches, I loops over all patches

! Set J7 as the index of the Ith patch in body J1
               J7    = NRBSITES*(J1-1) + I

! J Loop over patches in body J2
               DO J = NRBSITES/2 + 1, NRBSITES   ! for patch-antipatch, J only loops over antipatches
!               DO J = 1, NRBSITES                ! for only one kind of patches, J loops over all patches

! Set J8 as the index of Jth (anti)patch in body J2
                  J8 = NRBSITES*(J2-1) + J

! Find the angles between the orientation of patch J7 and NR, and between the orientation of (anti)patch J8 and NR
! Negative for DOTI, since RIJ is displacement from J2 to J1
                  DOTI =-DOT_PRODUCT(E(J7,:),NR(:))
                  DOTJ = DOT_PRODUCT(E(J8,:),NR(:))

! Calculate the values of PHI for the angular potential, and their derivatives wrt a change in ABSR
                  IF (DOTI < PAPCD) THEN
! Angle greater than patch width, so no attraction
                     PHII =-1.D0
                  ELSE
! Angle less than patch width, so some attraction
                     ARGI = ANGFAC*(DOTI-PAPCD)
                     PHII = -COS(ARGI) ! changed to minus
                  ENDIF 
                  IF (DOTJ < PAPCD) THEN
! Angle greater than patch width, so no attraction
                     PHIJ =-1.D0
                  ELSE
! Angle less than patch width, so some attraction
                     ARGJ = ANGFAC*(DOTJ-PAPCD)
                     PHIJ = -COS(ARGJ) ! changed to minus
                  ENDIF 

! Add the patch-antipatch attraction to the potential energy
                  ENERGY = ENERGY + 0.25D0*PAPEPS*(1.D0 + PHII)*(1.D0 + PHIJ)*WP

                  IF (GTEST) THEN

! Need to find the derivatives of the patch-antipatch attraction wrt translational and rotational coordinates
                      IF ((DOTI >= PAPCD) .AND. (DOTJ >= PAPCD)) THEN

! Calculate the derivates of the angle between the patch orientation and RIJ for J7 and J8 wrt a change in position of J1
                         DDOTIDR(:) = -DOTI*RIJ(:)/DSS2 - E(J7,:)/ABSR
                         DDOTJDR(:) = -DOTJ*RIJ(:)/DSS2 + E(J8,:)/ABSR

! Find the derivatives of the PHIs wrt a change in position of J1
                         DPHIIDR(:) = ANGFAC*SIN(ARGI)*DDOTIDR(:) ! changed to plus
                         DPHIJDR(:) = ANGFAC*SIN(ARGJ)*DDOTJDR(:) ! changed to plus

! Add the contribution to the gradient of a change in position of J1 due to J2, and vice versa
                         G(J3-2:J3) = G(J3-2:J3) + 0.25D0*PAPEPS*((1.D0 + PHIJ)*WP*DPHIIDR(:)                  &
     &                              + (1.D0 + PHII)*WP*DPHIJDR(:) +(1.D0+PHII)*(1.D0+PHIJ)*DWPDR*NR(:))
                         G(J4-2:J4) = G(J4-2:J4) - 0.25D0*PAPEPS*((1.D0 + PHIJ)*WP*DPHIIDR(:)                  &
     &                              + (1.D0 + PHII)*WP*DPHIJDR(:) +(1.D0+PHII)*(1.D0+PHIJ)*DWPDR*NR(:))

                      ENDIF

                      IF (DOTI >= PAPCD) THEN
! Find the derivatives of PHII wrt a change in each of the rotational coordinates of body I
                         DPHIIDPI(1) = -ANGFAC*SIN(ARGI)*DOT_PRODUCT(NR(:),DE1(J7,:)) ! changed to minus
                         DPHIIDPI(2) = -ANGFAC*SIN(ARGI)*DOT_PRODUCT(NR(:),DE2(J7,:)) ! changed to minus
                         DPHIIDPI(3) = -ANGFAC*SIN(ARGI)*DOT_PRODUCT(NR(:),DE3(J7,:)) ! changed to minus
! Add the contribution to the gradient of a change in orientation of J1 due to J2
! Simpler because WP does not depend on orientation, and PHIJ does not depend on the orientation of J1
                         G(J5-2:J5) = G(J5-2:J5) + 0.25D0*PAPEPS*(1.D0 + PHIJ)*WP*DPHIIDPI(:)
                      ENDIF

                      IF (DOTJ >= PAPCD) THEN
! Find the derivatives of PHIJ wrt a change in each of the rotational coordinates of body J
                         DPHIJDPJ(1) = ANGFAC*SIN(ARGJ)*DOT_PRODUCT(NR(:),DE1(J8,:)) ! changed to plus
                         DPHIJDPJ(2) = ANGFAC*SIN(ARGJ)*DOT_PRODUCT(NR(:),DE2(J8,:)) ! changed to plus
                         DPHIJDPJ(3) = ANGFAC*SIN(ARGJ)*DOT_PRODUCT(NR(:),DE3(J8,:)) ! changed to plus
! Add the contribution to the gradient of a change in orientation of J2 due to J1
! Simpler because WP does not depend on orientation, and PHII does not depend on the orientation of J2
                         G(J6-2:J6) = G(J6-2:J6) + 0.25D0*PAPEPS*(1.D0 + PHII)*WP*DPHIJDPJ(:)
                      ENDIF

                  ENDIF

               ENDDO

            ENDDO
 
         ENDDO

      ENDDO

      END SUBROUTINE PAP
 
!     ----------------------------------------------------------------------------------------------

! Defines pap potential by specifying the positions of patches and antipatches
! The first half of the vectors are patches, the second half anitpatches 
      SUBROUTINE DEFPAP()
! RBSTLA: vectors specifying the directions of patches and antipatches
! NRBSITES: the total number of patches and antipatches per body
      USE COMMONS, ONLY: RBSTLA, NRBSITES

      IMPLICIT NONE

      IF (NRBSITES == 2) THEN
! Create a patch along the positive z axis and an antipatch along the negative z axis
        RBSTLA(1,:) = (/ 0.D0, 0.D0, 1.D0/)
        RBSTLA(2,:) = (/ 0.D0, 0.D0,-1.D0/)

      ELSE IF (NRBSITES == 4) THEN
! Create a tetrahedral arrangement of patches and antipatches
        RBSTLA(1,:)= 1.D0/SQRT(3.D0)*(/  1.D0,  1.D0,  1.D0/)
        RBSTLA(2,:)= 1.D0/SQRT(3.D0)*(/ -1.D0, -1.D0,  1.D0/)
        RBSTLA(3,:)= 1.D0/SQRT(3.D0)*(/ -1.D0,  1.D0, -1.D0/)
        RBSTLA(4,:)= 1.D0/SQRT(3.D0)*(/  1.D0, -1.D0, -1.D0/)

      ENDIF

      END SUBROUTINE DEFPAP

!     ----------------------------------------------------------------------------------------------

! Writes positions to file in a format readable by XMakeMol
      SUBROUTINE VIEWPAP()

! NATOMS: twice the number of bodies
! NRBSITES: the number of patches and antipatches per body
! RBSTLA: the patch and antipatch directions in the body frame
! NSAVE: the number of different configurations to be written
      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSTLA, NSAVE
! Global variable from QMODULE:
! QMIN: the energies of the stored configurations
! FF: the steps at which the configurations were first found
! QMINP: two dimensional array storing the positions, then orientations, for each configuration
      USE QMODULE

      IMPLICIT NONE

! I, J1-7: indices and counters to iterate over
      INTEGER          :: I, J1, J2, J3, J5, J7

! RMI: the rotation matrix for a body
! DRMI: dummy variable
! P: the rotation vector for a body
! RBCOORDS: the position of a patch or antipatch for a body
! LFCTR: length factor, distance of patches from centre of body
      DOUBLE PRECISION :: RMI(3,3), DRMI(3,3), P(3), RBCOORDS(3), LFCTR

! GTEST: indicates whether to calculate gradients or not
      LOGICAL          :: GTEST

! Choose vale of LFCTR to look pretty
      LFCTR = 0.4D0

! Open writer to file 'pap.xyz'
      OPEN(UNIT=26, FILE='pap.xyz', STATUS='UNKNOWN')

! Gradients are not required
      GTEST = .FALSE.

! J1 Loop over each configuration
      DO J1 = 1, NSAVE

! Write the number of cores, patches and antipatches 
         WRITE(26,'(I6)') (NATOMS/2)*(NRBSITES+1)
! Write the energy and step at which the configuration was first found
         WRITE(26,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)

! J3 Loop over each body
         DO J3 = 1, NATOMS/2

! Set J5 as the index of the third translational coordinate of body J3, and J7 as the third rotational coordinate
            J5   = 3*J3
            J7   = 3*NATOMS/2 + J5
! Set P as the rotational coordinates
            P(:) = QMINP(J1,J7-2:J7)

! Find the rotation matrix
! Since GTEST is false, DRMI is a dummy variable
            CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

! Write the position of the cores of the bodies
! Display as an Oxygen atom
            WRITE(26,'(A4,3F20.10)') 'O', QMINP(J1,J5-2), QMINP(J1,J5-1), QMINP(J1,J5)

! J2 Loop over patches and antipatches within body J3
            DO J2 = 1, NRBSITES

! Calculate the position of the patch J2 in the lab frame
              RBCOORDS(1:3) = QMINP(J1,J5-2:J5) + LFCTR*MATMUL(RMI(:,:),RBSTLA(J2,:))

! For patches and antipatches different, use this block
               IF (J2 <= NRBSITES/2) THEN
! If J2 is a patch, display as a fluorine atom
                  WRITE(26,'(A4,3F20.10)') 'F', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ELSE
! If J2 is an antipatch, display as a carbon atom
                  WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
               ENDIF

! For patches and antipatches the same, use this block
! Display patches and antipatches both as fluorine atoms
!              WRITE(26,'(A4,3F20.10)') 'F', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)

            ENDDO

         ENDDO

      ENDDO

      CLOSE (UNIT=26)

      END SUBROUTINE VIEWPAP

