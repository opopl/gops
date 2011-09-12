!|gd351>

SUBROUTINE ASAOOSPOT (X, G, ENERGY, GTEST)

  USE COMMONS, ONLY: NATOMS

  IMPLICIT NONE
  
  INTEGER          :: J, J1, J2, J3
  DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
  DOUBLE PRECISION :: ENERGY, R, RSQ, ARG, V, DV
  DOUBLE PRECISION :: RI(NATOMS,3), RIJ(3)
  LOGICAL          :: GTEST
  DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D0
  DOUBLE PRECISION :: AT1, AS, ASSQ, HW1, VLV0, VLK0, VCV0, VCR0, VSA0, VSA1, VSA2, VSA3
  DOUBLE PRECISION :: SIGMAC, SIGMAP, Q, ZP, BETA, SIGMAPSQ, SIGMACSQ, RANGESQ, ONEPLUSQCB, PREFACTOR 
  

  !hard core approximation:
  !           r < 0.80000: ~ linear
  ! 0.80000 < r < AS     : ~ r^-12
  ! AS      < r < 0.99500: ~ ArcTanh[(1.0-r)*AT1] + Asakura-Oosawa(1.0)
  ! 0.99500 < r < 1.00000: hard core and attractive potential smoothly merged by polynomial of grade 3
  !attractive potential
  ! 1.00000 < r < 1.10000: Asakura-Oosawa potential
  ! 1.10000 < r          : 0 

  !parameters change when ZP, BETA, SIGMAC or SIGMAP are changed


  AT1 = 50.D0 !!higher AT1 -> steeper hardcore approx
  AS = 1.00001D0 - 1.D0/AT1 
  ASSQ = AS**2

  !parameters for smoothing function VS (polynomial of grade 3)
  !VS(0.995)  = VATH(0.995)
  !VS'(0.995) = VATH(0.995)
  !VS(1.00)   = VAO(1.00)
  !VS'(1.00)  = VAO'(1.00)
  !!VSA0 = -2.5864432352065337D6 !-4.0694693576136194D6
  !!VSA1 =  7.7858744209950390D6 ! 1.2254219988133736D7
  !!VSA2 = -7.8124381430060350D6 !-1.2300050910062172D7
  !!VSA3 =  2.6130061194594866D6! 4.1152994417840126D6

  HW1 = 0.5D0*LOG((1.D0+0.005D0*AT1)/(1.D0-0.005D0*AT1))

  VSA0 = -653154.4451401674D0 + (39800.00008137438D0*AT1)/(1.D0 - 0.000025D0*AT1**2) - 1.5880000032516249D7*HW1 
  VSA1 =  1.9660251800061788D6 - (119600.00024473826D0*AT1)/(1.D0 - 0.000025D0*AT1**2) + 4.776000009779455D7*HW1
  VSA2 = -1.9726060312274103D6 + (119800.00024535337D0*AT1)/(1.D0 - 0.000025D0*AT1**2) - 4.788000009804034D7*HW1
  VSA3 =  659734.4586033578D0 - (40000.00008198949D0*AT1)/(1.D0 - 0.000025D0*AT1**2) + 1.6000000032762045D7*HW1


  !!these have to be calculated when changing AT1

  !parameters for core repulsion VC(R) = VCV0 + (R-VCR0)^(-12)
  !VATH(R) = ArcTanh(-50*(R-1.00))+VAO(1.00) 
  !VC(AS) = VATH(AS)
  !VC'(AS) = VATH'(AS)
  !!VCV0 = -2.1925907471026258D3 !!value for AT1=100
  !!VCR0 =  4.633407007430528D-1 !!value for AT1=100
  VCV0 = -2.1917372506987317D3 !!value for AT1=50
  VCR0 =  4.533305685892997D-1 !!value for AT1=50
  !!VCV0 = -2191.137307077289D0  !!value for AT1=25
  !!VCR0 =  0.4333255033866304D0 !!value for AT1=25
 
  !!these have to be calculated when changing AT1

  !parameters for core repulsion VL(R) = VLV0 + VLK0*R
  !VL(0.8) = VC(0.8)
  !VL'(0.8) = VC'(0.8)
  !!VLV0 =     1.3921755730646929D7 !!value for AT1=100 
  !!VLK0 =  -1.6815246938095480D7   !!value for AT1=100  
  VLV0 =  9.5208471705023050D6    !!value for AT1=50
  VLK0 =  -1.1488917740863195D7   !!value for AT1=50
  !!VLV0 =  4.599326776700881D6     !!value for AT1=25 
  !!VLK0 =  -5.540284774711806D6    !!value for AT1=25 
 


  SIGMAC=1.D0
  SIGMAP=1.D-1

  Q=SIGMAP/SIGMAC

  ZP=1.D2
  BETA=1.D0

  SIGMACSQ=SIGMAC**2
  SIGMAPSQ=SIGMAP**2
  RANGESQ=(SIGMAC+SIGMAP)**2

  ONEPLUSQCB=(1.D0+Q)**3
  PREFACTOR=-PI*ZP*ONEPLUSQCB*(SIGMAP**3)/(6.D0*BETA*Q**3)


  ENERGY = 0.D0
  IF (GTEST) G(:) = 0.D0
  
  DO J1 = 1, NATOMS
    
    J3 = 3*J1
    RI(J1,1:3) = X(J3-2:J3)
       
  END DO


  DO J1 = 1, NATOMS

    DO J2 = J1+1, NATOMS

      RIJ(:) = RI(J1,:) - RI(J2,:)
      RSQ = DOT_PRODUCT(RIJ(:),RIJ(:))

      IF (RSQ.LE.0.64D0*SIGMACSQ) THEN
      !core - linear
        R = SQRT(RSQ)
        
        ENERGY = ENERGY + (VLV0 + VLK0 * R)  

        IF (GTEST) THEN
          !calculate derivatives
          DV = VLK0 / R   !factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 


      ELSEIF ((RSQ.GT.0.64D0*SIGMACSQ).AND.(RSQ.LE.ASSQ*SIGMACSQ)) THEN
      !core - R^-12
        R = SQRT(RSQ)
        
        ENERGY = ENERGY + (VCV0 + (R-VCR0)**(-12))
!        print*,'r^-12',(VCV0 + (R-VCR0)**(-12))  

        IF (GTEST) THEN
          !calculate derivatives
          DV = -12.D0 / ( (R-VCR0)**(13) * R )   !factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 

      ELSEIF ((RSQ.GT.ASSQ*SIGMACSQ).AND.(RSQ.LE.0.990025D0*SIGMACSQ)) THEN
      !core - ArcTanh
        R = SQRT(RSQ)
        ARG = (-AT1)*(R-1.D0)
        V = 0.5D0*LOG((1+ARG)/(1-ARG)) - 0.837758D0 ! = ArcTanh(ARG)      
!        print*,'ArcTanh',V      
 
        ENERGY = ENERGY + V

        IF (GTEST) THEN
          !calculate derivatives
          DV = -AT1 / ( (1.D0-ARG**2) * R )  !factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 

      ELSEIF ((RSQ.GT.0.990025D0*SIGMACSQ).AND.(RSQ.LE.SIGMACSQ)) THEN
      !smooth merging
        R = SQRT(RSQ)
        V = VSA0 + VSA1 * R + VSA2 * RSQ + VSA3 * RSQ*R

        ENERGY = ENERGY + V
!        print*,'merging',V    

        IF (GTEST) THEN
          !calculate derivatives
          DV = VSA1 / R + 2.D0 * VSA2 + 3.D0 * VSA3 * R!factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 

      ELSEIF ((RSQ.GT.SIGMACSQ).AND.(RSQ.LE.RANGESQ)) THEN
      !Asakura-Oosawa attraction
        R = SQRT(RSQ)
        V = PREFACTOR * (1.D0 - (3.D0*R)/(2.D0*SIGMAC*(1.D0+Q)) + (RSQ*R)/(2.D0*ONEPLUSQCB*SIGMAC**3))

        ENERGY = ENERGY + V
!        print*,'asaoos',R,V    

        IF (GTEST) THEN
          !calculate derivatives
          DV = PREFACTOR * ( - 3.D0/(R*2.D0*SIGMAC*(1.D0+Q)) + 3.D0*R/(2.D0*ONEPLUSQCB*SIGMAC**3))!factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 

      END IF

    END DO
  END DO

END SUBROUTINE ASAOOSPOT



!<gd351|
