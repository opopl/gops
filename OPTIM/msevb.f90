!   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
!   COPYRIGHT (C) 1999-2006 DAVID J. WALES
!   THIS FILE IS PART OF OPTIM.
!   
!   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
!   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
!   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
!   (AT YOUR OPTION) ANY LATER VERSION.
!   
!   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
!   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
!   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
!   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
!   
!   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
!   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
!   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
!
!###################################################################
!################# SUBROUTINE MSEVB ################################
!###################################################################

!###################################################################

! 02/04/04 CHANGED JACOBI ROUTINE FOR LAPACK DSYEVR
! INDIVIDUAL RESULTS ARE CONSISTENT TO 1E-12 BUT THERE SEEM TO BE DIFFERENCES IN MARKOV CHAINS NONE THE LESS

     SUBROUTINE MSEVB(NCOORDS, SYSTEMCOORDS, ASSIGNVBSTATES, ENERGY, ONLYASSIGNSTATES)
     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

!###################################################################
!
! PARAMETER SET USED IS FROM:
!    U. SCHMITT AND G. A. VOTH, J. CHEM. PHYS. 111, 1999, 9361
! AND NOT THE EARLIER SET, NOR THE LATER ONE, OR EVEN THE BUGGY ONE.
!
! CALCULATE THE MSEVB ENERGY OF A GIVEN PROTONATED WATER CLUSTER.
!
! STEPS: (I) ASSIGN THE POSSIBLE EVB STATES IN THE CLUSTER, (II)
! TRANSFER THE EVB BONDING ARRANGEMENTS INTO COORDINATE VECTORS
! XX(I,J) ETC. (III) CALCULATE THE HAMILTONIAN ELEMENTS BY MEANS
! OF ONLY THE NUCLEAR POSITIONS, AND FINALLY DIAGONALIZING THE
! HAMILTONIAN TO YIELD THE EIGENVECTORS AND EIGENVALUES.
!
!###################################################################

! SUBROUTINE ARGUMENTS

     INTEGER, INTENT(IN) :: NCOORDS  
     DOUBLE PRECISION, INTENT(IN) :: SYSTEMCOORDS(NCOORDS)
     LOGICAL, INTENT(IN) :: ASSIGNVBSTATES
     DOUBLE PRECISION, INTENT(OUT) :: ENERGY
     LOGICAL, OPTIONAL, INTENT(IN) :: ONLYASSIGNSTATES   ! ONLY ASSIGN THE VB STATES, DO NOTHING ELSE

! LOCAL VARIABLES

     INTEGER I,J,NROT

       DOUBLE PRECISION MIN

! VARIABLES FOR LAPACK ROUTINE

     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: D
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: V, E
!     INTEGER :: NUMEVFOUND
!     INTEGER :: ISUPPZ(2)
!     INTEGER :: LOCALSTATUS
!     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK
!     INTEGER, ALLOCATABLE, DIMENSION(:) :: IWORK

!     DOUBLE PRECISION :: DLAMCH ! EXTERNAL FUNCTION

! STORE THE ATOM COORDINATES IN PSIX, PSIY, PSIZ RATHER THAN THE SINGLE ARRAY COORDS

     J = 0

     DO I = 1,NCOORDS,3
        J = J + 1
        PSIX(J) = SYSTEMCOORDS(I)
        PSIY(J) = SYSTEMCOORDS(I+1)
        PSIZ(J) = SYSTEMCOORDS(I+2)
     ENDDO

! CALCULATE THE INTERATOMIC DISTANCES

     CALL CALCINTERATOMICDISTANCES(NATOMS,PSIX,PSIY,PSIZ,INTERATOMICR)

!###################################################################
! CHECK THAT THERE IS NO ROO VECTORS SMALLER IN MAGNITUDE THAN
! THE RCUT CUTOFF VALUE
! 
!      DO I = 3,NATOMS-4,3
!         DO J = I+3,NATOMS-1,3
!            DX = PSIX(I) - PSIX(J) 
!            DY = PSIY(I) - PSIY(J)   
!            DZ = PSIZ(I) - PSIZ(J)   
!            R = DSQRT(DX**2 + DY**2 + DZ**2)
!            IF (R.LT.RCUT) THEN
!            ENERGY = 1.0D03
!            RETURN
!            ENDIF
!         ENDDO
!      ENDDO
! 
!###################################################################
! ASSIGN THE VB STATES IF REQUIRED, ELSE USE THE ONES IN MEMORY
!
!###################################################################

      IF (ASSIGNVBSTATES) THEN
        CALL ASSIGNVB3
        IF (PRESENT(ONLYASSIGNSTATES).AND.ONLYASSIGNSTATES) RETURN
     ENDIF

! ALLOCATE THE RELEVANT ARRAYS

     ALLOCATE(XX(REDUCED_NUM_EIG,NATOMS), YY(REDUCED_NUM_EIG,NATOMS), ZZ(REDUCED_NUM_EIG,NATOMS))
     ALLOCATE(AIJ(REDUCED_NUM_EIG,REDUCED_NUM_EIG))
     ALLOCATE (D(REDUCED_NUM_EIG), V(REDUCED_NUM_EIG, REDUCED_NUM_EIG), E(REDUCED_NUM_EIG, REDUCED_NUM_EIG))
!     ALLOCATE (D(REDUCED_NUM_EIG), V(REDUCED_NUM_EIG, 1), E(REDUCED_NUM_EIG, REDUCED_NUM_EIG))

! VIJEXCH AND GRSTWFU ARE THE ONLY ARRAYS USED EXTERNALLY, ALL THE REST ARE DEALLOCATED AT THE END OF THE ROUTINE

     IF (ALLOCATED(VIJEXCH)) DEALLOCATE(VIJEXCH)
     ALLOCATE(VIJEXCH(REDUCED_NUM_EIG,REDUCED_NUM_EIG))

     IF (ALLOCATED(GRSTWFU)) DEALLOCATE(GRSTWFU)
     ALLOCATE (GRSTWFU(REDUCED_NUM_EIG))

! SOMETIMES WANT HAM FOR DEBUGGING

!     IF (ALLOCATED(HAM)) DEALLOCATE(HAM)
     ALLOCATE(HAM(REDUCED_NUM_EIG,REDUCED_NUM_EIG))

!###################################################################
! ASSIGN THE CONFIGURATIONS TO XX, THE COORDINATE VECTOR USED IN
! THE SUBSEQUENT CALCULATION OF THE HAMILTONIAN ELEMENTS
!
!###################################################################
 
     DO I = 1,REDUCED_NUM_EIG
        DO J = 1,NATOMS
           XX(I,J) = PSIX(ATMPL(I,J))
           YY(I,J) = PSIY(ATMPL(I,J))
           ZZ(I,J) = PSIZ(ATMPL(I,J))
        ENDDO
     ENDDO

!###################################################################
! BUILD THE HAMILTONIAN: CALL THE ROUTINES FOR CONSTRUCTING EACH 
! ELEMENT HIJ
!
!###################################################################

     CALL BILDHAM

       DO I = 1,REDUCED_NUM_EIG
          DO J = 1,REDUCED_NUM_EIG
             E(I,J) = HAM(I,J)
          ENDDO
       ENDDO

!###################################################################
! DIAGONALIZE THE HAMILTONIAN
!
! FIND THE GROUND STATE ENERGY AND EIGENVECTOR. NOTE THE MONITORING
! OF WHICH EVB STATE IS OCCUPIED IN THIS GEOMETRY
!
!###################################################################

      CALL JACOBI(E,REDUCED_NUM_EIG,D,V,NROT)

! D CONTAINS THE EIGENVALUES
! V CONTAINS THE EIGENVECTORS

!----FIND THE LOWEST EIGENVALUE AND PRINT OUT THE GEOMETRY
      MIN=1.0D10
      GROUND_STATE = 3
      DO I = 1,REDUCED_NUM_EIG
         IF (D(I).LT.MIN) THEN
            MIN = D(I)
            GROUND_STATE = I
         ENDIF 
       ENDDO

!----ASSIGN THE LOWEST ENERGY FROM THE LOWEST EIGENVALUE
      ENERGY = MIN

! RUN LAPACK ROUTINE WHICH CALCULATES ONLY THE LOWEST EIGENVALUE, RATHER THAN THE FULL SPECTRUM
! ONLY NEED ONE TRIANGLE
! ROUTINE OVERWRITES THE HAMILTONIAN WHICH IS WHY WE MUST PASS IT A COPY RATHER THAN THE ORIGINAL

!     DO I = 1,REDUCED_NUM_EIG
!        DO J = 1,I
!           E(I,J) = HAM(I,J)
!        ENDDO
!     ENDDO

!     ALLOCATE(WORK(33*REDUCED_NUM_EIG), IWORK(33*REDUCED_NUM_EIG))
     
!     CALL DSYEVR('V','I','L',REDUCED_NUM_EIG,E,REDUCED_NUM_EIG,0.0D0,0.0D0,1,1,DLAMCH('SAFE  MINIMUM'),NUMEVFOUND,
!    &              D,V,REDUCED_NUM_EIG,ISUPPZ,WORK,SIZE(WORK),IWORK,SIZE(IWORK),LOCALSTATUS)

!     DEALLOCATE(WORK,IWORK)

! CHECK FOR PROBLEMS WITH LAPACK ROUTINE

!     IF (LOCALSTATUS.NE.0 .OR. NUMEVFOUND.NE.1) THEN
!        PRINT *, 'ERROR IN DSYEVR ROUTINE: RETURN STATUS:', LOCALSTATUS, 'NUM EVS FOUND:', NUMEVFOUND
!        STOP
!     ENDIF

! ASSIGN ENERGY AND GROUND STATE WAVEFUNCTION

!     ENERGY = D(1)
!     GRSTWFU(:) = V(:,1)

!###################################################################
! ASSIGN THE GROUND STATE WAVEFUNCTION TO VARIABLE
!###################################################################
       DO I = 1,REDUCED_NUM_EIG
          GRSTWFU(I) = V(I,GROUND_STATE)
       ENDDO      

! DEALLOCATE MEMORY

     DEALLOCATE(XX, YY, ZZ, AIJ)
     DEALLOCATE(D, V, E)
     DEALLOCATE(HAM)

     RETURN
     END     

!###################################################################
!############### SUBROUTINE VH2OINTER #############################
!###################################################################

     SUBROUTINE VH2OINTER(VH2O_INTER)
     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

!################################################################### 
! CALCULATE THE WATER INTERMOLECULAR INTERACTION
! UNITS ARE KCAL/MOLE 
!
!################################################################### 

     DOUBLE PRECISION VH2O_INTER(REDUCED_NUM_EIG)
     DOUBLE PRECISION DX_OO,DY_OO,DZ_OO,ROO,ROO_6,ROO_12,LJ_TMP, &
                COULOMB_TMP,COULOMB_MONO,LJ_TERM,DX,DY,DZ, &
                RCOUL,QK,QJ,E_SQ,COULOMB,A,C
     INTEGER II,I,J,KAS,KK,START,COUNT,COUNT_SUB,JJ,KASKAS
     CHARACTER*1 DUMMY

!###################################################################
! BEGIN LOOP OVER THE EIGENSTATES
!
!################################################################### 
      DO 699 II = 1,REDUCED_NUM_EIG
!----REINITIALIZE THE VARIABLES:
        LJ_TMP = 0.0D0
        LJ_TERM = 0.0D0
        COULOMB = 0.0D0
        COULOMB_TMP = 0.0D0
        COULOMB_MONO = 0.0D0

!################################################################### 
! DETERMINE THE LJ COMPONENT OF THE WATER INTERMOLECULAR INTERACTION
! FOR THIS EIGENSTATE.
!
!################################################################### 
!     IF (QUICK) THEN
        DO I = 6,(NATOMS-4),3
           JJ = ATMPL(II,I) 
           DO J = I+3,(NATOMS-1),3 
              KASKAS = ATMPL(II,J)
              LJ_TMP = LJR(JJ,KASKAS)
              LJ_TERM = LJ_TMP + LJ_TERM
           ENDDO
        ENDDO

!     ELSE
!
!      DO I = 6,(NATOMS-4),3
!        DO J = I+3,(NATOMS-1),3
!           DX_OO = XX(II,I) - XX(II,J)
!           DY_OO = YY(II,I) - YY(II,J)
!           DZ_OO = ZZ(II,I) - ZZ(II,J)
!           ROO = DSQRT(DX_OO*DX_OO + DY_OO*DY_OO + DZ_OO*DZ_OO)
!            ROO = 3.1506D0/ROO
!           ROO_6 = ROO**6
!           ROO_12 = ROO_6*ROO_6
!              LJ_TMP = (ALJ*ROO_12) - (CLJ*ROO_6)
!            LJ_TMP = 4.0D0*(1.522D-01)*(ROO_12 - ROO_6)
!           LJ_TERM = LJ_TMP + LJ_TERM
!        ENDDO
!      ENDDO
!
!     ENDIF
!
!################################################################### 
! CALCULATE THE COULOMBIC COMPONENT OF THE WATER INTERMOLECULAR
! INTERACTION FOR THIS EIGENSTATE. 
!
!################################################################### 
     START = 5

       DO J = START,(NATOMS-3) 
          IF (J.EQ.(START+3)) THEN
               START = START + 3
           ENDIF 
          JJ = ATMPL(II,J)
          DO KAS = START+3,NATOMS  
               KASKAS = ATMPL(II,KAS)
             COULOMB_MONO = WATER_INTER_COULOMB(JJ,KASKAS)
             COULOMB_TMP = COULOMB_TMP + COULOMB_MONO            
          ENDDO
       ENDDO
  
!     DO J = START,(NATOMS-3)
!        IF (MOD(J,3).EQ.0) THEN
!           QJ = -0.834D0
!        ELSE
!          QJ = 0.417D0
!        ENDIF
!        IF (J.EQ.(START+3)) THEN
!           START = START + 3
!        ENDIF
!        DO KAS = START+3,NATOMS
!          DX = XX(II,J) - XX(II,KAS)
!           DY = YY(II,J) - YY(II,KAS)
!          DZ = ZZ(II,J) - ZZ(II,KAS)
!          RCOUL = DSQRT(DX*DX + DY*DY + DZ*DZ)
!          IF (MOD(KAS,3).EQ.0) THEN
!             QK = -0.834D0
!          ELSE
!             QK = 0.417D0
!          ENDIF
!              COULOMB_MONO = (QK*QJ*FOURPIEO)/(RCOUL)
!             COULOMB_TMP = COULOMB_TMP + COULOMB_MONO
!        ENDDO
!     ENDDO
!
!###################################################################
! FINALLY WORK OUT THE WATER INTERMOLECULAR INTERACTION FOR THIS
! EIGENSTATE
!
!###################################################################

     VH2O_INTER(II) = COULOMB_TMP + LJ_TERM

!###################################################################
! END OF EIGENSTATE

 699     ENDDO             

     RETURN 
     END 

!###################################################################
!################### VH2OINTRA SUBROUTINE ##########################
!###################################################################

     SUBROUTINE VH2OINTRA(VH2O_INTRA)
     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

!###################################################################
!
! DETERMINATION OF THE INTRAMOLECULAR POTENTIAL FOR H20
!
!################################################################### 

     DOUBLE PRECISION VH2O_INTRA(REDUCED_NUM_EIG)
     DOUBLE PRECISION DX_ROH1,DY_ROH1,DZ_ROH1,DX_ROH2,DY_ROH2,DZ_ROH2
     DOUBLE PRECISION ROH(3),ALPHA
     DOUBLE PRECISION VH2O_INTRA_TMP,VH2O_INTRA_R,VH2O_INTRA_ALPHA &
                ,VH2O_INTRA_RTMP,VH2O_INTRA_MONO
     INTEGER I,J,L,M,N
     CHARACTER*4 FORM
     DOUBLE PRECISION CHECKED_ACOS

!----INTIALIZE VARIABLES
     DO I = 1,REDUCED_NUM_EIG
        VH2O_INTRA(I) = 0.0D0
     ENDDO
        
!################################################################### 
! CALCULATE THE 2(O-H) BOND LENGTHS AND H-O-H BOND ANGLE IN EACH OF
! THE WATER MOLECULES OF THE EIGENSTATE CONSIDERED
!
!################################################################### 
     FORM = 'TIP3'
!----SUM OVER ALL EIGENSTATES
     DO 799 I = 1,REDUCED_NUM_EIG
        VH2O_INTRA_MONO = 0.0D00
!----SUM OVER ALL THE WATER MOLECULES
     DO N = 5,NATOMS-2,3
!----INTIALIZE VARIABLES
           ROH(1) = 0.0D00
           ROH(2) = 0.0D00
           ALPHA = 0.0D00
!----DETERMINE THE ROH DISTANCES IN WATER
        DX_ROH1 = PSIX(ATMPL(I,N)) - PSIX(ATMPL(I,N+1))
        DY_ROH1 = PSIY(ATMPL(I,N)) - PSIY(ATMPL(I,N+1))
           DZ_ROH1 = PSIZ(ATMPL(I,N)) - PSIZ(ATMPL(I,N+1))      
        ROH(1) = DSQRT(DX_ROH1*DX_ROH1 + DY_ROH1*DY_ROH1 + DZ_ROH1*DZ_ROH1)

        DX_ROH2 = PSIX(ATMPL(I,N+2)) - PSIX(ATMPL(I,N+1))
        DY_ROH2 = PSIY(ATMPL(I,N+2)) - PSIY(ATMPL(I,N+1))   
        DZ_ROH2 = PSIZ(ATMPL(I,N+2)) - PSIZ(ATMPL(I,N+1))   
           ROH(2) = DSQRT(DX_ROH2*DX_ROH2 + DY_ROH2*DY_ROH2 + DZ_ROH2*DZ_ROH2)

!----DETERMINE THE BENDING ANGLE, ALPHA
        ALPHA = (DX_ROH1*DX_ROH2 + DY_ROH1*DY_ROH2 + DZ_ROH1*DZ_ROH2)
        ALPHA = CHECKED_ACOS(ALPHA/(ROH(1)*ROH(2)))

!###################################################################
! NOW DETERMINE THE BENDING AND STRETCHING TERMS TOGETHER BY LOOPING 
! OVER THE 2 STRETCHED AND ONE ANGLE.
!
!###################################################################
        VH2O_INTRA_R = 0.0D00
        VH2O_INTRA_ALPHA = 0.0D00
        VH2O_INTRA_RTMP = 0.0D00

        DO J=1,2
           VH2O_INTRA_R = ROH(J) - H2OROHEQ
           VH2O_INTRA_R = 0.5D0*H2OKB*(VH2O_INTRA_R*VH2O_INTRA_R) 
              VH2O_INTRA_RTMP = VH2O_INTRA_RTMP + VH2O_INTRA_R
        ENDDO

        VH2O_INTRA_ALPHA = ALPHA-THETAEQ
        VH2O_INTRA_ALPHA = 0.5D0*H2OKTHETA*(VH2O_INTRA_ALPHA*VH2O_INTRA_ALPHA)
          
!################################################################### 
! DETERMINE THE CONTRIBUTION FROM THIS, SINGLE, WATER MOLECULE
!
!################################################################### 

        VH2O_INTRA_MONO = VH2O_INTRA_MONO + VH2O_INTRA_RTMP + VH2O_INTRA_ALPHA
 
!----CALC'ED THE ENERGY FOR 1 WATER IN 1 EIGENSTATE
!    FINISHED THE LOOP OVER A WATER MOLECULE
       ENDDO

!----FINALLY ASSIGN THE VH2O_INTRA FOR EIGENSTATE (I) |I>
        VH2O_INTRA(I) = VH2O_INTRA_MONO
 799     ENDDO

     RETURN
     END 

!#####################################################################
!################## VH3OINTRA SUBROUTINE #############################
!#####################################################################

     SUBROUTINE VH3OINTRA(VH3O_INTRA)
     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

!###################################################################
! DETERMINE THE CONTRIBUTION TO THE HAMILTONIAN FROM THE H3O+
! INTERMOLECULAR INTERACTION
!
! ONLY REQUIRES THE NUCLEAR COORDINATE VECTOR XX(I,J)
!
!###################################################################

     DOUBLE PRECISION VH3O_INTRA(REDUCED_NUM_EIG)
     DOUBLE PRECISION DX_ROH,DY_ROH,DZ_ROH,DX_HH,DY_HH,DZ_HH
     DOUBLE PRECISION ROH(3),ALPHA(3),HH(3)
     DOUBLE PRECISION VH3O_INTRA_TMP,VH3O_INTRA_R,VH3O_INTRA_ALPHA
     INTEGER I,J,KEG,L,M
     DOUBLE PRECISION CHECKED_ACOS

!###################################################################
! DETERMINE THE STRETCHING (3) AND BENDING (3) VECTORS
!
!###################################################################

     DO 899 I = 1,REDUCED_NUM_EIG
!----INTIALIZE VARIABLES
     VH3O_INTRA_TMP = 0.0D0
        DO KEG = 1,3
           ROH(KEG) = 0.0D0
           HH(KEG) = 0.0D0
           ALPHA(KEG) = 0.0D0
        ENDDO
!----DETERMINE THE ROH DISTANCES IN HYDRONIUM
        DX_ROH = XX(I,1) - XX(I,3)
        DY_ROH = YY(I,1) - YY(I,3)
        DZ_ROH = ZZ(I,1) - ZZ(I,3)
        ROH(1) = DSQRT(DX_ROH*DX_ROH + DY_ROH*DY_ROH + DZ_ROH*DZ_ROH)
        DX_ROH = XX(I,2) - XX(I,3)
           DY_ROH = YY(I,2) - YY(I,3)
        DZ_ROH = ZZ(I,2) - ZZ(I,3)
           ROH(2) = DSQRT(DX_ROH*DX_ROH + DY_ROH*DY_ROH + DZ_ROH*DZ_ROH)
        DX_ROH = XX(I,4) - XX(I,3)
           DY_ROH = YY(I,4) - YY(I,3)
           DZ_ROH = ZZ(I,4) - ZZ(I,3)
           ROH(3) = DSQRT(DX_ROH*DX_ROH + DY_ROH*DY_ROH + DZ_ROH*DZ_ROH)
!----DETERMINE THE BENDING ANGLES, ALPHA
!------(A)FIRST DETERMINE THE HH DISTANCES
        DX_HH = XX(I,1) - XX(I,2)
        DY_HH = YY(I,1) - YY(I,2)
           DZ_HH = ZZ(I,1) - ZZ(I,2)
        HH(1) = DSQRT(DX_HH*DX_HH + DY_HH*DY_HH + DZ_HH*DZ_HH)
        DX_HH = XX(I,2) - XX(I,4)
           DY_HH = YY(I,2) - YY(I,4)
           DZ_HH = ZZ(I,2) - ZZ(I,4)
           HH(2) = DSQRT(DX_HH*DX_HH + DY_HH*DY_HH + DZ_HH*DZ_HH)
        DX_HH = XX(I,4) - XX(I,1)
           DY_HH = YY(I,4) - YY(I,1)
           DZ_HH = ZZ(I,4) - ZZ(I,1)
           HH(3) = DSQRT(DX_HH*DX_HH + DY_HH*DY_HH + DZ_HH*DZ_HH)
!------(B)NOW DETERMINE THE ANGLES ALPHA
        ALPHA(1) = (-(HH(1)*HH(1)) + (ROH(1)*ROH(1)) + (ROH(2)*ROH(2)))
        ALPHA(1) = CHECKED_ACOS(ALPHA(1)/(2.0D0*ROH(1)*ROH(2)))

        ALPHA(2) = (-(HH(2)*HH(2)) + (ROH(2)*ROH(2)) + (ROH(3)*ROH(3)))
        ALPHA(2) =  CHECKED_ACOS(ALPHA(2)/(2.0D0*ROH(2)*ROH(3)))

        ALPHA(3) = (-(HH(3)*HH(3)) + (ROH(3)*ROH(3)) + (ROH(1)*ROH(1)))
        ALPHA(3) = CHECKED_ACOS(ALPHA(3)/(2.0D0*ROH(3)*ROH(1)))

!###################################################################
! NOW DETERMINE THE H3O+ INTERMOLECULAR ENERGY FOR THIS EIGENSTATE 
! FROM THE BENDING AND STRETCHING TERMS.
!
!###################################################################
        VH3O_INTRA_R = 0.0D0
        VH3O_INTRA_ALPHA = 0.0D0
         VH3O_INTRA_TMP = 0.0D0
      
          DO J=1,3
           VH3O_INTRA_R = (1.0D00 - DEXP(-1.0D0*AOHEQ*(ROH(J)-ROHEQ)))
           VH3O_INTRA_R = DOH * (VH3O_INTRA_R*VH3O_INTRA_R)
           VH3O_INTRA_ALPHA = ALPHA(J)-ALPHAEQ
           VH3O_INTRA_ALPHA = 0.5D0*KALPHA*(VH3O_INTRA_ALPHA* VH3O_INTRA_ALPHA)
           VH3O_INTRA_TMP = VH3O_INTRA_TMP + VH3O_INTRA_R + VH3O_INTRA_ALPHA
        ENDDO
!###################################################################
! FINALLY ASSIGN THE VH3O_INTRA FOR EIGENSTATE (I) |I>
!
!###################################################################
        VH3O_INTRA(I) = VH3O_INTRA_TMP

!###################################################################
! END OF LOOP OVER EIGENSTATES

 899     ENDDO

     RETURN
     END 

!###################################################################
!############### SUBROUTINE VINTER #################################
!###################################################################

     SUBROUTINE VINTER(V_INTER)
     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

!###################################################################
! CALCULATION OF THE INTERMOLECULAR INTERACTION BETWEEN H3O+ AND H2O
!
! COULOMB, REPULSION AND LENNARD-JONES TERMS
!
!###################################################################

     DOUBLE PRECISION V_INTER(REDUCED_NUM_EIG)
     DOUBLE PRECISION DX,DY,DZ,ROO,DXQQ,DYQQ,DZQQ,RQQ,QII,QJJ,ROO_6,ROO_12
     DOUBLE PRECISION LJ_TERM,REPULSE,COULOMB
     INTEGER I,II,JJ,KK,LL,III,JJJ
     DOUBLE PRECISION LJ_TMP,REPULSE_TMP,COULOMB_TMP,QTIP3P_HY,QTIP3P_OX
     CHARACTER*1 DUMMY

!###################################################################
! START OF LOOP OVER EIGENSTATES

     DO 599 LL = 1,REDUCED_NUM_EIG
        LJ_TERM = 0.0D00
!        LJ_TMP =0.0D0
        COULOMB = 0.0D0
!        COULOMB_TMP = 0.0D0
        REPULSE = 0.0D0
!        REPULSE_TMP = 0.0D0

!###################################################################
! DETERMINATION OF THE LJ AND REPLUSIVE TERMS
!
!###################################################################
     
        II = ATMPL(LL,3)     ! THE HYDRONIUM O ATOM IN THIS EIGENSTATE

!     IF (QUICK) THEN
        DO I = 6,(NATOMS-1),3
           JJ = ATMPL(LL,I)        ! THE CURRENT H2O ATOM
           REPULSE = REPULSE + REPULSE_INTER(II,JJ)
           LJ_TERM = LJ_TERM + LJ_INTER(II,JJ)
        ENDDO

!     ELSE
!
!     DO I = 6,(NATOMS-1),3 
!        DX = XX(LL,3) - XX(LL,I)
!         DY = YY(LL,3) - YY(LL,I)
!        DZ = ZZ(LL,3) - ZZ(LL,I)
!        ROO = DSQRT(DX*DX + DY*DY + DZ*DZ)
!        REPULSE_TMP = BIG_B*(1.0D0-DTANH(SMALL_B*(ROO - DOOEQ)))
!        ROO = SIGMA_MIX/ROO
!         ROO = ROO * ROO
!         ROO_6 = ROO**6 
!         ROO_12 = ROO**12
!        LJ_TMP = EPSILON_MIX*(ROO_12-ROO_6)
!        LJ_TERM = LJ_TERM + LJ_TMP
!        REPULSE = REPULSE + REPULSE_TMP
!     ENDDO
!
!     ENDIF

!###################################################################
! DETERMINATION OF THE COULOMBIC TERMS
!
!###################################################################

!      IF (QUICK) THEN
            DO II = 1,4
             III = ATMPL(LL,II)
             DO JJ = 5,NATOMS
                JJJ = ATMPL(LL,JJ)
                COULOMB = COULOMB + INTER_COULOMB(III,JJJ)
             ENDDO
          ENDDO
 
!      ELSE
!
!     DO II = 1,4
!        IF (II.EQ.3) THEN
!            QII = QINTERO
!        ELSE
!            QII = QINTERH
!        ENDIF
!        DO JJ = 5,NATOMS
!           DXQQ = XX(LL,II) - XX(LL,JJ) 
!           DYQQ = YY(LL,II) - YY(LL,JJ)
!           DZQQ = ZZ(LL,II) - ZZ(LL,JJ)
!           RQQ = DSQRT(DXQQ*DXQQ+DYQQ*DYQQ+DZQQ*DZQQ)
!           IF (MOD(JJ,3).EQ.0) THEN
!             QJJ = QTIP3P_O
!           ELSE
!             QJJ = QTIP3P_H
!           ENDIF
!           COULOMB_TMP = DAMPCHGE*(FOURPIEO*QII*QJJ)/(RQQ)
!           COULOMB = COULOMB + COULOMB_TMP
!        ENDDO
!     ENDDO
!
!      ENDIF
!
!###################################################################
! WORK OUT THE V_INTER() TERM FOR THIS EIGENSTATE FROM THE SUM OF 
! COULOMB, LJ AND  REPULSIVE TERMS

     V_INTER(LL) = LJ_TERM + COULOMB + REPULSE

!     PRINT *, 'EIGENSTATE:', LL
!     PRINT *, 'HYDRONIUM O:', ATMPL(LL,3)
!     PRINT *, 'LJ INTERACTION:', LJ_TERM
!     PRINT *, 'COULOMB ATTRACTION:', COULOMB
!     PRINT *, 'REPULSIVE TERM:', REPULSE

!     LJ_VALUES(LL)=LJ_TERM
!     COULOMB_VALUES(LL)=COULOMB
!     REP_VALUES(LL)=REPULSE

!###################################################################
! END OF LOOP OVER EIGENSTATES

 599     ENDDO

     RETURN
     END 

!#######################################################################
!################# SUBROUTINE BILDHAM ##################################
!####################################################################### 

     SUBROUTINE BILDHAM
     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

!###################################################################
! FORM THE NUMERICAL HAMILTONIAN FORMED BY THE NUCLEAR COORDINATES
!
!###################################################################

     INTEGER I,J
     DOUBLE PRECISION VH2O_INTER(REDUCED_NUM_EIG),VH3O_INTRA(REDUCED_NUM_EIG)
        DOUBLE PRECISION VH2O_INTRA(REDUCED_NUM_EIG),V_INTER(REDUCED_NUM_EIG)
      
!###################################################################
! FROM ONLY THE NUCLEAR COORDINATE VECTORS, XX(I,J) ETC, CALCULATE
! THE DIAGONAL AND THEN THE OFFDIAGONAL HAMILTONIAN ELEMENTS
!
! PASS BACK FROM THESE ROUTINES THE ARRAY OF DIAGONAL ELEMENT 
! COMPONENTS.
!
!###################################################################

     DO I = 1,REDUCED_NUM_EIG
           VH3O_INTRA(I) = 0.0D0
        VH2O_INTRA(I) = 0.0D0
        V_INTER(I) = 0.0D0
           VH2O_INTER(I) = 0.0D0
     ENDDO

     CALL VH3OINTRA(VH3O_INTRA)
     CALL VH2OINTRA(VH2O_INTRA)

! CALCULATE THE VARIOUS ENERGY COMPONENTS THAT ARE NEEDED

     CALL CALCATOMICINTERACTIONS
              
!----IS THERE MORE THAN ONE WATER MOLECULE?
        IF (NUM_EIG.GT.2) CALL VH2OINTER(VH2O_INTER)

     CALL VOFFDIAG2
     CALL VINTER(V_INTER)

!###################################################################
! BUILD THE HAMILTONIAN BY USING THE ELEMENT COMPONENTS JUST 
! DETERMINED
!
!###################################################################

     DO I = 1,REDUCED_NUM_EIG
         HAM(I,I) = VH3O_INTRA(I) + VH2O_INTRA(I) + V_INTER(I) + VH2O_INTER(I)
        IF (I.NE.REDUCED_NUM_EIG) THEN
           DO J = I+1,REDUCED_NUM_EIG
              HAM(I,J) = AIJ(I,J)
              HAM(J,I) = AIJ(I,J)
           ENDDO
        ENDIF

! DEBUGGING

!        H3O_INTRA_ENERGY(I)=VH3O_INTRA(I)
!        H2O_INTRA_ENERGY(I)=VH2O_INTRA(I)
!        H3O_H2O_ENERGY(I)=V_INTER(I)
!        H2O_H2O_ENERGY(I)=VH2O_INTER(I)

     ENDDO

!     DO I=1, REDUCED_NUM_EIG-1
!        DO J=I+1, REDUCED_NUM_EIG
!           OFF_DIAG_ENERGIES(I,J)=HAM(I,J)
!        ENDDO
!     ENDDO

     RETURN
     END

!#####################################################################

! SEPARATE OUT THE CALCULATION OF THE ENERGY TERMS BECAUSE WE ARE GOING REUSE THEM
! IN THE FORCE CALCULATIONS

     SUBROUTINE CALCATOMICINTERACTIONS ()

     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

! LOCAL VARIABLES

     INTEGER :: I, J
     DOUBLE PRECISION :: COULOMBPREFACTOR
     DOUBLE PRECISION :: QEXCH1, QWATER1, QH3O1, QEXCH2, QWATER2, QH3O2
     DOUBLE PRECISION :: R6,R12,R

     DOUBLE PRECISION, PARAMETER :: SIGMA6  = SIGMA_MIX**6
        DOUBLE PRECISION, PARAMETER :: TIP3P6  = H2OINTERSIGMA**6
     DOUBLE PRECISION, PARAMETER :: SIGMA12 = SIGMA6*SIGMA6
        DOUBLE PRECISION, PARAMETER :: TIP3P12 = TIP3P6*TIP3P6

! CALCULATE COULOMBIC INTERACTIONS

     ATOM_COULOMB = 0.0D0
     EACH_COULOMB = 0.0D0
     WATER_INTER_COULOMB = 0.0D0
     INTER_COULOMB = 0.0D0
     
     DO I = 1, NATOMS-1

        IF (MOD(I,3).EQ.0) THEN
           QEXCH1 = QEXCHO
           QWATER1 = QTIP3P_O
           QH3O1 = QINTERO
        ELSE
           QEXCH1 = QEXCHH
           QWATER1 = QTIP3P_H
           QH3O1 = QINTERH
        ENDIF

        DO J = I+1, NATOMS

           COULOMBPREFACTOR = FOURPIEO/INTERATOMICR(I,J)

           IF (MOD(J,3).EQ.0) THEN
           QEXCH2 = QEXCHO
           QWATER2 = QTIP3P_O
           QH3O2 = QINTERO
           ELSE
           QEXCH2 = QEXCHH
           QWATER2 = QTIP3P_H
           QH3O2 = QINTERH
           ENDIF

           EACH_COULOMB(I,J) = QEXCH1*QWATER2*COULOMBPREFACTOR
           ATOM_COULOMB(I) = ATOM_COULOMB(I) + EACH_COULOMB(I,J)
           EACH_COULOMB(J,I) = QEXCH2*QWATER1*COULOMBPREFACTOR
           ATOM_COULOMB(J) = ATOM_COULOMB(J) + EACH_COULOMB(J,I)

           WATER_INTER_COULOMB(I,J) = QWATER1*QWATER2*COULOMBPREFACTOR
           WATER_INTER_COULOMB(J,I) = WATER_INTER_COULOMB(I,J)

           INTER_COULOMB(I,J) = DAMPCHGE*QH3O1*QWATER2*COULOMBPREFACTOR
           INTER_COULOMB(J,I) = DAMPCHGE*QH3O2*QWATER1*COULOMBPREFACTOR      
        ENDDO
     ENDDO

! O-O TERMS

     DO I = 3, NATOMS-4, 3             ! LOOP OVER O ATOMS
        DO J = I+3, NATOMS-1, 3 
           R = INTERATOMICR(I,J) 
           R6 = 1.0D0/R**6
           R12 = R6*R6

! H2O-H2O O-O LJ TERM
           LJR(I,J) = 4.0D0*H2OINTEREPSILON*((TIP3P12*R12) - (TIP3P6*R6))
           LJR(J,I) = LJR(I,J)
           
! H2O-H3O LJ TERM
           LJ_INTER(I,J) = EPSILON_MIX*((SIGMA12*R12) - (SIGMA6*R6))
           LJ_INTER(J,I) = LJ_INTER(I,J)

! H2O-H3O REPULSION TERM
           REPULSE_INTER(I,J) = BIG_B*(1.0D0-DTANH(SMALL_B*(R - DOOEQ)))  
           REPULSE_INTER(J,I) = REPULSE_INTER(I,J)
        ENDDO
     ENDDO

     RETURN

     END SUBROUTINE

!#####################################################################
!################# SUBROUTINE JACOBI #################################
!#####################################################################
     
      SUBROUTINE JACOBI(A,N,D,V,NROT)
     IMPLICIT NONE

      INTEGER N,NROT,NMAX
      DOUBLE PRECISION A(N,N),D(N),V(N,N)
      PARAMETER (NMAX=500)
      INTEGER I,IP,IQ,J,RAC,MV
      DOUBLE PRECISION C,G,H,S,SM,T,TAU,THETA,TRESH,B(NMAX),Z(NMAX)

      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.0D0
11      CONTINUE
        V(IP,IP)=1.0D0
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.0D0
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.0D0
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.0D0)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2D0*SM/N**2
        ELSE
          TRESH=0.D0
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.D0*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP))).AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.D0
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5D0*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.D0)T=-T
              ENDIF
              C=1.D0/SQRT(1+T**2)
              S=T*C
              TAU=S/(1.D0+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.D0
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
23      CONTINUE
24    CONTINUE
      PRINT *, 'TOO MANY ITERATIONS IN JACOBI -- ABORTING '
     STOP

      RETURN
      END

!#####################################################################
!################## SUBROUTINE VOFFDIAG ##############################
!#####################################################################

     SUBROUTINE VOFFDIAG2
     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

!================================================================
!  R.CHRISTIE   10/11/98
!
! DETERMINATION OF THE OFF-DIAGONAL ELEMENTS OF THE HAMILTONIAN
!
!================================================================

     INTEGER I,J
     DOUBLE PRECISION ROH,ROO
     DOUBLE PRECISION ASYMM_COORD,F1,F2,F3

! GLOBAL STORAGE OF THE ZUNDEL TERMS FOR REUSE IN THE FORCE CALCULATIONS

     ZUNDEL_F = 0.0D0
     ZUNDEL_G = 0.0D0

     DO I = 1,REDUCED_NUM_EIG
        DO J = 1,REDUCED_NUM_EIG
           VIJEXCH(I,J) = 0.0D0
        ENDDO
     ENDDO

     IF (NUM_EIG.GT.2) CALL EXCHANGE2

!###################################################################
! CALCULATE SIMPLE MOLECULAR MECHANICS FUNCTIONS IN OFF DIAGONAL
! TERM. ONLY TWO VARIABLES: THE ASYMMETRIC STRETCH COORDINATE, AND
! (II) THE R(OO) VECTOR IN THE H5O2+ DIMER CONSIDERED AS 
! PARTICIPATING IN THE EXCHANGE CHARGE INTERACTION.

     DO I = 1,(REDUCED_NUM_EIG-1)
          DO J = (I+1),REDUCED_NUM_EIG

           IF (STATESINTERACT(I,J)) THEN

!           PRINT *, 'STATES', I, 'AND', J, 'INTERACT'

! ALREADY ASSIGNED A ZUNDEL SPECIES FOR THE INTERACTION OF THESE TWO 
! STATES 

           ROH = INTERATOMICR(ZUNDEL_SPECIES(I,J,4), ATMPL(I,3))
                                                                                       
!----DETERMINE THE R(OO) SEPARATION BETWEEN THE CENTRAL HYDRONIUM
!    IN STATE |I> AND THAT IN STATE |2>, AS IS NECESSARY IN THE
!    EXCHANGE TERM

           ROO = INTERATOMICR(ATMPL(I,3), ATMPL(J,3))

!##################################################################
! TWO VERSIONS OF DETERMINING THE ASYMMETRIC STRETCH COORDINATE
! ----------------------------------------------------------------
! 1. UDO SCHMITT PRE-PRINT VERSION OF CALCULATING THE ASYMMETRIC
!    STRETCH                                                    
!----DETERMINE THE UNIT VECTOR EOO
!             EOOX = (1.0D0/ROO)*OODX
!            EOOY = (1.0D0/ROO)*OODY
!            EOOZ = (1.0D0/ROO)*OODZ
!----DETERMINE THE ASYMMETRIC STRETCH TERMS QI
!             QX = ((OODX*0.5D0) - DX)*EOOX 
!             QY = ((OODY*0.5D0) - DY)*EOOY
!             QZ = ((OODZ*0.5D0) - DZ)*EOOZ
!----DETERMINE THE MAGNITUDE OF THE ASYMMETRIC STRETCH
!             ASYMM_COORD = (QX + QY + QZ)
!-----------------------------------------------------------------
! 2. PLAIN OLD METHOD FOR DETERMINING THE ASYMMETRIC STRETCH
!    COORDINATE                                              
           ASYMM_COORD = 0.5D0*ROO - ROH
!#################################################################          

           F1 = (1.00D0 + P*DEXP(-1.0D00*KEXP*(ROO-DOO)*(ROO-DOO)))
           F2 = 0.50D0*(1.0D0-DTANH(BETA*(ROO - BIG_ROOEQ )))
           F3 = 10.0D0*DEXP(-1.0D00*ALPHA_BLEH*(ROO - ROOEQ))
           ZUNDEL_F(I,J) = F1*(F2+F3)
           ZUNDEL_G(I,J) = DEXP(-1.0D00*GAMMA_MSEVB*ASYMM_COORD*ASYMM_COORD)

!           PRINT *, 'STATE 1:', I, 'STATE 2:', J
!           PRINT *, 'COMPONENTS OF THE OFF DIAGONAL CONTRIBUTION TO THE ENERGY:'
!           PRINT *, 'G:', G
!           PRINT *, 'F=F1(F2+F3)'
!           PRINT *, 'F1:', F1
!           PRINT *, 'F2:', F2
!           PRINT *, 'F3:', F3
!           PRINT *, 'F:', F
!           PRINT *, 'ROO:', ROO
!           PRINT *, 'ROH:', ROH
!           PRINT *, 'Q:', ASYMM_COORD
!           PRINT *, 'VIJEXCH:', VIJEXCH(I,J)

!----ASSIGN THE MATRIX ELEMENT
           AIJ(I,J) = (VIJ+VIJEXCH(I,J))*ZUNDEL_F(I,J)*ZUNDEL_G(I,J)
!----HERMETIAN MATRIX MIND
           AIJ(J,I) = AIJ(I,J)

           ELSE
           AIJ(I,J) = 0.0D0
           AIJ(J,I) = 0.0D0
           ENDIF
         ENDDO
      ENDDO

     RETURN 
     END

! ##################################################################################################

        SUBROUTINE EXCHANGE2

     USE MSEVB_COMMON

     IMPLICIT NONE
!---------------------------------------------------------------------
!  CALCULATION OF THE INTERACTION OF THE EXCHANGE CHARGE DISTRIBUTION 
!  FOR THE H5O2+ DIMERS RESULTING FROM <I|J> WITH THE REST OF THE 
!  TIP3P CHARGES ON THE REMAINING WATER MOLECULES
!
!---------------------------------------------------------------------

! AT THIS POINT ATOM_COULOMB(I) HOLDS THE SUM OF EXCHANGE INTERACTIONS WITH ALL OTHER ATOMS,
! WHERE I HOLDS THE EXCHANGE CHARGE AND THE OTHER ATOMS ARE ASSIGNED TIP3 CHARGES
! NEED TO THEREFORE SUBTRACT THE INTERACTIONS WITH THE OTHER MEMBERS OF THE HYDRONIUM SPECIES

     INTEGER I,J,K,L

     DO I = 1,(REDUCED_NUM_EIG-1)
        DO J = (I+1),REDUCED_NUM_EIG

           VIJEXCH(I,J) = 0.0D0

           IF (STATESINTERACT(I,J)) THEN

           DO K=1,7

! SUM OF EXCHANGE CHARGE INTERACTIONS WITH ALL OTHER ATOMS
              VIJEXCH(I,J) = VIJEXCH(I,J) + ATOM_COULOMB(ZUNDEL_SPECIES(I,J,K))

! SUBTRACT THE INTERACTIONS WITH THE OTHER ZUNDEL ATOMS
              DO L=1,7
                 IF (L.EQ.K) CYCLE
                 VIJEXCH(I,J) = VIJEXCH(I,J) - EACH_COULOMB(ZUNDEL_SPECIES(I,J,K),ZUNDEL_SPECIES(I,J,L))
              ENDDO
           ENDDO
           ENDIF
        ENDDO
     ENDDO

        RETURN
     END

!#############################################################################################

! J CHEM PHYS 2002 VERSION OF THIS ALGORITHM
! INCLUDING THE H-BOND ANGLE CUTOFF THAT THEY SAY DOESN`T MAKE ANY DIFFERENCE
! 04/12/2003 MODIFIED (CORRECTED) SO THAT EACH O ATOM CAN NOW FORM PART OF MORE THAN ONE
! HYDRONIUM SPECIES

! EXPERIMENTAL MIN ROH RATHER THAN ROH^2 CODE

        SUBROUTINE ASSIGNVB3

     USE COMMONS
     USE MSEVB_COMMON
     USE MSEVB_INTERFACES

     IMPLICIT NONE

     DOUBLE PRECISION MINHBONDANGLERAD,HBONDANGLE

     DOUBLE PRECISION MINSUMROH,CURRENTSUMROH     
     LOGICAL :: ASSIGNEDPIVOTSTATE

     INTEGER :: GENERATORSTATE(MAXNUMVBSTATES)  ! VB STATE I WAS GENERATED FROM GENERATORSTATE(I)

! HOLDS THE CLOSEST H ATOMS FOR THE OS AND THE CLOSEST O ATOMS FOR THE HS

     INTEGER PROXIMALATOMS(NATOMS,NUM_HYD)   
     INTEGER NEIGHBOURCOUNTS(NATOMS)

     INTEGER CURRENTO, CURRENTH

     INTEGER FIRSTVBSTATE, SECONDVBSTATE, CURRENTSHELLHCOUNT, NEXTSHELLHCOUNT, CURRENTSHELL
     INTEGER CURRENTEXCHANGEH, NEWEIGENO, OLDVBSTATE
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: EXCHANGEHS
     INTEGER, ALLOCATABLE :: FINGERPRINTS(:,:)
     LOGICAL NEWFINGERPRINT
     INTEGER O1POS,O2POS,CENTRALHPOS
     INTEGER FIRSTSTATEO, SECONDSTATEO

! LOOP VARIABLES

     INTEGER I,J,K,L,EIGEN_O

! CONVERT MINIMUM H BOND ANGLE TO RADIANS

     MINHBONDANGLERAD = (MINHBONDANGLE/180.0D0)*PI

! ORDER THE O ATOMS IN TERMS OF THE DISTANCE FROM EACH H, ONLY CONSIDER THOSE WHICH ARE LESS THAN
! THE CUTOFF DISTANCE APART

     NEIGHBOURCOUNTS = 0

     DO I = 1,NUM_HYD

        CURRENTH = I + INT((I-1)/2)

        DO J = 1,NUM_EIG
           CURRENTO = 3*J
           IF (INTERATOMICR(CURRENTO,CURRENTH).LT.MAXHBONDLENGTH) THEN
           NEIGHBOURCOUNTS(CURRENTH) = NEIGHBOURCOUNTS(CURRENTH) + 1

           DO K=1,NEIGHBOURCOUNTS(CURRENTH)
              IF (K.EQ.NEIGHBOURCOUNTS(CURRENTH)) THEN
                 PROXIMALATOMS(CURRENTH,K) = CURRENTO
                 EXIT
              ELSEIF (INTERATOMICR(CURRENTO,CURRENTH).LT.INTERATOMICR(PROXIMALATOMS(CURRENTH,K),CURRENTH)) THEN
                 DO L=NEIGHBOURCOUNTS(CURRENTH),K+1,-1
                 PROXIMALATOMS(CURRENTH,L)=PROXIMALATOMS(CURRENTH,L-1)
                 ENDDO
                 PROXIMALATOMS(CURRENTH,K) = CURRENTO
                 EXIT
              ENDIF
           ENDDO
           ENDIF
        ENDDO

! CHECK WHETHER WE HAVE ASSIGNED AT LEAST ONE O 

        IF (NEIGHBOURCOUNTS(CURRENTH).EQ.0) THEN
           PRINT *, 'WITHIN THE CUTOFF OF', MAXHBONDLENGTH, 'PROGRAM ASSIGNED ZERO O ATOMS TO H', CURRENTH
           PRINT *, 'EACH H MUST BE ASSIGNED AT LEAST ONE O'
           PRINT *, 'INCREASE CUTOFF OR EXAMINE YOUR GEOMETRY'

           DO J=1, NATOMS
           IF (MOD(J,3).EQ.0) THEN
              PRINT *, 'O', PSIX(J), PSIY(J), PSIZ(J)
           ELSE
              PRINT *, 'H', PSIX(J), PSIY(J), PSIZ(J)
           ENDIF
           ENDDO
           STOP
        ENDIF

     ENDDO

!     PRINT *, 'ASSIGNED OS:'

!     DO I=1, NUM_HYD
        
!        PRINT *, 'H:', HPOSITION(I)
!        DO J=1, ASSIGNEDOS(I)
!           PRINT *, J, 3*NEIGHBOURINGOS(I,J)
!        ENDDO
!     ENDDO

! TRY A SIMPLE VB ASSIGNMENT FIRST, IF THAT DOESN`T WORK TRY THE ITERATIVE PROCEDURE

     ASSIGNEDPIVOTSTATE = .FALSE.
     REDUCED_NUM_EIG = 0

     CALL TRYSIMPLEVBASSIGNMENT(PROXIMALATOMS, NEIGHBOURCOUNTS, REDUCED_NUM_EIG, ATMPL(1:2,:))

     IF (REDUCED_NUM_EIG.EQ.0) THEN

        CALL ASSIGNPIVOTVBSTATE(PROXIMALATOMS, NEIGHBOURCOUNTS, ASSIGNEDPIVOTSTATE, ATMPL(1,:))

        IF (ASSIGNEDPIVOTSTATE) THEN
           REDUCED_NUM_EIG = 1
        ELSE
           PRINT *, 'UNABLE TO ASSIGN PIVOT VB STATE'
           PRINT *, 'SYSTEM:'
           DO I = 1, NATOMS
           IF (MOD(I,3).EQ.0) THEN
              PRINT *, 'O', PSIX(I), PSIY(I), PSIZ(I)
           ELSE
              PRINT *, 'H', PSIX(I), PSIY(I), PSIZ(I)
           ENDIF
           ENDDO
           STOP   
        ENDIF
     ENDIF

     GENERATORSTATE(1:REDUCED_NUM_EIG) = -1    ! -1 INDICATES THAT THE STATES ARE PIVOT STATES
     
!     PRINT *, 'NUMBER OF PIVOT STATES ASSIGNED:', REDUCED_NUM_EIG
     
!     DO I=1,REDUCED_NUM_EIG
!        PRINT *, '*** STATE', I, '***'
!        DO J=1,NATOMS
!           PRINT *, J, ATMPL(I,J)
!        ENDDO
!     ENDDO

!     STOP

     ALLOCATE(FINGERPRINTS(MAXNUMVBSTATES,NATOMS))
     FINGERPRINTS = 0

     DO I=1,REDUCED_NUM_EIG

        FINGERPRINTS(I,ATMPL(I,1)) = ATMPL(I,3)
        FINGERPRINTS(I,ATMPL(I,2)) = ATMPL(I,3)
        FINGERPRINTS(I,ATMPL(I,4)) = ATMPL(I,3)
        
        DO J=6, (3*NUM_EIG), 3
           FINGERPRINTS(I,ATMPL(I,J-1))=ATMPL(I,J)
           FINGERPRINTS(I,ATMPL(I,J+1))=ATMPL(I,J)
        ENDDO
     ENDDO

! INITIALISE THE VB STATE INTERACTION MATRIX

     STATESINTERACT = .FALSE.

! ALLOCATE THE EXCHANGE HS ARRAY

     IF (SHELLSTOCOUNT.LT.2) THEN
        ALLOCATE(EXCHANGEHS((REDUCED_NUM_EIG+2),2))
     ELSE
        ALLOCATE(EXCHANGEHS((REDUCED_NUM_EIG+2)*3*(2**(SHELLSTOCOUNT-2)), 2))
     ENDIF

! RUN THROUGH THE HYDRATION SHELLS

     IF (REDUCED_NUM_EIG.GT.1) THEN

        CURRENTSHELLHCOUNT = 0

        FIRSTPIVOTSTATE: DO I = 1, 4

           IF (I.NE.3) THEN

! FOR SECOND AND SUBSEQUENT PIVOT STATES THE FOURTH ENTRY IS THE SYMMETRIC H
           
           DO SECONDVBSTATE = 2, REDUCED_NUM_EIG

              IF (ATMPL(SECONDVBSTATE,4).EQ.ATMPL(1,I)) THEN

                 CALL ASSIGNZUNDELSPECIES(1,SECONDVBSTATE,ATMPL(1,I),.FALSE.,MINHBONDANGLERAD)

                 EXCHANGEHS(CURRENTSHELLHCOUNT+1,1) = ATMPL(SECONDVBSTATE,1)
                 EXCHANGEHS(CURRENTSHELLHCOUNT+1,2) = SECONDVBSTATE
                 EXCHANGEHS(CURRENTSHELLHCOUNT+2,1) = ATMPL(SECONDVBSTATE,2)
                 EXCHANGEHS(CURRENTSHELLHCOUNT+2,2) = SECONDVBSTATE

                 CURRENTSHELLHCOUNT = CURRENTSHELLHCOUNT+2
                 CYCLE FIRSTPIVOTSTATE
              ENDIF
           ENDDO

           CURRENTSHELLHCOUNT = CURRENTSHELLHCOUNT+1
           EXCHANGEHS(CURRENTSHELLHCOUNT,1) = ATMPL(1,I)
           EXCHANGEHS(CURRENTSHELLHCOUNT,2) = 1           
           ENDIF
        ENDDO FIRSTPIVOTSTATE

     ELSE
        CURRENTSHELLHCOUNT = 3
        EXCHANGEHS(1,1) = ATMPL(REDUCED_NUM_EIG,1)
        EXCHANGEHS(1,2) = REDUCED_NUM_EIG
        EXCHANGEHS(2,1) = ATMPL(REDUCED_NUM_EIG,2)
        EXCHANGEHS(2,2) = REDUCED_NUM_EIG
        EXCHANGEHS(3,1) = ATMPL(REDUCED_NUM_EIG,4)
        EXCHANGEHS(3,2) = REDUCED_NUM_EIG
     ENDIF
     
! LOOP THROUGH HYDRATION SHELLS, 3 AT PRESENT (SET IN MSEVB_COMMON MODULE)
! ONLY ASSIGN NEW VB STATES FOR THE SPECIFIED NUMBER OF SHELLS

     DO CURRENTSHELL = 1, SHELLSTOCOUNT

!          PRINT *, 'SHELL', CURRENTSHELL, CURRENTSHELLHCOUNT

!          PRINT *, 'EXCHANGING HS IN THIS SHELL'
!          DO I = 1, CURRENTSHELLHCOUNT
!             PRINT *, I, EXCHANGEHS(I,1)
!          ENDDO

        NEXTSHELLHCOUNT = 0

        CURRENTSHELLLOOP: DO I = 1, CURRENTSHELLHCOUNT
           
           CURRENTEXCHANGEH = EXCHANGEHS(I,1)

           IF (NEIGHBOURCOUNTS(CURRENTEXCHANGEH).GT.1) THEN

!             PRINT *, 'CURRENTLY CONSIDERING EXCHANGE H', CURRENTEXCHANGEH
!             PRINT *, 'ASSIGNED O COUNT:', NEIGHBOURCOUNTS(CURRENTEXCHANGEH)

! AT MOST, EACH H CAN BE ASSIGNED TO TWO DIFFERENT HYDRONIUM SPECIES THEREFORE CHECK THE CLOSEST TWO OS
! IF THE H WAS ONLY ASSIGNED ONE O THEN THIS HYDRONIUM HAS ALREADY BEEN TAKEN CARE OF
        
           FIRSTVBSTATE = EXCHANGEHS(I,2)

           IF (ATMPL(FIRSTVBSTATE,3).EQ.PROXIMALATOMS(CURRENTEXCHANGEH,1)) THEN
              NEWEIGENO = PROXIMALATOMS(CURRENTEXCHANGEH,2)
           ELSE
              NEWEIGENO = PROXIMALATOMS(CURRENTEXCHANGEH,1)
           ENDIF

! CHECK THE H BOND ANGLE

           O2POS = ATMPL(FIRSTVBSTATE,3)

           HBONDANGLE = (PSIX(NEWEIGENO)-PSIX(CURRENTEXCHANGEH))*(PSIX(O2POS)-PSIX(CURRENTEXCHANGEH))+& 
                             (PSIY(NEWEIGENO)-PSIY(CURRENTEXCHANGEH))*(PSIY(O2POS)-PSIY(CURRENTEXCHANGEH))+&
                          (PSIZ(NEWEIGENO)-PSIZ(CURRENTEXCHANGEH))*(PSIZ(O2POS)-PSIZ(CURRENTEXCHANGEH))


           HBONDANGLE = DACOS(HBONDANGLE/(INTERATOMICR(NEWEIGENO,CURRENTEXCHANGEH) &
                              * INTERATOMICR(O2POS,CURRENTEXCHANGEH)))
              
!              PRINT *, 'BOND ANGLE: ', O1POS, CENTRALHPOS, O2POS, ((HBONDANGLE/PI)*180D0)
              
           IF (HBONDANGLE.GE.MINHBONDANGLERAD) THEN

              REDUCED_NUM_EIG = REDUCED_NUM_EIG + 1

! CHECK WHETHER WE HAVE SEEN THIS VB STATE BEFORE - CAN HAPPEN WITH CYCLES

              FINGERPRINTS(REDUCED_NUM_EIG,:)=FINGERPRINTS(FIRSTVBSTATE,:)
              FINGERPRINTS(REDUCED_NUM_EIG,CURRENTEXCHANGEH)=NEWEIGENO

              DO J=1,REDUCED_NUM_EIG-1
                 IF (J.NE.FIRSTVBSTATE) THEN

                 NEWFINGERPRINT = .FALSE.

                 DO K=1, NATOMS
                    IF (FINGERPRINTS(J,K).NE.FINGERPRINTS(REDUCED_NUM_EIG,K)) THEN
                    NEWFINGERPRINT = .TRUE.
                    EXIT
                    ENDIF
                 ENDDO

                 IF (.NOT.NEWFINGERPRINT) THEN
!                    PRINT *, 'SEEN THIS FINGERPRINT BEFORE', REDUCED_NUM_EIG, J
!                    DO K=1, NATOMS
!                    PRINT *, K, ATMPL(J,K)
!                    ENDDO
                    REDUCED_NUM_EIG = REDUCED_NUM_EIG - 1
                    CYCLE CURRENTSHELLLOOP                   
                 ENDIF
                 ENDIF
              ENDDO
              
              ATMPL(REDUCED_NUM_EIG, 3) = NEWEIGENO
              ATMPL(REDUCED_NUM_EIG, 4) = CURRENTEXCHANGEH                 

!                   PRINT *, 'NEW EIGEN O:', NEWEIGENO

              DO J=6, NATOMS, 3
                 IF (ATMPL(FIRSTVBSTATE,J).EQ.NEWEIGENO) THEN
                 ATMPL(REDUCED_NUM_EIG,1) = ATMPL(FIRSTVBSTATE,J-1)
                 ATMPL(REDUCED_NUM_EIG,2) = ATMPL(FIRSTVBSTATE,J+1)

! ADD THESE PROTONS TO THE POTENTIAL EXCHANGE LIST FOR THE NEXT SHELL

                 IF (CURRENTSHELL.NE.SHELLSTOCOUNT) THEN

                    NEXTSHELLHCOUNT = NEXTSHELLHCOUNT + 1
                    IF ((CURRENTSHELLHCOUNT+NEXTSHELLHCOUNT).GT.SIZE(EXCHANGEHS,1)) THEN
                    PRINT *, 'TOO MANY HS IN EXCHANGE H LIST'
                    PRINT *, 'INCREASE SIZE OF ARRAY'
                    DO K = 1, NATOMS
                       IF (MOD(K,3).EQ.0) THEN
                          PRINT *, 'O', PSIX(K), PSIY(K), PSIZ(K)
                       ELSE
                          PRINT *, 'H', PSIX(K), PSIY(K), PSIZ(K)
                       ENDIF
                    ENDDO
                    PRINT *, 'SIZE OF ARRAY:', SIZE(EXCHANGEHS,1)
                    PRINT *, 'CURRENT SHELL:', CURRENTSHELL
                    DO K=1,CURRENTSHELLHCOUNT
                       PRINT *, K, EXCHANGEHS(K,1), EXCHANGEHS(K,2)
                    ENDDO
                    PRINT *, 'NEXT SHELL:', (CURRENTSHELL+1)
                    DO K=CURRENTSHELLHCOUNT+1, CURRENTSHELLHCOUNT+NEXTSHELLHCOUNT-1
                       PRINT *, K, EXCHANGEHS(K,1), EXCHANGEHS(K,2)
                    ENDDO
                    STOP
                    ENDIF
                    EXCHANGEHS(CURRENTSHELLHCOUNT+NEXTSHELLHCOUNT,1) = ATMPL(REDUCED_NUM_EIG,1)
                    EXCHANGEHS(CURRENTSHELLHCOUNT+NEXTSHELLHCOUNT,2) = REDUCED_NUM_EIG
                 
                    NEXTSHELLHCOUNT = NEXTSHELLHCOUNT + 1
                    IF ((CURRENTSHELLHCOUNT+NEXTSHELLHCOUNT).GT.SIZE(EXCHANGEHS,1)) THEN
                    PRINT *, 'TOO MANY HS IN EXCHANGE H LIST'
                    PRINT *, 'INCREASE SIZE OF ARRAY'
                    DO K = 1, NATOMS
                       IF (MOD(K,3).EQ.0) THEN
                          PRINT *, 'O', PSIX(K), PSIY(K), PSIZ(K)
                       ELSE
                          PRINT *, 'H', PSIX(K), PSIY(K), PSIZ(K)
                       ENDIF
                    ENDDO
                    PRINT *, 'SIZE OF ARRAY:', SIZE(EXCHANGEHS,1)
                    PRINT *, 'CURRENT SHELL:', CURRENTSHELL
                    DO K=1,CURRENTSHELLHCOUNT
                       PRINT *, K, EXCHANGEHS(K,1), EXCHANGEHS(K,2)
                    ENDDO
                    PRINT *, 'NEXT SHELL:', (CURRENTSHELL+1)
                    DO K=CURRENTSHELLHCOUNT+1, CURRENTSHELLHCOUNT+NEXTSHELLHCOUNT-1
                       PRINT *, K, EXCHANGEHS(K,1), EXCHANGEHS(K,2)
                    ENDDO
                    STOP
                    ENDIF
                    EXCHANGEHS(CURRENTSHELLHCOUNT+NEXTSHELLHCOUNT,1) = ATMPL(REDUCED_NUM_EIG,2)
                    EXCHANGEHS(CURRENTSHELLHCOUNT+NEXTSHELLHCOUNT,2) = REDUCED_NUM_EIG
                 ENDIF
                 EXIT
                 ENDIF
              ENDDO

! ADD THE OLD HYDRONIUM SPECIES, MINUS THE EXCHANGE PROTON

              ATMPL(REDUCED_NUM_EIG,6) = ATMPL(FIRSTVBSTATE, 3)
           
              IF (ATMPL(FIRSTVBSTATE,1).EQ.CURRENTEXCHANGEH) THEN
                 ATMPL(REDUCED_NUM_EIG,5) = ATMPL(FIRSTVBSTATE,2)
                 ATMPL(REDUCED_NUM_EIG,7) = ATMPL(FIRSTVBSTATE,4)
              ELSE
                 ATMPL(REDUCED_NUM_EIG,5) = ATMPL(FIRSTVBSTATE,1)            
                 
                 IF (ATMPL(FIRSTVBSTATE,2).EQ.CURRENTEXCHANGEH) THEN
                 ATMPL(REDUCED_NUM_EIG,7) = ATMPL(FIRSTVBSTATE,4)                    
                 ELSE
                 ATMPL(REDUCED_NUM_EIG,7) = ATMPL(FIRSTVBSTATE,2)     
                 ENDIF
              ENDIF
              
! COPY ACROSS THE REMAINDER OF THE WATER SPECIES FROM THE OLD VB STATE, THESE WILL NOT HAVE CHANGED

              K = 9

              DO J = 6, NATOMS-1, 3
                 IF (ATMPL(FIRSTVBSTATE, J).EQ.ATMPL(REDUCED_NUM_EIG,3)) CYCLE
              
                 ATMPL(REDUCED_NUM_EIG, K-1) = ATMPL(FIRSTVBSTATE,J-1)
                 ATMPL(REDUCED_NUM_EIG, K) = ATMPL(FIRSTVBSTATE,J)
                 ATMPL(REDUCED_NUM_EIG, K+1) = ATMPL(FIRSTVBSTATE,J+1)
              
                 K = K + 3
              ENDDO

!              IF (.NOT.NEWFINGERPRINT) THEN
!                 PRINT *, 'NEW VB STATE:'
!                 DO K= 1, NATOMS
!                 PRINT *, K, ATMPL(REDUCED_NUM_EIG,K)
!                 ENDDO
!              ENDIF

! ASSIGN THE ZUNDEL SPECIES

              CALL ASSIGNZUNDELSPECIES (FIRSTVBSTATE,REDUCED_NUM_EIG,CURRENTEXCHANGEH,.TRUE.,MINHBONDANGLERAD)

 !                 PRINT *, 'NEW ZUNDEL SPECIES:'
 !                 DO K = 1, 7
 !                 PRINT *, K, ZUNDEL_SPECIES(FIRSTVBSTATE, REDUCED_NUM_EIG, K)
 !                 ENDDO

              GENERATORSTATE(REDUCED_NUM_EIG) = FIRSTVBSTATE

! CHECK WHETHER ANY OTHER VB STATES USE THIS FIRST VB STATE O AS THE HYDRONIUM WITH CURRENTEXCHANGEH AS ONE OF THE HS
! ASSIGN INTERACTIONS IF THEY DO

              DO J = 1, (REDUCED_NUM_EIG-1)
                 IF (J.NE.FIRSTVBSTATE .AND. ATMPL(J,3).EQ.O2POS .AND. &
                            (CURRENTEXCHANGEH.EQ.ATMPL(J,1) & 
                            .OR. CURRENTEXCHANGEH.EQ.ATMPL(J,2) .OR. CURRENTEXCHANGEH.EQ.ATMPL(J,4))) THEN

                 CALL ASSIGNZUNDELSPECIES(J,REDUCED_NUM_EIG,CURRENTEXCHANGEH,.TRUE.,MINHBONDANGLERAD)

!                     PRINT *, 'STATE', J, 'ALSO INTERACTS WITH THIS NEW STATE'
!                     PRINT *, 'NEW ZUNDEL SPECIES:'
!                     DO K = 1, 7
!                     PRINT *, K, ZUNDEL_SPECIES(J, REDUCED_NUM_EIG, K)
!                     ENDDO               
                 ENDIF
              ENDDO
              
           ENDIF
 
           ENDIF

        ENDDO CURRENTSHELLLOOP ! END OF LOOP OVER CURRENT SHELL

! UPDATE THE EXCHANGE H LIST

        DO I = 1, NEXTSHELLHCOUNT
           EXCHANGEHS(I,1) = EXCHANGEHS(I+CURRENTSHELLHCOUNT,1)
           EXCHANGEHS(I,2) = EXCHANGEHS(I+CURRENTSHELLHCOUNT,2)      
        ENDDO
       
        CURRENTSHELLHCOUNT = NEXTSHELLHCOUNT

     ENDDO

! DEALLOCATE MEMORY

     DEALLOCATE(EXCHANGEHS,FINGERPRINTS)

! CHECK FOR UNASSIGNED INTERACTIONS - THIS BIT WORKS IN TANDEM WITH THE IN LOOP CHECKING
! THE IN LOOP CHECKING TAKES CARE OF EXTRA INTERACTIONS WITH STATES WHICH WERE ASSIGNED BEFORE THE CURRENT ONE
! THIS BIT TAKES CARE OF EXTRA INTERACTIONS WITH STATES WHICH WERE ASSIGNED AFTER THE CURRENT ONE
! SEE JOURNAL ENTRY 06/04/04

     DO I = 2, REDUCED_NUM_EIG-1
        IF (GENERATORSTATE(I).NE.-1) THEN

           DO J = I+1, REDUCED_NUM_EIG
           IF (J.NE.GENERATORSTATE(I) .AND. (.NOT.STATESINTERACT(I,J)) .AND. &
                     ATMPL(J,3).EQ.ATMPL(GENERATORSTATE(I),3)) THEN

!              PRINT *, 'MISSED INTERACTION BETWEEN STATES', I, 'AND', J

! ATMPL(I,4) HOLDS THE EXCHANGE H WITH WHICH THIS STATE WAS INITIALLY GENERATED (EXCEPT FOR THE PIVOT STATE(S))
              CALL ASSIGNZUNDELSPECIES (I,J,ATMPL(I,4),.TRUE.,MINHBONDANGLERAD)

!              PRINT *, 'NEW ZUNDEL SPECIES'
!              DO K=1,7
!                 PRINT *, K, ZUNDEL_SPECIES(I,J,K)
!              ENDDO

           ENDIF
           ENDDO
        ENDIF
     ENDDO

     IF (REDUCED_NUM_EIG.GT.MAXNUMVBSTATES) THEN
        PRINT *, 'NUMBER OF VB STATES GENERATED EXCEEDS HARD LIMIT'
        STOP
     ENDIF

!      PRINT *, 'VB STATES ASSIGNED', REDUCED_NUM_EIG

!       DO I=1, REDUCED_NUM_EIG

!          PRINT *, '*** STATE', I, '***'

!          DO J=1, NATOMS
!             PRINT *, J, ATMPL(I,J)
!          ENDDO
!       ENDDO

!     PRINT *, 'EXCHANGING H LIST'

!     DO I=1, (CURRENTSHELLSTARTH + CURRENTSHELLHCOUNT - 1)
!        PRINT *, I, HPOSITION(EXCHANGEHS(I))
!     ENDDO

!       PRINT *, 'ASSIGNED ZUNDEL SPECIES'

!       DO I = 1, (REDUCED_NUM_EIG-1)
!          DO J = (I+1), REDUCED_NUM_EIG
!             IF (STATESINTERACT(I,J)) THEN
!             PRINT *, 'STATE 1:', I, 'STATE 2:', J
!             DO K=1,7
!                PRINT *, K, ZUNDEL_SPECIES(I,J,K)
!             ENDDO
!             ENDIF
!          ENDDO
!       ENDDO

!     STOP

     END SUBROUTINE

! ######################################################################################

! ASSIGNS THE EXCHANGE SPECIES BETWEEN TWO HYDRONIUM SPECIES, CHECKING THE H BOND ANGLE
! IF NECESSARY

     SUBROUTINE ASSIGNZUNDELSPECIES(FIRSTVBSTATE, SECONDVBSTATE, EXCHANGEH, CHECKEDHBOND, MINHBONDANGLERAD)

     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: FIRSTVBSTATE, SECONDVBSTATE, EXCHANGEH
     LOGICAL, INTENT(IN) :: CHECKEDHBOND
     DOUBLE PRECISION, INTENT(IN) :: MINHBONDANGLERAD
     
     INTEGER :: O1POS, O2POS, CENTRALHPOS
     DOUBLE PRECISION :: HBONDANGLE
     INTEGER :: LOWERINDEX, HIGHERINDEX

! MAKE SURE THAT THE FIRST VB STATE HAS THE LOWER INDEX - MATRICES ARE UPPER-TRIANGULAR

     IF (FIRSTVBSTATE.LT.SECONDVBSTATE) THEN
        LOWERINDEX = FIRSTVBSTATE
        HIGHERINDEX = SECONDVBSTATE
     ELSE
        LOWERINDEX = SECONDVBSTATE
        HIGHERINDEX = FIRSTVBSTATE
     ENDIF

     IF (.NOT.STATESINTERACT(LOWERINDEX,HIGHERINDEX)) THEN

! CHECK THE H BOND ANGLE IF NECESSARY

        IF (.NOT.CHECKEDHBOND) THEN
           O1POS = ATMPL(FIRSTVBSTATE,3)
           O2POS = ATMPL(SECONDVBSTATE,3)

           HBONDANGLE = (PSIX(O1POS)-PSIX(EXCHANGEH))*(PSIX(O2POS)-PSIX(EXCHANGEH)) + &
                          (PSIY(O1POS)-PSIY(EXCHANGEH))*(PSIY(O2POS)-PSIY(EXCHANGEH)) + &
                       (PSIZ(O1POS)-PSIZ(EXCHANGEH))*(PSIZ(O2POS)-PSIZ(EXCHANGEH))

           HBONDANGLE = DACOS(HBONDANGLE/(INTERATOMICR(O1POS,EXCHANGEH)*INTERATOMICR(O2POS,EXCHANGEH))) 
              
!              PRINT *, 'BOND ANGLE: ', O1POS, EXCHANGEH, O2POS, ((HBONDANGLE/PI)*180D0)
              
           IF (HBONDANGLE.LT.MINHBONDANGLERAD) RETURN
        ENDIF

        ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,3) = ATMPL(LOWERINDEX,3)  ! O
        ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,4) = EXCHANGEH            ! EXCHANGING H
        ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,6) = ATMPL(HIGHERINDEX,3) ! O
           
! SORT OUT THE OTHER PROTONS FROM THE FIRST HYDRONIUM SPECIES

        IF (ATMPL(LOWERINDEX,1).EQ.EXCHANGEH) THEN
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,1) = ATMPL(LOWERINDEX,2)         
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,2) = ATMPL(LOWERINDEX,4)
        ELSE
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,1) = ATMPL(LOWERINDEX,1)
              
           IF (ATMPL(LOWERINDEX,2).EQ.EXCHANGEH) THEN
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,2) = ATMPL(LOWERINDEX,4)
           ELSE
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,2) = ATMPL(LOWERINDEX,2)
           ENDIF
        ENDIF

! AND THOSE FROM THE SECOND HYDRONIUM SPECIES
              
        IF (ATMPL(HIGHERINDEX,1).EQ.EXCHANGEH) THEN
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,5) = ATMPL(HIGHERINDEX,2)    
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,7) = ATMPL(HIGHERINDEX,4)
        ELSE
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,5) = ATMPL(HIGHERINDEX,1)
           
           IF (ATMPL(HIGHERINDEX,2).EQ.EXCHANGEH) THEN
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,7) = ATMPL(HIGHERINDEX,4)
           ELSE
           ZUNDEL_SPECIES(LOWERINDEX,HIGHERINDEX,7) = ATMPL(HIGHERINDEX,2)
           ENDIF
        ENDIF

! SET THE INTERACTION LABEL

        STATESINTERACT(LOWERINDEX,HIGHERINDEX) = .TRUE.

     ENDIF

     RETURN

     END SUBROUTINE

! ###############################################################################################

! TRIES THE SIMPLEST POSSIBLE PIVOT VB ASSIGNMENT BY ASSIGNING EACH H TO THE NEAREST O
! AND CHECKING WHETHER OR NOT THAT IS VIABLE
! NOW CONSIDERS THE POSSIBILITY THAT THERE COULD BE TWO APPROXIMATELY EQUALLY GOOD 'PIVOT' STATES
! IN THE CASE WHERE THE LOCAL GEOMETRY REPRESENTS A ZUNDEL SPECIES

     SUBROUTINE TRYSIMPLEVBASSIGNMENT(PROXIMALATOMS, NEIGHBOURCOUNTS, NUMSTATESASSIGNED, PIVOTVBSTATES)

     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

! SUBROUTINE ARGUMENTS

     INTEGER, INTENT(IN) :: PROXIMALATOMS(NATOMS, NUM_HYD)
     INTEGER, INTENT(IN) :: NEIGHBOURCOUNTS(NATOMS)
     INTEGER, INTENT(OUT) :: NUMSTATESASSIGNED
     INTEGER, INTENT(OUT) :: PIVOTVBSTATES(2,NATOMS)

! LOCAL VARIABLES

     INTEGER :: CURRENTATOM, NEARESTOREF, INSERTREF
     INTEGER :: ASSIGNEDHCOUNTS(NUM_EIG)
     INTEGER :: ASSIGNEDHS(NUM_EIG,3)

     LOGICAL :: PIVOTOASSIGNED

     LOGICAL :: SYMMETRICH(NATOMS)
     INTEGER :: CHECKH, SHAREDH, PIVOTO1REF, PIVOTO2REF
     DOUBLE PRECISION, PARAMETER :: EQUIDISTANCECUTOFF = 1.2D0
     DOUBLE PRECISION :: DISTANCERATIO

! SEE IF WE CAN IDENTIFY A SYMMETRIC H ATOM

     SYMMETRICH = .FALSE.

     DO CURRENTATOM = 1, NATOMS
        IF ((MOD(CURRENTATOM,3).NE.0) .AND. (NEIGHBOURCOUNTS(CURRENTATOM).GT.1)) THEN

           DISTANCERATIO = INTERATOMICR(CURRENTATOM,PROXIMALATOMS(CURRENTATOM,2))/ &
                                  INTERATOMICR(CURRENTATOM,PROXIMALATOMS(CURRENTATOM,1))

           IF (DISTANCERATIO.LE.EQUIDISTANCECUTOFF) SYMMETRICH(CURRENTATOM) = .TRUE.
        ENDIF
     ENDDO

! FIND THE NEAREST O ATOM TO EACH H

     NUMSTATESASSIGNED = 1
     PIVOTOASSIGNED = .FALSE.
     ASSIGNEDHCOUNTS = 0
     PIVOTVBSTATES = 0

     DO CURRENTATOM = 1, NATOMS
        IF (MOD(CURRENTATOM,3).NE.0) THEN

           NEARESTOREF = PROXIMALATOMS(CURRENTATOM,1)/3

           IF (ASSIGNEDHCOUNTS(NEARESTOREF).EQ.3) THEN
           NUMSTATESASSIGNED = 0
           EXIT
           ELSE IF (ASSIGNEDHCOUNTS(NEARESTOREF).EQ.2) THEN
           IF (PIVOTOASSIGNED) THEN
              NUMSTATESASSIGNED = 0
              EXIT
           ELSE
              ASSIGNEDHCOUNTS(NEARESTOREF) = ASSIGNEDHCOUNTS(NEARESTOREF) + 1
              ASSIGNEDHS(NEARESTOREF,ASSIGNEDHCOUNTS(NEARESTOREF)) = CURRENTATOM
              PIVOTOASSIGNED = .TRUE.    
           ENDIF                
           ELSE
           ASSIGNEDHCOUNTS(NEARESTOREF) = ASSIGNEDHCOUNTS(NEARESTOREF) + 1
           ASSIGNEDHS(NEARESTOREF,ASSIGNEDHCOUNTS(NEARESTOREF)) = CURRENTATOM
           ENDIF
        ENDIF
     ENDDO

     IF (NUMSTATESASSIGNED.EQ.1) THEN

! FILL THE VB ASSIGNMENT

        INSERTREF = 6

        DO CURRENTATOM = 1, NUM_EIG
           IF (ASSIGNEDHCOUNTS(CURRENTATOM).EQ.3) THEN    ! H3O+ MOLECULE
           PIVOTVBSTATES(NUMSTATESASSIGNED,1) = ASSIGNEDHS(CURRENTATOM,1)
           PIVOTVBSTATES(NUMSTATESASSIGNED,2) = ASSIGNEDHS(CURRENTATOM,2)
           PIVOTVBSTATES(NUMSTATESASSIGNED,3) = CURRENTATOM*3
           PIVOTVBSTATES(NUMSTATESASSIGNED,4) = ASSIGNEDHS(CURRENTATOM,3)
           ELSE    ! NORMAL H2O MOLECULE
           PIVOTVBSTATES(NUMSTATESASSIGNED,INSERTREF-1) = ASSIGNEDHS(CURRENTATOM,1)
           PIVOTVBSTATES(NUMSTATESASSIGNED,INSERTREF) = CURRENTATOM*3
           PIVOTVBSTATES(NUMSTATESASSIGNED,INSERTREF+1) = ASSIGNEDHS(CURRENTATOM,2)

           INSERTREF = INSERTREF + 3
           ENDIF
        ENDDO

! CHECK IF THERE IS ANOTHER, ALMOST EQUALLY GOOD PIVOT STATE WHICH NEEDS ASSIGNING

        PIVOTO1REF = PIVOTVBSTATES(1,3)/3
        SHAREDH = -1

        DO CHECKH = 1, 3
           IF (SYMMETRICH(ASSIGNEDHS(PIVOTO1REF,CHECKH))) THEN
           IF (SHAREDH.EQ.-1) THEN
              SHAREDH = CHECKH
           ELSE 
              SHAREDH = -1
              EXIT
           ENDIF
           ENDIF
        ENDDO
        
        IF (SHAREDH.NE.-1) THEN

           PIVOTO2REF = PROXIMALATOMS(ASSIGNEDHS(PIVOTO1REF,SHAREDH),2)/3
           NUMSTATESASSIGNED = NUMSTATESASSIGNED + 1

! ASSIGN ANOTHER PIVOT STATE

           PIVOTVBSTATES(NUMSTATESASSIGNED,1) = ASSIGNEDHS(PIVOTO2REF,1)
           PIVOTVBSTATES(NUMSTATESASSIGNED,2) = ASSIGNEDHS(PIVOTO2REF,2)
           PIVOTVBSTATES(NUMSTATESASSIGNED,3) = PIVOTO2REF*3
           PIVOTVBSTATES(NUMSTATESASSIGNED,4) = ASSIGNEDHS(PIVOTO1REF,SHAREDH)

           PIVOTVBSTATES(NUMSTATESASSIGNED,6) = PIVOTVBSTATES(1,3)
           IF (SHAREDH.EQ.1) THEN
           PIVOTVBSTATES(NUMSTATESASSIGNED,5) = ASSIGNEDHS(PIVOTO1REF,2)
           PIVOTVBSTATES(NUMSTATESASSIGNED,7) = ASSIGNEDHS(PIVOTO1REF,3)
           ELSE
           PIVOTVBSTATES(NUMSTATESASSIGNED,5) = ASSIGNEDHS(PIVOTO1REF,1)

           IF (SHAREDH.EQ.2) THEN
              PIVOTVBSTATES(NUMSTATESASSIGNED,7) = ASSIGNEDHS(PIVOTO1REF,3)
           ELSE
              PIVOTVBSTATES(NUMSTATESASSIGNED,7) = ASSIGNEDHS(PIVOTO1REF,2)              
           ENDIF
           ENDIF

           INSERTREF = 9

           DO CURRENTATOM = 1, NUM_EIG             
           IF ((CURRENTATOM.NE.PIVOTO1REF) .AND. (CURRENTATOM.NE.PIVOTO2REF)) THEN
              PIVOTVBSTATES(NUMSTATESASSIGNED,INSERTREF-1) = ASSIGNEDHS(CURRENTATOM,1)
              PIVOTVBSTATES(NUMSTATESASSIGNED,INSERTREF) = CURRENTATOM*3
              PIVOTVBSTATES(NUMSTATESASSIGNED,INSERTREF+1) = ASSIGNEDHS(CURRENTATOM,2)
              INSERTREF = INSERTREF + 3
           ENDIF
           ENDDO
        ENDIF

     ENDIF

     RETURN

     END SUBROUTINE

! ###############################################################################################

! ASSIGN THE PIVOT VB STATE BY MINIMISING THE SUM OF BONDED OH DISTANCES
! THE PIVOT O IS IDENTIFIED AS THE O WITH THE 3 CLOSEST HS

     SUBROUTINE ASSIGNPIVOTVBSTATE(PROXIMALATOMS, NEIGHBOURCOUNTS, SUCCESS, PIVOTVBSTATE)

     USE COMMONS
     USE MSEVB_COMMON

     IMPLICIT NONE

! SUBROUTINE ARGUMENTS

     INTEGER, INTENT(IN OUT) :: PROXIMALATOMS(NATOMS, NUM_HYD)
     INTEGER, INTENT(IN OUT) :: NEIGHBOURCOUNTS(NATOMS)
     LOGICAL, INTENT(OUT) :: SUCCESS
     INTEGER, INTENT(OUT) :: PIVOTVBSTATE(NATOMS)

! LOCAL VARIABLES

     INTEGER :: I, J, K, L
     INTEGER :: CURRENTO, CURRENTH, PIVOTOXYGEN
     INTEGER :: ASSIGNEDOCOUNT
     INTEGER :: OLIST(NUM_EIG)

     DOUBLE PRECISION :: MINSUMROH, CURRENTSUMROH

     LOGICAL :: USEDH(NATOMS), GOINGUP, FOUNDANH
     INTEGER :: HCOUNTERS(NUM_EIG,3)

     SUCCESS = .TRUE.

! ORDER THE H ATOMS IN TERMS OF THE DISTANCE FROM EACH O, ONLY CONSIDER THOSE WHICH ARE LESS THAN
! THE CUTOFF DISTANCE APART

     DO I=1,NUM_EIG

        CURRENTO = 3*I
        NEIGHBOURCOUNTS(CURRENTO) = 0   

        DO J=1,NUM_HYD

           CURRENTH = J + INT((J-1)/2)

           IF (INTERATOMICR(CURRENTO,CURRENTH).LT.MAXHBONDLENGTH) THEN
           NEIGHBOURCOUNTS(CURRENTO) = NEIGHBOURCOUNTS(CURRENTO) + 1

           DO K=1,NEIGHBOURCOUNTS(CURRENTO)
              IF (K.EQ.NEIGHBOURCOUNTS(CURRENTO)) THEN
                 PROXIMALATOMS(CURRENTO,K) = CURRENTH
                 EXIT
              ELSEIF (INTERATOMICR(CURRENTO,CURRENTH).LT.INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,K))) THEN
                 DO L=NEIGHBOURCOUNTS(CURRENTO),K+1,-1
                 PROXIMALATOMS(CURRENTO,L)=PROXIMALATOMS(CURRENTO,L-1)
                 ENDDO
                 PROXIMALATOMS(CURRENTO,K) = CURRENTH
                 EXIT
              ENDIF
           ENDDO
           ENDIF
        ENDDO

! CHECK WHETHER WE HAVE ASSIGNED SUFFICIENT HS WITHIN THE CUTOFF TO FORM AT LEAST A WATER MOLECULE

        IF (NEIGHBOURCOUNTS(CURRENTO).LT.2) THEN
           PRINT *, 'WITHIN THE CUTOFF OF', MAXHBONDLENGTH, 'PROGRAM ASSIGNED', NEIGHBOURCOUNTS(CURRENTO), &
                       'AS POSSIBLY BONDED TO O ATOM', CURRENTO
           PRINT *, 'REQUIRE A MINIMUM OF 2 ASSIGNED HS TO FORM A WATER SPECIES'
           PRINT *, 'INCREASE CUTOFF OR EXAMINE YOUR GEOMETRY'
           SUCCESS = .FALSE.
           RETURN
        ENDIF
     ENDDO

! FIND THE PIVOT OXYGEN, THE O ATOM WITH THE THREE CLOSEST HS

     MINSUMROH = HUGE(MINSUMROH)
     PIVOTOXYGEN = -1

     DO CURRENTO = 3, (NATOMS-1), 3
        IF (NEIGHBOURCOUNTS(CURRENTO).GE.3) THEN

           CURRENTSUMROH = INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,1)) + &
                           INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,2)) + &
                              INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,3))

           IF (CURRENTSUMROH.LT.MINSUMROH) THEN
           PIVOTOXYGEN = CURRENTO
           MINSUMROH = CURRENTSUMROH
           ENDIF
        ENDIF
     ENDDO

! WE MUST HAVE ASSIGNED A PIVOT O BECAUSE EVERY H HAS BEEN ASSIGNED AT LEAST ONE O
! ASSIGN THE REST OF THIS VB STATE

     PIVOTVBSTATE = -1

! ORDER THE OTHER OXYGENS IN DECREASING PROXIMITY FROM THE HYDRONIUM O
! THIS SHOULD BE A MORE EFFICIENT WAY OF FINDING THE BEST OVERALL BOND ASSIGNMENT

     OLIST = -1
     OLIST(1) = PIVOTOXYGEN

     ASSIGNEDOCOUNT = 1

     DO I=1,NUM_EIG
        CURRENTO = 3*I

        IF (CURRENTO.EQ.PIVOTOXYGEN) CYCLE
           
        ASSIGNEDOCOUNT = ASSIGNEDOCOUNT + 1

        DO J=2,ASSIGNEDOCOUNT
           
           IF (OLIST(J).EQ.-1) THEN
           OLIST(J) = CURRENTO
           EXIT
           ELSEIF (INTERATOMICR(PIVOTOXYGEN,CURRENTO).LT.INTERATOMICR(PIVOTOXYGEN,OLIST(J))) THEN
           
           DO K=ASSIGNEDOCOUNT,J+1,-1
              OLIST(K) = OLIST(K-1)
           ENDDO

           OLIST(J) = CURRENTO
           EXIT                 
           ENDIF
        ENDDO
     ENDDO

     USEDH = .FALSE.

! SET THE COUNTERS TO THEIR INITIAL VALUES

     HCOUNTERS(1,1) = 1
     HCOUNTERS(1,2) = 2
     HCOUNTERS(1,3) = 3

     ASSIGNEDOCOUNT = 1
     GOINGUP = .TRUE.

! ASSIGN THE BEST VB SPECIES BY MINIMISING THE SUM OF THE DISTANCES BETWEEN BONDED ATOMS
 
! THREE HS WILL BE ASSOCIATED WITH THE PIVOT O AND TWO WITH EACH OF THE OTHERS

     MINSUMROH = HUGE(MINSUMROH)
     CURRENTSUMROH = 0.0D0
           
     INFINITELOOP: DO
        
        CURRENTO = OLIST(ASSIGNEDOCOUNT)

        IF (ASSIGNEDOCOUNT.EQ.1) THEN ! THE EIGEN O IN THIS STATE, WHICH HAS 3 ASSOCIATED HS

! CHECK IF WE ARRIVED HERE AS A RESULT OF INCREASING OR DECREASING ASSIGNOCOUNT
              
           IF (GOINGUP) THEN     ! INITIALISATION
           USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1))) = .TRUE.
           USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2))) = .TRUE.
           USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3))) = .TRUE.

           CURRENTSUMROH = INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1)))+&
                                    INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2)))+& 
                                 INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3)))

           ASSIGNEDOCOUNT = ASSIGNEDOCOUNT + 1

           ELSE          ! INCREMENT THE COUNTERS IF POSSIBLE, EXIT IF ALL ITERATIONS ARE DONE

           IF (HCOUNTERS(ASSIGNEDOCOUNT,3).EQ.NEIGHBOURCOUNTS(CURRENTO)) THEN
              
              IF (HCOUNTERS(ASSIGNEDOCOUNT,2).EQ.(NEIGHBOURCOUNTS(CURRENTO)-1)) THEN
                 
                 IF (HCOUNTERS(ASSIGNEDOCOUNT,1).EQ.(NEIGHBOURCOUNTS(CURRENTO)-2)) THEN
                 EXIT INFINITELOOP
                 ELSE
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1))) = .FALSE.
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2))) = .FALSE.
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3))) = .FALSE.
                 HCOUNTERS(ASSIGNEDOCOUNT,1) = HCOUNTERS(ASSIGNEDOCOUNT,1) + 1
                 HCOUNTERS(ASSIGNEDOCOUNT,2) = HCOUNTERS(ASSIGNEDOCOUNT,1) + 1
                 HCOUNTERS(ASSIGNEDOCOUNT,3) = HCOUNTERS(ASSIGNEDOCOUNT,2) + 1
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1))) = .TRUE.
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2))) = .TRUE.
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3))) = .TRUE.

                 CURRENTSUMROH = &
                               INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1))) + &
                               INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2))) + &
                               INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3)))
       
                 IF (CURRENTSUMROH.GE.MINSUMROH) THEN
                    GOINGUP = .FALSE.
                 ELSE
                    GOINGUP = .TRUE.
                 ENDIF
                 ENDIF

              ELSE 
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2))) = .FALSE.
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3))) = .FALSE.
                 HCOUNTERS(ASSIGNEDOCOUNT,2) = HCOUNTERS(ASSIGNEDOCOUNT,2) + 1
                 HCOUNTERS(ASSIGNEDOCOUNT,3) = HCOUNTERS(ASSIGNEDOCOUNT,2) + 1
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2))) = .TRUE.
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3))) = .TRUE.

                 CURRENTSUMROH = &
                            INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1))) + &
                            INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2))) + &
                            INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3)))

                 IF (CURRENTSUMROH.GE.MINSUMROH) THEN
                 GOINGUP = .FALSE.
                 ELSE
                 GOINGUP = .TRUE.
                 ENDIF
                 
              ENDIF
           ELSE 
              USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3))) = .FALSE.
              HCOUNTERS(ASSIGNEDOCOUNT,3) = HCOUNTERS(ASSIGNEDOCOUNT,3) + 1
              USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3))) = .TRUE.       

              CURRENTSUMROH = &
                         INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1))) + &
                         INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2))) + &
                         INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,3)))

              IF (CURRENTSUMROH.GE.MINSUMROH) THEN
                 GOINGUP = .FALSE.
              ELSE
                 GOINGUP = .TRUE.
              ENDIF
              
           ENDIF

           IF (GOINGUP) ASSIGNEDOCOUNT = ASSIGNEDOCOUNT + 1

           ENDIF

        ELSE               ! A WATER O, HAS 2 ASSOCIATED HS

! CHECK WHETHER WE HAVE ARRIVED HERE AS A RESULT OF INCREASING OR DECREASING ASSIGNEDOCOUNT

           IF (GOINGUP) THEN ! INITIALISATION

!           PRINT *, 'ADDING A NEW WATER O:', CURRENTO
!                 PRINT *, 'ASSIGNMENT TO DATE:'
!                 DO I=1,ASSIGNEDOCOUNT-1
!                 PRINT *, 'O:', (3*OLIST(I))
!                 IF (I.EQ.1) THEN
!                    PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                               'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2))),
!    &                             'H3:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,3)))
!                 ELSE
!                    PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                             'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2)))
!                 ENDIF
!                 ENDDO

!                 PRINT *, 'ROH_SQ TO DATE:', CURRENTSUMROHSQ

           FOUNDANH = .FALSE.

           FIRSTHLOOP: DO I=1, NEIGHBOURCOUNTS(CURRENTO)-1

              IF (.NOT.USEDH(PROXIMALATOMS(CURRENTO,I))) THEN

! IF ADDING THIS OH DISTANCE BRINGS THE TOTAL TO MORE THAN THE CURRENT MINIMUM THEN ANY SUBSEQUENT H WILL ALSO DO
! SO SINCE THEY ARE ARRANGED IN ASCENDING DISTANCE FROM THE O

                 IF ((CURRENTSUMROH+INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I))).GE.MINSUMROH) THEN
                 EXIT FIRSTHLOOP
                 ELSE

                 DO J = I+1, NEIGHBOURCOUNTS(CURRENTO)
                    IF (.NOT.USEDH(PROXIMALATOMS(CURRENTO,J))) THEN

                    IF ((CURRENTSUMROH + INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I)) + &
                                     INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,J))).GE.MINSUMROH) THEN
                       EXIT FIRSTHLOOP
                    ELSE
                       HCOUNTERS(ASSIGNEDOCOUNT,1) = I
                       USEDH(PROXIMALATOMS(CURRENTO,I)) = .TRUE.
                       HCOUNTERS(ASSIGNEDOCOUNT,2) = J
                       USEDH(PROXIMALATOMS(CURRENTO,J)) = .TRUE.
                       FOUNDANH = .TRUE.

                       CURRENTSUMROH = CURRENTSUMROH + &
                                        INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I)) + &
                                        INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,J))
                       EXIT FIRSTHLOOP
                    ENDIF
                    ENDIF
                 ENDDO
                 ENDIF
              ENDIF
           ENDDO FIRSTHLOOP

           IF (FOUNDANH) THEN

!              PRINT *, 'SUCCESSFULLY ALLOCATED NEW HS'
              IF (ASSIGNEDOCOUNT.EQ.NUM_EIG) THEN

                 PIVOTVBSTATE(1) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,1))
                 PIVOTVBSTATE(2) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,2))
                 PIVOTVBSTATE(3) = PIVOTOXYGEN
                 PIVOTVBSTATE(4) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,3))

                 I = 5

                 DO J = 2, NUM_EIG
                 PIVOTVBSTATE(I) = PROXIMALATOMS(OLIST(J),HCOUNTERS(J,1))
                 PIVOTVBSTATE(I+1) = OLIST(J)
                 PIVOTVBSTATE(I+2) = PROXIMALATOMS(OLIST(J),HCOUNTERS(J,2))
                 I = I + 3
                 ENDDO

                 MINSUMROH = CURRENTSUMROH

!                 PRINT *, 'NEW ASSIGNMENT:', MINSUMROH
!                    DO I = 1,NUM_EIG
!                    PRINT *, 'O:', (3*OLIST(I))
!                    IF (I.EQ.1) THEN
!                       PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                             'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2))),
!    &                             'H3:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,3)))
!                    ELSE
!                       PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                                'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2)))
!                    ENDIF
!                    ENDDO

                 GOINGUP = .FALSE.
              ELSE
                 ASSIGNEDOCOUNT = ASSIGNEDOCOUNT + 1
              ENDIF
           ELSE
!              PRINT *, 'COULD NOT ALLOCATE NEW HS, DROPPING A LEVEL'
              ASSIGNEDOCOUNT = ASSIGNEDOCOUNT - 1
              GOINGUP = .FALSE.
           ENDIF   
           ELSE 

!           PRINT *, 'ATTEMPTING TO INCREMENT O:', CURRENTO
                 
!                 PRINT *, 'ASSIGNMENT AT THIS POINT:'
                 
!                 DO I=1,ASSIGNEDOCOUNT
!                 PRINT *, 'O:', (3*OLIST(I))
!                 IF (I.EQ.1) THEN
!                    PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                               'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2))),
!    &                             'H3:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,3)))
!                 ELSE
!                    PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                             'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2)))
!                 ENDIF
!                 ENDDO

!                 PRINT *, 'CURRENT SUM ROH_SQ:', CURRENTSUMROHSQ

! INCREMENT THE CURRENT COUNTER IF POSSIBLE, ELSE DROP DOWN ANOTHER LEVEL

           CURRENTSUMROH = CURRENTSUMROH - &
                      INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2)))

           USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,2))) = .FALSE.

           IF (HCOUNTERS(ASSIGNEDOCOUNT,2).EQ.NEIGHBOURCOUNTS(CURRENTO)) THEN                 
              IF (HCOUNTERS(ASSIGNEDOCOUNT,1).EQ.(NEIGHBOURCOUNTS(CURRENTO)-1)) THEN                 

                 CURRENTSUMROH = CURRENTSUMROH - & 
                            INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1)))

                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1))) = .FALSE.
!                    HCOUNTERS(ASSIGNEDOCOUNT,1) = 1
!                    HCOUNTERS(ASSIGNEDOCOUNT,2) = 2
                 ASSIGNEDOCOUNT = ASSIGNEDOCOUNT - 1
              ELSE
                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1))) = .FALSE.
                 CURRENTSUMROH = CURRENTSUMROH - &
                            INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1)))

! FIND A SUITABLE H ATOM FOR THE FIRST H ATOM

                 FOUNDANH = .FALSE.

                 FIRSTHLOOP2: DO I = HCOUNTERS(ASSIGNEDOCOUNT,1)+1, NEIGHBOURCOUNTS(CURRENTO)-1
                    IF (.NOT.USEDH(PROXIMALATOMS(CURRENTO,I))) THEN

                    IF ((CURRENTSUMROH + INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I))).GE.MINSUMROH) THEN
                    EXIT FIRSTHLOOP2
                    ELSE
                    DO J = I+1, NEIGHBOURCOUNTS(CURRENTO) 
                       IF (.NOT.USEDH(PROXIMALATOMS(CURRENTO,J))) THEN
                         
                          IF ((CURRENTSUMROH + &
                                           INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I)) + &
                                           INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,J))).GE.MINSUMROH) THEN
                          EXIT FIRSTHLOOP2
                          ELSE
                          HCOUNTERS(ASSIGNEDOCOUNT,1) = I
                          USEDH(PROXIMALATOMS(CURRENTO,I)) = .TRUE.
                          HCOUNTERS(ASSIGNEDOCOUNT,2) = J
                          USEDH(PROXIMALATOMS(CURRENTO,J)) = .TRUE.
                          FOUNDANH = .TRUE.

                          CURRENTSUMROH = CURRENTSUMROH + &
                                              INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I)) + &
                                              INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,J))

                          EXIT FIRSTHLOOP2
                          ENDIF
                       ENDIF
                    ENDDO
                    ENDIF
                 ENDIF
                 ENDDO FIRSTHLOOP2

                 IF (FOUNDANH) THEN
                 IF (ASSIGNEDOCOUNT.EQ.NUM_EIG) THEN

                    PIVOTVBSTATE(1) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,1))
                    PIVOTVBSTATE(2) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,2))
                    PIVOTVBSTATE(3) = PIVOTOXYGEN
                    PIVOTVBSTATE(4) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,3))

                    I = 5
                       
                    DO J = 2, NUM_EIG
                    PIVOTVBSTATE(I) = PROXIMALATOMS(OLIST(J),HCOUNTERS(J,1))
                    PIVOTVBSTATE(I+1) = OLIST(J)
                    PIVOTVBSTATE(I+2) = PROXIMALATOMS(OLIST(J),HCOUNTERS(J,2))
                    I = I + 3
                    ENDDO

                    MINSUMROH = CURRENTSUMROH

!                    PRINT *, 'NEW ASSIGNMENT:', MINSUMROH
!                       DO I = 1,NUM_EIG
!                          PRINT *, 'O:', (3*OLIST(I))
!                          IF (I.EQ.1) THEN
!                          PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                                   'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2))),
!    &                                   'H3:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,3)))
!                          ELSE
!                          PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                                      'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2)))
!                          ENDIF
!                       ENDDO
                 ELSE
                    ASSIGNEDOCOUNT = ASSIGNEDOCOUNT + 1
                    GOINGUP = .TRUE.
                 ENDIF
                 ELSE 

! DROP DOWN A LEVEL
                 ASSIGNEDOCOUNT = ASSIGNEDOCOUNT - 1
                 ENDIF
              ENDIF
           ELSE 

! CHECK FOR AN UNUSED H ATOM
                 
              FOUNDANH = .FALSE.
              
              DO I = HCOUNTERS(ASSIGNEDOCOUNT,2)+1,NEIGHBOURCOUNTS(CURRENTO)
                 IF (.NOT.USEDH(PROXIMALATOMS(CURRENTO,I))) THEN

                 IF ((CURRENTSUMROH + INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I))).GE.MINSUMROH) THEN
                    EXIT
                 ELSE
                    HCOUNTERS(ASSIGNEDOCOUNT,2) = I
                    USEDH(PROXIMALATOMS(CURRENTO,I)) = .TRUE.
                    FOUNDANH = .TRUE.
                    CURRENTSUMROH = CURRENTSUMROH + INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I))
                    EXIT
                 ENDIF
                 ENDIF
              ENDDO

              IF (FOUNDANH) THEN

                 IF (ASSIGNEDOCOUNT.EQ.NUM_EIG) THEN

                 PIVOTVBSTATE(1) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,1))
                 PIVOTVBSTATE(2) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,2))
                 PIVOTVBSTATE(3) = PIVOTOXYGEN
                 PIVOTVBSTATE(4) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,3))

                 I = 5

                 DO J = 2, NUM_EIG
                    PIVOTVBSTATE(I) = PROXIMALATOMS(OLIST(J),HCOUNTERS(J,1))
                    PIVOTVBSTATE(I+1) = OLIST(J)
                    PIVOTVBSTATE(I+2) = PROXIMALATOMS(OLIST(J),HCOUNTERS(J,2))
                    I = I + 3
                 ENDDO

                 MINSUMROH = CURRENTSUMROH

!                 PRINT *, 'NEW ASSIGNMENT:', MINSUMROH
!                    DO I = 1,NUM_EIG
!                       PRINT *, 'O:', (3*OLIST(I))

!                       IF (I.EQ.1) THEN
!                          PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                                         'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2))),
!    &                                         'H3:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,3)))
!                       ELSE
!                          PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                                         'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2)))
!                       ENDIF
!                    ENDDO
                 ELSE
                 ASSIGNEDOCOUNT = ASSIGNEDOCOUNT + 1
                 GOINGUP = .TRUE.
                 ENDIF

! COULD NOT FIND AN ACCEPTABLE SECOND H, TRY RUNNING THROUGH THE FIRST H SELECTION IN TANDEM

              ELSE     

                 USEDH(PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1))) = .FALSE.
                 CURRENTSUMROH = CURRENTSUMROH - &
                            INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,HCOUNTERS(ASSIGNEDOCOUNT,1)))

                 FOUNDANH = .FALSE.

                 FIRSTHLOOP3: DO I = HCOUNTERS(ASSIGNEDOCOUNT,1)+1, NEIGHBOURCOUNTS(CURRENTO)-1
                    IF (.NOT.USEDH(PROXIMALATOMS(CURRENTO,I))) THEN

                    IF ((CURRENTSUMROH + INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I))).GE.MINSUMROH) THEN
                    EXIT FIRSTHLOOP3
                    ELSE

                    DO J = I+1, NEIGHBOURCOUNTS(CURRENTO) 
                       IF (.NOT.USEDH(PROXIMALATOMS(CURRENTO,J))) THEN

                          IF ((CURRENTSUMROH + &
                                           INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I)) + &
                                           INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,J))).GE.MINSUMROH) THEN
                          EXIT FIRSTHLOOP3
                          ELSE

                          HCOUNTERS(ASSIGNEDOCOUNT,1) = I
                          USEDH(PROXIMALATOMS(CURRENTO,I)) = .TRUE.
                          HCOUNTERS(ASSIGNEDOCOUNT,2) = J
                          USEDH(PROXIMALATOMS(CURRENTO,J)) = .TRUE.

                          CURRENTSUMROH = CURRENTSUMROH + &
                                              INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,I)) + &
                                              INTERATOMICR(CURRENTO,PROXIMALATOMS(CURRENTO,J))
                          FOUNDANH = .TRUE.
                          EXIT FIRSTHLOOP3
                          ENDIF
                       ENDIF
                    ENDDO
                    ENDIF
                 ENDIF
                 ENDDO FIRSTHLOOP3

                 IF (FOUNDANH) THEN
                 IF (ASSIGNEDOCOUNT.EQ.NUM_EIG) THEN

                    PIVOTVBSTATE(1) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,1))
                    PIVOTVBSTATE(2) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,2))
                    PIVOTVBSTATE(3) = PIVOTOXYGEN
                    PIVOTVBSTATE(4) = PROXIMALATOMS(PIVOTOXYGEN, HCOUNTERS(1,3))

                    I = 5

                    DO J = 2, NUM_EIG
                    PIVOTVBSTATE(I) = PROXIMALATOMS(OLIST(J),HCOUNTERS(J,1))
                    PIVOTVBSTATE(I+1) = OLIST(J)
                    PIVOTVBSTATE(I+2) = PROXIMALATOMS(OLIST(J),HCOUNTERS(J,2))
                    I = I + 3
                    ENDDO

                    MINSUMROH = CURRENTSUMROH

!                    PRINT *, 'NEW ASSIGNMENT:', MINSUMROH
!                       DO I = 1,NUM_EIG
!                          PRINT *, 'O:', (3*OLIST(I))
!                          IF (I.EQ.1) THEN
!                          PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                                    'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2))),
!    &                                      'H3:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,3)))
!                          ELSE
!                          PRINT *, 'H1:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,1))),
!    &                                      'H2:', HPOSITION(NEIGHBOURINGHS(OLIST(I),HCOUNTERS(I,2)))
!                          ENDIF
!                       ENDDO
                 ELSE
                    ASSIGNEDOCOUNT = ASSIGNEDOCOUNT + 1                 
                    GOINGUP = .TRUE.
                 ENDIF
                 ELSE 

! DROP DOWN A LEVEL

                 ASSIGNEDOCOUNT = ASSIGNEDOCOUNT - 1
                 ENDIF                 
              ENDIF
           ENDIF
           ENDIF
        ENDIF
     ENDDO INFINITELOOP

! CHECK ALSO THAT IT WAS POSSIBLE TO ASSIGN A COMPLETE VB STATE

     DO I=1,NATOMS
        IF (PIVOTVBSTATE(I).EQ.-1) THEN
           PRINT *, 'UNABLE TO ASSIGN A COMPLETE VB STATE FOR PIVOT OXYGEN', PIVOTOXYGEN
           PRINT *, 'CHECK GEOMETRY OR CUTOFF'
           SUCCESS = .FALSE.
           EXIT
        ENDIF
     ENDDO

     RETURN

     END SUBROUTINE

! ###############################################################################################

! CALCULATE INTERATOMIC DISTANCES - THESE WILL BE USED BY THE ENERGY/FORCE CALCULATIONS

     SUBROUTINE CALCINTERATOMICDISTANCES(NATOMS,XCOORDS,YCOORDS,ZCOORDS,SEPARATIONS)
     IMPLICIT NONE

! SUBROUTINE ARGUMENTS

     INTEGER, INTENT(IN) :: NATOMS
     DOUBLE PRECISION, INTENT (IN) :: XCOORDS(NATOMS), YCOORDS(NATOMS), ZCOORDS(NATOMS)
     DOUBLE PRECISION, DIMENSION(NATOMS,NATOMS), INTENT (OUT) :: SEPARATIONS

! LOCAL VARIABLES

     INTEGER ATOM1, ATOM2

! INITIALISE

     SEPARATIONS = 0.0D0

! CALCULATE THE UPPER TRIANGLE BUT EXPLICITLY FILL IN THE WHOLE MATRIX TO AVOID IN FUTURE HAVING TO WORK OUT
! WHICH ATOM NUMBER IS LOWER

     DO ATOM1 = 1, (NATOMS-1)
        DO ATOM2 = (ATOM1+1), NATOMS
           SEPARATIONS(ATOM1,ATOM2) = DSQRT((XCOORDS(ATOM1)-XCOORDS(ATOM2))**2 + &
                                               (YCOORDS(ATOM1)-YCOORDS(ATOM2))**2 + &
                                               (ZCOORDS(ATOM1)-ZCOORDS(ATOM2))**2)
           SEPARATIONS(ATOM2,ATOM1) = SEPARATIONS(ATOM1,ATOM2)
        ENDDO
     ENDDO

     RETURN

     END SUBROUTINE

! ##################################################################################################     

