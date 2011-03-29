MODULE PES_SHELL
  USE HB3B_COEF,ONLY:CONN
  IMPLICIT NONE

  ! DEFINE CONSTANTS
  REAL,PARAMETER::EMASS=1822.88848
  REAL,PARAMETER::AUANG=0.5291772083
  REAL,PARAMETER::AUCM=219474.6313710
  REAL,PARAMETER::AUKCAL=627.51

  ! WATER FUNCTIONS MAP
  INTEGER,SAVE::IWATERFCN
  ! 1: HBB12
  ! 2: HBB12+3B
  ! 3: HBB12+HB3B

CONTAINS
  !=================================================!
  ! THIS SUBROUTINE IS USED TO INITALIZE THE PES    !
  !=================================================!
  SUBROUTINE PES_INIT(NW, DNAME)
    INTEGER,INTENT(IN)::NW
    CHARACTER (LEN=*), INTENT(IN) :: DNAME
    !::::::::::::::::::::
    INTEGER::I
    INTEGER,DIMENSION(NW)::IDX_O

    ! 3-BODY POT
    CALL PES_INIT_3B(DNAME)

    ! DIMER POT
    CALL PREPOT()

    ! H-BONDED 3B
    ALLOCATE(CONN(NW*2,NW-1))
    IDX_O=(/(I,I=1,NW)/)
    DO I=1,NW
       IDX_O(I)=0
       CONN(I*2-1,:)=PACK(IDX_O,MASK=IDX_O.NE.0)
       CONN(I*2,:)=CONN(I*2-1,:)
       IDX_O(I)=I
    END DO

    RETURN
  END SUBROUTINE PES_INIT

  !==================================================!
  ! WATER POTENTIALS                                 !
  !==================================================!
  FUNCTION F(X) RESULT(POT)
    REAL,DIMENSION(:,:),INTENT(IN)::X
    REAL::POT
    ! ::::::::::::::::::::
    REAL,DIMENSION(3,1:SIZE(X,2))::XNEW,XXA
    REAL::P1,P2,P3
    INTEGER::NATM,NW,I

    NATM=SIZE(X,2); NW=NATM/3
    XNEW=X
    XXA=X*AUANG

    P1=0.D0; P2=0.D0; P3=0.D0
    SELECT CASE(IWATERFCN)
    CASE(1) ! 1: HBB12
       CALL POT12BHBB(NATM,XNEW,P2)
    CASE(2) ! 2: HBB12+3B
       CALL POT12BHBB(NATM,XNEW,P2)
       CALL POT3B(NATM,XNEW,P3)          
    CASE(3) ! 3: HBB12+HB3B
       CALL POT12BHBB(NATM,XNEW,P2)
       CALL POT_HB3B(NW,XXA,P3)
    END SELECT

    POT=P3+P2+P1

    RETURN
  END FUNCTION F
  
  !==================================================
  ! CALCULATE THE THREE BODY POTENTIALS
  !==================================================
  SUBROUTINE POT3B(NATM,XX,POT)
    INTEGER,INTENT(IN)::NATM
    REAL,DIMENSION(3,NATM),INTENT(IN)::XX
    REAL,INTENT(INOUT)::POT
    !::::::::::::::::::::
    REAL,DIMENSION(3,9)::X3
    REAL::E3
    REAL,EXTERNAL::FPES
    INTEGER::I,J,K,FO

    FO=NATM/3*2;E3=0.D0;POT=0.D0
    DO I=1,NATM/3-2
       X3(:,7)=XX(:,FO+I)          ! O1
       X3(:,1)=XX(:,I*2-1)         ! H1
       X3(:,2)=XX(:,I*2)           ! H1'
       DO J=I+1,NATM/3-1
          X3(:,8)=XX(:,FO+J)       ! O2
          X3(:,3)=XX(:,J*2-1)      ! H2
          X3(:,4)=XX(:,J*2)        ! H2'
          DO K=J+1,NATM/3
             X3(:,9)=XX(:,FO+K)       ! O3
             X3(:,5)=XX(:,K*2-1)      ! H3
             X3(:,6)=XX(:,K*2)        ! H3'
             E3=FPES(X3)      ! THREE-BODY CORRECTION
             POT=POT+E3
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE POT3B

END MODULE PES_SHELL
