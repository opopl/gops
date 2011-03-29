!==================================================
! CALCULATE THE MONOMER AND DIMER POTENTIALS, USING HBB ONLY
!==================================================
SUBROUTINE POT12BHBB(NATM,XX,POT)
  INTEGER,INTENT(IN)::NATM
  REAL,DIMENSION(3,NATM),INTENT(IN)::XX
  REAL,INTENT(INOUT)::POT
  !::::::::::::::::::::
  REAL,DIMENSION(6,3)::X1,X2
  REAL::E1,E2,VECT(3)
  INTEGER::I,J,FO

  FO=NATM/3*2
  DO I=1,NATM/3-1
     X2(1,:)=XX(:,FO+I)        ! O1
     X2(2,:)=XX(:,I*2-1)       ! H1
     X2(3,:)=XX(:,I*2)         ! H2
     DO J=I+1,NATM/3
        X2(4,:)=XX(:,FO+J)     ! O2
        X2(5,:)=XX(:,J*2-1)    ! H1
        X2(6,:)=XX(:,J*2)      ! H2
        VECT(:)=X2(4,:)-X2(1,:)
        X1=X2
        VECT(:)=VECT(:)*200.D0
        X1(4,:)=X2(4,:)+VECT(:)
        X1(5,:)=X2(5,:)+VECT(:)
        X1(6,:)=X2(6,:)+VECT(:)
        CALL CALCPOT(E2,X2) 
        CALL CALCPOT(E1,X1) 
        POT=POT+E2-E1*(1.D0-1.D0/DBLE(NATM/3-1))
     END DO
  END DO

  RETURN
END SUBROUTINE POT12BHBB
