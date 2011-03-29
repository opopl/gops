!==================================================!
! THIS SUBROUTINE CALCULATES 3-BODY USING SKINNER'S!
! FUNCTIONS.                                       !
!==================================================!  
SUBROUTINE POT_HB3B(NW,X,POT)
  USE HB3B_COEF
  INTEGER,INTENT(IN)::NW
  REAL,DIMENSION(3,NW*3),INTENT(IN)::X
  REAL,INTENT(OUT)::POT
  !::::::::::::::::::::
  REAL,DIMENSION(NW*2,NW)::RR,FF,GG,HH
  INTEGER::FO,I,J,K,JO,KO,IH1,IH2,JH1,JH2,KH1,KH2
  REAL::FA,FB,FC
  INTEGER::NA,NB,NC

  FA=0.D0;FB=0.D0;FC=0.D0
  FO=NW*2
  RR=0.D0

  ! H-BOND LENGTHES
  DO I=1,NW*2
     DO J=1,NW-1
        JO=CONN(I,J)
        RR(I,JO)=SQRT(SUM((X(:,I)-X(:,FO+JO))**2))
     END DO
  END DO

  ! MORSE_VARIABLES
  FF=EXP(-PARA_A2*RR)
  GG=EXP(-PARA_B2*RR)
  HH=EXP(-PARA_C2*RR)

  !TYPE-A
  NA=0;NB=0;NC=0;
  DO I=1,NW
     IH1=I*2-1
     IH2=I*2
     DO J=1,NW-2
        JO=CONN(IH1,J)
        JH1=JO*2-1
        JH2=JO*2
        DO K=J+1,NW-1
           KO=CONN(IH1,K)          
           KH1=KO*2-1
           KH2=KO*2
           ! TYPE-A
           FA=FA+FF(IH1,JO)*FF(IH2,KO)+FF(IH1,KO)*FF(IH2,JO)

           ! TYPE-B
           FB=FB+GG(IH1,JO)*GG(JH1,KO)+GG(IH2,JO)*GG(JH1,KO)&
                +GG(IH1,JO)*GG(JH2,KO)+GG(IH2,JO)*GG(JH2,KO)&
                +GG(IH1,KO)*GG(KH1,JO)+GG(IH2,KO)*GG(KH1,JO)&
                +GG(IH1,KO)*GG(KH2,JO)+GG(IH2,KO)*GG(KH2,JO)

           ! TYPE-C
           FC=FC+HH(JH1,I)*HH(KH1,I)+HH(JH1,I)*HH(KH2,I)&
                +HH(JH2,I)*HH(KH1,I)+HH(JH2,I)*HH(KH2,I)

!!$             NA=NA+2
!!$             NB=NB+8
!!$             NC=NC+4
        END DO
     END DO
  END DO

  POT=(FA*PARA_A1+FB*PARA_B1+FC*PARA_C1)*AUKJ

  RETURN
END SUBROUTINE POT_HB3B
