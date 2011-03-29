!=================================================!
! THIS SUBROUTINE IS USED TO INITALIZE THE        !
!   PES FROM BAS                                  !
!=================================================!
SUBROUTINE PES_INIT_3B(DNAME)
  USE PES, ONLY: PES0_INIT,PES1_INIT
  CHARACTER (LEN=*) :: DNAME
  CHARACTER (LEN=*), PARAMETER:: PES_X6Y3_SYSALL=&
       'X1 Y1 X2 X1Y1 Y2 X3 X2Y1 X1Y2 Y3 X4 X3Y1 X2Y2 X1Y3 '// &
       'X5 X4Y1 X3Y2 X2Y3 X6 X5Y1 X4Y2 X3Y3 X6Y1 X5Y2 X4Y3 X5Y2 X4Y3 '// &
       'X6Y2 X5Y3 X6Y3'

  CALL PES0_INIT (DIR=TRIM(DNAME))
  CALL PES1_INIT (PES_X6Y3_SYSALL)

  RETURN
END SUBROUTINE PES_INIT_3B
!=========================================================!
! THIS IS THE MAIN FUNCTION TO GET THE POTENTIAL          !
! ENERGY FROM FITTING CODE                                !
!=========================================================!
FUNCTION FPES(X) 
  USE PES,ONLY:PES_X6Y3_POT
  REAL,DIMENSION(0:2,0:8),INTENT(IN)::X
  REAL::FPES
  ! ::::::::::::::::::::

  FPES=PES_X6Y3_POT(X)

  RETURN
END FUNCTION FPES
