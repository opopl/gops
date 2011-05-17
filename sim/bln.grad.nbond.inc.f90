
! ..... Non-bonded interaction forces ..... {{{

        DO I = 1, N-2
          DO J = I+2, N
  
            RAD(1)=LEN_DR(I,J) ; 
            RAD(2)=RAD(1)**2 ; 
            RAD(7)=RAD(1)**7 ; 
            RAD(14)=RAD(7)**2 ; 
            RAD(8)=RAD(7)*RAD(1)
        
            DF = 2.0*AB(I,J,1)*S(12)/RAD(12) + AB(I,J,2)*S(6)/RAD(6)
            DF=-24.0*DF/RAD(2)

            FRR(1:3) = DF*DR(I,J,1:3) 
            FNB(I,1:3) = FRR(1:3) + FNB(I,1:3)
            FNB(J,1:3) = -FRR(1:3) + FNB(J,1:3)

          ENDDO
        ENDDO
! }}}

