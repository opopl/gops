            E(1:5)=0.0D0 
        
! non-bonded: E(2) {{{

        DO I = 1, N-2
          DO J = I+2, N
            RAD(6)=LEN_DR(I,J)**6; 
            RAD(12)=RAD(6)**2
            E(2) = E(2) + 4.0D0*AB(I,J,1)*S(12)/RAD(12)
            E(2) = E(2) + 4.0D0*AB(I,J,2)*S(6)/RAD(6)
          ENDDO
        ENDDO
! }}}
! bonded: E(3) {{{

        DO I = 1, N-1
          E(3) = E(3) + 0.5*RK_R*(B(I)-SIGMA)**2
        ENDDO
! }}}
! bond angles: E(4) {{{

        DO I = 2, N-1
          E(4) = E(4) + 0.5*RK_THETA*(ANG(I,1)-THETA_0)**2
        ENDDO
! }}}
! torsional angles: E(5) {{{

        DO I = 2, N-2
          E(5) = E(5) + CD(I,1)*(1.0 + COS(ANG(I,2)))
          E(5) = E(5) + CD(I,2)*(1.0 + COS(3.0*ANG(I,2)))
        ENDDO
! }}}

        ! Now the total energy
        E(1)=SUM(E(2:5))

