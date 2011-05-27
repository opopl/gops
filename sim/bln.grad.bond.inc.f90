! ... Bond interaction forces ... {{{

        DO I = 1, N-1
           RVAR = SIGMA/B(I) ; DF = RK_R*(1.0 - RVAR) 
           FRR(1:3)=DF*BVR(I,1:3)
           GB(I,1:3) = FRR(1:3)+GB(I,1:3)             
           GB(I+1,1:3) = -FRR(1:3)+GB(I+1,1:3)
        ENDDO
! }}}

