        
             F=FNB+FB+FBA+FTA 

        DO I = 1, N
          J = (I-1)*3
          GRAD(J+1) = -F(I,1)
          GRAD(J+2) = -F(I,2)
          GRAD(J+3) = -F(I,3)
        ENDDO

