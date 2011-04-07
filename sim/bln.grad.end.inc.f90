
             F=FNB+FB+FBA+FTA 

        DO I = 1, N
          J = (I-1)*3
          FQ(J+1) = -F(I,1)
          FQ(J+2) = -F(I,2)
          FQ(J+3) = -F(I,3)
        ENDDO

