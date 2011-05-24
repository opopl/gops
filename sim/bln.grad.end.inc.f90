        
             G=GNB+GB+GBA+GTA 

        DO I = 1, N
          J = (I-1)*3
          GRAD(J+1) = -G(I,1)
          GRAD(J+2) = -G(I,2)
          GRAD(J+3) = -G(I,3)
        ENDDO

