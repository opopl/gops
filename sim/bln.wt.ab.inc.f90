        DO I = 1, N
          DO J = 1, N
            IF (NTYPE(I) .EQ. 3 .OR. NTYPE(J) .EQ. 3)THEN
              AB(I,J,1:2) = (/ 1.0*EPSILON, 0.0D0  /) 
            ELSEIF (NTYPE(I) .EQ. 1 .AND. NTYPE(J) .EQ. 1)THEN
              AB(I,J,1:2)= EPSILON*(/ 1.0D0, -1.0D0 /) 
            ELSE
              AB(I,J,1:2) = EPSILON*2.0/3.0 
            ENDIF
          ENDDO
        ENDDO

