! Parameters for the dihedral angle potential: in array CD(:,1:2) {{{

        DO I = 1, N-3
          ICOUNT = 0

          DO J = 0,3
            IF(NTYPE(I+J) .EQ. 3)THEN
              ICOUNT = ICOUNT + 1
            ENDIF
          ENDDO

          IF (ICOUNT .GE. 2) THEN
            CD(I+1,1:2) = (/ 0.0, 0.2*EPSILON /)
          ELSE
            CD(I+1,1:2) = 1.2*EPSILON
        ENDIF

        ICOUNT = 0

        ENDDO

