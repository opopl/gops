! Parameters for the dihedral angle potential: in array CD(:,1:2) {{{

        DO I = 1, N-3
          ICOUNT = 0

          DO J = 0,3
            ! count the neutral residues next in the chain:
            ! for i,i+1,i+2,i+3
            IF(NTYPE(I+J) .EQ. 3)THEN
              ICOUNT = ICOUNT + 1
            ENDIF
          ENDDO

          IF (ICOUNT .GE. 2) THEN
            ! if there are two or more neutral residues involved 
            ! in the definition of the dihedral angle, then...
            CD(I+1,1:2) = (/ 0.0D0, 0.2D0*EPSILON /)
            ELSE
              ! ... otherwise
            CD(I+1,1:2) = 1.2*EPSILON
          ENDIF

        ICOUNT = 0

        ENDDO

