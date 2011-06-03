      ELSE IF (WORD.EQ.'ADDPATH') THEN
         ADDPATH=.TRUE.
         CALL READA(PATHNAME)
C
C  ADDPT determines whether we call ADDPERM to add permutational isomers of every stationary point to the
C  min.data and ts.data databases. Speeds up 2DLJ7! For other tagged atom situations the number of isomers
C  will equal the number of atoms for C1 symmetry, so this is probably a bad idea.
C
