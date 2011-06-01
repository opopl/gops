C  Calculated  RMSD between minima and native .
C
      ELSE IF (WORD.EQ.'AMHRMSD') THEN
         AMHRMSDT=.TRUE.
         CALL READI(WHICHMIN)
          PRINT '(A,I6,A,I6)','keywords> Calculate AMH RMSD  ',WHICHMIN
         IF (NITEMS.GT.2) THEN
          PRINT '(A)','keywords> ERROR - AMHRMSD '
          STOP
         ENDIF

