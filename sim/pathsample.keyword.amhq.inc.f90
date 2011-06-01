C
C  Calculated  Q between minima and native .
C
      ELSE IF (WORD.EQ.'AMHQ') THEN
         AMHQT=.TRUE.
         CALL READI(WHICHMIN)
          PRINT '(A,I6,A,I6)','keywords> Calculate AMH Q  ',WHICHMIN
         IF (NITEMS.GT.2) THEN
          PRINT '(A)','keywords> ERROR - AMHQ '
          STOP
         ENDIF
C
