C
C  Calculated Relative Qs between minima.
C
      ELSE IF (WORD.EQ.'AMHRELQ') THEN
         AMHRELQT=.TRUE.
         CALL READI(QRELONE)
         CALL READI(QRELTWO)
          PRINT '(A,I6,A,I6)','keywords> AMHRELQ min 1 ',QRELONE,' min 2',QRELTWO
         IF (NITEMS.GT.3) THEN
          PRINT '(A)','keywords> ERROR - AMHRELQ GT 2'
          STOP
         ENDIF


