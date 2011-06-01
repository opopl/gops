
C  Calculated  Q_contact  between minima and native .
C
      ELSE IF (WORD.EQ.'AMHQCONT') THEN
         AMHQCONTT=.TRUE.
         CALL READI(WHICHMIN)
         CALL READF(QCONTCUT)
          PRINT '(A,I6)','keywords> Calculate AMH Q CONTACT ',WHICHMIN
          PRINT '(A,G3.2)','keywords> Within a distance cutoff ',QCONTCUT
         IF (NITEMS.GT.3) THEN
          PRINT '(A)','keywords> ERROR - AMH Q CONTACT '
          STOP
         ENDIF
C
