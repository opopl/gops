C
C  Calculated Relative Contact Order
C
       ELSE IF (WORD.EQ.'AMH_RELCO') THEN
         AMH_RELCOT=.TRUE.
         CALL READI(WHICHMIN)
         CALL READF(RELCOCUT)
          PRINT '(A,I6)','keywords> Calculate AMH Relative Contact Order ',WHICHMIN
          PRINT '(A,G20.10)','keywords> Within a distance cutoff ',RELCOCUT

         IF (NITEMS.GT.3) THEN
          PRINT '(A)','keywords> ERROR - RELCO GT 2'
          STOP
         ENDIF


