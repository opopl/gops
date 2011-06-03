!
!  Load the minima.
!
      !op226 open min.data; load the minima {{{

      CALL INQF('min.data',YESNO)
      IF (YESNO) THEN
         CALL OPENF(UMINDATA,'O','min.data')
         NMIN=0
         J1=0
         DO 
            J1=J1+1
            IF (J1.GT.MAXMIN) CALL MINDOUBLE
            IF (DUMMYTST.AND.LOWESTFRQT) THEN
               READ(UMINDATA,*,END=30) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1),MINCURVE(J1),MINFRQ2(J1)
            ELSE
               READ(UMINDATA,*,END=30) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1)
            ENDIF
            NMIN=NMIN+1
         ENDDO
      ELSEIF (.NOT.(READMINT.OR.MERGEDBT)) THEN
         PRINT '(A)','setup> ERROR - no min.data file'
         STOP
      ELSE
         NMIN=0
      ENDIF

      !op226 }}}
 
