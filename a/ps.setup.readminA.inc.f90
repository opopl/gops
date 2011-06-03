
      !op226 open min.A; read in locations of minima {{{
      
      INQUIRE(FILE='min.A',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(1,FILE='min.A',STATUS='OLD')
         READ(1,*) NMINA
         ALLOCATE(LOCATIONA(NMINA)) ! not deallocated
         READ(1,*) LOCATIONA(1:NMINA)
         CLOSE(1)
      ELSEIF (.NOT.(READMINT.OR.MERGEDBT)) THEN
         PRINT '(A)','setup> ERROR - no min.A file'
         STOP
      ELSE
         NMINA=0
      ENDIF
      IF (DEBUG) THEN
         PRINT '(A,I6,A)','setup> there are ',NMINA,' A minima at locations:'
         PRINT '(10I6)',LOCATIONA(1:NMINA)
      ENDIF
      IF (PRINTT) WRITE(*,'(A,I6,A)') 'setup> locations read for ',NMINA,' min of type A'

      !op226 }}}


