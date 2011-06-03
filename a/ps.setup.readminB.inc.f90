      !op226 open min.B; read in locations of minima {{{

      INQUIRE(FILE='min.B',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(1,FILE='min.B',STATUS='OLD')
         READ(1,*) NMINB
         ALLOCATE(LOCATIONB(NMINB)) ! not deallocated
         READ(1,*) LOCATIONB(1:NMINB)
         CLOSE(1)
      ELSEIF (.NOT.(READMINT.OR.MERGEDBT)) THEN
         PRINT '(A)','setup> ERROR - no min.B file'
         STOP
      ELSE
         NMINB=0
      ENDIF
      IF (DEBUG) THEN
         PRINT '(A,I6,A)','setup> there are ',NMINB,' B minima at locations:'
         PRINT '(10I6)',LOCATIONB(1:NMINB)
      ENDIF
      IF (PRINTT) WRITE(*,'(A,I6,A)') 'setup> locations read for ',NMINB,' min of type B'

      !op226 }}}

