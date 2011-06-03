      !op226 IF (STARTFROMPATH) THEN {{{

      IF (STARTFROMPATH) THEN
!
!  Get all the necessary information about the A and B minima from the <PATHNAME> file.
!  The assumption is that we have just one A minimum and one B in this case.
!  Use GETNEWPATH or GETALLPATHS to do the bookkeeping.
!  Detect presence of existing min.data file and STOP if found to prevent overwriting.
!
         INQUIRE(FILE='min.data',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file min.data already exists. Will not overwrite.'
            STOP
         ENDIF
         CALL MYSYSTEM(STATUS,DEBUG,'cp ' // TRIM(ADJUSTL(PATHNAME)) // ' path.info')
         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='NEW')
         INQUIRE(FILE='ts.data',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file ts.data already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='NEW')
         IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
         IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
         IF (STARTTRIPLES) THEN
            CALL GETALLPATHS
         ELSE
            CALL GETNEWPATH(0,0)
         ENDIF
         INQUIRE(FILE='min.A',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file min.A already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(1,FILE='min.A',STATUS='NEW') 
         WRITE(1,'(I6)') 1
         WRITE(1,'(I6)') STARTMINA
         CLOSE(1)
         INQUIRE(FILE='min.B',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file min.B already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(1,FILE='min.B',STATUS='NEW') 
         WRITE(1,'(I6)') 1
         WRITE(1,'(I6)') STARTMINB
         CLOSE(1)
         CLOSE(UMINDATA)
         CLOSE(UTSDATA)
         PRINT '(A,I5,A,I5,A)','setup> The unique A and B minima are ',STARTMINA,' and ',STARTMINB,' respectively'
      ENDIF
      !op226 }}}

