
!  In this case we try to connect two minima, whose coordinates are in
!  files odata.start and odata.finish for a DIJINITSTARTT run. We set up min.data as 
!  well as min.A and min.B, with entries 1 1 and 1 2.
!
!  DIJINITCONTT specifies a continuation of an initial path run.
!

      !op226 IF (DIJINITT.OR.DIJINITFLYT) THEN {{{
      IF (DIJINITT.OR.DIJINITFLYT) THEN
         INQUIRE(FILE='min.data',EXIST=YESNO)
         IF (DIJINITCONTT) THEN
            IF (.NOT.YESNO) THEN
               PRINT '(A)','setup> ERROR - min.data must exist for a DIJINIT continuation run'
               STOP
            ENDIF
         ELSEIF (DIJINITSTARTT) THEN
            IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'cp min.data min.data.save')
            PRINT '(A)','setup> Creating new min.data file'
            OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN')
            INQUIRE(FILE='pairs.data',EXIST=YESNO)
            IF (YESNO) THEN
               PRINT '(A)','setup> Moving old pairs.data file'
               CALL MYSYSTEM(STATUS,DEBUG,'mv pairs.data pairs.data.old')
            ENDIF
            INQUIRE(FILE='commit.data',EXIST=YESNO)
            IF (YESNO) THEN
               PRINT '(A)','setup> Moving old commit.data file'
               CALL MYSYSTEM(STATUS,DEBUG,'mv commit.data commit.data.old')
            ENDIF
            INQUIRE(FILE='min.data.info.1',EXIST=YESNO)
            IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'rm min.data.info.1')
            CALL MYSYSTEM(STATUS,DEBUG,'cp odata.start odata.1')
            IF (CHARMMT) CALL MYSYSTEM(STATUS,DEBUG,'cp start.crd input.crd.1')
            IF (AMHT) CALL MYSYSTEM(STATUS,DEBUG,'cp AMHstart start.1')
            CALL MYSYSTEM(STATUS,DEBUG,TRIM(ADJUSTL(EXEC))//' 1 > output.start')
            INQUIRE(FILE='min.data.info.1',EXIST=YESNO)
            IF (YESNO) THEN
               OPEN(UNIT=1,FILE='min.data.info.1',STATUS='OLD')
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  READ(1,*) EMIN(1),FVIBMIN(1),HORDERMIN(1),IXMIN(1),IYMIN(1),IZMIN(1),MINCURVE(1),MINFRQ2(1)
                  WRITE(UMINDATA,'(2F20.10,I6,5F20.10)') EMIN(1), FVIBMIN(1), HORDERMIN(1),IXMIN(1),
     &                                                      IYMIN(1),IZMIN(1),MINCURVE(1),MINFRQ2(1)
               ELSE
                  READ(1,*) EMIN(1),FVIBMIN(1),HORDERMIN(1),IXMIN(1),IYMIN(1),IZMIN(1)
                  WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(1), FVIBMIN(1), HORDERMIN(1),IXMIN(1), IYMIN(1), IZMIN(1)
               ENDIF
               READ(1,*) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(UMIN,REC=1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               CLOSE(1)
            ELSE
               PRINT *, 'setup> ERROR - no min.data.info.1 found - check OPTIM output in output.start'
               STOP        
            ENDIF
            INQUIRE(FILE='min.data.info.2',EXIST=YESNO)
            IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'rm min.data.info.2')
            CALL MYSYSTEM(STATUS,DEBUG,'cp odata.finish odata.2')
            IF (CHARMMT) CALL MYSYSTEM(STATUS,DEBUG,'cp finish.crd input.crd.2')
            IF (AMHT) CALL MYSYSTEM(STATUS,DEBUG,'cp AMHfinish start.2')
            CALL MYSYSTEM(STATUS,DEBUG,TRIM(ADJUSTL(EXEC))//' 2 > output.finish')
            INQUIRE(FILE='min.data.info.2',EXIST=YESNO)
            IF (YESNO) THEN
               OPEN(UNIT=1,FILE='min.data.info.2',STATUS='OLD')
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  READ(1,*) EMIN(2),FVIBMIN(2),HORDERMIN(2),IXMIN(2),IYMIN(2),IZMIN(2),MINCURVE(2),MINFRQ2(2)
                  WRITE(UMINDATA,'(2F20.10,I6,5F20.10)') EMIN(2), FVIBMIN(2), HORDERMIN(2),IXMIN(2), 
     &                                                   IYMIN(2),IZMIN(2),MINCURVE(2),MINFRQ2(2)
               ELSE
                  READ(1,*) EMIN(2),FVIBMIN(2),HORDERMIN(2),IXMIN(2),IYMIN(2),IZMIN(2)
                  WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(2), FVIBMIN(2), HORDERMIN(2),IXMIN(2), IYMIN(2), IZMIN(2)
               ENDIF
               CLOSE(UMINDATA)
               READ(1,*) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(UMIN,REC=2) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               CLOSE(1)
            ELSE
               PRINT *, 'setup> ERROR - no min.data.info.2 found - check OPTIM output in output.finish'
               STOP
            ENDIF
            INQUIRE(FILE='min.A',EXIST=YESNO)
            IF (DIJINITSTARTT) THEN
               IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'cp min.A min.A.save')
               PRINT '(A)','setup> Creating new min.A file'
            ELSEIF (DIJINITCONTT) THEN
               IF (.NOT.YESNO) THEN
                  PRINT '(A)','setup> ERROR - min.A must exist for a DIJINIT continuation run'
                  STOP
               ENDIF
            ENDIF
            OPEN(1,FILE='min.A',STATUS='UNKNOWN')
            WRITE(1,'(I6)') 1
            WRITE(1,'(I6)') 1
            CLOSE(1)
            INQUIRE(FILE='min.B',EXIST=YESNO)
            IF (DIJINITSTARTT) THEN
               IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'cp min.B min.B.save')
               PRINT '(A)','setup> Creating new min.B file'
            ELSEIF (DIJINITCONTT) THEN
               IF (.NOT.YESNO) THEN
                  PRINT '(A)','setup> ERROR - min.B must exist for a DIJINIT continuation run'
                  STOP
               ENDIF
            ENDIF
            OPEN(1,FILE='min.B',STATUS='UNKNOWN')
            WRITE(1,'(I6)') 1
            WRITE(1,'(I6)') 2
            CLOSE(1)
            PRINT '(A)','setup> initial OPTIM jobs run for odata.start and odata.finish'

            IF (DUMMYTST) THEN
               READ(UMIN,REC=1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               READ(UMIN,REC=2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                          RMAT,.FALSE.)
               IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, 
     &                                                  DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               MINDISTMIN(1)=DISTANCE
               MINDISTMIN(2)=DISTANCE
!
!  Must create an entry in ts.data in this case.
!  ETS,FVIBTS,HORDERTS,PLUS,MINUS,IXTS,IYTS,IZTS
!
               INQUIRE(FILE='ts.data',EXIST=YESNO)
               IF (YESNO) THEN
                  PRINT '(A)','setup> - file ts.data already exists. Copying to ts.data.save'
                  CALL MYSYSTEM(STATUS,DEBUG,'mv ts.data ts.data.save')
               ENDIF
               OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='NEW') 

               IF (IMFRQT) THEN
                  PRINT '()',"setup> ERROR: can''t guess negative eigenvalue - don''t use DUMMYTS and IMFRQ"
                  STOP
               ENDIF
               WRITE(UTSDATA,'(2F20.10,3I10,3F20.10)') TSEGUESS(EMIN(1),EMIN(2),MINCURVE(1),MINCURVE(2),DISTANCE), 
     &               TSFVIBGUESS(EMIN(1),EMIN(2),FVIBMIN(1),FVIBMIN(2),MINFRQ2(1),MINFRQ2(2),NATOMS),1,1,2,1.0D0,1.0D0,1.0D0

               CLOSE(UTSDATA)
            ENDIF
         ENDIF
      ENDIF
      !op226 }}}

