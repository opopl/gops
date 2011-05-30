
      DO J1=1,NMIN
         TOPPOINTER(J1)=-1
      ENDDO

      CALL INQF('ts.data',YESNO)
      IF (YESNO.AND.DIJINITSTARTT) THEN
         IF (.NOT.DUMMYTST) THEN
            CALL MYSYSTEM(STATUS,DEBUG,'mv ts.data ts.data.save')
            PRINT '(A)','setup> Removing old ts.data file'
            YESNO=.FALSE.
         ENDIF
      ENDIF

      IF (YESNO) THEN
         CALL OPENF(UTSDATA,'ts.data','O')
         ! {{{
         J1=0
         DO 
            J1=J1+1
            IF (J1.GT.MAXTS) CALL TSDOUBLE
            IF (IMFRQT) THEN
               READ(UTSDATA,*,END=40) ETS(J1),FVIBTS(J1),HORDERTS(J1),PLUS(J1),MINUS(J1),IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
            ELSE
               READ(UTSDATA,*,END=40) ETS(J1),FVIBTS(J1),HORDERTS(J1),PLUS(J1),MINUS(J1),IXTS(J1),IYTS(J1),IZTS(J1)
            ENDIF
            IF ((PLUS(J1).GT.NMIN).OR.(MINUS(J1).GT.NMIN)) THEN
               PRINT '(A,I6,A)','setup> ERROR - minima specified for ts ',J1,' lie beyond those specified in min.data'
               PRINT '(A,2I6)','setup> plus and minus minima are ',PLUS(J1),MINUS(J1)
               STOP
            ENDIF

            IF (DUMMYTST.AND.(.NOT.NOPOINTS).AND.(NATTEMPT.GT.0)) THEN
               READ(UMIN,REC=PLUS(J1)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               READ(UMIN,REC=MINUS(J1)) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
                  CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                             RMAT,.FALSE.)
                  IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, 
     &                                                     DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               IF (DISTANCE.LT.MINDISTMIN(PLUS(J1))) MINDISTMIN(PLUS(J1))=DISTANCE
               IF (DISTANCE.LT.MINDISTMIN(MINUS(J1))) MINDISTMIN(MINUS(J1))=DISTANCE
            ENDIF
         ENDDO
40       CLOSE(UTSDATA) ! SAT need to reopen this file

         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION="APPEND",ACTION="READWRITE",FORM="FORMATTED") ! read used in Dijkstra
         NTS=J1-1
         IF (PRINTT) WRITE(*,'(A,I7,A)') 'setup> parameters read for ',NTS,' ts'

         IF (DIJKSTRAT .OR. KSHORTESTPATHST) THEN
           ! {{{
            INQUIRE(FILE='ts.attempts',EXIST=YESNO)
            TSATTEMPT(1:NTS)=0
            IF (YESNO) THEN
               PRINT '(A)','setup> Reading the number of searches for existing transition states from ts.attempts'
               OPEN(UNIT=1,FILE='ts.attempts',STATUS='UNKNOWN')
               J2=0
               DO J1=1,NTS
                  READ(1,'(I8)',END=51) TSATTEMPT(J1)
                  J2=J2+1
               ENDDO
51             CLOSE(1)
               IF (J2.LT.NTS) PRINT '(A)','setup> WARNING - end of file ts.attempts, remaining attempts set to zero'
            ENDIF
            ! }}}
         ENDIF

         DO J1=1,NTS
            POINTERP(J1)=-1
            POINTERM(J1)=-1
         ENDDO
         DO J1=NTS,1,-1
            IF (J1.GT.TOPPOINTER(PLUS(J1)))  TOPPOINTER(PLUS(J1))=J1
            IF (J1.GT.TOPPOINTER(MINUS(J1))) TOPPOINTER(MINUS(J1))=J1
            DO J2=J1-1,1,-1
               IF (PLUS(J2).EQ.PLUS(J1)) THEN
                  POINTERP(J1)=J2
                  GOTO 41
               ELSE IF (MINUS(J2).EQ.PLUS(J1)) THEN
                  POINTERP(J1)=J2
                  GOTO 41
               ENDIF
            ENDDO
41          CONTINUE
            DO J2=J1-1,1,-1
               IF (PLUS(J2).EQ.MINUS(J1)) THEN
                  POINTERM(J1)=J2
                  GOTO 42
               ELSE IF (MINUS(J2).EQ.MINUS(J1)) THEN
                  POINTERM(J1)=J2
                  GOTO 42
               ENDIF
            ENDDO
42          CONTINUE
         ENDDO
   
         IF ((.NOT.NOPOINTS).AND.(NATTEMPT.GT.0).AND.(.NOT.DUMMYTST)) THEN
           ! {{{
            IF (AMHT) THEN
              WRITE(*,*) 'setup> AVOIDING MOMENT OF INITERIA CALC FOR AMH'
            ELSE
               DO J1=1,NTS
                  READ(UTS,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
                  CALL INERTIAWRAPPER(LOCALPOINTS,NATOMS,angleAxis,IXM,IYM,IZM)
!                 IF (DEBUG) WRITE(*,'(A,I7,2F17.7,3I6,3F15.5)') 'setup> ',J1,ETS(J1),FVIBTS(J1),HORDERTS(J1),
!    1                                                            PLUS(J1),MINUS(J1),IXM,IYM,IZM
                  IF ((ABS(IXM-IXTS(J1)).GT.IDIFFTOL).OR.&
                     (ABS(IYM-IYTS(J1)).GT.IDIFFTOL).OR.&
                     (ABS(IZM-IZTS(J1)).GT.IDIFFTOL)) THEN
                     WRITE(*,'(A,I10)') 'setup> WARNING - principal moments of inertia do not agree with input for ts ',J1
                     WRITE(*,'(A)') 'values in ts.data:'
                     WRITE(*,'(3F20.10)') IXTS(J1),IYTS(J1),IZTS(J1)
                     WRITE(*,'(A)') 'values recalculated from points.ts:'
                     WRITE(*,'(3F20.10)') IXM,IYM,IZM
!                    STOP  
                  ENDIF
               ENDDO
            ENDIF
            ! }}}
         ENDIF
         ! }}}
      ELSE
        ! {{{
         WRITE(*,'(A)') 'setup> no transition states found'
         INQUIRE(FILE='ts.data',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file ts.data already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='NEW')
         NTS=0
         ! }}}
      ENDIF

