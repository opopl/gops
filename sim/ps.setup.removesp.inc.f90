
!
!  Procedure to remove selected stationary points specified by min.remove and ts.remove.
!  First line of each file gives the numbers of structures to remove.
!
      IF (REMOVESP) THEN 
         ! {{{
         OPEN(UNIT=1,FILE='min.remove',STATUS='OLD')
         READ(1,*) NMINREMOVE
         ALLOCATE(MINREMOVE(NMINREMOVE),MINPREV(NMIN))
         READ(1,*) MINREMOVE(1:NMINREMOVE)
         CLOSE(1)
         PRINT '(A)','setup> removing the following minima:'
         PRINT '(10I8)',MINREMOVE(1:NMINREMOVE)
         OPEN(UNIT=2,FILE='min.data.removed',STATUS='UNKNOWN')
         OPEN(UNIT=4,FILE='points.min.removed',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
         NDUMMY=0
         MINPREV(1:NMIN)=0
         minloop: DO J1=1,NMIN
            DO J2=1,NMINREMOVE
               IF (MINREMOVE(J2).EQ.J1) THEN
                  PRINT '(A,I8)','setup> removing minimum ',J1
                  CYCLE minloop
               ENDIF
            ENDDO
            PRINT '(A,I8)','setup> not removing minimum ',J1
            NDUMMY=NDUMMY+1
            MINPREV(J1)=NDUMMY
            WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J1), FVIBMIN(J1), HORDERMIN(J1), IXMIN(J1), IYMIN(J1), IZMIN(J1)
            READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            WRITE(4,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ENDDO minloop
         CLOSE(2)
         CLOSE(4)
!
! rewrite min.A and min.B in min.A.removed and min.B.removed since numbering may change.
!
         NDUMMY=0
         Aloop: DO J1=1,NMINA
            DO J2=1,NMINREMOVE
               IF (MINREMOVE(J2).EQ.LOCATIONA(J1)) THEN
                  PRINT '(A,I8)','setup> removing A minimum ',LOCATIONA(J1)
                  CYCLE Aloop
               ENDIF
            ENDDO
            NDUMMY=NDUMMY+1
         ENDDO Aloop
         OPEN(UNIT=2,FILE='min.A.removed',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all A minima removed'
            STOP
         ENDIF
         DO J1=1,NMINA
            IF (MINPREV(LOCATIONA(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONA(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         NDUMMY=0
         Bloop: DO J1=1,NMINB
            DO J2=1,NMINREMOVE
               IF (MINREMOVE(J2).EQ.LOCATIONB(J1)) THEN
                  PRINT '(A,I8)','setup> removing B minimum ',LOCATIONB(J1)
                  CYCLE Bloop
               ENDIF
            ENDDO
            NDUMMY=NDUMMY+1
         ENDDO Bloop
         OPEN(UNIT=2,FILE='min.B.removed',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all B minima removed'
            STOP
         ENDIF
         DO J1=1,NMINB
            IF (MINPREV(LOCATIONB(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONB(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         OPEN(UNIT=1,FILE='ts.remove',STATUS='OLD')
         READ(1,*) NTSREMOVE
         IF (NTSREMOVE.LE.0) GOTO 444
         ALLOCATE(TSREMOVE(NTSREMOVE))
         READ(1,*) TSREMOVE(1:NTSREMOVE)
         CLOSE(1)
         PRINT '(A)','setup> removing the following transition states:'
         PRINT '(10I8)',TSREMOVE(1:NTSREMOVE)
444      CONTINUE
         OPEN(UNIT=3,FILE='ts.data.removed',STATUS='UNKNOWN')
         OPEN(UNIT=5,FILE='points.ts.removed',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
         NDUMMY=0
         tsloop: DO J1=1,NTS
            DO J2=1,NTSREMOVE
               IF (TSREMOVE(J2).EQ.J1) CYCLE tsloop
            ENDDO
            IF (MINPREV(PLUS(J1)).EQ.0) THEN
               PRINT '(A,I8,A,I8,A)','setup> possible ERROR - transition state ',J1,' links minimum ',PLUS(J1), 
     &                               ' which has been removed - removing this ts'
               CYCLE tsloop
            ENDIF
            IF (MINPREV(MINUS(J1)).EQ.0) THEN
               PRINT '(A,I8,A,I8,A)','setup> possible ERROR - transition state ',J1,' links minimum ',MINUS(J1), 
     &                               ' which has been removed - removing this ts'
               CYCLE tsloop
            ENDIF
            NDUMMY=NDUMMY+1
            IF (IMFRQT) THEN
               WRITE(3,'(2F20.10,3I10,4F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
            ELSE
               WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1)
            ENDIF
            READ(UTS,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            WRITE(5,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ENDDO tsloop
         CLOSE(3); CLOSE(5)
         STOP
         ! }}}
      ENDIF

