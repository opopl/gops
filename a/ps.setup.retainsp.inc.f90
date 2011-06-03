
!  Procedure to retain selected stationary points specified by min.retain.
!  First line of each file gives the numbers of structures to retain.
!  All ts linking minima in the retain list are themselves retained.

      IF (RETAINSP) THEN
        ! {{{
         OPEN(UNIT=1,FILE='min.retain',STATUS='OLD')
         READ(1,*) NMINRETAIN
         ALLOCATE(MINRETAIN(NMINRETAIN),MINPREV(NMIN))
         READ(1,*) MINRETAIN(1:NMINRETAIN)
         CLOSE(1)
         PRINT '(A)','setup> retaining the following minima:'
         PRINT '(10I8)',MINRETAIN(1:NMINRETAIN)
         OPEN(UNIT=2,FILE='min.data.retained',STATUS='UNKNOWN')
         OPEN(UNIT=4,FILE='points.min.retained',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
         NDUMMY=0
         MINPREV(1:NMIN)=0
         ! minloop2 {{{
         minloop2: DO J1=1,NMIN
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.J1) THEN
                  PRINT '(A,I8)','setup> retaining minimum ',J1
                  NDUMMY=NDUMMY+1
                  MINPREV(J1)=NDUMMY
                  WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J1), FVIBMIN(J1), HORDERMIN(J1), IXMIN(J1), IYMIN(J1), IZMIN(J1)
                  READ(UMIN,REC=J1) (LOCALPOINTS(J3),J3=1,3*NATOMS)
                  WRITE(4,REC=NDUMMY) (LOCALPOINTS(J3),J3=1,3*NATOMS)
                  CYCLE minloop2
               ENDIF
            ENDDO
            PRINT '(A,I8)','setup> removing minimum ',J1
         ENDDO minloop2
         ! }}}
         CLOSE(2)
         CLOSE(4)
!
! rewrite min.A and min.B in min.A.retained and min.B.retained since numbering may change.
!
         NDUMMY=0
         Aloop2: DO J1=1,NMINA
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.LOCATIONA(J1)) THEN
                  NDUMMY=NDUMMY+1
                  PRINT '(A,I8)','setup> retaining A minimum ',LOCATIONA(J1)
                  CYCLE Aloop2
               ENDIF
            ENDDO
         ENDDO Aloop2
         OPEN(UNIT=2,FILE='min.A.retained',STATUS='UNKNOWN')
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
         ! Bloop2 {{{
         Bloop2: DO J1=1,NMINB
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.LOCATIONB(J1)) THEN
                  NDUMMY=NDUMMY+1
                  PRINT '(A,I8)','setup> retaining B minimum ',LOCATIONB(J1)
                  CYCLE Bloop2
               ENDIF
            ENDDO
         ENDDO Bloop2
         ! }}}
         OPEN(UNIT=2,FILE='min.B.retained',STATUS='UNKNOWN')
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

         OPEN(UNIT=3,FILE='ts.data.retained',STATUS='UNKNOWN')
         OPEN(UNIT=5,FILE='points.ts.retained',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
         NDUMMY=0
         ! tsloop2 {{{
         tsloop2: DO J1=1,NTS
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.PLUS(J1)) THEN
                  DO J3=1,NMINRETAIN
                     IF (MINRETAIN(J3).EQ.MINUS(J1)) THEN
                        NDUMMY=NDUMMY+1
                        PRINT '(A,I8,A,2I8)','setup> retaining ts ',J1,' connected to retained minima ',PLUS(J1),MINUS(J1)
                        IF (IMFRQT) THEN
                           WRITE(3,'(2F20.10,3I10,4F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
                        ELSE
                           WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1)
                        ENDIF
                        READ(UTS,REC=J1) (LOCALPOINTS(J4),J4=1,3*NATOMS)
                        WRITE(5,REC=NDUMMY) (LOCALPOINTS(J4),J4=1,3*NATOMS)
                        CYCLE tsloop2
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            PRINT '(A,I8)','setup> removing ts ',J1
         ENDDO tsloop2
         ! }}}
         CLOSE(3); CLOSE(5)
         PRINT '(A,I8,A)','setup> ',NDUMMY,' transition states retained'
         STOP
         ! }}}
      ENDIF

