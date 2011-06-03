
! Create a ts entry for DUMMYTS runs if there are minima that seem to lack such entries.

      IF (DUMMYTST) THEN
        ! {{{
         PRINT '(A)',' setup> shortest distances for local minima:'
         DO J1=1,NMIN
            IF (MINDISTMIN(J1).GT.HUGE(1.0D0)/1.0D1) THEN
               PRINT '(A,I8,A,G20.10)',' setup> in setup, minimum ',J1,' shortest distance unassigned'
               READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               DO J3=1,NMIN
                  IF (J3.EQ.J1) CYCLE
                  READ(UMIN,REC=J3) (LOCALPOINTS2(J2),J2=1,3*NATOMS)
                  CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                             RMAT,.FALSE.)
                  IF ((DISTANCE.LT.MINDISTMIN(J1)).OR.(DISTANCE.LT.MINDISTMIN(J3))) THEN
                     IF (DISTANCE.LT.MINDISTMIN(J1)) MINDISTMIN(J1)=DISTANCE
                     IF (DISTANCE.LT.MINDISTMIN(J3)) MINDISTMIN(J3)=DISTANCE
!
!  Must create an entry in ts.data in this case.
!  ETS,FVIBTS,HORDERTS,PLUS,MINUS,IXTS,IYTS,IZTS
!
                     IF (IMFRQT) THEN
                        PRINT '()', "setup> ERROR: can''t guess negative eigenvalue - don''t use DUMMYTS and IMFRQ"
                        STOP
                     ENDIF

                     WRITE(UTSDATA,'(2F20.10,3I10,3F20.10)') TSEGUESS(EMIN(J1),EMIN(J3),MINCURVE(J1),MINCURVE(J3),DISTANCE),
     &                           TSFVIBGUESS(EMIN(J1),EMIN(J3),FVIBMIN(J1),FVIBMIN(J3),MINFRQ2(J1),MINFRQ2(J3),NATOMS),
     &                           1,J3,J1,1.0D0,1.0D0,1.0D0
                     CALL FLUSH(UTSDATA,ISTAT)
                     NTS=NTS+1
                     IF (NTS.GT.MAXTS) CALL TSDOUBLE
                     ETS(NTS)=TSEGUESS(EMIN(J1),EMIN(J3),MINCURVE(J1),MINCURVE(J3),DISTANCE)
                     FVIBTS(NTS)=TSFVIBGUESS(EMIN(J1),EMIN(J3),FVIBMIN(J1),FVIBMIN(J3),MINFRQ2(J1),MINFRQ2(J3),NATOMS)
                     HORDERTS(NTS)=1
                     IXTS(NTS)=1.0D0
                     IYTS(NTS)=1.0D0
                     IZTS(NTS)=1.0D0
                     PLUS(NTS)=J3
                     MINUS(NTS)=J1
                     IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
                     IF (DEBUG) WRITE(*,'(A,I6,A)') 'setup> dummy ts ',NTS,' writing parameters to file ts.data'
!
!  Update ts pointers.
!
                     POINTERP(NTS)=-1
                     POINTERM(NTS)=-1
                     IF (TOPPOINTER(PLUS(NTS)).GT.0) POINTERP(NTS)=TOPPOINTER(PLUS(NTS))
                     IF (TOPPOINTER(MINUS(NTS)).GT.0) POINTERM(NTS)=TOPPOINTER(MINUS(NTS))
                     TOPPOINTER(PLUS(NTS))=NTS
                     TOPPOINTER(MINUS(NTS))=NTS
                  ENDIF
               ENDDO
            ENDIF
            PRINT '(A,I8,A,G20.10)',' setup> shortest distance for minimum ',J1,' is ',MINDISTMIN(J1)
         ENDDO
         ! }}}
      ENDIF

