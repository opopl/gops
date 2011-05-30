
! Read data for minima from min.data.info style file. Multiple minima are allowed.

!op226 IF (READMINT) THEN {{{
      IF (READMINT) THEN
         CALL INQF(MINNAME,YESNO)
         IF (YESNO) THEN
            CALL OPENF(1,'O',MINNAME)
            DO 
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  READ(1,*,END=130) NEWEMIN,NEWFVIBMIN,NEWHORDERMIN,NEWIXMIN,NEWIYMIN,NEWIZMIN,NEWMINCURVE,NEWMINFRQ2
               ELSE
                  READ(1,*,END=130) NEWEMIN,NEWFVIBMIN,NEWHORDERMIN,NEWIXMIN,NEWIYMIN,NEWIZMIN
               ENDIF
               READ(1,*) (LOCALPOINTS(J2),J2=1,3*NATOMS)
!
! Must check it is not an old minimum!
!
               DO J2=1,NMIN
                  DISTANCE=1.0D100
                  IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
                     READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
                     CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, 
     &                                DIST2,RIGIDBODY,RMAT,.FALSE.)
                  ENDIF
      
                  IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
                     PRINT '(A,I6)','setup> minimum is database minimum ',J2
                     IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-4) THEN
                        WRITE(*,'(A,F15.5,A,F15.5)') 'setup> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
                     ENDIF
                     IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
                        WRITE(*,'(A,I6,A,I6)') 'setup> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
                        NEWHORDERMIN=MAX(NEWHORDERMIN,HORDERMIN(J2))
                        WRITE(*,'(A,I6)') 'setup> using maximum value: ',NEWHORDERMIN
                     ENDIF
                     GOTO 140
                  ENDIF
               ENDDO
               NMIN=NMIN+1
               IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
               EMIN(NMIN)=NEWEMIN
               FVIBMIN(NMIN)=NEWFVIBMIN
               HORDERMIN(NMIN)=NEWHORDERMIN
               IXMIN(NMIN)=NEWIXMIN
               IYMIN(NMIN)=NEWIYMIN
               IZMIN(NMIN)=NEWIZMIN
               IF (DUMMYTST.AND.LOWESTFRQT) MINCURVE(NMIN)=NEWMINCURVE
               IF (DUMMYTST.AND.LOWESTFRQT) MINFRQ2(NMIN)=NEWMINFRQ2
               WRITE(*,'(A,I6,A)') 'setup> new minimum ',NMIN,
     &                ' writing parameters to file min.data and points to points.min'
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  WRITE(UMINDATA,'(2F20.10,I6,5F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), 
     &                                             IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN), MINCURVE(NMIN),MINFRQ2(NMIN)
               ELSE
                  WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), 
     &                                             IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
               ENDIF
               CALL FLUSH(UMINDATA,ISTAT)
               WRITE(UMIN,REC=NMIN) (LOCALPOINTS(J2),J2=1,3*NATOMS)
140            CONTINUE
            ENDDO
130         CLOSE(1)
         ELSE
            PRINT '(A)','setup> ERROR - no file ',TRIM(ADJUSTL(MINNAME))
         ENDIF
         STOP
      ENDIF
      !op226 }}}

