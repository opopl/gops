

!  Initialise PAIRDIST array for use in making an intial connection.
!  PAIRDIST should contain zero if the two minima are linked by a transition state.

! {{{
      IF (DIJINITT) THEN
        ! {{{
         IF (.NOT.INDEXCOSTFUNCTION) THEN
            DO J1=1,NMIN
               READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               PAIRDIST(J1*(J1+1)/2)=0.0D0
               DO J2=J1+1,NMIN
                  READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
                  CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                             RMAT,.FALSE.)
                  IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,
     &                                                     DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
                  PAIRDIST(J2*(J2-1)/2+J1)=DISTANCE
!                 PRINT '(A,3I6,G20.10)','J1,J2,INDEX,DISTANCE=',J1,J2,J2*(J2-1)/2+J1,DISTANCE
               ENDDO
               PRINT '(A,I8)','setup> Finished pair distance calculation for minimum ',J1
               CALL FLUSH(6,ISTAT)
            ENDDO
         ELSE
            DO J1=1,NMIN
               DO J2=J1+1,NMIN
                  PAIRDIST(J2*(J2-1)/2+J1)=1.0D0 ! this is not actually used !
               ENDDO
            ENDDO
         ENDIF
         DO J1=1,NPAIRDONE
            PAIRDIST(MAX(PAIR1(J1),PAIR2(J1))*(MAX(PAIR1(J1),PAIR2(J1))-1)/2+MIN(PAIR1(J1),PAIR2(J1)))=HUGE(1.0D0)
         ENDDO
         DO J1=1,NTS
! JMC n.b. don't apply the nconnmin criteria at this point, hence the huge(1) 's in place of NCONN() for the plus and minus minima.
            CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),HUGE(1),HUGE(1), 
     &                   PLUS(J1),MINUS(J1),.TRUE.,CUT_UNDERFLOW,DEADTS)
            IF (.NOT. DEADTS) THEN
               J2=MAX(PLUS(J1),MINUS(J1))
               J3=MIN(PLUS(J1),MINUS(J1))
               PAIRDIST(J2*(J2-1)/2+J3)=0.0D0
            ENDIF
         ENDDO
!
! Set the pairs for which connections have already been tried to infinite distance,
! so they are not tried again. Don;t overwrite zero distance settings for connections
! that have actually been found!
!
! }}}
      ENDIF
! }}}

