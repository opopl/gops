
IF (EXTRACTMINT) THEN
         OPEN(UNIT=1,FILE='extractedmin', STATUS='UNKNOWN')
         IF (WHICHMIN.EQ.-123) THEN
           ! {{{
            NDUMMY=0
            PRINT '(A)', 'setup> rewriting binary points.min file from extractedmin file'
            CALL MYSYSTEM(STATUS,DEBUG,'cp points.min points.min.save')
            DO 
               NDUMMY=NDUMMY+1
               READ(1,*,END=777) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(UMIN,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            ENDDO
777         NDUMMY=NDUMMY-1
            PRINT '(A,I6)','setup> number of minima extracted=',NDUMMY
            ! }}}
         ELSEIF (WHICHMIN.LE.0) THEN
         ! {{{
            NDUMMY=0
            PRINT '(A)', 'setup> extracting all minima '
            DO 
               NDUMMY=NDUMMY+1
               READ(UMIN,REC=NDUMMY,ERR=877) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
            ENDDO
877         NDUMMY=NDUMMY-1
            PRINT '(A,I6)','setup> number of minima extracted=',NDUMMY
            ! }}}
         ELSE
           ! {{{
            PRINT '(A,I6)', 'setup> extracting minimum ',WHICHMIN
            READ(UMIN,REC=WHICHMIN) (LOCALPOINTS(J2),J2=1,3*NATOMS)

!           IF (AMHT) THEN
!              CALL AMHDUMP(LOCALPOINTS,'amhmin.pdb')
!  
!              IF (AMHQT)THEN
!                 CALL AMHQ(WHICHMIN)
!              ENDIF
!  
!              IF (AMHQCONTT)THEN
!                 CALL AMHQCONT(WHICHMIN,QCONTCUT)
!              ENDIF
!  
!              IF (AMHRMSDT)THEN
!                 CALL AMHRMSD(WHICHMIN)
!              ENDIF
!  
!              IF (AMHRELQT)THEN
!                 CALL AMHRELQ(QRELONE, QRELTWO)
!              ENDIF
!  
!              IF (AMH_RELCOT)THEN
!                 CALL AMH_RELCO(WHICHMIN, RELCOCUT)
!              ENDIF
!  
!              IF (AMHALLATOMMINT)THEN
!                 CALL AMHALLATOMMIN
!              ENDIF
!           ELSE
               WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
               IF (CHARMMT) CALL CHARMMDUMP(LOCALPOINTS,'extractedmin.crd')
!           ENDIF
        ! }}}
         ENDIF
         CLOSE(1)
         STOP
      ENDIF

