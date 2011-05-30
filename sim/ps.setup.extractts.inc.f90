     IF (EXTRACTTST) THEN
         OPEN(UNIT=1,FILE='extractedts', STATUS='UNKNOWN')
         IF (WHICHTS.EQ.-123) THEN
           ! {{{
            NDUMMY=0
            PRINT '(A)', 'setup> rewriting binary points.ts file from extractedts file'
            CALL MYSYSTEM(STATUS,DEBUG,'cp points.ts points.ts.save')
            DO
               NDUMMY=NDUMMY+1
               READ(1,*,END=778) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(UTS,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            ENDDO
778         NDUMMY=NDUMMY-1
            PRINT '(A,I6)','setup> number of ts extracted=',NDUMMY
            ! }}}
         ELSEIF (WHICHTS.LE.0) THEN
         ! {{{
            NDUMMY=0
            PRINT '(A)', 'setup> extracting all ts '
            DO
               NDUMMY=NDUMMY+1
               READ(UTS,REC=NDUMMY,ERR=878) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
            ENDDO
878         NDUMMY=NDUMMY-1
            PRINT '(A,I6)','setup> number of ts extracted=',NDUMMY
            ! }}}
         ELSE
           ! {{{
            PRINT '(A,I6)', 'setup> extracting ts ',WHICHTS
            READ(UTS,REC=WHICHTS) (LOCALPOINTS(J2),J2=1,3*NATOMS)

!           IF (AMHT) THEN
!              CALL AMHDUMP(LOCALPOINTS,'amhts.pdb')
!           ELSE
               WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
               IF (CHARMMT) CALL CHARMMDUMP(LOCALPOINTS,'extractedts.crd')
!          ENDIF
               ! }}}
         ENDIF
         CLOSE(1)
         STOP
      ENDIF


