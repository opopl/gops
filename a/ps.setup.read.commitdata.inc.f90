      IF (NPFOLD.GT.0) THEN
         INQUIRE(FILE='commit.data',EXIST=YESNO)
         GPFOLD(1:NMIN)=0.0D0
         IF (YESNO) THEN
            PRINT '(A)','setup> Reading initial committor probabilities read from commit.data'
            OPEN(UNIT=1,FILE='commit.data',STATUS='OLD')
            J2=0
            DO J1=1,NMIN
               READ(1,*,END=110) GPFOLD(J1)
               J2=J2+1
            ENDDO 
110         CLOSE(1)
            IF (J2.LT.NMIN) PRINT '(A)','setup> WARNING - end of file commit.data, remaining probabilities set to zero'
         ELSE
            IF (DIRECTION.EQ.'AB') THEN ! GPFOLD is PFA
               DO J1=1,NMINA
                  GPFOLD(LOCATIONA(J1))=1.0D0
               ENDDO
            ELSE ! GPFOLD is PFB
               DO J1=1,NMINB
                  GPFOLD(LOCATIONB(J1))=1.0D0
               ENDDO
            ENDIF
            PRINT '(A)','setup> Initial committor probabilities set to 0 or 1'
!           PRINT '(6G20.10)',GPFOLD(1:NMIN)
         ENDIF
      ENDIF


