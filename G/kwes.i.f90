
      IF (END .OR. WORD .EQ. 'STOP') THEN

         IF (NPCOUNT.LT.NPAR) THEN
            DO J1=NPCOUNT+1,NPAR
               STEP(J1)=STEP(1)
               ASTEP(J1)=ASTEP(1)
               OSTEP(J1)=OSTEP(1)
               BLOCK(J1)=BLOCK(1)
            ENDDO
         ENDIF
        RETURN
      ENDIF


