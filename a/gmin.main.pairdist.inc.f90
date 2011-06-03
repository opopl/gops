
IF (PAIRDISTT) THEN
         PAIRDIST_FH=3000+MYNODE
         CALL OPENF(PAIRDIST_FH,">>","pairdists.dat")
         WRITE(PAIRDIST_FH,'(A10)',ADVANCE="NO") "Quench  "
         DO J1=1,NPAIRS
            WRITE(ATOM1,*) PAIRDIST(J1,1)
            WRITE(ATOM2,*) PAIRDIST(J1,2)
            WRITE(ATOMPAIR,*) TRIM(ADJUSTL(ATOM1))//"-"//TRIM(ADJUSTL(ATOM2))
            WRITE(PAIRDIST_FH,'(A10)',ADVANCE="NO") TRIM(ADJUSTL(ATOMPAIR))//"  " 
         ENDDO
         WRITE(LFH,'(A)') ""
      ENDIF

