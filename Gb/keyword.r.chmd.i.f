!op226> read in chmd.par{{{
! th368: 20-10-2009 Read parameter file containing CHARMM DYNAmics 
! parameters if either CHARMM/MD or CHARMM/NEWRESTART_MD was
! requested terminate if file is not found

         IF(CHMDT .OR. ( CHRMMT .AND. NEWRESTART_MD)) THEN

           INQUIRE(FILE='chmd.par',EXIST=YESNO)

           IF (YESNO) THEN
              OPEN(99,file='chmd.par')
              READ(99,'(A)') CHMDPAR
              CLOSE(99)
           ELSE
              WRITE(*,*) 'keywords> File chmd.par has to be provided.'
              STOP
           ENDIF
         ENDIF
! end th368: 20-10-2009
!op226> }}}

