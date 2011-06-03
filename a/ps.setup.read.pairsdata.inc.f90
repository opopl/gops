
!  Read in the pairs of minima previously searched in pairs.data exists.

!{{{
      ALLOCATE(PAIR1(MAXPAIRS),PAIR2(MAXPAIRS))
      INQUIRE(FILE='pairs.data',EXIST=YESNO)
      NPAIRDONE=0
      IF (YESNO) THEN
         OPEN(UNIT=1,FILE='pairs.data',STATUS='OLD')
         DO 
            NPAIRDONE=NPAIRDONE+1
            IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
            READ(1,*,END=120) PAIR1(NPAIRDONE), PAIR2(NPAIRDONE)
            IF (DEBUG) PRINT '(A,I8,A,2I8)','setup > previously searched pair number ',
     &                                      NPAIRDONE,' is ',PAIR1(NPAIRDONE), PAIR2(NPAIRDONE)
            IF ((PAIR1(NPAIRDONE).GT.NMIN).OR.(PAIR2(NPAIRDONE).GT.NMIN)) THEN
               PRINT '(A)','setup> ERROR *** minima specified in pairs.data do not exist in min.data'
               STOP
            ENDIF
         ENDDO
120      CLOSE(1)
         NPAIRDONE=NPAIRDONE-1
         PRINT '(A,I8,A)','setup> ',NPAIRDONE,' pairs of minima already searched read from pairs.data'
      ENDIF
! }}}

