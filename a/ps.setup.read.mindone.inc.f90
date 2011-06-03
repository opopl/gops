
!  Read in the minima previously searched in UNTRAP runs.
! {{{
      ALLOCATE(MINDONE(MAXDONE))
      CALL INQF('min.done',YESNO)
      NMINDONE=0
      IF (YESNO) THEN
         CALL OPENF(1,'O','min.done')
         DO 
            NMINDONE=NMINDONE+1
            IF (NMINDONE.GT.MAXDONE) CALL DONEDOUBLE
            READ(1,*,END=121) MINDONE(NMINDONE)
            IF (DEBUG) PRINT '(A,I8,A,2I8)','setup > previously searched minimum ',
     &                                      NMINDONE,' is ',MINDONE(NMINDONE)
         ENDDO
121      CLOSE(1)
         NMINDONE=NMINDONE-1
         PRINT '(A,I8,A)','setup> ',NMINDONE,' minima already searched read from min.done'
      ENDIF
! }}}

