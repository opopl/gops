
! If USEPAIRST is true then read the sequence of minima from file USEPAIRSFILE
! USEPAIRSFILE must be formatted as a single Epath file
IF (USEPAIRST) THEN
        ! {{{
         call openf(1,'O',USEPAIRSFILE)
         NUSEPAIRS=0
         DO
            READ(1,*,END=111) NDUMMY, DUMMY, NDUMMY
            NUSEPAIRS=NUSEPAIRS+1
         ENDDO
111      REWIND(1)
         PRINT '(A,A,A,I8,A,I8,A)','setup> Number of lines in file ',TRIM(ADJUSTL(USEPAIRSFILE)),' is ',NUSEPAIRS,' i.e. ',
     &           (NUSEPAIRS+1)/2,' minima'
         NUSEPAIRS=(NUSEPAIRS+1)/2
         ALLOCATE(USEPAIRSMIN(NUSEPAIRS))
         DO J1=1,NUSEPAIRS
            READ(1,*) NDUMMY, DUMMY, USEPAIRSMIN(J1)
            IF (J1.EQ.NUSEPAIRS) EXIT
            READ(1,*) NDUMMY, DUMMY, NDUMMY
         ENDDO
         CLOSE(1)
         PRINT '(A)','setup> Sequence of local minima:'
         PRINT '(15I8)',USEPAIRSMIN(1:NUSEPAIRS)
         ! }}}
      ENDIF

