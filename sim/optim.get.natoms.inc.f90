
200   CONTINUE

      IF (CASTEP) THEN
         FNAME=SYS(1:LSYS) // '.cell'
         WRITE(*,'(A,A)') ' getparams> Counting atoms in file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
         coordsloop: DO
            READ(7,'(A21)') DSTRING1
            CALL UPPERCASE(DSTRING1)
!           WRITE(*,'(A21)') DSTRING1
            IF (DSTRING1(1:16).EQ.'%BLOCK POSITIONS') THEN
               DO 
                  READ(7,'(A21)') DSTRING1
                  CALL UPPERCASE(DSTRING1)
!                 WRITE(*,*) DSTRING1
                  IF (DSTRING1(1:19).EQ.'%ENDBLOCK POSITIONS') EXIT coordsloop
                  NATOMS=NATOMS+1
               ENDDO
            ENDIF
         ENDDO coordsloop
!        WRITE(*,'(A,I4,A)') 'getparams> CASTEP run for ',NATOMS,' atoms'
         CLOSE(7)
      ELSEIF (ONETEP) THEN
         FNAME=SYS(1:LSYS) // '.dat'
         WRITE(*,'(A,A)') 'getparams> Counting atoms in file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
         coordsloop2: DO
            READ(7,'(A21)') DSTRING1
            CALL UPPERCASE(DSTRING1)
!           WRITE(*,'(A21)') DSTRING1
            IF (DSTRING1(1:16).EQ.'%BLOCK POSITIONS') THEN
               DO 
                  READ(7,'(A21)') DSTRING1
                  CALL UPPERCASE(DSTRING1)
!                 WRITE(*,*) DSTRING1
                  IF (DSTRING1(1:19).EQ.'%ENDBLOCK POSITIONS') EXIT coordsloop2
                  NATOMS=NATOMS+1
               ENDDO
            ENDIF
         ENDDO coordsloop2
         WRITE(*,'(A,I4,A)') 'getparams> ONETEP run for ',NATOMS,' atoms'
         CLOSE(7)
      ELSEIF (CP2K) THEN 
         FNAME=SYS(1:LSYS) // '.inp'
         WRITE(*,'(A,A)') ' getparams> Counting atoms in file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
         coordsloop3: DO
            READ(7,'(A21)') DSTRING1
            DSTRING1=ADJUSTL(DSTRING1)
            CALL UPPERCASE(DSTRING1)
!           WRITE(*,'(A21)') DSTRING1
            IF (DSTRING1(1:6).EQ.'&COORD') THEN
               DO 
                  READ(7,'(A21)') DSTRING1
                  DSTRING1=ADJUSTL(DSTRING1)
                  CALL UPPERCASE(DSTRING1)
!                 WRITE(*,*) DSTRING1
                  IF (DSTRING1(1:10).EQ.'&END COORD') EXIT coordsloop3
                  IF (DSTRING1(1:6).EQ.'SCALED') NATOMS=NATOMS-1
                  IF (DSTRING1(1:4).EQ.'UNIT') NATOMS=NATOMS-1
                  NATOMS=NATOMS+1
               ENDDO
            ENDIF
         ENDDO coordsloop3
         WRITE(*,'(A,I4,A)') ' getparams> CP2K run for ',NATOMS,' atoms'
         CLOSE(7) 
      ELSE IF (CPMD) THEN
         FNAME=SYS(1:LSYS)
         WRITE(*,'(A,A)') 'getparams> Reading coordinates from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
11       READ(7,'(A)') FNAME
         IF (FNAME(1:6).EQ.'&ATOMS') THEN
            J1=0
12          READ(7,'(A)') TSTRING
            IF (TSTRING(1:1).EQ.'*') THEN
               J1=J1+1
               READ(7,'(A)') FNAME
               READ(7,*) NTYPE(J1)
               DO J2=1,NTYPE(J1)
                  READ(7,*) DUMMY
               ENDDO
               NATOMS=NATOMS+NTYPE(J1)
               GOTO 12
            ELSE IF (TSTRING(1:1).EQ.' ') THEN
               GOTO 12
            ELSE IF (TSTRING(1:4).EQ.'&END') THEN
               GOTO 13
            ENDIF
         ELSE
            GOTO 11
         ENDIF
13       CONTINUE
      ELSE IF (VARIABLES) THEN
         NATOMS=0
110      CALL INPUT(END)
         IF (.NOT.END) THEN
            CALL READF(DUMMY)
            NATOMS=NATOMS+1
            GOTO 110
         ENDIF
      ELSE IF (RINGPOLYMERT) THEN
         NATOMS=0
111      CALL INPUT(END)
         IF (.NOT.END) THEN
            CALL READF(DUMMY)
            NATOMS=NATOMS+NITEMS
            GOTO 111
         ENDIF
      ELSE 
         NATOMS=0
300      CALL INPUT(END)
         IF (.NOT.END) THEN
            CALL READA(ZDUM)
            IF (ZDUM(1:2).EQ.'  ') GOTO 300
            CALL READF(DUMMY)
            CALL READF(DUMMY)
            CALL READF(DUMMY)
            NATOMS=NATOMS+1
            GOTO 300
         ENDIF
      ENDIF
400   CLOSE(5)


