 
      SUBROUTINE INPUT(END)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT

      LOGICAL BLANK, TCOMMA, END
      DATA CONCAT /'+++'/, LC /3/
      SAVE CONCAT, LC
 
 
 
 
 
 
 
      INTEGER LEVEL
      SAVE LEVEL, IR0

      DATA LEVEL /0/
      DATA SPACE /' '/, BRA /'('/, KET /')'/, COMMA /','/,
     *    SQUOTE /''''/, DQUOTE /'"'/

      ECHO=.FALSE.

      END=.FALSE.
10    M=1
20    LAST=M+79
      READ (IR,1001,END=900) (CHAR(I), I=M,LAST)
      IF (ECHO) PRINT '(1X,80A1)', (CHAR(I), I=M,LAST)
1001  FORMAT (80A1)
      LINE=LINE+1
30    IF (CHAR(LAST) .EQ. SPACE) THEN
        LAST=LAST-1
        IF (LAST .GT. 1) GO TO 30
      ENDIF
      IF (LC .GT. 0 .AND. LAST .GE. LC) THEN
        L=LC
        M=LAST
40      IF (CHAR(M) .EQ. CONCAT(L:L)) THEN
          IF (L .EQ. 1) THEN
            GO TO 20
          ENDIF
          L=L-1
          M=M-1
          GO TO 40
        ENDIF
      ENDIF

      L=0
50    L=L+1
      IF (CHAR(L) .EQ. SPACE .AND. L .LT. LAST) GO TO 50
      IF (CHAR(L) .NE. '#') GO TO 70
      W=' '
      I=0
52    I=I+1
      W(I:I)=CHAR(L)
      L=L+1
      IF (CHAR(L) .NE. SPACE) GO TO 52
60    L=L+1
      F=' '
      IF (L .GT. LAST) GO TO 67
      IF (CHAR(L) .EQ. SPACE) GO TO 60
      IF (CHAR(L) .EQ. SQUOTE .OR. CHAR(L) .EQ. DQUOTE) THEN
        TERM=CHAR(L)
        L=L+1
      ELSE
        TERM=SPACE
      ENDIF
      I=0
65    IF (L .GT. LAST .OR. CHAR(L) .EQ. TERM) GO TO 67
      I=I+1
      F(I:I)=CHAR(L)
      L=L+1
      GO TO 65
67    IF (W .EQ. '#INCLUDE') THEN
        IF (F .EQ. ' ') THEN
           PRINT '(/A)', 'No filename given in #include directive'
           STOP
         ENDIF
        IF (LEVEL .EQ. 0) IR0=IR
        LEVEL=LEVEL+1
        IR=90+LEVEL
        OPEN (IR,FILE=F,STATUS='UNKNOWN')
      ELSE IF (W .EQ. '#CONCAT') THEN
        LC=I
      ELSE IF (W .EQ. '#REVERT') THEN
        LEVEL=LEVEL-1
        IF (LEVEL .EQ. 0) THEN
          IR=IR0
        ELSE
          IR=90+LEVEL
        ENDIF
      ELSE IF (W .EQ. '#ECHO') THEN
        ECHO=.TRUE.
         IF (F .EQ. 'OFF') ECHO=.FALSE.
      ELSE
        GO TO 70
      ENDIF
      GO TO 10

70    ITEM=0
      NITEM=0
      NCR=0
      L=0
 
80    TCOMMA=.TRUE.
90    BLANK=.TRUE.
      TERM=SPACE
 
100   L=L+1
      IF (L .GT. LAST) GO TO 800
      IF (BLANK) GO TO 200
 
      TCOMMA=CHAR(L) .EQ. COMMA .AND. TERM .EQ. SPACE
      BLANK=TCOMMA .OR. CHAR(L) .EQ. TERM
      GO TO 100
 
200   TERM=SPACE
      BLANK=CHAR(L) .EQ. SPACE
      IF (BLANK) GO TO 100
      IF (CHAR(L) .EQ. BRA) GO TO 400
 
      NITEM=NITEM+1
      LOC(NITEM)=L
 
      IF (CHAR(L) .EQ. SQUOTE .OR. CHAR(L) .EQ. DQUOTE) TERM=CHAR(L)
 
      IF (CHAR(L) .NE. COMMA) GO TO 100
      IF (TCOMMA) LOC(NITEM)=0
      GO TO 80
 
400   NEST=1
410   L=L+1
      IF (L .GT. LAST) GO TO 800
      IF (CHAR(L) .EQ. BRA) NEST=NEST+1
      IF (CHAR(L) .EQ. KET) NEST=NEST-1
      IF (NEST .GT. 0) GO TO 410
      GO TO 90
 
800   IF (SKIPBL .AND. NITEM .EQ. 0) GO TO 10
      RETURN
 
900   IF (CAT) THEN
        PRINT '(A)', 'Apparently concatenating at end-of-file'
      ENDIF
      IF (LEVEL .GT. 0) THEN
        LEVEL=LEVEL-1
        IF (LEVEL .EQ. 0) THEN
          IR=IR0
        ELSE
          IR=90+LEVEL
        ENDIF
        GO TO 10
      ENDIF
      DO 950 L=1,LAST
950   CHAR(L)=SPACE
      ITEM=0
      NITEM=-1
      END=.TRUE.
      RETURN
      END
 
      BLOCK DATA INBLK

      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT

      DATA ITEM /0/, NITEM /0/, CHAR /800*' '/, LOC /80*0/, LINE /0/
      DATA SKIPBL /.FALSE./, CLEAR /.TRUE./, NCR /0/, NERROR /0/
      DATA IR /5/, ECHO /.FALSE./
      END
 
 
 
 
 
 
      SUBROUTINE READF(A)
      INTEGER P, STATE
      LOGICAL NULL, XNULL
      DOUBLE PRECISION A, AD, B, TEN, SIGN
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT

     &    S1, S2, DIGIT(10)
      DATA DIGIT /'1', '2', '3', '4', '5', '6', '7', '8', '9', '0'/
      DATA PLUS /'+'/, MINUS /'-'/, DOT /'.'/, SPACE /' '/, STAR /'*'/,
     *    D /'D'/, E /'E'/, COMMA /','/
      DATA TEN /10D0/

      IF (CLEAR) A=0D0
 
      IF (ITEM .GE. NITEM) GO TO 110
      ITEM=ITEM+1
      NCR=0
      P=LOC(ITEM)
      IF (P .EQ. 0) GO TO 110
 
      B=0D0
      SIGN=1D0
      NULL=.TRUE.
      XNULL=.FALSE.
      KXIMP=0
      KX=0
      KXSIGN=1
      STATE=1
      GO TO 12
 
10    P=P+1
      IF (P .GT. LAST) GO TO 100
 
12    C=CHAR(P)
      IF (C .EQ. SPACE .OR. C .EQ. COMMA) GO TO 100
 
      DO 20 I=1,10
      IF (C .EQ. DIGIT(I)) GO TO (70,72,71,80,81), STATE
20    CONTINUE
 
 
 
      IF (C .EQ. DOT) GO TO (30,30,99,99,99), STATE
      IF (C .EQ. PLUS) GO TO (42,52,52,52,99), STATE
      IF (C .EQ. MINUS) GO TO (41,51,51,51,99), STATE
      IF (C .EQ. D .OR. C .EQ. E) GO TO (99,60,60,99,99), STATE
 
99    IF (NERROR .LE. 1) THEN
        WRITE (6,'(1X,A,A,A)') 'Invalid character "', C,
     &                             '" while reading number.         '
        I2=MIN(LAST,P+20)
        I1=MAX(1,I2-70)
        S1=' '
        IF (I1 .GT. 1) S1='...'
        S2=' '
        IF (I2 .LT. LAST) S2='...'
        PRINT '(1X,A/1X,A,1X,80A,1X,A)',
     &           'Current input line is:', S1, (CHAR(I), I=I1,I2), S2
        PRINT '(5X,80A1)', (' ', I=I1,P), '*'
      ENDIF
      IF (NERROR .LE. 0) THEN
        STOP
      ELSE IF (NERROR .EQ. 1) THEN
        WRITE (6,'(1X,A)') 'Item replaced by zero.'
        A=0D0
        RETURN
      ELSE
        NERROR=-1
        A=0D0
        RETURN
      ENDIF
 
 
30    STATE=3
      GO TO 10
 
41    SIGN=-1D0
42    STATE=2
      GO TO 10
 
51    KXSIGN=-1
52    STATE=5
      IF (NULL) GO TO 99
      XNULL=.TRUE.
      GO TO 10
 
60    STATE=4
      GO TO 10
 
 
70    STATE=2
      GO TO 72
71    KXIMP=KXIMP-1
72    AD=I
      IF (I .EQ. 10) AD=0D0
      B=B*TEN+AD
      NULL=.FALSE.
      GO TO 10
 
 
80    STATE=5
81    IF (I .EQ. 10) I=0
      KX=KX*10+I
      XNULL=.FALSE.
      GO TO 10
 
 
100   IF (NULL .OR. XNULL) GO TO 99
      KX=KX*KXSIGN+KXIMP
      A=SIGN*B*(TEN**KX)
 
110   RETURN
      END
 
 
      SUBROUTINE INPUTI(I)
      DOUBLE PRECISION A
      A=I
      I=A
      RETURN
      END
 
 
 
 
 
 
      SUBROUTINE READK(C,M,N,K)
 
 
 
 
 
      INTEGER M(*), N(*), TP, P
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT
      DATA SPACE /' '/

      K=-1
      IF (ITEM .GE. NITEM) RETURN
      ITEM=ITEM+1
      NCR=0
      P=LOC(ITEM)
      IF (P .EQ. 0) RETURN
      K=0
      TP=1
      GO TO 20
 
10    P=P+1
20    IF (P .LE. LAST) GO TO 30
      IF (TP .EQ. 1) K=-1
      RETURN
 
30    CALL MYUPCASE(CHAR(P))
40    IF(CHAR(P) .EQ. C(TP)) GO TO 50
      TP=N(TP)
      IF (TP .EQ. 0) RETURN
      IF (TP .GT. 0) GO TO 40
      GO TO 70
 
50    TP=M(TP)
      IF (TP .GT. 0) GO TO 10
 
70    K=-TP
      RETURN
      END
 
 
      SUBROUTINE READA(M)
      IMPLICIT NONE
 
      INTEGER I, K, NCR, NERROR, ITEM, NITEM, LOC, IR, LAST, L, LINE
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT
      LOGICAL QUOTE
      DATA SPACE /' '/, COMMA /','/, SQUOTE /''''/, DQUOTE /'"'/


      IF (ITEM .GE. NITEM .AND. .NOT. CLEAR) RETURN
      K=0
      IF (ITEM .GE. NITEM) GOTO 65
      ITEM=ITEM+1
      NCR=0
      L=LOC(ITEM)
      IF (L .EQ. 0) RETURN
 
      QUOTE=(CHAR(L) .EQ. SQUOTE .OR. CHAR(L) .EQ. DQUOTE)
      IF (.NOT. QUOTE) THEN
         K=K+1
         GOTO 10
      ENDIF
      TERM=CHAR(L)
      QUOTE=.TRUE.
      L=L+1
      GO TO 11
 
10    CONTINUE
11    IF (L+K .GT. LAST) GO TO 50


      IF ((.NOT. QUOTE .AND. CHAR(L+K) .NE. SPACE
     *                .AND. CHAR(L+K) .NE. COMMA).OR.
     *    (QUOTE .AND. CHAR(L+K) .NE. TERM)) THEN
         K=K+1
         GOTO 10
      ENDIF
 
50    IF (K .EQ. 0) RETURN
      K=MIN0(K,LEN(M))
      DO 60 I=1,K
60    M(I:I)=CHAR(L+I-1)
65    DO 70 I=K+1,LEN(M)
70    M(I:I)=SPACE
      RETURN
      END
 
 
      SUBROUTINE READU(M)
 
      RETURN

      ENTRY READL(M)
      RETURN
      END
 
 
      SUBROUTINE READCH(M)
 

      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT
      M=' '
      IF (ITEM .GE. NITEM .AND. NCR .EQ. 0) RETURN
      IF (NCR .EQ. 0) ITEM=ITEM+1
      L=LOC(ITEM)
      M=CHAR(L+NCR)
      NCR=NCR+1
      RETURN
      END
 
 
      SUBROUTINE GETF(A)
      DOUBLE PRECISION A

      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT

      LOGICAL END
 
10    IF (ITEM .LT. NITEM) THEN
        RETURN
      ENDIF
      IF (END) THEN
        WRITE (6,1001)
1001    FORMAT ('0End of file while attempting to read a number')
        STOP
      ENDIF
      GO TO 10
 
      END
 
 
      SUBROUTINE GETS(S)
      DOUBLE PRECISION A
      REAL S
      A=S
      S=A
      RETURN
      END
 
 
      SUBROUTINE GETI(I)
      DOUBLE PRECISION A
      A=I
      I=A
      RETURN
      END
 
      SUBROUTINE GETA(M)

      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT

      LOGICAL END

10    IF (ITEM .LT. NITEM) THEN
        RETURN
      ENDIF
      IF (END) THEN
        WRITE (6,1001)
1001    FORMAT
     &    ('0End of file while attempting to read a character string')
        STOP
      ENDIF
      GO TO 10
 
      END
 
 
      SUBROUTINE READSNGL(S)
      DOUBLE PRECISION A
      REAL S
      A=S
      S=A
      RETURN
      END
 
 
      SUBROUTINE READI(I)
      DOUBLE PRECISION A
      A=I
      I=A
      RETURN
      END
 
 
      SUBROUTINE REREAD(K)
 
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT
 
      IF (K .LT. 0) ITEM=ITEM+K
      IF (K .EQ. 0) ITEM=ITEM-1
      IF (K .GT. 0) ITEM=K-1
      IF (ITEM .LT. 0) ITEM=0
      NCR=0
      RETURN
      END
 
 
      SUBROUTINE MYUPCASE(WORD)
 
      DATA UC /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA LC /'abcdefghijklmnopqrstuvwxyz'/
 
      DO 10 I=1,LEN(WORD)
      K=INDEX(LC,WORD(I:I))
      IF (K .NE. 0) WORD(I:I)=UC(K:K)
10    CONTINUE
      RETURN
 
      ENTRY LOCASE(WORD)
 
      DO 20 I=1,LEN(WORD)
      K=INDEX(UC,WORD(I:I))
      IF (K .NE. 0) WORD(I:I)=LC(K:K)
20    CONTINUE
      RETURN
 
      END

      SUBROUTINE REPORT(C,REFLCT)
      LOGICAL REFLCT

      LOGICAL SKIPBL, CLEAR, ECHO, CAT
     &                NERROR, IR, ECHO, LAST, CAT

      INTEGER PRINT
      LOGICAL SWITCH, ONLINE
      EQUIVALENCE (ONLINE, SWITCH(4))


      PRINT '(1X,A)', C
      IF (REFLCT) THEN
        L=LOC(ITEM)
        I2=MIN(LAST,L+20)
        I1=MAX(1,I2-70)
        S1=' '
        IF (I1 .GT. 1) S1='...'
        S2=' '
        IF (I2 .LT. LAST) S2='...'
        PRINT '(1X,A/1X,A,1X,80A,1X,A)',
     &           'Current input line is:', S1, (BUFF(I), I=I1,I2), S2
        PRINT '(4X,80A1)', (' ', I=I1,L), '*'
      ENDIF
      END
