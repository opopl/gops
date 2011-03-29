C GPL LICENSE INFO {{{
C   GMIN: A PROGRAM FOR FINDING GLOBAL MINIMA
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF GMIN.
C
C   GMIN IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   GMIN IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C  INPUT  E3   ( 11:41:40 TUESDAY 7 SEPTEMBER 1993 )
C }}}
 
C  DOXYGEN COMMENTS {{{
C> \NAME SUBROUTINE INPUT(LOGICAL END)
C>
C>  READ NEXT INPUT RECORD FROM UNIT IR INTO BUFFER, AND ANALYSE FOR DATA
C>  ITEMS, NULL ITEMS (ENCLOSED BY COMMAS), AND COMMENT (ENCLOSED IN
C>  PARENTHESES, WHICH MAY BE NESTED, AND OCCURRING WHERE SPACES MAY
C>  OCCUR).
C>
C>  E1  ATTEMPTS TO HANDLE STREAM-SWITCHING COMMANDS IN THE DATA.
C>      #INCLUDE FILE-NAME
C>         SWITCHES INPUT TO BE FROM THE FILE SPECIFIED;
C>      #REVERT
C>         (OR END-OF-FILE) REVERTS TO THE PREVIOUS FILE.
C>    ALSO
C>      #CONCAT "STRING"
C>         SETS THE CONCATENATION STRING; E.G.
C>         #CONCAT "\"
C> \PARAM ITEM    IS THE NUMBER OF THE LAST ITEM READ FROM THE BUFFER
C> \PARAM NITEM   IS THE NUMBER OF ITEMS IN THE BUFFER
C> \PARAM CHAR(I) IS THE ITH CHARACTER IN THE BUFFER
C> \PARAM LOC(I)  IS THE STARTING POSITION OF THE ITH ITEM
C> \PARAM LINE    IS THE NUMBER OF LINES (RECORDS) READ
C>
C>  IF \B SKIPBL IS SET TO \B TRUE, LINES CONTAINING NO ITEMS OTHER THAN
C>  COMMENT ARE SKIPPED AUTOMATICALLY, SO THAT INPUT WILL ALWAYS RETURN
C>  A LINE WITH AT LEAST ONE ITEM (POSSIBLY NULL).  IF SKIPBL IS FALSE
C>  (DEFAULT) NO SKIPPING OCCURS AND THE NEXT DATA LINE IS RETURNED
C>  REGARDLESS OF WHAT IT CONTAINS.
C>
C>  IF \B CLEAR IS SET TO \B TRUE (DEFAULT) THEN AN ATTEMPT TO READ MORE THAN
C>  NITEM ITEMS FROM A LINE WILL CAUSE ZEROS OR BLANKS TO BE RETURNED. IF
C>  CLEAR IS FALSE, THE VARIABLES IN QUESTION ARE LEFT UNALTERED.
C>
C>  NCR IS THE NUMBER OF SINGLE CHARACTERS WHICH HAVE BEEN READ FROM THE
C>  CURRENT ITEM.  IT IS ALWAYS ZERO UNLESS READCH HAS BEEN USED; IF NON-
C>  ZERO, ITEM REFERS TO THE CURRENT ITEM FROM WHICH CHARACTERS ARE BEING
C>  READ.
C>
C>  NERROR SPECIFIES THE TREATMENT OF ERRORS FOUND WHILE READING NUMBERS:
C>    0  HARD ERROR - PRINT MESSAGE AND STOP (DEFAULT).
C>    1  SOFT ERROR - PRINT MESSAGE AND RETURN ZERO.
C>    2  SOFT ERROR - NO MESSAGE, RETURN ZERO, SET NERROR TO -1.
C>  IF NERROR IS SET TO 2 IT IS IMPORTANT TO TEST AND RESET IT AFTER READING
C>  A NUMBER.
C>
C>  IR IS THE INPUT STREAM FROM WHICH THE DATA ARE TO BE READ (DEFAULT 5).
C>
C>  IF ECHO IS TRUE, THE INPUT LINE WILL BE REFLECTED TO THE SYSTEM OUTPUT
C>  STREAM.
C>  LAST GIVES THE POSITION OF THE LAST NON-BLANK CHARACTER IN THE LINE.
C>
C>  CONCAT IS THE LINE CONCATENATION STRING: IF IT IS FOUND AS THE LAST
C>  NON-BLANK CHARACTER STRING ON A LINE, THE FOLLOWING LINE IS READ AND
C>  APPENDED TO THE CURRENT LINE, OVERWRITING THE CONCATENATION STRING.
C>  THIS PROCEDURE IS REPEATED IF NECESSARY UNTIL A LINE IS FOUND THAT DOES
C>  NOT END WITH THE CONCATENATION STRING.
C> 
C>  THE CONCATENATION STRING MAY BE MODIFIED BY CHANGING THE DATA STATEMENT
C>  ABOVE. THE INTEGER LC MUST BE SET TO THE NUMBER OF NON-BLANK CHARACTERS
C>  IN THE CONCATENATION STRING. IF IT IS SET TO ZERO, NO CONCATENATION
C>  OCCURS. THE CONCATENATION STRING MAY ALSO BE CHANGED BY THE #CONCAT
C>  DIRECTIVE IN THE DATA FILE.
C> 
C }}}
      SUBROUTINE INPUT(END)
C  COMMENTS, DECLARATIONS {{{
C
C  READ NEXT INPUT RECORD FROM UNIT IR INTO BUFFER, AND ANALYSE FOR DATA
C  ITEMS, NULL ITEMS (ENCLOSED BY COMMAS), AND COMMENT (ENCLOSED IN
C  PARENTHESES, WHICH MAY BE NESTED, AND OCCURRING WHERE SPACES MAY
C  OCCUR).
C
C  E1  ATTEMPTS TO HANDLE STREAM-SWITCHING COMMANDS IN THE DATA.
C      #INCLUDE FILE-NAME
C         SWITCHES INPUT TO BE FROM THE FILE SPECIFIED;
C      #REVERT
C         (OR END-OF-FILE) REVERTS TO THE PREVIOUS FILE.
C    ALSO
C      #CONCAT "STRING"
C         SETS THE CONCATENATION STRING; E.G.
C         #CONCAT "\"
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      LOGICAL BLANK, TCOMMA, END
      CHARACTER(LEN=8) CONCAT
      DATA CONCAT /'+++'/, LC /3/
      SAVE CONCAT, LC
 
C     ITEM    IS THE NUMBER OF THE LAST ITEM READ FROM THE BUFFER
C     NITEM   IS THE NUMBER OF ITEMS IN THE BUFFER
C     CHAR(I) IS THE ITH CHARACTER IN THE BUFFER
C     LOC(I)  IS THE STARTING POSITION OF THE ITH ITEM
C     LINE    IS THE NUMBER OF LINES (RECORDS) READ
C
C  IF SKIPBL IS SET TO TRUE, LINES CONTAINING NO ITEMS OTHER THAN
C
C  COMMENT ARE SKIPPED AUTOMATICALLY, SO THAT INPUT WILL ALWAYS RETURN
C  A LINE WITH AT LEAST ONE ITEM (POSSIBLY NULL).  IF SKIPBL IS FALSE
C  (DEFAULT) NO SKIPPING OCCURS AND THE NEXT DATA LINE IS RETURNED
C  REGARDLESS OF WHAT IT CONTAINS.
C
C  IF CLEAR IS SET TO TRUE (DEFAULT) THEN AN ATTEMPT TO READ MORE THAN
C  NITEM ITEMS FROM A LINE WILL CAUSE ZEROS OR BLANKS TO BE RETURNED. IF
C  CLEAR IS FALSE, THE VARIABLES IN QUESTION ARE LEFT UNALTERED.
C
C  NCR IS THE NUMBER OF SINGLE CHARACTERS WHICH HAVE BEEN READ FROM THE
C  CURRENT ITEM.  IT IS ALWAYS ZERO UNLESS READCH HAS BEEN USED; IF NON-
C  ZERO, ITEM REFERS TO THE CURRENT ITEM FROM WHICH CHARACTERS ARE BEING
C  READ.
C
C  NERROR SPECIFIES THE TREATMENT OF ERRORS FOUND WHILE READING NUMBERS:
C    0  HARD ERROR - PRINT MESSAGE AND STOP (DEFAULT).
C    1  SOFT ERROR - PRINT MESSAGE AND RETURN ZERO.
C    2  SOFT ERROR - NO MESSAGE, RETURN ZERO, SET NERROR TO -1.
C  IF NERROR IS SET TO 2 IT IS IMPORTANT TO TEST AND RESET IT AFTER READING
C  A NUMBER.
 
C  IR IS THE INPUT STREAM FROM WHICH THE DATA ARE TO BE READ (DEFAULT 5).
 
C  IF ECHO IS TRUE, THE INPUT LINE WILL BE REFLECTED TO THE SYSTEM OUTPUT
C  STREAM.
 
C  LAST GIVES THE POSITION OF THE LAST NON-BLANK CHARACTER IN THE LINE.
 
C  CONCAT IS THE LINE CONCATENATION STRING: IF IT IS FOUND AS THE LAST
C  NON-BLANK CHARACTER STRING ON A LINE, THE FOLLOWING LINE IS READ AND
C  APPENDED TO THE CURRENT LINE, OVERWRITING THE CONCATENATION STRING.
C  THIS PROCEDURE IS REPEATED IF NECESSARY UNTIL A LINE IS FOUND THAT DOES
C  NOT END WITH THE CONCATENATION STRING.
 
C  THE CONCATENATION STRING MAY BE MODIFIED BY CHANGING THE DATA STATEMENT
C  ABOVE. THE INTEGER LC MUST BE SET TO THE NUMBER OF NON-BLANK CHARACTERS
C  IN THE CONCATENATION STRING. IF IT IS SET TO ZERO, NO CONCATENATION
C  OCCURS. THE CONCATENATION STRING MAY ALSO BE CHANGED BY THE #CONCAT
C  DIRECTIVE IN THE DATA FILE.
 
      INTEGER LEVEL
      SAVE LEVEL, IR0
      CHARACTER(LEN=80) W, F

      CHARACTER SPACE, BRA, KET, COMMA, SQUOTE, DQUOTE, TERM
      DATA LEVEL /0/
      DATA SPACE /' '/, BRA /'('/, KET /')'/, COMMA /','/,
     *    SQUOTE /''''/, DQUOTE /'"'/
C }}}

C {{{ 
      ECHO=.FALSE.

      END=.FALSE.
      CAT=.FALSE.
10    M=1
20    LAST=M+79
      READ (IR,1001,END=900) (CHAR(I), I=M,LAST)
      IF (ECHO) PRINT '(1X,80A1)', (CHAR(I), I=M,LAST)
1001  FORMAT (80A1)
      CAT=.FALSE.
      LINE=LINE+1
C  FIND LAST NON-BLANK CHARACTER
30    IF (CHAR(LAST) .EQ. SPACE) THEN
        LAST=LAST-1
        IF (LAST .GT. 1) GO TO 30
      ENDIF
C  LOOK FOR CONCATENATION STRING
      IF (LC .GT. 0 .AND. LAST .GE. LC) THEN
        L=LC
        M=LAST
40      IF (CHAR(M) .EQ. CONCAT(L:L)) THEN
          IF (L .EQ. 1) THEN
            CAT=.TRUE.
            GO TO 20
          ENDIF
          L=L-1
          M=M-1
          GO TO 40
        ENDIF
      ENDIF
      CAT=.FALSE.

C  LOGICAL LINE ASSEMBLED. FIRST LOOK FOR INPUT DIRECTIVES
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
      CALL MYUPCASE(W)
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
           PRINT '(/A)', 'NO FILENAME GIVEN IN #INCLUDE DIRECTIVE'
           STOP
         ENDIF
        IF (LEVEL .EQ. 0) IR0=IR
        LEVEL=LEVEL+1
        IR=90+LEVEL
        OPEN (IR,FILE=F,STATUS='UNKNOWN')
      ELSE IF (W .EQ. '#CONCAT') THEN
        CONCAT=F
        LC=I
      ELSE IF (W .EQ. '#REVERT') THEN
        CLOSE(IR)
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

C  ANALYSE INPUT 
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
 
C  READING THROUGH AN ITEM, SEEKING TERMINATOR.
C  (NEXT LINE OF CODE SUPPRESSED: COMMENT NOT ALLOWED WITHIN AN ITEM)
C     IF (CHAR(L) .EQ. BRA) GO TO 400
      TCOMMA=CHAR(L) .EQ. COMMA .AND. TERM .EQ. SPACE
      BLANK=TCOMMA .OR. CHAR(L) .EQ. TERM
      GO TO 100
 
C  LOOKING FOR NEXT ITEM
200   TERM=SPACE
      BLANK=CHAR(L) .EQ. SPACE
      IF (BLANK) GO TO 100
      IF (CHAR(L) .EQ. BRA) GO TO 400
 
C  ITEM FOUND
      NITEM=NITEM+1
      LOC(NITEM)=L
 
C  QUOTED STRING?
      IF (CHAR(L) .EQ. SQUOTE .OR. CHAR(L) .EQ. DQUOTE) TERM=CHAR(L)
 
C  NULL ITEM?
      IF (CHAR(L) .NE. COMMA) GO TO 100
      IF (TCOMMA) LOC(NITEM)=0
      GO TO 80
 
C  COMMENT (ENCLOSED IN PARENTHESES, POSSIBLY NESTED)
400   NEST=1
410   L=L+1
      IF (L .GT. LAST) GO TO 800
      IF (CHAR(L) .EQ. BRA) NEST=NEST+1
      IF (CHAR(L) .EQ. KET) NEST=NEST-1
      IF (NEST .GT. 0) GO TO 410
      GO TO 90
 
C  END OF CARD -- WAS IT BLANK?
800   IF (SKIPBL .AND. NITEM .EQ. 0) GO TO 10
      RETURN
 
C  END OF DATA
900   IF (CAT) THEN
        PRINT '(A)', 'APPARENTLY CONCATENATING AT END-OF-FILE'
        CALL REPORT('UNEXPECTED END OF DATA FILE',.TRUE.)
      ENDIF
      IF (LEVEL .GT. 0) THEN
C  REVERT TO PREVIOUS INPUT
        CLOSE(IR)
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
C }}} 
      END
 
C-----------------------------------------------------------------------
C> \NAME BLOCK DATA INBLK 
      BLOCK DATA INBLK

      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      DATA ITEM /0/, NITEM /0/, CHAR /800*' '/, LOC /80*0/, LINE /0/
      DATA SKIPBL /.FALSE./, CLEAR /.TRUE./, NCR /0/, NERROR /0/
      DATA IR /5/, ECHO /.FALSE./
      END
 
C-----------------------------------------------------------------------
 
C     SUBROUTINE STREAM(N)
 
C  SET THE INPUT STREAM FOR SUBSEQUENT DATA TO BE N.
 
C     LOGICAL SKIPBL, CLEAR, ECHO, CAT
C     COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
C    &                NERROR, IR, ECHO, LAST, CAT
C
C     IR=N
C     RETURN
C     END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READF(A)
C> \NAME READF
C> \BRIEF READ THE NEXT ITEM FROM THE BUFFER AS A REAL NUMBER
C> \PARAM A (DP)
C>
C DECLARATIONS {{{ 
      INTEGER P, STATE
      LOGICAL NULL, XNULL
      DOUBLE PRECISION A, AD, B, TEN, SIGN
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      CHARACTER C, PLUS, MINUS, DOT, SPACE, STAR, D, E, COMMA,
     &    S1, S2, DIGIT(10)
      DATA DIGIT /'1', '2', '3', '4', '5', '6', '7', '8', '9', '0'/
      DATA PLUS /'+'/, MINUS /'-'/, DOT /'.'/, SPACE /' '/, STAR /'*'/,
     *    D /'D'/, E /'E'/, COMMA /','/
      DATA TEN /10D0/
C }}}

C SUBROUTINE BODY {{{ 
      IF (CLEAR) A=0D0
 
C  IF THE ITEM IS NULL, OR IF NO ITEM REMAINS, A IS UNCHANGED
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
 
C  TERMINATORS
12    C=CHAR(P)
      IF (C .EQ. SPACE .OR. C .EQ. COMMA) GO TO 100
 
C  DIGITS
      DO 20 I=1,10
      IF (C .EQ. DIGIT(I)) GO TO (70,72,71,80,81), STATE
20    CONTINUE
 
C  STATE NUMBERS:  1  NOTHING READ YET
C                  2  READING DIGITS OF INTEGRAL PART
C                  3  READING DECIMAL FRACTION
C                  4  EXPECTING EXPONENT
C                  5  READING DIGITS OF EXPONENT
 
C  NUMBER CONTROL CHARACTERS -- BRANCH ON STATE
 
      IF (C .EQ. DOT) GO TO (30,30,99,99,99), STATE
      IF (C .EQ. PLUS) GO TO (42,52,52,52,99), STATE
      IF (C .EQ. MINUS) GO TO (41,51,51,51,99), STATE
      IF (C .EQ. D .OR. C .EQ. E) GO TO (99,60,60,99,99), STATE
 
C  ERROR
99    IF (NERROR .LE. 1) THEN
        WRITE (6,'(1X,A,A,A)') 'INVALID CHARACTER "', C,
     &                             '" WHILE READING NUMBER.         '
        I2=MIN(LAST,P+20)
        I1=MAX(1,I2-70)
        S1=' '
        IF (I1 .GT. 1) S1='...'
        S2=' '
        IF (I2 .LT. LAST) S2='...'
        PRINT '(1X,A/1X,A,1X,80A,1X,A)',
     &           'CURRENT INPUT LINE IS:', S1, (CHAR(I), I=I1,I2), S2
        PRINT '(5X,80A1)', (' ', I=I1,P), '*'
      ENDIF
      IF (NERROR .LE. 0) THEN
        STOP
      ELSE IF (NERROR .EQ. 1) THEN
        WRITE (6,'(1X,A)') 'ITEM REPLACED BY ZERO.'
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
C }}} 
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE INPUTI(I)
      DOUBLE PRECISION A
      A=I
      CALL READF(A)
      I=A
      RETURN
      END
 
C-----------------------------------------------------------------------
 
C DOXYGEN COMMENTS {{{ 
C> \NAME READK
C> 
C>  THE ARRAYS C, M AND N COMPRISE A TREE, EACH NODE OF WHICH CONSISTS OF
C>  A CHARACTER C(I) AND TWO INTEGERS M(I) AND N(I). INITIALLY I=1. THE
C>  CHARACTER C(I) IS TESTED AGAINST THE CURRENT CHARACTER IN THE INPUT.
C>  IF IT MATCHES, THE NEXT INPUT CHARACTER IS SELECTED AND M(I) IS
C>  EXAMINED; OTHERWISE N(I) IS EXAMINED. IF A POSITIVE VALUE IS FOUND,
C>  IT IS A POINTER TO THE NEXT NODE TO BE TESTED; IF NEGATIVE OR ZERO,
C>  IT IS MINUS THE KEY VALUE TO BE RETURNED IN K. IF A KEYWORD IS
C>  PRESENT IN THE INPUT BUT DOES NOT MATCH A KEYWORD IN THE TREE, THEN
C>  ZERO IS RETURNED; IF A NULL ITEM IS PRESENT OR IF THERE ARE NO MORE
C>  ITEMS ON THE LINE, THEN -1 IS RETURNED. NOTE THAT THIS ARRANGEMENT
C>  ALLOWS THE CODE IN THE TREE TO BE SET UP IN MANY WAYS TO SUIT THE
C>  REQUIREMENTS OF THE PROBLEM. ABBREVIATIONS MAY BE DEFINED
C>  ECONOMICALLY, AND THE KEYWORD CAN BE REQUIRED TO TERMINATE WITH A
C>  SPACE OR OTHER SYMBOL, OR IT CAN BE SELF-TERMINATING.
 
C>  EXAMPLE:  PRINT AND PRT ARE TO RETURN KEY 19, AND MUST BE TERMINATED
C>  BY SPACE.  P ALONE (I.E. NOT PRESENT AS THE FIRST LETTER OF PRINT OR
C>  PRT) IS TO RETURN 7.  INITIAL SPACES ARE TO BE IGNORED.  THE CODE
C>  TABLE MIGHT BEGIN
 
C>        1   2   3   4   5   6   7   8
 
C>   C   ' '  P   R   I   N   T  ' '  ...
C>   M    1   3   4   5   6   7  -19  ...
C>   N    2   8  -7   6   0   0   0   ...
 
C>  SUCH CODE TABLES CAN BE CREATED AUTOMATICALLY BY THE PROCEDURE
C>  ENCODE, WHICH PRODUCES THE DATA STATEMENTS REQUIRED TO INITIALIZE
C>  THE ARRAY.
C }}}
      SUBROUTINE READK(C,M,N,K)
C ORIGINAL COMMENTS {{{ 
C  THE ARRAYS C, M AND N COMPRISE A TREE, EACH NODE OF WHICH CONSISTS OF
C  A CHARACTER C(I) AND TWO INTEGERS M(I) AND N(I). INITIALLY I=1. THE
C  CHARACTER C(I) IS TESTED AGAINST THE CURRENT CHARACTER IN THE INPUT.
C  IF IT MATCHES, THE NEXT INPUT CHARACTER IS SELECTED AND M(I) IS
C  EXAMINED; OTHERWISE N(I) IS EXAMINED. IF A POSITIVE VALUE IS FOUND,
C  IT IS A POINTER TO THE NEXT NODE TO BE TESTED; IF NEGATIVE OR ZERO,
C  IT IS MINUS THE KEY VALUE TO BE RETURNED IN K. IF A KEYWORD IS
C  PRESENT IN THE INPUT BUT DOES NOT MATCH A KEYWORD IN THE TREE, THEN
C  ZERO IS RETURNED; IF A NULL ITEM IS PRESENT OR IF THERE ARE NO MORE
C  ITEMS ON THE LINE, THEN -1 IS RETURNED. NOTE THAT THIS ARRANGEMENT
C  ALLOWS THE CODE IN THE TREE TO BE SET UP IN MANY WAYS TO SUIT THE
C  REQUIREMENTS OF THE PROBLEM. ABBREVIATIONS MAY BE DEFINED
C  ECONOMICALLY, AND THE KEYWORD CAN BE REQUIRED TO TERMINATE WITH A
C  SPACE OR OTHER SYMBOL, OR IT CAN BE SELF-TERMINATING.
 
C  EXAMPLE:  PRINT AND PRT ARE TO RETURN KEY 19, AND MUST BE TERMINATED
C  BY SPACE.  P ALONE (I.E. NOT PRESENT AS THE FIRST LETTER OF PRINT OR
C  PRT) IS TO RETURN 7.  INITIAL SPACES ARE TO BE IGNORED.  THE CODE
C  TABLE MIGHT BEGIN
 
C        1   2   3   4   5   6   7   8
 
C   C   ' '  P   R   I   N   T  ' '  ...
C   M    1   3   4   5   6   7  -19  ...
C   N    2   8  -7   6   0   0   0   ...
 
C  SUCH CODE TABLES CAN BE CREATED AUTOMATICALLY BY THE PROCEDURE
C  ENCODE, WHICH PRODUCES THE DATA STATEMENTS REQUIRED TO INITIALIZE
C  THE ARRAY.
C }}}
 
C DECLARATIONS {{{ 
      INTEGER M(*), N(*), TP, P
      CHARACTER C(*), SPACE
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      DATA SPACE /' '/
C }}} 

C SUBROUTINE BODY {{{
      K=-1
      IF (ITEM .GE. NITEM) RETURN
      ITEM=ITEM+1
      NCR=0
      P=LOC(ITEM)
      IF (P .EQ. 0) RETURN
      K=0
      TP=1
      GO TO 20
 
C  ADVANCE CHARACTER POINTER
10    P=P+1
20    IF (P .LE. LAST) GO TO 30
C  END OF LINE
      IF (TP .EQ. 1) K=-1
      RETURN
 
30    CALL MYUPCASE(CHAR(P))
40    IF(CHAR(P) .EQ. C(TP)) GO TO 50
C  NO MATCH -- ADVANCE TREE POINTER
      TP=N(TP)
C  ZERO POINTER -- UNDECODEABLE
      IF (TP .EQ. 0) RETURN
C  POSITIVE VALUE -- TEST SAME INPUT CHARACTER AGAIN
      IF (TP .GT. 0) GO TO 40
C  NEGATIVE FAIL POINTER: THE INPUT CONTAINS A KEYWORD WHICH IS AN INITIAL
C  SUBSTRING OF SOME OTHER KEYWORD.
      GO TO 70
 
C  MATCHED -- ADVANCE TREE POINTER TO NEXT CHARACTER OF COMMAND
50    TP=M(TP)
      IF (TP .GT. 0) GO TO 10
 
C  LAST CHARACTER OF KEYWORD FOUND.
70    K=-TP
C }}}
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READA(M)
      IMPLICIT NONE
 
C  READ CHARACTERS AND PACK THEM INTO THE CHARACTER VARIABLE M.
C  IF THE FIRST CHARACTER IS A SINGLE OR DOUBLE QUOTE, THE STRING IS
C  TERMINATED BY A MATCHING QUOTE AND THE QUOTES ARE REMOVED; OTHERWISE
C  IT IS TERMINATED BY SPACE OR COMMA.  IF NECESSARY, M IS FILLED
C  OUT WITH SPACES.
C DECLARATIONS {{{ 
      CHARACTER SPACE, BLANK, COMMA, SQUOTE, DQUOTE, TERM
      CHARACTER*(*) M
      CHARACTER CHAR
      INTEGER I, K, NCR, NERROR, ITEM, NITEM, LOC, IR, LAST, L, LINE
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      LOGICAL QUOTE
      DATA SPACE /' '/, COMMA /','/, SQUOTE /''''/, DQUOTE /'"'/
C }}}

C SUBROUTINE BODY {{{ 

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
C
C  THE NEW SUN COMPILER REFUSES TO EXECUTE THE FOLLOWING LINE - HENCE REWRITTEN
C
C     K=K+1
C
11    IF (L+K .GT. LAST) GO TO 50

C     IF (.NOT. QUOTE .AND. CHAR(L+K) .NE. SPACE
C    *                .AND. CHAR(L+K) .NE. COMMA) GO TO 10
C     IF (QUOTE .AND. CHAR(L+K) .NE. TERM) GO TO 10

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
C }}}
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READU(M)
      CHARACTER*(*) M
 
      CALL READA(M)
      CALL MYUPCASE(M)
      RETURN

      ENTRY READL(M)
      CALL READA(M)
      CALL LOCASE(M)
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READCH(M)
 
C  READ A SINGLE CHARACTER FROM THE NEXT (OR THE CURRENT) DATA ITEM.
C  NO ACCOUNT IS TAKEN OF SPECIAL CHARACTERS.

C DECLARATIONS {{{ 
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      CHARACTER*(*) M
C }}} 
      M=' '
      IF (ITEM .GE. NITEM .AND. NCR .EQ. 0) RETURN
      IF (NCR .EQ. 0) ITEM=ITEM+1
      L=LOC(ITEM)
      M=CHAR(L+NCR)
      NCR=NCR+1
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE GETF(A)
C  READ THE NEXT ITEM AS A DOUBLE-PRECISION NUMBER, READING NEW DATA
C  RECORDS IF NECESSARY
C DECLARATIONS  {{{
      DOUBLE PRECISION A

      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      LOGICAL END
C }}}
 
10    IF (ITEM .LT. NITEM) THEN
        CALL READF(A)
        RETURN
      ENDIF
      CALL INPUT(END)
      IF (END) THEN
        WRITE (6,1001)
1001    FORMAT ('0END OF FILE WHILE ATTEMPTING TO READ A NUMBER')
        STOP
      ENDIF
      GO TO 10
 
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE GETS(S)
C  GET A SINGLE-PRECISION NUMBER
      DOUBLE PRECISION A
      REAL S
      A=S
      CALL GETF(A)
      S=A
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE GETI(I)
C  GET AN INTEGER
      DOUBLE PRECISION A
      A=I
      CALL GETF(A)
      I=A
      RETURN
      END
 
C-----------------------------------------------------------------------
C> \NAME GETA
C> \BRIEF GET A CHARACTER STRING 
      SUBROUTINE GETA(M)
C  GET A CHARACTER STRING
C DECLARATIONS {{{
      CHARACTER*(*) M

      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

C }}}
      LOGICAL END

10    IF (ITEM .LT. NITEM) THEN
        CALL READA(M)
        RETURN
      ENDIF
      CALL INPUT(END)
      IF (END) THEN
        WRITE (6,1001)
1001    FORMAT
     &    ('0END OF FILE WHILE ATTEMPTING TO READ A CHARACTER STRING')
        STOP
      ENDIF
      GO TO 10
 
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READSNGL(S)
C>  \NAME READSNGL
C>  \BRIEF READ A NUMBER FROM THE CURRENT RECORD INTO THE SINGLE-PRECISION
C>  \BRIEF VARIABLE S
C  READ A NUMBER FROM THE CURRENT RECORD INTO THE SINGLE-PRECISION
C  VARIABLE S
      DOUBLE PRECISION A
      REAL S
      A=S
      CALL READF(A)
      S=A
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READI(I)
C> \NAME READI
C> \BRIEF READ AN INTEGER FROM THE CURRENT RECORD
      DOUBLE PRECISION A
      A=I
      CALL READF(A)
      I=A
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE REREAD(K)
C> \NAME REREAD 
C>
C>  K>0  REREAD FROM ITEM K
C>  K<0  GO BACK |K| ITEMS
C>  K=0  SAME AS K=-1, I.E. REREAD LAST ITEM.
 
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
 
      IF (K .LT. 0) ITEM=ITEM+K
      IF (K .EQ. 0) ITEM=ITEM-1
      IF (K .GT. 0) ITEM=K-1
      IF (ITEM .LT. 0) ITEM=0
      NCR=0
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE MYUPCASE(WORD)
      CHARACTER WORD*(*)
 
      CHARACTER UC*26, LC*26
      DATA UC /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA LC /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
 
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
C  REPORT  B1  ( 16:10:41 THURSDAY 30 APRIL 1992 )

      SUBROUTINE REPORT(C,REFLCT)
C DECLARATIONS {{{
      CHARACTER C*(*)
      LOGICAL REFLCT

      CHARACTER BUFF
      COMMON /BUFFER/ BUFF(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      INTEGER PRINT
      LOGICAL SWITCH, ONLINE
      COMMON /TEST/ SWITCH(8), PRINT
      EQUIVALENCE (ONLINE, SWITCH(4))

      CHARACTER(LEN=3) S1, S2
C  }}}

      PRINT '(1X,A)', C
      IF (REFLCT) THEN
        L=LOC(ITEM)
        I2=MIN(LAST,L+20)
        I1=MAX(1,I2-70)
        S1=' '
        IF (I1 .GT. 1) S1='...'
        S2=' '
        IF (I2 .LT. LAST) S2='...'
C       PRINT '(1X,4(A,I3))',
C    &      'L =', L, '   LAST =', LAST, '   I1 =', I1, '   I2 =', I2
        PRINT '(1X,A/1X,A,1X,80A,1X,A)',
     &           'CURRENT INPUT LINE IS:', S1, (BUFF(I), I=I1,I2), S2
        PRINT '(4X,80A1)', (' ', I=I1,L), '*'
      ENDIF
C     IF (.NOT. ONLINE) STOP
      END
