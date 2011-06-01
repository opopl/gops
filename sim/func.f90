
MODULE FUNC

IMPLICIT NONE
SAVE

INTEGER :: NUMBER_OF_ATOMS

CONTAINS

! Doxygen RAD  {{{
!
!> @name RAD
!> @brief Add energy and gradient correction terms for the container separately.
!> @param X
!> @param V
!> @param ENERGY
!> @param GTEST
!
! }}}
      SUBROUTINE RAD(X,V,ENERGY,GTEST)
      ! {{{

      USE COMMONS

      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J3
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), 
     1                 ENERGY, DUMMYX, DUMMYY, DUMMYZ
      LOGICAL EVAP, EVAPREJECT
      COMMON /EV/ EVAP, EVAPREJECT

      IF (BSPT) RETURN ! container is accounted for in bspt by recounting previous configuration
      IF (PERIODIC) RETURN
      EVAP=.FALSE.
      EVAPREJECT=.FALSE.
      DO J1=1,NATOMS
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
         IF (DIST.GT.RADIUS) THEN
!           WRITE(MYUNIT,'(A,I5,5G20.10)') 'J1,DIST,RADIUS in rad = ',J1,DIST,RADIUS,X(J3-2),X(J3-1),X(J3)
            EVAP=.TRUE.
           ! IF (EVAP.AND.(BSWL.OR.BSPT)) then
           !     EVAPREJECT=.TRUE.
           !     IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') 'EVAP: atom, radius=',J1,SQRT(DIST)
           !     RETURN
           ! ENDIF

            IF (DEBUG)  WRITE(MYUNIT,'(A,2G20.10,L10)') 'rad> EVAP: atom, radius, EVAP=',J1,SQRT(DIST),EVAP
C           PRINT*,'EVAP: atom, radius=',J1,SQRT(DIST)
CC           ENERGY=ENERGY+1.0D5*(DIST-RADIUS)**2
CC           IF (GTEST.AND.(.NOT.(SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE))) THEN
CC              DUMMYX=1.0D5*4.0D0*(DIST-RADIUS)*X(J3-2)
CC              DUMMYY=1.0D5*4.0D0*(DIST-RADIUS)*X(J3-1)
CC              DUMMYZ=1.0D5*4.0D0*(DIST-RADIUS)*X(J3)
CC              V(J3-2)=V(J3-2)+DUMMYX
CC              V(J3-1)=V(J3-1)+DUMMYY
CC              V(J3)=V(J3)+DUMMYZ
CC           ENDIF

             DIST=(SQRT(RADIUS)-0.5D0)/SQRT(DIST)
             X(J3-2)=X(J3-2)*DIST
             X(J3-1)=X(J3-1)*DIST
             X(J3)=X(J3)*DIST
!            WRITE(MYUNIT,'(A,3G20.10)') 'rad> reset coords: ',X(J3-2),X(J3-1),X(J3)
C
C  Put it back in at the opposite end of a diameter
C
C           X(J3-2)=-X(J3-2)*0.8D0
C           X(J3-1)=-X(J3-1)*0.8D0
C           X(J3)=-X(J3)*0.8D0
         ENDIF
      ENDDO
      ! }}}
      RETURN
      END

!  For rigid-body angle-axis coordinates, just move the fixed site

      SUBROUTINE RADR(X,V,ENERGY,GTEST)
      ! {{{
      USE COMMONS
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J3
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), 
     1                 ENERGY, DUMMYX, DUMMYY, DUMMYZ
      LOGICAL EVAP, evapreject
      COMMON /EV/ EVAP, evapreject

      IF (PERIODIC) RETURN
      EVAP=.FALSE.
      DO J1=1,NATOMS/2
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C        WRITE(*,'(A,I6,5F15.5)') 'J1,DIST,coords,RADIUS in radr=',J1,DIST,X(J3-2),X(J3-1),X(J3),RADIUS
         IF (DIST.GT.RADIUS) THEN
            EVAP=.TRUE.
            WRITE(MYUNIT,'(A,I5,5G17.8)') 'EVAP: molecule, coords, dist, radius=',J1,X(J3-2),X(J3-1),X(J3),SQRT(DIST),SQRT(RADIUS)
C           IF (DEBUG) WRITE(*,'(A,I5,2G20.10)') 'EVAP: molecule, dist, radius=',J1,SQRT(DIST),SQRT(RADIUS)
            DIST=SQRT(RADIUS*0.9D0/DIST)
C           X(J3-2)=X(J3-2)*DIST
C           X(J3-1)=X(J3-1)*DIST
C           X(J3)=X(J3)*DIST
C
C  Put it back in at the opposite end of a diameter
C
            X(J3-2)=-X(J3-2)*0.8D0
            X(J3-1)=-X(J3-1)*0.8D0
            X(J3)=-X(J3)*0.8D0
         ENDIF
      ENDDO

      RETURN
      ! }}}
      END

SUBROUTINE MYSYSTEM(STATUS,DEBUG,JOBSTRING)
! {{{
USE PORFUNCS

IMPLICIT NONE

LOGICAL DEBUG
INTEGER STATUS
CHARACTER(LEN=*) JOBSTRING

IF (DEBUG) WRITE(*,'(A)') 'mysystem> '//trim(adjustl(jobstring)) 
CALL SYSTEM_SUBR(JOBSTRING,STATUS)

! IF (DEBUG) PRINT '(A,I6)','command '//JOBSTRING//' exit status=',STATUS
! IF (STATUS.NE.0) PRINT '(A,I8)','mysystem> WARNING - '//JOBSTRING//' exit status=',STATUS

RETURN
! }}}
END SUBROUTINE MYSYSTEM

SUBROUTINE INQF(FILENAME,YESNO)

LOGICAL,INTENT(OUT) :: YESNO
INTEGER,INTENT(IN) :: FILEHANDLE

FILENAME=TRIM(ADJUSTL(FILENAME))
INQUIRE(FILE=FILENAME,EXIST=YESNO)

END SUBROUTINE INQF

SUBROUTINE OPENF(FILEHANDLE,MODE,FILENAME)
! {{{

INTEGER FILEHANDLE
CHARACTER (LEN=*) MODE,FILENAME

FILENAME=TRIM(ADJUSTL(FILENAME))

SELECTCASE(MODE)
        CASE(">")
                OPEN(UNIT=FILEHANDLE,FILE=FILENAME,STATUS="UNKNOWN",FORM="FORMATTED")
        CASE("O")
                OPEN(UNIT=FILEHANDLE,FILE=FILENAME,STATUS="OLD")
        CASE("<")
                OPEN(UNIT=FILEHANDLE,FILE=FILENAME,STATUS="OLD",ACTION="READ")
        CASE(">>")
                OPEN(FILEHANDLE,FILE=FILENAME,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
        CASE("RW>>")
                OPEN(FILEHANDLE,FILE=FILENAME,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND',ACTION="READWRITE")
        case("DA")
                OPEN(FILEHANDLE,FILE=FILENAME,ACCESS="DIRECT",STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=8*3*NATOMS)
                
ENDSELECT
! }}}
ENDSUBROUTINE OPENF

SUBROUTINE MYRESET(NATOMS,NSEED)
! {{{
USE COMMONS,ONLY : COORDS,COORDSO,VAT,VATO

IMPLICIT NONE

INTEGER IA, NSEED, NATOMS

COORDS(1:NATOMS-NSEED)=COORDSO
VAT=VATO

RETURN
! }}}
END
 
SUBROUTINE TAKESTEP
      ! {{{

      USE COMMONS

      IMPLICIT NONE

      DO IA=1,NATOMS
         CALL GET_RND(RND,3,-1.0D0,1.0D0)
         COORDS(IA,1:3)=COORDS(IA,1:3)+STEP*RND(1:3)
      ENDDO
      
      RETURN
      ! }}}
END SUBROUTINE 

SUBROUTINE CENTRE2(R)
! {{{
USE COMMONS

IMPLICIT NONE

DOUBLE PRECISION RMASS(3), R(NATOMS,3)
INTEGER I,K

RMASS=SUM(R,DIM=1)/NATOMS

do K=1,3
R(:,K)=R(:,K)-RMASS(K)
ENDDO

IF (DEBUG) WRITE(LFH,'(A,3G20.10)') 'centre2> centre of mass reset to the origin from ',&                                                   (RMASS(I),I=1,3)
RETURN
! }}}
END

!> @name GSORT
!> @brief This subprogram performs a sort on the input data and
!> arranges it from smallest to biggest. The exchange-sort algorithm is used.
!
      SUBROUTINE GSORT(N,NATOMS)
! {{{

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NATOMS,N
      !DOUBLE PRECISION(:), INTENT(INOUT) ::
 
      INTEGER J1, L, J3, J2, NTEMP
      DOUBLE PRECISION TEMP, C
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (QMIN(L).GT.QMIN(J2)) L=J2
10       CONTINUE
         TEMP=QMIN(L)
         QMIN(L)=QMIN(J1)
         QMIN(J1)=TEMP
         NTEMP=FF(L)
         FF(L)=FF(J1)
         FF(J1)=NTEMP
         DO J2=1,3*NATOMS
            C=QMINP(L,J2)
            QMINP(L,J2)=QMINP(J1,J2)
            QMINP(J1,J2)=C
         ENDDO
      ENDDO

      RETURN
! }}}
      END

! doxygen - GET_RND {{{
!> @name         GET_RND
! 
!> @brief        Get an array of random numbers inside the interval [XMIN,XMAX]
!
!> @param[in]    N              - dimension of the array RND
!> @param[out]   RND            - the generated array of random numbers 
!> @param[in]    XMIN,XMAX      
!
! }}}
SUBROUTINE GET_RND(RND,N,XMIN,XMAX)
! {{{
IMPLICIT NONE

! random number vector
DOUBLE PRECISION :: RND(:),XMIN,XMAX,DX
! dimension of RND(:)
INTEGER N,I

DX=XMAX-XMIN

DO I=1,N 
        RND(I)=XMAX-DX*DPRAND()
ENDDO

RETURN
! }}}
END 

     SUBROUTINE SDPRND (ISEED)
        ! declarations {{{
        DOUBLE PRECISION XMOD, YMOD, POLY(101), OTHER, OFFSET, X
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0)
        INTEGER ISEED, INDEX, IX, IY, IZ, I
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
        ! }}}
        ! subroutine body {{{
!
!   ISEED should be set to an integer between 0 and 9999 inclusive;
!   a value of 0 will initialise the generator only if it has not
!   already been done.
!
        IF (INITAL .OR. ISEED .NE. 0) THEN
            INITAL = .FALSE.
        ELSE
            RETURN
        END IF
!
!   INDEX must be initialised to an integer between 1 and 101
!   inclusive, POLY(1...N) to integers between 0 and 1000009710
!   inclusive (not all 0), and OTHER to a non-negative proper fraction
!   with denominator 33554432.  It uses the Wichmann-Hill generator to
!   do this.
!
        IX = MOD(ABS(ISEED),10000)+1
        IY = 2*IX+1
        IZ = 3*IX+1
        DO 10 I = -10,101
            IF (I .GE. 1) POLY(I) = AINT(XMOD*X)
            IX = MOD(171*IX,30269)
            IY = MOD(172*IY,30307)
            IZ = MOD(170*IZ,30323)
            X = MOD(DBLE(IX)/30269.0D0+DBLE(IY)/30307.0D0+
     1        DBLE(IZ)/30323.0D0,1.0D0)
  10    CONTINUE
        OTHER = AINT(YMOD*X)/YMOD
        OFFSET = 1.0D0/YMOD
        INDEX = 1
        ! }}}
        END

FUNCTION DPRAND()
        ! DECLARATIONS {{{
        DOUBLE PRECISION XMOD, YMOD, XMOD2, XMOD4, TINY, POLY(101), DPRAND,
     1    OTHER, OFFSET, X, Y
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0,
     1    XMOD2 = 2000019422.0D0, XMOD4 = 4000038844.0D0,
     2    TINY = 1.0D-17)
        INTEGER INDEX, N
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
        ! }}}
        ! SUBROUTINE BODY {{{
C
C   THIS RETURNS A UNIFORM (0,1) RANDOM NUMBER, WITH EXTREMELY GOOD
C   UNIFORMITY PROPERTIES.  IT ASSUMES THAT DOUBLE PRECISION PROVIDES
C   AT LEAST 33 BITS OF ACCURACY, AND USES A POWER OF TWO BASE.
C
        IF (INITAL) THEN
            CALL SDPRND (0)
            INITAL = .FALSE.
        END IF
C
C   SEE [KNUTH] FOR WHY THIS IMPLEMENTS THE ALGORITHM DESCRIBED IN
C   THE PAPER.  NOTE THAT THIS CODE IS TUNED FOR MACHINES WITH FAST
C   DOUBLE PRECISION, BUT SLOW MULTIPLY AND DIVIDE; MANY, MANY OTHER
C   OPTIONS ARE POSSIBLE.
C
        N = INDEX-64
        IF (N .LE. 0) N = N+101
        X = POLY(INDEX)+POLY(INDEX)
        X = XMOD4-POLY(N)-POLY(N)-X-X-POLY(INDEX)
        IF (X .LT. 0.0D0) THEN
            IF (X .LT. -XMOD) X = X+XMOD2
            IF (X .LT. 0.0D0) X = X+XMOD
        ELSE
            IF (X .GE. XMOD2) THEN
                X = X-XMOD2
                IF (X .GE. XMOD) X = X-XMOD
            END IF
            IF (X .GE. XMOD) X = X-XMOD
        END IF
        POLY(INDEX) = X
        INDEX = INDEX+1
        IF (INDEX .GT. 101) INDEX = INDEX-101
C
C   ADD IN THE SECOND GENERATOR MODULO 1, AND FORCE TO BE NON-ZERO.
C   THE RESTRICTED RANGES LARGELY CANCEL THEMSELVES OUT.
C
   10   Y = 37.0D0*OTHER+OFFSET
        OTHER = Y-AINT(Y)
        IF (OTHER .EQ. 0.0D0) GO TO 10
        X = X/XMOD+OTHER
        IF (X .GE. 1.0D0) X = X-1.0D0
        DPRAND = X+TINY
        ! }}}
        END

! input {{{

      SUBROUTINE INPUT(END)
      ! {{{
      USE PORFUNCS
 
!  Read next input record from unit IR into buffer, and analyse for data
!  items, null items (enclosed by commas), and comment (enclosed in
!  parentheses, which may be nested, and occurring where spaces may
!  occur).

!  E1  Attempts to handle stream-switching commands in the data.
!      #include file-name
!         switches input to be from the file specified;
!      #revert
!         (or end-of-file) reverts to the previous file.
!    Also
!      #concat "string"
!         sets the concatenation string; e.g.
!         #concat "\"
 
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      LOGICAL BLANK, TCOMMA, END
      CHARACTER(LEN=8) CONCAT
      DATA CONCAT /'+++'/, LC /3/
      SAVE CONCAT, LC
 
!  ITEM    is the number of the last item read from the buffer
 
!  NITEM   is the number of items in the buffer
 
!  CHAR(I) is the Ith character in the buffer
 
!  LOC(I)  is the starting position of the Ith item
 
!  LINE    is the number of lines (records) read
 
!  If SKIPBL is set to TRUE, lines containing no items other than
!  comment are skipped automatically, so that INPUT will always return
!  a line with at least one item (possibly null).  If SKIPBL is FALSE
!  (default) no skipping occurs and the next data line is returned
!  regardless of what it contains.
 
!  If CLEAR is set to TRUE (default) then an attempt to read more than
!  NITEM items from a line will cause zeros or blanks to be returned. If
!  CLEAR is FALSE, the variables in question are left unaltered.
 
!  NCR is the number of single characters which have been read from the
!  current item.  It is always zero unless READCH has been used; if non-
!  zero, ITEM refers to the current item from which characters are being
!  read.
 
!  NERROR specifies the treatment of errors found while reading numbers:
!    0  Hard error - print message and stop (default).
!    1  Soft error - print message and return zero.
!    2  Soft error - no message, return zero, set NERROR to -1.
!  If NERROR is set to 2 it is important to test and reset it after reading
!  a number.
 
!  IR is the input stream from which the data are to be read (default 5).
 
!  If ECHO is TRUE, the input line will be reflected to the system output
!  stream.
 
!  LAST gives the position of the last non-blank character in the line.
 
!  CONCAT is the line concatenation string: if it is found as the last
!  non-blank character string on a line, the following line is read and
!  appended to the current line, overwriting the concatenation string.
!  This procedure is repeated if necessary until a line is found that does
!  not end with the concatenation string.
 
!  The concatenation string may be modified by changing the DATA statement
!  above. The integer LC must be set to the number of non-blank characters
!  in the concatenation string. If it is set to zero, no concatenation
!  occurs. The concatenation string may also be changed by the #concat
!  directive in the data file.
 
      INTEGER LEVEL
      SAVE LEVEL, IR0
      CHARACTER(LEN=132) W, F

      CHARACTER SPACE, BRA, KET, COMMA, SQUOTE, DQUOTE, TERM
      DATA LEVEL /0/
      DATA SPACE /' '/, BRA /'('/, KET /')'/, COMMA /','/,
     *    SQUOTE /''''/, DQUOTE /'"'/
 
      ECHO=.FALSE.

      END=.FALSE.
      CAT=.FALSE.
10    M=1
20    LAST=M+131
      READ (IR,1001,END=900) (CHAR(I), I=M,LAST)
      IF (ECHO) PRINT '(1X,132A1)', (CHAR(I), I=M,LAST)
1001  FORMAT (132A1)
      CAT=.FALSE.
      LINE=LINE+1
!  Find last non-blank character
30    IF (CHAR(LAST) .EQ. SPACE) THEN
        LAST=LAST-1
        IF (LAST .GT. 1) GO TO 30
      ENDIF
!  Look for concatenation string
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

!  Logical line assembled. First look for input directives
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
      CALL UPPERCASE(W)
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

!  Analyse input 
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
 
!  Reading through an item, seeking terminator.
!  (Next line of code suppressed: comment not allowed within an item)
!     IF (CHAR(L) .EQ. BRA) GO TO 400
      TCOMMA=CHAR(L) .EQ. COMMA .AND. TERM .EQ. SPACE
      BLANK=TCOMMA .OR. CHAR(L) .EQ. TERM
      GO TO 100
 
!  Looking for next item
200   TERM=SPACE
      BLANK=CHAR(L) .EQ. SPACE
      IF (BLANK) GO TO 100
      IF (CHAR(L) .EQ. BRA) GO TO 400
 
!  Item found
      NITEM=NITEM+1
      LOC(NITEM)=L
 
!  Quoted string?
      IF (CHAR(L) .EQ. SQUOTE .OR. CHAR(L) .EQ. DQUOTE) TERM=CHAR(L)
 
!  Null item?
      IF (CHAR(L) .NE. COMMA) GO TO 100
      IF (TCOMMA) LOC(NITEM)=0
      GO TO 80
 
!  Comment (enclosed in parentheses, possibly nested)
400   NEST=1
410   L=L+1
      IF (L .GT. LAST) GO TO 800
      IF (CHAR(L) .EQ. BRA) NEST=NEST+1
      IF (CHAR(L) .EQ. KET) NEST=NEST-1
      IF (NEST .GT. 0) GO TO 410
      GO TO 90
 
!  End of card -- was it blank?
800   IF (SKIPBL .AND. NITEM .EQ. 0) GO TO 10
      RETURN
 
!  End of data
900   IF (CAT) THEN
        PRINT '(A)', 'Apparently concatenating at end-of-file'
        CALL REPORT('Unexpected end of data file',.TRUE.)
      ENDIF
      IF (LEVEL .GT. 0) THEN
!  Revert to previous input
        CLOSE(IR)
        LEVEL=LEVEL-1
        IF (LEVEL .EQ. 0) THEN
          IR=IR0
        ELSE
          IR=90+LEVEL
        ENDIF
        GO TO 10
      ENDIF
      DO L=1,LAST
         CHAR(L)=SPACE
      ENDDO
      ITEM=0
      NITEM=-1
      END=.TRUE.
      RETURN
! }}} 
      END
 
      BLOCK DATA INBLK
      ! {{{

      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      DATA ITEM /0/, NITEM /0/, CHAR /800*' '/, LOC /132*0/, LINE /0/
      DATA SKIPBL /.FALSE./, CLEAR /.TRUE./, NCR /0/, NERROR /0/
      DATA IR /5/, ECHO /.FALSE./
      ! }}}
      END
 
      SUBROUTINE READF(A)
      ! {{{
      use porfuncs
 
!  Read the next item from the buffer as a real number
 
      INTEGER P, STATE
      LOGICAL NULL, XNULL
      DOUBLE PRECISION A, AD, B, TEN, SIGN
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      CHARACTER C, PLUS, MINUS, DOT, SPACE, STAR, D, E, COMMA,
     &    S1, S2, DIGIT(10)
      DATA DIGIT /'1', '2', '3', '4', '5', '6', '7', '8', '9', '0'/
      DATA PLUS /'+'/, MINUS /'-'/, DOT /'.'/, SPACE /' '/, STAR /'*'/,
     *    D /'D'/, E /'E'/, COMMA /','/
      DATA TEN /10D0/
 
      IF (CLEAR) A=0D0
 
!  If the item is null, or if no item remains, A is unchanged
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
 
!  Terminators
12    C=CHAR(P)
      IF (C .EQ. SPACE .OR. C .EQ. COMMA) GO TO 100
 
!  Digits
      DO 20 I=1,10
      IF (C .EQ. DIGIT(I)) GO TO (70,72,71,80,81), STATE
20    CONTINUE
 
!  State numbers:  1  nothing read yet
!                  2  reading digits of integral part
!                  3  reading decimal fraction
!                  4  expecting exponent
!                  5  reading digits of exponent
 
!  Number control characters -- branch on state
 
      IF (C .EQ. DOT) GO TO (30,30,99,99,99), STATE
      IF (C .EQ. PLUS) GO TO (42,52,52,52,99), STATE
      IF (C .EQ. MINUS) GO TO (41,51,51,51,99), STATE
      IF (C .EQ. D .OR. C .EQ. E) GO TO (99,60,60,99,99), STATE
 
!  Error
99    IF (NERROR .LE. 1) THEN
        WRITE (6,'(1X,A,A,A)') 'Invalid character "', C,
     &                             '" while reading number.'
        I2=MIN(LAST,P+20)
        I1=MAX(1,I2-70)
        S1=' '
        IF (I1 .GT. 1) S1='...'
        S2=' '
        IF (I2 .LT. LAST) S2='...'
        PRINT '(1X,A/1X,A,1X,132A,1X,A)',
     &           'Current input line is:', S1, (CHAR(I), I=I1,I2), S2
        PRINT '(5X,132A1)', (' ', I=I1,P), '*'
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
! }}} 
      END
 
      SUBROUTINE INPUTI(I)
      ! {{{
      DOUBLE PRECISION A
      A=I
      CALL READF(A)
      I=A
      RETURN
      ! }}}
      END
 
      SUBROUTINE READK(C,M,N,K)
! {{{ 
!  The arrays C, M and N comprise a tree, each node of which consists of
!  a character C(I) and two integers M(I) and N(I). Initially I=1. The
!  character C(I) is tested against the current character in the input.
!  If it matches, the next input character is selected and M(I) is
!  examined; otherwise N(I) is examined. If a positive value is found,
!  it is a pointer to the next node to be tested; if negative or zero,
!  it is minus the key value to be returned in K. If a keyword is
!  present in the input but does not match a keyword in the tree, then
!  zero is returned; if a null item is present or if there are no more
!  items on the line, then -1 is returned. Note that this arrangement
!  allows the code in the tree to be set up in many ways to suit the
!  requirements of the problem. Abbreviations may be defined
!  economically, and the keyword can be required to terminate with a
!  space or other symbol, or it can be self-terminating.
 
!  Example:  PRINT and PRT are to return key 19, and must be terminated
!  by space.  P alone (i.e. not present as the first letter of PRINT or
!  PRT) is to return 7.  Initial spaces are to be ignored.  The code
!  table might begin
 
!        1   2   3   4   5   6   7   8
 
!   C   ' '  P   R   I   N   T  ' '  ...
!   M    1   3   4   5   6   7  -19  ...
!   N    2   8  -7   6   0   0   0   ...
 
!  Such code tables can be created automatically by the procedure
!  ENCODE, which produces the DATA statements required to initialize
!  the array.
 
 
 
      INTEGER M(*), N(*), TP, P
      CHARACTER C(*), SPACE
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(132), LINE, SKIPBL, CLEAR, NCR,
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
 
!  Advance character pointer
10    P=P+1
20    IF (P .LE. LAST) GO TO 30
!  End of line
      IF (TP .EQ. 1) K=-1
      RETURN
 
30    CALL UPPERCASE(CHAR(P))
40    IF(CHAR(P) .EQ. C(TP)) GO TO 50
!  No match -- advance tree pointer
      TP=N(TP)
!  Zero pointer -- undecodeable
      IF (TP .EQ. 0) RETURN
!  Positive value -- test same input character again
      IF (TP .GT. 0) GO TO 40
!  Negative fail pointer: the input contains a keyword which is an initial
!  substring of some other keyword.
      GO TO 70
 
!  Matched -- advance tree pointer to next character of command
50    TP=M(TP)
      IF (TP .GT. 0) GO TO 10
 
!  Last character of keyword found.
70    K=-TP
      RETURN
      ! }}}
      END
 
      SUBROUTINE READA(M)
      ! {{{
      IMPLICIT NONE
 
!  Read characters and pack them into the character variable M.
!  If the first character is a single or double quote, the string is
!  terminated by a matching quote and the quotes are removed; otherwise
!  it is terminated by space or comma.  If necessary, M is filled
!  out with spaces.
 
      CHARACTER SPACE, BLANK, COMMA, SQUOTE, DQUOTE, TERM
      CHARACTER*(*) M
      CHARACTER CHAR
      INTEGER I, K, NCR, NERROR, ITEM, NITEM, LOC, IR, LAST, L, LINE
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(132), LINE, SKIPBL, CLEAR, NCR,
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
!
!  The new Sun compiler refuses to execute the following line - hence rewritten
!
!     K=K+1
!
11    IF (L+K .GT. LAST) GO TO 50

!     IF (.NOT. QUOTE .AND. CHAR(L+K) .NE. SPACE
!    *                .AND. CHAR(L+K) .NE. COMMA) GO TO 10
!     IF (QUOTE .AND. CHAR(L+K) .NE. TERM) GO TO 10

      IF ((.NOT. QUOTE .AND. CHAR(L+K) .NE. SPACE
     *                .AND. CHAR(L+K) .NE. COMMA).OR.
     *    (QUOTE .AND. CHAR(L+K) .NE. TERM)) THEN
         K=K+1
         GOTO 10
      ENDIF
 
50    IF (K .EQ. 0) RETURN
      K=MIN0(K,LEN(M))
      DO I=1,K
         M(I:I)=CHAR(L+I-1)
      ENDDO
65    DO I=K+1,LEN(M)
         M(I:I)=SPACE
      ENDDO
      RETURN
      ! }}}
      END
 
      SUBROUTINE READU(M)
      ! {{{
      CHARACTER*(*) M
 
      CALL READA(M)
      CALL UPPERCASE(M)
      RETURN

      ENTRY READL(M)
      CALL READA(M)
      CALL LOCASE(M)
      RETURN
      ! }}}
      END
 
      SUBROUTINE READCH(M)
! {{{ 
!  Read a single character from the next (or the current) data item.
!  No account is taken of special characters.
 
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      CHARACTER*(*) M
 
      M=' '
      IF (ITEM .GE. NITEM .AND. NCR .EQ. 0) RETURN
      IF (NCR .EQ. 0) ITEM=ITEM+1
      L=LOC(ITEM)
      M=CHAR(L+NCR)
      NCR=NCR+1
      RETURN
      ! }}}
      END
 
      SUBROUTINE GETF(A)
      ! {{{
      USE PORFUNCS
!  Read the next item as a double-precision number, reading new data
!  records if necessary
      DOUBLE PRECISION A

      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      LOGICAL END
 
10    IF (ITEM .LT. NITEM) THEN
        CALL READF(A)
        RETURN
      ENDIF
      CALL INPUT(END)
      IF (END) THEN
        WRITE (6,1001)
1001    FORMAT ('0End of file while attempting to read a number')
        STOP
      ENDIF
      GO TO 10
 
      ! }}}
      END
 
      SUBROUTINE GETS(S)
      ! {{{
!  Get a single-precision number
      DOUBLE PRECISION A
      REAL S
      A=S
      CALL GETF(A)
      S=A
      RETURN
      ! }}}
      END
 
      SUBROUTINE GETI(I)
      ! {{{
!  Get an integer
      DOUBLE PRECISION A
      A=I
      CALL GETF(A)
      I=A
      RETURN
      ! }}}
      END
 
      SUBROUTINE GETA(M)
      ! {{{
      USE PORFUNCS
!  Get a character string
      CHARACTER*(*) M

      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      LOGICAL END

10    IF (ITEM .LT. NITEM) THEN
        CALL READA(M)
        RETURN
      ENDIF
      CALL INPUT(END)
      IF (END) THEN
        WRITE (6,1001)
1001    FORMAT
     &    ('0End of file while attempting to read a character string')
        STOP
      ENDIF
      GO TO 10
      ! }}} 
      END
 
      SUBROUTINE READSNGL(S)
      ! {{{
!  Read a number from the current record into the single-precision
!  variable S
      DOUBLE PRECISION A
      REAL S
      A=S
      CALL READF(A)
      S=A
      RETURN
      ! }}}
      END
 
      SUBROUTINE READI(I)
      ! {{{
!  Read an integer from the current record
      DOUBLE PRECISION A
      A=I
      CALL READF(A)
      I=A
      RETURN
      ! }}}
      END
 
      SUBROUTINE REREAD(K)
      ! {{{ 
!  K>0  Reread from item K
!  K<0  Go back |K| items
!  K=0  Same as K=-1, i.e. reread last item.
 
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
 
      IF (K .LT. 0) ITEM=ITEM+K
      IF (K .EQ. 0) ITEM=ITEM-1
      IF (K .GT. 0) ITEM=K-1
      IF (ITEM .LT. 0) ITEM=0
      NCR=0
      RETURN
      ! }}}
      END
 
      SUBROUTINE UPPERCASE(WORD)
      ! {{{
      CHARACTER WORD*(*)
 
      CHARACTER UC*26, LC*26
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
      ! }}} 
      END

      SUBROUTINE REPORT(C,REFLCT)
      ! declarations {{{

      CHARACTER C*(*)
      LOGICAL REFLCT

      CHARACTER BUFF
      COMMON /BUFFER/ BUFF(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEMS, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT

      INTEGER PRINT
      LOGICAL SWITCH, ONLINE
      COMMON /TEST/ SWITCH(8), PRINT
      EQUIVALENCE (ONLINE, SWITCH(4))

      CHARACTER(LEN=3) S1, S2
      ! }}}
      ! subroutine body {{{

      PRINT '(1X,A)', C
      IF (REFLCT) THEN
        L=LOC(ITEM)
        I2=MIN(LAST,L+20)
        I1=MAX(1,I2-70)
        S1=' '
        IF (I1 .GT. 1) S1='...'
        S2=' '
        IF (I2 .LT. LAST) S2='...'
!       PRINT '(1X,4(A,I3))',
!    &      'L =', L, '   LAST =', LAST, '   I1 =', I1, '   I2 =', I2
        PRINT '(1X,A/1X,A,1X,132A,1X,A)',
     &           'Current input line is:', S1, (BUFF(I), I=I1,I2), S2
        PRINT '(4X,132A1)', (' ', I=I1,L), '*'
      ENDIF
!     IF (.NOT. ONLINE) STOP
      ! }}}
      END
! }}}

! dae 
! get coordinates from external file 'coords'
! which is CHARMM format but use standard fortran
! reading commands, unlike CHARMM which has a limit
! on total line length and does involved procedure with strings
!
      SUBROUTINE READREF(FNAME)
      ! {{{

      USE UTILS
      USE COMMONS

      IMPLICIT NONE
  
      CHARACTER(LEN=*) FNAME
      INTEGER I,I1
 
!     CHARACTER(LEN=4) RESLABEL(NATOMS),ATOMLABEL(NATOMS)
!     INTEGER RESNUMBER(NATOMS)
!     COMMON /CHARMMLABEL/ RESLABEL, ATOMLABEL, RESNUMBER

      OPEN(UNIT=19,FILE=FNAME,STATUS='OLD',IOSTAT=I)
      CALL OPENIOSTAT(I,FNAME,'READREF')
      READ (19,*)
      DO I=1,NATOMS
         READ (19,*) I1,RESNUMBER(I),RESLABEL(I),ATOMLABEL(I)
      ENDDO
  
      CLOSE(19)

      RETURN
      ! }}}
      END

     SUBROUTINE OPENIOSTAT(I,FNAME,SUBRID)
     ! {{{

          IMPLICIT NONE
          
          INTEGER,INTENT(IN) :: I
          CHARACTER(LEN=*),INTENT(IN) :: FNAME
          CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: SUBRID

          IF (I/=0) THEN
               PRINT *, 'ERROR WHILE OPENING FILE "'//TRIM(ADJUSTL(FNAME))//'"'
               IF (PRESENT(SUBRID)) PRINT *, 'STOP IN '//TRIM(ADJUSTL(SUBRID))
               STOP
          ENDIF
          ! }}}
     END SUBROUTINE OPENIOSTAT


ENDMODULE





