C   COPYRIGHT (C) 1992  N.M. MACLAREN
C   COPYRIGHT (C) 1992  THE UNIVERSITY OF CAMBRIDGE

C   THIS SOFTWARE MAY BE REPRODUCED AND USED FREELY, PROVIDED THAT ALL
C   USERS OF IT AGREE THAT THE COPYRIGHT HOLDERS ARE NOT LIABLE FOR ANY
C   DAMAGE OR INJURY CAUSED BY USE OF THIS SOFTWARE AND THAT THIS
C   CONDITION IS PASSED ONTO ALL SUBSEQUENT RECIPIENTS OF THE SOFTWARE,
C   WHETHER MODIFIED OR NOT.



        SUBROUTINE SDPRND (ISEED)
        DOUBLE PRECISION XMOD, YMOD, POLY(101), OTHER, OFFSET, X
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0)
        INTEGER ISEED, INDEX, IX, IY, IZ, I
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
C
C   ISEED SHOULD BE SET TO AN INTEGER BETWEEN 0 AND 9999 INCLUSIVE;
C   A VALUE OF 0 WILL INITIALISE THE GENERATOR ONLY IF IT HAS NOT
C   ALREADY BEEN DONE.
C
        IF (INITAL .OR. ISEED .NE. 0) THEN
            INITAL = .FALSE.
        ELSE
            RETURN
        END IF
C
C   INDEX MUST BE INITIALISED TO AN INTEGER BETWEEN 1 AND 101
C   INCLUSIVE, POLY(1...N) TO INTEGERS BETWEEN 0 AND 1000009710
C   INCLUSIVE (NOT ALL 0), AND OTHER TO A NON-NEGATIVE PROPER FRACTION
C   WITH DENOMINATOR 33554432.  IT USES THE WICHMANN-HILL GENERATOR TO
C   DO THIS.
C
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
        END

        DOUBLE PRECISION FUNCTION DPRAND()
        DOUBLE PRECISION XMOD, YMOD, XMOD2, XMOD4, TINY, POLY(101),
     1    OTHER, OFFSET, X, Y
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0,
     1    XMOD2 = 2000019422.0D0, XMOD4 = 4000038844.0D0,
     2    TINY = 1.0D-17)
        INTEGER INDEX, N
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
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
        END
