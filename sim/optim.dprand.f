
C   Copyright (C) 1992  N.M. Maclaren
C   Copyright (C) 1992  The University of Cambridge

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
C
C   ISEED should be set to an integer between 0 and 9999 inclusive;
C   a value of 0 will initialise the generator only if it has not
C   already been done.
C
        IF (INITAL .OR. ISEED .NE. 0) THEN
            INITAL = .FALSE.
        ELSE
            RETURN
        END IF
C
C   INDEX must be initialised to an integer between 1 and 101
C   inclusive, POLY(1...N) to integers between 0 and 1000009710
C   inclusive (not all 0), and OTHER to a non-negative proper fraction
C   with denominator 33554432.  It uses the Wichmann-Hill generator to
C   do this.
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
        ! }}}
        END

  
