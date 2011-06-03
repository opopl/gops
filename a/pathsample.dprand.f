
        FUNCTION DPRAND()
        ! declarations {{{
        DOUBLE PRECISION XMOD, YMOD, XMOD2, XMOD4, TINY, POLY(101),
     1    OTHER, OFFSET, X, Y, DPRAND
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0,
     1    XMOD2 = 2000019422.0D0, XMOD4 = 4000038844.0D0,
     2    TINY = 1.0D-17)
        INTEGER INDEX, N
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
        ! }}}
        ! subroutine body {{{
C
C   This returns a uniform (0,1) random number, with extremely good
C   uniformity properties.  It assumes that double precision provides
C   at least 33 bits of accuracy, and uses a power of two base.
C
        IF (INITAL) THEN
            CALL SDPRND (0)
            INITAL = .FALSE.
        END IF
C
C   See [Knuth] for why this implements the algorithm described in
C   the paper.  Note that this code is tuned for machines with fast
C   double precision, but slow multiply and divide; many, many other
C   options are possible.
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
C   Add in the second generator modulo 1, and force to be non-zero.
C   The restricted ranges largely cancel themselves out.
C
   10   Y = 37.0D0*OTHER+OFFSET
        OTHER = Y-AINT(Y)
        IF (OTHER .EQ. 0.0D0) GO TO 10
        X = X/XMOD+OTHER
        IF (X .GE. 1.0D0) X = X-1.0D0
        DPRAND = X+TINY
        ! }}}
        END
