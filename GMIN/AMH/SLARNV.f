      SUBROUTINE SLARNV( IDIST, ISEED, N, X )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     SEPTEMBER 30, 1994
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            IDIST, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            ISEED( 4 )
      REAL               X( N )
*     ..
*
*  PURPOSE
*  =======
*
*  SLARNV RETURNS A VECTOR OF N RANDOM DOUBLE PRECISION NUMBERS FROM A UNIFORM OR
*  NORMAL DISTRIBUTION.
*
*  ARGUMENTS
*  =========
*
*  IDIST   (INPUT) INTEGER
*          SPECIFIES THE DISTRIBUTION OF THE RANDOM NUMBERS:
*          = 1:  UNIFORM (0,1)
*          = 2:  UNIFORM (-1,1)
*          = 3:  NORMAL (0,1)
*
*  ISEED   (INPUT/OUTPUT) INTEGER ARRAY, DIMENSION (4)
*          ON ENTRY, THE SEED OF THE RANDOM NUMBER GENERATOR; THE ARRAY
*          ELEMENTS MUST BE BETWEEN 0 AND 4095, AND ISEED(4) MUST BE
*          ODD.
*          ON EXIT, THE SEED IS UPDATED.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF RANDOM NUMBERS TO BE GENERATED.
*
*  X       (OUTPUT) REAL ARRAY, DIMENSION (N)
*          THE GENERATED RANDOM NUMBERS.
*
*  FURTHER DETAILS
*  ===============
*
*  THIS ROUTINE CALLS THE AUXILIARY ROUTINE SLARUV TO GENERATE RANDOM
*  DOUBLE PRECISION NUMBERS FROM A UNIFORM (0,1) DISTRIBUTION, IN BATCHES OF UP TO
*  128 USING VECTORISABLE CODE. THE BOX-MULLER METHOD IS USED TO
*  TRANSFORM NUMBERS FROM A UNIFORM TO A NORMAL DISTRIBUTION.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      REAL               ONE, TWO
      PARAMETER          ( ONE = 1.0E+0, TWO = 2.0E+0 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      REAL               TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663E+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, IL, IL2, IV
*     ..
*     .. LOCAL ARRAYS ..
      REAL                 U( LV )
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          COS, LOG, MIN, SQRT
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           SLARUV
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      DO 40 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         IF( IDIST.EQ.3 ) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF
*
*        CALL SLARUV TO GENERATE IL2 NUMBERS FROM A UNIFORM (0,1)
*        DISTRIBUTION (IL2 <= LV)
*
         CALL SLARUV( ISEED, IL2, U )
*
         IF( IDIST.EQ.1 ) THEN
*
*           COPY GENERATED NUMBERS
*
            DO 10 I = 1, IL
               X( IV+I-1 ) = U( I )
   10       CONTINUE
         ELSE IF( IDIST.EQ.2 ) THEN
*
*           CONVERT GENERATED NUMBERS TO UNIFORM (-1,1) DISTRIBUTION
*
            DO 20 I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
   20       CONTINUE
         ELSE IF( IDIST.EQ.3 ) THEN
*
*           CONVERT GENERATED NUMBERS TO NORMAL (0,1) DISTRIBUTION
*
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )*
     $                       COS( TWOPI*U( 2*I ) )
   30       CONTINUE
         END IF
   40 CONTINUE
      RETURN
*
*     END OF SLARNV
*
      END
