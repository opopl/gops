      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1999
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGETRI COMPUTES THE INVERSE OF A MATRIX USING THE LU FACTORIZATION
*  COMPUTED BY DGETRF.
*
*  THIS METHOD INVERTS U AND THEN COMPUTES INV(A) BY SOLVING THE SYSTEM
*  INV(A)*L = INV(U) FOR INV(A).
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE FACTORS L AND U FROM THE FACTORIZATION
*          A = P*L*U AS COMPUTED BY DGETRF.
*          ON EXIT, IF INFO = 0, THE INVERSE OF THE ORIGINAL MATRIX A.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  IPIV    (INPUT) INTEGER ARRAY, DIMENSION (N)
*          THE PIVOT INDICES FROM DGETRF; FOR 1<=I<=N, ROW I OF THE
*          MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  WORK    (WORKSPACE/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          ON EXIT, IF INFO=0, THEN WORK(1) RETURNS THE OPTIMAL LWORK.
*
*  LWORK   (INPUT) INTEGER
*          THE DIMENSION OF THE ARRAY WORK.  LWORK >= MAX(1,N).
*          FOR OPTIMAL PERFORMANCE LWORK >= N*NB, WHERE NB IS
*          THE OPTIMAL BLOCKSIZE RETURNED BY ILAENV.
*
*          IF LWORK = -1, THEN A WORKSPACE QUERY IS ASSUMED; THE ROUTINE
*          ONLY CALCULATES THE OPTIMAL SIZE OF THE WORK ARRAY, RETURNS
*          THIS VALUE AS THE FIRST ENTRY OF THE WORK ARRAY, AND NO ERROR
*          MESSAGE RELATED TO LWORK IS ISSUED BY XERBLA.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0:  IF INFO = I, U(I,I) IS EXACTLY ZERO; THE MATRIX IS
*                SINGULAR AND ITS INVERSE COULD NOT BE COMPUTED.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 )
     $   RETURN
*
*     FORM INV(U).  IF INFO > 0 FROM DTRTRI, THEN U IS SINGULAR,
*     AND THE INVERSE IS NOT COMPUTED.
*
      CALL DTRTRI( 'UPPER', 'NON-UNIT', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     SOLVE THE EQUATION INV(A)*L = INV(U) FOR INV(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        USE UNBLOCKED CODE.
*
         DO 20 J = N, 1, -1
*
*           COPY CURRENT COLUMN OF L TO WORK AND REPLACE WITH ZEROS.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           COMPUTE CURRENT COLUMN OF INV(A).
*
            IF( J.LT.N )
     $         CALL DGEMV( 'NO TRANSPOSE', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        USE BLOCKED CODE.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           COPY CURRENT BLOCK COLUMN OF L TO WORK AND REPLACE WITH
*           ZEROS.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           COMPUTE CURRENT BLOCK COLUMN OF INV(A).
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL DTRSM( 'RIGHT', 'LOWER', 'NO TRANSPOSE', 'UNIT', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     APPLY COLUMN INTERCHANGES.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     END OF DGETRI
*
      END
      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DTRTI2 COMPUTES THE INVERSE OF A REAL UPPER OR LOWER TRIANGULAR
*  MATRIX.
*
*  THIS IS THE LEVEL 2 BLAS VERSION OF THE ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          SPECIFIES WHETHER THE MATRIX A IS UPPER OR LOWER TRIANGULAR.
*          = 'U':  UPPER TRIANGULAR
*          = 'L':  LOWER TRIANGULAR
*
*  DIAG    (INPUT) CHARACTER*1
*          SPECIFIES WHETHER OR NOT THE MATRIX A IS UNIT TRIANGULAR.
*          = 'N':  NON-UNIT TRIANGULAR
*          = 'U':  UNIT TRIANGULAR
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE TRIANGULAR MATRIX A.  IF UPLO = 'U', THE
*          LEADING N BY N UPPER TRIANGULAR PART OF THE ARRAY A CONTAINS
*          THE UPPER TRIANGULAR MATRIX, AND THE STRICTLY LOWER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF UPLO = 'L', THE
*          LEADING N BY N LOWER TRIANGULAR PART OF THE ARRAY A CONTAINS
*          THE LOWER TRIANGULAR MATRIX, AND THE STRICTLY UPPER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF DIAG = 'U', THE
*          DIAGONAL ELEMENTS OF A ARE ALSO NOT REFERENCED AND ARE
*          ASSUMED TO BE 1.
*
*          ON EXIT, THE (TRIANGULAR) INVERSE OF THE ORIGINAL MATRIX, IN
*          THE SAME STORAGE FORMAT.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DSCAL, DTRMV, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        COMPUTE INVERSE OF UPPER TRIANGULAR MATRIX.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           COMPUTE ELEMENTS 1:J-1 OF J-TH COLUMN.
*
            CALL DTRMV( 'UPPER', 'NO TRANSPOSE', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        COMPUTE INVERSE OF LOWER TRIANGULAR MATRIX.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              COMPUTE ELEMENTS J+1:N OF J-TH COLUMN.
*
               CALL DTRMV( 'LOWER', 'NO TRANSPOSE', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     END OF DTRTI2
*
      END
      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DTRTRI COMPUTES THE INVERSE OF A REAL UPPER OR LOWER TRIANGULAR
*  MATRIX A.
*
*  THIS IS THE LEVEL 3 BLAS VERSION OF THE ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          = 'U':  A IS UPPER TRIANGULAR;
*          = 'L':  A IS LOWER TRIANGULAR.
*
*  DIAG    (INPUT) CHARACTER*1
*          = 'N':  A IS NON-UNIT TRIANGULAR;
*          = 'U':  A IS UNIT TRIANGULAR.
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE TRIANGULAR MATRIX A.  IF UPLO = 'U', THE
*          LEADING N-BY-N UPPER TRIANGULAR PART OF THE ARRAY A CONTAINS
*          THE UPPER TRIANGULAR MATRIX, AND THE STRICTLY LOWER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF UPLO = 'L', THE
*          LEADING N-BY-N LOWER TRIANGULAR PART OF THE ARRAY A CONTAINS
*          THE LOWER TRIANGULAR MATRIX, AND THE STRICTLY UPPER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF DIAG = 'U', THE
*          DIAGONAL ELEMENTS OF A ARE ALSO NOT REFERENCED AND ARE
*          ASSUMED TO BE 1.
*          ON EXIT, THE (TRIANGULAR) INVERSE OF THE ORIGINAL MATRIX, IN
*          THE SAME STORAGE FORMAT.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0: IF INFO = I, A(I,I) IS EXACTLY ZERO.  THE TRIANGULAR
*               MATRIX IS SINGULAR AND ITS INVERSE CAN NOT BE COMPUTED.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRI', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 )
     $   RETURN
*
*     CHECK FOR SINGULARITY IF NON-UNIT.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     DETERMINE THE BLOCK SIZE FOR THIS ENVIRONMENT.
*
      NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        USE UNBLOCKED CODE
*
         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        USE BLOCKED CODE
*
         IF( UPPER ) THEN
*
*           COMPUTE INVERSE OF UPPER TRIANGULAR MATRIX
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              COMPUTE ROWS 1:J-1 OF CURRENT BLOCK COLUMN
*
               CALL DTRMM( 'LEFT', 'UPPER', 'NO TRANSPOSE', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL DTRSM( 'RIGHT', 'UPPER', 'NO TRANSPOSE', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              COMPUTE INVERSE OF CURRENT DIAGONAL BLOCK
*
               CALL DTRTI2( 'UPPER', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           COMPUTE INVERSE OF LOWER TRIANGULAR MATRIX
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 COMPUTE ROWS J+JB:N OF CURRENT BLOCK COLUMN
*
                  CALL DTRMM( 'LEFT', 'LOWER', 'NO TRANSPOSE', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL DTRSM( 'RIGHT', 'LOWER', 'NO TRANSPOSE', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              COMPUTE INVERSE OF CURRENT DIAGONAL BLOCK
*
               CALL DTRTI2( 'LOWER', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     END OF DTRTRI
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     SEPTEMBER 30, 1994
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          CA, CB
*     ..
*
*  PURPOSE
*  =======
*
*  LSAME RETURNS .TRUE. IF CA IS THE SAME LETTER AS CB REGARDLESS OF
*  CASE.
*
*  ARGUMENTS
*  =========
*
*  CA      (INPUT) CHARACTER*1
*  CB      (INPUT) CHARACTER*1
*          CA AND CB SPECIFY THE SINGLE CHARACTERS TO BE COMPARED.
*
* =====================================================================
*
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ICHAR
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST IF THE CHARACTERS ARE EQUAL
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     NOW TEST FOR EQUIVALENCE IF BOTH CHARACTERS ARE ALPHABETIC.
*
      ZCODE = ICHAR( 'Z' )
*
*     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
*     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
*     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
*     ICHAR('A') ON AN EBCDIC MACHINE.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
*        UPPER CASE 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
*        UPPER CASE 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
*        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     END OF LSAME
*
      END
