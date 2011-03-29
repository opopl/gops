      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGETF2 COMPUTES AN LU FACTORIZATION OF A GENERAL M-BY-N MATRIX A
*  USING PARTIAL PIVOTING WITH ROW INTERCHANGES.
*
*  THE FACTORIZATION HAS THE FORM
*     A = P * L * U
*  WHERE P IS A PERMUTATION MATRIX, L IS LOWER TRIANGULAR WITH UNIT
*  DIAGONAL ELEMENTS (LOWER TRAPEZOIDAL IF M > N), AND U IS UPPER
*  TRIANGULAR (UPPER TRAPEZOIDAL IF M < N).
*
*  THIS IS THE RIGHT-LOOKING LEVEL 2 BLAS VERSION OF THE ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX A.  M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE M BY N MATRIX TO BE FACTORED.
*          ON EXIT, THE FACTORS L AND U FROM THE FACTORIZATION
*          A = P*L*U; THE UNIT DIAGONAL ELEMENTS OF L ARE NOT STORED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (MIN(M,N))
*          THE PIVOT INDICES; FOR 1 <= I <= MIN(M,N), ROW I OF THE
*          MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0: IF INFO = K, U(K,K) IS EXACTLY ZERO. THE FACTORIZATION
*               HAS BEEN COMPLETED, BUT THE FACTOR U IS EXACTLY
*               SINGULAR, AND DIVISION BY ZERO WILL OCCUR IF IT IS USED
*               TO SOLVE A SYSTEM OF EQUATIONS.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            J, JP
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        FIND PIVOT AND TEST FOR SINGULARITY.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           APPLY THE INTERCHANGE TO COLUMNS 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           COMPUTE ELEMENTS J+1:M OF J-TH COLUMN.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           UPDATE TRAILING SUBMATRIX.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     END OF DGETF2
*
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGETRF COMPUTES AN LU FACTORIZATION OF A GENERAL M-BY-N MATRIX A
*  USING PARTIAL PIVOTING WITH ROW INTERCHANGES.
*
*  THE FACTORIZATION HAS THE FORM
*     A = P * L * U
*  WHERE P IS A PERMUTATION MATRIX, L IS LOWER TRIANGULAR WITH UNIT
*  DIAGONAL ELEMENTS (LOWER TRAPEZOIDAL IF M > N), AND U IS UPPER
*  TRIANGULAR (UPPER TRAPEZOIDAL IF M < N).
*
*  THIS IS THE RIGHT-LOOKING LEVEL 3 BLAS VERSION OF THE ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX A.  M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE M-BY-N MATRIX TO BE FACTORED.
*          ON EXIT, THE FACTORS L AND U FROM THE FACTORIZATION
*          A = P*L*U; THE UNIT DIAGONAL ELEMENTS OF L ARE NOT STORED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (MIN(M,N))
*          THE PIVOT INDICES; FOR 1 <= I <= MIN(M,N), ROW I OF THE
*          MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0:  IF INFO = I, U(I,I) IS EXACTLY ZERO. THE FACTORIZATION
*                HAS BEEN COMPLETED, BUT THE FACTOR U IS EXACTLY
*                SINGULAR, AND DIVISION BY ZERO WILL OCCUR IF IT IS USED
*                TO SOLVE A SYSTEM OF EQUATIONS.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     DETERMINE THE BLOCK SIZE FOR THIS ENVIRONMENT.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        USE UNBLOCKED CODE.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        USE BLOCKED CODE.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           FACTOR DIAGONAL AND SUBDIAGONAL BLOCKS AND TEST FOR EXACT
*           SINGULARITY.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           ADJUST INFO AND THE PIVOT INDICES.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           APPLY INTERCHANGES TO COLUMNS 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              APPLY INTERCHANGES TO COLUMNS J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              COMPUTE BLOCK ROW OF U.
*
               CALL DTRSM( 'LEFT', 'LOWER', 'NO TRANSPOSE', 'UNIT', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 UPDATE TRAILING SUBMATRIX.
*
                  CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     END OF DGETRF
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1999
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLASWP PERFORMS A SERIES OF ROW INTERCHANGES ON THE MATRIX A.
*  ONE ROW INTERCHANGE IS INITIATED FOR EACH OF ROWS K1 THROUGH K2 OF A.
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE MATRIX OF COLUMN DIMENSION N TO WHICH THE ROW
*          INTERCHANGES WILL BE APPLIED.
*          ON EXIT, THE PERMUTED MATRIX.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.
*
*  K1      (INPUT) INTEGER
*          THE FIRST ELEMENT OF IPIV FOR WHICH A ROW INTERCHANGE WILL
*          BE DONE.
*
*  K2      (INPUT) INTEGER
*          THE LAST ELEMENT OF IPIV FOR WHICH A ROW INTERCHANGE WILL
*          BE DONE.
*
*  IPIV    (INPUT) INTEGER ARRAY, DIMENSION (M*ABS(INCX))
*          THE VECTOR OF PIVOT INDICES.  ONLY THE ELEMENTS IN POSITIONS
*          K1 THROUGH K2 OF IPIV ARE ACCESSED.
*          IPIV(K) = L IMPLIES ROWS K AND L ARE TO BE INTERCHANGED.
*
*  INCX    (INPUT) INTEGER
*          THE INCREMENT BETWEEN SUCCESSIVE VALUES OF IPIV.  IF IPIV
*          IS NEGATIVE, THE PIVOTS ARE APPLIED IN REVERSE ORDER.
*
*  FURTHER DETAILS
*  ===============
*
*  MODIFIED BY
*   R. C. WHALEY, COMPUTER SCIENCE DEPT., UNIV. OF TENN., KNOXVILLE, USA
*
* =====================================================================
*
*     .. LOCAL SCALARS ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     INTERCHANGE ROW I WITH ROW IPIV(I) FOR EACH OF ROWS K1 THROUGH K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     END OF DLASWP
*
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1998
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  PURPOSE
*  =======
*
*  IEEECK IS CALLED FROM THE ILAENV TO VERIFY THAT INFINITY AND
*  POSSIBLY NAN ARITHMETIC IS SAFE (I.E. WILL NOT TRAP).
*
*  ARGUMENTS
*  =========
*
*  ISPEC   (INPUT) INTEGER
*          SPECIFIES WHETHER TO TEST JUST FOR INIFINITY ARITHMETIC
*          OR WHETHER TO TEST FOR INFINITY AND NAN ARITHMETIC.
*          = 0: VERIFY INFINITY ARITHMETIC ONLY.
*          = 1: VERIFY INFINITY AND NAN ARITHMETIC.
*
*  ZERO    (INPUT) REAL
*          MUST CONTAIN THE VALUE 0.0
*          THIS IS PASSED TO PREVENT THE COMPILER FROM OPTIMIZING
*          AWAY THIS CODE.
*
*  ONE     (INPUT) REAL
*          MUST CONTAIN THE VALUE 1.0
*          THIS IS PASSED TO PREVENT THE COMPILER FROM OPTIMIZING
*          AWAY THIS CODE.
*
*  RETURN VALUE:  INTEGER
*          = 0:  ARITHMETIC FAILED TO PRODUCE THE CORRECT ANSWERS
*          = 1:  ARITHMETIC PRODUCED THE CORRECT ANSWERS
*
*     .. LOCAL SCALARS ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. EXECUTABLE STATEMENTS ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     RETURN IF WE WERE ONLY ASKED TO CHECK INFINITY ARITHMETIC
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1999
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  PURPOSE
*  =======
*
*  ILAENV IS CALLED FROM THE LAPACK ROUTINES TO CHOOSE PROBLEM-DEPENDENT
*  PARAMETERS FOR THE LOCAL ENVIRONMENT.  SEE ISPEC FOR A DESCRIPTION OF
*  THE PARAMETERS.
*
*  THIS VERSION PROVIDES A SET OF PARAMETERS WHICH SHOULD GIVE GOOD,
*  BUT NOT OPTIMAL, PERFORMANCE ON MANY OF THE CURRENTLY AVAILABLE
*  COMPUTERS.  USERS ARE ENCOURAGED TO MODIFY THIS SUBROUTINE TO SET
*  THE TUNING PARAMETERS FOR THEIR PARTICULAR MACHINE USING THE OPTION
*  AND PROBLEM SIZE INFORMATION IN THE ARGUMENTS.
*
*  THIS ROUTINE WILL NOT FUNCTION CORRECTLY IF IT IS CONVERTED TO ALL
*  LOWER CASE.  CONVERTING IT TO ALL UPPER CASE IS ALLOWED.
*
*  ARGUMENTS
*  =========
*
*  ISPEC   (INPUT) INTEGER
*          SPECIFIES THE PARAMETER TO BE RETURNED AS THE VALUE OF
*          ILAENV.
*          = 1: THE OPTIMAL BLOCKSIZE; IF THIS VALUE IS 1, AN UNBLOCKED
*               ALGORITHM WILL GIVE THE BEST PERFORMANCE.
*          = 2: THE MINIMUM BLOCK SIZE FOR WHICH THE BLOCK ROUTINE
*               SHOULD BE USED; IF THE USABLE BLOCK SIZE IS LESS THAN
*               THIS VALUE, AN UNBLOCKED ROUTINE SHOULD BE USED.
*          = 3: THE CROSSOVER POINT (IN A BLOCK ROUTINE, FOR N LESS
*               THAN THIS VALUE, AN UNBLOCKED ROUTINE SHOULD BE USED)
*          = 4: THE NUMBER OF SHIFTS, USED IN THE NONSYMMETRIC
*               EIGENVALUE ROUTINES
*          = 5: THE MINIMUM COLUMN DIMENSION FOR BLOCKING TO BE USED;
*               RECTANGULAR BLOCKS MUST HAVE DIMENSION AT LEAST K BY M,
*               WHERE K IS GIVEN BY ILAENV(2,...) AND M BY ILAENV(5,...)
*          = 6: THE CROSSOVER POINT FOR THE SVD (WHEN REDUCING AN M BY N
*               MATRIX TO BIDIAGONAL FORM, IF MAX(M,N)/MIN(M,N) EXCEEDS
*               THIS VALUE, A QR FACTORIZATION IS USED FIRST TO REDUCE
*               THE MATRIX TO A TRIANGULAR FORM.)
*          = 7: THE NUMBER OF PROCESSORS
*          = 8: THE CROSSOVER POINT FOR THE MULTISHIFT QR AND QZ METHODS
*               FOR NONSYMMETRIC EIGENVALUE PROBLEMS.
*          = 9: MAXIMUM SIZE OF THE SUBPROBLEMS AT THE BOTTOM OF THE
*               COMPUTATION TREE IN THE DIVIDE-AND-CONQUER ALGORITHM
*               (USED BY XGELSD AND XGESDD)
*          =10: IEEE NAN ARITHMETIC CAN BE TRUSTED NOT TO TRAP
*          =11: INFINITY ARITHMETIC CAN BE TRUSTED NOT TO TRAP
*
*  NAME    (INPUT) CHARACTER*(*)
*          THE NAME OF THE CALLING SUBROUTINE, IN EITHER UPPER CASE OR
*          LOWER CASE.
*
*  OPTS    (INPUT) CHARACTER*(*)
*          THE CHARACTER OPTIONS TO THE SUBROUTINE NAME, CONCATENATED
*          INTO A SINGLE CHARACTER STRING.  FOR EXAMPLE, UPLO = 'U',
*          TRANS = 'T', AND DIAG = 'N' FOR A TRIANGULAR ROUTINE WOULD
*          BE SPECIFIED AS OPTS = 'UTN'.
*
*  N1      (INPUT) INTEGER
*  N2      (INPUT) INTEGER
*  N3      (INPUT) INTEGER
*  N4      (INPUT) INTEGER
*          PROBLEM DIMENSIONS FOR THE SUBROUTINE NAME; THESE MAY NOT ALL
*          BE REQUIRED.
*
* (ILAENV) (OUTPUT) INTEGER
*          >= 0: THE VALUE OF THE PARAMETER SPECIFIED BY ISPEC
*          < 0:  IF ILAENV = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE.
*
*  FURTHER DETAILS
*  ===============
*
*  THE FOLLOWING CONVENTIONS HAVE BEEN USED WHEN CALLING ILAENV FROM THE
*  LAPACK ROUTINES:
*  1)  OPTS IS A CONCATENATION OF ALL OF THE CHARACTER OPTIONS TO
*      SUBROUTINE NAME, IN THE SAME ORDER THAT THEY APPEAR IN THE
*      ARGUMENT LIST FOR NAME, EVEN IF THEY ARE NOT USED IN DETERMINING
*      THE VALUE OF THE PARAMETER SPECIFIED BY ISPEC.
*  2)  THE PROBLEM DIMENSIONS N1, N2, N3, N4 ARE SPECIFIED IN THE ORDER
*      THAT THEY APPEAR IN THE ARGUMENT LIST FOR NAME.  N1 IS USED
*      FIRST, N2 SECOND, AND SO ON, AND UNUSED PROBLEM DIMENSIONS ARE
*      PASSED A VALUE OF -1.
*  3)  THE PARAMETER VALUE RETURNED BY ILAENV IS CHECKED FOR VALIDITY IN
*      THE CALLING SUBROUTINE.  FOR EXAMPLE, ILAENV IS USED TO RETRIEVE
*      THE OPTIMAL BLOCKSIZE FOR STRTRI AS FOLLOWS:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. LOCAL SCALARS ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,
     $        1100 ) ISPEC
*
*     INVALID VALUE FOR ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     CONVERT NAME TO UPPER CASE IF THE FIRST CHARACTER IS LOWER CASE.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII CHARACTER SET
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC CHARACTER SET
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        PRIME MACHINES:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  BLOCK SIZE
*
*     IN THESE EXAMPLES, SEPARATE CODE IS PROVIDED FOR SETTING NB FOR
*     REAL AND COMPLEX.  WE ASSUME THAT NB WILL TAKE THE SAME VALUE IN
*     SINGLE OR DOUBLE PRECISION.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  MINIMUM BLOCK SIZE
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  CROSSOVER POINT
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  NUMBER OF SHIFTS (USED BY XHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  MINIMUM COLUMN DIMENSION (NOT USED)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  CROSSOVER POINT FOR SVD (USED BY XGELSS AND XGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  NUMBER OF PROCESSORS (NOT USED)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  CROSSOVER POINT FOR MULTISHIFT (USED BY XHSEQR)
*
      ILAENV = 50
      RETURN
*
  900 CONTINUE
*
*     ISPEC = 9:  MAXIMUM SIZE OF THE SUBPROBLEMS AT THE BOTTOM OF THE
*                 COMPUTATION TREE IN THE DIVIDE-AND-CONQUER ALGORITHM
*                 (USED BY XGELSD AND XGESDD)
*
      ILAENV = 25
      RETURN
*
 1000 CONTINUE
*
*     ISPEC = 10: IEEE NAN ARITHMETIC CAN BE TRUSTED NOT TO TRAP
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 ) 
      END IF
      RETURN
*
 1100 CONTINUE
*
*     ISPEC = 11: INFINITY ARITHMETIC CAN BE TRUSTED NOT TO TRAP
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 ) 
      END IF
      RETURN
*
*     END OF ILAENV
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 3.0) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     SEPTEMBER 30, 1994
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  PURPOSE
*  =======
*
*  XERBLA  IS AN ERROR HANDLER FOR THE LAPACK ROUTINES.
*  IT IS CALLED BY AN LAPACK ROUTINE IF AN INPUT PARAMETER HAS AN
*  INVALID VALUE.  A MESSAGE IS PRINTED AND EXECUTION STOPS.
*
*  INSTALLERS MAY CONSIDER MODIFYING THE STOP STATEMENT IN ORDER TO
*  CALL SYSTEM-SPECIFIC EXCEPTION-HANDLING FACILITIES.
*
*  ARGUMENTS
*  =========
*
*  SRNAME  (INPUT) CHARACTER*6
*          THE NAME OF THE ROUTINE WHICH CALLED XERBLA.
*
*  INFO    (INPUT) INTEGER
*          THE POSITION OF THE INVALID PARAMETER IN THE PARAMETER LIST
*          OF THE CALLING ROUTINE.
*
* =====================================================================
*
*     .. EXECUTABLE STATEMENTS ..
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** ON ENTRY TO ', A6, ' PARAMETER NUMBER ', I2, ' HAD ',
     $      'AN ILLEGAL VALUE' )
*
*     END OF XERBLA
*
      END
