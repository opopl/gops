C--------------------------------------------------------------------
C  SPARSE BLAS TOOLKIT INTERFACE ROUTINE:
C--------------------------------------------------------------------
      SUBROUTINE DCOOMM( TRANSA, M, N, K, ALPHA, DESCRA,
     *           VAL, INDX, JNDX, NNZ,
     *           B, LDB, BETA, C, LDC, WORK, LWORK)

C--------------------------------------------------------------------
C         ------------ BEGIN INTERFACE DESCRIPTION ------------
C   TOOLKIT INTERFACE:
C   DCOOMM -- COMPRESSED SPARSE ROW FORMAT MATRIX-MATRIX MULTIPLY
C  
C   C <- ALPHA A B + BETA C
C  
C   ARGUMENTS:
C  
C   INT TRANSA	INDICATES HOW TO OPERATE WITH THE SPARSE MATRIX
C  		0 : OPERATE WITH MATRIX
C  		1 : OPERATE WITH TRANSPOSE MATRIX
C  
C   INT M	NUMBER OF ROWS IN MATRIX A
C  
C   INT N	NUMBER OF COLUMNS IN MATRIX C
C  
C   INT K	NUMBER OF COLUMNS IN MATRIX A
C  
C   DOUBLE ALPHA SCALAR PARAMETER
C  
C   DOUBLE BETA  SCALAR PARAMETER
C  
C   INT DESCRA()	DESCRIPTOR ARGUMENT.  NINE ELEMENT INTEGER ARRAY
C  		DESCRA(1) MATRIX STRUCTURE
C  			0 : GENERAL
C  			1 : SYMMETRIC
C  			2 : HERMITIAN
C  			3 : TRIANGULAR
C  			4 : SKEW(ANTI)-SYMMETRIC
C  			5 : DIAGONAL
C  		DESCRA(2) UPPER/LOWER TRIANGULAR INDICATOR
C  			1 : LOWER
C  			2 : UPPER
C  		DESCRA(3) MAIN DIAGONAL TYPE
C  			0 : NON-UNIT
C  			1 : UNIT
C  		DESCRA(4) ARRAY BASE 
C  			0 : C/C++ COMPATIBLE
C  			1 : FORTRAN COMPATIBLE
C  		DESCRA(5) REPEATED INDICES?
C  			0 : UNKNOWN
C  			1 : NO REPEATED INDICES
C  
C  
C
C   DOUBLE VAL()  SCALAR ARRAY OF LENGTH NNZ CONTAINING MATRIX ENTRIES.
C  
C   INT INDX()    INTEGER ARRAY OF LENGTH NNZ CONTAINING ROW INDICES.
C
C   INT JNDX()    INTEGER ARRAY OF LENGTH NNZ CONTAINING COLUMN INDICES.
C
C   INT NNZ       NUMBER OF NON-ZERO ELEMENTS IN A.
C
C   DOUBLE B()    RECTANGULAR ARRAY WITH FIRST DIMENSION LDB.
C  
C   DOUBLE C()    RECTANGULAR ARRAY WITH FIRST DIMENSION LDC.
C  
C   DOUBLE WORK() SCRATCH ARRAY OF LENGTH LWORK.  LWORK SHOULD BE AT LEAST
C                 MAX(M,N)
C  
C       ------------ END INTERFACE DESCRIPTION --------------
C--------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER J1
C
C     INTERFACE VARIABLES:
C
      INTEGER TRANSA, M, N, K, LDB, LDC, LWORK
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION BETA
      INTEGER DESCRA(*)
      DOUBLE PRECISION B(*), C(*)
      DOUBLE PRECISION WORK(*)
C
C     FORMAT SPECIFIC INTERFACE VARIABLES:
C
      INTEGER NNZ
      INTEGER INDX(*), JNDX(*)
      DOUBLE PRECISION VAL(*)
C
C     LOCAL VARIABLES:
C
      INTEGER INFO
      CHARACTER TRANSPOSE
C
C     EXTERNALS:
C
      EXTERNAL XERBLA

C
C     TEST INPUT PARAMETERS:
C

      INFO = 0
      IF ( (TRANSA .NE. 0) .AND. (TRANSA .NE. 1) ) THEN
         INFO = 1
      ELSE IF ( M .LT. 0 ) THEN
         INFO = 2
      ELSE IF (N .LT. 0) THEN
         INFO = 3
      ELSE IF (K .LT. 0) THEN
         INFO = 4
      ELSE
        IF (TRANSA .EQ. 0) THEN
C         CHECK FOR CONSISTANT DIMENSIONS:
          IF ( LDB .LT. K ) THEN 
            INFO = 15
          ELSE IF (LDC .LT. M) THEN
            INFO = 18
          ENDIF
        ELSE IF (TRANSA .EQ. 1) THEN
C         CHECK FOR CONSISTANT DIMENSIONS:
          IF ( LDB .LT. M ) THEN 
            INFO = 15
          ELSE IF (LDC .LT. K) THEN
            INFO = 18
          ENDIF
        ENDIF
      ENDIF

      IF ( INFO .NE. 0 ) THEN
        CALL XERBLA('COOMM', INFO)
        RETURN
      ENDIF


      IF ( (DESCRA(1) .GE. 0 .AND. DESCRA(1) .LE. 5 ) .AND.
     *      ALPHA .EQ. 0.D0                                 ) THEN
C       QUICK RETURN AFTER SCALING:
        CALL DSCAL(M*N, BETA, C, 1)
        RETURN
       ENDIF
      
      TRANSPOSE = 'N'
      IF ( TRANSA .EQ. 1 ) TRANSPOSE = 'T'
 
C
C     CALL APPROPRIATE KERNEL SUBROUTINE:
C

      IF (DESCRA(1) .EQ. 0   .OR.
     *    DESCRA(1) .EQ. 3   .OR.
     *    DESCRA(1) .EQ. 5        ) THEN
C
C        GENERAL MATRIX MULTIPLY:
C
         IF (TRANSPOSE .EQ. 'N') THEN
           CALL DCOOMMGK( M, N, K, ALPHA,
     *       VAL, INDX, JNDX, NNZ,
     *       B, LDB, BETA, C, LDC, DESCRA(4))
         ELSE
           CALL DCOOMMGK( M, N, K, ALPHA,
     *       VAL, JNDX, INDX, NNZ,
     *       B, LDB, BETA, C, LDC, DESCRA(4))
         ENDIF
        RETURN
      ELSE IF (DESCRA(1) .EQ. 1  .OR. 
     *         DESCRA(1) .EQ. 2        ) THEN
C
C       SYMMETRIC/HERMITIAN  MATRIX MULTIPLY:
C
        CALL DCOOMMSK( M, N, K, ALPHA,
     *       VAL, INDX, JNDX, NNZ,
     *       B, LDB, BETA, C, LDC, DESCRA(4))
        RETURN
      ELSE IF (DESCRA(1) .EQ. 4 ) THEN
C
C        SKEW-SYMMETRIC MATRIX MULTIPLY:
C
         IF (TRANSPOSE .EQ. 'N') THEN
          CALL DCOOMMKK( M, N, K, ALPHA,
     *       VAL, INDX, JNDX, NNZ,
     *       B, LDB, BETA, C, LDC, DESCRA(4))
         ELSE
          CALL DCOOMMKK( M, N, K, ALPHA,
     *       VAL, JNDX, INDX, NNZ,
     *       B, LDB, BETA, C, LDC, DESCRA(4))
         ENDIF
        RETURN
      ELSE
        INFO = 6
      ENDIF
  
      IF ( INFO .NE. 0 ) THEN
        CALL XERBLA('COOMM', INFO)
        RETURN
      ENDIF
 
      RETURN
      END 



C--------------------------------------------------------------------
C  SPARSE BLAS KERNEL ROUTINE(S):
C--------------------------------------------------------------------
      SUBROUTINE DCOOMMGK( M, N, K, ALPHA, 
     *           VAL, INDX, JNDX, NNZ,
     *           B, LDB, BETA, C, LDC, BASE)
      IMPLICIT NONE
      INTEGER M, N, K, NNZ, LDB, LDC, BASE
      INTEGER L,J
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION BETA
      INTEGER INDX(*), JNDX(*)
      DOUBLE PRECISION VAL(*)
      DOUBLE PRECISION B(LDB,*), C(LDC,*)
C
C     SCALE C BY BETA:
C
C DAE 
C
      CALL DSCAL( N*LDC, BETA, C(1,1), 1)
  
      DO 5 L=1,N
         DO 10 J=1,NNZ
           C(INDX(J),L) = C(INDX(J),L)+ALPHA*B(JNDX(J),L)*VAL(J)
 10      CONTINUE
 5    CONTINUE

      RETURN
      END
         
      SUBROUTINE DCOOMMSK( M, N, K, ALPHA, VAL, INDX, JNDX, NNZ,
     *                          B, LDB, BETA, C, LDC, BASE)
      IMPLICIT NONE
      INTEGER M, N, K, NNZ, LDB, LDC, BASE
      INTEGER L,J
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION BETA
      INTEGER INDX(*), JNDX(*)
      DOUBLE PRECISION VAL(*)
      DOUBLE PRECISION B(LDB,*), C(LDC,*)
C
C     SCALE C BY BETA:
C
      CALL DSCAL( N*M, BETA, C(1,1), 1)

      DO 5 L=1,N
         DO 10 J=1,NNZ
           IF ( INDX(J) .NE. JNDX(J) ) THEN
              C(INDX(J),L) = C(INDX(J),L)+ALPHA*B(JNDX(J),L)*VAL(J)
              C(JNDX(J),L) = C(JNDX(J),L)+ALPHA*B(INDX(J),L)*VAL(J)
           ELSE
              C(INDX(J),L) = C(INDX(J),L)+ALPHA*B(JNDX(J),L)*VAL(J)
           ENDIF
 10      CONTINUE
 5    CONTINUE

      RETURN
      END
         
      SUBROUTINE DCOOMMKK( M, N, K, ALPHA, VAL, INDX, JNDX, NNZ,
     *                          B, LDB, BETA, C, LDC, BASE)
      IMPLICIT NONE
      INTEGER M, N, K, NNZ, LDB, LDC, BASE
      INTEGER L,J
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION BETA
      INTEGER INDX(*), JNDX(*)
      DOUBLE PRECISION VAL(*)
      DOUBLE PRECISION B(LDB,*), C(LDC,*)
C
C     SCALE C BY BETA:
C
      CALL DSCAL( N*M, BETA, C(1,1), 1)

      DO 5 L=1,N
         DO 10 J=1,NNZ
           IF ( INDX(J) .NE. JNDX(J) ) THEN
             C(INDX(J),L) = C(INDX(J),L)+ALPHA*B(JNDX(J),L)*VAL(J)
             C(JNDX(J),L) = C(JNDX(J),L)-ALPHA*B(INDX(J),L)*VAL(J)
           ENDIF
 10      CONTINUE
 5    CONTINUE

      RETURN
      END
         
