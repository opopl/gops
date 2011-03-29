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
C     INTERFACE TO SPJV.F FOR CALCULATING MINIMUM DISTANCE
C     OF TWO ATOMIC CONFIGURATIONS WITH RESPECT TO
C     PARTICLE PERMUTATIONS.
C     THE FUNCTION PERMDIST DETERMINES THE DISTANCE OR WEIGHT FUNCTION,
C     MINPERM IS THE MAIN ROUTINE.
C
C         TOMAS OPPELSTRUP, JUL 10, 2003
C         TOMASO@NADA.KTH.SE
C

C     THIS IS THE MAIN ROUTINE FOR MINIMUM DISTANCE CALCULATION.
C     GIVEN TWO COORDINATE VECTORS P,Q OF PARTICLES EACH, RETURN
C     THE MINIMUM DISTANCE IN DIST, AND THE PERMUTATION IN PERM.
C     PERM IS AN INTEGER VECTOR SUCH THAT
C       P(I) <--> Q(PERM(I))
C     I.E.
C       SUM(I=1,N) PERMDIST(P(I), Q(PERM(I))) == DIST
C
      SUBROUTINE MINPERM(N, P, Q, SX, SY, SZ, PBC, PERM, DIST, WORSTDIST, WORSTRADIUS)
      IMPLICIT NONE

C     INPUT
C       N  : SYSTEM SIZE
C       P,Q: COORDINATE VECTORS (N PARTICLES)
C       S  : BOX LENGTHS (OR DUMMY IF OPEN B.C.)
C       PBC: PERIODIC BOUNDARY CONDITIONS?
      INTEGER N
      DOUBLE PRECISION P(3*N), Q(3*N), S(3), SX, SY, SZ, WORSTDIST, WORSTRADIUS
      LOGICAL PBC

C     OUTPUT
C       PERM: PERMUTATION SO THAT P(I) <--> Q(PERM(I))
C       DIST: MINIMUM ATTAINABLE DISTANCE
C     WE HAVE
      INTEGER PERM(N)
      DOUBLE PRECISION DIST, DUMMY
      
C     PARAMETERS
C       SCALE : PRECISION
C       MAXNEI: MAXIMUM NUMBER OF CLOSEST NEIGHBOURS
      DOUBLE PRECISION SCALE
      INTEGER MAXNEI

      PARAMETER (SCALE = 1.0D6   )
      PARAMETER (MAXNEI = 60     )

C     INTERNAL VARIABLES
C     CC, KK, FIRST:
C       SPARSE MATRIX OF DISTANCES
C     FIRST(I):
C       BEGINNING OF ROW I IN DATA,INDEX VECTORS
C     KK(FIRST(I)..FIRST(I+1)-1):
C       COLUMN INDEXES OF EXISTING ELEMENTS IN ROW I
C     CC(FIRST(I)..FIRST(I+1)-1):
C       MATRIX ELEMENTS OF ROW I
      INTEGER*8 KK(N*MAXNEI), FIRST(N+1), X(N), Y(N)
      INTEGER*8 CC(N*MAXNEI), U(N), V(N), H
      INTEGER   M, I, J, K, L, L2, T, A, I3, J3
      INTEGER*8 N8, SZ8, D
      INTEGER NDONE, J1, J2

C     DISTANCE FUNCTION
      DOUBLE PRECISION PERMDIST

      S(1)=SX
      S(2)=SY
      S(3)=SZ
      M = MAXNEI
      IF(N .LE. MAXNEI) M = N
      SZ8 = M*N
      N8 = N

      DO I=0,N
         FIRST(I+1) = I*M + 1
      ENDDO
   
      IF(M .EQ. N) THEN
C     COMPUTE THE FULL MATRIX...

         DO I=1,N
            K = FIRST(I)-1
            DO J=1,N
               CC(K+J) = PERMDIST(P(3*I-2), Q(3*J-2), S, PBC)*SCALE
               KK(K+J) = J
C              WRITE(*,*) I, J, '-->', CC(K+J)
            ENDDO
         ENDDO
      ELSE
C     WE NEED TO STORE THE DISTANCES OF THE MAXNEI CLOSEEST NEIGHBORS
C     OF EACH PARTICLE. THE FOLLOWING BUILDS A HEAP TO KEEP TRACK OF
C     THE MAXNEI CLOSEST NEIGHBOURS SEEN SO FAR. IT MIGHT BE MORE
C     EFFICIENT TO USE QUICK-SELECT INSTEAD... (THIS IS DEFINATELY
C     TRUE IN THE LIMIT OF INFINITE SYSTEMS.)
        DO I=1,N
           K = FIRST(I)-1
           DO J=1,M
              CC(K+J) = PERMDIST(P(3*I-2), Q(3*J-2), S, PBC)*SCALE
              KK(K+J) = J
              L = J
10            IF(L .LE. 1) GOTO 11
              L2 = L/2
              IF(CC(K+L2) .LT. CC(K+L)) THEN
                 H = CC(K+L2)
                 CC(K+L2) = CC(K+L)
                 CC(K+L) = H
                 T = KK(K+L2)
                 KK(K+L2) = KK(K+L)
                 KK(K+L) = T
                 L = L2
                 GOTO 10
              ENDIF
11         ENDDO
           
           DO J=M+1,N
              D = PERMDIST(P(3*I-2), Q(3*J-2), S, PBC)*SCALE
              IF(D .LT. CC(K+1)) THEN
                 CC(K+1) = D
                 KK(K+1) = J
                 L = 1
20               L2 = 2*L
                 IF(L2+1 .GT. M) GOTO 21
                 IF(CC(K+L2+1) .GT. CC(K+L2)) THEN
                    A = K+L2+1
                 ELSE
                    A = K+L2
                 ENDIF
                 IF(CC(A) .GT. CC(K+L)) THEN
                    H = CC(A)
                    CC(A) = CC(K+L)
                    CC(K+L) = H
                    T = KK(A)
                    KK(A) = KK(K+L)
                    KK(K+L) = T
                    L = A-K
                    GOTO 20
                 ENDIF
21               IF (L2 .LE. M) THEN ! SPLIT IF STATEMENTS TO AVOID A SEGMENTATION FAULT
                    IF (CC(K+L2) .GT. CC(K+L)) THEN
                       H = CC(K+L2)
                       CC(K+L2) = CC(K+L)
                       CC(K+L) = H
                       T = KK(K+L2)
                       KK(K+L2) = KK(K+L)
                       KK(K+L) = T
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
!       PRINT '(A,I6,A)','ATOM ',I,' NEAREST NEIGHBOURS AND DISTANCES:'
!       PRINT '(20I6)',KK(M*(I-1)+1:M*I)
!       PRINT '(12I15)',CC(M*(I-1)+1:M*I)
           
        ENDDO    
!
! CREATE AND MAINTAIN AN ORDERED LIST, SMALLEST TO LARGEST FROM KK(M*(I-1)+1:M*I) FOR ATOM I.
! NOTE THAT THERE IS NO SYMMETRY WITH RESPECT TO EXCHANGE OF I AND J!
! THIS RUNS SLOWER THAN THE ABOVE HEAP ALGORITHM.
!
!          CC(1:N*M)=HUGE(1)
!           DO I=1,N
!              K=FIRST(I)-1
!              DO J=1,N
!                 D=PERMDIST(P(3*I-2), Q(3*J-2), S, PBC)*SCALE
!                 IF (D.GT.CC(K+M)) CYCLE
!                 DO J1=M-1,1,-1
!                    IF (D.LT.CC(K+J1)) THEN
!                       CC(K+J1+1)=CC(K+J1)
!                       KK(K+J1+1)=KK(K+J1)
!                    ELSE
!                       CC(K+J1+1)=D
!                       KK(K+J1+1)=J
!                       GOTO 112
!                    ENDIF
!                 ENDDO
! !
! !  IF WE REACH THIS POINT THEN WE NEED TO INSERT AT THE BEGINNING.
! !
!                 CC(K+1)=D
!                 KK(K+1)=J
! 112             CONTINUE
!              ENDDO
!           ENDDO
      ENDIF

C     CALL BIPARTITE MATCHING ROUTINE
      CALL JOVOSAP(N8, SZ8, CC, KK, FIRST, X, Y, U, V, H)

      IF(H .LT. 0) THEN
C     IF INITIAL GUESS CORRECT, DEDUCE SOLUTION DISTANCE
C     WHICH IS NOT DONE IN JOVOSAP
         H = 0
         DO I=1,N
            J = FIRST(I)
 30         IF (J.GT.N*MAXNEI) THEN
!              PRINT '(A,I6,A)','MINPERM> WARNING A - MATCHING FAILED'
               DO J1=1,N
                  PERM(J1)=J1
               ENDDO
               RETURN
            ENDIF
            IF(KK(J) .NE. X(I)) THEN
               J = J + 1
               GOTO 30
            ENDIF
            H = H + CC(J)
         ENDDO
      ENDIF

      DO I=1,N
         PERM(I) = X(I)
         IF (PERM(I).GT.N) PERM(I)=N
         IF (PERM(I).LT.1) PERM(I)=1
      ENDDO

      DIST = DBLE(H) / SCALE

      WORSTDIST=-1.0D0
      DO I=1,N
        DUMMY=(P(3*(I-1)+1)-Q(3*(PERM(I)-1)+1))**2+(P(3*(I-1)+2)-Q(3*(PERM(I)-1)+2))**2+(P(3*(I-1)+3)-Q(3*(PERM(I)-1)+3))**2
         IF (DUMMY.GT.WORSTDIST) THEN
            WORSTDIST=DUMMY 
            WORSTRADIUS=P(3*(I-1)+1)**2+P(3*(I-1)+2)**2+P(3*(I-1)+3)**2
         ENDIF
      ENDDO
      WORSTDIST=SQRT(WORSTDIST)
      WORSTRADIUS=MAX(SQRT(WORSTRADIUS),1.0D0)

      END
      
C     PERMDIST IS THE DISTANCE OR WEIGHT FUNCTION. IT IS CODED
C     SEPARATELY FOR CLARITY. JUST HOPE THAT THE COMPILER
C     KNOWS HOW TO TO DO PROPER INLINING!
C     INPUT
C       P,Q: COORDINATES
C       S  : BOXLENGTHS (OR DUMMY IF OPEN B.C.)
C       PBC: PERIODIC BOUNDARY CONDITIONS?

      DOUBLE PRECISION FUNCTION PERMDIST(P, Q, S, PBC)
      IMPLICIT NONE
      DOUBLE PRECISION P(3), Q(3), S(3)
      LOGICAL PBC

      DOUBLE PRECISION T, D
      INTEGER I

      D = 0.0D0
      IF (PBC) THEN
         DO I=1,3
            IF (S(I).NE.0.0D0) THEN
               T = Q(I) - P(I)
               T = T - S(I)*ANINT(T/S(I))
               D = D + T*T
            ENDIF
         ENDDO
      ELSE
         D= (Q(1) - P(1))**2+(Q(2) - P(2))**2+(Q(3) - P(3))**2
      ENDIF
      PERMDIST = D
      END

C     THE FOLLOWING ROUTINE PERFORMS WEIGHTED BIPARTITE MATCHING FOR
C     FOR A SPARSE NON-NEGATIVE INTEGER WEIGHT MATRIX.
C     THE ORIGINAL SOURCE IS
C         HTTP://WWW.MAGICLOGIC.COM/ASSIGNMENT.HTML
C     A PUBLICATION REFERENCE CAN BE FOUND ON THE ABOVE HOMEPAGE AND
C     IN A COMMENT BELOW
C     

      SUBROUTINE JOVOSAP(N,SZ,CC,KK,FIRST,X,Y,U,V,H)
      IMPLICIT NONE
      INTEGER*8 N, SZ
      INTEGER*8 CC(SZ),KK(SZ),FIRST(N+1),X(N),Y(N),U(N),V(N)
      INTEGER*8 H,CNT,L0,T,T0,TD,V0,VJ,DJ
      INTEGER*8 LAB(N),D(N),FREE(N),TODO(N)
      LOGICAL OK(N)
      INTEGER*8 J, I, J0, L, J1, MIN, K, I0
      INTEGER*8 BIGINT

C     I DON'T KNOW HOW TO MAKE G77 READ INTEGER*8 CONSTANTS/PARAMETERS.
C       PARAMETER (BIGINT = 10**12) DOES NOT WORK(!)
C     NOR DOES
C       PARAMETER (BIGINT = 1000000000000)
C     BUT THIS SEEMS TO BE OK:
      BIGINT = 10**9
      BIGINT = BIGINT * 1000

C
C THIS SUBROUTINE SOLVES THE SPARSE LINEAR ASSIGNMENT PROBLEM
C ACCORDING 
C
C   "A SHORTEST AUGMENTING PATH ALGORITHM FOR DENSE AND SPARSE LINEAR   
C    ASSIGNMENT PROBLEMS," COMPUTING 38, 325-340, 1987
C   
C   BY
C   
C   R. JONKER AND A. VOLGENANT, UNIVERSITY OF AMSTERDAM.
C
C
C INPUT PARAMETERS :
C N = NUMBER OF ROWS AND COLUMNS
C C = WEIGHT MATRIX
C
C OUTPUT PARAMETERS
C X = COL ASSIGNED TO ROW
C Y = ROW ASSIGNED TO COL
C U = DUAL ROW VARIABLE
C V = DUAL COLUMN VARIABLE
C H = VALUE OF OPTIMAL SOLUTION
C
C INITIALIZATION

C     NEXT LINE ADDED BY TOMASO@NADA.KTH.SE, TO ENABLE DETECTION
C     OF SOLUTIONS BEING EQUIVALENT TO THE INITIAL GUESS

C
C  IF Y(:) IS INITIALISED TO ZERO THEN WE SEE SEGMENTATION FAULTS IF 
C  A Y ELEMENT IS UNSET, ETC.
C

      Y(1:N) = 0
      X(1:N) = 0
      TODO(1:N)=0
      H = -1
      DO 10 J=1,N
         V(J)=BIGINT
   10 CONTINUE
      DO 20 I=1,N
         X(I)=0
         DO 15 T=FIRST(I),FIRST(I+1)-1
            J=KK(T)
            IF (CC(T).LT.V(J)) THEN
              V(J)=CC(T)
              Y(J)=I
            END IF
   15    CONTINUE
   20 CONTINUE
      DO 30 J=1,N
         J0=N-J+1
         I=Y(J0)
         IF (I.EQ.0) THEN
!           PRINT '(A,I6,A)','MINPERM> WARNING B - MATCHING FAILED'
            RETURN
         ENDIF
         IF (X(I).NE.0) THEN
           X(I)=-ABS(X(I))
           Y(J0)=0
         ELSE
           X(I)=J0
         END IF
   30 CONTINUE
      L=0
      DO 40 I=1,N
         IF (X(I).EQ.0) THEN
           L=L+1
           FREE(L)=I
           GOTO 40
         END IF
         IF (X(I).LT.0) THEN
           X(I)=-X(I)
         ELSE
           J1=X(I)
           MIN=BIGINT
           DO 31 T=FIRST(I),FIRST(I+1)-1
              J=KK(T)
              IF (J.EQ.J1) GOTO 31
              IF (CC(T)-V(J).LT.MIN) MIN=CC(T)-V(J)
   31      CONTINUE
           V(J1)=V(J1)-MIN
         END IF
   40 CONTINUE
C IMPROVE THE INITIAL SOLUTION
      CNT=0
      IF (L.EQ.0) GOTO 1000
   41 L0=L
      K=1
      L=0
   50 I=FREE(K)
      K=K+1
      V0=BIGINT
      VJ=BIGINT
      DO 42 T=FIRST(I),FIRST(I+1)-1
         J=KK(T)
         H=CC(T)-V(J)
         IF (H.LT.VJ) THEN
           IF (H.GE.V0) THEN
             VJ=H
             J1=J
           ELSE
             VJ=V0
             V0=H
             J1=J0
             J0=J
           END IF
         END IF
   42 CONTINUE
      I0=Y(J0)
      IF (V0.LT.VJ) THEN
        V(J0)=V(J0)-VJ+V0
      ELSE
        IF (I0.EQ.0) GOTO 43
        J0=J1
        I0=Y(J1)
      END IF
      IF (I0.EQ.0) GOTO 43
      IF (V0.LT.VJ) THEN
        K=K-1
        FREE(K)=I0
      ELSE
        L=L+1
        FREE(L)=I0
      END IF
   43 X(I)=J0
      Y(J0)=I
      IF (K.LE.L0) GOTO 50
      CNT=CNT+1
      IF ((L.GT.0).AND.(CNT.LT.2)) GOTO 41
C AUGMENTATION PART
      L0=L
      DO 90 L=1,L0
         DO 51 J=1,N
            OK(J)=.FALSE.
            D(J)=BIGINT
   51    CONTINUE
         MIN=BIGINT
         I0=FREE(L)
         TD=N
         DO 52 T=FIRST(I0),FIRST(I0+1)-1
            J=KK(T)
            DJ=CC(T)-V(J)
            D(J)=DJ
            LAB(J)=I0
            IF (DJ.LE.MIN) THEN
              IF (DJ.LT.MIN) THEN
                MIN=DJ
                K=1
                TODO(1)=J
              ELSE
                K=K+1
                TODO(K)=J
              END IF
            END IF
   52    CONTINUE
         DO 53 H=1,K
            J=TODO(H)
            IF (J.EQ.0) THEN
!              PRINT '(A,I6,A)','MINPERM> WARNING C - MATCHING FAILED'
               RETURN
            ENDIF
            IF (Y(J).EQ.0) GOTO 80
            OK(J)=.TRUE.
   53    CONTINUE
C REPEAT UNTIL A FREE ROW HAS BEEN FOUND
   60    IF (K.EQ.0) THEN
!           PRINT '(A,I6,A)','MINPERM> WARNING D - MATCHING FAILED'
            RETURN
         ENDIF
         J0=TODO(K)
         K=K-1
         I=Y(J0)
         TODO(TD)=J0
         TD=TD-1
         T0=FIRST(I)
         T=T0-1
   61    T=T+1
         IF (KK(T).NE.J0) GOTO 61
         H=CC(T)-V(J0)-MIN
         DO 62 T=T0,FIRST(I+1)-1
            J=KK(T)
            IF (.NOT. OK(J)) THEN
              VJ=CC(T)-H-V(J)
              IF (VJ.LT.D(J)) THEN
                D(J)=VJ
                LAB(J)=I
                IF (VJ.EQ.MIN) THEN
                  IF (Y(J).EQ.0) GOTO 70
                  K=K+1
                  TODO(K)=J
                  OK(J)=.TRUE.
                END IF
              END IF
            END IF
   62    CONTINUE
         IF (K.NE.0) GOTO 60
         MIN=BIGINT-1
         DO 63 J=1,N
            IF (D(J).LE.MIN) THEN
              IF (.NOT. OK(J)) THEN
                IF (D(J).LT.MIN) THEN
                  MIN=D(J)
                  K=1
                  TODO(1)=J
                ELSE
                  K=K+1
                  TODO(K)=J
                END IF
              END IF
            END IF
   63    CONTINUE
         DO 64 J0=1,K
            J=TODO(J0)
            IF (Y(J).EQ.0) GOTO 70
            OK(J)=.TRUE.
   64    CONTINUE
         GOTO 60
   70    IF (MIN.EQ.0) GOTO 80
         DO 71 K=TD+1,N
            J0=TODO(K)
            V(J0)=V(J0)+D(J0)-MIN
   71    CONTINUE
   80    I=LAB(J)
         Y(J)=I
         K=J
         J=X(I)
         X(I)=K
         IF (I0.NE.I) GOTO 80
   90 CONTINUE
      H=0
      DO 100 I=1,N
         J=X(I)
         T=FIRST(I)-1
  101    T=T+1
         IF (T.GT.SZ) THEN
            PRINT '(A,I6,A)','MINPERM> WARNING D - ATOM ',I,' NOT MATCHED - MAXIMUM NUMBER OF NEIGHBOURS TOO SMALL?'
            RETURN
         ENDIF
         IF (KK(T).NE.J) GOTO 101
         DJ=CC(T)
         U(I)=DJ-V(J)
         H=H+DJ
  100 CONTINUE

 1000 END
