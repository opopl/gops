C     THIS IS AN UNTESTED WRAPPER TO MY Z2-ROUTINE THAT
C     USES LINKED CELLS TO EFFICIENTLY FIND WHICH ATOMS
C     ARE WITHIN INTERACTION RANGE OF A GIVEN ATOM.

      SUBROUTINE Z2FASTER(NP,XX,LX,LY,LZ,ETOT,GG,STEST)

C     FAIRLY EFFICIENT Z2 EVALUATOR FOR LARGE SYSTEMS
C     INPUT
C       NP -- NUMBER OF PARTICLES
C       XX -- COORDINATE VECTOR XX({1,2,3},I) = {X,Y,Z}(I)
C              OK TO CALL WITH XX OF DIMENSION XX(3*NP),
C              IN WHICH CASE {X,Y,Z}(I) = XX(3*(I-1) + {1,2,3})
C       LX -- BOX LENGTH IN X-DIRECTION
C       LY -- BOX LENGTH IN Y-DIRECTION
C       LZ -- BOX LENGTH IN Z-DIRECTION
C     OUTPUT
C       ETOT -- ENERGY OF CONFIGURATION
C       GG   -- FORCE (SAME DIMENSIONS AS XX)
C
      IMPLICIT NONE
      INTEGER NP
      REAL*8 XX(3,NP),LX,LY,LZ,ETOT,GG(3,NP)
      LOGICAL STEST

      REAL*8 RCUT,E,DE,BOXL(2,3)
      INTEGER I,NX,NY,NZ
      PARAMETER ( RCUT =  2.644877840738571D0   )

      IF (STEST) THEN
         PRINT '(A)','Z2FASTER> SECOND DERIVATIVES NOT CODED'
         STOP
      ENDIF

      NX = MAX(3,INT(LX/(RCUT+0.01D0)))
      NY = MAX(3,INT(LY/(RCUT+0.01D0)))
      NZ = MAX(3,INT(LZ/(RCUT+0.01D0)))
      BOXL(1,1) = LX
      BOXL(1,2) = LY
      BOXL(1,3) = LZ
      DO I=1,3
         BOXL(2,I) = 1.0D0 / BOXL(1,I)
      ENDDO

      CALL FGCALC(NP,NX,NY,NZ,XX,BOXL,E,GG,DE)
      ETOT = E+DE
      END


C     DIVIDE SPACE IN BOXES, AND COMPUTE WHICH PARTICLE BELONG
C     TO WHICH BOX. BOXES ARE SUPPOSEDLY HAVE SIDES .GE.
C     POTENTIAL REANGE, SO ONLY NEIGHBORING BOXES NEED BE CONSIDERED
C     FOR POTENTIAL EVALUATION.
      SUBROUTINE BUILD_BOXES(NP,NX,NY,NZ,XX,BOXL,FIRST,NEXT)
C     DIVIDE THE UNIT CUBE SIMULATION BOX INTO NX*NY*NZ CELLS,
C     AND PUT ALL PARTICLES INTO THEIR RESPECTIVE CELL.
      IMPLICIT NONE
      INTEGER NP,NX, NY,NZ,I,J,K,NB(3),IDX(3)
      INTEGER FIRST(NX,NY,NZ), NEXT(NP)
      REAL*8 XX(3,NP),BOXL(2,3)

      NB(1) = NX
      NB(2) = NY
      NB(3) = NZ

      DO I=1,NP
         NEXT(I) = 0
      ENDDO

      DO I=1,NB(3)
         DO J=1,NB(2)
            DO K=1,NB(1)
               FIRST(K,J,I) = 0
            ENDDO
         ENDDO
      ENDDO
      DO I=1,NP
         DO J=1,3
            IDX(J) = INT(XX(J,I)*BOXL(2,J)*NB(J))
            IF(XX(J,I) .LT. 0D0) THEN
               IDX(J) = 1 + MOD(IDX(J)-1-NB(J)*(IDX(J)-1),NB(J))
            ELSE
               IDX(J) = 1 + MOD(IDX(J),NB(J))
            ENDIF
            IF(IDX(J).LT.1 .OR. IDX(J).GT.NB(J)) THEN
               WRITE(*,*) I,J,XX(J,I),IDX(J),NB(J),BOXL(1,J),BOXL(2,J)
               STOP
            ENDIF
         ENDDO
         NEXT(I) = FIRST(IDX(1),IDX(2),IDX(3))
         FIRST(IDX(1),IDX(2),IDX(3)) = I
      ENDDO
      END

C     DISTANCE ROUTINE, TAKES PERIODIC BOUNDARY CONDITIONS INTO ACCOUNT
      SUBROUTINE Z2DIST(XX,YY,BOXL,DR)
      IMPLICIT NONE
      REAL*8 XX(3),YY(3),BOXL(2,3),DR(3)

      INTEGER I
      REAL*8 T

C     COMPUTE DISTANCE BETWEEN PARTICLES P AND Q
      DO I=1,3
         T = YY(I) - XX(I)
         DR(I) = T - BOXL(1,I)*NINT(T*BOXL(2,I))
      ENDDO
      END      

C     G77 DOES NOT HAVE CEILING IMPLEMENTED; THIS CAN BE SUBSTITUTED
C      REAL*8 FUNCTION CEILING(X)
C      IMPLICIT NONE
C      REAL*8 X,Y
C      Y = INT(X)
C      IF(Y.NE.X .AND. X.GE.0D0) Y = Y + 1D0
C      CEILING = Y
C      END

      REAL*8 FUNCTION EPAIR(X,Y,BOXL)
      IMPLICIT NONE
      REAL*8 X(3),Y(3),BOXL(2,3),DR(3),R
      REAL*8 RINV,S1,S2,S3,A,B,KF,ALPH,SIG,POW,ST,RCUT,RCUT2,TWOK,BSIG
      PARAMETER ( A    =  1.04D0      )
      PARAMETER ( B    =  4200000.0D0 )
      PARAMETER ( KF   =  4.139D0     )
      PARAMETER ( ALPH =  0.33D0      )
      PARAMETER ( SIG  =  0.348D0     )
      PARAMETER ( POW  = -14.5D0      )
C     PARAMETER ( ST=0.133915D0 ) ! VALUES BELOW OPTIMIZED BY DAVID
C     PARAMETER ( RCUT=2.645D0  ) ! WALES FOR EXACT ZERO AT RCUT
      PARAMETER ( ST   =  0.1339154253770228D0  )
      PARAMETER ( RCUT =  2.644877840738571D0   )
      PARAMETER ( RCUT2 =  RCUT**2              )
      PARAMETER ( TWOK  =  2.0D0 * KF           )
C     PARAMETER ( BSIG = B/SIG**POW             )
      PARAMETER ( BSIG  =  0.94656039153728358D0)

      CALL Z2DIST(X,Y,BOXL,DR)
      R = DR(1)**2 + DR(2)**2 + DR(3)**2
      IF(R .LT. RCUT2) THEN
C     ANALYTICAL POTENTIAL
         R = SQRT(R)
         RINV = 1.0D0 / R
         S1 = A * EXP(ALPH * R) * RINV**3
         S2 = COS(TWOK * R)
         S3 = R**POW
         EPAIR = S1*S2 + BSIG*S3 + ST
      ELSE
         EPAIR = 0D0
      ENDIF
      END


C     COMPUTE ENERGY AND GRADIENT. POTENTIAL ENERGY IS
C     E+DE ON OUTPUT. TYPICALLY, DE/E <~ 1E-16, SO THIS
C     SUM HOLDS ENERGY WITH MORE RESOLUTION THAN A SINGLE
C     64-BUT FLOATING POINT NUMBER.
C     ALGORITHM:
C     - DIVIDE SPACE IN BOXES, AND DISTRIBUTE PARTICLES.
C     - GO THROUGH BOXES, AND USE CURRENT BOX AND NEIGHBORS
C       TO FIND ALL PARTICLES INTERACTING WITH PARTICLES IN
C       CURRENT BOX.
C     - SUM ENERGY CONTRIBUTIONS WITH KAHAN SUMMATION, PRODUCING
C       HIGHER PRECISION. FOR THIS TO WORK, THE CODE NEED BE
C       COMPILED WITHOUT AGGRESSIVE OPTIMIZATINO (NO CHANGES OF
C       SEMANTICS, REORDERING OF CALCULATIONS ALLOWED). ALSO,
C       FLOATING POINT RESULT MUST BE STORED BACK IN MEMORY AFTER
C       EACH EXPRESSION, SO THAT EXCESS PRECISION POTENTIALLY
C       AVAILABLE IN THE CPU IS GOTTEN RID OF BEFORE NEXT EXPRESSION
C       IS EVALUATED. FOR G77, THIS IS ENABLED BY -FFLOAT-STORE,
C       AND FOR THE INTEL FORTRAN COMPILER WITH -FLTCONSISTENCY
C       USUALLY OPTIMIZATION LEVEL -O2 IS SAFE.
      SUBROUTINE FGCALC(NP,NX,NY,NZ,XX,BOXL,E,GG,DE)
      IMPLICIT NONE
      INTEGER NP,NX,NY,NZ
      REAL*8 XX(3,NP),BOXL(2,3),E,GG(3,NP),DE
      
      INTEGER FIRST(NX,NY,NZ),NEXT(NP),NB(3),DIVEC(3,13)
      INTEGER I,J,K,P,Q,II,I2,J2,K2,L
      REAL*8 V,R,GIJ,DR(3),DV,V1,T,H

C     UNCOMMENT IF G77, ALSO UNCOMMENT CEILING IMPLEMENTATION ABOVE
C      EXTERNAL CEILING
C      REAL*8 CEILING

C     FOR USE WHEN COMPUTING POTENTIAL BY TABLE LOOKUP
C      REAL*8 F
C      INTEGER IDX
C      INCLUDE 'TZ2C.H'


      REAL*8 RINV,S1,S2,S3,A,B,KF,ALPH,SIG,POW,ST,RCUT,RCUT2,TWOK,BSIG
      PARAMETER ( A    =  1.04D0      )
      PARAMETER ( B    =  4200000.0D0 )
      PARAMETER ( KF   =  4.139D0     )
      PARAMETER ( ALPH =  0.33D0      )
      PARAMETER ( SIG  =  0.348D0     )
      PARAMETER ( POW  = -14.5D0      )
C     PARAMETER ( ST=0.133915D0 ) ! VALUES BELOW OPTIMIZED BY DAVID
C     PARAMETER ( RCUT=2.645D0  ) ! WALES FOR EXACT ZERO AT RCUT
      PARAMETER ( ST   =  0.1339154253770228D0  )
      PARAMETER ( RCUT =  2.644877840738571D0   )

      PARAMETER ( RCUT2 =  RCUT**2              )
      PARAMETER ( TWOK  =  2.0D0 * KF           )
C     PARAMETER ( BSIG = B/SIG**POW             )
      PARAMETER ( BSIG  =  0.94656039153728358D0)


C     OFFSETS TO 13 UPPLEFT BOXES SURROUNDING EACH CELL.
      DATA DIVEC / -1,-1,-1 , -1,-1, 0 , -1,-1, 1 , -1, 0,-1
     :     ,       -1, 0, 0 , -1, 0, 1 , -1, 1,-1 , -1, 1, 0
     :     ,       -1, 1, 1 ,  0,-1,-1 ,  0,-1, 0 ,  0,-1, 1
     :     ,        0, 0,-1 /

      NB(1) = NX
      NB(2) = NY
      NB(3) = NZ
C     TO ENSURE THERE ARE NO NEGATIVE RESULTS IN MODULUS
C     COMPUTATIONS LATER ON, AND TO ACCOMODATE TO THE FACT
C     THAT FORTRAN INDICES START WITH ONE. DIVEC IS RESTORED
C     TO ITS ORIGINAL STATE UPON EXIT...
      DO I=1,13
         DO J=1,3
            DIVEC(J,I) = DIVEC(J,I) + NB(J) - 1
         ENDDO
      ENDDO

      DO I=1,NP
         DO J=1,3
            GG(J,I) = 0D0
         ENDDO
      ENDDO
      V = 0D0
      DV = 0D0

      CALL BUILD_BOXES(NP,NB(1),NB(2),NB(3),XX,BOXL,FIRST,NEXT)
      
      DO I=1,NB(3)
         DO J=1,NB(2)
            DO K=1,NB(1)
               IF(FIRST(K,J,I) .GT. 0) THEN
C     BEGIN PAIRS OF PARTICLES IN DIFFERENT BOXES

C     LOOP OVER UPPER LEFT 13 OF THE 26 SURROUNDING BOXES
                  DO II=1,13
                     I2 = MOD(I+DIVEC(3,II),NB(3)) + 1
                     J2 = MOD(J+DIVEC(2,II),NB(2)) + 1
                     K2 = MOD(K+DIVEC(1,II),NB(1)) + 1

                     IF(I.EQ.I2.AND.J.EQ.J2.AND.K.EQ.K2) THEN
                        WRITE(*,*)
     :                       'NEWSTWEB.F: BOX DOUBLE COUNTING: ',
     :                       I,J,K,II
                        STOP
                     ENDIF

                     IF(FIRST(K2,J2,I2) .GT. 0) THEN
                        P = FIRST(K,J,I)
                        DO WHILE(P .GT. 0)
                           Q = FIRST(K2,J2,I2)
                           DO WHILE(Q .GT. 0)
                              CALL Z2DIST(XX(1,P),XX(1,Q),BOXL,DR)
                              R = DR(1)**2 + DR(2)**2 + DR(3)**2
                              IF(R .LT. RCUT2) THEN
C                             ANALYTICAL POTENTIAL
                                 R = SQRT(R)
                                 RINV = 1.0D0 / R
                                 S1 = A * EXP(ALPH * R) * RINV**3
                                 S2 = COS(TWOK * R)
                                 S3 = R**POW
               
C                             CONTRIBUTION TO POTENTIAL ENERGY PHI(R)
C     V = V + S1*S2 + BSIG*S3 + ST
                                 H = S1*S2 + BSIG*S3 + ST
C     KAHAN SUMMATION, BROKEN UP SO -FFLOAT-STORE CAN DO ITS JOB
C     V1 = V+H+DV;
C     D  = -(V1-V)+H+DV;
C     V = V1;
                                 T = V+H
                                 V1 = T+DV
                                 T = -(V1-V)
                                 T = T+H
                                 DV = T+DV
                                 V = V1
C                             CONTRIBUTION TO GRADIENT 1/R * DPHI(R)/DR
                                 GIJ = (S1*((ALPH-3.0D0*RINV)*S2 -
     :                                TWOK*SIN(TWOK*R))
     :                                + POW*BSIG*S3*RINV)*RINV

C!C                             POTENTIAL EVALUATED FROM TABLE 
C!C                             CONTRIBUTION TO POTENTIAL ENERGY PHI(R)
C!                                 IDX = CEILING(R*IDRTAB);
C!                                 F = R*IDRTAB - IDX;
C!                                 V = V + POT(IDX) + F*POT1(IDX)
C!
C!C                             CONTRIBUTION TO GRADIENT 1/R * DPHI(R)/DR
C!                                 GIJ = DPOT(IDX) + F*DPOT1(IDX)
                                 DO L=1,3
                                    GG(L,P) = GG(L,P) - GIJ*DR(L)
                                    GG(L,Q) = GG(L,Q) + GIJ*DR(L)
                                 ENDDO
                              ENDIF
                              Q = NEXT(Q)
                           ENDDO
                           P = NEXT(P)
                        ENDDO
                     ENDIF
                  ENDDO
C     END

C     BEGIN PARTICLE PAIRS WITHIN BOX K,J,I
                  P = FIRST(K,J,I)
                  DO WHILE(P .GT. 0)
                     Q = NEXT(P)
                     DO WHILE(Q .GT. 0)
                        CALL Z2DIST(XX(1,P),XX(1,Q),BOXL,DR)
                        R = DR(1)**2 + DR(2)**2 + DR(3)**2
                        IF(R .LT. RCUT2) THEN
C                          ANALYTICAL POTENTIAL
                           R = SQRT(R)
                           RINV = 1.0D0 / R
                           S1 = A * EXP(ALPH * R) * RINV**3
                           S2 = COS(TWOK * R)
                           S3 = R**POW
               
C                         CONTRIBUTION TO POTENTIAL ENERGY PHI(R)
C     V = V + S1*S2 + BSIG*S3 + ST
                           H = S1*S2 + BSIG*S3 + ST
C     KAHAN SUMMATION, BROKEN UP SO -FFLOAT-STORE CAN DO ITS JOB
C     V1 = V+H+DV;
C     D  = -(V1-V)+H+DV;
C     V = V1;
                           T = V+H
                           V1 = T+DV
                           T = -(V1-V)
                           T = T+H
                           DV = T+DV
                           V = V1

C                         CONTRIBUTION TO GRADIENT 1/R * DPHI(R)/DR
                                 GIJ = (S1*((ALPH-3.0D0*RINV)*S2 -
     :                                TWOK*SIN(TWOK*R))
     :                                + POW*BSIG*S3*RINV)*RINV
C!C                       POTENTIAL EVALUATED FROM TABLE 
C!C                       CONTRIBUTION TO POTENTIAL ENERGY PHI(R)
C!                           IDX = CEILING(R*IDRTAB);
C!                           F = R*IDRTAB - IDX;
C!                           V = V + POT(IDX) + F*POT1(IDX)
C!                           
C!C                          CONTRIBUTION TO GRADIENT 1/R * DPHI(R)/DR
C!                           GIJ = DPOT(IDX) + F*DPOT1(IDX)
                           DO L=1,3
                              GG(L,P) = GG(L,P) - GIJ*DR(L)
                              GG(L,Q) = GG(L,Q) + GIJ*DR(L)
                           ENDDO
                        ENDIF
                        Q = NEXT(Q)
                     ENDDO
                     P = NEXT(P)
                  ENDDO
C     END
               ENDIF
            ENDDO
         ENDDO
      ENDDO

C     V1 = V+DV
C     DV  = -(V1-V)+DV
C     V  = V1
      V1 = V+DV
      T = V1-V
      DV = -T+DV
      V = V1

      E = V
      DE = DV

C     RESTORE TO ORGINAL STATE, SINCE DIVEC MIGHT BE SAVED
C     BETWEEN CALLS AND NOT REINITIALIZED ON THE NEXT ENTRY...
      DO I=1,13
         DO J=1,3
            DIVEC(J,I) = DIVEC(J,I) - NB(J) + 1
         ENDDO
      ENDDO
      END
