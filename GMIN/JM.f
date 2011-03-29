C   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF OPTIM.
C
C   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C
C
C*************************************************************************
C
C  TWO- AND THREE-BODY TERMS IN THE ENERGY OF THE JM POTENTIAL WITH
C  CUTOFF - NOT PERIODIC. THIS MUST BE CALLED BEFORE JM2C
C  OR JM3C AS IT DOES SOME SETTING UP FOR THEM.
C                                        
C*************************************************************************
C
      SUBROUTINE JMEC(N,X,VNEW,POTEL,GTEST,STEST)
      IMPLICIT NONE

      INTEGER N, J1, J2, I, J, IJ
      LOGICAL YESNO, GTEST,STEST
      INTEGER M(N,N,3)
      DOUBLE PRECISION CUTOFF, X(3*N), POTEL, P2, P3, RAB, RHOAB, RBC, 
     1                 RAC, QQ2, RHOBC, RHOAC, QA, QB, QC, QQ3
      DOUBLE PRECISION DIST(N,N), SDIST(N,N), VNEW(3*N), 
     1              VEC(N,N,3), RH(N,N), NDUM(N,N),
     2              SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3

      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3, CUTOFF
      LOGICAL :: CALLED=.FALSE.
C
      SR2=1.414213562D0
      SR3=1.732050808D0
      SR6=2.449489743D0
      P2=0.0D0
      P3=0.0D0
      IF (.NOT.CALLED) THEN
         CUTOFF=1.0D10
         INQUIRE(FILE='JMPARAMS',EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            PRINT*,'DATA FILE JMPARAMS NOT FOUND - QUIT'
            STOP
         ELSE
            OPEN(UNIT=33,FILE='JMPARAMS',STATUS='OLD')
            READ(33,*) C0
            READ(33,*) C1
            READ(33,*) C2
            READ(33,*) C3
            READ(33,*) C4
            READ(33,*) C5
            READ(33,*) C6
            READ(33,*) C7
            READ(33,*) C8
            READ(33,*) C9
            READ(33,*) C10
            READ(33,*) RE
            READ(33,*) D
            READ(33,*) AN2
            READ(33,*) AN3
            READ(33,*,END=666) CUTOFF
666      CLOSE(33)
            CALLED=.TRUE.
         ENDIF
      ENDIF
C
C  CALCULATION OF CONNECTING VECTORS
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  CALCULATION OF DISTANCES:
C
      DO 20 J1=1,N
         DIST(J1,J1)=0.0D0
         SDIST(J1,J1)=0.0D0
         DO 10 J2=J1+1,N
            DIST(J1,J2)= VEC(J2,J1,1)**2 +
     1                   VEC(J2,J1,2)**2 +
     2                   VEC(J2,J1,3)**2 
            DIST(J1,J2)=DSQRT(DIST(J1,J2))
            DIST(J2,J1)=DIST(J1,J2)
            SDIST(J1,J2)=1/DIST(J1,J2)
            SDIST(J2,J1)=SDIST(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  CALCULATE THE ENERGY
C
      DO 22 I=1,N
         DO 23 J=1,N
            IF (I.NE.J.AND.DIST(I,J).LT.CUTOFF) THEN 
               RAB=DIST(I,J)
               RHOAB=(RAB-RE)/RE
               P2=P2-D*(1+AN2*RHOAB)*DEXP(-AN2*RHOAB)
               DO 24 IJ=1,N
                  IF ((IJ.NE.I).AND.(IJ.NE.J).AND.
     1                DIST(IJ,I).LT.CUTOFF.AND.
     2                DIST(IJ,J).LT.CUTOFF) THEN
                     RBC=DIST(IJ,J)
                     RAC=DIST(IJ,I)
                     RHOBC=(RBC-RE)/RE
                     RHOAC=(RAC-RE)/RE
                     QA=(RHOAB+RHOBC+RHOAC)/SR3
                     QQ2=(RHOBC-RHOAC)/SR2
                     QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                     QB=QQ2**2+QQ3**2
                     QC=QQ3**3-3*QQ3*QQ2**2
                     P3=P3+
     1                 (C0+C1*QA+C2*(QA**2)+C3*(QB)+C4*(QA**3)+
     2                  C5*QA*QB+C6*QC+C7*(QA**4) +
     3                  C8*(QA**2)*QB+C9*(QB**2)+C10*QA*QC)
     4                *D*DEXP(-AN3*QA)
                   ENDIF
24              CONTINUE
             ENDIF
23       CONTINUE
22    CONTINUE
      P3=P3/6.0D0
      P2=P2/2.0D0
      POTEL=P2+P3
      IF (GTEST.OR.STEST) THEN
         CALL JM2C(N, X, VNEW, DIST, SDIST, VEC, RH, M, NDUM,STEST)
         CALL JM3C(N, X, VNEW, DIST, SDIST, VEC, RH, M, NDUM,STEST)
      ENDIF
      RETURN
      END
C
C*************************************************************************
C
C  TWO BODY TERM OF MURRELL POTENTIAL - ANALYTIC DERIVATIVES WITH
C  PERIODIC BOUNDARY CONDITIONS
C
C*************************************************************************
C
      SUBROUTINE JM2C (N, X, V, RAB, RRAB, VEC, RHOAB, M, NDUM, STEST)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      INTEGER M(N,N,3)
      LOGICAL STEST
      DOUBLE PRECISION RAB(N,N), RRAB(N,N), 
     1              VEC(N,N,3), RHOAB(N,N), NDUM(N,N)

      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3, CUTOFF
      DOUBLE PRECISION CUTOFF, SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), ABY, ABX
C 
C
C  STORE DISTANCE MATRICES.
C
      DO 20 J1=1,N
         DO 10 J2=J1+1,N
            IF (RAB(J2,J1).GE.CUTOFF) THEN
               RRAB(J2,J1)=0.0D0
               RRAB(J1,J2)=0.0D0
            ENDIF
            RHOAB(J2,J1)=(RAB(J2,J1)-RE)/RE
            RHOAB(J1,J2)=RHOAB(J2,J1)
10       CONTINUE
20    CONTINUE
C
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
                  ABX=VEC(J4,J1,J2)
                  V(J3)=V(J3)
     1              +AN2*ABX*D*(RRAB(J4,J1)-
     2              (1+AN2*RHOAB(J4,J1))*RRAB(J4,J1)) 
     3              /(DEXP(AN2*RHOAB(J4,J1))*RE) 
30          CONTINUE
C           PRINT*,'J3,V(J3)=',J3,V(J3)
40       CONTINUE
50    CONTINUE
      RETURN
      END
C
C*************************************************************************
C
C  HERE WE CALCULATE THE ANALYTIC GRADIENT AND SECOND DERIVATIVES
C  FOR THE THREE-BODY TERM WITH PERIODIC BOUNDARY CONDITIONS.
C                                        
C*************************************************************************
C
      SUBROUTINE JM3C (N, X, V, R2, RR2, VEC, RH, M, NDUM, STEST)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5
      INTEGER M(N,N,3)
      LOGICAL STEST
      DOUBLE PRECISION R2(N,N), RR2(N,N), 
     1              VEC(N,N,3), RH(N,N), NDUM(N,N)
      COMMON /PARAMS/ SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3, CUTOFF
      DOUBLE PRECISION X(3*N), ABX, ACX, BCX, V(3*N),
     1                 TEMP, RRBC,
     2                 RAB, RRAB, RAC, RRAC, RBC,
     3                 ABY, ACY, BCY, RHOAB, RHOBC ,RHOAC ,QA, QQ2,
     4                 QQ3, QB, QC, 
     5                 SR2, SR3, SR6, C0, C1, C2, C3, C4,
     1                 C5, C6, C7, C8, C9, C10, RE, D, AN2, AN3
      DOUBLE PRECISION CUTOFF
C
C
C  FIRST THE GRADIENT.
 
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               IF (J3.NE.J1.AND.R2(J3,J1).LT.CUTOFF) THEN
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABX=VEC(J3,J1,J2)
               RHOAB=(RAB-RE)/RE
               DO 95 J4=J3+1,N
                  IF (J4.NE.J1.AND.R2(J4,J1).LT.CUTOFF
     1               .AND.(R2(J4,J3).LT.CUTOFF)) THEN
                  BCX=VEC(J4,J3,J2)
                  ACX=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1) 
                  RHOAC=(RAC-RE)/RE
                  RHOBC=(RBC-RE)/RE
                  QA=(RHOAB+RHOAC+RHOBC)/SR3
                  QQ3=(2*RHOAB-RHOBC-RHOAC)/SR6
                  QQ2=(RHOBC-RHOAC)/SR2
                  QB=(QQ2**2)+(QQ3**2)
                  QC=QQ3**3-3*QQ3*QQ2**2
                  TEMP=TEMP+
     1           D*((C1+3*C4*QA**2+4*C7*QA**3+C5*QB+2*QA*(C2+C8*QB)+
     2           C10*QC-AN3*(C0+C4*QA**3+C7*QA**4+C3*QB+C9*QB**2+
     3           QA**2*(C2+C8*QB)+C6*QC+QA*( C1+C5*QB+C10*QC)))*
     4           (-(ABX*RRAB/RE)-ACX*RRAC/RE)/SR3+
     5           (2*(C3+C5*QA+C8*QA**2+2*C9*QB)*
     6           (ACX*QQ2*RRAC/SR2+QQ3*(-2*ABX*RRAB+ACX*RRAC)/SR6)+
     7           (C6+C10*QA)*(-6*ACX*QQ2*QQ3*RRAC/SR2+
     8           (-3*QQ2**2+3*QQ3**2)*(-2*ABX*RRAB+ACX*RRAC)/SR6))/RE)/
     9           DEXP(AN3*QA)
               ENDIF
95             CONTINUE
            ENDIF
100         CONTINUE
            V(3*(J1-1)+J2)=V(3*(J1-1)+J2)+TEMP
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
      RETURN
      END
