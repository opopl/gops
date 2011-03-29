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
      SUBROUTINE MLATMIN(N,XC,RHO,BOXLX,BOXLY,BOXLZ,CUTOFF)
      IMPLICIT NONE
      INTEGER I, J, J1, NCOUNT, N, INFO
      DOUBLE PRECISION XC(3*N),CUTOFF,F1,GRAD(3),DIAG(3), TEMPA(9),
     1                 SECOND(3,3), XSAVE(3*N), BXSAVE, 
     2                 PSTEP(3), STEP(3), FOB(3), BYSAVE, BZSAVE,
     3                 RHO, BOXLX, BOXLY, BOXLZ, LP, LP1, LP2, STPMAG
      LOGICAL NRTEST

C
C  VALUE OF DIF IS THE ORDER OF MAGNITUDE TO WHICH THE LATTICE
C  CONSTANT CAN BE OPTIMISED. SETTING IT SMALLER THAN 10^(-7)
C  CAUSES NUMERICAL PROBLEMS ON THE DEC.
C
      NCOUNT=1
      BXSAVE=BOXLX
      BYSAVE=BOXLY
      BZSAVE=BOXLZ
      DO 20 J1=1,3*N
         XSAVE(J1)=XC(J1)
20    CONTINUE

10    DO 25 J1=1,3*N
         IF (MOD(J1,3).EQ.1) XC(J1)=XSAVE(J1)*BOXLX/BXSAVE
         IF (MOD(J1,3).EQ.2) XC(J1)=XSAVE(J1)*BOXLY/BYSAVE
         IF (MOD(J1,3).EQ.0) XC(J1)=XSAVE(J1)*BOXLZ/BZSAVE
25    CONTINUE

      PRINT*,'BOXLENGTH IN X IS ',BOXLX
      PRINT*,'BOXLENGTH IN Y IS ',BOXLY
      PRINT*,'BOXLENGTH IN Z IS ',BOXLZ

      CALL MLATDIFF(N,XC,F1,GRAD,SECOND,RHO,BOXLX,BOXLY,BOXLZ,CUTOFF)

C      PRINT*,'SECOND DERIVATIVE MATRIX IS ',(SECOND(1,J),J=1,3)
C      PRINT*,'SECOND DERIVATIVE MATRIX IS ',(SECOND(2,J),J=1,3)
C      PRINT*,'SECOND DERIVATIVE MATRIX IS ',(SECOND(3,J),J=1,3)

C     CALL TRED2(SECOND,3,3,DIAG,E,1)
C     CALL TQLI(DIAG,E,3,3,SECOND,1)

      CALL DSYEV('V','U',3,SECOND,3,DIAG,TEMPA,9,INFO)
      IF (INFO.NE.0) THEN
         PRINT*,'WARNING - INFO=',INFO,' IN DSYEV'
      ENDIF

      PRINT*
      PRINT*,'ENERGY FOR LATTICE CYCLE ',NCOUNT,' IS ',F1
      PRINT*,'GRADIENT WRT BOX LENGTH X =',GRAD(1)
      PRINT*,'GRADIENT WRT BOX LENGTH Y =',GRAD(2)
      PRINT*,'GRADIENT WRT BOX LENGTH Z =',GRAD(3)
      PRINT*,'EIGENVALUES OF SECOND DERIVATIVE MATRIX WRT BOX LENGTH='
      PRINT*,DIAG(1)
      PRINT*,DIAG(2)
      PRINT*,DIAG(3)
C      PRINT*,'EIGENVECTOR MATRIX IS ',(SECOND(1,J),J=1,3)
C      PRINT*,'EIGENVECTOR MATRIX IS ',(SECOND(2,J),J=1,3)
C      PRINT*,'EIGENVECTOR MATRIX IS ',(SECOND(3,J),J=1,3)
      NCOUNT=NCOUNT+1

C  CONVERT GRADIENT TO EIGENVECTOR BASIS

      DO 140 I=1,3
         FOB(I)=0.0D0
         DO 130 J = 1, 3
            FOB(I)=FOB(I)+GRAD(J)*SECOND(J,I)
130      CONTINUE
140   CONTINUE

C CALCULATE STEP IN EIGENVECTOR BASIS
      NRTEST=.FALSE.
C      NRTEST=.TRUE.
      DO 240 I=1,3
         PSTEP(I)=0.0D0
         LP1=DABS(DIAG(I))/2.0D0
         LP2=1.0D0 + 4.0D0*(FOB(I)/DIAG(I))**2
C EIGENVECTOR FOLLOWING STEP
         LP=LP1*(1.0D0+DSQRT(LP2))
C NEWTON RAPHSON STEP
         IF (NRTEST) LP=DIAG(I)
         PSTEP(I)=-FOB(I)/LP
240   CONTINUE

      DO 350 J=1,3
         STEP(J)=0.0D0
         DO 340 I=1,3
            STEP(J)=STEP(J)+PSTEP(I)*SECOND(J,I)
340      CONTINUE
350   CONTINUE

      STPMAG=DSQRT(STEP(1)**2+STEP(2)**2+STEP(3)**2)
      DO 370 I=1,3
        IF (STPMAG.GT.0.1D0) STEP(I)=STEP(I)*0.1D0/STPMAG
370   CONTINUE
C      PRINT*, 'SIZE OF LATTICE STEP IS ',STPMAG

      PRINT*,'STEP IN X DIRECTION', STEP(1)
      PRINT*,'STEP IN Y DIRECTION', STEP(2)
      PRINT*,'STEP IN Z DIRECTION', STEP(3)
      
      BOXLX=BOXLX+STEP(1)
      BOXLY=BOXLY+STEP(2)
      BOXLZ=BOXLZ+STEP(3)

      IF (STPMAG.GT.0.1D0) STPMAG=0.1D0

      IF (STPMAG.GT.1.0D-6) GOTO 10

      DO 70 J1=1,3*N
         IF (MOD(J1,3).EQ.1) XC(J1)=XSAVE(J1)*BOXLX/BXSAVE
         IF (MOD(J1,3).EQ.2) XC(J1)=XSAVE(J1)*BOXLY/BYSAVE
         IF (MOD(J1,3).EQ.0) XC(J1)=XSAVE(J1)*BOXLZ/BZSAVE
70    CONTINUE

      RETURN
      END
C
C*************************************************************************
C
C  SUBROUTINE MLATDIFF CALCULATES THE CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE WITH RESPECT TO THE BOX SIZE ANALYTICALLY FOR THE 
C  MORSE POTENTIAL.
C
C*************************************************************************
C
      SUBROUTINE MLATDIFF(N,X,P2,VLAT,ALAT,RHO,BOXLX,BOXLY,BOXLZ,CUTOFF)
      IMPLICIT NONE 
      INTEGER N, J1, J2, I, J, K, L, MX, MY, MZ 
      DOUBLE PRECISION X(3*N), RHO, BOXL(3),
     1                 DIST, R(N,N), P2, MEXP, VTEMP,
     2                 BOXLX, BOXLY, BOXLZ, VEC(N,N,3), 
     3                 CUTOFF, VLAT(3), ALAT(3,3), RCUT
C
C  DEAL WITH ATOMS LEAVING THE BOX:
C
      IF ((BOXLX.LT.BOXLY).AND.(BOXLX.LT.BOXLZ)) THEN
        RCUT=BOXLX*CUTOFF
      ELSE IF (BOXLY.LT.BOXLZ) THEN
        RCUT=BOXLY*CUTOFF
      ELSE
        RCUT=BOXLZ*CUTOFF
      ENDIF
      PRINT*,'CUTOFF USED = ',RCUT

      DO 41 J1=1,N
         X(3*(J1-1)+1)=X(3*(J1-1)+1) - BOXLX*DNINT(X(3*(J1-1)+1)/BOXLX)  
         X(3*(J1-1)+2)=X(3*(J1-1)+2) - BOXLY*DNINT(X(3*(J1-1)+2)/BOXLY)
         X(3*(J1-1)+3)=X(3*(J1-1)+3) - BOXLZ*DNINT(X(3*(J1-1)+3)/BOXLZ)
41    CONTINUE
C
C  CALCULATION OF CONNECTING VECTORS; TO IMPLEMENT THE PERIODIC
C  BOUNDARY CONDITIONS, THE SHORTEST VECTOR BETWEEN TWO ATOMS IS
C  USED:
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            MX=NINT(VEC(J2,J1,1)/BOXLX)
            MY=NINT(VEC(J2,J1,2)/BOXLY)
            MZ=NINT(VEC(J2,J1,3)/BOXLZ)
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * MX
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * MY
            VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * MZ
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  STORE DISTANCE MATRICES.
C
      P2=0.0D0
      DO 5 K=1,3
        VLAT(K)=0.0D0
        DO 4 L=1,3
          ALAT(K,L)=0.0D0
4       CONTINUE
5     CONTINUE
      BOXL(1)=BOXLX
      BOXL(2)=BOXLX
      BOXL(3)=BOXLZ
      
      DO 20 J1=1,N
         R(J1,J1)=0.0D0
         DO 10 J2=1,J1-1
            DIST=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 +
     1           VEC(J1,J2,3)**2
            R(J2,J1)=DSQRT(DIST)
            R(J1,J2)=R(J2,J1)
10       CONTINUE
20    CONTINUE
C
C  HERE WE CALCULATE THE ENERGY
C                                        
      DO 22 I=1,N-1
         DO 23 J=I+1,N
           IF (R(I,J).LT.RCUT) THEN
            MEXP=DEXP(RHO*(1.0D0-R(I,J)))
            P2=P2+MEXP*(MEXP-2.0D0)
            VTEMP=2.0D0*RHO*MEXP*(1.0D0-MEXP)/R(I,J)
            DO 21 K=1,3
             VLAT(K)=VLAT(K)+VTEMP*VEC(I,J,K)**2/BOXL(K)
             ALAT(K,K)=ALAT(K,K)+2.0D0*RHO*MEXP*VEC(I,J,K)**2/
     1                (R(I,J)*BOXL(K)**2)*((VEC(I,J,K)**2/R(I,J)**2-
     2                 1.0D0)*(MEXP-1.0D0)+(2.0D0*MEXP-1.0D0)*RHO*
     3                 VEC(I,J,K)**2/R(I,J))
             DO 18 L=K+1,3
             ALAT(K,L)=ALAT(K,L)+2.0D0*RHO*VEC(I,J,K)**2*VEC(I,J,L)**2*
     1                   MEXP/(BOXL(K)*BOXL(L)*R(I,J)**2)*(RHO*
     2                   (2.0D0*MEXP-1.0D0)+(MEXP-1.0D0)/R(I,J))
18           CONTINUE
21          CONTINUE
           ENDIF
23       CONTINUE
22    CONTINUE

      ALAT(2,1)=ALAT(1,2)
      ALAT(3,1)=ALAT(1,3)
      ALAT(3,2)=ALAT(2,3)

      RETURN
      END

