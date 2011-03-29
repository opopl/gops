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
C*******************************************************************************
C
C THIS SUBROUTINE CALCULATES THE AL/NI FARKAS ENERGY AND FIRST DERIVATIVE
C
C*******************************************************************************

      SUBROUTINE FARKAS(X,GRADFARK,EFARK,GRADT,NATOMS)
C      USE COMMONS
      IMPLICIT NONE
C     INTEGER AMAX

      INTEGER J1, J2, LFDEN, LFEMBED, LFPAIR, NATOMS
      DOUBLE PRECISION FDENX(500), FDENY(500), FDENY2(500),
     1                 FEMBEDX(500), FEMBEDY(500), FEMBEDY2(500),
     2                 FPAIRX(500), FPAIRY(500), FPAIRY2(500),
     3                 VPAIR, F, RHO(NATOMS), X(3*NATOMS), EFARK, 
     4                 RHOTEMP, DIST, RCUT, GRADFARK(3*NATOMS), DF(NATOMS), 
     5                 DVTEMP, DVTEMP1, DVTEMP2, DVTEMP3, DRHO,
     6                 DVPAIR, VPTEMP, RCUTSQ, DVPLO, RMIN
      LOGICAL GRADT
      COMMON /CFARKAS/ LFDEN, LFEMBED, LFPAIR, 
     1                FDENX, FDENY, FDENY2,
     2                FEMBEDX, FEMBEDY, FEMBEDY2,
     3                FPAIRX, FPAIRY, FPAIRY2

C CALCULATE ENERGY

      RCUT=FPAIRX(LFPAIR)
      RMIN=FPAIRX(1)
      RCUTSQ=RCUT**2
      CALL DSPLINT(FPAIRX,FPAIRY,FPAIRY2,LFPAIR,FPAIRX(1),DVPLO)
      EFARK=0.0D0
      DO 22 J1=1,NATOMS
         RHO(J1)=0.0D0
         VPAIR=0.0D0
         DO 23 J2=1,NATOMS
            IF (J1.NE.J2) THEN
               DIST=(X(3*(J2-1)+1)-X(3*(J1-1)+1))**2 +
     1              (X(3*(J2-1)+2)-X(3*(J1-1)+2))**2 +
     2              (X(3*(J2-1)+3)-X(3*(J1-1)+3))**2
               IF (DIST.LT.RCUTSQ) THEN
                 DIST=DSQRT(DIST)
                 IF (DIST.GT.RMIN) THEN
                   CALL SPLINT(FDENX,FDENY,FDENY2,LFDEN,DIST,
     1                         RHOTEMP)
                 ELSE
                   RHOTEMP=FDENY(1)
                 ENDIF
                 RHO(J1)=RHO(J1)+RHOTEMP
                 IF (J1.LT.J2) THEN
                   IF (DIST.GT.RMIN) THEN
                     CALL SPLINT(FPAIRX,FPAIRY,FPAIRY2,LFPAIR,
     1                         DIST,VPTEMP)
                   ELSE
                     VPTEMP=FPAIRY(1)+(DIST-RMIN)*DVPLO
C                     PRINT*,'ED: R OUTSIDE RANGE', DIST
                   ENDIF
                   VPAIR=VPAIR+VPTEMP
                 ENDIF
               ENDIF
            ENDIF
23       CONTINUE
         IF (RHO(J1).LT.FEMBEDX(LFEMBED)) THEN
           CALL SPLINT(FEMBEDX,FEMBEDY,FEMBEDY2,LFEMBED,RHO(J1),F)
         ELSE
C           PRINT*, 'ED: RHO OUTSIDE OF RANGE', RHO(J1)
           F=FEMBEDY(LFEMBED)
         ENDIF
         EFARK=EFARK+VPAIR+F
22    CONTINUE

C      PRINT*,EFARK

C CALCULATE GRADIENT

      IF (.NOT.GRADT) RETURN

      DO 30 J1=1,NATOMS
        IF (RHO(J1).LT.FEMBEDX(LFEMBED)) THEN
          CALL DSPLINT(FEMBEDX,FEMBEDY,FEMBEDY2,LFEMBED,RHO(J1),
     1               DF(J1))
        ELSE
C          PRINT*, 'ED: RHO OUTSIDE OF RANGE', RHO(J1)
          DF(J1)=0.0D0
        ENDIF
30    CONTINUE    
      
      DO 50 J1=1,NATOMS
         DVTEMP1=0.0D0
         DVTEMP2=0.0D0
         DVTEMP3=0.0D0
         DO 40 J2=1,NATOMS
            IF (J1.NE.J2) THEN
               DIST=( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     1              ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     2              ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2
               IF (DIST.LT.RCUTSQ) THEN
                 DIST=DSQRT(DIST)
                 IF (DIST.GT.RMIN) THEN
                   CALL DSPLINT(FPAIRX,FPAIRY,FPAIRY2,LFPAIR,DIST,
     1                        DVPAIR)
                   CALL DSPLINT(FDENX,FDENY,FDENY2,LFDEN,DIST,DRHO)
                 ELSE
C                   PRINT*, 'ED: R OUTSIDE OF RANGE', DIST
                   DVPAIR=DVPLO 
                   DRHO=0.0D0
                 ENDIF
                 DVTEMP=(DVPAIR+(DF(J1)+DF(J2))*DRHO)/DIST
                 DVTEMP1=DVTEMP1+DVTEMP*(X(3*(J1-1)+1)-X(3*(J2-1)+1))
                 DVTEMP2=DVTEMP2+DVTEMP*(X(3*(J1-1)+2)-X(3*(J2-1)+2))
                 DVTEMP3=DVTEMP3+DVTEMP*(X(3*(J1-1)+3)-X(3*(J2-1)+3))
               ENDIF
            ENDIF
40       CONTINUE
         GRADFARK(3*(J1-1)+1)=DVTEMP1
         GRADFARK(3*(J1-1)+2)=DVTEMP2
         GRADFARK(3*(J1-1)+3)=DVTEMP3
50    CONTINUE

      RETURN
      END

C*******************************************************************************
C
C THIS SUBROUTINE INITIALIZES THE ARRAYS OF 2ND DERIVATIVES FOR THE 
C CUBIC SPLINE INTERPOLATIONS FOR NI.
C
C*******************************************************************************

      SUBROUTINE NIINIT

      IMPLICIT NONE

      INTEGER I, LFDEN, LFEMBED, LFPAIR
      DOUBLE PRECISION FDENX(500), FDENY(500), FDENY2(500),
     1                 FEMBEDX(500), FEMBEDY(500), FEMBEDY2(500),
     2                 FPAIRX(500), FPAIRY(500), FPAIRY2(500),
     3                 DYLO, DYHI
      COMMON /CFARKAS/ LFDEN, LFEMBED, LFPAIR, 
     1                FDENX, FDENY, FDENY2,
     2                FEMBEDX, FEMBEDY, FEMBEDY2,
     3                FPAIRX, FPAIRY, FPAIRY2

      OPEN(UNIT=8,FILE='NI.DEN',STATUS='OLD')
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*) DYLO, DYHI
      READ(8,*)
      DO 20 I=1,500
        READ(8,*,END=30) FDENX(I), FDENY(I)
20    CONTINUE
      PRINT*,'WARNING: INCREASE THE DIMENSIONS IN NIINIT'
      STOP

30    LFDEN=I-1
      CLOSE(UNIT=8)
      PRINT*,'NI.DEN HAS ', LFDEN, ' ENTRIES'
      PRINT*,'FIRST DERIVATIVES AT END POINTS ', DYLO, DYHI

      CALL SPLINEGMIN(FDENX,FDENY,LFDEN,DYLO,DYHI,FDENY2)

      OPEN(UNIT=8,FILE='NI.PAIR',STATUS='OLD')
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*) DYLO, DYHI
      READ(8,*)
      DO 70 I=1,500
        READ(8,*,END=80) FPAIRX(I), FPAIRY(I)
70    CONTINUE
      PRINT*,'WARNING: INCREASE THE DIMENSIONS IN NIINIT'
      STOP

80    LFPAIR=I-1
      CLOSE(UNIT=8)
      PRINT*,'NI.PAIR HAS ', LFPAIR, ' ENTRIES'
      PRINT*,'SECOND DERIVATIVES AT END POINTS ', DYLO, DYHI

      CALL SPLINEGMIN(FPAIRX,FPAIRY,LFPAIR,DYLO,DYHI,FPAIRY2)

      OPEN(UNIT=8,FILE='NI.EMBED',STATUS='OLD')
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*) DYLO, DYHI
      READ(8,*)
      DO 120 I=1,500
        READ(8,*,END=130) FEMBEDX(I), FEMBEDY(I)
120   CONTINUE
      PRINT*,'WARNING: INCREASE 500'
      STOP

130   LFEMBED=I-1
      CLOSE(UNIT=8)
      PRINT*,'NI.EMBED HAS ', LFEMBED, ' ENTRIES'
      PRINT*,'FIRST DERIVATIVES AT END POINTS ',DYLO,DYHI

      CALL SPLINEGMIN(FEMBEDX,FEMBEDY,LFEMBED,DYLO,DYHI,
     1            FEMBEDY2)

      RETURN
      END

C*******************************************************************************
C
C THIS SUBROUTINE INITIALIZES THE ARRAYS OF 2ND DERIVATIVES FOR THE 
C CUBIC SPLINE INTERPOLATIONS FOR AL.
C
C*******************************************************************************

      SUBROUTINE ALINIT

      IMPLICIT NONE

      INTEGER I, LFDEN, LFEMBED, LFPAIR
      DOUBLE PRECISION FDENX(500), FDENY(500), FDENY2(500),
     1                 FEMBEDX(500), FEMBEDY(500), FEMBEDY2(500),
     2                 FPAIRX(500), FPAIRY(500), FPAIRY2(500), 
     3                 DYLO, DYHI
      COMMON /CFARKAS/ LFDEN, LFEMBED, LFPAIR, 
     1                FDENX, FDENY, FDENY2,
     2                FEMBEDX, FEMBEDY, FEMBEDY2,
     3                FPAIRX, FPAIRY, FPAIRY2

      OPEN(UNIT=8,FILE='AL.DEN',STATUS='OLD')
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*) DYLO, DYHI
      READ(8,*)
      DO 20 I=1,500
        READ(8,*,END=30) FDENX(I), FDENY(I)
20    CONTINUE
      PRINT*,'WARNING: INCREASE THE DIMENSIONS IN ALINIT'
      STOP

30    LFDEN=I-1
      CLOSE(UNIT=8)
      PRINT*,'AL.DEN HAS ', LFDEN, ' ENTRIES'
      PRINT*,'FIRST DERIVATIVES AT END POINTS ',DYLO,DYHI

      CALL SPLINEGMIN(FDENX,FDENY,LFDEN,DYLO,DYHI,FDENY2)

      OPEN(UNIT=8,FILE='AL.PAIR',STATUS='OLD')
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*) DYLO, DYHI
      READ(8,*)
      DO 70 I=1,500
        READ(8,*,END=80) FPAIRX(I), FPAIRY(I)
70    CONTINUE
      PRINT*,'WARNING: INCREASE THE DIMENSIONS IN ALINIT'
      STOP

80    LFPAIR=I-1
      CLOSE(UNIT=8)
      PRINT*,'AL.PAIR HAS ', LFPAIR, ' ENTRIES'
      PRINT*,'FIRST DERIVATIVES AT END POINTS ',DYLO,DYHI

      CALL SPLINEGMIN(FPAIRX,FPAIRY,LFPAIR,DYLO,DYHI,FPAIRY2)

      OPEN(UNIT=8,FILE='AL.EMBED',STATUS='OLD')
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*) DYLO, DYHI
      READ(8,*)
      DO 120 I=1,500
        READ(8,*,END=130) FEMBEDX(I), FEMBEDY(I)
120   CONTINUE
      PRINT*,'WARNING: INCREASE THE DIMENSIONS IN ALINIT'
      STOP

130   LFEMBED=I-1
      CLOSE(UNIT=8)
      PRINT*,'AL.EMBED HAS ', LFEMBED, ' ENTRIES'
      PRINT*,'FIRST DERIVATIVES AT END POINTS ',DYLO,DYHI

      CALL SPLINEGMIN(FEMBEDX,FEMBEDY,LFEMBED,DYLO,DYHI,FEMBEDY2)

      RETURN
      END

      SUBROUTINE SPLINEGMIN(X,Y,N,YP1,YPN,Y2)
      INTEGER N,NMAX
      DOUBLE PRECISION YP1,YPN,X(N),Y(N),Y2(N)
      PARAMETER (NMAX=5000)
      INTEGER I,K
      DOUBLE PRECISION P,QN,SIG,UN,U(NMAX)
      IF (YP1.GT..99D30) THEN
        Y2(1)=0.D0
        U(1)=0.D0
      ELSE
        Y2(1)=-0.5D0
        U(1)=(3.D0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.D0
        Y2(I)=(SIG-1.D0)/P
        U(I)=(6.D0*((Y(I+1)-Y(I))/(X(I+
     *1)-X(I))-(Y(I)-Y(I-1))/(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*
     *U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99D30) THEN
        QN=0.D0
        UN=0.D0
      ELSE
        QN=0.5D0
        UN=(3.D0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.D0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      INTEGER N
      DOUBLE PRECISION X,Y,XA(N),Y2A(N),YA(N)
      INTEGER KHI,KLO
      DOUBLE PRECISION A,B,H,POS

      IF ((X.LT.XA(1)).OR.(X.GT.XA(N))) THEN
        PRINT*, 'WARNING: X OUT OF RANGE IN SPLINT', X, XA(1),XA(N)
        STOP
      ENDIF

C WORKS IF INTERVALS EQUALLY SPACED 

      POS=1+(N-1)*(X-XA(1))/(XA(N)-XA(1))
      KLO=INT(POS)
      KHI=INT(POS)+1

C CHECK IF ON END POINT: UNLIKELY BUT POSSIBLE

      IF (KLO.EQ.N) THEN
        KLO=N-1
        KHI=N
      ENDIF

C INTERVAL BISECTION IF NOT
C      KLO=1
C      KHI=N
C1     IF (KHI-KLO.GT.1) THEN
C        K=(KHI+KLO)/2
C        IF(XA(K).GT.X)THEN
C          KHI=K
C        ELSE
C          KLO=K
C        ENDIF
C      GOTO 1
C      ENDIF

      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.D0) PRINT*, 'BAD XA INPUT IN SPLINT'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**
     *2)/6.D0
      RETURN
      END

      SUBROUTINE DSPLINT(XA,YA,Y2A,N,X,DYDX)
      INTEGER N
      DOUBLE PRECISION X,DYDX,XA(N),Y2A(N),YA(N)
      INTEGER KHI,KLO
      DOUBLE PRECISION A,B,H,POS

      IF ((X.LT.XA(1)).OR.(X.GT.XA(N))) THEN
        PRINT*, 'WARNING: X OUT OF RANGE IN DSPLINT', X, XA(1),XA(N)
        STOP
      ENDIF

C WORKS IF INTERVALS EQUALLY SPACED 

      POS=1+(N-1)*(X-XA(1))/(XA(N)-XA(1))
      KLO=INT(POS)
      KHI=INT(POS)+1

C CHECK IF ON END POINT: UNLIKELY BUT POSSIBLE

      IF (KLO.EQ.N) THEN
        KLO=N-1
        KHI=N
      ENDIF

C INTERVAL BISECTION IF NOT
C      KLO=1
C      KHI=N
C1     IF (KHI-KLO.GT.1) THEN
C        K=(KHI+KLO)/2
C        IF(XA(K).GT.X)THEN
C          KHI=K
C        ELSE
C          KLO=K
C        ENDIF
C      GOTO 1
C      ENDIF

      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.D0) PRINT*, 'BAD XA INPUT IN SPLINT'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      DYDX=(YA(KHI)-YA(KLO))/(XA(KHI)-XA(KLO))-
     1     (3*A**2-1)*(XA(KHI)-XA(KLO))*Y2A(KLO)/6.0D0+
     2     (3*B**2-1)*(XA(KHI)-XA(KLO))*Y2A(KHI)/6.0D0
      RETURN
      END
