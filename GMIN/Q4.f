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
C******************************************************************************
C
C  SUBROUTINE ORDER CALCULATES Q4 ONLY
C
C******************************************************************************

      SUBROUTINE ORDERQ4(NATOMS,POINTS,CURQ4)

      IMPLICIT DOUBLE PRECISION(A-H,P-Z)
      INTEGER J, K, NB, NATOMS
      DOUBLE PRECISION DX,DY,DZ,DIST,POINTS(3*NATOMS),
     1                 Q40,Q4R(4),Q4I(4), CURQ4

      Q40=0.0D0
      DO 10 J=1,4
        Q4R(J)=0.0D0
        Q4I(J)=0.0D0
10    CONTINUE

      NB=0
      DO 30 J=1,NATOMS
         DO 40 K=J+1,NATOMS
            DX=POINTS(3*(J-1)+1)-POINTS(3*(K-1)+1)
            DY=POINTS(3*(J-1)+2)-POINTS(3*(K-1)+2)
            DZ=POINTS(3*(J-1)+3)-POINTS(3*(K-1)+3)
            DIST=DSQRT(DX*DX+DY*DY+DZ*DZ)
            IF (DIST.LT.1.3909D0) THEN
               NB=NB+1
               CALL EVASH4(DX/DIST,DY/DIST,DZ/DIST,Q40,Q4R,Q4I)
            ENDIF
40       CONTINUE
30    CONTINUE

      IF (NB.EQ.0) PRINT*,'*** WARNING - NB=0'

      CURQ4 = Q4(NB,Q40,Q4R,Q4I)

      RETURN
      END

C
C**************************************************************
C     
      SUBROUTINE SHINIT
 
C *** CALCULATE COEFFICIENTS FOR Q4, Q6, W4 AND W6
C     NOTEBOOK JVD A92
 
      COMMON/Q4COEF/ Q4COEF
 
      INTEGER    M, I
      DOUBLE PRECISION       FAC, FACS
 
      DOUBLE PRECISION      Q4COEF(3,0:4)
 
C *** IT WOULD BE NICE INDEED TO KNOW THE EXACT FORMULA!
 
C     DATA      Q4COEF /   4.375,  -3.75,  0.375,
C    :                   -17.5,     7.5,   0.0,
C    :                   -52.5,    60.0,  -7.5,
C    :                  -105.0,     0.0,   0.0,
C    :                   105.0,     0.0,   0.0 /
 
C *** NOW MULTIPLY THE COEFFICIENTS WITH SQRT(2.0* (L-M)!/(L+M)! )

      Q4COEF(1,0)=4.375D0
      Q4COEF(2,0)=-3.75D0
      Q4COEF(3,0)=0.375D0
      Q4COEF(1,1)=-17.5D0
      Q4COEF(2,1)=7.5D0
      Q4COEF(3,1)=0.0D0
      Q4COEF(1,2)=-52.5D0
      Q4COEF(2,2)=60.0D0
      Q4COEF(3,2)=-7.5D0
      Q4COEF(1,3)=-105.0D0
      Q4COEF(2,3)=0.0D0
      Q4COEF(3,3)=0.0D0
      Q4COEF(1,4)=105.0D0
      Q4COEF(2,4)=0.0D0
      Q4COEF(3,4)=0.0D0
 
      DO 99 M=1, 4
         DO 98 I=1, 3
            FACS=DSQRT(2.0*FAC(4-M)/FAC(4+M))
            Q4COEF(I,M)=Q4COEF(I,M)*FACS
98       CONTINUE
99    CONTINUE
 
 
      DO 89 M=1, 6
         DO 88 I=1, 4
            FACS=DSQRT(2.0*FAC(6-M)/FAC(6+M))
C           Q6COEF(I,M)=Q6COEF(I,M)*FACS
88       CONTINUE
89    CONTINUE
 
      RETURN
      END
 
      SUBROUTINE EVASH4(X, Y, Z, Q0, QR, QI)
 
C *** EVALUATE THE SPHERICAL HARMONICS OF DEGREE 4 ********************
C
C     REAL*8   X,Y,Z    A VECTOR ON THE UNIT-SPHERE TO BE PROCESSED
C     DOUBLE PRECISION
C            Q0       ACCUMULATOR OF P40(COS(THETA))
C            QR(4)    ACCUMULATOR OF P4M(COS(THETA))*COS(PHI)
C            QI(4)    ACCUMULATOR OF P4M(COS(THETA))*SIN(PHI)
C
C     NOTE THAT THE ACTUAL SPHERICAL HARMONICS Y4M DIFFER FROM Q0
C     AND Q (= QR + I*QI) BY A FACTOR
C                     SQRT((9/8*PI)
C     THIS WILL BE TAKEN CARE OF IN FUNCTION Q4
C
C     ANGLES THETA AND PHI ARE DEFINED IN THE REGULAR WAY:
C
C         X = COS(PHI)SIN(THETA)
C         Y = SIN(PHI)SIN(THETA)
C         Z = COS(THETA)
C
C     NOTEBOOK JVD A38 AND A85
C
C *********************************************************************
 
      COMMON/Q4COEF/ Q4COEF
      DOUBLE PRECISION      Q4COEF(3,0:4)
 
      DOUBLE PRECISION      X, Y, Z
      DOUBLE PRECISION  Q0, QR(4), QI(4)
 
      DOUBLE PRECISION COSTH, COSTH2, COSTH4, SINTH, SINTH2, SCTH
      DOUBLE PRECISION TWOCPH, COSPHI(4), SINPHI(4), TEMP
 
      COSTH  = Z
      COSTH2 = COSTH  * COSTH
      COSTH4 = COSTH2 * COSTH2
      SINTH2 = 1. - COSTH2
      SINTH  = DSQRT(SINTH2)
      SCTH   = SINTH  * COSTH
 
 
C *** IS THETA = 0 ? THEN PHI IS IRRELEVANT *****************
 
      IF (SINTH .EQ. 0.) THEN
         COSPHI(1) = 1.
         SINPHI(1) = 1.
      ELSE
         COSPHI(1) = X/SINTH
         SINPHI(1) = Y/SINTH
      ENDIF
      TWOCPH    = 2.*COSPHI(1)
      COSPHI(2) = TWOCPH*COSPHI(1)-1.
      SINPHI(2) = TWOCPH*SINPHI(1)
      COSPHI(3) = TWOCPH*COSPHI(2)-COSPHI(1)
      SINPHI(3) = TWOCPH*SINPHI(2)-SINPHI(1)
      COSPHI(4) = TWOCPH*COSPHI(3)-COSPHI(2)
      SINPHI(4) = TWOCPH*SINPHI(3)-SINPHI(2)
 
C *** NOW THE SPHERICAL HARMONICS ARE CALCULATED AND SUMMED WITH Q
C
C     THIS PART OF THE SUBROUTINE WOULD HAVE BEEN EASIER TO UNDERSTAND
C     WHEN COMPLEX NUMBERS WOULD HAVE BEEN USED. HOWEVER, COMPLEX*16
C     VARIABLES ARE NOT INCLUDED IN THE FORTRAN 77 STANDARD.
 
      Q0    = Q0 + Q4COEF(1,0)*COSTH4+Q4COEF(2,0)*COSTH2
     :             +Q4COEF(3,0)
 
      TEMP  = (Q4COEF(1,1) * COSTH2 + Q4COEF(2,1)) * SCTH
      QR(1) = QR(1) + COSPHI(1)*TEMP
      QI(1) = QI(1) + SINPHI(1)*TEMP
 
      TEMP  = Q4COEF(1,2)*COSTH4+Q4COEF(2,2)*COSTH2+Q4COEF(3,2)
      QR(2) = QR(2) + COSPHI(2)*TEMP
      QI(2) = QI(2) + SINPHI(2)*TEMP
 
      TEMP  = Q4COEF(1,3) * SCTH * SINTH2
      QR(3) = QR(3) + COSPHI(3)*TEMP
      QI(3) = QI(3) + SINPHI(3)*TEMP
 
      TEMP  = Q4COEF(1,4) * SINTH2 * SINTH2
      QR(4) = QR(4) + COSPHI(4)*TEMP
      QI(4) = QI(4) + SINPHI(4)*TEMP
 
      RETURN
      END
 
 
      DOUBLE PRECISION FUNCTION Q4(NB, Q0, QR, QI)
 
C *** THE ORDER PARAMETER Q4 IS CALCULATED FROM Q0 AND Q *************
C
C     INTEGER  NB      NUMBER OF BONDS ACCUMULATED
C     DOUBLE PRECISION
C              Q0      ACCUMULATOR OF P40(COS(THETA))
C              QR(4)   ACCUMULATOR OF P4M(COS(THETA))*COS(PHI)
C              QI(4)   ACCUMULATOR OF P4M(COS(THETA))*SIN(PHI)
C
C     IN THE EXPRESSION FOR Q4, THE FACTOR PI*9/4 DISAPPEARS
C     THE SUMMATION HAS TO BE PERFORMED OVER THE SQUARES OF Q4M, WITH
C     M RUNNING FROM -4 TO +4. HOWEVER, Q4M AND Q4-M ARE EQUAL WHEN
C     SQUARED. SO THE SUMMATION IS DONE FOR POSITIVE M AND A FACTOR
C     2 IS INTRODUCED.
C
C ********************************************************************
 
      INTEGER   NB
      DOUBLE PRECISION  Q0, QR(4), QI(4)
 
      Q4 = DSQRT(    Q0 * Q0
     :          + ( QR(1)*QR(1) + QI(1)*QI(1) )
     :          + ( QR(2)*QR(2) + QI(2)*QI(2) )
     :          + ( QR(3)*QR(3) + QI(3)*QI(3) )
     :          + ( QR(4)*QR(4) + QI(4)*QI(4) ) ) / NB
 
      RETURN
      END
 
 
      DOUBLE PRECISION FUNCTION FAC(N)
 
      INTEGER I, N
 
      FAC = 1.
      IF(N.LT.2)RETURN
      DO 100 I = 2, N
         FAC = FAC*I
100   CONTINUE
      RETURN
      END
