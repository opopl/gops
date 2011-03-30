C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE FDM(P,GRAD,ENERGY,GRADT)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION P(3*NATOMS), ND5H, NTD, NOH, NIH, PI, ENERGY, 
     1                 GRAD(3*NATOMS), DUMMY1,
     2                 DUMMY2, X, Y, Z, R, DUMMY3, X2, Y2, X4, Y4, 
     3                 DUMMY4, Z2, Z4,
     4                 FAC, XMASS, YMASS, ZMASS, GRADX, GRADY, GRADZ
      INTEGER J1, J2, I, J4, J3
      LOGICAL GRADT
      PARAMETER (PI=3.141592654D0)
C     PARAMETER (ND5H=16.0D0*DSQRT(1.0D0/7.0D0)/(3.0D0*PI), 
C    1           NTD=8.0D0/PI,
C    2           NOH=64.0D0*DSQRT(2.0D0/4449.0D0)/PI,
C    3           NIH=256.0D0*DSQRT(2.0D0/105614964576169.0D0)/PI)
      PARAMETER (ND5H=2.015810523D0/PI,
     1           NTD=8.0D0/PI,
     2           NOH=1.356949761D0/PI,
     3           NIH=0.011140181D0/PI)

      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      DO 10 I=1,NATOMS
         XMASS=XMASS+P(3*(I-1)+1)
         YMASS=YMASS+P(3*(I-1)+2)
         ZMASS=ZMASS+P(3*(I-1)+3)
10    CONTINUE
      XMASS=XMASS/NATOMS
      YMASS=YMASS/NATOMS
      ZMASS=ZMASS/NATOMS

      IF (D5HT) THEN
         DUMMY1=0.0D0
         DUMMY2=FD5H*ND5H
         DO J1=1,NATOMS
            J2=3*J1
            X=P(J2-2)-XMASS
            Y=P(J2-1)-YMASS
            Z=P(J2)-ZMASS
            R=MAX(DSQRT(X*X+Y*Y+Z*Z),1.0D-10)
            FAC=DEXP(-EXPFAC*(R-EXPD)**2)*DSQRT(EXPFAC/PI)
            X=X/R
            Y=Y/R
            Z=Z/R
            X2=X*X
            Y2=Y*Y
            X4=X2*X2
            Y4=Y2*Y2
            DUMMY4=DUMMY2/R
            DUMMY3=X*(X4 - 10.0D0*X2*Y2 + 5.0D0*Y4)
            DUMMY1=DUMMY1 + DUMMY3*FAC
            IF (GRADT) THEN
               GRADX=FAC*(5.0D0*((X4-6.0D0*X2*Y2+Y4)-DUMMY3*X)
     1         -2.0D0*DUMMY3*EXPFAC*(R-EXPD)*X*R)*DUMMY4
               GRADY=FAC*(5.0D0*Y*(4.0D0*X*(-X2+Y2)-DUMMY3)
     1         -2.0D0*DUMMY3*EXPFAC*(R-EXPD)*Y*R)*DUMMY4
               GRADZ=FAC*(-5.0D0*DUMMY3*Z-2.0D0*DUMMY3*EXPFAC*(R-EXPD)
     1                  *Z*R)*DUMMY4 
               DO J3=1,NATOMS
                  J4=3*J3
                  GRAD(J4-2)=GRAD(J4-2)-GRADX/NATOMS
                  GRAD(J4-1)=GRAD(J4-1)-GRADY/NATOMS
                  GRAD(J4)=GRAD(J4)    -GRADZ/NATOMS
               ENDDO
               GRAD(J2-2)=GRAD(J2-2)+GRADX
               GRAD(J2-1)=GRAD(J2-1)+GRADY
               GRAD(J2)=  GRAD(J2)  +GRADZ
            ENDIF
         ENDDO
         ENERGY=ENERGY+DUMMY1*DUMMY2
      ENDIF
      IF (TDT) THEN
         DUMMY1=0.0D0
         DUMMY2=FTD*NTD
         DO J1=1,NATOMS
            J2=3*J1
            X=P(J2-2)-XMASS
            Y=P(J2-1)-YMASS
            Z=P(J2)-ZMASS
            R=MAX(DSQRT(X*X+Y*Y+Z*Z),1.0D-10)
            FAC=DEXP(-EXPFAC*(R-EXPD)**2)*DSQRT(EXPFAC/PI)
            X=X/R
            Y=Y/R
            Z=Z/R
            DUMMY4=DUMMY2/R
            DUMMY3=X*Y*Z
            DUMMY1=DUMMY1 + DUMMY3*FAC
            IF (GRADT) THEN
               GRADX=FAC*(Y*Z-3.0D0*DUMMY3*X-2.0D0*DUMMY3*EXPFAC
     1                  *(R-EXPD)*X*R)*DUMMY4
               GRADY=FAC*(X*Z-3.0D0*DUMMY3*Y-2.0D0*DUMMY3*EXPFAC
     1                  *(R-EXPD)*Y*R)*DUMMY4
               GRADZ=FAC*(X*Y-3.0D0*DUMMY3*Z-2.0D0*DUMMY3*EXPFAC
     1                  *(R-EXPD)*Z*R)*DUMMY4
               DO J3=1,NATOMS
                  J4=3*J3
                  GRAD(J4-2)=GRAD(J4-2)-GRADX/NATOMS
                  GRAD(J4-1)=GRAD(J4-1)-GRADY/NATOMS
                  GRAD(J4)=GRAD(J4)    -GRADZ/NATOMS
               ENDDO
               GRAD(J2-2)=GRAD(J2-2)+GRADX
               GRAD(J2-1)=GRAD(J2-1)+GRADY
               GRAD(J2)=  GRAD(J2)  +GRADZ
            ENDIF
         ENDDO
         ENERGY=ENERGY+DUMMY1*DUMMY2
      ENDIF
      IF (OHT) THEN
         DUMMY1=0.0D0
         DUMMY2=FOH*NOH
         DO J1=1,NATOMS
            J2=3*J1
            X=P(J2-2)-XMASS
            Y=P(J2-1)-YMASS
            Z=P(J2)-ZMASS
            R=MAX(DSQRT(X*X+Y*Y+Z*Z),1.0D-10)
            FAC=DEXP(-EXPFAC*(R-EXPD)**2)*DSQRT(EXPFAC/PI)
            X=X/R
            Y=Y/R
            Z=Z/R
            X2=X*X
            Y2=Y*Y
            Z2=Z*Z
            X4=X2*X2
            Y4=Y2*Y2
            Z4=Z2*Z2
            DUMMY4=DUMMY2/R
            DUMMY3=X4 + Y4 + Z4 - 3.0D0*(X2*Y2 + X2*Z2 + Y2*Z2)
            DUMMY1=DUMMY1 + DUMMY3*FAC
            IF (GRADT) THEN
               GRADX=FAC*(4.0D0*X*(X2-1.5D0*(Y2+Z2)-DUMMY3)
     1                  -2.0D0*DUMMY3*EXPFAC*(R-EXPD)*X*R)*DUMMY4
               GRADY=FAC*(4.0D0*Y*(Y2-1.5D0*(X2+Z2)-DUMMY3)
     1                  -2.0D0*DUMMY3*EXPFAC*(R-EXPD)*Y*R)*DUMMY4
               GRADZ=FAC*(4.0D0*Z*(Z2-1.5D0*(X2+Y2)-DUMMY3)
     1                  -2.0D0*DUMMY3*EXPFAC*(R-EXPD)*Z*R)*DUMMY4
               DO J3=1,NATOMS
                  J4=3*J3
                  GRAD(J4-2)=GRAD(J4-2)-GRADX/NATOMS
                  GRAD(J4-1)=GRAD(J4-1)-GRADY/NATOMS
                  GRAD(J4)=GRAD(J4)    -GRADZ/NATOMS
               ENDDO
               GRAD(J2-2)=GRAD(J2-2)+GRADX
               GRAD(J2-1)=GRAD(J2-1)+GRADY
               GRAD(J2)=  GRAD(J2)  +GRADZ
            ENDIF
         ENDDO
         ENERGY=ENERGY+DUMMY1*DUMMY2
      ENDIF
      IF (IHT) THEN
         DUMMY1=0.0D0
         DUMMY2=FIH*NIH
         DO J1=1,NATOMS
            J2=3*J1
            X=P(J2-2)-XMASS
            Y=P(J2-1)-YMASS
            Z=P(J2)-ZMASS
            R=MAX(DSQRT(X*X+Y*Y+Z*Z),1.0D-10)
            FAC=DEXP(-EXPFAC*(R-EXPD)**2)*DSQRT(EXPFAC/PI)
            X=X/R
            Y=Y/R
            Z=Z/R
            X2=X*X
            Y2=Y*Y
            Z2=Z*Z
            X4=X2*X2
            Y4=Y2*Y2
            Z4=Z2*Z2
            DUMMY4=DUMMY2/R
            DUMMY3=742500.0D0*X2*Y2*Z2 - 25574.0D0*(X4*Y2 + Y4*Z2 
     1                  + X2*Z4) +
     1                     131824.0D0*(X2*Y4 + X4*Z2 + Y2*Z4) 
     1                  + 8250.0D0*(X2*X4 + Y2*Y4 + Z2*Z4)
            DUMMY1=DUMMY1 + DUMMY3*FAC
            IF (GRADT) THEN
               GRADX=FAC*(2.0D0*X*(24750.0D0*X4 + 742500.0D0*Y2*Z2 
     1                  + 131824.0D0*(Y4 + 2*X2*Z2) -
     1                               25574.0D0*(2*X2*Y2 + Z4) 
     1                  - 3.0D0*DUMMY3)
     2                               -2.0D0*DUMMY3*EXPFAC*(R-EXPD)*X*R)
     1                  *DUMMY4
               GRADY=FAC*(2.0D0*Y*(24750.D0*Y4 + 742500.D0*X2*Z2 
     1                  - 25574.D0*(X4 + 2*Y2*Z2) +
     1                               131824.D0*(2*X2*Y2 + Z4) 
     1                  - 3.0D0*DUMMY3)
     2                               -2.0D0*DUMMY3*EXPFAC*(R-EXPD)*Y*R)
     1                  *DUMMY4
               GRADZ=FAC*(2.0D0*Z*(742500.0D0*X2*Y2 + 24750.0D0*Z4 
     1                  - 25574.0D0*(Y4 + 2*X2*Z2) +
     1                               131824.0D0*(X4 + 2*Y2*Z2) 
     1                  - 3.0D0*DUMMY3)
     2                               -2.0D0*DUMMY3*EXPFAC*(R-EXPD)*Z*R)
     1                  *DUMMY4
               DO J3=1,NATOMS
                  J4=3*J3
                  GRAD(J4-2)=GRAD(J4-2)-GRADX/NATOMS
                  GRAD(J4-1)=GRAD(J4-1)-GRADY/NATOMS
                  GRAD(J4)=GRAD(J4)    -GRADZ/NATOMS
               ENDDO
               GRAD(J2-2)=GRAD(J2-2)+GRADX
               GRAD(J2-1)=GRAD(J2-1)+GRADY
               GRAD(J2)=  GRAD(J2)  +GRADZ
            ENDIF
         ENDDO
         ENERGY=ENERGY+DUMMY1*DUMMY2
      ENDIF

      RETURN
      END
