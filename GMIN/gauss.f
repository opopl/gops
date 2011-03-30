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
      SUBROUTINE KEGEN
      USE COMMONS, ONLY: GAUSSKK, GAUSSEE, GMODES, NATOMS, GKSMALL
      IMPLICIT NONE
      INTEGER IDUM, J, I
      DOUBLE PRECISION PI, GASDEV, DPRAND, SMALLEST

      IDUM = 1313417
      PI = ATAN(1.0D0)*4.0D0
      DO J=1,GMODES
!        EE(J) = 2*PI*ran1(idum)  
         GAUSSEE(J) = 2.0D0*PI*DPRAND()  
         DO I=1,3*NATOMS
            GAUSSKK(I,J)=GASDEV(IDUM)
         ENDDO
      ENDDO
      DO I=1,3*NATOMS
         SMALLEST=1.0D100
         DO J=1,GMODES
            IF (DABS(GAUSSKK(I,J)).LT.SMALLEST) SMALLEST=GAUSSKK(I,J)
         ENDDO
         PRINT '(A,I8,G20.10)','I,smallest kk=',I,SMALLEST
!
!   GKSMALL(I) sets the length scale for the corresponding variable.
!   The efective range is 2*Pi/GKSMALL(I)
!
         GKSMALL(I)=SMALLEST
      ENDDO

      RETURN
      END

      FUNCTION gasdev(idum)
      IMPLICIT NONE
      INTEGER idum
      DOUBLE PRECISION gasdev, DPRAND
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
! 1       v1=2.*ran1(idum)-1.
!         v2=2.*ran1(idum)-1.
1       v1=2.*DPRAND()-1.
        v2=2.*DPRAND()-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

      SUBROUTINE GFIELD(X,V,FUNCVALUE,GTEST)
      USE COMMONS, ONLY: NATOMS, GMODES, GAUSSKK, GAUSSEE
      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J, JJ, I
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), FUNCVALUE, NORM
      DOUBLE PRECISION SPROD(GMODES), SUMF

      NORM=SQRT(2.0D0)/SQRT(1.0D0*GMODES)

      DO J=1,GMODES
         SPROD(J)=0.0D0
         DO JJ=1,3*NATOMS
            SPROD(J)=SPROD(J)+GAUSSKK(JJ,J)*X(JJ)
         ENDDO
      ENDDO

      SUMF=0.0D0
      DO J=1,GMODES
          SUMF=SUMF+COS(SPROD(J)+GAUSSEE(J))
      ENDDO
 
      FUNCVALUE=NORM*SUMF

      DO JJ=1,3*NATOMS
         SUMF=0.0D0
         DO J=1,GMODES
            SUMF=SUMF-GAUSSKK(JJ,J)*SIN(SPROD(J)+GAUSSEE(J))
         ENDDO
         V(JJ)=SUMF*NORM
      ENDDO
  
      RETURN
       
      END
