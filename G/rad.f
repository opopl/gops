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
C
C  Add energy and gradient correction terms for the container separately.
C
      SUBROUTINE RAD(X,V,ENERGY,GTEST)
      USE COMMONS
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J3
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), 
     1                 ENERGY, DUMMYX, DUMMYY, DUMMYZ
      LOGICAL EVAP, EVAPREJECT
      COMMON /EV/ EVAP, EVAPREJECT

      IF (BSPT) RETURN ! container is accounted for in bspt by recounting previous configuration
      IF (PERIODIC) RETURN
      EVAP=.FALSE.
      EVAPREJECT=.FALSE.
      DO J1=1,NATOMS
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
         IF (DIST.GT.RADIUS) THEN
!           WRITE(LFH,'(A,I5,5G20.10)') 'J1,DIST,RADIUS in rad = ',J1,DIST,RADIUS,X(J3-2),X(J3-1),X(J3)
            EVAP=.TRUE.
           ! IF (EVAP.AND.(BSWL.OR.BSPT)) then
           !     EVAPREJECT=.TRUE.
           !     IF (DEBUG) WRITE(LFH,'(A,2G20.10)') 'EVAP: atom, radius=',J1,SQRT(DIST)
           !     RETURN
           ! ENDIF

            IF (DEBUG)  WRITE(LFH,'(A,2G20.10,L10)') 'rad> EVAP: atom, radius, EVAP=',J1,SQRT(DIST),EVAP
C           PRINT*,'EVAP: atom, radius=',J1,SQRT(DIST)
CC           ENERGY=ENERGY+1.0D5*(DIST-RADIUS)**2
CC           IF (GTEST.AND.(.NOT.(SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE))) THEN
CC              DUMMYX=1.0D5*4.0D0*(DIST-RADIUS)*X(J3-2)
CC              DUMMYY=1.0D5*4.0D0*(DIST-RADIUS)*X(J3-1)
CC              DUMMYZ=1.0D5*4.0D0*(DIST-RADIUS)*X(J3)
CC              V(J3-2)=V(J3-2)+DUMMYX
CC              V(J3-1)=V(J3-1)+DUMMYY
CC              V(J3)=V(J3)+DUMMYZ
CC           ENDIF

             DIST=(SQRT(RADIUS)-0.5D0)/SQRT(DIST)
             X(J3-2)=X(J3-2)*DIST
             X(J3-1)=X(J3-1)*DIST
             X(J3)=X(J3)*DIST
!            WRITE(LFH,'(A,3G20.10)') 'rad> reset coords: ',X(J3-2),X(J3-1),X(J3)
C
C  Put it back in at the opposite end of a diameter
C
C           X(J3-2)=-X(J3-2)*0.8D0
C           X(J3-1)=-X(J3-1)*0.8D0
C           X(J3)=-X(J3)*0.8D0
         ENDIF
      ENDDO

      RETURN
      END
C
C  For rigid-body angle-axis coordinates, just move the fixed site
C
      SUBROUTINE RADR(X,V,ENERGY,GTEST)
      USE commons
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J3
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), 
     1                 ENERGY, DUMMYX, DUMMYY, DUMMYZ
      LOGICAL EVAP, evapreject
      COMMON /EV/ EVAP, evapreject

      IF (PERIODIC) RETURN
      EVAP=.FALSE.
      DO J1=1,NATOMS/2
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C        WRITE(*,'(A,I6,5F15.5)') 'J1,DIST,coords,RADIUS in radr=',J1,DIST,X(J3-2),X(J3-1),X(J3),RADIUS
         IF (DIST.GT.RADIUS) THEN
            EVAP=.TRUE.
            WRITE(LFH,'(A,I5,5G17.8)') 'EVAP: molecule, coords, dist, radius=',J1,X(J3-2),X(J3-1),X(J3),SQRT(DIST),SQRT(RADIUS)
C           IF (DEBUG) WRITE(*,'(A,I5,2G20.10)') 'EVAP: molecule, dist, radius=',J1,SQRT(DIST),SQRT(RADIUS)
            DIST=SQRT(RADIUS*0.9D0/DIST)
C           X(J3-2)=X(J3-2)*DIST
C           X(J3-1)=X(J3-1)*DIST
C           X(J3)=X(J3)*DIST
C
C  Put it back in at the opposite end of a diameter
C
            X(J3-2)=-X(J3-2)*0.8D0
            X(J3-1)=-X(J3-1)*0.8D0
            X(J3)=-X(J3)*0.8D0
         ENDIF
      ENDDO

      RETURN
      END
