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
C
C  ADD ENERGY AND GRADIENT CORRECTION TERMS FOR THE CONTAINER SEPARATELY.
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

      IF (BSPT) RETURN ! CONTAINER IS ACCOUNTED FOR IN BSPT BY RECOUNTING PREVIOUS CONFIGURATION
      IF (PERIODIC) RETURN
      EVAP=.FALSE.
      EVAPREJECT=.FALSE.
      DO J1=1,NATOMS
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
         IF (DIST.GT.RADIUS) THEN
!           WRITE(MYUNIT,'(A,I5,5G20.10)') 'J1,DIST,RADIUS IN RAD = ',J1,DIST,RADIUS,X(J3-2),X(J3-1),X(J3)
            EVAP=.TRUE.
           ! IF (EVAP.AND.(BSWL.OR.BSPT)) THEN
           !     EVAPREJECT=.TRUE.
           !     IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') 'EVAP: ATOM, RADIUS=',J1,SQRT(DIST)
           !     RETURN
           ! ENDIF

            IF (DEBUG)  WRITE(MYUNIT,'(A,2G20.10,L10)') 'RAD> EVAP: ATOM, RADIUS, EVAP=',J1,SQRT(DIST),EVAP
C           PRINT*,'EVAP: ATOM, RADIUS=',J1,SQRT(DIST)
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
!            WRITE(MYUNIT,'(A,3G20.10)') 'RAD> RESET COORDS: ',X(J3-2),X(J3-1),X(J3)
C
C  PUT IT BACK IN AT THE OPPOSITE END OF A DIAMETER
C
C           X(J3-2)=-X(J3-2)*0.8D0
C           X(J3-1)=-X(J3-1)*0.8D0
C           X(J3)=-X(J3)*0.8D0
         ENDIF
      ENDDO

      RETURN
      END
C
C  FOR RIGID-BODY ANGLE-AXIS COORDINATES, JUST MOVE THE FIXED SITE
C
      SUBROUTINE RADR(X,V,ENERGY,GTEST)
      USE COMMONS
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J3
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), 
     1                 ENERGY, DUMMYX, DUMMYY, DUMMYZ
      LOGICAL EVAP, EVAPREJECT
      COMMON /EV/ EVAP, EVAPREJECT

      IF (PERIODIC) RETURN
      EVAP=.FALSE.
      DO J1=1,NATOMS/2
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C        WRITE(*,'(A,I6,5F15.5)') 'J1,DIST,COORDS,RADIUS IN RADR=',J1,DIST,X(J3-2),X(J3-1),X(J3),RADIUS
         IF (DIST.GT.RADIUS) THEN
            EVAP=.TRUE.
            WRITE(MYUNIT,'(A,I5,5G17.8)') 'EVAP: MOLECULE, COORDS, DIST, RADIUS=',J1,X(J3-2),X(J3-1),X(J3),SQRT(DIST),SQRT(RADIUS)
C           IF (DEBUG) WRITE(*,'(A,I5,2G20.10)') 'EVAP: MOLECULE, DIST, RADIUS=',J1,SQRT(DIST),SQRT(RADIUS)
            DIST=SQRT(RADIUS*0.9D0/DIST)
C           X(J3-2)=X(J3-2)*DIST
C           X(J3-1)=X(J3-1)*DIST
C           X(J3)=X(J3)*DIST
C
C  PUT IT BACK IN AT THE OPPOSITE END OF A DIAMETER
C
            X(J3-2)=-X(J3-2)*0.8D0
            X(J3-1)=-X(J3-1)*0.8D0
            X(J3)=-X(J3)*0.8D0
         ENDIF
      ENDDO

      RETURN
      END
