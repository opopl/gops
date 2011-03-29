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
C*************************************************************************
C
C  SUBROUTINE MORSE CALCULATES THE CARTESIAN GRADIENT ANALYTICALLY FOR 
C  THE MORSE POTENTIAL.
C
C*************************************************************************
C
      SUBROUTINE MORSE(X,V,EMORSE,GTEST)
      USE COMMONS
      IMPLICIT NONE 
      LOGICAL GTEST
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), DIST, R, DUMMY,
     1                 RR(NATOMS,NATOMS), EMORSE, DUMMYX, 
     2                 DUMMYY, DUMMYZ, XMUL2
!     LOGICAL EVAP, EVAPREJECT
!     COMMON /EV/ EVAP, EVAPREJECT

!     EVAP=.FALSE.
      EMORSE=0.0D0
      DO J1=1,NATOMS
         VT(J1)=0.0D0
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
!           IF (DIST.GT.RADIUS) THEN
!              EVAP=.TRUE.
!              EMORSE=EMORSE+(DIST-RADIUS)**2
!           ENDIF
            RR(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=DSQRT((X(J3-2)-X(J4-2))**2 + (X(J3-1)-X(J4-1))**2 
     1                  + (X(J3)-X(J4))**2)
               R=DEXP(RHO-RHO*DIST)
               RR(J2,J1)=2.0D0*R*(R-1.0D0)/DIST
               RR(J1,J2)=RR(J2,J1)
               DUMMY=R*(R-2.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               EMORSE=EMORSE+DUMMY
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
!           IF (DIST.GT.RADIUS) THEN
!              EVAP=.TRUE.
!              EMORSE=EMORSE+(DIST-RADIUS)**2
!           ENDIF
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=DSQRT((X(J3-2)-X(J4-2))**2 + (X(J3-1)-X(J4-1))**2 
     1                  + (X(J3)-X(J4))**2)
               R=DEXP(RHO-RHO*DIST)
               DUMMY=R*(R-2.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               EMORSE=EMORSE+DUMMY
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.GTEST) RETURN
 
      DO J1=1,NATOMS
         J3=3*J1
         IF (SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE) THEN
            V(J3-2)=0.0D0
            V(J3-1)=0.0D0
            V(J3)=0.0D0
         ELSE
            DUMMYX=0.0D0
            DUMMYY=0.0D0
            DUMMYZ=0.0D0
            DO J4=1,NATOMS
               J2=3*J4
               XMUL2=RR(J4,J1)
               DUMMYX=DUMMYX-(X(J3-2)-X(J2-2))*XMUL2
               DUMMYY=DUMMYY-(X(J3-1)-X(J2-1))*XMUL2
               DUMMYZ=DUMMYZ-(X(J3)-X(J2))*XMUL2
            ENDDO
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
            IF (DIST.GT.RADIUS) THEN
               DUMMYX=DUMMYX+4.0D0*(DIST-RADIUS)*X(J3-2)
               DUMMYY=DUMMYY+4.0D0*(DIST-RADIUS)*X(J3-1)
               DUMMYZ=DUMMYZ+4.0D0*(DIST-RADIUS)*X(J3)
            ENDIF
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDIF
      ENDDO

      RETURN
      END
