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
C  ENERGY AND GRADIENT FOR LJ.
C
      SUBROUTINE LJ(X,V,ELJ,GTEST,SECT)
      USE COMMONS
      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
     1                 R6, ELJ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY
  
      IF (SECT) THEN
         CALL LJDIFF(NATOMS, X)
         RETURN
      ENDIF
      ELJ=0.0D0
      DO J1=1,NATOMS 
         VT(J1)=0.0D0
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            G(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6-1.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               ELJ=ELJ+DUMMY
               DIST=DIST*R6
               G(J2,J1)=-24.0D0*(2.0D0*R6-1.0D0)*DIST
               G(J1,J2)=G(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6-1.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               ELJ=ELJ+DUMMY
            ENDDO
         ENDDO
      ENDIF
      ELJ=ELJ*4.0D0

      IF (.NOT.GTEST) RETURN

      IF (SEEDT.AND.FREEZECORE) THEN
         DO J1=1,NATOMS-NSEED
            J3=3*J1
            DUMMYX=0.0D0
            DUMMYY=0.0D0
            DUMMYZ=0.0D0
            DO J4=1,NATOMS
               J2=3*J4
               XMUL2=G(J4,J1)
               DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
               DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
               DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
            ENDDO
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDDO
         V(3*(NATOMS-NSEED-1)+1:3*NATOMS)=0.0D0
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DUMMYX=0.0D0
            DUMMYY=0.0D0
            DUMMYZ=0.0D0
            DO J4=1,NATOMS
               J2=3*J4
               XMUL2=G(J4,J1)
               DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
               DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
               DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
            ENDDO
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDDO
      ENDIF

      RETURN
      END
C
C*************************************************************************
C
C  SUBROUTINE LJDIFF CALCULATES THE CARTESIAN SECOND
C  DERIVATIVE MATRIX ANALYTICALLY. REDUCED UNITS.
C
C*************************************************************************
C
      SUBROUTINE LJDIFF(N, X)
      IMPLICIT NONE
      INTEGER N, J1, J2,NATOMS
      DOUBLE PRECISION X(3*N), 
     1                 R2(N,N), 
     2                 R8(N,N), G(N,N),
     3                 R14(N,N), F(N,N)
C 
C  STORE DISTANCE MATRICES.
C

      DO J1=1,N
         R2(J1,J1)=0.0D0
         R8(J1,J1)=0.0D0
         R14(J1,J1)=0.0D0
         DO J2=J1+1,N
            R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1               +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2               +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R2(J2,J1)=1.0D0/R2(J2,J1)
            R8(J2,J1)=R2(J2,J1)**4
            R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
            R2(J1,J2)=R2(J2,J1)
         ENDDO
      ENDDO

      CALL LJS(G,F,R2,R14,R8,X,N)

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJS(G,F,R2,R14,R8,X,N)
      USE COMMONS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(NATOMS,NATOMS), R14(NATOMS,NATOMS), R8(NATOMS,NATOMS),
     1                 R2(NATOMS,NATOMS), F(NATOMS,NATOMS), 
     2                 X(3*NATOMS),DUMMY

      DO J1=1,N
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=672.0D0*R14(J2,J1)-192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
            G(J2,J1)=-24.0D0*(2.0D0*R14(J2,J1)-R8(J2,J1))
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
            ENDDO
            HESS(J3,J3)=DUMMY
         ENDDO
      ENDDO
C
C  NEXT ARE THE TERMS WHERE X_I AND X_J ARE ON THE SAME ATOM
C  BUT ARE DIFFERENT, E.G. Y AND Z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               DUMMY=0.0D0
               DO J5=1,N
                  DUMMY=DUMMY + F(J5,J1)*R2(J5,J1)* 
     1           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)) 
               ENDDO
               HESS(3*(J1-1)+J4,J3)=DUMMY
            ENDDO
         ENDDO
      ENDDO
C
C  CASE III, DIFFERENT ATOMS, SAME CARTESIAN COORDINATE.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
            ENDDO
         ENDDO
      ENDDO
C
C  CASE IV: DIFFERENT ATOMS AND DIFFERENT CARTESIAN COORDINATES.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  SYMMETRISE HESSIAN
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  ENERGY AND GRADIENT FOR LJ WITH A CUTOFF
C
      SUBROUTINE LJCUT(X,V,ELJ,GTEST,SECT)
      USE COMMONS
      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
     1                 R6, ELJ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY, CUTOFF2
      LOGICAL EVAP, EVAPREJECT
      COMMON /EV/ EVAP, EVAPREJECT

      CUTOFF2=CUTOFF**2
      IF (SECT) THEN
         CALL LJDIFF(NATOMS, X)
         RETURN
      ENDIF
C     EVAP=.FALSE.
      ELJ=0.0D0
      DO J1=1,NATOMS 
         VT(J1)=0.0D0
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C            IF ((DIST.GT.RADIUS).AND.(.NOT.LJMFT)) THEN
CC              IF (DEBUG) WRITE(*,'(A,I3,A,3G15.5,A,G15.5,A,G15.5)') 
CC    1                    'ATOM ',J1,' AT ',X(J3-2),X(J3-1),X(J3),' DIST=',DIST,' RADIUS=',RADIUS
C               EVAP=.TRUE.
C               ELJ=ELJ+10.0D2*(DIST-RADIUS)**2
C            ENDIF
            G(J1,J1)=0.0D0
            ATOM2: DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.GT.CUTOFF2) CYCLE ATOM2
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6-1.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               ELJ=ELJ+DUMMY
               DIST=DIST*R6
               G(J2,J1)=-24.0D0*(2.0D0*R6-1.0D0)*DIST
               G(J1,J2)=G(J2,J1)
            ENDDO ATOM2
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C            IF ((DIST.GT.RADIUS).AND.(.NOT.LJMFT)) THEN
C               EVAP=.TRUE.
CC              IF (DEBUG) WRITE(*,'(A,I3,A,3G15.5,A,G15.5,A,G15.5)') 
CC    1                    'ATOM ',J1,' AT ',X(J3-2),X(J3-1),X(J3),' DIST=',DIST,' RADIUS=',RADIUS
C               ELJ=ELJ+10.0D2*(DIST-RADIUS)**2
C            ENDIF
            ATOM22: DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.GT.CUTOFF2) CYCLE ATOM22
               DIST=1.0D0/DIST
               R6=DIST**3
               DUMMY=R6*(R6-1.0D0)
               VT(J1)=VT(J1)+DUMMY
               VT(J2)=VT(J2)+DUMMY
               ELJ=ELJ+DUMMY
            ENDDO ATOM22
         ENDDO
      ENDIF
C      IF (DEBUG.AND.EVAP) THEN
CC        WRITE(*,'(A)') 'AN ATOM HAS EVAPORATED - DUMPING COORDINATES'
CC        WRITE(40,'(I4)') NATOMS
CC        WRITE(40,'(A,I4,A,F15.5)') 'ENERGY AFTER EVAP=',ELJ
CC        WRITE(40,'(A2,3F20.10)') ('LJ ',X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3),J1=1,NATOMS)
C      ENDIF
      ELJ=ELJ*4.0D0

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
               XMUL2=G(J4,J1)
               DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
               DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
               DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
            ENDDO
C           DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C           IF ((DIST.GT.RADIUS).AND.(.NOT.LJMFT)) THEN
C              DUMMYX=DUMMYX+10.0D2*16.0D0*(DIST-RADIUS)*X(J3-2)
C              DUMMYY=DUMMYY+10.0D2*16.0D0*(DIST-RADIUS)*X(J3-1)
C              DUMMYZ=DUMMYZ+10.0D2*16.0D0*(DIST-RADIUS)*X(J3)
C           ENDIF
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDIF
      ENDDO

      RETURN
      END
