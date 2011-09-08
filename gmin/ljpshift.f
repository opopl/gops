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
C*************************************************************************
C
C  Subroutine LJPSHIFT calculates the energy, cartesian gradient and second
C  derivative matrix analytically for Lennard-Jones in reduced units
C  (epsilon=sigma=1) using a shifted, truncated potential.
C
C  Adapted for the binary LJ glass described by Sastry, Debenetti and
C  Stillinger, Nature, 393, 554, 1998. Atom types are A and B. The first
C  NTYPEA are A, the next NBTYPE=NATOMS-NTYPEA are B. epsilon and sigma for A are the
C  units of energy and distance, so we also need EPSAB, EPSAB, SIGAB and
C  SIGAA in these units. Sastry et al. density is 1.2 i.e. a box length
C  of 5.975206 for 256 atoms. 
C
C
C*************************************************************************
C
      SUBROUTINE LJPSHIFT(X, V, POTEL, GTEST, STEST)
      USE commons
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA
      DOUBLE PRECISION X(3*NATOMS), A(3*NATOMS,3*NATOMS), VEC1, VEC2, VEC3,
     1                 V(3*NATOMS), R2(NATOMS,NATOMS),  R6, R2DUM,
     2                 R8(NATOMS,NATOMS), G(NATOMS,NATOMS), EPSAB, EPSBB, SIGAB, SIGBB,
     3                 R14(NATOMS,NATOMS), F(NATOMS,NATOMS), 
     4                 POTEL, SIGAB6, SIGAB12, SIGBB6, SIGBB12,
     5                 XVEC(NATOMS,NATOMS,3), IRCUT2AA, IRCUT2AB, IRCUT2BB, SIGRCAA6, SIGRCAA12,
     6                 SIGRCAB6, SIGRCAB12, SIGRCBB6, SIGRCBB12, CONSTAA, CONSTBB, CONSTAB,
     7                 RCONSTAA, RCONSTAB, RCONSTBB, CUTAA, CUTAB, CUTBB
      LOGICAL GTEST, STEST
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /RCONST/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB

      N=NATOMS
      CUTAA=CUTOFF
      CUTAB=CUTOFF*SIGAB
      CUTBB=CUTOFF*SIGBB
      IRCUT2AA = 1.D0/CUTAA**2
      IRCUT2AB = 1.D0/CUTAB**2
      IRCUT2BB = 1.D0/CUTBB**2
      SIGAB6=SIGAB**6
      SIGAB12=SIGAB6**2
      SIGBB6=SIGBB**6
      SIGBB12=SIGBB6**2
      SIGRCAA6= 1.0D0/CUTAA**6
      SIGRCAA12=SIGRCAA6**2
      SIGRCAB6=SIGAB6/CUTAB**6
      SIGRCAB12=SIGRCAB6**2
      SIGRCBB6=SIGBB6/CUTBB**6
      SIGRCBB12=SIGRCBB6**2
      CONSTAA=4.0D0*SIGRCAA6-7.0D0*SIGRCAA12
      CONSTAB=4.0D0*SIGRCAB6-7.0D0*SIGRCAB12
      CONSTBB=4.0D0*SIGRCBB6-7.0D0*SIGRCBB12
      RCONSTAA=(6.0D0*SIGRCAA12-3.0D0*SIGRCAA6)/CUTAA**2
      RCONSTAB=(6.0D0*SIGRCAB12-3.0D0*SIGRCAB6)/CUTAB**2
      RCONSTBB=(6.0D0*SIGRCBB12-3.0D0*SIGRCBB6)/CUTBB**2
C
C  Work out cutoff for potential. Two particles interact if r<c, but
C  we will use the equivalent condition 1/r^2 > 1/c^2.
C
C  Deal with any atoms that have left the box.
C
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF
C
C  Calculate interatomic vectors using the minimum image convention.
C  VEC(i,j,alpha) is the alpha (x, y or z) component of the vector between
C  atoms i and j.
C
      POTEL=0.0D0
      IF (GTEST) THEN
         DO J1=1, N
            XVEC(J1,J1,1)=0.0D0
            XVEC(J1,J1,2)=0.0D0
            XVEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               XVEC(J2,J1,1)=X(J3+1)-X(J4+1)
               XVEC(J2,J1,2)=X(J3+2)-X(J4+2)
               XVEC(J2,J1,3)=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(XVEC(J2,J1,1)/BOXLX)
                  ANV(J2,J1,2)=NINT(XVEC(J2,J1,2)/BOXLY)
                  ANV(J2,J1,3)=NINT(XVEC(J2,J1,3)/BOXLZ)
               ENDIF
               XVEC(J2,J1,1)=XVEC(J2,J1,1)-BOXLX*ANV(J2,J1,1)
               XVEC(J2,J1,2)=XVEC(J2,J1,2)-BOXLY*ANV(J2,J1,2)
               XVEC(J2,J1,3)=XVEC(J2,J1,3)-BOXLZ*ANV(J2,J1,3)
               XVEC(J1,J2,1)=-XVEC(J2,J1,1)
               XVEC(J1,J2,2)=-XVEC(J2,J1,2)
               XVEC(J1,J2,3)=-XVEC(J2,J1,3)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,NTYPEA
               R2DUM=1.0D0/(XVEC(J1,J2,1)**2+XVEC(J1,J2,2)**2+XVEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2AA) THEN
                  POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
               ENDIF
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               R2DUM=1.0D0/(XVEC(J1,J2,1)**2+XVEC(J1,J2,2)**2+XVEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2AB) THEN
                  POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB) ! AB
               ENDIF
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2DUM=1.0D0/(XVEC(J1,J2,1)**2+XVEC(J1,J2,2)**2+XVEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2BB) THEN
                  POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
               ENDIF
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NTYPEA
            J3=3*(J1-1)
            DO J2=J1+1,NTYPEA
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC1/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC2/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
               ENDIF
               VEC1=VEC1-BOXLX*ANV(J2,J1,1)
               VEC2=VEC2-BOXLY*ANV(J2,J1,2)
               VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2AA) THEN
                  POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
               ENDIF
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            J3=3*(J1-1)
            DO J2=NTYPEA+1, N
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC1/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC2/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
               ENDIF
               VEC1=VEC1-BOXLX*ANV(J2,J1,1)
               VEC2=VEC2-BOXLY*ANV(J2,J1,2)
               VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2AB) THEN
                  POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB)  ! AB
               ENDIF
            ENDDO
         ENDDO
         DO J1=NTYPEA+1, N
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC1/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC2/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
               ENDIF
               VEC1=VEC1-BOXLX*ANV(J2,J1,1)
               VEC2=VEC2-BOXLY*ANV(J2,J1,2)
               VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2BB) THEN
                  POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      POTEL = POTEL * 4.D0
     
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
      CALL LJPSHIFTG(G,R2,R14,R8,V,XVEC,N,NATOMS,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)
      
      IF (.NOT.STEST) RETURN
      CALL LJPSHIFTS(G,F,R2,R14,R8,A,XVEC,N,NATOMS,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJPSHIFTG(G,R2,R14,R8,V,XVEC,N,NATOMS,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)
      IMPLICIT NONE
      INTEGER N, NATOMS, J1, J2, J3, J4, NTYPEA
      DOUBLE PRECISION G(NATOMS,NATOMS), R14(NATOMS,NATOMS), R8(NATOMS,NATOMS),
     1                 XVEC(NATOMS,NATOMS,3), V(3*NATOMS), R2(NATOMS,NATOMS),
     2                 EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,
     3                 RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
      COMMON /RCONST/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
C
C  Calculate the g tensor.
C
      DO J1=1,NTYPEA
         G(J1,J1)=0.0D0
         DO J2=J1+1,NTYPEA 
            G(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2AA) G(J2,J1)=-8.0D0      *(3.0D0*(2.0D0*R14(J2,J1)        -R8(J2,J1)       )-RCONSTAA) 
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      DO J1=1,NTYPEA
         DO J2=NTYPEA+1,N 
            G(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2AB) G(J2,J1)=-8.0D0*EPSAB*(3.0D0*(2.0D0*R14(J2,J1)*SIGAB12-R8(J2,J1)*SIGAB6)-RCONSTAB)
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      DO J1=NTYPEA+1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N 
            G(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2BB) G(J2,J1)=-8.0D0*EPSBB*(3.0D0*(2.0D0*R14(J2,J1)*SIGBB12-R8(J2,J1)*SIGBB6)-RCONSTBB)
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO J4=1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) V(J3)=V(J3)+G(J4,J1)*XVEC(J4,J1,J2)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) V(J3)=V(J3)+G(J4,J1)*XVEC(J4,J1,J2)
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) V(J3)=V(J3)+G(J4,J1)*XVEC(J4,J1,J2)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) V(J3)=V(J3)+G(J4,J1)*XVEC(J4,J1,J2)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPSHIFTS(G,F,R2,R14,R8,A,XVEC,N,NATOMS,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)
      IMPLICIT NONE
      INTEGER N, NATOMS, J1, J2, J3, J4, J5, J6, NTYPEA
      DOUBLE PRECISION G(NATOMS,NATOMS), R14(NATOMS,NATOMS), R8(NATOMS,NATOMS), IRCUT2AA, IRCUT2AB, IRCUT2BB,
     1                 XVEC(NATOMS,NATOMS,3), R2(NATOMS,NATOMS), RCONSTAA, RCONSTAB, RCONSTBB,
     2                 F(NATOMS,NATOMS), A(3*NATOMS,3*NATOMS), EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12
      COMMON /RCONST/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=0.0D0
            IF (J1.LE.NTYPEA) THEN
               IF (J2.LE.NTYPEA) THEN
                  IF (R2(J2,J1).GT.IRCUT2AA) F(J2,J1)=8.0D0*      (84.0D0*R14(J2,J1)-24.0D0*R8(J2,J1))  
               ELSE
                  IF (R2(J2,J1).GT.IRCUT2AB) F(J2,J1)=8.0D0*EPSAB*(84.0D0*R14(J2,J1)*SIGAB12-24.0D0*R8(J2,J1)*SIGAB6)
               ENDIF
            ELSE
               IF (R2(J2,J1).GT.IRCUT2BB) F(J2,J1)=   8.0D0*EPSBB*(84.0D0*R14(J2,J1)*SIGBB12-24.0D0*R8(J2,J1)*SIGBB6)
            ENDIF
C           F(J2,J1)=672.0D0*R14(J2,J1)-192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            A(J3,J3)=0.0D0
            DO J4=1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) A(J3,J3)=A(J3,J3)+F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)**2 + G(J4,J1)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) A(J3,J3)=A(J3,J3)+F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)**2 + G(J4,J1)
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) A(J3,J3)=A(J3,J3)+F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)**2 + G(J4,J1)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) A(J3,J3)=A(J3,J3)+F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)**2 + G(J4,J1)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               A(3*(J1-1)+J4,J3)=0.0D0
               DO J5=1,N
                  IF (J1.LE.NTYPEA) THEN
                     IF (J5.LE.NTYPEA) THEN
                        IF (R2(J5,J1).GT.IRCUT2AA) 
     1                      A(3*(J1-1)+J4,J3)=A(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*XVEC(J5,J1,J2)*XVEC(J5,J1,J4)
                     ELSE
                        IF (R2(J5,J1).GT.IRCUT2AB) 
     1                      A(3*(J1-1)+J4,J3)=A(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*XVEC(J5,J1,J2)*XVEC(J5,J1,J4)
                     ENDIF
                  ELSE
                     IF (J5.LE.NTYPEA) THEN
                        IF (R2(J5,J1).GT.IRCUT2AB)
     1                      A(3*(J1-1)+J4,J3)=A(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*XVEC(J5,J1,J2)*XVEC(J5,J1,J4)
                     ELSE
                        IF (R2(J5,J1).GT.IRCUT2BB) 
     1                      A(3*(J1-1)+J4,J3)=A(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*XVEC(J5,J1,J2)*XVEC(J5,J1,J4)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) THEN
                        A(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)**2-G(J4,J1)
                     ELSE
                        A(3*(J4-1)+J2,J3)=0.0D0
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        A(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)**2-G(J4,J1)
                     ELSE
                        A(3*(J4-1)+J2,J3)=0.0D0
                     ENDIF
                  ENDIF
               ELSE
                  IF (R2(J4,J1).GT.IRCUT2BB) THEN
                     A(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)**2-G(J4,J1)
                  ELSE
                     A(3*(J4-1)+J2,J3)=0.0D0
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  IF (J1.LE.NTYPEA) THEN
                     IF (J4.LE.NTYPEA) THEN
                        IF (R2(J4,J1).GT.IRCUT2AA) THEN
                           A(J6,J3)=-F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)*XVEC(J4,J1,J5)
                        ELSE
                           A(J6,J3)=0.0D0
                        ENDIF
                     ELSE
                        IF (R2(J4,J1).GT.IRCUT2AB) THEN
                           A(J6,J3)=-F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)*XVEC(J4,J1,J5)
                        ELSE
                           A(J6,J3)=0.0D0
                        ENDIF
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) THEN
                        A(J6,J3)=-F(J4,J1)*R2(J4,J1)*XVEC(J4,J1,J2)*XVEC(J4,J1,J5)
                     ELSE
                        A(J6,J3)=0.0D0
                     ENDIF
                  ENDIF
                  A(3*(J4-1)+J2,3*(J1-1)+J5)=A(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            A(J1,J2)=A(J2,J1)
         ENDDO
      ENDDO

      RETURN
      END

