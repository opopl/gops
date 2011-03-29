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
C  SUBROUTINE LJPSHIFT CALCULATES THE ENERGY, CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY FOR LENNARD-JONES IN REDUCED UNITS
C  (EPSILON=SIGMA=1) USING A SHIFTED, TRUNCATED POTENTIAL.
C
C  ADAPTED FOR THE BINARY LJ GLASS DESCRIBED BY SASTRY, DEBENETTI AND
C  STILLINGER, NATURE, 393, 554, 1998. ATOM TYPES ARE A AND B. THE FIRST
C  NTYPEA ARE A, THE NEXT NBTYPE=NATOMS-NTYPEA ARE B. EPSILON AND SIGMA FOR A ARE THE
C  UNITS OF ENERGY AND DISTANCE, SO WE ALSO NEED EPSAB, EPSAB, SIGAB AND
C  SIGAA IN THESE UNITS. SASTRY ET AL. DENSITY IS 1.2 I.E. A BOX LENGTH
C  OF 5.975206 FOR 256 ATOMS. 
C
C
C*************************************************************************
C
      SUBROUTINE LJPSHIFT(X, V, POTEL, GTEST, STEST)
      USE COMMONS
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
C  WORK OUT CUTOFF FOR POTENTIAL. TWO PARTICLES INTERACT IF R<C, BUT
C  WE WILL USE THE EQUIVALENT CONDITION 1/R^2 > 1/C^2.
C
C  DEAL WITH ANY ATOMS THAT HAVE LEFT THE BOX.
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
C  CALCULATE INTERATOMIC VECTORS USING THE MINIMUM IMAGE CONVENTION.
C  VEC(I,J,ALPHA) IS THE ALPHA (X, Y OR Z) COMPONENT OF THE VECTOR BETWEEN
C  ATOMS I AND J.
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
C  CALCULATE THE G TENSOR.
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
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
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
C  NEXT ARE THE TERMS WHERE X_I AND X_J ARE ON THE SAME ATOM
C  BUT ARE DIFFERENT, E.G. Y AND Z.
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
C  CASE III, DIFFERENT ATOMS, SAME CARTESIAN COORDINATE.
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
C  CASE IV: DIFFERENT ATOMS AND DIFFERENT CARTESIAN COORDINATES.
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
C  SYMMETRISE HESSIAN
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            A(J1,J2)=A(J2,J1)
         ENDDO
      ENDDO

      RETURN
      END

