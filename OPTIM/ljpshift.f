C   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF OPTIM.
C
C   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
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
C  SUBROUTINE LJPSHIFTBIN CALCULATES THE ENERGY, CARTESIAN GRADIENT AND SECOND
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
C  THE SHIFTING AND TRUNCATION IS AS DESCRIBED BY STODDARD AND FORD, PHYS.
C  REV. A, 8, 1504, 1973. ATOM TYPE 'LS'.
C
C
C*************************************************************************
C
      SUBROUTINE LJPSHIFTBIN(N, X, V, POTEL, BOXLX, BOXLY, BOXLZ, CUTOFF, GTEST, STEST, PTEST, BOXTEST)
      USE KEY
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA, ANV(N,N,3)
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3, VBOX(3), SBOX(3),
     1                 V(3*N), R2(N,N),  R6, R2DUM,
     2                 R8(N,N), G(N,N), EPSAB, EPSBB, SIGAB, SIGBB,
     3                 R14(N,N), F(N,N), 
     4                 POTEL, BOXLX, BOXLY, BOXLZ, SIGAB6, SIGAB12, SIGBB6, SIGBB12,
     5                 VEC(N,N,3), IRCUT2AA, IRCUT2AB, IRCUT2BB, CUTOFF, SIGRCAA6, SIGRCAA12,
     6                 SIGRCAB6, SIGRCAB12, SIGRCBB6, SIGRCBB12, CONSTAA, CONSTBB, CONSTAB,
     7                 RCONSTAA, RCONSTAB, RCONSTBB, CUTAA, CUTAB, CUTBB
      LOGICAL GTEST, STEST, PTEST, BOXTEST
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
      COMMON /BOXDERIVS/ VBOX, SBOX
C     COMMON /ANVCOMMON/ ANV

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
C     IF (PTEST) PRINT*,'CUTOFF USED = ',CUTOFF
C
C  DEAL WITH ANY ATOMS THAT HAVE LEFT THE BOX.
C
C     PRINT*,'IN LJPSHIFT FIXIMAGE,NORESET,GTEST=',FIXIMAGE,NORESET,GTEST
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
C           IF (ANINT(X(J2+1)/BOXLX).NE.0) PRINT*,'RESETTING X FOR J1,X,XNEW=',J1,X(J2+1),X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
C           IF (ANINT(X(J2+2)/BOXLX).NE.0) PRINT*,'RESETTING Y FOR J1,X,XNEW=',J1,X(J2+2),X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
C           IF (ANINT(X(J2+3)/BOXLX).NE.0) PRINT*,'RESETTING Z FOR J1,X,XNEW=',J1,X(J2+3),X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
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
      IF (GTEST.OR.BOXTEST) THEN
         DO J1=1, N
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(J2,J1,1)=X(J3+1)-X(J4+1)
               VEC(J2,J1,2)=X(J3+2)-X(J4+2)
               VEC(J2,J1,3)=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC(J2,J1,1)/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC(J2,J1,2)/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC(J2,J1,3)/BOXLZ)
               ENDIF
               VEC(J2,J1,1)=VEC(J2,J1,1)-BOXLX*ANV(J2,J1,1)
               VEC(J2,J1,2)=VEC(J2,J1,2)-BOXLY*ANV(J2,J1,2)
               VEC(J2,J1,3)=VEC(J2,J1,3)-BOXLZ*ANV(J2,J1,3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,NTYPEA
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2AA) THEN
                  POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
C                 IF (R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA.GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
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
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2AB) THEN
                  POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB) ! AB
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                         EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
C                 IF (EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB).GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                         EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
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
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2BB) THEN
                  POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB),POTEL
C                 IF (EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB).GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB),POTEL
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
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
C                 IF (R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA.GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
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
C              IF (J2.EQ.211) WRITE(*,'(A,I4,4F20.10)') 'J1,VEC12,VEC22,VEC32,R2DUM=',J1,VEC1**2,VEC2**2,VEC3**2,R2DUM
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2AB) THEN
                  POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB)  ! AB
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                       EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
C                 IF (EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB).GT.10.0D0)
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                       EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
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
C              IF (J2.EQ.211) WRITE(*,'(A,I4,4F20.10)') 'J1,VEC12,VEC22,VEC32,R2DUM=',J1,VEC1**2,VEC2**2,VEC3**2,R2DUM
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2BB) THEN
                  POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
C                 IF (EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB).GT.10.0D0)
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                 EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB),POTEL
               ENDIF
            ENDDO
         ENDDO
C        WRITE(*,'(A,6F20.10)') 'X53,X93=,',X(157),X(158),X(159),X(277),X(278),X(279)
      ENDIF
      POTEL = POTEL * 4.D0
     
      IF ((.NOT.GTEST).AND.(.NOT.STEST).AND.(.NOT.BOXTEST)) RETURN
      CALL LJPSHIFTGBIN(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,VBOX,BOXTEST,FRACTIONAL,ANV)
      IF (BOXTEST) THEN
         IF (FRACTIONAL) THEN
            VBOX(1)=VBOX(1)/BOXLX+PRESS*BOXLY*BOXLZ
            VBOX(2)=VBOX(2)/BOXLY+PRESS*BOXLX*BOXLZ
            VBOX(3)=VBOX(3)/BOXLZ+PRESS*BOXLY*BOXLX
         ELSE
            VBOX(1)=VBOX(1)+PRESS*BOXLY*BOXLZ
            VBOX(2)=VBOX(2)+PRESS*BOXLX*BOXLZ
            VBOX(3)=VBOX(3)+PRESS*BOXLY*BOXLX
         ENDIF
C        IF (BOXTEST) WRITE(*,'(A,3G20.10)') 'BOX LENGTH GRADIENT: ',VBOX(1),VBOX(2),VBOX(3)
      ENDIF
      
      IF ((.NOT.STEST).AND.(.NOT.BOXTEST)) RETURN
      CALL LJPSHIFTSBIN(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,SBOX,BOXTEST,FRACTIONAL,ANV)
      IF (BOXTEST) THEN
         IF (FRACTIONAL) THEN
            SBOX(1)=SBOX(1)/BOXLX**2
            SBOX(2)=SBOX(2)/BOXLY**2
            SBOX(3)=SBOX(3)/BOXLZ**2
         ENDIF
C        IF (BOXTEST) WRITE(*,'(A,3G20.10)') 'BOX LENGTH SECOND DERIVATIVES: ',SBOX(1),SBOX(2),SBOX(3)
      ENDIF

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJPSHIFTGBIN(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,VBOX,BOXTEST,
     1                        FRACTIONAL,ANV)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA, ANV(N,N,3)
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), VBOX(3),
     1                 VEC(N,N,3), V(3*N), R2(N,N),
     2                 EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,
     3                 RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
      LOGICAL BOXTEST,FRACTIONAL
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
C     COMMON /ANVCOMMON/ ANV
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
      IF (BOXTEST) GOTO 10
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO J4=1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C  BOX LENGTH DERIVATIVES
C
10    VBOX(1)=0.0D0
      VBOX(2)=0.0D0
      VBOX(3)=0.0D0
      IF (FRACTIONAL) THEN
         DO J1=1,N
            DO J4=J1+1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) THEN
                        VBOX(1)=VBOX(1)+G(J4,J1)*VEC(J4,J1,1)**2
                        VBOX(2)=VBOX(2)+G(J4,J1)*VEC(J4,J1,2)**2
                        VBOX(3)=VBOX(3)+G(J4,J1)*VEC(J4,J1,3)**2
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        VBOX(1)=VBOX(1)+G(J4,J1)*VEC(J4,J1,1)**2
                        VBOX(2)=VBOX(2)+G(J4,J1)*VEC(J4,J1,2)**2
                        VBOX(3)=VBOX(3)+G(J4,J1)*VEC(J4,J1,3)**2
                     ENDIF
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        VBOX(1)=VBOX(1)+G(J4,J1)*VEC(J4,J1,1)**2
                        VBOX(2)=VBOX(2)+G(J4,J1)*VEC(J4,J1,2)**2
                        VBOX(3)=VBOX(3)+G(J4,J1)*VEC(J4,J1,3)**2
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) THEN
                        VBOX(1)=VBOX(1)+G(J4,J1)*VEC(J4,J1,1)**2
                        VBOX(2)=VBOX(2)+G(J4,J1)*VEC(J4,J1,2)**2
                        VBOX(3)=VBOX(3)+G(J4,J1)*VEC(J4,J1,3)**2
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO J1=1,N
            DO J4=J1+1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) THEN
                        VBOX(1)=VBOX(1)-G(J4,J1)*VEC(J4,J1,1)*ANV(J4,J1,1)
                        VBOX(2)=VBOX(2)-G(J4,J1)*VEC(J4,J1,2)*ANV(J4,J1,2)
                        VBOX(3)=VBOX(3)-G(J4,J1)*VEC(J4,J1,3)*ANV(J4,J1,3)
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        VBOX(1)=VBOX(1)-G(J4,J1)*VEC(J4,J1,1)*ANV(J4,J1,1)
                        VBOX(2)=VBOX(2)-G(J4,J1)*VEC(J4,J1,2)*ANV(J4,J1,2)
                        VBOX(3)=VBOX(3)-G(J4,J1)*VEC(J4,J1,3)*ANV(J4,J1,3)
                     ENDIF
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        VBOX(1)=VBOX(1)-G(J4,J1)*VEC(J4,J1,1)*ANV(J4,J1,1)
                        VBOX(2)=VBOX(2)-G(J4,J1)*VEC(J4,J1,2)*ANV(J4,J1,2)
                        VBOX(3)=VBOX(3)-G(J4,J1)*VEC(J4,J1,3)*ANV(J4,J1,3)
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) THEN
                        VBOX(1)=VBOX(1)-G(J4,J1)*VEC(J4,J1,1)*ANV(J4,J1,1)
                        VBOX(2)=VBOX(2)-G(J4,J1)*VEC(J4,J1,2)*ANV(J4,J1,2)
                        VBOX(3)=VBOX(3)-G(J4,J1)*VEC(J4,J1,3)*ANV(J4,J1,3)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPSHIFTSBIN(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,SBOX,BOXTEST,
     1                        FRACTIONAL,ANV)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NTYPEA, ANV(N,N,3)
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), IRCUT2AA, IRCUT2AB, IRCUT2BB, SBOX(3),
     1                 VEC(N,N,3), R2(N,N), RCONSTAA, RCONSTAB, RCONSTBB,
     2                 F(N,N), EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12
      LOGICAL BOXTEST,FRACTIONAL
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
C     COMMON /ANVCOMMON/ ANV

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
      
      IF (BOXTEST) GOTO 10
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO J4=1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
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
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO J5=1,N
                  IF (J1.LE.NTYPEA) THEN
                     IF (J5.LE.NTYPEA) THEN
                        IF (R2(J5,J1).GT.IRCUT2AA) 
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ELSE
                        IF (R2(J5,J1).GT.IRCUT2AB) 
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ENDIF
                  ELSE
                     IF (J5.LE.NTYPEA) THEN
                        IF (R2(J5,J1).GT.IRCUT2AB)
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ELSE
                        IF (R2(J5,J1).GT.IRCUT2BB) 
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
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
                        HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
                     ELSE
                        HESS(3*(J4-1)+J2,J3)=0.0D0
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
                     ELSE
                        HESS(3*(J4-1)+J2,J3)=0.0D0
                     ENDIF
                  ENDIF
               ELSE
                  IF (R2(J4,J1).GT.IRCUT2BB) THEN
                     HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
                  ELSE
                     HESS(3*(J4-1)+J2,J3)=0.0D0
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
                           HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                        ELSE
                           HESS(J6,J3)=0.0D0
                        ENDIF
                     ELSE
                        IF (R2(J4,J1).GT.IRCUT2AB) THEN
                           HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                        ELSE
                           HESS(J6,J3)=0.0D0
                        ENDIF
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) THEN
                        HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                     ELSE
                        HESS(J6,J3)=0.0D0
                     ENDIF
                  ENDIF
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
C
C  BOX LENGTH DERIVATIVES
C
10    SBOX(1)=0.0D0
      SBOX(2)=0.0D0
      SBOX(3)=0.0D0
      IF (FRACTIONAL) THEN
         DO J1=1,N
            DO J4=J1+1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) THEN
                        SBOX(1)=SBOX(1)+(F(J4,J1)*VEC(J4,J1,1)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,1)**2
                        SBOX(2)=SBOX(2)+(F(J4,J1)*VEC(J4,J1,2)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,2)**2
                        SBOX(3)=SBOX(3)+(F(J4,J1)*VEC(J4,J1,3)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,3)**2
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        SBOX(1)=SBOX(1)+(F(J4,J1)*VEC(J4,J1,1)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,1)**2
                        SBOX(2)=SBOX(2)+(F(J4,J1)*VEC(J4,J1,2)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,2)**2
                        SBOX(3)=SBOX(3)+(F(J4,J1)*VEC(J4,J1,3)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,3)**2
                     ENDIF
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        SBOX(1)=SBOX(1)+(F(J4,J1)*VEC(J4,J1,1)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,1)**2
                        SBOX(2)=SBOX(2)+(F(J4,J1)*VEC(J4,J1,2)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,2)**2
                        SBOX(3)=SBOX(3)+(F(J4,J1)*VEC(J4,J1,3)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,3)**2
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) THEN
                        SBOX(1)=SBOX(1)+(F(J4,J1)*VEC(J4,J1,1)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,1)**2
                        SBOX(2)=SBOX(2)+(F(J4,J1)*VEC(J4,J1,2)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,2)**2
                        SBOX(3)=SBOX(3)+(F(J4,J1)*VEC(J4,J1,3)**2*R2(J4,J1)+G(J4,J1))*VEC(J4,J1,3)**2
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO J1=1,N
            DO J4=J1+1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) THEN
                        SBOX(1)=SBOX(1)+(F(J4,J1)*VEC(J4,J1,1)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,1)**2
                        SBOX(2)=SBOX(2)+(F(J4,J1)*VEC(J4,J1,2)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,2)**2
                        SBOX(3)=SBOX(3)+(F(J4,J1)*VEC(J4,J1,3)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,3)**2
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        SBOX(1)=SBOX(1)+(F(J4,J1)*VEC(J4,J1,1)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,1)**2
                        SBOX(2)=SBOX(2)+(F(J4,J1)*VEC(J4,J1,2)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,2)**2
                        SBOX(3)=SBOX(3)+(F(J4,J1)*VEC(J4,J1,3)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,3)**2
                     ENDIF
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        SBOX(1)=SBOX(1)+(F(J4,J1)*VEC(J4,J1,1)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,1)**2
                        SBOX(2)=SBOX(2)+(F(J4,J1)*VEC(J4,J1,2)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,2)**2
                        SBOX(3)=SBOX(3)+(F(J4,J1)*VEC(J4,J1,3)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,3)**2
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) THEN
                        SBOX(1)=SBOX(1)+(F(J4,J1)*VEC(J4,J1,1)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,1)**2
                        SBOX(2)=SBOX(2)+(F(J4,J1)*VEC(J4,J1,2)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,2)**2
                        SBOX(3)=SBOX(3)+(F(J4,J1)*VEC(J4,J1,3)**2*R2(J4,J1)+G(J4,J1))*ANV(J4,J1,3)**2
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END
C
C  SAME AS ABOVE WITHOUT PERIODIC BOUNDARIES.
C
      SUBROUTINE LJPSHIFTBINC(N, X, V, POTEL, CUTOFF, GTEST, STEST, PTEST)
      USE KEY
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3, 
     1                 V(3*N), R2(N,N),  R6, R2DUM,
     2                 R8(N,N), G(N,N), EPSAB, EPSBB, SIGAB, SIGBB,
     3                 R14(N,N), F(N,N), 
     4                 POTEL, SIGAB6, SIGAB12, SIGBB6, SIGBB12,
     5                 VEC(N,N,3), IRCUT2AA, IRCUT2AB, IRCUT2BB, CUTOFF, SIGRCAA6, SIGRCAA12,
     6                 SIGRCAB6, SIGRCAB12, SIGRCBB6, SIGRCBB12, CONSTAA, CONSTBB, CONSTAB,
     7                 RCONSTAA, RCONSTAB, RCONSTBB, CUTAA, CUTAB, CUTBB
      LOGICAL GTEST, STEST, PTEST
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB

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
C     IF (PTEST) PRINT*,'CUTOFF USED = ',CUTOFF
C
C  CALCULATE INTERATOMIC VECTORS 
C  VEC(I,J,ALPHA) IS THE ALPHA (X, Y OR Z) COMPONENT OF THE VECTOR BETWEEN
C  ATOMS I AND J.
C
      POTEL=0.0D0
      IF (GTEST) THEN
         DO J1=1, N
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(J2,J1,1)=X(J3+1)-X(J4+1)
               VEC(J2,J1,2)=X(J3+2)-X(J4+2)
               VEC(J2,J1,3)=X(J3+3)-X(J4+3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,NTYPEA
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2AA) THEN
                  POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
C                 IF (R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA.GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
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
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2AB) THEN
                  POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB) ! AB
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                         EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
C                 IF (EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB).GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                         EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
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
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2BB) THEN
                  POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB),POTEL
C                 IF (EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB).GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB),POTEL
               ENDIF
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
C        WRITE(*,'(A,6F20.10)') 'X53,X93=,',X(157),X(158),X(159),X(277),X(278),X(279)
      ELSE
         DO J1=1,NTYPEA
            J3=3*(J1-1)
            DO J2=J1+1,NTYPEA
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2AA) THEN
                  POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
C                 IF (R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA.GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
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
               R2DUM=VEC1**2+VEC2**2+VEC3**2
C              IF (J2.EQ.211) WRITE(*,'(A,I4,4F20.10)') 'J1,VEC12,VEC22,VEC32,R2DUM=',J1,VEC1**2,VEC2**2,VEC3**2,R2DUM
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2AB) THEN
                  POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB)  ! AB
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                       EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
C                 IF (EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB).GT.10.0D0)
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                       EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
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
               R2DUM=VEC1**2+VEC2**2+VEC3**2
C              IF (J2.EQ.211) WRITE(*,'(A,I4,4F20.10)') 'J1,VEC12,VEC22,VEC32,R2DUM=',J1,VEC1**2,VEC2**2,VEC3**2,R2DUM
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2BB) THEN
                  POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
C                 IF (EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB).GT.10.0D0)
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                 EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB),POTEL
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      POTEL = POTEL * 4.D0
     
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
      CALL LJPSHIFTGBINC(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)
      
      IF (.NOT.STEST) RETURN
      CALL LJPSHIFTSBINC(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJPSHIFTGBINC(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), 
     1                 VEC(N,N,3), V(3*N), R2(N,N),
     2                 EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,
     3                 RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
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
                     IF (R2(J4,J1).GT.IRCUT2AA) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPSHIFTSBINC(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NTYPEA
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), IRCUT2AA, IRCUT2AB, IRCUT2BB, 
     1                 VEC(N,N,3), R2(N,N), RCONSTAA, RCONSTAB, RCONSTBB,
     2                 F(N,N), EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB

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
            HESS(J3,J3)=0.0D0
            DO J4=1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
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
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO J5=1,N
                  IF (J1.LE.NTYPEA) THEN
                     IF (J5.LE.NTYPEA) THEN
                        IF (R2(J5,J1).GT.IRCUT2AA) 
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ELSE
                        IF (R2(J5,J1).GT.IRCUT2AB) 
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ENDIF
                  ELSE
                     IF (J5.LE.NTYPEA) THEN
                        IF (R2(J5,J1).GT.IRCUT2AB)
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ELSE
                        IF (R2(J5,J1).GT.IRCUT2BB) 
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
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
                        HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
                     ELSE
                        HESS(3*(J4-1)+J2,J3)=0.0D0
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
                     ELSE
                        HESS(3*(J4-1)+J2,J3)=0.0D0
                     ENDIF
                  ENDIF
               ELSE
                  IF (R2(J4,J1).GT.IRCUT2BB) THEN
                     HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
                  ELSE
                     HESS(3*(J4-1)+J2,J3)=0.0D0
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
                           HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                        ELSE
                           HESS(J6,J3)=0.0D0
                        ENDIF
                     ELSE
                        IF (R2(J4,J1).GT.IRCUT2AB) THEN
                           HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                        ELSE
                           HESS(J6,J3)=0.0D0
                        ENDIF
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) THEN
                        HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                     ELSE
                        HESS(J6,J3)=0.0D0
                     ENDIF
                  ENDIF
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
      SUBROUTINE LJPSHIFT(N, X, V, POTEL, BOXLX, BOXLY, BOXLZ, CUTOFF, GTEST, STEST, PTEST)
      USE KEY
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, ANV(N,N,3)
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3,
     1                 V(3*N), R2(N,N),  R6, R2DUM,
     2                 R8(N,N), G(N,N), 
     3                 R14(N,N), F(N,N), 
     4                 POTEL, BOXLX, BOXLY, BOXLZ, 
     5                 VEC(N,N,3), IRCUT2AA, CUTOFF, SIGRCAA6, SIGRCAA12,
     6                 CONSTAA, RCONSTAA, CUTAA
      LOGICAL GTEST, STEST, PTEST
      COMMON /RCONST/ RCONSTAA, IRCUT2AA
C     COMMON /ANVCOMMON/ ANV

      CUTAA=CUTOFF
      IRCUT2AA = 1.D0/CUTAA**2
      SIGRCAA6= 1.0D0/CUTAA**6
      SIGRCAA12=SIGRCAA6**2
      CONSTAA=4.0D0*SIGRCAA6-7.0D0*SIGRCAA12
      RCONSTAA=(6.0D0*SIGRCAA12-3.0D0*SIGRCAA6)/CUTAA**2
C
C  WORK OUT CUTOFF FOR POTENTIAL. TWO PARTICLES INTERACT IF R<C, BUT
C  WE WILL USE THE EQUIVALENT CONDITION 1/R^2 > 1/C^2.
C
C     IF (PTEST) PRINT*,'CUTOFF USED = ',CUTOFF
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
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(J2,J1,1)=X(J3+1)-X(J4+1)
               VEC(J2,J1,2)=X(J3+2)-X(J4+2)
               VEC(J2,J1,3)=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC(J2,J1,1)/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC(J2,J1,2)/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC(J2,J1,3)/BOXLZ)
               ENDIF
               VEC(J2,J1,1)=VEC(J2,J1,1)-BOXLX*ANV(J2,J1,1)
               VEC(J2,J1,2)=VEC(J2,J1,2)-BOXLY*ANV(J2,J1,2)
               VEC(J2,J1,3)=VEC(J2,J1,3)-BOXLZ*ANV(J2,J1,3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
            ENDDO
         ENDDO
         DO J1=1,N
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               IF (DOUBLET.AND.(J1.EQ.1).AND.(J2.EQ.2)) THEN
                  R2(J2,J1)=0.0D0 
                  R6=0.0D0
               ELSE
                  R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
                  R2(J2,J1)=R2DUM
                  R6=R2(J2,J1)**3
                  IF (R2DUM.GT.IRCUT2AA) THEN
                     POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA
                  ENDIF
               ENDIF
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,N
            J3=3*(J1-1)
            DO J2=J1+1,N
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
               IF (.NOT.(DOUBLET.AND.(J1.EQ.1).AND.(J2.EQ.2))) THEN
                  IF (R2DUM.GT.IRCUT2AA) THEN
                     POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA 
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      POTEL = POTEL * 4.D0
     
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
      CALL LJPSHIFTG(G,R2,R14,R8,V,VEC,N,DOUBLET)
      
      IF (.NOT.STEST) RETURN
      CALL LJPSHIFTS(G,F,R2,R14,R8,VEC,N)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJPSHIFTG(G,R2,R14,R8,V,VEC,N,DOUBLET)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 VEC(N,N,3), V(3*N), R2(N,N),
     3                 RCONSTAA, IRCUT2AA
      LOGICAL DOUBLET
      COMMON /RCONST/ RCONSTAA, IRCUT2AA
C
C  CALCULATE THE G TENSOR.
C
      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            G(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2AA) G(J2,J1)=-8.0D0      *(3.0D0*(2.0D0*R14(J2,J1)        -R8(J2,J1)       )-RCONSTAA) 
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      IF (DOUBLET) THEN
         G(1,2)=0.0D0
         G(2,1)=0.0D0
      ENDIF
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO J4=1,N
               IF (R2(J4,J1).GT.IRCUT2AA) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPSHIFTS(G,F,R2,R14,R8,VEC,N)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), IRCUT2AA, 
     1                 VEC(N,N,3), R2(N,N), RCONSTAA, 
     2                 F(N,N)
      COMMON /RCONST/ RCONSTAA, IRCUT2AA

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2AA) F(J2,J1)=8.0D0*      (84.0D0*R14(J2,J1)-24.0D0*R8(J2,J1))  
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO J4=1,N
               IF (R2(J4,J1).GT.IRCUT2AA) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
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
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO J5=1,N
                  IF (R2(J5,J1).GT.IRCUT2AA) 
     1          HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
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
               IF (R2(J4,J1).GT.IRCUT2AA) THEN
                  HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
               ELSE
                  HESS(3*(J4-1)+J2,J3)=0.0D0
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
                  IF (R2(J4,J1).GT.IRCUT2AA) THEN
                     HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                  ELSE
                     HESS(J6,J3)=0.0D0
                  ENDIF
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

      SUBROUTINE LJPSLATMIN(N,XC,BOXLX,BOXLY,BOXLZ,CUTOFF,V)
      IMPLICIT NONE
      INTEGER N, J1, NCOUNT
      DOUBLE PRECISION XC(3*N),F1,F2,F3,GRAD,SECOND,
     1                 XSAVE(3*N),BSAVE,V(3*N),
     2                 BOXLX,BOXLY,BOXLZ,CUTOFF,DIF,TEMP1

C
C  VALUE OF DIF IS THE ORDER OF MAGNITUDE TO WHICH THE LATTICE
C  CONSTANT CAN BE OPTIMISED. SETTING IT SMALLER THAN 10^(-7)
C  CAUSES NUMERICAL PROBLEMS ON THE DEC.
C
      DIF=1.0D-7
      BSAVE=BOXLX
      NCOUNT=1
      DO 20 J1=1,3*N
         XSAVE(J1)=XC(J1)
20    CONTINUE

10    TEMP1=BOXLX
      BOXLX=TEMP1+DIF
      BOXLY=TEMP1+DIF
      BOXLZ=TEMP1+DIF

      CALL LJPSHIFT(N, XC,  V, F1,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      BOXLX=TEMP1-DIF
      BOXLY=TEMP1-DIF
      BOXLZ=TEMP1-DIF
      DO 30 J1=1,3*N
         XC(J1)=XSAVE(J1)
30    CONTINUE
      CALL LJPSHIFT(N, XC,  V, F2,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      GRAD=(F1-F2)/(2.0D0*DIF)

      BOXLX=TEMP1
      BOXLY=TEMP1
      BOXLZ=TEMP1
      DO 40 J1=1,3*N
         XC(J1)=XSAVE(J1)
40    CONTINUE
      CALL LJPSHIFT(N, XC,  V, F1,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      BOXLX=TEMP1+2.0D0*DIF
      BOXLY=TEMP1+2.0D0*DIF
      BOXLZ=TEMP1+2.0D0*DIF
      DO 50 J1=1,3*N
         XC(J1)=XSAVE(J1)
50    CONTINUE
      CALL LJPSHIFT(N, XC,  V, F2,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      BOXLX=TEMP1-2.0D0*DIF
      BOXLY=TEMP1-2.0D0*DIF
      BOXLZ=TEMP1-2.0D0*DIF
      DO 60 J1=1,3*N
         XC(J1)=XSAVE(J1)
60    CONTINUE
      CALL LJPSHIFT(N, XC,  V, F3,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      SECOND=(F3+F2-2.0D0*F1)/(4.0D0*DIF*DIF)
      PRINT*,'ENERGY FOR LATTICE CYCLE ',NCOUNT,' IS ',F1
      PRINT*,'GRADIENT WRT BOX LENGTH=',GRAD
      PRINT*,'SECOND DERIVATIVE WRT BOX LENGTH=',SECOND
      NCOUNT=NCOUNT+1
      IF (DABS(GRAD/SECOND).GT.1.0D-4) THEN
         BOXLX=BOXLX-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         BOXLY=BOXLY-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         BOXLZ=BOXLZ-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         PRINT*,'STEP=',-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         GOTO 10
      ELSE
         BOXLX=BOXLX-GRAD/SECOND
         BOXLY=BOXLY-GRAD/SECOND
         BOXLZ=BOXLZ-GRAD/SECOND
         PRINT*,'STEP=',-GRAD/SECOND
         IF (DABS(GRAD/SECOND).GT.1.0D-6) GOTO 10
      ENDIF
      RETURN
      END
C
C*************************************************************************
C
C  SUBROUTINE LJPSHIFTBIN2 CALCULATES THE ENERGY, CARTESIAN GRADIENT AND SECOND
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
C  THIS ROUTINE USES THE SAME NEIGHBOURS WHEN FIXIMAGE IS SET TO TRY AND AVOID
C  DISCONTINUITIES IN THE SECOND DERIVATIVES.
C
C
C*************************************************************************
C
      SUBROUTINE LJPSHIFTBIN2(N, X, V, POTEL, BOXLX, BOXLY, BOXLZ, CUTOFF, GTEST, STEST, PTEST)
      USE KEY
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA, ANV(N,N,3)
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3,
     1                 V(3*N), R2(N,N),  R6, R2DUM,
     2                 R8(N,N), G(N,N), EPSAB, EPSBB, SIGAB, SIGBB,
     3                 R14(N,N), F(N,N), 
     4                 POTEL, BOXLX, BOXLY, BOXLZ, SIGAB6, SIGAB12, SIGBB6, SIGBB12,
     5                 VEC(N,N,3), IRCUT2AA, IRCUT2AB, IRCUT2BB, CUTOFF, SIGRCAA6, SIGRCAA12,
     6                 SIGRCAB6, SIGRCAB12, SIGRCBB6, SIGRCBB12, CONSTAA, CONSTBB, CONSTAB,
     7                 RCONSTAA, RCONSTAB, RCONSTBB, CUTAA, CUTAB, CUTBB
      LOGICAL GTEST, STEST, PTEST, CHECKCUT(N,N)
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
C     COMMON /ANVCOMMON/ ANV
C     SAVE CHECKCUT

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
      DO J1=1,N
         CHECKCUT(J1,J1)=.FALSE.
      ENDDO
C
C  WORK OUT CUTOFF FOR POTENTIAL. TWO PARTICLES INTERACT IF R<C, BUT
C  WE WILL USE THE EQUIVALENT CONDITION 1/R^2 > 1/C^2.
C
C     IF (PTEST) PRINT*,'CUTOFF USED = ',CUTOFF
C
C  DEAL WITH ANY ATOMS THAT HAVE LEFT THE BOX.
C
C     PRINT*,'IN LJPSHIFTBIN2 FIXIMAGE,NORESET=',FIXIMAGE,NORESET
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
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(J2,J1,1)=X(J3+1)-X(J4+1)
               VEC(J2,J1,2)=X(J3+2)-X(J4+2)
               VEC(J2,J1,3)=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC(J2,J1,1)/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC(J2,J1,2)/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC(J2,J1,3)/BOXLZ)
               ENDIF
               VEC(J2,J1,1)=VEC(J2,J1,1)-BOXLX*ANV(J2,J1,1)
               VEC(J2,J1,2)=VEC(J2,J1,2)-BOXLY*ANV(J2,J1,2)
               VEC(J2,J1,3)=VEC(J2,J1,3)-BOXLZ*ANV(J2,J1,3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,NTYPEA
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2AA) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2AB) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB) ! AB
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
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2BB) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
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
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2AA) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R6=R2DUM**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
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
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2AB) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R6=R2DUM**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB)  ! AB
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
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2BB) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R6=R2DUM**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
            ENDDO
         ENDDO
      ENDIF
      POTEL = POTEL * 4.D0
     
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
      CALL LJPSHIFTGBIN2(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,CHECKCUT)
      
      IF (.NOT.STEST) RETURN
      CALL LJPSHIFTSBIN2(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,CHECKCUT)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJPSHIFTGBIN2(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,CHECKCUT)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA
      LOGICAL CHECKCUT(N,N)
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 VEC(N,N,3), V(3*N), R2(N,N),
     2                 EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,
     3                 RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
C
C  CALCULATE THE G TENSOR.
C
      DO J1=1,NTYPEA
         G(J1,J1)=0.0D0
         DO J2=J1+1,NTYPEA 
            G(J2,J1)=0.0D0
            IF (CHECKCUT(J2,J1)) G(J2,J1)=-8.0D0      *(3.0D0*(2.0D0*R14(J2,J1)        -R8(J2,J1)       )-RCONSTAA) 
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      DO J1=1,NTYPEA
         DO J2=NTYPEA+1,N 
            G(J2,J1)=0.0D0
            IF (CHECKCUT(J2,J1)) G(J2,J1)=-8.0D0*EPSAB*(3.0D0*(2.0D0*R14(J2,J1)*SIGAB12-R8(J2,J1)*SIGAB6)-RCONSTAB)
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      DO J1=NTYPEA+1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N 
            G(J2,J1)=0.0D0
            IF (CHECKCUT(J2,J1)) G(J2,J1)=-8.0D0*EPSBB*(3.0D0*(2.0D0*R14(J2,J1)*SIGBB12-R8(J2,J1)*SIGBB6)-RCONSTBB)
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO J4=1,N
               IF (CHECKCUT(J4,J1)) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
            ENDDO
         ENDDO
      ENDDO
C     DO J1=1,10
C        WRITE(*,'(A,I6,F20.10)') 'J1,V=',J1,V(J1)
C     ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPSHIFTSBIN2(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,CHECKCUT)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NTYPEA
      LOGICAL CHECKCUT(N,N)
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), IRCUT2AA, IRCUT2AB, IRCUT2BB,
     1                 VEC(N,N,3), R2(N,N), RCONSTAA, RCONSTAB, RCONSTBB,
     2                 F(N,N), EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB

C     DO J1=1,N
C        DO J2=1,N
C           WRITE(*,'(A,2I6,L)') 'J1,J2,CHECK=',J1,J2,CHECKCUT(J1,J2)
C        ENDDO
C     ENDDO

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=0.0D0
            IF (J1.LE.NTYPEA) THEN
               IF (J2.LE.NTYPEA) THEN
                  IF (CHECKCUT(J2,J1)) F(J2,J1)=8.0D0*      (84.0D0*R14(J2,J1)-24.0D0*R8(J2,J1))  
               ELSE
                  IF (CHECKCUT(J2,J1)) F(J2,J1)=8.0D0*EPSAB*(84.0D0*R14(J2,J1)*SIGAB12-24.0D0*R8(J2,J1)*SIGAB6)
               ENDIF
            ELSE
               IF (CHECKCUT(J2,J1)) F(J2,J1)=   8.0D0*EPSBB*(84.0D0*R14(J2,J1)*SIGBB12-24.0D0*R8(J2,J1)*SIGBB6)
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
            HESS(J3,J3)=0.0D0
            DO J4=1,N
               IF (CHECKCUT(J4,J1)) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
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
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO J5=1,N
                  IF (CHECKCUT(J5,J1)) HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
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
               IF (CHECKCUT(J4,J1)) THEN
                  HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
               ELSE
                  HESS(3*(J4-1)+J2,J3)=0.0D0
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
                  IF (CHECKCUT(J4,J1)) THEN
                     HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                  ELSE
                     HESS(J6,J3)=0.0D0
                  ENDIF
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

C     DO J1=1,10
C        DO J2=J1,10
C           WRITE(*,'(A,2I6,F20.10)') 'J1,J2,A=',J1,J2,HESS(J1,J2)
C        ENDDO
C     ENDDO


      RETURN
      END
