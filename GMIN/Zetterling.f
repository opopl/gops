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
C FROM FZET@NADA.KTH.SE THU AUG 24 08:56:18 2000
C HELLO!
C I AM SORRY FOR THE DELAY. 
C I HAVE ATTACHED TWO POTENTIALS IN THE DESIRED FORMAT IN 
C THIS EMAIL, THE FIRST ONE IS THE GAMMA-BRASSY ONE, THE 
C SECOND ONE IS THE "DO-NOT-CRYSTALLIZE" ONE. IT IS OF 
C COURSE THE SAME POTENTIAL, BUT WITH DIFFERENT VALUES OF 
C THE PARAMETERS. THE ANALYTICAL EXPRESSION IS:
C 
C V=A*EXP(ALPHA*R)*COS(2*K_F*R)/R**3 + B*(R/SIGMA)**POW + CONST.
C 
C THE CONSTANT IS NEEDED TO HAVE V=0 AT THE CUTOFF.
C (THE CUTOFF IS OF COURSE ANOTHER PARAMETER...)
C 
C HERE IS THE FIRST ONE:

!-----------------------------------------------------------------------*
C
C  ENERGY AND GRADIENT FOR THE MODIFIED DZUGUTOV(ZETTERLING) POTENTIAL-1.
C
C  23/08-2000
C  FREDRIK ZETTERLING <FZET@PDC.KTH.SE>
C
! CORRESPONDENCE WITH THE PARAMETERS IN THE FILE "POTENTIAL.F"
! V     = GRAD
! EDZ   = EREAL
! GTEST = GRADT
!
! MEANING OF THE VARIABLES
! DIST - DISTANCE
! EDZ  - TOTAL INTERACTION ENERGY
! VT   - (3*N) VECTOR OF INTERACTION ENERGIES OF EACH PARTICLE
! G    - (3*N,3*N) MATRIX OF "FORCES" (FIRST DERIVATIVES OF THE POTENTIAL)/R
! V    - (3*N) VECTOR - GRADIENT (FORCE ACTING ON A PARTICLE)
! 
!-----------------------------------------------------------------------*
      SUBROUTINE Z1(X,V,EDZ,GTEST)
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4
      LOGICAL GTEST, DZT
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
     1                 EDZ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY,
     2                 DDUMMY, DDDUMMY, RCUT2, NEARD(NATOMS)
      INTEGER NEAREST(NATOMS)
      PARAMETER (RCUT2=7.0176544700968D0)
      COMMON /DZ/ DZT

      DZT=.FALSE.
      EDZ=0.0D0
      DO J1=1,NATOMS 
         VT(J1)=0.0D0
         NEARD(J1)=1.0D100
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            G(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
                     DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.LT.NEARD(J2)) THEN
                  NEARD(J2)=DIST
                  NEAREST(J2)=J1   
               ENDIF
               IF (DIST.LT.NEARD(J1)) THEN
                  NEARD(J1)=DIST
                  NEAREST(J1)=J2   
               ENDIF
               G(J2,J1)=0.0D0
               IF (DIST.LT.RCUT2) THEN 
                  CALL DERPHIZ1(DIST,DUMMY,DDUMMY,DDDUMMY,.TRUE.,.FALSE.)
                  VT(J1)=VT(J1)+DUMMY
                  VT(J2)=VT(J2)+DUMMY
                  EDZ=EDZ+DUMMY
                  G(J2,J1)=DDUMMY
               ENDIF
               G(J1,J2)=G(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.LT.NEARD(J2)) THEN
                  NEARD(J2)=DIST
                  NEAREST(J2)=J1   
               ENDIF
               IF (DIST.LT.NEARD(J1)) THEN
                  NEARD(J1)=DIST
                  NEAREST(J1)=J2   
               ENDIF
               IF (DIST.LT.RCUT2) THEN 
                  CALL DERPHIZ1(DIST,DUMMY,DDUMMY,DDDUMMY,.FALSE.,.FALSE.)
                  VT(J1)=VT(J1)+DUMMY
                  VT(J2)=VT(J2)+DUMMY
                  EDZ=EDZ+DUMMY
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.GTEST) GOTO 10

      DO J1=1,NATOMS
         J3=3*J1
         IF (SEEDT.AND.(J1.GT.NATOMS-NSEED)) THEN
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
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDIF
      ENDDO

10    DO J1=1,NATOMS
C        IF (VT(J1).EQ.0.0D0) THEN
         IF (NEARD(J1).GT.2.25D0) THEN
            DZT=.TRUE.
            X(3*(J1-1)+1)=X(3*(NEAREST(J1)-1)+1)+(X(3*(J1-1)+1)-X(3*(NEAREST(J1)-1)+1))/NEARD(J1)
            X(3*(J1-1)+2)=X(3*(NEAREST(J1)-1)+2)+(X(3*(J1-1)+2)-X(3*(NEAREST(J1)-1)+2))/NEARD(J1)
            X(3*(J1-1)+3)=X(3*(NEAREST(J1)-1)+3)+(X(3*(J1-1)+3)-X(3*(NEAREST(J1)-1)+3))/NEARD(J1)
C           PRINT*,'J1,VT,DIST=',J1,VT(J1),DIST
            NEARD(NEAREST(J1))=1.0D0
         ENDIF
      ENDDO

      RETURN
      END

!%%%%%%%
! VALUES OF THE POTENTIAL, FIRST, AND SECOND DERIVATIVES
! IF SECDER = .FALSE. DDPHI = 0 ON OUTPUT
!%%%%%%%
C  ***  THE DERIVATIVE IS DIVIDED BY R !
!-----------------------------------------------------------------------*
      SUBROUTINE DERPHIZ1(R2,PHI,DPHI,DDPHI,GTEST,SECDER)
      IMPLICIT NONE
      DOUBLE PRECISION R, R2, PHI, DPHI, DDPHI
      LOGICAL SECDER, GTEST

      DOUBLE PRECISION A,B,KF,ALPH,SIG,POW,ST,RCUT ! PAIR POTENTIAL PARAMETERS
      DOUBLE PRECISION DUMMY1,DUMMY2,DUMMY3

      PARAMETER ( A=1.58 )
      PARAMETER ( B=420000000.0 )
      PARAMETER ( KF=4.12 )
      PARAMETER ( ALPH=-0.22 )
      PARAMETER ( SIG=0.331 )
      PARAMETER ( POW=-18 )
C     PARAMETER ( ST=0.046822 )
C     PARAMETER ( RCUT=2.65D0 )
      PARAMETER ( ST=0.04682632412414856D0 )
      PARAMETER ( RCUT=2.649085591311991D0 )

      R = SQRT(R2)

      DUMMY1=B*(R/SIG)**POW
      DUMMY2=A*DEXP(ALPH*R)/R**3
      DUMMY3=DCOS(2*KF*R)
      PHI=ST+DUMMY2*DUMMY3+DUMMY1

      IF (GTEST) DPHI=(DUMMY2*(DUMMY3*(ALPH-3/R)-2*KF*DSIN(2*KF*R))+POW*DUMMY1/R)/R
            
      IF (SECDER) THEN 
         DDPHI=DBLE(POW)*DBLE(POW-1)*B*((R/SIG)**POW)/(R**2)
     1      +A*(DEXP(ALPH*R)/(R**3))*((ALPH-3/R)*((ALPH-3/R)
     2      *DCOS(2*KF*R)-4*KF*DSIN(2*KF*R))+(3/(R**2)-4*KF*KF)
     3      *DCOS(2*KF*R))
      ENDIF

      END                       ! SUBROUTINE DERPHI
!-----------------------------------------------------------------------*

C AND HERE IS THE SECOND ONE:

!-----------------------------------------------------------------------*
C
C  ENERGY AND GRADIENT FOR THE MODIFIED DZUGUTOV(ZETTERLING) POTENTIAL-2.
C
C  23/08-2000
C  FREDRIK ZETTERLING <FZET@PDC.KTH.SE>
C
! CORRESPONDENCE WITH THE PARAMETERS IN THE FILE "POTENTIAL.F"
! V     = GRAD
! EDZ   = EREAL
! GTEST = GRADT
!
! MEANING OF THE VARIABLES
! DIST - DISTANCE
! EDZ  - TOTAL INTERACTION ENERGY
! VT   - (3*N) VECTOR OF INTERACTION ENERGIES OF EACH PARTICLE
! G    - (3*N,3*N) MATRIX OF "FORCES" (FIRST DERIVATIVES OF THE POTENTIAL)/R
! V    - (3*N) VECTOR - GRADIENT (FORCE ACTING ON A PARTICLE)
! 
!-----------------------------------------------------------------------*
      SUBROUTINE Z2(X,V,EDZ,GTEST)
      USE COMMONS
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
     1                 EDZ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY,
     2                 DDUMMY, DDDUMMY, RCUT2, NEARD(NATOMS)
      LOGICAL DZT
      INTEGER NEAREST(NATOMS)
      PARAMETER (RCUT2=6.995378792429925D0)
      COMMON /DZ/ DZT

      DZT=.FALSE.
      EDZ=0.0D0
      DO J1=1,NATOMS 
         VT(J1)=0.0D0
         NEARD(J1)=1.0D100
      ENDDO
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            G(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
                     DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.LT.NEARD(J2)) THEN
                  NEARD(J2)=DIST
                  NEAREST(J2)=J1   
               ENDIF
               IF (DIST.LT.NEARD(J1)) THEN
                  NEARD(J1)=DIST
                  NEAREST(J1)=J2   
               ENDIF
               G(J2,J1)=0.0D0
               IF (DIST.LT.RCUT2) THEN
                  CALL DERPHIZ2(DIST,DUMMY,DDUMMY,DDDUMMY,.TRUE.,.FALSE.)
                  VT(J1)=VT(J1)+DUMMY
                  VT(J2)=VT(J2)+DUMMY
                  EDZ=EDZ+DUMMY
                  G(J2,J1)=DDUMMY
               ENDIF
               G(J1,J2)=G(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.LT.NEARD(J2)) THEN
                  NEARD(J2)=DIST
                  NEAREST(J2)=J1   
               ENDIF
               IF (DIST.LT.NEARD(J1)) THEN
                  NEARD(J1)=DIST
                  NEAREST(J1)=J2   
               ENDIF
               IF (DIST.LT.RCUT2) THEN
                  CALL DERPHIZ2(DIST,DUMMY,DDUMMY,DDDUMMY,.FALSE.,.FALSE.)
                  VT(J1)=VT(J1)+DUMMY
                  VT(J2)=VT(J2)+DUMMY
                  EDZ=EDZ+DUMMY
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.GTEST) GOTO 10

      DO J1=1,NATOMS
         J3=3*J1
         IF (SEEDT.AND.(J1.GT.NATOMS-NSEED)) THEN
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
            V(J3-2)=DUMMYX
            V(J3-1)=DUMMYY
            V(J3)=DUMMYZ
         ENDIF
      ENDDO

10    DO J1=1,NATOMS
C        IF (VT(J1).EQ.0.0D0) THEN
         IF (NEARD(J1).GT.2.25D0) THEN
            DZT=.TRUE.
            X(3*(J1-1)+1)=X(3*(NEAREST(J1)-1)+1)+(X(3*(J1-1)+1)-X(3*(NEAREST(J1)-1)+1))/NEARD(J1)
            X(3*(J1-1)+2)=X(3*(NEAREST(J1)-1)+2)+(X(3*(J1-1)+2)-X(3*(NEAREST(J1)-1)+2))/NEARD(J1)
            X(3*(J1-1)+3)=X(3*(NEAREST(J1)-1)+3)+(X(3*(J1-1)+3)-X(3*(NEAREST(J1)-1)+3))/NEARD(J1)
C           PRINT*,'J1,VT,DIST=',J1,VT(J1),DIST
            NEARD(NEAREST(J1))=1.0D0
         ENDIF
      ENDDO

      RETURN
      END

!%%%%%%%
! VALUES OF THE POTENTIAL, FIRST, AND SECOND DERIVATIVES
! IF SECDER = .FALSE. DDPHI = 0 ON OUTPUT
!%%%%%%%
!-----------------------------------------------------------------------*
      SUBROUTINE DERPHIZ2(R2,PHI,DPHI,DDPHI,GTEST,SECDER)
      IMPLICIT NONE
      DOUBLE PRECISION R, R2, PHI, DPHI, DDPHI
      LOGICAL SECDER, GTEST

      DOUBLE PRECISION A,B,KF,ALPH,SIG,POW,ST,RCUT ! PAIR POTENTIAL PARAMETERS
      DOUBLE PRECISION DUMMY1,DUMMY2,DUMMY3

      PARAMETER ( A=1.04D0 )
      PARAMETER ( B=4200000.0D0 )
      PARAMETER ( KF=4.139D0 )
      PARAMETER ( ALPH=0.33D0 )
      PARAMETER ( SIG=0.348D0 )
      PARAMETER ( POW=-14.5D0 )
C     PARAMETER ( ST=0.133915D0 )
C     PARAMETER ( RCUT=2.645D0 )
      PARAMETER ( ST=0.1339154253770228D0 )
      PARAMETER ( RCUT=2.644877840738571D0 )

      R = SQRT(R2)

      DUMMY1=B*(R/SIG)**POW
      DUMMY2=A*DEXP(ALPH*R)/R**3
      DUMMY3=DCOS(2*KF*R)
      PHI=ST+DUMMY2*DUMMY3+DUMMY1

      IF (GTEST) DPHI=(DUMMY2*(DUMMY3*(ALPH-3/R)-2*KF*DSIN(2*KF*R))+POW*DUMMY1/R)/R

      IF (SECDER) THEN 
            
         DDPHI=DBLE(POW)*DBLE(POW-1)*B*((R/SIG)**POW)/(R**2)
     1      +A*(DEXP(ALPH*R)/(R**3))*((ALPH-3/R)*((ALPH-3/R)
     2      *DCOS(2*KF*R)-4*KF*DSIN(2*KF*R))+(3/(R**2)-4*KF*KF)
     3      *DCOS(2*KF*R))

      ENDIF
  
      END                       ! SUBROUTINE DERPHI
!-----------------------------------------------------------------------*
