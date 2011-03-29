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
C  ENERGY AND DERIVATIVES OF THE BINARY SALT POTENTIAL DESCRIBED
C  BY WELCH ET AL, JCP, 94, 4980, 1976 AND PHILLIPS ET AL, JCP,
C  94, 4980, 1991.
C  ENERGY IS IN HARTREE, LENGTH IN BOHR.
C
      SUBROUTINE WEL(X, V, POTEL, GTEST, STEST)
      USE COMMONS
      IMPLICIT NONE 
      INTEGER N, I, J
C             , J1, K, L, J2, J3
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), POTEL, XQ(NATOMS), ALPHA(NATOMS),
     1                 AC(NATOMS,NATOMS), Q(NATOMS), RRD(NATOMS,NATOMS), RRD3(NATOMS,NATOMS),
     2                 XMU(3*NATOMS), XMAT(3*NATOMS,3*NATOMS), XMINV(3*NATOMS,3*NATOMS),
     3                 RRD5(NATOMS,NATOMS)
C                      , DIF, TEMP1, V1, V2, XMU1(3*NATOMS), XMU2(3*NATOMS)
      LOGICAL GTEST, STEST

      N=NATOMS
      DO I=1,N
         IF (ZSYM(I).EQ.'MI') THEN
            Q(I)=-1.0D0
         ELSE IF (ZSYM(I).EQ.'PL') THEN
            Q(I)=1.0D0
         ELSE
            WRITE(*,'(A)') ' ALL ATOMS MUST BE TYPE PL OR MI FOR WELCH'
            STOP
         ENDIF
      ENDDO

      DO I=1,N
         IF (ZSYM(I).EQ.'PL') THEN 
            XQ(I)=XQP
            ALPHA(I)=ALPHAP
         ELSE
            XQ(I)=XQM
            ALPHA(I)=ALPHAM
         ENDIF
         DO J=I,N
            IF (ZSYM(I).EQ.'PL') THEN 
               IF (ZSYM(J).EQ.'PL') THEN 
                  AC(I,J)=APP
                  AC(J,I)=APP
               ELSE IF (ZSYM(J).EQ.'MI') THEN 
                  AC(I,J)=APM
                  AC(J,I)=APM
               ENDIF
            ELSE IF (ZSYM(I).EQ.'MI') THEN
               IF (ZSYM(J).EQ.'MI') THEN 
                  AC(I,J)=AMM
                  AC(J,I)=AMM
               ELSE IF (ZSYM(J).EQ.'PL') THEN 
                  AC(I,J)=APM
                  AC(J,I)=APM
               ENDIF
            ENDIF 
         ENDDO
      ENDDO
C
C  CALCULATE THE ENERGY.
C
      CALL WENERGY(XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,POTEL,AC,XMAT,XMINV)

      IF (.NOT.GTEST) RETURN
C
C  GRADIENT.
C
      CALL WGRAD(XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,AC,V,XMAT,XMINV)

C     PRINT*,'NUMERICAL DERIVATIVES:'
C     DIF=1.0D-4
C     DO J1=1,3*N
C        TEMP1=X(J1)
C        X(J1)=X(J1)+DIF
C        CALL WENERGY(XMU1,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,V1,AC,XMAT,XMINV)
C        X(J1)=X(J1)-2.0D0*DIF
C        CALL WENERGY(XMU2,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,V2,AC,XMAT,XMINV)
C        V(J1)=(V1-V2)/(2.0D0*DIF)
C        WRITE(*,'(I3,F20.10)') J1,V(J1)
C        X(J1)=TEMP1
C     ENDDO

      IF (.NOT.STEST) RETURN
C
C  HESSIAN.
C
      RETURN
      END
C
C*******************************************************************************
C
C  ENERGY FOR THE WELCH POTENTIAL
C
      SUBROUTINE WENERGY(XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,POTEL,AC,XMAT,XMINV)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION ESELF, ECC, EREP, ECID, EIDID, POTEL, XMU(3*NATOMS), ALPHA(NATOMS), Q(NATOMS),
     1                 RRD(NATOMS,NATOMS), X(3*NATOMS), XQ(NATOMS), DUMMY, AC(NATOMS,NATOMS),
     2                 RRD3(NATOMS,NATOMS), XMAT(3*NATOMS,3*NATOMS), XMINV(3*NATOMS,3*NATOMS),
     3                 RRD5(NATOMS,NATOMS), RD, RHOL
      INTEGER I, J, N, K, L

      N=NATOMS
      RHOL=1.0D0/RHO
C
C  CALCULATE THE INTERPARTICLE DISTANCES.
C 
      DO I=1,N
         K=1+3*(I-1)         
         RRD(I,I)=0.0D0
         RRD3(I,I)=0.0D0
         RRD5(I,I)=0.0D0
         DO J=I+1,N
            L=K+3*(J-I)
            RD= DSQRT((X(K)-X(L))*(X(K)-X(L))+
     1                (X(K+1)-X(L+1))*(X(K+1)-X(L+1))+
     2                (X(K+2)-X(L+2))*(X(K+2)-X(L+2))) 
            RRD(I,J)=1.0D0/RD
            RRD3(I,J)=RRD(I,J)**3
            RRD5(I,J)=RRD3(I,J)*RRD(I,J)**2
            RRD(J,I)=RRD(I,J)
            RRD3(J,I)=RRD3(I,J)
            RRD5(J,I)=RRD5(I,J)
         ENDDO
      ENDDO
C
C  CALCULATE THE INDUCED DIPOLES.
C
      CALL DIP(X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMAT,XMINV)

      ESELF=0.0D0
      ECC=0.0D0
      EREP=0.0D0
      ECID=0.0D0
      EIDID=0.0D0

      DO I=1,N
C
C  INDUCED DIPOLE SELF-ENERGY
C
         ESELF=ESELF+(XMU(3*(I-1)+1)**2+XMU(3*(I-1)+2)**2+XMU(3*(I-1)+3)**2)/(2.0D0*ALPHA(I))
         DO J=I+1,N
C
C  CHARGE-CHARGE.
C
            ECC=ECC+Q(I)*Q(J)*RRD(J,I)
C
C  EXPONENTIAL REPULSION.
C
            DUMMY=(X(3*(J-1)+1)-X(3*(I-1)+1)+XMU(3*(I-1)+1)/XQ(I)-XMU(3*(J-1)+1)/XQ(J))**2
     1           +(X(3*(J-1)+2)-X(3*(I-1)+2)+XMU(3*(I-1)+2)/XQ(I)-XMU(3*(J-1)+2)/XQ(J))**2
     2           +(X(3*(J-1)+3)-X(3*(I-1)+3)+XMU(3*(I-1)+3)/XQ(I)-XMU(3*(J-1)+3)/XQ(J))**2
            DUMMY=SQRT(DUMMY)
            EREP=EREP+AC(J,I)*DEXP(-RHOL*DUMMY)
C
C  CHARGE-INDUCED DIPOLE TERM. INCLUDE I,J AND J,I.
C
            DUMMY=(((X(3*(J-1)+1)-X(3*(I-1)+1))*XMU(3*(I-1)+1)+
     1              (X(3*(J-1)+2)-X(3*(I-1)+2))*XMU(3*(I-1)+2)+
     2              (X(3*(J-1)+3)-X(3*(I-1)+3))*XMU(3*(I-1)+3))*Q(J)
     3           + ((X(3*(I-1)+1)-X(3*(J-1)+1))*XMU(3*(J-1)+1)+
     1              (X(3*(I-1)+2)-X(3*(J-1)+2))*XMU(3*(J-1)+2)+
     2              (X(3*(I-1)+3)-X(3*(J-1)+3))*XMU(3*(J-1)+3))*Q(I))*RRD3(J,I)
            ECID=ECID+DUMMY
C
C  INDUCED DIPOLE-INDUCED DIPOLE TERM.
C
            EIDID=EIDID+(XMU(3*(I-1)+1)*XMU(3*(J-1)+1)
     1                  +XMU(3*(I-1)+2)*XMU(3*(J-1)+2)
     2                  +XMU(3*(I-1)+3)*XMU(3*(J-1)+3))*RRD3(J,I)
            DUMMY=((X(3*(J-1)+1)-X(3*(I-1)+1))*XMU(3*(I-1)+1)+
     1             (X(3*(J-1)+2)-X(3*(I-1)+2))*XMU(3*(I-1)+2)+
     2             (X(3*(J-1)+3)-X(3*(I-1)+3))*XMU(3*(I-1)+3))
     3           *((X(3*(J-1)+1)-X(3*(I-1)+1))*XMU(3*(J-1)+1)+
     1             (X(3*(J-1)+2)-X(3*(I-1)+2))*XMU(3*(J-1)+2)+
     2             (X(3*(J-1)+3)-X(3*(I-1)+3))*XMU(3*(J-1)+3))
            EIDID=EIDID-3.0D0*DUMMY*RRD5(J,I)
         ENDDO
      ENDDO
      POTEL=ESELF+ECC+EREP+ECID+EIDID
C     WRITE(*,'(A,F20.10)') 'POTEL=',POTEL
C     WRITE(*,'(A,F20.10)') 'ESELF=',ESELF
C     WRITE(*,'(A,F20.10)') 'ECC=',ECC
C     WRITE(*,'(A,F20.10)') 'EREP=',EREP
C     WRITE(*,'(A,F20.10)') 'ECID=',ECID
C     WRITE(*,'(A,F20.10)') 'EIDID=',EIDID

      RETURN
      END
C
C***************************************************************
C
C  CALCULATE THE INDUCED DIPOLES BY MATRIX INVERSION.
C
      SUBROUTINE DIP(X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMAT,XMINV)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION ALPHA(NATOMS), RRD(NATOMS,NATOMS), RRD3(NATOMS,NATOMS), XMU(3*NATOMS),
     1                 XMAT(3*NATOMS,3*NATOMS), XMINV(3*NATOMS,3*NATOMS), DUMMY1, XVEC(3*NATOMS),
     2                 DUMMY2, DUMMY3, DUM, D, X(3*NATOMS), Q(NATOMS), 
     3                 RRD5(NATOMS,NATOMS), DET(2), ZWORK(3*NATOMS), RCOND
      INTEGER IPVT(3*NATOMS), J1, J2, J3, J4, NTEMP
      LOGICAL SFLAG
C
C  SET UP THE MATRIX AND VECTOR.
C
      DO J1=1,NATOMS
         J3=3*(J1-1)
         XMAT(J3+1,J3+1)=1.0D0
         XMAT(J3+1,J3+2)=0.0D0
         XMAT(J3+1,J3+3)=0.0D0
         XMAT(J3+2,J3+1)=0.0D0
         XMAT(J3+2,J3+2)=1.0D0
         XMAT(J3+2,J3+3)=0.0D0
         XMAT(J3+3,J3+1)=0.0D0
         XMAT(J3+3,J3+2)=0.0D0
         XMAT(J3+3,J3+3)=1.0D0
         DO J2=J1+1,NATOMS
            J4=3*(J2-1)
            DUMMY1=RRD3(J2,J1)
            DUMMY2=RRD5(J2,J1)
            XMAT(J4+1,J3+1)=3.0D0*(X(J4+1)-X(J3+1))**2 * DUMMY2 - DUMMY1
            XMAT(J4+1,J3+2)=3.0D0*(X(J4+1)-X(J3+1))*(X(J4+2)-X(J3+2))*DUMMY2
            XMAT(J4+1,J3+3)=3.0D0*(X(J4+1)-X(J3+1))*(X(J4+3)-X(J3+3))*DUMMY2
            XMAT(J4+2,J3+1)=XMAT(J4+1,J3+2)
            XMAT(J4+2,J3+2)=3.0D0*(X(J4+2)-X(J3+2))**2 * DUMMY2 - DUMMY1
            XMAT(J4+2,J3+3)=3.0D0*(X(J4+2)-X(J3+2))*(X(J4+3)-X(J3+3))*DUMMY2
            XMAT(J4+3,J3+1)=XMAT(J4+1,J3+3)
            XMAT(J4+3,J3+2)=XMAT(J4+2,J3+3)
            XMAT(J4+3,J3+3)=3.0D0*(X(J4+3)-X(J3+3))**2 * DUMMY2 - DUMMY1

            XMAT(J3+1,J4+1)=-XMAT(J4+1,J3+1)*ALPHA(J1)*0.5D0
            XMAT(J3+1,J4+2)=-XMAT(J4+1,J3+2)*ALPHA(J1)*0.5D0
            XMAT(J3+1,J4+3)=-XMAT(J4+1,J3+3)*ALPHA(J1)*0.5D0
            XMAT(J3+2,J4+1)=-XMAT(J4+2,J3+1)*ALPHA(J1)*0.5D0
            XMAT(J3+2,J4+2)=-XMAT(J4+2,J3+2)*ALPHA(J1)*0.5D0
            XMAT(J3+2,J4+3)=-XMAT(J4+2,J3+3)*ALPHA(J1)*0.5D0
            XMAT(J3+3,J4+1)=-XMAT(J4+3,J3+1)*ALPHA(J1)*0.5D0
            XMAT(J3+3,J4+2)=-XMAT(J4+3,J3+2)*ALPHA(J1)*0.5D0
            XMAT(J3+3,J4+3)=-XMAT(J4+3,J3+3)*ALPHA(J1)*0.5D0

            XMAT(J4+1,J3+1)=-XMAT(J4+1,J3+1)*ALPHA(J2)*0.5D0
            XMAT(J4+1,J3+2)=-XMAT(J4+1,J3+2)*ALPHA(J2)*0.5D0
            XMAT(J4+1,J3+3)=-XMAT(J4+1,J3+3)*ALPHA(J2)*0.5D0
            XMAT(J4+2,J3+1)=-XMAT(J4+2,J3+1)*ALPHA(J2)*0.5D0
            XMAT(J4+2,J3+2)=-XMAT(J4+2,J3+2)*ALPHA(J2)*0.5D0
            XMAT(J4+2,J3+3)=-XMAT(J4+2,J3+3)*ALPHA(J2)*0.5D0
            XMAT(J4+3,J3+1)=-XMAT(J4+3,J3+1)*ALPHA(J2)*0.5D0
            XMAT(J4+3,J3+2)=-XMAT(J4+3,J3+2)*ALPHA(J2)*0.5D0
            XMAT(J4+3,J3+3)=-XMAT(J4+3,J3+3)*ALPHA(J2)*0.5D0
         ENDDO
      ENDDO
      
C     PRINT*,'XMAT:'
C     DO J1=1,3*NATOMS
C        WRITE(*,'(6F12.8)') (XMAT(J1,J2),J2=1,3*NATOMS)
C     ENDDO
   
      DO J1=1,NATOMS
         DUMMY1=0.0D0
         DUMMY2=0.0D0
         DUMMY3=0.0D0
         DO J2=1,NATOMS
            DUM=RRD3(J2,J1)*Q(J2)
            DUMMY1=DUMMY1-(X(3*(J1-1)+1)-X(3*(J2-1)+1))*DUM
            DUMMY2=DUMMY2-(X(3*(J1-1)+2)-X(3*(J2-1)+2))*DUM
            DUMMY3=DUMMY3-(X(3*(J1-1)+3)-X(3*(J2-1)+3))*DUM
         ENDDO
         XVEC(3*(J1-1)+1)=-DUMMY1*ALPHA(J1)*0.5D0
         XVEC(3*(J1-1)+2)=-DUMMY2*ALPHA(J1)*0.5D0
         XVEC(3*(J1-1)+3)=-DUMMY3*ALPHA(J1)*0.5D0
      ENDDO

C     PRINT*,'XVEC:'
C     DO J1=1,3*NATOMS
C        WRITE(*,'(I3,F12.8)') J1,XVEC(J1)
C     ENDDO

C
C  SET XMINV TO XMAT.
C
      DO J1=1,3*NATOMS
         XMINV(J1,J1)=1.0D0
         DO J2=J1+1,3*NATOMS
            XMINV(J2,J1)=XMAT(J2,J1)
            XMINV(J1,J2)=XMAT(J1,J2)
         ENDDO
      ENDDO

      NTEMP=3*NATOMS
      CALL DGECO(XMINV,NTEMP,NTEMP,IPVT,RCOND,ZWORK)
      CALL DGEDI(XMINV,NTEMP,NTEMP,IPVT,DET,ZWORK,11)
      
      DO J1=1,3*NATOMS
         DUMMY1=0.0D0
         DO J2=1,3*NATOMS
            DUMMY1=DUMMY1+XMINV(J1,J2)*XVEC(J2)
         ENDDO
         XMU(J1)=DUMMY1
      ENDDO

      RETURN
      END
C
C***************************************************************
C
C  CALCULATE THE ANALYTIC DERIVATIVES OF THE INDUCED DIPOLES.
C
      SUBROUTINE DIPGRAD(X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMAT,XMINV,XMUGRAD)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION ALPHA(NATOMS), RRD(NATOMS,NATOMS), RRD3(NATOMS,NATOMS), XMU(3*NATOMS),
     1                 XMAT(3*NATOMS,3*NATOMS), XMINV(3*NATOMS,3*NATOMS), DUMMY, 
     2                 DUM, X(3*NATOMS), Q(NATOMS), VEC3(3*NATOMS), RRD5(NATOMS,NATOMS),
     3                 XMUGRAD(3*NATOMS,3*NATOMS), VEC1(3*NATOMS), VEC2(3*NATOMS)
      INTEGER J1, J2, J3, J4, J5, J6, J7
      
C
C  SET UP THE DERIVATIVE MATRIX AND VECTOR FOR ATOM J5 COMPONENT J6.
C
C     PRINT*,'ANALYTIC DERIVATIVES OF MU'
      DO J5=1,NATOMS           ! J5 = I (CAN BE A OR B)
         J4=3*(J5-1)
         DO J6=1,3             ! J6 = GAMMA
            DO J1=1,3*NATOMS
               VEC2(J1)=0.0D0
            ENDDO
            DO J7=1,3          ! FIRST INDEX OF T TENSOR = ALPHA
               DO J1=1,NATOMS  ! ATOM B
                  J3=3*(J1-1)
                  DUM=-15.0D0*RRD5(J1,J5)*ALPHA(J5)*(X(J4+J7)-X(J3+J7))*(X(J4+J6)-X(J3+J6))*RRD(J1,J5)**2
                  DO J2=1,3
                     VEC1(J3+J2)=(X(J4+J2)-X(J3+J2))*DUM
                  ENDDO
                  DUM=RRD5(J1,J5)*ALPHA(J5)
                  VEC1(J3+J7)=VEC1(J3+J7)+3.0D0*(X(J4+J6)-X(J3+J6))*DUM
                  VEC1(J3+J6)=VEC1(J3+J6)+3.0D0*(X(J4+J7)-X(J3+J7))*DUM
               ENDDO
               DUMMY=0.0D0
               DO J1=1,3*NATOMS
                  DUMMY=DUMMY+VEC1(J1)*XMU(J1)
               ENDDO
               VEC2(J4+J7)=DUMMY
            ENDDO
            DO J1=1,NATOMS  ! ATOM B
               J3=3*(J1-1)
               DUM=RRD5(J1,J5)
               VEC1(J3+1)=(X(J4+1)-X(J3+1))*DUM
               VEC1(J3+2)=(X(J4+2)-X(J3+2))*DUM
               VEC1(J3+3)=(X(J4+3)-X(J3+3))*DUM
            ENDDO
            DUMMY=0.0D0
            DO J1=1,3*NATOMS
               DUMMY=DUMMY+VEC1(J1)*XMU(J1)
            ENDDO
            VEC2(J4+J6)=VEC2(J4+J6)+DUMMY*3.0D0*ALPHA(J5)
   
            DO J7=1,3          ! FIRST INDEX OF T TENSOR = ALPHA
               DO J1=1,NATOMS  ! ATOM A THIS TIME
                  J3=3*(J1-1)
                  DUM=RRD5(J1,J5)*ALPHA(J1)
                  DUMMY=(X(J4+1)-X(J3+1))*XMU(J4+1)+(X(J4+2)-X(J3+2))*XMU(J4+2)+(X(J4+3)-X(J3+3))*XMU(J4+3)
                  VEC2(J3+J7)=VEC2(J3+J7)+3.0D0*DUM*((X(J4+J6)-X(J3+J6))*XMU(J4+J7)
     1                                              +(X(J4+J7)-X(J3+J7))*(XMU(J4+J6)
     2                                  -5.0D0*DUMMY*(X(J4+J6)-X(J3+J6))*RRD(J1,J5)**2))
               ENDDO
            ENDDO
            DO J1=1,NATOMS  ! ATOM A THIS TIME
               J3=3*(J1-1)
               DUM=RRD5(J1,J5)*ALPHA(J1)
               DO J2=1,3    ! BETA INDEX OF T(AI)(2) WHERE B OR I=J1
                  DUMMY=3.0D0*(X(J4+J2)-X(J3+J2))*DUM
                  VEC2(J3+J6)=VEC2(J3+J6)+DUMMY*XMU(J4+J2)
               ENDDO
            ENDDO
C
C  NOW VEC2 SHOULD CONTAIN (D M/D X) MU APART FROM THE FACTOR OF -1/2.
C
            DO J1=1,NATOMS   ! SECOND ATOM B
               J3=3*(J1-1)
               DUM=RRD3(J1,J5)*ALPHA(J1)*Q(J5)
               DO J2=1,3     ! ALPHA INDEX OF T(IB)(1)
                  VEC3(J3+J2)=-3.0D0*(X(J4+J2)-X(J3+J2))*(X(J4+J6)-X(J3+J6))*RRD(J1,J5)**2*DUM
               ENDDO
               VEC3(J3+J6)=VEC3(J3+J6)+DUM
            ENDDO

            DO J1=1,NATOMS   ! FIRST ATOM A
               J3=3*(J1-1)
               DUM=RRD3(J1,J5)*ALPHA(J5)*Q(J1)
               DO J2=1,3     ! ALPHA INDEX OF T(AI)(1)
                  VEC3(J4+J2)= VEC3(J4+J2)+3.0D0*(X(J4+J2)-X(J3+J2))*(X(J4+J6)-X(J3+J6))*RRD(J1,J5)**2*DUM
               ENDDO
               VEC3(J4+J6)=VEC3(J4+J6)-DUM
            ENDDO
C
C  NOW VEC3 - VEC2 SHOULD CONTAIN (D Y/D X) - (D M/D X) MU APART FROM A FACTOR OF -1/2.
C
            DO J1=1,3*NATOMS
               DUMMY=0.0D0
               DO J2=1,3*NATOMS
                  DUMMY=DUMMY-0.5D0*XMINV(J1,J2)*(VEC3(J2)-VEC2(J2))
               ENDDO
               XMUGRAD(J1,J4+J6)=DUMMY
C              WRITE(*,'(2I3,F20.10)') J1,J4+J6,XMUGRAD(J1,J4+J6)
            ENDDO

         ENDDO
      ENDDO

      RETURN
      END

C
C*******************************************************************************
C
C  ANALYTIC GRADIENT FOR THE WELCH POTENTIAL
C
      SUBROUTINE WGRAD(XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,AC,GRAD,XMAT,XMINV)
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION VSELF(3*NATOMS), VCC(3*NATOMS), VREP(3*NATOMS), VCID(3*NATOMS), VIDID(3*NATOMS), 
     1                 GRAD(3*NATOMS), XMU(3*NATOMS), ALPHA(NATOMS), Q(NATOMS),
     2                 RRD(NATOMS,NATOMS), X(3*NATOMS), XQ(NATOMS), DUMMY, AC(NATOMS,NATOMS),
     3                 XMUGRAD(3*NATOMS,3*NATOMS), XMAT(3*NATOMS,3*NATOMS),
     4                 XMINV(3*NATOMS,3*NATOMS), DUM, DUM2, VD(3), RDUM, RRD3(NATOMS,NATOMS), RRD5(NATOMS,NATOMS),
     5                 T1, T2, T3, RHOL
      INTEGER N, J1, J2, J3, J4, J5, J6, J7

      N=NATOMS
      RHOL=1.0D0/RHO
C
C  CALCULATE INDUCED DIPOLE DERIVATIVES.
C
      CALL DIPGRAD(X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMAT,XMINV,XMUGRAD)
C
C  CHARGE-CHARGE.
C
      DO J1=1,N
         DO J2=1,3
            DUMMY=0.0D0
            J4=3*(J1-1)+J2
            DO J3=1,N
               DUMMY=DUMMY-Q(J3)*RRD3(J3,J1)*(X(J4)-X(3*(J3-1)+J2))
            ENDDO
            VCC(J4)=Q(J1)*DUMMY
         ENDDO
      ENDDO
C
C  SELF-ENERGY.
C
      DO J1=1,N       ! THIS IS I
         DO J3=1,3    ! THIS IS ALPHA
            J2=3*(J1-1)+J3
            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+(XMU(J5+1)*XMUGRAD(J5+1,J2)
     1                     +XMU(J5+2)*XMUGRAD(J5+2,J2)
     2                     +XMU(J5+3)*XMUGRAD(J5+3,J2))/ALPHA(J4)
            ENDDO
            VSELF(J2)=DUMMY
         ENDDO
      ENDDO
C
C  CHARGE-INDUCED DIPOLE
C
      DO J1=1,N         ! ATOM I
         J2=3*(J1-1)
         DO J3=1,3      ! INDEX GAMMA
            DUMMY=0.0D0
            DO J4=1,N   ! A
               DUMMY=DUMMY+XMU(3*(J4-1)+J3)*RRD3(J4,J1)
            ENDDO
            VCID(J2+J3)=Q(J1)*DUMMY
            DUMMY=0.0D0
            DO J4=1,N   ! B
               DUMMY=DUMMY-Q(J4)*RRD3(J4,J1)
            ENDDO
            VCID(J2+J3)=VCID(J2+J3)+DUMMY*XMU(J2+J3)
            DUMMY=0.0D0
            DO J4=1,N   ! B
               J5=3*(J4-1)
               DUMMY=DUMMY+Q(J4)*(X(J2+J3)-X(3*(J4-1)+J3))*
     1               ((X(J5+1)-X(J2+1))*XMU(J2+1)
     2               +(X(J5+2)-X(J2+2))*XMU(J2+2)
     3               +(X(J5+3)-X(J2+3))*XMU(J2+3))*RRD5(J4,J1)
            ENDDO
            VCID(J2+J3)=VCID(J2+J3)-3.0D0*DUMMY
            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+(X(J2+J3)-X(J5+J3))*
     1               ((X(J2+1)-X(J5+1))*XMU(J5+1)
     2               +(X(J2+2)-X(J5+2))*XMU(J5+2)
     3               +(X(J2+3)-X(J5+3))*XMU(J5+3))*RRD5(J4,J1)
            ENDDO
            VCID(J2+J3)=VCID(J2+J3)-3.0D0*DUMMY*Q(J1)
            DUMMY=0.0D0
            DO J4=1,N    ! A
               J6=3*(J4-1)
               T1=XMUGRAD(J6+1,J2+J3)
               T2=XMUGRAD(J6+2,J2+J3)
               T3=XMUGRAD(J6+3,J2+J3)
               DO J5=1,N ! B
                  J7=3*(J5-1)
                  DUMMY=DUMMY+Q(J5)*((X(J7+1)-X(J6+1))*T1
     1                              +(X(J7+2)-X(J6+2))*T2
     2                              +(X(J7+3)-X(J6+3))*T3)*RRD3(J5,J4)
               ENDDO
            ENDDO
            VCID(J2+J3)=VCID(J2+J3)+DUMMY
         ENDDO
      ENDDO
C
C  INDUCED DIPOLE-INDUCED DIPOLE.
C 
      DO J1=1,N
         J2=3*(J1-1)
         DO J3=1,3
            DUMMY=0.0D0
            DO J4=1,N
               J6=3*(J4-1)
               T1=XMUGRAD(J6+1,J2+J3)
               T2=XMUGRAD(J6+2,J2+J3)
               T3=XMUGRAD(J6+3,J2+J3)
               DO J5=1,N
                  J7=3*(J5-1)
                  DUMMY=DUMMY+(XMU(J7+1)*T1+XMU(J7+2)*T2+XMU(J7+3)*T3)*RRD3(J5,J4)
               ENDDO
            ENDDO
            VIDID(J2+J3)=DUMMY

            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+(XMU(J2+1)*XMU(J5+1)+XMU(J2+2)*XMU(J5+2)+XMU(J2+3)*XMU(J5+3))
     1                    *(X(J2+J3)-X(J5+J3))*RRD5(J1,J4)
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)-3.0D0*DUMMY

            DUMMY=0.0D0
            DO J4=1,N     ! A
               J6=3*(J4-1)
               T1=XMUGRAD(J6+1,J2+J3)
               T2=XMUGRAD(J6+2,J2+J3)
               T3=XMUGRAD(J6+3,J2+J3)
               DO J5=1,N  ! B
                  J7=3*(J5-1)
                  DUMMY=DUMMY+(T1*(X(J7+1)-X(J6+1))
     1                        +T2*(X(J7+2)-X(J6+2))
     2                        +T3*(X(J7+3)-X(J6+3)))
     3                       *(XMU(J7+1)*(X(J7+1)-X(J6+1))
     4                        +XMU(J7+2)*(X(J7+2)-X(J6+2))
     5                        +XMU(J7+3)*(X(J7+3)-X(J6+3)))*RRD5(J5,J4)
               ENDDO
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)-3.0D0*DUMMY

            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+XMU(J5+J3)*(XMU(J2+1)*(X(J2+1)-X(J5+1))
     1                                +XMU(J2+2)*(X(J2+2)-X(J5+2))
     2                                +XMU(J2+3)*(X(J2+3)-X(J5+3)))*RRD5(J4,J1)
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)-3.0D0*DUMMY

            DUMMY=0.0D0
            DO J4=1,N
               J5=3*(J4-1)
               DUMMY=DUMMY+(XMU(J5+1)*(X(J2+1)-X(J5+1))
     1                     +XMU(J5+2)*(X(J2+2)-X(J5+2))
     2                     +XMU(J5+3)*(X(J2+3)-X(J5+3)))*RRD5(J4,J1)
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)-3.0D0*DUMMY*XMU(J2+J3)

            DUMMY=0.0D0
            DO J4=1,N  ! B
               J5=3*(J4-1)
               DUMMY=DUMMY+(XMU(J2+1)*(X(J2+1)-X(J5+1))
     1                     +XMU(J2+2)*(X(J2+2)-X(J5+2))
     2                     +XMU(J2+3)*(X(J2+3)-X(J5+3)))
     3                    *(XMU(J5+1)*(X(J2+1)-X(J5+1))
     4                     +XMU(J5+2)*(X(J2+2)-X(J5+2))
     5                     +XMU(J5+3)*(X(J2+3)-X(J5+3)))
     6                    *(X(J2+J3)-X(J5+J3))*RRD5(J4,J1)*RRD(J4,J1)**2
            ENDDO
            VIDID(J2+J3)=VIDID(J2+J3)+15.0D0*DUMMY
         ENDDO
      ENDDO
C
C  EXPONENTIAL REPULSION.
C
      DO J1=1,N
         J2=3*(J1-1)
         DO J3=1,3
            DUMMY=0.0D0
            DO J4=1,N     ! A
               J6=3*(J4-1)
               T1=XMUGRAD(J6+1,J2+J3)
               T2=XMUGRAD(J6+2,J2+J3)
               T3=XMUGRAD(J6+3,J2+J3)
               DO J5=J4+1,N  ! B
                  J7=3*(J5-1)
                  DUM=(X(J7+1)-X(J6+1)+XMU(J6+1)/XQ(J4)-XMU(J7+1)/XQ(J5))**2
     1               +(X(J7+2)-X(J6+2)+XMU(J6+2)/XQ(J4)-XMU(J7+2)/XQ(J5))**2
     2               +(X(J7+3)-X(J6+3)+XMU(J6+3)/XQ(J4)-XMU(J7+3)/XQ(J5))**2
                  DUM=SQRT(DUM)
                  RDUM=1.0D0/DUM
                  VD(1)=T1/XQ(J4)-XMUGRAD(J7+1,J2+J3)/XQ(J5)
                  VD(2)=T2/XQ(J4)-XMUGRAD(J7+2,J2+J3)/XQ(J5)
                  VD(3)=T3/XQ(J4)-XMUGRAD(J7+3,J2+J3)/XQ(J5)
                  IF (J1.EQ.J4) VD(J3)=VD(J3)-1.0D0
                  IF (J1.EQ.J5) VD(J3)=VD(J3)+1.0D0
                  DUM2=VD(1)*(X(J7+1)-X(J6+1)+XMU(J6+1)/XQ(J4)-XMU(J7+1)/XQ(J5))
     1                +VD(2)*(X(J7+2)-X(J6+2)+XMU(J6+2)/XQ(J4)-XMU(J7+2)/XQ(J5))
     2                +VD(3)*(X(J7+3)-X(J6+3)+XMU(J6+3)/XQ(J4)-XMU(J7+3)/XQ(J5))
                  DUMMY=DUMMY+AC(J5,J4)*DEXP(-RHOL*DUM)*DUM2*RDUM
               ENDDO
            ENDDO
            VREP(J2+J3)=-DUMMY*RHOL

         ENDDO
      ENDDO

      DO J1=1,3*N
         GRAD(J1)=VCC(J1)+VSELF(J1)+VREP(J1)+VCID(J1)+VIDID(J1)
      ENDDO
C     PRINT*,'ANALYTIC GRADIENT:'
C     WRITE(*,'(I3,F20.10)') (J1,GRAD(J1),J1=1,3*NATOMS)

      RETURN
      END

      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),DET(2),WORK(*)
C
C     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     DOUBLE PRECISION(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,DSWAP
C     FORTRAN DABS,MOD
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END

      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(*)
      DOUBLE PRECISION A(LDA,*),Z(*)
      DOUBLE PRECISION RCOND
C
C     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DGEFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,DSIGN
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END

      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END

