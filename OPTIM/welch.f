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
C  ENERGY AND DERIVATIVES OF THE BINARY SALT POTENTIAL DESCRIBED
C  BY WELCH ET AL, JCP, 94, 4980, 1976 AND PHILLIPS ET AL, JCP,
C  94, 4980, 1991.
C  ENERGY IS IN HARTREE, LENGTH IN BOHR.
C
      SUBROUTINE WEL(N, X, V, POTEL, APP, AMM, APM, RHO, XQP, XQM, ALPHAP, ALPHAM, ZSYM, GTEST, STEST)
      USE PORFUNCS
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, I, J, J1, J2
      DOUBLE PRECISION X(3*N), V(3*N), POTEL, XQ(N), ALPHA(N),
     1                 AC(N,N), Q(N), RRD(N,N), RRD3(N,N),EDUM,
     2                 XMU(3*N), XMAT(3*N,3*N), XMINV(3*N,3*N),DIF,TEMP1,
     3                 RRD5(N,N), APP, AMM, APM, RHO, XQP, XQM, ALPHAP, ALPHAM, FG1(3*N)
      CHARACTER(LEN=5) ZSYM(N)
      LOGICAL GTEST, STEST

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
      CALL WENERGY(N,XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,POTEL,AC,XMAT,XMINV,RHO)

      IF (.NOT.GTEST) RETURN
C
C  GRADIENT.
C
      CALL WGRAD(N,XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,AC,V,XMAT,XMINV,RHO)

      IF (.NOT.STEST) RETURN
C
C  HESSIAN.
C
      DIF=1.0D-3
      DO J1=1,3*N
         TEMP1=X(J1)
         X(J1)=X(J1)+DIF
         CALL WENERGY(N,XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,EDUM,AC,XMAT,XMINV,RHO)
         CALL WGRAD(N,XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,AC,FG1,XMAT,XMINV,RHO)
         X(J1)=X(J1)-DIF
         DO J2=J1,3*N
            HESS(J2,J1)=(FG1(J2)-V(J2))/DIF
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO

      RETURN
      END
C
C*******************************************************************************
C
C  ENERGY FOR THE WELCH POTENTIAL
C
      SUBROUTINE WENERGY(N,XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,POTEL,AC,XMAT,XMINV,RHO)
      IMPLICIT NONE
      INTEGER I, J, N, K, L
      DOUBLE PRECISION ESELF, ECC, EREP, ECID, EIDID, POTEL, XMU(3*N), ALPHA(N), Q(N),
     1                 RRD(N,N), X(3*N), XQ(N), DUMMY, AC(N,N),
     2                 RRD3(N,N), XMAT(3*N,3*N), XMINV(3*N,3*N),
     3                 RRD5(N,N), RD, RHOL, RHO

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
      CALL DIP(N,X,ALPHA,RRD3,RRD5,XMU,Q,XMAT,XMINV)

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
      SUBROUTINE DIP(NATOMS,X,ALPHA,RRD3,RRD5,XMU,Q,XMAT,XMINV)
      IMPLICIT NONE
      INTEGER NATOMS, IPVT(3*NATOMS), J1, J2, J3, J4, NTEMP
      DOUBLE PRECISION ALPHA(NATOMS), RRD3(NATOMS,NATOMS), XMU(3*NATOMS),
     1                 XMAT(3*NATOMS,3*NATOMS), XMINV(3*NATOMS,3*NATOMS), DUMMY1, XVEC(3*NATOMS),
     2                 DUMMY2, DUMMY3, DUM, D, X(3*NATOMS), Q(NATOMS), 
     3                 RRD5(NATOMS,NATOMS), DET(2), ZWORK(3*NATOMS), RCOND
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
C
C  SET XMINV TO XMAT.
C
      DO J1=1,3*NATOMS
         XMINV(J1,J1)=XMAT(J1,J1)
         DO J2=J1+1,3*NATOMS
            XMINV(J2,J1)=XMAT(J2,J1)
            XMINV(J1,J2)=XMAT(J1,J2)
         ENDDO
      ENDDO

      NTEMP=3*NATOMS
      CALL DGECO(XMINV,NTEMP,NTEMP,IPVT,RCOND,ZWORK)
      CALL DGEDI(XMINV,NTEMP,NTEMP,IPVT,DET,ZWORK,11)

C     IF (ABS(DET(2)).LT.-6) PRINT '(A,G20.10)',' WELCH> WARNING - MATRIX XMAT CLOSE TO SINGULARITY, DETERMINANT=',
C    &                        DET(1)*10.0**DET(2)

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
      SUBROUTINE DIPGRAD(NATOMS,X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMINV,XMUGRAD)
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, J5, J6, J7, NATOMS
      DOUBLE PRECISION ALPHA(NATOMS), RRD(NATOMS,NATOMS), RRD3(NATOMS,NATOMS), XMU(3*NATOMS),
     1                 XMINV(3*NATOMS,3*NATOMS), DUMMY, 
     2                 DUM, X(3*NATOMS), Q(NATOMS), VEC3(3*NATOMS), RRD5(NATOMS,NATOMS),
     3                 XMUGRAD(3*NATOMS,3*NATOMS), VEC1(3*NATOMS), VEC2(3*NATOMS)
      
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
      SUBROUTINE WGRAD(N,XMU,ALPHA,Q,RRD,RRD3,RRD5,X,XQ,AC,GRAD,XMAT,XMINV,RHO)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, J7
      DOUBLE PRECISION VSELF(3*N), VCC(3*N), VREP(3*N), VCID(3*N), VIDID(3*N), 
     1                 GRAD(3*N), XMU(3*N), ALPHA(N), Q(N),
     2                 RRD(N,N), X(3*N), XQ(N), DUMMY, AC(N,N),
     3                 XMUGRAD(3*N,3*N), XMAT(3*N,3*N),
     4                 XMINV(3*N,3*N), DUM, DUM2, VD(3), RDUM, RRD3(N,N), RRD5(N,N),
     5                 T1, T2, T3, RHOL, RHO

      RHOL=1.0D0/RHO
C
C  CALCULATE INDUCED DIPOLE DERIVATIVES.
C
      CALL DIPGRAD(N,X,ALPHA,RRD,RRD3,RRD5,XMU,Q,XMINV,XMUGRAD)
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
