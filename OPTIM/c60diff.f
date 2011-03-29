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
C  SUBROUTINE C60DIFF CALCULATES THE CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY. EPSILON IS SET TO UNITY AND
C  SIGMA = 1. USES F AND G TENSORS.
C
C*************************************************************************
C
      SUBROUTINE C60DIFF(N, X, V, ENERGY, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6 
      LOGICAL GTEST, STEST, RADMOVED, OVERLP
      DOUBLE PRECISION X(3*N), R2(N,N), D2DELTA,
     1                 V(3*N), F(N,N), ENERGY, TEMP, 
     2                 R(N,N), G(N,N), 
     3                 T2, T4, T3, T9, D, DSQ, D2, 
     4                 CM3, CP3, CM4, CP4, CM9, CP9, CM10, CP10, 
     5                 D2M, D2P, CM5, CP5, CM11, CP11, T6, T12 
      PARAMETER(D=1.023349668D0, DSQ=D*D, D2=2.0D0*D, D2DELTA=0.02D0) 
      COMMON /DISCON/ RADMOVED, OVERLP
C 
C  SIGMA=EPSILON=1 
C 
C  STORE DISTANCE MATRICES.  
C
      OVERLP=.FALSE.
10    CONTINUE
      DO 20 J1=1,N
         R(J1,J1)=0.0D0
         R2(J1,J1)=0.0D0
         DO J2=J1+1,N
            R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1               +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2               +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            R(J2,J1)=DSQRT(R2(J2,J1))
            R(J1,J2)=R(J2,J1)
            R2(J2,J1)=1.0D0/R2(J2,J1)
            R2(J1,J2)=R2(J2,J1)
         ENDDO
20    CONTINUE
C
C  CALCULATE THE ENERGY AND THE G AND F TENSORS.
C
      ENERGY=0.0D0
      DO 21 J1=1,N
         G(J1,J1)=0.0D0
         F(J1,J1)=0.0D0
         DO 22 J2=J1+1,N 
            TEMP=R(J2,J1)
            T2=1.0D0/(TEMP*TEMP)
            T3=T2/TEMP
            T4=T2*T2
            T6=T4*T2
            T9=T6*T3
            T12=T9*T3
            D2M=1.0D0/(-D2+TEMP)
            D2P=1.0D0/(D2+TEMP)
            CM3=D2M**3
            CP3=D2P**3
            CM4=CM3*D2M
            CP4=CP3*D2P
            CM5=CM4*D2M
            CP5=CP4*D2P
            CM9=CM5*CM4
            CP9=CP5*CP4
            CM10=CM9*D2M
            CP10=CP9*D2P
            CM11=CM10*D2M
            CP11=CP10*D2P
            G(J2,J1)=(4.0D0*((-10.0D0*(-2.0D0*T9+CM9+CP9)
     1       +75.0D0*(-2.0D0*T3+CM3+CP3))*T2
     2      +(10.0D0*(18.0D0*T9/TEMP-9.0D0*(CM10+CP10))
     3       -75.0D0*(6.0D0*T4-3.0D0*(CM4+CP4)))/TEMP))/(DSQ*TEMP)
            G(J1,J2)=G(J2,J1)

            F(J2,J1)=(-9600*T12+14400*T6
     1              +3600*(CM11+CP11-CM5-CP5)/TEMP
     3              +(1080*(CM10+CP10)-2700*(CM4+CP4))*T2
     5              +(120*(CM9+CP9)-900*(CM3+CP3))*T3)/DSQ
            F(J1,J2)=F(J2,J1)
C
C  FOR OVERLAPPING MOLECULES MAKE THE INTERACTION REPULSIVE!
C
            IF (TEMP.LT.D2+D2DELTA) THEN
               WRITE(*,'(A,G20.10,A,2I5,A)') ' WARNING DISTANCE = ',1.0D0/D2M+D2,' FOR ATOMS ',J1,J2,' IN C60DIFF, RESCALING'
               CALL RESCALEC6(X,N,D2+D2DELTA)
               OVERLP=.TRUE.
               GOTO 10
C              ENERGY=ENERGY-(10.0D0*(-2.0D0*T9+CM9+CP9)-75.0D0*(-2.0D0*T3+CM3+CP3))/TEMP
C              G(J2,J1)=-G(J2,J1)
C              G(J1,J2)=-G(J1,J2)
            ELSE
               ENERGY=ENERGY+(10.0D0*(-2.0D0*T9+CM9+CP9)-75.0D0*(-2.0D0*T3+CM3+CP3))/TEMP  
            ENDIF
22       CONTINUE
21    CONTINUE
      ENERGY=ENERGY*4.0D0/DSQ
      ENERGY=ENERGY/96.794820624738D0
C
C  FROM HERE ON DOWN THE CODE IS SYSTEM-INDEPENDENT!
C
C  FIRST THE GRADIENT.
C
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
      DO 13 J1=1,N
         DO 14 J2=1,3
            J3=3*(J1-1)+J2
            TEMP=0.0D0
            DO 15 J4=1,N
               TEMP=TEMP+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
15          CONTINUE
            V(J3)=TEMP
14       CONTINUE
13    CONTINUE
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      IF (.NOT.STEST) RETURN
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            TEMP=0.0D0
            DO 60 J4=1,N
               TEMP=TEMP+F(J4,J1)*R2(J4,J1)*
     1                 (X(J3)-X(3*(J4-1)+J2))**2 + G(J4,J1)   
60          CONTINUE
            HESS(J3,J3)=TEMP
70       CONTINUE
80    CONTINUE
C
C  NEXT ARE THE TERMS WHERE X_I AND X_J ARE ON THE SAME ATOM
C  BUT ARE DIFFERENT, E.G. Y AND Z.
C
      DO 120 J1=1,N
         DO 110 J2=1,3
            J3=3*(J1-1)+J2
            DO 100 J4=J2+1,3
               TEMP=0.0D0
               DO 90 J5=1,N
                  TEMP=TEMP+
     1            F(J5,J1)*R2(J5,J1)*
     2           (X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
90             CONTINUE
              HESS(3*(J1-1)+J4,J3)=TEMP
100         CONTINUE
110      CONTINUE
120   CONTINUE
C
C  CASE III, DIFFERENT ATOMS, SAME CARTESIAN COORDINATE.
C
      DO 150 J1=1,N
         DO 140 J2=1,3
            J3=3*(J1-1)+J2
            DO 130 J4=J1+1,N
               HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*
     1                           (X(J3)-X(3*(J4-1)+J2))**2-G(J4,J1) 
130         CONTINUE
140      CONTINUE
150   CONTINUE
C
C  CASE IV: DIFFERENT ATOMS AND DIFFERENT CARTESIAN COORDINATES.
C
      DO 180 J1=1,N
         DO 170 J2=1,3
            J3=3*(J1-1)+J2
            DO 160 J4=J1+1,N
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)
     1                    *(X(J3)-X(3*(J4-1)+J2))
     2                    *(X(3*(J1-1)+J5)-X(J6))
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
155            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
C
C  SYMMETRISE HESSIAN
C
      DO 200 J1=1,3*N
         DO 190 J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
190      CONTINUE
200   CONTINUE
      RETURN
      END

      SUBROUTINE RESCALEC6(X,NATOMS,D2)
      IMPLICIT NONE
      INTEGER NATOMS, J1, J2
      DOUBLE PRECISION X(3*NATOMS), D2, CX, CY, CZ, DX, DY, DZ, DIST

10    CONTINUE
      DO J1=1,NATOMS
         DO J2=J1+1,NATOMS
            DIST=SQRT((X(3*(J1-1)+1)-X(3*(J2-1)+1))**2+(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2+(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2)
            IF (DIST.LT.D2) THEN
               DX=(X(3*(J1-1)+1)-X(3*(J2-1)+1))*(D2+0.4D0)/DIST/2
               DY=(X(3*(J1-1)+2)-X(3*(J2-1)+2))*(D2+0.4D0)/DIST/2
               DZ=(X(3*(J1-1)+3)-X(3*(J2-1)+3))*(D2+0.4D0)/DIST/2
               CX=(X(3*(J1-1)+1)+X(3*(J2-1)+1))/2
               CY=(X(3*(J1-1)+2)+X(3*(J2-1)+2))/2
               CZ=(X(3*(J1-1)+3)+X(3*(J2-1)+3))/2

               X(3*(J1-1)+1)=CX+DX
               X(3*(J1-1)+2)=CY+DY
               X(3*(J1-1)+3)=CZ+DZ
               X(3*(J2-1)+1)=CX-DX
               X(3*(J2-1)+2)=CY-DY
               X(3*(J2-1)+3)=CZ-DZ
C              DIST=SQRT((X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
C    1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
C    2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2)
C              PRINT*,J1,J2,'SEPARATION RESET TO ',DIST
               GOTO 10
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END
