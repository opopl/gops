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
C  SUBROUTINE CAARDIFF CALCULATES THE CARTESIAN GRADIENT AND SECOND
C  DERIVATIVE MATRIX ANALYTICALLY. ATOMIC UNITS FOR CAAR_N CLUSTERS
C
C*************************************************************************
C
      SUBROUTINE CAARDIFF(N, X, V, ENERGY, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*N), ENERGY, R6,
     1                 V(3*N), R2(N,N), R2T,
     2                 R8(N,N), G(N,N),
     3                 R14(N,N), F(N,N)
      DOUBLE PRECISION R1(N),R1T
      DOUBLE PRECISION SIGAR6,EPSAR,CAAR0,CAAR1,CAAR2,CAAR3
C 
C  STORE DISTANCE MATRICES.
C

      SIGAR6=58808.93841
      EPSAR=0.00045935
      CAAR0=0.85*14.6732
      CAAR1=1.09357
      CAAR2=-0.85*606.872
      CAAR3=0.85*6.4D5

C      SIGAR6=58808.93841
C      EPSAR=0.00045935
C      CAAR0=12.472220
C      CAAR1=1.09357
C      CAAR2=-515.84120
C      CAAR3=5.55475D5

      ENERGY=0.0D0
      IF (GTEST) THEN
         DO J1=1,N
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2(J2,J1)=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1                  +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2                  +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
               IF (J1.EQ.1) THEN
                  R1(J2)=DSQRT(R2(J2,J1))
               ENDIF
               R2(J2,J1)=1.0D0/R2(J2,J1)
               R6=R2(J2,J1)**3
               IF (J1.EQ.1) THEN
                  ENERGY=ENERGY+CAAR0*DEXP(-CAAR1*R1(J2))
                  ENERGY=ENERGY+CAAR2*R6+CAAR3*R6*R6
               ELSE
                  ENERGY=ENERGY+4.D0*EPSAR*SIGAR6*R6*(SIGAR6*R6-1.D0)
               ENDIF
               R8(J2,J1)=R2(J2,J1)*R6
               R14(J2,J1)=R8(J2,J1)*R8(J2,J1)/R2(J2,J1)
               R2(J1,J2)=R2(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,N
            J3=3*(J1-1)
            DO J2=J1+1,N
               J4=3*(J2-1)
               R1T=0.D0
               R2T=(X(J3+1)-X(J4+1))**2+(X(J3+2)-X(J4+2))**2+(X(J3+3)-X(J4+3))**2
               IF (J1.EQ.1) R1T=DSQRT(R2T)
               R2T=1.0D0/R2T
               R6=R2T**3
               IF (J1.EQ.1) THEN
                  ENERGY=ENERGY+CAAR0*DEXP(-CAAR1*R1T)
                  ENERGY=ENERGY+CAAR2*R6+CAAR3*R6*R6
               ELSE
                  ENERGY=ENERGY+4.D0*EPSAR*SIGAR6*R6*(SIGAR6*R6-1.D0)
               ENDIF
            ENDDO
         ENDDO

      ENDIF

      ENERGY=ENERGY/EPSAR

      IF (.NOT.GTEST) RETURN
      CALL CAARG(G,R1,R14,R8,V,X,N)
      
      IF (.NOT.STEST) RETURN
      CALL CAARS(G,F,R1,R2,R14,R8,X,N)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE CAARG(G,R1,R14,R8,V,X,N)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 V(3*N), X(3*N), DUMMY
      DOUBLE PRECISION R1(N)
      DOUBLE PRECISION SIGAR6,EPSAR,CAAR0,CAAR1,CAAR2,CAAR3
C
C  CALCULATE THE G TENSOR.
C
      SIGAR6=58808.93841
      EPSAR=0.00045935
      CAAR0=0.85*14.6732
      CAAR1=1.09357
      CAAR2=-0.85*606.872
      CAAR3=0.85*6.4D5

C      SIGAR6=58808.93841
C      EPSAR=0.00045935
C      CAAR0=12.472220
C      CAAR1=1.09357
C      CAAR2=-515.84120
C      CAAR3=5.55475D5

      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            IF (J1.EQ.1) THEN
               G(J2,J1)=-CAAR1*CAAR0*DEXP(-CAAR1*R1(J2))/R1(J2)
               G(J2,J1)=G(J2,J1)-6.D0*CAAR2*R8(J2,J1)-12.D0*CAAR3*R14(J2,J1)
            ELSE
               G(J2,J1)=-24.D0*SIGAR6*EPSAR*(2.D0*SIGAR6*R14(J2,J1)-R8(J2,J1))
            ENDIF
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO

      DO J1=1,N
         DO J2=1,N
            G(J1,J2)=G(J1,J2)/EPSAR
         ENDDO
      ENDDO

C
C  FROM HERE ON DOWN THE CODE IS SYSTEM-INDEPENDENT!
C  FIRST CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DUMMY=0.0D0
            DO J4=1,N
               DUMMY=DUMMY+G(J4,J1)*(X(J3)-X(3*(J4-1)+J2))
            ENDDO
            V(J3)=DUMMY
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE CAARS(G,F,R1,R2,R14,R8,X,N)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, I
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 R2(N,N), F(N,N), 
     2                 X(3*N)
      DOUBLE PRECISION R1(N)
      DOUBLE PRECISION SIGAR6,EPSAR,CAAR0,CAAR1,CAAR2,CAAR3,SIGAR2
      DOUBLE PRECISION X0,Y0,Z0,XX,YY,ZZ,XY,YZ,XZ,XI,YI,ZI,XI2,YI2,ZI2,RI,RI2
      DOUBLE PRECISION W1,W11,WX1,WY1,WZ1,W2,WX2,WY2,WZ2,W3,WX3,WY3,WZ3
      DOUBLE PRECISION WXX,WYY,WZZ,WXY,WYZ,WXZ,WTOT,RIJ,XIJ,YIJ,ZIJ
      DOUBLE PRECISION W,S,WX,WY,WZ

      SIGAR2=38.887896
      SIGAR6=58808.93841
      EPSAR=0.00045935
      CAAR0=0.85*14.6732
      CAAR1=1.09357
      CAAR2=-0.85*606.872
      CAAR3=0.85*6.4D5

C      SIGAR6=58808.93841
C      EPSAR=0.00045935
C      CAAR0=12.472220
C      CAAR1=1.09357
C      CAAR2=-515.84120
C      CAAR3=5.55475D5

C     FIRST DO THE AR-AR PART

      DO J1=1,N*3
         DO J2=1,N*3
            HESS(J1,J2)=0.D0
         ENDDO
      ENDDO

C     AR-AR PART
      
      DO I=2,N
         XX=0.0
         YY=0.0
         ZZ=0.0
         XI=X(3*I-2)
         YI=X(3*I-1)
         ZI=X(3*I  )
         DO J1=2,N
            IF (I.NE.J1) THEN
               XIJ=XI-X(3*J1-2)
               YIJ=YI-X(3*J1-1)
               ZIJ=ZI-X(3*J1  )
                 RIJ=XIJ**2.+YIJ**2.+ZIJ**2.
                 W=2.*(SIGAR2/RIJ)**3-1.
                 S=(7.*(SIGAR2/RIJ)**3.-2.)/RIJ
                 WX=W-4.*S*XIJ**2.
                 WY=W-4.*S*YIJ**2.
                 WZ=W-4.*S*ZIJ**2.
                 WX=24.*EPSAR*WX*(SIGAR2**3)/(RIJ**4.)
                 WY=24.*EPSAR*WY*(SIGAR2**3)/(RIJ**4.)
                 WZ=24.*EPSAR*WZ*(SIGAR2**3)/(RIJ**4.)
                 XX=XX-WX
                 YY=YY-WY
                 ZZ=ZZ-WZ
                 
                 HESS(3*I-2,3*J1-2)=WX
                 HESS(3*I-1,3*J1-1)=WY
                 HESS(3*I,3*J1)=WZ
                 
                 HESS(3*J1,3*I)=HESS(3*I,3*J1)
                 HESS(3*J1-1,3*I-1)=HESS(3*I-1,3*J1-1)
                 HESS(3*J1-2,3*I-2)=HESS(3*I-2,3*J1-2)
                 
              ENDIF
           ENDDO
           HESS(3*I-2,3*I-2)=XX
           HESS(3*I-1,3*I-1)=YY
           HESS(3*I,3*I)=ZZ
        ENDDO

        DO I=2,N
           XY=0.0
           XZ=0.0
           YZ=0.0
           XI=X(3*I-2)
           YI=X(3*I-1)
           ZI=X(3*I)
           DO J1=2,N
              IF (I.NE.J1) THEN
                 XIJ=XI-X(3*J1-2)
                 YIJ=YI-X(3*J1-1)
                 ZIJ=ZI-X(3*J1)
                 RIJ=XIJ**2.+YIJ**2.+ZIJ**2.
                 S=SIGAR6*(7.*(SIGAR2/RIJ)**3.-2.)/(RIJ**5.)
                 WZ=EPSAR*96.*S*XIJ*YIJ
                 WY=EPSAR*96.*S*XIJ*ZIJ
                 WX=EPSAR*96.*S*YIJ*ZIJ
                 XY=XY+WZ
                 XZ=XZ+WY
                 YZ=YZ+WX

                 HESS(3*I-2,3*J1-1)=-WZ
                 HESS(3*I-2,3*J1)=-WY
                 HESS(3*I-1,3*J1)=-WX
                 HESS(3*J1,3*I-1)=HESS(3*I-1,3*J1)
                 HESS(3*J1,3*I-2)=HESS(3*I-2,3*J1)
                 HESS(3*J1-1,3*I-2)=HESS(3*I-2,3*J1-1)
                 
                 HESS(3*J1-2,3*I-1)=HESS(3*I-2,3*J1-1)
                 HESS(3*J1-2,3*I)=HESS(3*I-2,3*J1)
                 HESS(3*I,3*J1-2)=HESS(3*J1,3*I-2)
                 HESS(3*I-1,3*J1-2)=HESS(3*J1-1,3*I-2)
                 HESS(3*J1-1,3*I)=HESS(3*I-1,3*J1)
                 HESS(3*I,3*J1-1)=HESS(3*J1,3*I-1)
                 
              ENDIF
           ENDDO
           HESS(3*I-2,3*I-1)=XY
           HESS(3*I-1,3*I-2)=XY
           HESS(3*I-2,3*I)=XZ
           HESS(3*I,3*I-2)=XZ
           HESS(3*I-1,3*I)=YZ
           HESS(3*I,3*I-1)=YZ
        ENDDO

C     CA-AR PART

      X0=X(1)
      Y0=X(2)
      Z0=X(3)
      XX=0.D0
      YY=0.D0
      ZZ=0.D0
      XY=0.D0
      XZ=0.D0
      YZ=0.D0

      DO I=2,N
         XI=X0-X(3*I-2)
         YI=Y0-X(3*I-1)
         ZI=Z0-X(3*I  )
         XI2=XI**2
         YI2=YI**2
         ZI2=ZI**2
         RI2=XI2+YI2+ZI2
         XI2=XI2/RI2
         YI2=YI2/RI2
         ZI2=ZI2/RI2
         RI=DSQRT(RI2)

         W1=CAAR0*CAAR1*DEXP(-CAAR1*RI)
         W1=W1/RI
         W11=1.D0+CAAR1*RI
         WX1=W1*(XI2*W11-1.D0)
         WY1=W1*(YI2*W11-1.D0)
         WZ1=W1*(ZI2*W11-1.D0)

         W2=6.D0*CAAR2/RI2**4
         WX2=W2*(1.D0-8.D0*XI2)
         WY2=W2*(1.D0-8.D0*YI2)
         WZ2=W2*(1.D0-8.D0*ZI2)

         W3=12.D0*CAAR3/RI2**7
         WX3=W3*(1.D0-14.D0*XI2)
         WY3=W3*(1.D0-14.D0*YI2)
         WZ3=W3*(1.D0-14.D0*ZI2)

         WXX=WX1-WX2-WX3
         WYY=WY1-WY2-WY3
         WZZ=WZ1-WZ2-WZ3

C         WXX=-WXX
C         WYY=-WYY
C         WZZ=-WZZ

         XX=XX+WXX
         YY=YY+WYY
         ZZ=ZZ+WZZ

         HESS(1,3*I-2)=HESS(1,3*I-2)-WXX
         HESS(3*I-2,1)=HESS(3*I-2,1)-WXX
         HESS(2,3*I-1)=HESS(2,3*I-1)-WYY
         HESS(3*I-1,2)=HESS(3*I-1,2)-WYY
         HESS(3,3*I  )=HESS(3,3*I  )-WZZ
         HESS(3*I  ,3)=HESS(3*I  ,3)-WZZ

         HESS(3*I-2,3*I-2)=HESS(3*I-2,3*I-2)+WXX
         HESS(3*I-1,3*I-1)=HESS(3*I-1,3*I-1)+WYY
         HESS(3*I  ,3*I  )=HESS(3*I  ,3*I  )+WZZ

         W1=CAAR0*CAAR1*(CAAR1+1.D0/RI)*DEXP(-CAAR1*RI)/RI2
         W2=48.D0*CAAR2/RI2**5
         W3=168.D0*CAAR3/RI2**8
         WTOT=W1+W2+W3

         WXY=WTOT*XI*YI
         WXZ=WTOT*XI*ZI
         WYZ=WTOT*ZI*YI

         XY=XY+WXY
         YZ=YZ+WYZ
         XZ=XZ+WXZ

         HESS(1,3*I-1)=HESS(1,3*I-1)-WXY
         HESS(3*I-1,1)=HESS(3*I-1,1)-WXY
         HESS(2,3*I-2)=HESS(2,3*I-2)-WXY
         HESS(3*I-2,2)=HESS(3*I-2,2)-WXY

         HESS(1,3*I  )=HESS(1,3*I  )-WXZ
         HESS(3*I  ,1)=HESS(3*I  ,1)-WXZ
         HESS(3,3*I-2)=HESS(3,3*I-2)-WXZ
         HESS(3*I-2,3)=HESS(3*I-2,3)-WXZ

         HESS(2,3*I  )=HESS(2,3*I  )-WYZ
         HESS(3*I  ,2)=HESS(3*I  ,2)-WYZ
         HESS(3,3*I-1)=HESS(3,3*I-1)-WYZ
         HESS(3*I-1,3)=HESS(3*I-1,3)-WYZ

         HESS(3*I-2,3*I-1)=HESS(3*I-2,3*I-1)+WXY
         HESS(3*I-1,3*I-2)=HESS(3*I-1,3*I-2)+WXY
         HESS(3*I-1,3*I  )=HESS(3*I-1,3*I  )+WYZ
         HESS(3*I  ,3*I-1)=HESS(3*I  ,3*I-1)+WYZ
         HESS(3*I-2,3*I  )=HESS(3*I-2,3*I  )+WXZ
         HESS(3*I  ,3*I-2)=HESS(3*I  ,3*I-2)+WXZ

      ENDDO

      HESS(1,1)=XX
      HESS(2,2)=YY
      HESS(3,3)=ZZ
      HESS(1,2)=XY
      HESS(2,1)=XY
      HESS(1,3)=XZ
      HESS(3,1)=XZ
      HESS(2,3)=YZ
      HESS(3,2)=YZ

      DO J1=1,3*N
         DO J2=1,3*N
            HESS(J1,J2)=HESS(J1,J2)/EPSAR
         ENDDO
      ENDDO

      RETURN
      END
