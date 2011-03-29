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

           SUBROUTINE ASAR1(R,EX,EX1)
C          REAL*8  FUNCTION UTN(R,I,IP)
C  TERM X0G+ AR2 (AZIZ - JCP 99 (1993) 4518)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DATA E/143.235D0/,RM/3.757D0/, A/9.03228328D0/,B/-2.37132823D0/,
     * C6,C8,C10/1.09309955D0,.51568309D0,.32521242D0/,AA/8.73933927D4/
     *,C12,C14/.27818156D0,.31111959D0/,RO/1.107D0/
     *,NCALL/1/,EKAU/3.1668D-6/,RAU/0.52918D0/
      IF(NCALL.EQ.1) THEN
      RM=RM/RAU
      E=E*EKAU
      G61=2.1D0/6.D0
      G62=0.109D0/SQRT(6.D0)
      G81=2.1D0/8.D0
      G82=0.109D0/SQRT(8.D0)
      G101=2.1D-1
      G102=0.109D0/SQRT(10.D0)
      G121=2.1D0/12.D0
      G122=0.109D0/SQRT(12.D0)
      G141=0.15D0
      G142=0.109D0/SQRT(14.D0)
      NCALL=2
                     END IF
      X=R/RM
      X2=X*X
      X4=X2*X2
      X6=X4*X2
      X8=X4*X4
      X10=X6*X4
      X12=X6*X6
      X14=X8*X6
      C6X6=C6/X6
      C8X8=C8/X8
      C10X10=C10/X10
      C12X12=C12/X12
      C14X14=C14/X14
      RR=RO*R
      F=1.D0-RR**1.68*EXP(-0.78D0*RR)
      G6S=1.D0-EXP(-G61*RR-G62*RR*RR)
      G6D=G6S*G6S
      G6T=G6D*G6D
      G6=G6T*G6D
      G8S=1.D0-EXP(-G81*RR-G82*RR*RR)
      G8D=G8S*G8S
      G8T=G8D*G8D
      G8=G8T*G8T
      G10S=1.D0-EXP(-G101*RR-G102*RR*RR)
      G10D=G10S*G10S
      G10T=G10D*G10D
      G10=G10T*G10T*G10D
      G12S=1.D0-EXP(-G121*RR-G122*RR*RR)
      G12D=G12S*G12S
      G12T=G12D*G12D
      G12=G12T*G12T*G12T
      G14S=1.D0-EXP(-G141*RR-G142*RR*RR)
      G14D=G14S*G14S
      G14T=G14D*G14D
      G14O=G14T*G14T
      G14=G14O*G14T*G14D
      VSCF=AA*EXP((-A+B*X)*X)
      VDISP=-F*(C6X6*G6+C8X8*G8+C10X10*G10+C12X12*G12+C14X14*G14)

        EX=(VSCF+VDISP)*E

      VSCF1=VSCF*(-A+2.D0*B*X)
      F1=(F-1.D0)*(1.68D0/RR-0.78D0)
      G6P=6.D0*G6/G6S*(1.D0-G6S)*(G61+2.D0*G62*RR)
      G8P=8.D0*G8/G8S*(1.D0-G8S)*(G81+2.D0*G82*RR)
      G10P=10.D0*G10/G10S*(1.D0-G10S)*(G101+2.D0*G102*RR)
      G12P=12.D0*G12/G12S*(1.D0-G12S)*(G121+2.D0*G122*RR)
      G14P=14.D0*G14/G14S*(1.D0-G14S)*(G141+2.D0*G142*RR)
      VDISP1=RO*RM*(F1*VDISP/F
     &-F*(C6X6*G6P+C8X8*G8P+C10X10*G10P+C12X12*G12P+C14X14*G14P))
     &+F*(6.D0*C6X6*G6+8.D0*C8X8*G8+10.D0*C10X10*G10
     &+12.D0*C12X12*G12+14.D0*C14X14*G14)/X

        EX1=(VSCF1+VDISP1)*E/RM

      RETURN
      END



      SUBROUTINE ENERGYAR(R,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
      IMPLICIT NONE
      DOUBLE PRECISION R,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1
      DOUBLE PRECISION R0,UT1,AR2UT0
      SAVE

      R0=R
      ESU=AR2UT0(R0,1,UT1)
      ESU1=UT1
        ESG=AR2UT0(R0,2,UT1)  
      ESG1=UT1
        EPU=AR2UT0(R0,3,UT1)  
      EPU1=UT1
        EPG=AR2UT0(R0,4,UT1)  
      EPG1=UT1
      CALL ASAR1(R,EX,EX1)
      RETURN
      END

          DOUBLE PRECISION  FUNCTION AR2UT0(RR,I,UT1)
      IMPLICIT NONE
      INTEGER NB(8)
      INTEGER, DIMENSION(8) :: NB0=(/ 17,17,17,17,17,17,17,17 /)
      INTEGER :: NF=8
      INTEGER :: NCALL=1
      INTEGER NEND(8)
C    * ,L(8)
      DOUBLE PRECISION COEF(4,33,8),RR,UT1,SOAT,EAS,RMIN,RFNBI,AAE,AE,PPVALU
      INTEGER I,NF1,J,J1,JN,JJ,N1,NN,IBCBEG,K,NENDI,NBI,NB0I
C    * ,BREAK(33,8),S(17,6)
CCC  AR2* (SPIEGELMANN,1984)
      DOUBLE PRECISION, DIMENSION(136) ::  RF= (/
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15.,
     * 4.,4.5,4.6,4.7,4.8,5.,5.25,5.5,6.,6.5,7.,7.5,8.,9.,10.,12.,15. /)
      DOUBLE PRECISION, DIMENSION(136) :: ENF=(/
     *-.55165,-.57103,-.57143,-.57120,-.57051,-.56819,-.56444,-.56055,    
     *-.55401,-.55017,-.54784,-.54667,-.54615,-.54586,-.54582,-.54573,
     *-.54551,
     *-.51143,-.53534,-.53663,-.53730,-.53754,-.53725,-.53642,-.53591,    
     *-.53664,-.53912,-.54116,-.54277,-.54392,-.54518,-.54563,-.54572,
     *-.54551,
     *-.47586,-.49780,-.49872,-.50108,-.50994,-.51996,-.52839,-.53382,    
     *-.53961,-.54275,-.54402,-.54474,-.54518,-.54565,-.54578,-.54571,
     *-.54551,
     *-.50651,-.52661,-.52702,-.52676,-.52599,-.53231,-.53737,-.54033,    
     *-.54303,-.54456,-.54498,-.54524,-.54545,-.54572,-.54580,-.54571,
     *-.54551,
     *-.54822,-.56754,-.56792,-.56765,-.56691,-.56448,-.56055,-.55647,    
     *-.54951,-.54525,-.54261,-.54120,-.54048,-.53992,-.53971,-.53941,
     *-.53889,
     *-.50609,-.52941,-.53058,-.53113,-.53123,-.53060,-.52935,-.52849,    
     *-.52893,-.53143,-.53341,-.53499,-.53617,-.53760,-.53826,-.53866,
     *-.53889,
     *-.47526,-.49699,-.49769,-.49774,-.50045,-.51213,-.52052,-.52595,    
     *-.53176,-.53493,-.53627,-.53710,-.53768,-.53839,-.53873,-.53888,
     *-.53889,
     *-.50591,-.52605,-.52647,-.52620,-.52543,-.52757,-.53254,-.53525,    
     *-.53755,-.53869,-.53888,-.53900,-.53912,-.53931,-.53937,-.53924,
     *-.53889 /)
        SAVE

      AR2UT0=0.0
      UT1=0.0
      IF(I.GT.NF) RETURN
      IF(NCALL.GT.1) GO TO 5

      NCALL=2
      SOAT=SOAT/8065/27.212
      NB(1)=1
      NF1=NF-1
      DO 3 J=1,NF1
    3 NB(J+1)=NB(J)+NB0(J)
C     EAS=0.
      EAS=ENF(NB0(1))
      RMIN=0.
      DO 2 J=1,NF
      J1=NB(J)
      JN=J1+NB0(J)-1
      NEND(J)=JN
      RMIN=MAX(RMIN,RF(J1))
C     EAS=ENF(JN)
      DO 1 JJ=J1,JN
      ENF(JJ)=ENF(JJ)-EAS 
    1 CONTINUE
    2 CONTINUE

      DO 4 J=1,NF
      N1=NB(J)
      NN=NB0(J)
      IBCBEG=0
      IF(J.EQ.3) THEN
       IBCBEG=2
       COEF(2,1,J)=0.
               END IF
      COEF(2,NN,J)=0.
      DO 7 K=1,NN
    7 COEF(1,K,J)=ENF(N1+K-1)
C     L(J)=2*NB0(J)-1
      CALL CUBSPL(RF(N1),COEF(1,1,J),NN,IBCBEG,1)
C     CALL TAUTSP(RF(NB(J)),ENF(NB(J)),NB0(J),GAM,S,BREAK(1,J),
C    *                                  COEF(1,1,J),L(J),K,IFLAG)
    4 CONTINUE

    5 CONTINUE
        NENDI=NEND(I)
      IF(RR.LT.RF(NENDI)) THEN
        NBI=NB(I)
        NB0I=NB0(I)
        RFNBI=RF(NBI)
      IF(RR.LT.RFNBI) THEN
      AAE=ENF(NBI) +0.1
      AE=-PPVALU(RF(NBI),COEF(1,1,I),NB0I,4,RFNBI,1)/AAE
      AR2UT0=AAE*EXP(-AE*(RR-RFNBI)) -0.1
      UT1=-AE*(AR2UT0 +0.1)
                  ELSE
C     UT0=PPVALU(BREAK(1,I),COEF(1,1,I),L(I),4,RR,0)
      AR2UT0=PPVALU(RF(NBI),COEF(1,1,I),NB0I,4,RR,0)
      UT1=PPVALU(RF(NBI),COEF(1,1,I),NB0I,4,RR,1) 
                  END IF
                      ELSE
      AR2UT0=ENF(NENDI)
      UT1=0.0
                      END IF
      RETURN
        END


      SUBROUTINE RGNII(N,X,GRAD,EREAL,GRADT) !  ,H0,H1,EE,EV,W,NATOMS))
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
        DOUBLE PRECISION X(3,N),H0(9*N*(N-1)/2,9*N*(N-1)/2)
     &        ,H1(9*N*(N-1)/2,9*N*(N-1)/2)
     &        ,GRAD(N*3),EE(9*N*(N-1)/2),EV(9*N*(N-1)/2,9*N*(N-1)/2)
     &        ,W(9*N*(N-1)/2,1)
      LOGICAL GRADT
        INTEGER INFO
        DOUBLE PRECISION TEMPA(9*N)
      DATA IEV/1/, NCALL/1/, SHIFT/10.0D0/
      N3=9*N*(N-1)/2

        ID=1
        IT=0
      CALL HMATD(X,N,N3,H0,H1,H2,0)
      IF(NCALL.EQ.1.OR.IT.EQ.0) THEN
         IF(ID.EQ.0) THEN
C            CALL JACOB2(H0,EV,W(1,1),EE,W(1,2),N3,IEV)
                     ELSE
C       CALL TRED2(H0,N3,N3,EE,W(1,1),IEV)
C       CALL TQLI(EE,W(1,1),N3,N3,H0,IEV)

        CALL DSYEV('V','U',N3,H0,N3,EE,TEMPA,9*N,INFO)
        IF (INFO.NE.0) THEN
           PRINT*,'WARNING - INFO=',INFO,' IN DSYEV'
        ENDIF

      DO I=1,N3
      DO J=1,N3
      EV(I,J)=H0(I,J)
      ENDDO
      ENDDO
                     END IF
      CALL SORTV(EE,EV,N3,IEV)
      EREAL=EE(N3)
      NCALL=2
                                ELSE
      DO I=1,N3
      H0(I,I)=H0(I,I)-SHIFT
      ENDDO
C       CALL ITEIG(H0,EV(1,N3),EREAL,W(1,1),W(1,2),N3)
      EREAL=EREAL+SHIFT
                                END IF

      IF(GRADT) THEN

      DO I=1,N*3
      CALL HMATD(X,N,N3,H0,H1,H2,I)
      VHVI=0.0D0
      DO J=1,N3
      VH=0.0D0
      DO K=1,N3
      VH=VH+EV(K,N3)*H1(K,J)
      ENDDO
      VHVI=VHVI+VH*EV(J,N3)
      ENDDO
      GRAD(I)=VHVI
      ENDDO

              END IF
      RETURN
      END
C ------------------------------------------------------------
C      REAL*8 FUNCTION ENERGY(R0,I,IP)
C      REAL*8 R0
C      R=R0
C      ENERGY=UT(R,I,IP)
C      RETURN
C      END
      SUBROUTINE HMATD(X,N,N3,H0,H1,H2,IGRAD)
      USE COMMONS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C        COMMON/TEST/ITEST
      DOUBLE PRECISION X(3,N),H0(9*N*(N-1)/2,9*N*(N-1)/2)
     1               ,H1(9*N*(N-1)/2,9*N*(N-1)/2)

      DATA RAUA/0.52918D0/ ,RUNIT/1.D0/
C           ISU/1/,ISG/2/,IPU/3/,IPG/4/,IX/5/


C       IF (9*N*(N-1)/2.GT.3*MXATMS) THEN
C          PRINT*,'INCREASE MXATMS FOR HMATD'
C          STOP
C       ENDIF
        IF (NEON) RUNIT=RAUA
      
      IF(IGRAD.NE.0) GO TO 10

      NN9=9*N*(N-1)/2
      DO I=1,NN9
      DO J=1,NN9
      H0(I,J)=0.0D0
      ENDDO
      ENDDO

      EXS=0.0D0
      DO K=1,N-1
      DO L=K+1,N
      RKL=
     1      SQRT((X(1,L)-X(1,K))**2+(X(2,L)-X(2,K))**2+(X(3,L)-X(3,K))**2)
CC      EX=ENERGY(RKL,IX,0)
        IF (NEON) THEN
           EX=UTN(RKL/RAUA,5,0)
        ELSE
           CALL ASAR1(RKL,EX,EX1)
        ENDIF
      EXS=EXS+EX
      ENDDO
      ENDDO

      DO I=1,N-1
      DO J=I+1,N

      IJ = (2*N-I)*(I-1)/2 +J-I
      RIJ=
     1      SQRT((X(1,I)-X(1,J))**2+(X(2,I)-X(2,J))**2+(X(3,I)-X(3,J))**2)
CC      EXIJ=ENERGY(RIJ,IX,0)
        IF (NEON) THEN
           EXIJ=UTN(RIJ/RAUA,5,0)
        ELSE
           CALL ASAR1(RIJ,EXIJ,EXIJ1)
        ENDIF

      DO KI=1,3
      DO KJ=1,3

      IJK = 9*(IJ-1) +3*(KI-1) +KJ
      H0(IJK,IJK) = H0(IJK,IJK) +EXS -EXIJ +RUNIT/RIJ

      IF(J.EQ.N) GO TO 2
      DO J0=J+1,N

      DX=X(1,J0)-X(1,J)
      DY=X(2,J0)-X(2,J)
      DZ=X(3,J0)-X(3,J)
      DXY=SQRT(DX*DX+DY*DY)
      RJJ0=SQRT(DXY*DXY+DZ*DZ)
      CT=DZ/RJJ0
      ST=DXY/RJJ0
      CF=1.0D0
      SF=0.0D0
      IF(DXY.GT.0.0D0) THEN
      CF=DX/DXY
      SF=DY/DXY
                   END IF
        IF (NEON) THEN
           CALL FENERGY(RJJ0,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ELSE
         CALL ENERGYAR(RJJ0,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ENDIF
CC      ESU=ENERGY(RJJ0,ISU,0)
CC      ESG=ENERGY(RJJ0,ISG,0)
CC      EPU=ENERGY(RJJ0,IPU,0)
CC      EPG=ENERGY(RJJ0,IPG,0)
      ESS=0.5*(ESU+ESG)
      ESD=0.5*(ESU-ESG)
      EPS=0.5*(EPU+EPG)
      EPD=0.5*(EPU-EPG)
CC      EX=ENERGY(RJJ0,IX,0)
      
      IJ0 = (2*N-I)*(I-1)/2 +J0-I

      DO KJ0=1,3

      IJK0 = 9*(IJ0-1) +3*(KI-1) +KJ0

      IF(KJ.EQ.1) THEN
      IF(KJ0.EQ.1) THEN
      HS=ESS*(ST*CF)**2+EPS*((CT*CF)**2+SF**2)
      H0(IJK,IJK)=H0(IJK,IJK)+HS-EX
      H0(IJK0,IJK0)=H0(IJK0,IJK0)+HS-EX
      HD=ESD*(ST*CF)**2+EPD*((CT*CF)**2+SF**2)      
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KJ0.EQ.2) THEN
        HS=(ESS-EPS)*ST**2*SF*CF
        H0(IJK,IJK+1)=H0(IJK,IJK+1)+HS
        H0(IJK+1,IJK)=H0(IJK+1,IJK)+HS
      HD=(ESD-EPD)*ST**2*SF*CF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KJ0.EQ.3) THEN
        HS=(ESS-EPS)*ST*CT*CF
        H0(IJK,IJK+2)=H0(IJK,IJK+2)+HS
        H0(IJK+2,IJK)=H0(IJK+2,IJK)+HS
      HD=(ESD-EPD)*ST*CT*CF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
                    END IF
      IF(KJ.EQ.2) THEN
      IF(KJ0.EQ.2) THEN
      HS=ESS*(ST*SF)**2+EPS*((CT*SF)**2+CF**2)
      H0(IJK,IJK)=H0(IJK,IJK)+HS-EX
      H0(IJK0,IJK0)=H0(IJK0,IJK0)+HS-EX
      HD=ESD*(ST*SF)**2+EPD*((CT*SF)**2+CF**2)
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KJ0.EQ.1) THEN
        HS=(ESS-EPS)*ST**2*SF*CF
        H0(IJK0,IJK0+1)=H0(IJK0,IJK0+1)+HS
        H0(IJK0+1,IJK0)=H0(IJK0+1,IJK0)+HS
      HD=(ESD-EPD)*ST**2*SF*CF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KJ0.EQ.3) THEN
        HS=(ESS-EPS)*ST*CT*SF
        H0(IJK,IJK+1)=H0(IJK,IJK+1)+HS
        H0(IJK+1,IJK)=H0(IJK+1,IJK)+HS
      HD=(ESD-EPD)*ST*CT*SF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
                    END IF
      IF(KJ.EQ.3) THEN
      IF(KJ0.EQ.3) THEN
      HS=ESS*(CT)**2   +EPS*(ST)**2
      H0(IJK,IJK)=H0(IJK,IJK)+HS-EX
      H0(IJK0,IJK0)=H0(IJK0,IJK0)+HS-EX
      HD=ESD*(CT)**2   +EPD*(ST)**2
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KJ0.EQ.1) THEN
        HS=(ESS-EPS)*ST*CT*CF
        H0(IJK0,IJK0+2)=H0(IJK0,IJK0+2)+HS
        H0(IJK0+2,IJK0)=H0(IJK0+2,IJK0)+HS
      HD=(ESD-EPD)*ST*CT*CF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KJ0.EQ.2) THEN
        HS=(ESS-EPS)*ST*CT*SF
        H0(IJK0,IJK0+1)=H0(IJK0,IJK0+1)+HS
        H0(IJK0+1,IJK0)=H0(IJK0+1,IJK0)+HS
      HD=(ESD-EPD)*ST*CT*SF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
                    END IF
      ENDDO
      ENDDO

   2      CONTINUE
      DO I0=I+1,N
      IF(I0.EQ.J) GO TO 1

      DX=X(1,I0)-X(1,I)
      DY=X(2,I0)-X(2,I)
      DZ=X(3,I0)-X(3,I)
      DXY=SQRT(DX*DX+DY*DY)
      RII0=SQRT(DXY*DXY+DZ*DZ)
      CT=DZ/RII0
      ST=DXY/RII0
      CF=1.0D0
      SF=0.0D0
      IF(DXY.GT.0.0D0) THEN
      CF=DX/DXY
      SF=DY/DXY
                   END IF
        IF (NEON) THEN
           CALL FENERGY(RII0,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ELSE
         CALL ENERGYAR(RII0,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ENDIF
CC      ESU=ENERGY(RII0,ISU,0)
CC      ESG=ENERGY(RII0,ISG,0)
CC      EPU=ENERGY(RII0,IPU,0)
CC      EPG=ENERGY(RII0,IPG,0)
      ESS=0.5*(ESU+ESG)
      ESD=0.5*(ESU-ESG)
      EPS=0.5*(EPU+EPG)
      EPD=0.5*(EPU-EPG)
CC      EX=ENERGY(RII0,IX,0)
      
      IF(I0.LT.J) IJ0 = (2*N-I0)*(I0-1)/2 +J-I0
      IF(I0.GT.J) IJ0 = (2*N-J)*(J-1)/2 +I0-J

      DO KI0=1,3

      IF(I0.LT.J) IJK0 = 9*(IJ0-1) +3*(KI0-1) +KJ
      IF(I0.GT.J) IJK0 = 9*(IJ0-1) +3*(KJ-1) +KI0

      IF(KI.EQ.1) THEN
      IF(KI0.EQ.1) THEN
      HS=ESS*(ST*CF)**2+EPS*((CT*CF)**2+SF**2)
      H0(IJK,IJK)=H0(IJK,IJK)+HS-EX
      H0(IJK0,IJK0)=H0(IJK0,IJK0)+HS-EX
      HD=ESD*(ST*CF)**2+EPD*((CT*CF)**2+SF**2)      
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KI0.EQ.2) THEN
        HS=(ESS-EPS)*ST**2*SF*CF
        H0(IJK,IJK+3)=H0(IJK,IJK+3)+HS
        H0(IJK+3,IJK)=H0(IJK+3,IJK)+HS
      HD=(ESD-EPD)*ST**2*SF*CF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KI0.EQ.3) THEN
        HS=(ESS-EPS)*ST*CT*CF
        H0(IJK,IJK+6)=H0(IJK,IJK+6)+HS
        H0(IJK+6,IJK)=H0(IJK+6,IJK)+HS
      HD=(ESD-EPD)*ST*CT*CF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
                    END IF
      IF(KI.EQ.2) THEN
      IF(KI0.EQ.2) THEN
      HS=ESS*(ST*SF)**2+EPS*((CT*SF)**2+CF**2)
      H0(IJK,IJK)=H0(IJK,IJK)+HS-EX
      H0(IJK0,IJK0)=H0(IJK0,IJK0)+HS-EX
      HD=ESD*(ST*SF)**2+EPD*((CT*SF)**2+CF**2)
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KI0.EQ.1) THEN
        HS=(ESS-EPS)*ST**2*SF*CF
      IJKD=1
      IF(I0.LT.J) IJKD=3
        H0(IJK0,IJK0+IJKD)=H0(IJK0,IJK0+IJKD)+HS
        H0(IJK0+IJKD,IJK0)=H0(IJK0+IJKD,IJK0)+HS
      HD=(ESD-EPD)*ST**2*SF*CF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KI0.EQ.3) THEN
        HS=(ESS-EPS)*ST*CT*SF
        H0(IJK,IJK+3)=H0(IJK,IJK+3)+HS
        H0(IJK+3,IJK)=H0(IJK+3,IJK)+HS
      HD=(ESD-EPD)*ST*CT*SF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
                    END IF
      IF(KI.EQ.3) THEN
      IF(KI0.EQ.3) THEN
      HS=ESS*(CT)**2   +EPS*(ST)**2
      H0(IJK,IJK)=H0(IJK,IJK)+HS-EX
      H0(IJK0,IJK0)=H0(IJK0,IJK0)+HS-EX
      HD=ESD*(CT)**2   +EPD*(ST)**2
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KI0.EQ.1) THEN
        HS=(ESS-EPS)*ST*CT*CF
      IJKD=2
      IF(I0.LT.J) IJKD=6
        H0(IJK0,IJK0+IJKD)=H0(IJK0,IJK0+IJKD)+HS
        H0(IJK0+IJKD,IJK0)=H0(IJK0+IJKD,IJK0)+HS
      HD=(ESD-EPD)*ST*CT*CF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
      IF(KI0.EQ.2) THEN
        HS=(ESS-EPS)*ST*CT*SF
      IJKD=1
      IF(I0.LT.J) IJKD=3
        H0(IJK0,IJK0+IJKD)=H0(IJK0,IJK0+IJKD)+HS
        H0(IJK0+IJKD,IJK0)=H0(IJK0+IJKD,IJK0)+HS
      HD=(ESD-EPD)*ST*CT*SF
      H0(IJK,IJK0)=HD
      H0(IJK0,IJK)=HD
                    END IF
                    END IF
      ENDDO
   1      CONTINUE
      ENDDO

      ENDDO
      ENDDO

      ENDDO
      ENDDO

C ---------------------------------------------------------------------

  10      CONTINUE
      IF(IGRAD.LE.0) GO TO 20

      NN9=9*N*(N-1)/2
      DO I=1,NN9
      DO J=1,NN9
      H1(I,J)=0.0D0
      ENDDO
      ENDDO

      M=INT(IGRAD/3)+1
      IR=MOD(IGRAD,3)
      IF(IR.EQ.0) THEN
      IR=3
      M=M-1
                  END IF

      EXS1=0.0D0
      DO K=1,N-1
      DO L=K+1,N
      IF(L.EQ.M.OR.K.EQ.M) THEN
      RKL=
     1      SQRT((X(1,L)-X(1,K))**2+(X(2,L)-X(2,K))**2+(X(3,L)-X(3,K))**2)
CC      EX1= ENERGY(RKL,IX,1)*(X(IR,L)-X(IR,K))/RKL
        IF (NEON) THEN
        EX1= UTN(RKL/RAUA,5,1)/RAUA*(X(IR,L)-X(IR,K))/RKL
        ELSE
           CALL ASAR1(RKL,EX,EX1)
        EX1= EX1*(X(IR,L)-X(IR,K))/RKL
        ENDIF
      IF(K.EQ.M) EX1=-EX1
      EXS1=EXS1+EX1
                       END IF
      ENDDO
      ENDDO

      DO I=1,N-1
      DO J=I+1,N

      IJ = (2*N-I)*(I-1)/2 +J-I

      IF(I.EQ.M.OR.J.EQ.M) THEN

      RIJ=
     1      SQRT((X(1,I)-X(1,J))**2+(X(2,I)-X(2,J))**2+(X(3,I)-X(3,J))**2)
      DRIJDX=(X(IR,I)-X(IR,J))/RIJ
      IF(J.EQ.M) DRIJDX=-DRIJDX
CC      EXIJ1= ENERGY(RIJ,IX,1)*DRIJDX
        IF (NEON) THEN
        EXIJ1= UTN(RIJ/RAUA,5,1)/RAUA*DRIJDX
        ELSE
           CALL ASAR1(RIJ,EXIJ,EXIJ1)
        EXIJ1= EXIJ1*DRIJDX
        ENDIF

                           END IF

      DO KI=1,3
      DO KJ=1,3

      IJK = 9*(IJ-1) +3*(KI-1) +KJ
      H1(IJK,IJK) = H1(IJK,IJK) +EXS1

      IF(I.EQ.M.OR.J.EQ.M) THEN

      H1(IJK,IJK) = H1(IJK,IJK) -EXIJ1 -RUNIT/(RIJ*RIJ)*DRIJDX

                           END IF

      IF(J.EQ.N) GO TO 22
      DO J0=J+1,N

      IF(J.EQ.M.OR.J0.EQ.M) THEN

      DX=X(1,J0)-X(1,J)
      DY=X(2,J0)-X(2,J)
      DZ=X(3,J0)-X(3,J)
      DXY=SQRT(DX*DX+DY*DY)
      RJJ0=SQRT(DXY*DXY+DZ*DZ)
      CT=DZ/RJJ0
      ST=DXY/RJJ0
      CF=1.0D0
      SF=0.0D0
      IF(DXY.GT.0.0D0) THEN
      CF=DX/DXY
      SF=DY/DXY
                   END IF

      DRDX= (X(IR,J0)-X(IR,J))/RJJ0
      IF(J.EQ.M) DRDX=-DRDX

      CT1=-CT/RJJ0*DRDX
      IF(IR.EQ.3) THEN
      IF(M.EQ.J0) CT1=CT1+1.0D0/RJJ0
      IF(M.EQ.J) CT1=CT1-1.0D0/RJJ0
                END IF
      ST1=-ST/RJJ0*DRDX
      DXYDX=0.0D0
      IF(IR.EQ.1.OR.IR.EQ.2) THEN
C      DXYDX=1.0D0
      IF(DXY.GT.0.0D0) THEN
      DXYDX= (X(IR,J0)-X(IR,J))/DXY
      IF(J.EQ.M) DXYDX=-DXYDX
                   END IF
      ST1=ST1+DXYDX/RJJ0
                             END IF
      CF1=0.0D0
      SF1=0.0D0
      IF(DXY.GT.0.0D0) THEN
      IF(IR.EQ.1.OR.IR.EQ.2) THEN
      CF1=-CF/DXY*DXYDX
      IF(IR.EQ.1) THEN
      IF(M.EQ.J0) CF1=CF1+1.0D0/DXY
      IF(M.EQ.J) CF1=CF1-1.0D0/DXY
                    END IF
      SF1=-SF/DXY*DXYDX
      IF(IR.EQ.2) THEN
      IF(M.EQ.J0) SF1=SF1+1.0D0/DXY
      IF(M.EQ.J) SF1=SF1-1.0D0/DXY
                    END IF
                             END IF
                   END IF

        IF (NEON) THEN
           CALL FENERGY(RJJ0,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ELSE
         CALL ENERGYAR(RJJ0,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ENDIF
CC      ESU=ENERGY(RJJ0,ISU,0)
CC      ESG=ENERGY(RJJ0,ISG,0)
CC      EPU=ENERGY(RJJ0,IPU,0)
CC      EPG=ENERGY(RJJ0,IPG,0)
      ESS=0.5*(ESU+ESG)
      ESD=0.5*(ESU-ESG)
      EPS=0.5*(EPU+EPG)
      EPD=0.5*(EPU-EPG)
CC      EX=ENERGY(RJJ0,IX,0)

CC      ESU1=ENERGY(RJJ0,ISU,1)*DRDX
CC      ESG1=ENERGY(RJJ0,ISG,1)*DRDX
CC      EPU1=ENERGY(RJJ0,IPU,1)*DRDX
CC      EPG1=ENERGY(RJJ0,IPG,1)*DRDX
        ESU1=ESU1*DRDX
        ESG1=ESG1*DRDX
        EPU1=EPU1*DRDX
        EPG1=EPG1*DRDX
      ESS1=0.5*(ESU1+ESG1)
      ESD1=0.5*(ESU1-ESG1)
      EPS1=0.5*(EPU1+EPG1)
      EPD1=0.5*(EPU1-EPG1)
CC      EX1=ENERGY(RJJ0,IX,1)*DRDX
      EX1=EX1*DRDX
      
      IJ0 = (2*N-I)*(I-1)/2 +J0-I

      DO KJ0=1,3

      IJK0 = 9*(IJ0-1) +3*(KI-1) +KJ0

      IF(KJ.EQ.1) THEN
      IF(KJ0.EQ.1) THEN
           TRIG1=(ST*CF)**2
           TRIG2=(CT*CF)**2+SF**2
           TRIG3=2.0D0*ST*CF*(CF*ST1+ST*CF1)
           TRIG4=(2.0D0*CT*CF*(CF*CT1+CT*CF1)+2.0D0*SF*SF1)
CC      HS1=ESS1*(ST*CF)**2+EPS1*((CT*CF)**2+SF**2) +
CC     1          ESS*2.0D0*ST*CF*(CF*ST1+ST*CF1)+
CC     1      EPS*(2.0D0*CT*CF*(CF*CT1+CT*CF1)+2.0D0*SF*SF1)
        HS1= ESS1*TRIG1+EPS1*TRIG2+ESS*TRIG3+EPS*TRIG4
      H1(IJK,IJK)=H1(IJK,IJK)+HS1-EX1
      H1(IJK0,IJK0)=H1(IJK0,IJK0)+HS1-EX1
CC      HD1=ESD1*(ST*CF)**2+EPD1*((CT*CF)**2+SF**2) +
CC     1          ESD*2.0D0*ST*CF*(CF*ST1+ST*CF1)+
CC     1      EPD*(2.0D0*CT*CF*(CF*CT1+CT*CF1)+2.0D0*SF*SF1)
        HD1= ESD1*TRIG1 + EPD1*TRIG2 + ESD*TRIG3 + EPD*TRIG4
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KJ0.EQ.2) THEN
           TRIG1=ST**2*SF*CF
           TRIG2=2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1)
CC        HS1=(ESS1-EPS1)*ST**2*SF*CF +
CC     1      (ESS-EPS)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1)) 
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(IJK,IJK+1)=H1(IJK,IJK+1)+HS1
        H1(IJK+1,IJK)=H1(IJK+1,IJK)+HS1
CC      HD1=(ESD1-EPD1)*ST**2*SF*CF + 
CC     1      (ESD-EPD)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1))
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KJ0.EQ.3) THEN
           TRIG1=ST*CT*CF
           TRIG2=ST1*CT*CF+ST*CT1*CF+ST*CT*CF1
CC        HS1=(ESS1-EPS1)*ST*CT*CF + 
CC     1      (ESS-EPS)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(IJK,IJK+2)=H1(IJK,IJK+2)+HS1
        H1(IJK+2,IJK)=H1(IJK+2,IJK)+HS1
CC      HD1=(ESD1-EPD1)*ST*CT*CF + 
CC     1      (ESD-EPD)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
                    END IF

      IF(KJ.EQ.2) THEN
      IF(KJ0.EQ.2) THEN
           TRIG1=(ST*SF)**2
           TRIG2=(CT*SF)**2+CF**2
           TRIG3=2.0D0*ST*SF*(SF*ST1+ST*SF1)
           TRIG4=2.0D0*CT*SF*(SF*CT1+CT*SF1)+2.0D0*CF*CF1
CC      HS1=ESS1*(ST*SF)**2+EPS1*((CT*SF)**2+CF**2) +
CC     1          ESS*2.0D0*ST*SF*(SF*ST1+ST*SF1)+
CC     1      EPS*(2.0D0*CT*SF*(SF*CT1+CT*SF1)+2.0D0*CF*CF1)
        HS1=ESS1*TRIG1+EPS1*TRIG2+ESS*TRIG3+EPS*TRIG4
      H1(IJK,IJK)=H1(IJK,IJK)+HS1-EX1
      H1(IJK0,IJK0)=H1(IJK0,IJK0)+HS1-EX1
CC      HD1=ESD1*(ST*SF)**2+EPD1*((CT*SF)**2+CF**2) +
CC     1          ESD*2.0D0*ST*SF*(SF*ST1+ST*SF1)+
CC     1      EPD*(2.0D0*CT*SF*(SF*CT1+CT*SF1)+2.0D0*CF*CF1)
        HD1=ESD1*TRIG1+EPD1*TRIG2+ESD*TRIG3+EPD*TRIG4
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KJ0.EQ.1) THEN
           TRIG1=ST**2*SF*CF
           TRIG2=2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1)
CC        HS1=(ESS1-EPS1)*ST**2*SF*CF +
CC     1      (ESS-EPS)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1))
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(IJK0,IJK0+1)=H1(IJK0,IJK0+1)+HS1
        H1(IJK0+1,IJK0)=H1(IJK0+1,IJK0)+HS1
CC      HD1=(ESD1-EPD1)*ST**2*SF*CF + 
CC     1      (ESD-EPD)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1))
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KJ0.EQ.3) THEN
           TRIG1=ST*CT*SF
           TRIG2=ST1*CT*SF+ST*CT1*SF+ST*CT*SF1
CC        HS1=(ESS1-EPS1)*ST*CT*SF +
CC     1      (ESS-EPS)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(IJK,IJK+1)=H1(IJK,IJK+1)+HS1
        H1(IJK+1,IJK)=H1(IJK+1,IJK)+HS1
CC      HD1=(ESD1-EPD1)*ST*CT*SF + 
CC     1      (ESD-EPD)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
                    END IF

      IF(KJ.EQ.3) THEN
      IF(KJ0.EQ.3) THEN
           TRIG1=(CT)**2
           TRIG2=(ST)**2
           TRIG3=2.0D0*CT*CT1
           TRIG4=2.0D0*ST*ST1
CC      HS1=ESS1*(CT)**2   +EPS1*(ST)**2 +
CC     1      ESS*2.0D0*CT*CT1   +EPS*2.0D0*ST*ST1
        HS1=ESS1*TRIG1 + EPS1*TRIG2 + ESS*TRIG3 + EPS*TRIG4
      H1(IJK,IJK)=H1(IJK,IJK)+HS1-EX1
      H1(IJK0,IJK0)=H1(IJK0,IJK0)+HS1-EX1
CC      HD1=ESD1*(CT)**2   +EPD1*(ST)**2 +
CC     1      ESD*2.0D0*CT*CT1   +EPD*2.0D0*ST*ST1
        HD1=ESD1*TRIG1 + EPD1*TRIG2 + ESD*TRIG3 + EPD*TRIG4
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KJ0.EQ.1) THEN
           TRIG1=ST*CT*CF
           TRIG2=ST1*CT*CF+ST*CT1*CF+ST*CT*CF1
CC        HS1=(ESS1-EPS1)*ST*CT*CF + 
CC     1      (ESS-EPS)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(IJK0,IJK0+2)=H1(IJK0,IJK0+2)+HS1
        H1(IJK0+2,IJK0)=H1(IJK0+2,IJK0)+HS1
CC      HD1=(ESD1-EPD1)*ST*CT*CF + 
CC     1      (ESD-EPD)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KJ0.EQ.2) THEN
           TRIG1=ST*CT*SF
           TRIG2=ST1*CT*SF+ST*CT1*SF+ST*CT*SF1
CC        HS1=(ESS1-EPS1)*ST*CT*SF + 
CC     1      (ESS-EPS)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(IJK0,IJK0+1)=H1(IJK0,IJK0+1)+HS1
        H1(IJK0+1,IJK0)=H1(IJK0+1,IJK0)+HS1
CC      HD1=(ESD1-EPD1)*ST*CT*SF + 
CC     1      (ESD-EPD)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
                    END IF

      ENDDO
                       END IF
      ENDDO

  22      CONTINUE
      DO I0=I+1,N
      IF(I0.EQ.J) GO TO 11

      IF(I.EQ.M.OR.I0.EQ.M) THEN

      DX=X(1,I0)-X(1,I)
      DY=X(2,I0)-X(2,I)
      DZ=X(3,I0)-X(3,I)
      DXY=SQRT(DX*DX+DY*DY)
      RII0=SQRT(DXY*DXY+DZ*DZ)
      CT=DZ/RII0
      ST=DXY/RII0
      CF=1.0D0
      SF=0.0D0
      IF(DXY.GT.0.0D0) THEN
      CF=DX/DXY
      SF=DY/DXY
                   END IF

      DRDX= (X(IR,I0)-X(IR,I))/RII0
      IF(I.EQ.M) DRDX=-DRDX

      CT1=-CT/RII0*DRDX
      IF(IR.EQ.3) THEN
      IF(M.EQ.I0) CT1=CT1+1.0D0/RII0
      IF(M.EQ.I) CT1=CT1-1.0D0/RII0
                END IF
      ST1=-ST/RII0*DRDX
      DXYDX=0.0D0
      IF(IR.EQ.1.OR.IR.EQ.2) THEN
C      DXYDX=1.0D0
      IF(DXY.GT.0.0D0) THEN
      DXYDX= (X(IR,I0)-X(IR,I))/DXY
      IF(I.EQ.M) DXYDX=-DXYDX
                   END IF
      ST1=ST1+DXYDX/RII0
                             END IF
      CF1=0.0D0
      SF1=0.0D0
      IF(DXY.GT.0.0D0) THEN
      IF(IR.EQ.1.OR.IR.EQ.2) THEN
      CF1=-CF/DXY*DXYDX
      IF(IR.EQ.1) THEN
      IF(M.EQ.I0) CF1=CF1+1.0D0/DXY
      IF(M.EQ.I) CF1=CF1-1.0D0/DXY
                    END IF
      SF1=-SF/DXY*DXYDX
      IF(IR.EQ.2) THEN
      IF(M.EQ.I0) SF1=SF1+1.0D0/DXY
      IF(M.EQ.I) SF1=SF1-1.0D0/DXY
                    END IF
                             END IF
                   END IF

        IF (NEON) THEN
           CALL FENERGY(RII0,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ELSE
         CALL ENERGYAR(RII0,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ENDIF
CC      ESU=ENERGY(RII0,ISU,0)
CC      ESG=ENERGY(RII0,ISG,0)
CC      EPU=ENERGY(RII0,IPU,0)
CC      EPG=ENERGY(RII0,IPG,0)
      ESS=0.5*(ESU+ESG)
      ESD=0.5*(ESU-ESG)
      EPS=0.5*(EPU+EPG)
      EPD=0.5*(EPU-EPG)
CC      EX=ENERGY(RII0,IX,0)

CC      ESU1=ENERGY(RII0,ISU,1)*DRDX
CC      ESG1=ENERGY(RII0,ISG,1)*DRDX
CC      EPU1=ENERGY(RII0,IPU,1)*DRDX
CC      EPG1=ENERGY(RII0,IPG,1)*DRDX
        ESU1=ESU1*DRDX
        ESG1=ESG1*DRDX
        EPU1=EPU1*DRDX
        EPG1=EPG1*DRDX
      ESS1=0.5*(ESU1+ESG1)
      ESD1=0.5*(ESU1-ESG1)
      EPS1=0.5*(EPU1+EPG1)
      EPD1=0.5*(EPU1-EPG1)
CC      EX1=ENERGY(RII0,IX,1)*DRDX
      EX1=EX1*DRDX      

      IF(I0.LT.J) IJ0 = (2*N-I0)*(I0-1)/2 +J-I0
      IF(I0.GT.J) IJ0 = (2*N-J)*(J-1)/2 +I0-J

      DO KI0=1,3

      IF(I0.LT.J) IJK0 = 9*(IJ0-1) +3*(KI0-1) +KJ
      IF(I0.GT.J) IJK0 = 9*(IJ0-1) +3*(KJ-1) +KI0

      IF(KI.EQ.1) THEN
      IF(KI0.EQ.1) THEN
           TRIG1=(ST*CF)**2
           TRIG2=(CT*CF)**2+SF**2
           TRIG3=2.0D0*ST*CF*(CF*ST1+ST*CF1)
           TRIG4=(2.0D0*CT*CF*(CF*CT1+CT*CF1)+2.0D0*SF*SF1)
CC      HS1=ESS1*(ST*CF)**2+EPS1*((CT*CF)**2+SF**2) +
CC     1          ESS*2.0D0*ST*CF*(CF*ST1+ST*CF1)+
CC     1      EPS*(2.0D0*CT*CF*(CF*CT1+CT*CF1)+2.0D0*SF*SF1)
        HS1= ESS1*TRIG1+EPS1*TRIG2+ESS*TRIG3+EPS*TRIG4
      H1(IJK,IJK)=H1(IJK,IJK)+HS1-EX1
      H1(IJK0,IJK0)=H1(IJK0,IJK0)+HS1-EX1
CC      HD1=ESD1*(ST*CF)**2+EPD1*((CT*CF)**2+SF**2) +
CC     1          ESD*2.0D0*ST*CF*(CF*ST1+ST*CF1)+
CC     1      EPD*(2.0D0*CT*CF*(CF*CT1+CT*CF1)+2.0D0*SF*SF1)
        HD1= ESD1*TRIG1 + EPD1*TRIG2 + ESD*TRIG3 + EPD*TRIG4
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KI0.EQ.2) THEN
           TRIG1=ST**2*SF*CF
           TRIG2=2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1)
CC        HS1=(ESS1-EPS1)*ST**2*SF*CF +
CC     1      (ESS-EPS)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1)) 
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(IJK,IJK+3)=H1(IJK,IJK+3)+HS1
        H1(IJK+3,IJK)=H1(IJK+3,IJK)+HS1
CC      HD1=(ESD1-EPD1)*ST**2*SF*CF + 
CC     1      (ESD-EPD)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1))
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KI0.EQ.3) THEN
           TRIG1=ST*CT*CF
           TRIG2=ST1*CT*CF+ST*CT1*CF+ST*CT*CF1
CC        HS1=(ESS1-EPS1)*ST*CT*CF + 
CC     1      (ESS-EPS)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(IJK,IJK+6)=H1(IJK,IJK+6)+HS1
        H1(IJK+6,IJK)=H1(IJK+6,IJK)+HS1
CC      HD1=(ESD1-EPD1)*ST*CT*CF + 
CC     1      (ESD-EPD)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
                    END IF

      IF(KI.EQ.2) THEN
      IF(KI0.EQ.2) THEN
           TRIG1=(ST*SF)**2
           TRIG2=(CT*SF)**2+CF**2
           TRIG3=2.0D0*ST*SF*(SF*ST1+ST*SF1)
           TRIG4=2.0D0*CT*SF*(SF*CT1+CT*SF1)+2.0D0*CF*CF1
CC      HS1=ESS1*(ST*SF)**2+EPS1*((CT*SF)**2+CF**2) +
CC     1          ESS*2.0D0*ST*SF*(SF*ST1+ST*SF1)+
CC     1      EPS*(2.0D0*CT*SF*(SF*CT1+CT*SF1)+2.0D0*CF*CF1)
        HS1=ESS1*TRIG1+EPS1*TRIG2+ESS*TRIG3+EPS*TRIG4
      H1(IJK,IJK)=H1(IJK,IJK)+HS1-EX1
      H1(IJK0,IJK0)=H1(IJK0,IJK0)+HS1-EX1
CC      HD1=ESD1*(ST*SF)**2+EPD1*((CT*SF)**2+CF**2) +
CC     1          ESD*2.0D0*ST*SF*(SF*ST1+ST*SF1)+
CC     1      EPD*(2.0D0*CT*SF*(SF*CT1+CT*SF1)+2.0D0*CF*CF1)
        HD1=ESD1*TRIG1+EPD1*TRIG2+ESD*TRIG3+EPD*TRIG4
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KI0.EQ.1) THEN
           TRIG1=ST**2*SF*CF
           TRIG2=2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1)
CC        HS1=(ESS1-EPS1)*ST**2*SF*CF +
CC     1      (ESS-EPS)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1))
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
      IJKD=1
      IF(I0.LT.J) IJKD=3
        H1(IJK0,IJK0+IJKD)=H1(IJK0,IJK0+IJKD)+HS1
        H1(IJK0+IJKD,IJK0)=H1(IJK0+IJKD,IJK0)+HS1
CC      HD1=(ESD1-EPD1)*ST**2*SF*CF + 
CC     1      (ESD-EPD)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1))
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KI0.EQ.3) THEN
           TRIG1=ST*CT*SF
           TRIG2=ST1*CT*SF+ST*CT1*SF+ST*CT*SF1
CC        HS1=(ESS1-EPS1)*ST*CT*SF +
CC     1      (ESS-EPS)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(IJK,IJK+3)=H1(IJK,IJK+3)+HS1
        H1(IJK+3,IJK)=H1(IJK+3,IJK)+HS1
CC      HD1=(ESD1-EPD1)*ST*CT*SF + 
CC     1      (ESD-EPD)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
                    END IF

      IF(KI.EQ.3) THEN
      IF(KI0.EQ.3) THEN
           TRIG1=(CT)**2
           TRIG2=(ST)**2
           TRIG3=2.0D0*CT*CT1
           TRIG4=2.0D0*ST*ST1
CC      HS1=ESS1*(CT)**2   +EPS1*(ST)**2 +
CC     1      ESS*2.0D0*CT*CT1   +EPS*2.0D0*ST*ST1
        HS1=ESS1*TRIG1 + EPS1*TRIG2 + ESS*TRIG3 + EPS*TRIG4
      H1(IJK,IJK)=H1(IJK,IJK)+HS1-EX1
      H1(IJK0,IJK0)=H1(IJK0,IJK0)+HS1-EX1
CC      HD1=ESD1*(CT)**2   +EPD1*(ST)**2 +
CC     1      ESD*2.0D0*CT*CT1   +EPD*2.0D0*ST*ST1
        HD1=ESD1*TRIG1 + EPD1*TRIG2 + ESD*TRIG3 + EPD*TRIG4
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KI0.EQ.1) THEN
           TRIG1=ST*CT*CF
           TRIG2=ST1*CT*CF+ST*CT1*CF+ST*CT*CF1
CC        HS1=(ESS1-EPS1)*ST*CT*CF + 
CC     1      (ESS-EPS)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
      IJKD=2
      IF(I0.LT.J) IJKD=6
        H1(IJK0,IJK0+IJKD)=H1(IJK0,IJK0+IJKD)+HS1
        H1(IJK0+IJKD,IJK0)=H1(IJK0+IJKD,IJK0)+HS1
CC      HD1=(ESD1-EPD1)*ST*CT*CF + 
CC     1      (ESD-EPD)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
      IF(KI0.EQ.2) THEN
           TRIG1=ST*CT*SF
           TRIG2=ST1*CT*SF+ST*CT1*SF+ST*CT*SF1
CC        HS1=(ESS1-EPS1)*ST*CT*SF + 
CC     1      (ESS-EPS)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
      IJKD=1
      IF(I0.LT.J) IJKD=3
        H1(IJK0,IJK0+IJKD)=H1(IJK0,IJK0+IJKD)+HS1
        H1(IJK0+IJKD,IJK0)=H1(IJK0+IJKD,IJK0)+HS1
CC      HD1=(ESD1-EPD1)*ST*CT*SF + 
CC     1      (ESD-EPD)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(IJK,IJK0)=HD1
      H1(IJK0,IJK)=HD1
                    END IF
                    END IF

      ENDDO
                       END IF
  11      CONTINUE
      ENDDO

      ENDDO
      ENDDO

      ENDDO
      ENDDO

  20      CONTINUE

      RETURN
      END
      SUBROUTINE GRND(N,X,GRAD,EREAL,GRADT)
      USE COMMONS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(3,N),GRAD(N*3)
      LOGICAL GRADT

      SAVE
      DATA RAUA/0.52918D0/

      EXS=0.0D0
      DO K=1,N-1
      DO L=K+1,N
      RKL=
     1      SQRT((X(1,L)-X(1,K))**2+(X(2,L)-X(2,K))**2+(X(3,L)-X(3,K))**2)
        IF (NEON) THEN
           EX=UTN(RKL/RAUA,5,0)
        ELSE
           CALL ASAR1(RKL,EX,EX1)
        ENDIF
      EXS=EXS+EX
      ENDDO
      ENDDO
      EREAL=EXS


      IF(GRADT) THEN

      DO M=1,N
      DO IR=1,3
      IC=(M-1)*3+IR

      EXS1=0.0D0
      DO K=1,N-1
      DO L=K+1,N
      IF(L.EQ.M.OR.K.EQ.M) THEN
      RKL=
     1      SQRT((X(1,L)-X(1,K))**2+(X(2,L)-X(2,K))**2+(X(3,L)-X(3,K))**2)
        IF (NEON) THEN
      EX1= UTN(RKL/RAUA,5,1)/RAUA*(X(IR,L)-X(IR,K))/RKL
      IF(K.EQ.M) EX1=-EX1
        ELSE
           CALL ASAR1(RKL,EX,EX1)
        EX1= EX1*(X(IR,L)-X(IR,K))/RKL
        IF(K.EQ.M) EX1=-EX1
        ENDIF
      EXS1=EXS1+EX1
                       END IF
      ENDDO
      ENDDO

      GRAD(IC)=EXS1

      ENDDO
      ENDDO

                  END IF

      RETURN
      END
             SUBROUTINE CUBSPL(TAU,C,N,IBCBEG,IBCEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IBCBEG,IBCEND,N,I,J,L,M
      DOUBLE PRECISION C(4,N),TAU(N),DIVDF1,DIVDF3,DTAU,G
      L=N-1
      DO 10 M=2,N
      C(3,M)=TAU(M)-TAU(M-1)
   10 C(4,M)=(C(1,M)-C(1,M-1))/C(3,M)
      IF(IBCBEG-1) 11,15,16
   11 IF(N.GT.2) GO TO 12
      C(4,1)=1.D0
      C(3,1)=1.D0
      C(2,1)=2.D0*C(4,2)
      GO TO 25
   12 C(4,1)=C(3,3)
      C(3,1)=C(3,2)+C(3,3)
      C(2,1)=((C(3,2)+2.D0*C(3,1))*C(4,2)*C(3,3)+C(3,2)**2*C(4,3))
     1      /C(3,1)
      GO TO 19
   15 C(4,1)=1.D0
      C(3,1)=0.D0
      GO TO 18
   16 C(4,1)=2.D0
      C(3,1)=1.D0
      C(2,1)=3.D0*C(4,2)-C(3,2)/2.D0*C(2,1)
   18 IF(N.EQ.2) GO TO 25
   19 DO 20 M=2,L
      G=-C(3,M+1)/C(4,M-1)
      C(2,M)=G*C(2,M-1)+3.D0*(C(3,M)*C(4,M+1)+C(3,M+1)*C(4,M))
   20 C(4,M)=G*C(3,M-1)+2.D0*(C(3,M)+C(3,M+1))
      IF(IBCEND-1) 21,30,24
   21 IF(N.EQ.3.AND.IBCBEG.EQ.0) GO TO 22
      G=C(3,N-1)+C(3,N)
      C(2,N)=( (C(3,N)+2.D0*G)*C(4,N)*C(3,N-1)+C(3,N)**2*
     1         (C(1,N-1)-C(1,N-2))/C(3,N-1) )/G
      G=-G/C(4,N-1)
      C(4,N)=C(3,N-1)
      GO TO 29
   22 C(2,N)=2.D0*C(4,N)
      C(4,N)=1.D0
      GO TO 28
   24 C(2,N)=3.D0*C(4,N)+C(3,N)/2.D0*C(2,N)
      C(4,N)=2.D0
      GO TO 28
   25 IF(IBCEND-1) 26,30,24
   26 IF(IBCBEG.GT.0) GO TO 22
      C(2,N)=C(4,N)
      GO TO 30
   28 G=-1.D0/C(4,N-1)
   29 C(4,N)=G*C(3,N-1)+C(4,N)
      C(2,N)=(G*C(2,N-1)+C(2,N))/C(4,N)
   30 DO 40 J1=1,L
      J=L+1-J1
   40 C(2,J)=(C(2,J)-C(3,J)*C(2,J+1))/C(4,J)
      DO 50 I=2,N
      DTAU=C(3,I)
      DIVDF1=(C(1,I)-C(1,I-1))/DTAU
      DIVDF3=C(2,I-1)+C(2,I)-2.D0*DIVDF1
      C(3,I-1)=2.D0*(DIVDF1-C(2,I-1)-DIVDF3)/DTAU
   50 C(4,I-1)=(DIVDF3/DTAU)*(6.D0/DTAU)
      RETURN
      END
      SUBROUTINE INTERV( XT, LXT, X, LEFT, MFLAG)
COMPUTES LEFT = MAX( I, 1 .LE. I .LE. LXT .AND. XT(I) .LE. X)
C
C******  I N P U T  ******
C  XT  ... A REAL SEQUENCE, OF LENGTH  LXT, ASSUMED TO BE NONDECREASING
C  LXT ... NUMBER OF TERMS IN THE SEQUENCE  XT.
C  X   ... THE POINT WHOSE LOCATION WITH RESPECT TO THE SEQUENCE  XT
C          IS TO BE DETERMINED.
C
C******  O U T P U T  ******
C  LEFT, MFLAG ... BOTH INTEGER, WHOSE VALUE IS
C
C   1    -1      IF  XT(LXT) .LE. X .LT.  XT(1)
C   1     0      IF   XT(I)  .LE. X .LT. XT(I+1)
C  LXT    1      IF  XT(LXT) .LE. X
C
C     IN PARTICULAR, MFLAG = 0 IS THE 'USUAL' CASE. MFLAG .NE. 0
C     INDICATES THAT  X  LIES OUTSIDE THE HALFOPEN INTERVAL
C     XT(1) .LE. Y .LT. XT(LXT) . THE ASSYMETRIC TREATMENT OF THE
C     INTERVAL IS DUE TO THE DECISION TO MAKE ALL  PP  FUNCTIONS CONT-
C     INUOUS FROM THE RIGHT.
C
C******  M E T H O D  ******
C  THE PROGRAM IS DESIGNED TO BE EFFICIENT IN THE COMMON SITUATION THAT
C  IT IS CALLED REPEATEDLY, WITH  X  TAKEN FROM AN INCREASING OR DECREA-
C  SING SEQUENCE. THIS WILL HAPPEN, E.G., WHEN A PP FUNCTION IS TO BE
C  GRAPHED. THE FIRST GUESS FOR  LEFT  IS THEREFORE TAKEN TO BE THE VAL-
C  UE RETURNED AT THE PREVIOUS CALL AND STORED IN THE  L O C A L  VARIA
C  BLE  ILO . A FIRST CHECK ASCERTAIN THAT  ILO .LT. LXT (THIS IS NEC-
C  ESSARY SINCE THE PRESENT CALL MAY HAVE NOTHING TO DO WITH THE PREVI-
C  OUS CALL). THEN, IF  XT(ILO) .LE. X .LT. XT(ILO+1), WE SET  LEFT =
C  ILO  AND ARE DONE AFTER JUST THREE COMPARISONS.
C     OTHERWISE, WE REPEATEDLY DOUBLE THE DIFFERENCE  ISTEP = IHI - ILO
C  WHILE ALSO MOVING  ILO  IND  IHI  IN THE DIRECTION OF  X , UNTIL
C             XT(ILO) .LE. X .LT. XT(IHI) ,
C  AFTER WHICH WE USE BISECTION TO GET, IN ADDITION, ILO+1 = IHI .
C  LEFT = ILO IS THEN RETURNED.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LEFT,LXT,MFLAG,
     *        IHI,ILO,ISTEP,MIDDLE
      DOUBLE PRECISION  X,XT(LXT)
      SAVE
      DATA ILO /1/
C     SAVE ILO  (A VALID FORTRAN STATEMENT IN THE NEW 1977 STANDARD)
      IHI = ILO + 1
      IF (IHI .LT. LXT)                GO TO 20
        IF (X .GE. XT(LXT))            GO TO 110
          IF (LXT .LE. 1)              GO TO 90
            ILO = LXT - 1
            IHI = LXT
C
   20 IF (X .GE. XT(IHI))              GO TO 40
        IF (X .GE. XT(ILO))            GO TO 100
C      **** NOW X .LT. XT(ILO) . DECREASE  ILO  TO CAPTURE  X .
            ISTEP = 1
   31    IHI = ILO
         ILO = IHI - ISTEP
         IF (ILO .LE. 1)               GO TO 35
           IF (X .GE. XT(ILO))         GO TO 50
             ISTEP = ISTEP*2
                                       GO TO 31
   35 ILO = 1
      IF (X .LT. XT(1))                GO TO 90
                                       GO TO 50
C     **** NOW X .GE. XT(IHI) . INCREASE  IHI TO CAPTURE  X .
   40 ISTEP = 1
   41   ILO = IHI
        IHI = ILO + ISTEP
        IF (IHI .GE. LXT)              GO TO 45
          IF (X .LT. XT(IHI))          GO TO 50
            ISTEP =ISTEP*2
                                       GO TO 41
   45 IF (X .GE. XT(LXT))              GO TO 110
        IHI = LXT
C
C       **** NOW XT(ILO) .LE. X .LT. XT(IHI) . NPRROW THE INTERVAL.
   50 MIDDLE = (ILO + IHI)/2
      IF (MIDDLE .EQ. ILO)             GO TO 100
C     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO + 1 .
         IF (X .LT. XT(MIDDLE))         GO TO 53
           ILO = MIDDLE
                                        GO TO 50
   53      IHI = MIDDLE
                                        GO TO 50
C**** SET OUTPUT AND RETURN.
   90 MFLAG = -1
      LEFT  =  1
                                        GO TO 99999
  100 MFLAG =  0
      LEFT  = ILO
                                        GO TO 99999
  110 MFLAG =  1
      LEFT  = LXT
99999                                   RETURN
      END

      DOUBLE PRECISION FUNCTION PPVALU( BREAK, COEF, L, K, X, JDERIV)
C   CALLS INTERV
CALCULATES VALUE AT  X  OF JDERIV-TH DERIVATIVE OF  PP  FROM PP-REPR
C
C******  I N P U T  ******
C  BREAK, COEF, L, K ... FORMS THE PP-REPRESENTATION OF THE FUNCTION  F
C         TO BE EVALUATED.SPECIFICALLY, THE J-TH DERIVATIVE OF  F  IS
C         GIVEN BY
C
C   D(**J)F(X) = COEF(J+1,I) + H*(COEF(J+2,I) + H*( ... (COEF(K-1,I)+
C                            + H* COEF(  K,I) / (K-J-1))/(K-J-2) . )/2)
C
C   WITH  H = X - BREAK(I) , AND
C
C   I = MAX ( 1, MAX( J, BREAK(J) .LE. X, 1 .LE. J .LE. L ) ).
C
C  X      ... THE POINT AT WHICH TO EVALUATE.
C  JDERIV ... INTEGER GIVING THE ORDER OF THE DERIVATIVE TO BE EVALUATE-
C             ED. A S S U M E D  TO BE ZERO OR POSITIVE.
C
C******  O U T P U T  ******
C  PPVALU ... THE VALUE OF THE (JDERIV)-TH DERIVATIVE OF  F  AT  X.
C
C******  M E T H O D  ******
C     THE INTERVAL INDEX  I, APPROPRIATE FOR X, IS FOUND THROUGH A
C  CALL TO INTERV. THE FORMULA ABOVE FOR THE  JDERIV-TH DERIVATIVE
C  OF  F  IS THEN EVALUATED (BY NESTED MULTIPLICATION).
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      INTEGER JDERIV,K,L,
     *        I,M,NDUMMY,JDRP1
      DOUBLE PRECISION    BREAK(L),COEF(K,L)
C      DOUBLE PRECISION  FMMJDR,H,PPSUM,FDBLE,YD,X,FSNGL,YS
      PPSUM = 0.0
      FMMJDR = (FLOAT(K - JDERIV))
C     DERIVATIVES OF ORDER  K  OR HIGHER ARE IDENTICALLY ZERO.
      IF (FMMJDR .LE. 0.0)             GO TO 99
C     FIND INDEX  I  OF LARGEST BREAKPOINT TO THE LEFT OF  X.
      CALL INTERV( BREAK, L, X, I, NDUMMY)
C     EVALUATE JDERIV-TH DERIVATIVE OF  I-TH  POLYNOMIAL PIECE AT  X.
      H = X - BREAK(I)
      JDRP1 = JDERIV + 1
      M = K
      DO 10 NDUMMY=JDRP1,K
         PPSUM = (PPSUM/FMMJDR)*H + COEF(M,I)
         FMMJDR = FMMJDR - 1.0
   10    M = M - 1
   99 PPVALU = PPSUM
                                       RETURN
      END

        SUBROUTINE DIPOLES(X,N,H0,H1,IGRAD)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC      SAVE
      COMMON/DIPOL/ ADIP,RUNIT,REX

        DOUBLE PRECISION X(3,N),H0(N*3,N*3)
     &               ,H1(N*3,N*3)

      IF(IGRAD.GT.0) GO TO 10

      DO I=1,N

        X1I=X(1,I)/RUNIT
        X2I=X(2,I)/RUNIT
        X3I=X(3,I)/RUNIT
      EDDI=0.D0

      DO J=1,N-1
      IF(J.EQ.I) GO TO 1

      X1J=X(1,J)/RUNIT
        X2J=X(2,J)/RUNIT
        X3J=X(3,J)/RUNIT
      X1IJ=X1J-X1I
      X2IJ=X2J-X2I
      X3IJ=X3J-X3I
      RIJ2= X1IJ*X1IJ + X2IJ*X2IJ + X3IJ*X3IJ
      RIJ=SQRT(RIJ2)
      DIJ=ADIP/RIJ2

      DO K=J+1,N
      IF(K.EQ.I.OR.K.EQ.J) GO TO 2 

        X1K=X(1,K)/RUNIT
        X2K=X(2,K)/RUNIT
        X3K=X(3,K)/RUNIT
        X1IK=X1K-X1I
        X2IK=X2K-X2I
        X3IK=X3K-X3I
        RIK2= X1IK*X1IK + X2IK*X2IK + X3IK*X3IK        
        RIK=SQRT(RIK2) 
        DIK=ADIP/RIK2

        X1JK=X1K-X1J                                                   
        X2JK=X2K-X2J                                                   
        X3JK=X3K-X3J                                                   
        RJK2= X1JK*X1JK + X2JK*X2JK + X3JK*X3JK                                 
      RJK=SQRT(RJK2)

      IF(RJK.LT.0.5*REX) GO TO 2

      DIJDIK=DIJ*DIK*(X1IJ*X1IK + X2IJ*X2IK + X3IJ*X3IK)/(RIJ*RIK)
        DIJRJK=DIJ*(X1IJ*X1JK + X2IJ*X2JK + X3IJ*X3JK)/RIJ
        DIKRJK=DIK*(X1IK*X1JK + X2IK*X2JK + X3IK*X3JK)/RIK          

      EDDJKI=(RJK2*DIJDIK -3.D0*DIJRJK*DIKRJK)/(RJK2*RJK2*RJK)
      EDDI=EDDI+EDDJKI

  2     CONTINUE
        ENDDO

  1      CONTINUE
        ENDDO

        I3=3*(I-1)
        DO L=1,3
        I3L=I3+L
        H0(I3L,I3L)=H0(I3L,I3L)+EDDI
        ENDDO

      ENDDO

  10      CONTINUE
C----------------------------------------------------------------------
      IF(IGRAD.LE.0) GO TO 20

        M=INT(IGRAD/3)+1
        IR=MOD(IGRAD,3)
        IF(IR.EQ.0) THEN
        IR=3
        M=M-1
                    END IF

        DO I=1,N

        X1I=X(1,I)/RUNIT
        X2I=X(2,I)/RUNIT
        X3I=X(3,I)/RUNIT
        EDDI1=0.D0

        DO J=1,N-1
        IF(J.EQ.I) GO TO 11

        X1J=X(1,J)/RUNIT
        X2J=X(2,J)/RUNIT
        X3J=X(3,J)/RUNIT
        X1IJ=X1J-X1I
        X2IJ=X2J-X2I
        X3IJ=X3J-X3I
        RIJ2= X1IJ*X1IJ + X2IJ*X2IJ + X3IJ*X3IJ
        RIJ=SQRT(RIJ2)
        DIJ=ADIP/RIJ2

      RIJ1=0.D0
      DIJ1=0.D0
      IF(I.EQ.M.OR.J.EQ.M) THEN
      RIJ1=(X(IR,J)-X(IR,I))/RIJ/RUNIT
      IF(I.EQ.M) RIJ1=-RIJ1
      DIJ1=-2.D0*DIJ/RIJ*RIJ1
                       END IF

        DO K=J+1,N
        IF(K.EQ.I.OR.K.EQ.J) GO TO 22

        X1K=X(1,K)/RUNIT
        X2K=X(2,K)/RUNIT
        X3K=X(3,K)/RUNIT
        X1IK=X1K-X1I
        X2IK=X2K-X2I
        X3IK=X3K-X3I
        RIK2= X1IK*X1IK + X2IK*X2IK + X3IK*X3IK
        RIK=SQRT(RIK2)
        DIK=ADIP/RIK2

      RIK1=0.D0
      DIK1=0.D0
        IF(I.EQ.M.OR.K.EQ.M) THEN
        RIK1=(X(IR,K)-X(IR,I))/RIK/RUNIT
        IF(I.EQ.M) RIK1=-RIK1      
        DIK1=-2.D0*DIK/RIK*RIK1
                             END IF

        X1JK=X1K-X1J
        X2JK=X2K-X2J
        X3JK=X3K-X3J
        RJK2= X1JK*X1JK + X2JK*X2JK + X3JK*X3JK
        RJK=SQRT(RJK2)

      IF (RJK.LT.0.5*REX) GO TO 22

      XXJKI=X1IJ*X1IK + X2IJ*X2IK + X3IJ*X3IK
        DIJDIK=DIJ*DIK*XXJKI/(RIJ*RIK)
      XXIKJ=X1IJ*X1JK + X2IJ*X2JK + X3IJ*X3JK
        DIJRJK=DIJ*XXIKJ/RIJ
      XXIJK=X1IK*X1JK + X2IK*X2JK + X3IK*X3JK
        DIKRJK=DIK*XXIJK/RIK

      RJK1=0.D0
      IF(J.EQ.M.OR.K.EQ.M) THEN
        RJK1=(X(IR,K)-X(IR,J))/RJK/RUNIT
        IF(J.EQ.M) RJK1=-RJK1
                       END IF

      XXJKI1=0.D0
      DIJDIK1=0.D0
        IF(I.EQ.M.OR.J.EQ.M.OR.K.EQ.M) THEN
      IF(J.EQ.M) XXJKI1=(X(IR,K)-X(IR,I))/RUNIT
        IF(K.EQ.M) XXJKI1=(X(IR,J)-X(IR,I))/RUNIT
      IF(I.EQ.M) XXJKI1=-(X(IR,K)-X(IR,I)+X(IR,J)-X(IR,I))/RUNIT
      DIJDIK1=(DIJ1*DIK*XXJKI+DIJ*DIK1*XXJKI+DIJ*DIK*XXJKI1)/(RIJ*RIK)
     &             -DIJDIK*(RIJ1/RIJ+RIK1/RIK)
                                          END IF
        XXIKJ1=0.D0
        DIJRJK1=0.D0
        IF(I.EQ.M.OR.J.EQ.M.OR.K.EQ.M) THEN
        IF(I.EQ.M) XXIKJ1=-(X(IR,K)-X(IR,J))/RUNIT
        IF(K.EQ.M) XXIKJ1=(X(IR,J)-X(IR,I))/RUNIT
        IF(J.EQ.M) XXIKJ1=(X(IR,K)-X(IR,J)-(X(IR,J)-X(IR,I)))/RUNIT
        DIJRJK1=(DIJ1*XXIKJ+DIJ*XXIKJ1)/RIJ
     &         -DIJRJK*RIJ1/RIJ
                                       END IF
        XXIJK1=0.D0  
        DIKRJK1=0.D0
        IF(I.EQ.M.OR.J.EQ.M.OR.K.EQ.M) THEN
        IF(I.EQ.M) XXIJK1=-(X(IR,K)-X(IR,J))/RUNIT
        IF(J.EQ.M) XXIJK1=-(X(IR,K)-X(IR,I))/RUNIT
        IF(K.EQ.M) XXIJK1=(X(IR,K)-X(IR,J)+X(IR,K)-X(IR,I))/RUNIT 
        DIKRJK1=(DIK1*XXIJK+DIK*XXIJK1)/RIK                             
     &         -DIKRJK*RIK1/RIK                            
                                       END IF


        EDDJKI=(RJK2*DIJDIK -3.D0*DIJRJK*DIKRJK)/(RJK2*RJK2*RJK)
        EDDJKI1=(RJK2*DIJDIK1 -3.D0*(DIJRJK1*DIKRJK+DIJRJK*DIKRJK1)
     &             +2.D0*RJK*RJK1*DIJDIK)/(RJK2*RJK2*RJK)
     &             - 5.D0*EDDJKI/RJK*RJK1
        EDDI1=EDDI1+EDDJKI1

  22    CONTINUE
        ENDDO

  11    CONTINUE
        ENDDO

        I3=3*(I-1)
        DO L=1,3
        I3L=I3+L
        H1(I3L,I3L)=H1(I3L,I3L)+EDDI1/RUNIT
        ENDDO

        ENDDO

  20    CONTINUE

      RETURN
      END


C   2  2  2  2  2
C NE2+SU.DAT_QNE2+SG.DAT_QNE2+PU.DAT_QNE2+PG.DAT_QNE2X.DAT_Q
C NE2+SU.DAT  NE2+SG.DAT  NE2+PU.DAT  NE2+PG.DAT  NE2X.DAT   

      SUBROUTINE FENERGY(R,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON/DIPOL/ ADIP,RUNIT,REX

CC        REAL R0,UT,UT1
      SAVE
      DATA RAUA/0.52918D0/

      RUNIT=RAUA
      ADIP=2.68D0

        R0=R
        ESU=UT(R0,1,UT1)
        ESU1=UT1
        ESG=UT(R0,2,UT1)
        ESG1=UT1
        EPU=UT(R0,3,UT1)
        EPU1=UT1
        EPG=UT(R0,4,UT1)
        EPG1=UT1
      EX =UTN(R/RAUA,5,0)
      EX1=UTN(R/RAUA,5,1)/RAUA
        RETURN
        END

      DOUBLE PRECISION FUNCTION UT(RR,I0,UT1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON/POINTS/RU,EU, RG,EG, RPU,EPU, RPG,EPG, RX,EX,
     1              COEFU,COEFG,COEFPU,COEFPG,COEFX, IR1,IR2,IR3,IR4,IR5
      CHARACTER(LEN=12) NAME1,NAME2,NAME3,NAME4,NAME5
      DOUBLE PRECISION RU(50),EU(50), RG(50),EG(50), 
     1     RPU(50),EPU(50), RPG(50),EPG(50), RX(50),EX(50),
     1     COEFU(4,50), COEFG(4,50), COEFPU(4,50), COEFPG(4,50), 
     1     COEFX(4,50)
      DATA NC/1/,NC1/1/,NC2/1/,NC3/1/,NC4/1/,NC5/1/, RAUA/0.52918D0/
      IF(NC.EQ.1) THEN
      NC=2
C     OPEN(110,FILE='RG2+_3+N1.F')
!     OPEN(110,FILE='ENERGY.F')
      IND1=2  
      IND2=2  
      IND3=2  
      IND4=2  
      IND5=2
!     READ(110,'(2X,5I3)') IND1,IND2,IND3,IND4,IND5
CC      PRINT '(2X,5I3)', IND1,IND2,IND3,IND4,IND5
!     READ(110,'(2X,5A12)') NAME1,NAME2,NAME3,NAME4,NAME5
CC      PRINT '(2X,5A12)', NAME1,NAME2,NAME3,NAME4,NAME5
!      NAME1='NE2+SU.DAT'  
!      NAME2='NE2+SG.DAT'
!      NAME3='NE2+PU.DAT'
!      NAME4='NE2+PG.DAT'
       NAME1='NE2+SU.DAT_Q'  
       NAME2='NE2+SG.DAT_Q'
       NAME3='NE2+PU.DAT_Q'
       NAME4='NE2+PG.DAT_Q'
       NAME5='NE2X.DAT'
                  END IF
      UT=0.0
C
      IF(I0.EQ.1) THEN
      IF(IND1.EQ.0.OR.IND1.EQ.2) THEN
      IF(NC1.EQ.1) THEN
      NC1=2
      OPEN(11,FILE=NAME1,STATUS='OLD')
      IR1=0
   10 CONTINUE
      IR1=IR1+1
      READ(11,*,END=11) RU(IR1), EU(IR1)
CC      PRINT *, RU(IR1), EU(IR1)
      GO TO 10
   11 CONTINUE
      IR1=IR1-1
      IF(IND1.EQ.2) THEN
      EASU=EU(IR1)
      DO I=1,IR1
      COEFU(1,I)=EU(I) -EASU
      ENDDO
      COEFU(2,IR1)=0.0
      CALL CUBSPL(RU,COEFU,IR1,0,1)
                    END IF
                   END IF
      IF(RR.LT.RU(IR1)) THEN
      IF(IND1.EQ.0) THEN
      DO 1 I=1,IR1-1
      IF(.NOT.( RU(I).LE.RR.AND.RR.LE.RU(I+1) )) GO TO 1
      UT = EU(I) + (EU(I+1)-EU(I))/(RU(I+1)-RU(I))*(RR-RU(I)) -EASU
    1 CONTINUE
                    END IF
      IF(IND1.EQ.2) THEN
      UT=PPVALU(RU,COEFU,IR1,4,RR,0)
        UT1=PPVALU(RU,COEFU,IR1,4,RR,1) 
                END IF
                        ELSE
      UT=EU(IR1) -EASU
      UT1=0.0
                        END IF
                    END IF
      IF(IND1.EQ.1) UT=UTN(RR/RAUA,I0,0)
      RETURN
                 END IF
C
      IF(I0.EQ.2) THEN
      IF(IND2.EQ.0.OR.IND2.EQ.2) THEN
      IF(NC2.EQ.1) THEN
      NC2=2
      OPEN(22,FILE=NAME2,STATUS='OLD')
      IR2=0
   20 CONTINUE
      IR2=IR2+1
      READ(22,*,END=22) RG(IR2), EG(IR2)
CC      PRINT *, RG(IR2), EG(IR2)
      GO TO 20
   22 CONTINUE
      IR2=IR2-1
      IF(IND2.EQ.2) THEN
      EASG=EG(IR2)
      DO I=1,IR2
      COEFG(1,I)=EG(I) -EASG
      ENDDO
      COEFG(2,IR2)=0.0
      CALL CUBSPL(RG,COEFG,IR2,0,1)
                    END IF
                   END IF
      IF(RR.LT.RG(IR2)) THEN
      IF(IND2.EQ.0) THEN
      DO 2 I=1,IR2-1
      IF(.NOT.( RG(I).LE.RR.AND.RR.LE.RG(I+1) )) GO TO 2
      UT = EG(I) + (EG(I+1)-EG(I))/(RG(I+1)-RG(I))*(RR-RG(I)) -EASG
    2 CONTINUE
                    END IF
      IF(IND2.EQ.2) THEN
      UT=PPVALU(RG,COEFG,IR2,4,RR,0)
        UT1=PPVALU(RG,COEFG,IR2,4,RR,1) 
                END IF
                        ELSE
      UT=EG(IR2) -EASG
      UT1=0.0
                        END IF
                    END IF 
      IF(IND2.EQ.1) UT=UTN(RR/RAUA,I0,0)
      RETURN
                 END IF
C
      IF(I0.EQ.3) THEN
      IF(IND3.EQ.0.OR.IND3.EQ.2) THEN
      IF(NC3.EQ.1) THEN
      NC3=2
      OPEN(33,FILE=NAME3,STATUS='OLD')
      IR3=0
   30 CONTINUE
      IR3=IR3+1
      READ(33,*,END=33) RPU(IR3), EPU(IR3)
CC      PRINT *, RPU(IR3), EPU(IR3)
      GO TO 30
   33 CONTINUE
      IR3=IR3-1
      IF(IND3.EQ.2) THEN
      EASPU=EPU(IR3)
      DO I=1,IR3
      COEFPU(1,I)=EPU(I) -EASPU
      ENDDO
      COEFPU(2,IR3)=0.0
      CALL CUBSPL(RPU,COEFPU,IR3,0,1)
                    END IF
                   END IF
      IF(RR.LT.RPU(IR3)) THEN
      IF(IND3.EQ.0) THEN  
      DO 3 I=1,IR3-1
      IF(.NOT.( RPU(I).LE.RR.AND.RR.LE.RPU(I+1) )) GO TO 3
      UT = EPU(I) + (EPU(I+1)-EPU(I))/(RPU(I+1)-RPU(I))*(RR-RPU(I))
     1     -EASPU
    3 CONTINUE
                    END IF
      IF(IND3.EQ.2) THEN
      UT=PPVALU(RPU,COEFPU,IR3,4,RR,0) 
        UT1=PPVALU(RPU,COEFPU,IR3,4,RR,1) 
                END IF
                        ELSE
      UT=EPU(IR3) -EASPU
      UT1=0.0
                        END IF
                    END IF
      IF(IND3.EQ.1) UT=UTN(RR/RAUA,I0,0)
      RETURN
               END IF
C
      IF(I0.EQ.4) THEN
      IF(IND4.EQ.0.OR.IND4.EQ.2) THEN
      IF(NC4.EQ.1) THEN
      NC4=2
      OPEN(44,FILE=NAME4,STATUS='OLD')
      IR4=0
   40 CONTINUE
      IR4=IR4+1
      READ(44,*,END=44) RPG(IR4), EPG(IR4)
CC      PRINT *, RPG(IR4), EPG(IR4)
      GO TO 40
   44 CONTINUE
      IR4=IR4-1
      IF(IND4.EQ.2) THEN
      EASPG=EPG(IR4)
      DO I=1,IR4
      COEFPG(1,I)=EPG(I) -EASPG
      ENDDO
      COEFPG(2,IR4)=0.0
      CALL CUBSPL(RPG,COEFPG,IR4,0,1)
                    END IF
                   END IF
      IF(RR.LT.RPG(IR4)) THEN
      IF(IND4.EQ.0) THEN
      DO 4 I=1,IR4-1
      IF(.NOT.( RPG(I).LE.RR.AND.RR.LE.RPG(I+1) )) GO TO 4
      UT = EPG(I) + (EPG(I+1)-EPG(I))/(RPG(I+1)-RPG(I))*(RR-RPG(I))
     1     -EASPG
    4 CONTINUE
                    END IF
      IF(IND4.EQ.2) THEN
      UT=PPVALU(RPG,COEFPG,IR4,4,RR,0)
        UT1=PPVALU(RPG,COEFPG,IR4,4,RR,1)
                END IF
                        ELSE
      UT=EPG(IR4) -EASPG
      UT1=0.0
                        END IF
                    END IF 
      IF(IND4.EQ.1) UT=UTN(RR/RAUA,I0,0)
      RETURN
                 END IF
C
      IF(I0.EQ.5) THEN
      IF(IND5.EQ.0.OR.IND5.EQ.2) THEN
      IF(NC5.EQ.1) THEN
      NC5=2
      OPEN(55,FILE=NAME5,STATUS='OLD')
      IR5=0
   50 CONTINUE
      IR5=IR5+1
      READ(55,*,END=55) RX(IR5), EX(IR5)
CC      PRINT *, RX(IR5), EX(IR5)
      GO TO 50
   55 CONTINUE
      IR5=IR5-1
      IF(IND5.EQ.2) THEN
      EASX=EX(IR5)
      DO I=1,IR5
      COEFX(1,I)=EX(I) -EASX
      ENDDO
      COEFX(2,IR5)=0.0
      CALL CUBSPL(RX,COEFX,IR5,0,1)
                    END IF
                   END IF
      IF(RR.LT.RX(IR5)) THEN
      IF(IND5.EQ.0) THEN
      DO 5 I=1,IR5-1
      IF(.NOT.( RX(I).LE.RR.AND.RR.LE.RX(I+1) )) GO TO 5
      UT = EX(I) + (EX(I+1)-EX(I))/(RX(I+1)-RX(I))*(RR-RX(I)) -EASX
    5 CONTINUE
                    END IF
      IF(IND5.EQ.2) THEN
      UT=PPVALU(RX,COEFX,IR5,4,RR,0)
      UT1=PPVALU(RX,COEFX,IR5,4,RR,1)
                END IF
                        ELSE 
      UT=EX(IR5) -EASX
      UT1=0.0
                        END IF
                    END IF
      IF(IND5.EQ.1) UT=UTN(RR/RAUA,I0,0)
      RETURN
                 END IF
      RETURN 
      END


      SUBROUTINE RGNI(N,X,GRAD,EREAL,GRADT) !  ,H0,H1,EE,EV,W,NATOMS)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER INFO
        DOUBLE PRECISION TEMPA(9*N)

      PARAMETER (NEMAX=100)
        INTEGER IW(NEMAX,3),  IWORK(3*5*N), IFAIL(3*N)

        DOUBLE PRECISION X(3,N),H0(N*3,N*3)
     &        ,H1(N*3,N*3)
     &        ,GRAD(N*3),EE(N*3),EV(N*3,N*3),W(N*3,N*3),WORK(3*8*N)
      LOGICAL GRADT
      DATA IEV/1/, NCALL/1/, SHIFT/10.0D0/
     &  , ELOW/-1.D3/,EPS/1.D-6/,NSTEP/1000/,MESSAGE/6/,NWORK/777/
     &  , IFLAG/-2/ 
      SAVE

      N3=N*3
CC      IF(GRADT) IGRAD=1

        ID=1
        IT=0
C       IT=1
C       IT=2
      CALL HMAT(X,N,N3,H0,H1,H2,0)
      IF(NCALL.EQ.1.OR.IT.EQ.0) THEN
         IF(ID.EQ.0) THEN
C            CALL JACOB2(H0,EV,W(1,1),EE,W(1,2),N3,IEV)
                     ELSE
C       CALL TRED2(H0,N3,N3,EE,W(1,1),IEV)
C       CALL TQLI(EE,W(1,1),N3,N3,H0,IEV)

        CALL DSYEV('V','U',N3,H0,N3,EE,TEMPA,9*N,INFO)
        IF (INFO.NE.0) THEN
           PRINT*,'WARNING - INFO=',INFO,' IN DSYEV'
        ENDIF

C       ABSTOL=0.0D0
C
C  THIS IS LOOP IS ONLY NEEDED FOR DSYEVX WHEN ONLY ONE EIGENVALUE IS FOUND.
C  THE EV=H0 LOOP SHOULD BE COMMENTED IS DSYEVX IS USED.
C
C       DO J1=1,N3
C          EE(J1)=1.0D6
C       ENDDO
C       CALL DSYEVX('V','I','U',N3,H0,N3,DUMMY1,DUMMY2,1,1,ABSTOL,M,EE,EV,N3,WORK,3*8*N,IWORK,IFAIL,INFO)

        DO I=1,N3
        DO J=1,N3
        EV(I,J)=H0(I,J)
        ENDDO
        ENDDO
                     END IF
        CALL SORTV(EE,EV,N3,IEV)
        EREAL=EE(N3)

C      IF(IT.EQ.2) OPEN(NWORK,FILE='ITLAN.TMP')

      NCALL=2
                                ELSE
      IF(IT.EQ.1) THEN
      DO I=1,N3
      H0(I,I)=H0(I,I)-SHIFT
      ENDDO
C        CALL ITEIG(H0,EV(1,N3),EREAL,W(1,1),W(1,2),N3)
      EREAL=EREAL+SHIFT
                END IF

      IF(IT.EQ.2) THEN

        EHIGH=DMAX1(0.D0,2.D0*EREAL)
      DO ISTEP=1,NSTEP
C      CALL ITLANE(N3,ELOW,EHIGH,EPS,N,N3,NSTEP,MESSAGE,NWORK
C     &      ,IFLAG,EV,EV(1,N3),EE,IW,NE,W,W(1,2),IW(1,3),H1,H1(1,2))
      IF(IFLAG.EQ.0) GO TO 1
      IF(IFLAG.EQ.1) THEN
      DO I=1,N3
      HVI=0.D0
      DO J=1,N3
      HVI=HVI+H0(I,J)*EV(J,N3)
      ENDDO
      EV(I,1)=EV(I,1) +HVI
      ENDDO
                          END IF
        IF(IFLAG.EQ.2) THEN
        PRINT *,' ITLANE: IFLAG=', IFLAG
        STOP
                       END IF
      ENDDO
  1     EREAL=EE(1)
      IF(IEV.NE.0) THEN
      NW=0
      PRINT '(A)','ERROR, NE HAS NOT BEEN SET!'
      STOP
      DO I=1,NE
         NW=MAX0(NW,IW(1,I))
      ENDDO
C        CALL ITLANV(N3,NSTEP,MESSAGE,NWORK
C     &      ,EE,IW,1,W,W(1,2),N3,NW,JFLAG,EV(1,N3),H1,H1(1,2))
                 END IF
C     IF(JFLAG.GT.0) THEN
C        PRINT *,' ITLANV: JFLAG= ',JFLAG
C        STOP
C     END IF

                END IF

                                END IF

      IF(GRADT) THEN

      DO I=1,N3
      CALL HMAT(X,N,N3,H0,H1,H2,I)
      VHVI=0.0D0
      DO J=1,N3
      VH=0.0D0
      DO K=1,N3
CC      VH=VH+EV(K,N3)*H1(K,J,I)
      VH=VH+EV(K,N3)*H1(K,J)
      ENDDO
      VHVI=VHVI+VH*EV(J,N3)
      ENDDO
      GRAD(I)=VHVI
      ENDDO

              END IF
      RETURN
      END
C ------------------------------------------------------------
C      DOUBLE PRECISION FUNCTION ENERGY(R0,I,IP)
C      DOUBLE PRECISION R0
C      R=R0
C      ENERGY=UT(R,I,IP)
C      RETURN
C      END

              SUBROUTINE SORTV(EE,EV,N,IEV)
      DOUBLE PRECISION EE(N),EV(N,N),EM,EVNN
      DO 17 M=1,N-1
      NN=N+1-M
      EM=EE(1)
      IND=1
      DO 15 K=2,NN
      IF(EE(K).GE.EM) GO TO 15
      EM=EE(K)
      IND=K
   15 CONTINUE
      IF(IND.EQ.NN) GO TO 17
      EE(IND)=EE(NN)
      EE(NN)=EM
      IF(IEV.GT.0) THEN
      DO 16 K=1,N
      EVNN=EV(K,NN)
      EV(K,NN)=EV(K,IND)
      EV(K,IND)=EVNN
   16 CONTINUE
                   END IF
   17 CONTINUE
      RETURN
      END

      SUBROUTINE HMAT(X,N,N3,H0,H1,H2,IGRAD)
      USE COMMONS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C       COMMON/TEST/ITEST

      DOUBLE PRECISION X(3,N),H0(N*3,N*3)
     &               ,H1(N*3,N*3)

      SAVE
      DATA RAUA/0.52918D0/
CC     &      ,ISU/1/,ISG/2/,IPU/3/,IPG/4/,IX/5/

      IF(IGRAD.NE.0) GO TO 10

      DO I=1,N3
      DO J=1,N3
      H0(I,J)=0.0D0
      ENDDO
      ENDDO

      EXS=0.0D0
      DO K=1,N-1
      DO L=K+1,N
      RKL=
     &      SQRT((X(1,L)-X(1,K))**2+(X(2,L)-X(2,K))**2+(X(3,L)-X(3,K))**2)
CC      EX=ENERGY(RKL,IX,0)
        IF (NEON) THEN
           EX=UTN(RKL/RAUA,5,0)
        ELSE
           CALL ASAR1(RKL,EX,EX1)
        ENDIF
      EXS=EXS+EX
      ENDDO
      ENDDO

      DO I=1,N3-1
      H0(I,I)=H0(I,I)+EXS
      DO J=I+1,N3

      KK=INT(I/3)
      IF(MOD(I,3).GT.0) KK=KK+1
      LL=INT(J/3)
      IF(MOD(J,3).GT.0) LL=LL+1
      IF(LL.EQ.KK) GO TO 1

      DX=X(1,LL)-X(1,KK)
      DY=X(2,LL)-X(2,KK)
      DZ=X(3,LL)-X(3,KK)
      DXY=SQRT(DX*DX+DY*DY)
      RKL=SQRT(DXY*DXY+DZ*DZ)
      CT=DZ/RKL
      ST=DXY/RKL
      CF=1.0D0
      SF=0.0D0
      IF(DXY.GT.0.0D0) THEN
      CF=DX/DXY
      SF=DY/DXY
                   END IF
        IF (NEON) THEN
           CALL FENERGY(RKL,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ELSE
         CALL ENERGYAR(RKL,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ENDIF
CC      ESU=ENERGY(RKL,ISU,0)
CC      ESG=ENERGY(RKL,ISG,0)
CC      EPU=ENERGY(RKL,IPU,0)
CC      EPG=ENERGY(RKL,IPG,0)
      ESS=0.5D0*(ESU+ESG)
      ESD=0.5D0*(ESU-ESG)
      EPS=0.5D0*(EPU+EPG)
      EPD=0.5D0*(EPU-EPG)
CC      EX=ENERGY(RKL,IX,0)
      
      IF(MOD(I,3).EQ.1) THEN
      IF(MOD(J,3).EQ.1) THEN
         TRIG1=(ST*CF)**2
         TRIG2=(CT*CF)**2+SF**2
CC      HS=ESS*(ST*CF)**2+EPS*((CT*CF)**2+SF**2)
      HS=ESS*TRIG1+EPS*TRIG2
      H0(I,I)=H0(I,I)+HS-EX
      H0(J,J)=H0(J,J)+HS-EX
CC      HD=ESD*(ST*CF)**2+EPD*((CT*CF)**2+SF**2)      
      HD=ESD*TRIG1+EPD*TRIG2
      H0(I,J)=HD
      H0(J,I)=HD
                    END IF
      IF(MOD(J,3).EQ.2) THEN
         TRIG=ST**2*SF*CF
CC        HS=(ESS-EPS)*ST**2*SF*CF
        HS=(ESS-EPS)*TRIG
        H0(I,I+1)=H0(I,I+1)+HS
        H0(I+1,I)=H0(I+1,I)+HS
CC      HD=(ESD-EPD)*ST**2*SF*CF
      HD=(ESD-EPD)*TRIG
      H0(I,J)=HD
      H0(J,I)=HD
                    END IF
      IF(MOD(J,3).EQ.0) THEN
         TRIG=ST*CT*CF
CC        HS=(ESS-EPS)*ST*CT*CF
        HS=(ESS-EPS)*TRIG
        H0(I,I+2)=H0(I,I+2)+HS
        H0(I+2,I)=H0(I+2,I)+HS
CC      HD=(ESD-EPD)*ST*CT*CF
      HD=(ESD-EPD)*TRIG
      H0(I,J)=HD
      H0(J,I)=HD
                    END IF
                    END IF

      IF(MOD(I,3).EQ.2) THEN
      IF(MOD(J,3).EQ.2) THEN
         TRIG1=(ST*SF)**2
         TRIG2=(CT*SF)**2+CF**2
CC      HS=ESS*(ST*SF)**2+EPS*((CT*SF)**2+CF**2)
      HS=ESS*TRIG1+EPS*TRIG2
      H0(I,I)=H0(I,I)+HS-EX
      H0(J,J)=H0(J,J)+HS-EX
CC      HD=ESD*(ST*SF)**2+EPD*((CT*SF)**2+CF**2)
      HD=ESD*TRIG1+EPD*TRIG2
      H0(I,J)=HD
      H0(J,I)=HD
                    END IF
      IF(MOD(J,3).EQ.1) THEN
         TRIG=ST**2*SF*CF
CC        HS=(ESS-EPS)*ST**2*SF*CF
        HS=(ESS-EPS)*TRIG
        H0(J,J+1)=H0(J,J+1)+HS
        H0(J+1,J)=H0(J+1,J)+HS
CC      HD=(ESD-EPD)*ST**2*SF*CF
      HD=(ESD-EPD)*TRIG
      H0(I,J)=HD
      H0(J,I)=HD
                    END IF
      IF(MOD(J,3).EQ.0) THEN
         TRIG=ST*CT*SF
CC        HS=(ESS-EPS)*ST*CT*SF
        HS=(ESS-EPS)*TRIG
        H0(I,I+1)=H0(I,I+1)+HS
        H0(I+1,I)=H0(I+1,I)+HS
CC      HD=(ESD-EPD)*ST*CT*SF
      HD=(ESD-EPD)*TRIG
      H0(I,J)=HD
      H0(J,I)=HD
                    END IF
                    END IF

      IF(MOD(I,3).EQ.0) THEN
      IF(MOD(J,3).EQ.0) THEN
         TRIG1=(CT)**2
         TRIG2=(ST)**2
CC      HS=ESS*(CT)**2   +EPS*(ST)**2
      HS=ESS*TRIG1   +EPS*TRIG2
      H0(I,I)=H0(I,I)+HS-EX
      H0(J,J)=H0(J,J)+HS-EX
CC      HD=ESD*(CT)**2   +EPD*(ST)**2
      HD=ESD*TRIG1   +EPD*TRIG2
      H0(I,J)=HD
      H0(J,I)=HD
                    END IF
      IF(MOD(J,3).EQ.1) THEN
         TRIG=ST*CT*CF
CC        HS=(ESS-EPS)*ST*CT*CF
        HS=(ESS-EPS)*TRIG
        H0(J,J+2)=H0(J,J+2)+HS
        H0(J+2,J)=H0(J+2,J)+HS
CC      HD=(ESD-EPD)*ST*CT*CF
      HD=(ESD-EPD)*TRIG
      H0(I,J)=HD
      H0(J,I)=HD
                    END IF
      IF(MOD(J,3).EQ.2) THEN
         TRIG=ST*CT*SF
CC        HS=(ESS-EPS)*ST*CT*SF
        HS=(ESS-EPS)*TRIG
        H0(J,J+1)=H0(J,J+1)+HS
        H0(J+1,J)=H0(J+1,J)+HS
CC      HD=(ESD-EPD)*ST*CT*SF
      HD=(ESD-EPD)*TRIG
      H0(I,J)=HD
      H0(J,I)=HD
                    END IF
                    END IF

   1      CONTINUE
      ENDDO
      ENDDO
      H0(N3,N3)=H0(N3,N3)+EXS

      IF(DIPOLE) CALL DIPOLES(X,N,H0,H1,IGRAD)

C ---------------------------------------------------------------------

  10      CONTINUE
      IF(IGRAD.LE.0) GO TO 2

      DO I=1,N3
      DO J=1,N3
CC      DO K=1,N3
CC      H1(I,J,K)=0.0D0
CC      ENDDO
      H1(I,J)=0.0D0
      ENDDO
      ENDDO

CC      DO M=1,N
CC      DO IR=1,3
CC      IC=(M-1)*3+IR

      M=INT(IGRAD/3)+1
      IR=MOD(IGRAD,3)
      IF(IR.EQ.0) THEN
      IR=3
      M=M-1
                  END IF

      EXS1=0.0D0
      DO K=1,N-1
      DO L=K+1,N
      IF(L.EQ.M.OR.K.EQ.M) THEN
      RKL=
     &      SQRT((X(1,L)-X(1,K))**2+(X(2,L)-X(2,K))**2+(X(3,L)-X(3,K))**2)
CC      IF(L.EQ.M) EX1= ENERGY(RKL,IX,1)*(X(IR,L)-X(IR,K))/RKL
CC      IF(K.EQ.M) EX1=-ENERGY(RKL,IX,1)*(X(IR,L)-X(IR,K))/RKL
        IF (NEON) THEN
      IF(L.EQ.M) EX1= UTN(RKL/RAUA,5,1)/RAUA*(X(IR,L)-X(IR,K))/RKL
      IF(K.EQ.M) EX1=-UTN(RKL/RAUA,5,1)/RAUA*(X(IR,L)-X(IR,K))/RKL
        ELSE
           CALL ASAR1(RKL,EX,EX1)
        IF(L.EQ.M) EX1= EX1*(X(IR,L)-X(IR,K))/RKL
        IF(K.EQ.M) EX1=-EX1*(X(IR,L)-X(IR,K))/RKL
        ENDIF
      EXS1=EXS1+EX1
                       END IF
      ENDDO
      ENDDO

      DO I=1,N3-1
      H1(I,I)=H1(I,I)+EXS1
      DO J=I+1,N3

      KK=INT(I/3)
      IF(MOD(I,3).GT.0) KK=KK+1
      LL=INT(J/3)
      IF(MOD(J,3).GT.0) LL=LL+1
      IF(LL.EQ.KK) GO TO 11

      IF(LL.EQ.M.OR.KK.EQ.M) THEN

      DX=X(1,LL)-X(1,KK)
      DY=X(2,LL)-X(2,KK)
      DZ=X(3,LL)-X(3,KK)
      DXY=SQRT(DX*DX+DY*DY)
      RKL=SQRT(DXY*DXY+DZ*DZ)
      CT=DZ/RKL
      ST=DXY/RKL
      CF=1.0D0
      SF=0.0D0
      IF(DXY.GT.0.0D0) THEN
      CF=DX/DXY
      SF=DY/DXY
                   END IF

      IF(LL.EQ.M) DRDX= (X(IR,LL)-X(IR,KK))/RKL
      IF(KK.EQ.M) DRDX=-(X(IR,LL)-X(IR,KK))/RKL

      CT1=-CT/RKL*DRDX
      IF(IR.EQ.3) THEN
      IF(M.EQ.LL) CT1=CT1+1.0D0/RKL
      IF(M.EQ.KK) CT1=CT1-1.0D0/RKL
                END IF
      ST1=-ST/RKL*DRDX
      DXYDX=0.0D0
      IF(IR.EQ.1.OR.IR.EQ.2) THEN
C      DXYDX=1.0D0
      IF(DXY.GT.0.0D0) THEN
      IF(LL.EQ.M) DXYDX= (X(IR,LL)-X(IR,KK))/DXY
      IF(KK.EQ.M) DXYDX=-(X(IR,LL)-X(IR,KK))/DXY
                   END IF
      ST1=ST1+DXYDX/RKL
                             END IF
      CF1=0.0D0
      SF1=0.0D0
      IF(DXY.GT.0.0D0) THEN
      IF(IR.EQ.1.OR.IR.EQ.2) THEN
      CF1=-CF/DXY*DXYDX
      IF(IR.EQ.1) THEN
      IF(M.EQ.LL) CF1=CF1+1.0D0/DXY
      IF(M.EQ.KK) CF1=CF1-1.0D0/DXY
                    END IF
      SF1=-SF/DXY*DXYDX
      IF(IR.EQ.2) THEN
      IF(M.EQ.LL) SF1=SF1+1.0D0/DXY
      IF(M.EQ.KK) SF1=SF1-1.0D0/DXY
                    END IF
                             END IF
                   END IF

        IF (NEON) THEN
           CALL FENERGY(RKL,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ELSE
         CALL ENERGYAR(RKL,ESU,ESG,EPU,EPG,EX,ESU1,ESG1,EPU1,EPG1,EX1)
        ENDIF
CC      ESU=ENERGY(RKL,ISU,0)
CC      ESG=ENERGY(RKL,ISG,0)
CC      EPU=ENERGY(RKL,IPU,0)
CC      EPG=ENERGY(RKL,IPG,0)
      ESS=0.5D0*(ESU+ESG)
      ESD=0.5D0*(ESU-ESG)
      EPS=0.5D0*(EPU+EPG)
      EPD=0.5D0*(EPU-EPG)
CC      EX=ENERGY(RKL,IX,0)

CC      ESU1=ENERGY(RKL,ISU,1)*DRDX
CC      ESG1=ENERGY(RKL,ISG,1)*DRDX
CC      EPU1=ENERGY(RKL,IPU,1)*DRDX
CC      EPG1=ENERGY(RKL,IPG,1)*DRDX
      ESU1=ESU1*DRDX
      ESG1=ESG1*DRDX
      EPU1=EPU1*DRDX
      EPG1=EPG1*DRDX
      ESS1=0.5D0*(ESU1+ESG1)
      ESD1=0.5D0*(ESU1-ESG1)
      EPS1=0.5D0*(EPU1+EPG1)
      EPD1=0.5D0*(EPU1-EPG1)
CC      EX1=ENERGY(RKL,IX,1)*DRDX
      EX1=EX1*DRDX
      
      IF(MOD(I,3).EQ.1) THEN
      IF(MOD(J,3).EQ.1) THEN
         TRIG1=(ST*CF)**2
         TRIG2=(CT*CF)**2+SF**2
         TRIG3=2.0D0*ST*CF*(CF*ST1+ST*CF1)
         TRIG4=(2.0D0*CT*CF*(CF*CT1+CT*CF1)+2.0D0*SF*SF1)
C      HS=ESS*(ST*CF)**2+EPS*((CT*CF)**2+SF**2)
CC      HS1=ESS1*(ST*CF)**2+EPS1*((CT*CF)**2+SF**2) +
CC     &          ESS*2.0D0*ST*CF*(CF*ST1+ST*CF1)+
CC     &      EPS*(2.0D0*CT*CF*(CF*CT1+CT*CF1)+2.0D0*SF*SF1)
      HS1= ESS1*TRIG1+EPS1*TRIG2+ESS*TRIG3+EPS*TRIG4
      H1(I,I)=H1(I,I)+HS1-EX1
      H1(J,J)=H1(J,J)+HS1-EX1
C      HD=ESD*(ST*CF)**2+EPD*((CT*CF)**2+SF**2)      
CC      HD1=ESD1*(ST*CF)**2+EPD1*((CT*CF)**2+SF**2) +
CC     &          ESD*2.0D0*ST*CF*(CF*ST1+ST*CF1)+
CC     &      EPD*(2.0D0*CT*CF*(CF*CT1+CT*CF1)+2.0D0*SF*SF1)
      HD1= ESD1*TRIG1 + EPD1*TRIG2 + ESD*TRIG3 + EPD*TRIG4
      H1(I,J)=HD1
      H1(J,I)=HD1
                    END IF
      IF(MOD(J,3).EQ.2) THEN
         TRIG1=ST**2*SF*CF
         TRIG2=2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1)
C        HS=(ESS-EPS)*ST**2*SF*CF
CC        HS1=(ESS1-EPS1)*ST**2*SF*CF +
CC     &      (ESS-EPS)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1)) 
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(I,I+1)=H1(I,I+1)+HS1
        H1(I+1,I)=H1(I+1,I)+HS1
C      HD=(ESD-EPD)*ST**2*SF*CF
CC      HD1=(ESD1-EPD1)*ST**2*SF*CF + 
CC     &      (ESD-EPD)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1))
      HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(I,J)=HD1
      H1(J,I)=HD1
                    END IF
      IF(MOD(J,3).EQ.0) THEN
         TRIG1=ST*CT*CF
         TRIG2=ST1*CT*CF+ST*CT1*CF+ST*CT*CF1
C        HS=(ESS-EPS)*ST*CT*CF
CC        HS1=(ESS1-EPS1)*ST*CT*CF + 
CC     &      (ESS-EPS)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(I,I+2)=H1(I,I+2)+HS1
        H1(I+2,I)=H1(I+2,I)+HS1
C      HD=(ESD-EPD)*ST*CT*CF
CC      HD1=(ESD1-EPD1)*ST*CT*CF + 
CC     &      (ESD-EPD)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
      HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(I,J)=HD1
      H1(J,I)=HD1
                    END IF
                    END IF

      IF(MOD(I,3).EQ.2) THEN
      IF(MOD(J,3).EQ.2) THEN
         TRIG1=(ST*SF)**2
         TRIG2=(CT*SF)**2+CF**2
         TRIG3=2.0D0*ST*SF*(SF*ST1+ST*SF1)
         TRIG4=2.0D0*CT*SF*(SF*CT1+CT*SF1)+2.0D0*CF*CF1
C      HS=ESS*(ST*SF)**2+EPS*((CT*SF)**2+CF**2)
CC      HS1=ESS1*(ST*SF)**2+EPS1*((CT*SF)**2+CF**2) +
CC     &          ESS*2.0D0*ST*SF*(SF*ST1+ST*SF1)+
CC     &      EPS*(2.0D0*CT*SF*(SF*CT1+CT*SF1)+2.0D0*CF*CF1)
      HS1=ESS1*TRIG1+EPS1*TRIG2+ESS*TRIG3+EPS*TRIG4
      H1(I,I)=H1(I,I)+HS1-EX1
      H1(J,J)=H1(J,J)+HS1-EX1
C      HD=ESD*(ST*SF)**2+EPD*((CT*SF)**2+CF**2)
CC      HD1=ESD1*(ST*SF)**2+EPD1*((CT*SF)**2+CF**2) +
CC     &          ESD*2.0D0*ST*SF*(SF*ST1+ST*SF1)+
CC     &      EPD*(2.0D0*CT*SF*(SF*CT1+CT*SF1)+2.0D0*CF*CF1)
      HD1=ESD1*TRIG1+EPD1*TRIG2+ESD*TRIG3+EPD*TRIG4
      H1(I,J)=HD1
      H1(J,I)=HD1
                    END IF
      IF(MOD(J,3).EQ.1) THEN
         TRIG1=ST**2*SF*CF
         TRIG2=2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1)
C        HS=(ESS-EPS)*ST**2*SF*CF
CC        HS1=(ESS1-EPS1)*ST**2*SF*CF +
CC     &      (ESS-EPS)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1))
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(J,J+1)=H1(J,J+1)+HS1
        H1(J+1,J)=H1(J+1,J)+HS1
C      HD=(ESD-EPD)*ST**2*SF*CF
CC      HD1=(ESD1-EPD1)*ST**2*SF*CF + 
CC     &      (ESD-EPD)*(2.0D0*ST*ST1*SF*CF+ST**2*(SF1*CF+SF*CF1))
      HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(I,J)=HD1
      H1(J,I)=HD1
                    END IF
      IF(MOD(J,3).EQ.0) THEN
         TRIG1=ST*CT*SF
         TRIG2=ST1*CT*SF+ST*CT1*SF+ST*CT*SF1
C        HS=(ESS-EPS)*ST*CT*SF
CC        HS1=(ESS1-EPS1)*ST*CT*SF +
CC     &      (ESS-EPS)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(I,I+1)=H1(I,I+1)+HS1
        H1(I+1,I)=H1(I+1,I)+HS1
C      HD=(ESD-EPD)*ST*CT*SF
CC      HD1=(ESD1-EPD1)*ST*CT*SF + 
CC     &      (ESD-EPD)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
      HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(I,J)=HD1
      H1(J,I)=HD1
                    END IF
                    END IF

      IF(MOD(I,3).EQ.0) THEN
      IF(MOD(J,3).EQ.0) THEN
         TRIG1=(CT)**2
         TRIG2=(ST)**2
         TRIG3=2.0D0*CT*CT1
         TRIG4=2.0D0*ST*ST1
C      HS=ESS*(CT)**2   +EPS*(ST)**2
CC      HS1=ESS1*(CT)**2   +EPS1*(ST)**2 +
CC     &      ESS*2.0D0*CT*CT1   +EPS*2.0D0*ST*ST1
      HS1=ESS1*TRIG1 + EPS1*TRIG2 + ESS*TRIG3 + EPS*TRIG4
      H1(I,I)=H1(I,I)+HS1-EX1
      H1(J,J)=H1(J,J)+HS1-EX1
C      HD=ESD*(CT)**2   +EPD*(ST)**2
CC      HD1=ESD1*(CT)**2   +EPD1*(ST)**2 +
CC     &      ESD*2.0D0*CT*CT1   +EPD*2.0D0*ST*ST1
      HD1=ESD1*TRIG1 + EPD1*TRIG2 + ESD*TRIG3 + EPD*TRIG4
      H1(I,J)=HD1
      H1(J,I)=HD1
                    END IF
      IF(MOD(J,3).EQ.1) THEN
         TRIG1=ST*CT*CF
         TRIG2=ST1*CT*CF+ST*CT1*CF+ST*CT*CF1
C        HS=(ESS-EPS)*ST*CT*CF
CC        HS1=(ESS1-EPS1)*ST*CT*CF + 
CC     &      (ESS-EPS)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(J,J+2)=H1(J,J+2)+HS1
        H1(J+2,J)=H1(J+2,J)+HS1
C      HD=(ESD-EPD)*ST*CT*CF
CC      HD1=(ESD1-EPD1)*ST*CT*CF + 
CC     &      (ESD-EPD)*(ST1*CT*CF+ST*CT1*CF+ST*CT*CF1)
      HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(I,J)=HD1
      H1(J,I)=HD1
                    END IF
      IF(MOD(J,3).EQ.2) THEN
         TRIG1=ST*CT*SF
         TRIG2=ST1*CT*SF+ST*CT1*SF+ST*CT*SF1
C        HS=(ESS-EPS)*ST*CT*SF
CC        HS1=(ESS1-EPS1)*ST*CT*SF + 
CC     &      (ESS-EPS)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
        HS1=(ESS1-EPS1)*TRIG1 + (ESS-EPS)*TRIG2
        H1(J,J+1)=H1(J,J+1)+HS1
        H1(J+1,J)=H1(J+1,J)+HS1
C      HD=(ESD-EPD)*ST*CT*SF
CC      HD1=(ESD1-EPD1)*ST*CT*SF + 
CC     &      (ESD-EPD)*(ST1*CT*SF+ST*CT1*SF+ST*CT*SF1)
      HD1=(ESD1-EPD1)*TRIG1 + (ESD-EPD)*TRIG2
      H1(I,J)=HD1
      H1(J,I)=HD1
                    END IF
                    END IF

                       END IF

  11      CONTINUE
      ENDDO
      ENDDO
      H1(N3,N3)=H1(N3,N3)+EXS1

CC      ENDDO
CC      ENDDO

        IF(DIPOLE) CALL DIPOLES(X,N,H0,H1,IGRAD)

   2      CONTINUE

      RETURN
      END

CC           REAL*4 FUNCTION UTN(R,I,IP)
CC           UTN=0.0D0
CC           IF(I.NE.5) RETURN
CC           UTN=ASNE1(R,I,IP)
CC           RETURN
CC           END
CC                 FUNCTION ASNE1(R,I,IP)
                 DOUBLE PRECISION FUNCTION UTN(R,I,IP)
           IMPLICIT DOUBLE PRECISION (A-H,O-Z)
           SAVE
             COMMON/DIPOL/ ADIP,RUNIT,REX                                            
           DATA E/1.34D-4/,RM/5.841D0/, A/13.86434671D0/,B/-.12993822D0/, 
     & D/1.36D0/,
     * C6,C8,C10/1.21317545D0,.53222749D0,.24570703D0/, AA/8.9571795D5/
      REX=RM
           X=R/RM
           F=1.D0
           IF(X.LT.D) F=EXP(-(D/X-1.D0)**2)

           X2=X*X
           X4=X2*X2
           X6=X4*X2
           X8=X4*X4
           X10=X6*X4
           IF(IP.EQ.0) THEN
           V=AA*EXP((-A+B*X)*X)-F*(C6/X6+C8/X8+C10/X10)
CC           ASNE1=V*E
           UTN=V*E
           RETURN
                       END IF
           IF(IP.EQ.1) THEN
           F1=0.D0
           IF(X.LT.D) F1=F*2.D0*(D/X-1.D0)*D/X2
           V=AA*EXP((-A+B*X)*X) *(-A+2.D0*B*X)
     * +F*(6.D0*C6/X6+8.D0*C8/X8+10.D0*C10/X10)/X
           IF(X.LT.D) V=V-F1*(C6/X6+C8/X8+C10/X10)
CC           ASNE1=V*E/RM
           UTN=V*E/RM
                       END IF
           RETURN
           END

        SUBROUTINE RGNX2(N,X,GRAD,EREAL,GRADT)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DOUBLE PRECISION X(3,N),GRAD(N*3)
        LOGICAL GRADT
        SAVE

      DATA HRX2/0.995D0/, EAUEV/27.212D0/, RAUA/0.52918D0/

      EREAL=0.D0
        IF(N.LE.0) RETURN
      DO I=1,N*3
      GRAD(I)=0.D0
      ENDDO

        IF(N.GT.1) THEN
        DO K=1,N-1
        DO L=K+1,N
        RKL=
     1  SQRT((X(1,L)-X(1,K))**2+(X(2,L)-X(2,K))**2+(X(3,L)-X(3,K))**2)
        CALL ASHE1(RKL/RAUA,ERR,ERR1)
        EREAL=EREAL+ERR*EAUEV

      IF(GRADT) THEN
      DO M=1,3
      DEDXM=ERR1*EAUEV/RAUA*(X(M,L)-X(M,K))/RKL
      LGRAD=3*(L-1)+M
      GRAD(LGRAD)=GRAD(LGRAD)+DEDXM
        KGRAD=3*(K-1)+M
        GRAD(KGRAD)=GRAD(KGRAD)-DEDXM
      ENDDO
              END IF
        ENDDO
        ENDDO
                   END IF

      DO I=1,N
      D2XY=X(1,I)**2+X(2,I)**2
      RXA=SQRT(D2XY + (X(3,I)-HRX2)**2)
        RXB=SQRT(D2XY + (X(3,I)+HRX2)**2)                
      CXA=(X(3,I)-HRX2)/RXA
      CXB=(X(3,I)+HRX2)/RXB
      CALL EHECL2(RXA,CXA,ERXA,ERXA1,ERXA1C)
        CALL EHECL2(RXB,-CXB,ERXB,ERXB1,ERXB1C)
      EREAL=EREAL+ERXA+ERXB

      IF(GRADT) THEN
      DO M=1,3
      IF(M.LE.2) THEN
      DEADXM=(ERXA1 - ERXA1C*CXA/RXA)*X(M,I)/RXA
        DEBDXM=(ERXB1 + ERXB1C*CXB/RXB)*X(M,I)/RXB
               ELSE
        DEADXM=(ERXA1 - ERXA1C*CXA/RXA)*(X(M,I)-HRX2)/RXA + ERXA1C/RXA
        DEBDXM=(ERXB1 + ERXB1C*CXB/RXB)*(X(M,I)+HRX2)/RXB - ERXB1C/RXB
               END IF
      IGRAD=3*(I-1)+M
        GRAD(IGRAD)=GRAD(IGRAD)+DEADXM+DEBDXM
      ENDDO
              END IF
      ENDDO

      RETURN
      END


      SUBROUTINE EHECL2(R,C,E,E1,E1C)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/INIT/ES,EP,EPB
      DATA HRE/0.995D0/, NC/1/
      SAVE

      ES = HECL2UT0(R,1,ES1)
      EP = HECL2UT0(R,2,EP1)
        EPB= HECL2UT0(R,3,EPB1)

        RATIO=HRE/R
        A= 1.D0 + RATIO*RATIO
        B= 2.*RATIO
        A1= -2.D0*(A-1.D0)/R
        B1=-B/R

        C2= C*C
      S2= 1.D0 - C2
      ABC= A + B*C
      S02= S2/ABC
      S021= -S02/ABC*(A1 + B1*C)
      S021C= -(2.D0*C + S02*B)/ABC
      EPT=EPB
      EPT1=EPB1
      IF(S02.GT.0.5D0) THEN
      FACTOR= 4.D0*S02*(1.D0 - S02)
        FACTOR1= 4.D0*S021*(1.D0 - 2.D0*S02)
      FACTOR1C= 4.D0*S021C*(1.D0 - 2.D0*S02)
      EPT= EP + (EPB - EP)*FACTOR
      EPT1= EP1 + (EPB1 - EP1)*FACTOR + (EPB - EP)*FACTOR1
        EPT1C= (EPB - EP)*FACTOR1C 
                   END IF

      E = ES*C2 + EPT*S2 
      E1= ES1*C2 + EPT1*S2
      E1C= 2.D0*C*(ES - EPT) + EPT1C*S2

      RETURN
      END

      DOUBLE PRECISION FUNCTION HECL2UT0(RR,I0,UT1)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IR
      PARAMETER(IR=41)
      DOUBLE PRECISION, DIMENSION(41) :: RU=(/
     1      2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90,     
     1      3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90,
     1      4.00, 4.10, 4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90,
     1      5.00, 5.10, 5.20, 5.30, 5.40, 5.50, 5.60, 5.70, 5.80, 5.90,
     1      6.00 /)
      DOUBLE PRECISION, DIMENSION(41) :: EU=(/
     1       0.3609288,  0.2394266,  0.1526561,  0.0937677,  0.0559396,
     1        0.0322421,  0.0168441,  0.0070019,  0.0010727, -0.0022567,
     1      -0.0039128, -0.0046519, -0.0048347, -0.0046815, -0.0043490,
     1      -0.0039317, -0.0035008, -0.0030874, -0.0027041, -0.0023572,
     1      -0.0020523, -0.0017784, -0.0015446, -0.0013409, -0.0011662,
     1       -0.0010178, -0.0008888, -0.0007792, -0.0006841, -0.0006028,
     1      -0.0005304, -0.0004683, -0.0004153, -0.0003680, -0.0003261,
     1      -0.0002896, -0.0002588, -0.0002309, -0.0002080, -0.0001890,
     1      -0.0001690/)
      DOUBLE PRECISION, DIMENSION(41) :: EG=(/
     1       0.6226088,  0.4446164,  0.3092623,  0.2090408,  0.1372434,
     1       0.0879939,  0.0557593,  0.0348669,  0.0208967,  0.0116730,
     1       0.0057593,  0.0021143, -0.0000627, -0.0013357, -0.0020061,
     1      -0.0023021, -0.0023595, -0.0022811, -0.0021323, -0.0019467,
     1      -0.0017533, -0.0015611, -0.0013833, -0.0012201, -0.0010725,
     1      -0.0009412, -0.0008269, -0.0007237, -0.0006362, -0.0005610,
     1      -0.0004941, -0.0004352, -0.0003848, -0.0003403, -0.0003010,
     1      -0.0002675, -0.0002369, -0.0002114, -0.0001892, -0.0001696,
     1      -0.0001528/)
      DOUBLE PRECISION, DIMENSION(41) :: EPU=(/
     1       0.8823669,  0.6168400,  0.4294591,  0.2970549,  0.2044612,
     1       0.1390856,  0.0935656,  0.0623482,  0.0408527,  0.0261192,
     1       0.0161569,  0.0095486,  0.0051923,  0.0023638,  0.0005976,
     1      -0.0004888, -0.0011015, -0.0014176, -0.0015430, -0.0015577,
     1      -0.0015120, -0.0014252, -0.0013083, -0.0011809, -0.0010553,
     1      -0.0009371, -0.0008288, -0.0007346, -0.0006490, -0.0005784,
     1      -0.0005115, -0.0004544, -0.0004049, -0.0003628, -0.0003226,
     1      -0.0002862, -0.0002551, -0.0002274, -0.0001998, -0.0001736,
     1      -0.0001570/)
      DOUBLE PRECISION COEFU(4,IR), COEFG(4,IR), COEFPU(4,IR)
      DATA NC1/1/,NC2/1/,NC3/1/
      SAVE

      HECL2UT0=0.D0
      UT1=0.D0

      IF(I0.EQ.1) THEN
      IF(NC1.EQ.1) THEN
      NC1=2
      DO I=1,IR
      COEFU(1,I)=EU(I) 
      ENDDO
      COEFU(2,IR)=-6.D0*EU(IR)/RU(IR)
      CALL CUBSPL(RU,COEFU,IR,0,1)
        R0U=RU(1)
        AAU=EU(1)
        ALPHAU=-COEFU(2,1)/AAU
                   END IF
        IF(RR.LT.RU(IR)) THEN
        IF(RR.LT.R0U) THEN
        HECL2UT0=AAU*EXP(-ALPHAU*(RR-R0U))
      UT1=-ALPHAU*HECL2UT0
                        ELSE
      IND=INT(10.D0*(RR-R0U)) +1
      DR=RR-RU(IND)
      HECL2UT0=COEFU(1,IND)+DR*(COEFU(2,IND)
     1      +0.5D0*DR*(COEFU(3,IND)+DR*COEFU(4,IND)/3.D0))
        UT1=COEFU(2,IND) + DR*(COEFU(3,IND) + 0.5D0*DR*COEFU(4,IND))
                        END IF
                         ELSE
            HECL2UT0=EU(IR)*(RU(IR)/RR)**6
      UT1=-6.D0*HECL2UT0/RR
                         END IF
      RETURN
                  END IF

      IF(I0.EQ.2) THEN
      IF(NC2.EQ.1) THEN
      NC2=2
      DO I=1,IR
      COEFG(1,I)=EG(I)
      ENDDO
      COEFG(2,IR)=-6.D0*EG(IR)/RU(IR)
      CALL CUBSPL(RU,COEFG,IR,0,1)
        R0G=RU(1)
        AAG=EG(1)
        ALPHAG=-COEFG(2,1)/AAG
                   END IF
        IF(RR.LT.RU(IR)) THEN
        IF(RR.LT.R0G) THEN
        HECL2UT0=AAG*EXP(-ALPHAG*(RR-R0G))
        UT1=-ALPHAG*HECL2UT0
                        ELSE
        IND=INT(10.D0*(RR-R0G)) +1
        DR=RR-RU(IND)
        HECL2UT0=COEFG(1,IND)+DR*(COEFG(2,IND)
     1      +0.5D0*DR*(COEFG(3,IND)+DR*COEFG(4,IND)/3.D0))        
        UT1=COEFG(2,IND) + DR*(COEFG(3,IND) + 0.5D0*DR*COEFG(4,IND))
                        END IF
                         ELSE
        HECL2UT0=EG(IR)*(RU(IR)/RR)**6
        UT1=-6.D0*HECL2UT0/RR
                         END IF
      RETURN
                  END IF

      IF(I0.EQ.3) THEN
      IF(NC3.EQ.1) THEN
      NC3=2
      DO I=1,IR
      COEFPU(1,I)=EPU(I)
      ENDDO
      COEFPU(2,IR)=-6.D0*EPU(IR)/RU(IR)
      CALL CUBSPL(RU,COEFPU,IR,0,1)
        R0PU=RU(1)
        AAPU=EPU(1)
        ALPHAPU=-COEFPU(2,1)/AAPU
                   END IF
       IF(RR.LT.RU(IR)) THEN
       IF(RR.LT.R0PU) THEN
        HECL2UT0=AAPU*EXP(-ALPHAPU*(RR-R0PU))
        UT1=-ALPHAPU*HECL2UT0
                        ELSE
        IND=INT(10.D0*(RR-R0PU)) +1
        DR=RR-RU(IND)
        HECL2UT0=COEFPU(1,IND)+DR*(COEFPU(2,IND)
     1      +0.5D0*DR*(COEFPU(3,IND)+DR*COEFPU(4,IND)/3.D0))
        UT1=COEFPU(2,IND) + DR*(COEFPU(3,IND) + 0.5D0*DR*COEFPU(4,IND))
                        END IF
                         ELSE
        HECL2UT0=EPU(IR)*(RU(IR)/RR)**6
        UT1=-6.D0*HECL2UT0/RR                               
                         END IF
      RETURN
                  END IF
      RETURN 
      END

           SUBROUTINE ASHE1(R,EX,EX1)
C  POTENTIAL OF HE2 (AZIZ ET AL, 1987)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DATA E/10.948D0/,RM/2.963D0/, A/10.43329537D0/,B/-2.27965105D0/, D/1.4826D0/
     *    ,C6,C8,C10/1.36745214D0,0.42123807D0,0.17473318D0/, AA/1.8443101D5/
     *    ,NCALL/1/,EKAU/3.1668D-6/,RAU/0.52918D0/
      IF(NCALL.EQ.1) THEN
      RM=RM/RAU
      E=E*EKAU
      NCALL=2
                     END IF
      X=R/RM
        X2=X*X
      F=1.D0
      F1=0.D0
      IF(X.LT.D) THEN
      DX1=D/X-1.D0
      F=EXP(-DX1*DX1)
      F1=F*2.D0*DX1*D/X2
             END IF
      X6=X2*X2*X2
      VDW=(C6+(C8+C10/X2)/X2)/X6
      AEXABX=AA*EXP((-A+B*X)*X)
      V=AEXABX - F*VDW
      V1=(-A+2.D0*B*X)*AEXABX - F1*VDW + 
     1            F*(6.D0*C6+(8.D0*C8+10.D0*C10/X2)/X2)/X6/X
      EX=V*E
      EX1=V1*E/RM
      RETURN
      END
        SUBROUTINE RGNXY(N,X,GRAD,EREAL,GRADT,SOCOUPLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        SAVE
        DOUBLE PRECISION X(3,N),GRAD(N*3)
        LOGICAL GRADT,SOCOUPLE

      DATA EAUEV/27.212D0/, RAUA/0.52918D0/
     1      HSPLIT/0.0074267D0/

      EREAL=0.D0
        IF(N.LE.0) RETURN
      DO I=1,N*3
      GRAD(I)=0.D0
      ENDDO

        IF(N.GT.1) THEN
        DO K=1,N-1
        DO L=K+1,N
        RKL=
     1  SQRT((X(1,L)-X(1,K))**2+(X(2,L)-X(2,K))**2+(X(3,L)-X(3,K))**2)
        CALL ASAR1(RKL/RAUA,ERR,ERR1)
        EREAL=EREAL+ERR*EAUEV

      IF(GRADT) THEN
      DO M=1,3
      DEDXM=ERR1*EAUEV/RAUA*(X(M,L)-X(M,K))/RKL
      LGRAD=3*(L-1)+M
      GRAD(LGRAD)=GRAD(LGRAD)+DEDXM
        KGRAD=3*(K-1)+M
        GRAD(KGRAD)=GRAD(KGRAD)-DEDXM
      ENDDO
              END IF
        ENDDO
        ENDDO
                   END IF

      E1R=0.D0
      E2R=0.D0
      E12R=0.D0
      DO I=1,N

      X1I=X(1,I)
      X2I=X(2,I)
      X3I=X(3,I)
      DXY2=X1I*X1I+X2I*X2I
      RI=SQRT(DXY2+X3I*X3I)
      CI=X3I/RI
      CALL ARNOENERGY(RI,CI,E1RI,E1RI1,E1RI1C,E2RI,E2RI1,E2RI1C)
      DXY=SQRT(DXY2)
      CAI=1.D0
      SAI=0.D0
      IF(DXY.GT.0.D0) THEN
      CAI=X1I/DXY
      SAI=X2I/DXY
                  END IF
      DERI=E1RI-E2RI
      CAI2=CAI*CAI
      E1R=E1R +DERI*CAI2 +E2RI
      E2R=E2R -DERI*CAI2 +E1RI
      E12R=E12R +CAI*SAI*DERI

        ENDDO

        E12R2=E12R*E12R
      IF(SOCOUPLE) E12R2=E12R2 + HSPLIT*HSPLIT

      D = SQRT((E1R-E2R)**2 + 4.D0*E12R2)
      EREAL= EREAL + 0.5D0*(E1R+E2R + D)
      IF(SOCOUPLE) EREAL= EREAL - HSPLIT

      IF(GRADT) THEN

        DO I=1,N

      X1I=X(1,I)
      X2I=X(2,I)
      X3I=X(3,I)
      DXY2=X1I*X1I+X2I*X2I
      RI=SQRT(DXY2+X3I*X3I)
      CI=X3I/RI
        CALL ARNOENERGY(RI,CI,E1RI,E1RI1,E1RI1C,E2RI,E2RI1,E2RI1C)
      DXY=SQRT(DXY2)
      CAI=1.D0
      SAI=0.D0
        IF(DXY.GT.0.D0) THEN
      CAI=X1I/DXY
      SAI=X2I/DXY
                  END IF

      DERI1=E1RI1-E2RI1
      CAI2=CAI*CAI
      E1R1= DERI1*CAI2 + E2RI1
      E2R1=-DERI1*CAI2 + E1RI1
      E12R1= CAI*SAI*DERI1
      DEDRI=0.5D0*( E1R1 + E2R1 + 
     1            ((E1R-E2R)*(E1R1-E2R1) + 4.D0*E12R*E12R1)/D )

      DERI1C=E1RI1C-E2RI1C
      E1R1C= DERI1C*CAI2 + E2RI1C
      E2R1C=-DERI1C*CAI2 + E1RI1C
      E12R1C= CAI*SAI*DERI1C
        DEDCI=0.5D0*( E1R1C + E2R1C +
     1          ((E1R-E2R)*(E1R1C-E2R1C) + 4.D0*E12R*E12R1C)/D )

        DERI=E1RI-E2RI
      E1R1CA= 2.D0*CAI*DERI
CC      E2R1CA= -E1R1CA
      E12R1CA= 0.D0
      IF(SAI.NE.0.D0) E12R1CA= (SAI - CAI2/SAI)*DERI
        DEDCAI= ((E1R-E2R)*E1R1CA + 2.D0*E12R*E12R1CA)/D 

      DO M=1,3
      DRIDXM= X(M,I)/RI
      DCIDXM= -CI/RI*DRIDXM
      IF(M.EQ.3) DCIDXM= DCIDXM +1.D0/RI
      DCAIDXM= 0.D0
      IF(M.LE.2.AND.DXY.GT.0.D0) THEN
      DCAIDXM= -CAI*X(M,I)/DXY2
      IF(M.EQ.1) DCAIDXM= DCAIDXM +1.D0/DXY
               END IF
        DEDXM= DEDRI*DRIDXM + DEDCI*DCIDXM + DEDCAI*DCAIDXM
      IGRAD=3*(I-1)+M
        GRAD(IGRAD)=GRAD(IGRAD)+DEDXM
      ENDDO

        ENDDO

              END IF


      RETURN
      END

      SUBROUTINE ARNOENERGY(R,C1,V,VD,VDC,V1,V1D,V1DC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=9)
      SAVE
      DOUBLE PRECISION P(N),VI(N),VID(N),P1(N)

      IF(ABS(C1).GT.1.D0) C1=C1/ABS(C1)

      V0 = UT0(R,1,V0D)
      N1=N-1
      DO I=1,N1
        VI(I) = UT0(R,I+1,VDI)
      VID(I)=VDI
      ENDDO

        CALL POLEGN(N1,C1,P,P1)

        V = V0 
      VD = V0D
      VDC = 0.D0
      DO I=1,N1
      V = V + VI(I)*P(I)
      VD = VD + VID(I)*P(I)
      VDC = VDC + VI(I)*P1(I)
      ENDDO

        V0 = UT01(R,1,V0D)
      DO I=1,N1
      VI(I) = UT01(R,I+1,VDI)
      VID(I) = VDI
      ENDDO

      V1 = V0
      V1D = V0D
      V1DC = 0.D0
      DO I=1,N1
      V1 = V1 + VI(I)*P(I)
      V1D = V1D + VID(I)*P(I)
      V1DC = V1DC + VI(I)*P1(I)
      ENDDO

      RETURN
      END


      DOUBLE PRECISION FUNCTION UT0(RR,I0,UT0D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER(N0=24,N1=24,N2=24,N3=20,N4=20,N5=20,N6=20,N7=17,N8=18)

      DOUBLE PRECISION, DIMENSION(N0) :: RU=(/ 
     1 0.00, 0.25, 0.50, 0.75, 1.00,
     1 1.25, 1.50, 1.75, 2.00, 2.25,
     1 2.50, 2.75,
     1 3.00, 3.25, 3.50, 3.75, 4.00, 
     1 4.25, 4.50, 5.00, 5.50, 6.00, 
     1 7.00, 8.00/)

      DOUBLE PRECISION, DIMENSION(N0) :: EU=(/
     123.5757810,19.0536413,15.1520331,11.8249968, 9.0265728,
     1 6.7107999, 4.8317169, 3.3433641, 2.1997803, 1.3550050,
     1  .7630777,  .3780378,
     1  .1539246,  .0447773,  .0045563, -.0075913, -.0095564,
     1 -.0083419, -.0064951, -.0035558, -.0019332, -.0010942,
     1 -.0004037, -.0001733/)
      DOUBLE PRECISION, DIMENSION(N1) :: EG=(/
     1 2.9072499, 2.4029129, 1.9608503, 1.5769362, 1.2470459,
     1  .9670573,  .7328453,  .5402856,  .3852534,  .2636248,
     1  .1712753,  .1040805,
     1  .0579161,  .0286581,  .0122599,  .0046522,  .0013999, 
     1  .0001305, -.0002864, -.0003208, -.0001932, -.0001146,
     1 -.0000367, -.0000163/)
      DOUBLE PRECISION, DIMENSION(N2) :: EPU=(/
     136.1544630,29.2137453,23.2280250,18.1262749,13.8374705,
     110.2905828, 7.4145843, 5.1384474, 3.3911439, 2.1016460,
     1 1.1989267,  .6119579,
     1  .2697118,  .1011612,  .0353422,  .0104119,  .0016621,
     1 -.0009844, -.0014811, -.0010339, -.0005530, -.0002932,
     1 -.0000931, -.0000396/)
      DOUBLE PRECISION, DIMENSION(N3) :: E3=(/  
     1 5.7704509, 4.6198981, 3.6335576, 2.7987752, 2.1029035,
     1 1.5332943, 1.0772957,  .7222591,  .4555346,  .2644726,
     1  .1364226,  .0587354,
     1  .0187612,  .0038494,  .0012187,  .0002678, -.0000507,
     1 -.0001010, -.0000828, -.0000240/)
      DOUBLE PRECISION, DIMENSION(N4) :: E4=(/  
     1 6.4428066, 5.1725998, 4.0815667, 3.1560647, 2.3824570,
     1 1.7471023, 1.2363623,  .8365934,  .5341561,  .3154102,
     1  .1667152,  .0744307,
     1  .0249163,  .0045315, -.0003323, -.0010465, -.0007592,
     1 -.0004024, -.0001781, -.0000020/)
      DOUBLE PRECISION, DIMENSION(N5) :: E5=(/  
     1-4.6869956,-3.7386937,-2.9276822,-2.2432417,-1.6746429,
     1-1.2111622, -.8420752, -.5566573, -.3441831, -.1939286,
     1 -.0951689, -.0371794,
     1 -.0092352, -.0006112, -.0004883, -.0002925, -.0001497,
     1 -.0000855, -.0000426, -.0000077/)
      DOUBLE PRECISION, DIMENSION(N6) :: E6=(/ 
     1-4.1197525,-3.2840736,-2.5697842,-1.9673653,-1.4672967,
     1-1.0600520, -.7361093, -.4859449, -.3000351, -.1688565,
     1 -.0828862, -.0326009,
     1 -.0084771, -.0009921, -.0007325, -.0004486, -.0002404,
     1 -.0001257, -.0000563, -.0000099/)
      DOUBLE PRECISION, DIMENSION(N7) :: E7=(/  
     1 2.0150222, 1.6081385, 1.2600261,  .9661046,  .7217950,
     1  .5225175,  .3636947,  .2407501,  .1491035,  .0841764,
     1  .0413916,  .0161697,
     1  .0039326,  .0001014,  .0000209, -.0000055, -.0000049/)
      DOUBLE PRECISION, DIMENSION(N8) :: E8=(/
     1 2.8090209, 2.2393604, 1.7523957, 1.3416469, 1.0006275,
     1  .7228596,  .5018594,  .3311493,  .2042477,  .1146736,
     1  .0559460,  .0215842,
     1  .0051075,  .0000354, -.0000178, -.0000116, -.0000069,
     1 -.0000047/)

      DOUBLE PRECISION COEFU(4,N0), COEFG(4,N1), COEFPU(4,N2), COEF3(4,N3),
     1COEF4(4,N4), COEF5(4,N5), COEF6(4,N6), COEF7(4,N7), COEF8(4,N8)
      DATA 
     1 NC0/1/,NC1/1/,NC2/1/,NC3/1/,NC4/1/,NC5/1/,NC6/1/,NC7/1/,NC8/1/

      UT0=0.D0

      IF(I0.EQ.1) THEN
      IF(NC0.EQ.1) THEN
      NC0=2
      DO I=1,N0
      COEFU(1,I)=EU(I) 
      ENDDO
      COEFU(2,N0)=-6.D0*EU(N0)/RU(N0)
      CALL CUBSPL(RU,COEFU,N0,0,1)
        R0U=RU(1)
        AAU=EU(1)
        ALPHAU=-COEFU(2,1)/AAU
                   END IF

CC      RR=DMAX1(RR0,R0U)
        IF(RR.LT.RU(N0)) THEN
        IF(RR.LT.R0U) THEN
        UT0=AAU*EXP(-ALPHAU*(RR-R0U))
      UT0D=-ALPHAU*UT0
                        ELSE
      DO I=1,N0-1
      IND=I
      IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 10
      ENDDO
   10      CONTINUE
      DR=RR-RU(IND)
      UT0=COEFU(1,IND)+DR*(COEFU(2,IND)
     1      +0.5D0*DR*(COEFU(3,IND)+DR*COEFU(4,IND)/3.D0))
      UT0D=COEFU(2,IND) + DR*(COEFU(3,IND) + 0.5D0*DR*COEFU(4,IND))
                        END IF
                         ELSE
            UT0=EU(N0)*(RU(N0)/RR)**6
        UT0D=-6.D0*UT0/RR
                   END IF
      RETURN
                  END IF

      IF(I0.EQ.2) THEN
      IF(NC1.EQ.1) THEN
      NC1=2
      DO I=1,N1
      COEFG(1,I)=EG(I)
      ENDDO
      COEFG(2,N1)=-6.D0*EG(N1)/RU(N1)
      CALL CUBSPL(RU,COEFG,N1,0,1)
        R0G=RU(1)
        AAG=EG(1)
        ALPHAG=-COEFG(2,1)/AAG
                   END IF

CC        RR=DMAX1(RR0,R0G)
        IF(RR.LT.RU(N1)) THEN
        IF(RR.LT.R0G) THEN
        UT0=AAG*EXP(-ALPHAG*(RR-R0G))
      UT0D=-ALPHAG*UT0
                        ELSE
        DO I=1,N1-1
        IND=I
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 11
        ENDDO
   11   CONTINUE
        DR=RR-RU(IND)
        UT0=COEFG(1,IND)+DR*(COEFG(2,IND)
     1      +0.5D0*DR*(COEFG(3,IND)+DR*COEFG(4,IND)/3.D0))        
      UT0D=COEFG(2,IND) + DR*(COEFG(3,IND) + 0.5D0*DR*COEFG(4,IND))
                        END IF
                         ELSE
        UT0=EG(N1)*(RU(N1)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
      RETURN
                  END IF

      IF(I0.EQ.3) THEN
      IF(NC2.EQ.1) THEN
      NC2=2
      DO I=1,N2
      COEFPU(1,I)=EPU(I)
      ENDDO
      COEFPU(2,N2)=-6.D0*EPU(N2)/RU(N2)
      CALL CUBSPL(RU,COEFPU,N2,0,1)
        R0PU=RU(1)
        AAPU=EPU(1)
        ALPHAPU=-COEFPU(2,1)/AAPU
                   END IF
CC        RR=DMAX1(RR0,R0PU)
      IF(RR.LT.RU(N2)) THEN
       IF(RR.LT.R0PU) THEN
        UT0=AAPU*EXP(-ALPHAPU*(RR-R0PU))
      UT0D=-ALPHAPU*UT0
                        ELSE
        DO I=1,N2-1                     
        IND=I                     
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 12
        ENDDO                     
   12    CONTINUE                     
        DR=RR-RU(IND)
        UT0=COEFPU(1,IND)+DR*(COEFPU(2,IND)
     1      +0.5D0*DR*(COEFPU(3,IND)+DR*COEFPU(4,IND)/3.D0))
      UT0D=COEFPU(2,IND) + DR*(COEFPU(3,IND) + 0.5D0*DR*COEFPU(4,IND))
                        END IF
                         ELSE
        UT0=EPU(N2)*(RU(N2)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
      RETURN
                  END IF

      IF(I0.EQ.4) THEN
      IF(NC3.EQ.1) THEN
      NC3=2
      DO I=1,N3
      COEF3(1,I)=E3(I)
      ENDDO
      COEF3(2,N3)=-6.D0*E3(N3)/RU(N3)
      CALL CUBSPL(RU,COEF3,N3,0,1)
        R03=RU(1)
        AA3=E3(1)
        ALPHA3=-COEF3(2,1)/AA3
                   END IF

CC        RR=DMAX1(RR0,R03)
        IF(RR.LT.RU(N3)) THEN
        IF(RR.LT.R03) THEN
        UT0=AA3*EXP(-ALPHA3*(RR-R03))
      UT0D=-ALPHA3*UT0
                        ELSE
        DO I=1,N3-1                     
        IND=I                     
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 13
        ENDDO                     
   13   CONTINUE                     
        DR=RR-RU(IND)
        UT0=COEF3(1,IND)+DR*(COEF3(2,IND)
     1  +0.5D0*DR*(COEF3(3,IND)+DR*COEF3(4,IND)/3.D0))
      UT0D=COEF3(2,IND) + DR*(COEF3(3,IND) + 0.5D0*DR*COEF3(4,IND))
                        END IF
                         ELSE
        UT0=E3(N3)*(RU(N3)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
      RETURN
                  END IF

      IF(I0.EQ.5) THEN
      IF(NC4.EQ.1) THEN
      NC4=2
      DO I=1,N4
      COEF4(1,I)=E4(I)
      ENDDO
      COEF4(2,N4)=-6.D0*E4(N4)/RU(N4)
      CALL CUBSPL(RU,COEF4,N4,0,1)
        R04=RU(1)
        AA4=E4(1)
        ALPHA4=-COEF4(2,1)/AA4
                   END IF

CC        RR=DMAX1(RR0,R04)
        IF(RR.LT.RU(N4)) THEN
        IF(RR.LT.R04) THEN
        UT0=AA4*EXP(-ALPHA4*(RR-R04))
      UT0D=-ALPHA4*UT0
                        ELSE
        DO I=1,N4-1
        IND=I
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 14
        ENDDO
   14   CONTINUE 
        DR=RR-RU(IND)
        UT0=COEF4(1,IND)+DR*(COEF4(2,IND)
     1  +0.5D0*DR*(COEF4(3,IND)+DR*COEF4(4,IND)/3.D0))
      UT0D=COEF4(2,IND) + DR*(COEF4(3,IND) + 0.5D0*DR*COEF4(4,IND))
                        END IF
                         ELSE
        UT0=E4(N4)*(RU(N4)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
      RETURN
                  END IF

      IF(I0.EQ.6) THEN
      IF(NC5.EQ.1) THEN
      NC5=2
      DO I=1,N5
      COEF5(1,I)=E5(I)
      ENDDO
      COEF5(2,N5)=-6.D0*E5(N5)/RU(N5)
      CALL CUBSPL(RU,COEF5,N5,0,1)
        R05=RU(1)
        AA5=E5(1)
        ALPHA5=-COEF5(2,1)/AA5
                   END IF

CC        RR=DMAX1(RR0,R05)
        IF(RR.LT.RU(N5)) THEN
        IF(RR.LT.R05) THEN
        UT0=AA5*EXP(-ALPHA5*(RR-R05))
      UT0D=-ALPHA5*UT0
                        ELSE
        DO I=1,N5-1                     
        IND=I                     
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 15                      
        ENDDO                     
   15   CONTINUE                      
        DR=RR-RU(IND)
        UT0=COEF5(1,IND)+DR*(COEF5(2,IND)
     1  +0.5D0*DR*(COEF5(3,IND)+DR*COEF5(4,IND)/3.D0))
      UT0D=COEF5(2,IND) + DR*(COEF5(3,IND) + 0.5D0*DR*COEF5(4,IND))
                        END IF
                         ELSE
        UT0=E5(N5)*(RU(N5)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
      RETURN
                  END IF

      IF(I0.EQ.7) THEN
      IF(NC6.EQ.1) THEN
      NC6=2
      DO I=1,N6
      COEF6(1,I)=E6(I)
      ENDDO
      COEF6(2,N6)=-6.D0*E6(N6)/RU(N6)
      CALL CUBSPL(RU,COEF6,N6,0,1)
        R06=RU(1)
        AA6=E6(1)
        ALPHA6=-COEF6(2,1)/AA6
                   END IF

CC        RR=DMAX1(RR0,R06)
      IF(RR.LT.RU(N6)) THEN
        IF(RR.LT.R06) THEN
        UT0=AA6*EXP(-ALPHA6*(RR-R06))
      UT0D=-ALPHA6*UT0
                        ELSE
        DO I=1,N6-1                                          
        IND=I                                          
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 16
        ENDDO                                          
   16   CONTINUE                                           
        DR=RR-RU(IND)
        UT0=COEF6(1,IND)+DR*(COEF6(2,IND)
     1  +0.5D0*DR*(COEF6(3,IND)+DR*COEF6(4,IND)/3.D0))
      UT0D=COEF6(2,IND) + DR*(COEF6(3,IND) + 0.5D0*DR*COEF6(4,IND))
                        END IF
                         ELSE
        UT0=E6(N6)*(RU(N6)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
      RETURN
                  END IF

      IF(I0.EQ.8) THEN
      IF(NC7.EQ.1) THEN
      NC7=2
      DO I=1,N7
      COEF7(1,I)=E7(I)
      ENDDO
      COEF7(2,N7)=-6.D0*E7(N7)/RU(N7)
      CALL CUBSPL(RU,COEF7,N7,0,1)
        R07=RU(1)
        AA7=E7(1)
        ALPHA7=-COEF7(2,1)/AA7
                   END IF

CC      RR=DMAX1(RR0,R07)
      IF(RR.LT.RU(N7)) THEN
        IF(RR.LT.R07) THEN
        UT0=AA7*EXP(-ALPHA7*(RR-R07))
      UT0D=-ALPHA7*UT0
                        ELSE
        DO I=1,N7-1                                                             
        IND=I                                                               
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 17                     
        ENDDO                                                               
   17   CONTINUE                                                                
        DR=RR-RU(IND)
        UT0=COEF7(1,IND)+DR*(COEF7(2,IND)
     1  +0.5D0*DR*(COEF7(3,IND)+DR*COEF7(4,IND)/3.D0))
      UT0D=COEF7(2,IND) + DR*(COEF7(3,IND) + 0.5D0*DR*COEF7(4,IND))
                        END IF
                         ELSE
        UT0=E7(N7)*(RU(N7)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
      RETURN
                  END IF

      IF(I0.EQ.9) THEN
      IF(NC8.EQ.1) THEN
      NC8=2
      DO I=1,N8
      COEF8(1,I)=E8(I)
      ENDDO
      COEF8(2,N8)=-6.D0*E8(N8)/RU(N8)
      CALL CUBSPL(RU,COEF8,N8,0,1)
        R08=RU(1)
        AA8=E8(1)
        ALPHA8=-COEF8(2,1)/AA8
                   END IF

CC      RR=DMAX1(RR0,R08)
      IF(RR.LT.RU(N8)) THEN
        IF(RR.LT.R08) THEN
        UT0=AA8*EXP(-ALPHA8*(RR-R08))
      UT0D=-ALPHA8*UT0
                        ELSE
        DO I=1,N8-1  
        IND=I        
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 18                              
        ENDDO        
   18   CONTINUE     
        DR=RR-RU(IND)
        UT0=COEF8(1,IND)+DR*(COEF8(2,IND)
     1  +0.5D0*DR*(COEF8(3,IND)+DR*COEF8(4,IND)/3.D0))
      UT0D=COEF8(2,IND) + DR*(COEF8(3,IND) + 0.5D0*DR*COEF8(4,IND))
                        END IF
                         ELSE
        UT0=E8(N8)*(RU(N8)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
      RETURN
                  END IF

      RETURN 
      END

      SUBROUTINE POLEGN(N,X,P,P1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION P(N),P1(N)
      P(1)=X
      P1(1)=1.D0
      IF(N.GT.1) THEN
      PN0=1.D0
      PN=X
      PN0D=0.D0
      PND=1.D0
      N1=N-1
      DO I=1,N1
      RI=1.0D0*I
      PN1=((2.D0*RI+1.D0)*X*PN - RI*PN0)/(RI+1.D0)
      PN1D=((2.D0*RI+1.D0)*(PN+X*PND) - RI*PN0D)/(RI+1.D0)
      PN0=PN
      PN=PN1
      P(I+1)=PN1
      PN0D=PND
      PND=PN1D
      P1(I+1)=PN1D
      ENDDO
               END IF
      RETURN
      END

C------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION UT01(RR,I0,UT01D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      PARAMETER(N0=13,N1=12,N2=12,N3=10,N4=10,N5=8,N6=8,N7=6,N8=6)

      DOUBLE PRECISION, DIMENSION(N0) :: RU=(/ 
     1 3.00, 3.25, 3.50, 3.75, 4.00, 
     1 4.25, 4.50, 5.00, 5.50, 6.00, 
     1 7.00, 8.00, 10.0/)

      DOUBLE PRECISION, DIMENSION(N0) :: EU=(/
     1  .1220345,  .0289952, -.0023568, -.0104175, -.0105767,
     1 -.0086238, -.0065077, -.0034779, -.0018770, -.0010615,
     1 -.0003881, -.0001650, -.0000238/)
      DOUBLE PRECISION, DIMENSION(N1) :: EG=(/
     1  .0275582,  .0126120,  .0047627,  .0012348, -.0001393,
     1 -.0005516, -.0005758, -.0003482, -.0001767, -.0000959,
     1 -.0000267, -.0000039/)
      DOUBLE PRECISION, DIMENSION(N2) :: EPU=(/
     1  .2295942,  .0839463,  .0271378,  .0064960, -.0002063,
     1 -.0018637, -.0018941, -.0011290, -.0005799, -.0003031,
     1 -.0000955, -.0000272/)
      DOUBLE PRECISION, DIMENSION(N3) :: E3=(/  
     1  .0347188,  .0157193,  .0067817,  .0027936,  .0010881,
     1  .0004036,  .0001308, -.0000088, -.0000128, -.0000039/)
      DOUBLE PRECISION, DIMENSION(N4) :: E4=(/  
     1  .0799795,  .0316737,  .0122230,  .0045691,  .0016355,
     1  .0005520,  .0001638,  .0000002, -.0000114, -.0000065/)
      DOUBLE PRECISION, DIMENSION(N5) :: E5=(/  
     1  .0078080,  .0032611,  .0013291,  .0005546,  .0002360,
     1  .0000854,  .0000302,  .0000003/)
      DOUBLE PRECISION, DIMENSION(N6) :: E6=(/  
     1  .0121196,  .0044330,  .0016626,  .0006155,  .0002380,
     1  .0000723,  .0000263,  .0000035/)
      DOUBLE PRECISION, DIMENSION(N7) :: E7=(/  
     1  .0012897,  .0004052,  .0001378,  .0000390,  .0000099,
     1  .0000014/)
      DOUBLE PRECISION, DIMENSION(N8) :: E8=(/
     1  .0014554,  .0004653,  .0001494,  .0000507,  .0000087,
     1  .0000041/)

      DOUBLE PRECISION COEFU(4,N0), COEFG(4,N1), COEFPU(4,N2), COEF3(4,N3),
     1COEF4(4,N4), COEF5(4,N5), COEF6(4,N6), COEF7(4,N7), COEF8(4,N8)
      DATA 
     1 NC0/1/,NC1/1/,NC2/1/,NC3/1/,NC4/1/,NC5/1/,NC6/1/,NC7/1/,NC8/1/

      UT0=0.D0
      UT0D=0.D0

      IF(I0.EQ.1) THEN
      IF(NC0.EQ.1) THEN
      NC0=2
      DO I=1,N0
      COEFU(1,I)=EU(I) 
      ENDDO
      COEFU(2,N0)=-6.D0*EU(N0)/RU(N0)
      CALL CUBSPL(RU,COEFU,N0,0,1)
        R0U=RU(1)
        AAU=EU(1)
        ALPHAU=-COEFU(2,1)/AAU
                   END IF

CC      RR=DMAX1(RR0,R0U)
        IF(RR.LT.RU(N0)) THEN
        IF(RR.LT.R0U) THEN
        UT0=AAU*EXP(-ALPHAU*(RR-R0U))
      UT0D=-ALPHAU*UT0
                        ELSE
      DO I=1,N0-1
      IND=I
      IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 10
      ENDDO
   10      CONTINUE
      DR=RR-RU(IND)
      UT0=COEFU(1,IND)+DR*(COEFU(2,IND)
     1      +0.5D0*DR*(COEFU(3,IND)+DR*COEFU(4,IND)/3.D0))
      UT0D=COEFU(2,IND) + DR*(COEFU(3,IND) + 0.5D0*DR*COEFU(4,IND))
                        END IF
                         ELSE
            UT0=EU(N0)*(RU(N0)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
      UT01=UT0
      UT01D=UT0D
      RETURN
                  END IF

      IF(I0.EQ.2) THEN
      IF(NC1.EQ.1) THEN
      NC1=2
      DO I=1,N1
      COEFG(1,I)=EG(I)
      ENDDO
      COEFG(2,N1)=-6.D0*EG(N1)/RU(N1)
      CALL CUBSPL(RU,COEFG,N1,0,1)
        R0G=RU(1)
        AAG=EG(1)
        ALPHAG=-COEFG(2,1)/AAG
                   END IF

CC        RR=DMAX1(RR0,R0G)
        IF(RR.LT.RU(N1)) THEN
        IF(RR.LT.R0G) THEN
        UT0=AAG*EXP(-ALPHAG*(RR-R0G))
      UT0D=-ALPHAG*UT0
                        ELSE
        DO I=1,N1-1
        IND=I
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 11
        ENDDO
   11   CONTINUE
        DR=RR-RU(IND)
        UT0=COEFG(1,IND)+DR*(COEFG(2,IND)
     1      +0.5D0*DR*(COEFG(3,IND)+DR*COEFG(4,IND)/3.D0))        
      UT0D=COEFG(2,IND) + DR*(COEFG(3,IND) + 0.5D0*DR*COEFG(4,IND))
                        END IF
                         ELSE
        UT0=EG(N1)*(RU(N1)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
        UT01=UT0
      UT01D=UT0D
      RETURN
                  END IF

      IF(I0.EQ.3) THEN
      IF(NC2.EQ.1) THEN
      NC2=2
      DO I=1,N2
      COEFPU(1,I)=EPU(I)
      ENDDO
      COEFPU(2,N2)=-6.D0*EPU(N2)/RU(N2)
      CALL CUBSPL(RU,COEFPU,N2,0,1)
        R0PU=RU(1)
        AAPU=EPU(1)
        ALPHAPU=-COEFPU(2,1)/AAPU
                   END IF
CC        RR=DMAX1(RR0,R0PU)
      IF(RR.LT.RU(N2)) THEN
       IF(RR.LT.R0PU) THEN
        UT0=AAPU*EXP(-ALPHAPU*(RR-R0PU))
      UT0D=-ALPHAPU*UT0
                        ELSE
        DO I=1,N2-1                     
        IND=I                     
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 12
        ENDDO                     
   12    CONTINUE                     
        DR=RR-RU(IND)
        UT0=COEFPU(1,IND)+DR*(COEFPU(2,IND)
     1      +0.5D0*DR*(COEFPU(3,IND)+DR*COEFPU(4,IND)/3.D0))
      UT0D=COEFPU(2,IND) + DR*(COEFPU(3,IND) + 0.5D0*DR*COEFPU(4,IND))
                        END IF
                         ELSE
        UT0=EPU(N2)*(RU(N2)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
        UT01=UT0
      UT01D=UT0D
      RETURN
                  END IF

      IF(I0.EQ.4) THEN
      IF(NC3.EQ.1) THEN
      NC3=2
      DO I=1,N3
      COEF3(1,I)=E3(I)
      ENDDO
      COEF3(2,N3)=-6.D0*E3(N3)/RU(N3)
      CALL CUBSPL(RU,COEF3,N3,0,1)
        R03=RU(1)
        AA3=E3(1)
        ALPHA3=-COEF3(2,1)/AA3
                   END IF

CC        RR=DMAX1(RR0,R03)
        IF(RR.LT.RU(N3)) THEN
        IF(RR.LT.R03) THEN
        UT0=AA3*EXP(-ALPHA3*(RR-R03))
      UT0D=-ALPHA3*UT0
                        ELSE
        DO I=1,N3-1                     
        IND=I                     
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 13
        ENDDO                     
   13   CONTINUE                     
        DR=RR-RU(IND)
        UT0=COEF3(1,IND)+DR*(COEF3(2,IND)
     1  +0.5D0*DR*(COEF3(3,IND)+DR*COEF3(4,IND)/3.D0))
      UT0D=COEF3(2,IND) + DR*(COEF3(3,IND) + 0.5D0*DR*COEF3(4,IND))
                        END IF
                         ELSE
        UT0=E3(N3)*(RU(N3)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
        UT01=UT0
      UT01D=UT0D
      RETURN
                  END IF

      IF(I0.EQ.5) THEN
      IF(NC4.EQ.1) THEN
      NC4=2
      DO I=1,N4
      COEF4(1,I)=E4(I)
      ENDDO
      COEF4(2,N4)=-6.D0*E4(N4)/RU(N4)
      CALL CUBSPL(RU,COEF4,N4,0,1)
        R04=RU(1)
        AA4=E4(1)
        ALPHA4=-COEF4(2,1)/AA4
                   END IF

CC        RR=DMAX1(RR0,R04)
        IF(RR.LT.RU(N4)) THEN
        IF(RR.LT.R04) THEN
        UT0=AA4*EXP(-ALPHA4*(RR-R04))
      UT0D=-ALPHA4*UT0
                        ELSE
        DO I=1,N4-1
        IND=I
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 14
        ENDDO
   14   CONTINUE 
        DR=RR-RU(IND)
        UT0=COEF4(1,IND)+DR*(COEF4(2,IND)
     1  +0.5D0*DR*(COEF4(3,IND)+DR*COEF4(4,IND)/3.D0))
      UT0D=COEF4(2,IND) + DR*(COEF4(3,IND) + 0.5D0*DR*COEF4(4,IND))
                        END IF
                         ELSE
        UT0=E4(N4)*(RU(N4)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
        UT01=UT0
      UT01D=UT0D
      RETURN
                  END IF

      IF(I0.EQ.6) THEN
      IF(NC5.EQ.1) THEN
      NC5=2
      DO I=1,N5
      COEF5(1,I)=E5(I)
      ENDDO
      COEF5(2,N5)=-6.D0*E5(N5)/RU(N5)
      CALL CUBSPL(RU,COEF5,N5,0,1)
        R05=RU(1)
        AA5=E5(1)
        ALPHA5=-COEF5(2,1)/AA5
                   END IF

CC        RR=DMAX1(RR0,R05)
        IF(RR.LT.RU(N5)) THEN
        IF(RR.LT.R05) THEN
        UT0=AA5*EXP(-ALPHA5*(RR-R05))
      UT0D=-ALPHA5*UT0
                        ELSE
        DO I=1,N5-1                     
        IND=I                     
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 15                      
        ENDDO                     
   15   CONTINUE                      
        DR=RR-RU(IND)
        UT0=COEF5(1,IND)+DR*(COEF5(2,IND)
     1  +0.5D0*DR*(COEF5(3,IND)+DR*COEF5(4,IND)/3.D0))
      UT0D=COEF5(2,IND) + DR*(COEF5(3,IND) + 0.5D0*DR*COEF5(4,IND))
                        END IF
                         ELSE
        UT0=E5(N5)*(RU(N5)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
        UT01=UT0
      UT01D=UT0D
      RETURN
                  END IF

      IF(I0.EQ.7) THEN
      IF(NC6.EQ.1) THEN
      NC6=2
      DO I=1,N6
      COEF6(1,I)=E6(I)
      ENDDO
      COEF6(2,N6)=-6.D0*E6(N6)/RU(N6)
      CALL CUBSPL(RU,COEF6,N6,0,1)
        R06=RU(1)
        AA6=E6(1)
        ALPHA6=-COEF6(2,1)/AA6
                   END IF

CC        RR=DMAX1(RR0,R06)
      IF(RR.LT.RU(N6)) THEN
        IF(RR.LT.R06) THEN
        UT0=AA6*EXP(-ALPHA6*(RR-R06))
      UT0D=-ALPHA6*UT0
                        ELSE
        DO I=1,N6-1                                          
        IND=I                                          
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 16
        ENDDO                                          
   16   CONTINUE                                           
        DR=RR-RU(IND)
        UT0=COEF6(1,IND)+DR*(COEF6(2,IND)
     1  +0.5D0*DR*(COEF6(3,IND)+DR*COEF6(4,IND)/3.D0))
      UT0D=COEF6(2,IND) + DR*(COEF6(3,IND) + 0.5D0*DR*COEF6(4,IND))
                        END IF
                         ELSE
        UT0=E6(N6)*(RU(N6)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
        UT01=UT0
      UT01D=UT0D
      RETURN
                  END IF

      IF(I0.EQ.8) THEN
      IF(NC7.EQ.1) THEN
      NC7=2
      DO I=1,N7
      COEF7(1,I)=E7(I)
      ENDDO
      COEF7(2,N7)=-6.D0*E7(N7)/RU(N7)
      CALL CUBSPL(RU,COEF7,N7,0,1)
        R07=RU(1)
        AA7=E7(1)
        ALPHA7=-COEF7(2,1)/AA7
                   END IF

CC      RR=DMAX1(RR0,R07)
      IF(RR.LT.RU(N7)) THEN
        IF(RR.LT.R07) THEN
        UT0=AA7*EXP(-ALPHA7*(RR-R07))
      UT0D=-ALPHA7*UT0
                        ELSE
        DO I=1,N7-1                                                             
        IND=I                                                               
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 17                     
        ENDDO                                                               
   17   CONTINUE                                                                
        DR=RR-RU(IND)
        UT0=COEF7(1,IND)+DR*(COEF7(2,IND)
     1  +0.5D0*DR*(COEF7(3,IND)+DR*COEF7(4,IND)/3.D0))
      UT0D=COEF7(2,IND) + DR*(COEF7(3,IND) + 0.5D0*DR*COEF7(4,IND))
                        END IF
                         ELSE
        UT0=E7(N7)*(RU(N7)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
        UT01=UT0
      UT01D=UT0D
      RETURN
                  END IF

      IF(I0.EQ.9) THEN
      IF(NC8.EQ.1) THEN
      NC8=2
      DO I=1,N8
      COEF8(1,I)=E8(I)
      ENDDO
      COEF8(2,N8)=-6.D0*E8(N8)/RU(N8)
      CALL CUBSPL(RU,COEF8,N8,0,1)
        R08=RU(1)
        AA8=E8(1)
        ALPHA8=-COEF8(2,1)/AA8
                   END IF

CC      RR=DMAX1(RR0,R08)
      IF(RR.LT.RU(N8)) THEN
        IF(RR.LT.R08) THEN
        UT0=AA8*EXP(-ALPHA8*(RR-R08))
      UT0D=-ALPHA8*UT0
                        ELSE
        DO I=1,N8-1  
        IND=I        
        IF(RU(I).LE.RR.AND.RR.LT.RU(I+1)) GO TO 18                              
        ENDDO        
   18   CONTINUE     
        DR=RR-RU(IND)
        UT0=COEF8(1,IND)+DR*(COEF8(2,IND)
     1  +0.5D0*DR*(COEF8(3,IND)+DR*COEF8(4,IND)/3.D0))
      UT0D=COEF8(2,IND) + DR*(COEF8(3,IND) + 0.5D0*DR*COEF8(4,IND))
                        END IF
                         ELSE
        UT0=E8(N8)*(RU(N8)/RR)**6
      UT0D=-6.D0*UT0/RR
                         END IF
        UT01=UT0
      UT01D=UT0D
      RETURN
                  END IF

        UT01=UT0
      UT01D=UT0D
      RETURN 
      END


