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
C       TB POTENTIAL FOR NA, AG, LI CLUSTERS

        SUBROUTINE NATB(NSIZE, P, GRAD, E0, GTEST, GUIDET)
        IMPLICIT NONE
        INTEGER ION, NSIZE
        DOUBLE PRECISION, PARAMETER :: EPSILON=1D-10
        LOGICAL GTEST, STEST, GUIDET
        DOUBLE PRECISION P(3*NSIZE), GRAD(3*NSIZE), E0
        CHARACTER(LEN=2) METAL
        COMMON/ION/ION
        COMMON/METAL/METAL

        METAL='NA'
        ION=0

C       IF (METAL.EQ.'NA') THEN
C          PRINT *,'PARAMETERS OPTIMIZED FOR SODIUM'
C       ELSEIF (METAL.EQ.'AG') THEN
C          PRINT *,'PARAMETERS OPTIMIZED FOR SILVER'
C       ELSEIF (METAL.EQ.'LI') THEN
C          PRINT *,'PARAMETERS OPTIMIZED FOR LITHIUM'
C       ENDIF

C       IF (ION.EQ.1) THEN
C          PRINT *,'CATION'
C       ELSEIF (ION.EQ.0) THEN
C          PRINT *,'NEUTRAL'
C       ELSE
C          PRINT *,'ANION'
C       ENDIF

C       PRINT*,'GUIDET=',GUIDET
        IF (GUIDET) P(1:3*NSIZE)=P(1:3*NSIZE)*6.02D0
        CALL ENTOTS(P,GRAD,NSIZE,E0)
        IF (GUIDET) P(1:3*NSIZE)=P(1:3*NSIZE)/6.02D0
        IF (GUIDET) GRAD(1:3*NSIZE)=GRAD(1:3*NSIZE)*6.02D0

        END
C______________________________________________________________________

        SUBROUTINE ENTOTS(P,DP,NBAT,E0)
        IMPLICIT REAL*8(A-H,O-Z)
        PARAMETER(ANM=41901.052,RKB=3.165D-6,EREPMAX=100.D0)

        DIMENSION P(3*NBAT)
        DIMENSION GRAD(3*NBAT),DP(3*NBAT)

        DIMENSION HDER(3*NBAT,NBAT*(NBAT+1)/2),R(NBAT*(NBAT+1)/2)
        DIMENSION HA(NBAT*(NBAT+1)/2)
        DIMENSION VECT(NBAT*(NBAT+1)/2),RINV(NBAT*(NBAT+1)/2)
        DIMENSION H(NBAT*(NBAT+1)/2),DIF(3*NBAT*(NBAT+1)/2)
        DIMENSION DBSZ(3*NBAT*(NBAT+1)/2),UNIT(3*NBAT*(NBAT+1)/2)
        DIMENSION VALP(NBAT),VT(NBAT,NBAT),NUM(NBAT)
        DOUBLE PRECISION HMAT(NBAT,NBAT)
        INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NBAT)
        INTEGER IWORK(33*3*NBAT), INFO, ISTAT
        DOUBLE PRECISION WORK(33*3*NBAT), ABSTOL, DLAMCH
        INTEGER ION
        COMMON/ION/ION

        DIMENSION DI(3),DJ(3)

        LOGICAL MONOC
        CHARACTER*2 METAL
        COMMON/METAL/METAL

C        PARAMETERS

         LWORK=33*3*NBAT
         ILWORK=33*3*NBAT

        IF (METAL.EQ.'NA') THEN
           ESP=0.077307
           RMIN=3.
           RMAX=100.
        ELSEIF (METAL.EQ.'AG') THEN
           ESP=0.227857
           RMIN=4.
           RMAX=20.
        ELSEIF (METAL.EQ.'LI') THEN
           ESP=0.067908D0
           RMIN=0.
           RMAX=1000.
        ENDIF

        NALK=NBAT
        IDIAG=1
        DISTOR=0.D0
        HH=1D-6
        FTOL=1D-4
        ITMAX=100
        IPARAM=1
        ASZ=1
        NGR=0
        IPRTH=0
        IPRTG=0
        NBLEC=NBAT-ION
        APOL=0.

        NDIMH=NBAT ! CHANGED BY DJW

        DO IJ=1,NBAT
           H(IJ)=0.D0
           DO II=1,3*NBAT
              HDER(II,IJ)=0.D0
           ENDDO
        ENDDO
        DO I=1,3*NBAT
           GRAD(I)=0.D0
        ENDDO
        DENOM=1.D0/ESP 
        
        IJ=0
        IJQ=0
        DO 2680 I=1,NBAT
           DO 2682 J=1,I-1
              IJ=IJ+1
              IJX=IJQ+1
              IJY=IJQ+2
              IJZ=IJQ+3 
              DIF(IJX)=P(3*J-2)-P(3*I-2)
              DIF(IJY)=P(3*J-1)-P(3*I-1)
              DIF(IJZ)=P(3*J)-P(3*I)
              RIJ=DIF(IJX)*DIF(IJX)+DIF(IJY)*DIF(IJY)+DIF(IJZ)*DIF(IJZ)
              RIJ=DSQRT(RIJ)
              R(IJ)=RIJ
              RINV(IJ)=1.D0/RIJ
              UNIT(IJX)=DIF(IJX)*RINV(IJ)
              UNIT(IJY)=DIF(IJY)*RINV(IJ)
              UNIT(IJZ)=DIF(IJZ)*RINV(IJ)
 2682           IJQ=IJQ+3
           IJ=IJ+1
           IJQ=IJQ+3
 2680        CONTINUE
        
C-----------------------------
C       INITIALISATION 
C------------------------------

        DO 5500 IJ=1,NBAT*(NBAT+1)/2
           DO 5500 K=1,3*NBAT
 5500   HDER(K,IJ)=0.D0

C       
C-------------------------------------------------------
C       PRECALCUL ET STOCKAGE DES TERMES DE PERTURBATION
C-------------------------------------------------------
C       
        IJQ=0
        IJ=0
        DO 5700 I=1,NALK
           DO 5720 J=1,I-1
              IJ=IJ+1
              RIJ=R(IJ)
              IF ((RIJ.LT.RMAX).AND.(RIJ.GE.RMIN)) THEN
                 CALL FTSZ(RIJ,VAL,DER)
              ELSE IF (RIJ.LT.RMIN) THEN
                 IF (METAL.EQ.'NA') THEN
                    VAL=-3.724D-10*DEXP(3.32*RIJ)
                    DER=3.32*VAL
                 ELSEIF (METAL.EQ.'AG') THEN
                    VAL=-3.724D-10*DEXP(3.32*RIJ)
                    DER=3.32*VAL
                 ELSEIF (METAL.EQ.'LI') THEN
                    VAL=-3.724D-10*DEXP(3.32*RIJ)
                    DER=3.32*VAL
                 ENDIF
              ELSE
                 VAL=0.D0
                 DER=0.D0
              ENDIF
 1346              VECT(IJ)=VAL
C       
C       DERIVEE DE BSZ
C       
              RINV(IJ)=1.D0/RIJ 
              DO L=1,3
                 IJQL=IJQ+L
                 DBSZ(IJQL)=DER*UNIT(IJQL)
              ENDDO
              IJQ=IJQ+3
 5720           CONTINUE
           IJQ=IJQ+3
           IJ=IJ+1
 5700        CONTINUE 
C       
        
C----------------------------------------------
C       HAMILTONIEN ET GRADIENT DE L HAMILTONIEN
C-----------------------------------------------
        IJ=0
        III=0
        IJQ=0
        DO 30 I=1,NALK
           JJJ=0
C-------------------------------------------
C       ELEMENTS EXTRA-DIAGONAUX (BETA EFFECTIF)
C-------------------------------------------
           DO 32 J=1,I-1
              IJ=IJ+1
C       
C       HAMILTONIEN SS
C       
              RIJ=R(IJ)
              IF(RIJ.GE.RMAX) THEN
                 HSS=0.D0
                 DER=0.D0
              ELSE
                 IF (RIJ.GE.RMIN) THEN
                    CALL FTSS(RIJ,VAL,DER)
                 ELSE
                    IF (METAL.EQ.'NA') THEN
                       VAL=-3.743D-8*DEXP(2.544*RIJ)
                       DER=2.544*VAL
                    ELSEIF (METAL.EQ.'AG') THEN
                       VAL=-3.743D-8*DEXP(2.544*RIJ)
                       DER=2.544*VAL
                    ELSEIF (METAL.EQ.'LI') THEN
                       VAL=-3.743D-8*DEXP(2.544*RIJ)
                       DER=2.544*VAL
                    ENDIF
                 ENDIF
                 HSS=VAL
              ENDIF
C
C       GRADIENT DES TERMES SS
C
              IDER=0
              DO 2504 L=1,3 
                 IJQL=IJQ+L
                 DERQ=DER*UNIT(IJQL)
                 HDER(IDER+I,IJ)=-DERQ
                 HDER(IDER+J,IJ)=DERQ
                 IDER=IDER+NBAT
 2504              CONTINUE    
              DO 2660 L=1,3
                 DI(L)=0.D0
 2660              DJ(L)=0.D0
              KKK=0  
C       
C       HAMILTONIEN PERTURBE SP (TERMES A 3 CORPS)
C       
              HIJ=0.D0
              DO 2600 K=1,NALK
                 IF((K.EQ.J).OR.(K.EQ.I)) GO TO 2605
                 IF(I.GT.K) THEN
                    ISI=1
                    IK=III+K
                 ELSE
                    IK=KKK+I
                    ISI=-1
                 ENDIF
                 IF(J.GT.K) THEN
                    ISJ=1
                    JK=JJJ+K
                 ELSE
                    JK=KKK+J
                    ISJ=-1
                 ENDIF
                 IKQ=IK+IK+IK-3
                 JKQ=JK+JK+JK-3
                 SCAL=UNIT(IKQ+1)*UNIT(JKQ+1)+UNIT(IKQ+2)*UNIT(JKQ+2)
     +               +UNIT(IKQ+3)*UNIT(JKQ+3)
                 SCAL=SCAL*ISI*ISJ
                 BIJK=VECT(IK)*VECT(JK)
                 HIJK=BIJK*SCAL
                 HIJ=HIJ-HIJK
                 RIJK=RINV(IK)*RINV(JK)
C       
C       GRADIENT DES TERMES SP
C       
                 IDER=0
                 DO 2610 L=1,3
                    IKQL=IKQ+L
                    JKQL=JKQ+L
                    UIK=ISI*UNIT(IKQL)
                    UJK=ISJ*UNIT(JKQL)
                    AIJK=BIJK*RIJK
                    DII=-DBSZ(IKQL)*VECT(JK)*ISI*SCAL
     1                    -ISJ*DIF(JKQL)*AIJK
     2                    +HIJK*UIK*RINV(IK)
                    DJJ=-DBSZ(JKQL)*VECT(IK)*ISJ*SCAL
     1                 -ISI*DIF(IKQL)*AIJK
     2           +HIJK*UJK*RINV(JK)
                    DK=-DII-DJJ
                    DI(L)=DI(L)+DII
                    DJ(L)=DJ(L)+DJJ
                    HDER(IDER+K,IJ)=-DK*DENOM
 2610            IDER=IDER+NBAT
 2605                 KKK=KKK+K      
 2600              CONTINUE
              
 3030              FORMAT( 3(2X,I3),3(2X,F10.6))
              H(IJ)=HSS+HIJ*DENOM
              IDER=0 
              DO 2670 L=1,3
                 HDER(IDER+I,IJ)=HDER(IDER+I,IJ)-DI(L)*DENOM
                 HDER(IDER+J,IJ)=HDER(IDER+J,IJ)-DJ(L)*DENOM
 2670         IDER=IDER+NBAT
              JJJ=JJJ+J
              IJQ=IJQ+3
 32           CONTINUE
C       
C--------------------------------------------------------
C       ELEMENTS DIAGONAUX (REPULSION)
C----------------------------------------------------------
C       
           IJ=IJ+1
           KKK=0
           HII=0.D0 
           DO 223 K=1,NALK
              IF(K.EQ.I) GO TO 221
              IF(I.GT.K) THEN
                 IK=III+K
                 ISI=1
              ELSE
                 IK=KKK+I
                 ISI=-1
              ENDIF
              IKQ=IK+IK+IK-3
C       
C       HAMILTONIEN
C       
              VAL=0.D0 
              DER=0.D0
              RIK=R(IK)
              IF (RIK.LE.RMAX) THEN
                 IF (RIK.GE.RMIN) THEN
                    CALL FRSS(RIK,VAL,DER)
                 ELSE
                    VAL=1.D0
                    DER=0.D0
                 ENDIF
                 HII=HII+VAL
              ENDIF
C       
C       GRADIENT
C       
              IDER=0
              DO 2584 L=1,3
                 IKQL=IKQ+L
                 DERQ=ISI*DER*UNIT(IKQL)
                 HDER(K+IDER,IJ)=DERQ
                 HDER(I+IDER,IJ)=HDER(I+IDER,IJ)-DERQ
                 IDER=IDER+NBAT
 2584              CONTINUE
 221              KKK=KKK+K
 223           CONTINUE
           H(IJ)=HII
           III=III+I
           IJQ=IJQ+3
              
 30        CONTINUE

C--------------------------------
C       DIAGONALISATION DE H
C---------------------------------
C       
C       NOMBRE D ORBITALES OCCUPEES
C       
           NDOC=NBLEC/2
           NOC=NDOC
           MONOC=.FALSE.
           IF ((2*NDOC).NE.NBLEC) THEN
              MONOC=.TRUE. 
              NOC=NDOC+1
           ENDIF
C
           DO I=1,NALK
              NUM(I)=I
           ENDDO
C       
C       SAUVEGARDE DE L HAMILTONIEN DANS HA
C       
           DO IJ=1,NALK*(NALK+1)/2
              HA(IJ)=H(IJ)
           ENDDO
           IF(NALK.EQ.1) THEN
              VALP(1)=-H(1)
              VECT(1)=1.D0
           ENDIF
           IF(NALK.GT.1) THEN
C       
C       APPEL A GIVENS
C       
              IF(IDIAG.EQ.1) THEN
C       
C       HA CONTIENT LA DEMI MATRICE INFERIEURE
C       H CONTIENT LA DEMI MATRICE SUPERIEURE
C       H EST DETRUIT APRES LA DIAGONALISATION
C       
                 II=0
                 IJ=0
                 DO I=1,NALK
                    JJ=II
                    DO J=I,NALK
                       IJ=IJ+1
                       H(IJ)=HA(JJ+I)
                       JJ=JJ+J
                    ENDDO
                    II=II+I
                 ENDDO

C START NEW DJW

                 HMAT(1:NBAT,1:NBAT)=0.0D0
                 II=0
                 DO I=1,NALK
                    DO J=I,NALK
                       II=II+1
                       HMAT(I,J)=H(II)
                    ENDDO
                 ENDDO
                 
                 ABSTOL=DLAMCH('SAFE  MINIMUM')
                 CALL DSYEVR('V','I','U',NALK,HMAT,NBAT,0.0D0,1.0D0,1,NOC,ABSTOL,NFOUND,VALP,
     1                        VT,NBAT,ISUPPZ,WORK,LWORK, IWORK, ILWORK, INFO )
                 IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' IN DSYEVR'

C END NEW DJW
                 
                 DO 210 I=1,NOC
                    VALP(I)=-VALP(I)
 210             CONTINUE
              ENDIF
           ENDIF
C       
           DO I=1,NOC
              VALP(I)=-VALP(I)
              II=(NUM(I)-1)*NALK
C             DO J=1,NALK           ! DJW COMMENT
C                VT(J,I)=VECT(II+J) ! DJW COMMENT
C             END DO                ! DJW COMMENT
           END DO
C       
C-----------------------------------------------------
C     ENERGIE TOTALE
C-----------------------------------------------------
           VGR=0.
           EN=VGR
           DO 401 I=1,NDOC
              EN=EN+VALP(I)*2.0
 401           CONTINUE
           IF (MONOC) EN=EN+VALP(NOC)
C       
C---------------------------------------------
C       DENSITE ELECTRONIQUE
C----------------------------------------------
           IND=0
           E0=EN
C-------------------------------------
C       CALCUL DU GRADIENT DE L ENERGIE
C-------------------------------------
C       
C       A LA FIN GRAD CONTIENT LES DERIVEES DANS L ORDRE
C       (DE/DXI,I=1,NBAT),(DE/DYI,I=1,NBAT),(DE/DZI,I=1,NBAT)
C       
           
           IF(NALK.GT.0) THEN
              DO 2570 K=1,3*NBAT
                 GRD=0.D0
                 DO 2572 I=1,NDOC
                    DO 2585 L=1,NALK
                       DO 2585 M=1,NALK
                          IF(L.GT.M) THEN
                             ML=L*(L-1)/2+M
                          ELSE
                             ML=M*(M-1)/2+L
                          ENDIF
                          GRD=GRD+2.D0*VT(L,I)*VT(M,I)*HDER(K,ML)
 2585            CONTINUE
 2572              CONTINUE
              IF(MONOC) THEN 
                 DO 2586 L=1,NALK
                    DO 2586 M=1,NALK
                       IF(L.GT.M) THEN
                          ML=L*(L-1)/2+M
                       ELSE
                          ML=M*(M-1)/2+L
                       ENDIF
 2586                 GRD=GRD+VT(L,NOC)*VT(M,NOC)*HDER(K,ML)
              ENDIF
 2570              GRAD(K)=GRD
           ENDIF
           
           IJ=0
           EREP=0.0
           DO I=1,NBAT
              DO J=1,I-1
                 IJ=IJ+1
                 EREP=EREP+1./(R(IJ)**12)
              ENDDO
              IJ=IJ+1
           ENDDO
           
           EREP=EREPMAX*EREP
           
           DO I=1,NBAT
              DP(3*I-2)=GRAD(I)
              DP(3*I-1)=GRAD(I+NBAT)
              DP(3*I  )=GRAD(I+2*NBAT)
           ENDDO
           
           DO I=1,NBAT
              DO J=1,NBAT
                 IF (I.NE.J) THEN
                    XIJ=P(3*I-2)-P(3*J-2)
                    YIJ=P(3*I-1)-P(3*J-1)
                    ZIJ=P(3*I)-P(3*J)
                    RIJ2=XIJ**2+YIJ**2+ZIJ**2
                    RIJ14=RIJ2**7
                    DP(3*I-2)=DP(3*I-2)-EREPMAX*12.*XIJ/RIJ14
                    DP(3*I-1)=DP(3*I-1)-EREPMAX*12.*YIJ/RIJ14
                    DP(3*I  )=DP(3*I  )-EREPMAX*12.*ZIJ/RIJ14
                 ENDIF
              ENDDO
           ENDDO
           
           E0=E0+EREP
           
        RETURN
        END

C-----------------------------------------------------
        SUBROUTINE FTSS(R,F,DF)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(10),B(10),C(10)
        CHARACTER*2 METAL
        COMMON/METAL/METAL
        IF (METAL.EQ.'NA') THEN
           A(1)=-0.0000600345
           A(2)=0.7366633997
           A(3)=-21.2536277450
           B(1)=7.8227026687
           B(2)=5.0784652072
           B(3)=-5.7014389921
           C(1)=1.4334140314
           C(2)=2.7260385637
           C(3)=0.D0
           NT=3
           F=0.D0
           DF=0.D0
           DO I=1,NT
              FI=A(I)*(R**B(I))*DEXP(-R*C(I))
              DFI=(-C(I)+B(I)/R)*FI
              F=F+FI
              DF=DF+DFI
           ENDDO
        ELSEIF (METAL.EQ.'AG') THEN
           A(1)=-.0148920093/27.21
           B(1)= 9.4301889740
           C(1)=2.2599701079
           A(2)=-11.8760873508/27.21
           B(2)=-1.0347495105
           C(2)= .5509616893
           NT=2
           F=0.D0
           DF=0.D0
           DO I=1,NT
              FI=A(I)*(R**B(I))*DEXP(-R*C(I))
              DFI=(-C(I)+B(I)/R)*FI
              F=F+FI
              DF=DF+DFI
           ENDDO
        ELSEIF (METAL.EQ.'LI') THEN
           AA=-0.000195D0
           BB=1.65D0
           XN=8.D0
           FEXP=AA*DEXP(-BB*R)
           F=(R**XN)*FEXP
           DF=(XN*R**(XN-1.D0)-BB*(R**XN))*FEXP
        ENDIF
        RETURN
        END
C-----------------------------------------------------
        SUBROUTINE FTSZ(R,F,DF)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(10),B(10),C(10)
        CHARACTER*2 METAL
        COMMON/METAL/METAL
        IF (METAL.EQ.'NA') THEN
           A(1)=0.0002818237
           A(2)=-0.0025316896
           A(3)=67.9894048289
           B(1)=9.2162096989
           B(2)= 4.0958902363
           B(3)= -5.2076863938
           C(1)= 2.2554055362
           C(2)=.9408906266
           C(3)=.6157151927
           NT=3
           F=0.D0
           DF=0.D0
           DO I=1,NT
              FI=A(I)*(R**B(I))*DEXP(-R*C(I))
              DFI=(-C(I)+B(I)/R)*FI
              F=F+FI
              DF=DF+DFI
           ENDDO
        ELSEIF (METAL.EQ.'AG') THEN
           NT=2
           A(1)= -.0058542/27.21
           B(1)= 9.5986911821
           C(1)= 2.0836237596
           A(2)= -17.9409106/27.219
           B(2)=-3.5128659875
           C(2)=-.0000000007
           A(1)=-.0043951845/27.21
           A(2)=-13.469637022/27.219
           F=0.D0
           DF=0.D0
           DO I=1,NT
              FI=A(I)*(R**B(I))*DEXP(-R*C(I))
              DFI=(-C(I)+B(I)/R)*FI
              F=F+FI
              DF=DF+DFI
           ENDDO
        ELSEIF (METAL.EQ.'LI') THEN
           AA=0.0000091D0
           BB=1.80D0
           XN=10.D0
           FEXP=AA*DEXP(-BB*R)
           F=(R**XN)*FEXP
           DF=(XN*R**(XN-1.D0)-BB*(R**XN))*FEXP
        ENDIF
        
        RETURN
        END
C-----------------------------------
        SUBROUTINE FRSS(R,F,DF)
        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER*2 METAL
        COMMON/METAL/METAL
        IF (METAL.EQ.'NA') THEN
           A=1.6822175236
           C= 1.3007271463
           RE=5.5425217957
           C6= 10.4276970395
           X6= 5.8999900152
           FREP=A*DEXP(-C*R)
           DREP=-C*FREP
           FVDW=C6*R**(-X6)
           DVDW=-X6*FVDW/R
           ARG= (RE/R-1.D0)
           IF(R.GT.(RE)) THEN
              COUP=1.D0
              DCOUP=0.D0
           ELSE
              COUP=DEXP(-ARG**2)
              DARG=-RE/R**2
              DCOUP=-2.D0*COUP*ARG*DARG
           ENDIF
           F=FREP-COUP*FVDW
           DF=DREP-DVDW*COUP-DCOUP*FVDW
        ELSEIF (METAL.EQ.'AG') THEN
           A=8245.404377/27.21
           C=2.292537
           RE=16.726372
           C6=10679.191179/27.21
           X6=6. 
           FREP=A*DEXP(-C*R)
           DREP=-C*FREP
           FVDW=C6*R**(-X6)
           DVDW=-X6*FVDW/R
           ARG= (RE/R-1.D0)
           IF(R.GT.(RE)) THEN
              COUP=1.D0
              DCOUP=0.D0
           ELSE
              COUP=DEXP(-ARG**2)
              DARG=-RE/R**2
              DCOUP=-2.D0*COUP*ARG*DARG
           ENDIF
           F=FREP-COUP*FVDW
           DF=DREP-DVDW*COUP-DCOUP*FVDW
        ELSEIF (METAL.EQ.'LI') THEN
           A=1.6822175236D0
           C= 1.3007271463D0*5.80D0/5.05D0
           RE=5.5425217957D0*5.05D0/5.80D0
           X6= 5.8999900152D0
           C6= 10.4276970395D0*(5.05D0/5.80D0)**X6
           FREP=A*DEXP(-C*R)
           DREP=-C*FREP
           FVDW=C6*R**(-X6)
           DVDW=-X6*FVDW/R
           ARG= (RE/R-1.D0)
           IF(R.GT.RE) THEN
              COUP=1.D0
              DCOUP=0.D0
           ELSE
              COUP=DEXP(-ARG**2)
              DARG=-RE/R**2
              DCOUP=-2.D0*COUP*ARG*DARG
           ENDIF
           F=1.43D0*(FREP-COUP*FVDW)
           DF=1.43D0*(DREP-DVDW*COUP-DCOUP*FVDW)
        ENDIF
        
        RETURN
        END
