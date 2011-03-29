C     ------------------- RAMA ----------------------
      SUBROUTINE RAMA(IRES,JSTRT,JFINS,PRO_CORD,F_CORD,
     *                RAMASCL,RAMAPOT,NITCORD,CPRCORD)

C     --------------------------------------------------

C     RAMA FIGURES OUT CHIRAL FORCES 

      USE AMHGLOBALS,  ONLY:MAXSIZ,MAXCRD,APS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: PRO_CORD(MAXSIZ,3,MAXCRD),RAMASCL,
     *          NITCORD(MAXSIZ,3),CPRCORD(MAXSIZ,3)
      DOUBLE PRECISION, INTENT(OUT) :: F_CORD(MAXSIZ,3,MAXCRD),RAMAPOT(MAXSIZ)
      INTEGER, INTENT(IN):: IRES(MAXSIZ),JSTRT,JFINS
 
C
C
C        A        VECTOR FROM CA (PRO_CORD) TO C'
C        B        VECTOR FROM N TO CA
C        C        VECTOR FROM C' TO N
C
C        ADB        A DOT B
C        BDC        B DOT C
C        ADC        A DOT C
C
C        BXA        B CROSS A
C
C
C
C     SIGN : RETURNS THE MAGNITUDE OF ITS FIRST ARGUMENT WITH THE SIGN OF ITS
C            SECOND ARGUMENT.
C
C     ARGUMENT DECLARATIONS:

             DOUBLE PRECISION 
     *        A(3),B(3),C(3),PA(3),PB(3),PC(3),
     *              PA2,PB2,PC2,PADB,PADC,PBDC,PBXA(3),
     *              PCDBXA,PCHI,PQUAD,PFAC,
     *        ADIR(0:MAXSIZ+1,3,2),
     *        BDIR(0:MAXSIZ+1,3,2),W(6),
     *        CDIR(0:MAXSIZ+1,3,2),
     *        A2,B2,C2,ADB,ADC,BDC,BXA(3),CDBXA,
     *        CHI,QUAD,PHI(MAXSIZ),PSI(MAXSIZ),
     *        SLOPE(2),PROB(MAXSIZ),
     *        PKAPPA1,PKAPPA2,
     *        FAC,KAPPA1,KAPPA2,KAPPA,PKAPPA,
     *        CVAL1,CVAL2,CVAL3,NVAL1,NVAL2,NVAL3

C        VARIABLES FOR NEW POTENTIAL

            DOUBLE PRECISION PHID(6),PSID(6),F(6),SIGMA(6),S(6)
        INTEGER IPS

C     INTERNAL VARIABLES:
C        --- DO LOOP INDICES ---
         INTEGER I_AXIS,I_RES,I1,I2,I3

             DOUBLE PRECISION PI

        DATA PI /3.141592653589793238462643383279502884197D0/
        DATA W /1.3149D0,1.17016D0,1.29264D0,1.78596D0,1.0D0,1.0D0/
        DATA SIGMA /15.398D0,100.521D0,49.0954D0,419.123D0,0.0D0,0.0D0/
        DATA PHID /-2.051D0,-1.353D0,-1.265D0,-1.265D0,0.0D0,0.0D0/
        DATA PSID /2.138D0,2.4D0,-0.218D0,-0.929D0,0.0D0,0.0D0/

         CVAL1=0.444D0
         CVAL2=0.235D0
         CVAL3=0.321D0

         NVAL1=0.483D0
         NVAL2=0.703D0
         NVAL3=-0.186D0

C     --------------------- BEGIN -----------------------
C     TESTING FOR CHRIS

!  ZERO FORCE AND ENERGY
      RAMAPOT(:)=0.0D0
      F_CORD=0.0D0

        DO 600 I1=1,MAXSIZ
          PHI(I1)=0.0D0
          PSI(I1)=0.0D0
600        CONTINUE

        DO 610 I1=0,MAXSIZ
          DO 620 I2=1,3
            DO 630 I3=1,2
              ADIR(I1,I2,I3)=0.0D0
              BDIR(I1,I2,I3)=0.0D0
              CDIR(I1,I2,I3)=0.0D0
630            CONTINUE
620          CONTINUE
610        CONTINUE
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  PHI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DO 501 I_RES=JSTRT+1,JFINS-1
C       CALCULATE VECTORS        

            PA(1)= CPRCORD(I_RES,1)-PRO_CORD(I_RES,1,1)
            PA(2)= CPRCORD(I_RES,2)-PRO_CORD(I_RES,2,1)
            PA(3)= CPRCORD(I_RES,3)-PRO_CORD(I_RES,3,1)

            PB(1)= PRO_CORD(I_RES,1,1)-NITCORD(I_RES,1)
            PB(2)= PRO_CORD(I_RES,2,1)-NITCORD(I_RES,2)
            PB(3)= PRO_CORD(I_RES,3,1)-NITCORD(I_RES,3)

            PC(1)= NITCORD(I_RES,1)-CPRCORD(I_RES-1,1)
            PC(2)= NITCORD(I_RES,2)-CPRCORD(I_RES-1,2)
            PC(3)= NITCORD(I_RES,3)-CPRCORD(I_RES-1,3)

C     CALCULATE DOT PRODUCTS 

        PADB=PA(1)*PB(1)+PA(2)*PB(2)+PA(3)*PB(3)
        PBDC=PB(1)*PC(1)+PB(2)*PC(2)+PB(3)*PC(3)
        PADC=PA(1)*PC(1)+PA(2)*PC(2)+PA(3)*PC(3)

        PBXA(1)=PB(2)*PA(3)-PB(3)*PA(2)
        PBXA(2)=PB(3)*PA(1)-PB(1)*PA(3)
        PBXA(3)=PB(1)*PA(2)-PB(2)*PA(1)

        PCDBXA=PC(1)*PBXA(1)+PC(2)*PBXA(2)+PC(3)*PBXA(3)

        PA2=PA(1)*PA(1)+PA(2)*PA(2)+PA(3)*PA(3)
        PB2=PB(1)*PB(1)+PB(2)*PB(2)+PB(3)*PB(3)
        PC2=PC(1)*PC(1)+PC(2)*PC(2)+PC(3)*PC(3)

        PKAPPA1=DSQRT( PA2*PB2-PADB*PADB )
        PKAPPA2=DSQRT( PB2*PC2-PBDC*PBDC )
        PKAPPA=PKAPPA1*PKAPPA2

        PCHI=(PADB*PBDC-PADC*PB2)/PKAPPA
C        PCHI=DMAX1(DMIN1(PCHI,0.99999999D9),-0.99999999D9)

C        IF (PCHI.GT.1D0) PCHI=1D0
C        IF (PCHI.LT.-1D0) PCHI=-1D0

        IF (PCHI.GE.1.0D0)PCHI=0.9999
        IF (PCHI.LE.-1.0D0)PCHI=-0.9999

        PQUAD=DSIGN(1.0D0,PCDBXA)
        PHI(I_RES)=PQUAD*DACOS(PCHI)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  PSI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

            A(1)=  NITCORD(I_RES+1,1)-CPRCORD(I_RES,1)
            A(2)=  NITCORD(I_RES+1,2)-CPRCORD(I_RES,2)
            A(3)=  NITCORD(I_RES+1,3)-CPRCORD(I_RES,3)

            B(1)=  CPRCORD(I_RES,1)-PRO_CORD(I_RES,1,1)
            B(2)=  CPRCORD(I_RES,2)-PRO_CORD(I_RES,2,1)
            B(3)=  CPRCORD(I_RES,3)-PRO_CORD(I_RES,3,1)

            C(1)=  PRO_CORD(I_RES,1,1)-NITCORD(I_RES,1)
            C(2)=  PRO_CORD(I_RES,2,1)-NITCORD(I_RES,2)
            C(3)=  PRO_CORD(I_RES,3,1)-NITCORD(I_RES,3)

C     CALCULATE DOT PRODUCTS

        ADB=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
        BDC=B(1)*C(1)+B(2)*C(2)+B(3)*C(3)
        ADC=A(1)*C(1)+A(2)*C(2)+A(3)*C(3)

        BXA(1)=B(2)*A(3)-B(3)*A(2)
        BXA(2)=B(3)*A(1)-B(1)*A(3)
        BXA(3)=B(1)*A(2)-B(2)*A(1)

        CDBXA=C(1)*BXA(1)+C(2)*BXA(2)+C(3)*BXA(3)

        A2=A(1)*A(1)+A(2)*A(2)+A(3)*A(3)
        B2=B(1)*B(1)+B(2)*B(2)+B(3)*B(3)
        C2=C(1)*C(1)+C(2)*C(2)+C(3)*C(3)

        KAPPA1=DSQRT( A2*B2-ADB*ADB )
        KAPPA2=DSQRT( B2*C2-BDC*BDC )
        KAPPA=KAPPA1*KAPPA2

        CHI=(ADB*BDC-ADC*B2)/KAPPA
C        CHI=DMAX1(DMIN1(CHI,0.999999999999999D9),-0.999999999999999999D9)

C        IF (CHI.GT.1D0) CHI=1D0
C        IF (CHI.LT.-1D0) CHI=-1D0

        IF (CHI.GE.1.0D0)CHI=0.9999
        IF (CHI.LE.-1.0D0)CHI=-0.9999

        QUAD=DSIGN(1.0D0,CDBXA)
        PSI(I_RES)=QUAD*DACOS(CHI)

        IF (PSI(I_RES).GE. PI) PSI(I_RES)=-PSI(I_RES)
        IF (PHI(I_RES).GE. PI) PHI(I_RES)=-PHI(I_RES)
        IF (IRES(I_RES).NE.8) THEN
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEW POTENTIAL%%%%%%%%%%%%%%%%%%

        PROB(I_RES) = 0.0D0

        DO  IPS = 1,6
          F(IPS) = DEXP(-SIGMA(IPS)*(DCOS(PHI(I_RES) -
     *     PHID(IPS)) -1.0D0)**2)
           F(IPS) = F(IPS)
           S(IPS) = DEXP(-SIGMA(IPS)*(DCOS(PSI(I_RES) -
     *     PSID(IPS)) -1.0D0)**2)
           S(IPS) = S(IPS)
        PROB(I_RES) = PROB(I_RES) + W(IPS)*F(IPS)*S(IPS)*APS(I_RES,IPS)

        ENDDO   ! IPS

        RAMAPOT(I_RES) = - RAMASCL*PROB(I_RES)
        SLOPE(1) = 0.0D0
        SLOPE(2) = 0.0D0
        DO  IPS = 1,6
           SLOPE(1) = SLOPE(1) +APS(I_RES,IPS)*
     *     W(IPS)*RAMASCL*2.0D0*SIGMA(IPS)*F(IPS)*S(IPS)*
     *      (DCOS(PHI(I_RES)-PHID(IPS)) -1.0D0)*DSIN(PHI(I_RES) - PHID(IPS))
           SLOPE(2) = SLOPE(2) + APS(I_RES,IPS)*
     *     W(IPS)*RAMASCL*2.0D0*SIGMA(IPS)*F(IPS)*S(IPS)*
     *      (DCOS(PSI(I_RES)-PSID(IPS)) -1.0D0)*DSIN(PSI(I_RES) - PSID(IPS))

        ENDDO  ! IPS

        PFAC=-PQUAD*SLOPE(1)/DSQRT(1.0D0-PCHI*PCHI)
       FAC=-QUAD*SLOPE(2)/DSQRT(1.0D0-CHI*CHI)

        ADIR(I_RES,1,1)=PFAC*((PBDC*PB(1)-PB2*PC(1))/PKAPPA
     *         -PCHI*(PB2*PA(1)-PADB*PB(1))/PKAPPA1**2)
        ADIR(I_RES,2,1)=PFAC*((PBDC*PB(2)-PB2*PC(2))/PKAPPA
     *         -PCHI*(PB2*PA(2)-PADB*PB(2))/PKAPPA1**2)
        ADIR(I_RES,3,1)=PFAC*((PBDC*PB(3)-PB2*PC(3))/PKAPPA
     *         -PCHI*(PB2*PA(3)-PADB*PB(3))/PKAPPA1**2)

        BDIR(I_RES,1,1)=PFAC
     *         *((PADB*PC(1)+PBDC*PA(1)-2.0D0*PADC*PB(1))/PKAPPA
     *         -PCHI*(PA2*PB(1)-PADB*PA(1))/PKAPPA1**2
     *         -PCHI*(PC2*PB(1)-PBDC*PC(1))/PKAPPA2**2)
        BDIR(I_RES,2,1)=PFAC
     *         *((PADB*PC(2)+PBDC*PA(2)-2.0D0*PADC*PB(2))/PKAPPA
     *         -PCHI*(PA2*PB(2)-PADB*PA(2))/PKAPPA1**2
     *         -PCHI*(PC2*PB(2)-PBDC*PC(2))/PKAPPA2**2)
        BDIR(I_RES,3,1)=PFAC
     *         *((PADB*PC(3)+PBDC*PA(3)-2.0D0*PADC*PB(3))/PKAPPA
     *         -PCHI*(PA2*PB(3)-PADB*PA(3))/PKAPPA1**2
     *         -PCHI*(PC2*PB(3)-PBDC*PC(3))/PKAPPA2**2)

        CDIR(I_RES,1,1)=PFAC*((PADB*PB(1)-PB2*PA(1))/PKAPPA
     *         -PCHI*(PB2*PC(1)-PBDC*PB(1))/PKAPPA2**2)
        CDIR(I_RES,2,1)=PFAC*((PADB*PB(2)-PB2*PA(2))/PKAPPA
     *         -PCHI*(PB2*PC(2)-PBDC*PB(2))/PKAPPA2**2)
        CDIR(I_RES,3,1)=PFAC*((PADB*PB(3)-PB2*PA(3))/PKAPPA
     *         -PCHI*(PB2*PC(3)-PBDC*PB(3))/PKAPPA2**2)

        ADIR(I_RES,1,2)=FAC*((BDC*B(1)-B2*C(1))/KAPPA
     *         -CHI*(B2*A(1)-ADB*B(1))/KAPPA1**2)
        ADIR(I_RES,2,2)=FAC*((BDC*B(2)-B2*C(2))/KAPPA
     *         -CHI*(B2*A(2)-ADB*B(2))/KAPPA1**2)
        ADIR(I_RES,3,2)=FAC*((BDC*B(3)-B2*C(3))/KAPPA
     *         -CHI*(B2*A(3)-ADB*B(3))/KAPPA1**2)

        BDIR(I_RES,1,2)=FAC
     *         *((ADB*C(1)+BDC*A(1)-2.0D0*ADC*B(1))/KAPPA
     *         -CHI*(A2*B(1)-ADB*A(1))/KAPPA1**2
     *         -CHI*(C2*B(1)-BDC*C(1))/KAPPA2**2)
        BDIR(I_RES,2,2)=FAC
     *         *((ADB*C(2)+BDC*A(2)-2.0D0*ADC*B(2))/KAPPA
     *         -CHI*(A2*B(2)-ADB*A(2))/KAPPA1**2
     *         -CHI*(C2*B(2)-BDC*C(2))/KAPPA2**2)
        BDIR(I_RES,3,2)=FAC
     *         *((ADB*C(3)+BDC*A(3)-2.0D0*ADC*B(3))/KAPPA
     *         -CHI*(A2*B(3)-ADB*A(3))/KAPPA1**2
     *         -CHI*(C2*B(3)-BDC*C(3))/KAPPA2**2)

        CDIR(I_RES,1,2)=FAC*((ADB*B(1)-B2*A(1))/KAPPA
     *         -CHI*(B2*C(1)-BDC*B(1))/KAPPA2**2)
        CDIR(I_RES,2,2)=FAC*((ADB*B(2)-B2*A(2))/KAPPA
     *         -CHI*(B2*C(2)-BDC*B(2))/KAPPA2**2)
        CDIR(I_RES,3,2)=FAC*((ADB*B(3)-B2*A(3))/KAPPA
     *         -CHI*(B2*C(3)-BDC*B(3))/KAPPA2**2)

             ENDIF
501        CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       COMPUTE FORCES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DO 503 I_AXIS=1,3
          DO 504 I_RES=JSTRT,JFINS

      F_CORD(I_RES,I_AXIS,1)=F_CORD(I_RES,I_AXIS,1)
     * +CVAL2*ADIR(I_RES-1,I_AXIS,1)
     * -(CVAL2+CVAL3)*ADIR(I_RES,I_AXIS,1)
     * +(NVAL1+NVAL3)*BDIR(I_RES,I_AXIS,1)
     * -(NVAL1)*BDIR(I_RES+1,I_AXIS,1)
     * +(NVAL2-CVAL2)*CDIR(I_RES,I_AXIS,1)
     * +(NVAL1-CVAL1)*CDIR(I_RES+1,I_AXIS,1)
     * +(NVAL2-CVAL2)*ADIR(I_RES-1,I_AXIS,2)
     * +(NVAL1-CVAL1)*ADIR(I_RES,I_AXIS,2)
     * +CVAL2*BDIR(I_RES-1,I_AXIS,2)
     * -(CVAL2+CVAL3)*BDIR(I_RES,I_AXIS,2)
     * +(NVAL1+NVAL3)*CDIR(I_RES,I_AXIS,2)
     * -NVAL1*CDIR(I_RES+1,I_AXIS,2)

       F_CORD(I_RES,I_AXIS,3)=F_CORD(I_RES,I_AXIS,3)
     * +CVAL3*ADIR(I_RES,I_AXIS,1)
     * +(-NVAL3)*BDIR(I_RES+1,I_AXIS,1)
     * +(-CVAL3+NVAL3)*CDIR(I_RES+1,I_AXIS,1)
     * +(-CVAL3+NVAL3)*ADIR(I_RES,I_AXIS,2)
     * +CVAL3*BDIR(I_RES,I_AXIS,2)
     * +(-NVAL3)*CDIR(I_RES+1,I_AXIS,2)

504      CONTINUE
503    CONTINUE

C     ---------------------- DONE -----------------------

      RETURN
      END
