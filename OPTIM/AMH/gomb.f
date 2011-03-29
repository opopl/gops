
C     --------------------- GOMB ----------------------

      SUBROUTINE GOMB(  ETA,INDEX,
     *                  TEMPAV,
     *                  NGOMB,
     *                  A_TO_NMO,E,IBIASGAUSS,BIAS_AV,
     *                  BIAS_VAR,BIAS_PREFACTOR,BIAS_F,
     *                  IBIASPOLY,NBIASPOLY,BIASPOLY)

C     ---------------------------------------------------

C      GOMB (`GO MANY-BODY') EVALUATES THE A-FACTORS 
C      REQUIRED FOR THE GOMB IMPLEMENTATION
C
C     ---------------------------------------------------

      USE AMHGLOBALS,  ONLY: MAXCNT,AMHMAXSIZ,NMRES,NUMLNG,ILONG,
     *  FORSE,VPOTNT,RINCSQ,AMH_GSE,ERES,MAXR

      IMPLICIT NONE

       DOUBLE PRECISION, INTENT(IN):: ETA(1:MAXCNT,1:4),NGOMB,BIAS_AV,BIAS_VAR,
     *  BIAS_PREFACTOR,BIASPOLY(1:100)
       DOUBLE PRECISION, INTENT(OUT):: A_TO_NMO(1:AMHMAXSIZ,1:2),E(:,:),BIAS_F
      INTEGER, INTENT(IN):: INDEX(1:MAXCNT,1:4),NBIASPOLY
      LOGICAL, INTENT(IN):: TEMPAV,IBIASGAUSS,IBIASPOLY

C     INTERNAL VARIABLES:
         DOUBLE PRECISION A_TO_N(1:AMHMAXSIZ,1:2),M,C,E_ADDITION,ETA_PRIME
         DOUBLE PRECISION AGOMB1(1:AMHMAXSIZ,1:4),AGOMB2(1:AMHMAXSIZ,1:4),Q
         INTEGER I_IXN,ISIT1,ISIT2,I,J,I_TAB,I_RES,INDX

C     --------------------- BEGIN -----------------------

!     ZERO ENERGY AND FORCE FACTOR

      E=0.0D0
      BIAS_F=0.0D0

C     INITIALISE --- I LABELS RESIDUE, J TABLE

      DO 504 I=1,NMRES
        DO 501 J=1,4
          AGOMB1(I,J)=0.0D0
          AGOMB2(I,J)=0.0D0
501     ENDDO
504   ENDDO

         DO 506 I_TAB=1,4
           DO 505 I_IXN=1,NUMLNG(NMRES,I_TAB)
C  WANT TO WORK OUT AGOMB(I)=SUM_J(THETA_IJ) WHERE THETA_IJ
C  ARE THE PAIR ENERGIES. THIS IS SLIGHTLY COMPLICATED BY
C  THE FACT THAT I LABELS AN ATOM (CA OR CB) AND THINGS
C  ARE ANNOYINGLY STORED IN 4 TABLES. AGOMB IS THE QUANTITY
C  THAT WE TAKE MOD OF AND ^(NGOMB-1) IN THE I509 LOOP. THE
C  SECOND INDEX IN A_TO_NMO LABELS ALPHA/BETA (1/2)

             ISIT1=ILONG(I_IXN,1,I_TAB)
             ISIT2=ILONG(I_IXN,2,I_TAB)

             INDX=INDEX(I_IXN,I_TAB)
             IF (INDX.GT.MAXR) CYCLE !OUTSIDE RANGE OF FORCE TABLE

C  FOR ACCURACY CHOOSE INTERPOLATE POTENTIAL SO IT EXACTLY 
C  (HOPEFULLY) CORRESPONDS TO INTEGRATED FORCE. THIS
C  NATURALLY DEPENDS ON METHOD OF INTERPOLATING THE FORCE
         
             M=     (FORSE(INDX+1,I_IXN,I_TAB)-
     *                       FORSE(INDX,I_IXN,I_TAB))

             C=FORSE(INDX,I_IXN,I_TAB)
     *                                *FLOAT(INDX+1)
     *           -FORSE(INDX+1,I_IXN,I_TAB)
     *                                   *FLOAT(INDX)

             ETA_PRIME=1.0D0-ETA(I_IXN,I_TAB) 

             E_ADDITION=VPOTNT(INDX+1,I_IXN,I_TAB)+
     *      RINCSQ(I_IXN,I_TAB)*
     *      ETA_PRIME* (M*
     *      (FLOAT(INDX+1)*(FLOAT(INDX+1)-ETA_PRIME)+
     *       ETA_PRIME*ETA_PRIME/3.0D0)  +
     *      C* (FLOAT(INDX+1)-0.5D0*ETA_PRIME)  )

             AGOMB1(ISIT1,I_TAB)=AGOMB1(ISIT1,I_TAB) + 
     *         E_ADDITION

             AGOMB2(ISIT2,I_TAB)=AGOMB2(ISIT2,I_TAB) + 
     *         E_ADDITION          

  505     CONTINUE
506     CONTINUE

       DO 509 I_RES=1,NMRES
        A_TO_NMO(I_RES,1)  = ABS( ( AGOMB1(I_RES,1) + AGOMB1(I_RES,2)
     +        + AGOMB2(I_RES,1) + AGOMB2(I_RES,3) ) )**(NGOMB-1.0)
        A_TO_NMO(I_RES,2)  = ABS( ( AGOMB1(I_RES,3) + AGOMB1(I_RES,4)
     +        + AGOMB2(I_RES,2) + AGOMB2(I_RES,4) ) )**(NGOMB-1.0)
509    CONTINUE


      IF (TEMPAV.OR.IBIASGAUSS.OR.IBIASPOLY) THEN
        DO 510 I_RES=1,NMRES
          A_TO_N(I_RES,1)  =  A_TO_NMO(I_RES,1) * ( AGOMB1(I_RES,1)
     +      + AGOMB1(I_RES,2) + AGOMB2(I_RES,1) + AGOMB2(I_RES,3) )
          A_TO_N(I_RES,2)  =  A_TO_NMO(I_RES,2) * ( AGOMB1(I_RES,3)
     +      + AGOMB1(I_RES,4) + AGOMB2(I_RES,2) + AGOMB2(I_RES,4) )
        E(1,5)=E(1,5) + A_TO_N(I_RES,1) + A_TO_N(I_RES,2)
!       WRITE(SO,*) I_RES,A_TO_N(I_RES,1) + A_TO_N(I_RES,2)
510     CONTINUE
        E(1,5)=E(1,5)*0.5D0
        AMH_GSE=-4.0D0*ERES*FLOAT(NMRES)

C    THESE BELOW ARE BACKWARDS WAYS OF ADDING AN
C    UNMBRELLING POTENTIAL

        IF (IBIASGAUSS) THEN
          E(1,10)=BIAS_PREFACTOR*EXP(-(0.5D0/BIAS_VAR)*
     *                    (E(1,5)/AMH_GSE-BIAS_AV)**2)
          BIAS_F=-(E(1,10)/BIAS_VAR)*
     *                ((E(1,5)/AMH_GSE)-BIAS_AV)/AMH_GSE
        ELSEIF (IBIASPOLY) THEN
          Q=E(1,5)/AMH_GSE
          E(1,10)=0.0D0
          DO I=1,NBIASPOLY
            E(1,10)=E(1,10)+BIASPOLY(I)*Q**I
            BIAS_F=BIAS_F+(1.0D0/AMH_GSE)*(FLOAT(I))
     *                             *BIASPOLY(I)*Q**(I-1) 
          ENDDO
        ENDIF
      ENDIF

C     ---------------------- DONE -----------------------

      RETURN
      END
