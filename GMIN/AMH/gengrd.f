
C     --------------------- GENGRD ----------------------

      SUBROUTINE GENGRD(MAXCNT,ILONG,NUMLNG,NMRES,MAXS,MAXRIJ,XWORK,DELTA,DELTE,
     *                  MAXSIZ,DELTZ,RINC,IDIGNS,OARCHV,RSEP,RCUTAMH)

C     ---------------------------------------------------

C     GENGRD GENERATES THE R-GRID

C     ARGUMENTS:

C        MAXCNT- MAXIMUM NUMBER OF INTERACTIONS (I)
C        ILONG - LIST OF INTERACTING SITES (I)
C        NUMLNG- BREAKDOWN OF STRUCTURE OF ILONG (I)
C        NMRES - NUMBER OF TARGET SITES (I)
C        MAXS  - NUMBER OF R-GRID POINTS (I)
C        MAXRIJ- MAXIMUM GRID VALUE FOR EACH  (W)
C        XWORK - WORK ARRAY (W)
C        DELTA - GAUSSIAN WELL-WIDTH (I)
C        DELTE - EXPONENT FOR SCALING OF GAUSSIAN AS 
C                FUNCTION OF I-J (I)
C        MAXSIZ- MAXIMUM PROTEIN LENGTH (I)
C        DELTZ - DENOMINATOR IN GAUSSIAN FOR EACH INTERACTION
C                PAIR (O)
C        RINC  - R-GRID INCREMENT FOR EACH INTERACTION PAIR (O)


C       INDXT  - TOTAL NUMBER OF INTERACTIONS.
C       IDIFF  - DISTANCE IN NUMBER OF RESIDUES
C
C     ---------------------------------------------------

      USE AMHGLOBALS,  ONLY: N_LETTERS

C     SET REQUIRED PARAMETERS

      IMPLICIT NONE

C     ARGUMENT DECLARATIONS:

         LOGICAL IDIGNS

         INTEGER MAXSIZ,MAXCNT,ILONG(MAXCNT,2),NUMLNG(0:MAXSIZ),NMRES,MAXS,OARCHV

         DOUBLE PRECISION DELTZ(MAXCNT),RINC(MAXCNT),XWORK(MAXCNT),
     *        MAXRIJ(MAXCNT),DELTA,DELTE,RSEP,RCUTAMH
     
C     INTERNAL VARIABLES:

         INTEGER INDXT,IDIFF

C        --- DO LOOP INDICES ---

         INTEGER I_INDX
 
C        --- IMPLIED DO LOOP INDICES ---


         DOUBLE PRECISION DELTC,SCALR


C     --------------------- BEGIN -----------------------

C     --- DIAGNOSTICS ---

C      WRITE(OARCHV,145)MAXCNT,NUMLNG(NMRES),NMRES,MAXS,
C     *                 MAXSIZ,DELTA,DELTE
C  145 FORMAT(/'GENGRD: MAXCNT ',I5,' NUMLNG ',I5,
C     *        ' NMRES ',I3,' MAXS ',I4,' MAXSIZ ',I4,
C     *        ' DELTA ',1PE10.3,' DELTE ',1PE10.3)

C     --- END DIAGNOSTICS ---

C     INITIALIZE NUMBER OF INTERACTIONS

      INDXT=NUMLNG(NMRES)
      NUMLNG(0)=0

C     INITIALIZE MAXIMUM R_{IJ} 

      DO 500 I_INDX=1,INDXT
         IDIFF=ILONG(I_INDX,2) - ILONG(I_INDX,1)    ! DISTANCE IN SEQUENCE SPACE
         MAXRIJ(I_INDX)=MIN( RSEP*FLOAT(IDIFF), RCUTAMH )
  500 CONTINUE
C     SET ENDPOINTS SO THAT GAUSSIAN POTENTIAL IS CONTINUOUS,
C     I.E., V AND F AT ENDPOINTS ARE 'SMALL'

      DELTC=0.5D0/DELTA**2
      SCALR=(-DELTE*2)

      DO 501 I_INDX=1,INDXT
         DELTZ(I_INDX)=DELTC*(FLOAT(ILONG(I_INDX,2) - ILONG(I_INDX,1) ))**SCALR
        IF ( (ILONG(I_INDX,2) - ILONG(I_INDX,1) .GT. 4 ) .AND.
     *  (ILONG(I_INDX,2) - ILONG(I_INDX,1) .LT. 13).AND.N_LETTERS.EQ.4) THEN

         DELTZ(I_INDX)=2.0D0*(FLOAT(ILONG(I_INDX,2) - ILONG(I_INDX,1) ))**(-0.60D0) 
        ENDIF

  501 CONTINUE

C     SET UPPER BOUND FOR R-GRID

      DELTC=-LOG(1.0E-08)

C     SET ADDEND FOR DETERMINING MAXIMUM R-GRID VALUE

      DO 502 I_INDX=1,INDXT
         XWORK(I_INDX)=SQRT(DELTC/DELTZ(I_INDX))
  502 CONTINUE
C     SET ENDPOINT FOR EACH INTERACTION PAIR

      DO 503 I_INDX=1,INDXT
         MAXRIJ(I_INDX)=MAXRIJ(I_INDX) + XWORK(I_INDX)
  503 CONTINUE
C     INITIALIZE GRID INCREMENTS BASED ON MIN AND MAX
C     DISTANCE VALUES FOR EACH CONSTRAINT

      DO 504 I_INDX=1,INDXT
         RINC(I_INDX)=MAXRIJ(I_INDX)/FLOAT(MAXS-1)
  504 CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     --- DIAGNOSTICS ---
      IF( IDIGNS )THEN
C        LIST SPANNING INTERVAL FOR EACH | I - J |
C         WRITE(OARCHV,113)INDXT
C  113    FORMAT(/'B1 AND DELTZ:INDX ',I8)
C         DO 517 I517=1,MIN( 53,NMRES )
CC            IF( (ILONG(I517,1).EQ.4).AND.
CC     *          (ILONG(I517,2).EQ.88) )THEN 
CC            DO 518 I518=NUMLNG(I517-1)+1,
CC     *              MIN(NUMLNG(I517-1)+10,NUMLNG(I517))
C             DO 519 I519=1,2
C               IF( I519.EQ.1 )THEN
C                  I518=NUMLNG(I517-1)+1
C               ELSE
C                  I518=NUMLNG(I517)
C               ENDIF
C               WRITE(OARCHV,114)ILONG(I518,1),ILONG(I518,2),
C     *                          MAXRIJ(I518)-XWORK(I518),
C     *                          XWORK(I518),DELTZ(I518),
C     *                          RINC(I518)
C  114          FORMAT(2(I3,1X),4(1X,1PE10.3))
C  519       CONTINUE
CC  518       CONTINUE
CC            ENDIF
C  517    CONTINUE
      ENDIF
C     --- END DIAGNOSTICS ---
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     ---------------------- DONE -----------------------

      RETURN
      END
