 
C     --------------------- MDCTRL ----------------------
 
      SUBROUTINE MDCTRL(JSTRT,JFINS,TEMPAV,ISHKIT,
     *                  BDSHAK,TMPIDX,
     *                  TOTKE
     *                  )

C     ---------------------------------------------------

C     MDCTRL IS THE MASTER SUBROUTINE FOR THE MOLECULAR
C            DYNAMICS ALGORITHM OUTLINED IN RYCKAERT AND 
C            CICCOTTI, MOL. PHYS. 58, 1125-1136 (1986). 

C     ARGUMENTS:

C        JSTRT - FIRST NONFIXED SITE (I)
C        JFINS - LAST NONFIXED SITE (I)
C        TEMPAV- IF TEMPAV, THEN COMPUTE AVERAGES (I)
C        MAXSHK- MAXIMUM NUMBER OF SHAKE ITERATIONS (I)
C        TOLSHK- SHAKE TOLERANCE (I)
C        ISHKIT- TRACKS THE NUMBER OF SHAKE ITERATIONS (O)
C        BDSHK - TRUE IF SHAKE DOESN'T CONVERGE WITHIN
C                SPECIFIED NUMBER OF ITERATIONS (O)
C        TMPIDX- TEMPERATURE INDEX (I)

C     ---------------------------------------------------

      USE GLOBALS, ONLY: MAXTAB,NUMCRD,NUMPRO,QRCORD,WORK1,
     *  VELOCP,PRCORD,ZRCORD,AVEPP,MAXR,VPOTNT,
     *  FORSE,WORK4,PEXCLD,
     *  NUMLNG,NMRES,RINCINV,RINCSQ,OARCHV,ILONG,CRDIXN,
     *  IRES,WORK3,TIMSTP,EQDIST,HBOND,OXEXCLDV,
     *  I540,I1,
     *  AVEP,MOVANAL,TRCORD,SRCORD,BONDLN,TOLSHK,MAXSHK,MAXPRO,
     *  MAXCRD,I511,I518


      IMPLICIT NONE

C     ARGUMENT DECLARATIONS:

         LOGICAL TEMPAV,BDSHAK

         INTEGER JSTRT,JFINS,ISHKIT,TMPIDX

         REAL TOTKE(MAXPRO)
         

C     INTERNAL VARIABLES:

         LOGICAL SCL_CALL

         INTEGER I_AXIS, I_COORD, I_PRO, I_RES

         REAL STARHI,TRGENG(MAXTAB,3)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     --- EVERY STEP PART 1 OF 3 ---
C  NEEDED FOR PRINTING OUT EVERY STEP IN A TRAJECTORY
C
C         CHARACTER*10 NAME_A(10000)
C         CHARACTER*3 RES_TYPE(MAXSIZ)
C         CHARACTER*10 SAVE_NAME
C         INTEGER NAME_LENGTH, NL
C         INTEGER PDB, PDB_NUM, PN
C         EXTERNAL GET_RES_NAME
C         DATA PDB_NUM /0/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C     REQUIRED SUBROUTINES

         EXTERNAL FORCE,VERLET,
     *            SHKDRV


C     --------------------- BEGIN -----------------------

!         WRITE(SO,*) 'IN MDCTRL'

C     SAVE CURRENT COORDINATES

      DO 500 I_COORD=1, NUMCRD             ! CA, CA, O
         DO 524 I_PRO=1,NUMPRO
            DO 525 I_AXIS=1,3
               DO 526 I_RES=1,JFINS
                  QRCORD(I_RES,I_AXIS,I_PRO,I_COORD)=
     *            PRCORD(I_RES,I_AXIS,I_PRO,I_COORD)
  526          CONTINUE
  525       CONTINUE
  524    CONTINUE
  500 CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FIND FORCES AND POTENTIAL ENERGIES 
C       WRITE(SO,*) 'FORCE    '

      SCL_CALL=.FALSE.

      CALL FORCE(NUMPRO,
     *            PRCORD,
     *            ZRCORD,AVEPP,TEMPAV,
     *            MAXR,VPOTNT,
     *            FORSE,TRGENG,
     *            PEXCLD,
     *            NUMLNG,NMRES,RINCINV,RINCSQ,
     *            ILONG,
     *            CRDIXN,
     *            IRES,
     *            EQDIST,
     *            HBOND,OXEXCLDV,
     *            SCL_CALL)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         IF( TEMPAV ) THEN

C           RECORD TOTAL E FOR T AVERAGE

           DO 540 I540=1,NUMPRO
             DO 100 I1=1,50
               AVEP(I540,1,I1)=AVEPP(I540,1,I1)
100             CONTINUE
  540      CONTINUE

         ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     UPDATE COORDINATES

      IF (.NOT.MOVANAL) THEN

!       WRITE(SO,*) 'IN VERLET'

      CALL VERLET(NUMPRO,PRCORD,ZRCORD,
     *            TIMSTP,JSTRT,JFINS,TRCORD)

C     FIND SET OF FINAL COORDINATES CONSISTENT WITH
C     BOND CONSTRAINTS

            CALL SHKDRV(PRCORD,QRCORD,SRCORD,
     *            ZRCORD,BONDLN,JSTRT,JFINS,TOLSHK,
     *            MAXSHK,BDSHAK,NUMPRO,ISHKIT,
     *            MAXPRO,MAXCRD,NUMCRD,
     *            OARCHV,WORK1,WORK3,WORK4,IRES)

      ENDIF

C     IF BDSHAK, THEN SHAKE WAS INEFFECTIVE; RETURN

      IF( BDSHAK )THEN
         WRITE(OARCHV,144)
  144    FORMAT(/'BAD SHAKE POST VERLET')
         RETURN
      ENDIF
C     SET VELOCITY 

      STARHI=0.5/TIMSTP

      DO 514 I_COORD=1,NUMCRD
         DO 509 I_PRO=1,NUMPRO
            DO 510 I_AXIS=1,3
                  DO 511 I511=JSTRT,JFINS
                     VELOCP(I511,I_AXIS,I_PRO,I_COORD)=
     *              (PRCORD(I511,I_AXIS,I_PRO,I_COORD) - 
     *               TRCORD(I511,I_AXIS,I_PRO,I_COORD))*STARHI
  511             CONTINUE
               DO 518 I518=JSTRT,JFINS
                  TRCORD(I518,I_AXIS,I_PRO,I_COORD)=
     *            QRCORD(I518,I_AXIS,I_PRO,I_COORD)
  518          CONTINUE
  510       CONTINUE
  509    CONTINUE
  514 CONTINUE

      IF (TEMPAV) THEN
        TOTKE=0.0
        DO I_COORD=1,NUMCRD
        DO I_PRO=1,NUMPRO
        DO I_AXIS=1,3
        DO I511=JSTRT,JFINS
          TOTKE(I_PRO)=TOTKE(I_PRO)+VELOCP(I511,I_AXIS,I_PRO,I_COORD)**2
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        TOTKE=0.5*TOTKE
      ENDIF

C     ---------------------- DONE -----------------------

      RETURN
      END
