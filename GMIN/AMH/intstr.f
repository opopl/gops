
C     --------------------- INTSTR ----------------------

      SUBROUTINE INTSTR

C     ---------------------------------------------------

C     INTSTR GENERATES INITIAL STRUCTURES AND PERFORMS
C            A QUICK ANALYSIS ON THE RG AND LR/SR

C     ---------------------------------------------------

      USE AMHGLOBALS,  ONLY: IRESRT,MOVANAL,MAXSIZ,NMRES,MAXPRO,
     *  NUMPRO,MAXCRD,NUMCRD,PRCORD,OARCHV,IRES,ISEED_AMH,QUENCH

      IMPLICIT NONE

C     INTERNAL VARIABLES:

         INTEGER PROCNT

C     REQUIRED SUBROUTINES

         EXTERNAL RNDCOL,RESTRT    !   RESTRT

C     --------------------- BEGIN -----------------------

C     PROCNT TRACKS THE NUMBER OF INITIAL STRUCTURES
C     GENERATED

      PROCNT=0

      IF( ( IRESRT.AND.(.NOT.MOVANAL) ) .OR. QUENCH )THEN

         CALL RESTRT(PROCNT,MAXSIZ,NMRES,MAXPRO,
     *               NUMPRO,MAXCRD,NUMCRD,PRCORD,
     *               OARCHV)

      ENDIF

      IF( PROCNT.LT.NUMPRO )THEN

C        IF NUMBER OF INITIAL CONFIGURATIONS FOUND
C        UP TO THIS POINT IS LESS THAN THE REQUIRED
C        NUMBER, THEN FIND RANDOM COILS FOR
C        THE REMAINDER
         CALL RNDCOL(NMRES,MAXPRO,NUMPRO,
     *                IRES,PROCNT,NUMCRD,PRCORD,
     *                ISEED_AMH)

      ENDIF


C     ---------------------- DONE -----------------------

      RETURN
      END
