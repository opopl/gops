
C     --------------------- GENTAB ----------------------

      SUBROUTINE GENTAB

C     ---------------------------------------------------

C     GENTAB IS DRIVER ROUTINE FOR GENERATING REQUIRED TABLES

C     ---------------------------------------------------
C
C  SUBROUTINE LOGIC
C
C     * CALL TIME
C     * LOOP OVER NUMBER OF TABLES
C       * DETERMINE THE INTERACTING SUBSET OF RESIDUES
C         [ THIS MAY BE HISTORICAL ]
C       * CALL POTIND
C
C       * CALL GENGRD
C 
C       * IF TABLE = 1, PASSI = TRUE
C       * IF TABLE = NUMTAB, PASSF = TRUE
C
C       * CALL GASPOT
C
C     * CALL TIME
C     * CONSTRUCT HBOND TABLES
C
C  SUBROUTINE VARIABLES
C
C      RSEP

      USE AMHGLOBALS,  ONLY:SO,NUMTAB,CRDIXN,NMRES,
     *  IDIGNS,MAXSIZ,MAXCNT,NUMLNG,ILONG,OARCHV,MAXS,WORK6,WORK3,
     *  DELTA,DELTE,DELTZ,RINC,RCUTAMH,TARGET,BONDLN,NUMCRD,
     *  N_DIVS_MAX,TARG_DIST,MAXCRD,IBIASFILE,I_REP,MAXTAB,ISS_STRUCT,
     *  I_QBIAS_A,WIDTH_QEXP_A,DEL_R_A,Q_IJ_A,DQ_DR_IJ_A,I_BIAS_NATIVE_A,
     *  I_QBIAS_B,I_BIAS_NATIVE_B,TOTALSSNORM,NUMCONST_A, NUMCONST_B,
     *  DEL_R_B,Q_IJ_B,DQ_DR_IJ_B,WIDTH_QEXP_B,SS_A,SS_B,SS_DIST

      USE AMH_INTERFACES, ONLY:REP_CONTACT

      IMPLICIT NONE

C     INTERNAL VARIABLES:

         LOGICAL PASSI,PASSF

         INTEGER JBEGN,JEND, I_TAB, I_DIST,I_DIFF,TAB,I_IXN,ISIT1,ISIT2,NMRES_CHECK,IDUMMY,
     *           I_RES,I_RES1, I_RES2

         DOUBLE PRECISION RSEP, MAX_CONT,MAX_RIJ,
     *         SMALL,XDIFF,YDIFF,ZDIFF,
     *         BIAS_CORD(MAXSIZ,3,MAXCRD),
     *         XDIFFTEMP,YDIFFTEMP,ZDIFFTEMP

        INTEGER TEMPNRES,SEQ_A(500),SEQ_B(500),ICOORD,SA_ALPHA(500),SS_ALPHA(500),
     *          SA_BETA(500),SS_BETA(500)

        DOUBLE PRECISION Q_VAR_A,DEL_Q_A,Q_VAR_B,DEL_Q_B,CORDS(500,3,3)
        
        CHARACTER BLAH*38

C     REQUIRED SUBROUTINES

        EXTERNAL POTIND,GENGRD,GASPOT,EV_SET_UP

C     --------------------- BEGIN -----------------------
C        WRITE(SO,*) 'GENTAB'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     GENERATE POTENTIAL FROM MEMORY PROTEINS
C     CONSTRUCT POTENTIAL AND FORCE TABLES
C     FIND TABLES FOR ALPHA-ALPHA COORDINATES

      DO 500 I_TAB=1,NUMTAB

C        DETERMINE WHICH RESIDUES WILL INTERACT GIVEN
C        THE PROTEIN LENGTH AND 'MUTATION' REGION

         IF( (CRDIXN(I_TAB,1).EQ.CRDIXN(I_TAB,2)).AND.
     *       (CRDIXN(I_TAB,1).EQ.1) )THEN
            JBEGN=2
            JEND=NMRES
            RSEP=5.0D0
         ELSEIF( CRDIXN(I_TAB,1).NE.CRDIXN(I_TAB,2) )THEN
            JBEGN=1
            JEND=NMRES
             RSEP=5.5D0
         ELSE
            JBEGN=1
            JEND=NMRES
            RSEP=10.5D0
         ENDIF

         IDIGNS=.FALSE.
         CALL POTIND(MAXSIZ,MAXCNT,NUMLNG(0,I_TAB),ILONG(1,1,I_TAB),JBEGN,JEND,NMRES,I_TAB)

C        SET ANALYSIS PARAMETERS TO FALSE;
C        FULL TABLES ARE TO BE GENERATED

C        CONSTRUCT R-GRID FOR POTENTIAL

C         IF( I_TAB.EQ.1 )IDIGNS=.FALSE.
         CALL GENGRD(MAXCNT,ILONG(1,1,I_TAB),
     *               NUMLNG(0,I_TAB),NMRES,MAXS,
     *               WORK6,WORK3,DELTA,DELTE,
     *               MAXSIZ,DELTZ(1,I_TAB),
     *               RINC(1,I_TAB),IDIGNS,
     *               OARCHV,RSEP,RCUTAMH)
         IDIGNS=.FALSE.

C        SET PASSI TRUE IF THIS IS THE FIRST CALL
C        TO GASPOT; OTHERWISE SET PASSI TO FALSE

         IF( I_TAB.EQ.1 )THEN
            PASSI=.TRUE.
         ELSE
            PASSI=.FALSE.
         ENDIF

C        SET PASSF TRUE IF THIS IS THE LAS TCALL
C        TO GASPOT; OTHERWISE SET PASSF TO FALSE

         IF( I_TAB.EQ.NUMTAB )THEN
            PASSF=.TRUE.
         ELSE
            PASSF=.FALSE.
         ENDIF

C         IF( I_TAB.EQ.2 )IDIGNS=.TRUE.
C         WRITE(SO,*) 'CALLING GASPOT'         
        CALL GASPOT(MAXSIZ,TARGET, BONDLN,I_TAB,NUMCRD,CRDIXN(I_TAB,1),CRDIXN(I_TAB,2),PASSI,PASSF)

         IDIGNS=.FALSE.

500     CONTINUE


!*****************************************************
! SET UP EXCLUDED VOLUMES

      CALL EV_SET_UP

! SET UP REPLICA INTERACTION IF REQUIRED

      IF (I_REP) CALL REP_CONTACT(TARGET)
        
        OPEN(ISS_STRUCT,FILE='PARAMS/ALPHZ',STATUS='OLD')
                                                                                                         
        READ(ISS_STRUCT,991)BLAH
991     FORMAT(A38)
                                                                                                         
        IF (BLAH.EQ."THIS IS A NON-CONTINUOUS PROTEIN CHAIN")THEN
           WRITE(6,*) 'I FOUND A NON-CONTINUTOUS PROTEIN'
           STOP
        ENDIF
        READ(ISS_STRUCT,199)TEMPNRES
199     FORMAT(I5)
        READ (ISS_STRUCT,25)(SEQ_A(I_RES),I_RES=1,TEMPNRES)
25      FORMAT(25(I2,1X))
        
        DO 100 ICOORD=1,3
            READ (ISS_STRUCT,30)(CORDS(I_RES,ICOORD,1),I_RES=1,TEMPNRES)
30          FORMAT(8(F8.3,1X))
100     CONTINUE
        DO 101 ICOORD=1,3
            READ (ISS_STRUCT,30)(CORDS(I_RES,ICOORD,2),I_RES=1,TEMPNRES)
101     CONTINUE
        DO 102 ICOORD=1,3
            READ (ISS_STRUCT,30)(CORDS(I_RES,ICOORD,3),I_RES=1,TEMPNRES)
102     CONTINUE
        READ(ISS_STRUCT,25)(SA_ALPHA(I_RES),I_RES=1,TEMPNRES)
        READ(ISS_STRUCT,25)(SS_ALPHA(I_RES),I_RES=1,TEMPNRES)
        CLOSE (ISS_STRUCT)

C    TARG_DIST USED IN Q_BIAS_SEG ALPHA

         DO 2501 I_RES1=1, TEMPNRES-1
          DO  2502 I_RES2=I_RES1 +1, TEMPNRES
            XDIFFTEMP=CORDS(I_RES1,1,1) - CORDS(I_RES2,1,1)
            YDIFFTEMP=CORDS(I_RES1,2,1) - CORDS(I_RES2,2,1)
            ZDIFFTEMP=CORDS(I_RES1,3,1) - CORDS(I_RES2,3,1)
            SS_DIST(I_RES1,I_RES2,1)=DSQRT(XDIFFTEMP**2 +
     *                        YDIFFTEMP**2 + ZDIFFTEMP**2)
2502       ENDDO
2501      ENDDO

        OPEN(ISS_STRUCT,FILE='PARAMS/BETAZ',STATUS='OLD')
        READ(ISS_STRUCT,991)BLAH
        IF (BLAH.EQ."THIS IS A NON-CONTINUOUS PROTEIN CHAIN")THEN
           WRITE(6,*) 'I FOUND A NON-CONTINUTOUS PROTEIN'
           STOP
        ENDIF
        READ(ISS_STRUCT,199)TEMPNRES
        READ (ISS_STRUCT,25)(SEQ_B(I_RES),I_RES=1,TEMPNRES)
        DO 200 ICOORD=1,3
            READ (ISS_STRUCT,30)(CORDS(I_RES,ICOORD,1),I_RES=1,TEMPNRES)
200     CONTINUE
        DO 201 ICOORD=1,3
            READ (ISS_STRUCT,30)(CORDS(I_RES,ICOORD,2),I_RES=1,TEMPNRES)
201     CONTINUE
        DO 202 ICOORD=1,3
            READ (ISS_STRUCT,30)(CORDS(I_RES,ICOORD,3),I_RES=1,TEMPNRES)
202     CONTINUE
        READ(ISS_STRUCT,25)(SA_BETA(I_RES),I_RES=1,TEMPNRES)
        READ(ISS_STRUCT,25)(SS_BETA(I_RES),I_RES=1,TEMPNRES)
        CLOSE (ISS_STRUCT)

C    TARG_DIST USED IN Q_BIAS_SEG BETA
                                                                                                         
         DO 3501 I_RES1=1, TEMPNRES-1
          DO  3502 I_RES2=I_RES1 +1, TEMPNRES
                                                                                                         
            XDIFFTEMP=CORDS(I_RES1,1,1) - CORDS(I_RES2,1,1)
            YDIFFTEMP=CORDS(I_RES1,2,1) - CORDS(I_RES2,2,1)
            ZDIFFTEMP=CORDS(I_RES1,3,1) - CORDS(I_RES2,3,1)
            SS_DIST(I_RES1,I_RES2,2)=DSQRT(XDIFFTEMP**2 + YDIFFTEMP**2 + ZDIFFTEMP**2)
3502       ENDDO
3501      ENDDO


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      MAKE TABLE FOR CALCULATION OF Q, AND DQ/DR<IJ>
C      BOTH SCALED BY 1/ ( 0.5 (N-1)(N-2) )
C      I HAVE NOTES ON THIS 'A SPECIFIC Q-DEPENDENT POTENTIAL'

        MAX_CONT=1E-8    !MAX CONT OF Q<IJ> TO Q BEFORE IGNORE IT

       IF (I_QBIAS_A)THEN
       DO  I_DIFF = 1,  NMRES
         Q_VAR_A=(FLOAT(I_DIFF)**(2.0D0*WIDTH_QEXP_A))
         MAX_RIJ=DSQRT( -2.0D0*LOG(MAX_CONT)*Q_VAR_A)
         DEL_R_A(I_DIFF)=MAX_RIJ/N_DIVS_MAX
           DO I_DIST=0,N_DIVS_MAX-1
       IF (SS_A) THEN
       Q_IJ_A(I_DIST,I_DIFF)=
     *     DEXP(-0.5D0*(FLOAT(I_DIST)*DEL_R_A(I_DIFF))**2/
     *     Q_VAR_A)*(2.0D0/(FLOAT(TOTALSSNORM)))
       ENDIF

      IF (.NOT. SS_A) THEN
         Q_IJ_A(I_DIST,I_DIFF)=
     *     DEXP(-0.5D0*(FLOAT(I_DIST)*DEL_R_A(I_DIFF))**2/
     *    Q_VAR_A)* (2.0D0/(FLOAT((NUMCONST_A-1)*(NUMCONST_A-2))))
      ENDIF
       
       ENDDO ! DO   I_DIST=0,N_DIVS_MAX-1
!       WRITE(6,*)'NUMCONST_A  IN GENTAB ',NUMCONST_A
         
         Q_IJ_A(N_DIVS_MAX,I_DIFF)=0.0D0
!     USE EXPANSION TO CALC GRADIENT, FOR ACCURACY
         DO I_DIST=1,N_DIVS_MAX
           SMALL=-(0.5D0/Q_VAR_A)*(DEL_R_A(I_DIFF)**2)
     *                     *(2.0D0*FLOAT(I_DIST-1)+1.0D0)
         DEL_Q_A=Q_IJ_A(I_DIST-1,I_DIFF)*(
     *          +SMALL+(1.0D0/2.0D0)*SMALL**2+(1.0D0/6.0D0)*SMALL**3+
     *          (1.0D0/24.0D0)*SMALL**4+(1.0D0/120.0D0)*SMALL**5)
         DQ_DR_IJ_A(I_DIST,I_DIFF)= (DEL_Q_A/DEL_R_A(I_DIFF))

C         IF (I_DIFF.EQ.20) THEN   ! UNCOMMENT TO CHECK Q_IJ_A AND DQ_IJ_A/DR_IJ
C             WRITE(6,*) FLOAT(I_DIST)*DEL_R_A(I_DIFF),
C    *         Q_IJ_A(I_DIST-1,I_DIFF),DQ_DR_IJ(I_DIST,I_DIFF),
C    *        (Q_IJ_A(I_DIST,I_DIFF)-Q_IJ_A(I_DIST-1,I_DIFF))/DEL_R_A(I_DIFF)
C          ENDIF
           ENDDO  ! DO I_DIST=1,N_DIVS_MAX
       ENDDO   ! I
       ENDIF  ! IF (I_QBIAS_A)THEN

         IF (I_QBIAS_B)THEN

C  NORMAILIZATION OF Q CONSTRAINT OF SECONDARY STRUCTURE UNITS
C  OVER RIGHT LENGTH OF RESIDUES FOR PROPER Q NORMALIZATION

       DO  I_DIFF = 1,  NMRES
         Q_VAR_B=(FLOAT(I_DIFF)**(2.0D0*WIDTH_QEXP_B))
         MAX_RIJ=DSQRT( -2.0D0*LOG(MAX_CONT)*Q_VAR_B)
         DEL_R_B(I_DIFF)=MAX_RIJ/N_DIVS_MAX
         DO I_DIST=0,N_DIVS_MAX-1
        
         IF (.NOT. SS_B) THEN
       Q_IJ_B(I_DIST,I_DIFF)=
     * DEXP(-0.5D0*(FLOAT(I_DIST)*DEL_R_B(I_DIFF))**2/
     * Q_VAR_B ) *(2.0D0/(FLOAT((NUMCONST_B-1)*(NUMCONST_B-2))))
        ENDIF

      IF (SS_B) THEN
         Q_IJ_B(I_DIST,I_DIFF)=DEXP(-0.5D0*(FLOAT(I_DIST)*DEL_R_B(I_DIFF))**2/
     *Q_VAR_B)*(2.0D0/(FLOAT((TOTALSSNORM-1)*(TOTALSSNORM-2))))
        ENDIF
                                                                                                         
           ENDDO ! DO I_DIST=0,N_DIVS_MAX-1
                                                                                                         
         Q_IJ_B(N_DIVS_MAX,I_DIFF)=0.0D0
                                                                                                         
!     USE EXPANSION TO CALC GRADIENT, FOR ACCURACY
                                                                                                         
         DO I_DIST=1,N_DIVS_MAX
                                                                                                         
           SMALL=-(0.5D0/Q_VAR_B)*(DEL_R_B(I_DIFF)**2)*(2.0D0*FLOAT(I_DIST-1)+1.0D0)
                                                                                                         
           DEL_Q_B=Q_IJ_B(I_DIST-1,I_DIFF)*(
     *            +SMALL+(1.0D0/2.0D0)*SMALL**2+(1.0D0/6.0D0)*SMALL**3+
     *            (1.0D0/24.0D0)*SMALL**4+(1.0D0/120.0D0)*SMALL**5)


           DQ_DR_IJ_B(I_DIST,I_DIFF)= (DEL_Q_B/DEL_R_B(I_DIFF))
                                                                                                         
C         IF (I_DIFF.EQ.20) THEN    ! UNCOMMENT TO CHECK Q_IJ AND DQ_IJ/DR_IJ
C         WRITE(6,*) FLOAT(I_DIST)*DEL_R_B(I_DIFF),
C     *     Q_IJ_B(I_DIST-1,I_DIFF),DQ_DR_IJ_B(I_DIST,I_DIFF),
C     *   (Q_IJ_B(I_DIST,I_DIFF)-Q_IJ_B(I_DIST-1,I_DIFF))/DEL_R_B(I_DIFF)
C          ENDIF
                                                                                                         
               ENDDO  ! DO I_DIST=1,N_DIVS_MAX
           ENDDO   ! I
           ENDIF ! IF (I_QBIAS_B)THEN

C       ALSO CALCULATE SEPARATIONS IN THE TARGET, WHICH WILL
C       NEED LATER IN SUBROUTINE Q_BIAS
C       (TRANSFERRED THERE VIA GLOBALS)
C       --OR ALTERNATIVELY READ IN STRUCTURE TO BIAS TO

           IF ((.NOT.I_BIAS_NATIVE_A .AND.  I_QBIAS_A ).OR.
     *         (.NOT.I_BIAS_NATIVE_B .AND.  I_QBIAS_B )) THEN

C         IF (.NOT.I_BIAS_NATIVE_A .AND.  I_QBIAS_A )THEN

           OPEN(IBIASFILE,FILE='MOVIESEG_BIAS',STATUS='UNKNOWN',
     *                                            ACTION='READ')
           READ(IBIASFILE,1040) NMRES_CHECK,IDUMMY,IDUMMY,IDUMMY
1040       FORMAT(4(I8,1X))
           IF (NMRES_CHECK.NE.NMRES) THEN
             WRITE(SO,*) 'MOVIESEG_BIAS FILE WRONG',NMRES_CHECK,NMRES
             STOP
           ENDIF
           READ(IBIASFILE,*)
1050       FORMAT(3X,3(1X,F8.4),2(4X,3(1X,F8.4)))
           DO I_RES=1,NMRES
             READ(IBIASFILE,1050)  BIAS_CORD(I_RES,1,1),
     *       BIAS_CORD(I_RES,2,1),BIAS_CORD(I_RES,3,1),   
     *       BIAS_CORD(I_RES,1,2),BIAS_CORD(I_RES,2,2),   
     *       BIAS_CORD(I_RES,3,2),BIAS_CORD(I_RES,1,3),   
     *       BIAS_CORD(I_RES,2,3),BIAS_CORD(I_RES,3,3)
           ENDDO
           CLOSE(IBIASFILE)

          ENDIF 

        IF ((I_BIAS_NATIVE_A .AND.  I_QBIAS_A ).OR.
     *      (I_BIAS_NATIVE_B .AND.  I_QBIAS_B )) THEN

             BIAS_CORD=TARGET
             WRITE(6,*) 'TARGET  USED FOR CONSTRAINT'

         DO 1001 TAB=1,MAXTAB
         DO 501 I_IXN=1,NUMLNG(NMRES,TAB)

            ISIT1=ILONG(I_IXN,1,TAB)
            ISIT2=ILONG(I_IXN,2,TAB)

            XDIFF=BIAS_CORD(ISIT1,1,CRDIXN(TAB,1)) -
     *                  BIAS_CORD(ISIT2,1,CRDIXN(TAB,2))
            YDIFF=BIAS_CORD(ISIT1,2,CRDIXN(TAB,1)) -
     *                  BIAS_CORD(ISIT2,2,CRDIXN(TAB,2))
            ZDIFF=BIAS_CORD(ISIT1,3,CRDIXN(TAB,1)) -
     *                  BIAS_CORD(ISIT2,3,CRDIXN(TAB,2))

            TARG_DIST(I_IXN,TAB)=DSQRT(XDIFF**2 + YDIFF**2 + ZDIFF**2)

501      ENDDO
1001     ENDDO

       ENDIF
C     ---------------------- DONE -----------------------
C        WRITE(SO,*) 'LEAVING GENTAB'      RETURN
      END
