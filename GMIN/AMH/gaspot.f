C   --------------------- GASPOT ----------------------

      SUBROUTINE GASPOT(MAXSD,TARGET,BONDLN,CURTAB,NMCRD,CDTYP1,CDTYP2,PASSI,PASSF)

C     ---------------------------------------------------

C     GASPOT FIND UNGENERALIZED POTENTIAL

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

C     ARGUMENTS:

C        MAXSD - MAXIMUM PROTEIN LENGTH (I)
C        TARGET- TARGET PROTEIN (I)
C        BONDLN- TARGET PROTEIN'S BOND LENGTHS (O)
C        CURTAB- INDEX FOR CURRENT TABLE TO BE
C                CONSTRUCTED (I)
C        NMCRD- NUMBER OF COORDINATE TYPES (I)
C        CDTYP1- TPYE ID FOR COORDINATE SET 1 (I)
C        CDTYP2- TPYE ID FOR COORDINATE SET 2 (I)
C        PASSI - TRUE ON FIRST PASS THROUGH 
C                SUBROUTINE; OTHERWISE FALSE (I)[
C        PASSF - TRUE ON LAST PASS THROUGH 
C                SUBROUTINE; OTHERWISE FALSE(I)
C
C      VARIABLE NAMES LIFTED FROM GAUSSV.F
C        GAUSSP- GAUSSIAN AS A FUNCTION OF THE R-GRID
C                IN RGRID (O)
C        RGRID - GRID OF R POINTS FOR WHICH THE GAUSSIAN
C                IS TO BE COMPUTED (I)        
C      -----------------------------------------------
C
C
C
C        DIST_IJ        DISTANCE I TO J IN THE MEMORY
C        ID1
C     ---------------------------------------------------
C
C  NB:
C      PART OF THE REASON FOR THE LENGTH, AND CONFUSION IN THIS SUBROUTINE IS
C     HISTORICAL IN NATURE.  BEFORE HOMOLOGY MATCHING OF RESIDUES, EACH 
C     RESIDUE WAS COMPARED TO OTHER RESIDUES IN A SLIDING WINDOW,
C     FROM -N11 TO N11, AND FROM -N33 TO N33, WERE COMPARED.
C
C
C PROGRAM LOGIC
C -------------
C       RTLIST = 0
C
C
C509        LOOP OVER MEMORIES (25%)
C                CALL GETMEM (29%)
C               CALL HPROFL (HYDSCL
C                CALL HPROFL (HYDSCL
C
C                IF FIRST MEMORY (TARGET)  -->
C                       CALL PTTARG
C                       CONTINUE TO NEXT MEMORY (35%)
C
C                IF FIRST MEMORY PROTIEN;
C                        COMPARE NAMES PROTNM(1) AND PROTNM(0)
C
C               READ IN MATCH FILE
C
C519                LOOP OVER EACH INTERACTION (43%)
C
C
C                        DETERMINE CATEGORY OF INTERACTION (55%)
C                        CALL CHARGE
C                       IF CHARGE = 0, TO TO 519 (? NEXT IXN)
C                        CALL GAUSSV TO CREATE GAUSSIAN CURVE AROUND DIST_IJ
C
C                        CONTINUE TO NEXT INTERACTION (80%)
C
C      CONTINUE TO NEXT MEMORY (80%)
C
      USE AMHGLOBALS,  ONLY:SO, MAXS,MAXCNT,MAXMEM,NUMLNG,NMRES,
     *      IPROLST,VPOTNT,ICON,FORSE,DIST_CUT,
     *      IDIGNS,NUMMEM,PROTNM,IMEMRI,MAXRES,JRES,NUMCRD,
     *      MEMPRE,OARCHV,SA,IWORK,AMINOA,HYDSEQ,
     *      MAXSIZ,IRES,SHYDRO,TARPRE,NUMTAB,RAN_FORCE,
     *      RAN_FILE,IRAN,ISEED_AMH,S_RAN,LAMBDAR,ALLOW_NEG,ORAN,ILONG,
     *      I_ALT_PROX,RINC,SRCUT,ALT_PROX_CUT,RAN_MIN_SEQ_DIST,
     *      QCHRG,QCHRG2,DELTZ,WELSCL,R_RAN,
     *      RINCINV,RINCSQ,WORK2,FOUR_PLUS,
     *      YWORK,WORK8,HYDSCL,MIN_SEQ_SEP,I_3RD_IS_CONTACT,
     *      N_LETTERS,R_MIN,R_MAX,GLY_CON,I_V_TEST,TEST_SITE,
     *      IXN_FROM_SITE,MINMR,MAXMR,MAX_WELL,ALPHA_C_OF_N,
     *      AB_C_OF_N_OLD,AB_C_OF_N_NEW,NUM_WELL,I_CONTACT_ORDER,
     *      I_CONTACT_ORDER_MIN,I_CONTACT_ORDER_MAX,R_MIN_CONTACT_ORDER,
     *      R_MAX_CONTACT_ORDER,GAMMA_CONTACT_ORDER,IMAT,
     *      N_CONTACT_ORDER_TERMS,NUMSEQ,AVE_SEQ,ITARG_SEQ,
     *      AVE_SEQ_AMC,NUMSEQ_AMC,TGSEQUENCES_AMC,
     *      OVGASPOT,MAXSEQ,TGSEQUENCES,TARG_CONS,GO_CON,GO_CON_DIST

      IMPLICIT NONE


C     ARGUMENT DECLARATIONS:

         LOGICAL PASSI,PASSF,IPROX_RAN

         INTEGER MAXSD,CURTAB,CDTYP1,CDTYP2,NMCRD,IPROX,HELIX,I_CONTACT_TERM
     
         DOUBLE PRECISION TARGET(MAXSD,3,NMCRD),BONDLN(MAXSD,NMCRD)

C     INTERNAL VARIABLES:

         LOGICAL PASST,HELIX_A,HELIX_B
         DOUBLE PRECISION RGRID(MAXS), GAUSSP(MAXS),GAUSSQ(MAXS),RAN_NUM(1:MAXCNT),DELT_SAFE,THETA,RNMRES,
     *        THETA_DOT,FORCE_TERM(6),VPOTINT_TEST(1:MAXS+1),T_MIN,T_MAX,R_TEMP

         INTEGER INDXT,NMRSS(0:MAXMEM),MAXSZ,MAXS2, MATCH(MAXRES),J,I_R,II,
     *           I,ID1, ID2, ID3, ID4,IDUMMY,I_TEST_TAB,
     *           I_TEST_IXN,TAB_FOR_CONTACT,OPEN_STATUS,COUNT,I_MEM, I_IXN, I_RES,I_WELL

        INTEGER I514,I515,I544,I1,I2,I516,I507,I522,I588,I504,I523,I589,I503,I532,I518

         INTEGER ISEQ,CONNMRES,TEMP_MCP_PASS

         DOUBLE PRECISION RNORM,DELT2,TOTCHG(MAX_WELL),LONG_NFACTOR(MAX_WELL),DIST_IJ,M,C,DIST_PROX

         CHARACTER*5 PROFL,TARFL
         CHARACTER*42 CONFILE
         CHARACTER*33 PROFILE
         CHARACTER*36 MATFILE
         CHARACTER*10 CCOUNT1,CCOUNT2,CCOUNT3
         CHARACTER*30 BLAH
C     REQUIRED SUBROUTINES

         EXTERNAL GETMEM,GAUSSV,HPROFL,PTTARG,GAUSSW,NUM_TO_CHAR,SLARNV

C     --------------------- BEGIN -----------------------
C     SET VARIOUS VARIABLES

      INDXT=NUMLNG(NMRES,CURTAB)

C     REWIND DATA BASE FILE
      REWIND IPROLST

C     SET NUMBER POINTS TO BE ANALYZED OR SET THE TABLE SIZE

       MAXSZ=MAXS
       MAXS2=MAXS

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   SET LONG_NFACTOR WHICH IS SOME FUDGE FACTOR FOR THE
C   CONTACT POTENTIALS (IT PARTIALLY ALLOWS FOR THE EFFECT OF 
C   SIZE -- WHICH THE FRACTION OF CONTACTS AT A CERTAIN DISTANCE
C   DEPENDS ON). IT SHOULD BE THE SAME AS INCLUDED IN THE
C   OPTIMISATION. IN COREY'S CODE IT IS IN QCHRGMK.F INSTEAD.
C   JOHAN: ADDED LONG_NFACTOR FOR 2-WELLS 

        LONG_NFACTOR=1.0D0

C   START IF ALPHA_C_OF_N 
        IF ( ALPHA_C_OF_N ) THEN

       DO II=1,MAX_WELL,1
         LONG_NFACTOR(II)=1.0D0
       ENDDO
       RNMRES=REAL(NMRES)

        IF (NUM_WELL .EQ. 2) THEN
         LONG_NFACTOR(2)=1.0D0/(RNMRES*0.012D0 + 0.87D0)
        ELSEIF (NUM_WELL .EQ. 3) THEN
         LONG_NFACTOR(2)=1.0D0/(RNMRES*0.0065D0 + 0.87D0)
         LONG_NFACTOR(3)=1.0D0/(RNMRES*0.04187261D0 + 0.1256658D0)
        ELSEIF (NUM_WELL .EQ. 10) THEN
         LONG_NFACTOR(1)=1.0D0/(RNMRES*0.0008D0 + 0.09D0)
         LONG_NFACTOR(2)=1.0D0/(RNMRES*0.0009D0 + 0.16D0)
         LONG_NFACTOR(3)=1.0D0/(RNMRES*0.001D0 + 0.20D0)
         LONG_NFACTOR(4)=1.0D0/(RNMRES*0.003D0 + 0.29D0)
         LONG_NFACTOR(5)=1.0D0/(RNMRES*0.004D0 + 0.53D0)
         LONG_NFACTOR(6)=1.0D0/(RNMRES*0.004D0 + 0.76D0)
         LONG_NFACTOR(7)=1.0D0/(RNMRES*0.005D0 + 0.77D0)
         LONG_NFACTOR(8)=1.0D0/(RNMRES*0.005D0 + 0.94D0)
         LONG_NFACTOR(9)=1.0D0/(RNMRES*0.006D0 + 1.18D0)
         LONG_NFACTOR(10)=1.0D0/(RNMRES*0.013D0 + 1.8D0)
        ELSE
           WRITE(SO,*) 'UNSUPPORTED NO OF WELLS (ALPHA_C_OF_N)',NUM_WELL
           STOP
        ENDIF

       ENDIF
C   END IF ALPHA_C_OF_N 

       IF (AB_C_OF_N_NEW) THEN
       RNMRES=REAL(NMRES)
        IF (NUM_WELL .EQ. 2) THEN
          LONG_NFACTOR(1)= ( 0.035D0*RNMRES)/(RNMRES*0.043D0 + 1.0D0)
          LONG_NFACTOR(2)=( 0.07D0*RNMRES)/(RNMRES*0.023D0 + 1.0D0)
           LONG_NFACTOR(1)=1.0D0/(LONG_NFACTOR(1))
           LONG_NFACTOR(2)=1.0D0/(LONG_NFACTOR(2))
        ELSEIF (NUM_WELL .EQ. 3) THEN
          LONG_NFACTOR(1)= ( 0.0843467D0*RNMRES)/(RNMRES*0.0453928D0 + 1.0D0)
          LONG_NFACTOR(2)=( 0.0669808D0*RNMRES)/(RNMRES*0.025112D0 + 1.0D0)
          LONG_NFACTOR(3)=( 0.18665D0*RNMRES) /(RNMRES*0.0107983D0 + 1.0D0)
           LONG_NFACTOR(1)=1.0D0/(LONG_NFACTOR(1))
           LONG_NFACTOR(2)=1.0D0/(LONG_NFACTOR(2))
           LONG_NFACTOR(3)=1.0D0/(LONG_NFACTOR(3))
        ELSEIF (NUM_WELL .EQ. 5) THEN

         LONG_NFACTOR(1)=( 0.0297375D0*RNMRES) /(RNMRES*0.02977935D0 + 1.0D0)
         LONG_NFACTOR(2)=( 0.0389704D0*RNMRES) /(RNMRES*0.021101D0 + 1.0D0)
         LONG_NFACTOR(3)=( 0.0596751D0*RNMRES) /(RNMRES*0.0133269D0  + 1.0D0)
         LONG_NFACTOR(4)=( 0.0681322D0*RNMRES) /(RNMRES*0.0100256D0 + 1.0D0)
         LONG_NFACTOR(5)=( 0.0729201D0*RNMRES) /(RNMRES*0.00347563D0 + 1.0D0)

         LONG_NFACTOR(1)=1.0D0/(LONG_NFACTOR(1))
         LONG_NFACTOR(2)=1.0D0/(LONG_NFACTOR(2))
         LONG_NFACTOR(3)=1.0D0/(LONG_NFACTOR(3))
         LONG_NFACTOR(4)=1.0D0/(LONG_NFACTOR(4))
         LONG_NFACTOR(5)=1.0D0/(LONG_NFACTOR(5))

       ELSEIF (NUM_WELL .EQ. 10) THEN

        LONG_NFACTOR(1)=( 0.0785047D0*RNMRES) /(RNMRES*0.245032D0 + 1.0D0)
        LONG_NFACTOR(2)=( 0.0152761D0*RNMRES) /(RNMRES*0.0283803D0 + 1.0D0)
        LONG_NFACTOR(3)=( 0.205481D0*RNMRES) /(RNMRES*0.385813D0  + 1.0D0)
        LONG_NFACTOR(4)=( 0.0174765D0*RNMRES) /(RNMRES*0.0174638D0 + 1.0D0)
        LONG_NFACTOR(5)=( 0.0352685D0*RNMRES) /(RNMRES*0.0269838D0 + 1.0D0)
        LONG_NFACTOR(6)=( 0.0474026D0*RNMRES) /(RNMRES*0.0249249D0 + 1.0D0)
        LONG_NFACTOR(7)=( 0.18665D0 *RNMRES) /(RNMRES*0.0107983D0 + 1.0D0)
        LONG_NFACTOR(8)=( 0.0390303D0*RNMRES) /(RNMRES*0.0140943D0 + 1.0D0)
        LONG_NFACTOR(9)=( 0.0327411D0*RNMRES) /(RNMRES*0.00812347D0 + 1.0D0)
        LONG_NFACTOR(10)=( 0.0561461D0*RNMRES) /(RNMRES*0.00743991D0 + 1.0D0)

           LONG_NFACTOR(1)=1.0D0/(LONG_NFACTOR(1))
           LONG_NFACTOR(2)=1.0D0/(LONG_NFACTOR(2))
           LONG_NFACTOR(3)=1.0D0/(LONG_NFACTOR(3))
           LONG_NFACTOR(4)=1.0D0/(LONG_NFACTOR(4))
           LONG_NFACTOR(5)=1.0D0/(LONG_NFACTOR(5))
           LONG_NFACTOR(6)=1.0D0/(LONG_NFACTOR(6))
           LONG_NFACTOR(7)=1.0D0/(LONG_NFACTOR(7))
           LONG_NFACTOR(8)=1.0D0/(LONG_NFACTOR(8))
           LONG_NFACTOR(9)=1.0D0/(LONG_NFACTOR(9))
           LONG_NFACTOR(10)=1.0D0/(LONG_NFACTOR(10))

        ELSE
            WRITE(SO,*) 'NUM_WELL PROBLEM FOR AB_C_OF_N_NEW',NUM_WELL
            STOP 
        ENDIF
       ENDIF

       IF (AB_C_OF_N_OLD) THEN
       RNMRES=REAL(NMRES)
       LONG_NFACTOR(1)=1.0D0/(RNMRES*0.0015D0 + 1.94D0)
       LONG_NFACTOR(2)=1.0D0/(RNMRES*0.0032D0 + 1.83D0)
       LONG_NFACTOR(3)=1.0D0/(RNMRES*0.022D0 + 7.77D0)

       IF (NUM_WELL.NE.3) THEN 
         WRITE(SO,*) 'NUMWELL MUST BE 3 FOR AB_C_OF_N_*OLD*',NUM_WELL
         STOP
       ENDIF

       ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INITIALIZE POTENTIAL AND FORCE TABLE

        DO 514 I514=1,INDXT
          DO 515 I515=0,MAXS2+1
            VPOTNT(I515,I514,CURTAB)=0.0D0
            FORSE(I515,I514,CURTAB)=0.0D0
  515     CONTINUE
  514   CONTINUE

      IDIGNS=.FALSE.
      IF( PASSI )THEN
         PASST=.TRUE.
      ELSE
         PASST=.FALSE.
      ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   READ TARGET (WHICH IS CONSIDERED 0TH MEMORY)
         READ (IPROLST,1000)TARFL
1000     FORMAT(A5)
C         PROFILE='/HOME/MPRENTIS/AMH/PROTEINS/'//TARFL
CCC           PROFILE='/HOME/MPRENTIS/AMH/PROTEINS/'//TARFL

         PROTNM(0)=TARFL
         MATFILE='MATCH/'//TRIM(TARFL)//'/'//TARFL

          OPEN(IMEMRI,FILE='PROTEINS/'//TARFL,
     *                          STATUS='OLD',IOSTAT=OPEN_STATUS)

         IF (OPEN_STATUS.NE.0) THEN
           WRITE(6,*) 'FAILURE TO OPEN PROFILE 0TH MEM FILE '
           WRITE(6,*) 'ERROR NUMBER ',OPEN_STATUS
           STOP
         ENDIF

C        READ IN MEMORY PROTEIN COORDINATES, AND PRIMARY SEQUENCE AND CRYSTAL
C        SECONDARY STRUCTURE

!         WRITE(6,*)'CALL GETMEM '

         CALL GETMEM(PROTNM(0),NMRSS(0),MAXRES,JRES,IMEMRI,NUMCRD,YWORK,
     *               MEMPRE(1,CURTAB),OARCHV,SA,IWORK,PASST,0)

         CLOSE(IMEMRI)

         WORK8(1)=FLOAT(NMRSS(0))

C        SET HYDROPHOBICITY PROFILE

         CALL HPROFL(MAXRES,NMRSS(0),JRES,HYDSEQ(1,1,CURTAB),HYDSCL(0,1))
         CALL HPROFL(MAXRES,NMRSS(0),JRES,HYDSEQ(1,2,CURTAB),HYDSCL(0,2))

C   PLACE COORDINATES IN TARGET
C       AND SEQUENCE PROFILE IN TRES

               CALL PTTARG(MAXSIZ,NMRSS(0),NUMCRD,TARGET,MAXRES,YWORK,IRES,JRES,
     *                     SHYDRO(1,1,CURTAB),HYDSEQ(1,1,CURTAB),TARPRE(1,CURTAB),
     *                     MEMPRE(1,CURTAB),BONDLN,OARCHV,PASSI)


C  END OF READING TARGET INFO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     CONTACT ORDER TERM

       IF (I_CONTACT_ORDER.AND.CURTAB.EQ.4) THEN
         DO I_IXN=1,INDXT
           ID3=ILONG(I_IXN,1,CURTAB)
           ID4=ILONG(I_IXN,2,CURTAB)
           DO I_CONTACT_TERM=1,N_CONTACT_ORDER_TERMS
           IF (  (ID4-ID3.GE.I_CONTACT_ORDER_MIN(I_CONTACT_TERM)) .AND.
     *        (ID4-ID3.LE.I_CONTACT_ORDER_MAX(I_CONTACT_TERM)) .AND.
     *                    IRES(ID3).NE.8 .AND.
     *                        IRES(ID4).NE.8 )  THEN
              DO I_R=1,MAXSZ
                 R_TEMP=RINC(I_IXN,CURTAB)*REAL(I_R)
          T_MIN=TANH(7.0D0*(R_TEMP-R_MIN_CONTACT_ORDER(I_CONTACT_TERM)))
          T_MAX=TANH(7.0D0*(R_MAX_CONTACT_ORDER(I_CONTACT_TERM)-R_TEMP))
                 THETA = (GAMMA_CONTACT_ORDER(I_CONTACT_TERM,1)
     *              +GAMMA_CONTACT_ORDER(I_CONTACT_TERM,2)*
     *           ((ID4-ID3)**GAMMA_CONTACT_ORDER(I_CONTACT_TERM,3)) )*
     *              0.25D0*( 1.0D0+T_MIN )*( 1.0+T_MAX )
                 THETA_DOT=7.0D0*THETA*(T_MAX-T_MIN)

                 VPOTNT(I_R,I_IXN,CURTAB)=VPOTNT(I_R,I_IXN,CURTAB)+THETA
                 FORSE(I_R,I_IXN,CURTAB)=FORSE(I_R,I_IXN,CURTAB)-THETA_DOT/R_TEMP
              ENDDO

           ENDIF
           ENDDO
         ENDDO
       ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     READ IN SEQUENCES TO AVERAGE FORCE OVER, IF AVG_SEQ IS ON
C      IF (AVE_SEQ) THEN
C          WRITE(6,*)'TARGET_SEQUENCES'
C        OPEN(ITARG_SEQ,FILE='TARGET_SEQUENCES',STATUS='OLD')
C           READ(ITARG_SEQ,*)NUMSEQ
C           WRITE(6,*)'NUMSEQ'
C        DO ISEQ = 1,NUMSEQ
C          READ(ITARG_SEQ,999)(TGSEQUENCES(ID3,ISEQ),ID3=1,NMRES)
C          WRITE(6,999)(TGSEQUENCES(ID3,ISEQ),ID3=1,NMRES)
C999       FORMAT(25(I2,1X))
C        ENDDO
C        CLOSE(ITARG_SEQ)
C      ELSE
C        NUMSEQ = 1
C!        TGSEQUENCES(1:NMRES,1)=IRES(1:NMRES)
C      ENDIF
C          IF((AVE_SEQ) .AND. (NUMSEQ.GT.MAXSEQ) ) THEN
C            WRITE(6,*) 'NUMSEQ GREATER THAN MAXSEQ'
C            STOP
C          ENDIF



         IF (TARG_CONS) THEN
       CONFILE='/HOME/MPRENTIS/AMH/MD_INPUT/TARGCONS/'//TARFL
C               123456789012345678901234567890123456789012345678
          WRITE(6,*)'TARG_CONS DIRECTORY = ' , CONFILE
          OPEN(ICON,FILE=CONFILE,STATUS='OLD',IOSTAT=OPEN_STATUS)
               IF (OPEN_STATUS.NE.0) THEN
                 WRITE(6,*) 'FAILURE TO OPEN FILE ',PROFILE
                 WRITE(6,*) 'ERROR NUMBER ',OPEN_STATUS
                 STOP
               ENDIF
               READ(ICON,*)BLAH
               READ(ICON,102)CONNMRES
               READ(ICON,999)(TGSEQUENCES(I1,NUMSEQ),I1=1,NMRES)
999            FORMAT(25(I2,1X))
              WRITE(6,*)'TARG_CONS SEQUENCE'
               WRITE(6,999)(TGSEQUENCES(I1,NUMSEQ),I1=1,NMRES)
                                                                                
102            FORMAT (I5)
                WRITE(6,*)'CONNMRES =',CONNMRES
                IF( CONNMRES .GT.MAXRES )THEN
C                   WRITE(OARCHV,611)PROTNM,NMRSS,MAXRES
611                FORMAT('TARG_CONS ',A5,' TOO LARGE ',I4,
     *                   ' RESERVED SPACE ',I4)
                     WRITE(6,611)PROTNM,NMRSS,MAXRES
                     STOP
                ENDIF
                IF (CONNMRES .NE. NMRES )THEN
C                   WRITE(OARCHV,612)CONNMRES,NMRES
612                FORMAT('CONNMRES',I5,' .NE. NMRES ',I5,' 612 GASPOT')
                   WRITE(6,612)CONNMRES,NMRES
                   STOP
                ENDIF
          CLOSE(ICON)
        ENDIF ! IF (TARG_CONS)

         IF (CURTAB.EQ.1) THEN
C                   WRITE(6,*)'TARGET  SEQUENCE   ',TARG_CONS
C               DO 909 A2=1,NUMSEQ
C                    WRITE(6,999)(TGSEQUENCES(I1,A2),I1=1,NMRES)
C909            CONTINUE
         ENDIF  
                                                                          
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     LOOP OVER EACH MEMORY

      DO 509 I_MEM=1,NUMMEM
         READ (IPROLST,1000)PROFL
                                                                                
         PROTNM(I_MEM)=PROFL
CCC         MATFILE='/HOME/MPRENTIS/AMH/MATCH/'//TRIM(TARFL)//'/'//PROFL
C INSERT TARG_CONS SEQ PART HERE FOR MEMORIES
C  /HOME/MPRENTIS/AMH/CONSENS_SEQ

C           WRITE(6,*) 'MATFILE ',MATFILE
C           WRITE(6,*) 'PROFILE ',PROFILE 

          OPEN(IMEMRI,FILE='PROTEINS/'//PROFL,STATUS='OLD',IOSTAT=OPEN_STATUS)

         IF (OPEN_STATUS.NE.0) THEN
           WRITE(SO,*) ' PROTEINS ',PROFL 
           WRITE(SO,*) 'FAILURE TO OPEN PROTEIN 543 FILE '
           WRITE(SO,*) 'ERROR NUMBER ',OPEN_STATUS
           STOP
         ENDIF
 
C        READ IN MEMORY PROTEIN COORDINATES, AND PRIMARY SEQUENCE AND CRYSTAL 
C        SECONDARY STRUCTURE

         CALL GETMEM(PROTNM(I_MEM),NMRSS(I_MEM),MAXRES,JRES,IMEMRI,NUMCRD,YWORK,
     *               MEMPRE(1,CURTAB),OARCHV,SA,IWORK,PASST,I_MEM)

         CLOSE(IMEMRI)
 
         WORK8(I_MEM+1)=FLOAT(NMRSS(I_MEM))

C        SET HYDROPHOBICITY PROFILE

         CALL HPROFL(MAXRES,NMRSS(I_MEM),JRES,HYDSEQ(1,1,CURTAB),HYDSCL(0,1))
         CALL HPROFL(MAXRES,NMRSS(I_MEM),JRES,HYDSEQ(1,2,CURTAB),HYDSCL(0,2))

C        READ IN AND LABEL SECONDARY STRUCTURE BASED ON DATA BANK LABLES

C000000000000000000000000000000000000000000000000000000000000000000


         IF( PASSF.AND.(I_MEM.EQ.1) )THEN
C            WRITE(OARCHV,905)PROTNM(0),PROTNM(I_MEM)
  905       FORMAT(/'2ND STRUCTURES FOR ',A5,2X,A5)
            DO 544 I544=1,NMRSS(I_MEM)
               WRITE(OARCHV,904)I544,AMINOA(JRES(I544)),
     *                          (TARPRE(I544,I1),
     *                          (INT(SHYDRO(I544,I2,I1)),
     *                           I2=1,2),I1=1,NUMTAB)
  904          FORMAT(I4,1X,A3,2X,4(3(I2,1X),1X))
  544       CONTINUE
         ENDIF

CMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C  READ IN MATCH FILE
C
C          READ IN ALIGNMENTS
           DO 670 I1=1,MAXRES
             MATCH(I1)=0
670           CONTINUE
         OPEN(IMAT,FILE='MATCH/'//TRIM(TARFL)//'/'//PROFL,
     *                     STATUS='OLD',IOSTAT=OPEN_STATUS)

         IF (OPEN_STATUS.NE.0) THEN
           WRITE(SO,*) 'MATCH ', PROFL 
           WRITE(SO,*) 'FAILURE TO OPEN FILE MATCH 594 '
           WRITE(SO,*) 'ERROR NUMBER ',OPEN_STATUS
           STOP
         ENDIF

           READ (IMAT,*)
           READ (IMAT,1002)(MATCH(I_RES), I_RES=1,NMRES)
1002       FORMAT(10(I4))
          CLOSE (IMAT)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         IF (RAN_FORCE .AND. CURTAB.EQ.1) THEN
           IF (RAN_FILE) THEN
             OPEN(IRAN,FILE='RANDOM_IN',STATUS='OLD')
             READ(IRAN,*) IDUMMY
              IF (IDUMMY .NE. INDXT) THEN
                WRITE(SO,*) 'WRONG NO OF INTERACTIONS IN RANDOM_IN'
                WRITE(SO,*) IDUMMY,INDXT
                STOP
              ENDIF
              READ(IRAN,*) (RAN_NUM(I), I=1,INDXT)
           ELSE
             TEMP_MCP_PASS=3
             CALL SLARNV(TEMP_MCP_PASS,ISEED_AMH(1),INDXT,RAN_NUM) !SHOULD BE GAUSSIAN DIST MEAN 0, SD 1  
             DO J=1,INDXT
               RAN_NUM(J) = RAN_NUM(J)*S_RAN + LAMBDAR
               IF (RAN_NUM(J).LT.0.0D0.AND.(.NOT.ALLOW_NEG)) RAN_NUM(J)=0.0D0
             ENDDO
           ENDIF
           OPEN(ORAN,FILE='RANDOM_OUT',STATUS='NEW')
             WRITE(ORAN,*) INDXT
             WRITE(ORAN,*) (RAN_NUM(I), I=1,INDXT)
           CLOSE(ORAN)
         ENDIF


C        LOOP OVER EACH CONSTRAINT

         DO 519 I_IXN=1,INDXT

C            SET R-GRID
               
             DO 516 I516=1,MAXSZ
               RGRID(I516)=RINC(I_IXN,CURTAB)*FLOAT(I516)
516          CONTINUE

C      ID3, ID4 ARE THE I AND J INDICES (RESPECTIVELY) FOR THE TARGET PROTEIN

            ID3=ILONG(I_IXN,1,CURTAB)
            ID4=ILONG(I_IXN,2,CURTAB)

C      ID1, ID2 ARE THE I AND J INDICES (RESPECTIVELY) IN THE MEMORY PROTEIN
C      THAT ID3 AND ID4 ARE ALIGNED TO

            ID1=MATCH(ID3)
            ID2=MATCH(ID4)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C USE MEMORY SECONDARY STRUCTURAL INFO: CURRENTLY IN SOME ENCODINGS THE INTERACTIONS
C ARE GIVEN A DIFFERENT WEIGHT (IE GAMMA VALUE) IF *EITHER* ONE OF THE RESIDUES
C INVOLVED IS ALIGNED TO A HELICAL RESIDUE IN THE MEMORY
C I MADE THIS BIT OF THE CODE A BIT MORE LONG-WINDED THAN BEFORE TO AVOID LOOKING UP
C MEMPRE(ID1,CURTAB) WHEN ID1=0 (NON-ALIGNED) BECAUSE THIS IS OUTSIDE THE ARRAY.

              IF (ID1.EQ.0) THEN  !FIRST RESIDUE NOT ALIGNED => NOT ALIGNED TO HELIX
                HELIX_A=.FALSE.
              ELSE
                IF (MEMPRE(ID1,CURTAB) .EQ. 1) THEN 
                   HELIX_A=.TRUE.    !FIRST RESIDUE ALIGNED TO HELIX
                ELSE
                   HELIX_A=.FALSE.    !FIRST RESIDUE *NOT* ALIGNED TO HELIX
                ENDIF
              ENDIF
              IF (ID2.EQ.0) THEN  !SECOND RESIDUE NOT ALIGNED => NOT ALIGNED TO HELIX
                HELIX_B=.FALSE.
              ELSE
                IF (MEMPRE(ID2,CURTAB) .EQ. 1) THEN 
                   HELIX_B=.TRUE.    !SECOND RESIDUE ALIGNED TO HELIX
                ELSE
                   HELIX_B=.FALSE.    !SECOND RESIDUE *NOT* ALIGNED TO HELIX
                ENDIF
              ENDIF
              
              IF (HELIX_A.OR.HELIX_B) THEN
                 HELIX = 1
              ELSE
                 HELIX = 2
              ENDIF

              
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC             CHECK IF MATCH FILE CORRUPTED   CCCCCCCCC

              IF( (ID1.GT.NMRSS(I_MEM)).OR.(ID1.LT.0) ) THEN
                  WRITE(SO,*) 'MATCH FILE CORRUPTED'
                  WRITE(SO,*) 'MEMORY NUMBER',I_MEM,'NAME',PROTNM(I_MEM),'RESIDUE',ID3
                  WRITE(SO,*) 'ID1=',ID1
                  WRITE(SO,*) 'ID3=',ID3
                  WRITE(SO,*) (MATCH(I),I=1,NMRES)
                  STOP
              ENDIF
              IF( (ID2.GT.NMRSS(I_MEM)).OR.(ID2.LT.0).OR.
     *              ( (ID2-ID1.LE.0) .AND. (ID2*ID1.NE.0) ) )THEN
                  WRITE(SO,*) 'MATCH FILE CORRUPTED'
                  WRITE(SO,*) 'MEMORY NUMBER',I_MEM,'NAME',PROTNM(I_MEM),'RESIDUE',ID4
                  WRITE(SO,*) 'ID2=',ID2
                  WRITE(SO,*) 'ID1=',ID1
                  STOP
              ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


            IF ( IRES(ID3).EQ.8 .AND. IRES(ID4).EQ.8 ) THEN
               TAB_FOR_CONTACT=1
            ELSEIF ( IRES(ID3).EQ.8 .AND. IRES(ID4).NE.8 ) THEN
               TAB_FOR_CONTACT=2
            ELSEIF ( IRES(ID3).NE.8 .AND. IRES(ID4).EQ.8 ) THEN
               TAB_FOR_CONTACT=3
            ELSE
               TAB_FOR_CONTACT=4
            ENDIF                 !DECIDE TABLE FOR CONTACT POTENTIAL   
                                  !--USUALLY BETA-BETA (IE 4) BUT MUST
                                  ! ALLOW FOR GLYCINES
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCC REASONS TO EXIT LOOP OVER INTS (CONSIDER FIRST ONLY IF NOT CONTACT-INT)

            IF ( .NOT. ( ID4-ID3.GT.12 .AND. I_MEM.EQ.1 .AND.
     *                  CURTAB.EQ.TAB_FOR_CONTACT .AND. N_LETTERS.EQ.4
     *                  .AND. I_3RD_IS_CONTACT ) ) THEN


!             ONLY INTERACTIONS BETWEEN RESIDUES OF CERTAIN SEPARATION
              IF ((ID4-ID3).LT.MIN_SEQ_SEP) GOTO 519

C             NO POTENTIAL FOR BETA-CARBON ON GLYCINE 
               IF( (CDTYP1.EQ.2).AND.(IRES(ID3).EQ.8) )GO TO 519
               IF( (CDTYP2.EQ.2).AND.(IRES(ID4).EQ.8) )GO TO 519


!              IF NOT ALIGNED, GO TO NEXT INTERACTION
              IF ( (ID1.EQ.0) .OR. (ID2.EQ.0) ) GOTO 519

!             IF ALIGN TO GLYCINE AND ON BETA CARBON THEN SKIP
              IF( (CDTYP1.EQ.2).AND.(JRES(ID1).EQ.8) )GO TO 519
              IF( (CDTYP2.EQ.2).AND.(JRES(ID2).EQ.8) )GO TO 519

! EXIT LOOP IF CONTACT INTERACTION AND HAVE A GLYCINE AND GLY_CON IS FALSE
            ELSEIF (.NOT.GLY_CON) THEN 
              IF ( (IRES(ID3).EQ.8).OR.(IRES(ID4).EQ.8) ) GOTO 519
            ENDIF

CCCCCCCCCCC      END OF REASONS TO EXIT LOOP OVER INTS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


               DIST_IJ = DSQRT((YWORK(ID1,1,CDTYP1) -
     *                            YWORK(ID2,1,CDTYP2))**2 +
     *                           (YWORK(ID1,2,CDTYP1) -
     *                            YWORK(ID2,2,CDTYP2))**2 +
     *                           (YWORK(ID1,3,CDTYP1) -
     *                            YWORK(ID2,3,CDTYP2))**2 )

              DIST_PROX = DSQRT((YWORK(ID1,1,1) -
     *                            YWORK(ID2,1,1))**2 +
     *                           (YWORK(ID1,2,1) -
     *                            YWORK(ID2,2,1))**2 +
     *                           (YWORK(ID1,3,1) -
     *                            YWORK(ID2,3,1))**2 ) 

C              IF ( (CURTAB.EQ.4).AND.(I_IXN.EQ.101) ) THEN
C                WRITE(SO,*) 'DISTANCE IS',DIST_IJ
C                WRITE(SO,*) 'CDTYP1 AND 2 ARE',CDTYP1,CDTYP2
C              ENDIF


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SET CATEGORY OF INTERACTION BASED ON DISTANCE IN SEQ AND IN SPACE.

              IF (.NOT.I_ALT_PROX) THEN

                   IF (DIST_IJ.LT.SRCUT) THEN
                     IF ((ID4-ID3).LT.MINMR) THEN
                       IPROX=1
                     ELSE
                       IPROX=2
                     ENDIF
                   ELSE
                     IPROX=3
                   ENDIF

C      MIKE CHANGED THIS SO SHORT RANGE IN SEQ PAIRS ACTUALLY
C      GO INTO THE LONG RANGE IN SEQ AND SPACE CATEGORY IF THEY
C      ARE SEPARATED BY MORE THAN SRCUT ANGSTROMS
C      NOTE ALSO, THAT I APPEAR TO SET THIS FOR EACH TABLE, WHILE
C      COREY CHOSES CLASS BASED ONLY ON ALPHA-ALPHA DIST (+PROX)
C      THIS IS PROBABLY NOT A BIG DEAL, BUT WHEN RUNNING USING
C      HIS GAMMAS, SWITCH I_ALT_PROX FLAG TO TRUE AND WILL GET 
C      THE FOLLOWING PROXIMITY CLASSES


              ELSEIF (.NOT.I_3RD_IS_CONTACT) THEN
 
                   
                   IPROX = 4
                   IF ((ID4-ID3).LT.MINMR) IPROX=1
                   IF ( ((ID4-ID3).GE.MINMR)
     *             .AND. ((ID4-ID3).LE.MAXMR)) IPROX=2
                   IF ( ((ID4-ID3) .GT. MAXMR) .AND.
     *             (DIST_PROX .LT. ALT_PROX_CUT) ) IPROX = 3

              ELSE

                   IPROX = 4
                   IF ((ID4-ID3).LT.MINMR) IPROX=1
                   IF ( ((ID4-ID3).GE.MINMR)
     *             .AND. ((ID4-ID3).LE.MAXMR)) IPROX=2
                   IF ( ((ID4-ID3) .GT. MAXMR).AND.(I_MEM.EQ.1)
     *               .AND. CURTAB.EQ.TAB_FOR_CONTACT ) IPROX = 3

              ENDIF

C      RANDOM INTERACTIONS ARE FOR SEQUENCE SEPARATIONS OF
C      RAN_MIN_SEQ_DIST AND ABOVE

              IF ((ID4-ID3).LT.RAN_MIN_SEQ_DIST) THEN
                IPROX_RAN=.FALSE.
              ELSE
                IPROX_RAN=.TRUE.
              ENDIF   

              IF ( I_ALT_PROX.AND.(IPROX.EQ.4) ) GOTO 519



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C               FIND INDIVIDUAL CHARGES AND TOTAL CHARGE CONTRIBUTION
C     USED TO BE DONE IN SUBROUTINE CHARGE, BUT NOW DONE HERE
C     BECAUSE SO SHORT 

               TOTCHG=0.0D0
          IF (I_3RD_IS_CONTACT.AND.IPROX.EQ.3) THEN
           DO I_WELL=1,NUM_WELL
            DO ISEQ=1,NUMSEQ_AMC
              TOTCHG(I_WELL)= TOTCHG(I_WELL)-
     *        QCHRG(TGSEQUENCES_AMC(ID3,ISEQ),TGSEQUENCES_AMC(ID4,ISEQ),
     *        JRES(ID1),JRES(ID2),I_WELL+2)
                   ENDDO
              TOTCHG(I_WELL)=TOTCHG(I_WELL)*
     *                          LONG_NFACTOR(I_WELL)
                 ENDDO
               ELSE
              IF (FOUR_PLUS .AND. IPROX .LE. 2) THEN
          DO ISEQ=1,NUMSEQ_AMC
          TOTCHG(1)=TOTCHG(1)-
     *     QCHRG2(TGSEQUENCES_AMC(ID3,ISEQ),TGSEQUENCES_AMC(ID4,ISEQ),JRES(ID1),JRES(ID2),HELIX,IPROX)
          ENDDO
                 ELSE
          DO ISEQ=1,NUMSEQ_AMC
          TOTCHG(1)=TOTCHG(1)-
     *     QCHRG(TGSEQUENCES_AMC(ID3,ISEQ),TGSEQUENCES_AMC(ID4,ISEQ),JRES(ID1),JRES(ID2),IPROX)
          ENDDO
                 ENDIF
                 DO I_WELL=2,NUM_WELL
                   TOTCHG(I_WELL)=0.0D0
                 ENDDO
               ENDIF

	       IF (DIST_CUT) THEN
                IF ( (IPROX .EQ. 2) .AND.(DIST_IJ .GT. 6.5)) THEN
                  TOTCHG=0.0D0
                ENDIF
                ENDIF

               TOTCHG=TOTCHG/REAL(NUMSEQ_AMC)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C              FIND R-DEPENDENCE OF GAUSSIAN POTENTIAL
               DELT_SAFE=DELTZ(I_IXN,CURTAB)
               IF( WELSCL )THEN
                  CALL GAUSSV(MAXSZ,DIST_IJ,DELTZ(I_IXN,CURTAB),GAUSSP,RGRID)
         IF ( (GO_CON) .AND. (DIST_IJ .GT. GO_CON_DIST)) THEN
                   GAUSSP=0.0D0
          ENDIF

               ELSE
                  CALL GAUSSW(MAXSZ,DIST_IJ,0.25D0*DIST_IJ,GAUSSP,RGRID)
               ENDIF

               IF (RAN_FORCE.AND. (CURTAB.EQ.1).AND.(IPROX_RAN)) THEN
                  CALL GAUSSV(MAXSZ,R_RAN,DELT_SAFE,GAUSSQ,RGRID)
               ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C              COMBINE SEQUENCE-STRUCTURE TERMS

               DO 507 I507=1,MAXSZ

                 IF ( I_3RD_IS_CONTACT .AND.(IPROX .EQ. 3) ) THEN
                   GAUSSP(I507) = 0.0D0
                   DO I_WELL=1,NUM_WELL
                     GAUSSP(I507)= GAUSSP(I507)+TOTCHG(I_WELL)*
     *               0.25D0*((1.0D0 + TANH(7.0D0*(RGRID(I507)-R_MIN(I_WELL) )))
     *               *(1.0D0+TANH(7.0D0*(R_MAX(I_WELL)-RGRID(I507)))))
                   ENDDO
                 ELSE
                   GAUSSP(I507)=TOTCHG(1)*GAUSSP(I507)
                 ENDIF
 
  507          CONTINUE

C                 SET UP POTENTIAL
     
                    DO 522 I522=1,MAXSZ
                        VPOTNT(I522,I_IXN,CURTAB)=VPOTNT(I522,I_IXN,CURTAB) + GAUSSP(I522)
  522               CONTINUE

                  IF ((CURTAB.EQ.1) .AND. RAN_FORCE.AND. 
     *                                   (IPROX_RAN)) THEN
                    DO 588 I588=1,MAXS
                        GAUSSQ(I588)=GAUSSQ(I588)*RAN_NUM(I_IXN)
                        VPOTNT(I588,I_IXN,1)=VPOTNT(I588,I_IXN,1) +GAUSSQ(I588)
 588               ENDDO
                  ENDIF                       

C                 INCLUDE R-DEPENDENT FORCE TERM=
C                 -(1/R)*DV(R)/DR

                  DELT2=2.0D0*DELTZ(I_IXN,CURTAB)
                  DO 504 I504=1,MAXSZ

                    IF (I_3RD_IS_CONTACT.AND.IPROX .EQ. 3) THEN

                      GAUSSP(I504)=0.0D0
                      DO I_WELL=1,NUM_WELL
                        FORCE_TERM(1)=
     *                    TANH(7.0D0*(R_MAX(I_WELL)-RGRID(I504)))**2  ! G
                        FORCE_TERM(2)=
     *                    -TANH(7.0D0*(RGRID(I504)-R_MIN(I_WELL)))**2   ! U
                        FORCE_TERM(3)= 
     *                    TANH(7.0D0*(R_MAX(I_WELL)-RGRID(I504)))    !G
                        FORCE_TERM(4)= 
     *                    -TANH(7.0D0*(R_MAX(I_WELL)-RGRID(I504)))*    !GU^2
     *                      TANH(7.0D0*(RGRID(I504)-R_MIN(I_WELL)))**2
                        FORCE_TERM(5)= 
     *                    -TANH(7.0D0*(RGRID(I504)-R_MIN(I_WELL)))    !U
                        FORCE_TERM(6)= 
     *                    TANH(7.0D0*(RGRID(I504)-R_MIN(I_WELL)))*  !UG^2
     *                    TANH(7.0D0*(R_MAX(I_WELL)-RGRID(I504)))**2


                       GAUSSP(I504)=GAUSSP(I504)
     *                      -7.0D0*(TOTCHG(I_WELL)/4.0)* (
     *                      FORCE_TERM(1) + FORCE_TERM(2)
     *                      + FORCE_TERM(3) + FORCE_TERM(4) +
     *                      FORCE_TERM(5) + FORCE_TERM(6) )/RGRID(I504)
                      ENDDO

                    ELSE
                      GAUSSP(I504)=DELT2*GAUSSP(I504)*
     *                         (1.0D0 - (DIST_IJ/
     *                          RGRID(I504)) )

                    ENDIF

  504             CONTINUE

C                 SET UP FORCE TERM
    
                    DO 523 I523=1,MAXSZ
                      FORSE(I523,I_IXN,CURTAB)=FORSE(I523,I_IXN,CURTAB) + GAUSSP(I523)
  523               CONTINUE

                  IF (CURTAB.EQ.1 .AND. RAN_FORCE
     *                          .AND.(IPROX_RAN)) THEN
                    DO 589 I589=1,MAXS
                      GAUSSQ(I589)=DELT2*GAUSSQ(I589)*
     *                           (1.0D0 - (R_RAN/RGRID(I589)) )

                        FORSE(I589,I_IXN,1)=FORSE(I589,I_IXN,1) + GAUSSQ(I589)
                      
 589                ENDDO
                  ENDIF                       

C            END LOOP OVER INSERTIONS AND DELETIONS

  519    CONTINUE                        ! END OF LOOP OVER INTERACTIONS
  509 CONTINUE                                ! END OF LOOP OVER MEMORIES

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         PRINTS TEST POTENTIAL + FORCE (+INTEGRATED FORCE)
 
       IF (I_V_TEST) THEN

         I_TEST_TAB=CURTAB
         I_TEST_IXN=IXN_FROM_SITE(TEST_SITE(1),TEST_SITE(2),I_TEST_TAB)
         IF ( (TEST_SITE(1).NE.ILONG(I_TEST_IXN,1,I_TEST_TAB))
     *    .OR.(TEST_SITE(2).NE.ILONG(I_TEST_IXN,2,I_TEST_TAB)) )
     *      THEN
            WRITE(SO,*) 'ERROR IN GASPOT'
            WRITE(SO,*) 'TEST SITES FOR TABLE',I_TEST_TAB,'WRONG'
            WRITE(SO,*) 'SHOULD BE',TEST_SITE(1),TEST_SITE(2)
            WRITE(SO,*) 'FIND FROM LOOKING UP INTERACTION THEY ARE'
            WRITE(SO,*) ILONG(I_TEST_IXN,1,I_TEST_TAB),
     *                    ILONG(I_TEST_IXN,1,I_TEST_TAB)
            GOTO 1111
          ENDIF

         VPOTINT_TEST=0.0D0

         COUNT=TEST_SITE(1)
         CALL NUM_TO_CHAR(COUNT,CCOUNT1)
         COUNT=TEST_SITE(2)
         CALL NUM_TO_CHAR(COUNT,CCOUNT2)
         COUNT=I_TEST_TAB
         CALL NUM_TO_CHAR(COUNT,CCOUNT3)
         OPEN(OVGASPOT,FILE='V_GASPOT_'//TRIM(CCOUNT1)//'_'//
     *    TRIM(CCOUNT2)//'.'//TRIM(CCOUNT3),STATUS='UNKNOWN',
     *    ACTION='WRITE')

       WRITE(OVGASPOT,*) 'R,VPOTNT,FORSE, INT FORCE FOR IXN NUM,TABLE'
     *                                  ,I_TEST_IXN,I_TEST_TAB
         WRITE(OVGASPOT,*) 'SITES',ILONG(I_TEST_IXN,1,I_TEST_TAB)
     *                               ,ILONG(I_TEST_IXN,2,I_TEST_TAB)
         DO I=MAXS,1,-1
           VPOTINT_TEST(I)=VPOTINT_TEST(I+1)+
     *                         (1.0D0*RINC(I_TEST_IXN,I_TEST_TAB))*0.5D0*
     *                  (FORSE(I+1,I_TEST_IXN,I_TEST_TAB)*FLOAT(I+1)+
     *                    FORSE(I,I_TEST_IXN,I_TEST_TAB)*
     *                      FLOAT(I))*RINC(I_TEST_IXN,I_TEST_TAB)
         ENDDO
         DO I=1,MAXS+1
           WRITE(OVGASPOT,*) FLOAT(I)*RINC(I_TEST_IXN,I_TEST_TAB),
     *                   VPOTNT(I,I_TEST_IXN,I_TEST_TAB),
     *     (FORSE(I,I_TEST_IXN,I_TEST_TAB)*RINC(I_TEST_IXN,I_TEST_TAB))
     *                          *FLOAT(I),VPOTINT_TEST(I)
         ENDDO    
         CLOSE(OVGASPOT)
1111   ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DO I_IXN=1,INDXT
         VPOTNT(MAXS+1,I_IXN,CURTAB)=0.0D0


       DO I=MAXS,0,-1

         M=(1.0D0/RINC(I_IXN,CURTAB))*(FORSE(I+1,I_IXN,CURTAB)-FORSE(I,I_IXN,CURTAB))

         C=FORSE(I,I_IXN,CURTAB)*FLOAT(I+1)-FORSE(I+1,I_IXN,CURTAB)*FLOAT(I)

         VPOTNT(I,I_IXN,CURTAB)=VPOTNT(I+1,I_IXN,CURTAB)+
     *     M*RINC(I_IXN,CURTAB)*RINC(I_IXN,CURTAB)*RINC(I_IXN,CURTAB)*
     *      (FLOAT( (I+1)*(I+1)-I)-2.0D0/3.0D0)  +
     *     C*RINC(I_IXN,CURTAB)*RINC(I_IXN,CURTAB)*(FLOAT(I)+0.5D0)

       ENDDO
       ENDDO

C     WRITE OUT PROTEINS USED IN GENERATING POTENTIAL

      IF( PASSI )THEN
C         WRITE(OARCHV,177)PROTNM(0),INT(WORK8(1))
  177    FORMAT(/'TARGET PROTEIN ',A5,1X,I5)
C         WRITE(OARCHV,111)NUMMEM
  111    FORMAT('# PROTEINS USED IN GENERATING POTENTIAL ',I3)
C         WRITE(OARCHV,114)(PROTNM(I1),INT(WORK8(I1+1)),I1=1,NUMMEM)
  114    FORMAT(5(A5,1X,I5,2X))
      ENDIF

C     INVERT R-GRID INCREMENT

      DO 518 I518=1,INDXT
         RINCINV(I518,CURTAB)=1.0D0/RINC(I518,CURTAB)
         RINCSQ(I518,CURTAB)=RINC(I518,CURTAB)*
     *                              RINC(I518,CURTAB)
  518 CONTINUE

C     --- DIAGNOSTICS ---


         DO 534 I_IXN=1,INDXT

            ID1=ILONG(I_IXN,1,CURTAB)
            ID2=ILONG(I_IXN,2,CURTAB)

C           CHECK IF ENDPOINTS OF FORCE SMALL ENOUGH

C            IF( ABS(FORSE(MAXS,I_IXN,CURTAB)).GT.5.0E-5 )
C     *         WRITE(OARCHV,121)CURTAB,I_IXN,ID1,ID2,
C     *               FORSE(MAXS,I_IXN,CURTAB)
C  121          FORMAT('F > 5.0E-05:TABLE ',I2,' IXN ',I5,
C     *                ' BETWEEN ',2(I3,1X),' F(MAXS) ',
C     *                1PE10.3)

            IF( (ID1.EQ.0).AND.(ID2.LE.12).AND.(CDTYP1.EQ.2)  )THEN

C              CHECK IF POTENTIAL HAS BEEN CORRECTLY
C              GENERATED

C              GENERATE R-GRID

               DO 503 I503=1,MAXS
                 RGRID(I503)=FLOAT(I503)*RINC(I_IXN,CURTAB)
  503          CONTINUE

C               WRITE(OARCHV,110)CURTAB,CDTYP1,CDTYP2,ID1,ID2,
C     *                          RINCINV(I_IXN,CURTAB)
  110          FORMAT(/'TEST FORCE: TABLE NUMBER ',I3,
     *                ' COORD1 ',I1,' COORD2 ',I1,
     *                ' SITES',2(1X,I3),' RINCINV ',1PE10.3)
               DO 532 I532=2,MAXS-1

                  WORK2(I532)=0.5D0*(VPOTNT(I532-1,I_IXN,CURTAB)-
     *                             VPOTNT(I532+1,I_IXN,CURTAB))/
     *                            (RGRID(I532-1)-RGRID(I532))

                  IF( FORSE(I532,I_IXN,CURTAB)*RGRID(I532)
     *                    .NE.0.0D0 )
     *               RNORM=-WORK2(I532)/
     *                  (FORSE(I532,I_IXN,CURTAB)*RGRID(I532))

  112             FORMAT(I3,5(1X,1PE10.3))
  532          CONTINUE

            ENDIF

  534    CONTINUE

C     --- END DIAGNOSTICS ---

          WRITE(SO,*) 'END OF GASPOT '
C     ---------------------- DONE -----------------------

      RETURN
      END
