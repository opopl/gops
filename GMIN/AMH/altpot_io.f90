      SUBROUTINE DEFAULT_ALT()

        USE GLOBALS_ALT, ONLY : ALTPOTFLAG,KAPPA_ALT,TRESHOLD_ALT,KAPPA_WELL,KAPPA_OB, & 
             DEBUG_NUMER_FORCES,OUTFILE1_ALT, OUTPUT_STEP_ALT, OUTPUT_STEPSIZE_ALT, & 
             DO_SEND_OUTPUT,ONEBODYFLAG, ONEBODY_TYPE , TIMINGS_FILE_ALT, ACCUMULATED_TIME

        USE ALTPOT_INTERFACES, ONLY: INIT_ONEBODY

        IMPLICIT NONE
        INTEGER I

        TRESHOLD_ALT=0.70D0
        KAPPA_ALT=7.00D0
        KAPPA_WELL=7.00D0
        KAPPA_OB=7.00D0
        DEBUG_NUMER_FORCES=.FALSE.
        ALTPOTFLAG=0
        OUTPUT_STEP_ALT=0
        OUTPUT_STEPSIZE_ALT=10
        DO_SEND_OUTPUT=.FALSE.
        ONEBODYFLAG=0
        ONEBODY_TYPE=1
        ACCUMULATED_TIME=0.00D0


        CALL INIT_ONEBODY()

        DO I=1,10,1
           OUTFILE1_ALT(I)=102+I
        ENDDO

!        OPEN(UNIT=OUTFILE1_ALT(1),FILE='CHECK_ALT_FILE_1',STATUS='UNKNOWN')
!        OPEN(UNIT=OUTFILE1_ALT(2),FILE='CHECK_ALT_FILE_2',STATUS='UNKNOWN')
!        OPEN(UNIT=OUTFILE1_ALT(3),FILE='CHECK_ALT_FILE_3',STATUS='UNKNOWN')
!        OPEN(UNIT=OUTFILE1_ALT(4),FILE='CHECK_ALT_FILE_4',STATUS='UNKNOWN')
!        OPEN(UNIT=OUTFILE1_ALT(5),FILE='OB_DENSITY',STATUS='UNKNOWN')
!        OPEN(UNIT=OUTFILE1_ALT(6),FILE='OB_ENERGY',STATUS='UNKNOWN')
!        OPEN(UNIT=OUTFILE1_ALT(7),FILE='RESIDUE_DENSITY',STATUS='UNKNOWN')
!        OPEN(UNIT=TIMINGS_FILE_ALT,FILE='ALT_TIMINGS',STATUS='UNKNOWN')

        RETURN
      END SUBROUTINE DEFAULT_ALT

      ! READ_INPUT CALLS DEFAULT_ALT AND READ_ALTGAMMA AND READS INPUT, 
      ! SHOULD BE CALLED AFTER INITIL
      SUBROUTINE READ_INPUT_ALT()

        USE GLOBALS_ALT, ONLY : ALTPOTFLAG,KAPPA_ALT,TRESHOLD_ALT,KAPPA_WELL,ONEBODYFLAG &
             , ONEBODY_TYPE,KAPPA_OB,ALIMITS_OB, DEBUG_NUMER_FORCES

        USE ALTPOT_INTERFACES, ONLY: DEFAULT_ALT
        USE AMHGLOBALS,  ONLY : SO

        IMPLICIT NONE
        INTEGER NINPUT,IOER,N_OUTPUT,DUMI
        PARAMETER(NINPUT=130)
        PARAMETER(N_OUTPUT=131)
        CHARACTER STR200*200,INPUTFILE*200,PARAMCHECKFILE*200
        INPUTFILE='INPUT_AMH'
        PARAMCHECKFILE='AMW_OUTPUT'
        OPEN(UNIT=NINPUT,FILE=INPUTFILE,STATUS='OLD')
!        OPEN(UNIT=N_OUTPUT,FILE=PARAMCHECKFILE,STATUS='UNKNOWN')

        ! INITIALIZE ALTPOT VALUES
        CALL DEFAULT_ALT()


        !++++++ START INPUTLOOP
!        WRITE(6,*) 'START DEFAULT ALT'
100     READ(NINPUT,'(A)',END=999)STR200

        IF(INDEX(STR200,'ALTPOTFLAG').GT.0)THEN
           BACKSPACE(NINPUT)
           READ(NINPUT,*,IOSTAT=IOER,ERR=900) ALTPOTFLAG(1),ALTPOTFLAG(2)
        ELSEIF(INDEX(STR200,'KAPPA_ALT').GT.0)THEN
           BACKSPACE(NINPUT)
           READ(NINPUT,*,IOSTAT=IOER,ERR=900) KAPPA_ALT
        ELSEIF(INDEX(STR200,'TRESHOLD_ALT').GT.0)THEN
           BACKSPACE(NINPUT)
           READ(NINPUT,*,IOSTAT=IOER,ERR=900) TRESHOLD_ALT
        ELSEIF(INDEX(STR200,'KAPPA_WELL').GT.0)THEN
           BACKSPACE(NINPUT)
           READ(NINPUT,*,IOSTAT=IOER,ERR=900) KAPPA_WELL
        ELSEIF(INDEX(STR200,'ONEBODYFLAG').GT.0)THEN
           BACKSPACE(NINPUT)
           READ(NINPUT,*,IOSTAT=IOER,ERR=900) ONEBODYFLAG,ONEBODY_TYPE
        ELSEIF(INDEX(STR200,'KAPPA_OB').GT.0)THEN
           BACKSPACE(NINPUT)
           READ(NINPUT,*,IOSTAT=IOER,ERR=900) KAPPA_OB
        ELSEIF(INDEX(STR200,'DENSITYLIMITS_OB').GT.0)THEN
           BACKSPACE(NINPUT)
           READ(NINPUT,*,IOSTAT=IOER,ERR=900)ALIMITS_OB
        ELSEIF(INDEX(STR200,'DEBUG_NUMER_FORCES').GT.0)THEN
           BACKSPACE(NINPUT)
           READ(NINPUT,*,IOSTAT=IOER,ERR=900)DUMI
           IF(DUMI.GT.0)DEBUG_NUMER_FORCES=.TRUE.
        ENDIF
        GOTO 100
        !------ END INPUTLOOP
!        RETURN

        !  REPORT VALUES TO FILE PARAMCHECK_ALT:
999      WRITE(6,*)
!999     WRITE(N_OUTPUT,*)'ALTPOTFLAGS ARE : ',ALTPOTFLAG(1),ALTPOTFLAG(2)
!        WRITE(N_OUTPUT,*)'KAPPA_ALT IS : ',KAPPA_ALT
!        WRITE(N_OUTPUT,*)'TRESHOLD_ALT IS : ',TRESHOLD_ALT
!        WRITE(N_OUTPUT,*)'KAPPA_WELL IS : ',KAPPA_WELL
!1        WRITE(N_OUTPUT,*)'ONEBODYFLAG,ONEBODY_TYPE IS : ',ONEBODYFLAG,ONEBODY_TYPE
!        WRITE(N_OUTPUT,*)'KAPPA_OB IS : ',KAPPA_OB
!        WRITE(N_OUTPUT,*)'ALIMITS ARE : ',ALIMITS_OB

900     IF(IOER.NE.0)THEN
           WRITE(*,'(A)')'PROBLEM WITH READING INPUT SUBROUTINE: READ_INPUT_ALT'
        ENDIF

        CLOSE(NINPUT)
!        CLOSE(N_OUTPUT)

        RETURN


      END SUBROUTINE READ_INPUT_ALT

      !------------------------------------------------

      !------------------------------------

      SUBROUTINE READ_ALTGAMMA()

        ! READS IN ALTGAMMA AND ONEBODY_GAMMA FROM GAMMA.DAT
        ! FORMAT IN GAMMA.DAT IS EXPECTED TO BE 
        ! GAMMA1   GAMMA2  WELL  KEYWORD=ALTPOT
        ! ALTPOT IS TO BE USED FOR ANISOTROPIC POTENTIAL
        !      USE AMHGLOBALS,  ONLY : MAXSIZ,MAXCRD,IRES,R_MIN,R_MAX,SORT_NON_ADD,
        !            GAMMA_NON_ADD,CLASS_OF_RES_2

        USE AMHGLOBALS,  ONLY : SO,N_LETTERS_CON
        USE GLOBALS_ALT, ONLY : ALTGAMMA,ALTPOTFLAG,MAX_WELL_ALT,ONEBODY_GAMMA,&
               ONEBODYFLAG

        USE ALTPOT_INTERFACES, ONLY: LONGSCALE

        IMPLICIT NONE 

        INTEGER I,J,IWELL,ALTPOT_INDEX,DID_READ_GAMMA(MAX_WELL_ALT)
        INTEGER READ_STATUS,DID_READ_ONEBODY_GAMMA,ONEBODY_INDEX
            DOUBLE PRECISION X1,X2,X3
        CHARACTER  STR200*200,GAMMAFIL*200


        DID_READ_GAMMA=0
        DID_READ_ONEBODY_GAMMA=0

        GAMMAFIL='GAMMA.DAT'
        OPEN(UNIT=130,FILE=GAMMAFIL,IOSTAT=READ_STATUS,STATUS='OLD')
!        OPEN(UNIT=131,FILE='ALT_GAMMACHECK',IOSTAT=READ_STATUS,STATUS='UNKNOWN')
!        OPEN(UNIT=132,FILE='ONEBODY_GAMMACHECK',IOSTAT=READ_STATUS,STATUS='UNKNOWN')

        !  LOOP AND SEARCH FOR KEYWORD ALTPOT IN GAMMA.DAT 
        !                                              !!!!!!!!!!!!!!!!!!! ! START READING LOOP
100     READ (130,'(A)',ERR=99,END=999)STR200                
        ALTPOT_INDEX=INDEX(STR200,'ALTPOT')
        ONEBODY_INDEX=INDEX(STR200,'ONEBODY')
        IF(ALTPOT_INDEX.GT.0)THEN                   ! ALTPOT
           BACKSPACE(130)
           DO I=1,N_LETTERS_CON,1
              DO J=I,N_LETTERS_CON,1
                READ (130,*)X1,X2,IWELL
                ALTGAMMA(I,J,1,IWELL)=X1
                ALTGAMMA(I,J,2,IWELL)=X2
                ALTGAMMA(J,I,1,IWELL)=ALTGAMMA(I,J,1,IWELL)
                ALTGAMMA(J,I,2,IWELL)=ALTGAMMA(I,J,2,IWELL)
!                WRITE(131,*)ALTGAMMA(I,J,1,IWELL),ALTGAMMA(I,J,2,IWELL), &
!                IWELL
              ENDDO
           ENDDO
           DID_READ_GAMMA(IWELL)=1
        ENDIF

        IF(ONEBODY_INDEX.GT.0)THEN                   ! ONEBODY
           BACKSPACE(130)
           DO I=1,N_LETTERS_CON,1
              READ (130,*)X1,X2,X3
              ONEBODY_GAMMA(I,1)=X1
              ONEBODY_GAMMA(I,2)=X2
              ONEBODY_GAMMA(I,3)=X3
!              WRITE(132,*)ONEBODY_GAMMA(I,1),ONEBODY_GAMMA(I,2), & 
!                          ONEBODY_GAMMA(I,3),I
           ENDDO
           DID_READ_ONEBODY_GAMMA=1
        ENDIF

        GOTO 100                                            
        !                                            !!!!!!!!!!!!!!!!!! ! END READING LOOP

999     CLOSE (130)
!        CLOSE (131)
!        CLOSE (132)

        DO I=1,MAX_WELL_ALT    ! ERRORCHECK INPUT
           IF((ALTPOTFLAG(I).GT.0).AND.(DID_READ_GAMMA(I).EQ.0))THEN
              WRITE(*,*)'ALTPOTFLAG  ',I,' SET BUT NOT PRESENT IN GAMMA.DAT'
              STOP
           ENDIF
        ENDDO

        IF((ONEBODYFLAG.GT.0) .AND. (DID_READ_ONEBODY_GAMMA.EQ.0))THEN
         WRITE(SO,*)'ONEBODYFLAG SET BUT NOT PRESENT IN GAMMA.DAT WILL USE HP_SCALE AS ONEBODYGAMMAS'
           DO I=1,N_LETTERS_CON,1
              WRITE(132,*)ONEBODY_GAMMA(I,1),ONEBODY_GAMMA(I,2),ONEBODY_GAMMA(I,3),I
           ENDDO
        ENDIF

        CALL LONGSCALE()          ! THE LARGE N FUDGEFACTOR

        !                               ! ERROR 

99      IF(READ_STATUS.NE.0)THEN   
           WRITE(SO,*)'SOMETHING WRONG WITH INPUT ALTERNATIVE POTENTIAL'
           STOP
        ENDIF

        RETURN

      END SUBROUTINE READ_ALTGAMMA

      !--------------------------------------------- 
      
        SUBROUTINE FINALIZE_ALT()

        USE GLOBALS_ALT, ONLY : OUTFILE1_ALT,TIMINGS_FILE_ALT

        IMPLICIT NONE

!        CLOSE (OUTFILE1_ALT(1))
!        CLOSE (OUTFILE1_ALT(2))
!        CLOSE (OUTFILE1_ALT(3))
!        CLOSE (OUTFILE1_ALT(4))
!        CLOSE (OUTFILE1_ALT(5))
!        CLOSE (OUTFILE1_ALT(6))
!        CLOSE (OUTFILE1_ALT(7))
!        CLOSE (TIMINGS_FILE_ALT)

        RETURN
      END SUBROUTINE FINALIZE_ALT


      !....


      SUBROUTINE SEND_OUTPUT_ALT(A,E_ALT,NMRES,E_OB)

        USE AMHGLOBALS,  ONLY: MAXSIZ
        USE GLOBALS_ALT, ONLY : OUTPUT_STEP_ALT, OUTPUT_STEPSIZE_ALT,OUTFILE1_ALT,DO_SEND_OUTPUT &
             ,COUNT_ALT,T_ALT,MAX_WELL_ALT,TRESHOLD_ALT,OB_DENSITY,OB_DNS_COUNT &
             , MAX_LETTERS,AMINOACIDS,TIMINGS_FILE_ALT,ACCUMULATED_TIME
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: A(MAXSIZ),E_ALT(2,MAX_WELL_ALT),E_OB(3)
        INTEGER, INTENT(IN) ::  NMRES
        INTEGER N_IN,I
   
       DOUBLE PRECISION RNMRES,E_ALT_SUM


        OUTPUT_STEP_ALT=OUTPUT_STEP_ALT+1
        N_IN=0
        RNMRES=DBLE(NMRES)

        !      STR200="# STEP T  N_IN/NMRES  E_ALT[ 11   12    21   22] END"
        !      IF(OUTPUT_STEP_ALT.EQ.1)WRITE(OUTFILE1_ALT(1),'(A)')STR200
        !      IF(OUTPUT_STEP_ALT.EQ.1)WRITE(OUTFILE1_ALT(1),'(A)')STR200(1:MIN(200,INDEX('END',STR200)))
!        IF(OUTPUT_STEP_ALT.EQ.1)THEN
!           OPEN(UNIT=120,FILE='OUTPUT_GUIDE',STATUS='UNKNOWN')
!           WRITE(120,*)'GUIDE FOR CHECK_ALT_FILES'
!1           WRITE(120,*)'CHECK_ALT_FILE_1: STEP   T  N_IN/NMRES  E_ALT[ 11   21    12   22]  E_ALT_TOT'
!           !           WRITE(120,*)'CHECK_ALT_FILE_2: STEP   T  N_IN/NMRES  A[I]  '
!           CLOSE(120)
!        ENDIF

        IF(DO_SEND_OUTPUT)THEN
           DO I=1,NMRES,1
              IF(A(I).GT.TRESHOLD_ALT)N_IN=N_IN+1
           ENDDO
           E_ALT_SUM=E_ALT(1,1)+E_ALT(1,2)+E_ALT(2,1)+E_ALT(2,2)
!           WRITE(OUTFILE1_ALT(1),1000)COUNT_ALT,T_ALT,DBLE(N_IN)/RNMRES,E_ALT(1,1),E_ALT(2,1) &
!                ,E_ALT(1,2),E_ALT(2,2),E_ALT_SUM
           !         REWIND(OUTFILE1_ALT(2))
!           WRITE(OUTFILE1_ALT(2),1000)COUNT_ALT,T_ALT,(A(I),I=1,NMRES)
!           WRITE(TIMINGS_FILE_ALT,800) COUNT_ALT,T_ALT,(ACCUMULATED_TIME(I),I=1,4),(ACCUMULATED_TIME(I),I=10,13)
!800        FORMAT("COUNT T  TIMINGS:[ALTPOT][PP][OB][NUM] [ALTPOT_UTIL][PP1][PP2][PP3]", & 
!                                       I8,1X,F8.3,4X,4(F8.3,1X),2X,4(F8.3,1X))

        ENDIF

        !      IF((OB_DNS_COUNT.GT.1.0) .AND. DO_SEND_OUTPUT)WRITE(OUTFILE1_ALT(5),1000)COUNT_ALT,T_ALT,OB_DENSITY/OB_DNS_COUNT
!        IF((OB_DNS_COUNT.GT.1.0) .AND. DO_SEND_OUTPUT)THEN
!           REWIND(OUTFILE1_ALT(5))
!           DO I=1,MAX_LETTERS,1
!              WRITE(OUTFILE1_ALT(5),*)1.0,OB_DENSITY(I,1),AMINOACIDS(I)
!              WRITE(OUTFILE1_ALT(5),*)2.0,OB_DENSITY(I,2),AMINOACIDS(I)
!              WRITE(OUTFILE1_ALT(5),*)3.0,OB_DENSITY(I,3),AMINOACIDS(I)
!              WRITE(OUTFILE1_ALT(5),2000)
!           ENDDO
!        ENDIF
!        IF(DO_SEND_OUTPUT)WRITE(OUTFILE1_ALT(6),1000)COUNT_ALT,T_ALT,E_OB


        IF(MOD(OUTPUT_STEP_ALT,OUTPUT_STEPSIZE_ALT).NE.0)RETURN

        !      WRITE(OUTFILE1_ALT(1),60)

        !60    FORMAT(20(1X,E12.5))
1000    FORMAT(I8,2X,200(F8.3,2X))
2000    FORMAT()

        RETURN

      END SUBROUTINE SEND_OUTPUT_ALT

      SUBROUTINE LONGSCALE()

        USE AMHGLOBALS,  ONLY:SO,  AB_C_OF_N_OLD,AB_C_OF_N_NEW,MAX_WELL, & 
              ALPHA_C_OF_N,NUM_WELL,NMRES,N_LETTERS_CON
        USE GLOBALS_ALT, ONLY : ALTGAMMA

        IMPLICIT NONE
        INTEGER I,J
            DOUBLE PRECISION LONG_NFACTOR(MAX_WELL),RNMRES


        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !C   SET LONG_NFACTOR WHICH IS SOME FUDGE FACTOR FOR THE
        !C   CONTACT POTENTIALS (IT PARTIALLY ALLOWS FOR THE EFFECT OF 
        !C   SIZE -- WHICH THE FRACTION OF CONTACTS AT A CERTAIN DISTANCE
        !C   DEPENDS ON). IT SHOULD BE THE SAME AS INCLUDED IN THE
        !C   OPTIMISATION. IN COREY'S CODE IT IS IN QCHRGMK.F INSTEAD.
        !C   JOHAN: ADDED LONG_NFACTOR FOR 2-WELLS 

        LONG_NFACTOR=1.0D0
        RNMRES=DBLE(NMRES)

        !C   START IF ALPHA_C_OF_N 
        IF ( ALPHA_C_OF_N ) THEN

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
        !C   END IF ALPHA_C_OF_N 

        IF (AB_C_OF_N_NEW) THEN
           RNMRES=DBLE(NMRES)
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
            LONG_NFACTOR(7)=( 0.18665D0*RNMRES) /(RNMRES*0.0107983D0 + 1.0D0)
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
           RNMRES=DBLE(NMRES)
           LONG_NFACTOR(1)=1.0D0/(RNMRES*0.0015D0 + 1.94D0)
           LONG_NFACTOR(2)=1.0D0/(RNMRES*0.0032D0 + 1.83D0)
           LONG_NFACTOR(3)=1.0D0/(RNMRES*0.022D0 + 7.77D0)

           IF (NUM_WELL.NE.3) THEN 
              WRITE(SO,*) 'NUMWELL MUST BE 3 FOR AB_C_OF_N_*OLD*',NUM_WELL
              STOP
           ENDIF

        ENDIF

        DO I=1,N_LETTERS_CON,1
           DO J=I,N_LETTERS_CON,1
              ALTGAMMA(I,J,1,1)=ALTGAMMA(I,J,1,1)*LONG_NFACTOR(1)
              ALTGAMMA(I,J,2,1)=ALTGAMMA(I,J,2,1)*LONG_NFACTOR(1)
              ALTGAMMA(I,J,1,2)=ALTGAMMA(I,J,1,2)*LONG_NFACTOR(2)
              ALTGAMMA(I,J,2,2)=ALTGAMMA(I,J,2,2)*LONG_NFACTOR(2)
              ALTGAMMA(J,I,1,1)=ALTGAMMA(I,J,1,1)
              ALTGAMMA(J,I,2,1)=ALTGAMMA(I,J,2,1)
              ALTGAMMA(J,I,1,2)=ALTGAMMA(I,J,1,2)
              ALTGAMMA(J,I,2,2)=ALTGAMMA(I,J,2,2)
           ENDDO
        ENDDO

        RETURN
      END SUBROUTINE LONGSCALE
