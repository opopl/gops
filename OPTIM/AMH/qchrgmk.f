C     --------------------- QCHRGMK ----------------------

      SUBROUTINE QCHRGMK(QCHRG,QCHRG2,OARCHV,N_LETTERS)

C     ---------------------------------------------------
C     QCHRGMK MAKES TABLE QCHRG
C     ---------------------------------------------------

      USE AMHGLOBALS,  ONLY: SO,MAX_LETTERS,I_3RD_IS_CONTACT,FOUR_PLUS,
     *  I_NON_ADD_CONTACT,SORT_NON_ADD,
     *  GAMMA_NON_ADD,NGAMMA_NON_ADD,CLASS_2,MAX_WELL,NUM_WELL,
     *  PARA_HB,ANTI_HB,ANTI_NHB,I_IGNORE_CHAIN_DIRN,MISMATCH,
     *  PARA_ONE,ANTI_ONE,N_LETTERS_CON, SRPLR,IGAMMA,IHBOND

      USE GLOBALS_ALT, ONLY:ALTPOTFLAG ! FOR ALTERNATIVE POTENTIAL

      IMPLICIT NONE

C
C	I_RANGES	3; NEAR IN SEQ.   -->     NEAR IN SPACE
C                          FAR  IN SEQ.   AND     NEAR IN SPACE
C                          FAR  IN SEQ.   AND     FAR  IN SPACE

         INTEGER N_LETTERS

C     INTERNAL VARIABLES:
                                                                                                     
         DOUBLE PRECISION QCHRG(21,21,20,20,MAX_WELL+2),OPTGAMMA(512),QCHRG2(21,21,20,20,2,2),GAMMA_TEMP

         INTEGER I500,I501,I502,I503,I504,I_RANGES,NGAMMA,NGAMMA_CHECK,OPEN_STATUS,
     *           READ_STATUS,ISIT1,ISIT2,N_GAMMAS_PER_CON_WELL,IC1,OARCHV                               

C     TSHEN
      INTEGER MYJ,MYIP,MYJP,MYK,MYKP,MYWEL
      DOUBLE PRECISION MYCHCH, MYCHPO

        INTEGER CLASS(20),CLASS_CON(20),ICL,JCL,IPCL,JPCL,
     *       GAMMAPT,SORT_2(2,2,2,2),CLASS_4(20),CLASS_20(20),
     *       SORT(MAX_LETTERS,MAX_LETTERS,MAX_LETTERS,MAX_LETTERS),
     *       GAMMAPT_CONTACT,SORT_CONTACT(MAX_LETTERS,MAX_LETTERS),
     *       N_RANGES,N_GAMMAS_IN_RANGE(MAX_WELL+2),I_WELL,I_WELL2,
     *       I_TYPE,ICL_CON,JCL_CON,
     *  SORT_SYM(MAX_LETTERS,MAX_LETTERS,MAX_LETTERS,MAX_LETTERS),
     *  SORT_PLUS(MAX_LETTERS,MAX_LETTERS,MAX_LETTERS,MAX_LETTERS,2)

         INTEGER I, J, IP, JP,COUNT

        DATA CLASS_4  /1,3,2,2,4,2,2,1,3,4,4,3,4,4,1,1,1,4,4,4/
        DATA CLASS_20 /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20/
        DATA SORT_2 /1,5,6,16,7,3,14,12,8,15,2,11,13,10,9,4/
        CLASS_2=(/2,1,1,1,2,1,1,2,1,2,2,1,2,2,1,1,1,2,2,2/)

C  SET CLASSES FOR DIFFERENT TYPES OF CODE
C  DUE TO HISTORY, GAMMAS ARE GROUPED DIFFERENTLY
C  FOR THE TWO-LETTER CODE COMPARED TO MULTI-LETTER
C  CODES. THE TWO-LETTER CODE IS THUS TREATED ON A
C  SLIGHTLY DIFFERENT FOOTING.

          WRITE(6,*)'N_LETTERS ',N_LETTERS 

          IF (N_LETTERS.EQ.2) THEN
            DO I=1,20                  !TWO-LETTER CODE
              CLASS(I)=CLASS_2(I)
            ENDDO
          ELSEIF (N_LETTERS.EQ.4) THEN
            DO I=1,20                  !FOUR-LETTER CODE
              CLASS(I)=CLASS_4(I)
            ENDDO
          ENDIF
                                                                                                     
         IF (N_LETTERS_CON.EQ.2) THEN  !CLASSES MAY BE DIFFERENT FOR CONTACT INTERACTIONS
             CLASS_CON(1:20)=CLASS_2(1:20)
          ELSEIF (N_LETTERS_CON.EQ.4) THEN
             CLASS_CON(1:20)=CLASS_4(1:20)
          ELSEIF (N_LETTERS_CON.EQ.20) THEN
             CLASS_CON(1:20)=CLASS_20(1:20)
          ENDIF
        IF (N_LETTERS.EQ.2) THEN

          DO I=1,2                  ! TWO-LETTER CODE
          DO J=1,2
          DO IP=1,2
          DO JP=1,2
          SORT(I,J,IP,JP)=SORT_2(I,J,IP,JP)
          ENDDO
          ENDDO
          ENDDO
          ENDDO
 
        ELSE

C        SORT_SYM IS FOR THE CASE WHERE FOUR_PLUS IS FALSE AND CHAIN DIRECTION IS IGNORED
                                                                                                     
          COUNT = 1
             DO I = 1,N_LETTERS
             DO J = I,N_LETTERS   ! NOTE J STARTS AT I (AND JP AT IP) SO LOOPING OVER ALL
             DO IP = 1,N_LETTERS  ! DISTINGUISHABLE PAIRS IN BOTH TARGET AND MEMORY
             DO JP = IP,N_LETTERS  !(A TOTAL OF N_LETTERS+(N_LETTERS -CHOOSE-2) IN EACH
               SORT_SYM(I,J,IP,JP) = COUNT
               SORT_SYM(J,I,JP,IP) = COUNT  ! SWAPPING BOTH PAIRS (IE (IJ-I'J' -> JI-J'I')
                                            ! GIVES SAME GAMMA IF CHAIN DIRECTION CONSIDERED IRRELEVANT
               IF ( (I.NE.J) .AND. (IP.NE.JP) ) COUNT=COUNT+1 !IF BOTH PAIRS DIFFERENT THEN SWAPPING IDENTITIES IN JUST
               SORT_SYM(J,I,IP,JP) = COUNT            !ONE PAIR WILL GIVE A NEW GAMMA, SO INCREMENT
               SORT_SYM(I,J,JP,IP) = COUNT
               COUNT = COUNT +1
             ENDDO
             ENDDO
             ENDDO
             ENDDO


          COUNT = 1                 !MULTI-LETTER CODE
          DO 2 I = 1,N_LETTERS
          DO 4 J = 1,N_LETTERS
          DO 6 IP = 1,N_LETTERS
          DO 8 JP = 1,N_LETTERS

          SORT(I,J,IP,JP) = COUNT
          COUNT = COUNT +1
8         CONTINUE
6         CONTINUE
4         CONTINUE
2         CONTINUE

           IF (FOUR_PLUS) THEN           !INCLUDE SECONDARY STRUCTURAL INFORMATION

           IF (I_IGNORE_CHAIN_DIRN.AND.(.NOT.MISMATCH)) THEN  ! HERE CONSIDER THE ORDER OF RESIDUES IN CHAIN
                                                              !  NOT TO EFFECT THEIR INTERACTION ENERGY
             COUNT = 1      
             DO I = 1,N_LETTERS
             DO J = I,N_LETTERS   ! NOTE J STARTS AT I (AND JP AT IP) SO LOOPING OVER ALL 
             DO IP = 1,N_LETTERS  ! DISTINGUISHABLE PAIRS IN BOTH TARGET AND MEMORY
             DO JP = IP,N_LETTERS  !(A TOTAL OF N_LETTERS+(N_LETTERS -CHOOSE-2) IN EACH 
             DO IC1  = 1,2
               SORT_PLUS(I,J,IP,JP,IC1) = COUNT
               SORT_PLUS(J,I,JP,IP,IC1) = COUNT  ! SWAPPING BOTH PAIRS (IE (IJ-I'J' -> JI-J'I') 
                                                    ! GIVES SAME GAMMA IF CHAIN DIRECTION CONSIDERED IRRELEVANT
               IF ( (I.NE.J) .AND. (IP.NE.JP) ) COUNT=COUNT+1 !IF BOTH PAIRS DIFFERENT THEN SWAPPING IDENTITIES IN JUST
               SORT_PLUS(J,I,IP,JP,IC1) = COUNT            !ONE PAIR WILL GIVE A NEW GAMMA, SO INCREMENT
               SORT_PLUS(I,J,JP,IP,IC1) = COUNT
               COUNT = COUNT +1
             ENDDO
             ENDDO
             ENDDO
             ENDDO
             ENDDO


           ELSEIF (MISMATCH.AND.I_IGNORE_CHAIN_DIRN) THEN !HERE IGNORE NOT ONLY CHAIN DIRECTION BUT ALSO
                                                          !ASSUME PAIRS THAT FOR EACH KIND OF TARGET PAIR HAVE ONLY TWO
                                                          !KINDS OF ALIGNMENT: EXACTLY RIGHT AND NOT EXACTLY RIGHT (IE MISMATCHED)
                 SORT_PLUS=0
                 COUNT = 1
                 DO I = 1,N_LETTERS
                 DO J = I,N_LETTERS
                 DO IC1 = 1,2
                 SORT_PLUS(I,J,I,J,IC1)=COUNT
                 SORT_PLUS(J,I,J,I,IC1)=COUNT
                 COUNT = COUNT + 1
                 ENDDO
                 ENDDO
                 ENDDO

                 COUNT = 21

                 DO I = 1,N_LETTERS
                 DO J = I,N_LETTERS
                 DO IC1 = 1,2
                        DO IP = 1,4
                        DO JP = 1,4
                        IF (IP .NE. I .OR. JP .NE. J) THEN
                         SORT_PLUS(I,J,IP,JP,IC1)=COUNT
                         SORT_PLUS(J,I,JP,IP,IC1)=COUNT
                        ENDIF
                        ENDDO
                        ENDDO
                  COUNT = COUNT + 1
                  ENDDO
                  ENDDO
                  ENDDO

           ELSE

          COUNT = 1                 !MULTI-LETTER CODE
          DO I = 1,N_LETTERS
          DO J = 1,N_LETTERS
          DO IP = 1,N_LETTERS
          DO JP = 1,N_LETTERS
          DO IC1  = 1,2
          SORT_PLUS(I,J,IP,JP,IC1) = COUNT
          COUNT = COUNT +1
          ENDDO
          ENDDO
          ENDDO
          ENDDO
          ENDDO

           ENDIF
           ENDIF  ! FOUR_PLUS

          IF (I_3RD_IS_CONTACT) THEN
            COUNT = 1
            DO I=1,N_LETTERS_CON
            DO J=I,N_LETTERS_CON
            SORT_CONTACT(I,J)=COUNT 
            SORT_CONTACT(J,I)=COUNT
            COUNT=COUNT+1
            ENDDO
            ENDDO
          ENDIF

        ENDIF
        
        IF (I_NON_ADD_CONTACT) THEN
          COUNT=1
          DO I=1,2
          DO J=1,2
          DO I_WELL=1,2
          DO I_TYPE=1,2
          DO I_WELL2=1,2
            SORT_NON_ADD(I,J,I_WELL,I_TYPE,I_WELL2)= COUNT
            COUNT=COUNT+1
          ENDDO
          ENDDO
          ENDDO
          ENDDO
          ENDDO
        ENDIF
 
       N_GAMMAS_PER_CON_WELL=N_LETTERS_CON*(N_LETTERS_CON+1)/2

        IF (N_LETTERS.EQ.2) THEN
          NGAMMA=48
        ELSEIF ((N_LETTERS.EQ.4).AND.I_3RD_IS_CONTACT) THEN
          IF (FOUR_PLUS) THEN
            IF (I_IGNORE_CHAIN_DIRN.AND.(.NOT.MISMATCH)) THEN
              NGAMMA=544+NUM_WELL*N_GAMMAS_PER_CON_WELL
            ELSEIF (I_IGNORE_CHAIN_DIRN.AND.MISMATCH) THEN
              NGAMMA=80+NUM_WELL*N_GAMMAS_PER_CON_WELL
            ELSE
              NGAMMA=1024+NUM_WELL*N_GAMMAS_PER_CON_WELL
            ENDIF
          ELSE
            IF (I_IGNORE_CHAIN_DIRN.AND.(.NOT.MISMATCH)) THEN
            NGAMMA=272+NUM_WELL*N_GAMMAS_PER_CON_WELL
            ELSE
            NGAMMA=512+NUM_WELL*N_GAMMAS_PER_CON_WELL
            ENDIF
          ENDIF
        ELSEIF ((N_LETTERS.EQ.4).AND.(.NOT.I_3RD_IS_CONTACT)) THEN
          NGAMMA=768
        ENDIF

        IF (I_NON_ADD_CONTACT) NGAMMA=NGAMMA+NGAMMA_NON_ADD
    
        OPEN(IGAMMA,FILE='GAMMA.DAT',ACTION='READ',STATUS='OLD',IOSTAT=OPEN_STATUS)
        IF (OPEN_STATUS.NE.0) THEN
          WRITE(SO,*) 'FAILURE TO OPEN GAMMA FILE'
          WRITE(SO,*) 'ERROR NUMBER ',OPEN_STATUS
          STOP
        ENDIF

        READ(IGAMMA,*) NGAMMA_CHECK
        IF (NGAMMA_CHECK.NE.NGAMMA) THEN
          WRITE(SO,*) 'GAMMA COCK-UP',NGAMMA_CHECK,NGAMMA
          STOP
        ENDIF

           WRITE (OARCHV,*)
           WRITE (OARCHV,*) 'GAMMAS USED'
           WRITE (OARCHV,*) '-----------'


        IF (I_3RD_IS_CONTACT.AND.(N_LETTERS.EQ.4)) THEN
           N_RANGES=2+NUM_WELL
           IF (FOUR_PLUS) THEN
             IF (I_IGNORE_CHAIN_DIRN.AND.(.NOT.MISMATCH)) THEN
               N_GAMMAS_IN_RANGE(1)=272
               N_GAMMAS_IN_RANGE(2)=272
             ELSEIF (I_IGNORE_CHAIN_DIRN.AND.MISMATCH) THEN
               N_GAMMAS_IN_RANGE(1)=40
               N_GAMMAS_IN_RANGE(2)=40
             ELSE
               N_GAMMAS_IN_RANGE(1)=512
               N_GAMMAS_IN_RANGE(2)=512
             ENDIF
           ELSE

             IF (I_IGNORE_CHAIN_DIRN) THEN
               N_GAMMAS_IN_RANGE(1)= 136
               N_GAMMAS_IN_RANGE(2)= 136
             ELSE
               N_GAMMAS_IN_RANGE(1)= 256
               N_GAMMAS_IN_RANGE(2)= 256
             ENDIF

           ENDIF

          DO I = 1,NUM_WELL
             N_GAMMAS_IN_RANGE(2+I)=N_GAMMAS_PER_CON_WELL
           ENDDO

        ELSE 
           N_RANGES=3
           DO I=1,N_RANGES
             N_GAMMAS_IN_RANGE(I)=4**N_LETTERS ! EG 1 TO 16 FOR 2-LETTER CODE
           ENDDO
        ENDIF


        DO 510 I_RANGES=1,N_RANGES

        DO 504 I504=1,N_GAMMAS_IN_RANGE(I_RANGES)  
           READ (IGAMMA,*,IOSTAT=READ_STATUS,ERR=99) OPTGAMMA(I504)
99         IF (READ_STATUS.NE.0) THEN
             WRITE(SO,*)  'FAILURE READING FROM GAMMA FILE (1ST READ)'
             WRITE(SO,*) 'ERROR NUMBER ',READ_STATUS
             WRITE(SO,*) OPTGAMMA(I504)
             STOP
           ENDIF
           WRITE (OARCHV,444) OPTGAMMA(I504)
444        FORMAT(F8.4)

504    CONTINUE

        DO 500 I500=1,20
           ICL=CLASS(I500)
           ICL_CON=CLASS_CON(I500)

           DO 501 I501=1,20
             JCL=CLASS(I501)
             JCL_CON=CLASS_CON(I501)
              GAMMAPT_CONTACT=SORT_CONTACT(ICL_CON,JCL_CON)

              DO 502 I502=1,20
                 IPCL=CLASS(I502)
                 DO 503 I503=1,20
                   JPCL=CLASS(I503)

              IF ( I_RANGES .LE. 2) THEN
                 IF (FOUR_PLUS) THEN
                   DO IC1 = 1,2
                   GAMMAPT=SORT_PLUS(ICL,JCL,IPCL,JPCL,IC1)
                   QCHRG2(I500,I501,I502,I503,IC1,I_RANGES)=OPTGAMMA(GAMMAPT)
                   ENDDO
                 ELSE ! FOUR_PLUS FALSE
                      IF (I_IGNORE_CHAIN_DIRN) THEN
                       GAMMAPT=SORT_SYM(ICL,JCL,IPCL,JPCL)
                       GAMMA_TEMP=OPTGAMMA(GAMMAPT)
                      ELSE
                       GAMMAPT=SORT(ICL,JCL,IPCL,JPCL)
                       GAMMA_TEMP=OPTGAMMA(GAMMAPT)
                      ENDIF
                 ENDIF ! FOUR_PLUS
               ENDIF ! I_RANGES 2
                                                                                                     
                   IF (I_3RD_IS_CONTACT.AND.(N_LETTERS.EQ.4))THEN
                    IF(I_RANGES.GT.2) THEN
                       GAMMA_TEMP=OPTGAMMA(GAMMAPT_CONTACT)
                    ENDIF
                   ENDIF
C                  THIRD DISTANCE CLASS GO CONTACTS
              IF ((.NOT.I_3RD_IS_CONTACT).AND.(N_LETTERS.EQ.4))THEN
                 IF(I_RANGES.GT.2) THEN
                     GAMMA_TEMP=OPTGAMMA(GAMMAPT)
                 ENDIF

              ENDIF
                                                                                  
C   IF ALTERNATIVE POTENTIAL QCHRG FOR THE CORRESPONDING WELL
C   SHOULD BE SET TO 0
                  IF((ALTPOTFLAG(1).GT.0).AND.(I_RANGES.EQ.3))THEN
                      QCHRG(I500,I501,I502,I503,I_RANGES)=0.0D0
                  ELSEIF((ALTPOTFLAG(2).GT.0).AND.(I_RANGES.EQ.4))THEN
                      QCHRG(I500,I501,I502,I503,I_RANGES)=0.0D0
C   THESE ARE ADDED BY GAREGIN PAPOIAN, 4/8/04
C   MODIFYING MEDIUM-RANGE AMH INTERACTIONS
                  ELSEIF((ALTPOTFLAG(2).GT.0).AND.(I_RANGES.EQ.2))THEN
                      QCHRG(I500,I501,I502,I503,I_RANGES)=GAMMA_TEMP*SRPLR(2)
C   MODIFYING SHORT-RANGE AMH INTERACTIONS
                  ELSEIF((ALTPOTFLAG(2).GT.0).AND.(I_RANGES.EQ.1))THEN
                      QCHRG(I500,I501,I502,I503,I_RANGES)=GAMMA_TEMP*SRPLR(1)
C   END OF GAP
                  ELSE
                   QCHRG(I500,I501,I502,I503,I_RANGES)=GAMMA_TEMP


                  ENDIF
                                                                                       
503              CONTINUE
502            CONTINUE

501      CONTINUE
500     CONTINUE

510     CONTINUE

C----------TSHEN ADD THE 21TH ONE, JUST RELY ON THE RELATIONS
C NEED TO GAMMA(21,"1 - 20","1 - 20","1 - 20","1 - MAX_WELL+2")  AND ITS OPPO
C AND GAMMA2 BUT WHAT IN GAMMA2??
C FIRST TRY MYCHCH=1.2 AND MYCHPO=1.1
      MYCHCH=1.2D0*1.2D0
      MYCHPO=1.2D0
                
      DO MYJ=1,20
         DO MYIP=1,20
            DO MYJP=1,20
               DO MYWEL=1, MAX_WELL+2
                  IF((CLASS(MYJ).EQ.4).OR.(CLASS(MYIP).EQ.4).OR.(CLASS(MYJP).EQ.4)) THEN
                     QCHRG(21,MYJ,MYIP,MYJP,MYWEL)=QCHRG(7,MYJ,MYIP,MYJP,MYWEL)
                     QCHRG(MYJ,21,MYIP,MYJP,MYWEL)=QCHRG(MYJ,7,MYIP,MYJP,MYWEL)
                  ELSEIF ((CLASS(MYJ).EQ.1).OR.(CLASS(MYIP).EQ.1).OR.(CLASS(MYJP).EQ.1)) THEN
                     QCHRG(21,MYJ,MYIP,MYJP,MYWEL)=MYCHPO*QCHRG(7,MYJ,MYIP,MYJP,MYWEL)
                     QCHRG(MYJ,21,MYIP,MYJP,MYWEL)=MYCHPO*QCHRG(MYJ,7,MYIP,MYJP,MYWEL)
                  ELSE
                     QCHRG(21,MYJ,MYIP,MYJP,MYWEL)=MYCHCH*QCHRG(7,MYJ,MYIP,MYJP,MYWEL)
                     QCHRG(MYJ,21,MYIP,MYJP,MYWEL)=MYCHCH*QCHRG(MYJ,7,MYIP,MYJP,MYWEL)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      DO MYJ=1,20
         DO MYIP=1,20
            DO MYJP=1,20
               DO MYK=1, 2
                  DO MYKP=1,2
                     IF((CLASS(MYJ).EQ.4).OR.(CLASS(MYIP).EQ.4).OR.(CLASS(MYJP).EQ.4)) THEN
                        QCHRG2(21,MYJ,MYIP,MYJP,MYK,MYKP)=QCHRG2(7,MYJ,MYIP,MYJP,MYK,MYKP)
                        QCHRG2(MYJ,21,MYIP,MYJP,MYK,MYKP)=QCHRG2(MYJ,7,MYIP,MYJP,MYK,MYKP)
                     ELSEIF ((CLASS(MYJ).EQ.1).OR.(CLASS(MYIP).EQ.1).OR.(CLASS(MYJP).EQ.1)) THEN
                        QCHRG2(21,MYJ,MYIP,MYJP,MYK,MYKP)=MYCHPO*QCHRG2(7,MYJ,MYIP,MYJP,MYK,MYKP)
                        QCHRG2(MYJ,21,MYIP,MYJP,MYK,MYKP)=MYCHPO*QCHRG2(MYJ,7,MYIP,MYJP,MYK,MYKP)
                     ELSE
                        QCHRG2(21,MYJ,MYIP,MYJP,MYK,MYKP)=MYCHCH*QCHRG2(7,MYJ,MYIP,MYJP,MYK,MYKP)
                        QCHRG2(MYJ,21,MYIP,MYJP,MYK,MYKP)=MYCHCH*QCHRG2(MYJ,7,MYIP,MYJP,MYK,MYKP)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO MYIP=1,20
            DO MYJP=1,20
               DO MYWEL=1, MAX_WELL+2
                  IF((CLASS(MYIP).EQ.4).OR.(CLASS(MYJP).EQ.4)) THEN
                     QCHRG(21,21,MYIP,MYJP,MYWEL)=QCHRG(7,21,MYIP,MYJP,MYWEL)
                     QCHRG(21,21,MYIP,MYJP,MYWEL)=QCHRG(21,7,MYIP,MYJP,MYWEL)
                  ELSEIF ((CLASS(MYIP).EQ.1).OR.(CLASS(MYJP).EQ.1)) THEN
                     QCHRG(21,21,MYIP,MYJP,MYWEL)=MYCHPO*QCHRG(7,21,MYIP,MYJP,MYWEL)
                     QCHRG(21,21,MYIP,MYJP,MYWEL)=MYCHPO*QCHRG(21,7,MYIP,MYJP,MYWEL)
                  ELSE
                     QCHRG(21,21,MYIP,MYJP,MYWEL)=MYCHCH*QCHRG(7,21,MYIP,MYJP,MYWEL)
                     QCHRG(21,21,MYIP,MYJP,MYWEL)=MYCHCH*QCHRG(21,7,MYIP,MYJP,MYWEL)
                  ENDIF
               ENDDO

               DO MYK=1, 2
                  DO MYKP=1,2
                     IF((CLASS(MYIP).EQ.4).OR.(CLASS(MYJP).EQ.4)) THEN
                        QCHRG2(21,21,MYIP,MYJP,MYK,MYKP)=QCHRG2(7,21,MYIP,MYJP,MYK,MYKP)
                        QCHRG2(21,21,MYIP,MYJP,MYK,MYKP)=QCHRG2(21,7,MYIP,MYJP,MYK,MYKP)
                     ELSEIF ((CLASS(MYIP).EQ.1).OR.(CLASS(MYJP).EQ.1)) THEN
                        QCHRG2(21,21,MYIP,MYJP,MYK,MYKP)=MYCHPO*QCHRG2(7,21,MYIP,MYJP,MYK,MYKP)
                        QCHRG2(21,21,MYIP,MYJP,MYK,MYKP)=MYCHPO*QCHRG2(21,7,MYIP,MYJP,MYK,MYKP)
                     ELSE
                        QCHRG2(21,21,MYIP,MYJP,MYK,MYKP)=MYCHCH*QCHRG2(7,21,MYIP,MYJP,MYK,MYKP)
                        QCHRG2(21,21,MYIP,MYJP,MYK,MYKP)=MYCHCH*QCHRG2(21,7,MYIP,MYJP,MYK,MYKP)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

C-------END OF TSHEN MODI-----------

        IF (I_NON_ADD_CONTACT) THEN
          DO I=1,NGAMMA_NON_ADD
            READ (IGAMMA,*,IOSTAT=READ_STATUS,ERR=199) GAMMA_NON_ADD(I)
199            IF (READ_STATUS.NE.0) THEN
              WRITE(SO,*) 'IAILURE READING FROM GAMMA FILE (2ND READ)'
              WRITE(SO,*)'ERROR NUMBER ',READ_STATUS
              STOP
            ENDIF
          ENDDO
        ENDIF

        CLOSE (IGAMMA)

        OPEN(IHBOND,FILE='PARAMS/ANTI_HB',STATUS='OLD',IOSTAT=OPEN_STATUS)
        IF (OPEN_STATUS.NE.0) THEN
          WRITE(SO,*)'FAILURE TO OPEN ANTI_HB FILE'
          WRITE(SO,*) 'ERROR NUMBER ',OPEN_STATUS
          STOP
        ENDIF

        DO ISIT1 = 1,20
        READ(IHBOND,77)(ANTI_HB(ISIT1,ISIT2,1),ISIT2=1,20)
        ENDDO
        READ(IHBOND,*)
        DO ISIT1 = 1,20
        READ(IHBOND,77)(ANTI_HB(ISIT1,ISIT2,2),ISIT2=1,20)
        ENDDO
        CLOSE(IHBOND)

        OPEN(IHBOND,FILE='PARAMS/ANTI_NHB',STATUS='OLD',IOSTAT=OPEN_STATUS)
        IF (OPEN_STATUS.NE.0) THEN
          WRITE(SO,*)'FAILURE TO OPEN ANTI_NHB FILE'
          WRITE(SO,*) 'ERROR NUMBER ',OPEN_STATUS
          STOP
        ENDIF

        DO ISIT1 = 1,20
        READ(IHBOND,77)(ANTI_NHB(ISIT1,ISIT2,1),ISIT2=1,20)
        ENDDO
        READ(IHBOND,*)
        DO ISIT1 = 1,20
        READ(IHBOND,77)(ANTI_NHB(ISIT1,ISIT2,2),ISIT2=1,20)
        ENDDO
        CLOSE(IHBOND)

        OPEN(IHBOND,FILE='PARAMS/PARA_ONE',STATUS='OLD')
        DO ISIT1 = 1,20
        READ(IHBOND,*)PARA_ONE(ISIT1)
        ENDDO
        CLOSE(IHBOND)

        OPEN(IHBOND,FILE='PARAMS/ANTI_ONE',STATUS='OLD')
        DO ISIT1 = 1,20
        READ(IHBOND,*)ANTI_ONE(ISIT1)
        ENDDO
        CLOSE(IHBOND)

        OPEN(IHBOND,FILE='PARAMS/PARA_HB',
     *  STATUS='OLD',IOSTAT=OPEN_STATUS)
        IF (OPEN_STATUS.NE.0) THEN
          WRITE(SO,*) 'FAILURE TO OPEN PARA_HB FILE'
          WRITE(SO,*) 'ERROR NUMBER ',OPEN_STATUS
          STOP
        ENDIF

        DO ISIT1 = 1,20
        READ(IHBOND,77)(PARA_HB(ISIT1,ISIT2,1),ISIT2=1,20)
        ENDDO
        READ(IHBOND,*)
        DO ISIT1 = 1,20
        READ(IHBOND,77)(PARA_HB(ISIT1,ISIT2,2),ISIT2=1,20)
        ENDDO
        CLOSE(IHBOND)

77        FORMAT(20(F8.5,1X))

        RETURN
      END
