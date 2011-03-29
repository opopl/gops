
C     --------------------- RNDCOL ----------------------

      SUBROUTINE RNDCOL(NMRES,MAXPRO,NUMPRO,
     *                  IRES,PROCNT,NUMCRD,PRCORD,ISEED_AMH)

C     ---------------------------------------------------

C     RNDCOL GENERATE RANDOM COIL CONFIGURATIONS

C     ARGUMENTS:

C        AMHMAXSIZ- MAXIMUM TARGET PROTEIN SIZE (I)
C        MAXPRO- MAXIMUM NUMBER OF TRIAL STRUCTURES (I)
C        NUMPRO- ACTUAL NUMBER OF TRIAL STRUCTURES (I)
C        PROCNT- NUMBER OF RANDOM COIL CONFIGURATIONS TO
C                BE CONSTRUCTED (I)
C        NUMCRD- NUMBER OF COORDINATE TYPES (I)
C        PRCORD- RANDOM COIL STRUCURES (O)
C        BONDLN- BOND LENGTHS FOR COORDINATES TYPES (I)

C     ---------------------------------------------------

      USE AMHGLOBALS,  ONLY:AMHMAXSIZ,PDB,QUENCH_CRD,X_MCP,AMINOA

      IMPLICIT NONE

C     PASSED ARGUMENT DECLARATIONS:

      INTEGER NMRES,NUMPRO,PROCNT,NUMCRD,MAXPRO,IRES(AMHMAXSIZ),GLY_C,ISEED_AMH(4)

      DOUBLE PRECISION PRCORD(AMHMAXSIZ,3,MAXPRO,NUMCRD)

      LOGICAL TOUCH

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     INTERNAL VARIABLES
       DOUBLE PRECISION PI,PI2,DIST1,DIST2,
     *     NITCORD(AMHMAXSIZ,3),CPRCORD(AMHMAXSIZ,3),PHI,CPHI,SPHI,CPSI,SPSI,
     *     XPRIME(3),YPRIME(3),ZPRIME(3),XLEN,YLEN,ZLEN

         INTEGER TRY(AMHMAXSIZ,2), III, JJJ

         CHARACTER*3 RES_TYPE

C        --- DO LOOP INDICES ---

         INTEGER I501,I503,I506,I507,I1

C MY CHANGES BELOW   - SKIP
       DOUBLE PRECISION RAN_NUM_ARRAY(1),RAN_NUM

       EXTERNAL SLARNV

C     --------------------- BEGIN -----------------------
        PI=3.141592653589793238462643383279502884197
        PI2=2.0D0*PI
        DO 600 I1=1,AMHMAXSIZ
          TRY(I1,1)=0
          TRY(I1,2)=0
600        CONTINUE


        PROCNT=PROCNT + 1
        DO 503 I503=PROCNT,NUMPRO

C        GENERATE INITIAL CONFIGURATIONS FOR RESIDUE 1

        PRCORD(1,1,I503,1)=0.0D0
        PRCORD(1,2,I503,1)=0.0D0
        PRCORD(1,3,I503,1)=0.0D0
        CPRCORD(1,1)=0.0D0
        CPRCORD(1,2)=0.0D0
        CPRCORD(1,3)=1.53D0
        NITCORD(1,1)=0.0D0
        NITCORD(1,2)=-1.3859293D0
        NITCORD(1,3)=-0.4899999D0
        PRCORD(1,1,I503,2)=-1.2574047D0
        PRCORD(1,2,I503,2)=0.7259629D0
        PRCORD(1,3,I503,2)=-0.5133333D0
        
        IF (IRES(1).EQ.8) THEN
          PRCORD(1,1,I503,2)=PRCORD(1,1,I503,1)
          PRCORD(1,2,I503,2)=PRCORD(1,2,I503,1)
          PRCORD(1,3,I503,2)=PRCORD(1,3,I503,1)
        ENDIF

        I501=1

         DO WHILE (I501.LT.NMRES)

C           GENERATE RANDOM NUMBERS FOR X- AND
C           Y-COORDINATES

105      ZPRIME(1)=(CPRCORD(I501,1)-PRCORD(I501,1,I503,1))
         ZPRIME(2)=(CPRCORD(I501,2)-PRCORD(I501,2,I503,1))
         ZPRIME(3)=(CPRCORD(I501,3)-PRCORD(I501,3,I503,1))

         ZLEN=DSQRT(ZPRIME(1)*ZPRIME(1)+ZPRIME(2)*ZPRIME(2)+ZPRIME(3)*ZPRIME(3))

         ZPRIME(1)=ZPRIME(1)/ZLEN
         ZPRIME(2)=ZPRIME(2)/ZLEN
         ZPRIME(3)=ZPRIME(3)/ZLEN

         YPRIME(1)=(PRCORD(I501,1,I503,1)-NITCORD(I501,1)-ZPRIME(1)*0.4899999D0)
         YPRIME(2)=(PRCORD(I501,2,I503,1)-NITCORD(I501,2)-ZPRIME(2)*0.4899999D0)
         YPRIME(3)=(PRCORD(I501,3,I503,1)-NITCORD(I501,3)-ZPRIME(3)*0.4899999D0)

         YLEN=DSQRT(YPRIME(1)*YPRIME(1)+YPRIME(2)*YPRIME(2)+YPRIME(3)*YPRIME(3))

         YPRIME(1)=YPRIME(1)/YLEN
         YPRIME(2)=YPRIME(2)/YLEN
         YPRIME(3)=YPRIME(3)/YLEN

         XPRIME(1)=YPRIME(2)*ZPRIME(3)-YPRIME(3)*ZPRIME(2)
         XPRIME(2)=YPRIME(3)*ZPRIME(1)-YPRIME(1)*ZPRIME(3)
         XPRIME(3)=YPRIME(1)*ZPRIME(2)-YPRIME(2)*ZPRIME(1)

         XLEN=DSQRT(XPRIME(1)*XPRIME(1)+XPRIME(2)*XPRIME(2)+XPRIME(3)*XPRIME(3))

         XPRIME(1)=XPRIME(1)/XLEN
         XPRIME(2)=XPRIME(2)/XLEN
         XPRIME(3)=XPRIME(3)/XLEN

C 110               PHI=RAND()*1.92+1.57

         CALL SLARNV(1, ISEED_AMH,1,RAN_NUM_ARRAY)
               RAN_NUM=RAN_NUM_ARRAY(1)

110            PHI=RAN_NUM*1.92D0+1.57D0

               IF (PHI.GT.PI) PHI=PHI+2.09D0
               CPHI=DCOS(PHI)
               SPHI=DSIN(PHI)

               PRCORD(I501,1,I503,3)
     *                  =-1.0515796D0*SPHI*XPRIME(1)
     *                   +1.0515796D0*CPHI*YPRIME(1)
     *                   +0.6570998D0*ZPRIME(1)
     *                   +CPRCORD(I501,1)
               PRCORD(I501,2,I503,3)
     *                  =-1.0515796D0*SPHI*XPRIME(2)
     *                   +1.0515796D0*CPHI*YPRIME(2)
     *                   +0.6570998D0*ZPRIME(2)
     *                   +CPRCORD(I501,2)
               PRCORD(I501,3,I503,3)
     *                  =-1.0515796D0*SPHI*XPRIME(3)
     *                   +1.0515796D0*CPHI*YPRIME(3)
     *                   +0.6570998D0*ZPRIME(3)
     *                   +CPRCORD(I501,3)

               NITCORD(I501+1,1)
     *                   =1.20588D0*SPHI*XPRIME(1)
     *                   -1.20588D0*CPHI*YPRIME(1)
     *                   +0.5368923D0*ZPRIME(1)
     *                   +CPRCORD(I501,1)
               NITCORD(I501+1,2)
     *                   =1.20588D0*SPHI*XPRIME(2)
     *                   -1.20588D0*CPHI*YPRIME(2)
     *                   +0.5368923D0*ZPRIME(2)
     *                   +CPRCORD(I501,2)
               NITCORD(I501+1,3)
     *                   =1.20588D0*SPHI*XPRIME(3)
     *                   -1.20588D0*CPHI*YPRIME(3)
     *                   +0.5368923D0*ZPRIME(3)
     *                   +CPRCORD(I501,3)

               PRCORD(I501+1,1,I503,1)
     *                   =1.4358387D0*SPHI*XPRIME(1)
     *                   -1.4358387D0*CPHI*YPRIME(1)
     *                   +1.9887942D0*ZPRIME(1)
     *                   +CPRCORD(I501,1)
               PRCORD(I501+1,2,I503,1)
     *                   =1.4358387D0*SPHI*XPRIME(2)
     *                   -1.4358387D0*CPHI*YPRIME(2)
     *                   +1.9887942D0*ZPRIME(2)
     *                   +CPRCORD(I501,2)
               PRCORD(I501+1,3,I503,1)
     *                   =1.4358387D0*SPHI*XPRIME(3)
     *                   -1.4358387D0*CPHI*YPRIME(3)
     *                   +1.9887942D0*ZPRIME(3)
     *                   +CPRCORD(I501,3)

C      MAKE CHAIN SELF-AVOIDING

            TOUCH=.FALSE.

            DO 506 I506=1,I501-1

             DIST1=DSQRT( (PRCORD(I501+1,1,I503,1)-PRCORD(I506,1,I503,1))**2
     *                   +(PRCORD(I501+1,2,I503,1)-PRCORD(I506,2,I503,1))**2
     *                   +(PRCORD(I501+1,3,I503,1)-PRCORD(I506,3,I503,1))**2)

             DIST2=DSQRT( (PRCORD(I501+1,1,I503,1)-PRCORD(I506,1,I503,2))**2
     *                   +(PRCORD(I501+1,2,I503,1)-PRCORD(I506,2,I503,2))**2
     *                   +(PRCORD(I501+1,3,I503,1)-PRCORD(I506,3,I503,2))**2)

                IF ( (DIST1.LT.5.0D0).OR.(DIST2.LT.4.0D0) ) 
     *             TOUCH=.TRUE.

506              CONTINUE

                IF (TOUCH) THEN

                   TRY(I501,1)=TRY(I501,1)+1
                   IF (TRY(I501,1).LT.11) GOTO 110
                   
C           FAILED 10 TIMES-BACKTRACK TO PREVIOUS GOOD BOND

C                   WRITE (6,*) I501,1,TRY(I501,1)

                   DO 301 I501=I501-1,1,-1

                   TRY(I501+1,1)=0
                   TRY(I501,2)=TRY(I501,2)+1
                   IF (TRY(I501,2).LT.11) GOTO 115 

                   TRY(I501,2)=0
                   TRY(I501,1)=TRY(I501,1)+1
                   IF (TRY(I501,1).LT.11) GOTO 105
C                   WRITE (6,*) I501,1,TRY(I501,1)

301                CONTINUE        

                   STOP ' FAILED IN RNDCOL '

                ENDIF

115    ZPRIME(1)=(PRCORD(I501+1,1,I503,1)-NITCORD(I501+1,1))/1.47D0
       ZPRIME(2)=(PRCORD(I501+1,2,I503,1)-NITCORD(I501+1,2))/1.47D0
       ZPRIME(3)=(PRCORD(I501+1,3,I503,1)-NITCORD(I501+1,3))/1.47D0

       ZLEN=DSQRT(ZPRIME(1)*ZPRIME(1)+ZPRIME(2)*ZPRIME(2)+ZPRIME(3)*ZPRIME(3))

               ZPRIME(1)=ZPRIME(1)/ZLEN
               ZPRIME(2)=ZPRIME(2)/ZLEN
               ZPRIME(3)=ZPRIME(3)/ZLEN

               YPRIME(1)=(NITCORD(I501+1,1)-CPRCORD(I501,1)
     *                     -ZPRIME(1)*0.7189235D0)
     *                     /1.1070451D0
               YPRIME(2)=(NITCORD(I501+1,2)-CPRCORD(I501,2)
     *                     -ZPRIME(2)*0.7189235D0)
     *                     /1.1070451D0
               YPRIME(3)=(NITCORD(I501+1,3)-CPRCORD(I501,3)
     *                     -ZPRIME(3)*0.7189235D0)
     *                     /1.1070451D0

        YLEN=SQRT(YPRIME(1)*YPRIME(1)+YPRIME(2)*YPRIME(2)+YPRIME(3)*YPRIME(3))

        YPRIME(1)=YPRIME(1)/YLEN
        YPRIME(2)=YPRIME(2)/YLEN
        YPRIME(3)=YPRIME(3)/YLEN

        XPRIME(1)=YPRIME(2)*ZPRIME(3)-YPRIME(3)*ZPRIME(2)
        XPRIME(2)=YPRIME(3)*ZPRIME(1)-YPRIME(1)*ZPRIME(3)
        XPRIME(3)=YPRIME(1)*ZPRIME(2)-YPRIME(2)*ZPRIME(1)

        XLEN=SQRT(XPRIME(1)*XPRIME(1)+XPRIME(2)*XPRIME(2)+XPRIME(3)*XPRIME(3))

        XPRIME(1)=XPRIME(1)/XLEN
        XPRIME(2)=XPRIME(2)/XLEN
        XPRIME(3)=XPRIME(3)/XLEN

C  120               PHI=RAND()*1.75+3.49

               CALL SLARNV(1, ISEED_AMH,1,RAN_NUM_ARRAY(1))
               RAN_NUM=RAN_NUM_ARRAY(1)

120            PHI=RAN_NUM*1.75D0+3.49D0
               CPHI=DCOS(PHI)
               SPHI=DSIN(PHI)
               CPSI=DCOS(PHI-PI2/3.0D0)
               SPSI=DSIN(PHI-PI2/3.0D0)

               PRCORD(I501+1,1,I503,2)
     *                  =1.4519259D0*SPSI*XPRIME(1)-1.4519259D0*CPSI*YPRIME(1)
     *                  +0.5133333D0*ZPRIME(1)+PRCORD(I501+1,1,I503,1)
               PRCORD(I501+1,2,I503,2)
     *                  =1.4519259D0*SPSI*XPRIME(2)-1.4519259D0*CPSI*YPRIME(2)
     *                  +0.5133333D0*ZPRIME(2)+PRCORD(I501+1,2,I503,1)
               PRCORD(I501+1,3,I503,2)
     *                  =1.4519259D0*SPSI*XPRIME(3)-1.4519259D0*CPSI*YPRIME(3)
     *                  +0.5133333D0*ZPRIME(3)+PRCORD(I501+1,3,I503,1)

        IF (IRES(I501+1).EQ.8) THEN
          PRCORD(I501+1,1,I503,2)=PRCORD(I501+1,1,I503,1)
          PRCORD(I501+1,2,I503,2)=PRCORD(I501+1,2,I503,1)
          PRCORD(I501+1,3,I503,2)=PRCORD(I501+1,3,I503,1)
        ENDIF

C      MAKE CHAIN SELF-AVOIDING

            TOUCH=.FALSE.

            DO 507 I507=1,I501-1

                DIST1=DSQRT( (PRCORD(I501+1,1,I503,2)
     *                      -PRCORD(I507,1,I503,1))**2
     *                     +(PRCORD(I501+1,2,I503,2)
     *                      -PRCORD(I507,2,I503,1))**2
     *                     +(PRCORD(I501+1,3,I503,2)
     *                      -PRCORD(I507,3,I503,1))**2)

                DIST2=DSQRT( (PRCORD(I501+1,1,I503,2)
     *                      -PRCORD(I507,1,I503,2))**2
     *                     +(PRCORD(I501+1,2,I503,2)
     *                      -PRCORD(I507,2,I503,2))**2
     *                     +(PRCORD(I501+1,3,I503,2)
     *                      -PRCORD(I507,3,I503,2))**2)

                IF ( (DIST1.LT.4.0D0).OR.(DIST2.LT.3.0D0) ) 
     *             TOUCH=.TRUE.

507        CONTINUE

                IF (TOUCH) THEN

                   TRY(I501,2)=TRY(I501,2)+1
                   IF (TRY(I501,2).LT.11) GOTO 120

C           FAILED 10 TIMES-BACKTRACK TO PREVIOUS GOOD BOND

C                   WRITE (6,*) I501,2,TRY(I501,2)
                   TRY(I501,2)=0
                   TRY(I501,1)=TRY(I501,1)+1
                   IF (TRY(I501,1).LT.11) GOTO 105
C                   WRITE (6,*) I501,1,TRY(I501,1)

                   DO 302 I501=I501-1,1,-1

                   TRY(I501+1,1)=0
                   TRY(I501,2)=TRY(I501,2)+1
                   IF (TRY(I501,2).LT.11) GOTO 115
C                   WRITE (6,*) I501,2,TRY(I501,2)

                   TRY(I501,2)=0
                   TRY(I501,1)=TRY(I501,1)+1
                   IF (TRY(I501,1).LT.11) GOTO 105
                   WRITE (6,*) I501,1,TRY(I501,1)

302             CONTINUE

                   STOP ' FAILED IN RNDCOL '

                ENDIF

       CPRCORD(I501+1,1)=1.4424978D0*SPHI*XPRIME(1)-1.4424978D0*CPHI*YPRIME(1)
     *                  +0.51D0*ZPRIME(1)+PRCORD(I501+1,1,I503,1)
       CPRCORD(I501+1,2)=1.4424978D0*SPHI*XPRIME(2)-1.4424978D0*CPHI*YPRIME(2)
     *                  +0.51D0*ZPRIME(2)+PRCORD(I501+1,2,I503,1)
       CPRCORD(I501+1,3)=1.4424978D0*SPHI*XPRIME(3)-1.4424978D0*CPHI*YPRIME(3)
     *                  +0.51D0*ZPRIME(3)+PRCORD(I501+1,3,I503,1)

        I501=I501+1
       END DO

 503   CONTINUE

         OPEN(UNIT=PDB,FILE='0.PDB',FORM='FORMATTED',STATUS = 'UNKNOWN')

       DO 664 III = 1,NMRES
            RES_TYPE = 'ALA'
             IF (IRES(III) .EQ. 8) THEN
              RES_TYPE = 'GLY'
             END IF

CACACACACACACACACACACAC
       WRITE(PDB,665) III, RES_TYPE,III,(PRCORD(III, JJJ, 1, 1), JJJ =1,3), III
665    FORMAT('ATOM    ',I3,'  CA  ',A3,'   ',I3,'    ',3(F8.3),
     *             '  1.00  0.00      TPDB ',I3)
C'C'C'C'C'C'C'C'C'C'C'C'
       WRITE(PDB,671) III, RES_TYPE,III, (CPRCORD(III, JJJ),JJJ =1,3), III
671   FORMAT('ATOM    ',I3,'  C   ',A3,'   ',I3,'    ',3(F8.3),'  1.00  0.00      TPDB ',I3)
COOOOOOOOOOOOOOOOOOOOOOOO
       IF (III .LT. NMRES) THEN
        WRITE(PDB,670) III, RES_TYPE,III, (PRCORD(III, JJJ, 1,3),JJJ =1,3), III
670     FORMAT('ATOM    ',I3,'  O   ',A3,'   ',I3,'    ',3(F8.3),'  1.00  0.00      TPDB ',I3)
       END IF
CNNNNNNNNNNNNNNNNNNNNNNN
        WRITE(PDB,669) III, RES_TYPE, III, (NITCORD(III, JJJ),JJJ =1,3), III
669     FORMAT('ATOM    ',I3,'  N   ',A3,'   ',I3,'    ',3(F8.3),'  1.00  0.00      TPDB ',I3)

CBCBCBCBCBCBCBCBCBCBCBCB
       IF (IRES(III) .NE. 8) THEN
        WRITE(PDB,668) III, III, (PRCORD(III, JJJ, 1, 2),JJJ =1,3), III
668     FORMAT('ATOM    ',I3,'  CB  ALA   ',I3,'    ',3(F8.3),'  1.00  0.00      TPDB ',I3)
       END IF
CCCCCCCCCCCCCCCCCCCCCCCCC

664     CONTINUE
        CLOSE(PDB)

C  MCP STRUCTURE FOR GMIN

       GLY_C = 0

       DO 764 III = 1,NMRES

       IF (IRES(III).EQ.8) GLY_C = GLY_C +1

C  PRCORD(III, JJJ, 1, 1)  = PRCORD(RES ID , COORD ID, NUM PRO ID , ATOM ID)
C  CONVERT NMRES TO NUM OF ATOMS 
C  ASSUME X(CA(X,Y,Z),CB(X,Y,Z),O(X,Y,Z))  

       RES_TYPE = AMINOA(IRES(III))

       IF (IRES(III).EQ.8) THEN
        X_MCP(9*(III-1)+1-(GLY_C-1)*3) = (PRCORD(III, 1, 1, 1))   !  CA X
        X_MCP(9*(III-1)+2-(GLY_C-1)*3) = (PRCORD(III, 2, 1, 1))   !  CA Y
        X_MCP(9*(III-1)+3-(GLY_C-1)*3) = (PRCORD(III, 3, 1, 1))   !  CA Z
!       X_MCP(9*(III-1)+4) = (PRCORD(III, 1, 1, 2))   !  CB X
!       X_MCP(9*(III-1)+5) = (PRCORD(III, 2, 1, 2))   !  CB Y
!       X_MCP(9*(III-1)+6) = (PRCORD(III, 3, 1, 2))   !  CB Z
        X_MCP(9*(III-1)+4-(GLY_C-1)*3) = (PRCORD(III, 1, 1, 3))   !   O X
        X_MCP(9*(III-1)+5-(GLY_C-1)*3) = (PRCORD(III, 2, 1, 3))   !   O Y
        X_MCP(9*(III-1)+6-(GLY_C-1)*3) = (PRCORD(III, 3, 1, 3))   !   O Z
      ELSE
        X_MCP(9*(III-1)+1-GLY_C*3) = (PRCORD(III, 1, 1, 1))   !  CA X
        X_MCP(9*(III-1)+2-GLY_C*3) = (PRCORD(III, 2, 1, 1))   !  CA Y
        X_MCP(9*(III-1)+3-GLY_C*3) = (PRCORD(III, 3, 1, 1))   !  CA Z
        X_MCP(9*(III-1)+4-GLY_C*3) = (PRCORD(III, 1, 1, 2))   !  CB X
        X_MCP(9*(III-1)+5-GLY_C*3) = (PRCORD(III, 2, 1, 2))   !  CB Y
        X_MCP(9*(III-1)+6-GLY_C*3) = (PRCORD(III, 3, 1, 2))   !  CB Z
        X_MCP(9*(III-1)+7-GLY_C*3) = (PRCORD(III, 1, 1, 3))   !   O X
        X_MCP(9*(III-1)+8-GLY_C*3) = (PRCORD(III, 2, 1, 3))   !   O Y
        X_MCP(9*(III-1)+9-GLY_C*3) = (PRCORD(III, 3, 1, 3))   !   O Z
      ENDIF

CACACACACACACACACACACAC
      IF (IRES(III) .NE. 8) THEN
C        WRITE(6,765) III, RES_TYPE,III,X_MCP(9*(III-1)+1-GLY_C*3),X_MCP(9*(III-1)+2-GLY_C*3),X_MCP(9*(III-1)+3-GLY_C*3), III
765     FORMAT('ATOM    ',I3,'  CA  ',A3,'   ',I3,'    ',3(F8.3),'  1.00  0.00      TPDB ',I3)
      ENDIF

      IF (IRES(III) .EQ. 8) THEN
C      WRITE(6,865)III,RES_TYPE,III,X_MCP(9*(III-1)+1-(GLY_C-1)*3),X_MCP(9*(III-1)+2-(GLY_C-1)*3),X_MCP(9*(III-1)+3-(GLY_C-1)*3),III
865     FORMAT('ATOM    ',I3,'  XA  ',A3,'   ',I3,'    ',3(F8.3),'  1.00  0.00     TPDB ',I3)
      ENDIF

CBCBCBCBCBCBCBCBCBCBCBCB
       IF (IRES(III) .NE. 8) THEN
C        WRITE(6,768) III, RES_TYPE,III, X_MCP(9*(III-1)+4- GLY_C*3),X_MCP(9*(III-1)+5-GLY_C*3),X_MCP(9*(III-1)+6-GLY_C*3) , III
768     FORMAT('ATOM    ',I3,'  CB  ',A3,'   ',I3,'    ',3(F8.3),'  1.00  0.00      TPDB ',I3)
       END IF

COOOOOOOOOOOOOOOOOOOOOOOO
C       IF (III .LT. NMRES) THEN
      IF (IRES(III) .NE. 8) THEN
C        WRITE(6,770) III, RES_TYPE,III, X_MCP(9*(III-1)+7-GLY_C*3),X_MCP(9*(III-1)+8-GLY_C*3),X_MCP(9*(III-1)+9-GLY_C*3), III
770     FORMAT('ATOM    ',I3,'  O   ',A3,'   ',I3,'    ',3(F8.3),'  1.00  0.00      TPDB ',I3)
       ENDIF

      IF (IRES(III) .EQ. 8) THEN
C       WRITE(6,965)III,RES_TYPE,III,X_MCP(9*(III-1)+4-(GLY_C-1)*3),X_MCP(9*(III-1)+5-(GLY_C-1)*3),X_MCP(9*(III-1)+6-(GLY_C-1)*3),III
965     FORMAT('ATOM    ',I3,'  X0  ',A3,'   ',I3,'    ',3(F8.3),'  1.00  0.00      TPDB ',I3)
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCC

764     CONTINUE

C      WRITE(6,*)'OUT RNDCOL'

C       DO 364 III = 1,NMRES*3*3
C           WRITE(6,*)'COORD   RES',  X_MCP(III), III
C364    CONTINUE 

   
      QUENCH_CRD(:,:,:,:,1)=PRCORD

      RETURN
      END
