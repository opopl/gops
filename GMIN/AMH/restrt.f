
C     --------------------- INTSTR ----------------------

      SUBROUTINE RESTRT(PROCNT,MAXSIZ,NMRES,MAXPRO,
     *                  NUMPRO,MAXCRD,NUMCRD,PRCORD,
     *                  OARCHV)

C     ---------------------------------------------------

C     RESTRT READS IN PREVIOUSLY CALCULATED STRUCTURES

C     ---------------------------------------------------

      USE AMHGLOBALS,  ONLY:SO, QUENCH,NQUENCH,QUENCH_CRD,MAXMOV,TEMTUR,
     *               TEMTUR_QUENCH,ITGRD,X_MCP,IRES

      IMPLICIT NONE

C     ARGUMENT DECLARATIONS:

       INTEGER PROCNT,MAXSIZ,NMRES,MAXPRO,MAXCRD,NUMCRD,OARCHV,NUMPRO,I602
       CHARACTER(LEN=20) :: OTEMP 
       CHARACTER(LEN=20) :: OSTRING
       CHARACTER(LEN=2) :: SDUMMY
       

       DOUBLE PRECISION PRCORD(MAXSIZ,3,MAXPRO,MAXCRD),X,Y,Z

C     INTERNAL VARIABLES:

       INTEGER NUMPRR,NMRSS,NUMCRR,IDUMMY,GLY_C

C        --- DO LOOP INDICES ---

       INTEGER I500,I502,I503
 
C        --- IMPLIED DO LOOP INDICES ---

         INTEGER I2,III

C        ----  NEW ANNSCH STRUFF
         INTEGER GRID_NUM,INC,I501,I504

         DOUBLE PRECISION RINC,Q_TEMP(MAXMOV)

         RINC = 0.0D0
         NUMPRR = 1 

C     --------------------- BEGIN -----------------------

C     READ IN PREVIOUS PROTEINS 
C     Q_TEMP RECORDS QUENCHING TEMPERATURE
C     OPEN FILE CONTAINING COORDINATES

C     THE STRUCTURE FILE USED TO BE T.DAT 
C     I'VE CHANGED SO THAT FORMAT IS COMPATIBLE
C     WITH A SEGMENT OF A MOVIE FILE

      OPEN(UNIT=80,FILE='MOVIESEG',STATUS='OLD',FORM='FORMATTED')

C     READ IN PREVIOUS STRUCTURES

C      READ(80,67)NMRSS,NUMCRR,NUMPRR,NQUENCH
      WRITE(SO,*) 'QUENCH, NQUENCH=',QUENCH,NQUENCH
CC   67 FORMAT(3(I3,1X))
   67 FORMAT(4(I8,1X))

       IF (NUMPRR.GT.NUMPRO) THEN
          WRITE(SO,*) 'TOO MANY STRUCTURES IN MOVIESEG FILE'
          WRITE(SO,*) NUMPRR,NUMPRO
          STOP
       ENDIF

       IF (NQUENCH.GT.MAXMOV .AND. QUENCH) THEN
          WRITE(SO,*) 'TOO MANY STRUCTURES TO QUENCH',MAXMOV,QUENCH
          STOP
       ENDIF
  
        IF (.NOT. QUENCH) NQUENCH = 1

C      IF( NMRSS.NE.NMRES )THEN
C         WRITE(OARCHV,455)NMRSS,NMRES
C455      FORMAT(/'RESTRT: RESTART -- # RESIDUES NOT CONSISTENT ',2(I4,1X))
C         STOP
C      ENDIF

C     READ IN TRIAL STRUCTURES FOR NON-QUENCH RUN

      IF (.NOT. QUENCH) THEN

      DO 502 I502=1,NUMPRR

          READ(80,*)
          READ(80,*)
            DO 500 I500=1,NMRES
                  READ(80,76)
     *            (PRCORD(I500,I2,I502,1),I2=1,3),
     *            (PRCORD(I500,I2,I502,2),I2=1,3),
     *            (PRCORD(I500,I2,I502,3),I2=1,3)

                  WRITE(6,76)
     *            (PRCORD(I500,I2,I502,1),I2=1,3),
     *            (PRCORD(I500,I2,I502,2),I2=1,3),
     *            (PRCORD(I500,I2,I502,3),I2=1,3)
76     FORMAT(4X,3(F8.3,1X),4X,3(F8.3,1X),4X, 3(F8.3,1X))    
500       CONTINUE
C             I602=1

C            DO 500 I500=1,NMRES
C                  WRITE(6,*) 'I500',I500,NMRES
C                  READ(80,*) SDUMMY,X,Y,Z
C                  WRITE(6,*) SDUMMY,X,Y,Z
CC                  PRCORD(I500,1,I602,1)=X
C                  PRCORD(I500,2,I602,1)=Y
C                  PRCORD(I500,3,I602,1)=Z
C
C                  READ(80,*) SDUMMY,X,Y,Z
C                  WRITE(6,*) SDUMMY,X,Y,Z
C                  PRCORD(I500,1,I602,2)=X
C                  PRCORD(I500,2,I602,2)=Y
C                  PRCORD(I500,3,I602,2)=Z

C                  READ(80,*) SDUMMY,X,Y,Z
C                  WRITE(6,*) SDUMMY,X,Y,Z
C                  PRCORD(I500,1,I602,3)=X
C                  PRCORD(I500,2,I602,3)=Y
C                  PRCORD(I500,3,I602,3)=Z

C76                FORMAT(A2,4X,3(F20.10))
C500              CONTINUE

502     CONTINUE ! DO 502 I502=1,NUMPRR
 
       GLY_C = 0

       DO III = 1,NMRES

       IF (IRES(III).EQ.8) THEN
        X_MCP(9*(III-1)+1-(GLY_C)*3) = DBLE(PRCORD(III, 1, 1, 1))   !  CA X
        X_MCP(9*(III-1)+2-(GLY_C)*3) = DBLE(PRCORD(III, 2, 1, 1))   !  CA Y
        X_MCP(9*(III-1)+3-(GLY_C)*3) = DBLE(PRCORD(III, 3, 1, 1))   !  CA Z
C       X_MCP(9*(III-1)+4) = (PRCORD(III, 1, 1, 2))   !  CB X
C       X_MCP(9*(III-1)+5) = (PRCORD(III, 2, 1, 2))   !  CB Y
C       X_MCP(9*(III-1)+6) = (PRCORD(III, 3, 1, 2))   !  CB Z
        X_MCP(9*(III-1)+4-(GLY_C)*3) = DBLE(PRCORD(III, 1, 1, 3))   !   O X
        X_MCP(9*(III-1)+5-(GLY_C)*3) = DBLE(PRCORD(III, 2, 1, 3))   !   O Y
        X_MCP(9*(III-1)+6-(GLY_C)*3) = DBLE(PRCORD(III, 3, 1, 3))   !   O Z
        GLY_C = GLY_C +1
      ELSE
        X_MCP(9*(III-1)+1-GLY_C*3) = DBLE(PRCORD(III, 1, 1, 1))   !  CA X
        X_MCP(9*(III-1)+2-GLY_C*3) = DBLE(PRCORD(III, 2, 1, 1))   !  CA Y
        X_MCP(9*(III-1)+3-GLY_C*3) = DBLE(PRCORD(III, 3, 1, 1))   !  CA Z
        X_MCP(9*(III-1)+4-GLY_C*3) = DBLE(PRCORD(III, 1, 1, 2))   !  CB X
        X_MCP(9*(III-1)+5-GLY_C*3) = DBLE(PRCORD(III, 2, 1, 2))   !  CB Y
        X_MCP(9*(III-1)+6-GLY_C*3) = DBLE(PRCORD(III, 3, 1, 2))   !  CB Z
        X_MCP(9*(III-1)+7-GLY_C*3) = DBLE(PRCORD(III, 1, 1, 3))   !   O X
        X_MCP(9*(III-1)+8-GLY_C*3) = DBLE(PRCORD(III, 2, 1, 3))   !   O Y
        X_MCP(9*(III-1)+9-GLY_C*3) = DBLE(PRCORD(III, 3, 1, 3))   !   O Z
      ENDIF

      ENDDO

        QUENCH_CRD(:,:,:,:,1)=PRCORD
      
        ELSE
          DO I503=1, NQUENCH
           DO I502=1,NUMPRR
             READ(80,683)IDUMMY,IDUMMY,IDUMMY,Q_TEMP(I503),IDUMMY
683          FORMAT(3(I6,1X),F8.4,1X,I5,' STUCT SNAP T T TID')
C             WRITE(6,*)'QUENCH_T  TEMTUR I503' , Q_TEMP(I503)

               DO I500=1,NMRSS
C             WRITE(SO,*) I500,NMRSS,NUMPRR,NQUENCH,I503,TEMTUR(I503)

                  READ(80,76)
     *            (QUENCH_CRD(I500,I2,I502,1,I503),I2=1,3),
     *            (QUENCH_CRD(I500,I2,I502,2,I503),I2=1,3),
     *            (QUENCH_CRD(I500,I2,I502,3,I503),I2=1,3)
              ENDDO
            ENDDO
       ENDDO  ! I503 (LOOP OVER QUENCH)
         DO  I504=1, NQUENCH
C        SET UP NEW ANNEALLING SCHEDULE  
         INC=0
         DO 540 GRID_NUM=1,4
            IF( ITGRD(GRID_NUM).GT.0 )THEN
               RINC=Q_TEMP(I504)
               RINC=RINC/DBLE(ITGRD(GRID_NUM))
                DO 501 I501=1,ITGRD(GRID_NUM)
                  TEMTUR_QUENCH(INC+I501,I504)=Q_TEMP(I504) -
     *                             RINC*DBLE(I501-1)
C                  WRITE(6,*)'TEMTUR_QUENCH STEPS STRUCTURE',
C     *                 TEMTUR_QUENCH(INC+I501,I504),INC+I501,I504
  501          CONTINUE
               INC=INC + ITGRD(GRID_NUM)
               ENDIF  !  ITGRD(GRID_NUM)
540          CONTINUE  !   DO 540 GRID_NUM=1,4
             ENDDO   !  DO I504=1, NQUENCH
        ENDIF ! QUENCH

C     SET NUMBER IF INITIAL STRUCTURES TO NUMBER OF
C     STRUCTURES JUST READ IN

      PROCNT=NUMPRR

C     SEND MESSAGE ACKNOWLEDGING THAT 'OLD' STRUCTURES
C     ARE BEING USED AS STARTING STRUCTURES

CTEMPHACK
C      WRITE(OARCHV,450)NUMPRR
C  450 FORMAT(/'RESTART WITH ',I3,' PROTEINS')

       CLOSE(80)

C     ---------------------- DONE -----------------------

      RETURN
      END
