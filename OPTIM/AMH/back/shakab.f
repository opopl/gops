
C     --------------------- SHAKAB -----------------------

      SUBROUTINE SHAKAB(MAXSIZ,PRCORD,QRCORD,SRCORD,
     *                  BONDLN,NUMPRO,JSTRT,JFINS,ZRCORD,
     *                  MAXSHK,TOLSHK,BDSHAK,ISHKIT,
     *                  MAXPRO,MAXCRD,MAXCNT,WORK1,WORK3,
     *                  WORK4,OARCHV,IRES)

C     ---------------------------------------------------

C     SHAKAB ATTEMPTS TO SATISFY ALPHA-BETA DISTANCE 
C            CONSTRAINTS USING THE ITERATIVE PROCEDURE 
C            PROPOSED BY RYCKAERT ET AL., J. COMP. PHYS. 
C            23, 327-341 (1977).

C     ARGUMENTS:

C        MAXSIZ- MAXIMUM NUMBER OF PROTEIN RESIDUES (I)
C        PRCORD- UPDATED COORDINATES SATISFYING BOND
C                CONSTRAINTS (I,O)
C        QRCORD- PREVIOUS COORDINATES (I)
C        SRCORD- SCRATCH SPACE (I)
C        BONDLN- SPECIFIED BOND LENGTHS (I)
C        NUMPRO- NUMBER OF CONFIGURATIONS (I)
C        JSTRT - FIRST NONFIXED SITE (I)
C        JFINS - LAST NONFIXED SITE (I)
C        ZRCORD- WORK ARRAY (I)
C        MAXSHK- MAXIMUM NUMBER OF SHAKE ITERATIONS (I)
C        TOLSHK- SHAKE TOLERANCE (I)
C        BDSHAK- FLAG TO ALERT DRIVER TO A INEFFECTIVE
C                SHAKE (O)
C        ISHKIT- FLAG USED TO TRACK NUMBER OF SHAKE 
C                ITERATIONS (I)
C        MAXPRO- MAXIMUM NUMBER OF PROTEINS (I)
C        MAXCRD- MAXIMUM NUMBER OF COORDINATES (I)

C     ---------------------------------------------------

      IMPLICIT NONE

C     ARGUMENT DECLARATIONS:

         INTEGER MAXSIZ,JSTRT,JFINS,NUMPRO,
     *           MAXSHK,ISHKIT,MAXPRO,MAXCRD,MAXCNT,
     *           OARCHV,IRES(MAXSIZ)

         REAL PRCORD(MAXSIZ,3,MAXPRO,MAXCRD),
     *        QRCORD(MAXSIZ,3,MAXPRO,MAXCRD),
     *        SRCORD(MAXSIZ,3,MAXPRO),
     *        ZRCORD(MAXSIZ,3,MAXPRO),
     *        BONDLN(MAXSIZ,MAXCRD),TOLSHK,
     *        WORK1(MAXCNT),WORK3(MAXCNT),
     *        WORK4(MAXCNT)
     
C     INTERNAL VARIABLES:

         LOGICAL LBAD,BDSHAK

         CHARACTER*3 RES_TYPE

 
C        --- IMPLIED DO LOOP INDICES ---

        INTEGER I_AXIS, I_PRO, I_RES, I_SHK_ITER,III,JJJ

         REAL DIFF

C     REQUIRED SUBROUTINES

C     --------------------- BEGIN -----------------------

C     --- DIAGNOSTICS ---

C      IF( IDIGNS )THEN
C         
CC        ECHO SCALAR ARGUMENT ARGUMENTS
C

CC         CHECK THAT THERE IS ROOM FOR STORAGE OF R 
CC         IN SRCORD(1,.,.)
C
C          IF( JSTRT.LE.0 )THEN
C             WRITE(OARCHV,106)JSTRT
C  106        FORMAT(/'SHAKE:JSTRT TOO BIG ',I3)
C             STOP
C          ENDIF

C        PRINT OUT CURRENT AND PREVIOUS BOND LENGTHS

C         WRITE(OARCHV,105)
C  105    FORMAT(/'SHAKAB:PRCORD QRCORD')
C        DO 539 I_RES=JSTRT,JFINS
C           WORK1(1)=(PRCORD(I_RES,1,1,1)-
C    *                PRCORD(I_RES,1,1,2))**2 +
C    *               (PRCORD(I_RES,2,1,1)-
C    *                PRCORD(I_RES,2,1,2))**2 +
C    *               (PRCORD(I_RES,3,1,1)-
C    *                PRCORD(I_RES,3,1,2))**2 
C           WORK1(1)=SQRT(WORK1(1))

C           WORK1(2)=(QRCORD(I_RES,1,1,1)-
C    *                QRCORD(I_RES,1,1,2))**2 +
C    *               (QRCORD(I_RES,2,1,1)-
C    *                QRCORD(I_RES,2,1,2))**2 +
C    *               (QRCORD(I_RES,3,1,1)-
C    *                QRCORD(I_RES,3,1,2))**2 
C           WORK1(2)=SQRT(WORK1(2))

C            WRITE(OARCHV,104)I_RES,
C     *                      (PRCORD(I_RES,I1,1,1),I1=1,3),
C     *                      (PRCORD(I_RES,I1,1,2),I1=1,3)
C  104       FORMAT(I3,2(3(1X,1PE10.3),2X))

C            DIFF=ABS(BONDLN(I_RES,2)-WORK1(1))
C            IF( (DIFF.GT.1.0) )
C     *         WRITE(OARCHV,191)I_RES,WORK1(1),WORK1(2),
C     *                          BONDLN(I_RES,2),DIFF
C  191          FORMAT('FRONT SHAKE ',I3,1X,4(1PE12.5,1X))

C  539    CONTINUE

CC        CHECK IF NEXT-NEAREST NEIGHBOR CONSTRAINTS ARE 
CC        SATISFIED; IF NOT, THEN PRINT OUT OFFENDING SITES
C
C         CALL DIST2(MAXSIZ,JSTRT,JFINS,PRCORD,DIFF)
C
C      ENDIF

C     --- END DIAGNOSTICS ---

C     FIND THE R OF RYCKAERT ET AL.

      DO 521 I_PRO=1,NUMPRO
         DO 520 I_AXIS=1,3

C           COMPUTE R

            DO 513 I_RES=JSTRT,JFINS
               IF( IRES(I_RES).NE.8 )THEN
                 ZRCORD(I_RES,I_AXIS,I_PRO)=QRCORD(I_RES,I_AXIS,I_PRO,2)
     *                                 - QRCORD(I_RES,I_AXIS,I_PRO,1)
               ELSE
                  ZRCORD(I_RES,I_AXIS,I_PRO)=0.0
               ENDIF
  513       CONTINUE

  520    CONTINUE
  521 CONTINUE


        
C     --- DIAGNOSTICS ---

C     ECHO R

C       WRITE(OARCHV,811)
C 811   FORMAT(/'SHAKAB:  PROTEIN: ZRCORD')
C      DO 522 I_RES=JSTRT,JFINS
C         WRITE(30,135)I_RES,(ZRCORD(I_RES,I1,1),I1=1,3)
C  13     FORMAT(I3,1X,3(1PE10.3,1X))
C  52  CONTINUE

C     SET ARRAY WORK4 TO THE (BOND LENGTHS)**2

      DO 530 I_RES=JSTRT,JFINS
         WORK4(I_RES)=BONDLN(I_RES,2)**2
  530 CONTINUE

C     LOOP OVER NUMPRO PROTEIN CONFIGURATIONS

      DO 501 I_PRO=1,NUMPRO

C        SET BETA-GLYCINE TO ALPHA-GLYCINE COORDINATES

         DO 534 I_RES=JSTRT,JFINS
            IF( IRES(I_RES).EQ.8 )THEN
               PRCORD(I_RES,1,I_PRO,2)=PRCORD(I_RES,1,I_PRO,1)
               PRCORD(I_RES,2,I_PRO,2)=PRCORD(I_RES,2,I_PRO,1)
               PRCORD(I_RES,3,I_PRO,2)=PRCORD(I_RES,3,I_PRO,1)
            ENDIF
  534    CONTINUE

C        PERFORM MAXSHK ITERATIONS IN AN ATTEMPT TO SATISFY THE 
C         NEAREST-NEIGHBOR DISTANCE CONSTRAINTS
 
         DO 500 I_SHK_ITER=1,MAXSHK

C           PERFORM CHECKERBOARD BREAKUP, I.E., ANALYZE
C           CONSTRAINTS INDEPENDENTLY OF ONE ANOTHER

C           FIND R' OF RYCKAERT ET AL.

            DO 503 I_AXIS=1,3
               DO 502 I_RES=JSTRT,JFINS
                  SRCORD(I_RES,I_AXIS,I_PRO)=
     *            PRCORD(I_RES,I_AXIS,I_PRO,2) -
     *            PRCORD(I_RES,I_AXIS,I_PRO,1)
  502          CONTINUE
  503       CONTINUE

C           FIND |R'|**2

            DO 504 I_RES=JSTRT,JFINS
               WORK1(I_RES)=SRCORD(I_RES,1,I_PRO)**2 + 
     *                     SRCORD(I_RES,2,I_PRO)**2 +
     *                     SRCORD(I_RES,3,I_PRO)**2
  504       CONTINUE

C           FIND D**2 - R'**2

            DO 505 I_RES=JSTRT,JFINS
               WORK1(I_RES)=WORK4(I_RES) - WORK1(I_RES)
  505       CONTINUE

C           CHECK IF DONE

            DO 506 I_RES=JSTRT,JFINS

C              IF ALL CONSTRAINTS NOT SATISFIED, THEN CONTINUE;
C              OTHERWISE, CONSIDER NEXT CONFIGURATION

               IF( (ABS(WORK1(I_RES)).GT.TOLSHK) )THEN
                     GO TO 522
               ENDIF

  506       CONTINUE

C           DONE

C           INCREMENT VARIABLE USED TO TRACK THE NUMBER
C           OF REQUIRED ITERATIONS

            ISHKIT=ISHKIT + I_SHK_ITER - 1
            GO TO 501

  522       CONTINUE

C           FIND R.R'

            DO 507 I_RES=JSTRT,JFINS

               WORK3(I_RES)=ZRCORD(I_RES,1,I_PRO)*
     *                     SRCORD(I_RES,1,I_PRO) +
     *                     ZRCORD(I_RES,2,I_PRO)*
     *                     SRCORD(I_RES,2,I_PRO) +
     *                     ZRCORD(I_RES,3,I_PRO)*
     *                     SRCORD(I_RES,3,I_PRO)

  507       CONTINUE


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C           --- DIAGNOSTICS ---
C           PRINT WORK3 IF TOO SMALL
C            DO 561 I_RES=JSTRT,JFINS
C               IF( ABS(WORK3(I_RES)).LT.EPSILN )THEN
C                  DIFF=0.25*WORK1(I_RES)/WORK3(I_RES)
C                  WRITE(OARCHV,124)I_SHK_ITER,I_RES,WORK3(I_RES)
C                  WRITE(OARCHV,711)
C     *            (ZRCORD(I_RES,I1,I_PRO),I1=1,3),
C     *            (SRCORD(I_RES,I1,I_PRO),I1=1,3)
C  711             FORMAT('Z + S',6(1X,1PE10.3))
C  124             FORMAT('SHAKE ITER ',I4,' SITE ',I3,
C     *                   ' WORK3 ',1PE10.3,' G ',1PE10.3)
C               ENDIF
C  561       CONTINUE
C           --- END DIAGNOSTICS ---
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C           COMPUTE G 

            DO 509 I_RES=JSTRT,JFINS
               IF( IRES(I_RES).NE.8 )THEN
                  WORK3(I_RES)=0.25*WORK1(I_RES)/WORK3(I_RES)
               ELSE
                  WORK3(I_RES)=0.0
               ENDIF
  509       CONTINUE 



C           UPDATE COORDINATES

            DO 511 I_AXIS=1,3

               DO 510 I_RES=JSTRT,JFINS
                  PRCORD(I_RES,I_AXIS,I_PRO,2)=
     *            PRCORD(I_RES,I_AXIS,I_PRO,2) +
     *            ZRCORD(I_RES,I_AXIS,I_PRO)*WORK3(I_RES)
  510          CONTINUE

               DO 524 I_RES=JSTRT,JFINS
                  PRCORD(I_RES,I_AXIS,I_PRO,1)=
     *            PRCORD(I_RES,I_AXIS,I_PRO,1) -
     *            ZRCORD(I_RES,I_AXIS,I_PRO)*WORK3(I_RES)
  524          CONTINUE
 
  511       CONTINUE

  500    CONTINUE

  501 CONTINUE

      CONTINUE

C     --- DIAGNOSTICS ---

      LBAD=.FALSE.

C     DETERMINE WHICH IF ANY OF THE CONSTRAINTS
C     ARE NOT SATISFIED

      DO 512 I_PRO=1,NUMPRO
         DO 514 I_AXIS=1,3
            DO 523 I_RES=JSTRT,JFINS
               SRCORD(I_RES,I_AXIS,I_PRO)=PRCORD(I_RES,I_AXIS,I_PRO,2)
     *                              - PRCORD(I_RES,I_AXIS,I_PRO,1)
  523       CONTINUE
  514    CONTINUE
         DO 515 I_RES=JSTRT,JFINS
            WORK1(I_RES)=SRCORD(I_RES,1,I_PRO)**2 +
     *                  SRCORD(I_RES,2,I_PRO)**2 +
     *                  SRCORD(I_RES,3,I_PRO)**2
  515    CONTINUE
         DO 516 I_RES=JSTRT,JFINS
            WORK1(I_RES)=SQRT(WORK1(I_RES))
  516    CONTINUE


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO 517 I_RES=JSTRT,JFINS 
            DIFF=ABS(WORK1(I_RES) - BONDLN(I_RES,2))
            IF( DIFF.GT.TOLSHK )THEN
               LBAD=.TRUE.
               GO TO 519
            ENDIF
  517    CONTINUE
  519    CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

           
         IF( LBAD )THEN


C           --- ALL CONSTRAINTS NOT SATISFIED 
C               AFTER MAXSHK ITERATIONS ---

            WRITE(OARCHV,103)
  103       FORMAT(/'NOT ALL BONDS SATISFIED IN SHAKAB')

            DO 518 I_RES=JSTRT,JFINS 
               DIFF=ABS(WORK1(I_RES) - BONDLN(I_RES,2))
               IF( DIFF.GT.TOLSHK )THEN
                  WRITE(OARCHV,110)I_PRO,I_RES,WORK1(I_RES),
     *                             BONDLN(I_RES,2),DIFF
  110             FORMAT('PRO ',I2,' SITE ',I3,' D(CAL) ',
     *                   1PE10.3,' D(EXACT) ',1PE10.3,
     *                   ' DIFF ',1PE10.3)
               ENDIF
  518       CONTINUE
   
            BDSHAK=.TRUE.
            DO 533 I_AXIS=1,3
               DO 525 I_RES=JSTRT,JFINS
                  PRCORD(I_RES,I_AXIS,I_PRO,2)=
     *            QRCORD(I_RES,I_AXIS,I_PRO,2)
  525             CONTINUE
  533       CONTINUE
         ELSE

         ENDIF
  512 CONTINUE
 
C     --- END DIAGNOSTICS ---

C     ---------------------- DONE ----------------------

 
      RETURN
      END
