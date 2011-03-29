
C     --------------------- SHAKE -----------------------

      SUBROUTINE SHAKE(MAXSIZ,PRCORD,QRCORD,SRCORD,
     *                 BONDLN,NUMPRO,JSTRT,JFINS,ZRCORD,
     *                 MAXSHK,TOLSHK,BDSHAK,ISHKIT,
     *                 MAXPRO,MAXCRD,MAXCNT,WORK1,WORK3,
     *                 WORK4,OARCHV)

C     ---------------------------------------------------

C     SHAKE ATTEMPTS TO SATISFY NEAREST-NEIGHBOR DISTANCE
C           CONSTRAINTS USING ITERATIVE PROCEDURE 
C           PROPOSED BY RYCKAERT ET AL., J. COMP. PHYS. 
C           23, 327-341 (1977).

C     ARGUMENTS:

C        MAXSIZ- MAXIMUM NUMBER OF PROTEIN RESIDUES
C        PRCORD- UPDATED COORDINATES SATISFYING BOND
C                CONSTRAINTS
C        QRCORD- PREVIOUS COORDINATES
C        SRCORD- SCRATCH SPACE
C        BONDLN- SPECIFIED BOND LENGTHS
C        NUMPRO- NUMBER OF CONFIGURATIONS
C        JSTRT - FIRST NONFIXED SITE
C        JFINS - LAST NONFIXED SITE
C        ZRCORD- WORK ARRAY
C        MAXSHK- MAXIMUM NUMBER OF SHAKE ITERATIONS
C        TOLSHK- SHAKE TOLERANCE
C        BDSHAK- FLAG TO ALERT DRIVER TO A INEFFECTIVE
C                SHAKE
C        ISHKIT- FLAG USED TO TRACK NUMBER OF SHAKE 
C                ITERATIONS
C        MAXPRO- MAXIMUM NUMBER OF PROTEINS

C     ---------------------------------------------------

      IMPLICIT NONE

C     ARGUMENT DECLARATIONS:

         INTEGER MAXSIZ,JSTRT,JFINS,NUMPRO,
     *           MAXSHK,ISHKIT,MAXPRO,MAXCNT,OARCHV,
     *           MAXCRD

         REAL PRCORD(MAXSIZ,3,MAXPRO,MAXCRD),
     *        QRCORD(MAXSIZ,3,MAXPRO,MAXCRD),
     *        SRCORD(MAXSIZ,3,MAXPRO),
     *        ZRCORD(MAXSIZ,3,MAXPRO),
     *        BONDLN(MAXSIZ,MAXCRD),TOLSHK,
     *        WORK1(MAXCNT),WORK3(MAXCNT),
     *        WORK4(MAXCNT)
    
         CHARACTER*3 RES_TYPE
 
C     INTERNAL VARIABLES:

         LOGICAL LBAD,BDSHAK

         INTEGER ISTRT

C        --- DO LOOP INDICES ---


        INTEGER I_AXIS, I_MAX_SHK_ITER, I_ODD_EVEN, 
     *          I_PRO, I_PRO_2, I_RES,I1,III,JJJ


C        --- IMPLIED DO LOOP INDICES ---

         REAL DIFF

C     REQUIRED SUBROUTINES

C     --------------------- BEGIN -----------------------

C     --- DIAGNOSTICS ---

C        ECHO SCALAR ARGUMENT ARGUMENTS

C         WRITE(OARCHV,200)MAXCRD,MAXCNT
C  200    FORMAT('MAXCRD ',I3,' MAXCNT ',I5)
C
CCC         CHECK THAT THERE IS ROOM FOR STORAGE OF R 
CCC         IN SRCORD(1,.,.)
CC
CC          IF( JSTRT.LE.1 )THEN
CC             WRITE(OARCHV,106)JSTRT
CC  106        FORMAT(/'SHAKE:JSTRT TOO BIG ',I3)
CC             STOP
CC          ENDIF

C        PRINT OUT CURRENT AND PREVIOUS BOND LENGTHS

C         WRITE(OARCHV,105)
C  105    FORMAT(/'SHAKE:PRCORD QRCORD')
C         DO 539 I_RES=JSTRT,JFINS
C            WORK1(1)=(PRCORD(I_RES,1,1,1)-
C     *                PRCORD(I_RES-1,1,1,1))**2 +
C     *               (PRCORD(I_RES,2,1,1)-
C     *                PRCORD(I_RES-1,2,1,1))**2 +
C     *               (PRCORD(I_RES,3,1,1)-
C     *                PRCORD(I_RES-1,3,1,1))**2 
C            WORK1(1)=SQRT(WORK1(1))

C            WORK1(2)=(QRCORD(I_RES,1,1,1)-
C     *                QRCORD(I_RES-1,1,1,1))**2 +
C     *               (QRCORD(I_RES,2,1,1)-
C     *                QRCORD(I_RES-1,2,1,1))**2 +
C     *               (QRCORD(I_RES,3,1,1)-
C     *                QRCORD(I_RES-1,3,1,1))**2 
C            WORK1(2)=SQRT(WORK1(2))

C            WRITE(OARCHV,104)I_RES,
C     *                      (PRCORD(I_RES,I1,1,1),I1=1,3),
C     *                      (QRCORD(I_RES,I1,1,1),I1=1,3)
C  104       FORMAT(I3,2(3(1X,1PE10.3),2X))

C            DIFF=ABS(BONDLN(I_RES,1)-WORK1(1))
C            IF( (DIFF.GT.1.0).OR.(DIFF.LE.5.0E-05) )
C     *         WRITE(OARCHV,191)I_RES,WORK1(1),WORK1(2),
C     *                          BONDLN(I_RES,1),DIFF
C  191          FORMAT('FRONT SHAKE ',I3,1X,4(1PE12.5,1X))

C  539    CONTINUE

C      ENDIF
C     --- END DIAGNOSTICS ---
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FIND THE R OF RYCKAERT ET AL.
      DO 521 I_PRO=1,NUMPRO
         DO 520 I_AXIS=1,3

C           COMPUTE R

            DO 513 I_RES=JSTRT,JFINS
               ZRCORD(I_RES,I_AXIS,I_PRO)=
     *            QRCORD(I_RES,I_AXIS,I_PRO,1) -
     *            QRCORD(I_RES-1,I_AXIS,I_PRO,1)
  513       CONTINUE
  520    CONTINUE
  521 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     --- DIAGNOSTICS ---

C      IF( IDIGNS )THEN

C        ECHO R

C          WRITE(OARCHV,*)I_RES,' PROTEIN: SRCORD'
C         DO 522 I_RES=JSTRT,JFINS
C            WRITE(30,135)I_RES,(SRCORD(I_RES,I1),I1=1,3)
C  135       FORMAT(I3,1X,3(1PE10.3,1X))
C  522    CONTINUE
C
C      ENDIF

C     SET ARRAY WORK4 TO THE (BOND LENGTHS)**2

      DO 530 I_RES=JSTRT,JFINS
         WORK4(I_RES)=BONDLN(I_RES,1)**2
  530 CONTINUE

C     LOOP OVER NUMPRO PROTEIN CONFIGURATIONS

      DO 501 I_PRO=1,NUMPRO

C        PERFORM MAXSHK ITERATIONS IN AN ATTEMPT 
C        TO SATISFY THE NEAREST-NEIGHBOR DISTANCE 
C        CONSTRAINTS

         DO 500 I_MAX_SHK_ITER=1,MAXSHK

C           PERFORM CHECKERBOARD BREAKUP, I.E., ANALYZE
C           CONSTRAINTS INDEPENDENTLY OF ONE ANOTHER



            DO 508 I_ODD_EVEN=1,2

               IF( I_ODD_EVEN.EQ.1 )THEN

C                 'ODD' CONSTRAINTS

                  ISTRT=JSTRT

               ELSE

C                 'EVEN' CONSTRAINTS

                  ISTRT=JSTRT + 1

               ENDIF





C              FIND R' OF RYCKAERT ET AL.

               DO 503 I_AXIS=1,3
                  DO 502 I_RES=JSTRT,JFINS
                     SRCORD(I_RES,I_AXIS,I_PRO)=
     *               PRCORD(I_RES,I_AXIS,I_PRO,1) -
     *               PRCORD(I_RES-1,I_AXIS,I_PRO,1)
  502             CONTINUE
  503          CONTINUE

C              FIND |R'|**2

               DO 504 I_RES=JSTRT,JFINS
                     WORK1(I_RES)=SRCORD(I_RES,1,I_PRO)**2 + 
     *                           SRCORD(I_RES,2,I_PRO)**2 +
     *                           SRCORD(I_RES,3,I_PRO)**2
  504          CONTINUE

C              FIND D**2 - R'**2

               DO 505 I_RES=JSTRT,JFINS
C                  WRITE(30,402)I_RES,WORK1(I_RES),WORK4(I_RES),
C     *                         ABS(WORK4(I_RES)-WORK1(I_RES))
C  402             FORMAT(I3,3(1X,1PE10.3))
                  WORK1(I_RES)=WORK4(I_RES) - WORK1(I_RES)
  505          CONTINUE

C              CHECK IF DONE


               DO 506 I_RES=JSTRT,JFINS

C                 IF ALL CONSTRAINTS NOT SATISFIED, THEN CONTINUE;
C                 OTHERWISE, CONSIDER NEXT CONFIGURATION

                  IF( ABS(WORK1(I_RES)).GT.TOLSHK )THEN
                     GO TO 522
                  ENDIF

  506          CONTINUE

C              DONE

C              INCREMENT VARIABLE USED TO TRACK THE NUMBER
C              OF REQUIRED ITERATIONS

               ISHKIT=ISHKIT + I_MAX_SHK_ITER - 1
               GO TO 501

  522          CONTINUE

C              FIND R.R'
              DO 507 I_RES=ISTRT,JFINS,2

              WORK3(I_RES)=ZRCORD(I_RES,1,I_PRO)*SRCORD(I_RES,1,I_PRO)
     *                   + ZRCORD(I_RES,2,I_PRO)*SRCORD(I_RES,2,I_PRO)
     *                   + ZRCORD(I_RES,3,I_PRO)*SRCORD(I_RES,3,I_PRO)

  507          CONTINUE

C              --- DIAGNOSTICS ---

C               IF( IDIGNS )THEN
C 
CC                 PRINT WORK3 IF TOO SMALL
C
C                  DO 561 I_RES=ISTRT,JFINS,2
C                     IF( ABS(WORK3(I_RES)).LT.EPSILN )THEN
CC                        DIFF=0.25*WORK1(I_RES)/WORK3(I_RES)
C                        WRITE(OARCHV,124)I_MAX_SHK_ITER,I_RES,WORK3(I_RES)
C                        WRITE(OARCHV,711)
C     *                  (ZRCORD(I_RES,I1,I_PRO),I1=1,3),
C     *                  (SRCORD(I_RES,I1,I_PRO),I1=1,3)
C  711                   FORMAT('Z + S',6(1X,1PE10.3))
C  124                   FORMAT('SHAKE ITER ',I4,' SITE ',I3,
C     *                         ' WORK3 ',1PE10.3,' G ',1PE10.3)
C                     ENDIF
C  561             CONTINUE
C
C               ENDIF

C              --- END DIAGNOSTICS ---

C              COMPUTE G 

               DO 509 I_RES=ISTRT,JFINS,2
                  WORK3(I_RES)=0.25*WORK1(I_RES)/WORK3(I_RES)
  509          CONTINUE 

C              UPDATE COORDINATES

               DO 511 I_AXIS=1,3

                  IF( (ISTRT.GT.2).AND.(I_ODD_EVEN.EQ.1) )THEN

                     DO 510 I_RES=ISTRT+1,JFINS-1,2
                        PRCORD(I_RES,I_AXIS,I_PRO,1)=
     *                  PRCORD(I_RES,I_AXIS,I_PRO,1) -
     *                  ZRCORD(I_RES+1,I_AXIS,I_PRO)*WORK3(I_RES+1)
  510                CONTINUE

                  ELSE

                     DO 537 I_RES=ISTRT-1,JFINS-1,2
                        PRCORD(I_RES,I_AXIS,I_PRO,1)=
     *                  PRCORD(I_RES,I_AXIS,I_PRO,1) -
     *                  ZRCORD(I_RES+1,I_AXIS,I_PRO)*WORK3(I_RES+1)
  537                CONTINUE

                  ENDIF

                  DO 524 I_RES=ISTRT,JFINS,2
                     PRCORD(I_RES,I_AXIS,I_PRO,1)=
     *               PRCORD(I_RES,I_AXIS,I_PRO,1) +
     *               ZRCORD(I_RES,I_AXIS,I_PRO)*WORK3(I_RES)
  524             CONTINUE
 
                  IF( (ISTRT.GT.2).AND.(I_ODD_EVEN.EQ.1) )THEN
                     PRCORD(ISTRT,I_AXIS,I_PRO,1)=
     *               PRCORD(ISTRT,I_AXIS,I_PRO,1) +
     *               ZRCORD(ISTRT,I_AXIS,I_PRO)*WORK3(ISTRT)
                  ENDIF

  511          CONTINUE

C               IF( I_ODD_EVEN.EQ.1 )THEN 
C                  WRITE(30,918)I_MAX_SHK_ITER,ISTRT-1,
C     *            (PRCORD(ISTRT-1,I1,I_PRO,1),
C     *                    WORK8(I1),I1=1,3)
C  918             FORMAT(2(I3,1X),6(1PE10.3,1X))
C               ENDIF

  508       CONTINUE

  500    CONTINUE

  501 CONTINUE

      CONTINUE



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     --- DIAGNOSTICS ---

C      IF( IDIGNS )THEN

         LBAD=.FALSE.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        DETERMINE WHICH IF ANY OF THE CONSTRAINTS
C        ARE NOT SATISFIED

         DO 512 I_PRO_2=1,MIN(NUMPRO,2)
            DO 514 I_AXIS=1,3
               DO 523 I_RES=JSTRT,JFINS
                  SRCORD(I_RES,I_AXIS,I_PRO_2)=
     *               PRCORD(I_RES,I_AXIS,I_PRO_2,1) -
     *               PRCORD(I_RES-1,I_AXIS,I_PRO_2,1)
  523          CONTINUE
  514       CONTINUE



            DO 515 I_RES=JSTRT,JFINS
               WORK1(I_RES)=SRCORD(I_RES,1,I_PRO_2)**2 +
     *                     SRCORD(I_RES,2,I_PRO_2)**2 +
     *                     SRCORD(I_RES,3,I_PRO_2)**2
  515       CONTINUE

            DO 516 I_RES=JSTRT,JFINS
               WORK1(I_RES)=SQRT(WORK1(I_RES))
  516       CONTINUE

            DO 517 I_RES=JSTRT,JFINS 
               DIFF=ABS(WORK1(I_RES) - BONDLN(I_RES,1))
               IF( DIFF.GT.TOLSHK )THEN
                  LBAD=.TRUE.
                  GO TO 519
               ENDIF
  517       CONTINUE


  519       CONTINUE
           
            IF( LBAD )THEN

C              --- ALL CONSTRAINTS NOT SATISFIED 
C                  AFTER MAXSHK ITERATIONS ---

               WRITE(OARCHV,103)
  103          FORMAT(/'NOT ALL BONDS SATISFIED IN SHAKE')

               DO 518 I_RES=JSTRT,JFINS 
                  DIFF=ABS(WORK1(I_RES) - BONDLN(I_RES,1))
                  IF( DIFF.GT.TOLSHK )THEN
                     WRITE(OARCHV,110)I_PRO_2,I_RES,WORK1(I_RES),
     *                                BONDLN(I_RES,1),DIFF
  110                FORMAT('PRO ',I2,' SITE ',I3,' D(CAL) ',
     *                      1PE10.3,' D(EXACT) ',1PE10.3,
     *                      ' DIFF ',1PE10.3)
                  ENDIF
  518          CONTINUE
   
               BDSHAK=.TRUE.
               DO 533 I_AXIS=1,3
                  DO 525 I_RES=JSTRT,JFINS
                     PRCORD(I_RES,I_AXIS,I_PRO_2,1)=
     *               QRCORD(I_RES,I_AXIS,I_PRO_2,1)
  525             CONTINUE
  533          CONTINUE
            ELSE

C               WRITE(30,*)'SHAKE OK ',I_MAX_SHK_ITER-1

            ENDIF
  512    CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
C      ENDIF

C     --- END DIAGNOSTICS ---

C     ---------------------- DONE ----------------------
     
      RETURN
      END
