
C     --------------------- SHAKOX -----------------------

      SUBROUTINE SHAKOX(PRCORD,QRCORD,SRCORD,
     *                 NUMPRO,JSTRT,JFINS,
     *                 MAXSHK,TOLSHK,ISHKIT,
     *                 MAXPRO,MAXCRD)

C     ---------------------------------------------------

C     SHAKOX ATTEMPTS TO SATISFY NEAREST-NEIGHBOR DISTANCE
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
C        MAXSHK- MAXIMUM NUMBER OF SHAKE ITERATIONS
C        TOLSHK- SHAKE TOLERANCE
C                SHAKE
C        ISHKIT- FLAG USED TO TRACK NUMBER OF SHAKE 
C                ITERATIONS
C        MAXPRO- MAXIMUM NUMBER OF PROTEINS
C        I_AXIS        - INDEX OVER AXES
C        I_ODD_EVEN        - INDEX OVER COORDINATE TYPES CA AND CB
C        I_PRO        - INDEX OVER PROTEINS IN ENSEMBLE
C        I_RES        - INDEX OVER RESIDUES
C        I_SHAKE        - INDEX OVER SHAKE ITERATIONS
C
C
C
C     SKIP ADDITIONS
C       EQ_DIST_SQ        EQUILIBRIUM DISTANCE SQUARED
C       
C
C     ---------------------------------------------------

      USE GLOBALS, ONLY:SO,MAXSIZ,MAXCNT

      IMPLICIT NONE

C     ARGUMENT DECLARATIONS:

         INTEGER JSTRT,JFINS,NUMPRO,
     *           MAXSHK,ISHKIT,MAXPRO,
     *           MAXCRD,MAXPRO1

         PARAMETER (MAXPRO1=25)

         REAL PRCORD(MAXSIZ,3,MAXPRO,MAXCRD),
     *        QRCORD(MAXSIZ,3,MAXPRO,MAXCRD),
     *        SRCORD(MAXSIZ,3,MAXPRO),
     *        TOLSHK
     
C     INTERNAL VARIABLES:

         CHARACTER*3 RES_TYPE


         INTEGER I_AXIS, I_ODD_EVEN, I_PRO, I_RES, I_SHAKE,III,JJJ

         REAL MYST(MAXCNT), MYST2(MAXCNT)

C        --- IMPLIED DO LOOP INDICES ---

         REAL EQDIST(2),LNGTH(MAXSIZ,3,MAXPRO1,2)
         REAL EQ_DIST_SQ(2)

C     --------------------- BEGIN -----------------------

CCCCCCCCCCCCCCCCCCCCC


        EQDIST(1)=2.42677
        EQDIST(2)=2.82146
        EQ_DIST_SQ(1)=EQDIST(1)*EQDIST(1)
        EQ_DIST_SQ(2)=EQDIST(2)*EQDIST(2)

C     FIND THE R OF RYCKAERT ET AL.

      DO 521 I_PRO=1,NUMPRO
         DO 520 I_AXIS=1,3

C           COMPUTE R

          DO 513 I_RES=JSTRT,JFINS-1
            LNGTH(I_RES,I_AXIS,I_PRO,1)=
     *          QRCORD(I_RES,I_AXIS,I_PRO,3)-QRCORD(I_RES,I_AXIS,I_PRO,1)
513     CONTINUE

          DO 514 I_RES=JSTRT,JFINS-1
            LNGTH(I_RES,I_AXIS,I_PRO,2)=
     *          QRCORD(I_RES+1,I_AXIS,I_PRO,1)-QRCORD(I_RES,I_AXIS,I_PRO,3)
  514     CONTINUE

  520    CONTINUE
  521 CONTINUE


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C     LOOP OVER NUMPRO PROTEIN CONFIGURATIONS
      DO 501 I_PRO=1,NUMPRO

C        PERFORM MAXSHK ITERATIONS IN AN ATTEMPT TO SATISFY THE 
C         NEAREST-NEIGHBOR DISTANCE CONSTRAINTS

         DO 500 I_SHAKE=1,MAXSHK

C           PERFORM CHECKERBOARD BREAKUP, I.E., ANALYZE
C           CONSTRAINTS INDEPENDENTLY OF ONE ANOTHER

            DO 508 I_ODD_EVEN=1,2

C              FIND R' OF RYCKAERT ET AL.

            DO 503 I_AXIS=1,3
              DO 502 I_RES=JSTRT,JFINS-1
               IF (I_ODD_EVEN.EQ.1) THEN
                 SRCORD(I_RES,I_AXIS,I_PRO)=
     *           PRCORD(I_RES,I_AXIS,I_PRO,3)-PRCORD(I_RES,I_AXIS,I_PRO,1)
               ELSE
                 SRCORD(I_RES,I_AXIS,I_PRO)=
     *           PRCORD(I_RES+1,I_AXIS,I_PRO,1)-PRCORD(I_RES,I_AXIS,I_PRO,3) 
               ENDIF
502           CONTINUE
503         CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

C              FIND |R'|**2

               DO 504 I_RES=JSTRT,JFINS-1
                     MYST2(I_RES)=SRCORD(I_RES,1,I_PRO)**2 + 
     *                           SRCORD(I_RES,2,I_PRO)**2 +
     *                           SRCORD(I_RES,3,I_PRO)**2

C                     WRITE(SO,*) 'MYST2 ', MYST2(I_RES)
504          CONTINUE

C                     WRITE(SO,*) 'MYST2 ', MYST2(75)

C              FIND D**2 - R'**2

               DO 505 I_RES=JSTRT,JFINS-1
                  MYST2(I_RES)=EQ_DIST_SQ(I_ODD_EVEN) - MYST2(I_RES)
505          CONTINUE

C              CHECK IF DONE

               DO 506 I_RES=JSTRT,JFINS-1

C                 IF ALL CONSTRAINTS NOT SATISFIED, THEN CONTINUE;
C                 OTHERWISE, CONSIDER NEXT CONFIGURATION

                  IF( ABS(MYST2(I_RES)).GT.TOLSHK )THEN
                     GO TO 522
                  ENDIF

  506          CONTINUE

C              DONE

C              INCREMENT VARIABLE USED TO TRACK THE NUMBER
C              OF REQUIRED ITERATIONS

               ISHKIT=ISHKIT + I_SHAKE - 1

               GO TO 501  ! GO TO LOOP FOR THE NEXT PROTIEN

  522          CONTINUE

C              FIND R.R'

               DO 507 I_RES=JSTRT,JFINS-1

                  MYST(I_RES)=LNGTH(I_RES,1,I_PRO,I_ODD_EVEN)
     *                            *SRCORD(I_RES,1,I_PRO) +
     *                        LNGTH(I_RES,2,I_PRO,I_ODD_EVEN)
     *                            *SRCORD(I_RES,2,I_PRO) +
     *                        LNGTH(I_RES,3,I_PRO,I_ODD_EVEN)
     *                            *SRCORD(I_RES,3,I_PRO)

507          CONTINUE

C              COMPUTE G 

               DO 509 I_RES=JSTRT,JFINS-1
                  IF (MYST(I_RES) .NE. 0.0) THEN
                  MYST(I_RES)= (0.25*MYST2(I_RES)/MYST(I_RES))
                  END IF
  509          CONTINUE 

               DO 511 I_AXIS=1,3

                  DO 524 I_RES=JSTRT,JFINS-1
                    IF (I_ODD_EVEN.EQ.1) THEN
                      PRCORD(I_RES,I_AXIS,I_PRO,3)=
     *                PRCORD(I_RES,I_AXIS,I_PRO,3) +
     *                LNGTH(I_RES,I_AXIS,I_PRO,1)*MYST(I_RES)
 
                      PRCORD(I_RES,I_AXIS,I_PRO,1)=
     *                PRCORD(I_RES,I_AXIS,I_PRO,1) -
     *                LNGTH(I_RES,I_AXIS,I_PRO,1)*MYST(I_RES)
                    ELSE
                           PRCORD(I_RES+1,I_AXIS,I_PRO,1)=
     *                PRCORD(I_RES+1,I_AXIS,I_PRO,1) +
     *                LNGTH(I_RES,I_AXIS,I_PRO,2)*MYST(I_RES)

                      PRCORD(I_RES,I_AXIS,I_PRO,3)=
     *                PRCORD(I_RES,I_AXIS,I_PRO,3) -
     *                LNGTH(I_RES,I_AXIS,I_PRO,2)*MYST(I_RES) 
                    ENDIF
  524             CONTINUE

  511          CONTINUE   ! END OF LOOP OVER AXIS

  508       CONTINUE        ! END OF LOOP OVER ODD, EVEN

  500    CONTINUE        ! END OF LOOP OVER MAX SHAKE ITERATIONS

        WRITE(SO,*) 'FAILED TO PASS MAX ITERATIONS IN SHAKOX'
        STOP 'PASSED MAX ITERATIONS IN SHAKOX'

501     CONTINUE                ! END OF LOOP OVER PROTEINS
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CONTINUE

C     ---------------------- DONE ----------------------
     
      RETURN
      END
