
C     --------------------- SAVSTR ----------------------

      SUBROUTINE SAVSTR(NMRES,NUMPRO,MAXPRO,
     *            MAXCRD,PRCORD,IRES,SAVE_NAME,OCONV)

C     --------------------------------------------------

C     SAVSTR SAVE PROTEIN STRUCTURES, EG. FINAL.PDB

C     ARGUMENTS:
C        NMRES  - NUMBER OF RESIDUES (I)
C        MAXSIZ - MAXIMUM NUMBER OF RESIDUES (I)
C        NUMPRO - NUMBER OF TRIAL PROTEIN STRUCTURES (I)
C        MAXPRO - MAXIMUM NUMBER OF TRIAL PROTEIN 
C                 STRUCTURES (I)
C        MAXCRD - MAXIMUM NUMBER OF ATOMS/ RESIDUE (I)
C        PRCORD - TRIAL STRUCTURES (I)
C        OCONV  - UNIT ID FOR FILE TO WHICH STRUCTURES
C                 ARE TO BE WRITTEN (I)

C     ---------------------------------------------------

      USE AMHGLOBALS,  ONLY: MAXSIZ


      IMPLICIT NONE

C     ARGUMENT DECLARATIONS:

       INTEGER NMRES,NUMPRO,MAXPRO,
     *           MAXCRD,OCONV


        DOUBLE PRECISION PRCORD(MAXSIZ,3,MAXPRO,MAXCRD)

        INTEGER IRES(MAXSIZ)

C     INTERNAL VARIABLES:

        DOUBLE PRECISION CPRCORD(MAXSIZ,3), NITCORD(MAXSIZ,3)

        INTEGER I_PRO, I_RES, I_AXIS, ATOM_NO

        CHARACTER*10 SAVE_NAME
        INTEGER NL

        CHARACTER*3 RES_TYPE(MAXSIZ)

        EXTERNAL GET_RES_NAME

C     --------------------- BEGIN -----------------------


        DO 11 NL = 10, 1, -1

          IF (SAVE_NAME(NL:NL) .NE. ' ') THEN

             GO TO 12

          END IF

11      CONTINUE
12      CONTINUE


         OPEN(UNIT=OCONV,FILE=SAVE_NAME(1:NL)//'.PDB',
     *        STATUS='NEW',FORM='FORMATTED')


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CALCULTATE N AND C' POSITIONS

      DO 490 I_PRO=1,NUMPRO
          DO 491 I_RES=1,NMRES
              DO 492 I_AXIS = 1,3
              
        CPRCORD(I_RES,I_AXIS)=0.4436538*PRCORD(I_RES,I_AXIS,I_PRO,1)
     *                   +0.2352006*PRCORD(I_RES+1,I_AXIS,I_PRO,1)
     *                   +0.3211455*PRCORD(I_RES,I_AXIS,I_PRO,3)
        NITCORD(I_RES+1,I_AXIS)=0.4831806*PRCORD(I_RES,I_AXIS,I_PRO,1)
     *                   +0.7032820*PRCORD(I_RES+1,I_AXIS,I_PRO,1)
     *                   -0.1864626*PRCORD(I_RES,I_AXIS,I_PRO,3)


492           CONTINUE
491        CONTINUE
490     CONTINUE


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  WRITE OUT FINAL STRUCTURE

        ATOM_NO = 1

      DO 502 I_PRO=1,NUMPRO
            DO 500 I_RES=1,NMRES

            CALL GET_RES_NAME(IRES(I_RES), RES_TYPE(I_RES))

CNNNNNNNNNNNNNNN
            IF( I_RES .GT. 1) THEN            ! WRITE N POSITION

         WRITE(OCONV,665) ATOM_NO, RES_TYPE(I_RES), I_RES,
     *        (NITCORD(I_RES, I_AXIS), I_AXIS =1,3), ATOM_NO

665        FORMAT('ATOM    ',I3,'  N   ', A3, '   ',I3,'    ',3(F8.3),
     *             '  1.00  0.00      TPDB ',I3)

            ATOM_NO = ATOM_NO + 1


            END IF

CCACACACACCACACACACA
         WRITE(OCONV,666) ATOM_NO, RES_TYPE(I_RES), I_RES, 
     *        (PRCORD(I_RES, I_AXIS, 1, 1), I_AXIS =1,3), ATOM_NO
666        FORMAT('ATOM    ',I3,'  CA  ', A3, '   ',I3,'    ',3(F8.3),
     *             '  1.00  0.00      TPDB ',I3)
            ATOM_NO = ATOM_NO + 1


CC'C'C'C'C'C'C'C'C'C'

            IF( I_RES .LT. NMRES) THEN            ! WRITE C' POSITION

         WRITE(OCONV,667) ATOM_NO, RES_TYPE(I_RES), I_RES,
     *        (CPRCORD(I_RES, I_AXIS), I_AXIS =1,3), ATOM_NO

667        FORMAT('ATOM    ',I3,'  C   ', A3, '   ',I3,'    ',3(F8.3),
     *             '  1.00  0.00      TPDB ',I3)

            ATOM_NO = ATOM_NO + 1

COOOOOOOOOOOOO                                  ! WRITE O POSITION

         WRITE(OCONV,668) ATOM_NO, RES_TYPE(I_RES), I_RES,
     *        (PRCORD(I_RES, I_AXIS, 1,3), I_AXIS =1,3), ATOM_NO

668        FORMAT('ATOM    ',I3,'  O   ', A3, '   ',I3,'    ',3(F8.3),
     *             '  1.00  0.00      TPDB ',I3)

            ATOM_NO = ATOM_NO + 1

            END IF


             IF (IRES(I_RES) .NE. 8) THEN
         WRITE(OCONV,669) ATOM_NO, RES_TYPE(I_RES), I_RES,
     *        (PRCORD(I_RES, I_AXIS, 1, 2), I_AXIS =1,3), ATOM_NO
669        FORMAT('ATOM    ',I3,'  CB  ', A3, '   ',I3,'    ',3(F8.3),
     *             '  1.00  0.00      TPDB ',I3)
            ATOM_NO = ATOM_NO + 1

             END IF 



  500       CONTINUE
  502 CONTINUE


      CLOSE(OCONV)

C     ---------------------- DONE -----------------------

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        SUBROUTINE GET_RES_NAME(RES_NUMBER, RES_NAME)



        IMPLICIT NONE


        INTEGER RES_NUMBER

        CHARACTER*3 RES_NAME

        IF (RES_NUMBER .EQ. 1) THEN 
                RES_NAME =  "ALA" 
        ENDIF
        IF (RES_NUMBER .EQ. 2) THEN 
                RES_NAME =  "ARG" 
        ENDIF
        IF (RES_NUMBER .EQ. 3) THEN 
                RES_NAME =  "ASN" 
        ENDIF
        IF (RES_NUMBER .EQ. 4) THEN 
                RES_NAME =  "ASP" 
        ENDIF
        IF (RES_NUMBER .EQ. 5) THEN 
                RES_NAME =  "CYS" 
        ENDIF
        IF (RES_NUMBER .EQ. 6) THEN 
                RES_NAME =  "GLN" 
        ENDIF
        IF (RES_NUMBER .EQ. 7) THEN 
                RES_NAME =  "GLU" 
        ENDIF
        IF (RES_NUMBER .EQ. 8) THEN 
                RES_NAME =  "GLY" 
        ENDIF
        IF (RES_NUMBER .EQ. 9) THEN 
                RES_NAME =  "HIS" 
        ENDIF
        IF (RES_NUMBER .EQ. 10) THEN 
                RES_NAME =  "ILE" 
        ENDIF
        IF (RES_NUMBER .EQ. 11) THEN 
                RES_NAME =  "LEU" 
        ENDIF
        IF (RES_NUMBER .EQ. 12) THEN 
                RES_NAME =  "LYS" 
        ENDIF
        IF (RES_NUMBER .EQ. 13) THEN 
                RES_NAME =  "MET" 
        ENDIF
        IF (RES_NUMBER .EQ. 14) THEN 
                RES_NAME =  "PHE" 
        ENDIF
        IF (RES_NUMBER .EQ. 15) THEN 
                RES_NAME =  "PRO" 
        ENDIF
        IF (RES_NUMBER .EQ. 16) THEN 
                RES_NAME =  "SER" 
        ENDIF
        IF (RES_NUMBER .EQ. 17) THEN 
                RES_NAME =  "THR" 
        ENDIF
        IF (RES_NUMBER .EQ. 18) THEN 
                RES_NAME =  "TRP" 
        ENDIF
        IF (RES_NUMBER .EQ. 19) THEN 
                RES_NAME =  "TYR" 
        ENDIF
        IF (RES_NUMBER .EQ. 20) THEN 
                RES_NAME =  "VAL" 
        ENDIF


        IF ((RES_NUMBER .GT. 20) .OR. (RES_NUMBER .LT. 1)) THEN
           WRITE(*,*) 'RESIDUE OUT OF RANGE'
        END IF
             

        RETURN

        END   

