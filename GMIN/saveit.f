C   GMIN: A PROGRAM FOR FINDING GLOBAL MINIMA
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF GMIN.
C
C   GMIN IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   GMIN IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
      SUBROUTINE GSAVEIT(EREAL,P,NP)
      USE COMMONS
      USE QMODULE
      IMPLICIT NONE


      INTEGER J1, J2, J3, NP, NQTOT, NPCALL, CSMIT
      DOUBLE PRECISION EREAL,P(3*NATOMS), AVVAL, CSMRMS

      COMMON /TOT/ NQTOT
      COMMON /PCALL/ NPCALL
C
C  SAVE THE LOWEST NSAVE DISTINGUISHABLE CONFIGURATIONS.
C
C     WRITE(*,'(A,12E15.7)') 'EREAL,ECONV,QMIN',EREAL,ECONV,(QMIN(J1),J1=1,NSAVE)
      DO J1=1,NSAVE
         IF (DABS(EREAL-QMIN(J1)).LT.ECONV) THEN
C
C  THESE ARE PROBABLY THE SAME - BUT JUST TO MAKE SURE WE SAVE THE 
C  LOWEST.
C
            IF (EREAL.LT.QMIN(J1)) THEN
               QMIN(J1)=EREAL
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=P(J2)
               ENDDO
            ENDIF
            GOTO 10
         ENDIF
         IF (EREAL.LT.QMIN(J1)) THEN

            J2=NSAVE
20          CONTINUE
      
            IF (NSAVE.GT.1) THEN
               QMIN(J2)=QMIN(J2-1)
               FF(J2)=FF(J2-1)
               DO J3=1,3*NATOMS
                  QMINP(J2,J3)=QMINP(J2-1,J3)
               ENDDO

               J2=J2-1
               IF (J2.GE.J1+1) GOTO 20
            ENDIF

            QMIN(J1)=EREAL
            FF(J1)=NQ(NP)
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=P(J2)
            ENDDO

            GOTO 10
         ENDIF
      ENDDO

10    CONTINUE

      DO J1=1,NTARGETS
         IF (EREAL-TARGETS(J1).LT.ECONV) THEN
            IF (NPAR.LT.2) THEN
               WRITE(MYUNIT,'(2(A,I15),A)') 'SAVEIT> TARGET HIT AFTER ',NQTOT,' QUENCHES ',NPCALL,' FUNCTION CALLS.'
               WRITE(MYUNIT,'(2(A,F20.10))') 'SAVEIT> ENERGY=',EREAL,' TARGET=',TARGETS(J1)
            ELSE
               WRITE(MYUNIT,'(A,I1,A,I2,2(A,I15),A)') '[',NP,']SAVEIT> TARGET HIT IN PARALLEL RUN ',NP,
     1                    ' AFTER ',NQTOT,' QUENCHES ',NPCALL,' FUNCTION CALLS.'
               WRITE(MYUNIT,'(A,I1,2(A,F20.10))') '[',NP,']SAVEIT> ENERGY=',EREAL,' TARGET=',TARGETS(J1)
            ENDIF
            HIT=.TRUE.
         ENDIF
      ENDDO

      RETURN
      END
C
C  CSW34> NEW SUBROUTINE TO SAVE THE LOWEST NSAVEINTE CONFIGURATIONS INDEXED BY INTERACTION ENTHALPY BETWEEN PROTEIN AND LIGAND
C
      SUBROUTINE A9INTESAVEIT(INTEREAL,P,NP)
      USE COMMONS
      USE QMODULE
      IMPLICIT NONE


      INTEGER J1, J2, J3, NP, NQTOT, NPCALL
      DOUBLE PRECISION INTEREAL,P(3*NATOMS)

      COMMON /TOT/ NQTOT
      COMMON /PCALL/ NPCALL
C
C  SAVE THE LOWEST NSAVEINTE DISTINGUISHABLE CONFIGURATIONS.
C
C     WRITE(*,'(A,12E15.7)') 'INTEREAL,ECONV,QMIN',INTEREAL,ECONV,(INTEQMIN(J1),J1=1,NSAVEINTE)
      DO J1=1,NSAVEINTE
         IF (DABS(INTEREAL-INTEQMIN(J1)).LT.ECONV) THEN
C
C  THESE ARE PROBABLY THE SAME - BUT JUST TO MAKE SURE WE SAVE THE 
C  LOWEST.
C
            IF (INTEREAL.LT.INTEQMIN(J1)) THEN
               INTEQMIN(J1)=INTEREAL
               DO J2=1,3*NATOMS
                  INTEQMINP(J1,J2)=P(J2)
               ENDDO
            ENDIF
            GOTO 10
         ENDIF
         IF (INTEREAL.LT.INTEQMIN(J1)) THEN

            J2=NSAVEINTE
20          CONTINUE
      
            IF (NSAVEINTE.GT.1) THEN
               INTEQMIN(J2)=INTEQMIN(J2-1)
               INTEFF(J2)=INTEFF(J2-1)
               DO J3=1,3*NATOMS
                  INTEQMINP(J2,J3)=INTEQMINP(J2-1,J3)
               ENDDO

               J2=J2-1
               IF (J2.GE.J1+1) GOTO 20
            ENDIF

            INTEQMIN(J1)=INTEREAL
            INTEFF(J1)=NQ(NP)
            DO J2=1,3*NATOMS
               INTEQMINP(J1,J2)=P(J2)
            ENDDO

            GOTO 10
         ENDIF
      ENDDO

10    CONTINUE

      RETURN

      END SUBROUTINE A9INTESAVEIT

