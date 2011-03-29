C   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF OPTIM.
C
C   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C
C***********************************************************************
C
      SUBROUTINE ADM(Q)
      USE PORFUNCS
      USE COMMONS
      USE KEY
      IMPLICIT NONE
      INTEGER IEXIT, JBOT, NTIMES, ICOUNT, I, ICN, J
      DOUBLE PRECISION DIS, TOL, DIST, Q(3*NATOMS)
C
C  PRINTS LOWER DIAGONAL ATOMIC DISTANCE MATRIX 
C
      PARAMETER (TOL = 0.5D0)

       WRITE(6,1092)
1092   FORMAT(' INTERATOMIC DISTANCE MATRIX ')
       JBOT=1
       NTIMES=1+(NATOMS-1)/5
       DO 20 ICOUNT=1,NTIMES
          WRITE(6,*)
          WRITE(6,142)(ZSYM(ICN),ICN=JBOT,MIN(NATOMS,JBOT+4))
142       FORMAT(18X,A3,4(10X,A3))
          WRITE(6,144)(ICN,ICN=JBOT,MIN(NATOMS,JBOT+4))
144       FORMAT(16X,:'[',I3,']',4(8X,:'[',I3,']'))
          DO 10 I=JBOT,NATOMS
             WRITE(6,143)ZSYM(I),I,(DIST(Q(3*I-2),Q(3*J-2)),J=JBOT
     &      ,MIN(I,JBOT+4))
143          FORMAT(T3,A3,'[',I3,']',5(3X,F10.5))
10     CONTINUE
       JBOT=JBOT+5
20     CONTINUE
C
C CHECK DISTANCES TO SEE IF ANY ARE TOO SHORT.  CALL EXIT IF BELOW
C  A CERTAIN TOLERANCE.
C
       IF (DCHECK) THEN
          DO I=1,NATOMS
             DO J=I+1,NATOMS
                DIS=DIST(Q(3*I-2),Q(3*J-2))
                IF(DIS.LT.TOL)THEN
C
C IF ONE IS A DUMMY, GO ON.
C
                    IF(ATMASS(I).LT.1.D-3.OR.ATMASS(J).LT.1.D-3) GOTO 12
                    WRITE(6,231)I,J,DIS,TOL
231                 FORMAT(T3,' ATOMS ',I2,' AND ',I2,
     1               ' ARE RATHER CLOSE.',
     &               /,T3,' DISTANCE OF ',F6.4,' IS BELOW THRESHOLD OF '
     &                        ,F6.4,'.')
                    IEXIT=0
                ENDIF
12              CONTINUE
             ENDDO
          ENDDO
       ENDIF
       IF(IEXIT .EQ. 1) THEN
          WRITE (6,*) ' *INSPECT DISTANCE MATRIX AND CORRECT ',
     $       'Z-MATRIX.'
          STOP
       ENDIF
      RETURN
      END

