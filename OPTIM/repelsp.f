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
      SUBROUTINE REPELSP(COORDS,RMS,INEG,MFLAG)
      USE COMMONS
      USE KEY
      IMPLICIT NONE
      DOUBLE PRECISION COORDS(3*NATOMS), DUM(3*NATOMS), DIST, RMAT(3,3)
      DOUBLE PRECISION RMS, DUMMY
      INTEGER INEG
      LOGICAL NOSHIFTSAVE, NOITSAVE, MFLAG

C     DOUBLE PRECISION REPELPUSH
C     DOUBLE PRECISION REPELTS(3*NATOMS,100)
C     INTEGER NREPELTS, REPELFROM
C     LOGICAL REPELTST, REPEL
C     COMMON /OTS/ NREPELTS, REPELTST, REPELPUSH, REPEL, REPELFROM

      INTEGER J1, J2
      SAVE NOSHIFTSAVE, NOITSAVE

C     PRINT*,'NOIT=',NOIT
C     PRINT*,'REPEL,REPELFROM,RMS=',REPEL,REPELFROM,RMS
C     PRINT*,'INEG,NOSHIFT,NOIT=',INEG,NOSHIFT,NOIT

      IF (RMS.LT.1.0D-4) THEN
         IF (REPEL.AND.(INEG.EQ.2)) THEN
            DO J2=1,NOPT
               DUM(J2)=REPELTS(J2,REPELFROM)
            ENDDO
            CALL NEWMINDIST(COORDS,DUM,NATOMS,DIST,BULKT,TWOD,ZSYM(1),.FALSE.,RIGIDBODY,DEBUG,RMAT)
            DUMMY=0.0D0
            DO J1=1,NOPT
               DUM(J1)=COORDS(J1)-DUM(J1)
               DUMMY=DUMMY+DUM(J1)**2
            ENDDO
            DUMMY=1.0D0/SQRT(DUMMY)
            DO J1=1,NOPT
               COORDS(J1)=COORDS(J1)+REPELPUSH*DUM(J1)*DUMMY
            ENDDO
            REPEL=.FALSE.
            INDEXT=.FALSE.
            HINDEX=1
            NOIT=NOITSAVE
C           NOSHIFT=NOSHIFTSAVE
         ELSE IF ((.NOT.REPEL).AND.(INEG.EQ.1)) THEN
            DO J1=1,NREPELTS
               DO J2=1,NOPT
                  DUM(J2)=REPELTS(J2,J1)
               ENDDO
               CALL NEWMINDIST(COORDS,DUM,NATOMS,DIST,BULKT,TWOD,ZSYM(1),.FALSE.,RIGIDBODY,DEBUG,RMAT)
               IF (DIST.LT.0.05D0) THEN
                  WRITE(*,'(A,F15.5,A,I5)') ' SEARCH HAS MOVED TO A DISTANCE OF ',DIST,' FROM TRANSITION STATE ',J1
                  REPEL=.TRUE.
                  MFLAG=.FALSE.
                  REPELFROM=J1
                  INDEXT=.TRUE.
                  HINDEX=2
                  NOITSAVE=NOIT
                  NOIT=.TRUE.
C                 NOSHIFTSAVE=NOSHIFT
C
C     EIGENVALUE SHIFTING
C
C                 IF (NOSHIFT) THEN
C                    NOSHIFT=.FALSE.
                     DO J2=1,6
                        SHIFTL(J2)=0.0D0
                     ENDDO
                     IF (RTEST) THEN
                        IF (JZ.NE.0.0D0) THEN
                           WRITE(*,'(A,G20.10)') ' USING TRANS/Z ROTATION SHIFT=',SHIFTV
                           SHIFTL(3)=SHIFTV
                           SHIFTL(6)=SHIFTV
                        ELSE
                           WRITE(*,'(A,G20.10)') ' USING Z ROTATION SHIFT=',SHIFTV
                           SHIFTL(6)=SHIFTV
                        ENDIF
C                    ELSE IF (NOSHIFT) THEN
C                       WRITE(*,'(A)') 'NO SHIFTING'
                     ELSE IF (TWOD) THEN
                        IF (.NOT.BULKT) THEN
                           SHIFTL(1)=SHIFTV
                           SHIFTL(2)=SHIFTV
                           SHIFTL(6)=SHIFTV
                           WRITE(*,'(A)') ' X,Y TRANSLATIONAL AND Z ROTATIONAL SHIFTING'
                        ELSE
                           SHIFTL(1)=SHIFTV
                        SHIFTL(2)=SHIFTV
                           WRITE(*,'(A)') ' X,Y TRANSLATIONAL SHIFTING'
                        ENDIF
                     ELSE IF (EYTRAPT.OR.(ZSYM(NATOMS).EQ.'BE')) THEN
                        SHIFTL(4)=SHIFTV
                        SHIFTL(5)=SHIFTV
                        SHIFTL(6)=SHIFTV
                        WRITE(*,'(A,G20.10)') ' USING ROTATIONAL SHIFT=',SHIFTV
                     ELSE IF (BULKT) THEN
                        WRITE(*,'(A,G20.10)') ' USING TRANSLATIONAL SHIFT=',SHIFTV
                        DO J2=1,3
                           SHIFTL(J2)=SHIFTV
                        ENDDO
                     ELSE
                        WRITE(*,'(A,G20.10)') ' USING TRANSLATIONAL/ROTATIONAL SHIFT=',SHIFTV
                        DO J2=1,6
                           SHIFTL(J2)=SHIFTV
                        ENDDO
                     ENDIF
C                 ENDIF

                  RETURN
               ENDIF
C           DO J2=1,NOPT
C              REPELTS(J2,J1)=DUM(J2)
C           ENDDO
C           IF (DIST.LT.DMIN) THEN
C              DMIN=DIST
C              JMIN=J1
C           ENDIF
C           COMP=0.0D0
C           DO J2=1,NOPT
C              COMP=COMP+VEC(J2)*(COORDS(J2)-DUM(J2))/DIST
C           ENDDO
C           DO J2=1,NOPT
C              VEC(J2)=VEC(J2)-(1.0D0-AMOUNT)*COMP*(COORDS(J2)-DUM(J2))/DIST
C           ENDDO
C           WRITE(*,'(A,I6,F20.10)') 'J1,COMP=',J1,COMP
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END

