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
      SUBROUTINE SHIFTH(Q,UNDO,NOPT,NATOMS,ATMASS)
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, NOPT, NATOMS, NIADD
      DOUBLE PRECISION EV(3*NATOMS,6), RNORM(6), DUMMY, CMX, CMY, CMZ, TMASS,
     1                 AMASS(NATOMS), RMASS(NATOMS), Q(3*NATOMS), ATMASS(NATOMS)

      LOGICAL UNDO

C     PRINT*,'IN SHIFTH SHIFTL=',SHIFTL
C
C  IF MASST THEN THE COORDINATES, Q, HAVE ALREADY BEEN MULTIPLIED BY THE SQUARE ROOT OF THE MASS
C
      SHIFTED=.TRUE.

      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      TMASS=0.0D0
      NIADD=0
      IF (RIGIDBODY) THEN
         NATOMS=NATOMS/2
         NIADD=3*NATOMS
      ENDIF
C     PRINT*,'NATOMS,NIADD=',NATOMS,NIADD
      DO J2=1,NATOMS
         AMASS(J2)=1.0D0
         RMASS(J2)=1.0D0
         IF (MASST) AMASS(J2)=ATMASS(J2)
         IF (MASST) RMASS(J2)=SQRT(ATMASS(J2))
         TMASS=TMASS+AMASS(J2)
      ENDDO
      DO J2=1,NATOMS
         CMX=CMX+Q(3*(J2-1)+1)*RMASS(J2)
         CMY=CMY+Q(3*(J2-1)+2)*RMASS(J2)
         CMZ=CMZ+Q(3*(J2-1)+3)*RMASS(J2)
      ENDDO
      CMX=CMX/TMASS
      CMY=CMY/TMASS
      CMZ=CMZ/TMASS
C     PRINT*,'CMX,CMY,CMZ, NATOMS=',CMX,CMY,CMZ,NATOMS

      DO J1=1,6
         RNORM(J1)=0.0D0
         DO J2=1,3*NATOMS+NIADD
            EV(J2,J1)=0.0D0
         ENDDO
      ENDDO

      IF (MASST) THEN

         DO J1=1,NATOMS
            J3=3*(J1-1)
            DUMMY=SQRT(ATMASS(J1))

            EV(J3+1,1)=DUMMY
            RNORM(1)=RNORM(1)+DUMMY

            EV(J3+2,2)=DUMMY
            RNORM(2)=RNORM(2)+DUMMY

            EV(J3+3,3)=DUMMY
            RNORM(3)=RNORM(3)+DUMMY
   
            EV(J3+2,4)= Q(J3+3)*DUMMY
            EV(J3+3,4)=-Q(J3+2)*DUMMY
            RNORM(4)=RNORM(4)+(Q(J3+3)**2+Q(J3+2)**2)*DUMMY**2
     
            EV(J3+1,5)=-Q(J3+3)*DUMMY
            EV(J3+3,5)= Q(J3+1)*DUMMY
            RNORM(5)=RNORM(5)+(Q(J3+3)**2+Q(J3+1)**2)*DUMMY**2
    
            EV(J3+1,6)= Q(J3+2)*DUMMY
            EV(J3+2,6)=-Q(J3+1)*DUMMY
            RNORM(6)=RNORM(6)+(Q(J3+2)**2+Q(J3+1)**2)*DUMMY**2
         ENDDO

      ELSE

         DO J1=1,NATOMS
            J3=3*(J1-1)

            EV(J3+1,1)=1.0D0
            RNORM(1)=RNORM(1)+1.0D0

            EV(J3+2,2)=1.0D0
            RNORM(2)=RNORM(2)+1.0D0

            EV(J3+3,3)=1.0D0
            RNORM(3)=RNORM(3)+1.0D0
  
            EV(J3+2,4)= Q(J3+3)-CMZ
            EV(J3+3,4)=-Q(J3+2)+CMY
            RNORM(4)=RNORM(4)+(Q(J3+3)-CMZ)**2+(Q(J3+2)-CMY)**2

C           EV(J3+2,4)= Q(J3+3)
C           EV(J3+3,4)=-Q(J3+2)
C           RNORM(4)=RNORM(4)+Q(J3+3)**2+Q(J3+2)**2

            EV(J3+1,5)=-Q(J3+3)+CMZ
            EV(J3+3,5)= Q(J3+1)-CMX
            RNORM(5)=RNORM(5)+(Q(J3+3)-CMZ)**2+(Q(J3+1)-CMX)**2

C           EV(J3+1,5)=-Q(J3+3)
C           EV(J3+3,5)= Q(J3+1)
C           RNORM(5)=RNORM(5)+Q(J3+3)**2+Q(J3+1)**2

            EV(J3+1,6)= Q(J3+2)-CMY
            EV(J3+2,6)=-Q(J3+1)+CMX
            RNORM(6)=RNORM(6)+(Q(J3+2)-CMY)**2+(Q(J3+1)-CMX)**2

C           EV(J3+1,6)= Q(J3+2)
C           EV(J3+2,6)=-Q(J3+1)
C           RNORM(6)=RNORM(6)+Q(J3+2)**2+Q(J3+1)**2
         ENDDO

      ENDIF

      DO J1=1,6
         RNORM(J1)=1.0D0/SQRT(RNORM(J1))
      ENDDO

      DO J1=1,NATOMS
         J3=3*(J1-1)

         EV(J3+1,1)=EV(J3+1,1)*RNORM(1)
         EV(J3+2,2)=EV(J3+2,2)*RNORM(2)
         EV(J3+3,3)=EV(J3+3,3)*RNORM(3)
 
         EV(J3+2,4)=EV(J3+2,4)*RNORM(4)
         EV(J3+3,4)=EV(J3+3,4)*RNORM(4)

         EV(J3+1,5)=EV(J3+1,5)*RNORM(5)
         EV(J3+3,5)=EV(J3+3,5)*RNORM(5)

         EV(J3+1,6)=EV(J3+1,6)*RNORM(6)
         EV(J3+2,6)=EV(J3+2,6)*RNORM(6)
      ENDDO

C     PRINT*,'UNSHIFTED HESS:'
C     WRITE(*,'(6F15.5)') ((HESS(J1,J2),J1=1,NOPT),J2=1,NOPT)

      IF (FREEZE) THEN
         DO J1=1,NATOMS+NIADD/3
            IF (.NOT.FROZEN(J1))  CYCLE
            J3=3*(J1-1)
            HESS(J3+1,J3+1)=HESS(J3+1,J3+1)+SHIFTL(1)
            HESS(J3+2,J3+2)=HESS(J3+2,J3+2)+SHIFTL(1)
            HESS(J3+3,J3+3)=HESS(J3+3,J3+3)+SHIFTL(1)
         ENDDO
      ELSEIF (.NOT.UNDO) THEN
         DO J1=1,NATOMS+NIADD/3
            J3=3*(J1-1)
            DO J2=1,NATOMS+NIADD/3
               J4=3*(J2-1)
   
               HESS(J4+1,J3+1)=HESS(J4+1,J3+1)-SHIFTL(1)*EV(J4+1,1)*EV(J3+1,1)
     1                                        -SHIFTL(5)*EV(J4+1,5)*EV(J3+1,5)-SHIFTL(6)*EV(J4+1,6)*EV(J3+1,6)
               HESS(J4+1,J3+2)=HESS(J4+1,J3+2)-SHIFTL(6)*EV(J4+1,6)*EV(J3+2,6)
               HESS(J4+1,J3+3)=HESS(J4+1,J3+3)-SHIFTL(5)*EV(J4+1,5)*EV(J3+3,5)
               HESS(J4+2,J3+1)=HESS(J4+2,J3+1)-SHIFTL(6)*EV(J4+2,6)*EV(J3+1,6)
               HESS(J4+2,J3+2)=HESS(J4+2,J3+2)-SHIFTL(2)*EV(J4+2,2)*EV(J3+2,2)
     1                                        -SHIFTL(4)*EV(J4+2,4)*EV(J3+2,4)-SHIFTL(6)*EV(J4+2,6)*EV(J3+2,6)
               HESS(J4+2,J3+3)=HESS(J4+2,J3+3)-SHIFTL(4)*EV(J4+2,4)*EV(J3+3,4)
               HESS(J4+3,J3+1)=HESS(J4+3,J3+1)-SHIFTL(5)*EV(J4+3,5)*EV(J3+1,5)
               HESS(J4+3,J3+2)=HESS(J4+3,J3+2)-SHIFTL(4)*EV(J4+3,4)*EV(J3+2,4)
               HESS(J4+3,J3+3)=HESS(J4+3,J3+3)-SHIFTL(3)*EV(J4+3,3)*EV(J3+3,3)
     1                                        -SHIFTL(4)*EV(J4+3,4)*EV(J3+3,4)-SHIFTL(5)*EV(J4+3,5)*EV(J3+3,5)

            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS+NIADD/3
            J3=3*(J1-1)
            DO J2=1,NATOMS+NIADD/3
               J4=3*(J2-1)
               HESS(J4+1,J3+1)=HESS(J4+1,J3+1)+SHIFTL(1)*EV(J4+1,1)*EV(J3+1,1)
     1                                        +SHIFTL(5)*EV(J4+1,5)*EV(J3+1,5)+SHIFTL(6)*EV(J4+1,6)*EV(J3+1,6)
               HESS(J4+1,J3+2)=HESS(J4+1,J3+2)+SHIFTL(6)*EV(J4+1,6)*EV(J3+2,6)
               HESS(J4+1,J3+3)=HESS(J4+1,J3+3)+SHIFTL(5)*EV(J4+1,5)*EV(J3+3,5)
               HESS(J4+2,J3+1)=HESS(J4+2,J3+1)+SHIFTL(6)*EV(J4+2,6)*EV(J3+1,6)
               HESS(J4+2,J3+2)=HESS(J4+2,J3+2)+SHIFTL(2)*EV(J4+2,2)*EV(J3+2,2)
     1                                        +SHIFTL(4)*EV(J4+2,4)*EV(J3+2,4)+SHIFTL(6)*EV(J4+2,6)*EV(J3+2,6)
               HESS(J4+2,J3+3)=HESS(J4+2,J3+3)+SHIFTL(4)*EV(J4+2,4)*EV(J3+3,4)
               HESS(J4+3,J3+1)=HESS(J4+3,J3+1)+SHIFTL(5)*EV(J4+3,5)*EV(J3+1,5)
               HESS(J4+3,J3+2)=HESS(J4+3,J3+2)+SHIFTL(4)*EV(J4+3,4)*EV(J3+2,4)
               HESS(J4+3,J3+3)=HESS(J4+3,J3+3)+SHIFTL(3)*EV(J4+3,3)*EV(J3+3,3)
     1                                        +SHIFTL(4)*EV(J4+3,4)*EV(J3+3,4)+SHIFTL(5)*EV(J4+3,5)*EV(J3+3,5)

            ENDDO
         ENDDO
      ENDIF

      IF (TWOD) THEN
         DO J2=3,3*NATOMS,3
            HESS(J2,J2)=SHIFTV
         ENDDO
      ENDIF

C     PRINT*,'SHIFTED HESS:'
C     WRITE(*,'(6F15.5)') ((HESS(J1,J2),J1=1,NOPT),J2=1,NOPT)
C     PRINT*,'COORDINATES:'
C     WRITE(*,'(6F15.5)') (Q(J1),J1=1,NOPT)

      IF (RIGIDBODY) NATOMS=NATOMS*2

      RETURN
      END
