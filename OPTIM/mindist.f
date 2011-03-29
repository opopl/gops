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
C  FINDS THE MINIMUM DISTANCE BETWEEN TWO GEOMETRIES USING LBFGS
C  GEOMETRY R0 IS RESET TO R1 AFTER EACH OPTIMISATION STEP. GEOMETRY
C  IN R SHOULD NOT CHANGE, SO NEITHER SHOULD RA. RB IS RETURNED AS THE
C  CLOSEST GEOMETRY TO RA.
C
C JMC AS LONG AS ZSYM ISN'T 'W' (IN WHICH CASE MIND DOES SOMETHING SPECIAL) MIND
C DOESN'T CARE WHAT ATOMIC SYMBOL WE GIVE IT
C
      SUBROUTINE MINDIST(RA,RB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET)
      USE KEY,ONLY : FREEZE, PULLT, EFIELDT
      IMPLICIT NONE
      INTEGER J1, IG, I, ITER, J2, NCOUNT
      DOUBLE PRECISION DPRAND
      DOUBLE PRECISION P(3), DIST, DIST0, RG, RG0, DISTFUNC, 
     1                 MYROTMAT(3,3), 
     2                 OVEC(3), H1VEC(3), H2VEC(3), DSAVE, OMEGATOT(3,3), RA(*), RB(*)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: R, R0, R1, R1SAVE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RBSAVE
      INTEGER NSIZE, NATOMS
      LOGICAL BULKT, MFLAG, TWOD, AGAIN, PRESERVET
      COMMON /GEOM/ NSIZE
      COMMON /MINDOM/ MYROTMAT, OMEGATOT
      CHARACTER(LEN=5) ZUSE

      ALLOCATE(RBSAVE(3*NATOMS))
      RBSAVE(1:3*NATOMS)=RB(1:3*NATOMS)
C
C  INITIALISE ACCUMULATED ROTATION MATRIX TO THE IDENTITY.
C
      DO J1=1,3
         DO J2=1,3
            OMEGATOT(J2,J1)=0.0D0
         ENDDO
         OMEGATOT(J1,J1)=1.0D0
      ENDDO

      AGAIN=.FALSE.
      DSAVE=0.0D0
      IF (ZUSE(1:1).EQ.'W') THEN
         ALLOCATE(R(3,3*(NATOMS/2)),R0(3,3*(NATOMS/2)),R1(3,3*(NATOMS/2)),R1SAVE(3,3*(NATOMS/2)))
C        NSIZE=3*NATOMS
         NSIZE=3*(NATOMS/2)
C        DO J1=1,NATOMS
         DO J1=1,NATOMS/2
            CALL CONVERT(RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
C    1              RA(3*(NATOMS+J1-1)+1),RA(3*(NATOMS+J1-1)+2),RA(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
     1              RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3),OVEC,H1VEC,H2VEC)
            R(1,(J1-1)*3+1)=OVEC(1)
            R(2,(J1-1)*3+1)=OVEC(2)
            R(3,(J1-1)*3+1)=OVEC(3)
            R(1,(J1-1)*3+2)=H1VEC(1)
            R(2,(J1-1)*3+2)=H1VEC(2)
            R(3,(J1-1)*3+2)=H1VEC(3)
            R(1,(J1-1)*3+3)=H2VEC(1)
            R(2,(J1-1)*3+3)=H2VEC(2)
            R(3,(J1-1)*3+3)=H2VEC(3)
            CALL CONVERT(RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3),
C    1            RB(3*(NATOMS+J1-1)+1),RB(3*(NATOMS+J1-1)+2),RB(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
     1            RB(3*(NATOMS/2+J1-1)+1),RB(3*(NATOMS/2+J1-1)+2),RB(3*(NATOMS/2+J1-1)+3),OVEC,H1VEC,H2VEC)
            R0(1,(J1-1)*3+1)=OVEC(1)
            R0(2,(J1-1)*3+1)=OVEC(2)
            R0(3,(J1-1)*3+1)=OVEC(3)
            R0(1,(J1-1)*3+2)=H1VEC(1)
            R0(2,(J1-1)*3+2)=H1VEC(2)
            R0(3,(J1-1)*3+2)=H1VEC(3)
            R0(1,(J1-1)*3+3)=H2VEC(1)
            R0(2,(J1-1)*3+3)=H2VEC(2)
            R0(3,(J1-1)*3+3)=H2VEC(3)
         ENDDO
      ELSE
         ALLOCATE(R(3,NATOMS),R0(3,NATOMS),R1(3,NATOMS),R1SAVE(3,NATOMS))
         NSIZE=NATOMS
         DO J1=1,NSIZE
            R(1,J1)=RA(3*(J1-1)+1)
            R(2,J1)=RA(3*(J1-1)+2)
            R(3,J1)=RA(3*(J1-1)+3)
            R0(1,J1)=RB(3*(J1-1)+1)
            R0(2,J1)=RB(3*(J1-1)+2)
            R0(3,J1)=RB(3*(J1-1)+3)
         ENDDO
      ENDIF
C
C     PUT C.O.M. TO ORIGIN
C
      IF (.NOT.FREEZE) THEN
         DO IG=1,3
            RG=0.D0
            RG0=0.D0
            DO I=1,NSIZE
               RG=RG+R(IG,I)
               RG0=RG0+R0(IG,I)
            ENDDO
            RG=RG/NSIZE
            RG0=RG0/NSIZE
            DO I=1,NSIZE
               R(IG,I)=R(IG,I)-RG
               R0(IG,I)=R0(IG,I)-RG0
            ENDDO
         ENDDO
      ENDIF
C
C     INITIAL ANGLES
C
      P(1)=0.0D0
      P(2)=0.0D0
      P(3)=0.0D0
C     IF (TWOD) P(1)=0.0D0
C     IF (TWOD) P(2)=0.0D0
C
C     CALCULATE INITIAL DISTANCE
C
      NCOUNT=0
10    DIST0=DISTFUNC(P,R,R0,R1)
      DIST0=SQRT(DIST0)
      IF (BULKT) THEN
         DIST=DIST0
         IF (PRESERVET) RB(1:3*NATOMS)=RBSAVE(1:3*NATOMS)
         DEALLOCATE(R, R0, R1, R1SAVE, RBSAVE)
         RETURN
      ENDIF

      CALL MMYLBFGS(P,1.0D-7,MFLAG,DIST,500,ITER,R,R0,R1)
      DIST=DISTFUNC(P,R,R0,R1)
      DIST=SQRT(DIST)
      IF (MFLAG) THEN
C        WRITE(*,'(A,2F15.5,A,I6)') 'INITIAL AND FINAL DISTANCES:',DIST0,DIST,' ITERATIONS=',ITER
C        PRINT*
      ELSE
         NCOUNT=NCOUNT+1
         IF (NCOUNT.GT.0) THEN 
            PRINT*,'CONVERGENCE FAILURE IN MIND'
C           STOP
         ELSE
C           WRITE(*,'(A,2F15.5,A,I6,A,I6)') 'INITIAL AND FINAL DISTANCES:',DIST0,DIST,' ITERATIONS=',ITER,' NCOUNT=',NCOUNT
            DO J1=1,NSIZE
               R0(1,J1)=R1(1,J1)
               R0(2,J1)=R1(2,J1)
               R0(3,J1)=R1(3,J1)
            ENDDO
            IF (.NOT.(TWOD.OR.PULLT.OR.EFIELDT)) P(1)=2*(DPRAND()-0.5D0)/10.0D0
            IF (.NOT.(TWOD.OR.PULLT.OR.EFIELDT)) P(2)=2*(DPRAND()-0.5D0)/10.0D0
            P(3)=2*(DPRAND()-0.5D0)/10.0D0
            GOTO 10
         ENDIF
      ENDIF
C
C  THIS BLOCK ALLOWS THE SECOND GEOMETRY TO ROTATE OUT OF THE XY PLANE; NOT ALLOWED!
C
C     IF (TWOD) THEN
C        IF (AGAIN) THEN
C           WRITE(*,'(A,2F20.10)') 'DIST,DSAVE=',DIST,DSAVE
C           IF (DIST.GT.DSAVE) THEN
C              DIST=DSAVE
C              DO J1=1,NSIZE
C                 R1(1,J1)=R1SAVE(1,J1)
C                 R1(2,J1)=R1SAVE(2,J1)
C                 R1(3,J1)=R1SAVE(3,J1)
C                 R0(1,J1)=-R0(1,J1)
C              ENDDO
C              DO J1=1,3
C                 DO J2=1,3
C                    OSAVE(J2,J1)=OSAVESAVE(J2,J1)
C                 ENDDO
C              ENDDO
C           ENDIF
C        ELSE
C           AGAIN=.TRUE.
C           DSAVE=DIST
C           DO J1=1,NSIZE
C              R0(1,J1)=-R0(1,J1)
C              R1SAVE(1,J1)=R1(1,J1)
C              R1SAVE(2,J1)=R1(2,J1)
C              R1SAVE(3,J1)=R1(3,J1)
C           ENDDO
C           DO J1=1,3
C              DO J2=1,3
C                 OSAVE(J2,J1)=0.0D0
C                 OSAVESAVE(J2,J1)=OSAVE(J2,J1)
C              ENDDO
C              OSAVE(J1,J1)=1.0D0
C           ENDDO
C           GOTO 10
C        ENDIF
C     ENDIF

      IF (ZUSE(1:1).EQ.'W') THEN
C        DO J1=1,NATOMS  !  WCOMMENT
         DO J1=1,NATOMS/2
            OVEC(1)=R(1,(J1-1)*3+1)
            OVEC(2)=R(2,(J1-1)*3+1)
            OVEC(3)=R(3,(J1-1)*3+1)
            H1VEC(1)=R(1,(J1-1)*3+2)
            H1VEC(2)=R(2,(J1-1)*3+2)
            H1VEC(3)=R(3,(J1-1)*3+2)
            H2VEC(1)=R(1,(J1-1)*3+3)
            H2VEC(2)=R(2,(J1-1)*3+3)
            H2VEC(3)=R(3,(J1-1)*3+3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,OVEC =',J1,OVEC(1),OVEC(2),OVEC(3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,H1VEC=',J1,H1VEC(1),H1VEC(2),H1VEC(3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,H2VEC=',J1,H2VEC(1),H2VEC(2),H2VEC(3)
            CALL CONVERT2(OVEC,H1VEC,H2VEC,RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
     1                 RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3))
C    1                 RA(3*(NATOMS+J1-1)+1),RA(3*(NATOMS+J1-1)+2),RA(3*(NATOMS+J1-1)+3)) ! WCOMMENT
C           WRITE(*,'(A,I6,6F15.5)') 'J1,RA=',J1,RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
C    1              RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3)
            OVEC(1)=R1(1,(J1-1)*3+1)
            OVEC(2)=R1(2,(J1-1)*3+1)
            OVEC(3)=R1(3,(J1-1)*3+1)
            H1VEC(1)=R1(1,(J1-1)*3+2)
            H1VEC(2)=R1(2,(J1-1)*3+2)
            H1VEC(3)=R1(3,(J1-1)*3+2)
            H2VEC(1)=R1(1,(J1-1)*3+3)
            H2VEC(2)=R1(2,(J1-1)*3+3)
            H2VEC(3)=R1(3,(J1-1)*3+3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,OVEC =',J1,OVEC(1),OVEC(2),OVEC(3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,H1VEC=',J1,H1VEC(1),H1VEC(2),H1VEC(3)
C           WRITE(*,'(A,I6,3F15.5)') 'J1,H2VEC=',J1,H2VEC(1),H2VEC(2),H2VEC(3)
            CALL CONVERT2(OVEC,H1VEC,H2VEC,RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3),
C    1                 RB(3*(NATOMS+J1-1)+1),RB(3*(NATOMS+J1-1)+2),RB(3*(NATOMS+J1-1)+3))
     1                 RB(3*(NATOMS/2+J1-1)+1),RB(3*(NATOMS/2+J1-1)+2),RB(3*(NATOMS/2+J1-1)+3))
C           WRITE(*,'(A,I6,6F15.5)') 'J1,RA=',J1,RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
C    1              RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3)
         ENDDO
      ELSE
         DO J1=1,NSIZE
            RB(3*(J1-1)+1)=R1(1,J1)
            RB(3*(J1-1)+2)=R1(2,J1)
            RB(3*(J1-1)+3)=R1(3,J1)
            RA(3*(J1-1)+1)=R(1,J1)
            RA(3*(J1-1)+2)=R(2,J1)
            RA(3*(J1-1)+3)=R(3,J1)
         ENDDO
      ENDIF

      DO J1=1,3
         DO J2=1,3
            MYROTMAT(J2,J1)=OMEGATOT(J2,J1)
         ENDDO
      ENDDO

      IF (PRESERVET) RB(1:3*NATOMS)=RBSAVE(1:3*NATOMS)

      DEALLOCATE(R, R0, R1, R1SAVE, RBSAVE)

      RETURN
      END

C__________________________________________________________________________

        SUBROUTINE INVERSE(A,B)
        IMPLICIT NONE
        DOUBLE PRECISION A(3,3),B(3,3),C(2,2),D(3,3),DET
        INTEGER I, J, JJ, II

        DO I=1,3
           DO J=1,3
              DO JJ=1,J-1
                 DO II=1,I-1
                    C(II,JJ)=A(II,JJ)
                 ENDDO
                 DO II=I+1,3
                    C(II-1,JJ)=A(II,JJ)
                 ENDDO
              ENDDO
              
              DO JJ=J+1,3
                 DO II=1,I-1
                    C(II,JJ-1)=A(II,JJ)
                 ENDDO
                 DO II=I+1,3
                    C(II-1,JJ-1)=A(II,JJ)
                    
                 ENDDO
              ENDDO
              
              B(I,J)=C(1,1)*C(2,2)-C(1,2)*C(2,1)
              D(I,J)=B(I,J)*(-1)**(I+J)
           ENDDO
        ENDDO
        
        DET=A(1,1)*D(1,1)+A(1,2)*D(1,2)+A(1,3)*D(1,3)
        IF (DET.EQ.0.0D0) WRITE(*,'(A,G20.10)') 'ERROR, DETERMINANT IN MIND ROUTINE INVERSE=',DET
        
        DO I=1,3
           DO J=1,3
              B(I,J)=D(J,I)/DET
           ENDDO
        ENDDO
        
        RETURN
        END
C_____________________________________________________________________

      SUBROUTINE PROD(A,B,C)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3),B(3,3),C(3,3)
      INTEGER I, J

      DO I=1,3
         DO J=1,3
            C(I,J)=A(I,1)*B(1,J)+A(I,2)*B(2,J)+A(I,3)*B(3,J)
         ENDDO
      ENDDO

      RETURN
      END
C_____________________________________________________________________

      FUNCTION DISTFUNC(P,R,R0,R1)
      IMPLICIT NONE

      DOUBLE PRECISION R(3,*),R0(3,*),A(3,3),R1(3,*),DIST,DISTFUNC
      DOUBLE PRECISION MYROTMAT(3,3),ROTX(3,3),ROTY(3,3),ROTZ(3,3),OMEGATOT(3,3)
      DOUBLE PRECISION P(3)
      INTEGER NSIZE,I,J
      COMMON /MINDOM/ MYROTMAT, OMEGATOT

      COMMON/GEOM/NSIZE

      ROTX(1,1)=1.0D0
      ROTX(1,2)=0.0D0
      ROTX(2,1)=0.0D0
      ROTX(1,3)=0.0D0
      ROTX(3,1)=0.0D0
      ROTX(2,2)=DCOS(P(1))
      ROTX(3,3)=ROTX(2,2)
      ROTX(2,3)=DSIN(P(1))
      ROTX(3,2)=-ROTX(2,3)

      ROTY(1,2)=0.0D0
      ROTY(2,1)=0.0D0
      ROTY(2,3)=0.0D0
      ROTY(3,2)=0.0D0
      ROTY(2,2)=1.0D0
      ROTY(1,1)=DCOS(P(2))
      ROTY(3,3)=ROTY(1,1)
      ROTY(1,3)=-DSIN(P(2))
      ROTY(3,1)=-ROTY(1,3)

      ROTZ(1,3)=0.0D0
      ROTZ(3,1)=0.0D0
      ROTZ(2,3)=0.0D0
      ROTZ(3,2)=0.0D0
      ROTZ(3,3)=1.0D0
      ROTZ(1,1)=DCOS(P(3))
      ROTZ(2,2)=ROTZ(1,1)
      ROTZ(1,2)=DSIN(P(3))
      ROTZ(2,1)=-ROTZ(1,2)

      A(1,1)=ROTY(1,1)*ROTZ(1,1)+ROTY(1,2)*ROTZ(2,1)+ROTY(1,3)*ROTZ(3,1)
      A(1,2)=ROTY(1,1)*ROTZ(1,2)+ROTY(1,2)*ROTZ(2,2)+ROTY(1,3)*ROTZ(3,2)
      A(1,3)=ROTY(1,1)*ROTZ(1,3)+ROTY(1,2)*ROTZ(2,3)+ROTY(1,3)*ROTZ(3,3)
      A(2,1)=ROTY(2,1)*ROTZ(1,1)+ROTY(2,2)*ROTZ(2,1)+ROTY(2,3)*ROTZ(3,1)
      A(2,2)=ROTY(2,1)*ROTZ(1,2)+ROTY(2,2)*ROTZ(2,2)+ROTY(2,3)*ROTZ(3,2)
      A(2,3)=ROTY(2,1)*ROTZ(1,3)+ROTY(2,2)*ROTZ(2,3)+ROTY(2,3)*ROTZ(3,3)
      A(3,1)=ROTY(3,1)*ROTZ(1,1)+ROTY(3,2)*ROTZ(2,1)+ROTY(3,3)*ROTZ(3,1)
      A(3,2)=ROTY(3,1)*ROTZ(1,2)+ROTY(3,2)*ROTZ(2,2)+ROTY(3,3)*ROTZ(3,2)
      A(3,3)=ROTY(3,1)*ROTZ(1,3)+ROTY(3,2)*ROTZ(2,3)+ROTY(3,3)*ROTZ(3,3)

      MYROTMAT(1,1)=ROTX(1,1)*A(1,1)+ROTX(1,2)*A(2,1)+ROTX(1,3)*A(3,1)
      MYROTMAT(1,2)=ROTX(1,1)*A(1,2)+ROTX(1,2)*A(2,2)+ROTX(1,3)*A(3,2)
      MYROTMAT(1,3)=ROTX(1,1)*A(1,3)+ROTX(1,2)*A(2,3)+ROTX(1,3)*A(3,3)
      MYROTMAT(2,1)=ROTX(2,1)*A(1,1)+ROTX(2,2)*A(2,1)+ROTX(2,3)*A(3,1)
      MYROTMAT(2,2)=ROTX(2,1)*A(1,2)+ROTX(2,2)*A(2,2)+ROTX(2,3)*A(3,2)
      MYROTMAT(2,3)=ROTX(2,1)*A(1,3)+ROTX(2,2)*A(2,3)+ROTX(2,3)*A(3,3)
      MYROTMAT(3,1)=ROTX(3,1)*A(1,1)+ROTX(3,2)*A(2,1)+ROTX(3,3)*A(3,1)
      MYROTMAT(3,2)=ROTX(3,1)*A(1,2)+ROTX(3,2)*A(2,2)+ROTX(3,3)*A(3,2)
      MYROTMAT(3,3)=ROTX(3,1)*A(1,3)+ROTX(3,2)*A(2,3)+ROTX(3,3)*A(3,3)

C     CALL PROD(ROTY,ROTZ,A)
C     CALL PROD(ROTX,A,MYROTMAT)

      DO I=1,NSIZE
         DO J=1,3
            R1(J,I)=MYROTMAT(J,1)*R0(1,I)+MYROTMAT(J,2)*R0(2,I)+MYROTMAT(J,3)*R0(3,I)
         ENDDO
      ENDDO

      DIST=0.D0
      DO I=1,NSIZE
         DIST=DIST+(R(1,I)-R1(1,I))**2+(R(2,I)-R1(2,I))**2+(R(3,I)-R1(3,I))**2
C        WRITE(*,'(A,I5,2G20.10)') 'I,DIST,TOTAL=',I,(R(1,I)-R1(1,I))**2+(R(2,I)-R1(2,I))**2+(R(3,I)-R1(3,I))**2,DIST
      ENDDO

      DISTFUNC=DIST
      
      RETURN
      END
C_____________________________________________________________________

      SUBROUTINE DDISTFUNC(P,F,R,R0,R1)
      IMPLICIT NONE

      DOUBLE PRECISION R(3,*),R0(3,*),R1(3,*)
      DOUBLE PRECISION P(3),F(3)
      INTEGER NSIZE, I

      COMMON/GEOM/NSIZE

      F(1)=0.D0
      F(2)=0.D0
      F(3)=0.D0
      DO I=1,NSIZE
         F(1)=F(1)+R(2,I)*R1(3,I)-R(3,I)*R1(2,I)
         F(2)=F(2)+R(3,I)*R1(1,I)-R(1,I)*R1(3,I)
         F(3)=F(3)+R(1,I)*R1(2,I)-R(2,I)*R1(1,I)
      ENDDO
      F(1)=-2.D0*F(1)
      F(2)=-2.D0*F(2)
      F(3)=-2.D0*F(3)

      RETURN

      END
C_____________________________________________________________________

      SUBROUTINE DDDISTFUNC(P,H,R,R0,R1)
      IMPLICIT NONE

      DOUBLE PRECISION R(3,*),R0(3,*),R1(3,*)
      DOUBLE PRECISION P(3),H(3,3)
      INTEGER NSIZE, I, J

      COMMON/GEOM/NSIZE

      DO I=1,3
         DO J=1,3
            H(I,J)=0.D0
            H(I,J)=0.D0
            H(I,J)=0.D0
         ENDDO
      ENDDO
      DO I=1,NSIZE
         H(1,1)=H(1,1)+R(2,I)*R1(2,I)+R(3,I)*R1(3,I)
         H(2,2)=H(2,2)+R(1,I)*R1(1,I)+R(3,I)*R1(3,I)
         H(3,3)=H(3,3)+R(1,I)*R1(1,I)+R(2,I)*R1(2,I)
         H(2,1)=H(2,1)-R(2,I)*R1(1,I)
         H(3,1)=H(3,1)-R(3,I)*R1(1,I)
         H(3,2)=H(3,2)-R(3,I)*R1(2,I)
      ENDDO
      H(1,2)=H(2,1)
      H(1,3)=H(3,1)
      H(2,3)=H(3,2)

      DO I=1,3
         DO J=1,3
            H(I,J)=2.D0*H(I,J)
         ENDDO
      ENDDO

      RETURN

      END

      SUBROUTINE ROTGEOM(NATOMS,COORDS)
      IMPLICIT NONE
      INTEGER I, J, K, NATOMS
      DOUBLE PRECISION COORDS(*), R1, R0(3), MYROTMAT(3,3),OMEGATOT(3,3)
      COMMON /MINDOM/ MYROTMAT, OMEGATOT

      DO I=1,NATOMS
         R0(1)=COORDS(3*(I-1)+1)
         R0(2)=COORDS(3*(I-1)+2)
         R0(3)=COORDS(3*(I-1)+3)
         DO J=1,3
            R1=0.0D0
            DO K=1,3
               R1=R1+MYROTMAT(J,K)*R0(K)
            ENDDO
            COORDS(3*(I-1)+J)=R1
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE CONVERT(X,Y,Z,A,B,C,OVEC,H1VEC,H2VEC)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,A,B,C,OVEC(3),H1VEC(3),H2VEC(3),PI,OL,TH,PH,PS,SINA,SINB,SINC,COSA,COSB,COSC,SP12,SP13,SP22,
     1                 SP23,SP32,SP33,SY12,SY22,SY32,SZ13,SZ23,SZ33,ANGLE,HOLEN,RANGLE,OHZ,HYL,HZL
      DATA ANGLE,HOLEN/104.52D0,0.9572D0/

      PI=4.0D0*ATAN(1.0D0)
      RANGLE=PI*ANGLE/360.0D0
      OHZ=HOLEN*COS(RANGLE)
      HYL=HOLEN*SIN(RANGLE)
      HZL=16.0D0*OHZ/18.0D0
      OL=-OHZ+HZL
C
CCCCCC    A,B,C ARE EULER ANGLES
C
      TH=A
      PH=B
      PS=C
      SINA=SIN(TH)
      SINB=SIN(PH)
      SINC=SIN(PS)
      COSA=COS(TH)
      COSB=COS(PH)
      COSC=COS(PS)
      SP12=-(SINC*COSB+COSA*SINB*COSC)
      SP13=SINA*SINB
      SP22=-SINC*SINB+COSA*COSB*COSC
      SP23=-SINA*COSB
      SP32=SINA*COSC
      SP33=COSA
      SY12=SP12*HYL
      SY22=SP22*HYL
      SY32=SP32*HYL
      SZ13=SP13*HZL
      SZ23=SP23*HZL
      SZ33=SP33*HZL
C
CCCCC HYDROGEN POSITIONS
C
      H1VEC(1)=SY12+SZ13+X
      H1VEC(2)=SY22+SZ23+Y
      H1VEC(3)=SY32+SZ33+Z

      H2VEC(1)=-SY12+SZ13+X
      H2VEC(2)=-SY22+SZ23+Y
      H2VEC(3)=-SY32+SZ33+Z
C
CCCC OXYGEN POSITION
C
      OVEC(1)=SP13*OL   +X
      OVEC(2)=SP23*OL   +Y
      OVEC(3)=SP33*OL   +Z

      RETURN
      END

      SUBROUTINE CONVERT2(OVEC,H1VEC,H2VEC,X,Y,Z,A,B,C)
C
C  CONVERT H, H, O POSITIONS TO EULERS - NEEDS H, H, O
C  POSITIONS AND CENTRES OF MASS. HERE WE ALLOW FOR THE POSSIBILITY THAT THE
C  IDEAL RIGID WATER GEOMETRY MAY BE BROKEN AND TAKE THE BEST FIT.
C
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,Z,A,B,C,OVEC(3),H1VEC(3),H2VEC(3),PI,OL,SINA,SINB,SINC,COSA,COSB,COSC,SP12,SP13,SP22,
     1                 SP23,SP32,SP33,SY12,SY22,SY32,SZ13,SZ23,SZ33,ANGLE,HOLEN,RANGLE,OHZ,HYL,HZL,SP13T,PI2,
     2                 SP23T,SP33T,ASAVE,BSAVE,CSAVE,SY12T,SY22T,SY32T,SZ13T,SZ23T,SZ33T,DMIN,ABEST,BBEST,CBEST
      INTEGER JA,JB,JC
      DATA ANGLE,HOLEN/104.52D0,0.9572D0/

      PI=4.0D0*ATAN(1.0D0) 
      PI2=2.0D0*PI
      RANGLE=PI*ANGLE/360.0D0 
      OHZ=HOLEN*COS(RANGLE) 
      HYL=HOLEN*SIN(RANGLE) 
      HZL=16.0D0*OHZ/18.0D0 
      OL=-OHZ+HZL 
C
      X=(H1VEC(1)+H2VEC(1)+16.0D0*OVEC(1))/18.0D0
      Y=(H1VEC(2)+H2VEC(2)+16.0D0*OVEC(2))/18.0D0
      Z=(H1VEC(3)+H2VEC(3)+16.0D0*OVEC(3))/18.0D0


      COSA=((OVEC(3)-Z)/OL)
      IF (ABS(COSA).GT.1.0D0) THEN
         COSA=ABS(COSA)/COSA
         A=DACOS(-COSA)
      ELSE
         A=DACOS((OVEC(3)-Z)/OL)
      ENDIF
      SINA=DSIN(A)
      IF (SINA.EQ.0.0D0) SINA=1.0D-10

      COSB=-(OVEC(2)-Y)/(OL*SINA)
      IF (ABS(COSB).GT.1.0D0) THEN
         COSB=ABS(COSB)/COSB
         B=DACOS(-COSB)
      ELSE
         B=DACOS(-(OVEC(2)-Y)/(OL*SINA))
      ENDIF
      SINB=DSIN(B)

      COSC=(H1VEC(3)-H2VEC(3))/(2.0D0*HYL*SINA)
      IF (ABS(COSC).GT.1.0D0) THEN
         COSC=ABS(COSC)/COSC
         C=DACOS(-COSC)
      ELSE
         C=DACOS((H1VEC(3)-H2VEC(3))/(2.0D0*HYL*SINA))
      ENDIF
      SINC=DSIN(C)

C     RETURN

      SP13T=(OVEC(1)-X)/OL
      SP23T=(OVEC(2)-Y)/OL
      SP33T=(OVEC(3)-Z)/OL

      SY12T=(H1VEC(1)-H2VEC(1))/(2.0D0)
      SY22T=(H1VEC(2)-H2VEC(2))/(2.0D0)
      SY32T=(H1VEC(3)-H2VEC(3))/(2.0D0)

      SZ13T=(H1VEC(1)+H2VEC(1)-2.0D0*X)/(2.0D0)
      SZ23T=(H1VEC(2)+H2VEC(2)-2.0D0*Y)/(2.0D0)
      SZ33T=(H1VEC(3)+H2VEC(3)-2.0D0*Z)/(2.0D0)
 
      ASAVE=A
      BSAVE=B
      CSAVE=C

      DMIN=1.0D100
      DO JA=1,4
         IF (JA.EQ.1) A=ASAVE
         IF (JA.EQ.2) A=PI2-ASAVE
         IF (JA.EQ.3) A=PI-ASAVE
         IF (JA.EQ.4) A=PI+ASAVE
         SINA=SIN(A)
         COSA=COS(A)
         DO JB=1,4
            IF (JB.EQ.1) B=BSAVE
            IF (JB.EQ.2) B=PI2-BSAVE
            IF (JB.EQ.3) B=PI-BSAVE
            IF (JB.EQ.4) B=PI+BSAVE
            SINB=SIN(B)
            COSB=COS(B)
            DO JC=1,4
               IF (JC.EQ.1) C=CSAVE
               IF (JC.EQ.2) C=PI2-CSAVE
               IF (JC.EQ.3) C=PI-CSAVE
               IF (JC.EQ.4) C=PI+CSAVE
               SINC=SIN(C)
               COSC=COS(C)

               SP12=-(SINC*COSB+COSA*SINB*COSC)
               SP13=SINA*SINB
               SP22=-SINC*SINB+COSA*COSB*COSC
               SP23=-SINA*COSB
               SP32=SINA*COSC
               SP33=COSA
               SY12=SP12*HYL
               SY22=SP22*HYL
               SY32=SP32*HYL
               SZ13=SP13*HZL
               SZ23=SP23*HZL
               SZ33=SP33*HZL

               IF ( DABS(SP13T-SP13)+
     1              DABS(SP23T-SP23)+
     1              DABS(SP33T-SP33)+
     1              DABS(SY12T-SY12)+
     1              DABS(SY22T-SY22)+
     1              DABS(SY32T-SY32)+
     1              DABS(SZ13T-SZ13)+
     1              DABS(SZ23T-SZ23)+
     1              DABS(SZ33T-SZ33) .LT. DMIN) THEN
                   DMIN= DABS(SP13T-SP13)+
     1              DABS(SP23T-SP23)+
     1              DABS(SP33T-SP33)+
     1              DABS(SY12T-SY12)+
     1              DABS(SY22T-SY22)+
     1              DABS(SY32T-SY32)+
     1              DABS(SZ13T-SZ13)+
     1              DABS(SZ23T-SZ23)+
     1              DABS(SZ33T-SZ33) 
                  ABEST=A
                  BBEST=B
                  CBEST=C
                  IF (DMIN.LT.1.0D-10) GOTO 20
               ENDIF
            ENDDO
         ENDDO
      ENDDO

20    A=ABEST
      B=BBEST
      C=CBEST
C     IF (DMIN.GT.0.1D0) WRITE(*,'(A,F15.5)') 'WARNING, DEVIATION FROM RIGID BODY GEOMETRY DETECTED, BEST FIT IS ',DMIN

      RETURN
      END

C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** JULY 1990 ***
C
C        LINE SEARCH REMOVED PLUS SMALL MODIFICATIONS, DJW 2001
C
      SUBROUTINE MMYLBFGS(X,EPS,MFLAG,ENERGY,ITMAX,ITDONE,R,R0,R1)
      USE KEY
      USE PORFUNCS
      IMPLICIT NONE
      INTEGER N,M,J1,J2,ITMAX,ITDONE,NFAIL,J
      PARAMETER (N=3,M=5)  !  MMUPDATE IS ACTUALLY IGNORED
      DOUBLE PRECISION X(N),G(3),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT,OVERLAP,DISTFUNC
      DOUBLE PRECISION EPS,DUMMY1,ENERGY,ENEW,RMS,ALPHA,GSAVE(3)
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,MYROTMAT(3,3),OMEGATOT(3,3),OTEMP(3,3)
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE
      LOGICAL MFLAG, PTEST
      INTEGER NSIZE, ISTAT
      DOUBLE PRECISION R(3,NSIZE), R0(3,NSIZE), R1(3,NSIZE), DOT1, DOT2
      COMMON /GEOM/ NSIZE
      COMMON /MINDOM/ MYROTMAT, OMEGATOT
C
C  SGI APPEARS TO NEED THIS SAVE STATEMENT!
C
      SAVE

      PTEST=.TRUE.
      PTEST=.FALSE.
      ALPHA=1.0D0
      NFAIL=0
      ITER=0
      ITDONE=0
      ENERGY=DISTFUNC(X,R,R0,R1)
      CALL DDISTFUNC(X,GSAVE,R,R0,R1)
      IF (TWOD.OR.PULLT.OR.EFIELDT) GSAVE(1)=0.0D0
      IF (TWOD.OR.PULLT.OR.EFIELDT) GSAVE(2)=0.0D0
      RMS=SQRT((GSAVE(1)**2+GSAVE(2)**2+GSAVE(3)**2)/3)
      DO J1=1,N
         G(J1)=GSAVE(J1)
      ENDDO

      IF (PTEST) WRITE(*,'(A,2F20.10,A,I6,A)') ' DISTANCE AND RMS FORCE=',ENERGY,RMS,' AFTER ',ITDONE,' LBFGS STEPS'

10    CALL FLUSH(6,ISTAT)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
         MFLAG=.TRUE.
         IF (MFLAG) THEN
C           WRITE(*,'(A,F20.10)') ' DIAGONAL INVERSE HESSIAN ELEMENTS ARE NOW ',DIAG(1)
            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
C        WRITE(*,'(A,F20.10)') ' DIAGONAL INVERSE HESSIAN ELEMENTS ARE NOW ',DIAG(1)
         RETURN
      ENDIF

      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(*,240)
 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         DO I=1,N
            DIAG(I)=1.0D0
         ENDDO
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
         ISPT= N+2*M
         IYPT= ISPT+N*M
C
C  NR STEP FOR DIAGONAL INVERSE HESSIAN
C
         DO I=1,N
            W(ISPT+I)= -G(I)*DIAG(I)
            W(I)= -G(I)*DIAG(I)
         ENDDO
         GNORM= DSQRT(DDOT(N,G,1,G,1))
C
C  MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS.
C
         STP=MIN(1.0D0/GNORM,GNORM)
C        STP=1.0D0
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
         IF (YS.EQ.0.0D0) YS=1.0D0
C
C  UPDATE ESTIMATE OF DIAGONAL INVERSE HESSIAN ELEMENTS
C  WE DIVIDE BY BOTH YS AND YY AT DIFFERENT POINTS, SO
C  THEY HAD BETTER NOT BE ZERO!
C
         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
         IF (YY.EQ.0.0D0) YY=1.0D0
         DUMMY1=YS/YY
C        DUMMY1=ABS(YS/YY)
         DO I=1,N
            DIAG(I)= DUMMY1
         ENDDO
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: NOCEDAL, J. 1980,
C     "UPDATING QUASI-NEWTON MATRICES WITH LIMITED STORAGE",
C     MATHEMATICS OF COMPUTATION, VOL.24, NO.151, PP. 773-782.
C     ---------------------------------------------------------
C
         CP= POINT
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS
         DO I=1,N
            W(I)= -G(I)
         ENDDO
         CP= POINT
         DO I= 1,BOUND
            CP=CP-1
            IF (CP.EQ. -1)CP=M-1
            SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)= W(N+CP+1)*SQ
            CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
         ENDDO
        
         DO I=1,N
            W(I)=DIAG(I)*W(I)
         ENDDO

         DO I=1,BOUND
            YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
            BETA= W(N+CP+1)*YR
            INMC=N+M+CP+1
            BETA= W(INMC)-BETA
            ISCN=ISPT+CP*N
            CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
            CP=CP+1
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0
      ENDIF
C
C  STORE THE NEW SEARCH DIRECTION
C
      IF (ITER.GT.0) THEN
         DO I=1,N
            W(ISPT+POINT*N+I)= W(I)
         ENDDO
      ENDIF
      IF (TWOD.OR.PULLT.OR.EFIELDT) W(ISPT+POINT*N+1)=0.0D0
      IF (TWOD.OR.PULLT.OR.EFIELDT) W(ISPT+POINT*N+2)=0.0D0

C     OVERLAP=DDOT(N,G,1,W,1)/SQRT(DDOT(N,G,1,G,1)*DDOT(N,W,1,W,1))
      DOT1=SQRT(DDOT(N,G,1,G,1))
      DOT2=SQRT(DDOT(N,W,1,W,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

C     PRINT*,'OVERLAP,DIAG(1)=',OVERLAP,DIAG(1)
C     PRINT*,'G . G=',DDOT(N,G,1,G,1)
C     PRINT*,'W . W=',DDOT(N,W,1,W,1)
      IF (OVERLAP.GT.0.0D0) THEN
         IF (PTEST) PRINT*,'SEARCH DIRECTION HAS POSITIVE PROJECTION ONTO GRADIENT - REVERSING STEP'
C        IF (PTEST) PRINT*,'SEARCH DIRECTION HAS POSITIVE PROJECTION ONTO GRADIENT - RESET'
         DO I=1,N
C           W(I)=-W(I)
            W(ISPT+POINT*N+I)= -W(I)
         ENDDO
C        ITER=0
C        GOTO 10
      ENDIF

      DO I=1,N
         W(I)=G(I)
      ENDDO
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
      IF (STP*SLENGTH.GT.MAXMBFGS) STP=MAXMBFGS/SLENGTH
C
C  WE NOW HAVE THE PROPOSED STEP.
C
      DO J1=1,N
         X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
      ENDDO 
      NDECREASE=0
20    ENEW=DISTFUNC(X,R,R0,R1)

      IF (ENEW-ENERGY.LE.1.0D-3) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
      ELSE 
C
C  ENERGY INCREASED - TRY AGAIN WITH A SMALLER STEP SIZE
C
C        IF (STP*SLENGTH.LT.1.0D-10) THEN
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            PRINT*,' IN MMYLBFGS LBFGS STEP CANNOT FIND A LOWER DISTANCE, NFAIL=',NFAIL
            ITER=0  !  TRY RESETTING
            IF (NFAIL.GT.20) THEN
               PRINT*,' TOO MANY FAILURES - GIVE UP'
               RETURN
            ENDIF
            DO J1=1,N
               X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
            ENDDO 
            GOTO 30
         ENDIF
         DO J1=1,N
            X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
         ENDDO 
         NDECREASE=NDECREASE+1
         STP=STP/10.0D0
         IF (PTEST) 
     1    WRITE(*,'(A,F19.10,A,F16.10,A,F15.8)') ' DISTANCE INCREASED FROM ',ENERGY,' TO ',ENEW,' DECREASING STEP TO ',STP*SLENGTH
         GOTO 20
      ENDIF
C
C  RESET COMPARISON GEOMETRY TO R1 AND ROTATION ANGLES TO ZERO. ACCUMULATE THE OVERALL TRANSFORMATION MATRIX IN 
C  OMEGATOT.
C
      OTEMP(1,1)=MYROTMAT(1,1)*OMEGATOT(1,1)+MYROTMAT(1,2)*OMEGATOT(2,1)+MYROTMAT(1,3)*OMEGATOT(3,1)
      OTEMP(1,2)=MYROTMAT(1,1)*OMEGATOT(1,2)+MYROTMAT(1,2)*OMEGATOT(2,2)+MYROTMAT(1,3)*OMEGATOT(3,2)
      OTEMP(1,3)=MYROTMAT(1,1)*OMEGATOT(1,3)+MYROTMAT(1,2)*OMEGATOT(2,3)+MYROTMAT(1,3)*OMEGATOT(3,3)
      OTEMP(2,1)=MYROTMAT(2,1)*OMEGATOT(1,1)+MYROTMAT(2,2)*OMEGATOT(2,1)+MYROTMAT(2,3)*OMEGATOT(3,1)
      OTEMP(2,2)=MYROTMAT(2,1)*OMEGATOT(1,2)+MYROTMAT(2,2)*OMEGATOT(2,2)+MYROTMAT(2,3)*OMEGATOT(3,2)
      OTEMP(2,3)=MYROTMAT(2,1)*OMEGATOT(1,3)+MYROTMAT(2,2)*OMEGATOT(2,3)+MYROTMAT(2,3)*OMEGATOT(3,3)
      OTEMP(3,1)=MYROTMAT(3,1)*OMEGATOT(1,1)+MYROTMAT(3,2)*OMEGATOT(2,1)+MYROTMAT(3,3)*OMEGATOT(3,1)
      OTEMP(3,2)=MYROTMAT(3,1)*OMEGATOT(1,2)+MYROTMAT(3,2)*OMEGATOT(2,2)+MYROTMAT(3,3)*OMEGATOT(3,2)
      OTEMP(3,3)=MYROTMAT(3,1)*OMEGATOT(1,3)+MYROTMAT(3,2)*OMEGATOT(2,3)+MYROTMAT(3,3)*OMEGATOT(3,3)
      DO J1=1,3
         DO J2=1,3
            OMEGATOT(J2,J1)=OTEMP(J2,J1)
         ENDDO
      ENDDO
      X(1)=0.0D0
      X(2)=0.0D0
      X(3)=0.0D0
      DO I=1,NSIZE
         DO J=1,3
            R0(J,I)=R1(J,I)
         ENDDO
      ENDDO

      CALL DDISTFUNC(X,GSAVE,R,R0,R1)
      IF (TWOD.OR.PULLT.OR.EFIELDT) GSAVE(1)=0.0D0
      IF (TWOD.OR.PULLT.OR.EFIELDT) GSAVE(2)=0.0D0
      RMS=SQRT((GSAVE(1)**2+GSAVE(2)**2+GSAVE(3)**2)/3)
      DO J1=1,N
         G(J1)=GSAVE(J1)
      ENDDO
      IF (PTEST) WRITE(*,'(A,2F20.10,A,I6,A,G15.5)') ' DISTANCE AND RMS FORCE=',ENERGY,RMS,' AFTER ',ITDONE,
     1        ' LBFGS STEPS, STEP:',STP*SLENGTH
C
C     COMPUTE THE NEW STEP AND GRADIENT CHANGE
C
30    NPT=POINT*N
      DO I=1,N
         W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
         W(IYPT+NPT+I)= G(I)-W(I)
      ENDDO
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      GOTO 10

      RETURN
      END
