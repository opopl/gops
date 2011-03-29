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
      SUBROUTINE TWOEND(ENERGY,GRAD,VECS,Q)
      USE PORFUNCS
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE ZWK
      USE MODCHARMM
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, ITER, NS, J3, ITDONE
      DOUBLE PRECISION DIST, ENERGY, RMS, VECS(3*NATOMS), EVALMIN, EVALMAX, 
     1                 GRAD(3*NATOMS), SOVER, RMS2, EREAL, PROJP,
     2                 EP, STARTP(3*NATOMS), DELTAP, DELTA, RVEC(3*NATOMS), DUMMY1, PROJ, Q(3*NATOMS), VECL(3*NATOMS)
C
C  ASSIGN ENOUGH MEMORY TO WORK FOR A BLOCKSIZE OF 32 TO BE POSSIBLE.
C  THIS IS FOR DSYEVR.
C
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NATOMS)
      INTEGER IWORK(33*3*NATOMS), INFO
      DOUBLE PRECISION WORK(33*3*NATOMS),  ABSTOL, DIAG(3*NATOMS), ZMAT(3*NATOMS,3*NATOMS), DLAMCH
      LOGICAL PVFLAG, MFLAG, RESET, SWITCH,CONVERGED           !  IS A FULL ZMAT REALLY NEEDED?
      COMMON /PVF/ PVFLAG
C     COMMON /WORK/ ZMAT

      LWORK=33*3*NATOMS
      ILWORK=33*3*NATOMS
C
C  BULK PRINCIPAL IMAGE CODE ASSUMES THAT PARAM1, PARAM2 AND PARAM3 ARE THE BOX LENGTHS.
C
      OPEN(UNIT=87,FILE='CANDIDATE',STATUS='UNKNOWN')

      DO J1=1,NATOMS
         J3=3*(J1-1)
         VECS(J3+1)=FIN(J3+1)-Q(J3+1)
         VECS(J3+2)=FIN(J3+2)-Q(J3+2)
         VECS(J3+3)=FIN(J3+3)-Q(J3+3)
         LANV(J3+1)=0
         LANV(J3+2)=0
         LANV(J3+3)=0
         IF (BULKT) THEN
            LANV(J3+1)=NINT(VECS(J3+1)/PARAM1)
            LANV(J3+2)=NINT(VECS(J3+2)/PARAM2)
            LANV(J3+3)=NINT(VECS(J3+3)/PARAM3)
         ENDIF
         VECS(J3+1)=VECS(J3+1)-PARAM1*LANV(J3+1)
         VECS(J3+2)=VECS(J3+2)-PARAM2*LANV(J3+2)
         VECS(J3+3)=VECS(J3+3)-PARAM3*LANV(J3+3)
         START(J3+1)=Q(J3+1)
         START(J3+2)=Q(J3+2)
         START(J3+3)=Q(J3+3)
      ENDDO

C     DO J1=1,NOPT
C        LANV(J1)=0
C        IF (BULKT) LANV(J1)=NINT((FIN(J1)-Q(J1))/PARAM3)
C        VECS(J1)=FIN(J1)-Q(J1)-PARAM3*LANV(J1)
C        START(J1)=Q(J1)
C     ENDDO

      CALL ORTHOGOPT(VECS,Q,.TRUE.)
      ITER=1
C
C  INITIAL DISTANCE
C
      DIST=0.0D0
      DO J1=1,NOPT
         DIST=DIST+(FIN(J1)-START(J1)-PARAM3*LANV(J1))**2
C        WRITE(*,'(A,I6,2F20.10,I4,F20.10)') 'J1,FIN,START,LANV,DIST=',J1,FIN(J1),START(J1),LANV(J1),DIST
      ENDDO
      DIST=SQRT(DIST)
      EP=0.0D0
      PROJP=0.0D0
      DELTAP=0.0D0
      DO J1=1,NOPT
         STARTP(J1)=START(J1)
      ENDDO
C
C  MINIMISATION SUBJECT TO FORCE CONSTANT FORCE
C
10    SWITCH=.FALSE.
      DIST=0.0D0
      DO J1=1,NOPT
         DIST=DIST+(FIN(J1)-START(J1)-PARAM3*LANV(J1))**2
      ENDDO
      DIST=SQRT(DIST)

      NUP=1
      DO J1=1,NOPT
         ZWORK(J1,1)=ZMAT(J1,1)
      ENDDO
C
C  CALCULATE REQUIRED FORCE CONSTANT
C
      DUMMY1=0.0D0
      DO J1=1,NOPT
         RVEC(J1)=FIN(J1)-START(J1)
         DUMMY1=DUMMY1+RVEC(J1)**2
      ENDDO
      DUMMY1=1.0D0/SQRT(DUMMY1)
C     WRITE(*,'(A,F20.10)') 'DIST=',1.0D0/DUMMY1
      DO J1=1,NOPT
         RVEC(J1)=RVEC(J1)*DUMMY1
      ENDDO
      DUMMY1=0.0D0
      DO J1=1,NOPT
        DUMMY1=DUMMY1+RVEC(J1)*GRAD(J1)
      ENDDO
      PROJ=DUMMY1
C     PRINT*,'PROJECTION OF GRADIENT=',DUMMY1
C     FORCE=-DUMMY1-MAX(FINC*ABS(DUMMY1),FINC)
      FORCE=-DUMMY1-FINC-ABS(DUMMY1)*FINC
      IF (DUMMY1.LT.-FINC-ABS(DUMMY1)*FINC) FORCE=0.0D0
      RESET=.FALSE.
      IF (ITER.EQ.1) RESET=.TRUE.
      IF (UNRST.OR.CHRMMT) THEN
         PRINT*,'ERROR - TWOEND NOT COMPATIBLE WITH UNRES OR CHARMM'
         STOP
      ENDIF
      CALL MYLBFGS(NOPT,MUPDATE,START,.FALSE.,RMSTWO,MFLAG,ENERGY,RMS2,EREAL,RMS,NTWO,RESET,
     1             ITDONE,.FALSE.,GRAD,.FALSE.,.FALSE.)
      WRITE(23,'(3F20.10)') (START(J1),J1=1,NOPT)
      IF (NOHESS) THEN
         CALL BEIG(ITER,START,ENERGY,VECS,EVALMIN,NS,SOVER,.TRUE.,CONVERGED)
      ELSE 
         CALL POTENTIAL(START,ENERGY,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
         IF (NOIT) THEN
            IF (.NOT.VARIABLES) CALL SHIFTH(START,.TRUE.,NOPT,NATOMS,ATMASS)

            ABSTOL=DLAMCH('SAFE  MINIMUM')
            CALL DSYEVR('V','I','U',NOPT,HESS,SIZE(HESS,1),0.0D0,1.0D0,1,1,ABSTOL,NFOUND,DIAG,ZMAT,3*NATOMS,ISUPPZ,WORK,
     1                     LWORK, IWORK, ILWORK, INFO )
            IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' IN DSYEVR'
C           PRINT*,'OPTIMAL AND ACTUAL VALUES OF LWORK=',WORK(1),LWORK
C           PRINT*,'OPTIMAL AND ACTUAL VALUES OF ILWORK=',IWORK(1),ILWORK
            EVALMIN=DIAG(1)

C           WRITE(*,'(A,F20.10)') ' SMALLEST EIGENVALUE=',EVALMIN
            DO J1=1,NOPT
               VECS(J1)=ZMAT(J1,1)
            ENDDO
C           PRINT*,'EIGENVALUES:'
C           WRITE(*,'(6F15.5)') (DIAG(J1),J1=1,6)
         ELSE
            CALL ITEIG(ITER,START,VECS,EVALMIN,EVALMAX,NS,SOVER,.FALSE.,VECL,CONVERGED)
         ENDIF
      ENDIF
      DELTA=ENERGY-EP
      DELTA=PROJ-PROJP
      WRITE(*,'(I4,A,F15.7,A,F12.5,A,F15.5,A,F15.5,A,F12.5,A,F12.5,A,F12.5)') 
     1  ITER,' K=',FORCE,' D=',DIST,' E=',ENERGY,' PROJ=',PROJ,' DELTA=',DELTA,' EVALUE=',EVALMIN
      IF ((DELTA.LT.TWOEVAL).AND.(ITER.GT.2)) THEN
         SWITCH=.TRUE.
      ELSE
         DO J1=1,NOPT
            STARTP(J1)=START(J1)
         ENDDO
      ENDIF
      EP=ENERGY
      PROJP=PROJ
      DELTAP=DELTA
      WRITE(4,'(3F20.10)') (START(J1),J1=1,NOPT)
      IF (DIST.LT.0.05) THEN
         PRINT*,'FINAL GEOMETRY REACHED IN TWOENDS'
         CLOSE(87)
         STOP
         RETURN
      ENDIF
      IF (DUMPV) THEN
         IF (.NOT.ALLSTEPS) REWIND(44)
         WRITE(44,'(E20.10)') EVALMIN
         WRITE(44,'(3F20.10)') (VECS(J1),J1=1,NOPT)
      ENDIF
      IF (SWITCH) THEN
         WRITE(*,'(A)') ' CANDIDATE FOR TRANSITION STATE SEARCH LOCATED'
C        DO J1=1,NOPT
C           Q(J1)=START(J1)
C        ENDDO
C        TWOENDS=.FALSE.
         WRITE(87,'(3F20.10)') (START(J1),J1=1,NOPT)
C        RETURN
      ENDIF
      ITER=ITER+1
      IF (ITER.GT.NTWOITER) THEN
         PRINT*,'MAXIMUM STEPS EXCEEDED IN DOUBLE ENDED SEARCH'
         STOP
      ENDIF

      IF (BULKT) THEN
         DO J1=1,NATOMS
            J3=3*(J1-1)
            LANV(J3+1)=NINT((FIN(J3+1)-Q(J3+1))/PARAM1)
            LANV(J3+2)=NINT((FIN(J3+2)-Q(J3+2))/PARAM2)
            LANV(J3+3)=NINT((FIN(J3+3)-Q(J3+3))/PARAM3)
         ENDDO
      ENDIF
      GOTO 10

      CLOSE(87)
      RETURN
      END
