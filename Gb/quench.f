C  Conjugate gradient driver. 
C  CFLAG convergence test
C  CTEST checks for changes in chirality for AMBER runs
C
!> \param QTEST LOGICAL 
!> \param NP INTEGER
!> \param ITER INTEGER
!> \param TIME DOUBLE PRECISION
!> \param BRUN INTEGER
!> \param QDONE INTEGER
!> \param P DOUBLE PRECISION(3*NATOMS)
      SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,BRUN,QDONE,P)
!op226> Declarations {{{ 
      ! modules {{{
      USE COMMONS
      USE QMODULE
      USE PORFUNCS
      ! }}}
      IMPLICIT NONE
      ! subroutine {{{
      LOGICAL :: QTEST
      INTEGER ::    NP,ITER
      DOUBLE PRECISION :: TIME
      INTEGER ::    BRUN,QDONE
      DOUBLE PRECISION, DIMENSION(3*NATOMS) :: P
      ! }}}
      ! local {{{
      INTEGER I, J1, NSQSTEPS, IFLAG, NOPT, J2, NDUMMY, J3, CSMIT, J5
      DOUBLE PRECISION POTEL,EREAL,RBCOORDS(18),DIST, QE, QX, AVVAL, CSMRMS
      LOGICAL CFLAG, RES, COMPON, EVAPREJECT, PASS, FAIL
      DOUBLE PRECISION, DIMENSION(3*NATOMS) :: GRAD, DUM, TMPCOORDS
      DOUBLE PRECISION :: DUMMY, DISTMIN, SSAVE, DIST2, RMAT(3,3)
      DOUBLE PRECISION, PARAMETER :: HALFPI=1.570796327D0

      CHARACTER(LEN=80) DSTRING
      LOGICAL GUIDECHANGET, GUIDET, CSMDOGUIDET, DUMMYL
      ! }}}
      ! common {{{
      COMMON /MYPOT/ POTEL
      COMMON /CO/ COMPON
      COMMON /DMIN/ DISTMIN
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /EV/ EVAPREJECT
      COMMON /CSMAVVAL/ AVVAL, CSMRMS, CSMIT
      ! }}}
!op226>}}} 
      ! SUBROUTINE BODY {{{

C  TURN ON GUIDING POTENTIALS. THESE GET TURNED OFF IN POTENTIAL.F WHEN
C  THE RMS FORCE IS SMALL ENOUGH.
C
      NOPT=3*NATOMS

      IF (QTEST) THEN
         GMAX=FQMAX
      ELSE
         GMAX=SQMAX
      ENDIF

      QDONE=0

      DO I=1,3*NATOMS
         P(I)=COORDS(I)
      ENDDO

      COMPON=.FALSE.

      CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.)

      IF (EVAPREJECT) RETURN

      POTEL=EREAL

      IF (CFLAG) QDONE=1
      IF (.NOT.CFLAG) THEN
         IF (QTEST) THEN
            WRITE(LFH,'(A,I6,A)') 'WARNING - FINAL QUENCH ',NQ,'  DID NOT CONVERGE'
         ELSE
               WRITE(LFH,'(A,I7,A)') 'WARNING - QUENCH ',NQ,'  DID NOT CONVERGE'
         ENDIF
      ENDIF

      CALL MYCPU_TIME(TIME)

      RES=.FALSE.

C  IF EPSSPHERE IS NON-ZERO WE ARE PRESUMABLY DOING A CALCULATION OF THE 
C  ENERGY DENSITY OF LOCAL MINIMA. WE NEED TO KNOW THE MINIMUM DISTANCE
C  BETWEEN THE STARTING POINT AND THE QUENCHED MINIMA.
C
      IF ((EPSSPHERE.NE.0.0D0)) THEN
         DO J1=1,3*NATOMS
            DUM(J1)=COORDS(J1)
         ENDDO
C
C  DUM IS RETURNED IN THE CLOSEST ORIENTATION TO P; P SHOULD NOT CHANGE.
C  THIS IS NEARLY THE SAME MIND AS OPTIM! TO EXECUTE A RANDOM WALK WE MUST TAKE 
C  ANOTHER STEP AND MINIMISE UNTIL THE DISTANCE BETWEEN THE STARTING POINT
C  AND THE QUENCH MINIMUM IS LESS THAN EPSSPHERE.
C
         !CALL NEWMINDIST(P,DUM,NATOMS,DISTMIN,PERIODIC,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
      ENDIF
C
C  DEAL WITH EPSSPHERE SAMPLING.
C
      IF (EPSSPHERE.NE.0.0D0) THEN
         IF ((DISTMIN.GT.EPSSPHERE).OR.(ABS(EREAL-EPREV).LE.EDIFF)) THEN
            WRITE(LFH,'(A,F12.5,A,4F14.5)') 'STEP ',STEP,' EREAL, EPREV, DISTMIN, EPSSPHERE=',
     1                                     EREAL, EPREV, DISTMIN, EPSSPHERE
            DO J1=1,3*NATOMS
               COORDS(J1)=COORDSO(J1)
            ENDDO
            CALL TAKESTEP
             WRITE(LFH,'(A,G20.10)' ) 'RESEEDING STEP, MAXIMUM DISPLACEMENT RESET TO ',STEP
         ELSE
            WRITE(LFH,'(A,2F20.10)') 'VALID STEP, DISTMIN, EPSSPHERE=',DISTMIN, EPSSPHERE
         ENDIF
      ENDIF

      IF (.NOT.NORESET) THEN
         DO J1=1,3*(NATOMS-NSEED)
            COORDS(J1)=P(J1)
         ENDDO
         DO J1=1,NATOMS
            VAT(J1)=VT(J1)
         ENDDO
      ENDIF

C  CALLING CENTRE HERE WITHOUT AN EVAPORATION CHECK CAN PUT PARTICLES
C  OUTSIDE THE CONTAINER, AND MAKE A VALID STEP IN TAKESTEP IMPOSSIBLE.
C
C     PRINT*,'CALLING CENTRE FROM QUENCH'
C     IF ((.NOT.FIELDT).AND.(.NOT.SEEDT).AND.CENT) CALL CENTRE2(COORDS(1:3*NATOMS,NP))

      IF (DUMPT) THEN
        ! {{{
         IF (NCORE.GT.0) THEN
            WRITE(DUMPVUNIT,'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT,70) NQ, EREAL, RMS
            WRITE(DUMPXYZUNIT,80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NCORE)
            IF (NCORE.GT.0) WRITE(DUMPXYZUNIT,80) 
     &                     ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NCORE+1,NATOMS)
         ELSE
            WRITE(DUMPVUNIT,'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT,70) NQ, EREAL, RMS
70          FORMAT(1X,'QUENCH NUMBER ',I6,' FINAL ENERGY=',F20.10,' RMS FORCE=',E20.10)
            WRITE(DUMPXYZUNIT,80) ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NS)
            IF (NS.NE.0) WRITE(DUMPXYZUNIT,80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NS+1,NATOMS)
80          FORMAT(A2,3F20.10)
         ENDIF
         ! }}}
      ENDIF

      IF ((NQ.GE.NSSTOP).AND.SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
         WRITE(LFH,'(I6,A,G20.10)') NSSTOP,' QUENCHES COMPLETED, SETTING COORDINATES TO THE LOWEST MINIMUM, E=',QMIN(1)
         DO J1=1,3*NATOMS
            COORDS(J1)=QMINP(1,J1)
         ENDDO
         POTEL=QMIN(1)
         EREAL=POTEL
      ENDIF

      RETURN
      ! }}}
      END
