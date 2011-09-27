
      SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,BRUN,QDONE,P)
!op226> Declarations {{{ 
      use COMMONS
      USE QMODULE
      use porfuncs
      IMPLICIT NONE

      INTEGER I, J1, NSQSTEPS, NP, IFLAG, ITER, NOPT, J2, NDUMMY, J3, CSMIT
      DOUBLE PRECISION P(3*NATOMS),POTEL,TIME,EREAL,RBCOORDS(18),TMPCOORDS(3*NATOMS), DIST, QE, QX, AVVAL, CSMRMS
      LOGICAL QTEST, CFLAG, RES, COMPON, evapreject
      DOUBLE PRECISION  GRAD(3*NATOMS), DUMMY, DUM(3*NATOMS), DISTMIN, SSAVE, DIST2, RMAT(3,3)
C     DOUBLE PRECISION  WORK(60*NATOMS)
      DOUBLE PRECISION, PARAMETER :: HALFPI=1.570796327D0

      CHARACTER(LEN=80) DSTRING
      COMMON /MYPOT/ POTEL
      COMMON /CO/ COMPON
      COMMON /DMIN/ DISTMIN
      LOGICAL GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      common /ev/ evapreject
      DOUBLE PRECISION QSTART, QFINISH
      COMMON /Q4C/ QSTART, QFINISH
      COMMON /CSMAVVAL/ AVVAL, CSMRMS, CSMIT

      INTEGER BRUN, QDONE,ii

      NOPT=3*NATOMS
      IF (QTEST) THEN
         GMAX=CQMAX
      ELSE
         GMAX=BQMAX
      ENDIF

      QDONE=0
      DO I=1,3*NATOMS
         P(I)=COORDS(I,NP)
      ENDDO

      COMPON=.FALSE.
10    continue

      CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
      IF (EVAPREJECT) RETURN
      POTEL=EREAL

      IF (CFLAG) QDONE=1
      IF (.NOT.CFLAG) THEN
         IF (QTEST) THEN
            WRITE(MYUNIT,'(A,I6,A)') 'WARNING - Final Quench ',NQ(NP),'  did not converge'
         ELSE
            IF (NPAR.GT.1) THEN
               WRITE(MYUNIT,'(A,I7,A)') 'WARNING - Quench ',NQ(NP),'  did not converge'
            ELSE
               WRITE(MYUNIT,'(A,I7,A)') 'WARNING - Quench ',NQ(NP),'  did not converge'
            ENDIF
         ENDIF
      ENDIF

      CALL MYCPU_TIME(TIME)

      RES=.FALSE.

      IF (SAVEQ) CALL GSAVEIT(EREAL,P,NP)

      IF (.NOT.NORESET) THEN
         DO J1=1,3*(NATOMS-NSEED)
            COORDS(J1,NP)=P(J1)
         ENDDO
         DO J1=1,NATOMS
            VAT(J1,NP)=VT(J1)
         ENDDO
      ENDIF

      IF ((NQ(NP).GE.NSSTOP).AND.SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
         WRITE(MYUNIT,'(I6,A,G20.10)') NSSTOP,' quenches completed, setting coordinates to the lowest minimum, E=',QMIN(1)
         DO J1=1,3*NATOMS
            COORDS(J1,NP)=QMINP(1,J1)
         ENDDO
         POTEL=QMIN(1)
         EREAL=POTEL
      ENDIF

      RETURN
      END
