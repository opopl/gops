
      SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,BRUN,QDONE,P)
!op226> Declarations {{{ 

      USE COMMONS
      USE QMODULE
      USE PORFUNCS

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

C
C   sf344> gradually changing parameters to prevent dissociation of PY ellipsoids with repulsive sites 
C
      DOUBLE PRECISION epssave(3)

C
C  Data for the screen saver.
C
      INTEGER BRUN, QDONE,ii
!op226>}}} 
C
C  Turn on guiding potentials. These get turned off in potential.F when
C  the RMS force is small enough.
C
      SSAVE=STEP(NP)
C
C csw34 Reset the NFIX counter
C
      NFIX=0

C  QTEST is set for the final quenches with tighter convergence criteria.
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

      CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)
      POTEL=EREAL

      IF (CFLAG) QDONE=1

      IF (.NOT.CFLAG) THEN
         IF (QTEST) THEN
            WRITE(MYUNIT,'(A,I6,A)') 'WARNING - Final Quench ',NQ(NP),'  did not converge'
         ELSE
         ENDIF
      ENDIF

      CALL MYCPU_TIME(TIME)

      RES=.FALSE.


      IF (.NOT.NORESET) THEN
         DO J1=1,3*(NATOMS-NSEED)
            COORDS(J1,NP)=P(J1)
         ENDDO
         DO J1=1,NATOMS
            VAT(J1,NP)=VT(J1)
         ENDDO
      ENDIF

      RETURN

      END
