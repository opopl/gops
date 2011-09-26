
      ! Doxygen {{{
!  Conjugate gradient driver. 
!  CFLAG convergence test
!  CTEST checks for changes in chirality for AMBER runs
!
!> \param QTEST LOGICAL 
!> \param NP INTEGER
!> \param ITER INTEGER
!> \param TIME DOUBLE PRECISION
!> \param BRUN INTEGER
!> \param QDONE INTEGER
!> \param P DOUBLE PRECISION(3*NATOMS)
      ! }}}
      SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,BRUN,QDONE,P)
!op226> Declarations {{{ 
      ! mod {{{
      USE COMMONS
      USE F 
      USE V
      !USE MODAMBER9, ONLY : CISARRAY1, CISARRAY2, CHIARRAY1, CHIARRAY2
      !USE MODAMBER9, ONLY : SETCHIRAL, NOCISTRANSDNA, NOCISTRANSRNA
      !USE QMODULE
      USE PORFUNCS

      IMPLICIT NONE
      ! }}}
      ! sub {{{ 
      DOUBLE PRECISION, DIMENSION(:),INTENT(OUT) :: P   
      INTEGER NP,ITER,BRUN,QDONE
      LOGICAL QTEST
      DOUBLE PRECISION ::   TIME
      ! }}}
      ! loc {{{ 

      DOUBLE PRECISION ::   POTEL,EREAL
      DOUBLE PRECISION,DIMENSION(NR) ::GRAD,TMPCOORDS
      INTEGER ii
      INTEGER I, J1, NSQSTEPS, IFLAG, J2, NDUMMY, J3, CSMIT, J5
      DOUBLE PRECISION RBCOORDS(18), DIST, QE, QX, AVVAL, CSMRMS
      LOGICAL CFLAG, RES, COMPON, EVAPREJECT, PASS, FAIL
      DOUBLE PRECISION  DUMMY, DUM(3*NATOMS), DISTMIN, SSAVE, DIST2, RMAT(3,3)
      DOUBLE PRECISION, PARAMETER :: HALFPI=1.570796327D0
      CHARACTER(LEN=80) DSTRING
      LOGICAL GUIDECHANGET, GUIDET, CSMDOGUIDET
      LOGICAL  DUMMYL
!   sf344> gradually changing parameters to prevent dissociation of PY ellipsoids with repulsive sites 
      DOUBLE PRECISION epssave(3)
      ! }}}
      ! common {{{
      COMMON /MYPOT/ POTEL
      COMMON /CO/ COMPON
      COMMON /DMIN/ DISTMIN
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /EV/ EVAPREJECT
      ! }}}

!op226>}}} 
      ! subroutine body {{{

!  Turn on guiding potentials. These get turned off in potential.F when
!  the RMS force is small enough.

      SSAVE=STEP(NP)
11    continue

!  QTEST is set for the final quenches with tighter convergence criteria.

      IF (QTEST) THEN
         GMAX=FQMAX
      ELSE
         GMAX=SQMAX
      ENDIF

      QDONE=0
      P=COORDS(1:NR,NP)

10    CALL MYLBFGS(NR,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.,NP)

      IF (EVAPREJECT) RETURN
      POTEL=EREAL

      IF (CFLAG) QDONE=1
      IF (.NOT.CFLAG) THEN
         IF (QTEST) THEN
            WRITE(LFH,'(A,I6,A)') 'WARNING - Final Quench ',NQ(NP),'  did not converge'
         ELSE
            WRITE(LFH,'(A,I7,A)') 'WARNING - Quench ',NQ(NP),'  did not converge'
         ENDIF
      ENDIF

      CALL MYCPU_TIME(TIME)

      IF(SAVEQ)CALL GSAVEIT(EREAL,P,NP)

      !  Deal with EPSSPHERE sampling.!{{{

      IF (EPSSPHERE.NE.0.0D0) THEN
         IF ((DISTMIN.GT.EPSSPHERE).OR.(ABS(EREAL-EPREV(NP)).LE.ECONV)) THEN
            WRITE(LFH,'(A,F12.5,A,4F14.5)') 'step ',STEP(NP),' EREAL, EPREV, DISTMIN, EPSSPHERE=',&
     &                                     EREAL, EPREV(NP), DISTMIN, EPSSPHERE
            DO J1=1,3*NATOMS
               COORDS(J1,NP)=COORDSO(J1,NP)
            ENDDO
            CALL TAKESTEP(NP)
             WRITE(LFH,'(A,G20.10)' ) 'reseeding step, maximum displacement reset to ',STEP(NP)
            GOTO 11
         ELSE
            WRITE(LFH,'(A,2F20.10)') 'valid step, DISTMIN, EPSSPHERE=',DISTMIN, EPSSPHERE
         ENDIF
      ENDIF!}}}

!  NORESET true does not set the configuration point to the quench geometry
!  A relaxed frozen core does not get saved, but the lowest minima are saved
!  by GSAVEIT.
!
      IF (.NOT.NORESET) THEN
         DO J1=1,3*(NATOMS-NSEED)
            COORDS(J1,NP)=P(J1)
         ENDDO
         DO J1=1,NATOMS
            VAT(J1,NP)=VT(J1)
         ENDDO
      ENDIF

      RETURN
      ! }}}
      END SUBROUTINE QUENCH

