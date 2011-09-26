
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
!
!  Calling CENTRE here without an evaporation check can put particles
!  outside the container, and make a valid step in takestep impossible.
!
!     PRINT*,'Calling centre from quench'
!     IF ((.NOT.FIELDT).AND.(.NOT.SEEDT).AND.CENT) CALL CENTRE2(COORDS(1:3*NATOMS,NP))

      IF (DUMPT) THEN
     !    IF (ARNO) THEN
            !WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS+2
            !WRITE(DUMPXYZUNIT+NP,70) NP,NQ(NP),EREAL,RMS
            !WRITE(DUMPXYZUNIT+NP,'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            !WRITE(DUMPXYZUNIT+NP,'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            !IF (NS.NE.0) WRITE(DUMPXYZUNIT+NP,65) (P(I),I=1,3*(NATOMS-NS))
!65          FORMAT('AR ',3F20.10)
         !ELSE IF (TIP) THEN
            !WRITE(DUMPVUNIT-NP,'(1X,F20.10,E20.10)') EREAL, POTEL
            !WRITE(DUMPXYZUNIT+NP,'(I6)') (NATOMS/2)*3
            !WRITE(DUMPXYZUNIT+NP,70) NP,NQ(NP), EREAL, RMS
            !DO J2=1,NATOMS/2
               !CALL TIPIO(P(3*(J2-1)+1),P(3*(J2-1)+2),P(3*(J2-1)+3),&
     !&              P(3*(NATOMS/2+J2-1)+1),P(3*(NATOMS/2+J2-1)+2),P(3*(NATOMS/2+J2-1)+3),RBCOORDS)
               !WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
               !WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
               !WRITE(DUMPXYZUNIT+NP,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            !ENDDO
         !ELSE IF (CHRMMT) THEN
            !CALL CHARMMDUMP3(P)
            !CALL CHARMMDUMP2(P,DUMPXYZUNIT+NP) ! xyz
         !ELSEIF (NCORE(NP).GT.0) THEN
         IF (NCORE(NP).GT.0) THEN
            WRITE(DUMPVUNIT-NP,'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,70) NQ(NP), EREAL, RMS
!           WRITE(DUMPXYZUNIT+NP,80) ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NCORE(NP))
!           WRITE(DUMPXYZUNIT+NP,80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NCORE(NP)+1,NATOMS)
            WRITE(DUMPXYZUNIT+NP,80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NCORE(NP))
            IF (NCORE(NP).GT.0) WRITE(DUMPXYZUNIT+NP,80) &
     &                     ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NCORE(NP)+1,NATOMS)
         ELSE
            WRITE(DUMPVUNIT-NP,'(1X,F20.10,E20.10)') EREAL, POTEL
            WRITE(DUMPXYZUNIT+NP,'(I4)') NATOMS
            WRITE(DUMPXYZUNIT+NP,70) NQ(NP), EREAL, RMS
70          FORMAT(1X,'QUENCH NUMBER ',I6,' final energy=',F20.10,' RMS force=',E20.10)
            WRITE(DUMPXYZUNIT+NP,80) ('LA ',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=1,NATOMS-NS)
            IF (NS.NE.0) WRITE(DUMPXYZUNIT+NP,80) ('LB',P(3*(I-1)+1),P(3*(I-1)+2),P(3*(I-1)+3),I=NATOMS-NS+1,NATOMS)
80          FORMAT(A2,3F20.10)
         ENDIF
      ENDIF

      IF (SQUEEZET) THEN
         IF ((EREAL.GT.0.0D0).AND.(SQUEEZED.LT.1.0D0)) THEN
            SQUEEZED=2.0D0-SQUEEZED
            NSQSTEPS=NQ(NP)
         ELSE
            NSQSTEPS=100000
         ENDIF
         DO J1=1,3*NVEC
            VEC(J1)=VEC(J1)*SQUEEZED
         ENDDO
         IF (NQ(NP).GT.2*NSQSTEPS) SQUEEZET=.FALSE.
      ENDIF
    
      IF ((NQ(NP).GE.NSSTOP).AND.SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
         WRITE(LFH,'(I6,A,G20.10)') NSSTOP,' quenches completed, setting coordinates to the lowest minimum, E=',QMIN(1)
         DO J1=1,3*NATOMS
            COORDS(J1,NP)=QMINP(1,J1)
         ENDDO
         POTEL=QMIN(1)
         EREAL=POTEL
      ENDIF

      RETURN
      ! }}}
      END SUBROUTINE QUENCH

