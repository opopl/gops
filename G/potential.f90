
      SUBROUTINE POTENTIAL(X,GRAD,EREAL,GRADT,SECT)
!op226> Declarations {{{ 
      ! modules {{{
      USE COMMONS
      USE V
      USE F
      USE MODBLN
      USE PORFUNCS
      !}}}
      IMPLICIT NONE
      ! subroutine {{{
      !DOUBLE PRECISION, DIMENSION(:),INTENT(IN) :: X
      !DOUBLE PRECISION, DIMENSION(:),INTENT(OUT) :: GRAD
      DOUBLE PRECISION, DIMENSION(:),INTENT(IN) :: X
      DOUBLE PRECISION, DIMENSION(:)  :: GRAD
      DOUBLE PRECISION, INTENT(OUT) :: EREAL
      LOGICAL, INTENT(IN) :: GRADT, SECT       
      ! }}}
      ! local {{{
      INTEGER BRUN
      INTEGER NQTOT
      DOUBLE PRECISION :: EA(10)

      LOGICAL GUIDECHANGET, GUIDET, CSMDOGUIDET
      LOGICAL FTEST, EVAP, COMPON, YESNO, evapreject
      INTEGER J1, J2, J3, NPCALL, PERM(NATOMS), NPERM, NORBIT1, NORBIT2, CSMIT
      CHARACTER FNAME*80, DUMM*4
      DOUBLE PRECISION DUMMY2, GEMAX, XG, YG, ZG, RMAT(3,3), XD, YD, ZD, XTEMP(3*NATOMS)
      DOUBLE PRECISION :: GRADLJ(3*NATOMS), EREALLJ, GRADMF(3*NATOMS)
      DOUBLE PRECISION :: EREALMF, TERMLJ, TERMMF, GRADDUM1(3*NATOMS), AVVAL
      DOUBLE PRECISION :: GRADDUM2(3*NATOMS), EDUM1, EDUM2 
      DOUBLE PRECISION :: DUMMY(3*NATOMS), DIST2, WORSTRAD, GBDUMMY, QE
      DOUBLE PRECISION :: QX, DISTANCE, AA(3)
      DOUBLE PRECISION :: SAVECSMNORM, CSMRMS, CSMGRAD(3)
      DOUBLE PRECISION :: SAVECSMPMAT(3,3), SAVECSMIMAGES(3*NATOMS*CSMGPINDEX)
      CHARACTER(LEN=3) CSMGPSAVE
      INTEGER CSMGPINDEXSAVE
      DOUBLE PRECISION PTGPSAVE(3,3,2*CSMGPINDEX), CSMNORMSAVE, ENERGY, VNEW(3*NATOMS)

      LOGICAL SOCOUPLE
      DOUBLE PRECISION EPLUS, EMINUS, GRADDUM(3*NATOMS), DIFF

      ! }}}
      ! common {{{
      ! list of energies, EA(1) is the total one 
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /EA/ EA
      COMMON /CO/ COMPON
      COMMON /FAIL/ FTEST
      COMMON /EV/ EVAP, EVAPREJECT
      COMMON /PCALL/ NPCALL
      COMMON /CSMAVVAL/ AVVAL, CSMRMS, CSMIT
      COMMON /TOT/ NQTOT

      ! }}}
      ! }}}
      ! body {{{
!op226> potentials  {{{ 

      GUIDECHANGET=.FALSE.
      BRUN=0

      ! BRUN ==1 {{{
!
!  Test BRUN to see if we should stop if a screen saver is interrupted.
!  Need to save a restart file containing:
!  Current minimum in the Markov chain. COORDS
!  Number of steps done. NQTOT/NPAR should be close enough!
!  The current lowest minima. QMIN has the energies, QMINP has the points.
!  The current values of the temperature, acceptance ratio and step length,
!  TEMP(JP), ACCRAT(JP), STEP(JP), ASTEP(JP)
!  which can get changed dynamically.
!

      IF (BRUN.EQ.1) THEN
         WRITE(LFH,'(A)' ) 'dumping restart file ssdump'
         OPEN(UNIT=88,FILE='ssdump',STATUS='UNKNOWN')
         WRITE(88,'(3G20.10)') ((COORDS(J1,J2),J1=1,3*NATOMS),J2=1,NPAR)
         WRITE(88,'(I6)') NQTOT/NPAR, NPCALL
         WRITE(88,'(G20.10)') (QMIN(J1),J1=1,NSAVE)
         WRITE(88,'(3G20.10)') ((QMINP(J2,J1),J1=1,3*NATOMS),J2=1,NSAVE)
         WRITE(88,'(G20.10)') (TEMP(J1),J1=1,NPAR)
         WRITE(88,'(G20.10)') (ACCRAT(J1),J1=1,NPAR)
         WRITE(88,'(G20.10)') (STEP(J1),J1=1,NPAR)
         WRITE(88,'(G20.10)') (ASTEP(J1),J1=1,NPAR)
         WRITE(88,'(G20.10)') (OSTEP(J1),J1=1,NPAR)
         CALL SYSTEM('rm ssave')
         STOP
      ENDIF
      ! }}}

      NPCALL=NPCALL+1
      
10    CONTINUE

      ! IF...THEN {{{
      IF (MYBLNT) THEN 
         CALL EBLN(NATOMS,X,EA,GRAD,HESS,BLNTYPE,GRADT,SECT)
         EREAL=EA(1)
      ELSE IF (P46) THEN
         CALL P46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
      ELSE IF (G46) THEN
         CALL G46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
      ELSE IF (BLNT) THEN
         CALL BLN(X,GRAD,EREAL,GRADT)
      !ELSE
!
!  RAD must be called before the routine that calculates the potential or LBFGS
!  will get confused even if EVAP is set .TRUE. correctly.
!
!        CALL RAD(X,GRAD,EREAL,GRADT)
         !IF (EVAPREJECT) return
         !IF (CUTT) THEN
            !CALL LJCUT(X,GRAD,EREAL,GRADT,SECT)
         !ELSE
            !CALL LJ(X,GRAD,EREAL,GRADT,SECT)
         !ENDIF
      ENDIF
      ! }}}
! }}}
      ! Add Fields {{{

      IF (PULLT) THEN
         EREAL=EREAL-PFORCE*(X(3*(PATOM1-1)+3)-X(3*(PATOM2-1)+3))
         GRAD(3*(PATOM1-1)+3)=GRAD(3*(PATOM1-1)+3)-PFORCE
         GRAD(3*(PATOM2-1)+3)=GRAD(3*(PATOM2-1)+3)+PFORCE
      ENDIF
      ! }}}
      ! other {{{
      !IF (COMPON) CALL COMPRESS(X,GRAD,EREAL,GRADT)

      IF (GRADT.OR.CSMT) THEN
          IF (FREEZE) THEN
            DO J1=1,NATOMS
               IF (FROZEN(J1)) THEN
                  GRAD(3*(J1-1)+1)=0.0D0
                  GRAD(3*(J1-1)+2)=0.0D0
                  GRAD(3*(J1-1)+3)=0.0D0
               ENDIF
            ENDDO
         ENDIF
       ENDIF

         IF (SEEDT.AND.FREEZECORE) THEN
            DO J3=3*(NATOMS-NSEED)+1,3*NATOMS
               GRAD(J3)=0.0D0
            ENDDO
         ENDIF
         RMS=0.0D0
         DUMMY2=0.0D0

         ! }}}
         RMS=MAX(DSQRT(SUM(GRAD**2)/3*NATOMS),1.0D-100)

      RETURN
!op226>}}} 
      END SUBROUTINE POTENTIAL

      
