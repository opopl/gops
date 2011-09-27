
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
      INTEGER J1,J2,J3
      INTEGER NQTOT,NPCALL

      LOGICAL GUIDECHANGET, GUIDET, CSMDOGUIDET
      ! }}}
      ! common {{{
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /PCALL/ NPCALL
      COMMON /TOT/ NQTOT

      ! }}}
      ! }}}
      ! body {{{
!op226> potentials  {{{ 
      ! com {{{
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
         WRITE(88,'(I6)') NQTOT, NPCALL
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
      ! }}}
      !IF...THEN {{{
      EA(1:10)=0.0D0
      IF (MYBLNT) THEN 
         CALL EBLN(NATOMS,X,EA,GRAD,HESS,BLNTYPE,GRADT,SECT,EA_FH,.false.)
      ELSE IF (P46.OR.G46) THEN
         CALL EP46(EA_FH,DEB_BLN,X,NATOMS,GRAD,EA,GRADT,BLNTYPE); EREAL=EA(1)
         !CALL P46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
      !ELSE
      ENDIF
      ! }}}
! }}}
      ! Add Fields {{{

      IF (PULLT) THEN
         dE_fz=-PFORCE*(X(3*(PATOM1-1)+3)-X(3*(PATOM2-1)+3))
         EREAL=EREAL+dE_fz
         GRAD(3*(PATOM1-1)+3)=GRAD(3*(PATOM1-1)+3)-PFORCE
         GRAD(3*(PATOM2-1)+3)=GRAD(3*(PATOM2-1)+3)+PFORCE
      ENDIF
      EA(6)=dE_fz
      ! }}}
      ! other {{{
      
         IF (SEEDT.AND.FREEZECORE) THEN
            DO J3=3*(NATOMS-NSEED)+1,3*NATOMS
               GRAD(J3)=0.0D0
            ENDDO
         ENDIF

         ! }}}
         RMS=MAX(DSQRT(SUM(GRAD**2)/3*NATOMS),1.0D-100)

      RETURN
!op226>}}} 
      END SUBROUTINE POTENTIAL

      
