
      SUBROUTINE IOM
      ! modules {{{
      USE COMMONS
      USE V
      USE F
      USE PORFUNCS
      ! }}}
      IMPLICIT NONE
      ! loc {{{
      INTEGER J1,J2,JP
      LOGICAL EXISTS,YESNO
      INTEGER NPCALL,NQTOT
      DOUBLE PRECISION,DIMENSION(NATOMS,3) :: R
      ! }}}
      ! body {{{
      JP=1

      ! read in coordinates from C_FILE {{{
         OPEN(UNIT=7,FILE=C_FILE,STATUS='OLD')
         REWIND(7)
               DO J1=1,NATOMS
                   READ(7,*) R(J1,1:3)
               ENDDO
         CLOSE(7)
      !}}}

      !IF (.NOT.SEEDT.AND..NOT.AMHT) THEN
         WRITE(LFH,20) 
20       FORMAT('Initial coordinates:')
         WRITE(LFH,30) R
30       FORMAT(3F20.10)
      !ENDIF
      COORDS(1:NR,1)=PACK( R, .true. )

      CALL ED(LFH) 
      IF (P46) THEN
        WRITE(LFH,'(I4,A)') NATOMS,' atoms, Wild-type BLN model'
      ELSEIF (G46) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' atoms, Go-like BLN model'
      ELSEIF (BLNT) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' atoms, general BLN model'
      ENDIF
      CALL ED(LFH) 
      
      IF (RADIUS.EQ.0.0D0) THEN
         RADIUS=2.0D0+(3.0D0*NATOMS/17.77153175D0)**(1.0D0/3.0D0)
         IF (P46) THEN
            RADIUS=RADIUS*3.0D0
         ELSE 
            RADIUS=RADIUS*2.0D0**(1.0D0/6.0D0)
         ENDIF
      ENDIF

      RADIUS=RADIUS**2
      
      IF (NORESET) THEN
         WRITE(LFH,'(A)') 'Configuration will not be reset to quench geometry'
         IF (CENT) THEN
            WRITE(LFH,'(A)') 'WARNING CENTRE can lead to atoms leaving '
            WRITE(LFH,'(A)') 'the container after takestep when the centre of mass is moved.'
!           STOP
         ENDIF
      ELSE
         WRITE(LFH,'(A)') 'Configuration will be reset to quench geometry'
      ENDIF

      WRITE(LFH,'(A)') 'Sampling using Boltzmann weights'
      IF (CUTT) WRITE(LFH,'(A,F15.7)') 'Cutoff=',CUTOFF

      CALL ED(LFH)
        WRITE(LFH,'(A,/)') 'Minimization:'
        WRITE (LFH,'(A)') 'Nocedal LBFGS minimization'
        WRITE(LFH,'(A,I6)') 'Number of updates before reset in LBFGS=',MUPDATE
        WRITE(LFH,'(A,F20.10)') 'Maximum step size=',MAXBFGS
        WRITE(LFH,'(A,G12.4)') 'Guess for initial diagonal elements in LBFGS=',DGUESS
        WRITE(LFH, '(A,G20.10)') 'Maximum allowed energy rise during a minimisation=',MAXERISE
      CALL ED(LFH)
        WRITE(LFH,'(A,/)') 'Quenches:'
        WRITE(LFH,'(A,F15.10)') 'Sloppy quench tolerance for RMS gradient ',BQMAX
        IF (EFAC.NE.0.0D0) WRITE(LFH,'(A,F12.4)') 'Exponential factor for proposed steps=',EFAC
        WRITE(LFH,'(A,F15.10)') 'Final quench tolerance for RMS gradient ',CQMAX
        WRITE(LFH,'(A,F15.10)') 'Energy difference criterion for minima=',ECONV
        WRITE(LFH,'(A,I5,A,I5)') 'Maximum number of iterations: sloppy quenches ',MAXIT,' final quenches ',MAXIT2
      CALL ED(LFH)
        IF (PULLT) THEN
            WRITE(LFH,'(A,/)') 'Pulling:'
            WRITE(LFH,'(A,F15.10)') 'Force: ',PFORCE
            WRITE(LFH,'(A,F15.10)') 'Atom 1: ',PATOM1
            WRITE(LFH,'(A,F15.10)') 'Atom 2: ',PATOM2
            CALL ED(LFH)
        ENDIF

      WRITE(LFH,120) MCSTEPS(1), TFAC(1)
120   FORMAT(I9,' steps with temperature scaled by ',E15.8)

      IF (DEBUG) THEN
         WRITE(LFH,160) 
160      FORMAT('Debug printing is on')
      ENDIF
      IF (TARGET) THEN
         WRITE(LFH,'(A)',ADVANCE='NO') 'Target energies: '
         WRITE(LFH,'(F20.10)',ADVANCE='NO') (TARGETS(J1),J1=1,NTARGETS)
         WRITE(LFH,'(A)') ' '
      ENDIF
     ! ssdump {{{
!  Look for the file that contains interrupted screen saver restart information.
!  Current minimum in the Markov chain. COORDS
!  Number of steps done. NQTOT/NPAR should be close enough!
!  The current lowest minima. QMIN has the energies, QMINP has the points.
!  The current values of the temperature, acceptance ratio and step length,
!  TEMP(JP), ACCRAT(JP), STEP(JP), ASTEP(JP) and OSTEP(JP)
!  which can get changed dynamically.

      INQUIRE(FILE='ssdump',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(UNIT=88,FILE='ssdump',STATUS='UNKNOWN')
          WRITE(LFH,'(A)') 'reading dump information from file ssdump'
         READ(88,'(3G20.10)') ((COORDS(J1,J2),J1=1,3*NATOMS),J2=1,NPAR)
         READ(88,'(2I6)') NQTOT, NPCALL
         MCSTEPS(1)=MAX(MCSTEPS(1)-NQTOT*NPAR,1)
         NQTOT=NQTOT*NPAR
         READ(88,'(G20.10)') (QMIN(J1),J1=1,NSAVE)
         READ(88,'(3G20.10)') ((QMINP(J2,J1),J1=1,3*NATOMS),J2=1,NSAVE)
         READ(88,'(G20.10)') (TEMP(J1),J1=1,NPAR)
         READ(88,'(G20.10)') (ACCRAT(J1),J1=1,NPAR)
         READ(88,'(G20.10)') (STEP(J1),J1=1,NPAR)
         READ(88,'(G20.10)') (ASTEP(J1),J1=1,NPAR)
         READ(88,'(G20.10)') (OSTEP(J1),J1=1,NPAR)
         CLOSE(88)
         YESNO=.FALSE.
      ENDIF
! }}}

      RETURN
      ! }}}
      END
