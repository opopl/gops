
      SUBROUTINE IO1

      USE COMMONS
      USE QMODULE
      USE PORFUNCS

      IMPLICIT NONE

      LOGICAL END, YESNO 
      INTEGER :: NQTOT,NPCALL
      INTEGER ::    J1,J2,jp

      ! look ssdump
      COMMON /TOT/ NQTOT
      COMMON /PCALL/ NPCALL

         CLOSE(7)
         OPEN(UNIT=7,FILE=C_FILE,STATUS='OLD')
         REWIND(7)
               DO J1=1,NATOMS
                  J2=3*(J1-1)
                   READ(7,*) COORDS(J2+1), COORDS(J2+2), COORDS(J2+3)
               ENDDO
         CLOSE(7)
      
      IF (.NOT.SEEDT) THEN
         WRITE(LFH,20) 
20       FORMAT('Initial coordinates:')
               WRITE(LFH,30) (COORDS(J1),J1=1,3*NATOMS)
30             FORMAT(3F20.10)
      ENDIF

      IF (P46) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' 3-colour, 46 bead model polypeptide'
      ELSE IF (BLNT) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' bead BLN model'
      ENDIF

      IF (DEBUG.OR.CHECKMARKOVT) THEN 
         WRITE(LFH,'(A,I6,A)') 'io1> checking the energy of the saved coordinates in the chain'
      ENDIF

      IF (RADIUS.EQ.0.0D0) THEN
         RADIUS=2.0D0+(3.0D0*NATOMS/17.77153175D0)**(1.0D0/3.0D0)
         IF (P46) THEN
            RADIUS=RADIUS*3.0D0
         ELSE 
            RADIUS=RADIUS*2.0D0**(1.0D0/6.0D0)
         ENDIF
      ENDIF

      RADIUS=RADIUS**2

      WRITE(LFH,'(A,F15.10)') 'Sloppy quench tolerance for RMS gradient ',SQMAX

      ! check FIXBOTH STEPOUT FIXSTEP FIXTEMP {{{
            IF (FIXBOTH) THEN
               WRITE(LFH,'(A,I3,A,F12.4,A,2F12.4,A)') 'In run ',0,&
                & ' temperature=',TEMP,&
                & ' step size and angular threshold=', STEP,ASTEP,&
                ' all fixed'
            ELSE IF (FIXSTEP) THEN
               WRITE(LFH,'(A,I3,A,2F12.4)') 'In run ',0,&
                & ' step size and angular threshold fixed at ',&
                & STEP,ASTEP
               IF (.NOT.FIXTEMP) THEN
                  WRITE(LFH,'(A,F12.4,A,F12.4)') 'Temperature will be &
                    adjusted for acceptance ratio ',& 
                    ACCRAT,' initial value=',TEMP
               ELSE
                  WRITE(LFH,'(A,I1,A,G12.4)') 'In run ',0,' temperature fixed at ',TEMP
               ENDIF
            ELSE IF (STEPOUT) THEN
               WRITE(LFH,'(A,I3,A,2F12.4,A,2F12.4)') 'In run ',& 
                    & 0,' step size and angular &
                    threshold will be adjusted to &
                    escape from basins. Initial values=',&
                   STEP,ASTEP
               IF (.NOT.FIXTEMP) THEN
                  WRITE(LFH,'(A,F12.4,A,F12.4)') & 
     &                    'Temperature will be adjusted for acceptance ratio ',ACCRAT,' initial value=',TEMP
               ELSE
                  WRITE(LFH,'(A,I1,A,G12.4)') 'In run ',0,' temperature fixed at ',TEMP
               ENDIF
            ELSE 
               WRITE(LFH,'(A,I3,A,G12.4)') 'In run ',0,' temperature fixed at ',TEMP
               WRITE(LFH,'(A,F12.4,A,2F12.4)') 'Step size and angular threshold will be adjusted for acceptance ratio ',&
     &                ACCRAT,' initial values=',STEP,ASTEP
            ENDIF
        ! }}}

      ! NORESET {{{
      IF (NORESET.OR.BSPT) THEN
         WRITE(LFH,'(A)') 'Configuration will not be reset to quench geometry'
      ELSE
         WRITE(LFH,'(A)') 'Configuration will be reset to quench geometry'
      ENDIF
      ! }}}

      IF (CENT .AND. FIXCOM) THEN  
           WRITE(LFH,'(A)') "WARNING: " 
           write(LFH,'(a)') "keywords CENTRE (fixing centre of coordinates) "
           write(LFH,'(a)') "and FIXCOM (fixing centre of mass) are incompatible"
          STOP
      ENDIF

      WRITE(LFH,'(A)') 'Sampling using Boltzmann weights'

      ! lbfgs {{{
      IF (BFGS) THEN
         WRITE (LFH,'(A)') 'BFGS minimization'
      ELSE IF (LBFGST) THEN
         WRITE (LFH,'(A)') 'Nocedal LBFGS minimization'
         WRITE(LFH,'(A,I6)') 'Number of updates before reset in LBFGS=',MUPDATE
         WRITE(LFH,'(A,F20.10)') 'Maximum step size=',MAXBFGS
         WRITE(LFH,'(A,G12.4)') 'Guess for initial diagonal elements in LBFGS=',DGUESS
      ELSE
         WRITE (LFH,'(A)') 'Conjugate gradient minimization'
      ENDIF
      ! }}}

      WRITE(LFH,'(A,F15.10)') 'Final quench tolerance for RMS gradient ',FQMAX
      WRITE(LFH,'(A,F15.10)') 'Energy difference criterion for minima=',EDIFF
      WRITE(LFH,'(A,I5,A,I5)') 'Maximum number of iterations: sloppy quenches ',MAXIT,' final quenches ',MAXIT2

       WRITE(LFH,120) MCSTEPS, TFAC
120    FORMAT(I6,' steps with temperature scaled by ',E15.8)

      IF (DEBUG) THEN
         WRITE(LFH,160) 
160      FORMAT('Debug printing is on')
      ENDIF

      WRITE(LFH, '(A,G20.10)') 'Maximum allowed energy rise during a minimisation=',MAXERISE

      IF (TARGET) THEN
         WRITE(LFH,'(A)',ADVANCE='NO') 'Target energies: '
         WRITE(LFH,'(F20.10)',ADVANCE='NO') (TARGETS(J1),J1=1,NTARGETS)
         WRITE(LFH,'(A)') ' '
      ENDIF
                  
      ! ssdump  {{{
!  Look for the file that contains interrupted screen saver restart information.
!  Current minimum in the Markov chain. COORDS
!  Number of steps done. NQTOT/NPAR should be close enough!
!  The current lowest minima. QMIN has the energies, QMINP has the points.
!  The current values of the temperature, acceptance ratio and step length,
!  TEMP, ACCRAT, STEP, ASTEP and OSTEP
!  which can get changed dynamically.

      INQUIRE(FILE='ssdump',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(UNIT=88,FILE='ssdump',STATUS='UNKNOWN')
         WRITE(LFH,'(A)') 'reading dump information from file ssdump'
         READ(88,'(3G20.10)') (COORDS(J1),J1=1,3*NATOMS)
         READ(88,'(2I6)') NQTOT, NPCALL
         READ(88,'(G20.10)') (QMIN(J1),J1=1,NSAVE)
         READ(88,'(3G20.10)') ((QMINP(J2,J1),J1=1,3*NATOMS),J2=1,NSAVE)
         READ(88,'(G20.10)') TEMP
         READ(88,'(G20.10)') ACCRAT 
         READ(88,'(G20.10)') STEP
         READ(88,'(G20.10)') ASTEP
         READ(88,'(G20.10)') OSTEP
         CLOSE(88)
         YESNO=.FALSE.
      ENDIF
      ! }}}


      RETURN
      END
