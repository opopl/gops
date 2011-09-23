
      SUBROUTINE IOM
      ! modules {{{
      USE COMMONS
      USE V
      USE PORFUNCS
      ! }}}
      IMPLICIT NONE

      ! loc {{{
      INTEGER J1,J2,JP
      LOGICAL EXISTS,YESNO
      INTEGER NPCALL,NQTOT
     ! DOUBLE PRECISION VECMN, DUMMY, EPSAB, EPSBB, SIGAB, SIGBB
      !LOGICAL END, YESNO, EXISTS
      !INTEGER J1, J2, JP, NTYPEA, NTYPE(105), ISTAT
      !INTEGER MP, LP
      !INTEGER Iostatus
      !DOUBLE PRECISION GTOL,STPMIN,STPMAX
      !CHARACTER FNAME*80, TSTRING*80
      !COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      !COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      !INTEGER  NQTOT, NPCALL
      !COMMON /TOT/ NQTOT
      !COMMON /PCALL/ NPCALL
      !DOUBLE PRECISION EPS2, RAD, HEIGHT
      !COMMON /CAPS/ EPS2, RAD, HEIGHT
      !LOGICAL SKIPBL, CLEAR, ECHO, CAT
      !INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST
      !COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     !&                NERROR, IR, ECHO, LAST, CAT
      ! }}}

      !IF (DL_POLY) THEN
         !OPEN(UNIT=91,FILE='CONFIG',STATUS='OLD')
         !READ(91,'(A1)') DUMMY
         !READ(91,'(A1)') DUMMY
!13       READ(91,'(A1)',END=14) DUMMY
         !NATOMS=NATOMS+1
         !READ(91,*) COORDS(3*(NATOMS-1)+1,1),COORDS(3*(NATOMS-1)+2,1),COORDS(3*(NATOMS-1)+3,1)
         !READ(91,'(A1)') DUMMY
         !READ(91,'(A1)') DUMMY
!C        WRITE(LFH,'(3G20.10)') COORDS(3*(NATOMS-1)+1,1),COORDS(3*(NATOMS-1)+2,1),COORDS(3*(NATOMS-1)+3,1)
         !GOTO 13
!14       CONTINUE
         !CLOSE(91)
      !ELSEIF (AMHT) THEN
          !write(LFH,'(A)')'DUMMY    '
      !IF (.NOT.(AMBERT.OR.AMBER.OR.CPMD.OR.CHRMMT)) THEN
!         CLOSE(7)!{{{
         !OPEN(UNIT=7,FILE='coords',STATUS='OLD')
         !IF (MOD(NATOMS,NPAR).NE.0) THEN
            !WRITE(LFH,'(A,I8,A,I8)') 'Number of atoms in coords file=',NATOMS,
     !&               ' is not divisible by number of runs=',NPAR
            !STOP
         !ENDIF
         !NATOMS=NATOMS/NPAR
         !IF (PERMOPT) THEN
            !OPEN(UNIT=17,FILE='finish',STATUS='OLD')
            !READ(17,*) (FIN(J1),J1=1,3*NATOMS)
            !WRITE(LFH,'(A)') 'Target coordinates read from file finish'
         !ENDIF
         !IF (TSALLIST.AND.(QTSALLIS.EQ.0)) QTSALLIS=1.0D0+1.0D0/(3*NATOMS)
         !IF (DFTBT.OR.TOSI.OR.WELCH) THEN
            !IR=7
!C           REWIND(7) ! this seems to be needed now?
            !DO JP=1,NPAR
               !DO J1=1,NATOMS
                  !CALL INPUT(END)
                  !J2=3*(J1-1)
                  !CALL READA(ZSYM(J1))
                  !CALL MYUPCASE(ZSYM(J1))
                  !IF (ZSYM(J1).EQ.'H ') IATNUM(J1+1)=1
                  !IF (ZSYM(J1).EQ.'HE') IATNUM(J1+1)=2
                  !IF (ZSYM(J1).EQ.'LI') IATNUM(J1+1)=3
                  !IF (ZSYM(J1).EQ.'BE') IATNUM(J1+1)=4
                  !IF (ZSYM(J1).EQ.'B ') IATNUM(J1+1)=5
                  !IF (ZSYM(J1).EQ.'C ') IATNUM(J1+1)=6
                  !IF (ZSYM(J1).EQ.'N ') IATNUM(J1+1)=7
                  !IF (ZSYM(J1).EQ.'O ') IATNUM(J1+1)=8
                  !IF (ZSYM(J1).EQ.'F ') IATNUM(J1+1)=9
                  !IF (ZSYM(J1).EQ.'S ') IATNUM(J1+1)=18
                  !CALL READF(COORDS(J2+1,JP))
                  !CALL READF(COORDS(J2+2,JP))
                  !CALL READF(COORDS(J2+3,JP))
               !ENDDO
            !ENDDO
         !ELSE
            !rewind(7)
            !DO JP=1,NPAR
               !DO J1=1,NATOMS
                  !J2=3*(J1-1)
                   !READ(7,*) COORDS(J2+1,JP), COORDS(J2+2,JP), COORDS(J2+3,JP)
               !ENDDO
            !ENDDO
         !ENDIF
         !CLOSE(7)!}}}
!      ELSE IF (AMBERT) THEN
!         DO JP=1,NPAR
!             COORDS(:,JP)=COORDS(:,1)   ! we have already read the coordinates for AMBER runs into this array
!         END DO
      !ENDIF

      !IF (CPMD) THEN
         !FNAME=SYS(1:LSYS)
         !WRITE(LFH,'(A,A)') ' Reading coordinates from file ',FNAME
         !OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         !NATOMS=0
!111      READ(7,'(A)') FNAME
         !IF (FNAME(1:6).EQ.'&ATOMS') THEN
            !J1=0
!121         READ(7,'(A)') TSTRING
            !IF (TSTRING(1:1).EQ.'*') THEN
               !J1=J1+1
               !READ(7,'(A)') FNAME
               !READ(7,*) NTYPE(J1)
               !DO J2=1,NTYPE(J1)
                  !IATNUM(NATOMS+J2)=1
                  !ZSYM(NATOMS+J2)='CP'
                  !READ(7,*) COORDS(3*(NATOMS+J2-1)+1,1),COORDS(3*(NATOMS+J2-1)+2,1),COORDS(3*(NATOMS+J2-1)+3,1)
               !ENDDO
               !NATOMS=NATOMS+NTYPE(J1)
               !GOTO 121
            !ELSE IF (TSTRING(1:1).EQ.' ') THEN
               !GOTO 121
            !ELSE IF (TSTRING(1:4).EQ.'&END') THEN
               !GOTO 131
            !ENDIF
         !ELSE
            !GOTO 111
         !ENDIF

!131      CONTINUE
         !CLOSE(7)

         !CALL SYSTEM(' grep -c ANGSTROM ' // SYS(1:LSYS) // ' > temp')
         !OPEN(UNIT=7,FILE='temp',STATUS='OLD')
         !READ(7,*) J1
         !CLOSE(7)
         !IF (J1.EQ.1) THEN
            !WRITE(LFH,'(A)') ' Converting initial coordinates from Angstrom to Bohr'
            !DO J1=1,3*NATOMS
               !COORDS(J1,1)=COORDS(J1,1)*1.889726164D0
            !ENDDO
         !ENDIF
      !ENDIF

      !IF (RESIZET) THEN
         !WRITE(LFH,'(A,F15.7)') 'Multiplying coordinates by ',RESIZE
         !DO JP=1,NPAR
            !DO J1=1,3*NATOMS
               !COORDS(J1,JP)=COORDS(J1,JP)*RESIZE
            !ENDDO
         !ENDDO
      !ENDIF

         CLOSE(7)!{{{
         OPEN(UNIT=7,FILE='coords',STATUS='OLD')
         IF (MOD(NATOMS,NPAR).NE.0) THEN
            WRITE(LFH,'(A,I8,A,I8)') 'Number of atoms in coords file=',NATOMS,&
     &               ' is not divisible by number of runs=',NPAR
            STOP
         ENDIF
         NATOMS=NATOMS/NPAR
         IF (PERMOPT) THEN
            OPEN(UNIT=17,FILE='finish',STATUS='OLD')
            READ(17,*) (FIN(J1),J1=1,3*NATOMS)
            WRITE(LFH,'(A)') 'Target coordinates read from file finish'
         ENDIF
         REWIND(7)
            DO JP=1,NPAR
               DO J1=1,NATOMS
                  J2=3*(J1-1)
                   READ(7,*) COORDS(J2+1,JP), COORDS(J2+2,JP), COORDS(J2+3,JP)
               ENDDO
            ENDDO
         CLOSE(7)!}}}

      IF (.NOT.SEEDT.AND..NOT.AMHT) THEN
         WRITE(LFH,20) 
20       FORMAT('Initial coordinates:')
         DO JP=1,NPAR
               WRITE(LFH,30) (COORDS(J1,JP),J1=1,3*NATOMS)
30             FORMAT(3F20.10)
         ENDDO
      ENDIF

      IF (P46) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' 3-colour, 46 bead model polypeptide'
      ELSEIF (BLNT) THEN
         WRITE(LFH,'(I4,A)') NATOMS,' bead BLN model'
      ENDIF
     
      IF (DEBUG.OR.CHECKMARKOVT) WRITE(LFH,'(A,I6,A)') 'io1> checking the energy of the saved coordinates in the chain'
      IF (FIXCOM .AND. (.NOT. CHRMMT)) THEN
          INQUIRE(FILE='masses',EXIST=EXISTS)
          IF (EXISTS) THEN
              OPEN(UNIT=1978,FILE="masses")
              DO J2=1,NATOMS
                 READ(1978,'(F12.5)') MASSES(J2)
              ENDDO
              CLOSE(1978)
          ELSE
             WRITE(LFH,'(A)') 'WARNING, FIXCOM is specified, but "masses" file is not present: setting all masses to unity.'
             MASSES(1:NATOMS) = 1.0D0
          ENDIF
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
      
      DO J1=1,NPAR
         IF (JUMPMOVE(J1)) WRITE(LFH,'(A,I1,A,I1,A,I4,A)') &
     & 'Jump moves will be attempted from run ',J1,' to run ',JUMPTO(J1),' every ',JUMPINT(J1),' steps'
      ENDDO
      IF (NEWJUMP) WRITE(LFH,'(A,F12.3)')  &
     & 'Jumping based only on current energies (parallel tempering) attempt probability=',PNEWJUMP
      WRITE(LFH,'(A,F15.10)') 'Sloppy quench tolerance for RMS gradient ',BQMAX
      !IF ((.NOT.BSPT).AND.(.NOT.PTMC)) THEN
!         DO JP=1,NPAR
            !IF (FIXBOTH(JP)) THEN
               !WRITE(LFH,'(A,I3,A,F12.4,A,2F12.4,A)') 
     !& 'In run ',JP,' temperature=',TEMP(JP),' step size and angular threshold=',
     !& STEP(JP),ASTEP(JP),' all fixed'
               !IF (RIGID) WRITE(LFH,'(A,F12.4)') 'Orientational step size for rigid bodies also fixed at ',
     !& OSTEP(JP)
               !IF (FIXD) WRITE(LFH,'(A,I3,A,I4)') 'In run ',JP,' number of hard sphere collision moves fixed at ',
     !& NHSMOVE
            !ELSE IF (FIXSTEP(JP)) THEN
               !WRITE(LFH,'(A,I3,A,2F12.4)') 'In run ',JP,' step size and angular threshold fixed at ',
     !& STEP(JP),ASTEP(JP)
               !IF (RIGID) WRITE(LFH,'(A,F12.4)') 'Orientational step size for rigid bodies also fixed at ',
     !& OSTEP(JP)
               !IF (FIXD) WRITE(LFH,'(A,I3,A,I4)') 'In run ',JP,' number of hard sphere collision moves fixed at ',
     !& NHSMOVE
               !IF (.NOT.FIXTEMP(JP)) THEN
                  !WRITE(LFH,'(A,F12.4,A,F12.4)') 
     !& 'Temperature will be adjusted for acceptance ratio ',ACCRAT(JP),' initial value=',TEMP(JP)
               !ELSE
                  !WRITE(LFH,'(A,I1,A,G12.4)') 'In run ',JP,' temperature fixed at ',TEMP(JP)
               !ENDIF
            !ELSE IF (STEPOUT) THEN
               !WRITE(LFH,'(A,I3,A,2F12.4,A,2F12.4)') 
     !& 'In run ',JP,' step size and angular threshold will be adjusted to escape from basins. Initial values=',
     !& STEP(JP),ASTEP(JP)
               !IF (RIGID) WRITE(LFH,'(A,F12.4)') 'Orientational step size for rigid bodies initial value ',
     !& OSTEP(JP)
               !IF (FIXD) WRITE(LFH,'(A,I3,A,I4)') 
     !& 'In run ',JP,' number of hard sphere collision moves will be adjusted. Initial value=',NHSMOVE
               !IF (.NOT.FIXTEMP(JP)) THEN
                  !WRITE(LFH,'(A,F12.4,A,F12.4)') 
     !& 'Temperature will be adjusted for acceptance ratio ',ACCRAT(JP),' initial value=',TEMP(JP)
               !ELSE
                  !WRITE(LFH,'(A,I1,A,G12.4)') 'In run ',JP,' temperature fixed at ',TEMP(JP)
               !ENDIF
            !ELSE 
               !WRITE(LFH,'(A,I3,A,G12.4)') 'In run ',JP,' temperature fixed at ',TEMP(JP)
               !WRITE(LFH,'(A,F12.4,A,2F12.4)') 'Step size and angular threshold will be adjusted for acceptance ratio ',
     !& ACCRAT(JP),' initial values=',STEP(JP),ASTEP(JP)
               !IF (RIGID) WRITE(LFH,'(A,F12.4)') 'Orientational step size for rigid bodies initial value ',
     !& OSTEP(JP)
               !IF (FIXD) WRITE(LFH,'(A,I3,A,I4)') 
     !& 'In run ',JP,' number of hard sphere collision moves will be adjusted. Initial value=',NHSMOVE
            !ENDIF
         !ENDDO 
      !ENDIF
      IF (EFAC.NE.0.0D0) WRITE(LFH,'(A,F12.4)') 'Exponential factor for proposed steps=',EFAC
      IF (NORESET.OR.BSPT) THEN
         WRITE(LFH,'(A)') 'Configuration will not be reset to quench geometry'
         IF (CENT) THEN
            WRITE(LFH,'(A)') 'WARNING CENTRE can lead to atoms leaving '
            WRITE(LFH,'(A)') 'the container after takestep when the centre of mass is moved.'
!           STOP
         ENDIF
      ELSE
         WRITE(LFH,'(A)') 'Configuration will be reset to quench geometry'
      ENDIF
     ! IF (CENT .AND. FIXCOM) THEN
          !WRITE(LFH,'(A)') 'WARNING: keywords CENTRE (fixing centre of coordinates) and FIXCOM (fixing centre of mass) 
     !& are incompatible'
          !STOP
      !ENDIF

      WRITE(LFH,'(A)') 'Sampling using Boltzmann weights'
      !IF (PERIODIC) WRITE(LFH,'(A,3F15.7)') 'Periodic boundary conditions, box lengths: ',BOXLX,BOXLY,BOXLZ
      IF (CUTT) WRITE(LFH,'(A,F15.7)') 'Cutoff=',CUTOFF
       WRITE (LFH,'(A)') 'Nocedal LBFGS minimization'
       WRITE(LFH,'(A,I6)') 'Number of updates before reset in LBFGS=',MUPDATE
       WRITE(LFH,'(A,F20.10)') 'Maximum step size=',MAXBFGS
       WRITE(LFH,'(A,G12.4)') 'Guess for initial diagonal elements in LBFGS=',DGUESS
      WRITE(LFH,'(A,F15.10)') 'Final quench tolerance for RMS gradient ',CQMAX
      WRITE(LFH,'(A,F15.10)') 'Energy difference criterion for minima=',ECONV
      WRITE(LFH,'(A,I5,A,I5)') 'Maximum number of iterations: sloppy quenches ',MAXIT,' final quenches ',MAXIT2
      IF (.NOT.BSPT) THEN
         DO J1=1,NRUNS
            WRITE(LFH,120) J1, MCSTEPS(J1), TFAC(J1)
120         FORMAT('Run ',I3,': ',I9,' steps with temperature scaled by ',E15.8)
         ENDDO
      ENDIF
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
      !IF (RENORM) THEN
         !XMOVERENORM=MIN(XMOVERENORM,NATOMS*0.9D0)
         !IF (XMOVERENORM.GT.3.0D0) THEN
            !WRITE(LFH,'(A,F12.1,A)')  'Large steps of ',XMOVERENORM,' hard sphere type moves will be used:'
         !ELSE
            !WRITE(LFH,'(A,F15.5,A)')  'Large steps using maximum displacement ',XMOVERENORM,' will be used:'
         !ENDIF
         !WRITE(LFH,'(A,I6)')    '  Initial interval for large steps is ',NRENORM
         !WRITE(LFH,'(A,F12.5)') '  Temperature used in Metropolis test for acceptance of large steps is ',TRENORM
      !ENDIF
      !IF (RESTORET) THEN
         !WRITE(LFH,'(A,A)') 'Restoring GMIN run from file ',TRIM(ADJUSTL(DUMPFILE))
      !ENDIF
      !IF (NEWRESTART) THEN
         !IF (.NOT.AVOIDRESEEDT) THEN
            !WRITE(LFH,'(A,F12.5,A,I6,A)') 'Steps will be rejected and taboo list populated if the energy decrease < ',ECONV,
     !&                                ' within ',NRELAX,' steps'
         !ELSE
            !WRITE(LFH,'(A,F12.5,A,I6,A)') 'Runs will be reseeded if the energy does not decrease by at least ',ECONV,
     !&                                ' within ',NRELAX,' steps'
         !ENDIF
         !IF (NHSRESTART.GT.0)  WRITE(LFH, '(I6,A)') NHSRESTART,' hard sphere-type moves will be used to reseed'
      !ENDIF
      !IF (AVOID) THEN
         !IF (AVOIDRESEEDT) THEN
            !WRITE(LFH,'(A,F10.2,A,I6,A)') 'Runs will be reseeded if the current minimum comes within ',
     !& AVOIDDIST,' of up to ',MAXSAVE,' previous minima'
         !ELSE
            !WRITE(LFH,'(A,F10.2,A,I6,A)') 'Steps will be rejected if the current minimum comes within ',
     !& AVOIDDIST,' of up to ',MAXSAVE,' previous minima'
         !ENDIF
      !ENDIF
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
      END
