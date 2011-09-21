
      SUBROUTINE POTENTIAL(X,GRAD,EREAL,GRADT,SECT)
!op226> Declarations {{{ 

      USE COMMONS
      USE QMODULE
      USE PORFUNCS

      IMPLICIT NONE
      
      LOGICAL GRADT, SECT, EVAP, COMPON, YESNO, GUIDET, EVAPREJECT
      INTEGER J1, J2, J3, NPCALL
      DOUBLE PRECISION EREAL, GRAD(*), X(*), DUMMY2
      DOUBLE PRECISION ENERGY, VNEW(3*NATOMS)

      LOGICAL GUIDECHANGET
      COMMON /CO/ COMPON
      COMMON /EV/ EVAP, EVAPREJECT
      COMMON /PCALL/ NPCALL
      INTEGER BRUN

      INTEGER NQTOT
      COMMON /TOT/ NQTOT

      DOUBLE PRECISION EPLUS, EMINUS, GRADDUM(3*NATOMS), DIFF

Cop226> }}}
!op226>{{{ 
      BRUN=0

C  Test BRUN to see if we should stop if a screen saver is interrupted.
C  Need to save a restart file containing:
C  Current minimum in the Markov chain. COORDS
C  Number of steps done. NQTOT/NPAR should be close enough!
C  The current lowest minima. QMIN has the energies, QMINP has the points.
C  The current values of the temperature, acceptance ratio and step length,
C  TEMP(JP), ACCRAT(JP), STEP(JP), ASTEP(JP)
C  which can get changed dynamically.
C
      ! ssdump {{{

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

      IF (P46) THEN
         CALL P46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
      ELSE IF (G46) THEN
         CALL G46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
      ELSE IF (BLNT) THEN
         CALL BLN(X,GRAD,EREAL,GRADT)
      ENDIF

Cop226>}}} 

C  --------------- End of possible potentials - now add fields if required ------------------------------

      IF (PULLT) THEN
         EREAL=EREAL-PFORCE*(X(3*(PATOM1-1)+3)-X(3*(PATOM2-1)+3))
         GRAD(3*(PATOM1-1)+3)=GRAD(3*(PATOM1-1)+3)-PFORCE
         GRAD(3*(PATOM2-1)+3)=GRAD(3*(PATOM2-1)+3)+PFORCE
      ENDIF

      RMS=0.0D0
      DUMMY2=0.0D0
      DO J1=1,3*NATOMS
            DUMMY2=DUMMY2+GRAD(J1)**2 
      ENDDO
      RMS=MAX(DSQRT(DUMMY2/(3*NATOMS)),1.0D-100)

      RETURN
      END
