      SUBROUTINE POTENTIAL(X,GRAD,EREAL,GRADT,SECT)
      USE COMMONS
      USE QMODULE
      USE PERMU
      use porfuncs
      IMPLICIT NONE
      
      LOGICAL GRADT, FTEST, SECT, EVAP, COMPON, YESNO, GUIDET, evapreject
      INTEGER J1, J2, J3, NPCALL, PERM(NATOMS), NPERM, NORBIT1, NORBIT2, CSMIT
      DOUBLE PRECISION EREAL, GRAD(*), X(*), DUMMY2, GEMAX, XG, YG, ZG, RMAT(3,3), XD, YD, ZD, XTEMP(3*NATOMS),
     1                 GRADLJ(3*NATOMS), EREALLJ, GRADMF(3*NATOMS), EREALMF, TERMLJ, TERMMF, GRADDUM1(3*NATOMS), AVVAL,
     2                 GRADDUM2(3*NATOMS), EDUM1, EDUM2, DUMMY(3*NATOMS), DIST2, WORSTRAD, GBDUMMY, QE, QX, DISTANCE, AA(3),
     3                 SAVECSMNORM, CSMRMS, CSMGRAD(3), SAVECSMPMAT(3,3), SAVECSMIMAGES(3*NATOMS*CSMGPINDEX)
      INTEGER CSMGPINDEXSAVE
      DOUBLE PRECISION PTGPSAVE(3,3,2*CSMGPINDEX), CSMNORMSAVE, ENERGY, VNEW(3*NATOMS)

      LOGICAL SOCOUPLE, GUIDECHANGET, CSMDOGUIDET
      INTEGER BRUN

      INTEGER NQTOT

      DOUBLE PRECISION EPLUS, EMINUS, GRADDUM(3*NATOMS), DIFF



      GUIDECHANGET=.FALSE.
      BRUN=0

      IF (BRUN.EQ.1) THEN
         WRITE(MYUNIT,'(A)' ) 'dumping restart file ssdump'
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
         STOP
      ENDIF

      NPCALL=NPCALL+1
      
10    CONTINUE

      IF (MSORIGT) THEN
         IF (CUTT) THEN
         ELSE
         ENDIF
         IF (FTEST) THEN
            RETURN
         ENDIF
      ELSE IF (MSTRANST) THEN
         IF (FTEST) THEN
            RETURN
         ENDIF
      ELSE IF (FRAUSIT) THEN
         IF (FTEST) THEN
            RETURN
         ENDIF
      ELSE IF ((NEON.AND.PLUS).OR.(NEON.AND.STAR).OR.(ARGON.AND.PLUS).OR.(ARGON.AND.STAR)) THEN
      ELSE IF (TWOPLUS) THEN
      ELSE IF (GROUND) THEN
         IF (AXTELL) CALL AXT(NATOMS,X,GRAD,EREAL,GRADT,ZSTAR)
      ELSE IF (RGCL2) THEN
         IF (AXTELL) CALL AXT(NATOMS,X,GRAD,EREAL,GRADT,ZSTAR)
      ELSE IF (ARNO) THEN
         SOCOUPLE=.TRUE.

         IF (AXTELL) CALL AXT(NATOMS,X,GRAD,EREAL,GRADT,ZSTAR)
      ELSE IF (OTPT) THEN
      ELSE IF (LJMFT) THEN
      ELSE IF (MORSET) THEN
      ELSE IF (TOSI) THEN
         if (evapreject) return
      ELSE IF (WELCH) THEN
      ELSE IF (SCT) THEN
         IF (CPMD) EREAL=EREAL+1.0D6
      ELSE IF (ACKLANDT) THEN
      ELSE IF (FAL.OR.FNI) THEN
      ELSE IF (DFTBT) THEN
         IF (FTEST) THEN
            RETURN
         ENDIF
      ELSE IF (P46) THEN
      ELSE IF (G46) THEN
      ELSE IF (BLNT) THEN
      ELSE IF (CHAPERONINT) THEN
      ELSE IF (SW) THEN
         IF (PERIODIC) THEN
         ELSE
         ENDIF
      ELSE IF (AMHT) THEN


      ELSE IF (AMBERT) THEN
      ELSE IF (CHRMMT) THEN
      ELSE IF (DZTEST) THEN
      ELSE IF (ZETT1) THEN
      ELSE IF (ZETT2) THEN
      ELSE IF (PACHECO) THEN
         IF (AXTELL) CALL AXT(NATOMS,X,GRAD,EREAL,GRADT,ZSTAR)
      ELSE IF (MODEL1T) THEN
      ELSE IF (EAMLJT) THEN
      ELSE IF (PBGLUET) THEN
      ELSE IF (FST) THEN
      ELSE IF (GUPTAT) THEN
      ELSE IF (NATBT) THEN
      ELSE IF (ALGLUET) THEN
      ELSE IF (MGGLUET) THEN
      ELSE IF (EAMALT) THEN
      ELSE IF (WENZEL) THEN
      ELSE IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
         IF (CSMDOGUIDET) THEN ! we want to use a guiding group
            PTGPSAVE(1:3,1:3,1:2*CSMGPINDEX)=PTGP(1:3,1:3,1:2*CSMGPINDEX) ! before we change the index!
            PTGP(1:3,1:3,1:2*CSMGPINDEX)=PTGPGUIDE(1:3,1:3,1:2*CSMGPINDEX)
         ENDIF
         IF (PERMDIST) THEN
            DO J1=1,CSMGPINDEX
               XTEMP(1:3*NATOMS)=X(1:3*NATOMS)
            ENDDO
         ELSE
            DO J1=1,CSMGPINDEX
            ENDDO
         ENDIF

         IF (DEBUG) THEN
            SAVECSMIMAGES(1:3*NATOMS*CSMGPINDEX)=CSMIMAGES(1:3*NATOMS*CSMGPINDEX)
            
            DO J2=1,CSMGPINDEX
               DO J3=1,NATOMS
                  XTEMP(3*(J3-1)+1)=CSMPMAT(1,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) 
     &                             +CSMPMAT(1,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) 
     &                             +CSMPMAT(1,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                  XTEMP(3*(J3-1)+2)=CSMPMAT(2,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) 
     &                             +CSMPMAT(2,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) 
     &                             +CSMPMAT(2,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                  XTEMP(3*(J3-1)+3)=CSMPMAT(3,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) 
     &                             +CSMPMAT(3,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) 
     &                             +CSMPMAT(3,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
               ENDDO
            ENDDO
            IF (PERMDIST) THEN
               DO J1=1,CSMGPINDEX
                  XTEMP(1:3*NATOMS)=CSMAV(1:3*NATOMS)
               ENDDO
            ELSE
               DO J1=1,CSMGPINDEX
               ENDDO
            ENDIF
            AA(1)=0.0D0; AA(2)=0.0D0; AA(3)=6.283185307D0 ! should give an identity matrix
            SAVECSMPMAT(1:3,1:3)=CSMPMAT(1:3,1:3)
            SAVECSMNORM=CSMNORM
            DO J1=1,NATOMS
            ENDDO
            WRITE(MYUNIT,'(A,2G20.10)') 'potential> CSM values for reference structure and average=',EREAL,AVVAL
         ENDIF
         IF (CSMDOGUIDET) THEN ! undo guiding changes
            PTGP(1:3,1:3,1:2*CSMGPINDEX)=PTGPSAVE(1:3,1:3,1:2*CSMGPINDEX)
         ENDIF
      ELSE IF (PERMOPT) THEN
      ELSE IF (BLJCLUSTER) THEN
      ELSE IF (BINARY) THEN
         IF (SHIFTCUT) THEN
         ELSE
         ENDIF
      ELSE IF (LB2T) THEN
      ELSE IF (LJCOULT) THEN
      ELSE IF (STOCKT) THEN


      ELSE IF (DBPT) THEN

      ELSE IF (DBPTDT) THEN

      ELSE IF (DMBLMT) THEN

      ELSE IF (LINRODT) THEN

      ELSE IF (LWOTPT) THEN
         
          IF (PERIODIC) THEN
          ELSE
          ENDIF

      ELSE IF (NCAPT) THEN
 
      ELSE IF (NPAHT) THEN

      ELSE IF (NTIPT) THEN

      ELSEIF (PATCHY) THEN
      ELSEIF (ASAOOS) THEN ! this is not anisotropic

      ELSE IF (GBT) THEN

      ELSE IF (GBDT) THEN

      ELSE IF (GBDPT) THEN 

      ELSE IF (GEMT) THEN

      ELSE IF (PAHAT) THEN

      ELSE IF (PAHW99T) THEN

      ELSE IF (PAPT) THEN

      ELSE IF (MULTPAHAT) THEN

      ELSE IF (PYGT) THEN

      ELSE IF (PYGDPT) THEN
 

      ELSE IF (MSPYGT) THEN

      ELSE IF (MULTISITEPYT) THEN

      ELSE IF (MSTBINT) THEN

      ELSE IF (MSSTOCKT) THEN
          
      ELSE IF (SILANET) THEN

      ELSE IF (STOCKAAT) THEN

      ELSE IF (TDHDT) THEN


      ELSE IF (WATERDCT) THEN


      ELSE IF (WATERKZT) THEN


      ELSE IF (GAYBERNET) THEN
      ELSE IF (PARAMONOVT) THEN
      ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
      ELSE IF (LJCAPSIDT) THEN
      ELSE IF (STICKYT) THEN
      ELSE IF (TIP) THEN
      ELSE IF (QUADT) THEN
      ELSE IF (CAPSID) THEN
      ELSE IF (STRANDT) THEN
      ELSE IF (PAHT) THEN
      ELSE IF (DIFFRACTT) THEN
         WRITE(MYUNIT,'(A)') 'potential> diffract subroutine is commented'
         STOP

      ELSE IF (THOMSONT) THEN
         IF (ODDCHARGE.EQ.1.0D0) THEN
         ELSE
         ENDIF
      ELSE IF (GAUSST) THEN
      ELSE IF (QDT) THEN
      ELSE IF (QD2T) THEN
      ELSE IF (MULLERBROWNT) THEN
      ELSE IF (JMT) THEN
      ELSE IF (DF1T) THEN
      ELSE
         IF (EVAPREJECT) return
         IF (CUTT) THEN
         ELSE
         ENDIF
      ENDIF

      IF (PULLT) THEN
         EREAL=EREAL-PFORCE*(X(3*(PATOM1-1)+3)-X(3*(PATOM2-1)+3))
         GRAD(3*(PATOM1-1)+3)=GRAD(3*(PATOM1-1)+3)-PFORCE
         GRAD(3*(PATOM2-1)+3)=GRAD(3*(PATOM2-1)+3)+PFORCE
      ENDIF
      IF (FIELDT) THEN
         IF (CENT) THEN
         ELSE
         ENDIF
      ELSE IF (CPMD.AND.(NPCALL.EQ.1)) THEN
      ELSE IF (CPMD.AND.(.NOT.(SCT))) THEN
         INQUIRE(FILE='RESTART.1',EXIST=YESNO)
         OPEN(UNIT=8,FILE='newgeom',STATUS='UNKNOWN')
         DO J1=1,NATOMS
            WRITE(8,'(6F20.10)') X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3),0.0D0,0.0D0,0.0D0
         ENDDO
         IF (.NOT.YESNO) THEN
         ELSE
         ENDIF
     1               '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD call/" > temp')
         OPEN (UNIT=7,FILE='temp',STATUS='OLD')
         READ(7,'(A)') FNAME
         WRITE(MYUNIT,'(A)') FNAME
         OPEN (UNIT=7,FILE='ENERGY',STATUS='OLD')
         READ(7,*) EREAL, GEMAX
         IF (GEMAX.GT.1.0D-5) WRITE(MYUNIT,'(A,G15.5,A)') 'WARNING, GEMAX=',GEMAX,' CPMD wavefunction convergence suspect'
         OPEN(UNIT=7,FILE='GEOMETRY',STATUS='OLD')
         DO J1=1,NATOMS
            READ(7,*) GEMAX,GEMAX,GEMAX,GRAD(3*(J1-1)+1),GRAD(3*(J1-1)+2),GRAD(3*(J1-1)+3)
            GRAD(3*(J1-1)+1)=-GRAD(3*(J1-1)+1)
            GRAD(3*(J1-1)+2)=-GRAD(3*(J1-1)+2)
            GRAD(3*(J1-1)+3)=-GRAD(3*(J1-1)+3)
         ENDDO
      ELSE IF (DL_POLY) THEN
         OPEN (UNIT=91,FILE='STATIS',STATUS='OLD')
         READ(91,*) DUMM
         READ(91,*) DUMM
         READ(91,*) DUMM
         READ(91,*) EREAL
         WRITE(MYUNIT,'(A,G20.10)' ) 'EREAL=',EREAL
         OPEN (UNIT=91,FILE='STATIS',STATUS='OLD')
         READ(91,*) (GRAD(J3),J3=1,3*NATOMS)
         WRITE(91,'(3F20.10)') (GRAD(J3),J3=1,3*NATOMS)
      ENDIF

      IF (COMPON) CALL COMPRESS(X,GRAD,EREAL,GRADT)

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

         IF (SEEDT.AND.FREEZECORE) THEN
            DO J3=3*(NATOMS-NSEED)+1,3*NATOMS
               GRAD(J3)=0.0D0
            ENDDO
         ENDIF
         RMS=0.0D0
         DUMMY2=0.0D0
         IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
         ELSEIF (.NOT.THOMSONT) THEN
            DO J1=1,3*NATOMS
               DUMMY2=DUMMY2+GRAD(J1)**2 
            ENDDO
            RMS=MAX(DSQRT(DUMMY2/(3*NATOMS)),1.0D-100)
         ELSE
            DO J1=1,2*NATOMS
               DUMMY2=DUMMY2+GRAD(J1)**2 
            ENDDO
            RMS=MAX(DSQRT(DUMMY2/(2*NATOMS)),1.0D-100)
         ENDIF
         IF (NATBT.AND.(RMS.LT.GUIDECUT).AND.GUPTAT) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching off Gupta'
            GUPTAT=.FALSE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF (CSMDOGUIDET) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)') 'potential> switching off guiding point group'
            WRITE(MYUNIT,'(A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu guide      E=',
     1               EREAL,' steps=',CSMIT,' RMS=',CSMRMS
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF (CPMD.AND.(RMS.LT.GUIDECUT).AND.SCT) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching off Sutton-Chen'
            SCT=.FALSE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF (WELCH.AND.TOSI.AND.(RMS.LT.GUIDECUT)) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching off Tosi'
            TOSI=.FALSE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF ((ZETT1.OR.ZETT2).AND.(RMS.LT.GUIDECUT).AND.MORSET) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching off MORSE'
            MORSET=.FALSE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF (PACHECO.AND.(RMS.LT.GUIDECUT).AND.(.NOT.AXTELL)) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)' ) 'switching on AXTELL'
            AXTELL=.TRUE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF (PERCOLATET .AND. COMPON .AND. (RMS.LT.GUIDECUT)) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A, G20.10, A, G20.10)' )
     & 'Switching off compression, GUIDECUT=', GUIDECUT, '; RMS=', RMS
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
      ENDIF
 

      RETURN
      END
