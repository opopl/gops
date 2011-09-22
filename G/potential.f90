!op226>=================================== 
!op226> GPL License info {{{ 
!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!op226>}}} 
!op226>=================================== 
      SUBROUTINE POTENTIAL(X,GRAD,EREAL,GRADT,SECT)
!op226> Declarations {{{ 
      ! modules {{{
      USE COMMONS
      USE V
      USE MODBLN
      USE QMODULE
      USE PERMU
      USE PORFUNCS
      !}}}
      IMPLICIT NONE
      ! subroutine {{{
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: X
      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: GRAD
      DOUBLE PRECISION, INTENT(OUT) :: EREAL
      LOGICAL, intent(in) :: GRADT, SECT       
      ! }}}
      ! local {{{
      INTEGER BRUN
      INTEGER NQTOT
      DOUBLE PRECISION :: EA(10)
      COMMON /EA/ EA

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
      COMMON /CO/ COMPON
      COMMON /FAIL/ FTEST
      COMMON /EV/ EVAP, EVAPREJECT
      COMMON /PCALL/ NPCALL
      COMMON /CSMAVVAL/ AVVAL, CSMRMS, CSMIT
      COMMON /TOT/ NQTOT

      ! }}}
      ! }}}
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

      ! IF..ELSE chain {{{

      IF (MSORIGT) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         IF (CUTT) THEN
             CALL MSORIGC(NATOMS,X,GRAD,EREAL,GRADT)
         ELSE
             CALL MSORIG(NATOMS,X,GRAD,EREAL,GRADT)
         ENDIF
         IF (FTEST) THEN
            RETURN
         ENDIF
      ELSE IF (MSTRANST) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL MSTRANS97(NATOMS,X,GRAD,EREAL,GRADT)
         IF (FTEST) THEN
            RETURN
         ENDIF
      ELSE IF (FRAUSIT) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL FRAUSI(NATOMS,X,GRAD,EREAL,GRADT,ANGST,NATOMS)
         IF (FTEST) THEN
            RETURN
         ENDIF
!
!  DIM Ne^+, Ne*, Ar^+, Ar*
!
      ELSE IF ((NEON.AND.PLUS).OR.(NEON.AND.STAR).OR.(ARGON.AND.PLUS).OR.(ARGON.AND.STAR)) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
!        CALL RGNI(NATOMS,X,GRAD,EREAL,GRADT,h0,h1,ee,ev,w,NATOMS)
         CALL RGNI(NATOMS,X,GRAD,EREAL,GRADT)
!
!  DIM Ar^{2+}
!
      ELSE IF (TWOPLUS) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
!        CALL RGNII(NATOMS,X,GRAD,EREAL,GRADT,h0,h1,ee,ev,w,NATOMS)
         CALL RGNII(NATOMS,X,GRAD,EREAL,GRADT)
!
!  DIM Ar and Ne neutrals.
!
      ELSE IF (GROUND) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL GRND(NATOMS,X,GRAD,EREAL,GRADT)
         IF (AXTELL) CALL AXT(NATOMS,X,GRAD,EREAL,GRADT,ZSTAR)
      ELSE IF (RGCL2) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL RGNX2(NATOMS,X,GRAD,EREAL,GRADT)
         IF (AXTELL) CALL AXT(NATOMS,X,GRAD,EREAL,GRADT,ZSTAR)
      ELSE IF (ARNO) THEN
         SOCOUPLE=.TRUE.
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL RGNXY(NATOMS,X,GRAD,EREAL,GRADT,SOCOUPLE)
!        PRINT*,'SOCOUPLE=',SOCOUPLE
!        PRINT*,'POINTS:'
!        WRITE(*,'(3F20.10)') (X(J1),J1=1,3*NATOMS)
!        PRINT*,'Energy=',EREAL
!        PRINT*,'Gradient:'
!        WRITE(*,'(3F20.10)') (GRAD(J1),J1=1,3*NATOMS)
!        IF (GRADT) THEN
!           DISP=1.0D-4
!           DO J1=1,3*NATOMS
!              DUMMY2=X(J1)
!              X(J1)=X(J1)+DISP
!              CALL RGNXY(NATOMS,X,GRADNUM,V1,.FALSE.,SOCOUPLE)
!              X(J1)=X(J1)-2.0D0*DISP
!              CALL RGNXY(NATOMS,X,GRADNUM,V2,.FALSE.,SOCOUPLE)
!              GRADNUM(J1)=(V1-V2)/(2.0D0*DISP)
!              X(J1)=DUMMY2
!           ENDDO
!        ENDIF
!        PRINT*,'Numerical derivatives for displacement ',DISP
!        WRITE(*,'(3F20.10)') (GRADNUM(J1),J1=1,3*NATOMS)
!        PRINT*,'Analytic/numerical derivatives'
!        WRITE(*,'(3F20.10)') (GRAD(J1)/GRADNUM(J1),J1=1,3*NATOMS)
!        STOP

         IF (AXTELL) CALL AXT(NATOMS,X,GRAD,EREAL,GRADT,ZSTAR)
      ELSE IF (OTPT) THEN
!         CALL OTP(X,GRAD,EREAL,GRADT,SECT)
!         DIFF=1.0D-4
!         PRINT*,'analytic and numerical gradients:'
!         DO J1=1,3*NATOMS
!            X(J1)=X(J1)+DIFF
!            CALL OTP(X,GRADDUM,EPLUS,.FALSE.,.FALSE.)
!            X(J1)=X(J1)-2.0D0*DIFF
!            CALL OTP(X,GRADDUM,EMINUS,.FALSE.,.FALSE.)
!            X(J1)=X(J1)+DIFF
!C           IF ((ABS(GRAD(J1)).NE.0.0D0).AND.(100.0D0*(GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1).GT.1.0D0)) THEN
!               WRITE(*,'(I5,2F20.10)') J1,GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!C           ENDIF
!         ENDDO
      ELSE IF (LJMFT) THEN
!         CALL LJ(X,GRADLJ,EREALLJ,GRADT,SECT)
!         CALL MF(X,GRADMF,EREALMF,GRADT)
!C        CALL LJ(X,GRADLJ,EREALLJ,.FALSE.,.FALSE.)
!C        CALL MF(X,GRADMF,EREALMF,.FALSE.)
!         WRITE(*,'(A,G20.10)') 'radius=',RADIUS
!         PRINT*,'EREALLJ,EREALMF=',EREALLJ,EREALMF
!         EREAL=EREALLJ+EREALMF
!         TERMLJ=EREALLJ
!         TERMMF=EREALMF
!         DO J1=1,3*NATOMS
!C           PRINT*,'J1,GRADLJ,GRADMF=',J1,GRADLJ(J1),GRADMF(J1)
!            GRAD(J1)=GRADLJ(J1)+GRADMF(J1)
!         ENDDO
!C         DIFF=1.0D-6
!CC        PRINT*,'analytic and numerical gradients:'
!C
!C         DO J1=1,3*NATOMS
!C            X(J1)=X(J1)+DIFF
!C            CALL LJ(X,GRADDUM1,EDUM1,.FALSE.,.FALSE.)
!C            CALL MF(X,GRADDUM2,EDUM2,.FALSE.)
!C            EPLUS=EDUM1+EDUM2
!C
!C            X(J1)=X(J1)-2.0D0*DIFF
!C
!C            CALL LJ(X,GRADDUM1,EDUM1,.FALSE.,.FALSE.)
!C            CALL MF(X,GRADDUM2,EDUM2,.FALSE.)
!C
!C            EMINUS=EDUM1+EDUM2
!C
!CC            X(J1)=X(J1)+DIFF
!C            GRAD(J1)=(EPLUS-EMINUS)/(2.0D0*DIFF)
!C            WRITE(*,'(I5,4F20.10)') J1,GRAD(J1),EPLUS,EMINUS,ABS(EPLUS-EMINUS)
!C         ENDDO
!
      ELSE IF (MORSET) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL MORSE(X,GRAD,EREAL,GRADT)
      ELSE IF (TOSI) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         if (evapreject) return
         CALL TOSIFUMI(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (WELCH) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL WEL(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (SCT) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL SC(X,GRAD,EREAL,GRADT)
         IF (CPMD) EREAL=EREAL+1.0D6
      ELSE IF (ACKLANDT) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL ACK(X,GRAD,EREAL,GRADT)
      ELSE IF (FAL.OR.FNI) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL FARKAS(X,GRAD,EREAL,GRADT,NATOMS)
      ELSE IF (DFTBT) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
!        CALL SECDIFF(NATOMS,X,w)
!        STOP
         CALL DFTB(NATOMS,X,GRAD,EREAL,GRADT)
         IF (FTEST) THEN
            RETURN
         ENDIF
      ELSE IF (MYBLNT) THEN
         CALL EBLN(NATOMS,X,EA,GRAD,HESS,BLNTYPE,GRADT,SECT)
         EREAL=EA(1)
      ELSE IF (P46) THEN
         CALL P46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
      ELSE IF (G46) THEN
         CALL G46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
      ELSE IF (BLNT) THEN
         CALL BLN(X,GRAD,EREAL,GRADT)
      ELSE IF (CHAPERONINT) THEN
         CALL CHAPERONIN(X,GRAD,EREAL,GRADT)
      ELSE IF (SW) THEN
         IF (PERIODIC) THEN
            CALL SISW(X,GRAD,EREAL,GRADT)
         ELSE
            CALL RAD(X,GRAD,EREAL,GRADT)
            CALL SWTWO(X,GRAD,EREAL,GRADT)
            CALL SWTHREE(GRAD,EREAL,GRADT)
         ENDIF
!     ELSE IF (AMBER) THEN
!        CALL AMB(X,GRAD,EREAL,GRADT)
      ELSE IF (AMHT) THEN
         CALL WALESAMH_INTERFACE(X,GRAD,EREAL)
!         DIFF=1.0D-4
!         SUMMDIFF=0.D0
!         PRINT*,'analytic and numerical gradients:'
!         DO J1=1,3*NATOMS
!           WRITE(*,'(F20.10,2x,F20.10,2xI5)')X(J1),GRAD(J1),J1

!          X(J1)=X(J1)+DIFF
!          CALL scltab_wales(X,GRADDUM,EPLUS)
!          X(J1)=X(J1)-2.0D0*DIFF
!          CALL scltab_wales(X,GRADDUM,EMINUS)
!          X(J1)=X(J1)+DIFF
!           IF (100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1)).gt.0.0D0)THEN
!                WRITE(*,'(I5,3F15.8)') J1,GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),
!     1                          100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1))
!          ENDIF
!            SUMMDIFF = SUMMDIFF + (GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1)
!         GRAD(J1)=(EPLUS-EMINUS)/(2.0D0*DIFF)
!         ENDDO
!           WRITE(6,*)'SUMM DIFF ', SUMMDIFF

! sf344> addition
      ELSE IF (AMBERT) THEN
         call amberenergies(X,GRAD,EREAL,.false.,.false.)
      ELSE IF (CHRMMT) THEN
!        WRITE(LFH,'(A)') 'potential> coords:'
!        WRITE(LFH,'(3G20.10)') X(1:3*NATOMS)
         CALL OCHARMM(X,GRAD,EREAL,GRADT)
      ELSE IF (DZTEST) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL DZPOT(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (ZETT1) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL Z1(X,GRAD,EREAL,GRADT)
      ELSE IF (ZETT2) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL Z2(X,GRAD,EREAL,GRADT)
      ELSE IF (PACHECO) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL PRC60(X,GRAD,EREAL,GRADT)
         IF (AXTELL) CALL AXT(NATOMS,X,GRAD,EREAL,GRADT,ZSTAR)
      ELSE IF (MODEL1T) THEN
         CALL MODEL1(X,GRAD,EREAL,QE,QX)
      ELSE IF (EAMLJT) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL EAMLJ(X,GRAD,EREAL,GRADT)
      ELSE IF (PBGLUET) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL PBGLUE(X,GRAD,EREAL,GRADT)
      ELSE IF (FST) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL FS(X,GRAD,EREAL,GRADT)
      ELSE IF (GUPTAT) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL GUPTA(X,GRAD,EREAL,GRADT)
      ELSE IF (NATBT) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL NATB(NATOMS,X,GRAD,EREAL,GRADT,GUIDET)
!        DIFF=1.0D-3
!        PRINT*,'analytic and numerical gradients:'
!        DO J1=1,3*NATOMS
!           X(J1)=X(J1)+DIFF
!           CALL NATB(NATOMS,X,GRADDUM,EPLUS,.FALSE.)
!           X(J1)=X(J1)-2.0D0*DIFF
!           CALL NATB(NATOMS,X,GRADDUM,EMINUS,.FALSE.)
!           X(J1)=X(J1)+DIFF
!           IF (GRAD(J1).NE.0.0D0) WRITE(*,'(I5,3F20.10)') J1,GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),
!    1                          100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1))
!        ENDDO
!        STOP
      ELSE IF (ALGLUET) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL ALGLUE(X,GRAD,EREAL,GRADT)
      ELSE IF (MGGLUET) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL MGGLUE(X,GRAD,EREAL,GRADT)
      ELSE IF (EAMALT) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL EAMAL(X,GRAD,EREAL,GRADT)
      ELSE IF (WENZEL) THEN
         CALL WEN(X,GRAD,EREAL,GRADT,NATOMS)
      ELSE IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
!        IF (DEBUG) OPEN(UNIT=765,FILE='CSMrot.xyz',STATUS='UNKNOWN')
         IF (CSMDOGUIDET) THEN ! we want to use a guiding group
            PTGPSAVE(1:3,1:3,1:2*CSMGPINDEX)=PTGP(1:3,1:3,1:2*CSMGPINDEX) ! before we change the index!
            CSMNORMSAVE=CSMNORM
            CSMNORM=CSMGUIDENORM
            CSMGPSAVE=CSMGP
            CSMGP=CSMGUIDEGP
            CSMGPINDEXSAVE=CSMGPINDEX
            CSMGPINDEX=CSMGUIDEGPINDEX
            PTGP(1:3,1:3,1:2*CSMGPINDEX)=PTGPGUIDE(1:3,1:3,1:2*CSMGPINDEX)
         ENDIF
         IF (PERMDIST) THEN
            DO J1=1,CSMGPINDEX
               XTEMP(1:3*NATOMS)=X(1:3*NATOMS)
               CALL CSMROT(XTEMP,DUMMY,1,J1)
               CALL MINPERMDIST(XTEMP,DUMMY,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,EREAL,DIST2,RIGID,RMAT)
               CALL CSMROT(DUMMY,XTEMP,-1,J1) ! need to rotate the permuted rotated images back to the reference orientation
               CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=XTEMP(1:3*NATOMS)
            ENDDO
         ELSE
            DO J1=1,CSMGPINDEX
               CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=X(1:3*NATOMS)
            ENDDO
         ENDIF
         CALL CSMMIN(X,EREAL,CSMRMS,CSMIT)

         IF (DEBUG) THEN
!
!  Saving CSMIMAGES, CSMPMAT is necessary here because otherwise these quantities will
!  be different in finalq when we need to write out the results to file CSMav.xyz.
!
            CSMAV(1:3*NATOMS)=0.0D0
!              WRITE(LFH,'(A)') 'potential> CSMPMAT:'
!              WRITE(LFH,'(3G20.10)') CSMPMAT(1:3,1:3)
            SAVECSMIMAGES(1:3*NATOMS*CSMGPINDEX)=CSMIMAGES(1:3*NATOMS*CSMGPINDEX)
            
            DO J2=1,CSMGPINDEX
!
! Rotate permuted image to best orientation with CSMPMAT
! Then apply point group operation J2
!
               DO J3=1,NATOMS
                  XTEMP(3*(J3-1)+1)=CSMPMAT(1,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) &
     &                             +CSMPMAT(1,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) &
     &                             +CSMPMAT(1,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                  XTEMP(3*(J3-1)+2)=CSMPMAT(2,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) & 
     &                             +CSMPMAT(2,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) &
     &                             +CSMPMAT(2,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
                  XTEMP(3*(J3-1)+3)=CSMPMAT(3,1)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+1) &
     &                             +CSMPMAT(3,2)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+2) &
     &                             +CSMPMAT(3,3)*CSMIMAGES(3*NATOMS*(J2-1)+3*(J3-1)+3)
               ENDDO
!              WRITE(LFH,'(A)') 'potential> XTEMP:'
!              WRITE(LFH,'(3G20.10)') XTEMP(1:3*NATOMS)
               CALL CSMROT(XTEMP,DUMMY,1,J2)
               CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)+DUMMY(1:3*NATOMS)
            ENDDO
            CSMAV(1:3*NATOMS)=CSMAV(1:3*NATOMS)/CSMGPINDEX
!
!  Check the CSM for the averaged structure. It should be zero if this structure has the
!  right point group. Need to reset CSMIMAGES and CSMNORM temporarily.
!
            IF (PERMDIST) THEN
               DO J1=1,CSMGPINDEX
                  XTEMP(1:3*NATOMS)=CSMAV(1:3*NATOMS)
                  CALL CSMROT(XTEMP,DUMMY,1,J1)
                  CALL MINPERMDIST(XTEMP,DUMMY,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,DUMMY2,DIST2,RIGID,RMAT)
                  CALL CSMROT(DUMMY,XTEMP,-1,J1) ! need to rotate the permuted rotated images back to the reference orientation
                  CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=XTEMP(1:3*NATOMS)
               ENDDO
            ELSE
               DO J1=1,CSMGPINDEX
                  CSMIMAGES(1+3*NATOMS*(J1-1):3*NATOMS*J1)=CSMAV(1:3*NATOMS)
               ENDDO
            ENDIF
            AA(1)=0.0D0; AA(2)=0.0D0; AA(3)=6.283185307D0 ! should give an identity matrix
            SAVECSMPMAT(1:3,1:3)=CSMPMAT(1:3,1:3)
            SAVECSMNORM=CSMNORM
            CSMNORM=0.0D0
            DO J1=1,NATOMS
               CSMNORM=CSMNORM+CSMAV(3*(J1-1)+1)**2+CSMAV(3*(J1-1)+2)**2+CSMAV(3*(J1-1)+3)**2
            ENDDO
            CSMNORM=2*CSMGPINDEX*CSMNORM
            CALL CSMPOTGRAD(CSMAV,AA,AVVAL,.TRUE.,CSMGRAD)
            CSMNORM=SAVECSMNORM
            CSMPMAT(1:3,1:3)=SAVECSMPMAT(1:3,1:3)
            CSMIMAGES(1:3*NATOMS*CSMGPINDEX)=SAVECSMIMAGES(1:3*NATOMS*CSMGPINDEX)
            WRITE(LFH,'(A,2G20.10)') 'potential> CSM values for reference structure and average=',EREAL,AVVAL
         ENDIF
         IF (CSMDOGUIDET) THEN ! undo guiding changes
            CSMGP=CSMGPSAVE
            CSMGPINDEX=CSMGPINDEXSAVE
            PTGP(1:3,1:3,1:2*CSMGPINDEX)=PTGPSAVE(1:3,1:3,1:2*CSMGPINDEX)
            CSMNORM=CSMNORMSAVE
         ENDIF
      ELSE IF (PERMOPT) THEN
!
!  EREAL is the distance in this case
!
         CALL MINPERMDIST(FIN,X,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,PERIODIC,TWOD,EREAL,DIST2,RIGID,RMAT)
      ELSE IF (BLJCLUSTER) THEN
         CALL RAD(X,GRAD,EREAL,GRADT)
         CALL LJPSHIFTBINC(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (BINARY) THEN
         IF (SHIFTCUT) THEN
            CALL LJPSHIFT(X,GRAD,EREAL,GRADT,SECT)
         ELSE
            CALL LJPBIN(X,GRAD,EREAL,GRADT,SECT)
         ENDIF
      ELSE IF (LB2T) THEN
         CALL RADR(X,GRAD,EREAL,GRADT)
         CALL LB2(X,GRAD,EREAL,GRADT)
!         CALL LB2(X,GRAD,EREAL,GRADT,SECT)
!         DIFF=1.0D-4 
!         PRINT*,'analytic and numerical gradients:' 
!         DO J1=1,3*NATOMS
!            X(J1)=X(J1)+DIFF
!            CALL LB2(X,GRADDUM,EPLUS,.FALSE.,.FALSE.)
!            X(J1)=X(J1)-2.0D0*DIFF
!            CALL LB2(X,GRADDUM,EMINUS,.FALSE.,.FALSE.)
!            X(J1)=X(J1)+DIFF
!C           IF ((ABS(GRAD(J1)).NE.0.0D0).AND.(100.0D0*(GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1).GT.1.0D0)) THEN
!            WRITE(*,'(I5,2F20.10)') J1,GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!C           ENDIF
!          ENDDO
      ELSE IF (LJCOULT) THEN
         CALL RADR(X,GRAD,EREAL,GRADT)
         CALL LJCOUL(X,GRAD,EREAL,GRADT)
      ELSE IF (STOCKT) THEN
         CALL RADR(X,GRAD,EREAL,GRADT)
         CALL STOCK(X,GRAD,EREAL,GRADT)
!        DIFF=1.0D-4
!        PRINT*,'analytic and numerical gradients:'
!        DO J1=1,3*NATOMS
!           X(J1)=X(J1)+DIFF
!           CALL STOCK(X,GRADDUM,EPLUS,.FALSE.,.FALSE.)
!           X(J1)=X(J1)-2.0D0*DIFF
!           CALL STOCK(X,GRADDUM,EMINUS,.FALSE.,.FALSE.)
!           X(J1)=X(J1)+DIFF
!           IF (GRAD(J1).NE.0.0D0) WRITE(*,'(I5,3F20.10)') J1,GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),
!    1                          100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1))
!        ENDDO
!
!       Anisotropic potentials:

!     DC430 >

      ELSE IF (DBPT) THEN
          CALL DUMBBELLP (X, GRAD, EREAL, GRADT)

      ELSE IF (DBPTDT) THEN
          CALL DMBLTD (X, GRAD, EREAL, GRADT)

      ELSE IF (DMBLMT) THEN
          CALL DMBLMORSE (X, GRAD, EREAL, GRADT)

      ELSE IF (LINRODT) THEN
          CALL LINROD (X, GRAD, EREAL, GRADT)

      ELSE IF (LWOTPT) THEN
         
          IF (PERIODIC) THEN
             CALL LWOTP(X, GRAD, EREAL, GRADT)
          ELSE
             CALL LWOTPC(X, GRAD, EREAL, GRADT)
          ENDIF

      ELSE IF (NCAPT) THEN
          CALL NEWCAPSID (X, GRAD, EREAL, GRADT)
 
      ELSE IF (NPAHT) THEN
          CALL NEWPAH (X, GRAD, EREAL, GRADT)

      ELSE IF (NTIPT) THEN
          CALL NEWTIP (X, GRAD, EREAL, GRADT)

!|gd351>
      ELSEIF (PATCHY) THEN
          CALL PATCHYPOT (X, GRAD, EREAL, GRADT)
      ELSEIF (ASAOOS) THEN ! this is not anisotropic
          CALL ASAOOSPOT (X, GRAD, EREAL, GRADT)
!<gd351|

      ELSE IF (GBT) THEN
          CALL GBCALAMITIC (X, GRAD, EREAL, GRADT)

      ELSE IF (GBDT) THEN
          CALL GBDISCOTIC (X, GRAD, EREAL, GRADT)

      ELSE IF (GBDPT) THEN 
          CALL GBDDP (X, GRAD, EREAL, GRADT)

      ELSE IF (GEMT) THEN
          CALL GEM (X, GRAD, EREAL, GRADT)

      ELSE IF (PAHAT) THEN
          CALL PAHA (X, GRAD, EREAL, GRADT)

      ELSE IF (PAHW99T) THEN
          CALL PAHW99 (X, GRAD, EREAL, GRADT)

      ELSE IF (PAPT) THEN
          CALL PAP (X, GRAD, EREAL, GRADT)

      ELSE IF (MULTPAHAT) THEN
          CALL MULTPAHA (X, GRAD, EREAL, GRADT)

      ELSE IF (PYGT) THEN
          CALL PYG (X, GRAD, EREAL, GRADT)

      ELSE IF (PYGDPT) THEN
          CALL PYGDP (X, GRAD, EREAL, GRADT)
 
!      ELSE IF (MSGBT) THEN
!          CALL MSGAYBERNE (X, GRAD, EREAL, GRADT) 

      ELSE IF (MSPYGT) THEN
          CALL MULTISITEPY (X, GRAD, EREAL, GRADT)

      ELSE IF (MULTISITEPYT) THEN
          CALL MULTISITEPY2 (X, GRAD, EREAL, GRADT)

      ELSE IF (MSTBINT) THEN
          CALL MSTBIN (X, GRAD, EREAL, GRADT)

      ELSE IF (MSSTOCKT) THEN
          CALL MULTSTOCK (X, GRAD, EREAL, GRADT)
          
      ELSE IF (SILANET) THEN
          CALL SILANE (X, GRAD, EREAL, GRADT)

      ELSE IF (STOCKAAT) THEN
          CALL STOCKAA (X, GRAD, EREAL, GRADT)

      ELSE IF (TDHDT) THEN

          CALL TETRAHEDRA (X, GRAD, EREAL, GRADT)

      ELSE IF (WATERDCT) THEN

          CALL WATERPDC (X, GRAD, EREAL, GRADT)

      ELSE IF (WATERKZT) THEN

          CALL WATERPKZ (X, GRAD, EREAL, GRADT)

!    sf344>
      ELSE IF (GAYBERNET) THEN
!         CALL RADR(X,GRAD,EREAL,GRADT)
         CALL GAYBERNE(X,GRAD,EREAL,GRADT)
      ELSE IF (PARAMONOVT) THEN
!        CALL PARAMONOVNUMFIRSTDER(X,GRADT)
         CALL OLDPARAMONOV(X,GRAD,EREAL,GRADT) ! old PY routine
      ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
! call PES routine
!        CALL PYPES(X)      
!        STOP
        CALL PYGPERIODIC(X,GRAD,EREAL,GRADT)
!        CALL PARAMONOVNUMFIRSTDER(X,GRADT)
      ELSE IF (LJCAPSIDT) THEN
        CALL LJCAPSIDMODEL(X,GRAD,EREAL,GRADT)
      ELSE IF (STICKYT) THEN
         CALL RADR(X,GRAD,EREAL,GRADT)
         CALL STICKY(X,GRAD,EREAL,GRADT,.FALSE.)
!        DIFF=1.0D-3
!        PRINT*,'analytic and numerical gradients:'
!        DO J1=1,3*NATOMS
!           X(J1)=X(J1)+DIFF
!           CALL STICKY(X,GRADDUM,EPLUS,.FALSE.,.FALSE.)
!           X(J1)=X(J1)-2.0D0*DIFF
!           CALL STICKY(X,GRADDUM,EMINUS,.FALSE.,.FALSE.)
!           X(J1)=X(J1)+DIFF
!           IF (GRAD(J1).NE.0.0D0) WRITE(*,'(I5,3F20.10)') J1,GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF),
!    1                          100*ABS((GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1))
!        ENDDO
!        STOP
      ELSE IF (TIP) THEN
         CALL RADR(X,GRAD,EREAL,GRADT)
         CALL TIP4P(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (QUADT) THEN
         CALL RADR(X,GRAD,EREAL,GRADT)
         CALL QUAD(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (CAPSID) THEN
!        CALL RADC(X,GRAD,EREAL,GRADT)
         CALL FCAPSID(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (STRANDT) THEN
         CALL STRAND(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (PAHT) THEN
         CALL PAH(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (DIFFRACTT) THEN
         WRITE(LFH,'(A)') 'potential> diffract subroutine is commented'
         STOP
!        CALL DIFFRACT(X,GRAD,EREAL,GRADT,SECT,NATOMS)

!        DIFF=1.0D-4
!        PRINT*,'analytic and numerical gradients:'
!        DO J1=1,3*NATOMS
!           X(J1)=X(J1)+DIFF
!           CALL DIFFRACT(X,GRADDUM,EPLUS,.FALSE.,.FALSE.)
!           X(J1)=X(J1)-2.0D0*DIFF
!           CALL DIFFRACT(X,GRADDUM,EMINUS,.FALSE.,.FALSE.)
!           X(J1)=X(J1)+DIFF
!           IF ((ABS(GRAD(J1)).NE.0.0D0).AND.(100.0D0*(GRAD(J1)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GRAD(J1).GT.-1.0D0)) 
!    1         WRITE(*,'(I5,2F20.10)') J1,GRAD(J1),(EPLUS-EMINUS)/(2.0D0*DIFF)
!           PRINT '(A,3G20.10)','function: ',EREAL
!           PRINT '(A,3G20.10)','variables: ',X(1:3)
!           PRINT '(A,3G20.10)','gradient:  ',GRAD(1:3)
!        ENDDO
      ELSE IF (THOMSONT) THEN
         IF (ODDCHARGE.EQ.1.0D0) THEN
            CALL THOMSON(X,GRAD,EREAL,GRADT)
         ELSE
            CALL ODDTHOMSON(X,GRAD,EREAL,GRADT)
         ENDIF
      ELSE IF (GAUSST) THEN
         CALL GFIELD(X,GRAD,EREAL,GRADT)
      ELSE IF (QDT) THEN
         CALL QDTEST(X,GRAD,EREAL,GRADT)
      ELSE IF (QD2T) THEN
         CALL QDTEST2(X,GRAD,EREAL,GRADT)
      ELSE IF (MULLERBROWNT) THEN
         CALL MB(X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (JMT) THEN
         CALL JMEC(NATOMS,X,GRAD,EREAL,GRADT,SECT)
      ELSE IF (DF1T) THEN
         CALL DF1GRAD(X,NATOMS,GRAD,EREAL,GRADT,BOXLX,BOXLY)
      ELSE
!
!  RAD must be called before the routine that calculates the potential or LBFGS
!  will get confused even if EVAP is set .TRUE. correctly.
!
         CALL RAD(X,GRAD,EREAL,GRADT)
         IF (EVAPREJECT) return
         IF (CUTT) THEN
            CALL LJCUT(X,GRAD,EREAL,GRADT,SECT)
         ELSE
            CALL LJ(X,GRAD,EREAL,GRADT,SECT)
         ENDIF
      ENDIF
! }}}

!op226>}}} 
!
!  --------------- End of possible potentials - now add fields if required ------------------------------
      ! Fields {{{
!
      IF (PULLT) THEN
         EREAL=EREAL-PFORCE*(X(3*(PATOM1-1)+3)-X(3*(PATOM2-1)+3))
         GRAD(3*(PATOM1-1)+3)=GRAD(3*(PATOM1-1)+3)-PFORCE
         GRAD(3*(PATOM2-1)+3)=GRAD(3*(PATOM2-1)+3)+PFORCE
      ENDIF

      IF (FIELDT) THEN
         IF (CENT) THEN
            CALL FDM(X,GRAD,EREAL,GRADT)
         ELSE
            CALL FD(X,GRAD,EREAL,GRADT)
         ENDIF
      ELSE IF (CPMD.AND.(NPCALL.EQ.1)) THEN
         CALL SYSTEM(' sed -e "s/DUMMY/RESTART WAVEFUNCTION GEOFILE LATEST/" ' //  SYS(1:LSYS) // ' > temp ')
         CALL SYSTEM(' mv temp ' // SYS(1:LSYS) // '.restart')
      ELSE IF (CPMD.AND.(.NOT.(SCT))) THEN
         INQUIRE(FILE='RESTART.1',EXIST=YESNO)
         OPEN(UNIT=8,FILE='newgeom',STATUS='UNKNOWN')
         DO J1=1,NATOMS
            WRITE(8,'(6F20.10)') X(3*(J1-1)+1),X(3*(J1-1)+2),X(3*(J1-1)+3),0.0D0,0.0D0,0.0D0
         ENDDO
         CLOSE(8)
         CALL SYSTEM(' mv newgeom GEOMETRY ')
         CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.out ' // SYS(1:LSYS) // '.old.out >& /dev/null ')
         IF (.NOT.YESNO) THEN
            CALL SYSTEM(' ( cpmd.x.2 ' // SYS(1:LSYS) // ' > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
         ELSE
            CALL SYSTEM(' ( cpmd.x.2 ' // SYS(1:LSYS) // '.restart > ' // SYS(1:LSYS) // '.out ) >& /dev/null')
         ENDIF
         CALL SYSTEM('grep "CPU TIME" ' // SYS(1:LSYS) // &
     &               '.out | tail -1 | sed -e "s/ *CPU TIME/ CPU time for CPMD call/" > temp')
         OPEN (UNIT=7,FILE='temp',STATUS='OLD')
         READ(7,'(A)') FNAME
         WRITE(LFH,'(A)') FNAME
         CLOSE(7)
         OPEN (UNIT=7,FILE='ENERGY',STATUS='OLD')
         READ(7,*) EREAL, GEMAX
         CLOSE(7)
         IF (GEMAX.GT.1.0D-5) WRITE(LFH,'(A,G15.5,A)') 'WARNING, GEMAX=',GEMAX,' CPMD wavefunction convergence suspect'
         OPEN(UNIT=7,FILE='GEOMETRY',STATUS='OLD')
         DO J1=1,NATOMS
            READ(7,*) GEMAX,GEMAX,GEMAX,GRAD(3*(J1-1)+1),GRAD(3*(J1-1)+2),GRAD(3*(J1-1)+3)
!           WRITE(*,'(6F20.10)') GEMAX,GEMAX,GEMAX,GRAD(3*(J1-1)+1),GRAD(3*(J1-1)+2),GRAD(3*(J1-1)+3)
            GRAD(3*(J1-1)+1)=-GRAD(3*(J1-1)+1)
            GRAD(3*(J1-1)+2)=-GRAD(3*(J1-1)+2)
            GRAD(3*(J1-1)+3)=-GRAD(3*(J1-1)+3)
         ENDDO
         CLOSE(7)
      ELSE IF (DL_POLY) THEN
         CALL SYSTEM('DLPOLY.X > forces')
         OPEN (UNIT=91,FILE='STATIS',STATUS='OLD')
         READ(91,*) DUMM
         READ(91,*) DUMM
         READ(91,*) DUMM
         READ(91,*) EREAL
         WRITE(LFH,'(A,G20.10)' ) 'EREAL=',EREAL
         CLOSE(91)
         OPEN (UNIT=91,FILE='STATIS',STATUS='OLD')
         READ(91,*) (GRAD(J3),J3=1,3*NATOMS)
         CLOSE(91)
         WRITE(91,'(3F20.10)') (GRAD(J3),J3=1,3*NATOMS)
      ENDIF
!     IF (SQUEEZET) CALL SQUEEZE(X,GRAD,EREAL,GRADT)

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
!
!  Preserve centre of mass if required.
!
!        IF (CENT.AND.(.NOT.SEEDT)) THEN
!           XG=0.0D0
!           YG=0.0D0
!           ZG=0.0D0
!           DO J1=1,NATOMS-NSEED
!              J2=3*J1
!              XG=XG+GRAD(J2-2)
!              YG=YG+GRAD(J2-1)
!              ZG=ZG+GRAD(J2)
!           ENDDO
!           XG=XG/(NATOMS-NSEED)
!           YG=YG/(NATOMS-NSEED)
!           ZG=ZG/(NATOMS-NSEED)
!           DO J1=1,NATOMS-NSEED
!              J2=3*J1
!              GRAD(J2-2)=GRAD(J2-2)-XG
!              GRAD(J2-1)=GRAD(J2-1)-YG
!              GRAD(J2)=  GRAD(J2)-ZG
!           ENDDO
!        ENDIF
         RMS=0.0D0
         DUMMY2=0.0D0
         ! }}}

         RMS=MAX(DSQRT(SUM(GRAD**2)/3*NATOMS),1.0D-100)

! rest not-important {{{


         !IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
         !ELSEIF (.NOT.THOMSONT) THEN
            !DO J1=1,3*NATOMS
               !DUMMY2=DUMMY2+GRAD(J1)**2 
            !ENDDO
            !RMS=MAX(DSQRT(DUMMY2/(3*NATOMS)),1.0D-100)
         !ELSE
            !DO J1=1,2*NATOMS
               !DUMMY2=DUMMY2+GRAD(J1)**2 
            !ENDDO
            !RMS=MAX(DSQRT(DUMMY2/(2*NATOMS)),1.0D-100)
         !ENDIF

! Guiding potentials are set back to true in quench.f 
!
         IF (NATBT.AND.(RMS.LT.GUIDECUT).AND.GUPTAT) THEN
            IF (DEBUG) WRITE(LFH,'(A)' ) 'switching off Gupta'
            GUPTAT=.FALSE.
            GUIDECHANGET=.TRUE.
!
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
!
            GOTO 10
         ENDIF
!        IF ((CSMRMS.LT.GUIDECUT).AND.(CSMDOGUIDET)) THEN
         IF (CSMDOGUIDET) THEN
            IF (DEBUG) WRITE(LFH,'(A)') 'potential> switching off guiding point group'
            WRITE(LFH,'(A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu guide      E=',&
     &               EREAL,' steps=',CSMIT,' RMS=',CSMRMS
            GUIDECHANGET=.TRUE.
            CSMDOGUIDET=.FALSE.
!
!  The GOTO is needed here in case LBFGS sets CFLAG to TRUE, in which
!  case we stop with the wrong energy (but almost the right coordinates)
!
            GOTO 10
         ENDIF
         IF (CPMD.AND.(RMS.LT.GUIDECUT).AND.SCT) THEN
            IF (DEBUG) WRITE(LFH,'(A)' ) 'switching off Sutton-Chen'
            SCT=.FALSE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF (WELCH.AND.TOSI.AND.(RMS.LT.GUIDECUT)) THEN
            IF (DEBUG) WRITE(LFH,'(A)' ) 'switching off Tosi'
            TOSI=.FALSE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF ((ZETT1.OR.ZETT2).AND.(RMS.LT.GUIDECUT).AND.MORSET) THEN
            IF (DEBUG) WRITE(LFH,'(A)' ) 'switching off MORSE'
            MORSET=.FALSE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF (PACHECO.AND.(RMS.LT.GUIDECUT).AND.(.NOT.AXTELL)) THEN
            IF (DEBUG) WRITE(LFH,'(A)' ) 'switching on AXTELL'
            AXTELL=.TRUE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
         IF (PERCOLATET .AND. COMPON .AND. (RMS.LT.GUIDECUT)) THEN
            IF (DEBUG) WRITE(LFH,'(A, G20.10, A, G20.10)' ) &
     & 'Switching off compression, GUIDECUT=', GUIDECUT, '; RMS=', RMS
            COMPON=.FALSE.
            GUIDECHANGET=.TRUE.
            GOTO 10
         ENDIF
      ENDIF
!
!  When called from mylbfgs the dimension of X is only 2*NATOMS for the TWOD case.
!      
!     IF (TWOD) THEN
!        DO J1=1,3*NATOMS
!           J2=3*J1
!           X(J2)=0.0D0
!           IF (GRADT) GRAD(J2)=0.0D0
!        ENDDO
!     ENDIF
 
!     WRITE(*,'(A,3F20.10)') 'X,EREAL,RMS=',X(1),EREAL,RMS
!     WRITE(*,'(A,G20.10)') 'energy in potential=',EREAL
!     PRINT*,'coordinates:'
!     WRITE(*,'(3G20.10)') (X(J1),J1=1,3*NATOMS)
!     PRINT*,'gradient:'
!     WRITE(*,'(3G20.10)') (GRAD(J1),J1=1,3*NATOMS)
      ! }}}

      RETURN
      END
